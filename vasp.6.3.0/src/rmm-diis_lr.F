!#define debug
#include "symbol.inc"
MODULE rmm_diis_lr
  USE prec
CONTAINS
!************************ SUBROUTINE LINEAR_RESPONSE_DIIS **************
!
! this subroutine solves the linear response equation
!    ( H(0) - e(0) S(0) ) |phi(1)> = - |xi>
! where xi is usually calculated to be
!    |xi> = ( H(1) - e(0) S(1) ) |phi(1)> - e(1) S(0) |phi(0)>
! i.e. the perturbation resulting from a change of the Hamiltonian
!
! in principle there is a related variational principle that reads
! < phi(1) | xi > + < xi | phi(1) > + <phi(1)| H(0) - e(0) S(0) |phi(1)>
! which could be optimised as well, but this requires to constrain
! the wavefunctions phi(1) to observe certain orthogonality constraints
!
! in the present implementation an inverse iteration like algorithm
! is the basic step in the linear response solver 
! the routine is a variant of the rmm-diis.F routine
!
!  INFO%IALGO   determine type of preconditioning and the algorithm
!    8    TAP preconditioning
!    9    Jacobi like preconditioning
!    (TAP Teter Alan Payne is presently hardcoded)
!  WEIMIN  treshhold for total energy minimisation
!    is the fermiweight of a band < WEIMIN,
!    minimisation will break after a maximum of two iterations
!  EBREAK  absolut break condition
!    intra-band minimisation is stopped if DE is < EBREAK
!  DEPER   intra-band break condition (see below)
!  ICOUEV  number of intraband evalue minimisations
!  DESUM   total change of the variational quantity 
!  RMS     norm of residual vector
!  LRESET  reset the wavefunction array entirely
!
!***********************************************************************

  SUBROUTINE LINEAR_RESPONSE_DIIS(GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WXI,W0,WDES, &
       LMDIM,CDIJ,CQIJ, RMS,DESUM,ICOUEV, SV, CSHIFT, IU6, IU0, LRESET, IERROR, FBREAK_MEAN_IN)
#ifdef _OPENACC
    USE mopenacc
#endif
    USE prec

    USE wave
    USE wave_high
    USE base
    USE lattice
    USE mpimy
    USE mgrid

    USE nonl_high
    USE hamil
    USE constant
    USE wave_mpi
    USE string, ONLY: str
    USE tutor, ONLY: vtutor

    IMPLICIT NONE

    TYPE (grid_3d)     GRID
    TYPE (info_struct) INFO
    TYPE (latt)        LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W             ! LR of orbitals   ( H(0) - e(0) S(0) ) |phi(1)> = - |xi>
    TYPE (wavespin)    WXI           ! |xi>
    TYPE (wavespin)    W0            ! original, unpeturbed orbitals
    TYPE (wavedes)     WDES

    RGRID   SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ) ! local potential
    INTEGER LMDIM  
    OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    OVERLAP CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    REAL(q) DESUM                    ! total change of e(1) related to phi(1)
    REAL(q) RMS                      ! magnitude of the residual vector
    INTEGER ICOUEV                   ! number of H | phi> evaluations
    REAL(q) CSHIFT                   ! complex shift
    INTEGER IU0, IU6                 ! units for output
    LOGICAL LRESET                   ! reset W0
    INTEGER, OPTIONAL :: IERROR      ! return error code
    REAL(q), OPTIONAL :: FBREAK_MEAN_IN ! in: mean FBREAK stopping codition in previous iteration
                                     ! out: mean FBREAK codition in current iteration
    !----- local work arrays
    TYPE (wavedes1) WDES1            ! descriptor for one k-point
    TYPE (wavefun1), TARGET :: W1(WDES%NSIM) ! current wavefunction
    TYPE (wavefun1), TARGET :: WTMP(WDES%NSIM)
    TYPE (wavefuna), TARGET :: WOPT,WTMPA

    REAL(q), ALLOCATABLE, TARGET :: PRECON(:,:)
    GDEF, ALLOCATABLE :: CWORK1(:,:),CHAM(:,:),B(:),CHAM_(:,:),B_(:)
    INTEGER, ALLOCATABLE :: IPIV(:)

    INTEGER :: NSIM,NRES             ! number of bands treated simultaneously
    INTEGER :: LD                    ! leading dimension of CF array
    INTEGER :: NITER                 ! maximum iteration count
    INTEGER :: NODE_ME, IONODE
    INTEGER :: NB(WDES%NSIM)         ! contains a list of bands currently optimized
    REAL(q) :: EVALUE0(WDES%NSIM)    ! eigenvalue e(0)
    COMPLEX(q) :: EVALUE0_C(WDES%NSIM)! version for complex shift
    REAL(q) :: FBREAK(WDES%NSIM)     ! relative break criterion for that band
    INTEGER :: IT(WDES%NSIM)         ! current iteration for this band

    REAL(q) :: FNORM(WDES%NSIM)      ! norm of residual vector for each band
    REAL(q) :: ORTH(WDES%NSIM)       ! orthogonality condition for each band
    REAL(q) :: EVAR(WDES%NSIM)       ! variational quantity for each band
    REAL(q) :: EVARP(WDES%NSIM)      ! previous variational quantity for each band

    REAL(q) :: FNORM_, ORTH_, EVAR_, EVALUE0_

    REAL(q) :: SLOCAL                ! average local potential
    REAL(q) :: DE_ATT                ! 1/4 of the total bandwidth
    REAL(q) :: EKIN                  ! kinetic energy
    REAL(q) :: TRIAL                 ! trial step
    REAL(q) :: OTRIAL                ! optimal trial step
    LOGICAL :: LSTOP,LNEWB
    INTEGER :: NP, ISP, NK, NB_DONE, N, IDUMP, ISPINOR, NPRO, M, MM
    INTEGER :: I, ITER, IFAIL, N1, N2
    REAL(q) :: X, X2
    REAL(q) :: ESTART
    COMPLEX(q) :: C
    REAL(q) :: FBREAK_MEAN
    INTEGER :: ICOU_FBREAK
!$  INTEGER __omp_id
!$  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM

    PROFILING_START('linear_response_diis')

#ifdef _OPENACC_TST
    PUSH_ACC_EXEC_ON(.TRUE.)

    ACC_ASYNC_Q=ACC_ASYNC_ASYNC

    CALL ACC_COPYIN_TYPED_VAR(WDES)
    CALL ACC_COPYIN_TYPED_VAR(W)
    CALL ACC_COPYIN_TYPED_VAR(WXI)
    CALL ACC_COPYIN_TYPED_VAR(W0)
    CALL ACC_COPYIN_TYPED_VAR(NONL_S)
    CALL ACC_COPYIN_TYPED_VAR(NONLR_S)
    CALL ACC_COPYIN_TYPED_VAR(GRID)

!$ACC ENTER DATA COPYIN(SV,CDIJ,CQIJ) __IF_ASYNC__
#endif

    INFO%IALGO=8

    NODE_ME=0
    IONODE =0
#ifdef MPI
    NODE_ME=WDES%COMM%NODE_ME
    IONODE =WDES%COMM%IONODE
#endif
!=======================================================================
!  INITIALISATION:
! maximum  number of iterations
!=======================================================================
    IF (PRESENT(IERROR)) THEN
       IERROR=0
    ENDIF
    IF (PRESENT(FBREAK_MEAN_IN)) THEN
       FBREAK_MEAN=0
       ICOU_FBREAK=0
       IF (FBREAK_MEAN_IN==0.0_q) FBREAK_MEAN_IN=1E10_q
    ENDIF

    NSIM=WDES%NSIM
    ! at least 6 iterations are required for save convergence
    ! since there is no other backup algorithm, safety first
    NITER=MAX(INFO%NDAV,6)
    NRES =NITER

    RMS   =0
    DESUM =0
    ESTART=0
    ICOUEV=0
    SLOCAL=MINVAL(REAL(W0%CELTOT(1,1:WDES%NKPTS,1:WDES%ISPIN),q))

    TRIAL = 0.3_q

    ALLOCATE(PRECON(WDES%NRPLWV,NSIM),CWORK1(NRES,NRES),CHAM(NRES,NRES),B(NRES),IPIV(NRES),CHAM_(NRES,NRES),B_(NRES))

    LD=WDES%NRPLWV*NRES*2

!$ACC ENTER DATA CREATE(WDES1) __IF_ASYNC__
    CALL SETWDES(WDES,WDES1,0)
!$ACC ENTER DATA CREATE(WTMPA,WOPT) __IF_ASYNC__
    CALL NEWWAVA(WOPT, WDES1, NRES*2, NSIM)
    CALL NEWWAVA(WTMPA, WDES1, NSIM)
!$ACC ENTER DATA CREATE(WTMP(:),W1(:)) __IF_ASYNC__
    DO NP=1,NSIM
       CALL NEWWAV_R(W1(NP),WDES1)
    ENDDO

    DO NP=1,SIZE(WTMP)
       CALL NEWWAV(WTMP(NP),WDES1,.FALSE.)
    ENDDO

!$ACC ENTER DATA CREATE(PRECON,FNORM,ORTH,EVAR,CHAM,B,CWORK1) __IF_ASYNC__
!!!$ACC ENTER DATA CREATE(PRECON,FNORM,ORTH,EVAR) __IF_ASYNC__

!=======================================================================
    spin:    DO ISP=1,WDES%ISPIN
    kpoints: DO NK=1,WDES%NKPTS
#ifdef MPI
    IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
    PROFILING_START('kpoints')

    CALL SETWDES(WDES,WDES1,NK)
!=======================================================================
    DE_ATT=ABS(W0%CELTOT(WDES%NB_TOT,NK,ISP)-W0%CELTOT(1,NK,ISP))/2

    IF (INFO%LREAL) THEN
       CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES)
    ELSE
       CALL PHASE(WDES,NONL_S,NK)
    ENDIF

    NB=0          ! empty the list of bands, which are optimized currently
    NB_DONE=0     ! index the bands already optimised
!=======================================================================
    bands: DO
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
       !
       !  check the NB list, whether there is any empty slot
       !  fill in a not yet optimized wavefunction into the slot
       !
       IDUMP=0

       PROFILING_START('newband')

       newband: DO NP=1,NSIM
          LNEWB=.FALSE.

          IF (NB(NP)==0.AND.NB_DONE<WDES%NBANDS) THEN
             NB_DONE=NB_DONE+1; N=NB_DONE; LNEWB=.TRUE.
          ENDIF

          IF (LNEWB) THEN

             SET_ACC_ASYNC_Q(NP)

             NB(NP)    =N
             FBREAK(NP)=0
             IT(NP)    =0
#ifndef _OPENMP
#ifdef debug
             IDUMP=2
#endif
#ifdef MPI
             IF (NODE_ME /= IONODE) IDUMP=0
#endif
             IF (IDUMP>=2) WRITE(*,*)
             IF (IDUMP>=2) WRITE(*,'(I3,1X)',ADVANCE='NO') N
#endif
             ! copy eigen energy from CELEN
             EVALUE0(NP) =W0%CELEN(N,NK,ISP)
             EVALUE0_C(NP)=EVALUE0(NP) +CMPLX(0.0_q,2.0_q*CSHIFT,q)

             EVALUE0_=EVALUE0(NP)-SLOCAL
             ! calculate the preconditioning matrix
             CALL SETUP_PRECOND( ELEMENT(W0, WDES1, N, ISP), 8,  IDUMP, PRECON(1,NP), EVALUE0_, DE_ATT )

             IF (LRESET) CALL APPLY_PRECOND( ELEMENT( WXI, WDES1, N, ISP), ELEMENT( W, WDES1, N, ISP), PRECON(1,NP), -1.0_q)

             CALL SETWAV(W,W1(NP),WDES1,N,ISP)
!-----------------------------------------------------------------------
! FFT of the current trial wave function
!-----------------------------------------------------------------------
             CALL FFTWAV_W1(W1(NP))
             IF (LRESET) CALL W1_PROJ(W1(NP),NONLR_S,NONL_S)
          ENDIF
       ENDDO newband

       PROFILING_STOP('newband')

!=======================================================================
! if the NB list is now empty end the bands DO loop
!=======================================================================
       LSTOP=.TRUE.
       W1%LDO  =.FALSE.
       DO NP=1,NSIM
          IF ( NB(NP) /= 0 ) THEN
             LSTOP  =.FALSE.
             W1(NP)%LDO=.TRUE.  ! band not finished yet
             IT(NP) =IT(NP)+1   ! increase iteration count
          ENDIF
       ENDDO
       IF (LSTOP) THEN
          IF (IDUMP>=2) WRITE(*,*)
          EXIT bands
       ENDIF

!=======================================================================
! intra-band minimisation
!=======================================================================
!-----------------------------------------------------------------------
! calculate the vector (H(0)-e(0) S(0)) |phi(1)_opt >
!-----------------------------------------------------------------------
       !  residual vector temporarily in WTMPA
       !  to have uniform stride for result array
       IF (CSHIFT==0) THEN
          CALL HAMILTMU(WDES1, W1, NONLR_S, NONL_S, EVALUE0, &
               &     CDIJ, CQIJ, SV, ISP,  WTMPA)
       ELSE
          CALL HAMILTMU_C(WDES1, W1, NONLR_S, NONL_S, EVALUE0_C, &
               &     CDIJ, CQIJ, SV, ISP, WTMPA)
       ENDIF
#ifdef _OPENACC
! NOTE(mM): this should be absorbed into HAMILTMU(_C) in the near future
       ACC_ASYNC_Q=ACC_ASYNC_ASYNC
#endif

       PROFILING_START('i2')

       DO NP=1,NSIM
          IF (.NOT. W1(NP)%LDO) CYCLE

          SET_ACC_ASYNC_Q(NP)

          CALL TRUNCATE_HIGH_FREQUENCY_W1( ELEMENT( WTMPA, NP), .FALSE., INFO%ENINI)

          PROFILING_START('gi2')

#ifdef _OPENACC
          CALL KERNEL1(NP,FNORM(NP),ORTH(NP),EVAR(NP))
#else
          N=NB(NP); ITER=IT(NP)

          FNORM_ =0
          ORTH_  =0
          EVAR_  =0

!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(MM,C) REDUCTION(+:FNORM_,ORTH_,EVAR_)
          DO ISPINOR=0,WDES%NRSPINORS-1
             DO M=1,WDES1%NGVECTOR
                MM=M+ISPINOR*WDES1%NGVECTOR

                !  |R> = H(0)-epsilon S(0) |phi(1)> + | xi >
                C=WTMPA%CW(MM,NP)+WXI%CW(MM,N,NK,ISP)
                !   <R|R>
                FNORM_ =FNORM_+C*CONJG(C)
                !   <phi(0)| H(0)-e(0) S(0) |phi(1)> +  <phi(0)| xi >
                !   since xi is orthogonal to phi(0), and <phi(0)| H(0)-e(0) S(0) 
                !   is zero as well, ORTH_ should be zero
                ORTH_  =ORTH_+C*CONJG(W0%CW(MM,N,NK,ISP))
                !   variational quantity
                !   <phi(1)|xi> + c.c + <phi(1)| H(0)-e(0) S(0)|phi(1)>
                EVAR_  =EVAR_+2*W%CW(MM,N,NK,ISP)*CONJG(WXI%CW(MM,N,NK,ISP)) & 
                     +W%CW(MM,N,NK,ISP)*CONJG(WTMPA%CW(MM,NP))
             ENDDO
          ENDDO
!$OMP END PARALLEL DO

          CALLMPI( M_sum_3(WDES%COMM_INB, FNORM_, ORTH_, EVAR_))

          FNORM(NP)=FNORM_
          ORTH(NP) =ORTH_
          EVAR(NP) =EVAR_
#endif
          PROFILING_STOP('gi2')

#ifdef _OPENACC
!=======================================================================
!=======================================================================
! BEGIN: split loop in OpenACC version for performance reasons

       ENDDO

       DO NP=1,NSIM
          IF (.NOT. W1(NP)%LDO) CYCLE

          N=NB(NP); ITER=IT(NP)

          SET_ACC_ASYNC_Q(NP)
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)

! END: split loop in OpenACC version for performance reasons
!=======================================================================
!=======================================================================
#endif

#if !defined(_OPENMP) && !defined(_OPENACC)
          IF (IDUMP>=2) THEN
             WRITE(*,'(E9.2,"R")',ADVANCE='NO') SQRT(ABS(FNORM(NP)))
             WRITE(*,'(E9.2,"O")',ADVANCE='NO') ORTH(NP)
             WRITE(*,'(E9.2,"E")',ADVANCE='NO') EVAR(NP)
          ENDIF
#endif
          IF (ITER==1) THEN
             ! total norm of error vector at start
             RMS=RMS+WDES%RSPIN*WDES%WTKPT(NK)*W0%FERWE(N,NK,ISP)*SQRT(ABS(FNORM(NP)))/WDES%NB_TOT
             ESTART=ESTART+WDES%RSPIN*WDES%WTKPT(NK)*W0%FERWE(N,NK,ISP)* EVAR(NP)
          ELSE
             DESUM =DESUM +WDES%RSPIN*WDES%WTKPT(NK)*W0%FERWE(N,NK,ISP)*(EVAR(NP)-EVARP(NP))
          ENDIF
          ! store variational quantity
          EVARP(NP)=EVAR(NP)
          W%CELEN(N,NK,ISP)=EVAR(NP)

          ! norm of total error vector before start
          ! norm smaller than EBREAK stop  |e -e(app)| < | Residuum |
          IF (ABS(FNORM(NP))<INFO%EBREAK/10) THEN
#if !defined(_OPENMP) && !defined(_OPENACC)
             IF (IDUMP>=2) WRITE(*,'("X")',ADVANCE='NO')
#endif
             W1(NP)%LDO=.FALSE.
          ENDIF

          ! stop working on this band if ITER > NITER
          IF (ITER>NITER) W1(NP)%LDO=.FALSE.

          IF (.NOT. W1(NP)%LDO) CYCLE
!-----------------------------------------------------------------------
! fill current wavefunctions into work array WOPT%CW at position ITER
!-----------------------------------------------------------------------
          CALL W1_COPY(W1(NP), ELEMENT(WOPT, ITER, NP))
          CALL W1_COPY(ELEMENT(WTMPA, NP), ELEMENT(WOPT, NRES+ITER,NP))

          IF (ITER>1) THEN
             ! better conditioning for search
             ! w(iter-1)=w(iter)-w(iter-1)
             CALL W1_DSCAL( ELEMENT( WOPT, ITER-1, NP), -1.0_q)
             CALL W1_DAXPY( ELEMENT( WOPT, ITER, NP), 1.0_q, ELEMENT( WOPT, ITER-1, NP)) 

             ! gradient(iter-1)=gradient(iter)-gradient(iter-1)
             CALL W1_DSCAL( ELEMENT( WOPT, NRES+ITER-1, NP), -1.0_q)
             CALL W1_DAXPY( ELEMENT( WOPT, NRES+ITER, NP), 1.0_q, ELEMENT( WOPT, NRES+ITER-1, NP)) 
          ENDIF

#ifdef _OPENACC
!=======================================================================
!=======================================================================
! BEGIN: split loop in OpenACC version for performance reasons

       ENDDO

       DO NP=1,NSIM
          IF (.NOT. W1(NP)%LDO) CYCLE

          N=NB(NP); ITER=IT(NP)

          SET_ACC_ASYNC_Q(NP)

! END: split loop in OpenACC version for performance reasons
!=======================================================================
!=======================================================================
#endif

!***********************************************************************
! inverse interation step
! minimize
!    | ( H - e S) | phi(1) > + | xi > |^ 2  -> min
! in the yet available subspace spanned by the wavefunction stored in CF
! if one denotes these wavefunctions as phi(1)_j, and R_j=  (H - e S) phi(1)_j
! the following equation is obtained:
!  sum_ij  b_i* < R_i | R_j > b_j + sum_i b_i* <R_i | xi > + c.c. -> min
! or equivalently
!  sum_j  < R_i | R_j > b_j  = - <R_i | xi >
! the new optimized wavefunction is given by solving this linear
! equation for b
! the scalar product < | > can be evaluated with any metric
!***********************************************************************

!$ACC KERNELS PRESENT(CHAM,B) __IF_ASYNC__
!!!$ACC ENTER DATA CREATE(CHAM,B,CWORK1) __IF_ASYNC__
          CHAM=0 ; B=0
!$ACC END KERNELS

          !    A(n2,n1)=    < phi_n2 |  ( H - e S) ( H - e S)  | phi_n1 >
          builda: DO N1=1,ITER
             CALL W1_GEMV( one, ELEMENTS( WOPT, NRES+N1, NRES+ITER, NP),  ELEMENT( WOPT, NRES+N1, NP), zero, CWORK1(1,N1), 1)

!$ACC PARALLEL LOOP PRESENT(CHAM,CWORK1) __IF_ASYNC__
             DO N2=N1,ITER
                CHAM(N2,N1)=       GREAL(CWORK1(N2-N1+1,N1))
                CHAM(N1,N2)=GREAL(GCONJG(CWORK1(N2-N1+1,N1)))
             ENDDO
          ENDDO builda

          !     B(n1) =   - <R_n1 | xi >= - < phi_n1 | ( H - e S) | xi >
          CALL W1_GEMV( one, ELEMENTS( WOPT, NRES+1, NRES+ITER, NP), ELEMENT( WXI, WDES1, N, ISP), zero, B(1), 1)

!$ACC KERNELS PRESENT(B) __IF_ASYNC__
          B(1:ITER)=-GREAL(B(1:ITER))
!$ACC END KERNELS

!$ACC UPDATE SELF(CHAM,B) __IF_ASYNC__
!!!$ACC EXIT DATA COPYOUT(CHAM,B) DELETE(CWORK1) __IF_ASYNC__

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)

#if !defined(_OPENMP) && !defined(_OPENACC)
          CHAM_=CHAM
          B_ =B
#endif

          ! calculate the solution of sum_j CHAM(i,j) * X(j) = B(i)
          ! overwrite B by X
          PROFILING_START('ggetrf')
          CALL GGETRF( ITER, ITER, CHAM, NRES, IPIV, IFAIL )
          PROFILING_STOP('ggetrf')

          IF (IFAIL==0) THEN
             PROFILING_START('ggetrs')
             CALL GGETRS('N', ITER, 1, CHAM, NRES, IPIV, B, NRES, IFAIL)
             PROFILING_STOP('ggetrs')
          ENDIF

#if !defined(_OPENMP) && !defined(_OPENACC)
          IF (.FALSE.) THEN
             ! dump the matrix and the solution vector
             io_begin
             N2=MIN(10,ITER)
             WRITE(6,*)
             DO N1=1,N2
                WRITE(*,'("m",I3,8E14.7)')N1, CHAM_(N1,1:N2)
             ENDDO
             WRITE(*,'(A4,8E14.7)') 'b', B_(1:N2)

             WRITE(*,*)
             WRITE(*,'(A4,8E14.7)') 'e', B (1:N2)
             io_end
          ENDIF
#endif

          ! matrix singular and first iteration
          ! this usually means the trial vectors are zero; stop immedeatly
          IF (IFAIL/=0 .AND. ITER==1) THEN
             CALL vtutor%bug("LINEAR_RESPONSE_DIIS: matrix is zero, try to call with LRESET in the first iteration " &
                              // str(N) // " " // str(NK),__FILE__,__LINE__)
          ELSEIF (IFAIL/=0) THEN
             IF (PRESENT(IERROR)) THEN
                IERROR=IERROR+1
             ELSE
                IF (IU6>=0) WRITE(IU6,219) IFAIL,ITER,N
                IF (IU0>=0) WRITE(IU0,219) IFAIL,ITER,N
 219            FORMAT('WARNING in EDDRMM_LR: call to GGETRF failed, returncode =',I4,I2,I2)
             ENDIF
             !  try to save things somehow, goto next band
#if !defined(_OPENMP) && !defined(_OPENACC)
             IF (IDUMP>=2) WRITE(*,'("Z")',ADVANCE='NO')
#endif
             W1(NP)%LDO=.FALSE.
             CYCLE
          ENDIF

#if !defined(_OPENMP) && !defined(_OPENACC)
          IF (ITER==2 .AND. IDUMP==2) THEN
             ! write out 'optimal trial step' i.e step which would have minimized
             ! the residuum
             IF (ITER==2) THEN
                OTRIAL= REAL( 1+B(1)/B(2) ,KIND=q)
                WRITE(*,'(1X,F7.4,"o")',ADVANCE='NO') OTRIAL
             ENDIF
          ENDIF

          IF (IDUMP >= 3) THEN
             ! set CWORK1(1) to < xi | xi >
             C=W1_DOT( ELEMENT(WXI, WDES1, N, ISP) , ELEMENT(WXI, WDES1, N, ISP))

             DO N1=1,ITER
                DO N2=1,ITER
                   C=C+GCONJG(B(N2))*CHAM_(N2,N1)*B(N1)
                ENDDO
                C=C-B_(N1)*GCONJG(B(N1))-GCONJG(B_(N1))*B(N1)
             ENDDO
             ! residual after the step
             WRITE(*,'(1X,E9.2,"rs")',ADVANCE='NO') SQRT(ABS(C))
          ENDIF
#endif

!=======================================================================
! now performe the trial step (default TRIAL)
!=======================================================================
          ! W1=0
          CALL W1_DSCAL( W1(NP), 0.0_q)

          ! W1=W1 + B(I,1) *WOPT(I,NP)
          DO I=1,ITER
             CALL W1_GAXPY( ELEMENT(WOPT, I,NP), B(I), W1(NP))
          ENDDO

          ICOUEV=ICOUEV+1

! check for break condition
          IF ( (ITER>1 .AND. ABS(FNORM(NP))<FBREAK(NP)) .OR. ITER == NITER) THEN
             IF (IDUMP>=2) WRITE(*,'("B")',ADVANCE='NO')
             W1(NP)%LDO=.FALSE.
          ELSE
             ! trial step on wavefunction moving from the yet optimised wavefunction
             ! along the residual vector for that wavefunction
             !      - Precond { sum_i b_i ( H(0) - e(0) S(0)) |phi(1)_i> + xi }
             ! this is somewhat dangerous in the very last step 
             DO I=1,ITER
                CALL APPLY_PRECOND( ELEMENT(WOPT, NRES+I, NP), WTMP(NP), PRECON(1,NP))
                CALL W1_GAXPY( WTMP(NP), (-TRIAL*B(I)), W1(NP))
             ENDDO
             ! - Precond xi
             CALL ADD_PRECOND( ELEMENT(WXI, WDES1, N, ISP), W1(NP), PRECON(1,NP), -TRIAL)
          END IF
          ! transform the wave-function to real space
          CALL FFTWAV_W1(W1(NP))

       ENDDO

#ifdef _OPENACC
       CALL ACC_SYNC_ASYNC_Q(NSIM,W1(:)%LDO)
       ACC_ASYNC_Q=ACC_ASYNC_ASYNC
#endif

       PROFILING_STOP('i2')

       ! calculate results of projection operatores
       CALL W1_PROJALL(WDES1, W1, NONLR_S, NONL_S, NSIM)
!=======================================================================
! break of intra-band-minimisation
!=======================================================================
       i3: DO NP=1,NSIM
          N=NB(NP); ITER=IT(NP); IF (.NOT. W1(NP)%LDO) CYCLE i3

          IF (ITER==1) THEN
             FBREAK(NP)=ABS(FNORM(NP))*INFO%DEPER
             IF (PRESENT(FBREAK_MEAN_IN)) THEN
                FBREAK(NP)=MIN(ABS(FNORM(NP))*INFO%DEPER, FBREAK_MEAN_IN)
                FBREAK_MEAN=FBREAK_MEAN+FBREAK(NP)
                ICOU_FBREAK=ICOU_FBREAK+1
             ENDIF
          ENDIF
       ENDDO i3
       
       ! one band just finished ?, set NB(NP) also to 0 and finish everything
       DO NP=1,NSIM
          N=NB(NP)
          IF (.NOT. W1(NP)%LDO .AND. N /=0 ) THEN
             NB(NP)=0
             IF (IDUMP==10) WRITE(*,*)
          ENDIF
       ENDDO
!=======================================================================
! move onto the next Band
!=======================================================================

    ENDDO bands

    PROFILING_STOP('kpoints')

!=======================================================================
    ENDDO kpoints
    ENDDO spin
!=======================================================================

    IF (PRESENT(IERROR)) THEN
       CALLMPI( M_sum_i(WDES%COMM_INTER, IERROR ,1))
       CALLMPI( M_sum_i(WDES%COMM_KINTER, IERROR ,1))
    ENDIF
    IF (PRESENT(FBREAK_MEAN_IN)) THEN
       CALLMPI( M_sum_d(WDES%COMM_INTER, FBREAK_MEAN ,1))
       CALLMPI( M_sum_d(WDES%COMM_KINTER, FBREAK_MEAN ,1))
       CALLMPI( M_sum_i(WDES%COMM_INTER, ICOU_FBREAK ,1))
       CALLMPI( M_sum_i(WDES%COMM_KINTER, ICOU_FBREAK ,1))
       FBREAK_MEAN_IN=FBREAK_MEAN/ICOU_FBREAK
    ENDIF
    CALLMPI( M_sum_d(WDES%COMM_INTER, RMS, 1))
    CALLMPI( M_sum_d(WDES%COMM_KINTER, RMS, 1))

    CALLMPI( M_sum_d(WDES%COMM_INTER, DESUM, 1))
    CALLMPI( M_sum_d(WDES%COMM_KINTER, DESUM, 1))

    CALLMPI( M_sum_i(WDES%COMM_INTER, ICOUEV ,1))
    CALLMPI( M_sum_i(WDES%COMM_KINTER, ICOUEV ,1))

    DO NP=1,NSIM
       CALL DELWAV_R(W1(NP))
    ENDDO

    DO NP=1,SIZE(WTMP)
       CALL DELWAV(WTMP(NP),.FALSE.)
    ENDDO
!$ACC EXIT DATA DELETE(WTMP(:),W1(:)) __IF_ASYNC__

    CALL DELWAVA(WOPT)
    CALL DELWAVA(WTMPA)
!$ACC EXIT DATA DELETE(WOPT,WTMPA) __IF_ASYNC__
!$ACC EXIT DATA DELETE(PRECON,FNORM,ORTH,EVAR,CHAM,B,CWORK1) __IF_ASYNC__

    DEALLOCATE(PRECON,CWORK1,CHAM,B,IPIV,CHAM_,B_)

#ifdef _OPENACC
    CALL ACC_DELETE_TYPED_VAR(WDES1)
!!!$ACC WAIT IF(ACC_EXEC_ON)
!!    ACC_ASYNC_Q=ACC_ASYNC_SYNC
#endif

    ! WRITE(*,*) 'start energy',ESTART

#ifdef _OPENACC_TST
!$ACC UPDATE DEVICE(W%CELEN) __IF_ASYNC__
!$ACC EXIT DATA DELETE(SV,CDIJ,CQIJ) __IF_ASYNC__

    CALL ACC_DELETE_TYPED_VAR(GRID)
    CALL ACC_DELETE_TYPED_VAR(NONLR_S)
    CALL ACC_DELETE_TYPED_VAR(NONL_S)
    CALL ACC_DELETE_TYPED_VAR(W0)
    CALL ACC_DELETE_TYPED_VAR(WXI)
    CALL ACC_COPYOUT_TYPED_VAR(W)
    CALL ACC_DELETE_TYPED_VAR(WDES)

    POP_ACC_EXEC_ON
#endif

    PROFILING_STOP('linear_response_diis')
#ifdef _OPENACC
    CONTAINS

       SUBROUTINE KERNEL1(NP,FNORM,ORTH,EVAR)
       REAL(q) :: FNORM,ORTH,EVAR
       INTEGER :: NP
       ! local
       COMPLEX(q) :: C
       INTEGER :: N,ISPINOR,M,MM

       N=NB(NP)

!$ACC KERNELS PRESENT(FNORM,ORTH,EVAR) __IF_ASYNC__
       FNORM =0
       ORTH  =0
       EVAR  =0
!$ACC END KERNELS

!$ACC PARALLEL LOOP PRESENT(WDES1,WTMPA,WXI,W0,W,FNORM,ORTH,EVAR) COLLAPSE(2) PRIVATE(MM,C) REDUCTION(+:FNORM,ORTH,EVAR) __IF_ASYNC__
       DO ISPINOR=0,WDES%NRSPINORS-1
          DO M=1,WDES1%NGVECTOR
             MM=M+ISPINOR*WDES1%NGVECTOR
             !  |R> = H(0)-epsilon S(0) |phi(1)> + | xi >
             C=WTMPA%CW(MM,NP)+WXI%CW(MM,N,NK,ISP)
             !   <R|R>
             FNORM =FNORM+C*CONJG(C)
             !   <phi(0)| H(0)-e(0) S(0) |phi(1)> +  <phi(0)| xi >
             !   since xi is orthogonal to phi(0), and <phi(0)| H(0)-e(0) S(0)
             !   is zero as well, ORTH should be zero
             ORTH  =ORTH+C*CONJG(W0%CW(MM,N,NK,ISP))
             !   variational quantity
             !   <phi(1)|xi> + c.c + <phi(1)| H(0)-e(0) S(0)|phi(1)>
             EVAR  =EVAR+2*W%CW(MM,N,NK,ISP)*CONJG(WXI%CW(MM,N,NK,ISP)) &
                  +W%CW(MM,N,NK,ISP)*CONJG(WTMPA%CW(MM,NP))
          ENDDO
       ENDDO

       CALLMPI( M_sum_d(WDES%COMM_INB, FNORM, 1))
       CALLMPI( M_sum_d(WDES%COMM_INB, ORTH , 1))
       CALLMPI( M_sum_d(WDES%COMM_INB, EVAR , 1))
!$ACC UPDATE SELF(FNORM,ORTH,EVAR) __IF_ASYNC__
       RETURN
       END SUBROUTINE KERNEL1
#endif
  END SUBROUTINE LINEAR_RESPONSE_DIIS
END MODULE RMM_DIIS_LR
