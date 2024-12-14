!#define debug
#include "symbol.inc"
MODULE rmm_diis
  USE prec
  USE hamil
CONTAINS
!************************ SUBROUTINE EDDRMM *****************************
! RCS:  $Id: rmm-diis.F,v 1.7 2002/08/14 13:59:42 kresse Exp $
!
! this subroutine performes an optimization of the trial wavefunctions
! minimizing the expectation value of the Hamiltonian  i.e.
!     < phi | H |  phi >
! or the norm of the residual vector  i.e.
!    r^2 = < phi | H -e S | H - e S | phi >
! or using an inverse iteration method. In the last case
!    || ( H - e_initial S)  | phi > - | phi_initial > ||
! is optimized.
! The full name of the residual vector  minimization method
! is residual vector minimiziation method-
! direct inversion of the iterative subspace (RMM-DIIS)
!    see: D. M. Wood and A. Zunger, J. Phys. A, 1343 (1985)
!    and  P. Pulay,  Chem. Phys. Lett. 73, 393 (1980).
!
!
!  INFO%IALGO   determine type of preconditioning and the algorithm
!    0    inverse iteration         +  TAP preconditioning
!    6    rms-minimization          +  TAP preconditioning
!    7    rms-minimization          +  no preconditioning
!    8    precond rms-minimization  +  TAP preconditioning
!    9    precond rms-minimization  +  Jacobi like preconditioning
!    (TAP Teter Alan Payne)
!  parameters:
!  LDELAY=.TRUE.
!          steepest descent eigenvalue minimization
!          maximum number of steps is 2
!  WEIMIN  treshhold for total energy minimisation
!    is the fermiweight of a band < WEIMIN,
!    minimisation will break after a maximum of two iterations
!  EBREAK  absolut break condition
!    intra-band minimisation is stopped if DE is < EBREAK
!  DEPER   intra-band break condition (see below)
!  ICOUEV  number of intraband evalue minimisations
!  DESUM   total change in eigenvalues
!  RMS     norm of residual vector
!
!> @details @ref openmp :
!
!> @details @ref openacc :
!
!***********************************************************************

  SUBROUTINE EDDRMM(HAMILTONIAN,GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WDES, &
       LMDIM,CDIJ,CQIJ, RMS,DESUM,ICOUEV, SV,IU6,IU0, LDELAY, FBREAK_MEAN_IN)
#ifdef _OPENACC
    USE mopenacc
#endif
    USE prec
    USE base
    USE constant
    USE lattice
    USE mpimy
    USE mgrid
    USE nonl_high
    USE hamil_struct_def
    USE wave_high
    USE wave_mpi

    USE openmp, ONLY : omp_nthreads_acc

    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (ham_handle)  HAMILTONIAN
    TYPE (grid_3d)     GRID
    TYPE (info_struct) INFO
    TYPE (latt)        LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W
    TYPE (wavedes)     WDES

    LOGICAL LDELAY
    RGRID   SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ) ! local potential
    OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    REAL(q), OPTIONAL :: FBREAK_MEAN_IN ! in:  mean FBREAK stopping condition in previous iteration
                                        ! out: mean FBREAK condition in current iteration

  ! local variables
    TYPE (wavedes1)            WDES1          ! descriptor for one k-point
    TYPE (wavefun1), TARGET :: W1(WDES%NSIM)  ! current wavefunction
    TYPE (wavefun1), TARGET :: WTMP(WDES%NSIM)! temporary trial wavefunction

    REAL(q) R(INFO%NDAV)
    PARAMETER (LWORK=20)
    DIMENSION CWORK(LWORK*INFO%NDAV)
    DIMENSION RWORK(3*INFO%NDAV)

    TYPE(wavefuna), TARGET :: W_INI, WOPT, WTMPA
    REAL(q),   ALLOCATABLE :: PRECON(:,:)

    GDEF, ALLOCATABLE    :: CHAM(:,:),CTMP(:,:),CWORK1(:)
    GDEF, ALLOCATABLE    :: B(:),B_(:)
    INTEGER, ALLOCATABLE :: IPIV(:)

    INTEGER :: NB(WDES%NSIM)         ! contains a list of bands currently optimized
    REAL(q) :: EVALUE_INI(WDES%NSIM) ! eigenvalue of that band at the beginning
    REAL(q) :: EVALUE(WDES%NSIM)     ! eigenvalue during optimization
    REAL(q) :: DEIT(WDES%NSIM)       ! relative break criterion for that band
    REAL(q) :: IT(WDES%NSIM)         ! current iteration for this band
    REAL(q) :: FPRE(WDES%NSIM)       ! norm of residual vector for each band
    REAL(q) :: TRIAL(WDES%NSIM)      ! trial step for each band
    LOGICAL :: LNEWB,LSTOP
    LOGICAL :: LABORT(WDES%NSIM)     ! abort iteration on this band
    REAL(q) :: FBREAK_MEAN           ! required minimum error in the residual vector
    INTEGER :: ICOU_FBREAK
    TYPE (REDIS_PW_CTR), POINTER :: H_PW

    PROFILING_START('eddrmm')
#ifdef _OPENACC
    ACC_ASYNC_Q=ACC_ASYNC_ASYNC
#endif
    NODE_ME=0
    IONODE =0
#ifdef MPI
    NODE_ME=WDES%COMM%NODE_ME
    IONODE =WDES%COMM%IONODE
#endif
!=======================================================================
!  INITIALISATION:
! maximum  number of iterations
! NRES position where H - E S| trial vector > is stored
!=======================================================================
    IF (PRESENT(FBREAK_MEAN_IN)) THEN
       FBREAK_MEAN=0
       ICOU_FBREAK=0
    ! no minimum error required by calling routine
       IF (FBREAK_MEAN_IN==0.0) FBREAK_MEAN_IN=1E10
    ENDIF
    NSIM=WDES%NSIM
    NITER=INFO%NDAV
    IF (LDELAY) NITER=MIN(NITER,1)
    IF (LDELAY .AND. INFO%IALGO ==0) NITER=1
    NRES =INFO%NDAV

    DESUM =0
    RMS   =0
    ICOUEV=0

    SLOCAL=0
    DO I=1,GRID%RL%NP
       SLOCAL=SLOCAL+SV(I,1)
    ENDDO

    CALLMPI( M_sum_d(WDES%COMM_INB, SLOCAL, 1))
    SLOCAL=SLOCAL/GRID%NPLWV

    ALLOCATE(PRECON(WDES%NRPLWV,NSIM), &
   &         CHAM(NRES,NRES),CTMP(NRES,NRES),CWORK1(NRES), &
   &         B(NRES),B_(NRES),IPIV(NRES))

!$ACC ENTER DATA CREATE(WDES1) __IF_ASYNC__
    CALL SETWDES(WDES,WDES1,0)

!$ACC ENTER DATA CREATE(W_INI,WTMPA,WOPT) __IF_ASYNC__
    CALL NEWWAVA(W_INI, WDES1, NSIM)
    CALL NEWWAVA(WTMPA, WDES1, NSIM)
    CALL NEWWAVA(WOPT,  WDES1, NRES*2, NSIM)

!$ACC ENTER DATA CREATE(WTMP(:),W1(:)) __IF_ASYNC__
    DO NP=1,NSIM
       CALL NEWWAV(WTMP(NP), WDES1, .TRUE.)
       CALL NEWWAV_R(W1(NP), WDES1)
    ENDDO

    CTMP=0
    CHAM=0

!=======================================================================
! do we have to distribute the wavefunctions back ?
!=======================================================================
    IF (W%OVER_BAND) THEN
#ifdef MPI
       NCPU=WDES%COMM_INTER%NCPU ! number of procs involved in band dis.
#else
       NCPU=1
#endif
       NSTRIP=MIN(NSIM,WDES%NBANDS)*2
       CALL REDIS_PW_ALLOC(WDES, NSTRIP, H_PW)
    ENDIF

#ifdef _OPENACC
!$ACC ENTER DATA COPYIN(SV,CDIJ,CQIJ) CREATE(PRECON) __IF_ASYNC__
    CALL ACC_COPYIN_TYPED_VAR(GRID)
    CALL ACC_COPYIN_TYPED_VAR(NONLR_S)
    CALL ACC_COPYIN_TYPED_VAR(NONL_S)
#endif
!=======================================================================
    spin:    DO ISP=1,WDES%ISPIN
    kpoints: DO NK=1,WDES%NKPTS
#ifdef MPI
    IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
    CALL SETWDES(WDES,WDES1,NK)
!=======================================================================
    !  first initiate communication between bands
    IF (W%OVER_BAND) THEN
       DO N=1,NSTRIP
          CALL REDIS_PW_START(WDES, W%CW(1,N,NK,ISP), N, H_PW)
       ENDDO
    ENDIF
!=======================================================================
    DE_ATT=ABS(W%CELTOT(WDES%NB_TOTK(NK,ISP),NK,ISP)-W%CELTOT(1,NK,ISP))/4

    IF (INFO%LREAL) THEN
       CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES)
    ELSE
       CALL PHASE(WDES,NONL_S,NK)
    ENDIF

    NB=0          ! empty the list of bands, which are optimized currently
    NB_DONE=0     ! index the bands allready optimised
    IT=0
!=======================================================================
    bands: DO
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
       !
       !  check the NB list, whether there is any empty slot
       !  fill in a not yet optimized wavefunction into the slot
       !
       IDUMP=0

       PROFILING_START('newband')

DOACC !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(NP,N,IPREC,LNEWB) SCHEDULE(dynamic) &
DOACC !$OMP NUM_THREADS(omp_nthreads_acc)
       newband: DO NP=1,NSIM
          LNEWB=.FALSE.
DOACC !$OMP CRITICAL
          IF (NB(NP)==0.AND.NB_DONE<WDES%NBANDS) THEN
             NB_DONE=NB_DONE+1; N=NB_DONE; LNEWB=.TRUE.
          ENDIF
DOACC !$OMP END CRITICAL
          IF (LNEWB) THEN
             NB(NP)    =N
             LABORT(NP)=.FALSE.

             SET_ACC_ASYNC_Q(NP)

             IF (W%OVER_BAND) THEN
                CALL REDIS_PW_STOP (WDES, W%CW(1,N,NK,ISP), N, H_PW)
                IF (N+NSTRIP<=WDES%NBANDS) &
                     CALL REDIS_PW_START(WDES, W%CW(1,N+NSTRIP,NK,ISP), N+NSTRIP, H_PW)
             ENDIF

             CALL SETWAV(W,W1(NP),WDES1,N,ISP)  ! fill band N into W1(NP)

             CALL W1_COPY(W1(NP), ELEMENT(W_INI, NP))
#ifndef _OPENMP
#ifdef debug
             IDUMP=2
#endif
#ifdef MPI
             IF (NODE_ME /= IONODE) IDUMP=0
#endif
             IF (IDUMP==2) WRITE(*,'(I3,1X)',ADVANCE='NO') N
#endif
             ! start with FFT and the exact evaluation of the eigenenergy
             CALL FFTWAV_W1(W1(NP))

             IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
                CALL ECCP_TAU(WDES1,W1(NP),W1(NP),LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP),LATT_CUR,HAMILTONIAN%MU(:,ISP),W%CELEN(N,NK,ISP))
             ELSE
                CALL ECCP(WDES1,W1(NP),W1(NP),LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP), W%CELEN(N,NK,ISP))
             ENDIF
             EVALUE_INI(NP)=W%CELEN(N,NK,ISP)

             IF (IDUMP==2) WRITE(*,'(F9.4,"E")',ADVANCE='NO') REAL( W%CELEN(N,NK,ISP) ,KIND=q)

             ! calculate the preconditioning matrix
             CALL TRUNCATE_HIGH_FREQUENCY_W1( W1(NP), LDELAY, INFO%ENINI)
             IPREC=INFO%IALGO
             IF (LDELAY) IPREC=8

             CALL SETUP_PRECOND( W1(NP), IPREC, IDUMP, PRECON(1,NP), EVALUE_INI(NP)-SLOCAL, DE_ATT )

             DEIT(NP)=0
             IT(NP)  =0
          ENDIF
       ENDDO newband
DOACC !$OMP END PARALLEL DO

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
       IF (LSTOP) EXIT bands
!=======================================================================
! intra-band minimisation
!=======================================================================

       PROFILING_START('i1')

DOACC !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(NP,N,ITER) SCHEDULE(dynamic) &
DOACC !$OMP NUM_THREADS(omp_nthreads_acc)
       i1: DO NP=1,NSIM
          N=NB(NP); ITER=IT(NP); IF (.NOT. W1(NP)%LDO) CYCLE i1

          SET_ACC_ASYNC_Q(NP)

          ! fill current wavefunctions into work arrays WOPT at position (ITER,NP)
          CALL W1_COPY(W1(NP), ELEMENT(WOPT, ITER, NP))
          EVALUE(NP)=W%CELEN(N,NK,ISP)
       ENDDO i1
DOACC !$OMP END PARALLEL DO

       PROFILING_STOP('i1')

       !  store H-epsilon Q | psi > temporarily in upmost storage positions (2*NRES)
       !  to have uniform stride for result array see ELEMENTS_SECOND
       IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
          CALL HAMILTMU_TAU(WDES1, W1, NONLR_S, NONL_S, EVALUE_INI, &
         &     CDIJ, CQIJ, SV, LATT_CUR, HAMILTONIAN%MU, ISP, WTMPA )
       ELSE
          CALL HAMILTMU(WDES1,W1,NONLR_S,NONL_S,EVALUE_INI, &
         &     CDIJ,CQIJ, SV, ISP, WTMPA)
       ENDIF

       PROFILING_START('i2')

DOACC !$OMP PARALLEL DO NUM_THREADS(omp_nthreads_acc) DEFAULT(SHARED) SCHEDULE(dynamic) &
DOACC !$OMP PRIVATE(NP,N,ITER,FNORM,CTMP,CHAM,N1,CWORK1,N2,R,CWORK,RWORK,IFAIL,OTRIAL,B,B_,IPIV,C,I) &
DOACC !$OMP REDUCTION(+:RMS,FBREAK_MEAN,ICOU_FBREAK)
       i2: DO NP=1,NSIM
          N=NB(NP); ITER=IT(NP); IF (.NOT. W1(NP)%LDO) CYCLE i2

          SET_ACC_ASYNC_Q(NP)

          ! copy to proper storage position (NRES+ITER, NP)
          CALL W1_COPY(ELEMENT(WTMPA, NP), ELEMENT(WOPT, NRES+ITER, NP) )

          CALL TRUNCATE_HIGH_FREQUENCY_W1( ELEMENT( WOPT, NRES+ITER, NP), LDELAY, INFO%ENINI)
          CALL PW_NORM_WITH_METRIC_W1(ELEMENT( WOPT, NRES+ITER, NP), FNORM, FPRE(NP), PRECON(1,NP))

          IF (INFO%IALGO==6) FPRE(NP)=FNORM

          IF (IDUMP==2) WRITE(*,'(E9.2,"R")',ADVANCE='NO') SQRT(ABS(FNORM))

          IF (ITER==1) RMS=RMS+WDES%RSPIN*WDES%WTKPT(NK)*W%FERWE(N,NK,ISP)*SQRT(ABS(FNORM))/WDES%NB_TOT*WDES%NRSPINORS

          ! norm of total error vector before start
          ! norm smaller than EBREAK stop |e -e(app)| < | Residuum |
          IF (ABS(FNORM)<INFO%EBREAK) THEN
             W1(NP)%LDO=.FALSE.
             CYCLE i2
          ENDIF

          !----------------------------------------------------------------------
          optsubspace: IF (.NOT. LDELAY .AND. ITER > 1) THEN

             ! better conditioning for search
             ! w(iter-1)=w(iter)-w(iter-1)
             CALL W1_DSCAL( ELEMENT( WOPT, ITER-1, NP), -1.0_q)
             CALL W1_DAXPY( ELEMENT( WOPT, ITER, NP), 1.0_q, ELEMENT( WOPT, ITER-1, NP)) 

             ! gradient(iter-1)=gradient(iter)-gradient(iter-1)
             CALL W1_DSCAL( ELEMENT( WOPT, NRES+ITER-1, NP), -1.0_q)
             CALL W1_DAXPY( ELEMENT( WOPT, NRES+ITER, NP), 1.0_q, ELEMENT( WOPT, NRES+ITER-1, NP)) 

!***********************************************************************
! RMM-DIIS step
!
! minimize norm of residual vector in the subspace spanned by
! the set of wavefunctions stored in WOPT
! residual vectors are stored in WOPT starting at NRES
!
!***********************************************************************
             optsub: IF (INFO%IALGO /=0) THEN
!----------------------------------------------------------------------
! calculate  matrix INFO%IALGO=7 or 8
! CHAM(n2,n1)= <R(n2)| preconditioning |R(n1)>
! CTMP(n2,n1)= <phi(n2)| S |phi(n1)>
!----------------------------------------------------------------------
!$ACC ENTER DATA CREATE(CTMP,CHAM,CWORK1) __IF_ASYNC__
!$ACC KERNELS PRESENT(CTMP,CHAM) __IF_ASYNC__
                CTMP(1:NRES,1:NRES)=0
                CHAM(1:NRES,1:NRES)=0
!$ACC END KERNELS
                buildh: DO N1=1,ITER
                   IF (INFO%IALGO==8 .OR. INFO%IALGO==9) THEN
!                   IF (INFO%IALGO==8) THEN
                      CALL APPLY_PRECOND( ELEMENT(WOPT, NRES+N1,NP), WTMP(NP), PRECON(1,NP))
                   ELSE
                      CALL W1_COPY( ELEMENT(WOPT, NRES+N1,NP), WTMP(NP))
                   ENDIF

                   ! elements WOPT(NRES+N1: NRES+ITER, NP)
                   CALL W1_GEMV( one, ELEMENTS( WOPT, NRES+N1, NRES+ITER, NP),  WTMP(NP), zero, CWORK1, 1)

!$ACC PARALLEL LOOP PRESENT(CHAM,CWORK1) __IF_ASYNC__
                   DO N2=N1,ITER
                      CHAM(N2,N1)=      GREAL(CWORK1(N2-N1+1))
                      CHAM(N1,N2)=GREAL(GCONJG(CWORK1(N2-N1+1)))
                   ENDDO

                   ! elements WOPT(N1: ITER, NP)
                   CALL W1_GEMV( one, ELEMENTS( WOPT, N1, ITER, NP),  ELEMENT( WOPT, N1, NP), zero, CWORK1, 1, CQIJ)

!$ACC PARALLEL LOOP PRESENT(CTMP,CWORK1) __IF_ASYNC__
                   DO N2=N1,ITER
                      CTMP(N2,N1)=      GREAL(CWORK1(N2-N1+1))
                      CTMP(N1,N2)=GREAL(GCONJG(CWORK1(N2-N1+1)))
                   ENDDO
                ENDDO buildh
!$ACC EXIT DATA COPYOUT(CTMP,CHAM) DELETE(CWORK1) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)

#ifndef gammareal
                DO N1=1,ITER
                   IF (ABS(AIMAG(CHAM(N1,N1)))>1E-2_q) THEN
DOACC !$OMP CRITICAL (omp_wrt_stdout)
                      WRITE(*,*)'WARNING: Sub-Space-Matrix is not hermitian in rmm',AIMAG(CHAM(N1,N1))
DOACC !$OMP END CRITICAL (omp_wrt_stdout)
                   ENDIF
                   CHAM(N1,N1)= REAL( CHAM(N1,N1) ,KIND=q)
                ENDDO
#endif

                ! solve eigenvalue-problem and calculate lowest eigenvector
                ! this eigenvector corresponds to a minimal residuum
                ! CHAM(n1,n2) U(n2,1) = E(1) S(n1,n2)  U(n2,1)
#ifndef _OPENMP
                IF (.FALSE.) THEN
                   io_begin
                   NPL2=MIN(10,ITER)
                   WRITE(6,*)
                   DO N1=1,NPL2
                      WRITE(6,1)N1,(REAL( CHAM(N1,N2) ,KIND=q) ,N2=1,NPL2)
                   ENDDO
                   WRITE(6,*)
#ifndef gammareal
                   DO N1=1,NPL2
                      WRITE(6,3)N1,(AIMAG(CHAM(N1,N2)),N2=1,NPL2)
                   ENDDO
                   WRITE(6,*)
#endif
                   DO N1=1,NPL2
                      WRITE(6,1)N1,(REAL( CTMP(N1,N2) ,KIND=q) ,N2=1,NPL2)
                   ENDDO
                   WRITE(6,*)
#ifndef gammareal
                   DO N1=1,NPL2
                      WRITE(6,3)N1,(AIMAG(CTMP(N1,N2)),N2=1,NPL2)
                   ENDDO
                   WRITE(6,*)
#endif

1                  FORMAT(1I2,3X,20F9.5)
3                  FORMAT(1I2,3X,20E9.1)
                   io_end
                ENDIF
#endif
                !
                ! onion award of the year for IBM,
                ! who use a completely different DSYGV calling sequence
                !
#ifdef gammareal
#ifdef essl
                CALL DSYGV &
                     &  (1,CHAM,NRES,CTMP,NRES,R,CHAM,NRES, &
                     &           ITER,CWORK(1),LWORK*INFO%NDAV)
#else
                CALL DSYGV &
                     &  (1,'V','U',ITER,CHAM,NRES,CTMP,NRES,R, &
                     &           CWORK(1),LWORK*INFO%NDAV,IFAIL)
#endif
#else
                PROFILING_START('zhegv')

                CALL ZHEGV &
                     &  (1,'V','U',ITER,CHAM,NRES,CTMP,NRES,R, &
                     &           CWORK(1),LWORK*INFO%NDAV,RWORK(1),IFAIL)

                PROFILING_STOP('zhegv')
#endif
                ! just to be sure merge results from all nodes
                CALLMPI( M_bcast_g(WDES%COMM_INB, CHAM, NRES*NRES))

                IF (IFAIL/=0) THEN
DOACC !$OMP CRITICAL (omp_wrt_stdout)
                   IF (IU6>=0) &
                        WRITE(IU6,219) IFAIL,ITER,N
                   IF (IU0>=0) &
                        WRITE(IU0,219) IFAIL,ITER,N
DOACC !$OMP END CRITICAL (omp_wrt_stdout)
                   !  try to save things somehow, goto next band
                   W1(NP)%LDO=.FALSE.
                   CYCLE i2
                ENDIF

219             FORMAT('WARNING in EDDRMM: call to ZHEGV failed, returncode =',I4,1X,I2,1X,I6)
                FPRE(NP)=R(1)
#ifndef _OPENMP
                IF (IDUMP==2)  WRITE(*,'(E9.2,"P")',ADVANCE='NO') SQRT(ABS(FPRE(NP)))

                !     write out 'optimal trial step' i.e step which would have minimized
                !     the residuum
                IF (ITER==2 .AND. IDUMP==2) THEN
                   OTRIAL= REAL( 1+CHAM(1,1)/CHAM(2,1) ,KIND=q) *TRIAL(NP)
                   WRITE(*,'(1X,F7.4,"o")',ADVANCE='NO') OTRIAL
                ENDIF
#endif
                !     some heuristic for numerical accuracy problems
                !     small residuum and negative step -> stop immediately
                IF (ITER==2) THEN
                   OTRIAL= REAL( 1+CHAM(1,1)/CHAM(2,1) ,KIND=q) *TRIAL(NP)
                   IF (OTRIAL <0 .AND.  ABS(FPRE(NP))< 1E-9_q) THEN
#ifndef NO_NUMPROB_WARN
DOACC !$OMP CRITICAL (omp_wrt_stdout)
                      IF (IU0>=0) WRITE(IU0,'(" num prob ")',ADVANCE='NO')
DOACC !$OMP END CRITICAL (omp_wrt_stdout)
#endif
                      W1(NP)%LDO=.FALSE.
                      CYCLE i2
                   ENDIF
                ENDIF
#ifndef _OPENMP
                IF (.FALSE.) THEN
                   io_begin
                   NPL2=MIN(10,ITER)

                   WRITE(77,*)
                   DO N1=1,NPL2
                      WRITE(77,1)N1,R(N1),(REAL( CHAM(N2,N1) ,KIND=q) ,N2=1,NPL2)
                   ENDDO
                   WRITE(77,*)
#ifndef gammareal
                   DO N1=1,NPL2
                      WRITE(77,3)N1,R(N1),(AIMAG(CHAM(N2,N1)),N2=1,NPL2)
                   ENDDO
                   WRITE(77,*)
#endif
                   io_end
                ENDIF
#endif
             ELSE optsub
!***********************************************************************
! inverse interation step (usually not used)
! minimize
!    | ( H - e S) | phi > - | phi_ini > |^ 2  -> min
! in the yet available subspace spanned by the wavefunction stored in WOPT%CW
! if one denotes these wavefunctions as phi_j, and R_j=  (H - e S) phi_j
! the following equation is obtained:
!  sum_ij  b_i* < R_i | R_j > b_j - sum_i b_i*<R_i |phi_ini> + c.c. -> min
! or equivalently
!  sum_j  < R_i | R_j > b_j  = <R_i |phi_ini>
! the new optimized wavefunction is given by solving this linear
! equation for b
!
!***********************************************************************
!$ACC ENTER DATA CREATE(CWORK1,CHAM,B) __IF_ASYNC__
!$ACC KERNELS PRESENT(CHAM,B) __IF_ASYNC__
                CHAM=0
                B   =0
!$ACC END KERNELS
                !    A(n2,n1)=    < phi_n2 |  ( H - e S) ( H - e S)  | phi_n1 >
                builda: DO N1=1,ITER
                   CALL W1_GEMV( one, ELEMENTS( WOPT, NRES+N1, NRES+ITER, NP), ELEMENT( WOPT, NRES+N1, NP), zero, CWORK1, 1)
!$ACC PARALLEL LOOP PRESENT(CHAM,CWORK1) __IF_ASYNC__
                   DO N2=N1,ITER
                      CHAM(N2,N1)=       GREAL(CWORK1(N2-N1+1))
                      CHAM(N1,N2)=GREAL(GCONJG(CWORK1(N2-N1+1)))
                   ENDDO
                ENDDO builda

                ! B(n1) = < phi_n1 |  ( H - e S) | phi_ini >
                CALL W1_GEMV( one, ELEMENTS( WOPT, NRES+1, NRES+ITER, NP), ELEMENT( W_INI,  NP), zero, B(1), 1)

!$ACC EXIT DATA COPYOUT(CHAM,B) DELETE(CWORK1) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)

                CTMP=CHAM
                B_=B
                ! calculate the solution of sum_j CHAM(i,j) * X(j) = B(i)
                ! overwrite B by X
                PROFILING_START('ggetrf')

                CALL GGETRF( ITER, ITER, CHAM, NRES, IPIV, IFAIL )

                PROFILING_STOP('ggetrf')

                IF (IFAIL ==0) CALL GGETRS('N', ITER, 1, CHAM, NRES, IPIV, B, NRES, IFAIL)

#ifndef _OPENMP
                IF (.FALSE.) THEN
                   ! dump the matrix and the solution vector
                   io_begin
                   N2=MIN(10,ITER)
                   WRITE(6,*)
                   DO N1=1,N2
                      WRITE(6,'("m",I3,8E14.7)')N1, CTMP(N1,1:N2)
                   ENDDO
                   WRITE(6,'(A4,8E14.7)') 'b', B_(1:N2)

                   WRITE(6,*)
                   WRITE(6,'(A4,8E14.7)') 'e', B (1:N2)
                   io_end
                ENDIF
#endif
                IF (ITER == 1) B(1) = 1

                IF (IFAIL/=0) THEN
DOACC !$OMP CRITICAL (omp_wrt_stdout)
                   IF (IU6>=0) &
                        WRITE(IU6,219) IFAIL,ITER,N
                   IF (IU0>=0) &
                        WRITE(IU0,219) IFAIL,ITER,N
DOACC !$OMP END CRITICAL (omp_wrt_stdout)
                   !  try to save things somehow, goto next band
                   W1(NP)%LDO=.FALSE.
                   CYCLE i2
                ENDIF
#ifndef _OPENMP
                IF (ITER==2 .AND. IDUMP==2) THEN
                   ! write out 'optimal trial step' i.e step which would have minimized
                   ! the residuum
                   OTRIAL= REAL( 1+B(1)/B(2) ,KIND=q) *TRIAL(NP)
                   WRITE(*,'(1X,F7.4,"o")',ADVANCE='NO') OTRIAL
                ENDIF
#endif

                CHAM(1:ITER,1)=B(1:ITER)

#ifndef _OPENMP
                IF (IDUMP >= 3) THEN
                   ! set C to < xi | xi >
                   C=W1_DOT( ELEMENT(W_INI,NP) , ELEMENT(W_INI,NP))
                   DO N1=1,ITER
                      DO N2=1,ITER
                         C=C+GCONJG(B(N2))*CTMP(N2,N1)*B(N1)
                      ENDDO
                      C=C-B_(N1)*GCONJG(B(N1))-GCONJG(B_(N1))*B(N1)
                   ENDDO
                   ! residual after the step
                   WRITE(*,'(1X,E9.2,"rs")',ADVANCE='NO') SQRT(ABS(C))
                ENDIF
#endif
             ENDIF optsub
!=======================================================================
! now perform the trial step (ITER > 1 use previous trial step)
! but restrict trial step to 1.0
!=======================================================================
             ! break when relative error has decreased sufficiently
             IF (ABS(FPRE(NP))<DEIT(NP)) LABORT(NP)=.TRUE.
             ! break when finished
             IF (ITER == NITER) LABORT(NP)=.TRUE.
             ! quick break on empty orbitals
             IF (ABS(W%FERWE(N,NK,ISP))<INFO%WEIMIN .AND. ITER >= 2) LABORT(NP)=.TRUE.
             
             ! when FBREAK_MEAN_IN is not present, set LABORT to .FALSE. to recover old RMM-DIIS behaviour
             ! this means we finish with a trial step and NOT the vector that minimizes the residuum
             IF (.NOT. PRESENT(FBREAK_MEAN_IN)) LABORT(NP)=.FALSE.

             ! make trial step positive
             IF (TRIAL(NP)<0) TRIAL(NP)=ABS(TRIAL(NP))
             ! WTMP(NP)=0
             CALL W1_DSCAL( WTMP(NP), 0.0_q)

             ! WTMP=WTMP + CHAM(I,1) *WOPT(I,NP)
             DO I=1,ITER
                CALL W1_GAXPY( ELEMENT(WOPT, I,NP), CHAM(I,1), WTMP(NP))
             ENDDO

             IF (.NOT.LABORT(NP)) THEN
                ! perform the trial step on the orbital
                DO I=1,ITER
                   CALL APPLY_PRECOND( ELEMENT(WOPT, NRES+I,NP), W1(NP), PRECON(1,NP))
                   CALL W1_GAXPY( W1(NP), -CHAM(I,1)*TRIAL(NP), WTMP(NP))
                ENDDO
             ENDIF

             ! transform the wave-function to real space
             CALL FFTWAV_W1(WTMP(NP))
!***********************************************************************
! for first iteration (ITER=1) or LDELAY
! minimize energy starting from current wavefunctions stored in
! W1 along the current searchdirection stored in WOPT(NRES+ITER,NP)
!
!***********************************************************************
          ELSE optsubspace
             IF (ITER==1) THEN
                DEIT(NP)=ABS(FPRE(NP))*INFO%DEPER
                IF (PRESENT(FBREAK_MEAN_IN)) THEN
                   ! use previous average break condition as upper bound
                   DEIT(NP)=MIN(DEIT(NP), FBREAK_MEAN_IN)
                   ! determine new average break condition
                   FBREAK_MEAN=FBREAK_MEAN+MIN(DEIT(NP), FBREAK_MEAN_IN)
                   ICOU_FBREAK=ICOU_FBREAK+1
                ENDIF
             ENDIF
             IF (ITER == NITER) LABORT(NP)=.TRUE.

             ! trial vector in line minimization
             CALL APPLY_PRECOND( ELEMENT(WOPT, NRES+ITER,NP), WTMP(NP), PRECON(1,NP))

             CALL FFTWAV_W1(WTMP(NP))
          ENDIF optsubspace
       ENDDO i2
DOACC !$OMP END PARALLEL DO

#ifdef _OPENACC
       CALL ACC_SYNC_ASYNC_Q(NSIM,W1(:)%LDO)
       ACC_ASYNC_Q=ACC_ASYNC_ASYNC
#endif

       PROFILING_STOP('i2')

       ! calculate results of projection operatores
       WTMP%LDO=W1%LDO
       CALL W1_PROJALL(WDES1, WTMP, NONLR_S, NONL_S, NSIM)

!$ACC WAIT(ACC_ASYNC_ASYNC) IF(ACC_EXEC_ON)

       PROFILING_START('i3')

DOACC !$OMP PARALLEL DO PRIVATE(NP,N,ITER,A2,WSCAL,A1,DE,DECEL) SCHEDULE(dynamic) &
DOACC !$OMP NUM_THREADS(omp_nthreads_acc) REDUCTION(+:DESUM,ICOUEV)
       i3: DO NP=1,NSIM
          N=NB(NP); ITER=IT(NP); IF (.NOT. W1(NP)%LDO) CYCLE i3

          SET_ACC_ASYNC_Q(NP)

!***********************************************************************
!
! finish trial step of RMM-DIIS/ inverse iteration by copying WTMP to W1
!
!***********************************************************************
          mine: IF (.NOT. LDELAY .AND. ITER > 1) THEN

             CALL W1_COPY(WTMP(NP), W1(NP))
!***********************************************************************
! for LDELAY or ITER=1
! finish line minimization of <phi| H | phi> along trial direction
!
!***********************************************************************
          ELSE mine
             ! < g | phi > 1. order energy change
             A2=W1_DOT( ELEMENT( WOPT, NRES+ ITER, NP), WTMP(NP))
             A2=-A2*2
             CALL W1_DAXPY(WTMP(NP), -1.0_q, W1(NP))

             IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
                CALL ECCP_TAU(WDES1,W1(NP),W1(NP),LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP),LATT_CUR,HAMILTONIAN%MU(:,ISP),W%CELEN(N,NK,ISP))
             ELSE
                CALL ECCP(WDES1,W1(NP),W1(NP),LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP), W%CELEN(N,NK,ISP))
             ENDIF
             CALL CNORMA(W1(NP),CQIJ, ISP, WSCAL)

             W%CELEN(N,NK,ISP) =W%CELEN(N,NK,ISP)*WSCAL**2
             A1= W%CELEN(N,NK,ISP)-EVALUE(NP)

             ! quadratic interpolation to find the minimum
             TRIAL(NP)= -A2/(A1-A2)/2
             DE       = (A2+(A1-A2)*TRIAL(NP))*TRIAL(NP)

             ! avoid too large trial steps
#ifndef _OPENMP
             IF (IDUMP>=2) WRITE(*,'(1X,F7.4,"T")',ADVANCE='NO') TRIAL(NP)
#endif
             IF (.NOT. LDELAY) THEN
                IF (TRIAL(NP)>0 .AND.(TRIAL(NP) > 1)) TRIAL(NP)= 1
                IF (TRIAL(NP)>0 .AND.(TRIAL(NP)<0.1)) TRIAL(NP)=0.1
!TODO is there any use to allow for negative trial steps 
!I think right now it kind of sucks (in all but the first step ABS(TRIAL(NP)) is used
! so probably does not matter what is used here
                IF (TRIAL(NP)<0 .AND.(TRIAL(NP) <-1)) TRIAL(NP)=-1
                IF (TRIAL(NP)<0 .AND.(TRIAL(NP)>-.1)) TRIAL(NP)=-0.1
             ENDIF
#ifndef _OPENMP
             IF (IDUMP>=2) WRITE(*,'(1X,F7.4,"T")',ADVANCE='NO') TRIAL(NP)
#endif
             ! set W1 finally
             CALL W1_DAXPY(WTMP(NP), -(TRIAL(NP)-1), W1(NP))
          ENDIF mine
!=======================================================================
! common code
!=======================================================================
          ! inverse of norm
          CALL CNORMA(W1(NP),CQIJ, ISP, WSCAL)
          ! scale W1
          CALL W1_DSCAL(W1(NP), WSCAL)

          IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
             CALL ECCP_TAU(WDES1,W1(NP),W1(NP),LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP),LATT_CUR,HAMILTONIAN%MU(:,ISP),W%CELEN(N,NK,ISP))
          ELSE
             CALL ECCP(WDES1,W1(NP),W1(NP),LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP), W%CELEN(N,NK,ISP))
          ENDIF
          DECEL =W%CELEN(N,NK,ISP)-EVALUE(NP)
          DE    =DECEL
#ifndef _OPENMP
          IF (IDUMP==2) WRITE(*,'(E10.2,2H |)',ADVANCE='NO') DECEL
#endif
          DESUM =DESUM +WDES%RSPIN*WDES%WTKPT(NK)*W%FERWE(N,NK,ISP)*DECEL
          ICOUEV=ICOUEV+1
!=======================================================================
! break of intra-band-minimisation
! at the moment we performe a break of the intra-band minimization if
! ) DE is less then INFO%DEPER % of the change in the first minimization
!     of this band (relative breakcondition)
! ) DE less then INFO%EBREAK (absolut breakcondition)
! ) if unoccupied band break after 2. iteration
!=======================================================================
          DE=ABS(DE)

          IF (LABORT(NP)) W1(NP)%LDO=.FALSE.
          ! stop if change of the one-electron energy is smaller than EBREAK
          ! saves quite a lot and not unstable (for the tested systems)
          IF (DE<INFO%EBREAK) W1(NP)%LDO=.FALSE.
          
          IF (.NOT. PRESENT(FBREAK_MEAN_IN)) THEN
          ! old code behaviour, test for all break conditions again
             IF (ABS(FPRE(NP))<DEIT(NP)) W1(NP)%LDO=.FALSE.
             IF (ITER == NITER) W1(NP)%LDO=.FALSE.
             ! stop if change of the one-electron energy is smaller than EBREAK
             IF (DE<INFO%EBREAK) W1(NP)%LDO=.FALSE.
             IF (ABS(W%FERWE(N,NK,ISP))<INFO%WEIMIN .AND. &
                  &    ABS(W%FERWE(N,NK,ISP)*DE)<INFO%WEIMIN .AND. ITER >= 2) W1(NP)%LDO=.FALSE.
          ENDIF
       ENDDO i3
DOACC !$OMP END PARALLEL DO

#ifdef _OPENACC
       CALL ACC_SYNC_ASYNC_Q(NSIM,W1(:)%LDO)
       ACC_ASYNC_Q=ACC_ASYNC_ASYNC
#endif

       PROFILING_STOP('i3')

       ! one band just finished ?, set NB(NP) also to 0 and finish everything
       DO NP=1,NSIM
          N=NB(NP)
          IF (.NOT. W1(NP)%LDO .AND. N /=0 ) THEN
             NB(NP)=0
             IF (IDUMP==2)  WRITE(*,'(F9.4,2H q)')REAL( W%CELEN(N,NK,ISP) ,KIND=q)
             IF (IDUMP==10) WRITE(*,*)
          ENDIF
       ENDDO

    ENDDO bands

    ENDDO kpoints
    ENDDO spin

    !=======================================================================
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

    IF (W%OVER_BAND) THEN
       W%OVER_BAND=.FALSE.
       CALL REDIS_PW_DEALLOC(H_PW)
    ENDIF

    DO NP=1,NSIM
       CALL DELWAV_R(W1(NP))
       CALL DELWAV(WTMP(NP) ,.TRUE.)
    ENDDO
!$ACC EXIT DATA DELETE(WTMP(:),W1(:)) __IF_ASYNC__
    CALL DELWAVA(W_INI)
    CALL DELWAVA(WOPT)
    CALL DELWAVA(WTMPA)
!$ACC EXIT DATA DELETE(W_INI,WOPT,WTMPA) __IF_ASYNC__
!$ACC EXIT DATA DELETE(SV,CDIJ,CQIJ,PRECON) __IF_ASYNC__
    DEALLOCATE(PRECON,CHAM,CTMP,CWORK1,B,IPIV,B_)
#ifdef _OPENACC
    CALL ACC_DELETE_TYPED_VAR(NONL_S)
    CALL ACC_DELETE_TYPED_VAR(NONLR_S)
    CALL ACC_DELETE_TYPED_VAR(GRID)
    CALL ACC_DELETE_TYPED_VAR(WDES1)
!$ACC WAIT IF(ACC_EXEC_ON)
    ACC_ASYNC_Q=ACC_ASYNC_SYNC
#endif
    PROFILING_STOP('eddrmm')

    RETURN
  END SUBROUTINE EDDRMM

END MODULE rmm_diis
