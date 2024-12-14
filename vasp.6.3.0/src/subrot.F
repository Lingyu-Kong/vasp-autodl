!#define timing
#include "symbol.inc"
#ifndef _OPENACC
#define USE_ZHEEVX
#endif

MODULE subrot
  USE prec
  USE dfast
  USE hamil
  USE pead, ONLY : LUSEPEAD,LPEAD_NO_SCF,W_STORE,PEAD_EIGENVALUES,PEAD_ACC_ADD_CPROJ,PEAD_ACC_ADD_PW

  ! electron volts by which states < WDES%NBANDSLOW-1 should be shifted
  REAL(q), SAVE :: RSHIFT_NBLOW  = 250._q
CONTAINS
!************************ SUBROUTINE EDDIAG ****************************
! RCS:  $Id: subrot.F,v 1.10 2003/06/27 13:22:23 kresse Exp kresse $
!
! this subroutine calculates the electronic eigenvalues and
! optionally performes a  sub-space diagonalisation
! i.e. unitary transforms the wavefunctions so that the Hamiltonian
!  becomes diagonal in the subspace spanned by the wavefunctions
! IFLAG:
!  0 only eigenvalues (without  diagonalisation no sub-space matrix)
!  1 only eigenvalues and sub-space matrix (no diagonalisation)
!  2 eigenvalues using diagonalisation of sub-space matrix
!    do not rotate wavefunctions
! 12 eigenvalues using diagonalisation of sub-space matrix
!    no kinetic energy
!  3 eigenvalues and sub-space diagonalisation rotate wavefunctions
! 13 (for 13 no Jacobi algorithm is allowed)
! 23 eigenvalues and sub-space diagonaliation in the occupied 
!    manifold (occupancies > 1-1E-3)
!  4 eigenvalues and sub-space diagonalisation rotate wavefunctions
!    using Loewdin pertubation theory (conserves ordering)
!  5 eigenvalues and sub-space diagonalisation + orthogonalization
!    unfortunately this option turns out to slow down the 
!    convergence of IALGO=48
!
! notes on the matrix CHAM used internally
! it is defined as 
!  CHAM(n2,n1)   = < phi_n2 | H | phi_n1 >
! however the returned matrix CHAMHF is defined as the complex conjugated 
!  CHAMHF(n2,n1) = < phi_n2 | H | phi_n1 >* = < phi_n1 | H | phi_n2 >
! since this is what rot.F uses
!
!  written by gK
!  last update Sep 24 2004 massive cleanup
!  optional arguments:
!  NBANDS_MAX   maximum number of bands for which eigenvalues are
!               calculated (only supported for IFLAG=0)
!  TODO: Martijn new arguments
!
!> @details @ref openacc :
!
!***********************************************************************

  SUBROUTINE EDDIAG(HAMILTONIAN, &
       GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
       LMDIM,CDIJ,CQIJ,IFLAG,SV,T_INFO,P,IU0,EXHF, & 
       CHAMHF,LFIRST,LLAST,NKSTART,NKSTOP,NBANDS_MAX,EXHF_ACFDT)
#ifdef _OPENACC
    USE mopenacc
#endif
    USE prec
    USE wave_high
    USE lattice
    USE mpimy
    USE mgrid
    USE nonl_high
    USE hamil_struct_def
    USE constant
    USE jacobi
    USE scala
    USE main_mpi
    USE fock
    USE pseudo
    USE ini
    USE sym_grad
    USE fileio
    USE openmp, ONLY : omp_nthreads
#ifdef fock_dblbuf
    USE fock_dbl
#endif

    IMPLICIT NONE
    TYPE (ham_handle)  HAMILTONIAN
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W
    TYPE (wavedes)     WDES
    TYPE (symmetry) :: SYMM
    INTEGER LMDIM
    OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ),CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    INTEGER            IFLAG            ! determines mode of diagonalisation
    RGRID   SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ) ! local potential
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    INTEGER IU0
    REAL(q) EXHF
! if CHAMHF is present, CHAMHF is added to the subspace matrix at each iteration
! in addition calculation of Fock part is bypassed
    GDEF, OPTIONAL :: CHAMHF(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)
! if LFIRST is supplied CHAMHF is set to CONJG(CHAMHF)-CHAM
! if LLAST  is supplied CHAMHF is set to CONJG(CHAM)
    LOGICAL, OPTIONAL :: LFIRST, LLAST
    INTEGER, OPTIONAL :: NKSTART, NKSTOP     ! start k-point
    INTEGER, OPTIONAL :: NBANDS_MAX  ! maximum band index
    LOGICAL :: LSCAAWARE_LOCAL
    REAL(q), OPTIONAL :: EXHF_ACFDT
    
! local
    ! work arrays for ZHEEV (blocksize times number of bands)
    INTEGER, PARAMETER :: LWORK=32
    GDEF       CWRK(LWORK*WDES%NB_TOT)
    REAL(q)    R(WDES%NB_TOT)
#ifndef USE_ZHEEVX
    REAL(q)    RWORK(3*WDES%NB_TOT)
#else
    REAL(q)    RWORK(7*WDES%NB_TOT), ABSTOL, VL, VU
    INTEGER IWORK(5*WDES%NB_TOT), INFO(WDES%NB_TOT),IL, IU, NB_CALC
#endif
    ! work arrays (do max of 16 strips simultaneously)

    TYPE (wavedes1)    WDES1          ! descriptor for one k-point
    TYPE (wavefun1)    W1             ! current wavefunction
    TYPE (wavefuna)    WA             ! array to store wavefunction
    TYPE (wavespin)    WFOCK          ! array to store the Fock (exchange contribution)
    TYPE (wavefuna)    WNONL          ! array to hold non local part D * wave function character
    TYPE (wavefuna)    WOVL           ! array to hold non local part Q * wave function character
    TYPE (wavefuna)    WHAM           ! array to store accelerations for a selected block
    COMPLEX(q)         CDCHF          ! HF double counting energy
    INTEGER :: ICALL=0                ! number of calls
    GDEF,ALLOCATABLE,TARGET::  CHAM(:,:),COVL(:,:) ! Hamiltonian and overlap matrix
    TYPE (wavefun1)    WTMP(NSTRIP_STANDARD)

    ! array for asyncronous data redistribution
    TYPE (REDIS_PW_CTR),POINTER :: H_PW1, H_PW2

    LOGICAL, EXTERNAL :: USEFOCK_CONTRIBUTION

    INTEGER :: NB_TOT, NBANDS, NSTRIP, ISP, NK, N, NP, NPOS, NSTRIP_ACT, &
         NPOS_RED, NSTRIP_RED, IFAIL, MY_NKSTART, MY_NKSTOP, NSIM_LOCAL
!$  INTEGER NSIM_FOCK

    PROFILING_START('eddiag')

#ifdef _OPENACC_TST
!$ACC WAIT IF(ACC_EXEC_ON)
    PUSH_ACC_EXEC_ON(.TRUE.)

    ACC_ASYNC_Q = ACC_ASYNC_ASYNC

    CALL ACC_COPYIN_TYPED_VAR(WDES)
    CALL ACC_COPYIN_TYPED_VAR(W)
    CALL ACC_COPYIN_TYPED_VAR(P)
    CALL ACC_COPYIN_TYPED_VAR(GRID)
    CALL ACC_COPYIN_TYPED_VAR(NONL_S)
    CALL ACC_COPYIN_TYPED_VAR(NONLR_S)

!$ACC ENTER DATA COPYIN(SV,CDIJ,CQIJ,LATT_CUR) __IF_ASYNC__
!$ACC ENTER DATA COPYIN(CHAMHF) ASYNC(ACC_ASYNC_Q) IF(ACC_EXEC_ON.AND.PRESENT(CHAMHF))
!$ACC ENTER DATA COPYIN(HAMILTONIAN,HAMILTONIAN%MU) IF(ASSOCIATED(HAMILTONIAN%MU).AND.ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(HAMILTONIAN,HAMILTONIAN%AVEC) IF(ASSOCIATED(HAMILTONIAN%AVEC).AND.ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
#endif

    LscaAWARE_LOCAL=__IF_ACC_OFF__(LSCAAWARE.AND.(.NOT.(IFLAG==4)))
    IF (PRESENT(CHAMHF) .OR. IFLAG==1 .OR. IFLAG==23) LSCAAWARE_LOCAL=.FALSE.

    CDCHF=0         ! double counting HF
    IF (PRESENT(EXHF_ACFDT)) EXHF_ACFDT=0
    ICALL=ICALL+1

    NB_TOT=WDES%NB_TOT
    NBANDS=WDES%NBANDS
    IF (PRESENT(NBANDS_MAX) .AND. IFLAG == 0) NBANDS=NBANDS_MAX

    NSTRIP=NSTRIP_STANDARD
!$  NSTRIP=NSTRIP*omp_nthreads
 
    ! allocate work space
    IF (.NOT. LscaAWARE_LOCAL) THEN
       ALLOCATE(CHAM(NB_TOT,NB_TOT))
!$ACC ENTER DATA CREATE(CHAM) IF(ACC_EXEC_ON)
    ELSE
       CALL INIT_scala(WDES%COMM_KIN, NB_TOT)
       ALLOCATE(CHAM(SCALA_NP(),SCALA_NQ()))
    ENDIF
!$ACC ENTER DATA CREATE(WDES1) IF(ACC_EXEC_ON)
    CALL SETWDES(WDES,WDES1,0)
!$ACC ENTER DATA CREATE(WHAM,WNONL) IF(ACC_EXEC_ON)
    CALL NEWWAVA(WHAM, WDES1, NSTRIP)
    CALL NEWWAVA_PROJ(WNONL, WDES1)
!$ACC ENTER DATA CREATE(W1) IF(ACC_EXEC_ON)
    CALL NEWWAV_R(W1, WDES1)

    IF (IFLAG==5) THEN
       ALLOCATE(COVL(NB_TOT,NB_TOT))
!$ACC ENTER DATA CREATE(COVL,WOVL) IF(ACC_EXEC_ON)
       CALL NEWWAVA_PROJ(WOVL, WDES1)
    ENDIF

    IF (PRESENT(NKSTART)) THEN
       MY_NKSTART=NKSTART
    ELSE
       MY_NKSTART=1
    ENDIF
    IF (PRESENT(NKSTOP)) THEN
       MY_NKSTOP=NKSTOP
    ELSE
       MY_NKSTOP=WDES%NKPTS
    ENDIF

!=======================================================================
! start with HF part and store results WFOCK
!=======================================================================
    IF (USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF))) THEN

       CALL ALLOCW(WDES,WFOCK)
#ifdef _OPENACC
       CALL ACC_COPYIN_TYPED_VAR(WFOCK)
#endif
#ifndef fock_dblbuf
       NSIM_LOCAL=(W%WDES%NSIM*2+W%WDES%NB_PAR-1)/W%WDES%NB_PAR
       IF (LPEAD_NO_SCF()) THEN
          DO N=1,NSIM_LOCAL
             CALL NEWWAV(WTMP(N) , WDES1, .FALSE.)
          ENDDO
       ENDIF

       DO ISP=1,WDES%ISPIN
       DO NK=MY_NKSTART,MY_NKSTOP
#ifdef MPI
          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
          CALL SETWDES(WDES,WDES1,NK)

          IF (NONLR_S%LREAL) THEN
             CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES)
          ELSE
             CALL PHASE(WDES,NONL_S,NK)
          ENDIF

#ifndef _OPENMP
          DO NPOS=1,NBANDS,NSIM_LOCAL
             NSTRIP_ACT=MIN(NBANDS+1-NPOS,NSIM_LOCAL)
             IF (LPEAD_NO_SCF()) THEN
             ! this is a little bit tricky/dirty
             ! if local field effects are not taken into account in the pead,
             ! (i.e. Hamiltonian constructed from initial wavefunctions)
             ! W_STORE (original wavefunctions) must be passed to FOCK_ACC and WTMP 
             ! (wavefunction on which the action is calculated)
             ! must be supplied as an additional argument
                DO N=NPOS,NPOS+NSTRIP_ACT-1
                   NP=N-NPOS+1
                   CALL W1_COPY(ELEMENT(W,WDES1,N,ISP),WTMP(NP))
                ENDDO
                CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W_STORE,   &
                     NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIP_ACT, &
                     WFOCK%CW(:,NPOS:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF, WTMP, LSYMGRAD=LSYMGRAD)
             ELSE
                IF (PRESENT(EXHF_ACFDT)) THEN
                   CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W,  &
                        NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIP_ACT, &
                        WFOCK%CW(:,NPOS:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF,  LSYMGRAD=LSYMGRAD, EXHF_ACFDT=EXHF_ACFDT)
                ELSE
                   CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W,  &
                        NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIP_ACT, &
                        WFOCK%CW(:,NPOS:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF,  LSYMGRAD=LSYMGRAD)
                ENDIF
             ENDIF
          ENDDO
#else
          IF (LPEAD_NO_SCF()) THEN
             DO NPOS=1,NBANDS,NSIM_LOCAL
                NSTRIP_ACT=MIN(NBANDS+1-NPOS,NSIM_LOCAL)
                ! this is a little bit tricky/dirty
                ! if local field effects are not taken into account in the pead,
                ! (i.e. Hamiltonian constructed from initial wavefunctions)
                ! W_STORE (original wavefunctions) must be passed to FOCK_ACC and WTMP 
                ! (wavefunction on which the action is calculated)
                ! must be supplied as an additional argument
                DO N=NPOS,NPOS+NSTRIP_ACT-1
                   NP=N-NPOS+1
                   CALL W1_COPY(ELEMENT(W,WDES1,N,ISP),WTMP(NP))
                ENDDO
                CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W_STORE,   &
                     NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIP_ACT, &
                     WFOCK%CW(:,NPOS:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF, WTMP, &
                     LSYMGRAD=LSYMGRAD)
             ENDDO
          ELSE
             NSIM_FOCK=NBLOCK_FOCK
             DO NPOS=1,MIN(W%WDES%NB_TOTK(NK,ISP),NBANDS*W%WDES%NB_PAR),NSIM_FOCK
                NSTRIP_ACT=MIN(NSIM_FOCK,W%WDES%NB_TOTK(NK,ISP)-NPOS+1)
                IF (PRESENT(EXHF_ACFDT)) THEN
                   CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W, &
                        NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIP_ACT, &
                        WFOCK%CW(:,:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF, &
                        LSYMGRAD=LSYMGRAD, EXHF_ACFDT=EXHF_ACFDT, &
                        LCOMMENSURATE=.FALSE.)
                ELSE
                   CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W, &
                        NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIP_ACT, &
                        WFOCK%CW(:,:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF, &
                        LSYMGRAD=LSYMGRAD, LCOMMENSURATE=.FALSE.)
                ENDIF
             ENDDO
          ENDIF
#endif
       ENDDO
       ENDDO

       CALLMPI( M_sum_z(WDES%COMM_KIN,CDCHF,1))
       CALLMPI( M_sum_z(WDES%COMM_KINTER,CDCHF,1))
       EXHF=CDCHF
 
       IF (PRESENT(EXHF_ACFDT)) THEN
          CALLMPI( M_sum_d(WDES%COMM_KIN,EXHF_ACFDT,1))
          CALLMPI( M_sum_d(WDES%COMM_KINTER,EXHF_ACFDT,1))
       ENDIF
#else
       IF (LPEAD_NO_SCF()) THEN
          IF (PRESENT(EXHF_ACFDT)) THEN
             CALL FOCK_ALL_DBLBUF(WDES,W_STORE,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ, &
            &   EX=EXHF,EX_ACFDT=EXHF_ACFDT,NBMAX=NBANDS*W%WDES%NB_PAR,NKSTART=MY_NKSTART,NKSTOP=MY_NKSTOP, &
            &   LSYMGRAD=LSYMGRAD,XI=WFOCK,WP=W)
          ELSE
             CALL FOCK_ALL_DBLBUF(WDES,W_STORE,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ, &
            &   EX=EXHF,NBMAX=NBANDS*W%WDES%NB_PAR,NKSTART=MY_NKSTART,NKSTOP=MY_NKSTOP, &
            &   LSYMGRAD=LSYMGRAD,XI=WFOCK,WP=W)
          ENDIF
       ELSE
          IF (PRESENT(EXHF_ACFDT)) THEN
             CALL FOCK_ALL_DBLBUF(WDES,W,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ, &
            &   EX=EXHF,EX_ACFDT=EXHF_ACFDT,NBMAX=NBANDS*W%WDES%NB_PAR,NKSTART=MY_NKSTART,NKSTOP=MY_NKSTOP, &
            &   LSYMGRAD=LSYMGRAD,XI=WFOCK)
          ELSE
             CALL FOCK_ALL_DBLBUF(WDES,W,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ, &
            &   EX=EXHF,NBMAX=NBANDS*W%WDES%NB_PAR,NKSTART=MY_NKSTART,NKSTOP=MY_NKSTOP, &
            &   LSYMGRAD=LSYMGRAD,XI=WFOCK)
          ENDIF
       ENDIF
#endif
       IF (LSYMGRAD) &
            CALL APPLY_SMALL_SPACE_GROUP_OP( W, WFOCK, NONLR_S, NONL_S,P, T_INFO%NIONS, LATT_CUR, SYMM, CQIJ, .FALSE. , -1, MY_NKSTART)
    ENDIF
!=======================================================================
    spin:  DO ISP=1,WDES%ISPIN
    kpoint: DO NK=MY_NKSTART,MY_NKSTOP
#ifdef MPI
       IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
!=======================================================================
       IF (LscaAWARE_LOCAL) CALL INIT_scala(WDES%COMM_KIN, WDES%NB_TOTK(NK,ISP))

       CALL SETWDES(WDES,WDES1,NK)
!=======================================================================
!  IFLAG=0 calculate eigenvalues  only
!=======================================================================
       IF (IFLAG==0) THEN
          W%CELEN(:,NK,ISP)=0
          DO N=1,NBANDS
             ! transform wavefunction to real space
             ! and calculate eigenvalues calling ECCP, no redistribution !
             CALL SETWAV(W, W1, WDES1, N, ISP) ! allocation for W1%CR done above
             CALL FFTWAV_W1(W1)
!            IF (ASSOCIATED(HAMILTONIAN%AVEC)) THEN
!               CALL ECCP_VEC(WDES1,W1,W1,LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP),HAMILTONIAN%AVEC, W1%CELEN)
             IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
                CALL ECCP_TAU(WDES1,W1,W1,LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP),LATT_CUR,HAMILTONIAN%MU(:,ISP),W1%CELEN)    
             ELSE
                CALL ECCP(WDES1,W1,W1,LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP), W1%CELEN)
             ENDIF
             IF (USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF))) THEN
                W1%CELEN=W1%CELEN+W1_DOT( ELEMENT( W, WDES1, N, ISP), ELEMENT (WFOCK, WDES1, N, ISP))
             ENDIF

             W%CELEN(N,NK,ISP)=W1%CELEN
          ENDDO
          CALL PEAD_EIGENVALUES(W,NK,ISP)
          CYCLE kpoint
       ENDIF

       WA=ELEMENTS(W, WDES1, ISP)
!=======================================================================
!  IFLAG /= 0 calculate Hamiltonian CHAM
!=======================================================================
       !  caclulate D |cfin_n> (D = non local strength of PP)
       IF (WDES%DO_REDIS .AND. LASYNC) THEN
          CALL REDIS_PW_ALLOC(WDES, NSTRIP, H_PW1)
          CALL REDIS_PW_ALLOC(WDES, NSTRIP, H_PW2)
          DO NPOS=1,NSTRIP
             CALL REDIS_PW_START(WDES, WA%CW(1,NPOS), NPOS, H_PW1)
          ENDDO
       ENDIF

       CALL OVERL(WDES1,.TRUE.,LMDIM,CDIJ(1,1,1,ISP),WA%CPROJ(1,1),WNONL%CPROJ(1,1))
       DO N=1,NBANDS
          CALL PEAD_ACC_ADD_CPROJ(WNONL%CPROJ(:,N),N,NK,ISP)
       ENDDO

       IF (IFLAG==5) THEN
          CALL OVERL(WDES1,WDES1%LOVERL,LMDIM,CQIJ(1,1,1,ISP),WA%CPROJ(1,1),WOVL%CPROJ(1,1))
       ENDIF

       ! redistribute the wavefunction characters
       CALL REDISTRIBUTE_PROJ(WA)
       CALL REDISTRIBUTE_PROJ(WNONL)
       IF (IFLAG==5) CALL REDISTRIBUTE_PROJ(WOVL)
!$ACC KERNELS PRESENT(CHAM) __IF_ASYNC__
       CHAM=0
!$ACC END KERNELS
       strip: DO NPOS=1,NBANDS,NSTRIP
          NSTRIP_ACT=MIN(NBANDS+1-NPOS,NSTRIP)

          !  calculate V_{local} |phi> + T | phi >
          !  for a block containing NSTRIP wavefunctions

          ! set Fock contribution
          IF (USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF))) THEN
!$ACC KERNELS PRESENT(WHAM,WFOCK) __IF_ASYNC__
             WHAM%CW(:,1:NSTRIP_ACT)=WFOCK%CW(:,NPOS:NPOS+NSTRIP_ACT-1,NK,ISP)
!$ACC END KERNELS
          ENDIF

          DO N=NPOS,NPOS+NSTRIP_ACT-1
             NP=N-NPOS+1
             CALL SETWAV(W, W1, WDES1, N, ISP)
             CALL FFTWAV_W1(W1)
             IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
                CALL HAMILT_LOCAL_TAU(W1, SV, LATT_CUR, HAMILTONIAN%MU, ISP, WHAM%CW(:,NP), &
               &   USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF)), IFLAG/=12) 
             ELSE
                CALL HAMILT_LOCAL(W1, SV, ISP, WHAM%CW(:,NP), USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF)), IFLAG/=12)
             ENDIF
             CALL PEAD_ACC_ADD_PW(WHAM%CW(:,NP),N,NK,ISP)
             IF (WDES%DO_REDIS.AND. LASYNC) CALL REDIS_PW_START(WDES, WHAM%CW(1,NP), N, H_PW2)
          ENDDO
          ! redistribute wavefunctions
          ! after this redistributed up to and including 1...NPOS+NSTRIP_ACT
          IF (WDES%DO_REDIS) THEN
             IF (LASYNC) THEN
                DO N=NPOS,NPOS+NSTRIP_ACT-1
                   NP=N-NPOS+1
                   CALL REDIS_PW_STOP (WDES, WA%CW(1,N), N, H_PW1)
                   IF (N+NSTRIP<=NBANDS) &
                        CALL REDIS_PW_START(WDES, WA%CW(1,N+NSTRIP), N+NSTRIP, H_PW1)
                   CALL REDIS_PW_STOP (WDES, WHAM%CW(1,NP), N, H_PW2)
                ENDDO
             ELSE
                CALL REDISTRIBUTE_PW( ELEMENTS( WA, NPOS, NPOS-1+NSTRIP_ACT))
                CALL REDISTRIBUTE_PW( ELEMENTS( WHAM, 1, NSTRIP_ACT))
             ENDIF
          ENDIF

          NPOS_RED  =(NPOS-1)*WDES%NB_PAR+1
          NSTRIP_RED=NSTRIP_ACT*WDES%NB_PAR

          IF (.NOT. LscaAWARE_LOCAL) THEN
             CALL ORTH1('U', &
               WA%CW_RED(1,1),WHAM%CW(1,1),WA%CPROJ_RED(1,1), &
               WNONL%CPROJ_RED(1,NPOS_RED),NB_TOT, &
               NPOS_RED, NSTRIP_RED, WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CHAM(1,1))
          ELSE
             CALL ORTH1_DISTRI('U', &
               WA%CW_RED(1,1),WHAM%CW(1,1),WA%CPROJ_RED(1,1), &
               WNONL%CPROJ_RED(1,NPOS_RED),NB_TOT, &
               NPOS_RED, NSTRIP_RED, WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CHAM(1,1), & 
               WDES%COMM_KIN, WDES%NB_TOTK(NK,ISP))
          ENDIF
       ENDDO strip

       IF (WDES%DO_REDIS .AND. LASYNC) THEN
          CALL REDIS_PW_DEALLOC(H_PW1)
          CALL REDIS_PW_DEALLOC(H_PW2)
       ENDIF

       IF (.NOT. LscaAWARE_LOCAL) THEN
          CALLMPI( M_sum_g(WDES%COMM_KIN,CHAM(1,1),NB_TOT*NB_TOT))
          ! add lower triangle
!$ACC PARALLEL LOOP GANG PRESENT(CHAM) __IF_ASYNC__
          DO N=1,NB_TOT
!$ACC LOOP VECTOR
             DO NP=N+1,NB_TOT
                CHAM(NP,N)=GCONJG(CHAM(N,NP))
             ENDDO
          ENDDO
          IF (IFLAG==23) CALL RESTRICT_TO_OCCUPIED_ONLY(WDES%NB_TOTK(NK,ISP), CHAM, W%FERTOT(:,NK, ISP))
       ENDIF

       IF (PRESENT(CHAMHF).AND.PRESENT(LFIRST)) THEN
!$ACC KERNELS PRESENT(CHAMHF,CHAM) __IF_ASYNC__
          IF (LFIRST) THEN
             CHAMHF(:,:,NK,ISP)=GCONJG(CHAMHF(:,:,NK,ISP))-CHAM(:,:)
          ENDIF
          CHAM(:,:)=CHAM(:,:)+CHAMHF(:,:,NK,ISP)
!$ACC END KERNELS
       ENDIF
#ifdef debug
!$ACC UPDATE SELF(CHAM) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
       IF (IU0>=0) CALL DUMP_HAM( "Hamilton matrix subrot",WDES, CHAM)
#endif
!-----------------------------------------------------------------------
! calculate the overlap matrix
!-----------------------------------------------------------------------
       IF (IFLAG==5) THEN
!$ACC KERNELS PRESENT(COVL) __IF_ASYNC__
          COVL=(0._q,0._q)
!$ACC END KERNELS
          DO NPOS=1,NB_TOT-NSTRIP_STANDARD_GLOBAL,NSTRIP_STANDARD_GLOBAL
             CALL ORTH1('U',WA%CW_RED(1,1),WA%CW_RED(1,NPOS),WA%CPROJ_RED(1,1), &
                  WOVL%CPROJ_RED(1,NPOS),NB_TOT, &
                  NPOS,NSTRIP_STANDARD_GLOBAL,WDES1%NPL_RED,WDES1%NPRO_O_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,COVL(1,1))
          ENDDO

          CALL ORTH1('U',WA%CW_RED(1,1),WA%CW_RED(1,NPOS),WA%CPROJ_RED(1,1), &
               WOVL%CPROJ_RED(1,NPOS),NB_TOT, &
               NPOS,NB_TOT-NPOS+1,WDES1%NPL_RED,WDES1%NPRO_O_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,COVL(1,1))
          CALLMPI( M_sum_g(WDES%COMM_KIN,COVL(1,1),NB_TOT*NB_TOT))
#ifdef debug
!$ACC UPDATE SELF(COVL) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
          IF (IU0>=0) CALL DUMP_HAM( "Overlap matrix",WDES, COVL)
#endif
       ENDIF
!=======================================================================
! IFLAG =2
!  simply copy eigenvalues
!=======================================================================
       IF (IFLAG==12) THEN
!$ACC UPDATE SELF(CHAM) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
          IF (IU0>=0) CALL DUMP_HAM( "Hamilton matrix",WDES, CHAM)
          CALL vtutor%stopCode()
       ENDIF

       IF (.NOT. LscaAWARE_LOCAL) THEN
!$ACC PARALLEL LOOP PRESENT(WDES,W,CHAM) __IF_ASYNC__
          DO N=1,WDES%NB_TOTK(NK,ISP)
             W%CELTOT(N,NK,ISP)=CHAM(N,N)
          ENDDO
!$ACC UPDATE SELF(W%CELTOT(1:WDES%NB_TOTK(NK,ISP),NK,ISP)) __IF_ASYNC__
       ENDIF
!=======================================================================
! IFLAG =4 use Loewdin perturbation to get rotation matrix
! this preserves the ordering of the eigenvalues
! MIND: does not work for real matrices
!=======================================================================
       IF (IFLAG==4) THEN
          CALL LOEWDIN_DIAG(WDES%NB_TOTK(NK,ISP), NB_TOT, CHAM)

          CALL ORSP(WDES%NB_TOTK(NK,ISP), NB_TOT, NB_TOT, CHAM)
!test writegamma
!          IF (IU0>=0) CALL WRITEGAMMA(NK, WDES%NB_TOTK(NK,ISP), NB_TOT, CHAM, .TRUE.)
#ifdef debug
           IF (IU0>=0) CALL DUMP_HAM( "Loewdin",WDES, CHAM)
#endif
       ELSE
!=======================================================================
! IFLAG > 1 and IFLAG <4
! diagonalization of CHAM
! we have lots of choices for the parallel version
! this  makes things rather complicated
! to allow for reasonable simple programming, once the diagonalisation
! has been done I jump to line 100
!=======================================================================
          IF (IFLAG==1) GOTO 1000

#ifndef gammareal
          IF (.NOT. LscaAWARE_LOCAL) THEN
NOACC !$OMP PARALLEL DO DEFAULT(SHARED) &
NOACC !$OMP PRIVATE(N)
!$ACC KERNELS PRESENT(WDES,CHAM) __IF_ASYNC__
             DO N=1,WDES%NB_TOTK(NK,ISP)
                IF (ABS(AIMAG(CHAM(N,N)))>1E-2_q .AND. IU0>=0) THEN
NOACC !$OMP CRITICAL (omp_wrt_stdout)
NOACC              WRITE(IU0,'(A,I5,E14.3)')'WARNING in EDDIAG: sub space matrix is not hermitian',N,AIMAG(CHAM(N,N))
NOACC !$OMP END CRITICAL (omp_wrt_stdout)
                ENDIF
                CHAM(N,N)= REAL( CHAM(N,N) ,KIND=q)
             ENDDO
!$ACC END KERNELS
NOACC !$OMP END PARALLEL DO
          ELSE
             CALL BG_CHANGE_DIAGONALE(WDES%NB_TOTK(NK,ISP),CHAM(1,1),IU0)
          ENDIF
#endif

          !
          ! parallel versions
          ! if fast Jacobi method exists use it (T3D, T3E only)
          ! use the first line in that case
          CALL SHIFT_BANDS_BETWEEN_LOW_HIGH( WDES, WDES%NB_TOTK(NK,ISP), LscaAWARE_LOCAL, CHAM )

          IFAIL=0

#if defined(MPI)
          IF ( LJACOBI .AND. IFLAG /=13 ) THEN
             IF (IU0>=0) WRITE(IU0,*)'jacobi called'
             CALL jacDSSYEV(WDES%COMM_KIN, CHAM(1,1), R, NB_TOT)
             CALLMPI( M_sum_g(WDES%COMM_KIN, CHAM(1,1),NB_TOT*NB_TOT))
             CALLMPI( M_sum_g(WDES%COMM_KIN, R , NB_TOT))

             GOTO 100
          ENDIF

          ! use scaLAPACK if available in parallel version
          IF (__IF_ACC_OFF__(LscaLAPACK.AND.IFLAG/=5)) THEN
             IF (.NOT. LscaAWARE_LOCAL) THEN
                CALL pDSSYEX_ZHEEVX(WDES%COMM_KIN, CHAM(1,1), R,  NB_TOT, WDES%NB_TOTK(NK,ISP))
!  alternatively use divide and conquer, somewhat slower though
!                CALL pDSYEV_ZHEEVD(WDES%COMM_KIN, CHAM(1,1), R,  NB_TOT, WDES%NB_TOTK(NK,ISP))
                CALLMPI( M_sum_g(WDES%COMM_KIN, CHAM(1,1),NB_TOT*NB_TOT))
             ELSE
!                ALLOCATE(COVL(SCALA_NP(),SCALA_NQ())); COVL=CHAM
!                DO N=1,20000; CHAM=COVL;  CALL BG_pDSSYEX_ZHEEVX(WDES%COMM_KIN, CHAM(1,1), R,  WDES%NB_TOTK(NK,ISP));  ENDDO
!                STOP

                CALL BG_pDSSYEX_ZHEEVX(WDES%COMM_KIN, CHAM(1,1), R,  WDES%NB_TOTK(NK,ISP))
!  alternatively use divide and conquer, somewhat slower though
!               CALL BG_pDSSYED_ZHEEVD(WDES%COMM_KIN, CHAM(1,1), R,  WDES%NB_TOTK(NK,ISP))
             ENDIF
             stotmg
             GOTO 100
          ENDIF
#endif
          !
          !  seriell codes
          !
#ifdef  gammareal
          IF (IFLAG == 5) THEN
             CALL __DSYGV__ &
                  (1,'V','U',WDES%NB_TOTK(NK,ISP),CHAM(1,1),NB_TOT,COVL(1,1),NB_TOT, &
                  R,CWRK,LWORK*NB_TOT, IFAIL)
          ELSE
#ifndef USE_ZHEEVX
             CALL __DSYEV__ &
                  ('V','U',WDES%NB_TOTK(NK,ISP),CHAM(1,1),NB_TOT, &
                  R,CWRK,LWORK*NB_TOT, IFAIL)
#else
             ABSTOL=1E-10_q
             VL=0 ; VU=0 ; IL=0 ; IU=0
             ALLOCATE(COVL(NB_TOT,NB_TOT))
!$ACC ENTER DATA CREATE(COVL) __IF_ASYNC__
             CALL __DSYEVX__ &
                  ( 'V', 'A', 'U', WDES%NB_TOTK(NK,ISP), CHAM(1,1) , NB_TOT, VL, VU, IL, IU, &
                  ABSTOL , NB_CALC , R, COVL(1,1), NB_TOT, CWRK, &
                  LWORK*NB_TOT, RWORK, IWORK, INFO, IFAIL)
             CALL __DCOPY__(NB_TOT*NB_TOT,COVL,1,CHAM,1)
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC EXIT DATA DELETE(COVL) IF(ACC_EXEC_ON)
             DEALLOCATE(COVL)
#endif
          ENDIF
#else
          IF (IFLAG == 5) THEN
             CALL __ZHEGV__ &
                  (1,'V','U',WDES%NB_TOTK(NK,ISP),CHAM(1,1),NB_TOT,COVL(1,1),NB_TOT, &
                  R,CWRK,LWORK*NB_TOT, RWORK, IFAIL)
          ELSE
#ifndef USE_ZHEEVX
             CALL __ZHEEV__ &
                  ('V','U',WDES%NB_TOTK(NK,ISP),CHAM(1,1),NB_TOT, &
                  R,CWRK,LWORK*NB_TOT, RWORK, IFAIL)
#else
             ABSTOL=1E-10_q
             VL=0 ; VU=0 ; IL=0 ; IU=0
             ALLOCATE(COVL(NB_TOT,NB_TOT))
!$ACC ENTER DATA CREATE(COVL) __IF_ASYNC__
             CALL __ZHEEVX__ &
                  ( 'V', 'A', 'U', WDES%NB_TOTK(NK,ISP), CHAM(1,1) , NB_TOT, VL, VU, IL, IU, &
                  ABSTOL , NB_CALC , R, COVL(1,1), NB_TOT, CWRK, &
                  LWORK*NB_TOT, RWORK, IWORK, INFO, IFAIL)
             CALL __ZCOPY__(NB_TOT*NB_TOT,COVL,1,CHAM,1)
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC EXIT DATA DELETE(COVL) IF(ACC_EXEC_ON)
             DEALLOCATE(COVL)
#endif
          ENDIF
#endif
          ! T3D uses a global sum which does not guarantee to give the same results on all nodes
          ! the following line is required to make the code waterproof (we had problems)
          ! since we now use a propritary sum (see mpi.F) we should not require
          ! this broadcast anymore
          ! stotmg
          ! CALLMPI( M_bcast_g(WDES%COMM, CHAM(1,1), NB_TOT*NB_TOT))

100       CONTINUE
          stotmg

          IF (IFAIL/=0) THEN
             CALL vtutor%error("ERROR in EDDIAG: call to ZHEEV/ZHEEVX/DSYEV/DSYEVX failed! error code &
                &was " // str(IFAIL))
          ENDIF

          ! shift eigenvalues back
          IF (IFLAG==23) CALL SHIFT_UNOCCUPIED_BACK(WDES%NB_TOTK(NK,ISP), R, W%FERTOT(:,NK, ISP))
          CALL SHIFT_BACK_BANDS_BETWEEN_LOW_HIGH(WDES, WDES%NB_TOTK(NK,ISP), R )

          DO N=1,WDES%NB_TOTK(NK,ISP)
             W%CELTOT(N,NK,ISP)=R(N)
          ENDDO
       ENDIF
!=======================================================================
! IFLAG > 2
! rotate wavefunctions
!=======================================================================
       IF (IFLAG==2) GOTO 1000

       IF (.NOT. LscaAWARE_LOCAL) THEN
          CALL LINCOM('F',WA%CW_RED(:,:),WA%CPROJ_RED(:,:),CHAM(1,1), &
            WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), & 
            WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
            WA%CW_RED(:,:),WA%CPROJ_RED(:,:))
       ELSE
          CALL LINCOM_DISTRI('F',WA%CW_RED(1,1),WA%CPROJ_RED(1,1),CHAM(1,1), &
            WDES%NB_TOTK(NK,ISP), & 
            WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,NB_TOT, &
            WDES%COMM_KIN, NBLK )
       ENDIF

1000   CONTINUE
       
       !  back redistribution over bands
       IF (WDES%DO_REDIS) THEN
          CALL REDISTRIBUTE_PROJ( ELEMENTS( W, WDES1, ISP))
          IF (LASYNC) THEN
             W%OVER_BAND=.TRUE.
          ELSE
             CALL REDISTRIBUTE_PW( ELEMENTS( W, WDES1, ISP))
          ENDIF
          DWRITE "redis ok"
       ENDIF

       ! return Hamiltonian (IFLAG==1) or rotation matrix (IFLAG=3,4)
       IF (PRESENT(LLAST)) THEN       
          IF (LLAST) THEN
!$ACC KERNELS PRESENT(CHAMHF,CHAM) __IF_ASYNC__
             IF (IFLAG==1) CHAM(:,:)=GCONJG(CHAM(:,:))
             CHAMHF(:,:,NK,ISP)=CHAM(:,:)
!$ACC END KERNELS
          ENDIF
       ENDIF

    ENDDO kpoint
    ENDDO spin

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)

!    CALLMPI( M_sum_z(WDES%COMM_KIN,CDCHF,1))
!    CALLMPI( M_sum_z(WDES%COMM_KINTER,CDCHF,1))
!    EXHF=CDCHF
!
!    IF (PRESENT(EXHF_ACFDT)) THEN
!       CALLMPI( M_sum_d(WDES%COMM_KIN,EXHF_ACFDT,1))
!       CALLMPI( M_sum_d(WDES%COMM_KINTER,EXHF_ACFDT,1))
!    ENDIF

    ! need to correct CELTOT at this point so that it is correct on all nodes
    IF (IFLAG==0) THEN
        CALL MRG_CEL(WDES,W)
    ENDIF
#ifdef MPI
    IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
       CALL KPAR_SYNC_CELTOT(WDES,W)
    ENDIF

    IF (IFLAG>=3.AND.(LHFCALC.OR.LUSEPEAD()).AND.WDES%COMM_KINTER%NCPU.GT.1) THEN
       CALL KPAR_SYNC_FERTOT(WDES,W)
       CALL KPAR_SYNC_CELTOT(WDES,W)
       CALL KPAR_SYNC_AUXTOT(WDES,W)
       CALL KPAR_SYNC_WAVEFUNCTIONS(WDES,W)
    ENDIF
#endif
!$ACC UPDATE DEVICE(W%CELTOT) __IF_ASYNC__

    ! deallocation ...
!$ACC WAIT IF(ACC_EXEC_ON)
!$ACC EXIT DATA DELETE(CHAM) IF((.NOT.LscaAWARE_LOCAL).AND.ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
    DEALLOCATE(CHAM)

    CALL DELWAV_R(W1)
!$ACC EXIT DATA DELETE(W1) __IF_ASYNC__

    CALL DELWAVA(WHAM)
    CALL DELWAVA_PROJ(WNONL)
!$ACC EXIT DATA DELETE(WHAM,WNONL) __IF_ASYNC__

    IF (USEFOCK_CONTRIBUTION().AND.(.NOT.PRESENT(CHAMHF))) THEN 
#ifdef _OPENACC
       CALL ACC_DELETE_TYPED_VAR(WFOCK)
#endif
       CALL DEALLOCW(WFOCK)
       IF (LPEAD_NO_SCF()) THEN
          DO N=1,NSIM_LOCAL
             CALL DELWAV(WTMP(N), .FALSE.)
          ENDDO
       ENDIF
    ENDIF

    IF (IFLAG==5) THEN
!$ACC EXIT DATA DELETE(COVL) __IF_ASYNC__
       DEALLOCATE(COVL)
       CALL DELWAVA_PROJ(WOVL)
!$ACC EXIT DATA DELETE(WOVL) __IF_ASYNC__
    ENDIF

    IF (PRESENT(CHAMHF).AND.PRESENT(LFIRST)) THEN
       LFIRST=.FALSE.
    ENDIF

#ifdef _OPENACC
    CALL ACC_DELETE_TYPED_VAR(WDES1)
#endif

#ifdef _OPENACC_TST
    CALL ACC_DELETE_TYPED_VAR(WDES)
    CALL ACC_COPYOUT_TYPED_VAR(W)
    CALL ACC_DELETE_TYPED_VAR(P)
    CALL ACC_DELETE_TYPED_VAR(GRID)
    CALL ACC_DELETE_TYPED_VAR(NONL_S)
    CALL ACC_DELETE_TYPED_VAR(NONLR_S)

!$ACC EXIT DATA DELETE(SV,CDIJ,CQIJ,LATT_CUR) __IF_ASYNC__
!$ACC EXIT DATA COPYOUT(CHAMHF) IF(ACC_EXEC_ON.AND.PRESENT(CHAMHF)) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA DELETE(HAMILTONIAN%MU,HAMILTONIAN) IF(ASSOCIATED(HAMILTONIAN%MU).AND.ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA DELETE(HAMILTONIAN%AVEC,HAMILTONIAN) IF(ASSOCIATED(HAMILTONIAN%AVEC).AND.ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)

!$ACC WAIT IF(ACC_EXEC_ON)
    ACC_ASYNC_Q = ACC_ASYNC_SYNC

    POP_ACC_EXEC_ON
#endif

    PROFILING_STOP('eddiag')

    RETURN
    
    CONTAINS

!***********************************************************************
!
! The following subroutines are used to constrain the
! diagonalization to the occupied manifold
! here orbitals with an occupancy close to 1
!
!***********************************************************************

    SUBROUTINE RESTRICT_TO_OCCUPIED_ONLY(NB_TOT, CHAM, FERTOT)
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      INTEGER :: NB_TOT
      GDEF    :: CHAM(:,:)
      REAL(q) :: FERTOT(:)
      INTEGER :: NB_OCC, I, J

      ! seek first occupancy that differs from 1.00
      DO NB_OCC=1, NB_TOT
         IF (ABS((FERTOT(NB_OCC))-1.0_q)>1E-5_q) EXIT
      ENDDO

      ! now set the lower (upper) triangle starting with row NB_OCC to zero
      ! loop over row index
!$ACC PARALLEL LOOP GANG PRESENT(CHAM) __IF_ASYNC__
      DO I=NB_OCC, NB_TOT
         ! loop over column index
!$ACC LOOP VECTOR
         DO J=1,I-1
            CHAM(I,J)=0
            CHAM(J,I)=0
         ENDDO
         ! shift diagonal elements by 10 eV
         ! so that the diagonalization does not mix occupied and unoccupied manyfold
         ! this shift needs to be removed by  SHIFT_UNOCCUPIED_BACK
         CHAM(I,I)=CHAM(I,I)+10.0_q
      ENDDO

    END SUBROUTINE RESTRICT_TO_OCCUPIED_ONLY

    SUBROUTINE SHIFT_UNOCCUPIED_BACK(NB_TOT, R, FERTOT)
      INTEGER :: NB_TOT
      REAL(q) :: R(:)
      REAL(q) :: FERTOT(:)
      INTEGER :: NB_OCC, I

      ! seek first occupancy that differs from 1.00
      DO NB_OCC=1, NB_TOT
         IF (ABS((FERTOT(NB_OCC)-1.0_q))>1E-5_q) EXIT
      ENDDO

      ! now set the lower (upper) triangle starting with row NB_OCC to zero
      ! loop over row index
      DO I=NB_OCC, NB_TOT
         R(I)=R(I)-10.0_q
      ENDDO
    END SUBROUTINE SHIFT_UNOCCUPIED_BACK

  END SUBROUTINE EDDIAG

!***********************************************************************
!
!  for a selected k-points and spin component
! ) calculate onsite (one centre) contribution to the Hamilton matrix
!  if LOVERL is .TRUE.
!
!    CCORR(m,n) += <psi_k,n| beta_i> D_ij <beta_j | psi_k,m> 
!
! ) if the local potential SV is passed  down the Hamilton matrix
!
!    CCORR(m,n) += <psi_k,n|  SV | psi_k,m> 
!
!  is added as well
!  the calling routine must initialize CCORR to zero
!  
!***********************************************************************

  SUBROUTINE ONE_CENTER_BETWEEN_STATES(HAMILTONIAN, LATT_CUR, LOVERL, WDES, W, NK, ISP, LMDIM, &
       CDIJC, CCORR, SV)
    USE prec
    USE wave_high
    USE dfast
    USE lattice
    USE hamil_struct_def

    TYPE (ham_handle)  HAMILTONIAN
    TYPE (latt)        LATT_CUR
    LOGICAL LOVERL
    INTEGER NK, ISP, LMDIM
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    OVERLAP CDIJC(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)  ! one centre correction
    GDEF ::  CCORR(WDES%NB_TOT,WDES%NB_TOT)
    RGRID, OPTIONAL  :: SV(DIMREAL(WDES%GRID%MPLWV),WDES%NCDIJ)   ! local potential
  ! local
    GDEF   , POINTER :: CPROJ_RED(:,:)
    INTEGER NPOS
    INTEGER NSTRIP, NSTRIP_ACT, NPOS_RED, NSTRIP_RED
    TYPE (wavefuna)    WNONL          ! array to hold non local part D * wave function character
    TYPE (wavefuna)    WOVL           ! array to hold non local part
    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1             ! current wavefunction
    TYPE (wavefuna)    WHAM           ! store Hamiltonian times wavefunction
    TYPE (wavefuna)    WA             ! array to store wavefunction

    NSTRIP=NSTRIP_STANDARD
!-----------------------------------------------------------------------
! non local part
!-----------------------------------------------------------------------
    IF (LOVERL) THEN
       CALL SETWDES(WDES,WDES1,NK)

       CALL NEWWAVA_PROJ(WNONL, WDES1)
       IF (WDES%DO_REDIS) THEN
          CALL NEWWAVA_PROJ(WOVL, WDES1)
          CALL WA_COPY_CPROJ(ELEMENTS(W, WDES1, ISP), WOVL)
       ELSE
          WOVL=ELEMENTS(W, WDES1, ISP)
       ENDIF
       
       CALL OVERL(WDES1, .TRUE., LMDIM, CDIJC(1,1,1,ISP), WOVL%CPROJ(1,1), WNONL%CPROJ(1,1))
       
       IF (WDES%DO_REDIS) THEN
          CALL REDIS_PROJ(WDES1, WDES%NBANDS, WNONL%CPROJ(1,1))
          CALL REDIS_PROJ(WDES1, WDES%NBANDS, WOVL%CPROJ (1,1))
       ENDIF
       
       DO NPOS=1,WDES%NBANDS,NSTRIP
          NSTRIP_ACT=MIN(WDES%NBANDS+1-NPOS,NSTRIP)
          NPOS_RED  =(NPOS-1)*WDES%NB_PAR+1
          NSTRIP_RED=NSTRIP_ACT*WDES%NB_PAR
          
          CALL ORTH1('U', &
               W%CW(1,1,NK,ISP),W%CW(1,1,NK,ISP),WOVL%CPROJ(1,1), &
               WNONL%CPROJ_RED(1,NPOS_RED),WDES%NB_TOT, &
               NPOS_RED, NSTRIP_RED, 0,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CCORR(1,1))
!             attention              - 

          
       ENDDO
       CALL DELWAVA_PROJ(WNONL)
       IF (WDES%DO_REDIS) CALL DELWAVA_PROJ(WOVL)
    ENDIF
!-----------------------------------------------------------------------
! local part
!-----------------------------------------------------------------------
    IF (PRESENT(SV)) THEN
       ! allocate work space
       ALLOCATE(W1%CR(WDES%GRID%MPLWV*WDES%NRSPINORS))

       CALL SETWDES(WDES,WDES1,NK)
       CALL NEWWAVA(WHAM, WDES1, NSTRIP)

       WA=ELEMENTS(W, WDES1, ISP)

       strip: DO NPOS=1,WDES%NBANDS,NSTRIP
          NSTRIP_ACT=MIN(WDES%NBANDS+1-NPOS,NSTRIP)

          !  calculate V_{local} |phi> + T | phi >
          !  for a block containing NSTRIP wavefunctions
          DO N=NPOS,NPOS+NSTRIP_ACT-1
             NP=N-NPOS+1

             CALL SETWAV(W, W1, WDES1, N, ISP)
             CALL FFTWAV_W1(W1)
             IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
                CALL HAMILT_LOCAL_TAU(W1, SV, LATT_CUR, HAMILTONIAN%MU, ISP,  WHAM%CW(:,NP), .FALSE., .FALSE.)
             ELSE
                CALL HAMILT_LOCAL(W1, SV, ISP,  WHAM%CW(:,NP), .FALSE., .FALSE.)
             ENDIF
          ENDDO
          ! redistribute wavefunctions
          ! after this redistributed up to and including 1...NPOS+NSTRIP_ACT
          IF (WDES%DO_REDIS) THEN
             CALL REDISTRIBUTE_PW( ELEMENTS( WA, NPOS, NPOS-1+NSTRIP_ACT))
             CALL REDISTRIBUTE_PW( ELEMENTS( WHAM, 1, NSTRIP_ACT))
          ENDIF

          NPOS_RED  =(NPOS-1)*WDES%NB_PAR+1
          NSTRIP_RED=NSTRIP_ACT*WDES%NB_PAR

          CALL ORTH1('U', &
               WA%CW_RED(1,1),WHAM%CW(1,1),WA%CPROJ_RED(1,1), &
               WA%CPROJ_RED(1,NPOS_RED),WDES%NB_TOT, &
               NPOS_RED, NSTRIP_RED, WDES1%NPL_RED,0,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CCORR(1,1))
          !           attention                   ---

       ENDDO strip

       IF (WDES%DO_REDIS) CALL REDISTRIBUTE_PW( ELEMENTS( W, WDES1, ISP))

       ! deallocation ...
       DEALLOCATE(W1%CR)
       CALL DELWAVA(WHAM)

    ENDIF

    CALLMPI( M_sum_g(WDES%COMM_KIN,CCORR(1,1),WDES%NB_TOT*WDES%NB_TOT))

  END SUBROUTINE ONE_CENTER_BETWEEN_STATES


!***********************************************************************
!
! read in the file GAMMA and add the density matrix or Hamiltonian
! to the diagonal density matrix F (one-electron occupancies) or
! Hamiltonian (one-electron eigenvalues)
! then diagonalize the resulting matrix and rotate the
! one electron orbitals accordingly
! 
!> @details @ref openacc :
!> On entering this subroutine W%CW and W%CPROJ are updated on the host
!> (if present on the device) and then OpenACC execution is switched off.
!> On exit the OpenACC execution mode reverts to its original status,
!> and when mopenacc_struct_dev::openacc_exec_on = .TRUE. the arrays
!> W%CW and W%CPROJ are updated on the device (if present).
!
!***********************************************************************

  SUBROUTINE ADD_GAMMA_FROM_FILE( WDES, W, KPOINTS, NELECT, NUP_DOWN, LWRITE, IO )
#ifdef _OPENACC
    USE mopenacc
#endif
    USE prec
    USE wave_high
    USE mkpoints
    USE dfast
    USE base
    USE fileio
    IMPLICIT NONE
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    REAL (q)  NELECT, NUP_DOWN        ! number of electrons, and total spin
    TYPE (kpoints_struct) KPOINTS     ! k-points structure
    LOGICAL :: LWRITE                 ! write final eigenvalues to OUTCAR
    TYPE (in_struct)   IO
  ! local
    INTEGER NK, ISP, NB, LMDIM
    TYPE (wavedes1)    WDES1
    TYPE (wavefuna)    WA             ! array to store wavefunction
    GDEF, ALLOCATABLE :: CHAM(:,:)
    REAL(q), ALLOCATABLE :: CELTOT(:,:,:)
    INTEGER :: IMODE                  ! add to Hamiltonian (IMODE=2) or density matrix (IMODE=1)
  ! LAPACK    
    REAL(q)    R(WDES%NB_TOT)
    INTEGER :: IFAIL
    INTEGER, PARAMETER :: LWORK=32
    GDEF       CWRK(LWORK*WDES%NB_TOT)
    REAL(q)    RWORK(3*WDES%NB_TOT)

    PROFILING_START('add_gamma_from_file')

!$ACC UPDATE SELF(W%CW,W%CPROJ) IF_PRESENT IF(ACC_EXEC_ON)
    PUSH_ACC_EXEC_ON(.FALSE.)

    CALL OPENGAMMA
    CALL READGAMMA_HEAD( IMODE, WDES%NKPTS , WDES%NB_TOT, IO)

    ALLOCATE(CHAM(WDES%NB_TOT, WDES%NB_TOT))

  ! for IMODE = 2 we save the CELTOT locally
    IF (IMODE==2) THEN
       ALLOCATE(CELTOT(SIZE(W%CELTOT,1),SIZE(W%CELTOT,2),SIZE(W%CELTOT,3)))
       CELTOT=W%CELTOT
    ENDIF


    DO ISP=1,WDES%ISPIN
    DO NK=1,WDES%NKPTS

       CHAM=0
       CALL READGAMMA(NK, WDES%NB_TOTK(NK,ISP), SIZE(CHAM,1), CHAM,  IO)
#ifdef MPI
       IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
       ! add one-particle occupancies
       IF (IMODE==1) THEN
          DO NB=1, WDES%NB_TOTK(NK,ISP)
             CHAM(NB,NB)=CHAM(NB,NB)+W%FERTOT(NB,NK,ISP)
          ENDDO
       ! add one-particle eigenvalues
       ELSE
          DO NB=1, WDES%NB_TOTK(NK,ISP)
             CHAM(NB,NB)=CHAM(NB,NB)+W%CELTOT(NB,NK,ISP)
          ENDDO
       ENDIF
#ifdef debug
       IF (IO%IU0>=0) CALL DUMP_HAM( "Gamma",WDES, CHAM)
#endif
       ! diagonalize the one-particle matrix (or Hamiltonian matrix)
       ! do not use the X routines, iterative diagonalization is not save
       ! for density matrix
       IFAIL=0
       
       ! change sign of density matrix to sort occupied states as lowest states
       IF (IMODE==1) CHAM=-CHAM
#ifdef testgamma_matrix
       IF (NK==1 .AND. ISP==1) THEN
          WRITE(*,'(8F14.7)') (REAL(CHAM(NB,NB),q),NB=1, W%WDES%NB_TOT)
       ENDIF
#endif
#ifdef  gammareal
       CALL DSYEV &
            ('V','U',WDES%NB_TOTK(NK,ISP),CHAM(1,1),WDES%NB_TOT, &
            R,CWRK,LWORK*WDES%NB_TOT, IFAIL)
#else
       CALL ZHEEV &
            ('V','U',WDES%NB_TOTK(NK,ISP),CHAM(1,1),WDES%NB_TOT, &
            R,CWRK,LWORK*WDES%NB_TOT, RWORK,  IFAIL)
#endif
#ifdef debug
       IF (IO%IU0>=0) CALL DUMP_HAM( "U",WDES, CHAM)
#endif
       ! change sign of eigenvalues
       IF (IMODE==1) R=-R

       IF (IFAIL/=0) THEN
          CALL vtutor%error("ERROR in ADD_GAMMA_FROM_FILE: call to ZHEEV/ DSYEV failed! error code was &
             &" // str(IFAIL))
       ENDIF

       CALL SETWDES(WDES,WDES1,NK)

       WA=ELEMENTS(W, WDES1, ISP)

       !  distribution over plane wave coefficients
       IF (WDES%DO_REDIS) CALL REDISTRIBUTE_PROJ( ELEMENTS( W, WDES1, ISP))
       IF (WDES%DO_REDIS) CALL REDISTRIBUTE_PW( ELEMENTS( W, WDES1, ISP))

       CALL LINCOM('F',WA%CW_RED(:,:),WA%CPROJ_RED(:,:),CHAM(1,1), &
            WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), & 
            WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,WDES%NB_TOT, &
            WA%CW_RED(:,:),WA%CPROJ_RED(:,:))

       !  back redistribution over bands
       IF (WDES%DO_REDIS) CALL REDISTRIBUTE_PROJ( ELEMENTS( W, WDES1, ISP))
       IF (WDES%DO_REDIS) CALL REDISTRIBUTE_PW( ELEMENTS( W, WDES1, ISP))

       ! updated one-electron occupancies
       IF (IMODE==1) THEN
          W%FERTOT(1:WDES%NB_TOTK(NK,ISP),NK,ISP)=R(1:WDES%NB_TOTK(NK,ISP))
       ELSE
          ! copy R to CELTOT
          W%CELTOT(1:WDES%NB_TOTK(NK,ISP),NK,ISP)=R(1:WDES%NB_TOTK(NK,ISP))
       ENDIF
    ENDDO
    ENDDO

    CALL CLOSEGAMMA

    ! just in case sync everything back to all nodes if KPAR is used
    CALL KPAR_SYNC_ALL(WDES,W)

    ! IMODE==2, update one-electron occupancies now
    IF (IMODE==2) THEN
       ! update fermi-weights
#ifdef  testgamma_matrix
       WRITE(*,'(8F14.7)') REAL(W%CELTOT(:,1,1),q)
#endif       
       CALL DENSTA_SIMPLE(W, KPOINTS, NELECT, NUP_DOWN )
       IF (LWRITE) THEN
          IF (IO%IU6>=0) THEN
             WRITE(IO%IU6,*) 'eigenvalues after inclusion of GAMMA file'
          ENDIF
          CALL WRITE_EIGENVAL( W%WDES, W, IO%IU6)
       ENDIF
       ! restore old one-electron eigenvalues
       W%CELTOT=CELTOT
    ENDIF

    DEALLOCATE(CHAM)

    POP_ACC_EXEC_ON
!$ACC UPDATE DEVICE(W%CW,W%CPROJ) IF_PRESENT IF(ACC_EXEC_ON)

  END SUBROUTINE ADD_GAMMA_FROM_FILE


!************************ SUBROUTINE EDDIAG_EXACT **********************
!
! this subroutine performs a full diagonalization of the Hamiltonian
! right now it is stupidly implemented since
! the number of bands is increased for all k-points
! than EDDIAG is called
! and finally the number of bands is set back to the original
! value
! this requires a lot of storage but is still convenient for 
! GW and RPA calculations
!
!> @details @ref openacc :
!> On entering this subroutine W%CW and W%CPROJ are updated on the host
!> (if present on the device) and then OpenACC execution is switched off.
!> On exit the OpenACC execution mode reverts to its original status,
!> and when mopenacc_struct_dev::openacc_exec_on = .TRUE. the arrays
!> W%CW and W%CPROJ are updated on the device (if present).
!
!***********************************************************************

  SUBROUTINE EDDIAG_EXACT(HAMILTONIAN, &
       GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
       LMDIM,CDIJ,CQIJ,IFLAG,SV,T_INFO,P,IU0,IU6,EXHF,EXHF_ACFDT)
#ifdef _OPENACC
    USE mopenacc
#endif
    USE prec
    USE wave_high
    USE lattice
    USE mpimy
    USE mgrid
    USE nonl_high
    USE hamil_struct_def
    USE main_mpi
    USE pseudo
    USE poscar
    USE ini
    USE choleski
    USE fock
    USE scala
    IMPLICIT NONE
    TYPE (ham_handle)  HAMILTONIAN
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W
    TYPE (wavedes)     WDES
    TYPE (symmetry) ::   SYMM      
    INTEGER LMDIM
    OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ),CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    INTEGER            IFLAG            ! determines mode of diagonalisation
    RGRID   SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ) ! local potential
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    INTEGER IU0, IU6
    REAL(q) EXHF
    REAL(q) EXHF_ACFDT
  ! local
    INTEGER NB_TOT    ! maximum number of plane wave coefficients = number of bands
    TYPE (wavedes)     WDES_TMP
    TYPE (wavespin)    W_TMP
    INTEGER NK, DEGREES_OF_FREEDOM
#ifdef _OPENACC
    LOGICAL, EXTERNAL :: USEFOCK_CONTRIBUTION
#endif

    PROFILING_START('eddiag_exact')

  ! just make sure that data distribution is over bands
    CALL REDIS_PW_OVER_BANDS(WDES, W)
  ! are all bands calculated anyway 
    DEGREES_OF_FREEDOM=MAXVAL(WDES%NPLWKP_TOT)
    IF (WDES%LGAMMA) THEN
       DEGREES_OF_FREEDOM=DEGREES_OF_FREEDOM*2-1
    ENDIF

    IF (DEGREES_OF_FREEDOM<=WDES%NB_TOT) THEN
       IFLAG=3
       CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
            LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IU0,EXHF,EXHF_ACFDT=EXHF_ACFDT)
    ELSE

       NB_TOT=((DEGREES_OF_FREEDOM+WDES%NB_PAR-1)/WDES%NB_PAR)*WDES%NB_PAR
    
       WDES_TMP=WDES
       WDES_TMP%NB_TOT=NB_TOT
       WDES_TMP%NBANDS=NB_TOT/WDES%NB_PAR
       CALL INIT_SCALAAWARE( WDES_TMP%NB_TOT, WDES_TMP%NRPLWV, WDES_TMP%COMM_KIN )
       
       NULLIFY(WDES_TMP%NB_TOTK)
       ALLOCATE(WDES_TMP%NB_TOTK(WDES%NKDIM,2))
! set the maximum number of bands k-point dependent
       DO NK=1,WDES_TMP%NKPTS
          IF (WDES_TMP%LGAMMA) THEN
             WDES_TMP%NB_TOTK(NK,:)=MIN(WDES_TMP%NB_TOT,WDES_TMP%NPLWKP_TOT(NK)*2-1)
          ELSE
             WDES_TMP%NB_TOTK(NK,:)=MIN(WDES_TMP%NB_TOT,WDES_TMP%NPLWKP_TOT(NK))
          ENDIF
       ENDDO
       CALL RESETUP_FOCK_WDES(WDES_TMP, LATT_CUR, LATT_CUR, -1)

       CALL ALLOCW(WDES_TMP,W_TMP)

#ifdef _OPENACC
       CALL ACC_COPYIN_TYPED_VAR(WDES_TMP)
!$ACC ENTER DATA COPYIN(W_TMP) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA CREATE(W_TMP%CPTWFP, &
!$ACC&                  W_TMP%CPROJ, &
!$ACC&                  W_TMP%FERTOT, &
!$ACC&                  W_TMP%FERWE, &
!$ACC&                  W_TMP%AUXTOT, &
!$ACC&                  W_TMP%CELTOT, &
!$ACC&                  W_TMP%CELEN) ASYNC(ACC_ASYNC_Q)
       IF (ASSOCIATED(W_TMP%WDES)) THEN
          CALL ACC_COPYIN_TYPED_VAR(W_TMP%WDES)
          CALL ACC_ATTACH_ASYNC(W_TMP%WDES,ACC_ASYNC_Q)
       ENDIF
#endif

       CALL DUMP_ALLOCATE(IU6)

       W_TMP%FERTOT=0
       W_TMP%CELTOT=0
    ! copy data back to work array
!$ACC KERNELS PRESENT(W,W_TMP) __IF_ASYNC__
       W_TMP%CW(:,1:WDES%NBANDS,:,:)    =W%CW(:,1:WDES%NBANDS,:,:)
       W_TMP%CPROJ(:,1:WDES%NBANDS,:,:) =W%CPROJ(:,1:WDES%NBANDS,:,:)
!$ACC END KERNELS
       W_TMP%CELTOT(1:WDES%NB_TOT,:,:)=W%CELTOT(1:WDES%NB_TOT,:,:)
       W_TMP%FERTOT(1:WDES%NB_TOT,:,:)=W%FERTOT(1:WDES%NB_TOT,:,:)

    ! random initialization beyond WDES%NBANDS
       CALL WFINIT(WDES_TMP, W_TMP, 1E10_q, WDES%NB_TOT+1) ! ENINI=1E10 not cutoff restriction
!$ACC UPDATE DEVICE(W_TMP%CW(:,WDES%NBANDS+1:,:,:),W_TMP%CPROJ(:,WDES%NBANDS+1:,:,:)) __IF_ASYNC__

    ! get characters   
       CALL PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,W_TMP)
    ! orthogonalization
       CALL ORTHCH(WDES_TMP,W_TMP, WDES%LOVERL, LMDIM,CQIJ)
    ! and diagonalization
       IFLAG=3
!note(sm): Updating FERTOT because FOCK_ALL_DBLBUF needs FERWE in sync on GPU and CPYU
!$ACC UPDATE DEVICE(W_TMP%FERTOT) IF(ACC_EXEC_ON .AND. USEFOCK_CONTRIBUTION()) ASYNC(ACC_ASYNC_Q)
       CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W_TMP,WDES_TMP,SYMM, &
            LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IU0,EXHF,EXHF_ACFDT=EXHF_ACFDT)
    ! copy data back to original array

!$ACC KERNELS PRESENT(W,W_TMP) __IF_ASYNC__
       W%CW(:,1:WDES%NBANDS,:,:)    =W_TMP%CW(:,1:WDES%NBANDS,:,:)
       W%CPROJ(:,1:WDES%NBANDS,:,:) =W_TMP%CPROJ(:,1:WDES%NBANDS,:,:)
!$ACC END KERNELS
       W%CELTOT(1:WDES%NB_TOT,:,:)=W_TMP%CELTOT(1:WDES%NB_TOT,:,:)
       W%FERTOT(1:WDES%NB_TOT,:,:)=W_TMP%FERTOT(1:WDES%NB_TOT,:,:)

#ifdef _OPENACC
       CALL ACC_DELETE_TYPED_VAR(W_TMP)
       CALL ACC_DELETE_TYPED_VAR(WDES_TMP)
#endif
       CALL DEALLOCW(W_TMP)
       DEALLOCATE(WDES_TMP%NB_TOTK)

       CALL RESETUP_FOCK_WDES(WDES, LATT_CUR, LATT_CUR, -1)

    ENDIF

    PROFILING_STOP('eddiag_exact')

  END SUBROUTINE EDDIAG_EXACT


!***********************************************************************
!
! The following subroutines are used to constrain the diagonalization
! to the set of orbitals between WDES%NBANDSLOW and WDES%NBANDSHIGH
!
!***********************************************************************

    SUBROUTINE SHIFT_BANDS_BETWEEN_LOW_HIGH(WDES, NB_TOT, LscaAWARE_LOCAL, CHAM )
      USE prec
      USE wave_high
      USE scala

      TYPE (wavedes)     WDES
      INTEGER :: NB_TOT
      LOGICAL :: LscaAWARE_LOCAL  ! use scaLAPACK
      GDEF    :: CHAM(:,:)        ! Hamiltonian (either distributed or global array)
      
      INTEGER :: I
      PROFILING_START( 'shift_bands_between_low_high' )

      ! first check whether we have to restrict the number of orbitals
      IF (WDES%NBANDSLOW > 0) THEN
         IF (LscaAWARE_LOCAL) THEN
#ifdef scaLAPACK
            CALL SCA_SHIFT_BANDS_LOW_NOINT( NB_TOT, WDES%NBANDSLOW, CHAM(1,1),-RSHIFT_NBLOW, DESCSTD)
#endif
         ELSE
            ! set off diagonal elements below NBANDSLOW to 0
!$ACC KERNELS PRESENT(WDES,CHAM) __IF_ASYNC__
            DO I=1, MIN(WDES%NBANDSLOW-1,NB_TOT)
               CHAM(I,I+1:)=0
               CHAM(I+1:,I)=0
               ! shift diagonal elements by - RSHIFT_NBLOW eV
               ! so that the diagonalization does not mix occupied and unoccupied manyfold
               ! this shift needs to be removed by  SHIFT_UNOCCUPIED_BACK
               CHAM(I,I)=CHAM(I,I) - RSHIFT_NBLOW
            ENDDO
!$ACC END KERNELS
         ENDIF
      ENDIF

      ! first check whether we have to restrict the number of orbitals
      IF (WDES%NBANDSHIGH > 0) THEN
         IF (LscaAWARE_LOCAL) THEN
            WRITE(*,*) 'shift high bands lapack',NB_TOT
#ifdef scaLAPACK
            CALL SCA_SHIFT_BANDS_HIGH_NOINT( NB_TOT, WDES%NBANDSHIGH, CHAM(1,1), RSHIFT_NBLOW, DESCSTD)
#endif
         ELSE
            ! set off diagonal elements below NBANDSLOW to 0
!$ACC KERNELS PRESENT(WDES,CHAM) __IF_ASYNC__
            DO I=WDES%NBANDSHIGH+1, NB_TOT
               CHAM(I,:I-1)=0
               CHAM(:I-1,I)=0
               ! shift diagonal elements by +100 eV
               ! so that the diagonalization does not mix occupied and unoccupied manyfold
               ! this shift needs to be removed by  SHIFT_UNOCCUPIED_BACK
               CHAM(I,I)=CHAM(I,I)+RSHIFT_NBLOW
            ENDDO
!$ACC END KERNELS
         ENDIF
      ENDIF

      PROFILING_STOP( 'shift_bands_between_low_high' )
    END SUBROUTINE SHIFT_BANDS_BETWEEN_LOW_HIGH


    SUBROUTINE SHIFT_BACK_BANDS_BETWEEN_LOW_HIGH(WDES, NB_TOT, R )
      USE prec
      USE wave_high

      TYPE (wavedes)     WDES
      INTEGER :: NB_TOT
      REAL(q) :: R(:)          ! eigenvalues
      
      INTEGER :: I

      ! first check whether we have to restrict the number of orbitals
      IF (WDES%NBANDSLOW > 0) THEN
         ! set off diagonal elements below NBANDSLOW to 0
         DO I=1, MIN(WDES%NBANDSLOW-1, NB_TOT)
            R(I)=R(I)+RSHIFT_NBLOW
         ENDDO
      ENDIF

      ! first check whether we have to restrict the number of orbitals
      IF (WDES%NBANDSHIGH > 0) THEN
         ! set off diagonal elements below NBANDSLOW to 0
         DO I=WDES%NBANDSHIGH+1, NB_TOT
            R(I)=R(I)-RSHIFT_NBLOW
         ENDDO
      ENDIF

    END SUBROUTINE SHIFT_BACK_BANDS_BETWEEN_LOW_HIGH

END MODULE subrot




!***********************************************************************
!
! dump a "Hamilton matrix" between the calculated states
!
!***********************************************************************

  
  SUBROUTINE DUMP_HAM( STRING, WDES, CHAM )
    USE wave
    CHARACTER (LEN=*) :: STRING
    TYPE (wavedes)     WDES
    GDEF ::  CHAM(WDES%NB_TOT,WDES%NB_TOT)
    INTEGER N1, N2, NPL2
    INTEGER NB_TOT

    NB_TOT=WDES%NB_TOT

    WRITE(*,*) STRING

    NPL2=MIN(12,NB_TOT)
    DO N1=1,NPL2
       WRITE(*,1)N1,(REAL( CHAM(N1,N2) ,KIND=q) ,N2=1,NPL2)
    ENDDO
    WRITE(*,*)
#ifndef gammareal
    DO N1=1,NPL2
       WRITE(6,2)N1,(AIMAG( CHAM(N1,N2)),N2=1,NPL2)
    ENDDO
    WRITE(*,*)
#endif

1   FORMAT(1I2,3X,40F9.5)
!1   FORMAT(1I2,3X,40F14.9)
2   FORMAT(1I2,3X,40F9.5)

  END SUBROUTINE DUMP_HAM


!***********************************************************************
!
! dump a "Hamilton matrix" between the calculated states
! single precision version
!
!***********************************************************************

  
  SUBROUTINE DUMP_HAM_SINGLE( STRING, WDES, CHAM)
    USE wave
    CHARACTER (LEN=*) :: STRING
    TYPE (wavedes)     WDES
    GDEFS ::  CHAM(WDES%NB_TOT,WDES%NB_TOT)
    INTEGER N1, N2, NPL2
    INTEGER NB_TOT

    NB_TOT=WDES%NB_TOT

    WRITE(*,*) STRING
    NPL2=MIN(10,NB_TOT)
    DO N1=1,NPL2
       WRITE(*,1)N1,(REAL( CHAM(N1,N2) ,KIND=q) ,N2=1,NPL2)
    ENDDO
    WRITE(*,*)
#ifndef gammareal
    DO N1=1,NPL2
       WRITE(6,2)N1,(AIMAG( CHAM(N1,N2)),N2=1,NPL2)
    ENDDO
    WRITE(*,*)
#endif
1   FORMAT(1I2,3X,40F9.5)
2   FORMAT(1I2,3X,40F9.5)
!2   FORMAT(1I2,3X,40E9.1)

  END SUBROUTINE DUMP_HAM_SINGLE


!***********************************************************************
!
! dump a "Hamilton matrix" between the calculated states
!
!***********************************************************************

  
  SUBROUTINE DUMP_HAM_SELECTED( STRING, WDES, CHAM, NDIM, NBANDS)
    USE wave
    CHARACTER (LEN=*) :: STRING
    TYPE (wavedes)     WDES
    INTEGER NDIM
    INTEGER NBANDS
    GDEF ::  CHAM(NDIM,NDIM)
  ! local
    INTEGER N1, N2, NPL2
    INTEGER NB_TOT

    NB_TOT=WDES%NB_TOT

    WRITE(*,*) STRING
    NPL2=MIN(12,NBANDS)
    DO N1=1,NPL2
       WRITE(*,1)N1,(REAL( CHAM(N1,N2) ,KIND=q) ,N2=1,NPL2)
    ENDDO
    WRITE(*,*)
#ifndef gammareal
    DO N1=1,NPL2
       WRITE(6,2)N1,(AIMAG( CHAM(N1,N2)),N2=1,NPL2)
    ENDDO
    WRITE(*,*)
#endif
1   FORMAT(1I2,3X,24F9.5)
2   FORMAT(1I2,3X,24F9.5)

  END SUBROUTINE DUMP_HAM_SELECTED


!=======================================================================
!
! small routine to dump a distributed matrix
! the descriptor in DESCA must properly describe the matrix
! since RECON_SLICE calls check 
!
!=======================================================================

#if defined(MPI) && defined(scaLAPACK)
  SUBROUTINE DUMP_HAM_DISTRI( STRING, WDES, CHAM_DISTRI, NB_TOT, DESCA, IU)
    USE wave
    USE scala
    IMPLICIT NONE
    CHARACTER (LEN=*) :: STRING            ! string to dump
    TYPE (wavedes)    :: WDES              ! wave function descriptor
    GDEF              :: CHAM_DISTRI(*)    ! distributed matrix
    INTEGER           :: NB_TOT
    INTEGER           :: DESCA(*)          ! distributed matrix descriptor array
    INTEGER           :: IU                ! unit to write to (not dump for IU<0)
  ! local
    INTEGER, PARAMETER :: NDUMP=16
    GDEF,ALLOCATABLE  :: CHAM(:,:)
    INTEGER N1, N2
    INTEGER COLUMN_HIGH, COLUMN_LOW

    COLUMN_LOW=1

    COLUMN_HIGH=MIN(COLUMN_LOW+NDUMP-1,NB_TOT)

    ALLOCATE(CHAM(NB_TOT, COLUMN_HIGH-COLUMN_LOW+1))

    CALL RECON_SLICE(CHAM, NB_TOT, NB_TOT, CHAM_DISTRI,  DESCA, COLUMN_LOW, COLUMN_HIGH)

    CALLMPI( M_sum_g(WDES%COMM_KIN, CHAM(1,1), SIZE(CHAM)))

    IF (IU>=0) THEN
    WRITE(IU,*) STRING
    DO N1=1,NDUMP
       WRITE(IU,1)N1+COLUMN_LOW-1,(REAL( CHAM(N1,N2+COLUMN_LOW-1) ,KIND=q) ,N2=1,NDUMP)
    ENDDO
    WRITE(IU,*)
#ifndef gammareal
    DO N1=1,NDUMP
       WRITE(IU,2)N1+COLUMN_LOW-1,(AIMAG( CHAM(N1,N2+COLUMN_LOW-1)),N2=1,NDUMP)
    ENDDO
    WRITE(IU,*)
#endif
    ENDIF

    DEALLOCATE(CHAM)
!1   FORMAT(1I2,3X,40F9.5)
!2   FORMAT(1I2,3X,40F9.5)
1   FORMAT(1I2,3X,40F7.4)
2   FORMAT(1I2,3X,40F7.4)
!1   FORMAT(1I2,3X,40F14.9)

  END SUBROUTINE DUMP_HAM_DISTRI
#endif

!************************ SUBROUTINE ORSP   ****************************
!
! this subroutine perfomes a gram-schmidt orthogonalistion of a set
! of vectors (all elements on local node)
! the subroutine uses BLAS 3 calls
!
!***********************************************************************

  SUBROUTINE ORSP(NBANDS, NPL, NRPLWV, CPTWFP)
#ifdef _OPENACC
    USE mopenacc
#endif
    USE prec
    USE tutor, ONLY: vtutor
    IMPLICIT NONE

    INTEGER NBANDS
    INTEGER NPL
    INTEGER NRPLWV
    GDEF CPTWFP(NRPLWV,NBANDS)
  ! local
    GDEF CPRO(NBANDS)
    REAL(q) WFMAG
    INTEGER I,N

    GDEF, EXTERNAL :: __GDOTC__

    IF (NBANDS> NRPLWV) THEN
       CALL vtutor%bug("internal error in ORSP: leading dimension of matrix too small", __FILE__, __LINE__)
    ENDIF

!$ACC DATA CREATE(CPRO) __IF_ASYNC__
!$ACC KERNELS PRESENT(CPRO) __IF_ASYNC__
    CPRO=0
!$ACC END KERNELS
    DO N=1,NBANDS

       ! normalise the vector

       WFMAG= __GDOTC__ (NPL,CPTWFP(1,N),1,CPTWFP(1,N),1)
       CALL   __GDSCAL__(NPL,1/SQRT(WFMAG),CPTWFP(1,N),1)

       ! now orthogonalise all higher vectors to the
       ! present vector

       IF (NBANDS/=N ) THEN
          CALL __GGEMV__( trans, NPL , NBANDS-N ,one , CPTWFP(1,N+1), &
               &             NRPLWV, CPTWFP(1,N), 1 , zero ,  CPRO, 1)

!$ACC PARALLEL LOOP PRESENT(CPRO) __IF_ASYNC__
          DO I=1,NBANDS
             CPRO(I)=GCONJG(CPRO(I))
          ENDDO

          CALL __GGEMM__( 'N', 'T' , NPL , NBANDS-N , 1 , -one , &
               &             CPTWFP(1,N), NRPLWV , CPRO , NBANDS , &
               &             one , CPTWFP(1,N+1) , NRPLWV )
       ENDIF
    ENDDO
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!$ACC END DATA
    RETURN
  END SUBROUTINE ORSP


!***********************************************************************
!
! use Loewdin perturbation to determine a rotation matrix
! this preserves the ordering of the eigenvalues
! MIND: does not work for real matrices
!
!***********************************************************************

  SUBROUTINE LOEWDIN_DIAG(NB_TOT, NBDIM, CHAM)
#ifdef _OPENACC
    USE mopenacc
#endif
    USE prec
    USE tutor, ONLY: vtutor
    IMPLICIT NONE
    INTEGER NB_TOT, NBDIM
    GDEF :: CHAM(NBDIM, NB_TOT)
  ! local
    REAL(q), PARAMETER  :: DIFMAX=0.001_q
    REAL(q) DIFCEL
    INTEGER N1, N2
    COMPLEX(q) :: CROT
    REAL(q) :: FAKT

    IF (NB_TOT>NBDIM) THEN
       CALL vtutor%bug("internal error in LOEWDIN_DIAG: leading dimension of matrix too small", __FILE__, __LINE__)
    ENDIF
!!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CHAM) PRIVATE(DIFCEL,CROT,FAKT) __IF_ASYNC__
!$ACC PARALLEL LOOP PRESENT(CHAM) PRIVATE(N1,DIFCEL,CROT,FAKT) __IF_ASYNC__
    DO N2=1,NB_TOT
       DO N1=1,N2-1
          DIFCEL= REAL( CHAM(N2,N2)-CHAM(N1,N1) ,KIND=q)
          IF (ABS(DIFCEL)<DIFMAX) THEN
             CROT  =0
          ELSE
             CROT  =GCONJG(CHAM(N1,N2))/DIFCEL
             IF (ABS(CROT)>0.1_q) THEN
                FAKT= 0.1_q/ABS(CROT)
                CROT  = CROT*FAKT
             ENDIF
          ENDIF
          CHAM(N2,N1) =-CROT
          CHAM(N1,N2) =-GCONJG(CROT)
       ENDDO
    ENDDO
!$ACC PARALLEL LOOP PRESENT(CHAM) __IF_ASYNC__
    DO N1=1,NB_TOT
       CHAM(N1,N1)=1
    ENDDO
  END SUBROUTINE LOEWDIN_DIAG


!*******************************************************************
!  calculate the matrix elements of a local potential
!           
!  CHAM(i,j)= <psi_i,k| V |psi_j,k> = int psi_i,k*(r) V(r) psi_j,k(r)
!
! between states
! the argument CVTOT must be in real space
! the result is retured in CHAM
!
!*******************************************************************

  SUBROUTINE LOCAL_BETWEEN_STATES( HAMILTONIAN, W, LATT_CUR, P, T_INFO, IRDMAX, LMDIM, &
       GRID_SOFT, GRIDC, GRIDUS, SOFT_TO_C, C_TO_US, CVTOT, CHAM)
    USE prec
    USE wave_high
    USE lattice
    USE poscar
    USE pseudo
    USE pot
    USE pawm
    USE subrot
    USE hamil_struct_def
    USE us
    USE tutor, ONLY: vtutor
    
    TYPE (ham_handle)  HAMILTONIAN
    TYPE (wavespin)    W
    
    INTEGER  IRDMAX      ! allocation required for augmentation
    TYPE (latt)        LATT_CUR
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    TYPE (grid_3d)     GRIDC                  ! grid for potentials/charge
    TYPE (grid_3d)     GRID_SOFT              ! grid for soft chargedensity
    TYPE (grid_3d)     GRIDUS                 ! grid for augmentation
    TYPE (transit)     SOFT_TO_C              ! index table between GRID_SOFT and GRIDC
    TYPE (transit)     C_TO_US                ! index table between GRID_SOFT and GRIDC
    COMPLEX(q)  CVTOT(GRIDC%MPLWV,W%WDES%NCDIJ) ! local potential
    INTEGER LMDIM
    GDEF       CHAM(W%WDES%NB_TOT,W%WDES%NB_TOT,W%WDES%NKPTS,W%WDES%ISPIN)
  ! local
    INTEGER ISP, NK
    RGRID ::   SV(DIMREAL(W%WDES%GRID%MPLWV),W%WDES%NCDIJ)   ! local potential
    OVERLAP :: CDIJ(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    OVERLAP :: CQIJ(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    INTEGER IRDMAA
    REAL(q)  DISPL(3,T_INFO%NIONS)

#ifdef MPI
    IF (W%WDES%COMM_KINTER%NCPU.NE.1) THEN
!PK Trivial but callers must be adapted
       CALL vtutor%error("LOCAL_BETWEEN_STATES: KPAR>1 not tested (but seems ok), sorry.")
    END IF
#endif

    DISPL=0

    ! get  the non local strenght parameters
    CALL SETDIJ_(W%WDES, GRIDC, GRIDUS, C_TO_US, LATT_CUR, P, T_INFO, W%WDES%LOVERL, &
         LMDIM, CDIJ, CQIJ, CVTOT, .FALSE., IRDMAA, IRDMAX, DISPL)

    ! transform CVTOT to reciprocal space (required by SET_SV)
    DO ISP=1,W%WDES%NCDIJ
       CALL FFT_RC_SCALE(CVTOT(1,ISP),CVTOT(1,ISP),GRIDC)
    ENDDO

    ! now set SV from CVTOT
    CALL SET_SV( W%WDES%GRID, GRIDC, GRID_SOFT, W%WDES%COMM_INTER, SOFT_TO_C, W%WDES%NCDIJ, SV, CVTOT)
    
    CHAM=0

    DO ISP=1,W%WDES%NCDIJ
       DO NK=1,W%WDES%NKPTS
          CALL ONE_CENTER_BETWEEN_STATES( HAMILTONIAN, LATT_CUR, W%WDES%LOVERL, W%WDES, W, NK, ISP, LMDIM, &
               CDIJ, CHAM(1,1,NK,ISP), SV)
!         CALL DUMP_HAM( "Hamiltonian", W%WDES, CHAM(1,1,NK,ISP))

       ENDDO
    ENDDO
       

  END SUBROUTINE LOCAL_BETWEEN_STATES
