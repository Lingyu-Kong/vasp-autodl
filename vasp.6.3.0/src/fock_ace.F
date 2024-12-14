#include "symbol.inc"
!***********************************************************************
!
! This module implements the Adaptively Compressed Exchange (ACE) of:
!
!***********************************************************************

      MODULE fock_ace
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE wave
      USE fock_glb, ONLY : LHFCALC,AEXX,LFOCKACE

      IMPLICIT NONE

      PUBLIC :: LFOCK_ACE,CXIJ,WACE,FOCK_ACE_ALLOC_CXIJ,FOCK_ACE_ALLOCW, &
                FOCK_ACE_DEALLOCW,FOCK_ACE_ACC_PAW,FOCK_ACE_CONSTRUCT,FOCK_ACE_ACC_APPLY

      PRIVATE

      ! stores the action of the Fock exchange operator on the orbitals
      TYPE (wavespin) :: WACE 
      ! and a logical to track whether WACE has been allocated
      LOGICAL :: LALLOCATED=.FALSE.

      ! stores the Fock related part of the PAW strength parameters
      ! (computed in SET_DD_PAW)
      OVERLAP , ALLOCATABLE :: CXIJ(:,:,:,:)

      CONTAINS

!************************ SUBROUTINE FOCK_ACE_CONSTRUCT ****************
!
! Compute the Fock matrix:
!
!   X_ij = < \tilde psi_i | \tilde V_X | \tilde psi_j >
!
! its Choleski decomposition:
!
!   X = -LL^T
!
! and the inverse of L, A=L^-1.
!
! From these compute the following linear combination of the action
! of the Fock exchange on the orbitals
!
!   | \tilde X_i > = \sum_j \tilde V_X | \psi_j > (A^T)_ji
!
! The adaptively compressed representation of the Fock exchange
! operator acting on the pseudo orbitals \tilde V_X is then given by
!
!   \tilde V_ACE = -\sum_i | \tilde X_i > < \tilde X_i |
!
! (see FOCK_ACE_ACC_APPLY)
!
! On exit, the set | \tilde X_i > will be stored in WACE
!
!***********************************************************************

      SUBROUTINE FOCK_ACE_CONSTRUCT(WDES,W,NKSTART,NKSTOP)

      USE dfast
      USE scala
      USE wave_high, ONLY : ELEMENTS,REDISTRIBUTE_PW
      USE tutor, ONLY: vtutor
      USE string, ONLY : str
      USE openmp, ONLY : omp_nthreads

      TYPE (wavedes)  :: WDES
      TYPE (wavespin) :: W
      INTEGER, OPTIONAL :: NKSTART,NKSTOP

      ! local variables
      TYPE (wavedes1) :: WDESK
      TYPE (wavefuna) :: WA,WHAM

      INTEGER :: NPOS,NSTRIP,NPOS_RED,NSTRIP_RED
      INTEGER :: NB_TOT,NBANDS
      INTEGER :: ISP,NK
      INTEGER :: INFO,NODE_ME,IONODE
      INTEGER :: MY_NKSTART, MY_NKSTOP

      GDEF, ALLOCATABLE, TARGET :: CHAM(:,:)

      PROFILING_START('fock_ace_construct')

#ifdef MPI
      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE
#else
      NODE_ME=0
      IONODE=0
#endif
      NBANDS=WDES%NBANDS
      NB_TOT=WDES%NB_TOT

      NSTRIP=NSTRIP_STANDARD
!$    NSTRIP=NSTRIP*omp_nthreads

      MY_NKSTART=1; IF (PRESENT(NKSTART)) MY_NKSTART = NKSTART
      MY_NKSTOP=WDES%NKPTS; IF (PRESENT(NKSTOP)) MY_NKSTOP = NKSTOP

      IF (.NOT.__IF_ACC_OFF__(LscaAWARE)) THEN
         ALLOCATE(CHAM(WDES%NB_TOT,WDES%NB_TOT))
!$ACC ENTER DATA CREATE(CHAM) IF(ACC_EXEC_ON)
      ELSE
         CALL INIT_scala(WDES%COMM_KIN,WDES%NB_TOT)
         ALLOCATE(CHAM(SCALA_NP(),SCALA_NQ()))
      ENDIF

!$ACC ENTER DATA CREATE(WDESK) IF(ACC_EXEC_ON)

      spn: DO ISP=1,WDES%ISPIN
      kpt: DO NK=MY_NKSTART,MY_NKSTOP
#ifdef MPI
         IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE kpt
#endif
         IF (COUNT(.NOT.isEmpty(W%FERTOT(:,NK,ISP)))==0) THEN
!$ACC KERNELS PRESENT(WACE%CW) __IF_ASYNC__
            WACE%CW(:,:,NK,ISP)=0
!$ACC END KERNELS
            WACE%CPROJ(:,:,NK,ISP)=0
            CYCLE kpt
         ENDIF

         IF (__IF_ACC_OFF__(LscaAWARE)) CALL INIT_scala(WDES%COMM_KIN,WDES%NB_TOTK(NK,ISP))
         CALL SETWDES(WDES,WDESK,NK)

         WA  =ELEMENTS(W   ,WDESK,ISP)
         WHAM=ELEMENTS(WACE,WDESK,ISP)

         ! redistribute to "over-plane-wave" distribution
         IF (WDES%DO_REDIS) THEN
            CALL REDISTRIBUTE_PW(WA)
            CALL REDISTRIBUTE_PW(WHAM)
         ENDIF
!$ACC KERNELS PRESENT(CHAM) __IF_ASYNC__
         CHAM=0
!$ACC END KERNELS
         DO NPOS=1,NBANDS,NSTRIP
            NPOS_RED  =                 (NPOS-1)*WDES%NB_PAR+1
            NSTRIP_RED=MIN(NSTRIP,NBANDS-NPOS+1)*WDES%NB_PAR

            IF (.NOT.__IF_ACC_OFF__(LscaAWARE)) THEN
               CALL ORTH1('U', &
                  WA%CW_RED(1,1),WHAM%CW_RED(1,NPOS_RED),WA%CPROJ_RED(1,1),WHAM%CPROJ_RED(1,NPOS_RED), &
                  NB_TOT,NPOS_RED,NSTRIP_RED,WDESK%NPL_RED,0,WDESK%NRPLWV_RED,0,CHAM(1,1))
            ELSE
               CALL ORTH1_DISTRI('U', &
                  WA%CW_RED(1,1),WHAM%CW_RED(1,NPOS_RED),WA%CPROJ_RED(1,1),WHAM%CPROJ_RED(1,NPOS_RED), &
                  NB_TOT,NPOS_RED,NSTRIP_RED,WDESK%NPL_RED,0,WDESK%NRPLWV_RED,0,CHAM(1,1), & 
                  WDES%COMM_KIN,WDES%NB_TOTK(NK,ISP))
            ENDIF
         ENDDO

         IF (.NOT.__IF_ACC_OFF__(LscaAWARE)) THEN
            CALLMPI( M_sum_g(WDES%COMM_KIN,CHAM(1,1),NB_TOT*NB_TOT))
         ENDIF

#ifdef debug
         io_begin
!$ACC UPDATE SELF(CHAM) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
         CALL DUMP_HAM("<i|X|j>",WDES,CHAM)
         io_end
#endif
         ! X_ij = < \psi_i | V_X | \psi_j > is negative semidefinite,
         ! i.e., multiply with -1 before the Choleski decomposition
!$ACC KERNELS PRESENT(CHAM) __IF_ASYNC__
         CHAM=-CHAM
!$ACC END KERNELS
!=======================================================================
! Choleski-decomposition of the Fock matrix + inversion of the result
! calling LAPACK-routines ZPOTRF (decomposition) and ZTRTRI (inversion):
!=======================================================================
         PROFILING_START('choleski')

         IF (.NOT.__IF_ACC_OFF__(LscaAWARE)) THEN
            IF (__IF_ACC_OFF__(LscaLAPACK.AND.LscaLU)) THEN
               CALL pPOTRF_TRTRI(WDES%COMM_KIN,CHAM(1,1),NB_TOT,WDES%NB_TOTK(NK,ISP))
               CALLMPI( M_sum_g(WDES%COMM_KIN,CHAM(1,1),NB_TOT*NB_TOT))
            ELSE
               INFO=0
#ifdef gammareal
               CALL __DPOTRF__ &
#else
               CALL __ZPOTRF__ &
#endif
                    & ('U',WDES%NB_TOTK(NK,ISP),CHAM(1,1),NB_TOT,INFO)
               IF (INFO/=0) THEN
                  CALL vtutor%error("FOCK_ACE_CONSTRUCT: LAPACK routine ZPOTRF failed!\n" // &
                                    "kpoint: " // str(NK) // " spin: " // str(ISP))
               ENDIF
#ifdef gammareal
               CALL __DTRTRI__ &
#else
               CALL __ZTRTRI__ &
#endif
                    & ('U','N',WDES%NB_TOTK(NK,ISP),CHAM(1,1),NB_TOT,INFO)
               IF (INFO/=0) THEN
                  CALL vtutor%error("FOCK_ACE_CONSTRUCT: LAPACK routine ZTRTRI failed!\n" // &
                                    "kpoint: " // str(NK) // " spin: " // str(ISP))
               ENDIF
            ENDIF
         ELSE
            CALL BG_pPOTRF_TRTRI(CHAM(1,1),WDES%NB_TOTK(NK,ISP),INFO)
            IF (INFO/=0) THEN
               CALL vtutor%error("FOCK_ACE_CONSTRUCT: scaLAPACK routine ZPOTRF ZTRTRI failed!\n" // &
                                 "kpoint: " // str(NK) // " spin: " // str(ISP))
            ENDIF
         ENDIF

         PROFILING_STOP('choleski')
#ifdef debug
         io_begin
!$ACC UPDATE SELF(CHAM) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
         CALL DUMP_HAM("choleski",WDES,CHAM)
         io_end
#endif
!=======================================================================
! construct the set | X_i> = \sum_j (L^T)^-1_ji V_X | \tilde \psi_j >
! CHAM(j,i) = (L^T)^-1_ji
!=======================================================================
         IF (.NOT.__IF_ACC_OFF__(LscaAWARE)) THEN
            CALL LINCOM('U', &
                 WHAM%CW_RED,WHAM%CPROJ_RED,CHAM(1,1), &
                 WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), &
                 WDESK%NPL_RED,0,WDESK%NRPLWV_RED,0,NB_TOT, &
                 WHAM%CW_RED,WHAM%CPROJ_RED)
         ELSE
            CALL LINCOM_DISTRI('U', &
                 WHAM%CW_RED(1,1),WHAM%CPROJ_RED(1,1),CHAM(1,1), &
                 WDES%NB_TOTK(NK,ISP),WDESK%NPL_RED,0,WDESK%NRPLWV_RED,0,NB_TOT, &
                 WDES%COMM_KIN,NBLK)
         ENDIF

         ! back distribution of the orbitals to default distribution
         ! N.B.: | \tilde X_i > remains in "over-plane-wave" distribution
         IF (WDES%DO_REDIS) THEN
            CALL REDISTRIBUTE_PW(WA)
         ENDIF

      ENDDO kpt
      ENDDO spn

#ifdef _OPENACC
      CALL ACC_DELETE_TYPED_VAR(WDESK)
!$ACC EXIT DATA DELETE(CHAM) __IF_ASYNC__
#endif
      DEALLOCATE(CHAM)

      PROFILING_STOP('fock_ace_construct')

      RETURN
      END SUBROUTINE FOCK_ACE_CONSTRUCT


!************************ SUBROUTINE FOCK_ACE_ACC_APPLY ****************
!
! Compute the action of the adaptively compressed representation of
! the Fock exchange operator on a set of orbitals
!
!   \tilde V_ACE |\tilde \psi_i > =
!       -\sum_j | \tilde X_j > < \tilde X_j | \tilde \psi_i >
!
!***********************************************************************

      SUBROUTINE FOCK_ACE_ACC_APPLY(W1,ISP,CACC,CDCHF)

      USE wave_mpi
      USE wave_high, ONLY : ELEMENTS

      ! orbitals onto which V_ACE should act, | \tilde \psi_i >
      TYPE (wavefun1) :: W1(:)
      ! spin component
      INTEGER :: ISP
      ! resulting action, V_ACE | \tilde \psi_i >
      COMPLEX(q), TARGET :: CACC(:,:)

      ! contribution to Fock double counting energy
      COMPLEX(q) :: CDCHF

      ! local variables
      TYPE (wavedes1), POINTER :: WDESK
      TYPE (wavefuna) :: WA

      COMPLEX(q), ALLOCATABLE, TARGET :: CWIN(:,:)
      COMPLEX(q), POINTER :: CW_RED(:,:)

      GDEF, ALLOCATABLE   :: COVL(:,:)

      REAL(q) :: WEIGHT,WEIGHTK

      INTEGER :: NSTRIP,NB_ACE,NB,NG,NGVECTOR,NSTRIP_RED

      PROFILING_START('fock_ace_acc_apply')

      NSTRIP=SIZE(W1)

      ! sanity check
      IF (SIZE(CACC,2)/=NSTRIP) THEN
         CALL vtutor%error("FOCK_ACE_ACC_APPLY: ERROR: dimensions do not match: " // str(NSTRIP) // &
            " " // str(SIZE(CACC,2)))
      ENDIF

      WDESK => W1(1)%WDES1

      ! copy W1(i)%CW(:) to CWIN(:,i)
      ALLOCATE(CWIN(WDESK%NRPLWV,NSTRIP))
!$ACC ENTER DATA CREATE(CWIN) __IF_ASYNC__

      NGVECTOR=WDESK%NGVECTOR*WDESK%NRSPINORS
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CWIN,W1) __IF_ASYNC__
NOACC !$OMP PARALLEL SHARED(NSTRIP,NGVECTOR,CWIN,W1) PRIVATE(NB,NG)
      DO NB=1,NSTRIP
NOACC !$OMP DO
         DO NG=1,NGVECTOR
            CWIN(NG,NB)=W1(NB)%CW(NG)
         ENDDO
NOACC !$OMP END DO
      ENDDO
NOACC !$OMP END PARALLEL
#ifdef MPI
      NSTRIP_RED=NSTRIP*WDESK%COMM_INTER%NCPU
#else
      NSTRIP_RED=NSTRIP
#endif
      ! redistribute to "over-plane-waves"
      IF (WDESK%DO_REDIS) THEN
         CALL SET_WPOINTER(CW_RED,WDESK%NRPLWV_RED,NSTRIP_RED,CWIN(1,1))
         CALL REDIS_PW(WDESK,NSTRIP,CWIN(1,1))
      ELSE
         CW_RED => CWIN
      ENDIF

      WA=ELEMENTS(WACE,WDESK,ISP)

      NB_ACE=SIZE(WA%CW_RED,2)

      ALLOCATE(COVL(NB_ACE,NSTRIP_RED))
!$ACC ENTER DATA CREATE(COVL) __IF_ASYNC__

      ! compute COVL(i,j)=< \tilde X_i | \tilde \psi_j >,
      ! for i \in [1,NB_ACE] and j \in [1,NSTRIP_RED]
      CALL __GGEMM__(trans,'N', NB_ACE, NSTRIP_RED, m_ WDESK%NPL_RED, &
     &   one, WA%CW_RED(1,1), m_ WACE%WDES%NRPLWV_RED, CW_RED(1,1), m_ WDESK%NRPLWV_RED, &
     &   zero, COVL(1,1), NB_ACE)

!$ACC EXIT DATA DELETE(CWIN) __IF_ASYNC__
      DEALLOCATE(CWIN)

      CALLMPI( M_sum_g(WACE%WDES%COMM_KIN,COVL(1,1),NB_ACE*NSTRIP_RED))

      ! CACC(:,j) = \sum_i < COVL(i,j) WA(:,i)
      CALL __GGEMM__('N','N', m_ WDESK%NPL_RED, NSTRIP_RED, NB_ACE, &
     &   one, WA%CW_RED(1,1), m_ WACE%WDES%NRPLWV_RED, COVL(1,1), NB_ACE, &
     &   zero, CACC(1,1), m_ WDESK%NRPLWV_RED)

!$ACC EXIT DATA DELETE(COVL) __IF_ASYNC__
      DEALLOCATE(COVL)

      ! back distribution
      IF (WDESK%DO_REDIS) THEN
         CALL REDIS_PW(WDESK,NSTRIP,CACC(1,1))
      ENDIF

      ! calculate contribution to Fock double counting energy
      CDCHF=0
!$ACC ENTER DATA COPYIN(CDCHF) __IF_ASYNC__

      WEIGHTK=0.5_q*WDESK%WTKPT*WDESK%RSPIN
!$ACC PARALLEL LOOP COLLAPSE(force:2) PRESENT(CACC,W1,CDCHF) __IF_ASYNC__ &
!$ACC PRIVATE(WEIGHT) REDUCTION(+:CDCHF)
NOACC !$OMP PARALLEL SHARED(NSTRIP,W1,WEIGHTK,NGVECTOR,CACC) PRIVATE(NB,NG,WEIGHT) REDUCTION(+:CDCHF)
      DO NB=1,NSTRIP
         WEIGHT=W1(NB)%FERWE*WEIGHTK
NOACC !$OMP DO
         DO NG=1,NGVECTOR
            CACC(NG,NB)=-CACC(NG,NB)
            CDCHF=CDCHF-CONJG(W1(NB)%CW(NG))*CACC(NG,NB)*WEIGHT
         ENDDO
NOACC !$OMP END DO
      ENDDO
NOACC !$OMP END PARALLEL

!$ACC EXIT DATA COPYOUT(CDCHF) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)

      NULLIFY(WDESK)

      PROFILING_STOP('fock_ace_acc_apply')

      RETURN
      END SUBROUTINE FOCK_ACE_ACC_APPLY


!************************ SUBROUTINE FOCK_ACE_ACC_PAW ******************
!
! Add the PAW contribution:
!
!   \sum_ij | p_i > CXIJ_ij < p_j | \tilde \psi_l >
!
! to the action of the Fock exchange operator on the pseudo orbitals.
!
!***********************************************************************

      SUBROUTINE FOCK_ACE_ACC_PAW(WDES,W,LATT_CUR,NONL_S,NONLR_S)

      USE lattice, ONLY : latt
      USE nonl_high

      TYPE (wavedes)  :: WDES
      TYPE (wavespin) :: W

      TYPE (latt) :: LATT_CUR

      TYPE (nonl_struct ) :: NONL_S
      TYPE (nonlr_struct) :: NONLR_S

      ! local variables
      TYPE (wavedes1) :: WDESK

      INTEGER :: ISP,NK,NB,ISPINOR

      GDEF,       ALLOCATABLE :: CTMP1(:)
      COMPLEX(q), ALLOCATABLE :: CTMP2(:)

      ! quick exit if possible
      IF (.NOT.WDES%LOVERL) RETURN

      ! sanity check, CXIJ has be computed already
      IF (.NOT.ALLOCATED(CXIJ)) THEN
         CALL vtutor%error("FOCK_ACE_ACC_PAW: ERROR: CXIJ not computed yet.")
      ENDIF

      ALLOCATE(CTMP1(WDES%NPROD),CTMP2(WDES%GRID%MPLWV*WDES%NRSPINORS))

      spn: DO ISP=1,WDES%ISPIN
      kpt: DO NK=1,WDES%NKPTS
#ifdef MPI
         IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE kpt
#endif
         IF (NONLR_S%LREAL) THEN
            CALL PHASER(WDES%GRID,LATT_CUR,NONLR_S,NK,WDES)
         ELSE
            CALL PHASE(WDES,NONL_S,NK)
         ENDIF

         CALL SETWDES(WDES,WDESK,NK)

!$OMP PARALLEL DO PRIVATE(CTMP1,CTMP2)
         DO NB=1,WDES%NBANDS

            ! compute CTMP(i)=\sum_j CXIJ_ij < p_j | \tilde \psi_{nb,nk,isp} >
            CALL OVERL1(WDESK,SIZE(CXIJ,1),CXIJ(1,1,1,1),CXIJ(1,1,1,1),0._q,W%CPROJ(1,NB,NK,ISP),CTMP1(1))
 
            IF (NONLR_S%LREAL) THEN
               CTMP2=0
               CALL RACC0(NONLR_S,WDESK,CTMP1(1),CTMP2(1))
               DO ISPINOR=0,WDESK%NRSPINORS-1
                  CALL FFTEXT(WDESK%NGVECTOR,WDESK%NINDPW(1), &
                 &   CTMP2(1+ISPINOR*WDESK%GRID%MPLWV),WACE%CW(1+ISPINOR*WDESK%NGVECTOR,NB,NK,ISP),WDESK%GRID,.TRUE.)
               ENDDO
            ELSE
               CALL VNLAC0(NONL_S,WDESK,CTMP1(1),WACE%CW(1,NB,NK,ISP))
            ENDIF

         ENDDO
!$OMP END PARALLEL DO
      ENDDO kpt
      ENDDO spn

      DEALLOCATE(CTMP1,CTMP2)

      RETURN
      END SUBROUTINE FOCK_ACE_ACC_PAW


!************************ SUBROUTINE FOCK_ACE_ALLOC_CXIJ ***************
!
! Allocate space to store the Fock exchange part of the PAW
! strength parameters (in SET_DD_PAW): CXIJ
!
!***********************************************************************

      SUBROUTINE FOCK_ACE_ALLOC_CXIJ(LMDIM,NIONS,NCDIJ)
      INTEGER :: LMDIM,NIONS,NCDIJ

      ! if CXIJ is already allocated then check its sizes
      IF (ALLOCATED(CXIJ)) THEN
         IF ((SIZE(CXIJ,1)/=LMDIM).OR.(SIZE(CXIJ,2)/=LMDIM).OR. &
        &    (SIZE(CXIJ,3)/=NIONS).OR.(SIZE(CXIJ,4)/=NCDIJ)) &
        &    DEALLOCATE(CXIJ)
      ENDIF
      ! (re-)allocate CXIJ
      IF (.NOT.ALLOCATED(CXIJ)) ALLOCATE(CXIJ(LMDIM,LMDIM,NIONS,NCDIJ))
      ! and set to zero
      CXIJ=0

      END SUBROUTINE FOCK_ACE_ALLOC_CXIJ


!************************ SUBROUTINE FOCK_ACE_DEALLOC_CXIJ *************
!
! Deallocate CXIJ
!
!***********************************************************************

      SUBROUTINE FOCK_ACE_DEALLOC_CXIJ

      IF (ALLOCATED(CXIJ)) DEALLOCATE(CXIJ)

      RETURN
      END SUBROUTINE FOCK_ACE_DEALLOC_CXIJ


!************************ SUBROUTINE FOCK_ACE_ALLOCW *******************
!
! Allocate space to store the action of the Fock exchange on the
! pseudo orbitals: WACE
!
!***********************************************************************

      SUBROUTINE FOCK_ACE_ALLOCW(WDES)
      TYPE (wavedes) :: WDES
      ! local variables
      TYPE (wavefun) :: WTMP

      IF (LALLOCATED) RETURN

      CALL ALLOCW(WDES,WACE,WTMP,WTMP)
#ifdef _OPENACC
!$ACC ENTER DATA CREATE(WACE   ) __IF_ASYNC__
!$ACC ENTER DATA CREATE(WACE%CW) __IF_ASYNC__
      CALL ACC_COPYIN_TYPED_VAR(WACE%WDES)
#endif
      LALLOCATED=.TRUE.

      RETURN
      END SUBROUTINE FOCK_ACE_ALLOCW


!************************ SUBROUTINE FOCK_ACE_DEALLOCW *****************
!
! Deallocate WACE
!
!***********************************************************************

      SUBROUTINE FOCK_ACE_DEALLOCW

      IF (.NOT.LALLOCATED) RETURN

#ifdef _OPENACC
      CALL ACC_DELETE_TYPED_VAR(WACE%WDES)
!$ACC EXIT DATA DELETE(WACE%CW) __IF_ASYNC__
!$ACC EXIT DATA DELETE(WACE   ) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
#endif
      CALL DEALLOCW(WACE)

      LALLOCATED=.FALSE.

      RETURN
      END SUBROUTINE FOCK_ACE_DEALLOCW


!************************ FUNCTION LFOCK_ACE ***************************
!
! This function returns .TRUE. when the ACE is active
!
!***********************************************************************

      FUNCTION LFOCK_ACE()
      LOGICAL :: LFOCK_ACE
#ifndef _OPENACC
      LFOCK_ACE=LHFCALC.AND.(AEXX/=0).AND. LFOCKACE
#else
      ! in case ACC_EXEC_ON=.TRUE. those code paths that can use
      ! the ACE actually have to use it (for now)
      LFOCK_ACE=LHFCALC.AND.(AEXX/=0).AND.(LFOCKACE.OR.ACC_EXEC_ON)
#endif
      END FUNCTION LFOCK_ACE

      END MODULE fock_ace
