#include "symbol.inc"

MODULE full_kpoints
!*************************************************************************
!
! this module contains routines to support the "full-koint" grid
! usually VASP stores only the kpoints in the IRBZ, however in some
! cases it is required to use the full k-point set (e.g. Hartree Fock)
! the implemented routines determine the full set from the set of k-points
! in the IRBZ by applying all symmetry operations of the cell
! additionally subroutines for determining wavefunction at the full
! grid from the wavefunctions at the IRBZ are supplied
!
!*************************************************************************

  USE prec
  USE mkpoints_struct_def

  IMPLICIT NONE
!*************************************************************************
!
! interfaces to low level F77 routines
!
!*************************************************************************

  INTERFACE  
     SUBROUTINE ROTATE_VECTOR_ADD(LINV,VEC,ROTVEC,MMAX,LMAX,SL,P,FAKT)
       USE prec
       USE pseudo
    
       LOGICAL LINV                        ! inversion required (i.e. conjugation of final vector)
       GDEFS:: VEC                         ! initial vector
       GDEF :: ROTVEC                      ! final vector
       INTEGER MMAX                        ! first and second dimension of SL
       INTEGER LMAX                        ! final dimension of SL
       REAL(q) :: SL                       ! rotation matrix (always symmetric)
       TYPE (potcar) P                     ! pseudopotential descriptor
       REAL(q) FAKT
     END SUBROUTINE ROTATE_VECTOR_ADD

  END INTERFACE

  INTERFACE
     SUBROUTINE ROTATE_VECTOR(LINV,VEC,ROTVEC,MMAX,LMAX,SL,NPS,LPS)
!$ACC ROUTINE VECTOR
       USE prec
       USE pseudo
       LOGICAL LINV                        ! inversion required (i.e. conjugation of final vector)
       GDEF :: VEC                         ! initial vector
       GDEF :: ROTVEC                      ! final vector
       INTEGER MMAX                        ! first and second dimension of SL
       INTEGER LMAX                        ! final dimension of SL
       REAL(q) :: SL                       ! rotation matrix (always symmetric)
       INTEGER :: NPS                      ! number of pseudopotential channels
       INTEGER :: LPS                      ! l-quantum number of the pseudopotential channels
     END SUBROUTINE ROTATE_VECTOR
  END INTERFACE

  INTERFACE
     SUBROUTINE ROTATE_WAVE_ADD( NPL, C_ROT, C, CPHASE, NINDPW, LINV, LSHIFT, FAKT)
       USE prec
       USE mpimy
       USE mgrid
       TYPE (grid_3d)     GRID
       COMPLEX(qs) :: C
       COMPLEX(q)  :: C_ROT,CPHASE
       INTEGER       NINDPW
       LOGICAL :: LINV, LSHIFT
       REAL(q) FAKT
       INTEGER NPL
     END SUBROUTINE ROTATE_WAVE_ADD
  END INTERFACE

  INTERFACE
     SUBROUTINE ROTATE_WAVE( NPL, C_ROT, C, CPHASE, NINDPW, LINV, LSHIFT)
       USE prec
       USE mpimy
       USE mgrid
       TYPE (grid_3d)     GRID
       COMPLEX(q) :: C,C_ROT,CPHASE
       INTEGER       NINDPW
       LOGICAL       LINV, LSHIFT
       INTEGER NPL
     END SUBROUTINE ROTATE_WAVE
  END INTERFACE

  INTERFACE
     SUBROUTINE ROTATE_WAVE_BACK( NPL, C_ROT, C, CPHASE, NINDPW, LINV, LSHIFT)
       USE prec
       USE mpimy
       USE mgrid
       TYPE (grid_3d)     GRID
       COMPLEX(q) :: C,C_ROT,CPHASE
       INTEGER       NINDPW
       LOGICAL       LINV, LSHIFT
       INTEGER NPL
     END SUBROUTINE ROTATE_WAVE_BACK
  END INTERFACE

CONTAINS

!*************************************************************************
!
! if this subroutine is called the full k-point grid is used
! and stored in KPOINTS_FULL
!
!*************************************************************************

      SUBROUTINE USE_FULL_KPOINTS
        LFULL_KPOINTS=.TRUE.
      END SUBROUTINE USE_FULL_KPOINTS

!*************************************************************************
!
! if this subroutine is called time inversion symmetry is preferentially
! used to determine orbitals at -k
!
!*************************************************************************

      SUBROUTINE USE_TIME_INVERSION
        LTIME_INVERSION=.TRUE.
      END SUBROUTINE USE_TIME_INVERSION

!*************************************************************************
!
! CHECK_FULL_KPOINTS can be used to make sure that all entries in the
! KPOINTS_FULL structure are properly set up
!
!*************************************************************************

      SUBROUTINE CHECK_FULL_KPOINTS
        USE tutor, ONLY: vtutor
#ifdef _OPENACC
        USE mopenacc
#endif
        IF (.NOT.LFULL_KPOINTS .OR. .NOT. ASSOCIATED(KPOINTS_FULL)) THEN
           CALL vtutor%bug("internal error in CHECK_FULL_KPOINTS: KPOINTS_FULL not properly initialised", __FILE__, __LINE__)
        ENDIF
#ifdef _OPENACC
        IF (ACC_EXEC_ON.AND.(.NOT.ACC_IS_PRESENT(KPOINTS_FULL,1))) CALL ACC_COPYIN_TYPED_VAR(KPOINTS_FULL)
#endif
      END SUBROUTINE CHECK_FULL_KPOINTS


!******************** SUBROUTINE IBZKPT_HF *******************************
!
! subroutine IBZKPT_HF is deduced from ibzkpt
! it constructs the full k-point mesh and connects every point 
! with its generating point in the ibz and stores the symmetry 
! transformation
!
! it returns:        the structure KPOINTS_F
!
! by Robin Hirschl, 092003
!***********************************************************************

      SUBROUTINE IBZKPT_HF(LATT_CUR, KPOINTS, KPOINTS_F, NIONS, ROTMAP, MAGROT, ISYM, IU6,IU0)
      USE lattice
      USE mkpoints
      USE main_mpi
      USE string, ONLY: str
      USE tutor, ONLY: vtutor, isAlert, KPOINTS_HF
#ifndef oldsym
      USE spinsym
#endif
      IMPLICIT NONE
! passed structures and variables
      TYPE (latt)           LATT_CUR 
      TYPE (kpoints_struct) KPOINTS
      TYPE (skpoints_full)  KPOINTS_F
      INTEGER ISYM
      INTEGER               IU6,IU0
      INTEGER NIONS
      INTEGER ROTMAP(:,:,:)
      REAL(q) MAGROT(:,:)
! common symmetry variables
      INTEGER ISYMOP, NROT, IGRPOP, NROTK, INVMAP, NPCELL
      REAL(q) GTRANS, AP
      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL
! local variables and structures
      TYPE (skpoints_full) K_TMP
      REAL(q) TINY, V(3), VR(3), VT(3), ROP(3,3), GSQU, WTKPT
      INTEGER INVERS(9),IOP(3,3),DET
      INTEGER NKTMP,NKPTF,NKPT,NKPT_ZERO_WEIGHT,NK,NKF,NK2,NOP,I,J
      INTEGER N1, N2, N3
      LOGICAL LINV, LEXIST
! external routine
      LOGICAL,EXTERNAL :: SGREQL
      LOGICAL, SAVE :: LWARN=.TRUE.  ! k-point warning only once
! set data
      DATA TINY /1.E-6_q/, INVERS /-1,0,0,0,-1,0,0,0,-1/

! warnings from tutor
      INTEGER,PARAMETER :: NTUTOR=1
      REAL(q)     RTUT(NTUTOR),RDUM
      INTEGER  ITUT(NTUTOR),IDUM
      COMPLEX(q)  CDUM  ; LOGICAL  LDUM

      INTEGER :: IHNF(3,3),NKPT_REGULAR
      REAL(q) :: AK(3,3),TMP1(3),TMP2(3)
      LOGICAL :: LSKIP_IRREGULAR_POINTS

! mM: 03062020: for the moment KPOINTS%NKPX/Y/Z = -1 for generalized regular grids ...
!     LSKIP_IRREGULAR_POINTS=(KPOINTS%NKPX>0.AND.KPOINTS%NKPY>0.AND.KPOINTS%NKPZ>0 .AND. .NOT.KPOINTS%LKBLOWUP)
      LSKIP_IRREGULAR_POINTS=(KPOINTS%MODE>ModeUnknown.AND.KPOINTS%MODE<LineMode .AND. .NOT.KPOINTS%LKBLOWUP)

      IF (LSKIP_IRREGULAR_POINTS) THEN
         CALL RECIPS(1._q,KPOINTS%B(1,1),KPOINTS%B(1,2),KPOINTS%B(1,3),AK(1,1),AK(1,2),AK(1,3))
         CALL vtutor%warning('full set of k-points will be restricted to regular points only.')
      ENDIF

      IF (NIONS/=SIZE(ROTMAP,1)) CALL vtutor%bug("NIONS does not match ROTMAP",__FILE__,__LINE__)

! allocate the temp array K_TMP, the absolute maximum possible number
! of full k-points is 96*NKPT
      NKPT=KPOINTS%NKPTS
      NKPTF=96*NKPT
      ALLOCATE(K_TMP%VKPT(3,NKPTF),K_TMP%WTKPT(NKPTF), K_TMP%TRANS(3,NKPTF), K_TMP%LSHIFT(NKPTF), & 
           K_TMP%SPINFLIP(NKPTF), &
           K_TMP%NEQUIV(NKPTF),K_TMP%IROTOP(3,3,NKPTF),K_TMP%ISYMOP(3,3,NKPTF), &
           K_TMP%LINV(NKPTF),K_TMP%NOP(NKPTF),K_TMP%ROTMAP(NIONS,NKPTF))
#ifndef oldsym
      ALLOCATE(K_TMP%RSSYMOP(2,2,NKPTF))
#endif
!=======================================================================
! first copy kpoints of IRBZ onto first positions of K_TMP
!=======================================================================
      NKPT_ZERO_WEIGHT=0
      DO NK=1,NKPT
         IF (KPOINTS%WTKPT(NK)==0) NKPT_ZERO_WEIGHT=NKPT_ZERO_WEIGHT+1
         K_TMP%VKPT(1,NK)=KPOINTS%VKPT(1,NK)
         K_TMP%VKPT(2,NK)=KPOINTS%VKPT(2,NK)
         K_TMP%VKPT(3,NK)=KPOINTS%VKPT(3,NK)
         K_TMP%WTKPT(NK) =KPOINTS%WTKPT(NK)
         K_TMP%NEQUIV(NK)=NK
         K_TMP%IROTOP(1:3,1:3,NK)=0
         K_TMP%IROTOP(1,1,NK)=1
         K_TMP%IROTOP(2,2,NK)=1
         K_TMP%IROTOP(3,3,NK)=1
         K_TMP%ISYMOP(1:3,1:3,NK)=0
         K_TMP%ISYMOP(1,1,NK)=1
         K_TMP%ISYMOP(2,2,NK)=1
         K_TMP%ISYMOP(3,3,NK)=1
         K_TMP%TRANS(:,NK)=0._q
         K_TMP%SPINFLIP(NK)=0
         DO I=1,NIONS
            K_TMP%ROTMAP(I,NK)=I
         ENDDO
         K_TMP%LINV=.FALSE.
#ifndef oldsym
         K_TMP%RSSYMOP(1,1,NK)=CMPLX(1.0_q,0.0_q,KIND=q)
         K_TMP%RSSYMOP(2,2,NK)=CMPLX(1.0_q,0.0_q,KIND=q)
         K_TMP%RSSYMOP(1,2,NK)=zero
         K_TMP%RSSYMOP(2,1,NK)=zero
#endif
      ENDDO
      NKTMP=NKPT
!=======================================================================
! now do all point group operations with each k-point and check wether we
! generate a new k-point (in 1st BZ). If so, store it as well as the
! generating operation. By the way, check whether inversion is already one 
! of the sym ops.
!=======================================================================
 inversion: IF (.NOT. LTIME_INVERSION) THEN
      LINV=.FALSE.
      DO NK=1,NKPT
         ! zero weight k-points are not included in exchange kernel
         IF (KPOINTS%WTKPT(NK)==0) CYCLE
         symop1: DO NOP=1,NROTK
! test existence of inversion op
            IF (SGREQL(IGRPOP(1,1,NOP),INVERS)) LINV=.TRUE.
! copy symmetry op to real array
            ROP=IGRPOP(:,:,NOP)
! make new k-point and shift (for testing) it to 1st BZ
            VR(1)=KPOINTS%VKPT(1,NK)
            VR(2)=KPOINTS%VKPT(2,NK)
            VR(3)=KPOINTS%VKPT(3,NK) 
            V(1)=VR(1)*ROP(1,1)+VR(2)*ROP(2,1)+VR(3)*ROP(3,1)
            V(2)=VR(1)*ROP(1,2)+VR(2)*ROP(2,2)+VR(3)*ROP(3,2)
            V(3)=VR(1)*ROP(1,3)+VR(2)*ROP(2,3)+VR(3)*ROP(3,3)
! bring the point to the primitive cell
            VT(1)=MOD(V(1)+6.5_q,1._q)-0.5_q
            VT(2)=MOD(V(2)+6.5_q,1._q)-0.5_q
            VT(3)=MOD(V(3)+6.5_q,1._q)-0.5_q
            IF (LSKIP_IRREGULAR_POINTS) THEN
               ! First bring VT to cartesian coordinates
               TMP1=MATMUL(LATT_CUR%B,VT)
               ! then express it in the basis KPOINTS%B
               TMP2=MATMUL(TRANSPOSE(AK),TMP1)
               ! Check whether this points is an element of the
               ! (possibly shifted) k-lattice generated by KPOINTS%B
               TMP1=TMP2-KPOINTS%SHIFT
               ! and skip this point if it is not
               IF (ANY(MOD(ABS(TMP1-NINT(TMP1)),1._q)>TINY)) CYCLE symop1
            ENDIF
! test against previous k-points
            LEXIST=.FALSE.
            test1: DO NKF=1,NKTMP
               ! zero weight k-points are excluded
               ! since they are most likely manually entered into the k-point list
               ! and then one would need to adjust their weight
               IF(( K_TMP%WTKPT(NKF)/=0) .AND. &
                  ( ABS(MOD(VT(1)-K_TMP%VKPT(1,NKF)+6.5,1._q)-0.5_q)<TINY) .AND. &
                  ( ABS(MOD(VT(2)-K_TMP%VKPT(2,NKF)+6.5,1._q)-0.5_q)<TINY) .AND. &
                  ( ABS(MOD(VT(3)-K_TMP%VKPT(3,NKF)+6.5,1._q)-0.5_q)<TINY)) THEN
                  LEXIST=.TRUE.
                  EXIT test1
               ENDIF
            ENDDO test1
! add if it does not exist yet
            IF (.NOT. LEXIST) THEN
               NKTMP=NKTMP+1
               IF (NKTMP>SIZE(K_TMP%VKPT,2)) CALL vtutor%bug("K_TMP is not sufficiently large",__FILE__,__LINE__)
               K_TMP%VKPT(1,NKTMP)=V(1)
               K_TMP%VKPT(2,NKTMP)=V(2)
               K_TMP%VKPT(3,NKTMP)=V(3)
               K_TMP%WTKPT(NKTMP) =KPOINTS%WTKPT(NK)
               K_TMP%NEQUIV(NKTMP)=NK
               K_TMP%LINV(NKTMP)=.FALSE.
               K_TMP%IROTOP(:,:,NKTMP)=INT(ROP(:,:))
               K_TMP%ISYMOP(:,:,NKTMP)=ISYMOP(:,:,NOP)
               K_TMP%TRANS(:,NKTMP)=GTRANS(:,NOP)
               K_TMP%SPINFLIP(NKTMP) =0
               IF ( SIZE(MAGROT)>1) THEN
                  IF( MAGROT(NOP,1)==-1) THEN
                     K_TMP%SPINFLIP(NKTMP) =1
                  ENDIF
               ENDIF
               K_TMP%ROTMAP(:,NKTMP)=ROTMAP(:,NOP,1)
               K_TMP%NOP(NKTMP)=NOP
#ifndef oldsym
               K_TMP%RSSYMOP(:,:,NKTMP)=RSSYMOP(:,:,NOP)
#endif
            ENDIF
         ENDDO symop1
      ENDDO
!=======================================================================
! did not find LINV -> now we have to do it all over again with 
! all operators multiplied with INVERS
!=======================================================================
      IF (.NOT. LINV .AND. ISYM>=0) THEN
         DO NK=1,NKPT
            ! zero weight k-points are not included in exchange kernel
            IF (KPOINTS%WTKPT(NK)==0) CYCLE
            symop2: DO NOP=1,NROTK
               ! apply inversion symmetry to form to get IOP
               CALL SGRPRD(INVERS,IGRPOP(1,1,NOP),IOP(1,1))
               ROP=IOP  ! copy symmetry op to real array
! make new k-point and shift it (for testing) to 1st BZ
               VR(1)=KPOINTS%VKPT(1,NK)
               VR(2)=KPOINTS%VKPT(2,NK)
               VR(3)=KPOINTS%VKPT(3,NK)
               V(1)=VR(1)*ROP(1,1)+VR(2)*ROP(2,1)+VR(3)*ROP(3,1)
               V(2)=VR(1)*ROP(1,2)+VR(2)*ROP(2,2)+VR(3)*ROP(3,2)
               V(3)=VR(1)*ROP(1,3)+VR(2)*ROP(2,3)+VR(3)*ROP(3,3)
! bring the point to the primitive cell
               VT(1)=MOD(V(1)+6.5_q,1._q)-0.5_q
               VT(2)=MOD(V(2)+6.5_q,1._q)-0.5_q
               VT(3)=MOD(V(3)+6.5_q,1._q)-0.5_q
               IF (LSKIP_IRREGULAR_POINTS) THEN
                  ! First bring VT to cartesian coordinates
                  TMP1=MATMUL(LATT_CUR%B,VT)
                  ! then express it in the basis KPOINTS%B
                  TMP2=MATMUL(TRANSPOSE(AK),TMP1)
                  ! Check whether this points is an element of the
                  ! (possibly shifted) k-lattice generated by KPOINTS%B
                  TMP1=TMP2-KPOINTS%SHIFT
                  ! and skip this point if it is not
                  IF (ANY(MOD(ABS(TMP1-NINT(TMP1)),1._q)>TINY)) CYCLE symop2
               ENDIF
! test against previous k-points
               LEXIST=.FALSE.
               test2: DO NKF=1,NKTMP
               IF(( K_TMP%WTKPT(NKF)/=0) .AND. &
                  ( ABS(MOD(VT(1)-K_TMP%VKPT(1,NKF)+6.5,1._q)-0.5)<TINY) .AND. &
                  ( ABS(MOD(VT(2)-K_TMP%VKPT(2,NKF)+6.5,1._q)-0.5)<TINY) .AND. &
                  ( ABS(MOD(VT(3)-K_TMP%VKPT(3,NKF)+6.5,1._q)-0.5)<TINY)) THEN
                     LEXIST=.TRUE.
                     EXIT test2
                  ENDIF
               ENDDO test2
! add if it does not exist yet
               IF (.NOT. LEXIST) THEN
                  NKTMP=NKTMP+1
                  IF (NKTMP>SIZE(K_TMP%VKPT,2)) CALL vtutor%bug("K_TMP is not sufficiently large",__FILE__,__LINE__)
                  K_TMP%VKPT(1,NKTMP)=V(1)
                  K_TMP%VKPT(2,NKTMP)=V(2)
                  K_TMP%VKPT(3,NKTMP)=V(3)
                  K_TMP%WTKPT(NKTMP) =KPOINTS%WTKPT(NK)
                  K_TMP%NEQUIV(NKTMP)=NK
                  K_TMP%LINV(NKTMP)=.TRUE.
                  K_TMP%IROTOP(:,:,NKTMP)=INT(ROP(:,:))
                  K_TMP%ISYMOP(:,:,NKTMP)=ISYMOP(:,:,NOP)
                  K_TMP%TRANS(:,NKTMP)=GTRANS(:,NOP)
                  K_TMP%SPINFLIP(NKTMP) =0
                  IF ( SIZE(MAGROT)>1) THEN
                     IF( MAGROT(NOP,1)==-1) THEN
                        K_TMP%SPINFLIP(NKTMP) =1
                     ENDIF
                  ENDIF
                  K_TMP%ROTMAP(:,NKTMP)=ROTMAP(:,NOP,1)
                  K_TMP%NOP(NKTMP)=NOP
#ifndef oldsym
                  K_TMP%RSSYMOP(:,:,NKTMP)=RSSYMOP(:,:,NOP)
#endif
               ENDIF
            ENDDO symop2
         ENDDO
      ENDIF
      ELSE inversion
!=======================================================================
! this version uses time inversion symmetry preferentially
!=======================================================================
      DO NK=1,NKPT
         ! zero weight k-points are not included in exchange kernel
         IF (KPOINTS%WTKPT(NK)==0) CYCLE
         symop3: DO NOP=1,NROTK
! copy symmetry op to real array
            ROP=IGRPOP(:,:,NOP)
! make new k-point and shift (for testing) it to 1st BZ
            VR(1)=KPOINTS%VKPT(1,NK)
            VR(2)=KPOINTS%VKPT(2,NK)
            VR(3)=KPOINTS%VKPT(3,NK) 
            V(1)=VR(1)*ROP(1,1)+VR(2)*ROP(2,1)+VR(3)*ROP(3,1)
            V(2)=VR(1)*ROP(1,2)+VR(2)*ROP(2,2)+VR(3)*ROP(3,2)
            V(3)=VR(1)*ROP(1,3)+VR(2)*ROP(2,3)+VR(3)*ROP(3,3)
! bring the point to the primitive cell
            VT(1)=MOD(V(1)+6.5_q,1._q)-0.5_q
            VT(2)=MOD(V(2)+6.5_q,1._q)-0.5_q
            VT(3)=MOD(V(3)+6.5_q,1._q)-0.5_q
            IF (LSKIP_IRREGULAR_POINTS) THEN
               ! First bring VT to cartesian coordinates
               TMP1=MATMUL(LATT_CUR%B,VT)
               ! then express it in the basis KPOINTS%B
               TMP2=MATMUL(TRANSPOSE(AK),TMP1)
               ! Check whether this points is an element of the
               ! (possibly shifted) k-lattice generated by KPOINTS%B
               TMP1=TMP2-KPOINTS%SHIFT
               ! and skip this point if it is not
               IF (ANY(MOD(ABS(TMP1-NINT(TMP1)),1._q)>TINY)) GOTO 1000
            ENDIF
! test against previous k-points
            LEXIST=.FALSE.
            test3: DO NKF=1,NKTMP
               ! zero weight k-points are excluded
               ! since they are most likely manually entered into the k-point list
               ! and then one would need to adjust their weight
               IF(( K_TMP%WTKPT(NKF)/=0) .AND. &
                  ( ABS(MOD(VT(1)-K_TMP%VKPT(1,NKF)+6.5,1._q)-0.5_q)<TINY) .AND. &
                  ( ABS(MOD(VT(2)-K_TMP%VKPT(2,NKF)+6.5,1._q)-0.5_q)<TINY) .AND. &
                  ( ABS(MOD(VT(3)-K_TMP%VKPT(3,NKF)+6.5,1._q)-0.5_q)<TINY)) THEN
                  LEXIST=.TRUE.
                  EXIT test3
               ENDIF
            ENDDO test3
! add if it does not exist yet
            IF (.NOT. LEXIST) THEN
               NKTMP=NKTMP+1
               IF (NKTMP>SIZE(K_TMP%VKPT,2)) CALL vtutor%bug("K_TMP is not sufficiently large",__FILE__,__LINE__)
               K_TMP%VKPT(1,NKTMP)=V(1)
               K_TMP%VKPT(2,NKTMP)=V(2)
               K_TMP%VKPT(3,NKTMP)=V(3)
               K_TMP%WTKPT(NKTMP) =KPOINTS%WTKPT(NK)
               K_TMP%NEQUIV(NKTMP)=NK
               K_TMP%LINV(NKTMP)=.FALSE.
               K_TMP%IROTOP(:,:,NKTMP)=INT(ROP(:,:))
               K_TMP%ISYMOP(:,:,NKTMP)=ISYMOP(:,:,NOP)
               K_TMP%TRANS(:,NKTMP)=GTRANS(:,NOP)
               K_TMP%SPINFLIP(NKTMP) =0
               IF ( SIZE(MAGROT)>1) THEN
                  IF( MAGROT(NOP,1)==-1) THEN
                     K_TMP%SPINFLIP(NKTMP) =1
                  ENDIF
               ENDIF
               K_TMP%ROTMAP(:,NKTMP)=ROTMAP(:,NOP,1)
               K_TMP%NOP(NKTMP)=NOP
#ifndef oldsym
               K_TMP%RSSYMOP(:,:,NKTMP)=RSSYMOP(:,:,NOP)
#endif
            ENDIF
            ! apply inversion symmetry to form to get IOP
 1000       CALL SGRPRD(INVERS,IGRPOP(1,1,NOP),IOP(1,1))
            ROP=IOP  ! copy symmetry op to real array
! make new k-point and shift it (for testing) to 1st BZ
            VR(1)=KPOINTS%VKPT(1,NK)
            VR(2)=KPOINTS%VKPT(2,NK)
            VR(3)=KPOINTS%VKPT(3,NK)
            V(1)=VR(1)*ROP(1,1)+VR(2)*ROP(2,1)+VR(3)*ROP(3,1)
            V(2)=VR(1)*ROP(1,2)+VR(2)*ROP(2,2)+VR(3)*ROP(3,2)
            V(3)=VR(1)*ROP(1,3)+VR(2)*ROP(2,3)+VR(3)*ROP(3,3)
! bring the point to the primitive cell
            VT(1)=MOD(V(1)+6.5_q,1._q)-0.5_q
            VT(2)=MOD(V(2)+6.5_q,1._q)-0.5_q
            VT(3)=MOD(V(3)+6.5_q,1._q)-0.5_q
            IF (LSKIP_IRREGULAR_POINTS) THEN
               ! First bring VT to cartesian coordinates
               TMP1=MATMUL(LATT_CUR%B,VT)
               ! then express it in the basis KPOINTS%B
               TMP2=MATMUL(TRANSPOSE(AK),TMP1)
               ! Check whether this points is an element of the
               ! (possibly shifted) k-lattice generated by KPOINTS%B
               TMP1=TMP2-KPOINTS%SHIFT
               ! and skip this point if it is not
               IF (ANY(MOD(ABS(TMP1-NINT(TMP1)),1._q)>TINY)) CYCLE symop3
            ENDIF
! test against previous k-points
            LEXIST=.FALSE.
            test4: DO NKF=1,NKTMP
            IF(( K_TMP%WTKPT(NKF)/=0) .AND. &
               ( ABS(MOD(VT(1)-K_TMP%VKPT(1,NKF)+6.5,1._q)-0.5)<TINY) .AND. &
               ( ABS(MOD(VT(2)-K_TMP%VKPT(2,NKF)+6.5,1._q)-0.5)<TINY) .AND. &
               ( ABS(MOD(VT(3)-K_TMP%VKPT(3,NKF)+6.5,1._q)-0.5)<TINY)) THEN
                  LEXIST=.TRUE.
                  EXIT test4
               ENDIF
            ENDDO test4
! add if it does not exist yet
            IF (.NOT. LEXIST) THEN
               NKTMP=NKTMP+1
               IF (NKTMP>SIZE(K_TMP%VKPT,2)) CALL vtutor%bug("K_TMP is not sufficiently large",__FILE__,__LINE__)
               K_TMP%VKPT(1,NKTMP)=V(1)
               K_TMP%VKPT(2,NKTMP)=V(2)
               K_TMP%VKPT(3,NKTMP)=V(3)
               K_TMP%WTKPT(NKTMP) =KPOINTS%WTKPT(NK)
               K_TMP%NEQUIV(NKTMP)=NK
               K_TMP%LINV(NKTMP)=.TRUE.
               K_TMP%IROTOP(:,:,NKTMP)=INT(ROP(:,:))
               K_TMP%ISYMOP(:,:,NKTMP)=ISYMOP(:,:,NOP)
               K_TMP%TRANS(:,NKTMP)=GTRANS(:,NOP)
               K_TMP%SPINFLIP(NKTMP) =0
               IF ( SIZE(MAGROT)>1) THEN
                  IF( MAGROT(NOP,1)==-1) THEN
                     K_TMP%SPINFLIP(NKTMP) =1
                  ENDIF
               ENDIF
               K_TMP%ROTMAP(:,NKTMP)=ROTMAP(:,NOP,1)
               K_TMP%NOP(NKTMP)=NOP
#ifndef oldsym
               K_TMP%RSSYMOP(:,:,NKTMP)=RSSYMOP(:,:,NOP)
#endif
            ENDIF
         ENDDO symop3
      ENDDO

      END IF inversion
!=======================================================================
! set the weights properly
!=======================================================================
      ! loop over all original k-points in the IRZ
      DO NK=1,NKPT
         ! loop over all k-points in the full BZ
         I=0  ! I counts the number of equivalent k-points
         DO NKF=1,NKTMP
            IF (K_TMP%NEQUIV(NKF)==NK) I=I+1
         ENDDO
         DO NKF=1,NKTMP
            ! set new weight to weight of k-point in the IRZ divided by the number
            ! of equivalent points
            IF (K_TMP%NEQUIV(NKF)==NK) K_TMP%WTKPT(NKF)=KPOINTS%WTKPT(NK)/I
         ENDDO
      ENDDO
!=======================================================================
! allocate KPOINTS_F and copy values
!=======================================================================
      NULLIFY(KPOINTS_F%PHASE, KPOINTS_F%NG_INDEX)
      ALLOCATE(KPOINTS_F%VKPT(3,NKTMP),KPOINTS_F%WTKPT(NKTMP),KPOINTS_F%TRANS(3,NKTMP), & 
           KPOINTS_F%LSHIFT(NKTMP),KPOINTS_F%SPINFLIP(NKTMP), KPOINTS_F%NOP(NKTMP), &
           KPOINTS_F%NEQUIV(NKTMP),KPOINTS_F%IROTOP(3,3,NKTMP),KPOINTS_F%ISYMOP(3,3,NKTMP), &
           KPOINTS_F%LINV(NKTMP),KPOINTS_F%ROTMAP(NIONS,NKTMP))

      NULLIFY(KPOINTS_F%MAP_INTO_BZ)
     
      KPOINTS_F%NKPTS=NKTMP
#ifndef oldsym
      ALLOCATE(KPOINTS_F%RSSYMOP(2,2,NKTMP))
#endif
      DO NK=1,NKTMP
         KPOINTS_F%VKPT(1,NK)=K_TMP%VKPT(1,NK)
         KPOINTS_F%VKPT(2,NK)=K_TMP%VKPT(2,NK)
         KPOINTS_F%VKPT(3,NK)=K_TMP%VKPT(3,NK)
         KPOINTS_F%WTKPT(NK) =K_TMP%WTKPT(NK)
! set weights all equal (uncomment these lines to recover the old behaviour)
!         IF (K_TMP%WTKPT(NK)/=0) THEN
!            KPOINTS_F%WTKPT(NK)=1._q/REAL(NKTMP-NKPT_ZERO_WEIGHT,q)
!         ELSE
!            KPOINTS_F%WTKPT(NK)=0
!         ENDIF
         KPOINTS_F%NEQUIV(NK)=K_TMP%NEQUIV(NK)
         KPOINTS_F%LINV(NK)=K_TMP%LINV(NK)
         KPOINTS_F%IROTOP(:,:,NK)=K_TMP%IROTOP(:,:,NK)
         KPOINTS_F%ISYMOP(:,:,NK)=K_TMP%ISYMOP(:,:,NK)
         KPOINTS_F%TRANS(:,NK)=K_TMP%TRANS(:,NK)
         KPOINTS_F%LSHIFT(NK)=K_TMP%LSHIFT(NK)
         KPOINTS_F%SPINFLIP(NK)=K_TMP%SPINFLIP(NK)
         KPOINTS_F%NOP(NK)=K_TMP%NOP(NK)
         KPOINTS_F%ROTMAP(:,NK)=K_TMP%ROTMAP(:,NK)
#ifndef oldsym
         KPOINTS_F%RSSYMOP(:,:,NK)=K_TMP%RSSYMOP(:,:,NK)
#endif
      ENDDO
      IF (ABS(SUM(KPOINTS_F%WTKPT(:))-1)>1E-6_q) THEN
         CALL vtutor%bug("The k-point weights do not sum to one \n the (routine has &
            &been updated recently \n check lines around 'set weights all equal')\n"&
            &// str(KPOINTS_F%WTKPT),__FILE__,__LINE__)
      ENDIF

      KPOINTS_F%NKPX=KPOINTS%NKPX
      KPOINTS_F%NKPY=KPOINTS%NKPY
      KPOINTS_F%NKPZ=KPOINTS%NKPZ

      KPOINTS_F%B    =KPOINTS%B
      KPOINTS_F%SHIFT=KPOINTS%SHIFT

      KPOINTS_F%MODE =KPOINTS%MODE

      KPOINTS_F%NKPTS_NON_ZERO=KPOINTS_F%NKPTS-NKPT_ZERO_WEIGHT

      ! In case we know the generating lattice we can check whether application
      ! of the symops. of the lattice has generated any irregular points
      IF (KPOINTS%MODE>ModeUnknown.AND.KPOINTS%MODE<LineMode) THEN
         ! First determine the number of points in the regular grid
         CALL HNFORM(LATT_CUR%B,KPOINTS%B,IHNF); NKPT_REGULAR=IHNF(1,1)*IHNF(2,2)*IHNF(3,3)
         ! If there are more than NKPT_REGULAR points with non-zero weight then
         ! the list of k-points no longer represents a purely regular grid w.r.t.
         ! KPOINTS%B (and KPOINTS%SHIFT)
         IF (COUNT(KPOINTS_F%WTKPT>0)/=NKPT_REGULAR .AND. LWARN) THEN
            LWARN=.FALSE. ; CALL vtutor%write(isAlert, KPOINTS_HF)
         ENDIF
      ENDIF
!=======================================================================
!
! determine minimum distance between two k-points
!
!=======================================================================
      KPOINTS_F%VKPTMINDIST2=1
      DO NK=1,NKTMP
         ! phase shift (required for space group operations)
         IF ((ABS(KPOINTS_F%TRANS(1,NK))>TINY) .OR. &
             (ABS(KPOINTS_F%TRANS(2,NK))>TINY) .OR. &
             (ABS(KPOINTS_F%TRANS(3,NK))>TINY)) THEN
            KPOINTS_F%LSHIFT(NK)=.TRUE.
         ELSE
            KPOINTS_F%LSHIFT(NK)=.FALSE.
         ENDIF

         IF (KPOINTS_F%WTKPT(NK)==0)  CYCLE
         ! for regular grids we do not have to do the following more than once
         IF (KPOINTS_F%NKPX>0 .AND. KPOINTS_F%NKPY>0 .AND. KPOINTS_F%NKPZ>0 .AND. NK>1) CYCLE
         DO NKF=NK,NKTMP
            IF (KPOINTS_F%WTKPT(NKF)==0)  CYCLE
            VT(1)=(KPOINTS_F%VKPT(1,NK)-KPOINTS_F%VKPT(1,NKF))*LATT_CUR%B(1,1)+ &
                  (KPOINTS_F%VKPT(2,NK)-KPOINTS_F%VKPT(2,NKF))*LATT_CUR%B(1,2)+ &
                  (KPOINTS_F%VKPT(3,NK)-KPOINTS_F%VKPT(3,NKF))*LATT_CUR%B(1,3)
            VT(2)=(KPOINTS_F%VKPT(1,NK)-KPOINTS_F%VKPT(1,NKF))*LATT_CUR%B(2,1)+ &
                  (KPOINTS_F%VKPT(2,NK)-KPOINTS_F%VKPT(2,NKF))*LATT_CUR%B(2,2)+ &
                  (KPOINTS_F%VKPT(3,NK)-KPOINTS_F%VKPT(3,NKF))*LATT_CUR%B(2,3)
            VT(3)=(KPOINTS_F%VKPT(1,NK)-KPOINTS_F%VKPT(1,NKF))*LATT_CUR%B(3,1)+ &
                  (KPOINTS_F%VKPT(2,NK)-KPOINTS_F%VKPT(2,NKF))*LATT_CUR%B(3,2)+ &
                  (KPOINTS_F%VKPT(3,NK)-KPOINTS_F%VKPT(3,NKF))*LATT_CUR%B(3,3)
            ! loop over all boxes
            DO N1=-1,1
            DO N2=-1,1
            DO N3=-1,1
               V(1)=N1*LATT_CUR%B(1,1)+N2*LATT_CUR%B(1,2)+N3*LATT_CUR%B(1,3)
               V(2)=N1*LATT_CUR%B(2,1)+N2*LATT_CUR%B(2,2)+N3*LATT_CUR%B(2,3)
               V(3)=N1*LATT_CUR%B(3,1)+N2*LATT_CUR%B(3,2)+N3*LATT_CUR%B(3,3)
               IF (N1==0 .AND. N2==0 .AND. N3==0 .AND. NKF==NK) CYCLE

               GSQU=(VT(1)+V(1))**2+(VT(2)+V(2))**2+(VT(3)+V(3))**2
               KPOINTS_F%VKPTMINDIST2=MIN(GSQU,KPOINTS_F%VKPTMINDIST2)
               IF (ABS(KPOINTS_F%VKPTMINDIST2)<1E-8) THEN
                  CALL vtutor%error("Two equivalent k points were found when trying to generate the "//&
                                    "full brillouin zone: " // str(NK) // " " // str(NKF) //". "//&
                                    "Please check that the input k points contains a regular grid.")
               ENDIF
            ENDDO
            ENDDO
            ENDDO
         ENDDO
      ENDDO
! do some output
      IF (IU6>=0) THEN
         WRITE(IU6,'(A)') ' '
         WRITE(IU6,'(A)') ' Subroutine IBZKPT_HF returns following result:'
         WRITE(IU6,'(A)') ' =============================================='
         WRITE(IU6,'(A)') ' '
         WRITE(IU6,'(A,I6,A)') ' Found ',NKTMP,' k-points in 1st BZ'
         WRITE(IU6,'(A,I6,A)') ' the following ',KPOINTS_F%NKPTS,' k-points will be used (e.g. in the exchange kernel)'
         WRITE(IU6,'(A)') ' Following reciprocal coordinates:   # in IRBZ'
         DO NK=1,KPOINTS_F%NKPTS
            IF (KPOINTS_F%SPINFLIP(NK)==1) THEN
               WRITE(IU6,'(3F10.6,4X,F10.8,I4," t-inv ",L1," spinflp")') &
                    KPOINTS_F%VKPT(1,NK), &
                    KPOINTS_F%VKPT(2,NK), &
                    KPOINTS_F%VKPT(3,NK),KPOINTS_F%WTKPT(NK),KPOINTS_F%NEQUIV(NK),KPOINTS_F%LINV(NK)
            ELSE
               WRITE(IU6,'(3F10.6,4X,F10.8,I4," t-inv ",L1)') &
                    KPOINTS_F%VKPT(1,NK), &
                    KPOINTS_F%VKPT(2,NK), &
                    KPOINTS_F%VKPT(3,NK),KPOINTS_F%WTKPT(NK),KPOINTS_F%NEQUIV(NK),KPOINTS_F%LINV(NK)
            ENDIF
         ENDDO
      ENDIF

! deallocate tmp arrays
      DEALLOCATE(K_TMP%VKPT, K_TMP%TRANS, K_TMP%LSHIFT, K_TMP%SPINFLIP, &
           K_TMP%NEQUIV, K_TMP%IROTOP, K_TMP%ISYMOP, &
           K_TMP%LINV, K_TMP%NOP, K_TMP%ROTMAP)

#ifndef oldsym
      DEALLOCATE(K_TMP%RSSYMOP)
#endif
      END SUBROUTINE IBZKPT_HF


!************************ SUBROUTINE SET_FULL_KPOINTS ******************
!
! this routine resets the NKPTS_FOR_GEN_LAYOUT and the  VKPT entries 
! in the WDES array to the full k-point grid
! this is required to make sure that the GRID%RC structure contains
! all G vectors in reciprocal space that are required to perform an
! FFT from reciprocal to real space if the number of k-point changes
! dynamically
!
!*********************************************************************** 

  SUBROUTINE SET_FULL_KPOINTS(NKPTS_FOR_GEN_LAYOUT,VKPT)
    USE prec
    USE tutor, ONLY: vtutor
    IMPLICIT NONE 
    INTEGER NKPTS_FOR_GEN_LAYOUT
    REAL(q), POINTER :: VKPT(:,:)

    IF (.NOT. LFULL_KPOINTS) RETURN
    IF (.NOT. ASSOCIATED(KPOINTS_FULL)) THEN
       CALL vtutor%bug("internal error in SET_FULL_KPOINTS: KPOINTS_FULL not properly initialised", __FILE__, __LINE__)
    ENDIF

    NKPTS_FOR_GEN_LAYOUT=KPOINTS_FULL%NKPTS
    VKPT=>KPOINTS_FULL%VKPT

  END SUBROUTINE SET_FULL_KPOINTS


!******************* SUBROUTINE SETUP_FULL_KPOINTS *********************
!
! sets up the full k-point mesh and the transformation matrix from
! the kpoints in the IRBZ to the full mesh
! this is routine simply calls IBZKPT_HF
!
!***********************************************************************

  SUBROUTINE SETUP_FULL_KPOINTS(KPOINTS,LATT_CUR,NIONS,ROTMAP,MAGROT,ISYM,LINVERSION,IU6,IU0, LSYMGRAD)
#ifdef _OPENACC
    USE mopenacc
#endif
    USE prec
    USE mkpoints
    USE lattice
    USE tutor, ONLY: vtutor
    
    IMPLICIT NONE 
    
    TYPE (kpoints_struct) KPOINTS
    TYPE (latt) LATT_CUR
    INTEGER NIONS
    INTEGER ROTMAP(:,:,:)
    REAL(q) MAGROT(:,:)
    INTEGER ISYM
    LOGICAL LINVERSION
    INTEGER IU6,IU0
    LOGICAL, OPTIONAL :: LSYMGRAD
    INTEGER NK,NE
    REAL(q) VT(3) ! debugging reasons only

    LSYMGRAD_SAVE=.FALSE.
    IF (PRESENT(LSYMGRAD)) THEN
       LSYMGRAD_SAVE=LSYMGRAD
    ENDIF

    IF (.NOT. LFULL_KPOINTS) RETURN
    
#ifdef _OPENACC
    IF (ASSOCIATED(KPOINTS_FULL).AND.ACC_IS_PRESENT(KPOINTS_FULL,1)) THEN
       PUSH_ACC_EXEC_ON(.TRUE.)
       CALL ACC_DELETE_TYPED_VAR(KPOINTS_FULL)
       POP_ACC_EXEC_ON
    ENDIF
#endif

    NULLIFY(KPOINTS_FULL)
    ALLOCATE(KPOINTS_FULL)

    !check if KPOINTS is in line mode in which case we stop the code with an error message
    IF (KPOINTS%MODE==LineMode) THEN
        call vtutor%error('This calculation requires a uniform k mesh spanning the entire Brillouin zone, '//&
                          'and cannot be run using k points generated in line mode (KPOINTS file).')
    ENDIF

    CALL IBZKPT_HF(LATT_CUR, KPOINTS, KPOINTS_FULL, NIONS, ROTMAP, MAGROT, ISYM, IU6, IU0)
#ifdef debug
    101 FORMAT(2I4,2X,7(F8.5,1X),L)
    WRITE(*,*) "*****************************************"
    DO NK=1,KPOINTS_FULL%NKPTS
       NE=KPOINTS_FULL%NEQUIV(NK)
       VT(1)=KPOINTS%VKPT(1,NE)*KPOINTS_FULL%IROTOP(1,1,NK)+ &
            KPOINTS%VKPT(2,NE)*KPOINTS_FULL%IROTOP(2,1,NK)+ &
            KPOINTS%VKPT(3,NE)*KPOINTS_FULL%IROTOP(3,1,NK)
       VT(2)=KPOINTS%VKPT(1,NE)*KPOINTS_FULL%IROTOP(1,2,NK)+ &
            KPOINTS%VKPT(2,NE)*KPOINTS_FULL%IROTOP(2,2,NK)+ &
            KPOINTS%VKPT(3,NE)*KPOINTS_FULL%IROTOP(3,2,NK)
       VT(3)=KPOINTS%VKPT(1,NE)*KPOINTS_FULL%IROTOP(1,3,NK)+ &
            KPOINTS%VKPT(2,NE)*KPOINTS_FULL%IROTOP(2,3,NK)+ &
            KPOINTS%VKPT(3,NE)*KPOINTS_FULL%IROTOP(3,3,NK)
       WRITE(*,101) NK,NE,KPOINTS_FULL%VKPT(1,NK),KPOINTS_FULL%VKPT(2,NK), &
            KPOINTS_FULL%VKPT(3,NK),KPOINTS_FULL%TRANS(1,NK),KPOINTS_FULL%TRANS(2,NK), &
            KPOINTS_FULL%TRANS(3,NK),KPOINTS_FULL%WTKPT(NK), &
            KPOINTS_FULL%LINV(NK)
    ENDDO
    WRITE(*,*) "*****************************************"
#endif

    CALL SMALL_SPACE_GROUP_WEIGHTS( KPOINTS, KPOINTS_FULL, ROTMAP, MAGROT,  ISYM, LINVERSION, IU6)

#ifdef _OPENACC
    CALL ACC_COPYIN_TYPED_VAR(KPOINTS_FULL)
#endif
    RETURN
  END SUBROUTINE SETUP_FULL_KPOINTS


!******************* SUBROUTINE DEALLOC_FULL_KOINTS_STATIC *************
!
! deallocate static data structure KPOINTS_FULL as well as all
! related entities
!
!***********************************************************************

  SUBROUTINE DEALLOC_FULL_KPOINTS_STATIC
#ifdef _OPENACC
    USE mopenacc
#endif
    IF (.NOT. LFULL_KPOINTS) RETURN

    IF (ASSOCIATED(KPOINTS_FULL)) THEN
#ifdef _OPENACC
      IF (ACC_IS_PRESENT(KPOINTS_FULL,1)) THEN
         PUSH_ACC_EXEC_ON(.TRUE.)
         CALL ACC_DELETE_TYPED_VAR(KPOINTS_FULL)
         POP_ACC_EXEC_ON
      ENDIF
#endif
      DEALLOCATE(KPOINTS_FULL%VKPT,KPOINTS_FULL%WTKPT, &
           KPOINTS_FULL%TRANS,KPOINTS_FULL%LSHIFT, KPOINTS_FULL%SPINFLIP, KPOINTS_FULL%NOP, &
           KPOINTS_FULL%NEQUIV,KPOINTS_FULL%IROTOP,KPOINTS_FULL%ISYMOP, &
           KPOINTS_FULL%LINV,KPOINTS_FULL%ROTMAP)

       IF (ASSOCIATED(KPOINTS_FULL%PHASE)) THEN
          DEALLOCATE(KPOINTS_FULL%PHASE)
       ENDIF
       IF (ASSOCIATED(KPOINTS_FULL%NG_INDEX)) THEN
          DEALLOCATE(KPOINTS_FULL%NG_INDEX)
       ENDIF
       IF (ASSOCIATED(KPOINTS_FULL%MAP_INTO_BZ)) THEN
          DEALLOCATE(KPOINTS_FULL%MAP_INTO_BZ)
       ENDIF
#ifndef oldsym
       DEALLOCATE(KPOINTS_FULL%RSSYMOP)
#endif
       DEALLOCATE(KPOINTS_FULL)
    ENDIF
  END SUBROUTINE DEALLOC_FULL_KPOINTS_STATIC

!******************* SUBROUTINE DEALLOC_FULL_KOINTS ********************
!
! deallocate data structure KPOINTS_FULL handled down by the calling
! routine 
!
!***********************************************************************

  SUBROUTINE DEALLOC_FULL_KPOINTS(KPOINTS_FULL)
#ifdef _OPENACC
    USE mopenacc
#endif
    TYPE (skpoints_full), POINTER :: KPOINTS_FULL

    IF (.NOT. LFULL_KPOINTS) RETURN

    ! gK: I am not quite sure why this is required
    ! does this need to go or stay, I guess I will leave it in
    IF (ASSOCIATED(KPOINTS_FULL)) THEN
#ifdef _OPENACC
      IF (ACC_IS_PRESENT(KPOINTS_FULL,1)) THEN
         PUSH_ACC_EXEC_ON(.TRUE.)
         CALL ACC_DELETE_TYPED_VAR(KPOINTS_FULL)
         POP_ACC_EXEC_ON
      ENDIF
#endif
      DEALLOCATE(KPOINTS_FULL%VKPT,KPOINTS_FULL%WTKPT, &
           KPOINTS_FULL%TRANS,KPOINTS_FULL%LSHIFT, KPOINTS_FULL%SPINFLIP, KPOINTS_FULL%NOP, &
           KPOINTS_FULL%NEQUIV,KPOINTS_FULL%IROTOP,KPOINTS_FULL%ISYMOP, &
           KPOINTS_FULL%LINV,KPOINTS_FULL%ROTMAP)

       IF (ASSOCIATED(KPOINTS_FULL%PHASE)) THEN
          DEALLOCATE(KPOINTS_FULL%PHASE)
       ENDIF
       IF (ASSOCIATED(KPOINTS_FULL%NG_INDEX)) THEN
          DEALLOCATE(KPOINTS_FULL%NG_INDEX)
       ENDIF
       IF (ASSOCIATED(KPOINTS_FULL%MAP_INTO_BZ)) THEN
          DEALLOCATE(KPOINTS_FULL%MAP_INTO_BZ)
       ENDIF
#ifndef oldsym
       DEALLOCATE(KPOINTS_FULL%RSSYMOP)
#endif
       DEALLOCATE(KPOINTS_FULL)
    ENDIF
  END SUBROUTINE DEALLOC_FULL_KPOINTS

!*********************************************************************
!
! CHECK_GEN_LAYOUT_GRID checks whether two GRID structure are
! identical
!
!*********************************************************************


  SUBROUTINE CHECK_GEN_LAYOUT_GRID( GRID, GRID_NEW )
    USE mgrid
    USE string, ONLY: str
    USE tutor, ONLY: vtutor
    
    TYPE (grid_3d)     GRID
    TYPE (grid_3d)     GRID_NEW
  ! local
    INTEGER IND

    IF (GRID%RC%NCOL/=GRID_NEW%RC%NCOL)  THEN
       CALL vtutor%bug("internal error in RE_GEN_LAYOUT_GRID: GRID%RC%NCOL/=GRID_NEW%RC%NCOL " // &
          str(GRID%RC%NCOL) // " " // str(GRID_NEW%RC%NCOL), __FILE__, __LINE__)
    ENDIF
    
    DO IND=1,GRID%RC%NCOL
       IF (GRID%RC%I2(IND)/=GRID_NEW%RC%I2(IND)) THEN
          CALL vtutor%bug("internal error in RE_GEN_LAYOUT_GRID: GRID I2 changed " // str(IND) // &
             " " // str(GRID%RC%I2(IND)) // " " // str(GRID_NEW%RC%I2(IND)), __FILE__, __LINE__)
       ENDIF
       IF (GRID%RC%I3(IND)/=GRID_NEW%RC%I3(IND)) THEN
          CALL vtutor%bug("internal error in RE_GEN_LAYOUT_GRID: GRID I3 changed " // str(IND) // &
             " " // str(GRID%RC%I3(IND)) // " " // str(GRID_NEW%RC%I3(IND)), __FILE__, __LINE__)
       ENDIF
    ENDDO
  END SUBROUTINE CHECK_GEN_LAYOUT_GRID

!*********************************************************************
!
! CHECK_GEN_LAYOUT checks whether two wave descriptors posses identical 
! plane waves and an identical data layout
! this is required if both descriptors are used for "compressed"
! FFTs using FFTWAV and FFTEXT
!
!*********************************************************************

  SUBROUTINE CHECK_GEN_LAYOUT( WDES, WDES_NEW, NK_TEST)
    USE prec
    USE wave
    USE string, ONLY: str
    USE tutor, ONLY: vtutor
    
    TYPE (wavedes)     WDES
    TYPE (wavedes)     WDES_NEW
    INTEGER NK_TEST
  ! local
    INTEGER IND, NK

    ! check the number of G-vectors for the first NK_TEST k-points
    DO NK=1,NK_TEST
       IF (WDES%NGVECTOR(NK)/=WDES_NEW%NGVECTOR(NK)) THEN
          WRITE(0,*) 'internal error in CHECK_GEN_LAYOUT: number of G-vectors changed',NK,WDES%NGVECTOR(NK),WDES_NEW%NGVECTOR(NK)
!          STOP
       ENDIF
       DO IND=1,WDES%NGVECTOR(NK)
          IF (WDES_NEW%IGX(IND,NK)/=WDES%IGX(IND,NK) .OR. &
              WDES_NEW%IGY(IND,NK)/=WDES%IGY(IND,NK) .OR. &
              WDES_NEW%IGZ(IND,NK)/=WDES%IGZ(IND,NK)) THEN
              CALL vtutor%bug("internal error in CHECK_GEN_LAYOUT: G-vectors changed " // str(NK) &
                 // " " // str(IND) // " " // str(WDES_NEW%IGX(IND,NK)) // " " // str(WDES%IGX(IND,&
                 NK)) // " " // str(WDES_NEW%IGY(IND,NK)) // " " // str(WDES%IGY(IND,NK)) // " " // &
                 str(WDES_NEW%IGZ(IND,NK)) // " " // str(WDES%IGZ(IND,NK)), __FILE__, __LINE__)
           ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE CHECK_GEN_LAYOUT


!********************** SUBROUTINE SET_INDPW_FULL  *********************
!
! SET_INDPW_FULL sets the index array for wavefunctions belonging to
! k-points (k') that have symmetry-equivalent k-points in the IRBZ
! for a Bloch vector k' not located in the IRBZ 
! IGX, IGY and IGZ is set to to  G' = S G, where S is the symmetry 
! operation taking k to k' = S k
!
!  the following arrays in the WDES array are initialised as well:
!   NGVECTOR         number of vectors
!   IGX, IGY, IGZ    G'
!   NINDPW           index corresponding to G
!  furthermore
!   KPOINTS_F%PHASE  phase factor acquired by the wavefunction when 
!                    going from k to k' 
!                    this is related to space group operations
!
!  see also ROTATE_WAVE_FFT
!
!***********************************************************************

    SUBROUTINE SET_INDPW_FULL(GRID, WDES, KPOINTS_F)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE mgrid
      USE wave
      USE mpimy
      USE constant
      USE sym_prec
      USE main_mpi
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      TYPE (wavedes) :: WDES
      TYPE (skpoints_full), POINTER :: KPOINTS_F
      
      INTEGER, ALLOCATABLE :: IGRIDIND(:,:,:)
      INTEGER, ALLOCATABLE :: IGRIDIND_MERGED(:,:,:)
      ! IGX, Y and Z merged over all cores
      INTEGER, POINTER     :: IGX_MERGED(:)
      INTEGER, POINTER     :: IGY_MERGED(:)
      INTEGER, POINTER     :: IGZ_MERGED(:)
      INTEGER              :: NC,NI,NI_LOCAL_ORIG,N1,N2,N3,NK,NKI
      INTEGER              :: NG1I,NG2I,NG3I,NG1,NG2,NG3
      INTEGER              :: NGX,NGY,NGZ,NGVECI,NGVEC
      LOGICAL              :: LSHIFT
      REAL(q)              :: PHASE

#ifdef _OPENACC
      IF (ASSOCIATED(KPOINTS_F).AND.ACC_IS_PRESENT(KPOINTS_F,1)) THEN
         PUSH_ACC_EXEC_ON(.TRUE.)
         CALL ACC_DELETE_TYPED_VAR(KPOINTS_F)
         POP_ACC_EXEC_ON
      ENDIF
#endif

      IF (WDES%NKPTS == KPOINTS_F%NKPTS) THEN
         ALLOCATE(KPOINTS_F%PHASE(WDES%NGDIM,KPOINTS_F%NKPTS))
#ifdef _OPENACC
         CALL ACC_COPYIN_TYPED_VAR(KPOINTS_F)
#endif
         RETURN
      ENDIF
      IF (WDES%NKDIM<KPOINTS_F%NKPTS) THEN
         CALL vtutor%bug("internal error in SET_INDPW_FOCK: WDES%NKDIM is too small for HF " // &
            str(WDES%NKDIM) // " " // str(KPOINTS_F%NKPTS), __FILE__, __LINE__)
      ENDIF
      IF (SIZE(WDES%NINDPW,2)<KPOINTS_F%NKPTS .OR. SIZE(WDES%IGX,2)<KPOINTS_F%NKPTS ) THEN
         CALL vtutor%bug("internal error in SET_INDPW_FOCK: WDES arrays not properly allocated", __FILE__, __LINE__)
      ENDIF

! allocate arrays and set some values to zero
      NGX=GRID%NGX
      NGY=GRID%NGY
      NGZ=GRID%NGZ
      ALLOCATE(IGRIDIND(NGX,NGY,NGZ), IGRIDIND_MERGED(NGX,NGY,NGZ) )
      ALLOCATE(KPOINTS_F%PHASE(WDES%NGDIM,KPOINTS_F%NKPTS))
      ALLOCATE(KPOINTS_F%NG_INDEX(WDES%NGDIM,KPOINTS_F%NKPTS))

      WDES%NINDPW(:,WDES%NKPTS+1:KPOINTS_F%NKPTS)=0
      KPOINTS_F%PHASE=(1._q,0._q)
! first we set up an array that allows to determine the local storage index
! for a 3d index  (x, y, z)
      NI=0
      IGRIDIND=0
      col1: DO NC=1,GRID%RC%NCOL
         N2=GRID%RC%I2(NC) ; NG2=GRID%LPCTY(N2)
         N3=GRID%RC%I3(NC) ; NG3=GRID%LPCTZ(N3)
         row1: DO N1=1,GRID%RC%NROW
            NI=NI+1
            NG1=GRID%LPCTX(N1)
            IGRIDIND(MP(NG1,NGX),MP(NG2,NGY),MP(NG3,NGZ))=NI
         ENDDO row1
      ENDDO col1

! IGRIDIND_MERGED is only used to check that the rotated data goes somewhere
! on some node
      IGRIDIND_MERGED=IGRIDIND
      CALLMPI( M_sum_i(WDES%COMM_INB, IGRIDIND_MERGED , SIZE(IGRIDIND_MERGED)))
! now index the wavefunctions at the new k-point NK
! take index from old k-point NKI, rotate G-vector, enter new index
      DO NK=WDES%NKPTS+1,KPOINTS_F%NKPTS
         
         ! phase shift (required for space group operations)
         IF ((ABS(KPOINTS_F%TRANS(1,NK))>TINY) .OR. &
             (ABS(KPOINTS_F%TRANS(2,NK))>TINY) .OR. &
             (ABS(KPOINTS_F%TRANS(3,NK))>TINY)) THEN
            LSHIFT=.TRUE.
         ELSE
            LSHIFT=.FALSE.
         ENDIF
         
         NKI=KPOINTS_F%NEQUIV(NK)
         NGVECI=WDES%NGVECTOR(NKI)
         NGVEC=1

#ifdef MPI
         ! at this point we need to merge the "source" indices
         ! from all nodes
         IF (WDES%COMM_INB%NCPU > 1) THEN
            NGVECI=WDES%NPLWKP_TOT(NKI)/WDES%NRSPINORS
            ALLOCATE(IGX_MERGED(NGVECI),IGY_MERGED(NGVECI),IGZ_MERGED(NGVECI))
            ! copy the local indices into the work array
            IGX_MERGED=0
            IGY_MERGED=0
            IGZ_MERGED=0
            IGX_MERGED(WDES%NGVECTOR_POS(NKI):WDES%NGVECTOR_POS(NKI)+WDES%NGVECTOR(NKI)-1)=WDES%IGX(1:WDES%NGVECTOR(NKI),NKI)
            IGY_MERGED(WDES%NGVECTOR_POS(NKI):WDES%NGVECTOR_POS(NKI)+WDES%NGVECTOR(NKI)-1)=WDES%IGY(1:WDES%NGVECTOR(NKI),NKI)
            IGZ_MERGED(WDES%NGVECTOR_POS(NKI):WDES%NGVECTOR_POS(NKI)+WDES%NGVECTOR(NKI)-1)=WDES%IGZ(1:WDES%NGVECTOR(NKI),NKI)
            CALL M_sum_i(WDES%COMM_INB,IGX_MERGED ,SIZE(IGX_MERGED))
            CALL M_sum_i(WDES%COMM_INB,IGY_MERGED ,SIZE(IGY_MERGED))
            CALL M_sum_i(WDES%COMM_INB,IGZ_MERGED ,SIZE(IGZ_MERGED))
         ELSE
            IGX_MERGED=>WDES%IGX(:,NKI)
            IGY_MERGED=>WDES%IGY(:,NKI)
            IGZ_MERGED=>WDES%IGZ(:,NKI)
         ENDIF
#else
         IGX_MERGED=>WDES%IGX(:,NKI)
         IGY_MERGED=>WDES%IGY(:,NKI)
         IGZ_MERGED=>WDES%IGZ(:,NKI)
#endif
         ! loop over all reciprocal lattice vectors (on all nodes)
         DO NI=1,NGVECI
            NG1I=IGX_MERGED(NI)
            NG2I=IGY_MERGED(NI)
            NG3I=IGZ_MERGED(NI)
            ! apply symmetry operation taking point in IRZ into full BZ
            NG1=NG1I*KPOINTS_F%IROTOP(1,1,NK)+NG2I*KPOINTS_F%IROTOP(2,1,NK)+ &
              NG3I*KPOINTS_F%IROTOP(3,1,NK)   
            NG2=NG1I*KPOINTS_F%IROTOP(1,2,NK)+NG2I*KPOINTS_F%IROTOP(2,2,NK)+ &
              NG3I*KPOINTS_F%IROTOP(3,2,NK)   
            NG3=NG1I*KPOINTS_F%IROTOP(1,3,NK)+NG2I*KPOINTS_F%IROTOP(2,3,NK)+ &
              NG3I*KPOINTS_F%IROTOP(3,3,NK)

            ! determine whether rotated PW goes somewhere on some node
            ! if this is not the case, we are in trouble
            IF (IGRIDIND_MERGED(MP(NG1,NGX),MP(NG2,NGY),MP(NG3,NGZ))==0) THEN
               CALL vtutor%bug("internal error in SET_INDPW_FULL: G vector not found " // str(NK) &
                  // " " // str(NI) // " " // str(NG1) // " " // str(NG2) // " " // str(NG3) // " " &
                  // str(NG1I) // " " // str(NG2I) // " " // str(NG3I) // " " // str(MP(NG1,NGX)) // &
                  " " // str(MP(NG2,NGY)) // " " // str(MP(NG3,NGZ)), __FILE__, __LINE__)
            ENDIF

            ! does data go somewhere on local node
            ! if so store the positions
            IF (IGRIDIND(MP(NG1,NGX),MP(NG2,NGY),MP(NG3,NGZ))/=0) THEN
               IF ( NGVEC>SIZE(WDES%NINDPW,1)) THEN
                  CALL vtutor%bug("internal error in SET_INDPW_FULL: insufficient memory (see wave.F " &
                     // "safeguard) " // str(NGVEC) // " " // str(SIZE(WDES%NINDPW,1)), __FILE__, __LINE__)
               ENDIF

               WDES%IGX(NGVEC,NK)=NG1
               WDES%IGY(NGVEC,NK)=NG2
               WDES%IGZ(NGVEC,NK)=NG3

               WDES%NINDPW(NGVEC,NK)= &
                    IGRIDIND(MP(NG1,NGX),MP(NG2,NGY),MP(NG3,NGZ))
               ! store index, if NCPU=1 then NI = NGVEC
               KPOINTS_F%NG_INDEX(NGVEC,NK)=NI

               NGVEC=NGVEC+1
            ENDIF
            ! if the data originates from the local node  and shift is required
            ! store shift
            IF (LSHIFT) THEN
               NI_LOCAL_ORIG=NI-WDES%NGVECTOR_POS(NKI)+1
               IF (NI_LOCAL_ORIG >=1 .AND. NI_LOCAL_ORIG <= WDES%NGVECTOR(NKI)) THEN
                  PHASE=-TPI*(NG1*KPOINTS_F%TRANS(1,NK)+ &
                              NG2*KPOINTS_F%TRANS(2,NK)+ &
                              NG3*KPOINTS_F%TRANS(3,NK))
                  KPOINTS_F%PHASE(NI_LOCAL_ORIG,NK)=CMPLX(COS(PHASE),SIN(PHASE),KIND=q)
               ENDIF
            ENDIF
         ENDDO
         ! store number of local g-vectors for k in full BZ
         WDES%NGVECTOR(NK)=NGVEC-1
#ifdef MPI
         IF (WDES%COMM_INB%NCPU > 1) THEN
            DEALLOCATE(IGX_MERGED, IGY_MERGED, IGZ_MERGED)
         ENDIF

         NGVEC=WDES%NGVECTOR(NK)
         CALL M_sum_i(WDES%COMM_INB,NGVEC,1)
         IF (NGVEC/=WDES%NPLWKP_TOT(NKI)/WDES%NRSPINORS) THEN
            CALL vtutor%bug("internal error in SET_INDPW_FULL: number of final plane waves wrong "&
                // str(NK) // " " // str(NKI) // " " // str(NGVEC) // " " // str(WDES%NPLWKP_TOT(NKI)), &
                __FILE__, __LINE__)
         ENDIF
#endif
      ENDDO
      DEALLOCATE(IGRIDIND, IGRIDIND_MERGED)

#ifdef _OPENACC
      CALL ACC_COPYIN_TYPED_VAR(KPOINTS_F)
#endif
    END SUBROUTINE SET_INDPW_FULL


!**************************** FUNCTION MP ******************************
!
! small helper function, produces positive indices
!
!*********************************************************************** 

  FUNCTION MP(IND,MAXI)
    IMPLICIT NONE
    INTEGER MP,MAXI,IND
    IF (IND<=0) THEN
       MP=MAXI+IND
    ELSE
       MP=IND
    ENDIF
  END FUNCTION MP

!***********************************************************************
!
! search a particular k-point in the full list and return 
! the equivalent k-point in the IRBZ
!
!*********************************************************************** 


  FUNCTION KPOINT_IN_FULL_GRID(VKPT,KPOINTS_F)
    USE sym_prec
    USE string, ONLY: str
    USE tutor, ONLY: vtutor
    INTEGER :: KPOINT_IN_FULL_GRID
    REAL(q) VKPT(:)
    TYPE (skpoints_full) :: KPOINTS_F
 ! local
    REAL (q) :: X, Y, Z
    INTEGER :: IX, IY, IZ
    INTEGER NK, NK_TEST

    NK_TEST=-1
    IF ( ASSOCIATED(KPOINTS_F%MAP_INTO_BZ)) THEN
       X=MOD((VKPT(1)-KPOINTS_F%VKPT(1,1))*KPOINTS_F%NKPX+KPOINTS_F%NKPX*10, KPOINTS_F%NKPX*1.0_q)
       Y=MOD((VKPT(2)-KPOINTS_F%VKPT(2,1))*KPOINTS_F%NKPY+KPOINTS_F%NKPY*10, KPOINTS_F%NKPY*1.0_q)
       Z=MOD((VKPT(3)-KPOINTS_F%VKPT(3,1))*KPOINTS_F%NKPZ+KPOINTS_F%NKPZ*10, KPOINTS_F%NKPZ*1.0_q)

       IF (ABS(X-NINT(X)) < TINY .AND. ABS(Y-NINT(Y)) < TINY .AND. ABS(Z-NINT(Z)) < TINY) THEN
          IX=MOD( NINT(X), KPOINTS_F%NKPX)+1
          IY=MOD( NINT(Y), KPOINTS_F%NKPY)+1
          IZ=MOD( NINT(Z), KPOINTS_F%NKPZ)+1
          NK_TEST=KPOINTS_F%MAP_INTO_BZ(IX, IY, IZ)
          KPOINT_IN_FULL_GRID=NK_TEST
          RETURN
       ENDIF
    ENDIF

    DO NK=1,KPOINTS_F%NKPTS
       IF ( & 
         (ABS(MOD(VKPT(1)-KPOINTS_F%VKPT(1,NK)+10.5_q,1._q)-0.5_q)<TINY) .AND. &
         (ABS(MOD(VKPT(2)-KPOINTS_F%VKPT(2,NK)+10.5_q,1._q)-0.5_q)<TINY) .AND. &
         (ABS(MOD(VKPT(3)-KPOINTS_F%VKPT(3,NK)+10.5_q,1._q)-0.5_q)<TINY)) EXIT
    ENDDO

    IF (NK>KPOINTS_F%NKPTS) THEN
       CALL vtutor%bug("internal error in KPOINT_IN_FULL_GRID: can not find " // str(VKPT), __FILE__, __LINE__)
    ENDIF

    KPOINT_IN_FULL_GRID=NK

    IF (KPOINT_IN_FULL_GRID/= NK_TEST .AND. NK_TEST /= -1) THEN
       CALL vtutor%bug("internal error in VASP: KPOINT_IN_FULL_GRID is wrong " // str(NK_TEST) // &
          " " // str(KPOINT_IN_FULL_GRID) // " " // str(VKPT(1:3) - KPOINTS_F%VKPT(1:3,NK)), __FILE__, __LINE__)
    ENDIF
  END FUNCTION KPOINT_IN_FULL_GRID

!
! create a mapping for fast indexing of a k-point into KPOINTS_FULL
!
! mM: 02062020: this subroutine is broken for generalized regular grids ...

  SUBROUTINE SET_KPOINT_MAP_INTO_BZ( KPOINTS_F)
    USE sym_prec
    USE string, ONLY: str
    USE tutor, ONLY: vtutor
    TYPE (skpoints_full) :: KPOINTS_F
    REAL (q) :: X, Y, Z
    INTEGER :: IX, IY, IZ, NK

    IF (KPOINTS_F%NKPX*KPOINTS_F%NKPY*KPOINTS_F%NKPZ == KPOINTS_F%NKPTS_NON_ZERO ) THEN
       ALLOCATE(KPOINTS_F%MAP_INTO_BZ(KPOINTS_F%NKPX, KPOINTS_F%NKPY, KPOINTS_F%NKPZ))
       KPOINTS_F%MAP_INTO_BZ=-1

       DO NK=1,KPOINTS_F%NKPTS
          X=MOD((KPOINTS_F%VKPT(1,NK)-KPOINTS_F%VKPT(1,1))*KPOINTS_F%NKPX+KPOINTS_F%NKPX*10,KPOINTS_F%NKPX*1.0_q)
          Y=MOD((KPOINTS_F%VKPT(2,NK)-KPOINTS_F%VKPT(2,1))*KPOINTS_F%NKPY+KPOINTS_F%NKPY*10,KPOINTS_F%NKPY*1.0_q)
          Z=MOD((KPOINTS_F%VKPT(3,NK)-KPOINTS_F%VKPT(3,1))*KPOINTS_F%NKPZ+KPOINTS_F%NKPZ*10,KPOINTS_F%NKPZ*1.0_q)

          IF (ABS(X-NINT(X)) < TINY .AND. ABS(Y-NINT(Y)) < TINY .AND. ABS(Z-NINT(Z)) < TINY) THEN
             IX=MOD( NINT(X), KPOINTS_F%NKPX)+1
             IY=MOD( NINT(Y), KPOINTS_F%NKPY)+1
             IZ=MOD( NINT(Z), KPOINTS_F%NKPZ)+1
             IF (KPOINTS_F%MAP_INTO_BZ(IX, IY, IZ)/=-1) THEN
                CALL vtutor%bug("internal error in VASP: SET_KPOINT_MAP_INTO_BZ is already set " // &
                   str(KPOINTS_F%VKPT(1,NK) - KPOINTS_F%VKPT(1,1)) // " " // str(IX) // " " // str(IY) // " " // str(IZ), &
                   __FILE__, __LINE__)
             ENDIF

             KPOINTS_F%MAP_INTO_BZ(IX, IY, IZ) = NK
          ENDIF
       ENDDO

       DO IX=1,KPOINTS_F%NKPX
          DO IY=1,KPOINTS_F%NKPY
             DO IZ=1,KPOINTS_F%NKPZ
                IF (KPOINTS_F%MAP_INTO_BZ(IX, IY, IZ)==-1) THEN
                   CALL vtutor%bug("internal error in VASP: SET_KPOINT_MAP_INTO_BZ not all k-points " &
                      // "found " // str(IX) // " " // str(IY) // " " // str(IZ), __FILE__, __LINE__)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF

  END SUBROUTINE SET_KPOINT_MAP_INTO_BZ

!*********************************************************************** 
!
! small helper routine to identify identical k-points
! applies a minimum image convention to avoid problems at the
! boundary of the Brillouine zone
! returns .TRUE. if the k-points are identical
!
!*********************************************************************** 

  FUNCTION LIDENTICAL_KPOINT(VKPT1,VKPT2)
    USE sym_prec
    LOGICAL LIDENTICAL_KPOINT
    REAL(q) VKPT1(3),VKPT2(3)

    LIDENTICAL_KPOINT= &
         (ABS(MOD(VKPT1(1)-VKPT2(1)+10.5_q,1._q)-0.5_q)<TINY) .AND. &
         (ABS(MOD(VKPT1(2)-VKPT2(2)+10.5_q,1._q)-0.5_q)<TINY) .AND. &
         (ABS(MOD(VKPT1(3)-VKPT2(3)+10.5_q,1._q)-0.5_q)<TINY)

  END FUNCTION LIDENTICAL_KPOINT


!*********************************************************************
!
! determine the symmetry operations that leave the k-point NQ
! invariant
! then loop over second set NK2 and determine whether symmetry
! yields redundant pairs
!
!*********************************************************************

  SUBROUTINE SMALL_SPACE_GROUP_WEIGHTS( KPOINTS, KPOINTS_F, ROTMAP, MAGROT,  ISYM, LINVERSION, IU6)
    USE wave
    USE constant
    USE mkpoints
    USE string, ONLY: str
    USE tutor, ONLY: vtutor
! passed structures and variables
    TYPE (kpoints_struct)::  KPOINTS
    TYPE (skpoints_full) :: KPOINTS_F
    INTEGER ISYM
    INTEGER ROTMAP(:,:,:)
    REAL(q) MAGROT(:,:)
    LOGICAL LINVERSION
    INTEGER :: IU6
! external routine
    LOGICAL,EXTERNAL :: SGREQL
! common symmetry variables
    INTEGER ISYMOP, NROT, IGRPOP, NROTK, INVMAP, NPCELL
    REAL(q) GTRANS, AP
    COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
         &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL
! local variables
    REAL(q) :: TINY=1E-6, V(3), VR(3), ROP(3,3)
    INTEGER :: INVERS(9)=(/-1,0,0,0,-1,0,0,0,-1/),IOP(3,3)
    INTEGER NOP, N, NK1, NK2, NKF
    LOGICAL LINV, LEXIST

    INTEGER SMALL_GROUP_OP(48)
    REAL(q) SMALL_GROUP_ROP(3,3,48)
    INTEGER ITRANS(3,48)
    LOGICAL SMALL_GROUP_LINV(48)
    INTEGER NOP_SMALL_GROUP

    LOGICAL K_POINT_DONE(KPOINTS_F%NKPTS)

    IF (ALLOCATED(WEIGHT_K_POINT_PAIR_SMALL_GROUP)) THEN
       DEALLOCATE(WEIGHT_K_POINT_PAIR_SMALL_GROUP)
    ENDIF

    IF (.NOT. LSYMGRAD_SAVE)  RETURN

    ALLOCATE(WEIGHT_K_POINT_PAIR_SMALL_GROUP(KPOINTS%NKPTS, KPOINTS_F%NKPTS))

    IF (IU6>=0) THEN
       WRITE(IU6,'(A)') ' Determining symmetry at each k-point (small space group operations):'
    ENDIF

    WEIGHT_K_POINT_PAIR_SMALL_GROUP=0
!=======================================================================
! loop over all point group operations and check whether 
! symmetry leaves the k-point invariant
!=======================================================================
    DO NK1=1,KPOINTS%NKPTS
       NOP_SMALL_GROUP=0
       VR(1)=KPOINTS%VKPT(1,NK1)
       VR(2)=KPOINTS%VKPT(2,NK1)
       VR(3)=KPOINTS%VKPT(3,NK1)
       DO NOP=1,NROTK
          ! test existence of inversion op
          IF (SGREQL(IGRPOP(1,1,NOP),INVERS)) LINV=.TRUE.
          ! copy symmetry op to real array
          ROP=IGRPOP(:,:,NOP)
          
          ! generate new k-point S k
          V(1)=VR(1)*ROP(1,1)+VR(2)*ROP(2,1)+VR(3)*ROP(3,1)
          V(2)=VR(1)*ROP(1,2)+VR(2)*ROP(2,2)+VR(3)*ROP(3,2)
          V(3)=VR(1)*ROP(1,3)+VR(2)*ROP(2,3)+VR(3)*ROP(3,3)
          ! compare with original k-point
          IF ( ( ABS(MOD(VR(1)-V(1)+6.5,1._q)-0.5_q)<TINY) .AND. &
               ( ABS(MOD(VR(2)-V(2)+6.5,1._q)-0.5_q)<TINY) .AND. &
               ( ABS(MOD(VR(3)-V(3)+6.5,1._q)-0.5_q)<TINY)) THEN
             NOP_SMALL_GROUP=NOP_SMALL_GROUP+1
             ! round to next integer
             ITRANS(:,NOP_SMALL_GROUP)=CEILING(VR(:)-V(:)-0.5)
             SMALL_GROUP_OP (NOP_SMALL_GROUP)=NOP
             SMALL_GROUP_LINV(NOP_SMALL_GROUP)=.FALSE.
             SMALL_GROUP_ROP(:,:,NOP_SMALL_GROUP)=ROP
          ENDIF
       END DO
    ! same now applying time inversion symmetry
       IF (.NOT. LINV .AND. ISYM>=0 .AND. LINVERSION) THEN
          DO NOP=1,NROTK
             ! apply inversion symmetry to form to get IOP
             CALL SGRPRD(INVERS,IGRPOP(1,1,NOP),IOP(1,1))
             ! copy symmetry op to real array
             ROP=IOP
             V(1)=VR(1)*ROP(1,1)+VR(2)*ROP(2,1)+VR(3)*ROP(3,1)
             V(2)=VR(1)*ROP(1,2)+VR(2)*ROP(2,2)+VR(3)*ROP(3,2)
             V(3)=VR(1)*ROP(1,3)+VR(2)*ROP(2,3)+VR(3)*ROP(3,3)
             ! compare with original k-point
             IF ( ( ABS(MOD(VR(1)-V(1)+6.5,1._q)-0.5_q)<TINY) .AND. &
                  ( ABS(MOD(VR(2)-V(2)+6.5,1._q)-0.5_q)<TINY) .AND. &
                  ( ABS(MOD(VR(3)-V(3)+6.5,1._q)-0.5_q)<TINY)) THEN
                NOP_SMALL_GROUP=NOP_SMALL_GROUP+1
                ITRANS(:,NOP_SMALL_GROUP)=CEILING(VR(:)-V(:)-0.5)
                SMALL_GROUP_OP (NOP_SMALL_GROUP)=NOP
                SMALL_GROUP_LINV(NOP_SMALL_GROUP)=.TRUE.
                SMALL_GROUP_ROP(:,:,NOP_SMALL_GROUP)=ROP
             ENDIF
          END DO
       ENDIF
       IF (IU6>=0) WRITE(IU6,*) 'NK1=',NK1,' operations',NOP_SMALL_GROUP
!=======================================================================
! 
!=======================================================================
       K_POINT_DONE=.FALSE.
       DO NK2=1,KPOINTS_F%NKPTS
          ! zero weight k-point cycle
          IF (KPOINTS_F%WTKPT(NK2)==0) CYCLE
          ! already done this k-point, cycle
          IF (K_POINT_DONE(NK2)) CYCLE

          ! get coordinates of k-point
          VR(1)=KPOINTS_F%VKPT(1,NK2)
          VR(2)=KPOINTS_F%VKPT(2,NK2)
          VR(3)=KPOINTS_F%VKPT(3,NK2)

          ! loop over all small space group symmetry operations
          DO N=1,NOP_SMALL_GROUP
             NOP =SMALL_GROUP_OP (N)
             LINV=SMALL_GROUP_LINV(N)
             ROP =SMALL_GROUP_ROP(:,:,N)
             V(1)=VR(1)*ROP(1,1)+VR(2)*ROP(2,1)+VR(3)*ROP(3,1)
             V(2)=VR(1)*ROP(1,2)+VR(2)*ROP(2,2)+VR(3)*ROP(3,2)
             V(3)=VR(1)*ROP(1,3)+VR(2)*ROP(2,3)+VR(3)*ROP(3,3)
             LEXIST=.FALSE.
             ! match onto k-point in full grid
             DO NKF=1,KPOINTS_F%NKPTS
               ! zero weight k-points are again excluded
                IF(( KPOINTS_F%WTKPT(NKF)/=0) .AND. &
                  ( ABS(MOD(V(1)-KPOINTS_F%VKPT(1,NKF)+6.5,1._q)-0.5_q)<TINY) .AND. &
                  ( ABS(MOD(V(2)-KPOINTS_F%VKPT(2,NKF)+6.5,1._q)-0.5_q)<TINY) .AND. &
                  ( ABS(MOD(V(3)-KPOINTS_F%VKPT(3,NKF)+6.5,1._q)-0.5_q)<TINY)) THEN
                  ! if we hit a new k-point mark that we have hit it
                  ! and increase weight
                  IF (.NOT. K_POINT_DONE(NKF)) THEN
                     WEIGHT_K_POINT_PAIR_SMALL_GROUP(NK1,NK2)=WEIGHT_K_POINT_PAIR_SMALL_GROUP(NK1,NK2)+1
                     ! found this mark it so no need to do this one again
                     K_POINT_DONE(NKF)=.TRUE.
                  ENDIF
                  ! when we have found the corresponding k-point in the full grid
                  ! we are done
                  LEXIST=.TRUE.
                  EXIT
               ENDIF
            ENDDO
            IF (.NOT. LEXIST) THEN
               CALL vtutor%bug("internal error in SMALL_SPACE_GROUP_WEIGHTS: NK2 not found", __FILE__, __LINE__)
            ENDIF
         ENDDO
      ENDDO
!      WRITE(*,'(16F4.1)') WEIGHT_K_POINT_PAIR_SMALL_GROUP(NK1,:)
      IF (SUM(WEIGHT_K_POINT_PAIR_SMALL_GROUP(NK1,:)) /=KPOINTS_F%NKPTS_NON_ZERO) THEN
         CALL vtutor%bug(str(WEIGHT_K_POINT_PAIR_SMALL_GROUP(NK1,:)) // "\n internal error in " &
            // "SMALL_SPACE_GROUP_WEIGHTS: NK2 not found", __FILE__, __LINE__)
      ENDIF
         
   ENDDO

 END SUBROUTINE SMALL_SPACE_GROUP_WEIGHTS


!********************** SUBROUTINE ROTATE_WAVE     *********************
!
! determine the wavefunction (in real space and  the wave function character) 
! at the k-point WQ%WDES1%NK from the corresponding k-point in the IRBZ
! the wavefunction at the k-point in the IRBZ must be supplied in W_ORIG
! see also ROTATE_WAVE_CHARACTER and ROTATE_WAVE_FFT
!
!***********************************************************************
#ifndef gammareal
  SUBROUTINE W1_ROTATE_AND_FFT(W_NEW, W_ORIG, ROT_HANDLE, P, LATT_CUR, LSHIFT)
    USE wave_high
    USE pseudo
    USE lattice
#ifndef oldsym
    USE spinsym
#endif
    TYPE( rotation_handle), POINTER :: ROT_HANDLE
    TYPE (potcar)      P(:)
    TYPE (latt)        LATT_CUR
    TYPE (wavefun1)    W_NEW         ! new wavefunction 
    TYPE (wavefun1)    W_ORIG
    LOGICAL LSHIFT

    INTEGER NQ, ISPINOR

    NQ=W_NEW%WDES1%NK

    ! rotate wavefunction
    DO ISPINOR=0,W_NEW%WDES1%NRSPINORS-1
       CALL ROTATE_WAVE_FFT_MRG(W_ORIG%WDES1%COMM_INB, W_ORIG%WDES1%NGVECTOR, W_ORIG%WDES1%NGVECTOR_POS, W_ORIG%WDES1%NPL_TOT/W_ORIG%WDES1%NRSPINORS, KPOINTS_FULL%NG_INDEX(1,NQ),  & 
            W_NEW%WDES1%NGVECTOR, W_NEW%WDES1%NINDPW(1), &
            KPOINTS_FULL%PHASE(1,NQ),W_NEW%CR(1+ISPINOR*W_NEW%WDES1%GRID%MPLWV), &
            W_ORIG%CW(1+ISPINOR*W_ORIG%WDES1%NGVECTOR), &
            W_NEW%CW(1+ISPINOR*W_NEW%WDES1%NGVECTOR),W_NEW%WDES1%GRID,KPOINTS_FULL%LINV(NQ),LSHIFT)
    ENDDO
    ! rotate the wavefunction character
    IF (W_NEW%WDES1%LOVERL) &
         CALL ROTATE_WAVE_CHARACTER(ROT_HANDLE, P, LATT_CUR, W_NEW%WDES1, W_ORIG, W_NEW)

    ! test recalculate wave function character
    ! IF (W_NEW%WDES1%LOVERL) CALL W1_PROJ(W_NEW, NONLR_S, NONL_S)

#ifndef oldsym
    CALL ROTATE_WAVE_SPIN(W_NEW, KPOINTS_FULL%RSSYMOP(:,:,NQ))
#endif

  END SUBROUTINE W1_ROTATE_AND_FFT
#endif

#ifndef gammareal
  SUBROUTINE W1_ROTATE_AND_FFT_NO_PROJ(W_NEW, W_ORIG, LSHIFT)
    USE wave_high
    USE pseudo
    USE lattice
#ifndef oldsym
    USE spinsym
#endif
    TYPE (wavefun1)    W_NEW         ! new wavefunction 
    TYPE (wavefun1)    W_ORIG
    LOGICAL LSHIFT

    INTEGER NQ, ISPINOR

    NQ=W_NEW%WDES1%NK
    
    ! rotate wavefunction
    DO ISPINOR=0,W_NEW%WDES1%NRSPINORS-1
       CALL ROTATE_WAVE_FFT_MRG(W_ORIG%WDES1%COMM_INB, W_ORIG%WDES1%NGVECTOR, W_ORIG%WDES1%NGVECTOR_POS, W_ORIG%WDES1%NPL_TOT/W_ORIG%WDES1%NRSPINORS, KPOINTS_FULL%NG_INDEX(1,NQ), & 
            W_NEW%WDES1%NGVECTOR,W_NEW%WDES1%NINDPW(1), &
            KPOINTS_FULL%PHASE(1,NQ),W_NEW%CR(1+ISPINOR*W_NEW%WDES1%GRID%MPLWV), &
            W_ORIG%CW(1+ISPINOR*W_ORIG%WDES1%NGVECTOR), &
            W_NEW%CW(1+ISPINOR*W_NEW%WDES1%NGVECTOR),W_NEW%WDES1%GRID,KPOINTS_FULL%LINV(NQ),LSHIFT)
    ENDDO

#ifndef oldsym
    CALL ROTATE_WAVE_SPIN(W_NEW, KPOINTS_FULL%RSSYMOP(:,:,NQ))
#endif

  END SUBROUTINE W1_ROTATE_AND_FFT_NO_PROJ
#endif


!*********************************************************************
!
! subroutine to rotate the wavefunction character (W%CPROJ)
! from a kpoint in the IRBZ into a desired k-points WDES1%NK, which
! is not located in the IRBZ 
!
! the routine uses a rotation_handle (ROT_HANDLE) to cache intermediate 
! results and work arrays
! the handle is created on the fly upon the first call, and
! reinitialised whenever the WDES1%NK changes
!
!*********************************************************************

  SUBROUTINE ROTATE_WAVE_CHARACTER(ROT_HANDLE, P, LATT_CUR, WDES1, W, W_NEW)
#ifdef _OPENACC
!$ACC ROUTINE(ROTATE_VECTOR) VECTOR
    USE mopenacc
#endif
    USE wave_high
    USE pseudo
    USE lattice
    USE nonl_high
#ifndef oldsym
    USE spinsym
#endif
    TYPE( rotation_handle), POINTER :: ROT_HANDLE
    TYPE (potcar)      P(:)
    TYPE (wavedes1)    WDES1
    TYPE (latt)        LATT_CUR
    TYPE (wavefun1)    W, W_NEW
  ! local
    INTEGER :: NPRO, NT, NI, NI_GLBL, NPROP, NIP, ISPINOR
  ! local for setup
    GDEF :: CPROJ(WDES1%NPRO_TOT)

    CALL SET_ROT_HANDLE(P,LATT_CUR,WDES1,ROT_HANDLE)

!$ACC ENTER DATA CREATE(CPROJ) __IF_ASYNC__
    CALL MRG_PROJ(WDES1,CPROJ(1),W%CPROJ)

!$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(WDES1,KPOINTS_FULL,ROT_HANDLE,CPROJ,W_NEW) &
!$ACC PRIVATE(NT,NPRO,NI_GLBL,NIP,NPROP) __IF_ASYNC__
    DO ISPINOR=0,WDES1%NRSPINORS-1
       DO NI=1,WDES1%NIONS
          NT=WDES1%ITYP(NI)
          NPRO=WDES1%LMBASE(NI)+ISPINOR*(WDES1%NPRO/2)+1
#ifdef MPI
          ! global ionic index (see NI_GLBL function in wave.F)
          NI_GLBL=(NI-1)*WDES1%COMM_INB%NCPU+WDES1%COMM_INB%NODE_ME
#else
          NI_GLBL=NI
#endif
          NIP=KPOINTS_FULL%ROTMAP(NI_GLBL,WDES1%NK)
          NPROP=ROT_HANDLE%NPRO_NI(NIP)+ISPINOR*(WDES1%NPRO_TOT/2)

          CALL ROTATE_VECTOR( KPOINTS_FULL%LINV(WDES1%NK),CPROJ(NPROP),W_NEW%CPROJ(NPRO), &
               ROT_HANDLE%MMAX,ROT_HANDLE%LDIM,ROT_HANDLE%SL(1,1,0),ROT_HANDLE%NPS(NT),ROT_HANDLE%LPS(1,NT))
       ENDDO
    ENDDO

!$ACC EXIT DATA DELETE(CPROJ) __IF_ASYNC__

#ifndef oldsym
    CALL ROTATE_WAVE_CHARACTER_SPIN(WDES1, W_NEW%CPROJ, KPOINTS_FULL%RSSYMOP(:,:,WDES1%NK) )
#endif

  END SUBROUTINE ROTATE_WAVE_CHARACTER


!*********************************************************************
!
! subroutine to allocate a rotation handle, and initialize it
! for a specific k-point
!
!*********************************************************************

  SUBROUTINE SET_ROT_HANDLE(P,LATT_CUR,WDES1,ROT_HANDLE)
    USE wave_high
    USE pseudo
    USE lattice
    USE nonl_high
    USE asa

    TYPE( rotation_handle), POINTER :: ROT_HANDLE
    TYPE (potcar)      P(:)
    TYPE (wavedes1)    WDES1
    TYPE (latt)        LATT_CUR
 ! local
    INTEGER :: NPRO, NIS, NT, NI, NT_GLOBAL
    INTEGER, EXTERNAL :: MAXL

    INTEGER :: NITYP(SIZE(P)), NTYP, NIONS

    IF (.NOT. ASSOCIATED(ROT_HANDLE)) THEN
       ALLOCATE (ROT_HANDLE)

       ! old version: NTYP=WDES1%NTYP ; NITYP=WDES1%NITYP
       NTYP=SIZE(P)
       ! determine number of ions for each type
       NITYP=0
       DO NT=1, NTYP
          IF (WDES1%NITYP(NT)>0) THEN
             NITYP(WDES1%NT_GLOBAL(NT))=WDES1%NITYP(NT)
          ENDIF
       ENDDO
       CALLMPI( M_sum_i(WDES1%COMM_INB, NITYP , SIZE(NITYP)))
       NIONS=SUM(NITYP)
       
       ROT_HANDLE%LDIM=MAXL(SIZE(P),P(1))    ! maximum l quantum number
       ROT_HANDLE%MMAX=(2*ROT_HANDLE%LDIM+1) ! maximum m quantum number
       ALLOCATE (ROT_HANDLE%SL(ROT_HANDLE%MMAX,ROT_HANDLE%MMAX,0:ROT_HANDLE%LDIM), ROT_HANDLE%NPRO_NI(NIONS))

       ! setup index array NPRO_NI into globally summed CPROJ array
       ! this array is compatible to seriel version
       NPRO=1
       NIS=1
       DO NT=1, NTYP
          DO NI=NIS, NITYP(NT)+NIS-1
             ROT_HANDLE%NPRO_NI(NI)=NPRO
             NPRO= P(NT)%LMMAX+NPRO
          ENDDO
          NIS=NIS+NITYP(NT)
       ENDDO

       IF (NPRO-1 /= WDES1%NPRO_TOT/WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in ROTATE_WAVE_CHARACTER: NPRO is not correct " // &
             str(NPRO) // " " // str(WDES1%NPRO_TOT), __FILE__, __LINE__)
       ENDIF

       ! and copy some stuff from the potcar structure into ROT_HANDLE
       ALLOCATE(ROT_HANDLE%NPS(NTYP),ROT_HANDLE%LPS(P(1)%LDIM,NTYP))
       ROT_HANDLE%NPS=0; ROT_HANDLE%LPS=0
       DO NT=1,WDES1%NTYP
          NT_GLOBAL=WDES1%NT_GLOBAL(NT)
          ROT_HANDLE%NPS(NT)=P(NT_GLOBAL)%LMAX
          ROT_HANDLE%LPS(1:P(NT_GLOBAL)%LMAX,NT)=P(NT_GLOBAL)%LPS(1:P(NT_GLOBAL)%LMAX)
       ENDDO

       ! this signals the fact that ROT_HANDLE%SL has not been set yet
       ROT_HANDLE%NK=-1

!$ACC ENTER DATA COPYIN(ROT_HANDLE) __IF_ASYNC__
!$ACC ENTER DATA COPYIN(ROT_HANDLE%NPRO_NI,ROT_HANDLE%NPS,ROT_HANDLE%LPS) __IF_ASYNC__
!$ACC ENTER DATA CREATE(ROT_HANDLE%SL) __IF_ASYNC__
    ENDIF

    IF (ROT_HANDLE%NK/= WDES1%NK) THEN
       ROT_HANDLE%NK = WDES1%NK
       CALL SETUP_SYM_LL(ROT_HANDLE%MMAX,ROT_HANDLE%LDIM,KPOINTS_FULL%ISYMOP(:,:,WDES1%NK), & 
            ROT_HANDLE%SL,LATT_CUR%A,LATT_CUR%B)
!$ACC UPDATE DEVICE(ROT_HANDLE%NK,ROT_HANDLE%SL) IF_PRESENT __IF_ASYNC__
    ENDIF
  END SUBROUTINE SET_ROT_HANDLE

!*********************************************************************
!
! subroutine to deallocate a rotation handle
!
!*********************************************************************

  SUBROUTINE DEALLOCATE_ROT_HANDLE(ROT_HANDLE)
#ifdef _OPENACC
    USE mopenacc
#endif
    TYPE( rotation_handle), POINTER :: ROT_HANDLE

    IF (ASSOCIATED(ROT_HANDLE)) THEN
!$ACC EXIT DATA DELETE(ROT_HANDLE%SL,ROT_HANDLE%NPRO_NI,ROT_HANDLE%NPS,ROT_HANDLE%LPS,ROT_HANDLE) FINALIZE __IF_ASYNC__
       IF (ASSOCIATED(ROT_HANDLE%SL)) DEALLOCATE(ROT_HANDLE%SL)
       IF (ASSOCIATED(ROT_HANDLE%NPRO_NI)) DEALLOCATE(ROT_HANDLE%NPRO_NI)
       IF (ASSOCIATED(ROT_HANDLE%NPS)) DEALLOCATE(ROT_HANDLE%NPS)
       IF (ASSOCIATED(ROT_HANDLE%LPS)) DEALLOCATE(ROT_HANDLE%LPS)
       DEALLOCATE(ROT_HANDLE)
    ENDIF
    NULLIFY( ROT_HANDLE ) 
  
  END SUBROUTINE DEALLOCATE_ROT_HANDLE

END MODULE full_kpoints


!********************** SUBROUTINE COUNT_INDPW_FULL ********************
!
! This subroutine determines NRPLWV: the maximum number of plane waves
! per orbital each MPI-rank must be able to store for orbitals at
! k-points with NK = WDES%NKPTS+1, .. , KPOINTS_FULL%NKPTS, i.e., for
! orbitals at k-points of the full k-grid that are not part of the
! irreducible wedge of the 1BZ.
!
! N.B.: Currently COUNT_INDPW_FULL is only called from GEN_LAYOUT,
! and on entry NRPLWV is set to the maximum number of plane waves per
! orbital each MPI-rank must be able to store for orbitals at
! k-points that are part of the irreducible wedge of the 1BZ.
!
!***********************************************************************

    SUBROUTINE COUNT_INDPW_FULL(GRID, WDES, BI, NRPLWV)
      USE prec
      USE mgrid
      USE wave
      USE constant
      USE full_kpoints

      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      TYPE (wavedes) :: WDES

      REAL(q) :: BI(3,3)
      INTEGER :: NRPLWV

      INTEGER, ALLOCATABLE :: IGRIDIND(:,:,:)
      INTEGER, ALLOCATABLE :: IGRIDIND_MERGED(:,:,:)

      INTEGER, POINTER     :: IGX_MERGED(:)
      INTEGER, POINTER     :: IGY_MERGED(:)
      INTEGER, POINTER     :: IGZ_MERGED(:)

      INTEGER, ALLOCATABLE :: NGVECTOR_MERGED(:)

      INTEGER :: NI,NC,N1,N2,N3,NG1,NG2,NG3,NG1I,NG2I,NG3I
      INTEGER :: NGX,NGY,NGZ,NK,NKI,NGVEC,NGVECI
      REAL(q) :: G1,G2,G3,GIX,GIY,GIZ,ENERGI

#ifndef MPI
      ! in the serial version of the code this routine
      ! does not need to be called
      RETURN
#else
      ! exit if KPOINTS_FULL has not been set up
      IF (.NOT.LFULL_KPOINTS) RETURN
      ! exit if the number of k-points in the full grid equals
      ! the number of k-points in the irreducible grid
      IF (WDES%NKPTS == KPOINTS_FULL%NKPTS) RETURN
      ! exit if the plane wave coefficients of an orbital are
      ! all owned by a single MPI rank
      IF (WDES%COMM_INB%NCPU == 1) RETURN

      NGX=GRID%NGX
      NGY=GRID%NGY
      NGZ=GRID%NGZ

      ALLOCATE(IGRIDIND(NGX,NGY,NGZ),IGRIDIND_MERGED(NGX,NGY,NGZ))

! first we set up an array that allows to determine the local storage index for a 3d index (x, y, z)
      NI=0
      IGRIDIND=0
      col1: DO NC=1,GRID%RC%NCOL
         N2=GRID%RC%I2(NC) ; NG2=GRID%LPCTY(N2)
         N3=GRID%RC%I3(NC) ; NG3=GRID%LPCTZ(N3)
         row1: DO N1=1,GRID%RC%NROW
            NI=NI+1
            NG1=GRID%LPCTX(N1)
            IGRIDIND(MP(NG1,NGX),MP(NG2,NGY),MP(NG3,NGZ))=NI
         ENDDO row1
      ENDDO col1

! IGRIDIND_MERGED is only used to check that the rotated data goes somewher on some node
      IGRIDIND_MERGED=IGRIDIND
      CALLMPI( M_sum_i(WDES%COMM_INB, IGRIDIND_MERGED , SIZE(IGRIDIND_MERGED)))

! now index the wavefunctions at the new k-point NK
! take index from old k-point NKI, rotate G-vector, enter new index
      DO NK=WDES%NKPTS+1,KPOINTS_FULL%NKPTS

         NKI=KPOINTS_FULL%NEQUIV(NK)

         NGVEC=0
         NGVECI=0

         DO N3=1,GRID%NGZ_rd
            G3=(GRID%LPCTZ(N3)+WDES%VKPT(3,NK))
            NG3I=GRID%LPCTZ(N3)

            DO N2=1,GRID%NGY
               G2=(GRID%LPCTY(N2)+WDES%VKPT(2,NK))
               NG2I=GRID%LPCTY(N2)

               DO N1=1,GRID%NGX_rd
                  G1=(GRID%LPCTX(N1)+WDES%VKPT(1,NK))
                  NG1I=GRID%LPCTX(N1)

                  GIX= (G1*BI(1,1)+G2*BI(1,2)+G3*BI(1,3)) *TPI
                  GIY= (G1*BI(2,1)+G2*BI(2,2)+G3*BI(2,3)) *TPI
                  GIZ= (G1*BI(3,1)+G2*BI(3,2)+G3*BI(3,3)) *TPI
                  ENERGI=HSQDTM*((GIX**2)+(GIY**2)+(GIZ**2))

                  ! exclude some components for gamma-only version (C(G)=C*(-G))
                  IF (GRID%NGZ/=GRID%NGZ_rd) THEN
                     IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)<0) CYCLE
                     IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTX(N1)<0) CYCLE
                  ENDIF

                  ! is this part of the basis?
                  IF (ENERGI<WDES%ENMAX) THEN
                     ! apply symmetry operation taking point in IRZ into full BZ
                     NG1=NG1I*KPOINTS_FULL%IROTOP(1,1,NK)+NG2I*KPOINTS_FULL%IROTOP(2,1,NK)+NG3I*KPOINTS_FULL%IROTOP(3,1,NK)   
                     NG2=NG1I*KPOINTS_FULL%IROTOP(1,2,NK)+NG2I*KPOINTS_FULL%IROTOP(2,2,NK)+NG3I*KPOINTS_FULL%IROTOP(3,2,NK)   
                     NG3=NG1I*KPOINTS_FULL%IROTOP(1,3,NK)+NG2I*KPOINTS_FULL%IROTOP(2,3,NK)+NG3I*KPOINTS_FULL%IROTOP(3,3,NK)

                     NGVECI=NGVECI+1

                     ! does this data go on the local node?
                     IF (IGRIDIND(MP(NG1,NGX),MP(NG2,NGY),MP(NG3,NGZ))/=0) THEN
                        NGVEC=NGVEC+1
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         ! possibly increase the maximum number of plane wave coefficients stored locally
         NRPLWV=MAX(NRPLWV,NGVEC)
      ENDDO
 
      DEALLOCATE(IGRIDIND,IGRIDIND_MERGED)
#endif
    END SUBROUTINE COUNT_INDPW_FULL


!********************** SUBROUTINE ROTATE_WAVE_FFT *********************
!
!  this routine rotates a wavefunction (plane wave coefficients) 
!  from one k point in the IRBZ to another k'=S k and performs an FFT 
!  of the wavefunction to real space.
!  The trick is to know the index of the rotated plane wave G'= S G
!  this must be supplied in the NINDPW structure and is usually stored
!  in the WDES%NINDPW (see SET_INDPW_FULL).
!  Usually the wavefunction also needs to be multiplied by a phase factor
!  and possibly a conjugation is required (if time inversion is involved).
!  The difference to ROTATE_WAVE is, that ROTATE_WAVE uses index arrays
!  to store the wavefunctions in the proper position (as required by
!  GEN_LAYOUT and GEN_INDEX), whereas here the routine SET_INDPW_FULL
!  generates a "custom" layout different from the usual index arrays
!  used in VASP:
!  for k-vectors not located in the IRBZ the SET_INDPW_FULL sets
!  IGX, IGY and IGZ and NINDPW to G' = S G, 
!  where S is the symmetry operation taking k to k' = S k
!
!  Robin (09/2003) commented by gK
!
!***********************************************************************

  SUBROUTINE ROTATE_WAVE_FFT(COMM, &
       NPL,NINDPW,CPHASE,CR,C,C_ROT,GRID,LINV,LSHIFT)
    USE prec
    USE mpimy
    USE mgrid

    TYPE (communic) COMM
    TYPE (grid_3d)     GRID
    INTEGER NPL
    COMPLEX(q) :: C(NPL),C_ROT(NPL),CR(GRID%NPLWV),CPHASE(NPL)
    INTEGER       NINDPW(NPL)
    LOGICAL :: LINV,LSHIFT
    ! local
    INTEGER M

    ! first generate C_ROT
    ! i.e. the plane wave coefficients for the rotated wave function
    IF (LSHIFT .AND. LINV) THEN
       DO M=1,NPL
          C_ROT(M)=CONJG(C(M))*CPHASE(M)
       ENDDO
    ELSE IF (LSHIFT) THEN
       DO M=1,NPL
          C_ROT(M)=C(M)*CPHASE(M)
       ENDDO
    ELSE IF (LINV) THEN
       DO M=1,NPL
          C_ROT(M)=CONJG(C(M))
       ENDDO
    ELSE
       DO M=1,NPL
          C_ROT(M)=C(M)
       ENDDO
    ENDIF

!gK:  25.02.2012 GRID%NPLWV -> GRID%RC%NP (number of data points in rec. space)
    DO M=1,GRID%RC%NP
       CR(M)=(0.0_q,0.0_q)
    ENDDO

!DIR$ IVDEP
!OCL NOVREC
    DO M=1,NPL
       CR(NINDPW(M))=C_ROT(M)
    ENDDO
    ! in place fft from reciprocal to real space
    CALL FFT3D(CR,GRID,1)
    RETURN
  END SUBROUTINE ROTATE_WAVE_FFT


  SUBROUTINE ROTATE_WAVE_FFT_MRG(COMM, NPL_ORIG, NPL_POS, NPL_TOT, NG_INDEX, & 
       NPL,NINDPW,CPHASE,CR,C,C_ROT,GRID,LINV,LSHIFT)
#ifdef _OPENACC
    USE mopenacc_struct_def
#endif
    USE prec
    USE mpimy
    USE mgrid

    TYPE (communic) COMM
    TYPE (grid_3d)     GRID
    INTEGER NPL_ORIG  !  original number of plane waves (sender) 
    INTEGER NPL_POS   !  sum of NPL_ORIG up to (but not including) the current node
    INTEGER NPL_TOT   !  total number of plane waves
    INTEGER NG_INDEX(NPL_ORIG) ! position for PW coefficient in parallel version
    INTEGER NPL
    COMPLEX(q) ::  C(NPL),C_ROT(NPL),CR(GRID%MPLWV),CPHASE(NPL)
    INTEGER    ::  NINDPW(NPL)
    LOGICAL :: LINV,LSHIFT
    COMPLEX(q) :: C_TMP(NPL_TOT)
    ! local
    INTEGER M


    ! first generate C_ROT
    ! i.e. the plane wave coefficients for the rotated wave function
    IF (LSHIFT .AND. LINV) THEN
!$ACC PARALLEL LOOP PRESENT(C_ROT,C,CPHASE) __IF_ASYNC__
       DO M=1,NPL_ORIG
          C_ROT(M)=CONJG(C(M))*CPHASE(M)
       ENDDO
    ELSE IF (LSHIFT) THEN
!$ACC PARALLEL LOOP PRESENT(C_ROT,C,CPHASE) __IF_ASYNC__
       DO M=1,NPL_ORIG
          C_ROT(M)=C(M)*CPHASE(M)
       ENDDO
    ELSE IF (LINV) THEN
!$ACC PARALLEL LOOP PRESENT(C_ROT,C) __IF_ASYNC__
       DO M=1,NPL_ORIG
          C_ROT(M)=CONJG(C(M))
       ENDDO
    ELSE
!$ACC PARALLEL LOOP PRESENT(C_ROT,C) __IF_ASYNC__
       DO M=1,NPL_ORIG
          C_ROT(M)=C(M)
       ENDDO
    ENDIF
 
!$ACC PARALLEL LOOP PRESENT(CR) __IF_ASYNC__
    DO M=1,GRID%MPLWV
       CR(M)=(0.0_q,0.0_q)
    ENDDO

! now merge from all nodes using global sum
!$ACC ENTER DATA CREATE(C_TMP) __IF_ASYNC__
!$ACC KERNELS PRESENT(C_TMP) __IF_ASYNC__
    C_TMP=0
!$ACC END KERNELS
!$ACC PARALLEL LOOP PRESENT(C_TMP,C_ROT) __IF_ASYNC__
    DO M=1,NPL_ORIG
       C_TMP(NPL_POS+M-1)=C_ROT(M)
    ENDDO

    CALLMPI( M_sum_z(COMM, C_TMP , NPL_TOT))
!DIR$ IVDEP
!OCL NOVREC
!$ACC PARALLEL LOOP PRESENT(CR,NINDPW,C_TMP,NG_INDEX) __IF_ASYNC__
    DO M=1,NPL
       CR(NINDPW(M))=C_TMP(NG_INDEX(M))
    ENDDO
!$ACC EXIT DATA DELETE(C_TMP) __IF_ASYNC__
    ! in place fft from reciprocal to real space
    CALL FFT3D(CR,GRID,1)
    RETURN
  END SUBROUTINE ROTATE_WAVE_FFT_MRG

!********************** SUBROUTINE ROTATE_WAVE     *********************
!
! rotate a single wavefunction from one k-points to another one
! two routines are available
! the second one adds the results
!
!***********************************************************************


  SUBROUTINE ROTATE_WAVE( NPL, C_ROT, C, CPHASE, NINDPW, LINV, LSHIFT)
    USE prec
    USE mpimy
    USE mgrid
    IMPLICIT NONE

    TYPE (grid_3d)     GRID
    INTEGER NPL
    COMPLEX(q) :: C(NPL),C_ROT(NPL),CPHASE(NPL)
    INTEGER       NINDPW(NPL)
    LOGICAL :: LINV, LSHIFT
    ! local
    INTEGER M

    IF (LSHIFT .AND. LINV) THEN
       DO M=1,NPL
          C_ROT(NINDPW(M))=CONJG(C(M))*CPHASE(M)
       ENDDO
    ELSE IF (LSHIFT) THEN
       DO M=1,NPL
          C_ROT(NINDPW(M))=C(M)*CPHASE(M)
       ENDDO
    ELSE IF (LINV) THEN
       DO M=1,NPL
          C_ROT(NINDPW(M))=CONJG(C(M))
       ENDDO
    ELSE
       DO M=1,NPL
          C_ROT(NINDPW(M))=C(M)
       ENDDO
    ENDIF
    
    RETURN
  END SUBROUTINE ROTATE_WAVE

  SUBROUTINE ROTATE_WAVE_ADD( NPL, C_ROT, C, CPHASE, NINDPW, LINV, LSHIFT, FAKT)
    USE prec
    USE mpimy
    USE mgrid
    IMPLICIT NONE

    TYPE (grid_3d)     GRID
    INTEGER NPL
    COMPLEX(qs) :: C(NPL)
    COMPLEX(q)  :: C_ROT(NPL),CPHASE(NPL)
    INTEGER       NINDPW(NPL)
    LOGICAL :: LINV, LSHIFT
    REAL(q) FAKT
    ! local
    INTEGER M

    IF (LSHIFT .AND. LINV) THEN
       DO M=1,NPL
          C_ROT(NINDPW(M))=C_ROT(NINDPW(M))+CONJG(C(M))*CPHASE(M)*FAKT
       ENDDO
    ELSE IF (LSHIFT) THEN
       DO M=1,NPL
          C_ROT(NINDPW(M))=C_ROT(NINDPW(M))+C(M)*CPHASE(M)*FAKT
       ENDDO
    ELSE IF (LINV) THEN
       DO M=1,NPL
          C_ROT(NINDPW(M))=C_ROT(NINDPW(M))+CONJG(C(M))*FAKT
       ENDDO
    ELSE
       DO M=1,NPL
          C_ROT(NINDPW(M))=C_ROT(NINDPW(M))+C(M)*FAKT
       ENDDO
    ENDIF
    
    RETURN
  END SUBROUTINE ROTATE_WAVE_ADD

!***********************************************************************
! 
! back rotation of the wavefunction (inverse of ROTATE_WAVE_ADD)
! this is required to bring a wavefunction from the BZ to the IRBZ
! 
!***********************************************************************

  SUBROUTINE ROTATE_WAVE_BACK( NPL, C_ROT, C, CPHASE, NINDPW, LINV, LSHIFT)
    USE prec
    USE mpimy
    USE mgrid
    IMPLICIT NONE

    TYPE (grid_3d)     GRID
    INTEGER NPL
    COMPLEX(q)  :: C(NPL)
    COMPLEX(q)  :: C_ROT(NPL),CPHASE(NPL)
    INTEGER       NINDPW(NPL)
    LOGICAL :: LINV, LSHIFT
    ! local
    INTEGER M

    IF (LSHIFT .AND. LINV) THEN
       DO M=1,NPL
          C(M)=C(M)+CONJG(C_ROT(NINDPW(M)))*CPHASE(M)
       ENDDO
    ELSE IF (LSHIFT) THEN
       DO M=1,NPL
          C(M)=C(M)+C_ROT(NINDPW(M))*CONJG(CPHASE(M))
       ENDDO
    ELSE IF (LINV) THEN
       DO M=1,NPL
          C(M)=C(M)+CONJG(C_ROT(NINDPW(M)))
       ENDDO
    ELSE
       DO M=1,NPL
          C(M)=C(M)+C_ROT(NINDPW(M))
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE ROTATE_WAVE_BACK
 
#ifndef oldsym
  ! slight modification of rotate_wave_back, removed the c=c+x to c=x
  SUBROUTINE ROTATE_WAVE_BACK_NOADD( NPL, C_ROT, C, CPHASE, NINDPW, LINV, LSHIFT)
    USE prec
    USE mpimy
    USE mgrid
    IMPLICIT NONE

    TYPE (grid_3d)     GRID
    INTEGER NPL
    COMPLEX(q)  :: C(NPL)
    COMPLEX(q)  :: C_ROT(NPL),CPHASE(NPL)
    INTEGER       NINDPW(NPL)
    LOGICAL :: LINV, LSHIFT
    ! local
    INTEGER M

    IF (LSHIFT .AND. LINV) THEN
       DO M=1,NPL
          C(M)=CONJG(C_ROT(NINDPW(M)))*CPHASE(M)
       ENDDO
    ELSE IF (LSHIFT) THEN
       DO M=1,NPL
          C(M)=C_ROT(NINDPW(M))*CONJG(CPHASE(M))
       ENDDO
    ELSE IF (LINV) THEN
       DO M=1,NPL
          C(M)=CONJG(C_ROT(NINDPW(M)))
       ENDDO
    ELSE
       DO M=1,NPL
          C(M)=C_ROT(NINDPW(M))
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE ROTATE_WAVE_BACK_NOADD
#endif

!***********************************************************************
! 
! rotate a non-local potential or Hamilton matrix
! from one q-point to another
! according the VASP convention the second index transforms like G+q
! whereas the the first one transforms -(G+q) 
! [compare ADD_XI comments]
!
! what I have not yet tested is whether space group operations
! and the corresponding phase shifts (CPHASE) are properly 
! programmed
! I think I got it right, but maybe CPHASE needs to be
! conjugated in both lines
!
!***********************************************************************


  SUBROUTINE ROTATE_WPOT( NPL, W_ROT, W, NDIM, CPHASE, NINDPW, LINV, LSHIFT)
    USE prec
    USE mpimy
    USE mgrid
    IMPLICIT NONE

    TYPE (grid_3d)     GRID
    INTEGER NPL, NDIM
    COMPLEX(q) :: W(NDIM,NPL),W_ROT(NDIM,NPL)
    COMPLEX(q) :: CPHASE(NPL)
    INTEGER       NINDPW(NPL)
    LOGICAL :: LINV, LSHIFT
    ! local
    INTEGER M
    COMPLEX(q),ALLOCATABLE :: W_TMP(:,:)

    ALLOCATE(W_TMP(NDIM, NPL))

! first build the complex conjugated if required
    IF (LINV) THEN
       W_ROT=CONJG(W)
    ELSE
       W_ROT=W
    ENDIF
! second apply phase factors
    IF (LSHIFT) THEN
       DO M=1,NPL
          W_TMP(M,1:NPL)=W_ROT(M,1:NPL)*CONJG(CPHASE(M))
       ENDDO
       DO M=1,NPL
          W_ROT(1:NPL,M)=W_TMP(1:NPL,M)*CPHASE(M)
       ENDDO
    ENDIF
! third re-index rows and columns
    DO M=1,NPL
       W_TMP(NINDPW(M),1:NPL)=W_ROT(M,1:NPL)
    ENDDO
    DO M=1,NPL
       W_ROT(1:NPL,NINDPW(M))=W_TMP(1:NPL,M)
    ENDDO

    DEALLOCATE(W_TMP)
    
    RETURN
  END SUBROUTINE ROTATE_WPOT

!***********************************************************************
! 
! rotate a non-local potential or Hamilton matrix
! from one q-point in the BZ back into the IRZ
! this is the the inverse of the routine ROTATE_WPOT
! input:  W_ROT
! output: W
!
!***********************************************************************

  SUBROUTINE ROTATE_WPOT_BACK( NPL, W_ROT, W, NDIM, CPHASE, NINDPW, LINV, LSHIFT)
    USE prec
    USE mpimy
    USE mgrid
    IMPLICIT NONE

    TYPE (grid_3d)     GRID
    INTEGER NPL, NDIM
    COMPLEX(q) :: W(NDIM,NPL),W_ROT(NDIM,NPL)
    COMPLEX(q) :: CPHASE(NPL)
    INTEGER       NINDPW(NPL)
    LOGICAL :: LINV, LSHIFT
    ! local
    INTEGER M
    COMPLEX(q) :: W_TMP(NDIM,NPL)

! first re-index rows and columns
    DO M=1,NPL
       W_TMP(M,1:NPL)=W_ROT(NINDPW(M),1:NPL)
    ENDDO
    DO M=1,NPL
       W(1:NPL,M)=W_TMP(1:NPL,NINDPW(M))
    ENDDO
! second apply phase factors
    IF (LSHIFT) THEN
       DO M=1,NPL
          W_TMP(M,1:NPL)=W(M,1:NPL)*CPHASE(M)
       ENDDO
       DO M=1,NPL
          W(1:NPL,M)=W_TMP(1:NPL,M)*CONJG(CPHASE(M))
       ENDDO
    ENDIF
! final build the complex conjugated if required
    IF (LINV) THEN
       W=CONJG(W)
    ENDIF
    
    RETURN
  END SUBROUTINE ROTATE_WPOT_BACK

!***********************************************************************
! 
! rotate a vector (G+q) from one q-point to another
!
!***********************************************************************
  SUBROUTINE ROTATE_WPOT_DIAG( NPL, W_ROT, W, CPHASE, NINDPW, LINV, LSHIFT)
    USE prec
    USE mpimy
    USE mgrid
    IMPLICIT NONE

    TYPE (grid_3d)     GRID
    INTEGER NPL
    COMPLEX(q) :: W(NPL),W_ROT(NPL)
    COMPLEX(q) :: CPHASE(NPL)
    INTEGER       NINDPW(NPL)
    LOGICAL :: LINV, LSHIFT
    ! local
    INTEGER M
    COMPLEX(q),ALLOCATABLE :: W_TMP(:)

    ALLOCATE(W_TMP(NPL))

! first build the complex conjugated if required
    IF (LINV) THEN
       W_ROT=CONJG(W)
    ELSE
       W_ROT=W
    ENDIF
! second apply phase factors
    IF (LSHIFT) THEN
       DO M=1,NPL
          W_TMP(M)=W_ROT(M)*CONJG(CPHASE(M))
       ENDDO
       DO M=1,NPL
          W_ROT(M)=W_TMP(M)*CPHASE(M)
       ENDDO
    ENDIF
! third re-index rows and columns
    DO M=1,NPL
       W_TMP(NINDPW(M))=W_ROT(M)
    ENDDO
    DO M=1,NPL
       W_ROT(NINDPW(M))=W_TMP(M)
    ENDDO

    DEALLOCATE(W_TMP)
    
  END SUBROUTINE ROTATE_WPOT_DIAG

!***********************************************************************
! 
! subroutine to rotate one vector VEC (storing the wavefunction
! character for a specific atom)
! result is returned in ROTVEC
! 
!***********************************************************************

  SUBROUTINE ROTATE_VECTOR(LINV,VEC,ROTVEC,MMAX,LMAX,SL,NPS,LPS)
!$ACC ROUTINE VECTOR
    USE prec
    USE pseudo
    IMPLICIT NONE

    LOGICAL LINV                        ! inversion required (i.e. conjugation of final vector)
    GDEF :: VEC(*)                      ! initial vector
    GDEF :: ROTVEC(*)                   ! final vector
    INTEGER MMAX                        ! first and second dimension of SL
    INTEGER LMAX                        ! final dimension of SL
    REAL(q) :: SL(MMAX,MMAX,0:LMAX)     ! rotation matrix (allways symmetric)
    INTEGER :: NPS                      ! number of pseudopotential channels
    INTEGER :: LPS(NPS)                 ! l-quantum number of the pseudopotential channels
    ! local variables
    GDEF :: GTMP
    INTEGER CHANNEL,CHANNELS,IND,L,M,MP

    ! left hand transformation
    IND=0

    DO CHANNEL=1,NPS
       ! l-qantum number of this channel
       L=LPS(CHANNEL)
       ! rotate this l-block
!$ACC LOOP VECTOR PRIVATE(GTMP,MP)
       DO M=1,2*L+1
          GTMP=0
!$ACC LOOP REDUCTION(+:GTMP)
          DO MP=1,2*L+1
             GTMP=GTMP+SL(M,MP,L)*VEC(IND+MP)
          ENDDO
#ifdef _OPENACC
          IF (LINV) THEN
             ROTVEC(IND+M)=GCONJG(GTMP)
          ELSE
#endif
             ROTVEC(IND+M)=GTMP
#ifdef _OPENACC
          ENDIF
#endif
       ENDDO

       IND=IND+2*L+1
    ENDDO
#ifndef _OPENACC
    IF (LINV) THEN
       ROTVEC(1:IND)=GCONJG(ROTVEC(1:IND))
    ENDIF
#endif
  END SUBROUTINE ROTATE_VECTOR


  SUBROUTINE ROTATE_VECTOR_ADD(LINV,VEC,ROTVEC,MMAX,LMAX,SL,P,FAKT)
    USE prec
    USE pseudo
    IMPLICIT NONE
    
    LOGICAL LINV                        ! inversion required (i.e. conjugation of final vector)
    GDEFS:: VEC(*)                      ! initial vector
    GDEF :: ROTVEC(*)                   ! final vector
    INTEGER MMAX                        ! first and second dimension of SL
    INTEGER LMAX                        ! final dimension of SL
    REAL(q) :: SL(MMAX,MMAX,0:LMAX)     ! rotation matrix (allways symmetric)
    TYPE (potcar) P                     ! pseudopotential descriptor
    REAL(q) FAKT
 ! local variables
    INTEGER CHANNEL,CHANNELS,IND,L,M,MP
    
    CHANNELS=P%LMAX
    ! left hand transformation
    IND=0
      
    IF (LINV) THEN
       DO CHANNEL=1,CHANNELS
          ! l-qantum number of this channel
          L=P%LPS(CHANNEL)
          ! rotate this l-block
          DO M=1,(2*L+1)
             DO MP=1,(2*L+1)
                ROTVEC(IND+M)=ROTVEC(IND+M)+GCONJG(SL(M,MP,L)*VEC(IND+MP))*FAKT
             ENDDO
          ENDDO
          
          IND=IND+(2*L+1)
       ENDDO
    ELSE
       DO CHANNEL=1,CHANNELS
          ! l-qantum number of this channel
          L=P%LPS(CHANNEL)
          ! rotate this l-block
          DO M=1,(2*L+1)
             DO MP=1,(2*L+1)
                ROTVEC(IND+M)=ROTVEC(IND+M)+SL(M,MP,L)*VEC(IND+MP)*FAKT
             ENDDO
          ENDDO
          
          IND=IND+(2*L+1)
       ENDDO
    ENDIF

  END SUBROUTINE ROTATE_VECTOR_ADD
