
#include "symbol.inc"

!***********************************************************************
!
! this module implements a projection onto k-points
!
!***********************************************************************

MODULE mkproj
  USE prec
  USE poscar
  USE nonl_high

  IMPLICIT none

  LOGICAL, SAVE :: LKPROJ
  REAL(q), SAVE :: KPROJ_THRESHOLD

CONTAINS
!**********************************************************************
!
! read all variables related to exchange correlation treatment
! from the INCAR file
! this set both the variables in the setex module and
! in the local fock module
!
!**********************************************************************

  SUBROUTINE KPROJ_READER(IU5, IU0 )
      USE reader_tags
   
      IMPLICIT NONE
      INTEGER IU5, IU0
  ! local
      INTEGER IDUM, N, IERR
      REAL(q) RDUM
      COMPLEX(q) CDUM
      LOGICAL LOPEN,LDUM
      CHARACTER (1) :: CHARAC

      CALL OPEN_INCAR_IF_FOUND(IU5, LOPEN)

      LKPROJ=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LKPROJ', LKPROJ, IERR, WRITEXMLINCAR)
      
      KPROJ_THRESHOLD=1.E-2_q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'KPROJ_THRESHOLD', KPROJ_THRESHOLD, IERR, WRITEXMLINCAR)

      CALL CLOSE_INCAR_IF_FOUND(IU5)

    END SUBROUTINE KPROJ_READER

!***********************************************************************
! search a particular k-point and return
! the equivalent k-point in the array, if not found return -1
! this is equal to the routine KPOINT_IN_FULL_GRID, except 
! that this routine returns -1 and no error if point is not found
! in the KPOINTS_F array
!
!***********************************************************************

  FUNCTION FIND_KPOINT_IN_PC(N1,N2,N3,VKPT,PNK)
    !USE sym_prec
    USE lattice
    INTEGER :: FIND_KPOINT_IN_PC
    REAL(q) N1,N2,N3
    REAL(q),PARAMETER :: TINY=1E-8_q
    REAL(q) VKPT(:,:)
    INTEGER PNK
    ! local
    INTEGER IND
    DO IND=1,PNK
       IF ( &
         (ABS(MOD(N1-VKPT(1,IND)+10.5_q,1._q)-0.5_q)<TINY) .AND. &
         (ABS(MOD(N2-VKPT(2,IND)+10.5_q,1._q)-0.5_q)<TINY) .AND. &
         (ABS(MOD(N3-VKPT(3,IND)+10.5_q,1._q)-0.5_q)<TINY)) EXIT
    ENDDO

    IF (IND>PNK) THEN
       ! no kpoint found, set nk=-1
       IND=-1
    ENDIF

    FIND_KPOINT_IN_PC=IND
  END FUNCTION FIND_KPOINT_IN_PC

!**********************************************************************
!
!**********************************************************************
    SUBROUTINE KPROJ(IU5, IU0, IU6, GRID, LATT_CUR, W, SYM, CQIJ)
      USE lattice
      USE constant
      USE msymmetry
      USE base
      USE main_mpi
      USE vaspxml
      USE wave_high
      USE POSCAR, ONLY: RD_POSCAR_ALT
      USE pead, ONLY : LPEAD_SYM_RED
      USE tutor, ONLY: vtutor, isAlert, KProjection

      IMPLICIT NONE
      INTEGER IU5, IU0, IU6  ! units for IO
      TYPE (wavespin)    :: W          ! wavefunction
      TYPE (latt)        :: LATT_CUR 
      TYPE (grid_3d)     :: GRID
      TYPE (symmetry)    :: SYM
      OVERLAP            :: CQIJ(:,:,:,:) ! overlap operator
      ! local
      TYPE (latt)        :: LATT_PRIM
      TYPE (type_info)   :: T_INFO_PRIM
      TYPE (dynamics)    :: DYN_PRIM
      INTEGER :: N1, N2, N3, IND
      INTEGER :: NK, KPT, NGVECTOR, NGVECTORM, NKPTS_PRIM, NKPTS_IRZ
      INTEGER :: ISP, N, I
      REAL(q) :: TEMP
      REAL(q) :: G1, G2, G3, GIX, GIY, GIZ, GX, GY, GZ
      REAL(q) :: R1, R2, R3, K1, K2, K3, RS1, RS2, RS3
      REAL(q), ALLOCATABLE :: VKPT(:,:)      ! k-points in the BZ of primtive cell
      REAL(q), ALLOCATABLE :: VKPT_IRZ(:,:)   
      REAL(q), ALLOCATABLE :: WTKPT_IRZ(:)
      REAL(q), ALLOCATABLE :: KAPPA(:,:,:,:)
      INTEGER, ALLOCATABLE :: NKPTS_PRIM_INDEX(:,:), INDEX_IN_IRZ(:)
      TYPE (wavedes1)    WDES1
      INTEGER  ISYMOP,NROT,IGRPOP,NROTK,INVMAP,NPCELL,NB_GLOBAL,ISPINOR
      REAL(q)  GTRANS, AP, SUM_KAPPA
      COMMON /SYMM/   ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

      CALL KPROJ_READER(IU5, IU0)
      ! quick return if possible
      IF (.NOT. LKPROJ) RETURN

#ifdef MPI
      IF (W%WDES%COMM_INB%NCPU /=1) THEN
         CALL vtutor%write(isAlert, KProjection)
         RETURN
      ENDIF
#endif

      IF (IU0>=0) THEN
         WRITE(IU0,*) 'start k-points projection onto the first Brillouin zone of primitive cell.'
         WRITE(IU0,*) 'reading POSCAR.prim'
      ENDIF

      CALL RD_POSCAR_ALT("POSCAR.prim", LATT_PRIM, T_INFO_PRIM, DYN_PRIM, IU0, IU6)

      ALLOCATE(T_INFO_PRIM%ATOMOM(3*T_INFO_PRIM%NIOND))
      T_INFO_PRIM%ATOMOM=0

      IF (SYM%ISYM>0) THEN
         CALL INISYM(LATT_PRIM%A, DYN_PRIM%POSION, DYN_PRIM%VEL, T_INFO_PRIM%LSFOR, &
              T_INFO_PRIM%LSDYN,LPEAD_SYM_RED(),T_INFO_PRIM%NTYP,T_INFO_PRIM%NITYP,T_INFO_PRIM%NIOND, &
              SYM%PTRANS,SYM%NROT,SYM%NPTRANS,SYM%ROTMAP, &
              SYM%TAU,SYM%TAUROT,SYM%WRKROT, &
              SYM%INDROT,T_INFO_PRIM%ATOMOM,W%WDES%SAXIS,SYM%MAGROT,W%WDES%NCDIJ,IU6)
      ENDIF

      DEALLOCATE(T_INFO_PRIM%ATOMOM)

      NGVECTORM=1
      DO NK=1,W%WDES%NKPTS
         NGVECTOR=W%WDES%NGVECTOR(NK)
         NGVECTORM=MAX(NGVECTOR,NGVECTORM)
         !IF (IU0>=0) write(IU0,*) 'NK,NGVECTORM=', NK, NGVECTORM
      ENDDO

      ALLOCATE(VKPT(3,NGVECTORM*W%WDES%NKPTS),NKPTS_PRIM_INDEX(NGVECTORM,W%WDES%NKPTS))

      VKPT=0
      NKPTS_PRIM=0
!=======================================================================
! main loop over all special points
!=======================================================================
      kpoint: DO NK=1,W%WDES%NKPTS
         NGVECTOR=W%WDES%NGVECTOR(NK)
         CALL SETWDES(W%WDES,WDES1,NK)
!=======================================================================
! loop over all G-vectors in the basis at this k-point
!=======================================================================
         DO IND=1,NGVECTOR

            N1=MOD(W%WDES%IGX(IND,NK)+GRID%NGX,GRID%NGX)+1
            N2=MOD(W%WDES%IGY(IND,NK)+GRID%NGY,GRID%NGY)+1 
            N3=MOD(W%WDES%IGZ(IND,NK)+GRID%NGZ,GRID%NGZ)+1

! get G-vector of respective k-point and coefficient
            G1=(GRID%LPCTX(N1)+W%WDES%VKPT(1,NK))
            G2=(GRID%LPCTY(N2)+W%WDES%VKPT(2,NK))
            G3=(GRID%LPCTZ(N3)+W%WDES%VKPT(3,NK))
            GIX=(G1*LATT_CUR%B(1,1)+G2*LATT_CUR%B(1,2)+G3*LATT_CUR%B(1,3))
            GIY=(G1*LATT_CUR%B(2,1)+G2*LATT_CUR%B(2,2)+G3*LATT_CUR%B(2,3))
            GIZ=(G1*LATT_CUR%B(3,1)+G2*LATT_CUR%B(3,2)+G3*LATT_CUR%B(3,3))

!=======================================================================
! bring the G vectors into the first BZ of primitive cell
!=======================================================================
            R1=(GIX*LATT_PRIM%A(1,1)+GIY*LATT_PRIM%A(2,1)+GIZ*LATT_PRIM%A(3,1))
            R2=(GIX*LATT_PRIM%A(1,2)+GIY*LATT_PRIM%A(2,2)+GIZ*LATT_PRIM%A(3,2))
            R3=(GIX*LATT_PRIM%A(1,3)+GIY*LATT_PRIM%A(2,3)+GIZ*LATT_PRIM%A(3,3))

            K1=MOD(R1+1000.5_q,1._q)-0.5_q
            K2=MOD(R2+1000.5_q,1._q)-0.5_q
            K3=MOD(R3+1000.5_q,1._q)-0.5_q

!=======================================================================
! find the corresponding K and make a list of k-vectors in prim cell
!=======================================================================
            KPT=FIND_KPOINT_IN_PC(K1,K2,K3,VKPT,NKPTS_PRIM)

            ! k-point not found in the existing list
            ! add the k-point to the list
            IF (KPT<0) THEN
               NKPTS_PRIM=NKPTS_PRIM+1
               VKPT(1,NKPTS_PRIM)=K1
               VKPT(2,NKPTS_PRIM)=K2
               VKPT(3,NKPTS_PRIM)=K3
               KPT=NKPTS_PRIM
            ENDIF
            NKPTS_PRIM_INDEX(IND,NK)=KPT

         ENDDO
      ENDDO kpoint
      !IF (IU0>=0) WRITE(IU0,*) NKPTS_PRIM_INDEX
      IF (IU0>=0) WRITE(IU0,*) 'K-point list is generated..'
      IF (IU6>=0) WRITE(IU6,'(A,I6)') 'NKPTS_PRIM=', NKPTS_PRIM
      
      ALLOCATE(VKPT_IRZ(3,NKPTS_PRIM), INDEX_IN_IRZ(NKPTS_PRIM),WTKPT_IRZ(NKPTS_PRIM))

      CALL IBZKPT_LIST(LATT_PRIM, VKPT, NKPTS_PRIM, VKPT_IRZ, NKPTS_IRZ, &
               INDEX_IN_IRZ, SYM%ROTMAP, SYM%MAGROT, SYM%ISYM, IU6, IU0)
      
      WTKPT_IRZ=0
      DO NK=1,NKPTS_PRIM
         WTKPT_IRZ(INDEX_IN_IRZ(NK))=WTKPT_IRZ(INDEX_IN_IRZ(NK))+1
      ENDDO

      IF (IU0>=0) WRITE(IU0,*) 'Symmetry operation is done...'
      !IF (IU0>=0) WRITE(IU0,*) INDEX_IN_IRZ

      ALLOCATE(KAPPA(W%WDES%NB_TOT,W%WDES%NKPTS,W%WDES%ISPIN,NKPTS_IRZ))

      KAPPA=0
!=======================================================================
! loop over k-points and bands
!=======================================================================
      spin:    DO ISP=1,W%WDES%ISPIN
      kpoints: DO NK=1,W%WDES%NKPTS
         CALL SETWDES(W%WDES,WDES1,NK)
         NGVECTOR=W%WDES%NGVECTOR(NK)
      band: DO NB_GLOBAL=1,W%WDES%NB_TOT
         N=NB_LOCAL(NB_GLOBAL,WDES1)
         IF(N==0) CYCLE band
!=======================================================================
! take one wavefunction loop over all G vector and add Kappa 
! to the array Kappa that has the size of the number of k-points
!=======================================================================
         DO ISPINOR=0,WDES1%NRSPINORS-1
         DO IND=1,NGVECTOR
            KPT=INDEX_IN_IRZ(NKPTS_PRIM_INDEX(IND,NK))
            IF (KPT==0) THEN
               CALL vtutor%bug("internal error in KPROJ: G vector was not assigned to k-vector " // &
                  str(W%WDES%IGX(IND,NK)) // " " // str(W%WDES%IGY(IND,NK)) // " " // &
                  str(W%WDES%IGZ(IND,NK)), __FILE__, __LINE__)
            ENDIF
            TEMP=REAL(W%CPTWFP(IND+ISPINOR*NGVECTOR,N,NK,ISP)*CONJG(W%CPTWFP(IND+ISPINOR*NGVECTOR,N,NK,ISP)))
            KAPPA(NB_GLOBAL,NK,ISP,KPT)=KAPPA(NB_GLOBAL,NK,ISP,KPT)+TEMP
         ENDDO
         ENDDO
      ENDDO band
      ENDDO kpoints
      ENDDO spin
      CALLMPI( M_sum_d( W%WDES%COMM, KAPPA(1,1,1,1), SIZE(KAPPA)) )

!     IF (IU6>=0) THEN
!        WRITE(IU6,*) 'KPROJ_THRESHOLD=', KPROJ_THRESHOLD
!        WRITE(IU6,*) 
!        WRITE(IU6,*) '|K1|, eigenvalues, KAPPA'
!        WRITE(IU6,*) 
!     ENDIF
!     DO ISP=1,W%WDES%ISPIN
!     DO NK=1,W%WDES%NKPTS
!     DO KPT=1,NKPTS_IRZ
!     DO N=1,W%WDES%NB_TOT
!        SUM_KAPPA=SUM(KAPPA(N,NK,ISP,1:NKPTS_IRZ))
!        KAPPA(N,NK,ISP,1:NKPTS_IRZ)=KAPPA(N,NK,ISP,1:NKPTS_IRZ)/SUM_KAPPA
!
!     ! for Gamma-L direction 
!        IF (ABS(VKPT_IRZ(3,KPT))==ABS(VKPT_IRZ(2,KPT)) .AND. &
!            ABS(VKPT_IRZ(1,KPT))==ABS(VKPT_IRZ(2,KPT))) THEN
!           IF (KAPPA(N,NK,ISP,KPT)>KPROJ_THRESHOLD) THEN
!              IF (IU6>=0) THEN
!                 WRITE(IU6,'(3F12.6)') ABS(VKPT_IRZ(2,KPT)), &
!                                       REAL(W%CELTOT(N,NK,ISP)), &
!                                       KAPPA(N,NK,ISP,KPT)
!              ENDIF
!           ENDIF
!        ELSEIF (ABS(VKPT_IRZ(3,KPT))>1.E-3_q .AND. &
!            ABS(VKPT_IRZ(1,KPT))<1.E-4_q .AND. &
!            ABS(VKPT_IRZ(2,KPT))<1.E-4_q) THEN
!           IF (KAPPA(N,NK,ISP,KPT)>KPROJ_THRESHOLD) THEN
!              IF (IU6>=0) THEN
!                 WRITE(IU6,'(3F12.6)') ABS(VKPT_IRZ(3,KPT)), &
!                                       REAL(W%CELTOT(N,NK,ISP)), &
!                                       KAPPA(N,NK,ISP,KPT)
!              ENDIF
!           ENDIF
!        ELSEIF (ABS(VKPT_IRZ(2,KPT))>1.E-3_q .AND. &
!            ABS(VKPT_IRZ(1,KPT))<1.E-4_q .AND. &
!            ABS(VKPT_IRZ(3,KPT))<1.E-4_q) THEN
!           IF (KAPPA(N,NK,ISP,KPT)>KPROJ_THRESHOLD) THEN
!              IF (IU6>=0) THEN
!                 WRITE(IU6,'(3F12.6)') ABS(VKPT_IRZ(2,KPT)), &
!                                       REAL(W%CELTOT(N,NK,ISP)), &
!                                       KAPPA(N,NK,ISP,KPT)
!              ENDIF
!           ENDIF
!        ELSEIF (ABS(VKPT_IRZ(1,KPT))>1.E-3_q .AND. &
!            ABS(VKPT_IRZ(2,KPT))<1.E-4_q .AND. &
!            ABS(VKPT_IRZ(3,KPT))<1.E-4_q) THEN
!           IF (KAPPA(N,NK,ISP,KPT)>KPROJ_THRESHOLD) THEN
!              IF (IU6>=0) THEN
!                 WRITE(IU6,'(3F12.6)') ABS(VKPT_IRZ(1,KPT)), &
!                                       REAL(W%CELTOT(N,NK,ISP)), &
!                                       KAPPA(N,NK,ISP,KPT)
!              ENDIF
!           ENDIF
!        ENDIF
!     !! for Gamma-L direction 
!
!     ! for Gamma-X direction 
!     !   IF (ABS(VKPT_IRZ(3,KPT))<1.E-3_q .AND. &
!     !       ABS(VKPT_IRZ(1,KPT)-VKPT_IRZ(2,KPT))<1.E-3_q) THEN
!     !      IF (KAPPA(N,NK,ISP,KPT)>KPROJ_THRESHOLD) THEN
!     !         IF (IU6>=0) THEN
!     !            WRITE(IU6,'(3F12.6)') ABS(VKPT_IRZ(2,KPT)), &
!     !                                  REAL(W%CELTOT(N,NK,ISP)), &
!     !                                  KAPPA(N,NK,ISP,KPT)
!     !         ENDIF
!     !      ENDIF
!     !   ELSEIF (ABS(VKPT_IRZ(2,KPT))<1.E-3_q .AND. &
!     !       ABS(VKPT_IRZ(1,KPT)-VKPT_IRZ(3,KPT))<1.E-3_q) THEN
!     !      IF (KAPPA(N,NK,ISP,KPT)>KPROJ_THRESHOLD) THEN
!     !         IF (IU6>=0) THEN
!     !            WRITE(IU6,'(3F12.6)') ABS(VKPT_IRZ(3,KPT)), &
!     !                                  REAL(W%CELTOT(N,NK,ISP)), &
!     !                                  KAPPA(N,NK,ISP,KPT)
!     !         ENDIF
!     !      ENDIF
!     !   ELSEIF (ABS(VKPT_IRZ(1,KPT))<1.E-3_q .AND. &
!     !       ABS(VKPT_IRZ(2,KPT)-VKPT_IRZ(3,KPT))<1.E-3_q) THEN
!     !      IF (KAPPA(N,NK,ISP,KPT)>KPROJ_THRESHOLD) THEN
!     !         IF (IU6>=0) THEN
!     !            WRITE(IU6,'(3F12.6)') ABS(VKPT_IRZ(3,KPT)), &
!     !                                  REAL(W%CELTOT(N,NK,ISP)), &
!     !                                  KAPPA(N,NK,ISP,KPT)
!     !         ENDIF
!     !      ENDIF
!     !   ENDIF
!     ! for Gamma-X direction 
!     ENDDO
!     ENDDO
!     ENDDO
!     ENDDO
      
      IF (IU6>=0) THEN
         OPEN(UNIT=110,FILE=DIR_APP(1:DIR_LEN)//'PRJCAR',STATUS='UNKNOWN')
         WRITE(110,'(A)') 'Basis vectors reciprocal space of POSCAR.prim (units of 2pi):'
         WRITE(110,'(3F14.7)') LATT_PRIM%B
         WRITE(110,'(/A,I6)') 'number of k-points in IBZ of POSCAR.prim:',NKPTS_IRZ
         WRITE(110,'(/A)') "             b1            b2            b3      weight"
         DO KPT=1,NKPTS_IRZ
            WRITE(110,'(I4,3F14.7,3X,I4)') KPT,VKPT_IRZ(:,KPT),INT(WTKPT_IRZ(KPT))
         ENDDO

         DO ISP=1,W%WDES%ISPIN
            WRITE(110,'(/A,I4)') 'spin component:',ISP
            DO NK=1,W%WDES%NKPTS
               WRITE(110,'("k-point (associated with POSCAR):",I6,2X,"vkpt:",3F14.7,2X,"weight:",F14.7)') &
              &   NK,W%WDES%VKPT(:,NK),W%WDES%WTKPT(NK)
               DO N=1,W%WDES%NB_TOT
                  WRITE(110,'("band:",I6,2X,"energy:",F14.7)') N,REAL(W%CELTOT(N,NK,ISP),q)
                  I=0
                  DO KPT=1,NKPTS_IRZ
                     WRITE(110,'(E15.7)',ADVANCE='No') KAPPA(N,NK,ISP,KPT)
                     I=I+1; IF (MOD(I,10)==0) WRITE(110,*)
                  ENDDO 
                  IF (MOD(I,10)/=0) WRITE(110,*)
               ENDDO
            ENDDO
         ENDDO 
         CLOSE(110)
      ENDIF

      CALL XML_TAG("kpoints", comment="kpoints in IRZ of POSCAR.prim")
      CALL XML_KPOINTS_LIST(VKPT_IRZ(:,1:NKPTS_IRZ), WTKPT_IRZ(1:NKPTS_IRZ))
      CALL XML_CLOSE_TAG("kpoints")
      
      CALL XML_KPROJ(KAPPA, W%CELTOT, W%FERTOT, W%WDES%NB_TOT, W%WDES%NKPTS, W%WDES%ISPIN, NKPTS_IRZ)

      DEALLOCATE(VKPT, NKPTS_PRIM_INDEX, VKPT_IRZ, WTKPT_IRZ, INDEX_IN_IRZ)
      DEALLOCATE(KAPPA)

    END SUBROUTINE KPROJ


!******************** SUBROUTINE IBZKPT_LIST ***************************
!
! subroutine IBZKPT_LIST is deduced from ibzkpt
! it determines a list of k-points inside the IRZ by applying all
! symmetry operations to a supplied list of k-points
!
! it returns:        an array linking the k-points in IRZ 
!    to k-points in the full Brillouin zone
!
!***********************************************************************

      SUBROUTINE IBZKPT_LIST(LATT_CUR, VKPT, NKPTS, VKPT_IRZ, NKPTS_IRZ, INDEX_IN_IRZ, & 
           ROTMAP, MAGROT, ISYM, IU6, IU0)
      USE lattice
      USE mkpoints
      USE main_mpi

      IMPLICIT NONE
! passed structures and variables
      TYPE (latt)           LATT_CUR 
      INTEGER ISYM
      INTEGER               IU6,IU0
      INTEGER ROTMAP(:,:,:)
      REAL(q) :: VKPT(:,:)       ! k-points in the Brillouin zone
      REAL(q) :: VKPT_IRZ(:,:)   ! k-points in the IRZ
      INTEGER :: NKPTS           ! total number of k-points
      INTEGER :: NKPTS_IRZ       ! number of k-points in the IRZ
      INTEGER :: INDEX_IN_IRZ(:) ! index from each k-point in the Brillouin zone into IRZ
      REAL(q) MAGROT(:,:)
! common symmetry variables
      INTEGER ISYMOP, NROT, IGRPOP, NROTK, INVMAP, NPCELL
      REAL(q) GTRANS, AP
      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL
! local variables and structures
      REAL(q) TINY, V(3), VR(3), VT(3), ROP(3,3)
      REAL(q) KPT_WEIGHT(NKPTS)
      INTEGER INVERS(9),IOP(3,3)
      INTEGER NK,NKF,NOP
      LOGICAL LINV
! external routine
      LOGICAL,EXTERNAL :: SGREQL
      LOGICAL,ALLOCATABLE :: KPOINT_ALREADY_FOUND(:)
! set data
      DATA TINY /1.E-6_q/, INVERS /-1,0,0,0,-1,0,0,0,-1/

      ALLOCATE(KPOINT_ALREADY_FOUND(NKPTS))

      KPOINT_ALREADY_FOUND=.FALSE.
      NKPTS_IRZ=0
      INDEX_IN_IRZ=0
      LINV=.FALSE.
!=======================================================================
! now do all point group operations with each k-point and check wether we
! generate a new k-point (in 1st BZ). If so, store it as well as the
! generating operation. By the way, check whether inversion is already one 
! of the sym ops.
!=======================================================================
      DO NK=1,NKPTS
       IF (KPOINT_ALREADY_FOUND(NK)) THEN
          CYCLE
       ELSE
! store the k-point
         NKPTS_IRZ=NKPTS_IRZ+1
         IF (NKPTS_IRZ>SIZE(VKPT_IRZ,2)) THEN
            CALL vtutor%bug("IBZKPT_LIST: internal error.. VKPT_IRZ is not sufficiently large", __FILE__, __LINE__)
         ENDIF
         INDEX_IN_IRZ(NK)=NKPTS_IRZ
         VKPT_IRZ(1,NKPTS_IRZ)=VKPT(1,NK)
         VKPT_IRZ(2,NKPTS_IRZ)=VKPT(2,NK)
         VKPT_IRZ(3,NKPTS_IRZ)=VKPT(3,NK)             

         DO NOP=1,NROTK
! test existence of inversion op
            IF (SGREQL(IGRPOP(1,1,NOP),INVERS)) LINV=.TRUE.
! copy symmetry op to real array
            ROP=IGRPOP(:,:,NOP)
! make new k-point and shift (for testing) it to 1st BZ
            VR(1)=VKPT(1,NK)
            VR(2)=VKPT(2,NK)
            VR(3)=VKPT(3,NK) 
            V(1)=VR(1)*ROP(1,1)+VR(2)*ROP(2,1)+VR(3)*ROP(3,1)
            V(2)=VR(1)*ROP(1,2)+VR(2)*ROP(2,2)+VR(3)*ROP(3,2)
            V(3)=VR(1)*ROP(1,3)+VR(2)*ROP(2,3)+VR(3)*ROP(3,3)
! bring the point to the primitive cell
            VT(1)=MOD(V(1)+6.5_q,1._q)-0.5_q
            VT(2)=MOD(V(2)+6.5_q,1._q)-0.5_q
            VT(3)=MOD(V(3)+6.5_q,1._q)-0.5_q
! test against all other k-points
            test1: DO NKF=1,NKPTS
               IF(( ABS(MOD(VT(1)-VKPT(1,NKF)+6.5,1._q)-0.5_q)<TINY) .AND. &
                  ( ABS(MOD(VT(2)-VKPT(2,NKF)+6.5,1._q)-0.5_q)<TINY) .AND. &
                  ( ABS(MOD(VT(3)-VKPT(3,NKF)+6.5,1._q)-0.5_q)<TINY)) THEN
                  INDEX_IN_IRZ(NKF)=NKPTS_IRZ
                  KPOINT_ALREADY_FOUND(NKF)=.TRUE.
                  EXIT test1
               ENDIF
            ENDDO test1
         ENDDO

!=======================================================================
! did not find LINV -> now we have to do it all over again with 
! all operators multiplied with INVERS
!=======================================================================
         IF (.NOT. LINV .AND. ISYM>=0) THEN
            DO NOP=1,NROTK
               ! apply inversion symmetry to form to get IOP
               CALL SGRPRD(INVERS,IGRPOP(1,1,NOP),IOP(1,1))
               ROP=IOP  ! copy symmetry op to real array
! make new k-point and shift it (for testing) to 1st BZ
               VR(1)=VKPT(1,NK)
               VR(2)=VKPT(2,NK)
               VR(3)=VKPT(3,NK) 
               V(1)=VR(1)*ROP(1,1)+VR(2)*ROP(2,1)+VR(3)*ROP(3,1)
               V(2)=VR(1)*ROP(1,2)+VR(2)*ROP(2,2)+VR(3)*ROP(3,2)
               V(3)=VR(1)*ROP(1,3)+VR(2)*ROP(2,3)+VR(3)*ROP(3,3)
! bring the point to the primitive cell
               VT(1)=MOD(V(1)+6.5_q,1._q)-0.5_q
               VT(2)=MOD(V(2)+6.5_q,1._q)-0.5_q
               VT(3)=MOD(V(3)+6.5_q,1._q)-0.5_q
! test against all other k-points
               test2: DO NKF=1,NKPTS
                  IF(( ABS(MOD(VT(1)-VKPT(1,NKF)+6.5,1._q)-0.5)<TINY) .AND. &
                     ( ABS(MOD(VT(2)-VKPT(2,NKF)+6.5,1._q)-0.5)<TINY) .AND. &
                     ( ABS(MOD(VT(3)-VKPT(3,NKF)+6.5,1._q)-0.5)<TINY)) THEN
                     INDEX_IN_IRZ(NKF)=NKPTS_IRZ
                     KPOINT_ALREADY_FOUND(NKF)=.TRUE.
                     EXIT test2
                  ENDIF
               ENDDO test2
            ENDDO
         ENDIF
        ENDIF
      ENDDO


      KPT_WEIGHT=0
      DO NK=1,NKPTS
         KPT_WEIGHT(INDEX_IN_IRZ(NK))=KPT_WEIGHT(INDEX_IN_IRZ(NK))+1
      ENDDO

      DEALLOCATE(KPOINT_ALREADY_FOUND)

      IF (IU6>=0) THEN
         WRITE(IU6,*)
         WRITE(IU6,*) 'reciprocal lattice vectors of the primitive cell:'
         WRITE(IU6,'(3X,3F13.9)') LATT_CUR%B
         WRITE(IU6,*)
         WRITE(IU6,*) 'number of k-points in IRZ ',NKPTS_IRZ
         WRITE(IU6,*)' k-points in reciprocal lattice and weights: '
         DO NK=1,NKPTS_IRZ
            WRITE(IU6,'(1X,3F12.8,F12.3)') VKPT_IRZ(:,NK), KPT_WEIGHT(NK)
         ENDDO
      ENDIF
    END SUBROUTINE IBZKPT_LIST


END MODULE mkproj
