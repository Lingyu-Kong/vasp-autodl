#include "symbol.inc"
      MODULE POSCAR
      USE prec
      USE poscar_struct_def
      IMPLICIT NONE

      CONTAINS

#ifdef VASP_HDF5
      ! extract potential information out of poscar line
      SUBROUTINE EXTRACT_POTCAR_INFO(INPLIN, NTYP, TYPE, TYPEF, SHA256)
        CHARACTER (255), INTENT(IN) :: INPLIN
        CHARACTER (LEN=2), INTENT(OUT) :: TYPE(10000)
        CHARACTER (LEN=64), INTENT(OUT) :: SHA256(10000)
        CHARACTER (LEN=20), INTENT(OUT) :: TYPEF(10000)
        INTEGER, INTENT(IN) :: NTYP

        CHARACTER (LEN=20) :: ELEMENT
        CHARACTER (255) :: INPWRK
        INTEGER :: START_INDEX, POS, NI

        INPWRK = TRIM(INPLIN)

        START_INDEX = 1
        DO NI=1,NTYP
!          WRITE(*,*) 'INPWRK: ', INPWRK
          POS = INDEX(INPWRK(START_INDEX:), ' ')
          ELEMENT = INPWRK(START_INDEX:POS-1)
          INPWRK = TRIM(ADJUSTL(INPWRK(POS:)))
!          WRITE(*,*) 'Found Element: ', ELEMENT
          ! check for SHA256
          POS = INDEX(ELEMENT, '/')
          IF (POS /= 0) THEN
            TYPEF(NI) = ELEMENT(:POS-1)
            SHA256(NI) = ELEMENT(POS+1:)
          ELSE
            TYPEF(NI) = ELEMENT
            SHA256(NI) = ""
          END IF
          TYPE(NI) = TYPEF(NI)(1:2)
          IF (TYPE(NI)(2:2) == '_') THEN
            TYPE(NI)(2:2) = ' '
          ENDIF
          !WRITE(*,*)TYPE(NI), TYPEF(NI), SHA256(NI)
        END DO
      END SUBROUTINE
#endif
!=======================================================================
! RCS:  $Id: poscar.F,v 1.1 2000/11/15 08:13:54 kresse Exp $
!
!  Read UNIT=15: POSCAR file scan for total number of ions
!  and number of types
!  only T_INFO%NITYP is allocated at this point
!  all other arrays are allocated in RD_POSCAR
!=======================================================================
      SUBROUTINE RD_POSCAR_HEAD(LATT_CUR, T_INFO, &
     &           NIOND,NIONPD, NTYPD,NTYPPD, IU0, IU6)
      USE prec
      USE lattice
      USE main_mpi
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
#ifdef VASP_HDF5
      USE vhdf5_base
#endif
      IMPLICIT NONE

      INTEGER NIOND,NIONPD,NTYPPD,NTYPD
      INTEGER IU0,IU6

      TYPE (latt)::       LATT_CUR
      TYPE (type_info) :: T_INFO
      INTEGER NITYP(10000)          ! hard limit 10000 ions :->
      CHARACTER (LEN=2) TYPE(10000) ! type information
! temporary varibales
      CHARACTER (1)    CHARAC
      CHARACTER (255)  INPLIN,INPWRK
      INTEGER        NI,I,NT,NSCALE
      REAL(q)        SCALEX,SCALEY,SCALEZ
      INTEGER, EXTERNAL :: NITEMS
#ifdef VASP_HDF5
      INTEGER :: IERR, NUMB, POS
      INTEGER(HID_T) :: IH5INPOSCARGROUP_ID
!      CHARACTER(LEN=:), ALLOCATABLE :: FORM
      CHARACTER (LEN=64) :: SHA256(10000)
      CHARACTER (LEN=20) :: TYPEF(10000)
#endif

#ifdef VASP_HDF5
      INQUIRE(FILE=DIR_APP(1:DIR_LEN)//'POSCAR',EXIST=POSCAR_FOUND)
      IF (.NOT. POSCAR_FOUND) THEN
        IF (HDF5_FOUND) THEN
          IERR = VH5_GROUP_OPEN(IH5ININPUTGROUP_ID, SUBGRP_POSCAR, IH5INPOSCARGROUP_ID)
        ELSE
          GOTO 1000
        ENDIF
      ENDIF
#endif
#ifdef VASP_HDF5
      IF (POSCAR_FOUND) THEN
#endif
! Now extract from file POSCAR how many ion types we have ...
      OPEN(UNIT=15,FILE=DIR_APP(1:DIR_LEN)//'POSCAR',STATUS='OLD',ERR=1000)
      READ(15,'(A1)',ERR=147,END=147) CHARAC
! one scaling parameter or one for x, y and z
      READ(15,'(A)',ERR=147,END=147) INPLIN

! how many words/data items? --> number of ion types on file POSCAR!
      NSCALE=NITEMS(INPLIN,INPWRK,.TRUE.,'F')
      IF (NSCALE==1) THEN
         READ(INPLIN,*) LATT_CUR%SCALE
         SCALEX=1
         SCALEY=1
         SCALEZ=1
      ELSE IF (NSCALE==3) THEN
         LATT_CUR%SCALE=1
         READ(INPLIN,*) SCALEX,SCALEY,SCALEZ
      ELSE
         CALL vtutor%error("ERROR: there must be 1 or 3 items on line 2 of POSCAR")
      ENDIF

      DO I=1,3
        READ(15,*,ERR=147,END=147) LATT_CUR%A(1,I),LATT_CUR%A(2,I),LATT_CUR%A(3,I)
      ENDDO
#ifdef VASP_HDF5
      ELSE IF (HDF5_FOUND) THEN
        LATT_CUR%SCALE=1
        IERR = VH5_READ(IH5INPOSCARGROUP_ID,"scale_x", SCALEX)
        IERR = VH5_READ(IH5INPOSCARGROUP_ID,"scale_y", SCALEY)
        IERR = VH5_READ(IH5INPOSCARGROUP_ID,"scale_z", SCALEZ)
        IERR = VH5_READ(IH5INPOSCARGROUP_ID,"lattice_vectors", LATT_CUR%A); VH5_CHECK(IERR)
      ENDIF
#endif

      IF (LATT_CUR%SCALE<0._q) THEN
!----alternatively give a volume (=abs(scale)) and adjust the lengths of
!----the three lattice vectors to get the correct desired volume ... :
         CALL LATTIC(LATT_CUR)
         LATT_CUR%SCALE=(ABS(LATT_CUR%SCALE) &
     &                 / ABS(LATT_CUR%OMEGA))**(1._q/3._q)
      ENDIF
      
      LATT_CUR%A(1,:) =LATT_CUR%A(1,:)*SCALEX*LATT_CUR%SCALE
      LATT_CUR%A(2,:) =LATT_CUR%A(2,:)*SCALEY*LATT_CUR%SCALE
      LATT_CUR%A(3,:) =LATT_CUR%A(3,:)*SCALEZ*LATT_CUR%SCALE
         
      CALL LATTIC(LATT_CUR)

      IF (LATT_CUR%OMEGA<0) THEN
        CALL vtutor%error("ERROR: the triple product of the basis vectors is negative exchange two &
           &basis vectors")
      ENDIF
! 6th line, contains either the number of ions or their type
#ifdef VASP_HDF5
      IF (POSCAR_FOUND) THEN
#endif
      READ(15,'(A)',ERR=147,END=147) INPLIN
! how many words/data items? --> number of ion types on file POSCAR!
      READ(INPLIN,*,ERR=147,END=147) CHARAC
      IF (.NOT.(CHARAC>='0' .AND. CHARAC<='9')) THEN
         T_INFO%NTYP=NITEMS(INPLIN,INPWRK,.TRUE.,'A')
#ifdef VASP_HDF5
         CALL EXTRACT_POTCAR_INFO(INPLIN, T_INFO%NTYP, TYPE, TYPEF, SHA256)
#else
         READ(INPLIN,*,ERR=147,END=147) (TYPE(NI),NI=1,T_INFO%NTYP)
         DO NI=1,T_INFO%NTYP
            IF(TYPE(NI)(2:2) == '_') THEN
               TYPE(NI)(2:2) = ' '
            ENDIF
         ENDDO
#endif
         IF (IU0>=0) THEN
            WRITE(IU0,'(A,20A)') ' POSCAR found type information on POSCAR ',TYPE(1:T_INFO%NTYP)
         ENDIF

         READ(15,'(A)',ERR=147,END=147) INPLIN
         IF (T_INFO%NTYP/=NITEMS(INPLIN,INPWRK,.TRUE.,'I')) THEN
            CALL vtutor%error("ERROR: the type information is not consistent with the number of types")
         ENDIF
      ELSE
         T_INFO%NTYP=NITEMS(INPLIN,INPWRK,.TRUE.,'I')
         TYPE(1:T_INFO%NTYP+1)='  '
      ENDIF
! let me know how many ions
      READ(INPLIN,*,ERR=147,END=147) (NITYP(NI),NI=1,T_INFO%NTYP)
#ifdef VASP_HDF5
      ELSE
        IERR = VH5_READ(IH5INPOSCARGROUP_ID,"number_ion_types", NITYP, DIMENSIONS=T_INFO%NTYP)
        IERR = VH5_READ(IH5INPOSCARGROUP_ID,"ion_types", TYPEF)
        DO I=1,T_INFO%NTYP
          POS = INDEX(TYPEF(I), '_')
          IF (POS /= 0) THEN
            TYPE(I) = TYPEF(I)(1:POS-1)
          ELSE
            TYPE(I) = TYPEF(I)
          ENDIF
        ENDDO
        IERR = VH5_READ(IH5INPOSCARGROUP_ID,"ion_sha256", SHA256, DIMENSIONS=NUMB)
        IF (IERR == 0 .AND. IU0>=0) THEN
           WRITE(IU0,'(A,20A)') ' POSCAR group in vaspin.h5 found, type information: ',TYPE(1:T_INFO%NTYP)
         ELSE
           TYPE(1:T_INFO%NTYP+1)='  '
        ENDIF
        IF (IERR == 0 .AND. T_INFO%NTYP/=NUMB) THEN
           CALL vtutor%error("ERROR: the type information is not consistent with the number of types")
        ENDIF
      ENDIF
#endif
      T_INFO%NTYPP=T_INFO%NTYP
! how many ions do we have on file POSCAR ... ?
      T_INFO%NIONS=0
      DO NI=1,T_INFO%NTYP
         T_INFO%NIONS=T_INFO%NIONS+NITYP(NI)
      END DO

! there might be empty spheres scan for them

      T_INFO%NIONP=T_INFO%NIONS
      T_INFO%NTYPP=T_INFO%NTYP
#ifdef VASP_HDF5
      IF (POSCAR_FOUND) THEN
#endif
      READ(15,'(A1)',ERR=147,END=147) CHARAC
      IF ((CHARAC=='S').OR.(CHARAC=='s')) &
     &   READ(15,'(A1)',ERR=147,END=147) CHARAC
      DO NI=1,T_INFO%NIONS
         READ(15,'(A1)',ERR=147,END=147) CHARAC
      END DO

      READ(15,'(A1)',ERR=147,END=147) CHARAC
      IF ((CHARAC=='E').OR.(CHARAC=='e')) THEN
! this is also important for us ...
         READ(15,'(A)',ERR=147,END=147) INPLIN
! how many words/data items? --> number of empty sphere types!
         T_INFO%NTYPP=T_INFO%NTYPP+NITEMS(INPLIN,INPWRK,.TRUE.,'I')
         READ(INPLIN,*) (NITYP(NT),NT=T_INFO%NTYP+1,T_INFO%NTYPP)
         DO NT=T_INFO%NTYP+1,T_INFO%NTYPP
           T_INFO%NIONP=T_INFO%NIONP+NITYP(NT)
         ENDDO
         TYPE(T_INFO%NTYP+1:T_INFO%NTYPP)='  '
      ENDIF
! ... precise details later in the program ...
  147 REWIND 15
! set the require allocation parameters
#ifdef VASP_HDF5
      ELSE
        IERR = VH5_READ(IH5INPOSCARGROUP_ID,"number_empty_sphere_types", NITYP, DIMENSIONS=NUMB)
        IF (IERR==0) THEN
          T_INFO%NTYPP=T_INFO%NTYPP+NUMB
          DO NT=T_INFO%NTYP+1,T_INFO%NTYPP
            T_INFO%NIONP=T_INFO%NIONP+NITYP(NT)
          ENDDO
          TYPE(T_INFO%NTYP+1:T_INFO%NTYPP)='  '
        ENDIF
      ENDIF
#endif
      NIOND =T_INFO%NIONS
      NTYPD =T_INFO%NTYP
      NIONPD=T_INFO%NIONP
      NTYPPD=T_INFO%NTYPP
#ifdef VASP_HDF5
      ALLOCATE(T_INFO%NITYP(NTYPPD),T_INFO%TYPE(NTYPPD), T_INFO%SHA256(NTYPPD), T_INFO%TYPEF(NTYPPD))
#else
      ALLOCATE(T_INFO%NITYP(NTYPPD),T_INFO%TYPE(NTYPPD))
#endif
      T_INFO%NITYP(1:NTYPPD)=NITYP(1:NTYPPD)
      T_INFO%TYPE(1:NTYPPD) =TYPE (1:NTYPPD)

#ifdef VASP_HDF5
      T_INFO%SHA256(1:NTYPPD) = SHA256(1:NTYPPD)
      T_INFO%TYPEF(1:NTYPPD) = TYPEF(1:NTYPPD)

      IF (POSCAR_FOUND) THEN
#endif
      IF (IU0>=0) &
      WRITE(IU0,1) DIR_APP(1:DIR_LEN),NTYPPD,NIONPD
#ifdef VASP_HDF5
      ELSE
        IF (IU0>=0) &
        WRITE(IU0,2) DIR_APP(1:DIR_LEN),NTYPPD,NIONPD
      ENDIF
#endif

    1 FORMAT(' ',A,'POSCAR found : ',I2,' types and ',I7,' ions' )
    2 FORMAT(' ',A,'vaspin.h5 (POSCAR) found : ',I2,' types and ',I7,' ions' )
#ifdef VASP_HDF5
      IF (POSCAR_FOUND) THEN
#endif
      CLOSE(UNIT=15)
#ifdef VASP_HDF5
      ENDIF
#endif
      RETURN
 1000 CONTINUE
!
! all nodes report to unit 6 which have IU6 defined
! (guarantees that a sensible error message is allways written out)
!
      CALL vtutor%error("ERROR: the following files does not exist " // DIR_APP(1:DIR_LEN) // "POSCAR")
      END SUBROUTINE

!=======================================================================
!
!  Read UNIT=15: POSCAR Startjob and Continuation-job
!
!=======================================================================
      SUBROUTINE RD_POSCAR(LATT_CUR, T_INFO, DYN, &
     &           NIOND,NIONPD, NTYPD,NTYPPD, &
     &           IU0,IU6)
      USE prec
      USE lattice
      USE main_mpi
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
#ifdef VASP_HDF5
      use vhdf5_base
      use main_mpi, only: comm
#endif

      IMPLICIT NONE

      INTEGER NIOND,NIONPD,NTYPPD,NTYPD
      CHARACTER (255)  INPLIN,INPWRK
      INTEGER, EXTERNAL :: NITEMS
      TYPE (latt)::       LATT_CUR
      TYPE (type_info) :: T_INFO
      TYPE (dynamics)  :: DYN
      INTEGER IU0,IU6        ! io unit
! temporary
      CHARACTER(LEN=:), ALLOCATABLE :: SZNAM2
      CHARACTER (1)  CSEL
      CHARACTER (32) CSEL_LONG
      INTEGER I,NT,NI,NSCALE
      REAL(q) SCALEX,SCALEY,SCALEZ
      REAL(q) POTIMR

      INTEGER :: IBRAVS(3)

#ifdef VASP_HDF5
      INTEGER(HID_T) :: IH5OUTPOSCARGROUP_ID, IH5INPOSCARGROUP_ID
      INTEGER :: IH5ERR
      INTEGER :: NODE_ME, IONODE
      INTEGER :: IERR, NUMB
      LOGICAL :: LDUM
!      CHARACTER(LEN=:), ALLOCATABLE :: FORM
!      NODE_ME = COMM%NODE_ME
!      IONODE = COMM%IONODE

      IF (POSCAR_FOUND) THEN
#endif

      OPEN(UNIT=15,FILE=DIR_APP(1:DIR_LEN)//'POSCAR',STATUS='OLD')
#ifdef VASP_HDF5
      ELSE IF (HDF5_FOUND) THEN
        IERR = VH5_GROUP_OPEN(IH5ININPUTGROUP_ID, SUBGRP_POSCAR, IH5INPOSCARGROUP_ID)
      ENDIF
      IERR = VH5_GROUP_OPEN_OR_CREATE(IH5OUTINPUTGROUP_ID, SUBGRP_POSCAR, IH5OUTPOSCARGROUP_ID)
#endif
      IF (IU6>=0) WRITE(IU6,*)
!-----Basis vectors and scaling parameter ('lattice constant')
      T_INFO%SZNAM2=''
#ifdef VASP_HDF5
      IF (POSCAR_FOUND) THEN
#endif
      READ(15,'(A40)') T_INFO%SZNAM2
#ifdef VASP_HDF5
      ELSE IF (HDF5_FOUND) THEN
        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "system", SZNAM2)
        IF (ALLOCATED(SZNAM2)) T_INFO%SZNAM2 = SZNAM2
      ENDIF
      IERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, "system", T_INFO%SZNAM2, SKIP_PRESENT=.TRUE.)
#endif
IF (IU6>=0) WRITE(IU6,*)'POSCAR: ',TRIM(T_INFO%SZNAM2)
! one scaling parameter or one for x, y and z
#ifdef VASP_HDF5
      IF (POSCAR_FOUND) THEN
#endif
      READ(15,'(A)') INPLIN
! how many words/data items? --> number of ion types on file POSCAR!
      NSCALE=NITEMS(INPLIN,INPWRK,.TRUE.,'F')
      IF (NSCALE==1) THEN
         READ(INPLIN,*) LATT_CUR%SCALE
         SCALEX=1
         SCALEY=1
         SCALEZ=1
      ELSE IF (NSCALE==3) THEN
         LATT_CUR%SCALE=1
         READ(INPLIN,*) SCALEX,SCALEY,SCALEZ
      ELSE
         CALL vtutor%error("ERROR: there must be 1 or 3 items on line 2 of POSCAR")
      ENDIF

      DO I=1,3
        READ(15,*) LATT_CUR%A(1,I),LATT_CUR%A(2,I),LATT_CUR%A(3,I)
      ENDDO
#ifdef VASP_HDF5
      ELSE IF (HDF5_FOUND) THEN
        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "scale", LATT_CUR%SCALE)
        IF (IERR==0) THEN
          NSCALE=1
          SCALEX=1
          SCALEY=1
          SCALEZ=1
        ELSE
          NSCALE=3
          IERR = VH5_READ(IH5INPOSCARGROUP_ID, "scale_x", SCALEX)
          IERR = VH5_READ(IH5INPOSCARGROUP_ID, "scale_y", SCALEY)
          IERR = VH5_READ(IH5INPOSCARGROUP_ID, "scale_z", SCALEZ)
          LATT_CUR%SCALE=1
        ENDIF
        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "lattice_vectors", LATT_CUR%A)
      ENDIF
      IF (NSCALE==1) THEN
        IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, "scale", LATT_CUR%SCALE, SKIP_PRESENT=.TRUE.)
      ELSE IF (NSCALE==3) THEN
        IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, "scale_x", SCALEX, SKIP_PRESENT=.TRUE.)
        IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, "scale_y", SCALEY, SKIP_PRESENT=.TRUE.)
        IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, "scale_z", SCALEZ, SKIP_PRESENT=.TRUE.)
      ENDIF
      IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, "lattice_vectors", LATT_CUR%A, SKIP_PRESENT=.TRUE.)
#endif
      IF (LATT_CUR%SCALE<0._q) THEN
!----alternatively give a volume (=abs(scale)) and adjust the lengths of
!----the three lattice vectors to get the correct desired volume ... :
         CALL LATTIC(LATT_CUR)
         LATT_CUR%SCALE=(ABS(LATT_CUR%SCALE)  &
     &                 / ABS(LATT_CUR%OMEGA))**(1._q/3._q)
      ENDIF

      LATT_CUR%A(1,:) =LATT_CUR%A(1,:)*SCALEX*LATT_CUR%SCALE
      LATT_CUR%A(2,:) =LATT_CUR%A(2,:)*SCALEY*LATT_CUR%SCALE
      LATT_CUR%A(3,:) =LATT_CUR%A(3,:)*SCALEZ*LATT_CUR%SCALE

      CALL LATTIC(LATT_CUR)

      IF (LATT_CUR%OMEGA<0) &
        CALL vtutor%error("ERROR: the triple product of the basis vectors is negative exchange two basis vectors")

! sanity check (this should safe us trouble in IBZKPT lateron)
      IBRAVS = LATTCHK(LATT_CUR)

      T_INFO%NIOND =NIOND
      T_INFO%NIONPD=NIONPD
      T_INFO%NTYPD =NTYPD
      T_INFO%NTYPPD=NTYPPD
      ALLOCATE(T_INFO%LSFOR(3,NIOND),T_INFO%ITYP(NIOND))

      T_INFO%LSFOR=.TRUE.

!-----number of atoms per type
#ifdef VASP_HDF5
      IF (POSCAR_FOUND) THEN
#endif
      READ(15,'(A)') INPLIN
      READ(INPLIN,*) CSEL
      IF (.NOT.(CSEL>='0' .AND. CSEL<='9')) THEN
         ! construct format specifier
        !ALLOCATE(CHARACTER(2 + 4*(T_INFO%NTYP-1) + 1) :: FORM)
        !FORM(1:2)="("
        !DO NI=1,T_INFO%NTYP-1
        !  FORM(4*(NI-1)+2:4*NI+1) = "A,X,"
        !ENDDO
        !FORM(4*(T_INFO%NTYP-1)+2:4*(T_INFO%NTYP-1)+3) = "A)"
#ifdef VASP_HDF5
         IERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'ion_types', T_INFO%TYPEF(1:T_INFO%NTYP), SKIP_PRESENT=.TRUE.)
         IERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'ion_sha256', T_INFO%SHA256(1:T_INFO%NTYP), SKIP_PRESENT=.TRUE.)
#else
         READ(INPLIN,*) (T_INFO%TYPE(NI),NI=1,T_INFO%NTYP)
#endif
         READ(15,'(A)') INPLIN
         IF (T_INFO%NTYP/=NITEMS(INPLIN,INPWRK,.TRUE.,'I')) THEN
            CALL vtutor%error("ERROR: the type information is not consistent with the number of types")
         ENDIF
      ENDIF

      READ(INPLIN,*) (T_INFO%NITYP(NT),NT=1,T_INFO%NTYP)
#ifdef VASP_HDF5
      ELSE IF (HDF5_FOUND) THEN
        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "number_ion_types", T_INFO%NITYP, DIMENSIONS=T_INFO%NTYP)
        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "ion_types", T_INFO%TYPEF)
        DO I=1,T_INFO%NTYP
           T_INFO%TYPE(I) = T_INFO%TYPEF(I)
           IF (T_INFO%TYPE(I)(2:2) == '_') THEN
              T_INFO%TYPE(I)(2:2) = ' '
           ENDIF
        ENDDO
        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "ion_sha256", T_INFO%SHA256)
        IF (IERR ==0) THEN
         IERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'ion_types', T_INFO%TYPEF(1:T_INFO%NTYP), SKIP_PRESENT=.TRUE.)
         IERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'ion_sha256', T_INFO%SHA256(1:T_INFO%NTYP), SKIP_PRESENT=.TRUE.)
        ENDIF
      ENDIF
      IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'number_ion_types', T_INFO%NITYP(1:T_INFO%NTYP), SKIP_PRESENT=.TRUE.)
#endif

!---- Set up the table from which we get type of each ion
      NI=1
      DO NT=1,T_INFO%NTYP
      DO NI=NI,T_INFO%NITYP(NT)+NI-1
        T_INFO%ITYP(NI)=NT
      ENDDO
      ENDDO
!
!   positions
!
      T_INFO%NIONS=0
      DO NT=1,T_INFO%NTYP
      T_INFO%NIONS= T_INFO%NIONS+ T_INFO%NITYP(NT)
      ENDDO

      T_INFO%NIONP=T_INFO%NIONS

      IF (T_INFO%NIONS>NIOND) THEN
        CALL vtutor%error("ERROR: MAIN: increase NIOND " // str(T_INFO%NIONS))
      ENDIF

#ifdef VASP_HDF5
      IF (POSCAR_FOUND) THEN
#endif
      READ(15,*) CSEL_LONG
      ! remove whitespaces at the edges
      CALL STRIP( CSEL_LONG, I, 'A' )
      ! only first letter matters in the following 
      CSEL = CSEL_LONG(1:1)
      T_INFO%LSDYN=((CSEL=='S').OR.(CSEL=='s'))
      IF (T_INFO%LSDYN) READ(15,*) CSEL_LONG
      ! remove whitespaces at the edges
      CALL STRIP( CSEL_LONG, I, 'A' )
      ! only first letter matters in the following 
      CSEL = CSEL_LONG(1:1)
      IF (CSEL=='K'.OR.CSEL=='k'.OR. &
     &    CSEL=='C'.OR.CSEL=='c') THEN
        CSEL='K'
        IF (IU6>=0) &
        WRITE(IU6,*)' positions in cartesian coordinates'
        T_INFO%LDIRCO=.FALSE.
      ELSE
        IF (IU6>=0) &
        WRITE(IU6,*)' positions in direct lattice'
        T_INFO%LDIRCO=.TRUE.
      ENDIF
#ifdef VASP_HDF5
      ELSE IF (HDF5_FOUND) THEN
        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "selective_dynamics", T_INFO%LSDYN)
        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "direct_coordinates", T_INFO%LDIRCO)
        IF (T_INFO%LDIRCO) THEN
          CSEL='D'
          IF (IU6>=0) &
          WRITE(IU6,*)' positions in direct coordinates'
        ELSE
          CSEL='K'
          IF (IU6>=0) &
          WRITE(IU6,*)' positions in cartesian lattice'
        ENDIF
      ENDIF
      IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'selective_dynamics', T_INFO%LSDYN, SKIP_PRESENT=.TRUE.)
      IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'direct_coordinates', T_INFO%LDIRCO, SKIP_PRESENT=.TRUE.)
#endif

      ALLOCATE(DYN%POSION(3,NIONPD),DYN%POSIOC(3,NIONPD), &
     &         DYN%D2C(3,NIOND), &
     &         DYN%VEL(3,NIOND),DYN%D2(3,NIOND),DYN%D3(3,NIOND))

! alias T_INFO%POSION
      T_INFO%POSION => DYN%POSION
      DYN%POSION=0
      DYN%VEL   =0
      DYN%D2    =0
      DYN%D2C   =0
      DYN%D3    =0
!tb start
      LATT_CUR%AVEL=0
!tb end

#ifdef VASP_HDF5
      IF (POSCAR_FOUND) THEN
#endif
      DO NI=1,T_INFO%NIONS
      IF (T_INFO%LSDYN) THEN
      READ(15,*,ERR=400,END=400) DYN%POSION(1,NI),DYN%POSION(2,NI),DYN%POSION(3,NI), &
     &      T_INFO%LSFOR(1,NI),T_INFO%LSFOR(2,NI),T_INFO%LSFOR(3,NI)
      ELSE
      READ(15,*,ERR=400,END=400) DYN%POSION(1,NI),DYN%POSION(2,NI),DYN%POSION(3,NI)
      ENDIF
      ENDDO
#ifdef VASP_HDF5
      ELSE IF (HDF5_FOUND) THEN
        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "position_ions", DYN%POSION(:,1:T_INFO%NIONS))
        IF (T_INFO%LSDYN) THEN
          IERR = VH5_READ(IH5INPOSCARGROUP_ID, "selective_dynamics_ions", T_INFO%LSFOR)
        ENDIF
      ENDIF
      IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'position_ions', DYN%POSION(:,1:T_INFO%NIONS), SKIP_PRESENT=.TRUE.)
      IF (T_INFO%LSDYN) THEN
         IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'selective_dynamics_ions', T_INFO%LSFOR, SKIP_PRESENT=.TRUE.)
      ENDIF
#endif
      IF (CSEL=='K') THEN
        DYN%POSION(1,:)=LATT_CUR%SCALE*DYN%POSION(1,:)*SCALEX
        DYN%POSION(2,:)=LATT_CUR%SCALE*DYN%POSION(2,:)*SCALEY
        DYN%POSION(3,:)=LATT_CUR%SCALE*DYN%POSION(3,:)*SCALEZ
        
        CALL KARDIR(T_INFO%NIONS,DYN%POSION,LATT_CUR%B)
      ENDIF
      CALL TOPRIM(T_INFO%NIONS,DYN%POSION)
      DYN%POSIOC=DYN%POSION

      DYN%INIT=0
      DYN%SNOSE(1)=1
      LATT_CUR%INITlatv=-1

!
!   empty spheres
!
#ifdef VASP_HDF5
      IF (POSCAR_FOUND) THEN
#endif
      READ(15,'(A1)',ERR=424,END=410) CSEL
 ! read in new components and velocities of lattice vector components
      IF ((CSEL=='L').OR.(CSEL=='l')) THEN
        READ(15,*,ERR=4001,END=4001)  LATT_CUR%INITlatv
#ifdef VASP_HDF5
        IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'lattice_velocities_initialized', LATT_CUR%INITLATV, SKIP_PRESENT=.TRUE.)
#endif
        IF (  LATT_CUR%INITlatv>=0) THEN
          DO NI=1,3
            READ(15,'(3E16.8)',ERR=4001,END=4001) LATT_CUR%AVEL(1,NI),LATT_CUR%AVEL(2,NI),LATT_CUR%AVEL(3,NI)
          ENDDO
#ifdef VASP_HDF5
        IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'lattice_velocities', LATT_CUR%AVEL, SKIP_PRESENT=.TRUE.)
#endif
          LATT_CUR%AVEL=LATT_CUR%AVEL*DYN%POTIM
        ENDIF
        
        IF (  LATT_CUR%INITlatv>0) THEN
          DO NI=1,3
             READ(15,'(3E16.8)',ERR=4001,END=4001) LATT_CUR%A(1,NI),LATT_CUR%A(2,NI),LATT_CUR%A(3,NI)
          ENDDO
#ifdef VASP_HDF5
          IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'lattice_vectors_initialized', LATT_CUR%A, SKIP_PRESENT=.TRUE.)
#endif

        ENDIF
        CALL LATTIC(LATT_CUR)
        READ(15,'(A1)',ERR=424,END=410) CSEL
      ENDIF
#ifdef VASP_HDF5
      ELSE IF (HDF5_FOUND) THEN
        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "lattice_velocities_initialized", LATT_CUR%INITlatv)
        IF (IERR == 3) GOTO 424
        IF (  LATT_CUR%INITlatv>=0) THEN
          IERR = VH5_READ(IH5INPOSCARGROUP_ID, "lattice_velocities", LATT_CUR%AVEL)
          IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'lattice_velocities', LATT_CUR%AVEL, SKIP_PRESENT=.TRUE.)
          LATT_CUR%AVEL=LATT_CUR%AVEL*DYN%POTIM
        ENDIF
        IF (  LATT_CUR%INITlatv>0) THEN
          IERR = VH5_READ(IH5INPOSCARGROUP_ID, "lattice_vectors_initialized", LATT_CUR%A)
          IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'lattice_vectors_initialized', LATT_CUR%A, SKIP_PRESENT=.TRUE.)
        ENDIF
        CALL LATTIC(LATT_CUR)
      ENDIF
#endif

#ifdef VASP_HDF5
      IF (POSCAR_FOUND) THEN
#endif
  424 IF ((CSEL=='E').OR.(CSEL=='e')) THEN
        IF (T_INFO%NTYPP>NTYPPD) THEN
         CALL vtutor%error("ERROR: MAIN: increase NEMPTY " // str(T_INFO%NTYPP-T_INFO%NTYP))
        ENDIF
        READ(15,*,ERR=410,END=410) (T_INFO%NITYP(NT),NT=T_INFO%NTYP+1,T_INFO%NTYPP)
#ifdef VASP_HDF5
        IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'number_empty_sphere_types', T_INFO%NITYP(T_INFO%NTYP+1:T_INFO%NTYPP), SKIP_PRESENT=.TRUE.)
#endif
        DO NT=T_INFO%NTYP+1,T_INFO%NTYPP
          T_INFO%NIONP=T_INFO%NIONP+T_INFO%NITYP(NT)
        ENDDO
        IF (T_INFO%NIONP>NIONPD) THEN
         CALL vtutor%error("ERROR: MAIN: increase NEMPTY " // str(T_INFO%NIONP-T_INFO%NIONS))
        ENDIF
        T_INFO%NIONP=T_INFO%NIONP

        DO NI=T_INFO%NIONS+1,T_INFO%NIONP
         READ(15,*,ERR=410,END=410) &
     &      DYN%POSION(1,NI),DYN%POSION(2,NI),DYN%POSION(3,NI)
        ENDDO
#ifdef VASP_HDF5
        IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'position_empty_spheres', DYN%POSION(:, T_INFO%NIONS+1:T_INFO%NIONP), SKIP_PRESENT=.TRUE.)
#endif
        IF (.NOT.T_INFO%LDIRCO) THEN
          DO NI=T_INFO%NIONS+1,T_INFO%NIONP
            DYN%POSION(1,NI)=LATT_CUR%SCALE*DYN%POSION(1,NI)*SCALEX
            DYN%POSION(2,NI)=LATT_CUR%SCALE*DYN%POSION(2,NI)*SCALEY
            DYN%POSION(3,NI)=LATT_CUR%SCALE*DYN%POSION(3,NI)*SCALEZ
            CALL KARDIR(1,DYN%POSION(1:3,NI),LATT_CUR%B)
          ENDDO
        ENDIF
        DYN%POSIOC=DYN%POSION
        READ(15,'(A1)',ERR=425,END=410) CSEL
      ENDIF
#ifdef VASP_HDF5
      ELSE IF (HDF5_FOUND) THEN
        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "number_empty_sphere_types",  T_INFO%NITYP(T_INFO%NTYP+1:T_INFO%NTYPP))
        IF (IERR == 0) THEN
          IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'number_empty_sphere_types', T_INFO%NITYP(T_INFO%NTYP+1:T_INFO%NTYPP), SKIP_PRESENT=.TRUE.)
        ELSE IF (IERR == 3) THEN
          GOTO 410
        ENDIF
        DO NT=T_INFO%NTYP+1,T_INFO%NTYPP
          T_INFO%NIONP=T_INFO%NIONP+T_INFO%NITYP(NT)
        ENDDO
        IF (T_INFO%NIONP>NIONPD) THEN
         CALL vtutor%error("ERROR: MAIN: increase NEMPTY " // str(T_INFO%NIONP-T_INFO%NIONS))
        ENDIF
        T_INFO%NIONP=T_INFO%NIONP

        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "position_empty_spheres", DYN%POSION(:, T_INFO%NIONS+1:T_INFO%NIONP))
        IF (IERR == 0) THEN
          IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'position_empty_spheres', DYN%POSION(:, T_INFO%NIONS+1:T_INFO%NIONP), SKIP_PRESENT=.TRUE.)
        ELSE IF (IERR == 3) THEN
          GOTO 410
        ENDIF
        IF (.NOT.T_INFO%LDIRCO) THEN
          DO NI=T_INFO%NIONS+1,T_INFO%NIONP
            DYN%POSION(1,NI)=LATT_CUR%SCALE*DYN%POSION(1,NI)*SCALEX
            DYN%POSION(2,NI)=LATT_CUR%SCALE*DYN%POSION(2,NI)*SCALEY
            DYN%POSION(3,NI)=LATT_CUR%SCALE*DYN%POSION(3,NI)*SCALEZ
            CALL KARDIR(1,DYN%POSION(1:3,NI),LATT_CUR%B)
          ENDDO
        ENDIF
      DYN%POSIOC=DYN%POSION
      ENDIF
#endif


      DYN%POSIOC=DYN%POSION
#ifdef VASP_HDF5
      IF (POSCAR_FOUND) THEN
#endif
  425 IF (CSEL=='K'.OR.CSEL=='k'.OR.CSEL==' ' &
     &    .OR.CSEL=='C'.OR.CSEL=='c') THEN
        CSEL='K'
        IF (IU6>=0) &
        WRITE(IU6,*)' velocities in cartesian coordinates'
#ifdef VASP_HDF5
        IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'direct_coordinates_velocities', .FALSE., SKIP_PRESENT=.TRUE.)
#endif
      ELSE
        IF (IU6>=0) &
        WRITE(IU6,*)' velocities in direct lattice'
#ifdef VASP_HDF5
        IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'direct_coordinates_velocities', .TRUE., SKIP_PRESENT=.TRUE.)
#endif
      ENDIF
#ifdef VASP_HDF5
      ELSE IF (HDF5_FOUND) THEN
        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "direct_coordinates_velocities", LDUM)
        IF (IERR == 0) THEN
          IF (LDUM .EQV. .TRUE.) THEN
            CSEL="D"
          ELSE
            CSEL="K"
          ENDIF
        ENDIF
      ENDIF
#endif
!
!-----if we have velocities, read them in and transform from
!     cartesian coordinates to direct lattice
#ifdef VASP_HDF5
      IF (POSCAR_FOUND) THEN
#endif
      DO NI=1,T_INFO%NIONS
        READ(15,*,ERR=410,END=410)  &
     &            DYN%VEL(1,NI),DYN%VEL(2,NI),DYN%VEL(3,NI)
      ENDDO
#ifdef VASP_HDF5
      ELSE IF (HDF5_FOUND) THEN
        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "ion_velocities", DYN%VEL)
        IF (IERR == 3) THEN
          GOTO 410
        ENDIF
      ENDIF
      IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'ion_velocities', DYN%VEL, SKIP_PRESENT=.TRUE.)
#endif
      DO NI=1,T_INFO%NIONS
        !tb start
        IF (CSEL=='K' .AND. DYN%IBRION/=44 .AND. DYN%IBRION/=40) THEN
        !IF (CSEL=='K') THEN
        !tb end
        CALL  KARDIR(1,DYN%VEL(1:3,NI),LATT_CUR%B)
        DYN%VEL(1:3,NI)=DYN%VEL(1:3,NI)*DYN%POTIM
        ENDIF
      ENDDO
      
!
!-----try to read in predictor Coordinates
!
#ifdef VASP_HDF5
      IF (POSCAR_FOUND) THEN
#endif
      READ(15,*,ERR=430,END=430)
      READ(15,*,ERR=430,END=430) DYN%INIT
#ifdef VASP_HDF5
      ELSE IF (HDF5_FOUND) THEN
        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "dyn_init", DYN%INIT)
        IF (IERR == 3) THEN
          GOTO 430
        ENDIF
      ENDIF
      IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'dyn_init', DYN%INIT, SKIP_PRESENT=.TRUE.)
#endif
 
      !c initialization of ionic and latice DOFs must be consistent
      !c kill the job otherwise
      !c the folloving combinations of INITlatv and INIT make sense:
      !c     -1 -1 - just structure is read in
      !c     -1  0 - structure and ionic velocities
      !c      0 -1 - structure and lattice velocities
      !c      0  0 - structure and all velocities
      !c      1  1 - continuation of previous MD
      IF (DYN%IBRION==0 .AND. (DYN%ISIF>=3)) THEN
        IF ((DYN%INIT/=1 .AND.  LATT_CUR%INITlatv==1) .OR. (DYN%INIT==1 .AND.  LATT_CUR%INITlatv/=1)) GOTO 4002
      ENDIF

!-----if INIT is there and it is 1 we have predictor-coordinates on the
!-----file so we can start with them
      IF (DYN%INIT==0) GOTO 430
#ifdef VASP_HDF5
      IF (POSCAR_FOUND) THEN
#endif
      READ(15,*) POTIMR
#ifdef VASP_HDF5
      ELSE IF (HDF5_FOUND) THEN
        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "potim", POTIMR)
      ENDIF
      IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'potim', POTIMR, SKIP_PRESENT=.TRUE.)
#endif
      IF (POTIMR/=DYN%POTIM) THEN
        IF (IU6>=0) THEN
           WRITE(IU6,*)
           WRITE(IU6,*)' There are predictor-coordinates on the file.'
           WRITE(IU6,*)' we can''t use them due to change of POTIM!'
        ENDIF
        GOTO 430
      ENDIF

#ifdef VASP_HDF5
      IF (POSCAR_FOUND) THEN
#endif
!-----Read in Nose-Parameter
      READ(15,*) DYN%SNOSE
!-----Read in predictor-coordinates (always in direct lattice)
      READ(15,*) (DYN%POSION(1,NI),DYN%POSION(2,NI),DYN%POSION(3,NI),NI=1,T_INFO%NIONS)
      READ(15,*) (DYN%D2(1,NI),DYN%D2(2,NI),DYN%D2(3,NI),NI=1,T_INFO%NIONS)
      READ(15,*) (DYN%D3(1,NI),DYN%D3(2,NI),DYN%D3(3,NI),NI=1,T_INFO%NIONS)
#ifdef VASP_HDF5
      ELSE IF (HDF5_FOUND) THEN
        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "potim", POTIMR)
        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "predictor_coordinates", DYN%POSION(:,1:T_INFO%NIONS))
        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "predictor_coordinates_2", DYN%D2)
        IERR = VH5_READ(IH5INPOSCARGROUP_ID, "predictor_coordinates_3", DYN%D3)
      ENDIF
      IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'nose_thermostat', DYN%SNOSE, SKIP_PRESENT=.TRUE.)
      IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'predictor_coordinates', DYN%POSION(:,1:T_INFO%NIONS), SKIP_PRESENT=.TRUE.)
      IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'predictor_coordinates_2', DYN%D2, SKIP_PRESENT=.TRUE.)
      IH5ERR = VH5_WRITE(IH5OUTPOSCARGROUP_ID, 'predictor_coordinates_3', DYN%D3, SKIP_PRESENT=.TRUE.)
#endif

      IF (IU6>=0) THEN
         WRITE(IU6,*)
         WRITE(IU6,*)' Using predictor-coordinates on the file'
      ENDIF

#ifdef VASP_HDF5
      IF (POSCAR_FOUND) THEN
#endif
      CLOSE(UNIT=15)
#ifdef VASP_HDF5
      ENDIF
#endif
      RETURN
!-----------------------------------------------------------------------
!  Reading Inputfile 15 finished
!  if you end up at 430  INIT is set to 0,
!    INIT is used in the call to STEP  (predictors are not initialised)
!    in that way we tell STEP that it must initialize everything for us
!----------------------------------------------------------------------
  400 CONTINUE
      CALL vtutor%error("No initial positions read in")

  4001 CONTINUE
      CALL vtutor%error("No initial lattice velocities or vectors read in")

  4002 CONTINUE
      CALL vtutor%error("Inconsistent initialization of ion. and lat. DOFs")

  410 CONTINUE
      DYN%INIT=-1
      IF (IU6>=0) &
      WRITE(IU6,*)' No initial velocities read in'

      CLOSE(UNIT=15)
      !c initialized lattice DOFs and unitialized atomic DOFs makes no sense
      IF (DYN%IBRION==0 .AND. (DYN%ISIF==3 .OR.  DYN%ISIF==8 .OR. DYN%ISIF==9 .OR. DYN%ISIF==10)) THEN
        IF ( LATT_CUR%INITlatv==1 ) THEN
          CALL vtutor%error("Inconsistent initialization of ion. and lat. DOFs")
        ENDIF
      ENDIF

      RETURN

  430 DYN%INIT=0
      CLOSE(UNIT=15)
      !c using next lattice vectors and previous atomic positions makes no sense
      IF (DYN%IBRION==0 .AND. (DYN%ISIF==3 .OR.  DYN%ISIF==8 .OR. DYN%ISIF==9 .OR. DYN%ISIF==10)) THEN
        IF ( LATT_CUR%INITlatv==1 ) THEN
          CALL vtutor%error("Inconsistent initialization of ion. and lat. DOFs")
        ENDIF
      ENDIF
      RETURN

      END SUBROUTINE



!=======================================================================
!  Read a POSCAR file FILE_NAME
!  This is a thin wrapper around RD_POSCAR_ALT_(HEAD|FULL) which are
!  functional copies of RD_POSCAR_HEAD and RD_POSCAR, respetively.
!  At some point, these should be merged with the original routines.
!=======================================================================
      SUBROUTINE RD_POSCAR_ALT(FILE_NAME, LATT_ALT, T_INFO, DYN, IU0, IU6)
         USE prec
         USE lattice
         IMPLICIT NONE
   
         CHARACTER(*),     INTENT(IN)  :: FILE_NAME
         TYPE(latt),       INTENT(OUT) :: LATT_ALT
         TYPE(type_info),  INTENT(OUT) :: T_INFO
         TYPE(dynamics),   INTENT(OUT) :: DYN
         INTEGER,          INTENT(IN)  :: IU0
         INTEGER,          INTENT(IN)  :: IU6
   
         INTEGER NIOND, NIONPD, NTYPPD, NTYPD
   
         CALL RD_POSCAR_ALT_HEAD(FILE_NAME, LATT_ALT, T_INFO, NIOND, NIONPD, NTYPD, NTYPPD, IU0, IU6)
         CALL RD_POSCAR_ALT_FULL(FILE_NAME, LATT_ALT, T_INFO, DYN, NIOND, NIONPD, NTYPD, NTYPPD, IU0, IU6)
      END SUBROUTINE
   
!=======================================================================
!  Read a POSCAR file FILE_NAME and scan for total number of ions
!  and number of types
!  only T_INFO%NITYP is allocated at this point
!  all other arrays are allocated in RD_POSCAR_ALT_FULL
!=======================================================================
      SUBROUTINE RD_POSCAR_ALT_HEAD(FILE_NAME, LATT_ALT, T_INFO, &
                 NIOND, NIONPD, NTYPD, NTYPPD, IU0, IU6)
         USE prec
         USE lattice
         USE main_mpi
         USE tutor,  ONLY: vtutor
         USE string, ONLY: str
         IMPLICIT NONE
   
         CHARACTER(*),     INTENT(IN)  :: FILE_NAME
         TYPE (latt),      INTENT(OUT) :: LATT_ALT
         TYPE (type_info), INTENT(OUT) :: T_INFO
         INTEGER,          INTENT(OUT) :: NIOND
         INTEGER,          INTENT(OUT) :: NIONPD
         INTEGER,          INTENT(OUT) :: NTYPPD
         INTEGER,          INTENT(OUT) :: NTYPD
         INTEGER,          INTENT(IN)  :: IU0
         INTEGER,          INTENT(IN)  :: IU6
   
         INTEGER NITYP(10000)          ! hard limit 10000 ions :->
         CHARACTER (LEN=2) TYPE(10000) ! type information
   ! temporary varibales
         CHARACTER (1)    CHARAC
         CHARACTER (255)  INPLIN,INPWRK
         INTEGER        NI,I,NT,NSCALE
         REAL(q)        SCALEX,SCALEY,SCALEZ
         INTEGER, EXTERNAL :: NITEMS
   
   ! Now extract from file POSCAR how many ion types we have ...
         OPEN(UNIT=15,FILE=DIR_APP(1:DIR_LEN)//FILE_NAME,STATUS='OLD',ERR=1000)
   
         READ(15,'(A1)',ERR=147,END=147) CHARAC
   
   ! one scaling parameter or one for x, y and z
         READ(15,'(A)',ERR=147,END=147) INPLIN
   ! how many words/data items? --> number of ion types on file POSCAR!
         NSCALE=NITEMS(INPLIN,INPWRK,.TRUE.,'F')
         IF (NSCALE==1) THEN
            READ(INPLIN,*) LATT_ALT%SCALE
            SCALEX=1
            SCALEY=1
            SCALEZ=1
         ELSE IF (NSCALE==3) THEN
            LATT_ALT%SCALE=1
            READ(INPLIN,*) SCALEX,SCALEY,SCALEZ
         ELSE
            CALL vtutor%error("ERROR: there must be 1 or 3 items on line 2 of "//FILE_NAME)
         ENDIF
   
         DO I=1,3
           READ(15,*,ERR=147,END=147) LATT_ALT%A(1,I),LATT_ALT%A(2,I),LATT_ALT%A(3,I)
         ENDDO
   
         IF (LATT_ALT%SCALE<0._q) THEN
   !----alternatively give a volume (=abs(scale)) and adjust the lengths of
   !----the three lattice vectors to get the correct desired volume ... :
            CALL LATTIC(LATT_ALT)
            LATT_ALT%SCALE=(ABS(LATT_ALT%SCALE) &
        &                 / ABS(LATT_ALT%OMEGA))**(1._q/3._q)
         ENDIF
         
         LATT_ALT%A(1,:) =LATT_ALT%A(1,:)*SCALEX*LATT_ALT%SCALE
         LATT_ALT%A(2,:) =LATT_ALT%A(2,:)*SCALEY*LATT_ALT%SCALE
         LATT_ALT%A(3,:) =LATT_ALT%A(3,:)*SCALEZ*LATT_ALT%SCALE
            
         CALL LATTIC(LATT_ALT)
   
         IF (LATT_ALT%OMEGA<0) THEN
           CALL vtutor%error("ERROR: the triple product of the basis vectors is negative exchange two &
              &basis vectors")
         ENDIF
   
   ! 6th line, contains either the number of ions or their type
         READ(15,'(A)',ERR=147,END=147) INPLIN
   ! how many words/data items? --> number of ion types on file POSCAR!
         READ(INPLIN,*,ERR=147,END=147) CHARAC
         IF (.NOT.(CHARAC>='0' .AND. CHARAC<='9')) THEN
            T_INFO%NTYP=NITEMS(INPLIN,INPWRK,.TRUE.,'A')
            READ(INPLIN,*,ERR=147,END=147) (TYPE(NI),NI=1,T_INFO%NTYP)
            IF (IU0>=0) THEN
               WRITE(IU0,'(A,20A3)') FILE_NAME//' found type information on POSCAR ',TYPE(1:T_INFO%NTYP)
            ENDIF
   
            READ(15,'(A)',ERR=147,END=147) INPLIN
            IF (T_INFO%NTYP/=NITEMS(INPLIN,INPWRK,.TRUE.,'I')) THEN
               CALL vtutor%error("ERROR: the type information is not consistent with the number of types")
            ENDIF
         ELSE
            T_INFO%NTYP=NITEMS(INPLIN,INPWRK,.TRUE.,'I')
            TYPE(1:T_INFO%NTYP+1)='  '
         ENDIF
   
         T_INFO%NTYPP=T_INFO%NTYP
   ! let me know how many ions
         READ(INPLIN,*,ERR=147,END=147) (NITYP(NI),NI=1,T_INFO%NTYP)
   ! how many ions do we have on file POSCAR ... ?
         T_INFO%NIONS=0
         DO NI=1,T_INFO%NTYP
            T_INFO%NIONS=T_INFO%NIONS+NITYP(NI)
         END DO
   
   ! there might be empty spheres scan for them
   
         T_INFO%NIONP=T_INFO%NIONS
         T_INFO%NTYPP=T_INFO%NTYP
   
         READ(15,'(A1)',ERR=147,END=147) CHARAC
         IF ((CHARAC=='S').OR.(CHARAC=='s')) &
        &   READ(15,'(A1)',ERR=147,END=147) CHARAC
         DO NI=1,T_INFO%NIONS
            READ(15,'(A1)',ERR=147,END=147) CHARAC
         END DO
   
         READ(15,'(A1)',ERR=147,END=147) CHARAC
         IF ((CHARAC=='E').OR.(CHARAC=='e')) THEN
   ! this is also important for us ...
            READ(15,'(A)',ERR=147,END=147) INPLIN
   ! how many words/data items? --> number of empty sphere types!
            T_INFO%NTYPP=T_INFO%NTYPP+NITEMS(INPLIN,INPWRK,.TRUE.,'I')
            READ(INPLIN,*) (NITYP(NT),NT=T_INFO%NTYP+1,T_INFO%NTYPP)
            DO NT=T_INFO%NTYP+1,T_INFO%NTYPP
              T_INFO%NIONP=T_INFO%NIONP+NITYP(NT)
            ENDDO
            TYPE(T_INFO%NTYP+1:T_INFO%NTYPP)='  '
         ENDIF
   ! ... precise details later in the program ...
     147 REWIND 15
   ! set the require allocation parameters
   
         NIOND =T_INFO%NIONS
         NTYPD =T_INFO%NTYP
         NIONPD=T_INFO%NIONP
         NTYPPD=T_INFO%NTYPP
   
         ALLOCATE(T_INFO%NITYP(NTYPPD),T_INFO%TYPE(NTYPPD))
   
         T_INFO%NITYP(1:NTYPPD)=NITYP(1:NTYPPD)
         T_INFO%TYPE(1:NTYPPD) =TYPE (1:NTYPPD)
   
         IF (IU0>=0) &
         WRITE(IU0,1) DIR_APP(1:DIR_LEN)//FILE_NAME,NTYPPD,NIONPD
   
       1 FORMAT(' ',A,' found : ',I2,' types and ',I4,' ions' )
   
         CLOSE(UNIT=15)
         RETURN
    1000 CONTINUE
   !
   ! all nodes report to unit 6 which have IU6 defined
   ! (guarantees that a sensible error message is allways written out)
   !
         CALL vtutor%error("ERROR: the following file does not exist " // DIR_APP(1:DIR_LEN) // FILE_NAME)
      END SUBROUTINE
   
   !=======================================================================
   !
   !  Read full POSCAR file FILE_NAME
   !
   !=======================================================================
      SUBROUTINE RD_POSCAR_ALT_FULL(FILE_NAME, LATT_ALT, T_INFO, DYN, &
                 NIOND, NIONPD, NTYPD, NTYPPD, IU0, IU6)
         USE prec
         USE lattice
         USE main_mpi
         USE tutor,  ONLY: vtutor
         USE string, ONLY: str
         IMPLICIT NONE
   
         CHARACTER(*),     INTENT(IN)  :: FILE_NAME
         TYPE (latt),      INTENT(OUT) :: LATT_ALT
         TYPE (type_info), INTENT(OUT) :: T_INFO
         TYPE (dynamics),  INTENT(OUT) :: DYN
         INTEGER,          INTENT(IN)  :: NIOND
         INTEGER,          INTENT(IN)  :: NIONPD
         INTEGER,          INTENT(IN)  :: NTYPPD
         INTEGER,          INTENT(IN)  :: NTYPD
         INTEGER,          INTENT(IN)  :: IU0
         INTEGER,          INTENT(IN)  :: IU6
         
         CHARACTER (255)  INPLIN,INPWRK
         INTEGER, EXTERNAL :: NITEMS
   ! temporary
         CHARACTER (1)  CSEL
         INTEGER I,NT,NI,NSCALE
         REAL(q) SCALEX,SCALEY,SCALEZ
         REAL(q) POTIMR
   
         OPEN(UNIT=15,FILE=DIR_APP(1:DIR_LEN)//FILE_NAME,STATUS='OLD')
   
         IF (IU6>=0) WRITE(IU6,*)
   !-----Basis vectors and scaling parameter ('lattice constant')
         READ(15,'(A40)') T_INFO%SZNAM2
         IF (IU6>=0) WRITE(IU6,*) FILE_NAME//': ',T_INFO%SZNAM2
   
   ! one scaling parameter or one for x, y and z
         READ(15,'(A)') INPLIN
   ! how many words/data items? --> number of ion types on file POSCAR!
         NSCALE=NITEMS(INPLIN,INPWRK,.TRUE.,'F')
         IF (NSCALE==1) THEN
            READ(INPLIN,*) LATT_ALT%SCALE
            SCALEX=1
            SCALEY=1
            SCALEZ=1
         ELSE IF (NSCALE==3) THEN
            LATT_ALT%SCALE=1
            READ(INPLIN,*) SCALEX,SCALEY,SCALEZ
         ELSE
            CALL vtutor%error("ERROR: there must be 1 or 3 items on line 2 of "//FILE_NAME)
         ENDIF
   
         DO I=1,3
           READ(15,*) LATT_ALT%A(1,I),LATT_ALT%A(2,I),LATT_ALT%A(3,I)
         ENDDO
   
         IF (LATT_ALT%SCALE<0._q) THEN
   !----alternatively give a volume (=abs(scale)) and adjust the lengths of
   !----the three lattice vectors to get the correct desired volume ... :
            CALL LATTIC(LATT_ALT)
            LATT_ALT%SCALE=(ABS(LATT_ALT%SCALE)  &
        &                 / ABS(LATT_ALT%OMEGA))**(1._q/3._q)
         ENDIF
   
         LATT_ALT%A(1,:) =LATT_ALT%A(1,:)*SCALEX*LATT_ALT%SCALE
         LATT_ALT%A(2,:) =LATT_ALT%A(2,:)*SCALEY*LATT_ALT%SCALE
         LATT_ALT%A(3,:) =LATT_ALT%A(3,:)*SCALEZ*LATT_ALT%SCALE
   
         CALL LATTIC(LATT_ALT)
   
         IF (LATT_ALT%OMEGA<0) THEN
           CALL vtutor%error("ERROR: the triple product of the basis vectors is negative exchange two &
              &basis vectors")
         ENDIF
   
         T_INFO%NIOND =NIOND
         T_INFO%NIONPD=NIONPD
         T_INFO%NTYPD =NTYPD
         T_INFO%NTYPPD=NTYPPD
         ALLOCATE(T_INFO%LSFOR(3,NIOND),T_INFO%ITYP(NIOND))
   
         T_INFO%LSFOR=.TRUE.
   
   !-----number of atoms per type
         READ(15,'(A)') INPLIN
   
         READ(INPLIN,*) CSEL
         IF (.NOT.(CSEL>='0' .AND. CSEL<='9')) THEN
            READ(INPLIN,*) (T_INFO%TYPE(NI),NI=1,T_INFO%NTYP)
            READ(15,'(A)') INPLIN
            IF (T_INFO%NTYP/=NITEMS(INPLIN,INPWRK,.TRUE.,'I')) THEN
               CALL vtutor%error("ERROR: the type information is not consistent with the number of types")
            ENDIF
         ENDIF
   
         READ(INPLIN,*) (T_INFO%NITYP(NT),NT=1,T_INFO%NTYP)
   !---- Set up the table from which we get type of each ion
         NI=1
         DO NT=1,T_INFO%NTYP
         DO NI=NI,T_INFO%NITYP(NT)+NI-1
           T_INFO%ITYP(NI)=NT
         ENDDO
         ENDDO
   !
   !   positions
   !
         T_INFO%NIONS=0
         DO NT=1,T_INFO%NTYP
         T_INFO%NIONS= T_INFO%NIONS+ T_INFO%NITYP(NT)
         ENDDO
   
         T_INFO%NIONP=T_INFO%NIONS
   
         IF (T_INFO%NIONS>NIOND) THEN
           CALL vtutor%error("ERROR: MAIN: increase NIOND " // str(T_INFO%NIONS))
         ENDIF
   
         READ(15,'(A1)') CSEL
         T_INFO%LSDYN=((CSEL=='S').OR.(CSEL=='s'))
         IF (T_INFO%LSDYN) READ(15,'(A1)') CSEL
         IF (CSEL=='K'.OR.CSEL=='k'.OR. &
        &    CSEL=='C'.OR.CSEL=='c') THEN
           CSEL='K'
           IF (IU6>=0) &
           WRITE(IU6,*)' positions in cartesian coordinates'
   
           T_INFO%LDIRCO=.FALSE.
         ELSE
           IF (IU6>=0) &
           WRITE(IU6,*)' positions in direct lattice'
           T_INFO%LDIRCO=.TRUE.
         ENDIF
   
         ALLOCATE(DYN%POSION(3,NIONPD),DYN%POSIOC(3,NIONPD), &
        &         DYN%D2C(3,NIOND), &
        &         DYN%VEL(3,NIOND),DYN%D2(3,NIOND),DYN%D3(3,NIOND))
   
   ! alias T_INFO%POSION
         T_INFO%POSION => DYN%POSION
         DYN%POSION=0
         DYN%VEL   =0
         DYN%D2    =0
         DYN%D2C   =0
         DYN%D3    =0
   
         DO NI=1,T_INFO%NIONS
         IF (T_INFO%LSDYN) THEN
         READ(15,*,ERR=400,END=400) DYN%POSION(1,NI),DYN%POSION(2,NI),DYN%POSION(3,NI), &
        &      T_INFO%LSFOR(1,NI),T_INFO%LSFOR(2,NI),T_INFO%LSFOR(3,NI)
         ELSE
         READ(15,*,ERR=400,END=400) DYN%POSION(1,NI),DYN%POSION(2,NI),DYN%POSION(3,NI)
         ENDIF
         ENDDO
   
         IF (CSEL=='K') THEN
           DYN%POSION(1,:)=LATT_ALT%SCALE*DYN%POSION(1,:)*SCALEX
           DYN%POSION(2,:)=LATT_ALT%SCALE*DYN%POSION(2,:)*SCALEY
           DYN%POSION(3,:)=LATT_ALT%SCALE*DYN%POSION(3,:)*SCALEZ
           
           CALL KARDIR(T_INFO%NIONS,DYN%POSION,LATT_ALT%B)
         ENDIF
         CALL TOPRIM(T_INFO%NIONS,DYN%POSION)
         DYN%POSIOC=DYN%POSION
   
         DYN%INIT=0
         DYN%SNOSE(1)=1
   !
   !   empty spheres
   !
         READ(15,'(A1)',ERR=424,END=410) CSEL
     424 IF ((CSEL=='E').OR.(CSEL=='e')) THEN
           IF (T_INFO%NTYPP>NTYPPD) THEN
            CALL vtutor%error("ERROR: MAIN: increase NEMPTY " // str(T_INFO%NTYPP-T_INFO%NTYP))
           ENDIF
           READ(15,*,ERR=410,END=410) (T_INFO%NITYP(NT),NT=T_INFO%NTYP+1,T_INFO%NTYPP)
           DO NT=T_INFO%NTYP+1,T_INFO%NTYPP
             T_INFO%NIONP=T_INFO%NIONP+T_INFO%NITYP(NT)
           ENDDO
           IF (T_INFO%NIONP>NIONPD) THEN
            CALL vtutor%error("ERROR: MAIN: increase NEMPTY " // str(T_INFO%NIONP-T_INFO%NIONS))
           ENDIF
           T_INFO%NIONP=T_INFO%NIONP
   
           DO NI=T_INFO%NIONS+1,T_INFO%NIONP
            READ(15,*,ERR=410,END=410) &
        &      DYN%POSION(1,NI),DYN%POSION(2,NI),DYN%POSION(3,NI)
           ENDDO
           IF (.NOT.T_INFO%LDIRCO) THEN
             DO NI=T_INFO%NIONS+1,T_INFO%NIONP
               DYN%POSION(1,NI)=LATT_ALT%SCALE*DYN%POSION(1,NI)*SCALEX
               DYN%POSION(2,NI)=LATT_ALT%SCALE*DYN%POSION(2,NI)*SCALEY
               DYN%POSION(3,NI)=LATT_ALT%SCALE*DYN%POSION(3,NI)*SCALEZ
               CALL KARDIR(1,DYN%POSION(1:3,NI),LATT_ALT%B)
             ENDDO
           ENDIF
           READ(15,'(A1)',ERR=425,END=410) CSEL
         ENDIF
         DYN%POSIOC=DYN%POSION
   
     425 IF (CSEL=='K'.OR.CSEL=='k'.OR.CSEL==' ' &
        &    .OR.CSEL=='C'.OR.CSEL=='c') THEN
           CSEL='K'
           IF (IU6>=0) &
           WRITE(IU6,*)' velocities in cartesian coordinates'
         ELSE
           IF (IU6>=0) &
           WRITE(IU6,*)' velocities in direct lattice'
         ENDIF
   
   !
   !-----if we have velocities, read them in and transform from
   !     cartesian coordinates to direct lattice
         DO NI=1,T_INFO%NIONS
           READ(15,*,ERR=410,END=410)  &
        &            DYN%VEL(1,NI),DYN%VEL(2,NI),DYN%VEL(3,NI)
           !tb start
           IF (CSEL=='K' .AND. DYN%IBRION/=44 ) THEN
           !IF (CSEL=='K') THEN
           !tb end
           CALL  KARDIR(1,DYN%VEL(1:3,NI),LATT_ALT%B)
           DYN%VEL(1:3,NI)=DYN%VEL(1:3,NI)*DYN%POTIM
           ENDIF
         ENDDO
         
   !
   !-----try to read in predictor Coordinates
   !
         READ(15,*,ERR=430,END=430)
         READ(15,*,ERR=430,END=430) DYN%INIT
   
   !-----if INIT is there and it is 1 we have predictor-coordinates on the
   !-----file so we can start with them
         IF (DYN%INIT==0) GOTO 430
         READ(15,*) POTIMR
         IF (POTIMR/=DYN%POTIM) THEN
           IF (IU6>=0) THEN
              WRITE(IU6,*)
              WRITE(IU6,*)' There are predictor-coordinates on the file.'
              WRITE(IU6,*)' we can''t use them due to change of POTIM!'
           ENDIF
           GOTO 430
         ENDIF
   
   !-----Read in Nose-Parameter
         READ(15,*) DYN%SNOSE
   !-----Read in predictor-coordinates (always in direct lattice)
         READ(15,*) (DYN%POSION(1,NI),DYN%POSION(2,NI),DYN%POSION(3,NI),NI=1,T_INFO%NIONS)
         READ(15,*) (DYN%D2(1,NI),DYN%D2(2,NI),DYN%D2(3,NI),NI=1,T_INFO%NIONS)
         READ(15,*) (DYN%D3(1,NI),DYN%D3(2,NI),DYN%D3(3,NI),NI=1,T_INFO%NIONS)
         IF (IU6>=0) THEN
            WRITE(IU6,*)
            WRITE(IU6,*)' Using predictor-coordinates on the file'
         ENDIF
   
         CLOSE(UNIT=15)
         RETURN
   !-----------------------------------------------------------------------
   !  Reading Inputfile 15 finished
   !  if you end up at 430  INIT is set to 0,
   !    INIT is used in the call to STEP  (predictors are not initialised)
   !    in that way we tell STEP that it must initialize everything for us
   !----------------------------------------------------------------------
     400 CONTINUE
         CALL vtutor%error("No initial positions read in")
   
     410 CONTINUE
         DYN%INIT=-1
         IF (IU6>=0) &
         WRITE(IU6,*)' No initial velocities read in'
   
         CLOSE(UNIT=15)
         RETURN
   
     430 DYN%INIT=0
         CLOSE(UNIT=15)
         RETURN
   
      END SUBROUTINE
   

!*********************************************************************
!  
! this subroutine counts the dregress of freedom
!
!*********************************************************************

      SUBROUTINE COUNT_DEGREES_OF_FREEDOM( T_INFO, NDEGREES_OF_FREEDOM, &
              IU6, IU0, IBRION)
      USE prec
      TYPE (type_info) :: T_INFO
      COMPLEX (q) :: CDUM; REAL(q) :: RDUM ; LOGICAL :: LDUM
      INTEGER NDEGREES_OF_FREEDOM, IBRION, IU6, IU0, N, J

      IF ( T_INFO%LSDYN ) THEN
         NDEGREES_OF_FREEDOM=0
         DO N=1,T_INFO%NIONS
            DO J=1,3
               IF ( T_INFO%LSFOR(J,N)) &
                    NDEGREES_OF_FREEDOM=NDEGREES_OF_FREEDOM+1
            ENDDO
         ENDDO
      ELSE
         NDEGREES_OF_FREEDOM=3*T_INFO%NIONS-3
      ENDIF
! zero degrees of freedom, do not make me happy
! so in that case I simply set NDEGREES_OF_FREEDOM to 3*NIONS
! this avoids floating point exceptions in lots of places
      IF (NDEGREES_OF_FREEDOM==0) NDEGREES_OF_FREEDOM=3*T_INFO%NIONS
      
      END SUBROUTINE COUNT_DEGREES_OF_FREEDOM

!*************************** SYMVEL **********************************
!  
!  this subroutine removes any drift from the velocities
!  and warns the user if that the drift has been removed
!
!*********************************************************************

      SUBROUTINE SYMVEL_WARNING(NIONS, NTYP, ITYP,POMASS,V,IU6,IU0)
      USE prec
      USE lattice
      USE tutor, ONLY: vtutor, isAlert, CenterOfMassDrift

      INTEGER NTYP, NIONS, IU6, IU0, N, J, NT
      REAL(q) V(3,NIONS)
      INTEGER ITYP(NIONS)
      REAL(q) POMASS(NTYP)
      REAL(q) TMP(3),AVERAGE
      LOGICAL LWARN
      COMPLEX(q)  CDUM  ; LOGICAL  LDUM; REAL(q) RDUM; INTEGER IDUM


      AVERAGE=0
      TMP=0
!
      DO N=1,NIONS
         NT=ITYP(N)
         AVERAGE=AVERAGE+POMASS(NT)
         DO  J=1,3
            TMP(J)=TMP(J)+V(J,N)*POMASS(NT)
         ENDDO
      ENDDO

      LWARN=.FALSE.

      DO J=1,3
         IF ( ABS(TMP(J))> 1E-6) THEN
            LWARN=.TRUE.
         ENDIF
         TMP(J)=-TMP(J)/AVERAGE
      ENDDO

      DO N=1,NIONS
         DO J=1,3
            V(J,N)=V(J,N)+TMP(J)
         ENDDO
      ENDDO
      IF (LWARN) THEN
         CALL vtutor%write(isAlert, CenterOfMassDrift)
      ENDIF

      RETURN
      END SUBROUTINE


!*************************SUBROUTINE OUTPOS ****************************
!
!   write lattice parameters and positions to specified unit
!   use POSCAR compatibel format
!   LLONG specifies wether a long or short format is created
!   should be called only on IONODE !!
!
!***********************************************************************

      SUBROUTINE OUTPOS(IU, LLONG,SZNAM, T_INFO, SCALE, A, LSDYN, POSION)
      USE prec
      IMPLICIT NONE

      LOGICAL LLONG,LSDYN
      TYPE (type_info) :: T_INFO
      REAL(q) SCALE
      REAL(q) A(3,3)
      REAL(q) POSION(3,T_INFO%NIONS)
      CHARACTER (40) FORM
      CHARACTER (40) SZNAM
      INTEGER IU
! LOCAL
      INTEGER NT, NI, I
!-----direct lattice
      WRITE(IU,'(A40)') SZNAM

      WRITE(IU,*)  SCALE
      IF (LLONG) THEN
        FORM='(1X,3F22.16)'
      ELSE
        FORM='(1X,3F12.6)'
      ENDIF
      WRITE(IU,FORM) (A(1,I)/SCALE,A(2,I)/SCALE,A(3,I)/SCALE,I=1,3)

      IF (T_INFO%TYPE(1)/='  ') THEN
         WRITE(IU,'(20A5)') (T_INFO%TYPE(NT),NT=1,T_INFO%NTYP)
      ENDIF

      WRITE(IU,'(20I6)') (T_INFO%NITYP(NT),NT=1,T_INFO%NTYP)
      IF (LSDYN) WRITE(13,'(A18)') 'Selective dynamics'
      WRITE(IU,'(A6)')'Direct'

      IF (LSDYN) THEN
      IF (LLONG) THEN
        FORM='(3F20.16,3L4)'
      ELSE
        FORM='(3F10.6,3L2)'
      ENDIF
      ELSE
      IF (LLONG) THEN
        FORM='(3F20.16)'
      ELSE
        FORM='(3F10.6)'
      ENDIF
      ENDIF

      IF (LSDYN) THEN
          WRITE(IU,FORM) &
     &      (POSION(1,NI),POSION(2,NI),POSION(3,NI), &
     &       T_INFO%LSFOR(1,NI),T_INFO%LSFOR(2,NI),T_INFO%LSFOR(3,NI),NI=1,T_INFO%NIONS)
         ELSE
          WRITE(IU,FORM) &
     &      (POSION(1,NI),POSION(2,NI),POSION(3,NI),NI=1,T_INFO%NIONS)
      ENDIF
      IF (.NOT.LLONG) WRITE(IU,*)
      RETURN
      END SUBROUTINE

!=======================================================================
!
! read positions from stdin
!
!=======================================================================
      SUBROUTINE INPOS(LATT_CUR, T_INFO, DYN, IU6, IU0, LSTOP, MYCOMM)
      USE prec
      USE lattice
      USE main_mpi
      IMPLICIT NONE

      TYPE (latt)::       LATT_CUR
      TYPE (type_info) :: T_INFO
      TYPE (dynamics)  :: DYN
      INTEGER :: IU6, IU0
      LOGICAL :: LSTOP
      TYPE (communic) :: MYCOMM
    ! local
      INTEGER NI, IERR

      IF (IU6>=0) THEN
         IF (IU0>=0) WRITE(IU0,'(A)') 'POSITIONS: reading from stdin'
         DO NI=1,T_INFO%NIONS
            READ(*,*,  IOSTAT=IERR) DYN%POSION(1,NI), DYN%POSION(2,NI), DYN%POSION(3,NI)
            IF (IERR/=0) EXIT
         ENDDO

         IF (IERR/=0) THEN
            LSTOP=.TRUE.
         ELSE
            LSTOP=.FALSE.
         ENDIF
         WRITE(*,'(3F14.7)') DYN%POSION
         WRITE(*,'(3F14.7)') DYN%POSIOC
         CALLMPI( M_sum_i(MYCOMM, IERR, 1))
         IF (IERR==0) THEN
            CALLMPI( M_sum_d(MYCOMM, DYN%POSION(1,1), T_INFO%NIONS*3))
         ENDIF
         IF (IU0>=0) WRITE(IU0,'(A)') 'POSITIONS: read from stdin'
      ELSE
         IERR=0
         CALLMPI( M_sum_i(MYCOMM, IERR, 1))
         IF (IERR==0) THEN
            DYN%POSION(:,1:T_INFO%NIONS)=0
            CALLMPI( M_sum_d(MYCOMM, DYN%POSION(1,1), T_INFO%NIONS*3))
         ENDIF
      ENDIF
      END SUBROUTINE INPOS



!*************************SUBROUTINE OUTPOS_TRAIL  *********************
! write trailer for CONTCAR file
!
!    should be called only on IONODE !!
!***********************************************************************

      SUBROUTINE OUTPOS_TRAIL(IU,LOPEN, LATT_CUR, T_INFO, DYN)
      USE prec
      USE lattice
      IMPLICIT NONE

      INTEGER IU
      LOGICAL LOPEN
      TYPE (latt)::       LATT_CUR
      TYPE (type_info) :: T_INFO
      TYPE (dynamics)  :: DYN
! local variables
      INTEGER NT,NI
      REAL(q) :: VTMP(3)

      IF (DYN%IBRION==0 .AND. (DYN%ISIF==3 .OR. DYN%ISIF==8 .OR. DYN%ISIF==9 .OR. DYN%ISIF==10)) THEN
        !tb start new
        ! WRITE(IU,'(A)') 'Lattice velocities' 
        ! WRITE(IU,480) &
        !&      (LATT_CUR%AVEL(1,NI)/DYN%POTIM,LATT_CUR%AVEL(2,NI)/DYN%POTIM,LATT_CUR%AVEL(3,NI)/DYN%POTIM,NI=1,3)
        IF ( LATT_CUR%INITlatv==1) THEN
          WRITE(IU,'(A)') 'Lattice velocities and vectors'
          WRITE(IU,*)  LATT_CUR%INITlatv
          WRITE(IU,480) &
          &      (LATT_CUR%AVEL(1,NI)/DYN%POTIM,LATT_CUR%AVEL(2,NI)/DYN%POTIM,LATT_CUR%AVEL(3,NI)/DYN%POTIM,NI=1,3)
          WRITE(IU,480) &
          &      (LATT_CUR%A(1,NI),LATT_CUR%A(2,NI),LATT_CUR%A(3,NI),NI=1,3)
        ENDIF
        !tb end new
      ENDIF


      IF (T_INFO%NIONP>T_INFO%NIONS) THEN
         WRITE(IU,'(A)') 'Empty spheres'
         WRITE(IU,'(20I5)') (T_INFO%NITYP(NT),NT=T_INFO%NTYP+1,T_INFO%NTYPP)
         WRITE(IU,'(3F20.16)') &
     &      (DYN%POSION(1,NI),DYN%POSION(2,NI),DYN%POSION(3,NI),NI=T_INFO%NIONS+1,T_INFO%NIONP)
      ENDIF
      WRITE(IU,*)
!-----write out velocities
      DO NI=1,T_INFO%NIONS
      !tb start
        IF (DYN%IBRION/=44 .AND. DYN%IBRION/=40) THEN
          VTMP(1)=   DYN%VEL(1,NI)/DYN%POTIM
          VTMP(2)=   DYN%VEL(2,NI)/DYN%POTIM
          VTMP(3)=   DYN%VEL(3,NI)/DYN%POTIM
          CALL  DIRKAR(1,VTMP,LATT_CUR%A)
        ELSE
          VTMP(1)=   DYN%VEL(1,NI)
          VTMP(2)=   DYN%VEL(2,NI)
          VTMP(3)=   DYN%VEL(3,NI)
        END IF
        !VTMP(1)=   DYN%VEL(1,NI)/DYN%POTIM
        !VTMP(2)=   DYN%VEL(2,NI)/DYN%POTIM
        !VTMP(3)=   DYN%VEL(3,NI)/DYN%POTIM
        !CALL  DIRKAR(1,VTMP,LATT_CUR%A)
     !tb end
        WRITE(IU,480) VTMP(1),VTMP(2),VTMP(3)
      ENDDO
  480 FORMAT(3E16.8)
!-----if there was a call to STEP write out predictor-coordinates
      IF (DYN%INIT==1) THEN
      WRITE(IU,*)
      WRITE(IU,*) DYN%INIT
      WRITE(IU,*) DYN%POTIM
!-----write Nose-Parameter
      WRITE(IU,'(4E16.8)') DYN%SNOSE
      WRITE(IU,480) (DYN%POSION(1,NI),DYN%POSION(2,NI),DYN%POSION(3,NI),NI=1,T_INFO%NIONS)
      WRITE(IU,480) (DYN%D2(1,NI),DYN%D2(2,NI),DYN%D2(3,NI),NI=1,T_INFO%NIONS)
      WRITE(IU,480) (DYN%D3(1,NI),DYN%D3(2,NI),DYN%D3(3,NI),NI=1,T_INFO%NIONS)
      ENDIF

      IF (LOPEN) THEN
         CALL REOPEN(IU)
      ELSE
         REWIND IU
      ENDIF
      RETURN
      END SUBROUTINE

!***********************************************************************
!  write out initial header for XDATCAR
!***********************************************************************

      SUBROUTINE XDAT_HEAD(IU, T_INFO, A, SZNAM1)
      USE prec
      USE lattice
      IMPLICIT NONE

      INTEGER IU
      TYPE (type_info) :: T_INFO
      CHARACTER (40) SZNAM1
! local variables
      REAL(q) AOMEGA
      INTEGER NT
      INTEGER I
      CHARACTER (40) FORM
      REAL(q) :: A(3,3)
      INTEGER, PARAMETER :: SCALE=1
      LOGICAL, PARAMETER :: LLONG=.FALSE.

!
!     WRITE(IU,'(A40)') SZNAM1
!     IF (T_INFO%TYPE(1)/='  ') THEN
!        WRITE(IU,'(20A5)') (T_INFO%TYPE(NT),NT=1,T_INFO%NTYP)
!     ENDIF
      
      WRITE(IU,'(A40)') SZNAM1
      WRITE(IU,*)  SCALE
      IF (LLONG) THEN
        FORM='(1X,3F22.16)'
      ELSE
        FORM='(1X,3F12.6)'
      ENDIF
      WRITE(IU,FORM) (A(1,I)/SCALE,A(2,I)/SCALE,A(3,I)/SCALE,I=1,3)

      IF (T_INFO%TYPE(1)/='  ') THEN
         WRITE(IU,'(20A5)') (T_INFO%TYPE(NT),NT=1,T_INFO%NTYP)
      ENDIF

      WRITE(IU,'(20I6)') (T_INFO%NITYP(NT),NT=1,T_INFO%NTYP)

      RETURN
      END SUBROUTINE

!***********************************************************************
!
!  nearest neighboar table
!  (special wish of Roland Stumpf, and I think a good idea indeed)
!
!***********************************************************************

      SUBROUTINE NEAREST_NEIGHBOUR(IU6, IU0, T_INFO,L, RWIGS_)
      USE prec
      USE lattice
      USE tutor, ONLY: vtutor, isAlert, Position
      IMPLICIT NONE

      INTEGER IU6, IU0
      TYPE (latt)      :: L       ! lattice
      TYPE (type_info),TARGET :: T_INFO
      REAL(q) :: RWIGS_(T_INFO%NTYP)
! local variables
      INTEGER, PARAMETER :: MAXNEIG=100
      REAL(q),POINTER :: POSION(:,:)
      INTEGER I1,I2,I3,NII,NI,NIONS,NT,NTT,NSWP,NOUT,I,II,IND,IDUM,IDIOT
      REAL(q) D,DX,DY,DZ,RWIGS1,RWIGS2,DIS,SWP,RWIGS_MIN,RDUM
      COMPLEX(q) CDUM

      INTEGER NEIGT(T_INFO%NIONS,MAXNEIG),NEIGN(T_INFO%NIONS)
      REAL(q) DNEIG(T_INFO%NIONS,MAXNEIG)
      LOGICAL LWARN,LDUM
      INTEGER IREP
      REAL(q) :: RWIGS(T_INFO%NTYP)

      POSION => T_INFO%POSION
      LWARN=.FALSE.
!--------------------------------------------------------------------
! build up nearest neighbor table
!--------------------------------------------------------------------
      IF (IU6 < 0) RETURN

      NIONS=T_INFO%NIONS
      IF (NIONS<1000) THEN
         IREP=1
      ELSE
      ! large systems, does not make sense to go over replicated cells
         IREP=0
      ENDIF

      RWIGS_MIN=1000

      ! use RWIGS from INCAR, optionally 
      IF ( SUM( T_INFO%RWIGS ) > 0 ) THEN
         RWIGS= T_INFO%RWIGS 
      ELSE
         RWIGS= RWIGS_ 
      ENDIF

      DO NI=1,NIONS
         NT=T_INFO%ITYP(NI)
         RWIGS1=RWIGS(NT)*1.2
         IF (RWIGS1 <= 0.1) RWIGS1=1.0
         RWIGS_MIN=MIN(RWIGS1,RWIGS_MIN)

         IND=1
         I1=0
         I2=0
         I3=0

         DO I1=-IREP,IREP
         DO I2=-IREP,IREP
         DO I3=-IREP,IREP

            DO NII=1,NIONS
               NTT=T_INFO%ITYP(NII)
               RWIGS2=RWIGS(NTT)*1.2
               IF (RWIGS2 <= 0.1) RWIGS2=1.0

               DIS=RWIGS1+RWIGS2

               DX = MOD(POSION(1,NI)-POSION(1,NII)+10.5_q,1._q) -.5+I1
               DY = MOD(POSION(2,NI)-POSION(2,NII)+10.5_q,1._q) -.5+I2
               DZ = MOD(POSION(3,NI)-POSION(3,NII)+10.5_q,1._q) -.5+I3

               D =SQRT((DX*L%A(1,1)+DY*L%A(1,2)+DZ*L%A(1,3))**2 &
                  + (DX*L%A(2,1)+DY*L%A(2,2)+DZ*L%A(2,3))**2 &
                  + (DX*L%A(3,1)+DY*L%A(3,2)+DZ*L%A(3,3))**2)
               IF (NII/=NI .AND. D < DIS .AND. IND<=MAXNEIG) THEN
                  NEIGT(NI,IND)=NII
                  DNEIG(NI,IND)=D
                  IND=IND+1
               ENDIF
               !
               ! warn if distance between two ions is smaller than
               ! their averaged sum of their Wigner-Seitz radii
               IF ( NII/=NI .AND. D<DIS/2 ) LWARN = .TRUE.
            ENDDO
         ENDDO
         ENDDO
         ENDDO
         NEIGN(NI)=IND-1
      ENDDO
!--------------------------------------------------------------------
! sort by length
!--------------------------------------------------------------------
      WRITE(IU6,*) 'ion  position               nearest neighbor table'
      DO NI=1,NIONS
         DO I =1,NEIGN(NI)
         DO II=1,I-1
            IF (DNEIG(NI,I) < DNEIG(NI,II)) THEN
               NSWP       =NEIGT(NI,I)
               NEIGT(NI,I)=NEIGT(NI,II)
               NEIGT(NI,II)=NSWP
               SWP        =DNEIG(NI,I)
               DNEIG(NI,I)=DNEIG(NI,II)
               DNEIG(NI,II)=SWP
            ENDIF
         ENDDO
         ENDDO

         NOUT     =MIN(NEIGN(NI),16)
         IF (NIONS<1000) THEN
            ! short format
            WRITE(IU6,11) NI,POSION(:,NI),(NEIGT(NI,IND),DNEIG(NI,IND),IND=1,NOUT)
         ELSE
            ! long format
            WRITE(IU6,111) NI,POSION(:,NI),(NEIGT(NI,IND),DNEIG(NI,IND),IND=1,NOUT)
         ENDIF
 11      FORMAT(I4,3F7.3,'-',8(I4,F5.2),(/,26X,8(I4,F5.2)))
 111     FORMAT(I6,3F7.3,'-',8(I6,F5.2),(/,26X,8(I6,F5.2)))
      ENDDO
      WRITE(IU6,*)


      DO NI=1,NIONS
         IF (NEIGN(NI) > 0) THEN
            IF (DNEIG(NI,1) < RWIGS_MIN) LWARN=.TRUE.
         ENDIF
      ENDDO

      IDIOT=3
      IF (LWARN) THEN
         CALL vtutor%write(isAlert, Position)
      ENDIF

      END SUBROUTINE

      END MODULE
