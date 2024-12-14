#if defined(MPI) || defined(MPI_CHAIN)
#include "symbol.inc"

!****************************************************************************************************
! Module to handle files used in machine-learning calculations
! 1) ABCAR is an input file storing all the necessary ab initio data for generating force field.
! 2) INCAR is an input file that contains force field parameters.
! 3) ABNCAR is an output file that should contain all the renewed ab initio data.
! 4) REGCAR is an output file that contains the results of the regressions.
!****************************************************************************************************

      MODULE IOHANDLE
        USE ML_FF_CONSTANT
        USE ML_FF_PREC
        USE LATTICE_ML
        USE MPI_DATA
        USE PRINT2TERMINAL
        USE TIME_CHECK
        IMPLICIT NONE

        CONTAINS

!****************************************************************************************************
! Increment global file_stream counter 
!****************************************************************************************************

        SUBROUTINE ML_FILE_STREAM_INCREMENT(FILE_STREAM)
          IMPLICIT NONE
          INTEGER, INTENT(INOUT) :: FILE_STREAM
          FILE_STREAM = FILE_STREAM + 1
        END SUBROUTINE ML_FILE_STREAM_INCREMENT

!****************************************************************************************************
! Set ending of filenames for given machine learning instance
!****************************************************************************************************

        SUBROUTINE ML_SET_POSTFIX(INSTANCE,FILE_POSTFIX)
          IMPLICIT NONE
          INTEGER, INTENT(IN)           :: INSTANCE
          CHARACTER(LEN=4), INTENT(OUT) :: FILE_POSTFIX
          CHARACTER(LEN=1)              :: H1
          CHARACTER(LEN=2)              :: H2
          CHARACTER(LEN=3)              :: H3
          IF (INSTANCE.EQ.1) THEN
             FILE_POSTFIX=""
          ELSE
             IF (INSTANCE.LT.10) THEN
                WRITE(H1,"(I1)") INSTANCE
                FILE_POSTFIX="_" // TRIM(H1)
             ELSE IF (INSTANCE.GT.9.AND.INSTANCE.LT.100) THEN
                WRITE(H2,"(I2)") INSTANCE
                FILE_POSTFIX="_" // TRIM(H2)
             ELSE IF (INSTANCE.GT.99.AND.INSTANCE.LT.1000) THEN
                WRITE(H3,"(I3)") INSTANCE
                FILE_POSTFIX="_" // TRIM(H3)
             ENDIF
          ENDIF
        END SUBROUTINE ML_SET_POSTFIX

!****************************************************************************************************
! Set file name and code
!****************************************************************************************************

        SUBROUTINE SET_FILE (CFILE_IN,IFILE_IN,CFILE_OUT,IFILE_OUT)
          IMPLICIT NONE
          CHARACTER(len=*), INTENT(IN)  :: CFILE_IN
          CHARACTER(len=*), INTENT(OUT) :: CFILE_OUT
          INTEGER IFILE_IN
          INTEGER IFILE_OUT
          CFILE_OUT=CFILE_IN
          IFILE_OUT=IFILE_IN
        END SUBROUTINE SET_FILE

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Read element number written in ABCAR file.
!****************************************************************************************************

        SUBROUTINE READ_ABCAR_ELEMENT_ORDER1 (ABCAR,IU,MTYP_OLD)
          IMPLICIT NONE
! Input variables.
          CHARACTER(LEN=9), INTENT(IN)    :: ABCAR
          INTEGER         , INTENT(IN)    :: IU
! Input and ouput variables.
          INTEGER         , INTENT(INOUT) :: MTYP_OLD
! Local variables.
          INTEGER                         :: I
          INTEGER                         :: IERR
! Open file.
          OPEN(UNIT=IU,FILE=ABCAR,STATUS='OLD',IOSTAT=IERR)
! Check if error occured in opening file
          IF (IERR.NE.0) THEN
             CALL ERROR_PRINT2TERMINAL("Error in opening ML_AB file")
          ENDIF
! Read file.
          DO I=1, 8
             READ(IU,*)
          ENDDO
          READ(IU,*) MTYP_OLD
! Close file.
          CLOSE(IU)
        END SUBROUTINE READ_ABCAR_ELEMENT_ORDER1

!****************************************************************************************************
! Read element data written in ABCAR file.
!****************************************************************************************************

        SUBROUTINE READ_ABCAR_ELEMENT_ORDER2 (MB_INPUT,MCONF_INPUT,ABCAR,EATOM_OLD,IU,MIONS_OLD,MITYP_OLD,MTYP_OLD,POMASS_OLD,TYPE_OLD)
          IMPLICIT NONE
! Input variables.
          INTEGER         , INTENT(IN)    :: MB_INPUT
          INTEGER         , INTENT(IN)    :: MCONF_INPUT
          CHARACTER(LEN=9), INTENT(IN)    :: ABCAR
          INTEGER         , INTENT(IN)    :: IU
          INTEGER         , INTENT(IN)    :: MTYP_OLD
! Input and ouput variables.
          REAL(q)         , INTENT(INOUT) :: EATOM_OLD(:) !(1:MTYP_OLD)
          INTEGER         , INTENT(INOUT) :: MIONS_OLD
          INTEGER         , INTENT(INOUT) :: MITYP_OLD
          REAL(q)         , INTENT(INOUT) :: POMASS_OLD(:) !(1:MTYP_OLD)
          CHARACTER(LEN=2), INTENT(INOUT) :: TYPE_OLD(:) !(1:MTYP_OLD)
! Local variables
          INTEGER                         :: I
          INTEGER                         :: ICOL
          INTEGER                         :: IERR
          INTEGER                         :: IROW
          INTEGER                         :: MCONF_AB
          INTEGER                         :: NCOL
          INTEGER                         :: NMOD
          INTEGER                         :: NUM_BASIS_OLD(1:MTYP_OLD)
          INTEGER                         :: NWRITE
! Set constants.
          NWRITE=3
          NCOL=MTYP_OLD/NWRITE
          NMOD=MOD(MTYP_OLD,NWRITE)
! Open file.
          OPEN(UNIT=IU,FILE=ABCAR,STATUS='OLD',IOSTAT=IERR)
! Check if error occured in opening file
          IF (IERR.NE.0) THEN
             CALL ERROR_PRINT2TERMINAL("Error in opening ML_AB file")
          ENDIF
! Read the element data.
          DO I=1, 4
             READ(IU,*)
          ENDDO
          READ(IU,*) MCONF_AB
          DO I=1, 4
             READ(IU,*)
          ENDDO
          IF (MCONF_AB.GT.MCONF_INPUT) THEN
             CALL ERROR_PRINT2TERMINAL("ML_MCONF is smaller than the number &
                &of structures in the ML_AB file, please increase ML_MCONF.")
          ENDIF
          DO I=1, 3
             READ(IU,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(IU,*) (TYPE_OLD(NWRITE*(ICOL-1)+IROW),IROW=1,NWRITE)
          ENDDO
          IF(NMOD.NE.0) THEN
             READ(IU,*) (TYPE_OLD(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO I=1, 3
             READ(IU,*)
          ENDDO
          READ(IU,*) MIONS_OLD
          DO I=1, 3
             READ(IU,*)
          ENDDO
          READ(IU,*) MITYP_OLD
          DO I=1, 3
             READ(IU,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(IU,*) (EATOM_OLD(NWRITE*(ICOL-1)+IROW),IROW=1,NWRITE)
          ENDDO
          IF(NMOD.NE.0) THEN
             READ(IU,*) (EATOM_OLD(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO I=1, 3
             READ(IU,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(IU,*) (POMASS_OLD(NWRITE*(ICOL-1)+IROW),IROW=1,NWRITE)
          ENDDO
          IF(NMOD.NE.0) THEN
             READ(IU,*) (POMASS_OLD(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO I=1, 3
             READ(IU,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(IU,*) (NUM_BASIS_OLD(NWRITE*(ICOL-1)+IROW),IROW=1,NWRITE)
          ENDDO
          IF(NMOD.NE.0) THEN
             READ(IU,*) (NUM_BASIS_OLD(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO I=1, MTYP_OLD
             IF (NUM_BASIS_OLD(I).GT.MB_INPUT) THEN
                CALL ERROR_PRINT2TERMINAL("ML_MB is smaller than the number &
                   &of local configurations in the ML_AB file, please &
                   &increase ML_MB.")
             ENDIF
          ENDDO
! Convert units already immediately after reading in
          EATOM_OLD(1:MTYP_OLD)=EATOM_OLD(1:MTYP_OLD)/EUNIT
          POMASS_OLD(1:MTYP_OLD)=POMASS_OLD(1:MTYP_OLD)/MUNIT
! Close file.
          CLOSE(IU)
        END SUBROUTINE READ_ABCAR_ELEMENT_ORDER2
#endif

!****************************************************************************************************
! This routine determines whether the ML_FF file is in binary format or formatted string format.
!****************************************************************************************************
        SUBROUTINE CHECK_ML_FF_FORMAT(FF)
           USE ML_FF_STRUCT, ONLY: FF_PAR, ML_IO_WRITE
           
           ! Input and output variables
           TYPE (FF_PAR)                   :: FF
           ! Local variables
           CHARACTER                       :: HELP_CHAR
           INTEGER                         :: IERR

           ! Set binary format for ML_FF file by default
           FF%ML_FF_BINARY=.TRUE.
           ! Open ML_FF file as unformatted
           OPEN(UNIT=FF%IU1,FILE=FF%FFCAR,STATUS='OLD', &
                FORM='UNFORMATTED',ACCESS='STREAM',IOSTAT=IERR)
           ! READ first character
           READ(FF%IU1) HELP_CHAR
           IF (HELP_CHAR == "*") THEN
              FF%ML_FF_BINARY=.FALSE.
           ENDIF
           CLOSE(FF%IU1)
        END SUBROUTINE


!****************************************************************************************************
! Read element number written in FFCAR file.
!****************************************************************************************************

        SUBROUTINE READ_FF_CAR_DEF_ELEMENT_ORDER_FF_1(FFCAR,IU,MTYP_OLD)
          USE ML_FF_STRUCT, ONLY: ML_IO_WRITE
          IMPLICIT NONE
! Input variables.
          CHARACTER(LEN=9), INTENT(IN)    :: FFCAR
          INTEGER         , INTENT(IN)    :: IU
! Output variables.
          INTEGER         , INTENT(OUT)   :: MTYP_OLD
! Local variables
          INTEGER                         :: I
          INTEGER                         :: IERR
          INTEGER                         :: IHELP(1:3)
! Open file.
          OPEN(UNIT=IU,FILE=FFCAR,STATUS='OLD',FORM='UNFORMATTED',ACCESS='STREAM',IOSTAT=IERR)
! Check if error occured in opening file
          IF (IERR.NE.0) THEN
             CALL ERROR_PRINT2TERMINAL("Error in opening ML_FF file")
          ENDIF
! Do actual reading
          READ(IU) IHELP(1), IHELP(2), IHELP(3)
          READ(IU) MTYP_OLD
! Close file.
          CLOSE(IU)
        END SUBROUTINE READ_FF_CAR_DEF_ELEMENT_ORDER_FF_1

!****************************************************************************************************
! Read element data written in ABCAR file.
!****************************************************************************************************
        SUBROUTINE READ_FF_CAR_DEF_ELEMENT_ORDER_FF_2(FFCAR,IU,EATOM_OLD,MTYP_OLD,POMASS_OLD,TYPE_OLD)
          IMPLICIT NONE
! Input variables.
          CHARACTER(LEN=9), INTENT(IN)    :: FFCAR
          INTEGER         , INTENT(IN)    :: IU
          INTEGER         , INTENT(IN)    :: MTYP_OLD
! Input and ouput variables.
          REAL(q)         , INTENT(OUT)   :: EATOM_OLD(:) !(1:MTYP_OLD)
          REAL(q)         , INTENT(OUT)   :: POMASS_OLD(:) !(1:MTYP_OLD)
          CHARACTER(LEN=2), INTENT(OUT)   :: TYPE_OLD(:) !(1:MTYP_OLD)
! Local variables
          INTEGER                         :: I
          INTEGER                         :: ICOL
          INTEGER                         :: IERR
          INTEGER                         :: INTYP
          INTEGER                         :: IROW
          INTEGER                         :: NCOL
          INTEGER                         :: NMOD
          INTEGER                         :: NWRITE
          INTEGER                         :: IDUMMY
          INTEGER                         :: IHELP(1:3)
! Set constant
          NWRITE = 3
! Open file.
          OPEN(UNIT=IU,FILE=FFCAR,STATUS='OLD',FORM='UNFORMATTED',ACCESS='STREAM',IOSTAT=IERR)
! Check if error occured in opening file
          IF (IERR.NE.0) THEN
             CALL ERROR_PRINT2TERMINAL("Error in opening ML_FF file")
          ENDIF
! Do actual reading
          READ(IU) IHELP(1), IHELP(2), IHELP(3)
          READ(IU) IDUMMY
          NCOL=MTYP_OLD/NWRITE
          NMOD=MOD(MTYP_OLD,NWRITE)
          DO ICOL=1, NCOL, 1
             READ(IU) (TYPE_OLD(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(IU) (TYPE_OLD(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO ICOL=1, NCOL, 1
             READ(IU) (EATOM_OLD(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(IU) (EATOM_OLD(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO ICOL=1, NCOL
             READ(IU) (POMASS_OLD(NWRITE*(ICOL-1)+IROW),IROW=1,NWRITE)
          ENDDO
          IF(NMOD.NE.0) THEN
             READ(IU) (POMASS_OLD(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
! Unit conversion
          EATOM_OLD(1:MTYP_OLD)=EATOM_OLD(1:MTYP_OLD)/EUNIT
          POMASS_OLD(1:MTYP_OLD)=POMASS_OLD(1:MTYP_OLD)/MUNIT
! Close file.
          CLOSE(IU)
        END SUBROUTINE READ_FF_CAR_DEF_ELEMENT_ORDER_FF_2

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Read the number of configurations in the ABCAR file.
!****************************************************************************************************

        SUBROUTINE READ_ABCAR_SYSTEM_ORDER1 (ABCAR,IU,NCONF)
          IMPLICIT NONE
! Input variables.
          CHARACTER(LEN=9), INTENT(IN)    :: ABCAR
          INTEGER         , INTENT(IN)    :: IU
! Input and output variables.
          INTEGER         , INTENT(INOUT) :: NCONF
! Local variables.
          INTEGER                         :: I
          INTEGER                         :: IERR
! Open file.
          OPEN(UNIT=IU,FILE=ABCAR,STATUS='OLD',IOSTAT=IERR)
! Check if error occured in opening file
          IF (IERR.NE.0) THEN
             CALL ERROR_PRINT2TERMINAL("Error in opening ML_AB file")
          ENDIF
          DO I=1, 4
             READ(IU,*)
          ENDDO
          READ(IU,*) NCONF
! Close file.
          CLOSE(IU)
        END SUbROUTINE READ_ABCAR_SYSTEM_ORDER1

!****************************************************************************************************
! Read the system data in the ABCAR file.
!****************************************************************************************************

        SUBROUTINE READ_ABCAR_SYSTEM_ORDER2 (ABCAR,IU,NCONF,NSYS,SZNAM2_HELP2)
          IMPLICIT NONE
! Input variables.
          CHARACTER(LEN=9) , INTENT(IN)    :: ABCAR
          INTEGER          , INTENT(IN)    :: IU
          INTEGER          , INTENT(IN)    :: NCONF
! Input and output variables.
          INTEGER          , INTENT(INOUT) :: NSYS
          CHARACTER(LEN=40), INTENT(INOUT) :: SZNAM2_HELP2(:) !(1:NCONF+1)
! Local variables.
          INTEGER                          :: I
          INTEGER                          :: ICOL
          INTEGER                          :: ICONF
          INTEGER                          :: IERR
          INTEGER                          :: IFLAG
          INTEGER                          :: INTYP
          INTEGER                          :: IROW
          INTEGER                          :: ISYS
          INTEGER                          :: MTYP
          INTEGER          , ALLOCATABLE   :: NB(:)
          INTEGER                          :: NCOL
          INTEGER                          :: NIONS
          INTEGER                          :: NMOD
          INTEGER                          :: NTYP
          INTEGER                          :: NWRITE
          CHARACTER(LEN=40)                :: SZNAM2_HELP1
          REAL(q)                          :: VERSION
          CHARACTER(LEN=40)                :: READ_HELPER
! Set constants.
          NWRITE=3
! Open file.
          OPEN(UNIT=IU,FILE=ABCAR,STATUS='OLD',IOSTAT=IERR)
! Check if error occured in opening file
          IF (IERR.NE.0) THEN
             CALL ERROR_PRINT2TERMINAL("Error in opening ML_AB file")
          ENDIF
! Read the system data written in the ABCAR file.
          NSYS=0
          READ(IU,*) VERSION
          DO I=1, 7
             READ(IU,*)
          ENDDO
          READ(IU,*) MTYP
          ALLOCATE(NB(1:MTYP))
          NCOL=MTYP/NWRITE
          NMOD=MOD(MTYP,NWRITE)
          DO I=1, 3
             READ(IU,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(IU,*)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(IU,*)
          ENDIF
          DO I=1, 11
             READ(IU,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(IU,*)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(IU,*)
          ENDIF
          DO I=1, 3
             READ(IU,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(IU,*)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(IU,*)
          ENDIF
          DO I=1, 3
             READ(IU,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(IU,*) (NB(NWRITE*(ICOL-1)+IROW), IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(IU,*) (NB(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO INTYP=1, MTYP
             DO I=1, (NB(INTYP)+3)
                READ(IU,*)
             ENDDO
          ENDDO
          DO ICONF=1, NCONF
             DO I=1, 5
                READ(IU,*)
             ENDDO
             READ(IU,'(A40)') SZNAM2_HELP1
             SZNAM2_HELP1=TRIM(ADJUSTL(SZNAM2_HELP1))
             IF(NSYS.EQ.0) THEN
                NSYS=NSYS+1
                SZNAM2_HELP2(NSYS)=SZNAM2_HELP1
             ELSE
                IFLAG=1
                DO ISYS=1, NSYS
                   IF(SZNAM2_HELP2(ISYS).EQ.SZNAM2_HELP1) THEN
                      IFLAG=0
                      EXIT
                   ENDIF
                ENDDO
                IF(IFLAG.EQ.1) THEN
                   NSYS=NSYS+1
                   SZNAM2_HELP2(NSYS)=SZNAM2_HELP1
                ENDIF
             ENDIF
             DO I=1, 3
                READ(IU,*)
             ENDDO
             READ(IU,*) NTYP
             DO I=1, 3
                READ(IU,*)
             ENDDO
             READ(IU,*) NIONS
             DO I=1, 3
                READ(IU,*)
             ENDDO
             DO I=1, NTYP
                READ(IU,*)
             ENDDO
             ! Optional reading of CTIFOR
             READ(IU,*)
             READ(IU,*) READ_HELPER
             READ(IU,*)
             READ_HELPER=ADJUSTL(TRIM(READ_HELPER))
             IF (READ_HELPER(1:1).EQ.'C') THEN
                READ(IU,*)
                DO I=1, 3
                   READ(IU,*)
                ENDDO
             ENDIF
             ! End optional reading of CTIFOR
             DO I=1,3
                READ(IU,*)
             ENDDO
             IF (VERSION.GT.0.0_q) THEN
                DO I=1, 2*NIONS+20
                   READ(IU,*)
                ENDDO
             ELSE
                DO I=1, 3*NIONS+20
                   READ(IU,*)
                ENDDO
             ENDIF
          ENDDO
! Deallocate array.
          DEALLOCATE(NB)
! Close file.
          CLOSE(IU)
        END SUbROUTINE READ_ABCAR_SYSTEM_ORDER2

!****************************************************************************************************
! Read size written in ABCAR file.
!****************************************************************************************************

        SUBROUTINE READ_ABCAR_SIZE1 (ABCAR,IU,LMLMB_FF,MDATA_FFM,MTYP,NCONF)
          IMPLICIT NONE
! Input variables.
          CHARACTER(LEN=9), INTENT(IN)    :: ABCAR
          INTEGER         , INTENT(IN)    :: IU
          LOGICAL         , INTENT(IN)    :: LMLMB_FF
          INTEGER         , INTENT(IN)    :: MTYP
! In put and output variables.
          INTEGER         , INTENT(INOUT) :: MDATA_FFM
          INTEGER         , INTENT(INOUT) :: NCONF
! Local variables.
          INTEGER                         :: I
          INTEGER                         :: ICOL
          INTEGER                         :: ICONF
          INTEGER                         :: IERR
          INTEGER                         :: INTYP
          INTEGER                         :: IROW
          INTEGER                         :: MIONS
          INTEGER                         :: MITYP
          INTEGER                         :: MTYP_LOC
          INTEGER                         :: NB(1:MTYP)
          INTEGER                         :: NCOL
          INTEGER                         :: NMOD
          INTEGER                         :: NIONS
          INTEGER                         :: NTYP
          INTEGER                         :: NWRITE
          REAL(q)                         :: VERSION
          CHARACTER(LEN=40)               :: READ_HELPER
! Set constants.
          NWRITE=3
! Read ABCAR file.
          OPEN(UNIT=IU,FILE=ABCAR,STATUS='OLD',IOSTAT=IERR)
! Check if error occured in opening file
          IF (IERR.NE.0) THEN
             CALL ERROR_PRINT2TERMINAL("Error in opening ML_AB file")
          ENDIF
! Initialization.
          MDATA_FFM=0
! Read header.
          READ(IU,*) VERSION
          DO I=1, 3
             READ(IU,*)
          ENDDO
          READ(IU,*) NCONF
          DO I=1, 3
             READ(IU,*)
          ENDDO
          READ(IU,*) MTYP_LOC
! Set NCOL and NMOD.
          NCOL=MTYP_LOC/NWRITE
          NMOD=MOD(MTYP_LOC,NWRITE)
! Read ABCAR further
          DO I=1, 3
             READ(IU,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(IU,*)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(IU,*)
          ENDIF
          DO I=1, 3
             READ(IU,*)
          ENDDO
          READ(IU,*) MIONS
          DO I=1, 3
             READ(IU,*)
          ENDDO
          READ(IU,*) MITYP
          DO I=1, 3
             READ(IU,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(IU,*)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(IU,*)
          ENDIF
          DO I=1, 3
             READ(IU,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(IU,*)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(IU,*)
          ENDIF
! Read ABCAR file further
          DO I=1, 3
             READ(IU,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(IU,*) (NB(NWRITE*(ICOL-1)+IROW), IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(IU,*) (NB(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO INTYP=1, MTYP_LOC
             DO I=1, (NB(INTYP)+3)
                READ(IU,*)
             ENDDO
          ENDDO
          DO ICONF=1, NCONF
             DO I=1, 9
                READ(IU,*)
             ENDDO
             READ(IU,*) NTYP
             DO I=1, 3
                READ(IU,*)
             ENDDO
             READ(IU,*) NIONS
             IF (LMLMB_FF) THEN
                MDATA_FFM=MDATA_FFM+3*NIONS+7
             ENDIF
             DO I=1, 3
                READ(IU,*)
             ENDDO
             DO I=1, NTYP
                READ(IU,*)
             ENDDO
             ! Optional reading of CTIFOR
             READ(IU,*)
             READ(IU,*) READ_HELPER
             READ(IU,*)
             READ_HELPER=ADJUSTL(TRIM(READ_HELPER))
             IF (READ_HELPER(1:1).EQ.'C') THEN
                READ(IU,*)
                DO I=1, 3
                   READ(IU,*)
                ENDDO
             ENDIF
             ! End optional reading of CTIFOR
             DO I=1,3
                READ(IU,*)
             ENDDO
             IF (VERSION.GT.0.0_q) THEN
                DO I=1, 2*NIONS+20
                   READ(IU,*)
                ENDDO
             ELSE
                DO I=1, 3*NIONS+20
                   READ(IU,*)
                ENDDO
             ENDIF
          ENDDO
          CLOSE(IU)
        END SUBROUTINE READ_ABCAR_SIZE1

!****************************************************************************************************
! Read size data in ABCAR file
! CTIFOR_FF, IERR_FF are overwritten.
!****************************************************************************************************
        SUBROUTINE READ_ABCAR_SIZE2 (ML_INPUT_HANDLE,AB,FF,MITYP_NEW)
          USE ML_FF_STRUCT, ONLY : ML_FF_INPUT,ABINITIO,FF_PAR
          IMPLICIT NONE
! Input and output variables
          TYPE (ML_FF_INPUT)              :: ML_INPUT_HANDLE
          TYPE (ABINITIO)                 :: AB
          TYPE (FF_PAR)                   :: FF
! Input variables.
          INTEGER         , INTENT(IN)    :: MITYP_NEW
! Local variables.
          REAL(q),ALLOCATABLE             :: EATOM_TMP(:)
          INTEGER                         :: I
          INTEGER                         :: ICOL
          INTEGER                         :: ICONF
          INTEGER                         :: IERR
          INTEGER                         :: IFLAG
          INTEGER                         :: INTYP
          INTEGER                         :: IROW
          INTEGER                         :: ISYS
          INTEGER                         :: JNTYP
          INTEGER                         :: MIONS_TMP
          INTEGER                         :: MITYP_TMP
          INTEGER                         :: MTYP_TMP
          INTEGER,ALLOCATABLE             :: NB_TMP(:)
          INTEGER                         :: NCOL
          INTEGER                         :: NIONS_TMP
          INTEGER,ALLOCATABLE             :: NITYP_TMP2(:)
          INTEGER                         :: NMOD
          INTEGER                         :: NTYP_TMP
          INTEGER                         :: NWRITE
          REAL(q),ALLOCATABLE             :: POMASS_TMP(:)
          CHARACTER(LEN=40)               :: READ_HELPER
          CHARACTER(LEN=40)               :: SZNAM2_HELP1
          CHARACTER(LEN=40),ALLOCATABLE   :: SZNAM2_HELP2(:)
          CHARACTER(LEN=2),ALLOCATABLE    :: TYPE_TMP(:)
          REAL(q)                         :: VERSION
! Set constants
          NWRITE=3
!          AB%NSYS=FF%NSYS
!          AB%MTYP=FF%MTYP
          AB%NSYS=0
! Open AB%ABCAR file.
          OPEN(UNIT=AB%IU1,FILE=AB%ABCAR,STATUS='OLD',IOSTAT=IERR)
! Check if error occured in opening file
          IF (IERR.NE.0) THEN
             CALL ERROR_PRINT2TERMINAL("Error in opening ML_AB file")
          ENDIF
! Read header
          READ(AB%IU1,*) VERSION
          DO I=1, 3
             READ(AB%IU1,*)
          ENDDO
          READ(AB%IU1,*) AB%NCONF
          IF(AB%NCONF.GT.FF%MCONF) THEN
             CALL ERROR_PRINT2TERMINAL ('!!! MLFF : Set FF%MCONF &
                  larger than NCONF in AB%ABCAR file !!!')
          ENDIF
          DO I=1, 3
             READ(AB%IU1,*)
          ENDDO
          READ(AB%IU1,*) AB%MTYP
! Allocate helping arrays
          ALLOCATE(TYPE_TMP(AB%MTYP))
          ALLOCATE(EATOM_TMP(AB%MTYP))
          ALLOCATE(NB_TMP(AB%MTYP))
          ALLOCATE(NITYP_TMP2(AB%MTYP))
          ALLOCATE(POMASS_TMP(AB%MTYP))
          ALLOCATE(SZNAM2_HELP2(AB%NCONF))
! Set NCOL and NMOL to read elemeent data.
          NCOL=AB%MTYP/NWRITE
          NMOD=MOD(AB%MTYP,NWRITE)
! Read AB%ABCAR further
          DO I=1, 3
             READ(AB%IU1,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(AB%IU1,*) (TYPE_TMP(NWRITE*(ICOL-1)+IROW), IROW=1, NWRITE)
          ENDDO
          IF(NMOD.NE.0) THEN
             READ(AB%IU1,*) (TYPE_TMP(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO I=1, 3
             READ(AB%IU1,*)
          ENDDO
          READ(AB%IU1,*) MIONS_TMP
          DO I=1, 3
             READ(AB%IU1,*)
          ENDDO
          READ(AB%IU1,*) MITYP_TMP
          DO I=1, 3
             READ(AB%IU1,*)
          ENDDO
          ! Read atomic reference energies
          DO ICOL=1, NCOL
             READ(AB%IU1,*) (EATOM_TMP(NWRITE*(ICOL-1)+IROW), IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(AB%IU1,*) (EATOM_TMP(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO I=1, 3
             READ(AB%IU1,*)
          ENDDO
          ! Read atomic weights
          DO ICOL=1, NCOL
             READ(AB%IU1,*) (POMASS_TMP(NWRITE*(ICOL-1)+IROW), IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(AB%IU1,*) (POMASS_TMP(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO I=1, 3
             READ(AB%IU1,*)
          ENDDO
          ! Read number of local configurations
          DO ICOL=1, NCOL
             READ(AB%IU1,*) (NB_TMP(NWRITE*(ICOL-1)+IROW), IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(AB%IU1,*) (NB_TMP(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          ! Dummy read local configurations
          DO INTYP=1, AB%MTYP
             DO I=1, (NB_TMP(INTYP)+3)
                READ(AB%IU1,*)
             ENDDO
          ENDDO
          ! Loop over all configurations in ML_AB file
          DO ICONF=1, AB%NCONF
             DO I=1, 5
                READ(AB%IU1,*)
             ENDDO
             READ(AB%IU1,'(A40)') SZNAM2_HELP1
             SZNAM2_HELP1=TRIM(ADJUSTL(SZNAM2_HELP1))
             IF(AB%NSYS.EQ.0) THEN
                AB%NSYS=AB%NSYS+1
                SZNAM2_HELP2(AB%NSYS)=SZNAM2_HELP1
             ELSE
                IFLAG=1
                DO ISYS=1, AB%NSYS
                   IF(SZNAM2_HELP2(ISYS).EQ.SZNAM2_HELP1) THEN
                      IFLAG=0
                      EXIT
                   ENDIF
                ENDDO
                IF(IFLAG.EQ.1) THEN
                   AB%NSYS=AB%NSYS+1
                   SZNAM2_HELP2(AB%NSYS)=SZNAM2_HELP1
                ENDIF
             ENDIF
             DO I=1, 3
                READ(AB%IU1,*)
             ENDDO
             READ(AB%IU1,*) NTYP_TMP
             DO I=1, 3
                READ(AB%IU1,*)
             ENDDO
             READ(AB%IU1,*) NIONS_TMP
             DO I=1, 3
                READ(AB%IU1,*)
             ENDDO
             DO I=1, NTYP_TMP
                READ(AB%IU1,*)
             ENDDO
             ! Optional reading of CTIFOR
             READ(AB%IU1,*)
             READ(AB%IU1,*) READ_HELPER
             READ(AB%IU1,*)
             READ_HELPER=ADJUSTL(TRIM(READ_HELPER))
             IF (READ_HELPER(1:1).EQ.'C') THEN
                READ(AB%IU1,*)
                DO I=1, 3
                   READ(AB%IU1,*)
                ENDDO
             ENDIF
             ! End optional reading of CTIFOR
             DO I=1,3
                READ(AB%IU1,*)
             ENDDO
             IF (VERSION.GT.0.0_q) THEN
                DO I=1, 2*NIONS_TMP+20
                   READ(AB%IU1,*)
                ENDDO
             ELSE
                DO I=1, 3*NIONS_TMP+20
                   READ(AB%IU1,*)
                ENDDO
             ENDIF
          ENDDO
          CLOSE(AB%IU1)
! Define sizes
          AB%MCONF=FF%MCONF
          IF (FF%ISTART.NE.3) THEN
             AB%MIONS=MAX(MIONS_TMP,ML_INPUT_HANDLE%NIONS)
             AB%MITYP=MAX(MITYP_TMP,MITYP_NEW)
! Also check if current system is not inside AB file
             IFLAG=1
             DO ISYS=1, AB%NSYS
                IF(ML_INPUT_HANDLE%SZNAM2.EQ.SZNAM2_HELP2(ISYS)) THEN
                   IFLAG=0
                   EXIT
                ENDIF
             ENDDO
             IF(IFLAG.EQ.1) THEN
                AB%NSYS=AB%NSYS+1
             ENDIF
! Check if current systems types are among all types
             MTYP_TMP=AB%MTYP
             DO INTYP=1, ML_INPUT_HANDLE%NTYP
                IFLAG=1
                DO JNTYP=1, AB%MTYP
                   IF (TYPE_TMP(JNTYP).EQ.ML_INPUT_HANDLE%TYPE(INTYP)) THEN
                      IFLAG=0
                      EXIT
                   ENDIF
                ENDDO
                IF (IFLAG.EQ.1) THEN
                   MTYP_TMP=MTYP_TMP+1
                ENDIF
             ENDDO
             AB%MTYP=MTYP_TMP
! Obtain maximum total number of ions in systemfrom current structure and ML_AB file
             AB%MIONS=MAX(AB%MIONS,ML_INPUT_HANDLE%NIONS)
! Obtain maximum number of atoms per atom type
             DO INTYP=1, ML_INPUT_HANDLE%NTYP
                MITYP_TMP=MAX(MITYP_TMP,ML_INPUT_HANDLE%NITYP(INTYP))
             ENDDO
             AB%MITYP=MITYP_TMP
          ELSE
             AB%MIONS=MIONS_TMP
             AB%MITYP=MITYP_TMP
          ENDIF
! Deallocate helping arrays
          DEALLOCATE(TYPE_TMP)
          DEALLOCATE(EATOM_TMP)
          DEALLOCATE(NB_TMP)
          DEALLOCATE(NITYP_TMP2)
          DEALLOCATE(POMASS_TMP)
          DEALLOCATE(SZNAM2_HELP2)
        END SUBROUTINE READ_ABCAR_SIZE2

        SUBROUTINE READ_ABCAR_RESET_SIZES(AB,FF,FFM,ML_INPUT_HANDLE,TAG_LIST)
          USE ML_FF_STRUCT, ONLY: ABINITIO, FF_PAR, FFM_PAR, ML_FF_INPUT
          USE ML_FF_TAGLIST
          IMPLICIT NONE
! Input and output variables
          TYPE (ABINITIO)               :: AB
          TYPE (FF_PAR)                 :: FF
          TYPE (FFM_PAR)                :: FFM
          TYPE (ML_FF_INPUT)            :: ML_INPUT_HANDLE
          TYPE (TAGLIST), INTENT(INOUT) :: TAG_LIST
! Local variables
          INTEGER              :: I
          INTEGER              :: ICOL
          INTEGER              :: IERR
          INTEGER              :: IROW
          INTEGER              :: MCONF_HELP
          INTEGER              :: MTYP_TMP
          INTEGER, ALLOCATABLE :: NB_TMP(:) 
          INTEGER              :: NB_LARGEST_TMP
          INTEGER              :: NCOL
          INTEGER              :: NMOD
          INTEGER              :: NWRITE
          NWRITE = 3
          OPEN(UNIT=AB%IU1,FILE=AB%ABCAR,STATUS='OLD',IOSTAT=IERR)
          ! Check if error occured in opening file
          IF (IERR.NE.0) THEN
             CALL ERROR_PRINT2TERMINAL("Error in opening ML_AB file")
          ENDIF
          DO I=1,4
             READ(AB%IU1,*)
          ENDDO
          ! Read number of training configurations
          READ(AB%IU1,*) MCONF_HELP
          DO I=1,3
             READ(AB%IU1,*)
          ENDDO
          READ(AB%IU1,*) MTYP_TMP
          ! Set NCOL and NMOL to read elemeent data.
          NCOL=MTYP_TMP/NWRITE
          NMOD=MOD(MTYP_TMP,NWRITE)
          ! Read AB%ABCAR further
          DO I=1, 3
             READ(AB%IU1,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(AB%IU1,*)
          ENDDO
          IF(NMOD.NE.0) THEN
             READ(AB%IU1,*)
          ENDIF
          DO I=1, 3
             READ(AB%IU1,*)
          ENDDO
          READ(AB%IU1,*)
          DO I=1, 3
             READ(AB%IU1,*)
          ENDDO
          READ(AB%IU1,*)
          DO I=1, 3
             READ(AB%IU1,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(AB%IU1,*)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(AB%IU1,*) 
          ENDIF
          DO I=1, 3
             READ(AB%IU1,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(AB%IU1,*)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(AB%IU1,*)
          ENDIF
          DO I=1, 3
             READ(AB%IU1,*)
          ENDDO
          ! Do allocation of helping array for local configuratins
          ALLOCATE(NB_TMP(1:MTYP_TMP))
          ! Read number of local configurations for each element
          DO ICOL=1, NCOL
             READ(AB%IU1,*) (NB_TMP(NWRITE*(ICOL-1)+IROW), IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(AB%IU1,*) (NB_TMP(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          CLOSE(AB%IU1)
          ! If ML_MB was not set in the INCAR file then 
          ! change the default here
          IF (TAG_LIST%GET_STATE("ML_MB") == TL_DEFAULT) THEN
             CALL TAG_LIST%SET_STATE("ML_MB", TL_OVERRIDE)
             NB_LARGEST_TMP=NB_TMP(1)
             DO I=1, MTYP_TMP
                NB_LARGEST_TMP=MAX(NB_LARGEST_TMP,NB_TMP(I))
             ENDDO
             FFM%MB = NB_LARGEST_TMP+500
             FFM%MB_ALLOCATE=FFM%MB
             FFM%MB_INPUT=FFM%MB
          ENDIF
          ! If ML_MCONF was not set in the INCAR file then
          ! change the default here if needed
          IF (TAG_LIST%GET_STATE("ML_MCONF") == TL_DEFAULT) THEN
             CALL TAG_LIST%SET_STATE("ML_MCONF", TL_OVERRIDE)
             FF%MCONF=MCONF_HELP+500+FF%MCONF_NEW
          ENDIF
          DEALLOCATE(NB_TMP)
        END SUBROUTINE READ_ABCAR_RESET_SIZES
#endif

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Read data in ABCAR file
!****************************************************************************************************

        SUBROUTINE READ_ABCAR (AB,FF,A_AB,ABCAR,B_AB,EATOM,ISTART_FF,IU,ISYS_AB, &
                   LADD_AB,LB_AB,LB_FFM,LBHEAD_FFM, &
                   LNCONF_FF,LNIONS_AB,LNIONS_FF,LNITYP_AB,LNTYP_AB,LSYS_AB, &
                   MTYP_AB,MTYP_FF, &
                   NB_FFM,NBMAX_TMP,NBTOT_FFM,NCONF_AB,NCONF_SYS_AB,NIONS_AB,NITYP_AB,NITYP_ABN, &
                   NSYS_AB,NSYS_FF,NTYP_AB,NTYP_NEW, &
                   OMEGA_AB,POMASS_AB,POMASS_FF,POSION_AB,SZNAM2_AB,SZNAM2_ADD_AB,SZNAM2_FF, &
                   TIFOR_AB,TOTEN_AB,TSIF_AB,TYPE_AB,TYPE_FF)
          USE ML_FF_STRUCT, ONLY: ABINITIO, FF_PAR
          IMPLICIT NONE
! Input variables.
          TYPE (ABINITIO)                  :: AB
          TYPE (FF_PAR)                    :: FF
          CHARACTER(LEN=9) , INTENT(IN)    :: ABCAR
          INTEGER          , INTENT(IN)    :: ISTART_FF
          INTEGER          , INTENT(IN)    :: IU
          INTEGER          , INTENT(IN)    :: NCONF_AB
          INTEGER          , INTENT(IN)    :: MTYP_AB
          INTEGER          , INTENT(IN)    :: MTYP_FF
          INTEGER          , INTENT(IN)    :: NITYP_ABN(:) !(1:MTYP_ABN)
          INTEGER          , INTENT(IN)    :: NSYS_AB
          INTEGER          , INTENT(IN)    :: NSYS_FF
          INTEGER          , INTENT(IN)    :: NTYP_NEW
          REAL(q)          , INTENT(IN)    :: POMASS_FF(:) !(1:MTYP_FF)
          CHARACTER(LEN=40), INTENT(IN)    :: SZNAM2_FF(:) !(1:NSYS_FF)
          CHARACTER(LEN=2) , INTENT(IN)    :: TYPE_FF(:) !(1:MTYP_FF)
! Input and output variables.
          REAL(q)          , INTENT(INOUT) :: A_AB(:,:,:) !(1:3,1:3,1:MCONF_AB)
          REAL(q)          , INTENT(INOUT) :: B_AB(:,:,:) !(1:3,1:3,1:MCONF_AB)
          REAL(q)          , INTENT(INOUT) :: EATOM(:)
          INTEGER          , INTENT(INOUT) :: ISYS_AB(:) !(1:MCONF_AB)
          INTEGER          , INTENT(INOUT) :: LADD_AB(:,:) !(1:MTYP_AB,1:MCONF_AB)
          INTEGER          , INTENT(INOUT) :: LB_AB(:,:) !(1:MIONS_AB,1:MCONF_AB)
          INTEGER          , INTENT(INOUT) :: LB_FFM(:,:) !(1:MB_FFM,1:MTYP_FF)
          INTEGER          , INTENT(INOUT) :: LBHEAD_FFM(:) !(1:MTYP_FF)
          INTEGER          , INTENT(INOUT) :: LNCONF_FF(:,:) !(1:MB_FFM,1:MTYP_FF)
          INTEGER          , INTENT(INOUT) :: LNIONS_FF(:,:) !(1:MB_FFM,1:MTYP_FF)
          INTEGER          , INTENT(INOUT) :: LNIONS_AB(:,:,:) !(1:MITYP_AB,1:MTYP_AB,1:MCONF_AB)
          INTEGER          , INTENT(INOUT) :: LNITYP_AB(:,:) !(1:MIONS_AB,1:MCONF_AB)
          INTEGER          , INTENT(INOUT) :: LNTYP_AB(:,:) !(1:MIONS_AB,1:MCONF_AB)
          INTEGER          , INTENT(INOUT) :: LSYS_AB(:) !(1:MCONF_AB)
          INTEGER          , INTENT(INOUT) :: NB_FFM(:) !(1:MTYP_FF)
          INTEGER          , INTENT(INOUT) :: NBMAX_TMP(:) !(1:MTYP_FF)
          INTEGER          , INTENT(INOUT) :: NBTOT_FFM
          INTEGER          , INTENT(INOUT) :: NCONF_SYS_AB(:) !(1:NSYS_AB)
          INTEGER          , INTENT(INOUT) :: NIONS_AB(:) !(1:MCONF_AB)
          INTEGER          , INTENT(INOUT) :: NITYP_AB(:,:) !(1:MTYP_AB,1:MCONF_AB)
          INTEGER          , INTENT(INOUT) :: NTYP_AB(:) !(1:MCONF_AB)
          REAL(q)          , INTENT(INOUT) :: OMEGA_AB(:) !(1:MCONF_AB)
          REAL(q)          , INTENT(INOUT) :: POMASS_AB(:) !(1:MTYP_AB)
          REAL(q)          , INTENT(INOUT) :: POSION_AB(:,:,:) !(1:3,1:MIONS_AB,1:MCONF_AB)
          CHARACTER(LEN=40), INTENT(INOUT) :: SZNAM2_AB(:) !(1:MCONF_AB)
          CHARACTER(LEN=40), INTENT(INOUT) :: SZNAM2_ADD_AB(:) !(1:NSYS_AB)
          REAL(q)          , INTENT(INOUT) :: TIFOR_AB(:,:,:) !(1:3,1:MIONS_AB,1:MCONF_AB)
          REAL(q)          , INTENT(INOUT) :: TOTEN_AB(:) !(1:MCONF_AB)
          REAL(q)          , INTENT(INOUT) :: TSIF_AB(:,:,:) !(1:3,1:3,1:MCONF_AB)
          CHARACTER(LEN=2) , INTENT(INOUT) :: TYPE_AB(:) !(1:MTYP_AB)
! Local variables.
          REAL(q)                          :: EATOM_TMP(1:MTYP_FF)
          INTEGER                          :: I
          INTEGER                          :: IB
          INTEGER                          :: ICOL
          INTEGER                          :: ICONF
          INTEGER                          :: IERR
          INTEGER                          :: INIONS
          INTEGER                          :: INITYP
          INTEGER                          :: INTYP
          INTEGER                          :: INTYP0
          INTEGER                          :: IROW
          INTEGER                          :: ISYS
          INTEGER                          :: IXYZ
          INTEGER                          :: JNTYP0
          INTEGER                          :: JXYZ
          INTEGER                          :: MIONS_TMP
          INTEGER                          :: MITYP_TMP
          INTEGER                          :: MTYP_TMP
          INTEGER                          :: NB_TMP(1:MTYP_FF)
          INTEGER                          :: NCOL
          INTEGER                          :: NITYP_TMP2(1:MTYP_FF)
          INTEGER                          :: NMOD
          INTEGER                          :: NWRITE
          REAL(q)                          :: POMASS_TMP(1:MTYP_FF)
          CHARACTER(LEN=2)                 :: TYPE_TMP(1:MTYP_FF)
          REAL(q)                          :: VERSION
          CHARACTER(LEN=40)                :: READ_HELPER
          EATOM_TMP = 0.0_q
! Set constants
          NWRITE=3
! Initialize logical parameters
          NBMAX_TMP=0
! Open ABCAR file.
          OPEN(UNIT=IU,FILE=ABCAR,STATUS='OLD',IOSTAT=IERR)
! Check if error occured in opening file
          IF (IERR.NE.0) THEN
             CALL ERROR_PRINT2TERMINAL("Error in opening ML_AB file")
          ENDIF
! Read header
          READ(IU,*) VERSION
          DO I=1, 3
             READ(IU,*)
          ENDDO
          READ(IU,*)
          DO I=1, 3
             READ(IU,*)
          ENDDO
          READ(IU,*) MTYP_TMP
! Set NCOL and NMOL to read elemeent data.
          NCOL=MTYP_TMP/NWRITE
          NMOD=MOD(MTYP_TMP,NWRITE)
! Read ABCAR further
          DO I=1, 3
             READ(IU,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(IU,*) (TYPE_TMP(NWRITE*(ICOL-1)+IROW), IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(IU,*) (TYPE_TMP(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO I=1, 3
             READ(IU,*)
          ENDDO
          READ(IU,*) MIONS_TMP
          DO I=1, 3
             READ(IU,*)
          ENDDO
          READ(IU,*) MITYP_TMP
          DO I=1, 3
             READ(IU,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(IU,*) (EATOM_TMP(NWRITE*(ICOL-1)+IROW), IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(IU,*) (EATOM_TMP(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO I=1, 3
             READ(IU,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(IU,*) (POMASS_TMP(NWRITE*(ICOL-1)+IROW), IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(IU,*) (POMASS_TMP(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
! Configure the model system name.
          DO ISYS=1, NSYS_FF
             SZNAM2_ADD_AB(ISYS)=SZNAM2_FF(ISYS)
          ENDDO
! Set TYPE and POMASS.
          DO INTYP0=1, MTYP_AB
             POMASS_AB(INTYP0)=POMASS_FF(INTYP0)
             TYPE_AB(INTYP0)=TYPE_FF(INTYP0)
          ENDDO
! Read ABCAR file further
          DO I=1, 3
             READ(IU,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(IU,*) (NB_TMP(NWRITE*(ICOL-1)+IROW), IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(IU,*) (NB_TMP(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          NBTOT_FFM=0
          DO INTYP=1, MTYP_TMP
             DO JNTYP0=1, MTYP_FF
                IF(TYPE_TMP(INTYP).EQ.TYPE_FF(JNTYP0)) THEN
                   NB_FFM(JNTYP0)=NB_TMP(INTYP)
                   NBTOT_FFM=NBTOT_FFM+NB_TMP(INTYP)
                   DO I=1,  3
                      READ(IU,*)
                   ENDDO
                   DO IB=1, NB_TMP(INTYP)
                      READ(IU,*) LNCONF_FF(IB,JNTYP0),LNIONS_FF(IB,JNTYP0)
                      LB_AB(LNIONS_FF(IB,JNTYP0),LNCONF_FF(IB,JNTYP0))=IB
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
          LBHEAD_FFM(1)=0
          DO INTYP0=2, MTYP_FF
             LBHEAD_FFM(INTYP0)=LBHEAD_FFM(INTYP0-1)+NB_FFM(INTYP0-1)
          ENDDO
          NBTOT_FFM=0
          DO INTYP0=1, MTYP_FF
             DO IB=1, NB_FFM(INTYP0)
                NBTOT_FFM=NBTOT_FFM+1
                LB_FFM(IB,INTYP0)=NBTOT_FFM
             ENDDO
          ENDDO
          DO ICONF=1, NCONF_AB
             DO I=1, 5
                READ(IU,*)
             ENDDO
             READ(IU,'(A40)') SZNAM2_AB(ICONF)
             SZNAM2_AB(ICONF)=TRIM(ADJUSTL(SZNAM2_AB(ICONF)))
             DO ISYS=1, NSYS_AB
                IF(SZNAM2_AB(ICONF).EQ.SZNAM2_ADD_AB(ISYS)) THEN
                   ISYS_AB(ICONF)=ISYS
                   EXIT
                ENDIF
             ENDDO
             DO I=1, 3
                READ(IU,*)
             ENDDO
             READ(IU,*) NTYP_AB(ICONF)
             DO I=1, 3
                READ(IU,*)
             ENDDO
             READ(IU,*) NIONS_AB(ICONF)
             DO I=1, 3
                READ(IU,*)
             ENDDO
             DO INTYP=1, NTYP_AB(ICONF)
                READ(IU,*) TYPE_TMP(INTYP),NITYP_TMP2(INTYP)
                DO JNTYP0=1, MTYP_FF
                   IF(TYPE_TMP(INTYP) .EQ. TYPE_FF(JNTYP0)) THEN
                      LADD_AB(INTYP,ICONF)=JNTYP0
                      NITYP_AB(JNTYP0,ICONF)=NITYP_TMP2(INTYP)
                      NBMAX_TMP(JNTYP0)=NBMAX_TMP(JNTYP0)+NITYP_TMP2(INTYP)
                      EXIT
                   ENDIF
                ENDDO
             ENDDO
             INIONS=0
             DO INTYP=1, NTYP_AB(ICONF)
                INTYP0=LADD_AB(INTYP,ICONF)
                DO INITYP=1, NITYP_AB(INTYP0,ICONF)
                   INIONS=INIONS+1
                   LNITYP_AB(INIONS,ICONF)=INITYP
                   LNIONS_AB(INITYP,INTYP0,ICONF)=INIONS
                   LNTYP_AB(INIONS,ICONF)=INTYP0
                ENDDO
             ENDDO
             ! Optional reading of CTIFOR
             READ(IU,*)
             READ(IU,*) READ_HELPER
             READ(IU,*)
             READ_HELPER=ADJUSTL(TRIM(READ_HELPER))
             IF (READ_HELPER(1:1).EQ.'C') THEN
                READ(IU,*) AB%CTIFOR_ALLCONF_OLD(ICONF)
                AB%CTIFOR_ALLCONF_OLD(ICONF)=AB%CTIFOR_ALLCONF_OLD(ICONF)/FUNIT
                DO I=1, 3
                   READ(IU,*)
                ENDDO
                IF (ICONF.GT.1.AND.(.NOT.FF%LCTIFOR_READ_FROM_FILE)) THEN
                   CALL ERROR_PRINT2TERMINAL("Error CTIFOR in ML_AB")
                ENDIF
                FF%LCTIFOR_READ_FROM_FILE=.TRUE.
             ELSE
                IF (ICONF.GT.1.AND.FF%LCTIFOR_READ_FROM_FILE) THEN
                   CALL ERROR_PRINT2TERMINAL("Error CTIFOR in ML_AB")
                ENDIF
                FF%LCTIFOR_READ_FROM_FILE=.FALSE.
             ENDIF
             ! End optional reading of CTIFOR
             DO IXYZ=1, 3
                READ(IU,*) (A_AB(JXYZ,IXYZ,ICONF),JXYZ=1, 3)
             ENDDO
             DO IXYZ=1, 3
                DO JXYZ=1, 3
                   A_AB(JXYZ,IXYZ,ICONF)=A_AB(JXYZ,IXYZ,ICONF)/AUTOA
                ENDDO
             ENDDO
             CALL RECIPROCAL_LATTICE (A_AB(1:,1:,ICONF),B_AB(1:,1:,ICONF),OMEGA_AB(ICONF))
             DO I=1, 3
                READ(IU,*)
             ENDDO
             DO INIONS=1, NIONS_AB(ICONF)
                READ(IU,*) (POSION_AB(IXYZ,INIONS,ICONF),IXYZ=1, 3)
                DO IXYZ=1, 3
                   POSION_AB(IXYZ,INIONS,ICONF)=POSION_AB(IXYZ,INIONS,ICONF)/AUTOA
                ENDDO
             ENDDO
             DO I=1, 3
                READ(IU,*)
             ENDDO
             READ(IU,*) TOTEN_AB(ICONF)
             TOTEN_AB(ICONF)=TOTEN_AB(ICONF)/EUNIT
             DO I=1, 3
                READ(IU,*)
             ENDDO
             DO INIONS=1, NIONS_AB(ICONF)
                READ(IU,*) (TIFOR_AB(IXYZ,INIONS,ICONF),IXYZ=1, 3)
                DO IXYZ=1, 3
                   TIFOR_AB(IXYZ,INIONS,ICONF)=TIFOR_AB(IXYZ,INIONS,ICONF)/FUNIT
                ENDDO
             ENDDO
             DO I=1, 5
                READ(IU,*)
             ENDDO
             READ(IU,*) TSIF_AB(1,1,ICONF),TSIF_AB(2,2,ICONF),TSIF_AB(3,3,ICONF)
             DO I=1, 3
                READ(IU,*)
             ENDDO
             READ(IU,*) TSIF_AB(2,1,ICONF),TSIF_AB(3,2,ICONF),TSIF_AB(3,1,ICONF)
             DO IXYZ=1, 2
                DO JXYZ=IXYZ+1, 3
                   TSIF_AB(IXYZ,JXYZ,ICONF)=TSIF_AB(JXYZ,IXYZ,ICONF)
                ENDDO
             ENDDO
             DO IXYZ=1, 3
                DO JXYZ=1, 3
                   TSIF_AB(JXYZ,IXYZ,ICONF)=TSIF_AB(JXYZ,IXYZ,ICONF)/SUNIT
                ENDDO
             ENDDO
             IF (VERSION.EQ.(0.0_q)) THEN
                DO I=1, 3
                   READ(IU,*)
                ENDDO
                DO INIONS=1, NIONS_AB(ICONF)
                   READ(IU,*)
                ENDDO
             ENDIF
          ENDDO
          CLOSE(IU)
! If CTIFORs were not read from ML_AB file then no previous knowledge on 
! CTIFORs is used
          IF (.NOT.FF%LCTIFOR_READ_FROM_FILE) THEN
             FF%LCTIFOR_USE_FROM_FILE=.FALSE.
          ENDIF
! Check whether the system stored in the database is the same as the current system
! This is only relevant for ISTART_FF.NE.3
          IF (ISTART_FF.NE.3) THEN
             DO ICONF=1, NCONF_AB
                LSYS_AB(ICONF)=1
                IF(NTYP_AB(ICONF).NE.NTYP_NEW) THEN
                   LSYS_AB(ICONF)=0
                ENDIF
                IF(LSYS_AB(ICONF).EQ.1) THEN
                   DO INTYP=1, NTYP_AB(ICONF)
                      INTYP0=LADD_AB(INTYP,ICONF)
                      IF(NITYP_AB(INTYP0,ICONF).NE.NITYP_ABN(INTYP0)) THEN
                         LSYS_AB(ICONF)=0
                         EXIT
                      ENDIF
                   ENDDO
                ENDIF
             ENDDO
          ENDIF
! Check the number of configuration per model system.
          NCONF_SYS_AB=0
          DO ICONF=1, NCONF_AB
             NCONF_SYS_AB(ISYS_AB(ICONF))=NCONF_SYS_AB(ISYS_AB(ICONF))+1
          ENDDO
! Set the atomic reference energies read from ML_AB file
          DO INTYP=1, MTYP_FF
             EATOM(INTYP)=EATOM_TMP(INTYP)/EUNIT
          ENDDO
        END SUBROUTINE READ_ABCAR
#endif

!****************************************************************************************************
! Read FFCAR file and overwrite the FF parameters by the parameters in this file.
!****************************************************************************************************

        SUBROUTINE OVERWRITE_FF_PAR_BY_PAR_IN_FFCAR (TAG_LIST,FF,FFM)
          USE ML_FF_STRUCT, ONLY : FF_PAR, FFM_PAR, ML_IO_WRITE
          USE ML_FF_TAGLIST
          IMPLICIT NONE
! Input and output variables.
          TYPE (TAGLIST), INTENT(INOUT) :: TAG_LIST
          TYPE (FF_PAR)                 :: FF
          TYPE (FFM_PAR)                :: FFM
! Local variables.
          INTEGER                       :: I
          INTEGER                       :: IB
          INTEGER                       :: ICOL
          INTEGER                       :: IERR
          INTEGER                       :: INTYP
          INTEGER                       :: IROW
          INTEGER                       :: JNTYP
          INTEGER                       :: MTYP_LOC
          INTEGER                       :: NB_FFM(1:FF%MTYP)
          INTEGER                       :: NBTOT_FFM
          INTEGER                       :: NCOL
          INTEGER                       :: NMOD
          INTEGER, PARAMETER            :: NWRITE = 3
          CHARACTER(LEN=2)              :: TYPE_LOC(1:FF%MTYP)
          INTEGER                       :: IDUMMY(1:NWRITE)
          REAL(q)                       :: RDUMMY(1:NWRITE)
! Open FFCAR file.
!          OPEN(FF%IU1,FILE=FF%FFCAR,STATUS='OLD',IOSTAT=IERR)
          OPEN(FF%IU1,FILE=FF%FFCAR,STATUS='OLD',FORM='UNFORMATTED',ACCESS='STREAM',IOSTAT=IERR)
! Check if error occured in opening file
          IF (IERR.NE.0) THEN
             CALL ERROR_PRINT2TERMINAL("Error in opening ML_FF file")
          ENDIF
! Read the parameters.
          ! Read the version number
          READ(FF%IU1) FF%ML_FF_VERSION_READ(1), FF%ML_FF_VERSION_READ(2), FF%ML_FF_VERSION_READ(3)
          ! The number of element type
          READ(FF%IU1) MTYP_LOC
          IF(FF%MTYP.NE.MTYP_LOC) THEN
             CALL ERROR_PRINT2TERMINAL ('!!! MLFF : you need to generate a force field before adapting ISTART_FF=2 !!!')
          ENDIF
          NCOL=FF%MTYP/NWRITE
          NMOD=MOD(FF%MTYP,NWRITE)
          ! The element type
          DO ICOL=1, NCOL, 1
             READ(FF%IU1) (TYPE_LOC(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1) (TYPE_LOC(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO I=1, FF%MTYP
             IF(TYPE_LOC(I).NE.FF%TYPE(I)) THEN
                CALL ERROR_PRINT2TERMINAL ('!!! MLFF : you need to generate force field before adapting ISTART_FF=2 !!!')
             ENDIF
          ENDDO
          ! Reference atomic energies (eV)
          DO ICOL=1, NCOL, 1
             READ(FF%IU1) (FF%EATOM_REF(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1) (FF%EATOM_REF(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO INTYP=1, MTYP_LOC
             FF%EATOM_REF(INTYP)=FF%EATOM_REF(INTYP)/EUNIT
          ENDDO
          ! Atomic mass
          DO ICOL=1, NCOL, 1
             READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
          ENDIF
          ! Type of weighting
          READ(FF%IU1) FF%IWEIGHT
          ! Weights
          READ(FF%IU1) FF%WTOTEN,FF%WTIFOR,FF%WTSIF
          ! Mean square root error in energy (eV atom^-1)
          READ(FF%IU1) FF%STOTEN
          FF%STOTEN=FF%STOTEN/EUNIT
          ! Mean square root error in force (eV Angst^-1)
          READ(FF%IU1) FF%STIFOR
          FF%STIFOR=FF%STIFOR/FUNIT
          ! Mean square root error in stress tensor (kB)
          READ(FF%IU1) FF%STSIF
          FF%STSIF=FF%STSIF/SUNIT
! Read many-body parameters
          ! Many-body: true or false? - always true
          READ(FF%IU1) FF%LMLMB
          IF(FF%LMLMB) THEN
             ! Self-interaction correction?
             READ(FF%IU1) FFM%LSIC
             ! Type of similarity measure
             READ(FF%IU1) FFM%LSUPERVEC
             ! Weights for radial and angular descriptors
             READ(FF%IU1) FFM%W1,FFM%W2
             ! Type of cutoff function
             READ(FF%IU1) FFM%ICUT1,FFM%ICUT2
             ! Cutoff radius (Angst)
             READ(FF%IU1) FFM%RCUT1,FFM%RCUT2
             FFM%RCUT1=FFM%RCUT1/AUTOA
             FFM%RCUT2=FFM%RCUT2/AUTOA
             ! Broadening type
             READ(FF%IU1) FFM%IBROAD1,FFM%IBROAD2
             ! Broadening parameter
             READ(FF%IU1) FFM%SION1,FFM%SION2
             FFM%SION1=FFM%SION1/AUTOA
             FFM%SION2=FFM%SION2/AUTOA
             ! The maximum number of radial basis set
             READ(FF%IU1) FFM%MRB1,FFM%MRB2
             ! The number of grids for integrations
             READ(FF%IU1) FFM%NR1,FFM%NR2
             ! The number of spline grids 
             READ(FF%IU1) FFM%MSPL1,FFM%MSPL2
             FFM%NSPL1=FFM%MSPL1
             FFM%NSPL2=FFM%MSPL2
             ! Lmax for spherical harmonics
             READ(FF%IU1) FFM%LMAX1,FFM%LMAX2
             IF(FFM%W1.GT.0.0_q) THEN
                FFM%MMVAR1=FFM%MRB1*MTYP_LOC
                FFM%NNVAR1=FFM%MMVAR1
             ENDIF
             ! Polynomial parameters for SOAP
             READ(FF%IU1) FFM%NHYP1,FFM%NHYP2
             ! Normalization
             READ(FF%IU1) FFM%LNORM1,FFM%LNORM2
             ! Radial window
             READ(FF%IU1) FFM%LWINDOW1,FFM%LWINDOW2
             IF(FFM%LWINDOW1) THEN
                ! Radial descriptor window
                READ(FF%IU1) FFM%IWINDOW1
             ENDIF
             IF(FFM%LWINDOW2) THEN
                ! Angular descriptor window
                READ(FF%IU1) FFM%IWINDOW2
             ENDIF
             ! Angular filtering
             READ(FF%IU1) FFM%LAFILT2
             IF(FFM%LAFILT2) THEN
                ! Type of angular filtering
                READ(FF%IU1) FFM%IAFILT2
                IF(FFM%IAFILT2.EQ.2) THEN
                   ! Damping parameter
                   READ(FF%IU1) FFM%AFILT2
                ENDIF
             ENDIF
             ! Radial metric
             READ(FF%IU1) FFM%LMETRIC1,FFM%LMETRIC2
             IF(FFM%LMETRIC1) THEN
                ! Radial descriptor metric
                READ(FF%IU1) FFM%NMETRIC1,FFM%RMETRIC1
             ENDIF
             IF(FFM%LMETRIC2) THEN
                ! Angular descriptor metric
                READ(FF%IU1) FFM%NMETRIC2,FFM%RMETRIC2
             ENDIF
             ! Variable transformation
             READ(FF%IU1) FFM%LVARTRAN1,FFM%LVARTRAN2
             IF(FFM%LVARTRAN1) THEN
                ! Radial variable transformation
                READ(FF%IU1) FFM%NVARTRAN1
             ENDIF
             IF(FFM%LVARTRAN2) THEN
                ! Angular variable transformation
                READ(FF%IU1) FFM%NVARTRAN2
             ENDIF
             IF(FFM%W2.GT.0.0_q) THEN
! The number of radial basis sets for each angular component.
                NCOL=(FFM%LMAX2+1)/NWRITE
                NMOD=MOD(FFM%LMAX2+1,NWRITE)
                ! The number of radial basis sets for each angular channel
                DO ICOL=1, NCOL
                   READ(FF%IU1) (IDUMMY(IROW),IROW=1, NWRITE)
                ENDDO
                IF(NMOD .NE. 0) THEN
                   READ(FF%IU1) (IDUMMY(IROW),IROW=1, NMOD)
                ENDIF
! Lists of active descriptors.
                DO INTYP=1, MTYP_LOC
                   ! The number of active descriptors for XXX
                   READ(FF%IU1) FFM%NNVAR2(INTYP)
                   NCOL=(FFM%NNVAR2(INTYP))/NWRITE
                   NMOD=MOD(FFM%NNVAR2(INTYP),NWRITE)
                   ! List of active descriptors
                   DO ICOL=1, NCOL
                      READ(FF%IU1) (IDUMMY(IROW),IROW=1, NWRITE)
                   ENDDO
                   IF(NMOD .NE. 0) THEN
                      READ(FF%IU1) (IDUMMY(IROW),IROW=1, NMOD)
                   ENDIF
                ENDDO
                FFM%MMVAR2=FFM%NNVAR2(1)
                DO INTYP=2, MTYP_LOC, 1
                   FFM%MMVAR2=MAX(FFM%MMVAR2,FFM%NNVAR2(INTYP))
                ENDDO
! Lists of active descriptors for self-interaction correction if necessary.
                IF(FFM%LSIC) THEN
                   DO INTYP=1, MTYP_LOC
                      ! The number of active SIC descriptors for XXX 
                      READ(FF%IU1) FFM%NNVAR2_SIC(INTYP)
                      NCOL=(FFM%NNVAR2_SIC(INTYP))/NWRITE
                      NMOD=MOD(FFM%NNVAR2_SIC(INTYP),NWRITE)
                      ! List of active descriptors for SIC
                      DO ICOL=1, NCOL
                         READ(FF%IU1) (IDUMMY(IROW),IROW=1, NWRITE)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1) (IDUMMY(IROW),IROW=1, NMOD)
                      ENDIF
                   ENDDO
                   FFM%MMVAR2_SIC=FFM%NNVAR2_SIC(1)
                   DO INTYP=2, MTYP_LOC, 1
                      FFM%MMVAR2_SIC=MAX(FFM%MMVAR2_SIC,FFM%NNVAR2_SIC(INTYP))
                   ENDDO
                ENDIF
             ENDIF
             NCOL=MTYP_LOC/NWRITE
             NMOD=MOD(MTYP_LOC,NWRITE)
             ! The numbers of basis sets
             DO ICOL=1, NCOL
                READ(FF%IU1) (NB_FFM(NWRITE*(ICOL-1)+IROW), IROW=1, NWRITE)
             ENDDO
             IF(NMOD .NE. 0) THEN
                READ(FF%IU1) (NB_FFM(NWRITE*NCOL+IROW),IROW=1, NMOD)
             ENDIF
! Set the required size automatically to the largest required if
! ML_ISTART=2
! If ML_MB was not read in then set to the first to
! first entry in ML_FF file and then to the largest 
             IF (TAG_LIST%GET_STATE("ML_MB") == TL_DEFAULT) THEN
                CALL TAG_LIST%SET_STATE("ML_MB", TL_OVERRIDE)
                FFM%MB=NB_FFM(1)+10
             ENDIF
             DO INTYP=1, MTYP_LOC
                IF(NB_FFM(INTYP) .GT. FFM%MB) THEN
                    FFM%MB=NB_FFM(INTYP)+10
                ENDIF
             ENDDO
! Further reading
             ! Type of the method to scale the energy data
             READ(FF%IU1) FFM%ISCALE_TOTEN
             ! Average energy per atom (eV atom^-1)
             READ(FF%IU1) FFM%TOTENAV
             FFM%TOTENAV=FFM%TOTENAV/EUNIT
             DO INTYP=1, MTYP_LOC
                NCOL=NB_FFM(INTYP)/NWRITE
                NMOD=MOD(NB_FFM(INTYP),NWRITE)
                ! Regression coefficients (eV atom^-1)
                ! For XXX
                DO ICOL=1, NCOL
                   READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
                ENDDO
                IF(NMOD .NE. 0) THEN
                   READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
                ENDIF
             ENDDO
             IF(FFM%W1 .GT. 0.0_q) THEN
                NCOL=FFM%MMVAR1/NWRITE
                NMOD=MOD(FFM%MMVAR1,NWRITE)
                DO INTYP=1, MTYP_LOC
                   ! C00 (radial descriptor)
                   ! For XXX
                   DO IB=1, NB_FFM(INTYP)
                      DO ICOL=1, NCOL
                         READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF
             IF(FFM%W2 .GT. 0.0_q) THEN
                DO INTYP=1, MTYP_LOC
                   NCOL=FFM%NNVAR2(INTYP)/NWRITE
                   NMOD=MOD(FFM%NNVAR2(INTYP),NWRITE)
                   ! Cnlm (angular descriptor) (-)
                   ! For XXX
                   DO IB=1, NB_FFM(INTYP)
                      DO ICOL=1, NCOL
                         READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF
             ! Parameters for Bayesian error estimation
             ! Normalized noise parameter (-)     
             READ(FF%IU1) FFM%SIGV
             ! Normalized uncertainity parameter (-)        
             READ(FF%IU1) FFM%SIGW
             ! Variance in energies per atom in training data (eV atom^-1)     
             READ(FF%IU1) RDUMMY(1)
             ! Variance in forces in training data (eV Angst^-1)
             READ(FF%IU1) RDUMMY(1),RDUMMY(2),RDUMMY(3)
             ! Variance in stress tensors in training data (kB)      
             READ(FF%IU1) RDUMMY(1),RDUMMY(2),RDUMMY(3)
             READ(FF%IU1) RDUMMY(1),RDUMMY(2),RDUMMY(3)
             NBTOT_FFM=0
             DO INTYP=1, MTYP_LOC
                NBTOT_FFM=NBTOT_FFM+NB_FFM(INTYP)
             ENDDO
             NCOL=NBTOT_FFM/NWRITE
             NMOD=MOD(NBTOT_FFM,NWRITE)
             ! Covariance matrix (-)
             DO IB=1, NBTOT_FFM
                DO ICOL=1, NCOL
                   READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
                ENDDO
                IF(NMOD.NE.0) THEN
                   READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
                ENDIF
             ENDDO
          ENDIF
! Read charge parameters
! Keep this in for the moment to be compatible with old format (for LMLPC)
!          DO I=1, 3
!             READ(FF%IU1)
!          ENDDO
!          READ(FF%IU1) 
! Close FFCAR file.
          CLOSE(FF%IU1)
        END SUBROUTINE OVERWRITE_FF_PAR_BY_PAR_IN_FFCAR


!****************************************************************************************************
! Read FFCAR file.
!****************************************************************************************************

        SUBROUTINE READ_FFCAR (FF,FFM,PAR_SUP_HANDLE,ISYS_NEW)
          USE ML_FF_STRUCT, ONLY : FF_PAR, FFM_PAR, PARALLEL_SUPER
          IMPLICIT NONE
! Input variables.
          TYPE (FF_PAR)                   :: FF
          TYPE (FFM_PAR)                  :: FFM
          TYPE (PARALLEL_SUPER)           :: PAR_SUP_HANDLE
          INTEGER         , INTENT(IN)    :: ISYS_NEW
! Local variables.
          INTEGER                         :: I
          INTEGER                         :: IB
          INTEGER                         :: ICOL
          INTEGER                         :: IERR
          INTEGER                         :: INTYP
          INTEGER                         :: INTYP0
          INTEGER                         :: IRB
          INTEGER                         :: IROW
          INTEGER                         :: IVAR
          INTEGER                         :: JNTYP
          INTEGER                         :: JNTYP0
          INTEGER                         :: JJNTYP0
          INTEGER                         :: JRB
          INTEGER                         :: L
          INTEGER                         :: LVAR_HELP(1:FFM%MMVAR2,1:FF%MTYP)
          INTEGER                         :: LVAR_SIC_HELP(1:FFM%MMVAR2_SIC,1:FF%MTYP)
          INTEGER                         :: MTYP_LOC
          INTEGER                         :: NCOL
          INTEGER                         :: NMOD
          INTEGER                         :: NNVAR2_HELP1
          INTEGER                         :: NNVAR2_HELP2
          INTEGER                         :: NNVAR2_SIC_HELP1
          INTEGER                         :: NNVAR2_SIC_HELP2
          INTEGER                         :: NRB2(0:FFM%LMAX2)
          INTEGER                         :: NROW
          INTEGER, PARAMETER              :: NWRITE = 3
          CHARACTER(LEN=2)                :: TYPE_LOC(1:FF%MTYP)
          INTEGER                         :: IDUMMY(1:NWRITE)
          REAL(q)                         :: RDUMMY(1:NWRITE)
! Start profiling
          PROFILING_START('read_ffcar')
! Open FF%FFCAR file.
          OPEN(FF%IU1,FILE=FF%FFCAR,STATUS='OLD',FORM='UNFORMATTED',ACCESS='STREAM',IOSTAT=IERR)
! Check if error occured in opening file
          IF (IERR.NE.0) THEN
             CALL ERROR_PRINT2TERMINAL("Error in opening ML_FF file")
          ENDIF
! Read the parameters.
          ! Read the version number
          READ(FF%IU1) FF%ML_FF_VERSION_READ(1), FF%ML_FF_VERSION_READ(2), FF%ML_FF_VERSION_READ(3)
          ! The number of element type
          READ(FF%IU1) MTYP_LOC
          NCOL=MTYP_LOC/NWRITE
          NMOD=MOD(MTYP_LOC,NWRITE)
          ! The element type
          DO ICOL=1, NCOL, 1
             READ(FF%IU1) (TYPE_LOC(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1) (TYPE_LOC(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          ! Reference atomic energies (eV)
          DO ICOL=1, NCOL, 1
             READ(FF%IU1) (FF%EATOM_REF(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1) (FF%EATOM_REF(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO INTYP=1, MTYP_LOC
             FF%EATOM_REF(INTYP)=FF%EATOM_REF(INTYP)/EUNIT
          ENDDO
          ! Atomic mass
          DO ICOL=1, NCOL, 1
             READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE) 
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
          ENDIF
          ! Type of weighting
          READ(FF%IU1) FF%IWEIGHT
          ! Weights
          READ(FF%IU1) FF%WTOTEN,FF%WTIFOR,FF%WTSIF
          IF(FF%IWEIGHT.EQ.1) THEN
             FF%WTIFOR=FF%WTIFOR/FUNIT
             FF%WTOTEN=FF%WTOTEN/EUNIT
             FF%WTSIF=FF%WTSIF/SUNIT
          ENDIF
          ! Mean square root error in energy (eV atom^-1)
          READ(FF%IU1) FF%STOTEN
          FF%STOTEN=FF%STOTEN/EUNIT
          ! Mean square root error in force (eV Angst^-1)
          READ(FF%IU1) FF%STIFOR
          FF%STIFOR=FF%STIFOR/FUNIT
          ! Mean square root error in stress tensor (kB)
          READ(FF%IU1) FF%STSIF
          FF%STSIF=FF%STSIF/SUNIT
! Read many-body parameters
          ! Many-body: true or false?
          READ(FF%IU1) FF%LMLMB
          IF(FF%LMLMB) THEN
             ! Self-interaction correction?
             READ(FF%IU1) FFM%LSIC
             ! Type of similarity measure
             READ(FF%IU1) FFM%LSUPERVEC
             ! Weights for radial and angular descriptors
             READ(FF%IU1) FFM%W1,FFM%W2
             ! Type of cutoff function
             READ(FF%IU1) FFM%ICUT1,FFM%ICUT2
             ! Cutoff radius (Angst)
             READ(FF%IU1) FFM%RCUT1,FFM%RCUT2
             FFM%RCUT1=FFM%RCUT1/AUTOA
             FFM%RCUT2=FFM%RCUT2/AUTOA
             ! Broadening type
             READ(FF%IU1) FFM%IBROAD1,FFM%IBROAD2
             ! Broadening parameter
             READ(FF%IU1) FFM%SION1,FFM%SION2
             FFM%SION1=FFM%SION1/AUTOA
             FFM%SION2=FFM%SION2/AUTOA
            ! The maximum number of radial basis set
             READ(FF%IU1) IDUMMY(1),IDUMMY(2)
            ! The number of grids for integrations
             READ(FF%IU1) FFM%NR1,FFM%NR2
            ! The number of spline grids
             READ(FF%IU1) FFM%MSPL1,FFM%MSPL2
            ! Lmax for spherical harmonics
             READ(FF%IU1) IDUMMY(1),IDUMMY(2)
            ! Polynomial parameters for SOAP
             READ(FF%IU1) FFM%NHYP1,FFM%NHYP2
            ! Normalization 
             READ(FF%IU1) FFM%LNORM1,FFM%LNORM2
            ! Radial window
             READ(FF%IU1) FFM%LWINDOW1,FFM%LWINDOW2
             IF(FFM%LWINDOW1) THEN
               READ(FF%IU1) FFM%IWINDOW1
             ENDIF
             IF(FFM%LWINDOW2) THEN
               READ(FF%IU1) FFM%IWINDOW2
             ENDIF
            ! Angular filtering
             READ(FF%IU1) FFM%LAFILT2
             IF(FFM%LAFILT2) THEN
               ! Type of angular filtering
                READ(FF%IU1) FFM%IAFILT2
                IF(FFM%IAFILT2.EQ.2) THEN
                  ! Damping parameter
                   READ(FF%IU1) FFM%AFILT2
                ENDIF
             ENDIF
            ! Radial metric
             READ(FF%IU1) FFM%LMETRIC1,FFM%LMETRIC2
             IF(FFM%LMETRIC1) THEN
               READ(FF%IU1) FFM%NMETRIC1,FFM%RMETRIC1
             ENDIF
             IF(FFM%LMETRIC2) THEN
               READ(FF%IU1) FFM%NMETRIC2,FFM%RMETRIC2
             ENDIF
            ! Variable transformation
             READ(FF%IU1) FFM%LVARTRAN1,FFM%LVARTRAN2
             IF(FFM%LVARTRAN1) THEN
               READ(FF%IU1) FFM%NVARTRAN1
             ENDIF
             IF(FFM%LVARTRAN2) THEN
               READ(FF%IU1) FFM%NVARTRAN2
             ENDIF
             IF(FFM%W2.GT.0.0_q) THEN
! Read the number of radial basis sets for each angular components.
                NCOL=(FFM%LMAX2+1)/NWRITE
                NMOD=MOD(FFM%LMAX2+1,NWRITE)
               ! The number of radial basis sets for each angular channel
                DO ICOL=1, NCOL, 1
                   READ(FF%IU1) (NRB2(NWRITE*(ICOL-1)+IROW),IROW=0, NWRITE-1)
                ENDDO
                IF(NMOD .NE. 0) THEN
                   READ(FF%IU1) (NRB2(NWRITE*NCOL+IROW),IROW=0, NMOD-1)
                ENDIF
! Read the list of active descriptors.
                DO INTYP=1, MTYP_LOC
                  !  The number of active descriptors for XXX
                   READ(FF%IU1) FFM%NNVAR2(INTYP)
                   NCOL=(FFM%NNVAR2(INTYP))/NWRITE
                   NMOD=MOD(FFM%NNVAR2(INTYP),NWRITE)
                  ! List of active descriptors
                   DO ICOL=1, NCOL, 1
                      READ(FF%IU1) (LVAR_HELP(NWRITE*(ICOL-1)+IROW,INTYP),IROW=1, NWRITE)
                   ENDDO
                   IF(NMOD .NE. 0) THEN
                      READ(FF%IU1) (LVAR_HELP(NWRITE*NCOL+IROW,INTYP),IROW=1, NMOD)
                   ENDIF
                ENDDO
! Read the list of active descriptors for SIC.
                IF(FFM%LSIC) THEN
                   DO INTYP=1, MTYP_LOC
                     ! The number of active SIC descriptors for XXX
                      READ(FF%IU1) FFM%NNVAR2_SIC(INTYP)
                      NCOL=(FFM%NNVAR2_SIC(INTYP))/NWRITE
                      NMOD=MOD(FFM%NNVAR2_SIC(INTYP),NWRITE)
                     ! List of active descriptors for SIC
                      DO ICOL=1, NCOL, 1
                         READ(FF%IU1) (LVAR_SIC_HELP(NWRITE*(ICOL-1)+IROW,INTYP),IROW=1, NWRITE)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1) (LVAR_SIC_HELP(NWRITE*NCOL+IROW,INTYP),IROW=1, NMOD)
                      ENDIF
                   ENDDO
                ENDIF
! Set lists.
! Determination of the numbers of variables related to the angular descriptor.
                FFM%LVAR=0
                FFM%LFLAG_VAR=.FALSE.
                FFM%WVAR=0.0_q
                IF(FFM%LSIC) THEN
                   FFM%LVAR_SIC=0
                   FFM%LFLAG_VAR_SIC=.FALSE.
                ENDIF
                DO INTYP0=1, FF%MTYP
                   NNVAR2_HELP1=0
                   NNVAR2_HELP2=0
                   IF(FFM%LSIC) THEN
                      NNVAR2_SIC_HELP1=0
                      NNVAR2_SIC_HELP2=0
                   ENDIF
                   DO JNTYP0=1, FF%MTYP
                      DO JJNTYP0=1, FF%MTYP
                         DO L=0, FFM%LMAX2
                            DO IRB=1, NRB2(L)
                               DO JRB=IRB, NRB2(L)
                                  NNVAR2_HELP1=NNVAR2_HELP1+1
                                  DO IVAR=1, FFM%NNVAR2(INTYP0)
                                     IF(NNVAR2_HELP1.EQ.LVAR_HELP(IVAR,INTYP0)) THEN
                                        NNVAR2_HELP2=NNVAR2_HELP2+1
                                        FFM%LVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=NNVAR2_HELP2
                                        FFM%LFLAG_VAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=.TRUE.
                                        IF(JRB.EQ.IRB) THEN
                                           FFM%WVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=1.0_q
                                        ELSE
                                           FFM%WVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=SQRT(2.0_q)
                                        ENDIF
                                        EXIT
                                     ENDIF
                                  ENDDO
                                  IF(FFM%LSIC) THEN
                                     IF(JJNTYP0.EQ.JNTYP0) THEN
                                        NNVAR2_SIC_HELP1=NNVAR2_SIC_HELP1+1
                                        DO IVAR=1, FFM%NNVAR2_SIC(INTYP0)
                                           IF(NNVAR2_SIC_HELP1.EQ.LVAR_SIC_HELP(IVAR,INTYP0)) THEN
                                              NNVAR2_SIC_HELP2=NNVAR2_SIC_HELP2+1
                                              FFM%LVAR_SIC(JRB,IRB,L,JNTYP0,INTYP0)=NNVAR2_SIC_HELP2
                                              FFM%LFLAG_VAR_SIC(JRB,IRB,L,JNTYP0,INTYP0)=.TRUE.
                                              EXIT
                                           ENDIF
                                        ENDDO
                                     ENDIF
                                  ENDIF
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
! Read basis sets.
             NCOL=MTYP_LOC/NWRITE
             NMOD=MOD(MTYP_LOC,NWRITE)
            ! The numbers of basis sets
             DO ICOL=1, NCOL, 1
                READ(FF%IU1) (FFM%NB(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
             ENDDO
             IF(NMOD .NE. 0) THEN
                READ(FF%IU1) (FFM%NB(NWRITE*NCOL+IROW),IROW=1, NMOD)
             ENDIF
            ! Type of the method to scale the energy data            
             READ(FF%IU1) FFM%ISCALE_TOTEN
            ! Average energy per atom (eV atom^-1)
             READ(FF%IU1) FFM%TOTENAV
             FFM%TOTENAV=FFM%TOTENAV/EUNIT
            DO INTYP=1, MTYP_LOC
                NCOL=FFM%NB(INTYP)/NWRITE
                NMOD=MOD(FFM%NB(INTYP),NWRITE)
               ! Regression coefficients (eV atom^-1)
                ! For XXX
#ifdef scaLAPACK
                CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES_ALLOCATE(PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_WMAT_FFM)
                IF (PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_ROW.EQ.0.AND.PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_COL.EQ.0) THEN
                   DO ICOL=1, NCOL
                      READ(FF%IU1) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*(ICOL-1)+IROW,1),IROW=1, NWRITE)
                   ENDDO
                   IF(NMOD .NE. 0) THEN
                      READ(FF%IU1) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*NCOL+IROW,1),IROW=1, NMOD)
                   ENDIF
                ELSE
                   DO ICOL=1, NCOL
                      READ(FF%IU1)
                   ENDDO
                   IF(NMOD.NE.0) THEN
                      READ(FF%IU1)
                   ENDIF
                ENDIF
                CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_WMAT_FFM,FFM%WMAT(1:,1:,INTYP),PAR_SUP_HANDLE%CONTEXT_WORLD)
                CALL FINIT_SCALAPACK_DISTRIBUTE_TO_ALL_CORES(PAR_SUP_HANDLE)
#else
                DO ICOL=1, NCOL
                   READ(FF%IU1) (FFM%WMAT(NWRITE*(ICOL-1)+IROW,1,INTYP),IROW=1, NWRITE)
                ENDDO
                IF(NMOD .NE. 0) THEN
                   READ(FF%IU1) (FFM%WMAT(NWRITE*NCOL+IROW,1,INTYP),IROW=1, NMOD)
                ENDIF
#endif
             ENDDO
             FFM%WMAT=FFM%WMAT/EUNIT
             IF(FFM%W1.GT.0.0_q) THEN
               NCOL=FFM%MMVAR1/NWRITE
                NMOD=MOD(FFM%MMVAR1,NWRITE)
                DO INTYP=1, MTYP_LOC
                  ! C00 (radial descriptor) (-)
                   ! For XXX
#ifdef scaLAPACK
                   CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES_ALLOCATE (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_C00_FFM)
                   IF (PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_ROW.EQ.0.AND.PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_COL.EQ.0) THEN
                      DO IB=1, FFM%NB(INTYP)
                         DO ICOL=1, NCOL
                            READ(FF%IU1) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*(ICOL-1)+IROW,IB),IROW=1, NWRITE)
                         ENDDO
                         IF(NMOD .NE. 0) THEN
                            READ(FF%IU1) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*NCOL+IROW,IB),IROW=1, NMOD)
                         ENDIF
                      ENDDO
                   ELSE
                      DO IB=1, FFM%NB(INTYP)
                         DO ICOL=1, NCOL
                            READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
                         ENDDO
                         IF(NMOD .NE. 0) THEN
                            READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
                         ENDIF
                      ENDDO
                   ENDIF
                   CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_C00_FFM,FFM%C00(1:,1:,INTYP),PAR_SUP_HANDLE%CONTEXT_WORLD)
                   CALL FINIT_SCALAPACK_DISTRIBUTE_TO_ALL_CORES(PAR_SUP_HANDLE)
#else
                   DO IB=1, FFM%NB(INTYP)
                      DO ICOL=1, NCOL
                         READ(FF%IU1) (FFM%C00(NWRITE*(ICOL-1)+IROW,IB,INTYP),IROW=1, NWRITE)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1) (FFM%C00(NWRITE*NCOL+IROW,IB,INTYP),IROW=1, NMOD)
                      ENDIF
                   ENDDO
#endif
                ENDDO
             ENDIF
             IF(FFM%W2.GT.0.0_q) THEN
               DO INTYP=1, MTYP_LOC
                   NCOL=FFM%NNVAR2(INTYP)/NWRITE
                   NMOD=MOD(FFM%NNVAR2(INTYP),NWRITE)
                  ! Cnlm (angular descriptor) (-)
                   ! For XXX
#ifdef scaLAPACK
                   CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES_ALLOCATE (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_PS_FFM)
                   IF (PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_ROW.EQ.0.AND.PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_COL.EQ.0) THEN
                      DO IB=1, FFM%NB(INTYP)
                         DO ICOL=1, NCOL
                            READ(FF%IU1) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*(ICOL-1)+IROW,IB),IROW=1, NWRITE)
                         ENDDO
                         IF(NMOD .NE. 0) THEN
                            READ(FF%IU1) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*NCOL+IROW,IB),IROW=1, NMOD)
                         ENDIF
                      ENDDO
                   ELSE
                      DO IB=1, FFM%NB(INTYP)
                         DO ICOL=1, NCOL
                            READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
                         ENDDO
                         IF(NMOD .NE. 0) THEN
                            READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
                         ENDIF
                      ENDDO
                   ENDIF
                   CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_PS_FFM,FFM%PS(1:,1:,INTYP),PAR_SUP_HANDLE%CONTEXT_WORLD)
                   CALL FINIT_SCALAPACK_DISTRIBUTE_TO_ALL_CORES(PAR_SUP_HANDLE)
#else
                   DO IB=1, FFM%NB(INTYP)
                      DO ICOL=1, NCOL
                         READ(FF%IU1) (FFM%PS(NWRITE*(ICOL-1)+IROW,IB,INTYP),IROW=1, NWRITE)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1) (FFM%PS(NWRITE*NCOL+IROW,IB,INTYP),IROW=1, NMOD)
                      ENDIF
                   ENDDO
#endif
                ENDDO
             ENDIF
! Parameters for Bayesian error estimation
            ! Normalized noise parameter (-)
             READ(FF%IU1) FFM%SIGV
            ! Normalized uncertainity parameter (-)
             READ(FF%IU1) FFM%SIGW
! Here, FFM%SIG is assumed to be kBT. On the basis of the statistical
! dynamics, the variance should be given as sqrt(NIONS)*kBT.
            ! Variance in energies per atom in training data (eV atom^-1)
             READ(FF%IU1) FFM%SIG(1)
             FFM%SIG(1)=FFM%SIG(1)/EUNIT
            ! Variance in forces in training data (eV Angst^-1)
             READ(FF%IU1) FFM%SIG(2),FFM%SIG(3),FFM%SIG(4)
             FFM%SIG(2:4)=FFM%SIG(2:4)/FUNIT
            ! Variance in stress tensors in training data (kB)
             READ(FF%IU1) FFM%SIG(5),FFM%SIG(6),FFM%SIG(7)
             READ(FF%IU1) FFM%SIG(8),FFM%SIG(9),FFM%SIG(10)
             FFM%SIG(5:10)=FFM%SIG(5:10)/SUNIT
             FFM%SIG_SYS(1:10,ISYS_NEW)=FFM%SIG(1:10)
             FFM%NBTOT=0
             DO INTYP=1, MTYP_LOC
                FFM%NBTOT=FFM%NBTOT+FFM%NB(INTYP)
             ENDDO
             FFM%LBHEAD(1)=0
             DO INTYP=2, MTYP_LOC
                FFM%LBHEAD(INTYP)=FFM%LBHEAD(INTYP-1)+FFM%NB(INTYP-1)
             ENDDO
             NCOL=FFM%NBTOT/NWRITE
             NMOD=MOD(FFM%NBTOT,NWRITE)
            ! Covariance matrix (-)
#ifdef scaLAPACK
             CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES_ALLOCATE (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_CMAT_FFM)
             IF (PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_ROW.EQ.0.AND.PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_COL.EQ.0) THEN
                DO IB=1, FFM%NBTOT
                   DO ICOL=1, NCOL
                      READ(FF%IU1) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*(ICOL-1)+IROW,IB),IROW=1, NWRITE)
                   ENDDO
                   IF(NMOD .NE. 0) THEN
                      READ(FF%IU1) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*NCOL+IROW,IB),IROW=1, NMOD)
                   ENDIF
                ENDDO
             ELSE
                DO IB=1, FFM%NBTOT
                   DO ICOL=1, NCOL
                      READ(FF%IU1) (RDUMMY(IROW),IROW=1, NWRITE)
                   ENDDO
                   IF(NMOD .NE. 0) THEN
                      READ(FF%IU1) (RDUMMY(IROW),IROW=1, NMOD)
                   ENDIF
                ENDDO
             ENDIF
             CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_CMAT_FFM,FFM%CMAT(1:,1:),PAR_SUP_HANDLE%CONTEXT_WORLD)
             CALL FINIT_SCALAPACK_DISTRIBUTE_TO_ALL_CORES(PAR_SUP_HANDLE)
#else
             DO IB=1, FFM%NBTOT
                DO ICOL=1, NCOL
                   READ(FF%IU1) (FFM%CMAT(NWRITE*(ICOL-1)+IROW,IB),IROW=1, NWRITE)
                ENDDO
                IF(NMOD .NE. 0) THEN
                   READ(FF%IU1) (FFM%CMAT(NWRITE*NCOL+IROW,IB),IROW=1, NMOD)
                ENDIF
             ENDDO
#endif
          ENDIF
          CLOSE(FF%IU1)
! Stop profiling
          PROFILING_STOP('read_ffcar')
        END SUBROUTINE READ_FFCAR

!****************************************************************************************************
! ML_FF readers for ASCII (non binary) format
!****************************************************************************************************
! Here are the old ML_FF readers for compatibility reasons
        SUBROUTINE READ_FF_CAR_DEF_ELEMENT_ORDER_FF_1_NON_BINARY(FFCAR,IU,MTYP_OLD)
          USE ML_FF_STRUCT, ONLY: ML_IO_WRITE
          IMPLICIT NONE
! Input variables.
          CHARACTER(LEN=9), INTENT(IN)    :: FFCAR
          INTEGER         , INTENT(IN)    :: IU
! Output variables.
          INTEGER         , INTENT(OUT)   :: MTYP_OLD
! Local variables
          INTEGER                         :: I
          INTEGER                         :: IERR
! Open file.
          OPEN(UNIT=IU,FILE=FFCAR,STATUS='OLD',IOSTAT=IERR)
! Check if error occured in opening file
          IF (IERR.NE.0) THEN
             CALL ERROR_PRINT2TERMINAL("Error in opening ML_FF file")
          ENDIF
! Do actual reading
          DO I=1, 3
             READ(IU,*)
          ENDDO
          READ(IU,*) MTYP_OLD
! Close file.
          CLOSE(IU)
        END SUBROUTINE READ_FF_CAR_DEF_ELEMENT_ORDER_FF_1_NON_BINARY

        SUBROUTINE READ_FF_CAR_DEF_ELEMENT_ORDER_FF_2_NON_BINARY(FFCAR,IU,EATOM_OLD,MTYP_OLD,POMASS_OLD,TYPE_OLD)
          IMPLICIT NONE
! Input variables.
          CHARACTER(LEN=9), INTENT(IN)    :: FFCAR
          INTEGER         , INTENT(IN)    :: IU
          INTEGER         , INTENT(IN)    :: MTYP_OLD
! Input and ouput variables.
          REAL(q)         , INTENT(OUT)   :: EATOM_OLD(:) !(1:MTYP_OLD)
          REAL(q)         , INTENT(OUT)   :: POMASS_OLD(:) !(1:MTYP_OLD)
          CHARACTER(LEN=2), INTENT(OUT)   :: TYPE_OLD(:) !(1:MTYP_OLD)
! Local variables
          INTEGER                         :: I
          INTEGER                         :: ICOL
          INTEGER                         :: IERR
          INTEGER                         :: INTYP
          INTEGER                         :: IROW
          INTEGER                         :: NCOL
          INTEGER                         :: NMOD
          INTEGER                         :: NWRITE
! Set constant
          NWRITE = 3
! Open file.
          OPEN(UNIT=IU,FILE=FFCAR,STATUS='OLD',IOSTAT=IERR)
! Check if error occured in opening file
          IF (IERR.NE.0) THEN
             CALL ERROR_PRINT2TERMINAL("Error in opening ML_FF file")
          ENDIF
! Do actual reading
          DO I=1, 3
             READ(IU,*)
          ENDDO
          READ(IU,*)
          NCOL=MTYP_OLD/NWRITE
          NMOD=MOD(MTYP_OLD,NWRITE)
          DO I=1, 3
             READ(IU,*)
          ENDDO
          DO ICOL=1, NCOL, 1
             READ(IU,*) (TYPE_OLD(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(IU,*) (TYPE_OLD(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO I=1, 3
             READ(IU,*)
          ENDDO
          DO ICOL=1, NCOL, 1
             READ(IU,*) (EATOM_OLD(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(IU,*) (EATOM_OLD(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO I=1, 3
             READ(IU,*)
          ENDDO
          DO ICOL=1, NCOL
             READ(IU,*) (POMASS_OLD(NWRITE*(ICOL-1)+IROW),IROW=1,NWRITE)
          ENDDO
          IF(NMOD.NE.0) THEN
             READ(IU,*) (POMASS_OLD(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
! Unit conversion
          EATOM_OLD(1:MTYP_OLD)=EATOM_OLD(1:MTYP_OLD)/EUNIT
          POMASS_OLD(1:MTYP_OLD)=POMASS_OLD(1:MTYP_OLD)/MUNIT
! Close file.
          CLOSE(IU)
        END SUBROUTINE READ_FF_CAR_DEF_ELEMENT_ORDER_FF_2_NON_BINARY

        SUBROUTINE OVERWRITE_FF_PAR_BY_PAR_IN_FFCAR_NON_BINARY (TAG_LIST,FF,FFM)
          USE ML_FF_STRUCT, ONLY : FF_PAR, FFM_PAR
          USE ML_FF_TAGLIST
          IMPLICIT NONE
! Input and output variables.
          TYPE (TAGLIST), INTENT(INOUT) :: TAG_LIST
          TYPE (FF_PAR)                 :: FF
          TYPE (FFM_PAR)                :: FFM
! Local variables.
          INTEGER                       :: I
          INTEGER                       :: IB
          INTEGER                       :: ICOL
          INTEGER                       :: IERR
          INTEGER                       :: INTYP
          INTEGER                       :: IROW
          INTEGER                       :: JNTYP
          INTEGER                       :: MTYP_LOC
          INTEGER                       :: NB_FFM(1:FF%MTYP)
          INTEGER                       :: NBTOT_FFM
          INTEGER                       :: NCOL
          INTEGER                       :: NMOD
          INTEGER                       :: NWRITE
          CHARACTER(LEN=2)              :: TYPE_LOC(1:FF%MTYP)
! Set constant
          NWRITE=3
! Open FFCAR file.
          OPEN(FF%IU1,FILE=FF%FFCAR,STATUS='OLD',IOSTAT=IERR)
! Check if error occured in opening file
          IF (IERR.NE.0) THEN
             CALL ERROR_PRINT2TERMINAL("Error in opening ML_FF file")
          ENDIF
! Read the parameters.
          DO I=1, 3
             READ(FF%IU1,*)
          ENDDO
          READ(FF%IU1,*) MTYP_LOC
          IF(FF%MTYP.NE.MTYP_LOC) THEN
             CALL ERROR_PRINT2TERMINAL ('!!! MLFF : you need to generate a force field before adapting ISTART_FF=2 !!!')
          ENDIF
          NCOL=FF%MTYP/NWRITE
          NMOD=MOD(FF%MTYP,NWRITE)
          DO I=1, 3
             READ(FF%IU1,*)
          ENDDO
          DO ICOL=1, NCOL, 1
             READ(FF%IU1,*) (TYPE_LOC(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1,*) (TYPE_LOC(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO I=1, FF%MTYP
             IF(TYPE_LOC(I).NE.FF%TYPE(I)) THEN
                CALL ERROR_PRINT2TERMINAL ('!!! MLFF : you need to generate force field before adapting ISTART_FF=2 !!!')
             ENDIF
          ENDDO
          DO I=1, 3
             READ(FF%IU1,*)
          ENDDO
          DO ICOL=1, NCOL, 1
             READ(FF%IU1,*) (FF%EATOM_REF(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1,*) (FF%EATOM_REF(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO INTYP=1, MTYP_LOC
             FF%EATOM_REF(INTYP)=FF%EATOM_REF(INTYP)/EUNIT
          ENDDO
          DO I=1, 3
             READ(FF%IU1,*)
          ENDDO
          DO ICOL=1, NCOL, 1
             READ(FF%IU1,*)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1,*)
          ENDIF
          DO I=1, 5
             READ(FF%IU1,*)
          ENDDO
          READ(FF%IU1,*) FF%IWEIGHT
          DO I=1, 3
             READ(FF%IU1,*)
          ENDDO
          READ(FF%IU1,*) FF%WTOTEN,FF%WTIFOR,FF%WTSIF
          DO I=1, 5
             READ(FF%IU1,*)
          ENDDO
          READ(FF%IU1,*) FF%STOTEN
          FF%STOTEN=FF%STOTEN/EUNIT
          DO I=1, 3
             READ(FF%IU1,*)
          ENDDO
          READ(FF%IU1,*) FF%STIFOR
          FF%STIFOR=FF%STIFOR/FUNIT
          DO I=1, 3
             READ(FF%IU1,*)
          ENDDO
          READ(FF%IU1,*) FF%STSIF
          FF%STSIF=FF%STSIF/SUNIT
! Read many-body parameters
          DO I=1, 3
             READ(FF%IU1,*)
          ENDDO
          READ(FF%IU1,*) FF%LMLMB
          IF(FF%LMLMB) THEN
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%LSIC
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%LSUPERVEC
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%W1,FFM%W2
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%ICUT1,FFM%ICUT2
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%RCUT1,FFM%RCUT2
             FFM%RCUT1=FFM%RCUT1/AUTOA
             FFM%RCUT2=FFM%RCUT2/AUTOA
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%IBROAD1,FFM%IBROAD2
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%SION1,FFM%SION2
             FFM%SION1=FFM%SION1/AUTOA
             FFM%SION2=FFM%SION2/AUTOA
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%MRB1,FFM%MRB2
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%NR1,FFM%NR2
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%MSPL1,FFM%MSPL2
             FFM%NSPL1=FFM%MSPL1
             FFM%NSPL2=FFM%MSPL2
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%LMAX1,FFM%LMAX2
             IF(FFM%W1.GT.0.0_q) THEN
                FFM%MMVAR1=FFM%MRB1*MTYP_LOC
                FFM%NNVAR1=FFM%MMVAR1
             ENDIF
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%NHYP1,FFM%NHYP2
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%LNORM1,FFM%LNORM2
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%LWINDOW1,FFM%LWINDOW2
             IF(FFM%LWINDOW1) THEN
                DO I=1, 3
                   READ(FF%IU1,*)
                ENDDO
                READ(FF%IU1,*) FFM%IWINDOW1
             ENDIF
             IF(FFM%LWINDOW2) THEN
                DO I=1, 3
                   READ(FF%IU1,*)
                ENDDO
                READ(FF%IU1,*) FFM%IWINDOW2
             ENDIF
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%LAFILT2
             IF(FFM%LAFILT2) THEN
                DO I=1, 3
                   READ(FF%IU1,*)
                ENDDO
                READ(FF%IU1,*) FFM%IAFILT2
                IF(FFM%IAFILT2.EQ.2) THEN
                   DO I=1, 3
                      READ(FF%IU1,*)
                   ENDDO
                   READ(FF%IU1,*) FFM%AFILT2
                ENDIF
             ENDIF
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%LMETRIC1,FFM%LMETRIC2
             IF(FFM%LMETRIC1) THEN
                DO I=1, 3
                   READ(FF%IU1,*)
                ENDDO
                READ(FF%IU1,*) FFM%NMETRIC1,FFM%RMETRIC1
             ENDIF
             IF(FFM%LMETRIC2) THEN
                DO I=1, 3
                   READ(FF%IU1,*)
                ENDDO
                READ(FF%IU1,*) FFM%NMETRIC2,FFM%RMETRIC2
             ENDIF
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%LVARTRAN1,FFM%LVARTRAN2
             IF(FFM%LVARTRAN1) THEN
                DO I=1, 3
                   READ(FF%IU1,*)
                ENDDO
                READ(FF%IU1,*) FFM%NVARTRAN1
             ENDIF
             IF(FFM%LVARTRAN2) THEN
                DO I=1, 3
                   READ(FF%IU1,*)
                ENDDO
                READ(FF%IU1,*) FFM%NVARTRAN2
             ENDIF
             IF(FFM%W2.GT.0.0_q) THEN
! The number of radial basis sets for each angular component.
                DO I=1, 3
                   READ(FF%IU1,*)
                ENDDO
                NCOL=(FFM%LMAX2+1)/NWRITE
                NMOD=MOD(FFM%LMAX2+1,NWRITE)
                DO ICOL=1, NCOL
                   READ(FF%IU1,*)
                ENDDO
                IF(NMOD .NE. 0) THEN
                   READ(FF%IU1,*)
                ENDIF
! Lists of active descriptors.
                DO INTYP=1, MTYP_LOC
                   DO ICOL=1, 3
                      READ(FF%IU1,*)
                   ENDDO
                   READ(FF%IU1,*) FFM%NNVAR2(INTYP)
                   DO ICOL=1, 3
                      READ(FF%IU1,*)
                   ENDDO
                   NCOL=(FFM%NNVAR2(INTYP))/NWRITE
                   NMOD=MOD(FFM%NNVAR2(INTYP),NWRITE)
                   DO ICOL=1, NCOL
                      READ(FF%IU1,*)
                   ENDDO
                   IF(NMOD .NE. 0) THEN
                      READ(FF%IU1,*)
                   ENDIF
                ENDDO
                FFM%MMVAR2=FFM%NNVAR2(1)
                DO INTYP=2, MTYP_LOC, 1
                   FFM%MMVAR2=MAX(FFM%MMVAR2,FFM%NNVAR2(INTYP))
                ENDDO
! Lists of active descriptors for self-interaction correction if necessary.
                IF(FFM%LSIC) THEN
                   DO INTYP=1, MTYP_LOC
                      DO ICOL=1, 3
                         READ(FF%IU1,*)
                      ENDDO
                      READ(FF%IU1,*) FFM%NNVAR2_SIC(INTYP)
                      DO ICOL=1, 3
                         READ(FF%IU1,*)
                      ENDDO
                      NCOL=(FFM%NNVAR2_SIC(INTYP))/NWRITE
                      NMOD=MOD(FFM%NNVAR2_SIC(INTYP),NWRITE)
                      DO ICOL=1, NCOL
                         READ(FF%IU1,*)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1,*)
                      ENDIF
                   ENDDO
                   FFM%MMVAR2_SIC=FFM%NNVAR2_SIC(1)
                   DO INTYP=2, MTYP_LOC, 1
                      FFM%MMVAR2_SIC=MAX(FFM%MMVAR2_SIC,FFM%NNVAR2_SIC(INTYP))
                   ENDDO
                ENDIF
             ENDIF
             NCOL=MTYP_LOC/NWRITE
             NMOD=MOD(MTYP_LOC,NWRITE)
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             DO ICOL=1, NCOL
                READ(FF%IU1,*) (NB_FFM(NWRITE*(ICOL-1)+IROW), IROW=1, NWRITE)
             ENDDO
             IF(NMOD .NE. 0) THEN
                READ(FF%IU1,*) (NB_FFM(NWRITE*NCOL+IROW),IROW=1, NMOD)
             ENDIF
! Set the required size automatically to the largest required if
! ML_FF_ISTART=2
! If ML_MB was not read in then set to the first to
! first entry in ML_FF file and then to the largest 
             IF (TAG_LIST%GET_STATE("ML_MB") == TL_DEFAULT) THEN
                CALL TAG_LIST%SET_STATE("ML_MB", TL_OVERRIDE)
                FFM%MB=NB_FFM(1)+10
             ENDIF
             DO INTYP=1, MTYP_LOC
                IF(NB_FFM(INTYP) .GT. FFM%MB) THEN
                    FFM%MB=NB_FFM(INTYP)+10
                ENDIF
             ENDDO
! Further reading
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%ISCALE_TOTEN
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%TOTENAV
             FFM%TOTENAV=FFM%TOTENAV/EUNIT
             DO I=1, 2
                READ(FF%IU1,*)
             ENDDO
             DO INTYP=1, MTYP_LOC
                NCOL=NB_FFM(INTYP)/NWRITE
                NMOD=MOD(NB_FFM(INTYP),NWRITE)
                DO I=1, 3
                   READ(FF%IU1,*)
                ENDDO
                DO ICOL=1, NCOL
                   READ(FF%IU1,*)
                ENDDO
                IF(NMOD .NE. 0) THEN
                   READ(FF%IU1,*)
                ENDIF
             ENDDO
             IF(FFM%W1 .GT. 0.0_q) THEN
                DO I=1, 2
                   READ(FF%IU1,*)
                ENDDO
                NCOL=FFM%MMVAR1/NWRITE
                NMOD=MOD(FFM%MMVAR1,NWRITE)
                DO INTYP=1, MTYP_LOC
                   DO I=1, 3
                      READ(FF%IU1,*)
                   ENDDO
                   DO IB=1, NB_FFM(INTYP)
                      DO ICOL=1, NCOL
                         READ(FF%IU1,*)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1,*)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF
             IF(FFM%W2 .GT. 0.0_q) THEN
                DO I=1, 2
                   READ(FF%IU1,*)
                ENDDO
                DO INTYP=1, MTYP_LOC
                   NCOL=FFM%NNVAR2(INTYP)/NWRITE
                   NMOD=MOD(FFM%NNVAR2(INTYP),NWRITE)
                   DO I=1, 3
                      READ(FF%IU1,*)
                   ENDDO
                   DO IB=1, NB_FFM(INTYP)
                      DO ICOL=1, NCOL
                         READ(FF%IU1,*)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1,*)
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF
             DO I=1, 5
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%SIGV
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%SIGW
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*)
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*)
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*)
             READ(FF%IU1,*)
             NBTOT_FFM=0
             DO INTYP=1, MTYP_LOC
                NBTOT_FFM=NBTOT_FFM+NB_FFM(INTYP)
             ENDDO
             NCOL=NBTOT_FFM/NWRITE
             NMOD=MOD(NBTOT_FFM,NWRITE)
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             DO IB=1, NBTOT_FFM
                DO ICOL=1, NCOL
                   READ(FF%IU1,*)
                ENDDO
                IF(NMOD.NE.0) THEN
                   READ(FF%IU1,*)
                ENDIF
             ENDDO
          ENDIF
          CLOSE(FF%IU1)
        END SUBROUTINE OVERWRITE_FF_PAR_BY_PAR_IN_FFCAR_NON_BINARY


!****************************************************************************************************
! Read FFCAR file.
!****************************************************************************************************

        SUBROUTINE READ_FFCAR_NON_BINARY (FF,FFM,PAR_SUP_HANDLE,ISYS_NEW)
          USE ML_FF_STRUCT, ONLY : FF_PAR, FFM_PAR, PARALLEL_SUPER
          IMPLICIT NONE
! Input variables.
          TYPE (FF_PAR)                   :: FF
          TYPE (FFM_PAR)                  :: FFM
          TYPE (PARALLEL_SUPER)           :: PAR_SUP_HANDLE
          INTEGER         , INTENT(IN)    :: ISYS_NEW
! Local variables.
          INTEGER                         :: I
          INTEGER                         :: IB
          INTEGER                         :: ICOL
          INTEGER                         :: IERR
          INTEGER                         :: INTYP
          INTEGER                         :: INTYP0
          INTEGER                         :: IRB
          INTEGER                         :: IROW
          INTEGER                         :: IVAR
          INTEGER                         :: JNTYP
          INTEGER                         :: JNTYP0
          INTEGER                         :: JJNTYP0
          INTEGER                         :: JRB
          INTEGER                         :: L
          INTEGER                         :: LVAR_HELP(1:FFM%MMVAR2,1:FF%MTYP)
          INTEGER                         :: LVAR_SIC_HELP(1:FFM%MMVAR2_SIC,1:FF%MTYP)
          INTEGER                         :: MTYP_LOC
          INTEGER                         :: NCOL
          INTEGER                         :: NMOD
          INTEGER                         :: NNVAR2_HELP1
          INTEGER                         :: NNVAR2_HELP2
          INTEGER                         :: NNVAR2_SIC_HELP1
          INTEGER                         :: NNVAR2_SIC_HELP2
          INTEGER                         :: NRB2(0:FFM%LMAX2)
          INTEGER                         :: NROW
          INTEGER                         :: NWRITE
          CHARACTER(LEN=2)                :: TYPE_LOC(1:FF%MTYP)
! Start profiling
          PROFILING_START('read_ffcar')
! Open FF%FFCAR file.
          OPEN(FF%IU1,FILE=FF%FFCAR,STATUS='OLD',IOSTAT=IERR)
! Check if error occured in opening file
          IF (IERR.NE.0) THEN
             CALL ERROR_PRINT2TERMINAL("Error in opening ML_FF file")
          ENDIF
! Set constant
          NWRITE=3
! Read the parameters.
          DO I=1, 3
             READ(FF%IU1,*)
          ENDDO
          READ(FF%IU1,*) MTYP_LOC
          NCOL=MTYP_LOC/NWRITE
          NMOD=MOD(MTYP_LOC,NWRITE)
          DO I=1, 3
             READ(FF%IU1,*)
          ENDDO
          DO ICOL=1, NCOL, 1
             READ(FF%IU1,*) (TYPE_LOC(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1,*) (TYPE_LOC(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO I=1, 3
             READ(FF%IU1,*)
          ENDDO
          DO ICOL=1, NCOL, 1
             READ(FF%IU1,*) (FF%EATOM_REF(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1,*) (FF%EATOM_REF(NWRITE*NCOL+IROW),IROW=1, NMOD)
          ENDIF
          DO INTYP=1, MTYP_LOC
             FF%EATOM_REF(INTYP)=FF%EATOM_REF(INTYP)/EUNIT
          ENDDO
          DO I=1, 3
             READ(FF%IU1,*)
          ENDDO
          DO ICOL=1, NCOL, 1
             READ(FF%IU1,*)
          ENDDO
          IF(NMOD .NE. 0) THEN
             READ(FF%IU1,*)
          ENDIF
          DO I=1, 5
             READ(FF%IU1,*)
          ENDDO
          READ(FF%IU1,*) FF%IWEIGHT
          DO I=1, 3
             READ(FF%IU1,*)
          ENDDO
          READ(FF%IU1,*) FF%WTOTEN,FF%WTIFOR,FF%WTSIF
          IF(FF%IWEIGHT.EQ.1) THEN
             FF%WTIFOR=FF%WTIFOR/FUNIT
             FF%WTOTEN=FF%WTOTEN/EUNIT
             FF%WTSIF=FF%WTSIF/SUNIT
          ENDIF
          DO I=1, 5
             READ(FF%IU1,*)
          ENDDO
          READ(FF%IU1,*) FF%STOTEN
          FF%STOTEN=FF%STOTEN/EUNIT
          DO I=1, 3
             READ(FF%IU1,*)
          ENDDO
          READ(FF%IU1,*) FF%STIFOR
          FF%STIFOR=FF%STIFOR/FUNIT
          DO I=1, 3
             READ(FF%IU1,*)
          ENDDO
          READ(FF%IU1,*) FF%STSIF
          FF%STSIF=FF%STSIF/SUNIT
! Read many-body parameters
          DO I=1, 3
             READ(FF%IU1,*)
          ENDDO
          READ(FF%IU1,*) FF%LMLMB
          IF(FF%LMLMB) THEN
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%LSIC
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%LSUPERVEC
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%W1,FFM%W2
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%ICUT1,FFM%ICUT2
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%RCUT1,FFM%RCUT2
             FFM%RCUT1=FFM%RCUT1/AUTOA
             FFM%RCUT2=FFM%RCUT2/AUTOA
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%IBROAD1,FFM%IBROAD2
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%SION1,FFM%SION2
             FFM%SION1=FFM%SION1/AUTOA
             FFM%SION2=FFM%SION2/AUTOA
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*)
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%NR1,FFM%NR2
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%MSPL1,FFM%MSPL2
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*)
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%NHYP1,FFM%NHYP2
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%LNORM1,FFM%LNORM2
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%LWINDOW1,FFM%LWINDOW2
             IF(FFM%LWINDOW1) THEN
                DO I=1, 3
                   READ(FF%IU1,*)
                ENDDO
                READ(FF%IU1,*) FFM%IWINDOW1
             ENDIF
             IF(FFM%LWINDOW2) THEN
                DO I=1, 3
                   READ(FF%IU1,*)
                ENDDO
                READ(FF%IU1,*) FFM%IWINDOW2
             ENDIF
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%LAFILT2
             IF(FFM%LAFILT2) THEN
                DO I=1, 3
                   READ(FF%IU1,*)
                ENDDO
                READ(FF%IU1,*) FFM%IAFILT2
                IF(FFM%IAFILT2.EQ.2) THEN
                   DO I=1, 3
                      READ(FF%IU1,*)
                   ENDDO
                   READ(FF%IU1,*) FFM%AFILT2
                ENDIF
             ENDIF
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%LMETRIC1,FFM%LMETRIC2
             IF(FFM%LMETRIC1) THEN
                DO I=1, 3
                   READ(FF%IU1,*)
                ENDDO
                READ(FF%IU1,*) FFM%NMETRIC1,FFM%RMETRIC1
             ENDIF
             IF(FFM%LMETRIC2) THEN
                DO I=1, 3
                   READ(FF%IU1,*)
                ENDDO
                READ(FF%IU1,*) FFM%NMETRIC2,FFM%RMETRIC2
             ENDIF
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%LVARTRAN1,FFM%LVARTRAN2
             IF(FFM%LVARTRAN1) THEN
                DO I=1, 3
                   READ(FF%IU1,*)
                ENDDO
                READ(FF%IU1,*) FFM%NVARTRAN1
             ENDIF
             IF(FFM%LVARTRAN2) THEN
                DO I=1, 3
                   READ(FF%IU1,*)
                ENDDO
                READ(FF%IU1,*) FFM%NVARTRAN2
             ENDIF
             IF(FFM%W2.GT.0.0_q) THEN
! Read the number of radial basis sets for each angular components.
                NCOL=(FFM%LMAX2+1)/NWRITE
                NMOD=MOD(FFM%LMAX2+1,NWRITE)
                DO I=1, 3
                   READ(FF%IU1,*)
                ENDDO
                DO ICOL=1, NCOL, 1
                   READ(FF%IU1,*) (NRB2(NWRITE*(ICOL-1)+IROW),IROW=0, NWRITE-1)
                ENDDO
                IF(NMOD .NE. 0) THEN
                   READ(FF%IU1,*) (NRB2(NWRITE*NCOL+IROW),IROW=0, NMOD-1)
                ENDIF
! Read the list of active descriptors.
                DO INTYP=1, MTYP_LOC
                   DO I=1, 3
                      READ(FF%IU1,*)
                   ENDDO
                   READ(FF%IU1,*) FFM%NNVAR2(INTYP)
                   NCOL=(FFM%NNVAR2(INTYP))/NWRITE
                   NMOD=MOD(FFM%NNVAR2(INTYP),NWRITE)
                   DO I=1, 3
                      READ(FF%IU1,*)
                   ENDDO
                   DO ICOL=1, NCOL, 1
                      READ(FF%IU1,*) (LVAR_HELP(NWRITE*(ICOL-1)+IROW,INTYP),IROW=1, NWRITE)
                   ENDDO
                   IF(NMOD .NE. 0) THEN
                      READ(FF%IU1,*) (LVAR_HELP(NWRITE*NCOL+IROW,INTYP),IROW=1, NMOD)
                   ENDIF
                ENDDO
! Read the list of active descriptors for SIC.
                IF(FFM%LSIC) THEN
                   DO INTYP=1, MTYP_LOC
                      DO I=1, 3
                         READ(FF%IU1,*)
                      ENDDO
                      READ(FF%IU1,*) FFM%NNVAR2_SIC(INTYP)
                      NCOL=(FFM%NNVAR2_SIC(INTYP))/NWRITE
                      NMOD=MOD(FFM%NNVAR2_SIC(INTYP),NWRITE)
                      DO I=1, 3
                         READ(FF%IU1,*)
                      ENDDO
                      DO ICOL=1, NCOL, 1
                         READ(FF%IU1,*) (LVAR_SIC_HELP(NWRITE*(ICOL-1)+IROW,INTYP),IROW=1, NWRITE)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1,*) (LVAR_SIC_HELP(NWRITE*NCOL+IROW,INTYP),IROW=1, NMOD)
                      ENDIF
                   ENDDO
                ENDIF
! Set lists.
! Determination of the numbers of variables related to the angular descriptor.
                FFM%LVAR=0
                FFM%LFLAG_VAR=.FALSE.
                FFM%WVAR=0.0_q
                IF(FFM%LSIC) THEN
                   FFM%LVAR_SIC=0
                   FFM%LFLAG_VAR_SIC=.FALSE.
                ENDIF
                DO INTYP0=1, FF%MTYP
                   NNVAR2_HELP1=0
                   NNVAR2_HELP2=0
                   IF(FFM%LSIC) THEN
                      NNVAR2_SIC_HELP1=0
                      NNVAR2_SIC_HELP2=0
                   ENDIF
                   DO JNTYP0=1, FF%MTYP
                      DO JJNTYP0=1, FF%MTYP
                         DO L=0, FFM%LMAX2
                            DO IRB=1, NRB2(L)
                               DO JRB=IRB, NRB2(L)
                                  NNVAR2_HELP1=NNVAR2_HELP1+1
                                  DO IVAR=1, FFM%NNVAR2(INTYP0)
                                     IF(NNVAR2_HELP1.EQ.LVAR_HELP(IVAR,INTYP0)) THEN
                                        NNVAR2_HELP2=NNVAR2_HELP2+1
                                        FFM%LVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=NNVAR2_HELP2
                                        FFM%LFLAG_VAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=.TRUE.
                                        IF(JRB.EQ.IRB) THEN
                                           FFM%WVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=1.0_q
                                        ELSE
                                           FFM%WVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=SQRT(2.0_q)
                                        ENDIF
                                        EXIT
                                     ENDIF
                                  ENDDO
                                  IF(FFM%LSIC) THEN
                                     IF(JJNTYP0.EQ.JNTYP0) THEN
                                        NNVAR2_SIC_HELP1=NNVAR2_SIC_HELP1+1
                                        DO IVAR=1, FFM%NNVAR2_SIC(INTYP0)
                                           IF(NNVAR2_SIC_HELP1.EQ.LVAR_SIC_HELP(IVAR,INTYP0)) THEN
                                              NNVAR2_SIC_HELP2=NNVAR2_SIC_HELP2+1
                                              FFM%LVAR_SIC(JRB,IRB,L,JNTYP0,INTYP0)=NNVAR2_SIC_HELP2
                                              FFM%LFLAG_VAR_SIC(JRB,IRB,L,JNTYP0,INTYP0)=.TRUE.
                                              EXIT
                                           ENDIF
                                        ENDDO
                                     ENDIF
                                  ENDIF
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
! Read basis sets.
             NCOL=MTYP_LOC/NWRITE
             NMOD=MOD(MTYP_LOC,NWRITE)
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             DO ICOL=1, NCOL, 1
                READ(FF%IU1,*) (FFM%NB(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
             ENDDO
             IF(NMOD .NE. 0) THEN
                READ(FF%IU1,*) (FFM%NB(NWRITE*NCOL+IROW),IROW=1, NMOD)
             ENDIF
! Read the scaling method.
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%ISCALE_TOTEN
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%TOTENAV
             FFM%TOTENAV=FFM%TOTENAV/EUNIT
             DO I=1, 2
                READ(FF%IU1,*)
             ENDDO
             DO INTYP=1, MTYP_LOC
                NCOL=FFM%NB(INTYP)/NWRITE
                NMOD=MOD(FFM%NB(INTYP),NWRITE)
                DO I=1, 3
                   READ(FF%IU1,*)
                ENDDO
#ifdef scaLAPACK
                CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES_ALLOCATE(PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_WMAT_FFM)
                IF (PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_ROW.EQ.0.AND.PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_COL.EQ.0) THEN
                   DO ICOL=1, NCOL
                      READ(FF%IU1,*) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*(ICOL-1)+IROW,1),IROW=1, NWRITE)
                   ENDDO
                   IF(NMOD .NE. 0) THEN
                      READ(FF%IU1,*) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*NCOL+IROW,1),IROW=1, NMOD)
                   ENDIF
                ELSE
                   DO ICOL=1, NCOL
                      READ(FF%IU1,*)
                   ENDDO
                   IF(NMOD.NE.0) THEN
                      READ(FF%IU1,*)
                   ENDIF
                ENDIF
                CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_WMAT_FFM,FFM%WMAT(1:,1:,INTYP),PAR_SUP_HANDLE%CONTEXT_WORLD)
                CALL FINIT_SCALAPACK_DISTRIBUTE_TO_ALL_CORES(PAR_SUP_HANDLE)
#else
                DO ICOL=1, NCOL
                   READ(FF%IU1,*) (FFM%WMAT(NWRITE*(ICOL-1)+IROW,1,INTYP),IROW=1, NWRITE)
                ENDDO
                IF(NMOD .NE. 0) THEN
                   READ(FF%IU1,*) (FFM%WMAT(NWRITE*NCOL+IROW,1,INTYP),IROW=1, NMOD)
                ENDIF
#endif
             ENDDO
             FFM%WMAT=FFM%WMAT/EUNIT
             IF(FFM%W1.GT.0.0_q) THEN
                DO I=1, 2
                   READ(FF%IU1,*)
                ENDDO
                NCOL=FFM%MMVAR1/NWRITE
                NMOD=MOD(FFM%MMVAR1,NWRITE)
                DO INTYP=1, MTYP_LOC
                   DO I=1, 3
                      READ(FF%IU1,*)
                   ENDDO
#ifdef scaLAPACK
                   CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES_ALLOCATE (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_C00_FFM)
                   IF (PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_ROW.EQ.0.AND.PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_COL.EQ.0) THEN
                      DO IB=1, FFM%NB(INTYP)
                         DO ICOL=1, NCOL
                            READ(FF%IU1,*) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*(ICOL-1)+IROW,IB),IROW=1, NWRITE)
                         ENDDO
                         IF(NMOD .NE. 0) THEN
                            READ(FF%IU1,*) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*NCOL+IROW,IB),IROW=1, NMOD)
                         ENDIF
                      ENDDO
                   ELSE
                      DO IB=1, FFM%NB(INTYP)
                         DO ICOL=1, NCOL
                            READ(FF%IU1,*)
                         ENDDO
                         IF(NMOD .NE. 0) THEN
                            READ(FF%IU1,*)
                         ENDIF
                      ENDDO
                   ENDIF
                   CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_C00_FFM,FFM%C00(1:,1:,INTYP),PAR_SUP_HANDLE%CONTEXT_WORLD)
                   CALL FINIT_SCALAPACK_DISTRIBUTE_TO_ALL_CORES(PAR_SUP_HANDLE)
#else
                   DO IB=1, FFM%NB(INTYP)
                      DO ICOL=1, NCOL
                         READ(FF%IU1,*) (FFM%C00(NWRITE*(ICOL-1)+IROW,IB,INTYP),IROW=1, NWRITE)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1,*) (FFM%C00(NWRITE*NCOL+IROW,IB,INTYP),IROW=1, NMOD)
                      ENDIF
                   ENDDO
#endif
                ENDDO
             ENDIF
             IF(FFM%W2.GT.0.0_q) THEN
                DO I=1, 2
                   READ(FF%IU1,*)
                ENDDO
                DO INTYP=1, MTYP_LOC
                   NCOL=FFM%NNVAR2(INTYP)/NWRITE
                   NMOD=MOD(FFM%NNVAR2(INTYP),NWRITE)
                   DO I=1, 3
                      READ(FF%IU1,*)
                   ENDDO
#ifdef scaLAPACK
                   CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES_ALLOCATE (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_PS_FFM)
                   IF (PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_ROW.EQ.0.AND.PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_COL.EQ.0) THEN
                      DO IB=1, FFM%NB(INTYP)
                         DO ICOL=1, NCOL
                            READ(FF%IU1,*) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*(ICOL-1)+IROW,IB),IROW=1, NWRITE)
                         ENDDO
                         IF(NMOD .NE. 0) THEN
                            READ(FF%IU1,*) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*NCOL+IROW,IB),IROW=1, NMOD)
                         ENDIF
                      ENDDO
                   ELSE
                      DO IB=1, FFM%NB(INTYP)
                         DO ICOL=1, NCOL
                            READ(FF%IU1,*)
                         ENDDO
                         IF(NMOD .NE. 0) THEN
                            READ(FF%IU1,*)
                         ENDIF
                      ENDDO
                   ENDIF
                   CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_PS_FFM,FFM%PS(1:,1:,INTYP),PAR_SUP_HANDLE%CONTEXT_WORLD)
                   CALL FINIT_SCALAPACK_DISTRIBUTE_TO_ALL_CORES(PAR_SUP_HANDLE)
#else
                   DO IB=1, FFM%NB(INTYP)
                      DO ICOL=1, NCOL
                         READ(FF%IU1,*) (FFM%PS(NWRITE*(ICOL-1)+IROW,IB,INTYP),IROW=1, NWRITE)
                      ENDDO
                      IF(NMOD .NE. 0) THEN
                         READ(FF%IU1,*) (FFM%PS(NWRITE*NCOL+IROW,IB,INTYP),IROW=1, NMOD)
                      ENDIF
                   ENDDO
#endif
                ENDDO
             ENDIF
             DO I=1, 5
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%SIGV
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%SIGW
! Here, FFM%SIG is assumed to be kBT. On the basis of the statistical
! dynamics, the variance should be given as sqrt(NIONS)*kBT.
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%SIG(1)
             FFM%SIG(1)=FFM%SIG(1)/EUNIT
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%SIG(2),FFM%SIG(3),FFM%SIG(4)
             FFM%SIG(2:4)=FFM%SIG(2:4)/FUNIT
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
             READ(FF%IU1,*) FFM%SIG(5),FFM%SIG(6),FFM%SIG(7)
             READ(FF%IU1,*) FFM%SIG(8),FFM%SIG(9),FFM%SIG(10)
             FFM%SIG(5:10)=FFM%SIG(5:10)/SUNIT
             FFM%SIG_SYS(1:10,ISYS_NEW)=FFM%SIG(1:10)
             FFM%NBTOT=0
             DO INTYP=1, MTYP_LOC
                FFM%NBTOT=FFM%NBTOT+FFM%NB(INTYP)
             ENDDO
             FFM%LBHEAD(1)=0
             DO INTYP=2, MTYP_LOC
                FFM%LBHEAD(INTYP)=FFM%LBHEAD(INTYP-1)+FFM%NB(INTYP-1)
             ENDDO
             NCOL=FFM%NBTOT/NWRITE
             NMOD=MOD(FFM%NBTOT,NWRITE)
             DO I=1, 3
                READ(FF%IU1,*)
             ENDDO
#ifdef scaLAPACK
             CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES_ALLOCATE (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_CMAT_FFM)
             IF (PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_ROW.EQ.0.AND.PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_COL.EQ.0) THEN
                DO IB=1, FFM%NBTOT
                   DO ICOL=1, NCOL
                      READ(FF%IU1,*) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*(ICOL-1)+IROW,IB),IROW=1, NWRITE)
                   ENDDO
                   IF(NMOD .NE. 0) THEN
                      READ(FF%IU1,*) (PAR_SUP_HANDLE%GENERAL_SINGLE(NWRITE*NCOL+IROW,IB),IROW=1, NMOD)
                   ENDIF
                ENDDO
             ELSE
                DO IB=1, FFM%NBTOT
                   DO ICOL=1, NCOL
                      READ(FF%IU1,*)
                   ENDDO
                   IF(NMOD .NE. 0) THEN
                      READ(FF%IU1,*)
                   ENDIF
                ENDDO
             ENDIF
             CALL SCALAPACK_DISTRIBUTE_TO_ALL_CORES (PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_CMAT_FFM,FFM%CMAT(1:,1:),PAR_SUP_HANDLE%CONTEXT_WORLD)
             CALL FINIT_SCALAPACK_DISTRIBUTE_TO_ALL_CORES(PAR_SUP_HANDLE)
#else
             DO IB=1, FFM%NBTOT
                DO ICOL=1, NCOL
                   READ(FF%IU1,*) (FFM%CMAT(NWRITE*(ICOL-1)+IROW,IB),IROW=1, NWRITE)
                ENDDO
                IF(NMOD .NE. 0) THEN
                   READ(FF%IU1,*) (FFM%CMAT(NWRITE*NCOL+IROW,IB),IROW=1, NMOD)
                ENDIF
             ENDDO
#endif
          ENDIF
          CLOSE(FF%IU1)
! Stop profiling
          PROFILING_STOP('read_ffcar')
        END SUBROUTINE READ_FFCAR_NON_BINARY
!****************************************************************************************************


#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Output regression results
!****************************************************************************************************

        SUBROUTINE OUTPUT_REG (COMM_WORLD,DIR_APP,DIR_LEN,IU3,MCONF,MIONS,NCONF,NIONS,NSTEP,REGCAR, &
                   TIFOR,TIFOR_MB,TOTEN,TOTEN_MB, &
                   TSIF,TSIF_MB)
          USE ML_FF_STRUCT, ONLY : ML_MPI_PAR
          IMPLICIT NONE
! Input variables
          TYPE (ML_MPI_PAR)                :: COMM_WORLD
          CHARACTER(LEN=10), INTENT(IN) :: DIR_APP     ! Location of directory passed from VASP.
          INTEGER, INTENT(IN)           :: DIR_LEN     ! Length of DIR_APP
          INTEGER, INTENT(IN)           :: IU3
          INTEGER, INTENT(IN)           :: MCONF
          INTEGER, INTENT(IN)           :: MIONS
          INTEGER, INTENT(IN)           :: NCONF
          INTEGER, INTENT(IN)           :: NIONS(:) !(1:NCONF)
          INTEGER, INTENT(IN)           :: NSTEP
          CHARACTER(LEN=10), INTENT(IN) :: REGCAR
          REAL(q), INTENT(IN)           :: TIFOR(:,:,:) !(1:3,1:MIONS,1:MCONF)
          REAL(q), INTENT(IN)           :: TIFOR_MB(:,:,:) !(1:3,1:MIONS,1:MCONF)
          REAL(q), INTENT(IN)           :: TOTEN(:) !(1:MCONF)
          REAL(q), INTENT(IN)           :: TOTEN_MB(:) !(1:MCONF)
          REAL(q), INTENT(IN)           :: TSIF(:,:,:) !(1:3,1:3,1:MCONF)
          REAL(q), INTENT(IN)           :: TSIF_MB(:,:,:) !(1:3,1:3,1:MCONF)
! Local variables
          INTEGER ICONF
          INTEGER IERR
          INTEGER INIONS
          INTEGER IXYZ
          INTEGER JXYZ
! Output regression results
          IF(COMM_WORLD%NODE_ME.EQ.0) THEN
             OPEN(UNIT=IU3,FILE=DIR_APP(1:DIR_LEN)// &
                  REGCAR,STATUS='UNKNOWN',IOSTAT=IERR)
! Check if error occured in opening file
             IF (IERR.NE.0) THEN
                CALL ERROR_PRINT2TERMINAL("Error in opening ML_REG file")
             ENDIF
             WRITE(IU3,1) NSTEP
             WRITE(IU3,2)
             DO ICONF=1, NCONF
                WRITE(IU3,3) TOTEN(ICONF)*EUNIT,(TOTEN_MB(ICONF))*EUNIT
             ENDDO
             WRITE(IU3,4)
             DO ICONF=1, NCONF
                DO INIONS=1, NIONS(ICONF)
                   DO IXYZ=1, 3
                      WRITE(IU3,3) TIFOR(IXYZ,INIONS,ICONF)*FUNIT, &
                                    (TIFOR_MB(IXYZ,INIONS,ICONF))*FUNIT
                   ENDDO
                ENDDO
             ENDDO
             WRITE(IU3,5)
             DO ICONF=1, NCONF
                DO IXYZ=1, 3
                   DO JXYZ=IXYZ, 3
                   WRITE(IU3,3) TSIF(JXYZ,IXYZ,ICONF)*SUNIT, &
                                 (TSIF_MB(JXYZ,IXYZ,ICONF))*SUNIT
                   ENDDO
                ENDDO
             ENDDO
             CLOSE(IU3)
          ENDIF
! Format
1         FORMAT (100("*")/5X,'Results at'1X,I10,1X,'MD step.')
2         FORMAT (100("=")/5X,'Total energies (eV)'/100("-"))
3         FORMAT (5X,F12.6,1X,F12.6)
4         FORMAT (100("=")/5X,'Forces (eV ang.^-1)'/100("-"))
5         FORMAT (100("=")/5X,'Stress (kbar)'/100("-"))
6         FORMAT (100("=")/5X,'Charge (e)'/100("-"))
        END SUBROUTINE OUTPUT_REG
#endif

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Output new ab initio data
!****************************************************************************************************

        SUBROUTINE OUTPUT_ABNCAR (AB,FF,COMM_WORLD,A,ABNCAR,DIR_APP,DIR_LEN,EATOM_REF_FF,IU2,LADD,LNCONF,LNIONS,MB,MCONF,MITYP,MIONS,MTYP, &
                   NB,NCONF,NITYP,NIONS,NTYP,POMASS,POSION,SZNAM2,TIFOR,TOTEN,TSIF,TYPE)
          USE ML_FF_STRUCT, ONLY : ABINITIO, FF_PAR, ML_MPI_PAR
          IMPLICIT NONE
! Input variables
          TYPE (ABINITIO)               :: AB
          TYPE (FF_PAR)                 :: FF
          TYPE (ML_MPI_PAR)             :: COMM_WORLD
          REAL(q)          , INTENT(IN) :: A(:,:,:) !(1:3,1:3,1:MCONF)
          CHARACTER (LEN=10), INTENT(IN) :: ABNCAR
          CHARACTER(LEN=10), INTENT(IN) :: DIR_APP
          INTEGER          , INTENT(IN) :: DIR_LEN
          REAL(q)          , INTENT(IN) :: EATOM_REF_FF(:) !(1:MTYP)
          INTEGER          , INTENT(IN) :: IU2
          INTEGER          , INTENT(IN) :: LADD(:,:) !(1:MTYP,1:MCONF)
          INTEGER          , INTENT(IN) :: LNCONF(:,:) !(1:MB,1:MTYP)
          INTEGER          , INTENT(IN) :: LNIONS(:,:) !(1:MB,1:MTYP)
          INTEGER          , INTENT(IN) :: MB
          INTEGER          , INTENT(IN) :: MCONF
          INTEGER          , INTENT(IN) :: MITYP
          INTEGER          , INTENT(IN) :: MIONS
          INTEGER          , INTENT(IN) :: MTYP
          INTEGER          , INTENT(IN) :: NB(:) !(1:MTYP)
          INTEGER          , INTENT(IN) :: NCONF
          INTEGER          , INTENT(IN) :: NITYP(:,:) !(1:MTYP,1:MCONF)
          INTEGER          , INTENT(IN) :: NIONS(:) !(1:MCONF)
          INTEGER          , INTENT(IN) :: NTYP(:) !(1:MCONF)
          REAL(q)          , INTENT(IN) :: POMASS(:) !(1:MTYP)
          REAL(q)          , INTENT(IN) :: POSION(:,:,:) !(1:3,1:MIONS,1:MCONF)
          CHARACTER(LEN=40), INTENT(IN) :: SZNAM2(:) !(1:MCONF)
          REAL(q)          , INTENT(IN) :: TIFOR(:,:,:) !(1:3,1:MIONS,1:MCONF)
          REAL(q)          , INTENT(IN) :: TOTEN(:) !(1:MCONF)
          REAL(q)          , INTENT(IN) :: TSIF(:,:,:) !(1:3,1:3,1:MCONF)
          CHARACTER (LEN=2), INTENT(IN) :: TYPE(:) !(1:MTYP)
! Local variables
          INTEGER                       :: IB
          INTEGER                       :: ICOL
          INTEGER                       :: ICONF
          INTEGER                       :: IERR
          INTEGER                       :: INIONS
          INTEGER                       :: INTYP,INTYP0
          INTEGER                       :: INITYP
          INTEGER                       :: IROW
          INTEGER                       :: IXYZ,JXYZ
          INTEGER                       :: LNCONF_LOC(1:MB,1:MTYP)
          INTEGER                       :: LNIONS_LOC(1:MB,1:MTYP)
          INTEGER                       :: NB_LOC(1:MTYP)
          INTEGER                       :: NCOL
          INTEGER                       :: NMOD
          INTEGER                       :: NWRITE
! Set constant.
          NWRITE=3
          NCOL=MTYP/NWRITE
          NMOD=MOD(MTYP,NWRITE)
! Set basis set lists
          LNCONF_LOC(1:MB,1:MTYP)=LNCONF(1:MB,1:MTYP)
          LNIONS_LOC(1:MB,1:MTYP)=LNIONS(1:MB,1:MTYP)
          NB_LOC(1:MTYP)=NB(1:MTYP)
! Output all the collected ab initio data
          IF(COMM_WORLD%NODE_ME.EQ.0) THEN
             OPEN(UNIT=IU2,FILE=DIR_APP(1:DIR_LEN)//ABNCAR,STATUS='UNKNOWN',IOSTAT=IERR)
! Check if error occured in opening file
             IF (IERR.NE.0) THEN
                CALL ERROR_PRINT2TERMINAL("Error in opening ML_ABN file")
             ENDIF
             WRITE(IU2,*) "1.0 Version"
             WRITE(IU2,5) NCONF
             WRITE(IU2,6) MTYP
             WRITE(IU2,7)
             DO ICOL=1, NCOL
                WRITE(IU2,8) (TYPE(NWRITE*(ICOL-1)+IROW), IROW=1, NWRITE)
             ENDDO
             IF(NMOD .NE. 0) THEN
                WRITE(IU2,8) (TYPE(NWRITE*NCOL+IROW),IROW=1, NMOD)
             ENDIF
             WRITE(IU2,9) MIONS
             WRITE(IU2,10) MITYP
             WRITE(IU2,11)
             DO ICOL=1, NCOL
                WRITE(IU2,*) (EATOM_REF_FF(NWRITE*(ICOL-1)+IROW)*EUNIT, IROW=1, NWRITE)
             ENDDO
             IF(NMOD .NE. 0) THEN
                WRITE(IU2,*) (EATOM_REF_FF(NWRITE*NCOL+IROW)*EUNIT,IROW=1, NMOD)
             ENDIF
             WRITE(IU2,12)
             DO ICOL=1, NCOL
                WRITE(IU2,*) (POMASS(NWRITE*(ICOL-1)+IROW)*MUNIT, IROW=1, NWRITE)
             ENDDO
             IF(NMOD .NE. 0) THEN
                WRITE(IU2,*) (POMASS(NWRITE*NCOL+IROW)*MUNIT,IROW=1, NMOD)
             ENDIF
             WRITE(IU2,13)
             DO ICOL=1, NCOL
                WRITE(IU2,14) (NB_LOC(NWRITE*(ICOL-1)+IROW), IROW=1, NWRITE)
             ENDDO
             IF(NMOD .NE. 0) THEN
                WRITE(IU2,14) (NB_LOC(NWRITE*NCOL+IROW),IROW=1, NMOD)
             ENDIF
             DO INTYP0=1, MTYP
                WRITE(IU2,15) TYPE(INTYP0)
                DO IB=1, NB_LOC(INTYP0)
                   WRITE(IU2,16) LNCONF_LOC(IB,INTYP0),LNIONS_LOC(IB,INTYP0)
                ENDDO
             ENDDO
             DO ICONF=1, NCONF
                WRITE(IU2,17) ICONF
                WRITE(IU2,18) SZNAM2(ICONF)
                WRITE(IU2,19) NTYP(ICONF)
                WRITE(IU2,20) NIONS(ICONF)
                WRITE(IU2,21)
                DO INTYP=1, NTYP(ICONF)
                   INTYP0=LADD(INTYP,ICONF)
                   WRITE(IU2,22) TYPE(INTYP0),NITYP(INTYP0,ICONF)
                ENDDO
                IF (FF%LCTIFOR_WRITE_FROM_FILE) THEN
                   WRITE(IU2,30)
                   WRITE(IU2,*) AB%CTIFOR_ALLCONF_NEW(ICONF)*FUNIT
                ENDIF
                WRITE(IU2,23)
                DO IXYZ=1, 3
                   WRITE(IU2,*) (A(JXYZ,IXYZ,ICONF)*AUTOA,JXYZ=1, 3)
                ENDDO
                WRITE(IU2,24)
                DO INIONS=1, NIONS(ICONF)
                   WRITE(IU2,*) (POSION(IXYZ,INIONS,ICONF)*AUTOA,IXYZ=1, 3)
                ENDDO
                WRITE(IU2,25)
                WRITE(IU2,*) TOTEN(ICONF)*EUNIT
                WRITE(IU2,26)
                DO INIONS=1, NIONS(ICONF)
                   WRITE(IU2,*) (TIFOR(IXYZ,INIONS,ICONF)*FUNIT,IXYZ=1, 3)
                ENDDO
                WRITE(IU2,27)
                WRITE(IU2,*) TSIF(1,1,ICONF)*SUNIT,TSIF(2,2,ICONF)*SUNIT,TSIF(3,3,ICONF)*SUNIT
                WRITE(IU2,28)
                WRITE(IU2,*) TSIF(1,2,ICONF)*SUNIT,TSIF(3,2,ICONF)*SUNIT,TSIF(3,1,ICONF)*SUNIT
             ENDDO
             CLOSE(IU2)
          ENDIF
5         FORMAT (50("*")/5X,'The number of configurations'/50("-")/5X,I6)
6         FORMAT (50("*")/5X,'The maximum number of atom type'/50("-")/5X,I3)
7         FORMAT (50("*")/5X,'The atom types in the data file'/50("-"))
8         FORMAT (5X,A2,1X,A2,1X,A2)
9         FORMAT (50("*")/5X,'The maximum number of atoms per system'/50("-")/5X,I10)
10        FORMAT (50("*")/5X,'The maximum number of atoms per atom type'/50("-")/5X,I10)
11        FORMAT (50("*")/5X,'Reference atomic energy (eV)'/50("-"))
12        FORMAT (50("*")/5X,'Atomic mass'/50("-"))
13        FORMAT (50("*")/5X,'The numbers of basis sets per atom type'/50("-"))
14        FORMAT (5X,I5,1X,I5,1X,I5)
15        FORMAT (50("*")/5X,'Basis set for'1X,A2/50("-"))
16        FORMAT (5X,I6,1X,I6)
17        FORMAT (50("*")/5X,'Configuration num.'1X,I6)
18        FORMAT (50("=")/5X,'System name'/50("-")/5X,A40)
19        FORMAT (50("=")/5X,'The number of atom types'/50("-")/5X,I3)
20        FORMAT (50("=")/5X,'The number of atoms'/50("-")/5X,I6)
21        FORMAT (50("*")/5X,'Atom types and atom numbers'/50("-"))
22        FORMAT (5X,A2,1X,I6)
23        FORMAT (50("=")/5X,'Primitive lattice vectors (ang.)'/50("-"))
24        FORMAT (50("=")/5X,'Atomic positions (ang.)'/50("-"))
25        FORMAT (50("=")/5X,'Total energy (eV)'/50("-"))
26        FORMAT (50("=")/5X,'Forces (eV ang.^-1)'/50("-"))
27        FORMAT (50("=")/5X,'Stress (kbar)'/50("-")/ &
                          5X,'XX YY ZZ'/50("-"))
28        FORMAT (50("-")/5X,'XY YZ ZX'/50("-"))
29        FORMAT (50("=")/5X,'Charges (e)'/50("-"))
30        FORMAT (50("=")/5X,'CTIFOR'/50("-"))
        END SUBROUTINE OUTPUT_ABNCAR
#endif

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Output new FFCAR (FFNCAR) file.
!****************************************************************************************************
        SUBROUTINE OUT_FFNCAR (FF,FFM,PAR_SUP_HANDLE,DIR_APP,DIR_LEN)
          USE ML_FF_STRUCT, ONLY : FF_PAR,FFM_PAR,PARALLEL_SUPER
          IMPLICIT NONE
! Input variables
          TYPE (FF_PAR)                 :: FF
          TYPE (FFM_PAR)                :: FFM
          TYPE (PARALLEL_SUPER)         :: PAR_SUP_HANDLE
          CHARACTER(LEN=10), INTENT(IN) :: DIR_APP
          INTEGER          , INTENT(IN) :: DIR_LEN
! Local variables
          INTEGER                       :: IB
          INTEGER                       :: ICOL
          INTEGER                       :: IERR
          INTEGER                       :: INTYP0
          INTEGER                       :: JNTYP0
          INTEGER                       :: JJNTYP0
          INTEGER                       :: IRB
          INTEGER                       :: IROW
          INTEGER                       :: IVAR
          INTEGER                       :: JRB
          INTEGER                       :: L
          INTEGER                       :: LVAR_HELP(1:FFM%MMVAR2)
          INTEGER                       :: LVAR_SIC_HELP(1:FFM%MMVAR2_SIC)
          INTEGER                       :: NCOL
          INTEGER,PARAMETER             :: NROW = 3
          INTEGER                       :: NNVAR2_HELP1
          INTEGER                       :: NNVAR2_HELP2
! If FF%IFF=.TRUE., which means the presence of MLFF, output the FF parameters
          IF(FF%IFF) THEN
             PROFILING_START('out_ffncar')
! Open file
            OPEN(FF%IU2,FILE=DIR_APP(1:DIR_LEN)//FF%FFNCAR, &
                  STATUS='REPLACE',FORM='UNFORMATTED', &
                  ACCESS='STREAM',IOSTAT=IERR)
! Check if error occured in opening file
          IF (IERR.NE.0) THEN
             CALL ERROR_PRINT2TERMINAL("Error in opening ML_FFN file")
          ENDIF
! Output parameters
             IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
! Output the general parameters.
                ! Read the version number
                WRITE(FF%IU2) FF%ML_FF_VERSION_WRITE(1), FF%ML_FF_VERSION_WRITE(2), FF%ML_FF_VERSION_WRITE(3)
                ! The number of element type
                WRITE(FF%IU2) FF%MTYP 
                NCOL=FF%MTYP/NROW
                ! The element type
                DO ICOL=1, NCOL, 1
                   IROW=NROW*(ICOL-1)
                   WRITE(FF%IU2) (FF%TYPE(INTYP0), INTYP0=IROW+1,IROW+NROW)
                ENDDO
                IF(MOD(FF%MTYP,NROW) .GT. 0) THEN
                   WRITE(FF%IU2) (FF%TYPE(INTYP0), INTYP0=NCOL*NROW+1,NCOL*NROW+MOD(FF%MTYP,NROW))
                ENDIF
                ! Reference atomic energies (eV)
                NCOL=FF%MTYP/NROW
                DO ICOL=1, NCOL, 1
                   IROW=NROW*(ICOL-1)
                   WRITE(FF%IU2) (FF%EATOM_REF(INTYP0)*EUNIT, INTYP0=IROW+1,IROW+NROW)
                ENDDO
                IF(MOD(FF%MTYP,NROW) .GT. 0) THEN
                   WRITE(FF%IU2) (FF%EATOM_REF(INTYP0)*EUNIT, INTYP0=NCOL*NROW+1,NCOL*NROW+MOD(FF%MTYP,NROW))
                ENDIF
                ! Atomic mass
                NCOL=FF%MTYP/NROW
                DO ICOL=1, NCOL, 1
                   IROW=NROW*(ICOL-1)
                   WRITE(FF%IU2) (FF%POMASS(INTYP0)*MUNIT, INTYP0=IROW+1,IROW+NROW)
                ENDDO
                IF(MOD(FF%MTYP,NROW) .GT. 0) THEN
                   WRITE(FF%IU2) (FF%POMASS(INTYP0)*MUNIT, INTYP0=NCOL*NROW+1,NCOL*NROW+MOD(FF%MTYP,NROW))
                ENDIF
                ! Type of weighting
                WRITE(FF%IU2) FF%IWEIGHT
                IF(FF%IWEIGHT.EQ.1) THEN
                   WRITE(FF%IU2) FF%WTOTEN*EUNIT,FF%WTIFOR*FUNIT,FF%WTSIF*SUNIT
                ELSE
                   WRITE(FF%IU2) FF%WTOTEN,FF%WTIFOR,FF%WTSIF
                ENDIF
                ! Mean square root error in energy (eV atom^-1)
                WRITE(FF%IU2) FF%STOTEN*EUNIT
                ! Mean square root error in force (eV Angst^-1)
                WRITE(FF%IU2) FF%STIFOR*FUNIT
                ! Mean square root error in stress tensor (kB)
                WRITE(FF%IU2) FF%STSIF*SUNIT
! Write Many-body parameters.
                ! Many-body: true or false? - always true
                WRITE(FF%IU2) FF%LMLMB
             ENDIF
! Output the parameters for many-body interactions.
             IF(FF%LMLMB) THEN
                IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                   ! Self-interaction correction?
                   WRITE(FF%IU2) FFM%LSIC
                   ! Type of similarity measure
                   WRITE(FF%IU2) FFM%LSUPERVEC
                   ! Weights for radial and angular descriptors
                   WRITE(FF%IU2) FFM%W1,FFM%W2
                   ! Type of cutoff function
                   WRITE(FF%IU2) FFM%ICUT1,FFM%ICUT2
                   ! Cutoff radius (Angst)
                   WRITE(FF%IU2) FFM%RCUT1*AUTOA,FFM%RCUT2*AUTOA
                   ! Broadening type
                   WRITE(FF%IU2) FFM%IBROAD1,FFM%IBROAD2
                   ! Broadening parameter
                   WRITE(FF%IU2) FFM%SION1*AUTOA,FFM%SION2*AUTOA
                   ! The maximum number of radial basis set
                   WRITE(FF%IU2) FFM%MRB1,FFM%MRB2
                   ! The number of grids for integrations
                   WRITE(FF%IU2) FFM%NR1,FFM%NR2
                   ! The number of spline grids
                   WRITE(FF%IU2) FFM%NSPL1,FFM%NSPL2
                   ! Lmax for spherical harmonics
                   WRITE(FF%IU2) FFM%LMAX1,FFM%LMAX2
                   ! Polynomial parameters for SOAP
                   WRITE(FF%IU2) FFM%NHYP1,FFM%NHYP2
                   ! Normalization
                   WRITE(FF%IU2) FFM%LNORM1,FFM%LNORM2
                   ! Radial window
                   WRITE(FF%IU2) FFM%LWINDOW1,FFM%LWINDOW2
                   IF(FFM%LWINDOW1) THEN
                      ! Radial descriptor window
                      WRITE(FF%IU2) FFM%IWINDOW1
                   ENDIF
                   IF(FFM%LWINDOW2) THEN
                      ! Angular descriptor window
                      WRITE(FF%IU2) FFM%IWINDOW2
                   ENDIF
                   ! Angular filtering
                   WRITE(FF%IU2) FFM%LAFILT2
                   IF(FFM%LAFILT2) THEN
                      ! Type of angular filtering
                      WRITE(FF%IU2) FFM%IAFILT2
                      IF(FFM%IAFILT2.EQ.2) THEN
                         ! Damping parameter
                         WRITE(FF%IU2) FFM%AFILT2
                      ENDIF
                   ENDIF
                   ! Radial metric
                   WRITE(FF%IU2) FFM%LMETRIC1,FFM%LMETRIC2
                   IF(FFM%LMETRIC1) THEN
                      ! Radial descriptor metric
                      WRITE(FF%IU2) FFM%NMETRIC1,FFM%RMETRIC1*AUTOA
                   ENDIF
                   IF(FFM%LMETRIC2) THEN
                      ! Angular descriptor metric
                      WRITE(FF%IU2) FFM%NMETRIC2,FFM%RMETRIC2*AUTOA
                   ENDIF
                   ! Variable transformation
                   WRITE(FF%IU2) FFM%LVARTRAN1,FFM%LVARTRAN2
                   IF(FFM%LVARTRAN1) THEN
                      ! Radial variable transformation
                      WRITE(FF%IU2) FFM%NVARTRAN1
                   ENDIF
                   IF(FFM%LVARTRAN2) THEN
                      ! Angular variable transformation
                      WRITE(FF%IU2) FFM%NVARTRAN2
                   ENDIF
                   IF(FFM%W2.GT.0.0_q) THEN
                      NCOL=(FFM%LMAX2+1)/NROW
                      ! The number of radial basis sets for each angular channel
                      DO ICOL=1, NCOL, 1
                         IROW=NROW*(ICOL-1)
                         WRITE(FF%IU2) (FFM%NRB2(L), L=IROW,IROW+NROW-1)
                      ENDDO
                      IF(MOD(FFM%LMAX2+1,NROW) .GT. 0) THEN
                         WRITE(FF%IU2) (FFM%NRB2(L), L=NCOL*NROW,FFM%LMAX2)
                      ENDIF
                      DO INTYP0=1, FF%MTYP
                         LVAR_HELP=0
                         NNVAR2_HELP1=0
                         NNVAR2_HELP2=0
                         DO JNTYP0=1, FF%MTYP
                            DO JJNTYP0=1, FF%MTYP
                               DO L=0, FFM%LMAX2
                                  DO IRB=1, FFM%NRB2(L)
                                     DO JRB=IRB, FFM%NRB2(L)
                                        NNVAR2_HELP1=NNVAR2_HELP1+1
                                        IF(FFM%LFLAG_VAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)) THEN
                                           NNVAR2_HELP2=NNVAR2_HELP2+1
                                           LVAR_HELP(NNVAR2_HELP2)=NNVAR2_HELP1
                                        ENDIF
                                     ENDDO
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDDO
                         NCOL=NNVAR2_HELP2/NROW
                         ! The number of active descriptors for XXX
                         WRITE(FF%IU2) NNVAR2_HELP2
                         ! List of active descriptors
                         DO ICOL=1, NCOL, 1
                            IROW=NROW*(ICOL-1)
                            WRITE(FF%IU2) (LVAR_HELP(IVAR), IVAR=IROW+1,IROW+NROW)
                         ENDDO
                         IF(MOD(NNVAR2_HELP2,NROW) .GT. 0) THEN
                            WRITE(FF%IU2) (LVAR_HELP(IVAR), IVAR=NCOL*NROW+1,NNVAR2_HELP2)
                         ENDIF
                      ENDDO
                      ! Lists of active descriptors for self-interaction correction if necessary.
                      IF(FFM%LSIC) THEN
                         DO INTYP0=1, FF%MTYP
                            LVAR_SIC_HELP=0
                            NNVAR2_HELP1=0
                            NNVAR2_HELP2=0
                            DO JNTYP0=1, FF%MTYP
                               DO L=0, FFM%LMAX2
                                  DO IRB=1, FFM%NRB2(L)
                                     DO JRB=IRB, FFM%NRB2(L)
                                        NNVAR2_HELP1=NNVAR2_HELP1+1
                                        IF(FFM%LFLAG_VAR_SIC(JRB,IRB,L,JNTYP0,INTYP0)) THEN
                                           NNVAR2_HELP2=NNVAR2_HELP2+1
                                           LVAR_SIC_HELP(NNVAR2_HELP2)=NNVAR2_HELP1
                                        ENDIF
                                     ENDDO
                                  ENDDO
                               ENDDO
                            ENDDO
                            NCOL=NNVAR2_HELP2/NROW
                           ! The number of active SIC descriptors for XXX
                            WRITE(FF%IU2) NNVAR2_HELP2
                            ! List of active descriptors for SIC
                            DO ICOL=1, NCOL, 1
                               IROW=NROW*(ICOL-1)
                               WRITE(FF%IU2) (LVAR_SIC_HELP(IVAR), IVAR=IROW+1,IROW+NROW)
                            ENDDO
                            IF(MOD(NNVAR2_HELP2,NROW) .GT. 0) THEN
                               WRITE(FF%IU2) (LVAR_SIC_HELP(IVAR), IVAR=NCOL*NROW+1,NNVAR2_HELP2)
                            ENDIF
                         ENDDO
                      ENDIF
                   ENDIF
                   ! The numbers of basis sets
                   NCOL=FF%MTYP/NROW
                   DO ICOL=1, NCOL, 1
                      IROW=NROW*(ICOL-1)
                      WRITE(FF%IU2) (FFM%NB(INTYP0), INTYP0=IROW+1,IROW+NROW)
                   ENDDO
                   IF(MOD(FF%MTYP,NROW) .GT. 0) THEN
                      WRITE(FF%IU2) (FFM%NB(INTYP0), INTYP0=NCOL*NROW+1,FF%MTYP)
                   ENDIF
                   ! Type of the method to scale the energy data
                   WRITE(FF%IU2) FFM%ISCALE_TOTEN
                   ! Average energy per atom (eV atom^-1)
                   WRITE(FF%IU2) FFM%TOTENAV*EUNIT
                ENDIF
! Merge FFM%WMAT arrays if necessary
                DO INTYP0=1, FF%MTYP
#ifdef scaLAPACK
                   !Regression coefficients (eV atom^-1)
                   ! For XXX
                   CALL SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_WMAT_FFM,FFM%WMAT(1:,1:,INTYP0),PAR_SUP_HANDLE%CONTEXT_WORLD)
                   IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                      NCOL=FFM%NB(INTYP0)/NROW
                      DO ICOL=1, NCOL, 1
                         IROW=NROW*(ICOL-1)
                         WRITE(FF%IU2) (PAR_SUP_HANDLE%GENERAL_SINGLE(IB,1)*EUNIT,IB=IROW+1,IROW+NROW)
                      ENDDO
                      IF(MOD(FFM%NB(INTYP0),NROW) .GT. 0) THEN
                         WRITE(FF%IU2) (PAR_SUP_HANDLE%GENERAL_SINGLE(IB,1)*EUNIT,IB=NCOL*NROW+1,FFM%NB(INTYP0))
                      ENDIF
                   ENDIF
                   CALL FINIT_SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE)
#else
                   IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                      NCOL=FFM%NB(INTYP0)/NROW
                      DO ICOL=1, NCOL, 1
                         IROW=NROW*(ICOL-1)
                         WRITE(FF%IU2) (FFM%WMAT(IB,1,INTYP0)*EUNIT,IB=IROW+1,IROW+NROW)
                      ENDDO
                      IF(MOD(FFM%NB(INTYP0),NROW) .GT. 0) THEN
                         WRITE(FF%IU2) (FFM%WMAT(IB,1,INTYP0)*EUNIT,IB=NCOL*NROW+1,FFM%NB(INTYP0))
                      ENDIF
                   ENDIF
#endif
                ENDDO
                IF(FFM%W1.GT.0.0_q) THEN
                   NCOL=FFM%MMVAR1/NROW
                   DO INTYP0=1, FF%MTYP
                     ! C00 (radial descriptor) (-)
                      ! For XXX
! Merge FFM%C00 arrays if necessary
#ifdef scaLAPACK
                      CALL SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_C00_FFM,FFM%C00(1:,1:,INTYP0),PAR_SUP_HANDLE%CONTEXT_WORLD)
                      IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                         DO IB=1, FFM%NB(INTYP0)
                            DO ICOL=1, NCOL, 1
                               IROW=NROW*(ICOL-1)
                               WRITE(FF%IU2) (PAR_SUP_HANDLE%GENERAL_SINGLE(IVAR,IB),IVAR=IROW+1,IROW+NROW)
                            ENDDO
                            IF(MOD(FFM%MMVAR1,NROW) .GT. 0) THEN
                               WRITE(FF%IU2) (PAR_SUP_HANDLE%GENERAL_SINGLE(IVAR,IB),IVAR=NCOL*NROW+1,FFM%MMVAR1)
                            ENDIF
                         ENDDO
                      ENDIF
                      CALL FINIT_SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE)
#else
                      IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                         DO IB=1, FFM%NB(INTYP0)
                            DO ICOL=1, NCOL, 1
                               IROW=NROW*(ICOL-1)
                               WRITE(FF%IU2) (FFM%C00(IVAR,IB,INTYP0),IVAR=IROW+1,IROW+NROW)
                            ENDDO
                            IF(MOD(FFM%MMVAR1,NROW) .GT. 0) THEN
                               WRITE(FF%IU2) (FFM%C00(IVAR,IB,INTYP0),IVAR=NCOL*NROW+1,FFM%MMVAR1)
                            ENDIF
                         ENDDO
                      ENDIF
#endif
                   ENDDO
                ENDIF
! Output the angular descriptors.
! PS_SIC is not necessary because PS includes it.
                IF(FFM%W2.GT.0.0_q) THEN
                   DO INTYP0=1, FF%MTYP
                      NCOL=FFM%NNVAR2(INTYP0)/NROW
                      ! Cnlm (angular descriptor) (-)
                      ! For XXX
!Merge FFM%PS if necessary
#ifdef scaLAPACK
                      CALL SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_PS_FFM,FFM%PS(1:,1:,INTYP0),PAR_SUP_HANDLE%CONTEXT_WORLD)
                      IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                         DO IB=1, FFM%NB(INTYP0)
                            DO ICOL=1, NCOL, 1
                               IROW=NROW*(ICOL-1)
                               WRITE(FF%IU2) (PAR_SUP_HANDLE%GENERAL_SINGLE(IVAR,IB),IVAR=IROW+1,IROW+NROW)
                            ENDDO
                            IF(MOD(FFM%NNVAR2(INTYP0),NROW) .GT. 0) THEN
                               WRITE(FF%IU2) (PAR_SUP_HANDLE%GENERAL_SINGLE(IVAR,IB),IVAR=NCOL*NROW+1,FFM%NNVAR2(INTYP0))
                            ENDIF
                         ENDDO
                      ENDIF
                      CALL FINIT_SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE)
#else
                      IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                         DO IB=1, FFM%NB(INTYP0)
                            DO ICOL=1, NCOL, 1
                               IROW=NROW*(ICOL-1)
                               WRITE(FF%IU2) (FFM%PS(IVAR,IB,INTYP0),IVAR=IROW+1,IROW+NROW)
                            ENDDO
                            IF(MOD(FFM%MMVAR2,NROW) .GT. 0) THEN
                               WRITE(FF%IU2) (FFM%PS(IVAR,IB,INTYP0),IVAR=NCOL*NROW+1,FFM%NNVAR2(INTYP0))
                            ENDIF
                         ENDDO
                      ENDIF
#endif
                   ENDDO
                ENDIF
                NCOL=FFM%NBTOT/NROW
                IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                   ! Normalized noise parameter (-)
                   WRITE(FF%IU2) FFM%SIGV
                   ! Normalized uncertainity parameter (-)
                   WRITE(FF%IU2) FFM%SIGW
! In this output, we assume that SIGtot/sqrt(NIONS) is constant, where
! SIGtot denotes the variance of the total energy of the systm (not per
! atom).
! This assumption is on the basis of the fact that the variance of the total energy
! should be proportional to sqrt(NIONS)*kBT. This fact indicates that
! SIGtot/sqrt(NIONS)=kBT=const. Then, const.=
! (SIGatom*NIONS)/sqrt(NIONS)=SIGatom*sqrt(NIONS), where SIGatom
! corresponds to the FFM%SIG.
! But I believe that this treatment should be done only for energy.
! For force and stress tensors, I believe that their variances do not
! depend on the system size. Hence, I do not scale them.
                   ! Variance in energies per atom in training data (eV atom^-1)
                   WRITE(FF%IU2) FFM%SIG(1)*EUNIT
                   ! Variance in forces in training data (eV Angst^-1)
                   WRITE(FF%IU2) FFM%SIG(2)*FUNIT,FFM%SIG(3)*FUNIT,FFM%SIG(4)*FUNIT
                   ! Variance in stress tensors in training data (kB)
                   WRITE(FF%IU2) FFM%SIG(5)*SUNIT,FFM%SIG(6)*SUNIT,FFM%SIG(7)*SUNIT 
                   WRITE(FF%IU2) FFM%SIG(8)*SUNIT,FFM%SIG(9)*SUNIT,FFM%SIG(10)*SUNIT
                ENDIF
! Merge FFM%CMAT if necessary
#ifdef scaLAPACK
                ! Covariance matrix (-)
                CALL SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_CMAT_FFM,FFM%CMAT(1:,1:),PAR_SUP_HANDLE%CONTEXT_WORLD)
                IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                   DO IB=1, FFM%NBTOT
                      DO ICOL=1, NCOL, 1
                         IROW=NROW*(ICOL-1)
                         WRITE(FF%IU2) (PAR_SUP_HANDLE%GENERAL_SINGLE(IVAR,IB),IVAR=IROW+1,IROW+NROW)
                      ENDDO
                      IF(MOD(FFM%NBTOT,NROW) .GT. 0) THEN
                         WRITE(FF%IU2) (PAR_SUP_HANDLE%GENERAL_SINGLE(IVAR,IB),IVAR=NCOL*NROW+1,FFM%NBTOT)
                      ENDIF
                   ENDDO
                ENDIF
                CALL FINIT_SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE)
#else
                IF(PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
                   DO IB=1, FFM%NBTOT
                      DO ICOL=1, NCOL, 1
                         IROW=NROW*(ICOL-1)
                         WRITE(FF%IU2) (FFM%CMAT(IVAR,IB),IVAR=IROW+1,IROW+NROW)
                      ENDDO
                      IF(MOD(FFM%NBTOT,NROW) .GT. 0) THEN
                         WRITE(FF%IU2) (FFM%CMAT(IVAR,IB),IVAR=NCOL*NROW+1,FFM%NBTOT)
                      ENDIF
                   ENDDO
                ENDIF
#endif
             ENDIF
             PROFILING_STOP('out_ffncar')
          ENDIF
! Close file
          CLOSE(FF%IU2)
        END SUBROUTINE OUT_FFNCAR

#endif


#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Output histogram
!****************************************************************************************************

        SUBROUTINE OUT_HISTOGRAM (COMM_WORLD,AV,DIR_APP,DIR_LEN,HISCAR,IU4,NX,SIG,X,Y)
          USE ML_FF_STRUCT, ONLY : ML_MPI_PAR
          IMPLICIT NONE
! Input variables
          TYPE (ML_MPI_PAR)                :: COMM_WORLD
          REAL(q)          , INTENT(IN) :: AV(:) !(1:10)
          CHARACTER(LEN=10), INTENT(IN) ::  DIR_APP
          INTEGER          , INTENT(IN) ::  DIR_LEN
          CHARACTER (LEN=10), INTENT(IN) :: HISCAR
          INTEGER          , INTENT(IN) :: IU4
          INTEGER          , INTENT(IN) :: NX
          REAL(q)          , INTENT(IN) :: SIG(:) !(1:10)
          REAL(q)          , INTENT(IN) :: X(:) !(1:(NX+1))
          REAL(q)          , INTENT(IN) :: Y(:,:) !(1:(NX+1),1:10)
! Local variables
          INTEGER                       :: I
          INTEGER                       :: IERR
          INTEGER                       :: IX
          IF(COMM_WORLD%NODE_ME.EQ.0) THEN
! Open file
             OPEN(UNIT=IU4,FILE=DIR_APP(1:DIR_LEN)// &
                 HISCAR,STATUS='UNKNOWN',IOSTAT=IERR)
! Check if error occured in opening file
             IF (IERR.NE.0) THEN
                CALL ERROR_PRINT2TERMINAL("Error in opening ML_HIS file")
             ENDIF
! Output the calculated histograms
             WRITE(IU4,1) AV(1)*EUNIT
             WRITE(IU4,2) AV(2)*FUNIT,AV(3)*FUNIT,AV(4)*FUNIT
             WRITE(IU4,3) AV(5)*SUNIT/10.0_q,AV(8)*SUNIT/10.0_q,AV(10)*SUNIT/10.0_q, &
                          AV(6)*SUNIT/10.0_q,AV(7)*SUNIT/10.0_q,AV(9) *SUNIT/10.0_q
             WRITE(IU4,4) SIG(1)*EUNIT
             WRITE(IU4,5) SIG(2)*FUNIT,SIG(3)*FUNIT,SIG(4)*FUNIT
             WRITE(IU4,6) SIG(5)*SUNIT/10.0_q,SIG(8)*SUNIT/10.0_q,SIG(10)*SUNIT/10.0_q, &
                          SIG(6)*SUNIT/10.0_q,SIG(7)*SUNIT/10.0_q,SIG(9) *SUNIT/10.0_q
             WRITE(IU4,7)
             DO IX=1, NX
                WRITE(IU4,8) X(IX),Y(IX,1),(Y(IX,I),I=2, 4),Y(IX,5), &
                      Y(IX,8),Y(IX,10),Y(IX,6),Y(IX,7),Y(IX,9)
             ENDDO
! Close file
             CLOSE(IU4)
          ENDIF
! Format
1         FORMAT(100("*")/5X,'Average of energy (eV atom^-1)' &
                              /100("-")/5X,F16.8)
2         FORMAT(100("*")/5X,'Average of force (X Y Z) (eV ang^-1 &
                              atom^-1)'/100("-")/5X,F16.8,F16.8,F16.8)
3         FORMAT(100("*")/5X,'Average of stress (XX YY ZZ XY XZ YZ) &
                        (GPa atom^-1)'/100("-")/5X,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8)
4         FORMAT(100("*")/5X,'Variance of energy (eV atom^-1)' &
                              /100("-")/5X,F16.8)
5         FORMAT(100("*")/5X,'Variance of force (X Y Z) (eV ang^-1 &
                              atom^-1)'/100("-")/5X,F16.8,F16.8,F16.8)
6         FORMAT(100("*")/5X,'Variance of stress (XX YY ZZ XY XZ YZ) &
                  (GPa atom^-1)'/100("-")/5X,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8)
7         FORMAT(100("*")/5X,'Normalized histogram'/ &
                          5X,'X, HIS_E, HIS_F, HIS_S'/100("-"))
8         FORMAT(F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8, &
                F16.8,F16.8,F16.8,F16.8)
        END SUBROUTINE OUT_HISTOGRAM
#endif

!**********************************************************************
! Output local atomic energies
!***********************************************************************

        SUBROUTINE OUT_EATOM (COMM_WORLD,A,EATOM_KIN,EATOM_POT,IU,LADD,LNIONS,MIONS,MITYP,MTYP,NIONS,NITYP,NSTEP,NTYP,POSION,TYPE)
           USE ML_FF_STRUCT, ONLY : ML_MPI_PAR
           IMPLICIT NONE
! Global variables.
           TYPE (ML_MPI_PAR)               :: COMM_WORLD
           REAL(q)         , INTENT(IN) :: A(:,:) !(1:3,1:3)
           REAL(q)         , INTENT(IN) :: EATOM_KIN(:) !(1:MIONS)
           REAL(q)         , INTENT(IN) :: EATOM_POT(:) !(1:MIONS)
           INTEGER         , INTENT(IN) :: IU
           INTEGER         , INTENT(IN) :: LADD(:) !(1:MTYP)
           INTEGER         , INTENT(IN) :: LNIONS(:,:) !(1:MITYP,1:MTYP)
           INTEGER         , INTENT(IN) :: MIONS
           INTEGER         , INTENT(IN) :: MITYP
           INTEGER         , INTENT(IN) :: MTYP
           INTEGER         , INTENT(IN) :: NIONS
           INTEGER         , INTENT(IN) :: NITYP(:) !(1:MTYP)
           INTEGER         , INTENT(IN) :: NSTEP
           INTEGER         , INTENT(IN) :: NTYP
           REAL(q)         , INTENT(IN) :: POSION(:,:) !(1:3,1:MIONS)
           CHARACTER(LEN=2), INTENT(IN) :: TYPE(:) !(1:MTYP)
! Local variables.
           INTEGER                      :: INIONS
           INTEGER                      :: INITYP
           INTEGER                      :: INTYP
           INTEGER                      :: INTYP0
           INTEGER                      :: IXYZ
           INTEGER                      :: JXYZ
           IF(COMM_WORLD%NODE_ME.EQ.0) THEN
              WRITE(IU,1) NSTEP
              WRITE(IU,2)
              DO IXYZ=1, 3
                 WRITE(IU,3) (A(JXYZ,IXYZ)*AUTOA,JXYZ=1, 3)
              ENDDO
              WRITE(IU,4)
              DO INTYP=1, NTYP
                 INTYP0=LADD(INTYP)
                 DO INITYP=1, NITYP(INTYP0)
                    INIONS=LNIONS(INITYP,INTYP0)
                    WRITE(IU,5) TYPE(INTYP0),(POSION(IXYZ,INIONS)*AUTOA,IXYZ=1, 3),EATOM_KIN(INIONS)*EUNIT,EATOM_POT(INIONS)*EUNIT
                 ENDDO
              ENDDO
           ENDIF
1          FORMAT (100("*")/5X,'NSTEP=',1X,I10)
2          FORMAT (100("=")/5X,'Lattice constant (Angst)'/100("-"))
3          FORMAT (5X,F16.8,1X,F16.8,1X,F16.8)
4          FORMAT (100("=")/5X,'Atomic positions (Angst) and its & 
                  energy (eV)'/100("-"))
5          FORMAT (5X,A2,1X,F16.8,F16.8,F16.8,F16.8,F16.8)
        END SUBROUTINE OUT_EATOM

!****************************************************************************************************
! Output the heat flux.
!****************************************************************************************************

        SUBROUTINE OUT_HEAT (COMM_WORLD,IU,NSTEP,QHEAT)
           USE ML_FF_STRUCT, ONLY : ML_MPI_PAR
           IMPLICIT NONE
! Input variables.
           TYPE (ML_MPI_PAR)               :: COMM_WORLD
           INTEGER         , INTENT(IN) :: IU
           INTEGER         , INTENT(IN) :: NSTEP
           REAL(q)         , INTENT(IN) :: QHEAT(:) !(1:3)
! Local variables.
           INTEGER                      :: IXYZ
           IF(COMM_WORLD%NODE_ME.EQ.0) THEN
              WRITE(IU,1) NSTEP,(QHEAT(IXYZ)*EUNIT*AUTOA/TUNIT,IXYZ=1, 3)
           ENDIF
1          FORMAT (1X,'NSTEP=',I10,1X,'QXYZ=',E16.8e2,E16.8e2,E16.8e2)
        END SUBROUTINE OUT_HEAT

      END MODULE IOHANDLE
#endif
