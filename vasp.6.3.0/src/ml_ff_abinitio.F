#if defined(MPI) || defined(MPI_CHAIN)
#include "symbol.inc"
!**********************************************************************
! Module for handling ab initio data
!**********************************************************************

      MODULE ABINITIO_DATA
        USE ML_FF_CONSTANT
        USE LOGFILE
        USE IOHANDLE
        USE LATTICE_ML
        USE MEMORY_USAGE
        USE MPI_DATA
        USE ML_FF_PREC
        USE PRINT2TERMINAL
        USE TIME_CHECK
        
        IMPLICIT NONE

        CONTAINS

!****************************************************************************************************
! Initialize variables necessary for storing new ab initio data.
!****************************************************************************************************
        SUBROUTINE SET_ABN_PAR_ISTART3(FF,AB,ABN,MLH_TMP,MEM1,ICONF,NITYP_NEW,LALLOCATE_ARRAYS)
           USE ML_FF_STRUCT, ONLY: FF_PAR, ABINITIO, ABINITIONEW, ML_SUPER_TYPE, MEMORY, ML_IO_WRITE
           IMPLICIT NONE
           ! Input and output variables
           TYPE (FF_PAR)         :: FF
           TYPE (ABINITIO)       :: AB
           TYPE (ABINITIONEW)    :: ABN
           TYPE (ML_SUPER_TYPE)  :: MLH_TMP
           TYPE (MEMORY)         :: MEM1
           ! Input variables
           INTEGER, INTENT(IN)   :: ICONF
           LOGICAL, INTENT(IN)   :: LALLOCATE_ARRAYS
           INTEGER, INTENT(IN)   :: NITYP_NEW(:)
           ! Local variables
           INTEGER               :: INIONS
           INTEGER               :: INITYP
           INTEGER               :: INTYP
           INTEGER               :: ISYS
           ABN%SZNAM2=MLH_TMP%AB%SZNAM2(ICONF)
           ABN%NSYS=FF%NSYS           
           ABN%MCONF=FF%MCONF_NEW
           ABN%MITYP=NITYP_NEW(1)
           DO INTYP=2, MLH_TMP%AB%NTYP(ICONF)
              ABN%MITYP=MAX(ABN%MITYP,NITYP_NEW(INTYP))
           ENDDO
           ABN%MIONS=MLH_TMP%AB%NIONS(ICONF)
           ABN%MTYP=AB%MTYP
           ABN%NCONF=1
           ABN%NTYP=MLH_TMP%AB%NTYP(ICONF)
           ABN%NIONS=MLH_TMP%AB%NIONS(ICONF)
           ! Allocate arrays used for storing new ab initio data.
           IF (LALLOCATE_ARRAYS) THEN
              CALL ALLOCATE_ABN_ARRAY (MEM1,ABN)
           ENDIF
! Configure the system name.
           DO ISYS=1, FF%NSYS
              ABN%SZNAM2_ADD(ISYS)=FF%SZNAM2(ISYS)
           ENDDO
           DO ISYS=1, ABN%NSYS
              IF(ABN%SZNAM2_ADD(ISYS).EQ.ABN%SZNAM2) THEN
                 ABN%ISYS=ISYS
                 EXIT
              ENDIF
           ENDDO
! Configure element types, the number of atoms and a list storing order of element type.
           DO INTYP=1, ABN%MTYP_ALLOCATE
              ABN%LADD(INTYP)=MLH_TMP%AB%LADD(INTYP,ICONF)
              ABN%NITYP(INTYP)=MLH_TMP%AB%NITYP(INTYP,ICONF)
              ABN%POMASS(INTYP)=MLH_TMP%AB%POMASS(INTYP)
              ABN%TYPE(INTYP)=MLH_TMP%AB%TYPE(INTYP)
           ENDDO
! Configure list of atoms and element type.
           DO INITYP=1, ABN%MITYP_ALLOCATE
              DO INTYP=1, ABN%MTYP_ALLOCATE
                 ABN%LNIONS(INITYP,INTYP)=MLH_TMP%AB%LNIONS(INITYP,INTYP,ICONF)
              ENDDO
           ENDDO
           DO INIONS=1, ABN%MIONS_ALLOCATE
              ABN%LNITYP(INIONS)=MLH_TMP%AB%LNITYP(INIONS,ICONF)
              ABN%LNTYP(INIONS)=MLH_TMP%AB%LNTYP(INIONS,ICONF)
           ENDDO
        END SUBROUTINE SET_ABN_PAR_ISTART3

        SUBROUTINE SET_ABN_PAR(ABN,MEM1,ISTART_FF,MCONF_NEW_FF,MTYP_FF,NIONS_NEW,NITYP_NEW,NSYS_FF, &
                   NTYP_NEW,POMASS_NEW,SZNAM2_FF,SZNAM2_NEW,TYPE_FF,TYPE_NEW,LALLOCATE_ARRAYS)
           USE ML_FF_STRUCT, ONLY: ABINITIONEW, MEMORY, ML_IO_WRITE
           IMPLICIT NONE
! Input variables
           TYPE (ABINITIONEW)             :: ABN
           TYPE (MEMORY)                  :: MEM1
           INTEGER, INTENT(IN)            :: ISTART_FF
           LOGICAL, INTENT(IN)            :: LALLOCATE_ARRAYS
           INTEGER, INTENT(IN)            :: MCONF_NEW_FF
           INTEGER, INTENT(IN)            :: MTYP_FF
           INTEGER, INTENT(IN)            :: NIONS_NEW
           INTEGER, INTENT(IN)            :: NITYP_NEW(:)
           INTEGER, INTENT(IN)            :: NSYS_FF
           INTEGER, INTENT(IN)            :: NTYP_NEW
           REAL(q), INTENT(IN)            :: POMASS_NEW(:)
           CHARACTER (LEN=40), INTENT(IN) :: SZNAM2_FF(:)
           CHARACTER (LEN=40), INTENT(IN) :: SZNAM2_NEW
           CHARACTER (LEN=2) , INTENT(IN) :: TYPE_FF(:)
           CHARACTER (LEN=2) , INTENT(IN) :: TYPE_NEW(:)
! Local variables
           INTEGER                        :: INIONS
           INTEGER                        :: INITYP
           INTEGER                        :: INTYP
           INTEGER                        :: INTYP0
           INTEGER                        :: ISYS
           INTEGER                        :: JNTYP0
           INTEGER                        :: MITYP_NEW
! Copy the VASP parameters to parameters for new ab initio data.
           ABN%SZNAM2=SZNAM2_NEW
           ABN%NSYS=NSYS_FF
           ABN%MCONF=MCONF_NEW_FF
           MITYP_NEW=NITYP_NEW(1)
           DO INTYP=2, NTYP_NEW
              MITYP_NEW=MAX(MITYP_NEW,NITYP_NEW(INTYP))
           ENDDO
           ABN%MITYP=MITYP_NEW
           ABN%MIONS=NIONS_NEW
           ABN%MTYP=MTYP_FF
           ABN%NCONF=1
           ABN%NTYP=NTYP_NEW
           ABN%NIONS=NIONS_NEW
! We need to set sizes for array allocation
           IF (ISTART_FF.NE.3) THEN
              ABN%MTYP_ALLOCATE=ABN%MTYP
              ABN%MIONS_ALLOCATE=ABN%MIONS
              ABN%MITYP_ALLOCATE=ABN%MITYP
              ABN%NSYS_ALLOCATE=ABN%NSYS
           ENDIF
! Allocate arrays used for storing new ab initio data.
           IF (LALLOCATE_ARRAYS) THEN
              CALL ALLOCATE_ABN_ARRAY (MEM1,ABN)
           ENDIF
! Configure the system name.
           DO ISYS=1, NSYS_FF
              ABN%SZNAM2_ADD(ISYS)=SZNAM2_FF(ISYS)
           ENDDO
           DO ISYS=1, ABN%NSYS
              IF(ABN%SZNAM2_ADD(ISYS).EQ.ABN%SZNAM2) THEN
                 ABN%ISYS=ISYS
                 EXIT
              ENDIF
           ENDDO
! Configure element types, the number of atoms and a list storing order of element type.
           DO INTYP=1, ABN%NTYP
              DO JNTYP0=1, MTYP_FF
                 IF(TYPE_NEW(INTYP) .EQ. TYPE_FF(JNTYP0)) THEN
                    ABN%LADD(INTYP)=JNTYP0
                    ABN%NITYP(JNTYP0)=NITYP_NEW(INTYP)
!unit                    ABN%POMASS(JNTYP0)=POMASS_NEW(INTYP)/MUNIT
                    ABN%POMASS(JNTYP0)=POMASS_NEW(INTYP)
                    ABN%TYPE(JNTYP0)=TYPE_NEW(INTYP)
                    EXIT
                 ENDIF
              ENDDO
           ENDDO
! Configure list of atoms and element type.
           INIONS=0
           DO INTYP=1, ABN%NTYP
              INTYP0=ABN%LADD(INTYP)
              DO INITYP=1, ABN%NITYP(INTYP0)
                 INIONS=INIONS+1
                 ABN%LNIONS(INITYP,INTYP0)=INIONS
                 ABN%LNITYP(INIONS)=INITYP
                 ABN%LNTYP(INIONS)=INTYP0
              ENDDO
           ENDDO
        END SUBROUTINE SET_ABN_PAR

!****************************************************************************************************
! Subroutine to set parameters for storing ab initio data
!****************************************************************************************************


        SUBROUTINE SET_AB_PAR (ML_INPUT_HANDLE,FF,TIM,COMM_WORLD,MEM1,LOGF,AB,ABN,CTIFOR_FF,IFF_FF,ISTART_FF, &
                   LB_FFM,LBHEAD_FFM,LFMAT_FF,LGENFF_FF,LLIST_FF,LLIST_NEW_FF,LMLONLY_FF,LNCONF_FF,LNIONS_FF, &
                   LOVER_FF,LPS_FF,LSIG_FF,LSOAP_FF,MB_FFM,MBTOT_FFM,MCONF_FF, &
                   MTYP_FF,NB_FFM,NBTOT_FFM,NITYP_NEW,NIONS_NEW,NSYS_FF,NTYP_NEW, &
                   POMASS_FF,SZNAM2_FF,TYPE_FF)
          USE ML_FF_STRUCT, ONLY: ML_FF_INPUT, FF_PAR, TIME, ML_MPI_PAR, MEMORY, LOG_PAR, ABINITIO, ABINITIONEW, ML_IO_WRITE
          IMPLICIT NONE
! Input variables
          TYPE (ML_FF_INPUT)             :: ML_INPUT_HANDLE
          TYPE (FF_PAR)                  :: FF
          TYPE (TIME)                    :: TIM
          TYPE(ML_MPI_PAR)               :: COMM_WORLD
          TYPE (MEMORY)                  :: MEM1
          TYPE (LOG_PAR)                 :: LOGF
          TYPE (ABINITIO)                :: AB
          TYPE (ABINITIONEW)             :: ABN
          REAL(q), INTENT(INOUT)         :: CTIFOR_FF
          LOGICAL, INTENT(INOUT)         :: IFF_FF
          INTEGER, INTENT(IN)            :: ISTART_FF
          INTEGER, INTENT(INOUT)         :: LB_FFM(:,:) !(1:MB_FFM,1:MTYP_FF)
          INTEGER, INTENT(INOUT)         :: LBHEAD_FFM(:) !(1:MTYP_FF)
          LOGICAL, INTENT(INOUT)         :: LFMAT_FF
          LOGICAL, INTENT(INOUT)         :: LGENFF_FF
          LOGICAL, INTENT(INOUT)         :: LLIST_FF
          LOGICAL, INTENT(INOUT)         :: LLIST_NEW_FF
          LOGICAL, INTENT(IN)            :: LMLONLY_FF
          INTEGER, INTENT(INOUT)         :: LNCONF_FF(:,:) !(1:MB_FFM,1:MTYP_FF)
          INTEGER, INTENT(INOUT)         :: LNIONS_FF(:,:) !(1:MB_FFM,1:MTYP_FF)
          LOGICAL, INTENT(INOUT)         :: LOVER_FF
          LOGICAL, INTENT(INOUT)         :: LPS_FF
          LOGICAL, INTENT(INOUT)         :: LSIG_FF
          LOGICAL, INTENT(INOUT)         :: LSOAP_FF
          INTEGER, INTENT(IN)            :: MB_FFM
          INTEGER, INTENT(IN)            :: MBTOT_FFM
          INTEGER, INTENT(IN)            :: MCONF_FF
          INTEGER, INTENT(IN)            :: MTYP_FF
          INTEGER, INTENT(INOUT)         :: NB_FFM(:) !(1:MTYP_FF)
          INTEGER, INTENT(INOUT)         :: NBTOT_FFM
          INTEGER, INTENT(IN)            :: NITYP_NEW(:) !(1:NTYP_NEW)
          INTEGER, INTENT(IN)            :: NIONS_NEW
          INTEGER, INTENT(IN)            :: NSYS_FF
          INTEGER, INTENT(IN)            :: NTYP_NEW
          REAL(q), INTENT(IN)            :: POMASS_FF(:) !(1:MTYP_FF)
          CHARACTER (LEN=40), INTENT(IN) :: SZNAM2_FF(:) !(1:NSYS_FF)
          CHARACTER (LEN=2) , INTENT(IN) :: TYPE_FF(:) !(1:MTYP_FF)
! Local variables
          REAL(q), ALLOCATABLE           :: A_TMP(:,:,:)
          REAL(q), ALLOCATABLE           :: B_TMP(:,:,:)
          REAL(q)                        :: CPUTIM0
          REAL(q)                        :: CPUTIM1
          INTEGER, ALLOCATABLE           :: ISYS_TMP(:)
          INTEGER, ALLOCATABLE           :: LADD_TMP(:,:)
          INTEGER, ALLOCATABLE           :: LB_TMP(:,:)
          INTEGER, ALLOCATABLE           :: LNIONS_TMP(:,:,:)
          INTEGER, ALLOCATABLE           :: LNITYP_TMP(:,:)
          INTEGER, ALLOCATABLE           :: LNTYP_TMP(:,:)
          INTEGER, ALLOCATABLE           :: LSYS_TMP(:)
          INTEGER, ALLOCATABLE           :: MIONS_NEW
          INTEGER                        :: MITYP_NEW
          INTEGER, ALLOCATABLE           :: NBMAX_TMP(:)
          INTEGER, ALLOCATABLE           :: NIONS_TMP(:)
          INTEGER, ALLOCATABLE           :: NITYP_TMP(:,:)
          INTEGER, ALLOCATABLE           :: NTYP_TMP(:)
          REAL(q), ALLOCATABLE           :: OMEGA_TMP(:)
          REAL(q), ALLOCATABLE           :: POSION_TMP(:,:,:)
          REAL(q)                        :: RMEM_HELP
          CHARACTER(LEN=40), ALLOCATABLE :: SZNAM2_TMP(:)
          REAL(q), ALLOCATABLE           :: TOTEN_TMP(:)
          REAL(q), ALLOCATABLE           :: TIFOR_TMP(:,:,:)
          REAL(q), ALLOCATABLE           :: TSIF_TMP(:,:,:)
          REAL(q)                        :: VPUTIM0
          REAL(q)                        :: VPUTIM1
! Start profiling
          PROFILING_START('set_ab_par')
! Measure the initial times.
          CALL VTIME (VPUTIM0,CPUTIM0)
#ifndef ML_LEARNING_NOT_ACTIVE
! Allocate temporary an array.
          CALL ALLOCATE_TMP_ARRAY1
#endif

! Initialize or read variables necessary for storing ab initio data.
          IF (ISTART_FF .EQ. 0) THEN
#ifndef ML_LEARNING_NOT_ACTIVE
             CALL INIT_AB_ISTART0
#endif
          ELSE IF (ISTART_FF.EQ.1) THEN
#ifndef ML_LEARNING_NOT_ACTIVE
             CALL INIT_AB_ISTART1
#endif
          ELSE IF (ISTART_FF.EQ.2) THEN
             CALL INIT_AB_ISTART2
          ELSE IF (ISTART_FF.EQ.3) THEN
#ifndef ML_LEARNING_NOT_ACTIVE
             CALL INIT_AB_ISTART3
#endif
          ENDIF
#ifndef ML_LEARNING_NOT_ACTIVE
          IF(ISTART_FF /= 2) THEN
             CALL LOGSECTION_AB_INITIO_DATA(LOGF, AB, ABN)
          END IF
#endif
! Measure the final times and add them
          CALL VTIME (VPUTIM1,CPUTIM1)
          TIM%CPUTIM(1)=TIM%CPUTIM(1)+CPUTIM1-CPUTIM0
          TIM%VPUTIM(1)=TIM%VPUTIM(1)+VPUTIM1-VPUTIM0
#ifndef ML_LEARNING_NOT_ACTIVE
! Deallocate a temporary array
          CALL DEALLOCATE_TMP_ARRAY1
#endif
! Stop profiling
          PROFILING_STOP('set_ab_par')

          CONTAINS

#ifndef ML_LEARNING_NOT_ACTIVE
! Initialize variables storing ab initio data.

          SUBROUTINE INIT_AB_ISTART0
            IMPLICIT NONE
! Local variables
            INTEGER                        :: INTYP
            INTEGER                        :: INTYP0
            INTEGER                        :: ISYS
! The parameters were already set in SET_FF_PAR subroutine in ff.F.
! But for some of them are also set, here.
            LFMAT_FF     =.FALSE.
            LGENFF_FF    =.FALSE.
            LLIST_FF     =.TRUE.
            LLIST_NEW_FF =.TRUE.
            LOVER_FF     =.FALSE.
            LPS_FF       =.FALSE.
            LSIG_FF      =.FALSE.
            LSOAP_FF     =.FALSE.
            AB%NSYS=NSYS_FF
            AB%MIONS=NIONS_NEW
            AB%MTYP=NTYP_NEW
            MITYP_NEW=NITYP_NEW(1)
            DO INTYP=2, NTYP_NEW
               MITYP_NEW=MAX(MITYP_NEW,NITYP_NEW(INTYP))
            ENDDO
            AB%MITYP=MITYP_NEW
            AB%MCONF=MCONF_FF
            AB%NCONF=0
! Set LCTIFOR_USE_FROM_FILE to .FALSE. when ML_FF_ISTART=0
            FF%LCTIFOR_USE_FROM_FILE=.FALSE.
! If one starts from scratch
            FF%LCTIFOR_WRITE_FROM_FILE=.TRUE.
! We need to set array sizes for allocation
            IF (ISTART_FF.LT.3) THEN
               AB%MCONF_ALLOCATE=AB%MCONF
               AB%MIONS_ALLOCATE=AB%MIONS
               AB%MITYP_ALLOCATE=AB%MITYP
               AB%MTYP_ALLOCATE=AB%MTYP
               AB%NSYS_ALLOCATE=AB%NSYS
            ENDIF
            CALL ALLOCATE_AB_ARRAY (MEM1,AB)
            DO ISYS=1, NSYS_FF
               AB%SZNAM2_ADD(ISYS)=SZNAM2_FF(ISYS)
            ENDDO
            DO INTYP0=1, AB%MTYP
               AB%POMASS(INTYP0)=POMASS_FF(INTYP0)
               AB%TYPE(INTYP0)=TYPE_FF(INTYP0)
            ENDDO
            NBMAX_TMP=0
            LB_FFM=0
            LBHEAD_FFM=0
            NB_FFM=0
            NBTOT_FFM=0
          END SUBROUTINE INIT_AB_ISTART0
#endif

#ifndef ML_LEARNING_NOT_ACTIVE
! Initialization of ab initio parameters for ISTART_FF=1.

          SUBROUTINE INIT_AB_ISTART1
             IMPLICIT NONE
             INTEGER :: INTYP
             INTEGER :: INTYP0
             INTEGER :: ICONF
             INTEGER :: ISYS
             INTEGER :: IFLAG
             MITYP_NEW=NITYP_NEW(1)
             DO INTYP=2, NTYP_NEW
                MITYP_NEW=MAX(MITYP_NEW,NITYP_NEW(INTYP))
             ENDDO
             CALL READ_ABCAR_SIZE2 (ML_INPUT_HANDLE,AB,FF,MITYP_NEW)
! We need to set array sizes for allocation
             AB%MCONF_ALLOCATE=AB%MCONF
             AB%MIONS_ALLOCATE=AB%MIONS
             AB%MITYP_ALLOCATE=AB%MITYP
             AB%MTYP_ALLOCATE=AB%MTYP
             AB%NSYS_ALLOCATE=AB%NSYS
             CALL ALLOCATE_AB_ARRAY (MEM1,AB)
             CALL READ_ABCAR (AB,FF,AB%A,AB%ABCAR,AB%B,AB%EATOM,ISTART_FF,AB%IU1,AB%ISYS, &
                  AB%LADD,AB%LB,LB_FFM,LBHEAD_FFM, &
                  LNCONF_FF,AB%LNIONS,LNIONS_FF,AB%LNITYP,AB%LNTYP,AB%LSYS, &
                  AB%MTYP,MTYP_FF, &
                  NB_FFM,NBMAX_TMP,NBTOT_FFM,AB%NCONF,AB%NCONF_SYS,AB%NIONS,AB%NITYP,ABN%NITYP, &
                  AB%NSYS,NSYS_FF,AB%NTYP,NTYP_NEW,AB%OMEGA,AB%POMASS,POMASS_FF,AB%POSION, &
                  AB%SZNAM2,AB%SZNAM2_ADD,SZNAM2_FF,AB%TIFOR,AB%TOTEN,AB%TSIF,AB%TYPE,TYPE_FF)
! Copy CTIFORs from ML_AB file to array holding the old and new errors
             IF (FF%LCTIFOR_USE_FROM_FILE) THEN
                AB%CTIFOR_ALLCONF_NEW(1:AB%NCONF)=AB%CTIFOR_ALLCONF_OLD(1:AB%NCONF)
             ENDIF
! CTIFOR
             IF (FF%LCTIFOR_READ_FROM_FILE.AND.FF%LCTIFOR_USE_FROM_FILE) THEN
                FF%LCTIFOR_WRITE_FROM_FILE=.TRUE.
             ENDIF
! Set some important control variables here
             LFMAT_FF=.TRUE.
             LLIST_FF=.TRUE.
             LLIST_NEW_FF=.TRUE.
             LPS_FF=.TRUE.
             LSIG_FF=.TRUE.
             LSOAP_FF=.TRUE.
! Check the availability of ab initio data for force field generation.
! If the number of any basis set equals zero, the machine cannot generate the force field.
! In this situation, the machine does not execute the force field generation and the sparsification
! necessary for the force field generation.
             LGENFF_FF=.TRUE.
             LOVER_FF=.TRUE.
             DO INTYP0=1, MTYP_FF
                IF(NB_FFM(INTYP0).LE.1) THEN
                   LGENFF_FF=.FALSE.
                   LOVER_FF =.FALSE.
                ENDIF
             ENDDO
             CALL CHANGE_ORDER_AB
          END SUBROUTINE INIT_AB_ISTART1
#endif

#ifndef ML_LEARNING_NOT_ACTIVE
! Change the order of ab initio data to ease the restructuring of
! training data in sparsification process

          SUBROUTINE CHANGE_ORDER_AB
            IMPLICIT NONE
! Local variables
            INTEGER                        :: I
            INTEGER                        :: ICONF
            INTEGER                        :: INIONS
            INTEGER                        :: INTYP0
            INTEGER                        :: NCONF_TMP
            CALL ALLOCATE_TMP_ARRAY2
            NCONF_TMP=0
            DO ICONF=1, AB%NCONF
               IF(AB%LSYS(ICONF).EQ.0) THEN
                  NCONF_TMP=NCONF_TMP+1
                  A_TMP(1:3,1:3,NCONF_TMP)=AB%A(1:3,1:3,ICONF)
                  B_TMP(1:3,1:3,NCONF_TMP)=AB%B(1:3,1:3,ICONF)
                  ISYS_TMP(NCONF_TMP)=AB%ISYS(ICONF)
                  LADD_TMP(1:AB%MTYP,NCONF_TMP)=AB%LADD(1:AB%MTYP,ICONF)
                  LB_TMP(1:AB%MIONS,NCONF_TMP)=AB%LB(1:AB%MIONS,ICONF)
                  LNIONS_TMP(1:AB%MITYP,1:AB%MTYP,NCONF_TMP)=AB%LNIONS(1:AB%MITYP,1:AB%MTYP,ICONF)
                  LNITYP_TMP(1:AB%MIONS,NCONF_TMP)=AB%LNITYP(1:AB%MIONS,ICONF)
                  LNTYP_TMP(1:AB%MIONS,NCONF_TMP)=AB%LNTYP(1:AB%MIONS,ICONF)
                  LSYS_TMP(NCONF_TMP)=AB%LSYS(ICONF)
                  NIONS_TMP(NCONF_TMP)=AB%NIONS(ICONF)
                  NITYP_TMP(1:AB%MTYP,NCONF_TMP)=AB%NITYP(1:AB%MTYP,ICONF)
                  NTYP_TMP(NCONF_TMP)=AB%NTYP(ICONF)
                  OMEGA_TMP(NCONF_TMP)=AB%OMEGA(ICONF)
                  POSION_TMP(1:3,1:AB%MIONS,NCONF_TMP)=AB%POSION(1:3,1:AB%MIONS,ICONF)
                  SZNAM2_TMP(NCONF_TMP)=AB%SZNAM2(ICONF)
                  TOTEN_TMP(NCONF_TMP)=AB%TOTEN(ICONF)
                  TIFOR_TMP(1:3,1:AB%MIONS,NCONF_TMP)=AB%TIFOR(1:3,1:AB%MIONS,ICONF)
                  TSIF_TMP(1:3,1:3,NCONF_TMP)=AB%TSIF(1:3,1:3,ICONF)
               ENDIF
            ENDDO
            DO ICONF=1, AB%NCONF
               IF(AB%LSYS(ICONF).EQ.1) THEN
                  NCONF_TMP=NCONF_TMP+1
                  A_TMP(1:3,1:3,NCONF_TMP)=AB%A(1:3,1:3,ICONF)
                  B_TMP(1:3,1:3,NCONF_TMP)=AB%B(1:3,1:3,ICONF)
                  ISYS_TMP(NCONF_TMP)=AB%ISYS(ICONF)
                  LADD_TMP(1:AB%MTYP,NCONF_TMP)=AB%LADD(1:AB%MTYP,ICONF)
                  LB_TMP(1:AB%MIONS,NCONF_TMP)=AB%LB(1:AB%MIONS,ICONF)
                  LNIONS_TMP(1:AB%MITYP,1:AB%MTYP,NCONF_TMP)=AB%LNIONS(1:AB%MITYP,1:AB%MTYP,ICONF)
                  LNITYP_TMP(1:AB%MIONS,NCONF_TMP)=AB%LNITYP(1:AB%MIONS,ICONF)
                  LNTYP_TMP(1:AB%MIONS,NCONF_TMP)=AB%LNTYP(1:AB%MIONS,ICONF)
                  LSYS_TMP(NCONF_TMP)=AB%LSYS(ICONF)
                  NIONS_TMP(NCONF_TMP)=AB%NIONS(ICONF)
                  NITYP_TMP(1:AB%MTYP,NCONF_TMP)=AB%NITYP(1:AB%MTYP,ICONF)
                  NTYP_TMP(NCONF_TMP)=AB%NTYP(ICONF)
                  OMEGA_TMP(NCONF_TMP)=AB%OMEGA(ICONF)
                  POSION_TMP(1:3,1:AB%MIONS,NCONF_TMP)=AB%POSION(1:3,1:AB%MIONS,ICONF)
                  SZNAM2_TMP(NCONF_TMP)=AB%SZNAM2(ICONF)
                  TOTEN_TMP(NCONF_TMP)=AB%TOTEN(ICONF)
                  TIFOR_TMP(1:3,1:AB%MIONS,NCONF_TMP)=AB%TIFOR(1:3,1:AB%MIONS,ICONF)
                  TSIF_TMP(1:3,1:3,NCONF_TMP)=AB%TSIF(1:3,1:3,ICONF)
               ENDIF
            ENDDO
            DO ICONF=1, AB%NCONF
               AB%A(1:3,1:3,ICONF)=A_TMP(1:3,1:3,ICONF)
               AB%B(1:3,1:3,ICONF)=B_TMP(1:3,1:3,ICONF)
               AB%ISYS(ICONF)=ISYS_TMP(ICONF)
               AB%LADD(1:AB%MTYP,ICONF)=LADD_TMP(1:AB%MTYP,ICONF)
               AB%LB(1:AB%MIONS,ICONF)=LB_TMP(1:AB%MIONS,ICONF)
               AB%LNIONS(1:AB%MITYP,1:AB%MTYP,ICONF)=LNIONS_TMP(1:AB%MITYP,1:AB%MTYP,ICONF)
               AB%LNITYP(1:AB%MIONS,ICONF)=LNITYP_TMP(1:AB%MIONS,ICONF)
               AB%LNTYP(1:AB%MIONS,ICONF)=LNTYP_TMP(1:AB%MIONS,ICONF)
               AB%LSYS(ICONF)=LSYS_TMP(ICONF)
               AB%NIONS(ICONF)=NIONS_TMP(ICONF)
               AB%NITYP(1:AB%MTYP,ICONF)=NITYP_TMP(1:AB%MTYP,ICONF)
               AB%NTYP(ICONF)=NTYP_TMP(ICONF)
               AB%OMEGA(ICONF)=OMEGA_TMP(ICONF)
               AB%POSION(1:3,1:AB%MIONS,ICONF)=POSION_TMP(1:3,1:AB%MIONS,ICONF)
               AB%SZNAM2(ICONF)=SZNAM2_TMP(ICONF)
               AB%TOTEN(ICONF)=TOTEN_TMP(ICONF)
               AB%TIFOR(1:3,1:AB%MIONS,ICONF)=TIFOR_TMP(1:3,1:AB%MIONS,ICONF)
               AB%TSIF(1:3,1:3,ICONF)=TSIF_TMP(1:3,1:3,ICONF)
            ENDDO
! Here, the order of the basis sets is changed. Hence, it should be noted that 
! the order of design matrix elements is also changed from that in the previous calculations
! which provided the ABCAR file.
            LNCONF_FF(1:MB_FFM,1:MTYP_FF)=0
            LNIONS_FF(1:MB_FFM,1:MTYP_FF)=0
            NB_FFM(1:MTYP_FF)=0
            DO ICONF=1, AB%NCONF
               DO INIONS=1, AB%NIONS(ICONF)
                  IF(AB%LB(INIONS,ICONF).NE.0) THEN
                     INTYP0=AB%LNTYP(INIONS,ICONF)
                     NB_FFM(INTYP0)=NB_FFM(INTYP0)+1
                     LNCONF_FF(NB_FFM(INTYP0),INTYP0)=ICONF
                     LNIONS_FF(NB_FFM(INTYP0),INTYP0)=INIONS
                  ENDIF
               ENDDO
            ENDDO
            CALL DEALLOCATE_TMP_ARRAY2
          END SUBROUTINE CHANGE_ORDER_AB
#endif

! Initialization of ab initio parameters for ISTART_FF=2.
          SUBROUTINE INIT_AB_ISTART2
            IMPLICIT NONE
            IFF_FF      =.TRUE.
            LFMAT_FF    =.FALSE.
            LGENFF_FF   =.FALSE.
            LLIST_FF    =.FALSE.
            LLIST_NEW_FF=.TRUE.
            LOVER_FF    =.FALSE.
            LPS_FF      =.FALSE.
            LSIG_FF     =.FALSE.
            LSOAP_FF    =.FALSE.
          END SUBROUTINE INIT_AB_ISTART2

#ifndef ML_LEARNING_NOT_ACTIVE
! Initialization of ab initio parameters for ISTART_FF=1.
          SUBROUTINE INIT_AB_ISTART3
             IMPLICIT NONE
             INTEGER :: INTYP
             INTEGER :: INTYP0
             INTEGER :: ICONF
! We set MITYP_NEW and NIONS_NEW only for dummy reasons here
             MITYP_NEW=0
             CALL READ_ABCAR_SIZE2 (ML_INPUT_HANDLE,AB,FF,MITYP_NEW)
! We need to set array sizes for allocation
             AB%MCONF_ALLOCATE=AB%MCONF
             AB%MIONS_ALLOCATE=AB%MIONS
             AB%MITYP_ALLOCATE=AB%MITYP
             AB%MTYP_ALLOCATE=AB%MTYP
             AB%NSYS_ALLOCATE=AB%NSYS
             CALL ALLOCATE_AB_ARRAY (MEM1,AB)
! Read the data from the ML_AB file
! ABN%NITYP is not used if ISTART_FF=3
             CALL READ_ABCAR (AB,FF,AB%A,AB%ABCAR,AB%B,AB%EATOM,ISTART_FF,AB%IU1,AB%ISYS, &
                  AB%LADD,AB%LB,LB_FFM,LBHEAD_FFM, &
                  LNCONF_FF,AB%LNIONS,LNIONS_FF,AB%LNITYP,AB%LNTYP,AB%LSYS, &
                  AB%MTYP,MTYP_FF, &
                  NB_FFM,NBMAX_TMP,NBTOT_FFM,AB%NCONF,AB%NCONF_SYS,AB%NIONS,AB%NITYP,ABN%NITYP, &
                  AB%NSYS,NSYS_FF,AB%NTYP,NTYP_NEW,AB%OMEGA,AB%POMASS,POMASS_FF,AB%POSION, &
                  AB%SZNAM2,AB%SZNAM2_ADD,SZNAM2_FF,AB%TIFOR,AB%TOTEN,AB%TSIF,AB%TYPE,TYPE_FF)
! Copy CTIFOR for usage
             IF (FF%LCTIFOR_READ_FROM_FILE.AND.FF%LCTIFOR_USE_FROM_FILE) THEN
                AB%CTIFOR_ALLCONF_NEW(1:AB%NCONF)=AB%CTIFOR_ALLCONF_OLD(1:AB%NCONF)
                FF%LCTIFOR_WRITE_FROM_FILE=.TRUE.
             ENDIF
! Set some important control variables here
             LFMAT_FF=.TRUE.
             LLIST_FF=.TRUE.
             LLIST_NEW_FF=.TRUE.
             LPS_FF=.TRUE.
             LSIG_FF=.TRUE.
             LSOAP_FF=.TRUE.
! Check the availability of ab initio data for force field generation.
! If the number of any basis set equals zero, the machine cannot generate the force field.
! In this situation, the machine does not execute the force field generation and the sparsification
! necessary for the force field generation.
             LGENFF_FF=.TRUE.
             LOVER_FF=.TRUE.
             DO INTYP0=1, MTYP_FF
                IF(NB_FFM(INTYP0).LE.1) THEN
                   LGENFF_FF=.FALSE.
                   LOVER_FF =.FALSE.
                ENDIF
             ENDDO
! We set AB%LSYS to 1 everywhere since the new structures are all structures from existing data
             AB%LSYS=1
! MITYP_NEW is set to that of the first structure in the ab initio database
             MITYP_NEW=AB%NITYP(1,1)
             DO INTYP=2,AB%NTYP(1)
                MITYP_NEW=MAX(MITYP_NEW,AB%NITYP(INTYP,1))
             ENDDO
! For ISTART=1 we need not to call CHANGE_ORDER like for ISTART=1
          END SUBROUTINE INIT_AB_ISTART3
#endif

! Allocate temporary arrays

#ifndef ML_LEARNING_NOT_ACTIVE
          SUBROUTINE ALLOCATE_TMP_ARRAY1
            IMPLICIT NONE
            IF(ALLOCATED(NBMAX_TMP)) THEN
               DEALLOCATE(NBMAX_TMP)
            ENDIF
            ALLOCATE (NBMAX_TMP(1:MTYP_FF))
          END SUBROUTINE ALLOCATE_TMP_ARRAY1
#endif

#ifndef ML_LEARNING_NOT_ACTIVE
          SUBROUTINE ALLOCATE_TMP_ARRAY2
            IMPLICIT NONE
            RMEM_HELP=0.0_q
            IF(ALLOCATED(A_TMP)) THEN
               DEALLOCATE(A_TMP)
            ENDIF
            ALLOCATE(A_TMP(1:3,1:3,1:AB%NCONF))
            RMEM_HELP=RMEM_HELP+REAL(AB%NCONF,q)*9*8/1.0E+06_q
            IF(ALLOCATED(B_TMP)) THEN
               DEALLOCATE(B_TMP)
            ENDIF
            ALLOCATE(B_TMP(1:3,1:3,1:AB%NCONF))
            RMEM_HELP=RMEM_HELP+REAL(AB%NCONF,q)*9*8/1.0E+06_q
            IF(ALLOCATED(ISYS_TMP)) THEN
               DEALLOCATE(ISYS_TMP)
            ENDIF
            ALLOCATE(ISYS_TMP(1:AB%NCONF))
            RMEM_HELP=RMEM_HELP+REAL(AB%NCONF,q)/1.0E+06_q
            IF(ALLOCATED(LADD_TMP)) THEN
               DEALLOCATE(LADD_TMP)
            ENDIF
            ALLOCATE(LADD_TMP(1:AB%MTYP,1:AB%NCONF))
            RMEM_HELP=RMEM_HELP+REAL(AB%MTYP,q)*REAL(AB%NCONF,q)/1.0E+06_q
            IF(ALLOCATED(LB_TMP)) THEN
               DEALLOCATE(LB_TMP)
            ENDIF
            ALLOCATE(LB_TMP(1:AB%MIONS,1:AB%NCONF))
            RMEM_HELP=RMEM_HELP+REAL(AB%MIONS,q)*REAL(AB%NCONF,q)/1.0E+06_q
            IF(ALLOCATED(LNIONS_TMP)) THEN
               DEALLOCATE(LNIONS_TMP)
            ENDIF
            ALLOCATE(LNIONS_TMP(1:AB%MITYP,1:AB%MTYP,1:AB%NCONF))
            RMEM_HELP=RMEM_HELP+REAL(AB%MITYP,q)*REAL(AB%MTYP,q)*REAL(AB%NCONF,q)/1.0E+06_q
            IF(ALLOCATED(LNITYP_TMP)) THEN
               DEALLOCATE(LNITYP_TMP)
            ENDIF
            ALLOCATE(LNITYP_TMP(1:AB%MIONS,1:AB%NCONF))
            RMEM_HELP=RMEM_HELP+REAL(AB%MIONS,q)*REAL(AB%NCONF,q)/1.0E+06_q
            IF(ALLOCATED(LNTYP_TMP)) THEN
               DEALLOCATE(LNTYP_TMP)
            ENDIF
            ALLOCATE(LNTYP_TMP(1:AB%MIONS,1:AB%NCONF))
            RMEM_HELP=RMEM_HELP+REAL(AB%MIONS,q)*REAL(AB%NCONF,q)/1.0E+06_q
            IF(ALLOCATED(LSYS_TMP)) THEN
               DEALLOCATE(LSYS_TMP)
            ENDIF
            ALLOCATE(LSYS_TMP(1:AB%NCONF))
            RMEM_HELP=RMEM_HELP+REAL(AB%NCONF,q)/1.0E+06_q
            IF(ALLOCATED(NIONS_TMP)) THEN
               DEALLOCATE(NIONS_TMP)
            ENDIF
            ALLOCATE(NIONS_TMP(1:AB%NCONF))
            RMEM_HELP=RMEM_HELP+REAL(AB%NCONF,q)/1.0E+06_q
            IF(ALLOCATED(NITYP_TMP)) THEN
               DEALLOCATE(NITYP_TMP)
            ENDIF
            ALLOCATE(NITYP_TMP(1:AB%MTYP,1:AB%NCONF))
            RMEM_HELP=RMEM_HELP+REAL(AB%MTYP,q)*REAL(AB%NCONF,q)/1.0E+06_q
            IF(ALLOCATED(NTYP_TMP)) THEN
               DEALLOCATE(NTYP_TMP)
            ENDIF
            ALLOCATE(NTYP_TMP(1:AB%NCONF))
            RMEM_HELP=RMEM_HELP+REAL(AB%NCONF,q)/1.0E+06_q
            IF(ALLOCATED(OMEGA_TMP)) THEN
               DEALLOCATE(OMEGA_TMP)
            ENDIF
            ALLOCATE(OMEGA_TMP(1:AB%NCONF))
            RMEM_HELP=RMEM_HELP+REAL(AB%NCONF,q)*8/1.0E+06_q
            IF(ALLOCATED(POSION_TMP)) THEN
               DEALLOCATE(POSION_TMP)
            ENDIF
            ALLOCATE(POSION_TMP(1:3,1:AB%MIONS,1:AB%NCONF))
            RMEM_HELP=RMEM_HELP+REAL(AB%MIONS,q)*REAL(AB%NCONF,q)*3*8/1.0E+06_q
            IF(ALLOCAtED(SZNAM2_TMP)) THEN
               DEALLOCATE(SZNAM2_TMP)
            ENDIF
            ALLOCATE(SZNAM2_TMP(1:AB%NCONF))
            RMEM_HELP=RMEM_HELP+REAL(AB%NCONF,q)*8/1.0E+06_q
            IF(ALLOCATED(TOTEN_TMP)) THEN
               DEALLOCATE(TOTEN_TMP)
            ENDIF
            ALLOCATE(TOTEN_TMP(1:AB%NCONF))
            RMEM_HELP=RMEM_HELP+REAL(AB%NCONF,q)*8/1.0E+06_q
            IF(ALLOCATED(TIFOR_TMP)) THEN
               DEALLOCATE(TIFOR_TMP)
            ENDIF
            ALLOCATE(TIFOR_TMP(1:3,1:AB%MIONS,1:AB%NCONF))
            RMEM_HELP=RMEM_HELP+REAL(AB%MIONS,q)*REAL(AB%NCONF,q)*3*8/1.0E+06_q
            IF(ALLOCATED(TSIF_TMP)) THEN
               DEALLOCATE(TSIF_TMP)
            ENDIF
            ALLOCATE(TSIF_TMP(1:3,1:3,AB%NCONF))
            RMEM_HELP=RMEM_HELP+REAL(AB%NCONF,q)*9*8/1.0E+06_q
            MEM1%RMEM(5)=MAX(MEM1%RMEM(5),RMEM_HELP)
          END SUBROUTINE ALLOCATE_TMP_ARRAY2
#endif

! Deallocate temporary arrays

#ifndef ML_LEARNING_NOT_ACTIVE
          SUBROUTINE DEALLOCATE_TMP_ARRAY1
            IMPLICIT NONE
            IF(ALLOCATED(NBMAX_TMP)) THEN
               DEALLOCATE (NBMAX_TMP)
            ENDIF
          END SUBROUTINE DEALLOCATE_TMP_ARRAY1
#endif

#ifndef ML_LEARNING_NOT_ACTIVE
          SUBROUTINE DEALLOCATE_TMP_ARRAY2
            IMPLICIT NONE
            IF(ALLOCATED(A_TMP)) THEN
               DEALLOCATE(A_TMP)
            ENDIF
            IF(ALLOCATED(B_TMP)) THEN
               DEALLOCATE(B_TMP)
            ENDIF
            IF(ALLOCATED(ISYS_TMP)) THEN
               DEALLOCATE(ISYS_TMP)
            ENDIF
            IF(ALLOCATED(LADD_TMP)) THEN
               DEALLOCATE(LADD_TMP)
            ENDIF
            IF(ALLOCATED(LB_TMP)) THEN
               DEALLOCATE(LB_TMP)
            ENDIF
            IF(ALLOCATED(LNIONS_TMP)) THEN
               DEALLOCATE(LNIONS_TMP)
            ENDIF
            IF(ALLOCATED(LNITYP_TMP)) THEN
               DEALLOCATE(LNITYP_TMP)
            ENDIF
            IF(ALLOCATED(LNTYP_TMP)) THEN
               DEALLOCATE(LNTYP_TMP)
            ENDIF
            IF(ALLOCATED(LSYS_TMP)) THEN
               DEALLOCATE(LSYS_TMP)
            ENDIF
            IF(ALLOCATED(NIONS_TMP)) THEN
               DEALLOCATE(NIONS_TMP)
            ENDIF
            IF(ALLOCATED(NITYP_TMP)) THEN
               DEALLOCATE(NITYP_TMP)
            ENDIF
            IF(ALLOCATED(NTYP_TMP)) THEN
               DEALLOCATE(NTYP_TMP)
            ENDIF
            IF(ALLOCATED(OMEGA_TMP)) THEN
               DEALLOCATE(OMEGA_TMP)
            ENDIF
            IF(ALLOCATED(POSION_TMP)) THEN
               DEALLOCATE(POSION_TMP)
            ENDIF
            IF(ALLOCATED(SZNAM2_TMP)) THEN
               DEALLOCATE(SZNAM2_TMP)
            ENDIF
            IF(ALLOCATED(TOTEN_TMP)) THEN
               DEALLOCATE(TOTEN_TMP)
            ENDIF
            IF(ALLOCATED(TIFOR_TMP)) THEN
               DEALLOCATE(TIFOR_TMP)
            ENDIF
            IF(ALLOCATED(TSIF_TMP)) THEN
               DEALLOCATE(TSIF_TMP)
            ENDIF
          END SUBROUTINE DEALLOCATE_TMP_ARRAY2
#endif

        END SUBROUTINE  SET_AB_PAR


!****************************************************************************************************
! Configure new structural data
! (Not parallelized since it's fast)
!****************************************************************************************************

        SUBROUTINE SET_ABN_ST (TIM,ABN,A,B,OMEGA,POSION,POTIM,VEL)
        USE ML_FF_STRUCT, ONLY: TIME, ABINITIONEW
        IMPLICIT NONE
! Input variables
        TYPE (TIME)                   :: TIM
        TYPE (ABINITIONEW)            :: ABN
        REAL(q)          , INTENT(IN) :: A(:,:) !(1:3,1:3)
        REAL(q)          , INTENT(IN) :: B(:,:) !(1:3,1:3)
        REAL(q)          , INTENT(IN) :: OMEGA
        REAL(q)          , INTENT(IN) :: POSION(:,:) !(1:3,1:ABN%MIONS)
        REAL(q)          , INTENT(IN) :: POTIM
        REAL(q)          , INTENT(IN) :: VEL(:,:) !(1:3,1:ABN%MIONS)
! Local variables
        INTEGER                       :: INIONS,INTYP,INTYP0
        REAL(q)                       :: CPUTIM0,CPUTIM1
        REAL(q)                       :: VPUTIM0,VPUTIM1
! Start profiling
        PROFILING_START('set_abn_st')
! Measure the initial times
        CALL VTIME (VPUTIM0,CPUTIM0)
! Store new structural data
        ABN%POTIM=POTIM
        ABN%LSYS(ABN%NCONF)=1
        ABN%A(1:3,1:3,ABN%NCONF)=A(1:3,1:3)
        ABN%B(1:3,1:3,ABN%NCONF)=B(1:3,1:3)
        ABN%OMEGA(ABN%NCONF)=OMEGA
        DO INIONS=1, ABN%NIONS
           ABN%POSION(1:3,INIONS,ABN%NCONF)=POSION(1:3,INIONS)
              ABN%VEL(1:3,INIONS,ABN%NCONF)=VEL(1:3,INIONS)
        ENDDO
! Measure the final times and add them
        CALL VTIME (VPUTIM1,CPUTIM1)
        TIM%CPUTIM(1)=TIM%CPUTIM(1)+CPUTIM1-CPUTIM0
        TIM%VPUTIM(1)=TIM%VPUTIM(1)+VPUTIM1-VPUTIM0
! Stop profiling
        PROFILING_STOP('set_abn_st')
        END SUBROUTINE SET_ABN_ST

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Configure new energetical data
! (Partially parallelized by MPI)
!****************************************************************************************************

        SUBROUTINE SET_ABN_ENERGY (FF,TIM,COMM_WORLD,LOGF,ABN,IFF_FF,MTYP_FF,NB_FFM,NSTEP,&
                                   TIFOR,TOTEN,TSIF,TYPE_FF,DTIFOR_AV_MAE,DTIFOR_AV_RME,DTIFOR_MAX)
          USE ML_FF_STRUCT, ONLY: FF_PAR, TIME, ML_MPI_PAR, LOG_PAR, ABINITIONEW
          IMPLICIT NONE
! Input variables
          TYPE (FF_PAR)       :: FF
          TYPE (TIME)         :: TIM
          TYPE (ML_MPI_PAR)   :: COMM_WORLD
          TYPE (LOG_PAR)      :: LOGF
          TYPE (ABINITIONEW)  :: ABN
          LOGICAL, INTENT(IN) :: IFF_FF
          INTEGER, INTENT(IN) :: MTYP_FF
          INTEGER, INTENT(IN) :: NB_FFM(:) !(1:MTYP_FF)
          INTEGER, INTENT(IN) :: NSTEP
          REAL(q), INTENT(IN) :: TIFOR(:,:) !(1:3,1:ABN%MIONS)
          REAL(q), INTENT(IN) :: TOTEN
          REAL(q), INTENT(IN) :: TSIF(:,:) !(1:3,1:3)
          REAL(q), INTENT(INOUT) :: DTIFOR_AV_MAE
          REAL(q), INTENT(INOUT) :: DTIFOR_AV_RME
          REAL(q), INTENT(INOUT) :: DTIFOR_MAX
          CHARACTER(LEN=2), INTENT(IN) :: TYPE_FF(:) !(1:MTYP_FF)
! Local variables
          REAL(q)             :: CPUTIM0
          REAL(q)             :: CPUTIM1
          REAL(q)             :: DTIFOR(1:ABN%MIONS)
          REAL(q)             :: DTOTEN
          REAL(q)             :: DTSIF(1:3,1:3)
          REAL(q)             :: DTSIF_AV_MAE
          REAL(q)             :: DTSIF_AV_RME
          REAL(q)             :: DTSIF_MAX
!unit          REAL(q)             :: FAKT
          REAL(q)             :: VPUTIM0
          REAL(q)             :: VPUTIM1
! Start profiling
          PROFILING_START('set_abn_energy')
! Constant
!unit          FAKT=EVTOJ*1.0E+22_q/ABN%OMEGA(ABN%NCONF)/AUTOA**3
! Measure the initial times
          CALL VTIME (VPUTIM0,CPUTIM0)
          ABN%NSTEP(ABN%NCONF)=NSTEP
          ABN%TOTEN(ABN%NCONF)=TOTEN
          ABN%TSIF(1:3,1:3,ABN%NCONF)=TSIF(1:3,1:3)
          ABN%TIFOR(1:3,1:ABN%MIONS,ABN%NCONF)=TIFOR(1:3,1:ABN%MIONS)
! Set CTIFOR to store
          ABN%CTIFOR_ALLCONF(ABN%NCONF)=FF%CTIFOR
! Evaluate actuall errors and compare them with estimater errors
          IF(IFF_FF) THEN
             CALL CALCULATE_ERRORS
          ENDIF
! Measure the final times and add them
          CALL VTIME (VPUTIM1,CPUTIM1)
          TIM%CPUTIM(1)=TIM%CPUTIM(1)+CPUTIM1-CPUTIM0
          TIM%VPUTIM(1)=TIM%VPUTIM(1)+VPUTIM1-VPUTIM0
! Stop profiling
          PROFILING_STOP('set_abn_energy')

          CONTAINS

! Calculate actual errors

          SUBROUTINE CALCULATE_ERRORS
            IMPLICIT NONE
! Local variables
            INTEGER             :: IERR
            INTEGER             :: INIONS
            INTEGER             :: INIONS_TMP
            INTEGER             :: IXYZ
            INTEGER             :: JXYZ
! Calculate real error in energy.
            DTOTEN=ABS(ABN%TOTEN(ABN%NCONF)-ABN%TOTEN_MB(ABN%NCONF))/REAL(ABN%NIONS,q)
! Calculate real errors in forces.
            DTIFOR=0.0_q
            DTIFOR_AV_MAE=0.0_q
            DTIFOR_AV_RME=0.0_q
            DTIFOR_MAX=0.0_q
            DO INIONS_TMP=1, ABN%NIONS, COMM_WORLD%NCPU
               INIONS=INIONS_TMP+COMM_WORLD%NODE_ME
               IF(INIONS.LE.ABN%NIONS) THEN
                  DO IXYZ=1, 3
                     DTIFOR_AV_MAE =DTIFOR_AV_MAE +ABS(ABN%TIFOR(IXYZ,INIONS,ABN%NCONF)-ABN%TIFOR_MB(IXYZ,INIONS,ABN%NCONF))
                     DTIFOR_AV_RME =DTIFOR_AV_RME +(ABN%TIFOR(IXYZ,INIONS,ABN%NCONF)-ABN%TIFOR_MB(IXYZ,INIONS,ABN%NCONF))**2
                     DTIFOR(INIONS)=DTIFOR(INIONS)+(ABN%TIFOR(IXYZ,INIONS,ABN%NCONF)-ABN%TIFOR_MB(IXYZ,INIONS,ABN%NCONF))**2
                  ENDDO
                  DTIFOR(INIONS)=SQRT(DTIFOR(INIONS))
                  DTIFOR_MAX=MAX(DTIFOR_MAX,DTIFOR(INIONS))
               ENDIF
            ENDDO
            CALL MPI_ALLREDUCE(MPI_IN_PLACE,DTIFOR(1) ,ABN%MIONS,MPI_REAL8,MPI_SUM,COMM_WORLD%MPI_COMM,IERR)
            CALL MPI_ALLREDUCE(MPI_IN_PLACE,DTIFOR_AV_MAE ,1,MPI_REAL8,MPI_SUM,COMM_WORLD%MPI_COMM,IERR)
            CALL MPI_ALLREDUCE(MPI_IN_PLACE,DTIFOR_AV_RME ,1,MPI_REAL8,MPI_SUM,COMM_WORLD%MPI_COMM,IERR)
            CALL MPI_ALLREDUCE(MPI_IN_PLACE,DTIFOR_MAX,1,MPI_REAL8,MPI_MAX,COMM_WORLD%MPI_COMM,IERR)
            DTIFOR_AV_MAE=DTIFOR_AV_MAE/REAL(3*ABN%NIONS,q)
            DTIFOR_AV_RME=SQRT(DTIFOR_AV_RME/REAL(3*ABN%NIONS,q))
! Calculate real errors in stress-tensors.
            DTSIF=0.0_q
            DTSIF_AV_MAE=0_q
            DTSIF_AV_RME=0_q
            DTSIF_MAX=0_q
            DO IXYZ=1, 3
               DO JXYZ=IXYZ, 3
                  DTSIF(JXYZ,IXYZ)=ABS(ABN%TSIF(JXYZ,IXYZ,ABN%NCONF)-ABN%TSIF_MB(JXYZ,IXYZ,ABN%NCONF))
                  DTSIF_AV_MAE=DTSIF_AV_MAE+ABS(ABN%TSIF(JXYZ,IXYZ,ABN%NCONF)-ABN%TSIF_MB(JXYZ,IXYZ,ABN%NCONF))
                  DTSIF_AV_RME=DTSIF_AV_RME+(ABN%TSIF(JXYZ,IXYZ,ABN%NCONF)-ABN%TSIF_MB(JXYZ,IXYZ,ABN%NCONF))**2
                  DTSIF_MAX=MAX(DTSIF_MAX,ABS(ABN%TSIF(JXYZ,IXYZ,ABN%NCONF)-ABN%TSIF_MB(JXYZ,IXYZ,ABN%NCONF)))
               ENDDO
            ENDDO
            DTSIF_AV_MAE=DTSIF_AV_MAE/6.0_q
            DTSIF_AV_RME=SQRT(DTSIF_AV_RME/6.0_q)
          END SUBROUTINE CALCULATE_ERRORS

        END SUBROUTINE SET_ABN_ENERGY
#endif

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Add new ab initio data
! (Not parallelized since all data must be shared)
!****************************************************************************************************

        SUBROUTINE PUT_ABN_DATA (FF,PAR_SUP_HANDLE,TIM,AB,ABN,CTIFOR_FF,IERR,IFF_FF,IFF_OLD_FF,ISTART_FF, &
                   LABORT,LDOUB_FF,LFMAT_FF,LGENFF_FF,LLIST_FF,LOVER_FF,LPS_FF,LSAMPLE_FF,LSIG_FF,LSOAP_FF,NSTEP,SIG_FFM,FORCING_UPDATE)
          USE ML_FF_STRUCT, ONLY: FF_PAR,PARALLEL_SUPER,TIME, ABINITIO, ABINITIONEW
#ifdef ML_DEBUG
          USE PRINT2TERMINAL
#endif
          IMPLICIT NONE
! Input variables
          TYPE (FF_PAR)          :: FF
          TYPE (PARALLEL_SUPER)  :: PAR_SUP_HANDLE
          TYPE (TIME)            :: TIM
          TYPE (ABINITIO)        :: AB
          TYPE (ABINITIONEW)     :: ABN
          LOGICAL, INTENT(IN)    :: FORCING_UPDATE ! Forces to set all control variables to .TRUE. and STORE_ABN
          LOGICAL, INTENT(IN)    :: IFF_FF ! IFF specifies the presence or absence of FF.
          INTEGER, INTENT(IN)    :: ISTART_FF
          LOGICAL, INTENT(IN)    :: LABORT ! Logical variable provided from VASP. When
                         ! LABORT = .TRUE., some problems occur in
                         ! electronic structure calculation. Problems
                         ! can be no-convergence of electronic
                         ! structure, for example. In this case, data is
                         ! not stored.
          INTEGER, INTENT(IN)    :: NSTEP
          REAL(q), INTENT(IN)    :: SIG_FFM(:) !(1:10)
! Input or output variables.
          REAL(q), INTENT(INOUT)    :: CTIFOR_FF ! Criteria for force
          LOGICAL, INTENT(INOUT) :: IERR ! If IERR=0, force field is predicting energetics accurately.
          LOGICAL, INTENT(INOUT) :: IFF_OLD_FF ! Old IFF data.
          LOGICAL, INTENT(INOUT) :: LDOUB_FF ! If LDOUB=.TRUE., the sampling is stopped, and a new force field is generated.
          LOGICAL, INTENT(INOUT) :: LFMAT_FF ! If LFMAT_FF=.TRUE., the design metrix elements are calculated.
          LOGICAL, INTENT(INOUT) :: LGENFF_FF ! If LGENFF_FF=.TRUE., FF generation will be executed.
          LOGICAL, INTENT(INOUT) :: LLIST_FF ! If LLIST_FF=.TRUE., the machine generates lists storing the location of the design-matrix elements.
          LOGICAL, INTENT(INOUT) :: LOVER_FF ! If LOVER = .TRUE., sparsification will be executed
          LOGICAL, INTENT(INOUT) :: LPS_FF  ! If LPS = .TRUE., PS calculations will be executed
          LOGICAL, INTENT(INOUT) :: LSAMPLE_FF ! If LSAMPLE=.TRUE., sampling is executed.
          LOGICAL, INTENT(INOUT) :: LSIG_FF ! if LSIG_FF=.TRUE., the machine will calculate the variance of the training data.
          LOGICAL, INTENT(INOUT) :: LSOAP_FF ! If LSOAP = .TRUE., soap matrix elements are calculated.
! Local variables.
          REAL(q)                :: CPUTIM0
          REAL(q)                :: CPUTIM1
          INTEGER                :: INIONS
          INTEGER                :: ISTAT
          REAL(q)                :: VPUTIM0
          REAL(q)                :: VPUTIM1
#ifdef ML_DEBUG
! variables for debugging
          INTEGER                :: CHECKSUM
          INTEGER                :: CHECKSUM_HELP
#endif
! Start profiling.
          PROFILING_START('put_abn_data')
! Measure the initial times.
          CALL VTIME (VPUTIM0,CPUTIM0)
! Store the old IFF data.
          IFF_OLD_FF=IFF_FF
! Calculate Bayesian errors
          IERR=.FALSE.
          DO INIONS=1, ABN%NIONS
             IF(SQRT(ABN%BE(3*(INIONS-1)+2,ABN%NCONF)+ &
                      ABN%BE(3*(INIONS-1)+3,ABN%NCONF)+ &
                      ABN%BE(3*(INIONS-1)+4,ABN%NCONF)) .GT. CTIFOR_FF) THEN
                IERR=.TRUE.
                EXIT
             ENDIF
          ENDDO
#ifdef ML_DEBUG
          CHECKSUM_HELP=0
          IF (IERR) THEN
             CHECKSUM=1
          ELSE
             CHECKSUM=0
          ENDIF
          CALL MPI_ALLREDUCE(CHECKSUM_HELP,CHECKSUM,1,MPI_INTEGER,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,ISTAT)
          IF ((CHECKSUM*PAR_SUP_HANDLE%COMM_WORLD%NCPU).NE.CHECKSUM_HELP) THEN
             CALL ERROR_PRINT2TERMINAL("IERR is not the &
                  same on all cores in PUT_ABN_DATA")
          ENDIF
#endif
           CALL M_BCAST (PAR_SUP_HANDLE%COMM_WORLD,IERR,0)
! Add data?
! If LABORT=.TRUE., some problems occur in electronic structure
! calculations in VASP. In this case, the problematic error should not
! be added to the database. Thus, ABN%NCONF is not renewed, and the new
! ABN%NCONF-th ab initio data will be recalculated.
          IF(LABORT) THEN
             LFMAT_FF =.FALSE.
             LGENFF_FF=.FALSE.
             LLIST_FF =.FALSE.
             LOVER_FF =.FALSE.
             LPS_FF   =.FALSE.
             LSIG_FF  =.FALSE.
             LSOAP_FF =.FALSE.
! Otherwise, the ab initio data is correctly calculated. Thus, check its
! availability.
          ELSE
! Execute a following sampling method.----------------------------------------------
! 1) If the estimated error is larger than the product of CDOUB by the predetermined threshold
!    LDO_AB_INITIO is set as .TRUE. to execute the ab initio calculation, and LDOUB is set as 
!    .TRUE. to generate the force field immediately.
! 2) Otherwise, the following procedure is executed.
!         2-1) If the estimated error is larger than the predetermined threshold, and if sampling
!              is switched on (LSAMPLE=.TRUE.), LDO_AB_INITIO is set as .TRUE. to execute the 
!              ab initio calculation.
!         2-2) Otherwise, LDO_AB_INITIO is set as .FALSE. to skip the ab initio calculation.
!         2-3) Regarless 2-1) and 2-2), if the LSAMPLE is switched off, and if the number of MD step 
!              from the previous force field generation exceeds NDSTEP, 
!              LSAMPLE is set as .TURE. to switch on the sampling.
            IF(((IERR) .AND. (ABN%NCONF .EQ. ABN%MCONF) .AND. (LSAMPLE_FF)) .OR. &
                ((IERR) .AND. (FF%NSTEP .EQ. FF%NSW) .AND. (LSAMPLE_FF)) .OR. & 
                (LDOUB_FF) .OR.  &
               (.NOT. IFF_FF)) THEN
                CALL STORE_ABN
                LFMAT_FF =.TRUE.
                LGENFF_FF=.TRUE.
                LLIST_FF =.TRUE.
                LOVER_FF =.TRUE.
                LPS_FF   =.TRUE.
                LSIG_FF  =.TRUE.
                LSOAP_FF =.TRUE.
            ELSE
! If FORCING_UPDATE is true then all parameters are turned on
! This has usually the following two reasons:
! The last step in the ionic cycle is reached
! The last system of the same type in the ML_AB file is reached when ML_FF_ISTART=3
               IF (FORCING_UPDATE) THEN
                  CALL STORE_ABN
                  LFMAT_FF =.TRUE.
                  LGENFF_FF=.TRUE.
                  LLIST_FF =.TRUE.
                  LOVER_FF =.TRUE.
                  LPS_FF   =.TRUE.
                  LSIG_FF  =.TRUE.
                  LSOAP_FF =.TRUE.
               ELSE
                  LFMAT_FF =.FALSE.
                  LGENFF_FF=.FALSE.
                  LLIST_FF =.FALSE.
                  LOVER_FF =.FALSE.
                  LPS_FF   =.FALSE.
                  LSIG_FF  =.FALSE.
                  LSOAP_FF =.FALSE.
               ENDIF
            ENDIF
          ENDIF
! Measure the final times and add them
          CALL VTIME (VPUTIM1,CPUTIM1)
          TIM%CPUTIM(1)=TIM%CPUTIM(1)+CPUTIM1-CPUTIM0
          TIM%VPUTIM(1)=TIM%VPUTIM(1)+VPUTIM1-VPUTIM0
! Stop profiling
          PROFILING_STOP('put_abn_data')

          CONTAINS

          SUBROUTINE STORE_ABN
            IMPLICIT NONE
! Local variables
            INTEGER                :: ICONF
            INTEGER                :: NCONF_LOC
            NCONF_LOC=AB%NCONF
            DO ICONF=1, ABN%NCONF
               NCONF_LOC=NCONF_LOC+1
               AB%ISYS(NCONF_LOC)=ABN%ISYS
               AB%LADD(1:ABN%MTYP,NCONF_LOC)=ABN%LADD(1:ABN%MTYP)
               AB%LNITYP(1:ABN%MIONS,NCONF_LOC)=ABN%LNITYP(1:ABN%MIONS)
               AB%LNIONS(1:ABN%MITYP,1:ABN%MTYP,NCONF_LOC)=ABN%LNIONS(1:ABN%MITYP,1:ABN%MTYP)
               AB%LNTYP(1:ABN%MIONS,NCONF_LOC)=ABN%LNTYP(1:ABN%MIONS)
               AB%LSYS(NCONF_LOC)=ABN%LSYS(ICONF)
               AB%NITYP(1:ABN%MTYP,NCONF_LOC)=ABN%NITYP(1:ABN%MTYP)
               AB%NTYP(NCONF_LOC)=ABN%NTYP
               AB%NIONS(NCONF_LOC)=ABN%NIONS
               AB%A(1:3,1:3,NCONF_LOC)=ABN%A(1:3,1:3,ICONF)
               AB%B(1:3,1:3,NCONF_LOC)=ABN%B(1:3,1:3,ICONF)
               AB%OMEGA(NCONF_LOC)=ABN%OMEGA(ICONF)
               AB%POSION(1:3,1:ABN%MIONS,NCONF_LOC)=ABN%POSION(1:3,1:ABN%MIONS,ICONF)
               AB%SZNAM2(NCONF_LOC)=ABN%SZNAM2
               AB%TOTEN(NCONF_LOC)=ABN%TOTEN(ICONF)
               AB%TIFOR(1:3,1:ABN%MIONS,NCONF_LOC)=ABN%TIFOR(1:3,1:ABN%MIONS,ICONF)
               AB%TSIF(1:3,1:3,NCONF_LOC)=ABN%TSIF(1:3,1:3,ICONF)
               IF (FF%ISTART.LT.3) THEN
                  AB%CTIFOR_ALLCONF_NEW(NCONF_LOC)=ABN%CTIFOR_ALLCONF(ICONF)
               ENDIF
            ENDDO
          END SUBROUTINE STORE_ABN

        END SUBROUTINE PUT_ABN_DATA
#endif

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Calculate the histgram of data distribution
!****************************************************************************************************

        SUBROUTINE HISTOGRAM (PAR_SUP_HANDLE,AB,DIR_APP,DIR_LEN)
          USE ML_FF_STRUCT, ONLY: PARALLEL_SUPER, ABINITIO
          IMPLICIT NONE
          TYPE (PARALLEL_SUPER)         :: PAR_SUP_HANDLE
          TYPE (ABINITIO)               :: AB
! Input variables
          CHARACTER(LEN=10), INTENT(IN) :: DIR_APP     ! Location of directory passed from VASP.
          INTEGER, INTENT(IN)           :: DIR_LEN     ! Length of DIR_APP
! Local variables
          REAL(q) AV(1:10)
          REAL(q) DX
          INTEGER I
          INTEGER ICONF
          INTEGER INIONS
          INTEGER IX
          INTEGER IXYZ
          INTEGER JXYZ
          INTEGER KXYZ
          INTEGER NDATA(1:10)
          INTEGER, PARAMETER :: NX=20
          REAL(q) SIG(1:10)
          REAL(q) X(1:(NX+1))
          REAL(q) XMIN,XMAX
          REAL(q) Y(1:(NX+1),1:10)
! Calculate the number of data and their averages
          AV=0.0_q
          NDATA=0
          DO ICONF=1, AB%NCONF
             NDATA(1)=NDATA(1)+1
             AV(1)=AV(1)+AB%TOTEN(ICONF)/REAL(AB%NIONS(ICONF),q)
             DO INIONS=1, AB%NIONS(ICONF)
                DO IXYZ=1, 3
                   NDATA(IXYZ+1)=NDATA(IXYZ+1)+1
                   AV(IXYZ+1)=AV(IXYZ+1)+AB%TIFOR(IXYZ,INIONS,ICONF)/REAL(AB%NIONS(ICONF),q)
                ENDDO
             ENDDO
             KXYZ=0
             DO IXYZ=1, 3
                DO JXYZ=IXYZ, 3
                   KXYZ=KXYZ+1
                   NDATA(KXYZ+4)=NDATA(KXYZ+4)+1
                   AV(KXYZ+4)=AV(KXYZ+4)+AB%TSIF(JXYZ,IXYZ,ICONF)/REAL(AB%NIONS(ICONF),q)
                ENDDO
             ENDDO
          ENDDO
          DO I=1, 10
             AV(I)=AV(I)/REAL(NDATA(I),q)
          ENDDO
! Calculate variance of the training data
          SIG=0.0_q
          DO ICONF=1, AB%NCONF
             SIG(1)=SIG(1)+(AB%TOTEN(ICONF)/REAL(AB%NIONS(ICONF),q)-AV(1))**2
             DO INIONS=1, AB%NIONS(ICONF)
                DO IXYZ=1, 3
                   SIG(IXYZ+1)=SIG(IXYZ+1)+(AB%TIFOR(IXYZ,INIONS,ICONF)/REAL(AB%NIONS(ICONF),q)-AV(IXYZ+1))**2
                ENDDO
             ENDDO
             KXYZ=0
             DO IXYZ=1, 3
                DO JXYZ=IXYZ, 3
                   KXYZ=KXYZ+1
                   SIG(KXYZ+4)=SIG(KXYZ+4)+(AB%TSIF(JXYZ,IXYZ,ICONF)/REAL(AB%NIONS(ICONF),q)-AV(KXYZ+4))**2
                ENDDO
             ENDDO
          ENDDO
          DO I=1, 10
             SIG(I)=SQRT(SIG(I)/REAL(NDATA(I),q))
          ENDDO
! Determine the boundaries
          XMIN=-3.0_q
          XMAX= 3.0_q
! Determine the increment
          DX=(XMAX-XMIN)/REAL(NX,q)
! Determine the X-mesh
          DO IX=1, NX+1
             X(IX)=XMIN+(IX-1)*DX
          ENDDO
! Calculate the histogram
          Y=0.0_q
          DO ICONF=1, AB%NCONF
             IF(SIG(1) .GT. 0.0_q) THEN
                IX=INT(((AB%TOTEN(ICONF)/REAL(AB%NIONS(ICONF),q)-AV(1))/SIG(1)-XMIN)/DX)+1
                IF((IX .GE. 1) .AND. (IX .LE. NX+1)) THEN
                   Y(IX,1)=Y(IX,1)+1.0_q
                ENDIF
             ENDIF
             DO INIONS=1, AB%NIONS(ICONF)
                DO IXYZ=1, 3
                   IF(SIG(IXYZ+1) .GT. 0.0_q) THEN
                     IX=INT(((AB%TIFOR(IXYZ,INIONS,ICONF)/REAL(AB%NIONS(ICONF),q)-AV(IXYZ+1))/SIG(IXYZ+1)-XMIN)/DX)+1
                     IF((IX .GE. 1) .AND. (IX .LE. NX+1)) THEN
                        Y(IX,IXYZ+1)=Y(IX,IXYZ+1)+1.0_q
                     ENDIF
                   ENDIF
                ENDDO
             ENDDO
             KXYZ=0
             DO IXYZ=1, 3
                DO JXYZ=IXYZ, 3
                   KXYZ=KXYZ+1
                   IF(SIG(KXYZ+4) .GT. 0.0_q) THEN
                      IX=INT(((AB%TSIF(JXYZ,IXYZ,ICONF)/REAL(AB%NIONS(ICONF),q)-AV(KXYZ+4))/SIG(KXYZ+4)-XMIN)/DX)+1
                      IF((IX .GE. 1) .AND. (IX .LE. NX+1)) THEN
                         Y(IX,KXYZ+4)=Y(IX,KXYZ+4)+1.0_q
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
          DO I=1, 10
             DO IX=1, NX
                Y(IX,I)=Y(IX,I)/REAL(NDATA(I),q)
             ENDDO
          ENDDO
! Output histogra
          CALL OUT_HISTOGRAM (PAR_SUP_HANDLE%COMM_WORLD,AV,DIR_APP,DIR_LEN,AB%HISCAR,AB%IU4,NX,SIG,X,Y)
        END SUBROUTINE HISTOGRAM
#endif

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! In test calculations, where ab initio calculation is always executed
! to compare the ab initio data and machine-learning data,
! machine-learning data need to be restored when it is judged to be
! reliable in order to obtain machine-learning MD trajectory. For this
! purpose, this subroutine restore the machine-learning data.
! This subroutine is used in a following situation.
! 1) LTEST=.TRUE. (Test calculation is turned on)
! 2) LTRJ=.FALSE. (Machine-learning trajectory is desired)
! 3) IERR=.FALSE. (Force field is accurate enough)
!****************************************************************************************************

        SUBROUTINE RESTORE_ML_DATA (COMM_WORLD,LOGF,ABN,NSTEP)
          USE ML_FF_STRUCT, ONLY: ML_MPI_PAR,LOG_PAR, ABINITIONEW
          IMPLICIT NONE
! Input variables
          TYPE (ML_MPI_PAR)      :: COMM_WORLD
          TYPE (LOG_PAR)      :: LOGF
          TYPE (ABINITIONEW)  :: ABN
          INTEGER, INTENT(IN) :: NSTEP
! Start profiling
          PROFILING_START('restore_ml_data')
! Restore machine-learning energy, force and stress
          ABN%TOTEN(ABN%NCONF)=ABN%TOTEN_MB(ABN%NCONF)
          ABN%TIFOR(1:3,1:ABN%MIONS,ABN%NCONF)=ABN%TIFOR_MB(1:3,1:ABN%MIONS,ABN%NCONF)
          ABN%TSIF(1:3,1:3,ABN%NCONF)=ABN%TSIF_MB(1:3,1:3,ABN%NCONF)
! Stop profiling
          PROFILING_STOP('restore_ml_data')
        END SUBROUTINE RESTORE_ML_DATA
#endif

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Reset ABN%NCONF
!****************************************************************************************************

        SUBROUTINE RENEW_SAMPLING_PAR (FF,TIM,ABN,IERR,LABORT,NSTEP)
          USE ML_FF_STRUCT, ONLY: FF_PAR, TIME, ABINITIONEW, ML_MPI_PAR
          IMPLICIT NONE
! Input variables
          TYPE (FF_PAR)          :: FF
          TYPE (TIME)            :: TIM
          TYPE (ABINITIONEW)     :: ABN
          LOGICAL, INTENT(IN)    :: IERR
          LOGICAL, INTENT(IN)    :: LABORT
          INTEGER, INTENT(IN)    :: NSTEP
! Local variables
          REAL(q) :: CERR
          REAL(q) :: ERR
          REAL(q) :: VPUTIM0,CPUTIM0
          REAL(q) :: VPUTIM1,CPUTIM1
! Start profiling
          PROFILING_START('renew_sampling_par')
! Measure the initial times
          CALL VTIME (VPUTIM0,CPUTIM0)
! Add data?
! If LABORT=.TRUE., some problems occur in electronic structure
! calculations in VASP. In this case, the problematic error was not
! added to the database. Thus, ABN%NCONF is not renewed.
! Otherwise, execute following procedure:
! 1) If the estimated error is larger than the product of CDOUB by the predetermined threshold
!    LDO_AB_INITIO is set as .TRUE. to execute the ab initio calculation, and LDOUB is set as 
!    .TRUE. to generate the force field immediately.
! 2) Otherwise, the following procedure is executed.
!         2-1) If the estimated error is larger than the predetermined threshold, and if sampling
!              is switched on (LSAMPLE=.TRUE.), LDO_AB_INITIO is set as .TRUE. to execute the 
!              ab initio calculation.
!         2-2) Otherwise, LDO_AB_INITIO is set as .FALSE. to skip the ab initio calculation.
!         2-3) Regarless 2-1) and 2-2), if the LSAMPLE is switched off, and if the number of MD step 
!              from the previous force field generation exceeds NDSTEP, 
!              LSAMPLE is set as .TURE. to switch on the sampling.
          IF(.NOT.LABORT) THEN
             IF(((IERR).AND.(ABN%NCONF.EQ.ABN%MCONF).AND.(FF%LSAMPLE)) .OR. &
                (.NOT.FF%IFF_OLD).OR. &
                (FF%LDOUB)) THEN
                ABN%NCONF=1
                FF%LDOUB=.FALSE.
                FF%LSAMPLE=.FALSE.
             ELSE IF ((IERR).AND.(ABN%NCONF.LT.ABN%MCONF).AND.(FF%LSAMPLE)) THEN
                IF (FF%NSTEP.EQ.FF%NSW) THEN
                   ABN%NCONF=1
                   FF%LDOUB=.FALSE.
                   FF%LSAMPLE=.FALSE.
                ElSE
                   ABN%NCONF=ABN%NCONF+1
                   FF%LSAMPLE=.FALSE.
                ENDIF
             ELSE 
! When we read from database (FF%ISTART=3) then all structures should be kept
! so FF%LSAMPLE is always set as .TRUE..
                IF (FF%ISTART.EQ.3) THEN
                   IF (FF%LDISCARD_STRUCTURES_NOT_GIVING_BASIS) THEN
                      IF ((IERR).AND.(ABN%NCONF.LT.ABN%MCONF)) THEN
                         ABN%NCONF=ABN%NCONF+1
                      ENDIF
                   ELSE
                      IF ((ABN%NCONF.EQ.ABN%MCONF)) THEN
                         ABN%NCONF=1
                      ELSE
                         ABN%NCONF=ABN%NCONF+1
                      ENDIF
                   ENDIF
                   FF%LSAMPLE=.TRUE.
                ENDIF
              ENDIF
          ENDIF
          IF (FF%FORCING_UPDATE) THEN
             ABN%NCONF=1
             FF%FORCING_UPDATE = .FALSE.
          ENDIF
! Measure the final times and add them
          CALL VTIME (VPUTIM1,CPUTIM1)
          TIM%CPUTIM(1)=TIM%CPUTIM(1)+CPUTIM1-CPUTIM0
          TIM%VPUTIM(1)=TIM%VPUTIM(1)+VPUTIM1-VPUTIM0
! Stop profiling
          PROFILING_STOP('renew_sampling_par')
        END SUBROUTINE RENEW_SAMPLING_PAR
#endif

!****************************************************************************************************
! Allocate arrays for new ab initio data
! TIFOR, TOTEN, TSIF for each component of interactiion are allocated
! even when that component is turned off for simplifying the program
! structure.
!****************************************************************************************************

        SUBROUTINE ALLOCATE_ABN_ARRAY (MEM1,ABN)
          USE ML_FF_STRUCT, ONLY: MEMORY, ABINITIONEW
          IMPLICIT NONE
! Input variables
          TYPE (MEMORY)       :: MEM1
          TYPE (ABINITIONEW)  :: ABN
! Local variable
          REAL(q) RMEM_LOC
          RMEM_LOC=0.0_q
          IF(ALLOCATED(ABN%A)) THEN
             DEALLOCATE(ABN%A)
          ENDIF
          ALLOCATE(ABN%A(1:3,1:3,1:ABN%MCONF))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MCONF,q)*9*8/1.0E+06_q
          IF(ALLOCATED(ABN%B)) THEN
             DEALLOCATE(ABN%B)
          ENDIF
          ALLOCATE(ABN%B(1:3,1:3,1:ABN%MCONF))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MCONF,q)*9*8/1.0E+06_q
          IF(ALLOCATED(ABN%BE)) THEN
             DEALLOCATE(ABN%BE)
          ENDIF
          ALLOCATE(ABN%BE(1:(3*ABN%MIONS_ALLOCATE+7),1:ABN%MCONF))
          RMEM_LOC=RMEM_LOC+REAL((3*ABN%MIONS_ALLOCATE+7),q)*REAL(ABN%MCONF,q)/1.0E+06_q
          IF(ALLOCATED(ABN%BEAV)) THEN
             DEALLOCATE(ABN%BEAV)
          ENDIF
          ALLOCATE(ABN%BEAV(1:3,1:ABN%MCONF))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MCONF,q)/1.0E+06_q
          IF(ALLOCATED(ABN%BEMAX)) THEN
             DEALLOCATE(ABN%BEMAX)
          ENDIF
          ALLOCATE(ABN%BEMAX(1:3,1:ABN%MCONF))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MCONF,q)/1.0E+06_q
          IF(ALLOCATED(ABN%CTIFOR_ALLCONF)) THEN
             DEALLOCATE(ABN%CTIFOR_ALLCONF)
          ENDIF
          ALLOCATE(ABN%CTIFOR_ALLCONF(1:ABN%MCONF))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MCONF,q)/1.0E+06_q
          IF(ALLOCATED(ABN%EATOM_KIN)) THEN
             DEALLOCATE(ABN%EATOM_KIN)
          ENDIF
          ALLOCATE(ABN%EATOM_KIN(1:ABN%MIONS_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MIONS_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(ABN%EATOM_POT)) THEN
             DEALLOCATE(ABN%EATOM_POT)
          ENDIF
          ALLOCATE(ABN%EATOM_POT(1:ABN%MIONS_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MIONS_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(ABN%LADD)) THEN
             DEALLOCATE(ABN%LADD)
          ENDIF
          ALLOCATE(ABN%LADD(1:ABN%MTYP_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MTYP_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(ABN%LNITYP)) THEN
             DEALLOCATE(ABN%LNITYP)
          ENDIF
          ALLOCATE(ABN%LNITYP(1:ABN%MIONS_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MIONS_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(ABN%LNIONS)) THEN
             DEALLOCATE(ABN%LNIONS)
          ENDIF
          ALLOCATE(ABN%LNIONS(1:ABN%MITYP_ALLOCATE,1:ABN%MTYP_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MITYP_ALLOCATE,q)*REAL(ABN%MTYP_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(ABN%LNTYP)) THEN
             DEALLOCATE(ABN%LNTYP)
          ENDIF
          ALLOCATE(ABN%LNTYP(1:ABN%MIONS_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MIONS_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(ABN%LSYS)) THEN
             DEALLOCATE(ABN%LSYS)
          ENDIF
          ALLOCATE(ABN%LSYS(1:ABN%MCONF))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MCONF,q)/1.0E+06_q
          IF(ALLOCATED(ABN%NITYP)) THEN
             DEALLOCATE(ABN%NITYP)
          ENDIF
          ALLOCATE(ABN%NITYP(1:ABN%MTYP_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MTYP_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(ABN%NSTEP)) THEN
             DEALLOCATE(ABN%NSTEP)
          ENDIF
          ALLOCATE(ABN%NSTEP(1:ABN%MCONF))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MCONF,q)/1.0E+06_q
          IF(ALLOCATED(ABN%OMEGA)) THEN
             DEALLOCATE(ABN%OMEGA)
          ENDIF
          ALLOCATE(ABN%OMEGA(1:ABN%MCONF))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MCONF,q)*8/1.0E+06_q
          IF(ALLOCATED(ABN%POMASS)) THEN
             DEALLOCATE(ABN%POMASS)
          ENDIF
          ALLOCATE(ABN%POMASS(1:ABN%MTYP_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MTYP_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(ABN%POSION)) THEN
             DEALLOCATE(ABN%POSION)
          ENDIF
          ALLOCATE(ABN%POSION(1:3,1:ABN%MIONS_ALLOCATE,1:ABN%MCONF))
          RMEM_LOC=RMEM_LOC+REAl(ABN%MIONS_ALLOCATE,q)*3*REAL(ABN%MCONF,q)*8/1.0E+06_q
          IF(ALLOCATED(ABN%QHEAT)) THEN
             DEALLOCATE(ABN%QHEAT)
          ENDIF
          ALLOCATE(ABN%QHEAT(1:3,1:ABN%MCONF))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MCONF,q)*3*8/1.0E+06_q
          IF(ALLOCATED(ABN%SZNAM2_ADD)) THEN
             DEALLOCATE(ABN%SZNAM2_ADD)
          ENDIF
          ALLOCATE(ABN%SZNAM2_ADD(1:ABN%NSYS_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(ABN%NSYS_ALLOCATE,q)*40/1.0E+06_q
          IF(ALLOCATED(ABN%TIFOR)) THEN
             DEALLOCATE(ABN%TIFOR)
          ENDIF
          ALLOCATE(ABN%TIFOR(1:3,1:ABN%MIONS_ALLOCATE,1:ABN%MCONF))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MIONS_ALLOCATE,q)*3*REAL(ABN%MCONF,q)*8/1.0E+06_q
          IF(ALLOCATED(ABN%TIFOR_MB)) THEN
             DEALLOCATE(ABN%TIFOR_MB)
          ENDIF
          ALLOCATE(ABN%TIFOR_MB(1:3,1:ABN%MIONS_ALLOCATE,1:ABN%MCONF))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MIONS_ALLOCATE,q)*3*REAL(ABN%MCONF,q)*8/1.0E+06_q
          IF(ALLOCATED(ABN%TOTEN)) THEN
             DEALLOCATE(ABN%TOTEN)
          ENDIF
          ALLOCATE(ABN%TOTEN(1:ABN%MCONF))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MCONF,q)*8/1.0E+06_q
          IF(ALLOCATED(ABN%TOTEN_MB)) THEN
             DEALLOCATE(ABN%TOTEN_MB)
          ENDIF
          ALLOCATE(ABN%TOTEN_MB(1:ABN%MCONF))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MCONF,q)*8/1.0E+06_q
          IF(ALLOCATED(ABN%TSIF)) THEN
             DEALLOCATE(ABN%TSIF)
          ENDIF
          ALLOCATE(ABN%TSIF(1:3,1:3,1:ABN%MCONF))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MCONF,q)*9*8/1.0E+06_q
          IF(ALLOCATED(ABN%TSIF_MB)) THEN
             DEALLOCATE(ABN%TSIF_MB)
          ENDIF
          ALLOCATE(ABN%TSIF_MB(1:3,1:3,1:ABN%MCONF))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MCONF,q)*9*8/1.0E+06_q
          IF(ALLOCATED(ABN%TYPE)) THEN
             DEALLOCATE(ABN%TYPE)
          ENDIF
          ALLOCATE(ABN%TYPE(1:ABN%MTYP_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MCONF,q)*2/1.0E+06_q
          IF(ALLOCATED(ABN%VEL)) THEN
             DEALLOCATE(ABN%VEL)
          ENDIF
          ALLOCATE(ABN%VEL(1:3,1:ABN%MIONS_ALLOCATE,1:ABN%MCONF))
          RMEM_LOC=RMEM_LOC+REAL(ABN%MIONS_ALLOCATE,q)*3*REAL(ABN%MCONF,q)*8/1.0E+06_q
          MEM1%RMEM(3)=MAX(MEM1%RMEM(3),RMEM_LOC)
! Initialization of variables
          ABN%A=0.0_q
          ABN%B=0.0_q
          ABN%CTIFOR_ALLCONF=0.0_q
          ABN%EATOM_KIN(1:ABN%MIONS_ALLOCATE)=0.0_q
          ABN%EATOM_POT(1:ABN%MIONS_ALLOCATE)=0.0_q
          ABN%LADD=0
          ABN%LNITYP=0
          ABN%LNIONS=0
          ABN%LNTYP=0
          ABN%LSYS=0
          ABN%NITYP=0
          ABN%NSTEP=0
          ABN%OMEGA=0.0_q
          ABN%POSION=0.0_q
          ABN%QHEAT=0.0_q
          ABN%TIFOR=0.0_q
          ABN%TIFOR_MB=0.0_q
          ABN%TOTEN=0.0_q
          ABN%TOTEN_MB=0.0_q
          ABN%TSIF=0.0_q
          ABN%TSIF_MB=0.0_q
          ABN%TYPE=''
          ABN%VEL=0.0_q
        END SUBROUTINE ALLOCATE_ABN_ARRAY

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Allocate arrays for ab initio data
! TIFOR, TOTEN, TSIF for each component of interactiion are allocated
! even when that component is turned off for simplifying the program
! structure.
!****************************************************************************************************

        SUBROUTINE ALLOCATE_AB_ARRAY (MEM1,AB)
          USE ML_FF_STRUCT, ONLY: MEMORY, ABINITIO
          IMPLICIT NONE
! Input variables
          TYPE (MEMORY)       :: MEM1
          TYPE (ABINITIO)     :: AB
! Local variables
          REAL(q)             :: RMEM_LOC
! Initialization of memory usage
          RMEM_LOC=0.0_q
! Allocate
          IF(ALLOCATED(AB%A)) THEN
             DEALLOCATE(AB%A)
          ENDIF
          ALLOCATE(AB%A(1:3,1:3,1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+9*REAL(AB%MCONF_ALLOCATE,q)*8/1.0E+06_q
          IF(ALLOCATED(AB%B)) THEN
             DEALLOCATE(AB%B)
          ENDIF
          ALLOCATE(AB%B(1:3,1:3,1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+9*REAL(AB%MCONF_ALLOCATE,q)*8/1.0E+06_q
          IF(ALLOCATED(AB%CTIFOR_ALLCONF_OLD)) THEN
             DEALLOCATE(AB%CTIFOR_ALLCONF_OLD)
          ENDIF
          ALLOCATE(AB%CTIFOR_ALLCONF_OLD(1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MCONF_ALLOCATE,q)*8/1.0E+06_q
          IF(ALLOCATED(AB%CTIFOR_ALLCONF_NEW)) THEN
             DEALLOCATE(AB%CTIFOR_ALLCONF_NEW)
          ENDIF
          ALLOCATE(AB%CTIFOR_ALLCONF_NEW(1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MCONF_ALLOCATE,q)*8/1.0E+06_q
          IF(ALLOCATED(AB%ISYS)) THEN
             DEALLOCATE(AB%ISYS)
          ENDIF
          ALLOCATE(AB%ISYS(1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MCONF_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(AB%EATOM)) THEN
             DEALLOCATE(AB%EATOM)
          ENDIF
          ALLOCATE(AB%EATOM(1:AB%MTYP_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MTYP_ALLOCATE,q)*8/1.0E+06_q
          IF(ALLOCATED(AB%LADD)) THEN
             DEALLOCATE(AB%LADD)
          ENDIF
          ALLOCATE(AB%LADD(1:AB%MTYP_ALLOCATE,1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MTYP_ALLOCATE,q)*REAL(AB%MCONF_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(AB%LB)) THEN
             DEALLOCATE(AB%LB)
          ENDIF
          ALLOCATE(AB%LB(1:AB%MIONS_ALLOCATE,1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MIONS_ALLOCATE,q)*REAL(AB%MCONF_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(AB%LNITYP)) THEN
             DEALLOCATE(AB%LNITYP)
          ENDIF
          ALLOCATE(AB%LNITYP(1:AB%MIONS_ALLOCATE,1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MIONS_ALLOCATE,q)*REAL(AB%MCONF_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(AB%LNIONS)) THEN
             DEALLOCATE(AB%LNIONS)
          ENDIF
          ALLOCATE(AB%LNIONS(1:AB%MITYP_ALLOCATE,1:AB%MTYP_ALLOCATE,1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MITYP_ALLOCATE,q)*REAL(AB%MTYP_ALLOCATE,q)*REAL(AB%MCONF_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(AB%LNTYP)) THEN
             DEALLOCATE(AB%LNTYP)
          ENDIF
          ALLOCATE(AB%LNTYP(1:AB%MIONS_ALLOCATE,1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MIONS_ALLOCATE,q)*REAL(AB%MCONF_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(AB%LSYS)) THEN
             DEALLOCATE(AB%LSYS)
          ENDIF
          ALLOCATE(AB%LSYS(1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MCONF_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(AB%MITYP_CONF)) THEN
             DEALLOCATE(AB%MITYP_CONF)
          ENDIF
          ALLOCATE(AB%MITYP_CONF(1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MCONF_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(AB%NCONF_SYS)) THEN
             DEALLOCATE(AB%NCONF_SYS)
          ENDIF
          ALLOCATE(AB%NCONF_SYS(1:AB%NSYS_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%NSYS_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(AB%NITYP)) THEN
             DEALLOCATE(AB%NITYP)
          ENDIF
          ALLOCATE(AB%NITYP(1:AB%MTYP_ALLOCATE,1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MTYP_ALLOCATE,q)*REAL(AB%MCONF_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(AB%NIONS)) THEN
             DEALLOCATE(AB%NIONS)
          ENDIF
          ALLOCATE(AB%NIONS(1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MCONF_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(AB%NTYP)) THEN
             DEALLOCATE(AB%NTYP)
          ENDIF
          ALLOCATE(AB%NTYP(1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MCONF_ALLOCATE,q)/1.0E+06_q
          IF(ALLOCATED(AB%OMEGA)) THEN
             DEALLOCATE(AB%OMEGA)
          ENDIF
          ALLOCATE(AB%OMEGA(1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MCONF_ALLOCATE,q)*8/1.0E+06_q
          IF(ALLOCATED(AB%POMASS)) THEN
             DEALLOCATE(AB%POMASS)
          ENDIF
          ALLOCATE(AB%POMASS(1:AB%MTYP_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MTYP_ALLOCATE,q)*2/1.0E+06_q
          IF(ALLOCATED(AB%POSION)) THEN
             DEALLOCATE(AB%POSION)
          ENDIF
          ALLOCATE(AB%POSION(1:3,1:AB%MIONS_ALLOCATE,1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+3*REAL(AB%MIONS_ALLOCATE,q)*REAL(AB%MCONF_ALLOCATE,q)*8/1.0E+06_q
          IF(ALLOCATED(AB%SZNAM2)) THEN
             DEALLOCATE(AB%SZNAM2)
          ENDIF
          ALLOCATE(AB%SZNAM2(1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MCONF_ALLOCATE,q)*40/1.0E+06_q
          IF(ALLOCATED(AB%SZNAM2_ADD)) THEN
             DEALLOCATE(AB%SZNAM2_ADD)
          ENDIF
          ALLOCATE(AB%SZNAM2_ADD(1:AB%NSYS_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MCONF_ALLOCATE,q)*40/1.0E+06_q
          IF(ALLOCATED(AB%TIFOR)) THEN
             DEALLOCATE(AB%TIFOR)
          ENDIF
          ALLOCATE(AB%TIFOR(1:3,1:AB%MIONS_ALLOCATE,1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+3*REAL(AB%MIONS_ALLOCATE,q)*REAL(AB%MCONF_ALLOCATE,q)*8/1.0E+06_q
          IF(ALLOCATED(AB%TIFOR_MB)) THEN
             DEALLOCATE(AB%TIFOR_MB)
          ENDIF
          ALLOCATE(AB%TIFOR_MB(1:3,1:AB%MIONS_ALLOCATE,1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+3*REAL(AB%MIONS_ALLOCATE,q)*REAL(AB%MCONF_ALLOCATE,q)*8/1.0E+06_q
          IF(ALLOCATED(AB%TOTEN)) THEN
             DEALLOCATE(AB%TOTEN)
          ENDIF
          ALLOCATE(AB%TOTEN(1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MCONF_ALLOCATE,q)*8/1.0E+06_q
          IF(ALLOCATED(AB%TOTEN_MB)) THEN
             DEALLOCATE(AB%TOTEN_MB)
          ENDIF
          ALLOCATE(AB%TOTEN_MB(1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MCONF_ALLOCATE,q)*8/1.0E+06_q
          IF(ALLOCATED(AB%TSIF)) THEN
             DEALLOCATE(AB%TSIF)
          ENDIF
          ALLOCATE(AB%TSIF(1:3,1:3,1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MCONF_ALLOCATE,q)*9*8/1.0E+06_q
          IF(ALLOCATED(AB%TSIF_MB)) THEN
             DEALLOCATE(AB%TSIF_MB)
          ENDIF
          ALLOCATE(AB%TSIF_MB(1:3,1:3,1:AB%MCONF_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MCONF_ALLOCATE,q)*9*8/1.0E+06_q
          IF(ALLOCATED(AB%TYPE)) THEN
             DEALLOCATE(AB%TYPE)
          ENDIF
          ALLOCATE(AB%TYPE(1:AB%MTYP_ALLOCATE))
          RMEM_LOC=RMEM_LOC+REAL(AB%MTYP_ALLOCATE,q)*2/1.0E+06_q
          MEM1%RMEM(3)=MAX(MEM1%RMEM(3),RMEM_LOC)
! Initialization of variables
          AB%A=0.0_q
          AB%B=0.0_q
          AB%CTIFOR_ALLCONF_OLD=0.0_q
          AB%CTIFOR_ALLCONF_NEW=0.0_q
          AB%LADD=0
          AB%LB=0
          AB%LNITYP=0
          AB%LNIONS=0
          AB%LNTYP=0
          AB%LSYS=0
          AB%NCONF_SYS=0
          AB%NTYP=0
          AB%NITYP=0
          AB%NIONS=0
          AB%POMASS=0.0_q
          AB%POSION=0.0_q
          AB%TIFOR=0.0_q
          AB%TIFOR_MB=0.0_q
          AB%TOTEN=0.0_q
          AB%TOTEN_MB=0.0_q
          AB%TSIF=0.0_q
          AB%TSIF_MB=0.0_q
          AB%TYPE=''
          AB%OMEGA=0.0_q
        END SUBROUTINE ALLOCATE_AB_ARRAY
#endif

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Deallocate data storing the ab initio data
!****************************************************************************************************

        SUBROUTINE DEALLOCATE_AB_ARRAY(AB)
          USE ML_FF_STRUCT, ONLY: ABINITIO
          IMPLICIT NONE
          TYPE (ABINITIO) :: AB
          IF(ALLOCATED(AB%A)) THEN
             DEALLOCATE(AB%A)
          ENDIF
          IF(ALLOCATED(AB%B)) THEN
             DEALLOCATE(AB%B)
          ENDIF
          IF(ALLOCATED(AB%CTIFOR_ALLCONF_OLD)) THEN
             DEALLOCATE(AB%CTIFOR_ALLCONF_OLD)
          ENDIF
          IF(ALLOCATED(AB%CTIFOR_ALLCONF_NEW)) THEN
             DEALLOCATE(AB%CTIFOR_ALLCONF_NEW)
          ENDIF
          IF(ALLOCATED(AB%EATOM)) THEN
             DEALLOCATE(AB%EATOM)
          ENDIF
          IF(ALLOCATED(AB%ISYS)) THEN
             DEALLOCATE(AB%ISYS)
          ENDIF
          IF(ALLOCATED(AB%LADD)) THEN
             DEALLOCATE(AB%LADD)
          ENDIF
          IF(ALLOCATED(AB%LB)) THEN
             DEALLOCATE(AB%LB)
          ENDIF
          IF(ALLOCATED(AB%LNITYP)) THEN
             DEALLOCATE(AB%LNITYP)
          ENDIF
          IF(ALLOCATED(AB%LNIONS)) THEN
             DEALLOCATE(AB%LNIONS)
          ENDIF
          IF(ALLOCATED(AB%LNTYP)) THEN
             DEALLOCATE(AB%LNTYP)
          ENDIF
          IF(ALLOCATED(AB%LSYS)) THEN
             DEALLOCATE(AB%LSYS)
          ENDIF
          IF(ALLOCATED(AB%MITYP_CONF)) THEN
             DEALLOCATE(AB%MITYP_CONF)
          ENDIF
          IF(ALLOCATED(AB%NCONF_SYS)) THEN
             DEALLOCATE(AB%NCONF_SYS)
          ENDIF
          IF(ALLOCATED(AB%NITYP)) THEN
             DEALLOCATE(AB%NITYP)
          ENDIF
          IF(ALLOCATED(AB%NIONS)) THEN
             DEALLOCATE(AB%NIONS)
          ENDIF
          IF(ALLOCATED(AB%NTYP)) THEN
             DEALLOCATE(AB%NTYP)
          ENDIF
          IF(ALLOCATED(AB%OMEGA)) THEN
             DEALLOCATE(AB%OMEGA)
          ENDIF
          IF(ALLOCATED(AB%POMASS)) THEN
             DEALLOCATE(AB%POMASS)
          ENDIF
          IF(ALLOCATED(AB%POSION)) THEN
             DEALLOCATE(AB%POSION)
          ENDIF
          IF(ALLOCATED(AB%SZNAM2)) THEN
             DEALLOCATE(AB%SZNAM2)
          ENDIF
          IF(ALLOCATED(AB%SZNAM2_ADD)) THEN
             DEALLOCATE(AB%SZNAM2_ADD)
          ENDIF
          IF(ALLOCATED(AB%TIFOR)) THEN
             DEALLOCATE(AB%TIFOR)
          ENDIF
          IF(ALLOCATED(AB%TIFOR_MB)) THEN
             DEALLOCATE(AB%TIFOR_MB)
          ENDIF
          IF(ALLOCATED(AB%TOTEN)) THEN
             DEALLOCATE(AB%TOTEN)
          ENDIF
          IF(ALLOCATED(AB%TOTEN_MB)) THEN
             DEALLOCATE(AB%TOTEN_MB)
          ENDIF
          IF(ALLOCATED(AB%TSIF)) THEN
             DEALLOCATE(AB%TSIF)
          ENDIF
          IF(ALLOCATED(AB%TSIF_MB)) THEN
             DEALLOCATE(AB%TSIF_MB)
          ENDIF
          IF(ALLOCATED(AB%TYPE)) THEN
             DEALLOCATE(AB%TYPE)
          ENDIF
        END SUBROUTINE DEALLOCATE_AB_ARRAY
#endif

!****************************************************************************************************
! Deallocate data storing the new ab initio data
!****************************************************************************************************

        SUBROUTINE DEALLOCATE_ABN_ARRAY(ABN)
          USE ML_FF_STRUCT, ONLY: ABINITIONEW
          IMPLICIT NONE
! Input variables
          TYPE (ABINITIONEW)  :: ABN
          IF(ALLOCATED(ABN%A)) THEN
             DEALLOCATE(ABN%A)
          ENDIF
          IF(ALLOCATED(ABN%B)) THEN
             DEALLOCATE(ABN%B)
          ENDIF
          IF(ALLOCATED(ABN%BE)) THEN
             DEALLOCATE(ABN%BE)
          ENDIF
          IF(ALLOCATED(ABN%BEAV)) THEN
             DEALLOCATE(ABN%BEAV)
          ENDIF
          IF(ALLOCATED(ABN%BEMAX)) THEN
             DEALLOCATE(ABN%BEMAX)
          ENDIF
          IF(ALLOCATED(ABN%CTIFOR_ALLCONF)) THEN
             DEALLOCATE(ABN%CTIFOR_ALLCONF)
          ENDIF
          IF(ALLOCATED(ABN%EATOM_KIN)) THEN
             DEALLOCATE(ABN%EATOM_KIN)
          ENDIF
          IF(ALLOCATED(ABN%EATOM_POT)) THEN
             DEALLOCATE(ABN%EATOM_POT)
          ENDIF
          IF(ALLOCATED(ABN%LADD)) THEN
             DEALLOCATE(ABN%LADD)
          ENDIF
          IF(ALLOCATED(ABN%LNITYP)) THEN
             DEALLOCATE(ABN%LNITYP)
          ENDIF
          IF(ALLOCATED(ABN%LNIONS)) THEN
             DEALLOCATE(ABN%LNIONS)
          ENDIF
          IF(ALLOCATED(ABN%LNTYP)) THEN
             DEALLOCATE(ABN%LNTYP)
          ENDIF
          IF(ALLOCATED(ABN%LSYS)) THEN
             DEALLOCATE(ABN%LSYS)
          ENDIF
          IF(ALLOCATED(ABN%NITYP)) THEN
             DEALLOCATE(ABN%NITYP)
          ENDIF
          IF(ALLOCATED(ABN%NSTEP)) THEN
             DEALLOCATE(ABN%NSTEP)
          ENDIF
          IF(ALLOCATED(ABN%OMEGA)) THEN
             DEALLOCATE(ABN%OMEGA)
          ENDIF
          IF(ALLOCATED(ABN%POMASS)) THEN
             DEALLOCATE(ABN%POMASS)
          ENDIF
          IF(ALLOCATED(ABN%POSION)) THEN
             DEALLOCATE(ABN%POSION)
          ENDIF
          IF(ALLOCATED(ABN%QHEAT)) THEN
             DEALLOCATE(ABN%QHEAT)
          ENDIF
          IF(ALLOCATED(ABN%SZNAM2_ADD)) THEN
             DEALLOCATE(ABN%SZNAM2_ADD)
          ENDIF
          IF(ALLOCATED(ABN%TIFOR)) THEN
             DEALLOCATE(ABN%TIFOR)
          ENDIF
          IF(ALLOCATED(ABN%TIFOR_MB)) THEN
             DEALLOCATE(ABN%TIFOR_MB)
          ENDIF
          IF(ALLOCATED(ABN%TOTEN)) THEN
             DEALLOCATE(ABN%TOTEN)
          ENDIF
          IF(ALLOCATED(ABN%TOTEN_MB)) THEN
             DEALLOCATE(ABN%TOTEN_MB)
          ENDIF
          IF(ALLOCATED(ABN%TSIF)) THEN
             DEALLOCATE(ABN%TSIF)
          ENDIF
          IF(ALLOCATED(ABN%TSIF_MB)) THEN
             DEALLOCATE(ABN%TSIF_MB)
          ENDIF
          IF(ALLOCATED(ABN%TYPE)) THEN
             DEALLOCATE(ABN%TYPE)
          ENDIF
          IF(ALLOCATED(ABN%VEL)) THEN
             DEALLOCATE(ABN%VEL)
          ENDIF
        END SUBROUTINE DEALLOCATE_ABN_ARRAY

#ifndef ML_LEARNING_NOT_ACTIVE
        SUBROUTINE CONSTRUCT_AB_FOR_HELPING_STRUCTURE(ML_SUPER_HANDLE,MLH_TMP)
           USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE
           IMPLICIT NONE
           ! Input variables
           TYPE (ML_SUPER_TYPE) :: ML_SUPER_HANDLE
           ! Output variables
           TYPE (ML_SUPER_TYPE) :: MLH_TMP
! Set scalar values
           MLH_TMP%AB%MIONS=ML_SUPER_HANDLE%AB%MIONS
           MLH_TMP%AB%MTYP=ML_SUPER_HANDLE%AB%MTYP
           MLH_TMP%AB%NSYS=ML_SUPER_HANDLE%AB%NSYS
           MLH_TMP%AB%NCONF=ML_SUPER_HANDLE%AB%NCONF
           MLH_TMP%AB%MIONS_ALLOCATE=ML_SUPER_HANDLE%AB%MIONS_ALLOCATE
           MLH_TMP%AB%MTYP_ALLOCATE=ML_SUPER_HANDLE%AB%MTYP_ALLOCATE
           MLH_TMP%AB%MITYP_ALLOCATE=ML_SUPER_HANDLE%AB%MITYP_ALLOCATE
           MLH_TMP%AB%MCONF_ALLOCATE=ML_SUPER_HANDLE%AB%MCONF_ALLOCATE
           MLH_TMP%AB%NSYS_ALLOCATE=ML_SUPER_HANDLE%AB%NSYS_ALLOCATE
! Allocate arrays
           CALL ALLOCATE_AB_ARRAY(MLH_TMP%MEM1,MLH_TMP%AB)
! Add memory consumption for helping arrays
           ML_SUPER_HANDLE%MEM1%RMEM(3)=2*ML_SUPER_HANDLE%MEM1%RMEM(3)
! Set helping arrays
           MLH_TMP%AB%A(1:3,1:3,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%A(1:3,1:3,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%B(1:3,1:3,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%B(1:3,1:3,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%ISYS(1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%ISYS(1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%EATOM(1:ML_SUPER_HANDLE%AB%MTYP_ALLOCATE)=ML_SUPER_HANDLE%AB%EATOM(1:ML_SUPER_HANDLE%AB%MTYP_ALLOCATE)
           MLH_TMP%AB%LADD(1:ML_SUPER_HANDLE%AB%MTYP_ALLOCATE,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%LADD &
                          (1:ML_SUPER_HANDLE%AB%MTYP_ALLOCATE,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%LB(1:ML_SUPER_HANDLE%AB%MIONS_ALLOCATE,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%LB &
                        (1:ML_SUPER_HANDLE%AB%MIONS_ALLOCATE,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%LNITYP(1:ML_SUPER_HANDLE%AB%MIONS_ALLOCATE,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)= &
                      ML_SUPER_HANDLE%AB%LNITYP(1:ML_SUPER_HANDLE%AB%MIONS_ALLOCATE,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%LNIONS(1:ML_SUPER_HANDLE%AB%MITYP_ALLOCATE,1:ML_SUPER_HANDLE%AB%MTYP_ALLOCATE, & 
                             1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%LNIONS(1:ML_SUPER_HANDLE%AB%MITYP_ALLOCATE, &
                             1:ML_SUPER_HANDLE%AB%MTYP_ALLOCATE,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%LNTYP(1:ML_SUPER_HANDLE%AB%MIONS_ALLOCATE,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%LNTYP &
                           (1:ML_SUPER_HANDLE%AB%MIONS_ALLOCATE,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%LSYS(1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%LSYS(1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%NCONF_SYS(1:ML_SUPER_HANDLE%AB%NSYS_ALLOCATE)=ML_SUPER_HANDLE%AB%NCONF_SYS(1:ML_SUPER_HANDLE%AB%NSYS_ALLOCATE)
           MLH_TMP%AB%NITYP(1:ML_SUPER_HANDLE%AB%MTYP_ALLOCATE,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%NITYP &
                           (1:ML_SUPER_HANDLE%AB%MTYP_ALLOCATE,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%NIONS(1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%NIONS(1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%NTYP(1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%NTYP(1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%OMEGA(1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%OMEGA(1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%POMASS(1:ML_SUPER_HANDLE%AB%MTYP_ALLOCATE)=ML_SUPER_HANDLE%AB%POMASS(1:ML_SUPER_HANDLE%AB%MTYP_ALLOCATE)
           MLH_TMP%AB%POSION(1:3,1:ML_SUPER_HANDLE%AB%MIONS_ALLOCATE,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%POSION &
                            (1:3,1:ML_SUPER_HANDLE%AB%MIONS_ALLOCATE,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%SZNAM2(1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%SZNAM2(1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%SZNAM2_ADD(1:ML_SUPER_HANDLE%AB%NSYS_ALLOCATE)=ML_SUPER_HANDLE%AB%SZNAM2_ADD(1:ML_SUPER_HANDLE%AB%NSYS_ALLOCATE)
           MLH_TMP%AB%TIFOR(1:3,1:ML_SUPER_HANDLE%AB%MIONS_ALLOCATE,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%TIFOR &
                           (1:3,1:ML_SUPER_HANDLE%AB%MIONS_ALLOCATE,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%TIFOR_MB(1:3,1:ML_SUPER_HANDLE%AB%MIONS_ALLOCATE,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%TIFOR_MB &
                              (1:3,1:ML_SUPER_HANDLE%AB%MIONS_ALLOCATE,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%TOTEN(1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%TOTEN(1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%TOTEN_MB(1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%TOTEN_MB(1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%TSIF(1:3,1:3,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%TSIF(1:3,1:3,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%TSIF_MB(1:3,1:3,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)=ML_SUPER_HANDLE%AB%TSIF_MB(1:3,1:3,1:ML_SUPER_HANDLE%AB%MCONF_ALLOCATE)
           MLH_TMP%AB%TYPE(1:ML_SUPER_HANDLE%AB%MTYP_ALLOCATE)=ML_SUPER_HANDLE%AB%TYPE(1:ML_SUPER_HANDLE%AB%MTYP_ALLOCATE)
        END SUBROUTINE CONSTRUCT_AB_FOR_HELPING_STRUCTURE
#endif

      END MODULE ABINITIO_DATA
#endif
