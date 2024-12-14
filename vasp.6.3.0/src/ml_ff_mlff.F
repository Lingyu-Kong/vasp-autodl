#if defined(MPI) || defined(MPI_CHAIN)
#include "symbol.inc"

MODULE ml_main_subroutines
!****************************************************************************************************
! Module containing subroutines executing following procedures
! 1) Initializing machine-learning force field generation,
! 2) Providing estimations on energy, forces, stress and their errors,
! 3) Generating force field,
! 4) Finalizing machine-learning force field generation,
! The module is used in ml_interface.F file, and subroutines are called
! also from ml_interface.F file.
!****************************************************************************************************
        CONTAINS
!****************************************************************************************************
! Initialization of the machine-learning force field generations
!****************************************************************************************************

        SUBROUTINE MACHINE_LEARNING_INIT_LIB (ML_SUPER_HANDLE)
          USE ML_FF_PREC
          USE ML_FF_CONSTANT
          USE LOGFILE
          USE IOHANDLE
          USE MPI_DATA
          USE TIME_CHECK
          USE ABINITIO_DATA
          USE FORCE_FIELD
          USE MEMORY_USAGE
          USE RANDOM_NUMBER_GENERATOR
          USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE, ML_FILESTREAM_GLOBAL, ML_IO_WRITE
          USE PRINT2TERMINAL
          IMPLICIT NONE
          TYPE (ML_SUPER_TYPE) :: ML_SUPER_HANDLE
! Local variable
          TYPE (ML_SUPER_TYPE) :: MLH_TMP
          INTEGER              :: IERR
          CHARACTER(LEN=4)     :: FILE_POSTFIX
          REAL(q)              :: RDUM
! Initialize timing.
          CALL INIT_TIME(ML_SUPER_HANDLE%TIM)
! Initialize MPI for all MPI processes
          CALL MPI_MAKE_NEW_COMMUNICATOR(ML_SUPER_HANDLE%ML_INPUT_HANDLE%COMM_VASP,ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,0, &
                                         ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM, &
                                         ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD%NCPU,ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD%NODE_ME)
#ifdef use_shmem
! Make intra and inter node communicators
          CALL M_divide_intra_inter_node(ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD,ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD, &
           ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_intra_node_world,ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_inter_node_world)
#endif
#ifdef scaLAPACK
! Initialize Scalapack process grid for NxM processors
          CALL SCALAPACK_INIT(ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD%NCPU,ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD, &
                              ML_SUPER_HANDLE%ML_INPUT_HANDLE%NDIM_SCALAPACK_FF, &
                              ML_SUPER_HANDLE%PAR_SUP_HANDLE%CONTEXT_WORLD%NDIM, &
                              ML_SUPER_HANDLE%PAR_SUP_HANDLE%CONTEXT_WORLD%ICTXT,ML_SUPER_HANDLE%PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,  &
                              ML_SUPER_HANDLE%PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,ML_SUPER_HANDLE%PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW, &
                              ML_SUPER_HANDLE%PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL,ML_SUPER_HANDLE%PAR_SUP_HANDLE%CONTEXT_WORLD%NPROCS)
! Initialize Scalapack process grid for 1 processor (mainly used for testing)
          CALL SCALAPACK_INIT(1,ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NDIM_SCALAPACK_FF, &
                              ML_SUPER_HANDLE%PAR_SUP_HANDLE%CONTEXT_SINGLE%NDIM, &
                              ML_SUPER_HANDLE%PAR_SUP_HANDLE%CONTEXT_SINGLE%ICTXT, &
                              ML_SUPER_HANDLE%PAR_SUP_HANDLE%CONTEXT_SINGLE%NP_ROW,ML_SUPER_HANDLE%PAR_SUP_HANDLE%CONTEXT_SINGLE%NP_COL, &
                              ML_SUPER_HANDLE%PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_ROW, &
                              ML_SUPER_HANDLE%PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_COL,ML_SUPER_HANDLE%PAR_SUP_HANDLE%CONTEXT_SINGLE%NPROCS)
#endif
! Check here if input parameter combinations work together
          CALL CHECK_PARAMTER_COMPATIBILITY(ML_SUPER_HANDLE%FF)
! Set the postfix for the current instance of machine learning
          CALL ML_SET_POSTFIX(ML_SUPER_HANDLE%INSTANCE,FILE_POSTFIX)
! Define file names
#ifndef ML_LEARNING_NOT_ACTIVE
          IF (ML_SUPER_HANDLE%FF%ISTART.NE.2) THEN
             ! New ab-initio data file
             CALL ML_FILE_STREAM_INCREMENT(ML_FILESTREAM_GLOBAL)
             IF (ML_SUPER_HANDLE%INSTANCE.EQ.1) THEN
                CALL SET_FILE ('ML_ABN',ML_FILESTREAM_GLOBAL,ML_SUPER_HANDLE%AB%ABNCAR ,ML_SUPER_HANDLE%AB%IU2)
             ELSE
                CALL SET_FILE ('ML_ABN'//TRIM(FILE_POSTFIX),ML_FILESTREAM_GLOBAL,ML_SUPER_HANDLE%AB%ABNCAR ,ML_SUPER_HANDLE%AB%IU2)
             ENDIF
             ! Output file summarizing New FF parameters
             CALL ML_FILE_STREAM_INCREMENT(ML_FILESTREAM_GLOBAL)
             IF (ML_SUPER_HANDLE%INSTANCE.EQ.1) THEN
                CALL SET_FILE ('ML_FFN',ML_FILESTREAM_GLOBAL,ML_SUPER_HANDLE%FF%FFNCAR ,ML_SUPER_HANDLE%FF%IU2)
             ELSE
                CALL SET_FILE ('ML_FFN'//TRIM(FILE_POSTFIX),ML_FILESTREAM_GLOBAL,ML_SUPER_HANDLE%FF%FFNCAR ,ML_SUPER_HANDLE%FF%IU2)
             ENDIF
             ! Old ab-initio data file
             CALL ML_FILE_STREAM_INCREMENT(ML_FILESTREAM_GLOBAL)
             IF (ML_SUPER_HANDLE%INSTANCE.EQ.1) THEN
                CALL SET_FILE ('ML_AB',ML_FILESTREAM_GLOBAL,ML_SUPER_HANDLE%AB%ABCAR  ,ML_SUPER_HANDLE%AB%IU1)
             ELSE
                CALL SET_FILE ('ML_AB'//TRIM(FILE_POSTFIX),ML_FILESTREAM_GLOBAL,ML_SUPER_HANDLE%AB%ABCAR  ,ML_SUPER_HANDLE%AB%IU1)
             ENDIF
          ENDIF
#endif
  
          ! Output file summzrizing regression results
          CALL ML_FILE_STREAM_INCREMENT(ML_FILESTREAM_GLOBAL)
          IF (ML_SUPER_HANDLE%INSTANCE.EQ.1) THEN
             CALL SET_FILE ('ML_REG',ML_FILESTREAM_GLOBAL,ML_SUPER_HANDLE%AB%REGCAR ,ML_SUPER_HANDLE%AB%IU3)
          ELSE
             CALL SET_FILE ('ML_REG'//TRIM(FILE_POSTFIX),ML_FILESTREAM_GLOBAL,ML_SUPER_HANDLE%AB%REGCAR ,ML_SUPER_HANDLE%AB%IU3)
          ENDIF
          ! Output file summarizing the histogram data
          CALL ML_FILE_STREAM_INCREMENT(ML_FILESTREAM_GLOBAL)
          IF (ML_SUPER_HANDLE%INSTANCE.EQ.1) THEN
             CALL SET_FILE ('ML_HIS',ML_FILESTREAM_GLOBAL,ML_SUPER_HANDLE%AB%HISCAR ,ML_SUPER_HANDLE%AB%IU4)
          ELSE
             CALL SET_FILE ('ML_HIS'//TRIM(FILE_POSTFIX),ML_FILESTREAM_GLOBAL,ML_SUPER_HANDLE%AB%HISCAR ,ML_SUPER_HANDLE%AB%IU4)
          ENDIF
          ! Input file summarizing New FF parameters
          CALL ML_FILE_STREAM_INCREMENT(ML_FILESTREAM_GLOBAL)
          IF (ML_SUPER_HANDLE%INSTANCE.EQ.1) THEN
             CALL SET_FILE ('ML_FF',ML_FILESTREAM_GLOBAL,ML_SUPER_HANDLE%FF%FFCAR  ,ML_SUPER_HANDLE%FF%IU1)
          ELSE
             CALL SET_FILE ('ML_FF'//TRIM(FILE_POSTFIX),ML_FILESTREAM_GLOBAL,ML_SUPER_HANDLE%FF%FFCAR  ,ML_SUPER_HANDLE%FF%IU1)
          ENDIF
          ! Output file !summarizing several important information on ML calculations
          CALL ML_FILE_STREAM_INCREMENT(ML_FILESTREAM_GLOBAL)
          IF (ML_SUPER_HANDLE%INSTANCE.EQ.1) THEN
             CALL SET_FILE ('ML_LOGFILE',ML_FILESTREAM_GLOBAL,ML_SUPER_HANDLE%LOGF%ML_LOGFILE,ML_SUPER_HANDLE%LOGF%IU)
          ELSE
             CALL SET_FILE ('ML_LOGFILE'//TRIM(FILE_POSTFIX),ML_FILESTREAM_GLOBAL,ML_SUPER_HANDLE%LOGF%ML_LOGFILE,ML_SUPER_HANDLE%LOGF%IU)
          ENDIF
          ! Output file including local atomic energies.
          IF(ML_SUPER_HANDLE%FFM%LEATOM) THEN
             CALL ML_FILE_STREAM_INCREMENT(ML_FILESTREAM_GLOBAL)
             IF (ML_SUPER_HANDLE%INSTANCE.EQ.1) THEN
                CALL SET_FILE ('ML_EATOM',ML_FILESTREAM_GLOBAL,ML_SUPER_HANDLE%FF%ML_EATOM,ML_SUPER_HANDLE%FF%IU3)
             ELSE
                CALL SET_FILE ('ML_EATOM'//TRIM(FILE_POSTFIX),ML_FILESTREAM_GLOBAL,ML_SUPER_HANDLE%FF%ML_EATOM,ML_SUPER_HANDLE%FF%IU3)
             ENDIF
          ENDIF
          ! Output file including local heat flux.
          IF(ML_SUPER_HANDLE%FFM%LHEAT) THEN
             CALL ML_FILE_STREAM_INCREMENT(ML_FILESTREAM_GLOBAL)
             IF (ML_SUPER_HANDLE%INSTANCE.EQ.1) THEN
                CALL SET_FILE ('ML_HEAT',ML_FILESTREAM_GLOBAL,ML_SUPER_HANDLE%FF%ML_HEAT,ML_SUPER_HANDLE%FF%IU4)
             ELSE
                CALL SET_FILE ('ML_HEAT'//TRIM(FILE_POSTFIX),ML_FILESTREAM_GLOBAL,ML_SUPER_HANDLE%FF%ML_HEAT,ML_SUPER_HANDLE%FF%IU4)
             ENDIF
          ENDIF
! Open files.
          OPEN(UNIT=ML_SUPER_HANDLE%LOGF%IU,FILE=ML_SUPER_HANDLE%ML_INPUT_HANDLE%DIR_APP(1:ML_SUPER_HANDLE%ML_INPUT_HANDLE%DIR_LEN) &
                                            //ML_SUPER_HANDLE%LOGF%ML_LOGFILE,STATUS='UNKNOWN')
          IF (ML_IO_WRITE) WRITE(ML_SUPER_HANDLE%LOGF%IU, '(A)') "ML_LOGFILE version 0.1.0"
          IF(ML_SUPER_HANDLE%FFM%LEATOM) THEN
             OPEN(UNIT=ML_SUPER_HANDLE%FF%IU3,FILE=ML_SUPER_HANDLE%ML_INPUT_HANDLE%DIR_APP(1:ML_SUPER_HANDLE%ML_INPUT_HANDLE%DIR_LEN) &
                                                                           //ML_SUPER_HANDLE%FF%ML_EATOM,STATUS='UNKNOWN')
          ENDIF
          IF(ML_SUPER_HANDLE%FFM%LHEAT) THEN
             OPEN(UNIT=ML_SUPER_HANDLE%FF%IU4,FILE=ML_SUPER_HANDLE%ML_INPUT_HANDLE%DIR_APP(1:ML_SUPER_HANDLE%ML_INPUT_HANDLE%DIR_LEN) &
                                                                           //ML_SUPER_HANDLE%FF%ML_HEAT,STATUS='UNKNOWN')
          ENDIF
! Set Version number of ML_FF (ML_FFN file)
          ML_SUPER_HANDLE%FF%ML_FF_VERSION_WRITE(1) = 0
          ML_SUPER_HANDLE%FF%ML_FF_VERSION_WRITE(2) = 1
          ML_SUPER_HANDLE%FF%ML_FF_VERSION_WRITE(3) = 0
! Set contents of ML_LOGFILE.
          ML_SUPER_HANDLE%LOGF%NWRITE=ML_SUPER_HANDLE%ML_INPUT_HANDLE%NWRITE_FF
! Set this control variable to .FALSE.
          ML_SUPER_HANDLE%FF%LSECOND_TIME_ML_TO_VASP = .FALSE.
! Set another control variable to false
          ML_SUPER_HANDLE%FF%LSECOND_TIME_ML_TO_VASP = .FALSE.
! Random seed for machine learning is initialized here
          CALL RANE_ION(RDUM,PUT=ML_SUPER_HANDLE%FF%SEED(1:))
! Initialize a variable that is the same for both drivers
          ML_SUPER_HANDLE%FF%FORCING_UPDATE=.FALSE.
! Initialize counter for force field writing
          ML_SUPER_HANDLE%FF%NCOUNT_GENFF=0
! Do the actual initializations here
          IF (ML_SUPER_HANDLE%FF%ISTART.EQ.3) THEN
#ifndef ML_LEARNING_NOT_ACTIVE
             CALL ISTART_3_DRIVER(ML_SUPER_HANDLE,MLH_TMP)
#endif
          ELSE
             CALL ISTART_0_1_2_DRIVER(ML_SUPER_HANDLE)
          ENDIF
        END SUBROUTINE MACHINE_LEARNING_INIT_LIB

!****************************************************************************************************
! Subroutine to predict energy, forces and stress tensor by using
! existing force field.
!****************************************************************************************************

        SUBROUTINE ML_TO_VASP_MACHINE_LEARNING_LIB (ML_SUPER_HANDLE,A,B,POTIM,NSTEP,OMEGA,POSION,VEL,SZNAM2,TOTEN,TIFOR,TSIF,LDO_AB_INITIO)
          USE ABINITIO_DATA
          USE ML_FF_CONSTANT
          USE FORCE_FIELD
          USE MPI_DATA
          USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE, ML_IO_WRITE
          USE LOGFILE
          IMPLICIT NONE
          TYPE (ML_SUPER_TYPE) :: ML_SUPER_HANDLE
! Input variables
!          REAL(q)          , INTENT(IN) :: A(1:,1:) !(1:3,1:3)
          REAL(q)  :: A(1:,1:) !(1:3,1:3)
          REAL(q)          , INTENT(IN) :: B(1:,1:) !(1:3,1:3)
          INTEGER          , INTENT(IN) :: NSTEP
          REAL(q)          , INTENT(IN) :: POSION(:,:) !(1:3,1:ML_SUPER_HANDLE%ABN%NIONS)
          REAL(q)          , INTENT(IN) :: POTIM
          REAL(q)          , INTENT(IN) :: OMEGA
          CHARACTER(LEN=40), INTENT(IN) :: SZNAM2
          REAL(q)          , INTENT(IN) :: VEL(:,:) !(1:3,1:ML_SUPER_HANDLE%ABN%NIONS)
! Output variables
          LOGICAL          , INTENT(OUT) :: LDO_AB_INITIO
          REAL(q)          , INTENT(OUT) :: TOTEN
          REAL(q)          , INTENT(OUT) :: TIFOR(:,:) !(1:3,1:ML_SUPER_HANDLE%ABN%NIONS)
          REAL(q)          , INTENT(OUT) :: TSIF(:,:) !(1:3,1:3)
! Local variable
          INTEGER                        :: IERR
          INTEGER                        :: INIONS,IXYZ
          INTEGER                        :: NCONF_AB_DUMMY
          INTEGER                        :: NCONF_ABN_DUMMY
          INTEGER                        :: NCONF_SYS_ABN_DUMMY
          INTEGER                        :: ICONF

          ! Need to set MD step in FF because logging is using FF%NSTEP.
          ML_SUPER_HANDLE%FF%NSTEP = NSTEP
! Configure new ab initio data.
          CALL SET_ABN_ST (ML_SUPER_HANDLE%TIM,ML_SUPER_HANDLE%ABN,A(1:,1:),B(1:,1:),OMEGA,POSION,POTIM,VEL)
! Set FF%FORCING_UPDATE=.TRUE. if it is the last ionic step
          IF (ML_SUPER_HANDLE%FF%ISTART.EQ.2) THEN
             ML_SUPER_HANDLE%FF%FORCING_UPDATE=.FALSE.
          ENDIF
! Predict energy, force and stress.
! When ML_SUPER_HANDLE%FF%ISTART is set to 2, the program does not calculate NCONF_SYS
! from ML_FF_ABCAR file because it is not necessary for the force field 
! only calculations. In this case, they are simply set to 1.
          IF(ML_SUPER_HANDLE%FF%ISTART.EQ.2) THEN
             NCONF_AB_DUMMY=1
             NCONF_ABN_DUMMY=1
             NCONF_SYS_ABN_DUMMY=1
! Otherwise, actual data are used.
          ELSE
#ifndef ML_LEARNING_NOT_ACTIVE
             NCONF_AB_DUMMY=ML_SUPER_HANDLE%AB%NCONF
             NCONF_ABN_DUMMY=ML_SUPER_HANDLE%ABN%NCONF
             NCONF_SYS_ABN_DUMMY=ML_SUPER_HANDLE%AB%NCONF_SYS(ML_SUPER_HANDLE%ABN%ISYS)
#endif
          ENDIF
! Actual prediction of force stress and error, MANYBODY is the main routine.
          CALL PES_FF (ML_SUPER_HANDLE%AB,ML_SUPER_HANDLE%ABN,ML_SUPER_HANDLE%TIM,ML_SUPER_HANDLE%PAR_SUP_HANDLE, &
               ML_SUPER_HANDLE%MLSHM,ML_SUPER_HANDLE%MEM1,ML_SUPER_HANDLE%LOGF,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM, &
               LDO_AB_INITIO,NCONF_AB_DUMMY,NCONF_ABN_DUMMY,NCONF_SYS_ABN_DUMMY,NSTEP)
! Output atomic energies, if necessary.
          IF(ML_SUPER_HANDLE%FFM%LEATOM.AND.ML_SUPER_HANDLE%FF%IFF) THEN
             CALL OUT_EATOM (ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD, &
                  ML_SUPER_HANDLE%ABN%A(1:,1:,ML_SUPER_HANDLE%ABN%NCONF),ML_SUPER_HANDLE%ABN%EATOM_KIN(1:), &
                  ML_SUPER_HANDLE%ABN%EATOM_POT(1:),ML_SUPER_HANDLE%FF%IU3,ML_SUPER_HANDLE%ABN%LADD(1:),ML_SUPER_HANDLE%ABN%LNIONS(1:,1:), &
                  ML_SUPER_HANDLE%ABN%MIONS,ML_SUPER_HANDLE%ABN%MITYP,ML_SUPER_HANDLE%ABN%MTYP,ML_SUPER_HANDLE%ABN%NIONS, &
                  ML_SUPER_HANDLE%ABN%NITYP(1:),NSTEP,ML_SUPER_HANDLE%ABN%NTYP, &
                  ML_SUPER_HANDLE%ABN%POSION(1:,1:,ML_SUPER_HANDLE%ABN%NCONF),ML_SUPER_HANDLE%ABN%TYPE(1:))
          ENDIF
! Output heat flux, if necessary.
          IF(ML_SUPER_HANDLE%FFM%LHEAT.AND.ML_SUPER_HANDLE%FF%IFF) THEN
             CALL OUT_HEAT (ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD, &
                  ML_SUPER_HANDLE%FF%IU4,NSTEP,ML_SUPER_HANDLE%ABN%QHEAT(1:,ML_SUPER_HANDLE%ABN%NCONF))
          ENDIF
! Store the maximum bayesian error and change the criteria if necessary.
          IF (ML_SUPER_HANDLE%FF%ISTART.LT.3 .OR. .NOT.ML_SUPER_HANDLE%FF%LCTIFOR_USE_FROM_FILE) THEN
             CALL CHANGE_CRITERIA (ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD,ML_SUPER_HANDLE%LOGF,ML_SUPER_HANDLE%FF, &
                  ML_SUPER_HANDLE%ABN%BEMAX(1:,ML_SUPER_HANDLE%ABN%NCONF),ML_SUPER_HANDLE%ABN%NIONS)
          ENDIF
! Provide energy, force and stress tensor to VASP
          CALL GIVE_FROM_FF_TO_VASP (ML_SUPER_HANDLE%ABN%MIONS,ML_SUPER_HANDLE%ABN%OMEGA(ML_SUPER_HANDLE%ABN%NCONF), &
               ML_SUPER_HANDLE%ABN%TOTEN_MB(ML_SUPER_HANDLE%ABN%NCONF),ML_SUPER_HANDLE%ABN%TIFOR_MB(1:,1:,ML_SUPER_HANDLE%ABN%NCONF), &
               ML_SUPER_HANDLE%ABN%TSIF_MB(1:,1:,ML_SUPER_HANDLE%ABN%NCONF),TOTEN,TIFOR,TSIF)
! Output end statement of prediction
! Set here the control variable which controls whether VASP_TO_ML is run
! or not in this ionic iteration. This is only done if ML_TO_VASP is 
! called the first time. For the last ionic step, VASP_TO_ML is maybe
! called on data within the block size if present, even if no ab initio 
! step is performed in the last ionic step.
          IF (ML_SUPER_HANDLE%FF%ISTART.NE.2) THEN
             IF (.NOT.ML_SUPER_HANDLE%FF%LSECOND_TIME_ML_TO_VASP) THEN
                IF (ML_SUPER_HANDLE%FF%NSTEP.NE.ML_SUPER_HANDLE%FF%NSW) THEN
                   IF (LDO_AB_INITIO) THEN
                      ML_SUPER_HANDLE%FF%LRUN_VASP_TO_ML=.TRUE.
                   ELSE
                      ML_SUPER_HANDLE%FF%LRUN_VASP_TO_ML=.FALSE.
                   ENDIF
                ELSE
                   IF (LDO_AB_INITIO) THEN
                      ML_SUPER_HANDLE%FF%LRUN_VASP_TO_ML=.TRUE.
                   ELSE      
                      IF (ML_SUPER_HANDLE%ABN%NCONF.GT.1) THEN
                         ML_SUPER_HANDLE%FF%LRUN_VASP_TO_ML=.TRUE.
                      ELSE
                         ML_SUPER_HANDLE%FF%LRUN_VASP_TO_ML=.FALSE.
                      ENDIF
                   ENDIF
                ENDIF
             ENDIF
          ELSE
             ML_SUPER_HANDLE%FF%LRUN_VASP_TO_ML=.FALSE.
          ENDIF
          ! In case the prediction is accurate, i.e. LDO_AB_INITIO is .FALSE.
          ! we can write the status line here. Otherwise, when local
          ! configurations are collected, or actual training is performed, we
          ! have to print out the status line later (see
          ! VASP_TO_ML_MACHINE_LEARNING_LIB) because at this point not all
          ! information is available. Also, in case this is the second call
          ! (LSECOND...) we have to avoid the status output and only print the
          ! finalizing Bayesian error estimation output.
          ! We only write at this point if ISTART.NE.3
          IF ((ML_SUPER_HANDLE%FF%ISTART.NE.3) .AND. &
              .NOT. LDO_AB_INITIO .AND. &
              .NOT. ML_SUPER_HANDLE%FF%LSECOND_TIME_ML_TO_VASP .AND. &
              .NOT. ML_SUPER_HANDLE%FF%LRUN_VASP_TO_ML) THEN
             CALL LOGLINE_STATUS(ML_SUPER_HANDLE%FF,   &
                                 ML_SUPER_HANDLE%LOGF, &
                                 LDO_AB_INITIO)
          END IF
          ! The update of the threshold is print out in the second call
          ! of ML_TO_VASP is IUPDATE_CRITERIA=1 otherwise it's print out
          ! in the first call
          IF (ML_SUPER_HANDLE%FF%IUPDATE_CRITERIA.EQ.1) THEN
             IF (ML_SUPER_HANDLE%FF%LSECOND_TIME_ML_TO_VASP) THEN
                IF (ML_SUPER_HANDLE%FF%LLOG_THRESHOLD_HISTORY) THEN
                   CALL LOGLINE_THRHIST(ML_SUPER_HANDLE%FF, ML_SUPER_HANDLE%LOGF)
                   ! Now reset control variable to .FALSE.
                   ML_SUPER_HANDLE%FF%LLOG_THRESHOLD_HISTORY = .FALSE.
                END IF
                IF (ML_SUPER_HANDLE%FF%LLOG_THRESHOLD_UPDATE) THEN
                   CALL LOGLINE_THRUPD(ML_SUPER_HANDLE%FF, ML_SUPER_HANDLE%LOGF)
                   ! Now reset control variable to .FALSE.
                   ML_SUPER_HANDLE%FF%LLOG_THRESHOLD_UPDATE = .FALSE.
                END IF
                CALL LOGLINE_BEE(ML_SUPER_HANDLE%FF,   &
                                 ML_SUPER_HANDLE%LOGF, &
                                 ML_SUPER_HANDLE%ABN,  &
                                 .TRUE.)
             ENDIF
          ELSE IF (ML_SUPER_HANDLE%FF%IUPDATE_CRITERIA.EQ.0) THEN
             IF (.NOT.LDO_AB_INITIO.AND.ML_SUPER_HANDLE%FF%ISTART.NE.3) THEN
!             IF (.NOT.LDO_AB_INITIO) THEN
                IF (ML_SUPER_HANDLE%FF%LLOG_THRESHOLD_HISTORY) THEN
                   CALL LOGLINE_THRHIST(ML_SUPER_HANDLE%FF, ML_SUPER_HANDLE%LOGF)
                   ! Now reset control variable to .FALSE.
                   ML_SUPER_HANDLE%FF%LLOG_THRESHOLD_HISTORY = .FALSE.
                END IF
                IF (ML_SUPER_HANDLE%FF%LLOG_THRESHOLD_UPDATE) THEN
                   CALL LOGLINE_THRUPD(ML_SUPER_HANDLE%FF, ML_SUPER_HANDLE%LOGF)
                   ! Now reset control variable to .FALSE.
                   ML_SUPER_HANDLE%FF%LLOG_THRESHOLD_UPDATE = .FALSE.
                END IF
             ENDIF
          ENDIF
! The output of BEEF if we have an accurate force field has to 
! be done here after the update of the thresholds, this
! is especially relevant if ML_IUPDATE_CRITERIA=0
          IF (ML_SUPER_HANDLE%FF%ISTART.NE.3 .AND. &
              .NOT. LDO_AB_INITIO .AND. &
              .NOT. ML_SUPER_HANDLE%FF%LSECOND_TIME_ML_TO_VASP .AND. &
              .NOT. ML_SUPER_HANDLE%FF%LRUN_VASP_TO_ML) THEN
             CALL LOGLINE_BEE(ML_SUPER_HANDLE%FF,   &
                              ML_SUPER_HANDLE%LOGF, &
                              ML_SUPER_HANDLE%ABN,  &
                              .TRUE.)
             IF (ML_SUPER_HANDLE%FFM%LCOUPLE) THEN
                CALL LOGLINE_DCOUPLE(ML_SUPER_HANDLE%LOGF,     &
                                     ML_SUPER_HANDLE%FF%NSTEP, &
                                     ML_SUPER_HANDLE%FFM%DHDRCOUPLE)
             END IF
          ENDIF
! This variable is needed later to check whether new data is coming from 
! ab initio calculations or not
          IF (LDO_AB_INITIO) THEN
             ML_SUPER_HANDLE%FF%LABINITIO=.TRUE.
          ELSE
             ML_SUPER_HANDLE%FF%LABINITIO=.FALSE.
          ENDIF
! Set back the helping variable after calling ML_TO_VASP the second time for a given NSTEP 
          IF (ML_SUPER_HANDLE%FF%LSECOND_TIME_ML_TO_VASP) THEN
             ML_SUPER_HANDLE%FF%LSECOND_TIME_ML_TO_VASP=.FALSE.
          ENDIF
        END SUBROUTINE ML_TO_VASP_MACHINE_LEARNING_LIB

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Subroutine to get new ab initio data from VASP and generate new force
! field
!****************************************************************************************************

        SUBROUTINE VASP_TO_ML_MACHINE_LEARNING_LIB (ML_SUPER_HANDLE,A,B,DIR_APP,DIR_LEN,LABORT,NSTEP,OMEGA,POSION,SZNAM2,TOTEN,TIFOR,TSIF)
          USE ABINITIO_DATA
          USE ML_FF_CONSTANT
          USE FORCE_FIELD
          USE IOHANDLE
          USE MPI_DATA
          USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE, ML_IO_WRITE
          USE LOGFILE
          IMPLICIT NONE
          TYPE (ML_SUPER_TYPE) :: ML_SUPER_HANDLE
! Input variables
          REAL(q)          , INTENT(IN)    :: A(:,:) !(1:3,1:3)
          REAL(q)          , INTENT(IN)    :: B(:,:) !(1:3,1:3)
          CHARACTER(LEN=10), INTENT(IN)    :: DIR_APP     ! Location of directory passed from VASP.
          INTEGER          , INTENT(IN)    :: DIR_LEN     ! Length of DIR_APP
          LOGICAL          , INTENT(IN)    :: LABORT
          INTEGER          , INTENT(IN)    :: NSTEP
          REAL(q)          , INTENT(IN)    :: POSION(:,:) !(1:3,1:ML_SUPER_HANDLE%ABN%NIONS)
          REAL(q)          , INTENT(IN)    :: OMEGA
          CHARACTER(LEN=*), INTENT(IN)     :: SZNAM2
! Inout variables
          REAL(q)          , INTENT(INOUT) :: TOTEN
          REAL(q)          , INTENT(INOUT) :: TIFOR(:,:) !(1:3,1:ML_SUPER_HANDLE%ABN%NIONS)
          REAL(q)          , INTENT(INOUT) :: TSIF(:,:) !(1:3,1:3)
! Local variables
          LOGICAL                          :: IERR
          INTEGER                          :: INFO
          INTEGER                          :: NB_ALL(1:ML_SUPER_HANDLE%FF%MTYP)
          REAL(q)                          :: DTIFOR_AV_MAE
          REAL(q)                          :: DTIFOR_AV_RME
          REAL(q)                          :: DTIFOR_MAX
! Setting MD step
          ML_SUPER_HANDLE%FF%NSTEP=NSTEP
! Here we need to check wether it is the last ionic step in the ionic loop.
! Then we also need to check if there are structures within the new candidates to learn
! Only if there are some, then learning is executed.
          IF (ML_SUPER_HANDLE%FF%ISTART.NE.3) THEN
             IF (NSTEP.EQ.ML_SUPER_HANDLE%FF%NSW.AND.ML_SUPER_HANDLE%ABN%NCONF.GT.1) THEN
                ML_SUPER_HANDLE%FF%FORCING_UPDATE=.TRUE.
             ENDIF
          ENDIF
! Initialize DTIFORs here
          DTIFOR_AV_MAE=0.0_q
          DTIFOR_AV_RME=0.0_q
          DTIFOR_MAX=0.0_q
! Put new ab initio data
          CALL SET_ABN_ENERGY (ML_SUPER_HANDLE%FF, &
                               ML_SUPER_HANDLE%TIM,ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD, &
                               ML_SUPER_HANDLE%LOGF,ML_SUPER_HANDLE%ABN,ML_SUPER_HANDLE%FF%IFF,ML_SUPER_HANDLE%FF%MTYP, &
                               ML_SUPER_HANDLE%FFM%NB,NSTEP,TIFOR,TOTEN,TSIF,ML_SUPER_HANDLE%FF%TYPE, &
                               DTIFOR_AV_MAE, DTIFOR_AV_RME, DTIFOR_MAX)
! We need to dial NCONF back one value if no new ab initio data was added
          IF (ML_SUPER_HANDLE%FF%ISTART.NE.3) THEN
             IF (.NOT.ML_SUPER_HANDLE%FF%LABINITIO) THEN
                ML_SUPER_HANDLE%ABN%NCONF=ML_SUPER_HANDLE%ABN%NCONF-1
             ENDIF
          ENDIF
! Add new ab initio data to database
          CALL PUT_ABN_DATA (ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%PAR_SUP_HANDLE, &
                             ML_SUPER_HANDLE%TIM,ML_SUPER_HANDLE%AB,ML_SUPER_HANDLE%ABN,ML_SUPER_HANDLE%FF%CTIFOR,IERR, &
                             ML_SUPER_HANDLE%FF%IFF,ML_SUPER_HANDLE%FF%IFF_OLD,ML_SUPER_HANDLE%FF%ISTART, &
                             LABORT,ML_SUPER_HANDLE%FF%LDOUB,ML_SUPER_HANDLE%FF%LFMAT,ML_SUPER_HANDLE%FF%LGENFF, &
                             ML_SUPER_HANDLE%FF%LLIST,ML_SUPER_HANDLE%FF%LOVER,ML_SUPER_HANDLE%FF%LPS,ML_SUPER_HANDLE%FF%LSAMPLE, &
                             ML_SUPER_HANDLE%FF%LSIG,ML_SUPER_HANDLE%FF%LSOAP,NSTEP,ML_SUPER_HANDLE%FFM%SIG,ML_SUPER_HANDLE%FF%FORCING_UPDATE)
! Here we can assume that LDO_AB_INITIO is .TRUE. and all important
! flags (GENFF, LSAMPLE,...) were already set correctly. Hence, it is
! save to write the status line to the log file.
          CALL LOGLINE_STATUS(ML_SUPER_HANDLE%FF, ML_SUPER_HANDLE%LOGF, .TRUE.)
! If IUPDATE_CRITERIA.EQ.0 then we need to write the threshold here
          IF (ML_SUPER_HANDLE%FF%IUPDATE_CRITERIA.EQ.0) THEN
             IF (ML_SUPER_HANDLE%FF%LLOG_THRESHOLD_HISTORY) THEN
                CALL LOGLINE_THRHIST(ML_SUPER_HANDLE%FF, ML_SUPER_HANDLE%LOGF)
                ML_SUPER_HANDLE%FF%LLOG_THRESHOLD_HISTORY = .FALSE.
             END IF
             IF (ML_SUPER_HANDLE%FF%LLOG_THRESHOLD_UPDATE) THEN
                CALL LOGLINE_THRUPD(ML_SUPER_HANDLE%FF, ML_SUPER_HANDLE%LOGF)
                ML_SUPER_HANDLE%FF%LLOG_THRESHOLD_UPDATE = .FALSE.
             END IF
          ENDIF
! Add new angular descriptor if ML_SUPER_HANDLE%FF%LPS=.TRUE.
          IF(ML_SUPER_HANDLE%FF%LPS) THEN
             CALL PUT_PS (ML_SUPER_HANDLE%TIM, &
                  ML_SUPER_HANDLE%PAR_SUP_HANDLE,ML_SUPER_HANDLE%LOGF,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM, &
                  ML_SUPER_HANDLE%ABN%BE(1:,1:),ML_SUPER_HANDLE%AB%LADD(1:,ML_SUPER_HANDLE%AB%NCONF+1:), &
                  ML_SUPER_HANDLE%AB%LB(1:,ML_SUPER_HANDLE%AB%NCONF+1:),ML_SUPER_HANDLE%AB%LNIONS(1:,1:,ML_SUPER_HANDLE%AB%NCONF+1:), &
                  ML_SUPER_HANDLE%AB%LNTYP(1:,ML_SUPER_HANDLE%AB%NCONF+1:),ML_SUPER_HANDLE%AB%MCONF,ML_SUPER_HANDLE%ABN%MCONF, &
                  ML_SUPER_HANDLE%AB%MITYP,ML_SUPER_HANDLE%AB%MIONS,ML_SUPER_HANDLE%ABN%MIONS,ML_SUPER_HANDLE%AB%MTYP, &
                  ML_SUPER_HANDLE%AB%NCONF,ML_SUPER_HANDLE%ABN%NCONF,ML_SUPER_HANDLE%AB%NITYP(1:,ML_SUPER_HANDLE%AB%NCONF+1:), &
                  ML_SUPER_HANDLE%AB%NIONS(ML_SUPER_HANDLE%AB%NCONF+1),ML_SUPER_HANDLE%ABN%NIONS,ML_SUPER_HANDLE%ABN%NSTEP, &
                  ML_SUPER_HANDLE%AB%NTYP(ML_SUPER_HANDLE%AB%NCONF+1:))
          ENDIF
! Store angular descriptor for basis sets in following conditions.
! 1) ML_SUPER_HANDLE%FF%LPS=.TRUE. In this condition, the machine possesses angular descriptor,
!    and user specifies that the machine needs to store the angular descriptor in all processes.
! 2) ML_SUPER_HANDLE%FF%LMLONLY=.TRUE. (ML_SUPER_HANDLE%FF%ISTART=2). In this condition, the machine also possesses angular descriptor
!    which was read from ML_FF file.
          IF(ML_SUPER_HANDLE%FF%LPS.OR.ML_SUPER_HANDLE%FF%LMLONLY) THEN
             NB_ALL(1:ML_SUPER_HANDLE%FF%MTYP)=ML_SUPER_HANDLE%FFM%NB(1:ML_SUPER_HANDLE%FF%MTYP)+ML_SUPER_HANDLE%FFM%NB_NEW(1:ML_SUPER_HANDLE%FF%MTYP)
             CALL STORE_PS_IN_ALL_PROCESSES (ML_SUPER_HANDLE%PAR_SUP_HANDLE,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM,NB_ALL)
          ENDIF
! Add or calculate list.
          IF(ML_SUPER_HANDLE%FF%LLIST) THEN
             CALL PUT_LIST (ML_SUPER_HANDLE%TIM,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM, &
                            ML_SUPER_HANDLE%AB%LNTYP(1:,ML_SUPER_HANDLE%AB%NCONF+1:),ML_SUPER_HANDLE%AB%MCONF, &
                            ML_SUPER_HANDLE%ABN%MCONF,ML_SUPER_HANDLE%AB%MIONS,ML_SUPER_HANDLE%AB%NCONF, &
                            ML_SUPER_HANDLE%ABN%NCONF,ML_SUPER_HANDLE%AB%NIONS(ML_SUPER_HANDLE%AB%NCONF+1:))
          ENDIF
! Add or calculate SOAP matrix
          IF(ML_SUPER_HANDLE%FF%LSOAP) THEN
             CALL PUT_SOAP(ML_SUPER_HANDLE%TIM,ML_SUPER_HANDLE%PAR_SUP_HANDLE,ML_SUPER_HANDLE%MEM1,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM)
          ENDIF
! Add or calculate design matrix elements
          IF(ML_SUPER_HANDLE%FF%LFMAT) THEN
             CALL PUT_FMAT (ML_SUPER_HANDLE%AB,ML_SUPER_HANDLE%TIM,ML_SUPER_HANDLE%PAR_SUP_HANDLE, &
                  ML_SUPER_HANDLE%MLSHM,ML_SUPER_HANDLE%MEM1,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM, &
                  ML_SUPER_HANDLE%AB%A,ML_SUPER_HANDLE%AB%B,ML_SUPER_HANDLE%AB%LADD,ML_SUPER_HANDLE%AB%LNIONS, &
                  ML_SUPER_HANDLE%AB%LNITYP,ML_SUPER_HANDLE%AB%LNTYP,ML_SUPER_HANDLE%AB%MCONF,ML_SUPER_HANDLE%ABN%MCONF, &
                  ML_SUPER_HANDLE%AB%MIONS,ML_SUPER_HANDLE%AB%MITYP,ML_SUPER_HANDLE%AB%MTYP,ML_SUPER_HANDLE%AB%NCONF, &
                  ML_SUPER_HANDLE%ABN%NCONF,ML_SUPER_HANDLE%AB%NIONS,ML_SUPER_HANDLE%AB%NITYP,ML_SUPER_HANDLE%AB%NTYP, &
                  ML_SUPER_HANDLE%AB%OMEGA,ML_SUPER_HANDLE%AB%POSION)
          ENDIF
! Renew the size data if ML_SUPER_HANDLE%FF%LGENFF=.TRUE.
! In this condition, the machine tries to generate new force field.
! Sizes of array, therefore, must be changed.
          IF(ML_SUPER_HANDLE%FF%LGENFF) THEN
             CALL PUT_SIZE (ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM, &
                  ML_SUPER_HANDLE%ABN%ISYS,LABORT,ML_SUPER_HANDLE%ABN%MCONF,ML_SUPER_HANDLE%AB%NCONF, &
                  ML_SUPER_HANDLE%ABN%NCONF,ML_SUPER_HANDLE%AB%NCONF_SYS,NSTEP,ML_SUPER_HANDLE%AB%NSYS)
          ENDIF
! Calculate variance of the ab initio data
          IF(ML_SUPER_HANDLE%FF%LSIG) THEN
             CALL CAL_SIG (ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%AB%ISYS,ML_SUPER_HANDLE%AB%LADD,ML_SUPER_HANDLE%AB%MCONF, &
                  ML_SUPER_HANDLE%AB%MIONS,ML_SUPER_HANDLE%AB%MTYP,ML_SUPER_HANDLE%AB%NCONF, &
                  ML_SUPER_HANDLE%AB%NCONF_SYS,ML_SUPER_HANDLE%AB%NIONS,ML_SUPER_HANDLE%AB%NITYP, &
                  ML_SUPER_HANDLE%AB%NSYS,ML_SUPER_HANDLE%AB%NTYP,ML_SUPER_HANDLE%AB%TIFOR, &
                  ML_SUPER_HANDLE%AB%TOTEN,ML_SUPER_HANDLE%AB%TSIF)
          ENDIF
! Sparsification of stored basis set data if available
          IF(ML_SUPER_HANDLE%FF%LOVER) THEN
             CALL SPARS_ICUR (ML_SUPER_HANDLE%FF%ML_EPS_LOW,ML_SUPER_HANDLE%TIM, &
                              ML_SUPER_HANDLE%LAPACKWORK,ML_SUPER_HANDLE%PAR_SUP_HANDLE, &
                              ML_SUPER_HANDLE%MEM1,ML_SUPER_HANDLE%LOGF, &
                              ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM, &
                              ML_SUPER_HANDLE%AB%A,ML_SUPER_HANDLE%AB%B,ML_SUPER_HANDLE%AB%ISYS, &
                              ML_SUPER_HANDLE%AB%LADD,ML_SUPER_HANDLE%AB%LB,ML_SUPER_HANDLE%AB%LNITYP, &
                              ML_SUPER_HANDLE%AB%LNIONS,ML_SUPER_HANDLE%AB%LNTYP,ML_SUPER_HANDLE%AB%LSYS, &
                              ML_SUPER_HANDLE%AB%MCONF,ML_SUPER_HANDLE%AB%MITYP,ML_SUPER_HANDLE%AB%MIONS, &
                              ML_SUPER_HANDLE%AB%MTYP,ML_SUPER_HANDLE%AB%NCONF,ML_SUPER_HANDLE%AB%NCONF_SYS, &
                              ML_SUPER_HANDLE%AB%NITYP,ML_SUPER_HANDLE%AB%NIONS,NSTEP,ML_SUPER_HANDLE%AB%NSYS, &
                              ML_SUPER_HANDLE%AB%NTYP,ML_SUPER_HANDLE%AB%POSION, &
                              ML_SUPER_HANDLE%AB%TIFOR,ML_SUPER_HANDLE%AB%TOTEN,ML_SUPER_HANDLE%AB%TSIF, &
                              ML_SUPER_HANDLE%AB%TYPE,ML_SUPER_HANDLE%AB%OMEGA,ML_SUPER_HANDLE%AB%SZNAM2)
          ENDIF
! For upcoming MD step, this subroutine needs to be called again here in order to store all the restructured C00 and PS arrays.
          IF(ML_SUPER_HANDLE%FF%LPS.OR.ML_SUPER_HANDLE%FF%LMLONLY) THEN
             CALL STORE_PS_IN_ALL_PROCESSES(ML_SUPER_HANDLE%PAR_SUP_HANDLE,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM,ML_SUPER_HANDLE%FFM%NB)
          ENDIF
! Calculate variance of the ab initio data
          IF(ML_SUPER_HANDLE%FF%LSIG) THEN 
             CALL CAL_SIG (ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%AB%ISYS,ML_SUPER_HANDLE%AB%LADD,ML_SUPER_HANDLE%AB%MCONF, &
                           ML_SUPER_HANDLE%AB%MIONS,ML_SUPER_HANDLE%AB%MTYP,ML_SUPER_HANDLE%AB%NCONF, &
                           ML_SUPER_HANDLE%AB%NCONF_SYS,ML_SUPER_HANDLE%AB%NIONS,ML_SUPER_HANDLE%AB%NITYP, &
                           ML_SUPER_HANDLE%AB%NSYS,ML_SUPER_HANDLE%AB%NTYP,ML_SUPER_HANDLE%AB%TIFOR, &
                           ML_SUPER_HANDLE%AB%TOTEN,ML_SUPER_HANDLE%AB%TSIF)
          ENDIF
! Generation of force field from the stored ab initio data if available
          IF(ML_SUPER_HANDLE%FF%LGENFF) THEN
             CALL GEN_FF (ML_SUPER_HANDLE%AB,ML_SUPER_HANDLE%TIM, &
                          ML_SUPER_HANDLE%LAPACKWORK,ML_SUPER_HANDLE%PAR_SUP_HANDLE, &
                          ML_SUPER_HANDLE%MEM1,ML_SUPER_HANDLE%LOGF, &
                          ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM, &
                          ML_SUPER_HANDLE%AB%A,ML_SUPER_HANDLE%AB%B,ML_SUPER_HANDLE%AB%ISYS,ML_SUPER_HANDLE%AB%LADD, &
                          ML_SUPER_HANDLE%AB%LNITYP,ML_SUPER_HANDLE%AB%LNIONS,ML_SUPER_HANDLE%AB%LNTYP, &
                          ML_SUPER_HANDLE%AB%MCONF,ML_SUPER_HANDLE%AB%MITYP,ML_SUPER_HANDLE%AB%MIONS,ML_SUPER_HANDLE%AB%MTYP, &
                          ML_SUPER_HANDLE%AB%NCONF,ML_SUPER_HANDLE%AB%NCONF_SYS,ML_SUPER_HANDLE%AB%NITYP, &
                          ML_SUPER_HANDLE%AB%NIONS,NSTEP,ML_SUPER_HANDLE%AB%NSYS,ML_SUPER_HANDLE%AB%NTYP, &
                          ML_SUPER_HANDLE%AB%POSION,ML_SUPER_HANDLE%AB%OMEGA)
          ENDIF
! Do the rest of writing here
! Do first writing of BEE only if LGENFF=.TRUE. and IFF=.TRUE., otherwise write BEEF
          IF(ML_SUPER_HANDLE%FF%LGENFF) THEN
             CALL LOGLINE_STDAB_ERR(ML_SUPER_HANDLE%FF, ML_SUPER_HANDLE%LOGF)
             CALL LOGLINE_CFE(ML_SUPER_HANDLE%LOGF,     &
                              ML_SUPER_HANDLE%FF%NSTEP, &
                              DTIFOR_AV_RME,            &
                              DTIFOR_AV_MAE,            &
                              DTIFOR_MAX)
             CALL LOGLINE_LAST_ERR(ML_SUPER_HANDLE%FF, &
                                   ML_SUPER_HANDLE%LOGF)
          ENDIF
          IF (ML_SUPER_HANDLE%FFM%LCOUPLE) THEN
             CALL LOGLINE_DCOUPLE(ML_SUPER_HANDLE%LOGF,     &
                                  ML_SUPER_HANDLE%FF%NSTEP, &
                                  ML_SUPER_HANDLE%FFM%DHDRCOUPLE)
          END IF
          IF (ML_SUPER_HANDLE%FF%LGENFF) THEN
             IF (ML_SUPER_HANDLE%FF%IUPDATE_CRITERIA.EQ.0) THEN
                CALL LOGLINE_BEE(ML_SUPER_HANDLE%FF,   &
                                 ML_SUPER_HANDLE%LOGF, &
                                 ML_SUPER_HANDLE%ABN,  &
                                 .TRUE.)
             ELSE
                CALL LOGLINE_BEE(ML_SUPER_HANDLE%FF,   &
                                 ML_SUPER_HANDLE%LOGF, &
                                 ML_SUPER_HANDLE%ABN,  &
                                 .FALSE.)
             ENDIF
          ELSE
             CALL LOGLINE_BEE(ML_SUPER_HANDLE%FF,   &
                              ML_SUPER_HANDLE%LOGF, &
                              ML_SUPER_HANDLE%ABN,  &
                              .TRUE.)
          ENDIF
! If ML_SUPER_HANDLE%FF%IFF=.TRUE., make all processes own the regressin coefficients, inverse SOAP matrices and covariance matrix.
          IF(ML_SUPER_HANDLE%FF%IFF) THEN
             CALL STORE_WMAT_FFM_IN_ALL_PROCESS(ML_SUPER_HANDLE%PAR_SUP_HANDLE,ML_SUPER_HANDLE%MEM1,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM)
             CALL STORE_CMAT_FFM_IN_ALL_PROCESS(ML_SUPER_HANDLE%PAR_SUP_HANDLE,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM)
          ENDIF
! If ML_SUPER_HANDLE%FF%IFF=.TRUE., make all processes own the regressin coefficients, inverse SOAP matrices and covariance matrix.
! Calculate histogram if ML_SUPER_HANDLE%FF%LGENFF=.TRUE. and ML_SUPER_HANDLE%FF%IFF=.TRUE.
! In this condition, the machine tried to generate the force field, and the generation was successfully finished.
! Therefore, the program can evaluate the distributions of training data and regression results.
          IF(ML_SUPER_HANDLE%FF%LGENFF.AND.ML_SUPER_HANDLE%FF%IFF) THEN
             CALL HISTOGRAM (ML_SUPER_HANDLE%PAR_SUP_HANDLE,ML_SUPER_HANDLE%AB,DIR_APP,DIR_LEN)
          ENDIF
! Print out ML_FFN if counter reached block size 
          ! For ML_ISTART=2 we won't write out
          IF (ML_SUPER_HANDLE%FF%ISTART.NE.2) THEN
             IF (ML_SUPER_HANDLE%FF%LGENFF.AND.ML_SUPER_HANDLE%FF%IFF) THEN
                ML_SUPER_HANDLE%FF%NCOUNT_GENFF=ML_SUPER_HANDLE%FF%NCOUNT_GENFF+1
                IF (ML_SUPER_HANDLE%FF%NCOUNT_GENFF.GE.ML_SUPER_HANDLE%FF%NBLOCK_FFN) THEN
                   CALL OUT_FFNCAR (ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM,ML_SUPER_HANDLE%PAR_SUP_HANDLE,DIR_APP,DIR_LEN)
                   ML_SUPER_HANDLE%FF%NCOUNT_GENFF=0
                ENDIF
             ENDIF
          ENDIF
! Output results of regression if ML_SUPER_HANDLE%FF%IFF=.TRUE. and ML_SUPER_HANDLE%FF%LGENFF=.TRUE.
          IF(ML_SUPER_HANDLE%FF%LGENFF.AND.ML_SUPER_HANDLE%FF%IFF) THEN
             CALL OUTPUT_REG (ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD, &
                  DIR_APP,DIR_LEN,ML_SUPER_HANDLE%AB%IU3,ML_SUPER_HANDLE%AB%MCONF,ML_SUPER_HANDLE%AB%MIONS, &
                  ML_SUPER_HANDLE%AB%NCONF,ML_SUPER_HANDLE%AB%NIONS,0,ML_SUPER_HANDLE%AB%REGCAR, &
                  ML_SUPER_HANDLE%AB%TIFOR,ML_SUPER_HANDLE%AB%TIFOR_MB,ML_SUPER_HANDLE%AB%TOTEN,ML_SUPER_HANDLE%AB%TOTEN_MB, &
                  ML_SUPER_HANDLE%AB%TSIF,ML_SUPER_HANDLE%AB%TSIF_MB)
          ENDIF
! Output the ab initio data if ML_SUPER_HANDLE%FF%LGENFF=.TRUE. and ML_SUPER_HANDLE%FF%IFF=.TRUE.
          IF(ML_SUPER_HANDLE%FF%LGENFF.AND.ML_SUPER_HANDLE%FF%IFF) THEN
             CALL OUTPUT_ABNCAR (ML_SUPER_HANDLE%AB,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD, &
                  ML_SUPER_HANDLE%AB%A,ML_SUPER_HANDLE%AB%ABNCAR,DIR_APP,DIR_LEN,ML_SUPER_HANDLE%FF%EATOM_REF,ML_SUPER_HANDLE%AB%IU2, &
                  ML_SUPER_HANDLE%AB%LADD,ML_SUPER_HANDLE%FFM%LNCONF,ML_SUPER_HANDLE%FFM%LNIONS, &
                  ML_SUPER_HANDLE%FFM%MB,ML_SUPER_HANDLE%AB%MCONF,ML_SUPER_HANDLE%AB%MITYP, &
                  ML_SUPER_HANDLE%AB%MIONS,ML_SUPER_HANDLE%AB%MTYP,ML_SUPER_HANDLE%FFM%NB,ML_SUPER_HANDLE%AB%NCONF,ML_SUPER_HANDLE%AB%NITYP, &
                  ML_SUPER_HANDLE%AB%NIONS,ML_SUPER_HANDLE%AB%NTYP,ML_SUPER_HANDLE%AB%POMASS,ML_SUPER_HANDLE%AB%POSION, &
                  ML_SUPER_HANDLE%AB%SZNAM2,ML_SUPER_HANDLE%AB%TIFOR,ML_SUPER_HANDLE%AB%TOTEN,ML_SUPER_HANDLE%AB%TSIF, &
                  ML_SUPER_HANDLE%AB%TYPE)
          ENDIF
! If LTEST=.TRUE. and. LTRJ=.FALSE., and if the force field is accurate
! enough, machine-learning data should be restored to obtain
! machine-learning MD trajectory
          IF((ML_SUPER_HANDLE%FF%IFF).AND.(ML_SUPER_HANDLE%FF%LTEST).AND.(.NOT.ML_SUPER_HANDLE%FF%LTRJ).AND.(.NOT.IERR)) THEN
             CALL RESTORE_ML_DATA (ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD,ML_SUPER_HANDLE%LOGF,ML_SUPER_HANDLE%ABN,NSTEP)
             CALL GIVE_FROM_FF_TO_VASP (ML_SUPER_HANDLE%ABN%MIONS,ML_SUPER_HANDLE%ABN%OMEGA(ML_SUPER_HANDLE%ABN%NCONF), &
                  ML_SUPER_HANDLE%ABN%TOTEN_MB(ML_SUPER_HANDLE%ABN%NCONF),ML_SUPER_HANDLE%ABN%TIFOR_MB(1:,1:,ML_SUPER_HANDLE%ABN%NCONF), &
                  ML_SUPER_HANDLE%ABN%TSIF_MB(1:,1:,ML_SUPER_HANDLE%ABN%NCONF),TOTEN,TIFOR,TSIF)
          ENDIF
!TEST 
          IF (ML_SUPER_HANDLE%FF%LGENFF) THEN
             ML_SUPER_HANDLE%FFM%FMAT_NEW=0.0_q
          ENDIF
!TEST
! Reset ABN%NCONF as one.
          CALL RENEW_SAMPLING_PAR (ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%TIM, &
                                   ML_SUPER_HANDLE%ABN,IERR,LABORT,NSTEP) 
! This flag needs to be set to now whether ML_TO_VASP is called a second time for the same NSTEP
          IF (ML_SUPER_HANDLE%FF%LGENFF .AND. ML_SUPER_HANDLE%FF%IUPDATE_CRITERIA.EQ.1) THEN
             ML_SUPER_HANDLE%FF%LSECOND_TIME_ML_TO_VASP = .TRUE.
          ENDIF
        END SUBROUTINE VASP_TO_ML_MACHINE_LEARNING_LIB
#endif

!****************************************************************************************************
! Subroutine to finalize the machine-learning calculations
!****************************************************************************************************

        SUBROUTINE MACHINE_LEARNING_FINISH_LIB (ML_SUPER_HANDLE,DIR_APP,DIR_LEN,NSTEP)
          USE ABINITIO_DATA
          USE LOGFILE
          USE IOHANDLE
          USE FORCE_FIELD
          USE MEMORY_USAGE
          USE MPI_DATA
          USE TIME_CHECK
          USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE
          IMPLICIT NONE
          TYPE (ML_SUPER_TYPE) :: ML_SUPER_HANDLE
! Input variable
          CHARACTER(LEN=10), INTENT(IN) :: DIR_APP
          INTEGER          , INTENT(IN) :: DIR_LEN
          INTEGER          , INTENT(IN) :: NSTEP
#ifndef ML_LEARNING_NOT_ACTIVE
          IF(ML_SUPER_HANDLE%FF%LGENFF.AND.ML_SUPER_HANDLE%FF%IFF) THEN
             CALL OUTPUT_ABNCAR (ML_SUPER_HANDLE%AB,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD, &
                  ML_SUPER_HANDLE%AB%A,ML_SUPER_HANDLE%AB%ABNCAR,DIR_APP,DIR_LEN,ML_SUPER_HANDLE%FF%EATOM_REF, &
                  ML_SUPER_HANDLE%AB%IU2,ML_SUPER_HANDLE%AB%LADD,ML_SUPER_HANDLE%FFM%LNCONF,ML_SUPER_HANDLE%FFM%LNIONS, &
                  ML_SUPER_HANDLE%FFM%MB,ML_SUPER_HANDLE%AB%MCONF,ML_SUPER_HANDLE%AB%MITYP,ML_SUPER_HANDLE%AB%MIONS, &
                  ML_SUPER_HANDLE%AB%MTYP,ML_SUPER_HANDLE%FFM%NB,ML_SUPER_HANDLE%AB%NCONF, &
                  ML_SUPER_HANDLE%AB%NITYP,ML_SUPER_HANDLE%AB%NIONS,ML_SUPER_HANDLE%AB%NTYP,ML_SUPER_HANDLE%AB%POMASS, &
                  ML_SUPER_HANDLE%AB%POSION,ML_SUPER_HANDLE%AB%SZNAM2,ML_SUPER_HANDLE%AB%TIFOR,ML_SUPER_HANDLE%AB%TOTEN, &
                  ML_SUPER_HANDLE%AB%TSIF,ML_SUPER_HANDLE%AB%TYPE)
          ENDIF
          ! For ML_ISTART=2 we won't write out
          IF (ML_SUPER_HANDLE%FF%ISTART.NE.2) THEN
             IF(ML_SUPER_HANDLE%FF%IFF) THEN
                CALL OUT_FFNCAR (ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM,ML_SUPER_HANDLE%PAR_SUP_HANDLE,DIR_APP,DIR_LEN)
             ENDIF
          ENDIF
#endif
          CALL LOGSECTION_MEMORY(ML_SUPER_HANDLE%LOGF, 1, ML_SUPER_HANDLE%MEM1)
          CALL LOGSECTION_TIME(ML_SUPER_HANDLE%LOGF,       &
                               ML_SUPER_HANDLE%TIM%CPUTIM, &
                               ML_SUPER_HANDLE%TIM%VPUTIM)
                               
          CALL DEALLOCATE_FFARRAY(ML_SUPER_HANDLE%MLSHM,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM)
#ifndef ML_LEARNING_NOT_ACTIVE
          CALL DEALLOCATE_AB_ARRAY(ML_SUPER_HANDLE%AB)
#endif
          CALL DEALLOCATE_ABN_ARRAY(ML_SUPER_HANDLE%ABN)
! Close log file
          CLOSE(ML_SUPER_HANDLE%LOGF%IU)
          IF(ML_SUPER_HANDLE%FFM%LEATOM) THEN
             CLOSE(UNIT=ML_SUPER_HANDLE%FF%IU3)
          ENDIF
          IF(ML_SUPER_HANDLE%FFM%LHEAT) THEN
             CLOSE(UNIT=ML_SUPER_HANDLE%FF%IU4)
          ENDIF
        END SUBROUTINE MACHINE_LEARNING_FINISH_LIB

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Subroutine for construction of new ML_ABN file (with different basis sets) only from ML_AB file
!****************************************************************************************************
        SUBROUTINE ISTART_3_DRIVER(ML_SUPER_HANDLE,MLH_TMP)
          USE ML_FF_PREC
          USE ML_FF_CONSTANT
          USE LOGFILE
          USE IOHANDLE
          USE MPI_DATA
          USE TIME_CHECK
          USE ABINITIO_DATA
          USE FORCE_FIELD
          USE MEMORY_USAGE
          USE RANDOM_NUMBER_GENERATOR
          USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE, ML_FILESTREAM_GLOBAL, ML_IO_WRITE
          IMPLICIT NONE
! Input variables
          TYPE (ML_SUPER_TYPE) :: ML_SUPER_HANDLE
          TYPE (ML_SUPER_TYPE) :: MLH_TMP
! Local variables
          INTEGER              :: ICONF
          INTEGER              :: ICONF_SYS
          INTEGER              :: ISYS
          INTEGER              :: ITYPE
          LOGICAL              :: LDO_AB_INITIO
          INTEGER              :: NSTEP
          LOGICAL              :: REINITIALIZE
          REAL(q), ALLOCATABLE :: TIFOR_TMP(:,:)
          REAL(q)              :: TOTEN_TMP
          REAL(q)              :: TSIF_TMP(1:3,1:3)
          REAL(q), ALLOCATABLE :: VEL_TMP(:,:)
! Set variables to skip machine learning routines later
          ML_SUPER_HANDLE%SKIP_ML_TO_VASP = .TRUE.
          ML_SUPER_HANDLE%SKIP_VASP_TO_ML = .TRUE.
! Initialize NSTEP
          NSTEP=1
! We need to set the filestream for the helping types
          MLH_TMP%AB%ABCAR        = ML_SUPER_HANDLE%AB%ABCAR
          MLH_TMP%AB%IU1          = ML_SUPER_HANDLE%AB%IU1
          MLH_TMP%FF%FFCAR        = ML_SUPER_HANDLE%FF%FFCAR
          MLH_TMP%FF%IU1          = ML_SUPER_HANDLE%FF%IU1
          MLH_TMP%AB%ABNCAR       = ML_SUPER_HANDLE%AB%ABNCAR
          MLH_TMP%AB%IU2          = ML_SUPER_HANDLE%AB%IU2
          MLH_TMP%FF%FFNCAR       = ML_SUPER_HANDLE%FF%FFNCAR
          MLH_TMP%FF%IU2          = ML_SUPER_HANDLE%FF%IU2
          MLH_TMP%LOGF%ML_LOGFILE = ML_SUPER_HANDLE%LOGF%ML_LOGFILE
          MLH_TMP%LOGF%IU         = ML_SUPER_HANDLE%LOGF%IU
! Do the inititialization, partially inititialized with first structure
! but initialization of the other parts of first structure need to be
! also called later (INIT_FF_ISTART_3_SUBSEQUENT_CALLS)
          IF(ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD%NODE_ME.EQ.0) THEN
             WRITE(*,*) "Starting to select new local configurations from ML_AB file (ML_FF_ISTART=3):"
          ENDIF
          CALL INIT_FF_ISTART_3_FIRST_CALL(ML_SUPER_HANDLE,MLH_TMP) 
! Do first output of loop header here
          CALL LOGSECTION_LOOP_HEADER(ML_SUPER_HANDLE%FF, ML_SUPER_HANDLE%LOGF)
          ICONF=0
! Loop over training structures
          DO ISYS=1, MLH_TMP%AB%NSYS
             DO ICONF_SYS=1,MLH_TMP%AB%NCONF_SYS(ISYS)
                ICONF=ICONF+1
! We need to reinitialize if the number of elements or atoms change in the loop of the structures
                IF (ICONF_SYS.EQ.1) THEN
                   REINITIALIZE=.TRUE.
                ELSE
                   REINITIALIZE=.FALSE.
                   IF (MLH_TMP%AB%SZNAM2(ICONF).NE.MLH_TMP%AB%SZNAM2(ICONF-1)) REINITIALIZE=.TRUE.
                   IF (MLH_TMP%AB%NIONS(ICONF).NE.MLH_TMP%AB%NIONS(ICONF-1)) REINITIALIZE=.TRUE.
                   IF (MLH_TMP%AB%NTYP(ICONF).NE.MLH_TMP%AB%NTYP(ICONF-1)) THEN
                      REINITIALIZE=.TRUE.
                   ELSE
                      DO ITYPE=1,MLH_TMP%AB%NTYP(ICONF)
                        IF (MLH_TMP%AB%NITYP(ITYPE,ICONF).NE.MLH_TMP%AB%NITYP(ITYPE,ICONF-1)) REINITIALIZE=.TRUE.
                     ENDDO
                   ENDIF
                ENDIF
! We need to reinitialize the structure if the number of elements or atoms change
                IF (REINITIALIZE) THEN
                   CALL INIT_FF_ISTART_3_SUBSEQUENT_CALLS(ML_SUPER_HANDLE,MLH_TMP,ICONF)
! Since we don't have a force field after the change of atoms in the ab initio data
                   IF (ALLOCATED(TIFOR_TMP)) DEALLOCATE(TIFOR_TMP)
                   IF (ALLOCATED(VEL_TMP)) DEALLOCATE(VEL_TMP)
                   ALLOCATE(TIFOR_TMP(1:3,1:ML_SUPER_HANDLE%ABN%NIONS))
                   ALLOCATE(VEL_TMP(1:3,1:ML_SUPER_HANDLE%ABN%NIONS))
                ENDIF
! Force learning and creation of force field if it's the last structure for a given system
                IF (ICONF_SYS.EQ.MLH_TMP%AB%NCONF_SYS(ISYS)) THEN
                   ML_SUPER_HANDLE%FF%FORCING_UPDATE=.TRUE.
                ENDIF
! Forcing update if maximum number of new configurations is reached
                IF (ML_SUPER_HANDLE%ABN%NCONF.EQ.ML_SUPER_HANDLE%ABN%MCONF) THEN
                   ML_SUPER_HANDLE%FF%FORCING_UPDATE=.TRUE.
                ENDIF
! Set here CTIFOR if CTIFOR was supplied on the ML_AB file and wants to be used
                IF (ML_SUPER_HANDLE%FF%LCTIFOR_USE_FROM_FILE) THEN
                   ML_SUPER_HANDLE%FF%CTIFOR=ML_SUPER_HANDLE%AB%CTIFOR_ALLCONF_OLD(ICONF)
                ENDIF 
! Call ML_TO_VASP and VASP_TO_ML for the rest of the structures
                CALL ML_TO_VASP_MACHINE_LEARNING_LIB(ML_SUPER_HANDLE, MLH_TMP%AB%A(1:,1:,ICONF), MLH_TMP%AB%B(1:,1:,ICONF), &
                     ML_SUPER_HANDLE%ML_INPUT_HANDLE%POTIM, NSTEP, MLH_TMP%AB%OMEGA(ICONF), MLH_TMP%AB%POSION(1:,1:,ICONF), VEL_TMP, &
                     MLH_TMP%AB%SZNAM2(ICONF), TOTEN_TMP, TIFOR_TMP, TSIF_TMP, LDO_AB_INITIO)
                CALL VASP_TO_ML_MACHINE_LEARNING_LIB(ML_SUPER_HANDLE, MLH_TMP%AB%A(1:,1:,ICONF), MLH_TMP%AB%B(1:,1:,ICONF), &
                     ML_SUPER_HANDLE%ML_INPUT_HANDLE%DIR_APP, ML_SUPER_HANDLE%ML_INPUT_HANDLE%DIR_LEN,.FALSE.,NSTEP, &
                     MLH_TMP%AB%OMEGA(ICONF), MLH_TMP%AB%POSION(1:,1:,ICONF), MLH_TMP%AB%SZNAM2(ICONF), MLH_TMP%AB%TOTEN(ICONF), &
                     MLH_TMP%AB%TIFOR(1:,1:,ICONF), MLH_TMP%AB%TSIF(1:,1:,ICONF))
                IF (ML_IO_WRITE) THEN
                     WRITE(*,*) "Conf. ",ICONF,"done from",MLH_TMP%AB%NCONF,"configurations| FF_GEN:",ML_SUPER_HANDLE%FF%LGENFF
                ENDIF
                IF (ML_SUPER_HANDLE%FF%LGENFF.AND.ML_SUPER_HANDLE%FF%IUPDATE_CRITERIA.EQ.1) THEN
                   CALL ML_TO_VASP_MACHINE_LEARNING_LIB(ML_SUPER_HANDLE, MLH_TMP%AB%A(1:,1:,ICONF), MLH_TMP%AB%B(1:,1:,ICONF), &
                        ML_SUPER_HANDLE%ML_INPUT_HANDLE%POTIM, NSTEP, MLH_TMP%AB%OMEGA(ICONF), MLH_TMP%AB%POSION(1:,1:,ICONF), VEL_TMP, &
                        MLH_TMP%AB%SZNAM2(ICONF), TOTEN_TMP, TIFOR_TMP, TSIF_TMP, LDO_AB_INITIO)
                ENDIF
! Increment the step of the calculation
                NSTEP=NSTEP+1
             ENDDO
          ENDDO
! Deallocate some of the helping arrays
          CALL DEALLOCATE_AB_ARRAY(MLH_TMP%AB)
          IF (ALLOCATED(ML_SUPER_HANDLE%AB%MIONS_CONF)) DEALLOCATE(ML_SUPER_HANDLE%AB%MIONS_CONF)
          IF (ALLOCATED(ML_SUPER_HANDLE%AB%MTYP_CONF)) DEALLOCATE(ML_SUPER_HANDLE%AB%MTYP_CONF)
          IF (ALLOCATED(ML_SUPER_HANDLE%AB%MITYP_CONF)) DEALLOCATE(ML_SUPER_HANDLE%AB%MITYP_CONF)
          IF (ALLOCATED(TIFOR_TMP)) DEALLOCATE(TIFOR_TMP)
          IF (ALLOCATED(VEL_TMP)) DEALLOCATE(VEL_TMP)
        END SUBROUTINE ISTART_3_DRIVER

!****************************************************************************************************
! First initialization of run needed from ISTART_3_DRIVER
!****************************************************************************************************
        SUBROUTINE INIT_FF_ISTART_3_FIRST_CALL(ML_SUPER_HANDLE,MLH_TMP)
           USE ML_FF_PREC
           USE ML_FF_CONSTANT
           USE LOGFILE
           USE IOHANDLE
           USE MPI_DATA
           USE TIME_CHECK
           USE ABINITIO_DATA
           USE FORCE_FIELD
           USE MEMORY_USAGE
           USE RANDOM_NUMBER_GENERATOR
           USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE, ML_FILESTREAM_GLOBAL, ML_IO_WRITE
           IMPLICIT NONE
! Input and output variables
           TYPE (ML_SUPER_TYPE) :: ML_SUPER_HANDLE
           TYPE (ML_SUPER_TYPE) :: MLH_TMP
! Local variables
           INTEGER              :: IC
           INTEGER              :: IHELP
           INTEGER              :: ITYP
           INTEGER              :: JTYP
           CHARACTER (LEN=2), ALLOCATABLE :: TYPE_HELP(:)
           INTEGER :: ICONF 
!-----------------------------------------------------------------
! READ_ABCAR_RESET_SIZES
!-----------------------------------------------------------------
! Check size of FF%MCONF_ALLOCATE, FF%MB... from ML_AB file if ML_FF_ISTART=1,3
          CALL READ_ABCAR_RESET_SIZES(ML_SUPER_HANDLE%AB,ML_SUPER_HANDLE%FF, &
                                      ML_SUPER_HANDLE%FFM,ML_SUPER_HANDLE%ML_INPUT_HANDLE, &
                                      ML_SUPER_HANDLE%TAG_LIST)
!-----------------------------------------------------------------
! DEF_ELEMENT_ORDER_FF
!-----------------------------------------------------------------
! Define element order
          CALL DEF_ELEMENT_ORDER_FF (ML_SUPER_HANDLE%TIM, &
                                     ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM,ML_SUPER_HANDLE%FF%FFCAR, &
                                     ML_SUPER_HANDLE%AB%ABCAR, &
                                     ML_SUPER_HANDLE%ML_INPUT_HANDLE%EATOM_VASP, &
                                     ML_SUPER_HANDLE%AB%IU1,ML_SUPER_HANDLE%FF%IU1, &
                                     ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP, &
                                     ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP, &
                                     ML_SUPER_HANDLE%ML_INPUT_HANDLE%POMASS, &
                                     ML_SUPER_HANDLE%ML_INPUT_HANDLE%TYPE,.TRUE.,3)
!-----------------------------------------------------------------
! DEF_SYSTEM_ORDER_FF
!-----------------------------------------------------------------
! Define system order.
          CALL DEF_SYSTEM_ORDER_FF (ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM,ML_SUPER_HANDLE%AB%ABCAR, &
                                    ML_SUPER_HANDLE%AB%IU1,ML_SUPER_HANDLE%ML_INPUT_HANDLE%SZNAM2,.TRUE.,3)
!-----------------------------------------------------------------
! ESTIMATE_MEMORY
!-----------------------------------------------------------------
! Estimate necessary memory
          CALL ESTIMATE_MEMORY (ML_SUPER_HANDLE%TIM, &
               ML_SUPER_HANDLE%PAR_SUP_HANDLE,ML_SUPER_HANDLE%MEM0,ML_SUPER_HANDLE%MEM1, &
               ML_SUPER_HANDLE%LOGF,ML_SUPER_HANDLE%ML_INPUT_HANDLE%A,ML_SUPER_HANDLE%ML_INPUT_HANDLE%B, &
               ML_SUPER_HANDLE%FFM%IBROAD1,ML_SUPER_HANDLE%FFM%IBROAD2,ML_SUPER_HANDLE%FF%ISTART, &
               ML_SUPER_HANDLE%FFM%LCOUPLE,ML_SUPER_HANDLE%FFM%LMAX1,ML_SUPER_HANDLE%FFM%LMAX2, &
               ML_SUPER_HANDLE%FF%LMLMB,ML_SUPER_HANDLE%FF%LMLONLY,ML_SUPER_HANDLE%FFM%LSIC,ML_SUPER_HANDLE%FFM%LSUPERVEC, &
               ML_SUPER_HANDLE%FFM%MB,ML_SUPER_HANDLE%FF%MCONF,ML_SUPER_HANDLE%FF%MCONF_NEW,ML_SUPER_HANDLE%FF%MIONS, &
               ML_SUPER_HANDLE%FF%MITYP,ML_SUPER_HANDLE%FFM%MSPL1,ML_SUPER_HANDLE%FFM%MSPL2,ML_SUPER_HANDLE%FF%MTYP, &
               ML_SUPER_HANDLE%FFM%NHYP1,ML_SUPER_HANDLE%FFM%NHYP2,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP, &
               ML_SUPER_HANDLE%FFM%NR1,ML_SUPER_HANDLE%FFM%NR2,ML_SUPER_HANDLE%FFM%MRB1,ML_SUPER_HANDLE%FFM%MRB2, &
               ML_SUPER_HANDLE%FF%NSW,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP, &
               ML_SUPER_HANDLE%ML_INPUT_HANDLE%OMEGA,ML_SUPER_HANDLE%FFM%RCUT1,ML_SUPER_HANDLE%FFM%RCUT2, &
               ML_SUPER_HANDLE%FFM%W1,ML_SUPER_HANDLE%FFM%W2)
!-----------------------------------------------------------------
! SET_FF_PAR
!-----------------------------------------------------------------
! Set FF data parameters
          CALL SET_FF_PAR (ML_SUPER_HANDLE%TAG_LIST,ML_SUPER_HANDLE%TIM, &
                           ML_SUPER_HANDLE%PAR_SUP_HANDLE,ML_SUPER_HANDLE%MLSHM, &
                           ML_SUPER_HANDLE%MEM1,ML_SUPER_HANDLE%LOGF,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM, &
                           ML_SUPER_HANDLE%AB%ABCAR,ML_SUPER_HANDLE%AB%IU1, &
                           ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP, &
                           ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP,.TRUE.,3)
!-----------------------------------------------------------------
! SET_AB_PAR
!-----------------------------------------------------------------
! Before calling SET_AB_PAR ABN%NITYP needs to be allocated
          IF (ALLOCATED(ML_SUPER_HANDLE%ABN%NITYP)) DEALLOCATE(ML_SUPER_HANDLE%ABN%NITYP)
          ALLOCATE(ML_SUPER_HANDLE%ABN%NITYP(ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS))
! Set ab initio data parameters
          CALL SET_AB_PAR (ML_SUPER_HANDLE%ML_INPUT_HANDLE, &
               ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%TIM,ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD,ML_SUPER_HANDLE%MEM1, &
               ML_SUPER_HANDLE%LOGF,ML_SUPER_HANDLE%AB,ML_SUPER_HANDLE%ABN, &
               ML_SUPER_HANDLE%FF%CTIFOR,ML_SUPER_HANDLE%FF%IFF,3, &
               ML_SUPER_HANDLE%FFM%LB,ML_SUPER_HANDLE%FFM%LBHEAD,ML_SUPER_HANDLE%FF%LFMAT,ML_SUPER_HANDLE%FF%LGENFF, &
               ML_SUPER_HANDLE%FF%LLIST,ML_SUPER_HANDLE%FF%LLIST_NEW,ML_SUPER_HANDLE%FF%LMLONLY, &
               ML_SUPER_HANDLE%FFM%LNCONF,ML_SUPER_HANDLE%FFM%LNIONS,ML_SUPER_HANDLE%FF%LOVER,ML_SUPER_HANDLE%FF%LPS, &
               ML_SUPER_HANDLE%FF%LSIG,ML_SUPER_HANDLE%FF%LSOAP,ML_SUPER_HANDLE%FFM%MB,ML_SUPER_HANDLE%FFM%MBTOT, &
               ML_SUPER_HANDLE%FF%MCONF,ML_SUPER_HANDLE%FF%MTYP,ML_SUPER_HANDLE%FFM%NB, &
               ML_SUPER_HANDLE%FFM%NBTOT,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS,ML_SUPER_HANDLE%FF%NSYS, &
               ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP,ML_SUPER_HANDLE%FF%POMASS,ML_SUPER_HANDLE%FF%SZNAM2, &
               ML_SUPER_HANDLE%FF%TYPE)
! After that we deallocate ABN%NITYP again, it will be later allocated in SET_ABN_PAR
          IF (ALLOCATED(ML_SUPER_HANDLE%ABN%NITYP)) DEALLOCATE(ML_SUPER_HANDLE%ABN%NITYP)
!-----------------------------------------------------------------
! ADDITIONAL STUFF
!-----------------------------------------------------------------
! Set values for necessary variables for allocation
          ML_SUPER_HANDLE%ABN%MTYP_ALLOCATE=ML_SUPER_HANDLE%AB%MTYP
          ML_SUPER_HANDLE%ABN%MIONS_ALLOCATE=ML_SUPER_HANDLE%AB%MIONS
          ML_SUPER_HANDLE%ABN%MITYP_ALLOCATE=ML_SUPER_HANDLE%AB%MITYP
          ML_SUPER_HANDLE%ABN%NSYS_ALLOCATE=ML_SUPER_HANDLE%AB%NSYS
! First allocate some important variables
          IF (ALLOCATED(ML_SUPER_HANDLE%AB%MIONS_CONF)) DEALLOCATE(ML_SUPER_HANDLE%AB%MIONS_CONF)
          ALLOCATE(ML_SUPER_HANDLE%AB%MIONS_CONF(1:ML_SUPER_HANDLE%AB%MCONF))
          IF (ALLOCATED(ML_SUPER_HANDLE%AB%MTYP_CONF)) DEALLOCATE(ML_SUPER_HANDLE%AB%MTYP_CONF)
          ALLOCATE(ML_SUPER_HANDLE%AB%MTYP_CONF(1:ML_SUPER_HANDLE%AB%MCONF))
          IF (ALLOCATED(ML_SUPER_HANDLE%AB%MITYP_CONF)) DEALLOCATE(ML_SUPER_HANDLE%AB%MITYP_CONF)
          ALLOCATE(ML_SUPER_HANDLE%AB%MITYP_CONF(1:ML_SUPER_HANDLE%AB%MCONF))
! Set some variables that will be needed later
! Calculate the maximum number of elements, atoms and atoms per element for each configuration
! up to that configuration
          ! We also need to allocate a helping array that keeps track of the types
          ALLOCATE(TYPE_HELP(ML_SUPER_HANDLE%AB%MTYP_ALLOCATE))
          ML_SUPER_HANDLE%AB%MIONS_CONF=0
          ML_SUPER_HANDLE%AB%MTYP_CONF=0
          ML_SUPER_HANDLE%AB%MITYP_CONF=0
          ML_SUPER_HANDLE%AB%MIONS_CONF(1)=ML_SUPER_HANDLE%AB%NIONS(1)
          ML_SUPER_HANDLE%AB%MTYP_CONF(1)=ML_SUPER_HANDLE%AB%NTYP(1)
          ML_SUPER_HANDLE%AB%MITYP_CONF(1)=ML_SUPER_HANDLE%AB%NITYP(1,1)
          DO ITYP=1,ML_SUPER_HANDLE%AB%NTYP(1)
             TYPE_HELP(ITYP)=ML_SUPER_HANDLE%AB%TYPE(ML_SUPER_HANDLE%AB%LADD(ITYP,1))
             ML_SUPER_HANDLE%AB%MITYP_CONF(1)=MAX(ML_SUPER_HANDLE%AB%MITYP_CONF(1),ML_SUPER_HANDLE%AB%NITYP(ML_SUPER_HANDLE%AB%LADD(ITYP,1),1))
          ENDDO
          IF (ML_SUPER_HANDLE%AB%NCONF.GT.1) THEN
             DO IC=2,ML_SUPER_HANDLE%AB%NCONF
                ML_SUPER_HANDLE%AB%MIONS_CONF(IC)=MAX(ML_SUPER_HANDLE%AB%NIONS(IC),ML_SUPER_HANDLE%AB%MIONS_CONF(IC-1))
                ML_SUPER_HANDLE%AB%MTYP_CONF(IC)=ML_SUPER_HANDLE%AB%MTYP_CONF(IC-1)
                ML_SUPER_HANDLE%AB%MITYP_CONF(IC)=ML_SUPER_HANDLE%AB%MITYP_CONF(IC-1)
                DO ITYP=1,ML_SUPER_HANDLE%AB%NTYP(IC)
                   IHELP=0
                   DO JTYP=1,ML_SUPER_HANDLE%AB%MTYP_CONF(IC-1)
                      IF (ML_SUPER_HANDLE%AB%TYPE(ML_SUPER_HANDLE%AB%LADD(ITYP,IC)).EQ.TYPE_HELP(JTYP)) THEN
                         IHELP=1
                      ENDIF
                   ENDDO
                   IF (IHELP.EQ.0) THEN
                      ML_SUPER_HANDLE%AB%MTYP_CONF(IC)=ML_SUPER_HANDLE%AB%MTYP_CONF(IC)+1
                      TYPE_HELP(ML_SUPER_HANDLE%AB%MTYP_CONF(IC))=ML_SUPER_HANDLE%AB%TYPE(ML_SUPER_HANDLE%AB%LADD(ITYP,IC))
                   ENDIF
                   ML_SUPER_HANDLE%AB%MITYP_CONF(IC)=MAX(ML_SUPER_HANDLE%AB%NITYP(ML_SUPER_HANDLE%AB%LADD(ITYP,IC),IC),ML_SUPER_HANDLE%AB%MITYP_CONF(IC))
                ENDDO 
             ENDDO
          ENDIF
          ! Deallocate helping array
          DEALLOCATE(TYPE_HELP) 
! We need to back up some of the values here that maybe overwritten later
! Everything backed up into the helping structure MLH_TMP
          CALL CONSTRUCT_AB_FOR_HELPING_STRUCTURE(ML_SUPER_HANDLE,MLH_TMP)
! Here we need to set many values as if we would start from the beginning
          ML_SUPER_HANDLE%FF%NSYS=1
          ML_SUPER_HANDLE%AB%NSYS=1
          ML_SUPER_HANDLE%AB%NCONF=0
          ML_SUPER_HANDLE%FFM%NB=0 
          ML_SUPER_HANDLE%FFM%NBTOT=0
          ML_SUPER_HANDLE%FFM%NBTOT_NEW=0
          ML_SUPER_HANDLE%FFM%LB=0
          ML_SUPER_HANDLE%FFM%LBHEAD=0
          ML_SUPER_HANDLE%AB%NCONF_SYS=0
!-----------------------------------------------------------------
! SET_ABN_PAR
!-----------------------------------------------------------------
! Set parameters for new ab initio data
          ML_SUPER_HANDLE%ABN%MTYP_ALLOCATE=ML_SUPER_HANDLE%AB%MTYP_ALLOCATE
          ML_SUPER_HANDLE%ABN%MIONS_ALLOCATE=ML_SUPER_HANDLE%AB%MIONS_ALLOCATE
          ML_SUPER_HANDLE%ABN%MITYP_ALLOCATE=ML_SUPER_HANDLE%AB%MITYP_ALLOCATE
! Initialize data with the first ab initio data and do allocations
!          CALL SET_ABN_PAR(ML_SUPER_HANDLE%ABN,ML_SUPER_HANDLE%MEM1,3, &
!               ML_SUPER_HANDLE%FF%MCONF_NEW,ML_SUPER_HANDLE%FF%MTYP, &
!               ML_SUPER_HANDLE%AB%NIONS(1),ML_SUPER_HANDLE%AB%NITYP(:,1), &
!               ML_SUPER_HANDLE%FF%NSYS,ML_SUPER_HANDLE%AB%NTYP(1),ML_SUPER_HANDLE%FF%POMASS, &
!               ML_SUPER_HANDLE%FF%SZNAM2,ML_SUPER_HANDLE%AB%SZNAM2(1), &
!               ML_SUPER_HANDLE%FF%TYPE,ML_SUPER_HANDLE%FF%TYPE,.TRUE.)
          CALL SET_ABN_PAR_ISTART3(ML_SUPER_HANDLE%FF, &
                                   ML_SUPER_HANDLE%AB, &
                                   ML_SUPER_HANDLE%ABN, &
                                   MLH_TMP, &
                                   ML_SUPER_HANDLE%MEM1, &
                                   1, &
                                   ML_SUPER_HANDLE%AB%NITYP(:,1), &
                                   .TRUE.)
!-----------------------------------------------------------------
! INIT_LIST
!-----------------------------------------------------------------
! Initialize list of training data
! This takes care of FFM%NDATA, FFM%NTOTEN, FFM%NTIFOR, FFM%NTSIF, Etc.
          IF(ML_SUPER_HANDLE%FF%LLIST) THEN
             CALL INIT_LIST (ML_SUPER_HANDLE%TIM,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM, &
                             ML_SUPER_HANDLE%AB%NCONF,ML_SUPER_HANDLE%AB%NIONS)
          ENDIF
!-----------------------------------------------------------------
! INIT_LIST_NEW
!-----------------------------------------------------------------
! Initialize list of candidates of new training data.
          IF(ML_SUPER_HANDLE%FF%LLIST_NEW) THEN
             CALL INIT_LIST_NEW (ML_SUPER_HANDLE%TIM,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM,ML_SUPER_HANDLE%ABN%NIONS)
          ENDIF
        END SUBROUTINE INIT_FF_ISTART_3_FIRST_CALL

!****************************************************************************************************
! Subsequent initilializations needed by  ISTART_3_DRIVER if number of atoms or elements change
!****************************************************************************************************
        SUBROUTINE INIT_FF_ISTART_3_SUBSEQUENT_CALLS(ML_SUPER_HANDLE,MLH_TMP,ICONF)
           USE ML_FF_PREC
           USE ML_FF_CONSTANT
           USE LOGFILE
           USE IOHANDLE
           USE MPI_DATA
           USE TIME_CHECK
           USE ABINITIO_DATA
           USE FORCE_FIELD
           USE MEMORY_USAGE
           USE RANDOM_NUMBER_GENERATOR
           USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE, ML_FILESTREAM_GLOBAL
           IMPLICIT NONE
! Input variables
           TYPE (ML_SUPER_TYPE)           :: ML_SUPER_HANDLE
           TYPE (ML_SUPER_TYPE)           :: MLH_TMP
           INTEGER, INTENT(IN)            :: ICONF
! Local variables
           INTEGER                        :: IC
           INTEGER                        :: IFLAG
           INTEGER                        :: INTYP
           INTEGER                        :: INTYP0
           INTEGER                        :: IRB
           INTEGER                        :: ISYS
           INTEGER                        :: ITYP
           INTEGER                        :: JNTYP0
           INTEGER                        :: JJNTYP0
           INTEGER                        :: JRB
           INTEGER                        :: JTYP
           INTEGER                        :: L
           INTEGER                        :: MIONS_TMP
           INTEGER                        :: MITYP_TMP
           INTEGER                        :: MTYP_TMP
           INTEGER, ALLOCATABLE           :: NITYP_HELP(:)
           CHARACTER(LEN=40)              :: SZNAM2_HELP2(1:MLH_TMP%AB%NSYS)
           CHARACTER (LEN=2), ALLOCATABLE :: TYPE_HELP(:)
!-----------------------------------------------------------------
! Set important variables
!-----------------------------------------------------------------
! Set neccessary variables for the current configuration
           ML_SUPER_HANDLE%AB%MIONS=ML_SUPER_HANDLE%AB%MIONS_CONF(ICONF)
           ML_SUPER_HANDLE%AB%MTYP=ML_SUPER_HANDLE%AB%MTYP_CONF(ICONF)
           ML_SUPER_HANDLE%AB%MITYP=ML_SUPER_HANDLE%AB%MITYP_CONF(ICONF)
           ML_SUPER_HANDLE%FF%MIONS=ML_SUPER_HANDLE%AB%MIONS
           ML_SUPER_HANDLE%FF%MTYP=ML_SUPER_HANDLE%AB%MTYP
           ML_SUPER_HANDLE%FF%MITYP=ML_SUPER_HANDLE%AB%MITYP
!-----------------------------------------------------------------
! Set important variables
!-----------------------------------------------------------------
! Possibly update the number of systems (NSYS) here
           ML_SUPER_HANDLE%AB%NSYS=1
           SZNAM2_HELP2=""
           SZNAM2_HELP2(1)= ML_SUPER_HANDLE%AB%SZNAM2(1)
           DO IC=1,ICONF
              IFLAG=1
              DO ISYS=1,ML_SUPER_HANDLE%AB%NSYS
                 IF (ML_SUPER_HANDLE%AB%SZNAM2(IC).EQ.SZNAM2_HELP2(ISYS)) THEN
                    IFLAG=0
                    EXIT
                 ENDIF
              ENDDO
              IF (IFLAG.EQ.1) THEN
                 ML_SUPER_HANDLE%AB%NSYS=ML_SUPER_HANDLE%AB%NSYS+1
                 SZNAM2_HELP2(ML_SUPER_HANDLE%AB%NSYS)=ML_SUPER_HANDLE%AB%SZNAM2(IC)
              ENDIF
           ENDDO
! Also update the number of systems for FF variables
           ML_SUPER_HANDLE%FF%NSYS=ML_SUPER_HANDLE%AB%NSYS
!-----------------------------------------------------------------
! SET_FF_PAR - SET_FF_PAR_ARRAY_ISTART1
!-----------------------------------------------------------------
           ML_SUPER_HANDLE%FFM%NBTOT=0
           DO INTYP0=1, ML_SUPER_HANDLE%FF%MTYP
              ML_SUPER_HANDLE%FFM%NBTOT=ML_SUPER_HANDLE%FFM%NBTOT+ML_SUPER_HANDLE%FFM%NB(INTYP0)
           ENDDO
           ML_SUPER_HANDLE%FFM%LBHEAD(1)=0
           DO INTYP0=2, ML_SUPER_HANDLE%FF%MTYP
              ML_SUPER_HANDLE%FFM%LBHEAD(INTYP0)=ML_SUPER_HANDLE%FFM%LBHEAD(INTYP0-1)+ML_SUPER_HANDLE%FFM%NB(INTYP0-1)
           ENDDO
! Determination of the numbers of variables related to the radial descriptor.
           IF(ML_SUPER_HANDLE%FFM%W1.GT.0.0_q) THEN
              ML_SUPER_HANDLE%FFM%NNVAR1=ML_SUPER_HANDLE%FFM%MRB1*ML_SUPER_HANDLE%FF%MTYP
           ENDIF
! Determination of the numbers of variables related to the angular descriptor.
           IF(ML_SUPER_HANDLE%FFM%W2.GT.0.0_q) THEN
              ML_SUPER_HANDLE%FFM%NNVAR2=0
              ML_SUPER_HANDLE%FFM%LVAR=0
              ML_SUPER_HANDLE%FFM%LFLAG_VAR=.FALSE.
              ML_SUPER_HANDLE%FFM%WVAR=0.0_q
              IF(ML_SUPER_HANDLE%FFM%LSIC) THEN
                 ML_SUPER_HANDLE%FFM%NNVAR2_SIC=0
                 ML_SUPER_HANDLE%FFM%LVAR_SIC=0
                 ML_SUPER_HANDLE%FFM%LFLAG_VAR_SIC=.FALSE.
              ENDIF
              DO INTYP0=1, ML_SUPER_HANDLE%FF%MTYP
                 DO JNTYP0=1, ML_SUPER_HANDLE%FF%MTYP
                    DO JJNTYP0=1, ML_SUPER_HANDLE%FF%MTYP
                       DO L=0, ML_SUPER_HANDLE%FFM%LMAX2
                          DO IRB=1, ML_SUPER_HANDLE%FFM%NRB2(L)
                             DO JRB=IRB, ML_SUPER_HANDLE%FFM%NRB2(L)
                                ML_SUPER_HANDLE%FFM%NNVAR2(INTYP0)=ML_SUPER_HANDLE%FFM%NNVAR2(INTYP0)+1
                                ML_SUPER_HANDLE%FFM%LVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=ML_SUPER_HANDLE%FFM%NNVAR2(INTYP0)
                                ML_SUPER_HANDLE%FFM%LFLAG_VAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=.TRUE.
                                IF(JRB.EQ.IRB) THEN
                                   ML_SUPER_HANDLE%FFM%WVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=1.0_q
                                ELSE
                                   ML_SUPER_HANDLE%FFM%WVAR(JRB,IRB,L,JJNTYP0,JNTYP0,INTYP0)=SQRT(2.0_q)
                                ENDIF
                                IF(ML_SUPER_HANDLE%FFM%LSIC) THEN
                                   IF(JJNTYP0.EQ.JNTYP0) THEN
                                      ML_SUPER_HANDLE%FFM%NNVAR2_SIC(INTYP0)=ML_SUPER_HANDLE%FFM%NNVAR2_SIC(INTYP0)+1
                                      ML_SUPER_HANDLE%FFM%LVAR_SIC(JRB,IRB,L,JNTYP0,INTYP0)=ML_SUPER_HANDLE%FFM%NNVAR2_SIC(INTYP0)
                                      ML_SUPER_HANDLE%FFM%LFLAG_VAR_SIC(JRB,IRB,L,JNTYP0,INTYP0)=.TRUE.
                                   ENDIF
                                ENDIF
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
          
!-----------------------------------------------------------------
! SET_ABN_PAR
!-----------------------------------------------------------------
! Set parameters for new ab initio data, no allocation is done here
! We need a helping array for the types of the current structures, since it can be 
! different from all types. This is needed to determin the location in the arrays of the type.
           IF (ALLOCATED(TYPE_HELP)) DEALLOCATE(TYPE_HELP)
           IF (ALLOCATED(NITYP_HELP)) DEALLOCATE(NITYP_HELP)
           ALLOCATE(TYPE_HELP(1:MLH_TMP%AB%NTYP(ICONF)))
           ALLOCATE(NITYP_HELP(1:MLH_TMP%AB%NTYP(ICONF)))
           ITYP=0
           DO JTYP=1,MLH_TMP%AB%MTYP
              IF (MLH_TMP%AB%LADD(JTYP,ICONF).GT.0) THEN
                 IF (MLH_TMP%AB%NITYP(MLH_TMP%AB%LADD(JTYP,ICONF),ICONF).GT.0) THEN
                    ITYP=ITYP+1
                    TYPE_HELP(ITYP)=MLH_TMP%AB%TYPE(MLH_TMP%AB%LADD(JTYP,ICONF))
                    NITYP_HELP(ITYP)=MLH_TMP%AB%NITYP(MLH_TMP%AB%LADD(JTYP,ICONF),ICONF)
                 ENDIF
              ENDIF
           ENDDO
!           CALL SET_ABN_PAR(ML_SUPER_HANDLE%ABN,ML_SUPER_HANDLE%MEM1,3, &
!                ML_SUPER_HANDLE%FF%MCONF_NEW,ML_SUPER_HANDLE%FF%MTYP, &
!                MLH_TMP%AB%NIONS(ICONF),NITYP_HELP(:), &
!                ML_SUPER_HANDLE%FF%NSYS,MLH_TMP%AB%NTYP(ICONF),ML_SUPER_HANDLE%FF%POMASS, &
!                ML_SUPER_HANDLE%FF%SZNAM2,MLH_TMP%AB%SZNAM2(ICONF), &
!                ML_SUPER_HANDLE%FF%TYPE,TYPE_HELP,.FALSE.)
           CALL SET_ABN_PAR_ISTART3(ML_SUPER_HANDLE%FF, &
                                    ML_SUPER_HANDLE%AB, &
                                    ML_SUPER_HANDLE%ABN, &
                                    MLH_TMP, &
                                    ML_SUPER_HANDLE%MEM1, &
                                    ICONF, &
                                    NITYP_HELP(:), &
                                    .FALSE.)
! Deallocate helping arrays
           IF (ALLOCATED(TYPE_HELP)) DEALLOCATE(TYPE_HELP)
           IF (ALLOCATED(NITYP_HELP)) DEALLOCATE(NITYP_HELP)
!-----------------------------------------------------------------
! SET_ABN_PAR
!-----------------------------------------------------------------
           ML_SUPER_HANDLE%FF%LLIST_NEW=.TRUE.
           IF(ML_SUPER_HANDLE%FF%LLIST_NEW) THEN
              CALL INIT_LIST_NEW (ML_SUPER_HANDLE%TIM,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM,ML_SUPER_HANDLE%ABN%NIONS)
           ENDIF

!-----------------------------------------------------------------
! Set this variable to 0 overall otherwise some wrong 
! elements can be copied
!-----------------------------------------------------------------
           ML_SUPER_HANDLE%FFM%FMAT_NEW=0.0_q

!-----------------------------------------------------------------
! Setting control variables similar to ISTART=1
!-----------------------------------------------------------------
            ML_SUPER_HANDLE%FF%IFF          =.FALSE.
            ML_SUPER_HANDLE%FF%LDOUB        =.FALSE.
            ML_SUPER_HANDLE%FF%LFMAT     =.TRUE.
            ML_SUPER_HANDLE%FF%LGENFF    =.TRUE.
            ML_SUPER_HANDLE%FF%LLIST     =.TRUE.
            ML_SUPER_HANDLE%FF%LLIST_NEW =.TRUE.
            ML_SUPER_HANDLE%FF%LOVER     =.TRUE.
            ML_SUPER_HANDLE%FF%LPS       =.TRUE.
            ML_SUPER_HANDLE%FF%LSAMPLE   =.TRUE.
            ML_SUPER_HANDLE%FF%LSIG      =.TRUE.
            ML_SUPER_HANDLE%FF%LSOAP     =.TRUE.
!-----------------------------------------------------------------
! The power spectrum needs to be recalculated, escpecially
! if element types and numbers have changed.
! The first time we don't need to call this (ICONF.NE.1).
!-----------------------------------------------------------------
            IF ((ICONF.NE.1).AND.(ML_SUPER_HANDLE%FF%LPS)) THEN
               CALL INIT_PS (ML_SUPER_HANDLE%TIM,ML_SUPER_HANDLE%PAR_SUP_HANDLE, &
                    ML_SUPER_HANDLE%MLSHM,ML_SUPER_HANDLE%MEM1,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM, &
                    ML_SUPER_HANDLE%AB%A,ML_SUPER_HANDLE%AB%B, &
                    ML_SUPER_HANDLE%AB%LADD,ML_SUPER_HANDLE%AB%LNIONS,ML_SUPER_HANDLE%AB%LNTYP, &
                    ML_SUPER_HANDLE%AB%MCONF,ML_SUPER_HANDLE%AB%MIONS,ML_SUPER_HANDLE%AB%MITYP, &
                    ML_SUPER_HANDLE%AB%MTYP,ML_SUPER_HANDLE%AB%NCONF,ML_SUPER_HANDLE%AB%NIONS, &
                    ML_SUPER_HANDLE%AB%NITYP,ML_SUPER_HANDLE%AB%NTYP,ML_SUPER_HANDLE%AB%POSION)
               CALL STORE_PS_IN_ALL_PROCESSES (ML_SUPER_HANDLE%PAR_SUP_HANDLE,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM,ML_SUPER_HANDLE%FFM%NB)
            ENDIF

            IF(ML_SUPER_HANDLE%FF%LSOAP) THEN
               CALL INIT_SOAP (ML_SUPER_HANDLE%TIM,ML_SUPER_HANDLE%PAR_SUP_HANDLE,ML_SUPER_HANDLE%MEM1,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM)
            ENDIF
        END SUBROUTINE INIT_FF_ISTART_3_SUBSEQUENT_CALLS
#endif

!****************************************************************************************************
! Usual initialization driver
!****************************************************************************************************
        SUBROUTINE ISTART_0_1_2_DRIVER(ML_SUPER_HANDLE)
          USE ML_FF_PREC
          USE ML_FF_CONSTANT
          USE LOGFILE
          USE IOHANDLE
          USE MPI_DATA
          USE TIME_CHECK
          USE ABINITIO_DATA
          USE FORCE_FIELD
          USE MEMORY_USAGE
          USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE, ML_FILESTREAM_GLOBAL, ML_IO_WRITE
          IMPLICIT NONE
          TYPE (ML_SUPER_TYPE) :: ML_SUPER_HANDLE
          INTEGER :: ICONF, ITYP
#ifndef ML_LEARNING_NOT_ACTIVE
! Check size of FF%MCONF_ALLOCATE, FF%MB... from ML_AB file if ML_FF_ISTART=1,3
          IF (ML_SUPER_HANDLE%FF%ISTART.EQ.1) THEN
              CALL READ_ABCAR_RESET_SIZES(ML_SUPER_HANDLE%AB,ML_SUPER_HANDLE%FF, &
                                          ML_SUPER_HANDLE%FFM,ML_SUPER_HANDLE%ML_INPUT_HANDLE, &
                                          ML_SUPER_HANDLE%TAG_LIST)
          ENDIF  
#endif
! Check format type of ML_FF file for reading (compatibility to non binary format)
          IF (ML_SUPER_HANDLE%FF%ISTART.EQ.2) THEN
             CALL CHECK_ML_FF_FORMAT(ML_SUPER_HANDLE%FF)
             ! Set version number to 0.0.0 for no version number like in non binary format
             IF (.NOT.ML_SUPER_HANDLE%FF%ML_FF_BINARY) THEN
                ML_SUPER_HANDLE%FF%ML_FF_VERSION_READ(1)=0
                ML_SUPER_HANDLE%FF%ML_FF_VERSION_READ(2)=0
                ML_SUPER_HANDLE%FF%ML_FF_VERSION_READ(3)=0
             ENDIF
          ENDIF
! Define element order
          CALL DEF_ELEMENT_ORDER_FF (ML_SUPER_HANDLE%TIM, &
                                     ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM,ML_SUPER_HANDLE%FF%FFCAR, &
                                     ML_SUPER_HANDLE%AB%ABCAR, &
                                     ML_SUPER_HANDLE%ML_INPUT_HANDLE%EATOM_VASP, &
                                     ML_SUPER_HANDLE%AB%IU1,ML_SUPER_HANDLE%FF%IU1, &
                                     ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP, &
                                     ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP, &
                                     ML_SUPER_HANDLE%ML_INPUT_HANDLE%POMASS, &
                                     ML_SUPER_HANDLE%ML_INPUT_HANDLE%TYPE,.TRUE.,ML_SUPER_HANDLE%FF%ISTART)
! Define system order.
          CALL DEF_SYSTEM_ORDER_FF (ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM,ML_SUPER_HANDLE%AB%ABCAR, &
                                    ML_SUPER_HANDLE%AB%IU1,ML_SUPER_HANDLE%ML_INPUT_HANDLE%SZNAM2,.TRUE.,ML_SUPER_HANDLE%FF%ISTART)
! Estimate necessary memory
          CALL ESTIMATE_MEMORY (ML_SUPER_HANDLE%TIM, &
               ML_SUPER_HANDLE%PAR_SUP_HANDLE,ML_SUPER_HANDLE%MEM0,ML_SUPER_HANDLE%MEM1, &
               ML_SUPER_HANDLE%LOGF,ML_SUPER_HANDLE%ML_INPUT_HANDLE%A,ML_SUPER_HANDLE%ML_INPUT_HANDLE%B, &
               ML_SUPER_HANDLE%FFM%IBROAD1,ML_SUPER_HANDLE%FFM%IBROAD2,ML_SUPER_HANDLE%FF%ISTART, &
               ML_SUPER_HANDLE%FFM%LCOUPLE,ML_SUPER_HANDLE%FFM%LMAX1,ML_SUPER_HANDLE%FFM%LMAX2, &
               ML_SUPER_HANDLE%FF%LMLMB,ML_SUPER_HANDLE%FF%LMLONLY,ML_SUPER_HANDLE%FFM%LSIC,ML_SUPER_HANDLE%FFM%LSUPERVEC, &
               ML_SUPER_HANDLE%FFM%MB,ML_SUPER_HANDLE%FF%MCONF,ML_SUPER_HANDLE%FF%MCONF_NEW,ML_SUPER_HANDLE%FF%MIONS, &
               ML_SUPER_HANDLE%FF%MITYP,ML_SUPER_HANDLE%FFM%MSPL1,ML_SUPER_HANDLE%FFM%MSPL2,ML_SUPER_HANDLE%FF%MTYP, &
               ML_SUPER_HANDLE%FFM%NHYP1,ML_SUPER_HANDLE%FFM%NHYP2,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP, &
               ML_SUPER_HANDLE%FFM%NR1,ML_SUPER_HANDLE%FFM%NR2,ML_SUPER_HANDLE%FFM%MRB1,ML_SUPER_HANDLE%FFM%MRB2, &
               ML_SUPER_HANDLE%FF%NSW,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP, &
               ML_SUPER_HANDLE%ML_INPUT_HANDLE%OMEGA,ML_SUPER_HANDLE%FFM%RCUT1,ML_SUPER_HANDLE%FFM%RCUT2, &
               ML_SUPER_HANDLE%FFM%W1,ML_SUPER_HANDLE%FFM%W2)
! Set FF data parameters
          CALL SET_FF_PAR (ML_SUPER_HANDLE%TAG_LIST,ML_SUPER_HANDLE%TIM, &
                           ML_SUPER_HANDLE%PAR_SUP_HANDLE,ML_SUPER_HANDLE%MLSHM, &
                           ML_SUPER_HANDLE%MEM1,ML_SUPER_HANDLE%LOGF,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM, &
                           ML_SUPER_HANDLE%AB%ABCAR,ML_SUPER_HANDLE%AB%IU1, &
                           ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP, &
                           ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP,.TRUE.,ML_SUPER_HANDLE%FF%ISTART)
! Set parameters for new ab initio data
          CALL SET_ABN_PAR(ML_SUPER_HANDLE%ABN,ML_SUPER_HANDLE%MEM1,ML_SUPER_HANDLE%FF%ISTART, &
               ML_SUPER_HANDLE%FF%MCONF_NEW,ML_SUPER_HANDLE%FF%MTYP, &
               ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP, &
               ML_SUPER_HANDLE%FF%NSYS,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP,ML_SUPER_HANDLE%ML_INPUT_HANDLE%POMASS, &
               ML_SUPER_HANDLE%FF%SZNAM2,ML_SUPER_HANDLE%ML_INPUT_HANDLE%SZNAM2, &
               ML_SUPER_HANDLE%FF%TYPE,ML_SUPER_HANDLE%ML_INPUT_HANDLE%TYPE,.TRUE.)
! Set ab initio data parameters
          CALL SET_AB_PAR (ML_SUPER_HANDLE%ML_INPUT_HANDLE, &
               ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%TIM,ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD,ML_SUPER_HANDLE%MEM1, &
               ML_SUPER_HANDLE%LOGF,ML_SUPER_HANDLE%AB,ML_SUPER_HANDLE%ABN, &
               ML_SUPER_HANDLE%FF%CTIFOR,ML_SUPER_HANDLE%FF%IFF,ML_SUPER_HANDLE%FF%ISTART, &
               ML_SUPER_HANDLE%FFM%LB,ML_SUPER_HANDLE%FFM%LBHEAD,ML_SUPER_HANDLE%FF%LFMAT,ML_SUPER_HANDLE%FF%LGENFF, &
               ML_SUPER_HANDLE%FF%LLIST,ML_SUPER_HANDLE%FF%LLIST_NEW,ML_SUPER_HANDLE%FF%LMLONLY, &
               ML_SUPER_HANDLE%FFM%LNCONF,ML_SUPER_HANDLE%FFM%LNIONS,ML_SUPER_HANDLE%FF%LOVER,ML_SUPER_HANDLE%FF%LPS, &
               ML_SUPER_HANDLE%FF%LSIG,ML_SUPER_HANDLE%FF%LSOAP,ML_SUPER_HANDLE%FFM%MB,ML_SUPER_HANDLE%FFM%MBTOT, &
               ML_SUPER_HANDLE%FF%MCONF,ML_SUPER_HANDLE%FF%MTYP,ML_SUPER_HANDLE%FFM%NB, &
               ML_SUPER_HANDLE%FFM%NBTOT,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS,ML_SUPER_HANDLE%FF%NSYS, &
               ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP,ML_SUPER_HANDLE%FF%POMASS,ML_SUPER_HANDLE%FF%SZNAM2, &
               ML_SUPER_HANDLE%FF%TYPE)
#ifndef ML_LEARNING_NOT_ACTIVE
! Initialization of angular descriptor if ML_SUPER_HANDLE%FF%LPS=.TRUE.
          IF(ML_SUPER_HANDLE%FF%LPS) THEN
             CALL INIT_PS (ML_SUPER_HANDLE%TIM,ML_SUPER_HANDLE%PAR_SUP_HANDLE, &
                  ML_SUPER_HANDLE%MLSHM,ML_SUPER_HANDLE%MEM1,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM, &
                  ML_SUPER_HANDLE%AB%A,ML_SUPER_HANDLE%AB%B, &
                  ML_SUPER_HANDLE%AB%LADD,ML_SUPER_HANDLE%AB%LNIONS,ML_SUPER_HANDLE%AB%LNTYP, &
                  ML_SUPER_HANDLE%AB%MCONF,ML_SUPER_HANDLE%AB%MIONS,ML_SUPER_HANDLE%AB%MITYP, &
                  ML_SUPER_HANDLE%AB%MTYP,ML_SUPER_HANDLE%AB%NCONF,ML_SUPER_HANDLE%AB%NIONS, &
                  ML_SUPER_HANDLE%AB%NITYP,ML_SUPER_HANDLE%AB%NTYP,ML_SUPER_HANDLE%AB%POSION)
          ENDIF 
! Sparsification of descriptors and recomputations of angular descriptors.
! Recomputation is necessary for renormalizations.
          IF(ML_SUPER_HANDLE%FF%LPS.AND.ML_SUPER_HANDLE%FFM%LSPARSDES) THEN
             CALL SPARS_DES_CUR (ML_SUPER_HANDLE%LAPACKWORK,ML_SUPER_HANDLE%PAR_SUP_HANDLE,ML_SUPER_HANDLE%LOGF, &
                  ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM,ML_SUPER_HANDLE%FFM%NB)
             CALL INIT_PS (ML_SUPER_HANDLE%TIM,ML_SUPER_HANDLE%PAR_SUP_HANDLE, &
                  ML_SUPER_HANDLE%MLSHM,ML_SUPER_HANDLE%MEM1,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM, &
                  ML_SUPER_HANDLE%AB%A,ML_SUPER_HANDLE%AB%B,ML_SUPER_HANDLE%AB%LADD, &
                  ML_SUPER_HANDLE%AB%LNIONS,ML_SUPER_HANDLE%AB%LNTYP,ML_SUPER_HANDLE%AB%MCONF, &
                  ML_SUPER_HANDLE%AB%MIONS,ML_SUPER_HANDLE%AB%MITYP,ML_SUPER_HANDLE%AB%MTYP, &
                  ML_SUPER_HANDLE%AB%NCONF,ML_SUPER_HANDLE%AB%NIONS,ML_SUPER_HANDLE%AB%NITYP, &
                  ML_SUPER_HANDLE%AB%NTYP,ML_SUPER_HANDLE%AB%POSION)
          ENDIF
! Do first print out of loopsection header here
          CALL LOGSECTION_LOOP_HEADER(ML_SUPER_HANDLE%FF, ML_SUPER_HANDLE%LOGF)
#endif
! Store angular descriptor for basis sets in following conditions.
! 1) ML_SUPER_HANDLE%FF%LPS=.TRUE. In this condition, the machine possesses angular descriptor,
!    and user specifies that the machine needs to store the angular descriptor in all processes.
! 2) ML_SUPER_HANDLE%FF%LMLONLY=.TRUE. (ML_SUPER_HANDLE%FF%ISTART=2). In this condition, the machine also possesses angular descriptors
!    which was read from ML_FF file.
          IF(ML_SUPER_HANDLE%FF%LPS.OR.ML_SUPER_HANDLE%FF%LMLONLY) THEN
             CALL STORE_PS_IN_ALL_PROCESSES (ML_SUPER_HANDLE%PAR_SUP_HANDLE,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM,ML_SUPER_HANDLE%FFM%NB)
          ENDIF
#ifndef ML_LEARNING_NOT_ACTIVE
! If ML_SUPER_HANDLE%FF%LSOAP=.TRUE., initialize SOAP matrix.
          IF(ML_SUPER_HANDLE%FF%LSOAP) THEN
             CALL INIT_SOAP (ML_SUPER_HANDLE%TIM,ML_SUPER_HANDLE%PAR_SUP_HANDLE,ML_SUPER_HANDLE%MEM1,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM)
          ENDIF
! Initialize list of training data
          IF(ML_SUPER_HANDLE%FF%LLIST) THEN
             CALL INIT_LIST (ML_SUPER_HANDLE%TIM,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM, &
                             ML_SUPER_HANDLE%AB%NCONF,ML_SUPER_HANDLE%AB%NIONS)
          ENDIF
! Initialize list of candidates of new training data.
          IF(ML_SUPER_HANDLE%FF%LLIST_NEW) THEN
             CALL INIT_LIST_NEW (ML_SUPER_HANDLE%TIM,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM,ML_SUPER_HANDLE%ABN%NIONS)
          ENDIF
! Initialization of design matrix
          IF(ML_SUPER_HANDLE%FF%LFMAT) THEN
             CALL INIT_FMAT (ML_SUPER_HANDLE%AB,ML_SUPER_HANDLE%PAR_SUP_HANDLE, &
                  ML_SUPER_HANDLE%MLSHM,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM,ML_SUPER_HANDLE%AB%A, &
                  ML_SUPER_HANDLE%AB%B,ML_SUPER_HANDLE%AB%LADD, &
                  ML_SUPER_HANDLE%AB%LNITYP,ML_SUPER_HANDLE%AB%LNIONS,ML_SUPER_HANDLE%AB%LNTYP, &
                  ML_SUPER_HANDLE%AB%MCONF,ML_SUPER_HANDLE%AB%MITYP,ML_SUPER_HANDLE%AB%MIONS, &
                  ML_SUPER_HANDLE%AB%MTYP,ML_SUPER_HANDLE%AB%NCONF,ML_SUPER_HANDLE%AB%NITYP, &
                  ML_SUPER_HANDLE%AB%NIONS,ML_SUPER_HANDLE%AB%NTYP, &
                  ML_SUPER_HANDLE%AB%POSION,ML_SUPER_HANDLE%AB%OMEGA)
          ENDIF
! Output the first step 
          CALL LOGLINE_STATUS(ML_SUPER_HANDLE%FF, ML_SUPER_HANDLE%LOGF, .TRUE.)
! Calculate the variance of the ab initio data
          IF(ML_SUPER_HANDLE%FF%LSIG) THEN
             CALL CAL_SIG (ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%AB%ISYS,ML_SUPER_HANDLE%AB%LADD, &
                           ML_SUPER_HANDLE%AB%MCONF,ML_SUPER_HANDLE%AB%MIONS, &
                           ML_SUPER_HANDLE%AB%MTYP,ML_SUPER_HANDLE%AB%NCONF, &
                           ML_SUPER_HANDLE%AB%NCONF_SYS,ML_SUPER_HANDLE%AB%NIONS, &
                           ML_SUPER_HANDLE%AB%NITYP,ML_SUPER_HANDLE%AB%NSYS, &
                           ML_SUPER_HANDLE%AB%NTYP,ML_SUPER_HANDLE%AB%TIFOR, &
                           ML_SUPER_HANDLE%AB%TOTEN,ML_SUPER_HANDLE%AB%TSIF)
          ENDIF
! Calculate the variance of the ab initio data
! Sparsification of stored basis set data if available
          IF(ML_SUPER_HANDLE%FF%LOVER) THEN
             CALL SPARS_ICUR (ML_SUPER_HANDLE%FF%ML_EPS_LOW,ML_SUPER_HANDLE%TIM, &
                              ML_SUPER_HANDLE%LAPACKWORK,ML_SUPER_HANDLE%PAR_SUP_HANDLE, &
                              ML_SUPER_HANDLE%MEM1,ML_SUPER_HANDLE%LOGF, &
                              ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM, &
                              ML_SUPER_HANDLE%AB%A,ML_SUPER_HANDLE%AB%B, &
                              ML_SUPER_HANDLE%AB%ISYS,ML_SUPER_HANDLE%AB%LADD, &
                              ML_SUPER_HANDLE%AB%LB,ML_SUPER_HANDLE%AB%LNITYP, &
                              ML_SUPER_HANDLE%AB%LNIONS,ML_SUPER_HANDLE%AB%LNTYP, &
                              ML_SUPER_HANDLE%AB%LSYS,ML_SUPER_HANDLE%AB%MCONF, &
                              ML_SUPER_HANDLE%AB%MITYP,ML_SUPER_HANDLE%AB%MIONS, &
                              ML_SUPER_HANDLE%AB%MTYP,ML_SUPER_HANDLE%AB%NCONF, &
                              ML_SUPER_HANDLE%AB%NCONF_SYS,ML_SUPER_HANDLE%AB%NITYP, &
                              ML_SUPER_HANDLE%AB%NIONS,0,ML_SUPER_HANDLE%AB%NSYS, &
                              ML_SUPER_HANDLE%AB%NTYP,ML_SUPER_HANDLE%AB%POSION, &
                              ML_SUPER_HANDLE%AB%TIFOR, &
                              ML_SUPER_HANDLE%AB%TOTEN,ML_SUPER_HANDLE%AB%TSIF, &
                              ML_SUPER_HANDLE%AB%TYPE,ML_SUPER_HANDLE%AB%OMEGA, &
                              ML_SUPER_HANDLE%AB%SZNAM2)
          ENDIF
#endif
! For upcoming MD step, this subroutine needs to be called again here in order to store all the restructured C00 and PS arrays.
          IF(ML_SUPER_HANDLE%FF%LPS.OR.ML_SUPER_HANDLE%FF%LMLONLY) THEN
             CALL STORE_PS_IN_ALL_PROCESSES (ML_SUPER_HANDLE%PAR_SUP_HANDLE,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM,ML_SUPER_HANDLE%FFM%NB)
          ENDIF
! Sparsification of descriptors.

#ifndef ML_LEARNING_NOT_ACTIVE
! Calculate the variance of the ab initio data
          IF(ML_SUPER_HANDLE%FF%LSIG) THEN
             CALL CAL_SIG (ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%AB%ISYS,ML_SUPER_HANDLE%AB%LADD, &
                           ML_SUPER_HANDLE%AB%MCONF,ML_SUPER_HANDLE%AB%MIONS, &
                           ML_SUPER_HANDLE%AB%MTYP,ML_SUPER_HANDLE%AB%NCONF, &
                           ML_SUPER_HANDLE%AB%NCONF_SYS,ML_SUPER_HANDLE%AB%NIONS, &
                           ML_SUPER_HANDLE%AB%NITYP,ML_SUPER_HANDLE%AB%NSYS, &
                           ML_SUPER_HANDLE%AB%NTYP,ML_SUPER_HANDLE%AB%TIFOR, &
                           ML_SUPER_HANDLE%AB%TOTEN,ML_SUPER_HANDLE%AB%TSIF)
          ENDIF
! Generation of force field from the stored ab initio data if available
          IF(ML_SUPER_HANDLE%FF%LGENFF) THEN
             CALL GEN_FF (ML_SUPER_HANDLE%AB,ML_SUPER_HANDLE%TIM, &
                  ML_SUPER_HANDLE%LAPACKWORK,ML_SUPER_HANDLE%PAR_SUP_HANDLE, &
                  ML_SUPER_HANDLE%MEM1,ML_SUPER_HANDLE%LOGF, &
                  ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM, &
                  ML_SUPER_HANDLE%AB%A,ML_SUPER_HANDLE%AB%B, &
                  ML_SUPER_HANDLE%AB%ISYS,ML_SUPER_HANDLE%AB%LADD, &
                  ML_SUPER_HANDLE%AB%LNITYP,ML_SUPER_HANDLE%AB%LNIONS, &
                  ML_SUPER_HANDLE%AB%LNTYP,ML_SUPER_HANDLE%AB%MCONF, &
                  ML_SUPER_HANDLE%AB%MITYP,ML_SUPER_HANDLE%AB%MIONS, &
                  ML_SUPER_HANDLE%AB%MTYP,ML_SUPER_HANDLE%AB%NCONF, &
                  ML_SUPER_HANDLE%AB%NCONF_SYS,ML_SUPER_HANDLE%AB%NITYP, &
                  ML_SUPER_HANDLE%AB%NIONS,0,ML_SUPER_HANDLE%AB%NSYS, &
                  ML_SUPER_HANDLE%AB%NTYP,ML_SUPER_HANDLE%AB%POSION, &
                  ML_SUPER_HANDLE%AB%OMEGA)
          ENDIF
#endif
! If ML_SUPER_HANDLE%FF%IFF=.TRUE., make all processes own the regressin coefficients, inverse SOAP matrices and covariance matrix.
          IF(ML_SUPER_HANDLE%FF%IFF) THEN
             CALL STORE_WMAT_FFM_IN_ALL_PROCESS(ML_SUPER_HANDLE%PAR_SUP_HANDLE,ML_SUPER_HANDLE%MEM1,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM)
             CALL STORE_CMAT_FFM_IN_ALL_PROCESS(ML_SUPER_HANDLE%PAR_SUP_HANDLE,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%FFM)
          ENDIF
#ifndef ML_LEARNING_NOT_ACTIVE
! Calculate histogram if ML_SUPER_HANDLE%FF%LGENFF=.TRUE. and ML_SUPER_HANDLE%FF%IFF=.TRUE.
! In this condition, the machine tried to generate the force field, and the generation was successfully finished.
! Therefore, the program can evaluate the distributions of training data and regression results.
          IF(ML_SUPER_HANDLE%FF%LGENFF.AND.ML_SUPER_HANDLE%FF%IFF) THEN
             CALL HISTOGRAM (ML_SUPER_HANDLE%PAR_SUP_HANDLE, &
                  ML_SUPER_HANDLE%AB,ML_SUPER_HANDLE%ML_INPUT_HANDLE%DIR_APP,ML_SUPER_HANDLE%ML_INPUT_HANDLE%DIR_LEN)
          ENDIF
! Do the rest of output to the ML_LOGFILE
          IF (ML_SUPER_HANDLE%FF%LGENFF) THEN
             CALL LOGLINE_STDAB_ERR(ML_SUPER_HANDLE%FF, ML_SUPER_HANDLE%LOGF)
             CALL LOGLINE_BEE(ML_SUPER_HANDLE%FF,   &
                              ML_SUPER_HANDLE%LOGF, &
                              ML_SUPER_HANDLE%ABN,  &
                              .TRUE.)
          ENDIF
! Output results of regression if ML_SUPER_HANDLE%FF%IFF=.TRUE. and ML_SUPER_HANDLE%FF%LGENFF=.TRUE.
          IF(ML_SUPER_HANDLE%FF%LGENFF.AND.ML_SUPER_HANDLE%FF%IFF) THEN
             CALL OUTPUT_REG (ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD, &
                  ML_SUPER_HANDLE%ML_INPUT_HANDLE%DIR_APP,ML_SUPER_HANDLE%ML_INPUT_HANDLE%DIR_LEN, &
                  ML_SUPER_HANDLE%AB%IU3,ML_SUPER_HANDLE%AB%MCONF,ML_SUPER_HANDLE%AB%MIONS,ML_SUPER_HANDLE%AB%NCONF, &
                  ML_SUPER_HANDLE%AB%NIONS,0,ML_SUPER_HANDLE%AB%REGCAR,ML_SUPER_HANDLE%AB%TIFOR, &
                  ML_SUPER_HANDLE%AB%TIFOR_MB,ML_SUPER_HANDLE%AB%TOTEN,ML_SUPER_HANDLE%AB%TOTEN_MB, &
                  ML_SUPER_HANDLE%AB%TSIF,ML_SUPER_HANDLE%AB%TSIF_MB)
          ENDIF
! Output the ab initio data if ML_SUPER_HANDLE%FF%LGENFF=.TRUE. and ML_SUPER_HANDLE%FF%IFF=.TRUE.
          IF(ML_SUPER_HANDLE%FF%LGENFF.AND.ML_SUPER_HANDLE%FF%IFF) THEN
             CALL OUTPUT_ABNCAR (ML_SUPER_HANDLE%AB,ML_SUPER_HANDLE%FF,ML_SUPER_HANDLE%PAR_SUP_HANDLE%COMM_WORLD, &
                  ML_SUPER_HANDLE%AB%A,ML_SUPER_HANDLE%AB%ABNCAR,ML_SUPER_HANDLE%ML_INPUT_HANDLE%DIR_APP, &
                  ML_SUPER_HANDLE%ML_INPUT_HANDLE%DIR_LEN,ML_SUPER_HANDLE%FF%EATOM_REF,ML_SUPER_HANDLE%AB%IU2,ML_SUPER_HANDLE%AB%LADD, &
                  ML_SUPER_HANDLE%FFM%LNCONF,ML_SUPER_HANDLE%FFM%LNIONS,ML_SUPER_HANDLE%FFM%MB,ML_SUPER_HANDLE%AB%MCONF, &
                  ML_SUPER_HANDLE%AB%MITYP,ML_SUPER_HANDLE%AB%MIONS,ML_SUPER_HANDLE%AB%MTYP,ML_SUPER_HANDLE%FFM%NB,ML_SUPER_HANDLE%AB%NCONF, &
                  ML_SUPER_HANDLE%AB%NITYP,ML_SUPER_HANDLE%AB%NIONS,ML_SUPER_HANDLE%AB%NTYP,ML_SUPER_HANDLE%AB%POMASS, &
                  ML_SUPER_HANDLE%AB%POSION,ML_SUPER_HANDLE%AB%SZNAM2,ML_SUPER_HANDLE%AB%TIFOR,ML_SUPER_HANDLE%AB%TOTEN, &
                  ML_SUPER_HANDLE%AB%TSIF,ML_SUPER_HANDLE%AB%TYPE)
          ENDIF
#endif
        END SUBROUTINE ISTART_0_1_2_DRIVER



END MODULE ml_main_subroutines 

#else

MODULE ml_main_subroutines
        CONTAINS
        SUBROUTINE MACHINE_LEARNING_INIT_LIB (ML_SUPER_HANDLE)
        END SUBROUTINE MACHINE_LEARNING_INIT_LIB
        SUBROUTINE ML_TO_VASP_MACHINE_LEARNING_LIB (ML_SUPER_HANDLE,A,B,POTIM,NSTEP,OMEGA,POSION,VEL,SZNAM2,TOTEN,TIFOR,TSIF,LDO_AB_INITIO)
        END SUBROUTINE ML_TO_VASP_MACHINE_LEARNING_LIB
        SUBROUTINE VASP_TO_ML_MACHINE_LEARNING_LIB (ML_SUPER_HANDLE,A,B,DIR_APP,DIR_LEN,LABORT,NSTEP,OMEGA,POSION,SZNAM2,TOTEN,TIFOR,TSIF)
        END SUBROUTINE VASP_TO_ML_MACHINE_LEARNING_LIB
        SUBROUTINE MACHINE_LEARNING_FINISH_LIB (ML_SUPER_HANDLE,DIR_APP,DIR_LEN,NSTEP)
        END SUBROUTINE MACHINE_LEARNING_FINISH_LIB
END MODULE ml_main_subroutines 
#endif
