#include "symbol.inc"

!**********************************************************************
!
! this module is interfacing VASP to a machine learning code
! it is part of VASP and coverged by the general VASP 
! copy right
!
!**********************************************************************
MODULE ml_interface
  USE mpimy
  USE prec
  USE tutor, ONLY: vtutor, isError, RandomSeed, argument
#ifdef ML_AVAILABLE
  USE ml_main_subroutines
#endif
  IMPLICIT NONE
! this variable specifies whether the MLFF calculation is executed or not.
     LOGICAL :: ML_LMLFF
! this variable determines whether ab initio calculation is executed.
     LOGICAL :: LDO_AB_INITIO=.TRUE.
     LOGICAL :: LRUN_VASP_TO_ML=.FALSE.
! this variable specifies the results of the judgement of the necessity of ab initio calculations on new configurations.
! if ML_LMLFF=.FALSE., LDO_AB_INITIO_NEW=LDO_AB_INITIO=.TRUE.
     LOGICAL :: LDO_AB_INITIO_NEW=.TRUE.
! this variable determines whether the calculations are executed using only the
! force field or not.
     LOGICAL :: ML_FF_LMLONLY
  
  CONTAINS

!=======================================================================
!
! This subroutine reads variables necessary for the machine-learning force field
! generations from INCAR file.
!
!=======================================================================

   SUBROUTINE MACHINE_LEARNING_GENERAL_AND_INSTANCE_READER(IU5,IU0)
      USE reader_tags
#ifdef ML_AVAILABLE
      ! Warning: "tag" states are not set in this routine because the
      ! corresponding instance "ML_SUPER_TYPE%TAG_LIST" is not available here.
      USE ML_FF_TAGLIST, ONLY: TAGLIST_GET_ALT_TAGS
      USE ML_FF_STRUCT, ONLY: ML_FILESTREAM_GLOBAL, ML_TOTNUM_INSTANCES
#endif
      IMPLICIT NONE
      INTEGER, INTENT(IN)           :: IU5
      INTEGER, INTENT(IN)           :: IU0
      INTEGER                       :: IERR
      INTEGER                       :: N
      LOGICAL                       :: LOPEN
      CHARACTER(len=:), ALLOCATABLE :: ALT_TAGS(:)
      INTEGER                       :: I
      ! Set ML_LMLFF as .FALSE.
      ML_LMLFF=.FALSE.
      ! check if INCAR file exists
      CALL OPEN_INCAR_IF_FOUND(IU5, LOPEN)
      ! Read ML_LMLFF.
      !   ML_LMLFF = .FALSE.: No machine-learning force field generations
      !   ML_LMLFF = .TRUE. : Execute machine-learning force field
      !   generations
      ! Need to set alternatives manually here because the taglist module is
      ! not available here.
      ALT_TAGS = ['ML_LMLFF   ', 'ML_FF_LMLFF']
      DO I = 1, SIZE(ALT_TAGS)
         CALL PROCESS_INCAR(LOPEN,             &
                            IU0,               &
                            IU5,               &
                            TRIM(ALT_TAGS(I)), &
                            ML_LMLFF,          &
                            IERR,              &
                            WRITEXMLINCAR,     &
                            FOUNDNUMBER=N)
         ! Skip alternative tags if one is already found.
         IF (N > 0) EXIT
      END DO
#ifdef ML_AVAILABLE
      ! By default the number of instances is 1
      ML_TOTNUM_INSTANCES = 1
      ! At this point there is no need for multiple instances so we disable
      ! reading of this tag.
      !ALT_TAGS = TAGLIST_GET_ALT_TAGS('ML_TOTNUM_INSTANCES')
      !DO I = 1, SIZE(ALT_TAGS)
      !   CALL PROCESS_INCAR(LOPEN,               &
      !                      IU0,                 &
      !                      IU5,                 &
      !                      TRIM(ALT_TAGS(I)),   &
      !                      ML_TOTNUM_INSTANCES, &
      !                      IERR,                &
      !                      WRITEXMLINCAR,       &
      !                      FOUNDNUMBER=N)
      !   ! Skip alternative tags if one is already found.
      !   IF (N > 0) EXIT
      !END DO

      ! Global filestream variable
      ML_FILESTREAM_GLOBAL=500
      ALT_TAGS = TAGLIST_GET_ALT_TAGS('ML_FILESTREAM_START')
      DO I = 1, SIZE(ALT_TAGS)
         CALL PROCESS_INCAR(LOPEN,                &
                            IU0,                  &
                            IU5,                  &
                            TRIM(ALT_TAGS(I)),    &
                            ML_FILESTREAM_GLOBAL, &
                            IERR,                 &
                            WRITEXMLINCAR,        &
                            FOUNDNUMBER=N)
         ! Skip alternative tags if one is already found.
         IF (N > 0) EXIT
      END DO
#else
      IF (ML_LMLFF) THEN
         CALL vtutor%error("This VASP executable does not support machine &
                           &learning (ML_LMLFF=.TRUE.)! Try recompiling with &
                           &-DML_AVAILABLE")
      ENDIF
#endif
   END SUBROUTINE MACHINE_LEARNING_GENERAL_AND_INSTANCE_READER
#ifdef ML_AVAILABLE

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !> This subroutine reads variables necessary for the machine-learning force
   !> field generations from INCAR file.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE MACHINE_LEARNING_READER (TAG_LIST,         &
                                       FFM,              &
                                       FF,               &
                                       ML_INPUT_HANDLE,  &
                                       IU5,              &
                                       IU0,              &
                                       NTYP,             &
                                       IWAVPR,           &
                                       CURRENT_INSTANCE, &
                                       INSTANCES)
      USE reader_tags
      USE string, ONLY: str
      USE ml_ff_constant, ONLY: EUNIT, FUNIT, SUNIT, AUTOA
      USE ML_FF_TAGLIST
      USE ml_ff_struct, ONLY: FFM_PAR,     &
                              FF_PAR,      &
                              ML_FF_INPUT, &
                              ML_FILESTREAM_GLOBAL
      ! Input variables
      TYPE(TAGLIST), INTENT(INOUT)     :: TAG_LIST
      TYPE(FFM_PAR), INTENT(INOUT)     :: FFM
      TYPE(FF_PAR),  INTENT(INOUT)     :: FF
      TYPE(ML_FF_INPUT), INTENT(INOUT) :: ML_INPUT_HANDLE
      INTEGER,        INTENT(IN)       :: IU5
      INTEGER,        INTENT(IN)       :: IU0
      INTEGER,        INTENT(IN)       :: NTYP
      INTEGER,        INTENT(IN)       :: CURRENT_INSTANCE
      INTEGER,        INTENT(IN)       :: INSTANCES
      INTEGER,        INTENT(INOUT)    :: IWAVPR
      ! Local variables
      CHARACTER(len=:), ALLOCATABLE :: PREFIX
      INTEGER :: DTVALUES(8)
      INTEGER :: IDUMLONG
      INTEGER, PARAMETER :: SEED1_MAX = 900000000 
      INTEGER, PARAMETER :: SEED2_MAX = 1000000
      INTEGER :: IERR
      INTEGER :: IWAVPR_HELP
      LOGICAL :: LOPEN
      INTEGER :: M
      INTEGER :: N

      INTERFACE READ_TAG
         PROCEDURE READ_TAG_REAL
         PROCEDURE READ_TAG_INTEGER
         PROCEDURE READ_TAG_LOGICAL
      END INTERFACE READ_TAG

      INTERFACE READ_TAG_ARRAY
         PROCEDURE READ_TAG_ARRAY_REAL
         PROCEDURE READ_TAG_ARRAY_INTEGER
      END INTERFACE READ_TAG_ARRAY

#if defined(MPI) || defined(MPI_CHAIN)
      ! Initialize tag list (sets all internal states to 0, i.e. "unset").
      CALL TAG_LIST%INIT()
      ! Set the prefix for the reading
      IF (INSTANCES == 1) THEN
         PREFIX = trim("")
      ELSE
         PREFIX = trim("ML_INSTANCE_" // str(CURRENT_INSTANCE) // "/")
      ENDIF

      ! Set LMLONLY_FF to .FALSE. is later maybe change if ML_ISTART=2
      FF%LMLONLY = .FALSE.
      ! Allocate arrays necessary to read INCAR file.
      IF (ALLOCATED(ML_INPUT_HANDLE%EATOM_VASP)) THEN
         DEALLOCATE(ML_INPUT_HANDLE%EATOM_VASP)
      ENDIF
      ALLOCATE(ML_INPUT_HANDLE%EATOM_VASP(1:NTYP))
      ! check if INCAR file exists
      CALL OPEN_INCAR_IF_FOUND(IU5, LOPEN)
      ! Reading of the main varialbe ML_LMLFF is done in READER!!!
      ! If ML_LMLFF=.TRUE., the parameters related to MLFF are read in the
      ! INCAR file.
      IF(ML_LMLFF) THEN
         ! check here if IVAWPR was set
         CALL PROCESS_INCAR(LOPEN,             &
                            IU0,               &
                            IU5,               &
                            PREFIX// 'IWAVPR', &
                            IWAVPR_HELP,       &
                            IERR,              &
                            .FALSE.,           &
                            LCONTINUE=.TRUE.,  &
                            FOUNDNUMBER=N)
         IF (N == 0) THEN
            CALL PROCESS_INCAR(LOPEN,            &
                               IU0,              &
                               IU5,              &
                               'IWAVPR',         &
                               IWAVPR_HELP,      &
                               IERR,             &
                               .FALSE.,          &
                               LCONTINUE=.TRUE., &
                               FOUNDNUMBER=M)
            IF (M == 0) THEN
               IWAVPR = 11
            ENDIF
         ENDIF
         ! Read FF%CDOUB.
         CALL READ_TAG('ML_CDOUB', FF%CDOUB, 2.0_q)
         ! Read FF%CSIG that is a threshold for variance in the stored
         ! estimated errors used to determine the criteria when
         ! FF%ICRITERIA > 0.
         CALL READ_TAG('ML_CSIG', FF%CSIG, 4E-01_q)
         ! Read FF%CSLOPE that is a threshold for slope of the stored estimated
         ! errors used to determine the criteria when FF%ICRITERIA > 0.
         CALL READ_TAG('ML_CSLOPE', FF%CSLOPE, 2E-01_q)
         ! Read CTIFOR that corresponds to the criteria for the errors in
         ! forces estimated by the Bayesian error estimation method.
         CALL READ_TAG('ML_CTIFOR', FF%CTIFOR, 0.002_q)
         ! Read SCALE_CTIFOR which controls the selected number of local
         ! reference configurations
         CALL READ_TAG('ML_SCLC_CTIFOR',FF%SCALE_CTIFOR,1.0_q)
         ! Read LCTIFOR_USE_FROM_FILE, which decides to use previous
         ! information of CTIFOR from ML_AB file
         CALL READ_TAG('ML_LCTIFOR_USE_FROM_FILE',  &
                       FF%LCTIFOR_USE_FROM_FILE, &
                       .TRUE.)
         ! Read reference atomic energies, which should be energies of isolated
         ! atoms.
         CALL READ_TAG_ARRAY('ML_EATOM_REF',             &
                             ML_INPUT_HANDLE%EATOM_VASP, &
                             NTYP,                       &
                             0.0_q)
         ML_INPUT_HANDLE%EATOM_VASP(1:NTYP) = &
            ML_INPUT_HANDLE%EATOM_VASP(1:NTYP) / EUNIT
         ! Reading EPS_LOW
         CALL READ_TAG('ML_EPS_LOW', FF%ML_EPS_LOW, 1E-10_q)
         ! Reading EPS_REG
         CALL READ_TAG('ML_EPS_REG', FF%ML_EPS_REG, 1E-14_q)
         ! Read FF%IMAT_SPARS, which specifies the type of cutoff function.
         !    FF%IMAT_SPARS = 1 : Covariance matrix.
         !    FF%IMAT_SPARS = 2 : SOAP matrix.
         CALL READ_TAG('ML_IMAT_SPARS', FF%IMAT_SPARS, 1)
         ! Choose which algorithm is used for linear regression
         !    ML_IALGO_LINREG = 1 : Solve L2-norm (ridge regression)
         !    ML_IALGO_LINREG = 2 : Solve L1-norm with QR factorization
         !    ML_IALGO_LINREG = 3 : Solve L1-norm with truncated SVD
         !    ML_IALGO_LINREG = 4 : SVD + Tihonov regularization
         CALL READ_TAG('ML_IALGO_LINREG', FF%IALGO_LINREG, 1)
         ! Read FF%ISTART.
         !    ML_ISTART = 0 : Start the calculation without any previous data.
         !    ML_ISTART = 1 : Restart the calculation after reading the
         !                    previous database file.
         !    ML_ISTART = 2 : Only force field is used without learning.
         !    ML_ISTART = 3 : New force-field created only from ML_AB file.
#ifndef ML_LEARNING_NOT_ACTIVE
         CALL READ_TAG('ML_ISTART', FF%ISTART, 0)
#else
         CALL READ_TAG('ML_ISTART', FF%ISTART, 2)
#endif
#ifndef ML_LEARNING_NOT_ACTIVE
#else
         IF (FF%ISTART /= 2) THEN
            CALL vtutor%error("ERROR: Learning is disabled (source was &
                              &compiled with flag ML_LEARNING_NOT_ACTIVE). &
                              &Only ML_ISTART=2 is available, exiting...")
         ENDIF
#endif
         IF(FF%ISTART == 0 .OR. FF%ISTART == 1) THEN
            FF%LMLONLY = .FALSE.
         ELSE IF (FF%ISTART == 2) THEN
            FF%LMLONLY = .TRUE.
         ELSE IF (FF%ISTART == 3) THEN
            FF%LMLONLY = .FALSE.
         ENDIF
         ! Read FF%ISVD.
         CALL READ_TAG('ML_ISVD', FF%ISVD, 1)
         ! Read FF%IUPDATE_CRITERIA.
         CALL READ_TAG('ML_IUPDATE_CRITERIA', FF%IUPDATE_CRITERIA, 1)
         ! Read FF%LBASIS_DISCARD.
         !   If FF%LBASIS_DISCARD = .TRUE., the basis sets are thrown away if
         !   NB > FFM%MB.
         CALL READ_TAG('ML_LBASIS_DISCARD', FF%LBASIS_DISCARD, .FALSE.)
         ! Read LCCONF_DISCARD_FF.
         !   If FF%LCONF_DISCARD = .TRUE., the training sets are thrown away if
         !   NCONF > FF%MCONF.
         CALL READ_TAG('ML_LCONF_DISCARD', FF%LCONF_DISCARD, .FALSE.)
         ! Read FF%ICRITERIA.
         CALL READ_TAG('ML_ICRITERIA', FF%ICRITERIA, 1)
         ! Read FF%LMLMB.
         !   FF%LMLMB = .FALSE. : No many-body interactions. But currently,
         !   this program does not work when FF%LMLMB is set as .FALSE.
         !   FF%LMLMB = .TRUE.  : Execute force field generation with the
         !   many-body interactions.
         CALL READ_TAG('ML_LMLMB', FF%LMLMB, .TRUE.)
         ! Read FF%LERR.
         !   FF%LERR = .TRUE.  : Execute error estimations.
         !   FF%LERR = .FALSE. : Do not execute error estimations.
         CALL READ_TAG('ML_LERR', FF%LERR, .TRUE.)
         ! Read ML_INPUT_HANDLE%NDIM_SCALAPACK_FF
         CALL READ_TAG('ML_NDIM_SCALAPACK',            &
                       ML_INPUT_HANDLE%NDIM_SCALAPACK_FF, &
                       2)
         ! Read ML_INPUT_HANDLE%NWRITE_FF.
         CALL READ_TAG('ML_FF_NWRITE', ML_INPUT_HANDLE%NWRITE_FF, 2)
         ! Read FF%IWEIGHT
         CALL READ_TAG('ML_IWEIGHT', FF%IWEIGHT, 3)
         ! Tag to control wether ML_NMDINT is multiplied with a random
         ! number between 0 and 1
         ! If .TRUE. then NMDINT is varying between 1 and NMDINT
         CALL READ_TAG('ML_LNMDINT_RANDOM', FF%LNMDINT_RANDOM, .FALSE.)
         ! Read FF%LTEST.
         !   FF%LTEST = .TRUE.  : Execute ab initio calculations to get test
         !                        data.
         !   FF%LTEST = .FALSE. : Do not execute ab initio calculations to get
         !                        test data.
         CALL READ_TAG('ML_LTEST', FF%LTEST, .FALSE.)
         ! Read FF%LTRJ.
         !   FF%LTRJ = .TRUE.  : MD trajectory is always obtained using ab
         !                       initio forces and stress tensors.
         !   FF%LTRJ = .FALSE. : MD trajectory is obtained using generating
         !                       force field when the estimated error is small.
         !                       Otherwise, ab initio forces and stress tensors
         !                       are used.
         CALL READ_TAG('ML_LTRJ', FF%LTRJ, .FALSE.)
         ! READ FF%LTOTEN_SYSTEM
         !    FF%LTOTEN_SYSTEM = .TRUE. : Total energy of system instead
         !                                energy/atom is learned
         !    FF%LTOTEN_SYSTEM = .FALSE. : Energy/atom is learned
         CALL READ_TAG('ML_LTOTEN_SYSTEM', FF%LTOTEN_SYSTEM, .FALSE.)
         ! Read FF%MCONF that is the maximum number of configurations used as
         ! training data.
         CALL READ_TAG('ML_MCONF', FF%MCONF, 1500)
         ! Read FF%MCONF_NEW that is the maximum number of new configurations
         ! stored as candidates of training data temporally.
         CALL READ_TAG('ML_MCONF_NEW', FF%MCONF_NEW, 5)
         ! Read FF%MHIS that is the number of estimated errors stored in memory
         ! to judge the criteria.
         CALL READ_TAG('ML_MHIS', FF%MHIS, 10)
         ! Read FF%NMDINT that is the number of interval where machine do not
         ! collect any samples.
         IF (FF%ISTART == 3) THEN
            CALL READ_TAG('ML_NMDINT', FF%NMDINT, 1)
         ELSE
            CALL READ_TAG('ML_NMDINT', FF%NMDINT, 10)
         ENDIF
         ! Read FF%NTEST.
         ! If FF%LTEST = .TRUE., test data will be obtained at every NTEST MD
         ! step.
         IF(FF%LTEST) THEN
            CALL READ_TAG('ML_NTEST', FF%NTEST, 10)
         ENDIF
         ! Read FF%WTIFOR, FF%WTOTEN and FF%WTSIF
         ! in (eV Angst^-1), (eV atom^-1) and (kbar), respectively.
         IF(FF%IWEIGHT == 1) THEN
            CALL READ_TAG('ML_WTIFOR', FF%WTIFOR, 5E-02_q)
            CALL READ_TAG('ML_WTOTEN', FF%WTOTEN, 5E-03_q)
            CALL READ_TAG('ML_WTSIF', FF%WTSIF, 5.0_q)
         ELSE
            CALL READ_TAG('ML_WTIFOR', FF%WTIFOR, 1.0_q)
            CALL READ_TAG('ML_WTOTEN', FF%WTOTEN, 1.0_q)
            CALL READ_TAG('ML_WTSIF', FF%WTSIF, 1.0_q)
         ENDIF
         ! Read FF%CX (parameter x for criteria update).
         CALL READ_TAG('ML_CX', FF%CX, 0.0_q)
         ! If FF%LMLMB = .TRUE., read variables necessary for the many-body
         ! fitting.
         IF (FF%LMLMB) THEN
            FFM%LSIC = .TRUE.
            ! Read FFM%LAFILT2.
            CALL READ_TAG('ML_LAFILT2', FFM%LAFILT2, .TRUE.)
            IF (FFM%LAFILT2) THEN
               CALL READ_TAG('ML_IAFILT2', FFM%IAFILT2, 2)
               IF (FFM%IAFILT2 == 2) THEN
                  CALL READ_TAG('ML_AFILT2', FFM%AFILT2, 2E-03_q)
               ENDIF
            ENDIF
            ! Read ML_IBROADX (X=1 or 2).
            !  ! this variable determines the way to braoden the atomic distributions.
            ! ML_IBROADX=1 : No broadening.
            ! ML_IBROADX=2 : Element-independent user specified broadening.
            ! ML_IBROADX=3 : Element-dependent broadening
            ! Probably only ML_IBROADX=2 works!!!
            CALL READ_TAG('ML_IBROAD1', FFM%IBROAD1, 2)
            CALL READ_TAG('ML_IBROAD2', FFM%IBROAD2, 2)
            ! Read ML_ICUTX (X=1 or 2), which specifies the type of
            ! cutoff function.
            !   ML_ICUTX = 1 : Behler-Parinello type cutoff function.
            !   ML_ICUTX = 2 : Miwa-Ohno type cutoff function.
            CALL READ_TAG('ML_ICUT1', FFM%ICUT1, 1)
            CALL READ_TAG('ML_ICUT2', FFM%ICUT2, 1)
            ! this variable determines the method to calculate inverse of SOAP
            ! matrix.
            !   FFM%INVERSE_SOAP = 1 : Inverse matrix is calculated by using
            !                          eigen-values and eigen-vectors.
            !   FFM%INVERSE_SOAP = 2 : Inverse matrix is calculated by LU
            !                          factorization.
            CALL READ_TAG('ML_INVERSE_SOAP', FFM%INVERSE_SOAP, 2)
            ! Read FFM%IREG
            CALL READ_TAG('ML_IREG', FFM%IREG, 2)
            ! Read FFM%ISCALE_TOTEN
            CALL READ_TAG('ML_ISCALE_TOTEN', FFM%ISCALE_TOTEN, 2)
            ! Read FFM%LSUPERVEC
            CALL READ_TAG('ML_LSUPERVEC', FFM%LSUPERVEC, .TRUE.)
            ! Read FFM%LEATOM
            CALL READ_TAG('ML_LEATOM', FFM%LEATOM, .FALSE.)
            ! Read FFM%LHEAT
            CALL READ_TAG('ML_LHEAT', FFM%LHEAT, .FALSE.)
            ! Read ML_LMETRICX.
            CALL READ_TAG('ML_LMETRIC1', FFM%LMETRIC1, .FALSE.)
            CALL READ_TAG('ML_LMETRIC2', FFM%LMETRIC2, .FALSE.)
            ! Read FFM%LSPARSDES and related tags.
            CALL READ_TAG('ML_LSPARSDES', FFM%LSPARSDES, .FALSE.)
            CALL READ_TAG('ML_NRANK_SPARSDES', FFM%NRANK_SPARSDES, 5)
            CALL READ_TAG('ML_RDES_SPARSDES', FFM%RDES_SPARSDES, 0.5_q)
            ! Sparsification of descriptors is only allowed in ISTART = 1 mode,
            ! hence overwrite user setting if set incorrectly. Related tags can
            ! be left unchanged, they will only be used if LSPARSDES = .TRUE.
            IF (FF%ISTART /= 1 .AND. &
                TAG_LIST%GET_STATE('ML_LSPARSDES') /= TL_DEFAULT) THEN
               FFM%LSPARSDES = .FALSE.
               CALL TAG_LIST%SET_STATE('ML_LSPARSDES', TL_OVERRIDE) 
            END IF
            ! Read ML_LVARTRANX.
            CALL READ_TAG('ML_LVARTRAN1', FFM%LVARTRAN1, .FALSE.)
            CALL READ_TAG('ML_LVARTRAN2', FFM%LVARTRAN2, .FALSE.)
            ! Read ML_NMETRICX.
            IF (FFM%LMETRIC1) THEN
               CALL READ_TAG('ML_NMETRIC1', FFM%NMETRIC1, 6)
            ENDIF
            IF (FFM%LMETRIC2) THEN
               CALL READ_TAG('ML_NMETRIC2', FFM%NMETRIC2, 6)
            ENDIF
            ! Read ML_NVARTRANX.
            IF (FFM%LVARTRAN1) THEN
               CALL READ_TAG('ML_NVARTRAN1', FFM%NVARTRAN1, 6)
            ENDIF
            IF (FFM%LVARTRAN2) THEN
               CALL READ_TAG('ML_NVARTRAN2', FFM%NVARTRAN2, 6)
            ENDIF
            ! Read LWINDOWX_FFY (X=1, 2 or nothing, and Y=M or Q)
            !   ML_LWINDOWX_Y = .TRUE. : window function is used to smoothen
            !                               the expanded atomic distribution.
            !   ML_LWINDOWX_Y = .FLASE.: window function is not used.
            CALL READ_TAG('ML_LWINDOW1', FFM%LWINDOW1, .FALSE.)
            CALL READ_TAG('ML_LWINDOW2', FFM%LWINDOW2, .FALSE.)
            ! Read IWINDOWX_FFY (X=1, 2 or nothing, Y=M or Q) if LWINDOWX_FFY =
            ! .TRUE.
            IF (FFM%LWINDOW1) THEN
               CALL READ_TAG('ML_IWINDOW1', FFM%IWINDOW1, 6)
            ENDIF
            IF (FFM%LWINDOW2) THEN
               CALL READ_TAG('ML_IWINDOW2', FFM%IWINDOW2, 6)
            ENDIF
            ! Read FFM%LMAX2 that defines the maximum angular momentum number
            ! of spherical harmonics used to expand atomic distributions.
            IF (FFM%LAFILT2) THEN
               CALL READ_TAG('ML_LMAX2', FFM%LMAX2, 4)
            ELSE
               CALL READ_TAG('ML_LMAX2', FFM%LMAX2, 6)
            ENDIF
            ! Read FFM%LNORM1 that specifies the execution of normalization.
            CALL READ_TAG('ML_LNORM1', FFM%LNORM1, .TRUE.)
            ! Read FFM%LNORM2 that specifies the execution of normalization.
            CALL READ_TAG('ML_LNORM2', FFM%LNORM2, .TRUE.)
            ! Read FFM%MB that is the maximum number of basis sets used for
            ! describing many-body interactions.
            CALL READ_TAG('ML_MB', FFM%MB, 1500)
            ! Read ML_MSPLX (X=1 or 2) that is the number of radial grids
            ! for spline-interpolations.
            CALL READ_TAG('ML_MSPL1', FFM%MSPL1, 1000)
            CALL READ_TAG('ML_MSPL2', FFM%MSPL2, FFM%MSPL1)
            ! Read NHYP that corresponds to the polynomial parameters of SOAP
            ! 1-point kernel.
            IF (FFM%LSUPERVEC) THEN
               CALL READ_TAG('ML_NHYP', FFM%NHYP1, 4)
            ELSE
               CALL READ_TAG('ML_NHYP', FFM%NHYP1, 1)
            ENDIF
            ! Read NHYP that corresponds to the polynomial parameters of SOAP
            ! 2-point kernel.
            CALL READ_TAG('ML_NHYP2', FFM%NHYP2, 4)
            IF (FFM%LSUPERVEC) THEN
               FFM%NHYP2 = FFM%NHYP1
            ENDIF
            ! Read ML_NRX (X=1 or 2) that is the number of radial grids
            ! used to execute radial integrations to compute the many-body
            ! descriptors.
            CALL READ_TAG('ML_NR1', FFM%NR1, 1000)
            CALL READ_TAG('ML_NR2', FFM%NR2, FFM%NR1)
            ! Read FFM%W1 that
            CALL READ_TAG('ML_W1', FFM%W1, 0.1_q)
            ! Read FFM%W2 that
            ! We disable the reading of ML_W2 for the moment because
            ! everything is determined by ML_W1
            FFM%W2 = 1.0_q - FFM%W1
            !CALL READ_TAG('ML_W2', FFM%W2, 1.0_q - FFM%W1)
            ! Read ML_RCUTX (X=1 or 2) that is the cutoff radius used for
            ! calculating many-body descriptors.
            CALL READ_TAG('ML_RCUT1', FFM%RCUT1, 5.0_q)
            CALL READ_TAG('ML_RCUT2', FFM%RCUT2, FFM%RCUT1)
            IF (FFM%W1 == 0.0_q) THEN
               FFM%RCUT1 = FFM%RCUT2
            ENDIF
            IF (FFM%W2 == 0.0_q) THEN
               FFM%RCUT2 = FFM%RCUT1
            ENDIF
            ! Read ML_RMETRICX (X=1 or 2).
            IF (FFM%LMETRIC1) THEN
               CALL READ_TAG('ML_RMETRIC1', FFM%RMETRIC1, 1.0_q)
            ENDIF
            IF (FFM%LMETRIC2) THEN
               CALL READ_TAG('ML_RMETRIC2', FFM%RMETRIC2, 1.0_q)
            ENDIF
            ! Read FFM%SIGV0
            CALL READ_TAG('ML_SIGV0', FFM%SIGV0, 1.0_q)
            ! Read FFM%SIGW0
            CALL READ_TAG('ML_SIGW0', FFM%SIGW0, 1.0_q)
            ! Read ML_SIONX (X=1 or 2) that specifies the width of
            ! Gaussian functions used for broadening the atomic distributions.
            IF (FFM%IBROAD1 /= 1) THEN
               CALL READ_TAG('ML_SION1', FFM%SION1, 5.0E-01_q)
            ENDIF
            IF (FFM%IBROAD2 /= 1) THEN
               CALL READ_TAG('ML_SION2', FFM%SION2, FFM%SION1)
            ENDIF
            ! Read FFM%LCOUPLE that specifies whether the coupling parameter is
            ! introduced for the calculation of the chemical potential.
            CALL READ_TAG('ML_LCOUPLE', FFM%LCOUPLE, .FALSE.)
            IF (FFM%LCOUPLE) THEN
               CALL READ_TAG('ML_NATOM_COUPLED', FFM%NATOM_COUPLED, 0)
               IF (ALLOCATED(FFM%ICOUPLE)) THEN
                  DEALLOCATE(FFM%ICOUPLE)
               ENDIF
               ALLOCATE(FFM%ICOUPLE(1:FFM%NATOM_COUPLED))
               CALL READ_TAG_ARRAY('ML_ICOUPLE',      &
                                   FFM%ICOUPLE,       &
                                   FFM%NATOM_COUPLED, &
                                   0)
               IF (N > 0 .AND. (N /= FFM%NATOM_COUPLED)) THEN
                  CALL vtutor%error("ML_INTERFACE: ERROR: Number of atoms for &
                                    &coupling .NE. ML_NATOM_COUPLED.")
               ENDIF
               CALL READ_TAG('ML_RCOUPLE', FFM%RCOUPLE, 1.0_q)
            ELSE
               FFM%NATOM_COUPLED = 0
               CALL TAG_LIST%SET_STATE('ML_NATOM_COUPLED', TL_DEFAULT)
               IF(ALLOCATED(FFM%ICOUPLE)) THEN
                  DEALLOCATE(FFM%ICOUPLE)
               ENDIF
               ALLOCATE(FFM%ICOUPLE(1:FFM%NATOM_COUPLED))
               FFM%ICOUPLE = 0
               CALL TAG_LIST%SET_STATE('ML_ICOUPLE', TL_DEFAULT)
               FFM%RCOUPLE = 1.0_q
               CALL TAG_LIST%SET_STATE('ML_RCOUPLE', TL_DEFAULT)
            ENDIF
            ! Read ML_MRBX (X=1 or 2) that is the number of radial basis
            ! sets used to expand the atomic distributions.
            ! Automatically scale values of FFM%MRB2 if no value is supplied
            CALL READ_TAG('ML_MRB1', &
                          FFM%MRB1,  &
                          8)
            CALL READ_TAG('ML_MRB2', &
                          FFM%MRB2,  &
                          FFM%MRB1)
            ! Use random seed if available
            CALL READ_TAG_ARRAY('ML_RANDOM_SEED', FF%SEED, 3, 0)
            ! Decide whether configurations not giving basis functions will be
            ! thrown away
            IF (FF%ISTART == 3) THEN
               CALL READ_TAG('ML_LDISCARD_STRUCTURES_NOT_GIVING_BASIS', &
                             FF%LDISCARD_STRUCTURES_NOT_GIVING_BASIS,   &
                             .FALSE.)
            ENDIF
            ! Blocking factor for printing out ML_FFN file
            CALL READ_TAG('ML_NBLOCK', FF%NBLOCK_FFN, 50)
            ! Set some important stuff
            IF (IERR == 3) THEN
               CALL DATE_AND_TIME(VALUES=DTVALUES)
               ! with this choice we'r sure not to exceed
               ! SEED1_MAX (maximal value we can get here is 893581000)
               IDUMLONG = DTVALUES(3) * 24 * 60 * 60 * 300 &
                        + DTVALUES(5) * 60 * 60 * 1000 &
                        + DTVALUES(6) * 60 * 1000 &
                        + DTVALUES(7) * 1000 &
                        + DTVALUES(8)
               FF%SEED(1) = MOD(IDUMLONG, SEED1_MAX)
            ELSE IF ((FF%SEED(1) < 0) .OR. (FF%SEED(1) > SEED1_MAX) .OR. &
                     (FF%SEED(2) < 0) .OR. (FF%SEED(2) > SEED2_MAX) .OR. &
                     (FF%SEED(3) < 0)) THEN
               CALL vtutor%write(isError, RandomSeed, argument(ival = FF%SEED))     
            ENDIF
! Set some important variables
            FF%NHIS = 0
            FFM%LMAX1 = 0
            FFM%MB_INPUT = FFM%MB
            FFM%MB_ALLOCATE = FFM%MB
            FFM%NSPL1 = FFM%MSPL1
            FFM%NSPL2 = FFM%MSPL2

! Unit conversions
            IF (FF%IWEIGHT == 1) THEN
               FF%WTIFOR = FF%WTIFOR / FUNIT
               FF%WTOTEN = FF%WTOTEN / EUNIT
               FF%WTSIF = FF%WTSIF / SUNIT
            ENDIF
            FF%CTIFOR = FF%CTIFOR / FUNIT
            FFM%RCUT1 = FFM%RCUT1 / AUTOA
            FFM%RCUT2 = FFM%RCUT2 / AUTOA
            FFM%SION1 = FFM%SION1 / AUTOA
            FFM%SION2 = FFM%SION2 / AUTOA
            IF (FFM%LMETRIC1) THEN
               FFM%RMETRIC1 = FFM%RMETRIC1 / AUTOA
            ENDIF
            IF (FFM%LMETRIC2) THEN
               FFM%RMETRIC2 = FFM%RMETRIC2 / AUTOA
            ENDIF
         ELSE
            FFM%IAFILT2 = 0
            FFM%INVERSE_SOAP = 0
            FFM%IREG = 0
            FFM%LMAX1 = 0
            FFM%LMAX2 = 0
            FFM%MB = 0
            FFM%MB_INPUT = 0
            FFM%MB_ALLOCATE = 0
            FFM%MSPL1 = 0
            FFM%MSPL2 = 0
            FFM%NHYP1 = 0
            FFM%NHYP2 = 0
            FFM%NR1 = 0
            FFM%NR2 = 0
            FFM%NSPL1 = 0
            FFM%NSPL2 = 0
            FFM%MRB1 = 0
            FFM%MRB2 = 0
            FFM%NMETRIC1 = 0
            FFM%NMETRIC2 = 0
            FFM%NVARTRAN1 = 0
            FFM%NVARTRAN2 = 0
            FFM%RCUT1 = 0.0_q
            FFM%RCUT2 = 0.0_q
            FFM%RMETRIC1 = 0.0_q
            FFM%RMETRIC2 = 0.0_q
            FFM%SIGV0 = 0.0_q
            FFM%SIGW0 = 0.0_q
            FFM%SION1 = 0.0_q
            FFM%SION2 = 0.0_q
            FFM%W1 = 0.0_q
            FFM%W2 = 0.0_q
         ENDIF
      ENDIF
#endif

      CONTAINS

#define TYPE_REAL
#include "ml_interface_read_tag.inc"
#undef TYPE_REAL

#define TYPE_INTEGER
#include "ml_interface_read_tag.inc"
#undef TYPE_INTEGER

#define TYPE_LOGICAL
#include "ml_interface_read_tag.inc"
#undef TYPE_LOGICAL

#define TYPE_REAL
#include "ml_interface_read_tag_array.inc"
#undef TYPE_REAL

#define TYPE_INTEGER
#include "ml_interface_read_tag_array.inc"
#undef TYPE_INTEGER

   END SUBROUTINE MACHINE_LEARNING_READER


#if defined(MPI) || defined(MPI_CHAIN)
!=======================================================================
!
! Routines important for the handling of global machine learning handles.
!
!=======================================================================

     SUBROUTINE MACHINE_LEARNING_HANDLE_ALLOCATE(NUM_INSTANCES)
       USE ML_FF_STRUCT, ONLY: ML_SUPER_HANDLE_MAIN
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: NUM_INSTANCES
       IF (ALLOCATED(ML_SUPER_HANDLE_MAIN)) THEN
          DEALLOCATE(ML_SUPER_HANDLE_MAIN)
       ENDIF
       ALLOCATE(ML_SUPER_HANDLE_MAIN(NUM_INSTANCES))
     END SUBROUTINE MACHINE_LEARNING_HANDLE_ALLOCATE

     SUBROUTINE MACHINE_LEARNING_HANDLE_DEALLOCATE
       USE ML_FF_STRUCT, ONLY: ML_SUPER_HANDLE_MAIN
       IMPLICIT NONE
       IF (ALLOCATED(ML_SUPER_HANDLE_MAIN)) THEN
          DEALLOCATE(ML_SUPER_HANDLE_MAIN)
       ENDIF
     END SUBROUTINE MACHINE_LEARNING_HANDLE_DEALLOCATE

!=======================================================================
!
! the following routine is called once, before initialization is called
! to set a global COMM_WORLD communicator. This communicator should be
! only used very sparsely for emergency aborts and etc.
!
!=======================================================================
    SUBROUTINE MACHINE_LEARNING_SET_COMM_WORLD_GLOBAL(VASP_COMMUNICATOR)
      USE mpimy, ONLY : communic
      USE ML_FF_STRUCT, ONLY : ML_MPI_PAR, ML_COMM_WORLD_GLOBAL, ML_IO_WRITE
      IMPLICIT NONE
      TYPE (communic)   :: VASP_COMMUNICATOR
      ML_COMM_WORLD_GLOBAL%MPI_COMM = VASP_COMMUNICATOR%MPI_COMM
      ML_COMM_WORLD_GLOBAL%NCPU     = VASP_COMMUNICATOR%NCPU
      ML_COMM_WORLD_GLOBAL%NODE_ME     = VASP_COMMUNICATOR%NODE_ME
      ! set here which rank is allowed to write (currently only the first rank)
      IF (ML_COMM_WORLD_GLOBAL%NODE_ME.EQ.1) THEN
         ML_IO_WRITE = .TRUE.
      ELSE
         ML_IO_WRITE = .FALSE.
      ENDIF
    END SUBROUTINE MACHINE_LEARNING_SET_COMM_WORLD_GLOBAL

!=======================================================================
!
! the following routine is called once, after VASP has read
! the ionic positions + the Bravais lattice
! at this point the machine learning code is supposed to perform
! the basic setup
!
!=======================================================================

    SUBROUTINE MACHINE_LEARNING_INIT(ML_SUPER_HANDLE,DIR_APP,DIR_LEN,DYN, INFO, LATT_CUR, T_INFO,VASP_COMM)
      USE base
      USE lattice
      USE poscar
      USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE
      USE ML_FF_CONSTANT, ONLY: AUTOA, PI, MUNIT, TUNIT
      IMPLICIT NONE
      TYPE (ML_SUPER_TYPE) :: ML_SUPER_HANDLE
      TYPE (dynamics)    DYN
      TYPE (info_struct) INFO
      TYPE (latt)        LATT_CUR             ! Bravais lattice information
      TYPE (type_info)   T_INFO               ! type and position information
      CHARACTER(LEN=10), INTENT(IN) ::  DIR_APP
      INTEGER, INTENT(IN)           ::  DIR_LEN
      INTEGER, INTENT(IN)           ::  VASP_COMM

! One has to provide proper structure name in the POSCAR file otherwise the structure cannot be identified
! when reading the ML_AB file
      IF (T_INFO%SZNAM2.EQ."") THEN
         CALL vtutor%error("ML_INTERFACE: ERROR: No structure name was provided in the POSCAR file, exiting...")
      ENDIF
! Stop code if NSW>1 and ML_ISTART=3
! Give some variables from VASP to MLFF and allocate some necessary arrays
      ML_SUPER_HANDLE%FF%NSW = DYN%NSW 
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%A = LATT_CUR%A/AUTOA
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%B = LATT_CUR%B*AUTOA*2.0_q*PI
!      ML_SUPER_HANDLE%ML_INPUT_HANDLE%A = LATT_CUR%A
!      ML_SUPER_HANDLE%ML_INPUT_HANDLE%B = LATT_CUR%B
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%DIR_APP = DIR_APP
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%DIR_LEN = DIR_LEN
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%LABORT = INFO%LABORT
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%COMM_VASP = VASP_COMM
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS = T_INFO%NIONS
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%SZNAM2 = TRIM(ADJUSTL(T_INFO%SZNAM2))
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%OMEGA = LATT_CUR%OMEGA/AUTOA**3
!      ML_SUPER_HANDLE%ML_INPUT_HANDLE%OMEGA = LATT_CUR%OMEGA
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%POTIM = DYN%POTIM/TUNIT
!      ML_SUPER_HANDLE%ML_INPUT_HANDLE%POTIM = DYN%POTIM
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP = T_INFO%NTYP
      ALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP(1:ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP))
      ALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%POMASS(1:ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP))
      ALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TYPE(1:ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP))
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP = T_INFO%NITYP
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%POMASS = T_INFO%POMASS/MUNIT
!      ML_SUPER_HANDLE%ML_INPUT_HANDLE%POMASS = T_INFO%POMASS
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%TYPE = T_INFO%TYPE
! Set TIFOR_ML,TSIF_ML.
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TIFOR_ML)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TIFOR_ML)
      ENDIF
      ALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TIFOR_ML(1:3,1:ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS))
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TSIF_ML)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TSIF_ML)
      ENDIF
      ALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TSIF_ML(1:3,1:3))
! Set VEL_ML.
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML)
      ENDIF
      ALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML(1:3,1:ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS))
! Send variables to machine-learning library
      CALL MACHINE_LEARNING_INIT_LIB (ML_SUPER_HANDLE)
! We need to set this if ISTART=3
      IF (ML_SUPER_HANDLE%FF%ISTART.EQ.3) THEN
         LDO_AB_INITIO=.FALSE.
         LDO_AB_INITIO_NEW=.FALSE.
      ENDIF
    END SUBROUTINE MACHINE_LEARNING_INIT

!=======================================================================
!
! the following routine whenever a new set of positions
! is available to VASP the ML code, might be already able
! to calculate forces etc. and in that case returns the energy
! forces and the stress tensor, as well as an error estimate
! if LDO_AB_INITIO_NEW is set to be true, VASP will do the 
! calculations for the present positions, otherwise
! we will use the predicted data
!
!=======================================================================

    SUBROUTINE ML_TO_VASP_MACHINE_LEARNING(ML_SUPER_HANDLE, DYN, GRIDC, IO, LATT_CUR, LREMOVE_DRIFT, LWRITE_FORCE, &
               NIOND, NSTEP, SYMM, T_INFO, TOTEN, TIFOR, TSIF, VASP_COMM)
      USE base
      USE classicfields
      USE constant
      USE ini
      USE lattice
      USE poscar
      USE mgrid
      USE mpimy
      USE msymmetry
      USE poscar
      USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE
      USE ML_FF_CONSTANT, ONLY: MUNIT, TUNIT
      IMPLICIT NONE
!Input and output variables
      TYPE (ML_SUPER_TYPE)          :: ML_SUPER_HANDLE
      TYPE (dynamics)               :: DYN
      TYPE (grid_3d)                :: GRIDC             ! grid for potentials/charge
      TYPE (in_struct)              :: IO
      TYPE (latt), INTENT(IN)       :: LATT_CUR          ! Bravais lattice information
      TYPE(symmetry)                :: SYMM
      TYPE (type_info), INTENT(IN)  :: T_INFO            ! type and position information
      LOGICAL, INTENT(IN)           :: LWRITE_FORCE      ! If LWRITE_FORCE=.TRUE., write forces in OUTCAR file.
      INTEGER, INTENT(IN)           :: NIOND             ! The number of ions.
      INTEGER, INTENT(IN)           :: NSTEP             ! MD step
      REAL (q), INTENT(OUT)         :: TOTEN             ! total energy
      REAL (q), INTENT(OUT)         :: TIFOR(:,:)        !(1:3,1:NIOND) ! forces on all atoms
      REAL (q), INTENT(OUT)         :: TSIF(:,:)         !(1:3,1:3)      ! stress tensor
      INTEGER, INTENT(IN)           :: VASP_COMM
! Local variables related to symmetry.
      INTEGER  ISYMOP,NROTK
      LOGICAL LREMOVE_DRIFT
      COMMON /SYMM/   ISYMOP(3,3,48),NROTK
      REAL (q)         :: TOTEN_TMP             ! total energy
      REAL (q)         :: TIFOR_TMP(1:3,1:NIOND)        !(1:3,1:NIOND) ! forces on all atoms
      REAL (q)         :: TSIF_TMP(1:3,1:3)         !(1:3,1:3)      ! stress tensor
      LOGICAL :: LDO_AB_INITIO_NEW_TMP
! Local variables related to bias forces.
      REAL(q)   :: EV2KB
      REAL(q)   :: FAC
      INTEGER   :: I
      INTEGER   :: INIONS
      INTEGER   :: IXYZ
      INTEGER   :: JXYZ
      INTEGER   :: NOFFS
      INTEGER   :: NI
      INTEGER   :: NT
      REAL(q)   :: OFIELD_E                        ! energy contribution from order field
      REAL(q)   :: OFIELD_FOR(1:3,T_INFO%NIONS)    ! forces from order field
      REAL(q)   :: POSION_LOC(1:3,T_INFO%NIONS)
      REAL(q)   :: PRESS
      REAL(q)   :: pressure
      REAL(q)   :: TMP(1:6)
      REAL(q)   :: VEL(1:3)
      REAL(q)   :: VTMP(1:3)

! Do nothing and return if this flag is set to true (for instance ML_ISTART=3)
      IF (ML_SUPER_HANDLE%SKIP_ML_TO_VASP) THEN
         TOTEN=0.0_q
         TIFOR=0.0_q
         TSIF=0.0_q
         RETURN
      ENDIF

! Constants.
      EV2KB = EVTOJ*1E+22_q/LATT_CUR%OMEGA
      FAC  = AMTOKG* &
             1E5_q* &
             1E5_q* &
             1E30_q* &
             1E-8_q
! Give some variables from VASP to MLFF
      ML_SUPER_HANDLE%FF%NSW = DYN%NSW
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%A = LATT_CUR%A/AUTOA
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%B = LATT_CUR%B*AUTOA*2.0_q*PI
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%COMM_VASP = VASP_COMM
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%NIONS = T_INFO%NIONS
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP = T_INFO%NITYP
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%NSTEP = NSTEP
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%NTYP = T_INFO%NTYP
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%POMASS = T_INFO%POMASS/MUNIT
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%POTIM = DYN%POTIM/TUNIT
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%SZNAM2 = TRIM(ADJUSTL(T_INFO%SZNAM2))
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%TYPE = T_INFO%TYPE
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%OMEGA = LATT_CUR%OMEGA/AUTOA**3
! Calculate positions in cartesian coordinates
! Set velocity. At NSTEP>1, velocity at the next step is estimated by adding 
! current accerelation to current velocity.
      DO INIONS=1, T_INFO%NIONS
         POSION_LOC(1:3,INIONS)=T_INFO%POSION(1,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,1) + &
                                T_INFO%POSION(2,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,2) + &
                                T_INFO%POSION(3,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,3)
         ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML(1:3,INIONS)=DYN%VEL(1,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,1) + &
                                                            DYN%VEL(2,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,2) + &
                                                            DYN%VEL(3,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,3)
         ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML(1:3,INIONS)=ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML(1:3,INIONS) / &
                                                           ML_SUPER_HANDLE%ML_INPUT_HANDLE%POTIM
      ENDDO
! Predict energy, force and stress tensor by using machine-learning force field.
      IF (ML_SUPER_HANDLE%FF%LSECOND_TIME_ML_TO_VASP) THEN
         CALL ML_TO_VASP_MACHINE_LEARNING_LIB (ML_SUPER_HANDLE,ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:,1:),ML_SUPER_HANDLE%ML_INPUT_HANDLE%B(1:,1:), &
              ML_SUPER_HANDLE%ML_INPUT_HANDLE%POTIM,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NSTEP, &
              ML_SUPER_HANDLE%ML_INPUT_HANDLE%OMEGA,POSION_LOC,ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML, &
              ML_SUPER_HANDLE%ML_INPUT_HANDLE%SZNAM2,TOTEN_TMP,TIFOR_TMP,TSIF_TMP,LDO_AB_INITIO_NEW_TMP)
      ELSE
         CALL ML_TO_VASP_MACHINE_LEARNING_LIB (ML_SUPER_HANDLE,ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:,1:),ML_SUPER_HANDLE%ML_INPUT_HANDLE%B(1:,1:), &
              ML_SUPER_HANDLE%ML_INPUT_HANDLE%POTIM,ML_SUPER_HANDLE%ML_INPUT_HANDLE%NSTEP, &
              ML_SUPER_HANDLE%ML_INPUT_HANDLE%OMEGA,POSION_LOC,ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML, &
              ML_SUPER_HANDLE%ML_INPUT_HANDLE%SZNAM2,TOTEN,TIFOR,TSIF,LDO_AB_INITIO_NEW)
! Add the external order field, if necessary
         CALL START_TIMING("G")
         CALL OFIELD(GRIDC%COMM, IO, LATT_CUR, DYN, T_INFO, OFIELD_FOR, OFIELD_E)
         TIFOR=TIFOR+OFIELD_FOR
         TOTEN=TOTEN+OFIELD_E
         CALL STOP_TIMING("G",IO%IU6,'OFIELD')
! Remove drift.
         IF (DYN%IBRION/=0) THEN
! remove drift from the forces
            IF (LREMOVE_DRIFT) CALL SYMVEC(T_INFO%NIONS,TIFOR)
         ENDIF
! Output to OUTCAR file.
         IF (IO%IU6>=0) THEN
            PRESS=(TSIF(1,1)+TSIF(2,2)+TSIF(3,3))/3._q- &
                  DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA
            WRITE(IO%IU6,1)       (TSIF(IXYZ,IXYZ),IXYZ=1,3),      TSIF (1,2),     TSIF(2,3),     TSIF (3,1), &
                            (EV2KB*TSIF(IXYZ,IXYZ),IXYZ=1,3),EV2KB*TSIF(1,2),EV2KB*TSIF(2,3),EV2KB*TSIF(3,1), &
                            PRESS*EV2KB,DYN%PSTRESS
            TMP = 0._q
            NOFFS = 0
            DO NT=1,T_INFO%NTYP
               DO NI=1+NOFFS,T_INFO%NITYP(NT)+NOFFS
                  VEL(1) = DYN%VEL(1,NI)/DYN%POTIM
                  VEL(2) = DYN%VEL(2,NI)/DYN%POTIM
                  VEL(3) = DYN%VEL(3,NI)/DYN%POTIM
                  CALL  DIRKAR( 1, VEL, LATT_CUR%A )
                  TMP(1) = TMP(1) + VEL(1)*VEL(2) * T_INFO%POMASS(NT)
                  TMP(2) = TMP(2) + VEL(2)*VEL(3) * T_INFO%POMASS(NT)
                  TMP(3) = TMP(3) + VEL(3)*VEL(1) * T_INFO%POMASS(NT)
                  TMP(4) = TMP(4) + VEL(1)*VEL(1) * T_INFO%POMASS(NT)
                  TMP(5) = TMP(5) + VEL(2)*VEL(2) * T_INFO%POMASS(NT)
                  TMP(6) = TMP(6) + VEL(3)*VEL(3) * T_INFO%POMASS(NT)
               ENDDO
               NOFFS = NOFFS + T_INFO%NITYP(NT)
            ENDDO
            TMP = TMP/LATT_CUR%OMEGA*FAC         ! NOW TMP IS IN KB
            WRITE(IO%IU6,'(''  kinetic pressure (ideal gas correction) = '',F9.2,'' kB'')') (TMP(4)+TMP(5)+TMP(6))/3
            pressure=(TSIF(1,1)+TSIF(2,2)+TSIF(3,3))*EV2KB
            pressure=pressure + (TMP(4)+TMP(5)+TMP(6))
            pressure=pressure/3.
            WRITE(IO%IU6,'(''  total pressure  = '',F9.2,'' kB'')') pressure
            WRITE(IO%IU6,'(''  Total+kin. '',F9.3,5F12.3)') TSIF(1,1)*EV2KB + TMP(4),  &
                                                            TSIF(2,2)*EV2KB + TMP(5),  &
                                                            TSIF(3,3)*EV2KB + TMP(6),  &
                                                            TSIF(1,2)*EV2KB + TMP(1),  &
                                                            TSIF(2,3)*EV2KB + TMP(2),  &
                                                            TSIF(3,1)*EV2KB + TMP(3)
            WRITE(IO%IU6,2) LATT_CUR%OMEGA, &
               ((LATT_CUR%A(IXYZ,JXYZ),IXYZ=1,3),(LATT_CUR%B(IXYZ,JXYZ),IXYZ=1,3),JXYZ=1,3), &
               (LATT_CUR%ANORM(IXYZ),IXYZ=1,3),(LATT_CUR%BNORM(IXYZ),IXYZ=1,3)
            IF(LWRITE_FORCE) THEN
               WRITE(IO%IU6,3)
               DO INIONS=1, T_INFO%NIONS
                  VTMP(1)=T_INFO%POSION(1,INIONS)
                  VTMP(2)=T_INFO%POSION(2,INIONS)
                  VTMP(3)=T_INFO%POSION(3,INIONS)
                  CALL  DIRKAR(1,VTMP(1),LATT_CUR%A(1,1))
                  WRITE(IO%IU6,4) (VTMP(I),I=1,3),(TIFOR(I,INIONS),I=1,3)
               ENDDO
               VTMP(1:3)=0.0_q
               DO INIONS=1, T_INFO%NIONS
                  DO IXYZ=1, 3
                     VTMP(IXYZ)=VTMP(IXYZ)+TIFOR(IXYZ,INIONS)
                  ENDDO
               ENDDO
               WRITE(IO%IU6,5) VTMP
            ENDIF
            WRITE(IO%IU6,6)
            WRITE(IO%IU6,7) '  ML FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)',TOTEN,TOTEN,TOTEN
         ENDIF
! Add previously subtracted pressure contribution to total energy egain
         IF (DYN%PSTRESS/=0) THEN
            TOTEN=TOTEN+DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA
            IF (IO%IU6>=0) WRITE(IO%IU6,8) TOTEN,DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA
8           FORMAT ('  enthalpy is ML TOTEN    = ',F18.8,' eV   P V=',F18.8/)
         ENDIF
! Format
1        FORMAT(/'  ML FORCE on cell =-STRESS in cart. coord. units (eV/cell)'/ &
                 '  Direction',4X,'XX', 10X,'YY', 10X,'ZZ', 10X,'XY', 10X,'YZ', 10X,'ZX'/ &
                 '  --------------------------------------------------------------------------------------'/ &
                 '  Total: ',6F12.5/ &
                 '  in kB  ',6F12.5/, &
                 '  external pressure = ',F11.2,' kB', &
                 '  Pullay stress = ',F11.2,' kB'/)
2        FORMAT( '  volume of cell :  ',F10.2/ &
                 '      direct lattice vectors',17X,'reciprocal lattice vectors'/ &
                 3(2(3X,3F13.9)/)/ &
                 '  length of vectors'/ &
                  (2(3X,3F13.9)/) /)
3        FORMAT( '  POSITION    ',35X,'TOTAL-FORCE (eV/Angst) (ML)'/ &
                 ' -----------------------------------------------------------------------------------')
4        FORMAT((3F13.5,3X,3F14.6))
5        FORMAT( ' ----------------------------------------------', &
                 '-------------------------------------',/ &
                 '    total drift:      ',20X,3F14.6)
6        FORMAT (5X, //, &
                '----------------------------------------------------', &
                '----------------------------------------------------'//)
7        FORMAT(/ &
                A/ &
                '  ---------------------------------------------------'/ &
                '  free  energy ML TOTEN  = ',F18.8,' eV'// &
                '  ML energy  without entropy=',F18.8, &
                '  ML energy(sigma->0) =',F18.8/)
      ENDIF
    END SUBROUTINE ML_TO_VASP_MACHINE_LEARNING

#ifndef ML_LEARNING_NOT_ACTIVE
!=======================================================================
!
! the following routine passes the energies and forces
! calculated by VASP to the ML code
! note that if LDO_AB_INITIO is false, these are just the
! very same forces and energies as calculated by the ML code
! Not TOTEN, TIFOR, TSIF but TOTEN_ML, TIFOR_ML, TSIF_ML must be
! provided to machine-learning library.
!
!=======================================================================

    SUBROUTINE VASP_TO_ML_MACHINE_LEARNING( ML_SUPER_HANDLE, DIR_APP,DIR_LEN,INFO, LATT_CUR, NSTEP, T_INFO, DYN, IO, GRIDC, IONS, TOTEN, TIFOR, TSIF )
      USE constant
      USE base
      USE mgrid
      USE ini
      USE lattice
      USE poscar
      USE mpimy
      USE classicfields
      USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE
      USE ML_FF_CONSTANT, ONLY: MUNIT, TUNIT, EUNIT, FUNIT, SUNIT
      IMPLICIT NONE
! Input and output variables
      TYPE (ML_SUPER_TYPE)          :: ML_SUPER_HANDLE
      CHARACTER(LEN=10), INTENT(IN) :: DIR_APP
      INTEGER          , INTENT(IN) :: DIR_LEN
      TYPE (info_struct)            :: INFO                      ! Information
      TYPE (latt)      , INTENT(IN) :: LATT_CUR    ! Bravais lattice information
      INTEGER          , INTENT(IN) :: NSTEP       ! MD step
      TYPE (type_info) , INTENT(IN) :: T_INFO      ! type and position information
      TYPE (dynamics)  , INTENT(IN) :: DYN
      TYPE (in_struct) , INTENT(IN) :: IO
      TYPE (grid_3d)   , INTENT(IN) :: GRIDC
      INTEGER                       :: IONS
      REAL(q)                       :: TOTEN
      REAL(q)                       :: TIFOR(:,:) !(1:3,1:IONS)
      REAL(q)                       :: TSIF(:,:) !(1:3,1:3)
! LOCAL variables
      REAL(q)                       :: OFIELD_E                        ! energy contribution from orderfield
      REAL(q)                       :: OFIELD_FOR (1:3,T_INFO%NIONS)    ! forces from order field
      REAL(q)                       :: FAKT
      INTEGER                       :: INIONS
      LOGICAL                       :: LDO_AB_INITIO_HELP
      REAL(q)                       :: POSION_LOC(1:3,T_INFO%NIONS)
      REAL(q)                       :: TOTEN_TMP
      REAL(q)                       :: TIFOR_TMP(1:3,1:T_INFO%NIONS)
      REAL(q)                       :: TSIF_TMP(1:3,1:3)
      REAL(q)                       :: VEL_TMP(1:3,1:T_INFO%NIONS)

! Do nothing and return if this flag is set to true (for instance ML_ISTART=3)
      IF (ML_SUPER_HANDLE%SKIP_VASP_TO_ML) RETURN

! Initialize helping arrays
      VEL_TMP=0.0_q

      ! Conversions needed for input to INIT_LIB
      FAKT=EVTOJ*1.0E+22_q/LATT_CUR%OMEGA !/AUTOA**3

      ! subtract stress contribution of energy that was added in main.F
      IF (DYN%PSTRESS/=0) THEN
         TOTEN=TOTEN-DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA
      ENDIF

      ! subtract contribution of OFIELD that was added in force.F
      CALL START_TIMING("G")
      CALL OFIELD(GRIDC%COMM, IO, LATT_CUR, DYN, T_INFO, OFIELD_FOR, OFIELD_E) 
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%TOTEN_ML=TOTEN-OFIELD_E
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%TIFOR_ML(1:3,1:IONS)=TIFOR(1:3,1:IONS)-OFIELD_FOR(1:3,1:IONS)
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%TSIF_ML(1:3,1:3)=TSIF(1:3,1:3)
      CALL STOP_TIMING("G",IO%IU6,'OFIELD')

! We need to do unit conversions here
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%TOTEN_ML=ML_SUPER_HANDLE%ML_INPUT_HANDLE%TOTEN_ML/EUNIT
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%TIFOR_ML(1:3,1:IONS)=ML_SUPER_HANDLE%ML_INPUT_HANDLE%TIFOR_ML(1:3,1:IONS)/FUNIT
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%TSIF_ML(1:3,1:3)=ML_SUPER_HANDLE%ML_INPUT_HANDLE%TSIF_ML(1:3,1:3)*FAKT/SUNIT
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,1:3) = LATT_CUR%A(1:3,1:3)/AUTOA
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%B(1:3,1:3) = LATT_CUR%B(1:3,1:3)*AUTOA*2.0_q*PI
      ML_SUPER_HANDLE%ML_INPUT_HANDLE%OMEGA = LATT_CUR%OMEGA/AUTOA**3
      DO INIONS=1, T_INFO%NIONS
         POSION_LOC(1:3,INIONS)=T_INFO%POSION(1,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,1) + &
                                T_INFO%POSION(2,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,2) + &
                                T_INFO%POSION(3,INIONS)*ML_SUPER_HANDLE%ML_INPUT_HANDLE%A(1:3,3)
      ENDDO
! Do the actual calculations here
      CALL VASP_TO_ML_MACHINE_LEARNING_LIB (ML_SUPER_HANDLE, ML_SUPER_HANDLE%ML_INPUT_HANDLE%A, ML_SUPER_HANDLE%ML_INPUT_HANDLE%B, &
           DIR_APP,DIR_LEN,INFO%LABORT, NSTEP, ML_SUPER_HANDLE%ML_INPUT_HANDLE%OMEGA, POSION_LOC, &
           TRIM(ADJUSTL(T_INFO%SZNAM2)), ML_SUPER_HANDLE%ML_INPUT_HANDLE%TOTEN_ML, ML_SUPER_HANDLE%ML_INPUT_HANDLE%TIFOR_ML, ML_SUPER_HANDLE%ML_INPUT_HANDLE%TSIF_ML )
    END SUBROUTINE VASP_TO_ML_MACHINE_LEARNING
#endif

!=======================================================================
!
! the following routine allows the ML code to terminate and
! write all data
!
!=======================================================================

    SUBROUTINE MACHINE_LEARNING_FINISH (ML_SUPER_HANDLE,DIR_APP,DIR_LEN,NSTEP)
      USE mpimy
      USE ML_FF_STRUCT, ONLY: ML_SUPER_TYPE
      IMPLICIT NONE
      TYPE (ML_SUPER_TYPE) :: ML_SUPER_HANDLE
      CHARACTER(LEN=10), INTENT(IN) :: DIR_APP
      INTEGER          , INTENT(IN) :: DIR_LEN
      INTEGER          , INTENT(IN) :: NSTEP ! Final MD step
      CALL MACHINE_LEARNING_FINISH_LIB (ML_SUPER_HANDLE,DIR_APP,DIR_LEN,NSTEP)
! Deallocate NITYP, POMASS and TYPE
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%NITYP)
      ENDIF
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%POMASS)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%POMASS)
      ENDIF
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TYPE)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TYPE)
      ENDIF
! Deallocate TIFOR_ML and TSIF_ML.
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TIFOR_ML)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TIFOR_ML)
      ENDIF
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TSIF_ML)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%TSIF_ML)
      ENDIF
! Deallocate VEL_ML.
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%VEL_ML)
      ENDIF
! Deallocate arrays storing the energies of isolated atoms.
      IF(ALLOCATED(ML_SUPER_HANDLE%ML_INPUT_HANDLE%EATOM_VASP)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%ML_INPUT_HANDLE%EATOM_VASP)
      ENDIF
! Deallocate ML_ICOUPLE that stores the list of atoms, where the coupling parameter is introduced for the calculation of the chemical potential.
      IF(ALLOCATED(ML_SUPER_HANDLE%FFM%ICOUPLE)) THEN
         DEALLOCATE(ML_SUPER_HANDLE%FFM%ICOUPLE)
      ENDIF
    END SUBROUTINE MACHINE_LEARNING_FINISH
#endif
#endif
END MODULE ml_interface

