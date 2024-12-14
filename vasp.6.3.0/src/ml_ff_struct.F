#include "symbol.inc"
!************************************************************************
!
!  This module contains all the types used in the 
!  machine-learning force field part.
!
!***********************************************************************
      MODULE ML_FF_STRUCT
      USE ML_FF_PREC
      USE ML_FF_TAGLIST
#ifdef use_shmem
      USE mpi_f08, ONLY : MPI_Win
      USE iso_c_binding
#endif

!-------------------------------------------------------------------------
!
! This structure contains all the input variables read in and also some
! related variables.
!
!-------------------------------------------------------------------------
      TYPE ML_FF_INPUT
! real space primitive lattice vectors
         REAL(q) :: A(1:3,1:3)    
! reciprocal space lattice vectors
         REAL(q) :: B(1:3,1:3) 
! communicator passed from VASP. Usually COMM_WORLD.
         INTEGER :: COMM_VASP              
! location of directory passed from VASP
         CHARACTER(LEN=10) :: DIR_APP 
! length of dir_app
         INTEGER           :: DIR_LEN     
! this arrays specify the atomic energies in isolated states.
         REAL(q), ALLOCATABLE :: EATOM_VASP(:)
         LOGICAL :: LABORT
! the number of MD steps
         INTEGER :: NSTEP
! the number of element types
         INTEGER :: NTYP
! the number of atoms for each species
         INTEGER, ALLOCATABLE :: NITYP(:)
! the number of atoms
         INTEGER :: NIONS
! this variable specifies the dimension of scaLAPACK grid
         INTEGER :: NDIM_SCALAPACK_FF
! this variable specifies the contents of ML_LOGFILE.
         INTEGER :: NWRITE_FF
! cell volume
         REAL(q) :: OMEGA
! the mass of atoms
         REAL(q), ALLOCATABLE :: POMASS(:)
! the length of each MD step
         REAL(q) :: POTIM
! name of current model system
         CHARACTER(LEN=40) :: SZNAM2
! total energy, force and stress tensor provided by ab initio calculations.
         REAL(q)              :: TOTEN_ML
         REAL(q), ALLOCATABLE :: TIFOR_ML(:,:)
         REAL(q), ALLOCATABLE :: TSIF_ML(:,:)
! type information for each species
         CHARACTER(LEN=2), ALLOCATABLE :: TYPE(:)
! velocity.
         REAL(q), ALLOCATABLE :: VEL_ML(:,:)
! Weights used for scaling the training data.
         REAL(q) :: WTIFOR_FF
         REAL(q) :: WTOTEN_FF
         REAL(q) :: WTSIF_FF
      END TYPE

!-------------------------------------------------------------------------------------
!
! This structure contains all the variables for ab-initio stuff (mainly
! from the ML_AB file.
!
!-------------------------------------------------------------------------------------

      TYPE ABINITIO
         REAL(q)          , ALLOCATABLE :: A(:,:,:) ! Lattice vectors. 1d:,2d:,3d:
         CHARACTER (LEN=9)              :: ABCAR ! Ab initio data file name.
         CHARACTER (LEN=10)             :: ABNCAR ! New ab initio data file name.
         REAL(q)          , ALLOCATABLE :: B(:,:,:) ! Reciprocal lattice vectors
         REAL(q)          , ALLOCATABLE :: CTIFOR_ALLCONF_OLD(:) ! old CTIFOR (Bayesian error criterion) for each training configuration
         REAL(q)          , ALLOCATABLE :: CTIFOR_ALLCONF_NEW(:) ! new CTIFOR (Bayesian error criterion) for each training configuration
         REAL(q)          , ALLOCATABLE :: EATOM(:) ! Atomic reference energies read from ML_AB file
         CHARACTER (LEN=10)             :: HISCAR ! Histogram file
         INTEGER          , ALLOCATABLE :: ISYS(:) ! Model system number.
         INTEGER                        :: IU1 ! Code number of ab initio data file.
         INTEGER                        :: IU2 ! Code number of new ab initio data file.
         INTEGER                        :: IU3 ! Code number of output file summarizing the regression results
         INTEGER                        :: IU4 ! Code number of histogram file
         INTEGER          , ALLOCATABLE :: LADD(:,:) ! Address of ion type in FF
         INTEGER          , ALLOCATABLE :: LB(:,:) ! Label of basis set
         INTEGER          , ALLOCATABLE :: LNITYP(:,:) ! List of atom number
         INTEGER          , ALLOCATABLE :: LNIONS(:,:,:) ! List of total atom number
         INTEGER          , ALLOCATABLE :: LNTYP(:,:) ! List of element number
         INTEGER          , ALLOCATABLE :: LSYS(:) ! Variable to identify whether the stored configuration is the same as the system handled by VASP currently
         INTEGER                        :: MCONF ! Maximum number of configuration
         INTEGER                        :: MCONF_ALLOCATE
         INTEGER                        :: MITYP ! Maximum number of atoms per specie
         INTEGER                        :: MITYP_ALLOCATE ! Size of MITYP used for allocation
         INTEGER          , ALLOCATABLE :: MITYP_CONF(:)
         INTEGER                        :: MIONS ! Maximum number of atoms per system
         INTEGER                        :: MIONS_ALLOCATE ! Size of MIONS used for allocation
         INTEGER          , ALLOCATABLE :: MIONS_CONF(:)
         INTEGER                        :: NCONF ! Number of configuration
         INTEGER          , ALLOCATABLE :: NCONF_SYS(:) ! The number of configuration per model system
         INTEGER          , ALLOCATABLE :: NITYP(:,:) ! Number of atoms for each specie
         INTEGER          , ALLOCATABLE :: NIONS(:) ! Total number of atom
         INTEGER                        :: NSYS
         INTEGER                        :: NSYS_ALLOCATE
         INTEGER          , ALLOCATABLE :: NTYP(:) ! Number of species
         INTEGER                        :: MTYP ! Maximum number of species
         INTEGER                        :: MTYP_ALLOCATE ! Size of MTYP used for allocation
         INTEGER          , ALLOCATABLE :: MTYP_CONF(:)
         REAL(q)          , ALLOCATABLE :: POMASS(:) ! Atomic mass
         REAL(q)          , ALLOCATABLE :: POSION(:,:,:) ! Atom positions
         CHARACTER (LEN=10)             :: REGCAR ! Output file summarizing the regression results.
         CHARACTER(LEN=40), ALLOCATABLE :: SZNAM2(:) ! Name of model system.
         CHARACTER(LEN=40), ALLOCATABLE :: SZNAM2_ADD(:) ! Name of model system.
         REAL(q)          , ALLOCATABLE :: TIFOR(:,:,:) ! Ab initio atom forces
         REAL(q)          , ALLOCATABLE :: TIFOR_MB(:,:,:) ! Machine-learning many-body interaction force
         REAL(q)          , ALLOCATABLE :: TOTEN(:) ! Ab initio total energies
         REAL(q)          , ALLOCATABLE :: TOTEN_MB(:)  ! Machine-learning many-body interaction energy
         REAL(q)          , ALLOCATABLE :: TSIF(:,:,:) ! Ab initio stress
         REAL(q)          , ALLOCATABLE :: TSIF_MB(:,:,:) ! Machine-learning many-body interaction stress
         CHARACTER (LEN=2), ALLOCATABLE :: TYPE(:) ! Element symbols
         REAL(q)          , ALLOCATABLE :: OMEGA(:) ! Volme of cell
      END TYPE

      TYPE ABINITIONEW
         REAL(q)          , ALLOCATABLE :: A(:,:,:)
         REAL(q)          , ALLOCATABLE :: B(:,:,:)
         REAL(q)          , ALLOCATABLE :: BE(:,:)
         REAL(q)          , ALLOCATABLE :: BEAV(:,:)  ! The average error in force estimated by Bayesian theorem.
         REAL(q)          , ALLOCATABLE :: BEMAX(:,:) ! The maximum error in force estimated by Bayesian theorem.
         REAL(q)          , ALLOCATABLE :: CTIFOR_ALLCONF(:)
         REAL(q)          , ALLOCATABLE :: EATOM_KIN(:) ! Atomic kinetic energies. If FFM%LEATOM=.TRUE. or FFM%LHEAT=.TRUE., they are calculated.
         REAL(q)          , ALLOCATABLE :: EATOM_POT(:) ! Atomic potential energies. If FFM%LEATOM=.TRUE. or FFM%LHEAT=.TRUE., they are calculated.
         INTEGER                        :: ISYS ! Model system number.
         INTEGER          , ALLOCATABLE :: LADD(:)
         INTEGER          , ALLOCATABLE :: LNITYP(:)
         INTEGER          , ALLOCATABLE :: LNIONS(:,:)
         INTEGER          , ALLOCATABLE :: LNTYP(:)
         INTEGER          , ALLOCATABLE :: LSYS(:)
         INTEGER                        :: MCONF ! If the number of new configuration reaches this number, execute force field generation
         INTEGER                        :: MITYP
         INTEGER                        :: MITYP_ALLOCATE
         INTEGER                        :: MIONS
         INTEGER                        :: MIONS_ALLOCATE
         INTEGER                        :: MTYP
         INTEGER                        :: MTYP_ALLOCATE
         INTEGER                        :: NCONF ! The number number of new configuration
         INTEGER          , ALLOCATABLE :: NITYP(:)
         INTEGER                        :: NIONS
         INTEGER          , ALLOCATABLE :: NSTEP(:)
         INTEGER                        :: NSYS
         INTEGER                        :: NSYS_ALLOCATE
         INTEGER                        :: NTYP
         REAL(q)          , ALLOCATABLE :: POMASS(:)
         REAL(q)          , ALLOCATABLE :: POSION(:,:,:)
         REAL(q)                        :: POTIM
         REAL(q)          , ALLOCATABLE :: QHEAT(:,:)
         CHARACTER(LEN=40)              :: SZNAM2
         CHARACTER(LEN=40), ALLOCATABLE :: SZNAM2_ADD(:)
         REAL(q)          , ALLOCATABLE :: TIFOR(:,:,:)
         REAL(q)          , ALLOCATABLE :: TIFOR_MB(:,:,:)
         REAL(q)          , ALLOCATABLE :: TOTEN(:)
         REAL(q)          , ALLOCATABLE :: TOTEN_MB(:)
         REAL(q)          , ALLOCATABLE :: TSIF(:,:,:)
         REAL(q)          , ALLOCATABLE :: TSIF_MB(:,:,:)
         CHARACTER (LEN=2), ALLOCATABLE :: TYPE(:)
         REAL(q)          , ALLOCATABLE :: OMEGA(:)
         REAL(q)          , ALLOCATABLE :: VEL(:,:,:)
      ENDTYPE

!--------------------------------------------------------------------------------------
! 
! This is the handle for force field variables. 
! 
!--------------------------------------------------------------------------------------

      TYPE FF_PAR
         REAL(q)          , ALLOCATABLE :: BEHIS(:) ! History of Bayesian error.
         REAL(q)                        :: CDOUB ! If the maximum Bayesian error is larger than CDOUB*CTIFOR
                                                 ! sampling is stopped, and force field is calculated.
         REAL(q)                        :: CSIG ! A criteia for variance to judge the convergence.
         REAL(q)                        :: CSIG_LOG ! Current criteria sigma, for logging only.
         REAL(q)                        :: CSLOPE ! A criteria for slope to judge the convergence.
         REAL(q)                        :: CSLOPE_LOG ! Current criteria slope, for logging only.
         REAL(q)                        :: CTIFOR ! A criteria for force. This parameter is not used, now.
         REAL(q)                        :: CTIFOR_PREV ! The threshold before the last update, used for logging.
         REAL(q)                        :: SCALE_CTIFOR ! Scaling factor for the threshold.
         REAL(q)                        :: CX ! Parameter x for criteria update.
         REAL(q)          , ALLOCATABLE :: EATOM_REF(:) ! Reference energy for atoms.
         REAL(q)          , ALLOCATABLE :: ERRHIS(:)
         CHARACTER (LEN=9)              :: FFCAR ! Input file summarizing old FF parameters. If ISTART=2, this file is read.
         CHARACTER (LEN=10)             :: FFNCAR ! Output file summarizing new FF parameters
         LOGICAL                        :: FORCING_UPDATE ! Forces learning, creation of new force field and eventually ab initio calculation
         INTEGER                        :: IALGO_LINREG ! Algorithm for linear regression
                                                        ! IALGO_LINREG=1, L2-norm with ridge regression
                                                        ! IALGO_LINREG=2, L1-norm with QR factorization
                                                        ! IALGO_LINREG=3, L1-norm with truncated SVD
                                                        ! IALGO_LINREG=4, SVD + Tikhonov regularization
         INTEGER                        :: ICRITERIA ! Parameter controlling update of threshold for Bayesian learning 
         LOGICAL                        :: IFF ! Flag for judging the presence of FF.
         LOGICAL                        :: IFF_OLD ! Old IFF data.
         INTEGER                        :: IMAT_SPARS ! Type of matrix used the sparsification.
                                                      ! IMAT_SPARS=1: Covariance matrix.
                                                      ! IMAT_SPARS=2: SOAP matrix.
         LOGICAL                        :: ISPARS ! ISPARS=.TRUE. means the successful execution of sparsification.
         INTEGER                        :: ISTART ! ISTART determines the mode of operation:
                                                  !   ISTART = 0 : Start the calculation without any previous data.
                                                  !   ISTART = 1 : Restart the calculation after reading the previous database file.
                                                  !   ISTART = 2 : Only force field is used without learning.
                                                  !   ISTART = 3 : New force-field created only from ML_AB file.
         INTEGER                        :: IU1 ! File code number for the old FF file
         INTEGER                        :: IU2 ! File code number for the new FF file
         INTEGER                        :: IU3 ! File code number for the local atomic energies
         INTEGER                        :: IU4 ! File code number for the local heat flux.
         INTEGER                        :: ISVD ! Type of the SVD for the sparsification.
                                                ! ISVD=1: Diagonalization.
                                                ! ISVD=2: Singular value decomposition.
         INTEGER                        :: IUPDATE_CRITERIA ! Decides whether Bayesian criterion is updated in same step or next
         INTEGER                        :: IWEIGHT ! This variable specifies the way to scale the training data.
                                                   ! If IWEIGHT = 1, scale the training data by using constants, WTIFOR, WTOTEN and WTSIF that are 
                                                   !  provided from INCAR file.
                                                   ! If IWEIGHT = 2, scale the training data by using variances in all training data.
                                                   ! If IWEIGHT = 3, scale the training data by using averages of variances in system-dependent training data.
                                                   ! If IWEIGHT = 4, scale the training data by using variances in all training data. 
                                                   ! In this scaling, the number of structure datasets are also taken into account.
                                                   ! If IWEIGHT = 5, scale the training data by using averages of variances in 
                                                   ! system-dependent training data. In this scaling, the number of structure datasets are also taken into account.
         LOGICAL                        :: LABINITIO ! Shows whether ab initio calculations are provided this turn or not
         INTEGER          , ALLOCATABLE :: LADD(:) ! Address of element
         LOGICAL                        :: LBASIS_DISCARD ! If LBASIS_DISCARD=.TRUE., basis sets are discarded when its number exceeds MB_FFM.
         LOGICAL                        :: LCONF_DISCARD ! If LCONF_DISCARD=.TRUE., configurations that do not provide any local reference configurations are discarded.
         LOGICAL                        :: LCTIFOR_READ_FROM_FILE ! If LCTIFOR_READ_FROM_FILE=.TRUE., then CTIFORs from ML_AB file are used
         LOGICAL                        :: LCTIFOR_USE_FROM_FILE ! Use CTIFOR from previous data
         LOGICAL                        :: LCTIFOR_WRITE_FROM_FILE ! Write CTIFOR for each structure
         LOGICAL                        :: LDISCARD_STRUCTURES_NOT_GIVING_BASIS
         LOGICAL                        :: LDOUB ! If LDOUB=.TRUE., the sampling is stopped, and a new force field is generated.
         LOGICAL                        :: LFMAT ! If LFMAT=.TRUE., design matrix elements are calculated
         LOGICAL                        :: LGENFF ! LGENFF=.TRUE. means that force field generation is executed.
         LOGICAL                        :: LLIST ! If LLIST=.TRUE., the machine generates lists storing the locations of design-matrix elements.
         LOGICAL                        :: LLIST_NEW ! If LLIST_NEW=.TRUE., the machine generates lists for storing the locations of 
                                                     ! desitn-matrix elements related to the candidates of new training data.
         LOGICAL                        :: LLOG_THRESHOLD_HISTORY ! If .TRUE. write threshold history to log file on next occasion.
         LOGICAL                        :: LLOG_THRESHOLD_UPDATE ! If .TRUE. write threshold update info to log file on next occasion.
         LOGICAL                        :: LMLONLY ! If LMLONLY=.TRUE., only the predictions by the machine-learning force field will be executed.
         LOGICAL                        :: LERR ! If LMLONLY=.TRUE. and LERR=.FALSE., the error estimation is skipped to accerelate the ML-MD
         LOGICAL                        :: LMLMB ! If LMLMB=.TRUE., the program generates force fields including many-body terms.
                                                 ! Actually, current program needs LMLMB=.TRUE.
         LOGICAL                        :: LNMDINT_RANDOM ! Decids wether NMDINT times a random number between 0 and 1 is used instead of NMDINT
         LOGICAL                        :: LSIG ! If LSIG=.TRUE., the machine calculates the variance of training data.
         LOGICAL                        :: LOVER ! LOVER=.TRUE. means that sparsification is executed to prevent the overcompleteness.
         LOGICAL                        :: LPS ! LPS=.TRUE. means that angular descriptor is calculated.
         LOGICAL                        :: LRUN_VASP_TO_ML ! Decides whether ML_TO_VASP is run or not.
         LOGICAL                        :: LSAMPLE ! If LSAMPLE=.TRUE., sampling is executed.
         LOGICAL                        :: LSECOND_TIME_ML_TO_VASP ! Helping variable to prevent incrementation of control variables
         LOGICAL                        :: LSOAP ! If LSOAP=.TRUE., SOAP matrix elements are calculated.
         LOGICAL                        :: LTEST ! If LTEST=.TRUE., ab initio calculation is executed at every MD step regardless the estimated error.
         LOGICAL                        :: LTOTEN_SYSTEM ! IF LTOTEN_SYSTEM=.TRUE. total energy of system is learned instead of energy/atom
         LOGICAL                        :: LTRJ ! If LTRJ=.TRUE., MD is executed using ab initio forces.
         INTEGER                        :: MCONF ! The number of maximum training configurations
         INTEGER                        :: MCONF_INIT ! MCONF value provided by INCAR file
         INTEGER                        :: MCONF_NEW ! The number of maximum new ab initio configurations
         INTEGER                        :: MHIS ! The maximum number of history of Bayesian errors stored in the memory.
         INTEGER                        :: MITYP ! Maximum number of atoms for each species
         INTEGER                        :: MITYP_ALLOCATE ! Size of MITYP used for allocation
         INTEGER                        :: MIONS ! Maximum number of atoms for each systeem
         INTEGER                        :: MIONS_ALLOCATE ! Size for allocation of MIONS
         CHARACTER (LEN=12)             :: ML_EATOM ! Output file containing the local atomic energies.
         REAL(q)                        :: ML_EPS_LOW ! Error tolerance EPS_LOW in SPARS_ICUR
         REAL(q)                        :: ML_EPS_REG ! Parameter for regularization in Bayesian linear regression in BLEA_MB
         CHARACTER (LEN=11)             :: ML_HEAT ! Output file containing the local heat flux.
         LOGICAL                        :: ML_FF_BINARY ! Keeps track whether ML_FF file is in binary format or not
         INTEGER                        :: ML_FF_VERSION_READ(1:3) ! Version number of ML_FF file read from ML_FF file, x.y.z format (x: major, y: minor, z:patch)
         INTEGER                        :: ML_FF_VERSION_WRITE(1:3) ! Version number of ML_FF file for writing, x.y.z format (x: major, y: minor, z:patch)
         INTEGER                        :: MTYP ! The number of atom type in force field file
         INTEGER                        :: MTYP_ALLOCATE ! Size for MTYP used for allocation
         INTEGER                        :: NHIS ! The number of history of Bayesian errors stored in the memory.
         INTEGER                        :: NBLOCK_FFN ! Block size for which ML_FFN file is printed out
         INTEGER                        :: NCOUNT_GENFF ! Counter how often force field was updated since last print of ML_FFN
         INTEGER                        :: NMDINT ! Interval of no-sampling.
         INTEGER                        :: NMDINT_RANDOM ! NMDINT multiplied with a random number
         INTEGER                        :: NSTEP ! step number of the MD
         INTEGER                        :: NSTEP_FROM_GENFF ! The number of MD step from the previous FF generation.
         INTEGER                        :: NSTEP_FROM_SAMPLING ! The number of MD step from the previous sampling.
         INTEGER                        :: NSYS ! The number of different systems in the traiing dataset.
         INTEGER                        :: NSYS_ALLOCATE
         INTEGER                        :: NSW ! Number of force field steps
         INTEGER                        :: NTEST ! This parameter determines how often the test ab initio calculation will be executed.
         INTEGER                        :: NTIFOR ! The number of ab initio force data.
         INTEGER                        :: NTOTEN ! The number of ab initio energy data.
         INTEGER                        :: NTSIF ! The number of ab initio stress tensor data.
         REAL(q)          , ALLOCATABLE :: POMASS(:) ! Mass of atom
         INTEGER                        :: SEED(1:3) ! Random seed used for random numbers
         REAL(q)                        :: SIG(1:10) ! Variance in the ab initio data.
         REAL(q)          , ALLOCATABLE :: SIG_SYS(:,:) ! Variance in the ab initio data for individual system.
         REAL(q)                        :: STIFOR ! Root mean square error in predicted force.
         REAL(q)                        :: STIFOR_LAST ! Root mean square error in predicted force for last added configuration
         REAL(q)                        :: STOTEN ! Root mean square error in predicted total energy.
         REAL(q)                        :: STOTEN_LAST ! Root mean square error in predicted total energy for last added configuration
         REAL(q)                        :: STSIF  ! Root mean square error in predicted stress.
         REAL(q)                        :: STSIF_LAST !Root mean square error in predicted stress for last added configuration
         CHARACTER(LEN=40), ALLOCATABLE :: SZNAM2(:) ! Name of the model system in the training dataset.
         REAL(q)                        :: TIFORAV(1:3) ! Average total forces.
         REAL(q)                        :: TOTENAV ! Average total energies.
         REAL(q)                        :: TSIFAV(1:3,1:3) ! Average total stress tensors.
         REAL(q)                        :: WTIFOR ! Weight for force.
         REAL(q)                        :: WTOTEN ! Weight for energy.
         REAL(q)                        :: WTSIF ! Weight for stress.
         CHARACTER (LEN=2), ALLOCATABLE :: TYPE(:) ! Type of ions
      END TYPE


!--------------------------------------------------------------------------------------
! 
! Handle for many-body variables of force field. 
! 
!--------------------------------------------------------------------------------------

      TYPE FFM_PAR ! Parameters used for many-body interaction fitting.
         REAL(q)              :: AFILT2 ! Parameter for angular filtering
         REAL(q), ALLOCATABLE :: C00(:,:,:) ! Expansion coefficients with l=m=0 used as pairwise descriptors for basis set.
         REAL(q), POINTER, CONTIGUOUS :: C00_ALL(:,:) ! Non-distributed C00 array.
         REAL(q), ALLOCATABLE :: C00_NEW(:,:,:,:) ! New expansion coefficients with l=m=0 used as pairwise descriptors.
         REAL(q), ALLOCATABLE :: CMAT(:,:) ! Covariance matrix that can be used for Bayesian error estimation.
         REAL(q), POINTER, CONTIGUOUS :: CMAT_ALL(:,:) ! Covariance matrix that can be used for Bayesian error estimation.
         REAL(q)              :: DHDRCOUPLE ! Derivative of Hamiltonian with respect to the coupling parameter.
         REAL(q)              :: DR1 ! Increment for radial mesh used for calculations of 1-point correlation (radial) descriptor.
         REAL(q)              :: DR2 ! Increment for radial mesh used for calculations of 2-points correlation (angular) descriptor.
         REAL(q)              :: DRSPL1 ! Increment for radial mesh used in spline-interpolation for 1-point correlation (radial) descriptor.
         REAL(q)              :: DRSPL2 ! Increment for radial mesh used in spline-interpolation for 2-points correlation (angular) descriptor.
         REAL(q), ALLOCATABLE :: DMAT_FMAT(:) ! Scaling vector for FMAT.
         REAL(q), ALLOCATABLE :: DMAT_YMAT(:) ! Scaling vector for FMAT.
         REAL(q), ALLOCATABLE :: FMAT(:,:,:) ! Design matrix.
         REAL(q), ALLOCATABLE :: FMAT_NEW(:,:,:) ! Design matrix for new configurations.
         REAL(q), ALLOCATABLE :: FNL1(:,:,:,:,:) ! Spline-interpolated radial functions used for calculating 1-point correlation (radial) descriptor.
         REAL(q), ALLOCATABLE :: FNL2(:,:,:,:,:) ! Spline-interpolated radial functions used for calculating 2-points correlation (angular) descriptor.
         REAL(q), ALLOCATABLE :: FWIN1(:,:) ! Window function used for calculating 1-point correlation (radial) descriptor.
         REAL(q), ALLOCATABLE :: FWIN2(:,:) ! Window function used for calculating 2-points correlation (angular) descriptor.
         INTEGER              :: IAFILT2 ! Type of angular filtering.
         INTEGER              :: IBROAD1 ! Type of broadening method for 1-point correlation (radial) descriptor.
         INTEGER              :: IBROAD2 ! Type of broadening method for 2-points correlation (angular) descriptor.
         INTEGER, ALLOCATABLE :: ICOUPLE(:) ! The list of atoms where the coupling parameter is introduced for the calculatoin of the chemical potential.
         INTEGER              :: ICUT1 ! If ICUT1=1, use Behler-Parinello cutoff function.
                                       ! If ICUT1=2, use Miwa's cutoff function.
                                       ! If ICUT1=3, used Jinnouchi's cutoff function.
         INTEGER              :: ICUT2 ! If ICUT2=1, use Behler-Parinello cutoff function.
                                       ! If ICUT2=2, use Miwa's cutoff function.
                                       ! If ICUT2=3, use Jinnouchi's cutoff function.
         INTEGER              :: INVERSE_SOAP ! If INVERSE_SOAP=1, diagonalization method is used.
                                              ! If INVERSE_SOAP=2, LU factorization is used.
         INTEGER              :: IREG ! Type of method to optimize the regularization parameters
         INTEGER              :: ISCALE_TOTEN ! Type of the way to scale the total energy for the many-body term.
                                                 !  ISCALE_TOTEN=1: Scale to the total energies of isolated atoms.
                                                 !  ISCALE_TOTEN=2: Scale to the average of the training energy data.
         LOGICAL              :: LSUPERVEC ! Type of the similarity measure.
                                             ! ML_FF_LSUPERVEC_MB = .FALSE. : The similarity measure written in 
                                             !                                "Jinnouchi et al., Phys. Rev. B, 
                                             !                                "On-the-fly machine learning force field generation: Application to melting points" (2019).
                                             ! ML_FF_LSUPERVEC_MB = .TRUE. : A new similarity measure using super vector.
         INTEGER              :: IWINDOW1 ! Type of window function used for calculating 1-point correlation (radial) descriptor.
                                          !   IWINDOW1=1 : Hann window.
                                          !   IWINDOW1=2 : Hamming window.
                                          !   IWINDOW1=3 : Blackman window.
                                          !   IWINDOW1=4 : Nuttall window
                                          !   IWINDOW1=5 : Blackman-Nutall window
                                          !   IWINDOW1=6 : Blackman-Harris window
         INTEGER              :: IWINDOW2 ! Type of window function used for calculating 2-points correlation (angular) descriptor.
                                          !   IWINDOW2=1 : Hann window.
                                          !   IWINDOW2=2 : Hamming window.
                                          !   IWINDOW2=3 : Blackman window.
                                          !   IWINDOW2=4 : Nuttall window
                                          !   IWINDOW2=5 : Blackman-Nutall window
                                          !   IWINDOW2=6 : Blackman-Harris window
         LOGICAL              :: LAFILT2   ! If LAFILT2=.TRUE., angular filtering is executed.
         INTEGER, ALLOCATABLE :: LB(:,:) ! Label of basis set.
         INTEGER, ALLOCATABLE :: LBHEAD(:) ! List of begining of basis sets
         LOGICAL              :: LCOUPLE ! If LCOUPLE=.TRUE., the coupling parameter is introduced for the calculation of the chemical potential.
         LOGICAL, ALLOCATABLE :: LIONS_COUPLE(:) ! For coupled atoms, LIONS_COUPLE is set to .TRUE.
         INTEGER, ALLOCATABLE :: L_LNRB1(:) ! List of angular momentum number used for calculating 1-point correlation (radial) descriptor.
         INTEGER, ALLOCATABLE :: L_LNRB2(:) ! List of angular momentum number used for calculating 2-points correlation (angular) descriptor.
         LOGICAL              :: LEATOM ! If LEATOM=.TRUE., atomic energies are calculated.
         LOGICAL, ALLOCATABLE :: LFLAG_VAR(:,:,:,:,:,:)
         LOGICAL, ALLOCATABLE :: LFLAG_VAR_SIC(:,:,:,:,:)
         LOGICAL              :: LHEAT ! It LHEAT=.TRUE., calculate heat flux
         INTEGER              :: LMAX1 ! Maximum angular momentum number for 1-point correlation (radial) descriptor. This is set as zero.
         INTEGER              :: LMAX2 ! Maximum angular momentum number for 2-points correlation (angular) descriptor.
         LOGICAL              :: LMETRIC1 ! A logical variable that switches on or off the metric function used for 
                                          ! calculating 1-point correlation (radial) descriptor.
         LOGICAL              :: LMETRIC2 ! A logical variable that switches on or off the metric function used for calculating 
                                          ! 2-points correlation (angular) descriptor.
         INTEGER, ALLOCATABLE :: LNCONF(:,:) ! Maps the local reference configuration index (NB(TYPE)) and a given atom type, to which
                                             ! training structure it belongs to.
                                             ! First index: local reference configuration index (NB), depends also on a type
                                             ! Second index: Atom type
                                             ! Result: training structure number (usually ICONF)
         INTEGER, ALLOCATABLE :: LNIONS(:,:) ! Maps the local reference configuration index (NB(TYPE)) and a given atom type, to which
                                             ! atom it belongs to. This is linked to LNCONF, so the atom should be one of the atoms
                                             ! of the corresponding LNCONF training structure.
                                             ! First index: local reference configuration index (NB), depends also on a type
                                             ! Second index: Atom type
                                             ! Result: atom number (usually INIONS) belonging to ICONF (from LNCONF(NB,TYPE))
         LOGICAL              :: LNORM1 ! If LNORM1=.TRUE., normalization is executed for 1-point correlation (radial) descriptor.
         LOGICAL              :: LNORM2 ! If LNORM2=.TRUE., normalization is executed for 2-points correlation (angular) descriptor.
         LOGICAL              :: LSIC ! If LSIC=.TRUE., self-interaction is removed from the angular descriptor.
         LOGICAL              :: LSPARSDES ! If LSPARSDES=.TRUE., descriptors are sparsified. This can be .TRUE., only when ISTART=1.
         INTEGER, ALLOCATABLE :: LTIFOR(:,:,:) ! List of training force data provided by the stored ab initio data.
                                               ! For a specific training set element it gives the map to it's
                                               ! position in the design matrix for example.
                                               ! First index: cartesian direction of force.
                                               ! Second index: respective atom for force.
                                               ! Third index: training configuration.
         INTEGER, ALLOCATABLE :: LTIFOR_NEW(:,:,:) ! List of training force data provided by the new configurations.
         INTEGER, ALLOCATABLE :: LTOTEN(:) ! List of training energy data provided by the stored ab initio data.
                                           ! Similar as LTIFOR, read description of that.
                                           ! First dimension: training configuration.
         INTEGER, ALLOCATABLE :: LTOTEN_NEW(:) ! List of training energy data provided by the new configurations.
         INTEGER, ALLOCATABLE :: LTSIF(:,:,:) ! List of training stress tensor data provided by the stored ab initio data.
                                              ! Similar as LTIFOR, read description of that.
                                              ! First dimension: First tensor direction.
                                              ! Second dimension: Second tensor direction.
                                              ! Third index: training configuration.
         INTEGER, ALLOCATABLE :: LTSIF_NEW(:,:,:) ! List of training stress tensor data provided by the new configurations.
         INTEGER, ALLOCATABLE :: LVAR(:,:,:,:,:,:) ! Label of angular descriptor.
         INTEGER, ALLOCATABLE :: LVAR_SIC(:,:,:,:,:) ! Label of angular descriptor.
         LOGICAL              :: LVARTRAN1 ! If LVARTRAN1=.TRUE., the radial coordinate is transformed.
         LOGICAL              :: LVARTRAN2 ! If LVARTRAN1=.TRUE., the radial coordinate is transformed.
         LOGICAL              :: LWINDOW1 ! If LWINDOW1=.TRUE., window function is used to smoothen the expanded atomic 
                                          ! distribution in calculations of 1-point correlation (radial) descriptor.
         LOGICAL              :: LWINDOW2 ! If LWINDOW2=.TRUE., window function is used to smoothen the expanded atomic 
                                          ! distribution in calculations of 2-points correlation (angular) descriptor.
         INTEGER              :: MB ! Maximum number of basis sets provided by the stored ab initio data.
         INTEGER              :: MB_INPUT
         INTEGER              :: MB_ALLOCATE
         INTEGER              :: MB_NEW ! Maximum number of basis sets provided by the new configurations.
         INTEGER              :: MBTOT ! Maximum total number of basis sets provided by the stored ab initio data.
         INTEGER              :: MBTOT_ALLOCATE
         INTEGER              :: MBTOT_NEW ! Maximum total number of basis sets provided by the new configurations.
         INTEGER              :: MDATA ! Maximum number of training data provided by the stored ab initio data.
         INTEGER              :: MDATA_ALLOCATE
         INTEGER              :: MDATA_NEW ! Maximum number of training data provided by the new configurations.
         INTEGER              :: MDATA_NEW_ALLOCATE
         INTEGER              :: MLNR1 ! (LMAX1+1)*MR1
         INTEGER              :: MLNR2 ! (LMAX2+1)*MR2
         INTEGER              :: MLNRB1 ! (LMAX1+1)*MRB1
         INTEGER              :: MLNRB2 ! (LMAX2+1)*MRB2
         INTEGER              :: MRB1 ! Maximum number of radial basis sets used in calculations of 1-point correlation (radial) descriptor.
         INTEGER              :: MRB2 ! Maximum number of radial basis sets used in calculations of 2-points correlation (angular) descriptor.
         INTEGER              :: MSPL1 ! Maximum number of radial grids used in spline-interpolation for calculations of 1-point correlation (radial) descriptor.
         INTEGER              :: MSPL2 ! Maximum number of radial grids used in spline-interpolation for calculations of 2-points correlation (angular) descriptor.
         INTEGER              :: MMVAR1 ! Maximum total number of variables for 1-point descriptor.
         INTEGER              :: MMVAR1_ALLOCATE
         INTEGER              :: MMVAR2 ! Maximum total number of variables for 2-point descriptor.
         INTEGER              :: MMVAR2_ALLOCATE
         INTEGER              :: MMVAR2_SIC ! Maximum total number of variables for self-interaction terms in 2-point (angular) descriptor.
         INTEGER              :: MMVAR2_SIC_ALLOCATE
         INTEGER              :: NATOM_COUPLED ! The number of atom where the coupling parameter is introduced.
         INTEGER, ALLOCATABLE :: NB(:) ! Number of basis set for each species provided by the stored ab initio data.
         INTEGER, ALLOCATABLE :: NB_NEW(:) ! The number of basis set for each species provided by the new configurations.
         INTEGER              :: NBTOT ! The total number of basis sets provided by the stored ab initio data.
         INTEGER              :: NBTOT_NEW ! The total number of basis sets provided by the new configurations.
         INTEGER              :: NDATA ! The number of training data provided by the stored ab initio data.
         INTEGER              :: NDATA_NEW ! The number of training data provided by the new configurations.
         INTEGER              :: NHYP1 ! Hyper-parameter for SOAP 1-point kernel.
         INTEGER              :: NHYP2 ! Hyper-parameter for SOAP 2-point kernel.
         INTEGER              :: NLNRB1 ! SUM_L[NRB1(L)]
         INTEGER              :: NLNRB2 ! SUM_L[NRB2(L)]
         INTEGER              :: NMETRIC1 ! Polynomical parameter for metric function used for calculating 1-point correlation (radial) descriptor.
         INTEGER              :: NMETRIC2 ! Polynomical parameter for metric function used for calculating 2-points correlation (angular) descriptor.
         INTEGER              :: NR1 ! The number of radial mesh used for calculating 1-point correlation (radial) descriptor.
         INTEGER              :: NR2 ! The number of radial mesh used for calculating 2-points correlation (angular) descriptor.
         INTEGER              :: NRANK_SPARSDES ! Rank of the matrix used in the sparsification of the descriptors.
         INTEGER, ALLOCATABLE :: NRB1(:) ! The number of radial basis set used for calculating 1-point correlation (radial) descriptor.
         INTEGER, ALLOCATABLE :: NRB2(:) ! The number of radial basis set used for calculating 2-points correlation (angular) descriptor.
         INTEGER, ALLOCATABLE :: NRB_LNRB1(:) ! List of basis set number used for calculating 1-point correlation (radial) descriptor.
         INTEGER, ALLOCATABLE :: NRB_LNRB2(:) ! List of basis set number used for calculating 2-points correlation (angular) descriptor.
         INTEGER              :: NSPL1 ! The number of radial grids used for spline interpolation used for calculating 1-point correlation (radial) descriptor.
         INTEGER              :: NSPL2 ! The number of radial grids used for spline interpolation used for calculating 2-points correlation (angular) descriptor.
         INTEGER              :: NTIFOR ! The total number of training force data.
         INTEGER              :: NTOTEN ! The total number of training energy data.
         INTEGER              :: NTSIF ! The total number of training stress-tensor data.
         INTEGER              :: NNVAR1 ! Total number of variables for 1-point descriptor.
         INTEGER              :: NNVAR1_ALLOCATE
         INTEGER, ALLOCATABLE :: NNVAR2(:) ! Total number of variables for 2-point descriptor.
         INTEGER, ALLOCATABLE :: NNVAR2_SIC(:) ! Total number of variables for self-interaction 2-point (angular) descriptor.
         INTEGER              :: NVARTRAN1 ! Polynomial parameter used for radial variable transform.
         INTEGER              :: NVARTRAN2 ! Polynomial parameter used for radial variable transform.
         REAL(q), ALLOCATABLE :: PS(:,:,:) ! Angular descriptor for basis set.
         REAL(q), POINTER, CONTIGUOUS :: PS_ALL(:,:) ! Non-distributed angular descriptor.
         REAL(q), ALLOCATABLE :: PS_NEW(:,:,:,:) ! Angular descriptor for basis set.
         REAL(q), ALLOCATABLE :: QGRID1(:,:) ! Q-grid for spherical Bessel functions used as radial basis sets for 1-point correlation (radial) descriptor.
         REAL(q), ALLOCATABLE :: QGRID2(:,:) ! Q-grid for spherical Bessel functions used as radial basis sets for 2-points correlation (angular) descriptor.
         REAL(q)              :: QMAX1 ! Maximum of Q-grid for 1-point correlation (radial) descriptor.
         REAL(q)              :: QMAX2 ! Maximum of Q-grid for 2-points correlation (angular) descriptor.
         REAL(q), ALLOCATABLE :: RB1(:,:,:) ! Radial basis set for 1-point correlation (radial) descriptor.
         REAL(q), ALLOCATABLE :: RB2(:,:,:) ! Radial basis set for 2-points correlation (angular) descriptor.
         REAL(q), ALLOCATABLE :: RC1(:) ! Positions of gaussians for 1-point correlation (radial) descriptor.
         REAL(q), ALLOCATABLE :: RC2(:) ! Positions of gaussians for 2-points correlation (angular) descriptor.
         REAL(q)              :: RCOUPLE ! Coupliing parameter for the calculation of the chemical potential.
         REAL(q), ALLOCATABLE :: RCOUPLE_IONS(:) ! Coupling parameters for each ions.
         REAL(q), ALLOCATABLE :: RG1(:) ! Radial grid for 1-point correlation (radial) descriptor.
                                        ! When metric is introduced, this grid corresponds to the equispaced grid over R^N.
                                        ! Otherwise, this grid corresponds to the equispaced over R.
         REAL(q), ALLOCATABLE :: RG2(:) ! Radial grid for 2-points correlation (angular) descriptor.
                                        ! When metric is introduced, this grid corresponds to the equispaced grid over R^N.
                                        ! Otherwise, this grid corresponds to the equispaced over R.
         REAL(q)              :: RCUT1 ! Cutoff radius for 1-point correlation (radial) descriptor.
         REAL(q)              :: RCUT2 ! Cutoff radius 2-points correlation (angular) descriptor.
         REAL(q)              :: RDES_SPARSDES ! Ratio of the number of selecting descriptors to all.
         REAL(q)              :: RMETRIC1 ! A parameter for metric function for 1-point correlation (radial) descriptor.
         REAL(q)              :: RMETRIC2 ! A parameter for metric function for 2-points correlation (angular) descriptor.
         REAL(q), ALLOCATABLE :: RNORM1(:,:) ! Norm of Bessel function used for calculating 1-point correlation (radial) descriptor.
         REAL(q), ALLOCATABLE :: RNORM2(:,:) ! Norm of Bessel function used for calculating 2-points correlation (angular) descriptor.
         REAL(q)              :: SIG(1:10) ! Variations in first principles data.
         REAL(q), ALLOCATABLE :: SIG_SYS(:,:) ! Variations in first principles data for individual system.
         REAL(q)              :: SIGV ,SIGW  ! Pricision parameters determined by the evidence approximation.
         REAL(q)              :: SIGV0,SIGW0 ! (Initial) pricision parameters determined by the evidence approximation.
         REAL(q)              :: SION1 ! A parameter that specifies the broadening width for 1-point correlation (radial) descriptor.
         REAL(q)              :: SION2 ! A parameter that specifies the broadening width for 2-points correlation (angular) descriptor.
         REAL(q), ALLOCATABLE :: SOAP(:,:,:) ! SOAP matrix used for sparsification.
         REAL(q)              :: TIFORAV(1:3) ! Average total force given from many-body interactions.
         REAL(q)              :: TOTENAV ! Average total energy given from many-body interactions.
         REAL(q)              :: TSIFAV(1:3,1:3) ! Average total stress given from many-body interactions.
         REAL(q), ALLOCATABLE :: USQ1(:,:) ! Orthogonalizing coefficients for radial Gaussian basis sets.
         REAL(q), ALLOCATABLE :: USQ2(:,:) ! Orthogonalizing coefficients for radial Gaussian basis sets.
         REAL(q)              :: W1 ! Weight for 1-point descriptor in SOAP descriptors.
         REAL(q)              :: W2 ! Weight for 2-point descriptor in SOAP descriptors.
         REAL(q), ALLOCATABLE :: WION1(:) ! Broadening width used for calculating 1-point correlation (radial) descriptor.
         REAL(q), ALLOCATABLE :: WION2(:) ! Broadening width used for calculating 2-points correlation (angular) descriptor.
         REAL(q), ALLOCATABLE :: WMAT(:,:,:) ! Regression coefficients.
         REAL(q), ALLOCATABLE :: WMAT_ALL(:) ! Regression coefficients owned by all processes.
         REAL(q)              :: WR1 ! WR1 specifies the radial resolution of the radial spherical Bessel function basis sets.
         REAL(q)              :: WR2 ! WR2 specifies the radial resolution of the radial spherical Bessel function basis sets.
         REAL(q), ALLOCATABLE :: WVAR(:,:,:,:,:,:) ! Weight for angular descriptor.
         REAL(q), ALLOCATABLE :: YMAT(:,:) ! Data vector
      END TYPE


!--------------------------------------------------------------------------------------
! 
! Important variables used for the logfile.
! 
!--------------------------------------------------------------------------------------

      TYPE LOG_PAR
         INTEGER           :: IU
         CHARACTER(LEN=14) :: ML_LOGFILE
         INTEGER           :: NWRITE
      END TYPE


!--------------------------------------------------------------------------------------
! 
! Variables needed for gathering information about memory usage.
! 
!--------------------------------------------------------------------------------------

      TYPE MEMORY
         REAL(q) RMEM(1:11)
         REAL(q) RMEM_CMAT ! Memory used for storing covariance matrix of many-body interactions.
         REAL(q) RMEM_CMAT_ALL ! Memory used for storing covariance matrix of many-body interactions.
         REAL(q) RMEM_FMAT ! Memory used for storing design matrix of many-body interactions.
         REAL(q) RMEM_PS ! Memory used for storing angular descriptor.
         REAL(q) RMEM_PS_ALL ! Memory used for storing non-distributed angular descriptor.
         REAL(q) RMEM_RST ! Memory necessary for restructuring of the arrays in sparsification process
         REAL(q) RMEM_SOAP ! Memory used for storing inverse SOAP matrix necessary for error estimation
         REAL(q) RMEM_SOAP_ALL ! Memory used for storing inverse SOAP matrix necessary for error estimation
         REAL(q) RMEM_SVD ! Memory used for SVD decomposition
         REAL(q) RMEM_WMAT_ALL ! Memory used for storing non-distributed regression coefficients.
      END TYPE MEMORY


!--------------------------------------------------------------------------------------
! 
! A lot of parallel stuff. 
! 
!--------------------------------------------------------------------------------------
      TYPE ML_MPI_PAR
           INTEGER :: MPI_COMM     ! MPI_Communicator
           INTEGER :: NCPU         ! total number of proc in this communicator
           INTEGER :: NODE_ME      ! node id starting from 1 ... NCPU
           INTEGER :: COLOR        ! color of process within (sub)communicator
      END TYPE ML_MPI_PAR

      TYPE SCALAPACK_GRID
         INTEGER :: ICTXT        ! context number
         INTEGER :: NDIM         ! dimension of scaLAPACK grid
         INTEGER :: NP_ROW       ! number of rows in processor grid
         INTEGER :: NP_COL       ! number of columns in processor grid
         INTEGER :: ME_ROW       ! row coordinate in processor grid
         INTEGER :: ME_COL       ! column coordinate in processor grid
         INTEGER :: NPROCS       ! number of processors in context
      END TYPE

      TYPE PARALLEL_DESCRIPTOR
         INTEGER :: DESC(1:9)    ! descriptor number
         INTEGER :: M            ! number of rows in global array
         INTEGER :: N            ! number of columns in global array
         INTEGER :: MB           ! blocking factor for rows
         INTEGER :: NB           ! blocking factor for columns
         INTEGER :: IRSRC        ! process row of first row of matrix distribution
         INTEGER :: ICSRC        ! process column of first column of matrix distribution
         INTEGER :: LLD          ! Leading dimension of local array(MAX(1, LOCr(M))
         INTEGER :: INFO         ! info
         INTEGER :: LOC_M        ! local size of M specific to one process
         INTEGER :: LOC_N        ! local size of N specific to one process
      END TYPE

      TYPE SCALA_WORK
         INTEGER              :: INFO
         INTEGER, ALLOCATABLE :: IPIV(:)
         INTEGER, ALLOCATABLE :: IWORK(:)
         INTEGER              :: LOC_M
         INTEGER              :: LWORK
         INTEGER              :: LIWORK
         REAL(q), ALLOCATABLE :: WORK(:)
      END TYPE

#ifdef use_shmem
      TYPE M_shmem_desc
#ifndef sysv
         TYPE (MPI_Win)           :: shmem_win
#else
         INTEGER(KIND=c_int)      :: shmem_win
#endif
         TYPE (c_ptr)             :: shmem_ptr
         TYPE (ML_MPI_PAR), POINTER  :: comm
         INTEGER                  :: dim1=0, dim2=0, dim3=0
      END TYPE M_shmem_desc

      TYPE M_shmem
         TYPE (M_shmem_desc) :: desc
         COMPLEX(q), POINTER :: fptr_c(:) => NULL()
         REAL(q), POINTER    :: fptr_r(:) => NULL()
      END TYPE M_shmem

#endif


      TYPE PARALLEL_SUPER
         TYPE(ML_MPI_PAR)                :: COMM_WORLD
#ifdef use_shmem
         TYPE(ML_MPI_PAR)                :: COMM_inter_node_world
         TYPE(ML_MPI_PAR)                :: COMM_intra_node_world
#endif
         TYPE(ML_MPI_PAR)                :: COMM_NB
         TYPE(ML_MPI_PAR)                :: COMM_NIONS
         TYPE(ML_MPI_PAR)                :: COMM_NEIB
         TYPE(SCALAPACK_GRID)            :: CONTEXT_WORLD
         TYPE(SCALAPACK_GRID)            :: CONTEXT_SINGLE
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_A_SPARSDES_CUR
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_BE_ENERGY
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_BE_FORCE
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_BE_STRESS
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_C00
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_C00_FFM
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_C00_NEW_FFM
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_CMAT
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_CMAT0
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_CMAT_FFM
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_DC00
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_DCLM1
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_DCLM2
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_DPS
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_DSOAP1
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_DSOAP2
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_EATOM
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_ENERGY_HELP
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_FORCE_HELP
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_FMAT_FFM
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_FMAT_NEW_FFM
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_FMAT_FFM_TRANSPOSED
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_GMAT
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_GENERAL_SINGLE
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_HEAT1
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_HEAT2
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_PS
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_PS_FFM
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_PS_NEW_FFM
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_RBS
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_SOAP
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_SOAP_FFM
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_SOAP_NEW
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_SOAP_NEW_TRANSPOSED
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_SOAP_HELP
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_SOAP_TRANSPOSED
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_STRESS_HELP
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_V_SPARSDES_CUR
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_VEC_LOC_ENERGY
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_VEC_LOC_FORCE1
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_VEC_LOC_FORCE2
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_VEC_LOC_HEAT1
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_VEC_LOC_HEAT2
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_VEC_LOC_STRESS1
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_VEC_LOC_STRESS2
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_WMAT
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_WMAT_FFM
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_YMAT_FFM
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_ZMAT
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_ZMAT_HEAT
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_FMAT_HELP
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_FMAT_TRANS_HELP
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_SING_VT
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_SING_U
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_SING_SIGMA
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_SING_VSIG
         TYPE(PARALLEL_DESCRIPTOR)       :: DESC_SING_B1_SVD
         TYPE(SCALA_WORK)                :: SCAWORK
         REAL(q),ALLOCATABLE             :: GENERAL_SINGLE(:,:) ! Array mainly to reduce distributed 
      END TYPE PARALLEL_SUPER

!--------------------------------------------------------------------------------------
! 
! All important shared memory variables used in ml_ff_ff.F. 
! 
!--------------------------------------------------------------------------------------
      
      TYPE ML_SHMEM
         REAL(q),POINTER,CONTIGUOUS :: CLM1(:,:)
         REAL(q),POINTER,CONTIGUOUS :: CLM2(:,:)
         REAL(q),POINTER,CONTIGUOUS :: CLM_COUPLE_HELP1_MANYBODY(:,:)
         REAL(q),POINTER,CONTIGUOUS :: CLM_COUPLE_HELP2_MANYBODY(:,:)
         REAL(q),POINTER,CONTIGUOUS :: PS_SIC(:,:)
         REAL(q),POINTER,CONTIGUOUS :: PS_SIC_COUPLE(:,:)
         REAL(q),POINTER,CONTIGUOUS :: SOAP_ALL(:,:)
         REAL(q),POINTER,CONTIGUOUS :: SOAP_C00(:,:)
         REAL(q),POINTER,CONTIGUOUS :: SOAP_PS(:,:)
         INTEGER                    :: MCOL_CLM1_GLOBAL
         INTEGER                    :: MCOL_CLM2_GLOBAL
         INTEGER                    :: MCOL_CLM_COUPLE_HELP1_MANYBODY
         INTEGER                    :: MCOL_CLM_COUPLE_HELP2_MANYBODY
         INTEGER                    :: MCOL_PS_SIC_GLOBAL
         INTEGER                    :: MCOL_PS_SIC_COUPLE_GLOBAL
         INTEGER                    :: MCOL_SOAP_ALL_GLOBAL
         INTEGER                    :: MCOL_SOAP_C00_GLOBAL
         INTEGER                    :: MCOL_SOAP_PS_GLOBAL
         INTEGER                    :: MROW_CLM1_GLOBAL
         INTEGER                    :: MROW_CLM2_GLOBAL
         INTEGER                    :: MROW_CLM_COUPLE_HELP1_MANYBODY
         INTEGER                    :: MROW_CLM_COUPLE_HELP2_MANYBODY
         INTEGER                    :: MROW_PS_SIC_GLOBAL
         INTEGER                    :: MROW_PS_SIC_COUPLE_GLOBAL
         INTEGER                    :: MROW_SOAP_ALL_GLOBAL
         INTEGER                    :: MROW_SOAP_C00_GLOBAL
         INTEGER                    :: MROW_SOAP_PS_GLOBAL
#ifdef use_shmem
         TYPE (M_shmem) :: MLFF_SHM_2D_C00_ALL
         TYPE (M_shmem) :: MLFF_SHM_2D_PS_ALL
         TYPE (M_shmem) :: MLFF_SHM_2D_CLM1
         TYPE (M_shmem) :: MLFF_SHM_2D_CLM2
         TYPE (M_shmem) :: MLFF_SHM_2D_CLM_COUPLE_HELP1_MANYBODY
         TYPE (M_shmem) :: MLFF_SHM_2D_CLM_COUPLE_HELP2_MANYBODY
         TYPE (M_shmem) :: MLFF_SHM_2D_PS_SIC
         TYPE (M_shmem) :: MLFF_SHM_2D_PS_SIC_COUPLE
         TYPE (M_shmem) :: MLFF_SHM_2D_SOAP_ALL
         TYPE (M_shmem) :: MLFF_SHM_2D_SOAP_C00
         TYPE (M_shmem) :: MLFF_SHM_2D_SOAP_PS
         TYPE (M_shmem) :: MLFF_SHM_2D_CMAT_ALL
         LOGICAL        :: MLFF_SHMEM_ALLOC_C00_ALL=.FALSE.
         LOGICAL        :: MLFF_SHMEM_ALLOC_PS_ALL=.FALSE.
         LOGICAL        :: MLFF_SHMEM_ALLOC_CLM1=.FALSE.
         LOGICAL        :: MLFF_SHMEM_ALLOC_CLM2=.FALSE.
         LOGICAL        :: MLFF_SHMEM_ALLOC_CLM_COUPLE_HELP1_MANYBODY=.FALSE.
         LOGICAL        :: MLFF_SHMEM_ALLOC_CLM_COUPLE_HELP2_MANYBODY=.FALSE.
         LOGICAL        :: MLFF_SHMEM_ALLOC_PS_SIC=.FALSE.
         LOGICAL        :: MLFF_SHMEM_ALLOC_PS_SIC_COUPLE=.FALSE.
         LOGICAL        :: MLFF_SHMEM_ALLOC_SOAP_ALL=.FALSE.
         LOGICAL        :: MLFF_SHMEM_ALLOC_SOAP_C00=.FALSE.
         LOGICAL        :: MLFF_SHMEM_ALLOC_SOAP_PS=.FALSE.
         LOGICAL        :: MLFF_SHMEM_ALLOC_CMAT_ALL=.FALSE.
#else
         REAL(q),POINTER :: C00_ALL_NOSHMEM(:,:) => NULL()
         REAL(q),POINTER :: PS_ALL_NOSHMEM(:,:) => NULL()
         REAL(q),POINTER :: CLM1_NOSHMEM(:,:) => NULL()
         REAL(q),POINTER :: CLM2_NOSHMEM(:,:) => NULL()
         REAL(q),POINTER :: CLM_COUPLE_HELP1_MANYBODY_NOSHMEM(:,:) => NULL()
         REAL(q),POINTER :: CLM_COUPLE_HELP2_MANYBODY_NOSHMEM(:,:) => NULL()
         REAL(q),POINTER :: PS_SIC_NOSHMEM(:,:) => NULL()
         REAL(q),POINTER :: PS_SIC_COUPLE_NOSHMEM(:,:) => NULL()
         REAL(q),POINTER :: SOAP_ALL_NOSHMEM(:,:) => NULL()
         REAL(q),POINTER :: SOAP_C00_NOSHMEM(:,:) => NULL()
         REAL(q),POINTER :: SOAP_PS_NOSHMEM(:,:) => NULL()
         REAL(q),POINTER :: CMAT_ALL_NOSHMEM(:,:) => NULL()
#endif
         INTEGER                    :: ISTATUS_SHMEM
         INTEGER                    :: MROW_CLM1_GLOBAL_ALLOCATE
         INTEGER                    :: MCOL_CLM1_GLOBAL_ALLOCATE
         INTEGER                    :: MROW_CLM_COUPLE_HELP1_MANYBODY_ALLOCATE
         INTEGER                    :: MCOL_CLM_COUPLE_HELP1_MANYBODY_ALLOCATE
         INTEGER                    :: MROW_CLM2_GLOBAL_ALLOCATE
         INTEGER                    :: MCOL_CLM2_GLOBAL_ALLOCATE
         INTEGER                    :: MROW_CLM_COUPLE_HELP2_MANYBODY_ALLOCATE
         INTEGER                    :: MCOL_CLM_COUPLE_HELP2_MANYBODY_ALLOCATE
         INTEGER                    :: MROW_PS_SIC_GLOBAL_ALLOCATE
         INTEGER                    :: MCOL_PS_SIC_GLOBAL_ALLOCATE
         INTEGER                    :: MROW_PS_SIC_COUPLE_GLOBAL_ALLOCATE
         INTEGER                    :: MCOL_PS_SIC_COUPLE_GLOBAL_ALLOCATE
         INTEGER                    :: MROW_SOAP_ALL_GLOBAL_ALLOCATE
         INTEGER                    :: MCOL_SOAP_ALL_GLOBAL_ALLOCATE
         INTEGER                    :: MROW_SOAP_C00_GLOBAL_ALLOCATE
         INTEGER                    :: MCOL_SOAP_C00_GLOBAL_ALLOCATE
         INTEGER                    :: MROW_SOAP_PS_GLOBAL_ALLOCATE
         INTEGER                    :: MCOL_SOAP_PS_GLOBAL_ALLOCATE
         
      END TYPE ML_SHMEM
        
!--------------------------------------------------------------------------------------
! 
! Helping arrays neccessary for LAPACK operations.
! 
!--------------------------------------------------------------------------------------
      TYPE LAPACK_WORK
         INTEGER              :: INFO
         INTEGER, ALLOCATABLE :: IPIV(:)
         INTEGER              :: LWORK
         REAL(q), ALLOCATABLE :: WORK(:)
      END TYPE 

!--------------------------------------------------------------------------------------
! 
! Variables storing timing information.
! 
!--------------------------------------------------------------------------------------
      TYPE TIME
         REAL(q) CPUTIM(1:12)
         REAL(q) VPUTIM(1:12)
      END TYPE TIME

!--------------------------------------------------------------------------------------
! 
! This is the handle super type used from main. 
! 
!--------------------------------------------------------------------------------------
      TYPE ML_SUPER_TYPE
         TYPE (ABINITIO)       :: AB
         TYPE (ABINITIONEW)    :: ABN
         TYPE (ML_FF_INPUT)    :: ML_INPUT_HANDLE
         TYPE (FF_PAR)         :: FF
         TYPE (FFM_PAR)        :: FFM
         TYPE (TAGLIST)        :: TAG_LIST
         TYPE (LOG_PAR)        :: LOGF
         TYPE (MEMORY)         :: MEM0
         TYPE (MEMORY)         :: MEM1
         TYPE (ML_SHMEM)       :: MLSHM
         TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
         TYPE (LAPACK_WORK)    :: LAPACKWORK
         TYPE (TIME)           :: TIM
         INTEGER               :: INSTANCE
         LOGICAL               :: SKIP_ML_TO_VASP = .FALSE.
         LOGICAL               :: SKIP_VASP_TO_ML = .FALSE.
      END TYPE

! Global variables
      TYPE (ML_SUPER_TYPE),ALLOCATABLE :: ML_SUPER_HANDLE_MAIN(:)
      TYPE (ML_MPI_PAR)         :: ML_COMM_WORLD_GLOBAL
      LOGICAL                   :: ML_IO_WRITE = .FALSE.
      INTEGER                   :: ML_FILESTREAM_GLOBAL
      INTEGER                   :: ML_TOTNUM_INSTANCES

      END MODULE
