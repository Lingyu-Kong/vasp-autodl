!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Contains list of tags, descriptions and taglist derived type.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ml_ff_taglist

   implicit none
   private
   public :: taglist,                       &
             taglist_get_tag_width,         &
             taglist_description,           &
             taglist_get_all_tags,          &
             taglist_get_alt_tags,          &
             taglist_get_description_width, &
             taglist_state_symbol

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !> List of MLFF input tags and related convenience functions
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   enum, bind(c)
      enumerator :: TL_STATE     = 0
      enumerator :: TL_UNSET     = 1
      enumerator :: TL_DEFAULT   = 2
      enumerator :: TL_INCAR     = 3
      enumerator :: TL_INCAR_ALT = 4
      enumerator :: TL_OVERRIDE  = 5
   end enum
   ! Note: There is no way to make enums public without listing them.
   public :: TL_STATE,     &
             TL_UNSET,     &
             TL_DEFAULT,   &
             TL_INCAR,     &
             TL_INCAR_ALT, &
             TL_OVERRIDE

   type taglist

      private
      integer(kind(TL_STATE)), allocatable :: state(:)
      integer                , allocatable :: alt_index(:)

      contains

      procedure, public, pass(this) :: init => taglist_init
      procedure, public, pass(this) :: get_state => taglist_get_state
      procedure, public, pass(this) :: set_state => taglist_set_state
      procedure, public, pass(this) :: get_deprecated => taglist_get_deprecated
      procedure, public, pass(this) :: symbol => taglist_state_symbol_from_tag

   end type taglist

   !> Number of tags listed here.
   integer, parameter                  :: ntags = 94
   !> Maximum number of alternative tags listed here.
   integer, parameter                  :: nalt = 2
   !> Static list of all known MLFF tags.
   ! WARNING: If this list needs to be extended with more alternatives, make
   ! sure to fill it left to right, i.e. empty strings should only appear at
   ! the end of each line (e.g. "A", "B", " " and NOT "A", " ", "B"). The first
   ! entry is considered to be the "preferred" tag, all others are deprecated.
   character(len=*), parameter         :: tags(nalt, ntags) = reshape([                            &
      ! General, mode of operation:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "ML_LMLFF                                  ", "ML_FF_LMLFF                               ",  &
      "ML_ISTART                                 ", "ML_FF_ISTART                              ",  &

      ! Descriptors:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "ML_RCUT1                                  ", "ML_FF_RCUT1_MB                            ",  &
      "ML_SION1                                  ", "ML_FF_SION1_MB                            ",  &
      "ML_MRB1                                   ", "ML_FF_MRB1_MB                             ",  &
      "ML_RCUT2                                  ", "ML_FF_RCUT2_MB                            ",  &
      "ML_SION2                                  ", "ML_FF_SION2_MB                            ",  &
      "ML_MRB2                                   ", "ML_FF_MRB2_MB                             ",  &
      "ML_LMAX2                                  ", "ML_FF_LMAX2_MB                            ",  &
      "ML_LAFILT2                                ", "ML_FF_LAFILT2_MB                          ",  &
      "ML_AFILT2                                 ", "ML_FF_AFILT2_MB                           ",  &
      "ML_IAFILT2                                ", "ML_FF_IAFILT2_MB                          ",  &
      "ML_LSPARSDES                              ", "ML_FF_LSPARSDES_MB                        ",  &
      "ML_NRANK_SPARSDES                         ", "ML_FF_NRANK_SPARSDES_MB                   ",  &
      "ML_RDES_SPARSDES                          ", "ML_FF_RDES_SPARSDES_MB                    ",  &

      ! Kernel:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "ML_W1                                     ", "ML_FF_W1_MB                               ",  &
      "ML_NHYP                                   ", "ML_FF_NHYP1_MB                            ",  &

      ! Controlling Bayesian criterion:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "ML_ICRITERIA                              ", "ML_FF_ICRITERIA                           ",  &
      "ML_IUPDATE_CRITERIA                       ", "ML_FF_IUPDATE_CRITERIA                    ",  &
      "ML_LCTIFOR_USE_FROM_FILE                  ", "ML_FF_LCTIFOR_USE_FROM_FILE               ",  &
      "ML_CTIFOR                                 ", "ML_FF_CTIFOR                              ",  &
      "ML_SCLC_CTIFOR                            ", "ML_FF_SCLC_CTIFOR                         ",  &
      "ML_CSIG                                   ", "ML_FF_CSIG                                ",  &
      "ML_CSLOPE                                 ", "ML_FF_CSLOPE                              ",  &
      "ML_CX                                     ", "ML_FF_XMIX                                ",  &
      "ML_CDOUB                                  ", "ML_FF_CDOUB                               ",  &

      ! Sparsification and regression:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "ML_EPS_LOW                                ", "ML_FF_EPS_LOW                             ",  &
      "ML_EPS_REG                                ", "ML_FF_EPS_REG                             ",  &
      "ML_IALGO_LINREG                           ", "ML_FF_IALGO_LINREG                        ",  &
      "ML_ISVD                                   ", "ML_FF_ISVD                                ",  &
      "ML_IREG                                   ", "ML_FF_IREG_MB                             ",  &
      "ML_SIGV0                                  ", "ML_FF_SIGV0_MB                            ",  &
      "ML_SIGW0                                  ", "ML_FF_SIGW0_MB                            ",  &

      ! Weighting and sampling:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "ML_WTOTEN                                 ", "ML_FF_WTOTEN                              ",  &
      "ML_WTIFOR                                 ", "ML_FF_WTIFOR                              ",  &
      "ML_WTSIF                                  ", "ML_FF_WTSIF                               ",  &
      "ML_MHIS                                   ", "ML_FF_MHIS                                ",  &
      "ML_NMDINT                                 ", "ML_FF_NMDINT                              ",  &
      "ML_IWEIGHT                                ", "ML_FF_IWEIGHT                             ",  &

      ! Sizes of static arrays:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "ML_MB                                     ", "ML_FF_MB_MB                               ",  &
      "ML_MCONF                                  ", "ML_FF_MCONF                               ",  &
      "ML_MCONF_NEW                              ", "ML_FF_MCONF_NEW                           ",  &

      ! Special features:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "ML_LHEAT                                  ", "ML_FF_LHEAT_MB                            ",  &
      "ML_LCOUPLE                                ", "ML_FF_LCOUPLE_MB                          ",  &
      "ML_NATOM_COUPLED                          ", "ML_FF_NATOM_COUPLED_MB                    ",  &
      "ML_ICOUPLE                                ", "ML_FF_ICOUPLE_MB                          ",  &
      "ML_RCOUPLE                                ", "ML_FF_RCOUPLE_MB                          ",  &

      ! Reference energies:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "ML_ISCALE_TOTEN                           ", "ML_FF_ISCALE_TOTEN_MB                     ",  &
      "ML_EATOM_REF                              ", "ML_FF_EATOM                               ",  &

      ! Output options:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "ML_LEATOM                                 ", "ML_FF_LEATOM_MB                           ",  &

      ! Not important, description not used and not yet reviewed:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "ML_IBROAD1                                ", "ML_FF_IBROAD1_MB                          ",  &
      "ML_IBROAD2                                ", "ML_FF_IBROAD2_MB                          ",  &
      "ML_TOTNUM_INSTANCES                       ", "ML_FF_TOTNUM_INSTANCES                    ",  &
      "ML_FILESTREAM_START                       ", "ML_FF_FILESTREAM_START                    ",  &
      "ML_IMAT_SPARS                             ", "ML_FF_IMAT_SPARS                          ",  &
      "ML_LBASIS_DISCARD                         ", "ML_FF_LBASIS_DISCARD                      ",  &
      "ML_LCONF_DISCARD                          ", "ML_FF_LCONF_DISCARD                       ",  &
      "ML_LMLMB                                  ", "ML_FF_LMLMB                               ",  &
      "ML_LERR                                   ", "ML_FF_LERR                                ",  &
      "ML_NDIM_SCALAPACK                         ", "ML_FF_NDIM_SCALAPACK                      ",  &
      "ML_FF_NWRITE                              ", "                                          ",  & !  MAYBE DELETE, SET DEFAULT!?
      "ML_LNMDINT_RANDOM                         ", "ML_FF_LNMDINT_RANDOM                      ",  &
      "ML_LTEST                                  ", "ML_FF_LTEST                               ",  &
      "ML_LTRJ                                   ", "ML_FF_LTRJ                                ",  &
      "ML_LTOTEN_SYSTEM                          ", "ML_FF_LTOTEN_SYSTEM                       ",  &
      "ML_NTEST                                  ", "ML_FF_NTEST                               ",  &
      "ML_INVERSE_SOAP                           ", "ML_FF_INVERSE_SOAP_MB                     ",  &
      "ML_LSUPERVEC                              ", "ML_FF_LSUPERVEC_MB                        ",  &
      "ML_LMETRIC1                               ", "ML_FF_LMETRIC1_MB                         ",  &
      "ML_LMETRIC2                               ", "ML_FF_LMETRIC2_MB                         ",  &
      "ML_LVARTRAN1                              ", "ML_FF_LVARTRAN1_MB                        ",  &
      "ML_LVARTRAN2                              ", "ML_FF_LVARTRAN2_MB                        ",  &
      "ML_NMETRIC1                               ", "ML_FF_NMETRIC1_MB                         ",  &
      "ML_NMETRIC2                               ", "ML_FF_NMETRIC2_MB                         ",  &
      "ML_NVARTRAN1                              ", "ML_FF_NVARTRAN1_MB                        ",  &
      "ML_NVARTRAN2                              ", "ML_FF_NVARTRAN2_MB                        ",  &
      "ML_LWINDOW1                               ", "ML_FF_LWINDOW1_MB                         ",  &
      "ML_LWINDOW2                               ", "ML_FF_LWINDOW2_MB                         ",  &
      "ML_IWINDOW1                               ", "ML_FF_IWINDOW1_MB                         ",  &
      "ML_IWINDOW2                               ", "ML_FF_IWINDOW2_MB                         ",  &
      "ML_MSPL1                                  ", "ML_FF_MSPL1_MB                            ",  &
      "ML_MSPL2                                  ", "ML_FF_MSPL2_MB                            ",  &
      "ML_LNORM1                                 ", "ML_FF_LNORM1_MB                           ",  &
      "ML_LNORM2                                 ", "ML_FF_LNORM2_MB                           ",  &
      "ML_NHYP2                                  ", "ML_FF_NHYP2_MB                            ",  &
      "ML_NR1                                    ", "ML_FF_NR1_MB                              ",  &
      "ML_NR2                                    ", "ML_FF_NR2_MB                              ",  &
      "ML_RMETRIC1                               ", "ML_FF_RMETRIC1_MB                         ",  &
      "ML_RMETRIC2                               ", "ML_FF_RMETRIC2_MB                         ",  &
      "ML_RANDOM_SEED                            ", "ML_FF_RANDOM_SEED                         ",  &
      "ML_LDISCARD_STRUCTURES_NOT_GIVING_BASIS   ", "ML_FF_LDISCARD_STRUCTURES_NOT_GIVING_BASIS",  &
      "ML_ICUT1                                  ", "ML_FF_ICUT1_MB                            ",  &
      "ML_ICUT2                                  ", "ML_FF_ICUT2_MB                            ",  &
      "ML_NBLOCK                                 ", "ML_FF_NBLOCK_FFN                          "], &
      [nalt, ntags])
   !> Short description of each tag.
   character(len=*), parameter :: descriptions(ntags) = [                                                      &
      ! General, mode of operation:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "Main control tag to enable machine learning                                                          ", & ! ML_LMLFF
      "Machine learning operation mode                                                                      ", & ! ML_ISTART
      !" ML_FF_ISTART=0: Training from scratch. ML_FF_ISTART=1: Training continued from ML_AB file. ML_FF_ISTART=2: Force-field only. No learning. ML_FF_ISTART=3: Reselecting local configurations from existing ML_AB file.

      ! Descriptors:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "Cutoff radius of radial descriptors                                                                  ", & ! ML_RCUT1
      "Gaussian width for broadening the atomic distribution for radial descriptors                         ", & ! ML_SION1
      "Number of radial basis functions for atomic distribution for radial descriptors                      ", & ! ML_MRB1                                   
      "Cutoff radius of angular descriptors                                                                 ", & ! ML_RCUT2
      "Gaussian width for broadening the atomic distribution for angular descriptors                        ", & ! ML_SION2
      "Number of radial basis functions for atomic distribution for angular descriptors                     ", & ! ML_MRB2
      "Maximum angular momentum quantum number of spherical harmonics used to expand atomic distributions   ", & ! ML_LMAX2
      "Angular filtering enabled                                                                            ", & ! ML_LAFILT2
      "Angular filtering parameter a_FILT                                                                   ", & ! ML_AFILT2
      "Angular filtering function type                                                                      ", & ! ML_IAFILT2
      "Enable sparsification of angular descriptors                                                         ", & ! ML_LSPARSDES
      "Number of highest eigenvalues relevant in the sparsification algorithm of the angular descriptors    ", & ! ML_NRANK_SPARSDES
      "Desired ratio of selected to all descriptors resulting from the angular descriptor sparsification    ", & ! ML_RDES_SPARSDES

      ! Kernel:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "Weight of radial descriptors in the kernel (the angular counterpart is chosen so that the sum is 1.0)", & ! ML_W1
      "Power of the polynomial kernel                                                                       ", & ! ML_NHYP

      ! Controlling Bayesian criterion:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "Enable automatic updating of the Bayesian error estimation threshold during on-the-fly training      ", & ! ML_ICRITERIA
      "Decides whether update of threshold is done in the same MD step or the next MD step                  ", & ! ML_IUPDATE_CRITERIA
      "Read Bayesian error estimation threshold from ML_AB file                                             ", & ! ML_LCTIFOR_USE_FROM_FILE
      "Bayesian error estimation threshold (initial or static value depending on other settings)            ", & ! ML_CTIFOR
      "Scaling factor for ML_CTIFOR. The interval 0<ML_SCLC_CTIFOR<1 increases num. of local configurations ", & ! ML_SCLC_CTIFOR
      "Standard error parameter required for the automatic update of the Bayesian error estimation threshold", & ! ML_CSIG
      "Slope parameter required for the automatic update of the Bayesian error estimation threshold         ", & ! ML_CSLOPE
      "Additional parameter controlling the update of the Bayesian error estimation threshold               ", & ! ML_CX
      "Factor controlling the occurence of enforced first principle calculations                            ", & ! ML_CDOUB

      ! Sparsification and regression:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "Threshold for the CUR algorithm used in the sparsification of local reference configurations         ", & ! ML_EPS_LOW
      "Convergence criterion for the optimization of parameters within the Bayesian linear regression       ", & ! ML_EPS_REG
      "Linear regression algorithm                                                                          ", & ! ML_IALGO_LINREG
      ! ML_FF_IALGO_LINREG=1: Solve L2-norm (ridge regression). ML_IALGO_LINREG=2: Solve L1-norm with QR factorization. ML_INPUT_HANDLE%ML_FF_IALGO_LINREG = 3 : Solve L1-norm with truncated SVD. ML_INPUT_HANDLE%ML_FF_IALGO_LINREG = 4 : SVD + Tihonov regularization
      "Leverage scoring calculation mode (for sparsification of local configurations)                       ", & ! ML_ISVD
      ! ML_ISVD=1: The values are given from the eigenvalues provided by the diagonalization. ML_ISVD=2: These values are provided from the eigenvalues given from SVD
      "Time dependency of regularization parameters                                                         ", & ! ML_IREG
      ! ML_IREG_MB=1: The (initial) precision (ML_FF_SIGV0_MB) and noise (ML_FF_SIGW0_MB) parameters are kept constant. ML_IREG_MB=2: The parameters are optimized
      "Initial regularization parameter (noise, reversed and squared)                                       ", & ! ML_SIGV0
      "Initial regularization parameter (precision, reversed and squared)                                   ", & ! ML_SIGW0

      ! Weighting and sampling:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "Scaling weight for total energies in the training data                                               ", & ! ML_WTOTEN
      "Scaling weight for forces in the training data                                                       ", & ! ML_WTIFOR
      "Scaling weight for stresses in the training data                                                     ", & ! ML_WTISIF
      "Number of estimated errors stored to determine the threshold for the Bayesian error                  ", & ! ML_MHIS
      "Minimum number of MD steps between potential training samples                                        ", & ! ML_NMDINT
      "Weighting method for energies, forces and stresses                                                   ", & ! ML_IWEIGHT
      ! ML_IWEIGHT=1: The unnormalized energies, forces and stress tensor training data are divided by the weights determined by the flags ML_WTOTEN (eV/atom), ML_WTIFOR (eV/Angstrom) and ML_WTSIF (kBar), respectively.. ML_IWEIGHT=2: The training data are normalized by using their standard deviations. Then, the normalized energy, forces and stress tensor are multiplied by unitless ML_WTOTEN, ML_WTIFOR and ML_WTSIF, respectively. ML_IWEIGHT=3: Same as ML_FF_IWEIGHT=2 but the training data is divided into individual subsets.

      ! Sizes of static arrays:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!
      "Maximum number of local configurations stored in memory                                              ", & ! ML_MB
      ! This number is needed for allocation of static arrays, especially the number of columns in the design matrix.
      "Maximum number of structures stored for training                                                     ", & ! ML_MCONF
      ! This number is needed for allocation of static arrays, especially the number of rows in the design matrix.
      "Maximum number of configurations stored temporarily as training candidates                           ", & ! ML_MCONF_NEW

      ! Special features:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "Enable heat flux calculation (output written to ML_HEAT)                                             ", & ! ML_LHEAT
      "Enable thermodynamic integration (t.i.)                                                              ", & ! ML_LCOUPLE
      "Number of atoms whose interaction is controlled by the t.i. coupling parameter                       ", & ! ML_NATOM_COUPLED
      "List of atoms whose interaction is controlled by the t.i. coupling parameter                         ", & ! ML_ICOUPLE
      "Thermodynamic integration coupling parameter                                                         ", & ! ML_RCOUPLE

      ! Reference energies:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "Scaling mode of energies                                                                             ", & ! ML_ISCALE_TOTEN
      ! Specifies how to scale the energy data. ML_ISCALE_TOTEN=1: The total energy is scaled to the total energy of the isolated atoms given by ML_EATOM. ML_ISCALE_TOTEN=2: The total energy is scaled to the average of the training data.
      "List of reference energies of isolated atoms for each species in the system                          ", & ! ML_EATOM_REF
      ! Output options:
      "Enable output of kin. and pot. energy for each atom at each MD time step (written to ML_EATOM)       ", & ! ML_LEATOM

      ! Not important, description not used and not yet reviewed:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      "Specifies whether broadening of the atomic distribution is used for the radial descriptor            ", & ! ML_IBROAD1
      ! ML_IBROAD1=1: No broadening. ML_IBROAD1=2: Element-independent broadening.
      "Specifies whether broadening of the atomic distribution is used for the angular descriptor           ", & ! ML_IBROAD2
      ! ML_IBROAD2=1: No broadening. ML_IBROAD2=2: Element-independent broadening.
      "Specifies the number of machine learning instances                                                   ", & ! ML_TOTNUM_INSTANCES
      "Specifies the starting value for the filehandle numbers within Fortran                               ", & ! ML_FILESTREAM_START
      "Specifies the type of matrix for sparsification                                                      ", & ! ML_IMAT_SPARS
      ! ML_IMAT_SPARS=1: Covariance matrix. ML_IMAT_SPARS=2: SOAP matrix.
      "Discard local configurations if the maximum number in storage is exceeded                            ", & ! ML_LBASIS_DISCARD
      !Specifies whether the calculation is continued if number of local configurations (NB) is larger than ML_MB and basis sets are thrown away to satisfy NB<=ML_NB", & ! ML_LBASIS_DISCARD
      !"ML_LCONF_DISCARD                          ", Specifies whether the calculation is continued if number of training structures (NCONF) is larger than ML_MCONF and training structures are thrown away to satisfy NCONF<=ML_MCONF, or not.
      "Discard training structures if the maximum number in storage is exceeded                             ", & ! ML_LCONF_DISCARD
      ! Specifies whether the calculation is continued if number of training structures (NCONF) is larger than ML_MCONF and training structures are thrown away to satisfy NCONF<=ML_MCONF, or not.
      "Turns many body interaction on. Code won't work without this tag                                     ", & ! ML_LMLMB
      "Decides whether error estimations are turned on or not                                               ", & ! ML_LERR
      ! Turning off can possibly accelerate the MD when only force fields are used.
      "Specifies the dimension of scalapack grid. Should never be touched                                   ", & ! ML_NDIM_SCALAPACK
      "Deprecated, in past it specified contents of the ML_LOGFILE                                          ", & ! ML_NWRITE
      "Decides wether ML_NMDINT times a random number between 0 and 1 is used instead of ML_NMDINT          ", & ! ML_LNMDINT_RANDOM
      "Whether ab initio calculations are executed at every MD step regardless of the estimated error       ", & ! ML_LTEST
      "DESCRIPTION OF ML_LTRJ                                                                               ", & ! ML_LTRJ
      ! ML_LTRJ=.TRUE.: MD trajectory is always obtained using ab initio forces and stress tensors. ML_LTRJ=.FALSE.: MD trajectory is obtained using generated force field when the estimated error is small. Otherwise, ab initio forces and stress tensors are used.
      "DESCRIPTION OF ML_LTOTEN_SYSTEM                                                                      ", & !ML_LTOTEN_SYSTEM
      ! ML_LTOTEN_SYSTEM=.TRUE.:  Total energy of system is used in regression. ML_LTOTEN_SYSTEM=.FALSE.: Energy per atom is used in regression.
      "This parameter determines how often the test ab initio calculation will be executed                  ", & ! ML_NTEST
      ! The test data will be obtained at every ML_NTEST MD step.
      "Decides how the inverse matrix is calculated for ridge regression                                    ", & ! ML_INVERSE_SOAP
      ! ML_INVERSE_SOAP=1: Diagonalization. ML_INVERSE_SOAP=2: LU factorization.
      "Specifies whether super vector is used for Kernel or not                                             ", & ! ML_LSUPERVEC
      "Specifies whether metric function is used for radial descriptor or not                               ", & ! ML_LMETRIC1
      "Specifies whether metric function is used for angular descriptor or not                              ", & ! ML_LMETRIC2
      "Specifies whether variable transform is done for radial coordinate of radial descriptor or not       ", & ! ML_LVARTRAN1
      "Specifies whether variable transform is done for radial coordinate of angular descriptor or not      ", & ! ML_LVARTRAN2
      "Specifies the polynomial parameter for the metric function for the radial descriptor                 ", & ! ML_NMETRIC1
      "Specifies the polynomial parameter for the metric function for the angular descriptor                ", & ! ML_NMETRIC2
      "Polynomial parameter for the variable transform of the radial coordinate for the radial descriptor   ", & ! ML_NVARTRAN1
      "Polynomial parameter for the variable transform of the radial coordinate for the angular descriptor  ", & ! ML_NVARTRAN2
      "Specifies wheter a window function is used for the radial descriptor or not                          ", & ! ML_LWINDOW1
      "Specifies wheter a window function is used for the angular descriptor or not                         ", & ! ML_LWINDOW2
      "Specifies the type of window function for radial descriptor if ML_LWINDOW1=.TRUE.                    ", & ! ML_IWINDOW1
      ! ML_IWINDOW1=1: Hann window. ML_IWINDOW1=2: Hamming window. ML_IWINDOW1=3: Blackman window. ML_IWINDOW1=4: Nuttall window. ML_IWINDOW1=5: Blackman-Nutall. ML_IWINDOW1=6: Blackman-Harris.
      "Specifies the type of window function for angular descriptor if ML_LWINDOW2=.TRUE.                   ", & ! ML_IWINDOW2
      ! ML_IWINDOW2=1: Hann window. ML_IWINDOW2=2: Hamming window. ML_IWINDOW2=3: Blackman window. ML_IWINDOW2=4: Nuttall window. ML_IWINDOW2=5: Blackman-Nutall. ML_IWINDOW2=6: Blackman-Harris.
      "Maximum number of radial grid points used in spline-interpolation for radial descriptor              ", & ! ML_MSPL1
      "Maximum number of radial grid points used in spline-interpolation for angular descriptor             ", & ! ML_MSPL2
      "Specifies whether normalization is taken into account for radial descriptor                          ", & ! ML_LNORM1
      "Specifies whether normalization is taken into account for angular descriptor                         ", & ! ML_LNORM2
      "Sets the hyper parameter (power) of SOAP                                                             ", & ! ML_NHYP2
      "Number of radial mesh points for calculating radial descriptor                                       ", & ! ML_NR1
      "Number of radial mesh points for calculating angular descriptor                                      ", & ! ML_NR2
      "Parameter used for metric function for radial descriptor                                             ", & ! ML_RMETRIC1
      "Parameter used for metric function for angular descriptor                                            ", & ! ML_RMETRIC2
      "Sets the random seed for ML_LNMDINT_RANDOM                                                           ", & ! ML_RANDOM_SEED
      "Specifies whether structures not giving new local configurations are discarded for ML_ISTART=3       ", & ! ML_LDISCARD_STRUCTURES_NOT_GIVING_BASIS
      ! Specifies whether structures not giving new local configurations are thrown away as trainings structures for ML_ISTART=3
      "Specifies the type of cut off function for the radial descriptor                                     ", & ! ML_ICUT1
      ! ML_ICUT1=1: Behler-Parinello type cutoff function.. ML_ICUT1=2: Miwa-Ohno type cutoff function.
      "Specifies the type of cut off function for the angular descriptor                                    ", & ! ML_ICUT2
      ! ML_ICUT2=1: Behler-Parinello type cutoff function.. ML_ICUT2=2: Miwa-Ohno type cutoff function.
      "DESCRIPTION OF ML_NBLOCK                                                                             "]   ! ML_NBLOCK
   !> Maximum number of characters for tags.
   integer, parameter                  :: tag_width = len(tags)
   !> Maximum number of characters for tag description.
   integer, parameter                  :: description_width = len(descriptions)
   !> Tag state symbol width.
   integer, parameter                  :: state_symbol_width = 3

   contains

   subroutine taglist_init(this)

      class(taglist), intent(inout) :: this

      ! Initialize state array to TL_UNSET.
      allocate(this%state(ntags))
      this%state = TL_UNSET
      ! Initialize alternative index to 1 (preferred tag).
      allocate(this%alt_index(ntags))
      this%alt_index = 1

   end subroutine taglist_init

   function taglist_get_state(this, tag) result(state)

      class(taglist),   intent(in) :: this
      character(len=*), intent(in) :: tag
      integer(kind(TL_STATE))      :: state

      state = this%state(taglist_find_index(tag))

   end function taglist_get_state

   subroutine taglist_set_state(this, tag, in_state, alt_index)

      class(taglist),          intent(inout)          :: this
      character(len=*),        intent(in)             :: tag
      integer(kind(TL_STATE)), intent(in)             :: in_state
      integer,                 intent(in),   optional :: alt_index

      integer :: t

      t = taglist_find_index(tag)
      this%state(t) = in_state
      if (present(alt_index)) this%alt_index(t) = alt_index

   end subroutine taglist_set_state

   subroutine taglist_get_deprecated(this, list)

      class(taglist),           intent(in)                 :: this
      character(len=tag_width), intent(inout), allocatable :: list(:,:)

      integer :: i
      integer :: c

      ! Count how many deprecated tags were found.
      c = 0
      do i = 1, ntags
         if (this%alt_index(i) > 1) c = c + 1
      end do

      if (allocated(list)) deallocate(list)
      allocate(character(len=tag_width) :: list(c, 2))

      ! Fill list with (deprecated, preferred) tags.
      c = 1
      do i = 1, ntags
         if (this%alt_index(i) > 1) then
            list(c, 1) = tags(this%alt_index(i), i)
            list(c, 2) = tags(1, i)
            c = c + 1
         end if
      end do

   end subroutine taglist_get_deprecated

   function taglist_state_symbol_from_tag(this, tag) result(symbol)

      class(taglist),   intent(in)      :: this
      character(len=*), intent(in)      :: tag
      character(len=state_symbol_width) :: symbol

      symbol = taglist_state_symbol(this%get_state(tag))

   end function taglist_state_symbol_from_tag

   !!!!!!!!!!!!!!!!!!!!!!!!!
   ! Non-member procedures !
   !!!!!!!!!!!!!!!!!!!!!!!!!

   pure function taglist_get_tag_width() result(width)

      integer :: width

      width = tag_width

   end function taglist_get_tag_width

   function taglist_description(tag) result(description)

      character(len=*), intent(in)  :: tag
      character(len=:), allocatable :: description

      description = descriptions(taglist_find_index(tag))

   end function taglist_description

   function taglist_find_index(tag) result(tag_index)

      use iso_fortran_env, only: error_unit
      character(len=*), intent(in) :: tag
      integer(kind(TL_STATE))      :: tag_index

      integer :: t ! Tag loop index.
      integer :: a ! Alternative tag loop index.

      tag_loop: &
      do t = 1, ntags
         do a = 1, nalt
            ! Avoid comparing with potential empty alternative string, as this
            ! is asking for trouble...
            if (a > 1 .and. adjustl(trim(tags(a, t))) == "") cycle
            if (adjustl(trim(tag)) == adjustl(trim(tags(a, t)))) exit tag_loop
         end do
      end do tag_loop
      if (t >= ntags + 1) then
         write(error_unit, '(A,1X,A)') "ERROR: Unknown tag: ", tag
         stop 1
      end if
      tag_index = t

   end function taglist_find_index

   function taglist_get_all_tags() result(list)

      character(len=:), allocatable :: list(:)

      list = tags(1, :)

   end function taglist_get_all_tags

   function taglist_get_alt_tags(tag) result(list)

      character(len=*)              :: tag
      character(len=:), allocatable :: list(:)

      integer :: t ! Tag loop index.
      integer :: a ! Alternative tag loop index.
      integer :: n ! Number of non-empty alternative strings.
      integer :: c ! Counter.

      t = taglist_find_index(tag)
      n = 0
      do a = 1, nalt
         if (.not. (adjustl(trim(tags(a, t))) == "")) n = n + 1
      end do
      allocate(character(len=tag_width) :: list(n))
      c = 1
      do a = 1, nalt
         if (.not. (adjustl(trim(tags(a, t))) == "")) then
            list(c) = tags(a, t)
            c = c + 1
         end if
      end do

   end function taglist_get_alt_tags

   pure function taglist_get_description_width() result(width)

      integer :: width

      width = description_width

   end function taglist_get_description_width

   pure function taglist_state_symbol(state) result(symbol)

      integer(kind(TL_STATE)), intent(in) :: state
      character(len=state_symbol_width)   :: symbol

      select case(state)
         case(TL_UNSET)
            symbol = "(?)"
         case(TL_DEFAULT)
            symbol = "   "
         case(TL_INCAR)
            symbol = "(I)"
         case(TL_INCAR_ALT)
            symbol = "(i)"
         case(TL_OVERRIDE)
            symbol = "(!)"
      end select

   end function taglist_state_symbol

end module ml_ff_taglist
