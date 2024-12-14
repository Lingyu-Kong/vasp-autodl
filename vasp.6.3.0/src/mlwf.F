!#define debug
#include "symbol.inc"
!> @brief Module to produce maximally localized Wannier functions
MODULE mlwf
   USE prec
   USE mpimy
   USE locproj, ONLY : LPRJ_function
   IMPLICIT NONE

   !> Different available modes for obtaining intial projections (AMN)
   ENUM, BIND(C)
      ENUMERATOR :: PROJ_UNKNOWN, PROJ_SCDM, PROJ_WANNIER90, PROJ_LOCPROJ, PROJ_LOCPROJAUTO
   END ENUM

   !> Different available modes for obtainin unitary transformations (UMN)
   ENUM, BIND(C)
      ENUMERATOR :: TRAFO_UNKNOWN, TRAFO_ONESHOTSVD, TRAFO_WANNIER90_MODE, TRAFO_UMCO, TRAFO_SCDM
   END ENUM

   !> Different available modes for the smooth band cutoff
   ENUM, BIND(C)
      ENUMERATOR :: CUTOFF_NONE, CUTOFF_ERFC, CUTOFF_GAUSSIAN, CUTOFF_FERMI, CUTOFF_NUM_WANN
   END ENUM

   !> Type containing the information to compute maximally localized wannier functions.
   TYPE MLWF_type
      INTEGER :: ISPIN     !< Number of independent spin channels
      INTEGER :: NRSPINORS !< NUmber of spinor components
      INTEGER :: NB_TOT    !< Total number of bands
      INTEGER :: NUM_WANN  !< Number of desired wannier functions
      INTEGER :: NKPTS     !< Number of k-points for which the projections will be computed
      INTEGER :: NEXCLB    !< Number of bands to exclude
      CHARACTER(LEN=40) :: SZNAM1 !< Label read from the INCAR file
      CHARACTER(len=9) :: SEED_NAME='wannier90' !< name of the files for wannier90
      LOGICAL, ALLOCATABLE :: EXCLUDE_BAND(:) !< Bands to exclude
      TYPE(LPRJ_function), ALLOCATABLE :: LPRJ_functions(:) !< List of functions determining the local projections
      TYPE(communic) :: COMM
      INTEGER :: PROJ_MODE = PROJ_UNKNOWN !< Mode to use to determine the projections
      INTEGER :: TRAFO_MODE = TRAFO_UNKNOWN !< Mode to use to determine unitary wannier transformation
      ! related to the k-point grid
      INTEGER :: mp_grid(3)
      INTEGER, ALLOCATABLE :: NEQUIV(:)
      REAL(q), ALLOCATABLE :: kpt_latt(:,:)
      ! related to the structure
      REAL(q) :: real_lattice(3,3)
      REAL(q) :: recip_lattice(3,3)
      INTEGER :: num_atoms
      REAL(q), ALLOCATABLE :: atoms_cart(:,:)
#ifdef gammareal
      LOGICAL :: gamma_only=.TRUE.
#else
      LOGICAL :: gamma_only=.FALSE.
#endif
      CHARACTER(len=20), ALLOCATABLE :: atom_symbols(:)
      !> Contains the projections of the Kohn-Sham orbitals onto i
      !> localized orbitals < \tilde \psi_{nb,k} | S | RYlm_i >.
      !> The dimensions are (NBANDS,NK,ISP,NPROJ)
      GDEF, ALLOCATABLE :: LPRJ_COVL(:,:,:,:)
      !> Contains the projections of the Kohn-Sham orbitals onto IAOs
      !> The dimensions are (NBANDS,NK,ISP,NPROJ)
      GDEF, ALLOCATABLE :: LPRJ_COVL_IAO(:,:,:,:)
      ! disentanglement windows
      REAL(q) :: DIS_FROZ_MAX
      REAL(q) :: DIS_FROZ_MIN
      REAL(q) :: DIS_WIN_MAX
      REAL(q) :: DIS_WIN_MIN
      ! related to MMN file
      INTEGER :: NNTOT
      INTEGER, ALLOCATABLE :: NNLIST(:,:)
      INTEGER, ALLOCATABLE :: NNCELL(:,:,:)
      ! Input to wannier90
      COMPLEX(q), ALLOCATABLE :: M_matrix(:,:,:,:,:)
      COMPLEX(q), ALLOCATABLE :: A_matrix(:,:,:,:)
      REAL(q), ALLOCATABLE :: eigenvalues(:,:,:)
      ! wannier90_run variables: output
      COMPLEX(q), ALLOCATABLE :: U_matrix(:,:,:,:)
      COMPLEX(q), ALLOCATABLE :: U_matrix_opt(:,:,:,:)
      LOGICAL, ALLOCATABLE :: lwindow(:,:,:)
      LOGICAL, ALLOCATABLE :: lexclude_band(:)
      REAL(q), ALLOCATABLE :: wann_centres(:,:,:)
      REAL(q), ALLOCATABLE :: wann_spreads(:,:)
      REAL(q), ALLOCATABLE :: spread(:,:)
   END TYPE MLWF_type

   ! A global instance of maximally localized wavefunctions
   ! This is kept for compatiblity with old approach of defining module variables
   ! This type of usage should be discontinued in favour of explicitly passing the
   ! instance to the funcitons that need it
   TYPE(MLWF_type), TARGET, SAVE :: MLWF_GLOBAL

   ! write the wave functions to UNK files
   LOGICAL, PRIVATE, SAVE :: WRITE_UNK

   !> write the bloch functions as XSF file
   LOGICAL, PRIVATE, SAVE :: LWRITE_BLOCH_XSF=.FALSE.
   !> select which bloch functions to write as a XSF file
   INTEGER, ALLOCATABLE :: WRITE_BLOCH_XSF_IW(:)

   !> write the IAO functions as a wavecar file
   LOGICAL, PUBLIC, SAVE :: LWAVE_WANNIER = .FALSE.
   !> write the wannier functions as XSF file
   LOGICAL, PRIVATE, SAVE :: LWRITE_WANNIER_XSF = .FALSE.
   !> select which wannier functions to write as a XSF file
   INTEGER, ALLOCATABLE :: WRITE_WANNIER_XSF_IW(:)

   !> Write the LOCPROJ functions as a wavecar file
   LOGICAL, PUBLIC, SAVE :: LWAVE_LOCPROJ = .FALSE.
   !> write the locproj functions as XSF file
   LOGICAL, PRIVATE, SAVE :: LWRITE_LOCPROJ_XSF = .FALSE.
   ! select which locproj functions to write as a XSF file
   INTEGER, ALLOCATABLE :: WRITE_LOCPROJ_XSF_IW(:)

   !> Write the IAOS as a wavecar file
   LOGICAL, PRIVATE, SAVE :: LWAVE_IAO = .FALSE.
   !> write the IAO functions as XSF file
   LOGICAL, PRIVATE, SAVE :: LWRITE_IAO_XSF = .FALSE.
   !> select which IAO functions to write as a XSF file
   INTEGER, ALLOCATABLE :: WRITE_IAO_XSF_IW(:)

   !> read amn file instead of computing it
   LOGICAL, PRIVATE, SAVE :: READ_AMN

   ! Number of Wannier orbitals
   INTEGER, SAVE :: NUM_WANN = 0

   ! prepare the input for wannier90
   LOGICAL, SAVE :: LWANNIER90=.FALSE.
   ! setup and run wannier90 in library mode
   LOGICAL, SAVE :: LWANNIER90_RUN=.FALSE.
   ! guess the disentanglement window
   LOGICAL, SAVE :: LWANNIER90_AUTO_WINDOW=.FALSE.
   !> write U_nm(k) matrix
   LOGICAL, SAVE :: LWRITE_WANPROJ
   !> write the mmn and amn files when WANNIER90 runs in lib mode
   LOGICAL, SAVE :: WRITE_MMN_AMN

   ! related to the wannier input file
   LOGICAL,SAVE :: WANNIER90_WIN_LEGACY_MODE=.FALSE.
   INTEGER, PARAMETER :: WANNIER90_WIN_MAXLEN=8192
   CHARACTER(LEN=WANNIER90_WIN_MAXLEN),SAVE :: WANNIER90_WIN

   ! SCDM input parameters
   LOGICAL, SAVE :: LSCDM = .FALSE.

   !> Type containing the information for the cutoff
   !> (used by the SCDM method and with projections from LOCPROJ)
   TYPE cutoff_type
      INTEGER :: ISPIN
      INTEGER :: TYPE = CUTOFF_NONE !< Enum for the cutoff mode
      INTEGER :: NUM_WANN
      !> Parameter for center or erfc or gaussian for each spin channel
      REAL(q) :: MU(2)
      !> Parameter for broadening of erfc or gaussian for each spin channel
      REAL(q) :: SIGMA(2)
   END TYPE cutoff_type

   TYPE(cutoff_type),SAVE :: CUTOFF_DES

#ifdef VASP2WANNIER90
   INTERFACE
     subroutine wannier_setup( &
       seed_name_, mp_grid_, num_kpts_, real_lattice_, recip_lattice_, &
       kpt_latt_, num_bands_tot_, num_atoms_, atom_symbols_, atoms_cart_, &
       gamma_only_, spinors_, nntot_, nnlist_, nncell_, num_bands_, &
       num_wann_, proj_site_, proj_l_, proj_m_, proj_radial_, proj_z_, &
       proj_x_, proj_zona_, exclude_bands_ &
       , proj_s_, proj_s_qaxis_ &
       )
       implicit none
       integer, parameter :: dp = kind(1.0d0) ! statt 'USE w90_constants'
       integer, parameter :: num_nnmax=12     ! statt 'USE w90_parameters'
       character(len=*), intent(in) :: seed_name_
       integer, dimension(3), intent(in) :: mp_grid_
       integer, intent(in) :: num_kpts_
       real(kind=dp), dimension(3,3), intent(in) :: real_lattice_
       real(kind=dp), dimension(3,3), intent(in) :: recip_lattice_
       real(kind=dp), dimension(3,num_kpts_), intent(in) :: kpt_latt_
       integer, intent(in) :: num_bands_tot_
       integer, intent(in) :: num_atoms_
       character(len=*), dimension(num_atoms_), intent(in) :: atom_symbols_
       real(kind=dp), dimension(3,num_atoms_), intent(in) :: atoms_cart_
       logical, intent(in) :: gamma_only_
       logical, intent(in) :: spinors_
       integer, intent(out) :: nntot_
       integer, dimension(num_kpts_,num_nnmax), intent(out) :: nnlist_
       integer,dimension(3,num_kpts_,num_nnmax), intent(out) :: nncell_
       integer, intent(out) :: num_bands_
       integer, intent(out) :: num_wann_
       real(kind=dp), dimension(3,num_bands_tot_), intent(out) :: proj_site_
       integer, dimension(num_bands_tot_), intent(out) :: proj_l_
       integer, dimension(num_bands_tot_), intent(out) :: proj_m_
       integer, dimension(num_bands_tot_), intent(out) :: proj_radial_
       real(kind=dp), dimension(3,num_bands_tot_), intent(out) :: proj_z_
       real(kind=dp), dimension(3,num_bands_tot_), intent(out) :: proj_x_
       real(kind=dp), dimension(num_bands_tot_), intent(out) :: proj_zona_
       integer, dimension(num_bands_tot_), intent(out) :: exclude_bands_
       integer, dimension(num_bands_tot_), optional, intent(out) :: proj_s_
       real(kind=dp), dimension(3,num_bands_tot_), optional, intent(out) :: proj_s_qaxis_
     end subroutine wannier_setup
     !
     subroutine wannier_run( &
       seed_name_, mp_grid_, num_kpts_, real_lattice_, recip_lattice_, &
       kpt_latt_, num_bands_, num_wann_, nntot_, num_atoms_, atom_symbols_, &
       atoms_cart_, gamma_only_, M_matrix_, A_matrix_, eigenvalues_, &
       U_matrix_, U_matrix_opt_, lwindow_, wann_centres_, wann_spreads_, spread_)
       implicit none
       integer, parameter :: dp = kind(1.0d0) ! statt 'USE w90_constants'
       character(len=*), intent(in) :: seed_name_
       integer, dimension(3), intent(in) :: mp_grid_
       integer, intent(in) :: num_kpts_
       real(kind=dp), dimension(3,3), intent(in) :: real_lattice_
       real(kind=dp), dimension(3,3), intent(in) :: recip_lattice_
       real(kind=dp), dimension(3,num_kpts_), intent(in) :: kpt_latt_
       integer, intent(in) :: num_bands_
       integer, intent(in) :: num_wann_
       integer, intent(in) :: nntot_
       integer, intent(in) :: num_atoms_
       character(len=*), dimension(num_atoms_), intent(in) :: atom_symbols_
       real(kind=dp), dimension(3,num_atoms_), intent(in) :: atoms_cart_
       logical, intent(in) :: gamma_only_
       complex(kind=dp), dimension(num_bands_,num_bands_,nntot_,num_kpts_), intent(in) :: M_matrix_
       complex(kind=dp), dimension(num_bands_,num_wann_,num_kpts_), intent(in) :: A_matrix_
       real(kind=dp), dimension(num_bands_,num_kpts_), intent(in) :: eigenvalues_
       complex(kind=dp), dimension(num_wann_,num_wann_,num_kpts_), intent(out) :: U_matrix_
       complex(kind=dp), dimension(num_bands_,num_wann_,num_kpts_), optional, intent(out) :: U_matrix_opt_
       logical, dimension(num_bands_,num_kpts_), optional, intent(out) :: lwindow_
       real(kind=dp), dimension(3,num_wann_), optional, intent(out) :: wann_centres_
       real(kind=dp), dimension(num_wann_), optional, intent(out) :: wann_spreads_
       real(kind=dp), dimension(3), optional, intent(out) :: spread_
     end subroutine wannier_run
   END INTERFACE
#endif

   PUBLIC :: MLWF_READER
   PUBLIC :: MLWF_MAIN
   !PUBLIC :: MLWF_ROTATE_ORBITALS
   !PUBLIC :: MLWF_ROTATE_ORBITALS_FULLK
   !PUBLIC :: MLWF_ROTATE_ORBITALS_NOSYMM
   !PUBLIC :: MLWF_ROTATE_ORBITALS_NOSYMM_U
   !PUBLIC :: MLWF_DPSI_DK
   !PUBLIC :: MLWF_DHAM_DK
   PUBLIC :: MLWF_GET_U
   PUBLIC :: MLWF_ROTATE_HAM
   !
   PUBLIC :: GET_WAVE_FUNCTIONS
   PUBLIC :: GET_WAVE_FUNCTIONS_PHASE
   PUBLIC :: READ_WANPROJ_FILE

   CONTAINS

!******************** SUBROUTINE MLWF_READER ***************************
!> @brief Read variables related to producing Maximally localized Wannier functions
!***********************************************************************
 SUBROUTINE MLWF_READER(INFO,IU5,IU6,IU0)
   USE base
   USE vaspxml
   USE full_kpoints
   USE wave_high
   USE reader_tags
   USE main_mpi
   USE pead, ONLY : PEAD_REQUEST
   USE string, ONLY : LOWERCASE

   TYPE(info_struct) INFO
   INTEGER IU5,IU6,IU0
   ! local variables
   REAL(q) :: MU(2), SIGMA(2)
   INTEGER :: CUTOFF_TYPE
   INTEGER IPOS
   INTEGER IDUM, N, IERR
   REAL(q) RDUM
   LOGICAL LOPEN,LDUM
   CHARACTER (LEN=20) :: CUTOFF_TYPE_STR
   CHARACTER (LEN=255) :: BUFLINE
   LOGICAL  :: LWIN_FOUND
   INTEGER  :: IWIN_SIZE
   LOGICAL  :: NOCCUR
   EXTERNAL :: NOCCUR

   CALL OPEN_INCAR_IF_FOUND(IU5, LOPEN)

   CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LWANNIER90', LWANNIER90, IERR, WRITEXMLINCAR)

   CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LWANNIER90_RUN', LWANNIER90_RUN, IERR, WRITEXMLINCAR)
#ifndef VASP2WANNIER90
   IF (WANNIER90()) CALL vtutor%error('VASP was compiled without the VASP2WANNIER90 interface.')
#endif

   CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LWANNIER90_AUTO_WINDOW', LWANNIER90_AUTO_WINDOW, IERR, WRITEXMLINCAR)

   ! read wannier90 input file
   WANNIER90_WIN = ""
   CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'WANNIER90_WIN', WANNIER90_WIN, IDUM, IERR, WRITEXMLINCAR)
   ! remove " from the input
   DO N=1,LEN(WANNIER90_WIN)-1
     IF (WANNIER90_WIN(N:N) == '"') WANNIER90_WIN(N:N) = " "
   ENDDO
   ! if wannier90_win is not present then we fallback to legacy mode were the wannier library
   ! reads wannier90.win from the hard drive
   IF (WANNIER90().AND.IERR/=0) THEN
     INQUIRE(FILE=DIR_APP(1:DIR_LEN)//'wannier90.win',EXIST=LWIN_FOUND,SIZE=IWIN_SIZE)
     IF (LWIN_FOUND.AND.IWIN_SIZE>0) THEN
        WANNIER90_WIN_LEGACY_MODE = .TRUE.
        CALL vtutor%advice('In the current version of VASP you can specify the input to wannier90 &
                            in the INCAR file. To do so use the WANNIER90_WIN tag which is &
                            a string delimited by ", that that can span mutiple lines. &
                            VASP will then write the contents of WANNIER90_WIN to wannier90.win &
                            and automatically set some variables to be consistent with the present run. &
                            Presently you have not specified WANNIER90_WIN, therefore the interface &
                            to wannier90 will run in legacy mode. This means a barebones wannier90.win &
                            will be created only if this file is not already present.')
     ENDIF
   ENDIF

   ! Do we want to write UWAN files?
   LWRITE_WANPROJ=.FALSE.
   CALL PROCESS_INCAR(.FALSE., IU0, IU5, 'LWRITE_WANPROJ', LWRITE_WANPROJ, IERR, WRITEXMLINCAR)

   ! Number of Wannier functions
   CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NUM_WANN', NUM_WANN, IERR, WRITEXMLINCAR)

   ! Do we want to write UNK files?
   WRITE_UNK=.FALSE.
   CALL PROCESS_INCAR(.FALSE., IU0, IU5, 'LWRITE_UNK', WRITE_UNK, IERR, WRITEXMLINCAR)

   ! Do we want to write the bloch functions as XSF files?
   CALL PROCESS_INCAR(.FALSE., IU0, IU5, 'LWRITE_BLOCH_XSF', LWRITE_BLOCH_XSF, IERR, WRITEXMLINCAR)
   ! Select which bloch functions to write as XSF file
   CALL PROCESS_INCAR(.FALSE., IU0, IU5, 'WRITE_BLOCH_XSF_IW', IDUM, IERR, FOUNDNUMBER=N)
   IF (N/=0) THEN
      ALLOCATE(WRITE_BLOCH_XSF_IW(N))
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'WRITE_BLOCH_XSF_IW', WRITE_BLOCH_XSF_IW, N, IERR, WRITEXMLINCAR)
   ENDIF

   ! Write wavecar with WANNIER orbitals?
   CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LWAVE_WANNIER', LWAVE_WANNIER, IERR)
   ! Do we want to write the wannier functions as XSF files?
   CALL PROCESS_INCAR(.FALSE., IU0, IU5, 'LWRITE_WANNIER_XSF', LWRITE_WANNIER_XSF, IERR, WRITEXMLINCAR)
   ! Select which wannier functions to write as XSF file
   CALL PROCESS_INCAR(.FALSE., IU0, IU5, 'WRITE_WANNIER_XSF_IW', IDUM, IERR, FOUNDNUMBER=N)
   IF (N/=0) THEN
      ALLOCATE(WRITE_WANNIER_XSF_IW(N))
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'WRITE_WANNIER_XSF_IW', WRITE_WANNIER_XSF_IW, N, IERR, WRITEXMLINCAR)
   ENDIF

   ! Write wavecar with LOCPROJ orbitals?
   CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LWAVE_LOCPROJ', LWAVE_LOCPROJ, IERR)
   ! Do we want to write the locproj functions as XSF files?
   CALL PROCESS_INCAR(.FALSE., IU0, IU5, 'LWRITE_LOCPROJ_XSF', LWRITE_LOCPROJ_XSF, IERR, WRITEXMLINCAR)
   ! Select which LOCPROJ functions to write as XSF file
   CALL PROCESS_INCAR(.FALSE., IU0, IU5, 'WRITE_LOCPROJ_XSF_IW', IDUM, IERR, FOUNDNUMBER=N)
   IF (N/=0) THEN
      ALLOCATE(WRITE_LOCPROJ_XSF_IW(N))
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'WRITE_LOCPROJ_XSF_IW', WRITE_LOCPROJ_XSF_IW, N, IERR, WRITEXMLINCAR)
   ENDIF

   ! Write wavecar with IAO orbitals?
   CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LWAVE_IAO', LWAVE_IAO, IERR)
   ! Do we want to write the IAO functions as XSF files?
   CALL PROCESS_INCAR(.FALSE., IU0, IU5, 'LWRITE_IAO_XSF', LWRITE_IAO_XSF, IERR, WRITEXMLINCAR)
   ! Select which IAO functions to write as XSF file
   CALL PROCESS_INCAR(.FALSE., IU0, IU5, 'WRITE_IAO_XSF_IW', IDUM, IERR, FOUNDNUMBER=N)
   IF (N/=0) THEN
      ALLOCATE(WRITE_IAO_XSF_IW(N))
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'WRITE_IAO_XSF_IW', WRITE_IAO_XSF_IW, N, IERR, WRITEXMLINCAR)
   ENDIF

   ! Do we want to write the mmn and amn files even when WANNIER90 runs in library mode?
   WRITE_MMN_AMN=WANNIER90()
   CALL PROCESS_INCAR(.FALSE., IU0, IU5, 'LWRITE_MMN_AMN', WRITE_MMN_AMN, IERR, WRITEXMLINCAR)

   ! Do we want to read the amn file instead of
   ! computing it anew?
   READ_AMN=.FALSE.
   CALL PROCESS_INCAR(.FALSE., IU0, IU5, 'LREAD_AMN', READ_AMN, IERR, WRITEXMLINCAR)

   ! Use SCDM method?
   CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSCDM', LSCDM, IERR, WRITEXMLINCAR)
   ! Location of the cut-off function
   MU = HUGE(0._q)
   CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'CUTOFF_MU', MU, INFO%ISPIN, IERR, WRITEXMLINCAR, FOUNDNUMBER=IDUM)
   IF (IDUM < INFO%ISPIN) MU(2) = MU(1)
   ! Width of the cut-off function
   SIGMA = HUGE(0._q)
   CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'CUTOFF_SIGMA', SIGMA, INFO%ISPIN, IERR, WRITEXMLINCAR, FOUNDNUMBER=IDUM)
   IF (IDUM < INFO%ISPIN) SIGMA(2) = SIGMA(1)
   ! Type of the cut-off function
   CUTOFF_TYPE_STR = '--'
   CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'CUTOFF_TYPE', CUTOFF_TYPE_STR, 20, IERR, WRITEXMLINCAR, FOUNDNUMBER=IDUM)
   CUTOFF_TYPE = CUTOFF_NONE
   IF (CUTOFF_TYPE_STR /= '--') THEN
      IF (NOCCUR(LOWERCASE(CUTOFF_TYPE_STR),'erfc',-1))     CUTOFF_TYPE = CUTOFF_ERFC
      IF (NOCCUR(LOWERCASE(CUTOFF_TYPE_STR),'gaussian',-1)) CUTOFF_TYPE = CUTOFF_GAUSSIAN
      IF (NOCCUR(LOWERCASE(CUTOFF_TYPE_STR),'fermi',-1))    CUTOFF_TYPE = CUTOFF_FERMI
      IF (NOCCUR(LOWERCASE(CUTOFF_TYPE_STR),'num_wann',-1)) CUTOFF_TYPE = CUTOFF_NUM_WANN
   ENDIF
   CALL CUTOFF_DES_INIT(CUTOFF_DES,CUTOFF_TYPE,INFO%ISPIN,MU,SIGMA,NUM_WANN)

   ! Switch on the PEAD routines
   IF (MLWFRUN()) THEN
#ifdef gammareal
      ! for the gamma-only version we only need the pead routines when MMN is requested
      IF (WRITE_MMN_AMN.OR.WANNIER90RUN()) THEN
         CALL PEAD_REQUEST
      ENDIF
#else
      CALL PEAD_REQUEST
#endif
      CALL USE_FULL_KPOINTS
   ENDIF

   CALL CLOSE_INCAR_IF_FOUND(IU5)

 END SUBROUTINE MLWF_READER

!***********************************************************************
!
! WANNIER90
!
!***********************************************************************
 FUNCTION MLWFRUN()
   use locproj, only: LPRJ_functions
   LOGICAL MLWFRUN
   MLWFRUN = ALLOCATED(LPRJ_functions).OR.WANNIER90().OR.LSCDM
 END FUNCTION

 FUNCTION WANNIER90()
   LOGICAL WANNIER90
   WANNIER90=LWANNIER90.OR.LWANNIER90_RUN
 END FUNCTION WANNIER90

 FUNCTION WANNIER90RUN()
   LOGICAL WANNIER90RUN
   WANNIER90RUN=LWANNIER90_RUN
 END FUNCTION WANNIER90RUN

!***********************************************************************
!> @brief Main interface point between VASP and Wannier90
!***********************************************************************
 SUBROUTINE MLWF_MAIN(WDES,W,P,KPOINTS,CQIJ,T_INFO,LATT_CUR,INFO,IO,MLWF,WANNIER_RUN)
   USE ini
   USE base
   USE constant
   USE pseudo
   USE mkpoints
   USE poscar
   USE lattice
   USE wave_high
   USE radial
   USE full_kpoints
   USE locproj, ONLY : WANPROJFILE, LPRJ_FROM_PSEUDO, LOCPROJ_AUTO, &
                       LPRJ_functions, LPRJ_WRITE, LPRJ_WRITE_FUNCTIONS

   TYPE(wavedes) WDES
   TYPE(wavespin) W
   TYPE(potcar) P(:)
   TYPE(kpoints_struct) KPOINTS
   TYPE(type_info) T_INFO
   TYPE(latt) LATT_CUR
   TYPE(info_struct) INFO
   TYPE(in_struct) IO
   TYPE(MLWF_type),OPTIONAL,TARGET :: MLWF
   LOGICAL,OPTIONAL,INTENT(IN) :: WANNIER_RUN
   OVERLAP CQIJ(:,:,:,:)

   !local variables
   TYPE(MLWF_type),TARGET :: MLWF_TMP
   TYPE(MLWF_type),POINTER :: P_MLWF
   LOGICAL,PARAMETER :: WRITE_MMN = .TRUE.
   LOGICAL,PARAMETER :: WRITE_AMN = .TRUE.
   LOGICAL,PARAMETER :: WRITE_EIG = .TRUE.
   COMPLEX(q),ALLOCATABLE :: A_MATRIX_TMP(:,:,:,:)
   LOGICAL :: MY_LWANNIER90_RUN
   INTEGER :: NUM_WANN_TMP
   INTEGER :: IB,IBP,IW,IK,IS
   INTEGER :: NCDIJ, NKPTS, NB_TOT

   PROFILING_START('mlwf_main')

   IF (IO%IU0>=0) THEN
        WRITE(IO%IU6,*) 'MLWF driver (Wannier functions)'
        WRITE(IO%IU6,*) '==============================='
   ENDIF

   MY_LWANNIER90_RUN = LWANNIER90_RUN; IF (PRESENT(WANNIER_RUN)) MY_LWANNIER90_RUN = WANNIER_RUN
   P_MLWF => MLWF_TMP; IF (PRESENT(MLWF)) P_MLWF => MLWF

   ! We can get initial projections in different ways defined by the MLWF%PROJ_MODE variable
   IF (LSCDM)                     P_MLWF%PROJ_MODE = PROJ_SCDM
   IF (ALLOCATED(LPRJ_functions)) P_MLWF%PROJ_MODE = PROJ_LOCPROJ
   IF (LOCPROJ_AUTO)              P_MLWF%PROJ_MODE = PROJ_LOCPROJAUTO
   ! WANNIER90_MODE has to be determined inside MLWF_WANNIER90_SETUP
   SELECT CASE (P_MLWF%PROJ_MODE)
     ! LOCPROJ parser (read from global parser)
     CASE (PROJ_LOCPROJ)
         P_MLWF%NUM_WANN = SIZE(LPRJ_functions)
         ALLOCATE(P_MLWF%LPRJ_functions(size(LPRJ_functions)))
         P_MLWF%LPRJ_functions = LPRJ_functions
     ! Automatic projections from the pseudopotential
     CASE (PROJ_LOCPROJAUTO)
         P_MLWF%NUM_WANN = SIZE(LPRJ_functions)
         ALLOCATE(P_MLWF%LPRJ_functions(size(LPRJ_functions)))
         P_MLWF%LPRJ_functions = LPRJ_functions
     ! SCDM projections (default)
     CASE (PROJ_SCDM,PROJ_UNKNOWN)
         IF (NUM_WANN == 0) THEN
            ! LPRJ_functions is not modified here because NUM_WANN_TMP is passed
            CALL LPRJ_FROM_PSEUDO(T_INFO,P,WDES%LNONCOLLINEAR,P_MLWF%LPRJ_functions,NUM_WANN_TMP)
            P_MLWF%NUM_WANN = WDES%NB_TOT; IF (NUM_WANN_TMP < WDES%NB_TOT) P_MLWF%NUM_WANN = NUM_WANN_TMP
         ELSE
            P_MLWF%NUM_WANN = NUM_WANN
         ENDIF
   END SELECT
   ! Call Wannier90 setup
   CALL MLWF_WANNIER90_SETUP(P_MLWF,WDES,W,KPOINTS,P,CQIJ,T_INFO,LATT_CUR,INFO,IO)
   ! At this point the projection mode must be known
   IF (IO%IU0>=0) WRITE(IO%IU0,*) PROJ_STRING(P_MLWF%PROJ_MODE),' mode'
   IF (IO%IU6>=0) WRITE(IO%IU6,*) PROJ_STRING(P_MLWF%PROJ_MODE),' mode'
   IF (P_MLWF%PROJ_MODE==TRAFO_UNKNOWN) THEN
         CALL vtutor%alert("No projection mode was specified."//&
              " You can specify projections using wannier90 (NUM_WANN and WANNIER90_WIN tags),"//&
              " LOCPROJ (LOCPROJ tag) or the SCDM method (NUM_WANN and LSCDM tags).")
   END IF
   ! Compute MMN
   IF (WRITE_MMN_AMN.OR.MY_LWANNIER90_RUN) THEN
      CALL MLWF_WANNIER90_MMN(P_MLWF,WDES,W,P,CQIJ,T_INFO,LATT_CUR,WRITE_MMN_AMN,INFO,IO)
      IF (WRITE_MMN_AMN) CALL MLWF_WANNIER90_WRITE_MMN(P_MLWF,IO)
   ENDIF
   ! Write .eig file
   CALL MLWF_WANNIER90_EIG(P_MLWF,W,WRITE_MMN_AMN,IO)
   ! Estimate cutoff descriptor
   IF (P_MLWF%PROJ_MODE==PROJ_SCDM.OR.&
       P_MLWF%PROJ_MODE==PROJ_LOCPROJ.OR.&
       P_MLWF%PROJ_MODE==PROJ_LOCPROJAUTO.OR.&
       LWANNIER90_AUTO_WINDOW) THEN
      CALL CUTOFF_DES_ESTIMATE(CUTOFF_DES,WDES,W,P_MLWF%NUM_WANN,KPOINTS,INFO,IO)
      CALL CUTOFF_DES_PRINT(CUTOFF_DES,IO)
   ENDIF
   IF (LWANNIER90_AUTO_WINDOW) CALL ESTIMATE_AUTO_WINDOW(MLWF,W)
   ! Compute AMN
   IF (P_MLWF%PROJ_MODE==PROJ_SCDM) THEN
      CALL CALC_WANNIER_TRAFO_SCDM(W,P_MLWF%NUM_WANN,A_MATRIX_TMP,P,LATT_CUR,P_MLWF%NKPTS,P_MLWF%kpt_latt,IO)
      ALLOCATE(MLWF%A_matrix(MLWF%NB_TOT-MLWF%NEXCLB,MLWF%NUM_WANN,MLWF%NKPTS,MLWF%ISPIN))
      MLWF%A_matrix = CMPLX(0.0_q,0.0_q,q)
      ! Need to transfer the data to what wannier90 expects
      DO IW = 1, MLWF%NUM_WANN
         IBP = 0
         DO IB = 1, MLWF%NB_TOT
            IF (MLWF%EXCLUDE_BAND(IB)) CYCLE
            IBP = IBP + 1
            MLWF%A_matrix(IBP,IW,:,:) = A_MATRIX_TMP(IB,IW,:,:)
         ENDDO
      ENDDO
      DEALLOCATE(A_MATRIX_TMP)
   ELSE
      CALL LPRJ_WRITE_FUNCTIONS(P_MLWF%LPRJ_functions,MLWF%NRSPINORS==2,IO)
      CALL MLWF_CALC_AMN(P_MLWF,WDES,W,P,CQIJ,T_INFO,LATT_CUR,INFO,IO)
   ENDIF
   IF (WRITE_MMN_AMN.AND.P_MLWF%PROJ_MODE/=PROJ_UNKNOWN) CALL MLWF_WANNIER90_WRITE_AMN(MLWF,IO)
   ! At this point LPRJ_COVL conntains the projections as computed from LOCPROJ
   ! we can write them to disk
   IF (MLWF%PROJ_MODE==PROJ_LOCPROJAUTO.OR.MLWF%PROJ_MODE==PROJ_LOCPROJ) THEN
      CALL LPRJ_WRITE(MLWF%LPRJ_functions,MLWF%LPRJ_COVL,IO%IU6,IO%IU0,W)
   ENDIF
   ! Write .win file if we are not in legacy mode
   IF (IO%IU6>=0.AND..NOT.WANNIER90_WIN_LEGACY_MODE) CALL MLWF_WANNIER90_WIN(P_MLWF)
   ! If we do not call Wannier90 from VASP and we do not have any projections then we can exit early
   IF ((P_MLWF%PROJ_MODE==PROJ_UNKNOWN).AND.(.NOT.MY_LWANNIER90_RUN)) RETURN
   ! We might need UNK files for wannier90
   IF (WRITE_UNK) THEN
      CALL MLWF_WRITE_UNK_XSF(P_MLWF,W,P,LATT_CUR,WRITE_UNK,.FALSE.,IO)
   ENDIF
   ! run Wannier90
   IF (MY_LWANNIER90_RUN) THEN
      CALL MLWF_WANNIER90_RUN(P_MLWF,IO)
   ! or the one-shot method
   ! or unitary matrix contrained optimization (umco)
   ELSE
      CALL MLWF_TRAFO_RUN(P_MLWF,W,P,CQIJ,LATT_CUR,T_INFO,INFO,IO)
   ENDIF
   ! Final output
   IF (LWRITE_WANPROJ) CALL MLWF_WRITE_WANPROJ(P_MLWF,WANPROJFILE,IO)
   IF (LWRITE_WANNIER_XSF) THEN
      CALL MLWF_WRITE_UNK_XSF(P_MLWF,W,P,LATT_CUR,.FALSE.,LWRITE_WANNIER_XSF,IO)
   ENDIF
   ! Write wavecar file with wannier orbitals
   IF (LWAVE_WANNIER) CALL MLWF_WRITE_WAVECAR(P_MLWF,W,LATT_CUR,IO)
   ! write XSF with the BLOCH wavefuncitons
   IF (LWRITE_BLOCH_XSF) THEN
      IF (IO%IU0>=0) WRITE(IO%IU0,*) "Write Bloch XSF files"
      CALL WRITE_WAVESPIN_XSF('bloch',W,LATT_CUR,P,W%WDES%NB_TOT,&
      MLWF%atom_symbols,MLWF%atoms_cart,WRITE_BLOCH_XSF_IW)
   ENDIF
   CALL MLWF_WRITE_FINAL(MLWF,IO)

   PROFILING_STOP('mlwf_main')
 END SUBROUTINE MLWF_MAIN

 !> Free arrays in mlwf
 SUBROUTINE MLWF_FREE(MLWF)
   TYPE(MLWF_type) :: MLWF
   IF (ALLOCATED(MLWF%EXCLUDE_BAND)) DEALLOCATE(MLWF%EXCLUDE_BAND)
   IF (ALLOCATED(MLWF%LPRJ_functions)) DEALLOCATE(MLWF%LPRJ_functions)
   IF (ALLOCATED(MLWF%NEQUIV)) DEALLOCATE(MLWF%NEQUIV)
   IF (ALLOCATED(MLWF%kpt_latt)) DEALLOCATE(MLWF%kpt_latt)
   IF (ALLOCATED(MLWF%atoms_cart)) DEALLOCATE(MLWF%atoms_cart)
   IF (ALLOCATED(MLWF%atom_symbols)) DEALLOCATE(MLWF%atom_symbols)
   IF (ALLOCATED(MLWF%NNLIST)) DEALLOCATE(MLWF%NNLIST)
   IF (ALLOCATED(MLWF%NNCELL)) DEALLOCATE(MLWF%NNCELL)

   IF (ALLOCATED(MLWF%M_matrix)) DEALLOCATE(MLWF%M_matrix)
   IF (ALLOCATED(MLWF%A_matrix)) DEALLOCATE(MLWF%A_matrix)
   IF (ALLOCATED(MLWF%eigenvalues)) DEALLOCATE(MLWF%eigenvalues)

   IF (ALLOCATED(MLWF%U_matrix)) DEALLOCATE(MLWF%U_matrix)
   IF (ALLOCATED(MLWF%U_matrix_opt)) DEALLOCATE(MLWF%U_matrix_opt)

   IF (ALLOCATED(MLWF%LPRJ_COVL)) DEALLOCATE(MLWF%LPRJ_COVL)
   IF (ALLOCATED(MLWF%LPRJ_COVL_IAO)) DEALLOCATE(MLWF%LPRJ_COVL_IAO)

   IF (ALLOCATED(MLWF%lwindow)) DEALLOCATE(MLWF%lwindow)
   IF (ALLOCATED(MLWF%lexclude_band)) DEALLOCATE(MLWF%lexclude_band)
   IF (ALLOCATED(MLWF%wann_centres)) DEALLOCATE(MLWF%wann_centres)
   IF (ALLOCATED(MLWF%wann_spreads)) DEALLOCATE(MLWF%wann_spreads)
   IF (ALLOCATED(MLWF%spread)) DEALLOCATE(MLWF%spread)
 END SUBROUTINE MLWF_FREE

 !> Convert projection mode enumerator to a human readable string
 FUNCTION PROJ_STRING(PROJ_MODE) RESULT(STRING)
   INTEGER,INTENT(IN) :: PROJ_MODE
   CHARACTER(:),ALLOCATABLE :: STRING
   SELECT CASE (PROJ_MODE)
       CASE(PROJ_SCDM)
           STRING = 'SCDM'
       CASE(PROJ_LOCPROJ)
           STRING = 'LOCPROJ'
       CASE(PROJ_LOCPROJAUTO)
           STRING = 'LOCPROJ_AUTO'
       CASE(PROJ_WANNIER90)
           STRING = 'Wannier90'
       CASE(TRAFO_UNKNOWN)
           STRING = 'UNKOWN'
   END SELECT
 END FUNCTION PROJ_STRING

 !> Convert trafo mode enumerator to a human readable string
 FUNCTION TRAFO_STRING(TRAFO_MODE) RESULT(STRING)
   INTEGER,INTENT(IN) :: TRAFO_MODE
   CHARACTER(:),ALLOCATABLE :: STRING
   SELECT CASE (TRAFO_MODE)
       CASE(TRAFO_ONESHOTSVD)
           STRING = 'One-shot SVD method'
       CASE(TRAFO_WANNIER90_MODE)
           STRING = 'Wannier90 Libary'
       CASE(TRAFO_UMCO)
           STRING = 'Unitary matrix contrained optimization (UMCO)'
       CASE(TRAFO_SCDM)
           STRING = 'Selected columns of the density matrix (SCDM)'
       CASE(TRAFO_UNKNOWN)
           STRING = 'UNKOWN'
   END SELECT
 END FUNCTION TRAFO_STRING

 !> Write information about the method that was used to tranfer from AMN to UMN
 SUBROUTINE MLWF_WRITE_FINAL(MLWF,IO)
   use base, ONLY: in_struct
   TYPE(MLWF_type),INTENT(IN) :: MLWF
   TYPE(in_struct),INTENT(IN) :: IO

   IF (IO%IU0>=0) THEN
        WRITE(IO%IU6,'(A)') ''
        WRITE(IO%IU6,'(A)') ' MMN -> overlap matrix elements'
        WRITE(IO%IU6,'(A)') ' AMM -> projections of Bloch states onto localized states'
        WRITE(IO%IU6,'(A)') ' UMN -> unitary matrix transforming from Bloch states to Wannier states'
        WRITE(IO%IU6,'(A)') ''

        WRITE(IO%IU6,*) 'Parameters of MLWF datastructure'
        WRITE(IO%IU6,*) 'num_band = ',MLWF%NB_TOT
        WRITE(IO%IU6,*) 'num_wann = ',MLWF%NUM_WANN
        WRITE(IO%IU6,'(A)') ''
        WRITE(IO%IU6,'(A)') ' Computiation of AMN'
        WRITE(IO%IU6,*) 'proj = ',PROJ_STRING(MLWF%PROJ_MODE)
        WRITE(IO%IU6,'(A)') ''

        WRITE(IO%IU6,'(A)') ' Transformation from AMN to UMN'
        WRITE(IO%IU6,*) 'trafo = ',TRAFO_STRING(MLWF%TRAFO_MODE)
        WRITE(IO%IU6,'(A)') ''
        WRITE(IO%IU6,'(A)') '------------------------ end of MLWF driver reached ----------------------------------------------------'
   ENDIF
 END SUBROUTINE MLWF_WRITE_FINAL

 !> Initialize the cutoff descriptor with all the necessary information
 SUBROUTINE CUTOFF_DES_INIT(SELF,CUTOFF_ITYPE,ISPIN,MU,SIGMA,NUM_WANN)
   TYPE(cutoff_type) :: SELF
   INTEGER,INTENT(IN) :: ISPIN
   INTEGER,INTENT(IN) :: CUTOFF_ITYPE
   REAL(q),INTENT(IN) :: MU(2)
   REAL(q),INTENT(IN) :: SIGMA(2)
   INTEGER,INTENT(IN) :: NUM_WANN
   SELF%TYPE = CUTOFF_ITYPE
   SELF%ISPIN = ISPIN
   SELF%MU = MU
   SELF%SIGMA = SIGMA
   SELF%NUM_WANN = NUM_WANN
 END SUBROUTINE CUTOFF_DES_INIT

 !> @brief Estimate cutoff descriptor which contains the necessary information
 !> to use the SCDM method or to obtain unitary transformations from LOCPROJ
 SUBROUTINE CUTOFF_DES_ESTIMATE(CUTOFF_DES,WDES,W,NUM_WANN,KPOINTS,INFO,IO)

   USE base
   USE wave_high
   USE mkpoints

   TYPE(cutoff_type) :: CUTOFF_DES
   TYPE(wavedes) WDES
   TYPE(wavespin) W
   TYPE(kpoints_struct) KPOINTS
   TYPE(info_struct) INFO
   TYPE(in_struct) IO

   ! local variables
   INTEGER,PARAMETER :: MAX_ITER = 1000
   LOGICAL :: IS_METAL
   INTEGER :: N, IK, ISP, NUM_WANN
   REAL(q) :: SCDM_SIGMA_MIN, SCDM_SIGMA_MAX
   REAL(q) :: NELECT, NELECT_MIN, NELECT_MAX

   ! Needed for DENSTA.
   INTEGER, PARAMETER :: NEDOS=256
   REAL(q) :: ENTROPY, EFERMI, SIGMA, OCC
   REAL(q) :: DOS(NEDOS,W%WDES%ISPIN),DOSI(NEDOS,W%WDES%ISPIN)
   REAL(q) :: PAR(1,1,1,1,W%WDES%NCDIJ),DOSPAR(1,1,1,W%WDES%NCDIJ)
   TYPE (wavespin) W_tmp

   ! Simple heuristic to determine CUTOFF_MU and SCDM_SIGMA from NUM_WANN
   IF (CUTOFF_DES%TYPE==CUTOFF_NONE) THEN
      ! HM: Perhaps the erfc makes more sense in general even in the case of a metal
      CUTOFF_DES%TYPE = CUTOFF_ERFC!; IF (isMetal(W%FERTOT)) CUTOFF_DES%TYPE = CUTOFF_GAUSSIAN
      IF (WDES%NB_TOT<=NUM_WANN) THEN
          CUTOFF_DES%TYPE = CUTOFF_NONE
          CUTOFF_DES%NUM_WANN = WDES%NB_TOT
      ENDIF
  ENDIF

   IF (CUTOFF_DES%TYPE==CUTOFF_GAUSSIAN) THEN
      DO ISP=1,WDES%ISPIN
         IF (CUTOFF_DES%MU(ISP)==HUGE(0._q)) THEN
            CUTOFF_DES%MU(ISP) = W%EFERMI(ISP)
         ENDIF
         IF (CUTOFF_DES%SIGMA(ISP)==HUGE(0._q)) THEN
            CUTOFF_DES%SIGMA(ISP) = KPOINTS%SIGMA
            ! Determine SIGMA using bissection method such that the integral yields NUM_WANN
            SCDM_SIGMA_MAX = (KPOINTS%EMAX-KPOINTS%EMIN)*10
            SCDM_SIGMA_MIN = 0.1
            NELECT_MAX = 0
            NELECT_MIN = 0
            DO IK=1,W%WDES%NKPTS
               CUTOFF_DES%SIGMA(ISP) = SCDM_SIGMA_MAX
               NELECT_MAX = NELECT_MAX + SUM(CUTOFF_DES_GET_CUTOFF(CUTOFF_DES,ISP,W%CELEN(:,IK,ISP)))*W%WDES%WTKPT(IK)
               CUTOFF_DES%SIGMA(ISP) = SCDM_SIGMA_MIN
               NELECT_MIN = NELECT_MIN + SUM(CUTOFF_DES_GET_CUTOFF(CUTOFF_DES,ISP,W%CELEN(:,IK,ISP)))*W%WDES%WTKPT(IK)
            ENDDO

            ! Bissection method
            DO N=1,MAX_ITER
               CUTOFF_DES%SIGMA(ISP) = (SCDM_SIGMA_MAX+SCDM_SIGMA_MIN)/2
               NELECT = 0
               DO IK=1,W%WDES%NKPTS
                  NELECT = NELECT + SUM(CUTOFF_DES_GET_CUTOFF(CUTOFF_DES,ISP,W%CELEN(:,IK,ISP)))*W%WDES%WTKPT(IK)
               ENDDO
               IF (SIGN(1.0_q,NELECT-NUM_WANN*0.5) == SIGN(1.0_q,NELECT_MIN-NUM_WANN*0.5)) THEN
                  SCDM_SIGMA_MIN = CUTOFF_DES%SIGMA(ISP)
                  NELECT_MIN = NELECT
               ELSE
                  SCDM_SIGMA_MAX = CUTOFF_DES%SIGMA(ISP)
                  NELECT_MAX = NELECT
               ENDIF
               IF (ABS(NELECT-NUM_WANN*0.5)<1e-6) EXIT
            END DO
         ENDIF
      ENDDO
   ELSE IF (CUTOFF_DES%TYPE==CUTOFF_ERFC) THEN
      DO ISP=1,WDES%ISPIN
         IF (CUTOFF_DES%SIGMA(ISP)==HUGE(0._q)) THEN
            CUTOFF_DES%SIGMA(ISP) = 1
         ENDIF
         IF (CUTOFF_DES%MU(ISP)==HUGE(0._q)) THEN
            ! Some trickery is required.
            W_TMP=W
            ALLOCATE(W_TMP%FERTOT(W%WDES%NB_TOT, W%WDES%NKPTS, W%WDES%ISPIN))
            W_TMP%FERWE => W_TMP%FERTOT(W%WDES%NB_LOW : W%WDES%NB_TOT : W%WDES%NB_PAR, :, :)

            ! Calculate the partial occupancies using smearing.
            OCC = 2._q; IF (W%WDES%ISPIN==2.OR.W%WDES%NRSPINORS==2) OCC = 1._q
            CALL DENSTA(-1, -1, W%WDES, W_TMP, KPOINTS, NUM_WANN * 0.8 * OCC, &
               0._q, ENTROPY, EFERMI, CUTOFF_DES%SIGMA(ISP), .FALSE., &
               NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
            CUTOFF_DES%MU(ISP) = EFERMI
         ENDIF
      ENDDO
   ELSE IF (CUTOFF_DES%TYPE==CUTOFF_NUM_WANN) THEN
      DO ISP=1,WDES%ISPIN
         IF (CUTOFF_DES%SIGMA(ISP)==HUGE(0._q)) THEN
            CUTOFF_DES%SIGMA(ISP) = 1
         ENDIF
      ENDDO
   ENDIF

 END SUBROUTINE CUTOFF_DES_ESTIMATE

 !> @brief Print SCDM descriptor
 SUBROUTINE CUTOFF_DES_PRINT(CUTOFF_DES,IO)
   USE base

   TYPE(cutoff_type) :: CUTOFF_DES
   TYPE(in_struct) :: IO

   IF (IO%IU6>=0) THEN
      IF (CUTOFF_DES%TYPE/=CUTOFF_NONE) THEN
         WRITE(IO%IU6,'(2A)') ''
         WRITE(IO%IU6,'(2A)') ' Information about cutoff-function'
         WRITE(IO%IU6,'(2A)') ' ---------------------------------'
      ENDIF
      IF (CUTOFF_DES%TYPE==CUTOFF_ERFC)     WRITE(IO%IU6,'(2A)') ' CUTOFF_TYPE  = ','erfc'
      IF (CUTOFF_DES%TYPE==CUTOFF_GAUSSIAN) WRITE(IO%IU6,'(2A)') ' CUTOFF_TYPE  = ','gaussian'
      IF (CUTOFF_DES%TYPE==CUTOFF_FERMI)    WRITE(IO%IU6,'(2A)') ' CUTOFF_TYPE  = ','fermi'
      IF (CUTOFF_DES%TYPE==CUTOFF_NUM_WANN) WRITE(IO%IU6,'(2A)') ' CUTOFF_TYPE  = ','num_wann'
      IF (CUTOFF_DES%TYPE/=CUTOFF_NONE) THEN
         IF (CUTOFF_DES%ISPIN==1) THEN
            WRITE(IO%IU6,'(A,F6.2,A)') ' CUTOFF_MU    = ',CUTOFF_DES%MU(1),   ' eV'
            WRITE(IO%IU6,'(A,F6.2,A)') ' CUTOFF_SIGMA = ',CUTOFF_DES%SIGMA(1),' eV'
         ELSE
            WRITE(IO%IU6,'(A)')         '                  spin_up  spin_dw'
            WRITE(IO%IU6,'(A,2F6.2,A)') ' CUTOFF_MU    = ',CUTOFF_DES%MU,   '  eV'
            WRITE(IO%IU6,'(A,2F6.2,A)') ' CUTOFF_SIGMA = ',CUTOFF_DES%SIGMA,'  eV'
         ENDIF
         WRITE(IO%IU6,'(2A)') ''
      ENDIF
   ENDIF

 END SUBROUTINE CUTOFF_DES_PRINT

 !> Estimate Auto window for Wannier90 localization
 SUBROUTINE ESTIMATE_AUTO_WINDOW(MLWF,W)
   USE wave_high, ONLY: wavespin

   TYPE(mlwf_type),INTENT(INOUT) :: MLWF
   TYPE(wavespin),INTENT(IN) :: W

   INTEGER :: IK,IS,IB
   REAL(q) :: ELIM(2,MLWF%NB_TOT)

   ! set froz_max
   MLWF%DIS_FROZ_MAX = CUTOFF_DES%MU(1)

   ELIM(1,:) =  HUGE(1_q)
   ELIM(2,:) = -HUGE(1_q)
   ! if more frozen states than wannier states decrease froz_max until ok
   DO IS=1,MLWF%ISPIN
      DO IK=1,W%WDES%NKPTS
         DO IB=1,W%WDES%NB_TOT
            ELIM(1,IB) = MIN(ELIM(1,IB),REAL(W%CELTOT(IB,IK,IS)))
            ELIM(2,IB) = MAX(ELIM(2,IB),REAL(W%CELTOT(IB,IK,IS)))
         END DO
      END DO
   END DO
   MLWF%DIS_FROZ_MAX = ELIM(1,MLWF%NUM_WANN)-1e-5

 END SUBROUTINE ESTIMATE_AUTO_WINDOW

!> @brief Call Wannier90 to read projections
 SUBROUTINE MLWF_WANNIER90_SETUP(MLWF,WDES,W,KPOINTS,P,CQIJ,T_INFO,LATT_CUR,INFO,IO)
   USE base
   USE constant
   USE pseudo
   USE poscar
   USE lattice
   USE wave_high
   USE radial
   USE mkpoints_struct_def, ONLY : kpoints_struct
   USE full_kpoints, ONLY : KPOINTS_FULL, LIDENTICAL_KPOINT

   TYPE(mlwf_type),INTENT(INOUT) :: MLWF
   TYPE(wavedes),INTENT(IN) :: WDES
   TYPE(wavespin),INTENT(IN) :: W
   TYPE(kpoints_struct),INTENT(IN) :: KPOINTS
   TYPE(potcar),INTENT(IN) :: P(:)
   OVERLAP,INTENT(IN) :: CQIJ(:,:,:,:)
   TYPE(type_info),INTENT(IN) :: T_INFO
   TYPE(latt),INTENT(IN) :: LATT_CUR
   TYPE(info_struct),INTENT(IN) :: INFO
   TYPE(in_struct),INTENT(IN) :: IO

   ! local variables
   INTEGER, PARAMETER :: num_nnmax=12
   LOGICAL :: spinors
   LOGICAL :: LWIN_FOUND
   LOGICAL :: ZONA_WARNING
   INTEGER :: N, NI, NKP, M, IS, IK, NEXCLB
   INTEGER :: READ_NUM_WANN !< Number of wannier orbitals returned by call to wannier interface
   INTEGER :: ESTIMATED_NUM_WANN !< Estimated number of orbitals based on the radial types set to be different from zero
   INTEGER :: I1,I2,I3
   REAL(q) :: V(3),V1(3),SHIFT(3)
   CHARACTER(LEN=WANNIER90_WIN_MAXLEN) :: WANNIER90_WIN_FILTERED

   ! wannier90_setup variables: output
   INTEGER :: num_bands
   REAL(q), ALLOCATABLE :: proj_site(:,:)
   INTEGER, ALLOCATABLE :: proj_l(:)
   INTEGER, ALLOCATABLE :: proj_m(:)
   INTEGER, ALLOCATABLE :: proj_s(:)
   INTEGER, ALLOCATABLE :: proj_radial(:)
   REAL(q), ALLOCATABLE :: proj_z(:,:)
   REAL(q), ALLOCATABLE :: proj_x(:,:)
   REAL(q), ALLOCATABLE :: proj_zona(:)
   REAL(q), ALLOCATABLE :: proj_s_qaxisx(:,:)
   INTEGER, ALLOCATABLE :: exclude_bands(:)

   ! setup basic dimensions in MLWF
   MLWF%ISPIN = WDES%ISPIN
   MLWF%NRSPINORS = WDES%NRSPINORS
   MLWF%NB_TOT = WDES%NB_TOT
   MLWF%SZNAM1 = INFO%SZNAM1
#define use_kptsfull
#ifdef use_kptsfull
   MLWF%NKPTS = KPOINTS_FULL%NKPTS
#else
   MLWF%NKPTS = KPOINTS%NKPX*KPOINTS%NKPY*KPOINTS%NKPZ
#endif
   MLWF%COMM = WDES%COMM
   IF (MLWF%NUM_WANN == 0) THEN
      CALL vtutor%error('NUM_WANN must be different from zero when using projections from Wannier90.')
   ENDIF

   MLWF%mp_grid(1)=KPOINTS%NKPX
   MLWF%mp_grid(2)=KPOINTS%NKPY
   MLWF%mp_grid(3)=KPOINTS%NKPZ

   ALLOCATE(MLWF%kpt_latt(3,MLWF%NKPTS))
   ALLOCATE(MLWF%NEQUIV(MLWF%NKPTS))

#ifdef use_kptsfull
   MLWF%kpt_latt=KPOINTS_FULL%VKPT
   MLWF%NEQUIV  =KPOINTS_FULL%NEQUIV
#else
   N = 0
   SHIFT = [0.0_q,0.0_q,0.0_q]
   DO IK=1,KPOINTS_FULL%NKPTS
      outer: DO I3=0,KPOINTS%NKPZ-1
         DO I2=0,KPOINTS%NKPY-1
            DO I1=0,KPOINTS%NKPX-1
               ! k-point coordinates in the basis BK:
               V(1)=(REAL(I1,KIND=q)+SHIFT(1))/FLOAT(KPOINTS%NKPX)
               V(2)=(REAL(I2,KIND=q)+SHIFT(2))/FLOAT(KPOINTS%NKPY)
               V(3)=(REAL(I3,KIND=q)+SHIFT(3))/FLOAT(KPOINTS%NKPZ)
               !V1=MOD(V+6.5_q,1._q)-0.5_q
               IF (LIDENTICAL_KPOINT(KPOINTS_FULL%VKPT(:,IK),V)) THEN
                  N=N+1
                  MLWF%NEQUIV(N)     = KPOINTS_FULL%NEQUIV(IK)
                  MLWF%kpt_latt(:,N) = KPOINTS_FULL%VKPT(:,IK)
                  EXIT outer
               ENDIF
            END DO
         END DO
      END DO outer
   END DO
   IF (N/=MLWF%NKPTS) CALL vtutor%error('Failed to generate regular grid of k-points from generating basis.')
#endif

   MLWF%real_lattice(1,:)=LATT_CUR%A(:,1)
   MLWF%real_lattice(2,:)=LATT_CUR%A(:,2)
   MLWF%real_lattice(3,:)=LATT_CUR%A(:,3)

   MLWF%recip_lattice(1,:)=LATT_CUR%B(:,1)
   MLWF%recip_lattice(2,:)=LATT_CUR%B(:,2)
   MLWF%recip_lattice(3,:)=LATT_CUR%B(:,3)

   MLWF%recip_lattice=MLWF%recip_lattice*TPI

   MLWF%num_atoms=T_INFO%NIONS
   MLWF%DIS_FROZ_MIN = HUGE(MLWF%DIS_FROZ_MIN)
   MLWF%DIS_FROZ_MAX = HUGE(MLWF%DIS_FROZ_MAX)
   MLWF%DIS_WIN_MIN = HUGE(MLWF%DIS_WIN_MIN)
   MLWF%DIS_WIN_MAX = HUGE(MLWF%DIS_WIN_MAX)

   ALLOCATE(MLWF%atom_symbols(MLWF%num_atoms),MLWF%atoms_cart(3,MLWF%num_atoms))
   DO NI=1,T_INFO%NIONS
      MLWF%atom_symbols(NI)=T_INFO%TYPE(T_INFO%ITYP(NI))
!     WRITE(MLWF%atom_symbols(NI),'(A2,18X)') T_INFO%TYPE(T_INFO%ITYP(NI))

      MLWF%atoms_cart(1,NI)=LATT_CUR%A(1,1)*T_INFO%POSION(1,NI)+ &
                            LATT_CUR%A(1,2)*T_INFO%POSION(2,NI)+ &
                            LATT_CUR%A(1,3)*T_INFO%POSION(3,NI)
      MLWF%atoms_cart(2,NI)=LATT_CUR%A(2,1)*T_INFO%POSION(1,NI)+ &
                            LATT_CUR%A(2,2)*T_INFO%POSION(2,NI)+ &
                            LATT_CUR%A(2,3)*T_INFO%POSION(3,NI)
      MLWF%atoms_cart(3,NI)=LATT_CUR%A(3,1)*T_INFO%POSION(1,NI)+ &
                            LATT_CUR%A(3,2)*T_INFO%POSION(2,NI)+ &
                            LATT_CUR%A(3,3)*T_INFO%POSION(3,NI)
   ENDDO

   spinors=WDES%LNONCOLLINEAR

   ! A minimal wannier90.win file must exist; it must at least
   ! contain the keyword "num_wann"
   IF (IO%IU6>=0.AND..NOT.WANNIER90_WIN_LEGACY_MODE) CALL MLWF_WANNIER90_WIN(MLWF,LSETUP=.TRUE.)

   ! Initialize everything to zero
   MLWF%NNTOT=0; num_bands=0; num_wann=0

   ALLOCATE(MLWF%NNLIST(MLWF%NKPTS,NUM_NNMAX),MLWF%NNCELL(3,MLWF%NKPTS,NUM_NNMAX))
   MLWF%NNLIST=0; MLWF%NNCELL=0

   ALLOCATE(proj_site(3,MLWF%NB_TOT),proj_l(MLWF%NB_TOT),proj_m(MLWF%NB_TOT), &
  &   proj_radial(MLWF%NB_TOT),proj_z(3,MLWF%NB_TOT),proj_x(3,MLWF%NB_TOT), &
  &   proj_zona(MLWF%NB_TOT),exclude_bands(MLWF%NB_TOT),proj_s(MLWF%NB_TOT),&
  &   proj_s_qaxisx(3,MLWF%NB_TOT))
   proj_site=0; proj_l=0; proj_m=0; proj_radial=0; proj_z=0; proj_x=0; proj_zona=0; exclude_bands=0
   proj_s=0;proj_s_qaxisx=0
   READ_NUM_WANN = 0

   ! Only one node will do the actual work,
   ! otherwise all will write to wannier90.wout
#ifdef VASP2WANNIER90
   IF (IO%IU6>=0.AND.WANNIER90()) THEN
   WRITE(IO%IU6,*) 'Calling wannier_setup of wannier90 in library mode'
   WRITE(IO%IU0,*) 'Calling wannier_setup of wannier90 in library mode'
   CALL wannier_setup(MLWF%SEED_NAME,MLWF%mp_grid,MLWF%NKPTS,MLWF%real_lattice,MLWF%recip_lattice, &
  &                   MLWF%kpt_latt,MLWF%NB_TOT,MLWF%num_atoms,MLWF%atom_symbols,MLWF%atoms_cart, &
  &                   MLWF%gamma_only,spinors,MLWF%NNTOT,MLWF%NNLIST,MLWF%NNCELL,num_bands,READ_NUM_WANN, &
  &                   proj_site,proj_l,proj_m,proj_radial,proj_z,proj_x,proj_zona,exclude_bands, &
  &                   proj_s,proj_s_qaxisx)
   ELSE
   IF (IO%IU6>=0) num_bands = MLWF%NB_TOT
   ENDIF
#else
   IF (IO%IU6>=0) num_bands = MLWF%NB_TOT
#endif
   ! Now communicate the output to the other nodes
   CALLMPI( M_sum_i(MLWF%COMM,MLWF%NNTOT,1) )
   CALLMPI( M_sum_i(MLWF%COMM,MLWF%NNLIST,MLWF%NKPTS*num_nnmax) )
   CALLMPI( M_sum_i(MLWF%COMM,MLWF%NNCELL,3*MLWF%NKPTS*num_nnmax) )
   CALLMPI( M_sum_i(MLWF%COMM,num_bands,1) )
   CALLMPI( M_sum_i(MLWF%COMM,READ_NUM_WANN,1) )
   CALLMPI( M_sum_d(MLWF%COMM,proj_site,3*MLWF%NB_TOT) )
   CALLMPI( M_sum_i(MLWF%COMM,proj_l,MLWF%NB_TOT) )
   CALLMPI( M_sum_i(MLWF%COMM,proj_m,MLWF%NB_TOT) )
   CALLMPI( M_sum_i(MLWF%COMM,proj_radial,MLWF%NB_TOT) )
   CALLMPI( M_sum_d(MLWF%COMM,proj_z,3*MLWF%NB_TOT) )
   CALLMPI( M_sum_d(MLWF%COMM,proj_x,3*MLWF%NB_TOT) )
   CALLMPI( M_sum_d(MLWF%COMM,proj_zona,MLWF%NB_TOT) )
   CALLMPI( M_sum_i(MLWF%COMM,exclude_bands,MLWF%NB_TOT) )
   CALLMPI( M_sum_i(MLWF%COMM,proj_s,MLWF%NB_TOT) )
   CALLMPI( M_sum_i(MLWF%COMM,proj_s_qaxisx,3*MLWF%NB_TOT) )

   ALLOCATE(MLWF%EXCLUDE_BAND(MLWF%NB_TOT))
   MLWF%EXCLUDE_BAND=.FALSE. ; NEXCLB=0
   DO N=1,MLWF%NB_TOT
      M=exclude_bands(N)
      IF (M>MLWF%NB_TOT.OR.M<0) THEN
         CALL vtutor%error('MLWF_WANNIER90_SETUP: ERROR: exclude_bands seems corrupt'//str(N)//str(M))
      ELSEIF (M==0) THEN
         EXIT
      ELSE
         MLWF%EXCLUDE_BAND(M)=.TRUE.
         NEXCLB=NEXCLB+1
      ENDIF
   ENDDO
   MLWF%NEXCLB=NEXCLB

   IF ((MLWF%NB_TOT-NEXCLB)/=num_bands) THEN
      CALL vtutor%error("MLWF_WANNIER90_SETUP: ERROR: num_bands seems corrupt " // &
         str(MLWF%NB_TOT-NEXCLB) // " " // str(num_bands))
   ENDIF

   IF (num_bands<num_wann) THEN
      CALL vtutor%error("MLWF_WANNIER90_SETUP: ERROR: num_bands < num_wann " // str(num_bands) // " " // &
         str(num_wann))
   ENDIF

#ifdef debug
   IF (IO%IU6>=0) THEN
!     WRITE(*,*) num_bands_tot,num_bands,num_wann
!     WRITE(*,*) EXCLUDE_BAND
      DO NK=1,MLWF%NKPTS
         WRITE(*,'(I4,3F14.7)') NK,MLWF%kpt_latt(:,NK)
         DO NKP=1,MLWF%NNTOT
         WRITE(*,'(I4,3F14.7,3I4)') MLWF%nnlist(NK,NKP),MLWF%kpt_latt(:,MLWF%nnlist(NK,NKP)),MLWF%nncell(:,NK,NKP)
         ENDDO
         WRITE(*,*)
      ENDDO
   ENDIF
#endif

   ! try to check how many orbitals wanier90 returned
   ESTIMATED_NUM_WANN = COUNT(proj_radial>0.and.proj_radial<4)
   IF (ESTIMATED_NUM_WANN>0.AND.MLWF%PROJ_MODE==PROJ_UNKNOWN) THEN
      IF (ESTIMATED_NUM_WANN/=MLWF%NUM_WANN) THEN
         call vtutor%alert('The wannier interface is returning '//str(ESTIMATED_NUM_WANN)//&
                           ' projections while NUM_WANN is '//str(MLWF%num_wann)//&
                           '. I will continue by using NUM_WANN but you should probably '//&
                           'change NUM_WANN or the wannier90 input file.')
      ENDIF
      MLWF%PROJ_MODE = PROJ_WANNIER90
   ELSE
      RETURN
   ENDIF
   ! transfer what is known about the desired wannier projections onto LOCPROJ array
   ALLOCATE(MLWF%LPRJ_functions(MLWF%NUM_WANN))

   ZONA_WARNING=.TRUE.
   DO IS=1,MLWF%NUM_WANN
      MLWF%LPRJ_functions(IS)%R(1)=MOD(proj_site(1,IS)+100.5_q,1._q)-0.5_q
      MLWF%LPRJ_functions(IS)%R(2)=MOD(proj_site(2,IS)+100.5_q,1._q)-0.5_q
      MLWF%LPRJ_functions(IS)%R(3)=MOD(proj_site(3,IS)+100.5_q,1._q)-0.5_q

      MLWF%LPRJ_functions(IS)%n=proj_radial(IS)
      MLWF%LPRJ_functions(IS)%l=proj_l(IS)
      MLWF%LPRJ_functions(IS)%m=proj_m(IS)

      ! set spinor
      MLWF%LPRJ_functions(IS)%spinor = 1
      IF (proj_s(IS) == -1) MLWF%LPRJ_functions(IS)%spinor = 2
      MLWF%LPRJ_functions(IS)%spin_qaxis = proj_s_qaxisx(3,IS)

      ! set ylmrotation
      MLWF%LPRJ_functions(IS)%LROTYLM=.TRUE.
      MLWF%LPRJ_functions(IS)%PROJ_X = proj_x(:,IS)
      MLWF%LPRJ_functions(IS)%PROJ_Z = proj_z(:,IS)

      ! The default is set so as to keep compatibility with the previous versions
      if (abs(proj_zona(IS)-1.0_q)<1e-8)  then
        MLWF%LPRJ_functions(IS)%za=1.0_q/AUTOA
      else
        IF (ZONA_WARNING) THEN
          call vtutor%alert('The value of zona was manually set in the Wannier90 input file (see wannier90 manual).'//&
                            ' Note that in 5.4.4 < vasp < 6.3 this was assumed to be in Bohr^{-1} wich is not the case (it is defined in Angstroem^{-1}). '//&
                            ' If you want to compare the current results with results produced by one of these versions of vasp '//&
                            ' please multiply the values of zona used on those runs by '//str(1/AUTOA,'(F10.8)'))
          ZONA_WARNING=.FALSE.
        ENDIF
        MLWF%LPRJ_functions(IS)%za=proj_zona(IS)
      endif
      MLWF%LPRJ_functions(IS)%radial_type=3
      MLWF%LPRJ_functions(IS)%WANNIER90_ORBITAL_DEFINITIONS=.TRUE.
   ENDDO

 END SUBROUTINE MLWF_WANNIER90_SETUP

!***********************************************************************
!> @brief Compute .mmn file
!***********************************************************************
 SUBROUTINE MLWF_WANNIER90_MMN(MLWF,WDES,W,P,CQIJ,T_INFO,LATT_CUR,WRITE_MMN,INFO,IO)
   USE base
   USE constant
   USE pseudo
   USE poscar
   USE lattice
   USE wave_high
   USE radial
   USE pead, ONLY : LUSEPEAD, PEAD_CALC_OVERLAP

   TYPE(mlwf_type),INTENT(INOUT) :: MLWF
   TYPE(wavedes),INTENT(IN) :: WDES
   TYPE(wavespin),INTENT(IN) :: W
   TYPE(potcar),INTENT(IN) :: P(:)
   OVERLAP,INTENT(IN) :: CQIJ(:,:,:,:)
   TYPE(type_info),INTENT(IN) :: T_INFO
   TYPE(latt),INTENT(IN) :: LATT_CUR
   LOGICAL,INTENT(IN) :: WRITE_MMN !> Logical controlling whether to write .mmn file to disk
   TYPE(info_struct),INTENT(IN) :: INFO
   TYPE(in_struct),INTENT(IN) :: IO

   ! local variables
   REAL(q) KI(3),KJ(3)
   COMPLEX(q),ALLOCATABLE ::  S(:,:)
   INTEGER :: ISP, N, NBANDS, NI, M, MI, NKI, NKJ

   ! allocation
   ALLOCATE(S(WDES%NB_TOT,WDES%NB_TOT))
   IF (ALLOCATED(MLWF%M_matrix)) DEALLOCATE(MLWF%M_matrix)
   NBANDS = MLWF%NB_TOT-MLWF%NEXCLB
   ALLOCATE(MLWF%M_matrix(NBANDS,NBANDS,MLWF%nntot,MLWF%NKPTS,WDES%ISPIN))
   MLWF%M_matrix = CMPLX(0.0_q,0.0_q,q)

   PROFILING_START( 'mlwf_wannier90_mmn' )

   IF (IO%IU0>=0) WRITE(IO%IU0,*) 'Computing MMN (overlap matrix elements)'
   IF (IO%IU6>=0) WRITE(IO%IU6,*) 'Computing MMN (overlap matrix elements)'

   IF (.NOT.LUSEPEAD()) THEN
      call vtutor%bug('Pead routines must be initialized to compute MMN matrix elements for wannier90',__FILE__,__LINE__)
   ENDIF

   ! loop over spin
   spin: DO ISP=1,WDES%ISPIN
      ! runs over all k-points
      ki_loop: DO NKI=1,MLWF%NKPTS
#ifdef MPI
         IF (MOD(NKI-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
         KI(:)=MLWF%kpt_latt(:,NKI)
         ! runs over the number of nearest-neighbours of k_i
         kj_loop: DO NKJ=1,MLWF%nntot
            KJ(:)=MLWF%kpt_latt(:,MLWF%nnlist(NKI,NKJ))+REAL(MLWF%nncell(:,NKI,NKJ),q)
#ifdef debug
            IF (IO%IU6>=0) WRITE(*,'(2I4,3F14.7,2X,3F14.7)') NKI,NKJ,KI,KJ
#endif
            CALL PEAD_CALC_OVERLAP(W,KI,KJ,ISP,P,CQIJ,LATT_CUR,T_INFO,S,LQIJB=.TRUE.)

            NI=0
            MI=0
            n_loop: DO N=1,WDES%NB_TOT
               IF (MLWF%EXCLUDE_BAND(N)) CYCLE n_loop
               NI=NI+1; MI=0
               m_loop: DO M=1,WDES%NB_TOT
                  IF (MLWF%EXCLUDE_BAND(M)) CYCLE m_loop
                  MI=MI+1
                  MLWF%M_matrix(MI,NI,NKJ,NKI,ISP)=S(M,N)
               ENDDO m_loop
            ENDDO n_loop
            ! check consistency
            IF (NI/=MI.OR.(NI/=MLWF%NB_TOT-MLWF%NEXCLB)) THEN
               CALL vtutor%error("MLWF_WANNIER90_MMN: ERROR: num_bands seems corrupt " // &
                  str(MLWF%NB_TOT) // " " // str(NI))
            ENDIF
         ENDDO kj_loop
      ENDDO ki_loop
   ENDDO spin
   DEALLOCATE(S)
   CALLMPI( M_sum_z( WDES%COMM_KINTER, MLWF%M_matrix, SIZE(MLWF%M_matrix)))

   PROFILING_STOP( 'mlwf_wannier90_mmn' )
 END SUBROUTINE MLWF_WANNIER90_MMN

!***********************************************************************
!> Write the .mmn file.
!***********************************************************************
 SUBROUTINE MLWF_WANNIER90_WRITE_MMN(MLWF,IO)
   USE base
   USE constant
   USE pseudo
   USE poscar
   USE lattice
   USE wave_high
   USE radial
   USE pead, ONLY : PEAD_CALC_OVERLAP

   TYPE(mlwf_type),INTENT(INOUT) :: MLWF
   TYPE(in_struct),INTENT(IN) :: IO

   ! local variables
   REAL(q) KI(3),KJ(3)
   INTEGER :: ISP, N, NBANDS, NI, M, MI, NKI, NKJ

   NBANDS = MLWF%NB_TOT-MLWF%NEXCLB

   IF (IO%IU6<0) RETURN
   WRITE(IO%IU6,*) 'Writing MMN to '//TRIM(MLWF%SEED_NAME)//'.mmn'

   ! loop over spin
   spin: DO ISP=1,MLWF%ISPIN

      IF (MLWF%ISPIN==1) THEN
         OPEN(UNIT=99,FILE=MLWF%SEED_NAME//'.mmn',STATUS='REPLACE')
      ELSE
         OPEN(UNIT=99,FILE=MLWF%SEED_NAME//'.'//str(ISP)//'.mmn',STATUS='REPLACE')
      ENDIF
      WRITE(99,'(A)') 'File generated by VASP: '//MLWF%SZNAM1
      WRITE(99,'(3I12)') NBANDS,MLWF%NKPTS,MLWF%nntot

      ! runs over all k-points
      ki_loop: DO NKI=1,MLWF%NKPTS
         KI(:)=MLWF%kpt_latt(:,NKI)
         ! runs over the number of nearest-neighbours of k_i
         kj_loop: DO NKJ=1,MLWF%nntot
            WRITE(99,'(5I5)') NKI,MLWF%nnlist(NKI,NKJ),MLWF%nncell(:,NKI,NKJ)

            NI=0
            MI=0
            n_loop: DO N=1,MLWF%NB_TOT
               IF (MLWF%EXCLUDE_BAND(N)) CYCLE n_loop
               NI=NI+1; MI=0
               m_loop: DO M=1,MLWF%NB_TOT
                  IF (MLWF%EXCLUDE_BAND(M)) CYCLE m_loop
                  MI=MI+1
                  WRITE(99,'(2F18.12)') MLWF%M_matrix(MI,NI,NKJ,NKI,ISP)
               ENDDO m_loop
            ENDDO n_loop
            ! check consistency
            IF (NI/=MI.OR.(NI/=MLWF%NB_TOT-MLWF%NEXCLB)) THEN
               CALL vtutor%error("MLWF_WANNIER90_MMN: ERROR: num_bands seems corrupt " // &
                  str(MLWF%NB_TOT) // " " // str(NI))
            ENDIF
         ENDDO kj_loop
      ENDDO ki_loop
      CLOSE(99)
   ENDDO spin
   PROFILING_STOP( 'mlwf_wannier90_mmn' )
 END SUBROUTINE MLWF_WANNIER90_WRITE_MMN

!***********************************************************************
!> Compute the .amn projection matrices
!***********************************************************************
 SUBROUTINE MLWF_CALC_AMN(MLWF,WDES,W,P,CQIJ,T_INFO,LATT_CUR,INFO,IO)
   USE base
   USE constant
   USE pseudo
   USE poscar
   USE lattice
   USE wave_high
   USE radial
   USE fileio, ONLY : OUTWAV
   USE umco_pipek_mezey, ONLY : PM_PROJTYPE, PM_PROJTYPE_IAO
   USE locproj, ONLY : LPRJ_function, LPRJ_PROALL, LOCPROJ_IAO,&
                       LPRJ_CALC_INTRINSIC_ATOMIC_ORBITALS

   TYPE(mlwf_type),TARGET,INTENT(INOUT) :: MLWF
   TYPE(wavedes),INTENT(IN) :: WDES
   TYPE(wavespin),INTENT(IN) :: W
   TYPE(potcar),INTENT(IN) :: P(:)
   OVERLAP,INTENT(IN) :: CQIJ(:,:,:,:)
   TYPE(type_info),INTENT(IN) :: T_INFO
   TYPE(latt),INTENT(IN) :: LATT_CUR
   TYPE(info_struct),INTENT(IN) :: INFO
   TYPE(in_struct),INTENT(IN) :: IO

   ! Local variables
   INTEGER :: IW, ISP, IB, IBP
   LOGICAL :: LBUILD_TRIALFUNC
   GDEF,POINTER :: LPRJ_COVL(:,:,:,:)
   !> Contains the intrinsic atomic orbitals as constructed by
   !> the routine CALC_INTRINSIC_ATOMIC_ORBITALS
   TYPE(wavespin) :: W_IAO
   !> Contains the localized orbitals as constructed by
   !> the LPRJ_GEN_TRIALORBITALS routine according to the LOCPROJ tag.
   !> These are needed for the CALC_INTRINSIC_ATOMIC_ORBITALS routine.
   TYPE(wavespin) :: W_TRIALFUNC

   PROFILING_START( 'mlwf_calc_amn' )
   ALLOCATE(MLWF%A_matrix(MLWF%NB_TOT-MLWF%NEXCLB,MLWF%NUM_WANN,MLWF%NKPTS,MLWF%ISPIN))
   MLWF%A_matrix = CMPLX(0.0_q,0.0_q,q)

   IF (MLWF%PROJ_MODE == PROJ_UNKNOWN) RETURN

   IF (IO%IU0>=0) WRITE(IO%IU0,*) 'Computing AMN (projections onto localized orbitals)'
   IF (IO%IU6>=0) WRITE(IO%IU6,*) 'Computing AMN (projections onto localized orbitals)'

   IF (.NOT.ALLOCATED(MLWF%LPRJ_functions)) &
      CALL vtutor%bug('Projections onto localized orbitals are not present',__FILE__,__LINE__)

   ! Decide whether we need to keep the wavespin object with the locproj orbitals
   LBUILD_TRIALFUNC=LOCPROJ_IAO.OR.LWAVE_IAO.OR.LWAVE_LOCPROJ.OR.LWRITE_LOCPROJ_XSF.OR.LWRITE_IAO_XSF

   ! compute LPRJ_PROALL
   IF (LBUILD_TRIALFUNC) THEN
       CALL LPRJ_PROALL(W,WDES,P,CQIJ,LATT_CUR,MLWF%LPRJ_functions,MLWF%LPRJ_COVL,&
                        T_INFO,INFO,IO,VKPT=MLWF%kpt_latt,W_TRIALFUNC=W_TRIALFUNC)
   ELSE
       CALL LPRJ_PROALL(W,WDES,P,CQIJ,LATT_CUR,MLWF%LPRJ_functions,MLWF%LPRJ_COVL,&
                        T_INFO,INFO,IO,VKPT=MLWF%kpt_latt)
   ENDIF
   LPRJ_COVL=>MLWF%LPRJ_COVL

   ! calc IAOs if desired
   IF (LOCPROJ_IAO.OR.PM_PROJTYPE==PM_PROJTYPE_IAO) THEN
     CALL LPRJ_CALC_INTRINSIC_ATOMIC_ORBITALS(W, W_TRIALFUNC, MLWF%LPRJ_COVL, MLWF%LPRJ_COVL_IAO, W_IAO, LATT_CUR, IO, INFO, CQIJ)
     LPRJ_COVL=>MLWF%LPRJ_COVL_IAO
     IF (LWAVE_IAO) THEN
        IF (IO%IU0>=0) WRITE(IO%IU0,*) 'Writing WAVECAR.iao containing IAOs (intrinsic atomic orbitals)'
        CALL OUTWAV(IO, W_IAO%WDES, W_IAO, LATT_CUR, 0.0_q, 'iao')
     ENDIF
     IF (LWRITE_IAO_XSF) THEN
        IF (IO%IU0>=0) WRITE(IO%IU0,*) "Write IAO XSF files"
        CALL WRITE_WAVESPIN_XSF('iao',W_IAO,LATT_CUR,P,MLWF%num_atoms,&
        MLWF%atom_symbols,MLWF%atoms_cart,WRITE_IAO_XSF_IW)
     ENDIF
   ENDIF

   ! If writting the localized orbitals was requested
   IF (LBUILD_TRIALFUNC) THEN
     IF (LWAVE_LOCPROJ) THEN
        IF (IO%IU0>=0) WRITE(IO%IU0,*) 'Writing WAVECAR.loc containing Localized orbitals specified by LOCPROJ'
        CALL OUTWAV(IO, W_TRIALFUNC%WDES, W_TRIALFUNC, LATT_CUR, 0.0_q, 'loc')
     ENDIF
     IF (LWRITE_LOCPROJ_XSF) THEN
        IF (IO%IU0>=0) WRITE(IO%IU0,*) "Write LOCPROJ XSF files"
        CALL WRITE_WAVESPIN_XSF('locproj',W_TRIALFUNC,LATT_CUR,P,MLWF%num_atoms,&
        MLWF%atom_symbols,MLWF%atoms_cart,WRITE_LOCPROJ_XSF_IW)
     ENDIF
     ! no need for W_TRIALFUNC anymore
     CALL DEALLOCW(W_TRIALFUNC)
   ENDIF

   ! Need to transfer the data to what wannier90 expects
   DO IW = 1, MLWF%NUM_WANN
      IBP = 0
      DO IB = 1, MLWF%NB_TOT
         IF (MLWF%EXCLUDE_BAND(IB)) CYCLE
         IBP = IBP + 1
         MLWF%A_matrix(IBP,IW,:,:) = LPRJ_COVL(IB,:,:,IW)
      ENDDO
   ENDDO
   PROFILING_STOP( 'mlwf_calc_amn' )

 END SUBROUTINE MLWF_CALC_AMN

!***********************************************************************
!> @brief writes A_{mn}^k tranformation matrix between
!  Bloch and localized orbitals functions to file.
!  It consists of projections of the Bloch orbitals onto localized states
!
!  u^w_{nk} = sum_m A_{mn}^k u^b_{mk}
!
!  where:
!  u^w_{nk} is the periodic part of the localized orbitals
!  u^n_{mk} is the periodic part of the bloch orbitals
!
!  The format used is the one of the .amn files of wannier90
!***********************************************************************
 SUBROUTINE WRITE_AMN_FILE(FILENAME,NBANDS,NPROJ,NKPTS,AMN,SZNAM1,IO)
   USE base
   CHARACTER(LEN=*),INTENT(IN) :: FILENAME !> filename in which to write AMN
   INTEGER,INTENT(IN) :: NBANDS !> number of Bloch states
   INTEGER,INTENT(IN) :: NPROJ  !> number of localized states
   INTEGER,INTENT(IN) :: NKPTS  !> number of k-points
   COMPLEX(q),INTENT(IN) :: AMN(:,:,:) !> A_{mn}^k matrix in the form AMN(NBANDS,NPROJ,NKPTS)
   CHARACTER(LEN=*),INTENT(IN) :: SZNAM1 !> name of the system
   TYPE(in_struct),INTENT(IN) :: IO
   !local variables
   INTEGER :: NKI,N,M,UNIT

   WRITE(IO%IU6,*) 'Writing AMN to '//FILENAME//'.amn'

   OPEN(NEWUNIT=UNIT,FILE=FILENAME//'.amn',STATUS='REPLACE')
   WRITE(UNIT,'(A)') 'File generated by VASP: '//SZNAM1
   WRITE(UNIT,'(3I12)') NBANDS,NKPTS,NPROJ
   DO NKI=1,NKPTS
     DO N=1,NPROJ
       DO M=1,NBANDS
         WRITE(UNIT,'(3I5,2F18.12)') M,N,NKI,AMN(M,N,NKI)
       ENDDO
     ENDDO
   ENDDO
   CLOSE(UNIT)
 END SUBROUTINE WRITE_AMN_FILE

 !***********************************************************************
 !> @brief Write .amn file to disk
 !***********************************************************************
 SUBROUTINE MLWF_WANNIER90_WRITE_AMN(MLWF,IO)
   USE base
   USE string
   TYPE(mlwf_type),INTENT(INOUT) :: MLWF
   TYPE(in_struct),INTENT(IN) :: IO
   !local variables
   INTEGER :: ISP

   IF (IO%IU6>=0) THEN
      DO ISP=1,MLWF%ISPIN
         IF (MLWF%ISPIN==1) THEN
            CALL WRITE_AMN_FILE(MLWF%SEED_NAME,MLWF%NB_TOT-MLWF%NEXCLB,MLWF%NUM_WANN,MLWF%NKPTS,MLWF%A_matrix(:,:,:,ISP),MLWF%SZNAM1,IO)
         ELSE
            CALL WRITE_AMN_FILE(MLWF%SEED_NAME//'.'//str(ISP),MLWF%NB_TOT-MLWF%NEXCLB,MLWF%NUM_WANN,MLWF%NKPTS,MLWF%A_matrix(:,:,:,ISP),MLWF%SZNAM1,IO)
         ENDIF
      ENDDO
   ENDIF

 END SUBROUTINE MLWF_WANNIER90_WRITE_AMN

!***********************************************************************
!> @brief Write .eig file to be used by Wannier90
!***********************************************************************
 SUBROUTINE MLWF_WANNIER90_EIG(MLWF,W,WRITE_EIG,IO)
    USE base
    USE wave_high
    USE full_kpoints
    USE tutor, ONLY : vtutor

    TYPE(mlwf_type),INTENT(INOUT) :: MLWF
    TYPE(wavespin) :: W
    LOGICAL,INTENT(IN) :: WRITE_EIG !> Logical controlling whether to write .eig file to disk
    TYPE(in_struct),INTENT(IN) :: IO

    !local variables
    INTEGER :: ISP, NK, NKI, M, MI

    ! Check consistency
    CALL CHECK_FULL_KPOINTS
    IF (MLWF%NB_TOT /= W%WDES%NB_TOT) call vtutor%error('Inconsistent number of bounds found')

    ALLOCATE(MLWF%eigenvalues(MLWF%NB_TOT,MLWF%NKPTS,MLWF%ISPIN))

    IF (IO%IU6>=0) WRITE(IO%IU6,*) 'Writing EIG to '//TRIM(MLWF%SEED_NAME)//'.eig'

    DO ISP=1,MLWF%ISPIN
       IF (IO%IU6>=0.AND.WRITE_EIG) THEN
          IF (MLWF%ISPIN==1) THEN
             OPEN(UNIT=99,FILE=MLWF%SEED_NAME//'.eig',STATUS='REPLACE')
          ELSE
             OPEN(UNIT=99,FILE=MLWF%SEED_NAME//'.'//str(ISP)//'.eig',STATUS='REPLACE')
          ENDIF
       ENDIF

       DO NKI=1,MLWF%NKPTS
          NK=MLWF%NEQUIV(NKI)
          MI=0
          DO M=1,MLWF%NB_TOT
             IF (MLWF%EXCLUDE_BAND(M)) CYCLE
             MI=MI+1
             MLWF%eigenvalues(MI,NKI,ISP)=REAL(W%CELTOT(M,NK,ISP))
             IF (IO%IU6>0.AND.WRITE_EIG) WRITE(99,'(2I12,F22.12)') MI,NKI,REAL(W%CELTOT(M,NK,ISP))
          ENDDO
       ENDDO
       IF (IO%IU6>=0.AND.WRITE_EIG) CLOSE(99)
    ENDDO
 END SUBROUTINE MLWF_WANNIER90_EIG

!> @brief Write .win input file for Wannier90
 SUBROUTINE MLWF_WANNIER90_WIN(MLWF,LSETUP)
    USE string, ONLY : lowercase
    TYPE(mlwf_type),INTENT(INOUT) :: MLWF
    LOGICAL, OPTIONAL :: LSETUP
    ! local variables
    INTEGER :: NK, NI, IERR
    CHARACTER(LEN=WANNIER90_WIN_MAXLEN) :: WANNIER90_WIN_FILTERED
    ! for checking that the file was written to disk
    INTEGER,PARAMETER :: MAX_LINE_LENGTH = 255
    CHARACTER(LEN=MAX_LINE_LENGTH) :: BUFLIN
    LOGICAL  :: NOCCUR
    EXTERNAL :: NOCCUR

    CALL FILTER_WIN_FILE(MLWF,WANNIER90_WIN,WANNIER90_WIN_FILTERED)
    OPEN(UNIT=99,FILE=MLWF%SEED_NAME//'.win',STATUS='REPLACE')
    WRITE(99,'(A)') TRIM(WANNIER90_WIN_FILTERED)
    WRITE(99,'(A)') '# This part was generated automatically by VASP'

    ! write number of Bloch orbitals
    WRITE(99,'(A,I0)') 'num_bands = ', MLWF%NB_TOT-MLWF%NEXCLB

    ! write number of Wannier orbitals
    WRITE(99,'(A,I0)') 'num_wann = ', MLWF%NUM_WANN

    ! write bloch phases
    IF ((.NOT.PRESENT(LSETUP)).AND.(MLWF%NUM_WANN==0.OR.MLWF%PROJ_MODE==PROJ_UNKNOWN)) THEN
       WRITE(99, '(A)') 'use_bloch_phases = .true.'
    ENDIF

    ! write auto projections
    !IF (LSCDM) WRITE(99,'(A)') 'auto_projections = .true.'

    ! write spinors
    IF (MLWF%NRSPINORS==2) WRITE(99,'(A)') 'spinors = .true.'

    ! automatic disentanglement window
    IF (MLWF%DIS_WIN_MIN/=HUGE(MLWF%DIS_WIN_MIN)) THEN
       WRITE(99,'(A,F14.7)') 'dis_win_min =', MLWF%DIS_WIN_MIN
    ENDIF
    IF (MLWF%DIS_WIN_MAX/=HUGE(MLWF%DIS_WIN_MAX)) THEN
       WRITE(99,'(A,F14.7)') 'dis_win_max =', MLWF%DIS_WIN_MAX
    ENDIF
    IF (MLWF%DIS_FROZ_MIN/=HUGE(MLWF%DIS_FROZ_MIN)) THEN
       WRITE(99,'(A,F14.7)') 'dis_froz_min =', MLWF%DIS_FROZ_MIN
    ENDIF
    IF (MLWF%DIS_FROZ_MAX/=HUGE(MLWF%DIS_FROZ_MAX)) THEN
       WRITE(99,'(A,F14.7)') 'dis_froz_max =', MLWF%DIS_FROZ_MAX
    ENDIF

    ! write unit cell
    WRITE(99,'(A)') 'begin unit_cell_cart'
    WRITE(99,'(3F14.7)') MLWF%real_lattice(1,:)
    WRITE(99,'(3F14.7)') MLWF%real_lattice(2,:)
    WRITE(99,'(3F14.7)') MLWF%real_lattice(3,:)
    WRITE(99,'(A)')  'end unit_cell_cart'

    ! write atoms
    WRITE(99,'(A)') 'begin atoms_cart'
    DO NI=1,MLWF%num_atoms
       WRITE(99,'(A2,2X,3F14.7)') MLWF%atom_symbols(NI),MLWF%atoms_cart(:,NI)
    ENDDO
    WRITE(99,'(A)')  'end atoms_cart'

    ! write kpoints
    WRITE(99,'(A,3I6)') 'mp_grid =',MLWF%mp_grid
    WRITE(99,'(A)') 'begin kpoints'
    DO NK=1,MLWF%NKPTS
       WRITE(99,'(3F20.12)') MLWF%kpt_latt(:,NK)
    ENDDO
    WRITE(99,'(A)')  'end kpoints'
    CLOSE(99)

 END SUBROUTINE MLWF_WANNIER90_WIN

!***********************************************************************
!> @brief Call to Wannier90 interface
!***********************************************************************
 SUBROUTINE MLWF_WANNIER90_RUN(MLWF,IO)
    USE base
    USE string, ONLY : str

    TYPE(mlwf_type),INTENT(INOUT) :: MLWF
    TYPE(in_struct) :: IO

    ! local variables
    INTEGER :: ISP

    ALLOCATE(MLWF%U_matrix(MLWF%NUM_WANN,MLWF%NUM_WANN,MLWF%NKPTS,MLWF%ISPIN))
    ALLOCATE(MLWF%U_matrix_opt(MLWF%NB_TOT-MLWF%NEXCLB,MLWF%NUM_WANN,MLWF%NKPTS,MLWF%ISPIN))
    ALLOCATE(MLWF%lwindow(MLWF%NB_TOT-MLWF%NEXCLB,MLWF%NKPTS,MLWF%ISPIN))
    ALLOCATE(MLWF%wann_centres(3,MLWF%NUM_WANN,MLWF%ISPIN))
    ALLOCATE(MLWF%wann_spreads(MLWF%NUM_WANN,MLWF%ISPIN))
    ALLOCATE(MLWF%spread(3,MLWF%ISPIN))

    MLWF%U_matrix=0._q; MLWF%U_matrix_opt=0._q; MLWF%lwindow=.FALSE.
    MLWF%wann_centres=0._q; MLWF%wann_spreads=0._q; MLWF%spread=0._q

#ifdef VASP2WANNIER90
    ! Execute wannier_run only on the masternode
    IF (IO%IU6>=0) THEN
       IF (MLWF%ISPIN==1) THEN
          WRITE(IO%IU6,*) 'Calling wannier_run of wannier90 in library mode (check wannier90.wout)'
          WRITE(IO%IU0,*) 'Calling wannier_run of wannier90 in library mode (check wannier90.wout)'
          CALL wannier_run(MLWF%seed_name,MLWF%mp_grid,MLWF%NKPTS,MLWF%real_lattice,MLWF%recip_lattice, &
         &                 MLWF%kpt_latt,MLWF%NB_TOT-MLWF%NEXCLB,MLWF%NUM_WANN,MLWF%nntot,MLWF%num_atoms,MLWF%atom_symbols, &
         &                 MLWF%atoms_cart,MLWF%gamma_only,&
         &                 MLWF%M_matrix(:,:,:,:,1),MLWF%A_matrix(:,:,:,1),MLWF%eigenvalues(:,:,1), &
         &                 MLWF%U_matrix(:,:,:,1),MLWF%U_matrix_opt(:,:,:,1),MLWF%lwindow(:,:,1), &
         &                 MLWF%wann_centres(:,:,1),MLWF%wann_spreads(:,1),MLWF%spread(:,1))
      ELSE
          DO ISP=1,MLWF%ISPIN

          IF (ISP==1) WRITE(*,*) 'Running wannier90 in library mode for spin up   (check wannier90.wout)'
          IF (ISP==2) WRITE(*,*) 'Running wannier90 in library mode for spin down (check wannier90.wout)'

          ! copy wannier90.win to wannier90.(up/dn).win
          CALL COPYFILE(MLWF%seed_name//'.win',MLWF%seed_name//'.'//str(ISP)//'.win')

          ! run wannier90 for spin component ISP
          CALL wannier_run(MLWF%seed_name//'.'//str(ISP), &
         &                 MLWF%mp_grid,MLWF%NKPTS,MLWF%real_lattice,MLWF%recip_lattice, &
         &                 MLWF%kpt_latt,MLWF%NB_TOT-MLWF%NEXCLB,MLWF%NUM_WANN,MLWF%nntot,MLWF%num_atoms,MLWF%atom_symbols, &
         &                 MLWF%atoms_cart,MLWF%gamma_only,&
         &                 MLWF%M_matrix(:,:,:,:,ISP),MLWF%A_matrix(:,:,:,ISP),MLWF%eigenvalues(:,:,ISP), &
         &                 MLWF%U_matrix(:,:,:,ISP),MLWF%U_matrix_opt(:,:,:,ISP),MLWF%lwindow(:,:,ISP), &
         &                 MLWF%wann_centres(:,:,ISP),MLWF%wann_spreads(:,ISP),MLWF%spread(:,ISP))
          ENDDO
       ENDIF
    ENDIF
#endif

    ! Communicate the results to the other nodes
    CALLMPI( M_sum_z(MLWF%COMM,MLWF%U_matrix,MLWF%NUM_WANN*MLWF%NUM_WANN*MLWF%NKPTS*MLWF%ISPIN) )
    CALLMPI( M_sum_z(MLWF%COMM,MLWF%U_matrix_opt,MLWF%NB_TOT*MLWF%NUM_WANN*MLWF%NKPTS*MLWF%ISPIN) )
    CALLMPI( M_bcast_l(MLWF%COMM,MLWF%lwindow,MLWF%NB_TOT*MLWF%NKPTS*MLWF%ISPIN) )
    CALLMPI( M_sum_d(MLWF%COMM,MLWF%wann_centres,3*MLWF%NUM_WANN*MLWF%ISPIN) )
    CALLMPI( M_sum_d(MLWF%COMM,MLWF%wann_spreads,MLWF%NUM_WANN*MLWF%ISPIN) )
    CALLMPI( M_sum_d(MLWF%COMM,MLWF%spread,3*MLWF%ISPIN) )

    ! Set trafo method in mlwf structure
    MLWF%TRAFO_MODE = TRAFO_WANNIER90_MODE

 END SUBROUTINE MLWF_WANNIER90_RUN

!***********************************************************************
!> @brief Call to TRAFO routine to compute one-shot Wannier functions
!***********************************************************************
 SUBROUTINE MLWF_TRAFO_RUN(MLWF,W,P,CQIJ,LATT_CUR,T_INFO,INFO,IO)
    USE base
    USE lattice, ONLY: latt
    USE poscar_struct_def, ONLY: type_info
    USE pseudo_struct_def, ONLY: potcar
    USE wave_high, ONLY: wavespin
    USE umco_pipek_mezey, ONLY : PM_PROJTYPE
    USE umco, ONLY : IS_UMCO, UMCO_CALC_TRAFO

    TYPE(mlwf_type),INTENT(INOUT) :: MLWF
    TYPE(wavespin) :: W
    TYPE(potcar) :: P(:)
    OVERLAP      :: CQIJ(:,:,:,:)
    TYPE(latt) :: LATT_CUR
    TYPE(type_info) :: T_INFO
    TYPE(info_struct) :: INFO
    TYPE(in_struct) :: IO

    ! Local variables
    INTEGER :: IK,IS,IW
    GDEF,ALLOCATABLE :: T(:,:,:,:)

    PROFILING_START('mlwf_trafo_run')
    ALLOCATE(MLWF%U_matrix(MLWF%NUM_WANN,MLWF%NUM_WANN,MLWF%NKPTS,MLWF%ISPIN))
    ALLOCATE(MLWF%U_matrix_opt(MLWF%NB_TOT,MLWF%NUM_WANN,MLWF%NKPTS,MLWF%ISPIN))
    ALLOCATE(MLWF%lwindow(MLWF%NB_TOT,MLWF%NKPTS,MLWF%ISPIN))
    MLWF%lwindow = .TRUE.
    MLWF%U_matrix = CMPLX(0.0_q,0.0_q,KIND=q)
    DO IS=1,MLWF%ISPIN
       DO IK=1,MLWF%NKPTS
          DO IW=1,MLWF%NUM_WANN
             MLWF%U_matrix(IW,IW,IK,IS) = CMPLX(1.0_q,0.0_q,KIND=q)
          ENDDO
       ENDDO
    ENDDO

    SELECT CASE (MLWF%PROJ_MODE)
        ! In case the projections were generated with LOCPROJ we need to run TRAFO
        CASE(PROJ_LOCPROJ,PROJ_LOCPROJAUTO,PROJ_WANNIER90)
            IF (IS_UMCO()) THEN
                ! UMCO contributed by Tobias Schafer
                SELECT CASE(PM_PROJTYPE)
                    CASE (1) ! intrinsic atomic orbitals
                        CALL UMCO_CALC_TRAFO(MLWF%LPRJ_COVL_IAO, T, W, W%WDES%GRID, P, CQIJ, LATT_CUR, T_INFO, INFO, IO)
                    CASE (2) ! trial orbitals as defined by LOCPROJ
                        CALL UMCO_CALC_TRAFO(MLWF%LPRJ_COVL, T, W, W%WDES%GRID, P, CQIJ, LATT_CUR, T_INFO, INFO, IO)
                END SELECT
                ! TODO: the number of wannier orbitals for UMCO should be computed beforehand
                MLWF%num_wann=SIZE(T,1)
                MLWF%U_matrix_opt=CMPLX(0.0_q,0.0_q,KIND=q)
                MLWF%U_matrix_opt(:MLWF%num_wann,:MLWF%num_wann,:,:)=T(:,:,:,:)
                MLWF%TRAFO_MODE = TRAFO_UMCO
            ELSE
                ! One-shot SVD method
                CALL CALC_WANNIER_TRAFO_PROJ(T,W,MLWF)
                MLWF%U_matrix_opt=T
                MLWF%TRAFO_MODE = TRAFO_ONESHOTSVD
            ENDIF
            DEALLOCATE(T)
        ! In case SCDM was used, the projections are already unitary
        CASE(PROJ_SCDM)
            MLWF%U_matrix_opt=MLWF%A_matrix
            MLWF%TRAFO_MODE = TRAFO_SCDM
    END SELECT
    PROFILING_STOP('mlwf_trafo_run')

 END SUBROUTINE MLWF_TRAFO_RUN

!***********************************************************************
!>@brief Write WANPROJ file
!***********************************************************************
 SUBROUTINE MLWF_WRITE_WANPROJ(MLWF, WANPROJFILE, IO)
    USE base
    USE wave_high
    USE locproj, ONLY : WRITE_WANPROJ_FILE

    TYPE(mlwf_type),INTENT(INOUT) :: MLWF
    CHARACTER(LEN=*) :: WANPROJFILE
    TYPE(in_struct),INTENT(IN) :: IO
    ! local variables
    GDEF, ALLOCATABLE :: U(:,:,:,:)

    PROFILING_START('mlwf_write_wanproj')
    IF (IO%IU0>0) WRITE(IO%IU0,*)' Writing '//WANPROJFILE//' file'
    IF (IO%IU6>0) WRITE(IO%IU6,*)' Writing '//WANPROJFILE//' file'
    CALL MLWF_GET_U(MLWF,MLWF%NKPTS,MLWF%kpt_latt,U)
    CALL WRITE_WANPROJ_FILE(WANPROJFILE,MLWF%NB_TOT,MLWF%ISPIN,MLWF%NKPTS,MLWF%NUM_WANN,MLWF%kpt_latt,U)
    DEALLOCATE(U)
    PROFILING_STOP('mlwf_write_wanproj')

 END SUBROUTINE MLWF_WRITE_WANPROJ

!
!> Write WAVECAR with wannier orbitals computed by the MLWF module
!
 SUBROUTINE MLWF_WRITE_WAVECAR(SELF,W,LATT_CUR,IO)
   use poscar, ONLY : latt
   use base, ONLY : in_struct
   use wave_struct_def, ONLY : wavespin, wavedes
   use fileio, ONLY : OUTWAV
   use tutor, ONLY : vtutor
   use string, ONLY : str
   use locproj, ONLY : LPRJ_LINCOM
   use wave, ONLY : DEALLOCW, ALLOCW
   TYPE(mlwf_type),INTENT(IN) :: SELF
   TYPE(wavespin),INTENT(IN) :: W !< KS bloch states computed by VASP
   TYPE(latt),INTENT(IN) :: LATT_CUR !< Lattice datastructure
   TYPE(in_struct),INTENT(IN) :: IO
   ! local
   TYPE(wavespin) :: W_WANNIER
   GDEF,ALLOCATABLE :: U(:,:,:,:)
   INTEGER :: NPAR,NBANDS

   PROFILING_START('mlwf_write_wavecar')

   IF (SELF%num_wann>W%WDES%NB_TOT) THEN
      call vtutor%alert('Number of wannier states NUM_WANN='//str(SELF%NUM_WANN) //' is larger than number of bands NBANDS='//str(W%WDES%NB_TOT)//&
                        '. Writting the WAVECAR is not implemented in this case. Re-run the calculation with NBANDS>'//str(SELF%NUM_WANN))
      RETURN
   ENDIF

   ! create a temporary array with transformation matrix
   CALL ALLOCW(W%WDES,W_WANNIER)
   CALL MLWF_GET_U(SELF,SELF%NKPTS,SELF%kpt_latt,U)
   CALL LPRJ_LINCOM(W,U,W_WANNIER)
   CALL OUTWAV(IO, W_WANNIER%WDES, W_WANNIER, LATT_CUR, 0.0_q, 'wan')
   CALL DEALLOCW(W_WANNIER)
   PROFILING_STOP('mlwf_write_wavecar')
 END SUBROUTINE MLWF_WRITE_WAVECAR

!***********************************************************************
!> @brief Write UNK for Wannier90 or XSF files
!***********************************************************************
 SUBROUTINE MLWF_WRITE_UNK_XSF(MLWF, W, P, LATT_CUR, WRITE_UNK, WRITE_WANNIER_XSF, IO)
    USE base
    USE pseudo
    USE poscar
    USE lattice
    USE wave_high

    TYPE(mlwf_type),INTENT(INOUT) :: MLWF
    TYPE(wavespin) :: W
    TYPE(potcar),INTENT(IN) :: P(:)
    TYPE(latt),INTENT(IN) :: LATT_CUR
    LOGICAL,INTENT(IN) :: WRITE_UNK
    LOGICAL,INTENT(IN) :: WRITE_WANNIER_XSF
    TYPE(in_struct),INTENT(IN) :: IO

    ! local variables
    INTEGER :: NGX,NGY,NGZ,NX,NY,NZ,NXI,NYI,NZI,IND
    INTEGER :: NKI,NI,IW,ISP,ISPINOR,IB, NPWS
    INTEGER :: I
    INTEGER :: UNKUNIT, WANNUNIT
    LOGICAL :: LSKIP
    CHARACTER(LEN=11):: UNKFILE
    REAL(q) :: KI(3)
    GDEF, ALLOCATABLE :: U(:,:,:,:)
    GDEF,ALLOCATABLE :: WVFN_BANDS(:,:)
    GDEF,ALLOCATABLE :: WVFN_WANN(:,:,:)
    COMPLEX(q) :: UP,DW
    CHARACTER(LEN=50) :: FILENAME

    NGX=W%WDES%GRID%NGX
    NGY=W%WDES%GRID%NGY
    NGZ=W%WDES%GRID%NGZ

    PROFILING_START('mlwf_write_unk_xsf')
    IF (IO%IU0>=0.AND.WRITE_UNK)         WRITE(IO%IU0,*) "Write UNK files"
    IF (IO%IU0>=0.AND.WRITE_WANNIER_XSF) WRITE(IO%IU0,*) "Write wannier XSF files"
    IF (WRITE_WANNIER_XSF) THEN
       CALL MLWF_GET_U(MLWF,MLWF%NKPTS,MLWF%kpt_latt,U)
       ALLOCATE(WVFN_BANDS(NGX*NGY*NGZ*W%WDES%NRSPINORS,W%WDES%NB_TOT))
       ALLOCATE(WVFN_WANN(NGX*NGY*NGZ*W%WDES%NRSPINORS,MLWF%NUM_WANN,W%WDES%ISPIN))
       WVFN_WANN = 0
    ENDIF
    ! loop over spin
    DO ISP=1,MLWF%ISPIN
       DO NKI=1,MLWF%NKPTS
          KI = MLWF%kpt_latt(:,NKI)
          UNKUNIT = -1
          IF (IO%IU6>=0.AND.WRITE_UNK) THEN
             IF (W%WDES%NRSPINORS==2) THEN
                WRITE(UNKFILE,"('UNK',I5.5,'.NC')") NKI
             ELSE
                WRITE(UNKFILE,"('UNK',I5.5,'.',I1)") NKI,ISP
             ENDIF
             UNKUNIT=99
             OPEN(UNIT=UNKUNIT,FILE=UNKFILE,FORM='UNFORMATTED',STATUS='REPLACE')
             WRITE(UNKUNIT) NGX,NGY,NGZ,NKI,MLWF%NB_TOT
          ENDIF
          IF (WRITE_WANNIER_XSF) THEN
             CALL GET_WAVE_FUNCTIONS(W,KI,ISP,P,LATT_CUR,UNKUNIT,WVFN_BANDS)
             DO IW=1,MLWF%NUM_WANN
                DO IB=1,W%WDES%NB_TOT
                   IF (MLWF%EXCLUDE_BAND(IB)) CYCLE
                   WVFN_WANN(:,IW,ISP) = WVFN_WANN(:,IW,ISP)+WVFN_BANDS(:,IB)*U(IB,IW,NKI,ISP)
                ENDDO
             ENDDO
          ELSE
             CALL GET_WAVE_FUNCTIONS(W,KI,ISP,P,LATT_CUR,UNKUNIT,EXCLUDE_BAND=MLWF%EXCLUDE_BAND)
          ENDIF
          IF (IO%IU6>=0.AND.WRITE_UNK) CLOSE(UNKUNIT)
       ENDDO !end loop k-points
    ENDDO
    IF (WRITE_WANNIER_XSF) THEN
       ! Normalize WVFN_WANN
       WVFN_WANN = WVFN_WANN/MLWF%NKPTS
       ! Free some memory
       DEALLOCATE(WVFN_BANDS)
       DO ISP=1,MLWF%ISPIN
        DO IW=1,MLWF%NUM_WANN
           ! Only write the orbitals that were selected
           IF (ALLOCATED(WRITE_WANNIER_XSF_IW)) THEN
              LSKIP=.TRUE.
              DO I=1,SIZE(WRITE_WANNIER_XSF_IW)
                 IF (IW==WRITE_WANNIER_XSF_IW(I)) LSKIP=.FALSE.
              ENDDO
              IF (LSKIP) CYCLE
           ENDIF
           ! write XSF file from WVFN_WANNIER
           DO ISPINOR=1,W%WDES%NRSPINORS
              NPWS = (NGX*NGY*NGZ)*(ISPINOR-1)
              FILENAME = XSF_FILENAME('wannier',IW,ISP,ISPINOR,W%WDES%NRSPINORS,W%WDES%ISPIN)
              CALL WRITE_XSF(FILENAME,TRANSPOSE(MLWF%real_lattice),MLWF%num_atoms,&
                             MLWF%atom_symbols,MLWF%atoms_cart,NGX,NGY,NGZ,WVFN_WANN(NPWS+1:,IW,ISP))
           ENDDO ! loop over spinors
        ENDDO ! loop over wannier states
       ENDDO ! loop over spin states
       DEALLOCATE(WVFN_WANN)
    ENDIF
    PROFILING_STOP('mlwf_write_unk_xsf')
 END SUBROUTINE MLWF_WRITE_UNK_XSF

!
!> Build xsf filename based on the prefix spinor or spin index
!> spin refers to a spin polarized calculation (ISPIN=2)
!> spinor to a non-collinear calculation (LNONCOLLINEAR=.TRUE.)
!> even though physically speaking there is not distinction
!
 FUNCTION XSF_FILENAME(PREFIX,IB,ISP,ISPINOR,NRSPINORS,ISPIN) RESULT(FILENAME)
    use string, ONLY: str
    CHARACTER(LEN=*) :: PREFIX
    INTEGER :: IB !< index of the state
    INTEGER :: ISP !< spin index
    INTEGER :: ISPINOR !< spinor index
    INTEGER :: NRSPINORS !< number of spinors
    INTEGER :: ISPIN !< number of spin channels
    CHARACTER(LEN=50) :: FILENAME
    IF (NRSPINORS==1.AND.ISPIN==1) THEN
       WRITE(FILENAME,'(A)') TRIM(PREFIX)//"_"//str(IB)//".xsf"
    ELSEIF (ISPIN==2) THEN
       WRITE(FILENAME,'(A)') TRIM(PREFIX)//"_"//str(IB)//"_spin_"//str(ISP)//".xsf"
    ELSE
       WRITE(FILENAME,'(A)') TRIM(PREFIX)//"_"//str(IB)//"_spinor_"//str(ISPINOR)//".xsf"
    ENDIF
 END FUNCTION

!
!> Write a XSF file from real space information
!
 SUBROUTINE WRITE_XSF(FILENAME,LAT,NATOMS,ATOM_SYMBOLS,ATOM_CART_POS,NGX,NGY,NGZ,DAT)
    CHARACTER(LEN=*),INTENT(IN) :: FILENAME
    REAL(q),INTENT(IN) :: LAT(3,3)
    INTEGER,INTENT(IN) :: NATOMS
    CHARACTER(LEN=*),INTENT(IN) :: ATOM_SYMBOLS(NATOMS)
    REAL(q),INTENT(IN) :: ATOM_CART_POS(3,NATOMS)
    INTEGER,INTENT(IN) :: NGX,NGY,NGZ
    GDEF,INTENT(IN) :: DAT(:)
    ! local variables
    INTEGER :: NX,NY,NZ
    INTEGER :: NXI,NYI,NZI
    INTEGER :: NI,IND
    INTEGER :: WANNUNIT

    OPEN(NEWUNIT=WANNUNIT,FILE=TRIM(FILENAME),FORM='FORMATTED',STATUS='REPLACE')
    WRITE(WANNUNIT,*) '# XSF File generated by VASP'
    WRITE(WANNUNIT,*) 'CRYSTAL'
    WRITE(WANNUNIT,*) 'PRIMVEC'
    WRITE(WANNUNIT,'(3F14.7)') LAT(:,1)
    WRITE(WANNUNIT,'(3F14.7)') LAT(:,2)
    WRITE(WANNUNIT,'(3F14.7)') LAT(:,3)
    WRITE(WANNUNIT,*) 'CONVVEC'
    WRITE(WANNUNIT,'(3F14.7)') LAT(:,1)
    WRITE(WANNUNIT,'(3F14.7)') LAT(:,2)
    WRITE(WANNUNIT,'(3F14.7)') LAT(:,3)
    WRITE(WANNUNIT,*) 'PRIMCOORD'
    WRITE(WANNUNIT,*) NATOMS,1
    DO NI=1,NATOMS
       WRITE(WANNUNIT,'(A2,2X,3F14.7)') ATOM_SYMBOLS(NI),ATOM_CART_POS(:,NI)
    ENDDO
    ! write DATAGRID 3D
    WRITE(WANNUNIT,*)
    WRITE(WANNUNIT,*) 'BEGIN_BLOCK_DATAGRID_3D'
    WRITE(WANNUNIT,*) '3D_field'
    WRITE(WANNUNIT,*) 'BEGIN_DATAGRID_3D_UNKNOWN'
    WRITE(WANNUNIT,*) NGX+1,NGY+1,NGZ+1
    WRITE(WANNUNIT,'(3F14.7)') [0.0_q,0.0_q,0.0_q]
    WRITE(WANNUNIT,'(3F14.7)') LAT(:,1)
    WRITE(WANNUNIT,'(3F14.7)') LAT(:,2)
    WRITE(WANNUNIT,'(3F14.7)') LAT(:,3)

    ! A bit of complicated logic because XCrysden expects what they call
    ! 'general grids' instead of 'periodic grids'. For more details see:
    !  http://www.xcrysden.org/doc/XSF.html (retrieved 19/02/2020)
    DO NZ=1,NGZ+1
       NZI = NZ; IF (NZ>NGZ) NZI=1
       DO NY=1,NGY+1
          NYI = NY; IF (NY>NGY) NYI=1
          DO NX=1,NGX+1
             NXI = NX; IF (NX>NGX) NXI=1
             IND=NXI+(NYI-1)*NGX+(NZI-1)*NGX*NGY
             WRITE(WANNUNIT, '(E13.5)',ADVANCE='NO') REAL(DAT(IND))
          ENDDO
          WRITE(WANNUNIT,*)
       ENDDO
       WRITE(WANNUNIT,*)
    ENDDO
    WRITE(WANNUNIT,*) 'END_DATAGRID_3D'
    WRITE(WANNUNIT,*) 'END_BLOCK_DATAGRID_3D'
    CLOSE(WANNUNIT)
 END SUBROUTINE WRITE_XSF

!
!> Write the WFs at the gamma point to a XSF file
!
 SUBROUTINE WRITE_WAVESPIN_XSF(PREFIX,W,LATT_CUR,P,NATOMS,ATOM_SYMBOLS,ATOM_CART_POS,XSF_IW)
    use wave_struct_def, ONLY: wavespin
    use poscar_struct_def, ONLY: latt
    use pseudo_struct_def, ONLY: potcar
    CHARACTER(LEN=*) :: PREFIX
    TYPE(wavespin) :: W
    TYPE(latt),INTENT(IN) :: LATT_CUR
    TYPE(potcar),INTENT(IN) :: P(:)
    INTEGER,INTENT(IN) :: NATOMS
    CHARACTER(LEN=*),INTENT(IN) :: ATOM_SYMBOLS(NATOMS)
    REAL(q),INTENT(IN) :: ATOM_CART_POS(3,NATOMS)
    INTEGER,ALLOCATABLE,INTENT(IN) :: XSF_IW(:)
    ! local
    INTEGER :: I, IB, ISP
    INTEGER :: NGX, NGY, NGZ
    INTEGER :: ISPINOR
    INTEGER :: NPWS
    LOGICAL :: LSKIP
    GDEF, ALLOCATABLE :: WF_BANDS(:,:)
    CHARACTER(LEN=50) :: FILENAME
    NGX=W%WDES%GRID%NGX
    NGY=W%WDES%GRID%NGY
    NGZ=W%WDES%GRID%NGZ
    ALLOCATE(WF_BANDS(NGX*NGY*NGZ*W%WDES%NRSPINORS,W%WDES%NB_TOT))
    DO ISP=1,W%WDES%ISPIN
       CALL GET_WAVE_FUNCTIONS(W,[0.0_q,0.0_q,0.0_q],ISP,P,LATT_CUR,-1,WF_BANDS)
       DO IB=1,W%WDES%NB_TOT
          ! Only write the orbitals that were selected
          IF (ALLOCATED(XSF_IW)) THEN
             LSKIP=.TRUE.
             DO I=1,SIZE(XSF_IW)
                IF (IB==XSF_IW(I)) LSKIP=.FALSE.
             ENDDO
             IF (LSKIP) CYCLE
          ENDIF
          DO ISPINOR=1,W%WDES%NRSPINORS
             NPWS = (NGX*NGY*NGZ)*(ISPINOR-1)
             FILENAME = XSF_FILENAME(PREFIX,IB,ISP,ISPINOR,W%WDES%NRSPINORS,W%WDES%ISPIN)
             CALL WRITE_XSF(FILENAME,LATT_CUR%A,NATOMS,ATOM_SYMBOLS,ATOM_CART_POS,NGX,NGY,NGZ,WF_BANDS(NPWS+1:,IB))
          ENDDO ! loop over spinors
       ENDDO ! loop over bands
    ENDDO ! loop over spin
 END SUBROUTINE WRITE_WAVESPIN_XSF


#if 0
!******************** SUBROUTINE MLWF_ROTATE_ORBITALS ******************
!
!***********************************************************************
      SUBROUTINE MLWF_ROTATE_ORBITALS( &
     &   WDES,W,KPOINTS,GRID,T_INFO,P,NONL_S,SYMM,LATT_CUR,IO)
      USE base
      USE lattice
      USE pseudo
      USE poscar
      USE mgrid
      USE msymmetry
      USE nonl_high
      USE wave_high
      USE mkpoints
      USE full_kpoints
      USE kpoints_change
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      TYPE(kpoints_struct) KPOINTS
      TYPE(grid_3d) GRID
      TYPE(type_info) T_INFO
      TYPE(potcar) P(T_INFO%NTYP)
      TYPE(nonl_struct) NONL_S
      TYPE(symmetry) SYMM
      TYPE(latt) LATT_CUR
      TYPE(in_struct) IO
      ! local variables

      IF (SYMM%ISYM>=0.AND. .NOT.WDES%LGAMMA) THEN
         ! switch of symmetry
         CALL NOSYMM(LATT_CUR%A,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS,SYMM%PTRANS, &
        &   SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%MAGROT,WDES%ISPIN,IO%IU6)
         ! reread k-points with LINVERSION=.FALSE. to generate full mesh
         CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR,.FALSE., &
        &   T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IO%IU6,IO%IU0)
         CALL KPAR_SYNC_ALL(WDES,W)
         CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_CUR,-1, IO%IU0)
         CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
      ENDIF

      CALL MLWF_ROTATE_ORBITALS_FULLK(WDES,W)

      RETURN
      END SUBROUTINE MLWF_ROTATE_ORBITALS


!******************** SUBROUTINE MLWF_ROTATE_ORBITALS_FULLK ************
!
!***********************************************************************
      SUBROUTINE MLWF_ROTATE_ORBITALS_FULLK(WDES,W)
      USE dfast
      USE wave_high
      USE full_kpoints
      USE kpoints_change
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      ! local variables
      TYPE(wavedes1) WDES1
      TYPE(wavefuna) WA
      INTEGER I,J,K,IP,IPP
      INTEGER ISP,NK,NKP,NW,NB,NBEXCL,NBwin
      COMPLEX(q) CTMP
      GDEF, ALLOCATABLE :: U(:,:)

      CALL CHECK_FULL_KPOINTS

      IF ((.NOT.ALLOCATED(U_matrix)).OR.(.NOT.ALLOCATED(U_matrix_opt))) THEN
         CALL vtutor%error("MLWF_ROTATE_ORBITALS_FULLK: ERROR: rotation matrices not available")
      ENDIF

      NW=SIZE(U_matrix,1)
      NB=SIZE(U_matrix_opt,1)
      IF (SIZE(lexclude_band)/=WDES%NB_TOT) THEN
         CALL vtutor%error("MLWF_ROTATE_ORBITALS_FULLK: ERROR: lexclude_band array is corrupt: " // &
            str(SIZE(lexclude_band)) // " " // str(WDES%NB_TOT))
      ENDIF
      NBEXCL=0
      DO I=1,SIZE(lexclude_band)
         IF (lexclude_band(I)) NBEXCL=NBEXCL+1
      ENDDO
      IF ((NB+NBEXCL)/=WDES%NB_TOT) THEN
         CALL vtutor%error("MLWF_ROTATE_ORBITALS_FULLK: ERROR: inconsistent number of bands: " // &
             str(NB+NBEXCL)//' '//str(WDES%NB_TOT))
      ENDIF

      ALLOCATE(U(WDES%NB_TOT,WDES%NB_TOT))

      spin: DO ISP=1,WDES%ISPIN
         kpoints: DO NK=1,KPOINTS_FULL%NKPTS
            ! find the corresponding entry in KPOINTS_FULL_ORIG, that contains
            ! the set of k-points used in the computation of the rotation matrices
            DO NKP=1,KPOINTS_FULL_ORIG%NKPTS
               IF (LIDENTICAL_KPOINT(KPOINTS_FULL%VKPT(:,NK),KPOINTS_FULL_ORIG%VKPT(:,NKP))) exit
            ENDDO
            IF (NKP>KPOINTS_FULL_ORIG%NKPTS) THEN
               CALL vtutor%error("MLWF_ROTATE_ORBITALS_FULLK: ERROR: no matching k-point found in &
                  &KPOINTS_FULL_ORIG " // str(NK))
            ENDIF
            ! setup the effective rotation matrix U at the present k-point
            U=0._q; IP=0; IPP=0
            bands: DO I=1,WDES%NB_TOT
               IF (lexclude_band(I)) CYCLE bands
               IP=IP+1
               IF (.NOT.lwindow(IP,NKP,ISP)) CYCLE bands
               IPP=IPP+1
               DO J=1,NW
                  CTMP=0._q
                  DO K=1,NW
                     CTMP=CTMP+U_matrix_opt(IPP,K,NKP,ISP)*U_matrix(K,J,NKP,ISP)
                  ENDDO
                  U(I,J)=CTMP
               ENDDO
            ENDDO bands

            CALL SETWDES(WDES,WDES1,NK)
            WA=ELEMENTS(W,WDES1,ISP)

            ! redistribute over plane wave coefficients
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF

            ! build linear combinations of the wave functions in
            ! accordance with the effective rotation matrix U
            CALL LINCOM('F',WA%CW_RED,WA%CPROJ_RED,U(1,1), &
           &     WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), &
           &     WDES1%NPL_RED,WDES1%NPRO_RED,WDES%NRPLWV_RED,WDES%NPROD_RED,WDES%NB_TOT, &
           &     WA%CW_RED,WA%CPROJ_RED)

            ! redistribute back over bands
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF
         ENDDO kpoints
      ENDDO spin

      DEALLOCATE(U)

      RETURN
      END SUBROUTINE MLWF_ROTATE_ORBITALS_FULLK

!******************** SUBROUTINE MLWF_ROTATE_ORBITALS_NOSYMM ***********
!
!***********************************************************************
      SUBROUTINE MLWF_ROTATE_ORBITALS_NOSYMM(WDES,W)
      USE dfast
      USE wave_high
      USE full_kpoints
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      ! local variables
      TYPE(wavedes1) WDES1
      TYPE(wavefuna) WA
      INTEGER I,J,K,IP,IPP
      INTEGER ISP,NK,NKP,NW,NB,NBEXCL,NBwin
      COMPLEX(q) CTMP
      GDEF, ALLOCATABLE :: U(:,:)

      CALL CHECK_FULL_KPOINTS

      IF ((.NOT.ALLOCATED(U_matrix)).OR.(.NOT.ALLOCATED(U_matrix_opt))) THEN
         CALL vtutor%error("MLWF_ROTATE_ORBITALS_NOSYMM: ERROR: rotation matrices not available")
      ENDIF

      NW=SIZE(U_matrix,1)
      NB=SIZE(U_matrix_opt,1)
      IF (SIZE(lexclude_band)/=WDES%NB_TOT) THEN
         CALL vtutor%error("MLWF_ROTATE_ORBITALS_NOSYMM: ERROR: lexclude_band array is corrupt: " // &
            str(SIZE(lexclude_band)) // " " // str(WDES%NB_TOT))
      ENDIF
      NBEXCL=0
      DO I=1,SIZE(lexclude_band)
         IF (lexclude_band(I)) NBEXCL=NBEXCL+1
      ENDDO
      IF ((NB+NBEXCL)/=WDES%NB_TOT) THEN
         WRITE(*,*) 'MLWF_ROTATE_ORBITALS_NOSYMM: ERROR: inconsistent number of bands:', &
        &   NB+NBEXCL,WDES%NB_TOT
      ENDIF

      ALLOCATE(U(WDES%NB_TOT,WDES%NB_TOT))

!      WRITE(*,*) KPOINTS_FULL%NKPTS
!      WRITE(*,'(3F14.7)') KPOINTS_FULL%VKPT

      spin: DO ISP=1,WDES%ISPIN
         kpoints: DO NK=1,WDES%NKPTS
            ! find the corresponding entry in KPOINTS_FULL, that contains
            ! the set of k-points used in the computation of the rotation matrices
            DO NKP=1,KPOINTS_FULL%NKPTS
               IF (LIDENTICAL_KPOINT(WDES%VKPT(:,NK),KPOINTS_FULL%VKPT(:,NKP))) exit
            ENDDO
            IF (NKP>KPOINTS_FULL%NKPTS) THEN
               CALL vtutor%error("MLWF_ROTATE_ORBITALS: ERROR: no matching k-point found in &
                  &KPOINTS_FULL " // str(NK))
            ENDIF
            ! setup the effective rotation matrix U at the present k-point
            U=0._q; IP=0; IPP=0
            bands: DO I=1,WDES%NB_TOT
               IF (lexclude_band(I)) CYCLE bands
               IP=IP+1
               IF (.NOT.lwindow(IP,NKP,ISP)) CYCLE bands
               IPP=IPP+1
               DO J=1,NW
                  CTMP=0._q
                  DO K=1,NW
                     CTMP=CTMP+U_matrix_opt(IPP,K,NKP,ISP)*U_matrix(K,J,NKP,ISP)
                  ENDDO
                  U(I,J)=CTMP
               ENDDO
            ENDDO bands

            CALL SETWDES(WDES,WDES1,NK)
            WA=ELEMENTS(W,WDES1,ISP)

            ! redistribute over plane wave coefficients
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF

!           CALL ROTATE_CDER_NOINTER(WDES, U, SIZE(U,1), NK, ISP)

            ! build linear combinations of the wave functions in
            ! accordance with the effective rotation matrix U
            CALL LINCOM('F',WA%CW_RED,WA%CPROJ_RED,U(1,1), &
           &     WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), &
           &     WDES1%NPL_RED,WDES1%NPRO_RED,WDES%NRPLWV_RED,WDES%NPROD_RED,WDES%NB_TOT, &
           &     WA%CW_RED,WA%CPROJ_RED)

            ! redistribute back over bands
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF
         ENDDO kpoints
      ENDDO spin

      DEALLOCATE(U)

      RETURN
      END SUBROUTINE MLWF_ROTATE_ORBITALS_NOSYMM


      SUBROUTINE MLWF_ROTATE_ORBITALS_NOSYMM_U(WDES,W,U)
      USE wave_high
      USE full_kpoints
      USE dfast
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      GDEF :: U(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)
      ! local variables
      TYPE(wavedes1) WDES1
      TYPE(wavefuna) WA
      INTEGER I,J,K,IP,IPP
      INTEGER ISP,NK,NKP,NW,NB,NBEXCL,NBwin
      GDEF :: UU(WDES%NB_TOT,WDES%NB_TOT)

      spin: DO ISP=1,WDES%ISPIN
         kpoints: DO NK=1,WDES%NKPTS
            UU=U(:,:,NK,ISP)

            CALL SETWDES(WDES,WDES1,NK)
            WA=ELEMENTS(W,WDES1,ISP)

            ! redistribute over plane wave coefficients
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF

            !CALL ROTATE_CDER_NOINTER(WDES, UU, SIZE(UU,1), NK, ISP)

            ! build linear combinations of the wave functions in
            ! accordance with the effective rotation matrix U
            CALL LINCOM('F',WA%CW_RED,WA%CPROJ_RED,UU, &
           &     WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), &
           &     WDES1%NPL_RED,WDES1%NPRO_RED,WDES%NRPLWV_RED,WDES%NPROD_RED,WDES%NB_TOT, &
           &     WA%CW_RED,WA%CPROJ_RED)

            ! redistribute back over bands
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF
         ENDDO kpoints
      ENDDO spin

      RETURN
    END SUBROUTINE MLWF_ROTATE_ORBITALS_NOSYMM_U
#endif

#if 1
!******************** SUBROUTINE MLWF_ROTATE_CDER **********************
!> @brief Compute the optical matrix elements (CDER) and rotate to Wannier basis
!> This routine computes the optical matrix elements on the IBZ and rotates
!> them to the full Brillouin zone and then rotates them to the Wannier basis
!***********************************************************************
 SUBROUTINE MLWF_ROTATE_CDER(W,KPOINTS,KPOINTS_FBZ,U,GRID,GRIDC,GRIDUS,C_TO_US,IRDMAX,LMDIM,&
                             T_INFO,NONLR_S,NONL_S,P,SV,CQIJ,CDIJ,LATT_INI,LATT_CUR,NW,CHAM_WANN,SYMM,INFO,IO)
    USE constant
    USE base
    USE pseudo
    USE lattice
    USE mkpoints
    USE msymmetry
    USE nonl_high
    USE wave_high
    USE subrot_cluster
    USE kpoints_change
    USE mlr_optic
    USE mlrf_main
    USE choleski
    USE density_of_states, ONLY : DENINI
    USE poscar, ONLY : type_info
    TYPE (wavespin) W
    TYPE (latt) :: LATT_INI
    TYPE (latt) :: LATT_CUR
    TYPE (kpoints_struct) :: KPOINTS
    TYPE (kpoints_struct) :: KPOINTS_FBZ
    TYPE (potcar) :: P(:)
    TYPE (type_info) :: T_INFO
    TYPE (nonlr_struct)   :: NONLR_S
    TYPE (nonl_struct)    :: NONL_S
    TYPE (grid_3d) :: GRIDC
    TYPE (grid_3d) :: GRIDUS
    TYPE (grid_3d) :: GRID
    TYPE (transit) :: C_TO_US
    TYPE (info_struct) :: INFO
    TYPE (in_struct) :: IO
    TYPE (symmetry) :: SYMM
    INTEGER :: IRDMAX, LMDIM
    !> U_{nm}^k tranformation matrix between Bloch and Wannier functions
    GDEF, INTENT(IN) :: U(:,:,:,:)
    INTEGER,INTENT(IN) :: NW !< Number of wannier functions
    !> Optical matrix elements in the Wannier basis (NW,NW,KPOINTS_FBZ%NKPTS,W%WDES%ISPIN,3)
    GDEF, ALLOCATABLE, INTENT(OUT) :: CHAM_WANN(:,:,:,:,:)

    OVERLAP :: CQIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    OVERLAP :: CDIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    RGRID   :: SV(DIMREAL(GRID%MPLWV),W%WDES%NCDIJ)

    ! local variables
    INTEGER :: IK, IKF, IB, ISP, IDIR, NB
    TYPE(wavedes1), TARGET :: WDESK
    GDEFS,ALLOCATABLE :: CHAM_TMP(:,:,:,:,:) !< (W%WDES%NB_TOT,W%WDES%NB_TOT,KPOINTS_FBZ%NKPTS,W%WDES%ISPIN,3)
    GDEFS,ALLOCATABLE :: CHAM_FBZ(:,:,:) !< (W%WDES%NB_TOT,W%WDES%NB_TOT,3)
    GDEFS,ALLOCATABLE :: C(:,:) !< (NW,W%WDES%NB_TOT)
    GDEFS,ALLOCATABLE :: D(:,:) !< (W%WDES%NB_TOT,W%WDES%NB_TOT)
    GDEFS,ALLOCATABLE :: V(:,:) !< (W%WDES%NB_TOT,NW)
    GDEFS,ALLOCATABLE :: R(:,:) !< (NW,NW)
    TYPE (eigenf_cluster_pointer),POINTER :: DEG_CLUSTER(:,:)
    !> Wavefunction container for computation of optical matrix elements
    TYPE (wavespin) :: WDIFF
    LOGICAL :: LDONE
    INTEGER :: IB1,IB2
    REAL(q) :: E1,E2,EDIFF
    REAL(q) :: S(3,3)

    LOGICAL :: LINV
    INTEGER :: I,J,K,L
    INTEGER :: NC
    INTEGER :: IDIR_IBZ, IK_IBZ, ISP_IBZ

    CALL CHECK_FULL_KPOINTS

    ! Expand the WF to the full Brillouin zone
    !CALL RE_GEN_LAYOUT( GRID, W%WDES, KPOINTS_FBZ, LATT_CUR, LATT_INI,-1, IO%IU0)
    !CALL REALLOCATE_WAVE( W, GRID, W%WDES, NONL_S, T_INFO, P, LATT_CUR)

    ! Allocate wavefunction object to hold the derivative of the wavefunction
    NB = W%WDES%NB_TOT
    ALLOCATE(C(NB,NW))
    ALLOCATE(V(NB,NW))
    ALLOCATE(R(NW,NW))
    ALLOCATE(CHAM_TMP(NB,NB,KPOINTS_FBZ%NKPTS,W%WDES%ISPIN,3))
    ALLOCATE(CHAM_FBZ(NB,NB,3))
    ALLOCATE(CHAM_WANN(NW,NW,KPOINTS_FBZ%NKPTS,W%WDES%ISPIN,3))
    CALL ALLOCW(W%WDES,WDIFF)

    ! Find degenerate energies
    NULLIFY(DEG_CLUSTER)
    CALL FIND_DEG_CLUSTERS(W%WDES, W, DEG_CLUSTER)

    ! Compute optical matrix elements for each direction in the FBZ
    DO IDIR=1,3
       IF (IO%IU6>=0) WRITE(IO%IU0,*) 'Computing direction', IDIR
       WDIFF%CELTOT=0
       WDIFF%CW   =0
       WDIFF%CPROJ=0

       CALL FOCK_K_DER_ANALYT(KPOINTS, GRID, LATT_CUR, LATT_INI, &
            T_INFO,  NONLR_S, NONL_S, W, WDIFF, LMDIM, P, CQIJ, &
            SYMM, IDIR, LDONE, IO%IU0, IO%IU6)

       CALL LRF_RPHI0( P,NONLR_S,NONL_S,W,LATT_CUR, &
            T_INFO,INFO,IO,GRID,GRIDC,GRIDUS,C_TO_US,IRDMAX, &
            CDIJ,CQIJ,SV,LMDIM,DEG_CLUSTER, IDIR, WDIFF, LDONE, .TRUE.)

       ! Pull out the <u_n| r | u_m > for the calculation of the dielectric function
       CALL INPROD_W(WDIFF,W,CHAM_TMP(:,:,:,:,IDIR),INFO%LOVERL,IO%IU0)
    ENDDO

    ! Transform to Wannier basis using
    ! D_{MN} = \sum_{mn} (U_{mM})^T D_{mn} U_{nN}
    ! for each k-point, whith
    ! m,n -> runs over Bloch states
    ! M,N -> runs over Wannier states
    ! and ^T means complex transpose
    DO ISP=1,W%WDES%ISPIN
       DO IKF=1,KPOINTS_FBZ%NKPTS

          DO IK=1,KPOINTS_FULL%NKPTS
             IF (LIDENTICAL_KPOINT(KPOINTS_FBZ%VKPT(:,IKF),KPOINTS_FULL%VKPT(:,IK))) EXIT
          END DO

          ! Get indexes to perform transformation
          IK_IBZ=KPOINTS_FULL%NEQUIV(IK)
          ISP_IBZ=ISP
          IF (KPOINTS_FULL%SPINFLIP(IK)==1) ISP_IBZ=3-ISP
          LINV=KPOINTS_FULL%LINV(IK)

          ! determine transformation matrix in real space when going from IBZ to FBZ k-point
          S=0
          DO L=1,3
             DO K=1,3
                DO J=1,3
                   DO I=1,3
                      S(L,I)=S(L,I)+LATT_CUR%A(L,K)*KPOINTS_FULL%ISYMOP(J,K,IK)*LATT_CUR%B(I,J)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

          !apply rotation from IBZ to FBZ
          CHAM_FBZ = 0
          IF (LINV) THEN
             DO IDIR_IBZ=1,3
                DO IDIR=1,3
                   DO IB2=1,NB
                      DO IB1=1,NB
                         IF (IB1<=IB2) THEN
                            CHAM_FBZ(IB1,IB2,IDIR) = CHAM_FBZ(IB1,IB2,IDIR)-S(IDIR,IDIR_IBZ)*CHAM_TMP(IB2,IB1,IK_IBZ,ISP_IBZ,IDIR_IBZ)
                         ELSE
                            CHAM_FBZ(IB1,IB2,IDIR) = CHAM_FBZ(IB1,IB2,IDIR)-S(IDIR,IDIR_IBZ)*GCONJG(CHAM_TMP(IB1,IB2,IK_IBZ,ISP_IBZ,IDIR_IBZ))
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             DO IDIR_IBZ=1,3
                DO IDIR=1,3
                   CHAM_FBZ(:,:,IDIR) = CHAM_FBZ(:,:,IDIR)+S(IDIR,IDIR_IBZ)*CHAM_TMP(:,:,IK_IBZ,ISP_IBZ,IDIR_IBZ)
                ENDDO
             ENDDO
          ENDIF

          ! apply rotation to Wannier states
          DO IDIR=1,3
             V = U(:,:,IKF,ISP)
#ifdef gammareal
             CALL SGEMM('N', 'N', NB, NW, NB, 1.0, CHAM_FBZ(:,:,IDIR), NB, V, NW, 0.0, C, NB)
             CALL SGEMM('T', 'N', NW, NW, NB, 1.0, V, NB, C, NB, 0.0, R, NW)
#else
             CALL CGEMM('N', 'N', NB, NW, NB, CMPLX(1.0,0.0), CHAM_FBZ(:,:,IDIR), NB, V, NB, (0.0,0.0), C, NB)
             CALL CGEMM('C', 'N', NW, NW, NB, CMPLX(1.0,0.0), V, NB, C, NB, (0.0,0.0), R, NW)
#endif
             CHAM_WANN(:,:,IKF,ISP,IDIR) = R
          END DO !IDIR
       ENDDO !KPOINTS
    ENDDO !SPIN

    CALL DEALLOCW(WDIFF)
    CALL FREE_DEG_CLUSTERS(W%WDES, DEG_CLUSTER)

 END SUBROUTINE MLWF_ROTATE_CDER
#else
!******************** SUBROUTINE MLWF_ROTATE_CDER **********************
!> @brief Compute the optical matrix elements (CDER) and rotate to Wannier basis
!>! This routine expands the WFs to the full Brillouin zone then
!>! computes the optical matrix elements and rotates them to the Wannier basis
!***********************************************************************
 SUBROUTINE MLWF_ROTATE_CDER(W,KPOINTS,KPOINTS_FBZ,U,GRID,GRIDC,GRIDUS,C_TO_US,IRDMAX,LMDIM,&
                             T_INFO,NONLR_S,NONL_S,P,SV,CQIJ,CDIJ,LATT_INI,LATT_CUR,NW,CHAM_WANN,SYMM,INFO,IO)
    USE constant
    USE base
    USE pseudo
    USE lattice
    USE mkpoints
    USE msymmetry
    USE nonl_high
    USE wave_high
    USE subrot_cluster
    USE kpoints_change
    USE mlr_optic
    USE mlrf_main
    USE choleski
    USE density_of_states, ONLY : DENINI
    USE poscar, ONLY : type_info
    TYPE (wavespin) W
    TYPE (latt) :: LATT_INI
    TYPE (latt) :: LATT_CUR
    TYPE (kpoints_struct) :: KPOINTS
    TYPE (kpoints_struct) :: KPOINTS_FBZ
    TYPE (potcar) :: P(:)
    TYPE (type_info) :: T_INFO
    TYPE (nonlr_struct)   :: NONLR_S
    TYPE (nonl_struct)    :: NONL_S
    TYPE (grid_3d) :: GRIDC
    TYPE (grid_3d) :: GRIDUS
    TYPE (grid_3d) :: GRID
    TYPE (transit) :: C_TO_US
    TYPE (info_struct) :: INFO
    TYPE (in_struct) :: IO
    TYPE (symmetry) :: SYMM
    INTEGER :: IRDMAX, LMDIM
    !> U_{nm}^k tranformation matrix between Bloch and Wannier functions
    GDEF, INTENT(IN) :: U(:,:,:,:)
    INTEGER,INTENT(IN) :: NW !< Number of wannier functions
    !> Optical matrix elements in the Wannier basis (NW,NW,KPOINTS_FBZ%NKPTS,W%WDES%ISPIN,3)
    GDEF, ALLOCATABLE, INTENT(OUT) :: CHAM_WANN(:,:,:,:,:)

    OVERLAP :: CQIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    OVERLAP :: CDIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    RGRID   :: SV(DIMREAL(GRID%MPLWV),W%WDES%NCDIJ)

    ! local variables
    INTEGER :: IK, IB, ISP, IDIR, NB
    TYPE(wavedes1), TARGET :: WDESK
    GDEFS,ALLOCATABLE :: CHAM_TMP(:,:,:,:) !< (W%WDES%NB_TOT,W%WDES%NB_TOT,KPOINTS_FBZ%NKPTS,W%WDES%ISPIN)
    GDEFS,ALLOCATABLE :: C(:,:) !< (NW,W%WDES%NB_TOT)
    GDEFS,ALLOCATABLE :: D(:,:) !< (W%WDES%NB_TOT,W%WDES%NB_TOT)
    GDEFS,ALLOCATABLE :: V(:,:) !< (W%WDES%NB_TOT,NW)
    GDEFS,ALLOCATABLE :: R(:,:) !< (NW,NW)
    TYPE (eigenf_cluster_pointer),POINTER :: DEG_CLUSTER(:,:)
    !> Wavefunction container for computation of optical matrix elements
    TYPE (wavespin) :: WDIFF
    LOGICAL :: LDONE
    INTEGER :: IB1,IB2
    REAL(q) :: E1,E2,EDIFF

    CALL CHECK_FULL_KPOINTS

    ! Expand the WF to the full Brillouin zone
    CALL RE_GEN_LAYOUT( GRID, W%WDES, KPOINTS_FBZ, LATT_CUR, LATT_INI,-1, IO%IU0)
    CALL REALLOCATE_WAVE( W, GRID, W%WDES, NONL_S, T_INFO, P, LATT_CUR)

    ! Allocate wavefunction object to hold the derivative of the wavefunction
    NB = W%WDES%NB_TOT
    ALLOCATE(C(NB,NW))
    ALLOCATE(V(NB,NW))
    ALLOCATE(R(NW,NW))
    ALLOCATE(CHAM_TMP(NB,NB,KPOINTS_FBZ%NKPTS,W%WDES%ISPIN))
    ALLOCATE(CHAM_WANN(NW,NW,KPOINTS_FBZ%NKPTS,W%WDES%ISPIN,3))
    CALL ALLOCW(W%WDES,WDIFF)

    ! Find degenerate energies
    NULLIFY(DEG_CLUSTER)
    CALL FIND_DEG_CLUSTERS(W%WDES, W, DEG_CLUSTER)

    ! Compute optical matrix elements for each direction in the FBZ
    DO IDIR=1,3
       IF (IO%IU6>=0) WRITE(IO%IU0,*) 'Computing direction', IDIR
       WDIFF%CELTOT=0
       WDIFF%CW   =0
       WDIFF%CPROJ=0

       CALL FOCK_K_DER_ANALYT(KPOINTS, GRID, LATT_CUR, LATT_INI, &
            T_INFO,  NONLR_S, NONL_S, W, WDIFF, LMDIM, P, CQIJ, &
            SYMM, IDIR, LDONE, IO%IU0, IO%IU6)

       CALL LRF_RPHI0( P,NONLR_S,NONL_S,W,LATT_CUR, &
            T_INFO,INFO,IO,GRID,GRIDC,GRIDUS,C_TO_US,IRDMAX, &
            CDIJ,CQIJ,SV,LMDIM,DEG_CLUSTER, IDIR, WDIFF, LDONE, .TRUE.)

       ! Pull out the <u_n| r | u_m > for the calculation of the dielectric function
       CALL INPROD_W(WDIFF,W,CHAM_TMP(:,:,:,:),INFO%LOVERL,IO%IU0)

       ! Transform to Wannier basis using
       ! D_{MN} = \sum_{mn} (U_{mM})^T D_{mn} U_{nN}
       ! for each k-point, whith
       ! m,n -> runs over Bloch states
       ! M,N -> runs over Wannier states
       ! and ^T means complex transpose
       DO ISP=1,W%WDES%ISPIN
          DO IK=1,KPOINTS_FBZ%NKPTS
             V = U(:,:,IK,ISP)
#ifdef gammareal
             CALL SGEMM('N', 'N', NB, NB, NB, 1.0, CHAM_TMP(:,:,IK,ISP), W%WDES%NB, U(:,:,IK,ISP), NW, 0.0, C, NB)
             CALL SGEMM('T', 'N', NW, NB, NW, 1.0, U(:,:,IK,ISP), NW, C, NW, 0.0, R, NW)
             CHAM_WANN(:,:,IK,ISP,IDIR) = R
#else
             CALL CGEMM('N', 'N', NB, NW, NB, CMPLX(1.0,0.0), CHAM_TMP(:,:,IK,ISP), NB, V, NB, (0.0,0.0), C, NB)
             CALL CGEMM('C', 'N', NW, NW, NB, CMPLX(1.0,0.0), V, NB, C, NB, (0.0,0.0), R, NW)
             CHAM_WANN(:,:,IK,ISP,IDIR) = R
#endif
          ENDDO !KPOINTS
       ENDDO !SPIN
    ENDDO !IDIR

    CALL DEALLOCW(WDIFF)
    CALL FREE_DEG_CLUSTERS(W%WDES, DEG_CLUSTER)

    ! Return wavefunction to previous state
    CALL RE_GEN_LAYOUT( GRID, W%WDES, KPOINTS, LATT_CUR, LATT_CUR,-1, IO%IU0)
    CALL REALLOCATE_WAVE( W, GRID, W%WDES, NONL_S, T_INFO, P, LATT_CUR)

 END SUBROUTINE MLWF_ROTATE_CDER
#endif

#if 0
!******************** SUBROUTINE MLWF_DPSI_DK **************************
!> @brief Calculate the derivatives of the Wannier functions wrt. k
!***********************************************************************
      SUBROUTINE MLWF_DPSI_DK(W,KPOINTS,LATT_CUR,RPHI)
      USE constant
      USE pseudo
      USE lattice
      USE mkpoints
      USE full_kpoints
      USE wave_high
      USE twoelectron4o
      USE pead, ONLY : FAC
      TYPE(wavespin) W
      TYPE(latt) LATT_CUR
      TYPE(kpoints_struct) KPOINTS
      COMPLEX(qs) RPHI(W%WDES%NRPLWV,W%WDES%NBANDS,KPOINTS%NKPTS,W%WDES%ISPIN,3)
      ! local variables
      TYPE (wavedes1) WDESK,WDESIK
      TYPE (wavefun1) WK
      INTEGER ISP,IK,NK,IB,ISPINOR,IBG
      INTEGER IDIR,J,IDELTA,IORDER,ISGN
      REAL(q) K0(3),K(3),DK(3),KP(3)
      COMPLEX(q), ALLOCATABLE :: CPHASE(:)
      LOGICAL LPHASE

      IORDER=4

      DK(1)=1._q/REAL(KPOINTS%NKPX,KIND=q)
      DK(2)=1._q/REAL(KPOINTS%NKPY,KIND=q)
      DK(3)=1._q/REAL(KPOINTS%NKPZ,KIND=q)

      RPHI=0
      spin: DO ISP=1,W%WDES%ISPIN
      kpoint: DO IK=1,KPOINTS%NKPTS
         CALL SETWDES(W%WDES,WDESIK,IK)
         K0(:)=KPOINTS%VKPT(:,IK)

         dir: DO IDIR=1,3
         delta: DO IDELTA=1,IORDER
         sgn: DO ISGN=-1,1,2

            K(:)=K0(:)
            K(IDIR)=K0(IDIR)+ISGN*IDELTA*DK(IDIR)
            NK=KPOINT_IN_FULL_GRID(K,KPOINTS_FULL)
            KP(:)=KPOINTS_FULL%VKPT(:,NK)

            CALL SETWDES(W%WDES,WDESK,NK)
            CALL NEWWAV(WK,WDESK,.TRUE.)

            ALLOCATE(CPHASE(WDESK%GRID%MPLWV))
            CALL SETPHASE( KP(:)-K(:),WDESK%GRID,CPHASE,LPHASE)

!           bands: DO IB=1,MLWF_num_wann
            bands: DO IB=1,W%WDES%NBANDS
               ! global band index
#ifdef MPI
               IBG=(IB-1)*W%WDES%NB_PAR+W%WDES%COMM_INTER%NODE_ME
               IF(IBG>MLWF_num_wann) EXIT
#endif

               CALL W1_COPY(ELEMENT(W,WDESK,IB,ISP),WK)
               CALL FFTWAV_W1(WK)

               IF (LPHASE) CALL APPLY_PHASE(WDESK%GRID,CPHASE,WK%CR,WK%CR)

               WK%CR=WK%CR/WDESK%GRID%NPLWV
               DO ISPINOR=0,WDESK%NRSPINORS-1
                  CALL FFTEXT(WDESIK%NGVECTOR,WDESIK%NINDPW(1),WK%CR(1+ISPINOR*WDESK%GRID%MPLWV),WK%CW(1+ISPINOR*WDESIK%NGVECTOR),WDESK%GRID,.FALSE.)
               ENDDO

               RPHI(1:SIZE(WK%CW),IB,IK,ISP,IDIR)=RPHI(1:SIZE(WK%CW),IB,IK,ISP,IDIR)+CMPLX(ISGN*WK%CW(1:SIZE(WK%CW))*FAC(IDELTA,IORDER)/DK(IDIR)/TPI/2._q,KIND=qs)
!               DO J=1,3
!                  RPHI(1:SIZE(WK%CW),IB,IK,ISP,J)=RPHI(1:SIZE(WK%CW),IB,IK,ISP,J)+CMPLX(ISGN*LATT_CUR%A(J,IDIR)*WK%CW(1:SIZE(WK%CW))*FAC(IDELTA,IORDER)/DK(IDIR)/TPI/2._q,KIND=qs)
!               ENDDO
            ENDDO bands

            DEALLOCATE(CPHASE)
            CALL DELWAV(WK,.TRUE.)

         ENDDO sgn
         ENDDO delta
         ENDDO dir

      ENDDO kpoint
      ENDDO spin

      RETURN
      END SUBROUTINE MLWF_DPSI_DK
#endif

!******************** SUBROUTINE GET_WAVE_FUNCTIONS ******************
!> @brief get the wavefunctions in real space in a format compatible with the
!  UNK files needed by wannier90 to plot the Wannier functions.
!  Can write directly to file or return the WFs as an WVFN_BANDS(NGX*NGY*NGZ,NB_TOT) array
!***********************************************************************
 SUBROUTINE GET_WAVE_FUNCTIONS(W,K,ISP,P,LATT_CUR,IU,WVFN_BANDS,CELEN,EXCLUDE_BAND)
    USE pseudo
    USE lattice
    USE wave_high
    USE full_kpoints
    USE pead, ONLY : LUSEPEAD,WDES_FULL_PEAD,PEAD_WA_ROTATE
    TYPE(wavespin),INTENT(IN) :: W
    TYPE(potcar),INTENT(IN) :: P(:) !< Potcar datastructure
    TYPE(latt),INTENT(IN) :: LATT_CUR !< Current lattice datastructure
    INTEGER,INTENT(IN) :: ISP !< Spin index for which to obtain the real-space WF
    !> Unit of the file to write the real-space WF. Negative number to skip writting
    INTEGER,INTENT(IN) :: IU
    REAL(q),INTENT(IN) :: K(3)  !< Coordinates of the k-point for which to write the real-space WF
    LOGICAL,OPTIONAL :: EXCLUDE_BAND(W%WDES%NB_TOT)
    !> Optionally return the orbitals for this k-point in a single array
    GDEF,OPTIONAL,INTENT(OUT) :: WVFN_BANDS(:,:)
    !> Optionally return the eigenvalues of this array
    COMPLEX(q),OPTIONAL,INTENT(OUT) :: CELEN(:)
    ! local variables
    TYPE(wavespin) WP
    TYPE(wavefuna) WK
    TYPE(wavedes), POINTER :: WDES
    TYPE(wavedes1), TARGET :: WDESK
    TYPE (wavefun1), ALLOCATABLE :: WSTRIP(:),WCOLLECT(:)

    TYPE(rotation_handle), POINTER :: ROT_HANDLE

    INTEGER NK,NB,ISTRIP,ISPINOR,NPW,NPWS
    INTEGER NX,NY,NZ,NC,NGX,NGY,NGZ,IND,I,NWRITTEN
    INTEGER, PARAMETER :: NSTRIP=1

    GDEF, ALLOCATABLE :: WVFN(:)
    REAL(q) :: VKPT(3), GRID_POS(3)

    WP=W
    ALLOCATE(WDES)
    WDES=W%WDES; IF (LUSEPEAD()) WDES = WDES_FULL_PEAD
    WDES%NB_TOT = W%WDES%NB_TOT
    WDES%NBANDS = W%WDES%NBANDS
    WP%WDES=>WDES
    NGX=WP%WDES%GRID%NGX; NGY=WP%WDES%GRID%NGY; NGZ=WP%WDES%GRID%NGZ

    CALL CHECK_FULL_KPOINTS

    NULLIFY(ROT_HANDLE)

    ! search for kpoint k in BZ
    NK=KPOINT_IN_FULL_GRID(K,KPOINTS_FULL)
    CALL SETWDES(WP%WDES,WDESK,NK)
    IF (NK==KPOINTS_FULL%NEQUIV(NK)) THEN
       ! k is a kpoint in the IBZ
       WK=ELEMENTS(WP,WDESK,ISP)
    ELSE
       ! k is not a kpoint in the IBZ
       CALL NEWWAVA(WK,WDESK,WDESK%NBANDS)
       CALL PEAD_WA_ROTATE(WP,P,LATT_CUR,ISP,WK)
    ENDIF

    IF (PRESENT(CELEN)) CELEN = W%CELTOT(:,KPOINTS_FULL%NEQUIV(NK),ISP)
    ALLOCATE(WSTRIP(NSTRIP))

    ALLOCATE(WCOLLECT(NSTRIP*WP%WDES%NB_PAR))
    DO NB=1,NSTRIP*WP%WDES%NB_PAR
       CALL NEWWAV(WCOLLECT(NB),WDESK,.TRUE.)
    ENDDO

    NPW = NGX*NGY*NGZ
    ALLOCATE(WVFN(NPW*W%WDES%NRSPINORS))

    NWRITTEN=0

    DO NB=1,WDESK%NBANDS,NSTRIP
       DO ISTRIP=NB,MIN(WDESK%NBANDS,NB+NSTRIP-1)
          WSTRIP(ISTRIP-NB+1)=ELEMENT(WK,ISTRIP)
       ENDDO
       CALL W1_GATHER_W1(WP,MIN(WDESK%NBANDS-NB+1,NSTRIP),WSTRIP,WCOLLECT)
       DO ISTRIP=1,MIN(WDESK%NBANDS-NB+1,NSTRIP)*WP%WDES%NB_PAR
          ! skip this band?
          IF (PRESENT(EXCLUDE_BAND)) THEN
            IF (EXCLUDE_BAND((NB-1)*WP%WDES%NB_PAR+ISTRIP)) CYCLE
          ENDIF
          I=0; WVFN=0
          ! might seem strange to duplicate these loops,
          ! but I don't want this conditional statement nested inside
          IF (WP%WDES%GRID%RL%NFAST==1) THEN
             DO ISPINOR=1,WP%WDES%NRSPINORS
             NPWS = NPW*(ISPINOR-1)
             DO NC=1,WP%WDES%GRID%RL%NCOL
                NY=WP%WDES%GRID%RL%I2(NC)
                NZ=WP%WDES%GRID%RL%I3(NC)
                DO NX=1,WP%WDES%GRID%RL%NROW
                   IND=NX+(NY-1)*NGX+(NZ-1)*NGX*NGY
                   WVFN(IND+NPWS)=WCOLLECT(ISTRIP)%CR(NX+(NC-1)*WP%WDES%GRID%RL%NROW+NPWS)
                ENDDO
             ENDDO
             ENDDO
          ELSEIF (WP%WDES%GRID%RL%NFAST==3) THEN
             DO ISPINOR=1,WP%WDES%NRSPINORS
             NPWS = NPW*(ISPINOR-1)
             DO NC=1,WP%WDES%GRID%RL%NCOL
                NX=WP%WDES%GRID%RL%I2(NC)
                NY=WP%WDES%GRID%RL%I3(NC)
                DO NZ=1,WP%WDES%GRID%RL%NROW
                   IND=NX+(NY-1)*NGX+(NZ-1)*NGX*NGY
                   WVFN(IND+NPWS)=WCOLLECT(ISTRIP)%CR(NZ+(NC-1)*WP%WDES%GRID%RL%NROW+NPWS)
                ENDDO
             ENDDO
             ENDDO
          ELSE
             CALL vtutor%error("GET_WAVE_FUNCTIONS: ERROR: W1 grid not set")
          ENDIF
          CALLMPI(M_sum_g(WP%WDES%COMM_INB,WVFN(1),NGX*NGY*NGZ*W%WDES%NRSPINORS))

          ! write WVFN to file
#ifdef MPI
          IF (WP%WDES%COMM%NODE_ME==WP%WDES%COMM%IONODE.AND.IU>=0) THEN
#endif
#ifdef gammareal
             WRITE(IU) (CMPLX(WVFN(I),0.0_q),I=1,NGX*NGY*NGZ)
#else
             WRITE(IU) (WVFN(I),I=1,NGX*NGY*NGZ)
             IF (W%WDES%NRSPINORS==2) WRITE(IU) (WVFN(I+NPW),I=1,NGX*NGY*NGZ)
#endif
#ifdef MPI
          ENDIF
#endif

          ! return bands on a big array
          NWRITTEN=NWRITTEN+1
          IF (PRESENT(WVFN_BANDS)) WVFN_BANDS(:,NWRITTEN) = WVFN
       ENDDO
    ENDDO

    ! some deallocation to be done
    CALL DEALLOCATE_ROT_HANDLE(ROT_HANDLE)
    IF (NK/=KPOINTS_FULL%NEQUIV(NK)) CALL DELWAVA(WK)
    DO NB=1,NSTRIP*WP%WDES%NB_PAR
       CALL DELWAV(WCOLLECT(NB),.TRUE.)
    ENDDO
    DEALLOCATE(WSTRIP,WCOLLECT)
    DEALLOCATE(WVFN)
    DEALLOCATE(WDES)

    RETURN
 END SUBROUTINE GET_WAVE_FUNCTIONS

!******************** SUBROUTINE GET_WAVE_FUNCTIONS_PHASE **************
!> @brief get the phase of the wavefunction for this k-point
!***********************************************************************
 SUBROUTINE GET_WAVE_FUNCTIONS_PHASE(W,KPT,WVFN_PHASE,SHIFT)
    USE wave_high
    USE constant, ONLY: CITPI
    TYPE(wavespin),INTENT(IN) :: W
    REAL(q),INTENT(IN) :: KPT(3) !< Coordinates of the k-point
    !> Return the phase factor
    GDEF,INTENT(OUT) :: WVFN_PHASE(:)
    REAL(q),OPTIONAL :: SHIFT(3) !< Shift the real space points by this value
    ! local variables
    REAL(q) :: MY_SHIFT(3), GRID_POS(3)
    INTEGER :: NGX,NGY,NGZ
    INTEGER :: NX,NY,NZ,NC,IND

    NGX=W%WDES%GRID%NGX; NGY=W%WDES%GRID%NGY; NGZ=W%WDES%GRID%NGZ
    MY_SHIFT = [0.0_q, 0.0_q, 0.0_q]; IF (PRESENT(SHIFT)) MY_SHIFT = SHIFT

#ifdef gammareal
    WVFN_PHASE(:) = 0.0_q
#else
    IF (W%WDES%GRID%RL%NFAST==1) THEN
       DO NC=1,W%WDES%GRID%RL%NCOL
          NY=W%WDES%GRID%RL%I2(NC)
          NZ=W%WDES%GRID%RL%I3(NC)
          DO NX=1,W%WDES%GRID%RL%NROW
             IND=NX+(NY-1)*NGX+(NZ-1)*NGX*NGY
             GRID_POS(1) = REAL(NX-1)/NGX+MY_SHIFT(1)
             GRID_POS(2) = REAL(NY-1)/NGY+MY_SHIFT(2)
             GRID_POS(3) = REAL(NZ-1)/NGZ+MY_SHIFT(3)
             GRID_POS(1) = GRID_POS(1)-anint(GRID_POS(1))
             GRID_POS(2) = GRID_POS(2)-anint(GRID_POS(2))
             GRID_POS(3) = GRID_POS(3)-anint(GRID_POS(3))
             WVFN_PHASE(IND)=EXP(-CITPI * DOT_PRODUCT(GRID_POS, KPT))
          ENDDO
       ENDDO
    ELSEIF (W%WDES%GRID%RL%NFAST==3) THEN
       DO NC=1,W%WDES%GRID%RL%NCOL
          NX=W%WDES%GRID%RL%I2(NC)
          NY=W%WDES%GRID%RL%I3(NC)
          DO NZ=1,W%WDES%GRID%RL%NROW
             IND=NX+(NY-1)*NGX+(NZ-1)*NGX*NGY
             GRID_POS(1) = REAL(NX-1)/NGX+MY_SHIFT(1)
             GRID_POS(2) = REAL(NY-1)/NGY+MY_SHIFT(2)
             GRID_POS(3) = REAL(NZ-1)/NGZ+MY_SHIFT(3)
             GRID_POS(1) = GRID_POS(1)-anint(GRID_POS(1))
             GRID_POS(2) = GRID_POS(2)-anint(GRID_POS(2))
             GRID_POS(3) = GRID_POS(3)-anint(GRID_POS(3))
             WVFN_PHASE(IND)=EXP(-CITPI * DOT_PRODUCT(GRID_POS, KPT))
          ENDDO
       ENDDO
    ENDIF
#endif
 END SUBROUTINE GET_WAVE_FUNCTIONS_PHASE

!***********************************************************************
!> @brief Remove the kpoints, unit cell and atoms blocks as well as
!  spinors tag from the wannier90 input file
!***********************************************************************
 SUBROUTINE FILTER_WIN_FILE(MLWF,INPUT,OUTPUT)
    USE string, ONLY : lowercase
    TYPE(mlwf_type),INTENT(IN) :: MLWF
    CHARACTER(LEN=*),INTENT(IN)  :: INPUT
    CHARACTER(LEN=*),INTENT(OUT) :: OUTPUT
    ! Local variables
    INTEGER,PARAMETER :: MAX_LINE_LENGTH = 255
    CHARACTER(LEN=MAX_LINE_LENGTH) :: BUFLIN
    INTEGER :: I, NLINE, BUF_SIZE
    INTEGER :: INPOS, OLD_INPOS
    INTEGER :: OUTPOS, OLD_OUTPOS

    LOGICAL :: KPOINTS_BLOCK, CELL_BLOCK, ATOMS_BLOCK
    LOGICAL  :: NOCCUR
    EXTERNAL :: NOCCUR

    KPOINTS_BLOCK = .FALSE.
    CELL_BLOCK = .FALSE.
    ATOMS_BLOCK = .FALSE.
    BUF_SIZE = 100
    NLINE = 0
    INPOS = 1
    OUTPOS = 1
    OUTPUT = ""
    DO
      IF (INPOS == LEN(INPUT)) EXIT
      ! get the next line in the string
      OLD_INPOS = INPOS
      INPOS = INPOS + 1
      DO WHILE (INPOS<LEN(INPUT))
        IF (INPUT(INPOS:INPOS)==NEW_LINE('A')) EXIT
        INPOS = INPOS + 1
      END DO
      BUFLIN = INPUT(OLD_INPOS:INPOS-1) ! here I remove one line break
      ! filter unwanted information
      IF (NOCCUR(LOWERCASE(BUFLIN),'kpoints',-1)) THEN
        IF (NOCCUR(LOWERCASE(BUFLIN),'begin',-1)) KPOINTS_BLOCK = .true.
        IF (NOCCUR(LOWERCASE(BUFLIN),'end',-1)) KPOINTS_BLOCK = .false.
        CYCLE
      END IF
      IF (NOCCUR(LOWERCASE(BUFLIN),'unit_cell_cart',-1)) THEN
        IF (NOCCUR(LOWERCASE(BUFLIN),'begin',-1)) CELL_BLOCK = .true.
        IF (NOCCUR(LOWERCASE(BUFLIN),'end',-1)) CELL_BLOCK = .false.
        CYCLE
      END IF
      IF (NOCCUR(LOWERCASE(BUFLIN),'atoms_',-1)) THEN
        IF (NOCCUR(LOWERCASE(BUFLIN),'begin',-1)) ATOMS_BLOCK = .true.
        IF (NOCCUR(LOWERCASE(BUFLIN),'end',-1)) ATOMS_BLOCK = .false.
        CYCLE
      END IF
      IF (NOCCUR(LOWERCASE(BUFLIN),'mp_grid',-1)) CYCLE
      IF (NOCCUR(LOWERCASE(BUFLIN),'spinors',-1)) CYCLE
      IF (NOCCUR(LOWERCASE(BUFLIN),'num_bands',-1)) CYCLE
      IF (NOCCUR(LOWERCASE(BUFLIN),'num_wann',-1)) CYCLE
      IF (NOCCUR(LOWERCASE(BUFLIN),'dis_win_min',-1).AND.(MLWF%DIS_WIN_MIN/=HUGE(MLWF%DIS_WIN_MIN))) CYCLE
      IF (NOCCUR(LOWERCASE(BUFLIN),'dis_win_max',-1).AND.(MLWF%DIS_WIN_MAX/=HUGE(MLWF%DIS_WIN_MAX))) CYCLE
      IF (NOCCUR(LOWERCASE(BUFLIN),'dis_froz_min',-1).AND.(MLWF%DIS_FROZ_MIN/=HUGE(MLWF%DIS_FROZ_MIN))) CYCLE
      IF (NOCCUR(LOWERCASE(BUFLIN),'dis_froz_max',-1).AND.(MLWF%DIS_FROZ_MAX/=HUGE(MLWF%DIS_FROZ_MAX))) CYCLE
      IF (NOCCUR(LOWERCASE(BUFLIN),'use_bloch_phases',-1)) CYCLE
      IF (NOCCUR(LOWERCASE(BUFLIN),'auto_projections',-1)) CYCLE
      IF (NOCCUR(LOWERCASE(BUFLIN),'generated automatically by vasp',-1)) CYCLE
      !IF (TRIM(BUFLIN)=='') CYCLE
      IF (KPOINTS_BLOCK.OR.CELL_BLOCK.OR.ATOMS_BLOCK) CYCLE
      NLINE = NLINE + 1
      OLD_OUTPOS = OUTPOS
      OUTPOS = OUTPOS + LEN_TRIM(BUFLIN)
      OUTPUT(OLD_OUTPOS:OUTPOS) = TRIM(BUFLIN)
    ENDDO

 END SUBROUTINE FILTER_WIN_FILE

!******************** SUBROUTINE COPYFILE ******************************
!
!***********************************************************************
 SUBROUTINE COPYFILE(SRC,TRG)
    CHARACTER(*) :: SRC,TRG
    CHARACTER(LEN=1) :: CHARAC
    OPEN(UNIT=99,FILE=SRC,STATUS='OLD')
    OPEN(UNIT=98,FILE=TRG,STATUS='REPLACE')
 10 READ(99,'(A)',ADVANCE='No',EOR=20,END=30) CHARAC
    WRITE(98,'(A)',ADVANCE='No') CHARAC
    GOTO 10
 20 WRITE(98,*)
    GOTO 10
 30 CLOSE(98); CLOSE(99)
    RETURN
 END SUBROUTINE COPYFILE


#if 0
!************************* SETYLM_TST ********************************
!
!> @brief calculate the real spherical harmonics
! for a set of grid points up to LMAX
!
! YLM(:,1:LMAX) is set to Y_lm(hat x)
!
! hat x must be lying on the unit sphere
!
! written by Georg Kresse
!
!*********************************************************************

    SUBROUTINE SETYLM_TST(LYDIM,INDMAX,YLM,X,Y,Z)
      USE prec
      USE asa
      USE constant
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      INTEGER LYDIM           ! maximum L
      INTEGER INDMAX          ! number of points (X,Y,Z)
      REAL(q) YLM(:,:)        ! spherical harmonics
      REAL(q) X(:),Y(:),Z(:)  ! x,y and z coordinates

! local variables
      REAL(q) FAK
      INTEGER IND,LSET,LM,LP,LMINDX,ISTART,IEND,LNEW,L,M,MP,IC
!-----------------------------------------------------------------------
! runtime check of workspace
!-----------------------------------------------------------------------
      IF ( UBOUND(YLM,2) < (LYDIM+1)**2) THEN
         CALL vtutor%bug("internal ERROR: SETYLM, insufficient L workspace", __FILE__, __LINE__)
      ENDIF

      IF ( UBOUND(YLM,1) < INDMAX) THEN
         CALL vtutor%bug("internal ERROR: SETYLM, insufficient INDMAX workspace", __FILE__, __LINE__)
      ENDIF

      FAK=1/(2._q * SQRT(PI))
!-----------------------------------------------------------------------
! here is the code for L=0, hard coded
!-----------------------------------------------------------------------
      IF (LYDIM <0) GOTO 100
!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
        YLM(IND,1)=FAK
      ENDDO
!-----------------------------------------------------------------------
! here is the code for L=1, once again hard coded
!-----------------------------------------------------------------------
      IF (LYDIM <1) GOTO 100
!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
        YLM(IND,2)  = (FAK*SQRT(3._q))*Y(IND)
        YLM(IND,3)  = (FAK*SQRT(3._q))*Z(IND)
        YLM(IND,4)  = (FAK*SQRT(3._q))*X(IND)
      ENDDO
!-----------------------------------------------------------------------
! code for L=2,
!-----------------------------------------------------------------------
      IF (LYDIM <2) GOTO 100
!DIR$ IVDEP
!OCL NOVREC
      DO IND=1,INDMAX
        YLM(IND,5)= (FAK*SQRT(15._q))  *X(IND)*Y(IND)
        YLM(IND,6)= (FAK*SQRT(15._q))  *Y(IND)*Z(IND)
        YLM(IND,7)= (FAK*SQRT(5._q)/2._q)*(3*Z(IND)*Z(IND)-1)
        YLM(IND,8)= (FAK*SQRT(15._q))  *X(IND)*Z(IND)
        YLM(IND,9)= (FAK*SQRT(15._q)/2._q)*(X(IND)*X(IND)-Y(IND)*Y(IND))
      ENDDO
!-----------------------------------------------------------------------
! code for L=3,
!-----------------------------------------------------------------------
      IF (LYDIM <3) GOTO 100

      DO IND=1,INDMAX
        YLM(IND,10)=(FAK*SQRT(7._q)/2._q) *Z(IND)*(2*Z(IND)*Z(IND)-3*X(IND)*X(IND)-3*Y(IND)*Y(IND))
        YLM(IND,11)=(FAK*SQRT(35._q/2._q)/2._q) *Y(IND)*(3*X(IND)*X(IND)-Y(IND)*Y(IND))
        YLM(IND,12)=(FAK*SQRT(35._q/2._q)/2._q) *X(IND)*(X(IND)*X(IND)-3*Y(IND)*Y(IND))
        YLM(IND,13)=(FAK*SQRT(105._q)/2._q) *Z(IND)*(X(IND)*X(IND)-Y(IND)*Y(IND))
        YLM(IND,14)=(FAK*SQRT(105._q)) *X(IND)*Y(IND)*Z(IND)
        YLM(IND,15)=(FAK*SQRT(21._q/2._q)/2._q) *Y(IND)*(4*Z(IND)*Z(IND)-X(IND)*X(IND)-Y(IND)*Y(IND))
        YLM(IND,16)=(FAK*SQRT(21._q/2._q)/2._q) *X(IND)*(4*Z(IND)*Z(IND)-X(IND)*X(IND)-Y(IND)*Y(IND))
      ENDDO

 100  CONTINUE

    END SUBROUTINE SETYLM_TST

!******************** SUBROUTINE MLWF_DHAM_DK ******************
!> @brief Get dH/dk in wannier basis?
!***********************************************************************
      SUBROUTINE MLWF_DHAM_DK(W,KPOINTS,LATT_CUR,HAM,DHAM)
        USE constant
        USE pseudo
        USE lattice
        USE mkpoints
        USE full_kpoints
        USE wave_high
        USE pead, ONLY : FAC
        TYPE(wavespin) W
        TYPE(latt) LATT_CUR
        TYPE(kpoints_struct) KPOINTS
        !> derivative of the Hamiltonian in Wannier basis
        COMPLEX(q),INTENT(OUT) :: DHAM(W%WDES%NB_TOT,W%WDES%NB_TOT,KPOINTS_FULL%NKPTS,W%WDES%ISPIN,3)
        !> Hamiltonian in Wannier basis
        COMPLEX(q),INTENT(IN)  :: HAM(W%WDES%NB_TOT,W%WDES%NB_TOT,KPOINTS_FULL%NKPTS,W%WDES%ISPIN)
        ! local variables
        TYPE (wavedes1) WDESK
        TYPE (wavefun1) WK
        INTEGER ISP,IK,NK,IB,ISPINOR
        INTEGER IDIR,J,IDELTA,IORDER,ISGN
        REAL(q) K0(3),K(3),DK(3),KP(3)

        ! calculate derivative of Hamiltonian at each k-point
        IORDER=4

        DK(1)=1._q/REAL(KPOINTS_FULL%NKPX,KIND=q)
        DK(2)=1._q/REAL(KPOINTS_FULL%NKPY,KIND=q)
        DK(3)=1._q/REAL(KPOINTS_FULL%NKPZ,KIND=q)

        DHAM=0

        spin2: DO ISP=1,W%WDES%ISPIN
           kpoint2: DO IK=1,KPOINTS%NKPTS

              K0(:)=KPOINTS%VKPT(:,IK)

              dir: DO IDIR=1,3
                 delta: DO IDELTA=1,IORDER
                    sgn: DO ISGN=-1,1,2
                       K(:)=K0(:)
                       K(IDIR)=K0(IDIR)+ISGN*IDELTA*DK(IDIR)
                       NK=KPOINT_IN_FULL_GRID(K,KPOINTS_FULL)
                       KP(:)=KPOINTS_FULL%VKPT(:,NK)

                       ! non-cartesian coordinates !!!
                       DHAM(:,:,IK,ISP,IDIR)=DHAM(:,:,IK,ISP,IDIR)+CMPLX(ISGN*HAM(:,:,NK,ISP)*FAC(IDELTA,IORDER)/DK(IDIR)/TPI/2._q,KIND=q)
!!                       DO J=1,3
!!!                          DHAM(:,:,IK,ISP,J)=DHAM(:,:,IK,ISP,J)+CMPLX(ISGN*LATT_CUR%A(J,IDIR)*HAM(:,:,NK,ISP)*FAC(IDELTA,IORDER)/DK(IDIR)/TPI/2._q,KIND=q)
!!                          ! derivative w.r.t. to dimensionless coordinate
!!!                          DHAM(:,:,IK,ISP,J)=DHAM(:,:,IK,ISP,J)+CMPLX(ISGN*HAM(:,:,NK,ISP)*FAC(IDELTA,IORDER)/DK(IDIR)/TPI/2._q,KIND=q)
!!                       ENDDO

                    ENDDO sgn
                 ENDDO delta
              ENDDO dir

           ENDDO kpoint2
        ENDDO spin2

        RETURN
      END SUBROUTINE MLWF_DHAM_DK
#endif

!***********************************************************************
!> @brief Get U_(nm)^K rotation matrix computed after calling Wannier90
!***********************************************************************
 SUBROUTINE MLWF_GET_U(MLWF,NKPTS,VKPT,U)
    USE wave_high
    USE full_kpoints
    TYPE(MLWF_type) :: MLWF
    INTEGER :: NKPTS
    REAL(q) :: VKPT(3,NKPTS)
    !> U_(nm)^K rotation matrix (NB_TOT,NW,NKPTS,ISPIN)
    GDEF,ALLOCATABLE,INTENT(OUT) :: U(:,:,:,:)

    ! local variables
    INTEGER I,J,K,IP,IPP
    INTEGER ISP,NK,NKP,NBEXCL
    COMPLEX(q) CTMP

    IF ((.NOT.ALLOCATED(MLWF%U_matrix)).OR.(.NOT.ALLOCATED(MLWF%U_matrix_opt))) THEN
       CALL vtutor%error("MLWF_GET_U: ERROR: rotation matrices not available")
    ENDIF

    IF (SIZE(MLWF%EXCLUDE_BAND)/=MLWF%NB_TOT) THEN
       CALL vtutor%error("MLWF_GET_U: ERROR: exclude_band array is corrupt: " // &
          str(SIZE(MLWF%EXCLUDE_BAND)) // " " // str(MLWF%NB_TOT))
    ENDIF

    ! Allocate U matrix
    ALLOCATE(U(MLWF%NB_TOT,MLWF%NUM_WANN,MLWF%NKPTS,MLWF%ISPIN))

    spin: DO ISP=1,MLWF%ISPIN
       kpoint: DO NK=1,NKPTS
          ! find the corresponding entry in KPOINTS_FULL, that contains
          ! the set of k-points used in the computation of the rotation matrices
          DO NKP=1,MLWF%NKPTS
             IF (LIDENTICAL_KPOINT(VKPT(:,NK),MLWF%kpt_latt(:,NKP))) exit
          ENDDO

          IF (NKP>MLWF%NKPTS) THEN
             CALL vtutor%error("MLWF_GET_U: ERROR: no matching k-point found in &
                &KPOINTS_FULL " // str(NK) // " " // str(VKPT(:,NK)))
          ENDIF

          ! setup the effective rotation matrix U at the present k-point
          CALL MLWF_GET_U_IK_ISP(MLWF,NKP,ISP,U(:,:,NK,ISP))
       ENDDO kpoint
    ENDDO spin
 END SUBROUTINE MLWF_GET_U

!> @brief Setup the effective rotation matrix U at the present k-point and spin index
 SUBROUTINE MLWF_GET_U_IK_ISP(MLWF,NK,ISP,U)
     TYPE(MLWF_type) :: MLWF
     INTEGER :: NK, ISP
     INTEGER :: IP, IPP
     !> U_(nm) rotation matrix (MLWF%NB_TOT,MLWF%NUM_WANN)
     GDEF,INTENT(OUT) :: U(:,:)
     ! local variables
     INTEGER :: I,J,K
     COMPLEX(q) :: CTMP

     U=0._q; IP=0; IPP=0
     band: DO I=1,MLWF%NB_TOT
        IF (MLWF%EXCLUDE_BAND(I)) CYCLE band
        IP=IP+1
        IF (.NOT.MLWF%lwindow(IP,NK,ISP)) CYCLE band
        IPP=IPP+1
        DO J=1,MLWF%NUM_WANN
           CTMP=0._q
           DO K=1,MLWF%NUM_WANN
              CTMP=CTMP+MLWF%U_matrix_opt(IPP,K,NK,ISP)*MLWF%U_matrix(K,J,NK,ISP)
           ENDDO
           U(I,J)=CTMP
        ENDDO
     ENDDO band

 END SUBROUTINE MLWF_GET_U_IK_ISP

!***********************************************************************
!> @brief calculate Hamiltonian matrices in Wannier basis at each k-point
! setup matrix here
! The Hamiltonian in the Wannier band basis
!
! w .. Wannier bands
! b .. Bloch bands (cell-periodic part only)
! < w_m |H| w_n > =: H
! w_m = U_jm b_j
! eps = eps_j delta_ij
!
! therefore:
! H_mn = U*_jm U_in < b_j|H|b_i > = U*_jm U_in eps_j delta_ji = U*_jm U_jn eps_j
!***********************************************************************
 SUBROUTINE MLWF_ROTATE_HAM(WDES,W,NW,HAM,NEQUIV,U)
    USE wave_high
    USE full_kpoints
    TYPE(wavedes) WDES
    TYPE(wavespin) W
    INTEGER,INTENT(IN) :: NW !< Number of Wannier orbitals
    INTEGER,INTENT(IN) :: NEQUIV(WDES%NKPTS) !< Mapping from FBZ to IBZ
    !> Hamiltonian in the Wannier basis
    GDEF,INTENT(OUT) :: HAM(NW,NW,WDES%NKPTS,W%WDES%ISPIN)
    !> U_{nm}^k tranformation matrix between Bloch and Wannier functions
    GDEF,INTENT(IN) :: U(WDES%NB_TOT,NW,WDES%NKPTS,WDES%ISPIN)
    ! local variables
    INTEGER ISP,IK,B1,B2,IB,IKIBZ
    HAM=0
    spin1: DO ISP=1,WDES%ISPIN
       kpoint1: DO IK=1,WDES%NKPTS
          IKIBZ = NEQUIV(IK)

          DO B1=1,NW
             DO B2=1,NW
                DO IB=1,WDES%NB_TOT
                   HAM(B2,B1,IK,ISP)=HAM(B2,B1,IK,ISP)+GCONJG(U(IB,B2,IK,ISP))*W%CELTOT(IB,IKIBZ,ISP)*U(IB,B1,IK,ISP)
                END DO
             END DO
          END DO
       END DO kpoint1
    END DO spin1
    RETURN
 END SUBROUTINE MLWF_ROTATE_HAM

!**********************************************************************
!> @brief Read U_{nm}^k tranformation matrix between
!  Bloch and Wannier functions from WANPROJ
!**********************************************************************
 SUBROUTINE READ_WANPROJ_FILE(WANPROJFILE,WDES,IO,WANPROJ_EXISTS,NW,UWAN)
    USE base
    USE wave_high
    USE mkpoints
    USE constant
    USE full_kpoints
    USE fileio
    TYPE (wavedes)     WDES
    TYPE (in_struct)   IO
    CHARACTER(LEN=*),INTENT(IN) :: WANPROJFILE
    LOGICAL, INTENT(OUT) :: WANPROJ_EXISTS !< Whether a WANPROJ file was found
    INTEGER, INTENT(OUT) :: NW !< Number of wannier states read from file
    !> U_{nm}^k tranformation matrix between Bloch and Wannier functions
    GDEF, ALLOCATABLE, INTENT(OUT) :: UWAN(:,:,:,:)
    !local
    INTEGER            NKPTS,NB_TOT,ERROR
    INTEGER            ISPIN,IERR
    INTEGER            ISP,NKP,NK, I,J,I_,J_,K
    REAL(q)            VKPT(3)
    REAL(q)            RTMP1,RTMP2

    PROFILING_START( 'read_wanproj' )

    !if projection matrices are set up, return imediately
    IF ( ALLOCATED(UWAN)) THEN
       PROFILING_STOP( 'read_wanproj' )
       RETURN
    ENDIF

    INQUIRE(FILE=WANPROJFILE,EXIST=WANPROJ_EXISTS)

    NKP=0
    NW=0
    IERR=0
    !if file exist only root opens it and reads it
    IF ( WANPROJ_EXISTS ) THEN
       IF ( IO%IU0>=0 ) THEN
          WRITE(*,'(" Found ",A," file, reading...")')WANPROJFILE
          OPEN(UNIT=99,FILE=WANPROJFILE,STATUS='OLD')
          ! header is ignored by default
          READ(99,*,IOSTAT=ERROR)
          IF ( ERROR /=0 ) THEN
             call vtutor%error('READ_WANPROJ: could not read header')
          ENDIF
          !read header: number of k-points, bands and wannier states
          READ(99,*,IOSTAT=ERROR)ISPIN,NKPTS,NB_TOT,NW
          !maybe we need to kill the job here
          IF (NKPTS < WDES%NKPTS) THEN
             call vtutor%error('READ_WANPROJ: reports inconsistent NKPTS '//str(NKPTS)//' '//str(WDES%NKPTS))
          ENDIF
          IF (ISPIN < WDES%ISPIN) THEN
             call vtutor%error('READ_WANPROJ: reports inconsistent ISPIN '//str(ISPIN)//' '//str(WDES%ISPIN))
          ENDIF
          IF (NB_TOT <  WDES%NB_TOT) THEN
             call vtutor%error('READ_WANPROJ: reports inconsistent NBTOT '//str(NB_TOT)//' '//str(WDES%NB_TOT))
          ENDIF

          !write info on band reduction
          IF ( NB_TOT > WDES%NB_TOT ) &
             WRITE(*,'(" downsampling from ",I4," -> ",I4, " bands")')NB_TOT, WDES%NB_TOT
          !write info on k-points downsampling
          IF ( NKPTS > WDES%NKPTS ) &
             WRITE(*,'(" downsampling from ",I4," -> ",I4, " k-points")')NKPTS, WDES%NKPTS
          !write info on spin reduction
          IF ( ISPIN > WDES%ISPIN ) &
             WRITE(*,'(" downsampling from ",I4," -> ",I4, " spins")')ISPIN, WDES%ISPIN
          !skip the k-points list
          DO NK=1,NKPTS
             READ(99,*)I,VKPT(:)
          ENDDO
       ENDIF

       !communicate the number of Wannier states
       CALLMPI(M_sum_i(WDES%COMM,NW,1))

       ! Allocate UWAN with the correct dimensions
       ALLOCATE(UWAN(WDES%NB_TOT,NW,WDES%NKPTS,WDES%ISPIN))
       UWAN = 0

       IF ( IO%IU0>=0 ) THEN
          K=0
          scan: DO !scan the wann file
             K=K+1
             !read spin polarization and current k-point from file:
             READ(99,*,IOSTAT=ERROR)ISP,NB_TOT,VKPT(1:3)
             !look in the k-point list, if there is an entry for this point
             k_list: DO NK = 1, WDES%NKPTS
                !check if current k-point of list corresponds to k-point of file
                IF (.NOT.LIDENTICAL_KPOINT(WDES%VKPT(:,NK),VKPT(:))) CYCLE
                ! if an entry is found, store the corresponding matrix to UWAN
                NKP=NKP+1  !this is the total number of k-points found in list
                DO I =1, NB_TOT !bands
                   DO J =1,NW   !wannier states
                      READ(99,*,IOSTAT=ERROR)I_,J_,RTMP1,RTMP2
                      IF ( ERROR == 0 ) THEN
                         IF ( I_ > WDES%NB_TOT .OR. J_ > WDES%NB_TOT ) CYCLE
                         !store the matrix
#ifdef gammareal
                         UWAN(I_,J_,NK,ISP)=RTMP1
#else
                         UWAN(I_,J_,NK,ISP)=CMPLX(RTMP1,RTMP2,q)
#endif
                      ELSE
                         WRITE(IO%IU0,'(A,4I6)')' ERROR reading U_nk for',ISP,NK,I,J
                         EXIT scan
                      ENDIF
                   ENDDO
                ENDDO
                EXIT
             ENDDO k_list

             !if no point was found in k-points list, skip to next point
             IF ( NK > WDES%NKPTS )THEN
                DO I =1, NB_TOT
                   DO J =1,NW
                      READ(99,*,IOSTAT=ERROR)I_,J_,RTMP1,RTMP2
                      IF ( ERROR /= 0 ) THEN
                         WRITE(*,*)'ERROR (2)'
                         EXIT scan
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF
             !maybe we have enough points
             IF ( NKP==WDES%ISPIN*WDES%NKPTS)THEN
                WRITE(*,*)'Success.'
                EXIT scan
             ENDIF
             !also get out if all points in file have been done
             IF ( K==NKPTS*ISPIN) EXIT scan
          ENDDO  scan
          CLOSE(99)
       ENDIF

       !Terminate job eventually
       CALLMPI( M_bcast_i(WDES%COMM,NKP,1) )
       IF (NKP/=WDES%NKPTS*WDES%ISPIN) THEN
          CALL vtutor%error("ERROR: not all k-points could be found in " // str(NKP) // " " // &
             str(WDES%NKPTS))
       ENDIF

       !when everything went well, communicate the matrix to other nodes, and we are done
#ifdef gammareal
       CALLMPI(M_sum_d(WDES%COMM,UWAN,WDES%NKPTS*NW*WDES%NB_TOT*WDES%ISPIN))
#else
       CALLMPI(M_sum_g(WDES%COMM,UWAN,WDES%NKPTS*NW*WDES%NB_TOT*WDES%ISPIN))
#endif
    ENDIF
    PROFILING_STOP( 'read_wanproj' )
 END SUBROUTINE READ_WANPROJ_FILE

!******************** SUBROUTINE CALC_WANNIER_TRAFO_SCDM **********************
!
!> @brief Calculates the Wannier transformation matrices according to the SCDM method.
!
!******************************************************************************
 SUBROUTINE CALC_WANNIER_TRAFO_SCDM(W, NUM_WANN, T, P, LATT_CUR, NKPTS, VKPT, IO)
    USE base,               ONLY  : in_struct
    USE wave_struct_def,    ONLY  : wavespin
    USE pseudo_struct_def,  ONLY  : potcar
    USE LATTICE,            ONLY  : latt
    USE mathtools,          ONLY  : RRQR_DECOMPOSER, SV_DECOMPOSER
    USE tutor,              ONLY  : vtutor
    USE string,             ONLY  : str

    TYPE(wavespin),            INTENT(IN)  :: W
    INTEGER,                   INTENT(IN)  :: NUM_WANN
    COMPLEX(q), ALLOCATABLE,   INTENT(OUT) :: T(:, :, :, :)
    TYPE(potcar),              INTENT(IN)  :: P(:)
    TYPE(latt),                INTENT(IN)  :: LATT_CUR
    INTEGER,                   INTENT(IN)  :: NKPTS
    REAL(q),                   INTENT(IN)  :: VKPT(3, NKPTS)
    TYPE(in_struct),           INTENT(IN)  :: IO

    GDEF, ALLOCATABLE    :: PSI(:, :)
    GDEF, ALLOCATABLE    :: PSI_D_SCDM(:, :)
    GDEF, ALLOCATABLE    :: WVFN_PHASE(:)
    INTEGER              :: ISP
    INTEGER              :: K_IDX, K_IDX_IBZ
    INTEGER              :: IBZ, IB, IW
    INTEGER              :: GAMMA_IDX, J
    INTEGER, ALLOCATABLE :: PIVOT_UP(:),PIVOT_DW(:)
    REAL(q)              :: NORM
    REAL(q)              :: V(3)

    INTEGER              :: IDX_UP, IDX_DW, CURRENT_PIVOT, NP, IK
    REAL(q), ALLOCATABLE :: SIGMA(:)
    REAL(q), ALLOCATABLE :: CUTOFF_FACTOR(:)
    COMPLEX(q),ALLOCATABLE :: CELEN(:)
    TYPE(RRQR_DECOMPOSER)  :: RANK_REVEAL_QR
    TYPE(SV_DECOMPOSER)    :: SVD
    INTEGER                :: INFO

    IF (NUM_WANN > W%WDES%NB_TOT) CALL vtutor%error( &
       "CALC_WANNIER_TRAFO_SCDM: Number of Wannier functions must not be larger than number of bands.")

    IF (ALLOCATED(T)) DEALLOCATE(T)
    ALLOCATE(T(W%WDES%NB_TOT, NUM_WANN, NKPTS, W%WDES%ISPIN))
    T = 0

    ! Find Gamma point
    GAMMA_IDX = 0
    DO K_IDX = 1, W%WDES%NKPTS
       IF (W%WDES%AT_GAMMA(K_IDX)) THEN
          GAMMA_IDX = K_IDX
          EXIT
       ENDIF
    ENDDO

    ! If there is no Gamma, use the first available k-point
    IF (GAMMA_IDX == 0) THEN
       GAMMA_IDX = 1
       CALL vtutor%warning( &
          "The Gamma point was not found in the k-point grid. The SCDM method now uses the k-point " // &
          str(W%WDES%VKPT(:, GAMMA_IDX)) // " instead.")
    ENDIF

    NP = W%WDES%GRID%RL%NP
    ALLOCATE(PSI(NP*W%WDES%NRSPINORS,W%WDES%NB_TOT))
    ALLOCATE(CELEN(W%WDES%NB_TOT))
    ALLOCATE(WVFN_PHASE(NP*W%WDES%NRSPINORS))

    CALL RANK_REVEAL_QR%INIT(W%WDES%NB_TOT, NP * W%WDES%NRSPINORS, INFO)
    IF (INFO /= 0) CALL vtutor%error( &
      "CALC_WANNIER_TRAFO_SCDM: xGEQP3 query returned with error code " // str(INFO))

    CALL SVD%INIT('S',W%WDES%NB_TOT,NUM_WANN,INFO)
    IF (INFO /= 0) CALL vtutor%error( &
      "CALC_WANNIER_TRAFO_SCDM: xGESDD query returned with error code " // str(INFO))

    DO ISP = 1, W%WDES%ISPIN
       !  _______________________________
       ! |                               |
       ! |   Work on Gamma point first   |
       ! |_______________________________|

       ! Calculate real-spave wave function
       CALL GET_WAVE_FUNCTIONS(W,[0.0_q, 0.0_q, 0.0_q],ISP,P,LATT_CUR,-1,PSI,CELEN)
       DO IB=1,W%WDES%NB_TOT
          NORM = REAL(DOT_PRODUCT(PSI(:NP,IB), PSI(:NP,IB)))
          IF (W%WDES%NRSPINORS==2) NORM = NORM + REAL(DOT_PRODUCT(PSI(NP+1:,IB), PSI(NP+1:,IB)))
          PSI(:,IB) = PSI(:,IB) / SQRT(NORM)
       ENDDO

       ! Calculate SCDM matrix
       IF (ALLOCATED(PSI_D_SCDM)) DEALLOCATE(PSI_D_SCDM)
       ALLOCATE(PSI_D_SCDM(W%WDES%NB_TOT,NP*W%WDES%NRSPINORS))
       CUTOFF_FACTOR = CUTOFF_DES_GET_CUTOFF(CUTOFF_DES,ISP,CELEN)
       DO IB=1,W%WDES%NB_TOT
          PSI_D_SCDM(IB,:) = GCONJG(PSI(:,IB))*CUTOFF_FACTOR(IB)
       ENDDO

       ! Use rank-revealing QR decomposition to find a localized representation
       CALL RANK_REVEAL_QR%CALC(PSI_D_SCDM, INFO)
       IF (INFO /= 0) CALL vtutor%error( &
         "CALC_WANNIER_TRAFO_SCDM: xGEQP3 returned with error code " // str(INFO))

       ! In the case of non-collinear calculation we make sure the pivots correspond
       ! to the spin up and down components alternately (as defined in LPRJ_READER)
       !WRITE(*,*) PIVOT(:NUM_WANN) < NP
       IF (W%WDES%NRSPINORS==2) THEN
          ALLOCATE(PIVOT_UP(NP))
          ALLOCATE(PIVOT_DW(NP))
          IDX_UP = 1
          IDX_DW = 1
          ! Identify to which spinor component the pivots belong
          DO IW=1,W%WDES%GRID%RL%NP*2
             CURRENT_PIVOT = RANK_REVEAL_QR%PIVOT(IW)
             IF (CURRENT_PIVOT > W%WDES%GRID%RL%NP) THEN
                PIVOT_DW(IDX_DW) = CURRENT_PIVOT
                IDX_DW = IDX_DW + 1
             ELSE
                PIVOT_UP(IDX_UP) = CURRENT_PIVOT
                IDX_UP = IDX_UP + 1
             ENDIF
          ENDDO
          ! fill them in alternately in PIVOT (only up to NUM_WANN as only those are needed)
          DO IW=1,NUM_WANN/2
            RANK_REVEAL_QR%PIVOT((IW-1)*2+1) = PIVOT_UP(IW)
          ENDDO
          DO IW=1,NUM_WANN/2
            RANK_REVEAL_QR%PIVOT((IW-1)*2+2) = PIVOT_DW(IW)
          ENDDO
          DEALLOCATE(PIVOT_UP)
          DEALLOCATE(PIVOT_DW)
       ENDIF
       !WRITE(*,*) PIVOT(:NUM_WANN) < NP

       ! Calculate SCDM matrix
       DEALLOCATE(PSI_D_SCDM)
       ALLOCATE(PSI_D_SCDM(W%WDES%NB_TOT,NUM_WANN))
       CUTOFF_FACTOR = CUTOFF_DES_GET_CUTOFF(CUTOFF_DES,ISP,CELEN)
       DO IB=1,W%WDES%NB_TOT
          DO IW=1,NUM_WANN
             PSI_D_SCDM(IB,IW) = GCONJG(PSI(RANK_REVEAL_QR%PIVOT(IW),IB))*CUTOFF_FACTOR(IB)
          ENDDO
       ENDDO

       ! Orthonormalize via SVD
       CALL SVD%CALC(PSI_D_SCDM, INFO)
       IF (INFO /= 0) THEN
          CALL vtutor%error("CALC_WANNIER_TRAFO_PROJ: xGESDD returned with error code " // str(INFO))
       ENDIF
       ! Set transformation matrix (in Gamma-only version: real -> complex)
       T(:, :, GAMMA_IDX, ISP) = MATMUL(SVD%U_MAT, SVD%V_DAGGER_MAT)

       !  __________________________________
       ! |                                  |
       ! |   Loop over remaining k-points   |
       ! |__________________________________|
       DO K_IDX = 1, NKPTS

          ! Skip the gamma point as it was already computed
          IF (GAMMA_IDX==K_IDX) CYCLE

          ! Calculate real-spave wave function
          CALL GET_WAVE_FUNCTIONS(W,VKPT(:,K_IDX),ISP,P,LATT_CUR,-1,PSI,CELEN)
          DO IB=1,W%WDES%NB_TOT
             NORM = REAL(DOT_PRODUCT(PSI(:NP,IB), PSI(:NP,IB)))
             IF (W%WDES%NRSPINORS==2) NORM = NORM + REAL(DOT_PRODUCT(PSI(NP+1:,IB), PSI(NP+1:,IB)))
             PSI(:,IB) = PSI(:,IB) / SQRT(NORM)
          ENDDO

          CALL GET_WAVE_FUNCTIONS_PHASE(W,VKPT(:,K_IDX),WVFN_PHASE)

          ! Calculate SCDM matrix
          CUTOFF_FACTOR = CUTOFF_DES_GET_CUTOFF(CUTOFF_DES,ISP,CELEN)
          DO IB=1,W%WDES%NB_TOT
             DO IW=1,NUM_WANN
                J=RANK_REVEAL_QR%PIVOT(IW); IF ( J > NP) J = J-NP
                PSI_D_SCDM(IB,IW) = GCONJG(PSI(RANK_REVEAL_QR%PIVOT(IW),IB))*WVFN_PHASE(J)*CUTOFF_FACTOR(IB)
             ENDDO
          ENDDO

          ! Orthonormalize via SVD
          CALL SVD%CALC(PSI_D_SCDM, INFO)
          IF (INFO /= 0) THEN
             CALL vtutor%error( "CALC_WANNIER_TRAFO_PROJ: xGESDD returned with error code " // str(INFO))
          ENDIF
          ! Set transformation matrix (in Gamma-only version: real -> complex)
          T(:, :, K_IDX, ISP) = MATMUL(SVD%U_MAT, SVD%V_DAGGER_MAT)
       ENDDO
    ENDDO

 END SUBROUTINE CALC_WANNIER_TRAFO_SCDM

 !> @brief Get cutoff function for a set of energy states
 FUNCTION CUTOFF_DES_GET_CUTOFF(CUTOFF_DES,ISP,CELEN) RESULT(CUTOFF_FACTOR)
    TYPE(cutoff_type),INTENT(IN) :: CUTOFF_DES
    INTEGER,INTENT(IN) :: ISP
    COMPLEX(q),INTENT(IN) :: CELEN(:)
    REAL(q) :: CUTOFF_FACTOR(SIZE(CELEN))
    ! Local variables
    REAL(q) :: MU

    SELECT CASE(CUTOFF_DES%TYPE)
        CASE(CUTOFF_GAUSSIAN)
           CUTOFF_FACTOR = EXP(-(REAL(CELEN,q) - CUTOFF_DES%MU(ISP))**2 / CUTOFF_DES%SIGMA(ISP))
        CASE(CUTOFF_ERFC)
           CUTOFF_FACTOR = .5_q * ERFC((REAL(CELEN,q) - CUTOFF_DES%MU(ISP)) / CUTOFF_DES%SIGMA(ISP))
        CASE(CUTOFF_FERMI)
           CUTOFF_FACTOR = 1.0_q/(EXP((REAL(CELEN,q)-CUTOFF_DES%MU(ISP)))+1.0_q)
        CASE(CUTOFF_NUM_WANN)
           MU = REAL(CELEN(CUTOFF_DES%NUM_WANN),q)
           CUTOFF_FACTOR = .5_q * ERFC((REAL(CELEN,q) - MU) / CUTOFF_DES%SIGMA(ISP))
        CASE DEFAULT
           CUTOFF_FACTOR = 1._q
    END SELECT
 END FUNCTION CUTOFF_DES_GET_CUTOFF

 !> Compute transformation matrix to local orbitals using one-shot method
 !> starting from a non-unitary matrix of projections.
 !> This amounts to making the matrix "unitary" (it's still rectangular)
 !> by using a singular value decomposition
 ! _test_mani probably not yet implemented correctly for all k-points
 ! _test_mani SPINOR implementation missing
 SUBROUTINE CALC_WANNIER_TRAFO_PROJ(T, W, MLWF)
   USE wave_struct_def, ONLY: wavespin
   USE string,          ONLY: str
   USE tutor,           ONLY: vtutor
   USE full_kpoints,    ONLY: KPOINTS_FULL
   USE mathtools,       ONLY: SV_DECOMPOSER
   IMPLICIT none

   GDEF, ALLOCATABLE,   INTENT(OUT) :: T(:, :, :, :)
   TYPE(wavespin),      INTENT(IN)  :: W
   TYPE(MLWF_type),     INTENT(IN)  :: MLWF

   LOGICAL              :: DO_CUTOFF
   INTEGER              :: IK, ISP
   INTEGER              :: IKIBZ
   INTEGER              :: WAN_IDX
   TYPE(SV_DECOMPOSER)  :: SVD
   INTEGER              :: INFO
   REAL(q), ALLOCATABLE :: CUTOFF_FACTOR(:)

   !HM: The code still works correctly in this case, so I comment it
   !IF (MLWF%NB_TOT < MLWF%NUM_WANN) CALL vtutor%error( &
   !   "CALC_WANNIER_TRAFO_PROJ: Number of Wannier functions must not be larger than number of bands.")

   DO_CUTOFF = MLWF%NB_TOT > MLWF%NUM_WANN

   ALLOCATE(T(MLWF%NB_TOT, MLWF%NUM_WANN, MLWF%NKPTS, MLWF%ISPIN))
   T = MLWF%A_matrix

   CALL SVD%INIT("S", MLWF%NB_TOT, MLWF%NUM_WANN, INFO)
   IF (INFO /= 0) CALL vtutor%error( &
      "CALC_WANNIER_TRAFO_PROJ: xGESDD query returned with error code " // str(INFO))

   IF (DO_CUTOFF) ALLOCATE(CUTOFF_FACTOR(MLWF%NB_TOT))

   ! Loop over all spin and k components.
   DO ISP = 1, MLWF%ISPIN
      DO IK = 1, MLWF%NKPTS
         IF (DO_CUTOFF) THEN
            IKIBZ = KPOINTS_FULL%NEQUIV(IK)
            CUTOFF_FACTOR = CUTOFF_DES_GET_CUTOFF(CUTOFF_DES,ISP,W%CELTOT(:, IKIBZ, ISP))
            DO WAN_IDX = 1, MLWF%NUM_WANN
               T(:, WAN_IDX, IK, ISP) = T(:, WAN_IDX, IK, ISP) * CUTOFF_FACTOR
            ENDDO
         ENDIF

         CALL SVD%CALC(T(:, :, IK, ISP), INFO)
         IF (INFO /= 0) CALL vtutor%error( &
            "CALC_WANNIER_TRAFO_PROJ: xGESDD returned with error code " // str(INFO))

         T(:, :, IK, ISP) = MATMUL(SVD%U_MAT, SVD%V_DAGGER_MAT)
      ENDDO
   ENDDO
 ENDSUBROUTINE CALC_WANNIER_TRAFO_PROJ

END MODULE mlwf
