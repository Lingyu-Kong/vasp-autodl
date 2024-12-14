#include "symbol.inc"

MODULE elphon
   USE prec
   USE elphon_util
   USE base,               ONLY: in_struct
   USE tutor,              ONLY: vtutor
   USE string,             ONLY: str
   USE wave_struct_def,    ONLY: wavespin, &
                                 wavedes, &
                                 wavedes1
   USE msymmetry,          ONLY: PRIM_CELL_T
   USE poscar_struct_def,  ONLY: latt, &
                                 type_info
   USE wannier_mats,       ONLY: WANNIER_DESCRIPTOR_T, &
                                 WANNIER_MAT_SYMM_T, &
                                 WANNIER_MAT_MODES_T

! _test_mani UMKLAPP processes?
! _test_mani move phase calc to elphon_util?
! _test_mani check all USEs
! _test_mani check all unused
! _test_mani move hbar to top of module

   use lattice ! _test_mani for brillouin stuff.
   IMPLICIT NONE

   PRIVATE

   PUBLIC ELPHON_SETTINGS, ELPHON_INTERPOLATOR
   PUBLIC ELPH_READER, ELPH_CHECK_INPUT_CONSISTENCY
   PUBLIC ELPH_OVERLAP_WSWQ, ELPH_OUTPUT_STUFF, ELPH_IO

   ENUM, BIND(c)
      ENUMERATOR  :: METHOD_NONE
      ENUMERATOR  :: METHOD_LOCPROJ
      ENUMERATOR  :: METHOD_SCDM
      ENUMERATOR  :: METHOD_LOAD_PRIM
      ENUMERATOR  :: METHOD_LOAD_DATA
   END ENUM

   TYPE ELPHON_SETTINGS
      LOGICAL                    :: DO_WSWQ = .FALSE.
      LOGICAL                    :: WRITE_DATA = .FALSE.             !< Write elphon data to disk?
      LOGICAL                    :: WRITE_WANN_FUN = .FALSE.         !< Write Wannier functions as WAVECAR?
      LOGICAL                    :: ALL_ELECTRON = .FALSE.
      LOGICAL                    :: DO_POLAR = .FALSE.               !< Is code running in polar mode?
      LOGICAL                    :: DO_ZPR_DYNAMIC = .FALSE.         !< Include phonon frequency in ZPR?
      INTEGER                    :: BZ_SUBDIVISION = 1               !< How much the Brillouin Zone is subdivided during sampling.
      LOGICAL                    :: CALC_BANDSTRUCTURE = .FALSE.     !< Interpolate band structure using Wannier functions?
      LOGICAL                    :: CALC_ZPR = .FALSE.               !< Do we calculate the ZPR of the band gap?
      LOGICAL                    :: CALC_TEMPERATURE_GAP = .FALSE.   !< Do we calculate temperature-dependent band gap?
      REAL(q), ALLOCATABLE       :: BZ_POINTS(:, :)                  !< Calculation at these Brillouin-zone point (direct coordinates)
      REAL(q)                    :: POLAR_G_CUTOFF = 8._q            !< Cutoff in reciprocal space (multiple of longest G-vector)
      REAL(q)                    :: POLAR_EWALD_THRESHOLD = 1E-10_q  !< Threshold down to which Ewald summands contribute
      INTEGER(KIND(METHOD_NONE)) :: METHOD = METHOD_NONE             !< Method to construct Wannier orbitals
      INTEGER                    :: NUM_WANN_READ = 0                !< Number of Wannier orbitals to compute per primitive cell.
      INTEGER                    :: NUM_BANDS_MAX = 0                !< Maximum number of bands to include in the ZPR calculation.
      REAL(q)                    :: SCDM_MU = 0._q
      REAL(q)                    :: SCDM_SIGMA = 1._q
      LOGICAL                    :: EXTERNAL_INTERFACE = .FALSE.

      CONTAINS

      PROCEDURE,  PUBLIC   :: DO_INIT => SETTINGS_DO_INIT
      PROCEDURE,  PUBLIC   :: DO_PERTURB => SETTINGS_DO_PERTURB
      PROCEDURE,  PUBLIC   :: DO_CALC => SETTINGS_DO_CALC
      PROCEDURE,  PUBLIC   :: DO_LOAD_DATA => SETTINGS_DO_LOAD_DATA
   END TYPE

   ! Maybe useful for | Psi^(1) > calculation
   ! OVERLAP, ALLOCATABLE,         SAVE :: CQIJ_INIT(:, :, :, :)

   TYPE ELPHON_INTERPOLATOR
      LOGICAL, PRIVATE  :: IS_INITIALIZED = .FALSE.
      LOGICAL, PRIVATE  :: IS_COMPLETED = .FALSE.
      LOGICAL, PRIVATE  :: DO_POLAR = .FALSE.

      TYPE(wavespin),                           PRIVATE  :: W_WAN_INIT
      TYPE(WANNIER_DESCRIPTOR_T),               PRIVATE  :: WANNIER_DESCRIPTOR
      TYPE(WANNIER_MAT_SYMM_T), ALLOCATABLE,    PRIVATE  :: HAM_INIT(:)
      TYPE(WANNIER_MAT_MODES_T), ALLOCATABLE,   PRIVATE  :: COUPLING_WANN(:)
      TYPE(WANNIER_MAT_MODES_T), ALLOCATABLE,   PRIVATE  :: OVERLAP_PERT(:)

      OVERLAP, ALLOCATABLE,   PRIVATE  :: CDIJ_INIT(:, :, :, :)

      COMPLEX(q), ALLOCATABLE,   PRIVATE  :: COUPLING_WANN_LONG(:, :, :, :)
      REAL(q), ALLOCATABLE,      PRIVATE  :: EWALD_G_VECTORS(:, :)
      TYPE(POLAR_EWALD_T),       PRIVATE  :: POLAR_EWALD

      CONTAINS

      PROCEDURE,  PUBLIC   :: INIT => EPI_INIT
      FINAL                :: EPI_FINALIZE
      PROCEDURE,  PRIVATE  :: DEALLOC => EPI_DEALLOC
      PROCEDURE,  PUBLIC   :: CALC_COUPLING => EPI_CALC_COUPLING
      PROCEDURE,  PUBLIC   :: CALC_COUPLING_ATOMS => EPI_CALC_COUPLING_ATOMS
      PROCEDURE,  PUBLIC   :: PERTURB_WANNIER_PS => EPI_PERTURB_WANNIER_PS
      PROCEDURE,  PUBLIC   :: PERTURB_WANNIER_AE => EPI_PERTURB_WANNIER_AE

      PROCEDURE,  PUBLIC   :: WRITE_WANNIER_FUNCS => EPI_WRITE_WANNIER_FUNCS
      PROCEDURE,  PUBLIC   :: READ_WANNIER_FUNCS => EPI_READ_WANNIER_FUNCS
#ifdef VASP_HDF5
      PROCEDURE,  PUBLIC   :: WRITE_HDF5 => EPI_WRITE_HDF5
      PROCEDURE,  PUBLIC   :: READ_HDF5 => EPI_READ_HDF5
#endif

      PROCEDURE,  PUBLIC   :: INIT_POLAR => EPI_INIT_POLAR
      PROCEDURE,  PRIVATE  :: SETUP_COUPLING_WANN_LONG => EPI_SETUP_COUPLING_WANN_LONG
      PROCEDURE,  PRIVATE  :: ADD_COUPLING_LONG => EPI_ADD_COUPLING_LONG
      PROCEDURE,  PRIVATE  :: SUBTRACT_COUPLING_WANN_LONG => EPI_SUBTRACT_COUPLING_WANN_LONG
   END TYPE

   INTERFACE AVERAGE_DEGEN
      MODULE PROCEDURE AVERAGE_DEGEN_Z1
      MODULE PROCEDURE AVERAGE_DEGEN_Z2
      MODULE PROCEDURE AVERAGE_DEGEN_Z3
   END INTERFACE


! --------------------------------------
! --------------------------------------
! --------------------------------------
! --------------------------------------
! _test_mani
! Brillouin stuff! -> should either be moved, removed or merged.

   integer, parameter   :: num_neighbors = 26
   real(q), parameter   :: brill_tiny = 1e-10_q

   type brill_zone_t
   ! _test_mani protected?
      real(q), allocatable, public  :: qmesh(:, :)

      type(latt), private           :: lattice
      real(q), allocatable, private :: mid_point(:, :)
      real(q), allocatable, private :: mid_sqr(:)

   contains
      procedure,  public   :: init
      final                :: destructor

      procedure,  public   :: sample => sample_brill_zone
      procedure,  public   :: map_to_zone

      procedure,  private  :: generate_mid_points
      procedure,  private  :: check_point
      procedure,  private  :: calculate_shift
      procedure,  private  :: write_qmesh
   endtype brill_zone_t

CONTAINS


!                                 _____________
!________________________________/ ELPH_READER \________________________________
!
!> @brief Reads settings from the INCAR file
!>
!> # The following INCAR tags are available:
!>
!> ## LWSWQ
!> Calculates the overlop between current orbitals and orbitals stored in
!> file "WAVECAR.qqq"
!>
!> ## LELPH_BANDS
!> Interpolates the electronic band structure using Wannier orbitals
!>
!> ## LELPH_ZPR
!> Calculates the zero-point renormalization of the band gap
!>
!> ## LELPH_TEMP_GAP
!> Calculates the temperature-dependent band gap
!>
!> ## ELPH_METHOD
!> Wannierization method, see READ_ELPHON_METHOD
!>
!> ## NUM_WANN
!> Number of requested Wannier functions
!>
!> ## SCDM_MU
!> Location of the smooth cutoff for the SCDM method
!>
!> ## SCDM_SIGMA
!> Broadening of the smooth cutoff for the SCDM method
!>
!> ## ELPH_SUBDIV
!> Subdivision of q-point mesh used for integrals
!>
!> ## ELPH_MAX_BAND
!> Maximum band index to be included in the sum over states for the ZPR
!>
!> ## LELPH_ZPR_DYNAMIC
!> Includes phonon frequencies when calculating the ZPR
!>
!> ## ELPH_BZ_POINTS
!> ZPR is evaluated at these Brillouin-zone points. Points are given in
!> a space-separated list. Each 3-tuple is a point given in direct
!> coordinates of the reciprocal lattice of the primitive cell
!>
!> ## LELPH_WRITE_DATA
!> Writes electron-phonon data to HDF5
!>
!> ## LELPH_POLAR
!> Includes dipole-dipole corrections. For this mode, LPHON_POLAR must
!> also be set
!>
!> ## ELPH_G_CUTOFF
!> Works the same as PHON_G_CUTOFF but for electron-phonon calculations
!_______________________________________________________________________________
   SUBROUTINE ELPH_READER(INCAR, SETTINGS)
      USE incar_reader, ONLY: incar_file, &
                              process_incar, &
                              count_elements
      IMPLICIT NONE

      TYPE(incar_file),       INTENT(INOUT)  :: INCAR
      TYPE(ELPHON_SETTINGS),  INTENT(OUT)    :: SETTINGS

      INTEGER  :: NUM_ELEMENTS

      SETTINGS%METHOD = READ_ELPHON_METHOD(INCAR)

      ! For SCDM, we read additional parameters _test_mani move to wannier related stuff
      IF (SETTINGS%METHOD == METHOD_SCDM) THEN
         CALL process_incar(INCAR, "SCDM_MU", SETTINGS%SCDM_MU)
         CALL process_incar(INCAR, "SCDM_SIGMA", SETTINGS%SCDM_SIGMA)
      ENDIF

      CALL process_incar(INCAR, "LWSWQ", SETTINGS%DO_WSWQ)
      
      ! Do any kind of calculation?
      CALL process_incar(INCAR, "LELPH_BANDS", SETTINGS%CALC_BANDSTRUCTURE)
      CALL process_incar(INCAR, "LELPH_ZPR", SETTINGS%CALC_ZPR)
      CALL process_incar(INCAR, "LELPH_TEMP_GAP", SETTINGS%CALC_TEMPERATURE_GAP)
      CALL process_incar(INCAR, "LELPH_WRITE_DATA", SETTINGS%WRITE_DATA)
      CALL process_incar(INCAR, "LELPH_WRITE_WANN", SETTINGS%WRITE_WANN_FUN)

      ! Calculate the AE electron-phonon matrix element?
      CALL process_incar(INCAR, "LELPH_ALL_ELEC", SETTINGS%ALL_ELECTRON)

      ! Include corrections for polar materials?
      CALL process_incar(INCAR, "LELPH_POLAR", SETTINGS%DO_POLAR)
      CALL process_incar(INCAR, "ELPH_G_CUTOFF", SETTINGS%POLAR_G_CUTOFF)

      IF (SETTINGS%DO_INIT() .OR. SETTINGS%DO_PERTURB()) &
         CALL process_incar(INCAR, "NUM_WANN", SETTINGS%NUM_WANN_READ)

      IF (SETTINGS%DO_CALC()) THEN
         CALL process_incar(INCAR, "ELPH_SUBDIV", SETTINGS%BZ_SUBDIVISION)
         CALL process_incar(INCAR, "ELPH_MAX_BAND", SETTINGS%NUM_BANDS_MAX)
         CALL process_incar(INCAR, "LELPH_ZPR_DYNAMIC", SETTINGS%DO_ZPR_DYNAMIC)

         ! Read Brillouin-zone vectors
         NUM_ELEMENTS = count_elements(INCAR, "ELPH_BZ_POINTS")
         IF (NUM_ELEMENTS > 0) THEN
            IF (MODULO(NUM_ELEMENTS, 3) /= 0) CALL vtutor%error( &
               "The INCAR tag ELPH_BZ_POINTS only accepts a space-separated list of 3-vectors. " // &
               "Number of elements given: " // str(NUM_ELEMENTS))

            ALLOCATE(SETTINGS%BZ_POINTS(3, NUM_ELEMENTS / 3))
            CALL process_incar(INCAR, "ELPH_BZ_POINTS", SETTINGS%BZ_POINTS)
         ELSE
            ! Default to Gamma
            ALLOCATE(SETTINGS%BZ_POINTS(3, 1))
            SETTINGS%BZ_POINTS = 0
         ENDIF
      ENDIF

      ! Check consistency of input parameters
      IF (SETTINGS%WRITE_DATA .AND. SETTINGS%METHOD == METHOD_LOAD_DATA) CALL vtutor%error( &
         "LELPH_WRITE_DATA = .TRUE. cannot be used in conjunction with ELPH_METHOD = LOAD_DATA")

      ! Make sure a method is selected when we initialize
      IF (SETTINGS%METHOD == METHOD_NONE .AND. SETTINGS%DO_INIT()) CALL vtutor%error( &
         "Electron-phonon calculations were requested but ELPH_METHOD was not supplied.")
   END SUBROUTINE


!                             ____________________
!____________________________/ READ_ELPHON_METHOD \_____________________________
!
!> @brief Reads ELPH_METHOD tag from the INCAR file.
!>
!> # The following values are available:
!>
!> ## LOC_PROJ
!> Uses Local projections via LOCPROJ tag (see locproj.F)
!>
!> ## SCDM
!> Uses the SCDM method in the supercell (largly untested)
!>
!> ## LOAD_PRIM
!> Loads local orbitals from a previous primitive-cell calculation
!>
!> ## LOAD_DATA
!> Loads electron-phonon data from disk
!> 
!_______________________________________________________________________________
   INTEGER(KIND(METHOD_NONE)) FUNCTION READ_ELPHON_METHOD(INCAR) RESULT(METHOD)
      USE incar_reader, ONLY: incar_file, PROCESS_INCAR
      USE string,       ONLY: lowercase
      IMPLICIT NONE

      TYPE(incar_file), INTENT(INOUT) :: INCAR

      CHARACTER(:), ALLOCATABLE  :: READ_STRING

      METHOD = METHOD_NONE

      CALL PROCESS_INCAR(INCAR, 'ELPH_METHOD', READ_STRING)
      IF (.NOT.ALLOCATED(READ_STRING)) RETURN
      READ_STRING = lowercase(READ_STRING)

      SELECT CASE (READ_STRING)
         CASE ("loc_proj") ! _test_mani "locproj" broken
            METHOD = METHOD_LOCPROJ
         CASE ("scdm")
            METHOD = METHOD_SCDM
         CASE ("load_prim")
            METHOD = METHOD_LOAD_PRIM
         CASE ("load_data")
            METHOD = METHOD_LOAD_DATA

            ! Check to see if we have HDF5 support
#ifndef VASP_HDF5
            CALL vtutor%error('ELPH_METHOD = "' // READ_STRING // '" requires HDF5 support.')
#endif
         CASE DEFAULT
            CALL vtutor%error('ELPH_METHOD = "' // READ_STRING // '" was not recognized.')
      END SELECT
   END FUNCTION


   PURE LOGICAL FUNCTION SETTINGS_DO_INIT(THIS) RESULT(D)
      IMPLICIT NONE
      CLASS(ELPHON_SETTINGS), INTENT(IN)  :: THIS
      D = THIS%CALC_BANDSTRUCTURE .OR. THIS%WRITE_DATA .OR. THIS%DO_PERTURB()
      D = D .AND. .NOT. THIS%DO_LOAD_DATA()
   END FUNCTION

   PURE LOGICAL FUNCTION SETTINGS_DO_PERTURB(THIS) RESULT(D)
      IMPLICIT NONE
      CLASS(ELPHON_SETTINGS), INTENT(IN)  :: THIS
      D = THIS%CALC_ZPR .OR. THIS%CALC_TEMPERATURE_GAP .OR. THIS%WRITE_DATA .OR. THIS%DO_POLAR .OR. THIS%EXTERNAL_INTERFACE
      D = D .AND. .NOT. THIS%DO_LOAD_DATA()
   END FUNCTION

   PURE LOGICAL FUNCTION SETTINGS_DO_CALC(THIS) RESULT(D)
      IMPLICIT NONE
      CLASS(ELPHON_SETTINGS), INTENT(IN)  :: THIS
      D = THIS%CALC_BANDSTRUCTURE .OR. THIS%CALC_ZPR .OR. THIS%CALC_TEMPERATURE_GAP
   END FUNCTION

   PURE LOGICAL FUNCTION SETTINGS_DO_LOAD_DATA(THIS) RESULT(D)
      IMPLICIT NONE
      CLASS(ELPHON_SETTINGS), INTENT(IN)  :: THIS
      D = THIS%METHOD == METHOD_LOAD_DATA
   END FUNCTION


   SUBROUTINE ELPH_CHECK_INPUT_CONSISTENCY(INCAR, PHON_SETTINGS, ELPH_SETTINGS, POLAR_DATA)
      USE incar_reader, ONLY: incar_file, &
                              process_incar
      USE phonon,       ONLY: PHONON_SETTINGS
      IMPLICIT NONE

      TYPE(incar_file),       INTENT(INOUT)  :: INCAR
      TYPE(PHONON_SETTINGS),  INTENT(IN)     :: PHON_SETTINGS
      TYPE(ELPHON_SETTINGS),  INTENT(IN)     :: ELPH_SETTINGS
      TYPE(POLAR_DATA_T),     INTENT(IN)     :: POLAR_DATA

      INTEGER  :: IBRION
      INTEGER  :: IERR

      IF (ELPH_SETTINGS%DO_LOAD_DATA()) RETURN

      CALL process_incar(INCAR, "IBRION", IBRION, IERR)

      IF (IERR /= 0 .OR. IBRION /= 6) THEN
         IF (PHON_SETTINGS%DO_INIT().AND..NOT.PHON_SETTINGS%READ_FORCE_CONSTANTS) CALL vtutor%error( &
            "Phonon interpolations require finite differences using symmetry in the super cell. " // &
            "Set IBRION = 6 in the INCAR file.")
         IF (ELPH_SETTINGS%DO_PERTURB()) CALL vtutor%error( &
            "Electron-phonon interpolations require finite differences using symmetry in the super cell. " // &
            "Set IBRION = 6 in the INCAR file.")
      ENDIF

      IF (.NOT. ALLOCATED(POLAR_DATA%BORN_EFF_CHARGES)) THEN
         IF (PHON_SETTINGS%DO_POLAR) CALL vtutor%error( &
            "Polar phonon calculations require the Born effective charges and the dielectrict tensor. " // &
            "Please, supply them in the INCAR file.")
         IF (ELPH_SETTINGS%DO_POLAR) CALL vtutor%error( &
            "Polar electron-phonon calculations require the Born effective charges and the dielectrict tensor. " // &
            "Please, supply them in the INCAR file.")
      ENDIF
   END SUBROUTINE


!                                   _________
!__________________________________/ ELPH_IO \__________________________________
!
!> @brief Reads or writes electron-phonon data from/to disk
!>
!> Data is stored in an HDF5 file and can be used to perform Wannier
!> interpolation without having to recalculate the matrix elements.
!_______________________________________________________________________________
   SUBROUTINE ELPH_IO(EP_INTERPOL, SETTINGS, WDES, PRIM_CELL, NONLR_S, NONL_S, IO)
      USE nonlr_struct_def,   ONLY: nonlr_struct
      USE nonl_struct_def,    ONLY: nonl_struct
#ifdef VASP_HDF5
      USE hdf5
      USE vhdf5_base,         ONLY: vh5_file_create_or_overwrite, &
                                    vh5_file_close_writing, &
                                    vh5_file_open_read, &
                                    vh5_file_close, &
                                    vh5_error
      USE phonon,             ONLY: PHON_WRITE_HDF5, &
                                    PHON_READ_HDF5
#endif
      IMPLICIT NONE

      TYPE(ELPHON_INTERPOLATOR), INTENT(INOUT)  :: EP_INTERPOL
      TYPE(ELPHON_SETTINGS),     INTENT(IN)     :: SETTINGS
      TYPE(wavedes),             INTENT(IN)     :: WDES
      TYPE(PRIM_CELL_T),         INTENT(IN)     :: PRIM_CELL
      TYPE(nonlr_struct),        INTENT(IN)     :: NONLR_S
      TYPE(nonl_struct),         INTENT(IN)     :: NONL_S
      TYPE(in_struct),           INTENT(IN)     :: IO

      LOGICAL        :: FILE_EXISTS
#ifdef VASP_HDF5
      INTEGER(HID_T) :: FILE_ID
#endif

      ! Write data to file
      IF (SETTINGS%WRITE_DATA) THEN
         IF (SETTINGS%METHOD == METHOD_LOAD_DATA) CALL vtutor%bug( &
            "LELPH_WRITE_DATA = .TRUE. cannot be used in conjunction with ELPH_METHOD = LOAD_DATA", __FILE__, __LINE__)
         IF (.NOT. EP_INTERPOL%IS_INITIALIZED) CALL vtutor%bug("LELPH_WRITE_DATA but module not initialized!", __FILE__, __LINE__)
         IF (.NOT. EP_INTERPOL%IS_COMPLETED) CALL vtutor%bug("LELPH_WRITE_DATA but calculation not completed!", __FILE__, __LINE__)

#ifdef VASP_HDF5
         ! Write HDF5 data
         VH5_CHECK(vh5_file_create_or_overwrite('elphon_data.h5', FILE_ID))
         CALL EP_INTERPOL%WRITE_HDF5(FILE_ID, 'elphon')
         CALL PHON_WRITE_HDF5(FILE_ID, 'phonon')
         VH5_CHECK(vh5_file_close_writing(FILE_ID))
#else
         CALL vtutor%alert("LELPH_WRITE_DATA was set but HDF5 support is not enabled. No electron-phonon data is written.")
#endif
      ENDIF

      ! Write Wannier functions
      IF (SETTINGS%WRITE_WANN_FUN) CALL EP_INTERPOL%WRITE_WANNIER_FUNCS(IO)
      
      ! Read data from file
      IF (SETTINGS%METHOD == METHOD_LOAD_DATA) THEN
#ifdef VASP_HDF5
         ! Read HDF5 data
         INQUIRE(FILE = 'elphon_data.h5', EXIST = FILE_EXISTS)
         IF (.NOT. FILE_EXISTS) CALL vtutor%error( &
            "ELPH_IO: Attempt to read from 'elphon_data.h5' failed. File does not exist. " // &
            "Electron-phonon data must first be written by setting 'LELPH_WRITE_DATA' in the INCAR file.")

         VH5_CHECK(vh5_file_open_read('elphon_data.h5', FILE_ID))
         CALL EP_INTERPOL%READ_HDF5(PRIM_CELL, FILE_ID, 'elphon')
         CALL PHON_READ_HDF5(PRIM_CELL, FILE_ID, 'phonon')
         VH5_CHECK(vh5_file_close(FILE_ID))

         ! Read Wannier functions
         ! CALL EP_INTERPOL%READ_WANNIER_FUNCS(WDES, NONLR_S, NONL_S, IO)

         IF (SETTINGS%DO_CALC()) THEN
            CALL EP_INTERPOL%WANNIER_DESCRIPTOR%CALC_MIN_IMAGES()
            CALL EP_INTERPOL%WANNIER_DESCRIPTOR%CALC_MIN_IMAGES_TO_ATOM()
         ENDIF
#else
         CALL vtutor%error("ELPH_METHOD was set to 'LOAD_DATA' but HDF5 support is not enabled.")
#endif
      ENDIF
   END SUBROUTINE ELPH_IO


   ! _test_mani write wannier functions to disk.
   SUBROUTINE EPI_WRITE_WANNIER_FUNCS(THIS, IO)
      USE fileio, ONLY: OUTWAV
      IMPLICIT NONE

      CLASS(ELPHON_INTERPOLATOR),   INTENT(IN)  :: THIS
      TYPE(in_struct),              INTENT(IN)  :: IO

      CALL OUTWAV(IO, THIS%W_WAN_INIT%WDES, THIS%W_WAN_INIT, THIS%WANNIER_DESCRIPTOR%PRIM_CELL%SUPER_LATT, 0._q, 'wan')
   END SUBROUTINE


   ! _test_mani read wannier functions from disk.
   SUBROUTINE EPI_READ_WANNIER_FUNCS(THIS, WDES, NONLR_S, NONL_S, IO)
      USE nonlr_struct_def,   ONLY: nonlr_struct
      USE nonl_struct_def,    ONLY: nonl_struct
      USE fileio,             ONLY: INWAV_ALTERNATIVE
      USE nonl_high,          ONLY: PROALL
      IMPLICIT NONE

      CLASS(ELPHON_INTERPOLATOR),   INTENT(INOUT)  :: THIS
      TYPE(wavedes),                INTENT(IN)     :: WDES
      TYPE(nonlr_struct),           INTENT(IN)     :: NONLR_S
      TYPE(nonl_struct),            INTENT(IN)     :: NONL_S
      TYPE(in_struct),              INTENT(IN)     :: IO

      LOGICAL  :: LWREAD

      CALL INWAV_ALTERNATIVE(IO, WDES, THIS%W_WAN_INIT, WDES%GRID, THIS%WANNIER_DESCRIPTOR%PRIM_CELL%SUPER_LATT, LWREAD, 'wan')

      IF (.NOT.LWREAD) CALL vtutor%error('EPI_READ_WANNIER_FUNCS: Unable to read WAVECAR.wan.')

      CALL PROALL(WDES%GRID, THIS%WANNIER_DESCRIPTOR%PRIM_CELL%SUPER_LATT, NONLR_S, NONL_S, THIS%W_WAN_INIT)
   END SUBROUTINE


   ! _test_mani Serialization routine to write precalculated data to disk.
#ifdef VASP_HDF5
   SUBROUTINE EPI_WRITE_HDF5(THIS, LOC_ID, GROUP_PATH)
      USE hdf5
      USE vhdf5_base,   ONLY: vh5_error, &
                              vh5_write
      IMPLICIT NONE

      CLASS(ELPHON_INTERPOLATOR),   INTENT(IN)  :: THIS
      INTEGER(HID_T),               INTENT(IN)  :: LOC_ID
      CHARACTER(*),                 INTENT(IN)  :: GROUP_PATH

      INTEGER  :: N_SPIN, SPIN_IDX

      IF (.NOT. THIS%IS_INITIALIZED) CALL vtutor%bug("Module not initialized!", __FILE__, __LINE__)
      IF (.NOT. THIS%IS_COMPLETED) CALL vtutor%bug("Calculation not completed!", __FILE__, __LINE__)

      N_SPIN = SIZE(THIS%HAM_INIT)

      ! Write the Wannier descriptor
      CALL THIS%WANNIER_DESCRIPTOR%WRITE_HDF5(LOC_ID, GROUP_PATH // '/wannier_descriptor')
      ! Write number of spin components
      VH5_CHECK(vh5_write(LOC_ID, GROUP_PATH // '/num_spin', N_SPIN))
      ! Write the Wannier matrices
      DO SPIN_IDX = 1, N_SPIN
         CALL THIS%HAM_INIT(SPIN_IDX)%WRITE_HDF5(LOC_ID, GROUP_PATH // '/initial_hamiltonian/spin_' // str(SPIN_IDX))
         CALL THIS%COUPLING_WANN(SPIN_IDX)%WRITE_HDF5(LOC_ID, GROUP_PATH // '/coupling_wannier/spin_' // str(SPIN_IDX))
         CALL THIS%OVERLAP_PERT(SPIN_IDX)%WRITE_HDF5(LOC_ID, GROUP_PATH // '/perturbed_overlap/spin_' // str(SPIN_IDX))
      ENDDO
   END SUBROUTINE
   
   ! _test_mani Deserialization routine to fetch precalculated data from disk.
   SUBROUTINE EPI_READ_HDF5(THIS, PRIM_CELL, LOC_ID, GROUP_PATH)
      USE hdf5
      USE vhdf5_base,   ONLY: vh5_error, &
                              vh5_exists, &
                              vh5_read
      IMPLICIT NONE

      CLASS(ELPHON_INTERPOLATOR),   INTENT(INOUT)  :: THIS
      TYPE(PRIM_CELL_T),            INTENT(IN)     :: PRIM_CELL
      INTEGER(HID_T),               INTENT(IN)     :: LOC_ID
      CHARACTER(*),                 INTENT(IN)     :: GROUP_PATH

      INTEGER  :: SPIN_IDX, N_SPIN

      CALL THIS%DEALLOC()

      ! Read the Wannier descriptor
      CALL THIS%WANNIER_DESCRIPTOR%READ_HDF5(PRIM_CELL, LOC_ID, GROUP_PATH // '/wannier_descriptor')
      ! Read number of spin components
      VH5_CHECK(vh5_read(LOC_ID, GROUP_PATH // '/num_spin', N_SPIN))

      ! Setup allocation
      ALLOCATE(THIS%HAM_INIT(N_SPIN), THIS%COUPLING_WANN(N_SPIN), THIS%OVERLAP_PERT(N_SPIN))

      ! Read the Wannier matrices
      DO SPIN_IDX = 1, N_SPIN
         CALL THIS%HAM_INIT(SPIN_IDX)%READ_HDF5(THIS%WANNIER_DESCRIPTOR, LOC_ID, GROUP_PATH // '/initial_hamiltonian/spin_' // str(SPIN_IDX))
         CALL THIS%COUPLING_WANN(SPIN_IDX)%READ_HDF5(THIS%WANNIER_DESCRIPTOR, LOC_ID, GROUP_PATH // '/coupling_wannier/spin_' // str(SPIN_IDX))
         CALL THIS%OVERLAP_PERT(SPIN_IDX)%READ_HDF5(THIS%WANNIER_DESCRIPTOR, LOC_ID, GROUP_PATH // '/perturbed_overlap/spin_' // str(SPIN_IDX))
      ENDDO

      THIS%IS_INITIALIZED = .TRUE.
      THIS%IS_COMPLETED = .TRUE.
   END SUBROUTINE
#endif

   ! _test_mani interface code

   SUBROUTINE SERIALIZE_EXTERNAL_INTERFACE(PRIM_CELL, T_INFO, IO)
      IMPLICIT NONE

      TYPE(PRIM_CELL_T),   INTENT(IN)  :: PRIM_CELL
      TYPE(type_info),     INTENT(IN)  :: T_INFO
      TYPE(in_struct),     INTENT(IN)  :: IO

      INTEGER  :: UNIT

      OPEN(NEWUNIT = UNIT, FILE = 'binary_elphon_interface', STATUS = 'replace', FORM = 'unformatted')
      CALL PRIM_CELL%SERIALIZE(IO, UNIT)

      ! Type information
      WRITE(UNIT) SIZE(T_INFO%POMASS)
      WRITE(UNIT) T_INFO%POMASS(:)

      WRITE(UNIT) SIZE(T_INFO%ITYP)
      WRITE(UNIT) T_INFO%ITYP(:)

      CLOSE(UNIT)
   END SUBROUTINE


   ! _test_mani testing methods
   SUBROUTINE TEST_WANNIER_MAT_TRANS_SYMM(WANNIER_MAT, DESCRIPTOR, FILE_NAME, IO)
      IMPLICIT NONE

      GDEF,                         INTENT(IN)  :: WANNIER_MAT(:, :)
      TYPE(WANNIER_DESCRIPTOR_T),   INTENT(IN)  :: DESCRIPTOR
      CHARACTER(*),                 INTENT(IN)  :: FILE_NAME
      TYPE(in_struct),              INTENT(IN)  :: IO

      INTEGER  :: CELL_IDX_1, CELL_IDX_2, CELL_IDX_0
      INTEGER  :: UNIT
      LOGICAL  :: FOUND

      IF (IO%IU0 < 0) RETURN

      OPEN(NEWUNIT = UNIT, FILE = FILE_NAME, STATUS = 'replace')

      ASSOCIATE( &
         PC => DESCRIPTOR%PRIM_CELL, &
         N_ORBS => DESCRIPTOR%NUM_PRIM_ORBITALS)

         DO CELL_IDX_2 = 1, PC%NUM_CELLS
            DO CELL_IDX_1 = 1, PC%NUM_CELLS
               CELL_IDX_0 = PC%EQUIV_CELLS(CELL_IDX_1, CELL_IDX_2)

               WRITE(UNIT, '(F12.6, " ")', ADVANCE = 'no') MAXVAL(ABS( &
                  ABS(WANNIER_MAT( &
                     (CELL_IDX_1 - 1) * N_ORBS + 1 : CELL_IDX_1 * N_ORBS, &
                     (CELL_IDX_2 - 1) * N_ORBS + 1 : CELL_IDX_2 * N_ORBS)) - &
                  ABS(WANNIER_MAT( &
                     (CELL_IDX_0 - 1) * N_ORBS + 1 : CELL_IDX_0 * N_ORBS, 1 : N_ORBS))))
            ENDDO
            WRITE(UNIT, *)
         ENDDO
      END ASSOCIATE

      CLOSE(UNIT)
   END SUBROUTINE TEST_WANNIER_MAT_TRANS_SYMM

   LOGICAL FUNCTION READ_SINGLE_RECORD(UNIT, STRING) RESULT(SUCCESSFUL)
      IMPLICIT NONE

      INTEGER,                   INTENT(IN)  :: UNIT
      CHARACTER(:), ALLOCATABLE, INTENT(OUT) :: STRING

      CHARACTER(256) :: BUFFER
      INTEGER        :: SIZE

      STRING = ''

      DO
         READ(UNIT, '(A)', ADVANCE = 'no', SIZE = SIZE, EOR = 10, END = 20, ERR = 20) BUFFER
         STRING = STRING // BUFFER
      ENDDO

      10 CONTINUE
      STRING = STRING // BUFFER(1 : SIZE)
      SUCCESSFUL = .TRUE.
      RETURN

      20 CONTINUE
      SUCCESSFUL = .FALSE.
      RETURN
   END FUNCTION READ_SINGLE_RECORD


   SUBROUTINE READ_WANNIER_CENTERS(FILE_NAME, CENTERS)
      IMPLICIT NONE

      CHARACTER(*),           INTENT(IN)  :: FILE_NAME
      REAL(q), ALLOCATABLE,   INTENT(OUT) :: CENTERS(:, :)

      INTEGER                    :: UNIT
      CHARACTER(:), ALLOCATABLE  :: LINE
      INTEGER                    :: NUM_CENTERS
      INTEGER                    :: I
      CHARACTER(0)               :: TMP_STR

      OPEN(NEWUNIT = UNIT, FILE = FILE_NAME, STATUS = 'old')

      ! Read number of centers
      ! _test_mani exclude atomic centers
      IF (.NOT. READ_SINGLE_RECORD(UNIT, LINE)) CALL WANNIER_CENTER_ERROR()
      NUM_CENTERS = 0
      READ(LINE, *, ERR = 10, END = 10) NUM_CENTERS

      IF (NUM_CENTERS < 1) CALL WANNIER_CENTER_ERROR()
      ALLOCATE(CENTERS(3, NUM_CENTERS))

      ! Skip comment line
      IF (.NOT. READ_SINGLE_RECORD(UNIT, LINE)) CALL WANNIER_CENTER_ERROR()

      ! Read Wannier centers
      DO I = 1, NUM_CENTERS
         IF (.NOT. READ_SINGLE_RECORD(UNIT, LINE)) CALL WANNIER_CENTER_ERROR()
         READ(LINE, *, ERR = 10, END = 10) TMP_STR, CENTERS(:, I)
      ENDDO

      RETURN

      10 CONTINUE
      CALL WANNIER_CENTER_ERROR()

      CONTAINS

      ! _test_mani implement
      SUBROUTINE WANNIER_CENTER_ERROR()
         IMPLICIT NONE
         CALL vtutor%error("READ_WANNIER_CENTERS: Error reading Wannier centers from file " // FILE_NAME // ".")
      END SUBROUTINE WANNIER_CENTER_ERROR
   END SUBROUTINE READ_WANNIER_CENTERS


   SUBROUTINE READ_3D_FIELD_FROM_XSF(FILE_NAME, FIELD, GRID_LATTICE, GRID_OFFSET, PRIM_CELL, IO, PRIMITIVE_TRAFO)
      IMPLICIT NONE

      CHARACTER(*),           INTENT(IN)  :: FILE_NAME
      REAL(q), ALLOCATABLE,   INTENT(OUT) :: FIELD(:, :, :)
      REAL(q),                INTENT(OUT) :: GRID_LATTICE(3, 3)
      REAL(q),                INTENT(OUT) :: GRID_OFFSET(3)
      TYPE(PRIM_CELL_T),      INTENT(IN)  :: PRIM_CELL
      TYPE(in_struct),        INTENT(IN)  :: IO
      REAL(q), OPTIONAL,      INTENT(OUT) :: PRIMITIVE_TRAFO(3, 3)

      INTEGER                    :: UNIT
      CHARACTER(:), ALLOCATABLE  :: LINE
      INTEGER                    :: I
      REAL(q)                    :: PRIM_LATTICE(3, 3)
      REAL(q)                    :: PRIM_TRAFO(3, 3)
      INTEGER                    :: NUM_GRID_POINTS(3)
      LOGICAL                    :: READ_PRIM
      LOGICAL                    :: READ_FIELD

      READ_PRIM = .FALSE.
      READ_FIELD = .FALSE.

      OPEN(NEWUNIT = UNIT, FILE = FILE_NAME, STATUS = 'old')
      DO WHILE (READ_SINGLE_RECORD(UNIT, LINE))
         ! _test_mani lower/upper case
         LINE = TRIM(ADJUSTL(LINE))

         SELECT CASE (LINE)
         CASE("PRIMVEC")
            DO I = 1, 3
               IF (.NOT. READ_SINGLE_RECORD(UNIT, LINE)) CALL XSF_ERROR()
               READ(LINE, *, ERR = 10, END = 10) PRIM_LATTICE(:, I)
            ENDDO

            ! _test_mani improve trafo determination
            PRIM_TRAFO = MATMUL(PRIM_CELL%PRIM_LATT%A, MATRIX_INVERSE_3X3(PRIM_LATTICE))
            IF (ABS(ABS(DET_3X3(PRIM_TRAFO)) - 1) > 1E-5_q) &
               CALL vtutor%error("READ_3D_FIELD_FROM_XSF: Primitive cells do not match!")

            IF (PRESENT(PRIMITIVE_TRAFO)) PRIMITIVE_TRAFO = PRIM_TRAFO

            READ_PRIM = .TRUE.
         
         CASE("BEGIN_BLOCK_DATAGRID_3D")
            ! Early exit if primitive data is missing
            IF (.NOT. READ_PRIM) CALL XSF_ERROR()

            ! Skip comment
            IF (.NOT. READ_SINGLE_RECORD(UNIT, LINE)) CALL XSF_ERROR()

            ! Check if beginning of block is correct
            IF (.NOT. READ_SINGLE_RECORD(UNIT, LINE)) CALL XSF_ERROR()
            LINE = TRIM(ADJUSTL(LINE))
            IF (LEN(LINE) < 17 .OR. LINE(1 : 17) /= "BEGIN_DATAGRID_3D") CALL XSF_ERROR()

            ! Read number of grid points
            IF (.NOT. READ_SINGLE_RECORD(UNIT, LINE)) CALL XSF_ERROR()
            READ(LINE, *, ERR = 10, END = 10) NUM_GRID_POINTS(:)

            ! Read offset of grid origin
            IF (.NOT. READ_SINGLE_RECORD(UNIT, LINE)) CALL XSF_ERROR()
            READ(LINE, *, ERR = 10, END = 10) GRID_OFFSET(:)

            ! Read grid vectors
            DO I = 1, 3
               IF (.NOT. READ_SINGLE_RECORD(UNIT, LINE)) CALL XSF_ERROR()
               READ(LINE, *, ERR = 10, END = 10) GRID_LATTICE(:, I)
               GRID_LATTICE(:, I) = GRID_LATTICE(:, I) / (NUM_GRID_POINTS(I) - 1) * NUM_GRID_POINTS(I)
            ENDDO

            ! Transform lattice and offset to the VASP coordinate system
            GRID_LATTICE = MATMUL(PRIM_TRAFO, GRID_LATTICE)
            GRID_OFFSET = MATMUL(PRIM_TRAFO, GRID_OFFSET)

            ! Read field
            ALLOCATE(FIELD(NUM_GRID_POINTS(1), NUM_GRID_POINTS(2), NUM_GRID_POINTS(3)))
            READ(UNIT, *, ERR = 10, END = 10) FIELD(:, :, :)

            ! Check if end of block is correct
            IF (.NOT. READ_SINGLE_RECORD(UNIT, LINE)) CALL XSF_ERROR()
            LINE = TRIM(ADJUSTL(LINE))
            IF (LEN(LINE) < 15 .OR. LINE(1 : 15) /= "END_DATAGRID_3D") CALL XSF_ERROR()

            READ_FIELD = .TRUE.
            EXIT
         END SELECT
      ENDDO
      CLOSE(UNIT)

      IF (.NOT. READ_PRIM .OR. .NOT. READ_FIELD) CALL XSF_ERROR()

      RETURN

      10 CONTINUE
      CALL XSF_ERROR()

      CONTAINS

      ! _test_mani implement
      SUBROUTINE XSF_ERROR()
         IMPLICIT NONE
         CALL vtutor%error("READ_3D_FIELD_FROM_XSF: Error reading XSF file " // FILE_NAME // ".")
      END SUBROUTINE XSF_ERROR
   END SUBROUTINE READ_3D_FIELD_FROM_XSF


   SUBROUTINE SHIFT_GRID_TO_CENTER_FIELD(FIELD, GRID_OFFSET, GRID_LATTICE, FIELD_CENTER)
      IMPLICIT NONE

      REAL(q), INTENT(INOUT)  :: FIELD(:, :, :)
      REAL(q), INTENT(INOUT)  :: GRID_OFFSET(3)
      REAL(q), INTENT(IN)     :: GRID_LATTICE(3, 3)
      REAL(q), INTENT(IN)     :: FIELD_CENTER(3)

      REAL(q)  :: FIELD_CENTER_DIR(3)
      INTEGER  :: GRID_SHIFT(3)
      INTEGER  :: I

      FIELD_CENTER_DIR = MATMUL(MATRIX_INVERSE_3X3(GRID_LATTICE), FIELD_CENTER - GRID_OFFSET)
      GRID_SHIFT = NINT(.5_q - FIELD_CENTER_DIR)

      DO I = 1, 3
         FIELD = CSHIFT(FIELD, -GRID_SHIFT(I), I)
      ENDDO

      GRID_OFFSET = GRID_OFFSET - MATMUL(GRID_LATTICE, GRID_SHIFT)

   END SUBROUTINE SHIFT_GRID_TO_CENTER_FIELD


   SUBROUTINE CREATE_ORBITAL_FROM_FIELD(W1, LATT_CUR, NONL_S, NONLR_S, CQIJ, FIELD, GRID_LATTICE, GRID_OFFSET)
      USE wave_struct_def,    ONLY: wavefun1
      USE nonl_struct_def,    ONLY: nonl_struct
      USE nonlr_struct_def,   ONLY: nonlr_struct
      USE WAVE,               ONLY: NEWWAV, DELWAV
      USE wave_high,          ONLY: W1_COPY, CNORMN
      USE nonl_high,          ONLY: W1_PROJ
      IMPLICIT NONE

      TYPE(wavefun1)                   :: W1 ! _test_mani intent(inout) not working with elements()
      TYPE(latt),          INTENT(IN)  :: LATT_CUR
      TYPE(nonl_struct),   INTENT(IN)  :: NONL_S
      TYPE(nonlr_struct),  INTENT(IN)  :: NONLR_S
      OVERLAP,             INTENT(IN)  :: CQIJ(:, :, :,:)
      REAL(q),             INTENT(IN)  :: FIELD(:, :, :)
      REAL(q),             INTENT(IN)  :: GRID_LATTICE(3, 3)
      REAL(q),             INTENT(IN)  :: GRID_OFFSET(3)

      INTEGER              :: FIELD_GRID_SIZE(3)
      INTEGER              :: WAVE_GRID_SIZE(3)
      TYPE(wavefun1)       :: MY_W1
      COMPLEX(q), POINTER  :: CR_PTR(:, :, :)
      INTEGER              :: X, Y, Z
      INTEGER              :: SHIFT_X, SHIFT_Y, SHIFT_Z
      INTEGER              :: I1, I2, I3
      REAL(q)              :: GRID_TRAFO(3, 3)
      REAL(q)              :: FIELD_GRID_CENTER(3)
      REAL(q)              :: MIN_OFFSET(3)
      REAL(q)              :: POS_IN_WAVE_GRID(3), POS_IN_FIELD_GRID(3), POS_IN_BOX(3)
      REAL(q)              :: DIST_VEC(3)
      REAL(q)              :: DIST_SQR, MIN_DIST_SQR
      INTEGER              :: LOWER_BOUNDS(3), FIELD_IDX(3)
      REAL(q)              :: INTERPOL_BUFFER(0:1, 0:1, 0:1)
      REAL(q)              :: TMP

      IF (DET_3X3(GRID_LATTICE) > LATT_CUR%OMEGA + 1e-3_q) CALL vtutor%error("CREATE_ORBITAL_FROM_FIELD: Field grid too large.")

      FIELD_GRID_SIZE = SHAPE(FIELD)
      WAVE_GRID_SIZE = [W1%WDES1%GRID%NGX, W1%WDES1%GRID%NGY, W1%WDES1%GRID%NGZ]

      ! Get center of field grid relative to field-grid origin in super-cell coordinates
      FIELD_GRID_CENTER = MATMUL(MATMUL(GRID_LATTICE, [.5_q, .5_q, .5_q]), LATT_CUR%B)
      ! Get offset of field grid in super-cell coordinates
      MIN_OFFSET = MATMUL(GRID_OFFSET, LATT_CUR%B)
      ! Set offset such that field-grid center has shortest distance to super-cell origin
      MIN_OFFSET = MODULO(MIN_OFFSET + FIELD_GRID_CENTER + .5_q, [REAL(q):: 1, 1, 1]) - FIELD_GRID_CENTER - .5_q

      ! Get trafo between grid lattices
      GRID_TRAFO = MATMUL(MATRIX_INVERSE_3X3(GRID_LATTICE), LATT_CUR%A)

      ! Get offset in field-grid coordinates and scale by number of grid points
      MIN_OFFSET = MATMUL(GRID_TRAFO, MIN_OFFSET) * FIELD_GRID_SIZE

      ! Get trafo between actual grids
      GRID_TRAFO = GRID_TRAFO * DYAD(REAL(FIELD_GRID_SIZE, q), 1._q / WAVE_GRID_SIZE)

      CALL NEWWAV(MY_W1, W1%WDES1, .TRUE.)

      MY_W1%CR = 0
      CR_PTR( &
         0 : WAVE_GRID_SIZE(1) - 1, &
         0 : WAVE_GRID_SIZE(2) - 1, &
         0 : WAVE_GRID_SIZE(3) - 1) => MY_W1%CR

      DO Z = 0, WAVE_GRID_SIZE(3) - 1
         DO Y = 0, WAVE_GRID_SIZE(2) - 1
            DO X = 0, WAVE_GRID_SIZE(1) - 1

               MIN_DIST_SQR = HUGE(1._q)
               
               DO SHIFT_Z = -1, 0
                  DO SHIFT_Y = -1, 0
                     DO SHIFT_X = -1, 0
                        POS_IN_WAVE_GRID = [REAL(q):: X, Y, Z] + [REAL(q):: SHIFT_X, SHIFT_Y, SHIFT_Z] * WAVE_GRID_SIZE
                        POS_IN_FIELD_GRID = MATMUL(GRID_TRAFO, POS_IN_WAVE_GRID) - MIN_OFFSET
                        
                        LOWER_BOUNDS = FLOOR(POS_IN_FIELD_GRID)

                        IF (ANY(LOWER_BOUNDS + 1 <= 0 .OR. LOWER_BOUNDS + 1 > FIELD_GRID_SIZE)) CYCLE

                        DIST_VEC = MATMUL(GRID_LATTICE, .5_q * FIELD_GRID_SIZE - POS_IN_FIELD_GRID)
                        DIST_SQR = DOT_PRODUCT(DIST_VEC, DIST_VEC)
                        IF (DIST_SQR > MIN_DIST_SQR) CYCLE
                        MIN_DIST_SQR = DIST_SQR
                        
                        CR_PTR(X, Y, Z) = LERP_3D_0(FIELD, POS_IN_FIELD_GRID)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      ! CALL FFTEXT_WDES1(MY_W1%WDES1, MY_W1%CR, MY_W1%CPTWFP, .FALSE.) ! _test_mani can't use in GPU version
      CALL FFTEXT(MY_W1%WDES1%NGVECTOR, MY_W1%WDES1%NINDPW, MY_W1%CR, MY_W1%CPTWFP, MY_W1%WDES1%GRID, .FALSE.)

      CALL W1_COPY(MY_W1, W1)
      CALL DELWAV(MY_W1, .TRUE.)

      CALL W1_PROJ(W1, NONLR_S, NONL_S)
      CALL CNORMN(W1, CQIJ, 1, TMP)

   END SUBROUTINE CREATE_ORBITAL_FROM_FIELD


   PURE REAL(q) FUNCTION LERP_3D_0(FIELD, POS_IN_GRID) RESULT(VAL)
      IMPLICIT NONE

      REAL(q), INTENT(IN)  :: FIELD(:, :, :)
      REAL(q), INTENT(IN)  :: POS_IN_GRID(3)

      REAL(q)  :: BUFFER(0 : 1, 0 : 1, 0 : 1)
      INTEGER  :: GRID_SIZE(3)
      INTEGER  :: X, Y, Z
      INTEGER  :: LOWER_BOUNDS(3)
      INTEGER  :: FIELD_IDX(3)
      REAL(q)  :: POS_IN_BOX(3)

      GRID_SIZE = SHAPE(FIELD)
      LOWER_BOUNDS = FLOOR(POS_IN_GRID)
      POS_IN_BOX = POS_IN_GRID - LOWER_BOUNDS

      DO Z = 0, 1
         DO Y = 0, 1
            DO X = 0, 1
               FIELD_IDX = LOWER_BOUNDS + [X, Y, Z] + 1
               IF (ALL(FIELD_IDX > 0 .AND. FIELD_IDX <= GRID_SIZE)) THEN
                  BUFFER(X, Y, Z) = FIELD(FIELD_IDX(1), FIELD_IDX(2), FIELD_IDX(3))
               ELSE
                  BUFFER(X, Y, Z) = 0._q
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      BUFFER(:, :, 0) = BUFFER(:, :, 0) * (1._q - POS_IN_BOX(3)) + BUFFER(:, :, 1) * POS_IN_BOX(3)
      BUFFER(:, 0, 0) = BUFFER(:, 0, 0) * (1._q - POS_IN_BOX(2)) + BUFFER(:, 1, 0) * POS_IN_BOX(2)
      BUFFER(0, 0, 0) = BUFFER(0, 0, 0) * (1._q - POS_IN_BOX(1)) + BUFFER(1, 0, 0) * POS_IN_BOX(1)

      VAL = BUFFER(0, 0, 0)
   END FUNCTION LERP_3D_0


   SUBROUTINE IMPORT_WANNIER_FROM_PRIM(W, W_WANNIER, NUM_ORBS_IN_CELL, PRIM_CELL, LATT_CUR, NONL_S, NONLR_S, CQIJ, IO)
      USE nonl_struct_def,    ONLY: nonl_struct
      USE nonlr_struct_def,   ONLY: nonlr_struct
      USE wave,               ONLY: ALLOCW, &
                                    DEALLOCW, &
                                    SETWDES
      USE wave_high,          ONLY: ELEMENT, &
                                    ELEMENTS
      USE locproj,            ONLY: LPRJ_LINCOM
      IMPLICIT NONE

      TYPE(wavespin),      INTENT(IN)  :: W
      TYPE(wavespin),      INTENT(OUT) :: W_WANNIER
      INTEGER,             INTENT(IN)  :: NUM_ORBS_IN_CELL
      TYPE(PRIM_CELL_T),   INTENT(IN)  :: PRIM_CELL
      TYPE(latt),          INTENT(IN)  :: LATT_CUR
      TYPE(nonl_struct),   INTENT(IN)  :: NONL_S
      TYPE(nonlr_struct),  INTENT(IN)  :: NONLR_S
      OVERLAP,             INTENT(IN)  :: CQIJ(:, :, :, :)
      TYPE(in_struct),     INTENT(IN)  :: IO

      TYPE(wavedes1)             :: WDES1
      REAL(q), ALLOCATABLE       :: WANNIER_CENTERS(:, :)
      REAL(q), ALLOCATABLE       :: FIELD(:, :, :)
      REAL(q)                    :: GRID_LATTICE(3, 3)
      REAL(q)                    :: GRID_OFFSET(3)
      REAL(q)                    :: PRIM_TRAFO(3, 3)
      INTEGER                    :: FIELD_SIZE(3)
      CHARACTER(19)              :: FILE_NAME
      INTEGER                    :: ORB_IDX
      INTEGER                    :: BAND_IDX_GLOBAL, BAND_IDX_LOCAL

      CALL DEALLOCW(W_WANNIER)
      CALL ALLOCW(W%WDES, W_WANNIER)

      CALL SETWDES(W%WDES, WDES1, 1)

      CALL READ_WANNIER_CENTERS("wannier90_centres.xyz", WANNIER_CENTERS)

      DO ORB_IDX = 1, NUM_ORBS_IN_CELL
         DO BAND_IDX_LOCAL = 1, WDES1%NBANDS
            BAND_IDX_GLOBAL = (BAND_IDX_LOCAL - 1) * WDES1%NB_PAR + WDES1%NB_LOW
            IF (BAND_IDX_GLOBAL /= ORB_IDX) CYCLE

            WRITE(FILE_NAME, '(A, I5.5, A)') "wannier90_", ORB_IDX, ".xsf"

            CALL READ_3D_FIELD_FROM_XSF(FILE_NAME, FIELD, GRID_LATTICE, GRID_OFFSET, PRIM_CELL, IO, PRIM_TRAFO)

            WANNIER_CENTERS = MATMUL(PRIM_TRAFO, WANNIER_CENTERS)

            CALL SHIFT_GRID_TO_CENTER_FIELD(FIELD, GRID_OFFSET, GRID_LATTICE, WANNIER_CENTERS(:, ORB_IDX))
            FIELD_SIZE = SHAPE(FIELD)

            WRITE(*, *) "READ OK"

            CALL CREATE_ORBITAL_FROM_FIELD(ELEMENT(W_WANNIER, WDES1, BAND_IDX_LOCAL, 1), LATT_CUR, NONL_S, NONLR_S, CQIJ, FIELD, GRID_LATTICE, GRID_OFFSET)
            WRITE(*, *) "CREATE OK"
         ENDDO
      ENDDO
   END SUBROUTINE IMPORT_WANNIER_FROM_PRIM


   SUBROUTINE DUPLICATE_PRIM_ORBITALS(WA, NUM_ORBS_IN_CELL, PRIM_CELL)
      USE wave_struct_def, ONLY: wavefuna
      USE wave_high,       ONLY: REDISTRIBUTE_PW
      IMPLICIT NONE

      TYPE(wavefuna),      INTENT(INOUT)  :: WA
      INTEGER,             INTENT(IN)     :: NUM_ORBS_IN_CELL
      TYPE(PRIM_CELL_T),   INTENT(IN)     :: PRIM_CELL

      INTEGER  :: CELL_IDX

      ! _test_mani omg, find better method

      IF (WA%WDES1%DO_REDIS) CALL REDISTRIBUTE_PW(WA)

      DO CELL_IDX = 2, PRIM_CELL%NUM_CELLS
         WA%CW_RED(:, (CELL_IDX - 1) * NUM_ORBS_IN_CELL + 1 : CELL_IDX * NUM_ORBS_IN_CELL) = &
         WA%CW_RED(:, 1 : NUM_ORBS_IN_CELL)
      ENDDO

      IF (WA%WDES1%DO_REDIS) CALL REDISTRIBUTE_PW(WA)
   END SUBROUTINE DUPLICATE_PRIM_ORBITALS


   ! _test_mani rewrite to super_latt <-> prim_latt
   SUBROUTINE GET_GRID_POINTS_IN_PRIM_CELL(INDICES, GRID, LATT_CUR, PRIM_CELL)
      USE mgrid_struct_def,   ONLY: grid_3d
      IMPLICIT NONE

      INTEGER, ALLOCATABLE,   INTENT(OUT) :: INDICES(:)
      TYPE(grid_3d),          INTENT(IN)  :: GRID
      TYPE(latt),             INTENT(IN)  :: LATT_CUR
      TYPE(PRIM_CELL_T),      INTENT(IN)  :: PRIM_CELL

      REAL(q)              :: POINT_FROM_PC(3)
      REAL(q)              :: LATTICE_FROM_PC(3, 3)
      REAL(q)              :: TEST_POINT(3)
      INTEGER              :: X, Y, Z
      INTEGER              :: SHIFT_X, SHIFT_Y, SHIFT_Z
      INTEGER              :: CUR_INDEX
      INTEGER, ALLOCATABLE :: TMP_ARRAY(:)

      IF (ALLOCATED(INDICES)) DEALLOCATE(INDICES)
      ! Allocate array 10% larger than presumably required so reallocation is kept to a minimum.
      ALLOCATE(INDICES(INT(GRID%NPLWV * PRIM_CELL%PRIM_LATT%OMEGA * (1.1_q / LATT_CUR%OMEGA))))

      LATTICE_FROM_PC = MATMUL(TRANSPOSE(PRIM_CELL%PRIM_LATT%B), LATT_CUR%A)

      CUR_INDEX = 0

      DO Z = 0, GRID%NGZ - 1
         DO Y = 0, GRID%NGY - 1
            POINTS_LOOP: DO X = 0, GRID%NGX - 1

               POINT_FROM_PC = [REAL(q):: X, Y, Z] / [REAL(q):: GRID%NGX, GRID%NGY, GRID%NGZ]
               POINT_FROM_PC = MATMUL(MATMUL(LATT_CUR%A, POINT_FROM_PC), PRIM_CELL%PRIM_LATT%B)
      
               DO SHIFT_Z = -1, 0
                  DO SHIFT_Y = -1, 0
                     DO SHIFT_X = -1, 0

                        TEST_POINT = POINT_FROM_PC + MATMUL(LATTICE_FROM_PC, [REAL(q):: SHIFT_X, SHIFT_Y, SHIFT_Z])

                        ! _test_mani
                        TEST_POINT = TEST_POINT + 0._q / 8._q
                        
                        IF (ALL(TEST_POINT >= 0._q) .AND. ALL(TEST_POINT < 1._q)) THEN
                           ! Reallocate if necessary.
                           IF (CUR_INDEX == SIZE(INDICES)) THEN
                              ALLOCATE(TMP_ARRAY(CUR_INDEX + 64))
                              TMP_ARRAY(1 : CUR_INDEX) = INDICES
                              CALL MOVE_ALLOC(TMP_ARRAY, INDICES)
                           ENDIF

                           CUR_INDEX = CUR_INDEX + 1
                           INDICES(CUR_INDEX) = (GRID%RL%INDEX(Y, Z) - 1) * GRID%NGX + X + 1 ! _test_mani TODO: general case

                           CYCLE POINTS_LOOP
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO POINTS_LOOP
         ENDDO
      ENDDO
      
      ! Reallocate to correct size.
      IF (SIZE(INDICES) > CUR_INDEX) THEN
         ALLOCATE(TMP_ARRAY(CUR_INDEX))
         TMP_ARRAY(:) = INDICES(1 : CUR_INDEX)
         CALL MOVE_ALLOC(TMP_ARRAY, INDICES)
      ENDIF
      
   END SUBROUTINE GET_GRID_POINTS_IN_PRIM_CELL


   ! _test_mani TODO: general case.
   SUBROUTINE GATHER_PRIM_CELL_WAVE_FUNCTIONS(W, PRIM_CELL, LATT_CUR, WAVE_FUNC, BAND_INDICES)
      USE wave_struct_def, ONLY: wavefun1
      USE wave_high,       ONLY: FFTWAV_W1, &
                                 ELEMENT, &
                                 W1_COPY
      USE wave,            ONLY: SETWDES, &
                                 NEWWAV, &
                                 DELWAV
      IMPLICIT NONE

      TYPE(wavespin),         INTENT(IN)  :: W
      TYPE(PRIM_CELL_T),      INTENT(IN)  :: PRIM_CELL
      TYPE(latt),             INTENT(IN)  :: LATT_CUR
      GDEF, ALLOCATABLE,      INTENT(OUT) :: WAVE_FUNC(:, :)
      INTEGER, ALLOCATABLE,   INTENT(OUT) :: BAND_INDICES(:)

      TYPE(wavedes1)       :: WDES1
      TYPE(wavefun1)       :: W1
      INTEGER              :: NB_INC
      INTEGER, ALLOCATABLE :: BUFFER(:)
      INTEGER              :: NUM_POINTS
      INTEGER, ALLOCATABLE :: POINT_INDICES(:)
      INTEGER              :: PRIM_CELL_IDX
      INTEGER              :: I
      INTEGER              :: NB_LOW
      INTEGER              :: BAND_IDX

      CALL SETWDES(W%WDES, WDES1, 1)

      IF (WDES1%GRID%RL%NFAST /= 1) CALL vtutor%bug("GATHER_PRIM_CELL_WAVE_FUNCTIONS: NFAST /= 1 not implemented!", __FILE__, __LINE__)
#ifdef MPI
      IF (WDES1%GRID%COMM%NCPU /= 1) CALL vtutor%bug("GATHER_PRIM_CELL_WAVE_FUNCTIONS: NCPU /= 1 not implemented!", __FILE__, __LINE__)
#endif

      NB_INC = WDES1%NB_TOT
      ALLOCATE(BAND_INDICES(NB_INC))
      BAND_INDICES = [(I, I = 1, NB_INC)]

      ! _test_mani Improve MPI handling!
      NUM_POINTS = 0
#ifdef MPI
      IF (WDES1%COMM%NODE_ME == 1) THEN
#endif
         CALL GET_GRID_POINTS_IN_PRIM_CELL(POINT_INDICES, WDES1%GRID, LATT_CUR, PRIM_CELL)
         NUM_POINTS = SIZE(POINT_INDICES)
#ifdef MPI
      ENDIF

      CALL M_bcast_i(WDES1%COMM, NUM_POINTS, 1)
      IF (WDES1%COMM%NODE_ME /= 1) ALLOCATE(POINT_INDICES(NUM_POINTS))
      CALL M_bcast_i(WDES1%COMM, POINT_INDICES(1), NUM_POINTS)
#endif

      IF (ALLOCATED(WAVE_FUNC)) DEALLOCATE(WAVE_FUNC)
      ALLOCATE(WAVE_FUNC(NUM_POINTS, NB_INC))
      WAVE_FUNC = 0

      CALL NEWWAV(W1, WDES1, .true.)

      DO I = 1, NB_INC
         ! Check if we're on the correct node.
         NB_LOW = MODULO(BAND_INDICES(I) - 1, WDES1%NB_PAR) + 1
         IF (NB_LOW /= WDES1%NB_LOW) CYCLE
         ! Fetch local band index to index the array.
         BAND_IDX = (BAND_INDICES(I) - 1) / WDES1%NB_PAR + 1

         ! Copy the orbital and FFT to real-space wave function.
         CALL W1_COPY(ELEMENT(W, WDES1, BAND_IDX, 1), W1)
         CALL FFTWAV_W1(W1)

         ! Copy the relevant part of the wave function to the output array.
         WAVE_FUNC(:, I) = W1%CR(POINT_INDICES(:)) ! Imaginary part is zero in gamma-only version.
      ENDDO

      CALL DELWAV(W1, .TRUE.)

      CALLMPI(M_sum_g(W%WDES%COMM, WAVE_FUNC(1, 1), NUM_POINTS * NB_INC)) ! _test_mani comm?

   END SUBROUTINE GATHER_PRIM_CELL_WAVE_FUNCTIONS

   ! rewrite/remove
   SUBROUTINE GATHER_PRIM_CELL_WAVE_FUNCTIONS_3_27(W, PRIM_CELL, LATT_CUR, WAVE_FUNC, BAND_INDICES)
      USE wave_struct_def, ONLY: wavefun1
      USE wave_high,       ONLY: FFTWAV_W1, &
                                 ELEMENT, &
                                 W1_COPY
      USE wave,            ONLY: SETWDES, &
                                 NEWWAV, &
                                 DELWAV
      IMPLICIT NONE

      TYPE(wavespin),            INTENT(IN)  :: W
      TYPE(PRIM_CELL_T),         INTENT(IN)  :: PRIM_CELL
      TYPE(latt),                INTENT(IN)  :: LATT_CUR
      GDEF, ALLOCATABLE, TARGET, INTENT(OUT) :: WAVE_FUNC(:, :)
      INTEGER, ALLOCATABLE,      INTENT(OUT) :: BAND_INDICES(:)

      TYPE(wavedes1)       :: WDES1
      TYPE(wavefun1)       :: W1
      INTEGER              :: NB_INC
      INTEGER, ALLOCATABLE :: BUFFER(:)
      INTEGER              :: NUM_POINTS
      INTEGER, ALLOCATABLE :: POINT_INDICES(:)
      INTEGER              :: PRIM_CELL_IDX
      INTEGER              :: I
      INTEGER              :: NB_LOW
      INTEGER              :: BAND_IDX
      INTEGER              :: GRID_SIZE(3)
      REAL(q)              :: GRID_TRAFO(3, 3)
      REAL(q)              :: R_IN_GRID(3)
      REAL(q)              :: R_IN_BOX(3)
      GDEF, POINTER        :: WAVE_FUNC_PTR(:, :, :, :, :, :, :)
      COMPLEX(q), POINTER  :: CR_PTR(:, :, :)
      INTEGER              :: X, Y, Z
      INTEGER              :: SHIFT_X, SHIFT_Y, SHIFT_Z
      INTEGER              :: R_VEC(3)
      INTEGER              :: BOX_BOUNDS(3, 0 : 1)
      REAL(q)              :: INTERPOL_BUFFER(0 : 1, 0 : 1, 0 : 1)
      INTEGER              :: I1, I2, I3

      INTEGER, PARAMETER   :: NUM_POINTS_ALONG_PC(3) = [20, 20, 20]

      CALL SETWDES(W%WDES, WDES1, 1)

      IF (WDES1%GRID%RL%NFAST /= 1) CALL vtutor%error("GATHER_PRIM_CELL_WAVE_FUNCTIONS: NFAST /= 1 not implemented!")
#ifdef MPI
      IF (WDES1%GRID%COMM%NCPU /= 1) CALL vtutor%error("GATHER_PRIM_CELL_WAVE_FUNCTIONS: NCPU /= 1 not implemented!")
#endif

      NB_INC = WDES1%NB_TOT
      ALLOCATE(BAND_INDICES(NB_INC))
      BAND_INDICES = [(I, I = 1, NB_INC)]

      IF (ALLOCATED(WAVE_FUNC)) DEALLOCATE(WAVE_FUNC)
      ALLOCATE(WAVE_FUNC(PRODUCT(NUM_POINTS_ALONG_PC) * 27, NB_INC))
      WAVE_FUNC = 0

      ! Use rank remapping for easy indexing
      WAVE_FUNC_PTR( &
         0 : 2, 0 : 2, 0 : 2, &
         0 : NUM_POINTS_ALONG_PC(1) - 1, &
         0 : NUM_POINTS_ALONG_PC(2) - 1, &
         0 : NUM_POINTS_ALONG_PC(3) - 1, &
         1 : NB_INC) => WAVE_FUNC

      CALL NEWWAV(W1, WDES1, .true.)

      GRID_SIZE = [WDES1%GRID%NGX, WDES1%GRID%NGY, WDES1%GRID%NGZ]
      GRID_TRAFO = MATMUL(TRANSPOSE(LATT_CUR%B), PRIM_CELL%PRIM_LATT%A)
      GRID_TRAFO = GRID_TRAFO * DYAD(REAL(GRID_SIZE, q), 1._q / NUM_POINTS_ALONG_PC)

      DO I = 1, NB_INC
         ! Check if we're on the correct node.
         NB_LOW = MODULO(BAND_INDICES(I) - 1, WDES1%NB_PAR) + 1
         IF (NB_LOW /= WDES1%NB_LOW) CYCLE
         ! Fetch local band index to index the array.
         BAND_IDX = (BAND_INDICES(I) - 1) / WDES1%NB_PAR + 1

         ! Copy the orbital and FFT to real-space wave function.
         CALL W1_COPY(ELEMENT(W, WDES1, BAND_IDX, 1), W1)
         CALL FFTWAV_W1(W1)

         ! Use rank remapping for easy indexing
         CR_PTR(0 : GRID_SIZE(1) - 1, 0 : GRID_SIZE(2) - 1, 0 : GRID_SIZE(3) - 1) => W1%CR

         ! Loop over real-space grid points
         DO Z = 0, NUM_POINTS_ALONG_PC(3) - 1
            DO Y = 0, NUM_POINTS_ALONG_PC(2) - 1
               DO X = 0, NUM_POINTS_ALONG_PC(1) - 1

                  ! Loop over neighbor cells
                  DO SHIFT_Z = 0, 2
                     DO SHIFT_Y = 0, 2
                        DO SHIFT_X = 0, 2
                           R_VEC = [X, Y, Z] + (MODULO([SHIFT_X, SHIFT_Y, SHIFT_Z] + 1, 3) - 1) * NUM_POINTS_ALONG_PC

                           R_IN_GRID = MODULO(MATMUL(GRID_TRAFO, R_VEC), REAL(GRID_SIZE, q))
                           
                           BOX_BOUNDS(:, 0) = FLOOR(R_IN_GRID)
                           BOX_BOUNDS(:, 1) = MODULO(BOX_BOUNDS(:, 0) + 1, GRID_SIZE)

                           R_IN_BOX = R_IN_GRID - BOX_BOUNDS(:, 0)

                           DO I3 = 0, 1
                              DO I2 = 0, 1
                                 DO I1 = 0, 1
                                    INTERPOL_BUFFER(I1, I2, I3) = REAL(CR_PTR(BOX_BOUNDS(1, I1), BOX_BOUNDS(2, I2), BOX_BOUNDS(3, I3)))
                                 ENDDO
                              ENDDO
                           ENDDO

                           DO I2 = 0, 1
                              DO I1 = 0, 1
                                 INTERPOL_BUFFER(I1, I2, 0) = INTERPOL_BUFFER(I1, I2, 0) * (1._q - R_IN_BOX(1)) + INTERPOL_BUFFER(I1, I2, 1) * R_IN_BOX(1)
                              ENDDO
                           ENDDO

                           DO I1 = 0, 1
                              INTERPOL_BUFFER(I1, 0, 0) = INTERPOL_BUFFER(I1, 0, 0) * (1._q - R_IN_BOX(2)) + INTERPOL_BUFFER(I1, 1, 0) * R_IN_BOX(2)
                           ENDDO

                           WAVE_FUNC_PTR(SHIFT_X, SHIFT_Y, SHIFT_Z, X, Y, Z, I) = INTERPOL_BUFFER(1, 0, 0) * (1._q - R_IN_BOX(3)) + INTERPOL_BUFFER(0, 0, 0) * R_IN_BOX(3)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      CALL DELWAV(W1, .TRUE.)

      CALLMPI(M_sum_g(W%WDES%COMM, WAVE_FUNC(1, 1), SIZE(WAVE_FUNC))) ! _test_mani comm?

      ! do i = 1, size(WAVE_FUNC, 1)
      !    if (modulo(i - 1, 8) /= 0) WAVE_FUNC(i, :) = WAVE_FUNC(i, :) * 0.0000001
      ! enddo

   END SUBROUTINE GATHER_PRIM_CELL_WAVE_FUNCTIONS_3_27


   SUBROUTINE GATHER_PRIM_CELL_WAVE_FUNCTIONS_3(W, PRIM_CELL, LATT_CUR, WAVE_FUNC, BAND_INDICES)
      USE wave_struct_def, ONLY: wavefun1
      USE wave_high,       ONLY: FFTWAV_W1, &
                                 ELEMENT, &
                                 W1_COPY
      USE wave,            ONLY: SETWDES, &
                                 NEWWAV, &
                                 DELWAV
      IMPLICIT NONE

      TYPE(wavespin),            INTENT(IN)  :: W
      TYPE(PRIM_CELL_T),         INTENT(IN)  :: PRIM_CELL
      TYPE(latt),                INTENT(IN)  :: LATT_CUR
      GDEF, ALLOCATABLE, TARGET, INTENT(OUT) :: WAVE_FUNC(:, :)
      INTEGER, ALLOCATABLE,      INTENT(OUT) :: BAND_INDICES(:)

      TYPE(wavedes1)       :: WDES1
      TYPE(wavefun1)       :: W1
      INTEGER              :: NB_INC
      INTEGER, ALLOCATABLE :: BUFFER(:)
      INTEGER              :: NUM_POINTS
      INTEGER, ALLOCATABLE :: POINT_INDICES(:)
      INTEGER              :: PRIM_CELL_IDX
      INTEGER              :: I
      INTEGER              :: NB_LOW
      INTEGER              :: BAND_IDX
      INTEGER              :: GRID_SIZE(3)
      REAL(q)              :: GRID_TRAFO(3, 3)
      REAL(q)              :: R_IN_GRID(3)
      REAL(q)              :: R_IN_BOX(3)
      GDEF, POINTER        :: WAVE_FUNC_PTR(:, :, :, :, :, :, :)
      COMPLEX(q), POINTER  :: CR_PTR(:, :, :)
      INTEGER              :: X, Y, Z
      INTEGER              :: SHIFT_X, SHIFT_Y, SHIFT_Z
      INTEGER              :: R_VEC(3)
      INTEGER              :: BOX_BOUNDS(3, 0 : 1)
      REAL(q)              :: INTERPOL_BUFFER(0 : 1, 0 : 1, 0 : 1)
      INTEGER              :: I1, I2, I3

      INTEGER, PARAMETER   :: NUM_POINTS_ALONG_PC(3) = [20, 20, 20]

      CALL SETWDES(W%WDES, WDES1, 1)

      IF (WDES1%GRID%RL%NFAST /= 1) CALL vtutor%error("GATHER_PRIM_CELL_WAVE_FUNCTIONS: NFAST /= 1 not implemented!")
#ifdef MPI
      IF (WDES1%GRID%COMM%NCPU /= 1) CALL vtutor%error("GATHER_PRIM_CELL_WAVE_FUNCTIONS: NCPU /= 1 not implemented!")
#endif

      NB_INC = WDES1%NB_TOT
      ALLOCATE(BAND_INDICES(NB_INC))
      BAND_INDICES = [(I, I = 1, NB_INC)]

      IF (ALLOCATED(WAVE_FUNC)) DEALLOCATE(WAVE_FUNC)
      ALLOCATE(WAVE_FUNC(PRODUCT(NUM_POINTS_ALONG_PC) * 8, NB_INC))
      WAVE_FUNC = 0

      ! Use rank remapping for easy indexing
      WAVE_FUNC_PTR( &
         0 : 1, 0 : 1, 0 : 1, &
         0 : NUM_POINTS_ALONG_PC(1) - 1, &
         0 : NUM_POINTS_ALONG_PC(2) - 1, &
         0 : NUM_POINTS_ALONG_PC(3) - 1, &
         1 : NB_INC) => WAVE_FUNC

      CALL NEWWAV(W1, WDES1, .true.)

      GRID_SIZE = [WDES1%GRID%NGX, WDES1%GRID%NGY, WDES1%GRID%NGZ]
      GRID_TRAFO = MATMUL(TRANSPOSE(LATT_CUR%B), PRIM_CELL%PRIM_LATT%A)
      GRID_TRAFO = GRID_TRAFO * DYAD(REAL(GRID_SIZE, q), 1._q / NUM_POINTS_ALONG_PC)

      DO I = 1, NB_INC
         ! Check if we're on the correct node.
         NB_LOW = MODULO(BAND_INDICES(I) - 1, WDES1%NB_PAR) + 1
         IF (NB_LOW /= WDES1%NB_LOW) CYCLE
         ! Fetch local band index to index the array.
         BAND_IDX = (BAND_INDICES(I) - 1) / WDES1%NB_PAR + 1

         ! Copy the orbital and FFT to real-space wave function.
         CALL W1_COPY(ELEMENT(W, WDES1, BAND_IDX, 1), W1)
         CALL FFTWAV_W1(W1)

         ! Use rank remapping for easy indexing
         CR_PTR(0 : GRID_SIZE(1) - 1, 0 : GRID_SIZE(2) - 1, 0 : GRID_SIZE(3) - 1) => W1%CR

         ! Loop over real-space grid points
         DO Z = 0, NUM_POINTS_ALONG_PC(3) - 1
            DO Y = 0, NUM_POINTS_ALONG_PC(2) - 1
               DO X = 0, NUM_POINTS_ALONG_PC(1) - 1

                  ! Loop over neighbor cells
                  DO SHIFT_Z = 0, 1
                     DO SHIFT_Y = 0, 1
                        DO SHIFT_X = 0, 1
                           R_VEC = [X, Y, Z] + [SHIFT_X, SHIFT_Y, SHIFT_Z] * NUM_POINTS_ALONG_PC
                           R_VEC = MODULO(R_VEC + NUM_POINTS_ALONG_PC / 2, 2 * NUM_POINTS_ALONG_PC) - NUM_POINTS_ALONG_PC / 2

                           R_IN_GRID = MODULO(MATMUL(GRID_TRAFO, R_VEC), REAL(GRID_SIZE, q))
                           
                           BOX_BOUNDS(:, 0) = FLOOR(R_IN_GRID)
                           BOX_BOUNDS(:, 1) = MODULO(BOX_BOUNDS(:, 0) + 1, GRID_SIZE)

                           R_IN_BOX = R_IN_GRID - BOX_BOUNDS(:, 0)

                           DO I3 = 0, 1
                              DO I2 = 0, 1
                                 DO I1 = 0, 1
                                    INTERPOL_BUFFER(I1, I2, I3) = REAL(CR_PTR(BOX_BOUNDS(1, I1), BOX_BOUNDS(2, I2), BOX_BOUNDS(3, I3)))
                                 ENDDO
                              ENDDO
                           ENDDO

                           DO I2 = 0, 1
                              DO I1 = 0, 1
                                 INTERPOL_BUFFER(I1, I2, 0) = INTERPOL_BUFFER(I1, I2, 0) * (1._q - R_IN_BOX(1)) + INTERPOL_BUFFER(I1, I2, 1) * R_IN_BOX(1)
                              ENDDO
                           ENDDO

                           DO I1 = 0, 1
                              INTERPOL_BUFFER(I1, 0, 0) = INTERPOL_BUFFER(I1, 0, 0) * (1._q - R_IN_BOX(2)) + INTERPOL_BUFFER(I1, 1, 0) * R_IN_BOX(2)
                           ENDDO

                           WAVE_FUNC_PTR(SHIFT_X, SHIFT_Y, SHIFT_Z, X, Y, Z, I) = INTERPOL_BUFFER(1, 0, 0) * (1._q - R_IN_BOX(3)) + INTERPOL_BUFFER(0, 0, 0) * R_IN_BOX(3)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      CALL DELWAV(W1, .TRUE.)

      CALLMPI(M_sum_g(W%WDES%COMM, WAVE_FUNC(1, 1), SIZE(WAVE_FUNC))) ! _test_mani comm?

      ! do i = 1, size(WAVE_FUNC, 1)
      !    if (modulo(i - 1, 8) /= 0) WAVE_FUNC(i, :) = WAVE_FUNC(i, :) * 0.0000001
      ! enddo

   END SUBROUTINE GATHER_PRIM_CELL_WAVE_FUNCTIONS_3


   ! _test_mani doc: shift in direct coordinates
   SUBROUTINE SHIFT_ORBITAL(W1, SHIFT)
      USE wave_struct_def, ONLY: wavefun1
      USE CONSTANT
      IMPLICIT NONE

      TYPE(wavefun1)       :: W1 ! No intent(inout) so we can call with ELEMENT()
      REAL(q), INTENT(IN)  :: SHIFT(3)

      INTEGER  :: G_IDX
      REAL(q)  :: GK_DIR(3)
      
      DO G_IDX = 1, W1%WDES1%NGVECTOR
         GK_DIR = [REAL(q):: W1%WDES1%IGX(G_IDX), W1%WDES1%IGY(G_IDX), W1%WDES1%IGZ(G_IDX)] + W1%WDES1%VKPT
         W1%CPTWFP(G_IDX) = W1%CPTWFP(G_IDX) * EXP(-CITPI * DOT_PRODUCT(GK_DIR, SHIFT))
      ENDDO
   END SUBROUTINE SHIFT_ORBITAL


   ! I don't trust the compiler to do deallocations in a smart way, so I do all the deallocations manually as soon as possible.
   SUBROUTINE CALC_WANNIER_ORBITALS_SCDM(W, NUM_ORBS_IN_CELL, W_WANNIER, T, MU, SIGMA, PRIM_CELL, NONL_S, NONLR_S, CQIJ, IO)
      USE wave_struct_def,    ONLY: wavefun1
      USE nonl_struct_def,    ONLY: nonl_struct
      USE nonlr_struct_def,   ONLY: nonlr_struct
      USE nonl_high,          ONLY: PROALL
      USE wave,               ONLY: SETWDES, &
                                    ALLOCW, &
                                    DEALLOCW
      USE wave_high,          ONLY: ELEMENT, &
                                    ELEMENTS
      USE locproj,            ONLY: LPRJ_LINCOM
      USE mathtools,          ONLY: SV_DECOMPOSER
      IMPLICIT NONE

      TYPE(wavespin),      INTENT(IN)     :: W ! _test_mani
      INTEGER,             INTENT(INOUT)  :: NUM_ORBS_IN_CELL
      TYPE(wavespin),      INTENT(OUT)    :: W_WANNIER
      GDEF, ALLOCATABLE,   INTENT(OUT)    :: T(:, :, :, :)
      REAL(q),             INTENT(IN)     :: MU
      REAL(q),             INTENT(IN)     :: SIGMA
      TYPE(PRIM_CELL_T),   INTENT(IN)     :: PRIM_CELL
      TYPE(nonl_struct),   INTENT(IN)     :: NONL_S
      TYPE(nonlr_struct),  INTENT(IN)     :: NONLR_S
      OVERLAP,             INTENT(IN)     :: CQIJ(:, :, :, :)
      TYPE(in_struct),     INTENT(IN)     :: IO
      
      TYPE(wavedes1)       :: WDES1
      TYPE(wavefun1)       :: W1
      GDEF, ALLOCATABLE    :: PSI(:, :)
      GDEF, ALLOCATABLE    :: PSI_D(:, :)
      GDEF, ALLOCATABLE    :: SOVL(:, :)
      INTEGER, ALLOCATABLE :: BAND_INDICES(:)
      INTEGER              :: BAND_IDX
      INTEGER              :: BAND_IDX_GLOBAL
      INTEGER              :: NUM_GRID_POINTS
      INTEGER              :: NB_TOT
      INTEGER              :: NB_INC
      INTEGER              :: NUM_WAN
      INTEGER              :: CELL_IDX
      INTEGER              :: I, J
      REAL(q)              :: SHIFT(3)
      REAL(q)              :: CUTOFF_FACTOR

      INTEGER, ALLOCATABLE :: PIVOT(:)
      INTEGER              :: INFO

      TYPE(SV_DECOMPOSER)  :: SVD

      CALL SETWDES(W%WDES, WDES1, 1)

      NB_TOT = WDES1%NB_TOT

      ! _test_mani TEST
      CALL GATHER_PRIM_CELL_WAVE_FUNCTIONS_3(W, PRIM_CELL, PRIM_CELL%SUPER_LATT, PSI, BAND_INDICES)
      NUM_GRID_POINTS = SIZE(PSI, 1)
      NB_INC = SIZE(BAND_INDICES)
      
      ! Entanglement
      DO I = 1, NB_INC
         ! Gauss
         ! CUTOFF_FACTOR = (REAL(W%CELTOT(BAND_INDICES(I), 1, 1)) - SCDM_MU) / SCDM_SIGMA
         ! CUTOFF_FACTOR = EXP(-CUTOFF_FACTOR * CUTOFF_FACTOR)
         ! erfc
         CUTOFF_FACTOR = .5_q * ERFC((REAL(W%CELTOT(BAND_INDICES(I), 1, 1)) - MU) / SIGMA)

         ! IF (CUTOFF_FACTOR < 1e-10_q)

         PSI(:, I) = PSI(:, I) / SQRT(ABS(DOT_PRODUCT(PSI(:, I), PSI(:, I))))
         PSI(:, I) = PSI(:, I) * CUTOFF_FACTOR
      ENDDO

      ALLOCATE(PSI_D(NB_INC, NUM_GRID_POINTS))
      PSI_D = GCONJG(TRANSPOSE(PSI))
      ! DEALLOCATE(PSI) ! Save memory. _test_mani

      ! CALL CUSTOM_QR(PSI_D, PIVOT)
      CALL CUSTOM_QR_IMAGES(PSI_D, PIVOT, NUM_ORBS_IN_CELL, 8)

      IF (SIZE(PIVOT) == 0) CALL vtutor%error("CALC_WANNIER_ORBITALS_SCDM: " // &
         "The SCDM cutoff supresses all orbitals, please increase the parameters SCDM_MU and/or SCDM_SIGMA.")

      IF (SIZE(PIVOT) /= NUM_ORBS_IN_CELL) CALL vtutor%alert( &
         "The SCDM cutoff supresses some orbitals. New number of generated Wannier orbitals = " // &
         str(SIZE(PIVOT)) // ". To avoid this, increase the parameters SCDM_MU and/or SCDM_SIGMA.")

      ! Adapt the number of orbitals to pivot size
      NUM_ORBS_IN_CELL = SIZE(PIVOT)
      NUM_WAN = NUM_ORBS_IN_CELL * PRIM_CELL%NUM_CELLS

      PIVOT = (PIVOT - 1) * 8 + 1
      DEALLOCATE(PSI_D)
      ALLOCATE(PSI_D(NB_INC, NUM_ORBS_IN_CELL))

      ! Use pivot, then orthonormalization via SVD. _test_mani
      PSI_D(:, :) = GCONJG(TRANSPOSE(PSI(PIVOT(:), :)))
      CALL SVD%INIT('A', NB_INC, NUM_ORBS_IN_CELL, INFO) ! _test_mani handle errors
      CALL SVD%CALC(PSI_D, INFO) ! _test_mani handle errors

      ! Retrieve final transformation matrix.
      PSI_D = MATMUL(SVD%U_MAT(:, 1 : NUM_ORBS_IN_CELL), SVD%V_DAGGER_MAT) ! _test_mani
      CALL SVD%DEALLOC() ! Save memory

      ALLOCATE(T(NB_TOT, NUM_WAN, 1, 1))
      T = 0
      DO CELL_IDX = 1, PRIM_CELL%NUM_CELLS
         BAND_IDX = (CELL_IDX - 1) * NUM_ORBS_IN_CELL + 1
         T(BAND_INDICES(:), BAND_IDX : BAND_IDX + NUM_ORBS_IN_CELL - 1 , 1, 1) = PSI_D(:, 1 : NUM_ORBS_IN_CELL) ! _test_mani make pivot correct size
      ENDDO
      DEALLOCATE(PSI_D) ! Save memory.

      ! _test_mani lincom
      CALL DEALLOCW(W_WANNIER)
      CALL ALLOCW(W%WDES, W_WANNIER)
      CALL LPRJ_LINCOM(W, T, W_WANNIER)
      DEALLOCATE(T) ! Save memory.

      CALL POPULATE_SUPERCELL_WITH_ORBITALS(ELEMENTS(W_WANNIER, WDES1, 1), NUM_ORBS_IN_CELL, PRIM_CELL)
      CALL PROALL(W_WANNIER%WDES%GRID, PRIM_CELL%SUPER_LATT, NONLR_S, NONL_S, W_WANNIER)
      CALL CALC_ORTHONORMAL_TRAFO(T, W, W_WANNIER, NUM_WAN, CQIJ, IO, BAND_INDICES)

   END SUBROUTINE CALC_WANNIER_ORBITALS_SCDM


   SUBROUTINE CALC_WANNIER_ORBITALS_FROM_PRIM(W, NUM_ORBS_IN_CELL, W_WANNIER, T, PRIM_CELL, NONL_S, NONLR_S, CQIJ, IO)
      USE nonl_struct_def,    ONLY: nonl_struct
      USE nonlr_struct_def,   ONLY: nonlr_struct
      USE nonl_high,          ONLY: PROALL
      USE wave,               ONLY: SETWDES
      USE wave_high,          ONLY: ELEMENTS
      IMPLICIT NONE

      TYPE(wavespin),      INTENT(IN)  :: W
      INTEGER,             INTENT(IN)  :: NUM_ORBS_IN_CELL
      TYPE(wavespin),      INTENT(OUT) :: W_WANNIER
      GDEF, ALLOCATABLE,   INTENT(OUT) :: T(:, :, :, :)
      TYPE(PRIM_CELL_T),   INTENT(IN)  :: PRIM_CELL
      TYPE(nonl_struct),   INTENT(IN)  :: NONL_S
      TYPE(nonlr_struct),  INTENT(IN)  :: NONLR_S
      OVERLAP,             INTENT(IN)  :: CQIJ(:, :, :, :)
      TYPE(in_struct),     INTENT(IN)  :: IO
      
      TYPE(wavedes1)       :: WDES1
      INTEGER              :: NUM_WAN

      CALL SETWDES(W%WDES, WDES1, 1)
      NUM_WAN = NUM_ORBS_IN_CELL * PRIM_CELL%NUM_CELLS

      CALL IMPORT_WANNIER_FROM_PRIM(W, W_WANNIER, NUM_ORBS_IN_CELL, PRIM_CELL, PRIM_CELL%SUPER_LATT, NONL_S, NONLR_S, CQIJ, IO)
      CALL POPULATE_SUPERCELL_WITH_ORBITALS(ELEMENTS(W_WANNIER, WDES1, 1), NUM_ORBS_IN_CELL, PRIM_CELL)
      CALL PROALL(W_WANNIER%WDES%GRID, PRIM_CELL%SUPER_LATT, NONLR_S, NONL_S, W_WANNIER)
      CALL CALC_ORTHONORMAL_TRAFO(T, W, W_WANNIER, NUM_WAN, CQIJ, IO)
   END SUBROUTINE CALC_WANNIER_ORBITALS_FROM_PRIM


   SUBROUTINE CALC_ORTHONORMAL_TRAFO(T, W, W_WANNIER, NUM_WAN, CQIJ, IO, BAND_INDICES)
      USE wave,      ONLY: SETWDES
      USE wave_high, ONLY: ELEMENTS
      USE choleski,  ONLY: CALC_PAW_OVERLAP
      USE locproj,   ONLY: LPRJ_LINCOM
      USE mathtools, ONLY: SV_DECOMPOSER
      IMPLICIT NONE

      GDEF, ALLOCATABLE,   INTENT(OUT)    :: T(:, :, :, :)
      TYPE(wavespin),      INTENT(IN)     :: W
      TYPE(wavespin),      INTENT(INOUT)  :: W_WANNIER
      INTEGER,             INTENT(IN)     :: NUM_WAN
      OVERLAP,             INTENT(IN)     :: CQIJ(:, :, :, :)
      TYPE(in_struct),     INTENT(IN)     :: IO
      INTEGER, OPTIONAL,   INTENT(IN)     :: BAND_INDICES(:)

      TYPE(wavedes1)       :: WDES1
      INTEGER              :: NB_TOT
      INTEGER              :: K_IDX
      INTEGER              :: SPIN_IDX
      GDEF, ALLOCATABLE    :: SOVL(:, :)
      TYPE(SV_DECOMPOSER)  :: SVD
      INTEGER              :: INFO

      NB_TOT = W%WDES%NB_TOT

      ALLOCATE(SOVL(NB_TOT, NB_TOT))
      ALLOCATE(T(NB_TOT, NUM_WAN, W%WDES%NKPTS, W%WDES%ISPIN))
      T = 0

      CALL SVD%INIT('A', NB_TOT, NUM_WAN, INFO) ! _test_mani handle errors

      DO SPIN_IDX = 1, W%WDES%ISPIN
         DO K_IDX = 1, W%WDES%NKPTS
            CALL SETWDES(W%WDES, WDES1, K_IDX)

            CALL CALC_PAW_OVERLAP(ELEMENTS(W, WDES1, SPIN_IDX), ELEMENTS(W_WANNIER, WDES1, SPIN_IDX), CQIJ, SOVL)

            IF (PRESENT(BAND_INDICES)) THEN
               T(BAND_INDICES(:), :, K_IDX, SPIN_IDX) = SOVL(BAND_INDICES(:), 1 : NUM_WAN)
            ELSE
               T(:, :, K_IDX, SPIN_IDX) = SOVL(:, 1 : NUM_WAN)
            ENDIF

            ! Orthonormalization via SVD
            CALL SVD%CALC(T(:, :, 1, 1), INFO) ! _test_mani handle errors

            ! Retrieve final transformation matrix
            T(:, :, K_IDX, SPIN_IDX) = MATMUL(SVD%U_MAT(:, 1 : NUM_WAN), SVD%V_DAGGER_MAT)
         ENDDO
      ENDDO

      ! Perform transformation to receive orthonormal Wannier orbitals
      CALL LPRJ_LINCOM(W, T, W_WANNIER)
   END SUBROUTINE CALC_ORTHONORMAL_TRAFO


   SUBROUTINE POPULATE_SUPERCELL_WITH_ORBITALS(WA, NUM_ORBS_IN_CELL, PRIM_CELL)
      USE wave_struct_def, ONLY: wavefuna
      USE wave_high,       ONLY: ELEMENT
      IMPLICIT NONE

      TYPE(wavefuna)                      :: WA ! No intent(inout) so we can call with ELEMENT()
      INTEGER,             INTENT(IN)     :: NUM_ORBS_IN_CELL
      TYPE(PRIM_CELL_T),   INTENT(IN)     :: PRIM_CELL

      INTEGER        :: BAND_IDX
      INTEGER        :: BAND_IDX_GLOBAL
      INTEGER        :: NUM_ORBS_TOT
      INTEGER        :: CELL_IDX
      REAL(q)        :: SHIFT(3)

      CALL DUPLICATE_PRIM_ORBITALS(WA, NUM_ORBS_IN_CELL, PRIM_CELL)

      NUM_ORBS_TOT = NUM_ORBS_IN_CELL * PRIM_CELL%NUM_CELLS

      DO BAND_IDX = 1, WA%WDES1%NBANDS
         BAND_IDX_GLOBAL = (BAND_IDX - 1) * WA%WDES1%NB_PAR + WA%WDES1%NB_LOW
         IF (BAND_IDX_GLOBAL > NUM_ORBS_TOT) CYCLE

         ! Target cell is deduced from global band index (MPI)
         CELL_IDX = (BAND_IDX_GLOBAL - 1) / NUM_ORBS_IN_CELL + 1
         ! Calculate shift vector in direct coordinates
         SHIFT = PRIM_CELL%GET_POSION(CELL_IDX, 1) - PRIM_CELL%GET_POSION(1, 1)

         ! Shift orbital to target cell
         CALL SHIFT_ORBITAL(ELEMENT(WA, BAND_IDX), SHIFT)
      ENDDO
   END SUBROUTINE POPULATE_SUPERCELL_WITH_ORBITALS


   SUBROUTINE CUSTOM_QR_IMAGES(A, PIVOT, NUM_BLOCKS, BLOCK_SIZE)
      IMPLICIT NONE

      GDEF,                   INTENT(INOUT)  :: A(:, :)
      INTEGER, ALLOCATABLE,   INTENT(OUT)    :: PIVOT(:)
      INTEGER,                INTENT(IN)     :: NUM_BLOCKS
      INTEGER,                INTENT(IN)     :: BLOCK_SIZE

      INTEGER              :: NUM_BLOCKS_TOT
      INTEGER              :: BLOCK_IDX
      INTEGER              :: BLOCK_OFFSET
      INTEGER              :: CORNER_IDX
      INTEGER              :: NUM_ROWS, NUM_COLS
      INTEGER              :: COL_IDX, COL_IDX_2
      INTEGER              :: MAX_BLOCK
      INTEGER              :: MAX_COL
      REAL(q)              :: MAX_VAL, CUR_VAL
      REAL(q)              :: ALPHA
      REAL(q)              :: TAU
      INTEGER              :: TMP
      INTEGER, ALLOCATABLE :: TMP_PIVOT(:)
      GDEF, ALLOCATABLE    :: BUFFER(:)
      GDEF, ALLOCATABLE    :: Q_MAT(:, :)

      ! _test_mani check some stuff (num_blocks, etc.)

      NUM_ROWS = SIZE(A, 1)
      NUM_COLS = SIZE(A, 2)

      NUM_BLOCKS_TOT = NUM_COLS / BLOCK_SIZE

      ALLOCATE(PIVOT(NUM_BLOCKS_TOT))
      PIVOT = [(BLOCK_IDX, BLOCK_IDX = 1, NUM_BLOCKS_TOT)]

      ALLOCATE(BUFFER(NUM_COLS))
      ALLOCATE(Q_MAT(NUM_ROWS, NUM_ROWS))

      DO BLOCK_IDX = 1, NUM_BLOCKS
         CORNER_IDX = (BLOCK_IDX - 1) * BLOCK_SIZE + 1

         MAX_VAL = 0._q
         MAX_COL = 0

         DO COL_IDX = CORNER_IDX, NUM_COLS, BLOCK_SIZE
            CUR_VAL = ABS(DOT_PRODUCT(A(CORNER_IDX : NUM_ROWS, COL_IDX), A(CORNER_IDX : NUM_ROWS, COL_IDX)))

            IF (CUR_VAL > MAX_VAL) THEN
               MAX_VAL = CUR_VAL
               MAX_COL = COL_IDX
            ENDIF
         ENDDO

         ! If the maximum column is all zeros, we return early
         IF (MAX_COL == 0) EXIT

         MAX_BLOCK = (MAX_COL - 1) / BLOCK_SIZE + 1

         ! Swap pivot entries.
         TMP = PIVOT(BLOCK_IDX)
         PIVOT(BLOCK_IDX) = PIVOT(MAX_BLOCK)
         PIVOT(MAX_BLOCK) = TMP

         ! Iterate over all columns inside the block.
         DO BLOCK_OFFSET = 0, BLOCK_SIZE - 1
            COL_IDX = CORNER_IDX + BLOCK_OFFSET
            
            ! Swap columns. _test_mani dswap/zswap
            BUFFER(1 : NUM_ROWS) = A(:, COL_IDX)
            A(:, COL_IDX) = A(:, MAX_COL + BLOCK_OFFSET)
            A(:, MAX_COL + BLOCK_OFFSET) = BUFFER(1 : NUM_ROWS)

! _test_mani implement general case
#ifdef gammareal
            IF (COL_IDX < NUM_COLS) THEN
               CALL DLARFG(NUM_ROWS - COL_IDX + 1, A(COL_IDX, COL_IDX), A(COL_IDX + 1, COL_IDX), 1, TAU)
            ELSE
               CALL DLARFG(1, A(NUM_COLS, NUM_COLS), A(NUM_COLS, NUM_COLS), 1, TAU)
            ENDIF
            ALPHA = A(COL_IDX, COL_IDX)
            A(COL_IDX, COL_IDX) = 1._q
            CALL DLARF("L", NUM_ROWS - COL_IDX + 1, NUM_COLS - COL_IDX, A(COL_IDX, COL_IDX), 1, TAU, A(COL_IDX, COL_IDX + 1), NUM_ROWS, BUFFER)
            A(COL_IDX, COL_IDX) = ALPHA
#endif

            ! ! Perform Householder transformation.
            ! ! ALPHA = -SIGN(NORM2(A(COL_IDX : NUM_ROWS, COL_IDX)), A(COL_IDX, COL_IDX)) ! _test_mani
            ! ALPHA = NORM2(A(COL_IDX : NUM_ROWS, COL_IDX))
            ! BUFFER(1 : COL_IDX - 1) = 0
            ! BUFFER(COL_IDX) = A(COL_IDX, COL_IDX) - ALPHA
            ! BUFFER(COL_IDX + 1 : NUM_ROWS) = A(COL_IDX + 1 : NUM_ROWS, COL_IDX)

            ! BUFFER = BUFFER / NORM2(BUFFER)

            ! ! DO COL_IDX_2 = COL_IDX, NUM_COLS
            ! !    A(COL_IDX : NUM_ROWS, COL_IDX_2) = A(COL_IDX : NUM_ROWS, COL_IDX_2) - BUFFER(COL_IDX : NUM_ROWS) * ABS(DOT_PRODUCT(A(COL_IDX : NUM_ROWS, COL_IDX_2), BUFFER(COL_IDX : NUM_ROWS)))
            ! ! ENDDO

            ! Q_MAT = 0._q
            ! DO COL_IDX = 1, NUM_ROWS
            !    Q_MAT(COL_IDX, COL_IDX) = 1._q
            ! ENDDO
            ! Q_MAT = Q_MAT - 2 * DYAD(BUFFER, BUFFER)

            ! ! A = MATMUL(Q_MAT, A)

            ! DO COL_IDX = 1, NUM_COLS
            !    CALL GGEMV('N', NUM_ROWS, NUM_ROWS, one, Q_MAT, NUM_ROWS, A(:, COL_IDX), 1, zero, BUFFER, 1)
            !    A(:, COL_IDX) = BUFFER
            ! ENDDO
         ENDDO
      ENDDO

      ALLOCATE(TMP_PIVOT(BLOCK_IDX - 1))
      TMP_PIVOT = PIVOT(1 : BLOCK_IDX - 1)
      CALL MOVE_ALLOC(TMP_PIVOT, PIVOT)
   END SUBROUTINE CUSTOM_QR_IMAGES


   SUBROUTINE CUSTOM_QR(A, PIVOT)
      IMPLICIT NONE

      GDEF,                   INTENT(INOUT)  :: A(:, :)
      INTEGER, ALLOCATABLE,   INTENT(OUT)    :: PIVOT(:)

      INTEGER           :: NUM_ROWS, NUM_COLS
      INTEGER           :: ROW_IDX, COL_IDX
      INTEGER           :: MAX_COL
      REAL(q)           :: MAX_VAL, CUR_VAL
      REAL(q)           :: ALPHA
      INTEGER           :: TMP
      GDEF, ALLOCATABLE :: BUFFER(:)
      GDEF, ALLOCATABLE :: Q_MAT(:, :)

      NUM_ROWS = SIZE(A, 1)
      NUM_COLS = SIZE(A, 2)

      ALLOCATE(PIVOT(NUM_COLS))
      PIVOT = [(COL_IDX, COL_IDX = 1, NUM_COLS)]

      ALLOCATE(BUFFER(NUM_ROWS))
      ALLOCATE(Q_MAT(NUM_ROWS, NUM_ROWS))

      DO ROW_IDX = 1, NUM_ROWS
         MAX_VAL = 0._q
         MAX_COL = 0

         DO COL_IDX = ROW_IDX, NUM_COLS
            CUR_VAL = ABS(DOT_PRODUCT(A(ROW_IDX : NUM_ROWS, COL_IDX), A(ROW_IDX : NUM_ROWS, COL_IDX)))

            IF (CUR_VAL > MAX_VAL) THEN
               MAX_VAL = CUR_VAL
               MAX_COL = COL_IDX
            ENDIF
         ENDDO

         ! Swap pivot entries.
         TMP = PIVOT(ROW_IDX)
         PIVOT(ROW_IDX) = PIVOT(MAX_COL)
         PIVOT(MAX_COL) = TMP

         ! Swap columns.
         BUFFER = A(:, ROW_IDX)
         A(:, ROW_IDX) = A(:, MAX_COL)
         A(:, MAX_COL) = BUFFER

! _test_mani imlement general case
#ifdef gammareal
         ALPHA = -SIGN(NORM2(A(ROW_IDX : NUM_ROWS, ROW_IDX)), A(ROW_IDX, ROW_IDX))
         BUFFER(1 : ROW_IDX - 1) = 0
         BUFFER(ROW_IDX) = A(ROW_IDX, ROW_IDX) - ALPHA
         BUFFER(ROW_IDX + 1 : NUM_ROWS) = A(ROW_IDX + 1 : NUM_ROWS, ROW_IDX)

         BUFFER = BUFFER / NORM2(BUFFER)
         
         Q_MAT = 0._q
         DO COL_IDX = 1, NUM_ROWS
            Q_MAT(COL_IDX, COL_IDX) = 1._q
         ENDDO
         Q_MAT = Q_MAT - 2 * DYAD(BUFFER, BUFFER)

         A = MATMUL(Q_MAT, A)
#endif
      ENDDO
   END SUBROUTINE CUSTOM_QR


   SUBROUTINE WRITE_PRIM_WANN_FUNCS(W_WANN, DESC, T_INFO, IO)
      USE wave_struct_def,    ONLY: wavefun1
      USE wave,               ONLY: SETWDES, &
                                    NEWWAV, &
                                    DELWAV
      USE wave_high,          ONLY: ELEMENT, &
                                    W1_COPY, &
                                    FFTWAV_W1
      IMPLICIT NONE

      TYPE(wavespin),               INTENT(IN)  :: W_WANN
      TYPE(WANNIER_DESCRIPTOR_T),   INTENT(IN)  :: DESC
      TYPE(type_info),              INTENT(IN)  :: T_INFO
      TYPE(in_struct),              INTENT(IN)  :: IO

      TYPE(wavefun1)       :: W1
      TYPE(wavedes1)       :: WDES1
      INTEGER              :: UNIT
      INTEGER              :: BAND_IDX, BAND_IDX_GLOBAL
      INTEGER              :: NB_LOW
      COMPLEX(q), POINTER  :: CR_PTR(:, :, :)

      CALL SETWDES(W_WANN%WDES, WDES1, 1)
      CALL NEWWAV(W1, WDES1, .true.)

      ! DO BAND_IDX_GLOBAL = 1, WANNIER_DESCRIPTOR%NUM_PRIM_ORBITALS
      DO BAND_IDX_GLOBAL = DESC%NUM_PRIM_ORBITALS * 42 + 1, DESC%NUM_PRIM_ORBITALS * 43
         ! Check if we're on the correct node.
         NB_LOW = MODULO(BAND_IDX_GLOBAL - 1, WDES1%NB_PAR) + 1
         IF (NB_LOW /= WDES1%NB_LOW) CYCLE
         ! Fetch local band index to index the array.
         BAND_IDX = (BAND_IDX_GLOBAL - 1) / WDES1%NB_PAR + 1

         ! Copy the orbital and FFT to real-space wave function.
         CALL W1_COPY(ELEMENT(W_WANN, WDES1, BAND_IDX, 1), W1)
         CALL FFTWAV_W1(W1)

         ! Use rank remapping to represent 3D data
         CR_PTR(0 : WDES1%GRID%NGX - 1, 0 : WDES1%GRID%NGY - 1, 0 : WDES1%GRID%NGZ - 1) => W1%CR

         OPEN(NEWUNIT = UNIT, FILE = 'wannier90_' // str(BAND_IDX_GLOBAL - DESC%NUM_PRIM_ORBITALS * 42, '(I0.5)') // '.xsf', STATUS = 'replace')
         CALL WRITE_GRID_DATA_XSF(REAL(CR_PTR), DESC%PRIM_CELL, T_INFO, UNIT)
         CLOSE(UNIT)
      ENDDO

      CALL DELWAV(W1, .TRUE.)

      IF (IO%IU0 >= 0) THEN
         OPEN(NEWUNIT = UNIT, FILE = 'wannier90_centres.xyz', STATUS = 'replace')
         WRITE(UNIT, '(I0)') DESC%NUM_PRIM_ORBITALS
         WRITE(UNIT, '(A)') 'Written by VASP for testing purposes'
         WRITE(UNIT, '("X ", 3F12.6)') MATMUL(DESC%PRIM_CELL%SUPER_LATT%A, DESC%CENTERS)
         CLOSE(UNIT)
      ENDIF

      CALLMPI(M_barrier(W_WANN%WDES%COMM))

   END SUBROUTINE WRITE_PRIM_WANN_FUNCS


   SUBROUTINE WRITE_GRID_DATA_XSF(DATA, PRIM_CELL, T_INFO, UNIT)
      IMPLICIT NONE

      REAL(q),             INTENT(IN)  :: DATA(:, :, :)
      TYPE(PRIM_CELL_T),   INTENT(IN)  :: PRIM_CELL
      TYPE(type_info),     INTENT(IN)  :: T_INFO
      INTEGER,             INTENT(IN)  :: UNIT

      INTEGER  :: DATA_SHAPE(3)
      INTEGER  :: I

      DATA_SHAPE = SHAPE(DATA)

      WRITE(UNIT, '(A)') 'CRYSTAL', 'PRIMVEC'
      WRITE(UNIT, '(3F12.6)') PRIM_CELL%PRIM_LATT%A
      WRITE(UNIT, '(A)') 'CONVVEC'
      WRITE(UNIT, '(3F12.6)') PRIM_CELL%PRIM_LATT%A
      WRITE(UNIT, '(A)') 'PRIMCOORD'
      WRITE(UNIT, '(I0, " ", I0)') PRIM_CELL%NUM_PRIM_ATOMS, 1
      WRITE(UNIT, '(A6, " ", 3F12.6)') (T_INFO%TYPE(T_INFO%ITYP(PRIM_CELL%SUPER_IDX(1, I))), MATMUL(PRIM_CELL%SUPER_LATT%A, PRIM_CELL%GET_POSION(1, I)), I = 1, PRIM_CELL%NUM_PRIM_ATOMS)

      WRITE(UNIT, '(A)') 'BEGIN_BLOCK_DATAGRID_3D', '3D_field', 'BEGIN_DATAGRID_3D_UNKNOWN'
      WRITE(UNIT, '(3I8)') DATA_SHAPE
      ! WRITE(UNIT, '(3F12.6)') 0., 0., 0.
      WRITE(UNIT, '(3F12.6)') MATMUL(PRIM_CELL%SUPER_LATT%A, PRIM_CELL%GET_POSION(1, 1) - PRIM_CELL%GET_POSION(43, 1))
      WRITE(UNIT, '(3F12.6)') (PRIM_CELL%SUPER_LATT%A(:, I) / DATA_SHAPE(I) * (DATA_SHAPE(I) - 1), I = 1, 3)

      WRITE(UNIT, '(6E17.6E3)') DATA

      WRITE(UNIT, '(A)') 'END_DATAGRID_3D', 'END_BLOCK_DATAGRID_3D'
   END SUBROUTINE WRITE_GRID_DATA_XSF


   SUBROUTINE CALC_WANNIER_HAM_FROM_TRAFO(W_BLOCH, UK, HAM_WANNIER)
      IMPLICIT NONE

      TYPE (wavespin),     INTENT(IN)  :: W_BLOCH
      GDEF,                INTENT(IN)  :: UK(:, :, :, :)
      GDEF, ALLOCATABLE,   INTENT(OUT) :: HAM_WANNIER(:, :, :, :) ! Compatible with array returned by EDDIAG. _test_mani maybe use rank retargeting and rank 3 array?

      INTEGER  :: NB_TOT, NUM_WAN
      INTEGER  :: K_IDX, SPIN_IDX
      INTEGER  :: WAN_IDX_1, WAN_IDX_2

      NB_TOT = SIZE(UK, 1)
      NUM_WAN = SIZE(UK, 2)

      ! Setup Wannier-space Hamiltonian matrix.
      ALLOCATE(HAM_WANNIER(NUM_WAN, NUM_WAN, 1, W_BLOCH%WDES%ISPIN))
      HAM_WANNIER = 0

      ! Loop over all spin and k components.
      DO SPIN_IDX = 1, W_BLOCH%WDES%ISPIN
         DO K_IDX = 1, W_BLOCH%WDES%NKPTS

            ! Calculate the Hamiltonian matrix using a simple basis transform.
            DO WAN_IDX_2 = 1, NUM_WAN
               DO WAN_IDX_1 = 1, NUM_WAN
                  HAM_WANNIER(WAN_IDX_1, WAN_IDX_2, 1, SPIN_IDX) = &
                  HAM_WANNIER(WAN_IDX_1, WAN_IDX_2, 1, SPIN_IDX) + SUM( &
                  GCONJG(UK(:, WAN_IDX_1, K_IDX, SPIN_IDX)) * &
                  W_BLOCH%CELTOT(1 : NB_TOT, K_IDX, SPIN_IDX) * &
                  UK(:, WAN_IDX_2, K_IDX, SPIN_IDX))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDSUBROUTINE CALC_WANNIER_HAM_FROM_TRAFO


   ! _test_mani Helper routine ATM.
   SUBROUTINE ELPH_INTERPOLATE_BANDSTRUCTURE(INTERPOL, IO)
      USE mkpoints_struct_def,   ONLY: kpoints_struct
      USE mkpoints,              ONLY: RD_KPOINTS
      IMPLICIT NONE

      TYPE(ELPHON_INTERPOLATOR), INTENT(IN)  :: INTERPOL
      TYPE(in_struct),           INTENT(IN)  :: IO

      REAL(q), ALLOCATABLE    :: BAND_EV_K(:)
      COMPLEX(q), ALLOCATABLE :: TMP(:, :)
      TYPE(kpoints_struct)    :: Q_PATH
      REAL(q), ALLOCATABLE    :: PATH_QUANTITY(:, :)
      INTEGER                 :: PATH_LEN
      INTEGER                 :: Q_IDX

      ASSOCIATE( &
         PC => INTERPOL%WANNIER_DESCRIPTOR%PRIM_CELL, &
         HAM => INTERPOL%HAM_INIT(1))

         Q_PATH%LTET = .FALSE. ! _test_mani fix
         CALL RD_KPOINTS(Q_PATH, PC%PRIM_LATT, .FALSE., .FALSE., IO%IU6, IO%IU0, FILE_NAME = 'QPOINTS')
         PATH_LEN = Q_PATH%NKPTS

         ALLOCATE(PATH_QUANTITY(HAM%DESCRIPTOR%NUM_PRIM_ORBITALS, PATH_LEN))

         DO Q_IDX = 1, PATH_LEN
            CALL HAM%DIAG_AT_K(Q_PATH%VKPT(:, Q_IDX), BAND_EV_K, TMP)
            PATH_QUANTITY(:, Q_IDX) = BAND_EV_K
         ENDDO

         CALL WRITE_Q_DATA_TO_FILE(PATH_QUANTITY, Q_PATH%VKPT, 'debug_band_struct', PC%PRIM_LATT%B, IO)
      END ASSOCIATE
   END SUBROUTINE ELPH_INTERPOLATE_BANDSTRUCTURE


   ! _test_mani
   ! Calculates output e.g. band structure, el-phon matrix elements, ZPR of band gap, etc.
   SUBROUTINE ELPH_OUTPUT_STUFF(INTERPOL, SETTINGS, T_INFO, E_FERMI, IO)
      USE main_mpi
      USE chi_glb,   ONLY: SHIFT
      IMPLICIT NONE

      TYPE(ELPHON_INTERPOLATOR), INTENT(INOUT)  :: INTERPOL
      TYPE(ELPHON_SETTINGS),     INTENT(IN)     :: SETTINGS
      TYPE(type_info),           INTENT(IN)     :: T_INFO
      REAL(q),                   INTENT(IN)     :: E_FERMI
      TYPE(in_struct),           INTENT(IN)     :: IO

      ! locals
      INTEGER                 :: I, J
      REAL(q), ALLOCATABLE    :: ENERGY_SHIFT(:)
      TYPE(brill_zone_t)      :: BRILL_ZONE
      COMPLEX(q), ALLOCATABLE :: SELF_ENERGY_FM(:, :)
      REAL(q), ALLOCATABLE    :: SELF_ENERGY_DW(:)
      INTEGER                 :: UNIT
      REAL(q), ALLOCATABLE    :: ENERGIES(:)
      REAL(q)                 :: DELTA

      IF (.NOT. SETTINGS%DO_CALC()) RETURN
      
      IF (INTERPOL%IS_COMPLETED .AND. .NOT. INTERPOL%IS_INITIALIZED) CALL vtutor%bug("IS_COMPLETED but not IS_INITIALIZED", __FILE__, __LINE__)

      ASSOCIATE( &
         PC => INTERPOL%WANNIER_DESCRIPTOR%PRIM_CELL)

         IF (IO%IU0 >= 0) WRITE(*, *) "elphon start output"

         CALL PC%SETUP_ADDITIONAL_DATA()

         IF (SETTINGS%CALC_BANDSTRUCTURE) THEN
            IF (IO%IU0 >= 0) WRITE(*, *) "elphon bandstructure"

            IF (.NOT. INTERPOL%IS_INITIALIZED) CALL vtutor%bug( &
               "Bandstructure requested but module not initialized.", __FILE__, __LINE__)
            CALL ELPH_INTERPOLATE_BANDSTRUCTURE(INTERPOL, IO)
         ENDIF

         IF (IO%IU0 >= 0) WRITE(*, *) "elphon polar test"
         CALL ELPH_TEST_COUPLING(INTERPOL, T_INFO, IO)

         IF (SETTINGS%CALC_ZPR) THEN
            IF (IO%IU0 >= 0) WRITE(*, *) "elphon ZPR"

            IF (.NOT. INTERPOL%IS_COMPLETED) CALL vtutor%bug("ZPR requested but module not completed.", __FILE__, __LINE__)

            CALL BRILL_ZONE%INIT(PC%PRIM_LATT)
            CALL BRILL_ZONE%SAMPLE(PC%SUPER_LATT%B, SETTINGS%BZ_SUBDIVISION)
            CALLMPI(M_bcast_d(COMM, BRILL_ZONE%QMESH(1, 1), size(BRILL_ZONE%QMESH))) ! _test_mani Improve MPI handling

            IF (IO%IU0 >= 0) OPEN(NEWUNIT = UNIT, FILE = 'debug_renorm', STATUS = 'replace')
            DO J = 1, SIZE(SETTINGS%BZ_POINTS, 2)
               CALL ELPH_CALC_SELF_ENERGY_FAN_MIGDAL(SELF_ENERGY_FM, INTERPOL, T_INFO, SETTINGS%BZ_POINTS(:, J), BRILL_ZONE%QMESH, 1, 0._q, E_FERMI, SETTINGS%DO_ZPR_DYNAMIC, SHIFT, COMM=COMM, MAX_BAND=SETTINGS%NUM_BANDS_MAX)
               CALL ELPH_CALC_SELF_ENERGY_DEBYE_WALLER(SELF_ENERGY_DW, INTERPOL, T_INFO, SETTINGS%BZ_POINTS(:, J), BRILL_ZONE%QMESH, 1, 0._q, SHIFT, COMM, SETTINGS%NUM_BANDS_MAX)
               IF (IO%IU0 >= 0) THEN
                  IF (J == 1) WRITE(UNIT, '(*(A12))') 'k_x', 'k_y', 'k_z', ('band ' // str(I), I = 1, SIZE(SELF_ENERGY_FM, 2))
                  WRITE(UNIT, '(*(F12.6))') SETTINGS%BZ_POINTS(:, J), (REAL(SELF_ENERGY_FM(1, I)) + SELF_ENERGY_DW(I), I = 1, SIZE(SELF_ENERGY_FM, 2))
               ENDIF
            ENDDO
            IF (IO%IU0 >= 0) CLOSE(UNIT)

            ENERGIES = [REAL(q):: (-50._q + (32._q + 50._q) / 1000 * I, I = 0, 1000)]
            CALL ELPH_CALC_SELF_ENERGY_FAN_MIGDAL(SELF_ENERGY_FM, INTERPOL, T_INFO, SETTINGS%BZ_POINTS(:, 1), BRILL_ZONE%QMESH, 1, 0._q, E_FERMI, SETTINGS%DO_ZPR_DYNAMIC, SHIFT, ENERGIES, COMM)

            IF (IO%IU0 >= 0) THEN
               OPEN(NEWUNIT = UNIT, FILE = 'debug_renorm_real', STATUS = 'replace')
               WRITE(UNIT, '(*(F12.6))') ENERGIES
               DO J = 1, SIZE(SELF_ENERGY_FM, 2)
                  WRITE(UNIT, '(*(F12.6))') REAL(SELF_ENERGY_FM(:, J))
               ENDDO
               CLOSE(UNIT)
               OPEN(NEWUNIT = UNIT, FILE = 'debug_renorm_imag', STATUS = 'replace')
               WRITE(UNIT, '(*(F12.6))') ENERGIES
               DO J = 1, SIZE(SELF_ENERGY_FM, 2)
                  WRITE(UNIT, '(*(F12.6))') AIMAG(SELF_ENERGY_FM(:, J))
               ENDDO
               CLOSE(UNIT)
            ENDIF
         ENDIF

         IF (SETTINGS%CALC_TEMPERATURE_GAP) THEN
            IF (IO%IU0 >= 0) WRITE(*, *) "elphon temperature"

            IF (.NOT. INTERPOL%IS_COMPLETED) CALL vtutor%bug("Temperature dependence requested but module not completed.", __FILE__, __LINE__)

            IF (IO%IU0 >= 0) OPEN(NEWUNIT = UNIT, FILE = 'debug_temperature', STATUS = 'replace')
            DO J = 0, 7
               ! _test_mani implement adequately
               ! CALL ELPH_CALC_ENERGY_SHIFT(SETTINGS%BZ_POINTS(:, 1), 1, J * 100._q, T_INFO, PC, ENERGY_SHIFT)
               ! IF (IO%IU0 >= 0) WRITE(UNIT, '(F16.6, *(F16.6))') J * 100._q, ENERGY_SHIFT
            ENDDO

            IF (IO%IU0 >= 0) CLOSE(UNIT)
         ENDIF
      END ASSOCIATE
   END SUBROUTINE ELPH_OUTPUT_STUFF


   ! _test_mani helper routine ATM.
   SUBROUTINE ELPH_TEST_COUPLING(INTERPOL, T_INFO, IO)
      IMPLICIT NONE

      TYPE(ELPHON_INTERPOLATOR), INTENT(IN)  :: INTERPOL
      TYPE(type_info),           INTENT(IN)  :: T_INFO
      TYPE(in_struct),           INTENT(IN)  :: IO

      IF (.NOT. INTERPOL%IS_COMPLETED) RETURN
      CALL PLOT_COUPLING(INTERPOL%COUPLING_WANN(1), INTERPOL%HAM_INIT(1), [REAL(q):: 0, 0, 0], "debug_elph_element", .TRUE., T_INFO, IO)
      CALL WRITE_COUPLING(INTERPOL%COUPLING_WANN(1), INTERPOL%HAM_INIT(1), [REAL(q):: 0, 0, 0], "debug_elph_mat", T_INFO, IO, 4)
   END SUBROUTINE


   SUBROUTINE PLOT_COUPLING(WANNIER_MAT, HAM_INIT, K_VEC, FILE_NAME, AVERAGE, T_INFO, IO)
      USE mkpoints_struct_def,   ONLY: kpoints_struct
      USE mkpoints,              ONLY: RD_KPOINTS
      USE main_mpi,              ONLY: COMM_WORLD
      USE phonon,                ONLY: PHON_CALC_PHONONS
      IMPLICIT NONE

      TYPE(WANNIER_MAT_MODES_T), INTENT(IN)  :: WANNIER_MAT
      TYPE(WANNIER_MAT_SYMM_T),  INTENT(IN)  :: HAM_INIT
      REAL(q),                   INTENT(IN)  :: K_VEC(3)
      CHARACTER(*),              INTENT(IN)  :: FILE_NAME
      LOGICAL,                   INTENT(IN)  :: AVERAGE
      TYPE(type_info),           INTENT(IN)  :: T_INFO
      TYPE(in_struct),           INTENT(IN)  :: IO

      TYPE(kpoints_struct)             :: Q_PATH
      COMPLEX(q), ALLOCATABLE, TARGET  :: COUPLING(:, :, :)
      COMPLEX(q), POINTER              :: COUPLING_PTR(:, :)
      REAL(q), ALLOCATABLE             :: BAND_EV_K(:)
      REAL(q), ALLOCATABLE             :: FREQU(:)
      COMPLEX(q), ALLOCATABLE          :: TMP(:, :)
      REAL(q), ALLOCATABLE             :: PATH_QUANTITY(:, :, :)
      INTEGER                          :: PATH_LEN
      INTEGER                          :: Q_IDX
      INTEGER                          :: NUM_BANDS
      INTEGER                          :: BAND_IDX
      INTEGER                          :: Q_START, Q_END

      ASSOCIATE(PC => WANNIER_MAT%DESCRIPTOR%PRIM_CELL)

         ! Read the q-path from file.
         Q_PATH%LTET = .FALSE. ! _test_mani fix
         CALL RD_KPOINTS(Q_PATH, PC%PRIM_LATT, .FALSE., .FALSE., IO%IU6, IO%IU0, FILE_NAME = 'QPOINTS')
         PATH_LEN = Q_PATH%NKPTS
         NUM_BANDS = WANNIER_MAT%DESCRIPTOR%NUM_PRIM_ORBITALS

         ALLOCATE(PATH_QUANTITY(WANNIER_MAT%NUM_MODES, PATH_LEN, NUM_BANDS))
         PATH_QUANTITY = 0

         CALL GET_MPI_STRIP(COMM_WORLD, PATH_LEN, Q_START, Q_END)

         DO Q_IDX = Q_START, Q_END
            CALL WANNIER_MAT%TO_BLOCH(COUPLING, HAM_INIT, Q_PATH%VKPT(:, Q_IDX), K_VEC)

            COUPLING_PTR(1 : NUM_BANDS * NUM_BANDS, 1 : WANNIER_MAT%NUM_MODES) => COUPLING
            CALL ROTATE_WITH_MODES(COUPLING_PTR, Q_PATH%VKPT(:, Q_IDX), PC, T_INFO)

            ! Averaging over degenerate branches and bands.
            ! This is only used to output a fancy plot.
            IF (AVERAGE) THEN
               COUPLING = CONJG(COUPLING) * COUPLING
               CALL PHON_CALC_PHONONS(Q_PATH%VKPT(:, Q_IDX), PC, T_INFO, FREQU, TMP)
               CALL AVERAGE_DEGEN(FREQU, 3, COUPLING)
               CALL HAM_INIT%DIAG_AT_K([REAL(q):: 0, 0, 0], BAND_EV_K, TMP)
               CALL AVERAGE_DEGEN(BAND_EV_K, 2, COUPLING)
               ! CALL HAM_INIT%DIAG_AT_K([REAL(q):: 0, 0, 0] + Q_PATH%VKPT(:, Q_IDX), BAND_EV_K, TMP)
               ! CALL AVERAGE_DEGEN(BAND_EV_K, 1, COUPLING)
               COUPLING = SQRT(ABS(COUPLING))
            ENDIF

            DO BAND_IDX = 1, NUM_BANDS
               PATH_QUANTITY(:, Q_IDX, BAND_IDX) = COUPLING(BAND_IDX, BAND_IDX, :)
            ENDDO
         ENDDO

         CALLMPI(M_sum_d(COMM_WORLD, PATH_QUANTITY(1, 1, 1), SIZE(PATH_QUANTITY)))

         DO BAND_IDX = 1, NUM_BANDS
            CALL WRITE_Q_DATA_TO_FILE(PATH_QUANTITY(:, :, BAND_IDX), Q_PATH%VKPT, FILE_NAME // "_" // str(BAND_IDX), PC%PRIM_LATT%B, IO)
         ENDDO

      END ASSOCIATE

      ! Deallocate KPOINTS structure before it goes out of scope
      DEALLOCATE(Q_PATH%VKPT, Q_PATH%WTKPT, Q_PATH%IDTET)

   END SUBROUTINE PLOT_COUPLING


   SUBROUTINE WRITE_COUPLING(WANNIER_MAT, HAM_INIT, K_VEC, FILE_NAME, T_INFO, IO, MAX_BAND)
      USE poscar,    ONLY  : type_info
      USE mkpoints,  ONLY  : kpoints_struct, RD_KPOINTS
      USE main_mpi,  ONLY  : COMM_WORLD
      USE string,    ONLY  : str
      IMPLICIT NONE

      TYPE(WANNIER_MAT_MODES_T), INTENT(IN)  :: WANNIER_MAT
      TYPE(WANNIER_MAT_SYMM_T),  INTENT(IN)  :: HAM_INIT
      REAL(q),                   INTENT(IN)  :: K_VEC(3)
      CHARACTER(*),              INTENT(IN)  :: FILE_NAME
      TYPE(type_info),           INTENT(IN)  :: T_INFO
      TYPE(in_struct),           INTENT(IN)  :: IO
      INTEGER, OPTIONAL,         INTENT(IN)  :: MAX_BAND

      TYPE(kpoints_struct)             :: Q_PATH
      COMPLEX(q), ALLOCATABLE, TARGET  :: COUPLING(:, :, :)
      COMPLEX(q), POINTER              :: COUPLING_PTR(:, :)
      REAL(q), ALLOCATABLE             :: PATH_QUANTITY(:, :, :, :)
      INTEGER                          :: PATH_LEN
      INTEGER                          :: Q_IDX
      INTEGER                          :: MY_MAX_BAND
      INTEGER                          :: NUM_BANDS
      INTEGER                          :: BAND_IDX_1, BAND_IDX_2
      INTEGER                          :: Q_START, Q_END

      ASSOCIATE(PC => WANNIER_MAT%DESCRIPTOR%PRIM_CELL)

         ! Read the q-path from file.
         Q_PATH%LTET = .FALSE. ! _test_mani fix
         CALL RD_KPOINTS(Q_PATH, PC%PRIM_LATT, .FALSE., .FALSE., IO%IU6, IO%IU0, FILE_NAME = 'QPOINTS')
         PATH_LEN = Q_PATH%NKPTS
         NUM_BANDS = WANNIER_MAT%DESCRIPTOR%NUM_PRIM_ORBITALS

         MY_MAX_BAND = NUM_BANDS
         IF (PRESENT(MAX_BAND)) MY_MAX_BAND = MAX_BAND

         ALLOCATE(PATH_QUANTITY(WANNIER_MAT%NUM_MODES * 2, PATH_LEN, MY_MAX_BAND, MY_MAX_BAND))
         PATH_QUANTITY = 0

         CALL GET_MPI_STRIP(COMM_WORLD, PATH_LEN, Q_START, Q_END)

         DO Q_IDX = Q_START, Q_END
            CALL WANNIER_MAT%TO_BLOCH(COUPLING, HAM_INIT, Q_PATH%VKPT(:, Q_IDX), K_VEC)

            COUPLING_PTR(1 : NUM_BANDS * NUM_BANDS, 1 : WANNIER_MAT%NUM_MODES) => COUPLING
            CALL ROTATE_WITH_MODES(COUPLING_PTR, Q_PATH%VKPT(:, Q_IDX), PC, T_INFO)

            DO BAND_IDX_2 = 1, MY_MAX_BAND
               DO BAND_IDX_1 = 1, MY_MAX_BAND
                  PATH_QUANTITY(1::2, Q_IDX, BAND_IDX_1, BAND_IDX_2) = REAL(COUPLING(BAND_IDX_1, BAND_IDX_2, :))
                  PATH_QUANTITY(2::2, Q_IDX, BAND_IDX_1, BAND_IDX_2) = AIMAG(COUPLING(BAND_IDX_1, BAND_IDX_2, :))
               ENDDO
            ENDDO
         ENDDO

         CALLMPI(M_sum_d(COMM_WORLD, PATH_QUANTITY(1, 1, 1, 1), SIZE(PATH_QUANTITY)))

         DO BAND_IDX_2 = 1, MY_MAX_BAND
            DO BAND_IDX_1 = 1, MY_MAX_BAND
               CALL WRITE_Q_DATA_TO_FILE( &
                  PATH_QUANTITY(:, :, BAND_IDX_1, BAND_IDX_2), Q_PATH%VKPT, &
                  FILE_NAME // "_" // str(BAND_IDX_1) // "_" // str(BAND_IDX_2), PC%PRIM_LATT%B, IO)
            ENDDO
         ENDDO

      END ASSOCIATE

      ! Deallocate KPOINTS structure before it goes out of scope
      DEALLOCATE(Q_PATH%VKPT, Q_PATH%WTKPT, Q_PATH%IDTET)

   END SUBROUTINE WRITE_COUPLING


   ! _test_mani q must include 2pi
   ! _test_mani helper routine ATM.
   SUBROUTINE WRITE_Q_DATA_TO_FILE(Q_DATA, Q_PATH, FILENAME, LATT_VECTORS, IO)
      IMPLICIT NONE

      REAL(q),             INTENT(IN)  :: Q_DATA(:, :)
      REAL(q),             INTENT(IN)  :: Q_PATH(:, :)
      CHARACTER(*),        INTENT(IN)  :: FILENAME
      REAL(q),             INTENT(IN)  :: LATT_VECTORS(3, 3)
      TYPE(in_struct),     INTENT(IN)  :: IO

      INTEGER  :: Q_IDX
      INTEGER  :: PATH_LEN
      REAL(q)  :: Q_TOT
      REAL(q)  :: Q_CART(3)
      REAL(q)  :: Q_DIFF, Q_DIFF_OLD
      INTEGER  :: I

      CHARACTER(*), PARAMETER :: C_DELIM = '  '

      IF (IO%IU0 < 0) RETURN

      OPEN(UNIT = 1337, FILE = TRIM(FILENAME), STATUS = 'replace')

      ! Write a header line to the file.
      WRITE(1337, '(A)') &
         '#  NUM'       // C_DELIM // &
         '         q_x' // C_DELIM // &
         '         q_y' // C_DELIM // &
         '         q_z' // C_DELIM // &
         '     q_total' // C_DELIM // &
         ' Data ...'

      PATH_LEN = SIZE(Q_PATH, 2)
      IF (PATH_LEN > 1) THEN
         Q_CART = MATMUL(LATT_VECTORS, Q_PATH(:, 2) - Q_PATH(:, 1))
         Q_DIFF_OLD = SQRT(DOT_PRODUCT(Q_CART, Q_CART))
      ENDIF

      Q_TOT = 0

      DO Q_IDX = 1, PATH_LEN
         ! Write q-point to file
         WRITE(1337, '(I6, A, F12.6, A, F12.6, A, F12.6)', ADVANCE = 'no') &
            Q_IDX, C_DELIM, &
            Q_PATH(1, Q_IDX), C_DELIM, &
            Q_PATH(2, Q_IDX), C_DELIM, &
            Q_PATH(3, Q_IDX)

         ! Write total q-distance traveled
         IF (Q_IDX > 1) THEN
            Q_CART = MATMUL(LATT_VECTORS, Q_PATH(:, Q_IDX) - Q_PATH(:, Q_IDX - 1))
            Q_DIFF = SQRT(DOT_PRODUCT(Q_CART, Q_CART))

            IF (ABS(Q_DIFF) < 1E-7_q .OR. ABS(Q_DIFF - Q_DIFF_OLD) > 1E-7_q) THEN
               IF (Q_IDX < PATH_LEN) THEN
                  Q_CART = MATMUL(LATT_VECTORS, Q_PATH(:, Q_IDX + 1) - Q_PATH(:, Q_IDX))
                  Q_DIFF_OLD = SQRT(DOT_PRODUCT(Q_CART, Q_CART))
               ENDIF
            ELSE
               Q_TOT = Q_TOT + Q_DIFF
            ENDIF
         ENDIF
         WRITE(1337, '(A, F12.6)', ADVANCE = 'no') C_DELIM, Q_TOT

         ! Write couplings to file
         DO I = 1, SIZE(Q_DATA, 1)
            WRITE(1337, '(A, F16.6)', ADVANCE = 'no') C_DELIM, Q_DATA(I, Q_IDX)
         ENDDO
         WRITE(1337, *)
      ENDDO

      CLOSE(1337)
   END SUBROUTINE WRITE_Q_DATA_TO_FILE


   SUBROUTINE ELPH_CALC_SELF_ENERGY_FAN_MIGDAL(SELF_ENERGY, INTERPOL, T_INFO, K_VEC, Q_MESH, SPIN, TEMPERATURE, E_FERMI, DO_DYNAMIC, IMAG_SHIFT, ENERGIES, COMM, MAX_BAND)
      USE mpimy,     ONLY: communic
      USE phonon,    ONLY: PHON_CALC_PHONONS, &
                           PHON_FIND_TRANSLATIONAL_MODES
      USE constant
      IMPLICIT NONE

      COMPLEX(q), ALLOCATABLE,   INTENT(OUT)    :: SELF_ENERGY(:, :)
      TYPE(ELPHON_INTERPOLATOR), INTENT(INOUT)  :: INTERPOL
      TYPE(type_info),           INTENT(IN)     :: T_INFO
      REAL(q),                   INTENT(IN)     :: K_VEC(3)
      REAL(q),                   INTENT(IN)     :: Q_MESH(:, :)
      INTEGER,                   INTENT(IN)     :: SPIN
      REAL(q),                   INTENT(IN)     :: TEMPERATURE
      REAL(q),                   INTENT(IN)     :: E_FERMI
      LOGICAL,                   INTENT(IN)     :: DO_DYNAMIC
      REAL(q),                   INTENT(IN)     :: IMAG_SHIFT
      REAL(q), OPTIONAL,         INTENT(IN)     :: ENERGIES(:)
      TYPE(communic), OPTIONAL,  INTENT(IN)     :: COMM
      INTEGER, OPTIONAL,         INTENT(IN)     :: MAX_BAND

      INTEGER                 :: P_IDX
      INTEGER                 :: Q_IDX, Q_START, Q_END
      INTEGER                 :: MODE_IDX
      INTEGER                 :: BAND_IDX, BAND_IDX_2
      INTEGER                 :: MY_MAX_BAND

      REAL(q)                 :: Q_VEC(3)
      REAL(q), ALLOCATABLE    :: E_PHON(:)
      COMPLEX(q), ALLOCATABLE :: MODES(:, :)
      COMPLEX(q), ALLOCATABLE :: UK_D(:, :)
      COMPLEX(q), ALLOCATABLE :: PRIM_MAT(:, :, :)
      COMPLEX(q), ALLOCATABLE :: COUPLING_0(:, :, :), COUPLING_0S(:, :)
      REAL(q), ALLOCATABLE    :: E_ELEC_K(:), E_ELEC_KQ(:)
      COMPLEX(q), ALLOCATABLE :: COUPLING(:, :, :)
      REAL(q), ALLOCATABLE    :: FERMI_OCC(:), BOSE_OCC(:)
      REAL(q), ALLOCATABLE    :: MASSES(:)
      INTEGER                 :: TRANSLATIONAL_MODES(3)
      COMPLEX(q)              :: C_SHIFT

      COMPLEX(q), ALLOCATABLE, TARGET  :: S_MAT(:, :, :)
      COMPLEX(q), POINTER              :: S_PTR(:, :)

      REAL(q), PARAMETER      :: HBAR = 1.0545718E-24_q / SQRT(EVTOJ) / SQRT(AMTOKG)

      ! This should never happen
      IF (.NOT. INTERPOL%IS_COMPLETED) &
         CALL vtutor%bug("Internal error in ELPH_CALC_COUPLING: Missing Wannier-space data!", __FILE__, __LINE__)

      ASSOCIATE( &
         PC => INTERPOL%WANNIER_DESCRIPTOR%PRIM_CELL, &
         N_ORBS => INTERPOL%WANNIER_DESCRIPTOR%NUM_PRIM_ORBITALS, &
         N_MODES => INTERPOL%COUPLING_WANN(SPIN)%NUM_MODES)

         MY_MAX_BAND = N_ORBS
         IF (PRESENT(MAX_BAND)) THEN
            IF (MAX_BAND > 0) MY_MAX_BAND = MIN(N_ORBS, MAX_BAND)
         ENDIF

         C_SHIFT = CMPLX(0, IMAG_SHIFT, q)

         ! Setup the output matrix
         IF (PRESENT(ENERGIES)) THEN
            ALLOCATE(SELF_ENERGY(SIZE(ENERGIES), N_ORBS))
         ELSE
            ALLOCATE(SELF_ENERGY(1, N_ORBS))
         ENDIF
         SELF_ENERGY = 0

         ALLOCATE(FERMI_OCC(N_ORBS))
         ALLOCATE(BOSE_OCC(N_MODES))

         ! Diagonalize Hamiltonian at k-point
         CALL INTERPOL%HAM_INIT(SPIN)%DIAG_AT_K(K_VEC, E_ELEC_K, UK_D)

         IF (PRESENT(COMM)) THEN
            CALL GET_MPI_STRIP(COMM, SIZE(Q_MESH, 2), Q_START, Q_END)
         ELSE
            Q_START = 1
            Q_END = SIZE(Q_MESH, 2)
         ENDIF

         ! Sum over all q-points
         DO Q_IDX = Q_START, Q_END
            Q_VEC = Q_MESH(:, Q_IDX)

            ! _test_mani Minimize the number of phonon calculations by doing explicitly COUPLING_NO_MODES -> ROTATE_WITH_MODES
            ! Calculate coupling for < k + q | dq | k >
            CALL INTERPOL%CALC_COUPLING(COUPLING, Q_VEC, K_VEC, SPIN, T_INFO)
            ! Calculate eigenvalues at k + q
            CALL INTERPOL%HAM_INIT(SPIN)%DIAG_AT_K(K_VEC + Q_VEC, E_ELEC_KQ, MODES)
            ! Calculate phonons at q
            CALL PHON_CALC_PHONONS(Q_VEC, PC, T_INFO, E_PHON, MODES)
            E_PHON = E_PHON * HBAR

            IF (ALL(ABS(Q_VEC) < 1E-10)) THEN
               TRANSLATIONAL_MODES = PHON_FIND_TRANSLATIONAL_MODES(MODES)
            ELSE
               TRANSLATIONAL_MODES = 0
            ENDIF

            ! _test_mani dynamic test
            CALL INTERPOL%OVERLAP_PERT(SPIN)%TO_BLOCH(S_MAT, INTERPOL%HAM_INIT(SPIN), Q_VEC, K_VEC)
            S_PTR(1 : SIZE(S_MAT, 1) * SIZE(S_MAT, 2), 1 : SIZE(S_MAT, 3)) => S_MAT
            CALL ROTATE_WITH_MODES(S_PTR, Q_VEC, PC, T_INFO)

            FERMI_OCC = CALC_FERMI_DIRAC(TEMPERATURE, E_ELEC_KQ)
            BOSE_OCC = CALC_BOSE_EINSTEIN(TEMPERATURE, E_PHON)

            ! Sum over all phonon modes
            DO MODE_IDX = 1, N_MODES
               ! Skip the purely translational modes
               IF (ANY(MODE_IDX == TRANSLATIONAL_MODES)) CYCLE

               !  _________________________________________
               ! |                                         |
               ! |      Fan Migdal Term (regular part)     |
               ! |_________________________________________|

               DO BAND_IDX = 1, N_ORBS
                  IF (PRESENT(ENERGIES)) THEN
                     DO BAND_IDX_2 = 1, MY_MAX_BAND
                        SELF_ENERGY(:, BAND_IDX) = SELF_ENERGY(:, BAND_IDX) + &
                           ! CONJG(COUPLING(BAND_IDX_2, BAND_IDX, MODE_IDX)) * COUPLING(BAND_IDX_2, BAND_IDX, MODE_IDX) * &
                           CONJG(COUPLING(BAND_IDX_2, BAND_IDX, MODE_IDX) + (E_ELEC_K(BAND_IDX) - ENERGIES(:)) * S_MAT(BAND_IDX_2, BAND_IDX, MODE_IDX)) * &
                           (COUPLING(BAND_IDX_2, BAND_IDX, MODE_IDX) + (E_ELEC_K(BAND_IDX) - ENERGIES(:)) * S_MAT(BAND_IDX_2, BAND_IDX, MODE_IDX)) * &
                           (2._q * BOSE_OCC(MODE_IDX) + 1._q) / (ENERGIES(:) - E_ELEC_KQ(BAND_IDX_2) + C_SHIFT)
                     ENDDO
                  ELSE
                     IF (DO_DYNAMIC) THEN
                        SELF_ENERGY(1, BAND_IDX) = SELF_ENERGY(1, BAND_IDX) + SUM( &
                           ! CONJG(COUPLING(:MY_MAX_BAND, BAND_IDX, MODE_IDX)) * COUPLING(:MY_MAX_BAND, BAND_IDX, MODE_IDX) * ( &
                           CONJG(COUPLING(:, BAND_IDX, MODE_IDX)) * (COUPLING(:, BAND_IDX, MODE_IDX) + (SIGN(E_PHON(MODE_IDX), FERMI_OCC(:) - .5_q) - E_ELEC_K(BAND_IDX)) * S_MAT(:, BAND_IDX, MODE_IDX)) * ( &
                              (1._q - FERMI_OCC(:MY_MAX_BAND) + BOSE_OCC(MODE_IDX)) / (E_ELEC_K(BAND_IDX) - E_ELEC_KQ(:MY_MAX_BAND) - E_PHON(MODE_IDX) + C_SHIFT) + &
                              (FERMI_OCC(:MY_MAX_BAND) + BOSE_OCC(MODE_IDX)) / (E_ELEC_K(BAND_IDX) - E_ELEC_KQ(:MY_MAX_BAND) + E_PHON(MODE_IDX) + C_SHIFT) &
                           ) &
                        )
                     ELSE
                        SELF_ENERGY(1, BAND_IDX) = SELF_ENERGY(1, BAND_IDX) + SUM( &
                           CONJG(COUPLING(:MY_MAX_BAND, BAND_IDX, MODE_IDX)) * COUPLING(:MY_MAX_BAND, BAND_IDX, MODE_IDX) * &
                           (2._q * BOSE_OCC(MODE_IDX) + 1._q) / (E_ELEC_K(BAND_IDX) - E_ELEC_KQ(:MY_MAX_BAND) + C_SHIFT) &
                        )
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         ! Sum up the individual contributions to the energy shift from each node
         IF (PRESENT(COMM)) THEN
            CALLMPI(M_sum_z(COMM, SELF_ENERGY(1, 1), SIZE(SELF_ENERGY)))
         ENDIF

         ! Part of integration
         SELF_ENERGY = SELF_ENERGY / SIZE(Q_MESH, 2)

         !  _________________________________________
         ! |                                         |
         ! |      Fan Migdal Term (special part)     |
         ! |_________________________________________|

         ! Calculate coupling at q=0 without rotating phononic DOFs
         CALL INTERPOL%CALC_COUPLING_ATOMS(COUPLING_0, [0._q, 0._q, 0._q], K_VEC, SPIN)
         ALLOCATE(COUPLING_0S(N_ORBS, N_MODES))

         ! Calculate dS at q=0 without rotating phononic DOFs
         ! Fourier sum of dS
         CALL INTERPOL%OVERLAP_PERT(SPIN)%FOURIER_SUM(PRIM_MAT, K_VEC)
         ! Get < n,k | dS | n,k > by calculating the diagonal entries of (U * dS * U^Dagger)
         ! DOT_PRODUCT(a, b) takes SUM(CONJG(a) * b), hence UK_D needs transposing only
         DO MODE_IDX = 1, N_MODES
            DO BAND_IDX = 1, N_ORBS
               COUPLING_0S(BAND_IDX, MODE_IDX) = DOT_PRODUCT(UK_D(:, BAND_IDX), MATMUL(PRIM_MAT(:, :, MODE_IDX), UK_D(:, BAND_IDX)))
            ENDDO
         ENDDO

         ! The matrix elements can be obtained from COUPLING_0 and COUPLING_0S by rotating with phonon modes
         CALL PHON_CALC_PHONONS([0._q, 0._q, 0._q], PC, T_INFO, E_PHON, MODES)
         E_PHON = E_PHON * HBAR
         TRANSLATIONAL_MODES = PHON_FIND_TRANSLATIONAL_MODES(MODES)

         ! _test_mani array indexing
         ALLOCATE(MASSES(PC%NUM_PRIM_ATOMS))
         DO P_IDX = 1, PC%NUM_PRIM_ATOMS
            MASSES(P_IDX) = T_INFO%POMASS(T_INFO%ITYP(PC%SUPER_IDX(1, P_IDX)))
         ENDDO

         ! Rescale the phonon rotation matrix with the masses
         DO P_IDX = 1, PC%NUM_PRIM_ATOMS
            MODES(P_IDX * 3 - 2 : P_IDX * 3, :) = MODES(P_IDX * 3 - 2 : P_IDX * 3, :) / SQRT(MASSES(P_IDX))
         ENDDO

         BOSE_OCC = (2 * CALC_BOSE_EINSTEIN(TEMPERATURE, E_PHON) + 1) * HBAR * HBAR / (2 * E_PHON)

         DO MODE_IDX = 1, N_MODES
            ! Skip the purely translational modes
            ! _test_mani does this always work?! also, check the other instances of that condition
            IF (ANY(MODE_IDX == TRANSLATIONAL_MODES)) CYCLE

            ! _test_mani Performance
            DO BAND_IDX = 1, N_ORBS
               SELF_ENERGY(:, BAND_IDX) = SELF_ENERGY(:, BAND_IDX) - &
               SUM(COUPLING_0(BAND_IDX, BAND_IDX, :) * MODES(:, MODE_IDX)) * SUM(COUPLING_0S(BAND_IDX, :) * MODES(:, MODE_IDX)) * BOSE_OCC(MODE_IDX)
            ENDDO
         ENDDO
      END ASSOCIATE
   END SUBROUTINE ELPH_CALC_SELF_ENERGY_FAN_MIGDAL


   SUBROUTINE ELPH_CALC_SELF_ENERGY_DEBYE_WALLER(SELF_ENERGY, INTERPOL, T_INFO, K_VEC, Q_MESH, SPIN, TEMPERATURE, IMAG_SHIFT, COMM, MAX_BAND)
      USE mpimy,     ONLY: communic
      USE phonon,    ONLY: PHON_CALC_PHONONS, &
                           PHON_FIND_TRANSLATIONAL_MODES
      USE constant
      IMPLICIT NONE

      REAL(q), ALLOCATABLE,      INTENT(OUT)    :: SELF_ENERGY(:)
      TYPE(ELPHON_INTERPOLATOR), INTENT(INOUT)  :: INTERPOL
      TYPE(type_info),           INTENT(IN)     :: T_INFO
      REAL(q),                   INTENT(IN)     :: K_VEC(3)
      REAL(q),                   INTENT(IN)     :: Q_MESH(:, :)
      INTEGER,                   INTENT(IN)     :: SPIN
      REAL(q),                   INTENT(IN)     :: TEMPERATURE
      REAL(q),                   INTENT(IN)     :: IMAG_SHIFT
      TYPE(communic), OPTIONAL,  INTENT(IN)     :: COMM
      INTEGER, OPTIONAL,         INTENT(IN)     :: MAX_BAND

      INTEGER                 :: P_IDX_1, P_IDX_2
      INTEGER                 :: Q_IDX, Q_START, Q_END
      INTEGER                 :: MODE_IDX
      INTEGER                 :: MODE_IDX_1, MODE_IDX_2
      INTEGER                 :: BAND_IDX_1, BAND_IDX_2
      INTEGER                 :: CART_IDX_1, CART_IDX_2
      INTEGER                 :: MY_MAX_BAND

      REAL(q)                 :: Q_VEC(3)
      REAL(q), ALLOCATABLE    :: FREQU(:)
      COMPLEX(q), ALLOCATABLE :: MODES(:, :)
      COMPLEX(q), ALLOCATABLE :: UK_D(:, :)
      COMPLEX(q), ALLOCATABLE :: PRIM_MAT(:, :, :)
      COMPLEX(q), ALLOCATABLE :: COUPLING_0(:, :, :), COUPLING_0S(:, :)
      REAL(q), ALLOCATABLE    :: E_ELEC_K(:)
      REAL(q)                 :: OCC_FACTOR
      COMPLEX(q), ALLOCATABLE :: DW_FACTORS(:, :)
      COMPLEX(q)              :: DW_COUPLING
      REAL(q), ALLOCATABLE    :: MASSES(:)
      INTEGER                 :: TRANSLATIONAL_MODES(3)
      COMPLEX(q)              :: C_SHIFT

      REAL(q), PARAMETER      :: HBAR = 1.0545718E-24_q / SQRT(EVTOJ) / SQRT(AMTOKG)

      ! This should never happen
      IF (.NOT. INTERPOL%IS_COMPLETED) &
         CALL vtutor%bug("Internal error in ELPH_CALC_COUPLING: Missing Wannier-space data!", __FILE__, __LINE__)

      ASSOCIATE( &
         PC => INTERPOL%WANNIER_DESCRIPTOR%PRIM_CELL, &
         N_ORBS => INTERPOL%WANNIER_DESCRIPTOR%NUM_PRIM_ORBITALS, &
         N_MODES => INTERPOL%COUPLING_WANN(SPIN)%NUM_MODES)

         MY_MAX_BAND = N_ORBS
         IF (PRESENT(MAX_BAND)) THEN
            IF (MAX_BAND > 0) MY_MAX_BAND = MIN(N_ORBS, MAX_BAND)
         ENDIF

         C_SHIFT = CMPLX(0, IMAG_SHIFT, q)

         ! _test_mani array indexing
         ALLOCATE(MASSES(PC%NUM_PRIM_ATOMS))
         DO P_IDX_1 = 1, PC%NUM_PRIM_ATOMS
            MASSES(P_IDX_1) = T_INFO%POMASS(T_INFO%ITYP(PC%SUPER_IDX(1, P_IDX_1)))
         ENDDO

         ALLOCATE(DW_FACTORS(N_MODES, N_MODES))
         DW_FACTORS = 0

         ! Setup the output matrix
         ALLOCATE(SELF_ENERGY(N_ORBS))
         SELF_ENERGY = 0

         ! Diagonalize Hamiltonian at k-point
         CALL INTERPOL%HAM_INIT(SPIN)%DIAG_AT_K(K_VEC, E_ELEC_K, UK_D)

         IF (PRESENT(COMM)) THEN
            CALL GET_MPI_STRIP(COMM, SIZE(Q_MESH, 2), Q_START, Q_END)
         ELSE
            Q_START = 1
            Q_END = SIZE(Q_MESH, 2)
         ENDIF

         ! Sum over all q-points
         DO Q_IDX = Q_START, Q_END
            Q_VEC = Q_MESH(:, Q_IDX)

            ! Calculate phonons at q
            CALL PHON_CALC_PHONONS(Q_VEC, PC, T_INFO, FREQU, MODES)

            IF (ALL(ABS(Q_VEC) < 1E-10)) THEN
               TRANSLATIONAL_MODES = PHON_FIND_TRANSLATIONAL_MODES(MODES)
            ELSE
               TRANSLATIONAL_MODES = 0
            ENDIF

            ! Sum over all phonon modes
            DO MODE_IDX = 1, N_MODES
               ! Skip the purely translational modes
               IF (ANY(MODE_IDX == TRANSLATIONAL_MODES)) CYCLE

               ! Calculate the factor involving the Bose-Einstein distribution.
               OCC_FACTOR = 2 * CALC_BOSE_EINSTEIN(TEMPERATURE, FREQU(MODE_IDX) * HBAR) + 1

               !  ________________________________________________
               ! |                                                |
               ! |      Debye Waller Term (q-dependent part)      |
               ! |________________________________________________|

               ! Calculate the matrix involving the phonon modes (Theta)
               ! Include occupation factor and frequencies and sum over all modes
               DO MODE_IDX_2 = 1, N_MODES
                  P_IDX_2 = (MODE_IDX_2 + 2) / 3
                  CART_IDX_2 = MODULO(MODE_IDX_2 - 1, 3) + 1
                  DO MODE_IDX_1 = 1, N_MODES
                     P_IDX_1 = (MODE_IDX_1 + 2) / 3
                     CART_IDX_1 = MODULO(MODE_IDX_1 - 1, 3) + 1

                     DW_FACTORS(MODE_IDX_1, MODE_IDX_2) = DW_FACTORS(MODE_IDX_1, MODE_IDX_2) + &
                     (MODES(MODE_IDX_1, MODE_IDX) * CONJG(MODES((P_IDX_1 - 1) * 3 + CART_IDX_2, MODE_IDX)) / MASSES(P_IDX_1) + &
                     CONJG(MODES((P_IDX_2 - 1 ) * 3 + CART_IDX_1, MODE_IDX)) * MODES(MODE_IDX_2, MODE_IDX) / MASSES(P_IDX_2)) * &
                     OCC_FACTOR / FREQU(MODE_IDX)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         !  ______________________________________________
         ! |                                              |
         ! |      Debye Waller Term (remaining part)      |
         ! |______________________________________________|

         ! Calculate coupling at q=0 without rotating phononic DOFs
         CALL INTERPOL%CALC_COUPLING_ATOMS(COUPLING_0, [0._q, 0._q, 0._q], K_VEC, SPIN)
         ALLOCATE(COUPLING_0S(N_ORBS, N_MODES))

         ! Calculate dS at q=0 without rotating phononic DOFs
         ! Fourier sum of dS.
         CALL INTERPOL%OVERLAP_PERT(SPIN)%FOURIER_SUM(PRIM_MAT, K_VEC)
         ! Get < n,k | dS | n,k > by calculating the diagonal entries of (U * dS * U^Dagger)
         ! DOT_PRODUCT(a, b) takes SUM(CONJG(a) * b), hence UK_D needs transposing only
         DO MODE_IDX = 1, N_MODES
            DO BAND_IDX_1 = 1, N_ORBS
               COUPLING_0S(BAND_IDX_1, MODE_IDX) = DOT_PRODUCT(UK_D(:, BAND_IDX_1), MATMUL(PRIM_MAT(:, :, MODE_IDX), UK_D(:, BAND_IDX_1)))
            ENDDO
         ENDDO

         DW_FACTORS = DW_FACTORS * HBAR / 4._q

         ! Contract the DW factor matrix with the DW couplings
         DO MODE_IDX_2 = 1, N_MODES
            DO MODE_IDX_1 = 1, N_MODES
               DO BAND_IDX_2 = 1, N_ORBS
                  ! Special PAW contribution
                  DW_COUPLING = COUPLING_0(BAND_IDX_2, BAND_IDX_2, MODE_IDX_1) * COUPLING_0S(BAND_IDX_2, MODE_IDX_2)
                  ! Regular contribution
                  DW_COUPLING = DW_COUPLING - SUM( &
                     REAL(CONJG(COUPLING_0(:MY_MAX_BAND, BAND_IDX_2, MODE_IDX_1)) * COUPLING_0(:MY_MAX_BAND, BAND_IDX_2, MODE_IDX_2)) / &
                     (E_ELEC_K(BAND_IDX_2) - E_ELEC_K(:MY_MAX_BAND) + C_SHIFT) &
                  )

                  ! Add everything to the self energy
                  SELF_ENERGY(BAND_IDX_2) = SELF_ENERGY(BAND_IDX_2) + DW_COUPLING * DW_FACTORS(MODE_IDX_1, MODE_IDX_2)
               ENDDO
            ENDDO
         ENDDO

         ! Sum up the individual contributions to the energy shift from each node
         IF (PRESENT(COMM)) THEN
            CALLMPI(M_sum_d(COMM, SELF_ENERGY(1), N_ORBS))
         ENDIF

         ! Part of integration
         SELF_ENERGY = SELF_ENERGY / SIZE(Q_MESH, 2)
      END ASSOCIATE
   END SUBROUTINE


!                             ____________________
!____________________________/ CALC_BOSE_EINSTEIN \_____________________________
!
!> @brief Calculates the bosonic occupation number.
!>
!> \f[
!>    n = \frac{1}{\exp\qty(\frac{E}{k_\text{B} T}) - 1}
!> \f]
!>
!> @param[in] TEMPERATURE Temperature in Kelvin
!> @param[in] ENERGY Energy in eV
!> @return Occupation number
!_______________________________________________________________________________
   ELEMENTAL REAL(q) FUNCTION CALC_BOSE_EINSTEIN(TEMPERATURE, ENERGY) RESULT(OCC)
      USE constant
      IMPLICIT NONE

      REAL(q), INTENT(IN)  :: TEMPERATURE
      REAL(q), INTENT(IN)  :: ENERGY

      ! We run into potential problems at hard (negative) zero (IEEE -> -Infinity)
      IF (TEMPERATURE < TINY(1._q)) THEN
         OCC = 0._q
      ELSE
         ! General case
         OCC = 1._q / (EXP(ENERGY / (BOLKEV * TEMPERATURE)) - 1._q)
      ENDIF
   END FUNCTION CALC_BOSE_EINSTEIN


!                              __________________
!_____________________________/ CALC_FERMI_DIRAC \______________________________
!
!> @brief Calculates the fermionic occupation number.
!>
!> \f[
!>    f = \frac{1}{\exp\qty(\frac{E}{k_\text{B} T}) + 1}
!> \f]
!>
!> @param[in] TEMPERATURE Temperature in Kelvin
!> @param[in] ENERGY Energy in eV
!> @return Occupation number
!_______________________________________________________________________________
   ELEMENTAL REAL(q) FUNCTION CALC_FERMI_DIRAC(TEMPERATURE, ENERGY) RESULT(OCC)
      USE constant
      IMPLICIT NONE

      REAL(q), INTENT(IN)  :: TEMPERATURE
      REAL(q), INTENT(IN)  :: ENERGY

      ! We run into potential problems at hard (negative) zero (IEEE -> -Infinity)
      IF (TEMPERATURE < TINY(1._q)) THEN
         OCC = MERGE(0._q, 1._q, ENERGY > 0._q)
      ELSE
         ! General case
         OCC = 1._q / (EXP(ENERGY / (BOLKEV * TEMPERATURE)) + 1._q)
      ENDIF
   END FUNCTION CALC_FERMI_DIRAC


!                          ___________________________
!_________________________/ REORDER_LOCPROJ_FUNCTIONS \_________________________
!
!> @brief Rearranges columns of LOCPROJ projection matrix.
!>
!> The function indices as obtained by LOCPROJ follow the same pattern as
!> specified in the INCAR file. In order to be used by our Wannier
!> interpolation, the indices need to be rearranged such that all functions
!> (columns) are grouped first by cell index and then by primitive-atom index.
!>
!> @param[in,out] U LOCPROJ matrix to be rearranged
!> @param[in] MLWF_T Contains site information for local orbitals
!> @param[in] PRIM_CELL Primitive-cell information
!_______________________________________________________________________________
   SUBROUTINE REORDER_LOCPROJ_FUNCTIONS(U, MLWF_T, PRIM_CELL)
      USE mlwf,   ONLY: MLWF_type
      IMPLICIT NONE

      GDEF,                INTENT(INOUT)  :: U(:, :)
      TYPE(MLWF_type),     INTENT(IN)     :: MLWF_T
      TYPE(PRIM_CELL_T),   INTENT(IN)     :: PRIM_CELL

      INTEGER           :: CUR_IDX, FUN_IDX
      INTEGER           :: CELL_IDX, PRIM_IDX, S_IDX
      INTEGER           :: NUM_FUNS
      GDEF, ALLOCATABLE :: BUFFER(:, :)

      NUM_FUNS = SIZE(U, 2)
      ALLOCATE(BUFFER(SIZE(U, 1), NUM_FUNS))

      CUR_IDX = 0

      DO CELL_IDX = 1, PRIM_CELL%NUM_CELLS
         DO PRIM_IDX = 1, PRIM_CELL%NUM_PRIM_ATOMS
            S_IDX = PRIM_CELL%SUPER_IDX(CELL_IDX, PRIM_IDX)
            DO FUN_IDX = 1, NUM_FUNS
               IF (MLWF_T%LPRJ_functions(FUN_IDX)%poscar_site == S_IDX) THEN
                  CUR_IDX = CUR_IDX + 1
                  BUFFER(:, CUR_IDX) = U(:, FUN_IDX)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      
      U = BUFFER
   END SUBROUTINE


   SUBROUTINE TEST_HAM_MAT(HAM_MAT, ENERGIES, IO)
      IMPLICIT NONE

      GDEF,             INTENT(IN)  :: HAM_MAT(:, :)
      COMPLEX(q),       INTENT(IN)  :: ENERGIES(:)
      TYPE(in_struct),  INTENT(IN)  :: IO

      GDEF, ALLOCATABLE       :: EIGENVECTORS(:, :)
      INTEGER                 :: MAT_SIZE
      REAL(q), ALLOCATABLE    :: EIGENVALUES(:)
      REAL(q), ALLOCATABLE    :: WORK(:)
      INTEGER                 :: LWORK
      INTEGER                 :: IFAIL
      INTEGER                 :: UNIT
      INTEGER                 :: I

      IF (IO%IU0 < 0) RETURN

      MAT_SIZE = SIZE(HAM_MAT, 1)

      ALLOCATE(EIGENVECTORS(MAT_SIZE, MAT_SIZE), EIGENVALUES(MAT_SIZE))
      EIGENVECTORS = HAM_MAT

#ifdef gammareal
      ! Perform the diagonalization.
      LWORK = 3 * MAT_SIZE
      ALLOCATE(WORK(LWORK))
      CALL DSYEV('N', 'U', MAT_SIZE, EIGENVECTORS, MAT_SIZE, EIGENVALUES, WORK, LWORK, IFAIL)
      
      IF (IFAIL /= 0) CALL vtutor%error("TEST_HAM_MAT: DSYEV returned " // str(IFAIL) // ".")
#else
      CALL vtutor%bug("TEST_HAM_MAT: Complex not implemented!", __FILE__, __LINE__)
#endif
      OPEN(NEWUNIT = UNIT, FILE = 'debug_eigenvalues', STATUS = 'replace')
         WRITE(UNIT, '(2F12.6)') (REAL(ENERGIES(I)), EIGENVALUES(I), I = 1, MAT_SIZE)
      CLOSE(UNIT)
   END SUBROUTINE TEST_HAM_MAT


!                              ___________________
!_____________________________/ ELPH_INIT_WANNIER \_____________________________
!
!> @brief Begins the initialization of the electron-phonon module.
!>
!> A set of unperturbed Wannier functions is created from the ground-state data
!> and the unperturbed Hamiltonian matrix is calculated in Wannier space. The
!> routine needs to be called after electronic minimization has taken place but
!> before any ions are displaced. If initialization was successful,
!> ELPH_IS_INITIALIZED will return true. If, for any reason, a reset is
!> required, this routine can be called again to reinitialize the module.
!_______________________________________________________________________________
   SUBROUTINE EPI_INIT(THIS, SETTINGS, MLWF_T, W, CQIJ, CDIJ, P, PRIM_CELL, T_INFO, INFO, &
      GRID, GRIDC, GRIDUS, C_TO_US, NONL_S, NONLR_S, IRDMAX, IO)
      USE pseudo_struct_def,  ONLY: potcar
      USE base,               ONLY: info_struct
      USE mgrid_struct_def,   ONLY: grid_3d, &
                                    transit
      USE nonl_struct_def,    ONLY: nonl_struct
      USE nonlr_struct_def,   ONLY: nonlr_struct
      USE mlwf,               ONLY: MLWF_type, &
                                    CALC_WANNIER_TRAFO_PROJ
      USE locproj,            ONLY: LPRJ_LINCOM, LPRJ_IS_INITIALIZED, LPRJ_FUNCTIONS, &
                                    LPRJ_COVL, LPRJ_PROALL
      USE wave,               ONLY: ALLOCW, &
                                    DEALLOCW
      IMPLICIT NONE

      CLASS(ELPHON_INTERPOLATOR),   INTENT(INOUT)  :: THIS
      TYPE(ELPHON_SETTINGS),        INTENT(IN)     :: SETTINGS
      TYPE(MLWF_type),              INTENT(IN)     :: MLWF_T
      TYPE(wavespin),               INTENT(IN)     :: W
      OVERLAP,                      INTENT(IN)     :: CQIJ(:, :, :, :)
      OVERLAP,                      INTENT(IN)     :: CDIJ(:, :, :, :)
      TYPE(potcar),                 INTENT(INOUT)  :: P(:) ! NABLA array is set
      TYPE(PRIM_CELL_T), TARGET,    INTENT(INOUT)  :: PRIM_CELL
      TYPE(type_info),              INTENT(IN)     :: T_INFO
      TYPE(info_struct),            INTENT(IN)     :: INFO
      TYPE(grid_3d),                INTENT(IN)     :: GRID
      TYPE(grid_3d),                INTENT(IN)     :: GRIDC
      TYPE(grid_3d),                INTENT(IN)     :: GRIDUS
      TYPE(transit),                INTENT(IN)     :: C_TO_US
      TYPE(nonl_struct),            INTENT(IN)     :: NONL_S
      TYPE(nonlr_struct),           INTENT(IN)     :: NONLR_S
      INTEGER,                      INTENT(IN)     :: IRDMAX
      TYPE(in_struct),              INTENT(IN)     :: IO

      GDEF, ALLOCATABLE :: UK(:, :, :, :)
      GDEF, ALLOCATABLE :: HAM_INIT_MAT(:, :, :, :)
      INTEGER           :: SPIN_IDX
      INTEGER           :: NUM_ORBS_IN_CELL
      LOGICAL           :: FORCE_ATOM_CENTERED

      CALL THIS%DEALLOC()

      ! Allocate storage for Wannier function
      CALL ALLOCW(W%WDES, THIS%W_WAN_INIT)

      FORCE_ATOM_CENTERED = .FALSE.

      SELECT CASE (SETTINGS%METHOD)
         CASE (METHOD_SCDM)
            NUM_ORBS_IN_CELL = SETTINGS%NUM_WANN_READ
            CALL CALC_WANNIER_ORBITALS_SCDM(W, NUM_ORBS_IN_CELL, THIS%W_WAN_INIT, UK, SETTINGS%SCDM_MU, SETTINGS%SCDM_SIGMA, PRIM_CELL, NONL_S, NONLR_S, CQIJ, IO)
         CASE (METHOD_LOAD_PRIM)
            NUM_ORBS_IN_CELL = SETTINGS%NUM_WANN_READ
            CALL CALC_WANNIER_ORBITALS_FROM_PRIM(W, NUM_ORBS_IN_CELL, THIS%W_WAN_INIT, UK, PRIM_CELL, NONL_S, NONLR_S, CQIJ, IO)
         CASE (METHOD_LOCPROJ)
            IF (.NOT. ALLOCATED(MLWF_T%A_matrix) .OR. .NOT. ALLOCATED(MLWF_T%LPRJ_functions)) &
               CALL vtutor%error("ELPH_INIT_WANNIER: Projections not set up properly!")

            NUM_ORBS_IN_CELL = MLWF_T%NUM_WANN / PRIM_CELL%NUM_CELLS ! _test_mani LOCPROJ needs to be checked for translational invariance!

            CALL CALC_WANNIER_TRAFO_PROJ(UK, W, MLWF_T)

            ! Reorder the columns so they are grouped by site index
            DO SPIN_IDX = 1, W%WDES%ISPIN
               CALL REORDER_LOCPROJ_FUNCTIONS(UK(:, :, 1, SPIN_IDX), MLWF_T, PRIM_CELL)
            ENDDO

            ! Rotate to Wannier functions and store them in W_WAN_INIT.
            CALL LPRJ_LINCOM(W, UK, THIS%W_WAN_INIT)

            FORCE_ATOM_CENTERED = .TRUE.
         CASE DEFAULT
            CALL vtutor%bug("No elphon method selected!", __FILE__, __LINE__)
      END SELECT

      ! Setup Wannier-orbital descriptor.
      THIS%WANNIER_DESCRIPTOR = WANNIER_DESCRIPTOR_T( &
         THIS%W_WAN_INIT, NUM_ORBS_IN_CELL, PRIM_CELL, GRID, GRIDC, GRIDUS, &
         C_TO_US, T_INFO, P, SIZE(CQIJ, 1), INFO%LOVERL, IRDMAX, FORCE_ATOM_CENTERED)

      ! _test_mani export stuff
      ! IF (ELPH_METHOD == METHOD_LOCPROJ) THEN
      !    ALLOCATE(TEST_ARR(W_WAN_INIT%WDES%NRPLWV, 2))
      !    TEST_ARR(:, 1) = W_WAN_INIT%CPTWFP(:, 1, 1, 1)

      !    CALL WRITE_PRIM_WANN_FUNCS(THIS%W_WAN_INIT, THIS%WANNIER_DESCRIPTOR, T_INFO, IO)
      !    CALL CALC_WANNIER_ORBITALS_FROM_PRIM(W, NUM_ORBS_IN_CELL, W_WAN_INIT, UK, PRIM_CELL, NONL_S, NONLR_S, CQIJ, IO)

      !    TEST_ARR(:, 2) = W_WAN_INIT%CPTWFP(:, 1, 1, 1)

      !    OPEN(NEWUNIT = UNIT, FILE = 'debug_w_diff', STATUS = 'replace')
      !    WRITE(UNIT, '(3("(", F12.6, " ", F12.6, ")   "))') (TEST_ARR(I, 1), TEST_ARR(I, 2), TEST_ARR(I, 2) - TEST_ARR(I, 1), I = 1, W_WAN_INIT%WDES%NRPLWV)
      !    CLOSE(UNIT)
      !    WRITE(*, *) MAXVAL(ABS(TEST_ARR(:, 1) - TEST_ARR(:, 2)))
         
      !    WANNIER_DESCRIPTOR = WANNIER_DESCRIPTOR_T( &
      !       W_WAN_INIT, NUM_ORBS_IN_CELL, PRIM_CELL, GRID, GRIDC, GRIDUS, &
      !       C_TO_US, T_INFO, P, SIZE(CQIJ, 1), INFO%LOVERL, IRDMAX)
      ! ENDIF

      ! Calculate the unperturbed Wannier-space Hamiltonian matrix.
      CALL CALC_WANNIER_HAM_FROM_TRAFO(W, UK, HAM_INIT_MAT) ! _test_mani maybe replace with MLWF version

      ! Setup the initial Hamiltonian data structure.
      ALLOCATE(THIS%HAM_INIT(W%WDES%ISPIN))
      DO SPIN_IDX = 1, W%WDES%ISPIN
         THIS%HAM_INIT(SPIN_IDX) = WANNIER_MAT_SYMM_T(THIS%WANNIER_DESCRIPTOR, &
            HAM_INIT_MAT(1 : NUM_ORBS_IN_CELL * PRIM_CELL%NUM_CELLS, 1 : NUM_ORBS_IN_CELL, 1, SPIN_IDX))
      ENDDO

      ALLOCATE(THIS%CDIJ_INIT(SIZE(CDIJ, 1), SIZE(CDIJ, 2), SIZE(CDIJ, 3), SIZE(CDIJ, 4)))
      THIS%CDIJ_INIT = CDIJ

      IF (W%WDES%LOVERL) THEN
         ! Calculate derivative of partial waves
         ! CALL SET_NABIJ_AUG(P, SIZE(P))
         ! Alternative
         CALL CALC_NABIJ_2(P)
      ENDIF

      ! Module is initiliazed.
      THIS%IS_INITIALIZED = .TRUE.
      THIS%IS_COMPLETED = .FALSE.
   END SUBROUTINE


!                           ________________________
!__________________________/ EPI_PERTURB_WANNIER_PS \___________________________
!
!> @brief Calculates Wannier-space matrices after each displacement step.
!>
!> Calculates the electron-phonon and S-overlap matrix elements in the Wannier
!> representation using finite differences. It needs to be called after each
!> ionic displacmenet and the subsequent electronic minimization. Only after
!> the last call, i.e. after all displacements have been performed, is the
!> module properly setup and ready to use. This can be queried by a call to
!> ELPH_IS_COMPLETED.
!_______________________________________________________________________________
   SUBROUTINE EPI_PERTURB_WANNIER_PS(THIS, GRID, P, LATT_CUR, T_INFO, &
      HAMILTONIAN, CQIJ, CDIJ, SV, NONL_S, NONLR_S, SYMM, EXHF, IO)
      USE mgrid_struct_def,   ONLY: grid_3d
      USE pseudo_struct_def,  ONLY: potcar
      USE hamil_struct_def,   ONLY: ham_handle
      USE nonl_struct_def,    ONLY: nonl_struct
      USE nonlr_struct_def,   ONLY: nonlr_struct
      USE base,               ONLY: symmetry
      USE finite_differences, ONLY: GLOB_FINDIF_INFO
      USE wave,               ONLY: SETWDES
      USE nonl_high,          ONLY: PROALL
      USE wave_high,          ONLY: ELEMENTS
      USE subrot,             ONLY: EDDIAG
      USE nonl,               ONLY: PHASE
      USE choleski,           ONLY: CALC_PAW_OVERLAP
      IMPLICIT NONE

      CLASS(ELPHON_INTERPOLATOR),   INTENT(INOUT)  :: THIS
      TYPE(grid_3d),                INTENT(IN)     :: GRID
      TYPE(potcar),                 INTENT(IN)     :: P(:)
      TYPE(latt),                   INTENT(IN)     :: LATT_CUR
      TYPE(type_info),              INTENT(IN)     :: T_INFO

      TYPE(ham_handle),             INTENT(IN)     :: HAMILTONIAN
      OVERLAP,                      INTENT(IN)     :: CQIJ(:, :, :, :)
      OVERLAP,                      INTENT(IN)     :: CDIJ(:, :, :, :)
      RGRID,                        INTENT(IN)     :: SV(:, :)
      TYPE(nonl_struct),            INTENT(INOUT)  :: NONL_S ! _test_mani must be mutable!
      TYPE(nonlr_struct),           INTENT(INOUT)  :: NONLR_S ! _test_mani must be mutable!
      TYPE(symmetry),               INTENT(IN)     :: SYMM
      REAL(q),                      INTENT(IN)     :: EXHF

      TYPE(in_struct),              INTENT(IN)     :: IO

      ! locals
      GDEF, ALLOCATABLE, TARGET, SAVE  :: HAM_DIFF(:, :, :)
      GDEF, ALLOCATABLE, SAVE          :: OVERLAP_DIFF(:, :, :)
      GDEF, ALLOCATABLE                :: TMP_MAT(:, :)
      GDEF, POINTER                    :: HAM_PTR(:, :, :, :)
      INTEGER                          :: NB_TOT
      INTEGER                          :: SPIN_IDX, NUM_SPIN
      TYPE(wavedes1)                   :: WDES1
      REAL(q), POINTER                 :: POSION_PTR(:,:)
      INTEGER                          :: LMDIM
      GDEF                             :: FIN_DIFF_FACTOR
      ! Needed by EDDIAG.
      INTEGER                          :: IFLAG
      LOGICAL                          :: LFIRST, LLAST

      INTEGER, SAVE                    :: CUR_DOF = 0
      INTEGER, SAVE                    :: CUR_DISPL = 0

      ! If nothing has been set up, stop
      IF (.NOT. THIS%IS_INITIALIZED) CALL vtutor%bug("EPI object not initialized!", __FILE__, __LINE__)

      ! Check if finite difference has made progress.
      IF (GLOB_FINDIF_INFO%NUM_DOF == 0) RETURN
      IF (GLOB_FINDIF_INFO%CUR_DOF == CUR_DOF .AND. GLOB_FINDIF_INFO%CUR_DISPL == CUR_DISPL) RETURN

      CUR_DOF = GLOB_FINDIF_INFO%CUR_DOF
      CUR_DISPL = GLOB_FINDIF_INFO%CUR_DISPL

      IF (CUR_DISPL /= 1 .AND. CUR_DISPL /= 2) CALL vtutor%bug( &
         "DISPL > 2 not supported! " // str(CUR_DISPL), __FILE__, __LINE__)

      NB_TOT = THIS%W_WAN_INIT%WDES%NB_TOT
      NUM_SPIN = THIS%W_WAN_INIT%WDES%ISPIN
      LMDIM = SIZE(CQIJ, 1)
      CALL SETWDES(THIS%W_WAN_INIT%WDES, WDES1, 1)

      ! Allocate the perturbed Hamiltonian and S-overlap matrices.
      IF (.NOT.ALLOCATED(THIS%COUPLING_WANN)) THEN
         ALLOCATE(THIS%COUPLING_WANN(NUM_SPIN), THIS%OVERLAP_PERT(NUM_SPIN))
         DO SPIN_IDX = 1, NUM_SPIN
            THIS%COUPLING_WANN(SPIN_IDX) = WANNIER_MAT_MODES_T(THIS%WANNIER_DESCRIPTOR, GLOB_FINDIF_INFO%NUM_DOF)
            THIS%OVERLAP_PERT(SPIN_IDX) = WANNIER_MAT_MODES_T(THIS%WANNIER_DESCRIPTOR, GLOB_FINDIF_INFO%NUM_DOF)
         ENDDO

         ALLOCATE(HAM_DIFF(NB_TOT, NB_TOT, NUM_SPIN))
         ALLOCATE(OVERLAP_DIFF(NB_TOT, NB_TOT, NUM_SPIN))
      ENDIF

      ALLOCATE(TMP_MAT(NB_TOT, NB_TOT))

      ! _test_mani doesn't work for real-space projectors yet!
      ! Project initial Wannier functions onto perturbed projectors.
      CALL PROALL(GRID, LATT_CUR, NONLR_S, NONL_S, THIS%W_WAN_INIT)

      !  ___________________________________
      ! |                                   |
      ! |  Calculate perturbed Hamiltonian  |
      ! |        < W_m | dH | W_n >         |
      ! |___________________________________|

      ! Set EDDIAG parameters
      IFLAG = 1
      IF (CUR_DISPL == 1) THEN
         ! First displacement +\Delta
         HAM_DIFF = 0
         LFIRST = .TRUE.
         LLAST = .FALSE.
      ELSE
         ! Second displacement -\Delta
         LFIRST = .FALSE.
         LLAST = .TRUE.
      ENDIF
      
      HAM_PTR(1 : NB_TOT, 1 : NB_TOT, 1 : 1, 1 : NUM_SPIN) => HAM_DIFF

      ! This calculates dH(-\Delta) - dH(+\Delta), so negative of what we want.
      ! The minus sign is accounted for in the GEMM call below.
      CALL EDDIAG(   HAMILTONIAN, GRID, LATT_CUR, NONLR_S, NONL_S, THIS%W_WAN_INIT, THIS%W_WAN_INIT%WDES, SYMM, &
                     LMDIM, CDIJ, CQIJ, IFLAG, SV, T_INFO, P, IO%IU0, EXHF, HAM_PTR, LFIRST, LLAST)

      !  ___________________________________
      ! |                                   |
      ! |  Calculate perturbed PAW overlap  |
      ! |        < W_m | dS | W_n >         |
      ! |___________________________________|

      DO SPIN_IDX = 1, NUM_SPIN
         IF (CUR_DISPL == 1) THEN
            ! First displacement +\Delta
            CALL CALC_PAW_OVERLAP(ELEMENTS(THIS%W_WAN_INIT, WDES1, SPIN_IDX), ELEMENTS(THIS%W_WAN_INIT, WDES1, SPIN_IDX), CQIJ, OVERLAP_DIFF(:, :, SPIN_IDX))
         ELSE
            ! Second displacement -\Delta
            CALL CALC_PAW_OVERLAP(ELEMENTS(THIS%W_WAN_INIT, WDES1, SPIN_IDX), ELEMENTS(THIS%W_WAN_INIT, WDES1, SPIN_IDX), CQIJ, TMP_MAT)
            OVERLAP_DIFF(:, :, SPIN_IDX) = OVERLAP_DIFF(:, :, SPIN_IDX) - TMP_MAT

            ! Create Wannier matrix object
            CALL THIS%OVERLAP_PERT(SPIN_IDX)%SET_MODE(OVERLAP_DIFF(1 : THIS%HAM_INIT(SPIN_IDX)%MAT_SIZE, 1 : THIS%HAM_INIT(SPIN_IDX)%MAT_SIZE, SPIN_IDX), CUR_DOF)
            THIS%OVERLAP_PERT(SPIN_IDX)%ELEMENTS(:, :, :, :, CUR_DOF) = &
            THIS%OVERLAP_PERT(SPIN_IDX)%ELEMENTS(:, :, :, :, CUR_DOF) / (2._q * GLOB_FINDIF_INFO%CUR_STEP)
         ENDIF
      ENDDO

      !  ___________________________________________
      ! |                                           |
      ! |  Calculate Wannier-space coupling matrix  |
      ! |         < W_m | dH - dS H | W_n >         |
      ! |___________________________________________|

      IF (CUR_DISPL == 2) THEN
         DO SPIN_IDX = 1, NUM_SPIN
            ASSOCIATE(MAT_SIZE => THIS%HAM_INIT(SPIN_IDX)%MAT_SIZE)
               ! We need initial Wannier-space Hamiltonian matrix
               CALL THIS%HAM_INIT(SPIN_IDX)%UNPACK(TMP_MAT(1 : MAT_SIZE, 1 : MAT_SIZE))

               ! Takes care of negative sign in perturbed Hamiltonian
               FIN_DIFF_FACTOR = -1._q / (2._q * GLOB_FINDIF_INFO%CUR_STEP)

               ! Calculate coupling
               CALL GGEMM( 'N', 'N', &
                  MAT_SIZE, MAT_SIZE, MAT_SIZE, FIN_DIFF_FACTOR, OVERLAP_DIFF(:, :, SPIN_IDX), NB_TOT, &
                  TMP_MAT, NB_TOT, FIN_DIFF_FACTOR, HAM_DIFF(:, :, SPIN_IDX), NB_TOT)

               ! Create Wannier matrix object
               CALL THIS%COUPLING_WANN(SPIN_IDX)%SET_MODE(HAM_DIFF(1 : MAT_SIZE, 1 : MAT_SIZE, SPIN_IDX), CUR_DOF)
            END ASSOCIATE
         ENDDO
      ENDIF

      ASSOCIATE(PC => THIS%WANNIER_DESCRIPTOR%PRIM_CELL)
         ! Set NONL_S to initial positions.
         POSION_PTR => NONL_S%POSION
         NULLIFY(NONL_S%POSION)
         ALLOCATE(NONL_S%POSION(3, PC%NUM_SUPER_ATOMS))
         NONL_S%POSION = PC%POSION_SUPER
         CALL PHASE(THIS%W_WAN_INIT%WDES, NONL_S, 0)

         ! Restore initial projectors.
         CALL PROALL(GRID, PC%SUPER_LATT, NONLR_S, NONL_S, THIS%W_WAN_INIT)

         ! Reset NONL_S to perturbed positions.
         DEALLOCATE(NONL_S%POSION)
         NONL_S%POSION => POSION_PTR
         NULLIFY(POSION_PTR)
         CALL PHASE(THIS%W_WAN_INIT%WDES, NONL_S, 0)
      END ASSOCIATE

      IF (CUR_DOF == GLOB_FINDIF_INFO%NUM_DOF .AND. CUR_DISPL == GLOB_FINDIF_INFO%NUM_DISPL) THEN
         THIS%IS_COMPLETED = .TRUE.
         DEALLOCATE(HAM_DIFF, OVERLAP_DIFF)
      ENDIF
   END SUBROUTINE

   SUBROUTINE EPI_PERTURB_WANNIER_AE( THIS, GRID, P, LATT_CUR, T_INFO, &
      HAMILTONIAN, CQIJ, CDIJ, SV, NONL_S, NONLR_S, SYMM, EXHF, IO)
      USE mgrid_struct_def,   ONLY: grid_3d
      USE pseudo_struct_def,  ONLY: potcar
      USE hamil_struct_def,   ONLY: ham_handle
      USE nonl_struct_def,    ONLY: nonl_struct
      USE nonlr_struct_def,   ONLY: nonlr_struct
      USE base,               ONLY: symmetry
      USE finite_differences, ONLY: GLOB_FINDIF_INFO
      USE wave,               ONLY: SETWDES
      USE nonl_high,          ONLY: PROALL
      USE wave_high,          ONLY: ELEMENTS
      USE subrot,             ONLY: EDDIAG
      USE nonl,               ONLY: PHASE
      USE wave_struct_def,    ONLY: wavefuna
      USE choleski,           ONLY: CALC_PAW_OVERLAP
      IMPLICIT NONE

      CLASS(ELPHON_INTERPOLATOR),   INTENT(INOUT)  :: THIS
      TYPE(grid_3d),                INTENT(IN)     :: GRID
      TYPE(potcar),                 INTENT(IN)     :: P(:)
      TYPE(latt),                   INTENT(IN)     :: LATT_CUR
      TYPE(type_info),              INTENT(IN)     :: T_INFO

      TYPE(ham_handle),             INTENT(IN)     :: HAMILTONIAN
      OVERLAP,                      INTENT(IN)     :: CQIJ(:, :, :, :)
      OVERLAP,                      INTENT(IN)     :: CDIJ(:, :, :, :)
      RGRID,                        INTENT(IN)     :: SV(:, :) ! local potential
      TYPE(nonl_struct),            INTENT(INOUT)  :: NONL_S ! _test_mani must be mutable!
      TYPE(nonlr_struct),           INTENT(INOUT)  :: NONLR_S ! _test_mani must be mutable!

      TYPE(in_struct),              INTENT(IN)     :: IO

      TYPE(symmetry),               INTENT(IN)     :: SYMM
      REAL(q),                      INTENT(IN)     :: EXHF

      ! locals
      GDEF, ALLOCATABLE, TARGET, SAVE  :: HAM_DIFF(:, :, :)
      GDEF, ALLOCATABLE, SAVE          :: CQIJ_PROJ_DIFF(:, :, :)
      GDEF, ALLOCATABLE, SAVE          :: CDIJ_PROJ_DIFF(:, :, :)
      GDEF, ALLOCATABLE, SAVE          :: D_CDIJ_DIFF(:, :, :)
      GDEF, ALLOCATABLE                :: TMP_MAT(:, :), TMP_MAT_2(:, :), TMP_TMP(:, :), TOTAL_MAT(:, :)
      GDEF, POINTER                    :: HAM_PTR(:, :, :, :)
      OVERLAP, ALLOCATABLE             :: CQIJ_PERT(:, :, :, :)
      INTEGER                          :: NB_TOT
      INTEGER                          :: SPIN_IDX, NUM_SPIN
      TYPE(wavedes1)                   :: WDES1
      REAL(q), POINTER                 :: POSION_PTR(:,:)
      INTEGER                          :: LMDIM
      INTEGER                          :: ATOM_IDX
      INTEGER                          :: MAT_SIZE
      GDEF                             :: FIN_DIFF_FACTOR
      TYPE(wavefuna)                   :: WA_INIT_PROJ
      TYPE(wavefuna)                   :: WA_PERT_PROJ

      ! Needed by EDDIAG.
      INTEGER                          :: IFLAG
      LOGICAL                          :: LFIRST, LLAST

      INTEGER, SAVE                    :: CUR_DOF = 0
      INTEGER, SAVE                    :: CUR_DISPL = 0

      ! new
      GDEF, ALLOCATABLE, TARGET  :: CPROJ_PERT(:, :, :, :)
      GDEF, POINTER              :: CPROJ_INIT_PTR(:, :, :, :)
      OVERLAP, ALLOCATABLE       :: CQDIJ_ZERO(:, :, :, :)

      ! If nothing has been set up, stop
      IF (.NOT. THIS%IS_INITIALIZED) CALL vtutor%bug("Elphon not initialized!", __FILE__, __LINE__)

      ! Check if finite difference has made progress.
      IF (GLOB_FINDIF_INFO%NUM_DOF == 0) RETURN
      IF (GLOB_FINDIF_INFO%CUR_DOF == CUR_DOF .AND. GLOB_FINDIF_INFO%CUR_DISPL == CUR_DISPL) RETURN

      CUR_DOF = GLOB_FINDIF_INFO%CUR_DOF
      CUR_DISPL = GLOB_FINDIF_INFO%CUR_DISPL

      IF (CUR_DISPL /= 1 .AND. CUR_DISPL /= 2) CALL vtutor%bug( &
         "DISPL > 2 not supported! " // str(CUR_DISPL), __FILE__, __LINE__)

      NB_TOT = THIS%W_WAN_INIT%WDES%NB_TOT
      NUM_SPIN = THIS%W_WAN_INIT%WDES%ISPIN
      LMDIM = SIZE(CQIJ, 1)
      MAT_SIZE = THIS%HAM_INIT(1)%MAT_SIZE

      CALL SETWDES(THIS%W_WAN_INIT%WDES, WDES1, 1)

      ! _test_mani clean up
      ! Allocate the perturbed Hamiltonian and stuff
      IF (.NOT.ALLOCATED(THIS%COUPLING_WANN)) THEN
         ALLOCATE(THIS%COUPLING_WANN(NUM_SPIN), THIS%OVERLAP_PERT(NUM_SPIN))
         DO SPIN_IDX = 1, NUM_SPIN
            THIS%COUPLING_WANN(SPIN_IDX) = WANNIER_MAT_MODES_T(THIS%WANNIER_DESCRIPTOR, GLOB_FINDIF_INFO%NUM_DOF)
            THIS%OVERLAP_PERT(SPIN_IDX) = WANNIER_MAT_MODES_T(THIS%WANNIER_DESCRIPTOR, GLOB_FINDIF_INFO%NUM_DOF)
         ENDDO

         ALLOCATE(HAM_DIFF(NB_TOT, NB_TOT, NUM_SPIN))
         ALLOCATE(CQIJ_PROJ_DIFF(MAT_SIZE, MAT_SIZE, NUM_SPIN))
         ALLOCATE(CDIJ_PROJ_DIFF(MAT_SIZE, MAT_SIZE, NUM_SPIN))
         ALLOCATE(D_CDIJ_DIFF(MAT_SIZE, MAT_SIZE, NUM_SPIN))
      ENDIF

      ALLOCATE(TMP_MAT(NB_TOT, NB_TOT))

      ALLOCATE(CPROJ_PERT(SIZE(THIS%W_WAN_INIT%CPROJ, 1), SIZE(THIS%W_WAN_INIT%CPROJ, 2), SIZE(THIS%W_WAN_INIT%CPROJ, 3), SIZE(THIS%W_WAN_INIT%CPROJ, 4)))

      ALLOCATE(CQDIJ_ZERO(SIZE(CQIJ, 1), SIZE(CQIJ, 2), SIZE(CQIJ, 3), SIZE(CQIJ, 4)))
      CQDIJ_ZERO = 0

      ! _test_mani doesn't work for real-space projectors yet!
      ! Project initial Wannier functions onto perturbed projectors
      CPROJ_INIT_PTR => THIS%W_WAN_INIT%CPROJ
      THIS%W_WAN_INIT%CPROJ => CPROJ_PERT
      CALL PROALL(GRID, LATT_CUR, NONLR_S, NONL_S, THIS%W_WAN_INIT)
      
      ! _TEST_MANI
      THIS%W_WAN_INIT%CPROJ => CPROJ_INIT_PTR

      !  ___________________________________
      ! |                                   |
      ! |  Calculate perturbed Hamiltonian  |
      ! |        < W_m | dH | W_n >         |
      ! |___________________________________|

      ! Set EDDIAG parameters
      IFLAG = 1
      IF (CUR_DISPL == 1) THEN
         ! First displacement +\Delta
         HAM_DIFF = 0
         LFIRST = .TRUE.
         LLAST = .FALSE.
      ELSE
         ! Second displacement -\Delta
         LFIRST = .FALSE.
         LLAST = .TRUE.

         ! _test_mani move when cleaning up
         ALLOCATE(TMP_MAT_2(MAT_SIZE, MAT_SIZE), TMP_TMP(MAT_SIZE, MAT_SIZE), TOTAL_MAT(MAT_SIZE, MAT_SIZE))
         IF (THIS%W_WAN_INIT%WDES%LOVERL) THEN
            ! Calc R
            CALL SET_CQIJ_PERT(T_INFO, P, NUM_SPIN, GLOB_FINDIF_INFO%CUR_ATOM_IDX, GLOB_FINDIF_INFO%CUR_CART_IDX, CQIJ_PERT)
         ELSE
            ALLOCATE(CQIJ_PERT(1, 1, 1, 1))
            CQIJ_PERT = 0
         ENDIF
      ENDIF
      
      ! Use rank retargeting so we can omit k index
      HAM_PTR(1 : NB_TOT, 1 : NB_TOT, 1 : 1, 1 : NUM_SPIN) => HAM_DIFF

      CALL EDDIAG(   HAMILTONIAN, GRID, LATT_CUR, NONLR_S, NONL_S, THIS%W_WAN_INIT, THIS%W_WAN_INIT%WDES, SYMM, &
                     LMDIM, CQDIJ_ZERO, CQDIJ_ZERO, IFLAG, SV, T_INFO, P, IO%IU0, EXHF, HAM_PTR, LFIRST, LLAST)

      DO SPIN_IDX = 1, NUM_SPIN
         ! Set WA with initial projectors
         THIS%W_WAN_INIT%CPROJ => CPROJ_INIT_PTR
         WA_INIT_PROJ = ELEMENTS(THIS%W_WAN_INIT, WDES1, SPIN_IDX)
         ! Set WA with perturbed projectors
         THIS%W_WAN_INIT%CPROJ => CPROJ_PERT
         WA_PERT_PROJ = ELEMENTS(THIS%W_WAN_INIT, WDES1, SPIN_IDX)
         
         IF (CUR_DISPL == 1) THEN
            CQIJ_PROJ_DIFF(:, :, SPIN_IDX) = 0
            D_CDIJ_DIFF(:, :, SPIN_IDX) = 0
            CDIJ_PROJ_DIFF(:, :, SPIN_IDX) = 0
         ENDIF

         ! _test_mani clean up
         FIN_DIFF_FACTOR = SIGN(1._q, 1.5_q - CUR_DISPL) / (2._q * GLOB_FINDIF_INFO%CUR_STEP)

         ! derivative of projector with CQIJ (init = pert) using finite differences
         CALL CALC_PAW_OVERLAP(WA_INIT_PROJ, WA_PERT_PROJ, CQIJ, TMP_MAT, ONLY_AUG=.TRUE.)
         CQIJ_PROJ_DIFF(:, :, SPIN_IDX) = CQIJ_PROJ_DIFF(:, :, SPIN_IDX) + TMP_MAT(1 : MAT_SIZE, 1 : MAT_SIZE) * FIN_DIFF_FACTOR

         ! derivative of projector with CDIJ (init) using finite differences
         ! ^ could be done in one call but we don't use this in final implementation, just testing
         CALL CALC_PAW_OVERLAP(WA_INIT_PROJ, WA_PERT_PROJ, THIS%CDIJ_INIT, TMP_MAT, ONLY_AUG=.TRUE.)
         CDIJ_PROJ_DIFF(:, :, SPIN_IDX) = CDIJ_PROJ_DIFF(:, :, SPIN_IDX) + TMP_MAT(1 : MAT_SIZE, 1 : MAT_SIZE) * FIN_DIFF_FACTOR

         ! derivative of CDIJ using finite differences
         CALL CALC_PAW_OVERLAP(WA_INIT_PROJ, WA_INIT_PROJ, CDIJ, TMP_MAT, ONLY_AUG=.TRUE.)
         D_CDIJ_DIFF(:, :, SPIN_IDX) = D_CDIJ_DIFF(:, :, SPIN_IDX) + TMP_MAT(1 : MAT_SIZE, 1 : MAT_SIZE) * FIN_DIFF_FACTOR

         ! Calc all the components (pure CHAOS at this point)
         ! _test_mani clean up once working!

         IF (CUR_DISPL == 2) THEN
            ! Hamiltonian (V + D + P(D part))
            ! ONLY V NOW
            HAM_DIFF(1 : MAT_SIZE, 1 : MAT_SIZE, SPIN_IDX) = HAM_DIFF(1 : MAT_SIZE, 1 : MAT_SIZE, SPIN_IDX) * FIN_DIFF_FACTOR
            ! Add V + D + P(D part) to total g
            ! ONLY V NOW
            TOTAL_MAT = HAM_DIFF(1 : MAT_SIZE, 1 : MAT_SIZE, SPIN_IDX)

            ! -D
            ! HAM_DIFF(1 : MAT_SIZE, 1 : MAT_SIZE, SPIN_IDX) = HAM_DIFF(1 : MAT_SIZE, 1 : MAT_SIZE, SPIN_IDX) - D_CDIJ_DIFF(:, :, SPIN_IDX)
            TOTAL_MAT = TOTAL_MAT + D_CDIJ_DIFF(:, :, SPIN_IDX)
            ! -P(D part)
            ! HAM_DIFF(1 : MAT_SIZE, 1 : MAT_SIZE, SPIN_IDX) = HAM_DIFF(1 : MAT_SIZE, 1 : MAT_SIZE, SPIN_IDX) - CDIJ_PROJ_DIFF(:, :, SPIN_IDX) - GCONJG(TRANSPOSE(CDIJ_PROJ_DIFF(:, :, SPIN_IDX)))
            TOTAL_MAT = TOTAL_MAT + CDIJ_PROJ_DIFF(:, :, SPIN_IDX) + GCONJG(TRANSPOSE(CDIJ_PROJ_DIFF(:, :, SPIN_IDX)))

            ! _test_mani
            ! Mul P(Q part) with eigenvalues (matmul Wannier Hamiltonian)
            CALL THIS%HAM_INIT(SPIN_IDX)%UNPACK(TMP_MAT_2)
            CALL GGEMM( 'N', 'N', &
               MAT_SIZE, MAT_SIZE, MAT_SIZE, one, TMP_MAT_2, MAT_SIZE, &
               CQIJ_PROJ_DIFF(:, :, SPIN_IDX), MAT_SIZE, zero, TMP_MAT, NB_TOT)
            TMP_MAT_2 = -TMP_MAT(1 : MAT_SIZE, 1 : MAT_SIZE) - GCONJG(TRANSPOSE(TMP_MAT(1 : MAT_SIZE, 1 : MAT_SIZE)))

            ! Add P(Q part) to total g
            TOTAL_MAT = TOTAL_MAT + TMP_MAT_2

            ! Calc R
            CALL CALC_PAW_OVERLAP(WA_INIT_PROJ, WA_INIT_PROJ, CQIJ_PERT, TMP_MAT, ONLY_AUG=.TRUE.)

            ! _test_mani
            ! Mul R with eigenvalues (matmul Wannier Hamiltonian)
            CALL THIS%HAM_INIT(SPIN_IDX)%UNPACK(TMP_MAT_2)
            CALL GGEMM( 'N', 'N', &
               MAT_SIZE, MAT_SIZE, MAT_SIZE, one, TMP_MAT_2, MAT_SIZE, &
               TMP_MAT, NB_TOT, zero, TMP_TMP, MAT_SIZE)
            ! Mul from other side -> transpose (remember, anti hermitian)
            TMP_MAT_2 = -TMP_TMP - TRANSPOSE(TMP_TMP)
            TMP_MAT_2 = -TMP_MAT_2

            ! Add R to total g
            TOTAL_MAT = TOTAL_MAT + TMP_MAT_2

            ! Create Wannier matrix object (total g)
            CALL THIS%COUPLING_WANN(SPIN_IDX)%SET_MODE(TOTAL_MAT, CUR_DOF)
            ! Init and set overlap part to 0 so we can use the ZPR routine with this
            TOTAL_MAT = 0
            CALL THIS%OVERLAP_PERT(SPIN_IDX)%SET_MODE(TOTAL_MAT, CUR_DOF)
         ENDIF
      ENDDO

      ! Point to original memory
      THIS%W_WAN_INIT%CPROJ => CPROJ_INIT_PTR

      IF (CUR_DOF == GLOB_FINDIF_INFO%NUM_DOF .AND. CUR_DISPL == GLOB_FINDIF_INFO%NUM_DISPL) THEN
         THIS%IS_COMPLETED = .TRUE.
         ! Deallocate the SAVE arrays
         ! _test_mani add more after clean up
         DEALLOCATE(HAM_DIFF, CQIJ_PROJ_DIFF, CDIJ_PROJ_DIFF, D_CDIJ_DIFF)
      ENDIF
   END SUBROUTINE

   SUBROUTINE SET_CQIJ_PERT(T_INFO, P, NUM_SPIN, ATOM_IDX, CART_IDX, CQIJ_PERT)
      USE poscar, ONLY  : type_info
      USE pseudo, ONLY  : potcar
      IMPLICIT NONE

      TYPE(type_info),        INTENT(IN)  :: T_INFO
      TYPE(potcar),           INTENT(IN)  :: P(:)
      INTEGER,                INTENT(IN)  :: NUM_SPIN
      INTEGER,                INTENT(IN)  :: ATOM_IDX
      INTEGER,                INTENT(IN)  :: CART_IDX
      OVERLAP, ALLOCATABLE,   INTENT(OUT) :: CQIJ_PERT(:, :, :, :)

      INTEGER  :: LMDIM
      INTEGER  :: SPIN_IDX

      LMDIM = MAXVAL(P(:)%LMMAX)

      ALLOCATE(CQIJ_PERT(LMDIM, LMDIM, T_INFO%NIONS, NUM_SPIN))
      CQIJ_PERT = 0

      DO SPIN_IDX = 1, NUM_SPIN
         CQIJ_PERT(:, :, ATOM_IDX, SPIN_IDX) = -P(T_INFO%ITYP(ATOM_IDX))%NABLA(CART_IDX, :, :)
      ENDDO
   END SUBROUTINE SET_CQIJ_PERT

   ! _test_mani move somewhere else
   SUBROUTINE CALC_NABIJ_2(P)
      USE pseudo_struct_def,  ONLY: potcar
      USE asa,                ONLY: SETYLM_NABLA_YLM
      USE radial,             ONLY: GRAD_, & ! _test_mani use GRAD instead maybe?
                                    SIMPI
      IMPLICIT NONE

      TYPE(potcar), TARGET,   INTENT(INOUT)  :: P(:)

      INTEGER, EXTERNAL       :: MAXL1

      INTEGER                 :: TYPE_IDX
      INTEGER                 :: CART_IDX
      INTEGER                 :: NUM_R_POINTS, R_IDX
      INTEGER                 :: NUM_CHANNELS, CH1, CH2
      INTEGER                 :: MAX_L
      INTEGER                 :: L1, L2, M1, M2
      INTEGER                 :: LM1, LM2, LLM1, LLM2
      REAL(q)                 :: SUM
      REAL(q), ALLOCATABLE    :: DWAE(:), DWPS(:), TMP(:)
      REAL(q), ALLOCATABLE    :: YLM_NABLA_YLM(:,:,:), YLM_X_YLM(:,:,:)
      TYPE(potcar), POINTER   :: PP

      DO TYPE_IDX = 1, SIZE(P)
         PP => P(TYPE_IDX)
         NUM_CHANNELS = PP%LMAX
         MAX_L = MAXL1(PP)
         NUM_R_POINTS = PP%R%NMAX

         PP%NABLA = 0

         ALLOCATE(DWAE(NUM_R_POINTS), DWPS(NUM_R_POINTS), TMP(NUM_R_POINTS))

         ! ALLOCATE and setup nabla array
         ALLOCATE(YLM_NABLA_YLM((MAX_L+1)*(MAX_L+1), (MAX_L+1)*(MAX_L+1), 0:3))
         ALLOCATE(YLM_X_YLM((MAX_L+1)*(MAX_L+1), (MAX_L+1)*(MAX_L+1), 0:3))

         CALL SETYLM_NABLA_YLM(MAX_L, YLM_NABLA_YLM, YLM_X_YLM)

         LM2 = 1

         DO CH2 = 1, NUM_CHANNELS
            LM1 = 1
            L2 = PP%LPS(CH2)

            TMP = PP%WAE(:, CH2) / PP%R%R
            CALL GRAD_(PP%R, TMP, DWAE)
            DWAE = DWAE * PP%R%R

            TMP = PP%WPS(:, CH2) / PP%R%R
            CALL GRAD_(PP%R, TMP, DWPS)
            DWPS = DWPS * PP%R%R

            DO CH1 = 1, NUM_CHANNELS
               L1 = PP%LPS(CH1)

               DO M2 = 1, 2*L2 + 1
                  LLM2 = L2*L2 + M2

                  DO M1 = 1, 2*L1 + 1
                     LLM1 = L1*L1 + M1

                     DO CART_IDX = 1, 3

                        IF (ABS(YLM_X_YLM(LLM1, LLM2, CART_IDX)) > 1E-8_q .OR. &
                           ABS(YLM_NABLA_YLM(LLM1, LLM2, CART_IDX)) > 1E-8_q) THEN

                           DO R_IDX = 1, NUM_R_POINTS
                              TMP(R_IDX) = &
                                 PP%WAE(R_IDX, CH1) * ( DWAE(R_IDX) * YLM_X_YLM(LLM1, LLM2, CART_IDX) +         &
                                 YLM_NABLA_YLM(LLM1, LLM2, CART_IDX) / PP%R%R(R_IDX) * PP%WAE(R_IDX, CH2) ) -  &
                                 PP%WPS(R_IDX, CH1) * ( DWPS(R_IDX) * YLM_X_YLM(LLM1, LLM2, CART_IDX) +         &
                                 YLM_NABLA_YLM(LLM1, LLM2, CART_IDX) / PP%R%R(R_IDX) * PP%WPS(R_IDX, CH2) )
                           ENDDO

                           CALL SIMPI(PP%R,TMP,SUM)

                           PP%NABLA(CART_IDX, LM1+M1-1, LM2+M2-1) = &
                           PP%NABLA(CART_IDX, LM1+M1-1, LM2+M2-1) + SUM
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
               LM1 = LM1 + 2*L1 + 1
            ENDDO
            LM2 = LM2 + 2*L2 + 1
         ENDDO
         
         DEALLOCATE(YLM_NABLA_YLM, YLM_X_YLM)
         DEALLOCATE(DWAE, DWPS, TMP)

         DO LM2 = 1, PP%LMMAX
            DO LM1 = LM2, PP%LMMAX
               PP%NABLA(:, LM1, LM2) = (PP%NABLA(:, LM1, LM2) - PP%NABLA(:, LM2, LM1)) / 2
               PP%NABLA(:, LM2, LM1) = -PP%NABLA(:, LM1, LM2)
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE CALC_NABIJ_2

   ! _test_mani this might cause the original wavespin object (W) to become dysfunctional
   SUBROUTINE EPI_DEALLOC(THIS)
      USE wave,   ONLY: DEALLOCW
      IMPLICIT NONE

      CLASS(ELPHON_INTERPOLATOR),   INTENT(INOUT)  :: THIS

      CALL DEALLOCW(THIS%W_WAN_INIT)
      IF (ALLOCATED(THIS%HAM_INIT)) DEALLOCATE(THIS%HAM_INIT)
      IF (ALLOCATED(THIS%COUPLING_WANN)) DEALLOCATE(THIS%COUPLING_WANN)
      IF (ALLOCATED(THIS%OVERLAP_PERT)) DEALLOCATE(THIS%OVERLAP_PERT)
      IF (ALLOCATED(THIS%CDIJ_INIT)) DEALLOCATE(THIS%CDIJ_INIT)
      IF (ALLOCATED(THIS%COUPLING_WANN_LONG)) DEALLOCATE(THIS%COUPLING_WANN_LONG)
      IF (ALLOCATED(THIS%EWALD_G_VECTORS)) DEALLOCATE(THIS%EWALD_G_VECTORS)

      THIS%IS_INITIALIZED = .FALSE.
      THIS%IS_COMPLETED = .FALSE.
   END SUBROUTINE

   SUBROUTINE EPI_FINALIZE(THIS)
      IMPLICIT NONE
      TYPE(ELPHON_INTERPOLATOR), INTENT(INOUT)  :: THIS
      CALL THIS%DEALLOC()
   END SUBROUTINE


!                              ___________________
!_____________________________/ EPI_CALC_COUPLING \_____________________________
!
!> @brief Calculates the PAW electron-phonon matrix element.
!>
!> Calculates the quantity
!> \f[
!>    \tilde{g}_{mn,\nu}\qty(\vb k, \vb q) =
!>    \braket{
!>       \tilde{\Psi}_{m \vb k + \vb q} |
!>       \partial_{\nu \vb q} \tilde{H} -
!>       \varepsilon_{n \vb k} \partial_{\nu \vb q} \tilde{S} |
!>       \tilde{\Psi}_{n \vb k}}
!> \f]
!>
!> See Phys. Rev. B 101, 184302 for details
!>
!> @param[in,out] THIS  Instance reference
!> @param[out] COUPLING Calculated electron-phonon matrix
!> @param[in] Q_VEC     Phonon wave vector in direct coordinates of the
!>                      Brillouin zone of the primitive cell
!> @param[in] K_VEC     Electron Bloch vector in direct coordinates of the
!>                      Brillouin zone of the primitive cell
!> @param[in] SPIN      Spin index
!> @param[in] T_INFO    Atom-type information
!_______________________________________________________________________________
   SUBROUTINE EPI_CALC_COUPLING(THIS, COUPLING, Q_VEC, K_VEC, SPIN, T_INFO)
      IMPLICIT NONE

      CLASS(ELPHON_INTERPOLATOR),      INTENT(INOUT)  :: THIS
      COMPLEX(q), ALLOCATABLE, TARGET, INTENT(OUT)    :: COUPLING(:, :, :)
      REAL(q),                         INTENT(IN)     :: Q_VEC(3), K_VEC(3)
      INTEGER,                         INTENT(IN)     :: SPIN
      TYPE(type_info),                 INTENT(IN)     :: T_INFO

      COMPLEX(q), POINTER  :: COUPLING_PTR(:, :)

      ! This should never happen
      IF (.NOT. THIS%IS_COMPLETED) CALL vtutor%bug("EPI_CALC_COUPLING: Missing Wannier-space data!", __FILE__, __LINE__)

      ! Calculate coupling without rotating vibrational degrees of freedom
      CALL THIS%CALC_COUPLING_ATOMS(COUPLING, Q_VEC, K_VEC, SPIN)

      ! Rotate coupling with phonon eigenmodes
      COUPLING_PTR(1 : SIZE(COUPLING, 1) * SIZE(COUPLING, 2), 1 : SIZE(COUPLING, 3)) => COUPLING
      CALL ROTATE_WITH_MODES(COUPLING_PTR, Q_VEC, THIS%WANNIER_DESCRIPTOR%PRIM_CELL, T_INFO)
   END SUBROUTINE


!                           _________________________
!__________________________/ EPI_CALC_COUPLING_ATOMS \__________________________
!
!> @brief Calculates an intermediary PAW electron-phonon matrix element.
!>
!> Calculates the quantity
!> \f[
!>    \tilde{g}_{mn,\kappa \alpha}\qty(\vb k, \vb q) =
!>    \braket{
!>       \tilde{\Psi}_{m \vb k + \vb q} |
!>       \partial_{\kappa \alpha} \tilde{H} -
!>       \varepsilon_{n \vb k} \partial_{\kappa \alpha} \tilde{S} |
!>       \tilde{\Psi}_{n \vb k}}
!> \f]
!>
!> In order to obtain the electron-phonon matrix element in the phonon-mode
!> representation, this quantity still needs to be transformed with the matrix
!> of phonon modes.
!> 
!> See Phys. Rev. B 101, 184302 for details
!>
!> @param[in,out] THIS  Instance reference
!> @param[out] COUPLING Calculated electron-phonon matrix
!> @param[in] Q_VEC     Phonon wave vector in direct coordinates of the
!>                      Brillouin zone of the primitive cell
!> @param[in] K_VEC     Electron Bloch vector in direct coordinates of the
!>                      Brillouin zone of the primitive cell
!> @param[in] SPIN      Spin index
!_______________________________________________________________________________
   SUBROUTINE EPI_CALC_COUPLING_ATOMS(THIS, COUPLING, Q_VEC, K_VEC, SPIN)
      IMPLICIT NONE

      CLASS(ELPHON_INTERPOLATOR),   INTENT(INOUT)  :: THIS
      COMPLEX(q), ALLOCATABLE,      INTENT(OUT)    :: COUPLING(:, :, :)
      REAL(q),                      INTENT(IN)     :: Q_VEC(3), K_VEC(3)
      INTEGER,                      INTENT(IN)     :: SPIN

      ! This should never happen.
      IF (.NOT. THIS%IS_COMPLETED) CALL vtutor%bug("EPI_CALC_COUPLING_ATOMS: Missing Wannier-space data!", __FILE__, __LINE__)

      CALL THIS%COUPLING_WANN(SPIN)%TO_BLOCH(COUPLING, THIS%HAM_INIT(SPIN), Q_VEC, K_VEC)

      IF (THIS%DO_POLAR) THEN
         CALL THIS%SUBTRACT_COUPLING_WANN_LONG(COUPLING, Q_VEC, K_VEC, SPIN)
         CALL THIS%ADD_COUPLING_LONG(COUPLING, Q_VEC, K_VEC, SPIN)
      ENDIF
   END SUBROUTINE

!                              ___________________
!_____________________________/ ROTATE_WITH_MODES \_____________________________
!
!> @brief Rotates a matrix with the rescaled phonon-mode matrix.
!>
!> \f[ A'_{ab} = \sum_c A_{ac} u_{cb} \f]
!>
!> This is used to transform a Wannier-space matrix from an atomic-
!> displacement picture to a phonon-displacement picture. The index \f$c\f$ is
!> a composite index for primitive-atom index \f$\kappa\f$ and Cartesian
!> direction \f$\alpha\f$. \f$u_{cb}\f$ is the rescaled phonon-mode matrix.
!>
!> @param[in,out] ARRAY The array to be transformed
!> @param[in] Q_VEC     Phonon wave vector in direct coordinates of the
!>                      Brillouin zone of the primitive cell
!> @param[in] PRIM_CELL Primitive-cell information
!> @param[in] T_INFO    Atom-type information
!_______________________________________________________________________________
   SUBROUTINE ROTATE_WITH_MODES(ARRAY, Q_VEC, PRIM_CELL, T_INFO)
      USE phonon,    ONLY: PHON_CALC_PHONONS, &
                           PHON_FIND_TRANSLATIONAL_MODES
      USE constant
      IMPLICIT NONE

      COMPLEX(q),          INTENT(INOUT)  :: ARRAY(:, :)
      REAL(q),             INTENT(IN)     :: Q_VEC(3)
      TYPE(PRIM_CELL_T),   INTENT(IN)     :: PRIM_CELL
      TYPE(type_info),     INTENT(IN)     :: T_INFO

      INTEGER                 :: NUM_MODES
      INTEGER                 :: MODE_IDX
      INTEGER                 :: P_IDX
      REAL(q), ALLOCATABLE    :: FREQU(:)
      COMPLEX(q), ALLOCATABLE :: MODES(:, :)
      INTEGER                 :: TRANSLATIONAL_MODES(3)

      REAL(q), PARAMETER      :: HBAR = 1.0545718E-24_q / SQRT(EVTOJ) / SQRT(AMTOKG)

      ! Calculate phonons at q-point
      CALL PHON_CALC_PHONONS(Q_VEC, PRIM_CELL, T_INFO, FREQU, MODES)

      NUM_MODES = SIZE(MODES, 1)

      IF (ALL(ABS(Q_VEC) < 1e-10_q)) THEN
         TRANSLATIONAL_MODES = PHON_FIND_TRANSLATIONAL_MODES(MODES)
      ELSE
         TRANSLATIONAL_MODES = 0
      ENDIF

      ! Rescale the phonon eigenvectors by the masses
      DO P_IDX = 1, PRIM_CELL%NUM_PRIM_ATOMS
         MODES(P_IDX * 3 - 2 : P_IDX * 3, :) = MODES(P_IDX * 3 - 2 : P_IDX * 3, :) / &
         SQRT(T_INFO%POMASS(T_INFO%ITYP(PRIM_CELL%SUPER_IDX(1, P_IDX))))
      ENDDO

      ARRAY = MATMUL(ARRAY, MODES)

      ! Finally, rescale with frequency and attach correct units
      DO MODE_IDX = 1, NUM_MODES
         IF (ANY(MODE_IDX == TRANSLATIONAL_MODES)) THEN
            ! For translational modes at Gamma, set 0
            ARRAY(:, MODE_IDX) = 0._q
         ELSE IF (FREQU(MODE_IDX) < 0) THEN
            ! _test_mani
            ! For negative frequency, take complex SQRT, for now
            ARRAY(:, MODE_IDX) = ARRAY(:, MODE_IDX) * &
            SQRT(CMPLX(HBAR / (2 * FREQU(MODE_IDX)), 0, q))
         ELSE
            ! Prefactor to obtain eV
            ARRAY(:, MODE_IDX) = ARRAY(:, MODE_IDX) * &
            SQRT(HBAR / (2 * FREQU(MODE_IDX)))
         ENDIF
      ENDDO
   END SUBROUTINE ROTATE_WITH_MODES


!                              ___________________
!_____________________________/ ELPH_OVERLAP_WSWQ \_____________________________
!
!> @brief Computes the S-overlap between wavefuntions W and the wavefunctions
!> W^Q supplied in file WAVECAR.qqq.
!>
!> If all-electron corrections should be added to the matrix element,
!> additionally, the file POSCAR.qqq needs to be supplied that contains the
!> atomic configuration corresponding to the wavefunctions W^Q.
!>
!> @param W             Current wavefunctions
!> @param LATT_CUR      Current lattice
!> @param CQIJ          Augmentation charges
!> @param NONLR_S       Non-local structure (real space)
!> @param NONL_S        Non-local structure (reciprocal space)
!> @param P             Structure holding pseudo potential/PAW information
!> @param T_INFO        Structure holding atom type/position information
!> @param ALL_ELECTRON  Include all-electron corrections?
!> @param IO            IO unit info
!_______________________________________________________________________________
   SUBROUTINE ELPH_OVERLAP_WSWQ(W, LATT_CUR, CQIJ, NONLR_S, NONL_S, P, T_INFO, ALL_ELECTRON, IO)
      USE nonl_struct_def,    ONLY: nonl_struct
      USE nonlr_struct_def,   ONLY: nonlr_struct
      USE pseudo_struct_def,  ONLY: potcar
      USE nonl_high,          ONLY: PROALL
      USE nonl,               ONLY: PHASE
      USE wave,               ONLY: SETWDES, &
                                    DEALLOCW
      USE wave_high,          ONLY: ELEMENTS
      USE fileio,             ONLY: INWAV_ALTERNATIVE
      USE POSCAR,             ONLY: RD_POSCAR_ALT
      USE poscar_struct_def,  ONLY: dynamics
      USE choleski,           ONLY: CALC_PAW_OVERLAP
      IMPLICIT NONE

      TYPE(wavespin),      INTENT(INOUT)  :: W        ! must be mutable but is unchanged on exit
      TYPE(latt),          INTENT(IN)     :: LATT_CUR
      OVERLAP,             INTENT(IN)     :: CQIJ(:, :, :, :)
      TYPE(nonlr_struct),  INTENT(INOUT)  :: NONLR_S  ! must be mutable but is unchanged on exit
      TYPE(nonl_struct),   INTENT(INOUT)  :: NONL_S   ! must be mutable but is unchanged on exit
      TYPE(potcar),        INTENT(INOUT)  :: P(:)     ! NABLA array is allocated and setup
      TYPE(type_info),     INTENT(IN)     :: T_INFO
      LOGICAL,             INTENT(IN)     :: ALL_ELECTRON
      TYPE(in_struct),     INTENT(IN)     :: IO

      TYPE (wavespin)   :: WQ
      TYPE (wavedes1)   :: WDES1
      GDEF, ALLOCATABLE :: COVL(:,:)
      INTEGER           :: SPIN_IDX, K_IDX
      INTEGER           :: M, N
      INTEGER           :: ATOM_IDX, CART_IDX
      LOGICAL           :: LWREAD
      TYPE(LATT)        :: LATT_Q
      TYPE(type_info)   :: T_INFO_Q
      TYPE(dynamics)    :: DYN_Q
      INTEGER           :: UNIT

      REAL(q), POINTER     :: POSION_PTR(:, :)
      OVERLAP, ALLOCATABLE :: CRIJ(:, :, :, :)
      REAL(q)              :: DISP(3)

      CALL INWAV_ALTERNATIVE(IO, W%WDES, WQ, W%WDES%GRID, LATT_CUR, LWREAD, 'qqq')
      IF (.NOT.LWREAD) CALL vtutor%error('ELPH_OVERLAP_WSWQ: Unable to read file "WAVECAR.qqq"')

      IF (ALL_ELECTRON) THEN
         ! Read displaced structure from file
         CALL RD_POSCAR_ALT("POSCAR.qqq", LATT_Q, T_INFO_Q, DYN_Q, IO%IU0, IO%IU6)

         ! Set NONL_S to displaced positions
         POSION_PTR => NONL_S%POSION
         NONL_S%POSION => T_INFO_Q%POSION
         CALL PHASE(W%WDES, NONL_S, 0)

         ! Overlap between perturbed projectors |p^Q> and perturbed pseudo orbitals |W^Q>
         CALL PROALL(W%WDES%GRID, LATT_Q, NONLR_S, NONL_S, WQ)

         ! Reset NONL_S to equilibrium positions
         NONL_S%POSION => POSION_PTR
         CALL PHASE(W%WDES, NONL_S, 0)


         ! Calculate derivative of partial waves
         CALL CALC_NABIJ_2(P)

         ALLOCATE(CRIJ(SIZE(CQIJ, 1), SIZE(CQIJ, 2), SIZE(CQIJ, 3), SIZE(CQIJ, 4)))
         CRIJ = CQIJ

         ! Multiply with displacement pattern
         DO SPIN_IDX = 1, W%WDES%ISPIN
            DO ATOM_IDX = 1, SIZE(T_INFO_Q%POSION, 2)
               ! Get displacement vector in direct coordinates w.r.t. unperturbed supercell
               DISP = MATMUL(MATMUL(LATT_Q%A, T_INFO_Q%POSION(:, ATOM_IDX)), LATT_CUR%B) - T_INFO%POSION(:, ATOM_IDX)
               ! Bring desplacement vector to WS cell
               DISP = MODULO(DISP + 0.5_q, 1._q) - 0.5_q
               ! Transform to Cartesian coordinates
               DISP = MATMUL(LATT_CUR%A, DISP)

               DO CART_IDX = 1, 3
                  CRIJ(:, :, ATOM_IDX, SPIN_IDX) = &
                  CRIJ(:, :, ATOM_IDX, SPIN_IDX) + &
                  P(T_INFO%ITYP(ATOM_IDX))%NABLA(CART_IDX, :, :) * DISP(CART_IDX)
               ENDDO
            ENDDO
         ENDDO

         ! _test_mani maybe deallocate NABLA again if it was allocated?
      ELSE
         ! Overlap between initial projectors |p> and perturbed pseudo orbitals |W^Q>
         CALL PROALL(W%WDES%GRID, LATT_CUR, NONLR_S, NONL_S, WQ)
      ENDIF

      ! Storage for overlap matrices
      ALLOCATE(COVL(W%WDES%NB_TOT, W%WDES%NB_TOT))

      ! Generate output file
      IF (IO%IU0 >= 0) OPEN(NEWUNIT=UNIT, FILE='WSWQ', STATUS='REPLACE')

      DO SPIN_IDX = 1, W%WDES%ISPIN
         DO K_IDX = 1, W%WDES%NKPTS
            CALL SETWDES(W%WDES, WDES1, K_IDX)

            IF (ALL_ELECTRON) THEN
               ! <W|p> <f| |f^Q> <p^Q|W^Q>
               CALL CALC_PAW_OVERLAP(ELEMENTS(W, WDES1, SPIN_IDX), ELEMENTS(WQ, WDES1, SPIN_IDX), CRIJ, COVL)
            ELSE
               ! <W|S|W^Q>
               CALL CALC_PAW_OVERLAP(ELEMENTS(W, WDES1, SPIN_IDX), ELEMENTS(WQ, WDES1, SPIN_IDX), CQIJ, COVL)
            ENDIF

            ! Write transition matrix to file
            IF (IO%IU0 >= 0) THEN
               WRITE(UNIT, '(A, I1, A, I6)') '   spin=', SPIN_IDX, ', kpoint=', K_IDX

               WRITE(UNIT, '(A, I6, A, I6, A, 2F16.9) ') (('i=', M, ', j=', N, ' :', CMPLX(COVL(M, N), KIND=q), &
                  M = 1, W%WDES%NB_TOT), N = 1, W%WDES%NB_TOT)
            ENDIF
         ENDDO
      ENDDO

      IF (IO%IU0 >= 0) CLOSE(UNIT)

      ! Deallocate perturbed wavefunction
      CALL DEALLOCW(WQ)

      ! _test_mani make sure type_info and dynamics don't leak memory and are properly deallocated
      ! see issue #222
   END SUBROUTINE ELPH_OVERLAP_WSWQ


!                               ________________
!______________________________/ EPI_INIT_POLAR \_______________________________
!
!> @brief Initializes the polar electron-phonon code.
!>
!> The long-ranged polar electron-phonon matrix elements are setup in the
!> Wannier representation and stored for later use.
!>
!> @param[in,out] THIS     Instance reference
!> @param[in] SETTINGS     Electron-phonon settings
!> @param[in] POLAR_DATA   Born effective charges and dielectric tensor
!_______________________________________________________________________________
   SUBROUTINE EPI_INIT_POLAR(THIS, SETTINGS, POLAR_DATA)
      IMPLICIT NONE

      CLASS(ELPHON_INTERPOLATOR),   INTENT(INOUT)  :: THIS
      TYPE(ELPHON_SETTINGS),        INTENT(IN)     :: SETTINGS
      TYPE(POLAR_DATA_T),           INTENT(IN)     :: POLAR_DATA

      REAL(q)  :: ABS_G_CUTOFF

      IF (.NOT. THIS%IS_COMPLETED) CALL vtutor%bug("Electron-phonon code not completed!", __FILE__, __LINE__)

      ABS_G_CUTOFF = GET_G_CUTOFF(SETTINGS%POLAR_G_CUTOFF, THIS%WANNIER_DESCRIPTOR%PRIM_CELL%PRIM_LATT)
      CALL THIS%POLAR_EWALD%INIT(ABS_G_CUTOFF, SETTINGS%POLAR_EWALD_THRESHOLD, POLAR_DATA)
      CALL GENERATE_EWALD_G_VECTORS(THIS%POLAR_EWALD, THIS%WANNIER_DESCRIPTOR%PRIM_CELL%PRIM_LATT, THIS%EWALD_G_VECTORS)
      CALL THIS%SETUP_COUPLING_WANN_LONG(SETTINGS)
      THIS%DO_POLAR = .TRUE.
   END SUBROUTINE


!                        ______________________________
!_______________________/ EPI_SETUP_COUPLING_WANN_LONG \________________________
!
!> @brief Sets up the long-ranged polar electron-phonon matrix elements in the
!> Wannier representation.
!>
!> Performs a truncated Ewald sum (reciprocal part only) to calculate the long-
!> ranged part of the dipole contribution to the electron-phonon matrix in the
!> supercell Wannier representation.
!> It is stored for later use during the Wannier interpolation.
!>
!> @param[in,out] THIS     Instance reference
!> @param[in] SETTINGS     Electron-phonon settings
!_______________________________________________________________________________
   SUBROUTINE EPI_SETUP_COUPLING_WANN_LONG(THIS, SETTINGS)
      USE wannier_mats, ONLY: SETUP_PHASE_BUFFER
      USE constant
      IMPLICIT NONE

      CLASS(ELPHON_INTERPOLATOR),   INTENT(INOUT)  :: THIS
      TYPE(ELPHON_SETTINGS),        INTENT(IN)     :: SETTINGS

      REAL(q), ALLOCATABLE    :: G_VECTORS(:, :)
      REAL(q)                 :: G_CART(3)
      INTEGER                 :: G_IDX
      COMPLEX(q), ALLOCATABLE :: PHASE(:, :)
      INTEGER                 :: SPIN_IDX
      INTEGER                 :: P_IDX
      INTEGER                 :: MODE_IDX
      INTEGER                 :: CART_IDX

      ASSOCIATE( &
         PC => THIS%WANNIER_DESCRIPTOR%PRIM_CELL, &
         N_SPIN => SIZE(THIS%COUPLING_WANN), &
         N_PRIM => THIS%WANNIER_DESCRIPTOR%PRIM_CELL%NUM_PRIM_ATOMS, &
         N_CELLS => THIS%WANNIER_DESCRIPTOR%PRIM_CELL%NUM_CELLS, &
         N_ORBS => THIS%WANNIER_DESCRIPTOR%NUM_PRIM_ORBITALS, &
         EWALD => THIS%POLAR_EWALD)

         CALL GENERATE_EWALD_G_VECTORS(EWALD, PC%SUPER_LATT, G_VECTORS)

         ALLOCATE(THIS%COUPLING_WANN_LONG(N_ORBS, N_CELLS, N_PRIM * 3, N_SPIN))
         THIS%COUPLING_WANN_LONG = 0

         DO SPIN_IDX = 1, N_SPIN
            ! Loop over all G-vectors inside the Ewald cut-off except Gamma!
            DO G_IDX = 2, SIZE(G_VECTORS, 2)
               ! Already in Carthesian coordinates
               G_CART = G_VECTORS(:, G_IDX)

               ! Skip if Ewald suppression factor is too small
               IF (.NOT. EWALD%TRY_CALC(G_CART)) CYCLE

               DO P_IDX = 1, N_PRIM
                  CALL SETUP_PHASE_BUFFER(PHASE, THIS%WANNIER_DESCRIPTOR, P_IDX, G_CART)
                  PHASE = PHASE * PC%MIN_IMAGES(1, P_IDX, 1)%CALC_PRIM_PHASE(-G_CART)

                  DO CART_IDX = 1, 3
                     MODE_IDX = P_IDX * 3 - 3 + CART_IDX
                     THIS%COUPLING_WANN_LONG(:, :, MODE_IDX, SPIN_IDX) = &
                     THIS%COUPLING_WANN_LONG(:, :, MODE_IDX, SPIN_IDX) + &
                     EWALD%BORN_FACTOR(CART_IDX, P_IDX) * EWALD%DIELEC_FACTOR * EWALD%EXP_FACTOR * PHASE(:, :)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         THIS%COUPLING_WANN_LONG = THIS%COUPLING_WANN_LONG * (0._q, 1._q) * 4._q * PI * FELECT / (PC%SUPER_LATT%OMEGA)
      END ASSOCIATE
   END SUBROUTINE


!                       _________________________________
!______________________/ EPI_SUBTRACT_COUPLING_WANN_LONG \______________________
!
!> @brief Subtracts the long-ranged polar contribution from the electron-phonon
!> matrix.
!>
!> The previously setup long-ranged part is Wannier interpolated to the
!> primitive cell and subtracted from the supplied matrix. Note, that the
!> vibrational indices of the electron-phonon matrix must be in the atomic-
!> displacement picture for this to work (they are not yet rotated by the
!> phonon modes).
!>
!> We perform this interpolation by choice since it is the best trade of
!> between convenience of implemenation and computational speed.
!> This interpolation is a lot faster than the usual electron-phonon
!> interpolation since only a single sum over cells is used.
!>
!> @param[in] THIS         Instance reference
!> @param[in] SETTINGS     Electron-phonon settings
!> @param[in,out] COUPLING Electron-phonon matrix
!> @param[in] Q_VEC        Phonon wave vector in direct coordinates of the
!>                         Brillouin zone of the primitive cell
!> @param[in] K_VEC        Electron Bloch vector in direct coordinates of the
!>                         Brillouin zone of the primitive cell
!> @param[in] SPIN         Spin index
!_______________________________________________________________________________
   SUBROUTINE EPI_SUBTRACT_COUPLING_WANN_LONG(THIS, COUPLING, Q_VEC, K_VEC, SPIN)
      USE wannier_mats, ONLY: SETUP_PHASE_BUFFER
      IMPLICIT NONE

      CLASS(ELPHON_INTERPOLATOR),   INTENT(IN)     :: THIS
      COMPLEX(q),                   INTENT(INOUT)  :: COUPLING(:, :, :)
      REAL(q),                      INTENT(IN)     :: Q_VEC(3)
      REAL(q),                      INTENT(IN)     :: K_VEC(3)
      INTEGER,                      INTENT(IN)     :: SPIN

      COMPLEX(q), ALLOCATABLE             :: UK_D(:, :), UKQ_D(:, :)
      REAL(q), ALLOCATABLE                :: TMP_VEC(:)
      COMPLEX(q), ALLOCATABLE             :: BUFFER(:)
      COMPLEX(q), ALLOCATABLE             :: PHASE(:, :)
      INTEGER                             :: ORB_IDX_1, ORB_IDX_2
      INTEGER                             :: CELL_IDX
      INTEGER                             :: P_IDX
      INTEGER                             :: MODE_IDX
      INTEGER                             :: CART_IDX
      REAL(q)                             :: Q_CART(3)

      ASSOCIATE( &
         N_ORBS => THIS%WANNIER_DESCRIPTOR%NUM_PRIM_ORBITALS, &
         N_PRIM => THIS%WANNIER_DESCRIPTOR%PRIM_CELL%NUM_PRIM_ATOMS)

         ! Calculate the 'electronic' rotation matrices from Wannier to Bloch basis
         CALL THIS%HAM_INIT(SPIN)%DIAG_AT_K(K_VEC, TMP_VEC, UK_D)
         CALL THIS%HAM_INIT(SPIN)%DIAG_AT_K(K_VEC + Q_VEC, TMP_VEC, UKQ_D)

         ALLOCATE(BUFFER(N_ORBS))

         Q_CART = MATMUL(THIS%WANNIER_DESCRIPTOR%PRIM_CELL%PRIM_LATT%B, Q_VEC)

         DO P_IDX = 1, N_PRIM
            ! Pre-calculate phase factors
            CALL SETUP_PHASE_BUFFER(PHASE, THIS%WANNIER_DESCRIPTOR, P_IDX, -Q_CART)

            DO CART_IDX = 1, 3
               MODE_IDX = P_IDX * 3 - 3 + CART_IDX
               BUFFER = SUM(THIS%COUPLING_WANN_LONG(:, :, MODE_IDX, SPIN) * PHASE(:, :), 2)

               ! Rotate with Wannier rotation matrices
               ! U_{ma,k+q} g_{a} U^\dagger_{an,k}
               DO ORB_IDX_2 = 1, N_ORBS
                  DO ORB_IDX_1 = 1, N_ORBS
                     COUPLING(ORB_IDX_1, ORB_IDX_2, MODE_IDX) = &
                     COUPLING(ORB_IDX_1, ORB_IDX_2, MODE_IDX) - &
                     DOT_PRODUCT(UKQ_D(:, ORB_IDX_1), BUFFER(:) * UK_D(:, ORB_IDX_2))
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      END ASSOCIATE
   END SUBROUTINE


!                            _______________________
!___________________________/ EPI_ADD_COUPLING_LONG \___________________________
!
!> @brief Adds the long-ranged polar contribution to the electron-phonon matrix.
!>
!> The long-ranged part is calculated from an analytic expression and added to
!> the supplied Bloch-space matrix. Note, that the vibrational indices of the
!> electron- phonon matrix must be in the atomic-displacement picture for this
!> to work (they are not yet rotated by the phonon modes).
!>
!> @param[in,out] THIS     Instance reference
!> @param[in,out] COUPLING Electron-phonon matrix
!> @param[in] Q_VEC        Phonon wave vector in direct coordinates of the
!>                         Brillouin zone of the primitive cell
!> @param[in] K_VEC        Electron Bloch vector in direct coordinates of the
!>                         Brillouin zone of the primitive cell
!> @param[in] SPIN         Spin index
!_______________________________________________________________________________
   SUBROUTINE EPI_ADD_COUPLING_LONG(THIS, COUPLING, Q_VEC, K_VEC, SPIN)
      USE constant
      IMPLICIT NONE

      CLASS(ELPHON_INTERPOLATOR),   INTENT(INOUT)  :: THIS
      COMPLEX(q),                   INTENT(INOUT)  :: COUPLING(:, :, :)
      REAL(q),                      INTENT(IN)     :: Q_VEC(3)
      REAL(q),                      INTENT(IN)     :: K_VEC(3)
      INTEGER,                      INTENT(IN)     :: SPIN

      REAL(q)                    :: Q_CART(3), QG_CART(3)
      INTEGER                    :: G_IDX
      COMPLEX(q)                 :: PHASE
      COMPLEX(q), ALLOCATABLE    :: QG_FACTOR(:)
      INTEGER                    :: P_IDX
      INTEGER                    :: MODE_IDX
      INTEGER                    :: ORB_IDX_1, ORB_IDX_2
      COMPLEX(q), ALLOCATABLE    :: UKQ_D(:, :), UK_D(:, :)
      REAL(q), ALLOCATABLE       :: TMP_VEC(:)


      ASSOCIATE( &
         PC => THIS%WANNIER_DESCRIPTOR%PRIM_CELL, &
         N_ORBS => THIS%WANNIER_DESCRIPTOR%NUM_PRIM_ORBITALS, &
         N_MODES => THIS%COUPLING_WANN(SPIN)%NUM_MODES, &
         N_PRIM => THIS%WANNIER_DESCRIPTOR%PRIM_CELL%NUM_PRIM_ATOMS, &
         EWALD => THIS%POLAR_EWALD)

         ALLOCATE(QG_FACTOR(N_MODES))
         QG_FACTOR = 0

         ! Convert to Carthesian coordinates.
         Q_CART = MATMUL(PC%PRIM_LATT%B, Q_VEC)

         ! Loop over all G-vector inside the Ewald cut-off
         DO G_IDX = 1, SIZE(THIS%EWALD_G_VECTORS, 2) ! ! _test_mani
            ! G already in Carthesian coordinates
            QG_CART = Q_CART + THIS%EWALD_G_VECTORS(:, G_IDX)

            ! Only take non-zero q + G.
            IF (ALL(ABS(QG_CART) < 1E-10_q)) CYCLE

            ! Skip if Ewald suppression factor is too small
            IF (.NOT. EWALD%TRY_CALC(QG_CART)) CYCLE

            DO P_IDX = 1, N_PRIM
               ! Compute phase factor
               PHASE = PC%MIN_IMAGES(1, P_IDX, 1)%CALC_PRIM_PHASE(-QG_CART)

               ! Calculate prefactor
               QG_FACTOR(P_IDX * 3 - 2 : P_IDX * 3) = &
               QG_FACTOR(P_IDX * 3 - 2 : P_IDX * 3) + &
               EWALD%DIELEC_FACTOR * EWALD%BORN_FACTOR(:, P_IDX) * EWALD%EXP_FACTOR * PHASE
            ENDDO
         ENDDO

         QG_FACTOR = QG_FACTOR * &
         (0._q, 1._q) * 4._q * PI * FELECT / (PC%PRIM_LATT%OMEGA)

         CALL THIS%HAM_INIT(SPIN)%DIAG_AT_K(K_VEC, TMP_VEC, UK_D)
         CALL THIS%HAM_INIT(SPIN)%DIAG_AT_K(K_VEC + Q_VEC, TMP_VEC, UKQ_D)

         ! Rotate with Wannier rotation matrices
         ! U_{ma,k+q} U^\dagger_{an,k}
         DO MODE_IDX = 1, N_MODES
            DO ORB_IDX_2 = 1, N_ORBS
               DO ORB_IDX_1 = 1, N_ORBS
                  COUPLING(ORB_IDX_1, ORB_IDX_2, MODE_IDX) = &
                  COUPLING(ORB_IDX_1, ORB_IDX_2, MODE_IDX) + &
                  DOT_PRODUCT(UKQ_D(:, ORB_IDX_1), UK_D(:, ORB_IDX_2)) * QG_FACTOR(MODE_IDX)
               ENDDO
            ENDDO
         ENDDO
      END ASSOCIATE
   END SUBROUTINE


! _test_mani HELPER STUFF, MAINLY DEBUG!!

   SUBROUTINE AVERAGE_DEGEN_Z1(SPECTRUM, DIMEN, ARRAY)
      IMPLICIT NONE

      REAL(q),    INTENT(IN)     :: SPECTRUM(:)
      INTEGER,    INTENT(IN)     :: DIMEN
      COMPLEX(q), INTENT(INOUT)  :: ARRAY(:)

      ! locals
      INTEGER                    :: NUM_VALUES
      INTEGER, ALLOCATABLE       :: DEGEN_REF(:)
      REAL(q), ALLOCATABLE       :: TEMP_SPEC(:)
      INTEGER                    :: I, J

      NUM_VALUES = SIZE(ARRAY, DIMEN)

      ALLOCATE(DEGEN_REF(NUM_VALUES), TEMP_SPEC(NUM_VALUES))
      DEGEN_REF = 0

      TEMP_SPEC = SPECTRUM(1 : NUM_VALUES)
      TEMP_SPEC = SCALE(TEMP_SPEC, -EXPONENT(MAXVAL(ABS(TEMP_SPEC))))
      DO I = 1, NUM_VALUES
         IF (DEGEN_REF(I) == 0) THEN
            DO J = I + 1, NUM_VALUES
               IF (DEGEN_REF(J) == 0 .AND. ABS(TEMP_SPEC(I) - TEMP_SPEC(J)) < 1E-5_q) THEN
                  DEGEN_REF(I) = I
                  DEGEN_REF(J) = I
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      DO I = 1, NUM_VALUES
         IF (DEGEN_REF(I) == I) THEN
            ! Sum up all degenerate ARRAY values.
            DO J = I + 1, NUM_VALUES
               IF (DEGEN_REF(J) == I) THEN
                  ARRAY(I) = ARRAY(I) + ARRAY(J)
               ENDIF
            ENDDO
            ! Calculate mean value.
            ARRAY(I) = ARRAY(I) / COUNT(DEGEN_REF == I)
            ! Copy for each degenerate ARRAY value.
            DO J = I + 1, NUM_VALUES
               IF (DEGEN_REF(J) == I) THEN
                  ARRAY(J) = ARRAY(I)
               ENDIF
            ENDDO
         ENDIF
      ENDDO
   END SUBROUTINE AVERAGE_DEGEN_Z1


   SUBROUTINE AVERAGE_DEGEN_Z2(SPECTRUM, DIMEN, ARRAY)
      IMPLICIT NONE

      REAL(q),    INTENT(IN)     :: SPECTRUM(:)
      INTEGER,    INTENT(IN)     :: DIMEN
      COMPLEX(q), INTENT(INOUT)  :: ARRAY(:, :)

      ! locals
      INTEGER                    :: NUM_VALUES
      INTEGER, ALLOCATABLE       :: DEGEN_REF(:)
      REAL(q), ALLOCATABLE       :: TEMP_SPEC(:)
      INTEGER                    :: I, J

      NUM_VALUES = SIZE(ARRAY, DIMEN)

      ! IF (SIZE(SPECTRUM, 1) < NUM_VALUES) THEN
      !    CALL vtutor%error("Internal error in AVERAGE_DEGEN_Z1: Spectrum too small for specified &
      !       &matrix dimension: " // str(SIZE(SPECTRUM,1)) // " < " // str(NUM_VALUES) // ".")
      ! ENDIF

      ALLOCATE(DEGEN_REF(NUM_VALUES), TEMP_SPEC(NUM_VALUES))
      DEGEN_REF = 0

      TEMP_SPEC = SPECTRUM(1 : NUM_VALUES)
      TEMP_SPEC = SCALE(TEMP_SPEC, -EXPONENT(MAXVAL(ABS(TEMP_SPEC))))
      DO I = 1, NUM_VALUES
         IF (DEGEN_REF(I) == 0) THEN
            DO J = I + 1, NUM_VALUES
               IF (DEGEN_REF(J) == 0 .AND. ABS(TEMP_SPEC(I) - TEMP_SPEC(J)) < 1E-5_q) THEN
                  DEGEN_REF(I) = I
                  DEGEN_REF(J) = I
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      SELECT CASE(DIMEN)
      CASE(1)
         DO I = 1, NUM_VALUES
            IF (DEGEN_REF(I) == I) THEN
               ! Sum up all degenerate ARRAY values.
               DO J = I + 1, NUM_VALUES
                  IF (DEGEN_REF(J) == I) THEN
                     ARRAY(I, :) = ARRAY(I, :) + ARRAY(J, :)
                  ENDIF
               ENDDO
               ! Calculate mean value.
               ARRAY(I, :) = ARRAY(I, :) / COUNT(DEGEN_REF == I)
               ! Copy for each degenerate ARRAY value.
               DO J = I + 1, NUM_VALUES
                  IF (DEGEN_REF(J) == I) THEN
                     ARRAY(J, :) = ARRAY(I, :)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      CASE(2)
         DO I = 1, NUM_VALUES
            IF (DEGEN_REF(I) == I) THEN
               ! Sum up all degenerate ARRAY values.
               DO J = I + 1, NUM_VALUES
                  IF (DEGEN_REF(J) == I) THEN
                     ARRAY(:, I) = ARRAY(:, I) + ARRAY(:, J)
                  ENDIF
               ENDDO
               ! Calculate mean value.
               ARRAY(:, I) = ARRAY(:, I) / COUNT(DEGEN_REF == I)
               ! Copy for each degenerate ARRAY value.
               DO J = I + 1, NUM_VALUES
                  IF (DEGEN_REF(J) == I) THEN
                     ARRAY(:, J) = ARRAY(:, I)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      END SELECT
   END SUBROUTINE AVERAGE_DEGEN_Z2


   SUBROUTINE AVERAGE_DEGEN_Z3(SPECTRUM, DIMEN, ARRAY)
      IMPLICIT NONE

      REAL(q),    INTENT(IN)     :: SPECTRUM(:)
      INTEGER,    INTENT(IN)     :: DIMEN
      COMPLEX(q), INTENT(INOUT)  :: ARRAY(:, :, :)

      ! locals
      INTEGER                    :: NUM_VALUES
      INTEGER, ALLOCATABLE       :: DEGEN_REF(:)
      REAL(q), ALLOCATABLE       :: TEMP_SPEC(:)
      INTEGER                    :: I, J

      NUM_VALUES = SIZE(ARRAY, DIMEN)

      ! IF (SIZE(SPECTRUM, 1) < NUM_VALUES) THEN
      !    CALL vtutor%error("Internal error in AVERAGE_DEGEN_Z2: Spectrum too small for specified &
      !       &matrix dimension: " // str(SIZE(SPECTRUM,1)) // " < " // str(NUM_VALUES) // ".")
      ! ENDIF

      ALLOCATE(DEGEN_REF(NUM_VALUES), TEMP_SPEC(NUM_VALUES))
      DEGEN_REF = 0

      TEMP_SPEC = SPECTRUM(1 : NUM_VALUES)
      TEMP_SPEC = SCALE(TEMP_SPEC, -EXPONENT(MAXVAL(ABS(TEMP_SPEC))))
      DO I = 1, NUM_VALUES
         IF (DEGEN_REF(I) == 0) THEN
            DO J = I + 1, NUM_VALUES
               IF (DEGEN_REF(J) == 0 .AND. ABS(TEMP_SPEC(I) - TEMP_SPEC(J)) < 1E-5_q) THEN
                  DEGEN_REF(I) = I
                  DEGEN_REF(J) = I
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      SELECT CASE(DIMEN)
      CASE(1)
         DO I = 1, NUM_VALUES
            IF (DEGEN_REF(I) == I) THEN
               ! Sum up all degenerate ARRAY values.
               DO J = I + 1, NUM_VALUES
                  IF (DEGEN_REF(J) == I) THEN
                     ARRAY(I, :, :) = ARRAY(I, :, :) + ARRAY(J, :, :)
                  ENDIF
               ENDDO
               ! Calculate mean value.
               ARRAY(I, :, :) = ARRAY(I, :, :) / COUNT(DEGEN_REF == I)
               ! Copy for each degenerate ARRAY value.
               DO J = I + 1, NUM_VALUES
                  IF (DEGEN_REF(J) == I) THEN
                     ARRAY(J, :, :) = ARRAY(I, :, :)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      CASE(2)
         DO I = 1, NUM_VALUES
            IF (DEGEN_REF(I) == I) THEN
               ! Sum up all degenerate ARRAY values.
               DO J = I + 1, NUM_VALUES
                  IF (DEGEN_REF(J) == I) THEN
                     ARRAY(:, I, :) = ARRAY(:, I, :) + ARRAY(:, J, :)
                  ENDIF
               ENDDO
               ! Calculate mean value.
               ARRAY(:, I, :) = ARRAY(:, I, :) / COUNT(DEGEN_REF == I)
               ! Copy for each degenerate ARRAY value.
               DO J = I + 1, NUM_VALUES
                  IF (DEGEN_REF(J) == I) THEN
                     ARRAY(:, J, :) = ARRAY(:, I, :)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      CASE(3)
         DO I = 1, NUM_VALUES
            IF (DEGEN_REF(I) == I) THEN
               ! Sum up all degenerate ARRAY values.
               DO J = I + 1, NUM_VALUES
                  IF (DEGEN_REF(J) == I) THEN
                     ARRAY(:, :, I) = ARRAY(:, :, I) + ARRAY(:, :, J)
                  ENDIF
               ENDDO
               ! Calculate mean value.
               ARRAY(:, :, I) = ARRAY(:, :, I) / COUNT(DEGEN_REF == I)
               ! Copy for each degenerate ARRAY value.
               DO J = I + 1, NUM_VALUES
                  IF (DEGEN_REF(J) == I) THEN
                     ARRAY(:, :, J) = ARRAY(:, :, I)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      END SELECT
   END SUBROUTINE AVERAGE_DEGEN_Z3










! --------------------------------------
! --------------------------------------
! --------------------------------------
! --------------------------------------
! _test_mani
! Brillouin stuff! -> should either be moved, removed or merged.

!                                    ______
!  _________________________________/ init \__________________________________
! |                                                                           |
! |  Acts as a class constructor and must be called before any other          |
! |  routine.                                                                 |
! |___________________________________________________________________________|
! |                                                                           |
! |  lattice - Structure holding the lattice vectors.                         |
! |___________________________________________________________________________|
   subroutine init(this, lattice)
      implicit none

      class(brill_zone_t), intent(inout)  :: this
      type(latt),          intent(in)     :: lattice

      this%lattice%a = lattice%a
      call lattic(this%lattice)

      call this%generate_mid_points()
   END SUBROUTINE init

!                                 ____________
!  ______________________________/ destructor \_______________________________
! |                                                                           |
! |  The class destructor.                                                    |
! |___________________________________________________________________________|
   subroutine destructor(this)
      implicit none

      type(brill_zone_t), intent(inout)  :: this

      if (allocated(this%qmesh)) deallocate(this%qmesh)
      if (allocated(this%mid_point)) deallocate(this%mid_point)
      if (allocated(this%mid_sqr)) deallocate(this%mid_sqr)
   END SUBROUTINE destructor

!                             _____________________
!  __________________________/ generate_mid_points \__________________________
! |                                                                           |
! |  Generates points mid-way between neighboring reciprocal-lattice points.  |
! |___________________________________________________________________________|
   subroutine generate_mid_points(this)
      implicit none

      class(brill_zone_t), intent(inout)  :: this

      integer  :: x, y, z
      integer  :: num

      if (.not.allocated(this%mid_point)) allocate(this%mid_point(3, num_neighbors))
      if (.not.allocated(this%mid_sqr)) allocate(this%mid_sqr(num_neighbors))

      ! Generate all the neighboring mid-points.
      num = 0
      do x = -1, 1
      do y = -1, 1
      do z = -1, 1
         if (x /= 0 .or. y /= 0 .or. z /= 0) then
            num = num + 1
            this%mid_point(:, num) = matmul(this%lattice%b, [x, y, z]) / 2._q
            this%mid_sqr(num) = dot_product(this%mid_point(:, num), this%mid_point(:, num))
         endif
      enddo
      enddo
      enddo
   END SUBROUTINE generate_mid_points

!                                 _____________
!  ______________________________/ check_point \______________________________
! |                                                                           |
! |  This routine checks wether a given point in reciprocal space lies        |
! |  within the first BZ, on its boundary or outside.                         |
! |___________________________________________________________________________|
! |                                                                           |
! |  point     - Reciprocal point to check in Cartesian coordinates.          |
! |  is_inside - Returns true only if the point lies inside or on the         |
! |              boundary.                                                    |
! |  on_bounds - Returns true only if the point lies on the boundary.         |
! |___________________________________________________________________________|
   subroutine check_point(this, point, is_inside, on_bounds)
      implicit none

      class(brill_zone_t), intent(in)  :: this
      real(q),             intent(in)  :: point(3)
      logical,             intent(out) :: is_inside, on_bounds

      integer  :: j
      real(q)  :: dist_sqr

      is_inside = .true.
      on_bounds = .false.

      do j = 1, num_neighbors
         dist_sqr = dot_product(point, this%mid_point(:, j))

         ! Check if the point is inside the bz.
         if (this%mid_sqr(j) + brill_tiny < dist_sqr) then
            is_inside = .false.
            on_bounds = .false.
            exit
         endif

         ! Check if the point is on the bz boundaries.
         if (this%mid_sqr(j) - brill_tiny > dist_sqr) then
            on_bounds = .true.
         endif
      enddo
   END SUBROUTINE check_point

!                               _________________
!  ____________________________/ calculate_shift \____________________________
! |                                                                           |
! |  This routine calculates by how much the coordinate system of the         |
! |  generating lattice needs to be shifted in order to contain the entire    |
! |  first BZ of the target lattice within the parallelepiped spanned by the  |
! |  generator vectors.                                                       |
! |  This way, a regular sampling grid of the first BZ can more easily be     |
! |  described by a 3D index array.                                           |
! |___________________________________________________________________________|
! |                                                                           |
! |  gen_latt - The generating lattice used for sampling.                     |
! |  shift    - Integer Vector describing the actual shift as an integer      |
! |             multiple of the reciprocal generator lattice vectors.         |
! |___________________________________________________________________________|
   subroutine calculate_shift(this, gen_latt, shift)
      implicit none

      class(brill_zone_t), intent(in)  :: this
      type(latt),          intent(in)  :: gen_latt
      integer,             intent(out) :: shift(3)

      integer  :: i, j, k, x
      real(q)  :: equ_mat(3, 3), equ_vec(3)
      real(q)  :: mat_b(3, 3)
      real(q)  :: point(3)
      real(q)  :: det, inv_det
      logical  :: is_inside, on_bounds

      shift = 0
      
      ! Loop over all distinct triplets of neighbor points.
      ! In each sub-loop calculate part of the equation system.
      do i = 1, num_neighbors
         equ_mat(:, 3) = this%mid_point(:, i)
         equ_vec(3) = this%mid_sqr(i)
         do j = i + 1, num_neighbors
            equ_mat(:, 2) = this%mid_point(:, j) - this%mid_point(:, i)
            equ_vec(2) = this%mid_sqr(j) - this%mid_sqr(i)
            do k = j + 1, num_neighbors
               equ_mat(:, 1) = this%mid_point(:, k) - this%mid_point(:, j)
               equ_vec(1) = this%mid_sqr(k) - this%mid_sqr(j)

               ! Check if the equation system has a solution.
               det = det_3x3(equ_mat)
               if (abs(det) < brill_tiny) cycle
               inv_det = 1._q / det

               ! Solve the system using Cramer's rule.
               do x = 1, 3
                  mat_b = transpose(equ_mat)
                  mat_b(:, x) = equ_vec
                  point(x) = det_3x3(mat_b) * inv_det
               enddo

               ! We are only interested in points on the boundaries.
               call this%check_point(point, is_inside, on_bounds)
               if (.not.on_bounds) cycle

               ! Express the point in direct coordinates of the generators.
               point = matmul(point, gen_latt%b)

               ! Get the maximally negative coordinate per axis, rounded down.
               do x = 1, 3
                  shift(x) = min(shift(x), floor(point(x)))
               enddo
            enddo
         enddo
      enddo
   END SUBROUTINE calculate_shift


!                              ___________________
!  ___________________________/ sample_brill_zone \___________________________
! |                                                                           |
! |  This routine samples the first Brillouin zone using a uniform grid       |
! |  spanned by three linearly independent reciprocal lattice vectors         |
! |  (generator lattice).                                                     |
! |___________________________________________________________________________|
! |                                                                           |
! |  generators - Matrix containing the generator vectors. Same format as     |
! |               latt%a.                                                     |
! |  subdivide  - Subdivides the generator lattice along each axis, e.g.      |
! |               subdivision by 2 means 8 times as many samples.             |
! |___________________________________________________________________________|
   subroutine sample_brill_zone(this, generators, subdivide)
      use tutor, only: vtutor
      implicit none

      class(brill_zone_t), intent(inout)  :: this
      real(q),             intent(in)     :: generators(3, 3)
      integer,             intent(in)     :: subdivide

      ! locals
      type(latt)           :: gen_latt
      integer              :: shift(3)
      logical, allocatable :: check_grid(:, :, :)
      integer              :: grid_size(3)
      integer              :: x, y, z
      integer              :: i
      real(q)              :: point(3), tmp_point(3)
      logical              :: inside_bz, on_bounds
      integer              :: int_point(3)
      integer              :: num_points


real(q) :: subdiv_b(3, 3)
real(q) :: tmp_r
integer :: tmp_i
! Subdivide generator vectors.
! Lattice 'b' is used for the inverse transform.
subdiv_b = generators / subdivide
grid_size = huge(1)

do i = 1, 3
   tmp_point = matmul(subdiv_b(:, i), this%lattice%a)
   do x = 1, 3
      if (tmp_point(x) == 0) cycle
      tmp_r = abs(1._q / tmp_point(x))
      tmp_i = nint(tmp_r)
      if (abs(tmp_r - tmp_i) < 1e-5 .and. tmp_i < grid_size(i)) grid_size(i) = tmp_i
   ENDDO
enddo

if (any(grid_size == huge(1))) then
   call vtutor%stopCode()
endif

if (allocated(this%qmesh)) deallocate(this%qmesh)
allocate(this%qmesh(3, product(grid_size)))
i = 0

do z = 0, grid_size(3) - 1
   do y = 0, grid_size(2) - 1
      do x = 0, grid_size(1) - 1
         i = i + 1
         this%qmesh(:, i) = matmul(subdiv_b, [real(q):: x, y, z])
      enddo
   enddo
enddo

this%qmesh = matmul(transpose(this%lattice%a), this%qmesh)
this%qmesh = modulo(this%qmesh + .5_q, 1._q) - .5_q

return

      ! Subdivide generator vectors.
      ! Lattice 'b' is used for the inverse transform.
      gen_latt%a = generators / subdivide
      call lattic(gen_latt)

      ! Shift the generators so the whole BZ fits inside the positive 'octant'.
      call this%calculate_shift(gen_latt, shift)
      grid_size = abs(shift) * 2 + 1

      ! Use a bool grid to mark all points inside the BZ.
      allocate(check_grid(grid_size(1), grid_size(2), grid_size(3)))
      check_grid = .false.

      num_points = 0

      do z = 1, grid_size(1)
         do y = 1, grid_size(2)
            do x = 1, grid_size(3)
               ! If this point is a duplicate, skip.
               if (check_grid(x, y, z)) then
                  check_grid(x, y, z) = .false.
                  cycle
               endif

               point = matmul(gen_latt%a, [x, y, z] + shift - 1)
               call this%check_point(point, inside_bz, on_bounds)

               if (.not.inside_bz) cycle
               
               ! Add the point to the grid.
               check_grid(x, y, z) = .true.
               num_points = num_points + 1

               ! If the point is on the boundaries, search for duplicates.
               if (on_bounds) then
                  do i = 1, num_neighbors
                     tmp_point = point + 2 * this%mid_point(:, i)

                     call this%check_point(tmp_point, inside_bz, on_bounds)

                     ! Duplicates must lie on the boundaries.
                     if (.not.on_bounds) cycle

                     ! Get the point in direct generator coordinates.
                     tmp_point = matmul(tmp_point, gen_latt%b) - shift + 1._q
                     int_point = nint(tmp_point)

                     ! If the point is comensurable with the grid, mark it.
                     if (all(abs(tmp_point - int_point) < brill_tiny)) then
                        check_grid(int_point(1), int_point(2), int_point(3)) = .true.
                     endif
                  enddo
               endif
            enddo
         enddo
      enddo

      ! Allocate the q-point mesh to the correct size.
      if (allocated(this%qmesh)) deallocate(this%qmesh)
      allocate(this%qmesh(3, num_points))

      i = 0

      ! Fill the mesh with all points.
      do z = 1, grid_size(1)
         do y = 1, grid_size(2)
            do x = 1, grid_size(3)
               if (check_grid(x, y, z)) then
                  i = i + 1
                  this%qmesh(:, i) = matmul(gen_latt%a, [x, y, z] + shift - 1)
               endif
            enddo
         enddo
      enddo

      ! Convert q-points to direct coordinates w.r.t. reciprocal primitive cell.
      this%qmesh = matmul(transpose(this%lattice%a), this%qmesh)

      deallocate(check_grid)

   END SUBROUTINE sample_brill_zone


!                                 _____________
!  ______________________________/ write_qmesh \______________________________
! |                                                                           |
! |  This routine writes the generated q-point mesh to a file in readable     |
! |  form, e.g. for visualization purposes.                                   |
! |___________________________________________________________________________|
   subroutine write_qmesh(this)
      implicit none

      class(brill_zone_t), intent(in)  :: this

      integer  :: i
      real(q)  :: point(3)

      ! Output the q-mesh to a file for plotting etc.
      open(unit = 77767, file = 'debug_qmesh', status = 'replace')
      do i = 1, size(this%qmesh, 2)
         point = matmul(this%lattice%b, this%qmesh(:, i))
         write(77767, '(i6, a, f12.6, a, f12.6, a, f12.6)') &
      &  i, ',  ', point(1), ', ', point(2), ', ', point(3)
      enddo
      close(77767)
   END SUBROUTINE write_qmesh


!                              ___________________
!  ___________________________/ map_to_brill_zone \___________________________
! |                                                                           |
! |  This function maps any reciprocal-space point to an equivalent point     |
! |  inside the first Brillouin zone by shifting along reciprocal lattice     |
! |  vectors.                                                                 |
! |___________________________________________________________________________|
! |                                                                           |
! |  lattice     - Lattice structure used to describe the first Brillouin     |
! |                zone.                                                      |
! |  cart_qpoint - Point in reciprocal space in Cartesian coordinates.        |
! |  min_vec     - On return, point inside the first Brillouin zone.          |
! |___________________________________________________________________________|
   function map_to_brill_zone(lattice, cart_qpoint) result(min_vec)
      implicit none

      type(latt), intent(in)  :: lattice
      real(q), intent(in)     :: cart_qpoint(3)
      real(q)                 :: min_vec(3)

      ! locals
      integer  :: x, y, z
      real(q)  :: central_point(3)
      real(q)  :: test_vec(3)
      real(q)  :: test_sqr, min_sqr

      central_point = matmul(cart_qpoint, lattice%a)
      central_point = modulo(central_point + .5_q, 1._q) - .5_q

      ! Find cell containing q-point.
      min_sqr = huge(1._q)
      do x = -1, 1
      do y = -1, 1
      do z = -1, 1
         test_vec = matmul(lattice%b, central_point + [x, y, z])
         test_sqr = dot_product(test_vec, test_vec)
         if (test_sqr < min_sqr) then
            min_sqr = test_sqr
            min_vec = test_vec
         endif
      enddo
      enddo
      enddo
   endfunction map_to_brill_zone

!                                 _____________
!  ______________________________/ map_to_zone \______________________________
! |                                                                           |
! |  This function maps any reciprocal-space point to an equivalent point     |
! |  inside the first Brillouin zone by shifting along reciprocal lattice     |
! |  vectors.                                                                 |
! |___________________________________________________________________________|
! |                                                                           |
! |  cart_qpoint - Point in reciprocal space in Cartesian coordinates.        |
! |  min_vec     - On return, point inside the first Brillouin zone of the    |
! |                stored lattice.                                            |
! |___________________________________________________________________________|
   function map_to_zone(this, cart_qpoint) result(min_vec)
      implicit none

      class(brill_zone_t), intent(in)  :: this
      real(q), intent(in)              :: cart_qpoint(3)
      real(q)                          :: min_vec(3)

      min_vec = map_to_brill_zone(this%lattice, cart_qpoint)
   endfunction map_to_zone

END MODULE elphon
