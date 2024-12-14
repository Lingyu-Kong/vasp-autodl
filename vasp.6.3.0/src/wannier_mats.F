#include "symbol.inc"

MODULE wannier_mats
   USE prec
   USE elphon_util
   USE msymmetry, ONLY: MIN_IMAGE_T, &
                        PRIM_CELL_T
   USE tutor,     ONLY: vtutor
   USE string,    ONLY: str
   PRIVATE

   PUBLIC   WANNIER_DESCRIPTOR_T, WANNIER_MAT_SYMM_T, WANNIER_MAT_MODES_T, &
            SETUP_PHASE_BUFFER

!                            ______________________
!___________________________/ WANNIER_DESCRIPTOR_T \____________________________
!
!> @brief Descriptor for a set of Wannier orbitals.
!>
!> Contains information about the number of Wannier orbitals as well as their
!> spatial properties in the primitive and supercell.
!>
!> @var NUM_PRIM_ORBITALS  Number of orbitals in one primitive cell
!> @var PRIM_CELL          Pointer to the primitive-cell structure
!> @var CENTERS            Wannier centers in direct coordinates w.r.t. the
!>                         supercell
!> @var MIN_IMAGES         Holds minimum-image vectors w.r.t. Wannier functions
!>                         in the periodic supercell.
!> @var MIN_IMAGES_TO_ATOM Holds minimum-image vectors w.r.t. Wannier functions
!>                         in the periodic supercell and an atomic coordinate
!>                         in the primitive cell.
!_______________________________________________________________________________
   TYPE WANNIER_DESCRIPTOR_T
      INTEGER                          :: NUM_PRIM_ORBITALS
      TYPE(PRIM_CELL_T), POINTER       :: PRIM_CELL
      REAL(q), ALLOCATABLE             :: CENTERS(:, :)
      TYPE(MIN_IMAGE_T), ALLOCATABLE   :: MIN_IMAGES(:, :, :)
      TYPE(MIN_IMAGE_T), ALLOCATABLE   :: MIN_IMAGES_TO_ATOM(:, :, :)

      CONTAINS

#ifdef VASP_HDF5
      PROCEDURE,  PUBLIC   :: WRITE_HDF5 => WANNIER_DESCRIPTOR_WRITE_HDF5
      PROCEDURE,  PUBLIC   :: READ_HDF5 => WANNIER_DESCRIPTOR_READ_HDF5
#endif
      PROCEDURE,  PUBLIC   :: CALC_MIN_IMAGES => WANNIER_DESCRIPTOR_CALC_MIN_IMAGES
      PROCEDURE,  PUBLIC   :: CALC_MIN_IMAGES_TO_ATOM => WANNIER_DESCRIPTOR_CALC_MIN_IMAGES_TO_ATOM
      PROCEDURE,  PUBLIC   :: ENFORCE_ATOM_CENTERED => WANNIER_DESCRIPTOR_ENFORCE_ATOM_CENTERED
   END TYPE WANNIER_DESCRIPTOR_T

   INTERFACE WANNIER_DESCRIPTOR_T
      PROCEDURE   :: WANNIER_DESCRIPTOR_INIT
   END INTERFACE

!                             ____________________
!____________________________/ WANNIER_MAT_SYMM_T \_____________________________
!
!> @brief Represents a Wannier-space matrix in the periodic supercell.
!>
!> \f[ \braket{\mathcal{W}_{al} | \hat{A} | \mathcal{W}_{b0}} \f]
!>
!> The operator \f$\hat{A}\f$ needs to be lattice periodic since only a single
!> strip of the entire matrix is stored.
!>
!> @var MAT_SIZE     Size of the matrix along one of its dimensions. This is
!>                   equivalent to the number of Wannier states in the
!>                   supercell.
!> @var DESCRIPTOR   Pointer to Wannier descriptor
!> @var ELEMENTS     Matrix elements arranged in a memory layout optimized for
!>                   Fourier sums over cells: ELEMENTS(a, b, l)
!_______________________________________________________________________________
   TYPE WANNIER_MAT_SYMM_T
      INTEGER,                               PUBLIC   :: MAT_SIZE
      TYPE(WANNIER_DESCRIPTOR_T), POINTER,   PUBLIC   :: DESCRIPTOR
      GDEF, ALLOCATABLE,                     PUBLIC   :: ELEMENTS(:, :, :)

      CONTAINS

      ! Shorthand and generic interfaces to procedures.
      PROCEDURE,  PUBLIC   :: UNPACK => WANNIER_MAT_SYMM_UNPACK
      PROCEDURE,  PUBLIC   :: FOURIER_SUM => WANNIER_MAT_SYMM_FOURIER_SUM
      PROCEDURE,  PUBLIC   :: DIAG_AT_K => WANNIER_MAT_SYMM_DIAG_AT_K
#ifdef VASP_HDF5
      PROCEDURE,  PUBLIC   :: WRITE_HDF5 => WANNIER_MAT_SYMM_WRITE_HDF5
      PROCEDURE,  PUBLIC   :: READ_HDF5 => WANNIER_MAT_SYMM_READ_HDF5
#endif
   END TYPE

!                             _____________________
!____________________________/ WANNIER_MAT_MODES_T \____________________________
!
!> @brief Represents a perturbed Wannier-space matrix in the periodic supercell.
!>
!> \f[ \braket{\mathcal{W}_{al} | \partial_\nu \hat{A} | \mathcal{W}_{bl'}} \f]
!>
!> The perturbation is confined to one primitive cell with the number of modes
!> equal to the number of independent Cartesian displacements for each
!> primitive-cell atom (\nu = 3 N_p). Since the perturbation is localized, no
!> translational invariance is assumed and the matrix is stored in its entirety
!> w.r.t. the Wannier states (MAT_SIZE x MAT_SIZE).
!>
!> @var NUM_MODES    Number of modes/degrees of freedom
!> @var MAT_SIZE     Size of the matrix along one of its dimensions. This is
!>                   equivalent to the number of Wannier states in the
!>                   supercell.
!> @var DESCRIPTOR   Pointer to Wannier descriptor
!> @var ELEMENTS     Matrix elements arranged in a memory layout optimized for
!>                   Fourier sums over cells: ELEMENTS(a, b, l, l', \nu)
!_______________________________________________________________________________
   TYPE WANNIER_MAT_MODES_T
      INTEGER,                               PUBLIC   :: NUM_MODES
      INTEGER,                               PUBLIC   :: MAT_SIZE
      TYPE(WANNIER_DESCRIPTOR_T), POINTER,   PUBLIC   :: DESCRIPTOR
      GDEF, ALLOCATABLE,                     PUBLIC   :: ELEMENTS(:, :, :, :, :)

      CONTAINS

      ! Procedure names used for overloading. These are never used directly.
      PROCEDURE,  PRIVATE  :: UNPACK_1 => WANNIER_MAT_MODES_UNPACK_SINGLE_MODE
      PROCEDURE,  PRIVATE  :: UNPACK_2 => WANNIER_MAT_MODES_UNPACK_MODES

      ! Shorthand and generic interfaces to procedures.
      PROCEDURE,  PUBLIC   :: SET_MODE => WANNIER_MAT_MODES_SET_MODE
      GENERIC,    PUBLIC   :: UNPACK => UNPACK_1, UNPACK_2
      PROCEDURE,  PUBLIC   :: FOURIER_SUM => WANNIER_MAT_MODES_FOURIER_SUM
      PROCEDURE,  PUBLIC   :: TO_BLOCH => WANNIER_MAT_MODES_TO_BLOCH
#ifdef VASP_HDF5
      PROCEDURE,  PUBLIC   :: WRITE_HDF5 => WANNIER_MAT_MODES_WRITE_HDF5
      PROCEDURE,  PUBLIC   :: READ_HDF5 => WANNIER_MAT_MODES_READ_HDF5
#endif
   END TYPE

   ! WANNIER_MAT_SYMM_T Constructor
   INTERFACE WANNIER_MAT_SYMM_T
      PROCEDURE   :: WANNIER_MAT_SYMM_INIT
      PROCEDURE   :: WANNIER_MAT_SYMM_PACK
   END INTERFACE

   ! WANNIER_MAT_MODES_T Constructor
   INTERFACE WANNIER_MAT_MODES_T
      PROCEDURE   :: WANNIER_MAT_MODES_INIT
      PROCEDURE   :: WANNIER_MAT_MODES_PACK
   END INTERFACE

   
   CONTAINS

   SUBROUTINE CALC_WANNIER_SHIFT_MATRIX( &
      W_WAN, GRID, GRIDC, GRIDUS, C_TO_US, &
      LATT_CUR, T_INFO, P, LMDIM, LOVERL, IRDMAX, &
      G_INDICES, SPIN, SHIFT_MAT)

      USE wave_struct_def,    ONLY: wavespin
      USE mgrid_struct_def,   ONLY: grid_3d, transit
      USE poscar_struct_def,  ONLY: latt, type_info
      USE pseudo_struct_def,  ONLY: potcar
      USE WANNIER,            ONLY: WANNMAT
      IMPLICIT NONE

      TYPE(wavespin),            INTENT(IN)  :: W_WAN
      TYPE(grid_3d),             INTENT(IN)  :: GRID, GRIDC, GRIDUS
      TYPE(transit),             INTENT(IN)  :: C_TO_US
      TYPE(latt),                INTENT(IN)  :: LATT_CUR
      TYPE(type_info),           INTENT(IN)  :: T_INFO
      TYPE(potcar),              INTENT(IN)  :: P(T_INFO%NTYP)
      INTEGER,                   INTENT(IN)  :: LMDIM
      LOGICAL,                   INTENT(IN)  :: LOVERL
      INTEGER,                   INTENT(IN)  :: IRDMAX
      INTEGER,                   INTENT(IN)  :: G_INDICES(3)
      INTEGER,                   INTENT(IN)  :: SPIN
      COMPLEX(q), ALLOCATABLE,   INTENT(OUT) :: SHIFT_MAT(:, :)

      INTEGER              :: NUM_WAN
      REAL(q), ALLOCATABLE :: SHIFT_MAT_GAMMA(:, :)

      NUM_WAN = W_WAN%WDES%NB_TOT

      ALLOCATE(SHIFT_MAT(NUM_WAN, NUM_WAN))

#ifdef gammareal
      ! Used as storage for real/imaginary part in gamma-only version.
      ALLOCATE(SHIFT_MAT_GAMMA(NUM_WAN, NUM_WAN))

      ! Calculate < W_m | cos(G . r) | W_n >.
      CALL WANNMAT(  W_WAN, W_WAN%WDES, GRID, GRIDC, GRIDUS, C_TO_US, LATT_CUR, T_INFO, &
                     P, LMDIM, LOVERL, IRDMAX, G_INDICES, 1, SPIN, SHIFT_MAT_GAMMA, .FALSE.)

      ! First, store real part.
      SHIFT_MAT = SHIFT_MAT_GAMMA

      ! Calculate < W_m | sin(G . r) | W_n >.
      CALL WANNMAT(  W_WAN, W_WAN%WDES, GRID, GRIDC, GRIDUS, C_TO_US, LATT_CUR, T_INFO, &
                     P, LMDIM, LOVERL, IRDMAX, G_INDICES, 1, SPIN, SHIFT_MAT_GAMMA, .TRUE.)

      ! Then, add imaginary part.
      SHIFT_MAT = SHIFT_MAT + CMPLX(0._q, SHIFT_MAT_GAMMA, q)
#else
      ! Calculate < W_m | exp(i G . r) | W_n >.
      CALL WANNMAT(  W_WAN, W_WAN%WDES, GRID, GRIDC, GRIDUS, C_TO_US, LATT_CUR, T_INFO, &
                     P, LMDIM, LOVERL, IRDMAX, G_INDICES, 1, SPIN, SHIFT_MAT, .FALSE.)
#endif
   END SUBROUTINE CALC_WANNIER_SHIFT_MATRIX


   ! _test_mani direct coordinates
   SUBROUTINE CALC_WANNIER_CENTERS( &
      W_WAN, GRID, GRIDC, GRIDUS, C_TO_US, &
      LATT_CUR, T_INFO, P, LMDIM, LOVERL, IRDMAX, &
      WAN_CENTERS)

      USE wave_struct_def,    ONLY: wavespin, wavedes1
      USE mgrid_struct_def,   ONLY: grid_3d, transit
      USE poscar_struct_def,  ONLY: latt, type_info
      USE pseudo_struct_def,  ONLY: potcar
      USE WAVE,               ONLY: SETWDES
      USE CONSTANT,           ONLY: TPI
      IMPLICIT NONE

      TYPE(wavespin),   INTENT(IN)  :: W_WAN
      TYPE(grid_3d),    INTENT(IN)  :: GRID, GRIDC, GRIDUS
      TYPE(transit),    INTENT(IN)  :: C_TO_US
      TYPE(latt),       INTENT(IN)  :: LATT_CUR
      TYPE(type_info),  INTENT(IN)  :: T_INFO
      TYPE(potcar),     INTENT(IN)  :: P(T_INFO%NTYP)
      INTEGER,          INTENT(IN)  :: LMDIM
      LOGICAL,          INTENT(IN)  :: LOVERL
      INTEGER,          INTENT(IN)  :: IRDMAX
      REAL(q),          INTENT(OUT) :: WAN_CENTERS(:, :)

      INTEGER                 :: I
      INTEGER                 :: NUM_WAN, WAN_IDX
      TYPE(wavedes1)          :: WDES_K
      COMPLEX(q), ALLOCATABLE :: SHIFT_MAT(:, :)
      INTEGER                 :: G_INDICES(3)

      NUM_WAN = SIZE(WAN_CENTERS, 2)
      WAN_CENTERS = 0

      CALL SETWDES(W_WAN%WDES, WDES_K, 1)

      DO I = 1, 3
         G_INDICES = 0
         G_INDICES(I) = 1

         CALL CALC_WANNIER_SHIFT_MATRIX( &
            W_WAN, GRID, GRIDC, GRIDUS, C_TO_US, LATT_CUR, T_INFO, P, &
            LMDIM, LOVERL, IRDMAX, G_INDICES, 1, SHIFT_MAT)

         DO WAN_IDX = 1, NUM_WAN
            ! Yes, index order is swapped; avoids unnecessary allocations.
            WAN_CENTERS(I, WAN_IDX) = AIMAG(LOG(SHIFT_MAT(WAN_IDX, WAN_IDX))) / TPI
         ENDDO
      ENDDO
   END SUBROUTINE CALC_WANNIER_CENTERS


   SUBROUTINE WANNIER_DESCRIPTOR_CALC_MIN_IMAGES(THIS)
      IMPLICIT NONE

      CLASS(WANNIER_DESCRIPTOR_T),  INTENT(INOUT)  :: THIS

      REAL(q)  :: PRIM_DIST(3)
      INTEGER  :: ORB_IDX_1, ORB_IDX_2
      INTEGER  :: CELL_IDX_1

      IF (ALLOCATED(THIS%MIN_IMAGES)) DEALLOCATE(THIS%MIN_IMAGES)
      ALLOCATE(THIS%MIN_IMAGES(THIS%NUM_PRIM_ORBITALS, THIS%NUM_PRIM_ORBITALS, THIS%PRIM_CELL%NUM_CELLS))

      DO CELL_IDX_1 = 1, THIS%PRIM_CELL%NUM_CELLS
         DO ORB_IDX_2 = 1, THIS%NUM_PRIM_ORBITALS
            DO ORB_IDX_1 = 1, THIS%NUM_PRIM_ORBITALS
               PRIM_DIST = THIS%CENTERS(:, ORB_IDX_1) - THIS%CENTERS(:, ORB_IDX_2)
               THIS%MIN_IMAGES(ORB_IDX_1, ORB_IDX_2, CELL_IDX_1) = MIN_IMAGE_T(PRIM_DIST, CELL_IDX_1, THIS%PRIM_CELL)
            ENDDO
         ENDDO
      ENDDO
   END SUBROUTINE WANNIER_DESCRIPTOR_CALC_MIN_IMAGES


   SUBROUTINE WANNIER_DESCRIPTOR_CALC_MIN_IMAGES_TO_ATOM(THIS)
      IMPLICIT NONE

      CLASS(WANNIER_DESCRIPTOR_T),  INTENT(INOUT)  :: THIS

      REAL(q)  :: PRIM_DIST(3)
      INTEGER  :: ORB_IDX
      INTEGER  :: PRIM_IDX
      INTEGER  :: CELL_IDX

      IF (ALLOCATED(THIS%MIN_IMAGES_TO_ATOM)) DEALLOCATE(THIS%MIN_IMAGES_TO_ATOM)
      ALLOCATE(THIS%MIN_IMAGES_TO_ATOM(THIS%NUM_PRIM_ORBITALS, THIS%PRIM_CELL%NUM_CELLS, THIS%PRIM_CELL%NUM_PRIM_ATOMS))

      DO PRIM_IDX = 1, THIS%PRIM_CELL%NUM_PRIM_ATOMS
         DO CELL_IDX = 1, THIS%PRIM_CELL%NUM_CELLS
            DO ORB_IDX = 1, THIS%NUM_PRIM_ORBITALS
               PRIM_DIST = THIS%CENTERS(:, ORB_IDX) - THIS%PRIM_CELL%GET_POSION(1, PRIM_IDX)
               THIS%MIN_IMAGES_TO_ATOM(ORB_IDX, CELL_IDX, PRIM_IDX) = MIN_IMAGE_T(PRIM_DIST, CELL_IDX, THIS%PRIM_CELL)!, -1._q)
            ENDDO
         ENDDO
      ENDDO
   END SUBROUTINE WANNIER_DESCRIPTOR_CALC_MIN_IMAGES_TO_ATOM


   TYPE(WANNIER_DESCRIPTOR_T) FUNCTION WANNIER_DESCRIPTOR_INIT( &
      W_WANNIER, NUM_PRIM_ORBITALS, PRIM_CELL, GRID, GRIDC, GRIDUS, &
      C_TO_US, T_INFO, P, LMDIM, LOVERL, IRDMAX, FORCE_ATOM_CENTER) RESULT(DESCRIPTOR)
      
      USE wave_struct_def,    ONLY: wavespin
      USE mgrid_struct_def,   ONLY: grid_3d, transit
      USE poscar_struct_def,  ONLY: type_info
      USE pseudo_struct_def,  ONLY: potcar
      IMPLICIT NONE

      TYPE(wavespin),               INTENT(IN)     :: W_WANNIER
      INTEGER,                      INTENT(IN)     :: NUM_PRIM_ORBITALS
      TYPE(PRIM_CELL_T), TARGET,    INTENT(INOUT)  :: PRIM_CELL
      TYPE(grid_3d),                INTENT(IN)     :: GRID, GRIDC, GRIDUS
      TYPE(transit),                INTENT(IN)     :: C_TO_US
      TYPE(type_info),              INTENT(IN)     :: T_INFO
      TYPE(potcar),                 INTENT(IN)     :: P(T_INFO%NTYP)
      INTEGER,                      INTENT(IN)     :: LMDIM
      LOGICAL,                      INTENT(IN)     :: LOVERL
      INTEGER,                      INTENT(IN)     :: IRDMAX
      LOGICAL, OPTIONAL,            INTENT(IN)     :: FORCE_ATOM_CENTER

      LOGICAL  :: MY_FORCE

      CALL PRIM_CELL%SETUP_ADDITIONAL_DATA()

      MY_FORCE = .FALSE.
      IF (PRESENT(FORCE_ATOM_CENTER)) MY_FORCE = FORCE_ATOM_CENTER

      DESCRIPTOR%NUM_PRIM_ORBITALS = NUM_PRIM_ORBITALS
      DESCRIPTOR%PRIM_CELL => PRIM_CELL

      IF (ALLOCATED(DESCRIPTOR%CENTERS)) DEALLOCATE(DESCRIPTOR%CENTERS)
      ALLOCATE(DESCRIPTOR%CENTERS(3, DESCRIPTOR%NUM_PRIM_ORBITALS))

      CALL CALC_WANNIER_CENTERS( &
         W_WANNIER, GRID, GRIDC, GRIDUS, C_TO_US, PRIM_CELL%SUPER_LATT, &
         T_INFO, P, LMDIM, LOVERL, IRDMAX, DESCRIPTOR%CENTERS)

      IF (MY_FORCE) THEN
         CALL DESCRIPTOR%ENFORCE_ATOM_CENTERED()
      ELSE
         CALL DESCRIPTOR%CALC_MIN_IMAGES()
         CALL DESCRIPTOR%CALC_MIN_IMAGES_TO_ATOM()
      ENDIF
   END FUNCTION WANNIER_DESCRIPTOR_INIT


   SUBROUTINE WANNIER_DESCRIPTOR_ENFORCE_ATOM_CENTERED(THIS)
      IMPLICIT NONE

      CLASS(WANNIER_DESCRIPTOR_T),  INTENT(INOUT)  :: THIS

      INTEGER  :: ORB_IDX
      INTEGER  :: P_IDX
      REAL(q)  :: SHIFT(3), CLOSEST_SHIFT(3)
      REAL(q)  :: DIST, CLOSEST_DIST
      INTEGER  :: X, Y, Z

      DO ORB_IDX = 1, THIS%NUM_PRIM_ORBITALS
         CLOSEST_DIST = HUGE(1._q)
         DO P_IDX = 1, THIS%PRIM_CELL%NUM_PRIM_ATOMS
            DO Z = -1, 1
               DO Y = -1, 1
                  DO X = -1, 1
                     SHIFT = MATMUL( &
                        THIS%PRIM_CELL%SUPER_LATT%A, THIS%PRIM_CELL%GET_POSION(1, P_IDX) - &
                        THIS%CENTERS(:, ORB_IDX)) + &
                        MATMUL(THIS%PRIM_CELL%PRIM_LATT%A, [REAL(q):: X, Y, Z])
                     DIST = NORM2(SHIFT)
                     IF (DIST < CLOSEST_DIST) THEN
                        CLOSEST_DIST = DIST
                        CLOSEST_SHIFT = SHIFT
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         THIS%CENTERS(:, ORB_IDX) = THIS%CENTERS(:, ORB_IDX) + MATMUL(CLOSEST_SHIFT, THIS%PRIM_CELL%SUPER_LATT%B)
      ENDDO

      CALL THIS%CALC_MIN_IMAGES()
      CALL THIS%CALC_MIN_IMAGES_TO_ATOM()
   END SUBROUTINE


#ifdef VASP_HDF5
   SUBROUTINE WANNIER_DESCRIPTOR_WRITE_HDF5(THIS, LOC_ID, GROUP_PATH)
      USE hdf5,         ONLY: HID_T
      USE vhdf5_base,   ONLY: vh5_write, vh5_error
      IMPLICIT NONE
   
      CLASS(WANNIER_DESCRIPTOR_T),  INTENT(IN)  :: THIS
      INTEGER(HID_T),               INTENT(IN)  :: LOC_ID
      CHARACTER(*),                 INTENT(IN)  :: GROUP_PATH

      VH5_CHECK(vh5_write(LOC_ID, GROUP_PATH // "/num_prim_orbitals", THIS%NUM_PRIM_ORBITALS))
      VH5_CHECK(vh5_write(LOC_ID, GROUP_PATH // "/centers", THIS%CENTERS))
   END SUBROUTINE

   SUBROUTINE WANNIER_DESCRIPTOR_READ_HDF5(THIS, PRIM_CELL, LOC_ID, GROUP_PATH)
      USE hdf5,         ONLY: HID_T
      USE vhdf5_base,   ONLY: vh5_read, vh5_error
      IMPLICIT NONE
   
      CLASS(WANNIER_DESCRIPTOR_T),  INTENT(INOUT)  :: THIS
      TYPE(PRIM_CELL_T), TARGET,    INTENT(IN)     :: PRIM_CELL
      INTEGER(HID_T),               INTENT(IN)     :: LOC_ID
      CHARACTER(*),                 INTENT(IN)     :: GROUP_PATH
   
      THIS%PRIM_CELL => PRIM_CELL

      CALL THIS%PRIM_CELL%SETUP_ADDITIONAL_DATA()
   
      VH5_CHECK(vh5_read(LOC_ID, GROUP_PATH // "/num_prim_orbitals", THIS%NUM_PRIM_ORBITALS))
      IF (ALLOCATED(THIS%CENTERS)) DEALLOCATE(THIS%CENTERS)
      ALLOCATE(THIS%CENTERS(3, THIS%NUM_PRIM_ORBITALS))
      VH5_CHECK(vh5_read(LOC_ID, GROUP_PATH // "/centers", THIS%CENTERS))

      CALL THIS%CALC_MIN_IMAGES()
      CALL THIS%CALC_MIN_IMAGES_TO_ATOM()
   END SUBROUTINE
#endif


!                            _______________________
!___________________________/ WANNIER_MAT_SYMM_INIT \___________________________
!
!> @brief Basic constructor of the WANNIER_MAT_SYMM_T type.
!>
!> The matrix is setup with the correct size but is not initialized.
!>
!> @param[in] DESCRIPTOR   Wannier descriptor that holds primitive and super-
!>                         cell spatial information
!_______________________________________________________________________________
   TYPE(WANNIER_MAT_SYMM_T) FUNCTION WANNIER_MAT_SYMM_INIT(DESCRIPTOR) RESULT(WANNIER_MAT)
      IMPLICIT NONE

      TYPE(WANNIER_DESCRIPTOR_T), TARGET, INTENT(IN)  :: DESCRIPTOR

      IF (ALLOCATED(WANNIER_MAT%ELEMENTS)) DEALLOCATE(WANNIER_MAT%ELEMENTS)

      WANNIER_MAT%DESCRIPTOR => DESCRIPTOR
      WANNIER_MAT%MAT_SIZE = DESCRIPTOR%PRIM_CELL%NUM_CELLS * DESCRIPTOR%NUM_PRIM_ORBITALS

      ASSOCIATE(  N_ORBS => DESCRIPTOR%NUM_PRIM_ORBITALS, &
                  N_CELLS => DESCRIPTOR%PRIM_CELL%NUM_CELLS)
         ALLOCATE(WANNIER_MAT%ELEMENTS( N_ORBS, N_ORBS, N_CELLS))
      END ASSOCIATE
   END FUNCTION

!                            _______________________
!___________________________/ WANNIER_MAT_SYMM_PACK \___________________________
!
!> @brief Constructor for WANNIER_MAT_SYMM_T type with a cell-periodic Wannier-
!> space matrix.
!>
!> @param[in] DESCRIPTOR   Wannier descriptor that holds primitive and super-
!>                         cell spatial information
!> @param[in] MATRIX_STRIP One strip of Wannier-space matrix. Width of the strip
!>                         must equal number of Wannier orbitals in one
!>                         primitive cell
!_______________________________________________________________________________
   TYPE(WANNIER_MAT_SYMM_T) FUNCTION WANNIER_MAT_SYMM_PACK(DESCRIPTOR, MATRIX_STRIP) RESULT(WANNIER_MAT)
      IMPLICIT NONE

      TYPE(WANNIER_DESCRIPTOR_T),   INTENT(IN)  :: DESCRIPTOR
      GDEF,                         INTENT(IN)  :: MATRIX_STRIP(:, :)

      INTEGER  :: MAT_SHAPE(2)
      INTEGER  :: CELL_IDX

      WANNIER_MAT = WANNIER_MAT_SYMM_T(DESCRIPTOR)

      MAT_SHAPE = SHAPE(MATRIX_STRIP)

      IF (MAT_SHAPE(1) /= WANNIER_MAT%MAT_SIZE) CALL vtutor%bug( &
         "Dimension 1 of input matrix (" // str(MAT_SHAPE(1)) // ") differs from information " // &
         "found in descriptor (" // str(WANNIER_MAT%MAT_SIZE) // ").", __FILE__, __LINE__)

      IF (MAT_SHAPE(2) /= DESCRIPTOR%NUM_PRIM_ORBITALS) CALL vtutor%bug( &
         "Dimension 2 of input matrix (" // str(MAT_SHAPE(2)) // ") differs from information " // &
         "found in descriptor (" // str(DESCRIPTOR%NUM_PRIM_ORBITALS) // ").", __FILE__, __LINE__)

      ASSOCIATE(N_ORBS => WANNIER_MAT%DESCRIPTOR%NUM_PRIM_ORBITALS)
         DO CELL_IDX = 1, WANNIER_MAT%DESCRIPTOR%PRIM_CELL%NUM_CELLS
            WANNIER_MAT%ELEMENTS(:, :, CELL_IDX) = &
            MATRIX_STRIP((CELL_IDX - 1) * N_ORBS + 1 : CELL_IDX * N_ORBS, :)
         ENDDO
      END ASSOCIATE
   END FUNCTION

!                           _________________________
!__________________________/ WANNIER_MAT_SYMM_UNPACK \__________________________
!
!> @brief Unpacks Wannier matrix data to Array.
!>
!> Unpacks data stored in the optimized memory layout back into a square
!> Wannier-space matrix (MAT_SIZE x MAT_SIZE). Alternatively, only a single
!> strip can be unpacked by passing the correct shape (MAT_SIZE x
!> NUM_PRIM_ORBITALS). All other matrix shapes are illegal.
!>
!> @param[in] THIS         Instance reference
!> @param[out] MATRIX_OUT  Square matrix or matrix strip
!_______________________________________________________________________________
   SUBROUTINE WANNIER_MAT_SYMM_UNPACK(THIS, MATRIX_OUT)
      IMPLICIT NONE

      CLASS(WANNIER_MAT_SYMM_T), INTENT(IN)  :: THIS
      GDEF,                      INTENT(OUT) :: MATRIX_OUT(:, :)

      INTEGER  :: MAT_SHAPE(2)
      INTEGER  :: CELL_IDX_1, CELL_IDX_2

      MAT_SHAPE = SHAPE(MATRIX_OUT)

      IF (MAT_SHAPE(1) /= THIS%MAT_SIZE) CALL vtutor%bug( &
         "Dimension 1 of output matrix wrong: " // str(MAT_SHAPE(1)) // " /= " // str(THIS%MAT_SIZE), &
         __FILE__, __LINE__)

      ASSOCIATE(  N_ORBS => THIS%DESCRIPTOR%NUM_PRIM_ORBITALS, &
                  N_CELLS => THIS%DESCRIPTOR%PRIM_CELL%NUM_CELLS, &
                  EQU_CELLS => THIS%DESCRIPTOR%PRIM_CELL%EQUIV_CELLS)

         IF (MAT_SHAPE(2) == THIS%MAT_SIZE) THEN
            ! Return full MAT_SIZE x MAT_SIZE matrix
            DO CELL_IDX_2 = 1, N_CELLS
               DO CELL_IDX_1 = 1, N_CELLS
                  MATRIX_OUT( (CELL_IDX_1 - 1) * N_ORBS + 1 : CELL_IDX_1 * N_ORBS, &
                              (CELL_IDX_2 - 1) * N_ORBS + 1 : CELL_IDX_2 * N_ORBS) = &
                  THIS%ELEMENTS(:, :, EQU_CELLS(CELL_IDX_1, CELL_IDX_2))
               ENDDO
            ENDDO
         ELSE IF (MAT_SHAPE(2) == THIS%DESCRIPTOR%NUM_PRIM_ORBITALS) THEN
            ! Return MAT_SIZE x N_ORBS strip of matrix
            DO CELL_IDX_1 = 1, N_CELLS
               MATRIX_OUT((CELL_IDX_1 - 1) * N_ORBS + 1 : CELL_IDX_1 * N_ORBS, :) = &
               THIS%ELEMENTS(:, :, CELL_IDX_1)
            ENDDO
         ELSE
            CALL vtutor%bug("Dimension 2 of output matrix wrong: " // str(MAT_SHAPE(2)) // &
            " /= " // str(THIS%MAT_SIZE) // " and /= " // str(N_ORBS), __FILE__, __LINE__)
         ENDIF
      END ASSOCIATE
   END SUBROUTINE

!                        ______________________________
!_______________________/ WANNIER_MAT_SYMM_FOURIER_SUM \________________________
!
!> @brief Performs a Fourier sum over all cells.
!>
!> Transforms into the primitive cell at the given reciprocal-space vector.
!>
!> @param[in] THIS      Instance reference
!> @param[out] PRIM_MAT Primitive-cell matrix (bands x bands)
!> @param[in] K_VEC     Arbitrary k-vector in direct coordinates of the
!>                      Brillouin zone of the primitive cell
!_______________________________________________________________________________
   SUBROUTINE WANNIER_MAT_SYMM_FOURIER_SUM(THIS, PRIM_MAT, K_VEC)
      IMPLICIT NONE

      CLASS(WANNIER_MAT_SYMM_T), INTENT(IN)  :: THIS
      COMPLEX(q), ALLOCATABLE,   INTENT(OUT) :: PRIM_MAT(:, :)
      REAL(q),                   INTENT(IN)  :: K_VEC(3)

      INTEGER  :: CELL_IDX_1
      INTEGER  :: ORB_IDX_1, ORB_IDX_2
      INTEGER  :: NUM_ORBS
      REAL(q)  :: K_CART(3)

      NUM_ORBS = THIS%DESCRIPTOR%NUM_PRIM_ORBITALS

      ALLOCATE(PRIM_MAT(NUM_ORBS, NUM_ORBS))
      PRIM_MAT = 0

      ! Transform to Cartesian coordinates.
      ! _TEST_MANI DIR COORD
      K_CART = MATMUL(THIS%DESCRIPTOR%PRIM_CELL%PRIM_LATT%B, K_VEC)

      ! Loop over all primitive cells.
      DO CELL_IDX_1 = 1, THIS%DESCRIPTOR%PRIM_CELL%NUM_CELLS
         ! Double loop over all orbitals in each cell.
         DO ORB_IDX_2 = 1, NUM_ORBS
            DO ORB_IDX_1 = 1, NUM_ORBS
               PRIM_MAT(ORB_IDX_1, ORB_IDX_2) = &
               PRIM_MAT(ORB_IDX_1, ORB_IDX_2) + &
               THIS%ELEMENTS(ORB_IDX_1, ORB_IDX_2, CELL_IDX_1) * &
               THIS%DESCRIPTOR%MIN_IMAGES(ORB_IDX_1, ORB_IDX_2, CELL_IDX_1)%CALC_CELL_PHASE(-K_CART)
            ENDDO
         ENDDO
      ENDDO
   END SUBROUTINE

!                         ____________________________
!________________________/ WANNIER_MAT_SYMM_DIAG_AT_K \_________________________
!
!> @brief Performs a Fourier sum over all cells, followed by a diagonalization.
!>
!> Transforms into the primitive cell at the given reciprocal-space vector.
!> Afterwards, the matrix is diagonalized to retrieve eigenvalues and eigen-
!> vectors.
!>
!> @param[in] THIS         Instance reference
!> @param[in] K_VEC        Arbitrary k-vector in direct coordinates of the
!>                         Brillouin zone of the primitive cell
!> @param[out] EIGEN_VALS  Eigenvalues
!> @param[out] EIGEN_VECS  Eigenvectors
!_______________________________________________________________________________
   SUBROUTINE WANNIER_MAT_SYMM_DIAG_AT_K(THIS, K_VEC, EIGEN_VALS, EIGEN_VECS)
      IMPLICIT NONE

      CLASS(WANNIER_MAT_SYMM_T), INTENT(IN)  :: THIS
      REAL(q),                   INTENT(IN)  :: K_VEC(3)
      REAL(q), ALLOCATABLE,      INTENT(OUT) :: EIGEN_VALS(:)
      COMPLEX(q), ALLOCATABLE,   INTENT(OUT) :: EIGEN_VECS(:, :)

      ! Work arrays for ZHEEV.
      INTEGER                 :: LWORK, IFAIL
      COMPLEX(q), ALLOCATABLE :: CWORK(:)
      REAL(q), ALLOCATABLE    :: RWORK(:)

      ! Perform fourier sum to primitive cell.
      CALL THIS%FOURIER_SUM(EIGEN_VECS, K_VEC)

      ASSOCIATE(N_ORBS => THIS%DESCRIPTOR%NUM_PRIM_ORBITALS)
         ALLOCATE(EIGEN_VALS(N_ORBS))

         ! Perform the diagonalization.
         LWORK = 3 * N_ORBS
         ALLOCATE(CWORK(LWORK), RWORK(LWORK))
         CALL ZHEEV('V', 'U', N_ORBS, EIGEN_VECS, N_ORBS, EIGEN_VALS, CWORK, LWORK, RWORK, IFAIL)
      END ASSOCIATE

      IF (IFAIL /= 0) CALL vtutor%error("WANNIER_MAT_DIAG_AT_K: ZHEEV returned " // str(IFAIL))
   END SUBROUTINE

#ifdef VASP_HDF5
   SUBROUTINE WANNIER_MAT_SYMM_WRITE_HDF5(THIS, LOC_ID, GROUP_PATH)
      USE hdf5,         ONLY: HID_T
      USE vhdf5_base,   ONLY: vh5_write, vh5_error
      IMPLICIT NONE

      CLASS(WANNIER_MAT_SYMM_T), INTENT(IN)  :: THIS
      INTEGER(HID_T),            INTENT(IN)  :: LOC_ID
      CHARACTER(*),              INTENT(IN)  :: GROUP_PATH

      VH5_CHECK(vh5_write(LOC_ID, GROUP_PATH // "/matrix_size", THIS%MAT_SIZE))
      VH5_CHECK(vh5_write(LOC_ID, GROUP_PATH // "/elements", THIS%ELEMENTS))
   END SUBROUTINE

   SUBROUTINE WANNIER_MAT_SYMM_READ_HDF5(THIS, DESCRIPTOR, LOC_ID, GROUP_PATH)
      USE hdf5,         ONLY: HID_T
      USE vhdf5_base,   ONLY: vh5_read, vh5_error
      IMPLICIT NONE

      CLASS(WANNIER_MAT_SYMM_T),          INTENT(INOUT)  :: THIS
      TYPE(WANNIER_DESCRIPTOR_T), TARGET, INTENT(IN)     :: DESCRIPTOR
      INTEGER(HID_T),                     INTENT(IN)     :: LOC_ID
      CHARACTER(*),                       INTENT(IN)     :: GROUP_PATH

      THIS%DESCRIPTOR => DESCRIPTOR

      IF (ALLOCATED(THIS%ELEMENTS)) DEALLOCATE(THIS%ELEMENTS)

      VH5_CHECK(vh5_read(LOC_ID, GROUP_PATH // "/matrix_size", THIS%MAT_SIZE))

      ASSOCIATE(  N_ORBS => THIS%DESCRIPTOR%NUM_PRIM_ORBITALS, &
                  N_CELLS => THIS%DESCRIPTOR%PRIM_CELL%NUM_CELLS)
         ALLOCATE(THIS%ELEMENTS(N_ORBS, N_ORBS, N_CELLS))
      END ASSOCIATE

      VH5_CHECK(vh5_read(LOC_ID, GROUP_PATH // "/elements", THIS%ELEMENTS))
   END SUBROUTINE
#endif


!                           ________________________
!__________________________/ WANNIER_MAT_MODES_INIT \___________________________
!
!> @brief Basic constructor of the WANNIER_MAT_MODES_T type.
!>
!> The matrix is setup with the correct size but is not initialized.
!>
!> @param[in] DESCRIPTOR   Wannier descriptor that holds primitive and super-
!>                         cell spatial information
!> @param[in] NUM_MODES    Number of independent vibrational modes the matrix
!>                         holds
!_______________________________________________________________________________
   TYPE(WANNIER_MAT_MODES_T) FUNCTION WANNIER_MAT_MODES_INIT(DESCRIPTOR, NUM_MODES) RESULT(WANNIER_MAT)
      IMPLICIT NONE

      TYPE(WANNIER_DESCRIPTOR_T), TARGET, INTENT(IN)  :: DESCRIPTOR
      INTEGER,                            INTENT(IN)  :: NUM_MODES

      IF (ALLOCATED(WANNIER_MAT%ELEMENTS)) DEALLOCATE(WANNIER_MAT%ELEMENTS)

      WANNIER_MAT%DESCRIPTOR => DESCRIPTOR
      WANNIER_MAT%MAT_SIZE = DESCRIPTOR%PRIM_CELL%NUM_CELLS * DESCRIPTOR%NUM_PRIM_ORBITALS

      ! _test_mani check?
      WANNIER_MAT%NUM_MODES = NUM_MODES

      ASSOCIATE(  N_ORBS => DESCRIPTOR%NUM_PRIM_ORBITALS, &
                  N_CELLS => DESCRIPTOR%PRIM_CELL%NUM_CELLS)
         ALLOCATE(WANNIER_MAT%ELEMENTS(N_ORBS, N_ORBS, N_CELLS, N_CELLS, NUM_MODES))
      END ASSOCIATE
   END FUNCTION

!                           ________________________
!__________________________/ WANNIER_MAT_MODES_PACK \___________________________
!
!> @brief Constructor for WANNIER_MAT_MODES_T type with a perturbed Wannier-
!> space matrix.
!>
!> @param[in] DESCRIPTOR   Wannier descriptor that holds primitive and super-
!>                         cell spatial information
!> @param[in] MATRIX_IN    Wannier-space matrix. The first two dimensions label
!>                         the Wannier states (MAT_SIZE x MAT_SIZE). The last
!>                         dimension labels the vibrational modes (NUM_MODES).
!_______________________________________________________________________________
   TYPE(WANNIER_MAT_MODES_T) FUNCTION WANNIER_MAT_MODES_PACK(DESCRIPTOR, MATRIX_IN) RESULT(WANNIER_MAT)
      IMPLICIT NONE

      TYPE(WANNIER_DESCRIPTOR_T),   INTENT(IN)  :: DESCRIPTOR
      GDEF,                         INTENT(IN)  :: MATRIX_IN(:, :, :)

      INTEGER  :: MAT_SIZE
      INTEGER  :: CELL_IDX_1, CELL_IDX_2
      INTEGER  :: MODE_IDX

      MAT_SIZE = SIZE(MATRIX_IN, 1)

      IF (MAT_SIZE /= SIZE(MATRIX_IN, 2)) CALL vtutor%bug( &
         "Input matrix not square! " // str(SHAPE(MATRIX_IN)), __FILE__, __LINE__)

      WANNIER_MAT = WANNIER_MAT_MODES_T(DESCRIPTOR, SIZE(MATRIX_IN, 3))

      IF (MAT_SIZE /= WANNIER_MAT%MAT_SIZE) CALL vtutor%bug( &
         "Inconsistent matrix dimensions! " // str(MAT_SIZE) // " /= " // str(WANNIER_MAT%MAT_SIZE), &
         __FILE__, __LINE__)

      ASSOCIATE(  N_ORBS => WANNIER_MAT%DESCRIPTOR%NUM_PRIM_ORBITALS, &
                  N_CELLS => WANNIER_MAT%DESCRIPTOR%PRIM_CELL%NUM_CELLS)

         DO MODE_IDX = 1, WANNIER_MAT%NUM_MODES
            DO CELL_IDX_2 = 1, N_CELLS
               DO CELL_IDX_1 = 1, N_CELLS
                  WANNIER_MAT%ELEMENTS(:, :, CELL_IDX_1, CELL_IDX_2, MODE_IDX) = &
                  MATRIX_IN(  (CELL_IDX_1 - 1) * N_ORBS + 1 : CELL_IDX_1 * N_ORBS, &
                              (CELL_IDX_2 - 1) * N_ORBS + 1 : CELL_IDX_2 * N_ORBS, &
                              MODE_IDX)
               ENDDO
            ENDDO
         ENDDO
      END ASSOCIATE
   END FUNCTION

!                         ____________________________
!________________________/ WANNIER_MAT_MODES_SET_MODE \_________________________
!
!> @brief Allows to set a Wannier-space matrix explicitly for one individual
!> mode.
!>
!> @param[in,out] THIS  Instance reference
!> @param[in] MATRIX_IN Square Wannier-space matrix (MAT_SIZE x MAT_SIZE)
!> @param[in] MODE_IDX  Index of vibrational mode for which to set the matrix
!_______________________________________________________________________________
   SUBROUTINE WANNIER_MAT_MODES_SET_MODE(THIS, MATRIX_IN, MODE_IDX)
      IMPLICIT NONE

      CLASS(WANNIER_MAT_MODES_T),   INTENT(INOUT)  :: THIS
      GDEF,                         INTENT(IN)     :: MATRIX_IN(:, :)
      INTEGER,                      INTENT(IN)     :: MODE_IDX

      INTEGER  :: CELL_IDX_1, CELL_IDX_2

      IF (MODE_IDX < 1 .OR. MODE_IDX > THIS%NUM_MODES) CALL vtutor%bug( &
         "Mode index out of bounds! " // str(MODE_IDX), __FILE__, __LINE__)

      IF (ANY(THIS%MAT_SIZE /= SHAPE(MATRIX_IN))) CALL vtutor%bug( &
         "Inconsistent matrix dimensions! " // str(THIS%MAT_SIZE) // " /= " // &
         str(SHAPE(MATRIX_IN)), __FILE__, __LINE__)

      ASSOCIATE(  N_ORBS => THIS%DESCRIPTOR%NUM_PRIM_ORBITALS, &
                  N_CELLS => THIS%DESCRIPTOR%PRIM_CELL%NUM_CELLS)

         DO CELL_IDX_2 = 1, N_CELLS
            DO CELL_IDX_1 = 1, N_CELLS
               THIS%ELEMENTS(:, :, CELL_IDX_1, CELL_IDX_2, MODE_IDX) = &
               MATRIX_IN(  (CELL_IDX_1 - 1) * N_ORBS + 1 : CELL_IDX_1 * N_ORBS, &
                           (CELL_IDX_2 - 1) * N_ORBS + 1 : CELL_IDX_2 * N_ORBS)
            ENDDO
         ENDDO
      END ASSOCIATE
   END SUBROUTINE

!                    ______________________________________
!___________________/ WANNIER_MAT_MODES_UNPACK_SINGLE_MODE \____________________
!
!> @brief Unpacks Wannier matrix data to Array for single mode.
!>
!> Unpacks data stored in the optimized memory layout back into a square
!> Wannier-space matrix (MAT_SIZE x MAT_SIZE) for a single mode.
!>
!> @param[in] THIS         Instance reference
!> @param[out] MATRIX_OUT  Matrix to unpack into
!> @param[in] MODE_IDX     Index of vibrational mode for which to unpack
!_______________________________________________________________________________
   SUBROUTINE WANNIER_MAT_MODES_UNPACK_SINGLE_MODE(THIS, MATRIX_OUT, MODE_IDX)
      IMPLICIT NONE

      CLASS(WANNIER_MAT_MODES_T),   INTENT(IN)  :: THIS
      GDEF,                         INTENT(OUT) :: MATRIX_OUT(:, :)
      INTEGER,                      INTENT(IN)  :: MODE_IDX

      INTEGER  :: CELL_IDX_1, CELL_IDX_2
      INTEGER  :: I
      INTEGER  :: MAT_SHAPE(2)

      MAT_SHAPE = SHAPE(MATRIX_OUT)

      IF (MODE_IDX < 1 .OR. MODE_IDX > THIS%NUM_MODES) CALL vtutor%bug( &
         "Mode index out of bounds! " // str(MODE_IDX), __FILE__, __LINE__)

      DO I = 1, 2
         IF (MAT_SHAPE(I) /= THIS%MAT_SIZE) CALL vtutor%bug( &
            "Dimension " // str(I) // " of output matrix wrong: " // str(MAT_SHAPE(I)) // &
            " /= " // str(THIS%MAT_SIZE), __FILE__, __LINE__)
      ENDDO

      ASSOCIATE(  N_ORBS => THIS%DESCRIPTOR%NUM_PRIM_ORBITALS, &
                  N_CELLS => THIS%DESCRIPTOR%PRIM_CELL%NUM_CELLS)

         DO CELL_IDX_2 = 1, N_CELLS
            DO CELL_IDX_1 = 1, N_CELLS
               MATRIX_OUT( (CELL_IDX_1 - 1) * N_ORBS + 1 : CELL_IDX_1 * N_ORBS, &
                           (CELL_IDX_2 - 1) * N_ORBS + 1 : CELL_IDX_2 * N_ORBS) = &
               THIS%ELEMENTS(:, :, CELL_IDX_1, CELL_IDX_2, MODE_IDX)
            ENDDO
         ENDDO
      END ASSOCIATE
   END SUBROUTINE

!                       ________________________________
!______________________/ WANNIER_MAT_MODES_UNPACK_MODES \_______________________
!
!> @brief Unpacks Wannier matrix data to Array.
!>
!> @param[in] THIS         Instance reference
!> @param[out] MATRIX_OUT  Wannier-space matrix to unpack into. The first two
!>                         dimensions label the Wannier states (MAT_SIZE x
!>                         MAT_SIZE). Last dimension labels the vibrational
!>                         modes (NUM_MODES).
!_______________________________________________________________________________
   SUBROUTINE WANNIER_MAT_MODES_UNPACK_MODES(THIS, MATRIX_OUT)
      IMPLICIT NONE

      CLASS(WANNIER_MAT_MODES_T),   INTENT(IN)  :: THIS
      GDEF,                         INTENT(OUT) :: MATRIX_OUT(:, :, :)

      INTEGER  :: CELL_IDX_1, CELL_IDX_2
      INTEGER  :: I
      INTEGER  :: MAT_SHAPE(3)

      MAT_SHAPE = SHAPE(MATRIX_OUT)

      DO I = 1, 2
         IF (MAT_SHAPE(I) /= THIS%MAT_SIZE) CALL vtutor%bug( &
            "Dimension " // str(I) // " of output matrix wrong: " // str(MAT_SHAPE(I)) // &
            " /= " // str(THIS%MAT_SIZE), __FILE__, __LINE__)
      ENDDO

      IF (MAT_SHAPE(3) /= THIS%NUM_MODES) CALL vtutor%bug( &
         "Dimension 3 of output matrix wrong: " // str(MAT_SHAPE(3)) // &
         " /= " // str(THIS%NUM_MODES), __FILE__, __LINE__)

      ASSOCIATE(  N_ORBS => THIS%DESCRIPTOR%NUM_PRIM_ORBITALS, &
                  N_CELLS => THIS%DESCRIPTOR%PRIM_CELL%NUM_CELLS)
      
         DO CELL_IDX_2 = 1, N_CELLS
            DO CELL_IDX_1 = 1, N_CELLS
               MATRIX_OUT( (CELL_IDX_1 - 1) * N_ORBS + 1 : CELL_IDX_1 * N_ORBS, &
                           (CELL_IDX_2 - 1) * N_ORBS + 1 : CELL_IDX_2 * N_ORBS, :) = &
               THIS%ELEMENTS(:, :, CELL_IDX_1, CELL_IDX_2, :)
            ENDDO
         ENDDO
      END ASSOCIATE
   END SUBROUTINE

!                        _______________________________
!_______________________/ WANNIER_MAT_MODES_FOURIER_SUM \_______________________
!
!> @brief Performs a Fourier sum over all cells.
!>
!> Transforms into the primitive cell at the given reciprocal-lattice vectors.
!> The modes are left untouched, so are still in the atomic-displacements
!> picture.
!>
!> @param[in] THIS      Instance reference
!> @param[out] PRIM_MAT Primitive-cell matrix (bands, bands, modes)
!> @param[in] K_VEC     Reciprocal-space vector in direct coordinates
!>                      corresponding to the "Ket" state of the matrix element
!> @param[in] KP_VEC    Reciprocal-space vector in direct coordinates
!>                      corresponding to the "Bra" state of the matrix element
!>                      (default = conjg(K_VEC))
!_______________________________________________________________________________
   SUBROUTINE WANNIER_MAT_MODES_FOURIER_SUM(THIS, PRIM_MAT, K_VEC, KP_VEC)
      IMPLICIT NONE

      CLASS(WANNIER_MAT_MODES_T),   INTENT(IN)  :: THIS
      COMPLEX(q), ALLOCATABLE,      INTENT(OUT) :: PRIM_MAT(:, :, :)
      REAL(q),                      INTENT(IN)  :: K_VEC(3)
      REAL(q), OPTIONAL,            INTENT(IN)  :: KP_VEC(3)

      INTEGER                 :: NUM_ORBS
      INTEGER                 :: NUM_CELLS
      INTEGER                 :: CELL_IDX_1, CELL_IDX_2
      INTEGER                 :: ORB_IDX_1, ORB_IDX_2
      INTEGER                 :: PRIM_IDX
      INTEGER                 :: CART_IDX
      INTEGER                 :: MODE_IDX
      REAL(q)                 :: K_CART(3), KP_CART(3)
      LOGICAL                 :: USE_K, USE_KP
      COMPLEX(q), ALLOCATABLE :: PHASE_K(:, :), PHASE_KP(:, :)

      NUM_ORBS = THIS%DESCRIPTOR%NUM_PRIM_ORBITALS
      NUM_CELLS = THIS%DESCRIPTOR%PRIM_CELL%NUM_CELLS

      ALLOCATE(PRIM_MAT(NUM_ORBS, NUM_ORBS, THIS%NUM_MODES))
      PRIM_MAT = 0

      ALLOCATE(PHASE_K(NUM_ORBS, NUM_CELLS), PHASE_KP(NUM_ORBS, NUM_CELLS))

      ! Check which phase factors are needed.
      USE_K = ANY(ABS(K_VEC) > 1E-6_q)
      IF (PRESENT(KP_VEC)) THEN
         USE_KP = ANY(ABS(KP_VEC - K_VEC) > 1E-6_q)
      ELSE
         USE_KP = .FALSE.
      ENDIF

      ! Convert to Cartesian coordinates.
      IF (USE_K) K_CART = MATMUL(THIS%DESCRIPTOR%PRIM_CELL%PRIM_LATT%B, K_VEC)
      IF (USE_KP) KP_CART = MATMUL(THIS%DESCRIPTOR%PRIM_CELL%PRIM_LATT%B, KP_VEC)

      ! Loop over all primitive indices.
      DO PRIM_IDX = 1, THIS%DESCRIPTOR%PRIM_CELL%NUM_PRIM_ATOMS

         ! Pre-compute phase factor for k-point.
         IF (USE_K) THEN
            CALL SETUP_PHASE_BUFFER(PHASE_K, THIS%DESCRIPTOR, PRIM_IDX, K_CART)
         ELSE
            PHASE_K = 1
         ENDIF

         ! Pre-compute phase factor for kp-point.
         IF (USE_KP) THEN
            CALL SETUP_PHASE_BUFFER(PHASE_KP, THIS%DESCRIPTOR, PRIM_IDX, -KP_CART)
         ELSE
            PHASE_KP = CONJG(PHASE_K)
         ENDIF

         ! Loop over all Cartesian directions or just the requested one.
         DO CART_IDX = 1, 3
            ! Mode index of super matrix.
            MODE_IDX = PRIM_IDX * 3 - 3 + CART_IDX

            ! Double loop over all cells.
            DO CELL_IDX_2 = 1, NUM_CELLS
               DO CELL_IDX_1 = 1, NUM_CELLS

                  ! Double loop over all orbitals inside each cells.
                  DO ORB_IDX_2 = 1, NUM_ORBS
                     DO ORB_IDX_1 = 1, NUM_ORBS

                        PRIM_MAT(ORB_IDX_1, ORB_IDX_2, MODE_IDX) = &
                        PRIM_MAT(ORB_IDX_1, ORB_IDX_2, MODE_IDX) + &
                        THIS%ELEMENTS(ORB_IDX_1, ORB_IDX_2, CELL_IDX_1, CELL_IDX_2, MODE_IDX) * &
                        PHASE_KP(ORB_IDX_1, CELL_IDX_1) * PHASE_K(ORB_IDX_2, CELL_IDX_2)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   END SUBROUTINE

!                         ____________________________
!________________________/ WANNIER_MAT_MODES_TO_BLOCH \_________________________
!
!> @brief Transforms a Wannier-space matrix to the electronic Bloch picture.
!>
!> The WANNIER_MAT_MODES_T structure holding all vibrational modes is
!> transformed to a Bloch-space matrix via Wannier-Fourier interpolation. Only
!> the electronic part (band indices) is transformed, the vibrational degrees
!> of freedom remain untouched. The target k- and q-vectors can be arbitrary
!> vectors in the Brillouin zone.
!>
!> @param[in] THIS            Instance reference
!> @param[out] BLOCH_MAT      Output matrix (electronic Bloch picture)
!> @param[in] WANNIER_HAMIL   Wannier-space Hamiltonian matrix
!> @param[in] Q_VEC           Phonon wave vector in direct coordinates of the
!>                            Brillouin zone of the primitive cell
!> @param[in] K_VEC           Electron Bloch vector in direct coordinates of the
!>                            Brillouin zone of the primitive cell
!> @param[in] SPIN            Spin index
!_______________________________________________________________________________
   SUBROUTINE WANNIER_MAT_MODES_TO_BLOCH(THIS, BLOCH_MAT, WANNIER_HAMIL, Q_VEC, K_VEC)
      IMPLICIT NONE

      CLASS(WANNIER_MAT_MODES_T),   INTENT(IN)  :: THIS
      COMPLEX(q), ALLOCATABLE,      INTENT(OUT) :: BLOCH_MAT(:, :, :)
      TYPE(WANNIER_MAT_SYMM_T),     INTENT(IN)  :: WANNIER_HAMIL
      REAL(q),                      INTENT(IN)  :: Q_VEC(3)
      REAL(q),                      INTENT(IN)  :: K_VEC(3)

      INTEGER                 :: I
      COMPLEX(q), ALLOCATABLE :: UK_D(:, :), UKQ(:, :)
      REAL(q), ALLOCATABLE    :: ENERGY(:)

      ! Calculate the 'electronic' rotation matrices from Wannier to Bloch basis.
      CALL WANNIER_HAMIL%DIAG_AT_K(K_VEC, ENERGY, UK_D)
      CALL WANNIER_HAMIL%DIAG_AT_K(K_VEC + Q_VEC, ENERGY, UKQ)
      UKQ = CONJG(TRANSPOSE(UKQ))

      ! Perform Fourier sum to primitive cell.
      CALL THIS%FOURIER_SUM(BLOCH_MAT, K_VEC, K_VEC + Q_VEC)

      DO I = 1, THIS%NUM_MODES
         ! Rotate the primitive-cell matrix with the 'electronic' Wannier matrices.
         BLOCH_MAT(:, :, I) = MATMUL(MATMUL(UKQ, BLOCH_MAT(:, :, I)), UK_D)
      ENDDO
   END SUBROUTINE

!                             ____________________
!____________________________/ SETUP_PHASE_BUFFER \_____________________________
!
!> @brief Calculates the Fourier phase under a minimum-image constraint.
!>
!> Calculates \f$ e^{-i \vb k \vdot \vb R_l} \f$ under the constrant that
!> \f$ \vb R_l + \vb* \xi_a - \vb* \tau_\kappa \f$ is the minimum-image vector
!> in the periodic supercell.
!>
!> \f$k\f$ is a Bloch vector, \f$\vb R_l\f$ is the origin of unit cell \f$l\f$,
!> \f$\xi_a\f$ is the position of Wannier orbital \f$a\f$ and \f$\tau_\kappa\f$
!> is the position of primitive-cell atom \f$\kappa\f$ in unit cell 0. The
!> phase factor is calculated for each possible Wannier orbital \f$a\f$ in each
!> possible cell \f$l\f$.
!>
!> @param[in,out] PHASE_BUFFER   PHASE_BUFFER(a, l) constains the phase factor
!>                               corresponding to Wannier orbital a and cell l
!> @param[in] DESCRIPTOR         Wannier descriptor
!> @param[in] PRIM_IDX           Primitive index of atom involved in minimum-
!>                               image criterion
!> @param[in] K_CART             Bloch vector in Cartesian coordinates
!_______________________________________________________________________________
   SUBROUTINE SETUP_PHASE_BUFFER(PHASE_BUFFER, DESCRIPTOR, PRIM_IDX, K_CART)
      IMPLICIT NONE

      COMPLEX(q), ALLOCATABLE,      INTENT(OUT) :: PHASE_BUFFER(:, :)
      TYPE(WANNIER_DESCRIPTOR_T),   INTENT(IN)  :: DESCRIPTOR
      INTEGER,                      INTENT(IN)  :: PRIM_IDX
      REAL(q),                      INTENT(IN)  :: K_CART(3)

      INTEGER  :: NUM_ORBS, ORB_IDX
      INTEGER  :: NUM_CELLS, CELL_IDX

      NUM_ORBS = DESCRIPTOR%NUM_PRIM_ORBITALS
      NUM_CELLS = DESCRIPTOR%PRIM_CELL%NUM_CELLS

      ALLOCATE(PHASE_BUFFER(NUM_ORBS, NUM_CELLS))

      DO CELL_IDX = 1, NUM_CELLS
         DO ORB_IDX = 1, NUM_ORBS
            PHASE_BUFFER(ORB_IDX, CELL_IDX) = &
            DESCRIPTOR%MIN_IMAGES_TO_ATOM(ORB_IDX, CELL_IDX, PRIM_IDX)%CALC_CELL_PHASE(K_CART)
         ENDDO
      ENDDO
   END SUBROUTINE


#ifdef VASP_HDF5
   SUBROUTINE WANNIER_MAT_MODES_WRITE_HDF5(THIS, LOC_ID, GROUP_PATH)
      USE hdf5,         ONLY: HID_T
      USE vhdf5_base,   ONLY: vh5_write, vh5_error
      IMPLICIT NONE

      CLASS(WANNIER_MAT_MODES_T),   INTENT(IN)  :: THIS
      INTEGER(HID_T),               INTENT(IN)  :: LOC_ID
      CHARACTER(*),                 INTENT(IN)  :: GROUP_PATH

      VH5_CHECK(vh5_write(LOC_ID, GROUP_PATH // "/num_modes", THIS%NUM_MODES))
      VH5_CHECK(vh5_write(LOC_ID, GROUP_PATH // "/matrix_size", THIS%MAT_SIZE))
      VH5_CHECK(vh5_write(LOC_ID, GROUP_PATH // "/elements", THIS%ELEMENTS))
   END SUBROUTINE

   SUBROUTINE WANNIER_MAT_MODES_READ_HDF5(THIS, DESCRIPTOR, LOC_ID, GROUP_PATH)
      USE hdf5,         ONLY: HID_T
      USE vhdf5_base,   ONLY: vh5_read, vh5_error
      IMPLICIT NONE

      CLASS(WANNIER_MAT_MODES_T),         INTENT(INOUT)  :: THIS
      TYPE(WANNIER_DESCRIPTOR_T), TARGET, INTENT(IN)     :: DESCRIPTOR
      INTEGER(HID_T),                     INTENT(IN)     :: LOC_ID
      CHARACTER(*),                       INTENT(IN)     :: GROUP_PATH

      THIS%DESCRIPTOR => DESCRIPTOR

      IF (ALLOCATED(THIS%ELEMENTS)) DEALLOCATE(THIS%ELEMENTS)

      VH5_CHECK(vh5_read(LOC_ID, GROUP_PATH // "/num_modes", THIS%NUM_MODES))
      VH5_CHECK(vh5_read(LOC_ID, GROUP_PATH // "/matrix_size", THIS%MAT_SIZE))

      ASSOCIATE( &
         N_ORBS => THIS%DESCRIPTOR%NUM_PRIM_ORBITALS, &
         N_CELLS => THIS%DESCRIPTOR%PRIM_CELL%NUM_CELLS)

         ALLOCATE(THIS%ELEMENTS(N_ORBS, N_ORBS, N_CELLS, N_CELLS, THIS%NUM_MODES))
      END ASSOCIATE

      VH5_CHECK(vh5_read(LOC_ID, GROUP_PATH // "/elements", THIS%ELEMENTS))
   END SUBROUTINE
#endif

END MODULE wannier_mats