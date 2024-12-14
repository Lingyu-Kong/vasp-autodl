#include "symbol.inc"

MODULE elphon_util
   USE prec
   USE base,   ONLY: in_struct

   INTERFACE WRITE_MAT
      MODULE PROCEDURE WRITE_MATRIX_REAL
      MODULE PROCEDURE WRITE_MATRIX_COMPLEX
   END INTERFACE

!                                ______________
!_______________________________/ POLAR_DATA_T \________________________________
!
!> @brief Contains Born effective charges and dielectric tensor
!>
!> @var DIELECTRIC_TENSOR  Macroscopic static dielectric tensor
!> @var BORN_EFF_CHARGES   Born effective charges
!_______________________________________________________________________________
   TYPE POLAR_DATA_T
      REAL(q)              :: DIELECTRIC_TENSOR(3, 3)
      REAL(q), ALLOCATABLE :: BORN_EFF_CHARGES(:, :, :)

      CONTAINS

      PROCEDURE, PUBLIC :: CHECK => POLAR_DATA_CHECK
   END TYPE

!                                _______________
!_______________________________/ POLAR_EWALD_T \_______________________________
!
!> @brief Used to perform an Ewald sum for polar (electron-)phonon corrections.
!>
!> Mainly, provides access to the exponential Ewald cutoff for the reciprocal
!> part of the Ewald sum. For the calculation of polar corrections, the real-
!> space sum can safely be neglected.
!>
!> The G-vectors used in the sum need to be generated externally. A G-vector
!> can then be passed to the TRY_CALC routine to provide various useful
!> quantities for the polar Ewald sum.
!>
!> @var DATA            Contains Born effective charges and dielectric tensor
!> @var EXP_FACTOR      Exponential Ewald suppression factor
!> @var DIELEC_FACTOR   Factor involving dielectric tensor
!> @var BORN_FACTOR     Factor involving Born effective charges
!> @var G_CUTOFF        G-cutoff used to initialize the type
!> @var EWALD_PARAM     Ewald parameter
!> @var THRESHOLD       Ewald supression threshold
!_______________________________________________________________________________
   TYPE POLAR_EWALD_T
      TYPE(POLAR_DATA_T)   :: DATA

      REAL(q)              :: EXP_FACTOR
      REAL(q)              :: DIELEC_FACTOR
      REAL(q), ALLOCATABLE :: BORN_FACTOR(:, :)

      REAL(q)              :: G_CUTOFF
      REAL(q)              :: EWALD_PARAM
      REAL(q)              :: THRESHOLD

      CONTAINS

      PROCEDURE, PRIVATE   :: INIT_1 => POLAR_EWALD_INIT_BASE
      PROCEDURE, PRIVATE   :: INIT_2 => POLAR_EWALD_INIT_INDIV

      GENERIC, PUBLIC   :: INIT => INIT_1, INIT_2
      PROCEDURE, PUBLIC :: TRY_CALC => POLAR_EWALD_TRY_CALC
   END TYPE

   CONTAINS


!                              __________________
!_____________________________/ POLAR_DATA_CHECK \______________________________
!
!> @brief Checks if the Born-effective-charge array is correctly allocated
!>
!> @param[in] THIS      Instance reference
!> @param[in] PRIM_CELL Primitive-cell information
!_______________________________________________________________________________
   SUBROUTINE POLAR_DATA_CHECK(THIS, PRIM_CELL)
      USE msymmetry, ONLY: PRIM_CELL_T
      USE tutor,     ONLY: vtutor
      USE string,    ONLY: str
      IMPLICIT NONE

      CLASS(POLAR_DATA_T), INTENT(IN)  :: THIS
      TYPE(PRIM_CELL_T),   INTENT(IN)  :: PRIM_CELL

      INTEGER  :: N_PC, N_BORN

      IF (.NOT. ALLOCATED(THIS%BORN_EFF_CHARGES)) RETURN

      N_PC = PRIM_CELL%NUM_PRIM_ATOMS
      N_BORN = SIZE(THIS%BORN_EFF_CHARGES, 3)

      IF (N_BORN /= N_PC) CALL vtutor%error( &
         "The number of 3x3 matrices in the Born effective-charge tensor (" // str(N_BORN) // ") " // &
         "differs from the number of primitive-cell atoms (" // str(N_PC) // ").")
   END SUBROUTINE


!                            _______________________
!___________________________/ POLAR_EWALD_INIT_BASE \___________________________
!
!> @brief Initializes the Ewald parameter and other internals.
!>
!> G_CUTOFF is the radius of the sphere in reciprocal space that contains all
!> G-vectors used in the Ewald sum. Here, it is used to determine a suitable
!> Ewald parameter.
!>
!> The supplied threshold is the minimum value of the exponential Ewald
!> suppression factor before it is disregarded.
!>
!> @param[in,out] THIS     Instance reference
!> @param[in] G_CUTOFF     Cutoff in Cartesian G-space
!> @param[in] THRESHOLD    Ewald suppression threshold
!> @param[in] POLAR_DATA   Born effective charges and dielectric tensor
!_______________________________________________________________________________
   SUBROUTINE POLAR_EWALD_INIT_BASE(THIS, G_CUTOFF, THRESHOLD, POLAR_DATA)
      IMPLICIT NONE

      CLASS(POLAR_EWALD_T),   INTENT(INOUT)  :: THIS
      REAL(q),                INTENT(IN)     :: G_CUTOFF
      REAL(q),                INTENT(IN)     :: THRESHOLD
      TYPE(POLAR_DATA_T),     INTENT(IN)     :: POLAR_DATA

      THIS%G_CUTOFF = G_CUTOFF
      THIS%THRESHOLD = THRESHOLD
      THIS%DATA = POLAR_DATA
      THIS%EWALD_PARAM = GET_EWALD_PARAM(G_CUTOFF, THRESHOLD, THIS%DATA%DIELECTRIC_TENSOR)

      IF (ALLOCATED(THIS%BORN_FACTOR)) DEALLOCATE(THIS%BORN_FACTOR)
      ALLOCATE(THIS%BORN_FACTOR(3, SIZE(THIS%DATA%BORN_EFF_CHARGES)))
   END SUBROUTINE


!                           ________________________
!__________________________/ POLAR_EWALD_INIT_INDIV \___________________________
!
!> @brief Convenience overload for initializer.
!>
!> For details, look at the base initializer.
!>
!> @param[in,out] THIS           Instance reference
!> @param[in] G_CUTOFF           Cutoff in Cartesian G-space
!> @param[in] THRESHOLD          Ewald suppression threshold
!> @param[in] DIELECTRIC_TENSOR  Dielectric tensor
!> @param[in] BORN_EFF_CHARGES   Born effective charges
!_______________________________________________________________________________
   SUBROUTINE POLAR_EWALD_INIT_INDIV(THIS, G_CUTOFF, THRESHOLD, DIELECTRIC_TENSOR, BORN_EFF_CHARGES)
      IMPLICIT NONE

      CLASS(POLAR_EWALD_T),   INTENT(INOUT)  :: THIS
      REAL(q),                INTENT(IN)     :: G_CUTOFF
      REAL(q),                INTENT(IN)     :: THRESHOLD
      REAL(q),                INTENT(IN)     :: DIELECTRIC_TENSOR(3, 3)
      REAL(q),                INTENT(IN)     :: BORN_EFF_CHARGES(:, :, :)

      TYPE(POLAR_DATA_T)   :: COPY_DATA

      COPY_DATA%DIELECTRIC_TENSOR = DIELECTRIC_TENSOR
      COPY_DATA%BORN_EFF_CHARGES = BORN_EFF_CHARGES
      CALL THIS%INIT(G_CUTOFF, THRESHOLD, COPY_DATA)
   END SUBROUTINE


!                            ______________________
!___________________________/ POLAR_EWALD_TRY_CALC \____________________________
!
!> @brief Calculates Ewald-sum quantities for the supplied q- or G-vector.
!>
!> First, the Ewald suppression factor is calculated and compared against the
!> threshold. If it is too small, the routine stops and .false. is returned.
!> Otherwise, all quantities are calculated and .true. is returned.
!>
!> The calculated quantities are stored as member variables in the object.
!> Note that these quantities are not guaranteed to be properly set if this
!> routine returns falls, unless the optional FORCE parameter is set.
!>
!> @param[in,out] THIS  Instance reference
!> @param[in] QG_CART   q- or G-vector in Cartesian coordinates
!> @param[in] FORCE     If true, forces calculation of all quantities
!> @return              Whether the vector contributes to the Ewald sum
!_______________________________________________________________________________
   LOGICAL FUNCTION POLAR_EWALD_TRY_CALC(THIS, QG_CART, FORCE) RESULT(DOES_CONTRIBUTE)
      USE constant
      IMPLICIT NONE
      
      CLASS(POLAR_EWALD_T),   INTENT(INOUT)  :: THIS
      REAL(q),                INTENT(IN)     :: QG_CART(3)
      LOGICAL, OPTIONAL,      INTENT(IN)     :: FORCE

      INTEGER  :: P_IDX

      THIS%DIELEC_FACTOR = DOT_PRODUCT(QG_CART, MATMUL(THIS%DATA%DIELECTRIC_TENSOR, QG_CART))
      THIS%DIELEC_FACTOR = 1._q / (TPI * TPI * THIS%DIELEC_FACTOR)

      THIS%EXP_FACTOR = EXP(-1._q / (THIS%DIELEC_FACTOR * 4 * THIS%EWALD_PARAM**2))

      DOES_CONTRIBUTE = THIS%EXP_FACTOR >= THIS%THRESHOLD

      ! If below threshold, return early
      IF (PRESENT(FORCE)) THEN
         IF (.NOT. FORCE .AND. .NOT. DOES_CONTRIBUTE) RETURN
      ELSE IF (.NOT. DOES_CONTRIBUTE) THEN
         RETURN
      ENDIF

      DO P_IDX = 1, SIZE(THIS%BORN_FACTOR, 2)
         THIS%BORN_FACTOR(:, P_IDX) = TPI * MATMUL(QG_CART, THIS%DATA%BORN_EFF_CHARGES(:, :, P_IDX))
      ENDDO
   END FUNCTION


!                                ______________
!_______________________________/ GET_G_CUTOFF \________________________________
!
!> @brief Calculates G-space cutoff based on lattice geometry.
!>
!> The supplied cutoff is measured relative to the longest G-vector of the
!> supplied reciprocal lattice.
!>
!> @param[in] REL_G_CUTOFF Relative G-space cutoff
!> @param[in] G_LATT       Lattice (usually primitive cell)
!> @return                 Absolute G-space cutoff in Cartesian coordinates
!_______________________________________________________________________________
   REAL(q) PURE FUNCTION GET_G_CUTOFF(REL_G_CUTOFF, G_LATT) RESULT(G_CUTOFF)
      USE LATTICE,   ONLY: latt
      IMPLICIT NONE

      REAL(q),    INTENT(IN)  :: REL_G_CUTOFF
      TYPE(latt), INTENT(IN)  :: G_LATT

      INTEGER  :: I

      G_CUTOFF = 0
      DO I = 1, 3
         G_CUTOFF = MAX(G_CUTOFF, DOT_PRODUCT(G_LATT%B(:, I), G_LATT%B(:, I)))
      ENDDO
      G_CUTOFF = SQRT(G_CUTOFF) * REL_G_CUTOFF
   END FUNCTION


!                               _________________
!______________________________/ GET_EWALD_PARAM \______________________________
!
!> @brief Calculates a suitable Ewald parameter.
!>
!> The supplied G-space cutoff, \f$G_c\f$, being the radius of a sphere that
!> contains all the G-vectors used in the Ewald sum, should fulfill
!> \f[
!> e^\frac{-G \espilon G}{4 \alpha \alpha} = T
!> \f]
!> where \f$\alpha\f$ is the Ewald parameter, \f$\epsilon\f$ is the dielectric
!> constant and \f$T\f$ is the supplied threshold.
!> The optimal \f$\alpha\f$ is obtained by inverting this equation.
!>
!> @param[in] G_CUTOFF           G-space cutoff in Cartesian coordinates
!> @param[in] THRESHOLD          Dimensionless suppression threshold
!> @param[in] DIELECTRIC_TENSOR  Dielectric tensor
!> @return                       Suitable Ewald parameter
!_______________________________________________________________________________
   REAL(q) FUNCTION GET_EWALD_PARAM(G_CUTOFF, THRESHOLD, DIELECTRIC_TENSOR) RESULT(EWALD_PARAM)
      USE tutor,     ONLY: vtutor
      USE constant
      IMPLICIT NONE

      REAL(q), INTENT(IN)  :: G_CUTOFF
      REAL(q), INTENT(IN)  :: THRESHOLD
      REAL(q), INTENT(IN)  :: DIELECTRIC_TENSOR(3, 3)

      REAL(q)  :: DIELEC_COPY(3, 3)
      REAL(q)  :: EIGENVALS(3)
      REAL(q)  :: DIELEC_VAL
      REAL(q)  :: WORK(9)
      INTEGER  :: INFO

      ! Take maximum eigenvalue of dielectric tensor as measure
      CALL DSYEV("N", "U", 3, DIELECTRIC_TENSOR, 3, EIGENVALS, WORK, SIZE(WORK), INFO)
      DIELEC_VAL = MAXVAL(EIGENVALS)
      ! Fallback to maximum entry if diagonalization fails
      IF (INFO /= 0) DIELEC_VAL = MAXVAL(EIGENVALS)

      ! Determine Ewald parameter by inverting the exponential suppression
      EWALD_PARAM = PI * G_CUTOFF * SQRT(-MAXVAL(EIGENVALS) / (4 * LOG(THRESHOLD)))
   END FUNCTION


!                            _______________________
!___________________________/ SETUP_EWALD_G_VECTORS \___________________________
!
!> @brief Generate a list of all G-vectors to be used in the Ewald sum.
!>
!> The generated G-vectors are multiples of the reciprocal-lattice vectors of
!> the supplied lattice structure.
!> They are generated up to a maximum that is governed by the corresponding
!> exponential Ewald suppression factor.
!>
!> @param[in,out] POLAR_EWALD Polar Ewald structure
!> @param[in] GEN_LATTICE     Generator lattice defining the G-vectors
!> @param[out] G_VECTORS      List of G-vectors
!_______________________________________________________________________________
   SUBROUTINE GENERATE_EWALD_G_VECTORS(POLAR_EWALD, GEN_LATTICE, G_VECTORS)
      USE LATTICE,   ONLY: latt
      USE constant
      IMPLICIT NONE

      TYPE(POLAR_EWALD_T),    INTENT(INOUT)  :: POLAR_EWALD
      TYPE(latt),             INTENT(IN)     :: GEN_LATTICE
      REAL(q), ALLOCATABLE,   INTENT(OUT)    :: G_VECTORS(:, :)

      INTEGER              :: SHELL_SIZE
      INTEGER              :: X, Y, Z
      INTEGER              :: NUM_POINTS
      LOGICAL              :: FOUND_POINT
      REAL(q)              :: G_CART(3)
      INTEGER              :: G_POINT(3)
      INTEGER              :: ESTIMATED_POINTS

      ! Estimate the number of G vectors from the ratio between G-sphere and 1st BZ volumes
      ESTIMATED_POINTS = CEILING(4._q / 3._q * PI * POLAR_EWALD%G_CUTOFF**3 * GEN_LATTICE%OMEGA)

      ALLOCATE(G_VECTORS(3, ESTIMATED_POINTS))

      ! Add Gamma as first point
      G_VECTORS(:, 1) = [0, 0, 0]

      SHELL_SIZE = 0
      NUM_POINTS = 1
      FOUND_POINT = .TRUE.

      ! Keep increasing shell size until no more points contribute
      DO WHILE (FOUND_POINT)
         FOUND_POINT = .FALSE.
         SHELL_SIZE = SHELL_SIZE + 1

         ! Iterate over the current volume
         DO Z = -SHELL_SIZE, SHELL_SIZE
         DO Y = -SHELL_SIZE, SHELL_SIZE
         DO X = -SHELL_SIZE, SHELL_SIZE
            G_POINT = [X, Y, Z]

            ! Only take points on the shell
            IF (ANY(ABS(G_POINT) == SHELL_SIZE)) THEN
               ! Convert to Cartesian coordinates
               G_CART = MATMUL(GEN_LATTICE%B, REAL(G_POINT, q))

               ! Check if larger than threshold
               IF (POLAR_EWALD%TRY_CALC(G_CART)) THEN
                  NUM_POINTS = NUM_POINTS + 1
                  FOUND_POINT = .TRUE.

                  ! Grow G-vector array if needed
                  IF (NUM_POINTS > SIZE(G_VECTORS, 2)) &
                     CALL REALLOC_G_VECTORS(NUM_POINTS + MAX(31, INT(NUM_POINTS * 1.1)))

                  ! Add point
                  G_VECTORS(:, NUM_POINTS) = G_CART
               ENDIF
            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDDO

      ! Shrink G-vector array if needed
      IF (SIZE(G_VECTORS, 2) > NUM_POINTS) CALL REALLOC_G_VECTORS(NUM_POINTS)

      CONTAINS

      !> Simple helper routine to reallocate the G-vector array
      SUBROUTINE REALLOC_G_VECTORS(NEW_SIZE)
         IMPLICIT NONE

         INTEGER, INTENT(IN)  :: NEW_SIZE

         INTEGER              :: MIN_SIZE
         REAL(q), ALLOCATABLE :: TMP_ARRAY(:, :)

         MIN_SIZE = MIN(NEW_SIZE, SIZE(G_VECTORS, 2))

         ALLOCATE(TMP_ARRAY(3, NEW_SIZE))
         TMP_ARRAY(:, 1:MIN_SIZE) = G_VECTORS(:, 1:MIN_SIZE)
         CALL MOVE_ALLOC(TMP_ARRAY, G_VECTORS)
      END SUBROUTINE
   END SUBROUTINE


   SUBROUTINE SORT_PLOT_ARRAY(PLOT_ARR, LENGTH)
      IMPLICIT NONE

      REAL(q), INTENT(INOUT)  :: PLOT_ARR(:, :)
      INTEGER, INTENT(OUT)    :: LENGTH

      REAL(q), ALLOCATABLE :: REF_ARR(:, :)
      INTEGER              :: SORT_POS, I
      INTEGER              :: MIN_ELEMENT(1)
      LOGICAL, ALLOCATABLE :: MASK(:)

      ALLOCATE(REF_ARR(SIZE(PLOT_ARR, 1), 2), MASK(SIZE(PLOT_ARR, 1)))

      REF_ARR = PLOT_ARR
      MASK = .TRUE.
      SORT_POS = 0

      DO
         MIN_ELEMENT = MINLOC(REF_ARR(:, 1), MASK)
         IF (MIN_ELEMENT(1) == 0) EXIT

         DO I = 1, SIZE(REF_ARR, 1)
            IF (REF_ARR(I, 1) < REF_ARR(MIN_ELEMENT(1), 1) + 1E-5_q .AND. REF_ARR(I, 1) > REF_ARR(MIN_ELEMENT(1), 1) - 1E-5_q) THEN
               MASK(I) = .FALSE.
               IF (REF_ARR(I, 2) > REF_ARR(MIN_ELEMENT(1), 2)) THEN
                  MIN_ELEMENT(1) = I
               ENDIF
            ENDIF
         ENDDO

         SORT_POS = SORT_POS + 1
         PLOT_ARR(SORT_POS, :) = REF_ARR(MIN_ELEMENT(1), :)
      ENDDO

      LENGTH = SORT_POS

      DEALLOCATE(REF_ARR, MASK)
   END SUBROUTINE SORT_PLOT_ARRAY


   SUBROUTINE GET_MPI_STRIP(COMM, DATA_LENGTH, START, END)
      USE mpimy,  ONLY: communic
      IMPLICIT NONE

      TYPE(communic),   INTENT(IN)  :: COMM
      INTEGER,          INTENT(IN)  :: DATA_LENGTH
      INTEGER,          INTENT(OUT) :: START
      INTEGER,          INTENT(OUT) :: END

      INTEGER  :: STRIP, NSTRIPS
      INTEGER  :: SMALL_STRIP_LEN
      INTEGER  :: LARGE_STRIP_LEN
      INTEGER  :: NUM_LARGE_STRIPS
#ifdef MPI      
      NSTRIPS = COMM%NCPU
      STRIP = COMM%NODE_ME
#else
      NSTRIPS = 1
      STRIP = 1
#endif

      SMALL_STRIP_LEN = DATA_LENGTH / NSTRIPS
      LARGE_STRIP_LEN = SMALL_STRIP_LEN + 1
      NUM_LARGE_STRIPS = MODULO(DATA_LENGTH, NSTRIPS)

      IF (NUM_LARGE_STRIPS == 0) THEN
         END = STRIP * SMALL_STRIP_LEN
         START = END - SMALL_STRIP_LEN + 1
      ELSE IF (STRIP <= NUM_LARGE_STRIPS) THEN
         END = STRIP * LARGE_STRIP_LEN
         START = END - LARGE_STRIP_LEN + 1
      ELSE
         END = NUM_LARGE_STRIPS * LARGE_STRIP_LEN + (STRIP - NUM_LARGE_STRIPS) * SMALL_STRIP_LEN
         START = END - SMALL_STRIP_LEN + 1
      ENDIF

   END SUBROUTINE GET_MPI_STRIP


   SUBROUTINE WRITE_MATRIX_REAL(MAT, FORMAT_STR, FILE_NAME, IO)
      IMPLICIT NONE

      REAL(q),          INTENT(IN)  :: MAT(:, :)
      CHARACTER(*),     INTENT(IN)  :: FORMAT_STR
      CHARACTER(*),     INTENT(IN)  :: FILE_NAME
      TYPE(in_struct),  INTENT(IN)  :: IO

      INTEGER  :: UNIT
      INTEGER  :: I, J

      IF (IO%IU0 < 0) RETURN

      OPEN(NEWUNIT = UNIT, FILE = FILE_NAME, STATUS = 'replace')

      DO I = 1, SIZE(MAT, 1)
         DO J = 1, SIZE(MAT, 2)
            WRITE(UNIT, FORMAT_STR, ADVANCE = 'no') MAT(I, J)
            WRITE(UNIT, '(A)', ADVANCE = 'no') '   '
         ENDDO
         WRITE(UNIT, *)
      ENDDO
      CLOSE(UNIT)

   END SUBROUTINE WRITE_MATRIX_REAL


   SUBROUTINE WRITE_MATRIX_COMPLEX(MAT, FORMAT_STR, FILE_NAME, IO)
      IMPLICIT NONE

      COMPLEX(q),       INTENT(IN)  :: MAT(:, :)
      CHARACTER(*),     INTENT(IN)  :: FORMAT_STR
      CHARACTER(*),     INTENT(IN)  :: FILE_NAME
      TYPE(in_struct),  INTENT(IN)  :: IO

      INTEGER  :: UNIT_RE, UNIT_IM
      INTEGER  :: I, J

      IF (IO%IU0 < 0) RETURN

      OPEN(NEWUNIT = UNIT_RE, FILE = FILE_NAME // '_RE', STATUS = 'replace')
      OPEN(NEWUNIT = UNIT_IM, FILE = FILE_NAME // '_IM', STATUS = 'replace')

      DO I = 1, SIZE(MAT, 1)
         DO J = 1, SIZE(MAT, 2)
            WRITE(UNIT_RE, FORMAT_STR, ADVANCE = 'no') REAL(MAT(I, J))
            WRITE(UNIT_RE, '(A)', ADVANCE = 'no') '   '

            WRITE(UNIT_IM, FORMAT_STR, ADVANCE = 'no') AIMAG(MAT(I, J))
            WRITE(UNIT_IM, '(A)', ADVANCE = 'no') '   '
         ENDDO
         WRITE(UNIT_RE, *)
         WRITE(UNIT_IM, *)
      ENDDO
      CLOSE(UNIT_RE)
      CLOSE(UNIT_IM)

   END SUBROUTINE WRITE_MATRIX_COMPLEX

   PURE FUNCTION DYAD(A, B) RESULT(MAT)
      IMPLICIT NONE

      REAL(q), INTENT(IN)  :: A(:)
      REAL(q), INTENT(IN)  :: B(:)
      REAL(q)              :: MAT(SIZE(A), SIZE(B))

      INTEGER :: I, J

      DO J = 1, SIZE(B)
         DO I = 1, SIZE(A)
            MAT(I, J) = A(I) * B(J)
         ENDDO
      ENDDO
   END FUNCTION DYAD

!                                   _________
!__________________________________/ DET_3X3 \__________________________________
!
!> @brief Calculates determinant of a real 3x3 matrix.
!>
!> @param[in] A   Matrix
!> @return        Determinant
!_______________________________________________________________________________
   PURE FUNCTION DET_3X3(A) RESULT(DET)
      IMPLICIT NONE

      REAL(q), INTENT(IN)  :: A(3, 3)
      REAL(q)              :: DET

      DET   = A(1, 1) * A(2, 2) * A(3, 3) &
            + A(2, 1) * A(3, 2) * A(1, 3) &
            + A(1, 2) * A(2, 3) * A(3, 1) &
            - A(1, 3) * A(2, 2) * A(3, 1) &
            - A(2, 1) * A(1, 2) * A(3, 3) &
            - A(3, 2) * A(2, 3) * A(1, 1)
   ENDFUNCTION DET_3X3

!                             ____________________
!____________________________/ MATRIX_INVERSE_3X3 \_____________________________
!
!> @brief Calculates the inverse of a real 3x3 matrix.
!>
!> @param[in] A   Matrix
!> @return        Matrix inverse
!_______________________________________________________________________________
   PURE FUNCTION MATRIX_INVERSE_3X3(A) RESULT(B)
      IMPLICIT NONE

      REAL(q), INTENT(IN)  :: A(3, 3)
      REAL(q)              :: B(3, 3)
      REAL(q)              :: DET_INV

      ! Calculate inverse determinant.
      DET_INV = 1._q / DET_3X3(A)

      ! Calculate inverse of the matrix.
      B(1, 1) = +DET_INV * (A(2, 2) * A(3, 3) - A(2, 3) * A(3, 2))
      B(2, 1) = -DET_INV * (A(2, 1) * A(3, 3) - A(2, 3) * A(3, 1))
      B(3, 1) = +DET_INV * (A(2, 1) * A(3, 2) - A(2, 2) * A(3, 1))
      B(1, 2) = -DET_INV * (A(1, 2) * A(3, 3) - A(1, 3) * A(3, 2))
      B(2, 2) = +DET_INV * (A(1, 1) * A(3, 3) - A(1, 3) * A(3, 1))
      B(3, 2) = -DET_INV * (A(1, 1) * A(3, 2) - A(1, 2) * A(3, 1))
      B(1, 3) = +DET_INV * (A(1, 2) * A(2, 3) - A(1, 3) * A(2, 2))
      B(2, 3) = -DET_INV * (A(1, 1) * A(2, 3) - A(1, 3) * A(2, 1))
      B(3, 3) = +DET_INV * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1))
   END FUNCTION MATRIX_INVERSE_3X3

ENDMODULE





MODULE phonon
   USE prec
   USE base
   USE elphon_util
   USE tutor,              ONLY: vtutor
   USE string,             ONLY: str
   USE msymmetry,          ONLY: PRIM_CELL_T
   IMPLICIT NONE
   
   PRIVATE

   PUBLIC PHONON_SETTINGS, PHON_READER
   PUBLIC PHON_INIT, PHON_CALC_PHONONS, PHON_CALC_DISPERSION
   PUBLIC PHON_FIND_TRANSLATIONAL_MODES, PHON_INIT_POLAR

#ifdef VASP_HDF5
   PUBLIC PHON_WRITE_HDF5, PHON_READ_HDF5
#endif

   TYPE PHONON_SETTINGS
      LOGICAL  :: CALC_DISPERSION = .FALSE.        !< Calculate phonon frequencies and eigenvectors?
      LOGICAL  :: DO_POLAR = .FALSE.               !< Include dipole-dipole corrections?
      LOGICAL  :: READ_FORCE_CONSTANTS=.FALSE.     !< Read force constants from hdf5 file
      INTEGER  :: NWRITE = 1                       !< Determines how much output is written to OUTCAR
      REAL(q)  :: POLAR_G_CUTOFF = 8._q            !< Cutoff in reciprocal space (multiple of longest G-vector)
      REAL(q)  :: POLAR_EWALD_THRESHOLD = 1E-10_q  !< Threshold down to which Ewald summands contribute

      CONTAINS

      PROCEDURE,  PUBLIC   :: DO_INIT => SETTINGS_DO_INIT
      PROCEDURE,  PUBLIC   :: DO_CALC => SETTINGS_DO_CALC
   END TYPE

   REAL(q), ALLOCATABLE, SAVE             :: FORCE_CONST(:, :, :, :, :)       ! Force-constant matrix (1:3, 1:3, cell_idx_1, prim_idx_1, prim_idx_2)
   REAL(q), ALLOCATABLE, SAVE             :: FORCE_CONST_SHORT(:, :, :, :, :) ! Short-ranged force constants

   LOGICAL, SAVE                          :: DO_POLAR = .FALSE.
   TYPE(POLAR_EWALD_T), SAVE              :: POLAR_EWALD                      ! Facilitates Ewald sum for polar contribution
   REAL(q), ALLOCATABLE, SAVE             :: EWALD_G_VECTORS(:, :)            ! List of G-vectors for primitive-cell Ewald sum in Cartesian coordinates

CONTAINS


!                                 _____________
!________________________________/ PHON_READER \________________________________
!
!> @brief Reads settings from the INCAR file.
!>
!> # The following INCAR tags are available:
!>
!> ## LPHON_DISPERSION
!> Requests the calculation of the phonon dispersion along the q-point path
!> supplied in file QPOINTS (same format as KPOINTS).
!>
!> ## LPHON_POLAR
!> Includes dipole-dipole corrections. For this mode PHON_DIELECTRIC,
!> PHON_BORN_CHARGES and PHON_ATOM_INDICES must also be set.
!>
!> ## PHON_DIELECTRIC
!> Set dielectric tensor (3x3 values). The tensor is read in row by row.
!>
!> ## PHON_BORN_CHARGES
!> Set Born effective charges (3 x 3 x P_NUM values, P_NUM = number of atoms in
!> primitive cell). The tensor is read in row by row.
!>
!> ## PHON_G_CUTOFF
!> Cutoff radius in reciprocal space used to determine the number of G-vectors
!> involved in the Ewald sum in polar phonon calculations. It is specified as a
!> multiple of the longest reciprocal-lattice vector corresponding to the first
!> Brillouin zone of the primitive cell. This way, the Ewald sum is as
!> independent of the system as possible.
!>
!> ## PHON_NWRITE
!> Determines how much output is written to the OUTCAR file for phonon-
!> dispersion calculations. Positive numbers mean human readable output,
!> negative numbers writes in one-line format.
!>   2 = For each q-point, write the same as 1 and then the phonon modes
!>       with the displacement of each atom in the three cartesian directions per line
!>   1 = For each q-point, q-point coordinates are written in one line and
!>       the phonon frequencies are written one branch per line in different units.
!>   0 = No phonon output is written to OUTCAR
!>  -1 = For each q-point, only a single line is written containing q-point
!>       coordinates and frequencies.
!>  -2 = For each q-point, q-point coordinates and frequencies are written in
!>       separate blocks and frequencies are reported in different units.
!>  -3 = Like -2, but in addition, the phonon eigenvectors are written for each
!>       q-point.
!_______________________________________________________________________________
   SUBROUTINE PHON_READER(INCAR, SETTINGS, POLAR_DATA)
      USE incar_reader, ONLY: incar_file, &
                              process_incar, &
                              count_elements
      USE string,       ONLY: lowercase
      IMPLICIT NONE

      TYPE(incar_file),       INTENT(INOUT)  :: INCAR
      TYPE(PHONON_SETTINGS),  INTENT(OUT)    :: SETTINGS
      TYPE(POLAR_DATA_T),     INTENT(OUT)    :: POLAR_DATA

      INTEGER                    :: IERR
      INTEGER                    :: CHECK_IBRION
      INTEGER                    :: NUM_PRIM_ATOMS
      INTEGER                    :: CHECK_NUM
      INTEGER                    :: P_IDX
      CHARACTER(:), ALLOCATABLE  :: READ_STRING
      LOGICAL                    :: LOAD_DATA

      ! Calculate phonon dispersion?
      CALL process_incar(INCAR, "LPHON_DISPERSION", SETTINGS%CALC_DISPERSION)

      ! Include corrections for polar materials?
      CALL process_incar(INCAR, "LPHON_POLAR", SETTINGS%DO_POLAR)

      ! Read force constants from hdf5 file?
      CALL process_incar(INCAR, "LPHON_READ_FORCE_CONSTANTS", SETTINGS%READ_FORCE_CONSTANTS)

      ! Try to read the cut-off radius
      CALL process_incar(INCAR, "PHON_G_CUTOFF", SETTINGS%POLAR_G_CUTOFF)

      ! How much is written to OUTCAR?
      CALL process_incar(INCAR, "PHON_NWRITE", SETTINGS%NWRITE)

      ! Get number of elements in dielectric tensor
      CHECK_NUM = count_elements(INCAR, "PHON_DIELECTRIC")

      IF (CHECK_NUM > 0) THEN
         ! Check shape
         IF (CHECK_NUM /= 9) CALL vtutor%error('The dielectric tensor supplied via PHON_DIELECTRIC has the wrong number ' // &
            'of components: ' // str(CHECK_NUM) // '. Must be 9 (3x3).')

         ! Read dielectric tensor
         CALL process_incar(INCAR, "PHON_DIELECTRIC", POLAR_DATA%DIELECTRIC_TENSOR)
         POLAR_DATA%DIELECTRIC_TENSOR = TRANSPOSE(POLAR_DATA%DIELECTRIC_TENSOR)
      ENDIF

      ! Get number of elements in Born-effective-charge tensor
      CHECK_NUM = count_elements(INCAR, "PHON_BORN_CHARGES")

      IF (CHECK_NUM > 0) THEN
         ! Check shape
         NUM_PRIM_ATOMS = CHECK_NUM / 9
         IF (CHECK_NUM > 0 .AND. CHECK_NUM /= NUM_PRIM_ATOMS * 9) CALL vtutor%error('The Born effective-charge tensor supplied via PHON_BORN_CHARGES ' // &
            'has a number of components (' // str(CHECK_NUM) // ') that is not divisible by 9 (3x3 matrix for each primitive-cell atom).')

         ! Allocate and read Born-effective-charge tensor
         ALLOCATE(POLAR_DATA%BORN_EFF_CHARGES(3, 3, NUM_PRIM_ATOMS))
         CALL process_incar(INCAR, "PHON_BORN_CHARGES", POLAR_DATA%BORN_EFF_CHARGES)
         DO P_IDX = 1, NUM_PRIM_ATOMS
            POLAR_DATA%BORN_EFF_CHARGES(:, :, P_IDX) = TRANSPOSE(POLAR_DATA%BORN_EFF_CHARGES(:, :, P_IDX))
         ENDDO
      ENDIF
   END SUBROUTINE


   PURE LOGICAL FUNCTION SETTINGS_DO_INIT(THIS) RESULT(D)
      IMPLICIT NONE
      CLASS(PHONON_SETTINGS), INTENT(IN)  :: THIS
      D = THIS%CALC_DISPERSION
   END FUNCTION

   PURE LOGICAL FUNCTION SETTINGS_DO_CALC(THIS) RESULT(D)
      IMPLICIT NONE
      CLASS(PHONON_SETTINGS), INTENT(IN)  :: THIS
      D = THIS%CALC_DISPERSION
   END FUNCTION


   SUBROUTINE PHON_INIT(PRIM_CELL, FORCE_CONSTANTS, IO)
      IMPLICIT NONE

      TYPE(PRIM_CELL_T),   INTENT(IN)  :: PRIM_CELL
      REAL(q),             INTENT(IN)  :: FORCE_CONSTANTS(:, :)
      TYPE(in_struct),     INTENT(IN)  :: IO

      INTEGER              :: P_IDX_1, P_IDX_2
      INTEGER              :: S_IDX_1, S_IDX_2
      INTEGER              :: CELL_IDX_1
      REAL(q), ALLOCATABLE :: TMP_ARR(:, :, :)
      LOGICAL, SAVE        :: IS_INIT = .FALSE.

      IF (IS_INIT) CALL vtutor%bug('PHON_INIT called twice!', __FILE__, __LINE__)

      IF (ALLOCATED(FORCE_CONST)) DEALLOCATE(FORCE_CONST)
      ALLOCATE(FORCE_CONST(3, 3, PRIM_CELL%NUM_CELLS, PRIM_CELL%NUM_PRIM_ATOMS, PRIM_CELL%NUM_PRIM_ATOMS))

      ! Reorder force constants to optimized memory layout
      DO P_IDX_2 = 1, PRIM_CELL%NUM_PRIM_ATOMS
         S_IDX_2 = PRIM_CELL%SUPER_IDX(1, P_IDX_2)
         DO P_IDX_1 = 1, PRIM_CELL%NUM_PRIM_ATOMS
            DO CELL_IDX_1 = 1, PRIM_CELL%NUM_CELLS
               S_IDX_1 = PRIM_CELL%SUPER_IDX(CELL_IDX_1, P_IDX_1)

               ! _test_mani mass weigh
               FORCE_CONST(:, :, CELL_IDX_1, P_IDX_1, P_IDX_2) = &
               FORCE_CONSTANTS(3 * S_IDX_1 - 2 : 3 * S_IDX_1, 3 * S_IDX_2 - 2 : 3 * S_IDX_2)
            ENDDO
         ENDDO
      ENDDO

      CALL CORRECT_FORCE_CONSTANTS()

      IS_INIT = .TRUE.
   END SUBROUTINE PHON_INIT


!                           _________________________
!__________________________/ CORRECT_FORCE_CONSTANTS \__________________________
!
!> @brief Enforces the acoustic sum rule on the force constants.
!>
!> This is required to obtain correct acoustic frequencies at and close to the
!> Gamma point.
!_______________________________________________________________________________
   SUBROUTINE CORRECT_FORCE_CONSTANTS()
      IMPLICIT NONE

      INTEGER  :: NUM_PRIM, P_IDX_1, P_IDX_2
      INTEGER  :: NUM_CELLS, CELL_IDX_1

      NUM_PRIM = SIZE(FORCE_CONST, 5)
      NUM_CELLS = SIZE(FORCE_CONST, 3)

      DO P_IDX_2 = 1, NUM_PRIM
         FORCE_CONST(:, :, 1, P_IDX_2, P_IDX_2) = 0
         DO P_IDX_1 = 1, NUM_PRIM
            DO CELL_IDX_1 = 1, NUM_CELLS
               IF (P_IDX_1 /= P_IDX_2 .OR. CELL_IDX_1 /= 1) THEN
                  FORCE_CONST(:, :, 1, P_IDX_2, P_IDX_2) = &
                  FORCE_CONST(:, :, 1, P_IDX_2, P_IDX_2) - &
                  FORCE_CONST(:, :, CELL_IDX_1, P_IDX_1, P_IDX_2)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   END SUBROUTINE CORRECT_FORCE_CONSTANTS


!                               _________________
!______________________________/ PHON_INIT_POLAR \______________________________
!
!> @brief Initializes the polar phonon code.
!>
!> Sets up short-ranged force constants for later use in Fourier interpolation,
!> together with the relevant G-vectors used in the Ewald sum.
!>
!> @param[in] SETTINGS        Electron-phonon settings
!> @param[in,out] PRIM_CELL   Primitive-cell information
!> @param[in] POLAR_DATA      Born effective charges and dielectric tensor
!_______________________________________________________________________________
   SUBROUTINE PHON_INIT_POLAR(SETTINGS, PRIM_CELL, POLAR_DATA)
      IMPLICIT NONE

      TYPE(PHONON_SETTINGS),  INTENT(IN)     :: SETTINGS
      TYPE(PRIM_CELL_T),      INTENT(INOUT)  :: PRIM_CELL
      TYPE(POLAR_DATA_T),     INTENT(IN)     :: POLAR_DATA

      REAL(q)  :: ABS_G_CUTOFF

      ! Setup the minimum-distance vectors in the primitive cell if needed
      CALL PRIM_CELL%SETUP_ADDITIONAL_DATA()

      ABS_G_CUTOFF = GET_G_CUTOFF(SETTINGS%POLAR_G_CUTOFF, PRIM_CELL%PRIM_LATT)
      CALL POLAR_EWALD%INIT(ABS_G_CUTOFF, SETTINGS%POLAR_EWALD_THRESHOLD, POLAR_DATA)
      CALL GENERATE_EWALD_G_VECTORS(POLAR_EWALD, PRIM_CELL%PRIM_LATT, EWALD_G_VECTORS)
      CALL SETUP_FORCE_CONST_SHORT(PRIM_CELL)
      DO_POLAR = .TRUE.
   END SUBROUTINE


!                              ___________________
!_____________________________/ PHON_CALC_PHONONS \_____________________________
!
!> @brief Calculates phonon frequencies and modes.
!>
!> Performs a Fourier interpolation of the force-constant matrix to the
!> primitive cell, followed by a diagonalization of the dynamical matrix.
!> Imaginary frequencies are returned as negative values.
!>
!> @param[in] Q_VEC     Arbitrary vector in direct coordinates of the Brillouin
!>                      zone of the primitive cell
!> @param[in] PRIM_CELL Primitive-cell information
!> @param[in] T_INFO    Atom-type information
!> @param[out] FREQU    Phonon angular frequencies in internal units
!> @param[out] MODES    Normalized phonon modes (columns)
!_______________________________________________________________________________
   SUBROUTINE PHON_CALC_PHONONS(Q_VEC, PRIM_CELL, T_INFO, FREQU, MODES, DIRECTION)
      USE poscar_struct_def,  ONLY: type_info
      IMPLICIT NONE

      REAL(q),                   INTENT(IN)  :: Q_VEC(3)
      TYPE(PRIM_CELL_T),         INTENT(IN)  :: PRIM_CELL
      TYPE(type_info),           INTENT(IN)  :: T_INFO
      REAL(q), ALLOCATABLE,      INTENT(OUT) :: FREQU(:)
      COMPLEX(q), ALLOCATABLE,   INTENT(OUT) :: MODES(:, :)
      REAL(q), OPTIONAL,         INTENT(IN)  :: DIRECTION(3)

      ! locals
      COMPLEX(q), ALLOCATABLE    :: DYN_MAT(:, :, :, :)
      REAL(q), ALLOCATABLE       :: MASS_FACTOR(:)
      REAL(q)                    :: MY_DIRECTION(3)
      INTEGER                    :: P_IDX_1, P_IDX_2, P_NUM

      ! work arrays for ZHEEV
      INTEGER                 :: LWORK, IFAIL
      REAL(q)                 :: QNORM
      COMPLEX(q), ALLOCATABLE :: WORK(:)
      REAL(q), ALLOCATABLE    :: RWORK(:)

      MY_DIRECTION = [REAL(q):: 1, 0, 0]
      IF (PRESENT(DIRECTION)) THEN
         QNORM=SQRT(DIRECTION(1)**2+DIRECTION(2)**2+DIRECTION(3)**2)
         IF (QNORM>1e-6) MY_DIRECTION = DIRECTION/QNORM
      ENDIF
      IF (ALLOCATED(FREQU)) DEALLOCATE(FREQU)
      IF (ALLOCATED(MODES)) DEALLOCATE(MODES)

      P_NUM = PRIM_CELL%NUM_PRIM_ATOMS

      ALLOCATE(DYN_MAT(3, 3, P_NUM, P_NUM))

      ! Calculate dynamical matrix.
      IF (DO_POLAR) THEN

         ! Check if the polar code has been initialized.
         IF (.NOT. ALLOCATED(FORCE_CONST_SHORT)) CALL vtutor%bug("FORCE_CONST_SHORT not allocated!", __FILE__, __LINE__)

         ! Interpolate from short-ranged force constants.
         DYN_MAT = PHON_CALC_DYN_MAT_FROM_FC(PRIM_CELL, FORCE_CONST_SHORT, Q_VEC)
         ! Add the long-ranged part back up.
         DYN_MAT = DYN_MAT + CALC_DYN_MAT_LONG(PRIM_CELL, Q_VEC)
         ! Add the Gamma-point correction.
         IF (ALL(ABS(Q_VEC) < 1E-10_q)) THEN
            DYN_MAT = DYN_MAT + CALC_GAMMA_CORRECTION(PRIM_CELL, MY_DIRECTION)
         ENDIF
      ELSE
         ! Interpolate directly from the ab-initio force constants.
         DYN_MAT = PHON_CALC_DYN_MAT_FROM_FC(PRIM_CELL, FORCE_CONST, Q_VEC)
      ENDIF

      ! Set mass factors.
      ALLOCATE(MASS_FACTOR(P_NUM))
      DO P_IDX_1 = 1, P_NUM
         MASS_FACTOR(P_IDX_1) = 1._q / SQRT(T_INFO%POMASS(T_INFO%ITYP(PRIM_CELL%SUPER_IDX(1, P_IDX_1))))
      ENDDO

      ! Divide by square root of masses.
      DO P_IDX_2 = 1, P_NUM
         DO P_IDX_1 = 1, P_NUM
            DYN_MAT(:, :, P_IDX_1, P_IDX_2) = &
            DYN_MAT(:, :, P_IDX_1, P_IDX_2) * &
            MASS_FACTOR(P_IDX_1) * MASS_FACTOR(P_IDX_2)
         ENDDO
      ENDDO

      ALLOCATE(FREQU(3 * P_NUM), MODES(3 * P_NUM, 3 * P_NUM))

      ! Reorder array for diagonalization.
      DO P_IDX_2 = 1, P_NUM
         DO P_IDX_1 = 1, P_NUM
            MODES(3*P_IDX_1-2 : 3*P_IDX_1, 3*P_IDX_2-2 : 3*P_IDX_2) = DYN_MAT(:, :, P_IDX_1, P_IDX_2)
         ENDDO
      ENDDO

      DEALLOCATE(DYN_MAT)

      ! Allocate work arrays.
      LWORK = 6 * P_NUM - 1
      ALLOCATE(WORK(LWORK), RWORK(9 * P_NUM - 2))

      ! Calculate eigenvalues and eigenvectors
      CALL ZHEEV( &
      &  'V', 'U', 3 * P_NUM, MODES, 3 * P_NUM, &
      &  FREQU, WORK, LWORK, RWORK, IFAIL)

      DEALLOCATE(WORK, RWORK)

      IF (IFAIL /= 0) THEN
         CALL vtutor%bug("Internal error in PHON_CALC_PHONONS: Call to ZHEEV failed! Error code was " &
            // str(IFAIL), __FILE__, __LINE__)
      ENDIF

      ! Imaginary frequencies are represented as negative numbers.
      FREQU = SIGN(SQRT(ABS(FREQU)), FREQU)

! _test_mani Is this necessary?
! do i = 1, 3 * P_NUM
!    j = ((i - 1) / 3) * 3 + 1
!    if (abs(modes(j, i)) /= 0._q) then
!       modes(:, i) = modes(:, i) * exp((.0_q, -1._q) * atan2(aimag(modes(j, i)), real(modes(j, i))))
!    endif
! enddo

   END SUBROUTINE PHON_CALC_PHONONS


!                            _______________________
!___________________________/ CALC_GAMMA_CORRECTION \___________________________
!
!> @brief Calculates dipole-dipole correction at the Gamma point.
!>
!> The long-ranged dipole-dipole term is not defined *at* Gamma, so the best we
!> can do is to calculate the correction at an infinitesimal distance away from
!> Gamma in a certain direction.
!> The returned dynamic matrix is not mass weighted!
!>
!> @param[in] PRIM_CELL Primitive-cell information
!> @param[in] G_VEC     Direction vector (direct coordinates)
!> @return DYN_MAT      Dipole-dipole dynamical matrix at Gamma
!_______________________________________________________________________________
   FUNCTION CALC_GAMMA_CORRECTION(PRIM_CELL, G_VEC) RESULT(DYN_MAT)
      USE constant
      IMPLICIT NONE

      REAL(q), ALLOCATABLE :: DYN_MAT(:, :, :, :)

      TYPE(PRIM_CELL_T),   INTENT(IN)  :: PRIM_CELL
      REAL(q),             INTENT(IN)  :: G_VEC(3)

      REAL(q)  :: G_CART(3)
      INTEGER  :: CART_IDX_1, CART_IDX_2
      INTEGER  :: P_NUM
      INTEGER  :: P_IDX_1, P_IDX_2
      LOGICAL  :: L_TMP

      P_NUM = PRIM_CELL%NUM_PRIM_ATOMS

      ! Convert to Cartesian coordinates
      G_CART = MATMUL(PRIM_CELL%PRIM_LATT%B, G_VEC)

      ALLOCATE(DYN_MAT(3, 3, P_NUM, P_NUM))
      DYN_MAT = 0

      L_TMP = POLAR_EWALD%TRY_CALC(G_CART, .TRUE.)

      ! Calculate tensor product
      DO P_IDX_2 = 1, P_NUM
         DO P_IDX_1 = 1, P_NUM
            DO CART_IDX_2 = 1, 3
               DO CART_IDX_1 = 1, 3
                  DYN_MAT(CART_IDX_1, CART_IDX_2, P_IDX_1, P_IDX_2) = &
                  POLAR_EWALD%BORN_FACTOR(CART_IDX_1, P_IDX_1) * &
                  POLAR_EWALD%BORN_FACTOR(CART_IDX_2, P_IDX_2)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      ! Final result + correct units
      DYN_MAT = DYN_MAT * POLAR_EWALD%DIELEC_FACTOR
      DYN_MAT = DYN_MAT * 4 * PI / PRIM_CELL%PRIM_LATT%OMEGA * FELECT
   END FUNCTION CALC_GAMMA_CORRECTION


!                              ___________________
!_____________________________/ CALC_DYN_MAT_LONG \_____________________________
!
!> @brief Calculates the long-ranged dipole-dipole dynamical matrix.
!>
!> The long-ranged dipole-dipole contribution is calculated via an analytic
!> expression at the specified q-vector using an Ewald sum.
!> The returned dynamical matrix is not mass weighted!
!>
!> @param[in] PRIM_CELL Primitive-cell information
!> @param[in] Q_VEC     Arbitrary vector in direct coordinates of the Brillouin
!>                      zone of the primitive cell
!> @return DYN_MAT      Dipole-dipole dynamical matrix
!_______________________________________________________________________________
   FUNCTION CALC_DYN_MAT_LONG(PRIM_CELL, Q_VEC) RESULT(DYN_MAT)
      USE constant
      IMPLICIT NONE

      COMPLEX(q), ALLOCATABLE :: DYN_MAT(:, :, :, :)

      TYPE(PRIM_CELL_T),   INTENT(IN)  :: PRIM_CELL
      REAL(q),             INTENT(IN)  :: Q_VEC(3)

      REAL(q)              :: Q_CART(3), QG_CART(3)
      INTEGER              :: G_IDX
      INTEGER              :: P_NUM
      INTEGER              :: P_IDX_1, P_IDX_2
      INTEGER              :: CART_IDX_1, CART_IDX_2
      COMPLEX(q)           :: PHASE

      P_NUM = PRIM_CELL%NUM_PRIM_ATOMS

      IF (ALLOCATED(DYN_MAT)) DEALLOCATE(DYN_MAT)
      ALLOCATE(DYN_MAT(3, 3, P_NUM, P_NUM))

      DYN_MAT = 0

      ! Convert to Carthesian coordinates
      Q_CART = MATMUL(PRIM_CELL%PRIM_LATT%B, Q_VEC)

      ! Loop over all G-vectors inside the Ewald cut-off
      DO G_IDX = 1, SIZE(EWALD_G_VECTORS, 2)
         ! G already in Carthesian coordinates!
         QG_CART = Q_CART + EWALD_G_VECTORS(:, G_IDX)

         ! Only take non-zero q + G
         IF (ALL(ABS(QG_CART) < 1E-10_q)) CYCLE

         ! Skip if Ewald suppression factor is too small
         IF (.NOT. POLAR_EWALD%TRY_CALC(QG_CART)) CYCLE

         ! Double loop over primitive atoms
         DO P_IDX_2 = 1, P_NUM
            DO P_IDX_1 = 1, P_NUM
               PHASE = PRIM_CELL%MIN_IMAGES(1, P_IDX_1, P_IDX_2)%CALC_PRIM_PHASE(QG_CART)

               ! Add/Multiply everything up
               DO CART_IDX_2 = 1, 3
                  DO CART_IDX_1 = 1, 3
                     DYN_MAT(CART_IDX_1, CART_IDX_2, P_IDX_1, P_IDX_2) = &
                     DYN_MAT(CART_IDX_1, CART_IDX_2, P_IDX_1, P_IDX_2) + &
                     POLAR_EWALD%BORN_FACTOR(CART_IDX_1, P_IDX_1) * &
                     POLAR_EWALD%BORN_FACTOR(CART_IDX_2, P_IDX_2) * &
                     POLAR_EWALD%DIELEC_FACTOR * POLAR_EWALD%EXP_FACTOR * PHASE
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      ! Correct units
      DYN_MAT = DYN_MAT * 4._q * PI / PRIM_CELL%PRIM_LATT%OMEGA * FELECT

   END FUNCTION CALC_DYN_MAT_LONG


!                          ___________________________
!_________________________/ PHON_CALC_DYN_MAT_FROM_FC \_________________________
!
!> @brief Calculates the dynamical matrix via Fourier interpolation.
!>
!> Fourier interpolates the matrix of inter-atomic force constants to the
!> primitive cell at the given phonon wave vector.
!> The returned dynamical matrix is not mass weighted!
!>
!> @param[in] PRIM_CELL    Primitive-cell information
!> @param[in] FORCE_CONST  Force-constant matrix in supercell
!> @param[in] Q_VEC        Arbitrary vector in direct coordinates of the
!>                         Brillouin zone of the primitive cell
!> @return DYN_MAT         Dynamical matrix
!_______________________________________________________________________________
   FUNCTION PHON_CALC_DYN_MAT_FROM_FC(PRIM_CELL, FORCE_CONST, Q_VEC) RESULT(DYN_MAT)
      USE constant
      IMPLICIT NONE

      COMPLEX(q), ALLOCATABLE :: DYN_MAT(:, :, :, :)

      TYPE(PRIM_CELL_T),   INTENT(IN)  :: PRIM_CELL
      REAL(q),             INTENT(IN)  :: FORCE_CONST(:, :, :, :, :)
      REAL(q),             INTENT(IN)  :: Q_VEC(3)

      INTEGER     :: P_NUM
      INTEGER     :: P_IDX_1, P_IDX_2
      INTEGER     :: CELL_IDX_1
      REAL(q)     :: Q_CART(3)
      COMPLEX(q)  :: PHASE

      P_NUM = PRIM_CELL%NUM_PRIM_ATOMS

      IF (ALLOCATED(DYN_MAT)) DEALLOCATE(DYN_MAT)
      ALLOCATE(DYN_MAT(3, 3, P_NUM, P_NUM))
      DYN_MAT = 0

      ! Convert to Carthesian coordinates.
      Q_CART = MATMUL(PRIM_CELL%PRIM_LATT%B, Q_VEC)

      ! Double loop over primitive atoms.
      DO P_IDX_2 = 1, P_NUM
         DO P_IDX_1 = 1, P_NUM
            ! Loop over cells.
            DO CELL_IDX_1 = 1, PRIM_CELL%NUM_CELLS
               PHASE = PRIM_CELL%MIN_IMAGES(CELL_IDX_1, P_IDX_1, P_IDX_2)%CALC_CELL_PHASE(-Q_CART)

               ! Interpolate each 3x3 force-constant matrix to q.
               DYN_MAT(:, :, P_IDX_1, P_IDX_2) = &
               DYN_MAT(:, :, P_IDX_1, P_IDX_2) + &
               FORCE_CONST(:, :, CELL_IDX_1, P_IDX_1, P_IDX_2) * PHASE
            ENDDO
         ENDDO
      ENDDO

   END FUNCTION PHON_CALC_DYN_MAT_FROM_FC


!                        _______________________________
!_______________________/ PHON_FIND_TRANSLATIONAL_MODES \_______________________
!
!> @brief Finds translational-mode indices of supplied mode matrix.
!>
!> Translational modes only exist at the Gamma point. Therefore, any mode
!> matrix off Gamma does not have any. At Gamma, we look for real eigenvectors
!> that are parallel to each other.
!>
!> @return INDICES   Translational-mode indices in ascending order
!> @param[in] MODES  Phonon-mode matrix
!_______________________________________________________________________________
   FUNCTION PHON_FIND_TRANSLATIONAL_MODES(MODES) RESULT(INDICES)
      IMPLICIT NONE

      INTEGER  :: INDICES(3)

      COMPLEX(q), INTENT(IN)  :: MODES(:, :)

      REAL(q)              :: DOT_PROD
      REAL(q), ALLOCATABLE :: MEASURE(:)
      INTEGER              :: NUM_MODES, MODE_IDX
      INTEGER              :: NUM_ATOMS, ATOM_IDX
      INTEGER              :: CART_IDX

      NUM_MODES = SIZE(MODES, 2)
      NUM_ATOMS = NUM_MODES / 3

      INDICES = 0

      ! If there are only 3 modes, check imaginary part and return all of them
      IF (NUM_ATOMS == 1) THEN
         IF (ALL(ABS(AIMAG(MODES)) > 1E-5_q)) INDICES = [1, 2, 3]
         RETURN
      ENDIF

      ALLOCATE(MEASURE(NUM_MODES))
      MEASURE = 0

      DO MODE_IDX = 1, NUM_MODES
         ! Eigenvectors must be real at Gamma
         IF (ANY(ABS(AIMAG(MODES(:, MODE_IDX))) > 1E-5_q)) CYCLE

         MEASURE(MODE_IDX) = 1

         DO ATOM_IDX = 2, NUM_ATOMS
            DOT_PROD = NUM_ATOMS * DOT_PRODUCT( &
               REAL(MODES(1 : 3, MODE_IDX)), &
               REAL(MODES(ATOM_IDX * 3 - 2 : ATOM_IDX * 3, MODE_IDX)))

            IF (DOT_PROD < 0) DOT_PROD = 0

            MEASURE(MODE_IDX) = MEASURE(MODE_IDX) * DOT_PROD
         ENDDO
      ENDDO

      DO CART_IDX = 1, 3
         MODE_IDX = MAXLOC(MEASURE, 1)
         IF (MEASURE(MODE_IDX) <= 0) EXIT

         INDICES(CART_IDX) = MODE_IDX
         MEASURE(MODE_IDX) = 0
      ENDDO

      ! Sort indices in ascending order
      CALL SWAP_ORDER_INDICES(1, 2)
      CALL SWAP_ORDER_INDICES(1, 3)
      CALL SWAP_ORDER_INDICES(2, 3)

      CONTAINS

      SUBROUTINE SWAP_ORDER_INDICES(I1, I2)
         IMPLICIT NONE
         INTEGER, INTENT(IN)  :: I1, I2
         INTEGER  :: TMP

         IF (INDICES(I1) > INDICES(I2) .AND. INDICES(I2) /= 0) THEN
            TMP = INDICES(I2)
            INDICES(I2) = INDICES(I1)
            INDICES(I1) = TMP
         ENDIF
      END SUBROUTINE
   END FUNCTION PHON_FIND_TRANSLATIONAL_MODES


!                           _________________________
!__________________________/ SETUP_FORCE_CONST_SHORT \__________________________
!
!> @brief Sets up short-ranged force constants for use in Fourier interpolation.
!>
!> Short-ranged force constants are calculated by subtracting the long-ranged
!> dipole-dipole contribution from the ab-initio force constants. They can then
!> be readily used for Fourier interpolation to yield short-ranged dynamical
!> matrices.
!>
!> @param[in] PRIM_CELL Primitive-cell information
!_______________________________________________________________________________
   SUBROUTINE SETUP_FORCE_CONST_SHORT(PRIM_CELL)
      USE constant
      IMPLICIT NONE

      TYPE(PRIM_CELL_T),   INTENT(IN)  :: PRIM_CELL

      REAL(q), ALLOCATABLE :: G_VECTORS(:, :)
      REAL(q)              :: G_CART(3)
      INTEGER              :: G_IDX
      INTEGER              :: P_NUM
      INTEGER              :: P_IDX_1, P_IDX_2
      INTEGER              :: CELL_IDX_1
      INTEGER              :: CART_IDX_1, CART_IDX_2
      REAL(q)              :: PRE_FACTOR
      COMPLEX(q)           :: PHASE, PHASE_PRIM

      P_NUM = PRIM_CELL%NUM_PRIM_ATOMS

      IF (ALLOCATED(FORCE_CONST_SHORT)) DEALLOCATE(FORCE_CONST_SHORT)
      ALLOCATE(FORCE_CONST_SHORT(3, 3, PRIM_CELL%NUM_CELLS, P_NUM, P_NUM))
      FORCE_CONST_SHORT = FORCE_CONST

      CALL GENERATE_EWALD_G_VECTORS(POLAR_EWALD, PRIM_CELL%SUPER_LATT, G_VECTORS)

      ! Loop over all G-vector inside the Ewald cut-off except Gamma!
      DO G_IDX = 2, SIZE(G_VECTORS, 2)
         ! Already in Carthesian coordinates
         G_CART = G_VECTORS(:, G_IDX)

         ! Skip if Ewald suppression factor is too small
         IF (.NOT. POLAR_EWALD%TRY_CALC(G_CART)) CYCLE

         ! Calculate prefactor
         PRE_FACTOR = POLAR_EWALD%DIELEC_FACTOR * POLAR_EWALD%EXP_FACTOR
         PRE_FACTOR = PRE_FACTOR * 4._q * PI / PRIM_CELL%SUPER_LATT%OMEGA * FELECT

         ! Double loop over primitive atoms
         DO P_IDX_2 = 1, P_NUM
            DO P_IDX_1 = 1, P_NUM
               PHASE_PRIM = PRIM_CELL%MIN_IMAGES(1, P_IDX_1, P_IDX_2)%CALC_PRIM_PHASE(G_CART)
               DO CELL_IDX_1 = 1, PRIM_CELL%NUM_CELLS
                  PHASE = PHASE_PRIM * PRIM_CELL%MIN_IMAGES(CELL_IDX_1, P_IDX_1, P_IDX_2)%CALC_CELL_PHASE(G_CART)

                  ! Add/Multiply everything up
                  DO CART_IDX_2 = 1, 3
                     DO CART_IDX_1 = 1, 3
                        FORCE_CONST_SHORT(CART_IDX_1, CART_IDX_2, CELL_IDX_1, P_IDX_1, P_IDX_2) = &
                        FORCE_CONST_SHORT(CART_IDX_1, CART_IDX_2, CELL_IDX_1, P_IDX_1, P_IDX_2) - &
                        POLAR_EWALD%BORN_FACTOR(CART_IDX_1, P_IDX_1) * &
                        POLAR_EWALD%BORN_FACTOR(CART_IDX_2, P_IDX_2) * &
                        PRE_FACTOR * PHASE
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   END SUBROUTINE SETUP_FORCE_CONST_SHORT


!                         ____________________________
!________________________/ PHON_CALC_FC_FROM_DYN_MATS \_________________________
!
!> @brief Calulates supercell force constants from a set of dynamical matrices.
!>
!> Performs a Fourier transform of the dynamical matrices form the primitive
!> cell at each commensurate q-point to the supercell.
!> Currently unused and untested!
!>
!> @return FORCE_CONST  Supercell force constants
!> @param[in] PRIM_CELL Primitive-cell information
!> @param[in] DYN_MATS  Set of dynamical matrices for each q-point
!> @param[in] Q_MESH    Sef of commensurate q-points
!>
!> @todo If this is going to be used, improve it.
!_______________________________________________________________________________
   FUNCTION PHON_CALC_FC_FROM_DYN_MATS(PRIM_CELL, DYN_MATS, Q_MESH) RESULT(FORCE_CONST)
      USE constant
      IMPLICIT NONE

      REAL(q), ALLOCATABLE :: FORCE_CONST(:, :, :, :, :)

      TYPE(PRIM_CELL_T),   INTENT(IN)  :: PRIM_CELL
      COMPLEX(q),          INTENT(IN)  :: DYN_MATS(:, :, :, :, :)
      REAL(q),             INTENT(IN)  :: Q_MESH(:, :)

      COMPLEX(q)  :: FORCE_CONST_SUM(3, 3)
      INTEGER     :: NUM_Q_POINTS
      INTEGER     :: Q_IDX
      INTEGER     :: P_NUM
      INTEGER     :: P_IDX_1, P_IDX_2
      INTEGER     :: CELL_IDX_1
      COMPLEX(q)  :: PHASE
      REAL(q)     :: Q_CART(3)

      NUM_Q_POINTS = SIZE(DYN_MATS, 5)
      P_NUM = PRIM_CELL%NUM_PRIM_ATOMS

      IF (NUM_Q_POINTS /= SIZE(Q_MESH, 2)) THEN
         CALL vtutor%bug("Internal error in PHON_CALC_FC_FROM_DYN_MATS: Numbers of qpoints differ: " &
            // str(NUM_Q_POINTS) // " " // str(SIZE(Q_MESH,2)), __FILE__, __LINE__)
      ENDIF

      IF (ALLOCATED(FORCE_CONST)) DEALLOCATE(FORCE_CONST)

      ALLOCATE(FORCE_CONST(3, 3, PRIM_CELL%NUM_CELLS, P_NUM, P_NUM))

      ! Double loop over primitive atoms.
      DO P_IDX_2 = 1, P_NUM
         DO P_IDX_1 = 1, P_NUM
            ! Loop over cells.
            DO CELL_IDX_1 = 1, PRIM_CELL%NUM_CELLS
               FORCE_CONST_SUM = 0

               ! Loop over all supplied q-points.
               DO Q_IDX = 1, NUM_Q_POINTS
                  ! Convert to Carthesian coordinates.
                  Q_CART = MATMUL(PRIM_CELL%PRIM_LATT%B, Q_MESH(:, Q_IDX))

                  PHASE = PRIM_CELL%MIN_IMAGES(CELL_IDX_1, P_IDX_1, P_IDX_2)%CALC_CELL_PHASE(Q_CART)

                  ! Interpolate each 3x3 dynamical matrix to real space.
                  FORCE_CONST_SUM(:, :) = FORCE_CONST_SUM(:, :) + &
                  DYN_MATS(:, :, P_IDX_1, P_IDX_2, Q_IDX) * PHASE
               ENDDO

               ! If everything is ok, the force constants are real at this point.
               FORCE_CONST(:, :, CELL_IDX_1, P_IDX_1, P_IDX_2) = REAL(FORCE_CONST_SUM)
            ENDDO
         ENDDO
      ENDDO

      FORCE_CONST = FORCE_CONST / NUM_Q_POINTS

   END FUNCTION PHON_CALC_FC_FROM_DYN_MATS


!                            ______________________
!___________________________/ PHON_CALC_DISPERSION \____________________________
!
!> @brief Writes phonon frequencies along a path in reciprocal space.
!>
!> The set of q-points belonging to the path is generated from the QPOINTS file
!> (same format as KPOINTS).
!>
!> @param[in,out] PRIM_CELL   Primitive-cell information
!> @param[in] T_INFO          Atom-type information
!> @param[in] SETTINGS        Phonon settings
!> @param[in] IO              IO-unit information
!_______________________________________________________________________________
   SUBROUTINE PHON_CALC_DISPERSION(PRIM_CELL, T_INFO, SETTINGS, IO)
      USE poscar_struct_def,     ONLY: type_info
      USE mkpoints_struct_def,   ONLY: kpoints_struct
      USE mkpoints,              ONLY: RD_KPOINTS
      USE constant
      IMPLICIT NONE

      TYPE(PRIM_CELL_T),      INTENT(INOUT)  :: PRIM_CELL
      TYPE(type_info),        INTENT(IN)     :: T_INFO
      TYPE(PHONON_SETTINGS),  INTENT(IN)     :: SETTINGS
      TYPE(in_struct),        INTENT(IN)     :: IO

      ! locals
      INTEGER                 :: NUM_BRANCHES
      INTEGER                 :: Q_IDX
      TYPE(kpoints_struct)    :: Q_MESH
      REAL(q),    ALLOCATABLE :: FREQU(:, :)
      REAL(q),    ALLOCATABLE :: QDIRECTIONS(:,:)
      REAL(q)                 :: QDIRECTION(3), NEW_QDIRECTION(3)
      COMPLEX(q), ALLOCATABLE :: MODES(:, :, :)
      ! remove after refactor
      REAL(q),    ALLOCATABLE :: FREQU_TMP(:)
      COMPLEX(q), ALLOCATABLE :: MODES_TMP(:, :)

      NUM_BRANCHES = PRIM_CELL%NUM_PRIM_ATOMS * 3

      ! Setup the minimum-distance vectors in the primitive cell if needed.
      CALL PRIM_CELL%SETUP_ADDITIONAL_DATA()

      ! Read q-mesh from file for which the dispersion should be evaluated.
      Q_MESH%LTET = .FALSE. ! _test_mani fix
      CALL RD_KPOINTS(Q_MESH, PRIM_CELL%PRIM_LATT, .FALSE., .FALSE., IO%IU6, IO%IU0, FILE_NAME = 'QPOINTS')

      ALLOCATE(FREQU(NUM_BRANCHES, Q_MESH%NKPTS), MODES(NUM_BRANCHES, NUM_BRANCHES, Q_MESH%NKPTS))
      ALLOCATE(QDIRECTIONS(3,Q_MESH%NKPTS))

      ! Estimate q-directions
      ! TODO: propagate segments from QPOINTS file in the k-points structure and use them as directions
      QDIRECTION=[1.0_q,0.0_q,0.0_q]
      DO Q_IDX = 1, Q_MESH%NKPTS-1
         NEW_QDIRECTION = Q_MESH%VKPT(:, Q_IDX+1)-Q_MESH%VKPT(:, Q_IDX)
         IF (SQRT(SUM(NEW_QDIRECTION**2))>1e-6) QDIRECTION=NEW_QDIRECTION
         QDIRECTIONS(:,Q_IDX) = QDIRECTION
      ENDDO
      QDIRECTIONS(:,Q_MESH%NKPTS) = QDIRECTION

      ! _test_mani parallelize
      DO Q_IDX = 1, Q_MESH%NKPTS
         ! Calculate phonons for given q-point.
         ! _test_mani TODO: Refactor phonon routine to write into arrays instead of allocating them
         CALL PHON_CALC_PHONONS(Q_MESH%VKPT(:, Q_IDX), PRIM_CELL, T_INFO, FREQU_TMP, MODES_TMP, QDIRECTIONS(:,Q_IDX))
         FREQU(:, Q_IDX) = FREQU_TMP
         MODES(:, :, Q_IDX) = MODES_TMP
      ENDDO

      CALL PHON_WRITE_DATA(PRIM_CELL, Q_MESH, FREQU, MODES, SETTINGS%NWRITE, IO)

      ! CALL PHON_SHOW_FORCE_CONST_DECAY(T_INFO, LATT_CUR, IO)

      ! Deallocate KPOINTS structure before it goes out of scope.
      DEALLOCATE(Q_MESH%VKPT, Q_MESH%WTKPT, Q_MESH%IDTET)
   END SUBROUTINE PHON_CALC_DISPERSION


!                               _________________
!______________________________/ PHON_WRITE_DATA \______________________________
!
!> @brief Writes phonon frequencies and eigenvectors on a q-point mesh to disk.
!>
!> @param[in] PRIM_CELL Primitive-cell information
!> @param[in] Q_MESH    Mesh or path of q-points
!> @param[in] FREQU     Frequencies
!> @param[in] MODES     Eigenvectors
!> @param[in] NWRITE    Determines how much output is written to OUTCAR
!> @param[in] IO        IO-unit information
!_______________________________________________________________________________
   SUBROUTINE PHON_WRITE_DATA(PRIM_CELL, Q_MESH, FREQU, MODES, NWRITE, IO)
      USE mkpoints_struct_def,   ONLY: kpoints_struct
      USE finite_differences,    ONLY: CONVERT_FREQUENCY, &
                                       UNIT_INTERNAL, &
                                       UNIT_EV, &
                                       UNIT_HZ, &
                                       UNIT_2PI_HZ, &
                                       UNIT_INV_CM
#ifdef VASP_HDF5
      USE vhdf5,  ONLY: vh5_write_phonons, ih5outfileid
#endif
      USE constant
      IMPLICIT NONE

      TYPE(PRIM_CELL_T),      INTENT(IN)     :: PRIM_CELL
      TYPE(kpoints_struct),   INTENT(IN)     :: Q_MESH
      REAL(q),                INTENT(INOUT)  :: FREQU(:, :)
      COMPLEX(q),             INTENT(IN)     :: MODES(:, :, :)
      INTEGER,                INTENT(IN)     :: NWRITE
      TYPE(in_struct),        INTENT(IN)     :: IO

      INTEGER  :: Q_IDX
      INTEGER  :: I, J, IATOM
      REAL(q)  :: Q_TOT
      REAL(q)  :: Q_CART(3)
      REAL(q)  :: Q_DIFF, Q_DIFF_OLD
      INTEGER  :: NUM_MODES, NUM_ATOMS
      INTEGER  :: IU6

      CHARACTER(2), PARAMETER :: C_DELIM = '  '
      REAL(q),PARAMETER       :: PLANCK_SI = 6.626070040E-34_q
      REAL(q),PARAMETER       :: C_SI = 299792458_q

      IU6 = IO%IU6
      Q_TOT = 0
      NUM_MODES = SIZE(FREQU, 1)
      NUM_ATOMS = NUM_MODES/3

      IF (Q_MESH%NKPTS > 1) THEN
         Q_CART = MATMUL(PRIM_CELL%PRIM_LATT%B, Q_MESH%VKPT(:, 2) - Q_MESH%VKPT(:, 1))
         Q_DIFF_OLD = SQRT(DOT_PRODUCT(Q_CART, Q_CART))
      ENDIF

      IF (IO%IU0 >= 0) THEN
         ! Write header
         WRITE(IU6, '(/, A, //, A, /)') &
            '--------------------------------------------------------------------------------', &
            '                                    Phonons'

         ! Write phonon output to OUTCAR?
         IF (ABS(NWRITE) > 0) THEN
            ! Write header for compact list?
            IF (NWRITE == -1) WRITE(IU6, '(A)') &
               ' q-point No.' // C_DELIM // '         q_x' // C_DELIM // '         q_y' // C_DELIM // &
               '         q_z' // C_DELIM // 'q-path length' // C_DELIM // 'frequencies [THz]'

            ! Iterate over q-points
            DO Q_IDX = 1, Q_MESH%NKPTS
               ! Get total q-distance traveled.
               IF (Q_IDX > 1) THEN
                  Q_CART = MATMUL(PRIM_CELL%PRIM_LATT%B, Q_MESH%VKPT(:, Q_IDX) - Q_MESH%VKPT(:, Q_IDX - 1))
                  Q_DIFF = SQRT(DOT_PRODUCT(Q_CART, Q_CART))

                  IF (ABS(Q_DIFF) < 1E-7_q .OR. ABS(Q_DIFF - Q_DIFF_OLD) > 1E-7_q) THEN
                     IF (Q_IDX < Q_MESH%NKPTS) THEN
                        Q_CART = MATMUL(PRIM_CELL%PRIM_LATT%B, Q_MESH%VKPT(:, Q_IDX + 1) - Q_MESH%VKPT(:, Q_IDX))
                        Q_DIFF_OLD = SQRT(DOT_PRODUCT(Q_CART, Q_CART))
                     ENDIF
                  ELSE
                     Q_TOT = Q_TOT + Q_DIFF
                  ENDIF
               ENDIF

               ! Write single line for compact list?
               IF (NWRITE == -1) THEN
                  WRITE(IU6, '(A, I6, A, F12.6, A, F12.6, A, F12.6, A, F13.6, *(F14.6))') &
                     '      ', Q_IDX, C_DELIM, &
                     Q_MESH%VKPT(1, Q_IDX), C_DELIM, &
                     Q_MESH%VKPT(2, Q_IDX), C_DELIM, &
                     Q_MESH%VKPT(3, Q_IDX), C_DELIM, &
                     Q_TOT, CONVERT_FREQUENCY(FREQU(:, Q_IDX), UNIT_INTERNAL, UNIT_HZ) * 1E-12
               ENDIF

               ! Write qpoint info for human readable output
               IF (NWRITE > 0) THEN
                  WRITE(IU6, '(A)') &
                     ' q-point No.' // C_DELIM // '         q_x' // C_DELIM // '         q_y' // C_DELIM // &
                     '         q_z' // C_DELIM // 'q-path length'

                  WRITE(IU6, '(A, I6, A, F12.6, A, F12.6, A, F12.6, A, F13.6, *(F14.6))') &
                     '      ', Q_IDX, C_DELIM, &
                     Q_MESH%VKPT(1, Q_IDX), C_DELIM, &
                     Q_MESH%VKPT(2, Q_IDX), C_DELIM, &
                     Q_MESH%VKPT(3, Q_IDX), C_DELIM, &
                     Q_TOT
               ENDIF

               ! Write q-point block?
               IF (NWRITE <= -2) THEN
                  WRITE(IU6, '(A)') &
                     ' q-point No.' // C_DELIM // '         q_x' // C_DELIM // &
                     '         q_y' // C_DELIM // '         q_z' // C_DELIM // 'q-path length'
                  WRITE(IU6, '(A, I6, A, F12.6, A, F12.6, A, F12.6, A, F13.6)') &
                     '      ', Q_IDX, C_DELIM, &
                     Q_MESH%VKPT(1, Q_IDX), C_DELIM, &
                     Q_MESH%VKPT(2, Q_IDX), C_DELIM, &
                     Q_MESH%VKPT(3, Q_IDX), C_DELIM, &
                     Q_TOT
               ENDIF

               ! Write frequency block?
               IF (NWRITE <= -2) THEN
                  WRITE(IU6, '(A)') '  frequencies'
                  WRITE(IU6, '(A, *(F14.6))')   '   [THz]     ', CONVERT_FREQUENCY(FREQU(:, Q_IDX), UNIT_INTERNAL, UNIT_HZ) * 1E-12
                  WRITE(IU6, '(A, *(F14.6))')   '   [2Pi THz] ', CONVERT_FREQUENCY(FREQU(:, Q_IDX), UNIT_INTERNAL, UNIT_2PI_HZ) * 1E-12
                  WRITE(IU6, '(A, *(F14.6))')   '   [cm-1]    ', CONVERT_FREQUENCY(FREQU(:, Q_IDX), UNIT_INTERNAL, UNIT_INV_CM)
                  WRITE(IU6, '(A, *(F14.6))')   '   [meV]     ', CONVERT_FREQUENCY(FREQU(:, Q_IDX), UNIT_INTERNAL, UNIT_EV) * 1E3
               ENDIF

               ! Write eigenvector block?
               IF (NWRITE <= -3) THEN
                  WRITE(IU6, '(A)') '  eigenvectors'
                  DO I = 1, NUM_MODES
                     WRITE(IU6, '(A)', ADVANCE = 'no') '   '
                     DO J = 1, NUM_MODES
                        WRITE(IU6, '(F12.6, SP, F12.6, "i", A)', ADVANCE = 'no') MODES(I, J, Q_IDX), C_DELIM
                     ENDDO
                     WRITE(IU6, *)
                  ENDDO
               ENDIF

               ! Human readable output
               IF (NWRITE == 1) THEN
                  WRITE(IU6, '(A)')    ' frequencies'
                  WRITE(IU6, '(*(A))') 'branch index','         f[THz]', '      f[2PiTHz]', '        f[cm-1]','         f[meV]'
               ENDIF
               DO I=1,NUM_MODES

                  ! write phonon frequencies block
                  IF (NWRITE == 2) THEN
                     WRITE(IU6, '(A)')    ' frequencies'
                     WRITE(IU6, '(*(A))') 'branch index','         f[THz]', '      f[2PiTHz]', '        f[cm-1]','         f[meV]'
                  ENDIF

                  IF(NWRITE > 0) THEN
                     WRITE(IU6, '(I12)', ADVANCE = 'NO')   I
                     WRITE(IU6, '(*(F15.6))') FREQU(I,Q_IDX) * SQRT(EVTOJ / AMTOKG) * 1E10 / (1E12 * 2 * PI), &                   ! THz
                                              FREQU(I,Q_IDX) * SQRT(EVTOJ / AMTOKG) * 1E10 * 1E-12, &                             ! 2Pi THz
                                              FREQU(I,Q_IDX) * SQRT(EVTOJ / AMTOKG) * 1E10 / (1E2 * C_SI * 2 * PI), &             ! cm^(-1)
                                              FREQU(I,Q_IDX) * SQRT(EVTOJ / AMTOKG) * 1E10 * 1E3 * PLANCK_SI / (EVTOJ * 2 * PI)   ! meV
                  ENDIF

                  ! write phonon modes
                  IF (NWRITE == 2) THEN
                     WRITE(IU6, '(A)')    '       modes'
                     WRITE(IU6, '(*(A))') '  atom index','                        x',C_DELIM, &
                                                         '                        y',C_DELIM, &
                                                         '                        z'
                     DO IATOM=1, NUM_ATOMS
                        WRITE(IU6,'(I12)', ADVANCE = 'no') IATOM
                        DO J=1,3
                           WRITE(IU6, '(F12.6, SP, F12.6, "i", A)', ADVANCE = 'no') MODES(I, (IATOM-1)*3+J, Q_IDX), C_DELIM
                        ENDDO
                        WRITE(IU6, *)
                     ENDDO
                     WRITE(IU6, *)
                  ENDIF
               ENDDO

               ! Add a newline between q-point blocks
               IF (NWRITE <= -2 .OR. NWRITE > 0) WRITE(IU6, *)
            ENDDO
         ELSE
            ! Write short summary
            WRITE(IU6, *) "   Summary:"
            WRITE(IU6, *) "     ", Q_MESH%NKPTS, "q-points"
            WRITE(IU6, *) "     ", NUM_MODES, "frequencies per q-point"
         ENDIF

         ! Write footer
         WRITE(IU6, '(/, A, /)') '--------------------------------------------------------------------------------'
      ENDIF
      
      ! Write HDF5 output
#ifdef VASP_HDF5
      FREQU = CONVERT_FREQUENCY(FREQU, UNIT_INTERNAL, UNIT_HZ) * 1E-12
      CALL vh5_write_phonons(ih5outfileid, Q_MESH, FREQU, MODES)
      FREQU = CONVERT_FREQUENCY(FREQU, UNIT_HZ, UNIT_INTERNAL) * 1E12
#endif
   END SUBROUTINE PHON_WRITE_DATA


#ifdef VASP_HDF5
      SUBROUTINE PHON_WRITE_HDF5(LOC_ID, GROUP_PATH)
         USE hdf5
         USE vhdf5_base,   ONLY: vh5_error, vh5_write
         IMPLICIT NONE
   
         INTEGER(HID_T),   INTENT(IN)  :: LOC_ID
         CHARACTER(*),     INTENT(IN)  :: GROUP_PATH

         IF (.NOT. ALLOCATED(FORCE_CONST)) CALL vtutor%bug("Force constants not allocated!", __FILE__, __LINE__)
         VH5_CHECK(vh5_write(LOC_ID, GROUP_PATH // '/force_constants', FORCE_CONST))
      END SUBROUTINE PHON_WRITE_HDF5

      SUBROUTINE PHON_READ_HDF5(PRIM_CELL, LOC_ID, GROUP_PATH)
         USE hdf5
         USE vhdf5_base,   ONLY  : vh5_error &
                                 , vh5_exists &
                                 , vh5_read
         IMPLICIT NONE
   
         TYPE(PRIM_CELL_T),   INTENT(IN)  :: PRIM_CELL
         INTEGER(HID_T),      INTENT(IN)  :: LOC_ID
         CHARACTER(*),        INTENT(IN)  :: GROUP_PATH
   
         INTEGER  :: SPIN_IDX, NUM_SPIN

         IF (ALLOCATED(FORCE_CONST)) DEALLOCATE(FORCE_CONST)
         ALLOCATE(FORCE_CONST(3, 3, PRIM_CELL%NUM_CELLS, PRIM_CELL%NUM_PRIM_ATOMS, PRIM_CELL%NUM_PRIM_ATOMS))
         VH5_CHECK(vh5_read(LOC_ID, GROUP_PATH // '/force_constants', FORCE_CONST))

         CALL CORRECT_FORCE_CONSTANTS()
      END SUBROUTINE PHON_READ_HDF5
#endif


!                         _____________________________
!________________________/ PHON_SHOW_FORCE_CONST_DECAY \________________________
!
!> @brief Writes spatial decay of force-constant matrix to file.
!>
!> This is used for debugging purposes and can help determine if the supercell
!> is large enough for the Fourier interpolation to yield meaningful results.
!>
!> @param[in] PRIM_CELL Primitive-cell information
!> @param[in] IO        IO-unit information
!_______________________________________________________________________________
   SUBROUTINE PHON_SHOW_FORCE_CONST_DECAY(PRIM_CELL, IO)
      IMPLICIT NONE

      TYPE(PRIM_CELL_T),   INTENT(IN)  :: PRIM_CELL
      TYPE(in_struct),     INTENT(IN)  :: IO

      ! locals
      INTEGER              :: PRIM_ION
      INTEGER              :: S_NUM
      INTEGER              :: A1, A2, A3
      INTEGER              :: I, J
      REAL(q)              :: DIST, CUR_DIST
      REAL(q)              :: DIST_VEC(3)
      REAL(q), ALLOCATABLE :: DIST_ARR(:)
      REAL(q), ALLOCATABLE :: PLOT_ARR(:, :)
      INTEGER              :: DATA_LEN

      ! Generate arrays to hold nicely ordered distances and values.
      ALLOCATE(DIST_ARR(S_NUM), PLOT_ARR(S_NUM, 2))

      ! Set reference atom.
      PRIM_ION = PRIM_CELL%SUPER_IDX(1, 1)
      DIST_ARR(PRIM_ION) = 0._q

      ! Loop over all ions.
      DO I = 1, PRIM_CELL%NUM_SUPER_ATOMS
         ! This case was already handled.
         IF (I == PRIM_ION) CYCLE

         ! Get minimum-distance image.
         DIST = HUGE(0._q)
         ! _test_mani
         DO A1 = -2, 2
         DO A2 = -2, 2
         DO A3 = -2, 2
            DIST_VEC = PRIM_CELL%POSION_SUPER(:, I) - PRIM_CELL%POSION_SUPER(:, PRIM_ION) + [A1, A2, A3]
            DIST_VEC = MATMUL(PRIM_CELL%SUPER_LATT%A, DIST_VEC)
            CUR_DIST = DOT_PRODUCT(DIST_VEC, DIST_VEC)

            IF (CUR_DIST < DIST) DIST = CUR_DIST
         ENDDO
         ENDDO
         ENDDO

         ! Store the minimum distance in the array.
         DIST_ARR(I) = SQRT(DIST)
      ENDDO

      ! Write out spatial decay of the force constants.
      IF (IO%IU0 >= 0) THEN
         OPEN(UNIT = 1337, FILE = 'debug_fc_decay', STATUS = 'replace')

         ! First column is just distance.
         PLOT_ARR(:, 1) = DIST_ARR(:)
         ! Second column contains maximum values.
         DO I = 1, PRIM_CELL%NUM_SUPER_ATOMS
            PLOT_ARR(I, 2) = MAXVAL(ABS(FORCE_CONST(:, :, PRIM_CELL%CELL_IDX(I), PRIM_CELL%PRIM_IDX(I), 1)))
         ENDDO
         ! Sort the array based on distances.
         CALL SORT_PLOT_ARRAY(PLOT_ARR, DATA_LEN)
         ! Normalize to 1.
         ! PLOT_ARR(:, 2) = PLOT_ARR(:, 2) / MAXVAL(PLOT_ARR(:, 2))
         ! Write to file.
         WRITE(1337, '(F12.6, ", ", F12.6)') (PLOT_ARR(J, 1:2), J = 1, DATA_LEN)

         WRITE(1337, *); WRITE(1337, *)
         ! First column is just distance.
         PLOT_ARR(:, 1) = DIST_ARR(:)
         ! Second column contains maximum values.
         DO I = 1, PRIM_CELL%NUM_SUPER_ATOMS
            PLOT_ARR(I, 2) = MAXVAL(ABS(FORCE_CONST_SHORT(:, :, PRIM_CELL%CELL_IDX(I), PRIM_CELL%PRIM_IDX(I), 1)))
         ENDDO
         ! Sort the array based on distances.
         CALL SORT_PLOT_ARRAY(PLOT_ARR, DATA_LEN)
         ! Normalize to 1.
         ! PLOT_ARR(:, 2) = PLOT_ARR(:, 2) / MAXVAL(PLOT_ARR(:, 2))
         ! Write to file.
         WRITE(1337, '(F12.6, ", ", F12.6)') (PLOT_ARR(J, 1:2), J = 1, DATA_LEN)

         CLOSE(1337)
      ENDIF

   END SUBROUTINE PHON_SHOW_FORCE_CONST_DECAY

ENDMODULE phonon
