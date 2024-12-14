#include "symbol.inc"
MODULE msymmetry
  USE prec
  USE lattice
! TODO place COMMON block /SYMM/ in this module
! and include msymmetry everywhere where /SYMM/ is used
!

!************************* TYPE MIN_IMAGE_T ****************************
!
! Facilitates minimum-image convention between two entities (e.g.
! atoms, Wannier functions) in the periodic supercell.
! It is used in interpolation codes (phonons, electron-phonon) to
! calculate Fourier sums over primitive cells.
!
! The IMAGES array holds the inter-cell distance vectors that point
! to all primitive cells containing equivalent minimum images.
! PRIM_DIST contains the vector between the two entities in the
! primitive cell.
! WEIGHT is used in the evaluation of the Fourier phases. At the
! moment, it is simply 1 / IMAGE_COUNT.
!
!***********************************************************************
   TYPE MIN_IMAGE_T
      INTEGER, PRIVATE              :: IMAGE_COUNT = 0
      REAL(q), PRIVATE, ALLOCATABLE :: IMAGES(:, :)
      REAL(q), PRIVATE              :: PRIM_DIST(3)
      REAL(q), PRIVATE              :: WEIGHT

      CONTAINS
      PROCEDURE, PUBLIC :: COUNT => GET_MIN_IMAGE_COUNT
      PROCEDURE, PUBLIC :: GET_IMAGE => GET_MIN_IMAGE
      PROCEDURE, PUBLIC :: CALC_CELL_PHASE => CALC_MIN_IMAGE_CELL_PHASE
      PROCEDURE, PUBLIC :: CALC_PRIM_PHASE => CALC_MIN_IMAGE_PRIM_PHASE
      PROCEDURE, PUBLIC :: SERIALIZE => SERIALIZE_MIN_IMAGE
      PROCEDURE, PUBLIC :: DESERIALIZE => DESERIALIZE_MIN_IMAGE
   ENDTYPE MIN_IMAGE_T

   ! Overload type name with constructor
   INTERFACE MIN_IMAGE_T
      PROCEDURE   :: MIN_IMAGE_INIT
   END INTERFACE


!************************* TYPE PRIM_CELL_T ****************************
!
! Contains important information about the primitive cell and how it
! relates to the supercell.
!
! NUM_PRIM_ATOMS        - # of atoms in primitive cell.
! NUM_SUPER_ATOMS       - # of atoms in supercell.
! NUM_CELLS             - # of primitive cell.
! SUPER_LATT            - Lattice of supercell.
! PRIM_LATT             - Lattice of primitive cell.
! POSION_SUPER(:, i)    - Direct coordinates of supercell atom i w.r.t.
!                         supercell lattice.
! SUPER_IDX(l, k)       - Supercell atom index (= POSCAR index) of
!                         primitive atom k in cell l.
! PRIM_IDX(i)           - Primitive-atom index of supercell atom i.
! CELL_IDX(i)           - Cell index of supercell atom i.
! EQUIV_CELLS(l, k)     - Map from cell indices l and k to cell index m
!                         such that (R_l - R_k) = (R_m - R_0) within
!                         the periodic supercell.
! MIN_IMAGES(l, k, k')  - Minimum images between primitive atom k' in
!                         cell 1 and primitive atom k in cell l.
! PRIM_TO_SUPER(3, 3)   - Transformation matrix between direct
!                         coordinates from primitive to supercell.
! SUPER_TO_PRIM(3, 3)   - Inverse transformation of the above.
!
!***********************************************************************
   TYPE PRIM_CELL_T
      INTEGER                          :: NUM_PRIM_ATOMS
      INTEGER                          :: NUM_SUPER_ATOMS
      INTEGER                          :: NUM_CELLS
      TYPE(LATT)                       :: SUPER_LATT
      TYPE(LATT)                       :: PRIM_LATT
      REAL(q), ALLOCATABLE             :: POSION_SUPER(:, :)
      INTEGER, ALLOCATABLE             :: SUPER_IDX(:, :)
      INTEGER, ALLOCATABLE             :: PRIM_IDX(:)
      INTEGER, ALLOCATABLE             :: CELL_IDX(:)
      INTEGER, ALLOCATABLE             :: EQUIV_CELLS(:, :)
      TYPE(MIN_IMAGE_T), ALLOCATABLE   :: MIN_IMAGES(:, :, :)
      REAL(q)                          :: PRIM_TO_SUPER(3, 3)
      REAL(q)                          :: SUPER_TO_PRIM(3, 3)

      CONTAINS
      PROCEDURE,  PUBLIC   :: GET_POSION => PRIM_CELL_GET_POSION
      PROCEDURE,  PUBLIC   :: WRITE_POSCAR => PRIM_CELL_WRITE_POSCAR
      PROCEDURE,  PUBLIC   :: WRITE_OUTCAR => PRIM_CELL_WRITE_OUTCAR
      PROCEDURE,  PUBLIC   :: WRITE_XML => PRIM_CELL_WRITE_XML
      PROCEDURE,  PUBLIC   :: SERIALIZE => PRIM_CELL_SERIALIZE
      PROCEDURE,  PUBLIC   :: DESERIALIZE => PRIM_CELL_DESERIALIZE
      PROCEDURE,  PRIVATE  :: SETUP_SUPER_INDICES => PRIM_CELL_SETUP_SUPER_INDICES
      PROCEDURE,  PRIVATE  :: SETUP_PRIM_INDICES => PRIM_CELL_SETUP_PRIM_INDICES

      PROCEDURE,  PUBLIC   :: SETUP_ADDITIONAL_DATA => PRIM_CELL_SETUP_ADDITIONAL_DATA
      PROCEDURE,  PRIVATE  :: CALC_MIN_IMAGES => PRIM_CELL_CALC_MIN_IMAGES
      PROCEDURE,  PRIVATE  :: SETUP_EQUIV_CELLS => PRIM_CELL_SETUP_EQUIV_CELLS
#ifdef VASP_HDF5
      PROCEDURE,  PUBLIC   :: WRITE_HDF5 => PRIM_CELL_WRITE_HDF5
      PROCEDURE,  PUBLIC   :: READ_HDF5 => PRIM_CELL_READ_HDF5
#endif
   END TYPE PRIM_CELL_T

   ! Overload type name with constructor
   INTERFACE PRIM_CELL_T
      MODULE PROCEDURE   :: PRIM_CELL_INIT
   END INTERFACE


   CONTAINS

!********************** FUNCTION MIN_IMAGE_INIT ************************
!
! Constructor of MIN_IMAGE_T.
!
! PRIM_DISTANCE   - Distance between the two entities in direct
!                   coordinates w.r.t. the supercell lattice stored in
!                   PRIM_CELL.
! CELL_IDX        - Cell index for the second entity (entity 1 is
!                   always located in cell 1).
! PRIM_CELL       - Primitive-cell structure that holds cell info.
! TOLERANCE       - Optional; Max distance in Angstrom up to which two
!                   images are considered equivalent.
!                   Default = TINY from sym_prec
!
!***********************************************************************
   TYPE(MIN_IMAGE_T) FUNCTION MIN_IMAGE_INIT(PRIM_DISTANCE, CELL_IDX, PRIM_CELL, TOLERANCE) RESULT(MIN_IMAGE)
      USE LATTICE,   ONLY  : LATT
      USE sym_prec,  ONLY  : TINY
      IMPLICIT NONE
      REAL(q),             INTENT(IN)  :: PRIM_DISTANCE(3)
      INTEGER,             INTENT(IN)  :: CELL_IDX
      TYPE(PRIM_CELL_T),   INTENT(IN)  :: PRIM_CELL
      REAL(q), OPTIONAL,   INTENT(IN)  :: TOLERANCE

      INTEGER  :: X, Y, Z
      REAL(q)  :: DIST_VEC(3), CELL_DIST_VEC(3), SHIFT_VEC(3)
      REAL(q)  :: VEC_BUFFER(3, 8)
      INTEGER  :: NUM_IMAGES
      REAL(q)  :: MY_TOLERANCE
      REAL(q)  :: DIST, MIN_DIST

      MY_TOLERANCE = TINY
      IF (PRESENT(TOLERANCE)) MY_TOLERANCE = TOLERANCE

      MIN_IMAGE%PRIM_DIST = MODULO(PRIM_DISTANCE + .5_q, 1._q) - .5_q
      MIN_IMAGE%PRIM_DIST = MATMUL(PRIM_CELL%SUPER_LATT%A, PRIM_DISTANCE)

      ! Inter-cell distance vector in supercell
      CELL_DIST_VEC = PRIM_CELL%GET_POSION(CELL_IDX, 1) - PRIM_CELL%GET_POSION(1, 1)
      CELL_DIST_VEC = MODULO(CELL_DIST_VEC + .5_q, 1._q) - .5_q
      CELL_DIST_VEC = MATMUL(PRIM_CELL%SUPER_LATT%A, CELL_DIST_VEC)

      MIN_DIST = HUGE(1._q)

      ! Loop over all neighboring supercells
      DO Z = -1, 1
         DO Y = -1, 1
            DO X = -1, 1
               ! Supercell lattice shift in Cartesian coordinates
               SHIFT_VEC = MATMUL(PRIM_CELL%SUPER_LATT%A, [REAL(q):: X, Y, Z])
               ! Cartesian distance vector between the images
               DIST_VEC = CELL_DIST_VEC + MIN_IMAGE%PRIM_DIST + SHIFT_VEC
               DIST = NORM2(DIST_VEC)

               IF (ABS(DIST - MIN_DIST) < MY_TOLERANCE) THEN ! We already got this one!
                  NUM_IMAGES = NUM_IMAGES + 1
                  VEC_BUFFER(:, NUM_IMAGES) = CELL_DIST_VEC + SHIFT_VEC ! Only store inter-cell distance
                  MIN_DIST = DIST
               ELSE IF (DIST < MIN_DIST) THEN ! We found a new one!
                  NUM_IMAGES = 1
                  VEC_BUFFER(:, 1) = CELL_DIST_VEC + SHIFT_VEC ! Only store inter-cell distance
                  MIN_DIST = DIST
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      ! Allocate just the necessary space
      ALLOCATE(MIN_IMAGE%IMAGES(3, NUM_IMAGES))
      MIN_IMAGE%IMAGE_COUNT = NUM_IMAGES
      MIN_IMAGE%IMAGES = VEC_BUFFER(:, 1 : NUM_IMAGES)
      MIN_IMAGE%WEIGHT = 1._q / NUM_IMAGES
   END FUNCTION

!******************* FUNCTION GET_MIN_IMAGE_COUNT **********************
!
! Getter for number of images.
!
!***********************************************************************
   PURE INTEGER FUNCTION GET_MIN_IMAGE_COUNT(THIS) RESULT(COUNT)
      IMPLICIT NONE
      CLASS(MIN_IMAGE_T),  INTENT(IN)  :: THIS
      COUNT = THIS%IMAGE_COUNT
   END FUNCTION

!********************** FUNCTION GET_MIN_IMAGE *************************
!
! Getter for the inter-cell image vector with index INDEX.
!
!***********************************************************************
   PURE FUNCTION GET_MIN_IMAGE(THIS, INDEX) RESULT(VEC)
      IMPLICIT NONE
      REAL(q)                          :: VEC(3)
      CLASS(MIN_IMAGE_T),  INTENT(IN)  :: THIS
      INTEGER,             INTENT(IN)  :: INDEX
      VEC = THIS%IMAGES(:, INDEX)
   END FUNCTION

!**************** FUNCTION CALC_MIN_IMAGE_CELL_PHASE *******************
!
! Calculates the Fourier phase exp(2 pi i R . k) where R is an inter-
! cell minimum-image vector in the periodic supercell. Equivalent
! minimum-image vectors are handled correctly.
! 
! K_CART - The reciprocal-space vector k in Cartesian coordinates.
!
!***********************************************************************
   PURE COMPLEX(q) FUNCTION CALC_MIN_IMAGE_CELL_PHASE(THIS, K_CART) RESULT(PHASE)
      USE constant
      IMPLICIT NONE
      CLASS(MIN_IMAGE_T),  INTENT(IN)  :: THIS
      REAL(q),             INTENT(IN)  :: K_CART(3)

      INTEGER  :: IMAGE_IDX

      PHASE = 0
      DO IMAGE_IDX = 1, THIS%IMAGE_COUNT
         PHASE = PHASE + EXP(CITPI * DOT_PRODUCT(K_CART, THIS%IMAGES(:, IMAGE_IDX)))
      ENDDO
      PHASE = PHASE * THIS%WEIGHT
   END FUNCTION

!**************** FUNCTION CALC_MIN_IMAGE_PRIM_PHASE *******************
!
! Calculates the Fourier phase exp(2 pi i r . k) where r is the
! distance vector between the two entities in the primitive cell.
! 
! K_CART - The reciprocal-space vector k in Cartesian coordinates.
!
!***********************************************************************
   PURE COMPLEX(q) FUNCTION CALC_MIN_IMAGE_PRIM_PHASE(THIS, K_CART) RESULT(PHASE)
      USE constant
      IMPLICIT NONE
      CLASS(MIN_IMAGE_T),  INTENT(IN)  :: THIS
      REAL(q),             INTENT(IN)  :: K_CART(3)
      PHASE = EXP(CITPI * DOT_PRODUCT(K_CART, THIS%PRIM_DIST))
   END FUNCTION

! _test_mani serialization routine to be replaced with hdf5
   SUBROUTINE SERIALIZE_MIN_IMAGE(THIS, IO, OUT_UNIT)
      USE base,   ONLY  : in_struct
      IMPLICIT NONE
      CLASS(MIN_IMAGE_T),  INTENT(IN)  :: THIS
      TYPE(in_struct),     INTENT(IN)  :: IO
      INTEGER,             INTENT(IN)  :: OUT_UNIT

      IF (IO%IU0 < 0) RETURN

      WRITE(OUT_UNIT) THIS%IMAGE_COUNT
      WRITE(OUT_UNIT) THIS%IMAGES
      WRITE(OUT_UNIT) THIS%WEIGHT
   END SUBROUTINE

! _test_mani deserialization routine to be replaced with hdf5
   SUBROUTINE DESERIALIZE_MIN_IMAGE(THIS, IN_UNIT)
      IMPLICIT NONE
      CLASS(MIN_IMAGE_T),  INTENT(INOUT)  :: THIS
      INTEGER,             INTENT(IN)     :: IN_UNIT

      IF (ALLOCATED(THIS%IMAGES)) DEALLOCATE(THIS%IMAGES)

      READ(IN_UNIT) THIS%IMAGE_COUNT

      ALLOCATE(THIS%IMAGES(3, THIS%IMAGE_COUNT))

      READ(IN_UNIT) THIS%IMAGES
      READ(IN_UNIT) THIS%WEIGHT
   END SUBROUTINE

!********************** FUNCTION PRIM_CELL_INIT ************************
!
! Constructor of PRIM_CELL_T.
!
! LATT_CUR  - Supercell lattice.
! POSION    - Direct coordinates of supercell atoms w.r.t. supercell.
! SYMM      - Contains symmetry information.
!
!***********************************************************************
   TYPE(PRIM_CELL_T) FUNCTION PRIM_CELL_INIT(LATT_CUR, POSION, SYMM) RESULT(PRIM_CELL)
      USE lattice,   ONLY: LATT, LATTIC
      USE base,      ONLY: symmetry
      USE tutor,     ONLY: vtutor
      IMPLICIT NONE

      TYPE(LATT),       INTENT(IN)  :: LATT_CUR
      REAL(q),          INTENT(IN)  :: POSION(:, :)
      TYPE(symmetry),   INTENT(IN)  :: SYMM

      INTEGER  :: S_IDX

      ! Symmetry stuff
      INTEGER  :: ISYMOP, NROT, IGRPOP, NROTK, INVMAP, NPCELL
      REAL(q)  :: GTRANS, AP
      COMMON /SYMM/ ISYMOP(3, 3, 48), NROT, IGRPOP(3, 3, 48), NROTK, &
      GTRANS(3, 48), INVMAP(48), AP(3, 3), NPCELL

      IF (ALLOCATED(PRIM_CELL%POSION_SUPER)) DEALLOCATE(PRIM_CELL%POSION_SUPER)
      IF (ALLOCATED(PRIM_CELL%SUPER_IDX)) DEALLOCATE(PRIM_CELL%SUPER_IDX)
      IF (ALLOCATED(PRIM_CELL%PRIM_IDX)) DEALLOCATE(PRIM_CELL%PRIM_IDX)
      IF (ALLOCATED(PRIM_CELL%CELL_IDX)) DEALLOCATE(PRIM_CELL%CELL_IDX)
      IF (ALLOCATED(PRIM_CELL%MIN_IMAGES)) DEALLOCATE(PRIM_CELL%MIN_IMAGES)

      ! Set some variables for prim cell
      PRIM_CELL%NUM_SUPER_ATOMS = SIZE(POSION, 2)
      PRIM_CELL%NUM_PRIM_ATOMS = PRIM_CELL%NUM_SUPER_ATOMS / NPCELL
      PRIM_CELL%NUM_CELLS = NPCELL

      ! Set atom positions
      ALLOCATE(PRIM_CELL%POSION_SUPER(3, PRIM_CELL%NUM_SUPER_ATOMS))
      PRIM_CELL%POSION_SUPER = POSION

      ! Setup the primitive cell lattice
      PRIM_CELL%SUPER_LATT%A = LATT_CUR%A
      CALL LATTIC(PRIM_CELL%SUPER_LATT)
      PRIM_CELL%PRIM_LATT%A = AP
      PRIM_CELL%PRIM_LATT%SCALE = 1
      CALL LATTIC(PRIM_CELL%PRIM_LATT)

      ! Get trafo from super to primitive cell and vice versa
      PRIM_CELL%SUPER_TO_PRIM = MATMUL(TRANSPOSE(PRIM_CELL%PRIM_LATT%B), PRIM_CELL%SUPER_LATT%A)
      PRIM_CELL%PRIM_TO_SUPER = MATMUL(TRANSPOSE(PRIM_CELL%SUPER_LATT%B), PRIM_CELL%PRIM_LATT%A)

      IF (.NOT. PRIM_CELL%SETUP_SUPER_INDICES(SYMM)) THEN
         ! As a fallback, set the primitive cell to the supercell
         CALL vtutor%alert("Failed to determine mapping between supercell and primitive cell!")

         PRIM_CELL%NUM_PRIM_ATOMS = PRIM_CELL%NUM_SUPER_ATOMS
         PRIM_CELL%NUM_CELLS = 1

         DEALLOCATE(PRIM_CELL%SUPER_IDX)
         ALLOCATE(PRIM_CELL%SUPER_IDX(1, PRIM_CELL%NUM_SUPER_ATOMS))
         DO S_IDX = 1, PRIM_CELL%NUM_SUPER_ATOMS
            PRIM_CELL%SUPER_IDX(1, S_IDX) = S_IDX
         ENDDO

         PRIM_CELL%PRIM_LATT = PRIM_CELL%SUPER_LATT
         PRIM_CELL%SUPER_TO_PRIM = RESHAPE([REAL(q):: 1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
         PRIM_CELL%PRIM_TO_SUPER = PRIM_CELL%SUPER_TO_PRIM
      ENDIF

      CALL PRIM_CELL%SETUP_PRIM_INDICES()
   END FUNCTION

!********************* FUNCTION GET_ID_ROT_INDEX ***********************
!
! Obtains the index to the array ISYMOP that corresponds to the
! identity operator.
!
!***********************************************************************
   INTEGER FUNCTION GET_ID_ROT_INDEX() RESULT(IDX)
      USE tutor,  ONLY: vtutor
      IMPLICIT NONE

      INTEGER, PARAMETER   :: ID(3, 3) = RESHAPE([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])

      INTEGER  :: ISYMOP, NROT, IGRPOP, NROTK, INVMAP, NPCELL
      REAL(q)  :: GTRANS, AP
      COMMON /SYMM/ ISYMOP(3, 3, 48), NROT, IGRPOP(3, 3, 48), NROTK, &
      GTRANS(3, 48), INVMAP(48), AP(3, 3), NPCELL

      DO IDX = 1, NROTK
         IF (ALL(ISYMOP(:, :, IDX) == ID)) RETURN
      ENDDO

      CALL vtutor%bug("Could not find identity operation!", __FILE__, __LINE__)
   END FUNCTION

!********************* FUNCTION FIND_ATOM_REPLICA **********************
!
! Given a supercell atom, find all of its replicas in other primitive
! cells for a certain supercell rotation. Returns whether the search
! was successful (this allows silent failure/fallback and is a
! replacement for throwing an exception).
!
! ATOM_IDX  - Supercell atom index for which to find the replicas
! SYMM      - Contains symmetry information
! ROT_IDX   - Index that corresponds to the symmetry operation stored
!             in array ISYMOP
! REPLICA   - On return, the supercell indices to the replica
!
!***********************************************************************
   LOGICAL FUNCTION FIND_ATOM_REPLICA(ATOM_IDX, SYMM, ROT_IDX, REPLICA) RESULT(SUCCESS)
      USE base,   ONLY: symmetry
      USE tutor,  ONLY: vtutor
      IMPLICIT NONE

      INTEGER,          INTENT(IN)  :: ATOM_IDX
      TYPE(symmetry),   INTENT(IN)  :: SYMM
      INTEGER,          INTENT(IN)  :: ROT_IDX
      INTEGER,          INTENT(OUT) :: REPLICA(:)

      INTEGER  :: NUM_ATOMS, NUM_CELLS
      INTEGER  :: A_IDX
      INTEGER  :: TRANS_IDX
      INTEGER  :: CELL_IDX

      NUM_ATOMS = SIZE(SYMM%ROTMAP, 1)
      NUM_CELLS = SIZE(SYMM%ROTMAP, 3)

      SUCCESS = .TRUE.

      DO CELL_IDX = 1, NUM_CELLS
         DO A_IDX = 1, NUM_ATOMS
            IF (SYMM%ROTMAP(A_IDX, ROT_IDX, CELL_IDX) == ATOM_IDX) THEN
               REPLICA(CELL_IDX) = A_IDX
               EXIT
            ENDIF
         ENDDO
         IF (A_IDX > NUM_ATOMS) THEN
            SUCCESS = .FALSE.
            RETURN
         ENDIF
      ENDDO
   END FUNCTION

!************** SUBROUTINE PRIM_CELL_SETUP_SUPER_INDICES ***************
!
! Creates map from primitive and cell indices to supercell indices.
! Returns whether the operation was successful (this allows silent
! failure/fallback and is a replacement for throwing an exception).
!
!***********************************************************************
   LOGICAL FUNCTION PRIM_CELL_SETUP_SUPER_INDICES(THIS, SYMM) RESULT(SUCCESS)
      USE base,   ONLY: symmetry
      IMPLICIT NONE

      CLASS(PRIM_CELL_T),  INTENT(INOUT)  :: THIS
      TYPE(symmetry),      INTENT(IN)     :: SYMM

      LOGICAL, ALLOCATABLE :: ALREADY_ADDED(:)
      INTEGER              :: ID_IDX
      INTEGER              :: S_IDX
      INTEGER              :: P_IDX
      INTEGER              :: CELL_IDX
      REAL(q)              :: DIST_VEC(3)
      REAL(q)              :: DIST, MIN_DIST
      INTEGER              :: MIN_IDX

      IF (ALLOCATED(THIS%SUPER_IDX)) DEALLOCATE(THIS%SUPER_IDX)
      ALLOCATE(THIS%SUPER_IDX(THIS%NUM_CELLS, THIS%NUM_PRIM_ATOMS))
      THIS%SUPER_IDX = 0

      SUCCESS = .TRUE.

      IF (THIS%NUM_CELLS == 1) THEN
         DO S_IDX = 1, THIS%NUM_SUPER_ATOMS
            THIS%SUPER_IDX(1, S_IDX) = S_IDX
         ENDDO
         RETURN
      ENDIF

      ID_IDX = GET_ID_ROT_INDEX()

      ! Find replicas of first atom
      IF (.NOT. FIND_ATOM_REPLICA(1, SYMM, ID_IDX, THIS%SUPER_IDX(:, 1))) THEN
         SUCCESS = .FALSE.
         RETURN
      ENDIF

      ! Exlcude first atom + replicas from further searches
      ALLOCATE(ALREADY_ADDED(THIS%NUM_SUPER_ATOMS))
      ALREADY_ADDED = .FALSE.
      ALREADY_ADDED(THIS%SUPER_IDX(:, 1)) = .TRUE.
      P_IDX = 1

      ! Search for the next atom
      DO S_IDX = 2, THIS%NUM_SUPER_ATOMS
         ! Skip previously found atoms
         IF (ALREADY_ADDED(S_IDX)) CYCLE

         P_IDX = P_IDX + 1
         ! Find replicas of this atom
         ! The ordering might not correspond to the ones found previously
         IF (.NOT. FIND_ATOM_REPLICA(S_IDX, SYMM, ID_IDX, THIS%SUPER_IDX(:, P_IDX))) THEN
            SUCCESS = .FALSE.
            RETURN
         ENDIF

         ! From the replicas, find the atom that fits into the unit cell at (0, 0, 0)
         ! -> Simply compare distance to first atom
         MIN_DIST = HUGE(1._q)
         DO CELL_IDX = 1, THIS%NUM_CELLS
            DIST_VEC = THIS%GET_POSION(CELL_IDX, P_IDX) - THIS%GET_POSION(1, 1)
            DIST_VEC = MODULO(DIST_VEC + .5_q, 1._q) - .5_q
            DIST = DOT_PRODUCT(DIST_VEC, DIST_VEC)
            IF (DIST < MIN_DIST) THEN
               MIN_DIST = DIST
               MIN_IDX = CELL_IDX
            ENDIF
         ENDDO

         MIN_IDX = THIS%SUPER_IDX(MIN_IDX, P_IDX)
         ! Find the same replicas, but this time the ordering is consistent
         IF (.NOT. FIND_ATOM_REPLICA(MIN_IDX, SYMM, ID_IDX, THIS%SUPER_IDX(:, P_IDX))) THEN
            SUCCESS = .FALSE.
            RETURN
         ENDIF
         ! Exclude from further searches
         ALREADY_ADDED(THIS%SUPER_IDX(:, P_IDX)) = .TRUE.
      ENDDO
   END FUNCTION

!***************** SUBROUTINE PRIM_CELL_SETUP_INDICES ******************
!
! Creates map from supercell indices to primitive and cell indices.
!
!***********************************************************************
   SUBROUTINE PRIM_CELL_SETUP_PRIM_INDICES(THIS)
      IMPLICIT NONE
      CLASS(PRIM_CELL_T),  INTENT(INOUT)  :: THIS

      INTEGER  :: P_IDX
      INTEGER  :: S_IDX
      INTEGER  :: CELL_IDX

      IF (ALLOCATED(THIS%PRIM_IDX)) DEALLOCATE(THIS%PRIM_IDX)
      IF (ALLOCATED(THIS%CELL_IDX)) DEALLOCATE(THIS%CELL_IDX)
      ALLOCATE(THIS%PRIM_IDX(THIS%NUM_SUPER_ATOMS))
      ALLOCATE(THIS%CELL_IDX(THIS%NUM_SUPER_ATOMS))

      DO P_IDX = 1, THIS%NUM_PRIM_ATOMS
         DO CELL_IDX = 1, THIS%NUM_CELLS
            S_IDX = THIS%SUPER_IDX(CELL_IDX, P_IDX)
            THIS%PRIM_IDX(S_IDX) = P_IDX
            THIS%CELL_IDX(S_IDX) = CELL_IDX
         ENDDO
      ENDDO
   END SUBROUTINE


   SUBROUTINE PRIM_CELL_SETUP_ADDITIONAL_DATA(THIS, FORCE_RECALC)
      IMPLICIT NONE
      CLASS(PRIM_CELL_T),  INTENT(INOUT)  :: THIS
      LOGICAL, OPTIONAL,   INTENT(IN)     :: FORCE_RECALC

      LOGICAL  :: FORCE

      FORCE = .FALSE.
      IF (PRESENT(FORCE_RECALC)) FORCE = FORCE_RECALC

      IF (FORCE .OR. .NOT. ALLOCATED(THIS%MIN_IMAGES)) CALL THIS%CALC_MIN_IMAGES()
      IF (FORCE .OR. .NOT. ALLOCATED(THIS%EQUIV_CELLS)) CALL THIS%SETUP_EQUIV_CELLS()
   END SUBROUTINE


!*************** SUBROUTINE PRIM_CELL_SETUP_EQUIV_CELLS ****************
!
! Creates map from two cell indices l and l' to a cell index l'' such
! that (R_{l} - R_{l'}) = (R_{l''} - R_{0}) within the periodic
! supercell.
!
!***********************************************************************
   SUBROUTINE PRIM_CELL_SETUP_EQUIV_CELLS(THIS)
      USE sym_prec,  ONLY  : TINY
      USE tutor,     ONLY  : vtutor
      USE string,    ONLY  : str
      IMPLICIT NONE
      CLASS(PRIM_CELL_T),  INTENT(INOUT)  :: THIS

      INTEGER              :: CELL_IDX, CELL_IDX_1, CELL_IDX_2
      REAL(q)              :: DIST_VEC(3)
      REAL(q), ALLOCATABLE :: DIST_BUFFER(:, :)

      IF (ALLOCATED(THIS%EQUIV_CELLS)) DEALLOCATE(THIS%EQUIV_CELLS)
      ALLOCATE(THIS%EQUIV_CELLS(THIS%NUM_CELLS, THIS%NUM_CELLS))
      THIS%EQUIV_CELLS = 0

      ALLOCATE(DIST_BUFFER(3, THIS%NUM_CELLS))

      DO CELL_IDX = 1, THIS%NUM_CELLS
         THIS%EQUIV_CELLS(CELL_IDX, 1) = CELL_IDX
         DIST_BUFFER(:, CELL_IDX) = MODULO(THIS%GET_POSION(CELL_IDX, 1) - THIS%GET_POSION(1, 1), 1._q)
      ENDDO

      DO CELL_IDX_2 = 2, THIS%NUM_CELLS
         DO CELL_IDX_1 = 1, THIS%NUM_CELLS
            DIST_VEC = MODULO(THIS%GET_POSION(CELL_IDX_1, 1) - THIS%GET_POSION(CELL_IDX_2, 1), 1._q)

            DO CELL_IDX = 1, THIS%NUM_CELLS
               IF (ALL(ABS(DIST_VEC - DIST_BUFFER(:, CELL_IDX)) < TINY)) THEN
                  THIS%EQUIV_CELLS(CELL_IDX_1, CELL_IDX_2) = CELL_IDX
                  EXIT
               ENDIF
            ENDDO

            IF (THIS%EQUIV_CELLS(CELL_IDX_1, CELL_IDX_2) == 0) CALL vtutor%bug( &
               "Could not find equivalent cell for cell indices " // str(CELL_IDX_1) // " and " // &
               str(CELL_IDX_2), __FILE__, __LINE__)
         ENDDO
      ENDDO

   END SUBROUTINE

!******************** FUNCTION PRIM_CELL_GET_POSION ********************
!
! Getter for direct coordinates of primitive atom PRIM_IDX in cell
! CELL_IDX w.r.t. the supercell.
!
!***********************************************************************
   PURE FUNCTION PRIM_CELL_GET_POSION(THIS, CELL_IDX, PRIM_IDX) RESULT(POS)
      IMPLICIT NONE
      REAL(q)                          :: POS(3)
      CLASS(PRIM_CELL_T),  INTENT(IN)  :: THIS
      INTEGER,             INTENT(IN)  :: CELL_IDX
      INTEGER,             INTENT(IN)  :: PRIM_IDX
      POS = THIS%POSION_SUPER(:, THIS%SUPER_IDX(CELL_IDX, PRIM_IDX))
   END FUNCTION

!***************** FUNCTION PRIM_CELL_CALC_MIN_IMAGES ******************
!
! Sets up the minimum-image structures for all relevant pairs of atoms.
!
!***********************************************************************
   SUBROUTINE PRIM_CELL_CALC_MIN_IMAGES(THIS)
      IMPLICIT NONE

      CLASS(PRIM_CELL_T),  INTENT(INOUT)  :: THIS

      INTEGER  :: P_IDX_1, P_IDX_2
      INTEGER  :: CELL_IDX_1
      REAL(q)  :: PRIM_DIST(3)

      IF (ALLOCATED(THIS%MIN_IMAGES)) DEALLOCATE(THIS%MIN_IMAGES)
      ALLOCATE(THIS%MIN_IMAGES(THIS%NUM_CELLS, THIS%NUM_PRIM_ATOMS, THIS%NUM_PRIM_ATOMS))

      DO P_IDX_2 = 1, THIS%NUM_PRIM_ATOMS
         DO P_IDX_1 = 1, THIS%NUM_PRIM_ATOMS
            PRIM_DIST = THIS%GET_POSION(1, P_IDX_1) - THIS%GET_POSION(1, P_IDX_2)
            DO CELL_IDX_1 = 1, THIS%NUM_CELLS
               THIS%MIN_IMAGES(CELL_IDX_1, P_IDX_1, P_IDX_2) = MIN_IMAGE_T(PRIM_DIST, CELL_IDX_1, THIS)
            ENDDO
         ENDDO
      ENDDO
   END SUBROUTINE PRIM_CELL_CALC_MIN_IMAGES

!******************* FUNCTION PRIM_CELL_WRITE_OUTCAR *******************
!
! Writes primitive-cell information to OUTCAR.
!
!***********************************************************************
   SUBROUTINE PRIM_CELL_WRITE_OUTCAR(THIS, IO)
      USE base,   ONLY  : in_struct
      IMPLICIT NONE

      CLASS(PRIM_CELL_T),  INTENT(IN)  :: THIS
      TYPE(in_struct),     INTENT(IN)  :: IO

      REAL(q), ALLOCATABLE :: NICE_POS(:, :)
      INTEGER              :: I, J

      IF (IO%IU0 < 0) RETURN

      ALLOCATE(NICE_POS(3, THIS%NUM_PRIM_ATOMS))

      ! Get position of primitive atoms w.r.t. primitive cell
      NICE_POS = THIS%POSION_SUPER(:, THIS%SUPER_IDX(1, :))
      NICE_POS = MATMUL(THIS%SUPER_TO_PRIM, NICE_POS)
      ! Move positions to primitive-cell parallelepiped for output
      NICE_POS = MODULO(NICE_POS, 1._q)

      WRITE(IO%IU6, '(/, A, //, A, /)') &
         '----------------------------------------------------------------------------------------', &
         '                                     Primitive cell                                     '

      WRITE(IO%IU6, '(A, F12.4, /)') '  volume of cell : ', THIS%PRIM_LATT%OMEGA
      
      WRITE(IO%IU6, '(A, 20X, A)') '  direct lattice vectors', 'reciprocal lattice vectors'
      WRITE(IO%IU6, '(3(2(3X, 3F13.9), /))') ((THIS%PRIM_LATT%A(I, J), I = 1, 3), (THIS%PRIM_LATT%B(I, J), I = 1, 3), J = 1, 3)

      WRITE(IO%IU6, '(A)') '  length of vectors'
      WRITE(IO%IU6, '(2(3X, 3F13.9), /)') (THIS%PRIM_LATT%ANORM(I), I = 1, 3), (THIS%PRIM_LATT%BNORM(I), I = 1, 3)

      WRITE(IO%IU6, '(A)') '  position of ions in fractional coordinates (direct lattice)'

      WRITE(IO%IU6, '(3X, 3F13.9)') ((NICE_POS(I, J), I = 1, 3), J = 1, THIS%NUM_PRIM_ATOMS)

      WRITE(IO%IU6, '(/, A)') '  ion indices of the primitive-cell ions'
      WRITE(IO%IU6, '(A, 3X, A)') '   primitive index', 'ion index'
      WRITE(IO%IU6, '(13X, I5, 7X, I5)') (I, THIS%SUPER_IDX(1, I), I = 1, THIS%NUM_PRIM_ATOMS)

      WRITE(IO%IU6, '(/, A, /)') '----------------------------------------------------------------------------------------'
   END SUBROUTINE
   
!********************* FUNCTION PRIM_CELL_WRITE_XML ********************
!
! Writes primitive-cell information to xml.
!
!***********************************************************************
   SUBROUTINE PRIM_CELL_WRITE_XML(THIS, IO)
      USE vaspxml
      USE base,   ONLY  : in_struct
      IMPLICIT NONE

      CLASS(PRIM_CELL_T),  INTENT(IN)  :: THIS
      TYPE(in_struct),     INTENT(IN)  :: IO

      REAL(q), ALLOCATABLE :: NICE_POS(:, :)

      IF (IO%IU0 < 0) RETURN

      ALLOCATE(NICE_POS(3, THIS%NUM_PRIM_ATOMS))

      ! Get position of primitive atoms w.r.t. primitive cell
      NICE_POS = THIS%POSION_SUPER(:, THIS%SUPER_IDX(1, :))
      NICE_POS = MATMUL(THIS%SUPER_TO_PRIM, NICE_POS)
      ! Move positions to primitive-cell parallelepiped for output
      NICE_POS = MODULO(NICE_POS, 1._q)

      ! Write primitive cell information to XML.
      CALL XML_TAG("primitive_cell")

      ! First, crystal structure.
      CALL XML_TAG("structure", "primitive_cell")
      CALL XML_CRYSTAL(THIS%PRIM_LATT%A, THIS%PRIM_LATT%B, THIS%PRIM_LATT%OMEGA)
      CALL XML_POSITIONS(THIS%NUM_PRIM_ATOMS, NICE_POS)
      CALL XML_CLOSE_TAG("structure")

      ! Second, primitive atom indices.
      CALL XML_VECARRAY("primitive_index")
      CALL XML_ARRAY_INT(THIS%SUPER_IDX(1 : 1, :))
      CALL XML_CLOSE_TAG

      CALL XML_CLOSE_TAG
   ENDSUBROUTINE

!********************* FUNCTION PRIM_CELL_WRITE_POSCAR *****************
!
! Writes primitive-cell information to POSCAR.prim
!
!***********************************************************************
   SUBROUTINE PRIM_CELL_WRITE_POSCAR(THIS,ELEMENTS,T_INFO,IO)
      USE base,   ONLY  : in_struct
      USE poscar, ONLY  : type_info
      IMPLICIT NONE

      CLASS(PRIM_CELL_T),  INTENT(IN)  :: THIS
      CHARACTER(LEN=*),    INTENT(IN)  :: ELEMENTS(:)
      TYPE(type_info),     INTENT(IN)  :: T_INFO
      TYPE(in_struct),     INTENT(IN)  :: IO
      ! local
      INTEGER :: IUNIT
      INTEGER :: I, ISC, ITYP, IONTYP
      REAL(q), ALLOCATABLE :: NICE_POS(:, :)
      INTEGER, ALLOCATABLE :: ATOM_TYPES(:)
      INTEGER, ALLOCATABLE :: PRIM_NITYP(:)

      IF (IO%IU0 < 0) RETURN

      ALLOCATE(NICE_POS(3, THIS%NUM_PRIM_ATOMS))
      ALLOCATE(ATOM_TYPES(THIS%NUM_SUPER_ATOMS))
      ALLOCATE(PRIM_NITYP(T_INFO%NTYP))

       ! Get position of primitive atoms w.r.t. primitive cell
      NICE_POS = THIS%POSION_SUPER(:, THIS%SUPER_IDX(1, :))
      NICE_POS = MATMUL(THIS%SUPER_TO_PRIM, NICE_POS)
      ! Move positions to primitive-cell parallelepiped for output
      NICE_POS = MODULO(NICE_POS, 1._q)

      ! Write primitive cell information to file
      OPEN(NEWUNIT=IUNIT,FILE="POSCAR.prim")
      WRITE(IUNIT,'(A)') 'primitive cell generated by vasp'
      WRITE(IUNIT,'(A)') '1.0'

      ! First, crystal structure.
      WRITE(IUNIT,'(3F16.12)') THIS%PRIM_LATT%A(:,1)
      WRITE(IUNIT,'(3F16.12)') THIS%PRIM_LATT%A(:,2)
      WRITE(IUNIT,'(3F16.12)') THIS%PRIM_LATT%A(:,3)

      ! get array of atoms types in the supercell
      I=0
      DO ITYP=1,T_INFO%NTYP
        DO IONTYP=1,T_INFO%NITYP(ITYP)
           I=I+1
           ATOM_TYPES(I) = ITYP
        ENDDO
      ENDDO

      ! get number of each type in the primitive cell
      PRIM_NITYP = 0
      DO I=1, THIS%NUM_PRIM_ATOMS
         ITYP = ATOM_TYPES(THIS%SUPER_IDX(1, I))
         PRIM_NITYP(ITYP) = PRIM_NITYP(ITYP)+1
      ENDDO

      ! write elements strings
      WRITE(IUNIT,*) ELEMENTS
      WRITE(IUNIT,*) PRIM_NITYP

      ! write positions
      WRITE(IUNIT,'(A)') 'Direct'
      DO I=1, THIS%NUM_PRIM_ATOMS
        ISC = THIS%SUPER_IDX(1, I)
        WRITE(IUNIT,'(3F16.12)') NICE_POS(:,I)
      ENDDO
      CLOSE(IUNIT)

   END SUBROUTINE PRIM_CELL_WRITE_POSCAR

#ifdef VASP_HDF5
   SUBROUTINE PRIM_CELL_WRITE_HDF5(THIS, LOC_ID, GROUP_PATH)
      USE hdf5,         ONLY  : HID_T
      USE vhdf5_base,   ONLY  : vh5_write &
                              , vh5_error
      IMPLICIT NONE

      CLASS(PRIM_CELL_T),     INTENT(IN)  :: THIS
      INTEGER(HID_T),         INTENT(IN)  :: LOC_ID
      CHARACTER(*),           INTENT(IN)  :: GROUP_PATH

      VH5_CHECK(vh5_write(LOC_ID, GROUP_PATH // "num_prim_atoms", THIS%NUM_PRIM_ATOMS))
      VH5_CHECK(vh5_write(LOC_ID, GROUP_PATH // "num_super_atoms", THIS%NUM_SUPER_ATOMS))
      VH5_CHECK(vh5_write(LOC_ID, GROUP_PATH // "num_cells", THIS%NUM_CELLS))
      VH5_CHECK(vh5_write(LOC_ID, GROUP_PATH // "super_lattice", THIS%SUPER_LATT%A))
      VH5_CHECK(vh5_write(LOC_ID, GROUP_PATH // "prim_lattice", THIS%PRIM_LATT%A))
      VH5_CHECK(vh5_write(LOC_ID, GROUP_PATH // "direct_atom_coordinates_super", THIS%POSION_SUPER))
      VH5_CHECK(vh5_write(LOC_ID, GROUP_PATH // "super_atom_indices", THIS%SUPER_IDX))
      VH5_CHECK(vh5_write(LOC_ID, GROUP_PATH // "prim_atom_indices", THIS%PRIM_IDX))
      VH5_CHECK(vh5_write(LOC_ID, GROUP_PATH // "cell_indices", THIS%CELL_IDX))
      VH5_CHECK(vh5_write(LOC_ID, GROUP_PATH // "equivalent_cells", THIS%EQUIV_CELLS))
   END SUBROUTINE

   SUBROUTINE PRIM_CELL_READ_HDF5(THIS, LOC_ID, GROUP_PATH)
      USE hdf5,         ONLY  : HID_T
      USE vhdf5_base,   ONLY  : vh5_read &
                              , vh5_error
      USE LATTICE,      ONLY  : LATTIC
      IMPLICIT NONE

      CLASS(PRIM_CELL_T),  INTENT(INOUT)  :: THIS
      INTEGER(HID_T),         INTENT(IN)  :: LOC_ID
      CHARACTER(*),           INTENT(IN)  :: GROUP_PATH

      IF (ALLOCATED(THIS%POSION_SUPER)) DEALLOCATE(THIS%POSION_SUPER)
      IF (ALLOCATED(THIS%SUPER_IDX)) DEALLOCATE(THIS%SUPER_IDX)
      IF (ALLOCATED(THIS%PRIM_IDX)) DEALLOCATE(THIS%PRIM_IDX)
      IF (ALLOCATED(THIS%CELL_IDX)) DEALLOCATE(THIS%CELL_IDX)
      IF (ALLOCATED(THIS%MIN_IMAGES)) DEALLOCATE(THIS%MIN_IMAGES)

      VH5_CHECK(vh5_read(LOC_ID, GROUP_PATH // "num_prim_atoms", THIS%NUM_PRIM_ATOMS))
      VH5_CHECK(vh5_read(LOC_ID, GROUP_PATH // "num_super_atoms", THIS%NUM_SUPER_ATOMS))
      VH5_CHECK(vh5_read(LOC_ID, GROUP_PATH // "num_cells", THIS%NUM_CELLS))

      VH5_CHECK(vh5_read(LOC_ID, GROUP_PATH // "super_lattice", THIS%SUPER_LATT%A))
      VH5_CHECK(vh5_read(LOC_ID, GROUP_PATH // "prim_lattice", THIS%PRIM_LATT%A))

      CALL LATTIC(THIS%SUPER_LATT)
      CALL LATTIC(THIS%PRIM_LATT)

      ALLOCATE(THIS%POSION_SUPER(3, THIS%NUM_SUPER_ATOMS))
      ALLOCATE(THIS%SUPER_IDX(THIS%NUM_CELLS, THIS%NUM_PRIM_ATOMS))
      ALLOCATE(THIS%PRIM_IDX(THIS%NUM_SUPER_ATOMS))
      ALLOCATE(THIS%CELL_IDX(THIS%NUM_SUPER_ATOMS))
      ALLOCATE(THIS%EQUIV_CELLS(THIS%NUM_CELLS, THIS%NUM_CELLS))

      VH5_CHECK(vh5_read(LOC_ID, GROUP_PATH // "direct_atom_coordinates_super", THIS%POSION_SUPER))
      VH5_CHECK(vh5_read(LOC_ID, GROUP_PATH // "super_atom_indices", THIS%SUPER_IDX))
      VH5_CHECK(vh5_read(LOC_ID, GROUP_PATH // "prim_atom_indices", THIS%PRIM_IDX))
      VH5_CHECK(vh5_read(LOC_ID, GROUP_PATH // "cell_indices", THIS%CELL_IDX))
      VH5_CHECK(vh5_read(LOC_ID, GROUP_PATH // "equivalent_cells", THIS%EQUIV_CELLS))
   END SUBROUTINE
#endif


! _test_mani serialization routine to be replaced with hdf5
   SUBROUTINE PRIM_CELL_SERIALIZE(THIS, IO, OUT_UNIT)
      USE base,   ONLY  : in_struct
      IMPLICIT NONE

      CLASS(PRIM_CELL_T),  INTENT(IN)  :: THIS
      TYPE(in_struct),     INTENT(IN)  :: IO
      INTEGER,             INTENT(IN)  :: OUT_UNIT

      INTEGER  :: P_IDX_1, P_IDX_2
      INTEGER  :: CELL_IDX_1

      IF (IO%IU0 < 0) RETURN

      WRITE(OUT_UNIT) THIS%NUM_PRIM_ATOMS
      WRITE(OUT_UNIT) THIS%NUM_SUPER_ATOMS
      WRITE(OUT_UNIT) THIS%NUM_CELLS
      WRITE(OUT_UNIT) THIS%SUPER_LATT%A
      WRITE(OUT_UNIT) THIS%PRIM_LATT%A
      WRITE(OUT_UNIT) THIS%POSION_SUPER
      WRITE(OUT_UNIT) THIS%SUPER_IDX
      WRITE(OUT_UNIT) THIS%PRIM_IDX
      WRITE(OUT_UNIT) THIS%CELL_IDX
   END SUBROUTINE PRIM_CELL_SERIALIZE

! _test_mani deserialization routine to be replaced with hdf5
   SUBROUTINE PRIM_CELL_DESERIALIZE(THIS, IN_UNIT)
      USE LATTICE,   ONLY  : LATTIC
      IMPLICIT NONE

      CLASS(PRIM_CELL_T),  INTENT(INOUT)  :: THIS
      INTEGER,             INTENT(IN)     :: IN_UNIT

      INTEGER  :: P_IDX_1, P_IDX_2
      INTEGER  :: CELL_IDX_1

      IF (ALLOCATED(THIS%POSION_SUPER)) DEALLOCATE(THIS%POSION_SUPER)
      IF (ALLOCATED(THIS%SUPER_IDX)) DEALLOCATE(THIS%SUPER_IDX)
      IF (ALLOCATED(THIS%PRIM_IDX)) DEALLOCATE(THIS%PRIM_IDX)
      IF (ALLOCATED(THIS%CELL_IDX)) DEALLOCATE(THIS%CELL_IDX)
      IF (ALLOCATED(THIS%MIN_IMAGES)) DEALLOCATE(THIS%MIN_IMAGES)

      READ(IN_UNIT) THIS%NUM_PRIM_ATOMS
      READ(IN_UNIT) THIS%NUM_SUPER_ATOMS
      READ(IN_UNIT) THIS%NUM_CELLS

      READ(IN_UNIT) THIS%SUPER_LATT%A
      CALL LATTIC(THIS%SUPER_LATT)
      READ(IN_UNIT) THIS%PRIM_LATT%A
      CALL LATTIC(THIS%PRIM_LATT)

      ALLOCATE(THIS%POSION_SUPER(3, THIS%NUM_SUPER_ATOMS))
      ALLOCATE(THIS%SUPER_IDX(THIS%NUM_CELLS, THIS%NUM_PRIM_ATOMS))
      ALLOCATE(THIS%PRIM_IDX(THIS%NUM_SUPER_ATOMS))
      ALLOCATE(THIS%CELL_IDX(THIS%NUM_SUPER_ATOMS))

      READ(IN_UNIT) THIS%POSION_SUPER
      READ(IN_UNIT) THIS%SUPER_IDX
      READ(IN_UNIT) THIS%PRIM_IDX
      READ(IN_UNIT) THIS%CELL_IDX
   END SUBROUTINE PRIM_CELL_DESERIALIZE

   !> @brief Allocate arrays in the symm structure
   SUBROUTINE SYMM_ALLOC(SYMM,ISPIN,NIOND)
      USE BASE, ONLY: SYMMETRY
      TYPE (SYMMETRY),INTENT(INOUT) :: SYMM
      INTEGER,INTENT(IN) :: ISPIN
      INTEGER,INTENT(IN) :: NIOND
      ALLOCATE(SYMM%ROTMAP(NIOND,1,1))
      ALLOCATE(SYMM%TAU(NIOND,3))
      ALLOCATE(SYMM%TAUROT(NIOND,3))
      ALLOCATE(SYMM%WRKROT(3*(NIOND+2)))
      ALLOCATE(SYMM%PTRANS(NIOND+2,3))
      ALLOCATE(SYMM%INDROT(NIOND+2))
      IF (ISPIN==2) THEN
         ALLOCATE(SYMM%MAGROT(48,NIOND))
      ELSE
         ALLOCATE(SYMM%MAGROT(1,1))
      ENDIF
   END SUBROUTINE SYMM_ALLOC

   !> @brief Set the values of the symm structure from the common block
   SUBROUTINE SYMM_SET(SYMM)
      USE BASE, ONLY: SYMMETRY
      TYPE (SYMMETRY),INTENT(INOUT) :: SYMM
      INTEGER  ISYMOP,NROT,IGRPOP,NROTK,INVMAP,NPCELL
      REAL(q)  GTRANS,AP
      COMMON /SYMM/   ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL
      SYMM%ISYMOP=ISYMOP
      SYMM%IGRPOP=IGRPOP
      SYMM%GTRANS=GTRANS
      SYMM%INVMAP=INVMAP
      SYMM%AP=AP
      SYMM%NPCELL=NPCELL
   END SUBROUTINE SYMM_SET

!******************** SUBROUTINE INISYM ********************************
!
! subroutine inisym initialises and sets up all symmetry stuff needed by
! the charge-symmetrisation and the force-symmetrisation routines ...
! it sets up: the rotation matrices in direct space (ISYMOP)
!             the rotation matrices in reciprocal space (IGRPOP)
!             and a table indexing their inverse elements (INVMAP)
!             the non-trivial translations for each rotation (GTRANS)
!             a table (ROTMAP) which shows the connections between atoms
!             the number of primitive cells within a supercell (NPCELL)
!             and their coordinates (PTRANS) within the supercell
!             magnetization direction operations for each symmetry (MAGMOM)
! additionally some further output informations will be given on IU6
!
! this routine needs as input information: the lattice vectors (A)
!             the initial atomic positions (POSION)
!             the initial velocities (VEL)
!             the number of atomic species (NTYP)
!             the number of atoms of each species (NITYP)
!             some workarrays (TAU,TAUROT,WRKROT,INDROT)
!             the dimensioning parameter NIOND (max. no. of atoms)
!             the initial atomic magnetic moments (ATOMOM)
!             the number of spin components (ISPIN)
!             and the unit number where to write informations (IU6)
!             ... and good luck (and the hope that there are no bugs)!
!             MOST IMPORTANT: it needs coordinates etc. which have a
!                  precision better than the parameter TINY set in the
!                  SYMLIB or LATTLIB routines (I recommend at least
!                  10 valid digits - if possible take more ...)!!!!!!!
!                  Otherwise the programs might not accept equivalences
!                  between atomic coordinates, velocities etc. and will
!                  not find the correct full symmetry! If you have and
!                  if you want to avoid serious problems, please edit
!                  SYMLIB and LATTLIB and increase parameter TINY ... !
!
!***********************************************************************
#ifdef oldsym
    SUBROUTINE INISYM(A,POSION,VEL,LSFOR,LSDYN,NTYP,NITYP,NIOND, &
   &           PTRANS,NRTK,NPCLL,ROTMAP,TAU,TAUROT,WRKROT,INDROT,ATOMOM,SAXIS,MAGROT,ISPIN,IU6)
      USE prec
      USE random_seeded, ONLY: RANE

      IMPLICIT REAL(q) (A-H,O-Z)
      CHARACTER (4) GRPNAM
      LOGICAL LSDYN,LSFOR(3,NIOND)
      REAL(q) MAGROT(48,NIOND),SAXIS(3),ALPHA,BETA,ATOMOM_(3)
      INTEGER, POINTER ::  ROTMAP(:,:,:)

      DIMENSION A(3,3),POSION(3,NIOND),VEL(3,NIOND),NITYP(NTYP)
      DIMENSION TAU(NIOND,3),TAUROT(NIOND,3),WRKROT(3*(NIOND+2))
      DIMENSION PTRANS(NIOND+2,3),INDROT(NIOND+2),ATOMOM(3*NIOND)

      DIMENSION A1(3),A2(3),A3(3),IOPS(3,3,48),CELDIM(6),PDIM(6)
      DIMENSION COO1(3),COO2(3),B1(3),B2(3),B3(3),P1(3),P2(3),P3(3)

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

      GRPNAM='    '

!=======================================================================
! First me must detect the lattice type and find some set of lattice
! vectors which fulfills some standard crystallographic relations:
!=======================================================================
      DO 1 I=1,3
         A1(I)=A(I,1)
         A2(I)=A(I,2)
         A3(I)=A(I,3)
    1 ENDDO
      CALL LATTYP(A1,A2,A3,IBRAV,CELDIM,IU6)

!=======================================================================
! Find number of atoms and copy atomic positions to work array TAU
!=======================================================================
      NATOMS=0
      DO 2 I=1,NTYP
         NATOMS=NATOMS+NITYP(I)
    2 ENDDO
      DO 3 IA=1,NATOMS
         TAU(IA,1)=POSION(1,IA)
         TAU(IA,2)=POSION(2,IA)
         TAU(IA,3)=POSION(3,IA)
    3 ENDDO

!=======================================================================
! Here the first highlight: Do we have a primitive cell or not? This
! becomes important for symmetrisations on the Fourier grids (i.e. for
! the symmetrisation of the charge density) if we have no pure point
! symmetry because then we need the correct phase shifts when going from
! one grid point to a symmetry related grid point - but this requires
! the knowledge of  ALL  possible 'non-primitive' translations ... !
! And it is of course also needed for all other symmetrizations ...
!=======================================================================
      CALL PRICEL(IBRAV,CELDIM,A(1,1),A(1,2),A(1,3),TAU,P1,P2,P3,PTRANS, &
     &    NPCELL,IPTYP,PDIM,NTYP,1,NITYP,NIOND,TAUROT,INDROT,WRKROT,IU6)
      AP(1,1)=P1(1)
      AP(2,1)=P1(2)
      AP(3,1)=P1(3)
      AP(1,2)=P2(1)
      AP(2,2)=P2(2)
      AP(3,2)=P2(3)
      AP(1,3)=P3(1)
      AP(2,3)=P3(2)
      AP(3,3)=P3(3)

!=======================================================================
! The next step is to copy all coordinates into a temporary array and
! to transform them to the representation in the basis (A1,A2,A3):
!=======================================================================
      DO 4 IA=1,NATOMS
         COO1(1)=POSION(1,IA)
         COO1(2)=POSION(2,IA)
         COO1(3)=POSION(3,IA)
         CALL VECCON(COO1,COO2,1,A(1,1),A(1,2),A(1,3),A1,A2,A3)
         TAU(IA,1)=COO2(1)
         TAU(IA,2)=COO2(2)
         TAU(IA,3)=COO2(3)
    4 ENDDO

!=======================================================================
! Now try to find out which symmetry operations we have ... :
!=======================================================================
      IF (IU6>0) &
      WRITE(IU6,'(/2A)') 'Analysis of symmetry for initial ', &
     &                                         'positions (statically):'
      CALL SETGRP(ISYMOP,IOPS,GTRANS,NROT,NROTK,TAU,IBRAV,NOP, &
     &                      TAUROT,NTYP,1,NITYP,NIOND,INDROT,WRKROT,IU6)
      CALL PGROUP(ISYMOP,NROT,IPGIND,GRPNAM)
      IF (IU6>0) &
      WRITE(IU6,'(/4A)') 'The static configuration has the ', &
     &                                      'point symmetry ',GRPNAM,'.'
      IF (NROTK>NROT) THEN
         CALL PGROUP(ISYMOP,NROTK,IPGIND,GRPNAM)
         IF (IU6>0) &
         WRITE(IU6,*) 'The point group associated with ', &
     &                             'its full space group is ',GRPNAM,'.'
      ENDIF

!=======================================================================
! Now transform all back to the representation in our original basis (A)
! and set up a table showing the symmetry connections between the atoms
!=======================================================================
      CALL SGRCON(ISYMOP,ISYMOP,NROTK,A1,A2,A3,A(1,1),A(1,2),A(1,3))
      CALL VECCON(GTRANS,GTRANS,NROTK,A1,A2,A3,A(1,1),A(1,2),A(1,3))
      DO 5 IA=1,NATOMS
         TAU(IA,1)=POSION(1,IA)
         TAU(IA,2)=POSION(2,IA)
         TAU(IA,3)=POSION(3,IA)
    5 ENDDO
      IF (ASSOCIATED(ROTMAP)) THEN
         DEALLOCATE(ROTMAP)
      ENDIF
      ALLOCATE(ROTMAP(NIOND,NROTK,NPCELL))
      CALL POSMAP(TAU,ISYMOP,GTRANS,NROTK,PTRANS,NPCELL,NIOND,1, &
     &                                NTYP,NIOND,NITYP,ROTMAP(1,1,1),TAUROT)

!=======================================================================
! Now test the dynamical symmetry! Which of all the symmetry operations
! are compatible with the given initial velocities? (Final answer ...)!
! Regenerate the (possibly now invalid) map ROTMAP ...
!=======================================================================
      IF (IU6>0) &
      WRITE(IU6,'(/2A)') 'Analysis of symmetry for dynamics ', &
     &                             '(positions and initial velocities):'
      DO 6 IA=1,NATOMS
         TAU(IA,1)=VEL(1,IA)
         TAU(IA,2)=VEL(2,IA)
         TAU(IA,3)=VEL(3,IA)
    6 ENDDO
      CALL DYNSYM(TAU,ROTMAP(1,1,1),ISYMOP,GTRANS,NROTK,NROT,PTRANS, &
     &         NPCELL,1,NTYP,NIOND,NITYP,NIOND,TAUROT,WRKROT,INDROT,IU6)
      CALL PGROUP(ISYMOP,NROT,IPGIND,GRPNAM)
      IF (IU6>0) &
      WRITE(IU6,'(/4A)') 'The dynamic configuration has the ', &
     &                                      'point symmetry ',GRPNAM,'.'
      IF (NROTK>NROT) THEN
         CALL PGROUP(ISYMOP,NROTK,IPGIND,GRPNAM)
         IF (IU6>0) &
         WRITE(IU6,*) 'The point group associated with ', &
     &                             'its full space group is ',GRPNAM,'.'
      ENDIF
      DO 7 IA=1,NATOMS
         TAU(IA,1)=POSION(1,IA)
         TAU(IA,2)=POSION(2,IA)
         TAU(IA,3)=POSION(3,IA)
    7 ENDDO
      IF (ASSOCIATED(ROTMAP)) THEN
         DEALLOCATE(ROTMAP)
      ENDIF
      ALLOCATE(ROTMAP(NIOND,NROTK,NPCELL))
      CALL POSMAP(TAU,ISYMOP,GTRANS,NROTK,PTRANS,NPCELL,NIOND,1, &
     &                                NTYP,NIOND,NITYP,ROTMAP(1,1,1),TAUROT)

!=======================================================================
! Now test the symmetry constraints due to selective dynamics! Which of
! all the symmetry operations are compatible with the given constraints?
! (Mostly final answer)! Regenerate the (possibly now invalid) map ROTMAP
!=======================================================================
      IF (LSDYN) THEN
         IF (IU6>0) &
         WRITE(IU6,'(/2A)') 'Analysis of constrained symmetry ', &
     &                                'for selective dynamics:'
! Generate a random 'velocity field' ...
         DO IA=1,NATOMS
          DO I=1,3
            WRKROT(3*IA-3+I)=RANE()
          ENDDO
         ENDDO
! ... and symmetrize it according to the unconstrained symmetry:
         CALL VECSYM(WRKROT,ROTMAP(1,1,1),NTYP,NITYP,NIOND,TAUROT,TAU)
! Now set the fixed coordinates to zero ('break symmetry') ...
         DO IA=1,NATOMS
          DO I=1,3
            TAU(IA,I)=WRKROT(3*IA-3+I)
            IF (.NOT.LSFOR(I,IA)) TAU(IA,I)=0._q
          ENDDO
         ENDDO
! ... and check again with routine DYNSYM:
         CALL DYNSYM(TAU,ROTMAP(1,1,1),ISYMOP,GTRANS,NROTK,NROT,PTRANS, &
     &         NPCELL,1,NTYP,NIOND,NITYP,NIOND,TAUROT,WRKROT,INDROT,IU6)
         CALL PGROUP(ISYMOP,NROT,IPGIND,GRPNAM)
         IF (IU6>0) &
         WRITE(IU6,'(/4A)') 'The constrained configuration has the ', &
     &                                      'point symmetry ',GRPNAM,'.'
         IF (NROTK>NROT) THEN
            CALL PGROUP(ISYMOP,NROTK,IPGIND,GRPNAM)
            IF (IU6>0) &
            WRITE(IU6,*) 'The point group associated with ', &
     &                             'its full space group is ',GRPNAM,'.'
         ENDIF
         DO 10 IA=1,NATOMS
            TAU(IA,1)=POSION(1,IA)
            TAU(IA,2)=POSION(2,IA)
            TAU(IA,3)=POSION(3,IA)
   10    ENDDO
         IF (ASSOCIATED(ROTMAP)) THEN
            DEALLOCATE(ROTMAP)
         ENDIF
         ALLOCATE(ROTMAP(NIOND,NROTK,NPCELL))
         CALL POSMAP(TAU,ISYMOP,GTRANS,NROTK,PTRANS,NPCELL,NIOND,1, &
     &                                   NTYP,NIOND,NITYP,ROTMAP(1,1,1),TAUROT)
      ENDIF

!=======================================================================
! Finally test the symmetry constraints due to local magnetic moments
! all the symmetry operations are compatible with the given constraints?
! (Very final answer)! Regenerate the (possibly now invalid) map ROTMAP
!=======================================================================
      IF (ISPIN==2) THEN
         IF (IU6>0) &
         WRITE(IU6,'(/2A)') 'Analysis of magnetic symmetry:'
         CALL MAGSYM(ATOMOM,MAGROT,ROTMAP,ISYMOP,GTRANS,NROTK,NROT,PTRANS, &
     &         NPCELL,1,NTYP,NIOND,NITYP,NIOND,TAUROT,WRKROT,INDROT,IU6)
         CALL PGROUP(ISYMOP,NROT,IPGIND,GRPNAM)
         IF (IU6>0) &
         WRITE(IU6,'(/4A)') 'The magnetic configuration has the ', &
     &                                      'point symmetry ',GRPNAM,'.'
         IF (NROTK>NROT) THEN
            CALL PGROUP(ISYMOP,NROTK,IPGIND,GRPNAM)
            IF (IU6>0) &
            WRITE(IU6,*) 'The point group associated with ', &
     &                             'its full space group is ',GRPNAM,'.'
         ENDIF
         DO 11 IA=1,NATOMS
            TAU(IA,1)=POSION(1,IA)
            TAU(IA,2)=POSION(2,IA)
            TAU(IA,3)=POSION(3,IA)
   11    ENDDO
         IF (ASSOCIATED(ROTMAP)) THEN
            DEALLOCATE(ROTMAP)
         ENDIF
         ALLOCATE(ROTMAP(NIOND,NROTK,NPCELL))
         CALL POSMAP(TAU,ISYMOP,GTRANS,NROTK,PTRANS,NPCELL,NIOND,1, &
     &                                   NTYP,NIOND,NITYP,ROTMAP(1,1,1),TAUROT)
      ELSE IF (ISPIN==4) THEN
!
! this is a quick hack (Juergen what are we supposed to do here ???)
!
         IF (IU6>0) &
         WRITE(IU6,'(/2A)') 'Analysis of non collinear symmetry:'
         CALL RECIPS(1._q,A(1,1),A(1,2),A(1,3),B1,B2,B3)
         CALL EULER(SAXIS,ALPHA,BETA)
         DO 16 IA=1,NATOMS
! Rotate ATOMOM from saxis to the cartesian axes in which
! the integer rotation matrices are defined
            ATOMOM_(1)=COS(BETA)*COS(ALPHA)*ATOMOM(1+(IA-1)*3)- &
                 SIN(ALPHA)*ATOMOM(2+(IA-1)*3)+ &
                 SIN(BETA)*COS(ALPHA)*ATOMOM(3+(IA-1)*3)
            ATOMOM_(2)=COS(BETA)*SIN(ALPHA)*ATOMOM(1+(IA-1)*3)+ &
                 COS(ALPHA)*ATOMOM(2+(IA-1)*3)+ &
                 SIN(BETA)*SIN(ALPHA)*ATOMOM(3+(IA-1)*3)
            ATOMOM_(3)=-SIN(BETA)*ATOMOM(1+(IA-1)*3)+ &
                 COS(BETA)*ATOMOM(3+(IA-1)*3)
! Bring to direct coordinates
            TAU(IA,1)=ATOMOM_(1)*B1(1)+ATOMOM_(2)*B1(2)+ATOMOM_(3)*B1(3)
            TAU(IA,2)=ATOMOM_(1)*B2(1)+ATOMOM_(2)*B2(2)+ATOMOM_(3)*B2(3)
            TAU(IA,3)=ATOMOM_(1)*B3(1)+ATOMOM_(2)*B3(2)+ATOMOM_(3)*B3(3)
 16      ENDDO
         CALL DYNSYM(TAU,ROTMAP(1,1,1),ISYMOP,GTRANS,NROTK,NROT,PTRANS, &
     &        NPCELL,1,NTYP,NIOND,NITYP,NIOND,TAUROT,WRKROT,INDROT,IU6)
         CALL PGROUP(ISYMOP,NROT,IPGIND,GRPNAM)
         IF (IU6>0) &
         WRITE(IU6,'(/4A)') 'The non collinear magnetic moments have a ', &
     &        'point symmetry ',GRPNAM,'.'
         IF (NROTK>NROT) THEN
            CALL PGROUP(ISYMOP,NROTK,IPGIND,GRPNAM)
            IF (IU6>0) &
            WRITE(IU6,*) 'The point group associated with ', &
     &           'its full space group is ',GRPNAM,'.'
         ENDIF
         DO 17 IA=1,NATOMS
            TAU(IA,1)=POSION(1,IA)
            TAU(IA,2)=POSION(2,IA)
            TAU(IA,3)=POSION(3,IA)
 17      ENDDO
         IF (ASSOCIATED(ROTMAP)) THEN
            DEALLOCATE(ROTMAP)
         ENDIF
         ALLOCATE(ROTMAP(NIOND,NROTK,NPCELL))
         CALL POSMAP(TAU,ISYMOP,GTRANS,NROTK,PTRANS,NPCELL,NIOND,1, &
     &                                   NTYP,NIOND,NITYP,ROTMAP(1,1,1),TAUROT)

      ENDIF

!=======================================================================
! Finally set up the transformation matrices in reciprocal space and a
! map of the inverse elements of each (reciprocal space) group element
!======================================================================
      CALL RECIPS(1._q,A(1,1),A(1,2),A(1,3),B1,B2,B3)
      CALL SGRCON(ISYMOP,IGRPOP,NROTK,A(1,1),A(1,2),A(1,3),B1,B2,B3)
      CALL INVGRP(IGRPOP,NROTK,INVMAP)

!=======================================================================
! 'Just for fun': translate all primitive translations GTRANS into the
! generating primitive cell of the supercell (for reasons of beauty):
!=======================================================================
      IF (NPCELL>1) THEN
         DO 12 I=1,NROTK
! Cartesian coordinates of vector ...
            COO1(1)=GTRANS(1,I)*A(1,1)+GTRANS(2,I)*A(1,2)+ &
     &                                                GTRANS(3,I)*A(1,3)
            COO1(2)=GTRANS(1,I)*A(2,1)+GTRANS(2,I)*A(2,2)+ &
     &                                                GTRANS(3,I)*A(2,3)
            COO1(3)=GTRANS(1,I)*A(3,1)+GTRANS(2,I)*A(3,2)+ &
     &                                                GTRANS(3,I)*A(3,3)
! Coordinates in basis (P1,P2,P3) = generating cell ...
            CALL RECIPS(1._q,P1,P2,P3,B1,B2,B3)
            COO2(1)=COO1(1)*B1(1)+COO1(2)*B1(2)+COO1(3)*B1(3)
            COO2(2)=COO1(1)*B2(1)+COO1(2)*B2(2)+COO1(3)*B2(3)
            COO2(3)=COO1(1)*B3(1)+COO1(2)*B3(2)+COO1(3)*B3(3)
! Periodic boundary conditions (translation into the cell) ... : We use
! a shift of 5.E-10 to guarantee that 0. remains zero, also if we have
! (due to numerical errors) a small negative value (lets say -1.E-12 ?)
! which would give a coordinate of 0.999999999999. On the other hand
! values like 0.999999999999 will be treated like 1. and will be shifted
! to 0. (or -1.E-12 ...). It is somehow strange but useful ...
            COO2(1)=MOD(COO2(1)+6._q+5.E-10_q,1._q)-5.E-10_q
            COO2(2)=MOD(COO2(2)+6._q+5.E-10_q,1._q)-5.E-10_q
            COO2(3)=MOD(COO2(3)+6._q+5.E-10_q,1._q)-5.E-10_q
! Cartesian coordinates of translated vector ...
            COO1(1)=COO2(1)*P1(1)+COO2(2)*P2(1)+COO2(3)*P3(1)
            COO1(2)=COO2(1)*P1(2)+COO2(2)*P2(2)+COO2(3)*P3(2)
            COO1(3)=COO2(1)*P1(3)+COO2(2)*P2(3)+COO2(3)*P3(3)
! Coordinates in supercell basis (A) ...
            CALL RECIPS(1._q,A(1,1),A(1,2),A(1,3),B1,B2,B3)
            GTRANS(1,I)=COO1(1)*B1(1)+COO1(2)*B1(2)+COO1(3)*B1(3)
            GTRANS(2,I)=COO1(1)*B2(1)+COO1(2)*B2(2)+COO1(3)*B2(3)
            GTRANS(3,I)=COO1(1)*B3(1)+COO1(2)*B3(2)+COO1(3)*B3(3)
   12    ENDDO
!=======================================================================
! If we have a nonprimitive supercell update once again (and finally)
! the symmetry connection tables for the atomic position including all
! the "primitive" translations as additional symmetry operations ...
!=======================================================================
         DO 13 IA=1,NATOMS
            TAU(IA,1)=POSION(1,IA)
            TAU(IA,2)=POSION(2,IA)
            TAU(IA,3)=POSION(3,IA)
   13    ENDDO
         IF (ASSOCIATED(ROTMAP)) THEN
            DEALLOCATE(ROTMAP)
         ENDIF
         ALLOCATE(ROTMAP(NIOND,NROTK,NPCELL))
         CALL POSMAP(TAU,ISYMOP,GTRANS,NROTK,PTRANS,NPCELL,NIOND,1, &
     &                                   NTYP,NIOND,NITYP,ROTMAP(1,1,1),TAUROT)
      ENDIF

      NRTK=NROTK
      NPCLL=NPCELL

! YEAH! It was a hard job, but now we got all we need ...
      RETURN
    END SUBROUTINE INISYM
#else
#ifndef use_spglib
    SUBROUTINE INISYM(A,POSION,VEL,LSFOR,LSDYN,LPEAD,NTYP,NITYP,NIOND, &
   &           PTRANS,NRTK,NPCLL,ROTMAP,TAU,TAUROT,WRKROT,INDROT,ATOMOM,SAXIS,MAGROT,ISPIN,IU6)
      USE prec
      USE random_seeded, ONLY: RANE
! bexternal__
      USE bexternal
! bexternal__

      IMPLICIT REAL(q) (A-H,O-Z)
      CHARACTER (4) GRPNAM
      LOGICAL LSDYN,LSFOR(3,NIOND),LPEAD
      REAL(q) MAGROT(48,NIOND),SAXIS(3),ALPHA,BETA,ATOMOM_(3)
      INTEGER, POINTER ::  ROTMAP(:,:,:)

      DIMENSION A(3,3),POSION(3,NIOND),VEL(3,NIOND),NITYP(NTYP)
      DIMENSION TAU(NIOND,3),TAUROT(NIOND,3),WRKROT(3*(NIOND+2))
      DIMENSION PTRANS(NIOND+2,3),INDROT(NIOND+2),ATOMOM(3*NIOND)
      DIMENSION TMP(NIOND,3,5),VEC(NIOND,3,3),SPIN(NIOND)

      DIMENSION A1(3),A2(3),A3(3),IOPS(3,3,48),CELDIM(6),PDIM(6)
      DIMENSION COO1(3),COO2(3),B1(3),B2(3),B3(3),P1(3),P2(3),P3(3)

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

      PROFILING_START('inisym')

      GRPNAM='    '

!=======================================================================
! First me must detect the lattice type and find some set of lattice
! vectors which fulfills some standard crystallographic relations:
!=======================================================================
      DO 1 I=1,3
         A1(I)=A(I,1)
         A2(I)=A(I,2)
         A3(I)=A(I,3)
    1 ENDDO
      CALL LATTYP(A1,A2,A3,IBRAV,CELDIM,IU6)
!=======================================================================
! Find number of atoms and copy atomic positions to work array TMP(:,:,1)
!=======================================================================
      NATOMS=0
      DO 2 I=1,NTYP
         NATOMS=NATOMS+NITYP(I)
    2 ENDDO

      TMP=0._q
      DO 3 IA=1,NATOMS
         TMP(IA,1,1)=POSION(1,IA)
         TMP(IA,2,1)=POSION(2,IA)
         TMP(IA,3,1)=POSION(3,IA)
    3 ENDDO
!=======================================================================
! copy all coordinates into a temporary array and transform them to the
! representation in the basis (A1,A2,A3): store in TMP(:,:,2)
!=======================================================================
      DO 4 IA=1,NATOMS
         COO1(1)=POSION(1,IA)
         COO1(2)=POSION(2,IA)
         COO1(3)=POSION(3,IA)
         CALL VECCON(COO1,COO2,1,A(1,1),A(1,2),A(1,3),A1,A2,A3)
         TMP(IA,1,2)=COO2(1)
         TMP(IA,2,2)=COO2(2)
         TMP(IA,3,2)=COO2(3)
    4 ENDDO

!=======================================================================
! Here the first highlight: Do we have a primitive cell or not? This
! becomes important for symmetrisations on the Fourier grids (i.e. for
! the symmetrisation of the charge density) if we have no pure point
! symmetry because then we need the correct phase shifts when going from
! one grid point to a symmetry related grid point - but this requires
! the knowledge of  ALL  possible 'non-primitive' translations ... !
! And it is of course also needed for all other symmetrizations ...
!=======================================================================
      IF (IU6>0) &
      WRITE(IU6,'(/2(/A))') &
     &   'Analysis of symmetry for initial positions (statically):', &
     &   '====================================================================='

      TAU=TMP(:,:,1); VEC=0._q; SPIN=0._q

      CALL PRICELV(IBRAV,CELDIM,A(1,1),A(1,2),A(1,3),TAU,VEC,SPIN,P1,P2,P3,PTRANS, &
     &    NPCELL,IPTYP,PDIM,NTYP,1,NITYP,NIOND,3,TAUROT,INDROT,WRKROT,IU6)

      AP(1,1)=P1(1)
      AP(2,1)=P1(2)
      AP(3,1)=P1(3)
      AP(1,2)=P2(1)
      AP(2,2)=P2(2)
      AP(3,2)=P2(3)
      AP(1,3)=P3(1)
      AP(2,3)=P3(2)
      AP(3,3)=P3(3)

!=======================================================================
! Now try to find out which symmetry operations we have ... :
!=======================================================================
      TAU=TMP(:,:,2); VEC=0._q; SPIN=0._q

      CALL SETGRPV(ISYMOP,IOPS,GTRANS,NROT,NROTK,TAU,VEC,SPIN,IBRAV,NOP, &
     &                      TAUROT,NTYP,1,NITYP,NIOND,3,INDROT,WRKROT,IU6)

      CALL PGROUP(ISYMOP,NROT,IPGIND,GRPNAM)
      IF (IU6>0) &
      WRITE(IU6,'(/4A)') 'The static configuration has the ', &
     &                                      'point symmetry ',GRPNAM,'.'
      IF (NROTK>NROT) THEN
         CALL PGROUP(ISYMOP,NROTK,IPGIND,GRPNAM)
         IF (IU6>0) &
         WRITE(IU6,*) 'The point group associated with ', &
     &                             'its full space group is ',GRPNAM,'.'
      ENDIF

!=======================================================================
! Now test the dynamical symmetry! Which of all the symmetry operations
! are compatible with the given initial velocities? (Final answer ...)!
! Regenerate the (possibly now invalid) map ROTMAP ...
!=======================================================================
      IF (IU6>0) &
      WRITE(IU6,'(/2(/A))') &
     &   'Analysis of symmetry for dynamics (positions and initial velocities):', &
     &   '====================================================================='

! store the velocities in TMP(:,:,3)
      DO 6 IA=1,NATOMS
         TMP(IA,1,3)=VEL(1,IA)
         TMP(IA,2,3)=VEL(2,IA)
         TMP(IA,3,3)=VEL(3,IA)
    6 ENDDO

      TAU=TMP(:,:,1)
      VEC(:,:,1)=TMP(:,:,3)

      CALL PRICELV(IBRAV,CELDIM,A(1,1),A(1,2),A(1,3),TAU,VEC,SPIN,P1,P2,P3,PTRANS, &
     &    NPCELL,IPTYP,PDIM,NTYP,1,NITYP,NIOND,3,TAUROT,INDROT,WRKROT,IU6)

      AP(1,1)=P1(1)
      AP(2,1)=P1(2)
      AP(3,1)=P1(3)
      AP(1,2)=P2(1)
      AP(2,2)=P2(2)
      AP(3,2)=P2(3)
      AP(1,3)=P3(1)
      AP(2,3)=P3(2)
      AP(3,3)=P3(3)

      TAU=TMP(:,:,2)
      VEC(:,:,1)=TMP(:,:,3)
      DO IA=1,NATOMS
         CALL VECCON(VEC(IA,1:3,1),VEC(IA,1:3,1),1,A(1,1),A(1,2),A(1,3),A1,A2,A3)
      ENDDO

      CALL SETGRPV(ISYMOP,IOPS,GTRANS,NROT,NROTK,TAU,VEC,SPIN,IBRAV,NOP, &
     &                      TAUROT,NTYP,1,NITYP,NIOND,3,INDROT,WRKROT,IU6)
      CALL PGROUP(ISYMOP,NROT,IPGIND,GRPNAM)
      IF (IU6>0) &
      WRITE(IU6,'(/4A)') 'The dynamic configuration has the ', &
     &                                      'point symmetry ',GRPNAM,'.'
      IF (NROTK>NROT) THEN
         CALL PGROUP(ISYMOP,NROTK,IPGIND,GRPNAM)
         IF (IU6>0) &
         WRITE(IU6,*) 'The point group associated with ', &
     &                             'its full space group is ',GRPNAM,'.'
      ENDIF

      CALL SGRCON(ISYMOP,ISYMOP,NROTK,A1,A2,A3,A(1,1),A(1,2),A(1,3))
      CALL VECCON(GTRANS,GTRANS,NROTK,A1,A2,A3,A(1,1),A(1,2),A(1,3))

      TAU=TMP(:,:,1)
      IF (ASSOCIATED(ROTMAP)) THEN
         DEALLOCATE(ROTMAP)
      ENDIF
      ALLOCATE(ROTMAP(NIOND,NROTK,NPCELL))
      CALL POSMAP(TAU,ISYMOP,GTRANS,NROTK,PTRANS,NPCELL,NIOND,1, &
     &                                NTYP,NIOND,NITYP,ROTMAP(1,1,1),TAUROT)

!=======================================================================
! Now test the symmetry constraints due to selective dynamics! Which of
! all the symmetry operations are compatible with the given constraints?
! (Mostly final answer)! Regenerate the (possibly now invalid) map ROTMAP
!=======================================================================
      IF (LSDYN) THEN
         IF (IU6>0) &
         WRITE(IU6,'(/2(/A))') &
     &      'Analysis of constrained symmetry for selective dynamics:', &
     &      '====================================================================='

! Generate a random 'velocity field' ...
         DO IA=1,NATOMS
          DO I=1,3
            WRKROT(3*IA-3+I)=RANE()
          ENDDO
         ENDDO
! ... and symmetrize it according to the unconstrained symmetry:
         CALL VECSYM(WRKROT,ROTMAP(1,1,1),NTYP,NITYP,NIOND,TAUROT,TAU)
! Now set the fixed coordinates to zero ('break symmetry') ...
! and store the constrained velocity field in TMP(:,:,4)
         DO IA=1,NATOMS
          DO I=1,3
            TMP(IA,I,4)=WRKROT(3*IA-3+I)
            IF (.NOT.LSFOR(I,IA)) TMP(IA,I,4)=0._q
          ENDDO
         ENDDO

         TAU=TMP(:,:,1)
         VEC(:,:,1)=TMP(:,:,3); VEC(:,:,2)=TMP(:,:,4)

         CALL PRICELV(IBRAV,CELDIM,A(1,1),A(1,2),A(1,3),TAU,VEC,SPIN,P1,P2,P3,PTRANS, &
     &       NPCELL,IPTYP,PDIM,NTYP,1,NITYP,NIOND,3,TAUROT,INDROT,WRKROT,IU6)

         AP(1,1)=P1(1)
         AP(2,1)=P1(2)
         AP(3,1)=P1(3)
         AP(1,2)=P2(1)
         AP(2,2)=P2(2)
         AP(3,2)=P2(3)
         AP(1,3)=P3(1)
         AP(2,3)=P3(2)
         AP(3,3)=P3(3)

         TAU=TMP(:,:,2)
         VEC(:,:,1)=TMP(:,:,3); VEC(:,:,2)=TMP(:,:,4)
         DO IA=1,NATOMS
            CALL VECCON(VEC(IA,1:3,1),VEC(IA,1:3,1),1,A(1,1),A(1,2),A(1,3),A1,A2,A3)
            CALL VECCON(VEC(IA,1:3,2),VEC(IA,1:3,2),1,A(1,1),A(1,2),A(1,3),A1,A2,A3)
         ENDDO

         CALL SETGRPV(ISYMOP,IOPS,GTRANS,NROT,NROTK,TAU,VEC,SPIN,IBRAV,NOP, &
     &                         TAUROT,NTYP,1,NITYP,NIOND,3,INDROT,WRKROT,IU6)
         CALL PGROUP(ISYMOP,NROT,IPGIND,GRPNAM)
         IF (IU6>0) &
         WRITE(IU6,'(/4A)') 'The constrained configuration has the ', &
     &                                         'point symmetry ',GRPNAM,'.'
         IF (NROTK>NROT) THEN
            CALL PGROUP(ISYMOP,NROTK,IPGIND,GRPNAM)
            IF (IU6>0) &
            WRITE(IU6,*) 'The point group associated with ', &
     &                                'its full space group is ',GRPNAM,'.'
         ENDIF

         CALL SGRCON(ISYMOP,ISYMOP,NROTK,A1,A2,A3,A(1,1),A(1,2),A(1,3))
         CALL VECCON(GTRANS,GTRANS,NROTK,A1,A2,A3,A(1,1),A(1,2),A(1,3))
      ENDIF

!=======================================================================
! Test the symmetry constraints due to local magnetic moments! Are
! all symmetry operations compatible with the given constraints?
! Regenerate the (possibly now invalid) map ROTMAP
!=======================================================================
      IF (ISPIN==2) THEN
         IF (IU6>0) &
         WRITE(IU6,'(/2(/A))') &
     &       'Analysis of structural, dynamic, and magnetic symmetry:', &
     &       '====================================================================='

         CALL RECIPS(1._q,A(1,1),A(1,2),A(1,3),B1,B2,B3)
! bexternal__
         IF (LBEXTERNAL()) THEN
            DO IA=1,NATOMS
               ATOMOM(IA)=ATOMOM(IA)+BEXT(1)*100
            ENDDO
         ENDIF
! bexternal__
! Moments along cartesian z-axis, handling like non-collinear moments (see also
! below); bring them to lattice coordinates and store in TMP(:,:,5)
         DO IA=1,NATOMS
            TMP(IA,3,5)=ATOMOM(IA)
         ENDDO

         TAU=TMP(:,:,1)
         VEC(:,:,1)=TMP(:,:,3); VEC(:,:,2)=TMP(:,:,4)
         SPIN=TMP(:,3,5)

         CALL PRICELV(IBRAV,CELDIM,A(1,1),A(1,2),A(1,3),TAU,VEC,SPIN,P1,P2,P3,PTRANS, &
        &    NPCELL,IPTYP,PDIM,NTYP,1,NITYP,NIOND,3,TAUROT,INDROT,WRKROT,IU6)

         AP(1,1)=P1(1)
         AP(2,1)=P1(2)
         AP(3,1)=P1(3)
         AP(1,2)=P2(1)
         AP(2,2)=P2(2)
         AP(3,2)=P2(3)
         AP(1,3)=P3(1)
         AP(2,3)=P3(2)
         AP(3,3)=P3(3)

         TAU=TMP(:,:,2)
         VEC(:,:,1)=TMP(:,:,3); VEC(:,:,2)=TMP(:,:,4)
         DO IA=1,NATOMS
            CALL VECCON(VEC(IA,1:3,1),VEC(IA,1:3,1),1,A(1,1),A(1,2),A(1,3),A1,A2,A3)
            CALL VECCON(VEC(IA,1:3,2),VEC(IA,1:3,2),1,A(1,1),A(1,2),A(1,3),A1,A2,A3)
         ENDDO
         SPIN=TMP(:,3,5)

         CALL SETGRPV(ISYMOP,IOPS,GTRANS,NROT,NROTK,TAU,VEC,SPIN,IBRAV,NOP, &
     &                         TAUROT,NTYP,1,NITYP,NIOND,3,INDROT,WRKROT,IU6)
         CALL PGROUP(ISYMOP,NROT,IPGIND,GRPNAM)
         IF (IU6>0) &
         WRITE(IU6,'(/4A)') 'The overall configuration has the ', &
     &                                      'point symmetry ',GRPNAM,'.'
         IF (NROTK>NROT) THEN
            CALL PGROUP(ISYMOP,NROTK,IPGIND,GRPNAM)
            IF (IU6>0) &
            WRITE(IU6,*) 'The point group associated with ', &
     &                             'its full space group is ',GRPNAM,'.'
         ENDIF

         CALL SGRCON(ISYMOP,ISYMOP,NROTK,A1,A2,A3,A(1,1),A(1,2),A(1,3))
         CALL VECCON(GTRANS,GTRANS,NROTK,A1,A2,A3,A(1,1),A(1,2),A(1,3))
! bexternal__
         IF (LBEXTERNAL()) THEN
            DO IA=1,NATOMS
               ATOMOM(IA)=ATOMOM(IA)-BEXT(1)*100
            ENDDO
         ENDIF
! bexternal__
      ELSE IF (ISPIN==4) THEN
         IF (IU6>0) &
         WRITE(IU6,'(/2(/A))') &
     &      'Analysis of structural, dynamic, and magnetic symmetry:', &
     &      '====================================================================='

         CALL RECIPS(1._q,A(1,1),A(1,2),A(1,3),B1,B2,B3)
         CALL EULER(SAXIS,ALPHA,BETA)
! bexternal__
         IF (LBEXTERNAL()) THEN
            DO IA=1,NATOMS
               ATOMOM(1+(IA-1)*3)=ATOMOM(1+(IA-1)*3)+BEXT(1)*100
               ATOMOM(2+(IA-1)*3)=ATOMOM(2+(IA-1)*3)+BEXT(2)*100
               ATOMOM(3+(IA-1)*3)=ATOMOM(3+(IA-1)*3)+BEXT(3)*100
            ENDDO
         ENDIF
! bexternal__
         DO IA=1,NATOMS
! Rotate ATOMOM from saxis to the cartesian axes in which
! the integer rotation matrices are defined
            ATOMOM_(1)=COS(BETA)*COS(ALPHA)*ATOMOM(1+(IA-1)*3)- &
                 SIN(ALPHA)*ATOMOM(2+(IA-1)*3)+ &
                 SIN(BETA)*COS(ALPHA)*ATOMOM(3+(IA-1)*3)
            ATOMOM_(2)=COS(BETA)*SIN(ALPHA)*ATOMOM(1+(IA-1)*3)+ &
                 COS(ALPHA)*ATOMOM(2+(IA-1)*3)+ &
                 SIN(BETA)*SIN(ALPHA)*ATOMOM(3+(IA-1)*3)
            ATOMOM_(3)=-SIN(BETA)*ATOMOM(1+(IA-1)*3)+ &
                 COS(BETA)*ATOMOM(3+(IA-1)*3)
! Bring to direct coordinates and store in TMP(:,:,5)
            TMP(IA,1,5)=ATOMOM_(1)*B1(1)+ATOMOM_(2)*B1(2)+ATOMOM_(3)*B1(3)
            TMP(IA,2,5)=ATOMOM_(1)*B2(1)+ATOMOM_(2)*B2(2)+ATOMOM_(3)*B2(3)
            TMP(IA,3,5)=ATOMOM_(1)*B3(1)+ATOMOM_(2)*B3(2)+ATOMOM_(3)*B3(3)
         ENDDO

         TAU=TMP(:,:,1)
         VEC(:,:,1)=TMP(:,:,3); VEC(:,:,2)=TMP(:,:,4); VEC(:,:,3)=TMP(:,:,5)

         CALL PRICELV(IBRAV,CELDIM,A(1,1),A(1,2),A(1,3),TAU,VEC,SPIN,P1,P2,P3,PTRANS, &
        &    NPCELL,IPTYP,PDIM,NTYP,1,NITYP,NIOND,3,TAUROT,INDROT,WRKROT,IU6)

         AP(1,1)=P1(1)
         AP(2,1)=P1(2)
         AP(3,1)=P1(3)
         AP(1,2)=P2(1)
         AP(2,2)=P2(2)
         AP(3,2)=P2(3)
         AP(1,3)=P3(1)
         AP(2,3)=P3(2)
         AP(3,3)=P3(3)

         TAU=TMP(:,:,2)
         VEC(:,:,1)=TMP(:,:,3); VEC(:,:,2)=TMP(:,:,4); VEC(:,:,3)=TMP(:,:,5)
         DO IA=1,NATOMS
            CALL VECCON(VEC(IA,1:3,1),VEC(IA,1:3,1),1,A(1,1),A(1,2),A(1,3),A1,A2,A3)
            CALL VECCON(VEC(IA,1:3,2),VEC(IA,1:3,2),1,A(1,1),A(1,2),A(1,3),A1,A2,A3)
            CALL VECCON(VEC(IA,1:3,3),VEC(IA,1:3,3),1,A(1,1),A(1,2),A(1,3),A1,A2,A3)
         ENDDO

         CALL SETGRPV(ISYMOP,IOPS,GTRANS,NROT,NROTK,TAU,VEC,SPIN,IBRAV,NOP, &
     &                         TAUROT,NTYP,1,NITYP,NIOND,3,INDROT,WRKROT,IU6)
         CALL PGROUP(ISYMOP,NROT,IPGIND,GRPNAM)
         IF (IU6>0) &
         WRITE(IU6,'(/4A)') 'The overall configuration has the ', &
        &                                      'point symmetry ',GRPNAM,'.'
         IF (NROTK>NROT) THEN
            CALL PGROUP(ISYMOP,NROTK,IPGIND,GRPNAM)
            IF (IU6>0) &
            WRITE(IU6,*) 'The point group associated with ', &
        &                             'its full space group is ',GRPNAM,'.'
         ENDIF

         CALL SGRCON(ISYMOP,ISYMOP,NROTK,A1,A2,A3,A(1,1),A(1,2),A(1,3))
         CALL VECCON(GTRANS,GTRANS,NROTK,A1,A2,A3,A(1,1),A(1,2),A(1,3))
! bexternal__
         IF (LBEXTERNAL()) THEN
            DO IA=1,NATOMS
               ATOMOM(1+(IA-1)*3)=ATOMOM(1+(IA-1)*3)-BEXT(1)*100
               ATOMOM(2+(IA-1)*3)=ATOMOM(2+(IA-1)*3)-BEXT(2)*100
               ATOMOM(3+(IA-1)*3)=ATOMOM(3+(IA-1)*3)-BEXT(3)*100
            ENDDO
         ENDIF
! bexternal__
      ENDIF

!=======================================================================
! Set up the transformation matrices in reciprocal space
!=======================================================================
      CALL RECIPS(1._q,A(1,1),A(1,2),A(1,3),B1,B2,B3)
      CALL SGRCON(ISYMOP,IGRPOP,NROTK,A(1,1),A(1,2),A(1,3),B1,B2,B3)

!=======================================================================
! Remove all symmetry operation that do not conserve the finite
! difference stencil used in the PEAD method. (Mostly final answer)!
! Regenerate the (possibly now invalid) map ROTMAP.
!=======================================================================
      IF (LPEAD) THEN
         CALL STRCNSRV(REAL(RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/)),KIND=q),ISYMOP,IGRPOP,GTRANS,NROTK,NROT)

         CALL PGROUP(ISYMOP,NROT,IPGIND,GRPNAM)
         IF (IU6>0) &
         WRITE(IU6,'(/4A)') '(PEAD) The star conserving configuration has the ', &
        &                                      'point symmetry ',GRPNAM,'.'
         IF (NROTK>NROT) THEN
            CALL PGROUP(ISYMOP,NROTK,IPGIND,GRPNAM)
            IF (IU6>0) &
            WRITE(IU6,'(/4A)') '(PEAD) The point group associated with ', &
        &                             'its full space group is ',GRPNAM,'.'
         ENDIF
      ENDIF

!=======================================================================
! And finally, set up a map of the inverse elements of each (reciprocal
! space) group element
!=======================================================================
      CALL INVGRP(IGRPOP,NROTK,INVMAP)

!=======================================================================
! 'Just for fun': translate all primitive translations GTRANS into the
! generating primitive cell of the supercell (for reasons of beauty):
!=======================================================================
      DO 12 I=1,NROTK
! Cartesian coordinates of vector ...
         COO1(1)=GTRANS(1,I)*A(1,1)+GTRANS(2,I)*A(1,2)+ &
     &                                             GTRANS(3,I)*A(1,3)
         COO1(2)=GTRANS(1,I)*A(2,1)+GTRANS(2,I)*A(2,2)+ &
     &                                             GTRANS(3,I)*A(2,3)
         COO1(3)=GTRANS(1,I)*A(3,1)+GTRANS(2,I)*A(3,2)+ &
     &                                             GTRANS(3,I)*A(3,3)
! Coordinates in basis (P1,P2,P3) = generating cell ...
         CALL RECIPS(1._q,P1,P2,P3,B1,B2,B3)
         COO2(1)=COO1(1)*B1(1)+COO1(2)*B1(2)+COO1(3)*B1(3)
         COO2(2)=COO1(1)*B2(1)+COO1(2)*B2(2)+COO1(3)*B2(3)
         COO2(3)=COO1(1)*B3(1)+COO1(2)*B3(2)+COO1(3)*B3(3)
! Periodic boundary conditions (--> [-0.5,0.5) )
         COO2(1)=MOD(COO2(1)+100.5_q,1._q)-0.5_q
         COO2(2)=MOD(COO2(2)+100.5_q,1._q)-0.5_q
         COO2(3)=MOD(COO2(3)+100.5_q,1._q)-0.5_q
! Cartesian coordinates of translated vector ...
         COO1(1)=COO2(1)*P1(1)+COO2(2)*P2(1)+COO2(3)*P3(1)
         COO1(2)=COO2(1)*P1(2)+COO2(2)*P2(2)+COO2(3)*P3(2)
         COO1(3)=COO2(1)*P1(3)+COO2(2)*P2(3)+COO2(3)*P3(3)
! Coordinates in supercell basis (A) ...
         CALL RECIPS(1._q,A(1,1),A(1,2),A(1,3),B1,B2,B3)
         GTRANS(1,I)=COO1(1)*B1(1)+COO1(2)*B1(2)+COO1(3)*B1(3)
         GTRANS(2,I)=COO1(1)*B2(1)+COO1(2)*B2(2)+COO1(3)*B2(3)
         GTRANS(3,I)=COO1(1)*B3(1)+COO1(2)*B3(2)+COO1(3)*B3(3)
   12 ENDDO

!=======================================================================
! Update the symmetry connection table (ROTMAP) for the atomic
! positions and the table with spinflips (MAGROT)
!=======================================================================
      TAU=TMP(:,:,1)
      IF (ASSOCIATED(ROTMAP)) THEN
         DEALLOCATE(ROTMAP)
      ENDIF
      ALLOCATE(ROTMAP(NIOND,NROTK,NPCELL))
      CALL POSMAP(TAU,ISYMOP,GTRANS,NROTK,PTRANS,NPCELL,NIOND,1, &
     &                                NTYP,NIOND,NITYP,ROTMAP(1,1,1),TAUROT)

      IF (ISPIN==2) THEN
         NROTKP=NROTK; NROTP=NROT; NPCELLP=NPCELL

         CALL MAGSYM(ATOMOM,MAGROT,ROTMAP,ISYMOP,GTRANS,NROTKP,NROTP,PTRANS, &
     &         NPCELL,1,NTYP,NIOND,NITYP,NIOND,TAUROT,WRKROT,INDROT,-1)

         IF ((NROTKP/=NROTK).OR.(NROTP/=NROT).AND.(NPCELLP/=NPCELL)) THEN
            CALL ERROR('INISYM', 'Unable to resolve symmetry of collinear magnetic &
               &degrees of freedom.', 1)
         ENDIF
      ENDIF

      NRTK=NROTK
      NPCLL=NPCELL

! YEAH! It was a hard job, but now we got all we need ...
      IF (IU6>0) WRITE(IU6,50) NROTK,NROT,NPCELL
   50 FORMAT(//' Subroutine INISYM returns: Found ',I2, &
     &         ' space group operations'/,' (whereof ',I2, &
     &         ' operations are pure point group operations),'/, &
     &         ' and found ',I5,' ''primitive'' translations'/)

      PROFILING_STOP('inisym')

      RETURN
    END SUBROUTINE INISYM
#else
    SUBROUTINE INISYM(A,POSION,VEL,LSFOR,LSDYN,LPEAD,NTYP,NITYP,NIOND, &
   &           PTRANS,NRTK,NPCLL,ROTMAP,TAU,TAUROT,WRKROT,INDROT,ATOMOM,SAXIS,MAGROT,ISPIN,IU6)
      USE prec
      USE sym_prec
      USE spglib_f08
      USE iso_c_binding
      USE random_seeded, ONLY: RANE
! bexternal__
      USE bexternal, ONLY : LBEXTERNAL, BEXT
! bexternal__
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: NIOND, NTYP, ISPIN, IU6
      INTEGER, INTENT(IN)  :: NITYP(NTYP)
      REAL(q), INTENT(IN)  :: SAXIS(3)
      REAL(q), INTENT(IN)  :: A(3, 3)
      REAL(q), INTENT(IN)  :: POSION(3, NIOND)
      REAL(q), INTENT(IN)  :: VEL(3, NIOND)
      LOGICAL, INTENT(IN)  :: LSDYN
      LOGICAL, INTENT(IN)  :: LSFOR(3, NIOND)
      LOGICAL, INTENT(IN)  :: LPEAD
      INTEGER, INTENT(OUT) :: NRTK, NPCLL
      REAL(q), INTENT(OUT) :: PTRANS(NIOND+2, 3)
      REAL(q), INTENT(OUT) :: MAGROT(48, NIOND)

      REAL(q), INTENT(INOUT) :: ATOMOM(3 * NIOND)    ! Magnetic moments for collinear and non-collinear
      REAL(q), INTENT(INOUT) :: WRKROT(3*(NIOND+2))  ! Random numbers used for selective dynamics

      INTEGER, INTENT(IN)  :: INDROT(NIOND+2)        ! not used, only here to preserve the calling interface

      INTEGER, POINTER :: ROTMAP(:,:,:)

      INTEGER :: ISYMOP, NROT, IGRPOP, NROTK, INVMAP, NPCELL
      REAL(q) :: GTRANS, AP

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

      ! local variables
      INTEGER :: I, J, K
      INTEGER, ALLOCATABLE :: ISYMOP_TMP(:,:,:)
      REAL(q), ALLOCATABLE :: GTRANS_TMP(:,:)
      INTEGER, ALLOCATABLE :: SPIN_FLIPS(:)
      INTEGER, ALLOCATABLE :: ATOM_TYPES(:)
      INTEGER, ALLOCATABLE :: EQUIV_ATOMS(:)
      REAL(q) :: LATTICE(3,3)
      INTEGER :: NATOMS, NSYM_TMP

      REAL(q) :: TAU(NIOND,3), TAUROT(NIOND,3)

      CHARACTER(LEN=6) ::HermannMaguinSymbol

      INTEGER :: IRETURN

      REAL(q) :: CELDIM(6)
      INTEGER :: IPGIND,IBRAV
      CHARACTER(LEN=4) :: GRPNAM='    '

      TYPE(SpglibDataset) :: DATASET
      TYPE(SpglibSpacegroupType) :: SPGTYPE


      IF (IU6>0) WRITE(IU6,'(/A/,/X,A/)') &
         '----------------------------------------------------------------------------------------', &
         'INISYM: using SPGLIB'

      NATOMS = 0
      DO I = 1, NTYP
         NATOMS = NATOMS + NITYP(I)
      ENDDO

      ALLOCATE(ATOM_TYPES(NATOMS))
      ALLOCATE(EQUIV_ATOMS(NATOMS))

      ! Transpose because of memory order difference
      LATTICE = TRANSPOSE(A)

      ! As input, spglib requires an unique number for each atomic
      ! species for all atoms.
      K = 1
      DO I = 1, NTYP
         DO J = 1, NITYP(I)
            ATOM_TYPES(K) = I
            K = K + 1
         ENDDO
      ENDDO

      ! Symmetry search. For normal space group operations, spglib dataset is
      ! searched for getting more crystallographic information.

      DATASET = spg_get_dataset(LATTICE, POSION, ATOM_TYPES, NATOMS, TINY)

      AP = TRANSPOSE(DATASET % primitive_lattice)

      IF (IU6>0) THEN
         WRITE(IU6,'(/X,A,/X,A)') &
           'SPGLIB: symmetry analysis of lattice + basis:', &
           '------------------------------------------------------'
         WRITE(IU6,'(A,X,A17,6X,A,I4)') &
           ' International :', TRIM( DATASET % international_symbol ), 'Nr.', DATASET % spacegroup_number
         WRITE(IU6,'(A,X,A17,6X,A,I4)') &
           ' Hall          :', TRIM( DATASET % hall_symbol ),          'Nr.', DATASET % hall_number

         SPGTYPE = spg_get_spacegroup_type( DATASET % hall_number )

         WRITE(IU6,'(A,X,A17)') &
           ' Schoenflies   :', TRIM( SPGTYPE % Schoenflies )

         WRITE(IU6,'(/A,I3)') &
           ' Nr. of symops.:', DATASET % n_operations

         WRITE(IU6,'(/A/,3(/X,A,3(F15.10,A)))') &
           ' Primitive cell:', &
           'A1 = (', AP(1,1), ',', AP(2,1), ',', AP(3,1), ')', &
           'A2 = (', AP(1,2), ',', AP(2,2), ',', AP(3,2), ')', &
           'A3 = (', AP(1,3), ',', AP(2,3), ',', AP(3,3), ')'
      ENDIF

      NSYM_TMP = DATASET % n_operations

      ALLOCATE(ISYMOP_TMP(3, 3, NSYM_TMP))
      ALLOCATE(GTRANS_TMP(3, NSYM_TMP))
      ALLOCATE(SPIN_FLIPS(NSYM_TMP))

      ISYMOP_TMP = DATASET % rotations
      GTRANS_TMP = DATASET % translations

      CALL SPGLIB_VELOCITY(ISYMOP_TMP, GTRANS_TMP, NSYM_TMP, AP, &
           EQUIV_ATOMS, NATOMS, NIOND, LATTICE, POSION, ATOM_TYPES, VEL, IU6)

      IF (IU6>0) WRITE(IU6,'(//X,A,I3,X,A)') &
            'SPGLIB: ... accounting for velocities          ... ', NSYM_TMP,'symops. remaining.'

      IF (LSDYN) THEN
         ! Generate a random 'velocity field' ...
         DO J = 1, NATOMS
            DO I = 1, 3
               WRKROT(3 * J - 3 + I) = RANE()
            ENDDO
         ENDDO

         CALL SPGLIB_SELECTIVE_DYNAMICS(ISYMOP_TMP, GTRANS_TMP, NSYM_TMP, AP, &
              EQUIV_ATOMS, NATOMS, NIOND, LATTICE, POSION, ATOM_TYPES, &
              LSDYN, LSFOR, WRKROT, IU6)

         IF (IU6>0) WRITE(IU6,'(/X,A,I3,X,A)') &
            'SPGLIB: ... accounting for selective dynamics  ... ', NSYM_TMP,'symops. remaining.'
      ENDIF

      IF ((ISPIN == 2) .OR. (ISPIN == 4)) THEN
! bexternal__
         IF (LBEXTERNAL()) THEN
            IF (ISPIN==2) THEN
               DO I = 1, NATOMS
                  ATOMOM(I) = ATOMOM(I) + BEXT(1)*100
               ENDDO
            ELSE IF(ISPIN==4) THEN
               DO I = 1, NATOMS
                  ATOMOM(1 + (I-1) * 3) = ATOMOM(1 + (I-1) * 3) + BEXT(1)*100
                  ATOMOM(2 + (I-1) * 3) = ATOMOM(2 + (I-1) * 3) + BEXT(2)*100
                  ATOMOM(3 + (I-1) * 3) = ATOMOM(3 + (I-1) * 3) + BEXT(3)*100
               ENDDO
            ENDIF
         ENDIF
! bexternal__

         CALL SPGLIB_MAGNETIC(ISYMOP_TMP, GTRANS_TMP, SPIN_FLIPS, NSYM_TMP, &
              AP, EQUIV_ATOMS, NATOMS, NIOND, LATTICE, POSION, ATOM_TYPES, &
              ATOMOM, ISPIN, SAXIS, IU6)

! bexternal__
         IF (LBEXTERNAL()) THEN
            IF (ISPIN==2) THEN
               DO I = 1, NATOMS
                  ATOMOM(I) = ATOMOM(I) - BEXT(1)*100
               ENDDO
            ELSE IF(ISPIN==4) THEN
               DO I = 1, NATOMS
                  ATOMOM(1 + (I-1) * 3) = ATOMOM(1 + (I-1) * 3) - BEXT(1)*100
                  ATOMOM(2 + (I-1) * 3) = ATOMOM(2 + (I-1) * 3) - BEXT(2)*100
                  ATOMOM(3 + (I-1) * 3) = ATOMOM(3 + (I-1) * 3) - BEXT(3)*100
               ENDDO
            ENDIF
         ENDIF
! bexternal__

         IF (IU6>0) WRITE(IU6,'(/X,A,I3,X,A)') &
            'SPGLIB: ... accounting for magnetic structure  ... ', NSYM_TMP,'symops. remaining.'
      ENDIF

      CALL SPGLIB_TO_VASP(ISYMOP, GTRANS, PTRANS, IGRPOP, INVMAP, MAGROT, NROT, &
           NROTK, NPCELL, A, AP, ISYMOP_TMP, GTRANS_TMP, SPIN_FLIPS, NSYM_TMP, &
           NATOMS, NIOND, ISPIN, IU6)

      IF (LPEAD) THEN
         CALL STRCNSRV(REAL(RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/)),KIND=q),ISYMOP,IGRPOP,GTRANS,NROTK,NROT)

         IF (IU6>0) WRITE(IU6,'(/X,A,I3,X,A)') &
            'PEAD  : ... conserve finite difference stencil ... ', NROTK*NPCELL,'symops. remaining.'

         CALL INVGRP(IGRPOP,NROTK,INVMAP)
      ENDIF

      ! Generate ROTMAP
      DO I=1,NATOMS
         TAU(I,1)=POSION(1,I)
         TAU(I,2)=POSION(2,I)
         TAU(I,3)=POSION(3,I)
      ENDDO

      IF (ASSOCIATED(ROTMAP)) DEALLOCATE(ROTMAP)
      ALLOCATE(ROTMAP(NIOND, NROTK, NPCELL))

      CALL POSMAP(TAU, ISYMOP, GTRANS, NROTK, PTRANS, NPCELL, NIOND, 1, &
           & NTYP, NIOND, NITYP, ROTMAP(1,1,1), TAUROT)

!!      CALL SPGLIB_GET_PRIMITIVE(NATOMS, NIOND, NROTK, NPCELL, ROTMAP, POSION, A, AP, IU6)

      CALL SPGLIB_GET_POINTGROUP( ISYMOP, NROTK, IRETURN, HermannMaguinSymbol, -1)
      CALL PGROUP( ISYMOP, NROTK, IPGIND, GRPNAM)

      IF (IPGIND == IPOINTGROUP_SPGLIB2VASP(IRETURN)) THEN
         IF (IU6>0) WRITE(IU6,'(//X,A,/5X,A,A6,4X,A,X,A4)') &
            'Pointgroup associated with remaining symops.:','Hermann-Maguin (SPGLIB):',TRIM(HermannMaguinSymbol), &
            'Schoenflies (PGROUP):',GRPNAM
      ELSE
         IF (IU6>0) WRITE(IU6,'(/X,A,/5X,A,A6)') &
            'Pointgroup associated with remaining symops.:','Hermann-Maguin (SPGLIB):',TRIM(HermannMaguinSymbol)
      ENDIF

      IF (IU6>0) THEN
         WRITE(IU6,'(//X,A,/X,A)') &
           'LATTYP: primitive cell (standardized):', &
           '------------------------------------------------------'
         CALL LATTYP(AP(:,1),AP(:,2),AP(:,3),IBRAV,CELDIM,IU6)

         IF (NPCELL>1) THEN
            WRITE(IU6,'(/X,I5,X,A/)') NPCELL,'primitive cells build up your supercell'
         ELSE
            AP=A
!!            WRITE(IU6,'(/X,A/,3(/X,A,3(F15.10,A)))') &
!!              ' Lattice vectors:', &
!!              'A1 = (', AP(1,1), ',', AP(2,1), ',', AP(3,1), ')', &
!!              'A2 = (', AP(1,2), ',', AP(2,2), ',', AP(3,2), ')', &
!!              'A3 = (', AP(1,3), ',', AP(2,3), ',', AP(3,3), ')'
            WRITE(IU6,'(/5X,A/)') 'Original cell was already a primitive cell'
         ENDIF
      ENDIF

      DEALLOCATE(GTRANS_TMP)
      DEALLOCATE(ISYMOP_TMP)
      DEALLOCATE(SPIN_FLIPS)
      DEALLOCATE(ATOM_TYPES)
      DEALLOCATE(EQUIV_ATOMS)

      NRTK = NROTK
      NPCLL = NPCELL

      IF (IU6>0) WRITE(IU6,'(/X,A)') 'INISYM: done'

      RETURN
    END SUBROUTINE INISYM


    SUBROUTINE SPGLIB_VELOCITY(ISYMOP_TMP, GTRANS_TMP, NSYM_TMP, AP, &
           EQUIV_ATOMS, NATOMS, NIOND, LATTICE, POSION, ATOM_TYPES, VEL, IU6)
      USE sym_prec
      USE spglib_f08
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NATOMS, NIOND, IU6
      INTEGER, INTENT(IN) :: ATOM_TYPES(NATOMS)
      REAL(q), INTENT(IN) :: LATTICE(3, 3)  ! Basis vectors are given as row vectors.
      REAL(q), INTENT(IN) :: POSION(3, NIOND)
      REAL(q), INTENT(IN) :: VEL(3, NIOND)  ! Initial velocities of atoms
      INTEGER, INTENT(INOUT) :: NSYM_TMP
      INTEGER, INTENT(INOUT) :: ISYMOP_TMP(3, 3, NSYM_TMP)
      REAL(q), INTENT(INOUT) :: GTRANS_TMP(3, NSYM_TMP)
      REAL(q), INTENT(OUT) :: AP(3, 3)  ! Basis vectors are given as column vectors.
      INTEGER, INTENT(OUT) :: EQUIV_ATOMS(NATOMS)

      REAL(q), ALLOCATABLE :: TENSORS(:)
      REAL(q) :: PRIMITIVE_LATTICE(3,3) ! To store transpose(AP)
      REAL(q) :: VELOCITIES(3, NATOMS)
      INTEGER :: I
      INTEGER :: IS_MAGNETIC
      INTEGER :: TENSOR_RANK
      INTEGER, ALLOCATABLE :: DUMMY(:)

      IS_MAGNETIC = 0
      TENSOR_RANK = 1

      ALLOCATE(TENSORS(NATOMS * 3))

      VELOCITIES = MATMUL(TRANSPOSE(LATTICE), VEL(:, 1:NATOMS))
      DO I = 1, NATOMS
         TENSORS((1 + (I - 1) * 3):(3 + (I - 1) * 3)) = VELOCITIES(:, I)
      ENDDO

      NSYM_TMP = spg_get_symmetry_with_site_tensors( &
           & ISYMOP_TMP, GTRANS_TMP, EQUIV_ATOMS, PRIMITIVE_LATTICE, &
           & DUMMY, NSYM_TMP, LATTICE, POSION, ATOM_TYPES, TENSORS, &
           & TENSOR_RANK, NATOMS, IS_MAGNETIC, TINY)
      AP = TRANSPOSE(PRIMITIVE_LATTICE)

      DEALLOCATE(TENSORS)

      RETURN
    END SUBROUTINE SPGLIB_VELOCITY


    SUBROUTINE SPGLIB_SELECTIVE_DYNAMICS(ISYMOP_TMP, GTRANS_TMP, NSYM_TMP, AP, &
           EQUIV_ATOMS, NATOMS, NIOND, LATTICE, POSION, ATOM_TYPES, LSDYN, &
           LSFOR, WRKROT, IU6)
      USE sym_prec
      USE spglib_f08
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NATOMS, NIOND, IU6
      INTEGER, INTENT(IN) :: ATOM_TYPES(NATOMS)
      REAL(q), INTENT(IN) :: LATTICE(3, 3)  ! Basis vectors are given as row vectors.
      REAL(q), INTENT(IN) :: POSION(3, NIOND)
      LOGICAL, INTENT(IN) :: LSDYN
      LOGICAL, INTENT(IN) :: LSFOR(3, NIOND)
      REAL(q), INTENT(IN) :: WRKROT(3*(NIOND+2))  ! Random numbers used for selective dynamics
      INTEGER, INTENT(INOUT) :: NSYM_TMP
      INTEGER, INTENT(INOUT) :: ISYMOP_TMP(3, 3, NSYM_TMP)
      REAL(q), INTENT(INOUT) :: GTRANS_TMP(3, NSYM_TMP)
      REAL(q), INTENT(OUT) :: AP(3, 3)  ! Basis vectors are given as column vectors.
      INTEGER, INTENT(OUT) :: EQUIV_ATOMS(NATOMS)

      REAL(q), ALLOCATABLE :: TENSORS(:)
      REAL(q) :: PRIMITIVE_LATTICE(3,3) ! To store transpose(AP)
      REAL(q) :: VEL(3, NATOMS), VELSUM(3, NATOMS)  ! Direct coordinates
      REAL(q) :: VELOCITIES(3, NATOMS)  ! Cartesian coordinates
      REAL(q) :: POSION_ROTATED(3, NIOND)
      REAL(q) :: DIFF(3)
      REAL(q) :: TOLERANCE
      INTEGER :: I, J, K, L
      INTEGER :: IS_MAGNETIC
      INTEGER :: TENSOR_RANK
      INTEGER :: IND_MAP(NSYM_TMP, NATOMS)
      INTEGER, ALLOCATABLE :: DUMMY(:)

      IS_MAGNETIC = 0
      TENSOR_RANK = 1

      DO I = 1, NATOMS
         DO J = 1, 3
            VEL(J, I) = WRKROT((I - 1) * 3 + J)
            IF (.NOT. LSFOR(J, I)) VEL(J, I) = 0._q
         ENDDO
      ENDDO

      ! Atomic index mapping table by all space group operations in supercell
      ! (R_i,t_i)x_j = x_k --> IND_MAP(i, j) = k
      IND_MAP(:, :) = 0
      DO I = 1, NSYM_TMP
         POSION_ROTATED = MATMUL(TRANSPOSE(ISYMOP_TMP(:, :, I)), POSION)
         DO J = 1, NATOMS
            POSION_ROTATED(:, J) = POSION_ROTATED(:, J) + GTRANS_TMP(:, I)
         ENDDO
         DO J = 1, NATOMS  ! original atom

            ! Loop over K is to loosen tolerance to find DIFF(:) = 0 for
            ! slightly distorted case.
            TOLERANCE = TINY
            DO K = 1, 100
               DO L = 1, NATOMS  ! rotated atom
                  DIFF = POSION_ROTATED(:, L) - POSION(:, J)
                  DIFF = DIFF - NINT(DIFF)
                  DIFF = MATMUL(TRANSPOSE(LATTICE), DIFF)
                  IF ((ABS(DIFF(1)) < TOLERANCE) .AND. &
                       & (ABS(DIFF(2)) < TOLERANCE) .AND. &
                       & (ABS(DIFF(3)) < TOLERANCE)) THEN
                     IND_MAP(I, J) = L
                     EXIT
                  ENDIF
               ENDDO
               IF (IND_MAP(I, J) == 0) THEN
                  TOLERANCE = TOLERANCE * 1.05
               ELSE
                  EXIT
               ENDIF
            ENDDO

         ENDDO
      ENDDO

      ! A trivial consistency check of IND_MAP
      DO I = 1, NSYM_TMP
         DO J = 1, NATOMS  ! original atom
            IF (IND_MAP(I, J) == 0) THEN
               WRITE(IU6, *) 'warning: IND_MAP is wrongly assinged.'
            ENDIF
         ENDDO
      ENDDO

      ! Sum over all vectors generated by rotation
      VELSUM(:, :) = 0
      DO I = 1, NSYM_TMP
         DO J = 1, NATOMS  ! original atom
            K = IND_MAP(I, J)  ! rotated atom
            VELSUM(:, K) = VELSUM(:, K) + &
                 & MATMUL(TRANSPOSE(ISYMOP_TMP(:, :, I)), VEL(:, J))
         ENDDO
      ENDDO

      ! Break symmetry by setting zero at the elements corresponding to
      ! T or F in POSCAR
      DO I = 1, NATOMS
         DO J = 1, 3
            IF (.NOT. LSFOR(J, I)) VELSUM(J, I) = 0._q
         ENDDO
      ENDDO

      ALLOCATE(TENSORS(NATOMS * 3))

      ! Convert to Cartesian coordinates
      VELOCITIES = MATMUL(TRANSPOSE(LATTICE), VELSUM(:, 1:NATOMS))
      DO I = 1, NATOMS
         TENSORS((1 + (I - 1) * 3):(3 + (I - 1) * 3)) = VELOCITIES(:, I)
      ENDDO

      NSYM_TMP = spg_get_symmetry_with_site_tensors( &
           & ISYMOP_TMP, GTRANS_TMP, EQUIV_ATOMS, PRIMITIVE_LATTICE, &
           & DUMMY, NSYM_TMP, LATTICE, POSION, ATOM_TYPES, TENSORS, &
           & TENSOR_RANK, NATOMS, IS_MAGNETIC, TINY)
      AP = TRANSPOSE(PRIMITIVE_LATTICE)

      DEALLOCATE(TENSORS)

      RETURN
    END SUBROUTINE SPGLIB_SELECTIVE_DYNAMICS


    SUBROUTINE SPGLIB_MAGNETIC(ISYMOP_TMP, GTRANS_TMP, SPIN_FLIPS, NSYM_TMP, &
           AP, EQUIV_ATOMS, NATOMS, NIOND, LATTICE, POSION, ATOM_TYPES, &
           ATOMOM, ISPIN, SAXIS, IU6)
      USE sym_prec
      USE spglib_f08
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NATOMS, NIOND, ISPIN, IU6
      INTEGER, INTENT(IN) :: ATOM_TYPES(NATOMS)
      REAL(q), INTENT(IN) :: SAXIS(3)
      REAL(q), INTENT(IN) :: LATTICE(3,3)  ! Basis vectors are given as row vectors.
      REAL(q), INTENT(IN) :: POSION(3, NIOND)
      REAL(q), INTENT(IN) :: ATOMOM(3 * NIOND)  ! Magnetic moments for collinear and non-collinear
      INTEGER, INTENT(INOUT) :: NSYM_TMP
      INTEGER, INTENT(INOUT) :: ISYMOP_TMP(3, 3, NSYM_TMP)
      REAL(q), INTENT(INOUT) :: GTRANS_TMP(3, NSYM_TMP)
      REAL(q), INTENT(OUT) :: AP(3,3)  ! Basis vectors are given as column vectors.
      INTEGER, INTENT(OUT) :: EQUIV_ATOMS(NATOMS)
      INTEGER, INTENT(OUT) :: SPIN_FLIPS(NSYM_TMP)

      REAL(q), ALLOCATABLE :: TENSORS(:)
      REAL(q) :: ALPHA, BETA
      REAL(q) :: PRIMITIVE_LATTICE(3,3) ! To store transpose(AP)
      REAL(q) :: ATOMOM_(3 * NATOMS)
      INTEGER :: I, J
      INTEGER :: TENSOR_RANK
      INTEGER :: IS_MAGNETIC  ! 0 (false) or 1 (true)

      IF (ISPIN == 2) THEN
         TENSOR_RANK = 0
         IS_MAGNETIC = 1
         ALLOCATE(TENSORS(NATOMS))
         DO I = 1, NATOMS
            TENSORS(I) = ATOMOM(I)
         ENDDO
      ELSE IF (ISPIN == 4) THEN
         ALLOCATE(TENSORS(NATOMS * 3))

         ! ---- from INISYM ----
         CALL EULER(SAXIS, ALPHA, BETA)
         DO I = 1, NATOMS
            TENSORS(1 + (I - 1) * 3)=COS(BETA)*COS(ALPHA)*ATOMOM(1+(I-1)*3)- &
                 SIN(ALPHA)*ATOMOM(2+(I-1)*3)+ &
                 SIN(BETA)*COS(ALPHA)*ATOMOM(3+(I-1)*3)
            TENSORS(2 + (I - 1) * 3)=COS(BETA)*SIN(ALPHA)*ATOMOM(1+(I-1)*3)+ &
                 COS(ALPHA)*ATOMOM(2+(I-1)*3)+ &
                 SIN(BETA)*SIN(ALPHA)*ATOMOM(3+(I-1)*3)
            TENSORS(3 + (I - 1) * 3)=-SIN(BETA)*ATOMOM(1+(I-1)*3)+ &
                 COS(BETA)*ATOMOM(3+(I-1)*3)
         ENDDO
         ! ---- from INISYM ----

         TENSOR_RANK = 1
         IS_MAGNETIC = 1
      ENDIF

      ! spglib takes ATOMOM in cartesian coordinates.
      NSYM_TMP = spg_get_symmetry_with_site_tensors( &
           & ISYMOP_TMP, GTRANS_TMP, EQUIV_ATOMS, PRIMITIVE_LATTICE, &
           & SPIN_FLIPS, NSYM_TMP, LATTICE, POSION, ATOM_TYPES, TENSORS, &
           & TENSOR_RANK, NATOMS, IS_MAGNETIC, TINY)
      AP = TRANSPOSE(PRIMITIVE_LATTICE)

      DEALLOCATE(TENSORS)

      RETURN
    END SUBROUTINE SPGLIB_MAGNETIC


    SUBROUTINE SPGLIB_TO_VASP(ISYMOP, GTRANS, PTRANS, IGRPOP, INVMAP, MAGROT, &
           NROT, NROTK, NPCELL, A, AP, ISYMOP_TMP, GTRANS_TMP, SPIN_FLIPS, &
           NSYM_TMP, NATOMS, NIOND, ISPIN, IU6)

      USE sym_prec
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NIOND, NATOMS, NSYM_TMP, ISPIN, IU6
      REAL(q), INTENT(IN) :: A(3, 3), AP(3, 3)  ! Basis vectors are given as column vectors.
      INTEGER, INTENT(IN) :: ISYMOP_TMP(3, 3, NSYM_TMP)
      REAL(q), INTENT(IN) :: GTRANS_TMP(3, NSYM_TMP)
      INTEGER, INTENT(IN) :: SPIN_FLIPS(NSYM_TMP)
      REAL(q), INTENT(OUT) :: PTRANS(NIOND + 2, 3)
      INTEGER, INTENT(OUT) :: NROT, NROTK, NPCELL
      INTEGER, INTENT(OUT) :: ISYMOP(3, 3, 48), IGRPOP(3, 3, 48)
      INTEGER, INTENT(OUT) :: INVMAP(48)
      REAL(q), INTENT(OUT) :: GTRANS(3, 48)
      REAL(q), INTENT(OUT) :: MAGROT(48, NIOND)

      REAL(q) :: B1(3), B2(3), B3(3)

      CALL SPGLIB_OPERATIONS_IN_PRIMITIVE(ISYMOP, GTRANS, PTRANS, IGRPOP, &
           & INVMAP, MAGROT, NROT, NROTK, NPCELL, A, AP, ISYMOP_TMP, &
           & GTRANS_TMP, SPIN_FLIPS, NSYM_TMP, NATOMS, NIOND, ISPIN, IU6)

      CALL RECIPS(1._q, A(1,1), A(1,2), A(1,3), B1, B2, B3)
      CALL SGRCON(ISYMOP, IGRPOP, NROTK, A(1,1), A(1,2), A(1,3), B1, B2, B3)
      CALL INVGRP(IGRPOP, NROTK, INVMAP)

      RETURN
    END SUBROUTINE SPGLIB_TO_VASP


    SUBROUTINE SPGLIB_OPERATIONS_IN_PRIMITIVE(ISYMOP, GTRANS, PTRANS, IGRPOP, &
           INVMAP, MAGROT, NROT, NROTK, NPCELL, A, AP, ISYMOP_TMP, GTRANS_TMP, &
           SPIN_FLIPS, NSYM_TMP, NATOMS, NIOND, ISPIN, IU6)

      USE sym_prec
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NIOND, NATOMS, NSYM_TMP, ISPIN, IU6
      REAL(q), INTENT(IN) :: A(3, 3), AP(3, 3)  ! Basis vectors are given as column vectors.
      INTEGER, INTENT(IN) :: ISYMOP_TMP(3, 3, NSYM_TMP)
      REAL(q), INTENT(IN) :: GTRANS_TMP(3, NSYM_TMP)
      INTEGER, INTENT(IN) :: SPIN_FLIPS(NSYM_TMP)
      REAL(q), INTENT(OUT) :: PTRANS(NIOND + 2, 3)
      INTEGER, INTENT(OUT) :: NROT, NROTK, NPCELL
      INTEGER, INTENT(OUT) :: ISYMOP(3, 3, 48), IGRPOP(3, 3, 48)
      INTEGER, INTENT(OUT) :: INVMAP(48)
      REAL(q), INTENT(OUT) :: GTRANS(3, 48)
      REAL(q), INTENT(OUT) :: MAGROT(48, NIOND)

      INTEGER :: I, J, K, L
      INTEGER :: SYMIND_MAP(NSYM_TMP)
      INTEGER :: ITMP(3,3)
      LOGICAL :: IN_PRIM, IS_FOUND
      REAL(q) :: TOLERANCE
      REAL(q) :: DIFF(3), GTMP(3), RTMP
      REAL(q) :: BP(3, 3), TRANS_MAT(3, 3), INV_TRANS_MAT(3, 3)
      REAL(q) :: GTRANS_PRIM(3, NSYM_TMP)

      ! Count number of identity matrices in ISYMOPs that spglib returned.
      ! This number corresponds to NPCELL.
      ! If cross check is expected, NPCELL * NROTK
      NPCELL = 0
      DO I = 1, NSYM_TMP
         IF ((ISYMOP_TMP(1, 1, I) == 1) .AND. &
              & (ISYMOP_TMP(1, 2, I) == 0) .AND. &
              & (ISYMOP_TMP(1, 3, I) == 0) .AND. &
              & (ISYMOP_TMP(2, 1, I) == 0) .AND. &
              & (ISYMOP_TMP(2, 2, I) == 1) .AND. &
              & (ISYMOP_TMP(2, 3, I) == 0) .AND. &
              & (ISYMOP_TMP(3, 1, I) == 0) .AND. &
              & (ISYMOP_TMP(3, 2, I) == 0) .AND. &
              & (ISYMOP_TMP(3, 3, I) == 1)) THEN
            NPCELL = NPCELL + 1
            PTRANS(NPCELL, :) = GTRANS_TMP(:, I)
         ENDIF
      ENDDO

      ! Spglib returns all symmetry operations in supercell.
      ! To match what VASP requires, unique ISYMOPs are collected.
      NROTK = 0
      DO I = 1, NSYM_TMP
         IN_PRIM = .TRUE.
         DO J = 1, NROTK
            IF ((ISYMOP_TMP(1, 1, I) == ISYMOP(1, 1, J)) .AND. &
                 & (ISYMOP_TMP(1, 2, I) == ISYMOP(1, 2, J)) .AND. &
                 & (ISYMOP_TMP(1, 3, I) == ISYMOP(1, 3, J)) .AND. &
                 & (ISYMOP_TMP(2, 1, I) == ISYMOP(2, 1, J)) .AND. &
                 & (ISYMOP_TMP(2, 2, I) == ISYMOP(2, 2, J)) .AND. &
                 & (ISYMOP_TMP(2, 3, I) == ISYMOP(2, 3, J)) .AND. &
                 & (ISYMOP_TMP(3, 1, I) == ISYMOP(3, 1, J)) .AND. &
                 & (ISYMOP_TMP(3, 2, I) == ISYMOP(3, 2, J)) .AND. &
                 & (ISYMOP_TMP(3, 3, I) == ISYMOP(3, 3, J))) THEN
               IN_PRIM = .FALSE.
            ENDIF
         ENDDO
         IF (IN_PRIM) THEN
            NROTK = NROTK + 1
            ISYMOP(:, :, NROTK) = ISYMOP_TMP(:, :, I)
         ENDIF
      ENDDO

      ! Mapping all symmetry operations to those of unique rotation parts
      DO J = 1, NROTK
         DO I = 1, NSYM_TMP
            IF ((ISYMOP_TMP(1, 1, I) == ISYMOP(1, 1, J)) .AND. &
                 & (ISYMOP_TMP(1, 2, I) == ISYMOP(1, 2, J)) .AND. &
                 & (ISYMOP_TMP(1, 3, I) == ISYMOP(1, 3, J)) .AND. &
                 & (ISYMOP_TMP(2, 1, I) == ISYMOP(2, 1, J)) .AND. &
                 & (ISYMOP_TMP(2, 2, I) == ISYMOP(2, 2, J)) .AND. &
                 & (ISYMOP_TMP(2, 3, I) == ISYMOP(2, 3, J)) .AND. &
                 & (ISYMOP_TMP(3, 1, I) == ISYMOP(3, 1, J)) .AND. &
                 & (ISYMOP_TMP(3, 2, I) == ISYMOP(3, 2, J)) .AND. &
                 & (ISYMOP_TMP(3, 3, I) == ISYMOP(3, 3, J))) THEN
               SYMIND_MAP(I) = J
               EXIT
            ENDIF
         ENDDO
      ENDDO

      CALL RECIPS(1._q, AP(1,1), AP(1,2), AP(1,3), BP(1, 1), BP(1, 2), BP(1, 3))
      TRANS_MAT = MATMUL(TRANSPOSE(BP), A)
      CALL RECIPS(1._q, TRANS_MAT(1,1), TRANS_MAT(1,2), TRANS_MAT(1,3), &
           & INV_TRANS_MAT(1,1), INV_TRANS_MAT(1,2), INV_TRANS_MAT(1,3))
      INV_TRANS_MAT = TRANSPOSE(INV_TRANS_MAT)

      ! Transform vector part with respect to primitive basis vectors
      GTRANS_PRIM = MATMUL(TRANS_MAT, GTRANS_TMP)

      DO J = 1, NROTK
         DO I = 1, NSYM_TMP
            IF (SYMIND_MAP(I) == J) THEN
               GTRANS(:, J) = GTRANS_PRIM(:, I)
               EXIT
            ENDIF
         ENDDO
      ENDDO

      DO I = 1, NROTK
         DO J = 1, 3
            GTRANS(J, I) = GTRANS(J, I) - NINT(GTRANS(J, I))
            IF (GTRANS(J, I) < -TINY) THEN
               GTRANS(J, I) = GTRANS(J, I) + 1
            ENDIF
         ENDDO
      ENDDO

      ! Transform back vector part with respect to original basis vectors
      GTRANS(:, :NROTK) = MATMUL(INV_TRANS_MAT, GTRANS(:, :NROTK))

      ! NROT depends on origin of input crystal structure.
      ! So this is not very useful parameters.
      ! This below is just to reproduce that VASP does.
      NROT = 0
      DO I = 1, NROTK
         IF ((ABS(GTRANS(1, I)) < TINY) .AND. &
              & (ABS(GTRANS(2, I)) < TINY) .AND. &
              & (ABS(GTRANS(3, I)) < TINY)) THEN
            NROT = NROT + 1
         ENDIF
      ENDDO

      ! Set MAGROT for collinear magnetism
      IF (ISPIN == 2) THEN
         MAGROT(:, :) = 0
         DO I = 1, NROTK
            DO J = 1, NPCELL
               IS_FOUND = .FALSE.

               ! Loop over K is to loosen tolerance to find DIFF(:) = 0 for
               ! slightly distorted case.
               TOLERANCE = TINY
               DO K = 1, 100
                  DO L = 1, NSYM_TMP
                     IF (SYMIND_MAP(K) /= I) THEN
                        CYCLE
                     ENDIF
                     DIFF = GTRANS(:, I) + PTRANS(J, :) - GTRANS_TMP(:, L)
                     DIFF = DIFF - NINT(DIFF)
                     DIFF = matmul(A, DIFF)
                     IF ((ABS(DIFF(1)) < TOLERANCE) .AND. &
                          & (ABS(DIFF(2)) < TOLERANCE) .AND. &
                          & (ABS(DIFF(3)) < TOLERANCE)) THEN
                        MAGROT(I, J) = SPIN_FLIPS(L)
                        IS_FOUND = .TRUE.
                        EXIT
                     ENDIF
                     ! When not found, loosen tolerance.
                  ENDDO
                  IF (IS_FOUND) THEN
                     EXIT
                  ENDIF
                  TOLERANCE = TOLERANCE * 1.05
               ENDDO
            ENDDO
         ENDDO

         ! A trivial consistency check of MAGROT
         DO I = 1, NROTK
            DO J = 1, NPCELL
               IF (ABS(MAGROT(I, J)) < TINY) THEN
                  WRITE(IU6, *) 'warning: MAGROT is wrongly assigned.'
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      ! sort the symmetry operations, bring the NROT pure rotation
      ! operations to the front
      DO I = 1, NROT
         ! if this operation is a pure rotation operation cycle I
         IF (ALL(ABS(GTRANS(:, I)) < TINY)) CYCLE
         ! if not then swap it against an operation that is purely rotational
         IS_FOUND=.FALSE.
         DO J = NROT + 1, NROTK
            ! if operation J is a pure rotation operation, swap (I <--> J)
            IF (ALL(ABS(GTRANS(:, J)) < TINY)) THEN
               ! swap ISYMOP
               ITMP = ISYMOP(:, :, I); ISYMOP(:, :, I) = ISYMOP(:, :, J); ISYMOP(:, :, J) = ITMP
               ! swap GTRANS
               GTMP = GTRANS(:, I); GTRANS(:, I) = GTRANS(:, J); GTRANS(:, J) = GTMP
               ! and swap MAGROT
               IF (ISPIN==2) THEN
               DO K = 1, NPCELL
                  RTMP = MAGROT(I, K); MAGROT(I, K) = MAGROT(J, K); MAGROT(J, K) = RTMP
               ENDDO
               ENDIF
               ! and move to the next I
               IS_FOUND=.TRUE.; EXIT
            ENDIF
         ENDDO
         IF (.NOT.IS_FOUND) THEN
            ! this point should not be reached ...
            CALL vtutor%bug('error: there do not seem to be enough pure rotation operators.', __FILE__, __LINE__)
         ENDIF
      ENDDO

      RETURN
    END SUBROUTINE SPGLIB_OPERATIONS_IN_PRIMITIVE


    SUBROUTINE SPGLIB_GET_PRIMITIVE(NATOMS, NIOND, NROTK, NPCELL, ROTMAP, POSION, A, AP, IU6)

      USE sym_prec
      USE iso_c_binding
      USE spglib_f08
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NATOMS, NIOND, NROTK, NPCELL, IU6
      REAL(q), INTENT(IN) :: POSION(3, NIOND)
      REAL(q), INTENT(IN) :: A(3, 3)

      REAL(q), INTENT(OUT) :: AP(3,3)

      INTEGER :: ROTMAP(NIOND, NROTK, NPCELL)

      ! local variables
      REAL(q) :: LATTICE(3, 3)
      REAL(q) :: POSION_TMP(3, NATOMS)
      INTEGER :: ATOM_TYPES(NATOMS)
      INTEGER :: IRETURN

      TYPE(SpglibDataset) :: DATASET
      TYPE(SpglibSpacegroupType) :: SPGTYPE

      REAL(q) :: OMEGA, OMEGAP

      INTEGER :: I, J, K, ITYPE, IP

      LATTICE = TRANSPOSE(A)

      POSION_TMP(1:3, 1:NATOMS) = POSION(1:3, 1:NATOMS)

      ITYPE = 0 ; ATOM_TYPES = 0

      atoms: DO I = 1, NATOMS
         IF (ATOM_TYPES(I)/=0) CYCLE atoms

         ITYPE = ITYPE + 1
         ATOM_TYPES(I) = ITYPE

         DO J = 1, NROTK
         DO K = 1, NPCELL
            IP = ROTMAP( I, J, K )
            IF (ATOM_TYPES(IP) /= 0) THEN
              IF (ATOM_TYPES(IP) /= ITYPE) THEN
                 CALL vtutor%bug('SPGLIB_GET_PRIMITIVE: internal error: POSMAP corrupted', __FILE__, __LINE__)
              ENDIF
            ELSE
               ATOM_TYPES(IP) = ITYPE
            ENDIF
         ENDDO
         ENDDO

      ENDDO atoms

      IF (ANY(ATOM_TYPES == 0)) THEN
         CALL vtutor%bug('SPGLIB_GET_PRIMITIVE: internal error: failed to set up ATOM_TYPES', __FILE__, __LINE__)
      ENDIF

      DATASET = spg_get_dataset( LATTICE, POSION_TMP, ATOM_TYPES, NATOMS, TINY)

      AP=TRANSPOSE( Dataset % primitive_lattice)

!!      IRETURN = spg_standardize_cell( LATTICE, POSION_TMP, ATOM_TYPES, NATOMS, 1, 0, TINY)
!!      AP=TRANSPOSE(LATTICE)

      CALL CELVOL(A (1,1), A (1,2), A (1,3), OMEGA )
      CALL CELVOL(AP(1,1), AP(1,2), AP(1,3), OMEGAP)

      IF (NINT(OMEGA/OMEGAP)/=NPCELL) THEN
         CALL vtutor%bug('SPGLIB_GET_PRIMITIVE: internal error: volume of primitive cell is wrong' &
            // str(OMEGA/OMEGAP) // " " // str(NPCELL), __FILE__, __LINE__)
      ENDIF

      IF (IU6>0) THEN
         WRITE(IU6,'(//X,A,/X,A)') &
           'SPGLIB: symmetry analysis of effective configuration:', &
           '------------------------------------------------------'
         WRITE(IU6,'(A,X,A17,6X,A,I4)') &
           ' International :', TRIM( DATASET % international_symbol ), 'Nr.', DATASET % spacegroup_number
         WRITE(IU6,'(A,X,A17,6X,A,I4)') &
           ' Hall          :', TRIM( DATASET % hall_symbol ),          'Nr.', DATASET % hall_number

         SPGTYPE = spg_get_spacegroup_type( DATASET % hall_number )

         WRITE(IU6,'(A,X,A17)') &
           ' Schoenflies   :', TRIM( SPGTYPE % Schoenflies )

         WRITE(IU6,'(/A,I3)') &
           ' Nr. of symops.:', DATASET % n_operations

         WRITE(IU6,'(/A/,3(/X,A,3(F15.10,A)))') &
           ' Primitive cell:', &
           'A1 = (', AP(1,1), ',', AP(2,1), ',', AP(3,1), ')', &
           'A2 = (', AP(1,2), ',', AP(2,2), ',', AP(3,2), ')', &
           'A3 = (', AP(1,3), ',', AP(2,3), ',', AP(3,3), ')'
      ENDIF

      RETURN
    END SUBROUTINE SPGLIB_GET_PRIMITIVE


    SUBROUTINE SPGLIB_GET_POINTGROUP(ISYMOP, NSYM, IPGIND, HermannMaguinSymbol, IU6)
      USE iso_c_binding
      USE spglib_f08
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: ISYMOP(3, 3, NSYM), NSYM, IU6
      INTEGER, INTENT(OUT) :: IPGIND
      CHARACTER(LEN=6), INTENT(OUT) :: HermannMaguinSymbol
      ! local variables
      INTEGER :: I,ITMP(3,3)
      CHARACTER(KIND=C_CHAR) :: SYMBOL6(6)

      IPGIND = spg_get_pointgroup( SYMBOL6, ITMP, ISYMOP, NSYM)

      DO I = 1, SIZE(SYMBOL6)
         IF (SYMBOL6(I) == C_NULL_CHAR) THEN
            HermannMaguinSymbol(I:) = ' '
            EXIT
         ENDIF
         HermannMaguinSymbol(I:I) = SYMBOL6(I)
      ENDDO

      IF (IU6>0) WRITE(IU6,'(/X,A,2X,A)') 'SPGLIB_GET_POINTGROUP: associated pointgroup (Hermann-Maguin): ',HermannMaguinSymbol

      RETURN
    END SUBROUTINE SPGLIB_GET_POINTGROUP


    FUNCTION IPOINTGROUP_SPGLIB2VASP(I)
      INTEGER :: IND(32)=(/ 1, 2, 3, 4, 5, 6, 7, 8, &
        14, 15, 16, 17, 18, 19, 20,   9, 10, 11, 12, 13, &
        21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32 /)
      INTEGER :: I
      IF (I>0 .AND. I<33) THEN
         IPOINTGROUP_SPGLIB2VASP=IND(I)
      ELSE
         IPOINTGROUP_SPGLIB2VASP=0
      ENDIF
    END FUNCTION IPOINTGROUP_SPGLIB2VASP
#endif
#endif

!************************* SUBROUTINE WRTSYM ***************************
!
!***********************************************************************

      SUBROUTINE WRTSYM(NATOMS,NIOND,PTRANS,ROTMAP,MAGROT,ISPIN,IU6)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)

      INTEGER, POINTER :: ROTMAP(:,:,:)
      REAL(q), POINTER :: MAGROT(:,:)
      DIMENSION PTRANS(NIOND+2,3),INDROT(NIOND+2)
      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

      IF (IU6>0) THEN
         DO I=1,NROTK
            WRITE(IU6,'(//X,A,I4)') 'irot  :',I
            WRITE(IU6,'(X,A)') '--------------------------------------------------------------------'
            WRITE(IU6,'(X,A,3I4/,8X,3I4/,8X,3I4)') 'isymop:',ISYMOP(:,:,I)
            WRITE(IU6,'(/X,A,3F14.7)') 'gtrans:',GTRANS(:,I)
            DO J=1,NPCELL
!              WRITE(IU6,'(/A,I4)') 'NP=',J
               WRITE(IU6,'(/X,A,3F14.7)') 'ptrans:',PTRANS(J,:)
               IF (ISPIN==2) WRITE(IU6,'(/X,A,F14.7)') 'magrot:',MAGROT(I,J)
               WRITE(IU6,*)
               DO IA=1,NATOMS
                  INDROT(ROTMAP(IA,I,J))=IA
               ENDDO
               WRITE(IU6,'(X,A)') 'rotmap:'
               DO IA=1,NATOMS,5
!                 WRITE(IU6,'(A)',ADVANCE='No') "("
                  DO K=IA,MIN(IA+4,NATOMS)
!                    WRITE(IU6,'("(",I4,"<-",I4,")",2X)',ADVANCE='No')  K,ROTMAP(K,I,J)
                     WRITE(IU6,'(X,"(",I4,"->",I4,")",X)',ADVANCE='No')  K,INDROT(K)
                  ENDDO
                  WRITE(IU6,*)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE WRTSYM

!************************* SUBROUTINE NOSYMM ***************************
!
! subroutine nosymm should be called when no usage of any symmetry
! will be desired (initialise all stuff for point group C_1 ...)
! it sets up: the rotation matrices in direct space (ISYMOP)
!             the rotation matrices in reciprocal space (IGRPOP)
!             and a table indexing their inverse elements (INVMAP)
!             the non-trivial translations for each rotation (GTRANS)
!             a table (ROTMAP) which shows the connections between atoms
!             the number of primitive cells within a supercell (NPCELL)
!             and their coordinates (PTRANS) within the supercell
!     all will be initialised as if we had no symmetry (1 symmetry
!     operation = unity operator in both spaces, no non-trivial
!     translations, cell is a primitive cell ...)
! the routine needs the number of atomic species (NTYP)
!             the number of atoms of each species (NITYP)
!             the dimensioning parameter NIOND (max. no. of atoms)
!             and the unit number where to write informations (IU6)
!
!***********************************************************************

      SUBROUTINE NOSYMM(A,NTYP,NITYP,NIOND,PTRANS,NRTK,NPCLL,ROTMAP,MAGROT,ISPIN,IU6)
      USE prec


      IMPLICIT REAL(q) (A-H,O-Z)
      REAL(q) MAGROT(48,NIOND)
      INTEGER, POINTER :: ROTMAP(:,:,:)
      DIMENSION PTRANS(NIOND+2,3),NITYP(NTYP),A(3,3)

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL
      IF (IU6>0) THEN
      WRITE(IU6,'(/,2A,/2A,/)') 'IMPORTANT INFORMATION: All ', &
     &                  'symmetrisations will be switched off!', &
     &                  'NOSYMM: (Re-)initialisation of all symmetry ', &
     &                                      'stuff for point group C_1.'
      ENDIF
      NROT=1
      NROTK=1
      ISYMOP(1,1,1)=1
      ISYMOP(2,2,1)=1
      ISYMOP(3,3,1)=1
      ISYMOP(1,2,1)=0
      ISYMOP(2,1,1)=0
      ISYMOP(1,3,1)=0
      ISYMOP(3,1,1)=0
      ISYMOP(2,3,1)=0
      ISYMOP(3,2,1)=0
      IGRPOP(1,1,1)=1
      IGRPOP(2,2,1)=1
      IGRPOP(3,3,1)=1
      IGRPOP(1,2,1)=0
      IGRPOP(2,1,1)=0
      IGRPOP(1,3,1)=0
      IGRPOP(3,1,1)=0
      IGRPOP(2,3,1)=0
      IGRPOP(3,2,1)=0
      GTRANS(1,1)=0._q
      GTRANS(2,1)=0._q
      GTRANS(3,1)=0._q
      INVMAP(1)=1
      NATOMS=0
      DO I=1,NTYP
         NATOMS=NATOMS+NITYP(I)
      ENDDO
      DO I=1,NATOMS
         ROTMAP(I,1,1)=I
      ENDDO
      IF (ISPIN==2) THEN
         DO I=1,NATOMS
            MAGROT(1,I)=1._q
         ENDDO
      ENDIF
      NPCELL=1
      PTRANS(1,1)=0._q
      PTRANS(1,2)=0._q
      PTRANS(1,3)=0._q
      AP(1,1)=A(1,1)
      AP(2,1)=A(2,1)
      AP(3,1)=A(3,1)
      AP(1,2)=A(1,2)
      AP(2,2)=A(2,2)
      AP(3,2)=A(3,2)
      AP(1,3)=A(1,3)
      AP(2,3)=A(2,3)
      AP(3,3)=A(3,3)
      NRTK=NROTK
      NPCLL=NPCELL
      RETURN
      END SUBROUTINE NOSYMM


!******************** SUBROUTINE CHKTOL ********************************
!
! This routine calls LATTYP and INISYM with three different values for
! the symmetry tolerance parameter: TINY*SCALE, TINY, and TINY/SCALE
! When these calls detect different Bravais lattice types and/or a
! different number of symmetry operations a warning is printed.
!
!***********************************************************************

      SUBROUTINE CHKTOL(A,POSION,VEL,LSFOR,LSDYN,LPEAD,NTYP,NITYP,NIOND, &
     &           PTRANS,NRTK,NPCLL,ROTMAP,TAU,TAUROT,WRKROT,INDROT,ATOMOM,SAXIS,MAGROT,ISPIN,IU6,NWRITE)
      USE prec
      USE sym_prec
      USE tutor, ONLY : vtutor,isAlert,SymprecDependence

      IMPLICIT REAL(q) (A-H,O-Z)

      LOGICAL LSDYN,LSFOR(3,NIOND),LPEAD
      REAL(q) MAGROT(48,NIOND),SAXIS(3),ALPHA,BETA,ATOMOM_(3)
      INTEGER, POINTER ::  ROTMAP(:,:,:)

      DIMENSION A(3,3),POSION(3,NIOND),VEL(3,NIOND),NITYP(NTYP)
      DIMENSION TAU(NIOND,3),TAUROT(NIOND,3),WRKROT(3*(NIOND+2))
      DIMENSION PTRANS(NIOND+2,3),INDROT(NIOND+2),ATOMOM(3*NIOND)

      DIMENSION ATMP(3,3,3),IBRAV(3),CELDIM(6,3),NTMP(3,3)

      REAL(q), PARAMETER :: SCALE=10._q     

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

      TINY_SAVE=TINY

      TINY=TINY_SAVE*SCALE
      ATMP(:,:,1)=A
      CALL LATTYP(ATMP(1,1,1),ATMP(1,2,1),ATMP(1,3,1),IBRAV(1),CELDIM(1,1),-1)

      CALL INISYM(A,POSION,VEL,LSFOR,LSDYN,LPEAD,NTYP,NITYP,NIOND, &
     &   PTRANS,NRTK,NPCLL,ROTMAP,TAU,TAUROT,WRKROT,INDROT,ATOMOM,SAXIS,MAGROT,ISPIN,-1)
      NTMP(1,1)=NROTK; NTMP(2,1)=NROT; NTMP(3,1)=NPCELL

      TINY=TINY_SAVE
      ATMP(:,:,2)=A
      CALL LATTYP(ATMP(1,1,2),ATMP(1,2,2),ATMP(1,3,2),IBRAV(2),CELDIM(1,2),-1)

      CALL INISYM(A,POSION,VEL,LSFOR,LSDYN,LPEAD,NTYP,NITYP,NIOND, &
     &   PTRANS,NRTK,NPCLL,ROTMAP,TAU,TAUROT,WRKROT,INDROT,ATOMOM,SAXIS,MAGROT,ISPIN,-1)
      NTMP(1,2)=NROTK; NTMP(2,2)=NROT; NTMP(3,2)=NPCELL

      TINY=TINY_SAVE/SCALE
      ATMP(:,:,3)=A
      CALL LATTYP(ATMP(1,1,3),ATMP(1,2,3),ATMP(1,3,3),IBRAV(3),CELDIM(1,3),-1)

      CALL INISYM(A,POSION,VEL,LSFOR,LSDYN,LPEAD,NTYP,NITYP,NIOND, &
     &   PTRANS,NRTK,NPCLL,ROTMAP,TAU,TAUROT,WRKROT,INDROT,ATOMOM,SAXIS,MAGROT,ISPIN,-1)
      NTMP(1,3)=NROTK; NTMP(2,3)=NROT; NTMP(3,3)=NPCELL

      IF (((IBRAV(1)/=IBRAV(2)).OR.(IBRAV(2)/=IBRAV(3))).OR. &
          (ANY(NTMP(:,1)/=NTMP(:,2)).OR.ANY(NTMP(:,2)/=NTMP(:,3)))) THEN

         IF (IU6>=0.AND.NWRITE==3) THEN

            WRITE(IU6,'(/2X,A,E9.2)') 'SYMPREC =',TINY_SAVE*SCALE
            WRITE(IU6,'(2X,A)')       '------------------'
            CALL PRIBRAV(IBRAV(1),CELDIM(1,1),ATMP(1,1,1),ATMP(1,2,1),ATMP(1,3,1),IU6)
            WRITE(IU6,50) NTMP(:,1)

            WRITE(IU6,'( 2X,A,E9.2)') 'SYMPREC =',TINY_SAVE
            WRITE(IU6,'(2X,A)')       '------------------'
            CALL PRIBRAV(IBRAV(2),CELDIM(1,2),ATMP(1,1,2),ATMP(1,2,2),ATMP(1,3,2),IU6)
            WRITE(IU6,50) NTMP(:,2)

            WRITE(IU6,'( 2X,A,E9.2)') 'SYMPREC =',TINY_SAVE/SCALE
            WRITE(IU6,'(2X,A)')       '------------------'
            CALL PRIBRAV(IBRAV(3),CELDIM(1,3),ATMP(1,1,3),ATMP(1,2,3),ATMP(1,3,3),IU6)
            WRITE(IU6,50) NTMP(:,3)

         ENDIF

         CALL VTUTOR%WRITE(isAlert,SymprecDependence)

      ENDIF

      TINY=TINY_SAVE

   50 FORMAT(/' INISYM: Found ',I2,' space group operations (whereof ',I2,' pure point group operations),', &
     &       /'         and ',I4,' ''primitive'' translations'/)

      RETURN
      END SUBROUTINE CHKTOL


!************************* SUBROUTINE STRCNSRV *************************
!
!***********************************************************************

      SUBROUTINE STRCNSRV(A,ISYMOP,IGRPOP,GTRANS,NROTK,NROT)

      USE prec
      USE sym_prec

      IMPLICIT NONE

      INTEGER :: ISYMOP(3,3,48),IGRPOP(3,3,48)
      REAL(q) :: A(3,3),GTRANS(3,48)
      INTEGER :: NROTK,NROT

      INTEGER I,J,K,IS
      INTEGER INOT,INROT,ISORNOT,LIST(48)
      REAL(q) QST(3,6),RQST(3,6)
      INTEGER ISYMOPTMP(3,3,48),IGRPOPTMP(3,3,48)
      REAL(q) GTRANSTMP(3,48)
      LOGICAL ISMAPPED(6)

      REAL(q) V1,V2,V3

      INOT=0; INROT=0; LIST=1 

      ! loop on sym. op. 
      isym: DO IS=1,NROTK

         QST=0._q
         QST(1,1)= 1._q; QST(2,2)= 1._q; QST(3,3)= 1._q
         QST(1,4)=-1._q; QST(2,5)=-1._q; QST(3,6)=-1._q

         ! cartesian to direct 
         DO I=1,6 ! 3 vectors
            V1=QST(1,I)*A(1,1)+QST(2,I)*A(2,1)+QST(3,I)*A(3,1)
            V2=QST(1,I)*A(1,2)+QST(2,I)*A(2,2)+QST(3,I)*A(3,2)
            V3=QST(1,I)*A(1,3)+QST(2,I)*A(2,3)+QST(3,I)*A(3,3)
            QST(1,I)=V1
            QST(2,I)=V2
            QST(3,I)=V3
         ENDDO

         ! rotate the q-star
         RQST=0._q
         DO K=1,6 ! x y z direction of IQST
            DO I=1,3 ! component of IRQST
               DO J=1,3 ! summation index (row column) 
                  RQST(I,K)=RQST(I,K)+IGRPOP(J,I,IS)*QST(J,K)     
               ENDDO
            ENDDO
         ENDDO

         ISORNOT = 0
         idir: DO K=1,6 ! x y z direction for QST
            ISMAPPED=.FALSE.
            DO I=1,6 ! x y z direction for RQST 
               IF (ALL(ABS(QST(:,K)-RQST(:,I))<TINY)) ISMAPPED(I)=.TRUE.
#ifdef debug
               WRITE(*,'(A,3F9.5,A,3F9.5,A,L5)') 'original ',QST(:,K),' rotated ',RQST(:,I),' is it mapped ? ',ISMAPPED(I)
#endif
            ENDDO
            IF (ALL(.NOT.ISMAPPED)) ISORNOT=ISORNOT+1 
#ifdef debug
            WRITE(*,*) IS , K,ISORNOT,ISMAPPED
#endif
         ENDDO idir
         
         IF (ISORNOT>0) THEN
#ifdef debug
            IF (IU6>0) WRITE(*,'(A,I3,A)') 'SYM. OP.',IS,' is NOT SUITABLE FOR LR_NMR and set to zero'
#endif
            INOT=INOT+1
            ISYMOP(:,:,IS)=0
            IGRPOP(:,:,IS)=0
            GTRANS(:,IS)=0
            LIST(IS)=0
         ENDIF
      ENDDO isym
#ifdef debug
      IF (IU6>0) WRITE(*,*) 'NUMBER OF SYMMETRY OPERATIONS NOT SUITABLE',INOT
#endif
      ! reorganize the arrays
      K=0
      DO IS=1,NROTK
         IF (LIST(IS).EQ.1) THEN
            K=K+1
            ISYMOPTMP(:,:,K)=ISYMOP(:,:,IS)
            IGRPOPTMP(:,:,K)=IGRPOP(:,:,IS)
            GTRANSTMP(:,K)   =GTRANS(:,IS)
            IF ((GTRANS(1,IS)==0).AND.(GTRANS(2,IS)==0).AND.(GTRANS(3,IS)==0)) INROT=INROT+1
         ENDIF
      ENDDO

      ISYMOP(:,:,:)=ISYMOPTMP(:,:,:)
      IGRPOP(:,:,:)=IGRPOPTMP(:,:,:)
      GTRANS(:,:)  =GTRANSTMP(:,:)
      NROTK=NROTK-INOT
      NROT=INROT 

      RETURN
      END SUBROUTINE STRCNSRV


!************************* SUBROUTINE QSTAR_SYM_REDUCE *****************
!
!***********************************************************************
!!#define debug
      SUBROUTINE QSTAR_SYM_REDUCE(A,T_INFO,SYMM,IU6)
      USE prec
      USE base
      USE lattice
      USE poscar
      USE sym_prec
      IMPLICIT NONE
      TYPE(type_info) T_INFO
      TYPE(symmetry) SYMM
      REAL(q) A(3,3)
      INTEGER IU6
      ! local variables
      INTEGER I,J,K,IS
      INTEGER INOT,INROT,ISORNOT,LIST(48)
      REAL(q) QST(3,6),RQST(3,6)
      REAL(q) V1,V2,V3
      REAL(q) MODQST,MODRQST,PRODMOD,DOT
      INTEGER ISYMOPTMP(3,3,48),IGRPOPTMP(3,3,48)
      REAL(q) GTRANSTMP(3,48)
      LOGICAL ISMAPPED(6)
      INTEGER IPGIND
      CHARACTER(4) GRPNAM
      REAL(q) TAU(T_INFO%NIONS,3),TAUROT(T_INFO%NIONS,3)
      ! symmetry common block
      INTEGER  ISYMOP,NROT,IGRPOP,NROTK,INVMAP,NPCELL
      REAL(q)  GTRANS,AP
      COMMON /SYMM/   ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
           GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

#ifdef debug
      IF (IU6>0) THEN
      WRITE(*,'(A)') 'OLD ISYMOP'
      DO IS=1,NROTK 
         WRITE(*,'(A,I5)') 'SYM. OP',IS
         DO J=1,3
            WRITE(*,'(7I5)') (ISYMOP(I,J,IS),I=1,3), (IGRPOP(I,J,IS),I=1,3),INVMAP(IS)
         ENDDO
         WRITE(*,*) 'OLD GTRANS'
         WRITE(*,'(3F8.5)') (GTRANS(J,IS),J=1,3) 
      ENDDO  
      WRITE(*,*) 'NROT=',NROT,'NROTK=',NROTK
      ENDIF
#endif

      INOT=0; INROT=0; LIST=1 

      ! loop on sym. op. 
      isym: DO IS=1,NROTK
!        WRITE (*,'(A,I5)') 'SYMOP=',IS
         ! q-star
         QST=0._q
         QST(1,1)= 1._q; QST(2,2)= 1._q; QST(3,3)= 1._q
         QST(1,4)=-1._q; QST(2,5)=-1._q; QST(3,6)=-1._q

         ! cartesian to direct 
         DO I=1,6 ! 3 vectors
            V1=QST(1,I)*A(1,1)+QST(2,I)*A(2,1)+QST(3,I)*A(3,1)
            V2=QST(1,I)*A(1,2)+QST(2,I)*A(2,2)+QST(3,I)*A(3,2)
            V3=QST(1,I)*A(1,3)+QST(2,I)*A(2,3)+QST(3,I)*A(3,3)
            QST(1,I)=V1
            QST(2,I)=V2
            QST(3,I)=V3
         ENDDO

         ! rotate the q-star
         RQST=0._q
         DO K=1,6 ! x y z direction of IQST
            DO I=1,3 ! component of IRQST
               DO J=1,3 ! summation index (row column) 
                  RQST(I,K)=RQST(I,K)+IGRPOP(J,I,IS)*QST(J,K)     
               ENDDO
            ENDDO
         ENDDO

         ISORNOT = 0
         idir: DO K=1,6 ! x y z direction for QST
            ISMAPPED=.FALSE.
            DO I=1,6 ! x y z direction for RQST 
               IF (ALL(ABS(QST(:,K)-RQST(:,I))<TINY)) ISMAPPED(I)=.TRUE.
#ifdef debug
               WRITE(*,'(A,3F9.5,A,3F9.5,A,L5)') 'original ',QST(:,K),' rotated ',RQST(:,I),' is it mapped ? ',ISMAPPED(I)
#endif
            ENDDO
            IF (ALL(.NOT.ISMAPPED)) ISORNOT=ISORNOT+1 
#ifdef debug
            WRITE(*,*) IS , K,ISORNOT,ISMAPPED
#endif
         ENDDO idir
         
         IF (ISORNOT>0) THEN
#ifdef debug
            IF (IU6>0) WRITE(*,'(A,I3,A)') 'SYM. OP.',IS,' is NOT SUITABLE FOR LR_NMR and set to zero'
#endif
            INOT=INOT+1
            ISYMOP(:,:,IS)=0
            IGRPOP(:,:,IS)=0
            GTRANS(:,IS)=0
            LIST(IS)=0
         ENDIF
      ENDDO isym
#ifdef debug
      IF (IU6>0) WRITE(*,*) 'NUMBER OF SYMMETRY OPERATIONS NOT SUITABLE',INOT
#endif
      ! reorganize the arrays
      K=0
      DO IS=1,NROTK
         IF (LIST(IS).EQ.1) THEN
            K=K+1
            ISYMOPTMP(:,:,K)=ISYMOP(:,:,IS)
            IGRPOPTMP(:,:,K)=IGRPOP(:,:,IS)
            GTRANSTMP(:,K)   =GTRANS(:,IS)
            IF ((GTRANS(1,IS)==0).AND.(GTRANS(2,IS)==0).AND.(GTRANS(3,IS)==0)) INROT=INROT+1
         ENDIF
      ENDDO

      ISYMOP(:,:,:)=ISYMOPTMP(:,:,:)
      IGRPOP(:,:,:)=IGRPOPTMP(:,:,:)
      GTRANS(:,:)  =GTRANSTMP(:,:)
      NROTK=NROTK-INOT
      NROT=INROT 

!=======================================================================
! Now try to find out which symmetry operations we have here:
!=======================================================================
      CALL PGROUP(ISYMOP,NROT,IPGIND,GRPNAM)
      IF (IU6>0) &
      WRITE(IU6,'(/4A)') '(QSTAR_SYM_REDUCE) The star conserving configuration has the ', &
     &                                      'point symmetry ',GRPNAM,'.'
      IF (NROTK>NROT) THEN
         CALL PGROUP(ISYMOP,NROTK,IPGIND,GRPNAM)
         IF (IU6>0) &
         WRITE(IU6,'(/4A)') '(QSTAR_SYM_REDUCE) The point group associated with ', &
     &                             'its full space group is ',GRPNAM,'.'
      ENDIF
      DO I=1,T_INFO%NIONS
         TAU(I,1)=T_INFO%POSION(1,I)
         TAU(I,2)=T_INFO%POSION(2,I)
         TAU(I,3)=T_INFO%POSION(3,I)
      ENDDO
      IF (ASSOCIATED(SYMM%ROTMAP)) THEN
         DEALLOCATE(SYMM%ROTMAP)
      ENDIF
      ALLOCATE(SYMM%ROTMAP(T_INFO%NIONS,NROTK,NPCELL))
      CALL POSMAP(TAU,ISYMOP,GTRANS,NROTK,SYMM%PTRANS,NPCELL,T_INFO%NIONS,1, &
     &            T_INFO%NTYP,T_INFO%NIONS,T_INFO%NITYP,SYMM%ROTMAP(1,1,1),TAUROT)

!=======================================================================
! map of the inverse elements of each (reciprocal space) group element
!======================================================================
      CALL INVGRP(IGRPOP,NROTK,INVMAP)

      SYMM%NROT=NROTK

#ifdef debug
      IF (IU6>0) THEN
      WRITE(*,'(A)') 'NEW ISYMOP'
      DO IS=1,NROTK 
         WRITE(*,'(A,I5)') 'SYM. OP',IS
         DO J=1,3
            WRITE(*,'(7I5)') (ISYMOP(I,J,IS),I=1,3),(IGRPOP(I,J,IS),I=1,3),INVMAP(IS)
         ENDDO
         WRITE(*,*) 'NEW GTRANS'
         WRITE(*,'(3F8.5)') (GTRANS(J,IS),J=1,3) 
      ENDDO
      ENDIF
#endif
      RETURN
      END SUBROUTINE QSTAR_SYM_REDUCE
!!#undef debug

END MODULE msymmetry


!******************* SUBROUTINE RHOSYM *********************************
! RCS:  $Id: symmetry.F,v 1.8 2003/06/27 13:22:23 kresse Exp kresse $
!
! subroutine rhosym is an interface to the charge density symmetrization
! routines RHOSYR and RHOSYG which updates CHTOT (the reciprocal space
! charge density)
!
!***********************************************************************

      SUBROUTINE RHOSYM(CHTOT,GRIDC,PTRANS,NIOND,MAGROT,ISP)
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d) GRIDC

      COMPLEX(q)   CHTOT(GRIDC%RC%NP)
      REAL(q)   MAGROT(48,NIOND)
      DIMENSION PTRANS(NIOND+2,3)
! work arrays
      COMPLEX(q),ALLOCATABLE ::  CWORK(:),CWORK2(:)
      REAL(q),ALLOCATABLE    ::  WORK(:)

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL


! Trivial case: no symmetry ... :
      IF (NROTK==1) RETURN

      N =GRIDC%NGX_rd*GRIDC%NGY*GRIDC%NGZ_rd
      NP=GRIDC%NGX*GRIDC%NGY*GRIDC%NGZ
!-----------------------------------------------------------------------
! Pure point symmetry! we can use RHOSYR (real space symmetrization):
!-----------------------------------------------------------------------
      IF (NROT==NROTK .AND. NPCELL==1 .AND. .FALSE.)  THEN
         ALLOCATE(CWORK(N),WORK(NP))
!$ACC ENTER DATA CREATE(CWORK) __IF_ASYNC__

         CALL FFT3D(CHTOT,GRIDC,1)
         ! merge charge to CWORK (result is real)
         CALL MRG_GRID_RL(GRIDC, CWORK, CHTOT)
         CALL RHOSYR(CWORK,GRIDC%NGX,GRIDC%NGX,GRIDC%NGY,GRIDC%NGY,GRIDC%NGZ,GRIDC%NGZ,MAGROT(1:48,1:1),ISP,WORK)
         ! and back the data go to CHTOT
         CALL DIS_GRID_RL(GRIDC, CWORK, CHTOT, .FALSE.)
! FFT of the symmetrized real space density (CHTOT):
         CALL FFT_RC_SCALE(CHTOT,CHTOT,GRIDC)
         CALL SETUNB_COMPAT(CHTOT,GRIDC)

!$ACC EXIT DATA DELETE(CWORK) __IF_ASYNC__
         DEALLOCATE(CWORK,WORK)
      ELSE

!-----------------------------------------------------------------------
! General case including nontrivial translations! We have to work in
! reciprocal space, so we must use routine RHOSYG ... :
!-----------------------------------------------------------------------
         ALLOCATE(CWORK(N),CWORK2(N),WORK(NP))
!$ACC ENTER DATA CREATE(CWORK,CWORK2) __IF_ASYNC__

         ! merge charge to CWORK
         CALL MRG_GRID_RC(GRIDC, CWORK, CHTOT)
         CALL RHOSYG(CWORK2,CWORK,GRIDC%NGX_rd,GRIDC%NGX,GRIDC%NGY,GRIDC%NGY,&
               GRIDC%NGZ_rd,GRIDC%NGZ,PTRANS,NIOND,MAGROT,ISP,WORK)
         ! and back the data go to CHTOT
         CALL DIS_GRID_RC(GRIDC, CWORK2, CHTOT, .FALSE.)

!$ACC EXIT DATA DELETE(CWORK,CWORK2) __IF_ASYNC__
         DEALLOCATE(CWORK,CWORK2,WORK)
      ENDIF

      RETURN
      END SUBROUTINE


!******************* SUBROUTINE SYMFIELD *******************************
!
! Subroutine SYMFIELD is an interface to the magnetization
! symmetrization routine MAGSYG
!
!***********************************************************************

      SUBROUTINE SYMFIELD(MAGFLD,GRIDC,PTRANS,NIOND,MAGROT,SAXIS,LATT_CUR)

#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE prec
      USE mpimy
      USE mgrid
      USE lattice

      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d) GRIDC
      TYPE (latt) LATT_CUR

      COMPLEX(q) MAGFLD(GRIDC%MPLWV,3)

      REAL(q)   MAGROT(48,NIOND),SAXIS(3)
      DIMENSION PTRANS(NIOND+2,3)
! work arrays
      COMPLEX(q),ALLOCATABLE ::  CWORK(:),CWORK2(:)
      REAL(q),ALLOCATABLE    ::  WORK(:)

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL
#if defined(NGXhalf) || defined(NGZhalf)
      CALL ERROR('SYMFIELD', 'only supports full grid mode', 3)
#endif
! Trivial case: no symmetry ... :
      IF (NROTK==1) RETURN

      N =3*GRIDC%NGX*GRIDC%NGY*GRIDC%NGZ
      NP=GRIDC%NGX*GRIDC%NGY*GRIDC%NGZ

      ALLOCATE(CWORK(N),CWORK2(N),WORK(NP))
!$ACC ENTER DATA CREATE(CWORK,CWORK2) __IF_ASYNC__
      ! merge magnetization density to CWORK
      DO I=1,3
         CALL MRG_GRID_RC(GRIDC, CWORK(NP*(I-1)+1), MAGFLD(1,I))
      ENDDO
      CALL MAGSYG(CWORK2,CWORK,GRIDC%NGX,GRIDC%NGX,GRIDC%NGY,GRIDC%NGY,&
               GRIDC%NGZ,GRIDC%NGZ,PTRANS,NIOND,MAGROT,SAXIS,LATT_CUR,WORK)
      ! and back the data go to MAGFLD
      DO I=1,3
         CALL DIS_GRID_RC(GRIDC, CWORK2(NP*(I-1)+1), MAGFLD(1,I), .FALSE.)
      ENDDO
!$ACC EXIT DATA DELETE(CWORK,CWORK2) __IF_ASYNC__
      DEALLOCATE(CWORK,CWORK2,WORK)

      RETURN
      END


!******************* SUBROUTINE RHOSYG *********************************
!
! subroutine rhosyg symmetrizes the charge density in reciprocal space.
!
! this routine needs the charge density on the rec. space FFT-grid (RHO)
!                    the dimension parameters of array RHO (ID1,ID2,ID3)
!                    the number of FFT-grid points (NR1,NR2,NR3)
!                    all primitive cells lying in the supercell (PTRANS)
!                    the dimension parameter for array PTRANS (NIOND)
!                    and a workarray (DONE)
!
! N.B.: the only sensible settings are:
!       ID1=NGX_rd, ID2=NGY,    ID3=NGZ_rd
!       NR1=NGX   , NR2=NGY,    NR3=NGZ
!
! TODO: remove the arguments ID1, ID2, and ID3 entirely.
!
!***********************************************************************

      SUBROUTINE RHOSYG(RHOG,RHOGIN, &
     &              ID1,NR1,ID2,NR2,ID3,NR3,PTRANS,NIOND,MAGROT,ISP,DONE)
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE prec


      IMPLICIT REAL(q) (A-H,O-Z)
      COMPLEX(q) RHOG(ID1*ID2*ID3),RHOGIN(ID1*ID2*ID3)
      COMPLEX(q) SUM,PHASG,PHASE(48),RHOGRD
      REAL(q) MAGROT(48,NIOND)
      INTEGER G0(3),SG(3),GRPOP,TINDEX(4,48)
      DIMENSION PTRANS(NIOND+2,3),DONE(NR1*NR2*NR3)
      LOGICAL LCONJG

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,GRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

!$ACC UPDATE SELF(RHOGIN) __IF_ASYNC__
      TWOPI=8._q*ATAN(1._q)
!=======================================================================
! Mark all points as 'not yet visited' ...
!=======================================================================
      DO IR=1,NR1*NR2*NR3
        DONE(IR)=-1._q
      ENDDO
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!=======================================================================
! Pass over all FFT-grid points ...
!=======================================================================
      DO IR3=0,NR3-1
       G0(3)=MOD(IR3+NR3/2-1,NR3)-NR3/2+1
       DO IR2=0,NR2-1
         G0(2)=MOD(IR2+NR2/2-1,NR2)-NR2/2+1
         DO 80 IR1=0,NR1-1
! If the vector has already been visited take the next vector ... :
            IF (DONE(NR1*(NR2*IR3+IR2)+IR1+1)>0._q) GOTO 80
            G0(1)=MOD(IR1+NR1/2-1,NR1)-NR1/2+1
! New star ...
            SUM=(0._q,0._q)
! NG counts the number of G-vector within the symmetry star for G0 ...
            NG=0
            NSTAR=0
            DO 50 IROT=1,NROTK
!-----------------------------------------------------------------------
! Sum up the charges of all grid-points which are symmetry related to G0
!-----------------------------------------------------------------------
               IROTI=INVMAP(IROT)
               SG(1)=GRPOP(1,1,IROTI)*G0(1)+GRPOP(2,1,IROTI)*G0(2)+ &
     &                                            GRPOP(3,1,IROTI)*G0(3)
               SG(2)=GRPOP(1,2,IROTI)*G0(1)+GRPOP(2,2,IROTI)*G0(2)+ &
     &                                            GRPOP(3,2,IROTI)*G0(3)
               SG(3)=GRPOP(1,3,IROTI)*G0(1)+GRPOP(2,3,IROTI)*G0(2)+ &
     &                                            GRPOP(3,3,IROTI)*G0(3)

               IF ((SG(1)<-NR1/2+1).OR.(SG(1)>NR1/2).OR. &
     &             (SG(2)<-NR2/2+1).OR.(SG(2)>NR2/2).OR. &
     &             (SG(3)<-NR3/2+1).OR.(SG(3)>NR3/2)) CYCLE

               SG(1)=MOD(SG(1)+6*NR1,NR1)
               SG(2)=MOD(SG(2)+6*NR2,NR2)
               SG(3)=MOD(SG(3)+6*NR3,NR3)
               IGRID =NR1*(NR2*SG(3)+SG(2))+SG(1)+1
#ifdef NGXhalf
               IF (SG(1)<= (NR1/2)) THEN
                 IG1=SG(1)
                 IG2=SG(2)
                 IG3=SG(3)
                 LCONJG=.FALSE.
               ELSE
                 IG1=MOD(NR1-SG(1),NR1)
                 IG2=MOD(NR2-SG(2),NR2)
                 IG3=MOD(NR3-SG(3),NR3)
                 LCONJG=.TRUE.
               ENDIF
               IGRID2=ID1*(ID2*IG3+IG2)+IG1+1
#elif defined(NGZhalf)
               IF (SG(3)<= (NR3/2)) THEN
                 IG1=SG(1)
                 IG2=SG(2)
                 IG3=SG(3)
                 LCONJG=.FALSE.
               ELSE
                 IG1=MOD(NR1-SG(1),NR1)
                 IG2=MOD(NR2-SG(2),NR2)
                 IG3=MOD(NR3-SG(3),NR3)
                 LCONJG=.TRUE.
               ENDIF
               IGRID2=ID1*(ID2*IG3+IG2)+IG1+1
#elif defined(NGYhalf)
               CALL ERROR('RHOSYG', 'NGYhalf mode not supported', 2)
#else
               IGRID2=IGRID
               LCONJG=.FALSE.
#endif
!-----------------------------------------------------------------------
! Phase factors for the summation due to nonprimitive translation ... :
!-----------------------------------------------------------------------
               GPHASE=TWOPI*(FLOAT(G0(1))*GTRANS(1,IROT)+ &
     &          FLOAT(G0(2))*GTRANS(2,IROT)+FLOAT(G0(3))*GTRANS(3,IROT))
               SCALE=1._q
               IF (ISP==2) SCALE=MAGROT(IROT,1)
               PHASGR=COS(GPHASE)*SCALE
               PHASGI=SIN(GPHASE)*SCALE
! WARNING: If we use nonprimitive supercells we must include all
! 'primitive' translations of the generating unitcell which are
! nonprimitive translations of the supercell ... :
               DO 20 ITRANS=2,NPCELL
                  GPHASE= &
     &                  ((GTRANS(1,IROT)+PTRANS(ITRANS,1))*FLOAT(G0(1))+ &
     &                   (GTRANS(2,IROT)+PTRANS(ITRANS,2))*FLOAT(G0(2))+ &
     &                   (GTRANS(3,IROT)+PTRANS(ITRANS,3))*FLOAT(G0(3)))
                  SCALE=1._q
                  IF (ISP==2) SCALE=MAGROT(IROT,ITRANS)
                  PHASGR=PHASGR+(COS(TWOPI*GPHASE)*SCALE)
                  PHASGI=PHASGI+(SIN(TWOPI*GPHASE)*SCALE)
   20          CONTINUE
               PHASGR=PHASGR/FLOAT(NPCELL)
               PHASGI=PHASGI/FLOAT(NPCELL)
               PHASG= CMPLX( PHASGR , PHASGI ,KIND=q)
! Symmetrized charge density for vector G0 summed up here ... :
               RHOGRD=RHOGIN(IGRID2)
               IF (LCONJG) RHOGRD=CONJG(RHOGRD)
               SUM=SUM+RHOGRD*CONJG(PHASG)

               NSTAR=NSTAR+1
!=======================================================================
! We need also the phase factors for the re-filling of the star ...
!=======================================================================
! Store the equivalent points and count the number of equivalent points:
               IF (DONE(IGRID)<0) THEN
! Here we found a new member of the star, so increase the counter NG ...
                  NG=NG+1
! ... remember coordinates for this member (needed for re-filling) ...
                  TINDEX(1,NG)=SG(1)
                  TINDEX(2,NG)=SG(2)
                  TINDEX(3,NG)=SG(3)
! ... and set the counter for the number of transformations which give
! us this member of the symmetry star (here first transformation) ...
                  TINDEX(4,NG)=1
! Intialise phase-factor table for this member (fill in the phase factor
! due to the nontrivial translation of the symmetry operation which has
! generated this new member of the symmetry star related to G0:
                  PHASE(NG)=PHASG
               ELSE
! Here the rotated vector has already been generated by some other
! symmetry operation, so we have now to sum up some information ...
                  DO 30 IG=1,NG
! Find out which member we have hit ... :
                     IF ((SG(1)==TINDEX(1,IG)).AND. &
     &                           (SG(2)==TINDEX(2,IG)).AND. &
     &                                  (SG(3)==TINDEX(3,IG))) GOTO 40
   30             CONTINUE
                  CALL ERROR('RHOSYG', 'stars are not distinct, try to increase SYMPREC to e.g. 1E-4', NG)
! Here IG should point to the member which represents vector SG ...
   40             CONTINUE
! Count the number of symmetry operations which generate this vector:
                  TINDEX(4,IG)=TINDEX(4,IG)+1
! Add up the phase factor due to the nontrivial translation vector for
! this symmetry operation ... :
                  PHASE(IG)=PHASE(IG)+PHASG
               ENDIF
! Remember visited points ... :
               DONE(IGRID)=1._q
   50       CONTINUE
!=======================================================================
! Supply correct scaling to phase factors for re-filling ...
!=======================================================================
!OCL VDOPT(VL(8))
            DO 60 IG=1,NG
               PHASE(IG)=PHASE(IG)/FLOAT(TINDEX(4,IG))
   60       CONTINUE
! ... and correct scaling for SUM!
!           SUM=SUM/FLOAT(NROTK)
            SUM=SUM/FLOAT(NSTAR)
!=======================================================================
!  SUM contains now the symmetrized charge density at point G0.
!  Now fill the star of G  with this sum ... :
!=======================================================================
!OCL VDOPT(VL(8))
            DO 70 IG=1,NG
#ifdef NGXhalf
               IF (TINDEX(1,IG)> (NR1/2)) GOTO 70
#elif defined(NGZhalf)
               IF (TINDEX(3,IG)> (NR3/2)) GOTO 70
#endif
               IGRID= &
     &            ID1*(ID2* &
     &            (TINDEX(3,IG))+(TINDEX(2,IG)))+TINDEX(1,IG)+1
               RHOG(IGRID)=SUM*PHASE(IG)
   70       CONTINUE
   80    CONTINUE
       ENDDO
      ENDDO
!$ACC UPDATE DEVICE(RHOG) __IF_ASYNC__
      RETURN
      END


!******************* SUBROUTINE MAGSYG *********************************
!
! subroutine magsyg symmetrizes the magnetization density in reciprocal space.
!
! this routine needs the magnetization density on the rec. space FFT-grid (MAGG)
!                    the dimension parameters of array MAGG(ID1,ID2,ID3,3)
!                    the number of FFT-grid points (NR1,NR2,NR3)
!                    all primitive cells lying in the supercell (PTRANS)
!                    the dimension parameter for array PTRANS (NIOND)
!                    and a workarray (DONE)
!
!***********************************************************************

      SUBROUTINE MAGSYG(MAGG,MAGGIN, &
     &              ID1,NR1,ID2,NR2,ID3,NR3,PTRANS,NIOND, &
     &               MAGROT,SAXIS,LATT_CUR,DONE)

#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE prec
      USE lattice

      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (latt) LATT_CUR

      COMPLEX(q)  MAGG(ID1*ID2*ID3,3),MAGGIN(ID1*ID2*ID3,3)
      COMPLEX(q)  PHASG,PHASE(48),MAGX,MAGY,MAGZ,MAGX_,MAGY_,MAGZ_
      COMPLEX(q)  MAGTX,MAGTY,MAGTZ
      COMPLEX(q)  MAGG_TEMP(3),MAGG_(3)
      REAL(q)     MAGROT(48,NIOND),SAXIS(3),ALPHA,BETA
      INTEGER     G0(3),SG(3),GRPOP,TINDEX(5,48)
      DIMENSION   PTRANS(NIOND+2,3),DONE(ID1*ID2*ID3)

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,GRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL
      LOGICAL    LDUMP
#ifndef oldsym
      INTEGER    DET
#endif

!$ACC UPDATE SELF(MAGGIN) __IF_ASYNC__
      TWOPI=8._q*ATAN(1._q)
      CALL EULER(SAXIS,ALPHA,BETA)
!=======================================================================
! Mark all points as 'not yet visited' ...
!=======================================================================
      DO IR=1,ID1*ID2*ID3
        DONE(IR)=-1._q
      ENDDO
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!=======================================================================
! Pass over all FFT-grid points ...
!=======================================================================
      DO IR3=0,NR3-1
       G0(3)=MOD(IR3+NR3/2-1,NR3)-NR3/2+1
       DO IR2=0,NR2-1
         G0(2)=MOD(IR2+NR2/2-1,NR2)-NR2/2+1
         DO 80 IR1=0,NR1-1
! If the vector has already been visited take the next vector ... :
            IF (DONE(ID1*(ID2*IR3+IR2)+IR1+1)>0._q) GOTO 80
            G0(1)=MOD(IR1+NR1/2-1,NR1)-NR1/2+1
            LDUMP = .FALSE.
#ifdef debug
            IF (ABS(G0(1)) <= 1 .AND. ABS(G0(2)) <= 1 .AND. ABS(G0(3)) <= 1) LDUMP=.TRUE.
#ifndef oldsym
            LDUMP=.TRUE.
#endif
            IF (LDUMP) WRITE(*,*) '-------------------------------------------------------'
#endif
            MAGX=0 ; MAGY=0; MAGZ=0
! New star ...
! NG counts the number of G-vector within the symmetry star for G0 ...
            NG=0
            DO 50 IROT=1,NROTK
!-----------------------------------------------------------------------
! Seek out all grid-points which are symmetry related to G0
!-----------------------------------------------------------------------
               IROTI=INVMAP(IROT)
               SG(1)=GRPOP(1,1,IROTI)*G0(1)+GRPOP(2,1,IROTI)*G0(2)+ &
     &                                            GRPOP(3,1,IROTI)*G0(3)
               SG(2)=GRPOP(1,2,IROTI)*G0(1)+GRPOP(2,2,IROTI)*G0(2)+ &
     &                                            GRPOP(3,2,IROTI)*G0(3)
               SG(3)=GRPOP(1,3,IROTI)*G0(1)+GRPOP(2,3,IROTI)*G0(2)+ &
     &                                            GRPOP(3,3,IROTI)*G0(3)
               SG(1)=MOD(SG(1)+6*NR1,NR1)
               SG(2)=MOD(SG(2)+6*NR2,NR2)
               SG(3)=MOD(SG(3)+6*NR3,NR3)
               IGRID =ID1*(ID2*SG(3)+SG(2))+SG(1)+1

               IF (LDUMP) THEN
                  WRITE(*,'(3I3)') MOD(SG(1)+NR1/2,NR1)-NR1/2, &
                  MOD(SG(2)+NR2/2,NR2)-NR2/2,MOD(SG(3)+NR3/2,NR3)-NR3/2
                  WRITE(*,'("retrive  ",I6)') IGRID
#ifndef oldsym
                  WRITE(*,'("irot     ",I6)') IROT
#endif
               ENDIF
!-----------------------------------------------------------------------
! Phase factors for the summation due to nonprimitive translation ... :
!-----------------------------------------------------------------------
               GPHASE=TWOPI*(FLOAT(G0(1))*GTRANS(1,IROT)+ &
     &          FLOAT(G0(2))*GTRANS(2,IROT)+FLOAT(G0(3))*GTRANS(3,IROT))
               SCALE=1._q
               PHASGR=COS(GPHASE)*SCALE
               PHASGI=SIN(GPHASE)*SCALE
! WARNING: If we use nonprimitive supercells we must include all
! 'primitive' translations of the generating unitcell which are
! nonprimitive translations of the supercell ... :
               DO 20 ITRANS=2,NPCELL
                  GPHASE= &
     &                  ((GTRANS(1,IROT)+PTRANS(ITRANS,1))*FLOAT(G0(1))+ &
     &                   (GTRANS(2,IROT)+PTRANS(ITRANS,2))*FLOAT(G0(2))+ &
     &                   (GTRANS(3,IROT)+PTRANS(ITRANS,3))*FLOAT(G0(3)))
                  SCALE=1._q
                  PHASGR=PHASGR+(COS(TWOPI*GPHASE)*SCALE)
                  PHASGI=PHASGI+(SIN(TWOPI*GPHASE)*SCALE)
   20          CONTINUE
               PHASGR=PHASGR/FLOAT(NPCELL)
               PHASGI=PHASGI/FLOAT(NPCELL)
               PHASG= CMPLX( PHASGR , PHASGI ,KIND=q)
! Transform from "SAXIS basis" to the system of cartesian axes
! in which the integer rotation matrices are defined
               IF (LDUMP) WRITE(*,'(9X,6F14.7)') MAGGIN(IGRID,:)
               MAGG_TEMP(1)=COS(BETA)*COS(ALPHA)*MAGGIN(IGRID,1)- &
            &              SIN(ALPHA)*MAGGIN(IGRID,2)+ &
            &               SIN(BETA)*COS(ALPHA)*MAGGIN(IGRID,3)
               MAGG_TEMP(2)=COS(BETA)*SIN(ALPHA)*MAGGIN(IGRID,1)+ &
            &              COS(ALPHA)*MAGGIN(IGRID,2)+ &
            &               SIN(BETA)*SIN(ALPHA)*MAGGIN(IGRID,3)
               MAGG_TEMP(3)=-SIN(BETA)*MAGGIN(IGRID,1)+ &
            &               COS(BETA)*MAGGIN(IGRID,3)
               IF (LDUMP) WRITE(*,'(9X,6F14.7)') MAGG_TEMP
! Go to reciprocal coordinates
               CALL CKARDIR(1,MAGG_TEMP,LATT_CUR%A)
! Rotate to the orientation of G0
#ifdef oldsym
               MAGTX=GRPOP(1,1,IROTI)*MAGG_TEMP(1)+ &
              &       GRPOP(2,1,IROTI)*MAGG_TEMP(2)+ &
              &        GRPOP(3,1,IROTI)* MAGG_TEMP(3)
               MAGTY=GRPOP(1,2,IROTI)*MAGG_TEMP(1)+ &
              &       GRPOP(2,2,IROTI)*MAGG_TEMP(2)+ &
              &        GRPOP(3,2,IROTI)*MAGG_TEMP(3)
               MAGTZ=GRPOP(1,3,IROTI)*MAGG_TEMP(1)+ &
              &       GRPOP(2,3,IROTI)*MAGG_TEMP(2)+ &
              &        GRPOP(3,3,IROTI)*MAGG_TEMP(3)
#else
               MAGTX=GRPOP(1,1,IROT)*MAGG_TEMP(1)+ &
                    &       GRPOP(2,1,IROT)*MAGG_TEMP(2)+ &
                    &        GRPOP(3,1,IROT)* MAGG_TEMP(3)
               MAGTY=GRPOP(1,2,IROT)*MAGG_TEMP(1)+ &
                    &       GRPOP(2,2,IROT)*MAGG_TEMP(2)+ &
                    &        GRPOP(3,2,IROT)*MAGG_TEMP(3)
               MAGTZ=GRPOP(1,3,IROT)*MAGG_TEMP(1)+ &
                    &       GRPOP(2,3,IROT)*MAGG_TEMP(2)+ &
                    &        GRPOP(3,3,IROT)*MAGG_TEMP(3)

               DET=GRPOP(1,1,IROT)*GRPOP(2,2,IROT)*GRPOP(3,3,IROT)- &
                   GRPOP(1,1,IROT)*GRPOP(2,3,IROT)*GRPOP(3,2,IROT)+ &
                   GRPOP(1,2,IROT)*GRPOP(2,3,IROT)*GRPOP(3,1,IROT)- &
                   GRPOP(1,2,IROT)*GRPOP(2,1,IROT)*GRPOP(3,3,IROT)+ &
                   GRPOP(1,3,IROT)*GRPOP(2,1,IROT)*GRPOP(3,2,IROT)- &
                   GRPOP(1,3,IROT)*GRPOP(2,2,IROT)*GRPOP(3,1,IROT)
               MAGTX=DET*MAGTX
               MAGTY=DET*MAGTY
               MAGTZ=DET*MAGTZ
#endif
! And back to cartesian coordinates
               MAGG_TEMP(1)=MAGTX
               MAGG_TEMP(2)=MAGTY
               MAGG_TEMP(3)=MAGTZ
               CALL CDIRKAR(1,MAGG_TEMP,LATT_CUR%B)

               IF (LDUMP) WRITE(*,'(9X,6F14.7)') MAGG_TEMP
! And back to SAXIS representation
               MAGTX=COS(BETA)*COS(ALPHA)*MAGG_TEMP(1)+ &
            &       COS(BETA)*SIN(ALPHA)*MAGG_TEMP(2)- &
            &        SIN(BETA)*MAGG_TEMP(3)
               MAGTY=-SIN(ALPHA)*MAGG_TEMP(1)+ &
            &        COS(ALPHA)*MAGG_TEMP(2)
               MAGTZ=SIN(BETA)*COS(ALPHA)*MAGG_TEMP(1)+ &
            &       SIN(BETA)*SIN(ALPHA)*MAGG_TEMP(2)+ &
            &        COS(BETA)*MAGG_TEMP(3)

               IF (LDUMP) WRITE(*,'(9X,6F14.7)') MAGTX, MAGTY, MAGTZ
               IF (LDUMP) WRITE(*,'("  phase  ",6F14.7)') PHASG

! Symmetrized magnetization density for vector G0 summed up here ... :
               MAGX=MAGX+MAGTX*CONJG(PHASG)
               MAGY=MAGY+MAGTY*CONJG(PHASG)
               MAGZ=MAGZ+MAGTZ*CONJG(PHASG)
!=======================================================================
! We need also the phase factors for the re-filling of the star ...
!=======================================================================
! Store the equivalent points and count the number of equivalent points:
               IF (DONE(IGRID)<0) THEN
! Here we found a new member of the star, so increase the counter NG ...
                  NG=NG+1
! ... remember coordinates for this member (needed for re-filling) ...
                  TINDEX(1,NG)=SG(1)
                  TINDEX(2,NG)=SG(2)
                  TINDEX(3,NG)=SG(3)
! ... and set the counter for the number of transformations which give
! us this member of the symmetry star (here first transformation) ...
                  TINDEX(4,NG)=1
! ... and store the rotation that brings the rotated vector back to its
! original orientation as is needed when refilling the star
#ifdef oldsym
                  TINDEX(5,NG)=IROT
#else
                  TINDEX(5,NG)=IROTI
#endif
! Intialise phase-factor table for this member (fill in the phase factor
! due to the nontrivial translation of the symmetry operation which has
! generated this new member of the symmetry star related to G0:
                  PHASE(NG)=PHASG
               ELSE
! Here the rotated vector has already been generated by some other
! symmetry operation, so we have now to sum up some information ...
                  DO 30 IG=1,NG
! Find out which member we have hit ... :
                     IF ((SG(1)==TINDEX(1,IG)).AND. &
     &                           (SG(2)==TINDEX(2,IG)).AND. &
     &                                  (SG(3)==TINDEX(3,IG))) GOTO 40
   30             CONTINUE
! Here IG should point to the member which represents vector SG ...
   40             CONTINUE
! Count the number of symmetry operations which generate this vector:
                  TINDEX(4,IG)=TINDEX(4,IG)+1
! Add up the phase factor due to the nontrivial translation vector for
! this symmetry operation ... :
#ifdef oldsym
                  PHASE(IG)=PHASE(IG)+PHASG
#endif
               ENDIF
! Remember visited points ... :
               DONE(IGRID)=1._q
   50       CONTINUE
!=======================================================================
! Supply correct scaling to phase factors for re-filling ...
!=======================================================================
#ifdef oldsym
            DO 60 IG=1,NG
               PHASE(IG)=PHASE(IG)/FLOAT(TINDEX(4,IG))
   60       CONTINUE
#endif
! ... and correct scaling for SUM!
            MAGX=MAGX/FLOAT(NROTK)
            MAGY=MAGY/FLOAT(NROTK)
            MAGZ=MAGZ/FLOAT(NROTK)
            IF (LDUMP) WRITE(*,'(9X,6F14.7)') MAGX, MAGY, MAGZ
!=======================================================================
!  SUM contains now the symmetrized magnetization density at point G0.
!  Now fill the star of G  with this sum ... :
!=======================================================================
            DO 70 IG=1,NG
               IGRID =ID1*(ID2*TINDEX(3,IG)+TINDEX(2,IG))+TINDEX(1,IG)+1
               IF (LDUMP) WRITE(*,'("   store ",I6)') IGRID
#ifndef oldsym
               IF (LDUMP) WRITE(*,'("   phase ",2F14.7)') phase(ig)
               IF (LDUMP) WRITE(*,'("   irot ",I6)') tindex(5,ig)
#endif
               MAGX_=MAGX*PHASE(IG)
               MAGY_=MAGY*PHASE(IG)
               MAGZ_=MAGZ*PHASE(IG)
! Transform from "SAXIS basis" to the system of cartesian axes
! in which the integer rotation matrices are defined
               IF (LDUMP) WRITE(*,'(9X,6F14.7)') MAGG_,MAGY_,MAGZ_
               MAGG_TEMP(1)=COS(BETA)*COS(ALPHA)*MAGX_- &
            &              SIN(ALPHA)*MAGY_+ &
            &               SIN(BETA)*COS(ALPHA)*MAGZ_
               MAGG_TEMP(2)=COS(BETA)*SIN(ALPHA)*MAGX_+ &
            &              COS(ALPHA)*MAGY_+ &
            &               SIN(BETA)*SIN(ALPHA)*MAGZ_
               MAGG_TEMP(3)=-SIN(BETA)*MAGX_+ &
            &               COS(BETA)*MAGZ_
! Go to reciprocal coordinates
               IF (LDUMP) WRITE(*,'(9X,6F14.7)') MAGG_TEMP
               CALL CKARDIR(1,MAGG_TEMP,LATT_CUR%A)
! Rotate back to the original orientation
               MAGG_(1)=GRPOP(1,1,TINDEX(5,IG))*MAGG_TEMP(1)+ &
              &               GRPOP(2,1,TINDEX(5,IG))*MAGG_TEMP(2)+ &
              &                GRPOP(3,1,TINDEX(5,IG))*MAGG_TEMP(3)
               MAGG_(2)=GRPOP(1,2,TINDEX(5,IG))*MAGG_TEMP(1)+ &
              &               GRPOP(2,2,TINDEX(5,IG))*MAGG_TEMP(2)+ &
              &                GRPOP(3,2,TINDEX(5,IG))*MAGG_TEMP(3)
               MAGG_(3)=GRPOP(1,3,TINDEX(5,IG))*MAGG_TEMP(1)+ &
              &               GRPOP(2,3,TINDEX(5,IG))*MAGG_TEMP(2)+ &
              &                GRPOP(3,3,TINDEX(5,IG))*MAGG_TEMP(3)
#ifndef oldsym
               DET=GRPOP(1,1,TINDEX(5,IG))*&
                   GRPOP(2,2,TINDEX(5,IG))*&
                   GRPOP(3,3,TINDEX(5,IG)) - &
                   GRPOP(1,1,TINDEX(5,IG))*&
                   GRPOP(2,3,TINDEX(5,IG))*&
                   GRPOP(3,2,TINDEX(5,IG)) + &
                   GRPOP(1,2,TINDEX(5,IG))*&
                   GRPOP(2,3,TINDEX(5,IG))*&
                   GRPOP(3,1,TINDEX(5,IG)) - &
                   GRPOP(1,2,TINDEX(5,IG))*&
                   GRPOP(2,1,TINDEX(5,IG))*&
                   GRPOP(3,3,TINDEX(5,IG)) + &
                   GRPOP(1,3,TINDEX(5,IG))*&
                   GRPOP(2,1,TINDEX(5,IG))*&
                   GRPOP(3,2,TINDEX(5,IG)) - &
                   GRPOP(1,3,TINDEX(5,IG))*&
                   GRPOP(2,2,TINDEX(5,IG))*&
                   GRPOP(3,1,TINDEX(5,IG))
               MAGG_=MAGG_*DET
#endif
! And back to cartesian coordinates
               CALL CDIRKAR(1,MAGG_(1:3),LATT_CUR%B)
               IF (LDUMP) WRITE(*,'(9X,6F14.7)') MAGG_
! And back to SAXIS representation
               MAGG(IGRID,1)=COS(BETA)*COS(ALPHA)*MAGG_(1)+ &
            &               COS(BETA)*SIN(ALPHA)*MAGG_(2)- &
            &                SIN(BETA)*MAGG_(3)
               MAGG(IGRID,2)=-SIN(ALPHA)*MAGG_(1)+ &
            &                COS(ALPHA)*MAGG_(2)
               MAGG(IGRID,3)=SIN(BETA)*COS(ALPHA)*MAGG_(1)+ &
            &               SIN(BETA)*SIN(ALPHA)*MAGG_(2)+ &
            &                COS(BETA)*MAGG_(3)
            IF (LDUMP) WRITE(*,'(9X,6F14.7)') MAGG(IGRID,:)
   70       CONTINUE
   80    CONTINUE
      ENDDO
      ENDDO
!$ACC UPDATE DEVICE(MAGG) __IF_ASYNC__
      RETURN
      END

!**************** SUBROUTINE EULER  ************************************
!
!> @brief
!> Calculate the Euler angles (alpha,beta)_i
!> for the site-dependent spin directions (sx,sy,sz)_i.
!> beta  is the angle between the z axis and s
!> alpha is the angle between x and (sx,sy,0)
!>
!> These angles are used to determine a transformation that corresponds to
!> 1. rotation around y-axis of beta
!> 2. rotation around z-axis of alpha
!***********************************************************************

      SUBROUTINE EULER(S, ALPHA, BETA)

      USE prec
      USE constant

      IMPLICIT NONE
      REAL(q) T,T2,S(3), &
     &        ALPHA,BETA
! local

!       Normalisation of the spin directions
        T=SQRT(S(1)**2+S(2)**2+S(3)**2)

        IF (T<1.0E-10_q) THEN
          ALPHA=0
          BETA=0
        ELSE

          S(1)=S(1)/T
          S(2)=S(2)/T
          S(3)=S(3)/T

          T2=SQRT(S(1)**2+S(2)**2)
!       calculates the Euler angles (alpha,beta) for
!       each ions.
          IF(ABS(S(1))>=1.0E-10_q)THEN
            ALPHA=ATAN2(S(2),S(1))
          ELSE
!       in case sx(i)=0, avoids division by 0
            IF(ABS(S(2))>1.0E-10_q)THEN
              ALPHA=pi/2
            ENDIF
            IF(ABS(S(2))<1.0E-10_q)THEN
              ALPHA=0
            ENDIF
          ENDIF
          IF(ABS(S(3))>1.0E-10_q)THEN
            BETA=ATAN2(T2,S(3))
          ELSE
!       in case sz(i)=0, cos(beta(i))=0
            BETA=pi/2
          ENDIF
        ENDIF

      END SUBROUTINE EULER


!******************* SUBROUTINE RHOSYR *********************************
!
! subroutine rhosyr symmetrizes the charge density in direct space. This
! is easy and fast (it saves some FFTs) but it can only be used in the
! case of pure point symmetry - if we have nontrivial translations for
! certain symmetry elements we must use routine RHOSYG!
!
! this routine needs the charge density on the real space FFT-grid (RHO)
!                    the dimension parameters of array RHO (ID1,ID2,ID3)
!                    the number of FFT-grid points (NR1,NR2,NR3)
!                    and a workarray (DONE)
!
!***********************************************************************

      SUBROUTINE RHOSYR(RHO,ID1,NR1,ID2,NR2,ID3,NR3,MAGROT,ISP,DONE)
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE prec


      IMPLICIT REAL(q) (A-H,O-Z)
      REAL(q)   RHO(ID1*ID2*ID3),MAGROT(48)
      DIMENSION DONE(ID1,ID2,ID3)
      INTEGER R(3),SR(3),S

      COMMON /SYMM/ S(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

!$ACC UPDATE SELF(RHO) __IF_ASYNC__
!=======================================================================
! Mark all points as 'not yet visited' ...
!=======================================================================
      DO IR=1,ID1*ID2*ID3
         DONE(IR,1,1)=-1._q
      ENDDO
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
!=======================================================================
! Pass over all FFT-grid points ...
!=======================================================================
      DO IR3=1,NR3
       R(3)=IR3-1
       DO IR2=1,NR2
         R(2)=IR2-1
         DO 40 IR1=1,NR1
! If the vector has already been visited take the next vector ... :
            IF (DONE(IR1,IR2,IR3)>0._q) GOTO 40
            R(1)=IR1-1
! Note: in direct space the charge density is always real! So we use a
            SUM=0._q
!-----------------------------------------------------------------------
! Sum up the charges of all grid-points which are symmetry related to R
!-----------------------------------------------------------------------
            DO 20 IROT=1,NROT
               SCALE=1._q
               IF (ISP==2) SCALE=MAGROT(IROT)
               SR(1)=S(1,1,IROT)*R(1)+S(2,1,IROT)*R(2)+S(3,1,IROT)*R(3)
               SR(2)=S(1,2,IROT)*R(1)+S(2,2,IROT)*R(2)+S(3,2,IROT)*R(3)
               SR(3)=S(1,3,IROT)*R(1)+S(2,3,IROT)*R(2)+S(3,3,IROT)*R(3)
               SR(1)=MOD(SR(1)+6*NR1,NR1)
               SR(2)=MOD(SR(2)+6*NR2,NR2)
               SR(3)=MOD(SR(3)+6*NR3,NR3)
               IGRID=ID1*(ID2*SR(3)+SR(2))+SR(1)+1
               SUM=SUM+(RHO(IGRID)*SCALE)
   20       CONTINUE
! Correct scaling for SUM!
            SUM=SUM/FLOAT(NROT)
!-----------------------------------------------------------------------
!  SUM contains now the symmetrized charge density at point R.
!  Now fill the star of all symmetry related Rs  with this sum.
!-----------------------------------------------------------------------
            DO 30 IROT=1,NROT
               SCALE=1._q
               IF (ISP==2) SCALE=MAGROT(IROT)
               SR(1)=S(1,1,IROT)*R(1)+S(2,1,IROT)*R(2)+S(3,1,IROT)*R(3)
               SR(2)=S(1,2,IROT)*R(1)+S(2,2,IROT)*R(2)+S(3,2,IROT)*R(3)
               SR(3)=S(1,3,IROT)*R(1)+S(2,3,IROT)*R(2)+S(3,3,IROT)*R(3)
               SR(1)=MOD(SR(1)+6*NR1,NR1)
               SR(2)=MOD(SR(2)+6*NR2,NR2)
               SR(3)=MOD(SR(3)+6*NR3,NR3)
               IGRID=ID1*(ID2*SR(3)+SR(2))+SR(1)+1
! Remember the grid points which we have already visited ... :
               DONE(SR(1)+1,SR(2)+1,SR(3)+1)=1._q
               RHO(IGRID)=SUM*SCALE
   30       CONTINUE
   40    CONTINUE
      ENDDO
      ENDDO
!$ACC UPDATE DEVICE(RHO) __IF_ASYNC__
      RETURN
      END


!******************* SUBROUTINE POSSYM *********************************
!
! subroutine POSSYM is an interface to the lattlib-routine PSYM
! this routine is used to (re-)symmetrize the atomic positions
!
!***********************************************************************

      SUBROUTINE POSSYM(POSION,ROTMAP,PTRANS,NIOND,NTYP,NITYP,FROT)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      INTEGER ROTMAP(NIOND,NROTK,NPCELL)

      DIMENSION POSION(3,NIOND),PTRANS(NIOND+2,3),NITYP(NTYP)
      DIMENSION TAU(NIOND,3),TAUWRK(NIOND,3),FROT(3,NIOND)

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

      ! copy the positions to TAU
      DO I=1,NIOND
         TAU(I,:)=POSION(:,I)
      ENDDO

      CALL PSYM(TAU,ROTMAP,ISYMOP,GTRANS,PTRANS,NROTK,NPCELL,NIOND,1,NTYP,NITYP,TAUWRK)

      ! copy the symmetrized positions to FROT
      DO I=1,NIOND
         FROT(:,I)=TAUWRK(I,:)
      ENDDO

      RETURN
      END

!******************* SUBROUTINE FORSYM *********************************
!
! subroutine forsym is a simplified interface to lattlib-routine FSYM
! this routine serves generally for symmetrisation of a vector field
! F(IATOM) given in   c a r t e s i a n   coordinates.
!
! this routine requires the vector-field to be symmetrized (FION)
!             the connection map (ROTMAP) set up in routine INISYM
!             the number of atomic species (NTYP)
!             the number of atoms of each species (NITYP)
!             some workarrays (FROT,WRKROT)
!             the dimensioning parameter NIOND (max. no. of atoms)
!             and the lattice vectors (A)
!
! the input data will be replaced by the symmetric data on output
!
!***********************************************************************

      SUBROUTINE FORSYM(FION,ROTMAP,NTYP,NITYP,NIOND,FROT,WRKROT,A)
      USE prec


      IMPLICIT REAL(q) (A-H,O-Z)
      INTEGER ROTMAP(NIOND,NROTK,NPCELL)

      DIMENSION FION(3,NIOND),NITYP(NTYP),FROT(NIOND,3),WRKROT(NIOND,3)
      DIMENSION A(3,3),B1(3),B2(3),B3(3)

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

      CALL RECIPS(1._q,A(1,1),A(1,2),A(1,3),B1,B2,B3)
      CALL FSYM(FION,ROTMAP,ISYMOP,NROTK,NPCELL,1,NTYP,NIOND,NITYP, &
     &                        FROT,WRKROT,A(1,1),A(1,2),A(1,3),B1,B2,B3)

      RETURN
      END


!******************* SUBROUTINE VECSYM *********************************
!
! subroutine vecsym is a simplified interface to lattlib-routine FSYM
! this routine serves generally for symmetrisation of a vector field
! V(IATOM) given in   d i r e c t   lattice coordinates.
!
! this routine requires the vector-field to be symmetrized (FION)
!             the connection map (ROTMAP) set up in routine INISYM
!             the number of atomic species (NTYP)
!             the number of atoms of each species (NITYP)
!             some workarrays (FROT,WRKROT)
!             and the dimensioning parameter NIOND (max. no. of atoms)
!
! the input data will be replaced by the symmetric data on output
!
!***********************************************************************

      SUBROUTINE VECSYM(VEC,ROTMAP,NTYP,NITYP,NIOND,FROT,WRKROT)
      USE prec


      IMPLICIT REAL(q) (A-H,O-Z)
      INTEGER ROTMAP(NIOND,NROTK,NPCELL)

      DIMENSION VEC(3,NIOND),NITYP(NTYP),FROT(NIOND,3),WRKROT(NIOND,3)
      DIMENSION A1(3),A2(3),A3(3),B1(3),B2(3),B3(3)

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

      SAVE A1,A2,A3,B1,B2,B3
      DATA A1 /1._q,0._q,0._q/, A2 /0._q,1._q,0._q/, A3 /0._q,0._q,1._q/
      DATA B1 /1._q,0._q,0._q/, B2 /0._q,1._q,0._q/, B3 /0._q,0._q,1._q/

      CALL FSYM(VEC,ROTMAP,ISYMOP,NROTK,NPCELL,1,NTYP,NIOND,NITYP, &
     &                                    FROT,WRKROT,A1,A2,A3,B1,B2,B3)

      RETURN
      END


!******************* SUBROUTINE TENSYM *********************************
!
! subroutine tensym is a simplified interface to symlib-routine TSYM2
! this routine serves to symmetrize of a site specific tensor
! TION(IATOM) given in   c a r t e s i a n   coordinates.
!
! this routine requires the tensors to be symmetrized (TION)
!             the connection map (ROTMAP) set up in routine INISYM
!             the number of atomic species (NTYP)
!             the number of atoms of each species (NITYP)
!             some workarrays (TROT,WRKROT)
!             the dimensioning parameter NIOND (max. no. of atoms)
!             and the lattice vectors (A)
!
! the input data will be replaced by the symmetric data on output
!
!***********************************************************************

      SUBROUTINE TENSYM(TION,ROTMAP,NTYP,NITYP,NIOND,A)
      USE prec


      IMPLICIT REAL(q) (A-H,O-Z)
      INTEGER ROTMAP(NIOND,NROTK,NPCELL)

      DIMENSION TION(3,3,NIOND),NITYP(NTYP),TROT(NIOND,3,3),WRKROT(NIOND,3,3)
      DIMENSION A(3,3),B1(3),B2(3),B3(3)

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

      CALL RECIPS(1._q,A(1,1),A(1,2),A(1,3),B1,B2,B3)
      CALL TSYM2(TION,ROTMAP,ISYMOP,NROTK,NPCELL,1,NTYP,NIOND,NITYP, &
     &                        TROT,WRKROT,A(1,1),A(1,2),A(1,3),B1,B2,B3)

      RETURN
      END


!************************* SUBROUTINE COUNTER_SPIRAL *******************
!
! This subroutine is called in case one does spin spiral calculations.
! It rotates the initial magnetic moments clockwise in the xy
! plane, according to
!
!    M_rot_x = M_x cos(q.R) + M_y sin(q.R)
!    M_rot_y = M_y cos(q.R) - M_x sin(q.R)
!    M_rot_z = M_z
!
! where q is the spin spiral propagation vector, and R is the position
! to which the initial local moment refers (i.e. an atomic position).
!
! This counters the following counter clockwise rotation of the moments
! in the xy plane along the propagation direction of the spiral.
!
!***********************************************************************

      SUBROUTINE COUNTER_SPIRAL(QSPIRAL,NIONS,POSION,ATOMOM)

      USE prec
      USE constant

      IMPLICIT NONE

      INTEGER NIONS,I
      REAL(q) QSPIRAL(3)
      REAL(q) POSION(3,NIONS),ATOMOM(3*NIONS)
      REAL(q) QR,M_x,M_y

      DO I=1,NIONS
         QR=TPI*(QSPIRAL(1)*POSION(1,I)+QSPIRAL(2)*POSION(2,I)+QSPIRAL(3)*POSION(3,I))

         M_x=ATOMOM(1+(I-1)*3)*COS(QR)+ATOMOM(2+(I-1)*3)*SIN(QR)
         M_y=ATOMOM(2+(I-1)*3)*COS(QR)-ATOMOM(1+(I-1)*3)*SIN(QR)

         ATOMOM(1+(I-1)*3)=M_x
         ATOMOM(2+(I-1)*3)=M_y
      ENDDO
      END SUBROUTINE COUNTER_SPIRAL


      SUBROUTINE ERROR(subName, msg, errCode)
      USE prec
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: subName, msg
      INTEGER, INTENT(IN) :: errCode
      CALL vtutor%error('VERY BAD NEWS! internal error in subroutine' // subName // &
         ':' // msg // ' ' // str(errCode))
      END SUBROUTINE ERROR

      SUBROUTINE WARN(A,B,N)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      CHARACTER (LEN=*) A,B
      WRITE(*,1) A,B,N
    1 FORMAT(//' VERY BAD NEWS! internal error in subroutine',A,':'/,A,I8)
      END

!***********************************************************************
!
!   Routine POLSYM symmetrizes a vectorfield given in
!   cartesian coordinates according to the current symmetry operations
!   stored in the COMMON block /SYMM/
!
!   Input parameters:
!   -----------------
!
!      F(3)   contains on input the vector to be symmetrized
!             (given in crystal coordinates)
!
!      A(3,3) contains the lattice vectors defining the unit cell
!
!***********************************************************************

      SUBROUTINE POLSYM(F,A)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)

      REAL(q) F(3), A(3,3)


! symmetry related quantities (common block)
      INTEGER  ISYMOP,NROT,IGRPOP,NROTK,INVMAP,NPCELL
      REAL(q)  GTRANS,AP
      COMMON /SYMM/   ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
           GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

      CALL F1SYM(F,ISYMOP,NROTK,A)

      END SUBROUTINE


!******************* SUBROUTINE PARSYL *********************************
!
! Subroutine PARSYL symmetrizes l+site resolved partial charges
! or the corresponding local magnetic moments (when ISP=2),
! for l=0,..,LDIM
!
! These charges/moments must be supplied in the array:
!
!   VEC(1:LDIM+1,1:NIONS)
!
! On exit this array will contain the symmetrized result.

!***********************************************************************

      SUBROUTINE PARSYL(LDIM,NIONS,NIOND,NTYP,NITYP,VEC,ROTMAP,MAGROT,ISP)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)

      INTEGER :: LDIM,NIONS,NIOND,NTYP,ISP
      INTEGER :: NITYP(NTYP),ROTMAP(NIOND,NROTK,NPCELL)
      REAL(q) :: MAGROT(48,NPCELL)
      REAL(q) :: VEC(LDIM+1,NIONS)
! local variables
      INTEGER :: NROT,ITRANS,IA,IAP,IROT,IS,ISTART
      REAL(q) :: SCALE

      REAL(q), ALLOCATABLE :: TMP(:,:)

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

!-------------------------------------------------------------------
! allocate work array
!-------------------------------------------------------------------
      ALLOCATE(TMP(LDIM+1,NIONS))
      TMP=0
!-------------------------------------------------------------------
! do the symmetrization
!-------------------------------------------------------------------
      ISTART=1
      ! loop over all species
      DO IS=1,NTYP
        IF (IS>1) ISTART=ISTART+NITYP(IS-1)
        ! loop over all rotations
        DO IROT=1,NROTK
          ! loop over all ions
          DO IA=ISTART,NITYP(IS)+ISTART-1
            ! loop over all translations
            DO ITRANS=1,NPCELL
               ! destination atom
               IAP=ROTMAP(IA,IROT,ITRANS)
               SCALE=1._q
               IF (ISP==2) SCALE=MAGROT(IROT,ITRANS)
               TMP(:,IA)=TMP(:,IA)+VEC(:,IAP)*SCALE
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      ! divide final result by the number of translations and rotations
      SCALE=1._q/(NPCELL*NROTK)
      VEC=TMP*SCALE

      DEALLOCATE(TMP)

      END SUBROUTINE PARSYL


!******************* SUBROUTINE MAGSYL *********************************
!
! Subroutine MAGSYLM symmetrizes l+site resolved noncollinear
! local magnetic moments, for l=0,..,LDIM, and m=1,2*l+1
!
! These moments must be supplied in the array:
!
!   VEC(1:LDIM+1,1:NIONS,1:3)
!
! On exit this array will contain the symmetrized result.
!
!***********************************************************************

      SUBROUTINE MAGSYL(LDIM,NIONS,NIOND,NTYP,NITYP,VEC,ROTMAP,MAGROT,SAXIS,A,B)
      USE prec
      USE lattice
      IMPLICIT REAL(q) (A-H,O-Z)

      INTEGER :: LDIM,NIONS,NIOND,NTYP,ISP
      INTEGER :: NITYP(NTYP),ROTMAP(NIOND,NROTK,NPCELL)

      REAL(q) :: MAGROT(48,NPCELL),A(3,3),B(3,3)
      REAL(q) :: VEC(LDIM+1,NIONS,3),SAXIS(3)

! local variables
      INTEGER :: NROT,ITRANS,IA,IAP,IROT,IS,ISTART
      REAL(q) :: ALPHA,BETA,SCALE

      REAL(q), ALLOCATABLE :: TMP(:,:,:),TROTVEC1(:,:),TROTVEC2(:,:)
#ifndef oldsym
      INTEGER DET,TMPM(3,3)
#endif

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

      CALL EULER(SAXIS,ALPHA,BETA)

!-------------------------------------------------------------------
! allocate work arrays
!-------------------------------------------------------------------
      ALLOCATE(TMP(LDIM+1,NIONS,3),TROTVEC1(3,LDIM+1),TROTVEC2(3,LDIM+1))

      TMP=0
!-------------------------------------------------------------------
! do the symmetrization
!-------------------------------------------------------------------
      ISTART=1
      ! loop over all species
      DO IS=1,NTYP
        IF (IS>1) ISTART=ISTART+NITYP(IS-1)
        ! loop over all rotations
        DO IROT=1,NROTK
          ! loop over all ions
          DO IA=ISTART,NITYP(IS)+ISTART-1
            ! loop over all translations
            DO ITRANS=1,NPCELL
              ! destination atom
              IAP=ROTMAP(IA,IROT,ITRANS)
              SCALE=1._q

! Transform from "SAXIS basis" to the system of cartesian axes
              TROTVEC1(1,:)=COS(BETA)*COS(ALPHA)*VEC(:,IAP,1)- &
                   SIN(ALPHA)*VEC(:,IAP,2)+ &
                   SIN(BETA)*COS(ALPHA)*VEC(:,IAP,3)
              TROTVEC1(2,:)=COS(BETA)*SIN(ALPHA)*VEC(:,IAP,1)+ &
                   COS(ALPHA)*VEC(:,IAP,2)+ &
                   SIN(BETA)*SIN(ALPHA)*VEC(:,IAP,3)
              TROTVEC1(3,:)=-SIN(BETA)*VEC(:,IAP,1)+ &
                   COS(BETA)*VEC(:,IAP,3)
! Bring to direct coordinates in which the integer rotation matrices are defined
              CALL KARDIR(LDIM+1,TROTVEC1,B)
#ifndef oldsym
              TMPM=ISYMOP(:,:,IROT)
              DET=TMPM(1,1)*TMPM(2,2)*TMPM(3,3) - &
                  TMPM(1,1)*TMPM(2,3)*TMPM(3,2) + &
                  TMPM(1,2)*TMPM(2,3)*TMPM(3,1) - &
                  TMPM(1,2)*TMPM(2,1)*TMPM(3,3) + &
                  TMPM(1,3)*TMPM(2,1)*TMPM(3,2) - &
                  TMPM(1,3)*TMPM(2,2)*TMPM(3,1)
              IF (DET<0) TMPM=-TMPM

              TROTVEC2(1,:)= &
                   TROTVEC1(1,:)*TMPM(1,1)+ &
                   TROTVEC1(2,:)*TMPM(2,1)+ &
                   TROTVEC1(3,:)*TMPM(3,1)

              TROTVEC2(2,:)= &
                   TROTVEC1(1,:)*TMPM(1,2)+ &
                   TROTVEC1(2,:)*TMPM(2,2)+ &
                   TROTVEC1(3,:)*TMPM(3,2)

              TROTVEC2(3,:)= &
                   TROTVEC1(1,:)*TMPM(1,3)+ &
                   TROTVEC1(2,:)*TMPM(2,3)+ &
                   TROTVEC1(3,:)*TMPM(3,3)
#else
! Rotate in direct space
               TROTVEC2(1,:)=TROTVEC1(1,:)*S(1,1,IROT)+ &
                    TROTVEC1(2,:)*S(2,1,IROT)+ &
                    TROTVEC1(3,:)*S(3,1,IROT)

               TROTVEC2(2,:)=TROTVEC1(1,:)*S(1,2,IROT)+ &
                    TROTVEC1(2,:)*S(2,2,IROT)+ &
                    TROTVEC1(3,:)*S(3,2,IROT)

               TROTVEC2(3,:)=TROTVEC1(1,:)*S(1,3,IROT)+ &
                    TROTVEC1(2,:)*S(2,3,IROT)+ &
                    TROTVEC1(3,:)*S(3,3,IROT)
#endif
! bring TROTVEC2 to cartesian coordinates
              CALL DIRKAR(LDIM+1,TROTVEC2,A)
! And back to SAXIS representation
              TROTVEC1(1,:)=COS(BETA)*COS(ALPHA)*TROTVEC2(1,:)+ &
                   COS(BETA)*SIN(ALPHA)*TROTVEC2(2,:)- &
                   SIN(BETA)*TROTVEC2(3,:)
              TROTVEC1(2,:)=-SIN(ALPHA)*TROTVEC2(1,:)+ &
                   COS(ALPHA)*TROTVEC2(2,:)
              TROTVEC1(3,:)=SIN(BETA)*COS(ALPHA)*TROTVEC2(1,:)+ &
                   SIN(BETA)*SIN(ALPHA)*TROTVEC2(2,:)+ &
                   COS(BETA)*TROTVEC2(3,:)

              TMP(:,IA,1)=TMP(:,IA,1)+TROTVEC1(1,:)
              TMP(:,IA,2)=TMP(:,IA,2)+TROTVEC1(2,:)
              TMP(:,IA,3)=TMP(:,IA,3)+TROTVEC1(3,:)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      ! divide final result by the number of translations and rotations
      SCALE=1._q/(NPCELL*NROTK)

      VEC=TMP*SCALE

      DEALLOCATE(TMP,TROTVEC1,TROTVEC2)

      END SUBROUTINE MAGSYL


!******************* SUBROUTINE PARSYLM ********************************
!
! Subroutine PARSYLM symmetrizes l+m+site resolved partial charges
! or the corresponding local magnetic moments (when ISP=2),
! for l=0,..,LDIM, and m=1,2*l+1
!
! These charges/moments must be supplied in the array:
!
!   VEC(1:(LDIM+1)*(LDIM+1),1:NIONS)
!
! On exit this array will contain the symmetrized result.
!
!***********************************************************************

      SUBROUTINE PARSYLM(LDIM,NIONS,NIOND,NTYP,NITYP,VEC,ROTMAP,MAGROT,A,B,ISP)
      USE prec
      USE asa
      IMPLICIT REAL(q) (A-H,O-Z)

      INTEGER :: LDIM,NIONS,NIOND,NTYP,ISP
      INTEGER :: NITYP(NTYP),ROTMAP(NIOND,NROTK,NPCELL)

      REAL(q) :: MAGROT(48,NPCELL),A(3,3),B(3,3)
      REAL(q) :: VEC((LDIM+1)*(LDIM+1),NIONS)

! local variables
      INTEGER :: NROT,LMAX,LMDIM,ITRANS,IA,IAP,MMAX,IROT,IS,ISTART
      REAL(q) :: SCALE

      REAL(q), ALLOCATABLE :: TMP(:,:),ROTVEC(:,:)
      REAL(q), ALLOCATABLE :: SL(:,:,:)

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

!-------------------------------------------------------------------
! allocate work arrays
!-------------------------------------------------------------------
      MMAX=2*LDIM+1
      LMDIM=(LDIM+1)*(LDIM+1)
      ALLOCATE(TMP(LMDIM,NIONS),SL(MMAX,MMAX,0:LDIM),ROTVEC(LMDIM,NIONS))

      TMP=0
!-------------------------------------------------------------------
! do the symmetrization
!-------------------------------------------------------------------
      ISTART=1
      ! loop over all species
      DO IS=1,NTYP
        IF (IS>1) ISTART=ISTART+NITYP(IS-1)
        ! loop over all rotations
        DO IROT=1,NROTK
          ! setup rotation matrices for L=0,...,LMAX
          CALL SETUP_SYM_LL(MMAX,LDIM,ISYMOP(1,1,IROT),SL,A,B)
          ! loop over all ions
          DO IA=ISTART,NITYP(IS)+ISTART-1
            ! rotate the matrix and store result in ROTVEC
            CALL ROTATE_DIAGONAL(LDIM,VEC(1,IA),MMAX,LDIM,SL,ROTVEC(1,IA))
          ENDDO

          ! loop over all ions
          DO IA=ISTART,NITYP(IS)+ISTART-1
            ! loop over all translations
            DO ITRANS=1,NPCELL
              ! destination atom
              IAP=ROTMAP(IA,IROT,ITRANS)
              SCALE=1._q
              IF (ISP==2) SCALE=MAGROT(IROT,ITRANS)
              TMP(:,IA)=TMP(:,IA)+ROTVEC(:,IAP)*SCALE
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      ! divide final result by the number of translations and rotations
      SCALE=1._q/(NPCELL*NROTK)

      VEC=TMP*SCALE

      DEALLOCATE(TMP,ROTVEC,SL)

      END SUBROUTINE PARSYLM


!******************* SUBROUTINE MAGSYLM ********************************
!
! Subroutine MAGSYLM symmetrizes l+m+site resolved noncollinear
! local magnetic moments, for l=0,..,LDIM, and m=1,2*l+1
!
! These moments must be supplied in the array:
!
!   VEC(1:(LDIM+1)*(LDIM+1),1:NIONS,1:3)
!
! On exit this array will contain the symmetrized result.
!
!***********************************************************************

      SUBROUTINE MAGSYLM(LDIM,NIONS,NIOND,NTYP,NITYP,VEC,ROTMAP,MAGROT,SAXIS,A,B)
      USE prec
      USE asa
      USE lattice
      IMPLICIT REAL(q) (A-H,O-Z)

      INTEGER :: LDIM,NIONS,NIOND,NTYP,ISP
      INTEGER :: NITYP(NTYP),ROTMAP(NIOND,NROTK,NPCELL)

      REAL(q) :: MAGROT(48,NPCELL),A(3,3),B(3,3)
      REAL(q) :: VEC((LDIM+1)*(LDIM+1),NIONS,3),SAXIS(3)

! local variables
      INTEGER :: NROT,LMAX,LMDIM,ITRANS,IA,IAP,MMAX,IROT,IS,ISTART,IDIR
      REAL(q) :: ALPHA,BETA,SCALE

      REAL(q), ALLOCATABLE :: TMP(:,:,:),ROTVEC(:,:,:)
      REAL(q), ALLOCATABLE :: TROTVEC1(:,:),TROTVEC2(:,:)
      REAL(q), ALLOCATABLE :: SL(:,:,:)
#ifndef oldsym
      INTEGER DET,TMPM(3,3)
#endif

      COMMON /SYMM/ ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
     &                            GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

      CALL EULER(SAXIS,ALPHA,BETA)

!-------------------------------------------------------------------
! allocate work arrays
!-------------------------------------------------------------------
      MMAX=2*LDIM+1
      LMDIM=(LDIM+1)*(LDIM+1)
      ALLOCATE(TMP(LMDIM,NIONS,3),SL(MMAX,MMAX,0:LDIM),ROTVEC(LMDIM,NIONS,3))
      ALLOCATE(TROTVEC1(3,LMDIM),TROTVEC2(3,LMDIM))

      TMP=0
!-------------------------------------------------------------------
! do the symmetrization
!-------------------------------------------------------------------
      ISTART=1
      ! loop over all species
      DO IS=1,NTYP
        IF (IS>1) ISTART=ISTART+NITYP(IS-1)
        ! loop over all rotations
        DO IROT=1,NROTK
          ! setup rotation matrices for L=0,...,LMAX
          CALL SETUP_SYM_LL(MMAX,LDIM,ISYMOP(1,1,IROT),SL,A,B)
          ! loop over all ions
          DO IA=ISTART,NITYP(IS)+ISTART-1
            ! rotate the matrix and store result in ROTVEC
            DO IDIR=1,3
              CALL ROTATE_DIAGONAL(LDIM,VEC(1,IA,IDIR),MMAX,LDIM,SL,ROTVEC(1,IA,IDIR))
            ENDDO
          ENDDO

          ! loop over all ions
          DO IA=ISTART,NITYP(IS)+ISTART-1
            ! loop over all translations
            DO ITRANS=1,NPCELL
              ! destination atom
              IAP=ROTMAP(IA,IROT,ITRANS)
              SCALE=1._q

! Transform from "SAXIS basis" to the system of cartesian axes
              TROTVEC1(1,:)=COS(BETA)*COS(ALPHA)*ROTVEC(:,IAP,1)- &
                   SIN(ALPHA)*ROTVEC(:,IAP,2)+ &
                   SIN(BETA)*COS(ALPHA)*ROTVEC(:,IAP,3)
              TROTVEC1(2,:)=COS(BETA)*SIN(ALPHA)*ROTVEC(:,IAP,1)+ &
                   COS(ALPHA)*ROTVEC(:,IAP,2)+ &
                   SIN(BETA)*SIN(ALPHA)*ROTVEC(:,IAP,3)
              TROTVEC1(3,:)=-SIN(BETA)*ROTVEC(:,IAP,1)+ &
                   COS(BETA)*ROTVEC(:,IAP,3)
! Bring to direct coordinates in which the integer rotation matrices are defined
              CALL KARDIR(LMDIM,TROTVEC1,B)
#ifndef oldsym
              TMPM=ISYMOP(:,:,IROT)
              DET=TMPM(1,1)*TMPM(2,2)*TMPM(3,3) - &
                  TMPM(1,1)*TMPM(2,3)*TMPM(3,2) + &
                  TMPM(1,2)*TMPM(2,3)*TMPM(3,1) - &
                  TMPM(1,2)*TMPM(2,1)*TMPM(3,3) + &
                  TMPM(1,3)*TMPM(2,1)*TMPM(3,2) - &
                  TMPM(1,3)*TMPM(2,2)*TMPM(3,1)
              IF (DET<0) TMPM=-TMPM

              TROTVEC2(1,:)= &
                   TROTVEC1(1,:)*TMPM(1,1)+ &
                   TROTVEC1(2,:)*TMPM(2,1)+ &
                   TROTVEC1(3,:)*TMPM(3,1)

              TROTVEC2(2,:)= &
                   TROTVEC1(1,:)*TMPM(1,2)+ &
                   TROTVEC1(2,:)*TMPM(2,2)+ &
                   TROTVEC1(3,:)*TMPM(3,2)

              TROTVEC2(3,:)= &
                   TROTVEC1(1,:)*TMPM(1,3)+ &
                   TROTVEC1(2,:)*TMPM(2,3)+ &
                   TROTVEC1(3,:)*TMPM(3,3)
#else
! Rotate in direct space
               TROTVEC2(1,:)=TROTVEC1(1,:)*S(1,1,IROT)+ &
                    TROTVEC1(2,:)*S(2,1,IROT)+ &
                    TROTVEC1(3,:)*S(3,1,IROT)

               TROTVEC2(2,:)=TROTVEC1(1,:)*S(1,2,IROT)+ &
                    TROTVEC1(2,:)*S(2,2,IROT)+ &
                    TROTVEC1(3,:)*S(3,2,IROT)

               TROTVEC2(3,:)=TROTVEC1(1,:)*S(1,3,IROT)+ &
                    TROTVEC1(2,:)*S(2,3,IROT)+ &
                    TROTVEC1(3,:)*S(3,3,IROT)
#endif
! bring TROTVEC2 to cartesian coordinates
              CALL DIRKAR(LMDIM,TROTVEC2,A)
! And back to SAXIS representation
              TROTVEC1(1,:)=COS(BETA)*COS(ALPHA)*TROTVEC2(1,:)+ &
                   COS(BETA)*SIN(ALPHA)*TROTVEC2(2,:)- &
                   SIN(BETA)*TROTVEC2(3,:)
              TROTVEC1(2,:)=-SIN(ALPHA)*TROTVEC2(1,:)+ &
                   COS(ALPHA)*TROTVEC2(2,:)
              TROTVEC1(3,:)=SIN(BETA)*COS(ALPHA)*TROTVEC2(1,:)+ &
                   SIN(BETA)*SIN(ALPHA)*TROTVEC2(2,:)+ &
                   COS(BETA)*TROTVEC2(3,:)

              TMP(:,IA,1)=TMP(:,IA,1)+TROTVEC1(1,:)
              TMP(:,IA,2)=TMP(:,IA,2)+TROTVEC1(2,:)
              TMP(:,IA,3)=TMP(:,IA,3)+TROTVEC1(3,:)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      ! divide final result by the number of translations and rotations
      SCALE=1._q/(NPCELL*NROTK)

      VEC=TMP*SCALE

      DEALLOCATE(TMP,ROTVEC,SL,TROTVEC1,TROTVEC2)

      END SUBROUTINE MAGSYLM


!******************* SUBROUTINE ROTATE_DIAGONAL ************************
!
! Subroutine ROTATE_DIAGONAL rotates a vector of lm resolved partial
! charges (or magnetic moments), for l=0,LDIM and m=1,2l+1.
! The lm resolved partial charge/moments are supplied in:
!
!   VEC(1:(LDIM+1)*(LDIM+1))
!
! The storage layout of VEC is the standard one:
!
!    s, py, pz, px, dxy, dyz, dz2, dxz, dx2-y2  (.. etc ..)
!
! (for a complete list of components see for instance SPHPRO_FAST)
!
! The rotation matrix must be supplied in SL and has to be set up
! outside by a call to SETUP_SYM_LL.
!
! The result is returned in ROTVEC.
!
!***********************************************************************

      SUBROUTINE ROTATE_DIAGONAL(LDIM,VEC,MMAX,LMAX,SL,ROTVEC)
      USE prec
      REAL(q) :: VEC((LDIM+1)*(LDIM+1)),ROTVEC((LDIM+1)*(LDIM+1))
      REAL(q) :: SL(MMAX,MMAX,0:LMAX)
      INTEGER :: LDIM,MMAX,LMAX
! local variables
      REAL(q) :: TMP((LDIM+1)*(LDIM+1),(LDIM+1)*(LDIM+1))
      INTEGER :: IND,L,M,MP

      IND=0
      ROTVEC=0

      DO L=0,LDIM
        TMP=0
        ! left hand transformation
        DO M=1,(2*L+1)
        DO MP=1,(2*L+1)
          TMP(M,MP)=TMP(M,MP)+SL(M,MP,L)*VEC(IND+MP)
        ENDDO
        ENDDO
        ! right hand transformation
        DO M=1,(2*L+1)
        DO MP=1,(2*L+1)
          ROTVEC(IND+M)=ROTVEC(IND+M)+SL(M,MP,L)*TMP(M,MP)
        ENDDO
        ENDDO
        IND=IND+(2*L+1)
      ENDDO
      END SUBROUTINE ROTATE_DIAGONAL
