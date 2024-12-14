#include "symbol.inc"
MODULE mathtools
   USE prec
   USE string, ONLY: str
   USE tutor,  ONLY: vtutor
   IMPLICIT NONE 

   PRIVATE

   PUBLIC INVERT_REAL_MATRIX, SV_DECOMPOSER, RRQR_DECOMPOSER


!                                _______________
!_______________________________/ SV_DECOMPOSER \_______________________________
!
!> @brief Wrapper for DGESDD/ZGESDD LAPACK routines
!>
!> This type facilitates the functionality for performing a singular-value
!> decomposition multiple times without reallocating any of the buffers or work
!> arrays (for example inside a loop).
!>
!> The public type variables are used as storage for the output arrays.
!>
!> @var U_MAT           Left matrix in the SVD
!> @var V_DAGGER_MAT    Right matrix in the SVD
!> @var SINGULAR_VALS   Singular values
!_______________________________________________________________________________
   TYPE SV_DECOMPOSER
      GDEF, ALLOCATABLE,      PUBLIC   :: U_MAT(:, :)
      GDEF, ALLOCATABLE,      PUBLIC   :: V_DAGGER_MAT(:, :)
      REAL(q), ALLOCATABLE,   PUBLIC   :: SINGULAR_VALS(:)

      CHARACTER(1),  PRIVATE  :: MODE = ""
      INTEGER,       PRIVATE  :: DIM_1 = 0
      INTEGER,       PRIVATE  :: DIM_2 = 0
      INTEGER,       PRIVATE  :: MIN_DIM = 0
      INTEGER,       PRIVATE  :: LDU = 0
      INTEGER,       PRIVATE  :: LDVD = 0

      INTEGER,                PRIVATE  :: WORK_SIZE
      INTEGER, ALLOCATABLE,   PRIVATE  :: IWORK(:)
      GDEF, ALLOCATABLE,      PRIVATE  :: WORK(:)
      REAL(q), ALLOCATABLE,   PRIVATE  :: RWORK(:)

      CONTAINS

      PROCEDURE,  PUBLIC   :: INIT => SV_DECOMPOSER_INIT
      PROCEDURE,  PUBLIC   :: CALC => SV_DECOMPOSER_CALC
      PROCEDURE,  PUBLIC   :: DEALLOC => SV_DECOMPOSER_DEALLOC
   END TYPE


!                               _________________
!______________________________/ RRQR_DECOMPOSER \______________________________
!
!> @brief Wrapper for xGEQP3 LAPACK routines
!>
!> This type facilitates the functionality for performing a rank-revealing QR
!> decomposition multiple times without reallocating any of the buffers or work
!> arrays (for example inside a loop).
!>
!> The public type variables are used as storage for the output arrays.
!>
!> @var PIVOT  Pivot values
!> @var TAU    Scalar factors of the elementary reflectors
!_______________________________________________________________________________
   TYPE RRQR_DECOMPOSER
      INTEGER, ALLOCATABLE,   PUBLIC   :: PIVOT(:)
      GDEF, ALLOCATABLE,      PUBLIC   :: TAU(:)

      INTEGER,       PRIVATE  :: DIM_1 = 0
      INTEGER,       PRIVATE  :: DIM_2 = 0
      INTEGER,       PRIVATE  :: LDA = 0

      INTEGER,                PRIVATE  :: WORK_SIZE
      GDEF, ALLOCATABLE,      PRIVATE  :: WORK(:)
      REAL(q), ALLOCATABLE,   PRIVATE  :: RWORK(:)

      CONTAINS

      PROCEDURE,  PUBLIC   :: INIT => RRQR_DECOMPOSER_INIT
      PROCEDURE,  PUBLIC   :: CALC => RRQR_DECOMPOSER_CALC
      PROCEDURE,  PUBLIC   :: DEALLOC => RRQR_DECOMPOSER_DEALLOC
   END TYPE


   CONTAINS

!                             ____________________
!____________________________/ SV_DECOMPOSER_INIT \_____________________________
!
!> @brief Class constructor
!>
!> Internal buffers and work arrays are set up based on the supplied matrix
!> dimensions of a matrix A that is about to be singular-value decomposed by
!> this calculator.
!>
!> @param[in,out] THIS  Instance reference
!> @param[in] MODE      Same as `MODE` argument to DGESDD or ZGESDD
!> @param[in] DIM_1     SIZE(A, 1)
!> @param[in] DIM_2     SIZE(A, 2)
!> @param[out] INFO     `INFO` parameter of DGESDD/ZGESDD query
!_______________________________________________________________________________
   SUBROUTINE SV_DECOMPOSER_INIT(THIS, MODE, DIM_1, DIM_2, INFO)
      USE string, ONLY: uppercase
      IMPLICIT NONE

      CLASS(SV_DECOMPOSER),   INTENT(INOUT)  :: THIS
      CHARACTER(1),           INTENT(IN)     :: MODE
      INTEGER,                INTENT(IN)     :: DIM_1
      INTEGER,                INTENT(IN)     :: DIM_2
      INTEGER,                INTENT(OUT)    :: INFO

      GDEF  :: ASK_DIM
      GDEF  :: DUMMY_G

      CALL THIS%DEALLOC()

      THIS%MODE = uppercase(MODE)
      THIS%DIM_1 = DIM_1
      THIS%DIM_2 = DIM_2
      THIS%MIN_DIM = MIN(DIM_1, DIM_2)

      SELECT CASE (THIS%MODE)
      CASE("A")
         ALLOCATE(THIS%U_MAT(DIM_1, DIM_1), THIS%V_DAGGER_MAT(DIM_2, DIM_2))
         THIS%LDU = DIM_1
         THIS%LDVD = DIM_2
      CASE("S")
         ALLOCATE(THIS%U_MAT(DIM_1, THIS%MIN_DIM), THIS%V_DAGGER_MAT(THIS%MIN_DIM, DIM_2))
         THIS%LDU = DIM_1
         THIS%LDVD = THIS%MIN_DIM
      CASE("O")
         IF (DIM_1 >= DIM_2) THEN
            ALLOCATE(THIS%U_MAT(1, 1), THIS%V_DAGGER_MAT(DIM_2, DIM_2))
            THIS%LDU = 1
            THIS%LDVD = DIM_2
         ELSE
            ALLOCATE(THIS%U_MAT(DIM_1, DIM_1), THIS%V_DAGGER_MAT(1, 1))
            THIS%LDU = DIM_1
            THIS%LDVD = 1
         ENDIF
      CASE("N")
         ALLOCATE(THIS%U_MAT(1, 1), THIS%V_DAGGER_MAT(1, 1))
         THIS%LDU = 1
         THIS%LDVD = 1
      CASE DEFAULT
         CALL vtutor%bug("SV_DECOMPOSER_INIT: Mode '" // THIS%MODE // "' not implemented!", __FILE__, __LINE__)
      END SELECT

      ALLOCATE(THIS%SINGULAR_VALS(THIS%MIN_DIM))

      ! Allocate work arrays
      ALLOCATE(THIS%IWORK(8 * THIS%MIN_DIM))

      ! Query optimal size of work array
#ifdef gammareal
      CALL DGESDD(MODE, DIM_1, DIM_2, DUMMY_G, DIM_1, THIS%SINGULAR_VALS, &
                  THIS%U_MAT, THIS%LDU, THIS%V_DAGGER_MAT, THIS%LDVD, &
                  ASK_DIM, -1, THIS%IWORK, INFO)
#else
      ALLOCATE(THIS%RWORK(THIS%MIN_DIM * MAX(5 * (THIS%MIN_DIM + 1), 2 * (MAX(DIM_1, DIM_2) + THIS%MIN_DIM) + 1)))
      CALL ZGESDD(MODE, DIM_1, DIM_2, DUMMY_G, DIM_1, THIS%SINGULAR_VALS, &
                  THIS%U_MAT, THIS%LDU, THIS%V_DAGGER_MAT, THIS%LDVD, &
                  ASK_DIM, -1, THIS%RWORK, THIS%IWORK, INFO)
#endif

      IF (INFO /= 0) THEN
         CALL THIS%DEALLOC()
         RETURN
      ENDIF

      THIS%WORK_SIZE = INT(ASK_DIM)
      ALLOCATE(THIS%WORK(THIS%WORK_SIZE))
   END SUBROUTINE


!                             ____________________
!____________________________/ SV_DECOMPOSER_CALC \_____________________________
!
!> @brief Performs singular-value decomposition of matrix A.
!>
!> The return values are either returned in A or are stored in the public type
!> variables SINGULAR_VALS, U_MAT or V_DAGGER_MAT, in accordance with the
!> chosen mode. This is consistent with how the output is defined in DGESDD and
!> ZGESDD.
!>
!> @param[in,out] THIS  Instance reference
!> @param[in,out] A     Matrix A
!> @param[out] INFO     `INFO` parameter of DGESDD/ZGESDD
!_______________________________________________________________________________
   SUBROUTINE SV_DECOMPOSER_CALC(THIS, A, INFO)
      IMPLICIT NONE
   
      CLASS(SV_DECOMPOSER),   INTENT(INOUT)  :: THIS
      GDEF,                   INTENT(INOUT)  :: A(:, :) !! Content is destroyed on exit
      INTEGER,                INTENT(OUT)    :: INFO
      
      ! Perform singular value decomposition of A
#ifdef gammareal
      CALL DGESDD(THIS%MODE, THIS%DIM_1, THIS%DIM_2, A, THIS%DIM_1, THIS%SINGULAR_VALS, &
                  THIS%U_MAT, THIS%LDU, THIS%V_DAGGER_MAT, THIS%LDVD, &
                  THIS%WORK, THIS%WORK_SIZE, THIS%IWORK, INFO)
#else
      CALL ZGESDD(THIS%MODE, THIS%DIM_1, THIS%DIM_2, A, THIS%DIM_1, THIS%SINGULAR_VALS, &
                  THIS%U_MAT, THIS%LDU, THIS%V_DAGGER_MAT, THIS%LDVD, &
                  THIS%WORK, THIS%WORK_SIZE, THIS%RWORK, THIS%IWORK, INFO)
#endif
   END SUBROUTINE


!                            _______________________
!___________________________/ SV_DECOMPOSER_DEALLOC \___________________________
!
!> @brief Deallocates the buffer and work arrays.
!>
!> This can be used in case the allocated arrays are too large in between
!> operations. However, the object needs to be reinitialized for the next
!> operation. Deallocation happens at the end of scope regardless, so normally,
!> there shouldn't be a need to call this.
!>
!> @param[in,out] THIS  Instance reference
!_______________________________________________________________________________
   SUBROUTINE SV_DECOMPOSER_DEALLOC(THIS)
      IMPLICIT NONE
      CLASS(SV_DECOMPOSER),   INTENT(INOUT)  :: THIS

      IF (ALLOCATED(THIS%U_MAT)) DEALLOCATE(THIS%U_MAT)
      IF (ALLOCATED(THIS%V_DAGGER_MAT)) DEALLOCATE(THIS%V_DAGGER_MAT)
      IF (ALLOCATED(THIS%SINGULAR_VALS)) DEALLOCATE(THIS%SINGULAR_VALS)
      IF (ALLOCATED(THIS%IWORK)) DEALLOCATE(THIS%IWORK)
      IF (ALLOCATED(THIS%WORK)) DEALLOCATE(THIS%WORK)
      IF (ALLOCATED(THIS%RWORK)) DEALLOCATE(THIS%RWORK)
   END SUBROUTINE


!                            ______________________
!___________________________/ RRQR_DECOMPOSER_INIT \____________________________
!
!> @brief @brief Class constructor
!>
!> Internal buffers and work arrays are set up based on the supplied matrix
!> dimensions of a matrix A that is about to be QR decomposed by this
!> calculator.
!>
!> @param[in,out] THIS  Instance reference
!> @param[in] DIM_1     SIZE(A, 1)
!> @param[in] DIM_2     SIZE(A, 2)
!> @param[out] INFO     `INFO` parameter of xGEQP3 query
!_______________________________________________________________________________
   SUBROUTINE RRQR_DECOMPOSER_INIT(THIS, DIM_1, DIM_2, INFO)
      IMPLICIT NONE
      
      CLASS(RRQR_DECOMPOSER), INTENT(INOUT)  :: THIS
      INTEGER,                INTENT(IN)     :: DIM_1
      INTEGER,                INTENT(IN)     :: DIM_2
      INTEGER,                INTENT(OUT)    :: INFO

      GDEF  :: ASK_DIM
      GDEF  :: DUMMY_G

      CALL THIS%DEALLOC()

      THIS%DIM_1 = DIM_1
      THIS%DIM_2 = DIM_2

      ALLOCATE(THIS%PIVOT(DIM_2))
      THIS%PIVOT = 0

      ALLOCATE(THIS%TAU(MIN(DIM_1, DIM_2)))

      ! Query optimal size of work arrays
#ifdef gammareal
      CALL DGEQP3(DIM_1, DIM_2, DUMMY_G, DIM_1, THIS%PIVOT, THIS%TAU, ASK_DIM, -1, INFO)
#else
      ALLOCATE(THIS%RWORK(2 * DIM_2))
      CALL ZGEQP3(DIM_1, DIM_2, DUMMY_G, DIM_1, THIS%PIVOT, THIS%TAU, ASK_DIM, -1, THIS%RWORK, INFO)
#endif

      IF (INFO /= 0) THEN
         CALL THIS%DEALLOC()
         RETURN
      ENDIF

      THIS%WORK_SIZE = INT(ASK_DIM)
      ALLOCATE(THIS%WORK(THIS%WORK_SIZE))
   END SUBROUTINE


!                            ______________________
!___________________________/ RRQR_DECOMPOSER_CALC \____________________________
!
!> @brief Performs rank-revealing QR decomposition of matrix A.
!>
!> The return values are either returned in A or are stored in the public type
!> variables PIVOT and TAU. This is consistent with how the output is defined
!> in xGEQP3.
!>
!> @param[in,out] THIS  Instance reference
!> @param[in,out] A     Matrix A
!> @param[out] INFO     `INFO` parameter of xGEQP3
!_______________________________________________________________________________
   SUBROUTINE RRQR_DECOMPOSER_CALC(THIS, A, INFO)
      IMPLICIT NONE

      CLASS(RRQR_DECOMPOSER), INTENT(INOUT)  :: THIS
      GDEF,                   INTENT(INOUT)  :: A(:, :) !! Content is destroyed on exit
      INTEGER,                INTENT(OUT)    :: INFO

      ! Perform rank-revealing QR decomposition
#ifdef gammareal
      CALL DGEQP3(THIS%DIM_1, THIS%DIM_2, A, THIS%DIM_1, THIS%PIVOT, THIS%TAU, &
                  THIS%WORK, THIS%WORK_SIZE, INFO)
#else
      CALL ZGEQP3(THIS%DIM_1, THIS%DIM_2, A, THIS%DIM_1, THIS%PIVOT, THIS%TAU, &
                  THIS%WORK, THIS%WORK_SIZE, THIS%RWORK, INFO)
#endif
   END SUBROUTINE


!                           _________________________
!__________________________/ RRQR_DECOMPOSER_DEALLOC \__________________________
!
!> @brief Deallocates the buffer and work arrays.
!>
!> This can be used in case the allocated arrays are too large in between
!> operations. However, the object needs to be reinitialized for the next
!> operation. Deallocation happens at the end of scope regardless, so normally,
!> there shouldn't be a need to call this.
!>
!> @param[in,out] THIS
!_______________________________________________________________________________
   SUBROUTINE RRQR_DECOMPOSER_DEALLOC(THIS)
      IMPLICIT NONE
      CLASS(RRQR_DECOMPOSER), INTENT(INOUT)  :: THIS

      IF (ALLOCATED(THIS%PIVOT)) DEALLOCATE(THIS%PIVOT)
      IF (ALLOCATED(THIS%TAU)) DEALLOCATE(THIS%TAU)
      IF (ALLOCATED(THIS%WORK)) DEALLOCATE(THIS%WORK)
      IF (ALLOCATED(THIS%RWORK)) DEALLOCATE(THIS%RWORK)
   END SUBROUTINE


!____________________________________________________________________________
!!inverts a real matrix using LU decomposition
!____________________________________________________________________________
 
   SUBROUTINE INVERT_REAL_MATRIX( MATRIX, IU0 )
      REAL(q)      :: MATRIX(:,:)
      INTEGER      :: IU0
      !local 
      INTEGER             :: N
      REAL(q),ALLOCATABLE :: A(:,:)
      INTEGER,ALLOCATABLE :: IPIV(:)
      REAL(q),ALLOCATABLE :: WORK(:)
      INTEGER             :: IWORK
      INTEGER             :: INFO

      PROFILING_START('invert_real_matrix')

      !obtain size of matrix
      N=SIZE(MATRIX,1)

      !consistency check
      IF( N/= SIZE(MATRIX,2)) THEN
         CALL vtutor%bug("internal error in VASP: INVERT_REAL_MATRIX reports matrix not square " // &
            str(N) // " " // str(SIZE(MATRIX,2)), __FILE__, __LINE__)
      ENDIF

      !store MATRIX to A
      ALLOCATE(A(N,N))
      A(1:N,1:N)=MATRIX(1:N,1:N)

      !allocation of IPIV
      ALLOCATE(IPIV(N))

      !LU decomposition
      CALL DGETRF(N,N,A(1,1),N,IPIV,INFO)
      IF ( INFO/=0) THEN
         CALL vtutor%bug("internal error in VASP: INVERT_REAL_MATRIX reports " // str(INFO), __FILE__, __LINE__)
      ENDIF

      !obtain optimal working array sizes
      ALLOCATE(WORK(1))
      CALL DGETRI(N,A(1,1),N,IPIV,WORK,-1,INFO)
      IWORK=WORK(1)
      IF ( IWORK <= 0 ) THEN
         CALL vtutor%bug("internal error in VASP: INVERT_REAL_MATRIX reports illegal size from " &
            // "DGETRI " // str(IWORK), __FILE__, __LINE__)
      ENDIF
      !allocate working array with optimal size
      DEALLOCATE(WORK)
      ALLOCATE(WORK(IWORK))

      !actual inversion comes now
      CALL DGETRI(N,A(1,1),N,IPIV,WORK,IWORK,INFO)
      IF ( INFO/=0) THEN
         CALL vtutor%bug("internal error in VASP: INVERT_REAL_MATRIX reports DGETRI reports: " // &
            str(INFO), __FILE__, __LINE__)
      ENDIF

      !store result into MATRIX, deallocate and return
      MATRIX(1:N,1:N)=A(1:N,1:N)

      PROFILING_STOP('invert_real_matrix')
   END SUBROUTINE INVERT_REAL_MATRIX

END MODULE mathtools
