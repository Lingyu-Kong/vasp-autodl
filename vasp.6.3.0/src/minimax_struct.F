#include "symbol.inc"
MODULE minimax_struct
   USE prec
   USE base, ONLY: TOREAL
#ifdef qd_emulate
   USE qdmodule
#endif
   USE mpimy
   IMPLICIT NONE
!
!! general descriptor for a quadrature 
!
   TYPE quadrature_handle
      !! iidentification number 
      INTEGER          :: ITYPE=0
      !! grid id 
      INTEGER          :: GRID_ID=0
      !! number of quadrature points     
      INTEGER          :: N=0
      !! scaled interval [1,B] (or [0,B] ) for which quadrature error is minimized
      QDPREAL          :: A, B
      !! scaling factor
      QDPREAL          :: SCALING
      !! L_infty norm or max error of quadature
      QDPREAL          :: INFINITYNORM
      !! quadrature coefficients
      QDPREAL, POINTER :: C(:) => NULL()
      !! alternant of error function
      QDPREAL, POINTER :: X0(:) => NULL()
     
      !! function used to fit object function
      PROCEDURE( basis_function ), POINTER, NOPASS :: PHI2 => NULL()
      !! first derivative of basis function used for fit w.r.t. first argument
      PROCEDURE( basis_function ), POINTER, NOPASS :: D_PHI2_DZ => NULL()
      !! first derivative of basis function used for fit w.r.t. second argument
      PROCEDURE( basis_function ), POINTER, NOPASS :: D_PHI2_DL => NULL()
      !! second derivative of basis function used for fit w.r.t. first argument
      PROCEDURE( basis_function ), POINTER, NOPASS :: D2_PHI2_DZ2 => NULL()

      !! oject function and derivatives
      PROCEDURE( object_function ), POINTER, NOPASS :: F => NULL()
      PROCEDURE( object_function ), POINTER, NOPASS :: D_F_DZ => NULL()
      PROCEDURE( object_function ), POINTER, NOPASS :: D2_F_DZ2 => NULL()

      !! subroutines that set initial quadrature values
      PROCEDURE( init_quad_subroutine ), POINTER, NOPASS :: R1E3 => NULL()
      PROCEDURE( init_quad_subroutine ), POINTER, NOPASS :: R1E4 => NULL()
      PROCEDURE( init_quad_subroutine ), POINTER, NOPASS :: R1E5 => NULL()
      PROCEDURE( init_quad_subroutine ), POINTER, NOPASS :: R1E6 => NULL()
      !! tabulated coefficients in minimax_ini
      QDPREAL, POINTER :: CTAB(:) => NULL()

      !! basis used to fit function and derivatives
      PROCEDURE( basis_function ), POINTER, NOPASS :: PHI => NULL()
      !! dual basis 
      PROCEDURE( basis_function ), POINTER, NOPASS :: PSI => NULL()
      !! auxilary function used for transformation matrix for low frequency point
      PROCEDURE( object_function ), POINTER, NOPASS :: TRANS_TAYLOR => NULL()

      ! conjugate transformation functions
      !! basis used to fit function and derivatives
      PROCEDURE( basis_function ), POINTER, NOPASS :: PHI_CONJG => NULL()
      !! dual basis 
      PROCEDURE( basis_function ), POINTER, NOPASS :: PSI_CONJG => NULL()
      !! auxilary function used for transformation matrix for low frequency point
      PROCEDURE( object_function ), POINTER, NOPASS :: TRANS_TAYLOR_CONJG => NULL()

      !! vectorial function that fits the minimax coefficients 
      PROCEDURE( fit_coeff_function ), POINTER, NOPASS :: FITCOF => NULL()
      QDPREAL          :: R1 ! interval, where fit works
      QDPREAL          :: R2 !

      !! communicator used in internal minimax routines
      TYPE( communic ) :: COMM

   END TYPE quadrature_handle
!
!  unified imaginary grid descriptor for time and frequency domain
!  contains communicators and other local as well as global indices
!
   TYPE loop_des
      INTEGER :: NPOINTS                           ! total number of time/freq points tau
      INTEGER :: NPOINTSC                          ! loop counter, indexing the current grid  point handled
      INTEGER :: NPOINTS_IN_GROUP                  ! number of points in group
      INTEGER :: NPOINTS_IN_ROOT_GROUP             ! number of points in first group (always the largest)
      REAL(q) :: POINT_CURRENT                     ! imaginary time treated currently (only one at each step)
      REAL(q) :: POINTMIN_CURRENT                  ! current smallest tau 
      LOGICAL :: LDO_POINT_LOCAL                   ! skip GEMM calls?  
 
      TYPE (communic) :: COMM_BETWEEN_GROUPS       ! communicator between tau
                                                   ! processor groups
      TYPE (communic) :: COMM_IN_GROUP             ! communicator witin tau
                                                   ! processor group
      TYPE (communic) :: COMM                      ! global communicator
 
      REAL(q), POINTER :: POINTS_LOCAL(:)=>NULL()  ! local points treated by one processor group

      INTEGER, POINTER :: DISTRIBUTION(:,:)        ! information how points are distributed

      INTEGER, POINTER :: GROUP_OF_NODE(:)         ! processor group id
                                                   ! based on global id
      !inverse temperature (in eV), as default take the one for 300 K
      REAL(q) :: BETA = 38.68149168875351
   END TYPE loop_des
   
   ! define imaginary grid object
   TYPE imag_grid_handle
      ! number of grid points (NOMEGA)
      INTEGER          :: NOMEGA=0
      
      LOGICAL          :: LFINITE_TEMPERATURE = .FALSE.
      !inverse temperature (in eV), as default take the one for 300 K
      REAL(q)          :: BETA = 38.68149168875351_q
      !minimization interval (left boundary)
      REAL(q)          :: X1=0
      !minimization interval (right boundary)
      REAL(q)          :: X2=0
      !minimization interval (scaled right boundary)
      REAL(q)          :: R=0

      !this is used for CRPA 
      REAL(q) :: FACTOR = 1._q

      !integration error for Galitskii-Migdal energy
      REAL(q)          :: GM_INT_ERROR=0
      
      ! imaginary time points 
      REAL(q), POINTER :: TAU(:) => NULL( )
      ! imaginary time weights points 
      REAL(q), POINTER :: TAU_WEIGHT(:)=> NULL( )
      ! imaginary time integration error 
      REAL(q)          :: TAU_ERROR=0
      ! transformation weights from FER_RE to TAU=0+ point 
      REAL(q), POINTER :: TO_TAU0(:,:) => NULL( )
      REAL(q), POINTER :: TO_TAU0_ERROR(:) => NULL( )
      REAL(q), POINTER :: TO_TAU0_CONJG(:,:) => NULL( )
      REAL(q), POINTER :: TO_TAU0_CONJG_ERROR(:) => NULL( )
      ! transformation weights to TAU=0- are same as above with minus sign

      ! imaginary frequency points for real part of bosonic functions
      REAL(q), POINTER :: BOS_RE(:) => NULL( )
      ! weights for frequency integration of real part of bosonic functions
      REAL(q), POINTER :: BOS_RE_WEIGHT(:)=> NULL( )
      ! corresponding error 
      REAL(q)          :: BOS_RE_ERROR=0
      ! transformation matrix from TAU to BOS_RE grid
      REAL(q), POINTER :: TO_BOS_RE(:,:) => NULL( )
      ! error of transformation from TAU to BOS_RE grid points
      REAL(q), POINTER :: TO_BOS_RE_ERROR(:) => NULL( )
      ! transformation matrix from TAU to BOS_RE grid ( conjugate basis )
      REAL(q), POINTER :: TO_BOS_RE_CONJG(:,:) => NULL( )
      ! error of transformation from TAU to BOS_RE grid points
      REAL(q), POINTER :: TO_BOS_RE_CONJG_ERROR(:) => NULL( )

      ! imaginary frequency points for imaginary part of bosonic functions
      REAL(q), POINTER :: BOS_IM(:) => NULL( )
      ! weights for frequency integration of imaginary part of bosonic functions
      REAL(q), POINTER :: BOS_IM_WEIGHT(:)=> NULL( )
      ! corresponding error 
      REAL(q)          :: BOS_IM_ERROR=0
      ! transformation matrix from TAU to BOS_RE grid
      REAL(q), POINTER :: TO_BOS_IM(:,:) => NULL( )
      ! error of transformation from TAU to BOS_RE grid points
      REAL(q), POINTER :: TO_BOS_IM_ERROR(:) => NULL( )
      ! transformation matrix from TAU to BOS_RE grid (conjugate basis )
      REAL(q), POINTER :: TO_BOS_IM_CONJG(:,:) => NULL( )
      ! error of transformation from TAU to BOS_RE grid points
      REAL(q), POINTER :: TO_BOS_IM_CONJG_ERROR(:) => NULL( )

      ! imaginary frequency points for real part of fermionic functions
      REAL(q), POINTER :: FER_RE(:) => NULL( )
      ! weights for frequency integration of real part of fermionic functions
      REAL(q), POINTER :: FER_RE_WEIGHT(:)=> NULL( )
      ! corresponding error 
      REAL(q)          :: FER_RE_ERROR=0
      ! transformation matrix from TAU to FER_RE grid
      REAL(q), POINTER :: TO_FER_RE(:,:) => NULL( )
      ! error of transformation from TAU to FER_RE grid points
      REAL(q), POINTER :: TO_FER_RE_ERROR(:) => NULL( )
      ! transformation matrix from TAU to FER_RE grid (conjugate basis)
      REAL(q), POINTER :: TO_FER_RE_CONJG(:,:) => NULL( )
      ! error of transformation from TAU to FER_RE grid points
      REAL(q), POINTER :: TO_FER_RE_CONJG_ERROR(:) => NULL( )

      ! imaginary frequency points for imaginary part of fermionic functions
      REAL(q), POINTER :: FER_IM(:) => NULL( )
      ! weights for frequency integration of imaginary part of fermionic functions
      REAL(q), POINTER :: FER_IM_WEIGHT(:)=> NULL( )
      ! corresponding error 
      REAL(q)          :: FER_IM_ERROR=0
      ! transformation matrix from TAU to FER_RE grid
      REAL(q), POINTER :: TO_FER_IM(:,:) => NULL( )
      ! error of transformation from TAU to FER_RE grid points
      REAL(q), POINTER :: TO_FER_IM_ERROR(:) => NULL( )
      ! transformation matrix from TAU to FER_RE grid (conjugate basis)
      REAL(q), POINTER :: TO_FER_IM_CONJG(:,:) => NULL( )
      ! error of transformation from TAU to FER_RE grid points
      REAL(q), POINTER :: TO_FER_IM_CONJG_ERROR(:) => NULL( )
 
      ! quadrature handles, 5 distinct grids are possible 
      TYPE( quadrature_handle ) :: TIME
      ! grid for real part of bosonic correlation functions
      TYPE( quadrature_handle ) :: FREQ_BOS_RE
      ! grid for imaginary part of bosonic correlation functions
      TYPE( quadrature_handle ) :: FREQ_BOS_IM
      ! grid for real part of fermionic correlation functions
      TYPE( quadrature_handle ) :: FREQ_FER_RE
      ! grid for imaginary part of fermionic correlation functions
      TYPE( quadrature_handle ) :: FREQ_FER_IM

      ! loop descriptors for time grid
      TYPE(loop_des) :: T
      ! loop descriptors for bosonic grid
      TYPE(loop_des) :: B
      ! loop descriptors for fermionic grid
      TYPE(loop_des) :: F

      !! global communicator
      TYPE( communic ) :: COMM
      !! inter communicator between quadratures
      TYPE( communic ) :: COMM_INTER

   END TYPE imag_grid_handle


   ABSTRACT INTERFACE 
       ! interface for basis function call 
      FUNCTION OBJECT_FUNCTION( X )
         USE prec
#ifdef qd_emulate
         USE qdmodule
#endif
         IMPLICIT NONE
         QDPREAL,INTENT(IN)    :: X         ! argument X of error
         QDPREAL               :: OBJECT_FUNCTION

      END FUNCTION OBJECT_FUNCTION

       ! interface for basis function call 
      FUNCTION BASIS_FUNCTION( X, L )
         USE prec
#ifdef qd_emulate
         USE qdmodule
#endif
         IMPLICIT NONE
         QDPREAL,INTENT(IN)    :: X         ! argument X of error
         QDPREAL,INTENT(IN)    :: L         ! parameter
         QDPREAL               :: BASIS_FUNCTION

      END FUNCTION BASIS_FUNCTION

       ! interface for initialization of quadrature coefficients 
      SUBROUTINE INIT_QUAD_SUBROUTINE()
         USE prec
#ifdef qd_emulate
         USE qdmodule
#endif
         IMPLICIT NONE

      END SUBROUTINE INIT_QUAD_SUBROUTINE

       ! interface for function that calculates MM coefficients directly
      FUNCTION FIT_COEFF_FUNCTION( N, R )
         USE prec
         IMPLICIT NONE
         INTEGER,INTENT(IN)    :: N         ! fit order 
         REAL(q),INTENT(IN)    :: R         ! Interval border
         REAL(q)               :: FIT_COEFF_FUNCTION(4*N+4)
         ! polynomial coefficients of fit describing
         ! R-dependence of every quadrature coefficient and alternant
         REAL(q)               :: ALPHA( 24,100, 6 ) 

      END FUNCTION FIT_COEFF_FUNCTION

   END INTERFACE 

   CONTAINS

   FUNCTION ERROR_FUNCTION( QDT , Z ) 
      TYPE( quadrature_handle ) :: QDT
      QDPREAL :: Z
      QDPREAL :: ERROR_FUNCTION
      ! local 
      INTEGER :: I
      PROFILING_START( 'ERROR_FUNCTION' )

      ERROR_FUNCTION = QDT%F( Z )
      DO I = 1 , QDT%N 
         ERROR_FUNCTION = ERROR_FUNCTION - &
           QDT%C( I + QDT%N )*QDT%PHI2( Z, QDT%C(I) )
      ENDDO

      PROFILING_STOP( 'ERROR_FUNCTION' )
   END FUNCTION ERROR_FUNCTION

   FUNCTION D1_ERROR_FUNCTION( QDT , Z ) 
      TYPE( quadrature_handle ) :: QDT
      QDPREAL :: Z
      QDPREAL :: D1_ERROR_FUNCTION
      ! local 
      INTEGER :: I
      PROFILING_START( 'D1_ERROR_FUNCTION' )

      D1_ERROR_FUNCTION = QDT%D_F_DZ( Z )
      DO I = 1 , QDT%N 
         D1_ERROR_FUNCTION = D1_ERROR_FUNCTION - &
           QDT%C( I + QDT%N )*QDT%D_PHI2_DZ( Z, QDT%C(I) )
      ENDDO

      PROFILING_STOP( 'D1_ERROR_FUNCTION' )
   END FUNCTION D1_ERROR_FUNCTION

   FUNCTION D2_ERROR_FUNCTION( QDT , Z ) 
      TYPE( quadrature_handle ) :: QDT
      QDPREAL :: Z
      QDPREAL :: D2_ERROR_FUNCTION
      ! local 
      INTEGER :: I
      PROFILING_START( 'D2_ERROR_FUNCTION' )

      D2_ERROR_FUNCTION = QDT%D2_F_DZ2( Z )

      DO I = 1 , QDT%N 
         D2_ERROR_FUNCTION = D2_ERROR_FUNCTION - &
           QDT%C( I + QDT%N )*QDT%D2_PHI2_DZ2( Z, QDT%C(I) )
      ENDDO

      PROFILING_STOP( 'D2_ERROR_FUNCTION' )
   END FUNCTION D2_ERROR_FUNCTION

   FUNCTION GRAD_ERROR_FUNCTION( QDT , Z, J ) 
      TYPE( quadrature_handle ) :: QDT
      QDPREAL :: Z
      QDPREAL :: GRAD_ERROR_FUNCTION
      INTEGER :: J
      PROFILING_START( 'GRAD_ERROR_FUNCTION' )

      IF ( J <= QDT%N ) THEN
         GRAD_ERROR_FUNCTION = &
           -QDT%C( J + QDT%N )*QDT%D_PHI2_DL( Z, QDT%C(J) )
      ELSE
         GRAD_ERROR_FUNCTION = -QDT%PHI2( Z, QDT%C( J-QDT%N ) )
      ENDIF 

      PROFILING_STOP( 'GRAD_ERROR_FUNCTION' )
   END FUNCTION GRAD_ERROR_FUNCTION

!*******************************************************************************
!>helper,prints quadrature coefficients to file quad_IT.dat
!*******************************************************************************
   SUBROUTINE PRINT_QUADRATURE( QUADRATURE, IO, IT, FNAME )
      USE base
      USE prec
      TYPE( quadrature_handle )  :: QUADRATURE
      TYPE( in_struct )           :: IO
      ! local 
      INTEGER,INTENT(IN),OPTIONAL :: IT    ! add this number to file name
      CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: FNAME ! changes file name
      ! local 
      INTEGER                     :: I
      CHARACTER(LEN=4)            :: APP 

      PROFILING_START('PRINT_QUADRATURE')

      IF ( IO%IU0 >= 0 ) THEN
      IF (PRESENT(IT)) THEN
         WRITE(APP,'(I4)')IT
         APP=ADJUSTL(APP)
         IF ( PRESENT( FNAME ) ) THEN
            OPEN(UNIT=90,FILE=FNAME//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
         ELSE
            OPEN(UNIT=90,FILE='quad'//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
         ENDIF
      ELSE
         OPEN(UNIT=90,FILE='quad.dat',ACTION='WRITE',STATUS='REPLACE')
      ENDIF
      WRITE(90,'(2I4,3E13.6)')QUADRATURE%ITYPE,QUADRATURE%N,&
         TOPRINT(QUADRATURE%A), TOPRINT(QUADRATURE%B), TOPRINT(QUADRATURE%SCALING)
   
      !choose abscissas from coefficients
      DO I = 1, 2*QUADRATURE%N+1
         WRITE(90,'(2F45.32)')TOPRINT(QUADRATURE%C(I)),TOPRINT(QUADRATURE%X0(I))
      ENDDO

      CLOSE(90)
      ENDIF
      PROFILING_STOP('PRINT_QUADRATURE')
      RETURN
   END SUBROUTINE PRINT_QUADRATURE

!*******************************************************************************
!>helper,dumps quadrature coefficients to file 
!*******************************************************************************
   SUBROUTINE DUMP_QUADRATURE( QUADRATURE, IT, FNAME )
      USE prec
      USE base
      TYPE( quadrature_handle )  :: QUADRATURE
      ! local 
      INTEGER                    :: IT    
      CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: FNAME ! changes file name
      ! local 
      INTEGER                     :: I
      INTEGER                     :: IUNIT
      CHARACTER(LEN=4)            :: APP 

      PROFILING_START('DUMP_QUADRATURE')

      WRITE(APP,'(I4)')IT
      APP=ADJUSTL(APP)
      IF ( PRESENT( FNAME ) ) THEN
         OPEN(UNIT=IT,FILE=FNAME//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
      ELSE
         OPEN(UNIT=IT,FILE='quad'//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
      ENDIF
      WRITE(IT,'(2I4,3E13.6)')QUADRATURE%ITYPE,QUADRATURE%N,&
         TOPRINT(QUADRATURE%A), TOPRINT(QUADRATURE%B), TOPRINT(QUADRATURE%SCALING)
   
      !choose abscissas from coefficients
      DO I = 1, 2*QUADRATURE%N+1
         WRITE(IT,'(2F45.32)')TOPRINT(QUADRATURE%C(I)),TOPRINT(QUADRATURE%X0(I))
      ENDDO

      CLOSE(IT)
      PROFILING_STOP('DUMP_QUADRATURE')
      RETURN
   END SUBROUTINE DUMP_QUADRATURE

!*******************************************************************************
!>helper,dumps grids to file 
!*******************************************************************************
   SUBROUTINE DUMP_IGRID( IGRID, FNAME, INODE )
      USE prec
      TYPE( imag_grid_handle )  :: IGRID
      ! local 
      CHARACTER(LEN=*),INTENT(IN) :: FNAME ! changes file name
      INTEGER                     :: INODE
      ! local 
      INTEGER                     :: I,J
      INTEGER                     :: IT

      IF ( INODE == 1 ) THEN
         IT=100
         OPEN(UNIT=IT,FILE=FNAME//'.dat',ACTION='WRITE',STATUS='REPLACE')
         WRITE(IT,'(I4,L4,6E13.6)')IGRID%NOMEGA,IGRID%LFINITE_TEMPERATURE,&
            IGRID%BETA,IGRID%X1,IGRID%X2,IGRID%R,IGRID%FACTOR,IGRID%GM_INT_ERROR
      
         IF ( ASSOCIATED( IGRID%TAU ) ) THEN
         !dump tau grid
         WRITE(IT,'(A,E13.6)')'# TAU-GRID, Error=',IGRID%TAU_ERROR
         DO I = 1, IGRID%NOMEGA
            WRITE(IT,'(2F25.16)')IGRID%TAU(I),IGRID%TAU_WEIGHT(I)
         ENDDO
         ENDIF

         !dump BOS_RE grid
         IF ( ASSOCIATED( IGRID%BOS_RE ) ) THEN
         WRITE(IT,'(A,E13.6)')'# BOS_RE-GRID, Error=',IGRID%BOS_RE_ERROR
         DO I = 1, IGRID%NOMEGA
            WRITE(IT,'(2F25.16)')IGRID%BOS_RE(I),IGRID%BOS_RE_WEIGHT(I)
         ENDDO
         ENDIF
         IF ( ASSOCIATED( IGRID%TO_BOS_RE ) ) THEN
         !TO_BOS_RE transformation 
         IF ( IGRID%NOMEGA>6 ) THEN
            WRITE(IT,'(A)')'# TO_BOS_RE and TO_BOS_RE_CONJG transformation: '
            DO I = 1, IGRID%NOMEGA
            DO J = 1, IGRID%NOMEGA
               WRITE(IT,'(2I3,2F25.16)')I,J,IGRID%TO_BOS_RE(I,J),IGRID%TO_BOS_RE_CONJG(I,J)
            ENDDO
            ENDDO
         ELSE
            WRITE(IT,'(A)')'# TO_BOS_RE transformation: '
            DO I = 1, IGRID%NOMEGA
               WRITE(IT,'(6F25.16)')IGRID%TO_BOS_RE(I,:)
            ENDDO

            WRITE(IT,'(A)')'# TO_BOS_RE_CONJ transformation: '
            DO I = 1, IGRID%NOMEGA
               WRITE(IT,'(6F25.16)')IGRID%TO_BOS_RE_CONJG(I,:)
            ENDDO
         ENDIF
         WRITE(IT,'(A)')'# TO_BOS_RE_ERROR and TO_BOS_RE_CONJG_ERROR:'
         DO I = 1, IGRID%NOMEGA
            WRITE(IT,'(I3,2F25.16)')I,IGRID%TO_BOS_RE_ERROR(I),IGRID%TO_BOS_RE_CONJG_ERROR(I)
         ENDDO
         ENDIF

         IF ( ASSOCIATED( IGRID%TO_BOS_IM ) ) THEN
         !dump BOS_IM grid
         WRITE(IT,'(A,E13.6)')'# BOS_IM-GRID, Error=',IGRID%BOS_IM_ERROR
         DO I = 1, IGRID%NOMEGA
            WRITE(IT,'(2F25.16)')IGRID%BOS_IM(I),IGRID%BOS_IM_WEIGHT(I)
         ENDDO
         !TO_BOS_IM transformation 
         IF ( IGRID%NOMEGA>6 ) THEN
            WRITE(IT,'(A)')'# TO_BOS_IM and TO_BOS_IM_CONJG transformation: '
            DO I = 1, IGRID%NOMEGA
            DO J = 1, IGRID%NOMEGA
               WRITE(IT,'(2I3,2F25.16)')I,J,IGRID%TO_BOS_IM(I,J),IGRID%TO_BOS_IM_CONJG(I,J)
            ENDDO
            ENDDO
         ELSE
            WRITE(IT,'(A)')'# TO_BOS_IM transformation: '
            DO I = 1, IGRID%NOMEGA
               WRITE(IT,'(6F25.16)')IGRID%TO_BOS_IM(I,:)
            ENDDO

            WRITE(IT,'(A)')'# TO_BOS_IM_CONJ transformation: '
            DO I = 1, IGRID%NOMEGA
               WRITE(IT,'(6F25.16)')IGRID%TO_BOS_IM_CONJG(I,:)
            ENDDO
         ENDIF
         WRITE(IT,'(A)')'# TO_BOS_IM_ERROR and TO_BOS_IM_CONJG_ERROR:'
         DO I = 1, IGRID%NOMEGA
            WRITE(IT,'(I3,2F25.16)')I,IGRID%TO_BOS_IM_ERROR(I),IGRID%TO_BOS_IM_CONJG_ERROR(I)
         ENDDO
         ENDIF

         IF ( ASSOCIATED( IGRID%TO_FER_RE ) ) THEN
         !dump FER_RE grid
         WRITE(IT,'(A,E13.6)')'# FER_RE-GRID, Error=',IGRID%FER_RE_ERROR
         DO I = 1, IGRID%NOMEGA
            WRITE(IT,'(2F25.16)')IGRID%FER_RE(I),IGRID%FER_RE_WEIGHT(I)
         ENDDO
         !TO_FER_RE transformation 
         IF ( IGRID%NOMEGA>6 ) THEN
            WRITE(IT,'(A)')'# TO_FER_RE and TO_FER_RE_CONJG transformation: '
            DO I = 1, IGRID%NOMEGA
            DO J = 1, IGRID%NOMEGA
               WRITE(IT,'(2I3,2F25.16)')I,J,IGRID%TO_FER_RE(I,J),IGRID%TO_FER_RE_CONJG(I,J)
            ENDDO
            ENDDO
         ELSE
            WRITE(IT,'(A)')'# TO_FER_RE transformation: '
            DO I = 1, IGRID%NOMEGA
               WRITE(IT,'(6F25.16)')IGRID%TO_FER_RE(I,:)
            ENDDO

            WRITE(IT,'(A)')'# TO_FER_RE_CONJ transformation: '
            DO I = 1, IGRID%NOMEGA
               WRITE(IT,'(6F25.16)')IGRID%TO_FER_RE_CONJG(I,:)
            ENDDO
         ENDIF
         WRITE(IT,'(A)')'# TO_FER_RE_ERROR and TO_FER_RE_CONJG_ERROR:'
         DO I = 1, IGRID%NOMEGA
            WRITE(IT,'(I3,2F25.16)')I,IGRID%TO_FER_RE_ERROR(I),IGRID%TO_FER_RE_CONJG_ERROR(I)
         ENDDO
         ENDIF

         IF ( ASSOCIATED( IGRID%TO_FER_IM ) ) THEN
         !dump FER_IM grid
         WRITE(IT,'(A,E13.6)')'# FER_IM-GRID, Error=',IGRID%FER_IM_ERROR
         DO I = 1, IGRID%NOMEGA
            WRITE(IT,'(2F25.16)')IGRID%FER_IM(I),IGRID%FER_IM_WEIGHT(I)
         ENDDO
         !TO_FER_IM transformation 
         IF ( IGRID%NOMEGA>6 ) THEN
            WRITE(IT,'(A)')'# TO_FER_IM and TO_FER_IM_CONJG transformation: '
            DO I = 1, IGRID%NOMEGA
            DO J = 1, IGRID%NOMEGA
               WRITE(IT,'(2I3,2F25.16)')I,J,IGRID%TO_FER_IM(I,J),IGRID%TO_FER_IM_CONJG(I,J)
            ENDDO
            ENDDO
         ELSE
            WRITE(IT,'(A)')'# TO_FER_IM transformation: '
            DO I = 1, IGRID%NOMEGA
               WRITE(IT,'(6F25.16)')IGRID%TO_FER_IM(I,:)
            ENDDO

            WRITE(IT,'(A)')'# TO_FER_IM_CONJ transformation: '
            DO I = 1, IGRID%NOMEGA
               WRITE(IT,'(6F25.16)')IGRID%TO_FER_IM_CONJG(I,:)
            ENDDO
         ENDIF
         WRITE(IT,'(A)')'# TO_FER_IM_ERROR and TO_FER_IM_CONJG_ERROR:'
         DO I = 1, IGRID%NOMEGA
            WRITE(IT,'(I3,2F25.16)')I,IGRID%TO_FER_IM_ERROR(I),IGRID%TO_FER_IM_CONJG_ERROR(I)
         ENDDO
         ENDIF

      ENDIF

      CLOSE(IT)
   END SUBROUTINE DUMP_IGRID
   
!*******************************************************************************
!>helper,prints the nonlinear error function and its first two
! derivatives w.r.t. argument X
!*******************************************************************************
   SUBROUTINE PRINT_ERROR_FUNCTION( QUADRATURE, IO, IT, FNAME )
      USE base
      USE prec
      TYPE( quadrature_handle )  :: QUADRATURE
      TYPE( in_struct )           :: IO
      ! local 
      INTEGER,INTENT(IN),OPTIONAL :: IT    ! add this number to file name
      CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: FNAME ! changes file name
      ! local 
      INTEGER                     :: I,J,N
      QDPREAL                     :: X, X2, X1
      QDPREAL                     :: Y,D1Y,D2Y
      CHARACTER(LEN=4)            :: APP 
      INTEGER, PARAMETER          :: IMAX = 100

      PROFILING_START('PRINT_ERROR_FUNCTION')
   
      IF ( IO%IU0 >= 0 ) THEN
      IF (PRESENT(IT)) THEN
         WRITE(APP,'(I4)')IT
         APP=ADJUSTL(APP)
         IF ( PRESENT( FNAME ) ) THEN
            OPEN(UNIT=91,FILE=FNAME//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
         ELSE
            OPEN(UNIT=91,FILE='erf'//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
         ENDIF
      ELSE
         OPEN(UNIT=91,FILE='erf.dat',ACTION='WRITE',STATUS='REPLACE')
      ENDIF
      WRITE(91,*)'#',TOPRINT(QUADRATURE%A),TOPRINT(QUADRATURE%B)
      N = QUADRATURE%N
      DO J = 1, 2*N+1
         WRITE(91,'(A,2F50.32)')'# ',TOPRINT(QUADRATURE%C(J)),TOPRINT(QUADRATURE%X0(J))
      ENDDO
   
      !choose abscissas from coefficients
      X1 = QUADRATURE%A
      DO J = 1, N
         X2 = X1+QUADRATURE%C(J)*QUADRATURE%B/QUADRATURE%C(N)
         DO I=1,IMAX  !points in between two abscissas
            !equidistant abscissas between coefficients
            X=X1+((I-1)*(X2-X1))/IMAX
            IF( X > QUADRATURE%B) THEN
               X = QUADRATURE%B
               Y = ERROR_FUNCTION( QUADRATURE, X )
               D1Y = D1_ERROR_FUNCTION( QUADRATURE, X )
               D2Y = D2_ERROR_FUNCTION( QUADRATURE, X )
               WRITE(91,'(4E40.32)')TOPRINT(X),TOPRINT(Y),TOPRINT(D1Y),TOPRINT(D2Y)
               EXIT    
            ENDIF
            Y = ERROR_FUNCTION( QUADRATURE, X )
            D1Y = D1_ERROR_FUNCTION( QUADRATURE, X )
            D2Y = D2_ERROR_FUNCTION( QUADRATURE, X )
            WRITE(91,'(4E40.32)')TOPRINT(X),TOPRINT(Y),TOPRINT(D1Y),TOPRINT(D2Y)
         ENDDO
         X1 = X2
      ENDDO

      CLOSE(91)
      ENDIF 
      PROFILING_STOP('PRINT_ERROR_FUNCTION')
   END SUBROUTINE PRINT_ERROR_FUNCTION
   
!*******************************************************************************
!>helper,prints the nonlinear error function and its first two
! derivatives w.r.t. argument X
!*******************************************************************************
   SUBROUTINE PRINT_TRANS_ERROR_FUNCTION( QUAD_A, QUAD_B, ALPHA, OMEGA, &
      LCONJG, NODE_ME, IT )
      USE prec
      USE base
      TYPE( quadrature_handle )  :: QUAD_A
      TYPE( quadrature_handle )  :: QUAD_B
      QDPREAL                     :: ALPHA(:)
      QDPREAL                     :: OMEGA
      INTEGER                     :: NODE_ME 
      LOGICAL                     :: LCONJG
      INTEGER,INTENT(IN),OPTIONAL :: IT    ! add this number to file name
      ! local 
      INTEGER                     :: I, J, K
      QDPREAL                     :: X, X2, X1
      QDPREAL                     :: Y
      CHARACTER(LEN=4)            :: APP 
      INTEGER, PARAMETER          :: IMAX = 100

      PROFILING_START('PRINT_TRANS_ERROR_FUNCTION')
   
      IF (PRESENT(IT)) THEN
         WRITE(APP,'(I4)')IT
         APP=ADJUSTL(APP)
         IF ( LCONJG ) THEN
            OPEN(UNIT=100+NODE_ME,FILE='erf_w'//TRIM(APP)//'_conjg.dat',ACTION='WRITE',STATUS='REPLACE')
         ELSE
            OPEN(UNIT=100+NODE_ME,FILE='erf_w'//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
         ENDIF
      ELSE
         IF ( LCONJG ) THEN
            OPEN(UNIT=100+NODE_ME,FILE='erf_w_conjg.dat',ACTION='WRITE',STATUS='REPLACE')
         ELSE
            OPEN(UNIT=100+NODE_ME,FILE='erf_w.dat',ACTION='WRITE',STATUS='REPLACE')
         ENDIF
      ENDIF
      WRITE(100+NODE_ME,'(A,5E14.6)')'#',TOPRINT(QUAD_A%A),TOPRINT(QUAD_A%B),&
        TOPRINT(QUAD_B%A),TOPRINT(QUAD_B%B), TOPRINT( OMEGA )
   
      !choose abscissas from coefficients
      X1 = QUAD_B%A
      DO J = 1, QUAD_B%N
         IF ( QUAD_B%N > 1 .AND. ABS( QUAD_B%C(QUAD_B%N) ) > QDPP( 1.E-16 ) ) THEN
            X2 = X1+QUAD_B%C(J)*QUAD_B%B/QUAD_B%C(QUAD_B%N)
         ELSE
            X2 = X1+QUAD_B%B
         ENDIF
         DO I=1,IMAX  !points in between two abscissas
            !equidistant abscissas between coefficients
            X=X1+((I-1)*(X2-X1))/IMAX
            IF( X > QUAD_B%B) THEN
               X = QUAD_B%B
               IF ( LCONJG ) THEN
                  Y = QUAD_B%PHI_CONJG( X, OMEGA )
                  DO K = 1, QUAD_A%N
                     Y = Y - ALPHA( K )*QUAD_B%PSI_CONJG( X, QUAD_A%C( K ) )
                  ENDDO
               ELSE
                  Y = QUAD_B%PHI( X, OMEGA )
                  DO K = 1, QUAD_A%N
                     Y = Y - ALPHA( K )*QUAD_B%PSI( X, QUAD_A%C( K ) )
                  ENDDO
               ENDIF
               WRITE(100+NODE_ME,'(2E40.32)')TOPRINT(X),TOPRINT(Y)
               EXIT    
            ENDIF
            IF ( LCONJG ) THEN
               Y = QUAD_B%PHI_CONJG( X, OMEGA )
               DO K = 1, QUAD_A%N
                  Y = Y - ALPHA( K )*QUAD_B%PSI_CONJG( X, QUAD_A%C( K ) )
               ENDDO
            ELSE
               Y = QUAD_B%PHI( X, OMEGA )
               DO K = 1, QUAD_A%N
                  Y = Y - ALPHA( K )*QUAD_B%PSI( X, QUAD_A%C( K ) )
               ENDDO
            ENDIF
            WRITE(100+NODE_ME,'(2E40.32)')TOPRINT(X),TOPRINT(Y)
         ENDDO
         X1 = X2
      ENDDO
      CLOSE(100+NODE_ME)
      PROFILING_STOP('PRINT_TRANS_ERROR_FUNCTION')
   END SUBROUTINE PRINT_TRANS_ERROR_FUNCTION

!*******************************************************************************
!>helper,prints the nonlinear error function and its first two
! derivatives w.r.t. argument X
!*******************************************************************************
   SUBROUTINE PRINT_TRANSFORMATION_ERROR( OMEGA, ERRORS, N, LCONJG, NODE_ME, IT )
      USE prec
      USE base
      REAL(q)                     :: OMEGA(:)
      REAL(q)                     :: ERRORS(:)
      INTEGER                     :: N ! quadrature points come to separate file
      LOGICAL                     :: LCONJG
      INTEGER                     :: NODE_ME 
      INTEGER,INTENT(IN),OPTIONAL :: IT    ! add this number to file name
      ! local 
      INTEGER                     :: I
      CHARACTER(LEN=4)            :: APP 

      PROFILING_START('PRINT_TRANSFORMATION_ERROR')
      
      IF (PRESENT(IT)) THEN
         WRITE(APP,'(I4)')IT
         APP=ADJUSTL(APP)
         IF ( LCONJG ) THEN
            OPEN(UNIT=1000+NODE_ME,FILE='erf_trans'//TRIM(APP)//'_conjg.dat',ACTION='WRITE',STATUS='REPLACE')
            OPEN(UNIT=2000+NODE_ME,FILE='erf_trans_opt'//TRIM(APP)//'_conjg.dat',ACTION='WRITE',STATUS='REPLACE')
         ELSE
            OPEN(UNIT=1000+NODE_ME,FILE='erf_trans'//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
            OPEN(UNIT=2000+NODE_ME,FILE='erf_trans_opt'//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
         ENDIF
      ELSE
         IF ( LCONJG ) THEN
            OPEN(UNIT=1000+NODE_ME,FILE='erf_trans_conjg.dat',ACTION='WRITE',STATUS='REPLACE')
            OPEN(UNIT=2000+NODE_ME,FILE='erf_trans_opt_conjg.dat',ACTION='WRITE',STATUS='REPLACE')
         ELSE
            OPEN(UNIT=1000+NODE_ME,FILE='erf_trans.dat',ACTION='WRITE',STATUS='REPLACE')
            OPEN(UNIT=2000+NODE_ME,FILE='erf_trans_opt.dat',ACTION='WRITE',STATUS='REPLACE')
         ENDIF
      ENDIF
      !choose abscissas from coefficients
      DO I = 1, SIZE( ERRORS )-N
         WRITE(1000+NODE_ME,'( 2E20.10 )' )OMEGA(I),ERRORS(I)
      ENDDO
      CLOSE(1000+NODE_ME)

      !choose abscissas from coefficients
      DO I = SIZE(ERRORS)-N+1, SIZE( ERRORS )
         WRITE(2000+NODE_ME,'( 2E20.10 )' ) OMEGA(I), ERRORS(I)
      ENDDO
      CLOSE(2000+NODE_ME)

      PROFILING_STOP('PRINT_TRANSFORMATION_ERROR')

   END SUBROUTINE PRINT_TRANSFORMATION_ERROR

!*******************************************************************************
!>helper,prints Transformation Matrix to file
!*******************************************************************************
   SUBROUTINE PRINT_TRANSFORMATION_MATRIX( QUAD_A, QUAD_B, ALPHA, LCONJG, NODE_ME, IT )
      USE prec
      USE base
      TYPE( quadrature_handle )   :: QUAD_A
      TYPE( quadrature_handle )   :: QUAD_B
      REAL(q)                     :: ALPHA(:,:)
      INTEGER                     :: NODE_ME 
      LOGICAL                     :: LCONJG
      INTEGER,INTENT(IN),OPTIONAL :: IT    ! add this number to file name
      ! local 
      INTEGER                     :: I, J
      CHARACTER(LEN=4)            :: APP 
      INTEGER, PARAMETER          :: IMAX = 100

      PROFILING_START('PRINT_TRANSFORMATION_MATRIX')
   
      IF (PRESENT(IT)) THEN
         WRITE(APP,'(I4)')IT
         APP=ADJUSTL(APP)
         IF ( LCONJG ) THEN
            OPEN(UNIT=100+NODE_ME,FILE='trafo'//TRIM(APP)//'_conjg.dat',ACTION='WRITE',STATUS='REPLACE')
         ELSE
            OPEN(UNIT=100+NODE_ME,FILE='trafo'//TRIM(APP)//'.dat',ACTION='WRITE',STATUS='REPLACE')
         ENDIF
      ELSE
         IF ( LCONJG ) THEN
            OPEN(UNIT=100+NODE_ME,FILE='trafo_conjg.dat',ACTION='WRITE',STATUS='REPLACE')
         ELSE
            OPEN(UNIT=100+NODE_ME,FILE='trafo.dat',ACTION='WRITE',STATUS='REPLACE')
         ENDIF
      ENDIF
      WRITE(100+NODE_ME,'(A,4E14.6,L4)')'#',TOPRINT(QUAD_A%A),TOPRINT(QUAD_A%B),&
        TOPRINT(QUAD_B%A),TOPRINT(QUAD_B%B), LCONJG
      WRITE(100+NODE_ME,'(A)')'#     I      J              ALPHA(I,J)       '
   
      !choose abscissas from coefficients
      DO I = 1, QUAD_A%N
      DO J = 1, QUAD_B%N
         WRITE(100+NODE_ME,'(2I7,F25.16)')I,J,ALPHA(I,J)
      ENDDO
      ENDDO
      CLOSE(100+NODE_ME)
      PROFILING_STOP('PRINT_TRANSFORMATION_MATRIX')
   END SUBROUTINE PRINT_TRANSFORMATION_MATRIX

#ifdef MinimaxTest
!*******************************************************************************
!>helper, performs a simple integration test of quadrature
!>by calculating a model using the energies found in W
!*******************************************************************************

   SUBROUTINE MINIMAX_INTEGRATION_TEST( IMAG_GRIDS, W, EFERMI_, IO )
      USE prec
      USE base
      USE wave_struct_def
      USE minimax_functions2D
      TYPE(imag_grid_handle) :: IMAG_GRIDS 
      TYPE (wavespin)        :: W          ! wavefunction
      REAL(q)                :: EFERMI_
      TYPE (in_struct)       :: IO
      ! local 
      INTEGER                :: NK, NB, ISP
      INTEGER                :: I,J
      REAL(q)                :: GO_TAU(IMAG_GRIDS%NOMEGA)  ! Green's function in time
      REAL(q)                :: GU_TAU(IMAG_GRIDS%NOMEGA)  ! Green's function in time
      REAL(q)                :: CHIU_TAU(IMAG_GRIDS%NOMEGA)  ! Polarizability 
      REAL(q)                :: CHIO_TAU(IMAG_GRIDS%NOMEGA)  ! Polarizability 
      REAL(q)                :: SIGMAU_TAU(IMAG_GRIDS%NOMEGA)  ! Polarizability 
      REAL(q)                :: SIGMAO_TAU(IMAG_GRIDS%NOMEGA)  ! Polarizability 
      COMPLEX(q)             :: G_OMEGA(IMAG_GRIDS%NOMEGA)  ! Green's function in frequency
      COMPLEX(q)             :: G2_OMEGA(IMAG_GRIDS%NOMEGA)  ! Green's function in frequency
      COMPLEX(q)             :: CHI_OMEGA(IMAG_GRIDS%NOMEGA)  ! Polarizability 
      COMPLEX(q)             :: SIGMA_OMEGA(IMAG_GRIDS%NOMEGA)  ! Self-energy in time
      COMPLEX(q)             :: SIGMA2_OMEGA(IMAG_GRIDS%NOMEGA)  ! Self-energy in time
      REAL(q)                :: GM         ! Galitskii-Migdal energy
      REAL(q)                :: E, F, BETA, TAU, OMEGA
      INTEGER                :: NK2, NB2, ISP2, NSTATES
      REAL(q)                :: E2, F2
      REAL(q)                :: RE, IM
      REAL(q), POINTER       :: DELTA(:,:)
      INTEGER                :: ISPIN, NKPTS, NB_TOT
      REAL(q)                :: EFERMI, EMIN , EMAX 
      COMPLEX(q)             :: CELTOT(W%WDES%NB_TOT, W%WDES%NKPTS, W%WDES%ISPIN)
      REAL(q)                :: FERTOT(W%WDES%NB_TOT, W%WDES%NKPTS, W%WDES%ISPIN)
      ! initialize 
      BETA = IMAG_GRIDS%BETA
      EFERMI=EFERMI_
      EMIN=1E10
      EMAX=-1E10

      ISPIN = W%WDES%ISPIN
      NKPTS = W%WDES%NKPTS
      NB_TOT = W%WDES%NB_TOT 

      CELTOT=W%CELTOT
      FERTOT=W%FERTOT

!      ISPIN = 1
!      NKPTS=1
!      NB_TOT=W%WDES%NB_TOT
!      EFERMI=0
!      CALL MODEL( ISPIN, NKPTS, NB_TOT, CELTOT, FERTOT, BETA )
!
!
! Set up Green's function, in frequency and time 
! 
      GO_TAU = 0 
      GU_TAU = 0 
      G_OMEGA = 0 
      CHIU_TAU = 0 
      CHIO_TAU = 0 
      CHI_OMEGA = 0 
      NSTATES = ISPIN*NKPTS*NB_TOT
      IF ( IO%IU0>=0.AND.IO%NWRITE>3 ) &
         WRITE(*,*) "       TAU                GO(-TAU)            GU( TAU )     CHIU_TAU     CHIO_TAU"
      ! set up Correlation functions in time or frequency domain
      DO I = 1, IMAG_GRIDS%NOMEGA 
         TAU=IMAG_GRIDS%TAU( I ) 
         OMEGA=IMAG_GRIDS%FER_RE( I ) 
         DO ISP = 1, ISPIN
           DO NK = 1, NKPTS
              DO NB = 1, NB_TOT
                 ! set energy 
                 E = REAL( CELTOT(NB,NK,ISP) - EFERMI,q )
                 EMIN = MIN( E, EMIN ) 
                 EMAX = MAX( E, EMAX ) 
                 ! set fermi occupancy 
                 F = REAL( FERTOT(NB,NK,ISP), q)
                 ! green's function, work in -beta/2 and beta/2 domain
                 ! negative times
                 GO_TAU( I ) = GO_TAU( I ) + GOCC_WEIGHT_TIME( E, F, BETA,-TAU)/NSTATES
                 GU_TAU( I ) = GU_TAU( I ) + GUNO_WEIGHT_TIME( E, F, BETA, TAU)/NSTATES
                 ! in frequency domain one has
                 G_OMEGA( I ) = G_OMEGA( I ) + 1._q/(E-CMPLX(0._q,OMEGA,q))/NSTATES

                  DO ISP2 = 1, ISPIN
                    DO NK2 = 1, NKPTS
                       DO NB2 = 1, NB_TOT
                          ! set energy 
                          E2 = REAL( CELTOT(NB2,NK2,ISP2) - EFERMI,q )
                          ! set fermi occupancy 
                          F2 = REAL( FERTOT(NB2,NK2,ISP2), q)
                           
                          ! positive time chi
                          CHIU_TAU( I ) = CHIU_TAU( I ) &
                           + GUNO_WEIGHT_TIME( E, F, BETA, TAU)    &
                           * GOCC_WEIGHT_TIME( E2, F2, BETA,-TAU)/NSTATES/NSTATES
                          ! negative time chi
                          !CHIO_TAU( I ) = CHIO_TAU( I ) &
                          ! + GUNO_WEIGHT_TIME( E, F, BETA, BETA-TAU)    &
                          ! * GOCC_WEIGHT_TIME( E2, F2, BETA,TAU-BETA)
 
                          CHIO_TAU( I ) = CHIO_TAU( I ) &
                           + GUNO_WEIGHT_TIME( E, F, BETA, BETA-TAU)    &
                           * GOCC_WEIGHT_TIME( E2, F2, BETA,TAU-BETA)/NSTATES/NSTATES
 
                          ! positive frequencies
                          CHI_OMEGA( I ) = CHI_OMEGA( I ) &
                            - GG_OMEGA_WEIGHT( E, F, E2, F2, BETA, &
                            IMAG_GRIDS%BOS_RE(I), IMAG_GRIDS%BOS_RE(1) )/NSTATES/NSTATES
 
                       ENDDO
                    ENDDO
                  ENDDO


              ENDDO
           ENDDO
         ENDDO
         IF( IO%IU0>=0 .AND. IO%NWRITE>2 )WRITE(*,5)TAU, GO_TAU( I ) , GU_TAU( I), CHIO_TAU( I ), CHIU_TAU(I)
      ENDDO
IF ( IO%IU0>=0 ) WRITE(*,'(A,3F20.10)' )" EMIN/EMAX",BETA*EMIN, BETA*EMAX , IMAG_GRIDS%R

!
! set up polarizability in time and frequency 
! this also sets up the joined transition energies DELTA
! 
       CALL SETUP_SIGMA( W, DELTA, SIGMA_OMEGA, SIGMAO_TAU, SIGMAU_TAU, IMAG_GRIDS )
       IF ( IO%IU0>=0.AND.IO%NWRITE>3 ) &
         WRITE(*,*) "      OMEGA           G_OMEGA(OMEGA)*SIGMA_OMEGA(OMEGA)       GM "
       CALL INTEGRATE_GM( G_OMEGA, SIGMA_OMEGA, IMAG_GRIDS, GM ) 
       IF( IO%IU0>=0 )WRITE(*,'(" GM = ",I10, 3F25.12)')IMAG_GRIDS%NOMEGA, GM
       IF( IO%IU6>=0 )WRITE(IO%IU6,'(" GMTEST = ",I10, 3F25.12)')IMAG_GRIDS%NOMEGA, GM
!
! determine transformation error
! 
      E=0
      F=0
      E2=0
      F2=0
      G2_OMEGA = G_OMEGA
      G_OMEGA = 0 
      SIGMA2_OMEGA = SIGMA_OMEGA
      SIGMA_OMEGA = 0 
      ! perform transform to frequency domain
      IF ( IO%IU0>=0.AND.IO%NWRITE>3 ) &
         WRITE(*,*) "      OMEGA           RE:G_OMEGA(OMEGA)   IM:G_OMEGA(OMEGA)       FT "
      DO I = 1, IMAG_GRIDS%NOMEGA 
         OMEGA=IMAG_GRIDS%FER_RE( I ) 
         DO J = 1, IMAG_GRIDS%NOMEGA
            G_OMEGA( I ) = G_OMEGA( I ) + &
               ( GU_TAU(J) - GO_TAU(J) )*IMAG_GRIDS%TO_FER_RE( I, J ) 
            ! imaginary part 
            G_OMEGA( I ) = G_OMEGA( I )  + &
               (GU_TAU(J) + GO_TAU(J) )*IMAG_GRIDS%TO_FER_RE_CONJG( I, J)*CMPLX(0._q, 1._q,q)

            ! real part 
            SIGMA_OMEGA( I ) = SIGMA_OMEGA( I ) + &
               ( SIGMAU_TAU(J) - SIGMAO_TAU(J) )*IMAG_GRIDS%TO_FER_RE( I, J ) 
            ! imaginary part 
            SIGMA_OMEGA( I ) = SIGMA_OMEGA( I )+ &
               (SIGMAU_TAU(J) + SIGMAO_TAU(J) )*IMAG_GRIDS%TO_FER_RE_CONJG( I, J)*CMPLX(0._q, 1._q,q)

            ! polarizabilty 
!            SIGMA2_OMEGA( I ) = SIGMA2_OMEGA( I ) + (CHIU_TAU(J)+CHIO_TAU(J))/2*IMAG_GRIDS%TO_BOS_RE( I, J )
         ENDDO

         E = E + ABS( REAL( G2_OMEGA(I)-G_OMEGA(I) ,q))
         F = F + ABS( AIMAG( G2_OMEGA(I)-G_OMEGA(I) ) )

         E2 = E2 + ABS( REAL( SIGMA_OMEGA(I)-SIGMA2_OMEGA(I) ,q))
         F2 = F2 + ABS( AIMAG( SIGMA_OMEGA(I)-SIGMA2_OMEGA(I) ) )
         IF( IO%IU0>=0 .AND. IO%NWRITE>2 )WRITE(*,7)I,OMEGA, SIGMA_OMEGA(I), SIGMA2_OMEGA(I)
         !IF( IO%IU0>=0 .AND. IO%NWRITE>2 )WRITE(*,7)I,OMEGA, G_OMEGA(I), G2_OMEGA(I)
      ENDDO
      IF( IO%IU0>=0 )WRITE(*,6)" Total FT-Error (Re, Im)= ", E, F, E2,F2

      ! integrate
      CALL INTEGRATE_GM( G_OMEGA, SIGMA_OMEGA, IMAG_GRIDS, GM ) 
      IF( IO%IU0>=0 )WRITE(*,'(" GM FT= ",I10, 3F25.12)')IMAG_GRIDS%NOMEGA, GM
      IF( IO%IU6>=0 )WRITE(IO%IU6,'(" GMTESTFT = ",I10, 3F25.12)')IMAG_GRIDS%NOMEGA, GM


5     FORMAT( 8F20.10 )
6     FORMAT( A, 4E17.10 )
7     FORMAT( I6, 5F20.10 )

      CONTAINS
      ! occupied greens function weight for negative TAU 
      FUNCTION GOCC_WEIGHT_TIME( E, F, BETA, TAU )  
        REAL(q)    :: E,F,BETA, TAU
        REAL(q)    :: GOCC_WEIGHT_TIME
       
        ! if state if partially unoccupied 
        IF ( E > 0 ) THEN
           GOCC_WEIGHT_TIME=(1-F)*EXP(-E*(BETA+TAU))
        ! fully unoccupied
        ELSE
           GOCC_WEIGHT_TIME=F*EXP(-E*TAU)
        ENDIF
      END FUNCTION GOCC_WEIGHT_TIME

      ! unoccupied greens function for positive TAU
      FUNCTION GUNO_WEIGHT_TIME( E, F, BETA, TAU )  
        REAL(q)    :: E,F,BETA, TAU
        REAL(q)    :: GUNO_WEIGHT_TIME
       
        ! if state if partially occupied 
        IF ( E < 0 ) THEN
           GUNO_WEIGHT_TIME=F*EXP(E*(BETA-TAU))
        ! unoccupied or at fermi level
        ELSE
           GUNO_WEIGHT_TIME=(1-F)*EXP(-E*TAU)
        ENDIF
      END FUNCTION GUNO_WEIGHT_TIME

      ! set up chi in omega directly
      FUNCTION GG_OMEGA_WEIGHT(EU,FU, EO,FO, BETA, OMEGA, OMEGA0 )
        REAL(q) :: GG_OMEGA_WEIGHT
        REAL(q) :: WEIGHT
        REAL(q)    :: EU,FU, EO,FO, BETA, TAU, OMEGA0, D
        REAL(q) :: OMEGA
        ! local 
        ! integrates finite temperature time dependence in the interval [0,1] 
        ! with precision of 1E-16 factor 2 is missing for weights 
        REAL(q)    :: TIMEQUAD(6)
        DATA  TIMEQUAD(1:6)/&
        0.033620790194994629174196631993254414_q,& 
        0.168952138582706495872898244670068379_q,&
        0.380421654180154000801650227003847249_q,&
        0.042667387775989003284582423702886445_q,&
        0.090122670444184774241058732968667755_q,&
        0.117209941779825937979708783132082317_q/
        INTEGER    :: I
        INTEGER,PARAMETER :: N=3
    
        ! nummerically save way to evaluate the frequency dependence
        ! for finite temperature:
        !
        ! this is the time dependence
        ! X_uo(t) = f_o(1-f_u) e^{-D_uo t}, with D_uo = e_u-e_o 
        ! Fourier coefficient is obtained via integration
        ! int_0^beta dt X_uo(t) Cos( w t ) = 
        ! f_o(1-f_u) * [ D_uo/(D_uo^2 + w^2 )( 1-e^{-beta(e_u-e_o) } ) 
        !                + beta delta(D_uo) delta( w) ]
        ! 
        ! the first term becomes (including fermi weights): 
        ! D_uo/(D_uo^2 + w^2)*1/2( tanh( beta e_u/2 ) - tanh( beta e_o/2 ) )
        ! = D_uo/(D_uo^2 + w^2 )(f_o - f_u)
        ! 
        ! the first and second term are zero if w /= 0 and D_uo = 0 ,
        ! but non-zero if w = 0 and D_uo =0
         
        D = EU-EO
        IF ( D < 0 ) THEN
          GG_OMEGA_WEIGHT=0
          RETURN
        ENDIF
        ! metallic contributions 
        IF ( ABS( BETA*D ) < 1.E-4_q ) THEN
           ! for zero transition, there is a contribution only for w->0 from both terms
           ! which for finite temperature is the first frequency point OMEGA0
           IF ( ABS(OMEGA)-OMEGA0 < 1.E-6_q .OR. OMEGA==0) THEN
              ! beta* f_o*(1-f_u) comes from the first term  
              ! the same contribution from the second term 
              ! up to second order
              !WEIGHT = BETA*FO*(1-FO)/2*(2-EO*(2*FO-1)*BETA-EU*BETA*(1-2*FO))
              WEIGHT = 0 
              ! determine time weight and integrate to obtain omega = 0 value
              DO I = 1,  N 
                 WEIGHT = WEIGHT + ( 4 * BETA * TIMEQUAD(N+I) ) * &
                          GUNO_WEIGHT_TIME( EU, FU, BETA, TIMEQUAD(I)*BETA )*&
                          GOCC_WEIGHT_TIME( EO, FO, BETA,-TIMEQUAD(I)*BETA )
              ENDDO
           ENDIF
        ELSE 
           ! for non-zero transitions, one has following term 
           ! D_mn/(D_mn^2+w^2)(f_n-f_m) 
           WEIGHT = (D/(D**2+OMEGA**2))*(FO-FU) 
        ENDIF
        ! factor 2 missing for non self-exitations 
        IF ( D/=0 ) WEIGHT = WEIGHT * 2 
         
        GG_OMEGA_WEIGHT = -WEIGHT
      END FUNCTION GG_OMEGA_WEIGHT

      ! determine joined density of states 
      SUBROUTINE SETUP_CHI( W, DELTA, CHI_OMEGA, IMAG_GRIDS )
         TYPE (wavespin)         :: W          ! wavefunction
         REAL(q), POINTER        :: DELTA(:,:)
         TYPE(imag_grid_handle)  :: IMAG_GRIDS 
         COMPLEX(q),INTENT(INOUT):: CHI_OMEGA(:)  ! Green's function in frequency
         !local 
         INTEGER                :: ISP1, NK1, NB1
         INTEGER                :: N
         INTEGER                :: I, J, K
         INTEGER                :: NSTATES
         INTEGER                :: NSTATES2
         REAL(q),ALLOCATABLE    :: EPS(:), FER(:)
         
         ! number of states
         NSTATES = (ISPIN*W%WDES%NKPTS*NB_TOT)
         ALLOCATE( EPS( NSTATES ), FER( NSTATES ) ) 
         ! number of joined states
         NSTATES = NSTATES*(NSTATES+1)/2
         NSTATES2 = NSTATES*(NSTATES+1)/2
         ALLOCATE( DELTA( NSTATES, 2 ) ) 
         DELTA = 0 

         ! store states and occupancy in one array 
         N=0
         DO ISP1 = 1, ISPIN
           DO NK1 = 1, NKPTS
              DO NB1 = 1, NB_TOT
                 N=N+1
                 EPS( N ) = REAL( CELTOT( NB1, NK1, ISP1 ), q)
                 FER( N ) = REAL( FERTOT( NB1, NK1, ISP1 ), q)
              ENDDO
           ENDDO
         ENDDO
         ! this is the number of joined states
         NSTATES = N

         CHI_OMEGA = 0 
         ! now determine joined density 
         N=0
         DO J = 1, NSTATES
            DO I = 1, NSTATES
               ! upper triangle only
               IF ( I<J ) CYCLE 
               N=N+1
               DELTA( N, 1 ) = EPS( I )-EPS( J ) 
               DELTA( N, 2 ) = FER( I )*FER( J ) 
               DO K = 1, IMAG_GRIDS%NOMEGA 
                  IF ( DELTA( N, 1 ) == 0 ) THEN
                     IF ( IMAG_GRIDS%BOS_RE(K)==0 ) THEN

                     ENDIF

                  ELSE
                     CHI_OMEGA(K)=CHI_OMEGA(K) + DELTA(N,1)/(DELTA(N,1)**2+IMAG_GRIDS%BOS_RE(K)**2)/NSTATES*(FER(I)-FER(J))  
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         DEALLOCATE( EPS, FER )

      END SUBROUTINE SETUP_CHI

      ! determine joined density of states 
      SUBROUTINE SETUP_SIGMA( W, DELTA, SIGMA_OMEGA, SIGMAO_TAU, SIGMAU_TAU, IMAG_GRIDS )
         TYPE (wavespin)         :: W          ! wavefunction
         REAL(q), POINTER        :: DELTA(:,:)
         TYPE(imag_grid_handle)  :: IMAG_GRIDS 
         COMPLEX(q),INTENT(INOUT):: SIGMA_OMEGA(:)  ! Green's function in frequency
         REAL(q),INTENT(INOUT):: SIGMAO_TAU(:)  ! Green's function in frequency
         REAL(q),INTENT(INOUT):: SIGMAU_TAU(:)  ! Green's function in frequency
         !local 
         INTEGER                :: ISP1, NK1, NB1
         INTEGER                :: N
         INTEGER                :: I, J, K
         INTEGER                :: NSTATES
         INTEGER                :: NSTATES2
         REAL(q)                :: EMIN, EMAX 
         REAL(q),ALLOCATABLE    :: EPS(:), FER(:)
         EMIN = 1E10
         EMAX = -1E10
         
         ! number of states
         NSTATES = (ISPIN*NKPTS*NB_TOT)
         ALLOCATE( EPS( NSTATES ), FER( NSTATES ) ) 
         ! number of joined states
         NSTATES = NSTATES*(NSTATES+1)/2
         NSTATES2 = NSTATES*(NSTATES+1)/2
         ALLOCATE( DELTA( NSTATES, 2 ) ) 
         DELTA = 0 

         ! store states and occupancy in one array 
         N=0
         DO ISP1 = 1, ISPIN
           DO NK1 = 1, NKPTS
              DO NB1 = 1, NB_TOT
                 N=N+1
                 EPS( N ) = REAL( CELTOT( NB1, NK1, ISP1 ), q)
                 FER( N ) = REAL( FERTOT( NB1, NK1, ISP1 ), q)
              ENDDO
           ENDDO
         ENDDO
         ! this is the number of joined states
         NSTATES = N

         SIGMA_OMEGA = 0 
         SIGMAO_TAU = 0 
         SIGMAU_TAU = 0 
         ! now determine joined density 
         DO K = 1, IMAG_GRIDS%NOMEGA 
            N=0
            DO J = 1, NSTATES
               DO I = 1, NSTATES
               ! upper triangle only
               IF ( I<J ) CYCLE 
               N=N+1
               DELTA( N, 1 ) = EPS( I )-EPS( J ) 
               DELTA( N, 2 ) = FERMI( DELTA(N,1)*BETA )
               EMIN = MIN( EMIN, DELTA( N , 1 ) ) 
               EMAX = MAX( EMAX, DELTA( N , 1 ) ) 

               SIGMA_OMEGA(K)=SIGMA_OMEGA(K) + 1._q/(DELTA(N,1)-CMPLX(0._q,IMAG_GRIDS%FER_RE(K),q))/NSTATES  
               SIGMAO_TAU( K ) = SIGMAO_TAU( K ) + GOCC_WEIGHT_TIME( DELTA(N,1), DELTA(N,2), BETA,-IMAG_GRIDS%TAU(K))/NSTATES
               SIGMAU_TAU( K ) = SIGMAU_TAU( K ) + GUNO_WEIGHT_TIME( DELTA(N,1), DELTA(N,2), BETA, IMAG_GRIDS%TAU(K))/NSTATES
!IF( IO%IU0>=0 )WRITE(100,'("SIGMA",I10, 4F25.12)')N,DELTA(N,:)
              ENDDO
            ENDDO
            IF( IO%IU0>=0 .AND. IO%NWRITE>2)WRITE(*,'("SIGMA",I10, 4F25.12)')K, IMAG_GRIDS%FER_RE(K),SIGMAO_TAU(K),SIGMAU_TAU(K)
         ENDDO
         DEALLOCATE( EPS, FER )
         IF ( IO%IU0>=0 ) WRITE(*,'(A,3F20.10)' )" DELTAMIN/MAX",BETA*EMIN, BETA*EMAX , IMAG_GRIDS%R
      END SUBROUTINE SETUP_SIGMA

      ! define FERMI FUNCTION
      FUNCTION FERMI( X ) 
         REAL(q) :: X,Y 
         REAL(q) :: FERMI
         REAL(q), PARAMETER :: MAXE=100._q

         ! evaluate 
         IF ( X < 0 ) THEN
            IF ( ABS( X ) < MAXE ) THEN
                Y = EXP( X ) + 1._q 
            ELSE
                Y = EXP( -MAXE ) + 1._q 
            ENDIF 
         ELSE
            IF ( ABS( X ) < MAXE ) THEN
                Y = EXP( X ) + 1._q 
            ELSE
                Y = EXP( MAXE ) + 1._q 
            ENDIF 
         ENDIF
         FERMI = 1._q/Y
      END FUNCTION FERMI
     
      SUBROUTINE MODEL( ISPIN, NKPTS, NB_TOT, CELTOT, FERTOT, BETA )
         INTEGER   ,INTENT(INOUT) :: ISPIN, NKPTS, NB_TOT
         COMPLEX(q),INTENT(INOUT) :: CELTOT(NB_TOT, NKPTS, ISPIN)
         REAL(q)   ,INTENT(INOUT) :: FERTOT(NB_TOT, NKPTS, ISPIN)
         REAL(q)   ,INTENT(IN)    :: BETA 
         ! LOCAL 
         INTEGER                  :: ISP, NK, NB, N, K 
         REAL(q), PARAMETER       :: E_MIN = -10._q
         REAL(q), PARAMETER       :: E_MAX =  10._q
       
         NB_TOT = MAX( 2, NB_TOT ) 
         N = ISPIN*NKPTS*NB_TOT

         K = 0 
         DO ISP = 1, ISPIN
            DO NK = 1, NKPTS
               DO NB = 1, NB_TOT
                  K = K+1
                  CELTOT( NB, NK, ISP ) = E_MIN + (K-1)*(E_MAX-E_MIN)/(N-1)
                  FERTOT( NB, NK, ISP ) = FERMI( BETA*REAL(CELTOT( NB, NK, ISP),q) )  
               ENDDO
            ENDDO
         ENDDO

         CELTOT( NB_TOT, NKPTS, ISPIN ) = 0
         FERTOT( NB_TOT, NKPTS, ISPIN ) = FERMI( BETA*REAL(CELTOT( NB_TOT, NKPTS, ISPIN),q) )  

         CELTOT( NB_TOT-1, NKPTS, ISPIN ) = 0.5_q
         FERTOT( NB_TOT-1, NKPTS, ISPIN ) = FERMI( BETA*REAL(CELTOT( NB_TOT-1, NKPTS, ISPIN),q) )  
         CELTOT( NB_TOT-2, NKPTS, ISPIN ) = 0.5_q
         FERTOT( NB_TOT-2, NKPTS, ISPIN ) = FERMI( BETA*REAL(CELTOT( NB_TOT-2, NKPTS, ISPIN),q) )  
         CELTOT( NB_TOT-3, NKPTS, ISPIN ) = -0.1_q
         FERTOT( NB_TOT-3, NKPTS, ISPIN ) = FERMI( BETA*REAL(CELTOT( NB_TOT-3, NKPTS, ISPIN),q) )  
! 
!          K = 0 
!          DO ISP = 1, ISPIN
!             DO NK = 1, NKPTS
!                DO NB = 1, NB_TOT
!                   K = K+1
!                    WRITE(*,*)K, BETA*REAL(CELTOT( NB, NK, ISP ),q), FERTOT( NB, NK, ISP )
! 
!                ENDDO
!             ENDDO
!          ENDDO
! 
      END SUBROUTINE MODEL 

      SUBROUTINE  INTEGRATE_GM( G_OMEGA, SIGMA_OMEGA, IMAG_GRIDS, GM ) 
         INTEGER :: I 
         COMPLEX(q),INTENT(IN):: G_OMEGA(:)  ! Green's function in frequency
         COMPLEX(q),INTENT(IN):: SIGMA_OMEGA(:)  ! Green's function in frequency
         TYPE(imag_grid_handle)  :: IMAG_GRIDS 
         REAL(q)                :: GM
         GM=0
         DO I = 1, IMAG_GRIDS%NOMEGA
            OMEGA=IMAG_GRIDS%FER_RE(I)
            ! compute even part 
            ! (G+G*) x (S+S*)
            ! this part can be integrated with FER_RE grid
            GM=GM+ ( G_OMEGA(I) + CONJG( G_OMEGA(I) ) ) &
                 * ( SIGMA_OMEGA( I ) + CONJG( SIGMA_OMEGA(I) ) )&
                 * IMAG_GRIDS%FER_RE_WEIGHT( I )/4
            ! compute odd part 
            ! (G-G*) x (S-S*) 
            ! this part can be integrated with FER_IM grid
            GM=GM+ ( G_OMEGA(I) - CONJG( G_OMEGA(I) ) ) &
                 * ( SIGMA_OMEGA( I ) - CONJG( SIGMA_OMEGA(I) ) )&
                 * IMAG_GRIDS%FER_IM_WEIGHT( I )/4

            IF( IO%IU0>=0 .AND. IO%NWRITE>5 )WRITE(*,50)OMEGA, G_OMEGA(I)*SIGMA_OMEGA(I), GM
         ENDDO

50     FORMAT( 8F20.10 )
      END SUBROUTINE 

   END SUBROUTINE MINIMAX_INTEGRATION_TEST
#endif 
END MODULE minimax_struct
