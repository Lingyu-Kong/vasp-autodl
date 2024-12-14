#include "symbol.inc"

MODULE minimax
   USE prec
   USE base, ONLY : in_struct, TOREAL
#ifdef qd_emulate
   USE qdmodule
#endif
   USE minimax_struct
   USE minimax_varpro
   USE mpimy
   USE minimax_functions1D
   USE minimax_functions2D
   
   IMPLICIT NONE 
   PRIVATE
   PUBLIC :: SET_IMAG_GRID_HANDLE, &
             ALLOCATE_IMAG_GRID, &
             DEALLOCATE_IMAG_GRID_HANDLE, &
             LOCAL_INDEX_CYCLIC, &
             LOCAL_INDEX_TO_GLOBAL, &
             ADD_POINT_IMAG_GRID_HANDLE

#ifdef qd_emulate
#define ACC      qdreal('1.E-18')
#define ACC2     qdreal('1.E-12')
#define VERYTINY qdreal('1.E-28')
#define PIQ      qdpi()
#else
   !accuracy in minimization routines
   QDPREAL , PRIVATE            :: ACC=1.E-18_qd  
   !truncation threshold for SVD and accuracy in FT routine
   QDPREAL , PRIVATE            :: ACC2=1.E-12_qd   
   !tiny number for QDLUDCMP 
   QDPREAL , PRIVATE            :: VERYTINY=1.E-28_qd
#endif

   INTEGER, PRIVATE             :: ITMAXLSQ = 1500
   INTEGER, PRIVATE             :: NDATALSQ = 100

CONTAINS

!****************************************************************************
!>Main calling routine, caluclates imaginary time and frequency grid
!!and corresponding forward Fouier transformation matrix,i.e.
!!transformation from imaginary time to imanginary frequency  
!!optionally the sin transform is calculated 
!****************************************************************************

   SUBROUTINE SET_IMAG_GRID_HANDLE(X1, X2, NTAU, ITYP, SIGMA, IMAG_GRIDS, COMM, IO)
      USE prec 
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      REAL(q)                  :: X1, X2     ! energy interval 
      INTEGER                  :: NTAU       ! number of time/freq points
      INTEGER                  :: ITYP       ! type of grid requested
      REAL(q)                  :: SIGMA      ! temperature in K/eV
      TYPE(imag_grid_handle)   :: IMAG_GRIDS ! grid handle
      TYPE(communic)           :: COMM       ! communicator
      TYPE(in_struct)          :: IO         ! for output and input
      !local
      INTEGER                  :: NOMEGA     ! loacl number of frequency points
      QDPREAL                  :: RSTART     ! starting interval 
      TYPE( quadrature_handle) :: TAU_TMP    ! temporary quadrature in time

      PROFILING_START( 'SET_IMAG_GRID_HANDLE' )
      ! currently not more than 24 grid points supported   
      IF ( NTAU > 24 .AND. ITYP < 168 ) THEN 
         CALL vtutor%error("Currently not more than 24 grid points supported " // str(NTAU))
      ENDIF
      ! store requested minimization interval 
      IMAG_GRIDS%X1 = X1
      IMAG_GRIDS%X2 = X2
       
      IF (IO%IU6>=0) WRITE(IO%IU6,*)'set NWRITE>2 for more details about minimax calculation'
      IF (IO%IU0>=0) WRITE(IO%IU0,*)'set NWRITE>2 for more details about minimax calculation'
      ! initialize quadratures, this sets correct pointers to 
      ! basis functions as well as to the object functions of all quadratures
      CALL INITIALIZE_QUADRATURES( IMAG_GRIDS, NTAU, ITYP, SIGMA, COMM ) 

! the following is necessary if grids are computed via minimizer
      ! scale interval properly, such that only right boundary is relevant
      ! the left interval boundary is set properly
      CALL NORMALIZE_INTERVAL( IMAG_GRIDS, X1, X2, IO )
       
      ! obtain time grid first 
      IF( IMAG_GRIDS%TIME%N > 0 ) THEN
         CALL DETERMINE_QUADRATURE( IMAG_GRIDS%TIME, IMAG_GRIDS%R, IO ) 
      ENDIF

      ! obtain grid for real part of bosonic correlation function
      IF( IMAG_GRIDS%FREQ_BOS_RE%N > 0 ) THEN
         CALL DETERMINE_QUADRATURE( IMAG_GRIDS%FREQ_BOS_RE, IMAG_GRIDS%R,&
            IO, IMAG_GRIDS%COMM_INTER ) 
      ENDIF

      ! obtain grid for imaginary part of bosonic correlation function
      IF( IMAG_GRIDS%FREQ_BOS_IM%N > 0 ) THEN
         CALL DETERMINE_QUADRATURE( IMAG_GRIDS%FREQ_BOS_IM, IMAG_GRIDS%R,&
            IO, IMAG_GRIDS%COMM_INTER ) 
      ENDIF

      ! obtain grid for real part of fermionic correlation function
      IF( IMAG_GRIDS%FREQ_FER_RE%N > 0 ) THEN
         CALL DETERMINE_QUADRATURE( IMAG_GRIDS%FREQ_FER_RE, IMAG_GRIDS%R,&
            IO, IMAG_GRIDS%COMM_INTER ) 
      ENDIF

      ! obtain grid for imaginary part of fermionic correlation function
      IF( IMAG_GRIDS%FREQ_FER_IM%N > 0 ) THEN
         CALL DETERMINE_QUADRATURE( IMAG_GRIDS%FREQ_FER_IM, IMAG_GRIDS%R,&
            IO, IMAG_GRIDS%COMM_INTER ) 
      ENDIF

      ! sync quadratures is necessary 
      CALL SYNC_QUADRATURES( IMAG_GRIDS )
      ! -----------------------------------------------------
      ! Computation of transformation matrices 
      ! -----------------------------------------------------
      IF ( IMAG_GRIDS%TIME%N > 0 ) THEN
         ! from time to FREQ_BOS_RE grid 
         IF( IMAG_GRIDS%FREQ_BOS_RE%N >  0 ) THEN
            ! transformation to real correlation functions
            CALL  DETERMINE_TRANSFORMATION( IMAG_GRIDS%TIME, &
               IMAG_GRIDS%FREQ_BOS_RE, IMAG_GRIDS%TO_BOS_RE, &
               IMAG_GRIDS%TO_BOS_RE_ERROR, .FALSE.,&
               COMM, IO ) 
            ! transformation to imaginary correlation functions
            CALL  DETERMINE_TRANSFORMATION( IMAG_GRIDS%TIME, &
               IMAG_GRIDS%FREQ_BOS_RE, IMAG_GRIDS%TO_BOS_RE_CONJG, &
               IMAG_GRIDS%TO_BOS_RE_CONJG_ERROR, .TRUE.,&
               COMM, IO ) 
         ENDIF 

         ! from time to FREQ_BOS_IM grid
         IF( IMAG_GRIDS%FREQ_BOS_IM%N >  0 ) THEN
            ! transformation to real correlation functions
            CALL  DETERMINE_TRANSFORMATION( IMAG_GRIDS%TIME, &
               IMAG_GRIDS%FREQ_BOS_IM, IMAG_GRIDS%TO_BOS_IM, &
               IMAG_GRIDS%TO_BOS_IM_ERROR, .FALSE.,&
               COMM, IO ) 
            ! transformation to imaginary correlation functions
            CALL  DETERMINE_TRANSFORMATION( IMAG_GRIDS%TIME, &
               IMAG_GRIDS%FREQ_BOS_IM, IMAG_GRIDS%TO_BOS_IM_CONJG, &
               IMAG_GRIDS%TO_BOS_IM_CONJG_ERROR, .TRUE.,&
               COMM, IO ) 
         ENDIF 

         ! from time to FREQ_FER_RE grid
         IF( IMAG_GRIDS%FREQ_FER_RE%N >  0 ) THEN
            ! transformation to real correlation functions
            CALL  DETERMINE_TRANSFORMATION( IMAG_GRIDS%TIME, &
               IMAG_GRIDS%FREQ_FER_RE, IMAG_GRIDS%TO_FER_RE, &
               IMAG_GRIDS%TO_FER_RE_ERROR, .FALSE.,&
               COMM, IO ) 
            ! transformation to imaginary correlation functions
            CALL  DETERMINE_TRANSFORMATION( IMAG_GRIDS%TIME, &
               IMAG_GRIDS%FREQ_FER_RE, IMAG_GRIDS%TO_FER_RE_CONJG, &
               IMAG_GRIDS%TO_FER_RE_CONJG_ERROR, .TRUE.,&
               COMM, IO ) 
         ENDIF 

         ! from time to FREQ_FER_IM grid
         IF( IMAG_GRIDS%FREQ_FER_IM%N >  0 ) THEN
            ! transformation to real correlation functions
            CALL  DETERMINE_TRANSFORMATION( IMAG_GRIDS%TIME, &
               IMAG_GRIDS%FREQ_FER_IM, IMAG_GRIDS%TO_FER_IM, &
               IMAG_GRIDS%TO_FER_IM_ERROR, .FALSE.,&
               COMM, IO ) 
            ! transformation to imaginary correlation functions
            CALL  DETERMINE_TRANSFORMATION( IMAG_GRIDS%TIME, &
               IMAG_GRIDS%FREQ_FER_IM, IMAG_GRIDS%TO_FER_IM_CONJG, &
               IMAG_GRIDS%TO_FER_IM_CONJG_ERROR, .TRUE.,&
               COMM, IO ) 

            ! this is a dirty fix to reproduce the T=0 GWR results
            ! that conserve the number of particles for insulators
            ! inverse trafo to t=0 from frequency is not computed
            ! has effects in acont.F
            IF ( IMAG_GRIDS%LFINITE_TEMPERATURE ) THEN
               CALL COPY_QUADRATURE( IMAG_GRIDS%TIME, TAU_TMP, .TRUE. ) 
               ! transformation from frequency grid to tau=0+point
               CALL  DETERMINE_TRANSFORMATION( IMAG_GRIDS%FREQ_FER_RE, &
                  TAU_TMP, IMAG_GRIDS%TO_TAU0, &
                  IMAG_GRIDS%TO_TAU0_ERROR, .FALSE.,&
                  COMM, IO ) 
               ! transformation from frequency grid to tau=0+point
               CALL  DETERMINE_TRANSFORMATION( IMAG_GRIDS%FREQ_FER_RE, &
                  TAU_TMP, IMAG_GRIDS%TO_TAU0_CONJG, &
                  IMAG_GRIDS%TO_TAU0_CONJG_ERROR, .TRUE.,&
                  COMM, IO ) 
               CALL KILL_QUADRATURE( TAU_TMP ) 
            ENDIF 
         ENDIF 

      ENDIF
  
      ! copies quadrature handles to imag grid
      CALL SCALE_QUADRATURES( IMAG_GRIDS, IO ) 

      ! frees all quadrature handles in imaginary grid handle 
      CALL RELEASE_QUADRATURES( IMAG_GRIDS, IO ) 

      PROFILING_STOP('SET_IMAG_GRID_HANDLE')

      CONTAINS

      !*******************************************************************
      !>helper, translates requested interval to normalized one
      !*******************************************************************
      SUBROUTINE NORMALIZE_INTERVAL( GRIDS, X1, X2, IO ) 
         USE tutor, ONLY: vtutor, isAlert, OmegaGridOMEGAMIN
         TYPE( imag_grid_handle ) :: GRIDS
         REAL(q),INTENT(IN)       :: X1, X2     ! energy interval 
         TYPE(in_struct)          :: IO         ! for output and input
          
         PROFILING_START( 'NORMALIZE_INTERVAL' )
         ! first determine R and set X1, X2, the interval requested
         GRIDS%X1=X1
         GRIDS%X2=X2
         IF ( GRIDS%X1 < 0.00001_q ) THEN
            CALL vtutor%write(isAlert, OmegaGridOMEGAMIN)
            GRIDS%X1=MAX(0.00001_q,GRIDS%X1)
         ENDIF
   
         ! set A to very small value for finite temperature grids 
         IF ( GRIDS%LFINITE_TEMPERATURE ) THEN
            GRIDS%X1 = 0
            GRIDS%R=GRIDS%BETA*GRIDS%X2
         ELSE
            GRIDS%R=GRIDS%X2/GRIDS%X1
         ENDIF
   
         
         !determine, which coefficients to use
         IF (IO%IU6>=0) WRITE(IO%IU6,1)GRIDS%X1,GRIDS%X2
         IF (IO%IU0>=0) WRITE(IO%IU0,1)GRIDS%X1,GRIDS%X2
   
         IF ( IO%IU0 >=0 .AND. IO%NWRITE > 2 ) THEN
            IF ( GRIDS%LFINITE_TEMPERATURE ) THEN
               WRITE(IO%IU0, '(A,E10.3,A)')'  which corresponds to [0,',GRIDS%R,' ]'
            ELSE
               WRITE(IO%IU0, '(A,E10.3,A)')'  which corresponds to [1,',GRIDS%R,' ]'
            ENDIF
         ENDIF 
1        FORMAT(' quadrature errors minimized for energies in [',E10.3,',',E10.3,' ]')
         PROFILING_STOP( 'NORMALIZE_INTERVAL' )
      END SUBROUTINE NORMALIZE_INTERVAL


   END SUBROUTINE SET_IMAG_GRID_HANDLE

!****************************************************************************
! deacllocted grid handle 
!****************************************************************************
   SUBROUTINE DEALLOCATE_IMAG_GRID_HANDLE(GRIDS)
      TYPE( imag_grid_handle ) :: GRIDS
      PROFILING_START( 'DEALLOCATE_IMAG_GRID_HANDLE')
 
      ! release loop descriptors for time
      CALL RELEASE_LOOP_DES( GRIDS%T )     
      ! bosonic frequencies
      CALL RELEASE_LOOP_DES( GRIDS%B )     
      ! fermionic frequencies
      CALL RELEASE_LOOP_DES( GRIDS%F )     

      ! reset member variables of structure
      GRIDS%NOMEGA = 0 
      ! no finite temperature grid calculated
      GRIDS%LFINITE_TEMPERATURE = .FALSE.
      ! finite temperature value back to room temperature
      GRIDS%BETA = 38.68149168875351_q
      ! left minimization interval boundary
      GRIDS%X1 = 0 
      ! right minimization interval boundary
      GRIDS%X2 = 0 
      ! right minimization interval boundary (scaled)
      GRIDS%R = 0 
      ! auxilary factor for CRPAR 
      GRIDS%FACTOR = 1._q
      
      ! deallocate associated time points
      IF ( ASSOCIATED( GRIDS%TAU ) ) DEALLOCATE( GRIDS%TAU )
      IF ( ASSOCIATED( GRIDS%TAU_WEIGHT ) ) DEALLOCATE( GRIDS%TAU_WEIGHT )
      GRIDS%TAU_ERROR = 0 
      IF ( ASSOCIATED( GRIDS%TO_TAU0 ) ) DEALLOCATE( GRIDS%TO_TAU0 )
      IF ( ASSOCIATED( GRIDS%TO_TAU0_ERROR ) ) DEALLOCATE( GRIDS%TO_TAU0_ERROR )
      IF ( ASSOCIATED( GRIDS%TO_TAU0_CONJG ) ) DEALLOCATE( GRIDS%TO_TAU0_CONJG )
      IF ( ASSOCIATED( GRIDS%TO_TAU0_CONJG_ERROR ) ) DEALLOCATE( GRIDS%TO_TAU0_CONJG_ERROR )

      ! deallocate associated bosonic frequencies for real part 
      IF ( ASSOCIATED( GRIDS%BOS_RE ) ) DEALLOCATE( GRIDS%BOS_RE ) 
      IF ( ASSOCIATED( GRIDS%BOS_RE_WEIGHT ) ) DEALLOCATE( GRIDS%BOS_RE_WEIGHT ) 
      GRIDS%BOS_RE_ERROR = 0     
      IF ( ASSOCIATED( GRIDS%TO_BOS_RE ) ) DEALLOCATE( GRIDS%TO_BOS_RE ) 
      IF ( ASSOCIATED( GRIDS%TO_BOS_RE_ERROR ) ) DEALLOCATE( GRIDS%TO_BOS_RE_ERROR ) 
      IF ( ASSOCIATED( GRIDS%TO_BOS_RE_CONJG ) ) DEALLOCATE( GRIDS%TO_BOS_RE_CONJG ) 
      IF ( ASSOCIATED( GRIDS%TO_BOS_RE_CONJG_ERROR ) ) DEALLOCATE( GRIDS%TO_BOS_RE_CONJG_ERROR ) 
      
      ! deallocate associated bosonic frequencies for imaginary part 
      IF ( ASSOCIATED( GRIDS%BOS_IM ) ) DEALLOCATE( GRIDS%BOS_IM ) 
      IF ( ASSOCIATED( GRIDS%BOS_IM_WEIGHT ) ) DEALLOCATE( GRIDS%BOS_IM_WEIGHT ) 
      GRIDS%BOS_IM_ERROR = 0     
      IF ( ASSOCIATED( GRIDS%TO_BOS_IM ) ) DEALLOCATE( GRIDS%TO_BOS_IM ) 
      IF ( ASSOCIATED( GRIDS%TO_BOS_IM_ERROR ) ) DEALLOCATE( GRIDS%TO_BOS_IM_ERROR ) 
      IF ( ASSOCIATED( GRIDS%TO_BOS_IM_CONJG ) ) DEALLOCATE( GRIDS%TO_BOS_IM_CONJG ) 
      IF ( ASSOCIATED( GRIDS%TO_BOS_IM_CONJG_ERROR ) ) DEALLOCATE( GRIDS%TO_BOS_IM_CONJG_ERROR ) 

      ! deallocate associated fermionic frequencies for real part 
      IF ( ASSOCIATED( GRIDS%FER_RE ) ) DEALLOCATE( GRIDS%FER_RE ) 
      IF ( ASSOCIATED( GRIDS%FER_RE_WEIGHT ) ) DEALLOCATE( GRIDS%FER_RE_WEIGHT ) 
      GRIDS%FER_RE_ERROR = 0     
      IF ( ASSOCIATED( GRIDS%TO_FER_RE ) ) DEALLOCATE( GRIDS%TO_FER_RE ) 
      IF ( ASSOCIATED( GRIDS%TO_FER_RE_ERROR ) ) DEALLOCATE( GRIDS%TO_FER_RE_ERROR ) 
      IF ( ASSOCIATED( GRIDS%TO_FER_RE_CONJG ) ) DEALLOCATE( GRIDS%TO_FER_RE_CONJG ) 
      IF ( ASSOCIATED( GRIDS%TO_FER_RE_CONJG_ERROR ) ) DEALLOCATE( GRIDS%TO_FER_RE_CONJG_ERROR ) 

      ! deallocate associated fermionic frequencies for imaginary part 
      IF ( ASSOCIATED( GRIDS%FER_IM ) ) DEALLOCATE( GRIDS%FER_IM ) 
      IF ( ASSOCIATED( GRIDS%FER_IM_WEIGHT ) ) DEALLOCATE( GRIDS%FER_IM_WEIGHT ) 
      GRIDS%FER_IM_ERROR = 0     
      IF ( ASSOCIATED( GRIDS%TO_FER_IM ) ) DEALLOCATE( GRIDS%TO_FER_IM ) 
      IF ( ASSOCIATED( GRIDS%TO_FER_IM_ERROR ) ) DEALLOCATE( GRIDS%TO_FER_IM_ERROR ) 
      IF ( ASSOCIATED( GRIDS%TO_FER_IM_CONJG ) ) DEALLOCATE( GRIDS%TO_FER_IM_CONJG ) 
      IF ( ASSOCIATED( GRIDS%TO_FER_IM_CONJG_ERROR ) ) DEALLOCATE( GRIDS%TO_FER_IM_CONJG_ERROR ) 

      PROFILING_STOP( 'DEALLOCATE_IMAG_GRID_HANDLE')
      CONTAINS

      ! releases a loop_des type
      SUBROUTINE RELEASE_LOOP_DES( DES ) 
         TYPE( loop_des ) :: DES 
         PROFILING_START( 'release_loop_des' )
         IF ( ASSOCIATED( DES%POINTS_LOCAL ) ) DEALLOCATE( DES%POINTS_LOCAL )
         IF ( ASSOCIATED( DES%GROUP_OF_NODE ) ) DEALLOCATE( DES%GROUP_OF_NODE )
         IF ( ASSOCIATED( DES%DISTRIBUTION ) ) THEN
            DEALLOCATE( DES%DISTRIBUTION )
            CALLMPI( M_freec( DES%COMM_IN_GROUP ))
            CALLMPI( M_freec( DES%COMM_BETWEEN_GROUPS ))
         ENDIF
         PROFILING_STOP( 'release_loop_des' )
      END SUBROUTINE RELEASE_LOOP_DES
    
   END SUBROUTINE DEALLOCATE_IMAG_GRID_HANDLE  

!****************************************************************************
!>helper, initializes basic quadrature handles in imag_grid_handle
!****************************************************************************

   SUBROUTINE INITIALIZE_QUADRATURES( IMAG_GRIDS, NTAU, ITYP, SIGMA, COMM ) 
      USE minimax_functions1D
      USE minimax_functions2D
      USE minimax_ini
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      TYPE( imag_grid_handle ) :: IMAG_GRIDS
      INTEGER                  :: NTAU
      INTEGER                  :: ITYP 
      REAL(q)                  :: SIGMA 
      TYPE( communic )         :: COMM    ! global communicator 
      ! local                 
      INTEGER                  :: NGRIDS  ! counts the number of distinct grids requested

      PROFILING_START( 'INITIALIZE_QUADRATURES' ) 

      NGRIDS = 0 

      ! no finite temperature grid by default
      IMAG_GRIDS%LFINITE_TEMPERATURE = .FALSE.
      IMAG_GRIDS%NOMEGA = NTAU 
       
      ! initialize zero infinity norm 
      IMAG_GRIDS%TIME%INFINITYNORM= QDPP( 0 )
      IMAG_GRIDS%FREQ_BOS_RE%INFINITYNORM= QDPP( 0 )
      IMAG_GRIDS%FREQ_BOS_IM%INFINITYNORM= QDPP( 0 )
      IMAG_GRIDS%FREQ_FER_RE%INFINITYNORM= QDPP( 0 )
      IMAG_GRIDS%FREQ_FER_IM%INFINITYNORM= QDPP( 0 )
      IMAG_GRIDS%TIME%SCALING= QDPP( 1 )
      IMAG_GRIDS%FREQ_BOS_RE%SCALING= QDPP( 1 )
      IMAG_GRIDS%FREQ_BOS_IM%SCALING= QDPP( 1 )
      IMAG_GRIDS%FREQ_FER_RE%SCALING= QDPP( 1 )
      IMAG_GRIDS%FREQ_FER_IM%SCALING= QDPP( 1 )

      IMAG_GRIDS%TIME%R1= QDPP( 0 )
      IMAG_GRIDS%TIME%R2= QDPP( 0 )
      IMAG_GRIDS%FREQ_BOS_RE%R1= QDPP( 0 )
      IMAG_GRIDS%FREQ_BOS_RE%R2= QDPP( 0 )
      IMAG_GRIDS%FREQ_BOS_IM%R1= QDPP( 0 )
      IMAG_GRIDS%FREQ_BOS_IM%R2= QDPP( 0 )
      IMAG_GRIDS%FREQ_FER_RE%R1= QDPP( 0 )
      IMAG_GRIDS%FREQ_FER_RE%R2= QDPP( 0 )
      IMAG_GRIDS%FREQ_FER_IM%R1= QDPP( 0 )
      IMAG_GRIDS%FREQ_FER_IM%R2= QDPP( 0 )

      IMAG_GRIDS%BETA = 1/SIGMA
      IMAG_GRIDS%T%BETA = 1/SIGMA
      IMAG_GRIDS%B%BETA = 1/SIGMA
      IMAG_GRIDS%F%BETA = 1/SIGMA

      IF ( ITYP == 140 ) THEN
         ! requires frequency grid for real part of bosonic functions
         CALL INIT_T0_RE( NTAU, IMAG_GRIDS%FREQ_BOS_RE )

      ELSE IF ( ITYP == 141 ) THEN
         ! requires time grid 
         CALL INIT_T0_TIME( NTAU, IMAG_GRIDS%TIME ) 

      ELSE IF ( ITYP == 145 ) THEN
         ! requires time grid 
         CALL INIT_T0_TIME( NTAU, IMAG_GRIDS%TIME ) 
         ! requires BOS_RE
         CALL INIT_T0_RE( NTAU, IMAG_GRIDS%FREQ_BOS_RE )
         ! requires BOS_IM
!         CALL INIT_T0_IM( NTAU, IMAG_GRIDS%FREQ_BOS_IM )
         ! requires FER_RE
         CALL INIT_T0_RE2( NTAU, IMAG_GRIDS%FREQ_FER_RE )
         ! requires FER_IM
         CALL INIT_T0_IM( NTAU, IMAG_GRIDS%FREQ_FER_IM )

      ELSE IF ( ITYP == 150 ) THEN
         ! requires finite temperature BOS_RE
         CALL INIT_BOS_RE( NTAU, IMAG_GRIDS%FREQ_BOS_RE )

      ELSE IF ( ITYP == 151 ) THEN
         ! requires finite temperature time grid
         CALL INIT_T_TIME2( NTAU, IMAG_GRIDS%TIME )

      ELSE IF ( ITYP == 152 ) THEN
         ! requires BOS_RE
         CALL INIT_BOS_RE( NTAU, IMAG_GRIDS%FREQ_BOS_RE )
      ELSE IF ( ITYP == 153 ) THEN
         ! requires all finite temperature grids
         CALL INIT_T_TIME( NTAU, IMAG_GRIDS%TIME )
         ! requires BOS_RE
         !CALL INIT_BOS_RE_FOR_MP1( NTAU, IMAG_GRIDS%FREQ_BOS_RE )
         CALL INIT_BOS_RE( NTAU, IMAG_GRIDS%FREQ_BOS_RE )
         ! requires FER_RE
         CALL INIT_FER_RE_FOR_GM( NTAU, IMAG_GRIDS%FREQ_FER_RE )
         ! requires FER_IM
         CALL INIT_FER_IM( NTAU, IMAG_GRIDS%FREQ_FER_IM )
      ELSE IF ( ITYP == 154 ) THEN
         ! requires FER_IM
         !CALL INIT_FER_IM( NTAU, IMAG_GRIDS%FREQ_FER_IM )
         CALL INIT_FER_RE_FOR_GM( NTAU, IMAG_GRIDS%FREQ_FER_RE )

      ELSE IF ( ITYP == 155 ) THEN
         ! requires all finite temperature grids
         CALL INIT_T_TIME2( NTAU, IMAG_GRIDS%TIME )
         ! requires BOS_RE
         !CALL INIT_BOS_RE_FOR_MP1( NTAU, IMAG_GRIDS%FREQ_BOS_RE )
         CALL INIT_BOS_RE( NTAU, IMAG_GRIDS%FREQ_BOS_RE )
!         CALL INIT_BOS_IM( NTAU, IMAG_GRIDS%FREQ_BOS_IM )
         ! requires FER_RE
         CALL INIT_FER_RE_FOR_GM( NTAU, IMAG_GRIDS%FREQ_FER_RE )
         ! requires FER_IM
         !CALL INIT_FER_IM( NTAU, IMAG_GRIDS%FREQ_FER_IM )
         CALL INIT_FER_IM_FOR_GM( NTAU, IMAG_GRIDS%FREQ_FER_IM )

      ELSE IF ( ITYP == 156 ) THEN
         ! requires all finite temperature grids
         CALL INIT_T_TIME2( NTAU, IMAG_GRIDS%TIME )
         ! requires BOS_RE
         CALL INIT_BOS_RE( NTAU, IMAG_GRIDS%FREQ_BOS_RE )
         ! requires BOS_IM
!         CALL INIT_BOS_IM( NTAU, IMAG_GRIDS%FREQ_BOS_IM )
         ! requires FER_RE
         CALL INIT_FER_RE( NTAU, IMAG_GRIDS%FREQ_FER_RE )
         ! requires FER_IM
         CALL INIT_FER_IM( NTAU, IMAG_GRIDS%FREQ_FER_IM )

      ELSE IF ( ITYP == 157 ) THEN
         ! requires all finite temperature grids
         CALL INIT_T_TIME2( NTAU, IMAG_GRIDS%TIME )
         ! requires BOS_RE
         CALL INIT_BOS_RE( NTAU, IMAG_GRIDS%FREQ_BOS_RE )
         ! exchange RE and IM grid       
         ! requires FER_RE
         CALL INIT_FER_IM( NTAU, IMAG_GRIDS%FREQ_FER_RE )
         ! requires FER_IM
         CALL INIT_FER_RE( NTAU, IMAG_GRIDS%FREQ_FER_IM )

      ELSE IF ( ITYP == 158 ) THEN
         ! requires all finite temperature grids
         CALL INIT_T_TIME2( NTAU, IMAG_GRIDS%TIME )
         ! requires BOS_RE
         CALL INIT_BOS_RE_FOR_MP1( NTAU, IMAG_GRIDS%FREQ_BOS_RE )
         ! requires FER_RE
         CALL INIT_FER_RE_FOR_GM( NTAU, IMAG_GRIDS%FREQ_FER_RE )

      ELSE IF ( ITYP == 161 ) THEN
         CALL INIT_BOS_RE( NTAU, IMAG_GRIDS%FREQ_BOS_RE )
      ELSE IF ( ITYP == 162 ) THEN
         CALL INIT_FER_RE( NTAU, IMAG_GRIDS%FREQ_FER_RE )
      ELSE IF ( ITYP == 163 ) THEN
         CALL INIT_FER_RE_FOR_GM( NTAU, IMAG_GRIDS%FREQ_FER_RE )
      ! Matsubara grid 
      ELSE IF ( ITYP == 168 )  THEN
         ! requires finite temperature time grid
         CALL INIT_HYPERGEOM( NTAU, IMAG_GRIDS )
      ELSE IF ( ITYP == 169 )  THEN
         ! requires finite temperature time grid
         CALL INIT_MATSUBARA( NTAU, IMAG_GRIDS )

      ELSE
         CALL vtutor%error("ERROR, requested grid type not supported " // str(ITYP))
      ENDIF
  
      IF (  IMAG_GRIDS%LFINITE_TEMPERATURE .AND. &
         (SIGMA < 0.0001_q .OR. SIGMA > 10._q )) THEN
         CALL vtutor%error("ERROR, inverse temperature out of bounds " // str(SIGMA))
      ENDIF
      ! set temperature of grid
      IF ( IMAG_GRIDS%LFINITE_TEMPERATURE ) THEN
         IMAG_GRIDS%BETA = 1/SIGMA
      ENDIF 

#ifdef MPI
      IF ( COMM%NODE_ME == 1 ) THEN
         WRITE(*,'(A,I8)')' number of imaginary grid points requested:',NTAU
         WRITE(*,'(A,I3)')' number of distinct grids requested:',NGRIDS 
      ENDIF
#endif
   
      ! split communicators into groups 
      ! note that time grid is calculated by all Ranks
      ! frequency grid calculation might be distributed into groups 
      CALL SPLIT_COMMUNICATOR( COMM, NGRIDS, IMAG_GRIDS )
      PROFILING_STOP( 'INITIALIZE_QUADRATURES' ) 

      CONTAINS

      ! *********************************************
      ! intitializes T=0 time grid 
      ! *********************************************
      SUBROUTINE INIT_T0_TIME( NTAU, QUAD ) 
         INTEGER  :: NTAU
         TYPE( quadrature_handle )  :: QUAD
         PROFILING_START('init_t0_time')
         ! increase distinct grid counter
         NGRIDS = NGRIDS + 1 

         ! set number of basis functions
         QUAD%N = NTAU 
         ! set identifier
         QUAD%ITYPE       = 1
         ! no dual basis function is set for time quadrature

         ! basis function used for quadrature error
         QUAD%PHI2        => EXPF2
         ! first derivative w.r.t. first argument of basis function
         QUAD%D_PHI2_DZ   => D_EXPF2_DZ
         ! first derivative w.r.t. second argument of basis function
         QUAD%D_PHI2_DL   => D_EXPF2_DL
         ! second derivative w.r.t. first argument of basis function
         QUAD%D2_PHI2_DZ2 => D2_EXPF2_DZ2
         ! object function that is being fitted 
         QUAD%F           => INVERSE_Z
         ! first derivative of object function
         QUAD%D_F_DZ      => D1_INVERSE_Z
         ! second derivative of object function
         QUAD%D2_F_DZ2    => D2_INVERSE_Z
         ! quadrature scaling factor is band gap
         QUAD%SCALING     = QDPP(1)/QDPD( IMAG_GRIDS%X1 )
         ! grid id for parallelization
         QUAD%GRID_ID = NGRIDS-1

         ! set basis functions to exp(-z t) 
         QUAD%PHI         => EXPF
         ! set basis functions to ( x/(x^2+w^2) )
         QUAD%PSI         => UFREQ
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR=> COS_TAYLOR

         ! conjugate transformation error function
         ! set basis functions to ( w/(x^2+w^2) )
         QUAD%PHI_CONJG   => EXPF
         ! conjugate basis function set for time quadrature
         QUAD%PSI_CONJG   => VFREQ
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR_CONJG=> SIN_TAYLOR
        
         PROFILING_STOP('init_t0_time')
      END SUBROUTINE INIT_T0_TIME
 
      ! *********************************************
      ! intitializes error function T=0 real part
      ! *********************************************
      SUBROUTINE INIT_T0_RE( NOMEGA, QUAD ) 
         INTEGER                    :: NOMEGA
         TYPE( quadrature_handle )  :: QUAD 
         PROFILING_START('init_t0_re')
         ! increase distinct grid counter
         NGRIDS = NGRIDS + 1 

         ! set number of basis functions
         QUAD%N = NOMEGA
         ! set identifier
         QUAD%ITYPE       = 2
         ! basis function used for quadrature error
         QUAD%PHI2        => UFREQ2
         ! first derivative w.r.t. first argument of basis function
         QUAD%D_PHI2_DZ   => D_UFREQ2_DZ
         ! first derivative w.r.t. second argument of basis function
         QUAD%D_PHI2_DL   => D_UFREQ2_DL
         ! second derivative w.r.t. first argument of basis function
         QUAD%D2_PHI2_DZ2 => D2_UFREQ2_DZ2
         ! object function that is being fitted 
         QUAD%F           => INVERSE_Z
         ! first derivative of object function
         QUAD%D_F_DZ      => D1_INVERSE_Z
         ! second derivative of object function
         QUAD%D2_F_DZ2    => D2_INVERSE_Z
         ! quadrature scaling factor is band gap
         QUAD%SCALING     = QDPD( IMAG_GRIDS%X1 )
         ! grid id for parallelization
         QUAD%GRID_ID = NGRIDS-1

         ! set basis functions to ( x/(x^2+w^2) )
         QUAD%PHI         => UFREQ
         ! dual basis function set for time quadrature
         QUAD%PSI         => EXPF
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR=> COS_TAYLOR
         
         ! conjugate error
         ! set basis functions to ( w/(x^2+w^2) )
         QUAD%PHI_CONJG   => VFREQ
         ! dual basis function set for time quadrature
         QUAD%PSI_CONJG   => EXPF
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR_CONJG=> SIN_TAYLOR
        
         PROFILING_STOP('init_t0_re')
      END SUBROUTINE INIT_T0_RE

      ! *********************************************
      ! intitializes error function T=0 real part
      ! *********************************************
      SUBROUTINE INIT_T0_RE2( NOMEGA, QUAD ) 
         INTEGER                    :: NOMEGA
         TYPE( quadrature_handle )  :: QUAD 
         PROFILING_START('init_t0_re2')
         ! increase distinct grid counter
         NGRIDS = NGRIDS + 1 

         ! set number of basis functions
         QUAD%N = NOMEGA
         ! set identifier
         QUAD%ITYPE       = 20
         ! basis function used for quadrature error
         QUAD%PHI2        => UFREQ2
         ! first derivative w.r.t. first argument of basis function
         QUAD%D_PHI2_DZ   => D_UFREQ2_DZ
         ! first derivative w.r.t. second argument of basis function
         QUAD%D_PHI2_DL   => D_UFREQ2_DL
         ! second derivative w.r.t. first argument of basis function
         QUAD%D2_PHI2_DZ2 => D2_UFREQ2_DZ2
         ! object function that is being fitted 
         QUAD%F           => INVERSE_Z
         ! first derivative of object function
         QUAD%D_F_DZ      => D1_INVERSE_Z
         ! second derivative of object function
         QUAD%D2_F_DZ2    => D2_INVERSE_Z
         ! quadrature scaling factor is band gap
         QUAD%SCALING     = QDPD( IMAG_GRIDS%X1 )
         ! grid id for parallelization
         QUAD%GRID_ID = NGRIDS-1

         ! set basis functions to ( x/(x^2+w^2) )
         QUAD%PHI         => UFREQ
         ! dual basis function set for time quadrature
         QUAD%PSI         => EXPF
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR=> COS_TAYLOR
         
         ! conjugate error
         ! set basis functions to ( w/(x^2+w^2) )
         QUAD%PHI_CONJG   => VFREQ
         ! dual basis function set for time quadrature
         QUAD%PSI_CONJG   => EXPF
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR_CONJG=> SIN_TAYLOR
        
         PROFILING_STOP('init_t0_re2')
      END SUBROUTINE INIT_T0_RE2

      ! *********************************************
      ! intitializes error function T=0 imag. part
      ! *********************************************
      SUBROUTINE INIT_T0_IM( NOMEGA, QUAD ) 
         INTEGER                    :: NOMEGA
         TYPE( quadrature_handle )  :: QUAD 
         PROFILING_START('init_t0_im')
         ! increase distinct grid counter
         NGRIDS = NGRIDS + 1 
         ! set number of basis functions
         QUAD%N = NOMEGA
         ! set identifier
         QUAD%ITYPE       = 3
         ! basis function used for quadrature error
         QUAD%PHI2        => VFREQ2
         ! first derivative w.r.t. first argument of basis function
         QUAD%D_PHI2_DZ   => D_VFREQ2_DZ
         ! first derivative w.r.t. second argument of basis function
         QUAD%D_PHI2_DL   => D_VFREQ2_DL
         ! second derivative w.r.t. first argument of basis function
         QUAD%D2_PHI2_DZ2 => D2_VFREQ2_DZ2
         ! object function that is being fitted 
         QUAD%F           => INVERSE_Z
         ! first derivative of object function
         QUAD%D_F_DZ      => D1_INVERSE_Z
         ! second derivative of object function
         QUAD%D2_F_DZ2    => D2_INVERSE_Z
         ! quadrature scaling factor is band gap
         QUAD%SCALING     = QDPD( IMAG_GRIDS%X1 )
         ! grid id for parallelization
         QUAD%GRID_ID = NGRIDS-1

         ! set basis functions to ( w/(x^2+w^2) )
         QUAD%PHI         => VFREQ
         ! dual basis function set for time quadrature
         QUAD%PSI         => EXPF
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR=> SIN_TAYLOR

         ! conjugate transformation error function
         ! set basis functions to ( x/(x^2+w^2) )
         QUAD%PHI_CONJG   => UFREQ
         ! conjugate basis function set for time quadrature
         QUAD%PSI_CONJG   => EXPF
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR_CONJG=> COS_TAYLOR

         PROFILING_STOP('init_t0_im')
      END SUBROUTINE INIT_T0_IM

!
! T> 0 initiators
!
      ! *********************************************
      ! intitializes T>0 time grid 
      ! with phi as basis
      ! *********************************************
      SUBROUTINE INIT_T_TIME( NTAU, QUAD ) 
         INTEGER                    :: NTAU
         TYPE( quadrature_handle )  :: QUAD 
         PROFILING_START('init_t_time')

         ! increase distinct grid counter
         NGRIDS = NGRIDS + 1 
         ! flag this grid as a finite temperature grid 
         IMAG_GRIDS%LFINITE_TEMPERATURE = .TRUE.
         ! number of time grid points
         QUAD%N = NTAU 
         ! set identifier
         QUAD%ITYPE       = -10
         ! basis function used for quadrature error
         QUAD%PHI2        => UTIME
         ! first derivative w.r.t. first argument of basis function
         QUAD%D_PHI2_DZ   => D_UTIME_DZ
         ! first derivative w.r.t. second argument of basis function
         QUAD%D_PHI2_DL   => D_UTIME_DL
         ! second derivative w.r.t. first argument of basis function
         QUAD%D2_PHI2_DZ2 => D2_UTIME_DZ2
         ! object function that is being fitted 
         QUAD%F           => F_THZ_OVER_Z
         ! first derivative of object function
         QUAD%D_F_DZ      => D1_THZ_OVER_Z
         ! second derivative of object function
         QUAD%D2_F_DZ2    => D2_THZ_OVER_Z
         ! set initial coefficients
         QUAD%CTAB        => TIME_COSH
         ! calling function to set closest interval in table to requested one
!         QUAD%R1E3        => SET_COSH_R1E3
!         QUAD%R1E4        => SET_COSH_R1E4
!         QUAD%R1E5        => SET_COSH_R1E5
         QUAD%R1E6        => SET_COSH_R1E6
         ! quadrature scaling factor is temperature in eV (or inverse beta)
         QUAD%SCALING     = QDPP(1)/QDPD( SIGMA )
         ! grid id for parallelization
         QUAD%GRID_ID = NGRIDS-1

         ! set basis function to cosh(x/2(1-2t))/cosh(x/2)
         QUAD%PHI         => UTIME
         ! dual basis function is w/(x^2+w^2)
         QUAD%PSI         => VFREQ
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR=> COS_TAYLOR

         ! set basis function to sinh(x/2(1-2t))/cosh(x/2)
         QUAD%PHI_CONJG   => VTIME
         ! dual basis function is x/(x^2+w^2)
         QUAD%PSI_CONJG   => UFREQ
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR_CONJG=> SIN_TAYLOR
        
         PROFILING_STOP('init_t_time')
      END SUBROUTINE INIT_T_TIME

      ! *********************************************
      ! intitializes T>0 time grid 
      ! with phi^2 as basis
      ! *********************************************
      SUBROUTINE INIT_T_TIME2( NTAU, QUAD ) 
         USE minimax_dependence, ONLY:FIT_COEFF_COMP_COSH
         INTEGER                    :: NTAU
         TYPE( quadrature_handle )  :: QUAD 
         PROFILING_START('init_t_time2')

         ! increase distinct grid counter
         NGRIDS = NGRIDS + 1 
         ! flag this grid as a finite temperature grid 
         IMAG_GRIDS%LFINITE_TEMPERATURE = .TRUE.
         ! number of time grid points
         QUAD%N = NTAU 
         ! set identifier
         QUAD%ITYPE       = -1
         ! basis function used for quadrature error
         QUAD%PHI2        => UTIME2
         ! first derivative w.r.t. first argument of basis function
         QUAD%D_PHI2_DZ   => D_UTIME2_DZ
         ! first derivative w.r.t. second argument of basis function
         QUAD%D_PHI2_DL   => D_UTIME2_DL
         ! second derivative w.r.t. first argument of basis function
         QUAD%D2_PHI2_DZ2 => D2_UTIME2_DZ2
         ! object function that is being fitted 
         QUAD%F           => F_BETA_LRS
         ! first derivative of object function
         QUAD%D_F_DZ      => D1_F_BETA_LRS
         ! second derivative of object function
         QUAD%D2_F_DZ2    => D2_F_BETA_LRS
         ! set initial coefficients
         QUAD%CTAB        => TIME_COSH2
         ! calling function to set closest interval in table to requested one
         QUAD%R1E3        => SET_COSH2_R1E3
         QUAD%R1E4        => SET_COSH2_R1E4
         QUAD%R1E5        => SET_COSH2_R1E5
         QUAD%R1E6        => SET_COSH2_R1E6
         ! quadrature scaling factor is temperature in eV (or inverse beta)
         QUAD%SCALING     = QDPP(1)/QDPD( SIGMA )
         ! grid id for parallelization
         QUAD%GRID_ID = NGRIDS-1

         ! set basis function to cosh(x/2(1-2t))/cosh(x/2)
         QUAD%PHI         => UTIME
         ! dual basis function is w/(x^2+w^2)
         QUAD%PSI         => VFREQ
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR=> COS_TAYLOR

         ! set basis function to sinh(x/2(1-2t))/cosh(x/2)
         QUAD%PHI_CONJG   => VTIME
         ! dual basis function is x/(x^2+w^2)
         QUAD%PSI_CONJG   => UFREQ
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR_CONJG=> SIN_TAYLOR

         ! no fitting functions below order 6
         IF ( NTAU > 3 ) THEN
            ! associate function to calculate coefficients directly 
            QUAD%FITCOF => FIT_COEFF_COMP_COSH
         ENDIF
        
         PROFILING_STOP('init_t_time2')
      END SUBROUTINE INIT_T_TIME2

      ! *********************************************
      ! intitializes BOS_RE error function T>0 
      ! *********************************************
      SUBROUTINE INIT_BOS_RE( NOMEGA, QUAD ) 
         USE minimax_dependence, ONLY:FIT_COEFF_COMP_TANH2
         INTEGER                    :: NOMEGA
         TYPE( quadrature_handle )  :: QUAD 
         INTEGER                    :: R1(24),R2(24)
         DATA R1(6:24)/& 
         30,  & ! 6
         40,  & ! 7
         40,  & ! 8
         50,  & ! 9
         60,  & ! 10
         60,  & ! 11
         65,  & ! 12
         95,  & ! 13
         100, & ! 14
         105, & ! 15
         175, & ! 16
         190, & ! 17
         190, & ! 18
         230, & ! 19
         400, & ! 20
         400, & ! 21
         490, & ! 22
         590, & ! 23
         690/   ! 24
         DATA R2(6:24)/& 
         4E3,  &!  6
         4E3,  &!  7
         4E3,  &!  8
         4E3,  &!  9
         4E3,  &! 10
         4E3,  &! 11
         4E3,  &! 12
         3E4,  &! 13
         3E4,  &! 14
         3E4,  &! 15
         4E5,  &! 16
         4E5,  &! 17
         4E5,  &! 18
         5E5,  &! 19
         6E6,  &! 20
         6E6,  &! 21
         7E6,  &! 22
         7E6,  &! 23
         8E6/   ! 24
         PROFILING_START('init_bos_re')
         ! increase distinct grid counter
         NGRIDS = NGRIDS + 1 

         ! flag this grid as a finite temperature grid 
         IMAG_GRIDS%LFINITE_TEMPERATURE = .TRUE.
         ! quadrature contains no zero point, w=0 is added later on
         QUAD%N = NOMEGA - 1 
         ! set identifier
         QUAD%ITYPE       = -2
         ! basis function used for quadrature error
         QUAD%PHI2        => UTANH2
         ! first derivative w.r.t. first argument of basis function
         QUAD%D_PHI2_DZ   => D_UTANH2_DZ
         ! first derivative w.r.t. second argument of basis function
         QUAD%D_PHI2_DL   => D_UTANH2_DL
         ! second derivative w.r.t. first argument of basis function
         QUAD%D2_PHI2_DZ2 => D2_UTANH2_DZ2
         ! object function that is being fitted 
         QUAD%F           => F_BETA_LRS_REG
         ! first derivative of object function
         QUAD%D_F_DZ      => D1_F_BETA_LRS_REG
         ! second derivative of object function
         QUAD%D2_F_DZ2    => D2_F_BETA_LRS_REG
         ! set initial coefficients
         QUAD%CTAB        => FREQ_UTANH2
         ! calling function to set closest interval in table to requested one
         QUAD%R1E3        => SET_UTANH2_R1E3
         QUAD%R1E4        => SET_UTANH2_R1E4
         QUAD%R1E5        => SET_UTANH2_R1E5
         QUAD%R1E6        => SET_UTANH2_R1E6
         ! quadrature scaling factor is temperature in eV (or inverse beta)
         QUAD%SCALING     = QDPD( SIGMA )
         ! grid id for parallelization
         QUAD%GRID_ID = NGRIDS-1

         ! set basis function to x/(x^2+w^2) tanh(x/2)
         QUAD%PHI         => UTANH
         ! dual basis function is cosh(x/2(1-2))/cosh(x/2) 
         QUAD%PSI         => UTIME
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR=> COS_TAYLOR
         
         ! conjugate error
         ! set basis functions to ( w/(x^2+w^2) tanh(x/2) )
         QUAD%PHI_CONJG   => VTANH
         ! dual basis function is sinh(x/2(1-2))/cosh(x/2) 
         QUAD%PSI_CONJG   => VTIME
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR_CONJG=> SIN_TAYLOR

         ! no fitting functions below order 6
         IF ( QUAD%N > 2 ) THEN
            ! associate function to calculate coefficients directly 
            QUAD%FITCOF => FIT_COEFF_COMP_TANH2
         ENDIF
        
         PROFILING_STOP('init_bos_re')
      END SUBROUTINE INIT_BOS_RE

      ! *********************************************
      ! intitializes BOS_RE error function T>0 
      ! *********************************************
      SUBROUTINE INIT_BOS_RE_FOR_MP1( NOMEGA, QUAD ) 
         INTEGER                    :: NOMEGA
         TYPE( quadrature_handle )  :: QUAD 
         PROFILING_START('init_bos_re_for_mp1')
         ! increase distinct grid counter
         NGRIDS = NGRIDS + 1 

         ! flag this grid as a finite temperature grid 
         IMAG_GRIDS%LFINITE_TEMPERATURE = .TRUE.
         ! quadrature contains no zero point, w=0 is added later on
         QUAD%N = NOMEGA - 1 
         ! set identifier
         QUAD%ITYPE       = -20
         ! basis function used for quadrature error
         QUAD%PHI2        => UTANH
         ! first derivative w.r.t. first argument of basis function
         QUAD%D_PHI2_DZ   => D_UTANH_DZ
         ! first derivative w.r.t. second argument of basis function
         QUAD%D_PHI2_DL   => D_UTANH_DL
         ! second derivative w.r.t. first argument of basis function
         QUAD%D2_PHI2_DZ2 => D2_UTANH_DZ2
         ! object function that is being fitted 
         QUAD%F           => COTH_REG
         ! first derivative of object function
         QUAD%D_F_DZ      => D1_COTH_REG
         ! second derivative of object function
         QUAD%D2_F_DZ2    => D2_COTH_REG
         ! set initial coefficients
         QUAD%CTAB        => FREQ_UCOTH
         ! calling function to set closest interval in table to requested one
         QUAD%R1E3        => SET_UCOTH_R1E3
         QUAD%R1E4        => SET_UCOTH_R1E4
         QUAD%R1E5        => SET_UCOTH_R1E5
         QUAD%R1E6        => SET_UCOTH_R1E6
         ! quadrature scaling factor is temperature in eV (or inverse beta)
         QUAD%SCALING     = QDPD( SIGMA )
         ! grid id for parallelization
         QUAD%GRID_ID = NGRIDS-1

         ! set basis function to x/(x^2+w^2) tanh(x/2)
         QUAD%PHI         => UTANH
         ! dual basis function is cosh(x/2(1-2))/cosh(x/2) 
         QUAD%PSI         => UTIME
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR=> COS_TAYLOR
         
         ! conjugate error
         ! set basis functions to ( w/(x^2+w^2) tanh(x/2) )
         QUAD%PHI_CONJG   => VTANH
         ! dual basis function is sinh(x/2(1-2))/cosh(x/2) 
         QUAD%PSI_CONJG   => VTIME
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR_CONJG=> SIN_TAYLOR
        
         PROFILING_STOP('init_bos_re_for_mp1')
      END SUBROUTINE INIT_BOS_RE_FOR_MP1

      ! *********************************************
      ! intitializes T>0 freq grid for imaginary part
      ! *********************************************
      SUBROUTINE INIT_BOS_IM( NOMEGA, QUAD ) 
         INTEGER                    :: NOMEGA
         TYPE( quadrature_handle )  :: QUAD 
         PROFILING_START('init_bos_im')
         ! increase distinct grid counter
         NGRIDS = NGRIDS + 1 
         ! flag this grid as a finite temperature grid 
         IMAG_GRIDS%LFINITE_TEMPERATURE = .TRUE.
         ! quadrature contains no zero point
         QUAD%N = NOMEGA 
         ! set identifier
         QUAD%ITYPE       = -3
         ! basis function used for quadrature error
         QUAD%PHI2        => VTANH2
         ! first derivative w.r.t. first argument of basis function
         QUAD%D_PHI2_DZ   => D_VTANH2_DZ
         ! first derivative w.r.t. second argument of basis function
         QUAD%D_PHI2_DL   => D_VTANH2_DL
         ! second derivative w.r.t. first argument of basis function
         QUAD%D2_PHI2_DZ2 => D2_VTANH2_DZ2
         ! object function that is being fitted 
         QUAD%F           => F_BETA_RRS
         ! first derivative of object function
         QUAD%D_F_DZ      => D1_F_BETA_RRS
         ! second derivative of object function
         QUAD%D2_F_DZ2    => D2_F_BETA_RRS
         ! set initial coefficients
         QUAD%CTAB        => FREQ_VTANH2
         ! calling function to set closest interval in table to requested one
!         QUAD%R1E3        => SET_VTANH2_R1E3
!         QUAD%R1E4        => SET_VTANH2_R1E4
!         QUAD%R1E5        => SET_VTANH2_R1E5
         QUAD%R1E6        => SET_VTANH2_R1E6
         ! quadrature scaling factor is temperature in eV (or inverse beta)
         QUAD%SCALING     = QDPD( SIGMA )
         ! grid id for parallelization
         QUAD%GRID_ID = NGRIDS-1

         ! set basis function to w/(x^2+w^2) tanh(x/2)
         QUAD%PHI         => VTANH
         ! dual basis function is sinh(x/2(1-2))/cosh(x/2) 
         QUAD%PSI         => VTIME
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR=> SIN_TAYLOR

         ! conjugate error
         ! set basis function to x/(x^2+w^2) tanh(x/2)
         QUAD%PHI_CONJG   => UTANH
         ! dual basis function is cosh(x/2(1-2))/cosh(x/2) 
         QUAD%PSI_CONJG   => UTIME
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR_CONJG=> COS_TAYLOR

         PROFILING_STOP('init_bos_im')
      END SUBROUTINE INIT_BOS_IM

      ! *********************************************
      ! intitializes T>0 freq grid for real part
      ! *********************************************
      SUBROUTINE INIT_FER_RE( NOMEGA, QUAD ) 
         INTEGER                    :: NOMEGA
         TYPE( quadrature_handle )  :: QUAD 
         PROFILING_START('init_fer_re')
         ! increase distinct grid counter
         NGRIDS = NGRIDS + 1 

         ! flag this grid as a finite temperature grid 
         IMAG_GRIDS%LFINITE_TEMPERATURE = .TRUE.
         ! quadrature contains no zero point
         QUAD%N = NOMEGA 
         ! set identifier
         QUAD%ITYPE       = -4
         ! basis function used for quadrature error
         QUAD%PHI2        => UFREQ2
         ! first derivative w.r.t. first argument of basis function
         QUAD%D_PHI2_DZ   => D_UFREQ2_DZ
         ! first derivative w.r.t. second argument of basis function
         QUAD%D_PHI2_DL   => D_UFREQ2_DL
         ! second derivative w.r.t. first argument of basis function
         QUAD%D2_PHI2_DZ2 => D2_UFREQ2_DZ2
         ! object function that is being fitted 
         QUAD%F           => F_BETA_RRS
         ! first derivative of object function
         QUAD%D_F_DZ      => D1_F_BETA_RRS
         ! second derivative of object function
         QUAD%D2_F_DZ2    => D2_F_BETA_RRS
         ! set initial coefficients
         QUAD%CTAB        => FREQ_U2
         ! calling function to set closest interval in table to requested one
         ! QUAD%R1E3        => SET_U2_R1E3
         ! QUAD%R1E4        => SET_U2_R1E4
         ! QUAD%R1E5        => SET_U2_R1E5
         QUAD%R1E6        => SET_U2_R1E6
         ! quadrature scaling factor is temperature in eV (or inverse beta)
         QUAD%SCALING     = QDPD( SIGMA )
         ! grid id for parallelization
         QUAD%GRID_ID = NGRIDS-1

         ! set basis function to x/(x^2+w^2)
         QUAD%PHI         => UFREQ
         ! dual basis function is sinh(x/2(1-2))/cosh(x/2) 
         QUAD%PSI         => VTIME
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR=> COS_TAYLOR

         ! set basis function to w/(x^2+w^2)
         QUAD%PHI_CONJG   => VFREQ
         ! dual basis function is cosh(x/2(1-2))/cosh(x/2) 
         QUAD%PSI_CONJG   => UTIME
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR_CONJG=> SIN_TAYLOR
         
         PROFILING_STOP('init_fer_re')
      END SUBROUTINE INIT_FER_RE

      ! *********************************************
      ! intitializes T>0 freq grid for real part
      ! *********************************************
      SUBROUTINE INIT_FER_RE_FOR_GM( NOMEGA, QUAD ) 
         USE minimax_dependence, ONLY:FIT_COEFF_COMP_UFREQ
         INTEGER                    :: NOMEGA
         TYPE( quadrature_handle )  :: QUAD 
         INTEGER                    :: R1(24),R2(24)
         DATA R1(6:24)/& 
         30,  & ! 6
         40,  & ! 7
         40,  & ! 8
         50,  & ! 9
         60,  & ! 10
         60,  & ! 11
         65,  & ! 12
         95,  & ! 13
         100, & ! 14
         105, & ! 15
         175, & ! 16
         190, & ! 17
         190, & ! 18
         230, & ! 19
         400, & ! 20
         400, & ! 21
         490, & ! 22
         590, & ! 23
         690/   ! 24
         DATA R2(6:24)/& 
         4E3,  &!  6
         4E3,  &!  7
         4E3,  &!  8
         4E3,  &!  9
         4E3,  &! 10
         4E3,  &! 11
         4E3,  &! 12
         3E4,  &! 13
         3E4,  &! 14
         3E4,  &! 15
         4E5,  &! 16
         4E5,  &! 17
         4E5,  &! 18
         5E5,  &! 19
         6E6,  &! 20
         6E6,  &! 21
         7E6,  &! 22
         7E6,  &! 23
         8E6/   ! 24
         PROFILING_START('init_fer_re_for_gm')
         ! increase distinct grid counter
         NGRIDS = NGRIDS + 1 

         ! flag this grid as a finite temperature grid 
         IMAG_GRIDS%LFINITE_TEMPERATURE = .TRUE.
         ! quadrature contains no zero point
         QUAD%N = NOMEGA 
         ! set identifier
         QUAD%ITYPE       = -40
         ! basis function used for quadrature error
         QUAD%PHI2        => UFREQ
         ! first derivative w.r.t. first argument of basis function
         QUAD%D_PHI2_DZ   => D_UFREQ_DZ
         ! first derivative w.r.t. second argument of basis function
         QUAD%D_PHI2_DL   => D_UFREQ_DL
         ! second derivative w.r.t. first argument of basis function
         QUAD%D2_PHI2_DZ2 => D2_UFREQ_DZ2
         ! object function that is being fitted 
         QUAD%F           => F_TANH
         ! first derivative of object function
         QUAD%D_F_DZ      => D1_F_TANH
         ! second derivative of object function
         QUAD%D2_F_DZ2    => D2_F_TANH
         ! set initial coefficients
         QUAD%CTAB        => FREQ_COMP_TANH
         ! calling function to set closest interval in table to requested one
         QUAD%R1E6        => SET_COMP_TANH_R1E6
         QUAD%R1E5        => SET_COMP_TANH_R1E4
         QUAD%R1E4        => SET_COMP_TANH_R1E4
         ! quadrature scaling factor is temperature in eV (or inverse beta)
         QUAD%SCALING     = QDPD( SIGMA )
         ! grid id for parallelization
         QUAD%GRID_ID = NGRIDS-1

         ! set basis function to x/(x^2+w^2)
         QUAD%PHI         => UFREQ
         ! dual basis function is sinh(x/2(1-2))/cosh(x/2) 
         QUAD%PSI         => VTIME
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR=> COS_TAYLOR

         ! set basis function to w/(x^2+w^2)
         QUAD%PHI_CONJG   => VFREQ
         ! dual basis function is cosh(x/2(1-2))/cosh(x/2) 
         QUAD%PSI_CONJG   => UTIME
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR_CONJG=> SIN_TAYLOR

         ! no fitting functions below order 6
         IF ( NOMEGA > 3 ) THEN
            ! associate function to calculate coefficients directly 
            QUAD%FITCOF => FIT_COEFF_COMP_UFREQ
         ENDIF
         PROFILING_STOP('init_fer_re_for_gm')
      END SUBROUTINE INIT_FER_RE_FOR_GM

      ! *********************************************
      ! intitializes T>0 freq grid for imaginary part
      ! *********************************************
      SUBROUTINE INIT_FER_IM( NOMEGA, QUAD ) 
         INTEGER                    :: NOMEGA
         TYPE( quadrature_handle )  :: QUAD 
         PROFILING_START('init_fer_im')
         ! increase distinct grid counter
         NGRIDS = NGRIDS + 1 
         ! flag this grid as a finite temperature grid 
         IMAG_GRIDS%LFINITE_TEMPERATURE = .TRUE.
         ! quadrature contains no zero point
         QUAD%N = NOMEGA 
         ! set identifier
         QUAD%ITYPE       = -5
         ! basis function is w/(x^2+w^2)
         QUAD%PHI2        => VFREQ2
         ! first derivative w.r.t. first argument of basis function
         QUAD%D_PHI2_DZ   => D_VFREQ2_DZ
         ! first derivative w.r.t. second argument of basis function
         QUAD%D_PHI2_DL   => D_VFREQ2_DL
         ! second derivative w.r.t. first argument of basis function
         QUAD%D2_PHI2_DZ2 => D2_VFREQ2_DZ2
         ! object function that is being fitted 
         QUAD%F           => F_BETA_LRS
         ! first derivative of object function
         QUAD%D_F_DZ      => D1_F_BETA_LRS
         ! second derivative of object function
         QUAD%D2_F_DZ2    => D2_F_BETA_LRS
         ! set initial coefficients
         QUAD%CTAB        => FREQ_V2
         ! calling function to set closest interval in table to requested one
         ! QUAD%R1E3        => SET_V2_R1E3
         ! QUAD%R1E4        => SET_V2_R1E4
         ! QUAD%R1E5        => SET_V2_R1E5
         QUAD%R1E6        => SET_V2_R1E6
         ! quadrature scaling factor is temperature in eV (or inverse beta)
         QUAD%SCALING     = QDPD( SIGMA )
         ! grid id for parallelization
         QUAD%GRID_ID = NGRIDS-1

         ! set basis function to w/(x^2+w^2)
         QUAD%PHI         => VFREQ
         ! dual basis function is cosh(x/2(1-2))/cosh(x/2) 
         QUAD%PSI         => UTIME
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR=> SIN_TAYLOR

         ! set basis function to x/(x^2+w^2)
         QUAD%PHI_CONJG         => UFREQ
         ! dual basis function is sinh(x/2(1-2))/cosh(x/2) 
         QUAD%PSI_CONJG         => VTIME
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR_CONJG=> COS_TAYLOR

         PROFILING_STOP('init_fer_im')
      END SUBROUTINE INIT_FER_IM

      ! *********************************************
      ! intitializes T>0 freq grid for im part
      ! *********************************************
      SUBROUTINE INIT_FER_IM_FOR_GM( NOMEGA, QUAD ) 
         USE minimax_dependence, ONLY:FIT_COEFF_COMP_UFREQ
         INTEGER                    :: NOMEGA
         TYPE( quadrature_handle )  :: QUAD 
         INTEGER                    :: R1(24),R2(24)
         DATA R1(6:24)/& 
         30,  & ! 6
         40,  & ! 7
         40,  & ! 8
         50,  & ! 9
         60,  & ! 10
         60,  & ! 11
         65,  & ! 12
         95,  & ! 13
         100, & ! 14
         105, & ! 15
         175, & ! 16
         190, & ! 17
         190, & ! 18
         230, & ! 19
         400, & ! 20
         400, & ! 21
         490, & ! 22
         590, & ! 23
         690/   ! 24
         DATA R2(6:24)/& 
         4E3,  &!  6
         4E3,  &!  7
         4E3,  &!  8
         4E3,  &!  9
         4E3,  &! 10
         4E3,  &! 11
         4E3,  &! 12
         3E4,  &! 13
         3E4,  &! 14
         3E4,  &! 15
         4E5,  &! 16
         4E5,  &! 17
         4E5,  &! 18
         5E5,  &! 19
         6E6,  &! 20
         6E6,  &! 21
         7E6,  &! 22
         7E6,  &! 23
         8E6/   ! 24
         PROFILING_START('init_fer_im_for_gm')
         ! increase distinct grid counter
         NGRIDS = NGRIDS + 1 

         ! flag this grid as a finite temperature grid 
         IMAG_GRIDS%LFINITE_TEMPERATURE = .TRUE.
         ! quadrature contains no zero point
         QUAD%N = NOMEGA 
         ! set identifier
         QUAD%ITYPE       = -40
         ! basis function used for quadrature error
         QUAD%PHI2        => UFREQ
         ! first derivative w.r.t. first argument of basis function
         QUAD%D_PHI2_DZ   => D_UFREQ_DZ
         ! first derivative w.r.t. second argument of basis function
         QUAD%D_PHI2_DL   => D_UFREQ_DL
         ! second derivative w.r.t. first argument of basis function
         QUAD%D2_PHI2_DZ2 => D2_UFREQ_DZ2
         ! object function that is being fitted 
         QUAD%F           => F_TANH
         ! first derivative of object function
         QUAD%D_F_DZ      => D1_F_TANH
         ! second derivative of object function
         QUAD%D2_F_DZ2    => D2_F_TANH
         ! set initial coefficients
         QUAD%CTAB        => FREQ_COMP_TANH
         ! calling function to set closest interval in table to requested one
         QUAD%R1E6        => SET_COMP_TANH_R1E6
         ! quadrature scaling factor is temperature in eV (or inverse beta)
         QUAD%SCALING     = QDPD( SIGMA )
         ! grid id for parallelization
         QUAD%GRID_ID = NGRIDS-1

         ! set basis function to w/(x^2+w^2)
         QUAD%PHI         => VFREQ
         ! dual basis function is cosh(x/2(1-2))/cosh(x/2) 
         QUAD%PSI         => UTIME
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR=> SIN_TAYLOR

         ! set basis function to x/(x^2+w^2)
         QUAD%PHI_CONJG   => UFREQ
         ! dual basis function is sinh(x/2(1-2))/cosh(x/2) 
         QUAD%PSI_CONJG   => VTIME
         ! Taylor expansion for low frequencies that relate dual and basis 
         QUAD%TRANS_TAYLOR_CONJG=> COS_TAYLOR

         ! no fitting functions below order 6
         IF ( NOMEGA > 3 ) THEN
            ! associate function to calculate coefficients directly 
            QUAD%FITCOF => FIT_COEFF_COMP_UFREQ
         ENDIF
         
         PROFILING_STOP('init_fer_im_for_gm')
      END SUBROUTINE INIT_FER_IM_FOR_GM

      ! *********************************************
      ! split communicators into 4 groups
      ! *********************************************
       SUBROUTINE SPLIT_COMMUNICATOR( COMM, NGRIDS, GRIDS )
          TYPE( communic )         :: COMM
          INTEGER                  :: NGRIDS
          TYPE( imag_grid_handle ) :: GRIDS
          ! local 
          INTEGER                  :: I 
          INTEGER                  :: NGROUPS
          TYPE( communic )         :: CINTRA
         
          PROFILING_START( 'split_communicator' )
          ! number of groups is 
          NGROUPS = NGRIDS - 1 
          ! global communicator
          GRIDS%COMM = COMM 
          ! time grid is never distributed
          GRIDS%TIME%COMM = COMM
          ! initialize internal quadrature communicator to global one
          CINTRA = COMM
          ! BOS_RE grid is done by first group
          GRIDS%FREQ_BOS_RE%COMM = CINTRA
          ! BOS_IM grid is done by second group
          GRIDS%FREQ_BOS_IM%COMM = CINTRA
          ! FER_RE grid is done by third group
          GRIDS%FREQ_FER_RE%COMM = CINTRA
          ! FER_IM grid is done by third group
          GRIDS%FREQ_FER_IM%COMM = CINTRA
#ifdef MPI
          ! GRIDS%COMM_INTER%NCPU gives the number of groups
          GRIDS%COMM_INTER%NCPU = 1 
#endif 

          ! return here if only 1 grid calculated
          IF ( NGROUPS < 2 ) THEN
             PROFILING_STOP( 'split_communicator' )
             RETURN
          ENDIF

#if ( defined MPI && ! defined noQuadPrecision )
          ! in case only one grid is computed -> no distribution is done
          IF ( NGROUPS-1 < COMM%NCPU ) THEN
             ! split global communicator into groups based on MAP
             ! note that INTER communicator is not a proper communicator
             CALLMPI( M_divide_into_groups( COMM, GRIDS%COMM_INTER, CINTRA, NGROUPS) )
          ENDIF 
#endif
          ! BOS_RE grid is done by first group
          GRIDS%FREQ_BOS_RE%COMM = CINTRA
          ! BOS_IM grid is done by second group
          GRIDS%FREQ_BOS_IM%COMM = CINTRA
          ! FER_RE grid is done by third group
          GRIDS%FREQ_FER_RE%COMM = CINTRA
          ! FER_IM grid is done by third group
          GRIDS%FREQ_FER_IM%COMM = CINTRA

          PROFILING_STOP( 'split_communicator' )
       END SUBROUTINE SPLIT_COMMUNICATOR

      ! *********************************************
      ! intitializes all T>0 hypergeometric grid 
      ! of Ozaki
      ! *********************************************
      SUBROUTINE INIT_HYPERGEOM( NTAU, GRIDS) 
         INTEGER  :: NTAU
         TYPE( imag_grid_handle )  :: GRIDS
         PROFILING_START('init_hypergeom')

         ! set pointer to correct basis functions 
         CALL INIT_T_TIME2( NTAU, IMAG_GRIDS%TIME )        
         ! only two grids are allocated, because there is no time transform
         CALL INIT_BOS_RE( NTAU, IMAG_GRIDS%FREQ_BOS_RE )        
         CALL INIT_BOS_IM( NTAU, IMAG_GRIDS%FREQ_BOS_IM )        
         CALL INIT_FER_RE( NTAU, IMAG_GRIDS%FREQ_FER_RE )        
         CALL INIT_FER_IM( NTAU, IMAG_GRIDS%FREQ_FER_IM )        
         ! set number of grids 
         NGRIDS = 1 
         ! overwrite number of grid points 
!         IMAG_GRIDS%FREQ_BOS_IM%N = NTAU       
         IMAG_GRIDS%FREQ_BOS_RE%N = NTAU       
         ! shifted grid
         IMAG_GRIDS%TIME%ITYPE = -101
         ! even and odd bosonic functions have same grid 
         IMAG_GRIDS%FREQ_BOS_RE%ITYPE = -202
         IMAG_GRIDS%FREQ_BOS_IM%ITYPE = -202
         ! even and odd fermionic functions have same grid 
         IMAG_GRIDS%FREQ_FER_RE%ITYPE = -203
         IMAG_GRIDS%FREQ_FER_IM%ITYPE = -203

         ! set initial coefficients to null pointer
         IMAG_GRIDS%TIME%CTAB => NULL()
         IMAG_GRIDS%FREQ_BOS_RE%CTAB => NULL()
         IMAG_GRIDS%FREQ_BOS_IM%CTAB => NULL()
         IMAG_GRIDS%FREQ_FER_RE%CTAB => NULL()
         IMAG_GRIDS%FREQ_FER_IM%CTAB => NULL()
         ! calling function to set closest interval in table to requested one
         IMAG_GRIDS%TIME%R1E3 => NULL() 
         IMAG_GRIDS%TIME%R1E4 => NULL() 
         IMAG_GRIDS%TIME%R1E5 => NULL() 
         IMAG_GRIDS%TIME%R1E6 => NULL() 

         IMAG_GRIDS%FREQ_BOS_RE%R1E3 => NULL() 
         IMAG_GRIDS%FREQ_BOS_RE%R1E4 => NULL() 
         IMAG_GRIDS%FREQ_BOS_RE%R1E5 => NULL() 
         IMAG_GRIDS%FREQ_BOS_RE%R1E6 => NULL() 

         IMAG_GRIDS%FREQ_BOS_IM%R1E3 => NULL() 
         IMAG_GRIDS%FREQ_BOS_IM%R1E4 => NULL() 
         IMAG_GRIDS%FREQ_BOS_IM%R1E5 => NULL() 
         IMAG_GRIDS%FREQ_BOS_IM%R1E6 => NULL() 

         IMAG_GRIDS%FREQ_FER_RE%R1E3 => NULL() 
         IMAG_GRIDS%FREQ_FER_RE%R1E4 => NULL() 
         IMAG_GRIDS%FREQ_FER_RE%R1E5 => NULL() 
         IMAG_GRIDS%FREQ_FER_RE%R1E6 => NULL() 

         IMAG_GRIDS%FREQ_FER_IM%R1E3 => NULL() 
         IMAG_GRIDS%FREQ_FER_IM%R1E4 => NULL() 
         IMAG_GRIDS%FREQ_FER_IM%R1E5 => NULL() 
         IMAG_GRIDS%FREQ_FER_IM%R1E6 => NULL() 

         IMAG_GRIDS%TIME%FITCOF => NULL()
         IMAG_GRIDS%FREQ_BOS_RE%FITCOF => NULL()
         IMAG_GRIDS%FREQ_BOS_IM%FITCOF => NULL()
         IMAG_GRIDS%FREQ_FER_RE%FITCOF => NULL()
         IMAG_GRIDS%FREQ_FER_IM%FITCOF => NULL()

         PROFILING_STOP('init_hypergeom')
      END SUBROUTINE INIT_HYPERGEOM

! Matsubara intializer
      ! *********************************************
      ! intitializes all T>0 Matsubara grid 
      ! *********************************************
      SUBROUTINE INIT_MATSUBARA( NTAU, GRIDS) 
         INTEGER  :: NTAU
         TYPE( imag_grid_handle )  :: GRIDS
         PROFILING_START('init_matsubara')

         ! set pointer to correct basis functions 
         CALL INIT_T_TIME2( NTAU, IMAG_GRIDS%TIME )        
         CALL INIT_BOS_RE( NTAU, IMAG_GRIDS%FREQ_BOS_RE )        
!         CALL INIT_BOS_IM( NTAU, IMAG_GRIDS%FREQ_BOS_IM )        
         CALL INIT_FER_RE( NTAU, IMAG_GRIDS%FREQ_FER_RE )        
         CALL INIT_FER_RE( NTAU, IMAG_GRIDS%FREQ_FER_IM )        
         ! set number of grids 
         NGRIDS = 1 
         ! overwrite number of grid points 
!         IMAG_GRIDS%FREQ_BOS_IM%N = NTAU       
         IMAG_GRIDS%FREQ_BOS_RE%N = NTAU       
         ! shifted grid
         IMAG_GRIDS%TIME%ITYPE = -101
         ! even and odd bosonic functions have same grid 
         IMAG_GRIDS%FREQ_BOS_IM%ITYPE = -102
         IMAG_GRIDS%FREQ_BOS_RE%ITYPE = -102
         ! even and odd fermionic functions have same grid 
         IMAG_GRIDS%FREQ_FER_IM%ITYPE = -103
         IMAG_GRIDS%FREQ_FER_RE%ITYPE = -103

         ! set initial coefficients to null pointer
         IMAG_GRIDS%TIME%CTAB => NULL()
         IMAG_GRIDS%FREQ_BOS_RE%CTAB => NULL()
         IMAG_GRIDS%FREQ_BOS_IM%CTAB => NULL()
         IMAG_GRIDS%FREQ_FER_RE%CTAB => NULL()
         IMAG_GRIDS%FREQ_FER_IM%CTAB => NULL()
         ! calling function to set closest interval in table to requested one
         IMAG_GRIDS%TIME%R1E3 => NULL() 
         IMAG_GRIDS%TIME%R1E4 => NULL() 
         IMAG_GRIDS%TIME%R1E5 => NULL() 
         IMAG_GRIDS%TIME%R1E6 => NULL() 

         IMAG_GRIDS%FREQ_BOS_RE%R1E3 => NULL() 
         IMAG_GRIDS%FREQ_BOS_RE%R1E4 => NULL() 
         IMAG_GRIDS%FREQ_BOS_RE%R1E5 => NULL() 
         IMAG_GRIDS%FREQ_BOS_RE%R1E6 => NULL() 

         IMAG_GRIDS%FREQ_BOS_IM%R1E3 => NULL() 
         IMAG_GRIDS%FREQ_BOS_IM%R1E4 => NULL() 
         IMAG_GRIDS%FREQ_BOS_IM%R1E5 => NULL() 
         IMAG_GRIDS%FREQ_BOS_IM%R1E6 => NULL() 

         IMAG_GRIDS%FREQ_FER_RE%R1E3 => NULL() 
         IMAG_GRIDS%FREQ_FER_RE%R1E4 => NULL() 
         IMAG_GRIDS%FREQ_FER_RE%R1E5 => NULL() 
         IMAG_GRIDS%FREQ_FER_RE%R1E6 => NULL() 

         IMAG_GRIDS%FREQ_FER_IM%R1E3 => NULL() 
         IMAG_GRIDS%FREQ_FER_IM%R1E4 => NULL() 
         IMAG_GRIDS%FREQ_FER_IM%R1E5 => NULL() 
         IMAG_GRIDS%FREQ_FER_IM%R1E6 => NULL() 

         IMAG_GRIDS%TIME%FITCOF => NULL()
         IMAG_GRIDS%FREQ_BOS_RE%FITCOF => NULL()
         IMAG_GRIDS%FREQ_BOS_IM%FITCOF => NULL()
         IMAG_GRIDS%FREQ_FER_RE%FITCOF => NULL()
         IMAG_GRIDS%FREQ_FER_IM%FITCOF => NULL()

         PROFILING_STOP('init_matsubara')
      END SUBROUTINE INIT_MATSUBARA 

   END SUBROUTINE INITIALIZE_QUADRATURES

!****************************************************************************
!> DESCRIPTION:
!> calculates minimax quadrature for requested (scaled) interval [X1,X2]
!> A is implicitly set to 1 or 0 
!****************************************************************************
!#define PrintQuadrature
!#define PrintQuadratureMinimax
#ifdef MinimaxTest
!#define PrintQuadratureMinimax
#endif
   SUBROUTINE DETERMINE_QUADRATURE( QUADRATURE, B_REQUESTED, IO, &
      COMM_INTER ) 
      TYPE( quadrature_handle )  :: QUADRATURE 
      REAL(q)                    :: B_REQUESTED
      TYPE( in_struct )          :: IO
      TYPE( communic ), OPTIONAL :: COMM_INTER
      ! local variables
      QDPREAL                    :: B
      INTEGER                    :: NODE_ME = 1
      INTEGER                    :: IERR

      PROFILING_START( 'DETERMINE_QUADRATURE' )

#ifdef MPI
      NODE_ME = QUADRATURE%COMM%NODE_ME
#endif

      ! allocate quadrature 
      CALL ALLOCATE_QUADRATURE_HANDLE( QUADRATURE, IO ) 

#if ( defined MPI && ! defined noQuadPrecision )
      ! early return in case of distributed calculations
      IF ( PRESENT( COMM_INTER ) ) THEN
        IF ( QUADRATURE%GRID_ID /= COMM_INTER%NODE_ME .AND. COMM_INTER%NCPU/=1) THEN
           PROFILING_STOP( 'DETERMINE_QUADRATURE' )
           RETURN
        ENDIF
      ENDIF
#endif

      B = QDPD(B_REQUESTED)
      ! in case the MM solution can be calculated from 
      ! a polynomial fit 
      CALL COMPUTE_SOLUTION_FROM_FIT( QUADRATURE, B, IERR, IO) 

      IF ( IERR > 0 ) GOTO 999 
      IF ( IERR < 0 ) GOTO 998 

      ! to guarantee convergence of solver, the boundary has to be adjusted
      ! this routine also sets non-optimized grids
      ! such as Matsubara and the hypbergeometric Ozaki grid
      CALL ADJUST_BOUNDARY_FOR_CONVERGENCE( QUADRATURE, B, IO) 
      
      ! as first start obtain requested fitting order 
      ! this step might increase the minimization interval
      ! since convergence is prioritized upon requested 
      ! minimization interval
      CALL INCREASE_GRID_ORDER( QUADRATURE, IO ) 

998   CONTINUE
      ! requested minimization interval is B_REQUESTED
      B = QDPD( B_REQUESTED )

       
      ! if requested minimization interval is smaller than 
      ! interval needs to be decreased 
      ! obtain alternant, this checks if fit is OK
      CALL DECREASE_INTERVAL( QUADRATURE, B, IO ) 

      B_REQUESTED = B 
      ! if requested minimization interval is smaller than 
      ! interval needs to be decreased 
      ! obtain alternant, this checks if fit is OK
      CALL OBTAIN_LSQ_SOLUTION( QUADRATURE, IO ) 

999   CONTINUE 

      ! Solve Remez algorithm 
      CALL REMEZ_ALGORITHM( QUADRATURE, IO )

#ifdef PrintQuadratureMinimax
IF ( QUADRATURE%N<32 ) THEN
      CALL PRINT_ERROR_FUNCTION( QUADRATURE, IO, IT = QUADRATURE%ITYPE, &
         FNAME="erf-MM" )
      CALL PRINT_QUADRATURE( QUADRATURE, IO, IT = QUADRATURE%ITYPE, &
         FNAME="quad-MM" )
ENDIF
#endif

      IF ( NODE_ME == 1 .AND. IO%NWRITE>=2 ) THEN
         IF ( QUADRATURE%ITYPE == 1 ) THEN
            WRITE(*,1)TOPRINT(QUADRATURE%INFINITYNORM*QUADRATURE%SCALING)
         ELSE IF ( QUADRATURE%ITYPE == 2 ) THEN
            WRITE(*,2)TOPRINT(QUADRATURE%INFINITYNORM/QUADRATURE%SCALING)
         ELSE IF ( QUADRATURE%ITYPE == 3 ) THEN
            WRITE(*,3)TOPRINT(QUADRATURE%INFINITYNORM/QUADRATURE%SCALING)
         ELSE IF ( QUADRATURE%ITYPE == 20 ) THEN
            WRITE(*,20)TOPRINT(QUADRATURE%INFINITYNORM/QUADRATURE%SCALING)
         ELSE IF ( QUADRATURE%ITYPE == -1 ) THEN
            WRITE(*,11)TOPRINT(QUADRATURE%INFINITYNORM*QUADRATURE%SCALING)
         ELSE IF ( QUADRATURE%ITYPE == -2 ) THEN
            WRITE(*,12)TOPRINT(QUADRATURE%INFINITYNORM/QUADRATURE%SCALING)
         ELSE IF ( QUADRATURE%ITYPE == -3 ) THEN
            WRITE(*,13)TOPRINT(QUADRATURE%INFINITYNORM/QUADRATURE%SCALING)
         ELSE IF ( QUADRATURE%ITYPE == -4 .OR. QUADRATURE%ITYPE == -40 ) THEN
            WRITE(*,14)TOPRINT(QUADRATURE%INFINITYNORM/QUADRATURE%SCALING)
         ELSE IF ( QUADRATURE%ITYPE == -5 ) THEN
            WRITE(*,15)TOPRINT(QUADRATURE%INFINITYNORM/QUADRATURE%SCALING)
         ! non-optimized grids
         ELSE IF ( QUADRATURE%ITYPE <- 100  ) THEN
         ELSE
            WRITE(*,100)QUADRATURE%ITYPE, &
              TOPRINT(QUADRATURE%INFINITYNORM/QUADRATURE%SCALING)
         ENDIF
      ENDIF

      PROFILING_STOP( 'DETERMINE_QUADRATURE' )
1     FORMAT( ' time grid (T=0) determined with error:',E11.4 )
2     FORMAT( ' bosonic grid (T=0, re) determined with error:',E11.4 )
3     FORMAT( ' fermionic grid (T=0, im) determined with error:',E11.4 )
20    FORMAT( ' fermionic grid (T=0, re) determined with error:',E11.4 )
11    FORMAT( ' time grid (T>0) determined with error:',E11.4 )
12    FORMAT( ' bosonic grid (T>0, re) determined with error:',E11.4 )
13    FORMAT( ' bosonic grid (T>0, im) determined with error:',E11.4 )
14    FORMAT( ' fermionic grid (T>0, re) determined with error:',E11.4 )
15    FORMAT( ' fermionic grid (T>0, im) determined with error:',E11.4 )
100   FORMAT( ' type ',I3,' grid determined with error ',E11.4 )

      CONTAINS 
      !**********************************************************************
      !> allocates quadrature handle
      !**********************************************************************
      SUBROUTINE ALLOCATE_QUADRATURE_HANDLE( QUADRATURE, IO ) 
         USE string, ONLY: str
         USE tutor, ONLY: vtutor
         TYPE( quadrature_handle )  :: QUADRATURE 
         TYPE( in_struct )          :: IO
         ! local 
         INTEGER                    :: I 

         PROFILING_START( 'allocate_quadrature_handle' )

         IF ( QUADRATURE%N > 0 ) THEN
            NULLIFY(QUADRATURE%C)
            NULLIFY(QUADRATURE%X0)
            IF(.NOT.ASSOCIATED(QUADRATURE%C )) ALLOCATE(QUADRATURE%C (2*QUADRATURE%N+1)) 
            IF(.NOT.ASSOCIATED(QUADRATURE%X0)) ALLOCATE(QUADRATURE%X0(2*QUADRATURE%N+1)) 
            DO I = 1, 2*QUADRATURE%N+1
               QUADRATURE%C(I) = QDPP(0)
               QUADRATURE%X0(I) = QDPP(0)
            ENDDO
            QUADRATURE%INFINITYNORM = QDPP(0)
            QUADRATURE%A = QDPP(0)
            QUADRATURE%B = QDPP(0)
         ELSE
            ! error, since quadrature not flagged for allocation
            CALL vtutor%error("ERROR, quadrature not flagged for allocation " // str(QUADRATURE%ITYPE))
         ENDIF

         PROFILING_STOP( 'allocate_quadrature_handle' )
      END SUBROUTINE ALLOCATE_QUADRATURE_HANDLE

      !**********************************************************************
      !> this subroutine computes first guess for minimization 
      !> it increases minimization interval to obtain the fit for requested 
      !> fitting order QUADRATURE%N. 
      !**********************************************************************
      SUBROUTINE INCREASE_GRID_ORDER( QUADRATURE, IO ) 
         USE string, ONLY: str
         USE tutor, ONLY: vtutor
         TYPE( quadrature_handle )  :: QUADRATURE 
         TYPE( in_struct )          :: IO
         ! local variables
         INTEGER, PARAMETER         :: NMAX=300
         INTEGER, PARAMETER         :: IGUESS_MAX = 10
         QDPREAL, ALLOCATABLE       :: Y(:)
         QDPREAL, ALLOCATABLE       :: T(:,:)
         QDPREAL, ALLOCATABLE       :: A(:,:)
         QDPREAL, ALLOCATABLE       :: W(:)
         QDPREAL, ALLOCATABLE       :: ALF(:)
         QDPREAL, ALLOCATABLE       :: BETA(:)
         QDPREAL, ALLOCATABLE       :: ALFINIT(:)
         QDPREAL, ALLOCATABLE       :: BETAINIT(:)
         INTEGER                    :: I, IERR,  J
         INTEGER, PARAMETER         :: IVERBOSITY=3
         INTEGER                    :: IGUESS
         INTEGER                    :: NSTART=1
         INTEGER                    :: N
         INTEGER                    :: IV =1
         INTEGER                    :: P  =1
         INTEGER                    :: L  =1
         INTEGER                    :: NL =1
         LOGICAL                    :: LRETRY = .FALSE.
         !parameter for LSQ search algorithm 
#ifdef qd_emulate 
#define RANDLSQ  qdreal('50.')
#else
         QDPREAL , PARAMETER :: RANDLSQ =  50._qd
#endif
         PROFILING_START( 'increase_grid_order' )

         ! early return in case of non-opitmized grids
         IF ( QUADRATURE%ITYPE < -100 ) THEN
            PROFILING_STOP( 'increase_grid_order' )
            RETURN
         ENDIF 

         IF ( IO%IU0>=0 .AND. IO%NWRITE>IVERBOSITY ) THEN
            WRITE(*,'(A,I3,A,E10.4)')"  Looking for a LSQ solution for ITYPE=", &
               QUADRATURE%ITYPE, " B=", TOPRINT(QUADRATURE%B)
         ENDIF
         ! start from botton
         NL=1             !# number of non-linear parameters
         L=NL             !# of linear parameters
         P=NL             !# of non-vanishing non-linear partial derivatives
         IV=1             !# of linear independent variables 
   
         ALLOCATE(Y(NMAX),T(NMAX,1),A(NMAX,2*AMAX+2),W(NMAX))
DOESI    Y = 0; T = 0; A = 0;W = 0
         ALLOCATE(ALF(AMAX), BETA(BMAX))
DOESI    ALF = 0; BETA = 0
         ALLOCATE(ALFINIT(AMAX), BETAINIT(BMAX))
DOESI    ALFINIT = 0; BETAINIT = 0
   
         !initialize saved coefficients 
         DO I = 1, AMAX
            ALFINIT(I) =0 
         ENDDO
         DO I = 1, BMAX
           BETAINIT(I)=0 
         ENDDO
         DO I = 1, NMAX 
            Y(I) = 0
            W(I) = 0
            DO J= 1, 2*AMAX+2 
               A(I,J) = 0
            ENDDO
         ENDDO
   
         CALL SETUP_SAMPLING_POINTS( QUADRATURE, T, Y, W, IV, IO)
   
         ! skip procedure below if quadrature points have been found already 
         IF( QUADRATURE%INFINITYNORM /= 0 ) THEN
            IF ( IO%IU0 >=0 .AND. IO%NWRITE > IVERBOSITY ) WRITE(*,*)&
            ' first guess fitting coefficients found in table'
            !optionally, dump coefficients to OUTCAR
            IF ( IO%NWRITE>IVERBOSITY .AND. IO%IU6>=0) THEN
               WRITE(IO%IU6,*)' first guess fitting coefficients found in table'
               WRITE(IO%IU6,154)QUADRATURE%N
               DO I = 1, QUADRATURE%N
                  WRITE(IO%IU6,156)TOPRINT(QUADRATURE%C(I)),TOPRINT(QUADRATURE%C(I+QUADRATURE%N))
               ENDDO
            ENDIF
#ifdef PrintQuadrature
            CALL PRINT_ERROR_FUNCTION( QUADRATURE, IO, IT = QUADRATURE%ITYPE, &
               FNAME="erf-ini" )
            CALL PRINT_QUADRATURE( QUADRATURE, IO, IT = QUADRATURE%ITYPE, &
               FNAME="quad-ini" )
#endif
            GOTO 10 
         ENDIF
   
!!increase fitting order from NSTART 
!!to requested order stored in QUADRATURE%N
         blowup: DO N =NSTART, QUADRATURE%N
            NL=N  !# of current non-linear parameters
            L=NL  !# of linear parameters
            P=NL  !# of non-vanishing non-linear partial derivatives
            !use as starting guess coefficients of previous order
            IF ( N>1)THEN
               ALF =0
               BETA=0
               DO I = 1, N-1
                  ALF( I ) = ALFINIT( I )
                  BETA( I ) = BETAINIT( I )
               ENDDO
            ENDIF
            ! we need initial guesses for the parameters 
            guess: DO IGUESS =1, IGUESS_MAX
               ALF( N ) =((IGUESS-1)*RANDLSQ)/(IGUESS_MAX-1)
               BETA( N )=((IGUESS-1)*RANDLSQ)/(IGUESS_MAX-1)
               !print initial guess of non-linear coefficients 
               IF ( IO%IU6>=0 .AND. IO%NWRITE>IVERBOSITY ) THEN
                  IF ( N<QUADRATURE%N ) WRITE (IO%IU6,130) (TOPRINT(ALF( I )), I = 1, N )
                  WRITE (IO%IU6,140)
               ENDIF
               
               CALL VARPRO(N,N,NDATALSQ,NMAX,L+P+2,IV,T,Y,W,A,IO%NWRITE-2,ALF,BETA,ITMAXLSQ,&
                  QUADRATURE,IERR,IO%IU6,IO%NWRITE)
   
               !varpro converged
               IF (IERR>0) THEN
                  IF (LRETRY) THEN
                     IF ( IO%NWRITE>IVERBOSITY .AND. IO%IU6>=0 ) WRITE(IO%IU6,151)IGUESS
                     LRETRY=.FALSE.
                  ENDIF 
                  IF ( IO%NWRITE>IVERBOSITY .AND. IO%IU6>=0 ) WRITE(IO%IU6,152)NL,IERR
                  IF ( IO%NWRITE>IVERBOSITY .AND. IO%IU6>=0 ) WRITE(*,152)NL,IERR
                  !converged -> save coefficients for next higher order
                  ALFINIT =0
                  BETAINIT=0
                  DO I = 1, N
                     ALFINIT( I ) = ALF( I )
                     BETAINIT( I ) = BETA( I )
                  ENDDO
                  EXIT guess
               !varpro not converged 
               ELSE
                  IF (IO%IU6>=0 .AND. IO%NWRITE>IVERBOSITY ) WRITE (IO%IU6,'("                IERR =",2I7)') &
                     IERR,IGUESS
                  IF ( .NOT. LRETRY ) THEN
                     IF (IO%IU6>=0 .AND. IO%NWRITE>IVERBOSITY ) WRITE(IO%IU6,153)N
                  ENDIF 
                  LRETRY=.TRUE.
               ENDIF  
            ENDDO guess
   
            !truly converged, if less iterations than IGUESS_MAX needed 
            IF ( IGUESS <= IGUESS_MAX) THEN
               !optionally, dump coefficients to OUTCAR
               IF ( (N==QUADRATURE%N .AND. IO%NWRITE>IVERBOSITY) .AND. IO%IU6>=0) THEN
                  WRITE(IO%IU6,154)N
                  WRITE(IO%IU6,156)(TOPRINT(ALF(I)),TOPRINT(BETA(I)),I=1,N)
               ENDIF
   
               !save coefficients 
               DO I=1,N
                  QUADRATURE%C(I)=ABS(ALF(I))
                  QUADRATURE%C(I+N)=BETA(I)
               ENDDO
               ! update norm
               QUADRATURE%C(2*QUADRATURE%N+1) = ERROR_FUNCTION( QUADRATURE, QUADRATURE%B ) 
               QUADRATURE%INFINITYNORM = ABS( QUADRATURE%C(2*QUADRATURE%N+1) ) 
            ELSE
               CALL vtutor%error("INCREASE_GRID_ORDER failed to find a LSQ solutions for " // &
                  str(QUADRATURE%B))
            ENDIF
         ENDDO blowup
   
10       CONTINUE 

         DEALLOCATE(Y,T,A,W)
         DEALLOCATE(ALF, BETA)
         DEALLOCATE(ALFINIT, BETAINIT)
            
130      FORMAT ('0 INITIAL NONLINEAR PARAMETERS'/(4E20.7))
140      FORMAT ('0',50('*'))
151      FORMAT(/,'      success after ',I4,' attempts.')
152      FORMAT('     Order',I3' done after ',I6,' iterations')
153      FORMAT('          using random numbers for order',I3,' .')
154      FORMAT(' LSQ quadrature of ',I2,'th order found, abszissas and weights:')
156      FORMAT(2F45.32)
         PROFILING_STOP( 'increase_grid_order' )
   
      END SUBROUTINE INCREASE_GRID_ORDER

      !**********************************************************************
      !> decreases minimization interval from QUADRATURE%B to B_REQUESTED
      !> using VARPRO, if machine precission is reached minimization is 
      !> aborted.
      !**********************************************************************
!#define PrintQuadratureDecrease
      SUBROUTINE DECREASE_INTERVAL( QUADRATURE, B_REQUESTED, IO ) 
         TYPE( quadrature_handle )  :: QUADRATURE 
         QDPREAL                    :: B_REQUESTED
         TYPE( in_struct )          :: IO
         ! local variables
         INTEGER, PARAMETER         :: NMAX=300
         INTEGER, PARAMETER         :: IDOWN_MAX = 5000 
         INTEGER, PARAMETER         :: IVERBOSITY= 3
         QDPREAL                    :: Y(NMAX)
         QDPREAL                    :: T(NMAX,1)
         QDPREAL                    :: A(NMAX,2*AMAX+2)
         QDPREAL                    :: W(NMAX)
         QDPREAL                    :: ALF(AMAX)
         QDPREAL                    :: BETA(BMAX)
         INTEGER                    :: IERR =0
         INTEGER                    :: I,  J
         INTEGER                    :: IV =1
         INTEGER                    :: P  =1
         INTEGER                    :: L  =1
         INTEGER                    :: NL =1
         INTEGER                    :: IDOWN
         QDPREAL                    :: MINMAX(2)
         QDPREAL                    :: YZ 
         QDPREAL                    :: SSCALE

         PROFILING_START( 'decrease_interval' )

         ! early return in case of non-opitmized grids
         IF ( QUADRATURE%ITYPE < -100 ) THEN
            PROFILING_STOP( 'decrease_interval' )
            RETURN
         ENDIF 

         ! early exit 
         IF ( QUADRATURE%B < B_REQUESTED ) THEN
            PROFILING_STOP( 'decrease_interval' )
            RETURN
         ENDIF 
         IF ( IO%IU0>=0 .AND. IO%NWRITE>IVERBOSITY) THEN
            WRITE(*,'(A,2E12.4)')"  Decreasing interval ", &
               TOPRINT(QUADRATURE%B), TOPRINT(B_REQUESTED)
         ENDIF

         ! start from botton
         NL=QUADRATURE%N  !# number of non-linear parameters
         L=NL             !# of linear parameters
         P=NL             !# of non-vanishing non-linear partial derivatives
         IV=1             !# of linear independent variables 
   
         !initialize saved coefficients 
         DO I = 1, AMAX
            ALF(I) = 0 
            IF ( I <= QUADRATURE%N ) THEN
               ALF( I ) = QUADRATURE%C( I ) 
            ENDIF
         ENDDO
         DO I = 1, BMAX
            BETA(I) = 0 
            IF ( I <= QUADRATURE%N ) THEN
               BETA( I ) = QUADRATURE%C( I+QUADRATURE%N ) 
            ENDIF
         ENDDO
         DO I = 1, NMAX 
            Y(I) = 0
            W(I) = 0
            DO J= 1, 2*AMAX+2 
               A(I,J) = 0
            ENDDO
         ENDDO
        
         CALL ADAPT_CONVERGENCE_PARAMETER( QUADRATURE, ITMAXLSQ, SSCALE ) 

#ifdef PrintQuadratureDecrease
         ! this will use REMEZ algorithm from begining 
         IERR =1000
         CALL PRINT_QUADRATURE( QUADRATURE, IO, &
            IT=ABS( QUADRATURE%ITYPE )*100, FNAME="quad-MM_R" )
!         CALL PRINT_ERROR_FUNCTION( QUADRATURE, IO, &
!            IT = ABS( QUADRATURE%ITYPE )*100, FNAME="erf-MM_R")
#endif

         downsizing: DO IDOWN = 1, IDOWN_MAX
         IF ( IDOWN == 1 .AND. IO%NWRITE>IVERBOSITY .AND. IO%IU0 >=0 ) WRITE(*,1)
         IF ( B_REQUESTED < QUADRATURE%B ) THEN
            ! set start coefficients
            DO I = 1, QUADRATURE%N
               ALF( I ) = QUADRATURE%C( I ) 
               BETA( I ) = QUADRATURE%C( I+QUADRATURE%N ) 
            ENDDO
            CALL DETERMINE_ALTERNANT( QUADRATURE, IO )
            ! obtain minmax values
            MINMAX(1) = QDPP( 1E6 )
            MINMAX(2) =-QDPP( 1E6 )
            DO I = 1, 2*QUADRATURE%N+1
               YZ = ERROR_FUNCTION( QUADRATURE, QUADRATURE%X0( I ) )
               MINMAX(1) = MIN( MINMAX(1) , ABS( YZ ) )
               MINMAX(2) = MAX( MINMAX(2) , ABS( YZ ) )
            ENDDO
            QUADRATURE%C(2*QUADRATURE%N+1) = MINMAX(2)

            ! break condition:
            IF ( MINMAX(2) < QDPP(1E-24) ) THEN
               IF (IO%IU0>=0 .AND. IO%NWRITE>IVERBOSITY ) WRITE(*,*)&
               ' machine precision reached, aborting interval contraction'
               B_REQUESTED = QUADRATURE%B
               EXIT downsizing
            ENDIF

            CALL SET_RNEW( QUADRATURE, MINMAX, SSCALE, IERR, IVERBOSITY, IO )

            ! set sampling coefficients based on alternant         
            IF( QUADRATURE%ITYPE > 0 ) THEN
               CALL SETUP_SAMPLING_POINTS( QUADRATURE, T, Y, W, IV, IO )
            ELSE
               IF( IERR > 200 ) THEN
                  NDATALSQ=2*QUADRATURE%N+1
               ELSE
                  NDATALSQ=100
               ENDIF
               CALL SAMPLING_FROM_ALTERNANT( QUADRATURE, T, Y, W, IV, IO )
            ENDIF
   
            IF( IERR > 200 ) THEN
               CALL REMEZ_ALGORITHM( QUADRATURE, IO )
               SSCALE = QDPP(1)
            ELSE
               CALL VARPRO(QUADRATURE%N,QUADRATURE%N,NDATALSQ,NMAX,&
                  L+P+2,IV,T,Y,W,A,IO%NWRITE-2,ALF,BETA,ITMAXLSQ,&
                  QUADRATURE,IERR,IO%IU6,IO%NWRITE)
            ENDIF
      
            !varpro converged
            IF (IERR>0) THEN
               !converged -> save coefficients for next higher order
               DO I = 1, QUADRATURE%N
                  QUADRATURE%C( I ) = ALF( I )
                  QUADRATURE%C( I+QUADRATURE%N ) = BETA( I )
               ENDDO
               CALL ADAPT_SCALING( IERR, SSCALE )
               IF ( IO%NWRITE>IVERBOSITY ) THEN
                  IF( IO%IU0 >=0 )  WRITE(*,'( I8,F10.3)')IERR,TOPRINT(SSCALE)
                  IF( IO%IU6 >=0 )  WRITE(IO%IU6,'( I8,F10.3)')IERR,TOPRINT(SSCALE)
               ENDIF

            IF( IERR > 200 ) THEN
               CALL REMEZ_ALGORITHM( QUADRATURE, IO )
            ENDIF
#ifdef PrintQuadratureDecrease
            CALL PRINT_QUADRATURE( QUADRATURE, IO, &
               IT=ABS( QUADRATURE%ITYPE )*100+IDOWN, FNAME="quad-MM_R" )
!            CALL PRINT_ERROR_FUNCTION( QUADRATURE, IO, &
!               IT = ABS( QUADRATURE%ITYPE )*100+IDOWN , FNAME="erf-MM_R")
#endif
            ELSE
               WRITE(*,*)' '
               WRITE(*,*) 'ERROR, DECREASE_INTERVAL not converged',IERR, TOPRINT(QUADRATURE%B*SSCALE)
            ENDIF
         ELSE
            IF( IO%NWRITE>IVERBOSITY ) THEN
            IF( IO%IU0>=0  ) &
               WRITE(*,'(A,2E14.7)') '  Requested interval reached: ',&
                  TOPRINT(QUADRATURE%B), TOPRINT(B_REQUESTED)
            IF( IO%IU6>=0  ) &
               WRITE(IO%IU6,'(A,2E14.7)') '  Requested interval reached: ',&
                  TOPRINT(QUADRATURE%B), TOPRINT(B_REQUESTED)
            ENDIF
            QUADRATURE%B =  MAX( B_REQUESTED , QUADRATURE%B ) 
            EXIT downsizing
         ENDIF
         ENDDO downsizing
1        FORMAT ('        R_old       R_new       minimum     maximum      iters   damping' )
         PROFILING_STOP( 'decrease_interval' )
      END SUBROUTINE DECREASE_INTERVAL

      !**********************************************************************
      ! determines LSQ solution 
      !**********************************************************************
      SUBROUTINE OBTAIN_LSQ_SOLUTION( QUADRATURE, IO ) 
         TYPE( quadrature_handle )  :: QUADRATURE 
         TYPE( in_struct )          :: IO
         ! local variables
         INTEGER, PARAMETER         :: NMAX=300
         INTEGER, PARAMETER         :: IDOWN_MAX = 300 
         INTEGER, PARAMETER         :: IVERBOSITY= 3
         QDPREAL                    :: Y(NMAX)
         QDPREAL                    :: T(NMAX,1)
         QDPREAL                    :: A(NMAX,2*AMAX+2)
         QDPREAL                    :: W(NMAX)
         QDPREAL                    :: ALF(AMAX)
         QDPREAL                    :: BETA(BMAX)
         INTEGER                    :: I, IERR,  J
         INTEGER                    :: IV =1
         INTEGER                    :: P  =1
         INTEGER                    :: L  =1
         INTEGER                    :: NL =1
         QDPREAL                    :: MINMAX(2)
         QDPREAL                    :: YZ 
         INTEGER                    :: NODE_ME = 1
         PROFILING_START( 'obtain_lsq_solution' )
#ifdef MPI
         NODE_ME = QUADRATURE%COMM%NODE_ME
#endif
         ! early return in case of non-opitmized grids
         IF ( QUADRATURE%ITYPE < -100 ) THEN
            PROFILING_STOP( 'obtain_lsq_solution' )
            RETURN
         ENDIF 

         IF ( IO%IU0>=0 .AND. IO%NWRITE>IVERBOSITY ) THEN
            WRITE(*,'(A,E12.4)')"  Determining LSQ solution for B= ", TOPRINT(QUADRATURE%B)
         ENDIF

         ! start from botton
         NL=QUADRATURE%N  !# number of non-linear parameters
         L=NL             !# of linear parameters
         P=NL             !# of non-vanishing non-linear partial derivatives
         IV=1             !# of linear independent variables 
   
         !initialize saved coefficients 
         DO I = 1, AMAX
            ALF(I) = 0 
            IF ( I <= QUADRATURE%N ) THEN
               ALF( I ) = QUADRATURE%C( I ) 
            ENDIF
         ENDDO
         DO I = 1, BMAX
            BETA(I) = 0 
            IF ( I <= QUADRATURE%N ) THEN
               BETA( I ) = QUADRATURE%C( I+QUADRATURE%N ) 
            ENDIF
         ENDDO
         DO I = 1, NMAX 
            Y(I) = 0
            W(I) = 0
            DO J= 1, 2*AMAX+2 
               A(I,J) = 0
            ENDDO
         ENDDO
        
         CALL ADAPT_CONVERGENCE_PARAMETER( QUADRATURE, ITMAXLSQ ) 

         CALL DETERMINE_ALTERNANT( QUADRATURE, IO )
         ! obtain minmax values
         MINMAX(1) = QDPP( 1E6 )
         MINMAX(2) =- QDPP( 1E6 )
         DO I = 1, 2*QUADRATURE%N+1
            YZ = ERROR_FUNCTION( QUADRATURE, QUADRATURE%X0( I ) )
            MINMAX(1) = MIN( MINMAX(1) , ABS( YZ ) )
            MINMAX(2) = MAX( MINMAX(2) , ABS( YZ ) )
         ENDDO
         QUADRATURE%C(2*QUADRATURE%N+1) = MINMAX(2)
         ! set sampling coefficients based on alternant         
         IF( QUADRATURE%ITYPE > 0 ) THEN
            CALL SETUP_SAMPLING_POINTS( QUADRATURE, T, Y, W, IV, IO )
         ELSE
            CALL SAMPLING_FROM_ALTERNANT( QUADRATURE, T, Y, W, IV, IO )
         ENDIF
   
         CALL VARPRO(QUADRATURE%N,QUADRATURE%N,NDATALSQ,NMAX,&
            L+P+2,IV,T,Y,W,A,IO%NWRITE-2,ALF,BETA,ITMAXLSQ,&
            QUADRATURE,IERR,IO%IU6,IO%NWRITE)
      
         !varpro converged
         IF (IERR>0) THEN
            !converged -> save coefficients for next higher order
            DO I = 1, QUADRATURE%N
               QUADRATURE%C( I ) = ALF( I )
               QUADRATURE%C( I+QUADRATURE%N ) = BETA( I )
               IF ( IO%IU6>=0 .AND. IO%NWRITE>IVERBOSITY ) &
               WRITE(IO%IU6,'(2F25.16)') TOPRINT( ALF( I) ), TOPRINT( BETA( I ) )
            ENDDO
         ELSE
            IF ( NODE_ME == 1 ) &
               WRITE(*,*) 'WARNING, OBTAIN_LSQ_SOLUTION not converged',QUADRATURE%B
         ENDIF
         PROFILING_STOP( 'obtain_lsq_solution' )
      END SUBROUTINE OBTAIN_LSQ_SOLUTION

      !**********************************************************************
      ! might increase total number of allowed iterations for VARPRO 
      ! depending on the requested grid order 
      !**********************************************************************
      SUBROUTINE ADAPT_CONVERGENCE_PARAMETER( QUADRATURE, ITMAX, SSCALE ) 
         TYPE( quadrature_handle )  :: QUADRATURE 
         QDPREAL, OPTIONAL          :: SSCALE
         INTEGER                    :: ITMAX 
         PROFILING_START( 'adapt_convergence_parameter' )
        
         ! number of iterations for VARPRO
         IF ( QUADRATURE%N > 20 ) THEN
            ITMAX = 3000 
         ELSE
            ITMAX = 1500 
         ENDIF
          
         ! shinking scale for interval ( similar for POTIM )
         IF ( PRESENT( SSCALE ) ) THEN
            IF( QUADRATURE%ITYPE == -2 ) THEN
               SSCALE = QDPP(1)
            ELSE
               SSCALE = QDPP(0.95)
            ENDIF
         ENDIF

         PROFILING_STOP( 'adapt_convergence_parameter' )
      END SUBROUTINE ADAPT_CONVERGENCE_PARAMETER

      !**********************************************************************
      ! determines how strong interval is shrinked
      !**********************************************************************
      SUBROUTINE ADAPT_SCALING( ITER, SSCALE ) 
         INTEGER                    :: ITER
         QDPREAL                    :: SSCALE
         INTEGER, SAVE              :: ITER_OLD=-1
         PROFILING_START( 'adapt_scaling' )
        
         IF ( ITER_OLD < 0 ) THEN
            ! decrease SSCALE if quick convergence was reached
            IF ( ITER < 80 ) THEN
               SSCALE = SSCALE - QDPP(0.015)
            ELSE IF ( ITER < 100 ) THEN
               SSCALE = SSCALE - QDPP(0.01)
            ELSE IF ( ITER < 120 ) THEN
               SSCALE = SSCALE - QDPP(0.005)
            ELSE IF ( ITER < 140 ) THEN
               SSCALE = SSCALE - QDPP(0.0025)
            ELSE IF ( ITER < 170 ) THEN
               SSCALE = SSCALE - QDPP(0.0010)
            ELSE IF ( ITER < 190 ) THEN
               SSCALE = SSCALE - QDPP(0.0005)
            ELSE
               SSCALE = QDPP(1)
            ENDIF
         ELSE
            IF ( ITER_OLD > ITER ) THEN
               SSCALE = SSCALE - QDPP( 0.001 ) 
            ELSE IF ( ITER_OLD < ITER ) THEN
               SSCALE = SSCALE + QDPP( 0.0005 ) 
            ELSE IF ( ITER_OLD == ITER ) THEN
               SSCALE = SSCALE - QDPP( 0.0005 ) 
            ENDIF
            IF ( SSCALE > QDPP( 1 ) ) SSCALE = QDPP( 1 )
         ENDIF

         ITER_OLD = ITER 

         PROFILING_STOP( 'adapt_scaling' )
      END SUBROUTINE ADAPT_SCALING

      SUBROUTINE SET_RNEW( QUADRATURE, MINMAX, SCALING, IERR, IVERBOSITY, IO ) 
         TYPE( quadrature_handle )  :: QUADRATURE 
         QDPREAL                    :: MINMAX(2)
         QDPREAL                    :: SCALING
         INTEGER                    :: IVERBOSITY, IERR 
         TYPE( in_struct )          :: IO
         ! local 
         QDPREAL                    :: RNEW
         QDPREAL                    :: Y,X
         INTEGER                    :: LASTSIGN
         INTEGER                    :: NEWSIGN
         INTEGER                    :: NLOOP, I

         IF ( IERR <= 200 ) THEN
            RNEW = QUADRATURE%X0( 2*QUADRATURE%N )*SCALING
         ELSE
            ! remez algorithm used for shrinking, so using 
            ! last point that is of same sign as 2*QUADRATURE%N+1 
            LASTSIGN = SIGN_QD( ERROR_FUNCTION( QUADRATURE, &
               QUADRATURE%X0( 2*QUADRATURE%N+1) ) )

            NLOOP = 200 
10          CONTINUE

            sign_search: DO I = 1, NLOOP
               X = QUADRATURE%X0(2*QUADRATURE%N)+(I-1)*&
                (QUADRATURE%X0(2*QUADRATURE%N+1)-QUADRATURE%X0(2*QUADRATURE%N))/(NLOOP-1)
               NEWSIGN = SIGN_QD( ERROR_FUNCTION( QUADRATURE, X ) )
               IF ( NEWSIGN == LASTSIGN ) EXIT sign_search
            ENDDO sign_search

            IF ( I > NLOOP-1 ) THEN
               NLOOP = NLOOP * 2
               GOTO 10 
            ENDIF
            
            RNEW = X 
         ENDIF
         IF ( IO%NWRITE>IVERBOSITY ) THEN
            IF ( IO%IU0 >=0 ) &
            WRITE(*,2,ADVANCE="NO")TOPRINT( QUADRATURE%B ),&
               TOPRINT( RNEW ), TOPRINT( MINMAX(1) ),  TOPRINT( MINMAX(2) )
            IF ( IO%IU6 >=0 ) &
            WRITE(IO%IU6,3)TOPRINT( QUADRATURE%B ),&
               TOPRINT( RNEW ), TOPRINT( MINMAX(1) ),  TOPRINT( MINMAX(2) )
         ENDIF
         QUADRATURE%B = RNEW 

!2        FORMAT ('   ',4E12.4 )
2        FORMAT (' RNEW=',4E12.4 )
3        FORMAT (' RNEW=',4E12.4 )
   
      END SUBROUTINE SET_RNEW 

      !**********************************************************************
      !> this subroutine computes the minimax solution based using the 
      !> Remez algorithm
      !**********************************************************************
!#define PrintQuadratureMinimaxLoop
      SUBROUTINE REMEZ_ALGORITHM( QUADRATURE, IO, IERR ) 
         TYPE( quadrature_handle )  :: QUADRATURE 
         TYPE( in_struct )          :: IO
         INTEGER,OPTIONAL           :: IERR
         !parameter for Remez algorithm 
         INTEGER, PARAMETER         :: MAXREMEZ = 100
         INTEGER, PARAMETER         :: IVERBOSITY = 2
         INTEGER                    :: I
         INTEGER                    :: ITER
         INTEGER                    :: ITER_NL
         QDPREAL                    :: MINMAX(2)
         QDPREAL                    :: Y
         QDPREAL                    :: ACCURACY

         PROFILING_START( 'remez_algorithm' )

         ! early return in case of non-opitmized grids
         IF ( QUADRATURE%ITYPE < -100 ) THEN
            PROFILING_STOP( 'remez_algorithm' )
            RETURN
         ENDIF 

         IF ( IO%IU0>=0 .AND. IO%NWRITE>IVERBOSITY ) THEN
            WRITE(*,*)" "
            WRITE(*,'(A,I3,A,E10.4)')"  Looking for a minimax solution for ITYPE=", &
               QUADRATURE%ITYPE, " B=", TOPRINT(QUADRATURE%B)
         ENDIF

         CALL DETERMINE_ALTERNANT( QUADRATURE, IO )
         ! obtain minmax values
         MINMAX(1) = QDPP( 1E6 )
         MINMAX(2) =- QDPP( 1E6 )
         DO I = 1, 2*QUADRATURE%N+1
            Y = ERROR_FUNCTION( QUADRATURE, QUADRATURE%X0( I ) )
            MINMAX(1) = MIN( MINMAX(1) , ABS( Y ) )
            MINMAX(2) = MAX( MINMAX(2) , ABS( Y ) )
         ENDDO
         QUADRATURE%C(2*QUADRATURE%N+1) = MINMAX(2)
         ! required accuracy is two orders of magnitude lower
         ! than smallest error in interval 
         ACCURACY = MINMAX( 1 ) / 100

         ! print some info about convergence
         IF ( IO%IU0>=0 .AND. IO%NWRITE > IVERBOSITY ) THEN
            WRITE(*,1) TOPRINT(ACCURACY)
         ENDIF

         DO ITER=1, MAXREMEZ  !start the remez algorithm
#ifdef PrintQuadratureMinimaxLoop
            CALL PRINT_QUADRATURE( QUADRATURE, IO, &
               IT=ABS( QUADRATURE%ITYPE )*100+ITER, FNAME="quad-MML_R" )
            CALL PRINT_ERROR_FUNCTION( QUADRATURE, IO, &
               IT=ABS( QUADRATURE%ITYPE )*100+ITER, FNAME="erf-MML_R" )
#endif
            ! solve non linear system for new coefficients
            CALL SOLVE_NONLS( QUADRATURE, ITER_NL, IO )
         
            ! find alternant and maximum error 
            CALL DETERMINE_ALTERNANT( QUADRATURE, IO )
            ! obtain minmax values
            MINMAX(1) = QDPP( 1E6 )
            MINMAX(2) =- QDPP( 1E6 )
            DO I = 1, 2*QUADRATURE%N+1
               Y = ERROR_FUNCTION( QUADRATURE, QUADRATURE%X0( I ) )
               MINMAX(1) = MIN( MINMAX(1) , ABS( Y ) )
               MINMAX(2) = MAX( MINMAX(2) , ABS( Y ) )
            ENDDO
            QUADRATURE%INFINITYNORM = MINMAX(2)

            IF ( IO%IU0>=0 .AND. IO%NWRITE > IVERBOSITY ) THEN
               WRITE(*,2)ITER,TOPRINT( MINMAX(1) ), &
                  TOPRINT( MINMAX(2) ), &
                  TOPRINT( ABS(MINMAX(1)-MINMAX(2)) ), ITER_NL
            ENDIF
            ! algorithm converged
            IF ( ABS(MINMAX(1)-MINMAX(2)) < ACCURACY ) THEN
               EXIT 
            ENDIF
         ENDDO 

         IF ( ITER > MAXREMEZ ) THEN
            WRITE( *, *)' Warning Remez algorithm not converged, '//&
               'calculation continued with non-minimax solution'
         ENDIF
         
         IF ( PRESENT( IERR ) ) THEN 
            IERR = ITER 
         ENDIF
         PROFILING_STOP( 'remez_algorithm' )
1        FORMAT( '     iter      E_min       E_max       dError      iters_nl', E12.4)
2        FORMAT( '     ',I3,'     ',3E12.4,I10)
   
      END SUBROUTINE REMEZ_ALGORITHM
#undef PrintQuadratureMinimaxLoop

   END SUBROUTINE DETERMINE_QUADRATURE
#undef PrintQuadrature
#undef PrintQuadratureDecrease
#undef PrintQuadratureMinimax

   ! in case the MM solution can be calculated from 
   ! a polynomial fit 
!#define PrintMinimaxCheck
   SUBROUTINE COMPUTE_SOLUTION_FROM_FIT( QUADRATURE, B_REQUESTED, IERR, IO ) 
      TYPE( quadrature_handle )  :: QUADRATURE 
      QDPREAL                    :: B_REQUESTED
      INTEGER,INTENT(INOUT)      :: IERR
      TYPE( in_struct )          :: IO
      ! local                             
      INTEGER                    :: I
      REAL(q)                    :: LNR
      REAL(q)                    :: C(4*QUADRATURE%N+4)
      QDPREAL                    :: BOLD
      QDPREAL                    :: C1,C2
      LOGICAL                    :: LDECREASE=.FALSE.
      ! initialize ERROR 
      IERR = 0  
 
      ! in case fitting function is not associated, immediate return
      IF ( .NOT. ASSOCIATED( QUADRATURE%FITCOF ) ) THEN
         RETURN
      ENDIF 
      ! also return if requested R is not in interval
!      IF ( QUADRATURE%R1 > B_REQUESTED .OR. QUADRATURE%R2 < B_REQUESTED ) THEN
!         RETURN
!      ENDIF

      BOLD = QUADRATURE%B 
      ! otherwise continue 
      ! set B value of interval, 
      QUADRATURE%B = B_REQUESTED 

      ! compute coefficients and alternant from fit 
      LNR = TOREAL( LOG( B_REQUESTED ) )
      C = QUADRATURE%FITCOF( QUADRATURE%N, LNR  )

      ! check if LNR is in interval
      ! B_min = EXP( C( 4*N+3 ) ) and  B_max = EXP( C( 4*N+4 ) ) 
      ! fitting coefficients are valid for  B_min < B < B_max 

      IF ( C(4*QUADRATURE%N+3) > LNR .OR. C(4*QUADRATURE%N+4) < LNR ) THEN
         ! if not, set interval border to the one closest
         B_REQUESTED = QDPD( EXP( MIN( C(4*QUADRATURE%N+3), C(4*QUADRATURE%N+4) ) ) )
!WRITE(*,*)'BREQ', EXP( LNR ), B_REQUESTED, EXP( C(4*QUADRATURE%N+3 ) ), EXP( C(4*QUADRATURE%N+4)  ) 
         BOLD = QUADRATURE%B
         ! decrease of interval will be necessary 
         B_REQUESTED = BOLD 
         QUADRATURE%INFINITYNORM = QDPP( 0 ) 
         DO I = 1, 2*QUADRATURE%N+1
            QUADRATURE%C( I ) = QDPP( 0 )
            QUADRATURE%X0( I ) = QDPP( 0 ) 
         ENDDO
         IERR = 0 
         RETURN
      ENDIF

      DO I = 1, 2*QUADRATURE%N+1
         C1= EXP(C(I))
         IF ( C(2*QUADRATURE%N+1+I) < -99.99999999_q ) THEN
            C2= 0
         ELSE
            C2= EXP(C(2*QUADRATURE%N+1+I))
         ENDIF
         QUADRATURE%C( I ) = C1
         QUADRATURE%X0( I ) = C2
      ENDDO

      ! check if minimax solution can be found 
      CALL MINIMAX_CHECK( QUADRATURE, IO, IERR ) 
      ! set error 
      QUADRATURE%INFINITYNORM = ERROR_FUNCTION( QUADRATURE, B_REQUESTED )

      IF ( IERR > 0 ) THEN
         RETURN
      ENDIF

      ! if check went wrong, reset quadrature 
      QUADRATURE%B = BOLD
      QUADRATURE%INFINITYNORM = QDPP( 0 ) 
      DO I = 1, 2*QUADRATURE%N+1
         QUADRATURE%C( I ) = QDPP( 0 )
         QUADRATURE%X0( I ) = QDPP( 0 ) 
      ENDDO

      CONTAINS 

      SUBROUTINE MINIMAX_CHECK( QUADRATURE, IO, IERR ) 
         TYPE( quadrature_handle )  :: QUADRATURE 
         TYPE( in_struct )          :: IO
         INTEGER, INTENT(INOUT)     :: IERR
         !parameter for Remez algorithm 
         INTEGER, PARAMETER         :: MAXREMEZ = 100
         INTEGER, PARAMETER         :: IVERBOSITY = 4
         INTEGER                    :: I
         INTEGER                    :: ITER
         INTEGER                    :: ITER_NL
         QDPREAL                    :: MINMAX(2)
         QDPREAL                    :: Y
         QDPREAL                    :: ACCURACY

         ! obtain minmax values
         MINMAX(1) = QDPP( 1E6 )
         MINMAX(2) =- QDPP( 1E6 )
         DO I = 1, 2*QUADRATURE%N+1
            Y = ERROR_FUNCTION( QUADRATURE, QUADRATURE%X0( I ) )
            MINMAX(1) = MIN( MINMAX(1) , ABS( Y ) )
            MINMAX(2) = MAX( MINMAX(2) , ABS( Y ) )
         ENDDO
         QUADRATURE%C(2*QUADRATURE%N+1) = MINMAX(2)
         ! required accuracy is two orders of magnitude lower
         ! than smallest error in interval 
         !ACCURACY = MINMAX( 2 ) / 1000
         ACCURACY = VERYTINY

         ! print some info about convergence
         IF ( IO%IU0>=0 .AND. IO%NWRITE > IVERBOSITY ) THEN
            WRITE(*,1) TOPRINT(ACCURACY)
         ENDIF

         DO ITER=1, MAXREMEZ  !start the remez algorithm
            ! solve non linear system for new coefficients
            CALL SOLVE_NONLS( QUADRATURE, ITER_NL, IO )
#ifdef PrintMinimaxCheck
            CALL PRINT_QUADRATURE( QUADRATURE, IO, &
               IT=ABS( QUADRATURE%ITYPE )*100+ITER, FNAME="quad-check" )
            CALL PRINT_ERROR_FUNCTION( QUADRATURE, IO, &
               IT=ABS( QUADRATURE%ITYPE )*100+ITER, FNAME="erf-check_R" )
#endif
         
            ! no need to find alternant, correct one is already stored
            ! inside structure with sufficienty precission
!            CALL DETERMINE_ALTERNANT( QUADRATURE, IO )
            ! obtain minmax values
            MINMAX(1) = QDPP( 1E6 )
            MINMAX(2) =-QDPP( 1E6 )
            DO I = 1, 2*QUADRATURE%N+1
               Y = ERROR_FUNCTION( QUADRATURE, QUADRATURE%X0( I ) )
               MINMAX(1) = MIN( MINMAX(1) , ABS( Y ) )
               MINMAX(2) = MAX( MINMAX(2) , ABS( Y ) )
            ENDDO
            QUADRATURE%INFINITYNORM = MINMAX(2)

            IF ( IO%IU0>=0 .AND. IO%NWRITE > IVERBOSITY ) THEN
               WRITE(*,2)ITER,TOPRINT( MINMAX(1) ), &
                  TOPRINT( MINMAX(2) ), &
                  TOPRINT( ABS(MINMAX(1)-MINMAX(2)) ), ITER_NL
            ENDIF
            ! algorithm converged
            IF ( ABS(MINMAX(1)-MINMAX(2)) < ACCURACY ) THEN
               EXIT 
            ENDIF
         ENDDO 

         IERR = ITER 
         IF ( ITER > MAXREMEZ ) THEN
            WRITE( *, *)' Warning Remez algorithm not converged, '//&
               'calculation continued with non-minimax solution'
            IERR=0
         ENDIF
         
1        FORMAT( '     iter      E_min       E_max       dError      iters_nl', E12.4)
2        FORMAT( '     ',I3,'     ',3E12.4,I10)
   
      END SUBROUTINE MINIMAX_CHECK
     

   END SUBROUTINE COMPUTE_SOLUTION_FROM_FIT
#undef PrintMinimaxCheck

!****************************************************************************
!> DESCRIPTION:
!> determines transformation matrix from grid stored in QUAD_A to the 
!> one stored in QUAD_B
!****************************************************************************
#ifdef MinimaxTest
!#define PrintTransformationError
#endif
!#define PrintTransformationErrors
!#define PrintTransformationMatrix
!#define L2Norm
   SUBROUTINE DETERMINE_TRANSFORMATION( QUAD_A, QUAD_B, T, T_ERROR, LCONJG, &
      COMM, IO )
      USE minimax_functions2D, ONLY: UTIME, VTIME
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      TYPE( quadrature_handle )  :: QUAD_A
      TYPE( quadrature_handle )  :: QUAD_B
      REAL(q), POINTER           :: T(:,:)
      REAL(q), POINTER           :: T_ERROR(:)
      LOGICAL                    :: LCONJG
      TYPE( communic )           :: COMM 
      TYPE( in_struct )          :: IO
      ! local variables
      INTEGER                    :: I, J, K
      QDPREAL, ALLOCATABLE       :: XMAX(:)
      ! for LU decomposition
      INTEGER                    :: INFO          ! info variable
      INTEGER                    :: IPIV(QUAD_A%N)! pivot array  
      QDPREAL                    :: D ! auxillary variable for LU 
      INTEGER                    :: IND ! for loop blocking
      ! transformation arrays
      REAL(q), ALLOCATABLE       :: ALPHA(:,:) 
      REAL(q), ALLOCATABLE       :: ERRORS(:) 
      REAL(q), POINTER           :: OMEGA(:)
      QDPREAL                    :: OMEGA_K
      INTEGER                    :: NB 
      INTEGER                    :: MAXBISEC
      QDPREAL                    :: MAXERROR
      QDPREAL                    :: A(QUAD_A%N,QUAD_A%N) ! design matrix of fitting problem
      QDPREAL                    :: ALPHA_K(QUAD_A%N) ! function values and solution
      INTEGER                    :: NODE_ME = 1  
      INTEGER, PARAMETER         :: IVERBOSITY = 3
#ifdef L2Nornm
      INTEGER, PARAMETER         :: NSCALE = 10
#else
      INTEGER, PARAMETER         :: NSCALE = 1
#endif
      QDPREAL, ALLOCATABLE       :: X(:)
      PROCEDURE( basis_function ), POINTER :: QUAD_B_PHI => NULL() 
      PROCEDURE( basis_function ), POINTER :: QUAD_B_PSI => NULL() 
      PROCEDURE( object_function ), POINTER :: QUAD_B_TRANS_TAYLOR => NULL() 
    
      PROFILING_START( 'DETERMINE_TRANSFORMATION' )
#ifdef MPI 
      NODE_ME = COMM%NODE_ME
#endif
      IF ( LCONJG ) THEN    
         ! dual basis function must be assigned
         IF ( .NOT. ASSOCIATED( QUAD_B%PSI_CONJG )  ) THEN
            CALL vtutor%error("ERROR, dual basis function not assigned (c)")
         ENDIF
         ! basis function must be assigned
         IF ( .NOT. ASSOCIATED( QUAD_B%PHI_CONJG )  ) THEN
            CALL vtutor%error("ERROR, basis function not assigned (c)")
         ENDIF
         ! for low frequencies
         IF ( .NOT. ASSOCIATED( QUAD_B%TRANS_TAYLOR_CONJG )  ) THEN
            CALL vtutor%error("ERROR, dual basis function for small arguments not assigned (c)")
         ENDIF

         QUAD_B_PHI => QUAD_B%PHI_CONJG
         QUAD_B_PSI => QUAD_B%PSI_CONJG 
         QUAD_B_TRANS_TAYLOR => QUAD_B%TRANS_TAYLOR_CONJG
      ELSE
         ! dual basis function must be assigned
         IF ( .NOT. ASSOCIATED( QUAD_B%PSI )  ) THEN
            CALL vtutor%error("ERROR, dual basis function not assigned")
         ENDIF
         ! basis function must be assigned
         IF ( .NOT. ASSOCIATED( QUAD_B%PHI )  ) THEN
            CALL vtutor%error("ERROR, basis function not assigned")
         ENDIF
         ! for low frequencies
         IF ( .NOT. ASSOCIATED( QUAD_B%TRANS_TAYLOR )  ) THEN
            CALL vtutor%error("ERROR, dual basis function for small arguments not assigned")
         ENDIF

         QUAD_B_PHI => QUAD_B%PHI
         QUAD_B_PSI => QUAD_B%PSI
         QUAD_B_TRANS_TAYLOR => QUAD_B%TRANS_TAYLOR
      ENDIF
  
      ! allocate transformation array T and T_ERROR 
      CALL ALLOCATE_TRANSFORMATION( T, T_ERROR, QUAD_A, QUAD_B, &
         OMEGA, NB )

      ! early return in case of non-optimized Matsubara grid
      IF ( QUAD_A%ITYPE < -100 .AND. QUAD_B%ITYPE < -100  ) THEN
         CALL SET_DFT_MATRIX( T, T_ERROR, QUAD_A, QUAD_B )

         PROFILING_STOP( 'DETERMINE_TRANSFORMATION' )
         RETURN
      ENDIF

      MAXBISEC = SET_MAXBISEC( QUAD_A%INFINITYNORM )
      ALLOCATE( XMAX( MAXBISEC ) )  
      ALLOCATE( X( NSCALE*NDATALSQ ) )
      DO I = 1, MAXBISEC 
         XMAX(I) = QDPP( 0 )
      ENDDO
      ! allocate sampling points for minimization interval 
      CALL SAMPLING_FROM_ALTERNANT_FT( QUAD_A, NDATALSQ*NSCALE, X )
      CALL SAMPLING_FROM_ALTERNANT_FT( QUAD_A, MAXBISEC, XMAX )

      ! set up A_jl = \sum_i^NDATALSQ \psi_j(x_i)\psi_l(x_i)
      DO J = 1, QUAD_A%N
         DO K = 1, QUAD_A%N
            A( J, K )=0
#ifdef L2Norm
            A( J, K )=(X(1)-QUAD_A%A)*QUAD_B_PSI( X(1), QUAD_A%C(J) )*&
                                      QUAD_B_PSI( X(1), QUAD_A%C(K) )
            DO I = 2, NDATALSQ*NSCALE
                A( J, K ) = A( J, K ) + (X(I)-X(I-1))*QUAD_B_PSI( X(I), QUAD_A%C(J) )*&
                                        QUAD_B_PSI( X(I), QUAD_A%C(K) )
#else
            DO I = 1, NDATALSQ*NSCALE
                A( J, K ) = A( J, K ) + QUAD_B_PSI( X(I), QUAD_A%C(J) )*&
                                        QUAD_B_PSI( X(I), QUAD_A%C(K) )
#endif
            ENDDO
         ENDDO
      ENDDO

      ! obtain its LU decomposion and store U in A 
      CALL QDLUDCMP( A, QUAD_A%N, QUAD_A%N, IPIV, D, INFO )
      IF (INFO /= 0 ) THEN
         CALL vtutor%error("ERROR in DETERMINE_TRANSFORMATION: QDLUDCMP failed with: " // str(INFO))
      ENDIF 

      !allocate storage for transformation and its errors
      ALLOCATE( ALPHA( NB, QUAD_A%N ), ERRORS( NB ) )
      ALPHA = 0
      ERRORS = 0

      ! blocked frequency loop
      DO K = 1, NB
         OMEGA_K = QDPD( OMEGA( K ) )
#ifdef MPI
         ! block this loop with MPI 
         IND = LOCAL_INDEX_CYCLIC( COMM%NCPU, NB, COMM%NODE_ME, K )
         IF (IND < 0 ) CYCLE
#endif
         ! use time integration for small frequency points
         ! taylor expansion only if more than one point is computed
         IF ( ABS( OMEGA(K) ) < 0.01_q .AND. NB>1 ) THEN
            DO J = 1, QUAD_A%N 
               ALPHA_K( J ) = QDPP(2)*QUAD_A%C(J+QUAD_A%N)*QUAD_B_TRANS_TAYLOR( OMEGA_K*QUAD_A%C( J ) )
            ENDDO
            ! alternative scaling for finite temperature 
            ERRORS( K ) = QUAD_A%INFINITYNORM
         ELSE
            ! set up b vector for given frequency 
            DO J = 1, QUAD_A%N
               ALPHA_K( J ) = QDPP( 0 )

#ifdef L2Norm
               ALPHA_K( J ) = (QUAD_A%A-X(1))*QUAD_B_PHI( X(1), OMEGA_K)* &
                                    QUAD_B_PSI( X(1), QUAD_A%C(J) ) 
               DO I = 2, NDATALSQ*NSCALE
                  ALPHA_K( J ) = ALPHA_K( J ) + (X(I)-X(I-1))*QUAD_B_PHI( X(I), OMEGA_K)* &
                                    QUAD_B_PSI( X(I), QUAD_A%C(J) ) 
#else
               DO I = 1, NDATALSQ*NSCALE
                  ALPHA_K( J ) = ALPHA_K( J ) + QUAD_B_PHI( X(I), OMEGA_K)* &
                                    QUAD_B_PSI( X(I), QUAD_A%C(J) ) 
#endif

               ENDDO 
            ENDDO 

            !solve linear equation, where A is in upper triangluar form
            CALL QDLUBKSB( A, QUAD_A%N, QUAD_A%N, IPIV, ALPHA_K )
            
            ! find maximum value of error function
#ifdef L2Norm
            MAXERROR  = L2_NORM(QUAD_A, QUAD_B_PHI, QUAD_B_PSI, ALPHA_K, &
               XMAX, OMEGA_K,  MAXBISEC, IO%IU0)
#else
            MAXERROR  = CHEBYSHEV_NORM(QUAD_A, QUAD_B_PHI, QUAD_B_PSI, ALPHA_K, &
               XMAX, OMEGA_K,  MAXBISEC, IO%IU0)
#endif

            ERRORS(K)=MAXERROR
         ENDIF

         !save fourier matrix
         DO I=1,QUAD_A%N
            ALPHA( K, I ) = ALPHA_K( I )
         ENDDO

#ifdef PrintTransformationErrors
      ! this prints the error function of the transform 
      CALL PRINT_TRANS_ERROR_FUNCTION( QUAD_A, QUAD_B, ALPHA_K, OMEGA_K, &
         LCONJG, NODE_ME, IT = QUAD_B%ITYPE*100+SIGN(K,QUAD_B%ITYPE))
#endif
      ENDDO  !frequency loop 

      ! MPI communication 
      CALLMPI( M_sum_d(COMM, ALPHA, NB*QUAD_A%N))
      CALLMPI( M_sum_d(COMM, ERRORS, NB))

      IF ( IO%IU0 >=0 .AND. IO%NWRITE>IVERBOSITY ) THEN
         WRITE(*,1)QUAD_A%ITYPE,QUAD_B%ITYPE
         DO K = 1, NB 
            WRITE(*,2)K,OMEGA(K),ERRORS(K)
         ENDDO
#ifdef PrintTransformationError
         CALL PRINT_TRANSFORMATION_ERROR( OMEGA, ERRORS, QUAD_B%N, LCONJG, &
           NODE_ME, IT=QUAD_B%ITYPE)
#endif
      ENDIF

#ifdef PrintTransformationMatrix
      ! this prints the error function of the transform 
      CALL PRINT_TRANSFORMATION_MATRIX( QUAD_A, QUAD_B, ALPHA, &
         LCONJG, NODE_ME, IT = QUAD_B%ITYPE*100+SIGN(1,QUAD_B%ITYPE))
#endif

      ! store Transformation matrices 
      DO J = 1, QUAD_A%N
         DO K = 1, NB 
            T( K, J ) = ALPHA( K, J ) 
         ENDDO
      ENDDO
      DO K = 1, NB 
         T_ERROR( K) = ERRORS( K ) 
      ENDDO

      IF ( ALLOCATED( XMAX ) ) DEALLOCATE( XMAX )
      IF ( ALLOCATED( ERRORS) ) DEALLOCATE( ERRORS )
      IF ( ASSOCIATED( OMEGA ) ) DEALLOCATE( OMEGA )

      PROFILING_STOP( 'DETERMINE_TRANSFORMATION' )
1     FORMAT( '    Transformation Errors from grid ',I3,' to grid ',I3,/&
              '    Index       Omega_k       Error' )
2     FORMAT( '     ',I8,2E12.4 )
      CONTAINS

      ! ******************************************************************
      ! allocates transformation matrix T and corresponding error array 
      ! T_ERROR based on number of grid points found in QUAD_A and QUAD_B
      ! ******************************************************************
      SUBROUTINE ALLOCATE_TRANSFORMATION( T, T_ERROR, QUAD_A, QUAD_B, OMEGA, NB )
         USE string, ONLY: str
         USE tutor, ONLY: vtutor
         REAL(q), POINTER           :: T(:,:)
         REAL(q), POINTER           :: T_ERROR(:)
         TYPE( quadrature_handle )  :: QUAD_A
         TYPE( quadrature_handle )  :: QUAD_B
         REAL(q), POINTER           :: OMEGA(:)
         INTEGER                    :: NB
         ! local                    
         INTEGER                    :: I , J
#ifdef PrintTransformationError
         INTEGER                    :: IDBL
         INTEGER                    :: IDBL_MAX 
         QDPREAL                    :: X, X2, X1
         QDPREAL, ALLOCATABLE       :: X_I(:),X_S(:), X_M(:)
#endif 
         PROFILING_START('allocate_transformation')  

         ! quadratures have to be set 
         IF ( QUAD_A%N == 0 .OR. QUAD_B%N == 0 ) THEN
            CALL vtutor%error("ERROR, non-initialized handles passed to ALLOCATE_TRANSFORMATION " // &
               str(QUAD_A%N) // " " // str(QUAD_B%N))
         ENDIF    
         ! number of grid points in domains B
         NB = QUAD_B%N
         ! some bosonic grids require to add zero frequency point
         ! using the time integration quadrature weights
         IF( QUAD_B%ITYPE == -2 .OR. QUAD_B%ITYPE == -20 ) THEN
            NB = NB+1
         ENDIF 
#ifdef PrintTransformationError
         ! set number of sampling points dynamically based
         ! on frequency points
         ALLOCATE ( X_I( QUAD_B%N ) )
         DO I = 1, QUAD_B%N
            X_I( I ) = QUAD_B%B * QUAD_B%C( I ) / QUAD_B%C( QUAD_B%N )
         ENDDO
         IDBL_MAX = 20
         ! maximum number of doublications of sampling points 
         dbl_points: DO IDBL = 1, IDBL_MAX
            IF( ALLOCATED( X_S ) ) DEALLOCATE( X_S )
            ALLOCATE( X_S( SIZE( X_I ) ) )
            ! create sampling points based on intersection points
            X_S( 1 ) = QUAD_B%A + X_I( 1 ) / 2 
            DO I = 2, SIZE( X_I ) 
              X_S( I ) = QUAD_B%A + (X_I( I ) + X_I( I-1 )) / 2 
            ENDDO
            IF ( SIZE( X_S ) > 10000 ) EXIT

            ALLOCATE( X_M( 2*SIZE( X_I ) ) )
            DO I = 1, SIZE( X_I )
               X_M( 2*I ) = X_I( I ) 
               X_M( 2*I - 1 ) =  X_S( I ) - QUAD_A%A
            ENDDO
            IF ( ALLOCATED(X_I) ) DEALLOCATE( X_I )
            ALLOCATE( X_I( SIZE( X_M ) ) )
            DO I = 1, SIZE( X_M )
               X_I( I ) = X_M( I ) 
            ENDDO
            DEALLOCATE( X_M )
         ENDDO dbl_points

         NB = SIZE( X_S ) + QUAD_B%N 
         IF( IO%IU0>=0 .AND. IO%NWRITE>IVERBOSITY ) WRITE(*,*)"  used number of sampling points:",  NB
         IF( ALLOCATED( X_I ) ) DEALLOCATE( X_I )
         IF( ALLOCATED( X_M ) ) DEALLOCATE( X_M )
#endif       
         ! allocate storage for transformation matrices and its errors
         NULLIFY( T )
         NULLIFY( T_ERROR )
         NULLIFY( OMEGA )
         IF( .NOT. ASSOCIATED( T ) ) ALLOCATE( T( NB, QUAD_A%N ) )
         IF( .NOT. ASSOCIATED( T_ERROR ) ) ALLOCATE( T_ERROR( NB ) )
         IF( .NOT. ASSOCIATED( OMEGA ) ) ALLOCATE( OMEGA( NB ) )

         DO J = 1, NB
            DO I = 1, QUAD_A%N
               T(J,I) = QDPP( 0 )
            ENDDO
         ENDDO
         DO J = 1, NB
            T_ERROR(J) = QDPP( 0 )
         ENDDO

         ! set quadrature points for which transformation matrix 
         ! is calculated
         IF( QUAD_A%N == NB .OR. (QUAD_B%N == NB .AND. QUAD_B%N < QUAD_A%N) ) THEN
            IF (QUAD_B%ITYPE == -2 .OR. QUAD_B%ITYPE == -20 ) THEN
               OMEGA(1) = 0._q
               IF ( QUAD_B%N > 1 ) THEN
                  DO I = 2, NB
                     OMEGA( I ) = TOREAL( QUAD_B%C( I-1 ) )
                  ENDDO
               ENDIF
            ELSE 
               DO I = 1, NB
                  OMEGA( I ) = TOREAL( QUAD_B%C( I ) )
               ENDDO
            ENDIF
#ifdef PrintTransformationError
         ! this should never be executed, except for debuging
         ! fine sampling of frequency axis based on frequency abszissa
         ELSE 
            J=0
            DO I = 1, NB - QUAD_B%N 
               OMEGA( I ) = TOREAL( X_S( I ) )
               J=J+1
            ENDDO
            IF( ALLOCATED( X_S ) ) DEALLOCATE( X_S )
            ! last points correspond to optimized frequency quadrature
            DO I = 1, QUAD_B%N 
               OMEGA( J + I ) = TOREAL( QUAD_B%C( I ) )
            ENDDO
#else
         ELSE
            CALL vtutor%error("ERROR, ALLOCATE_TRANSFORMATION reports too many points " // str(NB) // &
               " " // str(QUAD_A%N))
#endif
         ENDIF

         PROFILING_STOP('allocate_transformation')  
      END SUBROUTINE ALLOCATE_TRANSFORMATION

      ! ******************************************************************
      ! allocates transformation matrix T and corresponding error array 
      ! T_ERROR based on number of grid points found in QUAD_A and QUAD_B
      ! ******************************************************************
      SUBROUTINE SET_DFT_MATRIX( T, T_ERROR, QUAD_A, QUAD_B )
         USE constant
         REAL(q)                    :: T(:,:)
         REAL(q)                    :: T_ERROR(:)
         TYPE( quadrature_handle )  :: QUAD_A
         TYPE( quadrature_handle )  :: QUAD_B
         ! local                    
         INTEGER                    :: I , J
         PROFILING_START('set_dft_matrix')  

         ! set Bosonic Fourier Transformation 
         DO I = 1, QUAD_A%N
            DO J = 1, QUAD_B%N
               IF ( LCONJG ) THEN    
                  T(J,I)= SIN( QUAD_B%C(J)*QUAD_A%C(I) )/QUAD_A%N
               ELSE
                  T(J,I)= COS( QUAD_B%C(J)*QUAD_A%C(I) )/QUAD_A%N
               ENDIF
            ENDDO
         ENDDO
         PROFILING_STOP('set_dft_matrix')  
      END SUBROUTINE SET_DFT_MATRIX

      ! ******************************************************************
      ! allocates transformation matrix T and corresponding error array 
      ! T_ERROR based on number of grid points found in QUAD_A and QUAD_B
      ! ******************************************************************
      SUBROUTINE SAMPLING_FROM_ALTERNANT_FT( QUADRATURE, NDATA, T )
         USE minimax_varpro, ONLY: PIKSRT2
         USE string, ONLY: str
         USE tutor, ONLY: vtutor
         TYPE( quadrature_handle )  :: QUADRATURE
         QDPREAL                    :: T(:)
         INTEGER                    :: NDATA
         !local 
         INTEGER                    :: I,K,L
         QDPREAL, ALLOCATABLE       :: XI(:)
 
         PROFILING_START( 'sampling_from_alternant_ft' ) 

         K = 2*QUADRATURE%N+1
         ALLOCATE( XI( K ) ) 
         DO I = 1, K 
            XI(I) = QUADRATURE%X0( I )
            T(I) = XI(I)
         ENDDO
         doubling: DO L = 1, 1000
            ! add points inbetween
            DO I = 1, 2*QUADRATURE%N-1
               IF ( I+K >NDATA ) EXIT doubling
               T( I + K ) = ( XI(I+1) + XI(I) ) / QDPP(2)
            ENDDO
            K = K + 2*QUADRATURE%N-1
            IF ( ALLOCATED( XI ) ) DEALLOCATE( XI ) 
            ALLOCATE( XI( K ) ) 
            XI( 1:K ) = T( 1:K )
            CALL PIKSRT2( K, XI ) 
         ENDDO doubling

         CALL PIKSRT2( NDATA, T(1:NDATA) ) 

         IF ( L > 1000 ) THEN
            CALL vtutor%error("error in SAMPLING_FROM_ALTERNANT_FT " // str(L))
         ENDIF
         IF ( ALLOCATED( XI ) ) DEALLOCATE( XI ) 

         PROFILING_STOP( 'sampling_from_alternant_ft' ) 
      END SUBROUTINE SAMPLING_FROM_ALTERNANT_FT

      !*******************************************************************
      ! sets the number of bisection points for sampling the axis
      ! used to determine the error of the transformation for a single 
      ! frequency point 
      !*******************************************************************
         FUNCTION SET_MAXBISEC( ERROR ) 
            USE prec
            IMPLICIT NONE
            INTEGER SET_MAXBISEC 
            QDPREAL :: ERROR
            PROFILING_START('set_maxbisec')
               
            IF ( ERROR > QDPP(1E-6) ) THEN
               SET_MAXBISEC = 200 
            ELSEIF ( ERROR > QDPP(1E-7) ) THEN
               SET_MAXBISEC = 300 
            ELSEIF ( ERROR > QDPP(1E-8) ) THEN
               SET_MAXBISEC = 400 
            ELSEIF ( ERROR > QDPP(1E-9) ) THEN
               SET_MAXBISEC = 400 
            ELSEIF ( ERROR > QDPP(1E-10) ) THEN
               SET_MAXBISEC = 500 
            ELSEIF ( ERROR > QDPP(1E-11) ) THEN
               SET_MAXBISEC = 600 
            ELSEIF ( ERROR > QDPP(1E-12) ) THEN
               SET_MAXBISEC = 700 
            ELSEIF ( ERROR > QDPP(1E-13) ) THEN
               SET_MAXBISEC = 800 
            ELSEIF ( ERROR > QDPP(1E-14) ) THEN
               SET_MAXBISEC = 900 
            ELSEIF ( ERROR > QDPP(1E-15) ) THEN
               SET_MAXBISEC = 1000 
            ELSE
               SET_MAXBISEC = 2000 
            ENDIF 

            PROFILING_STOP('set_maxbisec')
      
            RETURN
         ENDFUNCTION SET_MAXBISEC

      !*******************************************************************
      ! Calculates the chebyshev norm of frequency transformation function
      ! this function is defined as 
      ! \f$
      !   \eta( x, \vec{ \alpha }_k) = \phi_b( x, \omega_k ) -
      !   \sum_{j=1}^N \alpha_kj \phi_a( x, \tau_j )
      ! \f$
      !gj******************************************************************
      FUNCTION CHEBYSHEV_NORM( QUAD_A, QUAD_B_PHI, QUAD_B_PSI, ALPHA_K, XMAX,&
         OMEGA_K, MAXBISEC, IU0 )
         USE prec
         QDPREAL                     :: CHEBYSHEV_NORM
         TYPE( quadrature_handle )   :: QUAD_A
         PROCEDURE( basis_function ) :: QUAD_B_PHI
         PROCEDURE( basis_function ) :: QUAD_B_PSI
         QDPREAL                     :: ALPHA_K(:)
         QDPREAL                     :: XMAX(:)
         QDPREAL                     :: OMEGA_K
         INTEGER                     :: MAXBISEC
         INTEGER                     :: IU0
         !local variables 
         INTEGER                     :: I, J
         QDPREAL                     :: E_0, E_1
   
         PROFILING_START('chebyshev_norm')
   
         E_0 = QUAD_B_PHI( XMAX( 1 ), OMEGA_K )
         DO J = 1, QUAD_A%N
            E_0 = E_0 - ALPHA_K( J )*QUAD_B_PSI( XMAX( 1 ), QUAD_A%C( J ) )
         ENDDO
         E_0 = ABS( E_0 )
          
         DO I=2,MAXBISEC
            E_1 = QUAD_B_PHI( XMAX( I ), OMEGA_K )
            DO J = 1, QUAD_A%N
               E_1 = E_1 - ALPHA_K( J )*QUAD_B_PSI( XMAX( I ), QUAD_A%C( J ) )
            ENDDO
            E_0 = MAX( E_0, ABS( E_1 ) )
         ENDDO
      
         CHEBYSHEV_NORM=E_0
      
         PROFILING_STOP('chebyshev_norm')
      
      ENDFUNCTION CHEBYSHEV_NORM

      !*******************************************************************
      ! Calculates the L2 norm of frequency transformation function
      ! this function is defined as 
      ! \f$
      !   \eta( x, \vec{ \alpha }_k) = \phi_b( x, \omega_k ) -
      !   \sum_{j=1}^N \alpha_kj \phi_a( x, \tau_j )
      ! \f$
      !gj******************************************************************
      FUNCTION L2_NORM( QUAD_A, QUAD_B_PHI, QUAD_B_PSI, ALPHA_K, XMAX,&
         OMEGA_K, MAXBISEC, IU0 )
         USE prec
         QDPREAL                     :: L2_NORM
         TYPE( quadrature_handle )   :: QUAD_A
         PROCEDURE( basis_function ) :: QUAD_B_PHI
         PROCEDURE( basis_function ) :: QUAD_B_PSI
         QDPREAL                     :: ALPHA_K(:)
         QDPREAL                     :: XMAX(:)
         QDPREAL                     :: OMEGA_K
         INTEGER                     :: MAXBISEC
         INTEGER                     :: IU0
         !local variables 
         INTEGER                     :: I, J
         QDPREAL                     :: E_0
         QDPREAL                     :: L2
   
         PROFILING_START('l2_norm')
   
         E_0 = QUAD_B_PHI( XMAX( 1 ), OMEGA_K )
         DO J = 1, QUAD_A%N
            E_0 = E_0 - ALPHA_K( J )*QUAD_B_PSI( XMAX( 1 ), QUAD_A%C( J ) )
         ENDDO
         L2 = E_0*E_0*(XMAX(2)-XMAX(1))
          
         DO I=2,MAXBISEC-1
            E_0 = QUAD_B_PHI( XMAX( I ), OMEGA_K )
            DO J = 1, QUAD_A%N
               E_0 = E_0 - ALPHA_K( J )*QUAD_B_PSI( XMAX( I ), QUAD_A%C( J ) )
            ENDDO
            L2 = L2+(XMAX(I+1)-XMAX(I))*E_0*E_0
         ENDDO
      
         L2_NORM=L2
      
         PROFILING_STOP('l2_norm')
      
      ENDFUNCTION L2_NORM

   END SUBROUTINE DETERMINE_TRANSFORMATION

!****************************************************************************
!>helper, sets the boundary B of the requested interval such that convergence
!>is achieved for requested N, or tabulated coefficients are set if found 
!****************************************************************************

   SUBROUTINE ADJUST_BOUNDARY_FOR_CONVERGENCE( QUADRATURE, RSTART, IO ) 
      USE constant
      TYPE( quadrature_handle ) :: QUADRATURE 
      QDPREAL                   :: RSTART     ! energy interval 
      TYPE(in_struct)           :: IO         ! for output and input
      !maximum minimization interval for each order N
      REAL(q)                   :: RMAX
      INTEGER                   :: IA 
      PROFILING_START( 'ADJUST_BOUNDARY_FOR_CONVERGENCE' )

      ! in case tabulated coefficients are present 
      IF ( ASSOCIATED( QUADRATURE%CTAB ) )  THEN
         ! set tabulated coefficients
         CALL INIT_COEFFICIENTS( QUADRATURE, RSTART, IO ) 
         IA = 0
      ELSE IF ( QUADRATURE%ITYPE < -90 )  THEN
         ! Matsubara gri 
         CALL SET_QUADRATURE_NON_OPTIMIZED( QUADRATURE ) 
         IA = 0
      ELSE
         ! set empirically found value for right boundary value
         ! this should work for T=0 grids
         CALL SET_RSTART_EMPIRICAL(QUADRATURE%N,RSTART,IO%IU0)
         IA = 1
      ENDIF
      QUADRATURE%A = QDPD( IA ) 
      QUADRATURE%B = RSTART 

      PROFILING_STOP( 'ADJUST_BOUNDARY_FOR_CONVERGENCE' )
      CONTAINS 
      
      SUBROUTINE SET_RSTART_EMPIRICAL(N,RSTART,IU0)
         INTEGER                   :: N          !order 
         QDPREAL                   :: RSTART
         INTEGER                   :: IU0          
         !local
         REAL(q)                   :: RMAX
         REAL(q)                   :: R          !current minimization interval 
         INTEGER                   :: I, J         
         REAL(q)                   :: RNMAX(32)    
         DATA RNMAX(1:32)/&      
            8.667E+0,& ! 1  
            4.154E+1,& ! 2  
            1.468E+2,& ! 3 
            4.360E+2,& ! 4 
            1.153E+3,& ! 5 
            2.807E+3,& ! 6 
            6.373E+3,& ! 7 
            1.375E+4,& ! 8 
            2.839E+4,& ! 9 
            5.650E+4,& ! 10 
            1.089E+5,& ! 11 
            2.042E+6,& ! 12 
            3.737E+5,& ! 13 
            6.691E+5,& ! 14 
            1.175E+6,& ! 15 
            2.027E+6,& ! 16 
            3.440E+6,& ! 17 
            5.753E+6,& ! 18 
            9.491E+6,& ! 19 
            1.546E+7,& ! 20
            2.491E+7,& ! 21  
            3.969E+7,& ! 22 
            6.258E+7,& ! 23 
            9.776E+7,& ! 24 
            1.513E+8,& ! 25 
            2.325E+8,& ! 26 
            3.540E+8,& ! 27 
            5.353E+8,& ! 28 
            8.036E+8,& ! 29 
            1.198E+9,& ! 30
            1.775E+9,& ! 31
            2.614E+9/  ! 32
         REAL(q)                   :: CHECKPT(6) !subdivisions 
         INTEGER                   :: IVAL(6)    !supported orders of each subdivision
         !interval intersection (found empirically for moderate band gaps)
         DATA CHECKPT(1:6)/25._q, 50._q, 106._q, 221._q, 1005._q, 10005._q /
         !maximum number of grid points supported for intersection
         DATA IVAL(1:6)/ 12, 14, 16, 20, 24, 28/

         ! interval intersections for metals at low temperature 
         REAL(q) :: RMETAL(39)   !intervals that converge for cos sin and exp 
         DATA RMETAL(1:39)/&
         30005,   80005, 110005, 120501, 130005, 140005, 150005, 160005, 170005, 180005, &
         200005, 220005, 230501, 240005, 250005, 260005, 275551, 280005, 370005, 400501, &
         420005, 440005, 480005, 530005, 560005, 660005, 730005, 800005, 930005, 1000005,& 
         1200005, 1400005, 1600005, 2000005, 2200005, 2400005, 2600005, 3000005, 20000005/

         PROFILING_START('SET_RSTART_EMPIRICAL')
             
         R = RSTART

         !> For each N there is a maximum interval [1,RNMAX]. 
         !! Coefficients of [1,R>RNMAX] coincide with the coefficients of [1,RNMAX]
         R = MIN( R, RNMAX(N) ) 

         ! nothing to do if R is close to RNMAX( N )
         IF ( ABS( R - RNMAX( N ) )/MAX( R, RNMAX(N) )  <= 0.25_q) THEN
            RSTART = RNMAX( N ) 
            PROFILING_STOP( 'SET_RSTART_EMPIRICAL' )
            RETURN
         ENDIF

         !adjust RMAX according to R, use the next larger interval 
         !find correct RMAX 
         DO I = 1, SIZE(CHECKPT)
            RMAX = CHECKPT(I)
            IF ( R < RMAX) EXIT
         ENDDO
         ! check if requested R supports given N,
         IF ( N > IVAL(I) .AND. I<SIZE(CHECKPT) ) THEN
            ! if not go to next larger interval, which supports this N
            rmax_n: DO J=I+1, SIZE(CHECKPT)
               IF ( N<= IVAL(J) ) THEN
                  RMAX=CHECKPT(J)
                  EXIT rmax_n
               ENDIF
            ENDDO rmax_n
            I=J
         ENDIF 

         ! in case interval is larger than the ones listed above (often for metallic systems)
         ! go through following checkpoints defined below
         ! the minimization is started for following intervals, 
         ! where convergence was tested for cos, sin and exp with up to N=32 grid points
         !
         ! this should converge for all intervals with a band gap of > 0.0005 eV
         ! and a transition energy of up to 1500 eV 
         ! the smallest interval would be therefore be 
         ! [ 0.149, 1500 ] -> R ~ 10005
         ! the largest interval would be 
         ! [ 0.0005, 1500 ] -> R ~ 3E+6
         !

         IF ( I>SIZE(CHECKPT) .OR. N > IVAL( I )  )THEN
            ! check if R lies in between any of the following intervals
            search: DO I = 1, SIZE(RMETAL)-1
               ! if R is between two converged intervals, use larger one
               IF( RMETAL( I ) < R .AND. R < RMETAL( I + 1 ) ) THEN
                  RMAX = RMETAL( I+1 ) 
                  EXIT search
               ELSE IF ( ABS( RMETAL( I ) - R ) < ACC2 )  THEN
                  RMAX = RMETAL( I ) 
                  EXIT search
               ELSE IF ( ABS( RMETAL( I+1 ) - R ) < ACC2 ) THEN
                  RMAX = RMETAL( I+1 ) 
                  EXIT search 
               ENDIF
               RMAX = MAX( R, RMETAL(1) )
            ENDDO search

            ! nevertheless revert back to maximum grid defined by RNMAX(N) 
            ! this is the maximum error possible for a infinitesimal band gap 
            RMAX = MIN( RNMAX( N ) , RMAX ) 
            ! note that in principle we could start every minimization from RNMAX(N)
            ! for a given grid order N, but the procedure above accelerates 
            ! convergence since we start the minimization from a closer interval
         ENDIF
         RSTART = RMAX 

         PROFILING_STOP('SET_RSTART_EMPIRICAL')
      END SUBROUTINE SET_RSTART_EMPIRICAL

      !**********************************************************************
      !> direct calculation of quadrature
      !> this includes: bosonic, fermionic and hypergeometric matsubara grid
      !**********************************************************************
      SUBROUTINE SET_QUADRATURE_NON_OPTIMIZED( QUAD ) 
         USE string, ONLY: str
         USE tutor, ONLY: vtutor
         TYPE( quadrature_handle ) :: QUAD       
         ! local 
         INTEGER                  :: I 
         QDPREAL                  :: EA, EB

         PROFILING_START( 'set_quadrature_non_optimized' )

         ! uniform time grid 
         IF ( QUAD%ITYPE == -101 ) THEN
            DO I = 1, QUAD%N
               QUAD%C( I ) = QDPD(I-0.5)/QDPD(2*QUAD%N ) 
               QUAD%C( I+QUAD%N ) = QDPP(1)/(QDPD(2*QUAD%N))
            ENDDO

         ! Bosonic Matsubara grid 
         ELSE IF ( QUAD%ITYPE == -102 ) THEN
            DO I = 1, QUAD%N
               QUAD%C( I ) = QDPP(2)*PIQ*( QDPD( I ) - QDPP(1) ) 
               QUAD%C( I+QUAD%N ) = QDPP(2)
            ENDDO
            ! zero frequency point appears only once
            QUAD%C( 1+QUAD%N ) = QDPP(1)

         ! Fermionic Matsubara grid 
         ELSE IF ( QUAD%ITYPE == -103 ) THEN
            DO I = 1, QUAD%N
               QUAD%C( I ) = PIQ*( QDPP(2)*QDPD( I ) - QDPP(1) ) 
               QUAD%C( I+QUAD%N ) = QDPP(2)
            ENDDO
         ! bosonic Ozaki grid 
         ELSE IF ( QUAD%ITYPE == -202) THEN
            CALL SETUP_BOSONIC_GRID( QUAD%N, QUAD%C )  
         ! fermionic Ozaki grid 
         ELSE IF ( QUAD%ITYPE == -203) THEN
            CALL SETUP_FERMIONIC_GRID( QUAD%N, QUAD%C )
         ELSE
            CALL vtutor%bug("internal ERROR, SET_QUADRATURE_NON_OPTIMIZED reports: grid type unknown " &
               // str(QUAD%ITYPE), __FILE__, __LINE__)
         ENDIF        
        
         ! obtain maximum error 
         QUAD%INFINITYNORM = MAX( ABS( ERROR_FUNCTION( QUAD, QUAD%B+QDPD(0.0000001_q) ) ), &
                                  ABS( ERROR_FUNCTION( QUAD, QUAD%A ) ) )
         

         PROFILING_STOP( 'set_quadrature_non_optimized' )
      END SUBROUTINE SET_QUADRATURE_NON_OPTIMIZED

!*********************************************************************
!
! computes the bosonic grid 
!
!*********************************************************************

   SUBROUTINE SETUP_BOSONIC_GRID(NPOINTS, C)
     USE string, ONLY: str
     USE tutor, ONLY: vtutor
     INTEGER :: NPOINTS              ! total number of points wanted
     QDPREAL :: C(:)                 ! coefficients
     !local 
     INTEGER  :: M_DIM               ! dimension of matrix
     INTEGER  :: I, J                ! loop variables
     REAL(q), ALLOCATABLE :: D(:)   ! diagonal of tridiagonal matrix B
     REAL(q), ALLOCATABLE :: E(:)   ! sub-diagonal entries of B
     REAL(q), ALLOCATABLE :: Z(:,:) ! eigenvectors of tridiagonal matrix B
     INTEGER :: LDZ
     REAL(q), ALLOCATABLE :: WORK(:)! LAPACK working array
     INTEGER :: INFO                 ! LAPACK info variable
     REAL(q) :: AUX                 ! AUX
     REAL(q) :: TAUI                 ! AUX
     
     ! we need to approximate coth(x/2)+1/2-1/z by a continued fraction 
     ! it truns out this is 
     ! coth(x/2)+1/2-1/z = sum_i=1^NPOINTS R_i/(x-I P_i) + R_i/(x-I P_i)
     ! where P_i = 1/eig_i(B) and R_i = |<1|eigv(i)>|^2/(12 eig_i(B)^2)
     ! the matrix B is a symmetric tridiagonal matrix with zero diagonal 
     ! and 1/2{sqrt[(2n+1)(2n+3)] } as sub-diagonal
 
     M_DIM = 2*(NPOINTS-1)
     LDZ   = M_DIM
     ALLOCATE(D(M_DIM), E(M_DIM-1), Z(LDZ, M_DIM) )
     D=0
     E=0
     Z=0
     ALLOCATE(WORK(2*M_DIM-2))
     WORK=0

     !setup sub-diagonal
     DO I = 1, M_DIM-1
        E(I) = REAL(1,KIND=q)/(2*SQRT( REAL( (2*I+1)*(2*I+3 ) , KIND=q ) )) 
     ENDDO

     !call LAPACK 
     CALL DSTEV( 'V' , M_DIM, D, E, Z, LDZ, WORK, INFO )
     
     !stop in case something went wrong
     IF ( INFO /= 0 ) THEN
        CALL vtutor%error('ERROR during compution of bosonic frequencies, DSTEV reports:' // str(INFO))
     ENDIF
 
     !the first point is always 0 with residue 1 
     C(1) = QDPP(0)
     C(1+NPOINTS) = QDPP(1)

     !all others have R=2
     DO I = 1, M_DIM/2
        C( M_DIM/2-I+2 ) = QDPD( 1/D(M_DIM/2+I) )
        C( M_DIM/2-I+2 + NPOINTS ) = QDPD( 2*Z(1,M_DIM/2+I)**2/(12*D(M_DIM/2+I)**2) )
     ENDDO

     DEALLOCATE( WORK ) 
     DEALLOCATE( D, E, Z )

   ENDSUBROUTINE SETUP_BOSONIC_GRID

!*********************************************************************
!
! computes the fermionic grid 
!
!*********************************************************************

   SUBROUTINE SETUP_FERMIONIC_GRID(NPOINTS, C)
     USE string, ONLY: str
     USE tutor, ONLY: vtutor
     INTEGER :: NPOINTS              ! total number of points wanted
     QDPREAL :: C( : ) 
     !local 
     INTEGER  :: M_DIM               ! dimension of matrix
     INTEGER  :: I, J                ! loop variables
     REAL(q), ALLOCATABLE :: D(:)   ! diagonal of tridiagonal matrix B
     REAL(q), ALLOCATABLE :: E(:)   ! sub-diagonal entries of B
     REAL(q), ALLOCATABLE :: Z(:,:) ! eigenvectors of tridiagonal matrix B
     INTEGER :: LDZ
     REAL(q), ALLOCATABLE :: WORK(:)! LAPACK working array
     INTEGER :: INFO                 ! LAPACK info variable

     ! we need to approximate 1/2tanh(x/2) by a continued fraction 
     ! it truns out this is 
     ! 1/2tanh(x/2)-1/2 = sum_i=1^NPOINS R_i/(x-I P_i) + R_i/(x-I P_i)
     ! where P_i = 1/eig_i(B) and R_i = |<1|eigv(i)>|^2/(12 eig_i(B)^2)
     ! the matrix B is a symmetric tridiagonal matrix with zero diagonal 
     ! and 1/2{sqrt[(2n-1)(2n+1)] } as sub-diagonal
   
     M_DIM = 2*NPOINTS
     LDZ   = M_DIM
     ALLOCATE(D(M_DIM), E(M_DIM-1), Z(LDZ, M_DIM) )
     D=0
     E=0
     Z=0
 
     ALLOCATE(WORK(2*M_DIM-2))
     !setup sub-diagonal
     DO I = 1, M_DIM-1
        E(I) = REAL(1,KIND=q)/(2*SQRT( REAL( (2*I-1)*(2*I+1 ) , KIND=q ) )) 
     ENDDO

     WORK=0


     !call LAPACK 
     CALL DSTEV( 'V' , M_DIM, D, E, Z, LDZ, WORK, INFO )
     
     !stop in case something went wrong
     IF ( INFO /= 0 ) THEN
        CALL vtutor%error('ERROR during compution of bosonic frequencies, DSTEV reports:' // str(INFO))
     ENDIF
   
     !points appear in pairs with same residue, we store only the positive points
     DO I = 1, M_DIM/2
        C( M_DIM/2-I+1 ) = QDPD( 1/D(M_DIM/2+I) ) 
        C( M_DIM/2-I+1 + NPOINTS ) = 2*QDPD( Z(1,M_DIM/2+I)**2/(4*D(M_DIM/2+I)**2) )
     ENDDO

     DEALLOCATE( WORK ) 
     DEALLOCATE( D, E, Z )

   ENDSUBROUTINE SETUP_FERMIONIC_GRID


!**************************************************************************
!>sets inital quadrature_handle values based on requested interval length R
!>the routines used for initalization must be given in minimax_ini
!**************************************************************************
      SUBROUTINE INIT_COEFFICIENTS( QUADRATURE, R, IO ) 
         USE base
         USE prec
         USE minimax_ini
         USE tutor, ONLY: vtutor
         TYPE( quadrature_handle )  :: QUADRATURE ! quadrature
         QDPREAL                     :: R          ! requested interval boundary
         TYPE( in_struct )           :: IO
         ! local 
         INTEGER                     :: IERR
         QDPREAL                     :: RN, RM 
   
         RM = R 
         RN = R 
   
         ! quadrature must be initialized 
         IF ( QUADRATURE%N>0 ) THEN
            IF ( .NOT. ASSOCIATED( QUADRATURE%R1E6 ) ) THEN
               CALL vtutor%error('ERROR quadrature initialization routine R1E6 not associated')
            ENDIF
            ! find select RN of 1E6 interval COSH2
            CALL QUADRATURE%R1E6()
            CALL SET_RN_COEFFICIENTS( QUADRATURE, QUADRATURE%CTAB, RN, QUADRATURE%N, IERR )
            IF ( IERR ==-1 ) THEN
               CALL vtutor%error('ERROR, coefficients not found 1')
            ENDIF
            IF ( R > RN ) THEN
               RM=RN
               GOTO 10
            ENDIF
            ! RN holds current interval boundary 
   
            ! find select RM of 1E5 interval
            IF ( .NOT. ASSOCIATED( QUADRATURE%R1E5 ) ) THEN
               IF ( IO%IU0>=0 .AND. IO%NWRITE>2 ) &
               WRITE(*,*)' quadrature initialization routine R1E5 not found, will use R1E6'
               GOTO 10
            ENDIF
            CALL QUADRATURE%R1E5()
            CALL SET_RN_COEFFICIENTS( QUADRATURE, QUADRATURE%CTAB, RM, QUADRATURE%N, IERR )
            IF ( IERR ==-1 ) THEN
               CALL vtutor%error('ERROR, coefficients not found 2')
            ENDIF
            IF ( R > RM ) THEN
               CALL QUADRATURE%R1E5()
               CALL SET_RN_COEFFICIENTS( QUADRATURE, QUADRATURE%CTAB, RN, QUADRATURE%N, IERR )
               GOTO 10
            ELSE 
               RN = RM 
            ENDIF
   
            ! find select RM of 1E4 interval 
            IF ( .NOT. ASSOCIATED( QUADRATURE%R1E4 ) ) THEN
               IF ( IO%IU0>=0 .AND. IO%NWRITE>2 ) &
               WRITE(*,*)' quadrature initialization routine R1E5 not found, will use R1E5'
               GOTO 10
            ENDIF
            CALL QUADRATURE%R1E4()
            CALL SET_RN_COEFFICIENTS( QUADRATURE, QUADRATURE%CTAB, RM, QUADRATURE%N, IERR )
            IF ( IERR ==-1 ) THEN
               CALL vtutor%error('ERROR, coefficients not found 3')
            ENDIF
            IF ( R > RM ) THEN
               CALL QUADRATURE%R1E5()
               CALL SET_RN_COEFFICIENTS( QUADRATURE, QUADRATURE%CTAB, RN, QUADRATURE%N, IERR )
               GOTO 10
            ELSE 
               RN = RM 
            ENDIF
   
            ! find select RM of 1E4 interval
            IF ( .NOT. ASSOCIATED( QUADRATURE%R1E3 ) ) THEN
               IF ( IO%IU0>=0 .AND. IO%NWRITE>2 ) &
               WRITE(*,*)' quadrature initialization routine R1E3 not found, will use R1E5'
               GOTO 10
            ENDIF
            CALL QUADRATURE%R1E3()
            CALL SET_RN_COEFFICIENTS( QUADRATURE, QUADRATURE%CTAB, RM, QUADRATURE%N, IERR )
            IF ( IERR ==-1 ) THEN
               CALL vtutor%error('ERROR, coefficients not found 4')
            ENDIF
            IF ( R > RM ) THEN
               CALL QUADRATURE%R1E4()
               CALL SET_RN_COEFFICIENTS( QUADRATURE, QUADRATURE%CTAB, RN, QUADRATURE%N, IERR )
               GOTO 10
            ELSE 
               RN = RM 
            ENDIF
         ENDIF 
10    CONTINUE 
         QUADRATURE%B = RN 
         R = RN 
   
   
      END SUBROUTINE INIT_COEFFICIENTS

      SUBROUTINE SET_RN_COEFFICIENTS( QUADRATURE, CTAB, RN, N, IERR ) 
         TYPE( quadrature_handle ) :: QUADRATURE
         QDPREAL                    :: CTAB(:)
         QDPREAL                    :: RN    
         INTEGER                    :: N
         INTEGER                    :: IERR
         ! local                    
         INTEGER                    :: I, J 
   
         ! find out if requested order is available
         J = 0 
         DO I = 1 , N
            J = J + 2*I+1
         ENDDO
      
         IF ( SIZE( CTAB ) < J  ) THEN
            IERR = -1
         ELSE
            ! find position in array 
            J=1
            DO I = 1 , N-1
               J = J + 2*I+1
            ENDDO
            RN = CTAB(J)
            QUADRATURE%B = CTAB(J)
            
            DO I = 1 , 2*N
               QUADRATURE%C(I) = CTAB( J + I ) 
            ENDDO
            QUADRATURE%C(2*N+1) = ABS( ERROR_FUNCTION( QUADRATURE, QUADRATURE%B) )
            QUADRATURE%INFINITYNORM = ABS( QUADRATURE%C(2*N+1) ) 
            
            IERR = 1
         ENDIF 
      END SUBROUTINE SET_RN_COEFFICIENTS
   
      
   END SUBROUTINE ADJUST_BOUNDARY_FOR_CONVERGENCE

!****************************************************************************
!> helper: samples interval found in quadrature handle 
!****************************************************************************

   SUBROUTINE SETUP_SAMPLING_POINTS( QUADRATURE, T, Y, W, IV, IO )
      TYPE( quadrature_handle )             :: QUADRATURE
      QDPREAL                               :: T(:,:)
      QDPREAL                               :: Y(:)
      QDPREAL                               :: W(:)
      INTEGER                               :: IV
      TYPE( in_struct )                     :: IO
      !local 
      QDPREAL                               :: X1, X2
      INTEGER                               :: INU, NWRITE
      INTEGER I , J, L
      PROFILING_START( 'SETUP_SAMPLING_POINTS' )
     
      X1 = QUADRATURE%A
      X2 = QUADRATURE%B
      DO L =1, IV
         DO I =1, NDATALSQ
            T(I,L) = LOG(X2)+LOG(X2)*COS(((2*(2*NDATALSQ-I)+1)*PIQ)/(4*NDATALSQ))
            T(I,L) = EXP(T(I,L))
            T(I,L) = ( T(I,L)*( X2-X1 ) + ( X2*X1 - X2 ) )/( X2-QDPD(1) )
            W(I)   = QDPP(1)
            Y(I)   = QUADRATURE%F( T(I,L) )
            !print sample points only for IO%NWRITE >2
            IF ( IO%IU6>=0 .AND. IO%NWRITE > 2 ) THEN
               WRITE (IO%IU6,'(I5,2E16.7)') I, TOPRINT(T(I,L)), TOPRINT(Y(I))
            ENDIF
         ENDDO
      ENDDO
      PROFILING_STOP( 'SETUP_SAMPLING_POINTS' )
   END SUBROUTINE SETUP_SAMPLING_POINTS

!****************************************************************************
!> samples array T with NDATALSQ data points based on distribution in 
!! ALTERNANT
!****************************************************************************

   SUBROUTINE SAMPLING_FROM_ALTERNANT( QUADRATURE, T, Y, W, IV, IO )
      USE minimax_varpro, ONLY: PIKSRT2
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      TYPE( quadrature_handle )  :: QUADRATURE
      QDPREAL                    :: T(:,:)
      QDPREAL                    :: Y(:)
      QDPREAL                    :: W(:)
      INTEGER                    :: IV
      TYPE( in_struct )          :: IO
      !local 
      INTEGER                    :: I,K,L
      QDPREAL, ALLOCATABLE       :: XI(:)
 
      PROFILING_START( 'SAMPLING_FROM_ALTERNANT' ) 

      K = 2*QUADRATURE%N+1
      ALLOCATE( XI( K ) ) 
      XI(1:K) = QUADRATURE%X0( 1:K )
      T(1:K,IV) = XI(1:K)
      doubling: DO L = 1, 1000
         ! add points inbetween
         DO I = 1, 2*QUADRATURE%N-1
            IF ( I+K >NDATALSQ ) EXIT doubling
            T( I + K,IV ) = ( XI(I+1) + XI(I) ) / QDPP(2)
         ENDDO
         K = K + 2*QUADRATURE%N-1
         IF ( ALLOCATED( XI ) ) DEALLOCATE( XI ) 
         ALLOCATE( XI( K ) ) 
         XI( 1:K ) = T( 1:K,IV )
         CALL PIKSRT2( K, XI ) 
      ENDDO doubling

      CALL PIKSRT2( NDATALSQ, T(1:NDATALSQ,IV) ) 

      DO I = 1, NDATALSQ
         ! scale alternant to a smaller interval (2N extremum)
         ! smaller multiplication factor decreases interval even faster
         ! but is unstable 
         T(I,IV) = T(I,IV) * QDPP(1.00)*QUADRATURE%B/T(NDATALSQ,IV)
         Y(I) = QUADRATURE%F( T(I,IV) )
         W(I) = 1
      ENDDO
      IF ( L > 1000 ) THEN
         CALL vtutor%error("error in SAMPLING_FROM_ALTERNANT " // str(L))
      ENDIF

      PROFILING_STOP( 'SAMPLING_FROM_ALTERNANT' ) 
   END SUBROUTINE SAMPLING_FROM_ALTERNANT

!****************************************************************************
!> helper: determines alternant of error function
!****************************************************************************
!#define debug
   SUBROUTINE DETERMINE_ALTERNANT( QUADRATURE, IO )
      USE minimax_varpro, ONLY: PIKSRT2
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      TYPE( quadrature_handle ) :: QUADRATURE
      TYPE( in_struct )         :: IO
      !local 
      QDPREAL, ALLOCATABLE      :: X_I(:)
      QDPREAL, ALLOCATABLE      :: X_S(:)
      QDPREAL, ALLOCATABLE      :: X_0(:)
      QDPREAL, ALLOCATABLE      :: X0(:)
      QDPREAL, ALLOCATABLE      :: X_M(:)
      QDPREAL                   :: Y1,Y2
      INTEGER                   :: IDBL
      INTEGER                   :: IDBL_MAX
      INTEGER                   :: I,J
      INTEGER                   :: NALT
      QDPREAL                   :: Z 
      INTEGER                   :: NI, NS
      INTEGER,ALLOCATABLE       :: INFO(:)
      INTEGER                   :: NODE_ME = 1  
      INTEGER                   :: IND

      PROFILING_START( 'DETERMINE_ALTERNANT' )
      IF ( QUADRATURE%N < 1 ) THEN
         CALL vtutor%error("Wrong use of DETERMINE_ALTERNANT, quadrature_handle seems not to be set &
            &" // str(QUADRATURE%N))
      ENDIF
#ifdef MPI 
      NODE_ME = QUADRATURE%COMM%NODE_ME
#endif
      NALT = 0 
      ! allocate storage for alternant
      ALLOCATE ( X0( 2*QUADRATURE%N+1 ) )
      DO I = 1, 2*QUADRATURE%N+1
         X0(I) = 0
      ENDDO

      ! intialize intersection points 
      NI = QUADRATURE%N
      ALLOCATE ( X_I( NI ) )
      DO I = 1, NI
         X_I( I ) = QUADRATURE%B * &
            QUADRATURE%C( I ) / QUADRATURE%C( QUADRATURE%N )
      ENDDO

      Z = ERROR_FUNCTION(QUADRATURE,QUADRATURE%A)
      ! add left side of interval to alternant
      NALT = 0 
      IF ( ABS( Z ) > QDPP(1E-32) ) THEN
         X0( 1 ) = QUADRATURE%A
         NALT = 1
      ENDIF
      ! increase sampling points 
      IDBL_MAX = 10
      double_points: DO IDBL = 1, IDBL_MAX

         ALLOCATE( X_S( SIZE( X_I ) ) )
         ALLOCATE( X_0( SIZE( X_I ) ) )
         ALLOCATE( INFO( SIZE(X_0) ) )
         INFO=0
          
         ! create sampling points based on intersection points
         X_S( 1 ) = QUADRATURE%A + X_I( 1 ) / QDPP(2)
         DO I = 2, SIZE( X_I ) 
            X_S( I ) = QUADRATURE%A + (X_I( I ) + X_I( I-1 )) / QDPP(2)
         ENDDO

         !look for extremal points  
         root_search: DO I = 1, SIZE( X_S ) 
            X_0( I ) = QDPP( 0 )
#ifdef MPI
            ! block this loop with MPI 
            IND = LOCAL_INDEX_CYCLIC( QUADRATURE%COMM%NCPU, SIZE(X_S), &
               QUADRATURE%COMM%NODE_ME, I )
            IF (IND < 0 ) CYCLE root_search
#endif
            CALL RTSAFE( QUADRATURE, X_0( I ), X_S( I ), INFO(I) )
         ENDDO root_search

#ifdef MPI
         ! synchronization
         CALLMPI( M_sum_qd( QUADRATURE%COMM, X_0(1), SIZE( X_S ) ) )
         CALLMPI( M_sum_i( QUADRATURE%COMM, INFO(1), SIZE( X_S ) ) )
#endif

         DO I = 1, SIZE( X_S ) 
            IF ( INFO(I) == 0 ) CYCLE 
            ! compare with points in alternant array 
            ! might be first point
            IF ( NALT == 0 ) THEN
               NALT = NALT + 1 
               X0( NALT ) = X_0( I ) 
            ELSE
               ! points must differ
               compare:DO J = 1, NALT 
                  IF ( ABS( X0( J ) - X_0( I ) ) < QDPP(1E-8) ) THEN
                     EXIT compare
                  ENDIF
               ENDDO compare
  
               ! in case no double entries found, add to alternant
               IF ( J > NALT .AND. NALT < 2*QUADRATURE%N+1 ) THEN
                  CALL INJECT_X0( QUADRATURE%A, X0, X_0(I), NALT )
                  IF ( NALT == 2*QUADRATURE%N + 1 ) THEN
                     EXIT double_points
                  ENDIF
               ENDIF           
            ENDIF
         ENDDO

         ! check if all points have been found 
         IF ( NALT == 2*QUADRATURE%N ) THEN
            Y1 = ERROR_FUNCTION( QUADRATURE, X0( NALT ) ) 
            Y2 = ERROR_FUNCTION( QUADRATURE, QUADRATURE%B )
            ! add right border if signs differ
            IF ( SIGN_QD( Y1 ) /= SIGN_QD( Y2 ) ) THEN
               NALT = NALT + 1 
               X0( NALT ) = QUADRATURE%B 
               EXIT double_points
            ENDIF

         ! maybe all points have been found 
         ELSE IF ( NALT == 2*QUADRATURE%N + 1 ) THEN
            EXIT double_points
         ENDIF


         ! replace intersection points by current sampling points
         ! double sampling points but avoid previously tested sampling points
         ! in this way the interval is sampled more elegantly and efficiently
         ALLOCATE( X_M( 2*SIZE( X_I ) ) )
         DO I = 1, SIZE( X_I )
            X_M( 2*I ) = X_I( I ) 
            X_M( 2*I - 1 ) =  X_S( I ) - QUADRATURE%A
         ENDDO
         IF ( ALLOCATED(X_I) ) DEALLOCATE( X_I )
         ALLOCATE( X_I( SIZE( X_M ) ) )
         DO I = 1, SIZE( X_M )
            X_I( I ) = X_M( I ) 
         ENDDO
         DEALLOCATE( X_S )
         DEALLOCATE( X_0 )
         DEALLOCATE( X_M )
         DEALLOCATE( INFO )

      ENDDO double_points
      IF ( ALLOCATED( X_I ) ) DEALLOCATE( X_I )
      IF ( ALLOCATED( X_0 ) ) DEALLOCATE( X_0 )
      IF ( ALLOCATED( X_M ) ) DEALLOCATE( X_M )

      ! update alternant 
      DO I = 1, 2*QUADRATURE%N+1
         QUADRATURE%X0( I ) = X0( I ) 
#ifdef debug
IF ( I>2*QUADRATURE%N-1 )&
         WRITE(*,*)I,X0(I),ERROR_FUNCTION(QUADRATURE, X0(I) )
#endif
      ENDDO

      IF ( IDBL > IDBL_MAX ) THEN
         CALL vtutor%bug("internal ERROR, DETERMINE_ALTERNANT was not able to find alternant " // &
            str(NALT) // " " // str(SIZE(X_S)) // " " // str(IDBL_MAX), __FILE__, __LINE__)
      ENDIF 
      IF ( ALLOCATED( X_S ) ) DEALLOCATE( X_S )

      PROFILING_STOP( 'DETERMINE_ALTERNANT' )
      CONTAINS
   !*****************************************************************
   ! simple root finder for nonlinear error function
   !*****************************************************************
      SUBROUTINE RTSAFE( QUADRATURE, ROOT, GUESS, INFO )
         TYPE( quadrature_handle ) :: QUADRATURE
         QDPREAL                   :: ROOT
         QDPREAL                   :: GUESS
         INTEGER                   :: INFO
         ! local 
         INTEGER                   :: I 
         QDPREAL                   :: Y,D1Y, D2Y
         QDPREAL                   :: Z, DZ
         INTEGER, PARAMETER        :: MAXNEWT = 200
#ifdef qd_emulate 
#define XACC      qdreal('1.E-12')
#define EPS       qdreal('1.E-32')
#else
         QDPREAL, PARAMETER         :: XACC=1.E-12_qd  
         QDPREAL, PARAMETER         :: EPS =1.E-32_qd  
#endif
         PROFILING_START('RTSAFE')
         INFO=0
         Z=GUESS

         DO I=1,MAXNEWT
            D1Y = D1_ERROR_FUNCTION( QUADRATURE, Z )
            D2Y = D2_ERROR_FUNCTION( QUADRATURE, Z )
            DZ = D1Y / D2Y
            Z = Z - DZ
            IF ( ABS( DZ ) < XACC ) THEN
               ! must be in interval and non-zero function value
               Y = ERROR_FUNCTION( QUADRATURE, Z  )
               IF ( Z > QUADRATURE%A .AND. Z < QUADRATURE%B .AND. ABS( Y ) > EPS )THEN
                  INFO=I
               ENDIF
               ROOT = Z
               PROFILING_STOP('RTSAFE')
               RETURN
            ENDIF
         ENDDO
         
         INFO=0
         ROOT = Z
         PROFILING_STOP('RTSAFE')
      END SUBROUTINE RTSAFE

   !*****************************************************************
   ! adds alternant point to X0 array, sorts 
   !*****************************************************************
!#define debug2
      SUBROUTINE INJECT_X0( A, X0, X, NALT )
         USE string, ONLY: str
         USE tutor, ONLY: vtutor
         QDPREAL                    :: A
         QDPREAL                    :: X0(:)
         QDPREAL                    :: X
         INTEGER                    :: NALT
         ! local 
         QDPREAL,ALLOCATABLE        :: X1(:)
         INTEGER                    :: I, J, K

         PROFILING_START('INJECT_X0')
#ifdef debug2
   WRITE(*,*)'before',X,ALLOCATED(X1)
   DO I = 1, NALT 
     WRITE(*,*)I,X0(I)
   ENDDO
#endif

         inject: DO I = 1, NALT 
            IF ( X0( I ) < X ) THEN
               IF ( I < NALT ) THEN
                  IF (  X0( I + 1 ) > X0( I ) )THEN
                     IF ( X0( I + 1 ) > X ) THEN
                        ALLOCATE( X1( NALT-I ) )
                        DO J = I+1, NALT
                           X1( J-I ) = X0 ( J ) 
                        ENDDO
                        
                        X0( I+1 ) = X 

                        DO J = I+1, NALT
                           X0( J + 1 ) = X1( J-I ) 
                        ENDDO
                        DEALLOCATE( X1 )

                        EXIT inject

                     ENDIF
                  ENDIF
               ELSE 
                  X0( I+1 ) = X 
                  EXIT inject
               ENDIF
            ELSE IF ( X0( I ) > X .AND. I == 1 ) THEN
               ALLOCATE( X1( NALT ) )
               DO J = 1 , NALT
                  X1( J ) = X0 ( J ) 
               ENDDO
               X0( 1 ) = X

               DO J = 2 , NALT+1
                  X0( J ) = X1(J-1)
               ENDDO
               DEALLOCATE( X1 )
               EXIT inject
            ENDIF
         ENDDO inject 

         IF ( I > NALT ) THEN
            CALL vtutor%error("Error INJECT could not inject X")
         ENDIF

         NALT = NALT + 1

#ifdef debug2
   WRITE(*,*)'After', ALLOCATED( X1 )
   DO I = 1, NALT 
     WRITE(*,*)I,X0(I)
   ENDDO
   WRITE(*,*)
#endif
         PROFILING_STOP('INJECT_X0')
      END SUBROUTINE INJECT_X0
!#undef debug2

   END SUBROUTINE DETERMINE_ALTERNANT
#undef debug

!********************************************************************
!> stores scaled quadrature grid points to grid handle
!********************************************************************
   SUBROUTINE SCALE_QUADRATURES( GRIDS, IO ) 
      TYPE( imag_grid_handle) :: GRIDS
      TYPE( in_struct )       :: IO

      ! release quadrature for real part of bosonic functions
      IF ( GRIDS%FREQ_BOS_RE%N > 0 ) THEN
         CALL SCALE_QUADRATURE( GRIDS%FREQ_BOS_RE, GRIDS%BOS_RE, &
            GRIDS%BOS_RE_WEIGHT, GRIDS%BOS_RE_ERROR, IO ) 
      ENDIF  
      ! release quadrature for imaginary part of bosonic functions
      IF ( GRIDS%FREQ_BOS_IM%N > 0 ) THEN
         CALL SCALE_QUADRATURE( GRIDS%FREQ_BOS_IM, GRIDS%BOS_IM, &
            GRIDS%BOS_IM_WEIGHT, GRIDS%BOS_IM_ERROR, IO ) 
      ENDIF  
      ! release quadrature for re part of fermionic functions
      IF ( GRIDS%FREQ_FER_RE%N > 0 ) THEN
         CALL SCALE_QUADRATURE( GRIDS%FREQ_FER_RE, GRIDS%FER_RE, &
            GRIDS%FER_RE_WEIGHT, GRIDS%FER_RE_ERROR, IO ) 
      ENDIF  
      ! release quadrature for re part of fermionic functions
      IF ( GRIDS%FREQ_FER_IM%N > 0 ) THEN
         CALL SCALE_QUADRATURE( GRIDS%FREQ_FER_IM, GRIDS%FER_IM, &
            GRIDS%FER_IM_WEIGHT, GRIDS%FER_IM_ERROR, IO ) 
      ENDIF  

      ! release time quadrature 
      IF ( GRIDS%TIME%N > 0 ) THEN
         CALL SCALE_QUADRATURE( GRIDS%TIME, GRIDS%TAU, &
            GRIDS%TAU_WEIGHT, GRIDS%TAU_ERROR, IO ) 
         ! scale transformation matrix from time to BOS_RE grid
         IF ( GRIDS%FREQ_BOS_RE%N > 0 ) THEN
            CALL SCALE_TRANSFORMATION( GRIDS%FREQ_BOS_RE, GRIDS%BOS_RE, &
               GRIDS%TO_BOS_RE, GRIDS%TO_BOS_RE_ERROR, IO )
            CALL SCALE_TRANSFORMATION( GRIDS%FREQ_BOS_RE, GRIDS%BOS_RE, &
               GRIDS%TO_BOS_RE_CONJG, GRIDS%TO_BOS_RE_CONJG_ERROR, IO )
         ENDIF 
         ! scale transformation matrix from time to BOS_IM grid
         IF ( GRIDS%FREQ_BOS_IM%N > 0 ) THEN
            CALL SCALE_TRANSFORMATION( GRIDS%FREQ_BOS_IM, GRIDS%BOS_IM, &
               GRIDS%TO_BOS_IM, GRIDS%TO_BOS_IM_ERROR, IO )
            CALL SCALE_TRANSFORMATION( GRIDS%FREQ_BOS_IM, GRIDS%BOS_IM, &
               GRIDS%TO_BOS_IM_CONJG, GRIDS%TO_BOS_IM_CONJG_ERROR, IO )
         ENDIF 
         ! scale transformation matrix from time to FER_RE grid
         IF ( GRIDS%FREQ_FER_RE%N > 0 ) THEN
            CALL SCALE_TRANSFORMATION( GRIDS%FREQ_FER_RE, GRIDS%FER_RE, &
               GRIDS%TO_FER_RE, GRIDS%TO_FER_RE_ERROR, IO )
            CALL SCALE_TRANSFORMATION( GRIDS%FREQ_FER_RE, GRIDS%FER_RE, &
               GRIDS%TO_FER_RE_CONJG, GRIDS%TO_FER_RE_CONJG_ERROR, IO )
         ENDIF 
         ! scale transformation matrix from time to FER_IM grid
         IF ( GRIDS%FREQ_FER_IM%N > 0 ) THEN
            CALL SCALE_TRANSFORMATION( GRIDS%FREQ_FER_IM, GRIDS%FER_IM, &
               GRIDS%TO_FER_IM, GRIDS%TO_FER_IM_ERROR, IO )
            CALL SCALE_TRANSFORMATION( GRIDS%FREQ_FER_IM, GRIDS%FER_IM, &
               GRIDS%TO_FER_IM_CONJG, GRIDS%TO_FER_IM_CONJG_ERROR, IO )
         ENDIF 

         ! scale transformation matrix from time to FER_IM grid
         IF ( ASSOCIATED( GRIDS%TO_TAU0 )  ) THEN
            CALL SCALE_TRANSFORMATION( GRIDS%TIME, GRIDS%TAU, &
               GRIDS%TO_TAU0, GRIDS%TO_TAU0_ERROR, IO )
            CALL SCALE_TRANSFORMATION( GRIDS%TIME, GRIDS%TAU, &
               GRIDS%TO_TAU0_CONJG, GRIDS%TO_TAU0_CONJG_ERROR, IO )
         ENDIF 
      ENDIF  

      CONTAINS

      !***************************************************************
      !>copies quadrature handle to proper grid handle arrays
      !***************************************************************
      SUBROUTINE SCALE_QUADRATURE( QUAD, ABSZISSA, WEIGHT, GRID_ERROR, IO ) 
         TYPE( quadrature_handle ) :: QUAD
         REAL(q), POINTER          :: ABSZISSA(:)
         REAL(q), POINTER          :: WEIGHT(:)
         REAL(q)                   :: GRID_ERROR 
         TYPE( in_struct )         :: IO 
         !local 
         LOGICAL                   :: LWRITE=.FALSE.
         INTEGER                   :: I
         INTEGER                   :: M
         REAL(q)                   :: WEIGHT_SCALE=1._q

         ! number of grid points 
         M = QUAD%N
         ! artificially add zero point to grid 
         IF ( QUAD%ITYPE == -2 .OR. QUAD%ITYPE == -20 ) M = M+1 

         ! allocate abszissa array 
         NULLIFY( ABSZISSA ) 
         IF ( .NOT. ASSOCIATED( ABSZISSA ) ) ALLOCATE( ABSZISSA( M ) )
         ! allocate weight array 
         NULLIFY( WEIGHT ) 
         IF ( .NOT. ASSOCIATED( WEIGHT ) ) ALLOCATE( WEIGHT( M ) )

         ! add zero point
         IF ( M>QUAD%N ) THEN
            ABSZISSA( 1 ) = 0 
            WEIGHT( 1 ) = QDPP(1)*QUAD%SCALING
            DO I = 1, QUAD%N 
               ABSZISSA( I+1 ) = QUAD%C( I )*QUAD%SCALING
               WEIGHT( I+1 ) = QUAD%C( QUAD%N + I )*QUAD%SCALING
            ENDDO
         ELSE
            DO I = 1, QUAD%N 
               ABSZISSA( I ) = QUAD%C( I )*QUAD%SCALING
               WEIGHT( I ) = QUAD%C( QUAD%N + I )*QUAD%SCALING
            ENDDO
         ENDIF
         ! copy abszissas and weights, this weired loop should be able 
         ! to take into account if first point is zero point
 
         IF ( QUAD%ITYPE == 1 ) THEN
            WEIGHT_SCALE=0.5_q
         ELSE IF ( QUAD%ITYPE == 2 .OR.  QUAD%ITYPE == 3 .OR. QUAD%ITYPE == 20 ) THEN 
            WEIGHT_SCALE=0.25_q
         ELSE IF ( QUAD%ITYPE ==-1 ) THEN 
            WEIGHT_SCALE=2._q
         ELSE
            WEIGHT_SCALE=1._q
         ENDIF 

         DO I = 1, M
            WEIGHT( I ) = WEIGHT( I )*WEIGHT_SCALE 
         ENDDO

         ! scale quadrature error 
         ! error scaling is differs in frequency domain
         IF ( ABS( QUAD%ITYPE ) == 1 .OR. QUAD%ITYPE==-101) THEN
            GRID_ERROR = QUAD%INFINITYNORM * QUAD%SCALING
         ELSE
            GRID_ERROR = QUAD%INFINITYNORM / QUAD%SCALING
         ENDIF
         
   
         IF ( ABS( QUAD%ITYPE ) == 1 .OR. MOD( ABS( QUAD%ITYPE ), 2 ) == 0 ) THEN
            LWRITE = .TRUE. 
         ENDIF
         ! dump full grid for NWRITE>2
         IF ( IO%NWRITE > 2 ) THEN
             LWRITE = .TRUE. 
         ENDIF
         ! dump result to file
         IF( IO%IU6>=0 .AND. LWRITE ) THEN
            IF ( ABS( QUAD%ITYPE ) == 2 ) THEN
               WRITE( IO%IU6, '(A,E11.4)' ) ' bosonic grid (re):',GRID_ERROR
            ELSE IF ( QUAD%ITYPE == -3 ) THEN
               WRITE( IO%IU6, '(A,E11.4)' ) ' bosonic grid (im):',GRID_ERROR
            ELSE IF ( ABS( QUAD%ITYPE ) == 4 .OR. QUAD%ITYPE == -40 .OR.  QUAD%ITYPE == 20) THEN
               WRITE( IO%IU6, '(A,E11.4)' ) ' fermionic grid (re):',GRID_ERROR
            ELSE IF ( ABS( QUAD%ITYPE ) == 5 .OR. QUAD%ITYPE == 3 ) THEN
               WRITE( IO%IU6, '(A,E11.4)' ) ' fermionic grid (im):',GRID_ERROR
            ELSE IF ( ABS( QUAD%ITYPE ) == 1 ) THEN
               WRITE( IO%IU6, '(A,E11.4)' ) ' time grid:',GRID_ERROR
            ELSE IF ( QUAD%ITYPE == -101 ) THEN
               WRITE( IO%IU6, '(A)' ) ' uniform time grid:'
            ELSE IF ( QUAD%ITYPE == -102 ) THEN
               WRITE( IO%IU6, '(A)' ) ' bosonic Matsubara grid:'
            ELSE IF ( QUAD%ITYPE == -103 ) THEN
               WRITE( IO%IU6, '(A)' ) ' fermionic Matsubara grid:'
            ELSE IF ( QUAD%ITYPE == -202 ) THEN
               WRITE( IO%IU6, '(A)' ) ' hypergeometric bosonic Matsubara grid:'
            ELSE IF ( QUAD%ITYPE == -203 ) THEN
               WRITE( IO%IU6, '(A)' ) ' hypergeometric fermionic Matsubara grid:'
            ELSE
               WRITE( IO%IU6, '(I4,A)' ) QUAD%ITYPE,'-grid:'
            ENDIF
            WRITE( IO%IU6, '(5F13.4)' )ABSZISSA
         ENDIF 
         ! to std out
         IF( IO%IU0>=0 .AND. LWRITE )  THEN
            IF ( ABS( QUAD%ITYPE ) == 2 ) THEN
               WRITE( *, '(A,E11.4)' ) ' bosonic grid (re):',GRID_ERROR
            ELSE IF ( QUAD%ITYPE  == -3 ) THEN
               WRITE( *, '(A,E11.4)' ) ' bosonic grid (im):',GRID_ERROR
            ELSE IF ( ABS( QUAD%ITYPE ) == 4 .OR. QUAD%ITYPE == -40 .OR.  QUAD%ITYPE == 20 ) THEN
               WRITE( *, '(A,E11.4)' ) ' fermionic grid (re):',GRID_ERROR
            ELSE IF ( ABS( QUAD%ITYPE ) == 5 .OR. QUAD%ITYPE == 3 ) THEN
               WRITE( *, '(A,E11.4)' ) ' fermionic grid (im):',GRID_ERROR
            ELSE IF ( ABS( QUAD%ITYPE ) == 1 ) THEN
               WRITE( *, '(A,E11.4)' ) ' time grid:',GRID_ERROR
            ELSE IF ( QUAD%ITYPE == -101 ) THEN
               WRITE( *, '(A)' ) ' uniform time grid:'
            ELSE IF ( QUAD%ITYPE == -102 ) THEN
               WRITE( *, '(A)' ) ' bosonic Matsubara grid:'
            ELSE IF ( QUAD%ITYPE == -103 ) THEN
               WRITE( *, '(A)' ) ' fermionic Matsubara grid:'
            ELSE IF ( QUAD%ITYPE == -202 ) THEN
               WRITE( *, '(A)' ) ' hypergeomtric bosonic Matsubara grid:'
            ELSE IF ( QUAD%ITYPE == -203 ) THEN
               WRITE( *, '(A)' ) ' hypergeometric fermionic Matsubara grid:'
            ELSE
               WRITE( *, '(I4,A)' ) QUAD%ITYPE,'-grid:'
            ENDIF
            IF ( QUAD%N >  32 ) THEN 
               WRITE( *, '(4F13.4,A)' )ABSZISSA(1:4),'  ... '
               WRITE( *, '(A,4F13.4)' )'       ...   ',ABSZISSA(QUAD%N-3:QUAD%N)
            ELSE
               WRITE( *, '(5F13.4)' )ABSZISSA
            ENDIF
         ENDIF

         ! verbose mode 
         IF ( IO%NWRITE > 2 ) THEN
            IF( IO%IU6>=0 ) WRITE( IO%IU6, '(A,E11.4)' ) &
               '      abszissa       weights  scaling=',QUAD%SCALING
            DO I = 1, M
               IF( IO%IU6>=0 ) WRITE( IO%IU6, '(A,2E15.8)' )'     ',ABSZISSA(I),WEIGHT(I)
            ENDDO
         ENDIF

      END SUBROUTINE SCALE_QUADRATURE

      !***************************************************************
      !>scales transformation matrix 
      !***************************************************************
      SUBROUTINE SCALE_TRANSFORMATION( QUAD, POINTS, T, T_ERROR, IO) 
         TYPE( quadrature_handle ) :: QUAD
         REAL(q)                   :: POINTS(:)
         REAL(q)                   :: T(:,:)
         REAL(q)                   :: T_ERROR(:)
         TYPE( in_struct )         :: IO 
         !local 
         INTEGER                   :: I, J 
         REAL(q)                   :: EMAX
         REAL(q)                   :: WEIGHT_SCALE
         PROFILING_START( 'scale_transformation' )
    
         IF ( ABS(QUAD%ITYPE) == 2 ) THEN 
            WEIGHT_SCALE=2._q
         ! transformation for Fermions
         ELSE IF ( ABS(QUAD%ITYPE) == 103 ) THEN 
            WEIGHT_SCALE=0.5_q
         ELSE
            WEIGHT_SCALE=1._q
         ENDIF 
         
         ! scale transformation matrix 
         DO J = 1, SIZE( T, 2 ) 
            DO I = 1, SIZE( T, 1 ) 
               T(I,J) = T(I,J)/QUAD%SCALING*WEIGHT_SCALE
            ENDDO
         ENDDO

         IF( IO%IU6 >= 0 ) THEN
            IF ( ABS( QUAD%ITYPE ) == 2 ) THEN
               WRITE( IO%IU6, 20 ) 
            ELSE IF ( ABS( QUAD%ITYPE ) == 3 ) THEN
               WRITE( IO%IU6, 30 )
            ELSE IF ( ABS( QUAD%ITYPE ) == 4 .OR. QUAD%ITYPE == -40 .OR.  QUAD%ITYPE == 20) THEN
               WRITE( IO%IU6, 40 )
            ELSE IF ( ABS( QUAD%ITYPE ) == 5 ) THEN
               WRITE( IO%IU6, 50 )
            ELSE IF ( QUAD%ITYPE == -102 ) THEN
               WRITE( IO%IU6, 60 )
            ELSE IF ( QUAD%ITYPE == -103 ) THEN
               WRITE( IO%IU6, 70 )
            ELSE
               WRITE( IO%IU6, 10 ) QUAD%ITYPE
            ENDIF
            ! also dump transformation matrix 
            IF( IO%NWRITE > 2 ) THEN
               WRITE( IO%IU6, '(A)' )'  Transformation matrix (column major)'
               WRITE( IO%IU6, '(8E12.4)' ) T
            ENDIF 
         ENDIF

         IF( IO%IU6>=0 ) WRITE( IO%IU6, '(A,I3,A)' ) &
            '      point          trafo. error'
         
         EMAX = -1.E6_q
         ! scale errors
         DO I = 1, SIZE( T_ERROR )
            IF ( ABS( QUAD%ITYPE ) == 1 ) THEN
               T_ERROR(I) = T_ERROR(I)*QUAD%SCALING
            ELSE
               T_ERROR(I) = T_ERROR(I)/QUAD%SCALING
            ENDIF
            IF( IO%IU6>=0 ) WRITE( IO%IU6, '(A,2E15.8)' )'     ',POINTS(I),T_ERROR(I)
            EMAX = MAX( EMAX, T_ERROR( I ) )
         ENDDO

         ! dump result to stdou
         IF( IO%IU0>=0 .ANd. IO%NWRITE > 2 ) THEN
            IF ( ABS( QUAD%ITYPE ) == 2 ) THEN
               WRITE( *, 2 ) EMAX 
            ELSE IF ( ABS( QUAD%ITYPE ) == 3 ) THEN
               WRITE( *, 3 ) EMAX 
            ELSE IF ( ABS( QUAD%ITYPE ) == 4 .OR. QUAD%ITYPE == -40 .OR.  QUAD%ITYPE == 20) THEN
               WRITE( *, 4 ) EMAX 
            ELSE IF ( ABS( QUAD%ITYPE ) == 5 ) THEN
               WRITE( *, 5 ) EMAX 
            ELSE IF ( QUAD%ITYPE < -100) THEN
               ! dont write for Matsubara grid 
            ELSE
               WRITE( *, 1 ) QUAD%ITYPE, EMAX 
            ENDIF
         ENDIF 

         ! number of grid points 
         PROFILING_STOP( 'scale_transformation' )
1        FORMAT( ' Maximum transformation error to ',I4,' grid: ', E11.4 )
2        FORMAT( ' Maximum transformation error to bosonic (re) grid: ', E11.4 )
3        FORMAT( ' Maximum transformation error to bosonic (im) grid: ', E11.4 )
4        FORMAT( ' Maximum transformation error to fermion (re) grid: ', E11.4 )
5        FORMAT( ' Maximum transformation error to fermion (im) grid: ', E11.4 )
10       FORMAT( '  Transformation errors for each grid point of type', I3 )
20       FORMAT( '  Transformation errors for each bosonic (re) grid point' )
30       FORMAT( '  Transformation errors for each bosonic (im) grid point' )
40       FORMAT( '  Transformation errors for each fermionic (re) grid point' )
50       FORMAT( '  Transformation errors for each fermionic (im) grid point' )
60       FORMAT( '  Transformation errors for each bosonic Matsubara grid point' )
70       FORMAT( '  Transformation errors for each fermionic Matsubara grid point' )
      END SUBROUTINE SCALE_TRANSFORMATION

   END SUBROUTINE SCALE_QUADRATURES

!********************************************************************
!>frees all allocated quadratures in imaginary grid handle
!********************************************************************
   SUBROUTINE RELEASE_QUADRATURES( GRIDS, IO ) 
      TYPE( imag_grid_handle) :: GRIDS
      TYPE( in_struct )       :: IO
      PROFILING_START( 'RELEASE_QUADRATURES' )

      ! release time quadrature 
      IF ( GRIDS%TIME%N > 0 ) THEN
         CALL KILL_QUADRATURE( GRIDS%TIME ) 
      ENDIF  
      ! release quadrature for real part of bosonic functions
      IF ( GRIDS%FREQ_BOS_RE%N > 0 ) THEN
         CALL KILL_QUADRATURE( GRIDS%FREQ_BOS_RE ) 
      ENDIF  
      ! release quadrature for imaginary part of bosonic functions
      IF ( GRIDS%FREQ_BOS_IM%N > 0 ) THEN
         CALL KILL_QUADRATURE( GRIDS%FREQ_BOS_IM ) 
      ENDIF  
      ! release quadrature for re part of fermionic functions
      IF ( GRIDS%FREQ_FER_RE%N > 0 ) THEN
         CALL KILL_QUADRATURE( GRIDS%FREQ_FER_RE ) 
      ENDIF  
      ! release quadrature for re part of fermionic functions
      IF ( GRIDS%FREQ_FER_IM%N > 0 ) THEN
         CALL KILL_QUADRATURE( GRIDS%FREQ_FER_IM ) 
      ENDIF  

      PROFILING_STOP( 'RELEASE_QUADRATURES' )

   END SUBROUTINE RELEASE_QUADRATURES

   !***************************************************************
   !>destroys a quadrature handle QUAD
   !***************************************************************
   SUBROUTINE KILL_QUADRATURE( QUAD ) 
      TYPE( quadrature_handle ) :: QUAD
      PROFILING_START( 'KILL_QUADRATURE' )
       
      ! reset quadrature identifier
#ifdef CleanQuads
      QUAD%ITYPE = 0 
      ! reset quadrature order 
      QUAD%N = 0 
      ! reset interval boundaries
      QUAD%A = QDPP( 0 )
      QUAD%B = QDPP( 0 )
      ! reset scaling
      QUAD%SCALING = QDPP( 0 )
      ! reset infinity norm 
      QUAD%INFINITYNORM = QDPP( 0 )
#endif 

      ! release coefficients
      IF( ASSOCIATED(QUAD%C) ) DEALLOCATE( QUAD%C )
      ! release alternant
      IF( ASSOCIATED(QUAD%X0) ) DEALLOCATE( QUAD%X0 )
      ! release basis function
      IF( ASSOCIATED(QUAD%PHI) ) NULLIFY( QUAD%PHI )
      ! release dual basis
      IF( ASSOCIATED(QUAD%PSI) ) NULLIFY( QUAD%PSI )
      ! release basis function for quadrature error
      IF( ASSOCIATED(QUAD%PHI2) ) NULLIFY( QUAD%PHI2 )
      ! release first derivative of basis function for quadrature error
      IF( ASSOCIATED(QUAD%D_PHI2_DZ) ) NULLIFY( QUAD%D_PHI2_DZ )
      ! release first derivative of basis function for quadrature error
      IF( ASSOCIATED(QUAD%D_PHI2_DL) ) NULLIFY( QUAD%D_PHI2_DL )
      ! release first derivative of basis function for quadrature error
      IF( ASSOCIATED(QUAD%D2_PHI2_DZ2) ) NULLIFY( QUAD%D2_PHI2_DZ2 )
      ! release object function of quadrature 
      IF( ASSOCIATED(QUAD%F) ) NULLIFY( QUAD%F )
      ! release first derivative of object function of quadrature 
      IF( ASSOCIATED(QUAD%D_F_DZ) ) NULLIFY( QUAD%D_F_DZ )
      ! release second derivative of object function of quadrature 
      IF( ASSOCIATED(QUAD%D2_F_DZ2) ) NULLIFY( QUAD%D2_F_DZ2 )
      ! release auxilary function for low frequency expansion 
      IF( ASSOCIATED(QUAD%TRANS_TAYLOR) ) NULLIFY( QUAD%TRANS_TAYLOR )
      ! release initialization function for starting coefficients for R 1E3
      IF( ASSOCIATED(QUAD%R1E3) ) NULLIFY( QUAD%R1E3 )
      ! release initialization function for starting coefficients for R 1E4
      IF( ASSOCIATED(QUAD%R1E4) ) NULLIFY( QUAD%R1E4 )
      ! release initialization function for starting coefficients for R 1E5
      IF( ASSOCIATED(QUAD%R1E5) ) NULLIFY( QUAD%R1E5 )
      ! release initialization function for starting coefficients for R 1E6
      IF( ASSOCIATED(QUAD%R1E6) ) NULLIFY( QUAD%R1E6 )
      ! release initialized coefficients
      IF( ASSOCIATED(QUAD%CTAB) ) NULLIFY( QUAD%CTAB )

      IF( ASSOCIATED(QUAD%PHI_CONJG) ) NULLIFY( QUAD%PHI_CONJG )
      IF( ASSOCIATED(QUAD%PSI_CONJG) ) NULLIFY( QUAD%PSI_CONJG )
      IF( ASSOCIATED(QUAD%TRANS_TAYLOR_CONJG) ) NULLIFY( QUAD%TRANS_TAYLOR_CONJG )
      IF( ASSOCIATED(QUAD%FITCOF) ) NULLIFY( QUAD%FITCOF )

      PROFILING_STOP( 'KILL_QUADRATURE' )
   END SUBROUTINE KILL_QUADRATURE

!********************************************************************
!>adds one time and frequency point to GRIDS 
!********************************************************************

  SUBROUTINE ADD_POINT_IMAG_GRID_HANDLE(GRIDS )
     USE string, ONLY: str
     USE tutor, ONLY: vtutor
     TYPE(imag_grid_handle) :: GRIDS     !handle 
     TYPE(imag_grid_handle) :: GRIDS_NEW !handle 
   ! local
     INTEGER                :: NO        !number of frequency points

     PROFILING_START('add_point_imag_grid_handle')
     IF ( GRIDS%LFINITE_TEMPERATURE ) THEN
        CALL vtutor%error("adding point for T>0 not implemented yet, sorry")
     ENDIF

     !set number of grid points
     NO=GRIDS%NOMEGA
      
     CALL  ALLOCATE_IMAG_GRID_HANDLE_AS(GRIDS_NEW, GRIDS, NO+1)

     IF( ASSOCIATED( GRIDS%TAU ) ) THEN
        GRIDS_NEW%TAU(1:NO)=GRIDS%TAU(1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%TAU_WEIGHT ) ) THEN
        GRIDS_NEW%TAU_WEIGHT(1:NO)=GRIDS%TAU_WEIGHT(1:NO)
     ENDIF

     IF( ASSOCIATED( GRIDS%BOS_RE ) ) THEN
        GRIDS_NEW%BOS_RE(1:NO)=GRIDS%BOS_RE(1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%BOS_RE_WEIGHT ) ) THEN
        GRIDS_NEW%BOS_RE_WEIGHT(1:NO)=GRIDS%BOS_RE_WEIGHT(1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%TO_BOS_RE ) ) THEN
        GRIDS_NEW%TO_BOS_RE(1:NO,1:NO)=GRIDS%TO_BOS_RE(1:NO,1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%TO_BOS_RE_ERROR ) ) THEN
        GRIDS_NEW%TO_BOS_RE_ERROR(1:NO)=GRIDS%TO_BOS_RE_ERROR(1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%TO_BOS_RE_CONJG ) ) THEN
        GRIDS_NEW%TO_BOS_RE_CONJG(1:NO,1:NO)=GRIDS%TO_BOS_RE_CONJG(1:NO,1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%TO_BOS_RE_CONJG_ERROR ) ) THEN
        GRIDS_NEW%TO_BOS_RE_CONJG_ERROR(1:NO)=GRIDS%TO_BOS_RE_CONJG_ERROR(1:NO)
     ENDIF

     IF( ASSOCIATED( GRIDS%BOS_IM ) ) THEN
        GRIDS_NEW%BOS_IM(1:NO)=GRIDS%BOS_IM(1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%BOS_IM_WEIGHT ) ) THEN
        GRIDS_NEW%BOS_IM_WEIGHT(1:NO)=GRIDS%BOS_IM_WEIGHT(1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%TO_BOS_IM ) ) THEN
        GRIDS_NEW%TO_BOS_IM(1:NO,1:NO)=GRIDS%TO_BOS_IM(1:NO,1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%TO_BOS_IM_ERROR ) ) THEN
        GRIDS_NEW%TO_BOS_IM_ERROR(1:NO)=GRIDS%TO_BOS_IM_ERROR(1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%TO_BOS_IM_CONJG ) ) THEN
        GRIDS_NEW%TO_BOS_IM_CONJG(1:NO,1:NO)=GRIDS%TO_BOS_IM_CONJG(1:NO,1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%TO_BOS_IM_CONJG_ERROR ) ) THEN
        GRIDS_NEW%TO_BOS_IM_CONJG_ERROR(1:NO)=GRIDS%TO_BOS_IM_CONJG_ERROR(1:NO)
     ENDIF

     IF( ASSOCIATED( GRIDS%FER_RE ) ) THEN
        GRIDS_NEW%FER_RE(1:NO)=GRIDS%FER_RE(1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%FER_RE_WEIGHT ) ) THEN
        GRIDS_NEW%FER_RE_WEIGHT(1:NO)=GRIDS%FER_RE_WEIGHT(1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%TO_FER_RE ) ) THEN
        GRIDS_NEW%TO_FER_RE(1:NO,1:NO)=GRIDS%TO_FER_RE(1:NO,1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%TO_FER_RE_ERROR ) ) THEN
        GRIDS_NEW%TO_FER_RE_ERROR(1:NO)=GRIDS%TO_FER_RE_ERROR(1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%TO_FER_RE_CONJG ) ) THEN
        GRIDS_NEW%TO_FER_RE_CONJG(1:NO,1:NO)=GRIDS%TO_FER_RE_CONJG(1:NO,1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%TO_FER_RE_CONJG_ERROR ) ) THEN
        GRIDS_NEW%TO_FER_RE_CONJG_ERROR(1:NO)=GRIDS%TO_FER_RE_CONJG_ERROR(1:NO)
     ENDIF

     IF( ASSOCIATED( GRIDS%FER_IM ) ) THEN
        GRIDS_NEW%FER_IM(1:NO)=GRIDS%FER_IM(1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%FER_IM_WEIGHT ) ) THEN
        GRIDS_NEW%FER_IM_WEIGHT(1:NO)=GRIDS%FER_IM_WEIGHT(1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%TO_FER_IM ) ) THEN
        GRIDS_NEW%TO_FER_IM(1:NO,1:NO)=GRIDS%TO_FER_IM(1:NO,1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%TO_FER_IM_ERROR ) ) THEN
        GRIDS_NEW%TO_FER_IM_ERROR(1:NO)=GRIDS%TO_FER_IM_ERROR(1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%TO_FER_IM_CONJG ) ) THEN
        GRIDS_NEW%TO_FER_IM_CONJG(1:NO,1:NO)=GRIDS%TO_FER_IM_CONJG(1:NO,1:NO)
     ENDIF
     IF( ASSOCIATED( GRIDS%TO_FER_IM_CONJG_ERROR ) ) THEN
        GRIDS_NEW%TO_FER_IM_CONJG_ERROR(1:NO)=GRIDS%TO_FER_IM_CONJG_ERROR(1:NO)
     ENDIF

     ! safe call of this routine is done before distribution
     IF ( ASSOCIATED( GRIDS%T%DISTRIBUTION ) .OR. &
          ASSOCIATED( GRIDS%B%DISTRIBUTION ) .OR. &
          ASSOCIATED( GRIDS%F%DISTRIBUTION ) ) THEN
        CALL vtutor%bug("internal ERROR in VASP: ADD_POINT_IMAG_GRID_HANDLE cannot be called after " &
           // "distribution of grid points", __FILE__, __LINE__)
     ENDIF

     ! deallocate old GRID
     CALL RELEASE_IMAG_GRID_HANDLE(GRIDS)

     ! copy over remaining members
     GRIDS=GRIDS_NEW
     
     PROFILING_STOP('add_point_imag_grid_handle')

     CONTAINS
  
     SUBROUTINE ALLOCATE_IMAG_GRID_HANDLE_AS(GRIDS, GOLD, NOMEGA)
        use mpimy, ONLY : communic
        use base, ONLY : in_struct
        TYPE(imag_grid_handle)  :: GRIDS     !handle 
        TYPE(imag_grid_handle)  :: GOLD      ! old handle to check whats allocated
        INTEGER                 :: NOMEGA    !number of frequency points
        INTEGER                 :: I 

        PROFILING_START('ALLOCATE_IMAG_GRID_HANDLE_AS')

        NULLIFY( GRIDS%TAU ) 
        NULLIFY( GRIDS%TAU_WEIGHT ) 
        IF ( ASSOCIATED( GOLD%TAU ) ) THEN
           ALLOCATE( GRIDS%TAU( NOMEGA ) ) 
           GRIDS%TAU=0
        ENDIF
        IF ( ASSOCIATED( GOLD%TAU_WEIGHT ) ) THEN
           ALLOCATE( GRIDS%TAU_WEIGHT( NOMEGA ) ) 
           GRIDS%TAU_WEIGHT=0
        ENDIF

        NULLIFY( GRIDS%BOS_RE ) 
        NULLIFY( GRIDS%BOS_RE_WEIGHT ) 
        NULLIFY( GRIDS%TO_BOS_RE ) 
        NULLIFY( GRIDS%TO_BOS_RE_ERROR ) 
        NULLIFY( GRIDS%TO_BOS_RE_CONJG ) 
        NULLIFY( GRIDS%TO_BOS_RE_CONJG_ERROR ) 
        IF ( ASSOCIATED( GOLD%BOS_RE ) ) THEN
           ALLOCATE( GRIDS%BOS_RE( NOMEGA ) ) 
           GRIDS%BOS_RE=0
        ENDIF
        IF ( ASSOCIATED( GOLD%BOS_RE_WEIGHT ) ) THEN
           ALLOCATE( GRIDS%BOS_RE_WEIGHT( NOMEGA ) ) 
           GRIDS%BOS_RE_WEIGHT=0
        ENDIF
        IF ( ASSOCIATED( GOLD%TO_BOS_RE ) ) THEN
           ALLOCATE( GRIDS%TO_BOS_RE( NOMEGA, NOMEGA ) ) 
           GRIDS%TO_BOS_RE=0
        ENDIF
        IF ( ASSOCIATED( GOLD%TO_BOS_RE_ERROR) ) THEN
           ALLOCATE( GRIDS%TO_BOS_RE_ERROR( NOMEGA ) ) 
           GRIDS%TO_BOS_RE_ERROR=0
        ENDIF
        IF ( ASSOCIATED( GOLD%TO_BOS_RE_CONJG ) ) THEN
           ALLOCATE( GRIDS%TO_BOS_RE_CONJG( NOMEGA, NOMEGA ) ) 
           GRIDS%TO_BOS_RE_CONJG=0
        ENDIF
        IF ( ASSOCIATED( GOLD%TO_BOS_RE_CONJG_ERROR) ) THEN
           ALLOCATE( GRIDS%TO_BOS_RE_CONJG_ERROR( NOMEGA ) ) 
           GRIDS%TO_BOS_RE_CONJG_ERROR=0
        ENDIF

        NULLIFY( GRIDS%BOS_IM ) 
        NULLIFY( GRIDS%BOS_IM_WEIGHT ) 
        NULLIFY( GRIDS%TO_BOS_IM ) 
        NULLIFY( GRIDS%TO_BOS_IM_ERROR ) 
        NULLIFY( GRIDS%TO_BOS_IM_CONJG ) 
        NULLIFY( GRIDS%TO_BOS_IM_CONJG_ERROR ) 
        IF ( ASSOCIATED( GOLD%BOS_IM ) ) THEN
           ALLOCATE( GRIDS%BOS_IM( NOMEGA ) ) 
           GRIDS%BOS_IM=0
        ENDIF
        IF ( ASSOCIATED( GOLD%BOS_IM_WEIGHT ) ) THEN
           ALLOCATE( GRIDS%BOS_IM_WEIGHT( NOMEGA ) ) 
           GRIDS%BOS_IM_WEIGHT=0
        ENDIF
        IF ( ASSOCIATED( GOLD%TO_BOS_IM ) ) THEN
           ALLOCATE( GRIDS%TO_BOS_IM( NOMEGA, NOMEGA ) ) 
           GRIDS%TO_BOS_IM=0
        ENDIF
        IF ( ASSOCIATED( GOLD%TO_BOS_IM_ERROR) ) THEN
           ALLOCATE( GRIDS%TO_BOS_IM_ERROR( NOMEGA ) ) 
           GRIDS%TO_BOS_IM_ERROR=0
        ENDIF
        IF ( ASSOCIATED( GOLD%TO_BOS_IM_CONJG ) ) THEN
           ALLOCATE( GRIDS%TO_BOS_IM_CONJG( NOMEGA, NOMEGA ) ) 
           GRIDS%TO_BOS_IM_CONJG=0
        ENDIF
        IF ( ASSOCIATED( GOLD%TO_BOS_IM_CONJG_ERROR) ) THEN
           ALLOCATE( GRIDS%TO_BOS_IM_CONJG_ERROR( NOMEGA ) ) 
           GRIDS%TO_BOS_IM_CONJG_ERROR=0
        ENDIF

        NULLIFY( GRIDS%FER_RE ) 
        NULLIFY( GRIDS%FER_RE_WEIGHT ) 
        NULLIFY( GRIDS%TO_FER_RE ) 
        NULLIFY( GRIDS%TO_FER_RE_ERROR ) 
        NULLIFY( GRIDS%TO_FER_RE_CONJG ) 
        NULLIFY( GRIDS%TO_FER_RE_CONJG_ERROR ) 
        IF ( ASSOCIATED( GOLD%FER_RE ) ) THEN
           ALLOCATE( GRIDS%FER_RE( NOMEGA ) ) 
           GRIDS%FER_RE=0
        ENDIF
        IF ( ASSOCIATED( GOLD%FER_RE_WEIGHT ) ) THEN
           ALLOCATE( GRIDS%FER_RE_WEIGHT( NOMEGA ) ) 
           GRIDS%FER_RE_WEIGHT=0
        ENDIF
        IF ( ASSOCIATED( GOLD%TO_FER_RE ) ) THEN
           ALLOCATE( GRIDS%TO_FER_RE( NOMEGA, NOMEGA ) ) 
           GRIDS%TO_FER_RE=0
        ENDIF
        IF ( ASSOCIATED( GOLD%TO_FER_RE_ERROR) ) THEN
           ALLOCATE( GRIDS%TO_FER_RE_ERROR( NOMEGA ) ) 
           GRIDS%TO_FER_RE_ERROR=0
        ENDIF
        IF ( ASSOCIATED( GOLD%TO_FER_RE_CONJG ) ) THEN
           ALLOCATE( GRIDS%TO_FER_RE_CONJG( NOMEGA, NOMEGA ) ) 
           GRIDS%TO_FER_RE_CONJG=0
        ENDIF
        IF ( ASSOCIATED( GOLD%TO_FER_RE_CONJG_ERROR) ) THEN
           ALLOCATE( GRIDS%TO_FER_RE_CONJG_ERROR( NOMEGA ) ) 
           GRIDS%TO_FER_RE_CONJG_ERROR=0
        ENDIF

        NULLIFY( GRIDS%FER_IM ) 
        NULLIFY( GRIDS%FER_IM_WEIGHT ) 
        NULLIFY( GRIDS%TO_FER_IM ) 
        NULLIFY( GRIDS%TO_FER_IM_ERROR ) 
        NULLIFY( GRIDS%TO_FER_IM_CONJG ) 
        NULLIFY( GRIDS%TO_FER_IM_CONJG_ERROR ) 
        IF ( ASSOCIATED( GOLD%FER_IM ) ) THEN
           ALLOCATE( GRIDS%FER_IM( NOMEGA ) ) 
           GRIDS%FER_IM=0
        ENDIF
        IF ( ASSOCIATED( GOLD%FER_IM_WEIGHT ) ) THEN
           ALLOCATE( GRIDS%FER_IM_WEIGHT( NOMEGA ) ) 
           GRIDS%FER_IM_WEIGHT=0
        ENDIF
        IF ( ASSOCIATED( GOLD%TO_FER_IM ) ) THEN
           ALLOCATE( GRIDS%TO_FER_IM( NOMEGA, NOMEGA ) ) 
           GRIDS%TO_FER_IM=0
        ENDIF
        IF ( ASSOCIATED( GOLD%TO_FER_IM_ERROR) ) THEN
           ALLOCATE( GRIDS%TO_FER_IM_ERROR( NOMEGA ) ) 
           GRIDS%TO_FER_IM_ERROR=0
        ENDIF
        IF ( ASSOCIATED( GOLD%TO_FER_IM_CONJG ) ) THEN
           ALLOCATE( GRIDS%TO_FER_IM_CONJG( NOMEGA, NOMEGA ) ) 
           GRIDS%TO_FER_IM_CONJG=0
        ENDIF
        IF ( ASSOCIATED( GOLD%TO_FER_IM_CONJG_ERROR) ) THEN
           ALLOCATE( GRIDS%TO_FER_IM_CONJG_ERROR( NOMEGA ) ) 
           GRIDS%TO_FER_IM_CONJG_ERROR=0
        ENDIF

        PROFILING_STOP('ALLOCATE_IMAG_GRID_HANDLE_AS')
     END SUBROUTINE ALLOCATE_IMAG_GRID_HANDLE_AS

      SUBROUTINE RELEASE_IMAG_GRID_HANDLE(GRIDS)
         TYPE( imag_grid_handle ) :: GRIDS
         PROFILING_START( 'RELEASE_IMAG_GRID_HANDLE')
 
         ! deallocate associated time points
         IF ( ASSOCIATED( GRIDS%TAU ) ) DEALLOCATE( GRIDS%TAU )
         IF ( ASSOCIATED( GRIDS%TAU_WEIGHT ) ) DEALLOCATE( GRIDS%TAU_WEIGHT )
 
         ! deallocate associated bosonic frequencies for real part 
         IF ( ASSOCIATED( GRIDS%BOS_RE ) ) DEALLOCATE( GRIDS%BOS_RE ) 
         IF ( ASSOCIATED( GRIDS%BOS_RE_WEIGHT ) ) DEALLOCATE( GRIDS%BOS_RE_WEIGHT ) 
         IF ( ASSOCIATED( GRIDS%TO_BOS_RE ) ) DEALLOCATE( GRIDS%TO_BOS_RE ) 
         IF ( ASSOCIATED( GRIDS%TO_BOS_RE_ERROR ) ) DEALLOCATE( GRIDS%TO_BOS_RE_ERROR ) 
         IF ( ASSOCIATED( GRIDS%TO_BOS_RE_CONJG ) ) DEALLOCATE( GRIDS%TO_BOS_RE_CONJG ) 
         IF ( ASSOCIATED( GRIDS%TO_BOS_RE_CONJG_ERROR ) ) DEALLOCATE( GRIDS%TO_BOS_RE_CONJG_ERROR ) 
         
         ! deallocate associated bosonic frequencies for imaginary part 
         IF ( ASSOCIATED( GRIDS%BOS_IM ) ) DEALLOCATE( GRIDS%BOS_IM ) 
         IF ( ASSOCIATED( GRIDS%BOS_IM_WEIGHT ) ) DEALLOCATE( GRIDS%BOS_IM_WEIGHT ) 
         IF ( ASSOCIATED( GRIDS%TO_BOS_IM ) ) DEALLOCATE( GRIDS%TO_BOS_IM ) 
         IF ( ASSOCIATED( GRIDS%TO_BOS_IM_ERROR ) ) DEALLOCATE( GRIDS%TO_BOS_IM_ERROR ) 
         IF ( ASSOCIATED( GRIDS%TO_BOS_IM_CONJG ) ) DEALLOCATE( GRIDS%TO_BOS_IM_CONJG ) 
         IF ( ASSOCIATED( GRIDS%TO_BOS_IM_CONJG_ERROR ) ) DEALLOCATE( GRIDS%TO_BOS_IM_CONJG_ERROR ) 
 
         ! deallocate associated fermionic frequencies for real part 
         IF ( ASSOCIATED( GRIDS%FER_RE ) ) DEALLOCATE( GRIDS%FER_RE ) 
         IF ( ASSOCIATED( GRIDS%FER_RE_WEIGHT ) ) DEALLOCATE( GRIDS%FER_RE_WEIGHT ) 
         IF ( ASSOCIATED( GRIDS%TO_FER_RE ) ) DEALLOCATE( GRIDS%TO_FER_RE ) 
         IF ( ASSOCIATED( GRIDS%TO_FER_RE_ERROR ) ) DEALLOCATE( GRIDS%TO_FER_RE_ERROR ) 
         IF ( ASSOCIATED( GRIDS%TO_FER_RE_CONJG ) ) DEALLOCATE( GRIDS%TO_FER_RE_CONJG ) 
         IF ( ASSOCIATED( GRIDS%TO_FER_RE_CONJG_ERROR ) ) DEALLOCATE( GRIDS%TO_FER_RE_CONJG_ERROR ) 
 
         ! deallocate associated fermionic frequencies for imaginary part 
         IF ( ASSOCIATED( GRIDS%FER_IM ) ) DEALLOCATE( GRIDS%FER_IM ) 
         IF ( ASSOCIATED( GRIDS%FER_IM_WEIGHT ) ) DEALLOCATE( GRIDS%FER_IM_WEIGHT ) 
         IF ( ASSOCIATED( GRIDS%TO_FER_IM ) ) DEALLOCATE( GRIDS%TO_FER_IM ) 
         IF ( ASSOCIATED( GRIDS%TO_FER_IM_ERROR ) ) DEALLOCATE( GRIDS%TO_FER_IM_ERROR ) 
         IF ( ASSOCIATED( GRIDS%TO_FER_IM_CONJG ) ) DEALLOCATE( GRIDS%TO_FER_IM_CONJG ) 
         IF ( ASSOCIATED( GRIDS%TO_FER_IM_CONJG_ERROR ) ) DEALLOCATE( GRIDS%TO_FER_IM_CONJG_ERROR ) 
 
         PROFILING_STOP( 'RELEASE_IMAG_GRID_HANDLE')
       
      END SUBROUTINE RELEASE_IMAG_GRID_HANDLE  
   
   END SUBROUTINE ADD_POINT_IMAG_GRID_HANDLE

!**********************************************************************
!>  
!>   DESCRIPTION:
!    Solves following non linear system 
!>   \f$ 
!         f(x_i)-\sum_{j=1}^N \lambda_{j+N}\phi(x_i,\lambda_j)^2
!           +(-1)^{i-1} \lambda_{2N+1} = 0
!>   \f$ 
!
!**********************************************************************
!#define DumpJacobian
   SUBROUTINE SOLVE_NONLS( QUADRATURE, ITER, IO )
      USE prec
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      TYPE( quadrature_handle ) :: QUADRATURE 
      INTEGER, INTENT(INOUT)    :: ITER    !iterations needed to solution
      TYPE( in_struct )         :: IO
      !local
      INTEGER                   :: NL !number of variables
      QDPREAL                   :: J(2*QUADRATURE%N+1,2*QUADRATURE%N+1)  !Jacobian
      QDPREAL                   :: Y(2*QUADRATURE%N+1)                   !r.h.s. of system of equations
      INTEGER                   :: I,K 
      INTEGER, PARAMETER        :: MAXNEWTON = 5000
#ifdef qd_emulate 
#define DAMPNEWT      qdreal('0.35')
#define SUPERTINY     qdreal('1.E-33')
#else
      QDPREAL, PARAMETER   :: DAMPNEWT=0.35_qd  
      QDPREAL, PARAMETER   :: SUPERTINY=1.E-33_qd  
#endif
      PROFILING_START('SOLVE_NONLS')

      DO ITER=1,MAXNEWTON
#ifdef DumpJacobian
         IF( ITER < 3 ) THEN
            CALL PRINT_ERROR_FUNCTION( QUADRATURE, IO, IT = ITER+100 )
            CALL PRINT_QUADRATURE( QUADRATURE, IO, IT=ITER )
         ENDIF
#endif 
         NL = 2*QUADRATURE%N+1
         ! nonlinear system is solved iteratively using the Newton-Raphson 
         ! method. For this Jacobian w.r.t. the parameter C is needed
         DO I=1,NL
            !l.h.s.
            DO K=1,2*QUADRATURE%N
               J(I,K)=GRAD_ERROR_FUNCTION(QUADRATURE, QUADRATURE%X0( I ), K)
            ENDDO
            !r.h.s.
            IF ( MOD( I , 2 ) == 0 ) THEN
               J(I,NL)= QDPP(1)
               Y(I)=-ERROR_FUNCTION( QUADRATURE, QUADRATURE%X0( I ))-QUADRATURE%C(NL)
            ELSE
               J(I,NL)=-QDPP(1)
               Y(I)=-ERROR_FUNCTION( QUADRATURE, QUADRATURE%X0( I ))+QUADRATURE%C(NL)
            ENDIF
         ENDDO
#ifdef DumpJacobian
         IF ( ITER < 3 ) THEN
            DO I = 1, NL 
                WRITE(100+ITER,'(17F12.6)')J(I,:)
            ENDDO
            WRITE(100+ITER,'(17F12.6)')Y(:)
         ENDIF
#endif
         ! call linear solver for quadruple precision 
         ! on exit Y is solution vector X of the system J.X = Y 
         CALL SOLVE_LINEAR_SYSTEM_QD( NL, J, Y )  

         ! add damped solution vector to previous vector
         DO I=1, NL
            QUADRATURE%C(I)=QUADRATURE%C(I)+DAMPNEWT*Y(I)
            ! weights must be positive for a proper minimax solution
            IF ( (QUADRATURE%C(I)<0) .AND. (I>QUADRATURE%N .AND. I/=NL) ) THEN 
               CALL vtutor%error(str(ITER) // "Error in SOLVE_NONLS, weight is negative: " // str(I) &
                  // " " // str(QUADRATURE%C(I)) // " " // str(DAMPNEWT*Y(I)))
            ENDIF
         ENDDO  
          
         ! solution found
         !IF ( DOT_PRODUCT_MPR(Y,Y) < VERYTINY )  EXIT 
         !IF ( DOT_PRODUCT_MPR(Y,Y) < VERYTINY .OR.  Y(NL) < VERYTINY)  EXIT 
         IF ( DOT_PRODUCT_MPR(Y,Y) < VERYTINY .OR.  ABS(Y(NL)) < SUPERTINY)  EXIT 

!IF ( ITER < 200 ) THEN
!WRITE(100,'(I6,2E40.30)')ITER, DOT_PRODUCT_MPR(Y,Y)
!WRITE(100,'(I6,2E40.30)')( I, Y(I),Y(I+QUADRATURE%N), I=1,QUADRATURE%N )
!WRITE(100,'(I6,E40.30)')2*QUADRATURE%N+1, Y( 2*QUADRATURE%N+1)
!ELSE
!STOP
!ENDIF

      ENDDO 
!STOP

      IF ( ITER > MAXNEWTON ) THEN
         CALL vtutor%error("Error in SOLVE_NONLS: Newton-Raphson not converged " // str(MAXNEWTON))
      ENDIF 

      
      PROFILING_STOP('SOLVE_NONLS')

   END SUBROUTINE SOLVE_NONLS
#undef DumpJacobian


!**********************************************************************
!>  
!>   DESCRIPTION:
!>   Solves the linear system 
!>   \f$ 
!>       J \cdot \vec{x} = \vec{y}
!>   \f$ 
!>
!>   array J(:,:) is destroyed on exit
!>   
!>   @param[in]  J  (coefficient matrix of dimension NL x NL )
!>   @param[in]  Y  (r.h.s. vector \f$\vec{y}\f$ )
!>   @param[out] Y  (solution \f$\vec{x}\f$)
!>   
!**********************************************************************
   SUBROUTINE SOLVE_LINEAR_SYSTEM_QD( NL, J, Y  ) 
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      INTEGER              :: NL
      QDPREAL              :: J(:,:)
      QDPREAL              :: Y(:)
      !local
      QDPREAL              :: D             
      INTEGER              :: IPIV(NL)       ! for QDLUDCMP
      INTEGER              :: INFO           ! info code

      PROFILING_START( 'SOLVE_LINEAR_SYSTEM_QD' ) 
      ! LU decomposition
      CALL QDLUDCMP( J, NL, NL, IPIV, D, INFO )
      IF (INFO /= 0 ) THEN
         CALL vtutor%error("ERROR, SOLVE_LINEAR_SYSTEM_QD failed with code: " // str(INFO))
      ENDIF 
      !solve linear equation 
      CALL QDLUBKSB( J, NL, NL, IPIV, Y)

      PROFILING_STOP( 'SOLVE_LINEAR_SYSTEM_QD' ) 

   END SUBROUTINE SOLVE_LINEAR_SYSTEM_QD

!**********************************************************************
!
! LU decomposition for quadruple precision  (NR)
!
! Given a matrix a(1:n,1:n) , with physical dimension 
! np by np , this routine replaces it by the 
! LU decomposition of a rowwise permutation of itself. 
! a and n are input. a is output, arranged as in 
! equation (2.3.14) of NR above; indx(1:n) is an output vector 
! that records the row permutation effected by the partial 
! pivoting; d is output as ±1 depending on whether the number 
! of row interchanges was even or odd, respectively. 
! This routine is used in combination with qdlubksb to solve 
! linear equations
!
!**********************************************************************
   SUBROUTINE QDLUDCMP(A,N,NP,INDX,D, INFO)
      IMPLICIT NONE
      QDPREAL  :: A(NP,NP)
      INTEGER  :: N
      INTEGER  :: NP
      INTEGER  :: INDX(N)
      QDPREAL  :: D
      INTEGER  :: INFO 
      !local 
      INTEGER, PARAMETER  :: NMAX=500
      INTEGER I,IMAX,J,K
      QDPREAL  AAMAX,DUM,SUM,VV(NMAX) 

      PROFILING_START('QDLUDCMP')
 
      !vv stores the implicit scaling of each row.
      D=1
        
      !no row interchanges yet.
      !loop over rows to get the implicit scaling information.
      DO I=1,N
         AAMAX=0
         DO J=1,N
            IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
         ENDDO
         IF ( AAMAX .EQ. 0 ) THEN
            !singular matrix in qdludcmp no nonzero largest element.
            INFO=-I 
            PROFILING_STOP('QDLUDCMP')
            RETURN
         ENDIF  
         !save the scaling.
         VV(I)=1/AAMAX
      ENDDO

      !this is the loop over columns of crout’s method.
      DO J=1,N
         !this is equation (2.3.12) except for i = j.
         DO I=1,J-1
            SUM=A(I,J)
            DO K=1,I-1
               SUM=SUM-A(I,K)*A(K,J)
            ENDDO
            A(I,J)=SUM
         ENDDO
         !initialize for the search for largest pivot element.
         AAMAX=0

         !this is i = j of equation (2.3.12) and i = j + 1 . . . n
         !of equation (2.3.13).
         DO I=J,N
            SUM=A(I,J)
            DO K=1,J-1
               SUM=SUM-A(I,K)*A(K,J)
            ENDDO
            A(I,J)=SUM
            !figure of merit for the pivot.
            DUM=VV(I)*ABS(SUM)
            !is it better than the best so far?
            IF (DUM.GE.AAMAX) THEN
               IMAX=I
               AAMAX=DUM
            ENDIF
         ENDDO

         !do we need to interchange rows?
         IF (J.NE.IMAX)THEN
            !YES, DO SO...
            DO K=1,N
               DUM=A(IMAX,K)
               A(IMAX,K)=A(J,K)
               A(J,K)=DUM
            ENDDO
            !...and change the parity of d.
            D=-D
            !also interchange the scale factor.
            VV(IMAX)=VV(J)
         ENDIF

         INDX(J)=IMAX
         !if the pivot element is zero the matrix is singular (at least to the precision of the al-
         !gorithm). for some applications on singular matrices, it is desirable to substitute tiny
         !for zero.
         IF( A(J,J) .EQ. 0 ) A(J,J)=VERYTINY

         !now, finally, divide by the pivot element.
         IF(J.NE.N)THEN
            DUM=1/A(J,J)
            DO I=J+1,N
                A(I,J)=A(I,J)*DUM
            ENDDO
         ENDIF
      ENDDO
 
      INFO=0
      PROFILING_STOP('QDLUDCMP')
   END SUBROUTINE QDLUDCMP

!**********************************************************************
! quadruple precision of DGETRS
!
! Solves the set of n linear equations A · X = B. Here A is 
! input, not as the matrix A but rather as its LU 
! decomposition, determined by the routine ludcmp. indx is 
! input as the permutation vector returned by ludcmp. 
! b(1:n) is input as the right-hand side vector B, and 
! returns with the solution vector X. a , n , np, and 
! indx are not modified by this routine and can be left 
! in place for successive calls with different right-hand 
! sides b. This routine takes into account the possibility 
! that b will begin with many zero elements, so it is 
! efficient for use in matrix inversion.
!**********************************************************************
   SUBROUTINE QDLUBKSB(A,N,NP,INDX,B)
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: NP
      INTEGER, INTENT(IN) :: INDX(N)
      QDPREAL  A(NP,NP)
      QDPREAL  B(N)
      !local 
      INTEGER I,II,J,LL
      QDPREAL  SUM
   
      PROFILING_START('QDLUBKSB')
   
      II=0
      !when ii is set to a positive value, it will become the in-
      !dex of the first nonvanishing element of b. we now do
      !the forward substitution, equation (2.3.6). the only new
      !wrinkle is to unscramble the permutation as we go.
      DO I=1,N
         LL=INDX(I)
         SUM=B(LL)
         B(LL)=B(I)
         IF (II.NE.0)THEN
            DO J=II,I-1
               SUM=SUM-A(I,J)*B(J)
            ENDDO
         !a nonzero element was encountered, so from now on we will
         !have to do the sums in the loop above.
         ELSEIF (SUM.NE.0) THEN
            II=I
         ENDIF
         B(I)=SUM
      ENDDO
   
      !now we do the backsubstitution, equation (2.3.7).
      DO I=N,1,-1
         SUM=B(I)
         DO J=I+1,N
            SUM=SUM-A(I,J)*B(J)
         ENDDO 
         !store a component of the solution vector x.
         B(I)=SUM/A(I,I)
      ENDDO
      !all done!
   
      PROFILING_STOP('QDLUBKSB')
   END SUBROUTINE QDLUBKSB


!*******************************************************************
!>computes local index of a cyclic distribution of an array with NDA 
!!distributed among NCPU nodes
!*******************************************************************
      
   FUNCTION LOCAL_INDEX_CYCLIC( NCPU, NDAT, NODE_ME, M )
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INTEGER  :: NODE_ME
      INTEGER  :: M, L
      INTEGER  :: NCPU
      INTEGER  :: NDAT
      INTEGER  :: IND
      INTEGER  :: LOCAL_INDEX_CYCLIC
      INTEGER  :: IBLOCK
#ifndef MPI
      LOCAL_INDEX_CYCLIC = M
      RETURN
#else
      IF ( NCPU == 1 ) THEN
         LOCAL_INDEX_CYCLIC = M
         RETURN
      ENDIF

      ! number of ranks have to be at least 1 
      IF ( NCPU < 1  ) THEN
         CALL vtutor%bug('internal ERROR in LOCAL_INDEX_CYCLIC: NCPU < 1 ' &
            // str(NCPU), __FILE__, __LINE__)
      ENDIF
      ! in case there are more ranks than data available 
      IBLOCK = MOD( M-1, NCPU-1 )+1        ! local index array
      IND = FLOOR( (M-1)/REAL(NCPU-1) )+1  ! process id handling M
      ! if out of bounds
      IF ( IBLOCK /= NODE_ME ) THEN
         IND = -1000
      ENDIF
      LOCAL_INDEX_CYCLIC = IND
#endif
   END FUNCTION LOCAL_INDEX_CYCLIC

!*******************************************************************
!>determines global point N_GLOBAL from the local point N_LOCAL
!!and the group id ID.
!!distribution is assumed block cyclic (round robin)
!*******************************************************************

   FUNCTION LOCAL_INDEX_TO_GLOBAL( NLOCAL, ID, NG, NTOT )
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INTEGER :: LOCAL_INDEX_TO_GLOBAL   !global point
      INTEGER :: NLOCAL            !local point in group
      INTEGER :: ID                !number of group
      INTEGER :: NG                !total number of groups
      INTEGER, OPTIONAL :: NTOT    !size of global array

      !the global index is 
      LOCAL_INDEX_TO_GLOBAL = ID + ( NLOCAL - 1 ) * NG

      IF ( PRESENT(NTOT)) THEN
         IF ( LOCAL_INDEX_TO_GLOBAL > NTOT ) THEN
            CALL vtutor%error("ERROR in LOCAL_INDEX_TO_GLOBAL: global index is out bounds " // &
               str(LOCAL_INDEX_TO_GLOBAL) // " " // str(NTOT) // " " //&
               str( ID ) // " " // str( NLOCAL ) // " " // str( NG ) )  
         ENDIF
      ENDIF
   END FUNCTION LOCAL_INDEX_TO_GLOBAL

!*******************************************************************
! allocates an imaginary grid handle
!*******************************************************************
     SUBROUTINE ALLOCATE_IMAG_GRID(GRIDS, NOMEGA)
        use mpimy, ONLY : communic
        TYPE(imag_grid_handle)  :: GRIDS     !handle 
        INTEGER                 :: NOMEGA    !number of frequency points

        NULLIFY( GRIDS%BOS_RE ) 
        NULLIFY( GRIDS%BOS_RE_WEIGHT ) 
        NULLIFY( GRIDS%TO_BOS_RE ) 
        NULLIFY( GRIDS%TO_BOS_RE_ERROR ) 
        NULLIFY( GRIDS%TO_BOS_RE_CONJG ) 
        NULLIFY( GRIDS%TO_BOS_RE_CONJG_ERROR ) 
        IF ( .NOT.ASSOCIATED( GRIDS%BOS_RE ) ) THEN
           ALLOCATE( GRIDS%BOS_RE( NOMEGA ) ) 
           GRIDS%BOS_RE=0
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%BOS_RE_WEIGHT ) ) THEN
           ALLOCATE( GRIDS%BOS_RE_WEIGHT( NOMEGA ) ) 
           GRIDS%BOS_RE_WEIGHT=0
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_BOS_RE ) ) THEN
           ALLOCATE( GRIDS%TO_BOS_RE( NOMEGA, NOMEGA ) ) 
           CALL SET_IDENTITY( GRIDS%TO_BOS_RE ) 
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_BOS_RE_ERROR) ) THEN
           ALLOCATE( GRIDS%TO_BOS_RE_ERROR( NOMEGA ) ) 
           GRIDS%TO_BOS_RE_ERROR=0
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_BOS_RE_CONJG ) ) THEN
           ALLOCATE( GRIDS%TO_BOS_RE_CONJG( NOMEGA, NOMEGA ) ) 
           CALL SET_IDENTITY( GRIDS%TO_BOS_RE_CONJG ) 
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_BOS_RE_CONJG_ERROR) ) THEN
           ALLOCATE( GRIDS%TO_BOS_RE_CONJG_ERROR( NOMEGA ) ) 
           GRIDS%TO_BOS_RE_CONJG_ERROR=0
        ENDIF

        NULLIFY( GRIDS%BOS_IM ) 
        NULLIFY( GRIDS%BOS_IM_WEIGHT ) 
        NULLIFY( GRIDS%TO_BOS_IM ) 
        NULLIFY( GRIDS%TO_BOS_IM_ERROR ) 
        NULLIFY( GRIDS%TO_BOS_IM_CONJG ) 
        NULLIFY( GRIDS%TO_BOS_IM_CONJG_ERROR ) 
        IF ( .NOT.ASSOCIATED( GRIDS%BOS_IM ) ) THEN
           ALLOCATE( GRIDS%BOS_IM( NOMEGA ) ) 
           GRIDS%BOS_IM=0
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%BOS_IM_WEIGHT ) ) THEN
           ALLOCATE( GRIDS%BOS_IM_WEIGHT( NOMEGA ) ) 
           GRIDS%BOS_IM_WEIGHT=0
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_BOS_IM ) ) THEN
           ALLOCATE( GRIDS%TO_BOS_IM( NOMEGA, NOMEGA ) ) 
           CALL SET_IDENTITY( GRIDS%TO_BOS_IM ) 
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_BOS_IM_ERROR) ) THEN
           ALLOCATE( GRIDS%TO_BOS_IM_ERROR( NOMEGA ) ) 
           GRIDS%TO_BOS_IM_ERROR=0
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_BOS_IM_CONJG ) ) THEN
           ALLOCATE( GRIDS%TO_BOS_IM_CONJG( NOMEGA, NOMEGA ) ) 
           CALL SET_IDENTITY( GRIDS%TO_BOS_IM_CONJG ) 
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_BOS_IM_CONJG_ERROR) ) THEN
           ALLOCATE( GRIDS%TO_BOS_IM_CONJG_ERROR( NOMEGA ) ) 
           GRIDS%TO_BOS_IM_CONJG_ERROR=0
        ENDIF

        NULLIFY( GRIDS%FER_RE ) 
        NULLIFY( GRIDS%FER_RE_WEIGHT ) 
        NULLIFY( GRIDS%TO_FER_RE ) 
        NULLIFY( GRIDS%TO_FER_RE_ERROR ) 
        NULLIFY( GRIDS%TO_FER_RE_CONJG ) 
        NULLIFY( GRIDS%TO_FER_RE_CONJG_ERROR ) 
        IF ( .NOT.ASSOCIATED( GRIDS%FER_RE ) ) THEN
           ALLOCATE( GRIDS%FER_RE( NOMEGA ) ) 
           GRIDS%FER_RE=0
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%FER_RE_WEIGHT ) ) THEN
           ALLOCATE( GRIDS%FER_RE_WEIGHT( NOMEGA ) ) 
           GRIDS%FER_RE_WEIGHT=0
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_FER_RE ) ) THEN
           ALLOCATE( GRIDS%TO_FER_RE( NOMEGA, NOMEGA ) ) 
           CALL SET_IDENTITY( GRIDS%TO_FER_RE ) 
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_FER_RE_ERROR) ) THEN
           ALLOCATE( GRIDS%TO_FER_RE_ERROR( NOMEGA ) ) 
           GRIDS%TO_FER_RE_ERROR=0
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_FER_RE_CONJG ) ) THEN
           ALLOCATE( GRIDS%TO_FER_RE_CONJG( NOMEGA, NOMEGA ) ) 
           CALL SET_IDENTITY( GRIDS%TO_FER_RE_CONJG ) 
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_FER_RE_CONJG_ERROR) ) THEN
           ALLOCATE( GRIDS%TO_FER_RE_CONJG_ERROR( NOMEGA ) ) 
           GRIDS%TO_FER_RE_CONJG_ERROR=0
        ENDIF

        NULLIFY( GRIDS%FER_IM ) 
        NULLIFY( GRIDS%FER_IM_WEIGHT ) 
        NULLIFY( GRIDS%TO_FER_IM ) 
        NULLIFY( GRIDS%TO_FER_IM_ERROR ) 
        NULLIFY( GRIDS%TO_FER_IM_CONJG ) 
        NULLIFY( GRIDS%TO_FER_IM_CONJG_ERROR ) 
        IF ( .NOT.ASSOCIATED( GRIDS%FER_IM ) ) THEN
           ALLOCATE( GRIDS%FER_IM( NOMEGA ) ) 
           GRIDS%FER_IM=0
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%FER_IM_WEIGHT ) ) THEN
           ALLOCATE( GRIDS%FER_IM_WEIGHT( NOMEGA ) ) 
           GRIDS%FER_IM_WEIGHT=0
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_FER_IM ) ) THEN
           ALLOCATE( GRIDS%TO_FER_IM( NOMEGA, NOMEGA ) ) 
           CALL SET_IDENTITY( GRIDS%TO_FER_IM ) 
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_FER_IM_ERROR) ) THEN
           ALLOCATE( GRIDS%TO_FER_IM_ERROR( NOMEGA ) ) 
           GRIDS%TO_FER_IM_ERROR=0
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_FER_IM_CONJG ) ) THEN
           ALLOCATE( GRIDS%TO_FER_IM_CONJG( NOMEGA, NOMEGA ) ) 
           CALL SET_IDENTITY( GRIDS%TO_FER_IM_CONJG ) 
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_FER_IM_CONJG_ERROR) ) THEN
           ALLOCATE( GRIDS%TO_FER_IM_CONJG_ERROR( NOMEGA ) ) 
           GRIDS%TO_FER_IM_CONJG_ERROR=0
        ENDIF

        NULLIFY( GRIDS%TAU ) 
        NULLIFY( GRIDS%TAU_WEIGHT ) 
        NULLIFY( GRIDS%TO_TAU0 )
        NULLIFY( GRIDS%TO_TAU0_ERROR ) 
        NULLIFY( GRIDS%TO_TAU0_CONJG )
        NULLIFY( GRIDS%TO_TAU0_CONJG_ERROR )
        IF ( .NOT.ASSOCIATED( GRIDS%TAU ) ) THEN
           ALLOCATE( GRIDS%TAU( NOMEGA ) ) 
           GRIDS%TAU=0
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TAU_WEIGHT ) ) THEN
           ALLOCATE( GRIDS%TAU_WEIGHT( NOMEGA ) ) 
           GRIDS%TAU_WEIGHT=0
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_TAU0 ) ) THEN
           ALLOCATE( GRIDS%TO_TAU0( NOMEGA, NOMEGA ) ) 
           CALL SET_IDENTITY( GRIDS%TO_TAU0 ) 
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_TAU0_ERROR) ) THEN
           ALLOCATE( GRIDS%TO_TAU0_ERROR( NOMEGA ) ) 
           GRIDS%TO_TAU0_ERROR=0
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_TAU0_CONJG ) ) THEN
           ALLOCATE( GRIDS%TO_TAU0_CONJG( NOMEGA, NOMEGA ) ) 
           CALL SET_IDENTITY( GRIDS%TO_TAU0_CONJG ) 
        ENDIF
        IF ( .NOT.ASSOCIATED( GRIDS%TO_TAU0_CONJG_ERROR) ) THEN
           ALLOCATE( GRIDS%TO_TAU0_CONJG_ERROR( NOMEGA ) ) 
           GRIDS%TO_TAU0_CONJG_ERROR=0
        ENDIF

        CONTAINS 
        SUBROUTINE SET_IDENTITY( MAT ) 
           USE string, ONLY: STR
           USE tutor, ONLY: VTUTOR
           REAL(q) :: MAT( : , : )  
           !local 
           INTEGER :: I 

           MAT = 0 
           IF ( SIZE( MAT, 1 ) /= SIZE( MAT, 2 ) ) THEN
              CALL VTUTOR%ERROR( " SET_IDENTITY reports inconsistent array"//&
               " sizes: "//STR( SIZE( MAT, 1 ) )//" vs "//STR( SIZE(MAT, 2))) 
           ENDIF 

           DO I = 1, SIZE( MAT , 1 ) 
              MAT( I, I ) = 1._q            
           ENDDO

        END SUBROUTINE SET_IDENTITY

     END SUBROUTINE ALLOCATE_IMAG_GRID

!*******************************************************************
!> helper: synchronizes quaratures 
!*******************************************************************
   SUBROUTINE SYNC_QUADRATURES( GRIDS ) 
      TYPE( imag_grid_handle ) :: GRIDS 
#ifdef MPI
      PROFILING_START( 'SYNC_QUADRATURES' )
      ! sync quadratures, if they have been calculated and set 
      IF( GRIDS%FREQ_BOS_RE%N > 0 ) THEN
         CALL SYNC_QUADRATURE( GRIDS%COMM_INTER, GRIDS%COMM, GRIDS%FREQ_BOS_RE ) 
      ENDIF
      IF( GRIDS%FREQ_BOS_IM%N > 0 ) THEN
         CALL SYNC_QUADRATURE( GRIDS%COMM_INTER, GRIDS%COMM, GRIDS%FREQ_BOS_IM ) 
      ENDIF
      IF( GRIDS%FREQ_FER_RE%N > 0 ) THEN
         CALL SYNC_QUADRATURE( GRIDS%COMM_INTER, GRIDS%COMM, GRIDS%FREQ_FER_RE ) 
      ENDIF
      IF( GRIDS%FREQ_FER_IM%N > 0 ) THEN
         CALL SYNC_QUADRATURE( GRIDS%COMM_INTER, GRIDS%COMM, GRIDS%FREQ_FER_IM ) 
      ENDIF
      
      PROFILING_STOP( 'SYNC_QUADRATURES' )

      CONTAINS
      SUBROUTINE SYNC_QUADRATURE( COMM_INTER, COMM, QUAD ) 
         USE string, ONLY: str
         USE tutor, ONLY: vtutor
         TYPE( communic )          :: COMM
         TYPE( communic )          :: COMM_INTER
         TYPE( quadrature_handle ) :: QUAD
         INTEGER                   :: ID
         PROFILING_START( 'sync_quadrature' )

         ! quick return if no distribution has been performed
         IF ( COMM_INTER%NCPU == 1) THEN
            PROFILING_STOP( 'sync_quadrature' )
            RETURN
         ENDIF

         ! otherwise broadcast results between groups
         ! the group master rank sends the result to the others
         ID = 0 
         IF ( COMM_INTER%NODE_ME == QUAD%GRID_ID ) THEN
            IF ( QUAD%COMM%NODE_ME == 1 ) ID = COMM%NODE_ME
         ENDIF    
         CALLMPI( M_sum_i( COMM, ID, 1 ) )
         IF ( ID < 1 .OR. ID>COMM%NCPU ) THEN
            CALL vtutor%bug("internal ERROR, SYNC_QUADRATURE reports ID out of scope " // &
               str(QUAD%ITYPE) // " " // str(ID) // " " // str(COMM%NCPU), __FILE__, __LINE__)
         ENDIF
#ifdef debug
         CALL DUMP_QUADRATURE( QUAD, ABS(QUAD%ITYPE)+100*COMM%NODE_ME, &
         FNAME="before" )
#endif
         ! broadcast coefficients
         CALLMPI( M_bcast_i_from( COMM, QUAD%ITYPE, 1, ID ))
         CALLMPI( M_bcast_i_from( COMM, QUAD%N, 1, ID ))
         CALLMPI( M_bcast_qd_from( COMM, QUAD%A, 1, ID ))
         CALLMPI( M_bcast_qd_from( COMM, QUAD%B, 1, ID ))
         CALLMPI( M_bcast_qd_from( COMM, QUAD%SCALING, 1, ID ))
         CALLMPI( M_bcast_qd_from( COMM, QUAD%INFINITYNORM, 1, ID )) 
         CALLMPI( M_bcast_qd_from( COMM, QUAD%C(1), 2*QUAD%N+1, ID ))
         CALLMPI( M_bcast_qd_from( COMM, QUAD%X0(1), 2*QUAD%N+1, ID )) 

#ifdef debug
         CALL DUMP_QUADRATURE( QUAD, ABS(QUAD%ITYPE)+100*COMM%NODE_ME, &
         FNAME="after" )
#endif
      PROFILING_STOP( 'sync_quadrature' )
   END SUBROUTINE SYNC_QUADRATURE
#endif
      END SUBROUTINE SYNC_QUADRATURES

!*******************************************************************
! copies one quadrature to another quadrature b
! if b is not associated it allocates b
!*******************************************************************
   SUBROUTINE COPY_QUADRATURE( QUAD_A, QUAD_B, LSINGLE_POINT )
      TYPE( quadrature_handle )  :: QUAD_A
      TYPE( quadrature_handle )  :: QUAD_B
      LOGICAL                    :: LSINGLE_POINT
      !local
      INTEGER                    :: I 

      ! decide if single point should be allocated 
      ! used for creating a dummy quadrature for the calculation 
      ! of the tau->0 limit weights 
      QUAD_B%N = QUAD_A%N
      IF ( LSINGLE_POINT )  QUAD_B%N = 1 

      NULLIFY( QUAD_B%C ) 
      IF ( ASSOCIATED( QUAD_B%C ) ) DEALLOCATE( QUAD_B%C )
      IF(.NOT.ASSOCIATED(QUAD_B%C )) ALLOCATE(QUAD_B%C(2*QUAD_B%N+1)) 
      NULLIFY( QUAD_B%X0 ) 
      IF ( ASSOCIATED( QUAD_B%X0 ) ) DEALLOCATE( QUAD_B%X0 )
      IF(.NOT.ASSOCIATED(QUAD_B%X0 )) ALLOCATE(QUAD_B%X0(2*QUAD_B%N+1)) 
      DO I = 1, 2*QUAD_B%N+1
         QUAD_B%C(I) = QDPP(0)
         QUAD_B%X0(I) = QDPP(0)
      ENDDO
      QUAD_B%INFINITYNORM = QDPP(0)
      QUAD_B%A = QUAD_A%A
      QUAD_B%B = QUAD_A%B
      QUAD_B%ITYPE = QUAD_A%ITYPE

      ! in case full quadrature is copied
      IF ( QUAD_B%N == QUAD_A%N ) THEN
         DO I = 1, 2*QUAD_B%N+1
            QUAD_B%C(I) = QUAD_A%C(I)
            QUAD_B%X0(I) = QUAD_A%X0(I)
         ENDDO
         QUAD_B%INFINITYNORM = QUAD_A%INFINITYNORM
      ENDIF

      ! point to the same functions 
      IF( ASSOCIATED(QUAD_A%PHI) )               QUAD_B%PHI                => QUAD_A%PHI               
      IF( ASSOCIATED(QUAD_A%PSI) )               QUAD_B%PSI                => QUAD_A%PSI
      IF( ASSOCIATED(QUAD_A%D_PHI2_DZ) )         QUAD_B%D_PHI2_DZ          => QUAD_A%D_PHI2_DZ
      IF( ASSOCIATED(QUAD_A%D_PHI2_DL) )         QUAD_B%D_PHI2_DL          => QUAD_A%D_PHI2_DL
      IF( ASSOCIATED(QUAD_A%D2_PHI2_DZ2) )       QUAD_B%D2_PHI2_DZ2        => QUAD_A%D2_PHI2_DZ2 
      IF( ASSOCIATED(QUAD_A%F) )                 QUAD_B%F                  => QUAD_A%F
      IF( ASSOCIATED(QUAD_A%D_F_DZ) )            QUAD_B%D_F_DZ             => QUAD_A%D_F_DZ
      IF( ASSOCIATED(QUAD_A%D2_F_DZ2) )          QUAD_B%D2_F_DZ2           => QUAD_A%D2_F_DZ2
      IF( ASSOCIATED(QUAD_A%TRANS_TAYLOR) )      QUAD_B%TRANS_TAYLOR       => QUAD_A%TRANS_TAYLOR 
      IF( ASSOCIATED(QUAD_A%R1E3) )              QUAD_B%R1E3               => QUAD_A%R1E3  
      IF( ASSOCIATED(QUAD_A%R1E4) )              QUAD_B%R1E4               => QUAD_A%R1E4  
      IF( ASSOCIATED(QUAD_A%R1E5) )              QUAD_B%R1E5               => QUAD_A%R1E5  
      IF( ASSOCIATED(QUAD_A%R1E6) )              QUAD_B%R1E6               => QUAD_A%R1E6  
      IF( ASSOCIATED(QUAD_A%CTAB) )              QUAD_B%CTAB               => QUAD_A%CTAB  
      IF( ASSOCIATED(QUAD_A%PHI_CONJG) )         QUAD_B%PHI_CONJG          => QUAD_A%PHI_CONJG  
      IF( ASSOCIATED(QUAD_A%PSI_CONJG) )         QUAD_B%PSI_CONJG          => QUAD_A%PSI_CONJG  
      IF( ASSOCIATED(QUAD_A%TRANS_TAYLOR_CONJG) )QUAD_B%TRANS_TAYLOR_CONJG => QUAD_A%TRANS_TAYLOR_CONJG
      IF( ASSOCIATED(QUAD_A%FITCOF) )            QUAD_B%FITCOF             => QUAD_A%FITCOF
   END SUBROUTINE COPY_QUADRATURE 

END MODULE minimax
