#include "symbol.inc"

!*****************************************************************************
!
!  Unitary Matrix Constrained Optimization (UMCO)
!
!  The UMCO module find a unitary transformation that maximizes/minimimzes
!  some cost function (e.g. orbital localization)
!  The conjugate gradient algorithm and notation is used from
!  Lehtola et al., JCTC 2013, 9, 5365 (dx.doi.org/10.1021/ct400793q)
! 
!  01/2021: birth of this module (Tobias Schäfer)
!
!*****************************************************************************


!> @brief Module to provide functional derivatives of Pipek-Mezey functional
MODULE umco_pipek_mezey
  USE base
  USE locproj, ONLY: LPRJ_NUM_ORBITALS_ON_ATOM, LPRJ_GET_NUM_WAN, LPRJ_GET_SITE
#ifdef PROFILING
    USE profiling
#endif
  IMPLICIT NONE

  PUBLIC :: UMCO_PM_INIT, &
            UMCO_PM_DEALLOC, &
            UMCO_PM_EUCLIDEAN_DERIVATIVE, &
            UMCO_PM_COSTFUNCTION, &
            PM_PROJTYPE, PM_PROJTYPE_NONE, PM_PROJTYPE_IAO, PM_PROJTYPE_LOCPROJ

  ENUM, BIND(C)
    ENUMERATOR :: PM_PROJTYPE_NONE, PM_PROJTYPE_IAO, PM_PROJTYPE_LOCPROJ
  END ENUM

  PRIVATE

  INTEGER :: NATOMS  ! number of atoms for PM functional 
  INTEGER :: NWF     ! number of wannier functions to construct 
  INTEGER :: POWER   ! the exponent in the cost PM function

  ! how to build the projector on atomic charges?
  ! PM_PROJTYPE_NONE = don't use
  ! PM_PROJTYPE_IAO = use intrinsic atomic orbitals => produces intrinsic bonding oritbals
  ! PM_PROJTYPE_LOCPROJ = use trials functions as defined by LOCPROJ
  INTEGER :: PM_PROJTYPE

  ! overlap with atmoic functions: GOVL(i:mu:A) = <i|mu_A>
  GDEF, ALLOCATABLE :: GOVL(:,:,:) ! overlap matrix per atom

  CONTAINS

  !> Initialize Pipek-Mezey functional
  SUBROUTINE UMCO_PM_INIT(LPRJ_COVL, ISP, IWF_LOW, IWF_HIGH, T_INFO, IO)
    USE wave
    USE lattice
    USE pseudo
    USE poscar_struct_def
    USE reader_tags

    GDEF,INTENT(IN)     :: LPRJ_COVL(:,:,:,:)
    INTEGER             :: ISP
    INTEGER             :: IWF_LOW
    INTEGER             :: IWF_HIGH
    TYPE(type_info)     :: T_INFO
    TYPE(in_struct)     :: IO
    ! local
    INTEGER           :: NAO_TOT ! number of atomic orbitals (total) 
    INTEGER           :: IK, I, J
    INTEGER           :: MAX_NAO1
    INTEGER           :: IATOM, IAO, IAO1
    REAL(q)           :: SCALEF
    ! test
    GDEF, ALLOCATABLE    :: AD(:,:,:)
    INTEGER              :: MU, NU
    REAL(q), ALLOCATABLE :: REVL(:)  ! real eigenvalues 
    GDEF, ALLOCATABLE    :: WORK(:)  ! work array
    REAL(q), ALLOCATABLE :: RWORK(:) ! work array
    INTEGER              :: IERR
    INTEGER              :: LWORK    ! size of work array
    GDEF                 :: WORK_DIM
    ! debug
    REAL(q)              :: SUMME
    GDEF                 :: MAN
    INTEGER              :: FUNIT

    ! gamma only  at the moment
    IK = 1
    
    ! default value for exponent
    POWER = 2

    ! get number of orbitals to transform
    NWF = IWF_HIGH-IWF_LOW+1 

    ! get number of atomic functions from LOCPROJ
    NAO_TOT = LPRJ_GET_NUM_WAN() 

    ! store total number of atoms (ions)
    NATOMS = T_INFO%NIONP

    ! find max number of local functions per atom
    MAX_NAO1=0
    DO I = 1, NATOMS
      MAX_NAO1 = MAX(MAX_NAO1, LPRJ_NUM_ORBITALS_ON_ATOM(I))
    ENDDO

    ! allocate and fill the projector
    ALLOCATE(GOVL(NWF,MAX_NAO1,NATOMS))
    GOVL = zero 

    ! store overlaps per atom
    DO IATOM = 1, NATOMS
      IAO1 = 1 ! reset local atomic function index
      DO IAO = 1, NAO_TOT
        IF(LPRJ_GET_SITE(IAO) == IATOM) THEN
          GOVL(:,IAO1,IATOM) = LPRJ_COVL(IWF_LOW:IWF_HIGH, IK, ISP, IAO)
          IAO1 = IAO1 + 1
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE UMCO_PM_INIT



  !*****************************************************************************
  !
  ! Provides the functional derivative of the Pipek-Mezey (PM) cost function L
  ! with respect to the complex conjugate unitary matrix:
  !
  ! dL/du*(s,t) = p * sum_A [gX](s,t) * ( [X^H X](t,t) )**(p-1)
  !
  ! where p is the exponentiation in the PM functional, and
  ! X = g^H U
  ! g = overlap matrix (on the atom A)
  ! U = the unitary matrix
  !
  !*****************************************************************************
  SUBROUTINE UMCO_PM_EUCLIDEAN_DERIVATIVE(U, GAM)
    GDEF, intent(in)  :: U(NWF,NWF)   ! the unitary trafo matrix
    GDEF, intent(out) :: GAM(NWF,NWF) ! the euclidean derivative
    ! local
    INTEGER           :: IATOM, NAO1, NMU, MU, I, J
    GDEF, ALLOCATABLE :: GWORK(:,:)
    GDEF, ALLOCATABLE :: GWORK2(:,:)
    GDEF, ALLOCATABLE :: GWORKD(:)

    PROFILING_START("UMCO_PM_EUCLIDEAN_DERIVATIVE")

    ALLOCATE(GWORKD(NWF))
    ALLOCATE(GWORK2(NWF,NWF))

    GAM = zero

    DO IATOM = 1, NATOMS

      NMU = LPRJ_NUM_ORBITALS_ON_ATOM(IATOM)
      ALLOCATE(GWORK(NMU,NWF))
    
      !GWORK = GCONJG(TRANSPOSE(GOVL(:,1:NMU,IATOM)))
      !GWORK = MATMUL(GWORK,U)
      ! GWORK = GOVL^H*U
      CALL GGEMM('C', 'N', NMU, NWF, NWF, one, GOVL(:,1:NMU,IATOM), NWF, U, NWF, zero, GWORK, NMU)

      GWORKD = zero
      DO J = 1, NWF
        DO MU = 1, NMU
          GWORKD(J) = GWORKD(J) + ABS(GWORK(MU,J))**2.0_q
        ENDDO
        GWORKD(J) = GWORKD(J)**(POWER-1)
      ENDDO

      !GWORK2 = MATMUL(GOVL(:,1:NMU,IATOM),GWORK)
      ! GWORK2 = GOVL*GWORK
      CALL GGEMM('N', 'N', NWF, NWF, NMU, one, GOVL(:,1:NMU,IATOM), NWF, GWORK, NMU, zero, GWORK2, NWF)
      DO J = 1, NWF
        GWORK2(:,J) = GWORK2(:,J) * GWORKD(J)
      ENDDO

      GAM = GAM + GWORK2

      DEALLOCATE(GWORK)
    ENDDO

    GAM = GAM * POWER

    !! for ISPIN=1 
    !GAM = GAM * 2

    DEALLOCATE(GWORK2)
    DEALLOCATE(GWORKD)

    PROFILING_STOP("UMCO_PM_EUCLIDEAN_DERIVATIVE")

  END SUBROUTINE UMCO_PM_EUCLIDEAN_DERIVATIVE


  FUNCTION UMCO_PM_COSTFUNCTION(U) RESULT(L)
    GDEF, intent(in)  :: U(NWF,NWF)   ! the unitary trafo matrix
    REAL(q) :: L ! the cost function evaluated at U
    ! local
    INTEGER           :: IATOM, NAO1, NMU, MU, I, J
    REAL(q)           :: DIAG
    GDEF, ALLOCATABLE :: GWORK(:,:)

    PROFILING_START("UMCO_PM_COSTFUNCTION")

    L = 0.0_q

    DO IATOM = 1, NATOMS

      NMU = LPRJ_NUM_ORBITALS_ON_ATOM(IATOM)
      ALLOCATE(GWORK(NMU,NWF))
    
      GWORK = GCONJG(TRANSPOSE(GOVL(:,1:NMU,IATOM)))
      GWORK = MATMUL(GWORK,U)

      DO J = 1, NWF
        DIAG = 0.0_q
        DO MU = 1, NMU
          DIAG = DIAG + ABS(GWORK(MU,J))**2.0_q
        ENDDO
        L = L + DIAG**POWER
      ENDDO

      DEALLOCATE(GWORK)
    ENDDO

    PROFILING_STOP("UMCO_PM_COSTFUNCTION")

  END FUNCTION UMCO_PM_COSTFUNCTION


  SUBROUTINE UMCO_PM_DEALLOC
    IF (ALLOCATED(GOVL)) DEALLOCATE(GOVL)
  END SUBROUTINE UMCO_PM_DEALLOC

END MODULE umco_pipek_mezey




!> @brief Module to find unitary trafo of orbitals to max/minimize a functional
MODULE umco
  USE base
  USE locproj
  USE umco_pipek_mezey
  USE string, ONLY: str
  USE tutor, ONLY: vtutor
#ifdef PROFILING
    USE profiling
#endif
  IMPLICIT NONE

  PUBLIC :: IS_UMCO ! check if mocule is active
  PUBLIC :: UMCO_READER
  PUBLIC :: UMCO_CALC_TRAFO
  PUBLIC :: UMCO_UTRAFO

  GDEF, ALLOCATABLE :: UMCO_UTRAFO(:,:,:,:) ! dimensions are (band,band,kpoint,spin)



  PRIVATE

  ! the unitary transformation matrix
  REAL(q)           :: SGN    ! exp(SGN*mu*H) for update of unitary
  INTEGER           :: UMCO_NBLOW(2)  ! lower band for which to perform UMCO (for spin up and down)
  INTEGER           :: UMCO_NBHIGH(2) ! upper band for which to perform UMCO (for spin up and down)
  INTEGER           :: NWF    ! number of wannier functions to construct 
  INTEGER           :: UPOWER ! order of U in cost function
  INTEGER           :: POL    ! order of polynomial for line search algo
  LOGICAL           :: LUMCO
  REAL(q)           :: GRADIENT_NORM_STOP ! break condition for gradient norm 

  CONTAINS


  SUBROUTINE UMCO_READER(INFO, IO)
    USE reader_tags
    TYPE(info_struct) :: INFO
    Type(in_struct)   :: IO
    ! local
    INTEGER :: N, NUMB
    INTEGER :: IERR

    LUMCO = .FALSE.
    CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'LUMCO', LUMCO, IERR, WRITEXMLINCAR)

    ! read nblow
    CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'UMCO_NBLOW', N, IERR, FOUNDNUMBER=NUMB)
    IF (NUMB==0) THEN
      UMCO_NBLOW = -1
    ELSE
      CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'UMCO_NBLOW', UMCO_NBLOW, MIN(NUMB,2), IERR)
      IF (NUMB==1) UMCO_NBLOW(2)=UMCO_NBLOW(1)
    ENDIF

    ! read nbhigh
    CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'UMCO_NBHIGH', N, IERR, FOUNDNUMBER=NUMB)
    IF (NUMB==0) THEN
      UMCO_NBHIGH = -1
    ELSE
      CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'UMCO_NBHIGH', UMCO_NBHIGH, MIN(NUMB,2), IERR)
      IF (NUMB==1) UMCO_NBHIGH(2)=UMCO_NBHIGH(1)
    ENDIF

    ! default value for projector on atomic charges
    PM_PROJTYPE = PM_PROJTYPE_NONE
    CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'UMCO_PM_PROJTYPE', PM_PROJTYPE, IERR, WRITEXMLINCAR)
    IF (N==1) PM_PROJTYPE = PM_PROJTYPE_IAO
    IF (N==2) PM_PROJTYPE = PM_PROJTYPE_LOCPROJ

    GRADIENT_NORM_STOP = 1E-10_q
    CALL PROCESS_INCAR(IO%LOPEN, IO%IU0, IO%IU5, 'UMCO_GRAD_BREAK', GRADIENT_NORM_STOP, IERR, WRITEXMLINCAR)


  END SUBROUTINE UMCO_READER



  FUNCTION IS_UMCO() RESULT(IS_LUMCO)
    LOGICAL :: IS_LUMCO
    IS_LUMCO = LUMCO
  END FUNCTION IS_UMCO

  !> Get the maximum number of wannier orbitals produced by UMCO module
  FUNCTION UMCO_GET_NUM_WANN(W) RESULT(NW_MAX)
    use wave_struct_def, ONLY : wavespin
    TYPE(wavespin)      :: W
    ! local
    INTEGER :: NW_MAX
    INTEGER :: IK
    INTEGER :: ISP, IWF_LOW, IWF_HIGH

    IK=1

    ! determine dimensions for U matrix
    NW_MAX=0
    DO ISP=1,W%WDES%ISPIN
      IWF_LOW  = 1
      IWF_HIGH = NINT(SUM(W%FERTOT(:, IK, ISP)))
      ! overwrite with user setting
      IF (UMCO_NBLOW(ISP)  .GT. 0) IWF_LOW  = UMCO_NBLOW(ISP)
      IF (UMCO_NBHIGH(ISP) .GT. 0) IWF_HIGH = UMCO_NBHIGH(ISP)
      NW_MAX = MAX(IWF_HIGH,NW_MAX)
    ENDDO

  END FUNCTION UMCO_GET_NUM_WANN

  !> Use unitary matrix constrained optimization to minimize the Pipek-Mezey functional
  !> and obtain localized orbitals
  SUBROUTINE UMCO_CALC_TRAFO(LPRJ_COVL, UMCO_UTRAFO, W, GRID, P, CQIJ, LATT_CUR, T_INFO, INFO, IO)
    USE wave
    USE wave_high
    USE dfast
    USE fileio
    USE lattice
    USE pseudo
    USE poscar_struct_def
    USE constant

    GDEF,INTENT(IN) :: LPRJ_COVL(:,:,:,:)
    !> final transformation matrix between the bloch states and localized orbitals
    GDEF, ALLOCATABLE, INTENT(OUT):: UMCO_UTRAFO(:,:,:,:)
    TYPE(wavespin)      :: W
    TYPE(grid_3d)       :: GRID
    TYPE(potcar)        :: P(:)
    OVERLAP             :: CQIJ(:,:,:,:)
    TYPE(latt)          :: LATT_CUR
    TYPE(type_info)     :: T_INFO
    TYPE(info_struct)   :: INFO
    TYPE(in_struct)     :: IO
    ! local
    TYPE(wavefuna)          :: WA
    TYPE(wavedes1)          :: WDESK
    INTEGER                 :: ISP, IK, ITER, MAXITER, MINITER, I, J, N
    INTEGER                 :: IWF_LOW, IWF_HIGH
    REAL(q)                 :: CGPR_FACTOR, GRADIENT_NORM
    GDEF, ALLOCATABLE       :: U(:,:)
    GDEF, ALLOCATABLE       :: GWORK(:,:)
    COMPLEX(q), ALLOCATABLE :: CWORK(:,:)
    REAL(q), ALLOCATABLE    :: RWORK(:,:), RWORK2(:)
    ! to store Riemannian derivative of cost function of previous iteration step
    GDEF, ALLOCATABLE       :: PRIOR_RIEM_DERV(:,:) 
    GDEF, ALLOCATABLE       :: DJDUH(:,:)  ! eclidean derivative of cost function
    GDEF, ALLOCATABLE       :: H(:,:)      ! ascent direction
    COMPLEX(q), ALLOCATABLE :: V(:,:)      ! unitary to diagonalize H
    COMPLEX(q), ALLOCATABLE :: EVL(:)      ! eigenvalues of H
    REAL(q)                 :: MAX_ABS_EVL ! largest absolute value in EVL
    REAL(q), ALLOCATABLE    :: MU(:)       ! MU array for polynom. line search
    REAL(q)                 :: MUROOT
    REAL(q), ALLOCATABLE    :: DJDMU(:)    ! deriv. of cost function wrt mu
    COMPLEX(q), ALLOCATABLE :: EXPDIAG(:)  ! exp(mu EVL)
    GDEF, ALLOCATABLE       :: U2(:,:)     ! updated unitary U2=exp(SGN*mu*H)U
    INTEGER, ALLOCATABLE    :: IPIV(:)     ! pivot indices
    INTEGER                 :: IERR
    INTEGER                 :: WINNER_NODE
    INTEGER                 :: NODE_ME, NCPU
    INTEGER                 :: NW_MAX
    INTEGER                 :: IB1,IB2
    REAL(q)                 :: WORK_DIM, DUMMY
    REAL(q), ALLOCATABLE    :: WI(:), WR(:) ! real and imag part of eigenvals
    REAL(q), ALLOCATABLE    :: FERTOT(:,:,:)
    REAL(q)                 :: CURR_COSTFUNC, PREV_COSTFUNC, DIFF_COSTFUNC
    REAL(q), ALLOcATABLE    :: COSTFUNC_ON_NODE(:)
    ! debug
    GDEF :: t1,t2
    GDEF, ALLOCATABLE       :: HAMIL(:,:)
    INTEGER                 :: FUNIT

#ifndef gammareal
    CALL vtutor%error('UMCO_CALC_TRAFO only works in vasp_gam')
#endif

    PROFILING_START('UMCO_CALC_TRAFO')

    ! status output
    IF(IO%IU0 .GT. 0) THEN
      WRITE(IO%IU0,*) "UMCO mode" 
    ENDIF

    NODE_ME = 1
    NCPU = 1
#ifdef MPI
    NODE_ME = W%WDES%COMM%NODE_ME
    NCPU = W%WDES%COMM%NCPU
#endif

    ! backup FERTOT 
    ALLOCATE(FERTOT(W%WDES%NB_TOT, W%WDES%NKPTS, W%WDES%ISPIN))
    FERTOT = W%FERTOT

    ! gamma only 
    IK = 1

    NW_MAX = UMCO_GET_NUM_WANN(W)

    ! allocate storage for transformation matrix
    IF(ALLOCATED(UMCO_UTRAFO)) DEALLOCATE(UMCO_UTRAFO)
    ALLOCATE(UMCO_UTRAFO(NW_MAX,NW_MAX,1,W%WDES%ISPIN))
    UMCO_UTRAFO = zero

    ! set reasonable value for max iterations
    MINITER = 8
    MAXITER = 100000

    ! FUNCTIONAL SPECIFIC SETTINGS
    UPOWER = 4   ! PM at the moment
    POL = 5      ! order of polynomial (for line search)
    SGN = +1.0_q ! PM -> maximize cost function (+1 = maximize, -1=mimize)


    ! init the random number generator for random initial trafo matrix
    CALL INIT_RND_GENERATOR(NODE_ME)

    DO ISP = 1, W%WDES%ISPIN

      ! we allow to specify a band window  (default = occupied)
      IWF_LOW = 1
      IWF_HIGH = NINT(SUM(W%FERTOT(:, IK, ISP)))
      ! overwrite with user setting
      IF (UMCO_NBLOW(ISP)  .GT. 0) IWF_LOW  = UMCO_NBLOW(ISP)
      IF (UMCO_NBHIGH(ISP) .GT. 0) IWF_HIGH = UMCO_NBHIGH(ISP)
      NWF = IWF_HIGH-IWF_LOW+1
  
      IF(NWF .LT. 2) THEN
        CALL vtutor%alert("module UMCO: just one band given, nothing to &
                           optimize... so I do nothing and return.")
      ENDIF

      ! init Pipek-Mezey
      CALL UMCO_PM_INIT(LPRJ_COVL, ISP, IWF_LOW, IWF_HIGH, T_INFO, IO)

      ! output file 
      FUNIT=700+NODE_ME
      OPEN(unit = FUNIT,file = "UMCOiter_ISP" // str(ISP) // "_NODE" // str(NODE_ME), FORM='FORMATTED', access='stream', STATUS='REPLACE')

      ! output
      IF(IO%IU0 .GT. 0) THEN
        WRITE(IO%IU0,'(A,I1,A,I1)') " working on spin ", ISP, " of ", W%WDES%ISPIN
        WRITE(IO%IU0,*) "  ITER   GRAD_NORM     DIFF    CGPR_FACTOR      MUROOT         COST_FUNCTION   NODE"
      ENDIF
      WRITE(FUNIT,'(A,I1,A,I1)') " working on spin ", ISP, " of ", W%WDES%ISPIN
      WRITE(FUNIT,*) "  ITER   GRAD_NORM     DIFF    CGPR_FACTOR      MUROOT         COST_FUNCTION"
  
      ALLOCATE(GWORK(NWF,NWF))
      ALLOCATE(PRIOR_RIEM_DERV(NWF,NWF))
      ALLOCATE(DJDUH(NWF,NWF))
      ALLOCATE(H(NWF,NWF))
      GWORK = zero
      PRIOR_RIEM_DERV = zero
      H = zero
  
      ! allocations for iteration loop
      ALLOCATE(DJDMU(POL+1))
      ALLOCATE(V(NWF,NWF), EVL(NWF))
      ALLOCATE(CWORK(NWF,NWF))
      ALLOCATE(EXPDIAG(NWF))
      ALLOCATE(U2(NWF,NWF))
      ALLOCATE(MU(POL+1))
      ALLOCATE(RWORK(POL,POL))
      ALLOCATE(IPIV(POL))
      ALLOCATE(WR(POL), WI(POL))

      ALLOCATE(U(NWF,NWF))

      CALL UMCO_RANDOM_UNITARY(U)
      !U=0
      !DO I = 1, NWF
      !  U(I,I) = 1.0_q
      !ENDDO

      CURR_COSTFUNC = 0.0_q
      PREV_COSTFUNC = 0.0_q
      DIFF_COSTFUNC = 0.0_q

      PROFILING_START('UMCO_ITER')

      ! the main iteration loop
      iterations: &
      DO ITER = 1, MAXITER
  
        ! Calc Euclidean derivative
        CALL UMCO_PM_EUCLIDEAN_DERIVATIVE(U, DJDUH)
  
        ! Calc Riemannian derivative
        GWORK = MATMUL(DJDUH,GCONJG(TRANSPOSE(U))) - MATMUL(U,GCONJG(TRANSPOSE(DJDUH)))
  
        ! calc the gradient norm as a stop criterion
        GRADIENT_NORM = 0.5_q * REAL(GFROBENIUS_PRODUCT(GWORK,GWORK), q) 
  
        ! Calc Polak-Ribière-Polyak update factor for conjugate gradient algorithm
        IF(ITER .GT. 1) THEN
          CGPR_FACTOR = UMCO_CGPR_FACTOR(GWORK, PRIOR_RIEM_DERV)
        ELSE 
          CGPR_FACTOR = 0.0_q
        ENDIF
        PRIOR_RIEM_DERV = GWORK
  
        ! update ascent direction
        H = GWORK + CGPR_FACTOR * H
  
        ! check if set-back of H is necessary
        IF( (REAL(GFROBENIUS_PRODUCT(H, GWORK), q) .LT. 0.0_q) & 
             .OR. ( (MODULO(ITER-1,NWF) .EQ. 0 ) .AND. (ITER .GT. 2) ) ) THEN
          H = GWORK
          CGPR_FACTOR = 0.0_q
        ENDIF
  
        ! 
        ! begin: calc optimal step size (polynomial line search approximation)
        !
        
        ! calc DJDMU for mu=0
        CALL UMCO_CALC_DJDMU(DJDUH, U, H, DJDMU(POL+1))
        ! diagonalize H and store unitary in V and imaginary eigenvalues in EVL 
        CALL UMCO_DIAG_H(H, V, EVL)
        MAX_ABS_EVL = MAXVAL(ABS(EVL))
  
        ! define mu sampling points
        MU(POL+1)=0.0_q
        MU(1)=TPI/(POL*UPOWER*MAX_ABS_EVL)
        DO N=2,POL
          MU(N) = N*MU(1)
        ENDDO
      
        ! construct exp(SGN*mu*H) for mu(1)
        EXPDIAG(:) = EXP( SGN * EVL(:) * MU(1) )
        DO J = 1, NWF
          CWORK(:,J) = V(:,J) * EXPDIAG(J)
        ENDDO
        CWORK = MATMUL(CWORK,CONJG(TRANSPOSE(V)))
  
        ! construct exp(SING*mu*H)*U for mu(1)
#ifdef gammareal
        U2 = MATMUL(REAL(CWORK,KIND=q),U)
#else
        U2 = MATMUL(CWORK,U)
#endif
  
        DO N=1,POL
          ! Calc Euclidean derivative for current mu
          CALL UMCO_PM_EUCLIDEAN_DERIVATIVE(U2, DJDUH)
  
          ! calc DJDMU for current mu
          CALL UMCO_CALC_DJDMU(DJDUH, U2, H, DJDMU(N))
  
          ! construct exp(SING*mu*H)*U  for MU(N)
#ifdef gammareal
          U2 = MATMUL(REAL(CWORK,KIND=q),U2) 
#else
          U2 = MATMUL(CWORK,U2)
#endif
        ENDDO
  
        ! construct polynomial coefficients (polynomial fit)
        ! build matrix to invert
        ! ToDo better do this with DSYSV instead of matrix inversion
        DO J=1,POL
          DO I=1,POL
            RWORK(I,J) = MU(I)**J
          ENDDO
        ENDDO
        CALL DGETRF(POL, POL, RWORK, POL, IPIV, IERR) ! LU decomposition
        CALL DGETRI(POL, RWORK, POL, IPIV, WORK_DIM, -1, IERR) ! query inversion
        I = INT(WORK_DIM)
        ALLOCATE(RWORK2(I))
        CALL DGETRI(POL, RWORK, POL, IPIV, RWORK2, I, IERR) ! perform inversion
        DEALLOCATE(RWORK2)
        ALLOCATE(RWORK2(POL))
        DO N=1,POL
          RWORK2(N) = DJDMU(N) - DJDMU(POL+1)
        ENDDO
        RWORK2 = MATMUL(RWORK,RWORK2) ! final coefficients
  
        ! calculate roots of polynomial via companion matrix (use RWORK)
        RWORK = 0.0_q
        DO I=1,POL-1
          RWORK(I,I+1) = 1.0_q
        ENDDO
        RWORK(POL,1) = -1.0_q * DJDMU(POL+1) / RWORK2(POL)
        DO J=2,POL
          RWORK(POL,J) = -1.0_q * RWORK2(J-1) / RWORK2(POL)
        ENDDO
        DEALLOCATE(RWORK2)
  
        CALL DGEEV('N', 'N', POL, RWORK, POL, WR, WI, DUMMY, POL, & ! query diag.
                   DUMMY, POL, WORK_DIM, -1, IERR)
        I = INT(WORK_DIM)
        ALLOCATE(RWORK2(I))
        CALL DGEEV('N', 'N', POL, RWORK, POL, WR, WI, DUMMY, POL, & ! perf. diag.
                   DUMMY, POL, RWORK2, I, IERR)
        DEALLOCATE(RWORK2)
  
        ! find smalles positive real eigenvalue 
        ! (equivalent to smalles positive real root of polynomial)
        MUROOT=1D+10 ! large number
        DO I=1,POL
          IF( (WI(I) .GT. +1E-8_q) .OR. (WR(I) .LT. 0.0_q) ) CYCLE
          MUROOT=MIN(MUROOT,WR(I))
        ENDDO
  
        ! 
        ! end: calc optimal step size (polynomial line search approximation)
        !
  
        ! construct updated unitary for optimal mu, i.e.
        ! U -> exp(SGN * MUROOT * H) U
        EXPDIAG(:) = EXP( SGN * EVL(:) * MUROOT )
        DO J = 1, NWF
          CWORK(:,J) = V(:,J) * EXPDIAG(J)
        ENDDO
        CWORK = MATMUL(CWORK,CONJG(TRANSPOSE(V)))
#ifdef gammareal
        U = MATMUL(REAL(CWORK,KIND=q),U)
#else
        U = MATMUL(CWORK,U)
#endif

        CURR_COSTFUNC = UMCO_PM_COSTFUNCTION(U)
        DIFF_COSTFUNC = CURR_COSTFUNC - PREV_COSTFUNC
        PREV_COSTFUNC = CURR_COSTFUNC
  
        ! iteration output
        IF(IO%IU0 .GT. 0) THEN
          IF((ITER .GT. 1) .AND. (CGPR_FACTOR .EQ. 0.0_q) ) THEN
            WRITE(IO%IU0,'(I7,4E12.4,E22.14,I7,A)') ITER, GRADIENT_NORM, DIFF_COSTFUNC, CGPR_FACTOR, &
                                   MUROOT, CURR_COSTFUNC, NODE_ME, " (reset CG mixing)"
          ELSE 
            WRITE(IO%IU0,'(I7,4E12.4,E22.14,I7)') ITER, GRADIENT_NORM, DIFF_COSTFUNC, CGPR_FACTOR, &
                                   MUROOT, CURR_COSTFUNC, NODE_ME
          ENDIF
        ENDIF
        IF((ITER .GT. 1) .AND. (CGPR_FACTOR .EQ. 0.0_q) ) THEN
          WRITE(FUNIT,'(I7,4E12.4,E22.14,A)') ITER, GRADIENT_NORM, DIFF_COSTFUNC, CGPR_FACTOR, &
                                MUROOT, CURR_COSTFUNC, " (reset CG mixing)"
        ELSE 
          WRITE(FUNIT,'(I7,4E12.4,E22.14)') ITER, GRADIENT_NORM, DIFF_COSTFUNC, CGPR_FACTOR, &
                                 MUROOT, CURR_COSTFUNC
        ENDIF
  
        ! check if iterations are done
        IF((GRADIENT_NORM .LT. GRADIENT_NORM_STOP) .AND. (ITER .GT. MINITER)) THEN
          IF(IO%IU0 .GT. 0) THEN
            WRITE(IO%IU0,'(A,E22.14)') " converged value of cost function at main mpi rank", CURR_COSTFUNC
            IF(NCPU .GT. 1) WRITE(IO%IU0,'(A)') " waiting until all mpi ranks converged..."
          ENDIF
          WRITE(FUNIT,'(A,E22.14)') " converged value of cost function at main mpi rank", CURR_COSTFUNC
          EXIT
        ENDIF
  
      ENDDO iterations

      PROFILING_STOP('UMCO_ITER')

      CLOSE(FUNIT) ! close the iteration output file

      ! gather the result of cost function from all cores in order
      ! to check which core found the lowest/highest optimum
      ALLOCATE(COSTFUNC_ON_NODE(NCPU))
      COSTFUNC_ON_NODE(:) = 0.0_q
      COSTFUNC_ON_NODE(NODE_ME) = CURR_COSTFUNC
      CALLMPI(M_sum_g(W%WDES%COMM, COSTFUNC_ON_NODE(1), SIZE(COSTFUNC_ON_NODE)))
      DIFF_COSTFUNC = MAXVAL(COSTFUNC_ON_NODE) - MINVAL(COSTFUNC_ON_NODE)
      IF(ABS(DIFF_COSTFUNC) .GT. 10*GRADIENT_NORM_STOP) THEN
        CALL vtutor%alert("The individual mpi ranks converged to different extrema. Compare the converged extrama in the UMCOiter_ISP*_NODE* output files.")
      ENDIF
      IF(SGN .GT. 0) THEN
        WINNER_NODE = MAXLOC(COSTFUNC_ON_NODE,1)
      ELSE
        WINNER_NODE = MINLOC(COSTFUNC_ON_NODE,1)
      ENDIF
      IF(IO%IU0 .GT. 0) THEN
        WRITE(IO%IU0,'(A,I6,A)')   " we pick the result from node", WINNER_NODE, " where the highest/lowest maximum/minimum"
        WRITE(IO%IU0,'(A,E22.14)') " of the cost function was found:", COSTFUNC_ON_NODE(WINNER_NODE)
        WRITE(IO%IU0,'(A,E12.4)')  " Largest deviation of found extrema: " , DIFF_COSTFUNC
      ENDIF
      IF(IO%IU6 .GT. 0) THEN
        WRITE(IO%IU6,'(A,I6,A)')   " we pick the result from node", WINNER_NODE, " where the highest/lowest maximum/minimum"
        WRITE(IO%IU6,'(A,E22.14)') " of the cost function was found:", COSTFUNC_ON_NODE(WINNER_NODE)
        WRITE(IO%IU6,'(A,E12.4)')  " Largest deviation of found extrema: " , DIFF_COSTFUNC
      ENDIF
      DEALLOCATE(COSTFUNC_ON_NODE)

      ! share best unitary with all cores
      IF(NODE_ME .NE. WINNER_NODE) THEN
        U = zero
      ENDIF
      CALLMPI(M_sum_g(W%WDES%COMM, U(1,1), SIZE(U)))

      ! store trafo in global array
      DO IB1=IWF_LOW,IWF_HIGH
         DO IB2=IWF_LOW,IWF_HIGH
            UMCO_UTRAFO(IB2,IB1,IK,ISP) = U(IB2-IWF_LOW+1,IB1-IWF_LOW+1)
         ENDDO
      ENDDO


#ifdef tschaef_verbose
      ! HACK: write out hamiltonian in wannier basis
      IF(IO%IU0 .GT. 0) THEN
        ALLOCATE(HAMIL(NWF,NWF))
        DO J = 1, NWF
          HAMIL(:,J) = W%CELTOT(IWF_LOW:IWF_HIGH,1,ISP) * U(:,J)
        ENDDO
        HAMIL = MATMUL(GCONJG(TRANSPOSE(U)),HAMIL)
        FUNIT=700
        OPEN(unit = FUNIT,file = "hamil_ISP" // str(ISP) // ".dat", FORM='FORMATTED', access='stream', STATUS='REPLACE')
        WRITE(FUNIT,'(I7)') NWF
        DO I = 1, NWF
          DO J = 1, NWF
#ifdef gammareal
            WRITE(FUNIT,'(2I7,E26.16)') I, J, HAMIL(I,J)
#else
            WRITE(FUNIT,'(2I7,2E26.16)') I, J, REAL(HAMIL(I,J), KIND=q), AIMAG(HAMIL(I,J))
#endif
          ENDDO
        ENDDO
        CLOSE(FUNIT)
        DEALLOCATE(HAMIL)
      ENDIF
#endif


      DEALLOCATE(U)
      CALL UMCO_PM_DEALLOC()

      DEALLOCATE(GWORK)
      DEALLOCATE(PRIOR_RIEM_DERV)
      DEALLOCATE(DJDUH)
      DEALLOCATE(H)
      DEALLOCATE(DJDMU)
      DEALLOCATE(V,EVL)
      DEALLOCATE(CWORK)
      DEALLOCATE(EXPDIAG)
      DEALLOCATE(U2)
      DEALLOCATE(MU)
      DEALLOCATE(RWORK)
      DEALLOCATE(IPIV)
      DEALLOCATE(WR,WI)

    ENDDO ! ISPIN

    ! restore
    W%FERTOT = FERTOT
    DEALLOCATE(FERTOT)

    PROFILING_STOP('UMCO_CALC_TRAFO')

  END SUBROUTINE UMCO_CALC_TRAFO






  !********************** INIT_RND_GENERATOR  ******************************
  ! init the random number generator
  !*************************************************************************
  SUBROUTINE INIT_RND_GENERATOR(NODE_ME)
    INTEGER :: NODE_ME
    !local
    INTEGER :: I
    INTEGER :: RND_SEED_SIZE    
    INTEGER :: IRND_INIT
    INTEGER, DIMENSION(:), ALLOCATABLE :: RND_SEED    ! random seed


    CALL SYSTEM_CLOCK(COUNT=IRND_INIT)
    IRND_INIT = IRND_INIT + NODE_ME * 7 ! init nodes individually

    ! initialze random generator with system time
    CALL RANDOM_SEED(SIZE = RND_SEED_SIZE)
    ALLOCATE(RND_SEED(RND_SEED_SIZE))
    DO I = 1, RND_SEED_SIZE
      RND_SEED(I) = IRND_INIT + 97 * I   ! choose the seed values
    END DO
    CALL RANDOM_SEED(PUT = RND_SEED)
    DEALLOCATE(RND_SEED)
  ENDSUBROUTINE INIT_RND_GENERATOR



  !********************** UMCO_RANDOM_UNITARY  ******************************
  ! create a random unitary matrix
  !*************************************************************************
  SUBROUTINE UMCO_RANDOM_UNITARY(U)
    USE constant

    GDEF :: U(NWF,NWF)
    ! local 
    INTEGER              :: I, J
    INTEGER, ALLOCATABLE :: JPVT(:)
    GDEF, ALLOCATABLE    :: TAU(:)
    GDEF, ALLOCATABLE    :: WORK(:)
    GDEF, ALLOCATABLE    :: RWORK(:)
    GDEF, ALLOCATABLE    :: RDIAG(:)
    INTEGER              :: LWORK
    INTEGER              :: INFO
    REAL(q)              :: RNDM(NWF,NWF)
    ! debig
    REAL(q)              :: ALPHA(NWF)

    PROFILING_START('UMCO_RANDOM_UNITARY')

    U = (0.0_q, 0.0_q)
    
    ! generate random matrix
#ifdef gammareal
    CALL RANDOM_NUMBER(RNDM)
    U = U + ( RNDM - 0.5_q ) 
#else
    CALL RANDOM_NUMBER(RNDM)
    U = U + ( RNDM - 0.5_q ) * (1.0_q, 0.0_q)
    CALL RANDOM_NUMBER(RNDM)
    U = U + ( RNDM - 0.5_q ) * (0.0_q, 1.0_q)
#endif

    !
    ! calculate QR decomposition of U and use Q as new U
    !
    ALLOCATE(JPVT(NWF))
    ALLOCATE(TAU(NWF))
    ALLOCATE(WORK(1))
#ifdef gammareal
    CALL DGEQP3(NWF, NWF, U, NWF, JPVT, TAU, WORK, -1, INFO)
#else
    ALLOCATE(RWORK(2*NWF))
    CALL ZGEQP3(NWF, NWF, U, NWF, JPVT, TAU, WORK, -1, RWORK, INFO)
#endif
    LWORK = INT(WORK(1))
    DEALLOCATE(WORK)
    ALLOCATE(WORK(LWORK))
#ifdef gammareal
    CALL DGEQP3(NWF, NWF, U, NWF, JPVT, TAU, WORK, LWORK, INFO)
#else
    CALL ZGEQP3(NWF, NWF, U, NWF, JPVT, TAU, WORK, LWORK, RWORK, INFO)
    DEALLOCATE(RWORK)
#endif

    ALLOCATE(RDIAG(NWF))
    DO I = 1, NWF
      RDIAG(I) = U(I,I) / ABS(U(I,I))
    ENDDO
    
#ifdef gammareal
    CALL DORGQR(NWF, NWF, NWF, U, NWF, TAU, WORK, -1, INFO) 
#else
    CALL ZUNGQR(NWF, NWF, NWF, U, NWF, TAU, WORK, -1, INFO) 
#endif
    LWORK = INT(WORK(1))
    DEALLOCATE(WORK)
    ALLOCATE(WORK(LWORK))
#ifdef gammareal
    CALL DORGQR(NWF, NWF, NWF, U, NWF, TAU, WORK, LWORK, INFO) 
#else
    CALL ZUNGQR(NWF, NWF, NWF, U, NWF, TAU, WORK, LWORK, INFO) 
#endif

    DO J = 1, NWF
      U(:,J) = U(:,J) * RDIAG(J)
    ENDDO

    DEALLOCATE(RDIAG)
    DEALLOCATE(WORK)
    DEALLOCATE(TAU)
    DEALLOCATE(JPVT)

!    U = (0.0_q, 0.0_q)
!    DO J = 1, NWF
!      U(J,J) = 1.0_q
!    ENDDO
!
!#ifdef gammareal
!    ALPHA=0.0_q
!#else
!    CALL RANDOM_NUMBER(ALPHA)
!#endif
!    ALPHA=ALPHA*TPI
!
!    U(:,:) = (0.0_q, 0.0_q)
!    DO I = 1, NWF
!      U(I,I) = EXP(ALPHA(I)*(0.0_q,1.0_q))
!    ENDDO

    PROFILING_STOP('UMCO_RANDOM_UNITARY')
  ENDSUBROUTINE UMCO_RANDOM_UNITARY



#if 0
  ! HM NOTE: I comment this routine for the moment since it causes a circular dependency

  !********************** UMCO_INIT_SCDM_UNITARY  ******************************
  ! initialize the unitary with the SCDM rotation
  !*************************************************************************
  SUBROUTINE UMCO_INIT_SCDM_UNITARY(U, W, P, LATT_CUR, NKPTS, IO)
    USE wave
    USE lattice
    USE poscar_struct_def
    USE pseudo_struct_def
    USE mlwf
    USE wave_high
    USE dfast
    USE fileio
    GDEF            :: U(NWF,NWF)
    TYPE(wavespin)  :: W
    TYPE(potcar)    :: P(:)
    TYPE(latt)      :: LATT_CUR
    INTEGER         :: NKPTS
    TYPE(in_struct) :: IO
    ! local 
    COMPLEX(q), ALLOCATABLE :: T(:,:,:,:)
    INTEGER                 :: IK, ISP
    REAL(q)                 :: VKPT(3, NKPTS)
    ! debug
    TYPE(wavefuna)          :: WA
    TYPE(wavedes1)          :: WDESK

    PROFILING_START('UMCO_INIT_SCDM_UNITARY')

    ! this hast do be an argument later
    VKPT(:,1) = [0.0_q, 0.0_q, 0.0_q];


    IK  = 1
    ISP = 1
    
    CALL CALC_WANNIER_TRAFO_SCDM(W, NWF, T, P, LATT_CUR, NKPTS, VKPT, IO)

#ifdef gammareal
    U(:,:) = REAL(T(1:NWF,1:NWF,IK,ISP), KIND=q)
#else
    U(:,:) = T(1:NWF,1:NWF,IK,ISP)
#endif

    !DO ISP = 1, NWF
    !  DO IK = 1, NWF
    !    WRITE(*,'(2I5,E14.7)') ISP, IK, U(ISP,IK)
    !  ENDDO
    !ENDDO
    !STOP

    IK  = 1
    ISP = 1

    !CALL SETWDES(W%WDES, WDESK, IK)
    !WA = ELEMENTS(W, WDESK, ISP)

    !IF (W%WDES%DO_REDIS) THEN
    !   CALL REDISTRIBUTE_PROJ(WA)
    !   CALL REDISTRIBUTE_PW(WA)
    !ENDIF

    !CALL LINCOM('F', &
    !            WA%CW_RED(:,1:NWF), &
    !            WA%CPROJ_RED(:,1:NWF), &
    !            U(:,:), &
    !            NWF, &
    !            NWF, &
    !            WDESK%NPL_RED, &
    !            WDESK%NPRO_RED, &
    !            WDESK%NRPLWV_RED, &
    !            WDESK%NPROD_RED, &
    !            NWF, &
    !            WA%CW_RED(:,1:NWF), &
    !            WA%CPROJ_RED(:,1:NWF))

    !IF (W%WDES%DO_REDIS) THEN
    !   CALL REDISTRIBUTE_PROJ(WA)
    !   CALL REDISTRIBUTE_PW(WA)
    !ENDIF

    !CALL OUTWAV(IO, W%WDES, W, LATT_CUR, 0.0_q, 'scd')

    PROFILING_STOP('UMCO_INIT_SCDM_UNITARY')


  END SUBROUTINE UMCO_INIT_SCDM_UNITARY
#endif




  !********************** UMCO_DIAG_H  ******************************
  ! diagonalizes the skew-hermitian matrix H and stores the
  ! the unitary diagoanlization in W and the purely imaginary
  ! eigenvalues in EVL
  !******************************************************************
  SUBROUTINE UMCO_DIAG_H(H, V, EVL)
    GDEF, intent(in)        :: H(NWF,NWF)
    COMPLEX(q), intent(out) :: V(NWF,NWF)
    COMPLEX(q), intent(out) :: EVL(NWF)
    ! local
    COMPLEX(q), ALLOCATABLE :: M(:,:)   ! hermitian matrix, later eigenvectors
    REAL(q), ALLOCATABLE    :: REVL(:)  ! real eigenvalues of M
    COMPLEX(q), ALLOCATABLE :: WORK(:)  ! work array
    REAL(q), ALLOCATABLE    :: RWORK(:) ! work array
    INTEGER                 :: INFO
    INTEGER                 :: LWORK    ! size of work array
    COMPLEX(q)              :: WORK_DIM
    ! debug
    INTEGER :: I,J

    PROFILING_START('UMCO_DIAG_H')

    ALLOCATE(M(NWF,NWF), REVL(NWF))

    ! create hermitian M from skew hermitian H via M=iH
    M=(0.0_q,1.0_q)*H

    ALLOCATE(RWORK(3*NWF-2))
    ! query optimal work arrays for diagonalization
    CALL ZHEEV('V', 'L', NWF, M, NWF, REVL,  &
               WORK_DIM, -1, RWORK, INFO)
    LWORK =INT(REAL(WORK_DIM, KIND=q))

    ! perform diagonalization
    ALLOCATE(WORK(LWORK))
    CALL ZHEEV('V', 'L', NWF, M, NWF, REVL,  &
               WORK, LWORK, RWORK, INFO)

    !DO J=1,NWF
    !  DO I=1,NWF
    !    V(I,J) = M(I,J)
    !    WRITE(*,*) I,J, V(I,J)
    !  ENDDO
    !ENDDO
    V=M
    EVL=(0.0_q, -1.0_q)*REVL

    DEALLOCATE(WORK)
    DEALLOCATE(RWORK)
    DEALLOCATE(REVL, M)

    PROFILING_STOP('UMCO_DIAG_H')

  END SUBROUTINE UMCO_DIAG_H



  !******************** UMCO_CALC_DJDMU  ************************************
  ! calculates the derivative of the composition of the cost function 
  ! with a rotated unitary, J(exp(SGN*mu*H)), with respect to mu (DJDMU).
  ! As arguments we need the euclidean derivative of J wrt U^dagger (DJDUH),
  ! the rotated unitary U2=exp(SGN*mu*H)U and the descent direction H
  !**************************************************************************
  SUBROUTINE UMCO_CALC_DJDMU(DJDUH,U2,H, DJDMU)
    GDEF, intent(in)     :: DJDUH(NWF,NWF)
    GDEF, intent(in)     :: U2(NWF,NWF)
    GDEF, intent(in)     :: H(NWF,NWF)
    REAL(q), intent(out) :: DJDMU
    ! local
    GDEF, ALLOCATABLE :: GWORK(:,:)
    INTEGER           :: I

    PROFILING_START('UMCO_CALC_DJDMU')

    ALLOCATE(GWORK(NWF,NWF))

    GWORK = MATMUL(H,U2)
    GWORK = GCONJG(TRANSPOSE(GWORK))
    !GWORK = MATMUL(U2,GWORK)
    GWORK = MATMUL(DJDUH,GWORK)

    DJDMU = 0.0_q

    ! calculate trace 
    DO I=1,NWF
      DJDMU = DJDMU + 2 * SGN * REAL(GWORK(I,I), KIND=q)
    ENDDO

    DEALLOCATE(GWORK)

    PROFILING_STOP('UMCO_CALC_DJDMU')
  END SUBROUTINE UMCO_CALC_DJDMU



  !********************** UMCO_CGPR_FACTOR  ******************************
  ! calc the Polak-Ribière-Polyak update factor for conjugate gradient algorithm
  !*************************************************************************
  FUNCTION UMCO_CGPR_FACTOR(CURRENT_RIEM_DERV, PRIOR_RIEM_DERV) RESULT(CGPR_FACTOR)
    GDEF    :: CURRENT_RIEM_DERV(NWF,NWF)
    GDEF    :: PRIOR_RIEM_DERV(NWF,NWF)
    REAL(q) :: CGPR_FACTOR

    CGPR_FACTOR = REAL( &
      GFROBENIUS_PRODUCT(CURRENT_RIEM_DERV,CURRENT_RIEM_DERV-PRIOR_RIEM_DERV) &
      / GFROBENIUS_PRODUCT(PRIOR_RIEM_DERV,PRIOR_RIEM_DERV), KIND=q)

  END FUNCTION UMCO_CGPR_FACTOR



  !********************** GFROBENIUS_PRODUCT  ******************************
  ! calc the Polak-Ribière-Polyak update factor for conjugate gradient algorithm
  !*************************************************************************
  FUNCTION GFROBENIUS_PRODUCT(A,B) RESULT(S)
    GDEF :: A(:,:)
    GDEF :: B(:,:)
    GDEF :: S
    ! local
    INTEGER I,J
    GDEF, ALLOCATABLE :: C(:,:)

    PROFILING_START("GFROBENIUS_PRODUCT")

    S = zero

    IF( (SIZE(A,2) .NE. SIZE(B,1)) .OR. (SIZE(A,1) .NE. SIZE(B,2)) ) THEN
      WRITE(*,*) " GFROBENIUS_PRODUCT: wrong dimension"
    ENDIF

    ALLOCATE(C(SIZE(B,2),SIZE(B,1)))

    C = GCONJG(B)

    DO J = 1, SIZE(A,2)
      DO I = 1, SIZE(A,1)
        S = S + A(I,J)*C(I,J)
      ENDDO
    ENDDO
    PROFILING_STOP("GFROBENIUS_PRODUCT")
  END FUNCTION GFROBENIUS_PRODUCT



END MODULE umco
