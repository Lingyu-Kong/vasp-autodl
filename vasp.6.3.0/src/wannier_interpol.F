#include "symbol.inc"

!************************ MODULE WAVE_INTERPOLATE **********************
!>
!> @brief
!> Interpolate the band structure to an arbitrary k-point grid
!> using a localized basis (normally wannier)
!> This routine was writted by Henrique Miranda mostly based on
!> a previous implementation present in VASP done by Gerald Jordan
!>
!***********************************************************************

MODULE wannier_interpolation
  USE prec
  USE lr_helper
  USE twoelectron4o
  USE mkpoints
  USE full_kpoints
  USE wave_high
  USE lattice
  USE poscar
  USE mpimy
  USE mgrid
  USE nonl_high
  USE base
  USE pseudo
  USE kpoints_change
  USE hamil_struct_def
  USE msymmetry
  USE fileio
  IMPLICIT NONE

  !> base type to perform wannier interpolation
  TYPE wannier_interpolator
     TYPE (kpoints_struct) :: KPOINTS_COARSE_FULL !< original k-points in the FBZ
     TYPE (kpoints_struct) :: KPOINTS_INTER       !< kpoint array of IBZ of interpolated grid
     TYPE (wavedes) :: WDES_COARSE_FULL           !< descriptor for initial wave function in the FBZ
     TYPE (wavedes) :: WDES_INTER                 !< descriptor for interpolated wave function
     INTEGER :: ISPIN = 0                         !< number of spins
     INTEGER :: NW = 0                            !< number of wannier functions
     GDEF, ALLOCATABLE :: U(:,:,:,:)              !< wannier rotation matrices
     GDEF, ALLOCATABLE :: HAM(:,:,:,:)            !< Hamiltonian on coarse grid
     GDEF, ALLOCATABLE :: CDER(:,:,:,:,:)         !< Optical matrix elements on coarse grid in the Wannier basis
     INTEGER,ALLOCATABLE :: NEQUIV(:)             !< Mapping between FBZ and IBZ of coarse grid
  END TYPE wannier_interpolator

  ! Routines
  PUBLIC :: WANNIER_INTERPOL_READER
  PUBLIC :: WANNIER_INTERPOL_SETUP_KPOINTS
  PUBLIC :: WANNIER_INTERPOL_FROM_WANPROJ
  PUBLIC :: WANNIER_INTERPOL_FROM_MLWF
  PUBLIC :: WANNIER_INTERPOL_SETUP_HAM
  PUBLIC :: WANNIER_INTERPOL_SETUP_CDER
  PUBLIC :: WANNIER_INTERPOLATE
  PUBLIC :: WANNIER_INTERPOL_FREE
  ! Variables
  PUBLIC :: LWANNIERINTERPOL

  LOGICAL, SAVE :: LINTPOL_KPATH !< Compute eigenvalues and eigenvectors along a path
  LOGICAL, SAVE :: LINTPOL_OPTICS !< Compute optical properties
  LOGICAL, SAVE :: LINTPOL_VEL !< Compute the band-velocities
  LOGICAL, SAVE :: LINTPOL_DOS !< Compute the density of states
  LOGICAL, SAVE :: LWANNIERINTERPOL=.FALSE.!< Whether the wannier interpolator routines should be called

CONTAINS

  !> @brief Read and store input variables related to the wannier interpolation
  SUBROUTINE WANNIER_INTERPOL_READER(IO)
    USE base
    USE crpa, ONLY: LCRPAPLOT
    USE reader_tags

    TYPE (in_struct) IO
    ! local variables
    INTEGER IERR
    LOGICAL LOPEN

    CALL OPEN_INCAR_IF_FOUND(IO%IU5, LOPEN)

    ! interpolate bandstructure on a given k-grid
    INQUIRE(FILE="KPOINTS_WAN",EXIST=LINTPOL_KPATH)
    CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LINTPOL_KPATH', LINTPOL_KPATH, IERR, WRITEXMLINCAR)
    IF( LCRPAPLOT ) THEN
        IF (IO%IU0>=0) WRITE(IO%IU0,*)'Wannier interpolation of k-path set by LCRPAPLOT'
        LINTPOL_KPATH = .TRUE.
    ENDIF 

    ! Calculate optical properties
    LINTPOL_OPTICS=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LINTPOL_OPTICS', LINTPOL_OPTICS, IERR, WRITEXMLINCAR)

    ! Calculate conductivity
    LINTPOL_VEL=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LINTPOL_VEL', LINTPOL_VEL, IERR, WRITEXMLINCAR)

    ! Calculate density of states
    LINTPOL_DOS=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LINTPOL_DOS', LINTPOL_DOS, IERR, WRITEXMLINCAR)

    LWANNIERINTERPOL=LWANNIERINTERPOL .OR. LINTPOL_KPATH .OR. LINTPOL_OPTICS .OR. LINTPOL_VEL .OR. LINTPOL_DOS

    CALL CLOSE_INCAR_IF_FOUND(IO%IU5)

    IF (LWANNIERINTERPOL) THEN
#ifdef gammareal
       CALL vtutor%error("WANNIER_INTERPOL_READER: ERROR: Wannier interpolation is not available for &
          &gamma-only.")
#endif
    ENDIF

  END SUBROUTINE WANNIER_INTERPOL_READER

  !************* SUBROUTINE WANNIER_INTERPOL_SETUP_KPOINTS ****************
  ! 
  !> @brief Setup k-points to be used in Wannier interpolation
  ! 
  !************************************************************************
  SUBROUTINE WANNIER_INTERPOL_SETUP_KPOINTS(WANN_INTERP,KPOINTS,LATT_CUR,SYMM,W)
    USE prec
    USE wave_high
    USE full_kpoints
    USE tutor, ONLY: vtutor
    TYPE (wannier_interpolator) WANN_INTERP
    TYPE (kpoints_struct) KPOINTS
    TYPE (latt)        LATT_CUR
    TYPE (symmetry) SYMM
    TYPE (wavespin) :: W
    !    local
    INTEGER :: NKPTS_FULL

    CALL CHECK_FULL_KPOINTS

    ! initialize KPOINTS_INTER where interpolated quantities will be computed
    CALL COPY_KPOINTS(KPOINTS,WANN_INTERP%KPOINTS_INTER)
#ifdef oldsym
    CALL RD_KPOINTS(WANN_INTERP%KPOINTS_INTER, LATT_CUR, &
         SYMM%ISYM>=0.AND..NOT. W%WDES%LSORBIT.AND..NOT. W%WDES%LSPIRAL, &
         SYMM%ISYM<0, -1, -1, FILE_NAME="KPOINTS_WAN")
#else
    CALL RD_KPOINTS(WANN_INTERP%KPOINTS_INTER, LATT_CUR, &
         SYMM%ISYM>=0.AND..NOT. W%WDES%LNONCOLLINEAR, SYMM%ISYM<0, &
         -1, -1, FILE_NAME="KPOINTS_WAN")
#endif

    ! initialize wavefunction descriptor for interpolation
    WANN_INTERP%WDES_INTER=W%WDES
    WANN_INTERP%WDES_INTER%NKPTS= WANN_INTERP%KPOINTS_INTER%NKPTS
    WANN_INTERP%WDES_INTER%VKPT =>WANN_INTERP%KPOINTS_INTER%VKPT
    WANN_INTERP%WDES_INTER%WTKPT=>WANN_INTERP%KPOINTS_INTER%WTKPT

    ! intialize FBZ for coarse grid 
    NKPTS_FULL = KPOINTS_FULL%NKPTS
    ALLOCATE(WANN_INTERP%NEQUIV(NKPTS_FULL))
    WANN_INTERP%NEQUIV = KPOINTS_FULL%NEQUIV

    ! set KPOINTS_COARSE_FULL
    CALL COPY_KPOINTS(KPOINTS,WANN_INTERP%KPOINTS_COARSE_FULL)
    DEALLOCATE(WANN_INTERP%KPOINTS_COARSE_FULL%VKPT)
    DEALLOCATE(WANN_INTERP%KPOINTS_COARSE_FULL%WTKPT)
    ALLOCATE(WANN_INTERP%KPOINTS_COARSE_FULL%VKPT(3,NKPTS_FULL))
    ALLOCATE(WANN_INTERP%KPOINTS_COARSE_FULL%WTKPT(NKPTS_FULL))
    WANN_INTERP%KPOINTS_COARSE_FULL%NKPTS = NKPTS_FULL
    WANN_INTERP%KPOINTS_COARSE_FULL%VKPT  = KPOINTS_FULL%VKPT
    WANN_INTERP%KPOINTS_COARSE_FULL%WTKPT = 1.0/NKPTS_FULL

    ! set wdes
    WANN_INTERP%WDES_COARSE_FULL = W%WDES
    WANN_INTERP%WDES_COARSE_FULL%NKPTS= WANN_INTERP%KPOINTS_COARSE_FULL%NKPTS
    WANN_INTERP%WDES_COARSE_FULL%VKPT =>WANN_INTERP%KPOINTS_COARSE_FULL%VKPT
    WANN_INTERP%WDES_COARSE_FULL%WTKPT=>WANN_INTERP%KPOINTS_COARSE_FULL%WTKPT

  END SUBROUTINE WANNIER_INTERPOL_SETUP_KPOINTS


  !************* SUBROUTINE WANNIER_INTERPOL_FROM_WANPROJ ****************
  ! 
  !> @brief Read transformation matrices for Wannier orbitals from WANPROJ file
  !> and store in WANN_INTERP%U
  ! 
  !************************************************************************
  SUBROUTINE WANNIER_INTERPOL_FROM_WANPROJ(WANN_INTERP,IO)
    USE wave_high
    USE mlwf, ONLY: MLWF_type, READ_WANPROJ_FILE
    USE locproj, ONLY: WANPROJFILE
    TYPE (wannier_interpolator) :: WANN_INTERP
    TYPE (wavespin),TARGET :: W
    TYPE (in_struct)   IO
    ! local
    LOGICAL :: LWANPROJ
    INTEGER :: NW, NKPTS

    ! check if we can read the rotation matrices
    CALL READ_WANPROJ_FILE(WANPROJFILE,WANN_INTERP%WDES_COARSE_FULL,IO,LWANPROJ,WANN_INTERP%NW,WANN_INTERP%U)

    IF(.NOT.LWANPROJ)THEN
       CALL vtutor%error( WANPROJFILE//' file not found' )
    END IF

    ! initialize WDES_INTER to fit KPOINTS_INTER
    NW = WANN_INTERP%NW
    NKPTS = WANN_INTERP%KPOINTS_INTER%NKPTS
    WANN_INTERP%WDES_INTER%NB_TOT = NW
    WANN_INTERP%WDES_INTER%NBANDS = NW/WANN_INTERP%WDES_INTER%NB_PAR
#ifdef MPI
    ! Normally VASP ensures that the number of bands is divisible by the number of nodes
    ! but this is not possible here. We cannot add more bands to have a number divisible
    ! by the number of CPUS so we just distribute the leftover bands amongst the CPUs.
    ! Otherwise the DOS routine will fail to determine the fermi energy
    IF (WANN_INTERP%WDES_INTER%COMM%NODE_ME <= MOD(NW,WANN_INTERP%WDES_INTER%NB_PAR)) THEN
        WANN_INTERP%WDES_INTER%NBANDS = WANN_INTERP%WDES_INTER%NBANDS + 1
    ENDIF
#endif
    NULLIFY(WANN_INTERP%WDES_INTER%NB_TOTK)
    ALLOCATE(WANN_INTERP%WDES_INTER%NB_TOTK(NKPTS,WANN_INTERP%ISPIN))
    WANN_INTERP%WDES_INTER%NB_TOTK = NW

  END SUBROUTINE WANNIER_INTERPOL_FROM_WANPROJ


  !************* SUBROUTINE WANNIER_INTERPOL_FROM_MLWF ****************
  !
  !> @brief Read transformation matrices for Wannier orbitals from MLWF structure
  !
  !************************************************************************
  SUBROUTINE WANNIER_INTERPOL_FROM_MLWF(WANN_INTERP,MLWF_STRUCT,IO)
    USE wave_high
    USE mlwf, ONLY: MLWF_type, MLWF_GET_U
    USE locproj, ONLY: WANPROJFILE
    TYPE (wannier_interpolator) :: WANN_INTERP
    TYPE (mlwf_type) :: MLWF_STRUCT
    TYPE (in_struct)   IO
    ! local
    LOGICAL :: LWANPROJ
    INTEGER :: NW, NKPTS

    CALL MLWF_GET_U(MLWF_STRUCT,MLWF_STRUCT%NKPTS,MLWF_STRUCT%kpt_latt,WANN_INTERP%U)
    NW = MLWF_STRUCT%NUM_WANN
    WANN_INTERP%NW = MLWF_STRUCT%NUM_WANN

    ! initialize WDES_INTER to fit KPOINTS_INTER
    NKPTS = WANN_INTERP%KPOINTS_INTER%NKPTS
    WANN_INTERP%WDES_INTER%NB_TOT = NW
    WANN_INTERP%WDES_INTER%NBANDS = NW/WANN_INTERP%WDES_INTER%NB_PAR
#ifdef MPI
    ! Normally VASP ensures that the number of bands is divisible by the number of nodes
    ! but this is not possible here. We cannot add more bands to have a number divisible
    ! by the number of CPUS so we just distribute the leftover bands amongst the CPUs.
    ! Otherwise the DOS routine will fail to determine the fermi energy
    IF (WANN_INTERP%WDES_INTER%COMM%NODE_ME <= MOD(NW,WANN_INTERP%WDES_INTER%NB_PAR)) THEN
        WANN_INTERP%WDES_INTER%NBANDS = WANN_INTERP%WDES_INTER%NBANDS + 1
    ENDIF
#endif
    NULLIFY(WANN_INTERP%WDES_INTER%NB_TOTK)
    ALLOCATE(WANN_INTERP%WDES_INTER%NB_TOTK(NKPTS,WANN_INTERP%ISPIN))
    WANN_INTERP%WDES_INTER%NB_TOTK = NW

  END SUBROUTINE WANNIER_INTERPOL_FROM_MLWF


  !************* SUBROUTINE WANNIER_INTERPOL_SETUP_HAM ********************
  !
  !> @brief Setup Wannier Hamiltonian for interpolation
  !
  !************************************************************************
  SUBROUTINE WANNIER_INTERPOL_SETUP_HAM(WANN_INTERP,W,IO)
    USE wave_high
    USE xi
    USE mlwf
    TYPE (wannier_interpolator) WANN_INTERP
    TYPE (wavespin), TARGET :: W
    TYPE (in_struct)   IO

    ! Apply scissor operator
    IF(SCISSOR/=0)THEN
       CALL APPLY_SCISSOR(W,SCISSOR)
    END IF

    ! Compute Hamiltonian in the Wannier representation
    IF (IO%IU0>=0) WRITE(IO%IU0,*) 'setup Wannier Hamiltonian'
    ALLOCATE(WANN_INTERP%HAM(WANN_INTERP%NW,WANN_INTERP%NW,WANN_INTERP%KPOINTS_COARSE_FULL%NKPTS,W%WDES%ISPIN))
    CALL MLWF_ROTATE_HAM(WANN_INTERP%WDES_COARSE_FULL,W,WANN_INTERP%NW,WANN_INTERP%HAM,WANN_INTERP%NEQUIV,WANN_INTERP%U)

  END SUBROUTINE WANNIER_INTERPOL_SETUP_HAM

  !************* SUBROUTINE WANNIER_INTERPOL_SETUP_CDER ********************
  !
  !> @brief Setup interpolation of optical matrix elements
  !> CDER is computed in the FBZ, rotated to the Wannier basis and 
  !> stored in WANN_INTERP%CDER
  !
  !************************************************************************
  SUBROUTINE WANNIER_INTERPOL_SETUP_CDER(WANN_INTERP,W,KPOINTS,GRID,GRIDC,GRIDUS,C_TO_US,IRDMAX,LMDIM,T_INFO,&
                          NONLR_S,NONL_S,P,SV,CQIJ,CDIJ, &
                          LATT_INI,LATT_CUR,SYMM,INFO,IO)
    USE wave_high
    USE xi
    USE mlwf
    TYPE (wannier_interpolator) WANN_INTERP
    TYPE (wavespin) W
    TYPE (latt) :: LATT_INI
    TYPE (latt) :: LATT_CUR
    TYPE (kpoints_struct) :: KPOINTS
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
    OVERLAP :: CQIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    OVERLAP :: CDIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    RGRID   :: SV(DIMREAL(GRID%MPLWV),W%WDES%NCDIJ)

    IF (.NOT.(LINTPOL_VEL.OR.LINTPOL_OPTICS)) RETURN

    CALL MLWF_ROTATE_CDER(W,KPOINTS,WANN_INTERP%KPOINTS_COARSE_FULL,WANN_INTERP%U,&
                          GRID,GRIDC,GRIDUS,C_TO_US,IRDMAX,LMDIM,T_INFO,&
                          NONLR_S,NONL_S,P,SV,CQIJ,CDIJ, &
                          LATT_INI,LATT_CUR,WANN_INTERP%NW,WANN_INTERP%CDER,SYMM,INFO,IO)

  END SUBROUTINE WANNIER_INTERPOL_SETUP_CDER

  !************* SUBROUTINE HAM_DIAG **************************************
  ! 
  !> @brief Diagonalize Hamiltonian matrix
  !
  !************************************************************************
  SUBROUTINE HAM_DIAG(NBANDS,NKPTS,NSPIN,HAM,EIG)
    USE string, ONLY: str
    USE tutor, ONLY: vtutor
    INTEGER,INTENT(IN) :: NBANDS,NKPTS,NSPIN
    GDEF,INTENT(IN) :: HAM(NBANDS,NBANDS,NKPTS,NSPIN)
    REAL(q),INTENT(IN) :: EIG(NBANDS,NKPTS,NSPIN)
    ! local variables
    INTEGER :: ISP,IK,IFAIL
    ! work arrays for ZHEEV (blocksize times number of bands)
    INTEGER, PARAMETER :: LWORK=32
    GDEF       CWRK(LWORK*NBANDS)
    REAL(q)    RWORK(3*NBANDS)

    spin: DO ISP=1,NSPIN
       kpoint: DO IK=1,NKPTS

#ifdef gammareal
          CALL DSYEV('V','U',NBANDS,HAM(:,:,IK,ISP),NBANDS,EIG(:,IK,ISP),CWRK,LWORK*NBANDS,RWORK,IFAIL)
#else 
          CALL ZHEEV('V','U',NBANDS,HAM(:,:,IK,ISP),NBANDS,EIG(:,IK,ISP),CWRK,LWORK*NBANDS,RWORK,IFAIL)
#endif
          IF(IFAIL/=0) CALL vtutor%error("HAM_DIAG: error in ZHEEV: " // str(IFAIL))

       END DO kpoint
    END DO spin

  END SUBROUTINE HAM_DIAG

  !************* SUBROUTINE WANNIER_INTERPOLATE_KPATH *********************
  !
  !> @brief Interpolate energies along a path of k-points
  !
  !************************************************************************
  SUBROUTINE WANNIER_INTERPOLATE(WANN_INTERP,W,KPOINTS,SYMM,LATT_INI,LATT_CUR,NEDOS,T_INFO,INFO,IO)
    USE prec
    USE ini
    USE wave_high
    USE lattice
    USE poscar
    USE nonl_high
    USE base
    USE pseudo
    USE constant
    USE choleski
    USE subrot
    USE mlwf
    USE tutor, ONLY: vtutor
    USE mlrf_main, ONLY: OMEGAMAX_OPTIC, RTIME
    USE mlr_optic, ONLY: MAX_ENERGY_UNOCCUPIED, WPLASMON, CONDUCTIVITY_ENERGY_RESOLVED, &
                         LRF_EPSILON_COMPUTE, LR_OPTIC_OUTPUT, NBANDS_CDER
#ifdef VASP_HDF5
    USE vhdf5
#endif

    TYPE (wannier_interpolator) :: WANN_INTERP
    TYPE (wavespin) :: W
    TYPE (kpoints_struct) KPOINTS
    TYPE (latt),INTENT(IN) :: LATT_INI
    TYPE (latt),INTENT(IN) :: LATT_CUR
    TYPE (type_info)    T_INFO
    TYPE (symmetry) SYMM
    INTEGER :: NEDOS !< Number of points for the computation of DOS 
    TYPE (info_struct) INFO
    TYPE (in_struct)   IO
    ! local variables
    INTEGER :: I,IK,ISP,NW,IW,IDIR
    INTEGER :: NKPTS
    REAL(q) :: EDIFF
    GDEF, ALLOCATABLE :: HAM_KPATH(:,:,:,:)
    COMPLEX(q), ALLOCATABLE :: HAMR(:,:,:,:,:,:)
    !COMPLEX(q), ALLOCATABLE :: DHDK_KPATH(:,:,:,:,:)
    REAL(q), ALLOCATABLE :: EIG_KPATH(:,:,:)
    REAL(q), ALLOCATABLE :: VEL_KPATH(:,:,:,:)
    GDEF, ALLOCATABLE :: CHAM_KPATH(:,:,:,:)
    GDEFS, ALLOCATABLE :: CHAMS_KPATH(:,:,:,:,:)
    GDEFS :: C(WANN_INTERP%NW,WANN_INTERP%NW), D(WANN_INTERP%NW,WANN_INTERP%NW), U(WANN_INTERP%NW,WANN_INTERP%NW)

    !> Wavefunction container for interpolated wavefunction
    TYPE (wavespin) :: W_INTER

    ! Final computations
    REAL(q) :: EMAX, EMAX_COND
    REAL(q),ALLOCATABLE :: COND_ENERGY(:,:,:,:), EDOS(:)
    COMPLEX(q), ALLOCATABLE  :: EPSDD(:,:,:)
    COMPLEX(q), ALLOCATABLE  :: EPSJJ(:,:,:)
    REAL(q)    :: CON(3,3), BIMAG(3,3), DELTAE, EFERMI, ENTROPY
    COMPLEX(q) :: WPLASMA_INTER(3,3)
    REAL(q) :: PAR(1,1,1,1,W%WDES%NCDIJ)
    REAL(q) :: DOSPAR(1,1,1,W%WDES%NCDIJ)
    REAL(q) :: DOS(NEDOS,W%WDES%ISPIN),DOSI(NEDOS,W%WDES%ISPIN)
    INTEGER :: IB1,IB2
    REAL(q) :: E1,E2,OCC
    CHARACTER(LEN=:), ALLOCATABLE :: SUBGROUP

    EMAX = 0.0_q
    EMAX_COND = 0.0_q
    DELTAE = 0.0_q
    EFERMI = 0.0_q
    ENTROPY = 0.0_q
    E1 = 0.0_q
    E2 = 0.0_q
    OCC = 0.0_q
    CON = 0.0_q
    BIMAG = 0.0_q
    WPLASMA_INTER = (0.0_q, 0.0_q)
    PAR = 0.0_q
    DOSPAR = 0.0_q
    DOS = 0.0_q
    DOSI = 0.0_q

    ! Allocate a container for interpolated eigenvalues
    CALL ALLOCW_NOPLANEWAVE(WANN_INTERP%WDES_INTER, W_INTER)

    ! Allocate memory
    NW = WANN_INTERP%NW
    NKPTS = WANN_INTERP%KPOINTS_INTER%NKPTS
    ALLOCATE(HAM_KPATH(NW,NW,NKPTS,W%WDES%ISPIN))
    ALLOCATE(EIG_KPATH(NW,NKPTS,W%WDES%ISPIN))
    ALLOCATE(CHAM_KPATH(NW,NW,NKPTS,W%WDES%ISPIN))
    ALLOCATE(CHAMS_KPATH(NW,NW,NKPTS,W%WDES%ISPIN,3))
    ALLOCATE(VEL_KPATH(NW,NKPTS,W%WDES%ISPIN,3))

    ! HAMILTONIAN INTERPOLATION
    CALL FOURIER_INTERPOL_KTOR(WANN_INTERP%KPOINTS_COARSE_FULL,WANN_INTERP%HAM,HAMR)
    CALL FOURIER_INTERPOL_RTOK(HAMR,WANN_INTERP%KPOINTS_INTER%VKPT,HAM_KPATH,LATT_CUR)
    DEALLOCATE(HAMR)

    ! HAMILTONIAN DIAGONALIZATION
    CALL HAM_DIAG(WANN_INTERP%NW,NKPTS,W%WDES%ISPIN,HAM_KPATH,EIG_KPATH)

    DO ISP=1,W%WDES%ISPIN
       DO IK=1,NKPTS
          W_INTER%CELTOT(:NW,IK,ISP) = EIG_KPATH(:NW,IK,ISP)
       ENDDO
    ENDDO    

#if 0
    !-----------------------------------------------------
    ! Transform dH/dk to cartesian coordinates
    !
    ! d/dk_i = dk_m/dk_i d/dk_m   where i=cartesian, m=direct coordinate
    ! \vec{k} = k_m \vec{b}^m
    ! 2pi k_m = \vec{a}^m . \vec{k} = a^m_i k_i
    ! Therefore: dk_m/dk_i = a^m_i/2pi
    ALLOCATE(DHDK_INT_CART(NW,NW,3))
    DO ISP=1,W%WDES%ISPIN
       DO IK=1,NKPTS
          DHDK_INT_CART=0
          DO IDIR_CAR=1,3
             DO IDIR_RED=1,3
                DHDK_INT_CART(:,:,IDIR_CAR)=DHDK_INT_CART(:,:,IDIR_CAR)+&
                                            DHDK_KPATH(:,:,IK,ISP,IDIR_RED)*LATT_CUR%A(IDIR_CAR,IDIR_RED)/TPI
             END DO
          END DO
          DHDK_KPATH(:,:,IK,ISP,:) = DHDK_INT_CART
       END DO
    END DO
    DEALLOCATE(DHDK_INT_CART)
#endif

    IF (LINTPOL_OPTICS.OR.LINTPOL_VEL) THEN
       DO IDIR=1,3
          IF (IO%IU0>=0) WRITE(IO%IU0,*) 'Interpolating direction', IDIR
          ! interpolation of the optical matrix elements
          CALL FOURIER_INTERPOL_KTOR(WANN_INTERP%KPOINTS_COARSE_FULL,WANN_INTERP%CDER(:,:,:,:,IDIR),HAMR)
          CALL FOURIER_INTERPOL_RTOK(HAMR,WANN_INTERP%KPOINTS_INTER%VKPT,CHAM_KPATH,LATT_CUR)
          DEALLOCATE(HAMR)
          ! compute sandwich with wannier states to transform to Hamiltonian gauge again
          DO ISP=1,W%WDES%ISPIN
             DO IK=1,NKPTS
                ! transform A
                U = HAM_KPATH(:,:,IK,ISP)
                D = CHAM_KPATH(:,:,IK,ISP)
                CALL CGEMM('N', 'N', NW, NW, NW, CMPLX(1.0,0.0), D, NW, U, NW, CMPLX(0.0,0.0), C, NW)
                CALL CGEMM('C', 'N', NW, NW, NW, CMPLX(1.0,0.0), U, NW, C, NW, CMPLX(0.0,0.0), D, NW)
                ! Save velocities
                DO IW=1,NW
                   VEL_KPATH(IW,IK,ISP,IDIR) = REAL(D(IW,IW),q)
                ENDDO
                DO IB2=1,NW
                   E2 = REAL(EIG_KPATH(IB2,IK,ISP),q)
                   DO IB1=1,NW
                      E1 = REAL(EIG_KPATH(IB1,IK,ISP),q)
                      EDIFF = E2-E1
                      IF (ABS(EDIFF)>1E-10_q) THEN
                         CHAMS_KPATH(IB1,IB2,IK,ISP,IDIR) = D(IB1,IB2)/EDIFF
                      ELSE
                         CHAMS_KPATH(IB1,IB2,IK,ISP,IDIR) = 0
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF


    ! Compute IP optical spectra and conductivity on dense grid
    !-----------------------------------------------------------------
    ! This one is very nasty
    NBANDS_CDER = WANN_INTERP%WDES_INTER%NB_TOT

    ! Set energy bounds for DOS and conductivity calculation
    EMAX=MAX_ENERGY_UNOCCUPIED(W%WDES,W)*1.2; IF (OMEGAMAX_OPTIC/=-1) EMAX=OMEGAMAX_OPTIC
    EMAX_COND=5; IF (OMEGAMAX_OPTIC/=-1) EMAX_COND=OMEGAMAX_OPTIC
    DELTAE=EMAX/(NEDOS-1)

    ! Related to eigenvalues only
    !-----------------------------------------------
    IF (LINTPOL_DOS) THEN
       OCC = 2._q; IF (W%WDES%ISPIN==2.OR.W%WDES%NRSPINORS==2) OCC = 1._q
       IF (INFO%NELECT>=WANN_INTERP%NW*OCC) CALL vtutor%error('Not enough wannier states to compute the fermi energy')
       CALL DENSTA(IO%IU0, IO%IU6, WANN_INTERP%WDES_INTER, W_INTER, WANN_INTERP%KPOINTS_INTER, &
            INFO%NELECT, INFO%NUP_DOWN, ENTROPY, EFERMI, WANN_INTERP%KPOINTS_INTER%SIGMA, &
            .TRUE., NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)

       CALL XML_DOS(EFERMI, WANN_INTERP%KPOINTS_INTER%EMIN, WANN_INTERP%KPOINTS_INTER%EMAX, .FALSE., &
            DOS, DOSI, DOSPAR, NEDOS, 1, 1, WANN_INTERP%WDES_INTER%NCDIJ, &
            comment='kpoints_wan')
    ENDIF

    CALL XML_TAG("eigenvalues_interp",comment='wannier')
    CALL XML_TAG("kpoints")
    CALL XML_KPOINTS_LIST(WANN_INTERP%KPOINTS_INTER%VKPT, WANN_INTERP%KPOINTS_INTER%WTKPT)
    CALL XML_CLOSE_TAG("kpoints")
    CALL XML_EIGENVALUES_EXT(REAL(W_INTER%CELTOT),1, NW, WANN_INTERP%KPOINTS_INTER%NKPTS, W%WDES%ISPIN)
    CALL XML_CLOSE_TAG("eigenvalues_interp")

#ifdef VASP_HDF5
    CALL VH5_WRITE_DOS(IH5OUTFILEID, WANN_INTERP%WDES_INTER, WANN_INTERP%KPOINTS_INTER, DOS, DOSI, DOSPAR, EFERMI, T_INFO%NIONP, -1, SUBGROUP="electron_dos_kpoints_wan")
    CALL VH5_WRITE_EIGENVAL(IH5OUTFILEID, WANN_INTERP%WDES_INTER, W_INTER, WANN_INTERP%KPOINTS_INTER, SUBGROUP="electron_eigenvalues_kpoints_wan")
#endif

    ALLOCATE(COND_ENERGY(NEDOS,3,3,W%WDES%NCDIJ))
    ALLOCATE(EDOS(NEDOS))

    DO I=1,NEDOS
       EDOS(I)=EFERMI-EMAX_COND+(EMAX_COND*2)/(NEDOS-1)*(I-1)
    ENDDO

    !CALL VH5_WRITE_WAVEDER(IH5OUTFILEID, CHAMS_KPATH, WANN_INTERP%KPOINTS_INTER%VKPT)
    IF (LINTPOL_OPTICS) THEN
       ALLOCATE(EPSDD(NEDOS,3,3),EPSJJ(NEDOS,3,3))
       CALL LRF_EPSILON_COMPUTE(WANN_INTERP%WDES_INTER, W_INTER, LATT_CUR, LATT_INI, INFO, IO, WANN_INTERP%KPOINTS_INTER, SYMM, &
                                EFERMI, DELTAE, EMAX, EMAX_COND, NEDOS, &
                                VEL_KPATH, CHAMS_KPATH, &
                                WPLASMON, CON, BIMAG, WPLASMA_INTER, EPSDD, EPSJJ, COND_ENERGY )

       CALL LR_OPTIC_OUTPUT(EPSDD,EPSJJ,WPLASMA_INTER,CON,BIMAG,NEDOS,DELTAE,KPOINTS,LATT_CUR,INFO,IO)

#ifdef VASP_HDF5
       SUBGROUP = "linear_response_kpoints_wan"
       CALL VH5_WRITE_DIELECTRIC_ENERGIES(IH5OUTFILEID, DELTAE, NEDOS, SUBGROUP)
       CALL VH5_WRITE_DIELECTRIC_DYNAMIC(IH5OUTFILEID, EPSDD, "density_density", SUBGROUP)
       CALL VH5_WRITE_DIELECTRIC_DYNAMIC(IH5OUTFILEID, EPSJJ, "current_current", SUBGROUP)
       CALL VH5_WRITE_CONDUCTIVITY_DYNAMIC(IH5OUTFILEID, NEDOS, EDOS, COND_ENERGY, EFERMI)
#endif
       DEALLOCATE(EPSDD,EPSJJ)
    ELSEIF (LINTPOL_VEL.AND.(LINTPOL_OPTICS)) THEN
       CALL CONDUCTIVITY_ENERGY_RESOLVED( WANN_INTERP%WDES_INTER, W_INTER, VEL_KPATH, EFERMI, &
                                        NEDOS, EMAX_COND, COND_ENERGY, INFO%NELECT, INFO%NUP_DOWN, &
                                        WANN_INTERP%KPOINTS_INTER, SYMM, LATT_CUR, RTIME)
#ifdef VASP_HDF5
       CALL VH5_WRITE_CONDUCTIVITY_DYNAMIC(IH5OUTFILEID, NEDOS, EDOS, COND_ENERGY, EFERMI, SUBGROUP='conductivity_kpoints_wan')
#endif
    ENDIF
    DEALLOCATE(COND_ENERGY,EDOS)

    ! Free memory
    DEALLOCATE(HAM_KPATH,CHAM_KPATH,CHAMS_KPATH,EIG_KPATH,VEL_KPATH)

  END SUBROUTINE WANNIER_INTERPOLATE

  !>@brief Write the eigenvalues and velocities to a file for using with gnuplot
  SUBROUTINE WANNIER_WRITE_GNUPLOT(WANN_INTERP,LATT_CUR,EIG_KPATH,VEL_KPATH,IO)
    USE constant, ONLY: TPI
    TYPE (wannier_interpolator) :: WANN_INTERP
    TYPE (latt),INTENT(IN) :: LATT_CUR
    REAL(q),INTENT(IN) :: EIG_KPATH(:,:,:)
    REAL(q),INTENT(IN) :: VEL_KPATH(:,:,:,:)
    TYPE (in_struct),INTENT(IN) :: IO
    ! local variables
    REAL(q) :: VEL(3),DK(3)
    REAL(q) :: YMIN, YMAX, XSTART, XSTOP, XX
    INTEGER :: IK, ISP, IW

    ! Write results in bands_int.dat file
    IF(IO%IU0>=0) THEN
       OPEN(UNIT=199,FILE='bands_int.dat',FORM='FORMATTED',STATUS='REPLACE')
       DO ISP=1,WANN_INTERP%ISPIN
          DO IW=1,WANN_INTERP%NW
             XX = 0
             XSTART = XX
             DO IK=1,WANN_INTERP%KPOINTS_INTER%NKPTS
                IF (LINTPOL_OPTICS.OR.LINTPOL_VEL) THEN
                   VEL = ABS(VEL_KPATH(IW,IK,ISP,:))
                   WRITE(199,'(5E14.6)') XX,REAL(EIG_KPATH(IW,IK,ISP),q),VEL(1),VEL(2),VEL(3)
                ELSE
                   WRITE(199,'(5E14.6)') XX,REAL(EIG_KPATH(IW,IK,ISP),q)
                ENDIF
                IF (IK<WANN_INTERP%KPOINTS_INTER%NKPTS) THEN
                   DK = WANN_INTERP%KPOINTS_INTER%VKPT(:,IK+1) - WANN_INTERP%KPOINTS_INTER%VKPT(:,IK)
                   CALL DIRKAR(1,DK,LATT_CUR%B)
                   XX = XX + SQRT(SUM((TPI*DK)**2))
                ENDIF
             END DO
             WRITE(199,*)
             XSTOP = XX
          END DO
          WRITE(199,*)
          WRITE(199,*)
       END DO
       CLOSE(199)
    ENDIF

    ! Make a nice gnuplot plot file
    IF ( IO%IU0>=0 ) THEN
       YMIN = MINVAL( EIG_KPATH )
       YMAX = MAXVAL( EIG_KPATH )
       !increase values of YMIN an YMAX to 10% of energy range
       XX=ABS(YMAX-YMIN)*0.1_q
       YMIN=YMIN-XX
       YMAX=YMAX+XX

       OPEN(UNIT=199,FILE='bands_int.gnu',FORM='FORMATTED',STATUS='REPLACE')
       !write header
       WRITE(199,'(A)')'set nokey'
       WRITE(199,'("set xrange [",F12.6,":",F12.6,"]")')XSTART,XSTOP
       WRITE(199,'("set yrange [",F12.6,":",F12.6,"]")')YMIN,YMAX
       WRITE(199,'(A)')' plot "bands_int.dat" w l'
       CLOSE(199)
    ENDIF

  END SUBROUTINE WANNIER_WRITE_GNUPLOT

  !************* SUBROUTINE KPOINTS_SUPERGRID *****************************
  ! 
  ! create k-points grid suitable for Fourier interpolation
  ! fastest index kx, then ky, then kz
  ! kx,ky,kz in (-0.5,0.5] but ordered starting with 0
  ! i.e. [0,0.5] + (-0.5,0)
  !
  !************************************************************************
  SUBROUTINE KPOINTS_SUPERGRID(M1,M2,M3,VKPT)
    INTEGER :: M1,M2,M3
    REAL(q) :: VKPT(3,M1*M2*M3)
    !
    REAL(q) :: TINY=1.E-6_q,D(3)
    INTEGER :: I,I1,I2,I3

    I=0
    D(1)=1.0_q/REAL(M1,q)
    D(2)=1.0_q/REAL(M2,q)
    D(3)=1.0_q/REAL(M3,q)
    DO I3=0,M3-1
       DO I2=0,M2-1
          DO I1=0,M1-1
             I=I+1
             ! map values into interval (-0.5,0.5]
             ! 1. map interval (-0.5,0.5] to (0,1]
             ! 2. then apply MOD function to x
             ! 3. map back to (-0.5,0.5]
             ! MOD( x+ 0.5-eps, 1.0 ) - (0.5-eps)
             VKPT(:,I)=(/ MOD(D(1)*I1+0.5_q-TINY,1.0_q)-(0.5_q-TINY),&
                  MOD(D(2)*I2+0.5_q-TINY,1.0_q)-(0.5_q-TINY),&
                  MOD(D(3)*I3+0.5_q-TINY,1.0_q)-(0.5_q-TINY) &
                  /)
          END DO
       END DO
    END DO
  END SUBROUTINE KPOINTS_SUPERGRID


  !************* SUBROUTINE KPOINT_IN_GRID ********************************
  ! 
  ! Search for 3-vector in array (not limited to kpoints_struct)
  !
  !************************************************************************
  FUNCTION KPOINT_IN_GRID(VKPT,GRID)
    INTEGER :: KPOINT_IN_GRID
    REAL(q) :: VKPT(3), GRID(:,:)
    REAL(q),PARAMETER :: TINY=1E-8_q
    ! local
    INTEGER NK
    DO NK=1,SIZE(GRID,2)
       IF ( ALL( ABS(MOD(VKPT(:)-GRID(:,NK)+10.5_q,1._q)-0.5_q)<TINY ) ) EXIT
    ENDDO

    IF (NK>SIZE(GRID,2)) THEN
       ! no kpoint found, set nk=-1
       NK=-1
    ENDIF

    KPOINT_IN_GRID=NK
  END FUNCTION KPOINT_IN_GRID

  !******************* SUBROUTINE FOURIER_INTERPOL_KTOR ******************************
  !
  !> @brief Transform the from reciprocal (K) to real space (R) using FFT
  !
  !**********************************************************************************
  SUBROUTINE FOURIER_INTERPOL_KTOR(KPOINTS,HAM,HAMR)
    USE profiling
    USE base
    USE mkpoints
    USE lattice
    USE constant
#ifndef __NEC__
    include 'fftw3.f'
#else // __NEC__
    include 'aslfftw3.f'
#endif // __NEC__
    TYPE (kpoints_struct) :: KPOINTS !< Kpoints data structure
    GDEF :: HAM(:,:,:,:) !< Hamiltonian in reciprocal space
    COMPLEX(q),ALLOCATABLE :: HAMR(:,:,:,:,:,:) !< Hamiltonian in real space
    ! local variables
    INTEGER :: N1,N2,N3,NB1,NB2,NKP,NSP,&
         I1,I2,I3,IK,JK,IB1,IB2,ISP
    INTEGER(qi8) :: PLANF
    COMPLEX(q), ALLOCATABLE :: C(:,:,:)
    REAL(q), ALLOCATABLE :: VKPT(:,:)
    INTEGER, ALLOCATABLE :: KINDEX(:)

    PROFILING_START("fourier_interpol_ktor")

    ! coarse grid
    N1=KPOINTS%NKPX
    N2=KPOINTS%NKPY
    N3=KPOINTS%NKPZ

    NB1=SIZE(HAM,1)
    NB2=SIZE(HAM,2)
    NKP=KPOINTS%NKPTS
    NSP=SIZE(HAM,4)
    ALLOCATE(C(0:N1-1,0:N2-1,0:N3-1))
    ALLOCATE(HAMR(NB1,NB2,N1,N2,N3,NSP))

    ! create index array for k-points
    ! such that k-points are ordered by increasing x, then y, then z coordinate
    ! KINDEX(i) is the _original_ position of the i-th ordered k-point
    ! create ordered kpoints grid for Fourier interpolation
    ALLOCATE(KINDEX(NKP),VKPT(3,N1*N2*N3))
    CALL KPOINTS_SUPERGRID(N1,N2,N3,VKPT)
    DO IK=1,NKP
       ! find corresponding kpoints in original grid
       JK=KPOINT_IN_GRID(VKPT(:,IK),KPOINTS%VKPT)
       KINDEX(IK)=JK
    END DO
    DEALLOCATE(VKPT)

    !
    !  K -> R
    !
    ! first call: create a plan
    ! either using FFTW_MEASURE (takes a lot of time)
    ! or using FFTW_ESTIMATE (faster but FFTs less efficient)
    ! Neither of them actually performs a FFT, FFTW_MEASURE overwrites the arrays
    call dfftw_plan_dft_3d(PLANF,N1,N2,N3,C,C,FFTW_FORWARD,FFTW_ESTIMATE)

    DO ISP=1,NSP
       DO IB1=1,NB1
          DO IB2=1,NB2
             ! sort HAM in the same order as k-points
             IK=1
             DO I3=0,N3-1
                DO I2=0,N2-1
                   DO I1=0,N1-1
                      C(I1,I2,I3)=HAM(IB1,IB2,KINDEX(IK),ISP)
                      IK=IK+1
                   END DO
                END DO
             END DO

             call dfftw_execute_dft(PLANF,C,C)

             ! store HAM in real space
             DO I3=0,N3-1
                DO I2=0,N2-1
                   DO I1=0,N1-1
                      HAMR(IB1,IB2,I1+1,I2+1,I3+1,ISP) = C(I1,I2,I3)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    call dfftw_destroy_plan(PLANF)
    DEALLOCATE(C,KINDEX)

    PROFILING_STOP("fourier_interpol_ktor")

  END SUBROUTINE FOURIER_INTERPOL_KTOR

  !> @brief Transform from real space (R) to reciprocal space (K) using slow FFT
  SUBROUTINE FOURIER_INTERPOL_RTOK(HAMR,VKPT_INT,HAM_INT,LAT,DHDK_INT)
    USE profiling
    USE base
    USE mkpoints
    USE lattice
    USE constant
    COMPLEX(q) :: HAMR(:,:,:,:,:,:) !< Hamiltonian in real space
    REAL(q) :: VKPT_INT(:,:) !< K-points onto which to interpolate the Hamiltonian (when KINTER==0)
    GDEF :: HAM_INT(:,:,:,:) !< Interpolated Hamiltonian
    TYPE (latt) :: LAT !< lattice type
    GDEF, OPTIONAL :: DHDK_INT(:,:,:,:,:) !< Calculate derivatives
    ! local variables
    INTEGER :: N1,N2,N3,NSP,&
         I1,I2,I3,IK,JK,NEQ,ISP,&
         NKP_INT,IW1,IW2,IW3
    INTEGER :: NB1,NB2
    COMPLEX(q) :: PHASE
    INTEGER :: DISPL(3),DISPL_ALL(3,10)
    REAL(q) :: DISPL_CART(3),DIST,DIST_FOUND
    REAL(q) :: TINY=1.D-8
    LOGICAL :: DO_DHDK

    DO_DHDK=.FALSE.; IF(PRESENT(DHDK_INT)) DO_DHDK=.TRUE.

    PROFILING_START("fourier_interpol_rtok")

    NKP_INT=SIZE(VKPT_INT,2)

    IF(SIZE(HAM_INT,3)/=NKP_INT)THEN
       CALL vtutor%error("FOURIER_INTERPOL: HAM_INT has not the correct k-point size " // &
          str(SIZE(HAM_INT,3)) // " " // str(NKP_INT))
    END IF
    HAM_INT = 0
    IF(PRESENT(DHDK_INT))THEN
       IF(SIZE(DHDK_INT,3)/=NKP_INT)THEN
          CALL vtutor%error("FOURIER_INTERPOL: DHDK_INT has not the correct k-point size " // &
             str(SIZE(DHDK_INT,3)) // " " // str(NKP_INT))
       END IF
       DHDK_INT=0
    END IF

    NB1 = SIZE(HAMR,1)
    NB2 = SIZE(HAMR,2)
    N1 = SIZE(HAMR,3)
    N2 = SIZE(HAMR,4)
    N3 = SIZE(HAMR,5)
    NSP= SIZE(HAMR,6)

    !
    !  R -> K
    !
    ! transform again to reciprocal space
    DO I3=0,N3-1
       DO I2=0,N2-1
          DO I1=0,N1-1
             !
             ! search the equivalent displacement with the shortest distance
             DIST_FOUND=1E10
             DO IW3=-N3,N3,N3
                DO IW2=-N2,N2,N2
                   DO IW1=-N1,N1,N1
                      ! displacement in direct coordinates
                      ! The actual value of R1 is [(I1-1)+IW1]/N1*(N1*\vec{a1})
                      ! N1 cancels out, so we can use LAT%A directly in DIRKAR
                      ! analogously for R2 and R3
                      DISPL(1)=I1+IW1
                      DISPL(2)=I2+IW2
                      DISPL(3)=I3+IW3
                      ! R-point displacement in cartesian coordinates
                      DISPL_CART=DISPL
                      CALL DIRKAR(1, DISPL_CART(1), LAT%A)
                      DIST=SUM(DISPL_CART**2)
                      ! Here we try to find the minimum distance
                      ! if found we reset the list of equivalent displacements
                      IF (DIST<=DIST_FOUND) THEN
                         DIST_FOUND=DIST
                         NEQ = 0
                      ENDIF
                      ! if the current distance is the same as the minimum distance
                      ! then we start to build the list of equivalent displacements
                      IF (ABS(DIST-DIST_FOUND)<TINY) THEN
                         NEQ=NEQ+1
                         DISPL_ALL(:,NEQ)=DISPL
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             !
             ! find the corresponding points in the fine grid
             ! and distribute the value in C among these equivalent points evenly
             !
             DO ISP=1,NSP
                DO IK=1,NKP_INT
                   PHASE = 0
                   DO JK=1,NEQ
                      PHASE = PHASE+EXP( CITPI*SUM(VKPT_INT(:,IK)*DISPL_ALL(:,JK)) )/REAL(NEQ,q)
                   END DO
                  
                   HAM_INT(:,:,IK,ISP)=HAM_INT(:,:,IK,ISP)+HAMR(:,:,I1+1,I2+1,I3+1,ISP)*PHASE
                   IF (.NOT.DO_DHDK) CYCLE
                   DHDK_INT(:,:,IK,ISP,1) = DHDK_INT(:,:,IK,ISP,1)+CITPI*DISPL_ALL(1,JK)*HAMR(:,:,I1+1,I2+1,I3+1,ISP)*PHASE
                   DHDK_INT(:,:,IK,ISP,2) = DHDK_INT(:,:,IK,ISP,2)+CITPI*DISPL_ALL(2,JK)*HAMR(:,:,I1+1,I2+1,I3+1,ISP)*PHASE
                   DHDK_INT(:,:,IK,ISP,3) = DHDK_INT(:,:,IK,ISP,3)+CITPI*DISPL_ALL(3,JK)*HAMR(:,:,I1+1,I2+1,I3+1,ISP)*PHASE
                ENDDO !IK
             ENDDO !ISP
          ENDDO !I3
       ENDDO !I2
    ENDDO !I3

    HAM_INT=HAM_INT/REAL(N1*N2*N3,q)
    IF (DO_DHDK) DHDK_INT=DHDK_INT/REAL(N1*N2*N3,q)

    PROFILING_STOP("fourier_interpol_rtok")

  END SUBROUTINE FOURIER_INTERPOL_RTOK

  !> @brief Free wannier interpolation object
  SUBROUTINE WANNIER_INTERPOL_FREE(WANN_INTERP)
    TYPE (wannier_interpolator) :: WANN_INTERP
    DEALLOCATE(WANN_INTERP%KPOINTS_COARSE_FULL%VKPT)
    DEALLOCATE(WANN_INTERP%KPOINTS_COARSE_FULL%WTKPT)
    DEALLOCATE(WANN_INTERP%KPOINTS_COARSE_FULL%IDTET)
    DEALLOCATE(WANN_INTERP%KPOINTS_INTER%VKPT)
    DEALLOCATE(WANN_INTERP%KPOINTS_INTER%WTKPT)
    DEALLOCATE(WANN_INTERP%KPOINTS_INTER%IDTET)
    ! The following two WDES are part of the WANN_INTERP structure
    ! but they should nit be deallocated.
    ! They are made of the main WDES + kpoints from the two k-point structures above
    !CALL DEALLOCWDES(WANN_INTERP%WDES_COARSE_FULL)
    !CALL DEALLOCWDES(WANN_INTERP%WDES_INTER)
    IF (ALLOCATED(WANN_INTERP%U)) DEALLOCATE(WANN_INTERP%U)
    IF (ALLOCATED(WANN_INTERP%HAM)) DEALLOCATE(WANN_INTERP%HAM)
    IF (ALLOCATED(WANN_INTERP%CDER)) DEALLOCATE(WANN_INTERP%CDER)
    IF (ALLOCATED(WANN_INTERP%NEQUIV)) DEALLOCATE(WANN_INTERP%NEQUIV)
  END SUBROUTINE WANNIER_INTERPOL_FREE

END MODULE wannier_interpolation
