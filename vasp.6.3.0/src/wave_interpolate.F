#include "symbol.inc"

!************************ MODULE WAVE_INTERPOLATE **********************
!>
!> @brief
!> this module provides functions to interpolate the band structure
!> to an arbitrary k-point grid
!> This routine was writted by Henrique Miranda
!> Based on routine INTERPOLATE_BAND_STR by gK and
!> INTERPOLATE_BANDSTR_GRDSHFT by eFL
!>
!> The function INTERPOLATE_WAVE_INIT doubles the size of the
!> original wavefunctions using REALLOC_WAVE
!> TODO: 1. compute the first derivative with respect to k and
!>       store in wave_interpolator object and use it to interpolate
!>       2. use FFT to interpolate WF to adjacent k-points
!>
!> The function INTERPOLATE_WAVE_KPOINTS can be called interpolating
!> the wavefunctions to a specified batch of k-points
!> only the new wavefunctions are optimized
!>
!> INTERPOLATE_WAVE_FREE restores the underlying wavefunctions
!> to the original state
!
!> Consequences:
!> The INTERPOLATE_WAVE_INIT will reallocate the original wavefunctions
!> After INTERPOLATE_WAVE_FREE the original WFs are restored
!> Should also have some consequences in the FOCK module
!>
!***********************************************************************

MODULE wave_interpolate

  USE prec
  USE wave_high
  USE lattice
  USE poscar
  USE mpimy
  USE mgrid
  USE nonl_high
  USE base
  USE pseudo
  USE kpoints_change
  USE constant
  USE choleski
  USE subrot
  USE david
  USE rmm_diis
  USE hamil_struct_def
  USE tutor, ONLY: vtutor

  IMPLICIT NONE

  !> Mode to use for optional k-point grid
  INTEGER :: KPOINTS_OPT_MODE=1
  !> Logical to activate use of optional k-point grid
  LOGICAL :: LKPOINTS_OPT=.FALSE.

  TYPE wave_interpolator
     !> Number of k-points in the original wavefunction datastructure
     INTEGER :: NKPTS_ORIG
     !> Max number of k-points in every batch
     INTEGER :: MAX_NKPTS_BATCH
     !> Number of k-points in batch
     INTEGER :: NKPTS_BATCH
     !> Coordinates of k-points in batch
     REAL(q),ALLOCATABLE :: KPOINTS_BATCH(:,:)
     !> index of k-points in batch
     INTEGER,ALLOCATABLE :: KPOINTS_INDEX(:)
     !> option to decide with method is used to find intial guess for the wavefunctions
     !>  =1 Use WF in close k-point as initial guess
     !>  =2 Use derivative of the wavefunction to interpolate to new k-points
     INTEGER :: OPTION
     !> Kpoints structure
     TYPE (kpoints_struct) KPOINTS
     !> Map reallocated WF to original NGVECTOR
     INTEGER,ALLOCATABLE :: NGVECTOR(:)
     !> Map reallocated WF to original NINDPW
     INTEGER,ALLOCATABLE :: NINDPW(:,:)
     !> Wavefunction contained for computation of derivatives of the WFs
     TYPE (wavespin) :: WDIFF3(3)
     LOGICAL :: HAS_WAVEDER = .FALSE.
     !> Save status of LFOCKACE
     LOGICAL :: LFOCKACE_SAVE
  END TYPE wave_interpolator

  PUBLIC :: WAVE_INTERPOLATOR_INIT
  PUBLIC :: WAVE_INTERPOLATOR_KPOINTS
  PUBLIC :: WAVE_INTERPOLATOR_FREE

CONTAINS

!***********************************************************************
!> Read input variables controling the behaviour of the wave_interpolator
!***********************************************************************
SUBROUTINE WAVE_INTERPOLATOR_READER(IU0,IU5,IU6)
  USE vaspxml
  USE base
  USE reader_tags
  USE pead, ONLY: PEAD_REQUEST

  INTEGER IU5   ! input device (usually INCAR)
  INTEGER IU0   ! stderr
  INTEGER IU6   ! stdout

  ! local
  INTEGER IERR
  REAL(q) RDUM
  COMPLEX(q) CDUM
  LOGICAL LOPEN,LDUM

  !
  ! interpolation in k using linear response routines
  !
  KPOINTS_OPT_MODE=1
  CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'KPOINTS_OPT_MODE', KPOINTS_OPT_MODE, IERR, WRITEXMLINCAR)
  IF (KPOINTS_OPT_MODE/10==1.OR.KPOINTS_OPT_MODE/10==2) CALL PEAD_REQUEST

  !
  ! additional finer grid
  !
  INQUIRE(FILE="KPOINTS_OPT",EXIST=LKPOINTS_OPT)
  CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LKPOINTS_OPT', LKPOINTS_OPT, IERR, WRITEXMLINCAR)

END SUBROUTINE WAVE_INTERPOLATOR_READER

!***********************************************************************
!> Write information about wave_interpolator to the OUTCAR file
!***********************************************************************
SUBROUTINE WAVE_INTERPOLATOR_WRITE_OUTCAR(IU6)
  INTEGER,INTENT(IN) :: IU6               !> output unit

  IF (IU6>=0) THEN
     WRITE(IU6,*)          ' Optional k-point grid parameters'
     WRITE(IU6,'(A,L4,A)') '   LKPOINTS_OPT  =  ',LKPOINTS_OPT,     '    use optional k-point grid'
     WRITE(IU6,'(A,I6,A)') '   KPOINTS_OPT_MODE=',KPOINTS_OPT_MODE, '    mode for optional k-point grid'
     WRITE(IU6,*)
  ENDIF

END SUBROUTINE WAVE_INTERPOLATOR_WRITE_OUTCAR

!> Write information about the optinal k-points file to the xml file
SUBROUTINE XML_WRITE_WAVE_INTERPOLATOR()
  USE vaspxml

  !local varaibles
  LOGICAL :: LDUM
  INTEGER :: IDUM
  REAL(q) :: RDUM
  COMPLEX(q)  :: CDUM
  CHARACTER(1) :: CHARAC

  CALL XML_TAG("separator","optional k-points parameters")
  CALL XML_INCAR('KPOINTS_OPT_MODE','I',KPOINTS_OPT_MODE,RDUM,CDUM,LDUM,CHARAC,1)
  CALL XML_INCAR('LKPOINTS_OPT','L',IDUM,RDUM,CDUM,LKPOINTS_OPT,CHARAC,1)
  CALL XML_CLOSE_TAG
END SUBROUTINE XML_WRITE_WAVE_INTERPOLATOR

!***********************************************************************
!>
!> @brief Compute the derivative of the wavefunction with respect to k
!> this can be later on used to improve the intial guess of the wavefunction
!>
!***********************************************************************
SUBROUTINE WAVE_INTERPOLATOR_WAVEDER(WAVE_INTERP,W,KPOINTS,GRID,LATT_INI,LATT_CUR,T_INFO,&
                  NONLR_S,NONL_S,LMDIM,GRIDC,GRIDUS,C_TO_US,IRDMAX,P,CQIJ,CDIJ,SV,SYMM,INFO,IO)
  USE lattice, ONLY: latt
  USE subrot_cluster
  USE mlr_optic, ONLY: FOCK_K_DER_ANALYT
  USE mlrf_main, ONLY: LRF_RPHI0
  TYPE (wave_interpolator) WAVE_INTERP
  TYPE (wavespin) :: W
  TYPE (kpoints_struct) :: KPOINTS
  TYPE (grid_3d)        :: GRID
  TYPE (latt)           :: LATT_INI, LATT_CUR
  TYPE (potcar)         :: P(:)
  TYPE (grid_3d)        :: GRIDC
  TYPE (grid_3d)        :: GRIDUS
  TYPE (transit)        :: C_TO_US
  TYPE (info_struct)    :: INFO
  TYPE (type_info)      :: T_INFO
  TYPE (nonlr_struct)   :: NONLR_S
  TYPE (nonl_struct)    :: NONL_S
  TYPE (in_struct)      :: IO
  TYPE (symmetry)       :: SYMM
  INTEGER               :: LMDIM, IRDMAX
  OVERLAP :: CQIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
  OVERLAP :: CDIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
  RGRID   :: SV(DIMREAL(GRID%MPLWV),W%WDES%NCDIJ)
  ! local variables
  TYPE (eigenf_cluster_pointer),POINTER :: DEG_CLUSTER(:,:)
  INTEGER :: IDIR
  LOGICAL :: LDONE

  ! Find degenerate energies
  NULLIFY(DEG_CLUSTER)
  CALL FIND_DEG_CLUSTERS(W%WDES, W, DEG_CLUSTER)

  ! Compute derivative of WFs
  DO IDIR=1,3
     CALL ALLOCW(W%WDES,WAVE_INTERP%WDIFF3(IDIR))
     WAVE_INTERP%WDIFF3(IDIR)%CW   =0
     WAVE_INTERP%WDIFF3(IDIR)%CPROJ=0

     ! We use linear response only (to allow arbitrary grids)
     ! TODO: The derivatives for the original points are always recomputed
     CALL FOCK_K_DER_ANALYT(KPOINTS, GRID, LATT_CUR, LATT_INI, &
          T_INFO,  NONLR_S, NONL_S, W, WAVE_INTERP%WDIFF3(IDIR), LMDIM, P, CQIJ, &
          SYMM, IDIR, LDONE, IO%IU0, IO%IU6)

     CALL LRF_RPHI0( P,NONLR_S,NONL_S,W,LATT_CUR, &
          T_INFO,INFO,IO,GRID,GRIDC,GRIDUS,C_TO_US,IRDMAX, &
          CDIJ,CQIJ,SV,LMDIM,DEG_CLUSTER, IDIR, WAVE_INTERP%WDIFF3(IDIR), LDONE)
  ENDDO
  CALL FREE_DEG_CLUSTERS(W%WDES, DEG_CLUSTER)
  WAVE_INTERP%HAS_WAVEDER = .TRUE.

END SUBROUTINE WAVE_INTERPOLATOR_WAVEDER

!***********************************************************************
!>
!> @brief Initialize the wavefunction interpolator structure
!> increase the size of the wavefunction
!> in the old part the original wavefunctions are stored
!> the new part of the wavefunction has the size of the batches
!> of k-points that will be used
!>
!***********************************************************************

SUBROUTINE WAVE_INTERPOLATOR_INIT(WAVE_INTERP, KPOINTS, GRID, &
                  LATT_CUR, LATT_INI, T_INFO, NONL_S, W, P, &
                  SYMM, IO, OPTION)

  USE fock_glb, ONLY : LFOCKACE
  TYPE (wave_interpolator) WAVE_INTERP
  TYPE (grid_3d)     GRID
  TYPE (kpoints_struct) KPOINTS
  TYPE (latt)        LATT_CUR
  TYPE (latt)        LATT_INI
  TYPE (type_info)   T_INFO
  TYPE (nonl_struct) NONL_S
  TYPE (wavespin), TARGET :: W
  TYPE (symmetry) SYMM
  TYPE (potcar)      P(:)
  TYPE (wavespin) :: WDIFF3(3)
  TYPE (in_struct) :: IO
  INTEGER OPTION

!    local
  INTEGER IK
  INTEGER NKPTS_ORIG, MAX_NKPTS_BATCH
  REAL(q),ALLOCATABLE :: FAKE_VKPT(:,:)
  REAL(q),ALLOCATABLE,TARGET :: TMP_VKPT(:,:), TMP_WTKPT(:)

  CALL CHECK_FULL_KPOINTS

  ! original number of k-points
  NKPTS_ORIG      = W%WDES%NKPTS
  ! can change this to any number to configure the size of the batches
  MAX_NKPTS_BATCH = NKPTS_ORIG

  IF (IO%IU0>=0) WRITE(IO%IU0,*) "Start KPOINTS_OPT (optional k-point list driver)"

  ! Save information to map the G vectors from the IBZ
  ALLOCATE(WAVE_INTERP%NGVECTOR(SIZE(W%WDES%NGVECTOR,1)))
  ALLOCATE(WAVE_INTERP%NINDPW(SIZE(W%WDES%NINDPW,1),SIZE(W%WDES%NINDPW,2)))
  WAVE_INTERP%NGVECTOR = W%WDES%NGVECTOR
  WAVE_INTERP%NINDPW = W%WDES%NINDPW

  ! Resize k-point structure
  CALL COPY_KPOINTS(KPOINTS,WAVE_INTERP%KPOINTS)
  WAVE_INTERP%KPOINTS%NKDIM = NKPTS_ORIG+MAX_NKPTS_BATCH
  WAVE_INTERP%KPOINTS%NKPTS = NKPTS_ORIG+MAX_NKPTS_BATCH
  ALLOCATE(TMP_VKPT(3,WAVE_INTERP%KPOINTS%NKPTS))
  ALLOCATE(TMP_WTKPT(WAVE_INTERP%KPOINTS%NKPTS))
  TMP_VKPT  = 0_q; TMP_VKPT(:,:NKPTS_ORIG) = WAVE_INTERP%KPOINTS%VKPT
  TMP_WTKPT = 0_q; TMP_WTKPT(:NKPTS_ORIG)  = WAVE_INTERP%KPOINTS%WTKPT
  DEALLOCATE(WAVE_INTERP%KPOINTS%WTKPT)
  DEALLOCATE(WAVE_INTERP%KPOINTS%VKPT)
  WAVE_INTERP%KPOINTS%VKPT =>TMP_VKPT
  WAVE_INTERP%KPOINTS%WTKPT=>TMP_WTKPT
  CALL SETUP_FULL_KPOINTS(WAVE_INTERP%KPOINTS,LATT_CUR, T_INFO%NIONS, SYMM%ROTMAP, SYMM%MAGROT, SYMM%ISYM, &
       SYMM%ISYM>=0.AND..NOT.W%WDES%LNONCOLLINEAR, -1,IO%IU0)

  ! save important dimensions in the data structure 
  WAVE_INTERP%NKPTS_ORIG      = NKPTS_ORIG
  WAVE_INTERP%MAX_NKPTS_BATCH = MAX_NKPTS_BATCH
  IF (OPTION/10==2.AND..NOT.WAVE_INTERP%HAS_WAVEDER) THEN
     CALL vtutor%bug('To use WAVEDER to obtain better estimatives one must first call WAVE_INTERPOLATOR_WAVEDER',__FILE__,__LINE__)
  ENDIF
  WAVE_INTERP%OPTION          = OPTION

  ALLOCATE(WAVE_INTERP%KPOINTS_BATCH(3,WAVE_INTERP%MAX_NKPTS_BATCH))
  ALLOCATE(WAVE_INTERP%KPOINTS_INDEX(WAVE_INTERP%MAX_NKPTS_BATCH))
  WAVE_INTERP%NKPTS_BATCH = 0

  CALL KPAR_SYNC_ALL(W%WDES,W)
  CALL RE_GEN_LAYOUT( GRID, W%WDES, WAVE_INTERP%KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
  CALL REALLOCATE_WAVE( W, GRID, W%WDES, NONL_S, T_INFO, P, LATT_CUR, NKSTOP=WAVE_INTERP%NKPTS_ORIG)
  IF (WAVE_INTERP%OPTION/10==2) THEN
     CALL REALLOCATE_WAVE( WAVE_INTERP%WDIFF3(1), GRID, W%WDES, NONL_S, T_INFO, P, LATT_CUR, NKSTOP=WAVE_INTERP%NKPTS_ORIG)
     CALL REALLOCATE_WAVE( WAVE_INTERP%WDIFF3(2), GRID, W%WDES, NONL_S, T_INFO, P, LATT_CUR, NKSTOP=WAVE_INTERP%NKPTS_ORIG)
     CALL REALLOCATE_WAVE( WAVE_INTERP%WDIFF3(3), GRID, W%WDES, NONL_S, T_INFO, P, LATT_CUR, NKSTOP=WAVE_INTERP%NKPTS_ORIG)
  ENDIF
  ! Save status of LFOCKACE
  WAVE_INTERP%LFOCKACE_SAVE = LFOCKACE
  LFOCKACE = .FALSE.

END SUBROUTINE WAVE_INTERPOLATOR_INIT


!***********************************************************************
!>
!> @brief Interpolate the wavefunctions to a new set of k-points
!> The k-point are passed in the WAVE_INTERP datastructure:
!> WAVE_INTERP%NKPTS_BATCH        - the number of k-points in the batch
!> WAVE_INTERP%KPOINTS_BATCH(3,:) - the coordinates of the k-points
!> and must be set before calling this routine
!>
!***********************************************************************

SUBROUTINE WAVE_INTERPOLATOR_KPOINTS(WAVE_INTERP, HAMILTONIAN, GRID, &
                  LATT_CUR, LATT_INI, T_INFO, NONLR_S, NONL_S, W, LMDIM, P, SV, &
                  CQIJ, CDIJ, SYMM, INFO, IU6, IU0)

  USE density_of_states, ONLY: DENINI
  USE pead, ONLY : WDES_FULL_PEAD, PEAD_WA_ROTATE, PEAD_DPSI_DK_1K
  USE fock_ace, ONLY : LFOCK_ACE,FOCK_ACE_ALLOCW,WACE,FOCK_ACE_CONSTRUCT,FOCK_ACE_DEALLOCW
#ifdef fock_dblbuf
  USE fock_dbl, ONLY : FOCK_ALL_DBLBUF
#endif

  ! arguments 
  TYPE (wave_interpolator) WAVE_INTERP
  TYPE (ham_handle)  HAMILTONIAN
  TYPE (grid_3d)     GRID
  TYPE (latt)        LATT_CUR
  TYPE (latt)        LATT_INI
  TYPE (type_info)   T_INFO
  TYPE (nonlr_struct)NONLR_S
  TYPE (nonl_struct) NONL_S
  TYPE (potcar)      P(:)
  TYPE (wavespin), TARGET :: W
  TYPE (info_struct)    :: INFO
  INTEGER IU6, IU0
  INTEGER LMDIM
  OVERLAP         CQIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
  OVERLAP         CDIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
  TYPE (symmetry) SYMM
  RGRID           SV(DIMREAL(GRID%MPLWV),W%WDES%NCDIJ)

  ! local
  TYPE (energy)   :: E
  TYPE (wavespin) :: W_NEW
  TYPE (wavespin) :: DW_NEW(3)
  TYPE (wavefun)  :: WUP
  TYPE (wavefun)  :: WDW

  TYPE(wavefuna)  :: WK_FBZ, WK_IBZ, WK, DWK_FBZ(3)
  TYPE(wavedes1)  :: WDESK_FBZ, WDESK_IBZ, WDESK
  TYPE(wavedes)   :: WDES_FBZ
  COMPLEX(q),ALLOCATABLE :: CR(:)

  INTEGER IDIR, IDIR_FBZ
  INTEGER I,J,K,L,IDX
  INTEGER ISP, N, M, IK, IKFULL, IKIBZ, CLOSEST_IKFULL
  INTEGER NKPTS_ORIG, IELM, NSIM, NB, ICOUEV
  REAL(q) DIST_CART, MIN_DIST, EXHF, DESUM1, RMS, OLD_TOTEN, EININ
  INTEGER,ALLOCATABLE :: K_MAP_FROM_IBZ(:,:)
  REAL(q) VKPT1(3),VKPT2(3),DVKPT(3),DVKPT_TMP(3)
  REAL(q) S(3,3)
  CHARACTER(LEN=5) :: STR

  COMPLEX(q),ALLOCATABLE :: CWL(:), CWR(:)
  COMPLEX(q),ALLOCATABLE :: CW_SAVE(:,:)
  GDEF,ALLOCATABLE :: CPROJ_SAVE(:,:)

  NKPTS_ORIG = WAVE_INTERP%NKPTS_ORIG

  ! Consistency check
  IF (WAVE_INTERP%NKPTS_BATCH > WAVE_INTERP%MAX_NKPTS_BATCH) CALL VTUTOR%ERROR("Number of k-points in batch is larger than MAX_NKPTS_BATCH")
  IF (WAVE_INTERP%NKPTS_BATCH < 1)                           CALL VTUTOR%ERROR("Number of k-points in batch is smaller than 1")

  SELECT CASE (WAVE_INTERP%OPTION/10)
     ! Use random wavefunctions
     CASE(0)
        ! Replace the k-points with k-points that we want to interpolate to
        DO IK=1,WAVE_INTERP%NKPTS_BATCH
           W%WDES%VKPT(:,NKPTS_ORIG+IK)=WAVE_INTERP%KPOINTS_BATCH(:,IK)
        ENDDO
        CALL KPAR_SYNC_ALL(W%WDES,W)
        CALL RE_GEN_LAYOUT( GRID, W%WDES, WAVE_INTERP%KPOINTS, LATT_CUR, LATT_INI, IU6, IU0)
#if 1
        ! We need to map the orbitals as changing the k-points might change the number of planewaves
        ! and a mapping of the old to new g-vectors is needed even for the original k-points
        CALL REALLOCATE_WAVE( W, GRID, W%WDES, NONL_S, T_INFO, P, LATT_CUR, NKSTOP=NKPTS_ORIG)
#else
        CALL ALLOCW(W%WDES,W_NEW,WUP,WDW)

        ! copy wavefunction coeficients and projectors
        DO ISP=1,W%WDES%ISPIN
           DO IK=1,NKPTS_ORIG
              DO N=1,W%WDES%NBANDS
                 W_NEW%CW(:,N,IK,ISP)   = W%CW(:,N,IK,ISP)
                 W_NEW%CPROJ(:,N,IK,ISP)= W%CPROJ(:,N,IK,ISP)
                 W_NEW%CELTOT(N,IK,ISP) = W%CELTOT(N,IK,ISP)
                 W_NEW%FERTOT(N,IK,ISP) = W%FERTOT(N,IK,ISP)
                 W_NEW%AUXTOT(N,IK,ISP) = W%AUXTOT(N,IK,ISP)
              ENDDO
           ENDDO
        ENDDO

        CALL DEALLOCW(W)
        W=W_NEW

        IF (NONL_S%LRECIP) THEN
           CALL NONL_DEALLOC(NONL_S)
           CALL NONL_ALLOC (NONL_S,T_INFO,P,W%WDES, .FALSE.)

           CALL SPHER(GRID, NONL_S, P, W%WDES, LATT_CUR, 1)
           CALL PHASE(W%WDES,NONL_S,0)
        ENDIF
#endif

        DO ISP=1,W%WDES%ISPIN
            CALL DENINI(W%FERTOT(:,:,ISP),W%WDES%NB_TOT,W%WDES%NKPTS,INFO%NELECT,W%WDES%LNONCOLLINEAR, NKSTART=NKPTS_ORIG+1)
        ENDDO
        ! initialize with random numbers
        EININ = INFO%ENINI !INFO%ENMAX
        CALL WFINIT(W%WDES, W, EININ, NKSTART=NKPTS_ORIG+1)

     ! Use a closest WF as intial guess
     CASE(1,2)
        ! To initialize the new wavefunctions with a good initial guess
        ! we find the closest k-point from the coarse grid in the FBZ
        ! Loop over all the points to be interpolated
        ALLOCATE(K_MAP_FROM_IBZ(2,WAVE_INTERP%NKPTS_BATCH))
        DO IK=1,WAVE_INTERP%NKPTS_BATCH
           VKPT1 = WAVE_INTERP%KPOINTS_BATCH(:,IK)
           ! Find closest point in full grid using naive algorithm
           ! TODO: fast nearest-neighbour search
           MIN_DIST = 1e8
           CLOSEST_IKFULL = -1
           DO IKFULL=1,KPOINTS_FULL_ORIG%NKPTS
              VKPT2 = KPOINTS_FULL_ORIG%VKPT(:,IKFULL)
              DVKPT = VKPT1-VKPT2
              CALL DIRKAR(1, DVKPT(1), LATT_CUR%B)
              DIST_CART = DVKPT(1)**2+DVKPT(2)**2+DVKPT(3)**2
              IF (DIST_CART>MIN_DIST) CYCLE
              MIN_DIST = DIST_CART
              CLOSEST_IKFULL = IKFULL
           END DO
           K_MAP_FROM_IBZ(1,IK) = CLOSEST_IKFULL
           K_MAP_FROM_IBZ(2,IK) = KPOINTS_FULL_ORIG%NEQUIV(CLOSEST_IKFULL)
        END DO

        ! Replace the k-points with k-points that we want to interpolate to
        DO IK=1,WAVE_INTERP%NKPTS_BATCH
           W%WDES%VKPT(:,NKPTS_ORIG+IK)=WAVE_INTERP%KPOINTS_BATCH(:,IK)
        ENDDO
        CALL KPAR_SYNC_ALL(W%WDES,W)
        CALL RE_GEN_LAYOUT(GRID, W%WDES, WAVE_INTERP%KPOINTS, LATT_CUR, LATT_INI, IU6, IU0)
        CALL REALLOCATE_WAVE( W, GRID, W%WDES, NONL_S, T_INFO, P, LATT_CUR, NKSTOP=NKPTS_ORIG)
        IF (WAVE_INTERP%OPTION/10==2) THEN
           CALL REALLOCATE_WAVE( WAVE_INTERP%WDIFF3(1), GRID, W%WDES, NONL_S, T_INFO, P, LATT_CUR, NKSTOP=NKPTS_ORIG)
           CALL REALLOCATE_WAVE( WAVE_INTERP%WDIFF3(2), GRID, W%WDES, NONL_S, T_INFO, P, LATT_CUR, NKSTOP=NKPTS_ORIG)
           CALL REALLOCATE_WAVE( WAVE_INTERP%WDIFF3(3), GRID, W%WDES, NONL_S, T_INFO, P, LATT_CUR, NKSTOP=NKPTS_ORIG)
        ENDIF

        ! to map from the adjacent k-point in the FBZ to the current point
        ! we use Fourier tranformation in real space
        ALLOCATE(CR(W%WDES%GRID%RL%NP))
        W_NEW=W
        W_NEW%WDES=>WDES_FULL_PEAD
        IF (WAVE_INTERP%OPTION/10==2) THEN
           ALLOCATE(CWR(GRID%NPLWV))
           ALLOCATE(CWL(W%WDES%NRPLWV))
           ALLOCATE(CW_SAVE(SIZE(W_NEW%CW,1),SIZE(W_NEW%CW,2)))
           ALLOCATE(CPROJ_SAVE(SIZE(W_NEW%CPROJ,1),SIZE(W_NEW%CPROJ,2)))
           DO IDIR=1,3
              DW_NEW(IDIR)=WAVE_INTERP%WDIFF3(IDIR)
              DW_NEW(IDIR)%WDES=>WDES_FULL_PEAD
           ENDDO
        ENDIF
        DO ISP=1,W%WDES%ISPIN
           DO IK=1,WAVE_INTERP%NKPTS_BATCH
              IKFULL = K_MAP_FROM_IBZ(1,IK)
              IKIBZ  = K_MAP_FROM_IBZ(2,IK)
              ! Create WDESK
              CALL SETWDES(W%WDES,WDESK,NKPTS_ORIG+IK)
              ! Create WDESK_FBZ
              CALL SETWDES(W_NEW%WDES,WDESK_FBZ,IKFULL)
              WDESK_FBZ%NRPLWV = W%WDES%NRPLWV
              IF (IKFULL==IKIBZ) THEN
                 ! k is a kpoint in the IBZ
                 WK_FBZ=ELEMENTS(W_NEW,WDESK_FBZ,ISP)
              ELSE
                 ! k is not a kpoint in the IBZ
                 CALL NEWWAVA(WK_FBZ,WDESK_FBZ,WDESK_FBZ%NBANDS)
                 CALL PEAD_WA_ROTATE(W_NEW,P,LATT_CUR,ISP,WK_FBZ)
              ENDIF
              IF (WAVE_INTERP%OPTION/10==2) THEN
                 ! Rotate all directions from IBZ to FBZ
                 DO IDIR=1,3
                    IF (IKFULL==IKIBZ) THEN
                       ! k is a kpoint in the IBZ
                       DWK_FBZ(IDIR)=ELEMENTS(DW_NEW(IDIR),WDESK_FBZ,ISP)
                    ELSE
                       ! k is not a kpoint in the IBZ
                       CALL NEWWAVA(DWK_FBZ(IDIR),WDESK_FBZ,WDESK_FBZ%NBANDS)
                       CALL PEAD_WA_ROTATE(DW_NEW(IDIR),P,LATT_CUR,ISP,DWK_FBZ(IDIR))
                    ENDIF
                 ENDDO
                 ! determine transformation matrix in real space when going from IBZ to FBZ k-point
                 S=0
                 DO L=1,3
                    DO K=1,3
                       DO J=1,3
                          DO I=1,3
                             S(L,I)=S(L,I)+LATT_CUR%A(L,K)*&
                                    KPOINTS_FULL%ISYMOP(J,K,IKFULL)*LATT_CUR%B(I,J)
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
                 ! k-point difference
                 DVKPT = WDESK%VKPT - WDESK_FBZ%VKPT
                 CALL KARDIR(1,DVKPT,LATT_CUR%B)
                 ! transform the k-point difference
                 DVKPT_TMP = DVKPT
                 DVKPT = MATMUL(S,DVKPT_TMP)
              ENDIF
              ! for all the bands
              DO NB=1,W%WDES%NBANDS
                 IF (WAVE_INTERP%OPTION/10==2) THEN
                    !CWL = WK_FBZ%CW(:,NB)
                    CWL = WK_FBZ%CW(:,NB) + CITPI*( &
                                            DVKPT(1)*DWK_FBZ(1)%CW(:,NB) + &
                                            DVKPT(2)*DWK_FBZ(2)%CW(:,NB) + &
                                            DVKPT(3)*DWK_FBZ(3)%CW(:,NB) )
                    ! rec (fbz point) -> real
                    !CALL FFTWAV(WDESK_FBZ,CR,WK_FBZ%CW(:,NB))
                    CALL FFTWAV(WDESK_FBZ%NGVECTOR,WDESK_FBZ%NINDPW,CR,CWL,WDESK_FBZ%GRID)
                    ! real -> rec (point that we want)
                    CALL FFTEXT(WDESK%NGVECTOR,WDESK%NINDPW,CR,W%CW(:,NB,NKPTS_ORIG+IK,ISP),WDESK%GRID,.FALSE.)
                    ! For the projectors we just copy from the FBZ
                    !W%CPROJ(:,NB,NKPTS_ORIG+IK,ISP) = WK_FBZ%CPROJ(:,NB)
                    W%CPROJ(:,NB,NKPTS_ORIG+IK,ISP) = WK_FBZ%CPROJ(:,NB) + CITPI*( &
                                                      DVKPT(1)*DWK_FBZ(1)%CPROJ(:,NB) + &
                                                      DVKPT(2)*DWK_FBZ(2)%CPROJ(:,NB) + &
                                                      DVKPT(3)*DWK_FBZ(3)%CPROJ(:,NB) )
                    W%CELTOT(NB,NKPTS_ORIG+IK,ISP) = W%CELTOT(NB,IKIBZ,ISP) + CITPI*( &
                                                     DVKPT(1)*DW_NEW(1)%CELTOT(NB,IKIBZ,ISP) + &
                                                     DVKPT(2)*DW_NEW(2)%CELTOT(NB,IKIBZ,ISP) + &
                                                     DVKPT(3)*DW_NEW(3)%CELTOT(NB,IKIBZ,ISP) )
                 ELSE
                    ! rec (fbz point) -> real
                    CALL FFTWAV(WDESK_FBZ%NGVECTOR,WDESK_FBZ%NINDPW,CR,WK_FBZ%CW(:,NB),WDESK%GRID)
                    ! real -> rec (point that we want)
                    CALL FFTEXT(WDESK%NGVECTOR,WDESK%NINDPW,CR,W%CW(:,NB,NKPTS_ORIG+IK,ISP),WDESK%GRID,.FALSE.)
                    ! For the projectors we just copy from the FBZ
                    W%CPROJ(:,NB,NKPTS_ORIG+IK,ISP)= WK_FBZ%CPROJ(:,NB)
                    W%CELTOT(NB,NKPTS_ORIG+IK,ISP) = W%CELTOT(NB,IKIBZ,ISP)
                 ENDIF
                 ! Copy the rest of the information
                 !W%CELTOT(NB,NKPTS_ORIG+IK,ISP) = W%CELTOT(NB,IKIBZ,ISP)
                 W%FERTOT(NB,NKPTS_ORIG+IK,ISP) = W%FERTOT(NB,IKIBZ,ISP)
                 W%AUXTOT(NB,NKPTS_ORIG+IK,ISP) = W%AUXTOT(NB,IKIBZ,ISP)
              ENDDO
              IF (IK/=KPOINTS_FULL%NEQUIV(IK)) THEN
                 CALL DELWAVA(WK_FBZ)
                 IF (WAVE_INTERP%OPTION/10==2) THEN
                    CALL DELWAVA(DWK_FBZ(1))
                    CALL DELWAVA(DWK_FBZ(2))
                    CALL DELWAVA(DWK_FBZ(3))
                 ENDIF
              ENDIF
           ENDDO ! kpoints batch
        ENDDO ! spin index
        DEALLOCATE(CR)
        DEALLOCATE(K_MAP_FROM_IBZ)

        DO ISP=1,W%WDES%ISPIN
            CALL DENINI(W%FERTOT(:,:,ISP),W%WDES%NB_TOT,W%WDES%NKPTS,INFO%NELECT,W%WDES%LNONCOLLINEAR, NKSTART=NKPTS_ORIG+1)
        ENDDO

        ! initialize with random numbers beyond batch
        EININ = INFO%ENINI !INFO%ENMAX
        CALL WFINIT(W%WDES, W, EININ, NKSTART=NKPTS_ORIG+WAVE_INTERP%NKPTS_BATCH+1)
  END SELECT

  ! Prepare call to diag by orthonormalizing the orbitals
  CALL SET_DATAKE(W%WDES, LATT_CUR%B)
  IF (NONLR_S%LREAL) THEN
     CALL RSPHER(GRID,NONLR_S,LATT_CUR)
  ELSE
     CALL SPHER(GRID,NONL_S,P,W%WDES,LATT_CUR, 1)
  ENDIF
  CALL PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,W)
  CALL ORTHCH(W%WDES,W, W%WDES%LOVERL, LMDIM, CQIJ)

  IF (LFOCK_ACE()) THEN
    CALL FOCK_ACE_ALLOCW(W%WDES)
  ENDIF

  ! Call diagonalization
  NSIM=W%WDES%NSIM*2
#ifdef MPI
  NSIM=((W%WDES%NSIM*2+W%WDES%COMM_INTER%NCPU-1)/ &
         W%WDES%COMM_INTER%NCPU)*W%WDES%COMM_INTER%NCPU
#endif

  OLD_TOTEN =BANDSTRUCTURE_ENERGY_W1(W%WDES, W, NKPTS_ORIG+1)
  IF (IU0>=0) THEN
    WRITE(17,*)  '      N       E                     dE             ncg'
    WRITE(IU0,*) '      N       E                     dE             ncg'
  END IF

  DO IELM=1,INFO%NELM

#ifdef fock_dblbuf
     IF (LFOCK_ACE()) THEN
       CALL FOCK_ALL_DBLBUF(W%WDES,W,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ, &
                            NKSTART=NKPTS_ORIG+1,EX=E%EXHF,EX_ACFDT=E%EXHF_ACFDT,XI=WACE)
  !    CALL FOCK_ACE_ACC_PAW(W%WDES,W,LATT_CUR,NONL_S,NONLR_S)
       CALL FOCK_ACE_CONSTRUCT(W%WDES,W,NKSTART=NKPTS_ORIG+1)
     ENDIF
#endif

     SELECT CASE (MOD(WAVE_INTERP%OPTION,10))
        CASE (0)
          ! Lazy mode: do nothing!
        CASE (1)
             CALL EDDAV(HAMILTONIAN,P, GRID,INFO,LATT_CUR,NONLR_S,NONL_S, &
                       W,W%WDES, NSIM, LMDIM,CDIJ,CQIJ, RMS, DESUM1,ICOUEV, SV, &
                       E%EXHF, IU6,IU0, LDELAY=.FALSE., LSUBROTI=INFO%LDIAG, LHF=.TRUE.,&
                       LEMPTY=.FALSE.,NKSTART=NKPTS_ORIG+1)
             STR='DAV: '
     END SELECT

     E%EBANDSTR=BANDSTRUCTURE_ENERGY_W1(W%WDES, W, NKPTS_ORIG+1)

     IF (IU0>=0) THEN
        WRITE(17,  '(A,I3,A,E20.12,A,E12.5,A,I6)') STR,IELM,'   ',E%EBANDSTR,'   ',E%EBANDSTR-OLD_TOTEN,'   ',ICOUEV
        WRITE(IU0, '(A,I3,A,E20.12,A,E12.5,A,I6)') STR,IELM,'   ',E%EBANDSTR,'   ',E%EBANDSTR-OLD_TOTEN,'   ',ICOUEV
     END IF

     IF (ABS(OLD_TOTEN-E%EBANDSTR) < ABS(INFO%EDIFF) ) EXIT
     OLD_TOTEN = E%EBANDSTR

  ENDDO

  IF (LFOCK_ACE()) CALL FOCK_ACE_DEALLOCW

  CONTAINS
  ! The difference in this routine is that we don't use the weights of the k-points
  ! to compute the change of the band-structure. This is necessary as the correct
  ! weights of the original points should be mantained for Hybrid calculations
  ! TODO: Check whether to remove also FETOT as we might want to also converge empty states
  FUNCTION BANDSTRUCTURE_ENERGY_W1(WDES, W, NKSTART) RESULT(E)
    IMPLICIT NONE
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W      ! wavefunction
    REAL(q)            E
    INTEGER            NKSTART
    INTEGER            ISP, NK, NB

    E=0
    DO ISP=1,WDES%ISPIN
       DO NK=NKSTART,WDES%NKPTS
#ifdef MPI
          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
          DO NB=1,WDES%NB_TOT
             E=E+WDES%RSPIN* REAL( W%CELTOT(NB,NK,ISP) ,KIND=q) *W%FERTOT(NB,NK,ISP)
          ENDDO
       ENDDO
    ENDDO
    CALLMPI( M_sum_d(WDES%COMM_KINTER, E, 1))
  END FUNCTION BANDSTRUCTURE_ENERGY_W1

END SUBROUTINE WAVE_INTERPOLATOR_KPOINTS


!***********************************************************************
!>
!> @brief Restore wavefunctions to previous state and free interpolator structure
!>
!***********************************************************************

SUBROUTINE WAVE_INTERPOLATOR_FREE(WAVE_INTERP, KPOINTS, GRID, &
                  LATT_CUR, LATT_INI, T_INFO, NONL_S, W, P, SYMM) 

  USE fock_glb, ONLY : LFOCKACE
  TYPE (wave_interpolator) WAVE_INTERP
  TYPE (kpoints_struct) KPOINTS
  TYPE (grid_3d)     GRID
  TYPE (latt)        LATT_CUR
  TYPE (latt)        LATT_INI
  TYPE (type_info)   T_INFO
  TYPE (nonl_struct) NONL_S
  TYPE (potcar)      P(:)
  TYPE (wavespin), TARGET :: W
  TYPE (symmetry) SYMM
  ! local varaibles
  INTEGER :: IDIR

  !TODO: no need to re-read the k-points
#ifdef oldsym
  CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
       SYMM%ISYM>=0.AND..NOT.W%WDES%LSORBIT.AND..NOT.W%WDES%LSPIRAL, &
       T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,-1,-1)
#else
  CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
       SYMM%ISYM>=0.AND..NOT.W%WDES%LNONCOLLINEAR, &
       T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,-1,-1)
#endif
  CALL RE_GEN_LAYOUT( GRID, W%WDES, KPOINTS, LATT_CUR, LATT_INI, -1, -1)
  CALL REALLOCATE_WAVE( W, GRID, W%WDES, NONL_S, T_INFO, P, LATT_CUR)
  DEALLOCATE(WAVE_INTERP%KPOINTS_BATCH)
  DEALLOCATE(WAVE_INTERP%KPOINTS_INDEX)
  IF (WAVE_INTERP%OPTION/10==2) THEN
     DO IDIR=1,3
         CALL DEALLOCW(WAVE_INTERP%WDIFF3(IDIR))
     ENDDO
  ENDIF
  ! restore status of LFOCKACE
  LFOCKACE = WAVE_INTERP%LFOCKACE_SAVE

END SUBROUTINE WAVE_INTERPOLATOR_FREE

END MODULE wave_interpolate
