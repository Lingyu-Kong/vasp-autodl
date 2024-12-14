#include "symbol.inc"
MODULE mgrid_struct_def
  USE prec
  USE mpimy

  TYPE layout
     INTEGER NP                   ! local number of grid points
     INTEGER NALLOC               ! allocation required
     INTEGER NFAST                ! which index is fast (1-x, 2-y, 3-z)
     INTEGER NCOL                 ! number of columns in the grid
     INTEGER NROW                 ! number of elements in each column
     INTEGER,POINTER :: I2(:)     ! y/z/x-index of each column
     INTEGER,POINTER :: I3(:)     ! z/x/y-index of each column
     INTEGER,POINTER :: INDEX(:,:)! column index for each yz, zx or xy pair
  END TYPE layout

  TYPE grid_map
     INTEGER,POINTER :: N(:)      ! number of elements send by each node
     INTEGER,POINTER :: PTR(:)    ! sum_j=1,I N(j)
     INTEGER,POINTER :: RMT(:)    ! remote address for send  (shmem t3d)
     INTEGER,POINTER :: TBL(:)    ! address of each element
     ! inverse transformation (i.e. receiver information)
     INTEGER,POINTER :: NI(:)
     INTEGER,POINTER :: PTRI(:)
     INTEGER,POINTER :: TBLI(:)
     INTEGER,POINTER :: RMTI(:)
     LOGICAL LOCAL                ! all information is local
     LOGICAL LOCAL_COPY           ! no data redistribution required
  END TYPE grid_map

  TYPE grid_3d
     !only  GRID
     INTEGER NGX,NGY,NGZ          ! number of grid points in x,y,z
     INTEGER NGX_rd,NGY_rd,NGZ_rd ! in the complex mode the _rd values
     ! are equal to NGX, NGY, NGZ
     ! if a real to complex FFT is used only half of the data
     ! are stored in one direction and the corresponding NG?_rd is 
     ! set to (NG?+1)/2
     INTEGER NPLWV                ! total number of grid points NGX*NGY*NGZ
     INTEGER MPLWV                ! allocation in complex words required to do in place FFT's
     INTEGER NGPTAR(3)            ! equivalent to /(NGX,NGY,NGZ/)
     INTEGER,POINTER :: LPCTX(:)  ! loop counters in x
     INTEGER,POINTER :: LPCTY(:)  ! loop counters in y
     INTEGER,POINTER :: LPCTZ(:)  ! loop counters in z
     ! loop counters, in which the unbalanced contribution is zeroed
     INTEGER,POINTER :: LPCTX_(:) ! loop counters in x
     INTEGER,POINTER :: LPCTY_(:) ! loop counters in y
     INTEGER,POINTER :: LPCTZ_(:) ! loop counters in z
     ! reciprocal space layout (x is always  fast index)
     TYPE(layout)    :: RC
     ! intermediate layout (y is always  fast index, used only in parallel version)
     TYPE(layout)    :: IN
     ! real space layout   (x or z is the fast index)
     TYPE(layout)    :: RL
     ! real space layout for FFT  (x or z is the fast index)
     ! this structure is usually equivalent to RL (and hence points to RL)
     ! only if the serial version is used for the FFT of wavefunctions,
     ! the structure differs from RL for the GRID_SOFT structure
     ! in this case, the FFT has z as fast index as required for the parallel FFT
     ! but the RL structure has x as fast index to be compatible to the FFT
     ! of wavefunctions
     TYPE(layout), POINTER  :: RL_FFT
     ! information only required for real space representation
     INTEGER NGZ_complex          ! number of grid points for z fast
     ! mapping for parallel version
     TYPE(grid_map)  :: RC_IN     ! recip -> intermediate
     TYPE(grid_map)  :: IN_RL     ! intermediate -> real space
     TYPE(communic), POINTER :: COMM,COMM_KIN,COMM_KINTER
     LOGICAL         :: LREAL     ! are data stored as complex or real numbers in real space
     LOGICAL         :: REAL2CPLX ! real to complex FFT
     REAL(q), POINTER:: FFTSCA(:,:) ! scaling factors for real to complex fft (wavefunction fft) in reciprocal space (only defined for PW in cutoff-sphere)
     REAL(q), POINTER:: FFTWEIGHT(:)! for real space version, weight of each coefficient in reciprocal space for all coefficients
     ! plane waves with N1=0 (NX=0) within the cutoff sphere correspond to two 
     ! G vectors on the 3D FFT grid
     ! IND_IN_SPHERE stores the PW indices (within the sphere) with N1=0
     ! NINDPWCONJG stores the position where the conjugated coefficient needs to be stored
     INTEGER, POINTER:: NINDPWCONJG(:)
     INTEGER, POINTER:: IND_IN_SPHERE(:)
  END TYPE grid_3d
  ! comments:
  ! REAL2CPLX determines wether a real to complex FFT is use
  ! LREAL     determines wether the data in real space are 
  !           stored in real a complex array
  !
  ! the decision whether a serial of parallel FFT is performed
  ! is presently decided by the GRID%RL%NFAST tag
  ! if GRID%RL%NFAST =1 -> serial FFT
  ! if GRID%RL%NFAST =3 -> parallel FFT

  !
  ! transition table used to go from a large to a small grid
  ! or vice versa
  !
  TYPE transit
     INTEGER,POINTER :: IND1(:)   ! fast index transition table
     INTEGER,POINTER :: INDCOL(:) ! column to column transition table
  END TYPE transit


  ! If LPLANE_WISE is set, the data are distributed in real and reciprocal
  ! space plane by plane i.e. one processor holds all elements of
  ! a plane with a specific x index
  ! this  reduces the communication in the FFT considerably
  ! the default for LPLANE_WISE can be set in this file (see below),
  ! or using the flag LPLANE in the INCAR reader
#ifdef plane_wise
  LOGICAL :: LPLANE_WISE=.TRUE.
#else
  LOGICAL :: LPLANE_WISE=.FALSE.
#endif

  ! compatibility modus to vasp.4.4
  ! the flag determine among other things whether the charge at unbalanced
  ! lattice vectors NGX/2+1 are zeroed or not
  LOGICAL :: LCOMPAT

  ! determines the effort FFTW will invest in constructing its "wisdom" (calls to FFTMAKEPLAN)
  ! FFTW_PLAN_EFFORT=0 -> FFTW_ESTIMATE, FFTW_PLAN_EFFORT=1 -> FFTW_MEASURE (default).
  INTEGER :: FFTW_PLAN_EFFORT=1

END MODULE
