#include "symbol.inc"

MODULE fock_glb

  USE prec
  USE wave_struct_def
  USE nonlr_struct_def
#ifdef shmem_bcast_buffer
  USE mpi_shmem
#endif
#ifdef _OPENMP
  USE omp_lib_kinds, ONLY : OMP_LOCK_KIND
#endif
  IMPLICIT NONE

  !> use HF Hamiltonian
  LOGICAL,SAVE :: LHFCALC=.FALSE.

  !> precision for HF type calculations
  !> L  low    (coarse grid for HF, normal augmentation charge)
  !> F  fast   (coarse grid for HF, soft augmentation charge)
  !> N  normal (normal grid for HF, normal augmentation charge)
  !> A  accurate (normal grid for HF, soft augmentation charge)
  CHARACTER(12), SAVE ::    SZPRECFOCK='      '
  REAL(q) :: ENCUTFOCK  !< test version
                        !< apply spherical cutoff for HF part
                        !< recommended default ENCUT for PRECFOCK = Fast

  !> use EXX-optimized effective potential method
  !> usually this implies that the HF Hamiltonian
  !> is not apllied in the calculations of H phi
  INTEGER,SAVE :: EXXOEP=0

  !> this flag forces the use the the Fock Hamiltonian when calculating
  !> H phi
  !> throughout all routines (even if EXXOEP not zero)
  LOGICAL,SAVE :: LHFCALC_FORCE=.FALSE.

  !> LSYMGRAD determines whether symmetry is used to restore
  !> the proper symmetry of the gradient
  LOGICAL,SAVE :: LSYMGRAD=.FALSE.

  !> amount of exact exchange (the remainder is treated in the LDA framework)
  !> amount of exact exchange (the remainder is treated in the LDA framework)
  !> throughout all routines
  REAL(q),SAVE :: AEXX=0

  !> use Hartree Fock only in one center AE terms
  LOGICAL,SAVE :: LHFONE=.FALSE.

  !> screening length in screened exchange methods
  REAL(q),SAVE :: HFSCREEN=0

  !> screening length in screened exchange methods for correlation
  REAL(q),SAVE :: HFSCREENC=0

  !> cutoff for HF kernel applicable to 0 and 3 dimensional systems (molecules)
  REAL(q),SAVE :: HFRCUT=0

  !> smoothed cutoff function applied in correlated calculations
  REAL(q),SAVE :: HFRCUT_SMOOTH=1._q/200._q

  !> parameter for the determination of the convergence corrections
  REAL(q),SAVE :: HFALPHA

  !> this variable allows to reduce the exchange operator to a coarse
  !> grid q which is reduced in the first second and third direction
  !> by the values supplied by the user
  INTEGER, SAVE :: NKREDX=1, NKREDY=1, NKREDZ=1

  !> this variable allows to shift the grid on which the Fock
  !> operator is evaluated
  !> works only in combination with at least one NKRED set to 2
  LOGICAL :: HFKIDENT=.FALSE.

  !> this variable allows to shift the grid on which the Fock
  !> operator is evaluated
  !> works only in combination with at least one NKRED set to 2
  LOGICAL :: SHIFTRED=.FALSE.

  !> this variable allows to use only every second point in the HF
  !> grid, gamma point is included
  !> sort of similar to NKRED, but resulting in a net downsampling by 2
  LOGICAL :: EVENONLY =.FALSE.

  !> similar to the above flag, but it allows to use odd grid points only
  !> that is gamma is not included
  LOGICAL :: ODDONLY =.FALSE.

  !> maximum L quantum number to which the charge density
  !> is augmented on the plane wave grid
  !> this has quite some performance impact, but unfortunately
  !> LMAX_FOCK can hardly be ever decreased without noticeable
  !> loss in precision
  INTEGER, SAVE :: LMAX_FOCK

  !> maximum L quantum number to which the charge density
  !> is accurately augmented on the plane wave grid
  !> VASP usually restores all moments (monopole, dipole, quadrupole etc.)
  !> on the plane wave grid).
  !> In some routines like GW, RPA correlation, the one center terms
  !> are not implemented and one would like to restore the AE charge density
  !> very accurately on the plane wave grid.
  !> This has very significant performance impact, and since the one
  !> center terms correct for any difference between PS and AE
  !> wavefunction, use this flag only for GW type calculations
  !> NMAX_FOCKAE determines how many functions are used for each
  !> channel. Only up to two are supported
  INTEGER, SAVE :: LMAX_FOCKAE
  INTEGER, SAVE :: NMAX_FOCKAE
  !> this flag forces the use of LMAXFOCKAE_IN_DFT in every part of the code
  !> including DFT
  LOGICAL, SAVE :: LFOCKAEDFT=.FALSE.

  REAL(q),ALLOCATABLE, SAVE :: QMAX_FOCKAE(:)
    
  !> LRSCOR=.TRUE.  Range separated LDA correlation
  !> LRSCOR=.FALSE. Complete LDA correlation
  LOGICAL       :: LRSCOR      

  !> LRHFCALC long range HF interaction only
  !> the default is short range HF interaction only
  !> the variable should be identical to LUSE_LONGRANGE_HF in xclib
  LOGICAL       :: LRHFCALC

  !> temporarily switch off the HF treatment
  !
  LOGICAL       :: LSTACK_FOCK , LHFCALC_STACK

  !> L_THOMAS_FERMI Thomas Fermi screening in HF exchange
  !> should be identical to LUSE_THOMAS_FERMI in fock.F
  LOGICAL        :: L_THOMAS_FERMI

  !> LUSE_MODEL_HF short range full HF
  !> long range a fraction of the full HF specified by AEXX
  LOGICAL       :: L_MODEL_HF

  !> calculate four orbital integrals
  INTEGER        :: FOURORBIT
  REAL(q), SAVE  :: ENCUT4O

  !> this matrix stores the transformation matrix that
  !> is required to go from the occupancies RHO(lm, l'm') to
  !> the onsite L dependent occupancies RHO(LM)
  !> the matrix differs for each type
  REAL(q), POINTER, SAVE, CONTIGUOUS :: TRANS_MATRIX_FOCK(:,:,:,:) 

  !> this structures stores the data required to perform
  !> a fast augmentation of the pseudo charge density within
  !> the HF related routines
  TYPE (nonlr_struct), SAVE, TARGET :: FAST_AUG_FOCK

  !> structure to store the data layout of the one-center charge densities
  !> in the HF routines
  !> this must match the FAST_AUG_FOCK data layout
  TYPE (wavedes1), SAVE :: AUG_DES

  !> this structure stores the wavefunctions descriptor
  !> for the HF related routines
  !> usually it is equivalent to WDES, but it is possible
  !> to use a coarser FFT grid for all HF related routines
  TYPE (wavedes),SAVE, POINTER :: WDES_FOCK

  !> WDES_FOCK differs from WDES
  LOGICAL :: LWDES_FOCK=.FALSE.

  !> this structure stores the related grid structure
  !> for the 3d-FFT of orbitals
  !> usually simply a pointer to GRID
  !> in the parallel cases, this grid might contain
  !> only the minimal required grid points to reduce computational cost
  !> NOTE: for Gamma only in real space the arrays must be declared as COMPLEX
  !> as for all orbital FFT's in VASP
  TYPE (grid_3d), POINTER :: GRID_FOCK

  !> store the 3d structure to perform a 3d FFT of charges and potentials
  !> within HF related routines
  !> stores all grid points in real and reciprocal space
  !> NOTE: for Gamma only in real space the arrays must be declared as GDEF
  !> see also mgrid.F GEN_GRID_HF for further comments
  TYPE (grid_3d) :: GRIDHF

  !> leading dimensions of CDIJ, or CRHODE
  INTEGER :: LMDIM_FOCK

  !> use model GW
  INTEGER :: MODEL_GW=0

  !> static dielectric constant for GW
  REAL(q) MODEL_GW_EPS0

  !> parameter alpha for model GW
  !> although the Bechstedt et al. Solid State Comm. 84, 765 suggests to use 1.5_q
  !> it seems most Bechstedt people use 1.0
  REAL(q) :: MODEL_GW_ALPHA=1.0_q

  !> FSG_STORE stores the convergence correction
  !> array is set up by the SETUP_FOCK routines
  !> for each k-point one value is stored
  !
  REAL(q), POINTER :: FSG_STORE(:)

  !> NBANDSGWLOW_FOCK allows to exclude low energy states from the HF
  !> Hamiltonian
  !
  INTEGER :: NBANDSGWLOW_FOCK

  !> blocking factor in FOCK_ACC and FOCK_FORCE
  !
  INTEGER :: NBLOCK_FOCK=64

  !> this variable determines whether we use the Adeptively Compressed Exchange (ACE)
  LOGICAL, SAVE :: LFOCKACE=.FALSE.

  !> this handle in principle allows to store all quantities
  !> required to calculate the accelerations on the wavefunction in HF case
  !
  TYPE fock_handle
     TYPE (wavedes) :: WDES
     INTEGER  :: LMDIM
     COMPLEX(q),POINTER :: CXI(:,:,:)   !< stores acceleration in real space for each band
     GDEF, POINTER :: CKAPPA(:,:,:)     !< stores NL accelerations
     GDEF, POINTER, CONTIGUOUS :: CDLM(:)
     GDEF, POINTER :: CDIJ(:,:,:,:)
     TYPE (wavedes1) :: WDESK           ! 
     TYPE (wavefun1) :: W1              !< really only a dummy W1%CPROJ=> CDLM
     TYPE (wavefun1) :: WQ              !< wavefunction

!$   ! locks for CXI and CKAPPA (needed in STORE_GW_ACC)
!$   INTEGER(KIND=OMP_LOCK_KIND), POINTER :: OMP_LCK1(:),OMP_LCK2(:)
  END TYPE fock_handle

#ifdef shmem_bcast_buffer
  !
  !> these shared memory segments replace parts of the on-node
  !> communication in FOCK_ACC
  !
  TYPE (M_shmem) :: FOCK_SHM_BCST1,FOCK_SHM_BCST2,FOCK_SHM_BCST3
  LOGICAL, SAVE  :: LFOCK_SHM_ALLOC=.FALSE.
#endif

END MODULE fock_glb
