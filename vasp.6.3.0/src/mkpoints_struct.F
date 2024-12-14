#include "symbol.inc"
MODULE mkpoints_struct_def
  USE prec
! logical variable to determine whether a full k-point grid is required
  LOGICAL,SAVE :: LFULL_KPOINTS=.FALSE.
! preferably use time inversion symmetry to generate orbitals at k-points -k
  LOGICAL,SAVE :: LTIME_INVERSION=.FALSE.

! classifies the mode of k-points generation (set in RD_KPOINTS)
  ENUM, BIND(C)
     ENUMERATOR :: ModeUnknown, MonkhorstPack, kSpacing, RkLength, GeneratingLattice, LineMode, ExplicitList
  END ENUM

!
!  structure to store the k-point grid in the full wedge of the Brillouin zone
!
  TYPE skpoints_full
     INTEGER NKPTS                   ! actual number of k-points
     INTEGER NKPTS_NON_ZERO          ! number of k-points with non zero weight
     REAL(q),POINTER :: VKPT(:,:)    ! coordinate of k-point
     REAL(q),POINTER :: WTKPT(:)     ! symmetry weight-factor for each k-point
     INTEGER,POINTER :: NEQUIV(:)    ! equivlist full kp-ibzkpt
     INTEGER,POINTER :: IROTOP(:,:,:)! integer rotation matrix in reciprocal space 
     INTEGER,POINTER :: ISYMOP(:,:,:)! integer rotation matrix in real space
     REAL(q),POINTER :: TRANS(:,:)   ! possible nontrivial translation
     LOGICAL,POINTER :: LSHIFT(:)    ! translation required
     INTEGER,POINTER :: SPINFLIP(:)  ! spin flip required
     LOGICAL,POINTER :: LINV(:)      ! mirror point?
     COMPLEX(q),POINTER:: PHASE(:,:) ! phase shift for PW components 
     INTEGER,POINTER :: NG_INDEX(:,:)! used in parallel version to reindex global PW coefficients to local coefficients
     INTEGER,POINTER :: NOP(:)       ! index of symmetry operation
     INTEGER,POINTER :: ROTMAP(:,:)  ! map indexing atoms that are taken into each other when 
                                     ! the symmetry operation is applied
     INTEGER :: NKPX, NKPY, NKPZ     ! integer division along rec. lattice vectors for generating the mesh
     REAL(q) :: B(3,3)               ! generating lattice vector
     REAL(q) :: SHIFT(3)             ! shift of the generating lattice
     REAL(q) :: VKPTMINDIST2         ! squared minimal distance between two k-points 
     INTEGER :: MODE=ModeUnknown     ! stores the k-points generation format specification
#ifndef oldsym
     COMPLEX(q),POINTER :: RSSYMOP(:,:,:) ! Spin rotation matrix in real space  
#endif
     INTEGER, POINTER :: MAP_INTO_BZ(:,:,:)
  END TYPE skpoints_full

! pointer to the currently used k-points
  TYPE (skpoints_full), SAVE, POINTER :: KPOINTS_FULL=>NULL()
!
! the rotation handle allows a quick rotation of the wavefunction character
! it hashes all required entries
!
  TYPE rotation_handle
     INTEGER :: LDIM                 ! set to MAXL(T_INFO%NTYP,P)
     INTEGER :: MMAX                 ! set to 2*LMDIM+1
     INTEGER :: NK                   ! k-points to which one needs to rotate
     REAL(q),POINTER :: SL(:,:,:) => NULL()   ! rotation matrix
     INTEGER, POINTER :: NPRO_NI(:) => NULL() ! pointer to the position where the wavefunction character is stored for 
                                     ! a particular ion
     INTEGER, POINTER :: NPS(:) => NULL() ! number of channels for each pseudo-potential type (= P%LMAX)
     INTEGER, POINTER :: LPS(:,:)=> NULL() ! l-quantum number of all pseudo-potential channels
  END TYPE rotation_handle

  REAL(q), ALLOCATABLE, SAVE :: WEIGHT_K_POINT_PAIR_SMALL_GROUP(:,:)
  LOGICAL :: LSYMGRAD_SAVE=.FALSE.

!
!  structure required for kpoints generation
!
  TYPE kpoints_struct
!only  KPOINTS
     INTEGER :: NKDIM             ! maximal number of k-points
     INTEGER :: NKPTS             ! actual number of k-points
     INTEGER :: NKPTS_NON_ZERO    ! number of k-points with non zero weight
     ! it is usually save to set this to NKPTS, then some routines might do extra work
     ! NKPTS_NON_ZERO is presently only handled and used by the GW routines
     REAL(q),POINTER:: VKPT(:,:)  ! coordinate of k-point
     REAL(q),POINTER:: WTKPT(:)   ! symmetry weight-factor for each k-point
     INTEGER,POINTER:: IDTET(:,:) ! link list for tetrahedron
     REAL(q) :: VOLWGT            ! volume weight for each tetrahedron
     INTEGER :: NTET              ! number of tetrahedrons
     LOGICAL :: LTET              ! use tetrahedron method ?
     INTEGER :: ISMEAR            ! type of smearing
     REAL(q) :: SIGMA             ! type of smearing
     REAL(q) :: EMIN              ! minimal E for DOS
     REAL(q) :: EMAX              ! maximal E for DOS
     REAL(q) :: EFERMI            ! maximal E for DOS
     INTEGER :: NKPX,NKPY,NKPZ    ! integer division along rec. lattice vectors for generating the mesh
     REAL(q) :: B(3,3)            ! generating lattice for k-point mesh
     REAL(q) :: SHIFT(3)          ! shift of the generating lattice
     REAL(q) :: SPACING           ! spacing of k-point grid in A-1
     LOGICAL :: LGAMMA            ! gamma point is included
     LOGICAL :: LKBLOWUP          ! include "blow-up phase" in automatic k-points generation
     INTEGER :: MODE=ModeUnknown  ! stores the k-points generation format specification
     CHARACTER*40  SZNAMK         ! name of k-points file
    
  END TYPE kpoints_struct

!
!  structure to store all quantities that are required to go
!  from one k-point grid in the IRZ to another one 
!  i.e. if symmetry is reduced
!
  TYPE skpoints_trans
     INTEGER NKPTS                   ! not used
     LOGICAL,POINTER   :: LINV(:)    ! mirror point?
     LOGICAL,POINTER   :: LSHIFT(:)  ! phase shift required
     LOGICAL,POINTER   :: LSHIFT_G(:)! G vectors shifted
     COMPLEX(q),POINTER:: CPHASE(:,:)! phase shift for PW components 
     INTEGER,POINTER   :: NK_OLD(:)  ! index of original k-point in the IRZ
     INTEGER,POINTER   :: NINDPW(:,:)! index of plane wave component after application of symmetry
     INTEGER,POINTER  ::ISYMOP(:,:,:)! rotation matrix from k-point in the IRZ to full k-point
     INTEGER,POINTER   :: ROTMAP(:,:)! map indexing atoms that are taken into each other when 
                                     ! the symmetry operation is applied
     INTEGER,POINTER   :: SPINFLIP(:)! spin flip required
#ifndef oldsym
     COMPLEX(q),POINTER :: RSSYMOP(:,:,:) ! Spin rotation matrix in real space
#endif
  END TYPE skpoints_trans

! gK: the previous two variables are global. This raises major issues if symmetry is lowered stepwise.
! This occurs for the combination of finite differences and RPA calculations
! to resolve this issue, the current KPOINTS_ORIG can be pushed (PUSH_KPOINTS_ORIG) onto a stack and later
! restored using pop (POP_KPOINTS_ORIG)
! to implement a stack  
! I need a stack structure, since F90 does not have arrays of pointers
  TYPE kpoints_stack_type
     TYPE (kpoints_struct), POINTER :: KPOINTS_ORIG => NULL()
     TYPE (skpoints_full), POINTER  :: KPOINTS_FULL_ORIG => NULL()
  END TYPE kpoints_stack_type
  
!-----hard limits for k-point generation package
!     NTETD  is the initial number of tetrahedra which can be
!            stored when using the tetrahedron integration.
!            When needed this is increased in TETIRR
  INTEGER, PARAMETER :: NKDIMD=20000,NTETD=90000
      
! pointer to the currently used k-points
  TYPE (kpoints_struct),SAVE,POINTER :: KPOINTS_     

! logical variable to determine whether also all k-points
! which are the difference between two other k-points are used
! in the full k-point mesh (this is required for GW calculations)
  LOGICAL, SAVE :: LSHIFT_KPOINTS=.FALSE.

! logical variable that tracks whether KBLOWUP was read from INCAR or not
  LOGICAL, SAVE :: LKBLOWUP_WAS_READ=.FALSE.

! pointer to the k-points generated initially
! (it is possible the k-points change during a run, if the
!  e.g. the symmetry is lowered)
  TYPE (kpoints_struct),SAVE,POINTER :: KPOINTS_ORIG

! pointer to the full k-point grid generated initially
  TYPE (skpoints_full),SAVE,POINTER :: KPOINTS_FULL_ORIG

END MODULE mkpoints_struct_def
