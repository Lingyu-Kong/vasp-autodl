#include "symbol.inc"
MODULE nonl_struct_def
  USE prec
  !
  !  structure required to support non local projection operators in recip space
  !
  TYPE nonl_struct
     !only NONL_S
     LOGICAL LRECIP                ! structure set up ?
     INTEGER NTYP                  ! number of types
     INTEGER NIONS                 ! number of ions
     INTEGER NK                    ! kpoint for which CREXP is set up
     INTEGER SELECTED_ION          ! allows to generate a projector for a single ion
     INTEGER, POINTER :: NITYP(:)  ! number of ions for each type
     INTEGER, POINTER :: LMMAX(:)  ! max l-quantum number for each type
     LOGICAL LSPIRAL               ! do we want to calculate spin spirals
!    REAL(q), POINTER, CONTIGUOUS ::QPROJ(:,:,:,:,:) ! projectors in reciprocal space for each k-point
     REAL(q), POINTER ::QPROJ(:,:,:,:,:)
     COMPLEX(q),POINTER, CONTIGUOUS ::CREXP(:,:)  ! phase factor exp (i (G+k) R(ion))
     REAL(q),  POINTER  ::POSION(:,:) ! positions (required for setup)
     REAL(q),  POINTER  ::VKPT_SHIFT(:,:) => NULL()
                                   ! k-point shift for each ion
     COMPLEX(q),POINTER ::CQFAK(:,:)  ! i^l

     ! these arrays allow loops over all ions to be restructured from
     ! nested "over-types" + "over-ions-of-type" to straightforward
     ! loops "over-all-ions" (used in many subroutines).
     INTEGER, POINTER :: ITYP(:)
     INTEGER, POINTER :: LMBASE(:)
  END TYPE nonl_struct
END MODULE nonl_struct_def

MODULE nonlr_struct_def
  USE prec
  USE mpimy
  USE mpi_shmem
  !
  !  structures required to support non local projection operators in real space
  !

  TYPE nonlr_proj
     REAL(q), POINTER :: PSPRNL(:,:,:)
     INTEGER, POINTER :: LPS(:)
  END TYPE nonlr_proj

  !
  !> @brief
  !> The data structure of non local projection operators in real space
  !
  !> @details @ref openmp :
  !> Under OpenMP the members nonlr_struct::nlimax, nonlr_struct::nli,
  !> nonlr_struct::rproj, nonlr_struct::crrexp, and nonlr_struct::nlibase
  !> acquire and additional dimension.
  !
  TYPE nonlr_struct
     !only NONLR_S
     LOGICAL LREAL                 ! structure set up ?
     INTEGER NTYP                  ! number of types
     INTEGER NIONS                 ! number of ions
     INTEGER SELECTED_ION          ! allows to generate a projector for a single ion
     INTEGER IRMAX                 ! maximum number points in sphere
     INTEGER IRALLOC               ! size for allocation =IRMAX*LMDIM*NIONS
     INTEGER NK                    ! kpoint for which CRREXP is set up
     INTEGER, POINTER :: NITYP(:)  ! number of ions for each type
     INTEGER, POINTER :: ITYP(:)   ! type for each ion
     INTEGER, POINTER :: LMAX(:)   ! max l-quantum number for each type
     INTEGER, POINTER :: LMMAX(:)  ! number lmn-quantum numbers for each type
     INTEGER, POINTER ::CHANNELS(:)! number of ln-quantum for each type
     REAL(q), POINTER :: PSRMAX(:) ! real space cutoff
     REAL(q), POINTER :: RSMOOTH(:)! radius for smoothing the projectors around each point
     REAL(q), POINTER :: POSION(:,:) ! positions (required for setup)
     REAL(q), POINTER :: VKPT_SHIFT(:,:)  ! k-point shift for each ion
     TYPE(nonlr_proj), POINTER :: BETA(:) ! a set of structures containing pointers to 
     LOGICAL LSPIRAL               ! do we want to calculate spin spirals

     INTEGER, POINTER :: NLIMAX(:     __NOACC_omp_arg(:)) ! maximum index for each ion
     INTEGER, POINTER :: NLI   (:,:   __NOACC_omp_arg(:)) ! index for gridpoints
     REAL(qn),POINTER, CONTIGUOUS :: RPROJ (:     __NOACC_omp_arg(:)) ! projectors on real space grid
     COMPLEX(q),POINTER, CONTIGUOUS::CRREXP(:,:,: __NOACC_omp_arg(:)) ! phase factor exp (i k (R(ion)-r(grid)))

     ! these arrays allow loops over all ions to be restructured from
     ! nested "over-types" + "over-ions-of-type" to straightforward
     ! loops "over-all-ions" (used in many subroutines).
     INTEGER, POINTER :: LMBASE(:)
     INTEGER, POINTER :: NLIBASE(: __NOACC_omp_arg(:))
#ifdef shmem_rproj
     ! memory segment to store the projectors on the real space grid,
     ! shared between the MPI ranks of COMM_shmem (intra node)
     TYPE (communic), POINTER :: COMM_shmem => NULL()
     TYPE (M_shmem)           :: RPROJ_shmem
#endif
  END TYPE nonlr_struct

  TYPE smoothing_handle
     INTEGER :: N                  ! number of grid points
     REAL(q), POINTER :: WEIGHT(:) ! weight of each grid point
     REAL(q), POINTER :: X1(:), X2(:), X3(:) ! positions of additional grid points in fractional coordinates
  END TYPE smoothing_handle
END MODULE nonlr_struct_def
