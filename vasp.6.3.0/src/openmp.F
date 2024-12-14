#define use_omp_blas
#include "symbol.inc"
      MODULE openmp
      ! overall threading strategy

      !> @ref openmp :
      !> number of threads to be used in \"general\":
      !> by default this is set to the maximum number of available threads,
      !> openmp::omp_nthreads=omp_get_max_threads()
      INTEGER, SAVE :: omp_nthreads   =1

      ! specific routines

      !> @ref openmp :
      !> number of threads used in ::m_alltoall_d_omp
      !> (experimental and unused for now).
      INTEGER, SAVE :: omp_nthreads_alltoall=1

      !> @ref openmp :
      !> number of threads used in the domain decomposition of the
      !> real space projection operators.
      !> Default: openmp::omp_nthreads_nonlr_rspace=openmp::omp_nthreads
      INTEGER, SAVE :: omp_nthreads_nonlr_rspace=1

      !> @ref openmp :
      !> number of threads to be used in connection with OpenACC:
      !> in most instances this involves a distribution over orbitals
      INTEGER, SAVE :: omp_nthreads_acc=1

      !> @ref openmp :
      !> if set to true the domain of the real space projection operators
      !> is decomposed over (x,y)-planes, in a round-robin fashion, if set
      !> to false the projectors are distributed over (z,y)-columns.
      LOGICAL, SAVE :: omp_nonlr_planewise=.TRUE.

      !> @ref openmp :
      !> if set to true ::fftmakeplan or ::fftmakeplan_mpi will call
      !> dfftw_init_threads and set openmp::omp_dfftw_init_threads=.FALSE.
      LOGICAL, SAVE :: omp_dfftw_init_threads=.TRUE.

      CONTAINS

!> @details @ref openmp :
!> this subroutine reads
!>   openmp::omp_nthreads
!>   openmp::omp_nthreads_alltoall
!>   openmp::omp_nthreads_nonlr_rspace
!>   openmp::omp_nonlr_planewise
!>
!> When none of these tags are specified in the INCAR file,
!>   openmp::omp_nthreads              = omp_get_max_threads() (i.e., OMP_NUM_THREADS),
!>   openmp::omp_nthreads_nonlr_rspace = openmp::omp_nthreads
!
      SUBROUTINE INIT_OMP(IO)
      USE prec
      USE base
      USE reader_tags
      TYPE (in_struct) IO
      ! local
      INTEGER IDUM, N, IERR
      REAL(q) RDUM
      COMPLEX(q) CDUM
      LOGICAL LDUM
      CHARACTER(1) CHARAC
#ifdef _OPENMP
      INTEGER, EXTERNAL :: OMP_GET_MAX_THREADS


      ! prevent nesting in openMP
!!#ifdef _OPENMP45
      CALL OMP_SET_MAX_ACTIVE_LEVELS(1)
!!#else
!!      CALL OMP_SET_NESTED(.FALSE.)
!!#endif
 
      omp_nthreads=OMP_GET_MAX_THREADS()

      IF (omp_nthreads<2) THEN
         omp_nthreads=1 ; RETURN
      ENDIF

      omp_nthreads_alltoall=1
      CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5,'NTHREADS_A2A',omp_nthreads_alltoall,IERR)

      IF (omp_nthreads_alltoall<1)            omp_nthreads_alltoall=1
      IF (omp_nthreads_alltoall>omp_nthreads) omp_nthreads_alltoall=omp_nthreads

      omp_nthreads_nonlr_rspace=omp_nthreads
      CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5,'NTHREADS_RPROJ',omp_nthreads_nonlr_rspace,IERR)

      IF (omp_nthreads_nonlr_rspace<1.OR.omp_nthreads_nonlr_rspace>omp_nthreads)  &
     &   omp_nthreads_nonlr_rspace=omp_nthreads
#ifdef _OPENACC
      omp_nthreads_nonlr_rspace=1
#endif
#endif
      RETURN
      END SUBROUTINE INIT_OMP

      END MODULE openmp

      SUBROUTINE OMP_DCOPY(N,DX,INCX,DY,INCY)
      USE prec
      INTEGER :: N,INCX,INCY
      REAL(q) :: DX(*),DY(*)
#ifdef use_omp_blas
      IF (INCX==1.AND.INCY==1) THEN
      ! both increments equal one
!DIR$ IVDEP
!OCL NOVREC
!$OMP PARALLEL DO
         DO I=1,N
            DY(I)=DX(I)
         ENDDO
!$OMP END PARALLEL DO
         RETURN
      ENDIF
      IF (INCX==2.AND.INCY==1) THEN
      ! specialize complex -> real copy
!DIR$ IVDEP
!OCL NOVREC
!$OMP PARALLEL DO
         DO I=1,N
            II=2*I-1
            DY(I)=DX(II)
         ENDDO
!$OMP END PARALLEL DO
         RETURN
      ENDIF
      ! one or both increments unequal one
!$OMP PARALLEL DO PRIVATE(IX,IY)
      DO I=1,N
         IX=(I-1)*INCX+1
         IY=(I-1)*INCY+1
         DY(IY)=DX(IX)
      ENDDO
!$OMP END PARALLEL DO
      RETURN
#endif
      CALL DCOPY(N,DX(1),INCX,DY(1),INCY)

      END SUBROUTINE OMP_DCOPY


      SUBROUTINE OMP_DSCAL(N,DA,DX,INCX)
      USE prec
      INTEGER :: N,INCX
      REAL(q) :: DX(*)
      REAL(q) :: DA
      IF (DA==1._q) RETURN
#ifdef use_omp_blas
      IF (INCX==1) THEN
!DIR$ IVDEP
!OCL NOVREC
!$OMP PARALLEL DO
         DO I=1,N
            DX(I)=DX(I)*DA
         ENDDO
!$OMP END PARALLEL DO
      ELSE
!$OMP PARALLEL DO PRIVATE(IX)
         DO I=1,N
            IX=(I-1)*INCX+1
            DX(IX)=DX(IX)*DA
         ENDDO
!$OMP END PARALLEL DO
      ENDIF
      RETURN
#endif
      CALL DSCAL(N,DA,DX(1),INCX)

      END SUBROUTINE OMP_DSCAL


      SUBROUTINE OMP_DSCPY(N,DA,DX,INCX,DY,INCY)
      USE prec
      INTEGER :: N,INCX,INCY
      REAL(q) :: DX(*),DY(*)
      REAL(q) :: DA
      IF (DA==1._q) THEN
         CALL OMP_DCOPY(N,DX(1),INCX,DY(1),INCY)
      ELSE
#ifdef use_omp_blas
         IF (INCX==1.AND.INCY==1) THEN
         ! both increments equal one
!DIR$ IVDEP
!OCL NOVREC
!$OMP PARALLEL DO
            DO I=1,N
               DY(I)=DX(I)*DA
            ENDDO
!$OMP END PARALLEL DO
            RETURN
         ENDIF
         IF (INCX==2.AND.INCY==1) THEN
         ! specialze complex -> real
!DIR$ IVDEP
!OCL NOVREC
!$OMP PARALLEL DO
            DO I=1,N
               II=2*I-1
               DY(I)=DX(II)*DA
            ENDDO
!$OMP END PARALLEL DO
            RETURN
         ENDIF
         ! one or both increments unequal one
!$OMP PARALLEL DO PRIVATE(IX,IY)
         DO I=1,N
            IX=(I-1)*INCX+1
            IY=(I-1)*INCY+1
            DY(IY)=DX(IX)*DA
         ENDDO
!$OMP END PARALLEL DO
         RETURN
#else
         CALL DCOPY(N,DX,INCX,DY,INCY)
         CALL DSCAL(N,DA,DY,INCY)
#endif
      ENDIF
      END SUBROUTINE OMP_DSCPY


      SUBROUTINE OMP_DAXPY(N,DA,DX,INCX,DY,INCY)
      USE prec
      INTEGER :: N,INCX,INCY
      REAL(q) :: DX(*),DY(*)
      REAL(q) :: DA
      IF (DA==0._q) RETURN
#ifdef use_omp_blas 
      IF (INCX==1.AND.INCY==1) THEN
      ! both increments equal one
         IF (DA==1._q) THEN
!DIR$ IVDEP
!OCL NOVREC
!$OMP PARALLEL DO
            DO I=1,N
               DY(I)=DY(I)+DX(I)
            ENDDO
!$OMP END PARALLEL DO
         ELSE
!DIR$ IVDEP
!OCL NOVREC
!$OMP PARALLEL DO
            DO I=1,N
               DY(I)=DY(I)+DX(I)*DA
            ENDDO
!$OMP END PARALLEL DO
         ENDIF
         RETURN
      ENDIF
      IF (INCX==2.AND.INCY==1) THEN
      ! specialize complex -> real
         IF (DA==1._q) THEN
!DIR$ IVDEP
!OCL NOVREC
!$OMP PARALLEL DO
            DO I=1,N
               II=2*I-1
               DY(I)=DY(I)+DX(II)
            ENDDO
!$OMP END PARALLEL DO
         ELSE
!DIR$ IVDEP
!OCL NOVREC
!$OMP PARALLEL DO
            DO I=1,N
               II=2*I-1
               DY(I)=DY(I)+DX(II)*DA
            ENDDO
!$OMP END PARALLEL DO
         ENDIF
         RETURN
      ENDIF
      ! one or both increments unequal one
      IF (DA==1._q) THEN
!$OMP PARALLEL DO PRIVATE(IX,IY)
         DO I=1,N
            IX=(I-1)*INCX+1
            IY=(I-1)*INCY+1
            DY(IY)=DY(IY)+DX(IX)
         ENDDO
!$OMP END PARALLEL DO
      ELSE
!$OMP PARALLEL DO PRIVATE(IX,IY)
         DO I=1,N
            IX=(I-1)*INCX+1
            IY=(I-1)*INCY+1
            DY(IY)=DY(IY)+DX(IX)*DA
         ENDDO
!$OMP END PARALLEL DO
      ENDIF
      RETURN
#endif
      CALL DAXPY(N,DA,DX(1),INCX,DY(1),INCY)

      END SUBROUTINE OMP_DAXPY

!> @page openmp OpenMP
!! @tableofcontents
!! @section general All changes
!! All datatypes, subroutines, and functions that have been changed
!! to enable parallelization under OpenMP.
!!
!! \li ::apply_gfac
!! \li ::apply_gfac_2_encuts
!! \li ::apply_gfac_der
!! \li ::apply_gfac_encuts
!! \li ::apply_gfac_exchange
!! \li ::apply_gfac_exchange_2
!! \li augfast::calc_dllmm_trans
!! \li augfast::depsum_two_bands_rholm_trace
!! \li augfast::depsum_vector
!! \li augfast::setup_aug_des
!! \li broyden::brmix
!! \li chi_base::add_responsefunction_cache
!! \li chi_base::add_responsefunction_int
!! \li chi_base::add_xi
!! \li chi_base::add_xi_spectral
!! \li chi_base::allocate_responsefun_cache
!! \li chi_base::clean_responsefunction_cache
!! \li chi_base::clean_responsefunction_int
!! \li chi_base::deallocate_responsefun
!! \li chi_base::responsefunction
!! \li chi_base::screened_two_electron_cached
!! \li choleski::orthch
!! \li ::crrexp_mul_work_gadd
!! \li ::crexp_mul_wave
!! \li ::complex_add
!! \li ::cmplx_cmplx_cmplx_mul
!! \li ::cmplx_cmplx_real_mul
!! \li ::cmplx_real_real_mul
!! \li ::dllmm_kernel
!! \li ::eccp_nl
!! \li ::eccp_nl_fock
!! \li ::exchange_gfac
!! \li ::fexcg_
!! \li ::fftbas
!! \li ::fftbas_plan_mpi
!! \li ::fftbrc
!! \li ::fftext_mpi
!! \li ::ffthamil
!! \li ::fftmakeplan
!! \li ::fftmakeplan_mpi
!! \li ::fftwav
!! \li ::fftwav_mpi
!! \li fock::allocate_fock_handle
!! \li fock::deallocate_fock_handle
!! \li fock::fock_acc
!! \li fock::fock_charge_mu
!! \li fock::fock_force
!! \li fock::fock_handle
!! \li fock::xc_fock_reader
!! \li fock_dbl::fock_all_dblbuf
!! \li force::forhar
!! \li force::forloc
!! \li ::ggaall_grid
!! \li greens_real_space::transpose_g_r_response
!! \li greens_real_space::transpose_g_proj_response
!! \li greens_real_space::transpose_r_g
!! \li greens_real_space::transpose_r_proj
!! \li greens_real_space_k::calculate_response_super
!! \li greens_real_space_k::fft_g_super
!! \li hamil::eccp
!! \li hamil::eccp_tau
!! \li hamil::eccp_vec
!! \li hamil::hamiltmu
!! \li hamil::hamiltmu_c
!! \li hamil::hamiltmu_tau
!! \li hamil::hamiltmu_vec
!! \li hamil::setup_precond
!! \li hamil::simple_precond
!! \li hamil_lrf::hamiltmu_commutator
!! \li hamil_lrf::hamiltmu_lrf
!! \li hamil_lrf::lrf_commutator
!! \li hamil_lrf::lrf_hamil
!! \li ::kinhamil
!! \li ::kinhamil_c
!! \li ::kinhamil_tau
!! \li ::kinhamil_vec
!! \li main_mpi::init_mpi
!! \li main_mpi::wrt_distr
!! \li ::map_backward
!! \li ::map_forward
!! \li ::map_gather
!! \li ::map_scatter
!! \li mpimy::m_init
!! \li nonl::phase
!! \li nonl::proj1
!! \li nonl::projxyz
!! \li nonl::projxyz_wa
!! \li nonl::spher
!! \li nonl::strenl
!! \li nonlr::fornlr
!! \li nonlr::nonlr_alloc
!! \li nonlr::nonlr_alloc_crrexp
!! \li nonlr::nonlr_dealloc
!! \li nonlr::nonlr_set_single_ion
!! \li nonlr::phaser
!! \li nonlr::phaserr
!! \li nonlr::phaser_hf
!! \li nonlr::real_optlay
!! \li nonlr::real_optlay_grid
!! \li nonlr::rnlpr
!! \li nonlr::rpro1
!! \li nonlr::rpro1_hf
!! \li nonlr::rpromu
!! \li nonlr::rpromu_hf
!! \li nonlr::rspher_all
!! \li nonlr_struct_def::nonlr_struct
!! \li nonl_high::w1_projall
!! \li openmp::init_omp
!! \li ::overl
!! \li ::overl1
!! \li ::overl1_c
!! \li ::overl_fock
!! \li pawfock::coloumb_4term
!! \li pawfock::coloumb_4term_ps
!! \li pawfock_inter::ntyp_slater
!! \li pawfock_inter::s
!! \li pawm::set_dd_paw
!! \li ::pw_charge
!! \li ::pw_charge_trace
!! \li ::pw_norm_with_metric
!! \li ::racc0
!! \li ::racc0mu
!! \li ::racc0mu_hf
!! \li ::racc0_hf
!! \li radial::rad_lda_xc
!! \li ::real_add
!! \li ::real_cmplx_cmplx_mul
!! \li ::real_real_cmplx_mul
!! \li ::real_real_real_mul
!! \li ::rholm_kernel
!! \li ::rholm_kernel_aux
!! \li ::rholm_kernel_dgemm
!! \li ::rholm_one_center_kernel
!! \li rmm_diis::eddrmm
!! \li rmm_diis_lr::linear_response_diis
!! \li rot::edwav
!! \li rs_greensfunc_kernel::set_rsgf_all
!! \li rs_greensfunc_kernel::rs_coulomb_green_func
!! \li rs_greensfunc_kernel::rs_greensfunc_type
!! \li screened_2e::determine_slot
!! \li screened_2e::determine_slot_inter
!! \li screened_2e::determine_slot_inter_weight
!! \li screened_2e::integrate_w_2e_simple
!! \li screened_2e::integrate_w_2e_spectral
!! \li screened_2e::integrate_w_2e_spectral_imag
!! \li screened_2e::qp_shift
!! \li screened_2e::qp_shift_pade
!! \li ::setdij_
!! \li ::start_profiling
!! \li ::stop_profiling
!! \li subrot::eddiag
!! \li ::truncate_high_frequency_one
!! \li twoelectron4o::apply_phase
!! \li us::depsum
!! \li us::fordep
!! \li ::vhamil
!! \li ::vhamil_trace
!! \li wave::delwav_omp
!! \li wave::newwav_omp
!! \li wave_cacher::allocate_cacher
!! \li wave_cacher::deallocate_cacher
!! \li wave_cacher::find_new_cacher
!! \li wave_cacher::remove_cacher
!! \li wave_cacher::store_cacher
!! \li wave_cacher::store_gw_acc
!! \li wave_high::cnorma
!! \li wave_high::w1_copy
!! \li wave_high::w1_daxpy
!! \li wave_high::w1_dot
!! \li wave_high::w1_dscal
!! \li wave_high::w1_gaxpy
!! \li ::work_mul_crexp
!! \li xi::calculate_xi
!!
!! @section simple Simple loop parallelism
!! In the following cases only basic loop parallelism under OpenMP was added.
!! Sometimes in conjunction with a reduction operation, and sometimes
!! conditional upon a certain loop size.
!! These loops often run over all points in a real or reciprocal space grid,
!! or over all basis vectors in an orbital.
!!
!! \li ::apply_gfac
!! \li ::apply_gfac_2_encuts
!! \li ::apply_gfac_der
!! \li ::apply_gfac_encuts
!! \li ::apply_gfac_exchange
!! \li ::apply_gfac_exchange_2
!! \li augfast::depsum_vector
!! \li broyden::brmix
!! \li ::crrexp_mul_work_gadd
!! \li ::crexp_mul_wave
!! \li ::complex_add
!! \li ::dllmm_kernel
!! \li ::cmplx_cmplx_cmplx_mul
!! \li ::cmplx_cmplx_real_mul
!! \li ::cmplx_real_real_mul
!! \li ::eccp_nl
!! \li ::eccp_nl_fock
!! \li ::exchange_gfac
!! \li ::fexcg_
!! \li ::ffthamil
!! \li force::forhar
!! \li force::forloc
!! \li ::ggaall_grid
!! \li ::kinhamil
!! \li nonl::phase
!! \li nonl::strenl
!! \li pawfock::coloumb_4term
!! \li pawfock::coloumb_4term_ps
!! \li ::pw_charge
!! \li ::pw_charge_trace
!! \li ::pw_norm_with_metric
!! \li radial::rad_lda_xc
!! \li ::real_add
!! \li ::real_cmplx_cmplx_mul
!! \li ::real_real_cmplx_mul
!! \li ::real_real_real_mul
!! \li ::rholm_kernel
!! \li ::rholm_kernel_aux
!! \li ::rholm_kernel_dgemm
!! \li ::rholm_one_center_kernel
!! \li rs_greensfunc_kernel::rs_coulomb_green_func
!! \li ::truncate_high_frequency_one
!! \li twoelectron4o::apply_phase
!! \li ::vhamil
!! \li ::vhamil_trace
!! \li ::work_mul_crexp
