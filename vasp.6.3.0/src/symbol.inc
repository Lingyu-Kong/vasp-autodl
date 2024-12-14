!-------- to be costumized by user (usually done in the makefile)-------
!#define vector              compile for vector machine
!#define essl                use ESSL instead of LAPACK
!#define single_BLAS         use single prec. BLAS

!#define wNGXhalf            gamma only wavefunctions (X-red)
!#define wNGZhalf            gamma only wavefunctions (Z-red)

!#define NGXhalf             charge stored in REAL array (X-red)
!#define NGZhalf             charge stored in REAL array (Z-red)
!#define NOZTRMM             replace ZTRMM by ZGEMM
!#define MPI                 compile for parallel machine with MPI
!------------- end of user part --------------------------------
#if defined(USELIBXC5) && ! defined(USELIBXC)
#define USELIBXC
#endif
#ifdef single_BLAS
!
!   use single precission BLAS and LAPACK calls (usually on Cray)
!
#define ZGEMM  CGEMM
#define ZGEMV  CGEMV
#define ZDSCAL CSSCAL
#define ZDOTC  CDOTC
#define ZAXPY  CAXPY
#define ZHEEV  CHEEV
#define ZHEGV  CHEGV
#define ZTRMM  CTRMM
#define ZTRTRI CTRTRI
#define ZPOTRF CPOTRF
#define PZPOTRF PCPOTRF
#define PZTRTRI PCTRTRI
#define PZHEEVX PCHEEVX
#define ZCOPY   CCOPY
#define ZGETRF  CGETRF
#define ZGETRS  CGETRS
#define ZGEFA   CGEFA
#define ZGECO   CGECO
#define ZGEDI   CGEDI
#define ZHEEVX  CHEEVX

#define DGEMM  SGEMM
#define DGEMV  SGEMV
#define DGEGV  SGEGV
#define DDOT   SDOT
#define DSCAL  SSCAL
#define DAXPY  SAXPY
#define DSYEV  SSYEV
#define DSYGV  SSYGV
#define DTRMM  STRMM
#define DTRTRI STRTRI
#define DPOTRF SPOTRF
#define DGETRF SGETRF
#define DGETRS SGETRS
#define PDPOTRF PSPOTRF
#define PDTRTRI PSTRTRI
#define PDSYEVX PSSYEVX
#define DCOPY   SCOPY
#define DSYEVX  SSYEVX

#endif
#ifdef NGZhalf
!
!   charge density: half grid mode Z direction
!
#define realmode
#define  zm(n)             (n/2+1)
#define  xm(n)              n
#define  xmin(n)          (-(n/2-1))
#define  SETN1fromN1P     N1=MOD(N1P+xmin(GRIDC%NGX)+GRIDC%NGX,GRIDC%NGX)
#define  SETFACT1         FACTM=1
#define  SETFACT          IF (N3 /= 1) FACTM=2
#define  SETSQRT          IF (N3 /= 1) FACTM=SQRT(2._q)
#define  SETISQRT         IF (N3 /= 1) FACTM=1./SQRT(2._q)
#define  MULFACT          FACTM*
#elif defined(NGXhalf)
!
!   charge density: half grid mode X direction
!
#define realmode
#define  xm(n)             (n/2+1)
#define  xmin(n)            0
#define  zm(n)              n
#define  SETN1fromN1P     N1=N1P
#define  SETFACT1         FACTM=1
#define  SETFACT          IF (N1 /= 1) FACTM=2
#define  SETSQRT          IF (N1 /= 1) FACTM=SQRT(2._q)
#define  SETISQRT         IF (N1 /= 1) FACTM=1./SQRT(2._q)
#define  MULFACT          FACTM*
#else
!
!   charge density: full grid mode
!
#define  zm(n)              n
#define  xm(n)              n
#define  xmin(n)          (-(n/2-1))
#define  SETN1fromN1P      N1=MOD(N1P+xmin(GRIDC%NGX)+GRIDC%NGX,GRIDC%NGX)
#define  SETFACT1          FACTM=1
#define  SETFACT
#define  SETSQRT
#define  SETISQRT
#define  MULFACT
#endif
#ifdef realmode
!
!   charge density real
!
#define RGRID             REAL(q)
#define OVERLAP           REAL(q)
#define DIMREAL(n)        n*2
#define CHGCONJG(n)       (n)
#else
!
!   charge density complex
!
#define RGRID            COMPLEX(q)
#define OVERLAP          COMPLEX(q)
#define DIMREAL(n)       n
#define CHGCONJG(n)      CONJG(n)
#endif

#ifdef wNGZhalf
!
!   wavefunctions: half grid mode for Z direction
!
#define gammareal
#elif defined(wNGXhalf)
!
!   wavefunctions: half grid mode for X direction
!
#define gammareal
#else
!
!   wavefunctions: full grid mode
!
#endif
#ifdef  gammareal
!
!   wavefunctions real (gamma only)
!
#define GREAL(n)  REAL(n,KIND=q)
#define GCMPLX(n,m,q)  (n)
#define GDEF      REAL(q)
#define GDEFS     REAL(qs)
#define GDEFN     REAL(qn)
#define M_sum_g   M_sum_d
#define M_sumf_g  M_sumf_d
#define M_sum_g_giant M_sum_d_giant
#define M_bcast_g M_bcast_d
#define M_bcast_g_from M_bcast_d_from
#define M_allgatherv_g M_allgatherv_d
#define CPROJ     GPROJ
#define GCONJG(n) (n)
#define GAIMAG(n) 0._q
#define GGEMM     DGEMM
#define GGEMV     DGEMV
#define GGETRF    DGETRF
#define GGETRS    DGETRS
#define GDOTC     DDOT
#define GDSCAL    DSCAL
#define GAXPY     DAXPY
#define GDOT      DDOT
#define PGGEMM    PDGEMM
#define PGGEADD   PDGEADD
#define zero      0._q
#define one       1._q
#define m_        2*
#define ndata     1
#define ndata8    1_qi8
#define trans     'T'
#define gammaarg   NPLWKC,NINDPC,WVSCAL,WVSCAI
#define gammaargNK NPLWKC(NK),NINDPC(1,NK),WVSCAL(1,NK),WVSCAI(1,NK)
#else
!
!   wavefunctions complex
!
#define GREAL(n)  (n)
#define GCMPLX(n,m,q)  CMPLX(n,m,KIND=q)
#define GDEF      COMPLEX(q)
#define GDEFS     COMPLEX(qs)
#define GDEFN     COMPLEX(qn)
#define M_sumf_g  M_sumf_z
#define M_sum_g   M_sum_z
#define M_sum_g_giant M_sum_z_giant
#define M_bcast_g M_bcast_z
#define M_bcast_g_from M_bcast_z_from
#define M_allgatherv_g M_allgatherv_z
#define GCONJG(n) CONJG(n)
#define GAIMAG(n) AIMAG(n)
#define GGEMM     ZGEMM
#define GGEMV     ZGEMV
#define GGETRF    ZGETRF
#define GGETRS    ZGETRS
#define GDOTC     ZDOTC
#define GDSCAL    ZDSCAL
#define GAXPY     ZAXPY
#define GDOT      ZDOTC
#define PGGEMM    PZGEMM
#define PGGEADD   PZGEADD
#define zero      (0._q,0._q)
#define one       (1._q,0._q)
#define trans     'C'
#define m_
#define ndata     2
#define ndata8    2_qi8
#define gammaarg
#define gammaargNK
#endif
!
!   common definitions
!
#define ztrans     'C'
#define dtrans     'T'
#define zzero      (0._q,0._q)
#define zone       (1._q,0._q)
#define dzero      0._q
#define done       1._q

#if defined(MPI)
!
!   mpi parallel macros
!
#define CALLMPI(x) CALL x
#define CALLMPI_C(x) CALL x
#define STOP       CALL M_stop('VASP aborting ...'); stop

#define io_begin IF (NODE_ME==IONODE) THEN
#define do_io    IF (NODE_ME==IONODE)
#define io_end   ENDIF
#define FFTGRIDPLAN FFTGRIDPLAN_MPI
#define FFTMAKEPLAN FFTMAKEPLAN_MPI
#define FFT3D       FFT3D_MPI
#define FFTCHK      FFTCHK_MPI
#define FFT3RC      FFT3RC_MPI
#define FFTEXT      FFTEXT_MPI
#define FFTWAV      FFTWAV_MPI
#define FFTINI      FFTINI_MPI
#elif defined(MPI_CHAIN)
!
!   mpi nudged chain version
!
#define CALLMPI(x)
#define CALLMPI_C(x) CALL x
#define STOP       CALL M_stop('VASP aborting ...'); stop

#define io_begin IF (NODE_ME==IONODE) THEN
#define do_io    IF (NODE_ME==IONODE)
#define io_end   ENDIF
#else
#define CALLMPI(x)
#define CALLMPI_C(x)
#define nodes
#define do_io
#define io_begin
#define io_end
#undef use_shmem
#undef shmem_rproj
#undef shmem_bcast_buffer
#endif
#ifdef F90_T3D
!
!   T3D macros
!
!#define IVDEP
#endif
#if defined(T3D_SMA)
!
!   T3D SHMEM communication
!
#endif
#ifdef CRAY_MPP
!
!   CRAY_MPP macros
!
#define MPI_double_precision MPI_real
#define MPI_double_complex   MPI_complex
#endif
#ifdef USE_ERF
!
!  use instrinsic ERF
!
#define ERRF(x) ERF(x)
#define ERRFC(x) ERFC(x)
#endif

#define dwrite DWRITE
#ifdef debug
!
!  debugging primitives
!
#define DWRITE WRITE(0,*)
#define DWRITE0 do_io WRITE(0,*)
#else
#define DWRITE  !
#define DWRITE0 !
#endif

#define CW CPTWFP

#ifdef nonlr_single
!
! single precision non local projection operators  
!
#define RPROJ   RPROJS
#define qn      qs
#define DGEMVn  SGEMV
#define DGEMMn  SGEMM
#else
#define qn      q
#define DGEMVn  DGEMV
#define DGEMMn  DGEMM
#endif


#if  ! (defined(NGXhalf)  || defined(NGYhalf) ||  defined(NGZhalf))
#define RHOLM_complex
#endif

#ifdef dotiming
!
!  some timing primitives
!
#define usetmg REAL(q) :: TV,TC,TV0,TC0
#define statmg CALLMPI( MPI_barrier( WDES%COMM%MPI_COMM, ierror )) ;  CALL VTIME(TV0,TC0)
#define stotmg CALLMPI( MPI_barrier( WDES%COMM%MPI_COMM, ierror )) ;  CALL VTIME(TV,TC) ; IF (NODE_ME==IONODE) WRITE(*,'(A,2F10.4)') 'timing ',TV-TV0,TC-TC0 ; CALL VTIME(TV0,TC0)
#define stotmgt(X) CALLMPI( MPI_barrier( WDES%COMM%MPI_COMM, ierror )) ;  CALL VTIME(TV,TC) ; IF (NODE_ME==IONODE) WRITE(*,'(A,2F10.4)') X,TV-TV0,TC-TC0 ; CALL VTIME(TV0,TC0)
#define addtmg(T1,T2) CALLMPI(MPI_barrier(WDES%COMM%MPI_COMM,ierror)); CALL VTIME(TV,TC); T1=T1+TV-TV0;T2=T2+TC-TC0; TV0=TV;TC0=TC
#define addtmga(T1,T2) CALL VTIME(TV,TC); T1=T1+TV-TV0;T2=T2+TC-TC0; TV0=TV;TC0=TC
#define wrttmg(X,T1,T2) IF (NODE_ME==IONODE) WRITE(*,'(A,2F10.4)') X,T1,T2
#else
#define usetmg 
#define statmg 
#define stotmg 
#define stotmgt(X)
#define addtmg(T1,T2)
#define addtmga(T1,T2)
#define wrttmg(X,T1,T2)
#endif
!
! OpenMP macros
!
#ifdef _OPENMP45
#define _OPENMP40
#define _OPENMP
#endif
#ifdef _OPENMP40
#define _OPENMP
#endif
#ifdef _OPENMP
#define __omp_id               omp_id
#define __omp_get_max_threads  ,OMP_GET_MAX_THREADS()
#define __omp_get_num_threads  ,OMP_GET_NUM_THREADS()
#define __omp_get_thread_num   ,OMP_GET_THREAD_NUM()+1
#define __omp_arg(x)           ,x
#define __omp_arg_id           __omp_arg(__omp_id)
#define __omp_add(x)           x
#else
#define __omp_id               1
#define __omp_get_max_threads
#define __omp_get_num_threads
#define __omp_get_thread_num
#define __omp_arg(x)
#define __omp_arg_id
#define __omp_add(x)
#endif
!
! profiling macros
!
#ifdef PROFILING
#define PROFILING_START(x)     CALL START_PROFILING(x)
#define PROFILING_STOP(x)      CALL STOP_PROFILING(x)
#else
#define PROFILING_START(x)
#define PROFILING_STOP(x)
#endif
!
! shmem macros
!
#ifdef shmem_rproj
#define use_shmem
#endif
#ifdef shmem_bcast_buffer
#define use_shmem
#endif
!
! quadruple precision 
!
#ifdef qd_emulate
#define QDPREAL                TYPE(qd_real)
#define QDPCOMPLEX             TYPE(qd_complex)
#define QDPD(x)                qdreal(REAL(x,KIND=q))
#define TOPRINT(x)             TOREAL(x)
#define QDPC(x,y )             to_qdc_qd2( x,y ) 
#define QDPC2QDPREAL( x )      to_qd_qdc( x ) 
#ifdef __GNUC__
#define QDPP(x)                qdreal('x')
#define QDFP(x)                qdreal('x')
#else
#define QDPP(x)                qdreal(#x)
#define QDFP(x)                qdreal(#x)
#endif
#else
#define TOPRINT(x)             x
#define QDPC2QDPREAL( x )      REAL( x, qd ) 
#define QDPREAL                REAL(qd)
#define QDPCOMPLEX             COMPLEX(qd)
#define QDPD(x)                REAL(x,KIND=qd)
#define QDPP(x)                REAL(x,KIND=qd)
#define QDPC(x,y)              CMPLX(x,y,KIND=qd)
#define IDENTITYMACRO(x)       x 
#define QDFP(x)                IDENTITYMACRO(x)_qd
#endif
!
! for the VASP2WANNIER90 interface
!
#ifdef VASP2WANNIER90v2
#define VASP2WANNIER90
#endif
!
! CUDA includes
!
#ifdef CUDA_GPU
#include "cuda_fft.inc"
#include "cuda_macros.inc"
#include "cuda_profiling.inc"
#endif
!
! SIMD related definitions
!
#include "simd.inc"
!
! memalign macros
!
#if   defined(memalign64)
#define memalign
#define __align_BTS            64
#define __align_DBL            8
#define __align_DIM(x)         ((x+__align_DBL-1)/__align_DBL)*__align_DBL
#elif defined(memalign32)
#define memalign
#define __align_BTS            32
#define __align_DBL            4
#define __align_DIM(x)         ((x+__align_DBL-1)/__align_DBL)*__align_DBL
#else
#define __align_DIM(x)         x
#endif

!
! Macros for PGI/NV HPC compilers version specific code
!
#ifdef __PGIC__
#define NV_MAJOR_VER __PGIC__
#define NV_MINOR_VER __PGIC_MINOR__
#else
#ifdef __PGIF90__
#define NV_MAJOR_VER __PGIF90__
#define NV_MINOR_VER __PGIF90_MINOR__
#elif __NVCOMPILER
#define NV_MAJOR_VER __NVCOMPILER_MAJOR__
#define NV_MINOR_VER __NVCOMPILER_MINOR__  
#endif
#endif

#ifdef NV_MAJOR_VER
#define PGI_BEFORE_20_11 (NV_MAJOR_VER < 20 || (NV_MAJOR_VER == 20 && NV_MINOR_VER < 11))
#define PGI_BEFORE_21_1 (NV_MAJOR_VER < 21 || (NV_MAJOR_VER == 21 && NV_MINOR_VER < 1))
#define PGI_BEFORE_21_2 (NV_MAJOR_VER < 21 || (NV_MAJOR_VER == 21 && NV_MINOR_VER < 2))
#define PGI_BEFORE_XX_X 1
#define PGI_BETWEEN_20_9_AND_21_3 ((NV_MAJOR_VER == 20 && NV_MINOR_VER > 9) || (NV_MAJOR_VER == 21 && NV_MINOR_VER < 3))
#define PGI_BETWEEN_20_9_AND_XX_X ((NV_MAJOR_VER == 20 && NV_MINOR_VER > 9) || NV_MAJOR_VER > 20)
#else
#define PGI_BEFORE_20_11 0
#define PGI_BEFORE_21_1 0
#define PGI_BEFORE_21_2 0
#define PGI_BEFORE_XX_X 0
#define PGI_BETWEEN_20_9_AND_21_3 0
#define PGI_BETWEEN_20_9_AND_XX_X 0
#endif

!
! OpenACC macros
!
#ifdef _OPENACC
#define _OPENACC_TST
#define NOACC                  !!
#define DOACC
#define PUSH_ACC_EXEC_ON(x)    CALL ACC_PUSH_EXEC_ON(x)
#define POP_ACC_EXEC_ON        CALL ACC_POP_EXEC_ON
#if PGI_BEFORE_21_1
#define SET_ACC_ASYNC_Q(x)     CALL ACC_SET_ASYNC_Q(ACC_ASYNC_ASYNC)
#else
#define SET_ACC_ASYNC_Q(x)     CALL ACC_SET_ASYNC_Q(x)
#endif
#define __IF_ACC_OFF__(x)      (x.AND.(.NOT.ACC_EXEC_ON))
#define __IF_ASYNC__           IF(ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
#define __DCOPY__              ACC_DCOPY
#define __ZCOPY__              ACC_ZCOPY
#define __DAXPY__              ACC_DAXPY
#define __ZAXPY__              ACC_ZAXPY
#define __DSCAL__              ACC_DSCAL
#define __ZSCAL__              ACC_ZSCAL
#define __DDOT__               ACC_DDOT
#define __ZDOTC__              ACC_ZDOTC
#define __DGEMV__              ACC_DGEMV
#define __ZGEMV__              ACC_ZGEMV
#define __DGEMM__              ACC_DGEMM
#define __ZGEMM__              ACC_ZGEMM
#define __DTRMM__              ACC_DTRMM
#define __ZTRMM__              ACC_ZTRMM
#define __DPOTRF__             ACC_DPOTRF
#define __ZPOTRF__             ACC_ZPOTRF
#define __DGETRF__             ACC_DGETRF
#define __ZGETRF__             ACC_ZGETRF
#define __DGETRS__             ACC_DGETRS
#define __ZGETRS__             ACC_ZGETRS
#define __DTRTRI__             ACC_DTRTRI
#define __ZTRTRI__             ACC_ZTRTRI
#define __DSYGV__              ACC_DSYGV
#define __DSYEV__              ACC_DSYEV
#define __DSYEVX__             ACC_DSYEVX
#define __ZHEGV__              ACC_ZHEGV
#define __ZHEEV__              ACC_ZHEEV
#define __ZHEEVX__             ACC_ZHEEVX
#define __ZDSCAL__             ACC_ZDSCAL
#ifdef gammareal
#define __GCOPY__              ACC_DCOPY
#define __GGEMV__              ACC_DGEMV
#define __GGEMM__              ACC_DGEMM
#define __GDOTC__              ACC_DDOT
#define __GDSCAL__             ACC_DSCAL
#else
#define __GCOPY__              ACC_ZCOPY
#define __GGEMV__              ACC_ZGEMV
#define __GGEMM__              ACC_ZGEMM
#define __GDOTC__              ACC_ZDOTC
#define __GDSCAL__             ACC_ZDSCAL
#endif
#define __TUTOR__(x, y, u, v)  CALL ACC_VTUTOR_LOG(x, y)
#else
#define NOACC
#define DOACC                  !!
#define PUSH_ACC_EXEC_ON(x)
#define POP_ACC_EXEC_ON
#define SET_ACC_ASYNC_Q(x)
#define __IF_ACC_OFF__(x)      x
#define __IF_ASYNC__
#define __DCOPY__              DCOPY
#define __ZCOPY__              ZCOPY
#define __DAXPY__              DAXPY
#define __ZAXPY__              ZAXPY
#define __DSCAL__              DSCAL
#define __ZSCAL__              ZSCAL
#define __DDOT__               DDOT
#define __ZDOTC__              ZDOTC
#define __DGEMV__              DGEMV
#define __ZGEMV__              ZGEMV
#define __DGEMM__              DGEMM
#define __ZGEMM__              ZGEMM
#define __DTRMM__              DTRMM
#define __ZTRMM__              ZTRMM
#define __DPOTRF__             DPOTRF
#define __ZPOTRF__             ZPOTRF
#define __DGETRF__             DGETRF
#define __ZGETRF__             ZGETRF
#define __DGETRS__             DGETRS
#define __ZGETRS__             ZGETRS
#define __DTRTRI__             DTRTRI
#define __ZTRTRI__             ZTRTRI
#define __DSYGV__              DSYGV
#define __DSYEV__              DSYEV
#define __DSYEVX__             DSYEVX
#define __ZHEGV__              ZHEGV
#define __ZHEEV__              ZHEEV
#define __ZHEEVX__             ZHEEVX
#define __ZDSCAL__             ZDSCAL
#ifdef gammareal
#define __GCOPY__              DCOPY
#define __GGEMV__              DGEMV
#define __GGEMM__              DGEMM
#define __GDOTC__              DDOT
#define __GDSCAL__             DSCAL
#else
#define __GCOPY__              ZCOPY
#define __GGEMV__              ZGEMV
#define __GGEMM__              ZGEMM
#define __GDOTC__              ZDOTC
#define __GDSCAL__             ZDSCAL
#endif
#define __TUTOR__(x, y, u, v)  CALL vtutor%write(x, y, filename = u, linenumber = v)
#endif
!
! combined OpenMP and OpenACC macros
!
#if defined(_OPENMP) && ! defined(_OPENACC)
#define __NOACC_omp_arg(x)     ,x
#else
#define __NOACC_omp_arg(x)
#endif
!
! routines replaced in LAPACK >=3.6 
!
#ifdef LAPACK36
#define DGEGV DGGEV
#endif

!
! Macros for HDF5 error check
!
#ifdef VASP_HDF5
#define VH5_CHECK(e) call vh5_error(e,__FILE__,__LINE__)
#else
#define VH5_CHECK(e)
#endif

!
! Macros for GNU version specific code
!
#ifdef __GNUC__
#define GNU_BEFORE_7_5 (__GNUC__ < 7 || (__GNUC__ == 7 && __GNUC_MINOR__ < 5))
#endif

!
! For machine learning
!
#if defined(MPI) || defined(MPI_CHAIN)
#define ML_AVAILABLE
#endif

!
! Extra safe initializations (+ overflow protections)
!
#ifdef extra_safe_init
#define DOESI
#else
#define DOESI !
#endif
