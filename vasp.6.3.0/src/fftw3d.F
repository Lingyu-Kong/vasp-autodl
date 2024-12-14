#define use_r2c_fftw
#define inversion_through_mid_point
#include "fft3dbatched.F"
#include "fft3dsimple.F"
#ifndef MINLOOP
#define MINLOOP 1
#endif
#ifndef CACHE_SIZE

#if defined vector
#define CACHE_SIZE 0
#else
#define CACHE_SIZE 32768
#endif

#endif
!*************************************************************************
!*        Fast Fourier Transform for VASP using the FFTW package.        *
!*                                                                       *
!* This routine is just to give a good interface between VASP and FFTW.  *
!* FFTW3D gives an identical interface to the FFTBAS routine, apart from *
!* that the complex array is placed last instead of first. This is just  *
!* to prevent alignment errors. I.e. if the c-array would be too large   *
!* it is simply cut off with this arrangement (the other way around might*
!* cause trouble on the Macintosh.)                                      *
!* The Fortran calls simply makes use of the FFTW library written in     *
!* C/C++ (see also http://www.fftw.org). The result is considerably      *
!* faster then the fft3dfurth/fftbas system for all array sizes I have   *
!* tested. For small arrays (n<32) you get an improvement of several     *
!* 100%, and for larger you still get an improvement in the range of 50%.*
!* Those who like to use the FFTW should get the package from the url    *
!* above, compile it and link with the library. For the Mac the Absoft   *
!* C/C++ compiler makes a good job even if you have to create your own   *
!* your own makefile.                                                    *
!* It should also be noted that FFTW is just as conservative about mem-  *
!* ory use as FFTBAS, possibly even more so.                             *
!*                                                                       *
!*                                   Ph.D. Stefan Mankefors Mars 3 2000  *
!*                                         sem@fy.chalmers.se            *
!*                                                                       *
!* NOTE: When you compile the library you have to configure the Fortran  *
!* conventions used, i.e. 'upper case', 'lowercase' etc. This is done at *
!* the end of the 'config.h' file where you define *one* of the the conv-*
!* ventions used for Fortran. ('config.h' comes with all F77 conventions *
!* undefined.) Please note that this might be machine dependent.         *
!*                                                                       *
!* NOTE2: The real to complex FFT might also be exchanged in a similar   *
!* way. I have not tested this, however, and the gain in time is very    *
!* slim since this type of FFT is used very little by the VASP - at least*
!* as far as the benchmark runs goes (I do not have the experience yet to*
!* claim anything else.) Hence it is a question of an additional library *
!* against a gain of perhaps 1%. I have choosen not to use FFTW for this.*
!* Please observe that this means that fft3dlib still is needed in the   *
!* make process of VASP.                                                 *
!*************************************************************************

!> @details @ref openmp :
!> under OpenMP this subroutine calls dfftw_plan_with_nthreads.
!> Note: this subroutine may be called from within an OpenMP parallel
!> region as well: threadsafety is ensured by an explicit CRITICAL region
!> (VASP_FFT_PLAN_CREATE_DESTROY).
!
      SUBROUTINE FFTBAS(C,GRID,ISIGN)
#ifdef _OPENACC
      USE mopenacc
      USE mcufft
#endif
      USE prec
      USE mgrid

      USE iso_c_binding

      IMPLICIT NONE

#ifndef __NEC__
      include 'fftw3.f'
#else // __NEC__
      include 'aslfftw3.f'
#endif // __NEC__

      TYPE (grid_3d) GRID

      COMPLEX(q) :: C(*)
      INTEGER    :: ISIGN

      ! local variables
      TYPE(c_ptr) PLAN
      INTEGER(c_int) EFFORT

#if defined(fftw_cache_plans) || defined(__NEC_TUNE__)
      TYPE(c_ptr), SAVE :: PLAN_FWD = c_null_ptr, PLAN_BWD = c_null_ptr
      INTEGER, SAVE :: NFFT_FWD, NFFT_BWD
#endif

!$    INTEGER THREADS
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS

#ifdef _OPENACC
      TYPE(acc_cufft_plan) :: ACC_PLAN
      INTEGER :: IERR

      IF (ACC_IS_PRESENT(C).AND.ACC_EXEC_ON) THEN
         IERR=ACC_CUFFT_GETPLAN(GRID,INT(ACC_ASYNC_Q),1,ACC_PLAN)
         IF (IERR/=0) THEN
            CALL vtutor%error("FFTBAS: no cuFFT plan found for ACC_ASYNC_Q= " // str(INT(ACC_ASYNC_Q)))
         ENDIF
!$ACC HOST_DATA USE_DEVICE(C)
         IF (ISIGN.LE.0) THEN
            IERR=cufftExecZ2Z(ACC_PLAN%CUFFT_3D_Z2Z,C,C,CUFFT_FORWARD)
         ELSE
            IERR=cufftExecZ2Z(ACC_PLAN%CUFFT_3D_Z2Z,C,C,CUFFT_INVERSE)
         ENDIF
!$ACC END HOST_DATA
         RETURN
      ENDIF
#endif
      EFFORT=FFTW_ESTIMATE
#ifdef use_fftw_plan_effort
      IF (FFTW_PLAN_EFFORT==1) EFFORT=FFTW_MEASURE
#endif
!$!=======================================================================
!$!  initialise openMP FFT, has to be done here and not in main.F in
!$!  in general, since the FFTs are called outside and inside openMP
!$!  parallel regions.
!$!=======================================================================
!$OMP PARALLEL SHARED(THREADS)
!$OMP MASTER
!$    THREADS=OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!$    CALL dfftw_plan_with_nthreads(THREADS) 
#if !defined(fftw_cache_plans) && !defined(__NEC_TUNE__)
      IF (ISIGN.LE.0) THEN
         CALL dfftw_plan_dft_3d(PLAN,GRID%NGX,GRID%NGY,GRID%NGZ,C,C,FFTW_FORWARD, EFFORT)
      ELSE
         CALL dfftw_plan_dft_3d(PLAN,GRID%NGX,GRID%NGY,GRID%NGZ,C,C,FFTW_BACKWARD,EFFORT)
      ENDIF
#else // fftw_cache_plans or __NEC_TUNE__
      IF (ISIGN.LE.0) THEN
!!         write(*,*) 'forward:',plan_fwd,nfft_fwd,grid%nplwv,grid%ngx,grid%ngy,grid%ngz
         IF (c_associated(PLAN_FWD) .AND. NFFT_FWD == GRID%NPLWV) THEN
!!            write(*,*) 'using plan_fwd:',plan_fwd
            PLAN=PLAN_FWD
         ELSE
            IF (c_associated(PLAN_FWD)) CALL dfftw_destroy_plan(PLAN_FWD)
            CALL dfftw_plan_dft_3d(PLAN,GRID%NGX,GRID%NGY,GRID%NGZ,C,C,FFTW_FORWARD,EFFORT)
            PLAN_FWD=PLAN
            NFFT_FWD=GRID%NPLWV
!!            write(*,*) 'generated plan_fwd:',plan_fwd,nfft_fwd
         ENDIF
      ELSE
!!         write(*,*) 'bckward:',plan_bwd,nfft_bwd,grid%nplwv,grid%ngx,grid%ngy,grid%ngz
         IF (c_associated(PLAN_BWD) .AND. NFFT_BWD == GRID%NPLWV) THEN
!!            write(*,*) 'using plan_bck:',plan_bwd
            PLAN=PLAN_BWD
         ELSE
            IF (c_associated(PLAN_BWD)) CALL dfftw_destroy_plan(PLAN_BWD)
            CALL dfftw_plan_dft_3d(PLAN,GRID%NGX,GRID%NGY,GRID%NGZ,C,C,FFTW_BACKWARD,EFFORT)
            PLAN_BWD=PLAN
            NFFT_BWD=GRID%NPLWV
!!            write(*,*) 'generated plan_bck:',plan_bwd,nfft_bwd
         ENDIF
      ENDIF
#endif // fftw_cache_plans or __NEC_TUNE__
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

      CALL dfftw_execute_dft(PLAN,C,C)

#if !defined(fftw_cache_plans) && !defined(__NEC_TUNE__)
!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
      CALL dfftw_destroy_plan(PLAN)
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
#endif // fftw_cache_plans or __NEC_TUNE__

      RETURN
      END SUBROUTINE


      SUBROUTINE FFTGRIDPLAN(GRID)
      USE prec
      USE mgrid

      IMPLICIT NONE

      TYPE (grid_3d) GRID

      ! local variables
      COMPLEX(q), ALLOCATABLE :: CWORK(:)
      INTEGER, PARAMETER :: NADD=1024, NTIMES=4
      INTEGER :: I

      ALLOCATE(CWORK(GRID%MPLWV+NADD))

      DO I=1,NTIMES
         CALL INIDAT(GRID%MPLWV+NADD,CWORK)
         CALL FFTMAKEPLAN(CWORK(I),GRID)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP MASTER
!$       CALL FFTMAKEPLAN_MPI(CWORK(I),GRID)
!$OMP END MASTER
!$OMP END PARALLEL
      ENDDO

      DEALLOCATE(CWORK)

      RETURN
      END SUBROUTINE FFTGRIDPLAN


!> @details @ref openmp :
!> if openmp::omp_dfftw_init_threads=.TRUE. this subroutine calls
!> dfftw_init_threads, and sets openmp::omp_dfftw_init_threads=.FALSE. to
!> make sure dfftw_init_threads is called only once.
!
      SUBROUTINE FFTMAKEPLAN(c,grid)

      USE prec
      USE mgrid
      USE openmp, ONLY : omp_dfftw_init_threads

      USE iso_c_binding

      IMPLICIT NONE

#ifndef __NEC__
      include 'fftw3.f'
#else // __NEC__
      include 'aslfftw3.f'
#endif // __NEC__

      TYPE (grid_3d) grid
      COMPLEX(q) c(*), cdummy

      ! local variables
      INTEGER nx,ny,nz
      TYPE(c_ptr) plan
      INTEGER(c_int) EFFORT

!$    INTEGER IERROR, THREADS
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS

      EFFORT=FFTW_MEASURE
      IF (FFTW_PLAN_EFFORT==0) EFFORT=FFTW_ESTIMATE
      IF (FFTW_PLAN_EFFORT==1) EFFORT=FFTW_MEASURE

!$    IF (omp_dfftw_init_threads) THEN
!$       CALL dfftw_init_threads(ierror)
!$       omp_dfftw_init_threads=.FALSE.
!$    ENDIF

      nx=grid%ngx
      ny=grid%ngy
      nz=grid%ngz

!$OMP PARALLEL SHARED(THREADS)
!$OMP MASTER
!$    THREADS=OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!$    CALL dfftw_plan_with_nthreads(THREADS)
      IF (.NOT.grid%real2cplx) THEN
         CALL dfftw_plan_dft_3d(plan,nx,ny,nz, &
              c, c, & 
              FFTW_FORWARD, EFFORT)
         CALL dfftw_destroy_plan(plan)
         CALL dfftw_plan_dft_3d(plan,nx,ny,nz,&
              c, c, & 
              FFTW_BACKWARD, EFFORT)
         CALL dfftw_destroy_plan(plan)
      ELSE
         CALL dfftw_plan_dft_r2c_3d(plan,nx,ny,nz,&
              c,c, &
              EFFORT)
         CALL dfftw_destroy_plan(plan)
         CALL dfftw_plan_dft_c2r_3d(plan,nx,ny,nz,&
              c,c, &
              EFFORT)
         CALL dfftw_destroy_plan(plan)
      ENDIF
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

      RETURN
      END SUBROUTINE FFTMAKEPLAN


!> @details @ref openmp :
!> under OpenMP this subroutine calls dfftw_plan_with_nthreads.
!> Note: this subroutine may be called from within an OpenMP parallel
!> region as well: threadsafety is ensured by an explicit CRITICAL region
!> (VASP_FFT_PLAN_CREATE_DESTROY).
!
#ifdef use_r2c_fftw
      SUBROUTINE FFTBRC(C,grid,ISIGN)
#ifdef _OPENACC
      USE mopenacc
      USE mcufft
#endif
      USE prec
      USE mgrid

      USE iso_c_binding

      IMPLICIT NONE

#ifndef __NEC__
      include 'fftw3.f'
#else // __NEC__
      include 'aslfftw3.f'
#endif // __NEC__
      TYPE (grid_3d) grid
      COMPLEX(q) C(*)
      INTEGER ISIGN
      ! local variables
      TYPE(c_ptr) plan
      INTEGER(c_int) EFFORT
!$    INTEGER THREADS, IERROR
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS
#ifdef _OPENACC
      TYPE(acc_cufft_plan) :: ACC_PLAN
      REAL(q), POINTER  :: A(:)
      INTEGER :: IERR

      IF (ACC_IS_PRESENT(C).AND.ACC_EXEC_ON) THEN
         IERR=ACC_CUFFT_GETPLAN(GRID,INT(ACC_ASYNC_Q),1,ACC_PLAN)
         CALL c_f_pointer(c_loc(C),A,[2*GRID%MPLWV])
!$ACC HOST_DATA USE_DEVICE(C,A)
         IF (ISIGN.LE.0) THEN
            IERR=cufftExecD2Z(ACC_PLAN%CUFFT_3D_D2Z,A,C)
         ELSE
            IERR=cufftExecZ2D(ACC_PLAN%CUFFT_3D_Z2D,C,A)
         ENDIF
!$ACC END HOST_DATA
         RETURN
      ENDIF
#endif
      EFFORT=FFTW_ESTIMATE
#ifdef use_fftw_plan_effort
      IF (FFTW_PLAN_EFFORT==1) EFFORT=FFTW_MEASURE
#endif
!$!=======================================================================
!$!  initialise openMP FFT, has to be done here and not in main.F in
!$!  in general, since the FFTs are called outside and inside openMP
!$!  parallel regions.
!$!=======================================================================
!$OMP PARALLEL SHARED(THREADS)
!$OMP MASTER
!$    THREADS=OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL
!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!$    CALL dfftw_plan_with_nthreads(THREADS) 
      IF (ISIGN.LE.0) THEN
         CALL dfftw_plan_dft_r2c_3d(plan,grid%ngx,grid%ngy,grid%ngz,C,C,EFFORT)
      ELSE
         CALL dfftw_plan_dft_c2r_3d(plan,grid%ngx,grid%ngy,grid%ngz,C,C,EFFORT)
      ENDIF
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
      IF (ISIGN.LE.0) THEN
         CALL dfftw_execute_dft_r2c(plan,C,C)
      ELSE
         CALL dfftw_execute_dft_c2r(plan,C,C)
      ENDIF
!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
      CALL dfftw_destroy_plan(plan)
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
      RETURN
      END SUBROUTINE FFTBRC
#endif

!=======================================================================
!   this routine returns the next correct setting for the
!   three dimensional FFT
!=======================================================================

      SUBROUTINE FFTCHK(NFFT)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION NFFT(3)
      LOGICAL FFTCH1

      DO 100 IND=1,3
  200 CONTINUE
        IF (FFTCH1(NFFT(IND))) GOTO 100
        NFFT(IND)=NFFT(IND)+1
        GOTO 200
  100 CONTINUE
      END

      LOGICAL FUNCTION FFTCH1(NIN)
      USE prec
      LOGICAL, EXTERNAL :: FFTCHK_FURTH

      FFTCH1=FFTCHK_FURTH(NIN)
      END FUNCTION

      MODULE fft_private
      USE prec
      REAL(q),POINTER,SAVE ::  WORK(:)
      END MODULE

#ifndef use_r2c_fftw
!=======================================================================
!   generic   3-d fast fourier transformation
!   written by Juergen Furthmueller
!   performes the 3-d real to complex FFT
!
!     +1  q->r   vr= sum(q) vq exp(+iqr)    complex to real
!     -1  r->q   vq= sum(r) vr exp(-iqr)    real to complex
!
!=======================================================================
      SUBROUTINE FFTBRC(A,grid,ISIGN)
      USE prec
      USE mgrid
      USE fft_private
      USE smart_allocate
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE (grid_3d) grid

      COMPLEX(q) A(*)
      INTEGER, PARAMETER :: NMAXM=4096 ! propably sufficient forever
      DIMENSION TRIGX(2*NMAXM),TRIGY(2*NMAXM),TRIGZ(2*NMAXM)
      DIMENSION IFAC(19,3)
      SAVE TRIGX,TRIGY,TRIGZ,IFAC,NXO,NYO,NZO
      DATA NXO /0/, NYO /0/, NZO /0/

      NX=grid%ngx
      NY=grid%ngy
      NZ=grid%ngz
      NMAX=MAX(NX,NY)
      NMAX=MAX(NMAX,NZ)

#if CACHE_SIZE==0
#if defined vector
      NALLOC=4*NX*NY*NZ
#else
      NALLOC=4*NX*MAX(NY,NZ)
#endif
#else
      NALLOC=MAX(CACHE_SIZE, 4*MINLOOP*MAX(NX,NY,NZ))
#endif
      CALL SMART_ALLOCATE_REAL(WORK,NALLOC)

      IF (NX>NMAXM) THEN
         CALL vtutor%error("FFT3DFURTH: Increase NMAXM to " // str(NMAX))
      ENDIF
      IF (NY>NMAXM) THEN
         CALL vtutor%error("FFT3DFURTH: Increase NMAXM to " // str(NMAX))
      ENDIF
      IF (NZ>NMAXM) THEN
         CALL vtutor%error("FFT3DFURTH: Increase NMAXM to " // str(NMAX))
      ENDIF
! Initialize FFT if necessary (changes of mesh size, first call)!
      IF ((NX/=NXO).OR.(NY/=NYO).OR.(NZ/=NZO)) THEN
#if defined vector
      IF ((CACHE_SIZE==0).OR.(CACHE_SIZE>=4*NX*NY*NZ)) THEN
         CALL FFTR3V(A,NX+2,NX,NY,NY,NZ,NZ,WORK,IFAC, &
     &               TRIGX,TRIGY,TRIGZ,0,0,IERR)
      ELSE
         CALL FFTCRN(A,NX+2,NX,NY,NY,NZ,NZ,WORK,IFAC, &
     &               TRIGX,TRIGY,TRIGZ,0,0,IERR,CACHE_SIZE)
      ENDIF
#else
         CALL FFTCRN(A,NX+2,NX,NY,NY,NZ,NZ,WORK,IFAC, &
     &               TRIGX,TRIGY,TRIGZ,0,0,IERR,CACHE_SIZE)
#endif
         IF (IERR/=0) THEN
            CALL vtutor%error("INIT FFT3D: IERR = " // str(IERR))
         ENDIF
! Remember last mesh sizes!
         NXO=NX
         NYO=NY
         NZO=NZ
      END IF
! Do the transformation!
#if defined vector
      IF ((CACHE_SIZE==0).OR.(CACHE_SIZE>=4*NX*NY*NZ)) THEN
      CALL FFTR3V(A,NX+2,NX,NY,NY,NZ,NZ,WORK,IFAC, &
     &            TRIGX,TRIGY,TRIGZ,ISIGN,-ISIGN,IERR)
      ELSE
      CALL FFTCRN(A,NX+2,NX,NY,NY,NZ,NZ,WORK,IFAC, &
     &            TRIGX,TRIGY,TRIGZ,ISIGN,-ISIGN,IERR,CACHE_SIZE)
      ENDIF
#else
      CALL FFTCRN(A,NX+2,NX,NY,NY,NZ,NZ,WORK,IFAC, &
     &            TRIGX,TRIGY,TRIGZ,ISIGN,-ISIGN,IERR,CACHE_SIZE)
#endif
      IF (IERR/=0) THEN
         CALL vtutor%error("FFT3D: IERR = " // str(IERR))
      ENDIF

      RETURN
      END
#endif
