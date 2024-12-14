#include "symbol.inc"
#ifndef CACHE_SIZE
#define CACHE_SIZE 8000
#endif
#ifndef MINLOOP
#define MINLOOP 1
#endif

#undef FFTGRIDPLAN
#undef FFTMAKEPLAN
#undef FFT3D 
#undef FFT3DR 
#undef FFTCHK
#undef FFT3RC
#undef FFTEXT
#undef FFTWAV
#undef FFTINI


    SUBROUTINE FFTGRIDPLAN_MPI(GRID)
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
         CALL FFTMAKEPLAN_MPI(CWORK(I),GRID)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP MASTER
!$       CALL FFTMAKEPLAN_MPI(CWORK(I),GRID)
!$OMP END MASTER
!$OMP END PARALLEL
      ENDDO

      DEALLOCATE(CWORK)

      RETURN
    END SUBROUTINE FFTGRIDPLAN_MPI


!===============================================================================
!
!  FFTW requires a plan therefore this new calling interface is provided
!  which calls the FFTBAS and FFTBRC routine for generating these
!  plans
!
!> @details @ref openmp :
!> if openmp::omp_dfftw_init_threads=.TRUE. this subroutine calls
!> dfftw_init_threads, and sets openmp::omp_dfftw_init_threads=.FALSE. to
!> make sure dfftw_init_threads is called only once.
!
!===============================================================================

    SUBROUTINE FFTMAKEPLAN_MPI(A,GRID)
      USE prec
      USE mgrid
      USE openmp, ONLY : omp_dfftw_init_threads

      USE iso_c_binding

#ifndef __NEC__
      include 'fftw3.f'
#else // __NEC__
      include 'aslfftw3.f'
#endif // __NEC__

      TYPE (grid_3d) GRID
      REAL(q) A(*)

      ! local variables
      TYPE(c_ptr) :: planx, plany, planz
      INTEGER :: NX,NY,NZ,IDX,IDY,IDZ

      INTEGER(c_int) :: EFFORT

!$    INTEGER IERROR,THREADS
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS

!serFFT
      IF (GRID%RL%NFAST==1 .AND. GRID%RL_FFT%NFAST==1) THEN
         CALL FFTMAKEPLAN(A,GRID)
         RETURN
      ENDIF
!serFFTend

      EFFORT=FFTW_MEASURE
      IF (FFTW_PLAN_EFFORT==0) EFFORT=FFTW_ESTIMATE
      IF (FFTW_PLAN_EFFORT==1) EFFORT=FFTW_MEASURE

!$    IF (omp_dfftw_init_threads) THEN
!$       CALL dfftw_init_threads(ierror)
!$       omp_dfftw_init_threads=.FALSE.
!$    ENDIF

      NX=GRID%NGPTAR(1)
      NY=GRID%NGPTAR(2)
      NZ=GRID%NGPTAR(3)

      IDX=NX
      IDY=NY
      IDZ=NZ

!$OMP PARALLEL SHARED(THREADS)
!$OMP MASTER
!$    THREADS=OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!$    CALL dfftw_plan_with_nthreads(THREADS)
      IF (GRID%RC%NCOL > 0) THEN
         CALL dfftw_plan_many_dft(planx, 1, NX , GRID%RC%NCOL, &
                             A(1), NX, 1 , IDX, &
                             A(1), NX, 1 , IDX, &
                             FFTW_BACKWARD, EFFORT)
         CALL dfftw_destroy_plan(planx)
      ENDIF

      IF (GRID%IN%NCOL > 0) THEN
         CALL dfftw_plan_many_dft(plany, 1, NY , GRID%IN%NCOL, &
                             A(1), NY, GRID%IN%NCOL, 1 , &
                             A(1), NY, GRID%IN%NCOL, 1 , &
                             FFTW_BACKWARD, EFFORT)
         CALL dfftw_destroy_plan(plany)
      ENDIF

      IF (GRID%RL_FFT%NCOL > 0) THEN
         IF (NZ/2+1==GRID%NGZ_rd) THEN
           CALL dfftw_plan_many_dft_c2r(planz, 1, NZ , GRID%RL_FFT%NCOL, &
                             A(1), NZ, 1, (IDZ+2)/2 , &
                             A(1), NZ, 1, IDZ+2 , &
                             EFFORT)
         ELSE
           CALL dfftw_plan_many_dft(planz, 1, NZ , GRID%RL_FFT%NCOL, &
                             A(1), NZ, 1, IDZ , &
                             A(1), NZ, 1, IDZ , &
                             FFTW_BACKWARD, EFFORT)
         ENDIF
         CALL dfftw_destroy_plan(planz)
      ENDIF

      IF (GRID%RL_FFT%NCOL > 0) THEN
         IF (NZ/2+1==GRID%NGZ_rd) THEN
           CALL dfftw_plan_many_dft_r2c(planz, 1, NZ , GRID%RL_FFT%NCOL, &
                             A(1), NZ, 1, IDZ+2 , &
                             A(1), NZ, 1, (IDZ+2)/2 , &
                             EFFORT)
         ELSE
           CALL dfftw_plan_many_dft(planz, 1, NZ , GRID%RL_FFT%NCOL, &
                             A(1), NZ, 1, IDZ , &
                             A(1), NZ, 1, IDZ , &
                             FFTW_FORWARD, EFFORT)
         ENDIF
         CALL dfftw_destroy_plan(planz)
      ENDIF

      IF (GRID%IN%NCOL > 0) THEN
         CALL dfftw_plan_many_dft(plany, 1, NY , GRID%IN%NCOL, &
                             A(1), NY, GRID%IN%NCOL, 1 , &
                             A(1), NY, GRID%IN%NCOL, 1 , &
                             FFTW_FORWARD, EFFORT)
         CALL dfftw_destroy_plan(plany)
      ENDIF

      IF (GRID%RC%NCOL > 0) THEN
         CALL dfftw_plan_many_dft(planx, 1, NX , GRID%RC%NCOL, &
                             A(1), NX, 1 , IDX, &
                             A(1), NX, 1 , IDX, &
                             FFTW_FORWARD, EFFORT)
         CALL dfftw_destroy_plan(planx)
      ENDIF
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

      RETURN
    END SUBROUTINE


!-----------------------------------------------------------------------
! RCS:  $Id: fftmpi.F,v 1.3 2002/08/14 13:59:38 kresse Exp $
!
!   3-d parallel fast fourier transformation using fftw
!   written by Georg Kresse
!
!     +1  q->r   vr= sum(q) vq exp(+iqr)
!     -1  r->q   vq= sum(r) vr exp(-iqr)
! 
!   the FFTBAS_PLAN routine performs both the complex to complex, and
!   complex to real FFT
!
!   the FFTBAS routine is the calling interface for the 
!    complex, complex FFT 
!   whereas FFTBRC is the calling interface for complex to real FFT
!  
!=======================================================================

!
!  this subroutine calls the FFTBAS_PLAN routine  with FFTW_ESTIMATE
!

    SUBROUTINE FFTBAS_MPI(A,GRID,ISIGN)
      USE prec
      USE mgrid

      USE iso_c_binding

#ifndef __NEC__
      include 'fftw3.f'
#else // __NEC__
      include 'aslfftw3.f'
#endif // __NEC__

      TYPE (grid_3d) GRID
      REAL(q) A(*)
      INTEGER ISIGN   !  direction of fft
      ! local variables
      INTEGER(c_int) EFFORT

      PROFILING_START('fftbas_mpi')

      EFFORT=FFTW_ESTIMATE
#ifdef use_fftw_plan_effort
      IF (FFTW_PLAN_EFFORT==1) EFFORT=FFTW_MEASURE
#endif
      CALL FFTBAS_PLAN_MPI(A,GRID,ISIGN,EFFORT)

      PROFILING_STOP('fftbas_mpi')

    END SUBROUTINE


!> @details @ref openmp :
!> under OpenMP this subroutine calls dfftw_plan_with_nthreads.
!> Note: this subroutine may be called from within an OpenMP parallel
!> region as well: threadsafety is ensured by an explicit CRITICAL region
!> (VASP_FFT_PLAN_CREATE_DESTROY).
!
    SUBROUTINE FFTBAS_PLAN_MPI(A,GRID,ISIGN,EFFORT)
#ifdef _OPENACC
      USE mopenacc
      USE mcufft
#endif
      USE prec
      USE smart_allocate
      USE mpimy
      USE mgrid

      USE iso_c_binding

      IMPLICIT NONE

#ifndef __NEC__
      include 'fftw3.f'
#else // __NEC__
      include 'aslfftw3.f'
#endif // __NEC__

      TYPE (grid_3d) GRID
      REAL(q) A(*)
      INTEGER ISIGN          !  direction of fft
      INTEGER(c_int) EFFORT  !  planning effort
      ! local variables
      COMPLEX(q),POINTER,SAVE ::  RCVBUF(:),SNDBUF(:)
      TYPE(c_ptr) :: planx, plany, planz
      INTEGER :: NX,NY,NZ,IDX,IDY,IDZ

!$    INTEGER THREADS
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS
#ifdef _OPENACC
      COMPLEX(q), POINTER  :: C(:)
      TYPE(acc_cufft_plan) :: ACC_PLAN
      INTEGER :: IERR
#endif
      PROFILING_START('fftbas_plan_mpi')

      NX=GRID%NGPTAR(1)
      NY=GRID%NGPTAR(2)
      NZ=GRID%NGPTAR(3)

      CALL SMART_ALLOCATE_COMPLEX(RCVBUF,GRID%MPLWV)
      CALL SMART_ALLOCATE_COMPLEX(SNDBUF,GRID%MPLWV)

      IDX=NX
      IDY=NY
      IDZ=NZ

#ifdef _OPENACC
      IF (ACC_IS_PRESENT(A).AND.ACC_EXEC_ON) THEN
!$ACC DATA CREATE(RCVBUF,SNDBUF) ASYNC(ACC_ASYNC_Q)
         CALL c_f_pointer(c_loc(A),C,[GRID%MPLWV])

         IERR=ACC_CUFFT_GETPLAN(GRID,INT(ACC_ASYNC_Q),1,ACC_PLAN)
         IF (IERR/=0) THEN
            CALL vtutor%error("FFTBAS_PLAN_MPI: no cuFFT plan found for ACC_ASYNC_Q= " // &
               str(INT(ACC_ASYNC_Q)))
         ENDIF
!=======================================================================
! do the transformation forward (q->r)
!=======================================================================
         IF (ISIGN==1) THEN
! transformation along first dimension:
!$ACC HOST_DATA USE_DEVICE(C)
            IF (GRID%RC%NCOL>0) IERR=IERR+cufftExecZ2Z(ACC_PLAN%CUFFT_RC_Z2Z,C,C,CUFFT_INVERSE)
!$ACC END HOST_DATA
            CALL MAP_FORWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)
! transformation along second dimension:
!$ACC HOST_DATA USE_DEVICE(C)
            IF (GRID%IN%NCOL>0) IERR=IERR+cufftExecZ2Z(ACC_PLAN%CUFFT_IN_Z2Z,C,C,CUFFT_INVERSE)
!$ACC END HOST_DATA
            CALL MAP_FORWARD(A(1), GRID%RL_FFT%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)
! transformation along third dimension:
            IF (NZ/2+1==GRID%NGZ_rd) THEN
!$ACC HOST_DATA USE_DEVICE(C,A)
               IF (GRID%RL_FFT%NCOL>0) IERR=IERR+cufftExecZ2D(ACC_PLAN%CUFFT_RL_Z2D,C,A)
!$ACC END HOST_DATA
            ELSE
!$ACC HOST_DATA USE_DEVICE(C)
               IF (GRID%RL_FFT%NCOL>0) IERR=IERR+cufftExecZ2Z(ACC_PLAN%CUFFT_RL_Z2Z,C,C,CUFFT_INVERSE)
!$ACC END HOST_DATA
            ENDIF
!=======================================================================
! do the transformation backward (r->q)
!=======================================================================
         ELSE
! transformation along third dimension:
            IF (NZ/2+1==GRID%NGZ_rd) THEN
!$ACC HOST_DATA USE_DEVICE(C,A)
               IF (GRID%RL_FFT%NCOL>0) IERR=IERR+cufftExecD2Z(ACC_PLAN%CUFFT_RL_D2Z,A,C)
!$ACC END HOST_DATA
            ELSE
!$ACC HOST_DATA USE_DEVICE(C)
               IF (GRID%RL_FFT%NCOL>0) IERR=IERR+cufftExecZ2Z(ACC_PLAN%CUFFT_RL_Z2Z,C,C,CUFFT_FORWARD)
!$ACC END HOST_DATA
            ENDIF
            CALL MAP_BACKWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)
! transformation along second dimension:
!$ACC HOST_DATA USE_DEVICE(C)
            IF (GRID%IN%NCOL>0) IERR=IERR+cufftExecZ2Z(ACC_PLAN%CUFFT_IN_Z2Z,C,C,CUFFT_FORWARD)
!$ACC END HOST_DATA
            CALL MAP_BACKWARD(A(1), GRID%RC%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)
! transformation along first dimension:
!$ACC HOST_DATA USE_DEVICE(C)
            IF (GRID%RC%NCOL>0) IERR=IERR+cufftExecZ2Z(ACC_PLAN%CUFFT_RC_Z2Z,C,C,CUFFT_FORWARD)
!$ACC END HOST_DATA
         ENDIF
!$ACC END DATA
         PROFILING_STOP('fftbas_plan_mpi')
         RETURN
      ENDIF
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
!=======================================================================
! create plans
!=======================================================================
!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!$    CALL dfftw_plan_with_nthreads(THREADS) 
      IF (ISIGN==1) THEN
         IF (GRID%RC%NCOL > 0) THEN
            CALL dfftw_plan_many_dft(planx, 1, NX , GRID%RC%NCOL, &
                                A(1), NX, 1 , IDX, &
                                A(1), NX, 1 , IDX, &
                                FFTW_BACKWARD, EFFORT)
         ENDIF

         IF (GRID%IN%NCOL > 0) THEN
            CALL dfftw_plan_many_dft(plany, 1, NY , GRID%IN%NCOL, &
                                A(1), NY, GRID%IN%NCOL, 1 , &
                                A(1), NY, GRID%IN%NCOL, 1 , &
                                FFTW_BACKWARD, EFFORT)
         ENDIF

         IF (GRID%RL_FFT%NCOL > 0) THEN
            IF (NZ/2+1==GRID%NGZ_rd) THEN
!              WRITE(*,*) 'detected real to complex'
              CALL dfftw_plan_many_dft_c2r(planz, 1, NZ , GRID%RL_FFT%NCOL, &
                                A(1), NZ, 1, (IDZ+2)/2 , &
                                A(1), NZ, 1, IDZ+2 , &
                                EFFORT)
            ELSE
!              WRITE(*,*) 'complex to complex'
              CALL dfftw_plan_many_dft(planz, 1, NZ , GRID%RL_FFT%NCOL, &
                                A(1), NZ, 1, IDZ , &
                                A(1), NZ, 1, IDZ , &
                                FFTW_BACKWARD, EFFORT)
            ENDIF
         ENDIF
      ELSE
         IF (GRID%RL_FFT%NCOL > 0) THEN
            IF (NZ/2+1==GRID%NGZ_rd) THEN
!              WRITE(*,*) 'detected inverse real to complex'
              CALL dfftw_plan_many_dft_r2c(planz, 1, NZ , GRID%RL_FFT%NCOL, &
                                A(1), NZ, 1, IDZ+2 , &
                                A(1), NZ, 1, (IDZ+2)/2 , &
                                EFFORT)
            ELSE
!              WRITE(*,*) 'detected inverse complex to complex'
              CALL dfftw_plan_many_dft(planz, 1, NZ , GRID%RL_FFT%NCOL, &
                                A(1), NZ, 1, IDZ , &
                                A(1), NZ, 1, IDZ , &
                                FFTW_FORWARD, EFFORT)
            ENDIF
         ENDIF

         IF (GRID%IN%NCOL > 0) THEN
            CALL dfftw_plan_many_dft(plany, 1, NY , GRID%IN%NCOL, &
                                A(1), NY, GRID%IN%NCOL, 1 , &
                                A(1), NY, GRID%IN%NCOL, 1 , &
                                FFTW_FORWARD, EFFORT)
         ENDIF

         IF (GRID%RC%NCOL > 0) THEN
            CALL dfftw_plan_many_dft(planx, 1, NX , GRID%RC%NCOL, &
                                A(1), NX, 1 , IDX, &
                                A(1), NX, 1 , IDX, &
                                FFTW_FORWARD, EFFORT)
         ENDIF
      ENDIF
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!=======================================================================
! do the transformation forward (q->r)
!=======================================================================
      IF (ISIGN==1) THEN
! transformation along first dimension:
         PROFILING_START('dfftw_execute')
         IF (GRID%RC%NCOL>0) THEN
            CALL dfftw_execute_dft(planx,A(1),A(1))
         ENDIF
         PROFILING_STOP('dfftw_execute')
         CALL MAP_FORWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)
! transformation along second dimension:
         PROFILING_START('dfftw_execute')
         IF (GRID%IN%NCOL>0) THEN
            CALL dfftw_execute_dft(plany,A(1),A(1))
         ENDIF
         PROFILING_STOP('dfftw_execute')
         CALL MAP_FORWARD(A(1), GRID%RL_FFT%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)
! transformation along third dimension:
         PROFILING_START('dfftw_execute')
         IF (GRID%RL_FFT%NCOL>0) THEN
            IF (NZ/2+1==GRID%NGZ_rd) THEN
               CALL dfftw_execute_dft_c2r(planz,A(1),A(1))
            ELSE
               CALL dfftw_execute_dft(planz,A(1),A(1))
            ENDIF
         ENDIF
         PROFILING_STOP('dfftw_execute')
!=======================================================================
! do the transformation backward (r->q)
!=======================================================================
      ELSE
! transformation along third dimension:
         PROFILING_START('dfftw_execute')
         IF (GRID%RL_FFT%NCOL>0) THEN
            IF (NZ/2+1==GRID%NGZ_rd) THEN
               CALL dfftw_execute_dft_r2c(planz,A(1),A(1))
            ELSE
               CALL dfftw_execute_dft(planz,A(1),A(1))
            ENDIF
         ENDIF
         PROFILING_STOP('dfftw_execute')
         CALL MAP_BACKWARD(A(1), GRID%IN%NALLOC, SNDBUF(1), RCVBUF(1), GRID%IN_RL, GRID%COMM)
! transformation along second dimension:
         PROFILING_START('dfftw_execute')
         IF (GRID%IN%NCOL>0) THEN
            CALL  dfftw_execute_dft(plany,A(1),A(1))
         ENDIF
         PROFILING_STOP('dfftw_execute')
         CALL MAP_BACKWARD(A(1), GRID%RC%NALLOC, SNDBUF(1), RCVBUF(1), GRID%RC_IN, GRID%COMM)
! transformation along first dimension:
         PROFILING_START('dfftw_execute')
         IF (GRID%RC%NCOL>0) THEN
            CALL dfftw_execute_dft(planx,A(1),A(1))
         ENDIF
         PROFILING_STOP('dfftw_execute')
      ENDIF

!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
      IF (GRID%RC%NCOL > 0) THEN
         call dfftw_destroy_plan(planx)
      ENDIF

      IF (GRID%IN%NCOL > 0) THEN
         call dfftw_destroy_plan(plany)
      ENDIF

      IF (GRID%RL_FFT%NCOL > 0) THEN
         call dfftw_destroy_plan(planz)
      ENDIF
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

      PROFILING_STOP('fftbas_plan_mpi')

      RETURN
    END SUBROUTINE FFTBAS_PLAN_MPI


!=======================================================================
!   3-d parallel real to complex fast fourier transformation using 
!   fftw-kernels
!   communication routines and set of communication routines
!   in fftmpi_map.F written by Georg Kresse
!
!     +1  q->r   vr= sum(q) vq exp(+iqr)
!     -1  r->q   vq= sum(r) vr exp(-iqr)
!
!=======================================================================

!
!  this subroutine calls the FFTBAS_PLAN_MPI routine  with FFTW_ESTIMATE
!  the FFTBAS_PLAN_MPI routine detects automatically real to complex
!  transforms and handles them accordingly
!
    SUBROUTINE FFTBRC_MPI(A,GRID,ISIGN)
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
      REAL(q) A(*)
      INTEGER ISIGN   !  direction of fft
      ! local variables
      INTEGER(c_int) EFFORT

      PROFILING_START('fftbrc_mpi')

      EFFORT=FFTW_ESTIMATE
#ifdef use_fftw_plan_effort
      IF (FFTW_PLAN_EFFORT==1) EFFORT=FFTW_MEASURE
#endif
      CALL FFTBAS_PLAN_MPI(A,GRID,ISIGN,EFFORT)

      PROFILING_STOP('fftbrc_mpi')

    END SUBROUTINE FFTBRC_MPI


!***********************************************************************
! FROM HERE ON THE ROUTINES ARE IDENTICAL TO fftmpi.F
!***********************************************************************


!************************* SUBROUTINE FFTINI ***************************
!
!  if necessary this routine performes initialization
!  for FFTWAV and FFTEXT
!  usually this is only necessary for the Gamma point only
!  1-kpoint version
!
!   FFTSCA(.,1) is the scaling factor for extracting the wavefunction
!               from the FFT grid (FFTEXT)
!   FFTSCA(.,2) is the scaling factor for puting the wavefunction on
!               the grid
!***********************************************************************

    SUBROUTINE  FFTINI_MPI(NINDPW,NPLWKP,NKPTS,NRPLW,GRID)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (grid_3d)  GRID
      DIMENSION NPLWKP(NKPTS)
      DIMENSION NINDPW(NRPLW,NKPTS)

      IF (GRID%REAL2CPLX) THEN
         IF (GRID%RL%NFAST==1) THEN
            CALL FFTINI(NINDPW,NPLWKP,NKPTS,NRPLW,GRID)
            RETURN
         ENDIF
         
         IF (NKPTS>1) THEN
            CALL vtutor%error("FFT3D: real version works only for 1 k-point")
         ENDIF
         
         NK=1
         NPL=NPLWKP(NK)
         NULLIFY(GRID%FFTSCA)
         ALLOCATE(GRID%FFTSCA(NPL,2))
         
         DO N=1,NPL
            IND=NINDPW(N,NK)
            N1= MOD((IND-1),GRID%RC%NROW)+1
            NC= (IND-1)/GRID%RC%NROW+1
            N2= GRID%RC%I2(NC)
            N3= GRID%RC%I3(NC)
            
            FACTM=SQRT(2._q)
            IF (N1==1 .AND. N2==1 .AND. N3==1) FACTM=1
            GRID%FFTSCA(N,1)= FACTM
            GRID%FFTSCA(N,2)= 1/FACTM
            ! this statment is required
            ! because for z==0 only half of the FFT components are set
            ! upon calling FFTWAV
            IF (N3==1) GRID%FFTSCA(N,2)=FACTM
         ENDDO
      END IF
      RETURN
    END SUBROUTINE FFTINI_MPI

!************************* SUBROUTINE FFTWAV ***************************
!
!  this subroutine transforms a wavefunction C defined  within  the
!  cutoff-sphere to real space CR
! MIND:
! for the real version (gamma point only) it is assumed
! that the wavefunctions at NGZ != 0
! are multiplied by a factor sqrt(2) on the linear grid
! this factor has to be removed before the FFT transformation !
! (scaling with   FFTSCA(M,2))
!
!> @details @ref openmp :
!> the loops that fill CR with the reciprocal space components of
!> an orbital are distributed over all available OpenMP threads.
!
!***********************************************************************

    SUBROUTINE FFTWAV_MPI(NPL,NINDPW,CR,C,GRID)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      TYPE (grid_3d)     GRID
      COMPLEX(q):: C(NPL), CR(GRID%MPLWV)
      DIMENSION NINDPW(NPL)

      PROFILING_START('fftwav_mpi')

      IF (GRID%RL%NFAST==1) THEN
         CALL FFTWAV(NPL,NINDPW,CR,C,GRID)
         PROFILING_STOP('fftwav_mpi')
         RETURN
      ENDIF

      IF (GRID%LREAL) THEN
!DIR$ IVDEP
!OCL NOVREC
!$ACC PARALLEL LOOP PRESENT(CR) IF (ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(M) SHARED(CR,GRID)
         DO M=1,GRID%RL%NCOL*GRID%NGZ/2
            CR(M)=0.0_q
         ENDDO
NOACC !$OMP END PARALLEL DO
      ELSE
!DIR$ IVDEP
!OCL NOVREC
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(M) SHARED(CR,GRID)
!$ACC PARALLEL LOOP PRESENT(CR) IF (ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
         DO M=1,GRID%RL%NCOL*GRID%NGZ
            CR(M)=0.0_q
         ENDDO
NOACC !$OMP END PARALLEL DO
      ENDIF

      IF (GRID%REAL2CPLX) THEN
!DIR$ IVDEP
!OCL NOVREC
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(M) SHARED(NPL,C,CR,NINDPW,GRID)
!$ACC PARALLEL LOOP PRESENT(CR,C,NINDPW,GRID%FFTSCA) IF (ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
         DO M=1,NPL
            CR(NINDPW(M))=C(M)*GRID%FFTSCA(M,2)
         ENDDO
NOACC !$OMP END PARALLEL DO
      ELSE
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(M) SHARED(NPL,C,CR,NINDPW)
!$ACC PARALLEL LOOP PRESENT(CR,C,NINDPW) IF (ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
         DO M=1,NPL
            CR(NINDPW(M))=C(M)
         ENDDO
NOACC !$OMP END PARALLEL DO
      ENDIF
      CALL FFT3D_MPI(CR,GRID,1)

      PROFILING_STOP('fftwav_mpi')

      RETURN
    END SUBROUTINE FFTWAV_MPI

!************************* SUBROUTINE FFTEXT ***************************
!
! this subroutine performes a FFT to reciprocal space and extracts data
! from the FFT-mesh
! MIND:
! for the real version (gamma point only) it is assumed
! that the wavefunctions at NGX != 0
! are multiplied by a factor sqrt(2) on the linear grid
! this factor has to be applied after the FFT transformation !
!  (scaling with   FFTSCA(M))
!
!> @details @ref openmp :
!> the loops that extract the reciprocal space components of an orbital
!> from the full FFT grid are distributed over all available OpenMP
!> threads.
!
!***********************************************************************

    SUBROUTINE FFTEXT_MPI(NPL,NINDPW,CR,C,GRID,LADD)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      DIMENSION C(NPL),CR(GRID%MPLWV)
      DIMENSION NINDPW(NPL)
      LOGICAL   LADD

      PROFILING_START('fftext_mpi')

      CALL FFT3D_MPI(CR,GRID,-1)

      IF (LADD .AND. GRID%REAL2CPLX) THEN
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID%FFTSCA) IF (ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(M) SHARED(NPL,C,CR,NINDPW, GRID)
         DO M=1,NPL
            C(M)=C(M)+CR(NINDPW(M))*GRID%FFTSCA(M,1)
         ENDDO
NOACC !$OMP END PARALLEL DO
      ELSE IF (LADD .AND. .NOT. GRID%REAL2CPLX) THEN
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW) IF (ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(M) SHARED(NPL,CR,NINDPW) REDUCTION(+:C)
         DO M=1,NPL
            C(M)=C(M)+CR(NINDPW(M))
         ENDDO
NOACC !$OMP END PARALLEL DO
      ELSE IF (GRID%REAL2CPLX) THEN
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID%FFTSCA) IF (ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(M) SHARED(NPL,C,CR,NINDPW,GRID)
        DO M=1,NPL
          C(M)=CR(NINDPW(M))*GRID%FFTSCA(M,1)
        ENDDO
NOACC !$OMP END PARALLEL DO
      ELSE
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW) IF (ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
!DIR$ IVDEP
!OCL NOVREC
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(M) SHARED(NPL,C,CR,NINDPW)
        DO M=1,NPL
          C(M)=CR(NINDPW(M))
        ENDDO
NOACC !$OMP END PARALLEL DO
      ENDIF

      PROFILING_STOP('fftext_mpi')

      RETURN
    END SUBROUTINE FFTEXT_MPI


!===============================================================================
!
!    3-d fast fourier transform (possibly real to complex and vice versa)
!    for chardensities and potentials
!     +1  q->r   vr= sum(q) vq exp(+iqr)    (might be complex to real)
!     -1  r->q   vq= sum(r) vr exp(-iqr)    (might be real to complex)
!
!===============================================================================

    SUBROUTINE FFT3D_MPI(C,GRID,ISN)
      USE prec
      USE mpimy
      USE mgrid
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      
      TYPE (grid_3d)   GRID
      REAL(q) C(*)

      PROFILING_START('fft3d_mpi')

      NX=GRID%NGPTAR(1)
      NY=GRID%NGPTAR(2)
      NZ=GRID%NGPTAR(3)

!-------------------------------------------------------------------------------
! use serial FFT
!-------------------------------------------------------------------------------
      IF (GRID%RL%NFAST==1 .AND. GRID%RL_FFT%NFAST==1) THEN
         CALL FFT3D(C, GRID, ISN)
      ELSE IF (GRID%RL%NFAST==1) THEN
!-------------------------------------------------------------------------------
! parallel FFT with serial data layout  (GRID%RL%NFAST==1)
!-------------------------------------------------------------------------------
       ! complex to complex case
         IF (.NOT. GRID%REAL2CPLX) THEN
            IF (.NOT. (NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
               CALL vtutor%bug("internal error 1 in FFT3D_MPI: something not properly set " // &
                  str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) &
                  // " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
                  str(GRID%NGZ_rd), __FILE__, __LINE__)
            ENDIF
            
            IF (ISN==1) THEN
               ! q->r FFT
               CALL FFTBAS_MPI(C,GRID,ISN)
               ! bring real space data from parallel to serial data layout
               IF (GRID%COMM%NODE_ME==1) THEN
                  CALL FFTPAR_TO_SER(GRID%NGX, GRID%NGY, GRID%NGZ, C)
               ENDIF
            ELSE
               ! bring real space data from serial to parallel data layout
               IF (GRID%COMM%NODE_ME==1) THEN
                  CALL FFTSER_TO_PAR(GRID%NGX, GRID%NGY, GRID%NGZ, C)
               ENDIF
               ! r->q FFT
               CALL FFTBAS_MPI(C,GRID,ISN)
            ENDIF
            
            ! complex to real case
         ELSE IF (GRID%LREAL) THEN
            IF (.NOT. (NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ/2+1==GRID%NGZ_rd) ) THEN
               CALL vtutor%bug("internal error 2 in FFT3D_MPI: something not properly set " // &
                  str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) &
                  // " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
                  str(GRID%NGZ_rd), __FILE__, __LINE__)
            ENDIF

            !  in real space the first dimension in VASP is NGZ (REAL data)
            !  but the FFT requires NGZ+2 (real data)
            !  therefore some data movement is required
            IF (ISN==1) THEN
               ! q->r FFT
               CALL FFTBRC_MPI(C,GRID,ISN)
               ! x-lines (go from stride NZ+2 to NZ)
               CALL RESTRIDE_Q2R(GRID%RL_FFT%NCOL,NZ,C)
               ! bring real space data from parallel to serial data layout
               IF (GRID%COMM%NODE_ME==1) THEN
                  CALL FFTPAR_TO_SER_REAL(GRID%NGX, GRID%NGY, GRID%NGZ, C)
               ENDIF
            ELSE
               ! bring real space data from serial to parallel data layout
               IF (GRID%COMM%NODE_ME==1) THEN
                  CALL FFTSER_TO_PAR_REAL(GRID%NGX, GRID%NGY, GRID%NGZ, C)
               ENDIF
               ! x-lines (go from stride NZ to NZ+2)
               CALL RESTRIDE_R2Q(GRID%RL_FFT%NCOL,NZ,C)
               ! r->q FFT
               CALL FFTBRC_MPI(C,GRID,ISN)
            ENDIF
         ELSE
            CALL vtutor%error("ERROR in FFT3D_MPI: this version does not support the required half &
               &grid mode \n " // str(NX) // " " // str(NY) // " " // str(NZ) // "\n " // &
               str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // str(GRID%NGZ_rd))
         ENDIF
!-------------------------------------------------------------------------------
!  complex parallel FFT
!-------------------------------------------------------------------------------
      ELSE IF (.NOT. GRID%REAL2CPLX) THEN
         IF (.NOT. (NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("internal error 3 in FFT3D_MPI: something not properly set " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) // &
               " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
               str(GRID%NGZ_rd), __FILE__, __LINE__)
         ENDIF
         CALL FFTBAS_MPI(C,GRID,ISN)
!-------------------------------------------------------------------------------
!  real to complex parallel FFT
!-------------------------------------------------------------------------------
      ELSE IF (GRID%LREAL) THEN
         IF (.NOT.(NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ/2+1==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("internal error 4 in FFT3D_MPI: something not properly set " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) // &
               " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
               str(GRID%NGZ_rd), __FILE__, __LINE__)
         ENDIF
         
         !  in real space the first dimension in VASP is NGZ (REAL data)
         !  but the FFT requires NGZ+2 (real data)
         !  therefore some data movement is required
         IF (ISN==1) THEN
            ! q->r FFT
            CALL FFTBRC_MPI(C,GRID,ISN)
            ! concat z-lines (go from stride NZ+2 to NZ)
            CALL RESTRIDE_Q2R(GRID%RL%NCOL,NZ,C)
         ELSE
            ! concat z-lines (go from stride NZ to NZ+2)
            CALL RESTRIDE_R2Q(GRID%RL%NCOL,NZ,C)
            ! r->q FFT
            CALL FFTBRC_MPI(C,GRID,ISN)
         ENDIF
!-------------------------------------------------------------------------------
!  real to complex parallel FFT with complex storage layout in real space
!-------------------------------------------------------------------------------
      ELSE
         IF (.NOT.(NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ/2+1==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("internal error 5 in FFT3D_MPI: something not properly set " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // str(NX) // " " // str(NY) // &
               " " // str(NZ) // "\n " // str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // &
               str(GRID%NGZ_rd), __FILE__, __LINE__)
         ENDIF

         IF (ISN==1) THEN
            ! q->r FFT
            CALL FFTBRC_MPI(C,GRID,ISN)
            ! concat z-lines (go from stride NZ+2 to NZ)
            CALL RESTRIDE_Q2R_CMPLX(GRID%RL%NCOL,NZ,C)
         ELSE
            ! z-lines (go from complex stride NZ to real stride NZ+2)
            CALL RESTRIDE_R2Q_CMPLX(GRID%RL%NCOL,NZ,C)
            ! r->q FFT
            CALL FFTBRC_MPI(C,GRID,ISN)
         ENDIF
      ENDIF

      PROFILING_STOP('fft3d_mpi')

    END SUBROUTINE FFT3D_MPI


!=======================================================================
!   this routine returns the next correct setting for the
!   three dimensional FFT
!=======================================================================

    SUBROUTINE FFTCHK_MPI(NFFT)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION NFFT(3)
      LOGICAL FFTCH1_MPI

      DO IND=1,3
200      CONTINUE
         IF (FFTCH1_MPI(NFFT(IND))) CYCLE
         NFFT(IND)=NFFT(IND)+1
         GOTO 200
100   ENDDO
    END SUBROUTINE FFTCHK_MPI
    
    LOGICAL FUNCTION FFTCH1_MPI(NIN)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      PARAMETER (NFACT=4)
      DIMENSION IFACT(NFACT),NCOUNT(NFACT)
      DATA      IFACT /2,3,5,7/
      N=NIN
      DO 100 I=1,NFACT
         NCOUNT(I)=0
120      NEXT=N/IFACT(I)
         IF (NEXT*IFACT(I)==N) THEN
            N=NEXT
            NCOUNT(I)=NCOUNT(I)+1
            GOTO 120
         ENDIF
100   ENDDO
      IF (N==1 .AND. (NCOUNT(1)/=0)) &
           &  THEN
         FFTCH1_MPI=.TRUE.
      ELSE
         FFTCH1_MPI=.FALSE.
      ENDIF
      RETURN
    END FUNCTION FFTCH1_MPI

!=======================================================================
! 
! change data layout from parallel to serial data layout
! and vice versa for complex and real arrays
! operates usually in real space
!
!=======================================================================


    SUBROUTINE FFTPAR_TO_SER(NGX, NGY, NGZ, CORIG)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec

      INTEGER NGX, NGY, NGZ
      COMPLEX(q) :: CORIG(NGX*NGY*NGZ)
    ! local
      INTEGER IX, IY, IZ
      COMPLEX(q) :: C(NGX*NGY*NGZ)

!$ACC ENTER DATA CREATE(C) __IF_ASYNC__
!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CORIG) __IF_ASYNC__
      DO IX=0,NGX-1
         DO IY=0,NGY-1
!DIR$ IVDEP
!OCL NOVREC
            DO IZ=0,NGZ-1
               ! C(IX,IY,IZ)=CORIG(IZ,IX,IY)
               C(1+IX+NGX*(IY+NGY*IZ))=CORIG(1+IZ+NGZ*(IX+NGX*IY))
            ENDDO
         ENDDO
      ENDDO
!DIR$ IVDEP
!OCL NOVREC
!$ACC PARALLEL LOOP PRESENT(C,CORIG) __IF_ASYNC__
      DO IX=1,NGX*NGY*NGZ
         CORIG(IX)=C(IX)
      ENDDO
!$ACC EXIT DATA DELETE(C) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
    END SUBROUTINE FFTPAR_TO_SER

    SUBROUTINE FFTPAR_TO_SER_REAL(NGX, NGY, NGZ, CORIG)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec

      INTEGER NGX, NGY, NGZ
      REAL(q) :: CORIG(NGX*NGY*NGZ)
    ! local
      INTEGER IX, IY, IZ
      REAL(q) :: C(NGX*NGY*NGZ)

!$ACC ENTER DATA CREATE(C) __IF_ASYNC__
!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CORIG) __IF_ASYNC__
      DO IX=0,NGX-1
         DO IY=0,NGY-1
!DIR$ IVDEP
!OCL NOVREC
            DO IZ=0,NGZ-1
               ! C(IX,IY,IZ)=CORIG(IZ,IX,IY)
               C(1+IX+NGX*(IY+NGY*IZ))=CORIG(1+IZ+NGZ*(IX+NGX*IY))
            ENDDO
         ENDDO
      ENDDO
!DIR$ IVDEP
!OCL NOVREC
!$ACC PARALLEL LOOP PRESENT(C,CORIG) __IF_ASYNC__
      DO IX=1,NGX*NGY*NGZ
         CORIG(IX)=C(IX)
      ENDDO
!$ACC EXIT DATA DELETE(C) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
    END SUBROUTINE FFTPAR_TO_SER_REAL


    SUBROUTINE FFTSER_TO_PAR(NGX, NGY, NGZ, CORIG)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec

      INTEGER NGX, NGY, NGZ
      COMPLEX(q) :: CORIG(NGX*NGY*NGZ)
    ! local
      INTEGER IX, IY, IZ
      COMPLEX(q) :: C(NGX*NGY*NGZ)

!$ACC ENTER DATA CREATE(C) __IF_ASYNC__
!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CORIG) __IF_ASYNC__
      DO IX=0,NGX-1
         DO IY=0,NGY-1
!DIR$ IVDEP
!OCL NOVREC
            DO IZ=0,NGZ-1
               ! C(IZ,IX,IY)=CORIG(IX,IY,IZ)
               C(1+IZ+NGZ*(IX+NGX*IY))=CORIG(1+IX+NGX*(IY+NGY*IZ))
            ENDDO
         ENDDO
      ENDDO
!DIR$ IVDEP
!OCL NOVREC
!$ACC PARALLEL LOOP PRESENT(C,CORIG) __IF_ASYNC__
      DO IX=1,NGX*NGY*NGZ
         CORIG(IX)=C(IX)
      ENDDO
!$ACC EXIT DATA DELETE(C) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
    END SUBROUTINE FFTSER_TO_PAR

    SUBROUTINE FFTSER_TO_PAR_REAL(NGX, NGY, NGZ, CORIG)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec

      INTEGER NGX, NGY, NGZ
      REAL(q) :: CORIG(NGX*NGY*NGZ)
    ! local
      INTEGER IX, IY, IZ
      REAL(q) :: C(NGX*NGY*NGZ)

!$ACC ENTER DATA CREATE(C) __IF_ASYNC__
!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(C,CORIG) __IF_ASYNC__
      DO IX=0,NGX-1
         DO IY=0,NGY-1
!DIR$ IVDEP
!OCL NOVREC
            DO IZ=0,NGZ-1
               ! C(IZ,IX,IY)=CORIG(IX,IY,IZ)
               C(1+IZ+NGZ*(IX+NGX*IY))=CORIG(1+IX+NGX*(IY+NGY*IZ))
            ENDDO
         ENDDO
      ENDDO
!DIR$ IVDEP
!OCL NOVREC
!$ACC PARALLEL LOOP PRESENT(C,CORIG) __IF_ASYNC__
      DO IX=1,NGX*NGY*NGZ
         CORIG(IX)=C(IX)
      ENDDO
!$ACC EXIT DATA DELETE(C) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
    END SUBROUTINE FFTSER_TO_PAR_REAL

