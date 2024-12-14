#include "symbol.inc"

#undef FFTGRIDPLAN
#undef FFTMAKEPLAN
#undef FFT3D
#undef FFT3DR
#undef FFTCHK
#undef FFT3RC
#undef FFTEXT
#undef FFTWAV
#undef FFTINI

!****************** SUBROUTINE FFTWAV_MU *******************************
!
!***********************************************************************

    SUBROUTINE FFTWAV_MU(NPL,N,NINDPW,CR,LDR,C,LDC,GRID)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE mgrid

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: CR(LDR,*),C(LDC,*)
      INTEGER    :: NINDPW(NPL)
      INTEGER    :: NPL, N, LDR, LDC

      ! local variables
      INTEGER :: I,J

#ifdef _OPENACC
      LOGICAL :: ACC_ACTIVE
      ACC_ACTIVE=ACC_IS_PRESENT(NINDPW).AND.ACC_IS_PRESENT(CR).AND. &
          ACC_IS_PRESENT(C).AND.ACC_IS_PRESENT(GRID,1).AND.ACC_EXEC_ON
#endif

      PROFILING_START('fftwav_mu')

      ! zero all elements on the grid
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,GRID) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,GRID,CR) PRIVATE(I,J)
      DO I=1,N
         DO J=1,GRID%NGX_rd*GRID%NGY_rd*GRID%NGZ_rd
            CR(J,I)=(0.0_q,0.0_q)
         ENDDO
      ENDDO
NOACC !$OMP END PARALLEL DO

      ! now fill in non zero elements from
      ! within the radial cutoff sphere
      IF (GRID%REAL2CPLX) THEN

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,NINDPW,GRID) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,CR,NINDPW,C,GRID) PRIVATE(I,J)
         DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
            DO J=1,NPL
               CR(NINDPW(J),I)=C(J,I)*GRID%FFTSCA(J,2)
            ENDDO
         ENDDO
NOACC !$OMP END PARALLEL DO
#ifdef inversion_through_mid_point
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,GRID) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,CR,C,GRID) PRIVATE(I,J)
         DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
            DO J=1,SIZE(GRID%IND_IN_SPHERE)
               CR(GRID%NINDPWCONJG(J),I)=CONJG(C(GRID%IND_IN_SPHERE(J),I))*GRID%FFTSCA(GRID%IND_IN_SPHERE(J),2)
            ENDDO
         ENDDO
NOACC !$OMP END PARALLEL DO
#endif

      ELSE

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,NINDPW) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,CR,NINDPW,C) PRIVATE(I,J)
         DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
            DO J=1,NPL
               CR(NINDPW(J),I)=C(J,I)
            ENDDO
         ENDDO
NOACC !$OMP END PARALLEL DO

      ENDIF

      CALL FFT3D_MU(N,CR,LDR,GRID,1)

      PROFILING_STOP('fftwav_mu')

    END SUBROUTINE FFTWAV_MU


!****************** SUBROUTINE FFTWAV_USEINV_MU ************************
!
!***********************************************************************

    SUBROUTINE FFTWAV_USEINV_MU(NPL,N,NINDPW,NINDPW_INV,FFTSCA,CR,LDR,C,LDC,GRID)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE mgrid

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: CR(LDR,*),C(LDC,*)
      REAL(q)    :: FFTSCA(NPL)
      INTEGER    :: NINDPW(NPL),NINDPW_INV(NPL)
      INTEGER    :: NPL, N, LDR, LDC

      ! local variables
      INTEGER :: I,J

#ifdef _OPENACC
      LOGICAL :: ACC_ACTIVE
      ACC_ACTIVE=ACC_IS_PRESENT(NINDPW).AND.ACC_IS_PRESENT(NINDPW_INV).AND.ACC_IS_PRESENT(FFTSCA).AND. &
          ACC_IS_PRESENT(CR).AND.ACC_IS_PRESENT(C).AND.ACC_IS_PRESENT(GRID,1).AND.ACC_EXEC_ON
#endif

      PROFILING_START('fftwav_useinv_mu')

      ! zero all elements on the grid
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,GRID) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,GRID,CR) PRIVATE(I,J)
      DO I=1,N
         DO J=1,GRID%NGX_rd*GRID%NGY_rd*GRID%NGZ_rd
            CR(J,I)=(0.0_q,0.0_q)
         ENDDO
      ENDDO
NOACC !$OMP END PARALLEL DO

      ! now fill in non zero elements from
      ! within the radial cutoff sphere
      IF (GRID%REAL2CPLX) THEN

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,NINDPW,GRID) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,CR,NINDPW,C,GRID) PRIVATE(I,J)
         DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
            DO J=1,NPL
               CR(NINDPW(J),I)=C(J,I)*GRID%FFTSCA(J,2)
            ENDDO
         ENDDO
NOACC !$OMP END PARALLEL DO
#ifdef inversion_through_mid_point
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,GRID) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,CR,C,GRID) PRIVATE(I,J)
         DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
            DO J=1,SIZE(GRID%IND_IN_SPHERE)
               CR(GRID%NINDPWCONJG(J),I)=CONJG(C(GRID%IND_IN_SPHERE(J),I))*GRID%FFTSCA(GRID%IND_IN_SPHERE(J),2)
            ENDDO
         ENDDO
NOACC !$OMP END PARALLEL DO
#endif

      ELSE

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,NINDPW_INV,FFTSCA) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,CR,NINDPW_INV,C,FFTSCA) PRIVATE(I,J)
         DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
            DO J=1,NPL
               CR(NINDPW_INV(J),I)=CONJG(C(J,I))*FFTSCA(J)
            ENDDO
         ENDDO
NOACC !$OMP END PARALLEL DO

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CR,C,NINDPW,FFTSCA) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,CR,NINDPW,C,FFTSCA) PRIVATE(I,J)
         DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
            DO J=1,NPL
               CR(NINDPW(J),I)=C(J,I)*FFTSCA(J)
            ENDDO
         ENDDO
NOACC !$OMP END PARALLEL DO

      ENDIF

      CALL FFT3D_MU(N,CR,LDR,GRID,1)

      PROFILING_STOP('fftwav_useinv_mu')

    END SUBROUTINE FFTWAV_USEINV_MU


!****************** SUBROUTINE FFTEXT_MU *******************************
!
!***********************************************************************

    SUBROUTINE FFTEXT_MU(NPL,N,NINDPW,CR,LDR,C,LDC,GRID,LADD)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE mgrid

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: CR(LDR,*),C(LDC,*)
      INTEGER    :: NINDPW(NPL)
      INTEGER    :: NPL, N, LDR, LDC
      LOGICAL    :: LADD

      ! local variables
      INTEGER :: I,J

#ifdef _OPENACC
      LOGICAL :: ACC_ACTIVE
      ACC_ACTIVE=ACC_IS_PRESENT(NINDPW).AND.ACC_IS_PRESENT(CR).AND. &
          ACC_IS_PRESENT(C).AND.ACC_IS_PRESENT(GRID,1).AND.ACC_EXEC_ON
#endif

      PROFILING_START('fftext_mu')

      CALL FFT3D_MU(N,CR,LDR,GRID,-1)

      IF (GRID%REAL2CPLX) THEN
         IF (LADD) THEN
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CR,NINDPW,GRID) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW,GRID) PRIVATE(I,J)
            DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
               DO J=1,NPL
                  C(J,I)=C(J,I)+CR(NINDPW(J),I)*GRID%FFTSCA(J,1)
               ENDDO
            ENDDO
NOACC !$OMP END PARALLEL DO
         ELSE
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW,GRID) PRIVATE(I,J)
            DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
               DO J=1,NPL
                  C(J,I)=CR(NINDPW(J),I)*GRID%FFTSCA(J,1)
               ENDDO
            ENDDO
NOACC !$OMP END PARALLEL DO
         ENDIF
      ELSE
         IF (LADD) THEN
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CR,NINDPW) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW) PRIVATE(I,J)
            DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
               DO J=1,NPL
                  C(J,I)=C(J,I)+CR(NINDPW(J),I)
               ENDDO
            ENDDO
NOACC !$OMP END PARALLEL DO
         ELSE
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW) PRIVATE(I,J)
            DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
               DO J=1,NPL
                  C(J,I)=CR(NINDPW(J),I)
               ENDDO
            ENDDO
NOACC !$OMP END PARALLEL DO
         ENDIF
      ENDIF

      PROFILING_STOP('fftext_mu')

    END SUBROUTINE FFTEXT_MU


!****************** SUBROUTINE FFTEXT_USEINV_MU ************************
!
!***********************************************************************

    SUBROUTINE FFTEXT_USEINV_MU(NPL,N,NINDPW,FFTSCA,CR,LDR,C,LDC,GRID,LADD)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE mgrid

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: CR(LDR,*),C(LDC,*)
      REAL(q)    :: FFTSCA(NPL)
      INTEGER    :: NINDPW(NPL)
      INTEGER    :: NPL, N, LDR, LDC
      LOGICAL    :: LADD

      ! local variables
      INTEGER :: I,J

#ifdef _OPENACC
      LOGICAL :: ACC_ACTIVE
      ACC_ACTIVE=ACC_IS_PRESENT(NINDPW).AND.ACC_IS_PRESENT(FFTSCA).AND.ACC_IS_PRESENT(CR).AND. &
          ACC_IS_PRESENT(C).AND.ACC_IS_PRESENT(GRID,1).AND.ACC_EXEC_ON
#endif

      PROFILING_START('fftext_useinv_mu')

      CALL FFT3D_MU(N,CR,LDR,GRID,-1)

      IF (GRID%REAL2CPLX) THEN
         IF (LADD) THEN
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CR,NINDPW,GRID) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW,GRID) PRIVATE(I,J)
            DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
               DO J=1,NPL
                  C(J,I)=C(J,I)+CR(NINDPW(J),I)*GRID%FFTSCA(J,1)
               ENDDO
            ENDDO
NOACC !$OMP END PARALLEL DO
         ELSE
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,GRID) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW,GRID) PRIVATE(I,J)
            DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
               DO J=1,NPL
                  C(J,I)=CR(NINDPW(J),I)*GRID%FFTSCA(J,1)
               ENDDO
            ENDDO
NOACC !$OMP END PARALLEL DO
         ENDIF
      ELSE
         IF (LADD) THEN
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(C,CR,NINDPW,FFTSCA) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW,FFTSCA) PRIVATE(I,J)
            DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
               DO J=1,NPL
                  C(J,I)=C(J,I)+CR(NINDPW(J),I)*FFTSCA(J)
               ENDDO
            ENDDO
NOACC !$OMP END PARALLEL DO
         ELSE
!$ACC PARALLEL LOOP PRESENT(C,CR,NINDPW,FFTSCA) IF (ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(N,NPL,C,CR,NINDPW,FFTSCA) PRIVATE(I,J)
            DO I=1,N
!DIR$ IVDEP
!OCL NOVREC
               DO J=1,NPL
                  C(J,I)=CR(NINDPW(J),I)*FFTSCA(J)
               ENDDO
            ENDDO
NOACC !$OMP END PARALLEL DO
         ENDIF
      ENDIF

      PROFILING_STOP('fftext_useinv_mu')

    END SUBROUTINE FFTEXT_USEINV_MU


!****************** SUBROUTINE FFT3D_MU ********************************
!
!***********************************************************************

    SUBROUTINE FFT3D_MU(N,C,LDC,GRID,ISN)
      USE prec
      USE mgrid

      IMPLICIT NONE

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: C(LDC,*)
      INTEGER :: N,LDC,ISN

      ! local variables
      INTEGER :: NX,NY,NZ

      NX=GRID%NGPTAR(1)
      NY=GRID%NGPTAR(2)
      NZ=GRID%NGPTAR(3)

      IF (.NOT. GRID%REAL2CPLX) THEN
!-------------------------------------------------------------------------------
!  complex to complex FFTs
!-------------------------------------------------------------------------------
         IF (.NOT. (NX==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("FFT3D_MU: grid dimensions not properly set (1) " // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // &
               str(NX) // " " // str(NY) // " " // str(NZ) // "\n " // &
               str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " // str(GRID%NGZ_rd), &
               __FILE__, __LINE__)
         ENDIF

         IF  (.NOT. GRID%LREAL ) THEN

            CALL FFTBAS_MU(N,C,LDC,GRID,ISN)

         ELSE

            CALL vtutor%bug("FFT3D_MU: batched ffts for REAL2CPLX = " // str(GRID%REAL2CPLX) // &
               " and LREAL = " // str(GRID%LREAL) //" not supported yet",__FILE__,__LINE__)

         ENDIF
      ELSE
!-------------------------------------------------------------------------------
!  real to complex FFTs
!-------------------------------------------------------------------------------
         IF (.NOT. (NX/2+1==GRID%NGX_rd .AND. NY==GRID%NGY_rd .AND. NZ==GRID%NGZ_rd) ) THEN
            CALL vtutor%bug("FFT3D_MU: grid dimensions not properly set (2)" // &
               str(GRID%LREAL) // " " // str(GRID%REAL2CPLX) // "\n " // &
               str(NX) // " " // str(NY) // " " // str(NZ) // "\n " // &
               str(GRID%NGX_rd) // " " // str(GRID%NGY_rd) // " " //  str(GRID%NGZ_rd), &
                __FILE__, __LINE__)
         ENDIF

         IF (.NOT. GRID%LREAL) THEN
!
! the data in real space is stored as an array of complex numbers (with Im(c)=0)
!
            IF (ISN==1) THEN
               ! q -> r FFT
               CALL FFTBRC_MU(N,C,LDC,GRID,ISN)
               ! contract columns along x (go from stride NX+2 to stride NX)
               ! and store real valued transform as a array of complex numbers
               ! (with Im(c) = 0)
               CALL RESTRIDE_Q2R_CMPLX_MU(N,NY*NZ,NX,C,LDC)
            ELSE
               ! expand columns along x (go from stride NX to stride NX+2)
               ! and store array of complex numbers with Im(c)=0 as array 
               ! of reals
               CALL RESTRIDE_R2Q_CMPLX_MU(N,NY*NZ,NX,C,LDC)
               ! r -> q FFT
               CALL FFTBRC_MU(N,C,LDC,GRID,ISN)
            ENDIF

         ELSE
!
! the data in real space is stored as an array of real numbers
!
            IF (ISN==1) THEN
               ! q -> r FFT
               CALL FFTBRC_MU(N,C,LDC,GRID,ISN)
               ! contract columns along x (go from stride NX+2 to stride NX)
               CALL RESTRIDE_Q2R_MU(N,NY*NZ,NX,C,LDC)
            ELSE
               ! expand columns along x (go from stride NX to stride NX+2)
               CALL RESTRIDE_R2Q_MU(N,NY*NZ,NX,C,LDC)
               ! q -> r FFT
               CALL FFTBRC_MU(N,C,LDC,GRID,ISN)
            ENDIF

         ENDIF
      ENDIF

    END SUBROUTINE FFT3D_MU


!****************** SUBROUTINE RESTRIDE_Q2R_CMPLX_MU *******************
!
! The output of an inplace complex-to-real FFT on a grid of NX*NY*NZ
! points is an ( NX+2,NY,NZ ) array of real numbers (where the first
! dimension is padded with zeros.
!
! RESTRIDE_Q2R_CMPLX_MU does two things:
!
! ) The padding is removed, i.e., the data is re-stored as an
!   ( NX,NY,NZ ) array.
! ) The data is stored as an array of complex numbers (with Im(c)=0).
!
!***********************************************************************

    SUBROUTINE RESTRIDE_Q2R_CMPLX_MU(N,NCOL,NX,C,LDC)
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE prec
      IMPLICIT NONE
      INTEGER :: N,NCOL,NX,LDC
      REAL(q) :: C(2*LDC,*)
      ! local variables
      INTEGER :: ICOL,NDEST,NSRC,I,J

#ifdef _OPENACC
      REAL(q), ALLOCATABLE :: CTMP(:,:)

      ALLOCATE(CTMP(NCOL*(NX+2),N))
!$ACC ENTER DATA CREATE(CTMP) __IF_ASYNC__
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CTMP,C) __IF_ASYNC__
      DO J=1,N
         DO I=1,(NCOL-1)*(NX+2)+NX
            CTMP(I,J)=C(I,J)
         ENDDO
      ENDDO
#endif

!$ACC PARALLEL LOOP COLLAPSE(force:3) PRESENT(C,CTMP) PRIVATE(NDEST,NSRC) __IF_ASYNC__
NOACC !$OMP PARALLEL DO PRIVATE(J,ICOL,NDEST,NSRC,I)
      DO J=1,N
         DO ICOL=NCOL-1,0,-1
            NDEST=ICOL* NX*2
            NSRC =ICOL*(NX+2)
            DO I=NX,1,-1
NOACC          C(NDEST+I*2-1,J)=C   (NSRC+I,J)
DOACC          C(NDEST+I*2-1,J)=CTMP(NSRC+I,J)
               C(NDEST+I*2  ,J)=0
            ENDDO
         ENDDO
      ENDDO
NOACC !$OMP END PARALLEL DO

#ifdef _OPENACC
!$ACC EXIT DATA DELETE(CTMP) __IF_ASYNC__
      DEALLOCATE(CTMP)
#endif
    END SUBROUTINE RESTRIDE_Q2R_CMPLX_MU


!****************** SUBROUTINE RESTRIDE_R2Q_CMPLX_MU *******************
!
! RESTRIDE_R2Q_CMPLX_MU adresses two peculiarities:

! ) The input of an inplace real-to-complex FFT on a grid of NX*NY*NZ
!   points is an ( NX+2,NY,NZ ) array of real numbers (where the first
!   dimension is padded with zeros.
!
! ) In many cases the data (even though it is real) in real space is
!   stored as an array of complex numbers (with Im(c)=0)
!
! RESTRIDE_R2Q_CMPLX_MU does two things:
!
! ) Padding is added, i.e., the data is re-stored as an ( NX+2,NY,NZ )
!   array.
! ) The data is stored as an array of real numbers.
!
!***********************************************************************

    SUBROUTINE RESTRIDE_R2Q_CMPLX_MU(N,NCOL,NX,C,LDC)
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE prec
      IMPLICIT NONE
      INTEGER :: N,NCOL,NX,LDC
      REAL(q) :: C(2*LDC,*)
      ! local variables
      INTEGER :: ICOL,NDEST,NSRC,I,J

#ifdef _OPENACC
      REAL(q), ALLOCATABLE :: CTMP(:,:)

      ALLOCATE(CTMP(2*NCOL*NX,N))
!$ACC ENTER DATA CREATE(CTMP) __IF_ASYNC__
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CTMP,C) __IF_ASYNC__
      DO J=1,N
         DO I=1,2*NCOL*NX
            CTMP(I,J)=C(I,J)
         ENDDO
      ENDDO
#endif

!$ACC PARALLEL LOOP COLLAPSE(force:3) PRESENT(C,CTMP) PRIVATE(NSRC,NDEST) __IF_ASYNC__
NOACC !$OMP PARALLEL DO PRIVATE(J,ICOL,NSRC,NDEST,I)
      DO J=1,N
         DO ICOL=0,NCOL-1
            NSRC =ICOL* NX*2
            NDEST=ICOL*(NX+2)
            DO I=1,NX
NOACC          C(NDEST+I,J)=C   (NSRC+I*2-1,J)
DOACC          C(NDEST+I,J)=CTMP(NSRC+I*2-1,J)
            ENDDO
         ENDDO
      ENDDO
NOACC !$OMP END PARALLEL DO

#ifdef _OPENACC
!$ACC EXIT DATA DELETE(CTMP) __IF_ASYNC__
      DEALLOCATE(CTMP)
#endif
    END SUBROUTINE RESTRIDE_R2Q_CMPLX_MU


!****************** SUBROUTINE RESTRIDE_Q2R_MU *************************
!
! The output of an inplace complex-to-real FFT on a grid of NX*NY*NZ
! points is an ( NX+2,NY,NZ ) array of real numbers (where the first
! dimension is padded with zeros.
!
! RESTRIDE_Q2R_MU removes the padding, i.e., the data is re-stored as
! an ( NX,NY,NZ ) array.
!
!***********************************************************************

    SUBROUTINE RESTRIDE_Q2R_MU(N,NCOL,NX,C,LDC)
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE prec
      IMPLICIT NONE
      INTEGER :: N,NCOL,NX,LDC
      REAL(q) :: C(2*LDC,*)
      ! local variables
      INTEGER :: ICOL,NDEST,NSRC,I,J

#ifdef _OPENACC
      REAL(q), ALLOCATABLE :: CTMP(:,:)

      ALLOCATE(CTMP(NCOL*(NX+2),N))
!$ACC ENTER DATA CREATE(CTMP) __IF_ASYNC__
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CTMP,C) __IF_ASYNC__
      DO J=1,N
         DO I=1,(NCOL-1)*(NX+2)+NX
            CTMP(I,J)=C(I,J)
         ENDDO
      ENDDO
#endif

!$ACC PARALLEL LOOP COLLAPSE(force:3) PRESENT(C,CTMP) PRIVATE(NDEST,NSRC) __IF_ASYNC__
NOACC !$OMP PARALLEL DO PRIVATE(J,ICOL,NDEST,NSRC,I)
      DO J=1,N
         DO ICOL=1,NCOL-1
            NDEST=ICOL* NX
            NSRC =ICOL*(NX+2)
            DO I=1,NX
NOACC          C(NDEST+I,J)=C   (NSRC+I,J)
DOACC          C(NDEST+I,J)=CTMP(NSRC+I,J)
            ENDDO
         ENDDO
      ENDDO
NOACC !$OMP END PARALLEL DO

#ifdef _OPENACC
!$ACC EXIT DATA DELETE(CTMP) __IF_ASYNC__
      DEALLOCATE(CTMP)
#endif
    END SUBROUTINE RESTRIDE_Q2R_MU


!****************** SUBROUTINE RESTRIDE_R2Q_MU *************************
!
! The input of an inplace real-to-complex FFT on a grid of NX*NY*NZ
! points is an ( NX+2,NY,NZ ) array of real numbers (where the first
! dimension is padded with zeros.
!
! RESTRIDE_R2Q_MU add this padding, i.e., the data is re-stored as an
! ( NX+2,NY,NZ ) array.
!
!***********************************************************************

    SUBROUTINE RESTRIDE_R2Q_MU(N,NCOL,NX,C,LDC)
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE prec
      IMPLICIT NONE
      INTEGER :: N,NCOL,NX,LDC
      REAL(q) :: C(2*LDC,*)
      ! local variables
      INTEGER :: ICOL,NDEST,NSRC,I,J

#ifdef _OPENACC
      REAL(q), ALLOCATABLE :: CTMP(:,:)

      ALLOCATE(CTMP(NCOL*NX,N))
!$ACC ENTER DATA CREATE(CTMP) __IF_ASYNC__
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CTMP,C) __IF_ASYNC__
      DO J=1,N
         DO I=1,NCOL*NX
            CTMP(I,J)=C(I,J)
         ENDDO
      ENDDO
#endif

!$ACC PARALLEL LOOP COLLAPSE(force:3) PRESENT(C,CTMP) PRIVATE(NSRC,NDEST) __IF_ASYNC__
NOACC !$OMP PARALLEL DO PRIVATE(J,ICOL,NSRC,NDEST,I)
      DO J=1,N
         DO ICOL=NCOL-1,1,-1
            NSRC =ICOL* NX
            NDEST=ICOL*(NX+2)
            DO I=NX,1,-1
NOACC          C(NDEST+I,J)=C   (NSRC+I,J)
DOACC          C(NDEST+I,J)=CTMP(NSRC+I,J)
            ENDDO
         ENDDO
      ENDDO
NOACC !$OMP END PARALLEL DO

#ifdef _OPENACC
!$ACC EXIT DATA DELETE(CTMP) __IF_ASYNC__
      DEALLOCATE(CTMP)
#endif
    END SUBROUTINE RESTRIDE_R2Q_MU


!****************** SUBROUTINE FFTBAS_MU *******************************
!
!***********************************************************************

    SUBROUTINE FFTBAS_MU(N,C,LDC,GRID,ISN)
#ifdef _OPENACC
      USE mopenacc
      USE mcufft
#endif
      USE prec
      USE mgrid

      USE iso_c_binding

      IMPLICIT NONE

#ifndef __NEC__
      INCLUDE 'fftw3.f'
#else // __NEC__
      include 'aslfftw3.f'
#endif // __NEC__

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: C(*)
      INTEGER :: N,LDC,ISN

      ! local variables
      TYPE(c_ptr) :: PLAN

      INTEGER(c_int) :: DIRECTION,EFFORT

!$    INTEGER :: THREADS
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS

#ifdef _OPENACC
!=======================================================================
! Call cuFFT
!=======================================================================
      TYPE(acc_cufft_plan) :: ACC_PLAN
      INTEGER :: IERR

      IF (ACC_IS_PRESENT(C).AND.ACC_EXEC_ON) THEN
!
! Find corresponding plan
!
         IERR=ACC_CUFFT_GETPLAN(GRID,INT(ACC_ASYNC_Q),N,ACC_PLAN)

! Plan not found? Exit in error
         IF (IERR/=0) &
            CALL vtutor%bug("FFTBAS_MU: no cuFFT plan found for ACC_ASYNC_Q= " // &
               str(INT(ACC_ASYNC_Q)) // " and batch size N= " // str(N),__FILE__,__LINE__)

!
! Execute plan
!
#if PGI_BETWEEN_20_9_AND_XX_X
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_PLAN%LFORCE_ACC_WAIT)
#endif
!$ACC HOST_DATA USE_DEVICE(C)
         IF (ISN.LE.0) THEN
            IERR=cufftExecZ2Z(ACC_PLAN%CUFFT_3D_Z2Z,C,C,CUFFT_FORWARD)
         ELSE
            IERR=cufftExecZ2Z(ACC_PLAN%CUFFT_3D_Z2Z,C,C,CUFFT_INVERSE)
         ENDIF
!$ACC END HOST_DATA
         RETURN
      ENDIF
#endif

!=======================================================================
! Call FFTW
!=======================================================================
      IF (ISN.LE.0) THEN
         DIRECTION=FFTW_FORWARD
      ELSE
         DIRECTION=FFTW_BACKWARD
      ENDIF

      EFFORT=FFTW_ESTIMATE
#ifdef use_fftw_plan_effort
      IF (FFTW_PLAN_EFFORT==1) EFFORT=FFTW_MEASURE
#endif

!
! How many OpenMP threads are available?
!
!$OMP PARALLEL SHARED(THREADS)
!$OMP MASTER
!$    THREADS=OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

!
! Create batched 3d fft plan
!
!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!$    CALL dfftw_plan_with_nthreads(THREADS)

      CALL dfftw_plan_many_dft(PLAN, 3, GRID%NGPTAR, N, &
                               C(1), GRID%NGPTAR, 1, LDC, &
                               C(1), GRID%NGPTAR, 1, LDC, &
                               DIRECTION, EFFORT)
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

!
! Execute plan
!
      CALL dfftw_execute_dft(PLAN,C(1),C(1))

!
! Destroy plan
!
!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
      CALL dfftw_destroy_plan(PLAN)
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

      RETURN
    END SUBROUTINE FFTBAS_MU


!****************** SUBROUTINE FFTBRC_MU *******************************
!
!***********************************************************************

    SUBROUTINE FFTBRC_MU(N,C,LDC,GRID,ISN)
#ifdef _OPENACC
      USE mopenacc
      USE mcufft
#endif
      USE prec
      USE mgrid

      USE iso_c_binding

      IMPLICIT NONE

#ifndef __NEC__
      INCLUDE 'fftw3.f'
#else // __NEC__
      include 'aslfftw3.f'
#endif // __NEC__

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: C(*)
      INTEGER :: N,LDC,ISN

      ! local variables
      TYPE(c_ptr) :: PLAN

      INTEGER(c_int) :: EFFORT

      REAL(q), POINTER :: R(:)

!$    INTEGER :: THREADS
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS

#ifdef _OPENACC
!=======================================================================
! Call cuFFT
!=======================================================================
      TYPE(acc_cufft_plan) :: ACC_PLAN
      INTEGER :: IERR

      IF (ACC_IS_PRESENT(C).AND.ACC_EXEC_ON) THEN
!
! Find corresponding plan
!
         IERR=ACC_CUFFT_GETPLAN(GRID,INT(ACC_ASYNC_Q),N,ACC_PLAN)

! Plan not found? Exit in error
         IF (IERR/=0) &
            CALL vtutor%bug("FFTBAS_MU: no cuFFT plan found for ACC_ASYNC_Q= " // &
               str(INT(ACC_ASYNC_Q)) // " and batch size N= " // str(N),__FILE__,__LINE__)

!
! Connect a real pointer to the data buffer as well
        CALL c_f_pointer(c_loc(C),R,[2*LDC*N])

!
! Execute plan
!

!$ACC HOST_DATA USE_DEVICE(C,R)
         IF (ISN.LE.0) THEN
            IERR=cufftExecD2Z(ACC_PLAN%CUFFT_3D_D2Z,R,C)
         ELSE
            IERR=cufftExecZ2D(ACC_PLAN%CUFFT_3D_Z2D,C,R)
         ENDIF
!$ACC END HOST_DATA

! cleanup
         NULLIFY(R)

         RETURN
      ENDIF
#endif

!=======================================================================
! Call FFTW
!=======================================================================
      EFFORT=FFTW_ESTIMATE
#ifdef use_fftw_plan_effort
      IF (FFTW_PLAN_EFFORT==1) EFFORT=FFTW_MEASURE
#endif

!
! How many OpenMP threads are available?
!
!$OMP PARALLEL SHARED(THREADS)
!$OMP MASTER
!$    THREADS=OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

!
! Create and execute batched 3d r2c or c2r fft plan
!
!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!$    CALL dfftw_plan_with_nthreads(THREADS)

      IF (ISN.LE.0) THEN
         CALL dfftw_plan_many_dft_r2c(PLAN, 3, GRID%NGPTAR, N, &
                                      C(1), [GRID%NGX+2 , GRID%NGY   , GRID%NGZ   ], 1, 2* LDC, &
                                      C(1), [GRID%NGX_rd, GRID%NGY_rd, GRID%NGZ_rd], 1,    LDC, &
                                      EFFORT)
         CALL dfftw_execute_dft_r2c(PLAN,C(1),C(1))
      ELSE
         CALL dfftw_plan_many_dft_c2r(PLAN, 3, GRID%NGPTAR, N, &
                                      C(1), [GRID%NGX_rd, GRID%NGY_rd, GRID%NGZ_rd], 1,    LDC, &
                                      C(1), [GRID%NGX+2 , GRID%NGY   , GRID%NGZ   ], 1, 2* LDC, &
                                      EFFORT)
         CALL dfftw_execute_dft_c2r(PLAN,C(1),C(1))
      ENDIF
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

!
! Destroy plan
!
!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
      CALL dfftw_destroy_plan(PLAN)
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

      RETURN
    END SUBROUTINE FFTBRC_MU


!****************** SUBROUTINE FFTMAKEPLAN_MU **************************
!
!***********************************************************************

    SUBROUTINE FFTMAKEPLAN_MU(N,C,LDC,GRID)
      USE prec
      USE mgrid
      USE openmp, ONLY : omp_dfftw_init_threads

      USE iso_c_binding

      IMPLICIT NONE

#ifndef __NEC__
      INCLUDE 'fftw3.f'
#else // __NEC__
      include 'aslfftw3.f'
#endif // __NEC__

      TYPE(grid_3d) :: GRID

      COMPLEX(q) :: C(*)
      INTEGER :: N,LDC

      ! local variables
      TYPE(c_ptr) :: PLAN
!!      INTEGER, POINTER :: INULLPTR => NULL()

      INTEGER(c_int) :: EFFORT

!$    INTEGER :: THREADS, ierror
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS

      PROFILING_START('fftmakeplan_mu')

      EFFORT=FFTW_MEASURE
      IF (FFTW_PLAN_EFFORT==0) EFFORT=FFTW_ESTIMATE
      IF (FFTW_PLAN_EFFORT==1) EFFORT=FFTW_MEASURE

!$    IF (omp_dfftw_init_threads) THEN
!$       CALL dfftw_init_threads(ierror)
!$       omp_dfftw_init_threads=.FALSE.
!$    ENDIF

!$OMP PARALLEL SHARED(THREADS)
!$OMP MASTER
!$    THREADS=OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

!$OMP CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)
!$    CALL dfftw_plan_with_nthreads(THREADS)

      IF (.NOT.GRID%REAL2CPLX) THEN
         CALL dfftw_plan_many_dft(PLAN, 3, GRID%NGPTAR, N, &
                                  C, GRID%NGPTAR, 1, LDC, &
                                  C, GRID%NGPTAR, 1, LDC, &
                                  FFTW_FORWARD, EFFORT)
         CALL dfftw_destroy_plan(PLAN)

         CALL dfftw_plan_many_dft(PLAN, 3, GRID%NGPTAR, N, &
                                  C, GRID%NGPTAR, 1, LDC, &
                                  C, GRID%NGPTAR, 1, LDC, &
                                  FFTW_BACKWARD, EFFORT)
         CALL dfftw_destroy_plan(PLAN)
      ELSE
         CALL dfftw_plan_many_dft_r2c(PLAN, 3, GRID%NGPTAR, N, &
                                      C(1), [GRID%NGX+2 , GRID%NGY   , GRID%NGZ   ], 1, 2* LDC, &
                                      C(1), [GRID%NGX_rd, GRID%NGY_rd, GRID%NGZ_rd], 1,    LDC, &
                                      EFFORT)
         CALL dfftw_destroy_plan(PLAN)

         CALL dfftw_plan_many_dft_c2r(PLAN, 3, GRID%NGPTAR, N, &
                                      C(1), [GRID%NGX_rd, GRID%NGY_rd, GRID%NGZ_rd], 1,    LDC, &
                                      C(1), [GRID%NGX+2 , GRID%NGY   , GRID%NGZ   ], 1, 2* LDC, &
                                      EFFORT)
         CALL dfftw_destroy_plan(PLAN)
      ENDIF
!$OMP END CRITICAL (VASP_FFT_PLAN_CREATE_DESTROY)

      PROFILING_STOP('fftmakeplan_mu')

      RETURN
    END SUBROUTINE FFTMAKEPLAN_MU
