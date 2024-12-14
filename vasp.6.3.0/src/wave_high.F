#include "symbol.inc"
#define boundary_check
#define MPI_barrier_after_bcast
!***********************************************************************
!
! this module implements high level routines to operate on wavefunctions
! care is taken to avoid indexing pointer arrays since
! this incures performance penalties on all Intel compilers 
! (but particularly on the efc compiler)
! (constructs such as W1%CW(I,..) are avoided)
! instead low level F77 routines are used, which are implemented
! at the end of the routine
!
!***********************************************************************

MODULE wave_high
  USE prec
  USE wave
  USE wave_mpi
!***********************************************************************
!
! interfaces for functions that accept only F77 style arrays
! these F77 routines would cause a problem if the first element
! of a pointer array is passed to the F77 routines
!
!***********************************************************************

  INTERFACE
     SUBROUTINE ARRAY_TO_W1( W1, C, CPROJ)
       USE wave
       TYPE (wavefun1)    W1
       COMPLEX(q):: C
       GDEF, OPTIONAL :: CPROJ
     END SUBROUTINE ARRAY_TO_W1
  END INTERFACE
  
  INTERFACE
     SUBROUTINE W1_TO_ARRAY( W1, C, CPROJ)
       USE wave
       TYPE (wavefun1)    W1
       COMPLEX(q):: C
       GDEF, OPTIONAL :: CPROJ
     END SUBROUTINE W1_TO_ARRAY
  END INTERFACE

  INTERFACE
     SUBROUTINE ECCP_NL(LMDIM,LMMAXC,CDIJ,CPROJ1,CPROJ2,CNL)
!$ACC ROUTINE VECTOR
       USE prec
       GDEF      CNL
       INTEGER LMDIM, LMMAXC
       OVERLAP CDIJ
       GDEF CPROJ1,CPROJ2
     END SUBROUTINE ECCP_NL
  END INTERFACE

  INTERFACE
     SUBROUTINE OVERL(WDES1, LOVERL, LMDIM, CQIJ, CPROF, CRESUL)
       USE wave
       TYPE (wavedes1) WDES1
       LOGICAL LOVERL
       INTEGER LMDIM
       OVERLAP CQIJ,CDIJ
       GDEF CRESUL,CPROF
     END SUBROUTINE OVERL
  END INTERFACE

  INTERFACE
     SUBROUTINE OVERL1(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF,CRESUL)
       USE wave
       TYPE (wavedes1) WDES1
       INTEGER LMDIM
       OVERLAP CQIJ,CDIJ
       REAL(q) :: EVALUE
       GDEF CRESUL,CPROF
     END SUBROUTINE OVERL1
  END INTERFACE

  INTERFACE
     SUBROUTINE OVERL1_C(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF,CRESUL)
       USE wave
       TYPE (wavedes1) WDES1
       INTEGER LMDIM
       OVERLAP CQIJ,CDIJ
       COMPLEX(q) EVALUE
       GDEF CRESUL,CPROF
     END SUBROUTINE OVERL1_C
  END INTERFACE

  INTERFACE
     SUBROUTINE OVERL1_CCDIJ(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF,CRESUL)
       USE wave
       TYPE (wavedes1) WDES1
       INTEGER LMDIM
       COMPLEX(q) CQIJ,CDIJ
       REAL(q) :: EVALUE
       GDEF CRESUL,CPROF
     END SUBROUTINE OVERL1_CCDIJ
  END INTERFACE


  INTERFACE ELEMENT
     MODULE PROCEDURE W1_FROM_W
     MODULE PROCEDURE W1_FROM_WA
     MODULE PROCEDURE W1_FROM_WA2
  END INTERFACE

  INTERFACE ELEMENTS
     MODULE PROCEDURE WA_FROM_W
     MODULE PROCEDURE WA_FROM_WA
     MODULE PROCEDURE WA_FROM_WA2
  END INTERFACE

  INTERFACE REDISTRIBUTE_PROJ
     MODULE PROCEDURE W_REDIS_PROJ
     MODULE PROCEDURE WA_REDIS_PROJ
     MODULE PROCEDURE W1_REDIS_PROJ
  END INTERFACE

  INTERFACE REDISTRIBUTE_PW
     MODULE PROCEDURE W_REDIS_PW
     MODULE PROCEDURE WA_REDIS_PW
     MODULE PROCEDURE W1_REDIS_PW
  END INTERFACE

  INTERFACE ASSIGNMENT
     MODULE PROCEDURE W1_COPY_REVERSE_ARG
  END INTERFACE
  
  CONTAINS

!***********************************************************************
!
! FFT of a wavefunction to real space
!      
!***********************************************************************

    SUBROUTINE FFTWAV_W1( W1)
      IMPLICIT NONE
      INTEGER ISPINOR
      TYPE (wavefun1)    W1

      IF (W1%WDES1%NK==0) THEN
         CALL vtutor%bug("internal error in FFTWAV_W1: NK is set to zero: WDES not set up properly", __FILE__, __LINE__)
      ENDIF
      IF (W1%WDES1%LUSEINV) THEN
         IF (.NOT. ASSOCIATED(W1%WDES1%FFTSCA) .OR. .NOT. ASSOCIATED(W1%WDES1%NINDPW_INV)) THEN
            CALL vtutor%bug("internal error in FFTWAV_W1: FFTSCA is not associated and LUSEINV is &
               &set " // str(ASSOCIATED(W1%WDES1%FFTSCA)) // " " // str(ASSOCIATED(W1%WDES1%NINDPW_INV)), &
               __FILE__, __LINE__)
         ENDIF
         DO ISPINOR=0,W1%WDES1%NRSPINORS-1
            CALL FFTWAV_USEINV(W1%WDES1%NGVECTOR, W1%WDES1%NINDPW(1), W1%WDES1%NINDPW_INV(1), W1%WDES1%FFTSCA(1,2), & 
                 W1%CR(1+ISPINOR*W1%WDES1%GRID%MPLWV), & 
                 W1%CW(1+ISPINOR*W1%WDES1%NGVECTOR),W1%WDES1%GRID)
         ENDDO
      ELSE
         DO ISPINOR=0,W1%WDES1%NRSPINORS-1
            CALL FFTWAV(W1%WDES1%NGVECTOR, W1%WDES1%NINDPW(1), & 
                 W1%CR(1+ISPINOR*W1%WDES1%GRID%MPLWV), & 
                 W1%CW(1+ISPINOR*W1%WDES1%NGVECTOR),W1%WDES1%GRID)
         ENDDO
      ENDIF
    END SUBROUTINE FFTWAV_W1

!
! FFT of a wavefunction to real space
! similar version as above, however, receives WDES1 and explicit "WORK" arrays CR and CW
! this routine is new in vasp.6 and should consistently replace
! FFTWAV
! 

    SUBROUTINE FFTWAV_WDES1( WDES1, CR, CW)
      IMPLICIT NONE
      TYPE (wavedes1)    WDES1          ! descriptor for one k-point
      COMPLEX(q) :: CR(WDES1%GRID%RL%NP)  ! input:  orbitals in real space
      COMPLEX(q) :: CW(WDES1%NRPLWV)      ! output: orbitals in reciprocal space

      IF (WDES1%NK==0) THEN
         CALL vtutor%bug("internal error in FFTWAV_WDES1: NK is set to zero: WDES not set up properly", __FILE__, __LINE__)
      ENDIF
      IF (WDES1%LUSEINV) THEN
         IF (.NOT. ASSOCIATED(WDES1%FFTSCA) .OR. .NOT. ASSOCIATED(WDES1%NINDPW_INV)) THEN
            CALL vtutor%bug("internal error in FFTWAV_WDES1: FFTSCA is not associated and LUSEINV &
               &is set " // str(ASSOCIATED(WDES1%FFTSCA)) // " " // str(ASSOCIATED(WDES1%NINDPW_INV)), &
               __FILE__, __LINE__)
         ENDIF
         CALL FFTWAV_USEINV(WDES1%NGVECTOR, WDES1%NINDPW(1), WDES1%NINDPW_INV(1), WDES1%FFTSCA(1,2), & 
              CR(1),CW(1),WDES1%GRID)
      ELSE
         CALL FFTWAV(WDES1%NGVECTOR, WDES1%NINDPW(1), & 
              CR(1),CW(1),WDES1%GRID)
      ENDIF
    END SUBROUTINE FFTWAV_WDES1


!***********************************************************************
!
! FFT of a wavefunction to real space a WDES1 must be supplied
! as well as the input and output array
! this routine is new in vasp.6 and should consistently replace
! FFTEXT
!      
!***********************************************************************

    SUBROUTINE FFTEXT_WDES1(WDES1, CR, CW, LADD)
      IMPLICIT NONE
      INTEGER ISPINOR
      TYPE (wavedes1)    WDES1          ! descriptor for one k-point
      COMPLEX(q) :: CR(WDES1%GRID%RL%NP)! input:  orbitals in real space
      COMPLEX(q) :: CW(WDES1%NRPLWV)    ! output: orbitals in reciprocal space
      LOGICAL    :: LADD                ! add results to CW and does not clear CW beforehand

      IF (WDES1%NK==0) THEN
         CALL vtutor%bug("internal error in FFTWAV_WDES1: NK is set to zero: WDES not set up properly", __FILE__, __LINE__)
      ENDIF
      IF (WDES1%LUSEINV) THEN
         IF (.NOT. ASSOCIATED(WDES1%FFTSCA) .OR. .NOT. ASSOCIATED(WDES1%NINDPW_INV)) THEN
            CALL vtutor%bug("internal error in FFTEXT_WDES1: FFTSCA is not associated and LUSEINV &
               &is set " // str(ASSOCIATED(WDES1%FFTSCA)) // " " // str(ASSOCIATED(WDES1%NINDPW_INV)), &
               __FILE__, __LINE__)
         ENDIF
         CALL FFTEXT_USEINV(WDES1%NGVECTOR, WDES1%NINDPW(1), WDES1%FFTSCA(1,1), & 
              CR(1), CW(1), WDES1%GRID, LADD)
      ELSE
         CALL FFTEXT(WDES1%NGVECTOR, WDES1%NINDPW(1), & 
              CR(1), CW(1), WDES1%GRID, LADD)
      ENDIF
    END SUBROUTINE FFTEXT_WDES1


!***********************************************************************
!
! This routine accomplishes the same thing as the assignment:
!
!   W1_LHS = W1_RHS
!
! would, but allows for the use of "fast"-memory (in the sense that
! the fm_* arrays are excluded from the assignment).
!
!***********************************************************************

  SUBROUTINE W1_ASSIGN(W1_LHS,W1_RHS)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1_RHS
    TYPE (wavefun1) :: W1_LHS

    W1_LHS%WDES1  =>W1_RHS%WDES1
    W1_LHS%FERWE  = W1_RHS%FERWE
    W1_LHS%AUX    = W1_RHS%AUX
    W1_LHS%CELEN  = W1_RHS%CELEN
    W1_LHS%NB     = W1_RHS%NB
    W1_LHS%ISP    = W1_RHS%ISP  
    W1_LHS%LDO    = W1_RHS%LDO

    W1_LHS%CPTWFP =>W1_RHS%CPTWFP
    W1_LHS%CPROJ  =>W1_RHS%CPROJ    
    W1_LHS%CR     =>W1_RHS%CR

  END SUBROUTINE W1_ASSIGN


!***********************************************************************
!
! copy a W1 structure
!  W2 = W1  (W1 -> W2)
! the argument arrangement is similar to  DCOPY, ZCOPY in BLAS level 1
! for syntactic sugar the reverse operation is also present
!
!> @details @ref openmp :
!> BLAS1 calls were replaced by explicit OMP PARALLEL loops
!> (for some reason BLAS1 call thread really badly).
!
!***********************************************************************

  SUBROUTINE W1_COPY( W1, W2)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1    
    TYPE (wavefun1) :: W2
    INTEGER I

    PROFILING_START('w1_copy')

#if defined(_OPENACC) || ! defined(_OPENMP)
#if defined (use_fastbcopy)
    CALL FASTBCOPY8(W1%WDES1%NRPLWV*2 ,W1%CW(1), W2%CW(1))
#else
    CALL __ZCOPY__( W1%WDES1%NRPLWV, W1%CW(1), 1, W2%CW(1), 1)
#endif

    IF (W1%WDES1%LGAMMA) THEN
       CALL __DCOPY__( W1%WDES1%NPROD, W1%CPROJ(1), 1, W2%CPROJ(1), 1)
    ELSE
       CALL __ZCOPY__( W1%WDES1%NPROD, W1%CPROJ(1), 1, W2%CPROJ(1), 1)
    ENDIF

    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("W1_COPY: real space allocation is not correct",__FILE__,__LINE__)
       ENDIF
       IF (SIZE(W1%CR) /=SIZE(W2%CR)) THEN
          CALL vtutor%bug("W1_COPY: real space allocation is different " // &
             str(SIZE(W1%CR)) // " " // str(SIZE(W2%CR)),__FILE__,__LINE__)
       ENDIF
       CALL __ZCOPY__( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS, W1%CR(1), 1, W2%CR(1), 1)
    END IF
#else
!$OMP PARALLEL SHARED(W1,W2) PRIVATE(I)
!DIR$ IVDEP
!$OMP DO
    DO I=1,W1%WDES1%NRPLWV
       W2%CW(I)=W1%CW(I)
    ENDDO
!$OMP END DO
!DIR$ IVDEP
!$OMP DO
    DO I=1,W1%WDES1%NPROD
       W2%CPROJ(I)=W1%CPROJ(I)
    ENDDO
!$OMP END DO
    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("W1_COPY: real space allocation is not correct",__FILE__,__LINE__)
       ENDIF
       IF (SIZE(W1%CR) /=SIZE(W2%CR)) THEN
          CALL vtutor%bug("W1_COPY: real space allocation is different " // &
             str(SIZE(W1%CR)) // " " // str(SIZE(W2%CR)),__FILE__,__LINE__)
       ENDIF
!DIR$ IVDEP
!$OMP DO
       DO I=1,W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS
          W2%CR(I)=W1%CR(I)
       ENDDO
!$OMP END DO
    END IF
!$OMP END PARALLEL
#endif
    PROFILING_STOP('w1_copy')

  END SUBROUTINE W1_COPY

  SUBROUTINE W1_COPY_NOCR( W1, W2)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1    
    TYPE (wavefun1) :: W2

#if defined (use_fastbcopy)
    CALL FASTBCOPY8(W1%WDES1%NRPLWV*2 ,W1%CW(1), W2%CW(1))
#else
    CALL ZCOPY( W1%WDES1%NRPLWV, W1%CW(1), 1, W2%CW(1), 1)
#endif

    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( W1%WDES1%NPROD, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       CALL ZCOPY( W1%WDES1%NPROD, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF
  END SUBROUTINE W1_COPY_NOCR

  SUBROUTINE W1_COPY_CPROJ( W1, W2)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1
    TYPE (wavefun1) :: W2
    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( W1%WDES1%NPROD, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       CALL ZCOPY( W1%WDES1%NPROD, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF
  END SUBROUTINE W1_COPY_CPROJ

  SUBROUTINE W1_COPY_REVERSE_ARG( W2, W1 )
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1
    TYPE (wavefun1) :: W2

#if defined (use_fastbcopy)
    CALL FASTBCOPY8(W1%WDES1%NRPLWV*2 ,W1%CW(1), W2%CW(1))
#else
    CALL ZCOPY( W1%WDES1%NRPLWV, W1%CW(1), 1, W2%CW(1), 1)
#endif
    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( W1%WDES1%NPROD, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       CALL ZCOPY( W1%WDES1%NPROD, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF

    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is not correct", __FILE__, __LINE__)
       ENDIF
       IF (SIZE(W1%CR) /=SIZE(W2%CR)) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is different " // &
             str(SIZE(W1%CR)) // " " // str(SIZE(W2%CR)), __FILE__, __LINE__)
       ENDIF

       CALL ZCOPY( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS, W1%CR(1), 1, W2%CR(1), 1)
    END IF

    
  END SUBROUTINE W1_COPY_REVERSE_ARG


!***********************************************************************
!
! copy a WA structure
!  W2 = W1  (W1 -> W2)
! the argument arrangement is similar to  DCOPY, ZCOPY in BLAS level 1
! the redistributed wavefunctions pointers in the destination
! are also properly set
!
!***********************************************************************

  SUBROUTINE WA_COPY( W1, W2 )
    IMPLICIT NONE
    TYPE (wavefuna), INTENT(IN) :: W1    
    TYPE (wavefuna) :: W2

    IF (SIZE(W1%CW) /= SIZE(W2%CW)) THEN
       CALL vtutor%bug("internal error in WA_COPY: size mismatch in CW " // str(SIZE(W1%CW)) // &
          " " // str(SIZE(W2%CW)), __FILE__, __LINE__)
    ENDIF

    IF (SIZE(W1%CPROJ) /= SIZE(W2%CPROJ)) THEN
       CALL vtutor%bug("internal error in WA_COPY: size mismatch in CPROJ " // str(SIZE(W1%CPROJ)) &
          // " " // str(SIZE(W2%CPROJ)), __FILE__, __LINE__)
    ENDIF

    CALL ZCOPY( SIZE(W1%CW), W1%CW(1,1), 1, W2%CW(1,1), 1)

    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( SIZE(W1%CPROJ), W1%CPROJ(1,1), 1,  W2%CPROJ(1,1), 1)
    ELSE
       CALL ZCOPY( SIZE(W1%CPROJ), W1%CPROJ(1,1), 1,  W2%CPROJ(1,1), 1)
    ENDIF
    ! remember WDES1
    W2%WDES1 =>W1%WDES1
    ! one dimensional indexing assumed
    W2%FIRST_DIM=0
    ! remember spin index
    W2%ISP =W1%ISP
    ! set redistributed wavefunction indices
    IF (W2%WDES1%DO_REDIS) THEN
       CALL SET_WPOINTER(W2%CW_RED,    W2%WDES1%NRPLWV_RED, W2%WDES1%NB_TOT, W2%CW(1,1))
       CALL SET_GPOINTER(W2%CPROJ_RED, W2%WDES1%NPROD_RED,  W2%WDES1%NB_TOT, W2%CPROJ(1,1))
    ELSE
       W2%CW_RED=>W2%CPTWFP
       W2%CPROJ_RED=>W2%CPROJ
    ENDIF


  END SUBROUTINE WA_COPY


  SUBROUTINE WA_COPY_CPROJ( W1, W2 )
    IMPLICIT NONE
    TYPE (wavefuna), INTENT(IN) :: W1    
    TYPE (wavefuna) :: W2

    IF (SIZE(W1%CPROJ) /= SIZE(W2%CPROJ)) THEN
       CALL vtutor%bug("internal error in WA_COPY: size mismatch in CPROJ " // str(SIZE(W1%CPROJ)) &
          // " " // str(SIZE(W2%CPROJ)), __FILE__, __LINE__)
    ENDIF

    IF (W1%WDES1%LGAMMA) THEN
       CALL DCOPY( SIZE(W1%CPROJ), W1%CPROJ(1,1), 1,  W2%CPROJ(1,1), 1)
    ELSE
       CALL ZCOPY( SIZE(W1%CPROJ), W1%CPROJ(1,1), 1,  W2%CPROJ(1,1), 1)
      ENDIF
    ! remember WDES1
    W2%WDES1 =>W1%WDES1
    ! one dimensional indexing assumed
    W2%FIRST_DIM=0
    ! remember spin index
    W2%ISP =W1%ISP
    ! set redistributed wavefunction indices
    IF (W2%WDES1%DO_REDIS) THEN
       CALL SET_GPOINTER(W2%CPROJ_RED, W2%WDES1%NPROD_RED,  W2%WDES1%NB_TOT, W2%CPROJ(1,1))
    ELSE
       W2%CPROJ_RED=>W2%CPROJ
    ENDIF


  END SUBROUTINE WA_COPY_CPROJ


!***********************************************************************
!
! redistribute the wavefunction character for a W1 or WA array
!
!***********************************************************************

  SUBROUTINE W1_REDIS_PROJ( W1)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1

    IF (W1%WDES1%DO_REDIS) CALL REDIS_PROJ(W1%WDES1, 1, W1%CPROJ(1))
  END SUBROUTINE W1_REDIS_PROJ

  SUBROUTINE WA_REDIS_PROJ( WA)
    IMPLICIT NONE
    TYPE (wavefuna), INTENT(IN) :: WA

    IF (WA%WDES1%DO_REDIS) CALL REDIS_PROJ(WA%WDES1, SIZE(WA%CPROJ,2), WA%CPROJ(1,1))
  END SUBROUTINE WA_REDIS_PROJ

  SUBROUTINE W_REDIS_PROJ( W)
    IMPLICIT NONE
    TYPE (wavespin), INTENT(IN) :: W
    TYPE (wavedes1)    WDES1          ! descriptor for one k-point
    INTEGER :: K1, ISP

    IF (W%WDES%DO_REDIS) THEN
       DO K1=1,W%WDES%NKPTS
#ifdef MPI
          IF (MOD(K1-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
          CALL SETWDES(W%WDES,WDES1,K1)
          DO ISP=1,W%WDES%ISPIN
             CALL WA_REDIS_PROJ( ELEMENTS( W, WDES1, ISP))
          ENDDO
       ENDDO
    ENDIF
  END SUBROUTINE W_REDIS_PROJ



!***********************************************************************
!
! redistribute the plane wave coefficients for a W1 or WA array
!
!***********************************************************************

  SUBROUTINE W1_REDIS_PW( W1)
    IMPLICIT NONE
    TYPE (wavefun1), INTENT(IN) :: W1

    IF (W1%WDES1%DO_REDIS) CALL REDIS_PW(W1%WDES1, 1, W1%CW(1))
  END SUBROUTINE W1_REDIS_PW

  SUBROUTINE WA_REDIS_PW( WA)
    IMPLICIT NONE
    TYPE (wavefuna), INTENT(IN) :: WA

    IF (WA%WDES1%DO_REDIS) CALL REDIS_PW(WA%WDES1, SIZE(WA%CW,2), WA%CW(1,1))
  END SUBROUTINE WA_REDIS_PW


  SUBROUTINE W_REDIS_PW( W)
    IMPLICIT NONE
    TYPE (wavespin), INTENT(IN) :: W
    TYPE (wavedes1)    WDES1          ! descriptor for one k-point
    INTEGER :: K1, ISP

    IF (W%WDES%DO_REDIS) THEN
       DO K1=1,W%WDES%NKPTS
#ifdef MPI
          IF (MOD(K1-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
          CALL SETWDES(W%WDES,WDES1,K1)
          DO ISP=1,W%WDES%ISPIN
             CALL WA_REDIS_PW( ELEMENTS( W, WDES1, ISP))
          ENDDO
       ENDDO
    ENDIF
  END SUBROUTINE W_REDIS_PW


!***********************************************************************
!
! update of vector
!  W2 = W1*a + W2
!
!> @details @ref openmp :
!> BLAS1 calls were replaced by explicit OMP PARALLEL loops
!> (for some reason BLAS1 call thread really badly).
!
!***********************************************************************

  SUBROUTINE W1_DAXPY( W1, SCALE, W2)
    IMPLICIT NONE
    REAL(q) SCALE
    TYPE (wavefun1)    W1, W2
!$  INTEGER I

    PROFILING_START('w1_daxpy')

#if defined(_OPENACC) || ! defined( _OPENMP)
    CALL __DAXPY__( W1%WDES1%NPL*2, SCALE, W1%CW(1), 1, W2%CW(1), 1)

    IF (W1%WDES1%LGAMMA) THEN
       CALL __DAXPY__( W1%WDES1%NPRO, SCALE, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       CALL __DAXPY__( W1%WDES1%NPRO*2, SCALE, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF

    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is not correct", __FILE__, __LINE__)
       ENDIF

       ! real space wavefunction complex, SCALE real
       CALL __DAXPY__( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS*2, SCALE, W1%CR(1), 1, W2%CR(1), 1)
    END IF
#else
!$OMP PARALLEL SHARED(W1,W2,SCALE) PRIVATE(I)
!DIR$ IVDEP
!$OMP DO
    DO I=1,W1%WDES1%NPL
       W2%CW(I)=W2%CW(I)+SCALE*W1%CW(I)
    ENDDO
!$OMP END DO NOWAIT
!DIR$ IVDEP
!$OMP DO
    DO I=1,W1%WDES1%NPRO
       W2%CPROJ(I)=W2%CPROJ(I)+SCALE*W1%CPROJ(I)
    ENDDO
!$OMP END DO NOWAIT
    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is not correct", __FILE__, __LINE__)
       ENDIF
!DIR$ IVDEP
!$OMP DO
    DO I=1,W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS
       W2%CR(I)=W2%CR(I)+SCALE*W1%CR(I)
    ENDDO
!$OMP END DO NOWAIT
    END IF
!$OMP END PARALLEL
#endif

    PROFILING_STOP('w1_daxpy')

  END SUBROUTINE W1_DAXPY

!
!> @details @ref openmp :
!> BLAS1 calls were replaced by explicit OMP PARALLEL loops
!> (for some reason BLAS1 call thread really badly).
!
  SUBROUTINE W1_GAXPY( W1, SCALE, W2)
    IMPLICIT NONE
    GDEF SCALE
    TYPE (wavefun1)    W1, W2
!$  INTEGER I

    PROFILING_START('w1_gaxpy')

#if defined(_OPENACC) || ! defined( _OPENMP)
    IF (W1%WDES1%LGAMMA) THEN
       ! wavefunction complex, SCALE real
       CALL __DAXPY__( W1%WDES1%NPL*2, SCALE, W1%CW(1), 1, W2%CW(1), 1)
       CALL __DAXPY__( W1%WDES1%NPRO, SCALE, W1%CPROJ(1), 1,  W2%CPROJ(1), 1)
    ELSE
       ! wavefunction complex, SCALE complex
       CALL __ZAXPY__( W1%WDES1%NPL, SCALE, W1%CW(1), 1, W2%CW(1), 1)
       CALL __ZAXPY__( W1%WDES1%NPRO, SCALE, W1%CPROJ(1),  1, W2%CPROJ(1), 1)
    ENDIF

    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is not correct", __FILE__, __LINE__)
       ENDIF

       IF (W1%WDES1%LGAMMA) THEN
          ! real space wavefunction complex, SCALE real
          CALL __DAXPY__( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS*2, SCALE, W1%CR(1), 1, W2%CR(1), 1)
       ELSE
          ! real space wavefunction complex, SCALE complex
          CALL __ZAXPY__( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS, SCALE, W1%CR(1), 1, W2%CR(1), 1)
       ENDIF
             
    END IF
#else
!$OMP PARALLEL SHARED(W1,W2,SCALE) PRIVATE(I)
!DIR$ IVDEP
!$OMP DO
       DO I=1,W1%WDES1%NPL
          W2%CW(I)=W2%CW(I)+SCALE*W1%CW(I)
       ENDDO
!$OMP END DO NOWAIT
!DIR$ IVDEP
!$OMP DO
       DO I=1,W1%WDES1%NPRO
          W2%CPROJ(I)=W2%CPROJ(I)+SCALE*W1%CPROJ(I)
       ENDDO
!$OMP END DO NOWAIT
    IF (ASSOCIATED(W1%CR) .AND. ASSOCIATED(W2%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is not correct", __FILE__, __LINE__)
       ENDIF
!DIR$ IVDEP
!$OMP DO
       DO I=1,W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS
          W2%CR(I)=W2%CR(I)+SCALE*W1%CR(I)
       ENDDO
!$OMP END DO NOWAIT
    END IF
!$OMP END PARALLEL
#endif

    PROFILING_STOP('w1_gaxpy')

  END SUBROUTINE W1_GAXPY

!***********************************************************************
!
! update of vector
!  W1 = W1*a
!
!> @details @ref openmp :
!> BLAS1 calls were replaced by explicit OMP PARALLEL loops
!> (for some reason BLAS1 call thread really badly).
!
!***********************************************************************


  SUBROUTINE W1_DSCAL( W1, SCALE)
    IMPLICIT NONE
    REAL(q) SCALE
    TYPE (wavefun1)    W1
!$  INTEGER I

    PROFILING_START('w1_dscal')

    ! since this function can be used to zero out an array we operate
    ! on all elements (dimension) and not only on those that are
    ! actually used
#if defined(_OPENACC) || ! defined(_OPENMP)
    CALL __DSCAL__( W1%WDES1%NRPLWV*2, SCALE, W1%CW(1), 1)

    IF (W1%WDES1%LGAMMA) THEN
       CALL __DSCAL__( W1%WDES1%NPROD, SCALE, W1%CPROJ(1), 1)
    ELSE
       CALL __DSCAL__( W1%WDES1%NPROD*2, SCALE, W1%CPROJ(1),  1)
    ENDIF

    IF (ASSOCIATED(W1%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is not correct", __FILE__, __LINE__)
       ENDIF
       ! real space wavefunction complex, SCALE real
       CALL __DSCAL__( W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS*2, SCALE, W1%CR(1), 1)
    END IF
#else
!$OMP PARALLEL SHARED(W1,SCALE) PRIVATE(I)
!DIR$ IVDEP
!$OMP DO
     DO I=1,W1%WDES1%NRPLWV
        W1%CW(I)=SCALE*W1%CW(I)
     ENDDO 
!$OMP END DO NOWAIT
!DIR$ IVDEP
!$OMP DO
     DO I=1,W1%WDES1%NPROD
        W1%CPROJ(I)=SCALE*W1%CPROJ(I)
     ENDDO 
!$OMP END DO NOWAIT
    IF (ASSOCIATED(W1%CR)) THEN
       IF (SIZE(W1%CR) /=W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in W1_COPY: real space allocation is not correct", __FILE__, __LINE__)
       ENDIF
!DIR$ IVDEP
!$OMP DO
     DO I=1,W1%WDES1%GRID%MPLWV*W1%WDES1%NRSPINORS
        W1%CR(I)=SCALE*W1%CR(I)
     ENDDO 
!$OMP END DO NOWAIT
    END IF
!$OMP END PARALLEL
#endif

    PROFILING_STOP('w1_dscal')

  END SUBROUTINE W1_DSCAL


!***********************************************************************
!
! calculate the dot product between two wavefunctions
!  C=   W1^* x W2
! this is a substitue for the routine CINDPROD but mind
! the W1 and W2 are interchanged
!
!> @details @ref openmp :
!> a ZDOTC call is replaced by an explicit OMP PARALLEL loop
!> (for some reason BLAS1 call thread really badly).
!> The nested loop over \"types\" + \"ions-of-type\" is replaced by
!> a loop over \"all ions\", that is distributed over all available
!> threads.
!
!***********************************************************************

  FUNCTION W1_DOT( W1, W2, CQIJ) RESULT (C)
#ifdef _OPENACC
    USE mopenacc
#endif
    IMPLICIT NONE
    TYPE (wavefun1)   :: W1, W2
    GDEF              :: C
    OVERLAP, OPTIONAL :: CQIJ(:,:,:,:) ! optional overlap operator
  ! local    
    GDEF    :: CNL, CTMP
    INTEGER :: LMDIM, NPRO, NPRO_, ISPINOR, ISPINOR_, LMMAXC, NT, NI

    REAL(q),    EXTERNAL :: __DDOT__
    COMPLEX(q), EXTERNAL :: __ZDOTC__

!$  INTEGER :: I
!$ACC ROUTINE(ECCP_NL) VECTOR

    PROFILING_START('w1_dot')

#if defined(_OPENACC) || ! defined(_OPENMP)
    IF (W1%WDES1%LGAMMA) THEN
       C=__DDOT__(2*W1%WDES1%NPL,W1%CW(1),1,W2%CW(1),1)
    ELSE
       C=__ZDOTC__( W1%WDES1%NPL,W1%CW(1),1,W2%CW(1),1)
    ENDIF
#else
    C=0
!$OMP PARALLEL DO PRIVATE(I) SHARED(W1,W2) REDUCTION(+:C)
    DO I=1,W1%WDES1%NPL
       C=C+CONJG(W1%CW(I))*W2%CW(I)
    ENDDO
!$OMP END PARALLEL DO
#endif
    IF (PRESENT(CQIJ) .AND. W1%WDES1%LOVERL .AND.  W1%WDES1%NPROD>0 ) THEN
       CNL=0
       LMDIM=SIZE(CQIJ,1)
#ifdef _OPENACC
!$ACC ENTER DATA COPYIN(W1,W2) __IF_ASYNC__
!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(CQIJ,W1,W1%WDES1,W1%CPROJ,W2,W2%CPROJ) REDUCTION(+:CNL) &
!$ACC PRIVATE(NT,LMMAXC,NPRO,NPRO_,CTMP) __IF_ASYNC__
#else
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP SHARED(LMDIM,W1,W2,CQIJ) &
!$OMP PRIVATE(ISPINOR,ISPINOR_,NI,NT,LMMAXC,NPRO,NPRO_,CTMP) &
!$OMP REDUCTION(+:CNL)
#endif
       spinor: DO ISPINOR=0,W1%WDES1%NRSPINORS-1
       DO ISPINOR_=0,W1%WDES1%NRSPINORS-1
          DO NI=1,W1%WDES1%NIONS
             NT=W1%WDES1%ITYP(NI)
             LMMAXC=W1%WDES1%LMMAX(NT)
             IF (LMMAXC==0) CYCLE
             NPRO =ISPINOR *(W1%WDES1%NPRO/2)+W1%WDES1%LMBASE(NI)
             NPRO_=ISPINOR_*(W1%WDES1%NPRO/2)+W1%WDES1%LMBASE(NI)

             CTMP=0; CALL ECCP_NL(LMDIM,LMMAXC,CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W2%CPROJ(NPRO_+1),W1%CPROJ(NPRO+1),CTMP)
             CNL=CNL+CTMP
          ENDDO
       ENDDO
       ENDDO spinor
#ifdef _OPENACC
!$ACC EXIT DATA DELETE(W1,W2) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
#else
!$OMP END PARALLEL DO
#endif
       C=C+CNL
    ENDIF
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
    CALLMPI( M_sum_g(W1%WDES1%COMM_INB, C, 1))

    PROFILING_STOP('w1_dot')
      
  END FUNCTION W1_DOT


!***********************************************************************
!
! calculate the inproduct between one wavefunction and 
! a set of wavefunctions with possible scales and add the result
! to a third vector
!  C=   WA^* x W1 * SCALEA + C * SCALEC
!     
!***********************************************************************

  SUBROUTINE W1_GEMV( SCALEA, WA, W1, SCALEC, C, NC , CQIJ)
#ifdef _OPENACC
    USE mopenacc
#endif
    IMPLICIT NONE
    GDEF :: SCALEA                  ! scaleing constant for WA x W1
    TYPE (wavefuna)    WA
    TYPE (wavefun1)    W1
    GDEF :: SCALEC                  ! scaling constant for C
    GDEF :: C(*)                    ! result
    INTEGER NC                      ! stride of C
    OVERLAP, OPTIONAL :: CQIJ(:,:,:,:) ! optional overlap operator
  ! local
    INTEGER :: LMDIM, I, N
    GDEF CRESUL(WA%WDES1%NPRO)
    GDEF CTMP(SIZE(WA%CW,2))

    PROFILING_START('w1_gemv')

    N=SIZE(WA%CW,2)

!$ACC ENTER DATA CREATE(CTMP) __IF_ASYNC__
    IF (W1%WDES1%LGAMMA) THEN
       CALL __DGEMV__( dtrans, 2* WA%WDES1%NPL, N, done , WA%CW(1,1) , &
            &           2* WA%WDES1%NRPLWV, W1%CW(1) , 1 , dzero,  CTMP(1), 1)
    ELSE
       CALL __ZGEMV__( ztrans,  WA%WDES1%NPL, N, zone , WA%CW(1,1) , &
            &            WA%WDES1%NRPLWV, W1%CW(1) , 1 , zzero,  CTMP(1), 1)
    ENDIF

    IF (PRESENT(CQIJ) .AND. WA%WDES1%LOVERL .AND.  WA%WDES1%NPROD>0 ) THEN
       LMDIM = SIZE(CQIJ,1)
!$ACC ENTER DATA CREATE(CRESUL) __IF_ASYNC__
       CALL OVERL1(WA%WDES1, LMDIM, CQIJ(1,1,1,1), CQIJ(1,1,1,1), 0.0_q, WA%CPROJ(1,1), CRESUL(1))
       IF (W1%WDES1%LGAMMA) THEN
          CALL __DGEMV__( dtrans, WA%WDES1%NPRO,  N, done ,WA%CPROJ(1,1), &
               &           WA%WDES1%NPROD, CRESUL , 1 , done,  CTMP(1), 1)
       ELSE
          CALL __ZGEMV__( ztrans, WA%WDES1%NPRO,  N, zone ,WA%CPROJ(1,1), &
               &           WA%WDES1%NPROD, CRESUL , 1 , zone,  CTMP(1), 1)
       ENDIF
!$ACC EXIT DATA DELETE(CRESUL) __IF_ASYNC__
    ENDIF
    CALLMPI( M_sum_g(WA%WDES1%COMM_INB, CTMP, N))

    IF (SCALEC==0) THEN
!$ACC ENTER DATA CREATE(C) __IF_ASYNC__
!$ACC PARALLEL LOOP PRESENT(C,CTMP) __IF_ASYNC__
       DO I=0,N-1
          C(I*NC+1)=CTMP(I+1)*SCALEA
       ENDDO
    ELSE
!$ACC ENTER DATA COPYIN(C) __IF_ASYNC__
!$ACC PARALLEL LOOP PRESENT(C,CTMP) __IF_ASYNC__
       DO I=0,N-1
          C(I*NC+1)=C(I*NC+1)*SCALEC+CTMP(I+1)*SCALEA
       ENDDO
    ENDIF
!$ACC EXIT DATA COPYOUT(C) DELETE(CTMP) __IF_ASYNC__

    PROFILING_STOP('w1_gemv')

  END SUBROUTINE W1_GEMV


!***********************************************************************
!
! W1 descriptor from W array
! alternative to SETWAV (returns a W1 descriptor)
! it is somewhat slimmed down to optimize performance
! and returns only the wavefunction and wavefunction character pointers
! the full version remains SETWAV
! 
!***********************************************************************

  FUNCTION W1_FROM_W( W, WDES1, NB, ISP) RESULT (W1)
    IMPLICIT NONE
    INTEGER NB, ISP
    TYPE (wavespin) W
    TYPE (wavefun1) W1
    TYPE (wavedes1), TARGET :: WDES1
    INTEGER NK
    
    NK=WDES1%NK
#ifdef boundary_check
    IF (NB<=0 .OR. NB> SIZE(W%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_W: bounds exceed " // str(NB) // " " // &
          str(SIZE(W%CPTWFP,2)), __FILE__, __LINE__)
    ENDIF
#endif    
    W1%CPTWFP=>W%CPTWFP(:,NB,NK,ISP)
    W1%CPROJ =>W%CPROJ(:,NB,NK,ISP)
    W1%WDES1 => WDES1
    NULLIFY(W1%CR)
    W1%LDO=.TRUE.
  END FUNCTION W1_FROM_W

!***********************************************************************
!
! W1 descriptor from WA descriptor
!     
!***********************************************************************

  FUNCTION W1_FROM_WA( WA, N1) RESULT (W1)
    IMPLICIT NONE
    INTEGER N1
    TYPE (wavefun1) W1
    TYPE (wavefuna) WA
#ifdef boundary_check
    IF (N1<=0 .OR. N1> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1) // " " // &
          str(SIZE(WA%CPTWFP,2)), __FILE__, __LINE__)
    ENDIF
#endif    
    W1%CPTWFP=>WA%CPTWFP(:,N1)
    W1%CPROJ =>WA%CPROJ(:,N1)
    W1%WDES1 =>WA%WDES1

    IF (ASSOCIATED(WA%CR)) THEN
       W1%CR =>WA%CR(:,N1)
    ELSE
       NULLIFY(W1%CR)
    ENDIF

    W1%LDO=.TRUE.
  END FUNCTION W1_FROM_WA

!***********************************************************************
!
! subindex a WA array
! return WA(N1:N2)
!     
!***********************************************************************

  FUNCTION WA_FROM_WA( WA, N1, N2) RESULT (W1)
    IMPLICIT NONE
    INTEGER N1, N2
    TYPE (wavefuna) W1
    TYPE (wavefuna) WA
#ifdef boundary_check
    IF (N1<=0 .OR. N1> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in WA_FROM_WA: bounds exceed " // str(N1) // " " // &
          str(SIZE(WA%CPTWFP,2)), __FILE__, __LINE__)
    ENDIF
    IF (N2<N1 .OR. N2> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in WA_FROM_WA: bounds exceed " // str(N1) // " " // str(N2) &
          // " " // str(SIZE(WA%CPTWFP,2)), __FILE__, __LINE__)
    ENDIF
#endif    
    W1%CPTWFP=>WA%CPTWFP(:,N1:N2)
    W1%CPROJ =>WA%CPROJ(:,N1:N2)
    W1%WDES1 =>WA%WDES1
  END FUNCTION WA_FROM_WA

!***********************************************************************
!
! one element from a WA array using two indices
! return WA(N1,N2)
!     
!***********************************************************************

  FUNCTION W1_FROM_WA2( WA, N1, N2) RESULT (W1)
    IMPLICIT NONE
    INTEGER N1, N2
    TYPE (wavefun1) W1
    TYPE (wavefuna) WA
#ifdef boundary_check
    IF (N1<=0 .OR. N1> WA%FIRST_DIM) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1) // " " // &
          str(WA%FIRST_DIM), __FILE__, __LINE__)
    ENDIF
    IF (N1+ (N2-1)*WA%FIRST_DIM<=0 .OR. N1+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1+(N2-1)) // " " // &
          str(SIZE(WA%CPTWFP,2)), __FILE__, __LINE__)
    ENDIF
#endif
    W1%CPTWFP=>WA%CPTWFP(:,N1+ (N2-1)*WA%FIRST_DIM)
    W1%CPROJ =>WA%CPROJ (:,N1+ (N2-1)*WA%FIRST_DIM)
    W1%WDES1 =>WA%WDES1
    NULLIFY(W1%CR)
    W1%LDO=.TRUE.
  END FUNCTION W1_FROM_WA2

!***********************************************************************
!
! subindex a WA array using two indices
! return WA(N1:N12,N2)
!     
!***********************************************************************

  FUNCTION WA_FROM_WA2( WA, N1, N12, N2) RESULT (W1)
    IMPLICIT NONE
    INTEGER N1, N12, N2
    TYPE (wavefuna) W1
    TYPE (wavefuna) WA
#ifdef boundary_check
    IF (N1<=0 .OR. N1> WA%FIRST_DIM) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1) // " " // &
          str(WA%FIRST_DIM), __FILE__, __LINE__)
    ENDIF
    IF (N1+ (N2-1)*WA%FIRST_DIM<=0 .OR. N1+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1+(N2-1)) // " " // &
          str(SIZE(WA%CPTWFP,2)), __FILE__, __LINE__)
    ENDIF
    IF (N12<N1 .OR. N12> WA%FIRST_DIM) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1) // " " // str(N12) &
          // " " // str(WA%FIRST_DIM), __FILE__, __LINE__)
    ENDIF
    IF (N12+ (N2-1)*WA%FIRST_DIM<N1+ (N2-1)*WA%FIRST_DIM .OR. N12+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N12+(N2-1)) // " " // &
          str(SIZE(WA%CPTWFP,2)), __FILE__, __LINE__)
    ENDIF

#endif
    W1%CPTWFP=>WA%CPTWFP(:,N1+ (N2-1)*WA%FIRST_DIM:N12+ (N2-1)*WA%FIRST_DIM)
    W1%CPROJ =>WA%CPROJ (:,N1+ (N2-1)*WA%FIRST_DIM:N12+ (N2-1)*WA%FIRST_DIM)
    W1%WDES1 =>WA%WDES1
  END FUNCTION WA_FROM_WA2


!***********************************************************************
!
! get a WA structure from a wavefunction array
! this is equivalent to the low level routine SETWAVA
!
!***********************************************************************

  FUNCTION WA_FROM_W( W, WDES1, ISP ) RESULT (WA)
    IMPLICIT NONE
    TYPE (wavespin), INTENT(IN) :: W
    TYPE (wavedes1), TARGET ::  WDES1
    INTEGER ISP
    TYPE (wavefuna) ::  WA

    WA%CPTWFP=>W%CPTWFP(:,:,WDES1%NK,ISP)
    WA%CPROJ =>W%CPROJ(:,:,WDES1%NK,ISP)
    WA%FERWE =>W%FERWE(:,WDES1%NK,ISP)
    WA%AUX   =>W%AUX  (:,WDES1%NK,ISP)
    WA%CELEN =>W%CELEN(:,WDES1%NK,ISP)
    WA%WDES1 =>WDES1
    ! remember spin index
    WA%ISP =ISP
    ! set redistributed wavefunction indices
    IF (WDES1%DO_REDIS) THEN
       CALL SET_WPOINTER(WA%CW_RED,    WDES1%NRPLWV_RED, W%WDES%NB_TOT, WA%CW(1,1))
       CALL SET_GPOINTER(WA%CPROJ_RED, WDES1%NPROD_RED,  W%WDES%NB_TOT, WA%CPROJ(1,1))
    ELSE
       WA%CW_RED=>WA%CPTWFP
       WA%CPROJ_RED=>WA%CPROJ
    ENDIF
  END FUNCTION WA_FROM_W


!=======================================================================
!  create storage for a wavefunction array WA
!=======================================================================

  SUBROUTINE NEWWAVA(WA, WDES1, NDIM, NDIM2, ALLOC_REAL)
#ifdef _OPENACC
    USE mopenacc
#endif
    IMPLICIT NONE
    TYPE (wavefuna), TARGET :: WA
    TYPE (wavedes1), TARGET :: WDES1
    INTEGER :: NDIM
    INTEGER, OPTIONAL :: NDIM2
    LOGICAL, OPTIONAL :: ALLOC_REAL
    ! local variables
    INTEGER :: NDIMTOT

    NDIMTOT=NDIM
    WA%FIRST_DIM=0

    IF (PRESENT(NDIM2)) THEN
       NDIMTOT=NDIM*NDIM2
       WA%FIRST_DIM=NDIM
    ENDIF

    ALLOCATE(WA%CPTWFP(WDES1%NRPLWV,NDIMTOT), WA%CPROJ(WDES1%NPROD,NDIMTOT), &
             WA%CELEN(NDIMTOT), WA%FERWE(NDIMTOT), WA%AUX(NDIMTOT))
!$ACC ENTER DATA CREATE(WA%CPTWFP,WA%CPROJ) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)

DOESI WA%CPTWFP=0; WA%CPROJ=0

    ! set redistributed wavefunction indices
    IF (WDES1%DO_REDIS) THEN
       CALL SET_WPOINTER(WA%CW_RED,    WDES1%NRPLWV_RED, NDIM*WDES1%NB_PAR, WA%CW(1,1))
       CALL SET_GPOINTER(WA%CPROJ_RED, WDES1%NPROD_RED,  NDIM*WDES1%NB_PAR, WA%CPROJ(1,1))
    ELSE
       WA%CW_RED   =>WA%CPTWFP
       WA%CPROJ_RED=>WA%CPROJ
    ENDIF
!$ACC ENTER DATA CREATE(WA%CW_RED,WA%CPROJ_RED) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)

    ! optionally allocate space for the orbitals in real space
    IF (PRESENT(ALLOC_REAL)) THEN
       IF (ALLOC_REAL) THEN
          ALLOCATE(WA%CR(WDES1%GRID%MPLWV*WDES1%NRSPINORS,NDIMTOT))
!$ACC ENTER DATA CREATE(WA%CR) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
       ENDIF
    ENDIF

    WA%ISP=-1

    WA%WDES1=>WDES1
#ifdef _OPENACC
    IF (ACC_EXEC_ON.AND.ACC_IS_PRESENT(WA,1)) THEN
!$ACC UPDATE DEVICE(WA%FIRST_DIM,WA%ISP) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(WA%WDES1) ASYNC(ACC_ASYNC_Q)
!!      CALL ACC_COPYIN_TYPED_VAR(WA%WDES1)
!!      CALL ACC_ATTACH_ASYNC(WA%WDES1,ACC_ASYNC_Q)
    ENDIF
#endif
  END SUBROUTINE NEWWAVA


!=======================================================================
!  create storage for one wavefunction W array
!  to store non local part only
!=======================================================================

  SUBROUTINE NEWWAVA_PROJ(WA, WDES1, NDIM)
#ifdef _OPENACC
    USE mopenacc
#endif
    USE prec
    USE wave_mpi
    IMPLICIT NONE
    TYPE (wavefuna), TARGET :: WA
    TYPE (wavedes1), TARGET :: WDES1
    INTEGER, OPTIONAL :: NDIM
    IF (PRESENT(NDIM)) THEN
       ALLOCATE(WA%CPROJ(WDES1%NPROD,NDIM))
    ELSE
       ALLOCATE(WA%CPROJ(WDES1%NPROD,WDES1%NBANDS))
    ENDIF
!$ACC ENTER DATA CREATE(WA%CPROJ) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
    WA%FIRST_DIM=0
    WA%ISP=-1
    NULLIFY(WA%CPTWFP)
    NULLIFY(WA%CW_RED)
    NULLIFY(WA%FERWE)
    NULLIFY(WA%AUX  )
    NULLIFY(WA%CELEN)
    ! set redistributed projection indices
    IF (WDES1%DO_REDIS) THEN
       IF (PRESENT(NDIM)) THEN
          CALL SET_GPOINTER(WA%CPROJ_RED, WDES1%NPROD_RED,  NDIM*WDES1%NB_PAR, WA%CPROJ(1,1))
       ELSE
          CALL SET_GPOINTER(WA%CPROJ_RED, WDES1%NPROD_RED,  WDES1%NB_TOT, WA%CPROJ(1,1))
       ENDIF
    ELSE
       WA%CPROJ_RED=>WA%CPROJ
    ENDIF
!$ACC ENTER DATA CREATE(WA%CPROJ_RED) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
    WA%WDES1=>WDES1
#ifdef _OPENACC
    IF (ACC_EXEC_ON.AND.ACC_IS_PRESENT(WA,1)) THEN
!$ACC UPDATE DEVICE(WA%FIRST_DIM,WA%ISP) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(WA%WDES1) ASYNC(ACC_ASYNC_Q)
!!      CALL ACC_COPYIN_TYPED_VAR(WA%WDES1)
!!      CALL ACC_ATTACH_ASYNC(WA%WDES1,ACC_ASYNC_Q)
    ENDIF
#endif
  END SUBROUTINE NEWWAVA_PROJ


!=======================================================================
!  create storage for an array of single orbital structures W1(:)
!=======================================================================

  SUBROUTINE NEWWAV_ARRAY(W1,WDES1,ALLOC_REAL,WA)
#ifdef _OPENACC
    USE mopenacc
#endif
    IMPLICIT NONE

    TYPE (wavefun1) :: W1(:)
    TYPE (wavedes1) :: WDES1
    TYPE (wavefuna), OPTIONAL :: WA

    LOGICAL :: ALLOC_REAL

    ! local variables
    INTEGER :: NDIM,I

    NDIM=SIZE(W1)

    IF (PRESENT(WA)) THEN
!
! Allocate a wavefuna structure WA to create storage space for the pointer
! members of the array of single orbital structures W1(:).
! This guarantees that the respective orbital data arrays W1(:)%CPTWFP and
! W1(:)%CR will be laid out contiguously in memory.
!
!$ACC ENTER DATA CREATE(WA) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(W1(:),1)) ASYNC(ACC_ASYNC_Q)
       ! allocate a wavefuna structure WA of size NDIM
       CALL NEWWAVA(WA, WDES1, NDIM, ALLOC_REAL=ALLOC_REAL)
       ! connect the W1 members to WA
       DO I=1,NDIM
          ! connect pointers
          W1(I) = ELEMENT(WA, I)
!$ACC ENTER DATA COPYIN(W1(I)%CPTWFP,W1(I)%CPROJ) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(W1(I),1)) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(W1(I)%CR) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(W1(I),1).AND.ALLOC_REAL) ASYNC(ACC_ASYNC_Q)
!!!$ACC ENTER DATA CREATE(W1(I)%CPTWFP,W1(I)%CPROJ) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(W1(I),1)) ASYNC(ACC_ASYNC_Q)
!!!$ACC ENTER DATA CREATE(W1(I)%CR) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(W1(I),1).AND.ALLOC_REAL) ASYNC(ACC_ASYNC_Q)

          ! set other members of W1 (see NEWWAV)
          W1(I)%FERWE= 0
          W1(I)%CELEN= 0
          W1(I)%NB   =-1
          W1(I)%ISP  =-1
          W1(I)%LDO  = .TRUE.

#ifdef _OPENACC
          IF (ACC_EXEC_ON.AND.ACC_IS_PRESENT(W1(I),1)) THEN
!$ACC UPDATE DEVICE(W1(I)%FERWE,W1(I)%CELEN,W1(I)%NB,W1(I)%ISP,W1(I)%LDO) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(W1(I)%WDES1) ASYNC(ACC_ASYNC_Q)
          ENDIF
#endif
       ENDDO
    ELSE
!
! In case WA is not present we fall back to the use of NEWWAV.
!
       DO I=1,NDIM
          CALL NEWWAV(W1(I), WDES1, ALLOC_REAL)
       ENDDO
    ENDIF
  END SUBROUTINE NEWWAV_ARRAY


!=======================================================================
!  destroy storage for an array of single orbital structures W1(:)
!  that was possibly created using a wavefuna structure
!=======================================================================

  SUBROUTINE DELWAV_ARRAY(W1,DEALLOC_REAL,WA)
#ifdef _OPENACC
    USE mopenacc
#endif
    IMPLICIT NONE

    TYPE (wavefun1) :: W1(:)
    TYPE (wavefuna), OPTIONAL :: WA

    LOGICAL :: DEALLOC_REAL

    ! local variables
    INTEGER :: NDIM,I

    NDIM=SIZE(W1)

    IF (PRESENT(WA)) THEN
       IF (NDIM/=SIZE(WA%CPTWFP,2)) CALL vtutor%bug("DELWAV_ARRAY: partial deallocation is not supported",__FILE__,__LINE__)
       ! deallocate WA structure
       CALL DELWAVA(WA)
!$ACC EXIT DATA DELETE(WA) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
       ! nullify W1 pointer members
       DO I=1,NDIM
!$ACC EXIT DATA DELETE(W1(I)%CPTWFP,W1(I)%CPROJ,W1(I)%WDES1) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(W1(I),1)) ASYNC(ACC_ASYNC_Q)
          NULLIFY(W1(I)%CPTWFP,W1(I)%CPROJ,W1(I)%WDES1)
          IF (DEALLOC_REAL) THEN
!$ACC EXIT DATA DELETE(W1(I)%CR) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(W1(I),1)) ASYNC(ACC_ASYNC_Q)
             NULLIFY(W1(I)%CR)
          ENDIF
       ENDDO
    ELSE
       DO I=1,NDIM
          CALL DELWAV(W1(I), DEALLOC_REAL)
       ENDDO
    ENDIF
  END SUBROUTINE DELWAV_ARRAY


!=======================================================================
!  Connect the pointer members of a wavefuna structure to the members
!  of an array of single orbital structures
!  Mind: the latter has to have been created using NEWWAV_ARRAY, to
!  ensure that things like W1(:)%CPTWFP are contiguous in memory
!
!  WA_CONNECT will try to make sure that the members of WA and W1(:)
!  are properly mapped and signal success by returning WA_CONNECT=0 or 1
!  The latter means that WA%CR is associated and one may perform
!  Fast-Fourier-Transfers with this structure.
!
!=======================================================================

  FUNCTION WA_CONNECT(W1,WA)
#ifdef _OPENACC
    USE mopenacc
#endif
    USE iso_c_binding
    IMPLICIT NONE
    TYPE (wavefun1) :: W1(:)
    TYPE (wavefuna) :: WA

    INTEGER :: WA_CONNECT

    ! local variables
    INTEGER :: I

    CALL c_f_pointer(c_loc(W1(1)%CPTWFP(1)),WA%CPTWFP,[SIZE(W1(1)%CPTWFP),SIZE(W1)])
    CALL c_f_pointer(c_loc(W1(1)%CPROJ (1)),WA%CPROJ ,[SIZE(W1(1)%CPROJ ),SIZE(W1)])

    IF (ASSOCIATED(W1(1)%CR)) CALL c_f_pointer(c_loc(W1(1)%CR(1)),WA%CR,[SIZE(W1(1)%CR),SIZE(W1)])

    WA%WDES1 => W1(1)%WDES1

    ! test consistency
    WA_CONNECT=0
    DO I=1,SIZE(W1)
       IF (.NOT.c_associated(c_loc(WA%CPTWFP(1,I)),c_loc(W1(I)%CPTWFP(1)))) &
          CALL vtutor%bug("WA_CONNECT: target appears to be non-contiguous",__FILE__,__LINE__)
    ENDDO

    ! signal the fact that we can do FFTs with this structure
    IF (WA_CONNECT==0.AND.ASSOCIATED(WA%CR)) WA_CONNECT=1

!$ACC ENTER DATA COPYIN(WA) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(W1,1)) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(WA%CPTWFP,WA%CPROJ,WA%WDES1) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(W1,1)) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(WA%CR) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(W1,1).AND.ASSOCIATED(WA%CR)) ASYNC(ACC_ASYNC_Q)
  END FUNCTION WA_CONNECT


!=======================================================================
!  nullify the pointer members of a wavefuna structure
!=======================================================================

  SUBROUTINE WA_NULLIFY(WA)
#ifdef _OPENACC
    USE mopenacc
#endif
    IMPLICIT NONE
    TYPE (wavefuna) :: WA
!$ACC EXIT DATA DELETE(WA%CR) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(WA,1).AND.ASSOCIATED(WA%CR)) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA DELETE(WA%CPTWFP,WA%CPROJ,WA%WDES1) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA DELETE(WA) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(WA,1))  ASYNC(ACC_ASYNC_Q)
    NULLIFY(WA%CPTWFP,WA%CPROJ,WA%WDES1); IF (ASSOCIATED(WA%CR)) NULLIFY(WA%CR)
  END SUBROUTINE WA_NULLIFY


!=======================================================================
!  batched FFT (q -> r) of a group of orbitals
!=======================================================================

  SUBROUTINE FFTWAV_WA(WA)
    IMPLICIT NONE
    TYPE (wavefuna) :: WA
    ! local variables
    INTEGER :: ISPINOR

    IF (WA%WDES1%NK==0) &
       CALL vtutor%bug("FFTWAV_WA: NK is set to zero: WDES not set up properly", __FILE__, __LINE__)

    IF (.NOT.ASSOCIATED(WA%CR)) &
       CALL vtutor%bug("FFTWAV_WA: can not perform FFT: CR is not allocated", __FILE__, __LINE__)

    IF (WA%WDES1%LUSEINV) THEN
       DO ISPINOR=0,WA%WDES1%NRSPINORS-1
          CALL FFTWAV_USEINV_MU(WA%WDES1%NGVECTOR,SIZE(WA%CW,2),WA%WDES1%NINDPW(1), &
                                WA%WDES1%NINDPW_INV(1),WA%WDES1%FFTSCA(1,2), &
                                WA%CR(1+ISPINOR*WA%WDES1%GRID%MPLWV,1),SIZE(WA%CR,1), &
                                WA%CW(1+ISPINOR*WA%WDES1%NGVECTOR,1),SIZE(WA%CW,1), &
                                WA%WDES1%GRID)
       ENDDO
    ELSE
       DO ISPINOR=0,WA%WDES1%NRSPINORS-1
          CALL FFTWAV_MU(WA%WDES1%NGVECTOR,SIZE(WA%CW,2),WA%WDES1%NINDPW(1), &
                         WA%CR(1+ISPINOR*WA%WDES1%GRID%MPLWV,1),SIZE(WA%CR,1), &
                         WA%CW(1+ISPINOR*WA%WDES1%NGVECTOR,1),SIZE(WA%CW,1), &
                         WA%WDES1%GRID)
       ENDDO
    ENDIF
  END SUBROUTINE FFTWAV_WA


!=======================================================================
!  batched FFT (r -> q) of a group of orbitals
!=======================================================================

  SUBROUTINE FFTEXT_WA(WA,LADD)
    IMPLICIT NONE
    TYPE (wavefuna) :: WA
    LOGICAL :: LADD
    ! local variables
    INTEGER :: ISPINOR

    IF (WA%WDES1%NK==0) &
       CALL vtutor%bug("FFTEXT_WA: NK is set to zero: WDES not set up properly", __FILE__, __LINE__)

    IF (.NOT.ASSOCIATED(WA%CR)) &
       CALL vtutor%bug("FFTEXT_WA: can not perform FFT: CR is not allocated", __FILE__, __LINE__)

    IF (WA%WDES1%LUSEINV) THEN
       DO ISPINOR=0,WA%WDES1%NRSPINORS-1
          CALL FFTEXT_USEINV_MU(WA%WDES1%NGVECTOR,SIZE(WA%CW,2),WA%WDES1%NINDPW(1),WA%WDES1%FFTSCA(1,1), &
                                WA%CR(1+ISPINOR*WA%WDES1%GRID%MPLWV,1),SIZE(WA%CR,1), &
                                WA%CW(1+ISPINOR*WA%WDES1%NGVECTOR,1),SIZE(WA%CW,1), &
                                WA%WDES1%GRID,LADD)
       ENDDO
    ELSE
       DO ISPINOR=0,WA%WDES1%NRSPINORS-1
          CALL FFTEXT_MU(WA%WDES1%NGVECTOR,SIZE(WA%CW,2),WA%WDES1%NINDPW(1), &
                         WA%CR(1+ISPINOR*WA%WDES1%GRID%MPLWV,1),SIZE(WA%CR,1), &
                         WA%CW(1+ISPINOR*WA%WDES1%NGVECTOR,1),SIZE(WA%CW,1), &
                         WA%WDES1%GRID,LADD)
       ENDDO
    ENDIF
  END SUBROUTINE FFTEXT_WA


!=======================================================================
!  set one single wavefunction array (WA) from an array of wavefunctions
!=======================================================================

  SUBROUTINE SETWAVA(W, WA, WDES1, ISP)
    USE prec
    IMPLICIT NONE
    INTEGER NB,ISP
    TYPE (wavespin) W
    TYPE (wavefuna) WA
    TYPE (wavedes1), TARGET :: WDES1
    INTEGER NK

    NK=WDES1%NK

    WA%CPTWFP=>W%CPTWFP(:,:,NK,ISP)
    WA%CPROJ =>W%CPROJ(:,:,NK,ISP)
    WA%FERWE =>W%FERWE(:,NK,ISP)
    WA%AUX   =>W%AUX  (:,NK,ISP)
    WA%CELEN =>W%CELEN(:,NK,ISP)
    WA%WDES1 =>WDES1
    ! one dimensional indexing assumed
    WA%FIRST_DIM=0
    ! remember spin index
    WA%ISP =ISP
    ! set redistributed wavefunction indices
    IF (WDES1%DO_REDIS) THEN
       CALL SET_WPOINTER(WA%CW_RED,    WA%WDES1%NRPLWV_RED, W%WDES%NB_TOT, WA%CW(1,1))
       CALL SET_GPOINTER(WA%CPROJ_RED, WA%WDES1%NPROD_RED,  W%WDES%NB_TOT, WA%CPROJ(1,1))
    ELSE
       WA%CW_RED=>WA%CPTWFP
       WA%CPROJ_RED=>WA%CPROJ
    ENDIF

  END SUBROUTINE SETWAVA

!=======================================================================
!  destroy storage for a wavefunctionarray WA
!=======================================================================

      SUBROUTINE DELWAVA(WA)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      IMPLICIT NONE
      TYPE (wavefuna) WA

      IF (.NOT. ASSOCIATED(WA%CPTWFP) .OR. .NOT. ASSOCIATED(WA%FERWE) .OR. &
          .NOT. ASSOCIATED(WA%CELEN)  .OR. .NOT. ASSOCIATED(WA%AUX)) THEN
         CALL vtutor%bug("internal error in DELWAVA: not all enities are associated, try DELWAVA_PROJ", __FILE__, __LINE__)
      ENDIF

!$ACC EXIT DATA DELETE(WA%CPTWFP,WA%CPROJ,WA%CW_RED,WA%CPROJ_RED) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
      DEALLOCATE(WA%CPTWFP,WA%CPROJ,WA%CELEN, WA%FERWE, WA%AUX)
      NULLIFY(WA%CW_RED,WA%CPROJ_RED)

      IF (ASSOCIATED(WA%CR)) THEN
!$ACC EXIT DATA DELETE(WA%CR) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
         DEALLOCATE(WA%CR)
      ENDIF

#ifdef _OPENACC
      IF (ACC_EXEC_ON.AND.ACC_IS_PRESENT(WA,1)) THEN
!$ACC EXIT DATA DELETE(WA%WDES1) ASYNC(ACC_ASYNC_Q)
!!        CALL ACC_DELETE_TYPED_VAR(WA%WDES1)
!!        CALL ACC_DETACH_ASYNC(WA%WDES1,ACC_ASYNC_Q)
      ENDIF
#endif
      NULLIFY(WA%WDES1)
!!$ACC EXIT DATA DELETE(WA) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
      END SUBROUTINE DELWAVA


!=======================================================================
!  destroy storage for a wavefunctionarray WA
!  wavefunction character  only
!=======================================================================

      SUBROUTINE DELWAVA_PROJ(WA)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      IMPLICIT NONE
      TYPE (wavefuna) WA

      IF (ASSOCIATED(WA%CPTWFP) .OR. ASSOCIATED(WA%FERWE) .OR. &
          ASSOCIATED(WA%CELEN) .OR. ASSOCIATED(WA%AUX)) THEN
         CALL vtutor%bug("internal error in DELWAVA_PROJ: enities are associated, try DELWAVA", __FILE__, __LINE__)
      ENDIF
!$ACC EXIT DATA DELETE(WA%CPROJ,WA%CPROJ_RED) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
      DEALLOCATE(WA%CPROJ)
      NULLIFY(WA%CPROJ_RED)
#ifdef _OPENACC
      IF (ACC_EXEC_ON.AND.ACC_IS_PRESENT(WA,1)) THEN
!$ACC EXIT DATA DELETE(WA%WDES1) ASYNC(ACC_ASYNC_Q)
!!        CALL ACC_DELETE_TYPED_VAR(WA%WDES1)
!!        CALL ACC_DETACH_ASYNC(WA%WDES1,ACC_ASYNC_Q)
      ENDIF
#endif
      NULLIFY(WA%WDES1)
!!$ACC EXIT DATA DELETE(WA) IF(ACC_EXEC_ON.AND.ACC_IS_PRESENT(WA,1)) ASYNC(ACC_ASYNC_Q)
    END SUBROUTINE DELWAVA_PROJ


!***********************************************************************
!
! index the wavefunction array or character array in WA
!     
!***********************************************************************

  FUNCTION PCW( WA, N1, N2)
    INTEGER N1, N2
    TYPE (wavefuna) WA
    COMPLEX(q), POINTER :: PCW(:)
#ifdef boundary_check
    IF (N1<=0 .OR. N1> WA%FIRST_DIM) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1) // " " // &
          str(WA%FIRST_DIM), __FILE__, __LINE__)
    ENDIF
    IF (N1+ (N2-1)*WA%FIRST_DIM<=0 .OR. N1+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1+(N2-1)) // " " // &
          str(SIZE(WA%CPTWFP,2)), __FILE__, __LINE__)
    ENDIF
#endif
    PCW=>WA%CPTWFP(:,N1+ (N2-1)*WA%FIRST_DIM)
  END FUNCTION PCW

  FUNCTION PCPROJ( WA, N1, N2)
    INTEGER N1, N2
    TYPE (wavefuna) WA
    GDEF, POINTER :: PCPROJ(:)
#ifdef boundary_check
    IF (N1<=0 .OR. N1> WA%FIRST_DIM) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1) // " " // &
          str(WA%FIRST_DIM), __FILE__, __LINE__)
    ENDIF
    IF (N1+ (N2-1)*WA%FIRST_DIM<=0 .OR. N1+ (N2-1)*WA%FIRST_DIM> SIZE(WA%CPTWFP,2)) THEN
       CALL vtutor%bug("internal error in W1_FROM_WA: bounds exceed " // str(N1+(N2-1)) // " " // &
          str(SIZE(WA%CPTWFP,2)), __FILE__, __LINE__)
    ENDIF
#endif
    PCPROJ=>WA%CPROJ(:,N1+ (N2-1)*WA%FIRST_DIM)
  END FUNCTION PCPROJ



!************************* SUBROUTINE ORTHON ***************************
!
! orthogonalize a wavefunction W1 to all other bands
! including the current band
! the subroutine uses BLAS 3 calls,
!
!***********************************************************************
    
  SUBROUTINE ORTHON(NK, W, W1, CQIJ, ISP)
    IMPLICIT NONE

    INTEGER NK
    TYPE (wavespin)   W
    TYPE (wavefun1)   W1
    LOGICAL LOVERL
    OVERLAP CQIJ(:,:,:,:)
    INTEGER ISP
  ! local
    GDEF :: CPRO(W%WDES%NBANDS),CWORK(W%WDES%NPRO)
    REAL(q) :: WFMAG
    INTEGER :: I

#ifdef MPI
    IF (W%WDES%COMM_KIN%NCPU /= W%WDES%COMM_INB%NCPU) THEN
       CALL vtutor%bug("internal error: ORTHON does not support band-par.", __FILE__, __LINE__)
    ENDIF
#endif
    IF (W1%WDES1%LOVERL) THEN
       CALL OVERL1(W1%WDES1, SIZE(CQIJ,1),CQIJ(1,1,1,ISP),CQIJ(1,1,1,ISP), 0.0_q, W1%CPROJ(1),CWORK(1))
    ENDIF

    CALL GGEMV( trans , m_ W%WDES%NPLWKP(NK) , W%WDES%NBANDS ,one , W%CW(1,1,NK,ISP), &
         &             m_ W%WDES%NRPLWV, W1%CW(1) , 1 , zero ,  CPRO(1), 1)

    IF (W1%WDES1%LOVERL) THEN
       IF (W%WDES%NPRO /= 0) &
            CALL GGEMV( trans ,  W%WDES%NPRO , W%WDES%NBANDS ,one , W%CPROJ(1,1,NK,ISP) , &
            W%WDES%NPROD, CWORK(1), 1 , one ,  CPRO(1), 1)
    ENDIF

    CALLMPI( M_sum_g(W%WDES%COMM_KIN, CPRO(1), W%WDES%NBANDS))

    CALL GGEMM( 'N', 'N' , m_ W%WDES%NPLWKP(NK) , 1 , W%WDES%NBANDS , -one , &
         W%CW(1,1,NK,ISP), m_ W%WDES%NRPLWV , CPRO(1) , W%WDES%NBANDS , &
         one , W1%CW(1) , m_ W%WDES%NRPLWV )

    IF (W%WDES%NPRO /= 0) &
         CALL GGEMM( 'N', 'N' ,  W%WDES%NPRO , 1 , W%WDES%NBANDS  , -one , &
         W%CPROJ(1,1,NK,ISP) ,  W%WDES%NPROD , CPRO(1) , W%WDES%NBANDS , &
         one , W1%CPROJ(1) ,  W%WDES%NPROD  )

  END SUBROUTINE ORTHON

  SUBROUTINE ORTHON1P(NK, W, W1, CQIJ, ISP, NB)
    IMPLICIT NONE

    INTEGER NK
    TYPE (wavespin)   W
    TYPE (wavefun1)   W1(:)
    LOGICAL LOVERL
    OVERLAP CQIJ(:,:,:,:)
    INTEGER ISP,NB(:)
  ! local
    GDEF :: CPRO(W%WDES%NBANDS,SIZE(W1)),CWORK(W%WDES%NPRO,SIZE(W1))
    GDEF,ALLOCATABLE :: CW1(:,:),CPROJ1(:,:)
    REAL(q) :: WFMAG
    INTEGER :: I,BW

    BW=128

#ifdef MPI
    IF (W%WDES%COMM_KIN%NCPU /= W%WDES%COMM_INB%NCPU) THEN
       CALL vtutor%bug("internal error: ORTHON does not support band-par.", __FILE__, __LINE__)
    ENDIF
#endif

    ALLOCATE(CW1(SIZE(W1(1)%CW,1),SIZE(W1)))
    DO I=1,SIZE(W1)
       CW1(:,I)=W1(I)%CW
    END DO

    DO I=1,SIZE(W1)
       IF (W1(I)%WDES1%LOVERL)THEN
          CALL OVERL1(W1(I)%WDES1, SIZE(CQIJ,1),CQIJ(1,1,1,ISP),CQIJ(1,1,1,ISP), 0.0_q, W1(I)%CPROJ(1),CWORK(1,I))
       ELSE
          CWORK(:,I)=0.0_q
       END IF
    END DO

    CALL GGEMM( trans , 'N', W%WDES%NBANDS, size(w1), m_ W%WDES%NPLWKP(NK), &
         & one , W%CW(1,1,NK,ISP), m_ W%WDES%NRPLWV, CW1(1,1) , m_ W%WDES%NRPLWV, zero ,  CPRO(1,1), W%WDES%NBANDS)


    IF (W%WDES%NPRO /= 0) &
         CALL GGEMM( trans , 'N', W%WDES%NBANDS , size(w1), W%WDES%NPRO , one ,&
         & W%CPROJ(1,1,NK,ISP) , W%WDES%NPROD, CWORK(1,1), W%WDES%NPRO, one,  CPRO(1,1), W%WDES%NBANDS)

    CALLMPI( M_sum_g(W%WDES%COMM_KIN, CPRO(1,1), W%WDES%NBANDS*SIZE(W1)))

    DO I=1,SIZE(W1)
       CPRO(NB(I),I)=0.d0
       IF(NB(I)-BW>=1)CPRO(1:NB(I)-BW,I)=0.d0
       IF(NB(I)+BW<=SIZE(CPRO,1))CPRO(NB(I)+BW:SIZE(CPRO,1),I)=0.d0
    END DO

    CALL GGEMM( 'N', 'N' , m_ W%WDES%NPLWKP(NK) , SIZE(W1) , W%WDES%NBANDS , -one , &
         W%CW(1,1,NK,ISP), m_ W%WDES%NRPLWV , CPRO(1,1) , W%WDES%NBANDS , &
         one , CW1(1,1) , m_ W%WDES%NRPLWV )

    IF (W%WDES%NPRO /= 0)THEN
       DO I=1,SIZE(W1)
          CWORK(:,i)=W1(i)%CPROJ(:)
        END DO
       CALL GGEMM( 'N', 'N' ,  W%WDES%NPRO , SIZE(W1) , W%WDES%NBANDS  , -one , &
            W%CPROJ(1,1,NK,ISP) ,  W%WDES%NPROD , CPRO(1,1) , W%WDES%NBANDS , &
            one , CWORK(1,1) ,  W%WDES%NPROD  )
    end IF

    DO I=1,SIZE(W1)
       W1(I)%CW=CW1(:,I)
       W1(I)%CPROJ=CWORK(:,I)
    END DO

    DEALLOCATE(CW1)

  END SUBROUTINE ORTHON1P

!************************* SUBROUTINE CNORMN  **************************
!
! this subroutine normalises a wavefunction
! subroutine is not important for performance
!
!***********************************************************************

  SUBROUTINE CNORMN(W, CQIJ, ISP, WSCAL)
#ifdef _OPENACC
    USE mopenacc
#endif
    IMPLICIT NONE
    TYPE (wavefun1) W

    OVERLAP CQIJ(:,:,:,:)
    INTEGER ISP
    REAL(q) WSCAL
  ! local
    GDEF    CP
    REAL(q) WFMAG
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, LMMAXC, NPRO2, NPRO2_

    COMPLEX(q), EXTERNAL ::  __ZDOTC__
!$ACC ROUTINE(ECCP_NL) VECTOR

    WFMAG=__ZDOTC__(W%WDES1%NPL,W%CW(1),1,W%CW(1),1)
!=======================================================================
! if necessary caclulate <w| P |w>
!=======================================================================
    IF (W%WDES1%LOVERL) THEN
       CP  =0

!$ACC PARALLEL LOOP COLLAPSE(force:3) PRESENT(W,CQIJ) REDUCTION(+:CP) &
!$ACC PRIVATE(NPRO,NPRO_,NPRO2,NPRO2_,NT,LMMAXC) DEFAULT(none) __IF_ASYNC__
       spinor: DO ISPINOR=0,W%WDES1%NRSPINORS-1
          DO ISPINOR_=0,W%WDES1%NRSPINORS-1

             NPRO =ISPINOR *(W%WDES1%NPRO/2)
             NPRO_=ISPINOR_*(W%WDES1%NPRO/2)

             DO NI=1,W%WDES1%NIONS
                NT=W%WDES1%ITYP(NI)
                LMMAXC=W%WDES1%LMMAX(NT)
                IF (LMMAXC==0) CYCLE

                NPRO2 =W%WDES1%LMBASE(NI)+NPRO
                NPRO2_=W%WDES1%LMBASE(NI)+NPRO_
                CALL ECCP_NL(SIZE(CQIJ,1),LMMAXC,CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W%CPROJ(NPRO2_+1),W%CPROJ(NPRO2+1),CP)
             ENDDO
          ENDDO
       ENDDO spinor

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
       WFMAG=WFMAG+CP
    ENDIF

    CALLMPI( M_sum_d(W%WDES1%COMM_INB, WFMAG, 1))

    !-----check that it is non-zero
    IF(WFMAG<=0) THEN
!=======================================================================
! if it is smaller zero write a warning
!=======================================================================
#ifdef MPI
       IF (W%WDES1%COMM_INB%NODE_ME == W%WDES1%COMM_INB%IONODE) THEN
#endif
          WRITE(*,*)'WARNING: CNORMN: search vector ill defined'
#ifdef MPI
       ENDIF
#endif
       WSCAL= -1._q/SQRT(-WFMAG)
    ELSE
       WSCAL= 1._q/SQRT(WFMAG)
    ENDIF
    CALL __ZDSCAL__( W%WDES1%NPL ,WSCAL,W%CW(1),1)
    CALL __GDSCAL__( W%WDES1%NPRO,WSCAL,W%CPROJ(1),1)

  END SUBROUTINE CNORMN


!************************* SUBROUTINE CNORMN_REAL **********************
!
! performs operations on real space part of wavefunction
! after a call to CPROJCN and CNORMN
!
!***********************************************************************

  SUBROUTINE CNORMN_REAL(W, W1, ISP, WSCAL, CSCPD )
#ifdef _OPENACC
    USE mopenacc_struct_def
#endif
    IMPLICIT NONE
    TYPE (wavefun1)    W, W1

    INTEGER   ISP
    REAL(q) WSCAL
    COMPLEX(q) :: CSCPD
  ! local
    INTEGER ISPINOR, K, KK

!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(W,W%WDES1,W%WDES1%GRID,W%WDES1%GRID%RL,W1) PRIVATE(KK) __IF_ASYNC__
    DO ISPINOR=0,W%WDES1%NRSPINORS-1
       DO K=1,W%WDES1%GRID%RL%NP
          KK=K+ISPINOR*W%WDES1%GRID%MPLWV
          W%CR(KK)=(W%CR(KK)-CSCPD*W1%CR(KK))*WSCAL
       ENDDO
    ENDDO

  END SUBROUTINE CNORMN_REAL


!************************* SUBROUTINE CNORMA  **************************
!
! this subroutine calculates the norm of a wavefunction
!
!> @details @ref openmp :
!> a ZDOTC call is replaced by an explicit OMP PARALLEL loop
!> (for some reason BLAS1 call thread really badly).
!> The nested loop over \"types\" + \"ions-of-type\" is replaced by
!> a loop over \"all ions\", that is distributed over all available
!> threads.
!
!***********************************************************************

  SUBROUTINE CNORMA(W, CQIJ, ISP, WSCAL)
#ifdef _OPENACC
    USE mopenacc
#endif
    IMPLICIT NONE

    TYPE (wavefun1) W
    OVERLAP :: CQIJ(:,:,:,:)
    INTEGER :: ISP
    REAL(q) :: WSCAL
  ! local
    GDEF    :: CP, CTMP
    REAL(q) :: WFMAG
    INTEGER :: ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NI, LMMAXC, LMDIM

    REAL(q), EXTERNAL ::  __DDOT__

!$  INTEGER I
!$ACC ROUTINE(ECCP_NL) VECTOR

    PROFILING_START('cnorma')

#if defined(_OPENACC) || ! defined(_OPENMP)
!      WFMAG=ZDOTC(W%WDES1%NPL,W%CW(1),1,W%CW(1),1)
       WFMAG=__DDOT__(W%WDES1%NPL*2,W%CW(1),1,W%CW(1),1)
#else
       WFMAG=0
!$OMP PARALLEL DO SIMD SHARED(W) REDUCTION(+:WFMAG)
       DO I=1,W%WDES1%NPL
          WFMAG=WFMAG+CONJG(W%CW(I))*W%CW(I)
       ENDDO
!$OMP END PARALLEL DO SIMD
#endif
!=======================================================================
! if necessary caclulate <w| P |w>
!=======================================================================
    IF (W%WDES1%LOVERL) THEN

       CP=0
       LMDIM=SIZE(CQIJ,1)

       PROFILING_START('eccp_nl')

#ifdef _OPENACC
!$ACC ENTER DATA COPYIN(W) __IF_ASYNC__
!$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(CQIJ,W) REDUCTION(+:CP) &
!$ACC PRIVATE(NT,LMMAXC,NPRO,NPRO_,CTMP) __IF_ASYNC__
#else
!$OMP PARALLEL DO COLLAPSE(3) SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP PRIVATE(ISPINOR,ISPINOR_,NI,NT,LMMAXC,NPRO,NPRO_,CTMP) &
!$OMP SHARED(W,CQIJ,LMDIM,ISP) REDUCTION(+:CP)
#endif
       spinor: DO ISPINOR=0,W%WDES1%NRSPINORS-1
          DO ISPINOR_=0,W%WDES1%NRSPINORS-1
             DO NI=1,W%WDES1%NIONS
                NT=W%WDES1%ITYP(NI)
                LMMAXC=W%WDES1%LMMAX(NT)
                IF (LMMAXC==0) CYCLE
                NPRO =ISPINOR *(W%WDES1%NPRO/2)+W%WDES1%LMBASE(NI)
                NPRO_=ISPINOR_*(W%WDES1%NPRO/2)+W%WDES1%LMBASE(NI)

                CTMP=0; CALL ECCP_NL(LMDIM,LMMAXC,CQIJ(1,1,NI,ISP+ISPINOR_+2*ISPINOR),W%CPROJ(NPRO_+1),W%CPROJ(NPRO+1),CTMP)
                CP=CP+CTMP
             ENDDO
          ENDDO
       ENDDO spinor
#ifdef _OPENACC
!$ACC EXIT DATA DELETE(W) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
#else
!$OMP END PARALLEL DO
#endif
       PROFILING_STOP('eccp_nl')

       WFMAG=WFMAG+CP
    ENDIF
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
    CALLMPI( M_sum_d(W%WDES1%COMM_INB, WFMAG, 1))

    !-----check that it is non-zero
    IF(WFMAG<=0) THEN
!=======================================================================
! if it is smaller zero write a warning
!=======================================================================
#ifdef MPI
       IF (W%WDES1%COMM_INB%NODE_ME == W%WDES1%COMM_INB%IONODE) THEN
#endif
          WRITE(*,*)'WARNING: CNORMN: search vector ill defined'
#ifdef MPI
       ENDIF
#endif
       WSCAL= -1._q/SQRT(-WFMAG)
    ELSE
       WSCAL= 1._q/SQRT(WFMAG)
    ENDIF

    PROFILING_STOP('cnorma')

  END SUBROUTINE CNORMA

!************************* SUBROUTINE CINPROD  *************************
!
! this subroutine calculates the inproduct between two wavefunctions
!  <W2 | S | W1>      = W2(G)*  W1(G) +  \sum_i W2_i* Q_ij W1_j
!***********************************************************************

  SUBROUTINE CINPROD(W1,W2,CQIJ,CWFMAG)
    IMPLICIT NONE

    TYPE (wavefun1)    W1
    TYPE (wavefun1)    W2

    OVERLAP   CQIJ(:,:,:,:)
    COMPLEX(q) CWFMAG
  ! local
    GDEF      CP
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, LMMAXC
    COMPLEX(q), EXTERNAL :: ZDOTC
    REAL(q), EXTERNAL ::  DDOT

    CWFMAG=ZDOTC(W1%WDES1%NPL,W2%CW(1),1,W1%CW(1),1)
    CP  =0
!=======================================================================
! if necessary caclulate <w| P |w>
!=======================================================================
    IF (W1%WDES1%LOVERL) THEN
       NPRO=0

       spinor: DO ISPINOR=0,W1%WDES1%NRSPINORS-1
          DO ISPINOR_=0,W1%WDES1%NRSPINORS-1

             NPRO =ISPINOR *(W1%WDES1%NPRO/2)
             NPRO_=ISPINOR_*(W1%WDES1%NPRO/2)

             NIS =1
             DO NT=1,W1%WDES1%NTYP
                LMMAXC=W1%WDES1%LMMAX(NT)
                IF (LMMAXC==0) GOTO 230

                DO NI=NIS,W1%WDES1%NITYP(NT)+NIS-1
                   CALL ECCP_NL(SIZE(CQIJ,1),LMMAXC,CQIJ(1,1,NI,1+ISPINOR_+2*ISPINOR),W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CP)
                   NPRO = LMMAXC+NPRO
                   NPRO_= LMMAXC+NPRO_
                ENDDO
230             NIS = NIS+W1%WDES1%NITYP(NT)
             ENDDO
          ENDDO
       ENDDO spinor

       CWFMAG=CWFMAG+CP
    ENDIF

    CALLMPI( M_sum_d(W1%WDES1%COMM_INB, CWFMAG, 2))
  END SUBROUTINE CINPROD


!************************* SUBROUTINE PROJCN ***************************
!
! this subroutine projects out from one wavefunction
! CF another wavefunction  CPRO
! subroutine is not important for performance
!
!***********************************************************************

  SUBROUTINE PROJCN(W1, W2, CQIJ, ISP, CSCPD)
#ifdef _OPENACC
    USE mopenacc
#endif
    IMPLICIT NONE

    TYPE (wavefun1)    W1,W2
    GDEF    CADD
    OVERLAP CQIJ(:,:,:,:)
    INTEGER :: ISP
    COMPLEX(q) :: CSCPD
    ! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, LMMAXC
    COMPLEX(q), EXTERNAL ::  __ZDOTC__
!$ACC ROUTINE(ECCP_NL) VECTOR

    CSCPD= GREAL(__ZDOTC__(W1%WDES1%NPL,W2%CW(1),1,W1%CW(1),1))
!=======================================================================
! if necessary caclulate <p| P |w>
!=======================================================================
    IF (W1%WDES1%LOVERL) THEN
       CADD=0
!$ACC PARALLEL LOOP COLLAPSE(3) DEFAULT(NONE) PRESENT(W1,W1%WDES1,W2,CQIJ) REDUCTION(+:CADD) &
!$ACC PRIVATE(NPRO,NPRO_,NT,LMMAXC) __IF_ASYNC__
       spinor: DO ISPINOR=0,W1%WDES1%NRSPINORS-1
          DO ISPINOR_=0,W1%WDES1%NRSPINORS-1
             DO NI=1,W1%WDES1%NIONS
                   NT=W1%WDES1%ITYP(NI)
                   LMMAXC=W1%WDES1%LMMAX(NT)
                   IF (LMMAXC==0) CYCLE

                   NPRO =ISPINOR *(W1%WDES1%NPRO/2)+W1%WDES1%LMBASE(NI)
                   NPRO_=ISPINOR_*(W1%WDES1%NPRO/2)+W1%WDES1%LMBASE(NI)
                   CALL ECCP_NL(SIZE(CQIJ,1),LMMAXC,CQIJ(1,1,NI,ISP+ISPINOR_+2*ISPINOR),W1%CPROJ(NPRO_+1),W2%CPROJ(NPRO+1),CADD)
             ENDDO
          ENDDO
       ENDDO spinor

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
       CSCPD=GREAL(CSCPD+CADD)
    ENDIF
!=======================================================================
! performe orthogonalisations
!=======================================================================
    CALLMPI( M_sum_z(W1%WDES1%COMM_INB, CSCPD, 1))

    CALL __ZAXPY__(W1%WDES1%NPL ,-CSCPD,W2%CW(1)   ,1,W1%CW(1)   ,1)
    IF (W1%WDES1%LGAMMA) THEN
       CALL __DAXPY__(W1%WDES1%NPRO,-CSCPD,W2%CPROJ(1),1,W1%CPROJ(1),1)
    ELSE
       CALL __ZAXPY__(W1%WDES1%NPRO,-CSCPD,W2%CPROJ(1),1,W1%CPROJ(1),1)
    ENDIF

    RETURN
  END SUBROUTINE PROJCN


!************************ SUBROUTINE W1_GATHER ************************
!
! This subroutine gathers a set of wavefunctions starting 
! from band NB1 until NB2 to all nodes
! 
!**********************************************************************

  SUBROUTINE W1_GATHER( W, NB1, NB2, ISP, W1)
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed

  ! local
    INTEGER :: NN, N, NLOC, NCPU

    NCPU=W%WDES%NB_PAR

    DO N=NB1,NB2
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       CALL W1_COPY( ELEMENT( W, W1(NN)%WDES1, N, ISP), W1(NN) )
       CALL FFTWAV_W1( W1(NN))
    ENDDO

    NLOC=(NB2-NB1+1)*W%WDES%NB_PAR

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO NN=1,NLOC
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CW(1), &
               SIZE(W1(NN)%CW),MOD(NN-1,NCPU)+1)
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(NN-1,NCPU)+1)
#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#endif
       ENDDO

#ifdef MPI_barrier_after_bcast
      CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif
    ENDIF
#endif
  END SUBROUTINE W1_GATHER

!
!  new version that uses MPI_Allgather
!  this seems to be even slower then the version above
!
  SUBROUTINE W1_GATHER_NEW( W, NB1, NB2, ISP, W1)
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed

    COMPLEX(q) :: CWBUFF(W%WDES%NRPLWV,W%WDES%NB_PAR)        ! receive buffer
    GDEF       :: CPROJBUFF(W%WDES%NPROD,W%WDES%NB_PAR)      ! receive buffer for projectors
    COMPLEX(q) :: CRBUFF(W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,W%WDES%NB_PAR) ! receive buffer for real space part
  ! local
    INTEGER :: NN, N, NLOC, NCPU
    INTEGER :: ierr          ! error for mpi

    NCPU=W%WDES%NB_PAR
    
    IF (SIZE(W1(1)%CW) /= W%WDES%NRPLWV) THEN
       CALL vtutor%bug("internal error in W1_GATHER: sizes inconsistent " // str(SIZE(W1(1)%CW)), __FILE__, __LINE__)
    ENDIF

    DO N=NB1,NB2
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       CALL W1_COPY( ELEMENT( W, W1(NN)%WDES1, N, ISP), W1(NN) )
       CALL FFTWAV_W1( W1(NN))

    ! distribute W1 to all nodes
#ifdef MPI
       IF (W%WDES%DO_REDIS) THEN
          CALL MPI_Allgather(W1(NN)%CR(1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS, MPI_double_complex, CRBUFF(1,1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS, & 
               MPI_double_complex, W%WDES%COMM_INTER%MPI_COMM, ierr)
          IF ( ierr /= MPI_success ) CALL vtutor%error('ERROR: M_Allgather returns' // str(ierror))

          CALL MPI_Allgather(W1(NN)%CW(1), W%WDES%NRPLWV , MPI_double_complex, CWBUFF(1,1), W%WDES%NRPLWV , & 
               MPI_double_complex, W%WDES%COMM_INTER%MPI_COMM, ierr)
          IF ( ierr /= MPI_success ) CALL vtutor%error('ERROR: M_Allgather returns' // str(ierror))

#ifdef gammareal
          IF (W%WDES%LOVERL) CALL MPI_Allgather(W1(NN)%CPROJ(1), W%WDES%NPROD, MPI_double_precision, &
               CPROJBUFF(1,1), W%WDES%NPROD, MPI_double_precision, W%WDES%COMM_INTER%MPI_COMM, ierr)
#else 
          IF (W%WDES%LOVERL) CALL MPI_Allgather(W1(NN)%CPROJ(1), W%WDES%NPROD, MPI_double_complex, &
               CPROJBUFF(1,1), W%WDES%NPROD, MPI_double_complex, W%WDES%COMM_INTER%MPI_COMM, ierr)
#endif
          IF ( ierr /= MPI_success ) CALL vtutor%error('ERROR: M_Allgather returns' // str(ierror))

          ! store results back in correct storage position
          DO NN=1,W%WDES%NB_PAR
             W1((N-NB1)*NCPU+NN)%CW(:)   =CWBUFF(:,NN)
             W1((N-NB1)*NCPU+NN)%CR(:)   =CRBUFF(:,NN)
             W1((N-NB1)*NCPU+NN)%CPROJ(:)=CPROJBUFF(:,NN)
          ENDDO
       END IF
#endif
    ENDDO
  END SUBROUTINE W1_GATHER_NEW

!************************ SUBROUTINE W1_GATHER ************************
!
! This subroutine gathers a set of wavefunctions starting 
! from band NB1 until NB2 to all nodes
! compared to the previous routine only up to NLOC bands are collected
! 
!**********************************************************************

  SUBROUTINE W1_GATHER_N( W, NB1, NB2, ISP, W1, NLOC)
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed
    INTEGER :: NLOC          ! total number of bands to be collected

  ! local
    INTEGER :: NN, N, NCPU

    NCPU=W%WDES%NB_PAR

    DO N=NB1,NB2
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       IF (NN>NLOC) EXIT
       CALL W1_COPY( ELEMENT( W, W1(NN)%WDES1, N, ISP), W1(NN) )
       CALL FFTWAV_W1( W1(NN))
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO NN=1,NLOC
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CW(1), &
               SIZE(W1(NN)%CW),MOD(NN-1,NCPU)+1)
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(NN-1,NCPU)+1)
#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#endif
       ENDDO

#ifdef MPI_barrier_after_bcast
      CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif

    ENDIF
#endif

  END SUBROUTINE W1_GATHER_N


!************************ SUBROUTINE W1_GATHER_GLB ********************
!
! This subroutine gathers a set of wavefunctions starting 
! from band NB1 until NB2 to all nodes
! compared to the previous version the global instead of local band 
! indices are supplied
! 
!**********************************************************************

  SUBROUTINE W1_GATHER_GLB( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed

  ! local
    INTEGER :: N_INTO_TOT, N, NCPU
    INTEGER :: ierror

    PROFILING_START('w1_gather_glb')

    NCPU=W%WDES%NB_PAR
    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(N_INTO_TOT-NB1+1)%WDES1, N, ISP), W1(N_INTO_TOT-NB1+1) )
          CALL FFTWAV_W1( W1(N_INTO_TOT-NB1+1))
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CW(1), &
               SIZE(W1(N)%CW),MOD(N_INTO_TOT-1,NCPU)+1)
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(N_INTO_TOT-1,NCPU)+1)
#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#endif
       ENDDO

#ifdef MPI_barrier_after_bcast
      CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif

    ENDIF
#endif
    PROFILING_STOP('w1_gather_glb')

  END SUBROUTINE W1_GATHER_GLB

  SUBROUTINE W1_GATHER_GLB_( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed

  ! local
    INTEGER :: N_INTO_TOT, N, NB_LOCAL

    PROFILING_START('w1_gather_glb_')

    DO N=NB1,NB2
       N_INTO_TOT=N-NB1+1
       IF (MOD(N-1,W%WDES%NB_PAR)+1==W%WDES%NB_LOW) THEN
          NB_LOCAL=1+(N-1)/W%WDES%NB_PAR
          CALL W1_COPY(ELEMENT(W,W1(N_INTO_TOT)%WDES1,NB_LOCAL,ISP),W1(N_INTO_TOT))
          CALL FFTWAV_W1(W1(N_INTO_TOT))
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO N=NB1,NB2
          N_INTO_TOT=N-NB1+1
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N_INTO_TOT)%CW(1), &
               SIZE(W1(N_INTO_TOT)%CW),MOD(N-1,W%WDES%NB_PAR)+1)         
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N_INTO_TOT)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(N-1,W%WDES%NB_PAR)+1)
#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N_INTO_TOT)%CPROJ(1), &
               W%WDES%NPROD,MOD(N-1,W%WDES%NB_PAR)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(N_INTO_TOT)%CPROJ(1), &
               W%WDES%NPROD,MOD(N-1,W%WDES%NB_PAR)+1)
#endif
       ENDDO
#ifdef MPI_barrier_after_bcast
       CALLMPI( M_barrier(W%WDES%COMM_INTER))
#endif
    ENDIF
#endif
    PROFILING_STOP('w1_gather_glb_')

  END SUBROUTINE W1_GATHER_GLB_


!************************ SUBROUTINE W1_GATHER_GLB_NOFFT ********************
! identical to W1_GATHER_GLB but does not FFT to real space 
!**********************************************************************

  SUBROUTINE W1_GATHER_GLB_NOFFT( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed

  ! local
    INTEGER :: N_INTO_TOT, N, NCPU
    INTEGER :: ierror

    PROFILING_START('w1_gather_glb_nofft')

    NCPU=W%WDES%NB_PAR
    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(N_INTO_TOT-NB1+1)%WDES1, N, ISP), W1(N_INTO_TOT-NB1+1) )
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CW(1), &
               SIZE(W1(N)%CW),MOD(N_INTO_TOT-1,NCPU)+1)
#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#endif
       ENDDO

#ifdef MPI_barrier_after_bcast
      CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif

    ENDIF
#endif
    PROFILING_STOP('w1_gather_glb_nofft')

  END SUBROUTINE W1_GATHER_GLB_NOFFT


!************************ SUBROUTINE W1_IGATHER_GLB *******************
!
! Gathers a set of wavefunctions starting from band NB1 until NB2
! (global band indices) to all nodes using non-blocking bcast_from.
!
! In case the code is compiled with -Dshmem_bcast_buffer communication
! will be between ranks within COMM_inter_node, i.e., only those ranks
! within COMM_INTER that are NOT on the same physical node will talk to
! eachother. MPI communication between ranks within COMM_INTER that
! reside on the same node is not necessary since they access a common
! shared memory segment.
!
! In case the code is NOT compiled with -Dshmem_bcast_buffer, the
! communication will be between all ranks in COMM_INTER.
! 
!**********************************************************************
#ifdef shmem_bcast_buffer
  SUBROUTINE W1_IGATHER_GLB( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed

  ! local
    INTEGER :: NI, N, NB_LOCAL

    INTEGER :: IDO(NB2-NB1+1)

    INTEGER :: nrequests,requests(3*(NB2-NB1+1))

    PROFILING_START('w1_igather_glb')

    CALL M_barrier(W%WDES%COMM_intra_node)

    DO N=NB1,NB2
       NI=N-NB1+1
       IF (MOD(N-1,W%WDES%NB_PAR)+1==W%WDES%NB_LOW) THEN
          NB_LOCAL=1+(N-1)/W%WDES%NB_PAR
          CALL W1_COPY(ELEMENT(W,W1(NI)%WDES1,NB_LOCAL,ISP),W1(NI))
          CALL FFTWAV_W1(W1(NI))
          IDO(NI)=W%WDES%COMM_inter_node%NODE_ME
       ELSE
          IDO(NI)=0
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%COMM_inter_node%NCPU>1) THEN

       PROFILING_START('stage_1')

       CALLMPI( M_sum_i(W%WDES%COMM_inter_node,IDO,NB2-NB1+1))

       nrequests=0
       DO N=NB1,NB2
          NI=N-NB1+1
          IF (IDO(NI)>0) THEN
             nrequests=nrequests+1
             CALL M_ibcast_z_from(W%WDES%COMM_inter_node,W1(NI)%CW(1), &
            &     SIZE(W1(NI)%CW),IDO(NI),requests(nrequests))
             nrequests=nrequests+1
             CALL M_ibcast_z_from(W%WDES%COMM_inter_node,W1(NI)%CR(1), &
            &     SIZE(W1(NI)%CR),IDO(NI),requests(nrequests))
             IF (W%WDES%LOVERL) THEN
                nrequests=nrequests+1
#ifndef gammareal
                CALL M_ibcast_z_from(W%WDES%COMM_inter_node,W1(NI)%CPROJ(1), &
               &     SIZE(W1(NI)%CPROJ),IDO(NI),requests(nrequests))
#else
                CALL M_ibcast_d_from(W%WDES%COMM_inter_node,W1(NI)%CPROJ(1), &
               &     SIZE(W1(NI)%CPROJ),IDO(NI),requests(nrequests))
#endif
             ENDIF
          ENDIF
       ENDDO

       CALL M_waitall(nrequests,requests(1))

       PROFILING_STOP('stage_1')

    ENDIF

    CALL M_barrier(W%WDES%COMM_intra_node)
!!#ifndef shmem_bcast_buffer
!!!   CALL M_barrier(W%WDES%COMM_INTER)
!!
!!    IF (W%WDES%COMM_intra_node%NCPU>1) THEN
!!       DO N=1,NB2-NB1+1
!!          IF (IDO(N)>0) IDO(N)=W%WDES%COMM_intra_node%NODE_ME
!!       ENDDO
!!
!!       CALLMPI( M_sum_i(W%WDES%COMM_intra_node,IDO,NB2-NB1+1))
!!
!!       PROFILING_START('stage_2')
!!
!!       nrequests=0
!!       DO N=NB1,NB2
!!          NI=N-NB1+1
!!          IF (IDO(NI)>0) THEN
!!             nrequests=nrequests+1
!!             CALL M_ibcast_z_from(W%WDES%COMM_intra_node,W1(NI)%CW(1), &
!!            &     SIZE(W1(NI)%CW),IDO(NI),requests(nrequests))
!!             nrequests=nrequests+1
!!             CALL M_ibcast_z_from(W%WDES%COMM_intra_node,W1(NI)%CR(1), &
!!            &     SIZE(W1(NI)%CR),IDO(NI),requests(nrequests))
!!             IF (W%WDES%LOVERL) THEN
!!                nrequests=nrequests+1
!!#ifndef gammareal
!!                CALL M_ibcast_z_from(W%WDES%COMM_intra_node,W1(NI)%CPROJ(1), &
!!               &     SIZE(W1(NI)%CPROJ),IDO(NI),requests(nrequests))
!!#else
!!                CALL M_ibcast_d_from(W%WDES%COMM_intra_node,W1(NI)%CPROJ(1), &
!!               &     SIZE(W1(NI)%CPROJ),IDO(NI),requests(nrequests))
!!#endif
!!             ENDIF
!!          ENDIF
!!       ENDDO
!!
!!       CALL M_waitall(nrequests,requests(1))
!!
!!       PROFILING_STOP('stage_2')
!!
!!    ENDIF
!!#endif
#endif
    PROFILING_STOP('w1_igather_glb')

  END SUBROUTINE W1_IGATHER_GLB
#else
  SUBROUTINE W1_IGATHER_GLB( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed

  ! local
    INTEGER :: NI, N, NB_LOCAL
    INTEGER :: ierror

    INTEGER :: nrequests,requests(3*(NB2-NB1+1))

    PROFILING_START('w1_igather_glb')

    DO N=NB1,NB2
       NI=N-NB1+1
       IF (MOD(N-1,W%WDES%NB_PAR)+1==W%WDES%NB_LOW) THEN
          NB_LOCAL=1+(N-1)/W%WDES%NB_PAR
          CALL W1_COPY(ELEMENT(W,W1(NI)%WDES1,NB_LOCAL,ISP),W1(NI))
          CALL FFTWAV_W1(W1(NI))
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%COMM_INTER%NCPU>1) THEN

       nrequests=0
       DO N=NB1,NB2
          NI=N-NB1+1

          nrequests=nrequests+1
          CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CW(1), &
         &     SIZE(W1(NI)%CW),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
          nrequests=nrequests+1
          CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CR(1), &
         &     SIZE(W1(NI)%CR),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
          IF (W%WDES%LOVERL) THEN
             nrequests=nrequests+1
#ifndef gammareal
             CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CPROJ(1), &
            &     SIZE(W1(NI)%CPROJ),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
#else
             CALL M_ibcast_d_from(W%WDES%COMM_INTER,W1(NI)%CPROJ(1), &
            &     SIZE(W1(NI)%CPROJ),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
#endif
          ENDIF
       ENDDO

       CALL M_waitall(nrequests,requests(1))

    ENDIF
#endif
    PROFILING_STOP('w1_igather_glb')

  END SUBROUTINE W1_IGATHER_GLB
#endif

!************************ SUBROUTINE W1_REDUCE_GLB ********************
!
! Take the sum of W1(i)%CW, W1(i)%CPROJ, for i=1,NB2-NB1+1 over all
! ranks in COMM_inter_node, using reduce_to.
! The result is reduced onto the MPI-rank in COMM_inter_node that would
! normally own (part of) the band with the global index N=NB1,...,NB2.
!
!**********************************************************************

  SUBROUTINE W1_REDUCE_GLB(WDES1,W1,NB1,NB2)
    IMPLICIT NONE
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1(:) ! array into which the merge is performed
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band

  ! local
    INTEGER :: N_INTO_TOT, N
    INTEGER :: ierror

    INTEGER :: IDO(NB2-NB1+1)

    PROFILING_START('w1_reduce_glb')
#ifdef MPI
    DO N=NB1,NB2
       N_INTO_TOT=N-NB1+1
       IF (MOD(N-1,WDES1%NB_PAR)+1==WDES1%NB_LOW) THEN
          IDO(N_INTO_TOT)=WDES1%COMM_inter_node%NODE_ME
       ELSE
          IDO(N_INTO_TOT)=0
       ENDIF
    ENDDO

    CALLMPI( M_sum_i(WDES1%COMM_inter_node,IDO,NB2-NB1+1))

    ! reduce W1
    IF (WDES1%DO_REDIS) THEN

       DO N=NB1,NB2
          N_INTO_TOT=N-NB1+1
          IF (IDO(N_INTO_TOT)>0) THEN
             CALL M_reduce_z_to(WDES1%COMM_inter_node,W1(N_INTO_TOT)%CW(1), &
            &     SIZE(W1(N_INTO_TOT)%CW),IDO(N_INTO_TOT))
#ifndef gammareal
             IF (WDES1%LOVERL) &
                CALL M_reduce_z_to(WDES1%COMM_inter_node,W1(N_INTO_TOT)%CPROJ(1), &
            &     SIZE(W1(N_INTO_TOT)%CPROJ),IDO(N_INTO_TOT))
#else
             IF (WDES1%LOVERL) &
                CALL M_reduce_d_to(WDES1%COMM_inter_node,W1(N_INTO_TOT)%CPROJ(1), &
            &     SIZE(W1(N_INTO_TOT)%CPROJ),IDO(N_INTO_TOT))
#endif
          ENDIF
       ENDDO

    ENDIF
#endif
    PROFILING_STOP('w1_reduce_glb')

  END SUBROUTINE W1_REDUCE_GLB


!************************ SUBROUTINE W1_IREDUCE_GLB *******************
!
! Take the sum of W1(i)%CW, W1(i)%CPROJ, for i=1,NB2-NB1+1 over all
! ranks in COMM_INTER, using non-blocking reduce_to.
! The result is reduced onto the MPI-rank in COMM_INTER that would
! normally own (part of) the band with the global index N=NB1,...,NB2.
!
!**********************************************************************

  SUBROUTINE W1_IREDUCE_GLB(WDES1,W1,NB1,NB2)
    IMPLICIT NONE
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1(:) ! array into which the merge is performed
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band

  ! local
    INTEGER :: NI, N
    INTEGER :: requests(2*(NB2-NB1+1))
    INTEGER :: nrequests

    PROFILING_START('w1_ireduce_glb')

    ! reduce W1
#ifdef MPI
    IF (WDES1%DO_REDIS) THEN

       nrequests=0
       DO N=NB1,NB2
          NI=N-NB1+1

          nrequests=nrequests+1
          CALL M_ireduce_z_to(WDES1%COMM_INTER,W1(NI)%CW(1), &
         &     SIZE(W1(NI)%CW),MOD(N-1,WDES1%NB_PAR)+1,requests(nrequests))

          IF (WDES1%LOVERL) THEN
             nrequests=nrequests+1
#ifndef gammareal
             CALL M_ireduce_z_to(WDES1%COMM_INTER,W1(NI)%CPROJ(1), &
            &     SIZE(W1(NI)%CPROJ),MOD(N-1,WDES1%NB_PAR)+1,requests(nrequests))
#else
             CALL M_ireduce_d_to(WDES1%COMM_INTER,W1(NI)%CPROJ(1), &
            &     SIZE(W1(NI)%CPROJ),MOD(N-1,WDES1%NB_PAR)+1,requests(nrequests))
#endif
          ENDIF
       ENDDO

       CALL M_waitall(nrequests,requests(1))

    ENDIF
#endif
    PROFILING_STOP('w1_ireduce_glb')

  END SUBROUTINE W1_IREDUCE_GLB

!************************ SUBROUTINE W1_GATHER_GLB_NOCR ***************
!
! This subroutine gathers a set of wavefunctions starting 
! from band NB1 until NB2 to all nodes
! compared to the previous version the global instead of local band 
! indices are supplied and the real space part is not
! communicated
! 
!**********************************************************************

  SUBROUTINE W1_GATHER_GLB_NOCR( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed

  ! local
    INTEGER :: N_INTO_TOT, N, NCPU

    PROFILING_START('w1_gather_glb_nocr')

    NCPU=W%WDES%NB_PAR
    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(N_INTO_TOT-NB1+1)%WDES1, N, ISP), W1(N_INTO_TOT-NB1+1) )
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CW(1), &
               SIZE(W1(N)%CW),MOD(N_INTO_TOT-1,NCPU)+1)
#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#endif
       ENDDO
    ENDIF
#endif
    PROFILING_STOP('w1_gather_glb_nocr')

  END SUBROUTINE W1_GATHER_GLB_NOCR

#ifdef use_shmem
!************************ SUBROUTINE W1_GATHER_GLB_NOCR_SHMEM *********
!
! This subroutine gathers a set of wavefunctions starting 
! from band NB1 until NB2 to all nodes
! compared to the previous version the global instead of local band 
! indices are supplied
! 
!**********************************************************************

  SUBROUTINE W1_GATHER_GLB_NOCR_SHMEM(W,NB1,NB2,ISP,W1)
    IMPLICIT NONE
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1) :: W1(:) ! array into which the merge is performed
    ! local variables
    COMPLEX(q) :: CWORK(W1(1)%WDES1%NRPLWV) 
    INTEGER :: N_INTO_TOT,N,NCPU

    NCPU=W%WDES%NB_PAR
    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(N_INTO_TOT-NB1+1)%WDES1, N, ISP), W1(N_INTO_TOT-NB1+1) )
!         CALL FFTWAV_W1( W1(N_INTO_TOT-NB1+1))
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          ! the node that owns this band stores it onto CWORK
          IF (MOD(N_INTO_TOT-1,NCPU)+1==W%WDES%COMM_INTER%NODE_ME) CWORK=W1(N)%CW
          ! and broadcasts CWORK to all other nodes
          CALL M_bcast_z_from(W%WDES%COMM_INTER,CWORK(1), &
               SIZE(CWORK),MOD(N_INTO_TOT-1,NCPU)+1)
          ! but only the root nodes in COMM_SHMEM need to store CWORK back into W1
          IF (MOD(N_INTO_TOT-1,NCPU)+1/=W%WDES%COMM_INTER%NODE_ME &
             .AND.W%WDES%COMM_SHMEM%NODE_ME==1) W1(N)%CW=CWORK
#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(N)%CPROJ(1), &
               W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#endif
       ENDDO
    ENDIF
#endif

  END SUBROUTINE W1_GATHER_GLB_NOCR_SHMEM


!************************ SUBROUTINE W1_GATHER_GLB_ALLK_SHMEM *********
!
! Gathers a set of wavefunctions including band NB1 to NB2 over
! all k-points to all nodes using shared memory
!
!**********************************************************************

  SUBROUTINE W1_GATHER_GLB_ALLK_SHMEM( W, NB1, NB2, ISP, WF)
    IMPLICIT NONE
    TYPE (wavespin) W          ! wavefunction
    INTEGER :: NB1             ! starting band
    INTEGER :: NB2             ! final band
    INTEGER :: ISP             ! spin
    TYPE (wavefun1) :: WF(:,:) ! array into which the merge is performed
    ! local variables
    TYPE(wavefun1):: WAUX  
    TYPE(wavedes1), TARGET :: WDESAUX
    COMPLEX(q) :: CWORK(W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)
    GDEF CPROJ(W%WDES%NPROD)
    INTEGER :: N_INTO_TOT,N,NCPU,IK,NKPTS

    CALL SETWDES(W%WDES,WDESAUX,0)
    CALL NEWWAV(WAUX,WDESAUX,.TRUE.)
    NCPU=W%WDES%NB_PAR
    NKPTS=W%WDES%NKPTS

    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          DO IK=1,NKPTS
             ! set wave descriptor to current k-point for FFT
             CALL SETWDES(W%WDES,WDESAUX,IK)
             CALL W1_COPY_NOCR( ELEMENT( W, WDESAUX, N, ISP), WAUX)
             CALL FFTWAV_W1(WAUX)
             CALL ZCOPY( WDESAUX%GRID%MPLWV*WDESAUX%NRSPINORS, WAUX%CR(1), 1, WF(N_INTO_TOT-NB1+1,IK)%CR(1), 1)
             IF (W%WDES%LGAMMA) THEN
                CALL DCOPY( W%WDES%NPROD, WAUX%CPROJ(1), 1,  WF(N_INTO_TOT-NB1+1,IK)%CPROJ(1), 1)
             ELSE
                CALL ZCOPY( W%WDES%NPROD, WAUX%CPROJ(1), 1,  WF(N_INTO_TOT-NB1+1,IK)%CPROJ(1), 1)
             ENDIF
          END DO
       ENDIF
    ENDDO
    CALL DELWAV(WAUX,.TRUE.)
    ! distribute WF to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          DO IK=1,NKPTS
             ! the node that owns this band stores it onto CWORK
             IF (MOD(N_INTO_TOT-1,NCPU)+1==W%WDES%COMM_INTER%NODE_ME) THEN
                CWORK(1:SIZE(WF(N,IK)%CR))=WF(N,IK)%CR(:); CPROJ(1:SIZE(WF(N,IK)%CPROJ))=WF(N,IK)%CPROJ(:)
             ENDIF
             ! and broadcasts CWORK to all other nodes
             CALL M_bcast_z_from(W%WDES%COMM_INTER,CWORK(1),SIZE(CWORK),MOD(N_INTO_TOT-1,NCPU)+1)
             IF (W%WDES%LOVERL) THEN
#ifdef gammareal
                CALL M_bcast_d_from(W%WDES%COMM_INTER,CPROJ(1),SIZE(CPROJ),MOD(N_INTO_TOT-1,NCPU)+1)
#else
                CALL M_bcast_z_from(W%WDES%COMM_INTER,CPROJ(1),SIZE(CPROJ),MOD(N_INTO_TOT-1,NCPU)+1)         
#endif
             ENDIF
             ! but only the root nodes in COMM_SHMEM need to store CWORK back into WF
             IF (MOD(N_INTO_TOT-1,NCPU)+1/=W%WDES%COMM_INTER%NODE_ME.AND.W%WDES%COMM_SHMEM%NODE_ME==1) THEN
                WF(N,IK)%CR=CWORK(1:SIZE(WF(N,IK)%CR)); WF(N,IK)%CPROJ(:)=CPROJ(1:SIZE(WF(N,IK)%CPROJ))
             ENDIF
          ENDDO
       ENDDO
#ifdef MPI_barrier_after_bcast
       CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif
    ENDIF
#endif
  END SUBROUTINE W1_GATHER_GLB_ALLK_SHMEM
#endif

!************************ SUBROUTINE W1_GATHER_GLB_ALLK ***************
!
! Gathers a set of wavefunctions including band NB1 to NB2 over
! all k-points to all nodes
!
!**********************************************************************

  SUBROUTINE W1_GATHER_GLB_ALLK( W, NB1, NB2, ISP, WF)
    IMPLICIT NONE
    TYPE (wavespin) W          ! wavefunction
    INTEGER :: NB1             ! starting band
    INTEGER :: NB2             ! final band
    INTEGER :: ISP             ! spin
    TYPE (wavefun1) :: WF(:,:) ! array into which the merge is performed
    ! local
    TYPE (wavefun1):: WAUX  
    TYPE(wavedes1), TARGET :: WDESAUX
    INTEGER :: N_INTO_TOT,N,NCPU,IK,NKPTS


    CALL SETWDES(W%WDES,WDESAUX,0)
    CALL NEWWAV(WAUX,WDESAUX,.TRUE.)
    NCPU=W%WDES%NB_PAR
    NKPTS=W%WDES%NKPTS

    ! W1 contains bands NB1 to NB2 for all k-points
    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          DO IK=1,NKPTS
             ! set wave descriptor to current k-point for FFT
             CALL SETWDES(W%WDES,WDESAUX,IK)
             CALL W1_COPY_NOCR( ELEMENT( W, WDESAUX, N, ISP), WAUX)
             CALL FFTWAV_W1(WAUX)
             CALL ZCOPY( WDESAUX%GRID%MPLWV*WDESAUX%NRSPINORS, WAUX%CR(1), 1, WF(N_INTO_TOT-NB1+1,IK)%CR(1), 1)
             IF (W%WDES%LGAMMA) THEN
                CALL DCOPY( W%WDES%NPROD, WAUX%CPROJ(1), 1,  WF(N_INTO_TOT-NB1+1,IK)%CPROJ(1), 1)
             ELSE
                CALL ZCOPY( W%WDES%NPROD, WAUX%CPROJ(1), 1,  WF(N_INTO_TOT-NB1+1,IK)%CPROJ(1), 1)
             ENDIF
          END DO
       ENDIF
    ENDDO
    CALL DELWAV(WAUX,.TRUE.)
    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO N_INTO_TOT=NB1,NB2
          N=N_INTO_TOT-NB1+1
          DO IK=1,NKPTS
             CALL M_bcast_z_from(W%WDES%COMM_INTER,WF(N,IK)%CR(1), &
                  W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(N_INTO_TOT-1,NCPU)+1)
#ifndef gammareal
             IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,WF(N,IK)%CPROJ(1), &
                  W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#else
             IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,WF(N,IK)%CPROJ(1), &
                  W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#endif
          ENDDO
       END DO

#ifdef MPI_barrier_after_bcast
       CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif

    ENDIF
#endif
  END SUBROUTINE W1_GATHER_GLB_ALLK


!************************ SUBROUTINE W1_GATHER ************************
!
! This subroutine gathers a set of orbitals
! from band NB1 until NB2 to all nodes
! compared to previous versions the set in collected
! into work arrays 
! 
!**********************************************************************

  SUBROUTINE W1_GATHER_ARRAY( W, NB1, NB2, ISP, W1, CR, CPROJ)
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1     ! array for FFT 
    GDEF    :: CR(:,:)       ! collected real space orbitals
    GDEF    :: CPROJ(:,:)    ! collected projected orbitals

  ! local
    INTEGER :: NN, N, NLOC, NCPU

    IF (W%WDES%NRSPINORS/=1) THEN
       CALL vtutor%bug("internal error in W1_GATHER_ARRAY: at the moment spinors are not supported", __FILE__, __LINE__)
    ENDIF

    NCPU=W%WDES%NB_PAR

    DO N=NB1,NB2
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       CALL W1_COPY( ELEMENT( W, W1%WDES1, N, ISP), W1 )
       CALL FFTWAV_W1( W1)

       IF (NN> SIZE(CR,2) .OR. NN >SIZE(CPROJ,2)) THEN
          CALL vtutor%bug("internal error in W1_GATHER_ARRAY: bound exceed " // str(NN) // " " // &
             str(SIZE(CR,2)) // " " // str(SIZE(CPROJ,2)), __FILE__, __LINE__)
       ENDIF
       IF (SIZE(CR,1)/=W1%WDES1%GRID%NPLWV) THEN
          WRITE(*,*) 'internal error in W1_GATHER_ARRAY: size mismatch ',SIZE(CR,1),W1%WDES1%GRID%NPLWV
       ENDIF

       CR(1:W1%WDES1%GRID%RL%NP, NN)=W1%CR(1:W1%WDES1%GRID%RL%NP)
       ! pad CR with zeros (just in case)
       CR(W1%WDES1%GRID%RL%NP+1:SIZE(CR,1), NN)=0
       CPROJ(:,NN)=W1%CPROJ(:)

    ! distribute W1 to all nodes
#ifdef MPI
       IF (W%WDES%DO_REDIS) THEN
        DO NN=1,NCPU
#ifndef gammareal
          CALL M_bcast_z_from(W%WDES%COMM_INTER,CR(1,NN+(N-NB1)*NCPU), &
               SIZE(CR,1),MOD(NN-1,NCPU)+1)
#else
          CALL M_bcast_d_from(W%WDES%COMM_INTER,CR(1,NN+(N-NB1)*NCPU), &
               SIZE(CR,1),MOD(NN-1,NCPU)+1)
#endif

#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,CPROJ(1,NN+(N-NB1)*NCPU), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,CPROJ(1,NN+(N-NB1)*NCPU), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#endif
         ENDDO
       END IF
#endif
      ENDDO

  END SUBROUTINE W1_GATHER_ARRAY


!************************ SUBROUTINE W1_GATHER ************************
!
! This subroutine gathers a set of orbitals
! from band NB1 until NB2 to all nodes
! this version collects the place wave coefficients and the
! CPROJ coefficients
! 
!**********************************************************************

  SUBROUTINE W1_GATHER_ARRAY_RECIPROCAL( W, NB1, NB2, ISP, W1, CG, CPROJ)
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1     ! array for FFT 
    COMPLEX(q) :: CG(:,:)    ! collected real space orbitals
    GDEF       :: CPROJ(:,:) ! collected projected orbitals

    COMPLEX(q) :: CWBUFF(SIZE(CG,1))        ! send buffer
    GDEF       :: CPROJBUFF(SIZE(CPROJ,1))  ! send buffer for projectors

  ! local
    INTEGER :: NN, N, NLOC, NCPU
    INTEGER :: ierr          ! error for mpi

    IF (W%WDES%NRSPINORS/=1) THEN
       CALL vtutor%bug("internal error in W1_GATHER_ARRAY: at the moment spinors are not supported", __FILE__, __LINE__)
    ENDIF
    IF (SIZE(CPROJ,1) > SIZE(W1%CPROJ,1)) THEN
       CALL vtutor%bug("internal error in W1_GATHER_ARRAY_RECIPROCAL: CPROJ size inconsistent " // &
          str(SIZE(CPROJ,1)) // " " // str(SIZE(W1%CPROJ,1)), __FILE__, __LINE__)
    ENDIF

    NCPU=W%WDES%NB_PAR

    ! NB1 and NB2 are local indices for collecting the bands
    DO N=NB1,NB2
       !NB_LOW is the off-set of the local node
       !NN is index into CG array
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       !this copies orbital N (reciprocal and real (if allocated) part to W1)
       !these are local indices
       CALL W1_COPY( ELEMENT( W, W1%WDES1, N, ISP), W1 )

       IF (NN> SIZE(CG,2) .OR. NN >SIZE(CPROJ,2)) THEN
          CALL vtutor%bug("internal error in W1_GATHER_ARRAY_RECIPROCAL: bound exceed " // str(NN) &
             // " " // str(SIZE(CG,2)) // " " // str(SIZE(CPROJ,2)), __FILE__, __LINE__)
       ENDIF
  
       ! copy data over to return arrays: CG and CPROJ
       CG(1:W1%WDES1%NPL ,NN)=W1%CW(1:W1%WDES1%NPL)
       CG(W1%WDES1%NPL+1:,NN)=0.0_q                 ! pad with zero

       CPROJ(:,NN)=W1%CPROJ(1:SIZE(CPROJ,1))
      ! distribute W1 to all nodes, CG(1,(N-NB)*NCPU)
#ifdef MPI
       ! MPI_IN_PLACE was tested but found to be much slower (code has been removed)
       IF (W%WDES%DO_REDIS) THEN
          CPROJBUFF(:)=CPROJ(:,NN)
          CWBUFF(:)   =CG(:,NN)
          CALL MPI_Allgather(CWBUFF(1), SIZE(CWBUFF), MPI_double_complex, CG(1,(N-NB1)*NCPU+1), SIZE(CG,1), & 
               MPI_double_complex, W%WDES%COMM_INTER%MPI_COMM, ierr)
          IF ( ierr /= MPI_success ) CALL vtutor%error('ERROR: M_Allgather returns' // str(ierror))
#ifdef gammareal
          IF (W%WDES%LOVERL) CALL MPI_Allgather(CPROJBUFF(1), SIZE(CPROJBUFF), MPI_double_precision, &
               CPROJ(1,(N-NB1)*NCPU+1), SIZE(CPROJ,1), MPI_double_precision, W%WDES%COMM_INTER%MPI_COMM, ierr)
#else 
          IF (W%WDES%LOVERL) CALL MPI_Allgather(CPROJBUFF(1), SIZE(CPROJBUFF), MPI_double_complex,&
               CPROJ(1,(N-NB1)*NCPU+1), SIZE(CPROJ,1), MPI_double_complex, W%WDES%COMM_INTER%MPI_COMM, ierr)
#endif
          IF ( ierr /= MPI_success ) CALL vtutor%error('ERROR: M_Allgather returns' // str(ierror))
       END IF
#endif
    ENDDO

  END SUBROUTINE W1_GATHER_ARRAY_RECIPROCAL

!
! old version, uses bcast and is much slower
!
  SUBROUTINE W1_GATHER_ARRAY_RECIPROCAL_OLD( W, NB1, NB2, ISP, W1, CG, CPROJ)
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1     ! array for FFT 
    COMPLEX(q) :: CG(:,:)    ! collected real space orbitals
    GDEF       :: CPROJ(:,:) ! collected projected orbitals

  ! local
    INTEGER :: NN, N, NLOC, NCPU

    IF (W%WDES%NRSPINORS/=1) THEN
       CALL vtutor%bug("internal error in W1_GATHER_ARRAY: at the moment spinors are not supported", __FILE__, __LINE__)
    ENDIF

    NCPU=W%WDES%NB_PAR

    DO N=NB1,NB2
       !NB_LOW will be the off-set of the local node
       NN=(N-NB1)*NCPU+W%WDES%NB_LOW
       !this copies orbital N (reciprocal and real (if allocated) part to W1)
       !these are local indices
       CALL W1_COPY( ELEMENT( W, W1%WDES1, N, ISP), W1 )

       IF (NN> SIZE(CG,2) .OR. NN >SIZE(CPROJ,2)) THEN
          CALL vtutor%bug("internal error in W1_GATHER_ARRAY: bound exceed " // str(NN) // " " // &
             str(SIZE(CG,2)) // " " // str(SIZE(CPROJ,2)), __FILE__, __LINE__)
       ENDIF
  
       !copy the local data to the position with global index in the strip
       CG(1:W1%WDES1%NPL, NN)=W1%CW(1:W1%WDES1%NPL)
       CG(W1%WDES1%NPL+1: SIZE(CG,1), NN)=0
       !the first index CPROJ size can be smaller than the second one (less cpus/tau than total number)
       !but this copies only the data for valid indices in the first array (or?)
       !!!test
       !IF (SIZE(CPROJ,1).LT.SIZE(W1%CPROJ,1)) WRITE (*,*) NN, W1%CPROJ(SIZE(CPROJ,1)+1:)
       IF (SIZE(CPROJ,1).GT.SIZE(W1%CPROJ,1)) THEN
         WRITE (*,*) 'internal error in W1_GATHER_ARRAY: Sizes inconsistent ', SIZE(CPROJ,1), SIZE(W1%CPROJ,1)
       ENDIF
       CPROJ(:,NN)=W1%CPROJ(1:SIZE(CPROJ,1))
       !IF (N.eq.NB1) WRITE (*,*) 'in gather array', SIZE(CPROJ,1), SIZE(CPROJ,2), SIZE(W1%CPROJ,1)
       !!!test

    ! distribute W1 to all nodes
#ifdef MPI
       IF (W%WDES%DO_REDIS) THEN
       DO NN=1,NCPU
          !so broadcast the data from node holding index NN to all other nodes 
          !I wonder how efficient this is, might be all right 10000 coefficients is 160kB of data, 
          !a bit on the low side but fine
          CALL M_bcast_z_from(W%WDES%COMM_INTER,CG(1,NN+(N-NB1)*NCPU), &
               SIZE(CG,1),MOD(NN-1,NCPU)+1)
#ifndef gammareal
        ! OK, the point is that the W%WDES%NPROD array can be longer than the size of CPROJ
        ! since both are a multiple of NCPU but once the total number and in the second case per tau
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,CPROJ(1,NN+(N-NB1)*NCPU), &
               SIZE(CPROJ,1),MOD(NN-1,NCPU)+1)

#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,CPROJ(1,NN+(N-NB1)*NCPU), &
               SIZE(CPROJ,1),MOD(NN-1,NCPU)+1)
#endif
       ENDDO
       END IF
#endif
    ENDDO

  END SUBROUTINE W1_GATHER_ARRAY_RECIPROCAL_OLD


!************************ SUBROUTINE W1_GATHER_DISTR ******************
!
! This subroutine gathers a set of wavefunctions starting 
! from band NB1 until NB2 to all nodes
! compared to the previous version the global instead of local band 
! indices are supplied
! 
!**********************************************************************

  SUBROUTINE W1_GATHER_DISTR( W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed
    ! local variables
    INTEGER :: N_INTO_TOT, N, NCPU
    INTEGER :: NBMIN, NBMAX
    TYPE (wavefun1):: WTMP
    
    NCPU=W%WDES%NB_PAR
    
    CALL NEWWAV(WTMP, W1(1)%WDES1, .FALSE.)
        
    ! establish global band index interval
    NBMIN=NB1
    NBMAX=NB2
#ifdef MPI
    DO N=1,NCPU
       CALL M_bcast_i_from(W%WDES%COMM_INTER, NBMIN, 1, n)
       CALL M_bcast_i_from(W%WDES%COMM_INTER, NBMAX, 1, n)
       IF (NB1<NBMIN) NBMIN=NB1
       IF (NB2>NBMAX) NBMAX=NB2
    ENDDO
#endif
    
    DO N_INTO_TOT=NBMIN,NBMAX
    
       ! if band resides on this node copy it to WTMP
       IF (MOD(N_INTO_TOT-W%WDES%NB_LOW,NCPU)==0) THEN
          N=(N_INTO_TOT-W%WDES%NB_LOW)/NCPU+1
          CALL W1_COPY( ELEMENT( W, WTMP%WDES1, N, ISP), WTMP )
       ENDIF

#ifdef MPI
       ! broadcast WTMP from the node where it resides
       CALL M_bcast_z_from(W%WDES%COMM_INTER,WTMP%CW(1), &
            SIZE(WTMP%CW),MOD(N_INTO_TOT-1,NCPU)+1)
#ifndef gammareal
       IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,WTMP%CPROJ(1), &
            W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#else
       IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,WTMP%CPROJ(1), &
            W%WDES%NPROD,MOD(N_INTO_TOT-1,NCPU)+1)
#endif
#endif
       ! if this band is targeted to reside on this node
       ! copy WTMP to W1(N_INTO_TOT-NB1+1)
       IF (N_INTO_TOT>=NB1 .AND. N_INTO_TOT<=NB2) THEN
          N=N_INTO_TOT-NB1+1
          CALL W1_COPY(WTMP,W1(N))
       ENDIF
    ENDDO
    
    ! FFT to real space
    DO N=1,(NB2-NB1)+1
       CALL FFTWAV_W1(W1(N))
    ENDDO

    CALL DELWAV(WTMP, .FALSE.)

  END SUBROUTINE W1_GATHER_DISTR


!************************ SUBROUTINE W1_GATHER_KSEL *******************
!
! This subroutine gathers a set of wavefunctions starting 
! from band NB1 until NB2 and distributes the data over k in a round
! robin fashion
! the data distribution is based on an index (k-point index) supplied as
! the last argument
!
!**********************************************************************

  SUBROUTINE W1_GATHER_KSEL( W, NB1, NB2, ISP, W1, NK)
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed
    INTEGER :: NK            ! k-points intex

   ! local
    TYPE (wavefun1):: WTMP   ! temporary for FFT
    INTEGER :: NN, N, NLOC, NCPU

    NCPU=W%WDES%NB_PAR

    CALL NEWWAV(WTMP, W1(1)%WDES1, .TRUE.)

    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(1)%WDES1, N, ISP), WTMP )
          CALL FFTWAV_W1( WTMP)
       ENDIF
#ifdef MPI
       IF (MOD(NK-1,NCPU)+1 ==  W%WDES%NB_LOW) THEN
          ! receive from all other nodes and local copy
          DO NN=1,NCPU
             N_INTO_TOT=(N-1)*NCPU+NN
             IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
                IF (NN==W%WDES%NB_LOW) THEN
                   CALL W1_COPY( WTMP, W1(N_INTO_TOT-NB1+1) )
                ELSE
!                   WRITE(*,*) W%WDES%NB_LOW, 'receive', N_INTO_TOT-NB1+1, 'from' , NN
                   CALL M_recv_z(W%WDES%COMM_INTER, NN , & 
                        W1(N_INTO_TOT-NB1+1)%CR(1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)
#ifndef gammareal
                   IF (W%WDES%LOVERL) CALL M_recv_z(W%WDES%COMM_INTER, NN, & 
                        W1(N_INTO_TOT-NB1+1)%CPROJ(1), W%WDES%NPROD)
#else
                   IF (W%WDES%LOVERL) CALL M_recv_d(W%WDES%COMM_INTER, NN, & 
                        W1(N_INTO_TOT-NB1+1)%CPROJ(1), W%WDES%NPROD)
#endif
                ENDIF
             ENDIF
          END DO
       ELSE
          IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
!             WRITE(*,*) W%WDES%NB_LOW,'send', N_INTO_TOT, 'to' , MOD(NK-1,NCPU)+1
             CALL M_send_z(W%WDES%COMM_INTER, MOD(NK-1,NCPU)+1, & 
                  WTMP%CR(1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)
#ifndef gammareal
             IF (W%WDES%LOVERL) CALL M_send_z(W%WDES%COMM_INTER, MOD(NK-1,NCPU)+1, & 
                  WTMP%CPROJ(1), W%WDES%NPROD)
#else
             IF (W%WDES%LOVERL) CALL M_send_d(W%WDES%COMM_INTER, MOD(NK-1,NCPU)+1, & 
                  WTMP%CPROJ(1), W%WDES%NPROD)
#endif
          ENDIF
       ENDIF
#endif
    ENDDO

    CALL DELWAV(WTMP, .TRUE.)

  END SUBROUTINE W1_GATHER_KSEL


!************************ SUBROUTINE W1_GATHER_KNODESEL ***************
!
! This subroutine is similar to W1_GATHER_KSEL, but at variance with
! it, one explicitly select the node on which you gather the wave
! fucntion (instead of distributing k-points in round robin fashion).
!
!**********************************************************************

  SUBROUTINE W1_GATHER_KNODESEL( W, NB1, NB2, ISP, W1, NODE)
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed
    INTEGER :: NODE          ! which node receives

   ! local
    TYPE (wavefun1):: WTMP   ! temporary for FFT
    INTEGER :: NN, N, NLOC, NCPU

    NCPU=W%WDES%NB_PAR

    CALL NEWWAV(WTMP, W1(1)%WDES1, .TRUE.)

    DO N=(NB1-1)/W%WDES%NB_PAR+1,(NB2-1)/W%WDES%NB_PAR+1
       N_INTO_TOT=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
          CALL W1_COPY( ELEMENT( W, W1(1)%WDES1, N, ISP), WTMP )
          CALL FFTWAV_W1( WTMP)
       ENDIF
#ifdef MPI
       IF (NODE==W%WDES%NB_LOW) THEN
          ! receive from all other nodes and local copy
          DO NN=1,NCPU
             N_INTO_TOT=(N-1)*NCPU+NN
             IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
                IF (NN==W%WDES%NB_LOW) THEN
                   CALL W1_COPY( WTMP, W1(N_INTO_TOT-NB1+1) )
                ELSE
                   CALL M_recv_z(W%WDES%COMM_INTER, NN , & 
                        W1(N_INTO_TOT-NB1+1)%CR(1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)
#ifndef gammareal
                   IF (W%WDES%LOVERL) CALL M_recv_z(W%WDES%COMM_INTER, NN, & 
                        W1(N_INTO_TOT-NB1+1)%CPROJ(1), W%WDES%NPROD)
#else
                   IF (W%WDES%LOVERL) CALL M_recv_d(W%WDES%COMM_INTER, NN, & 
                        W1(N_INTO_TOT-NB1+1)%CPROJ(1), W%WDES%NPROD)
#endif
                ENDIF
             ENDIF
          END DO
       ELSE
          IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
             CALL M_send_z(W%WDES%COMM_INTER,NODE, & 
                  WTMP%CR(1), W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)
#ifndef gammareal
             IF (W%WDES%LOVERL) CALL M_send_z(W%WDES%COMM_INTER,NODE, & 
                  WTMP%CPROJ(1), W%WDES%NPROD)
#else
             IF (W%WDES%LOVERL) CALL M_send_d(W%WDES%COMM_INTER,NODE, & 
                  WTMP%CPROJ(1), W%WDES%NPROD)
#endif
          ENDIF
       ENDIF
#endif
    ENDDO
    CALL DELWAV(WTMP, .TRUE.)
  END SUBROUTINE W1_GATHER_KNODESEL


!************************ SUBROUTINE W1_GATHER_W1 *********************
!
! This subroutine gathers from a wavefunction array
! instead of W
! 
!**********************************************************************

  SUBROUTINE W1_GATHER_W1( W, NB2, W1_ORIG, W1)
    TYPE (wavespin):: W       ! wavefunction
    TYPE (wavefun1):: W1_ORIG(:) ! wavefunction array to be merged
    TYPE (wavefun1):: W1(:)   ! array into which the merge is performed
    INTEGER :: NB2            ! final band

  ! local
    INTEGER :: NN, N, NLOC, NCPU

    NCPU=W%WDES%NB_PAR
    IF (SIZE(W1_ORIG)< NB2) THEN
       CALL vtutor%bug("internal error in W1_GATHER_W1: W1_ORIG is not sufficiently large", __FILE__, __LINE__)
    ENDIF
       
    DO N=1,NB2
       NN=(N-1)*NCPU+W%WDES%NB_LOW
       CALL W1_COPY_NOCR( W1_ORIG(N), W1(NN) )
       CALL FFTWAV_W1( W1(NN))
    ENDDO

    NLOC=NB2*W%WDES%NB_PAR

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO NN=1,NLOC
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(NN-1,NCPU)+1)
#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(NN)%CPROJ(1), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#endif
       ENDDO

#ifdef MPI_barrier_after_bcast
      CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif

    ENDIF
#endif

  END SUBROUTINE W1_GATHER_W1


!************************ SUBROUTINE W1_GATHER_STRIP ******************
!
! This subroutine gathers from a wavefunction array instead of W.
! Collected are the bands W1_ORIG(N): NB1<=NB_LOW+(N-1)*NB_PAR<=NB2,
! into W1(1:NB2-NB1+1)
! 
!**********************************************************************

  SUBROUTINE W1_GATHER_STRIP( W, NB1, NB2, W1_ORIG, W1)
    TYPE (wavespin):: W       ! wavefunction
    TYPE (wavefun1):: W1_ORIG(:) ! wavefunction array to be merged
    TYPE (wavefun1):: W1(:)   ! array into which the merge is performed
    INTEGER :: NB1, NB2       ! final band

  ! local
    INTEGER :: NN, N, NCPU

    PROFILING_START('w1_gather_strip')

    NCPU=W%WDES%NB_PAR
    IF (SIZE(W1_ORIG)*NCPU<NB2) THEN
       CALL vtutor%bug("internal error in W1_GATHER_STRIP: W1_ORIG is not sufficiently large", __FILE__, __LINE__)
    ENDIF

    IF (SIZE(W1)<NB2-NB1+1) THEN
       CALL vtutor%bug("internal error in W1_GATHER_STRIP: W1 is not sufficiently large", __FILE__, __LINE__)
    ENDIF
       
    DO N=1,SIZE(W1_ORIG)
       NN=(N-1)*NCPU+W%WDES%NB_LOW
       IF (NN<NB1.OR.NN>NB2) CYCLE
       CALL W1_COPY_NOCR( W1_ORIG(N), W1(NN-NB1+1) )
       CALL FFTWAV_W1( W1(NN-NB1+1))
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%DO_REDIS) THEN
       DO NN=NB1,NB2
          CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN-NB1+1)%CR(1), &
               W%WDES%GRID%MPLWV*W%WDES%NRSPINORS,MOD(NN-1,NCPU)+1)
#ifndef gammareal
          IF (W%WDES%LOVERL) CALL M_bcast_z_from(W%WDES%COMM_INTER,W1(NN-NB1+1)%CPROJ(1), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#else
          IF (W%WDES%LOVERL) CALL M_bcast_d_from(W%WDES%COMM_INTER,W1(NN-NB1+1)%CPROJ(1), &
               W%WDES%NPROD,MOD(NN-1,NCPU)+1)
#endif
       ENDDO

#ifdef MPI_barrier_after_bcast
      CALLMPI( M_barrier( W%WDES%COMM_INTER ))
#endif

    ENDIF
#endif
    PROFILING_STOP('w1_gather_strip')

  END SUBROUTINE W1_GATHER_STRIP


!************************ SUBROUTINE W1_IGATHER_STRIP *****************
!
! Gathers a set of wavefunctions starting from band NB1 until NB2
! (global band indices) to all nodes using non-blocking bcast_from.
! This subroutine gathers from a wavefunction array W1_ORIG instead
! of W.
!
! In case the code is compiled with -Dshmem_bcast_buffer communication
! will be between ranks within COMM_inter_node, i.e., only those ranks
! within COMM_INTER that are NOT on the same physical node will talk to
! eachother. MPI communication between ranks within COMM_INTER that
! reside on the same node is not necessary since they access a common
! shared memory segment.
!
! In case the code is NOT compiled with -Dshmem_bcast_buffer, the
! communication will be between all ranks in COMM_INTER.
! 
!**********************************************************************
#ifdef shmem_bcast_buffer
  SUBROUTINE W1_IGATHER_STRIP( W, NB1, NB2, W1_ORIG, W1)
    TYPE (wavespin):: W          ! wavefunction
    TYPE (wavefun1):: W1_ORIG(:) ! wavefunction array to be merged
    TYPE (wavefun1):: W1(:)      ! array into which the merge is performed
    INTEGER :: NB1,NB2           ! first and final band

  ! local
    INTEGER :: NI,N,NB_LOCAL
    INTEGER :: ierror

    INTEGER :: IDO(NB2-NB1+1)

    INTEGER :: nrequests,requests(3*(NB2-NB1+1))


    PROFILING_START('w1_igather_strip')

    IF (SIZE(W1_ORIG)*W%WDES%NB_PAR<NB2) THEN
       CALL vtutor%bug("internal error in W1_IGATHER_STRIP: W1_ORIG is not sufficiently large", __FILE__, __LINE__)
    ENDIF

    IF (SIZE(W1)<NB2-NB1+1) THEN
       CALL vtutor%bug("internal error in W1_IGATHER_STRIP: W1 is not sufficiently large", __FILE__, __LINE__)
    ENDIF
       
    CALL M_barrier(W%WDES%COMM_intra_node)

    DO N=NB1,NB2
       NI=N-NB1+1
       IF (MOD(N-1,W%WDES%NB_PAR)+1==W%WDES%NB_LOW) THEN
          NB_LOCAL=1+(N-1)/W%WDES%NB_PAR
          CALL W1_COPY_NOCR(W1_ORIG(NB_LOCAL),W1(NI)); CALL FFTWAV_W1(W1(NI))
          IDO(NI)=W%WDES%COMM_inter_node%NODE_ME
       ELSE
          IDO(NI)=0
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%COMM_inter_node%NCPU>1) THEN

       CALLMPI( M_sum_i(W%WDES%COMM_inter_node,IDO,NB2-NB1+1))

       nrequests=0
       DO N=NB1,NB2
          NI=N-NB1+1
          IF (IDO(NI)>0) THEN
             nrequests=nrequests+1
             CALL M_ibcast_z_from(W%WDES%COMM_inter_node,W1(NI)%CW(1), &
            &     SIZE(W1(NI)%CW),IDO(NI),requests(nrequests))
             nrequests=nrequests+1
             CALL M_ibcast_z_from(W%WDES%COMM_inter_node,W1(NI)%CR(1), &
            &     SIZE(W1(NI)%CR),IDO(NI),requests(nrequests))
             IF (W%WDES%LOVERL) THEN
                nrequests=nrequests+1
#ifndef gammareal
                CALL M_ibcast_z_from(W%WDES%COMM_inter_node,W1(NI)%CPROJ(1), &
               &     SIZE(W1(NI)%CPROJ),IDO(NI),requests(nrequests))
#else
                CALL M_ibcast_d_from(W%WDES%COMM_inter_node,W1(NI)%CPROJ(1), &
               &     SIZE(W1(NI)%CPROJ),IDO(NI),requests(nrequests))
#endif
             ENDIF
          ENDIF
       ENDDO

       CALL M_waitall(nrequests,requests(1))

    ENDIF

    CALL M_barrier(W%WDES%COMM_intra_node)
!!#ifndef shmem_bcast_buffer
!!    IF (W%WDES%COMM_intra_node%NCPU>1) THEN
!!       DO N=1,NB2-NB1+1
!!          IF (IDO(N)>0) IDO(N)=W%WDES%COMM_intra_node%NODE_ME
!!       ENDDO
!!
!!       CALLMPI( M_sum_i(W%WDES%COMM_intra_node,IDO,NB2-NB1+1))
!!
!!       nrequests=0
!!       DO N=NB1,NB2
!!          NI=N-NB1+1
!!          IF (IDO(NI)>0) THEN
!!             nrequests=nrequests+1
!!             CALL M_ibcast_z_from(W%WDES%COMM_intra_node,W1(NI)%CW(1), &
!!            &     SIZE(W1(NI)%CW),IDO(NI),requests(nrequests))
!!             nrequests=nrequests+1
!!             CALL M_ibcast_z_from(W%WDES%COMM_intra_node,W1(NI)%CR(1), &
!!            &     SIZE(W1(NI)%CR),IDO(NI),requests(nrequests))
!!             IF (W%WDES%LOVERL) THEN
!!                nrequests=nrequests+1
!!#ifndef gammareal
!!                CALL M_ibcast_z_from(W%WDES%COMM_intra_node,W1(NI)%CPROJ(1), &
!!               &     SIZE(W1(NI)%CPROJ),IDO(NI),requests(nrequests))
!!#else
!!                CALL M_ibcast_d_from(W%WDES%COMM_intra_node,W1(NI)%CPROJ(1), &
!!               &     SIZE(W1(NI)%CPROJ),IDO(NI),requests(nrequests))
!!#endif
!!             ENDIF
!!          ENDIF
!!       ENDDO
!!
!!       CALL M_waitall(nrequests,requests(1))
!!
!!    ENDIF
!!#endif
#endif
    PROFILING_STOP('w1_igather_strip')

  END SUBROUTINE W1_IGATHER_STRIP
#else
  SUBROUTINE W1_IGATHER_STRIP( W, NB1, NB2, W1_ORIG, W1)
    TYPE (wavespin):: W          ! wavefunction
    TYPE (wavefun1):: W1_ORIG(:) ! wavefunction array to be merged
    TYPE (wavefun1):: W1(:)      ! array into which the merge is performed
    INTEGER :: NB1,NB2           ! first and final band

  ! local
    INTEGER :: NI,N,NB_LOCAL
    INTEGER :: ierror

    INTEGER :: nrequests,requests(3*(NB2-NB1+1))


    PROFILING_START('w1_igather_strip')

    IF (SIZE(W1_ORIG)*W%WDES%NB_PAR<NB2) THEN
       CALL vtutor%bug("internal error in W1_IGATHER_STRIP: W1_ORIG is not sufficiently large", __FILE__, __LINE__)
    ENDIF

    IF (SIZE(W1)<NB2-NB1+1) THEN
       CALL vtutor%bug("internal error in W1_IGATHER_STRIP: W1 is not sufficiently large", __FILE__, __LINE__)
    ENDIF
       
    DO N=NB1,NB2
       NI=N-NB1+1
       IF (MOD(N-1,W%WDES%NB_PAR)+1==W%WDES%NB_LOW) THEN
          NB_LOCAL=1+(N-1)/W%WDES%NB_PAR
          CALL W1_COPY_NOCR(W1_ORIG(NB_LOCAL),W1(NI))
          CALL FFTWAV_W1(W1(NI))
       ENDIF
    ENDDO

    ! distribute W1 to all nodes
#ifdef MPI
    IF (W%WDES%COMM_INTER%NCPU>1) THEN

       nrequests=0
       DO N=NB1,NB2
          NI=N-NB1+1

          nrequests=nrequests+1
          CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CW(1), &
         &     SIZE(W1(NI)%CW),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
          nrequests=nrequests+1
          CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CR(1), &
         &     SIZE(W1(NI)%CR),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
          IF (W%WDES%LOVERL) THEN
             nrequests=nrequests+1
#ifndef gammareal
             CALL M_ibcast_z_from(W%WDES%COMM_INTER,W1(NI)%CPROJ(1), &
            &     SIZE(W1(NI)%CPROJ),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
#else
             CALL M_ibcast_d_from(W%WDES%COMM_INTER,W1(NI)%CPROJ(1), &
            &     SIZE(W1(NI)%CPROJ),MOD(N-1,W%WDES%NB_PAR)+1,requests(nrequests))
#endif
          ENDIF
       ENDDO

       CALL M_waitall(nrequests,requests(1))

    ENDIF
#endif
    PROFILING_STOP('w1_igather_strip')

  END SUBROUTINE W1_IGATHER_STRIP
#endif

!************************* SUBROUTINE WVREAL_PRECISE *******************
!
! this subroutine forces the wavefunction to be real at the Gamma-point
! it is required for the gamma point only mode
! to avoid that small non real components develop
! this version is exact and works through an FFT to real space
! and the forcing the wavefunction to be real
! the routine is required if only subspace rotations are performed
! since the subspace rotation routine can not force wavefunctions
! to become real
!
! if LORBITALREAL is .TRUE., the routine also tries to make the
! orbitals real (in real space) by calculating the orbital
!   phi(r) = u_k(r) e^ikr
! taking the real part, storing the coefficients back
! this is possible at special k-points such as points at the BZ boundary
! e.g. k=(+-0.5,+-0.5,+-0.5)
!      k=(+-0.5,+-0.5,  0)
!      k=(+-0.5,  0  ,  0)    and permutations thereof
! the routine will not work properly at any other k-points
! (and is in fact by-passed in this case)
! ideally no subspace rotation should be called after calling this
! routine (though orthogonalization is fine)
!
!***********************************************************************

  SUBROUTINE WVREAL_PRECISE(W)
#ifdef _OPENACC
    USE mopenacc
#endif
    USE constant
    USE random_seeded, ONLY: RANE
    IMPLICIT NONE
    TYPE (wavespin) W
  ! local
    INTEGER :: NK, ISP, NB, ISPINOR, N
    TYPE (wavedes1)    WDES1          ! descriptor for one k-point
    TYPE (wavefun1)    W1             ! current wavefunction
    LOGICAL :: LPHASE
    COMPLEX(q) :: CPHASE(W%WDES%GRID%MPLWV)
    COMPLEX(q) :: AVERAGE_PHASE

! for spinors it is not (always) possible to make wave functions real-valued
! (for that reason we need to use the "_ncl" of VASP even for Gamma-only ...)
     IF (W%WDES%LNONCOLLINEAR) RETURN

#ifdef gammareal
!$ACC ENTER DATA CREATE(W1,WDES1) IF(ACC_EXEC_ON)
    CALL SETWDES(W%WDES,WDES1,0)
    CALL NEWWAV(W1, WDES1, .TRUE.)
    !
    ! search for position of G=0 zero element
    !
    DO NK  =1,W%WDES%NKPTS
#ifdef MPI
       IF (MOD(NK-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
       IF (ABS(SUM(W%WDES%VKPT(:,NK)*W%WDES%VKPT(:,NK)))< 1E-12) THEN
          CALL SETWDES(W%WDES,WDES1,NK)
          DO ISP=1,W%WDES%ISPIN
             DO NB=1,W%WDES%NBANDS
                CALL W1_COPY(ELEMENT(W, WDES1, NB, ISP), W1)
                CALL FFTWAV_W1(W1)
                DO ISPINOR =0,WDES1%NRSPINORS-1
!$ACC PARALLEL LOOP PRESENT(W1,W1%CR,WDES1,WDES1%GRID,WDES1%GRID%MPLWV) __IF_ASYNC__
                   DO N=1,WDES1%GRID%RL%NP
                      W1%CR(N+ISPINOR*WDES1%GRID%MPLWV)=REAL(W1%CR(N+ISPINOR*WDES1%GRID%MPLWV),q) & 
                          *(1.0_q/WDES1%GRID%NPLWV)
                   ENDDO
                   CALL FFTEXT(WDES1%NGVECTOR, WDES1%NINDPW(1), W1%CR(1+ISPINOR*WDES1%GRID%MPLWV),W%CW(1+ISPINOR*WDES1%NGVECTOR,NB,NK,ISP),WDES1%GRID,.FALSE.)
                ENDDO
             ENDDO
          ENDDO
       END IF
    ENDDO

    CALL DELWAV(W1, .TRUE.)
#ifdef _OPENACC
!$ACC EXIT DATA DELETE(W1) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)
    CALL ACC_DELETE_TYPED_VAR(WDES1)
#endif
#else
    IF (W%WDES%LORBITALREAL) THEN
!$ACC ENTER DATA CREATE(W1,WDES1) IF(ACC_EXEC_ON)
    CALL SETWDES(W%WDES,WDES1,0)
    CALL NEWWAV(W1, WDES1, .TRUE.)
!
! force orbitals to be real for non-gamma point only method
!
    DO NK  =1,W%WDES%NKPTS
#ifdef MPI
       IF (MOD(NK-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
       IF (ABS(MOD(W%WDES%VKPT(1,NK)*2+100,1._q))< 1E-6 .AND. &
           ABS(MOD(W%WDES%VKPT(2,NK)*2+100,1._q))< 1E-6 .AND. &
           ABS(MOD(W%WDES%VKPT(3,NK)*2+100,1._q))< 1E-6) THEN
       CALL SETPHASE_WVREAL(W%WDES%VKPT(:,NK), W%WDES%GRID, CPHASE, LPHASE)
!$ACC ENTER DATA COPYIN(CPHASE) IF(ACC_EXEC_ON .AND. LPHASE) ASYNC(ACC_ASYNC_Q)
       
       CALL SETWDES(W%WDES,WDES1,NK)
       DO ISP=1,W%WDES%ISPIN
          DO NB=1,W%WDES%NBANDS
             CALL W1_COPY(ELEMENT(W, WDES1, NB, ISP), W1)
             CALL FFTWAV_W1(W1)
             ! chose random phase in 1st and 4th quadrant, avoid proximity to +-i
             AVERAGE_PHASE=EXP(CITPI*(RANE()-0.5)/2.1)
             CALLMPI (M_bcast_z( WDES1%COMM_INB, AVERAGE_PHASE, 1))

             DO ISPINOR =0,WDES1%NRSPINORS-1
                IF (LPHASE) THEN
!$ACC PARALLEL LOOP PRESENT(W1,W1%CR,WDES1,WDES1%GRID,WDES1%GRID%MPLWV,CPHASE) &
!$ACC& __IF_ASYNC__
                   DO N=1,WDES1%GRID%RL%NP
                      W1%CR(N+ISPINOR*WDES1%GRID%MPLWV)=CONJG(CPHASE(N))*REAL(W1%CR(N+ISPINOR*WDES1%GRID%MPLWV)*CPHASE(N)*AVERAGE_PHASE,q) & 
                      *(1.0_q/WDES1%GRID%NPLWV)
                   ENDDO
                ELSE
!$ACC PARALLEL LOOP PRESENT(W1,W1%CR,WDES1,WDES1%GRID,WDES1%GRID%MPLWV) __IF_ASYNC__
                   DO N=1,WDES1%GRID%RL%NP
                      W1%CR(N+ISPINOR*WDES1%GRID%MPLWV)=REAL(W1%CR(N+ISPINOR*WDES1%GRID%MPLWV)*AVERAGE_PHASE,q) & 
                      *(1.0_q/WDES1%GRID%NPLWV)
                   ENDDO
                ENDIF

                CALL FFTEXT(WDES1%NGVECTOR, WDES1%NINDPW(1), W1%CR(1+ISPINOR*WDES1%GRID%MPLWV),W%CW(1+ISPINOR*WDES1%NGVECTOR,NB,NK,ISP),WDES1%GRID,.FALSE.)
             ENDDO
          ENDDO
       ENDDO
!$ACC EXIT DATA DELETE(CPHASE) IF(ACC_EXEC_ON .AND. LPHASE) ASYNC(ACC_ASYNC_Q)
       ENDIF
    ENDDO

    CALL DELWAV(W1, .TRUE.)
#ifdef _OPENACC
!$ACC EXIT DATA DELETE(W1) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)
    CALL ACC_DELETE_TYPED_VAR(WDES1)
#endif
    ENDIF
#endif

    RETURN
  END SUBROUTINE WVREAL_PRECISE


  SUBROUTINE SETPHASE_WVREAL(VKPT, GRID, CPHASE, LPHASE)
    USE constant

    IMPLICIT NONE

    REAL(q) :: VKPT(3)
    TYPE (grid_3d) GRID
    COMPLEX(q) :: CPHASE(GRID%MPLWV)
    LOGICAL LPHASE
  ! local
    REAL(q),PARAMETER :: TINY=1E-6_q
    REAL(q) F1, F2, F3
    INTEGER NC, N, IND
    COMPLEX(q) C, CD, CSUM
    CSUM=0
#ifdef gammareal
    LPHASE=.FALSE.
    RETURN
#endif
    
    IF (ABS(VKPT(1))>TINY .OR. ABS(VKPT(2))>TINY .OR. ABS(VKPT(3))>TINY) THEN
       LPHASE=.TRUE.
       F1=TPI/GRID%NGX*VKPT(1)
       F2=TPI/GRID%NGY*VKPT(2)
       F3=TPI/GRID%NGZ*VKPT(3)

       IF (GRID%RL%NFAST==3) THEN
          CD=EXP(CMPLX(0,F3,q))
          IND=0
          DO NC=1,GRID%RL%NCOL
             C=EXP(CMPLX(0,F1*(GRID%RL%I2(NC)-1)+F2*(GRID%RL%I3(NC)-1),q))
             DO N=1,GRID%RL%NROW
                IND=IND+1
                CPHASE(IND)=C
                C=C*CD
             ENDDO
          ENDDO
       ELSE
          CD=EXP(CMPLX(0,F1,q))
          IND=0
          DO NC=1,GRID%RL%NCOL
             C=EXP(CMPLX(0,F2*(GRID%RL%I2(NC)-1)+F3*(GRID%RL%I3(NC)-1),q))
             DO N=1,GRID%RL%NROW
                IND=IND+1
                CPHASE(IND)=C
                CSUM=CSUM+C
                C=C*CD
             ENDDO
          ENDDO
       ENDIF
       LPHASE=.TRUE.
    ELSE
       LPHASE=.FALSE.
    ENDIF
  END SUBROUTINE SETPHASE_WVREAL

END MODULE wave_high


!***********************************************************************
!
! assign a a W1 structure an initial value from a wavefunction 
! array
! the plane wave coefficients and the wave function character (optional)
! must be supplied
! for performance reason no runtime checking on anything is performed
!      
!***********************************************************************

  SUBROUTINE ARRAY_TO_W1( W1, C, CPROJ)
    USE prec
    USE wave
    IMPLICIT NONE
    TYPE (wavefun1)    W1
    COMPLEX(q):: C(*)
    GDEF, OPTIONAL :: CPROJ(*)

#if defined (use_fastbcopy)
    CALL FASTBCOPY8(W1%WDES1%NRPLWV*2 ,C, W1%CW(1))
#else
    CALL ZCOPY( W1%WDES1%NRPLWV, C, 1, W1%CW(1), 1)
#endif
    IF (PRESENT(CPROJ)) THEN
       IF (W1%WDES1%LGAMMA) THEN
          CALL DCOPY( W1%WDES1%NPROD, CPROJ(1), 1, W1%CPROJ(1), 1)
       ELSE
          CALL ZCOPY( W1%WDES1%NPROD, CPROJ(1), 1, W1%CPROJ(1), 1)
       ENDIF
    END IF


  END SUBROUTINE ARRAY_TO_W1


!***********************************************************************
!
! assign a W1 structure to an array (reverse of the previous operation)
!      
!***********************************************************************

  SUBROUTINE W1_TO_ARRAY( W1, C,  CPROJ)
    USE prec
    USE wave
    IMPLICIT NONE
    TYPE (wavefun1)    W1
    COMPLEX(q) :: C(*)
    GDEF, OPTIONAL :: CPROJ(*)

#if defined (use_fastbcopy)
    CALL FASTBCOPY8(W1%WDES1%NRPLWV*2 , W1%CW(1), C)
#else
    CALL ZCOPY( W1%WDES1%NRPLWV, W1%CW(1), 1, C, 1)
#endif
    IF (PRESENT(CPROJ)) THEN
       IF (W1%WDES1%LGAMMA) THEN
          CALL DCOPY( W1%WDES1%NPROD,  W1%CPROJ(1), 1, CPROJ(1), 1)
       ELSE
          CALL ZCOPY( W1%WDES1%NPROD,  W1%CPROJ(1), 1, CPROJ(1), 1)
       ENDIF
    END IF

  END SUBROUTINE W1_TO_ARRAY


!************************* SUBROUTINE ECCP_NL   ************************
!
! this subroutine calculates the expectation value of <c|H|cp>
! where c and cp are two wavefunctions; non local part only
! for one ion only
! I have put this in a separate routine because optimization
! is than easier
!
!***********************************************************************

  SUBROUTINE ECCP_NL(LMDIM,LMMAXC,CDIJ,CPROJ1,CPROJ2,CNL)
!$ACC ROUTINE VECTOR
    USE prec
    IMPLICIT NONE
    GDEF      CNL
    INTEGER LMDIM, LMMAXC
    OVERLAP CDIJ(LMDIM,LMDIM)
    GDEF CPROJ1(LMMAXC),CPROJ2(LMMAXC)
 ! local
    INTEGER L, LP

!   PROFILING_START('eccp_nl')

!DIR$ IVDEP
!OCL NOVREL
!$ACC LOOP VECTOR COLLAPSE(2) REDUCTION(+:CNL)
    DO L=1,LMMAXC
       DO LP=1,LMMAXC
          CNL=CNL+CDIJ(LP,L)*CPROJ1(LP)*GCONJG(CPROJ2(L))
       ENDDO
    ENDDO

!   PROFILING_STOP('eccp_nl')

  END SUBROUTINE ECCP_NL


!************************* SUBROUTINE OVERL ***************************
!
! F77 low level routine
! calculate the result of the overlap-operator acting onto a set of
! wave function characters
!
!' CRESUL^n_N,nlm = sum_n'l'm' D_N,n'l'm',nlm CPROF^n_N,n'l'm'  
!  CRESUL^n_N,i   = sum_j D_N,j,i CPROF^n_N,j  
!
! n is the band index, N the ion index, nlm=i the index for the 
! one centre partial waves
!
!> @details @ref openmp :
!> the loop over the bands owned locally by a particular MPI-rank
!> is distributed over all available OpenMP threads.
!
!**********************************************************************

  SUBROUTINE OVERL(WDES1, LOVERL, LMDIM, CQIJ, CPROF, CRESUL)
!$ACC ROUTINE(SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX) SEQ
#ifdef _OPENACC
    USE mopenacc
#endif
    USE wave
    IMPLICIT NONE
    TYPE (wavedes1) WDES1
    LOGICAL LOVERL
    INTEGER LMDIM
    OVERLAP CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    GDEF CRESUL(WDES1%NPROD,WDES1%NBANDS),CPROF(WDES1%NPROD,WDES1%NBANDS)
  ! local
    GDEF CTMP
    INTEGER NB,NP,ISPINOR,ISPINOR_,NPRO,NPRO_,NT,NI,LMMAXC,L,LP

    PROFILING_START('overl')

    IF (LOVERL) THEN
#ifdef _OPENACC
!$ACC KERNELS PRESENT(CRESUL)  __IF_ASYNC__
       CRESUL(1:WDES1%NPRO,1:WDES1%NBANDS)=0
!$ACC END KERNELS
#endif
NOACC !$OMP PARALLEL DO DEFAULT(NONE) &
NOACC !$OMP SHARED(WDES1,CRESUL,CQIJ,CPROF) &
NOACC !$OMP PRIVATE(NB,NP,ISPINOR,ISPINOR_,NPRO,NPRO_,NI,NT,LMMAXC,L,LP,CTMP)
!$ACC PARALLEL LOOP PRESENT(WDES1,CQIJ,CPROF,CRESUL) COLLAPSE(4) GANG __IF_ASYNC__ &
!$ACC PRIVATE(NP,ISPINOR,ISPINOR_,NI,NT,LMMAXC,NPRO,NPRO_,L,LP,CTMP)
       bands: DO NB=1,WDES1%NBANDS
NOACC     DO NP=1,WDES1%NPRO
NOACC        CRESUL(NP,NB)=0
NOACC     ENDDO

          spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
          DO ISPINOR_=0,WDES1%NRSPINORS-1
             DO NI=1,WDES1%NIONS
                NT=WDES1%ITYP(NI)
                LMMAXC=WDES1%LMMAX(NT)
                IF (LMMAXC==0) CYCLE
                NPRO =WDES1%LMBASE(NI)+ISPINOR *(WDES1%NPRO/2)
                NPRO_=WDES1%LMBASE(NI)+ISPINOR_*(WDES1%NPRO/2)

!$ACC LOOP
!NEC$ select_vector
                DO L =1,LMMAXC
#ifndef __NEC__
                   CTMP=0
!DIR$ IVDEP
!OCL NOVREC
!$ACC LOOP VECTOR REDUCTION(+:CTMP)
                   DO LP=1,LMMAXC
                      CTMP=CTMP+CQIJ(LP,L,NI,1+ISPINOR_+2*ISPINOR)*CPROF(LP+NPRO_,NB)
                   ENDDO
#if !defined(_OPENACC) || defined(gammareal)
!$ACC ATOMIC UPDATE
                   CRESUL(L+NPRO,NB)=CRESUL(L+NPRO,NB)+CTMP
#else
                   CALL SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX(CRESUL(L+NPRO,NB),CTMP)
#endif
#else // __NEC__
!NEC$ outerloop_unroll(8)
                   DO LP=1,LMMAXC
                      CRESUL(L+NPRO,NB)=CRESUL(L+NPRO,NB)+CQIJ(LP,L,NI,1+ISPINOR_+2*ISPINOR)*CPROF(LP+NPRO_,NB)
                   ENDDO
#endif // __NEC__
                ENDDO
             ENDDO
          ENDDO
          ENDDO spinor

       ENDDO bands
NOACC !$OMP END PARALLEL DO
    ENDIF

    PROFILING_STOP('overl')

  END SUBROUTINE OVERL


!************************* SUBROUTINE OVERL1 **************************
!
! F77 low level routine
! calculate the result of the overlap-operator acting onto one
! wave function character
!
!  CRESUL_N,nlm = sum_n'l'm' (D_N,n'l'm',nlm-e Q_N,n'l'm',nlm) CPROF_N,n'l'm'  
!  CRESUL_N,i   = sum_j (D_N,j,i -e Q_N,j,i) CPROF_N,j  
!
! N the ion index, nlm=i the index for the one centre partial waves
!
!> @details @ref openmp :
!> the nested loops over \"types\" + \"ions-of-type\" are replaced
!> by a single loop over \"all ions\" that is distributed over
!> all available threads.
!
!**********************************************************************

!
! scheduled for removal to be replaced by OVERL1_
!
  SUBROUTINE OVERL1(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF, CRESUL)
#ifdef _OPENACC
    USE mopenacc_struct_def
#endif
    USE wave
    IMPLICIT NONE

    TYPE (wavedes1) WDES1
    INTEGER LMDIM
    OVERLAP CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
         CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    REAL(q) :: EVALUE
    GDEF CRESUL(WDES1%NPRO),CPROF(WDES1%NPRO)
  ! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NPRO2, NPRO2_, NT, NI, LMMAXC, L, LP
    GDEF CTMP

    PROFILING_START('overl1')

!$ACC KERNELS PRESENT(CRESUL) __IF_ASYNC__
    CRESUL=0
!$ACC END KERNELS

!$ACC ENTER DATA COPYIN(EVALUE) __IF_ASYNC__

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
    DO ISPINOR_=0,WDES1%NRSPINORS-1

       NPRO =ISPINOR *(WDES1%NPRO/2)
       NPRO_=ISPINOR_*(WDES1%NPRO/2)
#ifdef _OPENACC
!$ACC PARALLEL LOOP GANG VECTOR_LENGTH(32) __IF_ASYNC__ &
!$ACC& PRESENT(WDES1,CDIJ,CQIJ,CPROF,CRESUL,EVALUE) PRIVATE(NT,LMMAXC,NPRO2,NPRO2_,L,LP,CTMP)
#else
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP SHARED(NPRO,NPRO_,WDES1,EVALUE,ISPINOR,ISPINOR_,CDIJ,CQIJ,CPROF,CRESUL) &
!$OMP PRIVATE(NI,NT,LMMAXC,NPRO2,NPRO2_,L,LP,CTMP)
#endif
!NEC$ novector
       DO NI=1,WDES1%NIONS
          NT=WDES1%ITYP(NI)
          LMMAXC=WDES1%LMMAX(NT)
          IF (LMMAXC==0) CYCLE
          NPRO2 =WDES1%LMBASE(NI)+NPRO
          NPRO2_=WDES1%LMBASE(NI)+NPRO_
          IF (EVALUE==0) THEN
!$ACC LOOP VECTOR PRIVATE(CTMP)
!NEC$ select_vector
             DO L =1,LMMAXC
#ifndef __NEC__
                CTMP=0
!DIR$ IVDEP
!OCL NOVREC
                DO LP=1,LMMAXC
                   CTMP=CTMP + &
                        CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)*CPROF(LP+NPRO2_)
                ENDDO
                CRESUL(L+NPRO2)=CRESUL(L+NPRO2)+CTMP
#else // __NEC__
!NEC$ outerloop_unroll(8)
                DO LP=1,LMMAXC
                   CRESUL(L+NPRO2)=CRESUL(L+NPRO2) + &
                        CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)*CPROF(LP+NPRO2_)
                ENDDO
#endif // __NEC__
             ENDDO
          ELSE
!$ACC LOOP VECTOR PRIVATE(CTMP)
!NEC$ select_vector
             DO L =1,LMMAXC
#ifndef __NEC__
                CTMP=0
!DIR$ IVDEP
!OCL NOVREC
                DO LP=1,LMMAXC
                   CTMP=CTMP+ &
                       (CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)- &
                        EVALUE*CQIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)) * CPROF(LP+NPRO2_)
                ENDDO
                CRESUL(L+NPRO2)=CRESUL(L+NPRO2)+CTMP
#else // __NEC__
!NEC$ outerloop_unroll(8)
                DO LP=1,LMMAXC
                   CRESUL(L+NPRO2)=CRESUL(L+NPRO2)+ &
                       (CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)- &
                        EVALUE*CQIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)) * CPROF(LP+NPRO2_)
                ENDDO
#endif // __NEC__
             ENDDO
          ENDIF
       ENDDO
#ifndef _OPENACC
!$OMP END PARALLEL DO
#endif
    ENDDO
    ENDDO spinor

!$ACC EXIT DATA DELETE(EVALUE) __IF_ASYNC__

    PROFILING_STOP('overl1')

  END SUBROUTINE OVERL1

#ifdef _OPENACC
  SUBROUTINE OVERL_MU(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, W1, CRESUL, NSIM)
    USE mopenacc_struct_def
    USE wave
    IMPLICIT NONE

    TYPE (wavedes1) WDES1
    TYPE (wavefun1) W1(NSIM)
    INTEGER LMDIM, NSIM
    OVERLAP CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
         CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    REAL(q) :: EVALUE(NSIM)
    GDEF CRESUL(WDES1%NPROD,NSIM)
  ! local
    LOGICAL :: LDO(NSIM)
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NPRO2, NPRO2_, NT, NI, LMMAXC, L, LP, NP
    GDEF CTMP

    PROFILING_START('overl_mu')

!$ACC KERNELS PRESENT(CRESUL) __IF_ASYNC__
    CRESUL=0
!$ACC END KERNELS

    LDO=W1%LDO

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
    DO ISPINOR_=0,WDES1%NRSPINORS-1

       NPRO =ISPINOR *(WDES1%NPRO/2)
       NPRO_=ISPINOR_*(WDES1%NPRO/2)

!$ACC PARALLEL LOOP GANG COLLAPSE(force:2) VECTOR_LENGTH(32) __IF_ASYNC__ &
!$ACC PRIVATE(NT,LMMAXC,NPRO2,NPRO2_) PRESENT(WDES1,CDIJ,CQIJ,W1,CRESUL,EVALUE)
       DO NP=1,NSIM
          IF(.NOT.LDO(NP)) CYCLE
          DO NI=1,WDES1%NIONS
             NT=WDES1%ITYP(NI)
             LMMAXC=WDES1%LMMAX(NT)
             IF (LMMAXC==0) CYCLE
             NPRO2 =WDES1%LMBASE(NI)+NPRO
             NPRO2_=WDES1%LMBASE(NI)+NPRO_
             IF (EVALUE(NP)==0) THEN
!$ACC LOOP VECTOR PRIVATE(CTMP,LP)
                DO L =1,LMMAXC
                   CTMP=0
                   DO LP=1,LMMAXC
                      CTMP=CTMP + &
                           CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)*W1(NP)%CPROJ(LP+NPRO2_)
                   ENDDO
                   CRESUL(L+NPRO2,NP)=CRESUL(L+NPRO2,NP)+CTMP
                ENDDO
             ELSE
!$ACC LOOP VECTOR PRIVATE(CTMP,LP)
                DO L =1,LMMAXC
                   CTMP=0
                   DO LP=1,LMMAXC
                      CTMP=CTMP+ &
                          (CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)- &
                           EVALUE(NP)*CQIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)) * W1(NP)%CPROJ(LP+NPRO2_)
                   ENDDO
                   CRESUL(L+NPRO2,NP)=CRESUL(L+NPRO2,NP)+CTMP
                ENDDO
             ENDIF
          ENDDO
       ENDDO

    ENDDO
    ENDDO spinor

    PROFILING_STOP('overl_mu')

  END SUBROUTINE OVERL_MU
#endif

!
! identical to previous version but with complex EVALUE
!
!> @details @ref openmp :
!> the nested loops over \"types\" + \"ions-of-type\" are replaced
!> by a single loop over \"all ions\" that is distributed over
!> all available threads.
!
  SUBROUTINE OVERL1_C(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF,CRESUL)
#ifdef _OPENACC
    USE mopenacc_struct_def
#endif
    USE wave
    IMPLICIT NONE

    TYPE (wavedes1) WDES1
    INTEGER LMDIM
    OVERLAP CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
         CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    COMPLEX(q) EVALUE
    GDEF CRESUL(WDES1%NPRO),CPROF(WDES1%NPRO)
  ! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NPRO2, NPRO2_, NT, NI, LMMAXC, L, LP
    GDEF CTMP

    PROFILING_START('overl1_c')

!$ACC KERNELS PRESENT(CRESUL) __IF_ASYNC__
    CRESUL=0
!$ACC END KERNELS

!$ACC ENTER DATA COPYIN(EVALUE) __IF_ASYNC__

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
    DO ISPINOR_=0,WDES1%NRSPINORS-1

       NPRO =ISPINOR *(WDES1%NPRO/2)
       NPRO_=ISPINOR_*(WDES1%NPRO/2)
#ifdef _OPENACC
!$ACC PARALLEL LOOP GANG __IF_ASYNC__ &
!$ACC PRESENT(WDES1,CDIJ,CQIJ,CPROF,CRESUL,EVALUE) PRIVATE(NT,LMMAXC,NPRO2,NPRO2_,L,LP,CTMP)
#else
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP SHARED(NPRO,NPRO_,WDES1,EVALUE,ISPINOR,ISPINOR_,CDIJ,CQIJ,CPROF,CRESUL) &
!$OMP PRIVATE(NI,NT,LMMAXC,NPRO2,NPRO2_,L,LP,CTMP)
#endif
       DO NI=1,WDES1%NIONS
          NT=WDES1%ITYP(NI)
          LMMAXC=WDES1%LMMAX(NT)
          IF (LMMAXC==0) CYCLE
          NPRO2 =WDES1%LMBASE(NI)+NPRO
          NPRO2_=WDES1%LMBASE(NI)+NPRO_
          DO L =1,LMMAXC
             CTMP=0
!DIR$ IVDEP
!OCL NOVREC
!$ACC LOOP VECTOR REDUCTION(+:CTMP)
             DO LP=1,LMMAXC
                CTMP=CTMP+(CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)- &
                     EVALUE*CQIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)) * CPROF(LP+NPRO2_)
             ENDDO
             CRESUL(L+NPRO2)=CRESUL(L+NPRO2)+CTMP
          ENDDO
       ENDDO
#ifndef _OPENACC
!$OMP END PARALLEL DO
#endif
    ENDDO
    ENDDO spinor

!$ACC EXIT DATA DELETE(EVALUE) __IF_ASYNC__

    PROFILING_STOP('overl1_c')

  END SUBROUTINE OVERL1_C

!
! indentical to OVERL1 but with CDIJ and CQIJ always complex
!
  SUBROUTINE OVERL1_CCDIJ(WDES1, LMDIM, CDIJ, CQIJ, EVALUE, CPROF,CRESUL)
    USE wave
    IMPLICIT NONE

    TYPE (wavedes1) WDES1
    INTEGER LMDIM
    COMPLEX(q) CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
         CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    REAL(q) :: EVALUE
    GDEF CRESUL(WDES1%NPRO),CPROF(WDES1%NPRO)
  ! local
    INTEGER ISPINOR, ISPINOR_, NPRO, NPRO_, NT, NIS, NI, LMMAXC, L, LP

    CRESUL=0
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
    DO ISPINOR_=0,WDES1%NRSPINORS-1

       NPRO =ISPINOR *(WDES1%NPRO/2)
       NPRO_=ISPINOR_*(WDES1%NPRO/2)

       NIS =1
       DO NT=1,WDES1%NTYP
          LMMAXC=WDES1%LMMAX(NT)
          IF (LMMAXC/=0) THEN
             DO NI=NIS,WDES1%NITYP(NT)+NIS-1
                IF (EVALUE==0) THEN
                   DO L =1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
                   DO LP=1,LMMAXC
                      CRESUL(L+NPRO)=CRESUL(L+NPRO)+ &
                           CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)*CPROF(LP+NPRO_)
                   ENDDO
                   ENDDO
                ELSE
                   DO L =1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
                   DO LP=1,LMMAXC
                      CRESUL(L+NPRO)=CRESUL(L+NPRO)+ &
                           (CDIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)- &
                           EVALUE*CQIJ(LP,L,NI,ISPINOR_+2*ISPINOR+1)) * CPROF(LP+NPRO_)
                   ENDDO
                   ENDDO
                ENDIF
                NPRO = LMMAXC+NPRO
                NPRO_= LMMAXC+NPRO_
             ENDDO
          ENDIF
          NIS = NIS+WDES1%NITYP(NT)
       ENDDO
    ENDDO
    ENDDO spinor

  END SUBROUTINE OVERL1_CCDIJ
