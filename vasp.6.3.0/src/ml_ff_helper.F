#if defined(MPI) || defined(MPI_CHAIN)


!******************************************************************************************
! Module containing necessary arrays for lapack within the machine
! learning force field code.
!******************************************************************************************
      MODULE LAPACK_DATA

        USE ML_FF_PREC
        IMPLICIT NONE

        CONTAINS

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DSYEV_INIT
! This subroutine sets LWORK as -1 and allocate WORK for calculating the minimum size of 
! work arrays
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DSYEV_INIT(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           LAPACKWORK%LWORK=-1
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
           ALLOCATE(LAPACKWORK%WORK(1))
        END SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DSYEV_INIT

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DSYEV_FIN
! This subroutine deallocate work array that was used to determine the necessary size of 
! work arrays.
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DSYEV_FIN(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           LAPACKWORK%LWORK=INT(LAPACKWORK%WORK(1))
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
        END SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DSYEV_FIN

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_DSYEV_INIT
! This subroutine sets LWORK as -1 and allocate WORK for calculating the minimum size of 
! work arrays
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_DSYEV_INIT(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
           ALLOCATE(LAPACKWORK%WORK(LAPACKWORK%LWORK))
        END SUBROUTINE LAPACK_WORK_ARRAY_DSYEV_INIT

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_DSYEV_FIN
! This subroutine deallocate work array that was used to determine the necessary size of 
! work arrays.
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_DSYEV_FIN(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
        END SUBROUTINE LAPACK_WORK_ARRAY_DSYEV_FIN

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGESVD_INIT
! This subroutine sets LWORK as -1 and allocate WORK for calculating the minimum size of 
! work arrays
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGESVD_INIT(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           LAPACKWORK%LWORK=-1
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
           ALLOCATE(LAPACKWORK%WORK(1))
        END SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGESVD_INIT

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGESVD_FIN
! This subroutine deallocate work array that was used to determine the necessary size of 
! work arrays.
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGESVD_FIN(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           LAPACKWORK%LWORK=INT(LAPACKWORK%WORK(1))
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
        END SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGESVD_FIN

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_DGESVD_INIT
! This subroutine sets LWORK as -1 and allocate WORK for calculating the minimum size of 
! work arrays
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_DGESVD_INIT(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
           ALLOCATE(LAPACKWORK%WORK(LAPACKWORK%LWORK))
        END SUBROUTINE LAPACK_WORK_ARRAY_DGESVD_INIT

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_DSYEV_FIN
! This subroutine deallocate work array that was used to determine the necessary size of 
! work arrays.
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_DGESVD_FIN(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
        END SUBROUTINE LAPACK_WORK_ARRAY_DGESVD_FIN

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_DGETRF_INIT
! This subroutine sets LAPACKWORK%IPIV for LU factorization
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_DGETRF_INIT (LAPACKWORK,M)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           INTEGER, INTENT(IN) :: M
           IF (ALLOCATED(LAPACKWORK%IPIV)) THEN
             DEALLOCATE(LAPACKWORK%IPIV)
           ENDIF
           ALLOCATE(LAPACKWORK%IPIV(M))
        END SUBROUTINE LAPACK_WORK_ARRAY_DGETRF_INIT

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_DGETRF_FIN
! This subroutine deallocates LAPACKWORK%IPIV
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_DGETRF_FIN (LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           IF (ALLOCATED(LAPACKWORK%IPIV)) THEN
             DEALLOCATE(LAPACKWORK%IPIV)
           ENDIF
        END SUBROUTINE LAPACK_WORK_ARRAY_DGETRF_FIN

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGETRI_INIT
! This subroutine allocates arrays necessary for determining the optimal sizes of 
! working arrays for DGETRI
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGETRI_INIT(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           LAPACKWORK%LWORK=-1
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
           ALLOCATE(LAPACKWORK%WORK(1:1))
        END SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGETRI_INIT

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGETRI_FIN
! This subroutine deallocates arrays used for determining the optimal sizes of 
! working arrays for DGETRI
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGETRI_FIN(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           LAPACKWORK%LWORK=INT(LAPACKWORK%WORK(1))
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
        END SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGETRI_FIN

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_DGETRI_INIT
! This subroutine allocates arrays necessary for DGETRI
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_DGETRI_INIT(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
           ALLOCATE(LAPACKWORK%WORK(1:LAPACKWORK%LWORK))
        END SUBROUTINE LAPACK_WORK_ARRAY_DGETRI_INIT

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_DGETRI_FIN
! This subroutine deallocates arrays used in DGETRI
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_DGETRI_FIN(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
              DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
        END SUBROUTINE LAPACK_WORK_ARRAY_DGETRI_FIN

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGELSY_INIT
! This subroutine sets LWORK as -1 and allocate WORK for calculating the minimum size of 
! work arrays
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGELSY_INIT(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           LAPACKWORK%LWORK=-1
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
           ALLOCATE(LAPACKWORK%WORK(1))
        END SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGELSY_INIT

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGELSY_FIN
! This subroutine deallocate work array that was used to determine the necessary size of 
! work arrays.
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGELSY_FIN(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           LAPACKWORK%LWORK=INT(LAPACKWORK%WORK(1))
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
        END SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGELSY_FIN

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_DGELSY_INIT
! This subroutine sets LWORK as -1 and allocate WORK for calculating the minimum size of 
! work arrays
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_DGELSY_INIT(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
           ALLOCATE(LAPACKWORK%WORK(LAPACKWORK%LWORK))
        END SUBROUTINE LAPACK_WORK_ARRAY_DGELSY_INIT

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_DGELSY_FIN
! This subroutine deallocate work array that was used to determine the necessary size of 
! work arrays.
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_DGELSY_FIN(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
        END SUBROUTINE LAPACK_WORK_ARRAY_DGELSY_FIN

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGELSS_INIT
! This subroutine sets LWORK as -1 and allocate WORK for calculating the minimum size of 
! work arrays
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGELSS_INIT(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           LAPACKWORK%LWORK=-1
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
           ALLOCATE(LAPACKWORK%WORK(1))
        END SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGELSS_INIT

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGELSS_FIN
! This subroutine deallocate work array that was used to determine the necessary size of 
! work arrays.
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGELSS_FIN(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           LAPACKWORK%LWORK=INT(LAPACKWORK%WORK(1))
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
        END SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGELSS_FIN

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_DGELSS_INIT
! This subroutine sets LWORK as -1 and allocate WORK for calculating the minimum size of 
! work arrays
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_DGELSS_INIT(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
           ALLOCATE(LAPACKWORK%WORK(LAPACKWORK%LWORK))
        END SUBROUTINE LAPACK_WORK_ARRAY_DGELSS_INIT

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_DGELSS_FIN
! This subroutine deallocate work array that was used to determine the necessary size of 
! work arrays.
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_DGELSS_FIN(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
        END SUBROUTINE LAPACK_WORK_ARRAY_DGELSS_FIN

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGELSD_INIT
! This subroutine sets LWORK as -1 and allocate WORK for calculating the minimum size of 
! work arrays
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGELSD_INIT(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           LAPACKWORK%LWORK=-1
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
           ALLOCATE(LAPACKWORK%WORK(1))
           IF (ALLOCATED(LAPACKWORK%IPIV)) THEN
             DEALLOCATE(LAPACKWORK%IPIV)
           ENDIF
           ALLOCATE(LAPACKWORK%IPIV(1))
        END SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGELSD_INIT

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGELSD_FIN
! This subroutine deallocate work array that was used to determine the necessary size of 
! work arrays.
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGELSD_FIN(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           LAPACKWORK%LWORK=INT(LAPACKWORK%WORK(1))
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
        END SUBROUTINE LAPACK_WORK_ARRAY_SIZE_DGELSD_FIN

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_DGELSD_INIT
! This subroutine sets LWORK as -1 and allocate WORK for calculating the minimum size of 
! work arrays
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_DGELSD_INIT(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           INTEGER LIWORK
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
           ALLOCATE(LAPACKWORK%WORK(LAPACKWORK%LWORK))
           LIWORK = INT(LAPACKWORK%IPIV(1))
           IF (ALLOCATED(LAPACKWORK%IPIV)) THEN
             DEALLOCATE(LAPACKWORK%IPIV)
           ENDIF
           ALLOCATE(LAPACKWORK%IPIV(LIWORK))
        END SUBROUTINE LAPACK_WORK_ARRAY_DGELSD_INIT

!******************************************************************************************
! SUBROUTINE LAPACK_WORK_ARRAY_DGELD_FIN
! This subroutine deallocate work array that was used to determine the necessary size of 
! work arrays.
!******************************************************************************************
        SUBROUTINE LAPACK_WORK_ARRAY_DGELSD_FIN(LAPACKWORK)
           USE ML_FF_STRUCT, ONLY: LAPACK_WORK
           IMPLICIT NONE
           TYPE (LAPACK_WORK) LAPACKWORK
           IF (ALLOCATED(LAPACKWORK%WORK)) THEN
             DEALLOCATE(LAPACKWORK%WORK)
           ENDIF
           IF (ALLOCATED(LAPACKWORK%IPIV)) THEN
             DEALLOCATE(LAPACKWORK%IPIV)
           ENDIF
        END SUBROUTINE LAPACK_WORK_ARRAY_DGELSD_FIN


      END MODULE LAPACK_DATA


!****************************************************************************************************
! Module to measure the time. Definitions of the times are listed below.
! 1 : Time for reading files, setting parameters and so on
! 2 : Time for calculating Coulomb descriptor and design matrix
! 3 : Time for calculating pairwise descriptor and design matrix
! 4 : Time for calculating many-body descriptor and design matrix
! 5 : Time for sparsification
! 6 : Time for regression of Coulomb interactions
! 7 : Time for regression of pairwise interactions
! 8 : Time for regression of many-body interactions
! 9 : Time for predicting Coulomb interactions
! 10: Time for predicting pairwise interactions
! 11: Time for predicting many-body interactions
! 12: Total time
!****************************************************************************************************

      MODULE TIME_CHECK
        USE ML_FF_PREC
        IMPLICIT NONE

        CONTAINS

!****************************************************************************************************
! Initialization of times
!****************************************************************************************************

        SUBROUTINE INIT_TIME(TIM)
        USE ML_FF_STRUCT, ONLY : TIME
        IMPLICIT NONE
        TYPE (TIME) :: TIM
        TIM%CPUTIM=0.0_q
        TIM%VPUTIM=0.0_q
        END SUBROUTINE INIT_TIME

!****************************************************************************************************
! Timing routine
!****************************************************************************************************
        SUBROUTINE VTIME(VPUTIME,CPUTIME)
        IMPLICIT NONE
        REAL(q) :: VPUTIME
        REAL(q) :: CPUTIME
        INTEGER :: TIME_COUNT
        INTEGER :: TIME_COUNT_RATE
        INTEGER :: TIME_COUNT_MAX
        CALL CPU_TIME(CPUTIME)
        CALL SYSTEM_CLOCK(TIME_COUNT,TIME_COUNT_RATE,TIME_COUNT_MAX)
        VPUTIME=REAL(TIME_COUNT,q)/REAL(TIME_COUNT_RATE,q)
        END SUBROUTINE

      END MODULE TIME_CHECK


!****************************************************************************************************
! This module contains all the random number generators
!****************************************************************************************************
      MODULE RANDOM_NUMBER_GENERATOR
         USE ML_FF_PREC
         IMPLICIT NONE
 
         CONTAINS
        
         SUBROUTINE RMARIN_ION(IJ,KL,U,C,CD,CM,I97,J97,TEST)
         USE ML_FF_PREC
         USE PRINT2TERMINAL
         IMPLICIT NONE
         INTEGER :: I,J,K,L,M,IJ,KL,II,JJ
         REAL(q) :: S,T
         REAL(q) U(:) !(97)
         LOGICAL:: TEST
         REAL(q) ::  C, CD, CM
         INTEGER :: I97, J97

         IF ( (IJ<0) .OR. (IJ>31328) .OR. &
            (KL<0) .OR. (KL>30081) ) THEN
            CALL ERROR_PRINT2TERMINAL ("The first random number seed must have a value between 0 and 31328 \n The &
                 second seed must have a value between 0 and 30081")
         ENDIF

         I = MOD(IJ/177, 177) + 2
         J = MOD(IJ    , 177) + 2
         K = MOD(KL/169, 178) + 1
         L = MOD(KL,     169)

         DO II = 1, 97
            S = 0.0_q
            T = 0.5_q
            DO jj = 1, 24
               M = MOD(MOD(I*J, 179)*K, 179)
               I = J
               J = K
               K = M
               L = MOD(53*L+1, 169)
               IF (MOD(L*M, 64) >= 32) THEN
                  S = S + T
               ENDIF
               T = 0.5_q * T
            ENDDO
            U(II) = S
         ENDDO

         C = 362436.0_q / 16777216.0_q
         CD = 7654321.0_q / 16777216.0_q
         CM = 16777213.0_q /16777216.0_q

         I97 = 97
         J97 = 33

         TEST = .TRUE.
         END SUBROUTINE RMARIN_ION

         SUBROUTINE RANMAR_ION(RDUMMY,U, C, CD, CM, I97, J97, TEST)
         USE ML_FF_PREC
         IMPLICIT NONE
         REAL(q) U(:) !(97)
         REAL(q) C, CD, CM, RDUMMY
         REAL(q) :: UNI
         INTEGER I97, J97
         LOGICAL TEST

         INTEGER IVEC

         IF (.NOT. TEST) THEN
            WRITE(*,*) 'Call the init routine (RMARIN_ION) &
                         before calling RANMAR_ION!'
            WRITE(*,*) 'Initializing now with built-in seeds &
                        1802 and 9373 ...'
            CALL RMARIN_ION(1802,9373,U,C,CD,CM,I97,J97,TEST)
         ENDIF

         UNI = U(I97) - U(J97)
         IF ( UNI < 0.0_q ) UNI = UNI + 1.0_q
         U(I97) = UNI
         I97 = I97 - 1
         IF (I97 == 0) I97 = 97
         J97 = J97 - 1
         IF (J97 == 0) J97 = 97
         C = C - CD
         IF ( C < 0.0_q ) C = C + CM
         UNI = UNI - C
         IF ( UNI < 0.0_q ) UNI = UNI + 1.0_q
         RDUMMY = UNI
         END SUBROUTINE RANMAR_ION

         SUBROUTINE RANE_ION(RNG_ION,PUT,GET)
         USE ML_FF_PREC
         IMPLICIT NONE
         INTEGER, OPTIONAL :: PUT(:),GET(:)
         INTEGER,PARAMETER :: IJ_MAX=31328
         INTEGER,PARAMETER ::KL_MAX=30081
         INTEGER,PARAMETER :: SEED1_MAX=900000000
         INTEGER,SAVE :: RNG_STATE1=54217137 !IJ*(KL_MAX+1)+KL
         INTEGER,SAVE :: RNG_STATE2=0            !c number of rannd. numbers generated
         INTEGER,SAVE :: RNG_STATE3=0            !c up to this point is:
         INTEGER,PARAMETER :: COUNTLARGE=1000000 !c RNG_STATE3*COUNTLARGE+RNG_STATE2
         !INTEGER,PARAMETER :: COUNTLARGE=1000 !c RNG_STATE3*COUNTLARGE+RNG_STATE2
         REAL(q) :: rng_ion,rdummy
         INTEGER :: IMODE,i,j
         LOGICAL,SAVE :: LFIRST=.TRUE.,LINIT=.FALSE.,TEST=.FALSE.
         REAL(q),SAVE :: U(97)
         REAL(q),SAVE :: C, CD, CM
         INTEGER,SAVE :: I97, J97
         INTEGER, SAVE :: IJ, KL,ICALL

         IF (LFIRST) THEN
            ICALL=0
            IJ=1802
            KL=9373
            C = 362436.0_q / 16777216.0_q
            CD = 7654321.0_q / 16777216.0_q
            CM = 16777213.0_q /16777216.0_q
            I97 = 97
            J97 = 33
            TEST=.FALSE.
            LINIT=.TRUE.
         ENDIF

         !c just return the RNG seed (and initialize if needed)
         IF (PRESENT(GET)) THEN
            IMODE=2
            GET(1)=RNG_STATE1
            GET(2)=RNG_STATE2
            GET(3)=RNG_STATE3
            !write(*,*) "get?"
            IF ((.NOT. LFIRST)) RETURN
         !c reseed and reinitialise if SEED is provided
         ELSE IF (PRESENT(PUT)) THEN
            IMODE=1
            !write(*,*) 'put?',SEED
            !c ensure that SEED(1) doesn't exceed
            !c max. sensible value
            RNG_STATE1=MOD(PUT(1),SEED1_MAX)
            RNG_STATE2=MOD(PUT(2),COUNTLARGE)
            RNG_STATE3=PUT(3)+PUT(2)/COUNTLARGE

            !c regenerate pair of seeds used originally
            IJ=RNG_STATE1/(KL_MAX+1)
            KL=RNG_STATE1-(KL_MAX+1)*IJ
            LINIT=.TRUE.
         ELSE
            IMODE=0
         ENDIF

         IF (LINIT) THEN
            CALL RMARIN_ION(IJ,KL,U,C,CD,CM,I97,J97,TEST)

            !c skip RNG_STATE3*COUNTLONG+RNG_STATE2 numbers 
            !c (presumably already used in a previous run) 
            !c on initialization
            DO i=1,RNG_STATE3
               DO j=1, COUNTLARGE
                  CALL RANMAR_ION(rdummy,U, C, CD, CM, I97, J97, TEST)
               ENDDO
            ENDDO
            DO i=1,RNG_STATE2
               CALL RANMAR_ION(rdummy,U, C, CD, CM, I97, J97, TEST)
            ENDDO
            LINIT=.FALSE.
            !write(*,*) "rng_init_finish"
         ENDIF

         !c return random number if no special action is taken
         IF (IMODE==0) THEN
            ! write(*,*) "rng_num?",RNG_STATE1,RNG_STATE2
            RNG_STATE2=RNG_STATE2+1
            IF (RNG_STATE2>=COUNTLARGE) THEN
               RNG_STATE2=RNG_STATE2-COUNTLARGE
               RNG_STATE3=RNG_STATE3+1
            ENDIF
            CALL RANMAR_ION(RNG_ION,U, C, CD, CM, I97, J97, TEST)
         ENDIF

         IF (LFIRST) LFIRST=.FALSE.
         END SUBROUTINE RANE_ION
      
      END MODULE RANDOM_NUMBER_GENERATOR

#ifdef ml_debug
!-----------------------------------------------------------------
! The following contains subroutines to write out important 
! distributed and non-distributed arrays. The data is written out
! to a file, where the offset of the unit is defined by IUNIT_INPUT
!-----------------------------------------------------------------
      MODULE MLFF_TESTING
         USE ML_FF_PREC
         USE MPI_DATA
         IMPLICIT NONE
      
         CONTAINS 
       
         SUBROUTINE PRINT_ORDERED_FMAT(PAR_SUP_HANDLE,FF,FFM,IUNIT_INPUT)
            USE ML_FF_STRUCT, ONLY: PARALLEL_SUPER, FF_PAR, FFM_PAR, ML_IO_WRITE
            IMPLICIT NONE
            ! Input variables
            TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
            TYPE (FF_PAR)         :: FF
            TYPE (FFM_PAR)        :: FFM
            INTEGER               :: IUNIT_INPUT
            ! Local variables
            INTEGER :: IHELP
            INTEGER :: IUNIT
            INTEGER :: JHELP
            INTEGER :: JJMAX
            INTEGER :: JJMIN
            INTEGER :: KHELP
            INTEGER :: KKMAX
            INTEGER :: KKMIN
            INTEGER :: ORDERED_REF_CONF(1:FFM%MBTOT,1:FF%MTYP)
       
            ! make reverse list for writing out in an ordered way
            CALL MAKE_ORDERED_LIST(PAR_SUP_HANDLE,FF,FFM,ORDERED_REF_CONF)
 
            ! write fmat
            DO IHELP=1,FF%MTYP
               IUNIT=IUNIT_INPUT+IHELP
               CALL SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_FMAT_FFM,FFM%FMAT(:,:,IHELP),PAR_SUP_HANDLE%CONTEXT_WORLD)
               IF (ML_IO_WRITE) THEN
                  KKMIN=1
                  KKMAX=FFM%NDATA
                  JJMIN=1
                  JJMAX=FFM%NB(IHELP)
                  IF (IHELP.EQ.1) THEN
                  ELSEIF (IHELP.EQ.2) THEN
                  ELSE
                  ENDIF
                  DO JHELP=JJMIN,JJMAX
                     DO KHELP=KKMIN,KKMAX
                        write(IUNIT,*) IHELP, JHELP, ORDERED_REF_CONF(JHELP,IHELP), KHELP, PAR_SUP_HANDLE%GENERAL_SINGLE(ORDERED_REF_CONF(JHELP,IHELP),KHELP)
                     ENDDO
                  ENDDO
               ENDIF
               CALL FINIT_SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE)
            ENDDO
         END SUBROUTINE PRINT_ORDERED_FMAT

         SUBROUTINE PRINT_YMAT(PAR_SUP_HANDLE,FF,FFM,IUNIT_INPUT)
            USE ML_FF_STRUCT, ONLY: PARALLEL_SUPER, FF_PAR, FFM_PAR, ML_IO_WRITE
            IMPLICIT NONE
            ! Input variables
            TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
            TYPE (FF_PAR)         :: FF
            TYPE (FFM_PAR)        :: FFM
            INTEGER               :: IUNIT_INPUT
            ! Local variables
            INTEGER :: IHELP
            INTEGER :: IUNIT
            INTEGER :: JHELP
            INTEGER :: JJMAX
            INTEGER :: JJMIN
            INTEGER :: KHELP
 
            ! Make IO unit
            IUNIT=IUNIT_INPUT
            ! Reduce YMAT to one core
            CALL SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE,PAR_SUP_HANDLE%DESC_YMAT_FFM,FFM%YMAT(:,:),PAR_SUP_HANDLE%CONTEXT_WORLD)
            IF (ML_IO_WRITE) THEN
               DO IHELP=1,FFM%NDATA
                  write(IUNIT,*) IHELP, FFM%NDATA, FFM%NDATA, FFM%NDATA, PAR_SUP_HANDLE%GENERAL_SINGLE(IHELP,1)
               ENDDO
            ENDIF
            CALL FINIT_SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE)
         END SUBROUTINE PRINT_YMAT

         SUBROUTINE PRINT_3D_ARRAY_DISTRIBUTED(ABN,PAR_SUP_HANDLE,FF,FFM,IUNIT_INPUT,ARRAY_3D,DESC_ARRAY_3D,INTYP0)
            USE ML_FF_STRUCT, ONLY: ABINITIONEW, PARALLEL_SUPER, FF_PAR, FFM_PAR, PARALLEL_DESCRIPTOR, ML_IO_WRITE
            IMPLICIT NONE
            ! Input variables
            TYPE (ABINITIONEW)         :: ABN
            TYPE (PARALLEL_SUPER)      :: PAR_SUP_HANDLE
            TYPE (FF_PAR)              :: FF
            TYPE (FFM_PAR)             :: FFM
            INTEGER                    :: IUNIT_INPUT
            REAL(q)                    :: ARRAY_3D(:,:,:)
            TYPE (PARALLEL_DESCRIPTOR) :: DESC_ARRAY_3D
            INTEGER                    :: INTYP0
            ! Local variables
            INTEGER :: IB
            INTEGER :: ICOUNT
            INTEGER :: IHELP
            INTEGER :: IERR
            INTEGER :: INIONS
            INTEGER :: INITYP
            INTEGER :: INITYP_TMP
            INTEGER :: IUNIT
            INTEGER :: JHELP
            INTEGER :: JJMAX
            INTEGER :: JJMIN
            INTEGER :: KHELP
            INTEGER :: KKMAX
            INTEGER :: KKMIN
            INTEGER :: LOCVAR_COL
            INTEGER :: ORDERED_REF_CONF(1:FFM%MBTOT,1:FF%MTYP)
            REAL(q) :: SOAP_ALL(1:FF%MITYP,1:FFM%MBTOT)
            REAL(q) :: SOAP_ALL2(1:FF%MITYP,1:FFM%MBTOT)

            ! make reverse list for writing out in an ordered way
            CALL MAKE_ORDERED_LIST(PAR_SUP_HANDLE,FF,FFM,ORDERED_REF_CONF)

            IUNIT=IUNIT_INPUT+INTYP0
            CALL SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE,DESC_ARRAY_3D,ARRAY_3D(:,:,INTYP0),PAR_SUP_HANDLE%CONTEXT_WORLD)
            IF (ML_IO_WRITE) THEN
               KKMIN=1
               KKMAX=FFM%NNVAR2(INTYP0)
               JJMIN=1
               JJMAX=FFM%NB(INTYP0)
               DO JHELP=JJMIN,JJMAX
                  DO KHELP=KKMIN,KKMAX
                     write(IUNIT,*) JHELP, KHELP, ORDERED_REF_CONF(JHELP,INTYP0), INTYP0, PAR_SUP_HANDLE%GENERAL_SINGLE(KHELP,ORDERED_REF_CONF(JHELP,INTYP0))
                  ENDDO
               ENDDO
            ENDIF
            CALL FINIT_SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE)
         END SUBROUTINE PRINT_3D_ARRAY_DISTRIBUTED

         SUBROUTINE MAKE_ORDERED_LIST(PAR_SUP_HANDLE,FF,FFM,ORDERED_REF_CONF)
            USE ML_FF_STRUCT, ONLY: PARALLEL_SUPER, FF_PAR, FFM_PAR, ML_IO_WRITE
            IMPLICIT NONE
            ! Input variables
            TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
            TYPE (FF_PAR)         :: FF
            TYPE (FFM_PAR)        :: FFM
            ! Input and output variables
            INTEGER :: ORDERED_REF_CONF(:,:)
            ! Local variables
            INTEGER :: CURRENT_CONF
            INTEGER :: CURRENT_ION
            INTEGER :: CURRENT_NB
            INTEGER :: ICOUNT
            INTEGER :: IHELP
            INTEGER :: JHELP
            INTEGER :: KHELP
            INTEGER :: LHELP
            LOGICAL :: L_INLIST
            DO IHELP=1,FF%MTYP
               ! Find the first element of list
               ICOUNT=0
               CURRENT_CONF=FFM%LNCONF(1,IHELP)
               CURRENT_ION=FFM%LNIONS(1,IHELP)
               CURRENT_NB=1
               DO JHELP=1,FFM%NB(IHELP)
                  IF (FFM%LNCONF(JHELP,IHELP).LT.CURRENT_CONF) THEN
                     CURRENT_CONF=FFM%LNCONF(JHELP,IHELP)
                     CURRENT_ION=FFM%LNIONS(JHELP,IHELP)
                     CURRENT_NB=JHELP
                  ELSEIF (FFM%LNCONF(JHELP,IHELP).EQ.CURRENT_CONF) THEN
                     IF (FFM%LNIONS(JHELP,IHELP).LT.CURRENT_ION) THEN
                        CURRENT_CONF=FFM%LNCONF(JHELP,IHELP)
                        CURRENT_ION=FFM%LNIONS(JHELP,IHELP)
                        CURRENT_NB=JHELP
                     ENDIF
                  ENDIF
               ENDDO
               ICOUNT=ICOUNT+1
               ORDERED_REF_CONF(ICOUNT,IHELP)=CURRENT_NB
               ! Do the rest of the list
               DO JHELP=2, FFM%NB(IHELP)
                  ! Find next starting entry
                  DO KHELP=1,FFM%NB(IHELP)
                     L_INLIST=.FALSE.
                     DO LHELP=1,ICOUNT
                        IF (ORDERED_REF_CONF(LHELP,IHELP).EQ.KHELP) THEN
                           L_INLIST=.TRUE.
                        ENDIF
                     ENDDO
                     IF (.NOT.L_INLIST) THEN
                        CURRENT_CONF=FFM%LNCONF(KHELP,IHELP)
                        CURRENT_ION=FFM%LNIONS(KHELP,IHELP)
                        CURRENT_NB=KHELP
                        EXIT
                     ENDIF
                  ENDDO
                  ! check if any of the other elements is smaller
                  DO KHELP=1,FFM%NB(IHELP)
                     L_INLIST=.FALSE.
                     DO LHELP=1,ICOUNT
                        IF (ORDERED_REF_CONF(LHELP,IHELP).EQ.KHELP) THEN
                           L_INLIST=.TRUE.
                        ENDIF
                     ENDDO
                     IF (L_INLIST) THEN
                        CYCLE
                     ELSE
                        IF (FFM%LNCONF(KHELP,IHELP).LT.CURRENT_CONF) THEN
                           CURRENT_CONF=FFM%LNCONF(KHELP,IHELP)
                           CURRENT_ION=FFM%LNIONS(KHELP,IHELP)
                           CURRENT_NB=KHELP
                        ELSEIF (FFM%LNCONF(KHELP,IHELP).EQ.CURRENT_CONF) THEN
                           IF (FFM%LNIONS(KHELP,IHELP).LT.CURRENT_ION) THEN
                              CURRENT_CONF=FFM%LNCONF(KHELP,IHELP)
                              CURRENT_ION=FFM%LNIONS(KHELP,IHELP)
                              CURRENT_NB=KHELP
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
                  ICOUNT=ICOUNT+1
                  ORDERED_REF_CONF(ICOUNT,IHELP)=CURRENT_NB
               ENDDO
            ENDDO
         END SUBROUTINE MAKE_ORDERED_LIST

      END MODULE MLFF_TESTING
#endif
#endif
