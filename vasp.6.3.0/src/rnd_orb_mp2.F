!#define dotiming
#include "symbol.inc"

#ifdef MPI



!***********************************************************************
!
! this module provides all necessary statistical calculations like
! mean value, variance, etc... 
!
!***********************************************************************
MODULE rnd_stat
  USE prec
  IMPLICIT none

  INTEGER(qi8), PRIVATE :: N_DI_SAMPLES             ! number of direct sample
  INTEGER(qi8), PRIVATE :: N_EX_SAMPLES             ! number of exchange sample
  INTEGER, PRIVATE :: NENCUTS                    ! number of cutoffs
  REAL(q), PRIVATE, ALLOCATABLE :: DI_MEAN(:)    ! mean value of direct samples for each cutoff
  REAL(q), PRIVATE, ALLOCATABLE :: EX_MEAN(:)    ! mean value of exchange samples for each cutoff
  REAL(q), PRIVATE :: DI_VAR                     ! var of direct samples
  REAL(q), PRIVATE :: EX_VAR                     ! var of exchange samples
  REAL(q), PRIVATE :: DI_ERR                     ! error of direct samples
  REAL(q), PRIVATE :: EX_ERR                     ! error of exchange samples
  REAL(q), PRIVATE, ALLOCATABLE :: DI_MP2(:,:,:)   ! total direct energy for each cutoff and tau and q point
  REAL(q), PRIVATE, ALLOCATABLE :: EX_MP2(:,:,:)   ! total exchange energy for each cutoff and tau and q point
  REAL(q), PRIVATE, ALLOCATABLE :: DI_ERRS(:,:,:)  ! total direct errors for each cutoff and tau and q point
  REAL(q), PRIVATE, ALLOCATABLE :: EX_ERRS(:,:,:)  ! total exchange errors for each cutoff and tau and q point
  INTEGER, PRIVATE :: MAX_TAU                    ! total number of tau points
  INTEGER, PRIVATE :: MAX_Q                      ! total number of q-pointtotal number of q-pointss
  INTEGER, PRIVATE :: CURRENT_TAU
  INTEGER, PRIVATE :: CURRENT_Q
  ! file units
  INTEGER, PRIVATE :: DI_MEAN_FILEUNIT 
  INTEGER, PRIVATE :: EX_MEAN_FILEUNIT 
  INTEGER, PRIVATE :: DI_ERR_FILEUNIT
  INTEGER, PRIVATE :: EX_ERR_FILEUNIT
  INTEGER, PRIVATE :: DI_VAR_FILEUNIT
  INTEGER, PRIVATE :: EX_VAR_FILEUNIT

  CONTAINS
  
  !***************************** INIT_STATISTICS *******************************
  ! initializes the module (mem allocation)
  !*****************************************************************************
  SUBROUTINE INIT_STATISTICS(SET_MAX_TAU, SET_MAX_Q, N_MP2_ENCUTS)
    USE ini
    INTEGER :: SET_MAX_TAU
    INTEGER :: SET_MAX_Q
    INTEGER :: N_MP2_ENCUTS

    NENCUTS = N_MP2_ENCUTS

    MAX_TAU = SET_MAX_TAU
    MAX_Q = SET_MAX_Q

    ALLOCATE( DI_MEAN(NENCUTS), &
              EX_MEAN(NENCUTS), &
              DI_MP2(NENCUTS, MAX_Q, MAX_TAU),   &
              EX_MP2(NENCUTS, MAX_Q, MAX_TAU),   &
              DI_ERRS(NENCUTS, MAX_Q, MAX_TAU),   &
              EX_ERRS(NENCUTS, MAX_Q, MAX_TAU) )

    DI_MP2 = 0
    EX_MP2 = 0
    DI_ERRS = 0
    EX_ERRS = 0
  END SUBROUTINE INIT_STATISTICS 



  !***************************** DEALLOCATE_STATISTICS  *******************************
  ! for mem deallocation
  !*******************************************************************************
  SUBROUTINE DEALLOCATE_STATISTICS()  
    DEALLOCATE(DI_MEAN, EX_MEAN, DI_MP2, EX_MP2, DI_ERRS, EX_ERRS)
  END SUBROUTINE DEALLOCATE_STATISTICS  



  !***************************** FLUSH_STATISTICS  *******************************
  ! 
  ! Set q and tau-point and set all variables to 0 (except DI_MP2 and EX_MP2)
  !
  !*******************************************************************************
  SUBROUTINE FLUSH_STATISTICS(TAU, Q)  
    INTEGER :: TAU
    INTEGER :: Q
    CHARACTER(LEN = 5) :: NQ_STR
    CHARACTER(LEN = 5) :: NTAU_STR

    CURRENT_TAU = TAU
    CURRENT_Q = Q

    DI_MEAN = 0 
    EX_MEAN = 0 
    DI_VAR  = 0 
    EX_VAR  = 0 
    DI_ERR  = 0 
    EX_ERR  = 0 
    N_DI_SAMPLES = 0
    N_EX_SAMPLES = 0

    DI_MEAN_FILEUNIT = 421
    EX_MEAN_FILEUNIT = 422
    DI_ERR_FILEUNIT = 431
    EX_ERR_FILEUNIT = 432
    DI_VAR_FILEUNIT = 441
    EX_VAR_FILEUNIT = 442
    CLOSE(DI_MEAN_FILEUNIT)
    CLOSE(EX_MEAN_FILEUNIT)
    CLOSE(DI_ERR_FILEUNIT)
    CLOSE(EX_ERR_FILEUNIT)
    CLOSE(DI_VAR_FILEUNIT)
    CLOSE(EX_VAR_FILEUNIT)
    WRITE(NQ_STR, '(I5.5)') CURRENT_Q
    WRITE(NTAU_STR, '(I5.5)') CURRENT_TAU
    OPEN(UNIT = DI_MEAN_FILEUNIT, FILE = "di_mp2_NQ" // TRIM(NQ_STR) // "_TAU" // TRIM(NTAU_STR), STATUS = "REPLACE")
    OPEN(UNIT = EX_MEAN_FILEUNIT, FILE = "ex_mp2_NQ" // TRIM(NQ_STR) // "_TAU" // TRIM(NTAU_STR), STATUS = "REPLACE")
    OPEN(UNIT = DI_ERR_FILEUNIT, FILE = "di_err_NQ" // TRIM(NQ_STR) // "_TAU" // TRIM(NTAU_STR), STATUS = "REPLACE")
    OPEN(UNIT = EX_ERR_FILEUNIT, FILE = "ex_err_NQ" // TRIM(NQ_STR) // "_TAU" // TRIM(NTAU_STR), STATUS = "REPLACE")
    OPEN(UNIT = DI_VAR_FILEUNIT, FILE = "di_var_NQ" // TRIM(NQ_STR) // "_TAU" // TRIM(NTAU_STR), STATUS = "REPLACE")
    OPEN(UNIT = EX_VAR_FILEUNIT, FILE = "ex_var_NQ" // TRIM(NQ_STR) // "_TAU" // TRIM(NTAU_STR), STATUS = "REPLACE")
  END SUBROUTINE FLUSH_STATISTICS  



  !***************************** STATISTICS_TO_OUTCAR  *******************************
  ! 
  ! write statistics to OUTCAR file
  !
  !*******************************************************************************
  SUBROUTINE STATISTICS_TO_OUTCAR(IU6)  
    INTEGER :: IU6

    WRITE(IU6,*) ""
    WRITE(IU6,'(A,I5,A,I3)') " MP2 statistics for NQ=",CURRENT_Q, " and TAU=", CURRENT_TAU 
    WRITE(IU6,'(A)')         "   direct term"
    WRITE(IU6,'(A,F16.6)')   "     sample average     :", DI_MEAN(1) 
    WRITE(IU6,'(A,F16.6)')   "     statistical error  :", DI_ERR
    WRITE(IU6,'(A,F16.6)')   "     sample variance    :", DI_VAR
    WRITE(IU6,'(A,I12)')     "     sample size        :", N_DI_SAMPLES
    WRITE(IU6,'(A)')         "   exchange term"
    WRITE(IU6,'(A,F16.6)')   "     sample average     :", EX_MEAN(1) 
    WRITE(IU6,'(A,F16.6)')   "     statistical error  :", EX_ERR
    WRITE(IU6,'(A,F16.6)')   "     sample variance    :", EX_VAR
    WRITE(IU6,'(A,I12)')     "     sample size        :", N_EX_SAMPLES
  END SUBROUTINE STATISTICS_TO_OUTCAR


  !***************************** CLOSE_STATISTICS  *******************************
  ! closes the output file
  !*******************************************************************************
  SUBROUTINE CLOSE_STATISTICS(IU0, IU6, ENCUTS, ENCUTSOFTS)  
    INTEGER :: IU0, IU6
    REAL(q) :: ENCUTS(NENCUTS)
    REAL(q) :: ENCUTSOFTS(NENCUTS)
    ! local
    REAL(q) :: TOT_DI_MP2(NENCUTS), TOT_EX_MP2(NENCUTS), TOT_MP2(NENCUTS)
    REAL(q) :: TOT_DI_ERR, TOT_EX_ERR
    INTEGER :: NTAU, NQ, I
    REAL(q) :: EXTRAPOLATED_DIRECT_MP2
    REAL(q) :: EXTRAPOLATED_EXCHANGE_MP2

    TOT_DI_MP2 = 0
    TOT_EX_MP2 = 0
    TOT_DI_ERR = 0
    TOT_EX_ERR = 0

    DO NTAU = 1, MAX_TAU
      DO NQ = 1, MAX_Q
        TOT_DI_MP2(:) = TOT_DI_MP2(:) + DI_MP2(:, NQ, NTAU)
        TOT_EX_MP2(:) = TOT_EX_MP2(:) + EX_MP2(:, NQ, NTAU)
        TOT_DI_ERR = TOT_DI_ERR + DI_ERRS(1, NQ, NTAU)**2
        TOT_EX_ERR = TOT_EX_ERR + EX_ERRS(1, NQ, NTAU)**2
      ENDDO
    ENDDO
    TOT_DI_ERR = SQRT(TOT_DI_ERR)
    TOT_EX_ERR = SQRT(TOT_EX_ERR)

    TOT_MP2 = TOT_DI_MP2 + TOT_EX_MP2

    WRITE(IU0,'(A)') ""
    WRITE(IU0,'(A)') " MP2 result"
    WRITE(IU0,'(A)') " =========="
    WRITE(IU0,'(A,F16.6,A,F16.6)') "   direct MP2 energy  : ", TOT_DI_MP2(1), " +/- ", TOT_DI_ERR
    WRITE(IU0,'(A,F16.6,A,F16.6)') "   exchange MP2 energy: ", TOT_EX_MP2(1), " +/- ", TOT_EX_ERR
    WRITE(IU0,'(A,F16.6,A,F16.6)') "   total MP2 energy   : ", TOT_MP2(1), " +/- ", SQRT(TOT_DI_ERR**2+TOT_EX_ERR**2)
    WRITE(IU0,'(A)') ""

    WRITE(IU6,'(A)') ""
    WRITE(IU6,'(A)') " MP2 result"
    WRITE(IU6,'(A)') " =========="
    WRITE(IU6,'(A,F16.6,A,F16.6)') "   direct MP2 energy  : ", TOT_DI_MP2(1), " +/- ", TOT_DI_ERR
    WRITE(IU6,'(A,F16.6,A,F16.6)') "   exchange MP2 energy: ", TOT_EX_MP2(1), " +/- ", TOT_EX_ERR
    WRITE(IU6,'(A,F16.6,A,F16.6)') "   total MP2 energy   : ", TOT_MP2(1), " +/- ", SQRT(TOT_DI_ERR**2+TOT_EX_ERR**2)
    WRITE(IU6,'(A)') ""

    CALL LIN_REG_ENERGIES(TOT_DI_MP2, ENCUTS, EXTRAPOLATED_DIRECT_MP2)
    CALL LIN_REG_ENERGIES(TOT_EX_MP2, ENCUTS, EXTRAPOLATED_EXCHANGE_MP2)

    WRITE(IU6,*) " "
    WRITE(IU6,*) " -------------------------------------------------------------------------"
    WRITE(IU6,'(5A)') "      ENCUTGW", "  ENCUTGWSOFT", "      direct MP2", "    exchange MP2", "       total MP2"
    WRITE(IU6,*) " -------------------------------------------------------------------------"
    DO I = NENCUTS, 1, -1
      WRITE(IU6,'(2F13.3,3F16.8)') ENCUTS(I), ENCUTSOFTS(I), TOT_DI_MP2(I), TOT_EX_MP2(I), TOT_MP2(I)
    ENDDO
    WRITE(IU6,'(A)') "      linear regression"
    WRITE(IU6,'(A,3F16.8)') "      converged values    ", EXTRAPOLATED_DIRECT_MP2, EXTRAPOLATED_EXCHANGE_MP2, &
                                                          EXTRAPOLATED_DIRECT_MP2+EXTRAPOLATED_EXCHANGE_MP2

    !DO I = 1, NENCUTS
    !  WRITE(IU6,*) "ENCUT", I, TOT_DI_MP2(I)+TOT_EX_MP2(I)
    !ENDDO
    !WRITE(IU6,*) "converged", EXTRAPOLATED_ENERGY

    CLOSE(DI_MEAN_FILEUNIT)
    CLOSE(EX_MEAN_FILEUNIT)
    CLOSE(DI_VAR_FILEUNIT)
    CLOSE(EX_VAR_FILEUNIT)
    CLOSE(DI_ERR_FILEUNIT)
    CLOSE(EX_ERR_FILEUNIT)

  END SUBROUTINE CLOSE_STATISTICS  



  !****************************** LIN_REG_ENERGIES *****************************
  !
  ! Perform a linear regression of the input data: y = M*x + B, where 
  ! x = ENCUTS^(-3/2) and y = ENERGIES. 
  !
  ! This is used to extrapolate ENCUT -> infinity
  ! Because: x(ENCUT -> infinity) = 0 => y -> b = extrapolated energy
  !
  ! see: http://dx.doi.org/10.1103/PhysRevB.90.075125
  !
  !***********************************************************************************
  SUBROUTINE LIN_REG_ENERGIES(ENERGIES, ENCUTS, EXTRAPOLATED_ENERGY)
    REAL(q) :: ENERGIES(:)         ! calculated energies, each for a given ENCUT
    REAL(q) :: ENCUTS(:)           ! all ENCUTs
    REAL(q) :: EXTRAPOLATED_ENERGY    ! the resulting converged energy
    ! local
    INTEGER :: I
    INTEGER :: N                            ! number of sample points
    REAL(q) :: SUM_X, SUM_X2, SUM_Y, SUM_XY ! sums of the samples
    REAL(q) :: M, B                         ! parameters of y = M*x + B
  
    N = SIZE(ENERGIES)
    SUM_X  = 0
    SUM_X2 = 0
    SUM_Y  = 0
    SUM_XY = 0
  
    ! perform linear regression
    DO I = 1, N
      SUM_X  = SUM_X  + ENCUTS(I)**(-3.0_q/2)
      SUM_X2 = SUM_X2 + ENCUTS(I)**(-3.0_q)    ! ENCUTS^3 corresponds to x^2
      SUM_Y  = SUM_Y  + ENERGIES(I)
      SUM_XY = SUM_XY + ENCUTS(I)**(-3.0_q/2) * ENERGIES(I)
    ENDDO
    M = ( N*SUM_XY - SUM_X*SUM_Y) / ( N*SUM_X2 - SUM_X**2 )
    B = ( SUM_Y - M*SUM_X ) / N
  
    EXTRAPOLATED_ENERGY = B
  
  END SUBROUTINE LIN_REG_ENERGIES



  !***************************** ADD_DI_SAMPLES *******************************
  ! Adds direct samples and calcualtes all statistical quantities
  !************************************************************************
  SUBROUTINE ADD_DI_SAMPLES(DMP2_SAMPLES, NSMPLS_PER_CPU, NCPU)
    COMPLEX(q) :: DMP2_SAMPLES(:,:,:)
    INTEGER :: NSMPLS_PER_CPU
    INTEGER :: NCPU
    ! local
    INTEGER :: I, J
    ! new calcualtion
    REAL(q) :: DI_MEAN_OLD(NENCUTS)
    REAL(q) :: THIS_SAMPLE(NENCUTS)

#ifdef PROFILING
    CALL START_PROFILING('ADD_SAMPLES')
#endif

    DO J = 1, NCPU
      DO I = 1, NSMPLS_PER_CPU
        N_DI_SAMPLES = N_DI_SAMPLES + 1

        THIS_SAMPLE(:) = REAL(DMP2_SAMPLES(I,:,J),q)

        !! dum sample ---------------------------------------
        !WRITE(300,'(F20.10)') THIS_SAMPLE(1)
        !!---------------------------------------------------

        IF(N_DI_SAMPLES < 2) THEN
          DI_MEAN(:) = THIS_SAMPLE(:)
          DI_VAR = 0
        ELSE
          DI_MEAN_OLD(:) = DI_MEAN(:)
          DI_MEAN(:) = DI_MEAN_OLD(:) + (THIS_SAMPLE(:) - DI_MEAN_OLD(:)) / (N_DI_SAMPLES*1.0_q)

          ! Welfords algorithm to calculate variance
          DI_VAR = (N_DI_SAMPLES - 2.0_q) * DI_VAR + (THIS_SAMPLE(1) - DI_MEAN_OLD(1)) * (THIS_SAMPLE(1) - DI_MEAN(1))
          DI_VAR = DI_VAR / (N_DI_SAMPLES - 1.0_q)
        ENDIF
      ENDDO
    ENDDO

    DI_MP2(:, CURRENT_Q, CURRENT_TAU) = DI_MEAN(:)

    DI_ERR = SQRT(DI_VAR / N_DI_SAMPLES)
    DI_ERRS(1, CURRENT_Q, CURRENT_TAU) = DI_ERR

#ifdef PROFILING
    CALL STOP_PROFILING('ADD_SAMPLES')
#endif

  END SUBROUTINE ADD_DI_SAMPLES

  !***************************** ADD_EX_SAMPLE *******************************
  ! Adds exchange samples and calcualtes all statistical quantities
  !************************************************************************
  SUBROUTINE ADD_EX_SAMPLES(XMP2_SAMPLES, NSMPLS_PER_CPU, NCPU)
    COMPLEX(q) :: XMP2_SAMPLES(:,:,:)
    INTEGER :: NSMPLS_PER_CPU
    INTEGER :: NCPU
    ! local
    INTEGER :: I, J
    REAL(q) :: EX_MEAN_OLD(NENCUTS)
    REAL(q) :: THIS_SAMPLE(NENCUTS)

#ifdef PROFILING
    CALL START_PROFILING('ADD_SAMPLES')
#endif

    DO J = 1, NCPU
      DO I = 1, NSMPLS_PER_CPU
        N_EX_SAMPLES = N_EX_SAMPLES + 1

        THIS_SAMPLE(:) = REAL(XMP2_SAMPLES(I,:,J),q)
        
        !! dum sample ---------------------------------------
        !WRITE(400,'(F20.10)') THIS_SAMPLE(1)
        !!---------------------------------------------------

        IF(N_EX_SAMPLES < 2) THEN
          EX_MEAN(:) = THIS_SAMPLE(:)
          EX_VAR = 0
        ELSE
          EX_MEAN_OLD(:) = EX_MEAN(:)
          EX_MEAN(:) = EX_MEAN_OLD(:) + (THIS_SAMPLE(:) - EX_MEAN_OLD(:)) / (N_EX_SAMPLES*1.0_q)

          ! Welfords algorithm to calculate variance
          EX_VAR = (N_EX_SAMPLES - 2.0_q) * EX_VAR + (THIS_SAMPLE(1) - EX_MEAN_OLD(1)) * (THIS_SAMPLE(1) - EX_MEAN(1))
          EX_VAR = EX_VAR / (N_EX_SAMPLES - 1.0_q)
        ENDIF
      ENDDO
    ENDDO

    EX_MP2(:, CURRENT_Q, CURRENT_TAU) = EX_MEAN(:)

    EX_ERR = SQRT(EX_VAR / N_EX_SAMPLES)
    EX_ERRS(1, CURRENT_Q, CURRENT_TAU) = EX_ERR

#ifdef PROFILING
    CALL STOP_PROFILING('ADD_SAMPLES')
#endif

  END SUBROUTINE ADD_EX_SAMPLES



  !***************************** GET_STAT_ERRORS  *******************************
  !
  ! returns the statistical error of the direct and exchange term
  ! error = sqrt(var / #samples)
  ! 
  !******************************************************************************
  SUBROUTINE GET_STAT_ERRORS(DI_ERROR, EX_ERROR)
    REAL(q) :: DI_ERROR
    REAL(q) :: EX_ERROR

    DI_ERROR = DI_ERR
    EX_ERROR = EX_ERR

  END SUBROUTINE GET_STAT_ERRORS  



  !***************************** DUMP_DI_STATISTICS  *******************************
  ! writes the direct statistics to the output files
  !******************************************************************************
  SUBROUTINE DUMP_DI_STATISTICS(IU0)  
    INTEGER :: IU0
    IF (IU0 >= 0) THEN
      WRITE(DI_MEAN_FILEUNIT, '(I12,2F20.10)') N_DI_SAMPLES, DI_MP2(1, CURRENT_Q, CURRENT_TAU)
      WRITE(DI_ERR_FILEUNIT, '(I12,F20.10)') N_DI_SAMPLES, DI_ERR
      WRITE(DI_VAR_FILEUNIT, '(I12,2E20.10)') N_DI_SAMPLES, DI_VAR
    ENDIF
  END SUBROUTINE DUMP_DI_STATISTICS  
  !***************************** DUMP_EX_STATISTICS  *******************************
  ! writes the direct statistics to the output files
  !******************************************************************************
  SUBROUTINE DUMP_EX_STATISTICS(IU0)  
    INTEGER :: IU0
    IF (IU0 >= 0) THEN
      WRITE(EX_MEAN_FILEUNIT, '(I12,F20.10)') N_EX_SAMPLES, EX_MP2(1, CURRENT_Q, CURRENT_TAU)
      WRITE(EX_ERR_FILEUNIT, '(I12,F20.10)') N_EX_SAMPLES, EX_ERR
      WRITE(EX_VAR_FILEUNIT, '(I12,E20.10)') N_EX_SAMPLES, EX_VAR
    ENDIF
  END SUBROUTINE DUMP_EX_STATISTICS  

END MODULE rnd_stat






!***********************************************************************
!
! this module implements the stochastic laplace transformed MP2 
! based on an approach by Baer and Neuhauser 
! (dx.doi.org/10.1021/ct300946j , J. Chem. Theory Comput. 2013, 9, 24−27)
!
!***********************************************************************
MODULE stoch_lt_mp2
  USE prec
  USE fock
  USE chi_base
  USE lattice
  USE full_kpoints
  USE paw
  USE sym_prec
#ifdef PROFILING
  USE profiling
#endif
  IMPLICIT none

  ! different ENCUTs and ENCUTSOFTs for extrapolation
  INTEGER, PARAMETER, PRIVATE :: N_MP2_ENCUTS = 8             ! number of different ENCUTs
  REAL(q), PARAMETER, PRIVATE :: ENCUT_DIVISOR = 1.05_q       ! ENCUTS(n+1) = ENCUTS(n) / ENCUT_DIVISOR
  REAL(q), ALLOCATABLE, PRIVATE :: ENCUTS(:)                  ! different ENCUTs
  REAL(q), ALLOCATABLE, PRIVATE :: ENCUTSOFTS(:)              ! different ENCUTSOFTs
  REAL(q), PRIVATE :: ENCUT_LOW                               ! lower cutoff for the low energy MP2 contibution

CONTAINS
    
  !****************************** RND_ORBITALS_MP2 ******************************
  !
  ! This is the main routine that calculates the MP2 energy using stochastic
  ! orbitals as described in https://arxiv.org/abs/1711.05709
  !
  ! Additionally it can used to calculate only the high energy contributions to
  ! the MP2 correlation energy by setting NB_LOW > 0
  !
  !*************************************************************************************
  SUBROUTINE RND_ORBITALS_MP2(LATT_CUR, W, LMDIM, CQIJ, NONLR_S, NONL_S, P, TAU, TAUWEIGHT, WGW, &
                              ENCUT, ENCUTSOFT, IU0, IU5, IU6)
    USE ini
    USE reader_tags
    USE nonl_high
    USE wave_high
    USE mpimy
    USE mgrid
    USE constant
    USE pseudo
    USE kpoints_change, ONLY : KPOINTS_ORIG
    USE rnd_stat
    IMPLICIT NONE

    TYPE (latt) :: LATT_CUR
    TYPE (wavespin) ::  W
    INTEGER :: LMDIM
    OVERLAP CQIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ) ! not used
    TYPE (nonlr_struct) :: NONLR_S
    TYPE (nonl_struct) :: NONL_S
    TYPE (potcar) :: P(NONLR_S%NTYP)
    REAL(q) :: TAU(:)                    ! tau-points for the evaluation of the tau-integral
    REAL(q) :: TAUWEIGHT(:)              ! tau-weights for Gauss integration
    TYPE (wavedes), POINTER :: WGW       
    REAL(q) :: ENCUT, ENCUTSOFT          ! cutoff for outer sums over G (actually only the potential is set to 0)
    INTEGER :: IU0, IU5, IU6
    ! local
    INTEGER :: I, IERR
    INTEGER :: NQ, NK, NPL
    INTEGER :: NTAU
    INTEGER :: LAST_FILLED
    INTEGER :: NB_LOW                                ! number of bands for the lower band cutoff (NB_LOW < NBANDS) 
    INTEGER :: NCYCLES                               ! number of outer loop cycles
    INTEGER :: ICYCLE
    INTEGER :: NTHETA                                ! number of random basis vectors
    INTEGER :: NTHETA_EXTENDED                       ! NTHETA_EXTENDED = NTHETA * W%WDES%COMM_INTER 
    INTEGER :: NCORREL_THETA                         ! number of THETA for correlated sampling (should be divisior of NTHETA)
    INTEGER :: NCORREL_THETA_BLOCKS                  ! = NTHETA / NCORREL_THETA
    INTEGER :: ICORREL_THETA_BLOCKS                  ! iteration variable for NCORREL_THETA_BLOCKS
    INTEGER :: CORREL_THETA_START                    ! start and...
    INTEGER :: CORREL_THETA_STOP                     ! ...stop of current NTHETA block
    INTEGER :: REACHED_EX_ACCURACY, REACHED_DI_ACCURACY
    INTEGER :: ROOT_NODE                             ! MPI root node
    TYPE (wavedes)  :: WDESOCC                       ! descriptor to store occupied orbitals
    TYPE (wavedes)  :: WDESRNDOCC1, WDESRNDOCC2      ! descriptor to store random orbitals
    TYPE (wavedes)  :: WDESRNDVIRT1, WDESRNDVIRT2    ! descriptor to store random orbitals
    TYPE (wavedes)  :: WDESRNDVIRT1_LOW, WDESRNDVIRT2_LOW    ! descriptor to store random orbitals (for lower NBANDS)
    TYPE (wavedes1) :: WDESK 
    TYPE (wavespin) :: WOCC                          ! occupied orbitals |i>
    TYPE (wavespin) :: WRNDOCC1, WRNDOCC2            ! stochastic occupied orbitals |kappa> and |lambda> 
    TYPE (wavespin) :: WRNDVIRT1, WRNDVIRT2          ! stochastic virtual orbitals |alpha> and |beta>
    TYPE (wavespin) :: WRNDVIRT1_LOW, WRNDVIRT2_LOW  ! stochastic virtual orbitals with lower NBANDS (NB_LOW)
    TYPE (wavedes1)    WGWQ                          ! wave descriptor for outer q-point loop
    !COMPLEX(q), ALLOCATABLE :: OVERLAPS1(:,:)       ! stores <kappa_theta  | G | alpha_theta>
    !COMPLEX(q), ALLOCATABLE :: OVERLAPS2(:,:)       ! stores <lambda_theta | G | beta_theta>
    !COMPLEX(q), ALLOCATABLE :: OVERLAPS3(:,:)       ! stores <alpha_theta | G | lambda_theta>
    !COMPLEX(q), ALLOCATABLE :: OVERLAPS4(:,:)       ! stores <beta_theta  | G | kappa_theta>
    COMPLEX(q), ALLOCATABLE :: DMP2_SAMPLES(:,:)     ! buffer to store sample of direct mp2 energy
    COMPLEX(q), ALLOCATABLE :: XMP2_SAMPLES(:,:)     ! buffer to store sample of exchange mp2 energy
    COMPLEX(q), ALLOCATABLE :: DMP2_SAMPLES_GATH(:,:,:)   ! buffer to gather samples of direct mp2 energy from all cores
    COMPLEX(q), ALLOCATABLE :: XMP2_SAMPLES_GATH(:,:,:)   ! buffer to gather sample of exchange mp2 energy from all cores
    COMPLEX(q) :: EDRCT_TOT, EEXCHG_TOT              ! total direct and exchange MP2 energy
    REAL(q) :: DI_ERR, EX_ERR, ESTOP
    LOGICAL :: CHANGED_CORREL_THETA
    REAL(q) :: ENCUTGW_LOW                           ! lower than ENCUTGW
    INTEGER :: SYS_CLOCK
    LOGICAL :: LOPEN_INCAR
    REAL(q) :: EFERMI(W%WDES%ISPIN)

#ifdef PROFILING
    CALL START_PROFILING('RND_ORBITALS_MP2')
#endif

    EFERMI = W%EFERMI
    ! read the flags from INCAR  -----------------------------------------------------
    CALL OPEN_INCAR_IF_FOUND(IU5, LOPEN_INCAR)

    NTHETA = -1 
    CALL PROCESS_INCAR(LOPEN_INCAR, IU0, IU5, 'NSTORB', NTHETA, IERR, WRITEXMLINCAR)
    IF(NTHETA < 0) THEN
      CALL vtutor%error("ERROR: number of stochastic orbitals (NSTORB) not set or illegal value in INCAR!")
    ENDIF
    
    NCORREL_THETA = NTHETA ! standard value
    CALL PROCESS_INCAR(LOPEN_INCAR, IU0, IU5, 'NSCORREL', NCORREL_THETA, IERR, WRITEXMLINCAR)
    IF(NCORREL_THETA < 0 .OR. MOD(NTHETA,NCORREL_THETA) > 0) THEN
      CALL vtutor%error("ERROR: number of stochastic orbitals to correlate (NSCORREL) not set or &
         &illegal value in INCAR!")
    ENDIF
    
    ESTOP = 0.05 ! standard value
    CALL PROCESS_INCAR(LOPEN_INCAR, IU0, IU5, 'ESTOP', ESTOP, IERR, WRITEXMLINCAR)
    IF(ESTOP < 0.0) THEN
      CALL vtutor%error("ERROR: illegal value for ESTOP in INCAR!")
    ENDIF

    NCYCLES=-1
    CALL PROCESS_INCAR(LOPEN_INCAR, IU0, IU5, 'NCYCLES', NCYCLES, IERR, WRITEXMLINCAR)

    ENCUT_LOW = 0 ! standard value
    CALL PROCESS_INCAR(LOPEN_INCAR, IU0, IU5, 'ENCUTGWLOW', ENCUT_LOW, IERR, WRITEXMLINCAR)
    IF(ENCUT_LOW < 0 .OR. ENCUT_LOW > ENCUT) THEN
      CALL vtutor%error("ERROR: illegal value for ENCUTGWLOW in INCAR!")
    ENDIF

    NB_LOW = W%WDES%NB_TOT * (ENCUT_LOW/ENCUT)**(1.5) ! standard value
    CALL PROCESS_INCAR(LOPEN_INCAR, IU0, IU5, 'NBANDSLOW', NB_LOW, IERR, WRITEXMLINCAR)
    IF(NB_LOW < 0 .OR. NB_LOW > W%WDES%NB_TOT) THEN
      CALL vtutor%error("ERROR: illegal value for NBANDSLOW in INCAR!")
    ENDIF

    CALL CLOSE_INCAR_IF_FOUND(IU5)

    IF(NCYCLES < 1 .AND. ESTOP < 0.0000001) THEN
      WRITE(IU0,'(A)') " WARNING: no meaningfull stop criterion found for the stochastic LTMP2 calculation! Will probably run forever..."
    ENDIF
    !--------------------------------------------------------------------------------------

    !status output
    IF (IU0>0) THEN
       WRITE(IU0,'(A)') " " 
       WRITE(IU0,'(A)') " Starting the stochastic LTMP2 calculation" 
       WRITE(IU0,'(A)') " -------------------------------" 
    ENDIF

    IF(IU6>0) THEN
      WRITE(IU6,'(A)') ""
      WRITE(IU6,'(A)') " Settings for stochastic orbitals:"
      WRITE(IU6,'(A)') " ================================="
      WRITE(IU6,'(A18,I8)') " NTHETA = ", NTHETA
      WRITE(IU6,'(A18,I8)') " NCORREL_THETA = ", NCORREL_THETA
      WRITE(IU6,'(A18,E12.3)') " ESTOP = ", ESTOP
      WRITE(IU6,'(A18,I8)') " NCYCLES = ", NCYCLES
      WRITE(IU6,'(A18)') ""
      IF(NB_LOW > 0) THEN
        WRITE(IU6,'(A)') " NOTE: ONLY THE HIGH ENERGY MP2 CONTRIBUTIONS ARE CALCULATED"
        WRITE(IU6,'(A18,I8)') " NB_LOW = ", NB_LOW
      WRITE(IU6,'(A18,E12.3)') " ENCUT_LOW = ", ENCUT_LOW
      ENDIF
    ENDIF
    IF(IU0>0 .AND. NB_LOW > 0) THEN
      WRITE(IU0,'(A)') " NOTE: ONLY THE HIGH ENERGY MP2 CONTRIBUTIONS ARE CALCULATED"
    ENDIF

    NCORREL_THETA_BLOCKS = NTHETA / NCORREL_THETA

    ! set root node as the output node
    ROOT_NODE = 0
    IF(IU6>0) ROOT_NODE = W%WDES%COMM%NODE_ME
    CALLMPI(M_sum_i(W%WDES%COMM, ROOT_NODE, 1))

    ! allocate memory for the different ENCUTs
    ALLOCATE(ENCUTS(N_MP2_ENCUTS))
    ALLOCATE(ENCUTSOFTS(N_MP2_ENCUTS))
    ! calculate different ENCUTs
    DO I = 1, N_MP2_ENCUTS
      ENCUTS(I) = ENCUT / ( ENCUT_DIVISOR**(I-1) )
      ENCUTSOFTS(I) = ENCUTSOFT / ( ENCUT_DIVISOR**(I-1) )
    ENDDO

    ! first calculate the last filled and first empty orbital indices among all k-points and spins to make live simple
    LAST_FILLED = 0
    DO I = 1, W%WDES%ISPIN
      DO NK = 1, W%WDES%NKPTS
        LAST_FILLED = MAX(LAST_FILLED_XI_NOMOD(W, NK, I), LAST_FILLED)
      ENDDO
    ENDDO

    WDESOCC = W%WDES
    WDESRNDOCC1 = W%WDES
    WDESRNDOCC2 = W%WDES
    WDESRNDVIRT1 = W%WDES
    WDESRNDVIRT2 = W%WDES

    ! round occupied bands to next larger value dividable by NB_PAR
    WDESOCC%NBANDS = (LAST_FILLED + WDESOCC%NB_PAR - 1) / WDESOCC%NB_PAR
    WDESOCC%NB_TOT = WDESOCC%NBANDS * WDESOCC%NB_PAR

    ! since the orbitals can be distributed over several mpi ranks, we need to blow up NTHETA (then 
    ! called NTHETA_EXTENDED) in order to have NTHETA stochastic orbitals at each mpi rank
    NTHETA_EXTENDED = NTHETA * W%WDES%COMM_INTER%NCPU

    ! set global number of random bands to NTHETA_EXTENDED
    !WDESRNDOCC1%NBANDS  = (NTHETA_EXTENDED + WDESRNDOCC1%NB_PAR - 1) / WDESRNDOCC1%NB_PAR
    WDESRNDOCC1%NBANDS  = NTHETA
    WDESRNDOCC1%NB_TOT  = NTHETA_EXTENDED
    WDESRNDOCC2%NBANDS  = NTHETA
    WDESRNDOCC2%NB_TOT  = NTHETA_EXTENDED
    WDESRNDVIRT1%NBANDS = NTHETA
    WDESRNDVIRT1%NB_TOT = NTHETA_EXTENDED
    WDESRNDVIRT2%NBANDS = NTHETA
    WDESRNDVIRT2%NB_TOT = NTHETA_EXTENDED

    ! now create orbital array for the occupied and stochastic orbitals
    CALL ALLOCW(WDESOCC, WOCC)           ! occupied orbitals
    CALL ALLOCW(WDESRNDOCC1, WRNDOCC1)   ! random occupied orbitals
    CALL ALLOCW(WDESRNDOCC2, WRNDOCC2)   ! random occupied orbitals
    CALL ALLOCW(WDESRNDVIRT1, WRNDVIRT1) ! random virtual orbitals
    CALL ALLOCW(WDESRNDVIRT2, WRNDVIRT2) ! random virtual orbitals

    ! do the same for the low energy contribution (if desired)
    IF(NB_LOW > 0) THEN
      WDESRNDVIRT1_LOW = W%WDES
      WDESRNDVIRT2_LOW = W%WDES
      WDESRNDVIRT1_LOW%NBANDS = NTHETA
      WDESRNDVIRT1_LOW%NB_TOT = NTHETA_EXTENDED
      WDESRNDVIRT2_LOW%NBANDS = NTHETA
      WDESRNDVIRT2_LOW%NB_TOT = NTHETA_EXTENDED
      CALL ALLOCW(WDESRNDVIRT1_LOW, WRNDVIRT1_LOW) ! random virtual orbitals
      CALL ALLOCW(WDESRNDVIRT2_LOW, WRNDVIRT2_LOW) ! random virtual orbitals
    ENDIF

    ! copy original orbitals into the occupied wave function array
    WOCC%CW(:,1:WOCC%WDES%NBANDS,:,:)    = W%CW(:,1:WOCC%WDES%NBANDS,:,:)
    WOCC%CPROJ(:,1:WOCC%WDES%NBANDS,:,:) = W%CPROJ(:,1:WOCC%WDES%NBANDS,:,:)
    WOCC%CELTOT(1:WOCC%WDES%NB_TOT,:,:)  = W%CELTOT(1:WOCC%WDES%NB_TOT,:,:)
    WOCC%FERTOT(1:WOCC%WDES%NB_TOT,:,:)  = W%FERTOT(1:WOCC%WDES%NB_TOT,:,:)
    WOCC%AUXTOT(1:WOCC%WDES%NB_TOT,:,:)  = W%AUXTOT(1:WOCC%WDES%NB_TOT,:,:)

    ! redistribute original orbitals over plane wave coefficients 
    ! necessary for CALC_RND_ORBTIALS: there we need all orbitals at each core 
    DO I = 1, W%WDES%ISPIN
      DO NK = 1, W%WDES%NKPTS
        ! occupied
        CALL SETWDES(WOCC%WDES, WDESK, NK)            
        CALL REDIS_PW(WDESK, WOCC%WDES%NBANDS, WOCC%CW(1, 1, NK, I))
        CALL REDIS_PROJ(WDESK, WOCC%WDES%NBANDS, WOCC%CPROJ(1, 1, NK, I))
        ! unoccupied
        CALL SETWDES(W%WDES, WDESK, NK)            
        CALL REDIS_PW(WDESK, W%WDES%NBANDS, W%CW(1, 1, NK, I))
        CALL REDIS_PROJ(WDESK, W%WDES%NBANDS, W%CPROJ(1, 1, NK, I))
      ENDDO
    ENDDO 

    ! init the random number generator
    ! it is EXTREMELY IMPORTANT that 
    !   (1) SYS_CLOCK is the same for all nodes and 
    !   (2) the second argument of INIT_RND_GENERATOR is the same for all nodes in COMM_KIN
    SYS_CLOCK = 0
    IF(W%WDES%COMM%NODE_ME == 1) CALL SYSTEM_CLOCK(COUNT=SYS_CLOCK)
    CALLMPI(M_sum_i(W%WDES%COMM, SYS_CLOCK, 1))
    CALL INIT_RND_GENERATOR(SYS_CLOCK, W%WDES%COMM_KINTER%NODE_ME)

    EDRCT_TOT = (0.0, 0.0)
    EEXCHG_TOT = (0.0, 0.0)

    ! allocate memory for sample buffer
    ALLOCATE(DMP2_SAMPLES(NCORREL_THETA**2, N_MP2_ENCUTS), &
             XMP2_SAMPLES(NCORREL_THETA**2, N_MP2_ENCUTS))
    CALL REGISTER_ALLOCATE(2*16._q*SIZE(DMP2_SAMPLES), "samples")

    ! allocate memory for sample gather buffer
    IF(IU6>0) THEN
      ALLOCATE(DMP2_SAMPLES_GATH(NCORREL_THETA**2, N_MP2_ENCUTS, W%WDES%COMM%NCPU), &
               XMP2_SAMPLES_GATH(NCORREL_THETA**2, N_MP2_ENCUTS, W%WDES%COMM%NCPU))
      CALL REGISTER_ALLOCATE(2*16._q*SIZE(DMP2_SAMPLES_GATH), "smpl_gath")
    ENDIF

    !OPEN(UNIT = 549, FILE = 'rnd', STATUS = 'old', ACTION = 'read')

    IF(IU6>0) CALL INIT_STATISTICS(SIZE(TAU), KPOINTS_ORIG%NKPTS, N_MP2_ENCUTS)

    CHANGED_CORREL_THETA = .FALSE.

    ! outer q-point loop (loop over all q-points in the ireducible wedge) (UNUSED)
    qpt: & 
    DO NQ = 1, KPOINTS_ORIG%NKPTS
      CALL SETWDES(WGW, WGWQ, NQ)

      NPL = WGWQ%NPL

      ! print warning if NTHETA is < 2 (because this is unefficient)
      IF(IU0>0 .AND. NTHETA < 2) THEN
        WRITE(*,*) " " 
        WRITE(*,*) " !! WARNING !! NTHETA is too small:" 
        WRITE(*,*) "  -> BLAS routines are not efficient in this case." 
        WRITE(*,*) " " 
      ENDIF

      tau_loop: &
      DO NTAU= 1, SIZE(TAU)

        IF(IU0>0) THEN
          WRITE(*,*) "  Working on NTAU", NTAU 
        ENDIF

        ! for the last few tau points we should decrease NCORREL_THETA by a factor of 0.25
        ! since the contribution to the MP2 energy converges very fast for large tau
        IF(1.0*NTAU/SIZE(TAU) >= 0.75 .AND. NCORREL_THETA >= 15 .AND. ESTOP > 0.0 .AND. .NOT. CHANGED_CORREL_THETA) THEN
          CHANGED_CORREL_THETA = .TRUE.
          NCORREL_THETA = NCORREL_THETA / 5
          NCORREL_THETA_BLOCKS = NTHETA / NCORREL_THETA
          DEALLOCATE(DMP2_SAMPLES, XMP2_SAMPLES)
          ALLOCATE(DMP2_SAMPLES(NCORREL_THETA**2, N_MP2_ENCUTS), XMP2_SAMPLES(NCORREL_THETA**2, N_MP2_ENCUTS))
          IF(IU6>0) THEN
            DEALLOCATE(DMP2_SAMPLES_GATH, XMP2_SAMPLES_GATH)
            ALLOCATE(DMP2_SAMPLES_GATH(NCORREL_THETA**2, N_MP2_ENCUTS, W%WDES%COMM%NCPU), &
                     XMP2_SAMPLES_GATH(NCORREL_THETA**2, N_MP2_ENCUTS, W%WDES%COMM%NCPU))
            WRITE(IU6,*) " "
            WRITE(IU6,'(A,I8)') " changed NCORREL_THETA to ", NCORREL_THETA
            WRITE(IU6,*) " "
          ENDIF
        ENDIF

        REACHED_DI_ACCURACY = 0
        REACHED_EX_ACCURACY = 0

        IF(IU6>0) CALL FLUSH_STATISTICS(NTAU, NQ)

        ICYCLE = 0
        DO ! loop until ICYCLE > NCYCLES or until the desired statistical error is reached
          ICYCLE = ICYCLE + 1
          IF(ICYCLE > NCYCLES .AND. NCYCLES > 0) EXIT

          ! calculate the stochastic occupied and virtual orbitals and store in WRNDOCC and WRNDVIRT
          CALL CALC_RND_ORBITALS(W, WOCC, NTHETA, NTHETA_EXTENDED, NQ, TAU(NTAU), EFERMI,  WRNDOCC1, WRNDOCC2, &
                                 WRNDVIRT1, WRNDVIRT2, NB_LOW, WRNDVIRT1_LOW, WRNDVIRT2_LOW)

          ! loop through the NTHETA_CORREL BLOCK
          DO ICORREL_THETA_BLOCKS = 1, NCORREL_THETA_BLOCKS

            ! calc first and last NTHETA index for thiy cycle
            CORREL_THETA_START = NCORREL_THETA * (ICORREL_THETA_BLOCKS-1) + 1
            CORREL_THETA_STOP  = CORREL_THETA_START + NCORREL_THETA - 1

            CALL CALC_DIRECT_AND_EXCHANGE_MP2(WRNDOCC1, WRNDOCC2, WRNDVIRT1, WRNDVIRT2, LATT_CUR, WGWQ, NTHETA, &
                                              NCORREL_THETA, CORREL_THETA_START, CORREL_THETA_STOP, W%WDES%RSPIN, &
                                              REACHED_DI_ACCURACY, REACHED_EX_ACCURACY, DMP2_SAMPLES, XMP2_SAMPLES, &
                                              NB_LOW, WRNDVIRT1_LOW, WRNDVIRT2_LOW)

            ! dump memory usage
            IF(NTAU == 1 .AND. ICYCLE == 1 .AND. ICORREL_THETA_BLOCKS == 1) CALL DUMP_ALLOCATE(IU6)

            DMP2_SAMPLES(:,:) = DMP2_SAMPLES(:,:) * TAUWEIGHT(NTAU) * KPOINTS_ORIG%WTKPT(NQ) 
            XMP2_SAMPLES(:,:) = XMP2_SAMPLES(:,:) * TAUWEIGHT(NTAU) * KPOINTS_ORIG%WTKPT(NQ) 

            CALLMPI(MPI_GATHER(DMP2_SAMPLES, SIZE(DMP2_SAMPLES), MPI_double_complex, DMP2_SAMPLES_GATH, SIZE(DMP2_SAMPLES), MPI_double_complex, ROOT_NODE-1, W%WDES%COMM%MPI_COMM, IERR))
            IF(IERR > 0 .AND. IU0>0) THEN
              CALL vtutor%error("ERROR in RND_ORBITALS_MP2: MPI_GATHER failed")
            ENDIF

            CALLMPI(MPI_GATHER(XMP2_SAMPLES, SIZE(XMP2_SAMPLES), MPI_double_complex, XMP2_SAMPLES_GATH, SIZE(XMP2_SAMPLES), MPI_double_complex, ROOT_NODE-1, W%WDES%COMM%MPI_COMM, IERR))
            IF(IERR > 0 .AND. IU0>0) THEN
              CALL vtutor%error("ERROR in RND_ORBITALS_MP2: MPI_GATHER failed")
            ENDIF

            IF(IU6>0) THEN
              IF(REACHED_DI_ACCURACY == 0) THEN
                CALL ADD_DI_SAMPLES(DMP2_SAMPLES_GATH, NCORREL_THETA**2, W%WDES%COMM%NCPU)
                CALL DUMP_DI_STATISTICS(IU0)
              ENDIF
              IF(REACHED_EX_ACCURACY == 0) THEN
                CALL ADD_EX_SAMPLES(XMP2_SAMPLES_GATH, NCORREL_THETA**2, W%WDES%COMM%NCPU)
                CALL DUMP_EX_STATISTICS(IU0)
              ENDIF

              ! check if accuracy is reached
              CALL GET_STAT_ERRORS(DI_ERR, EX_ERR)
              IF(DI_ERR < ESTOP .AND. 1_8 * ICYCLE*CORREL_THETA_STOP*W%WDES%COMM%NCPU > 1000) REACHED_DI_ACCURACY = 1
              IF(EX_ERR < ESTOP .AND. 1_8 * ICYCLE*CORREL_THETA_STOP*W%WDES%COMM%NCPU > 1000) REACHED_EX_ACCURACY = 1
            ELSE
              ! all other mpi ranks should have zero otherwise M_sum_i produces nonsense
              REACHED_DI_ACCURACY = 0
              REACHED_EX_ACCURACY = 0
            ENDIF

            CALLMPI(M_sum_i(W%WDES%COMM, REACHED_DI_ACCURACY, 1))
            CALLMPI(M_sum_i(W%WDES%COMM, REACHED_EX_ACCURACY, 1))
            IF(REACHED_DI_ACCURACY == 1 .AND. REACHED_EX_ACCURACY == 1) EXIT
          
          ENDDO

          IF(REACHED_DI_ACCURACY == 1 .AND. REACHED_EX_ACCURACY == 1) EXIT

        ENDDO

        IF(IU6>0) CALL STATISTICS_TO_OUTCAR(IU6)

      ENDDO tau_loop

      !DEALLOCATE(OVERLAPS1, OVERLAPS2, OVERLAPS3, OVERLAPS4)

    ENDDO qpt

#ifdef PROFILING
    CALL STOP_PROFILING('RND_ORBITALS_MP2')
#endif

    ! close the files for the statistic output
    IF(IU6>0) THEN
      CALL CLOSE_STATISTICS(IU0, IU6, ENCUTS, ENCUTSOFTS) ! write results to terminal and OUTCAR (also exrapolates MP2 result)
      CALL DEALLOCATE_STATISTICS()
    ENDIF

    IF(IU6>0) DEALLOCATE(DMP2_SAMPLES_GATH, XMP2_SAMPLES_GATH)
    DEALLOCATE(DMP2_SAMPLES, XMP2_SAMPLES)

    CALL DEALLOCW(WOCC)
    CALL DEALLOCW(WRNDOCC1)
    CALL DEALLOCW(WRNDOCC2)
    CALL DEALLOCW(WRNDVIRT1)
    CALL DEALLOCW(WRNDVIRT2)

    IF(NB_LOW > 0) THEN
      CALL DEALLOCW(WRNDVIRT1_LOW)
      CALL DEALLOCW(WRNDVIRT2_LOW)
    ENDIF

    DEALLOCATE(ENCUTS)
    DEALLOCATE(ENCUTSOFTS)

    CALL XML_CLOSE_TAG("calculation") 
    CALL DUMP_FINAL_TIMING(IU6)
    CALL STOP_XML
#ifdef PROFILING
     ! started in INIT_PROFILING
     PROFILING_STOP('total_time')
     IF (IU0>=0) THEN
        CALL PRINT_FULL_PROFILE(IU6,PFLAT=.FALSE.)
     ENDIF
#endif

    CALL vtutor%stopCode()
  END SUBROUTINE RND_ORBITALS_MP2
  


  !********************** INIT_RND_GENERATOR  ******************************
  ! init the random number generator
  !*************************************************************************
  SUBROUTINE INIT_RND_GENERATOR(FIXVAL, VARVAL)
    INTEGER :: FIXVAL         ! a number to initialize the random generator which all nodes share
    INTEGER :: VARVAL         ! a number to individually initialize the random generator for each node
    !local
    INTEGER :: I
    INTEGER :: RND_SEED_SIZE    ! number of integers for the random seed
    INTEGER, DIMENSION(:), ALLOCATABLE :: RND_SEED    ! random seed

    ! initialze random generator with system time
    CALL RANDOM_SEED(SIZE = RND_SEED_SIZE)
    ALLOCATE(RND_SEED(RND_SEED_SIZE))
    DO I = 1, RND_SEED_SIZE
      RND_SEED(I) = FIXVAL + 37 * VARVAL + 997 * I   ! choose the seed values
    END DO
    CALL RANDOM_SEED(PUT = RND_SEED)
    DEALLOCATE(RND_SEED)
  ENDSUBROUTINE INIT_RND_GENERATOR



  !********************** CALC_RND_ORBITALS  ******************************
  ! 
  ! Calculates the stochastic orbitals 
  ! |kappa_theta>  = WRNDOCC1  = sum_i^occ  r1(i,theta) e^(+epsilon_i tau/2) |i>
  ! |lambda_theta> = WRNDOCC2  = sum_i^occ  r2(i,theta) e^(+epsilon_i tau/2) |i>
  ! |alpha_theta>  = WRNDVIRT1 = sum_a^virt r3(a,theta) e^(-epsilon_a tau/2) |a>
  ! |beta_theta>   = WRNDVIRT2 = sum_a^virt r4(a,theta) e^(-epsilon_a tau/2) |a>
  !
  ! and if passed then also the stochastic virtual orbitals for a smaller number of
  ! bands are calculated
  ! |alpha_low_theta>  = WRNDVIRT1 = sum_a^lowvirt r3(a,theta) e^(-epsilon_a tau/2) |a>
  ! |beta_low_theta>   = WRNDVIRT2 = sum_a^lowvirt r4(a,theta) e^(-epsilon_a tau/2) |a>
  ! where NBLOW = lowvirt
  ! 
  !************************************************************************
  SUBROUTINE CALC_RND_ORBITALS(W, WOCC, NTHETA, NTHETA_EXTENDED, NQ, TAU, EFERMI, WRNDOCC1, WRNDOCC2, WRNDVIRT1, WRNDVIRT2, & 
                               NBLOW, WRNDVIRT1_LOW, WRNDVIRT2_LOW)
    USE ini
    TYPE (wavespin) ::  W
    TYPE (wavespin) ::  WOCC
    INTEGER :: NTHETA
    INTEGER :: NTHETA_EXTENDED       ! NTHETA_EXTENDED = NTHETA * W%WDES%COMM_INTER%NCPU
    INTEGER :: NQ
    REAL(q) :: TAU
    REAL(q) :: EFERMI(W%WDES%ISPIN)
    TYPE (wavespin) ::  WRNDOCC1
    TYPE (wavespin) ::  WRNDOCC2
    TYPE (wavespin) ::  WRNDVIRT1
    TYPE (wavespin) ::  WRNDVIRT2
    INTEGER, OPTIONAL ::  NBLOW
    TYPE (wavespin), OPTIONAL ::  WRNDVIRT1_LOW
    TYPE (wavespin), OPTIONAL ::  WRNDVIRT2_LOW
    ! local
    INTEGER NB_ALL, NB_OCC
    INTEGER ::  NB_LOW
    INTEGER ITHETA, I, A
    INTEGER :: ISP, K1, K2
    GDEF, ALLOCATABLE :: RND_MAT_OCC1(:,:), RND_MAT_VIRT1(:,:) ! random matrices
    GDEF, ALLOCATABLE :: RND_MAT_OCC2(:,:), RND_MAT_VIRT2(:,:) ! random matrices
#ifndef gammareal
    REAL(q), ALLOCATABLE :: RND_MAT_AUX1(:,:), RND_MAT_AUX2(:,:) ! auxiliary matrices to generate complex random arrays
#endif
    REAL(q) :: FAC
    TYPE (wavefuna) :: WA, WI                                     ! wavefun for a virtual/unoccupied orbital
    TYPE (wavefuna) :: WRNDA1, WRNDA2, WRNDI1, WRNDI2             ! wavefun for a randm virtual/occupied orbital
    TYPE (wavefuna) :: WRNDA1_LOW, WRNDA2_LOW                     ! wavefun for a randm virtual orbital for the lower NBANDS
    TYPE (wavedes1) :: WDESK1, WDESK2

#ifdef PROFILING
    CALL START_PROFILING('CALC_RND_ORBITALS')
#endif
  
    NB_LOW = -1 
    IF(PRESENT(NBLOW)) NB_LOW = NBLOW 

    ISP = 1 ! (SPIN IGNORED)
    K1 = 1 ! GAMMA ONLY
    K2 = 1 ! GAMMA ONLY
        
    CALL SETWDES(WOCC%WDES, WDESK1, K1) ! GAMMA ONLY
    CALL SETWDES(W%WDES, WDESK2, K2) ! GAMMA ONLY

    NB_ALL = W%WDES%NB_TOT
    NB_OCC = WOCC%WDES%NB_TOT

    ! allocate memory for the random matrices
    ALLOCATE( RND_MAT_OCC1(NB_OCC, NTHETA_EXTENDED), RND_MAT_VIRT1(NB_ALL, NTHETA_EXTENDED), &
              RND_MAT_OCC2(NB_OCC, NTHETA_EXTENDED), RND_MAT_VIRT2(NB_ALL, NTHETA_EXTENDED) )
    CALL REGISTER_ALLOCATE(16.0_q / m_ (2*SIZE(RND_MAT_OCC1)+2*SIZE(RND_MAT_VIRT1)), "rand_mat")

    RND_MAT_OCC1 = 0
    RND_MAT_OCC2 = 0
    RND_MAT_VIRT1 = 0
    RND_MAT_VIRT2 = 0

    
    ! generate random numbers:
    ! inside COMM_KIN the random numbers are equal due to the equal seed. this
    ! is necessary since the random matrix has to be equal inside COMM_KIN and
    ! the redundant generation of random numbers avoids communication.
#ifdef PROFILING
    CALL START_PROFILING('RANDOM_NUMBER')
#endif
#ifdef gammareal
    CALL RANDOM_NUMBER(RND_MAT_OCC1)
    CALL RANDOM_NUMBER(RND_MAT_OCC2)
    CALL RANDOM_NUMBER(RND_MAT_VIRT1)
    CALL RANDOM_NUMBER(RND_MAT_VIRT2)
#else
    ALLOCATE(RND_MAT_AUX1(NB_OCC, NTHETA_EXTENDED), RND_MAT_AUX2(NB_OCC, NTHETA_EXTENDED))
    CALL RANDOM_NUMBER(RND_MAT_AUX1)
    CALL RANDOM_NUMBER(RND_MAT_AUX2)
    RND_MAT_OCC1 = CMPLX(RND_MAT_AUX1, RND_MAT_AUX2,q)
    CALL RANDOM_NUMBER(RND_MAT_AUX1)
    CALL RANDOM_NUMBER(RND_MAT_AUX2)
    RND_MAT_OCC2 = CMPLX(RND_MAT_AUX1, RND_MAT_AUX2,q)
    DEALLOCATE(RND_MAT_AUX1, RND_MAT_AUX2)
    ALLOCATE(RND_MAT_AUX1(NB_ALL, NTHETA_EXTENDED), RND_MAT_AUX2(NB_ALL, NTHETA_EXTENDED))
    CALL RANDOM_NUMBER(RND_MAT_AUX1)
    CALL RANDOM_NUMBER(RND_MAT_AUX2)
    RND_MAT_VIRT1 = CMPLX(RND_MAT_AUX1, RND_MAT_AUX2,q)
    CALL RANDOM_NUMBER(RND_MAT_AUX1)
    CALL RANDOM_NUMBER(RND_MAT_AUX2)
    RND_MAT_VIRT2 = CMPLX(RND_MAT_AUX1, RND_MAT_AUX2,q)
    DEALLOCATE(RND_MAT_AUX1, RND_MAT_AUX2)
#endif
#ifdef PROFILING
    CALL STOP_PROFILING('RANDOM_NUMBER')
#endif

#ifdef gammareal
    RND_MAT_OCC1(:,:)  = (RND_MAT_OCC1(:,:)  - 0.5_q) * 2.0_q * SQRT(3.0_q) 
    RND_MAT_OCC2(:,:)  = (RND_MAT_OCC2(:,:)  - 0.5_q) * 2.0_q * SQRT(3.0_q) 
    RND_MAT_VIRT1(:,:) = (RND_MAT_VIRT1(:,:) - 0.5_q) * 2.0_q * SQRT(3.0_q) 
    RND_MAT_VIRT2(:,:) = (RND_MAT_VIRT2(:,:) - 0.5_q) * 2.0_q * SQRT(3.0_q) 
#else
    RND_MAT_OCC1(:,:)  = (RND_MAT_OCC1(:,:)  - CMPLX(0.5_q, 0.5_q,q)) * 2.0_q * SQRT(3.0_q/2.0_q) 
    RND_MAT_OCC2(:,:)  = (RND_MAT_OCC2(:,:)  - CMPLX(0.5_q, 0.5_q,q)) * 2.0_q * SQRT(3.0_q/2.0_q) 
    RND_MAT_VIRT1(:,:) = (RND_MAT_VIRT1(:,:) - CMPLX(0.5_q, 0.5_q,q)) * 2.0_q * SQRT(3.0_q/2.0_q) 
    RND_MAT_VIRT2(:,:) = (RND_MAT_VIRT2(:,:) - CMPLX(0.5_q, 0.5_q,q)) * 2.0_q * SQRT(3.0_q/2.0_q) 
#endif

    ! multiply the orbital energies ( exp(+/- epsilon tau / 2) ) into the random matrices (time evolution)
    DO I = 1, NB_OCC
      FAC = EXP(+TAU * 0.5_q * ( REAL(W%CELTOT(I,K1,ISP) - EFERMI(ISP) ) ) ) * W%FERTOT(I,K1,ISP)
      RND_MAT_OCC1(I,:)  = RND_MAT_OCC1(I,:) * FAC
      RND_MAT_OCC2(I,:)  = RND_MAT_OCC2(I,:) * FAC
    ENDDO
    DO A = 1, NB_ALL
      IF (ABS(1- W%FERTOT(A,K2,ISP)) < 1E-8) THEN
        RND_MAT_VIRT1(A,:) = 0
        RND_MAT_VIRT2(A,:) = 0
      ELSE
        FAC = EXP(-TAU * 0.5_q * ( REAL(W%CELTOT(A,K2,ISP) - EFERMI(ISP) ) ) ) * (1 - W%FERTOT(A,K2,ISP) )
        RND_MAT_VIRT1(A,:) = RND_MAT_VIRT1(A,:) * FAC
        RND_MAT_VIRT2(A,:) = RND_MAT_VIRT2(A,:) * FAC
      ENDIF
    ENDDO

    ! created subarrays for the considered k-points and spin  
    WI = ELEMENTS(WOCC, WDESK1, ISP)
    WA = ELEMENTS(W, WDESK2, ISP)
    WRNDI1 = ELEMENTS(WRNDOCC1, WDESK1, ISP)
    WRNDI2 = ELEMENTS(WRNDOCC2, WDESK1, ISP)
    WRNDA1 = ELEMENTS(WRNDVIRT1, WDESK2, ISP)
    WRNDA2 = ELEMENTS(WRNDVIRT2, WDESK2, ISP)

    ! calculate random occupied orbitals |kappa> and |lambda>: WRNDI%CPTWFP = WI%CPTWFP * RND_MAT_OCC
    CALL GGEMM('N', 'N', m_ WDESK1%NPL_RED, WRNDOCC1%WDES%NB_TOT, WDESK1%NB_TOT, one, WI%CW_RED(1,1), &
               m_ WDESK1%NRPLWV_RED, RND_MAT_OCC1, NB_OCC, zero, WRNDI1%CW_RED(1,1), m_ WDESK1%NRPLWV_RED)
    CALL GGEMM('N', 'N', WDESK1%NPRO_RED, WRNDOCC1%WDES%NB_TOT, WDESK1%NB_TOT, one, WI%CPROJ_RED(1,1), &
               WDESK1%NPROD_RED, RND_MAT_OCC1, NB_OCC, zero, WRNDI1%CPROJ_RED(1,1), WDESK1%NPROD_RED)

    CALL GGEMM('N', 'N', m_ WDESK1%NPL_RED, WRNDOCC2%WDES%NB_TOT, WDESK1%NB_TOT, one, WI%CW_RED(1,1), &
               m_ WDESK1%NRPLWV_RED, RND_MAT_OCC2, NB_OCC, zero, WRNDI2%CW_RED(1,1), m_ WDESK1%NRPLWV_RED)
    CALL GGEMM('N', 'N', WDESK1%NPRO_RED, WRNDOCC2%WDES%NB_TOT, WDESK1%NB_TOT, one, WI%CPROJ_RED(1,1), &
               WDESK1%NPROD_RED, RND_MAT_OCC2, NB_OCC, zero, WRNDI2%CPROJ_RED(1,1), WDESK1%NPROD_RED)

    ! calculate random virtual orbitals |alpha> and |beta>: WRNDA%CPTWFP = WA%CPTWFP * RND_MAT_VIRT
    CALL GGEMM('N', 'N', m_ WDESK2%NPL_RED, WRNDVIRT1%WDES%NB_TOT, WDESK2%NB_TOT, one, WA%CW_RED(1,1), &
               m_ WDESK2%NRPLWV_RED, RND_MAT_VIRT1, NB_ALL, zero, WRNDA1%CW_RED(1,1), m_ WDESK2%NRPLWV_RED)
    CALL GGEMM('N', 'N', WDESK2%NPRO_RED, WRNDVIRT1%WDES%NB_TOT, WDESK2%NB_TOT, one, WA%CPROJ_RED(1,1), &
               WDESK2%NPROD_RED, RND_MAT_VIRT1, NB_ALL, zero, WRNDA1%CPROJ_RED(1,1), WDESK2%NPROD_RED)

    CALL GGEMM('N', 'N', m_ WDESK2%NPL_RED, WRNDVIRT2%WDES%NB_TOT, WDESK2%NB_TOT, one, WA%CW_RED(1,1), &
               m_ WDESK2%NRPLWV_RED, RND_MAT_VIRT2, NB_ALL, zero, WRNDA2%CW_RED(1,1), m_ WDESK2%NRPLWV_RED)
    CALL GGEMM('N', 'N', WDESK2%NPRO_RED, WRNDVIRT2%WDES%NB_TOT, WDESK2%NB_TOT, one, WA%CPROJ_RED(1,1), &
               WDESK2%NPROD_RED, RND_MAT_VIRT2, NB_ALL, zero, WRNDA2%CPROJ_RED(1,1), WDESK2%NPROD_RED)

#ifdef PROFILING
    CALL START_PROFILING('REDISTRIBUTE_ORBITALS')
#endif
    CALL REDISTRIBUTE_PROJ(WRNDI1)
    CALL REDISTRIBUTE_PROJ(WRNDI2)
    CALL REDISTRIBUTE_PROJ(WRNDA1)
    CALL REDISTRIBUTE_PROJ(WRNDA2)
    CALL REDISTRIBUTE_PW(WRNDI1)
    CALL REDISTRIBUTE_PW(WRNDI2)
    CALL REDISTRIBUTE_PW(WRNDA1)
    CALL REDISTRIBUTE_PW(WRNDA2)
#ifdef PROFILING
    CALL STOP_PROFILING('REDISTRIBUTE_ORBITALS')
#endif

    ! do the same for the stochastic virtuals with the lower number of NBANDS (if desired)
    IF( NB_LOW > 0) THEN
      WRNDA1_LOW = ELEMENTS(WRNDVIRT1_LOW, WDESK2, ISP)
      WRNDA2_LOW = ELEMENTS(WRNDVIRT2_LOW, WDESK2, ISP)
      CALL GGEMM('N', 'N', m_ WDESK2%NPL_RED, WRNDVIRT1_LOW%WDES%NB_TOT, NB_LOW, one, WA%CW_RED(1,1), &
                 m_ WDESK2%NRPLWV_RED, RND_MAT_VIRT1, NB_ALL, zero, WRNDA1_LOW%CW_RED(1,1), m_ WDESK2%NRPLWV_RED)
      CALL GGEMM('N', 'N', WDESK2%NPRO_RED, WRNDVIRT1_LOW%WDES%NB_TOT, NB_LOW, one, WA%CPROJ_RED(1,1), &
                 WDESK2%NPROD_RED, RND_MAT_VIRT1, NB_ALL, zero, WRNDA1_LOW%CPROJ_RED(1,1), WDESK2%NPROD_RED)
      CALL GGEMM('N', 'N', m_ WDESK2%NPL_RED, WRNDVIRT2_LOW%WDES%NB_TOT, NB_LOW, one, WA%CW_RED(1,1), &
                 m_ WDESK2%NRPLWV_RED, RND_MAT_VIRT2, NB_ALL, zero, WRNDA2_LOW%CW_RED(1,1), m_ WDESK2%NRPLWV_RED)
      CALL GGEMM('N', 'N', WDESK2%NPRO_RED, WRNDVIRT2_LOW%WDES%NB_TOT, NB_LOW, one, WA%CPROJ_RED(1,1), &
                 WDESK2%NPROD_RED, RND_MAT_VIRT2, NB_ALL, zero, WRNDA2_LOW%CPROJ_RED(1,1), WDESK2%NPROD_RED)
#ifdef PROFILING
      CALL START_PROFILING('REDISTRIBUTE_ORBITALS')
#endif
      CALL REDISTRIBUTE_PROJ(WRNDA1_LOW)
      CALL REDISTRIBUTE_PROJ(WRNDA2_LOW)
      CALL REDISTRIBUTE_PW(WRNDA1_LOW)
      CALL REDISTRIBUTE_PW(WRNDA2_LOW)
#ifdef PROFILING
    CALL STOP_PROFILING('REDISTRIBUTE_ORBITALS')
#endif
    ENDIF

    DEALLOCATE(RND_MAT_OCC1, RND_MAT_VIRT1, RND_MAT_OCC2, RND_MAT_VIRT2)

#ifdef PROFILING
    CALL STOP_PROFILING('CALC_RND_ORBITALS')
#endif
  END SUBROUTINE CALC_RND_ORBITALS




  !************************* CALC_DIRECT_AND_EXCHANGE_MP2 ******************************
  !
  ! This routines calculates the direct and exchange MP2 energy contribution given
  ! the set of stochastic orbitals WRNDOCC1, WRNDOCC2, WRNDVIRT1, and WRNDVIRT2.
  ! In the two loops over THETA1 and THETA2 the relevant overlaps densities,
  ! OVERLAPS1-4, are calculated and multiplied by the Coulomb kernel (4pi/G^2).
  ! Hence, for each combination of THETA1 and THETA2 the following is calcualted:
  !
  ! d = sum_G 4pi/G^2 OVERLAPS1(G) OVERLAPS2(G)*
  ! x = sum_G 4pi/G^2 OVERLAPS3(G) OVERLAPS4(G)*
  !
  ! DMP2 = -2|d|^2
  ! XMP2 = d*x
  !
  ! Furthermore, if NBLOW > 0 is passed, then only the high energy
  ! contributions to the MP2 energy are calculated
  !
  !************************************************************************************
  SUBROUTINE CALC_DIRECT_AND_EXCHANGE_MP2(WRNDOCC1, WRNDOCC2, WRNDVIRT1, WRNDVIRT2, LATT_CUR, WGWQ, NTHETA, &
                                          NCORREL_THETA, CORREL_THETA_START, CORREL_THETA_STOP, RSPIN,  &
                                          REACHED_DI_ACCURACY, REACHED_EX_ACCURACY, EDRCT, EEXCHG, NBLOW, &
                                          WRNDVIRT1_LOW, WRNDVIRT2_LOW)
    USE ini
    TYPE (wavespin) :: WRNDOCC1
    TYPE (wavespin) :: WRNDOCC2
    TYPE (wavespin) :: WRNDVIRT1
    TYPE (wavespin) :: WRNDVIRT2
    TYPE (latt) :: LATT_CUR
    TYPE (wavedes1) :: WGWQ     
    INTEGER :: NTHETA
    INTEGER :: NCORREL_THETA
    INTEGER :: CORREL_THETA_START
    INTEGER :: CORREL_THETA_STOP
    REAL(q) :: RSPIN
    INTEGER :: REACHED_DI_ACCURACY
    INTEGER :: REACHED_EX_ACCURACY
    COMPLEX(q) :: EDRCT(NCORREL_THETA**2, N_MP2_ENCUTS)  ! the resulting direct MP2 contribution
    COMPLEX(q) :: EEXCHG(NCORREL_THETA**2, N_MP2_ENCUTS) ! the resulting exchange MP2 contribution
    INTEGER, OPTIONAL ::  NBLOW
    TYPE (wavespin), OPTIONAL ::  WRNDVIRT1_LOW
    TYPE (wavespin), OPTIONAL ::  WRNDVIRT2_LOW
    ! local
    INTEGER           ::  NB_LOW
    TYPE (wavefun1), ALLOCATABLE :: W1OCC1(:)    ! for random occupied orbitals in real space
    TYPE (wavefun1), ALLOCATABLE :: W1OCC2(:)    ! for random occupied orbitals in real space
    TYPE (wavefun1), ALLOCATABLE :: W1VIRT1(:)   ! for random virtual orbitals in real space
    TYPE (wavefun1), ALLOCATABLE :: W1VIRT2(:)   ! for random virtual orbitals in real space
    TYPE (wavespin) :: WHF                       ! temp wavespin with WDES_FOCK
    TYPE (wavedes1), TARGET :: WDESK             ! for gamma only we need only one wavedes1
    INTEGER :: K, ISP
    INTEGER :: N, T, T1, T2
    REAL(q) :: FSG
    REAL(q) :: POTFAK(WGWQ%NPL, N_MP2_ENCUTS)
    GDEF, ALLOCATABLE :: GWORK(:,:)              ! work array for overlap density
    GDEF, ALLOCATABLE :: CRHOLM(:,:)             ! augmentation occupancy matrix
    COMPLEX(q), ALLOCATABLE :: OVERLAPS1(:,:)
    COMPLEX(q), ALLOCATABLE :: OVERLAPS2(:,:)
    COMPLEX(q), ALLOCATABLE :: OVERLAPS3(:,:)
    COMPLEX(q), ALLOCATABLE :: OVERLAPS4(:,:)
    COMPLEX(q) :: ETMP1(N_MP2_ENCUTS)
    COMPLEX(q) :: ETMP2(N_MP2_ENCUTS)

#ifdef PROFILING
    CALL START_PROFILING('CALC_DIRECT_AND_EXCHANGE_MP2')
#endif

    NB_LOW = -1 
    IF(PRESENT(NBLOW)) NB_LOW = NBLOW 

    K = 1   ! GAMMA ONLY
    ISP = 1 ! NO SPIN POL.

    CALL CHECK_FULL_KPOINTS

    ! fourier transform all stochastic orbitals to real space, and start with...
    ! ...WRNDOCC1 to W1OCC1
    WHF = WRNDOCC1
    WHF%WDES => WDES_FOCK
    WHF%WDES%NBANDS = WRNDOCC1%WDES%NBANDS
    WHF%WDES%NB_TOT = WRNDOCC1%WDES%NB_TOT
    ALLOCATE(W1OCC1(NCORREL_THETA))
    CALL SETWDES(WHF%WDES, WDESK, 0)
    DO N = 1, NCORREL_THETA
      CALL NEWWAV(W1OCC1(N) , WDESK, .TRUE.)
    ENDDO
    CALL SETWDES(WHF%WDES, WDESK, K)
    DO N = 1, NCORREL_THETA
      T = CORREL_THETA_START + N - 1
      CALL W1_COPY(ELEMENT(WHF, WDESK, T, ISP), W1OCC1(N))
      CALL FFTWAV_W1(W1OCC1(N))
    ENDDO

    ! ...WRNDOCC2 to W1OCC2
    WHF = WRNDOCC2
    WHF%WDES => WDES_FOCK
    WHF%WDES%NBANDS = WRNDOCC2%WDES%NBANDS
    WHF%WDES%NB_TOT = WRNDOCC2%WDES%NB_TOT
    ALLOCATE(W1OCC2(NCORREL_THETA))
    CALL SETWDES(WHF%WDES, WDESK, 0)
    DO N = 1, NCORREL_THETA
      CALL NEWWAV(W1OCC2(N) , WDESK, .TRUE.)
    ENDDO
    CALL SETWDES(WHF%WDES, WDESK, K)
    DO N = 1, NCORREL_THETA
      T = CORREL_THETA_START + N - 1
      CALL W1_COPY(ELEMENT(WHF, WDESK, T, ISP), W1OCC2(N))
      CALL FFTWAV_W1(W1OCC2(N))
    ENDDO

    ! ...WRNDVIRT1 to W1VIRT1
    WHF = WRNDVIRT1
    WHF%WDES => WDES_FOCK
    WHF%WDES%NBANDS = WRNDVIRT1%WDES%NBANDS
    WHF%WDES%NB_TOT = WRNDVIRT1%WDES%NB_TOT
    ALLOCATE(W1VIRT1(NCORREL_THETA))
    CALL SETWDES(WHF%WDES, WDESK, 0)
    DO N = 1, NCORREL_THETA
      CALL NEWWAV(W1VIRT1(N) , WDESK, .TRUE.)
    ENDDO
    CALL SETWDES(WHF%WDES, WDESK, K)
    DO N = 1, NCORREL_THETA
      T = CORREL_THETA_START + N - 1
      CALL W1_COPY(ELEMENT(WHF, WDESK, T, ISP), W1VIRT1(N))
      CALL FFTWAV_W1(W1VIRT1(N))
    ENDDO

    ! ...WRNDVIRT2 to W1VIRT2
    WHF = WRNDVIRT2
    WHF%WDES => WDES_FOCK
    WHF%WDES%NBANDS = WRNDVIRT2%WDES%NBANDS
    WHF%WDES%NB_TOT = WRNDVIRT2%WDES%NB_TOT
    ALLOCATE(W1VIRT2(NCORREL_THETA))
    CALL SETWDES(WHF%WDES, WDESK, 0)
    DO N = 1, NCORREL_THETA
      CALL NEWWAV(W1VIRT2(N) , WDESK, .TRUE.)
    ENDDO
    CALL SETWDES(WHF%WDES, WDESK, K)
    DO N = 1, NCORREL_THETA
      T = CORREL_THETA_START + N - 1
      CALL W1_COPY(ELEMENT(WHF, WDESK, T, ISP), W1VIRT2(N))
      CALL FFTWAV_W1(W1VIRT2(N))
    ENDDO


    ! allocate work arrays and overlap densities
    ALLOCATE( GWORK(m_ GRIDHF%MPLWV, NCORREL_THETA), CRHOLM(AUG_DES%NPROD*WHF%WDES%NRSPINORS, NCORREL_THETA), &
              OVERLAPS1(WGWQ%NPL, NCORREL_THETA), &
              OVERLAPS2(WGWQ%NPL, NCORREL_THETA), &
              OVERLAPS3(WGWQ%NPL, NCORREL_THETA), &
              OVERLAPS4(WGWQ%NPL, NCORREL_THETA))

    ! register memory of real space orbtials
#ifdef gammareal
       CALL REGISTER_ALLOCATE(4*NCORREL_THETA*(16._q*SIZE(W1OCC1(1)%CPTWFP,KIND=selected_int_kind(15))          &
                                              + 8._q*SIZE(W1OCC1(1)%CPROJ,KIND=selected_int_kind(15))           &
                                              +16._q*SIZE(W1OCC1(1)%CR,KIND=selected_int_kind(15))), "realsp")
       CALL REGISTER_ALLOCATE(8.0_q*SIZE(GWORK)+8.0_q*SIZE(CRHOLM)+4*16.0_q*SIZE(OVERLAPS1), "ovrlps")
#else
       CALL REGISTER_ALLOCATE(4*NCORREL_THETA*(16._q*SIZE(W1OCC1(1)%CPTWFP,KIND=selected_int_kind(15))          &
                                              +16._q*SIZE(W1OCC1(1)%CPROJ,KIND=selected_int_kind(15))           &
                                              +16._q*SIZE(W1OCC1(1)%CR,KIND=selected_int_kind(15))), "realsp")
       CALL REGISTER_ALLOCATE(16.0_q*SIZE(GWORK)+16.0_q*SIZE(CRHOLM)+4*16.0_q*SIZE(OVERLAPS1), "ovrlps")
#endif

    FSG = 0 ! at the moment no non-analytic for G = 0
    CALL RND_MP2_SET_GFAC_VKPT(WGWQ, LATT_CUR, FSG, POTFAK)
   
    ! calc OVERLAPS1 = < kappa_theta | G | alpha_theta >
    DO T = 1, NCORREL_THETA
      CALL FOCK_CHARGE(W1OCC1(T), W1VIRT1(T), GWORK(:,T), CRHOLM(:,T))
    ENDDO
    DO T = 1, NCORREL_THETA
      CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), GWORK(1,T), OVERLAPS1(1,T), WGWQ%GRID, .FALSE.)
    ENDDO
    OVERLAPS1 = OVERLAPS1 / (GRIDHF%NPLWV * 1.0_q)

    ! calc OVERLAPS2 = < lambda_theta | G | beta_theta >
    DO T = 1, NCORREL_THETA
      CALL FOCK_CHARGE(W1VIRT2(T), W1OCC2(T), GWORK(:,T), CRHOLM(:,T))
    ENDDO
    DO T = 1, NCORREL_THETA
      CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), GWORK(1,T), OVERLAPS2(1,T), WGWQ%GRID, .FALSE.)
    ENDDO
    OVERLAPS2 = OVERLAPS2 / (GRIDHF%NPLWV * 1.0_q)

    theta2: &
    DO T2 = 1, NCORREL_THETA

        IF(REACHED_EX_ACCURACY == 0) THEN
          ! calc OVERLAPS3 = < alpha_theta1 | G | lambda_theta2 >
          CALL FOCK_CHARGE_MU(W1VIRT1, W1OCC2(T2), GWORK, CRHOLM)
          DO T1 = 1, NCORREL_THETA
            CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), GWORK(1,T1), OVERLAPS3(1,T1), WGWQ%GRID, .FALSE.)
          ENDDO
          OVERLAPS3 = OVERLAPS3 / (GRIDHF%NPLWV * 1.0_q)

          ! calc OVERLAPS4 = < kappa_theta1 | G | beta_theta2 >
          CALL FOCK_CHARGE_MU(W1OCC1, W1VIRT2(T2), GWORK, CRHOLM)
          DO T1 = 1, NCORREL_THETA
            CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), GWORK(1,T1), OVERLAPS4(1,T1), WGWQ%GRID, .FALSE.)
          ENDDO
          OVERLAPS4 = OVERLAPS4 / (GRIDHF%NPLWV * 1.0_q)
        ENDIF

        theta1: &
        DO T1 = 1, NCORREL_THETA

          CALL APPLY_GFAC_2(OVERLAPS1(:,T1), OVERLAPS2(:,T2), POTFAK, WGWQ%NPL, ETMP1)
          IF(REACHED_EX_ACCURACY == 0) &
            CALL APPLY_GFAC_2(OVERLAPS3(:,T1), OVERLAPS4(:,T1), POTFAK, WGWQ%NPL, ETMP2)

#ifdef gammareal
          ETMP1 = REAL(ETMP1,q)
          IF(REACHED_EX_ACCURACY == 0) &
            ETMP2 = REAL(ETMP2,q)
#endif
          EDRCT( (T1-1)*NCORREL_THETA+T2,:) = -ABS(ETMP1(:))**2 
          IF(REACHED_EX_ACCURACY == 0) &
            EEXCHG((T1-1)*NCORREL_THETA+T2,:) = ETMP1(:) * ETMP2(:)

        ENDDO theta1
    ENDDO theta2

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! if NB_LOW > 0 is passed, then we want to calculate only the high energy
    ! contribution to MP2 -> calculate and substract the low energy MP2 contribution
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IF( NB_LOW > 0) THEN

      ! set the Coulomb kernel using the ENCUT_LOW
      CALL RND_MP2_SET_GFAC_VKPT(WGWQ, LATT_CUR, FSG, POTFAK, LOW_CUT=.TRUE.)

      DO N = 1, NCORREL_THETA
         CALL DELWAV(W1VIRT1(N), .TRUE.)
         CALL DELWAV(W1VIRT2(N), .TRUE.)
      ENDDO

      ! fourier transform virtual low band stochastic orbitals to real space, and start with...
      ! ...WRNDVIRT1_LOW to W1VIRT1
      WHF = WRNDVIRT1_LOW
      WHF%WDES => WDES_FOCK
      WHF%WDES%NBANDS = WRNDVIRT1_LOW%WDES%NBANDS
      WHF%WDES%NB_TOT = WRNDVIRT1_LOW%WDES%NB_TOT
      CALL SETWDES(WHF%WDES, WDESK, 0)
      DO N = 1, NCORREL_THETA
        CALL NEWWAV(W1VIRT1(N) , WDESK, .TRUE.)
      ENDDO
      CALL SETWDES(WHF%WDES, WDESK, K)
      DO N = 1, NCORREL_THETA
        T = CORREL_THETA_START + N - 1
        CALL W1_COPY(ELEMENT(WHF, WDESK, T, ISP), W1VIRT1(N))
        CALL FFTWAV_W1(W1VIRT1(N))
      ENDDO

      ! ...WRNDVIRT2_LOW to W1VIRT2
      WHF = WRNDVIRT2_LOW
      WHF%WDES => WDES_FOCK
      WHF%WDES%NBANDS = WRNDVIRT2_LOW%WDES%NBANDS
      WHF%WDES%NB_TOT = WRNDVIRT2_LOW%WDES%NB_TOT
      CALL SETWDES(WHF%WDES, WDESK, 0)
      DO N = 1, NCORREL_THETA
        CALL NEWWAV(W1VIRT2(N) , WDESK, .TRUE.)
      ENDDO
      CALL SETWDES(WHF%WDES, WDESK, K)
      DO N = 1, NCORREL_THETA
        T = CORREL_THETA_START + N - 1
        CALL W1_COPY(ELEMENT(WHF, WDESK, T, ISP), W1VIRT2(N))
        CALL FFTWAV_W1(W1VIRT2(N))
      ENDDO

      ! calc OVERLAPS1 = < kappa_theta | G | alpha_theta >
      DO T = 1, NCORREL_THETA
        CALL FOCK_CHARGE(W1OCC1(T), W1VIRT1(T), GWORK(:,T), CRHOLM(:,T))
      ENDDO
      DO T = 1, NCORREL_THETA
        CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), GWORK(1,T), OVERLAPS1(1,T), WGWQ%GRID, .FALSE.)
      ENDDO
      OVERLAPS1 = OVERLAPS1 / (GRIDHF%NPLWV * 1.0_q)

      ! calc OVERLAPS2 = < lambda_theta | G | beta_theta >
      DO T = 1, NCORREL_THETA
        CALL FOCK_CHARGE(W1VIRT2(T), W1OCC2(T), GWORK(:,T), CRHOLM(:,T))
      ENDDO
      DO T = 1, NCORREL_THETA
        CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), GWORK(1,T), OVERLAPS2(1,T), WGWQ%GRID, .FALSE.)
      ENDDO
      OVERLAPS2 = OVERLAPS2 / (GRIDHF%NPLWV * 1.0_q)

      theta2_: &
      DO T2 = 1, NCORREL_THETA

          IF(REACHED_EX_ACCURACY == 0) THEN
            ! calc OVERLAPS3 = < alpha_theta1 | G | lambda_theta2 >
            CALL FOCK_CHARGE_MU(W1VIRT1, W1OCC2(T2), GWORK, CRHOLM)
            DO T1 = 1, NCORREL_THETA
              CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), GWORK(1,T1), OVERLAPS3(1,T1), WGWQ%GRID, .FALSE.)
            ENDDO
            OVERLAPS3 = OVERLAPS3 / (GRIDHF%NPLWV * 1.0_q)

            ! calc OVERLAPS4 = < kappa_theta1 | G | beta_theta2 >
            CALL FOCK_CHARGE_MU(W1OCC1, W1VIRT2(T2), GWORK, CRHOLM)
            DO T1 = 1, NCORREL_THETA
              CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), GWORK(1,T1), OVERLAPS4(1,T1), WGWQ%GRID, .FALSE.)
            ENDDO
            OVERLAPS4 = OVERLAPS4 / (GRIDHF%NPLWV * 1.0_q)
          ENDIF

          theta1_: &
          DO T1 = 1, NCORREL_THETA

            CALL APPLY_GFAC_2(OVERLAPS1(:,T1), OVERLAPS2(:,T2), POTFAK, WGWQ%NPL, ETMP1)
            IF(REACHED_EX_ACCURACY == 0) &
              CALL APPLY_GFAC_2(OVERLAPS3(:,T1), OVERLAPS4(:,T1), POTFAK, WGWQ%NPL, ETMP2)

#ifdef gammareal
            ETMP1 = REAL(ETMP1,q)
            IF(REACHED_EX_ACCURACY == 0) &
              ETMP2 = REAL(ETMP2,q)
#endif
            EDRCT( (T1-1)*NCORREL_THETA+T2,:) = EDRCT( (T1-1)*NCORREL_THETA+T2,:) + ABS(ETMP1(1))**2 
            IF(REACHED_EX_ACCURACY == 0) &
              EEXCHG((T1-1)*NCORREL_THETA+T2,:) = EEXCHG((T1-1)*NCORREL_THETA+T2,:) - ETMP1(1) * ETMP2(1)

          ENDDO theta1_
      ENDDO theta2_
    ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    EDRCT   = EDRCT  * RSPIN**2 / 2.0_q
    EEXCHG = EEXCHG * RSPIN   / 2.0_q

    !EEXCHG = EDRCT + EEXCHG

    ! deallocate
    DEALLOCATE(GWORK, CRHOLM, OVERLAPS1, OVERLAPS2, OVERLAPS3, OVERLAPS4)
    DO N = 1, NCORREL_THETA
       CALL DELWAV(W1OCC1(N), .TRUE.)
       CALL DELWAV(W1OCC2(N), .TRUE.)
       CALL DELWAV(W1VIRT1(N), .TRUE.)
       CALL DELWAV(W1VIRT2(N), .TRUE.)
    ENDDO
    DEALLOCATE(W1OCC1)
    DEALLOCATE(W1OCC2)
    DEALLOCATE(W1VIRT1)
    DEALLOCATE(W1VIRT2)

#ifdef PROFILING
    CALL STOP_PROFILING('CALC_DIRECT_AND_EXCHANGE_MP2')
#endif

  END SUBROUTINE CALC_DIRECT_AND_EXCHANGE_MP2




  !********************** RND_MP2_SET_GFAC_VKPT  ******************************
  !
  ! This routine provides the array POTFAK containing
  ! 4 pi e^2 / (G+q)^2
  ! on the plan-wave grid given by WGWQ multiplied by 
  ! a smooth cosinus cutoff between ENCUTSOFT and ENCUT
  !
  ! additionally set the zero grid point value to the convergence
  ! correction FSG if G+q = 0 
  !
  !********************************************************************
  SUBROUTINE RND_MP2_SET_GFAC_VKPT(WGWQ, LATT_CUR, FSG, POTFAK, LOW_CUT)
    USE constant

    TYPE(wavedes1) :: WGWQ
    TYPE(latt) :: LATT_CUR
    REAL(q) :: FSG
    REAL(q) :: POTFAK(WGWQ%NPL, N_MP2_ENCUTS)
    LOGICAL, OPTIONAL :: LOW_CUT    ! switch if only the low energy MP2 contributions are considered
    ! local
    INTEGER :: NG, I
    REAL(q) :: GX, GY, GZ, GSQU
    REAL(q) :: DQX, DQY, DQZ
    REAL(q) :: EN
    REAL(q) :: SCALE
    LOGICAL :: LOW_CUTOFF              ! internally used switch

    PROFILING_START( 'rnd_mp2_set_gfac_vkpt' )

    LOW_CUTOFF = .FALSE.
    IF( PRESENT( LOW_CUT ) ) LOW_CUTOFF = LOW_CUT 

    SCALE = EDEPS / LATT_CUR%OMEGA / TPI**2

    DQX = WGWQ%VKPT(1) * LATT_CUR%B(1,1) + WGWQ%VKPT(2) * LATT_CUR%B(1,2) + WGWQ%VKPT(3) * LATT_CUR%B(1,3)
    DQY = WGWQ%VKPT(1) * LATT_CUR%B(2,1) + WGWQ%VKPT(2) * LATT_CUR%B(2,2) + WGWQ%VKPT(3) * LATT_CUR%B(2,3)
    DQZ = WGWQ%VKPT(1) * LATT_CUR%B(3,1) + WGWQ%VKPT(2) * LATT_CUR%B(3,2) + WGWQ%VKPT(3) * LATT_CUR%B(3,3)

    ! check if only the low energy contributions are considered (ENCUT_LOW)...

    IF( LOW_CUTOFF ) THEN
      POTFAK = 0
      DO NG = 1, WGWQ%NPL

        ! calculate (G+q)^2
        GX = WGWQ%IGX(NG)*LATT_CUR%B(1,1) + WGWQ%IGY(NG)*LATT_CUR%B(1,2) + WGWQ%IGZ(NG)*LATT_CUR%B(1,3)
        GY = WGWQ%IGX(NG)*LATT_CUR%B(2,1) + WGWQ%IGY(NG)*LATT_CUR%B(2,2) + WGWQ%IGZ(NG)*LATT_CUR%B(2,3)
        GZ = WGWQ%IGX(NG)*LATT_CUR%B(3,1) + WGWQ%IGY(NG)*LATT_CUR%B(3,2) + WGWQ%IGZ(NG)*LATT_CUR%B(3,3)
        GSQU = (GX+DQX)**2 + (GY+DQY)**2 + (GZ+DQZ)**2

        ! set FSG if G+q = 0
        IF (GSQU < G2ZERO) THEN
          POTFAK(NG,:) = FSG
          CYCLE
        ENDIF
      
        ! calculate corresponding energy to GSQU
        EN = HSQDTM * GSQU * TPI**2 

        IF (EN <= ENCUT_LOW) THEN
          POTFAK(NG,1) = SCALE / GSQU 
        ELSE 
          POTFAK(NG,1) = 0
        ENDIF
      ENDDO
    ! ...or if the full Coulomb kernel is desired (ENCUT and ENCUTSOFT)
    ELSE
      DO NG = 1, WGWQ%NPL

        ! calculate (G+q)^2
        GX = WGWQ%IGX(NG)*LATT_CUR%B(1,1) + WGWQ%IGY(NG)*LATT_CUR%B(1,2) + WGWQ%IGZ(NG)*LATT_CUR%B(1,3)
        GY = WGWQ%IGX(NG)*LATT_CUR%B(2,1) + WGWQ%IGY(NG)*LATT_CUR%B(2,2) + WGWQ%IGZ(NG)*LATT_CUR%B(2,3)
        GZ = WGWQ%IGX(NG)*LATT_CUR%B(3,1) + WGWQ%IGY(NG)*LATT_CUR%B(3,2) + WGWQ%IGZ(NG)*LATT_CUR%B(3,3)
        GSQU = (GX+DQX)**2 + (GY+DQY)**2 + (GZ+DQZ)**2

        ! set FSG if G+q = 0
        IF (GSQU < G2ZERO) THEN
          POTFAK(NG,:) = FSG
          CYCLE
        ENDIF
      
        ! calculate corresponding energy to GSQU
        EN = HSQDTM * GSQU * TPI**2 

        DO I = 1, N_MP2_ENCUTS
          ! check in which cutoff area we are
          IF (EN <= ENCUTSOFTS(I)) THEN
            POTFAK(NG,I) = SCALE / GSQU ! set 4 pi e^2 / (G+q)^2 if we are below the smooth cutoff threshold
          ELSE IF (EN > ENCUTS(I)) THEN
            POTFAK(NG,I) = 0
          ELSE
            POTFAK(NG,I) = ( 1 + COS((EN-ENCUTSOFTS(I))/(ENCUTS(I)-ENCUTSOFTS(I))*PI) ) * SCALE / ( 2 * GSQU )
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    PROFILING_STOP( 'rnd_mp2_set_gfac_vkpt' )
  END SUBROUTINE RND_MP2_SET_GFAC_VKPT 




  !********************** ORBITALS_TO_REAL_SPACE  ******************************
  !
  !********************************************************************
  SUBROUTINE ORBITALS_TO_REAL_SPACE(W, NB1, NB2, ISP, W1)
    IMPLICIT NONE
    TYPE (wavespin) W        ! wavefunction
    INTEGER :: NB1           ! starting band
    INTEGER :: NB2           ! final band
    INTEGER :: ISP           ! spin
    TYPE (wavefun1):: W1(:)  ! array into which the merge is performed
    ! local
    INTEGER :: N_INTO_TOT, N, NCPU

    PROFILING_START('RND_ORBITALS_TO_REAL_SPACE')

    NCPU = W%WDES%NB_PAR
    DO N = (NB1-1)/W%WDES%NB_PAR+1, (NB2-1)/W%WDES%NB_PAR+1
      N_INTO_TOT = (N-1)*NCPU + W%WDES%NB_LOW
      IF (NB1<=N_INTO_TOT .AND. N_INTO_TOT<=NB2) THEN
        CALL W1_COPY(ELEMENT(W, W1(N_INTO_TOT-NB1+1)%WDES1, N, ISP), W1(N_INTO_TOT-NB1+1) )
        CALL FFTWAV_W1(W1(N_INTO_TOT-NB1+1))
      ENDIF
    ENDDO

    PROFILING_STOP('RND_ORBITALS_TO_REAL_SPACE')

  END SUBROUTINE ORBITALS_TO_REAL_SPACE




  !********************** APPLY_GFAC_2  ******************************
  !
  ! This routine calculates
  ! sum_G POTFAK(G) CWORK(G) CWORKG(G)^*
  ! and stores the result in EN
  !
  !********************************************************************
  SUBROUTINE APPLY_GFAC_2(CWORK, CWORKG, POTFAK, NPL, EN, LOW_CUT)
    COMPLEX(q)  :: CWORK(NPL)
    COMPLEX(q)  :: CWORKG(NPL)
    REAL(q)     :: POTFAK(NPL,N_MP2_ENCUTS)
    INTEGER     :: NPL
    COMPLEX(q)  :: EN(N_MP2_ENCUTS)
    LOGICAL, OPTIONAL :: LOW_CUT       ! switch if only the low energy MP2 contributions are considered
    ! local
    INTEGER :: NG, I
    LOGICAL :: LOW_CUTOFF              ! internally used switch

#ifdef PROFILING
    CALL START_PROFILING('APPLY_GFAC_2')
#endif

    EN = 0
    LOW_CUTOFF = .FALSE.
    IF( PRESENT( LOW_CUT ) ) LOW_CUTOFF = LOW_CUT 

    ! check if only the low energy contributions are considered (ENCUT_LOW)...
    IF(LOW_CUTOFF) THEN
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(NG) &
!$OMP REDUCTION(+:EN)
      DO NG = 1, NPL
        EN(1) = EN(1) + CWORK(NG) * DCONJG( CWORKG(NG) ) * POTFAK(NG,1)
      ENDDO
!$OMP END PARALLEL DO
    ! ...or if the full Coulomb kernel is desired (ENCUT and ENCUTSOFT)
    ELSE  
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(NG) &
!$OMP REDUCTION(+:EN)
      DO NG = 1, NPL
        DO I = N_MP2_ENCUTS, 1, -1
          ! check if POTFAKs are equal for different cutoffs (if so treat all POTFAKs equally)
          IF (ABS(POTFAK(NG,I) - POTFAK(NG,1)) < 1.E-12_q ) THEN
            EN(1:I) = EN(1:I) + CWORK(NG) * DCONJG( CWORKG(NG) ) * POTFAK(NG,I)
            EXIT
          ! if POTFAKs differ treat POTFAK(I) individually
          ELSE
            EN(I) = EN(I) + CWORK(NG) * DCONJG( CWORKG(NG) ) * POTFAK(NG,I)
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
    ENDIF

#ifdef PROFILING
    CALL STOP_PROFILING('APPLY_GFAC_2')
#endif

  END SUBROUTINE APPLY_GFAC_2



END MODULE stoch_lt_mp2
#else
MODULE stoch_lt_mp2
  CONTAINS
  SUBROUTINE stoch_lt_mp2_dummy
      WRITE(*,*)'Im a DEC compiler so I need this line'
  END SUBROUTINE
END MODULE stoch_lt_mp2
#endif





