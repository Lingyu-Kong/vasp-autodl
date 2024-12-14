!#define dotiming
#include "symbol.inc"

! we either have the gammareal or the time_inversion mode
#ifndef gammareal
#define time_inversion
#endif
#ifdef MPI
!***********************************************************************
!
! This module distributes and stores the overlap densities for
! the MP2 calculation
!
!***********************************************************************
MODULE overlap_densities
  USE prec
  USE mpimy
  USE ini
  USE WAVE
  IMPLICIT NONE
  
  TYPE (communic), PRIVATE :: COMM_OVERLAP   ! communicator for distribution of plane wave coefficients
  INTEGER, PRIVATE :: G_PER_CORE             ! number of G indices per core 
  INTEGER, PRIVATE :: START_G_INDEX          ! global start G index on local core
  INTEGER, PRIVATE :: FINAL_G_INDEX          ! global final G index on local core

  ! OVERLAP_DENS stores OVERLAP_DENS(G, i, a, k, sigma) = < u_a,k | e^-iGr e^iK(k+q)r | u_i,T(k+q) > 
  ! note that the occupied indices i are GLOBAL indices whereas the unoccupied a are LOCAL !
  COMPLEX(q), PRIVATE, ALLOCATABLE, TARGET :: OVERLAP_DENS(:,:,:,:,:) 

  CONTAINS

  
  ! sets the communicator for the parallelization over G (first argument)
  SUBROUTINE OVERLAP_SET_COMM(COMM)
    TYPE (communic) COMM
    COMM_OVERLAP = COMM
  END SUBROUTINE OVERLAP_SET_COMM




  !****************************** OVERLAP_DISTRIBUTE_G  ******************************
  !
  ! establishes the parallelization strategy of the G vectors over the cores
  ! 
  !***********************************************************************************
  SUBROUTINE OVERLAP_DISTRIBUTE_G(WGWQ, IU0, IU6)
    TYPE (wavedes1) :: WGWQ
    INTEGER :: IU0, IU6
    INTEGER :: G_OVERRUN

    G_PER_CORE = (WGWQ%NPL + COMM_OVERLAP%NCPU - 1) / COMM_OVERLAP%NCPU 
    START_G_INDEX = G_PER_CORE * (COMM_OVERLAP%NODE_ME - 1) + 1
    FINAL_G_INDEX = G_PER_CORE *  COMM_OVERLAP%NODE_ME 
    
    G_OVERRUN = FINAL_G_INDEX - WGWQ%NPL ! calculate the overrun (how many useless local G indices there are)
    IF (G_OVERRUN < 0) G_OVERRUN = 0

    CALLMPI(M_sum_i(COMM_OVERLAP, G_OVERRUN, 1))

    !IF (IU0>=0) THEN
    !  WRITE(IU0,*) " "
    !  WRITE(IU0,*) " G_PER_CORE, WGWQ%NPL, COMM_OVERLAP%NCPU"
    !  WRITE(IU0,*) G_PER_CORE, WGWQ%NPL, COMM_OVERLAP%NCPU
    !  WRITE(IU0,'(A,F6.2,A)') " Processor power NOT contributing to plane waves inside of the cutoff (plane wave parallelization inefficiency):", 100.0_q * G_OVERRUN / WGWQ%NPL, "%"
    !ENDIF
    !IF (IU6>=0) THEN
    !  WRITE(IU6,*) " "
    !  WRITE(IU6,*) " G_PER_CORE, WGWQ%NPL, COMM_OVERLAP%NCPU"
    !  WRITE(IU6,*) G_PER_CORE, WGWQ%NPL, COMM_OVERLAP%NCPU
    !  WRITE(IU6,'(A,F6.2,A)') " Processor power NOT contributing to plane waves inside of the cutoff (plane wave parallelization inefficiency):", 100.0_q * G_OVERRUN / WGWQ%NPL, "%"
    !ENDIF
  END SUBROUTINE OVERLAP_DISTRIBUTE_G




  !****************************** OVERLAP_ALLOCATE  ******************************
  !
  ! allocate memory or the overlap densities
  ! 
  !*******************************************************************************
  SUBROUTINE OVERLAP_ALLOCATE(LAST_FILLED, NKPTS, NBANDS, ISPIN, IU6)
    INTEGER :: NBANDS      ! number of orbitals treated locally on each node
    INTEGER :: LAST_FILLED ! total number (global) of occupied orbitals
    INTEGER :: NKPTS       ! number of k-points in BZ
    INTEGER :: ISPIN       ! number of spins
    INTEGER :: IU6

    ALLOCATE(OVERLAP_DENS(G_PER_CORE, LAST_FILLED, NBANDS, NKPTS, ISPIN))
    CALL REGISTER_ALLOCATE(16._q*SIZE(OVERLAP_DENS), "overlap")
    CALL DUMP_ALLOCATE(IU6)
    OVERLAP_DENS = (0, 0)
  END SUBROUTINE OVERLAP_ALLOCATE




  !****************************** OVERLAP_DEALLOCATE  ******************************
  !
  ! deallocates the memory for the overlap densites
  ! 
  !***********************************************************************************
  SUBROUTINE OVERLAP_DEALLOCATE(IU6)
    INTEGER :: IU6

    CALL DEREGISTER_ALLOCATE(16._q*SIZE(OVERLAP_DENS), "overlap")
    DEALLOCATE(OVERLAP_DENS)
    !CALL DUMP_ALLOCATE(IU6)
  END SUBROUTINE OVERLAP_DEALLOCATE




  !****************************** OVERLAP_STORE ******************************
  !
  ! this routine gets an overlap density in reciprocal space and picks and stores  
  ! all G associated with this core 
  ! 
  !***************************************************************************
  SUBROUTINE OVERLAP_STORE(OVRLP, NI, NA, K2, ISP)
    COMPLEX(q) :: OVRLP(:)    ! overlap density in reciprocal space
    INTEGER :: NI             ! GLOBAL occupied index
    INTEGER :: NA             ! LOCAL unoccupied index 
    INTEGER :: K2, ISP        ! k-point and spin
    ! local 
    INTEGER :: TRUNC_START_G
    INTEGER :: TRUNC_FINAL_G
    INTEGER :: SORTED_NG, LOCAL_NG
    IF ( NI > SIZE(OVERLAP_DENS,2) .OR. NA > SIZE(OVERLAP_DENS,3) .OR. K2 > SIZE(OVERLAP_DENS,4)) THEN
       CALL vtutor%bug("internal error in OVERLAP_STORE: bounds exceeded " // str(NI) // " " // &
          str(SIZE(OVERLAP_DENS,2)) // " " // str(NA) // " " // str(SIZE(OVERLAP_DENS,3)) // " " // &
          str(K2) // " " // str(SIZE(OVERLAP_DENS,4)), __FILE__, __LINE__)
    ENDIF
    ! ToDO tS: rethink this more carefully
    TRUNC_START_G = MIN(START_G_INDEX, SIZE(OVRLP))
    TRUNC_FINAL_G = MIN(FINAL_G_INDEX, SIZE(OVRLP))
    !TRUNC_START_G = START_G_INDEX
    !TRUNC_FINAL_G = FINAL_G_INDEX
    OVERLAP_DENS(1:(TRUNC_FINAL_G-TRUNC_START_G+1), NI, NA, K2, ISP) = OVRLP(TRUNC_START_G : TRUNC_FINAL_G)
  END SUBROUTINE OVERLAP_STORE




  !****************************** OVERLAP_GET_GLOBAL_G ******************************
  !
  ! returns the global G index given a local G index
  ! 
  !*************************************************************************************
  FUNCTION OVERLAP_GET_GLOBAL_G(NG_LOCAL)
    INTEGER :: NG_LOCAL
    INTEGER :: OVERLAP_GET_GLOBAL_G
    IF(NG_LOCAL > G_PER_CORE  .OR. NG_LOCAL < 1) THEN
      OVERLAP_GET_GLOBAL_G = -1
    ELSE
      OVERLAP_GET_GLOBAL_G = START_G_INDEX + NG_LOCAL - 1
    ENDIF
  END FUNCTION OVERLAP_GET_GLOBAL_G




  ! returns G_PER_CORE
  FUNCTION OVERLAP_GET_G_PER_CORE()
    INTEGER :: OVERLAP_GET_G_PER_CORE
    OVERLAP_GET_G_PER_CORE = G_PER_CORE
  END FUNCTION OVERLAP_GET_G_PER_CORE




  !****************************** OVERLAP_RETURN_I ******************************
  !
  ! returns the stored overlap density for all GLOBAL occupied states i, for given 
  ! LOCAL unoccupied orbital a, and given spin and given k2
  ! 
  !******************************************************************************
  FUNCTION OVERLAP_RETURN_I(NG_LOCAL, N2, K2, ISP)
    COMPLEX (q), POINTER :: OVERLAP_RETURN_I(:)
    INTEGER :: NG_LOCAL
    INTEGER :: N2       ! LOCAL unoccupied index a
    INTEGER :: K2, ISP

    OVERLAP_RETURN_I=>OVERLAP_DENS(NG_LOCAL,:, N2, K2, ISP)
  END FUNCTION OVERLAP_RETURN_I




  !****************************** FIND_GLOBAL_G_INDEX ******************************
  !
  ! returns the VASP index for G given a 3D reciprocal vector
  ! 
  !*********************************************************************************
  FUNCTION FIND_GLOBAL_G_INDEX(IGX, IGY, IGZ, WGWQ)
    INTEGER :: FIND_GLOBAL_G_INDEX
    INTEGER :: IGX, IGY, IGZ
    TYPE (wavedes1)    WGWQ
    ! local
    INTEGER :: N

    FIND_GLOBAL_G_INDEX = 0
    DO N = 1, WGWQ%NPL
      IF (WGWQ%IGX(N) == IGX .AND. WGWQ%IGY(N) == IGY .AND. WGWQ%IGZ(N) == IGZ) THEN
        FIND_GLOBAL_G_INDEX = N
      ENDIF
    ENDDO
    
    IF (FIND_GLOBAL_G_INDEX == 0) THEN
      CALL vtutor%bug("internal error in FIND_GLOBAL_G_INDEX: wave vector not found " // str(IGX) &
         // " " // str(IGY) // " " // str(IGZ), __FILE__, __LINE__)
    ENDIF
  END FUNCTION FIND_GLOBAL_G_INDEX

END MODULE overlap_densities






!***********************************************************************
!
! this module implements the laplace transformed MP2
! based on an Schmidt decomposition
!
!***********************************************************************
MODULE lt_mp2
  USE prec
  USE fock
  USE chi_base
  USE lattice
  USE full_kpoints
  USE paw
  USE overlap_densities
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

CONTAINS
    
  !****************************** LAPLACE_TRANSFORMED_MP2 ******************************
  !
  ! This is the main routine that calculates the MP2 energy
  ! 
  !*************************************************************************************
  SUBROUTINE LAPLACE_TRANSFORMED_MP2(LATT_CUR, W, LMDIM, CQIJ, NONLR_S, NONL_S, P, TAU, TAUWEIGHT, WGW, ENCUT, ENCUTSOFT, IU0, IU6)
    USE ini
    USE nonl_high
    USE wave_high
    USE mpimy
    USE mgrid
    USE constant
    USE pseudo
    USE kpoints_change, ONLY : KPOINTS_ORIG
    USE tutor, ONLY: vtutor, isAlert, BSEantires
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
    REAL(q) :: ENCUT, ENCUTSOFT                   ! cutoff for outer sums over G (actually only the potential is set to 0)
    INTEGER :: IU0, IU6
    ! local
    REAL(q) :: E_EXCH(N_MP2_ENCUTS)            ! exchange energy
    REAL(q) :: EXCHANGE_TOTAL(N_MP2_ENCUTS) 
    REAL(q) :: EXCHANGE_CONVERGED
    REAL(q) :: E_DIRECT(N_MP2_ENCUTS)          ! direct energy
    REAL(q) :: DIRECT_TOTAL(N_MP2_ENCUTS)
    REAL(q) :: DIRECT_CONVERGED
    INTEGER :: I, NK, NTAU                     ! loop variables
    INTEGER :: LAST_FILLED                     ! last orbital treated as filled
    INTEGER :: NG_LOCAL, NG                    ! G loop variables
    INTEGER :: NQ                              ! difference vector for which MP2 is calculated (outer q-point loop)
    TYPE (wavedes1)    WGWQ                    ! wave descriptor for outer q-point loop
    TYPE (wavedes1)    WDESK                   ! temporary wavedescriptor for NK
    TYPE (wavedes)  :: WDESOCC                 ! descriptor to store filled orbitals
    TYPE (wavespin) :: WT                      ! Schmidt decomposed states |w>
    TYPE (wavespin) :: WOCC                    ! occupied orbitals |i>
    INTEGER :: LREAL
    REAL (q) :: SCALE, GX, GY, GZ, GSQU
    REAL (q) :: NFLOAT, NFFT                   ! counter
    REAL (q) :: B_OVERRUN_PERCENT              ! artificial increase of the occupied bands in per cent
    
#ifdef PROFILING
    CALL START_PROFILING('LAPLACE_TRANSFORMED_MP2')
#endif

    !status output
    IF (IU0 >= 0) THEN
       WRITE(IU0,'(A)') " " 
       WRITE(IU0,'(A)') " Starting the LT-MP2 calculation" 
       WRITE(IU0,'(A)') " -------------------------------" 
    ENDIF

    ! allocate memory for the different ENCUTs
    ALLOCATE(ENCUTS(N_MP2_ENCUTS))
    ALLOCATE(ENCUTSOFTS(N_MP2_ENCUTS))
    ! calculate different ENCUTs
    DO I = 1, N_MP2_ENCUTS
      ENCUTS(I) = ENCUT / ( ENCUT_DIVISOR**(I-1) )
      ENCUTSOFTS(I) = ENCUTSOFT / ( ENCUT_DIVISOR**(I-1) )
    ENDDO


    ! set to communictor for the overlap densities (first argument)
    CALL OVERLAP_SET_COMM(W%WDES%COMM_KINTER) ! here we missuse the k-parallelization flag KPAR

    ! first calculate the last filled and first empty orbital indices 
    ! among all k-points and spins to make live simple
    LAST_FILLED = 0
    DO I = 1, W%WDES%ISPIN
      DO NK = 1,W%WDES%NKPTS
        LAST_FILLED = MAX(LAST_FILLED_XI_NOMOD(W, NK, I), LAST_FILLED)
      ENDDO
    ENDDO

    WDESOCC=W%WDES

    ! round WDESOCC%NBANDS to next larger value dividable by WDESOCC%NB_PAR
    WDESOCC%NBANDS = (LAST_FILLED + WDESOCC%NB_PAR - 1) / WDESOCC%NB_PAR
    WDESOCC%NB_TOT = WDESOCC%NBANDS * WDESOCC%NB_PAR

    ! calculate inefficieny of parallelization over bands
    CALL B_PARALLEL_EFFICIENCY(WDESOCC%NBANDS, LAST_FILLED, WDESOCC%NB_PAR, B_OVERRUN_PERCENT, IU0, IU6)

    IF (IU6>=0) WRITE(IU6,'(/A,I6 /)') ' Bands included in the LTMP2 as occupied', LAST_FILLED

    ! now create orbital array for the original orbitals
    ! but restrict to occupied manyfold
    CALL ALLOCW(WDESOCC, WOCC) ! occupied orbitals
    CALL ALLOCW(WDESOCC, WT)   ! Schmidt decomposed states
#ifndef gammareal
    ! need LORBITALREAL to be true, print warning otherwise
    IF (.NOT. W%WDES%LORBITALREAL) THEN
       CALL vtutor%write(isAlert, BSEantires)
    ENDIF
#endif

    ! copy original orbitals into this "occupied" wave function array
    WOCC%CW(:,1:WOCC%WDES%NBANDS,:,:)    = W%CW(:,1:WOCC%WDES%NBANDS,:,:)
    WOCC%CPROJ(:,1:WOCC%WDES%NBANDS,:,:) = W%CPROJ(:,1:WOCC%WDES%NBANDS,:,:)
    WOCC%CELTOT(1:WOCC%WDES%NB_TOT,:,:)  = W%CELTOT(1:WOCC%WDES%NB_TOT,:,:)
    WOCC%FERTOT(1:WOCC%WDES%NB_TOT,:,:)  = W%FERTOT(1:WOCC%WDES%NB_TOT,:,:)
    WOCC%AUXTOT(1:WOCC%WDES%NB_TOT,:,:)  = W%AUXTOT(1:WOCC%WDES%NB_TOT,:,:)

    ! scaling for Coulomb kernel
    SCALE = EDEPS / (LATT_CUR%OMEGA * TPI**2)

    EXCHANGE_TOTAL = 0
    DIRECT_TOTAL=0

    ! outer q-point loop
    ! loop over all q-points in the ireducible wedge
    qpt: & 
    DO NQ = 1, KPOINTS_ORIG%NKPTS

      !status output
      IF (IU0 >= 0) THEN
         WRITE(IU0,'(A)') " " 
         WRITE(IU0,'(A,1I3,A,1I3)') " Working on NQ = ", NQ, " of ", KPOINTS_ORIG%NKPTS 
      ENDIF

      CALL SETWDES(WGW, WGWQ, NQ)

      ! determine data distribution of G vectors among groups for the overlap densities
      CALL OVERLAP_DISTRIBUTE_G(WGWQ, IU0, IU6)
      
      ! calculate how much cpu power is for nothing due to parallelization over plane waves
      CALL G_PARALLEL_EFFICIENCY(OVERLAP_GET_G_PER_CORE(), WGWQ%NPL, W%WDES%COMM_KINTER%NCPU, NQ, B_OVERRUN_PERCENT, IU0, IU6)

      IF (NQ == 1) CALL START_TIMING("ovrlp") ! measure time only for the first q-point

      ! allocate required storage for overlap densities
      CALL OVERLAP_ALLOCATE(LAST_FILLED, W%WDES%NKPTS, W%WDES%NBANDS, W%WDES%ISPIN, IU6)

      !status output
      IF (IU0 >= 0) THEN
         WRITE(IU0,'(A$)') "    Calculating the overlap densities..." 
      ENDIF

      ! calculate overlap density  for all spins and distribute over all cores
      CALL CALCULATE_OVERLAP_DENSITY(W, LATT_CUR, LAST_FILLED, NQ, WGWQ, W%WDES%ISPIN, P, IU6)

      IF (NQ == 1) CALL STOP_TIMING("ovrlp", IU6) ! measure time only for the first q-point

      ! status output
      IF (IU0 >= 0) THEN
         WRITE(IU0,'(A)') " done" 
      ENDIF

      ! redistribute original orbitals over plane wave coefficients 
      ! necessary for SCHMIDT_DECOMPOSE: there we need all orbitals at each core 
      ! note that we DON NOT redistribute the occupied orbitals WOCC%CW and WOCC%CPROJ
      DO I = 1, W%WDES%ISPIN
        DO NK = 1, W%WDES%NKPTS
          CALL SETWDES(W%WDES, WDESK, NK)            
          CALL REDIS_PW(WDESK, W%WDES%NBANDS, W%CW(1, 1, NK, I))
          CALL REDIS_PROJ(WDESK, W%WDES%NBANDS, W%CPROJ(1, 1, NK, I))
        ENDDO
      ENDDO 

      NFLOAT = 0
      NFFT = 0

      !status output
      IF (IU0 >= 0) THEN
         WRITE(IU0,'(A)') "    Calculating the MP2 contribution for NQ" 
      ENDIF

      tau_loop: &
      DO NTAU=1,SIZE(TAU)
       
        !status output
        IF (IU0 >= 0) THEN
           WRITE(IU0,'(A,1I3,A,1I3)') "       Working on TAU = ", NTAU, " of ", SIZE(TAU)
        ENDIF

        ! outer G loop
        DO NG_LOCAL = 1, OVERLAP_GET_G_PER_CORE()

          ! calculate (G+q)^2
          NG = OVERLAP_GET_GLOBAL_G(NG_LOCAL)  
          IF(NG > WGWQ%NPL .OR. NG < 1) CYCLE
          GX = (WGWQ%IGX(NG)+WGWQ%VKPT(1))*LATT_CUR%B(1,1) &
             + (WGWQ%IGY(NG)+WGWQ%VKPT(2))*LATT_CUR%B(1,2) & 
             + (WGWQ%IGZ(NG)+WGWQ%VKPT(3))*LATT_CUR%B(1,3)
          GY = (WGWQ%IGX(NG)+WGWQ%VKPT(1))*LATT_CUR%B(2,1) &
             + (WGWQ%IGY(NG)+WGWQ%VKPT(2))*LATT_CUR%B(2,2) &
             + (WGWQ%IGZ(NG)+WGWQ%VKPT(3))*LATT_CUR%B(2,3)
          GZ = (WGWQ%IGX(NG)+WGWQ%VKPT(1))*LATT_CUR%B(3,1) &
             + (WGWQ%IGY(NG)+WGWQ%VKPT(2))*LATT_CUR%B(3,2) &
             + (WGWQ%IGZ(NG)+WGWQ%VKPT(3))*LATT_CUR%B(3,3)
          GSQU = GX**2 + GY**2 + GZ**2

          ! first calculate the Schmidt decomposed states (WT) then
          ! evaluate the direct and exchange diagram for this G, q and tau
#if defined(gammareal)
          DO LREAL = 0, 1
            CALL SCHMIDT_DECOMPOSE(W, WT, WDESOCC, LAST_FILLED, NQ, NG_LOCAL, TAU(NTAU), W%WDES%ISPIN, NFLOAT, LREAL == 0 )
            CALL CALCULATE_DIRECT(WOCC, WT, LATT_CUR, E_DIRECT, ENCUTS, ENCUTSOFTS, WGWQ, NQ, W%WDES%ISPIN, NFFT, P)
            CALL CALCULATE_EXCHANGE(WOCC, WT, LAST_FILLED, LATT_CUR, E_EXCH, ENCUTS, ENCUTSOFTS, WGWQ, NQ, W%WDES%ISPIN, NFFT, P)
            
            ! loop over all different ENCUT and ENCUTSOFT
            DO I = 1, N_MP2_ENCUTS
              EXCHANGE_TOTAL(I) = EXCHANGE_TOTAL(I) - E_EXCH(I) &
                                * COULOMB_ATTENUATED(GSQU, SCALE, ENCUTSOFTS(I), ENCUTS(I)) &
                                * TAUWEIGHT(NTAU) * W%WDES%RSPIN * KPOINTS_ORIG%WTKPT(NQ)
              DIRECT_TOTAL(I) =  DIRECT_TOTAL(I) - E_DIRECT(I) &
                              * COULOMB_ATTENUATED(GSQU, SCALE, ENCUTSOFTS(I), ENCUTS(I))&
                              * TAUWEIGHT(NTAU) * KPOINTS_ORIG%WTKPT(NQ)
            ENDDO
          ENDDO
#elif defined(time_inversion)
          CALL SCHMIDT_DECOMPOSE(W, WT, WDESOCC, LAST_FILLED, NQ, NG_LOCAL, TAU(NTAU), W%WDES%ISPIN, NFLOAT)
          CALL CALCULATE_DIRECT(WOCC, WT, LATT_CUR, E_DIRECT, ENCUTS, ENCUTSOFTS, WGWQ, NQ, W%WDES%ISPIN, NFFT, P)
          CALL CALCULATE_EXCHANGE(WOCC, WT, LAST_FILLED, LATT_CUR, E_EXCH, ENCUTS, ENCUTSOFTS, WGWQ, NQ, W%WDES%ISPIN, NFFT, P)

          ! loop over all different ENCUT and ENCUTSOFT
          DO I = 1, N_MP2_ENCUTS
            EXCHANGE_TOTAL(I) = EXCHANGE_TOTAL(I) - E_EXCH(I) &
                              * COULOMB_ATTENUATED(GSQU, SCALE, ENCUTSOFTS(I), ENCUTS(I)) &
                              * TAUWEIGHT(NTAU) * W%WDES%RSPIN * KPOINTS_ORIG%WTKPT(NQ)
            DIRECT_TOTAL(I) =  DIRECT_TOTAL(I) - E_DIRECT(I) &
                            * COULOMB_ATTENUATED(GSQU, SCALE, ENCUTSOFTS(I), ENCUTS(I)) &
                            * TAUWEIGHT(NTAU) * KPOINTS_ORIG%WTKPT(NQ)
          ENDDO
#endif
        ENDDO

        !CALL STOP_TIMING("G", IU6, "LTMP2")
        CALL WFORCE(IU6)

      ENDDO tau_loop

      CALL OVERLAP_DEALLOCATE(IU6)
      
      ! redistribute original orbitals back over bands (required for next call of CALCULATE_OVERLAP_DENSITY)
      ! note that we DON NOT redistribute the occupied orbitals WOCC%CW and WOCC%CPROJ
      DO I = 1, W%WDES%ISPIN
        DO NK = 1, W%WDES%NKPTS
          CALL SETWDES(W%WDES, WDESK, NK)            
          CALL REDIS_PW(WDESK, W%WDES%NBANDS, W%CW(1, 1, NK, I))
          CALL REDIS_PROJ(WDESK, W%WDES%NBANDS, W%CPROJ(1, 1, NK, I))
        ENDDO
      ENDDO

    ENDDO qpt

    CALL DEALLOCW(WOCC)
    CALL DEALLOCW(WT)

    ! collect results from each core
    CALLMPI(M_sum_d(W%WDES%COMM_KINTER, EXCHANGE_TOTAL, N_MP2_ENCUTS))
    CALLMPI(M_sum_d(W%WDES%COMM_KINTER, DIRECT_TOTAL, N_MP2_ENCUTS))

    ! perform linear regression to extrapolate ENCUT -> infinity
    CALL LIN_REG_CONVERGE_LTMP2(REAL(EXCHANGE_TOTAL, KIND=q), ENCUTS, EXCHANGE_CONVERGED)
    CALL LIN_REG_CONVERGE_LTMP2(REAL(DIRECT_TOTAL, KIND=q), ENCUTS, DIRECT_CONVERGED)

#ifdef PROFILING
    CALL STOP_PROFILING('LAPLACE_TRANSFORMED_MP2')
#endif

    ! output
    IF (IU6>=0) THEN
      WRITE(IU6,*) " "
      WRITE(IU6,*) " -------------------------------------------------------------------------"
      WRITE(IU6,'(5A)') "      ENCUTGW", "  ENCUTGWSOFT", "      direct MP2", "    exchange MP2", "       total MP2"
      WRITE(IU6,*) " -------------------------------------------------------------------------"
      DO I = N_MP2_ENCUTS, 1, -1
        WRITE(IU6,'(2F13.3,3F16.8)') ENCUTS(I), ENCUTSOFTS(I), REAL(DIRECT_TOTAL(I)), REAL(EXCHANGE_TOTAL(I)), &
                                     REAL(DIRECT_TOTAL(I)+EXCHANGE_TOTAL(I))
      ENDDO
      WRITE(IU6,'(A)') "      linear regression"
      WRITE(IU6,'(A,3F16.8)') "      converged values    ", DIRECT_CONVERGED, EXCHANGE_CONVERGED, DIRECT_CONVERGED+EXCHANGE_CONVERGED
    ENDIF
    IF (IU0>=0) THEN
      WRITE(IU0,*) " "
      WRITE(IU0,*) " -------------------------------------------------------------------------"
      WRITE(IU0,'(5A)') "      ENCUTGW", "  ENCUTGWSOFT", "      direct MP2", "    exchange MP2", "       total MP2"
      WRITE(IU0,*) " -------------------------------------------------------------------------"
      DO I = N_MP2_ENCUTS, 1, -1
        WRITE(IU0,'(2F13.3,3F16.8)') ENCUTS(I), ENCUTSOFTS(I), REAL(DIRECT_TOTAL(I)), REAL(EXCHANGE_TOTAL(I)), &
                                     REAL(DIRECT_TOTAL(I)+EXCHANGE_TOTAL(I))
      ENDDO
      WRITE(IU0,'(A)') "      linear regression"
      WRITE(IU0,'(A,3F16.8)') "      converged values    ", DIRECT_CONVERGED, EXCHANGE_CONVERGED, DIRECT_CONVERGED+EXCHANGE_CONVERGED
    ENDIF

    DEALLOCATE(ENCUTS)
    DEALLOCATE(ENCUTSOFTS)

    !CALLMPI( M_sum_d(W%WDES%COMM,NFLOAT,1))
    !CALLMPI( M_sum_d(W%WDES%COMM,NFFT,1))
    !IF (IU6>=0) WRITE(IU6,"(' number of  BLAS level 3 operations ',F10.2,' Gflops',F10.0,' kfft')") NFLOAT/1E9, NFFT/1E3
    !IF (IU0>=0) WRITE(IU0,"(' number of  BLAS level 3 operations ',F10.2,' Gflops',F10.0,' kfft')") NFLOAT/1E9, NFFT/1E3

     ! close XML file, dump profiling and stop
     CALL XML_CLOSE_TAG("calculation") ! for whatever reason this is still open
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

  END SUBROUTINE LAPLACE_TRANSFORMED_MP2




  !****************************** G_PARALLEL_EFFICIENCY ******************************
  !
  ! The efficiency of the parallelization over plane waves is estimated (how many 
  ! useless G above the cutoff are calculated?) 
  ! 
  !********************************************************************************
  SUBROUTINE G_PARALLEL_EFFICIENCY(G_PER_CORE, NPL, G_PAR, NQ, B_OVERRUN_PERCENT, IU0, IU6)
    INTEGER :: G_PER_CORE        ! plane waves per core
    INTEGER :: NPL               ! total number of plane waves
    INTEGER :: G_PAR             ! number of G vectors treatet in parallel
    INTEGER :: NQ                ! current q-point (outer k-point loop)
    REAL(q) :: B_OVERRUN_PERCENT ! artificial increase of the occupied bands in per cent
    INTEGER :: IU0, IU6
    ! local
    INTEGER :: G_OVERRUN                  ! how many G vectors are too much
    REAL(q) :: G_OVERRUN_PERCENT 
    REAL(q) :: COMP_TIME_INCREASE_PERCENT ! increase of computation time due to inefficiency of parallelization of
                                          ! plane waves AND bands

    G_OVERRUN = G_PAR * G_PER_CORE - NPL
    G_OVERRUN_PERCENT = 100.0_q * G_OVERRUN / NPL

    COMP_TIME_INCREASE_PERCENT = G_OVERRUN_PERCENT + B_OVERRUN_PERCENT + B_OVERRUN_PERCENT*G_OVERRUN_PERCENT / 100 
    !COMP_TIME_INCREASE_PERCENT = G_OVERRUN_PERCENT + 2*B_OVERRUN_PERCENT &
    !                           + (2*B_OVERRUN_PERCENT*G_OVERRUN_PERCENT + B_OVERRUN_PERCENT**2) / 100 &
    !                           + (B_OVERRUN_PERCENT**2*G_OVERRUN_PERCENT) / 10000


    IF(COMP_TIME_INCREASE_PERCENT >= 5.0) THEN
      IF(IU0 >= 0) THEN
        WRITE(IU0, '(/A)') "    Warning! The inefficiency of the parallelization is not negligible (change KPAR to solve this problem)."
        WRITE(IU0, '(A,F6.2,A)')  "      Plane wave parallelization: effectless CPU power at this NQ: ", G_OVERRUN_PERCENT, "%"
        WRITE(IU0, '(A,F6.2,A)')  "      Band parallelization: artificial increase of last occupied band: ", B_OVERRUN_PERCENT, "%"
        WRITE(IU0, '(A,F6.2,A/)') "      Increase of computation time due to parallelization inefficiency: ", COMP_TIME_INCREASE_PERCENT, "%"
      ENDIF
      IF(IU6 >= 0) THEN
        WRITE(IU6, '(/A,I3,A)') " Warning! At NQ = ", NQ, &
                           ", the inefficiency of the parallelization is not negligible (change KPAR to solve this problem)."
        WRITE(IU6, '(A,F6.2,A)')  "   Plane wave parallelization: effectless CPU power at this NQ: ", G_OVERRUN_PERCENT, "%"
        WRITE(IU6, '(A,F6.2,A)')  "   Band parallelization: artificial increase of last occupied band: ", B_OVERRUN_PERCENT, "%"
        WRITE(IU6, '(A,F6.2,A/)') "   Increase of computation time due to parallelization inefficiency: ", COMP_TIME_INCREASE_PERCENT, "%"
      ENDIF
    ELSE IF (IU6 >= 0) THEN
      WRITE(IU6, '(/A,I3)') " Inefficiency of parallelization over plane waves at NQ = ", NQ
      WRITE(IU6, '(A,F6.2,A/)') "   Effectless CPU power: ", G_OVERRUN_PERCENT, "%"
    ENDIF
  END SUBROUTINE G_PARALLEL_EFFICIENCY




  !****************************** B_PARALLEL_EFFICIENCY ******************************
  !
  ! The efficiency of the parallelization over occupied bands is estimated (is the  
  ! number of occupied bands increased artificially?)
  ! 
  !********************************************************************************
  SUBROUTINE B_PARALLEL_EFFICIENCY(NB_PER_CORE, LAST_FILLED, NB_PAR, B_OVERRUN_PERCENT, IU0, IU6)
    INTEGER :: NB_PER_CORE       ! bands per core
    INTEGER :: LAST_FILLED       ! total number of occupied bands
    INTEGER :: NB_PAR            ! number of bands vectors treatet in parallel
    REAL(q) :: B_OVERRUN_PERCENT ! artificial increase of the occupied bands in per cent
    INTEGER :: IU0, IU6
    ! local
    INTEGER :: B_OVERRUN         ! number of excessive occupied bands

    B_OVERRUN = NB_PAR * NB_PER_CORE - LAST_FILLED
    B_OVERRUN_PERCENT = 100.0_q * B_OVERRUN / LAST_FILLED

    IF(IU6 >= 0) THEN
      WRITE(IU6, '(/A)') " Inefficiency of parallelization over occupied bands:"
      WRITE(IU6, '(A,F6.2,A/)') "   Artificial increase of last occupied band: ", B_OVERRUN_PERCENT, "%"
    ENDIF

    !IF(IU0 >= 0) THEN
    !  IF(B_OVERRUN_PERCENT > 2.5) THEN
    !    WRITE(IU0, '(/A)') " Warning: The inefficiency of the parallelization over occupied bands is not negligible:"
    !    WRITE(IU0, '(A,F6.2,A/)') "   Artificial increase of last occupied band: ", B_OVERRUN_PERCENT, "%"
    !  ENDIF
    !ENDIF
  END SUBROUTINE B_PARALLEL_EFFICIENCY



  !****************************** COULOMB_ATTENUATED ******************************
  !
  ! returns Coulomb kernel 
  !
  ! 4 pi e^2 / G^2
  !
  ! for one given G multiplied by a smooth cosinus cutoff between ENCUTSOFT and ENCUT
  ! 
  !********************************************************************************
  FUNCTION COULOMB_ATTENUATED(GSQU, SCALE, ENCUTSOFT, ENCUT)
    USE constant

    REAL(q) :: COULOMB_ATTENUATED
    REAL(q) :: GSQU
    REAL(q) :: SCALE
    REAL(q) :: ENCUTSOFT, ENCUT
    ! local
    REAL(q) :: EN

    IF (ABS(GSQU) < G2ZERO) THEN
      !COULOMB_ATTENUATED = FSG_STORE(1)
      COULOMB_ATTENUATED = 0
      RETURN
    ENDIF
    
    ! calculate corresponding energy to GSQU
    EN = HSQDTM * GSQU * TPI**2 
    
    ! return 4 pi e^2 / G^2 if we are below the smooth cutoff threshold
    IF(ENCUTSOFT >= 0 .AND. EN <= ENCUTSOFT) THEN 
     COULOMB_ATTENUATED = SCALE / GSQU 

    ! return 0 if we are above the cutoff 
    ELSE IF(ENCUT >= 0 .AND. EN >= ENCUT) THEN
      COULOMB_ATTENUATED = 0

    ! apply smooth cutoff since we are between ENCUTSOFT and ENCUT 
    ELSE
      COULOMB_ATTENUATED = ( 1 + COS((EN-ENCUTSOFT)/(ENCUT-ENCUTSOFT)*PI) ) * SCALE / ( 2 * GSQU ) 
    ENDIF
  END FUNCTION COULOMB_ATTENUATED




  !****************************** CALCULATE_OVERLAP_DENSITY  ******************************
  !
  ! this routine calculates the overlap densities < u_a,k2 | e^-iGr e^iK(k2+q)r | u_i,T(k2+q) > 
  ! in reciprocal space for all a, i, k2, and G and given q
  ! 
  ! This routine can make no use of the parallelization over plane-waves. However, we have to assume that the
  ! number of groups parallelizing the plane-waves is quite high. Thus this routine would calculate a lot of
  ! redundant information. Hence, we need to implement a further parallelization over WHF%WDES%NBANDS.
  ! (see NA_PER_CORE)
  !  
  !*******************************************************************************
  SUBROUTINE CALCULATE_OVERLAP_DENSITY(W, LATT_CUR, LAST_FILLED, NQ, WGWQ, ISPIN, P, IU6)
    TYPE (wavespin) :: W                       ! all orbitals
    TYPE (latt) :: LATT_CUR
    INTEGER :: LAST_FILLED                     ! number of occupied states
    INTEGER :: NQ                              ! given q-point
    TYPE (wavedes1)    WGWQ                    ! wavedescriptor of the given q-point NQ
    INTEGER :: ISPIN                           ! number of spins
    TYPE (potcar) :: P(:)
    INTEGER :: IU6
    ! local variables
    TYPE (wavedes1), TARGET :: WDESK1, WDESK2, WDESK2_IRZ
    TYPE (wavefun1) :: W2                      ! unoccupied orbitals
    TYPE (wavefun1),ALLOCATABLE :: WOCC1(:)    ! occupied orbitals
    INTEGER K1, K2                             ! K1 = T(K2 + q)
    INTEGER N, N_, NA, NA_, NP, NGLB, ISP, ISP_IRZ
    LOGICAL LSHIFT
    COMPLEX(q),ALLOCATABLE :: OVERLAP_I(:,:)   ! overlap density in reciprocal space for fixed i 
    GDEF, ALLOCATABLE :: GWORK(:,:)            ! work array for overlap density
    GDEF, ALLOCATABLE :: CRHOLM(:,:)           ! augmentation occupancy matrix
    TYPE(rotation_handle), POINTER :: ROT_HANDLE
    TYPE (wavespin) WHF                        ! like W but with WDES => WDES_FOCK
    INTEGER IBLK,NBLK,NBSTART,NBSTOP           ! blocking variables
    LOGICAL LPHASE
    COMPLEX(q) :: CPHASE(WGWQ%GRID%MPLWV)      ! phase factor in real space
    INTEGER :: NA_PER_CORE, NA_START, NA_FINAL ! additional parallelization over WHF%WDES%NBANDS
    INTEGER :: CPU                             ! CPU LOOP
    ! test: dump orbitals and overlaps
    !INTEGER(qi8) :: J
    !REAL(q) :: GX, GY, GZ, GSQU, ENEG

#ifdef PROFILING
    CALL START_PROFILING('CALCULATE_OVERLAP_DENSITY')
#endif
    IF(AEXX==0) THEN
#ifdef PROFILING
      CALL STOP_PROFILING('CALCULATE_OVERLAP_DENSITY')
#endif
      RETURN
    ENDIF

    ! use temporarily another WDES
    WHF = W
    WHF%WDES => WDES_FOCK
    CALL CHECK_FULL_KPOINTS

    NULLIFY(ROT_HANDLE)

    ! initialize blocking of occupied orbitals
    NBLK = ((WHF%WDES%NSIM*2+WHF%WDES%NB_PAR-1)/WHF%WDES%NB_PAR)*W%WDES%NB_PAR
    !IF (NBLOCK_FOCK>0) NBLK = MIN(NBLK, NBLOCK_FOCK)  ! deactivated since NBLOCK_FOCK is to small if
                                                       ! OMP parallelization is activated
    NBLK = MIN(NBLK, 64) ! instead manually limit to 64
    
    ! allocate memory
    ALLOCATE( OVERLAP_I(WGWQ%NPL, WHF%WDES%COMM_KINTER%NCPU),  &
              GWORK(m_ GRIDHF%MPLWV,NBLK),                     &
              CRHOLM(AUG_DES%NPROD*WHF%WDES%NRSPINORS,NBLK) )
             
    CALL REGISTER_ALLOCATE(16._q*SIZE(OVERLAP_I), "ovrlp_i")
    CALL DUMP_ALLOCATE(IU6)

    CALL SETWDES(WHF%WDES, WDESK2, 0)
    CALL NEWWAV(W2, WDESK2,.TRUE.)
    CALL SETWDES(WHF%WDES, WDESK1, 0)
    ALLOCATE(WOCC1(NBLK))
    DO N=1, NBLK
      CALL NEWWAV(WOCC1(N) , WDESK1, .TRUE.)
    ENDDO

    ! This routine can make no use of the parallelization over plane-waves. However, we have to assume that the
    ! number of groups parallelizing the plane-waves is quite high. Thus this routine would calculate a lot of
    ! redundant information. Hence, we need to implement a further parallelization over WHF%WDES%NBANDS.
    NA_PER_CORE = (WHF%WDES%NBANDS + WHF%WDES%COMM_KINTER%NCPU - 1) / WHF%WDES%COMM_KINTER%NCPU
    NA_START = NA_PER_CORE * (WHF%WDES%COMM_KINTER%NODE_ME - 1) + 1
    NA_FINAL = NA_PER_CORE *  WHF%WDES%COMM_KINTER%NODE_ME 

    ! k-loop where we have the relationship k1 = T(k2+q) or k2 = T(k1-q) respectively
    kpoint: &
    DO K1=1, WHF%WDES%NKPTS
    
      CALL SETWDES(WHF%WDES, WDESK1, K1)
      
      ! calculate k2 = T(k1-q)
      K2 = KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K1) - KPOINTS_FULL%VKPT(:,NQ), KPOINTS_FULL)

      IF(KPOINTS_FULL%WTKPT(K2)==0 .OR.                                    &
        (HFKIDENT.AND.SKIP_THIS_KPOINT_IN_FOCK(WHF%WDES%VKPT(:,K2))) .OR. &
        (.NOT.HFKIDENT.AND.SKIP_THIS_KPOINT_IN_FOCK(KPOINTS_FULL%VKPT(:,K2)-WHF%WDES%VKPT(:,K1)))) CYCLE

      ! set phase factor for augmentation charges
      CALL PHASER_HF(GRIDHF, LATT_CUR, FAST_AUG_FOCK, WHF%WDES%VKPT(:,K1)-WHF%WDES%VKPT(:,K2))

      ! set correct phase: k1-k2-q = K(k2+q)
      CALL SETPHASE(WHF%WDES%VKPT(:,K1)-WHF%WDES%VKPT(:,K2)-KPOINTS_FULL%VKPT(:,NQ), WGWQ%GRID, CPHASE, LPHASE)

      spin: &
      DO ISP = 1, ISPIN

        block: &
        DO IBLK=1, LAST_FILLED, NBLK

          NGLB = MIN(LAST_FILLED-IBLK+1, NBLK)
          ! global index of the first and last band in this block
          NBSTART = IBLK
          NBSTOP  = NBSTART+NGLB-1

          ! fourier transform the bands for which the HF exchange need to be calculated
          ! to real space, then gather to WOCC1
          CALL W1_GATHER_GLB(WHF, NBSTART, NBSTOP, ISP, WOCC1)

          CALL SETWDES(WHF%WDES, WDESK2, K2)
          CALL SETWDES(WHF%WDES, WDESK2_IRZ, KPOINTS_FULL%NEQUIV(K2))

          ISP_IRZ = ISP
          IF (KPOINTS_FULL%SPINFLIP(K2) == 1) THEN
            ISP_IRZ = 3 - ISP
          ENDIF

          ! loop over all bands present on the local core and which are part of the 
          ! additional parallelization  (actually only the virtual bands are necessary 
          ! but including the occupied bands does no harm)
          unocc: &
          DO NA = NA_START, NA_FINAL
            
            ! check if this band exists at all
            IF(NA <= WHF%WDES%NBANDS) THEN
              IF (K2 <= WHF%WDES%NKPTS) THEN
                CALL W1_COPY(ELEMENT(WHF, WDESK2, NA, ISP), W2)
                CALL FFTWAV_W1(W2)
#ifndef gammareal
              ELSE
                ! symmetry must be considered if the wavefunctions for this
                ! k-point K2 are not stored in W
                LSHIFT=.FALSE.
                IF ((ABS(KPOINTS_FULL%TRANS(1,K2))>TINY) .OR. &
                    (ABS(KPOINTS_FULL%TRANS(2,K2))>TINY) .OR. &
                    (ABS(KPOINTS_FULL%TRANS(3,K2))>TINY)) LSHIFT=.TRUE.
                CALL W1_ROTATE_AND_FFT(W2, ELEMENT(WHF, WDESK2_IRZ, NA, ISP_IRZ), ROT_HANDLE, P, LATT_CUR, LSHIFT)
#endif        
              ENDIF

              ! calculate overlap density < u_a,k2 | r > < r | u_i,T(k2+q) >
              CALL FOCK_CHARGE_MU(WOCC1(1:NGLB), W2, GWORK, CRHOLM)
            ELSE
              GWORK = 0 ! the band did not exists
            ENDIF

            occ: &
            DO N = 1, NGLB
              N_ = NBSTART+N-1        ! global storage index of present band
#ifndef gammareal
              IF (LPHASE) CALL APPLY_PHASE(WGWQ%GRID, CPHASE(1), GWORK(1,N), GWORK(1,N))
#endif      
              OVERLAP_I = 0
              
              ! go to reciprocal space and store result in OVERLAP_I (use FFT's for k-point q)
              CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), GWORK(1,N), OVERLAP_I(1, WHF%WDES%COMM_KINTER%NODE_ME), WGWQ%GRID, .FALSE.)
              
              ! sum all overlaps of parallelization together
              ! note: each core in WHF%WDES%COMM_KINTER needs ALL unoccupied bands for the decomposition later
              CALLMPI(M_sum_z(WHF%WDES%COMM_KINTER, OVERLAP_I(1,1), SIZE(OVERLAP_I)))
              
              ! divide by number of grid points
              OVERLAP_I=OVERLAP_I / GRIDHF%NPLWV
              
              ! store the results of all cpus of the additional parallelization
              DO CPU = 1, WHF%WDES%COMM_KINTER%NCPU
                ! calculate the unoccupied index NA_ (note that NA only runs from NA_START to NA_FINAL but we need all bands)
                NA_ = (CPU - WHF%WDES%COMM_KINTER%NODE_ME) * NA_PER_CORE + NA 
                
                IF (NA_ > WHF%WDES%NBANDS) CYCLE
                
                ! store the overlap
                CALL OVERLAP_STORE(OVERLAP_I(:, CPU), N_, NA_, K2, ISP) ! N_ is the GLOBAL i index and NA_ is the LOCAL a index
              ENDDO

            ENDDO occ
          ENDDO unocc
        ENDDO block
      ENDDO spin
    ENDDO kpoint

    DEALLOCATE(OVERLAP_I, GWORK, CRHOLM)
    CALL DEREGISTER_ALLOCATE(16._q*SIZE(OVERLAP_I), "ovrlp_i")
    CALL DEALLOCATE_ROT_HANDLE(ROT_HANDLE)
    CALL DELWAV(W2,.TRUE.)
    DO N=1,NBLK
       CALL DELWAV(WOCC1(N) ,.TRUE.)
    ENDDO
    DEALLOCATE(WOCC1)

#ifdef PROFILING
    CALL STOP_PROFILING('CALCULATE_OVERLAP_DENSITY')
#endif

  END SUBROUTINE CALCULATE_OVERLAP_DENSITY



  
  !***************************** SCHMIDT_DECOMPOSE *******************************
  !
  ! calculates the Schmidt decomposed states | w_i,T(k1-q) > for all i and k1,
  ! for given G, tau and q
  ! 
  !*******************************************************************************
  SUBROUTINE SCHMIDT_DECOMPOSE(W, WT, WDESOCC, LAST_FILLED, NQ, NG_LOCAL, TAU, ISPIN, NFLOAT, LREAL)
    INTEGER :: NG_LOCAL            ! given G vector
    REAL(q) :: TAU                 ! given tau
    TYPE (wavespin) :: W           ! all orbitals
    TYPE (wavespin) :: WT          ! resulting Schmidt decomposed orbitals |w>
    TYPE (wavedes) :: WDESOCC      ! wavedescriptor for occupied indices
    INTEGER :: LAST_FILLED         ! number of occupied states
    INTEGER :: NQ                  ! given q-point
    INTEGER :: ISPIN               ! number of spins
    LOGICAL, OPTIONAL :: LREAL    
    REAL(q) :: NFLOAT 
    ! local
    INTEGER :: K1, N1, K2, N2, N2_GLOBAL, ISP 
    GDEF, ALLOCATABLE :: CL(:,:)       ! matrix with coefficients (overlaps) for Schmidt decomposition  
    TYPE (wavedes1) :: WDESK2
    TYPE (wavefuna)    WA              ! unoccupied orbital
    TYPE (wavefuna)    WTA             ! Schmidt decomposed state

#ifdef PROFILING
    CALL START_PROFILING('SCHMIDT_DECOMPOSE')
#endif

    ! allocate CL(a,i) and create space for WDESOCC%NB_TOT occupied bands for WT to match NPAR
    ! note that CL(a,i) for i > LAST_FILLED will be set to zero, hence WT for i > LAST_FILLED 
    ! will be zero
    ALLOCATE(CL(W%WDES%NB_TOT, WDESOCC%NB_TOT))
    DO ISP = 1, ISPIN
      ! loop over all k-points
      DO K2 = 1, W%WDES%NKPTS
        CALL SETWDES(W%WDES, WDESK2, K2)

        K1 = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,K2)+KPOINTS_FULL%VKPT(:,NQ), KPOINTS_FULL)

        CL = 0

        ! calculate the coefficient matrix
        ! loop over unoccupied bands stored locally
        DO N2 = 1, W%WDES%NBANDS
          N2_GLOBAL=(N2-1)*W%WDES%NB_PAR+W%WDES%NB_LOW
#ifdef gammareal
          IF (PRESENT (LREAL) .AND. LREAL) THEN
            CL(N2_GLOBAL,1:LAST_FILLED)= REAL(OVERLAP_RETURN_I(NG_LOCAL, N2, K2, ISP),q)
          ELSE
            CL(N2_GLOBAL,1:LAST_FILLED)= AIMAG(OVERLAP_RETURN_I(NG_LOCAL, N2, K2, ISP))
          ENDIF
#else
          CL(N2_GLOBAL,1:LAST_FILLED)= OVERLAP_RETURN_I(NG_LOCAL, N2, K2, ISP)
#endif  
          CL(N2_GLOBAL,LAST_FILLED+1: WDESOCC%NB_TOT)=0
          ! loop over occupied indices
          DO N1 = 1, LAST_FILLED
            CL(N2_GLOBAL,N1) = CL(N2_GLOBAL,N1) * EXP(TAU * ( REAL(W%CELTOT(N1,K1,ISP) - W%CELTOT(N2_GLOBAL,K2,ISP),q) )) &
                               * W%FERTOT(N1,K1,ISP)*(1- W%FERTOT(N2_GLOBAL,K2,ISP))
          ENDDO
        ENDDO

        ! gather results from all cores
        CALLMPI(M_sum_g(W%WDES%COMM_KIN, CL(1,1), SIZE(CL)))

        ! created subarrays for the considered k-points K2 and spin ISPIN 
        WA = ELEMENTS(W, WDESK2, ISP)
        WTA = ELEMENTS(WT, WDESK2, ISP)

        ! linear combine orbitals into WTA
        ! matrix mulitplication: WTA%CW_RED = WA%CW_RED * CL
        CALL GGEMM('N', 'N', m_ WDESK2%NPL_RED, WDESOCC%NB_TOT, WDESK2%NB_TOTK(ISP), one, WA%CW_RED(1,1), &
                   m_ WDESK2%NRPLWV_RED, CL, SIZE(CL,1), zero, WTA%CW_RED(1,1), m_ WDESK2%NRPLWV_RED)
        NFLOAT = NFLOAT + m_ WDESK2%NPL_RED * WDESOCC%NB_TOT * WDESK2%NB_TOTK(ISP) * ndata * ndata

        CALL GGEMM('N', 'N', WDESK2%NPRO_RED, WDESOCC%NB_TOT, WDESK2%NB_TOTK(ISP), one, WA%CPROJ_RED(1,1), &
                   WDESK2%NPROD_RED, CL, SIZE(CL,1), zero, WTA%CPROJ_RED(1,1), WDESK2%NPROD_RED)
        NFLOAT = NFLOAT + WDESK2%NPRO_RED * WDESOCC%NB_TOT * WDESK2%NB_TOTK(ISP) * ndata * ndata

        ! redistribute orbitals in WT over bands
        CALL REDISTRIBUTE_PROJ(WTA)
        CALL REDISTRIBUTE_PW(WTA)
      ENDDO
    ENDDO

    DEALLOCATE(CL)

#ifdef PROFILING
    CALL STOP_PROFILING('SCHMIDT_DECOMPOSE')
#endif

  END SUBROUTINE SCHMIDT_DECOMPOSE




  !***************************** CALCULATE_DIRECT  *******************************
  !
  ! calculate the MP2 direct contribution for the current G, q and tau
  !
  ! E_DIRECT = sum_G' v(G'-q) | sum_k1 sum_i  < u_i,k1 | e^-iG'r e^iK(k1-q)r | w_i,T(k1-q) > |^2
  ! 
  !*******************************************************************************
  SUBROUTINE CALCULATE_DIRECT(WOCC, WT, LATT_CUR, E_DIRECT, ENCUTS, ENCUTSOFTS, WGWQ, NQ, ISPIN, NFFT, P)
    TYPE (wavespin) :: WOCC           ! occupied orbitals
    TYPE (wavespin) :: WT             ! Schmidt decomposed orbitals
    TYPE (latt) :: LATT_CUR
    REAL(q) :: E_DIRECT(N_MP2_ENCUTS) ! resulting direct energy
    TYPE (wavedes1)    WGWQ           ! wavedescriptor for given q-point
    INTEGER :: NQ                     ! given q-point
    INTEGER :: ISPIN                  ! number of spins
    REAL(q) :: ENCUTS(N_MP2_ENCUTS)   ! energy cutoffs on Coulomb kernel
    REAL(q) :: ENCUTSOFTS(N_MP2_ENCUTS) ! smooth energy cutoffs on Coulomb kernel
    REAL(q) :: NFFT
    TYPE (potcar) :: P(:)
    ! local variables
    TYPE (wavedes1), TARGET :: WDESK1, WDESK2  !
    TYPE (wavefun1) :: WOCC1                   ! occupied states | u_i,k1 >
    TYPE (wavefun1) :: WT2                     ! Schmidt decomposed states | w_i,T(k1-q) >
    REAL(q) :: FSG                             ! singularity correction setting V(G=0) to a finite value
    INTEGER :: K1, K2, NI, ISP
    GDEF,      ALLOCATABLE :: GWORK(:)         ! work array for overlap densities
    GDEF,      ALLOCATABLE :: GWORKSUM(:)      ! work array for overlap densities
    GDEF,      ALLOCATABLE :: CRHOLM(:)        ! augmentation occupancy matrix
    REAL(q),   ALLOCATABLE :: POTFAKS(:,:)      ! Coulomb potential v(G-q) in reciprocal space
    TYPE (wavespin)        :: WOCCHF, WTHF 
    LOGICAL LPHASE
    COMPLEX(q) :: CPHASE(WGWQ%GRID%MPLWV)      ! phase factors

    IF (AEXX==0) THEN
       RETURN
    ENDIF

#ifdef PROFILING
    CALL START_PROFILING('CALCULATE_DIRECT')
#endif

    ! use temporarily another WDES
    WOCCHF = WOCC
    WOCCHF%WDES => WDES_FOCK
    WTHF = WT
    WTHF%WDES => WDES_FOCK

    CALL CHECK_FULL_KPOINTS

    ! allocate memory
    ALLOCATE( GWORK(m_ GRIDHF%MPLWV), &
              GWORKSUM(m_ GRIDHF%MPLWV), &
              CRHOLM(AUG_DES%NPROD*WOCCHF%WDES%NRSPINORS), &
              POTFAKS(GRIDHF%MPLWV, N_MP2_ENCUTS) )

    ! initialize orbitals and wavedescriptors
    CALL SETWDES(WOCCHF%WDES,WDESK1,0)
    CALL SETWDES(WOCCHF%WDES,WDESK2,0)
    CALL NEWWAV(WOCC1, WDESK1, .TRUE.)
    CALL NEWWAV(WT2, WDESK2, .TRUE.)

    ! average electrostatic potential for k=k' and n=n'
    !FSG=FSG_STORE(1)
    FSG=0 ! divergence correction; set to 0 for the time being

    ! set up coulomb potentials in reciprocal space, each for another ENCUT and ENCUTSOFT
    CALL LTMP2_SET_GFAC_VKPT_ENCUTS(GRIDHF, LATT_CUR, -KPOINTS_FULL%VKPT(:,NQ), NQ, FSG, POTFAKS, ENCUTS, ENCUTSOFTS)

    GWORK=0
    GWORKSUM=0
    
    DO K1=1,WOCCHF%WDES%NKPTS
      CALL SETWDES(WOCCHF%WDES, WDESK1, K1)
      
      ! k2 = T(k1-q)
      K2 = KPOINT_IN_FULL_GRID(WOCC%WDES%VKPT(:,K1)-KPOINTS_FULL%VKPT(:,NQ), KPOINTS_FULL)
      CALL SETWDES(WTHF%WDES, WDESK2, K2)

#ifndef gammareal
      ! set phase factors 
      CALL PHASER_HF(GRIDHF, LATT_CUR, FAST_AUG_FOCK, WOCC%WDES%VKPT(:,K2)-WOCC%WDES%VKPT(:,K1))
      ! set phase factor e^iK(k1-q)r where K(k1-q) = k2-k1-q
      CALL SETPHASE(WOCC%WDES%VKPT(:,K2)-WOCC%WDES%VKPT(:,K1)+KPOINTS_FULL%VKPT(:,NQ), WGWQ%GRID, CPHASE, LPHASE)
#endif

      DO ISP = 1, ISPIN
        DO NI = 1, WOCC%WDES%NBANDS ! corresponds to sum over i
          CALL W1_COPY(ELEMENT(WOCCHF, WDESK1, NI, ISP), WOCC1)
          CALL W1_COPY(ELEMENT(WTHF,   WDESK2, NI, ISP), WT2)

          CALL FFTWAV_W1(WOCC1)
          CALL FFTWAV_W1(WT2)
          NFFT = NFFT + 2

          ! calculate overlap density 
          ! < u_i,k1 | r > < r | w_i,T(k1-q) >
          CALL FOCK_CHARGE(WT2, WOCC1, GWORK, CRHOLM)

#ifndef gammareal
          IF (LPHASE) CALL APPLY_PHASE_GDEF( GRIDHF, CPHASE(1), GWORK(1), GWORK(1) )
#endif

          GWORKSUM = GWORKSUM + GWORK * WOCCHF%FERWE(NI, K1, ISP) * WDESK1%RSPIN
        ENDDO
      ENDDO
    ENDDO
     
    ! go back to reciprocal space
    CALL FFT3D(GWORKSUM,GRIDHF,-1)
    NFFT = NFFT + 1

    CALLMPI(M_sum_g(WOCCHF%WDES%COMM_KIN, GWORKSUM , SIZE(GWORKSUM)))

    IF (MCALPHA/=0) THEN
      ! with finite size corrections:
      CALL vtutor%bug("internal error in VASP: LAPLACE_TRANSFORMED_MP2 does not support MCALPHA/=0", __FILE__, __LINE__)
    ELSE
      ! sum_G': apply coulomb potentials for each cutoff  
      CALL APPLY_GFAC_2_ENCUTS(GRIDHF, GWORKSUM(1), GWORKSUM(1), POTFAKS(1,1), E_DIRECT(1), ENCUTS(1), N_MP2_ENCUTS)
      E_DIRECT = E_DIRECT * (0.5_q/GRIDHF%NPLWV) * WDESK1%WTKPT  ! divide by grid points
    ENDIF

    DEALLOCATE(GWORK, GWORKSUM, CRHOLM, POTFAKS)
    CALL DELWAV(WOCC1,.TRUE.)
    CALL DELWAV(WT2,.TRUE.)

#ifdef PROFILING
    CALL STOP_PROFILING('CALCULATE_DIRECT')
#endif
  END SUBROUTINE CALCULATE_DIRECT



  
  !****************************** CALCULATE_EXCHANGE ******************************
  !
  ! calculate the MP2 exchange contribution for the current G, q and tau
  !
  ! E_EXCH = sum_k1,k2 sum_ij sun_G' v(G'+T(k1+k2-q))) < phi*_i,k1 | exp(-iG'r) exp(-T(k1+k2-q)r) | w_j,T(k2-q) >* 
  !                                                    < phi*_j,k2 | exp(-iG'r) exp(-T(k1+k2-q)r) | w_i,T(k1-q) >
  ! 
  !********************************************************************************
  SUBROUTINE CALCULATE_EXCHANGE(WOCC, WT, LAST_FILLED, LATT_CUR, E_EXCH, ENCUTS, ENCUTSOFTS, WGWQ, NQ, ISPIN, NFFT, P)

    TYPE (wavespin) :: WOCC           ! occupied states
    TYPE (wavespin) :: WT             ! schmidt decomposed states
    INTEGER :: LAST_FILLED
    TYPE (latt) :: LATT_CUR
    REAL(q) :: E_EXCH(N_MP2_ENCUTS)   ! accumulated exchange energy
    TYPE (wavedes1)    WGWQ           ! wave descriptor of q-point
    INTEGER :: NQ                     ! q-point
    INTEGER :: ISPIN                  ! number of spins
    REAL(q) :: ENCUTS(N_MP2_ENCUTS)   ! energy cutoffs of Coulomb kernel
    REAL(q) :: ENCUTSOFTS(N_MP2_ENCUTS) ! smooth energy cutoffs of Coulomb kernel
    REAL(q) :: NFFT                   ! FFT counter
    TYPE (potcar) :: P(:)
    ! local
    TYPE (wavespin) WOCCHF, WTHF                      ! same as WOCC and WT but with HF wave descriptor
    TYPE (wavefun1) :: WOCC1                          ! occupied state
    TYPE (wavefun1) :: WT1                            ! schmidt state
    INTEGER :: NBLK, N, IBLK, NGLB, NBSTART, NBSTOP   ! variables for blocking 
    TYPE (wavefun1), ALLOCATABLE :: WOCCGATHER(:)     ! blocked occupied states
    TYPE (wavefun1), ALLOCATABLE :: WTGATHER(:)       ! blocked schmidt states
    REAL(q) :: FSG                                    ! singularity correction setting V(G=0) to a finite value
    TYPE(rotation_handle), POINTER :: ROT_HANDLE
    GDEF, ALLOCATABLE :: GWORK1(:,:)                  ! buffer for < phi*_i,k1 | exp(-iG'r | w_j,T(k2-q) >*
    GDEF, ALLOCATABLE :: GWORK2(:,:)                  ! buffer for < phi*_j,k2 | exp(-iG'r | w_i,T(k1-q) >
    GDEF, ALLOCATABLE :: GWORK_SUM(:)  
    GDEF, ALLOCATABLE :: CRHOLM(:,:)                  ! augmentation occupancy matrix
    REAL(q), ALLOCATABLE :: POTFAKS(:,:)                ! potential in reciprocal space
    TYPE (wavedes1), TARGET :: WDESK1, WDESK2         ! k1, k2
    TYPE (wavedes1), TARGET :: WDESQ1, WDESQ2         ! q1 = k1-q , q2 = k2-q
    !TYPE (wavedes1), TARGET :: WDESQ3                ! q3 = q-k1-k2 
    TYPE (wavedes1), TARGET :: WDESK1_IRZ, WDESQ1_IRZ ! k1, q1 in the irreducible BZ
    INTEGER :: ISP, ISP_IRZ                           ! spin loop
    INTEGER :: K1, K2, Q1, Q2, Q3                     ! k-point loops
    INTEGER :: NI, NJ_GLOBAL, NJ                      ! band loop
    COMPLEX(q) :: CPHASE1(WGWQ%GRID%MPLWV)            ! phase shift
    COMPLEX(q) :: CPHASE2(WGWQ%GRID%MPLWV)            ! phase shift
    REAL(q) :: EXX(N_MP2_ENCUTS)                   ! contr. to exchange energy
    LOGICAL :: LSHIFT, LPHASE1, LPHASE2
#ifdef gammareal
    INTEGER :: NGLB_TRIMMED, NI_GLOBAL                            
#endif

    IF (AEXX==0) THEN
      RETURN
    ENDIF

#ifdef PROFILING
    CALL START_PROFILING('CALCULATE_EXCHANGE')
#endif

    ! use temporarily another WDES
    WOCCHF = WOCC
    WOCCHF%WDES => WDES_FOCK
    WTHF = WT
    WTHF%WDES => WDES_FOCK

    CALL CHECK_FULL_KPOINTS

    NULLIFY(ROT_HANDLE)

    ! allocate memory, we have to do the action on nstripn bands
    NBLK = ((WOCCHF%WDES%NSIM*2 + WOCCHF%WDES%NB_PAR - 1) / WOCCHF%WDES%NB_PAR) * WOCCHF%WDES%NB_PAR
    !IF (NBLOCK_FOCK > 0) NBLK = MIN(NBLK, NBLOCK_FOCK)    ! deactivated since NBLOCK_FOCK is to small if
                                                           ! OMP parallelization is activated
    NBLK = MIN(NBLK, 64) ! instead manually limit to 64

    ALLOCATE( GWORK1(m_ GRIDHF%MPLWV, NBLK),                       &
              GWORK2(m_ GRIDHF%MPLWV, NBLK),                       &
              GWORK_SUM(m_ GRIDHF%MPLWV),                       &
              CRHOLM(AUG_DES%NPROD * WOCCHF%WDES%NRSPINORS, NBLK), &
              POTFAKS(GRIDHF%MPLWV, N_MP2_ENCUTS) )

    ! allocate memory for wave functions
    CALL SETWDES(WOCCHF%WDES, WDESK1, 0)
    CALL NEWWAV(WOCC1, WDESK1, .TRUE.)

    CALL SETWDES(WTHF%WDES, WDESQ1, 0)
    CALL NEWWAV(WT1, WDESQ1, .TRUE.)

    CALL SETWDES(WOCCHF%WDES, WDESK2, 0)
    CALL SETWDES(WTHF%WDES, WDESQ2, 0)
    ALLOCATE(WOCCGATHER(NBLK))
    ALLOCATE(WTGATHER(NBLK))
    DO N = 1, NBLK
      CALL NEWWAV(WOCCGATHER(N) , WDESK2, .TRUE.)
      CALL NEWWAV(WTGATHER(N) , WDESQ2, .TRUE.)
    ENDDO

    ! average electrostatic potential 
    !FSG = FSG_STORE(1) 
    FSG = 0 ! set to zero for the time beeing

    E_EXCH = 0

    spin: &
    DO ISP = 1, ISPIN

      block: &
      DO IBLK = 1, LAST_FILLED, NBLK
        NGLB = MIN(LAST_FILLED-IBLK+1, NBLK)
        ! global index of the first and last band in this block
        NBSTART = IBLK
        NBSTOP  = NBSTART + NGLB - 1

        ! irgendwo muss noch dieses ganze IF Zeug fuer die k-Punkte
        ! und das IF mit KPOINTS_FULL%NEQUIV      
        ! wo und warum WKPT und RSPIN?

            
        kpoint2: &
        DO K2 = 1, WOCCHF%WDES%NKPTS
          CALL SETWDES(WOCCHF%WDES, WDESK2, K2)
          Q2 = KPOINT_IN_FULL_GRID(WOCC%WDES%VKPT(:,K2) - KPOINTS_FULL%VKPT(:,NQ), KPOINTS_FULL)
          CALL SETWDES(WTHF%WDES, WDESQ2, Q2)

          ! get WOCCGATHER and WTGATHER
          CALL W1_GATHER_GLB(WOCCHF, NBSTART, NBSTOP, ISP, WOCCGATHER)
          CALL W1_GATHER_GLB(WTHF, NBSTART, NBSTOP, ISP, WTGATHER)
          NFFT = NFFT + 2*NGLB

          kpoint1: &
          DO K1 = 1, WOCCHF%WDES%NKPTS
            CALL SETWDES(WOCCHF%WDES, WDESK1, K1)
            CALL SETWDES(WOCCHF%WDES, WDESK1_IRZ, KPOINTS_FULL%NEQUIV(K1))
            Q1 = KPOINT_IN_FULL_GRID(WOCC%WDES%VKPT(:,K1) - KPOINTS_FULL%VKPT(:,NQ), KPOINTS_FULL)
            CALL SETWDES(WTHF%WDES, WDESQ1, Q1)
            CALL SETWDES(WOCCHF%WDES, WDESQ1_IRZ, KPOINTS_FULL%NEQUIV(Q1))
            Q3 = KPOINT_IN_FULL_GRID(WOCC%WDES%VKPT(:,K1) + WOCC%WDES%VKPT(:,K2) - KPOINTS_FULL%VKPT(:,NQ), KPOINTS_FULL)

            ! set up coulomb potentials in reciprocal space, each for another ENCUT and ENCUTSOFT
            CALL LTMP2_SET_GFAC_VKPT_ENCUTS(GRIDHF, LATT_CUR, KPOINTS_FULL%VKPT(:,Q3), Q3, FSG, POTFAKS, ENCUTS, ENCUTSOFTS)

            GWORK_SUM(:) = 0

            jband:&  ! corresponds to sum over j
            DO NJ = 1, WOCC%WDES%NBANDS
              ! get global band index 
              NJ_GLOBAL = (NJ-1) * WOCCHF%WDES%NB_PAR + WOCCHF%WDES%NB_LOW
#ifdef gammareal
              IF(NJ_GLOBAL < NBSTART) CYCLE jband
              NGLB_TRIMMED = MIN(NGLB, NJ_GLOBAL - NBSTART + 1) ! reset NGLB in the gamma-only version
#endif    
              IF(NJ_GLOBAL > LAST_FILLED) CYCLE jband
              IF(NJ_GLOBAL < NBANDSGWLOW_FOCK) CYCLE jband

              ! get WOCC1
              ! symmetry must be considered if the wavefunctions for this k-point are not stored in W
              IF (K1 <= WOCCHF%WDES%NKPTS) THEN
                CALL W1_COPY(ELEMENT(WOCCHF, WDESK1, NJ, ISP), WOCC1)
                CALL FFTWAV_W1(WOCC1)
                NFFT = NFFT + 1
#ifndef gammareal
              ELSE
                LSHIFT=.FALSE.
                IF ((ABS(KPOINTS_FULL%TRANS(1, K1)) > TINY) .OR. &
                    (ABS(KPOINTS_FULL%TRANS(2, K1)) > TINY) .OR. &
                    (ABS(KPOINTS_FULL%TRANS(3, K1)) > TINY)) LSHIFT=.TRUE.
                ISP_IRZ = ISP
                IF (KPOINTS_FULL%SPINFLIP(K1) == 1) THEN
                  ISP_IRZ = 3-ISP
                ENDIF
                CALL W1_ROTATE_AND_FFT(WOCC1, ELEMENT(WOCCHF, WDESK1_IRZ, NJ, ISP_IRZ), ROT_HANDLE, P, LATT_CUR, LSHIFT)
                NFFT = NFFT + 1
#endif
              ENDIF
              !complex conjugate 
              WOCC1%CR = CONJG(WOCC1%CR)
              WOCC1%CPROJ = GCONJG(WOCC1%CPROJ)

              ! get WT1
              ! symmetry must be considered if the wavefunctions for this k-point are not stored in WT
              IF (Q1 <= WTHF%WDES%NKPTS) THEN
                CALL W1_COPY(ELEMENT(WTHF, WDESQ1, NJ, ISP), WT1)
                CALL FFTWAV_W1(WT1)
                NFFT = NFFT + 1
#ifndef gammareal
              ELSE
                LSHIFT=.FALSE.
                IF ((ABS(KPOINTS_FULL%TRANS(1, Q1)) > TINY) .OR. &
                    (ABS(KPOINTS_FULL%TRANS(2, Q1)) > TINY) .OR. &
                    (ABS(KPOINTS_FULL%TRANS(3, Q1)) > TINY)) LSHIFT=.TRUE.
                ISP_IRZ = ISP
                IF (KPOINTS_FULL%SPINFLIP(Q1) == 1) THEN
                  ISP_IRZ = 3-ISP
                ENDIF
                CALL W1_ROTATE_AND_FFT(WT1, ELEMENT(WTHF, WDESQ1_IRZ, NJ, ISP_IRZ), ROT_HANDLE, P, LATT_CUR, LSHIFT)
                NFFT = NFFT + 1
#endif
              ENDIF
              !complex conjugate 
              WT1%CR = CONJG(WT1%CR)
              WT1%CPROJ = GCONJG(WT1%CPROJ)

              GWORK1 = 0
              GWORK2 = 0 

              ! calculate overlaps in real space (most expensive step when tested with C8 and C16)
#ifdef gammareal
              ! < r | u_j > < r | w_i >
              CALL FOCK_CHARGE_MU(WOCCGATHER(1:NGLB_TRIMMED), WT1, GWORK1(:,1:NGLB_TRIMMED), CRHOLM) ! quartic scaling 
              ! < r | u_i > < r | w_j >
              CALL FOCK_CHARGE_MU(WTGATHER(1:NGLB_TRIMMED), WOCC1, GWORK2(:,1:NGLB_TRIMMED), CRHOLM) ! quartic scaling 
#else
              ! set correct phase
              CALL PHASER_HF(GRIDHF, LATT_CUR, FAST_AUG_FOCK, WOCC%WDES%VKPT(:,K2) + KPOINTS_FULL%VKPT(:,Q1))
              CALL SETPHASE(KPOINTS_FULL%VKPT(:,Q1) + WOCC%WDES%VKPT(:,K2) - KPOINTS_FULL%VKPT(:,Q3), WGWQ%GRID, CPHASE1, LPHASE1)
              ! < r | u_j,k2 > < r | w_i,T(k1-q) >
              CALL FOCK_CHARGE_MU(WOCCGATHER(1:NGLB), WT1, GWORK1(:,1:NGLB), CRHOLM) ! quartic scaling

              ! set correct phase
              CALL PHASER_HF(GRIDHF, LATT_CUR, FAST_AUG_FOCK, WOCC%WDES%VKPT(:,K1) + KPOINTS_FULL%VKPT(:,Q2))
              CALL SETPHASE(KPOINTS_FULL%VKPT(:,Q2) + WOCC%WDES%VKPT(:,K1) - KPOINTS_FULL%VKPT(:,Q3), WGWQ%GRID, CPHASE2, LPHASE2)
              ! < r | w_j,k1 > < r | u_i,T(k2-q) > 
              CALL FOCK_CHARGE_MU(WTGATHER(1:NGLB), WOCC1, GWORK2(:,1:NGLB), CRHOLM) ! quartic scaling
#endif

              ! sum over i
#ifdef gammareal
              DO NI = 1, NGLB_TRIMMED
                ! FFT back to reciprocal space 
                CALL FFT3D(GWORK1(1,NI), GRIDHF, -1) ! quartic scaling
                CALL FFT3D(GWORK2(1,NI), GRIDHF, -1) ! quartic scaling
                NFFT = NFFT + 2
                
                NI_GLOBAL = NBSTART + NI - 1 ! gloabl i index
                IF(NI_GLOBAL /= NJ_GLOBAL) THEN
                  ! multiply overlaps and add to GWORK_SUM ! quartic scaling
                  CALL ADD_OVERLAPS(GRIDHF, GWORK1(1,NI), GWORK2(1,NI), GWORK_SUM(1), 2.0_q)
                ELSE
                  ! multiply overlaps and add to GWORK_SUM ! quartic scaling
                  CALL ADD_OVERLAPS(GRIDHF, GWORK1(1,NI), GWORK2(1,NI), GWORK_SUM(1), 1.0_q)
                ENDIF
              ENDDO ! NI loop
#else
              DO NI = 1, NGLB
                ! apply phase factors
                IF (LPHASE1) CALL APPLY_PHASE_GDEF(GRIDHF, CPHASE1(1), GWORK1(1,NI), GWORK1(1,NI))
                IF (LPHASE2) CALL APPLY_PHASE_GDEF(GRIDHF, CPHASE2(1), GWORK2(1,NI), GWORK2(1,NI))

                ! FFT back to reciprocal space 
                CALL FFT3D(GWORK1(1,NI), GRIDHF, -1) ! quartic scaling
                CALL FFT3D(GWORK2(1,NI), GRIDHF, -1) ! quartic scaling
                NFFT = NFFT + 2
                
                ! multiply overlaps and add to GWORK_SUM ! quartic scaling
                CALL ADD_OVERLAPS(GRIDHF, GWORK1(1,NI), GWORK2(1,NI), GWORK_SUM(1), 1.0_q)
              ENDDO ! NI loop
#endif
            ENDDO jband

            ! apply coulomb potential for each cutoff
            CALL APPLY_GFAC_ENCUTS(GRIDHF, GWORK_SUM(1), POTFAKS(1,1), EXX(1), ENCUTS(1), N_MP2_ENCUTS)
            E_EXCH = E_EXCH - EXX * (0.5_q/GRIDHF%NPLWV) * WDESK1%WTKPT  ! divide by grid points
          
          ENDDO kpoint1
        ENDDO kpoint2
      ENDDO block
    ENDDO spin

    CALLMPI(M_sum_d(WOCCHF%WDES%COMM_KIN, E_EXCH, N_MP2_ENCUTS))

    DEALLOCATE(GWORK1, GWORK2, GWORK_SUM, CRHOLM, POTFAKS) 
    CALL DEALLOCATE_ROT_HANDLE(ROT_HANDLE)
    CALL DELWAV(WOCC1, .TRUE.)
    CALL DELWAV(WT1, .TRUE.)
    DO N = 1, NBLK
      CALL DELWAV(WOCCGATHER(N), .TRUE.)
      CALL DELWAV(WTGATHER(N), .TRUE.)
    ENDDO
    DEALLOCATE(WOCCGATHER, WTGATHER)

#ifdef PROFILING
    CALL STOP_PROFILING('CALCULATE_EXCHANGE')
#endif

  END SUBROUTINE CALCULATE_EXCHANGE




  !****************************** LTMP2_SET_GFAC_VKPT_ENCUTS *****************************
  !
  ! set up the Coulomb kernels (POTFAK) for each cutoff (ENCUTS)
  !
  ! 4 pi e^2 / (G+q)^2
  !
  ! multiplied by a smooth cosinus cutoff between ENCUTSOFT and ENCUT
  !
  ! additionally set the zero grid point value to the convergence
  ! correction FSG if G+q = 0 
  !
  ! for maximal performance the routine includes a weighting factor 1/NRPLWV
  ! and furthermore includes a q-point weighting factor
  ! KPOINTS_FULL%WTKPT(NQ)*NKREDX*NKREDY*NKREDZ
  !
  !*****************************************************************************
  SUBROUTINE LTMP2_SET_GFAC_VKPT_ENCUTS(GRID, LATT_CUR, VKPT, NQ, FSG, POTFAKS, ENCUTS, ENCUTSOFTS)
    USE constant
    USE mgrid
    USE nonl_high

    TYPE (grid_3d) GRID
    TYPE (latt) LATT_CUR
    REAL(q) :: VKPT(3)            ! coordinate of k-point in reciprocal lattice units
    INTEGER :: NQ                 ! k-point at which contribution to exchange kernel is evaluated
    REAL(q) :: FSG
    REAL(q) :: POTFAKS(GRID%MPLWV, N_MP2_ENCUTS) ! the resulting Coulomb kernel
    REAL(q) :: ENCUTS(N_MP2_ENCUTS)
    REAL(q) :: ENCUTSOFTS(N_MP2_ENCUTS)
    ! local
    INTEGER :: I, NI, NC, N1, N2, N3
    REAL(q) :: DQX, DQY, DQZ, GX, GY, GZ, GSQU, SCALE, EN

    CALL CHECK_FULL_KPOINTS

    ! all k-points in the full grid have equal weight equal to first one
    SCALE = KPOINTS_FULL%WTKPT(NQ) * NKREDX * NKREDY * NKREDZ * EDEPS / LATT_CUR%OMEGA / TPI**2

    DQX = VKPT(1) * LATT_CUR%B(1,1) + VKPT(2) * LATT_CUR%B(1,2) + VKPT(3) * LATT_CUR%B(1,3)
    DQY = VKPT(1) * LATT_CUR%B(2,1) + VKPT(2) * LATT_CUR%B(2,2) + VKPT(3) * LATT_CUR%B(2,3)
    DQZ = VKPT(1) * LATT_CUR%B(3,1) + VKPT(2) * LATT_CUR%B(3,2) + VKPT(3) * LATT_CUR%B(3,3)

    NI = 0
    col: DO NC = 1, GRID%RC%NCOL
      N2 = GRID%RC%I2(NC)
      N3 = GRID%RC%I3(NC)
      row: DO N1 = 1, GRID%RC%NROW
        NI = NI + 1
        GX = (GRID%LPCTX(N1) * LATT_CUR%B(1,1) + GRID%LPCTY(N2) * LATT_CUR%B(1,2) + GRID%LPCTZ(N3) * LATT_CUR%B(1,3))
        GY = (GRID%LPCTX(N1) * LATT_CUR%B(2,1) + GRID%LPCTY(N2) * LATT_CUR%B(2,2) + GRID%LPCTZ(N3) * LATT_CUR%B(2,3))
        GZ = (GRID%LPCTX(N1) * LATT_CUR%B(3,1) + GRID%LPCTY(N2) * LATT_CUR%B(3,2) + GRID%LPCTZ(N3) * LATT_CUR%B(3,3))
        GSQU = (DQX+GX)**2 + (DQY+GY)**2 + (DQZ+GZ)**2

        ! calculate corresponding energy to GSQU
        EN = HSQDTM * GSQU * TPI**2 

        ! set FSG if G+q = 0
        IF (GSQU < KPOINTS_FULL%VKPTMINDIST2 / 40) THEN
          POTFAKS(NI,:) = FSG
          CYCLE
        ENDIF

        ! loop through all cutoffs
        DO I = N_MP2_ENCUTS, 1, -1
          ! check if this GSQU is below ENCUTSOFT for all cutoffs (if so treat all POTFAKSS equally)
          IF (EN <= ENCUTSOFTS(I)) THEN
            POTFAKS(NI, 1:I) = SCALE / GSQU ! set 4 pi e^2 / (G+q)^2 if we are below the smooth cutoff threshold
            EXIT
          ELSE
            IF (EN >= ENCUTS(I)) THEN
              POTFAKS(NI, I) = 0 ! set 0 if we are above the cutoff 
            ELSE
              ! apply smooth cutoff since we are between ENCUTSOFT and ENCUT 
              POTFAKS(NI, I) = ( 1 + COS((EN-ENCUTSOFTS(I))/(ENCUTS(I)-ENCUTSOFTS(I))*PI) ) * SCALE / ( 2 * GSQU )
            ENDIF
          ENDIF
        ENDDO

      ENDDO row
    ENDDO col

    DO I = 1, N_MP2_ENCUTS
      CALL SETUNB_REAL(POTFAKS(:,I), GRID)
    ENDDO

    IF (L_MODEL_HF) THEN
      DO NC = 1, GRID%RC%NP
        POTFAKS(NC,:) = 1.0_q * POTFAKS(NC,:) / GRID%NPLWV 
      ENDDO
    ELSE
      DO NC=1,GRID%RC%NP
        POTFAKS(NC,:) = POTFAKS(NC,:) * AEXX / GRID%NPLWV
      ENDDO
    ENDIF

  END SUBROUTINE LTMP2_SET_GFAC_VKPT_ENCUTS




  !****************************** LIN_REG_CONVERGE_LTMP2 *****************************
  !
  ! Perform a linear regression of the input data: y = M*x + B, where 
  ! x = ENCUTS^(-3/2) and y = ENERGIES. 
  !
  ! This is used to extrapolate ENCUT -> infinity
  ! Because: x(ENCUT -> infinity) = 0 => y -> b = converged energy
  !
  ! see: http://dx.doi.org/10.1103/PhysRevB.90.075125
  !
  !***********************************************************************************
  SUBROUTINE LIN_REG_CONVERGE_LTMP2(ENERGIES, ENCUTS, CONVERGED_ENERGY)
    REAL(q) :: ENERGIES(:)         ! calculated energies, each for a given ENCUT
    REAL(q) :: ENCUTS(:)           ! all ENCUTs
    REAL(q) :: CONVERGED_ENERGY   ! the resulting converged energy
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

    CONVERGED_ENERGY = B

  END SUBROUTINE LIN_REG_CONVERGE_LTMP2

END MODULE lt_mp2
#else
MODULE lt_mp2
  CONTAINS
  SUBROUTINE lt_mp2_dummy
      WRITE(*,*)'Im a DEC compiler so I need this line'
  END SUBROUTINE
END MODULE lt_mp2
#endif



  !****************************** APPLY_GFAC_2_ENCUTS ****************************** 
  !
  !  This subroutine multiplies the overlap densities CWORK and CONJG(CWORKG) by the 
  !  potential kernel POTFAKS. This is done for different cutoffs stored in ENCUTS.
  !
  !*********************************************************************************
  SUBROUTINE APPLY_GFAC_2_ENCUTS(GRID, CWORK, CWORKG, POTFAKS, ENERGIES, ENCUTS, N_MP2_ENCUTS)
    USE mgrid
    USE tutor, ONLY: vtutor
    IMPLICIT NONE

    TYPE (grid_3d) :: GRID
    COMPLEX(q)  :: CWORK(GRID%MPLWV)
    COMPLEX(q)  :: CWORKG(GRID%MPLWV)
    REAL(q)     :: POTFAKS(GRID%MPLWV, N_MP2_ENCUTS)
    REAL(q)  :: ENERGIES(N_MP2_ENCUTS)  ! results of the routine
    REAL(q) :: ENCUTS(N_MP2_ENCUTS)
    INTEGER :: N_MP2_ENCUTS
    ! LOCAL
    INTEGER NP, I

    ENERGIES=0
    IF (GRID%REAL2CPLX) THEN
      IF (.NOT. ASSOCIATED(GRID%FFTWEIGHT)) THEN
        CALL vtutor%bug("internal error in APPLY_GFAC_2_ENCUTS: FFTWEIGHT is not associated", __FILE__, __LINE__)
      ENDIF
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(NP) &
!$OMP REDUCTION(+:ENERGIES)
      DO NP=1, GRID%RC%NP
        ! loop through all cutoffs
        DO I = N_MP2_ENCUTS, 1, -1
          ! check if POTFAKSS are equal for different cutoffs (if so treat all POTFAKSS equally)
          IF ( ABS(POTFAKS(NP,I) - POTFAKS(NP,1)) < 1.E-8_q ) THEN
            ENERGIES(1:I) = ENERGIES(1:I) + POTFAKS(NP,1) * CWORK(NP) * CONJG(CWORKG(NP)) * GRID%FFTWEIGHT(NP)
            EXIT
          ! if POTFAKSS differ treat POTFAKS(I) individually
          ELSE
            ENERGIES(I) = ENERGIES(I) + POTFAKS(NP,I) * CWORK(NP) * CONJG(CWORKG(NP)) * GRID%FFTWEIGHT(NP) 
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
    ELSE
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(NP) &
!$OMP REDUCTION(+:ENERGIES)
      DO NP=1, GRID%RC%NP
        ! loop through all cutoffs
        DO I = N_MP2_ENCUTS, 1, -1
          ! check if POTFAKSS are equal for different cutoffs (if so treat all POTFAKSS equally)
          IF ( ABS(POTFAKS(NP,I) - POTFAKS(NP,1)) < 1.E-8_q ) THEN
            ENERGIES(1:I) = ENERGIES(1:I) + POTFAKS(NP,1) * CWORK(NP) * CONJG(CWORKG(NP)) 
            EXIT
          ! if POTFAKSS differ treat POTFAKS(I) individually
          ELSE
            ENERGIES(I) = ENERGIES(I) + POTFAKS(NP,I) * CWORK(NP) * CONJG(CWORKG(NP)) 
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
    ENDIF

  END SUBROUTINE APPLY_GFAC_2_ENCUTS




  !****************************** APPLY_GFAC_ENCUTS ****************************** 
  !
  !  This subroutine multiplies the overlap densitiy CWORK with the potential 
  !  kernel POTFAKS. This is done for different cutoffs stored in ENCUTS.
  !
  !*********************************************************************************
  SUBROUTINE APPLY_GFAC_ENCUTS(GRID, CWORK, POTFAKS, ENERGIES, ENCUTS, N_MP2_ENCUTS)
    USE mgrid
    USE tutor, ONLY: vtutor
    IMPLICIT NONE

    TYPE (grid_3d) :: GRID
    COMPLEX(q)  :: CWORK(GRID%MPLWV)
    REAL(q) :: POTFAKS(GRID%MPLWV, N_MP2_ENCUTS)
    REAL(q) :: ENERGIES(N_MP2_ENCUTS)  ! results of the routine
    REAL(q) :: ENCUTS(N_MP2_ENCUTS)
    INTEGER :: N_MP2_ENCUTS
    ! LOCAL
    INTEGER NP, I

    ENERGIES=0
    IF (GRID%REAL2CPLX) THEN
      IF (.NOT. ASSOCIATED(GRID%FFTWEIGHT)) THEN
        CALL vtutor%bug("internal error in APPLY_GFAC_ENCUTS: FFTWEIGHT is not associated", __FILE__, __LINE__)
      ENDIF
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(NP) &
!$OMP REDUCTION(+:ENERGIES)
      DO NP=1, GRID%RC%NP
        ! loop through all cutoffs
        DO I = N_MP2_ENCUTS, 1, -1
          ! check if POTFAKSS are equal for different cutoffs (if so treat all POTFAKSS equally)
          IF ( ABS(POTFAKS(NP,I) - POTFAKS(NP,1)) < 1.E-8_q ) THEN
            ENERGIES(1:I) = ENERGIES(1:I) + POTFAKS(NP,1) * CWORK(NP) * GRID%FFTWEIGHT(NP)
            EXIT
          ! if POTFAKSS differ treat POTFAKS(I) individually
          ELSE
            ENERGIES(I) = ENERGIES(I) + POTFAKS(NP,I) * CWORK(NP) * GRID%FFTWEIGHT(NP) 
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
    ELSE
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(NP) &
!$OMP REDUCTION(+:ENERGIES)
      DO NP=1, GRID%RC%NP
        ! loop through all cutoffs
        DO I = N_MP2_ENCUTS, 1, -1
          ! check if POTFAKSS are equal for different cutoffs (if so treat all POTFAKSS equally)
          IF ( ABS(POTFAKS(NP,I) - POTFAKS(NP,1)) < 1.E-8_q ) THEN
            ENERGIES(1:I) = ENERGIES(1:I) + POTFAKS(NP,1) * CWORK(NP)  
            EXIT
          ! if POTFAKSS differ treat POTFAKS(I) individually
          ELSE
            ENERGIES(I) = ENERGIES(I) + POTFAKS(NP,I) * CWORK(NP)  
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
    ENDIF

  END SUBROUTINE APPLY_GFAC_ENCUTS




  !****************************** ADD_OVERLAPS ****************************** 
  !
  ! Multiply two overlaps CWORK1*CONJG(CWORK2)*FACTOR and add it to CWORK_SUM
  ! This function is used to overcome a GDEF/COMPLEX problem
  !
  !*********************************************************************************
  SUBROUTINE ADD_OVERLAPS(GRID, CWORK1, CWORK2, CWORK_SUM, FACTOR)
    USE mgrid
    IMPLICIT NONE

    TYPE (grid_3d) :: GRID
    ! complex arrays of size GRID%MPLWV = real arrays of size 2*GRID%MPLWV
    REAL(q) :: CWORK1(2*GRID%MPLWV) 
    REAL(q) :: CWORK2(2*GRID%MPLWV)
    REAL(q) :: CWORK_SUM(2*GRID%MPLWV)
    REAL(q) :: FACTOR

    ! calculate CWORK1*CONJG(CWORK2)*FACTOR and ignore the imaginary part of the result
    ! note that Re(z1*z2^*) = Re(z1)*Re(z2) + Im(z1)*Im(z2)
    CWORK_SUM(1::2) = CWORK_SUM(1::2) + FACTOR * ( CWORK1(1::2)*CWORK2(1::2) + CWORK1(2::2)*CWORK2(2::2) ) 

  END SUBROUTINE ADD_OVERLAPS





!   SUBROUTINE MY_APPLY_GFAC_EXCHANGE_2(GRID, CWORK_SUM, CWORK, CWORKG, POTFAK, EXCHANGE, FACTOR)
!
!     USE prec
!     USE mgrid
!     IMPLICIT NONE
!
!     TYPE (grid_3d) GRID
!     INTEGER NP
!     REAL(q)     :: POTFAK(GRID%MPLWV)
!     !COMPLEX(q)  :: CWORK_SUM(GRID%MPLWV)
!     !COMPLEX(q)  :: CWORK(GRID%MPLWV)
!     !COMPLEX(q)  :: CWORKG(GRID%MPLWV)
!     REAL(q)  :: CWORK_SUM(2*GRID%MPLWV)
!     REAL(q)  :: CWORK(2*GRID%MPLWV)
!     REAL(q)  :: CWORKG(2*GRID%MPLWV)
!     COMPLEX(q)  :: EXCHANGE
!     REAL(q)     :: FACTOR
!     ! multiply by 4 pi e^2/G^2 and divide by # of gridpoints to obtain potential
!     EXCHANGE=0
!     IF (GRID%REAL2CPLX) THEN
!        IF (.NOT. ASSOCIATED(GRID%FFTWEIGHT)) THEN
!           WRITE(0,*) 'internal error in APPLY_GFAC_EXCHANGE: FFTWEIGHT is not associated'
!           STOP
!        ENDIF
!!$OMP PARALLEL DO DEFAULT(SHARED) &
!!$OMP PRIVATE(NP) &
!!$OMP REDUCTION(+:CWORK_SUM)
!        DO NP=1,2*GRID%RC%NP,2
!           !EXCHANGE=EXCHANGE+POTFAK(NP)*CWORK(NP)*CONJG(CWORKG(NP))*GRID%FFTWEIGHT(NP)
!           !CWORK_SUM(NP)=CWORK_SUM(NP) + FACTOR*(CWORK(NP)*CONJG(CWORKG(NP)))
!           !CWORK_SUM(NP)=CWORK_SUM(NP) + CWORK(NP)*CONJG(CWORKG(NP))
!           CWORK_SUM(NP)=CWORK_SUM(NP) + FACTOR* (CWORK(NP)*CWORKG(NP) + CWORK(NP+1)*CWORKG(NP+1) )
!        ENDDO
!        !CWORK_SUM(1::2) = CWORK_SUM(1::2) + FACTOR * ( CWORK(1::2)*CWORKG(1::2) + CWORK(2::2)*CWORKG(2::2) )
!!$OMP END PARALLEL DO
!     ELSE
!!$OMP PARALLEL DO DEFAULT(SHARED) &
!!$OMP PRIVATE(NP) &
!!$OMP REDUCTION(+:CWORK_SUM)
!        DO NP=1,2*GRID%RC%NP,2
!           !EXCHANGE=EXCHANGE+POTFAK(NP)*CWORK(NP)*CONJG(CWORKG(NP))
!           !CWORK_SUM(NP)=CWORK_SUM(NP) + FACTOR*(CWORK(NP)*CONJG(CWORKG(NP)))
!           !CWORK_SUM(NP)=CWORK_SUM(NP) + CWORK(NP)*CONJG(CWORKG(NP))
!           CWORK_SUM(NP)=CWORK_SUM(NP) + FACTOR* (CWORK(NP)*CWORKG(NP) + CWORK(NP+1)*CWORKG(NP+1) )
!        ENDDO
!        !CWORK_SUM(1::2) = CWORK_SUM(1::2) + FACTOR * ( CWORK(1::2)*CWORKG(1::2) + CWORK(2::2)*CWORKG(2::2) )
!!$OMP END PARALLEL DO
!     ENDIF
!   END SUBROUTINE MY_APPLY_GFAC_EXCHANGE_2
