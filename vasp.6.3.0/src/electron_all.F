#include "symbol.inc"
!**********************************************************************
! RCS:  $Id: electron.F,v 1.12 2003/06/27 13:22:15 kresse Exp kresse $
!
! subroutine for performing electronic minimization in VASP
! this version optimises all bands simultaneously using
! a conjugate gradient, or a damped molecular dynamics algorithm
!
!**********************************************************************

      SUBROUTINE ELMIN_ALL( &
          HAMILTONIAN,KINEDEN, &
          P,WDES,NONLR_S,NONL_S,W,W_F,W_G,LATT_CUR,LATT_INI, &
          T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI,CHF,CHAM,ECONV,XCSIF, &
          NSTEP,LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI,LDIMP,LMDIMP)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE charge
      USE pseudo
      USE lattice
      USE steep
      USE us
      USE pot
      USE force
      USE fileio
      USE nonl_high
      USE rmm_diis
      USE david
      USE ini
      USE ebs
      USE rot
      USE wave_high
      USE choleski
      USE mwavpre
      USE mwavpre_noio
      USE msphpro
      USE broyden
      USE msymmetry
      USE subrot
      USE melf
      USE base
      USE mpimy
      USE mgrid
      USE mkpoints
      USE constant
      USE setexm
      USE poscar
      USE wave
      USE hamil
      USE pawm
      USE cl
      USE vaspxml
      USE mdipol
      USE Constrained_M_modular
      USE LDAPLUSU_MODULE
      USE ini
      USE fock
      USE gw_model
      USE hamil_struct_def
      USE meta
      USE pead
      USE electron_common, ONLY: checkAbort, testNumberOfStep
      USE tutor, ONLY: vtutor, isAlert, ALGO_A_ISMEAR, ModelGW
! solvation__
      USE solvation
! solvation__
      USE fock_ace, ONLY : LFOCK_ACE,FOCK_ACE_ALLOCW,WACE,FOCK_ACE_CONSTRUCT,FOCK_ACE_DEALLOCW
#ifdef fock_dblbuf
      USE fock_dbl, ONLY : FOCK_ALL_DBLBUF
#endif
      IMPLICIT NONE
!=======================================================================
!  structures
!=======================================================================
      TYPE (ham_handle)  HAMILTONIAN
      TYPE (tau_handle)  KINEDEN
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      TYPE (wavedes)     WDES
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (wavespin)    W          ! wavefunction
      TYPE (wavespin)    W_F        ! wavefunction for all bands simultaneous
      TYPE (wavespin)    W_G        ! same as above
      TYPE (latt)        LATT_CUR
      TYPE (dynamics)    DYN
      TYPE (info_struct) INFO
      TYPE (in_struct)   IO
      TYPE (mixing)      MIX
      TYPE (kpoints_struct) KPOINTS
      TYPE (symmetry)    SYMM
      TYPE (grid_3d)     GRID       ! grid for wavefunctions
      TYPE (grid_3d)     GRID_SOFT  ! grid for soft chargedensity
      TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
      TYPE (grid_3d)     GRIDUS     ! temporary grid in us.F
      TYPE (grid_3d)     GRIDB      ! Broyden grid
      TYPE (transit)     B_TO_C     ! index table between GRIDB and GRIDC
      TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
      TYPE (transit)     SOFT_TO_C  ! index table between GRID_SOFT and GRIDC
      TYPE (energy)      E,E2
      TYPE (latt)        LATT_INI
     
      INTEGER NSTEP,LMDIM,IRDMAX,NEDOS
      REAL(q) TOTEN,EFERMI

      COMPLEX(q)  CHTOT(GRIDC%MPLWV,WDES%NCDIJ) ! charge-density in real / reciprocal space
      COMPLEX(q)  CHTOTL(GRIDC%MPLWV,WDES%NCDIJ)! old charge-density
      RGRID       DENCOR(GRIDC%RL%NP)           ! partial core
      COMPLEX(q)  CVTOT(GRIDC%MPLWV,WDES%NCDIJ) ! local potential
      COMPLEX(q)  CSTRF(GRIDC%MPLWV,T_INFO%NTYP)

!   augmentation related quantities
      OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
               CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
               CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
!  paw sphere charge density
      INTEGER N_MIX_PAW
      REAL(q)  RHOLM(N_MIX_PAW,WDES%NCDIJ),RHOLM_LAST(N_MIX_PAW,WDES%NCDIJ)
!  charge-density and potential on soft grid
      COMPLEX(q)  CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
      RGRID       SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ)
!  density of states
      REAL(q)    DOS(NEDOS,WDES%ISPIN),DOSI(NEDOS,WDES%ISPIN)
!  all-band-simultaneous-update arrays
      GDEF       CHF(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN), &
                 CHAM(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)
      REAL(q) :: XCSIF(3,3)

! local 
      REAL(q) :: TOTENL=0
      REAL(q) :: DESUM1,DESUM(INFO%NELM*4),ECONV
      INTEGER :: IONODE, NODE_ME
      INTEGER :: NELM, NSIM, IRET, ITMP, IERROR
!  needed temporary for aspherical GGA calculation 
      OVERLAP,ALLOCATABLE ::  CDIJ_TMP(:,:,:,:)
!  local l-projected wavefunction characters (not really used here)
      REAL(q)    PAR(1,1,1,1,WDES%NCDIJ),DOSPAR(1,1,1,WDES%NCDIJ)

      REAL(q), EXTERNAL :: RHO0
      INTEGER N,ISP,ICONJU,IROT,ICEL,I,II,IRDMAA, &
              IERR,IDUM,IFLAG,ICOUEV,ICOUEV2,NN,NORDER,IERRBR,L,LP, &
              NCOUNT
      REAL(q) BTRIAL,RDUM,RMS,ORT,TOTEN2,RMS2,RMST, &
              WEIGHT,BETATO,DESUM2,RMSC,RMSP
      REAL(q) RHOAUG(WDES%NCDIJ),RHOTOT(WDES%NCDIJ)
      COMPLEX(q) CDUM
      CHARACTER (LEN=1) CHARAC
      CHARACTER (LEN=3) CRET
      LOGICAL LABORT_WITHOUT_CONV
! parameters for FAST_SPHPRO
      INTEGER :: LDIMP,LMDIMP
      REAL(q),ALLOCATABLE:: PAR_DUMMY(:,:,:,:,:)
! dynamical contribution to GW 
      REAL(q) :: E_DFT(WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)

      PROFILING_START('elmin_all')

      E_DFT=W%CELTOT

      IONODE=0
      NODE_ME=0
#ifdef MPI
      IONODE  = WDES%COMM%IONODE
      NODE_ME = WDES%COMM%NODE_ME
#endif

      IF (KPOINTS%ISMEAR <-1) &
         CALL vtutor%write(isAlert, ALGO_A_ISMEAR)

      IF (INFO%IALGO /= 3 .AND. INFO%IALGO /= 4 .AND. INFO%IALGO /= 5 .AND. MODEL_GW/=0) &
         CALL vtutor%write(isAlert, ModelGW)

      NELM=INFO%NELMALL

      ! to make timing more sensefull syncronize now
      CALLMPI( M_barrier( WDES%COMM ))
      CALL START_TIMING("LOOP")

      io_begin
      IF (INFO%LONESW) THEN
         IF (IO%IU0>=0) WRITE(IO%IU0,141)
         WRITE(17,141)
 141     FORMAT('       N       E                     dE             ' &
        &      ,'d eps       ncg     rms          ort')
      ELSE
         IF (IO%IU0>=0) WRITE(IO%IU0,142)
         WRITE(17,142)
 142     FORMAT('       N       E                     dE             ' &
        &      ,'d eps       ncg     rms          rms(c)')
      ENDIF
      io_end

      DESUM1=0
      INFO%LMIX=.FALSE.

 130  FORMAT (5X, //, &
     &'----------------------------------------------------', &
     &'----------------------------------------------------'//)

 140  FORMAT (5X, //, &
     &'--------------------------------------- Iteration ', &
     &I6,'(',I4,')  ---------------------------------------'//)

      DWRITE0 'electronall entered'

      DOS=0
      DOSI=0
      CALL DIPOL_RESET()


!=======================================================================
!  reset some arrays to zero to avoid NaN errors
!=======================================================================
      CHF=0
      W_F%CW=0
      W_G%CW=0
      W_F%FERWE=0
      W_G%FERWE=0
      CHTOTL=CHTOT
      RHOLM_LAST=RHOLM

      W%AUX=1
      ! force an update of local potential 
      IF ( ABS(INFO%NELMDL) >=3 .AND. LHFCALC .AND. .NOT. INFO%LCHCON) INFO%LPOTOK=.FALSE.

#ifdef _OPENACC
      PUSH_ACC_EXEC_ON(.TRUE.)

      ACC_ASYNC_Q=ACC_ASYNC_SYNC

      CALL ACC_COPYIN_TYPED_VAR(WDES)
      CALL ACC_COPYIN_TYPED_VAR(W)
      CALL ACC_COPYIN_TYPED_VAR(W_F)
      CALL ACC_COPYIN_TYPED_VAR(W_G)
      CALL ACC_COPYIN_TYPED_VAR(NONL_S)
      CALL ACC_COPYIN_TYPED_VAR(NONLR_S)

      CALL ACC_COPYIN_TYPED_VAR(GRID)
      CALL ACC_COPYIN_TYPED_VAR(GRIDC)
      CALL ACC_COPYIN_TYPED_VAR(GRID_SOFT)
      CALL ACC_COPYIN_TYPED_VAR(SOFT_TO_C)
      CALL ACC_COPYIN_TYPED_VAR(C_TO_US)
      CALL ACC_COPYIN_TYPED_VAR(P)

!$ACC ENTER DATA COPYIN(SV,CDIJ,CQIJ) __IF_ASYNC__
!$ACC ENTER DATA COPYIN(DENCOR,CSTRF,LATT_CUR) __IF_ASYNC__
!$ACC ENTER DATA COPYIN(HAMILTONIAN) __IF_ASYNC__
!$ACC ENTER DATA COPYIN(HAMILTONIAN%MU) IF(ASSOCIATED(HAMILTONIAN%MU).AND.ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(HAMILTONIAN%AVEC) IF(ASSOCIATED(HAMILTONIAN%AVEC).AND.ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(KINEDEN,KINEDEN%TAUC) IF(ASSOCIATED(KINEDEN%TAUC).AND.ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(KINEDEN,KINEDEN%TAU) IF(ASSOCIATED(KINEDEN%TAU).AND.ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(CHF,CHAM) __IF_ASYNC__
!$ACC WAIT IF(ACC_EXEC_ON)
#endif
!=======================================================================
      electron: DO N=1,NELM

      ICOUEV=0

      CALL XML_TAG("scstep")
!======================================================================
      io_begin
      WRITE(IO%IU6,140) NSTEP,N
      io_end
!=======================================================================
! if recalculation of total local potential is necessary (INFO%LPOTOK=.F.)
! call POTLOK: the subroutine calculates
! ) the hartree potential from the electronic  charge density
! ) the exchange correlation potential
! ) and the total lokal potential
!  in addition all double counting correction and forces are calculated
! &
! call SETDIJ
! calculates the Integral of the depletion charges * local potential
! and sets CDIJ
!=======================================================================
      CALL START_TIMING("G")

      CALL WVREAL(WDES,GRID,W) ! only for gamma some action

      ! calculate <psi_n| rho | psi_m>
      CALL GW_MODEL_SET_RHOMAT(P, W, LATT_CUR, &
          T_INFO, GRIDC, GRIDUS, GRID_SOFT, SOFT_TO_C, C_TO_US, &
          IRDMAX, LMDIM, DENCOR, CHTOT )

      IF (.NOT. INFO%LPOTOK) THEN
         ! in delay phase, only DFT Hamiltonian is applied, so calculate full DFT XC potential
         IF ( ABS(INFO%NELMDL) >=3 .AND. N < ABS(INFO%NELMDL) .AND. LHFCALC) THEN
            ! Set all HF exchange contributions to zero
            CALL PUSH_XC_TYPE(P(1)%LEXCH,1._q,1._q,1._q,1._q,0._q)
            ! Initialize xc tables
            CALL SETUP_LDA_XC(1, -1, -1,IO%IDIOT)
         ENDIF
         CALL POTLOK(GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES, &
                     INFO,P,T_INFO,E,LATT_CUR, &
                     CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF)

         CALL POTLOK_METAGGA(KINEDEN, &
                     GRID,GRIDC,GRID_SOFT,WDES%COMM_INTER,WDES,INFO,P,T_INFO,E,LATT_CUR, &
                     CHDEN,CHTOT,DENCOR,CVTOT,SV,HAMILTONIAN%MUTOT,HAMILTONIAN%MU,SOFT_TO_C,XCSIF)
         CALL STOP_TIMING("G",IO%IU6,'POTLOK')
         DWRITE0 'potlok is ok'

         CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
                     LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

         CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
            WDES%NCDIJ, LMDIM, CDIJ(1,1,1,1),  RHOLM, CRHODE(1,1,1,1), &
            E,  LMETA=.FALSE., LASPH=INFO%LASPH, LCOREL= .FALSE.  )

         IF (  ABS(INFO%NELMDL) >=3 .AND. N < ABS(INFO%NELMDL) .AND. LHFCALC) THEN
            ! Restore the original situation
            CALL POP_XC_TYPE
            ! reset the xc tables
            IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
               CALL SETUP_LDA_XC(2,-1,-1,IO%IDIOT)
            ELSE
               CALL SETUP_LDA_XC(1,-1,-1,IO%IDIOT)
            ENDIF
         ENDIF
      
         CALL UPDATE_CMBJ(GRIDC,T_INFO,LATT_CUR,IO%IU6)

         ! write LDA+U occupancy matrices
         IF (USELDApU()) CALL LDAPLUSU_PRINTOCC(WDES,T_INFO%NIONS,T_INFO%ITYP,IO%IU6)

         CALL STOP_TIMING("G",IO%IU6,'SETDIJ')
         DWRITE0 'setdij is ok'

         INFO%LPOTOK=.TRUE.
      ENDIF
!=======================================================================
! use
! STEEPEST DESCENT/ CONJUGATE GRADIENT/ DAMPED EQUATION OF MOTION/ 
!  QUICKMIN
! UPDATE ALL BANDS SIMULTANEOUSLY
! the energy is exact and evaluated from the input potential
! this part of the code is influenced by NELMDL:
! during the delay (N<= INFO%NELMDL) the chargedensity is not updated
! and in each step the wavfunctions are diagonalized exactly
! steepest descent approach is used (conjugation seems to be unreliable)
!=======================================================================
      CALL START_TIMING("G")
!-----------------------------------------------------------------------
!  PREDICTOR STEP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  first step and code for delay
!-----------------------------------------------------------------------
      ICONJU=0
      IF (INFO%LSUBROT) THEN
         IF (INFO%LDIAG) THEN
            ! new algorithm
            IROT=7
            ICEL=2
         ELSE
            ! subspace preconditioning
            IROT=4
            ICEL=1
         ENDIF
      ELSE
         IF (INFO%LDIAG) THEN
            ! supersede by new algorithm
            IROT=6
            ICEL=2
         ELSE
            ! no subspace preconditioning
            IROT=3
            ICEL=1
         ENDIF
      ENDIF
      IF ((NINT(INFO%NELECT)==WDES%NB_TOT .AND. WDES%LNONCOLLINEAR) .OR.  & 
          (NINT(INFO%NELECT)==WDES%NB_TOT*2 .AND. .NOT. WDES%LNONCOLLINEAR)) THEN
         IROT=1
         ICEL=0
      ENDIF
      ! pead is applicable to insulators only, and unoccupied bands are "switched off"
      ! using NB_TOTK, hence no subspace rotation applied
      IF (LPEAD_EFIELD_SWITCHED_ON()) THEN
          IROT=1
          ICEL=0
      ENDIF

      ! HF and all bands simultaneous update 
      ! use more bands at a time (requires more storage but is faster)
      NSIM=WDES%NSIM*2
#ifdef MPI
      ! MPI dividable by WDES%COMM_INTER%NCPU
      NSIM=((NSIM+WDES%COMM_INTER%NCPU-1)/WDES%COMM_INTER%NCPU)*WDES%COMM_INTER%NCPU
#endif
      IF ( N <= ABS(INFO%NELMDL)) THEN
         IF (ABS(INFO%NELMDL) >=3 .AND. N<ABS(INFO%NELMDL) .AND. LHFCALC) THEN
         ! gK: if delay is larger 3, use DFT Hamiltonian in the delay phase
            CALL EDDAV(HAMILTONIAN,P, GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WDES, NSIM, &
                 LMDIM,CDIJ,CQIJ, RMS,DESUM1,ICOUEV, SV, E%EXHF, IO%IU6,IO%IU0, &
                 LDELAY=.FALSE.,  LSUBROTI=.TRUE. ,  LEMPTY=.FALSE. )
            IF (N+1==ABS(INFO%NELMDL) .AND. .NOT.INFO%LCHCON) INFO%LPOTOK=.FALSE.
         ELSE
#ifdef fock_dblbuf
            IF (LFOCK_ACE()) THEN
               CALL FOCK_ACE_ALLOCW(WDES)
               CALL FOCK_ALL_DBLBUF(WDES,W,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,CQIJ, &
              &   EX=E%EXHF,EX_ACFDT=E%EXHF_ACFDT,XI=WACE)
               CALL FOCK_ACE_CONSTRUCT(WDES,W)
            ENDIF
#endif
            CALL EDDAV(HAMILTONIAN,P, GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WDES, NSIM, &
                 LMDIM,CDIJ,CQIJ, RMS,DESUM1,ICOUEV, SV, E%EXHF, IO%IU6,IO%IU0, &
                 LDELAY=.FALSE.,  LSUBROTI=.TRUE. ,  LEMPTY=.FALSE. , LHF= .TRUE.)
#ifdef fock_dblbuf
            CALL FOCK_ACE_DEALLOCW
#endif
         ENDIF

         DWRITE0 'eddav is ok'
         ! set the W_F%CELTOT from which partial occupation-numbers are calculated
         W_F%CELTOT(:,:,:)=W%CELTOT(:,:,:)
!$ACC UPDATE DEVICE(W_F%CELTOT) __IF_ASYNC__
         CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
              INFO%NUP_DOWN, E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
              NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)

         IFLAG=0  ! signal that no corrector step is performed
      ELSE
         ! do not conjugate during delay
         IF (N> ABS(INFO%NELMDL)+1 .AND. N/=1) ICONJU=1
         ! calculate (possible) energy contribution due to finite electric field
         CALL PEAD_EDOTP(W,P,CQIJ,LATT_CUR,T_INFO,E)

         BTRIAL=INFO%TIME

!$ACC UPDATE DEVICE(W%FERTOT) IF(ACC_EXEC_ON)
         CALL EDWAV(HAMILTONIAN,KINEDEN, &
             INFO,ICONJU,IROT,ICEL,IRET,IO,BTRIAL,EFERMI, &
             ORT,RMS, E, TOTEN,TOTEN2,DESUM1, &
             GRID,KPOINTS,LATT_CUR,NONLR_S,NONL_S, T_INFO, P, W,WDES,W_F,W_G, &
             LMDIM, NSIM, CQIJ,CDIJ,CHAM,CHF,SV, &
             SYMM,GRID_SOFT,GRIDC,GRIDB,GRIDUS, &
             C_TO_US,B_TO_C,SOFT_TO_C,DENCOR,CSTRF, &
             MIX,N_MIX_PAW,IRDMAX,CHTOTL,RHOLM_LAST)

         DWRITE0 'edwav is ok'
         IFLAG=ICONJU
      ENDIF
!---- calculate old band structure energy
      E%EBANDSTR=BANDSTRUCTURE_ENERGY(WDES, W)
!---- old total energy
      TOTEN=E%EBANDSTR+E%DENC+E%XCENC+E%TEWEN+E%PSCENC+E%EENTROPY+E%PAWPS+E%PAWAE+INFO%EALLAT+E%EXHF+E%EDOTP+Ediel_sol
!---- for constrained moment calculations
      IF (M_CONSTRAINED()) TOTEN=TOTEN+E_CONSTRAINT()
      io_begin
      CALL WRITE_CONSTRAINED_M(17,.FALSE.)
      io_end

      DESUM(N)=TOTEN-TOTENL
      ECONV=DESUM(N)

 305  FORMAT(A3,': ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5, &
     &       I6,'  ',2E10.3)
 306  FORMAT('DAV: ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5, &
     &       I6,'  ',2E10.3)

      io_begin
      IF (N <= ABS(INFO%NELMDL)) THEN
         WRITE(17,306)   N,TOTEN,DESUM(N),DESUM1,WDES%NB_TOT*WDES%NKPTS,RMS
         IF (IO%IU0>=0) &
         WRITE(IO%IU0,306)  N,TOTEN,DESUM(N),DESUM1,WDES%NB_TOT*WDES%NKPTS,RMS
      ELSE
         CRET='SDA'
         SELECT CASE (IRET)
         CASE (1) 
            CRET='CGA'
         CASE (3)
            CRET='DMP'
         CASE (4)
            CRET='QMN'
         END SELECT
         WRITE(17,305)   CRET,N,TOTEN,DESUM(N),DESUM1,WDES%NB_TOT*WDES%NKPTS,RMS,ORT
         IF (IO%IU0>=0) &
         WRITE(IO%IU0,305) CRET,N,TOTEN,DESUM(N),DESUM1,WDES%NB_TOT*WDES%NKPTS,RMS,ORT
      ENDIF
      io_end
      CALL STOP_TIMING("G",IO%IU6,'TRIAL ')

!-----------------------------------------------------------------------
!   CORRECTOR STEPS if necessary
!-----------------------------------------------------------------------
      cor: DO
      IF (IFLAG==0) EXIT cor

      E2=E
      IF (ICEL==0) THEN
!---- new partial occupancies from CELTOT
         CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
              INFO%NUP_DOWN, E2%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
              NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
      ELSE IF (ICEL==1.OR. ICEL==2) THEN
!---- new partial occupations from W_F%CELTOT (requires some fiddling)
         ! use W_G temporarily
         W_G%CELTOT=W%CELTOT;  W%CELTOT=W_F%CELTOT
         CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
              INFO%NUP_DOWN, E2%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
              NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
         W%CELTOT=W_G%CELTOT ! restore W%CELTOT
      ENDIF

!---- update of charge if necessary
      IF (.NOT.INFO%LCHCON .AND. N > ABS(INFO%NELMDL)) THEN

         CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
              GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
              LATT_CUR, P, SYMM, T_INFO, &
              CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

         CALL SET_KINEDEN(GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM, &
              T_INFO%NIONS,W,WDES,KINEDEN)      

         CALL POTLOK(GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES, &
                  INFO,P,T_INFO,E2,LATT_CUR, &
                  CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF)

         CALL POTLOK_METAGGA(KINEDEN, &
                  GRID,GRIDC,GRID_SOFT,WDES%COMM_INTER,WDES,INFO,P,T_INFO,E2,LATT_CUR, &
                  CHDEN,CHTOT,DENCOR,CVTOT,SV,HAMILTONIAN%MUTOT,HAMILTONIAN%MU,SOFT_TO_C,XCSIF)

         CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
                  LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

         CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
            WDES%NCDIJ, LMDIM, CDIJ(1,1,1,1),  RHOLM, CRHODE(1,1,1,1), &
            E2,  LMETA=.FALSE., LASPH=INFO%LASPH, LCOREL= .FALSE.  )

         CALL UPDATE_CMBJ(GRIDC,T_INFO,LATT_CUR,IO%IU6)

      ENDIF

!---- recalculate eigenvalues, if not already calculated
      IFLAG=0
!$ACC KERNELS PRESENT(W) __IF_ASYNC__
      W%CELTOT=0
!$ACC END KERNELS

      CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
           LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IO%IU0,E2%EXHF)

!$ACC UPDATE SELF(W%CW,W%CPROJ,W%CELTOT) WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
      CALL PEAD_EDOTP(W,P,CQIJ,LATT_CUR,T_INFO,E2)

      E2%EBANDSTR= BANDSTRUCTURE_ENERGY(WDES, W)
      TOTEN2=E2%EBANDSTR+E2%DENC+E2%XCENC+E2%TEWEN+E2%PSCENC+E2%EENTROPY+E2%PAWPS+E2%PAWAE+INFO%EALLAT+E2%EXHF+E2%EDOTP+Ediel_sol
!---- in case of constrained moment calculations
      IF (M_CONSTRAINED()) TOTEN2=TOTEN2+E_CONSTRAINT()
      io_begin
      CALL WRITE_CONSTRAINED_M(17,.FALSE.)
      io_end

!---- correct the trial step
!     for testing set IFLAG to 10 than several trial steps are done
      IFLAG=2

      DWRITE0 'edwav call '
      CALL EDWAV(HAMILTONIAN,KINEDEN, &
             INFO,IFLAG,IROT,ICEL,IRET,IO,BTRIAL,EFERMI, &
             ORT,RMS, E, TOTEN,TOTEN2,DESUM1, &
             GRID,KPOINTS,LATT_CUR,NONLR_S,NONL_S, T_INFO, P, W,WDES,W_F,W_G, &
             LMDIM, NSIM, CQIJ,CDIJ,CHAM,CHF,SV, &
             SYMM,GRID_SOFT,GRIDC,GRIDB,GRIDUS, &
             C_TO_US,B_TO_C,SOFT_TO_C,DENCOR,CSTRF, &
             MIX,N_MIX_PAW,IRDMAX,CHTOTL,RHOLM_LAST)
      DWRITE0 'edwav is ok'

      ENDDO cor

!-----------------------------------------------------------------------
!   end of loop update fermi-weights
!   and check for abort condition
!-----------------------------------------------------------------------
      IF (ICEL==0) THEN
!---- new partial occupancies from CELTOT
         CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
              INFO%NUP_DOWN, E2%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
              NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
      ELSE IF (ICEL==1 .OR. ICEL==2) THEN
!---- new partial occupancies from CELTOT_F (see above)
        W_G%CELTOT=W%CELTOT;  W%CELTOT=W_F%CELTOT
        CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
                INFO%NUP_DOWN, E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
                NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
        W%CELTOT=W_G%CELTOT ! restore W%CELTOT
      ENDIF


      CALL STOP_TIMING("G",IO%IU6,'CORREC')
!
!---- test for break condition
!
      INFO%LABORT=.FALSE.
!---- energy converged
!     DESUM1 is the expected energy difference between the current wavefunctions
!     and the previous one, it is not yet known
!     DESUM(N) contains the previous change of energy
!     both must be smaller then the threshhold (requiring that the energy
!     does not change for two consecutive steps)
      IF(ABS(DESUM1)<INFO%EDIFF .AND. (ABS(DESUM(N))<INFO%EDIFF.OR. IS_MODEL_GW())) INFO%LABORT=.TRUE.
      IF(LPEAD_ABORT(INFO,RMS)) INFO%LABORT=.TRUE.
      CALL testNumberOfStep(N, NELM, INFO, LABORT_WITHOUT_CONV)

!      IF (INFO%LABORT.AND. &
!         ((NINT(INFO%NELECT)==WDES%NB_TOT .AND. WDES%LNONCOLLINEAR) .OR. &
!          (NINT(INFO%NELECT)==WDES%NB_TOT*2 .AND. .NOT. WDES%LNONCOLLINEAR))) THEN
! test
!     IF (INFO%LABORT) THEN
!        TOTEN=TOTEN+PEAD_EDOTP(W,P,CQIJ,LATT_CUR,T_INFO,E,LFORCE=.TRUE.)
!     ENDIF
! test

      IF (INFO%LABORT .AND. (INFO%LDIAG.OR.LUSEPEAD())) THEN
!$ACC UPDATE SELF(W%CW,W%CPROJ,W%CELTOT) WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
        CALL PEAD_ACC_CALC_ALL(W,KPOINTS,P,CQIJ,LATT_CUR,T_INFO)

        IF (IO%IU0>=0) WRITE(IO%IU0,*)'final diagonalization'
! the final diagonalization is desirable even for IROT=3 (used for LDIAG = .FALSE.)
! however, presently for IROT=3  and LDIAG = .FALSE. the code does not  get here
        IF (IROT==3 .OR. IROT==4 .OR. IROT==5 .OR. IROT==6 .OR. IROT==7 ) THEN
           IFLAG=23  ! rotate only in fully occupied many-fold
           CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
                LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IO%IU0,E%EXHF)
        ELSE
           IFLAG=3   ! rotate in entire sub-space of orbitals
           CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
                LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IO%IU0,E%EXHF)

           CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
                INFO%NUP_DOWN, E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
                NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
        ENDIF

        CALL STOP_TIMING("G",IO%IU6,'EDDIAG')
      ENDIF
! ======================================================================
! If the end of the electronic loop is reached
! calculate accurate initial state core level shifts
! if required
! ======================================================================
      IF (INFO%LABORT .AND. ACCURATE_CORE_LEVEL_SHIFTS()) THEN

         ALLOCATE(CDIJ_TMP(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ))
         CDIJ_TMP=CDIJ

         CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
              WDES%NCDIJ, LMDIM, CDIJ_TMP(1,1,1,1), RHOLM, CRHODE, &
              E, LMETA= .FALSE. , LASPH=INFO%LASPH , LCOREL=.TRUE. )
         DEALLOCATE(CDIJ_TMP)
      ENDIF
!========================= subroutine CHSP  ============================
! if charge density is updated
!  ) first copy current charge to CHTOTL
!  ) set  INFO%LPOTOK to .F. this requires a recalculation of the local pot.
!  ) set INFO%LMIX to .T.
!  ) call subroutine CHSP+ DEPLE to generate the new charge density
!  ) then performe mixing
! MIND:
! ) if delay is selected  do not update
! ) if convergence corrections to forces are calculated do not update charge
!   in last iteration
!=======================================================================

      INFO%LMIX=.FALSE.
      MIX%NEIG=0

      IF (.NOT.INFO%LCHCON.AND..NOT.(INFO%LABORT.AND.INFO%LCORR).AND.N>=ABS(INFO%NELMDL)) THEN
         CALL START_TIMING("G")

         IF (N<=ABS(INFO%NELMDL)) THEN
            PUSH_ACC_EXEC_ON(.FALSE.)
            DO ISP=1,WDES%NCDIJ
               CALL RC_ADD(CHTOT(1,ISP),1.0_q,CHTOT(1,ISP),0.0_q,CHTOTL(1,ISP),GRIDC)
            ENDDO
            RHOLM_LAST=RHOLM
            POP_ACC_EXEC_ON
         ENDIF

         CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
              GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
              LATT_CUR, P, SYMM, T_INFO, &
              CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

         CALL SET_KINEDEN(GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM, &
              T_INFO%NIONS,W,WDES,KINEDEN)

         INFO%LPOTOK=.FALSE.

         CALL STOP_TIMING("G",IO%IU6,'CHARGE')

!=======================================================================
! possibly we will mix here let's see (does not seem to work)
!=======================================================================
         IF (MIX%IMIX==2) THEN
            CALL START_TIMING("G")
            INFO%LMIX=.TRUE.
!---- simple mixing ... :
            CALL MIX_SIMPLE(GRIDC,MIX,WDES%NCDIJ, CHTOT,CHTOTL, &
            N_MIX_PAW, RHOLM, RHOLM_LAST, LATT_CUR%B, LATT_CUR%OMEGA, RMST)

            CALL STOP_TIMING("G",IO%IU6,'MIXING')
            DWRITE0 "mixing is ok"
         ENDIF

!---- end of charge update
      ENDIF

      CALL START_TIMING("G")
      IF (W%OVER_BAND) THEN
         CALL REDIS_PW_OVER_BANDS(WDES, W)
         CALL STOP_TIMING("G",IO%IU6,'REDIS')
      ENDIF
!=======================================================================
! total time used for this step
!=======================================================================
      CALL SEPERATOR_TIMING(IO%IU6)
      CALL STOP_TIMING("LOOP",IO%IU6,XMLTAG='total')
!=======================================================================
!  important write statements
!=======================================================================

 2440 FORMAT(/' eigenvalue-minimisations  :',I6,/ &
     &       ' total energy-change (2. order) :',E14.7,'  (',E14.7,')')

 200  FORMAT(' number of electron ',F15.7,' magnetization ',3F15.7)
 201  FORMAT(' augmentation part  ',F15.7,' magnetization ',3F15.7)

      NORDER=0
      IF (KPOINTS%ISMEAR>=0) NORDER=KPOINTS%ISMEAR

      DO I=1,WDES%NCDIJ
         RHOTOT(I)=RHO0(GRIDC, CHTOT(1,I))
         RHOAUG(I)=RHOTOT(I)-RHO0(GRID_SOFT, CHDEN(1,I))
      END DO

      io_begin

    ! iteration counts
      WRITE(IO%IU6,2440) ICOUEV,DESUM(N),DESUM1

    ! charge density
      WRITE(IO%IU6,200) RHOTOT
      IF (INFO%LOVERL) THEN
        WRITE(IO%IU6,201) RHOAUG
      ENDIF
    ! dipol moment
      IF (DIP%LCOR_DIP) CALL WRITE_DIP(IO%IU6)

      io1: IF (IO%NWRITE>=2 .OR. (NSTEP==1)) THEN
       ! energy
         WRITE(IO%IU6,7240) E%PSCENC,E%TEWEN,E%DENC,E%EXHF,E%XCENC,E%PAWPS,E%PAWAE, &
                            E%EENTROPY,E%EBANDSTR,INFO%EALLAT,TOTEN, &
                            TOTEN-E%EENTROPY,TOTEN-E%EENTROPY/(2+NORDER)

         IF (LHFCALC) WRITE(IO%IU6,'( "  exchange ACFDT corr.  = ",F18.8,"  see jH, gK, PRB 81, 115126")') E%EXHF_ACFDT

 7240    FORMAT(/ &
        &        ' Free energy of the ion-electron system (eV)'/ &
        &        '  ---------------------------------------------------'/ &
        &        '  alpha Z        PSCENC = ',F18.8/ &
        &        '  Ewald energy   TEWEN  = ',F18.8/ &
        &        '  -Hartree energ DENC   = ',F18.8/ &
        &        '  -exchange      EXHF   = ',F18.8/ &
        &        '  -V(xc)+E(xc)   XCENC  = ',F18.8/ &
        &        '  PAW double counting   = ',2F18.8/ &
        &        '  entropy T*S    EENTRO = ',F18.8/ &
        &        '  eigenvalues    EBANDS = ',F18.8/ &
        &        '  atomic energy  EATOM  = ',F18.8/ &
        &        '  ---------------------------------------------------'/ &
        &        '  free energy    TOTEN  = ',F18.8,' eV'// &
        &        '  energy without entropy =',F18.8, &
        &        '  energy(sigma->0) =',F18.8)
 72612   FORMAT(//&
        &        '  METAGGA EXCHANGE AND CORRELATION (eV)'/ &
        &        '  ---------------------------------------------------'/ &
        &        '  LDA+GGA E(xc)  EXCG   = ',F18.6/ &
        &        '  LDA+GGA PAW    PS : AE= ',2F18.6/ &
        &        '  core xc             AE= ',1F18.6/ &
        &        '  metaGGA E(xc)  EXCM   = ',F18.6/ &
        &        '  metaGGA PAW    PS : AE= ',2F18.6/ &
        &        '  metaGGA core xc     AE= ',1F18.6/ &
        &        '  ---------------------------------------------------'/ &
        &        '  METAGGA result:'/ &
        &        '  free  energy   TOTEN  = ',F18.6,' eV'// &
        &        '  energy  without entropy=',F18.6, &
        &        '  energy(sigma->0) =',F16.6)
      ELSE io1
         WRITE(IO%IU6,7242) TOTEN,TOTEN-E%EENTROPY
 7242    FORMAT(/'  free energy = ',E20.12, &
        &        '  energy without entropy= ',E20.12)

      ENDIF io1

!     too slow on many servers nowadays
!     IF (IO%LOPEN) CALL WFORCE(IO%IU6)
!     IF (IO%LOPEN) CALL WFORCE(17)
      WRITE(IO%IU6,130)
      io_end
!=======================================================================
!  perform some additional write statments if required
!=======================================================================
!-----Eigenvalues and weights
      IF (((NSTEP==1.OR.NSTEP==DYN%NSW).AND.INFO%LABORT).OR.(IO%NWRITE>=1.AND.INFO%LABORT).OR.IO%NWRITE>=3) THEN
         CALL KPAR_SYNC_CELTOT(WDES,W)

         CALL MRG_AUX(WDES,W)

         IF (IO%IU6>=0) &
            WRITE(IO%IU6,'("  average scaling for gradient ",F8.4)') SUM(W%AUXTOT)/SIZE(W%AUXTOT)

         ! calculate the core level shifts
         IF (INFO%LOVERL) THEN
            CALL CL_SHIFT_PW( GRIDC, LATT_CUR, IRDMAX,  &
            T_INFO, P, WDES%NCDIJ, CVTOT, MAX(INFO%ENAUG,INFO%ENMAX), IO%IU6)
         ELSE
            do_io WRITE(*,*) " **** core level shifts not calculated ****"
         ENDIF

         io_begin
         CALL RHOAT0(P,T_INFO, BETATO,LATT_CUR%OMEGA)

         WRITE(IO%IU6,2202) EFERMI,REAL( E%CVZERO ,KIND=q) ,E%PSCENC/INFO%NELECT+BETATO
 2202    FORMAT(' E-fermi : ', F8.4,'     XC(G=0): ',F8.4, &
        &       '     alpha+bet :',F8.4/)

         CALL WRITE_EIGENVAL( WDES, W, IO%IU6)
! at this point the CELTOT contains the eigenvalues
! <phi_i | T+ V_H + V_COH + V_SEX | phi_i>
! you need to call a subroutine that calculates the dynamical correction
! prints them out and writes out the corrected eigenvalues *including dynamical corrections
         CALL GWDYNSM1_IJ(W, WDES, KPOINTS, LATT_CUR, T_INFO, INFO, IO, E_DFT)

!-----Charge-density along one line
         WRITE(IO%IU6,130)
         DO I=1,WDES%NCDIJ
            WRITE(IO%IU6,*)'soft charge-density along one line, spin component',I
            WRITE(IO%IU6,'(10(6X,I4))') (II,II=0,9)
            CALL WRT_RC_LINE(IO%IU6,GRID_SOFT, CHDEN(1,I))
            IF (INFO%LOVERL) THEN
               WRITE(IO%IU6,*)'total charge-density along one line'
               CALL WRT_RC_LINE(IO%IU6,GRIDC, CHTOT(1,I))
            ENDIF
            WRITE(IO%IU6,*)
         ENDDO
!-----pseudopotential strength and augmentation charge
         DO I=1,WDES%NCDIJ
            WRITE(IO%IU6,*) 'pseudopotential strength for first ion, spin component:',I
            DO LP=1,P(1)%LMMAX
               WRITE(IO%IU6,'(16(F7.3,1X))') (CDIJ(L,LP,1,I),L=1,MIN(8,P(1)%LMMAX))
!                    (REAL(CDIJ(L,LP,1,I),q),L=1,MIN(16,P(1)%LMMAX))
            ENDDO
         ENDDO

         IF (INFO%LOVERL) THEN
            DO I=1,WDES%NCDIJ
               WRITE(IO%IU6,*) 'total augmentation occupancy for first ion, spin component:',I
               DO LP=1,P(1)%LMMAX
                  WRITE(IO%IU6,'(16(F7.3,1X))') (REAL(CRHODE(L,LP,1,I),q),L=1,MIN(16,P(1)%LMMAX))
               ENDDO
            ENDDO
         ENDIF
         io_end

      ENDIF
!=======================================================================
!  xml related output
!=======================================================================
      CALL XML_TAG("energy")
      IF (INFO%LABORT .OR. N==1) THEN
         CALL XML_TAG_REAL("alphaZ",E%PSCENC)
         CALL XML_TAG_REAL("ewald", E%TEWEN)
         CALL XML_TAG_REAL("hartreedc",E%DENC)
         CALL XML_TAG_REAL("XCdc",E%XCENC)
         CALL XML_TAG_REAL("pawpsdc",E%PAWPS)
         CALL XML_TAG_REAL("pawaedc",E%PAWAE)
         CALL XML_TAG_REAL("eentropy",E%EENTROPY)
         CALL XML_TAG_REAL("bandstr",E%EBANDSTR)
         CALL XML_TAG_REAL("atom",INFO%EALLAT)
         CALL XML_ENERGY(TOTEN, TOTEN-E%EENTROPY, TOTEN-E%EENTROPY/(2+NORDER))
      ELSE
         CALL XML_ENERGY(TOTEN, TOTEN-E%EENTROPY, TOTEN-E%EENTROPY/(2+NORDER))
      ENDIF
      CALL XML_CLOSE_TAG

!      IF (INFO%LABORT .AND. INFO%LMETAGGA) THEN
!         CALL XML_TAG("metagga")
!         CALL XML_TAG_REAL("e_fr_energy",E%TOTENMGGA)
!         CALL XML_TAG_REAL("e_wo_entrp", E%TOTENMGGA-E%EENTROPY)
!         CALL XML_TAG_REAL("e_0_energy", E%TOTENMGGA-E%EENTROPY/(2+NORDER))
!         CALL XML_CLOSE_TAG
!      ENDIF

      CALL XML_CLOSE_TAG("scstep")
!======================== end of loop ENDLSC ===========================
! This is the end of the selfconsistent calculation loop
!=======================================================================
      CALL checkAbort(LABORT_WITHOUT_CONV, IO, W%WDES%COMM, INFO)
      IF (INFO%LABORT.OR.INFO%LSOFT) EXIT electron
      TOTENL=TOTEN

      ENDDO electron

#ifdef _OPENACC
!$ACC WAIT IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(W%CW,W%CPROJ,W%CELTOT) IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(SV) IF(ACC_EXEC_ON)

!$ACC EXIT DATA DELETE(CHF,CHAM) IF(ACC_EXEC_ON)
!$ACC EXIT DATA DELETE(SV,CDIJ,CQIJ) IF(ACC_EXEC_ON)
!$ACC EXIT DATA DELETE(DENCOR,CSTRF,LATT_CUR) IF(ACC_EXEC_ON)
!$ACC EXIT DATA COPYOUT(KINEDEN%TAU,KINEDEN) IF(ASSOCIATED(KINEDEN%TAU).AND.ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA COPYOUT(KINEDEN%TAUC,KINEDEN) IF(ASSOCIATED(KINEDEN%TAUC).AND.ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA DELETE(HAMILTONIAN%AVEC) IF(ASSOCIATED(HAMILTONIAN%AVEC).AND.ACC_EXEC_ON)
!$ACC EXIT DATA DELETE(HAMILTONIAN%MU) IF(ASSOCIATED(HAMILTONIAN%MU).AND.ACC_EXEC_ON)
!$ACC EXIT DATA DELETE(HAMILTONIAN) IF(ACC_EXEC_ON)

      CALL ACC_DELETE_TYPED_VAR(P)
      CALL ACC_DELETE_TYPED_VAR(C_TO_US)
      CALL ACC_DELETE_TYPED_VAR(SOFT_TO_C)
      CALL ACC_DELETE_TYPED_VAR(GRID_SOFT)
      CALL ACC_DELETE_TYPED_VAR(GRIDC)
      CALL ACC_DELETE_TYPED_VAR(GRID)

      CALL ACC_DELETE_TYPED_VAR(NONLR_S)
      CALL ACC_DELETE_TYPED_VAR(NONL_S)
      CALL ACC_DELETE_TYPED_VAR(W_G)
      CALL ACC_DELETE_TYPED_VAR(W_F)
      CALL ACC_DELETE_TYPED_VAR(W)
      CALL ACC_DELETE_TYPED_VAR(WDES)

      POP_ACC_EXEC_ON
#endif
!
! calculate dipol corrections now
!
      IF ( DIP%IDIPCO >0 ) THEN
         IF (.NOT. DIP%LCOR_DIP) THEN
            CALL CDIPOL_CHTOT_REC(GRIDC, LATT_CUR,P,T_INFO, &
            CHTOT,CSTRF,CVTOT, WDES%NCDIJ, INFO%NELECT )

            CALL WRITE_DIP(IO%IU6)
            IF (IO%IU6>0) THEN
               WRITE(IO%IU6,*)
               WRITE(IO%IU6,*) &
               " *************** adding dipol energy to TOTEN NOW **************** "
            ENDIF
         TOTEN=TOTEN+ DIP%ECORR
         ENDIF
      ENDIF

      ! notify calling routine whether convergence has been reached
      INFO%LABORT=LABORT_WITHOUT_CONV

      DWRITE0 'electronall left'

      PROFILING_STOP('elmin_all')

      RETURN
      END SUBROUTINE
