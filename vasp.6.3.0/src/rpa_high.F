#include "symbol.inc"

MODULE rpa_high
  USE rpa_force
  CONTAINS
!*************** SUBROUTINE EDDIAG_EXACT_UPDATE_NBANDS *****************
!
! this subroutine performs a full diagonalization of the Hamiltonian
! and keeps the band number to the large number so determined
!
!***********************************************************************

  SUBROUTINE EDDIAG_EXACT_UPDATE_NBANDS(HAMILTONIAN,KINEDEN, &
       GRID,GRID_SOFT,GRIDC,GRIDB,GRIDUS,C_TO_US,SOFT_TO_C,B_TO_C,E, &
       CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
       IRDMAX,CRHODE,MIX,N_MIX_PAW,RHOLM,RHOLM_LAST,CHDEN, &
       LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
       LMDIM,CDIJ,CQIJ,SV,VDW_SET,T_INFO, DYN, P, IO, INFO, &
       XCSIF, EWSIF, TSIF, EWIFOR, TIFOR, PRESS, TOTEN, KPOINTS, EFERMI, NEDOS, NBANDSEXACT )

    USE prec
    USE wave_high
    USE lattice
    USE mpimy
    USE mgrid
    USE nonl_high
    USE hamil_struct_def
    USE main_mpi
    USE pseudo
    USE poscar
    USE ini
    USE choleski
    USE fock
    USE scala
    USE setexm
    USE meta
    USE us
    USE pawm
    USE vdwd4, ONLY: vdw_settings
    USE mkpoints
    USE mlr_optic
    USE subrot
    USE pead
    IMPLICIT NONE
    TYPE (ham_handle)  HAMILTONIAN
    TYPE (tau_handle)  KINEDEN
    TYPE (latt)        LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W
    TYPE (wavedes)     WDES
    TYPE (symmetry)    SYMM
    INTEGER LMDIM
    OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ),CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    OVERLAP CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    TYPE (type_info)   T_INFO
    TYPE (kpoints_struct) :: KPOINTS
    TYPE (dynamics)    DYN
    TYPE (grid_3d)     GRID
    TYPE (grid_3d)     GRID_SOFT  ! grid for soft chargedensity
    TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
    TYPE (grid_3d)     GRIDB      ! grid for Broyden mixer
    TYPE (grid_3d)     GRIDUS     ! temporary grid in us.F
    TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
    TYPE (transit)     B_TO_C     ! index table between GRIDB and GRIDC
    TYPE (transit)     SOFT_TO_C  ! index table between GRID_SOFT and GRIDC
    TYPE (vdw_settings), INTENT(IN) :: VDW_SET
    COMPLEX(q)  CHTOT(GRIDC%MPLWV,WDES%NCDIJ) ! charge-density in real / reciprocal space
    COMPLEX(q)  CHTOTL(GRIDC%MPLWV,WDES%NCDIJ)! old charge-density
    RGRID       DENCOR(GRIDC%RL%NP)           ! partial core
    COMPLEX(q)  CVTOT(GRIDC%MPLWV,WDES%NCDIJ) ! local potential
    COMPLEX(q)  CSTRF(GRIDC%MPLWV,T_INFO%NTYP)! structure factor
    RGRID   SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ) ! local potential
    COMPLEX(q)  CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
    INTEGER            IRDMAX, N_MIX_PAW
    REAL(q)  RHOLM(N_MIX_PAW,WDES%NCDIJ)
    REAL(q)  RHOLM_LAST(N_MIX_PAW,WDES%NCDIJ)
    TYPE (potcar)      P(T_INFO%NTYP)
    TYPE (energy)      E
    TYPE (in_struct)   IO
    TYPE (mixing)      MIX
    TYPE (info_struct) INFO
    REAL(q)   XCSIF(3,3)                      ! stress stemming from XC
    REAL(q)   EWSIF(3,3)                      ! stress from Ewald contribution
    REAL(q)   TSIF(3,3)                       ! total stress (set by routine)
    REAL(q)   EWIFOR(3,T_INFO%NIOND)          ! ewald force
    REAL(q)   TIFOR(3,T_INFO%NIOND)           ! total force (set by routine)
    REAL(q)   PRESS                           ! external pressure
    REAL(q)   TOTEN
    INTEGER   NBANDSEXACT                     ! bands to be used in exact diagonalization
    REAL(q)   EFERMI                          ! Fermi-energy
    INTEGER   NEDOS
  ! local
    INTEGER NB_TOT            ! maximum number of plane wave coefficients = number of bands
    TYPE (wavespin)    W_TMP
    INTEGER ISP, NK, NB, DEGREES_OF_FREEDOM
    INTEGER            IFLAG  ! determines mode of diagonalisation
    INTEGER NB_TOT_OLD, NBANDS_OLD
    LOGICAL :: LCORR_TMP
    INTEGER,ALLOCATABLE :: NB_TOTK(:,:)

  ! just make sure that data distribution is over bands
    CALL REDIS_PW_OVER_BANDS(WDES, W)

    ! recalculate charge density and  kinetic energy density
    CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
         GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
         LATT_CUR, P, SYMM, T_INFO, &
         CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

    CALL SET_KINEDEN(GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM, &
         T_INFO%NIONS,W,WDES,KINEDEN)      

  ! calculate local potential
    CALL UPDATE_POT

    LCORR_TMP=INFO%LCORR
    INFO%LCORR=.FALSE.

  ! we include the contributions \sum_ij H^HF_ij <phi_i| dp/dR Q |phi_j> + c.c,
  ! where i and j are both occupied,
  ! in the routine FORNL2 and FORNLR2 in rpa_force
  ! by setting CELTOT to 0, these contributions are bypassed here
    W%CELTOT=0
  ! Actually that contributions would be correct only if the
  ! HF Hamiltonian were diagonalized beforehand; this is obviously not the case yet
  ! furthermore LREMOVE_DRIFT, is .false. i.e. no poking with forces
    IF (IO%IU6>=0) WRITE(IO%IU6,*) "HF-forces evaluated using KS orbitals"
    IF (IO%IU6>=0) WRITE(IO%IU6,*) "-------------------------------------"
    CALL FORCE_AND_STRESS( &
         KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR, &
         T_INFO,T_INFO,DYN,INFO,IO,MIX,SYMM,GRID,GRID_SOFT, &
         GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C, &
         CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
         CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
         CHDEN,SV,VDW_SET, &
         LMDIM, IRDMAX, .TRUE., &
         DYN%ISIF/=0, DYN%ISIF/=0,.FALSE.,  &
         XCSIF, EWSIF, TSIF, EWIFOR, TIFOR, PRESS, TOTEN, KPOINTS )
    ! force routine updated NONLR_S, so reset it
    IF (INFO%LREAL) THEN
       CALL RSPHER(GRID,NONLR_S,LATT_CUR)
    ENDIF

    INFO%LCORR=LCORR_TMP

    PROFILING_START('eddiag_exact_update_nbands')

  ! the first step is to restore the original XC-functional (as read from POTCAR/INCAR)
  ! in order to use the appropriate Hamiltonian
    CALL POP_XC_TYPE
    IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
       CALL SETUP_LDA_XC(2,-1,-1,IO%IDIOT)
    ELSE
       CALL SETUP_LDA_XC(1,-1,-1,IO%IDIOT)
    ENDIF
  ! now update the PAW one-center terms to current functional
    CALL SET_PAW_ATOM_POT( P , T_INFO, WDES%LOVERL, LMDIM, INFO%EALLAT, INFO%LMETAGGA, IO%IU6 )
    ! and restore the convergence corrections
    CALL SET_FSG_STORE(GRIDHF, LATT_CUR, WDES)
    CALL UPDATE_POT

    DEGREES_OF_FREEDOM=MAXVAL(WDES%NPLWKP_TOT)
    IF (WDES%LGAMMA) THEN
       DEGREES_OF_FREEDOM=DEGREES_OF_FREEDOM*2-1
    ENDIF

    IF (DEGREES_OF_FREEDOM<=WDES%NB_TOT) THEN
  ! smaller than already included bands
  ! restore DFT eigenvalue
       IFLAG=3
       CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
            LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IO%IU0,E%EXHF,EXHF_ACFDT=E%EXHF_ACFDT)
    ELSE
       NB_TOT_OLD=WDES%NB_TOT
       NBANDS_OLD=WDES%NBANDS
       NB_TOT=((DEGREES_OF_FREEDOM+WDES%NB_PAR-1)/WDES%NB_PAR)*WDES%NB_PAR
    
       WDES%NB_TOT=NB_TOT
       WDES%NBANDS=NB_TOT/WDES%NB_PAR
       CALL INIT_SCALAAWARE( WDES%NB_TOT, WDES%NRPLWV, WDES%COMM_KIN )
       
! set the maximum number of bands k-point dependent
       DO NK=1,WDES%NKPTS
          IF (WDES%LGAMMA) THEN
             WDES%NB_TOTK(NK,:)=MIN(WDES%NB_TOT,WDES%NPLWKP_TOT(NK)*2-1)
          ELSE
             WDES%NB_TOTK(NK,:)=MIN(WDES%NB_TOT,WDES%NPLWKP_TOT(NK))
          ENDIF
       ENDDO
       CALL RESETUP_FOCK_WDES(WDES, LATT_CUR, LATT_CUR, -1)

    ! copy W to new W_TMP    
       CALL ALLOCW(WDES, W_TMP)
       W_TMP%FERTOT(:,:,:)=0
       W_TMP%CELTOT(:,:,:)=0
    ! copy data back to work array
       W_TMP%CW(:,1:NBANDS_OLD,:,:)    =W%CW(:,1:NBANDS_OLD,:,:)
       W_TMP%CPROJ(:,1:NBANDS_OLD,:,:) =W%CPROJ(:,1:NBANDS_OLD,:,:)
       W_TMP%CELTOT(1:NB_TOT_OLD,:,:)=W%CELTOT(1:NB_TOT_OLD,:,:)
       W_TMP%FERTOT(1:NB_TOT_OLD,:,:)=W%FERTOT(1:NB_TOT_OLD,:,:)
       CALL DEALLOCW(W)
    ! random initialization beyond WDES%NBANDS
       CALL WFINIT(WDES, W_TMP, 1E10_q, NB_TOT_OLD+1) ! ENINI=1E10 not cutoff restriction

    ! get characters
       CALL PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,W_TMP)
    ! orthogonalization
       CALL ORTHCH(WDES,W_TMP, WDES%LOVERL, LMDIM,CQIJ)
    ! and diagonalization
       IFLAG=3
       CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W_TMP,WDES,SYMM, &
            LMDIM,CDIJ,CQIJ, IFLAG,SV,T_INFO,P,IO%IU0,E%EXHF,EXHF_ACFDT=E%EXHF_ACFDT)
    ! set W=W_TMP
       W=W_TMP
    ! now if NBANDSEXACT is set, decrease the number of bands to NBANDSEXACT
       IF (NBANDSEXACT>=0 .AND. NBANDSEXACT< NB_TOT) THEN
          NBANDSEXACT=MAX(NB_TOT_OLD, NBANDSEXACT)  ! make sure user did not use stupid values
          NB_TOT=((NBANDSEXACT+WDES%NB_PAR-1)/WDES%NB_PAR)*WDES%NB_PAR
    
          WDES%NB_TOT=NB_TOT
          WDES%NBANDS=NB_TOT/WDES%NB_PAR
          CALL INIT_SCALAAWARE( WDES%NB_TOT, WDES%NRPLWV, WDES%COMM_KIN )

          ! set the maximum number of bands k-point dependent
          DO NK=1,WDES%NKPTS
             IF (WDES%LGAMMA) THEN
                WDES%NB_TOTK(NK,:)=MIN(WDES%NB_TOT,WDES%NPLWKP_TOT(NK)*2-1)
             ELSE
                WDES%NB_TOTK(NK,:)=MIN(WDES%NB_TOT,WDES%NPLWKP_TOT(NK))
             ENDIF
          ENDDO
          CALL RESETUP_FOCK_WDES(WDES, LATT_CUR, LATT_CUR, -1)
          
          CALL ALLOCW(WDES, W_TMP)
          ! copy data to W_TMP
          W_TMP%CW(:,1:WDES%NBANDS,:,:)    =W%CW(:,1:WDES%NBANDS,:,:)
          W_TMP%CPROJ(:,1:WDES%NBANDS,:,:) =W%CPROJ(:,1:WDES%NBANDS,:,:)
          W_TMP%CELTOT(1:WDES%NB_TOT,:,:)=W%CELTOT(1:WDES%NB_TOT,:,:)
          W_TMP%FERTOT(1:WDES%NB_TOT,:,:)=W%FERTOT(1:WDES%NB_TOT,:,:)
          CALL DEALLOCW(W)
          W=W_TMP
       ENDIF
       CALL DUMP_ALLOCATE(IO%IU6)
    ENDIF

    E%EBANDSTR=BANDSTRUCTURE_ENERGY(WDES, W)

    IF (IO%IU6>=0) THEN
       WRITE(IO%IU6,7240) " Exact diagonalization of KS Hamiltonian yields:", & 
            E%PSCENC,E%TEWEN,E%DENC,E%EXHF,E%XCENC,E%PAWPS,E%PAWAE, &
            E%EENTROPY,E%EBANDSTR,INFO%EALLAT,'  free energy    TOTEN  = ', &
            E%EBANDSTR+E%DENC+E%XCENC+E%TEWEN+E%PSCENC+E%EENTROPY+E%PAWPS+E%PAWAE+INFO%EALLAT+E%EXHF
    ENDIF

    IF (IO%LOPTICS) THEN
       CALL START_TIMING("G")
       CALL PEAD_RESETUP_WDES(WDES, GRID, KPOINTS, LATT_CUR, LATT_CUR, IO)
       CALL LR_OPTIC( &
            P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_CUR, &
            T_INFO,INFO,IO,KPOINTS,SYMM,GRID,GRID_SOFT, &
            GRIDC,GRIDUS,C_TO_US,SOFT_TO_C, &
            CHTOT,DENCOR,CVTOT,CSTRF, &
            CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM, &
            CHDEN,SV,LMDIM,IRDMAX,EFERMI,NEDOS, &
            LSTORE=.TRUE., LPOT=.FALSE.)
       CALL STOP_TIMING("G",IO%IU6,'OPTICS')
    ENDIF


    ! now we restore the HF xc-correlation functional
    CALL PUSH_XC_TYPE_FOR_GW
    IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
       CALL SETUP_LDA_XC(2,-1,-1,IO%IDIOT)
    ELSE
       CALL SETUP_LDA_XC(1,-1,-1,IO%IDIOT)
    ENDIF
    ! now update the PAW one-center terms to current functional
    CALL SET_PAW_ATOM_POT( P , T_INFO, WDES%LOVERL, LMDIM, INFO%EALLAT, INFO%LMETAGGA, IO%IU6 )

    ! and restore the convergence corrections
    CALL SET_FSG_STORE(GRIDHF, LATT_CUR, WDES)


7240 FORMAT(/ &
              A,/ &
              ' Free energy of the ion-electron system (eV)', / &
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
     &        A,F18.8,' eV' )

    CALL WRITE_EIGENVAL_NBANDS( WDES, W, IO%IU6, WDES%NB_TOT)

    PROFILING_STOP('eddiag_exact_update_nbands')

  CONTAINS 

    SUBROUTINE UPDATE_POT
      USE pot
      USE pawm
      USE morbitalmag
      USE us
      REAL(q) :: XCSIF(3,3)
      INTEGER :: IRDMAA

      CALL POTLOK(GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES, &
           INFO,P,T_INFO,E,LATT_CUR, &
           CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF)
      
      CALL POTLOK_METAGGA(KINEDEN, &
           GRID,GRIDC,GRID_SOFT,WDES%COMM_INTER,WDES,INFO,P,T_INFO,E,LATT_CUR, &
           CHDEN,CHTOT,DENCOR,CVTOT,SV,HAMILTONIAN%MUTOT,HAMILTONIAN%MU,SOFT_TO_C,XCSIF)
      
      CALL VECTORPOT(GRID, GRIDC, GRID_SOFT, SOFT_TO_C,  WDES%COMM_INTER, & 
           LATT_CUR, T_INFO%POSION, HAMILTONIAN%AVEC, HAMILTONIAN%AVTOT)
      
      
      CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
           LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)
      
      CALL SETDIJ_AVEC(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
           LMDIM,CDIJ,HAMILTONIAN%AVTOT, NONLR_S, NONL_S, IRDMAX)
      
      CALL SET_DD_MAGATOM(WDES, T_INFO, P, LMDIM, CDIJ)
      
      CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
           WDES%NCDIJ, LMDIM, CDIJ(1,1,1,1),  RHOLM, CRHODE(1,1,1,1), &
           E,  LMETA=.FALSE., LASPH=INFO%LASPH, LCOREL= .FALSE.  )
      
      CALL UPDATE_CMBJ(GRIDC,T_INFO,LATT_CUR,IO%IU6)
      
    END SUBROUTINE UPDATE_POT

  END SUBROUTINE EDDIAG_EXACT_UPDATE_NBANDS

!********************** SUBROUTINE UPDATE_CDIJ *************************
!
!  Update CDIJ to the current functional that is applied
!
!***********************************************************************

  SUBROUTINE UPDATE_CDIJ(HAMILTONIAN,KINEDEN,GRID,GRID_SOFT, &
           GRIDC,GRIDUS,C_TO_US,SOFT_TO_C,E,CHTOT,DENCOR, &
           CVTOT,CSTRF,IRDMAX,CRHODE,N_MIX_PAW,RHOLM,CHDEN, &
           LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM,LMDIM,CDIJ, &
           CQIJ,SV,T_INFO,P,IO,INFO,XCSIF)

    USE prec
    USE wave_high
    USE lattice
    USE mpimy
    USE mgrid
    USE nonl_high
    USE hamil_struct_def
    USE main_mpi
    USE pseudo
    USE poscar
    USE ini
    USE choleski
    USE fock
    USE scala
    USE setexm
    USE meta
    USE us
    USE pawm
    IMPLICIT NONE
    TYPE (ham_handle)  HAMILTONIAN
    TYPE (tau_handle)  KINEDEN
    TYPE (latt)        LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W
    TYPE (wavedes)     WDES
    TYPE (symmetry)    SYMM
    INTEGER LMDIM
    OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ),CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    OVERLAP CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    TYPE (type_info)   T_INFO
    TYPE (grid_3d)     GRID
    TYPE (grid_3d)     GRID_SOFT  ! grid for soft chargedensity
    TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
    TYPE (grid_3d)     GRIDUS     ! temporary grid in us.F
    TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
    TYPE (transit)     SOFT_TO_C  ! index table between GRID_SOFT and GRIDC
    COMPLEX(q)  CHTOT(GRIDC%MPLWV,WDES%NCDIJ) ! charge-density in real / reciprocal space
    RGRID       DENCOR(GRIDC%RL%NP)           ! partial core
    COMPLEX(q)  CVTOT(GRIDC%MPLWV,WDES%NCDIJ) ! local potential
    COMPLEX(q)  CSTRF(GRIDC%MPLWV,T_INFO%NTYP)! structure factor
    RGRID   SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ) ! local potential
    COMPLEX(q)  CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
    INTEGER            IRDMAX, N_MIX_PAW
    REAL(q)  RHOLM(N_MIX_PAW,WDES%NCDIJ)
    TYPE (potcar)      P(T_INFO%NTYP)
    TYPE (energy)      E
    TYPE (in_struct)   IO
    TYPE (info_struct) INFO
    REAL(q)   XCSIF(3,3)                      ! stress stemming from XC
    INTEGER NK

    PROFILING_START('update_cdij')
  ! the first step is to restore the original XC-functional (as read from POTCAR/INCAR)
  ! in order to use the appropriate Hamiltonian
    CALL POP_XC_TYPE
    IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
       CALL SETUP_LDA_XC(2,-1,-1,IO%IDIOT)
    ELSE
       CALL SETUP_LDA_XC(1,-1,-1,IO%IDIOT)
    ENDIF
  ! now update the PAW one-center terms to current functional
    CALL SET_PAW_ATOM_POT( P , T_INFO, WDES%LOVERL, LMDIM, INFO%EALLAT, INFO%LMETAGGA, IO%IU6 )
    ! and restore the convergence corrections
    CALL SET_FSG_STORE(GRIDHF, LATT_CUR, WDES)

  ! just make sure that data distribution is over bands
    CALL REDIS_PW_OVER_BANDS(WDES, W)

    ! recalculate charge density and  kinetic energy density
    CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
         GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
         LATT_CUR, P, SYMM, T_INFO, &
         CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

    CALL SET_KINEDEN(GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM, &
         T_INFO%NIONS,W,WDES,KINEDEN)      
    CALL UPDATE_POT

    PROFILING_STOP('update_cdij')

  CONTAINS 

    SUBROUTINE UPDATE_POT
      USE pot
      USE pawm
      USE morbitalmag
      USE us
      REAL(q) :: XCSIF(3,3)
      INTEGER :: IRDMAA

      CALL POTLOK(GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES, &
           INFO,P,T_INFO,E,LATT_CUR, &
           CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF)
      
      CALL POTLOK_METAGGA(KINEDEN, &
           GRID,GRIDC,GRID_SOFT,WDES%COMM_INTER,WDES,INFO,P,T_INFO,E,LATT_CUR, &
           CHDEN,CHTOT,DENCOR,CVTOT,SV,HAMILTONIAN%MUTOT,HAMILTONIAN%MU,SOFT_TO_C,XCSIF)
      
      CALL VECTORPOT(GRID, GRIDC, GRID_SOFT, SOFT_TO_C,  WDES%COMM_INTER, & 
           LATT_CUR, T_INFO%POSION, HAMILTONIAN%AVEC, HAMILTONIAN%AVTOT)
      
      
      CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
           LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)
      
      CALL SETDIJ_AVEC(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
           LMDIM,CDIJ,HAMILTONIAN%AVTOT, NONLR_S, NONL_S, IRDMAX)
      
      CALL SET_DD_MAGATOM(WDES, T_INFO, P, LMDIM, CDIJ)
      
      CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
           WDES%NCDIJ, LMDIM, CDIJ(1,1,1,1),  RHOLM, CRHODE(1,1,1,1), &
           E,  LMETA=.FALSE., LASPH=INFO%LASPH, LCOREL= .FALSE.  )
      
      CALL UPDATE_CMBJ(GRIDC,T_INFO,LATT_CUR,IO%IU6)
      
    END SUBROUTINE UPDATE_POT

  END SUBROUTINE UPDATE_CDIJ

!***********************************************************************
!
! restore the original WDES, as well as the original
! number of orbitals
! this routine closely follows EDDIAG_EXACT_UPDATE_NBANDS
!
!***********************************************************************

  SUBROUTINE RESTORE_NBANDS( WDES, W, LATT_CUR, INFO, LMDIM, T_INFO, P, IO)
    USE wave_high
    USE fock
    USE lattice
    USE setexm
    USE scala
    USE pseudo
    USE poscar
    USE pawm
    IMPLICIT NONE
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    TYPE (latt)        LATT_CUR
    TYPE (info_struct) INFO
    INTEGER            LMDIM
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    TYPE (in_struct)   IO
  ! local variables
    TYPE (wavespin)    W_TMP
    INTEGER NK

    IF (.NOT. LINIT_RPA_FORCE) RETURN

    WDES%NB_TOT=WDES_GROUNDSTATE%NB_TOT
    WDES%NBANDS=WDES_GROUNDSTATE%NBANDS

    CALL INIT_SCALAAWARE( WDES%NB_TOT, WDES%NRPLWV, WDES%COMM_KIN )

    DO NK=1,WDES%NKPTS
       IF (WDES%LGAMMA) THEN
          WDES%NB_TOTK(NK,:)=MIN(WDES%NB_TOT,WDES%NPLWKP_TOT(NK)*2-1)
       ELSE
          WDES%NB_TOTK(NK,:)=MIN(WDES%NB_TOT,WDES%NPLWKP_TOT(NK))
       ENDIF
    ENDDO

    CALL RESETUP_FOCK_WDES(WDES, LATT_CUR, LATT_CUR, -1)

    CALL ALLOCW(WDES, W_TMP)

    W_TMP%CW(:,1:WDES%NBANDS,:,:)    =W%CW(:,1:WDES%NBANDS,:,:)
    W_TMP%CPROJ(:,1:WDES%NBANDS,:,:) =W%CPROJ(:,1:WDES%NBANDS,:,:)
    W_TMP%CELTOT(1:WDES%NB_TOT,:,:)=W%CELTOT(1:WDES%NB_TOT,:,:)
    W_TMP%FERTOT(1:WDES%NB_TOT,:,:)=W%FERTOT(1:WDES%NB_TOT,:,:)

    CALL DEALLOCW(W)
    W=W_TMP

!    CALL POP_XC_TYPE
!    IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
!       CALL SETUP_LDA_XC(2,IO%IU6,IO%IU0,IO%IDIOT)
!    ELSE
!       CALL SETUP_LDA_XC(1,IO%IU6,IO%IU0,IO%IDIOT)
!    ENDIF
!    AEXX    =1.0-LDAX
!    HFSCREEN=LDASCREEN
  ! now update the PAW one-center terms to current functional
!    CALL SET_PAW_ATOM_POT( P , T_INFO, WDES%LOVERL, LMDIM, INFO%EALLAT, INFO%LMETAGGA, IO%IU6 )

    ! and restore the convergence corrections
    CALL SET_FSG_STORE(GRIDHF, LATT_CUR, WDES)
    
  END SUBROUTINE RESTORE_NBANDS

END MODULE
