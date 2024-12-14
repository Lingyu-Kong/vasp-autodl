#include "symbol.inc"

!*********************************************************************
!
! this module implements the calculation of the response/NATURAL
! function from the Green's function
!   G(r',r, tau ) G* + (r',r, tau)
! where G is the Green's function for the unoccupied states
!   G (r',r, tau) = sum_a <r'|a> <a| r> exp ( -(e_a- e_F) tau )
! and G+ is the Green's function for the occupied states
!   G+(r',r, tau) = sum_a <r'|i> <i| r> exp (  (e_i- e_F) tau )
! 
! The k-dependent code is selected with ALGO = RPARK, GWRK, ...
! and uses convolution in k-space resulting and is quite slow
! compared to chi_super. However, less memory is required 
!
! note that the resulting response functions are conjugated
! compared to the definition used in chi_base.F
!
!*********************************************************************

MODULE chi_GG
  USE prec
  USE greens_real_space
  USE GG_base
  USE greens_orbital
  USE rpa_force
  USE chi_glb 
  USE lt_mp2      
  USE stoch_lt_mp2 
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: CALCULATE_XI_REAL_GWRK

CONTAINS


!***********************************************************************
!
! main calling routine
! TODO: this is presently a simple copy of CALCULATE_XI_REAL
!   obviously it might be possible to remove some arguments
!   although this is not very likely
!
!
! TODO gK: very few routine handle multiple time points at a time
!  but there are so few, that the support is really broken, 
!  and it would be cleaner to remove it alltogehter
!  all Greensfunction like arrays are allocated for multiple TAU
!  most likely best to change this top down
!
!
! notes on KPAR parallelization
! step 1: all WDES%COMM replaced by WDES%COMM_KIN
!         WDES%COMM_INTER equals WDES%COMM_KIN if 
!         parallelization over orbitals is not selected 
!
!***********************************************************************

  SUBROUTINE CALCULATE_XI_REAL_GWRK( & 
          HAMILTONIAN, P, WDES, NONLR_S, NONL_S, W, LATT_CUR, LATT_INI, &
          GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, E, & 
          CHTOT, CHTOTL, DENCOR, CVTOT, CSTRF, IRDMAX, &
          T_INFO, DYN, INFO, IO, KPOINTS, SYMM, MIX, & 
          LMDIM, CQIJ, CDIJ, CRHODE, N_MIX_PAW, RHOLM, RHOLM_LAST, CHDEN, SV, & 
          EFERMI_IN, TOTEN, NEDOS, DOS, DOSI )
    USE base
    USE pseudo
    USE mpimy
    USE constant
    USE poscar
    USE wave
    USE wave_high
    USE ini
    USE fileio, ONLY: OUTWAV
    USE hamil_struct_def
    USE subrot
    USE acfdt_GG
    USE minimax_struct, ONLY : imag_grid_handle
    USE minimax
    USE pade_fit
    USE dmft
    USE pead, ONLY : PEAD_RESETUP_WDES
    USE mymath, ONLY : GREATEST_COMMON_DIVISOR
    USE tutor, ONLY: vtutor, isError, isAlert, ACFDTRnotSupported, FermiFiniteTemp, &
        OmegaGrid140, TauParCPU, OmegaPar, TauPar, &
        chiGG_NOmega, OmegaGridOmegaTL, argument, & 
        ! OEP related
        LPOTOK  
    ! required for OEP
    USE us, ONLY: SET_CHARGE
    USE mlr_optic, ONLY: LR_OPTIC_DONE
    USE density_of_states, ONLY: DENSTA_DER
#ifdef VASP_HDF5
    USE vhdf5
#endif
    IMPLICIT NONE
! structures
    TYPE (ham_handle)  HAMILTONIAN
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    TYPE (wavedes)     WDES
    TYPE (nonlr_struct)NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W
    TYPE (latt)        LATT_CUR, LATT_INI
    TYPE (grid_3d)     GRID       ! grid for wavefunctions
    TYPE (grid_3d)     GRID_SOFT  ! grid for soft chargedensity
    TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
    TYPE (grid_3d)     GRIDUS     ! temporary grid in us.F
    TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
    TYPE (transit)     SOFT_TO_C  ! index table between GRID_SOFT and GRIDC
    TYPE (energy)      E
    COMPLEX(q) CHTOT(GRIDC%MPLWV,WDES%NCDIJ) ! charge density
    COMPLEX(q) CHTOTL(GRIDC%MPLWV,WDES%NCDIJ)! old charge-density
    RGRID      DENCOR(GRIDC%RL%NP)
    COMPLEX(q)  CVTOT(GRIDC%MPLWV,WDES%NCDIJ) ! local potential
    COMPLEX(q)  CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
    INTEGER     IRDMAX
    TYPE (dynamics)    DYN
    TYPE (info_struct) INFO
    TYPE (in_struct)   IO
    TYPE (kpoints_struct) KPOINTS
    TYPE (symmetry)    SYMM
    TYPE (mixing)      MIX
    INTEGER  LMDIM
    OVERLAP  CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    OVERLAP  CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
!  paw sphere charge density
    INTEGER N_MIX_PAW
    REAL(q)  RHOLM(N_MIX_PAW,WDES%NCDIJ),RHOLM_LAST(N_MIX_PAW,WDES%NCDIJ)
!  charge-density and potential on soft grid
    COMPLEX(q)  CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
    RGRID       SV(DIMREAL(WDES%GRID%MPLWV),WDES%NCDIJ)
    REAL(q)     EFERMI_IN
    REAL(q)     EFERMI(WDES%ISPIN)
    REAL(q)     TOTEN
    INTEGER     NEDOS
    REAL(q)     DOS(NEDOS,WDES%ISPIN),DOSI(NEDOS,WDES%ISPIN)
#ifdef VASP_HDF5
    INTEGER IH5ERR
#endif
#ifndef scaLAPACK
    CALL vtutor%write(isError, ACFDTRnotSupported)
  END SUBROUTINE CALCULATE_XI_REAL_GWRK
#else

! local variables
  ! we need two descriptors, one for the time (tau) and one for the energy-frequency (nu) domain
  ! these two are required since in time domain the data are distributed over all cores
  ! and one imaginary time is usually done at time
  ! whereas in frequency domain a subset of cores can work on one frequency 
  ! this improves parallel efficiency but keeps the storage demand during the
  ! most memory intensive step at bay
    TYPE (greensfdes) :: GDES_TAU              ! Green function in time  domain
    TYPE (greensfdes) :: GDES                  ! Green function in freq. domain

    TYPE (correlation), POINTER :: COR
    LOGICAL LGAMMA
    TYPE (responsefunction), POINTER :: CHI_TAU(:)  ! response function in tau domain
    TYPE (responsefunction), POINTER :: CHI         ! response function in nu domain
    TYPE (responsefunction), POINTER :: CHI_QPT(:)  ! response function in nu domain
    TYPE (responsefunction), POINTER :: CHI0(:)  ! response function at w=0
    !!!need to figure out which one we need in the WORK state
    
    TYPE (greensf), POINTER :: GU(:)           ! Green function for unoccupied orbitals in imag. time  G( tau)
    TYPE (greensf), POINTER :: GO(:)           ! Green function for occupied orbitals in imag. time G(-tau)
    TYPE (greensf), POINTER :: SIGMAU(:)       ! Self-energy G( tau) W(tau)  (invokes only GU)
    TYPE (greensf), POINTER :: SIGMAO(:)       ! Self-energy G(-tau) W(tau)  (invokes only GO)

    TYPE (wavedes), POINTER :: WGW
    TYPE (wavedes1)         :: WGWQ
    INTEGER                 :: NQ, NQ_COUNTER, NK1, NK2, ISP
    TYPE (grid_3d), POINTER :: GRIDWGW
    TYPE (screened_2e_handle) :: S2E

    TYPE (imag_grid_handle) :: IMAG_GRIDS      ! imaginary grids handle
    TYPE (wavespin) WHF
    INTEGER :: N1, N2, NCPU, ierror
    REAL(q)    PAR(1,1,1,1,WDES%NCDIJ),DOSPAR(1,1,1,WDES%NCDIJ)
    GDEF, POINTER       :: SIGMAO_DIAG(:)      ! diagonal of <i| Sigma_o(tau) | i> of self-energy
    GDEF, POINTER       :: SIGMAU_DIAG(:)      ! diagonal of <i| Sigma_u(tau) | i> of self-energy
    GDEF, POINTER, CONTIGUOUS :: SIGMAO_MAT(:,:,:,:) ! matrix <i| Sigma_o(tau) | a> stored distributed
    GDEF, POINTER, CONTIGUOUS :: SIGMAU_MAT(:,:,:,:) ! matrix <i| Sigma_u(tau) | a> stored distributed
    GDEF, POINTER, CONTIGUOUS :: CHAM_MAT(:,:,:,:)   ! Hamilton matrix T + V_local + V_x or density matrix
    GDEF, POINTER, CONTIGUOUS :: GAMMA(:,:,:,:)      ! density matrix for one k-point
    GDEF, POINTER, CONTIGUOUS :: CORR_MAT(:,:,:,:)   ! correlation contribution to density matrix
    GDEF, POINTER, CONTIGUOUS :: SIGMA_MAT(:,:,:,:)=>NULL()  ! int G_0(w) (\Sigma(w)+H_HF) G_0(w) w dw
    GDEF, POINTER, CONTIGUOUS :: GO_MAT(:,:,:,:)     ! matrix <i| G_o(tau) | a> stored distributed
    GDEF, POINTER, CONTIGUOUS :: GU_MAT(:,:,:,:)     ! matrix <i| G_u(tau) | a> stored distributed
    COMPLEX(q), POINTER, CONTIGUOUS :: SIGMAW_MAT(:,:,:,:) ! frequency dependent self-energy 
    COMPLEX(q), POINTER :: SIGMA_COS(:,:,:,:)  ! cosine transformed self-energy
    COMPLEX(q), POINTER :: SIGMA_SIN(:,:,:,:)  ! sin transformed self-energy
    GDEF, POINTER       :: SIGMAO_TAU(:,:,:,:)=>NULL() ! diagonal of <i| Sigma_o(tau) | i> of self-energy
    GDEF, POINTER       :: SIGMAU_TAU(:,:,:,:)=>NULL() ! diagonal of <i| Sigma_u(tau) | i> of self-energy
    INTEGER         :: NTAU_ROOT
    INTEGER :: NOMEGA_GRID                     ! number of frequency points in grids
    INTEGER :: I, J
    INTEGER :: IU6_TMP, IU6_MEM
    REAL(q) :: TAU_TMP                            ! auxillary: temporary storage
                                               ! for local time point
    !variables for head and wings
    TYPE( head_handle) :: HEADHAND             ! handle for the long-wave limit of chi
    CHARACTER(LEN=16):: TEXT
    TYPE (skpoints_trans) ::  RES_KPTS_TRANS
                                               ! handle to generate response function at k-points outside the IRZ
    INTEGER :: NQ_IRZ                          ! q-point in IRZ
    INTEGER :: NKPTS_IRZ                       ! total number of k-points in the IRZ
    TYPE (wavespin)  WMEAN                     ! mean field Hamiltonian
                                               ! identical to W, except for CELTOT and FERTOT
    TYPE (wavespin)  WHF_IN_KS                 ! diagonals of HF Hamiltonian in KS basis
    REAL(q) :: EFERMI_MEAN(WDES%ISPIN)         ! Fermi level for mean field Hamiltonian
    REAL(q) :: EFERMI_SAVE(WDES%ISPIN)         ! Fermi level for mean field Hamiltonian
    REAL(q) :: EFERMI_VAR                      ! fermi energy variance 
    REAL(q) :: EFERMI_SLOPE(2)                 ! slope of Fermi level
   ! convergence correction for screened Coulomb potential and self-energy
    REAL(q) :: FSG0
    COMPLEX(q):: MP2SUM, MP2SUM2, CSINGLES
    COMPLEX(q):: E_MP2
    REAL(q) :: XCSIF(3,3)         ! temporary (stress tensor)
    INTEGER :: IRDMAA             ! temporary
    INTEGER :: NELM               ! loop counter in GW part
    INTEGER :: NELM_HF            ! loop counter in HF part
    GDEF    :: R(WDES%NB_TOT)
    REAL(q) :: RR(WDES%NB_TOT)
    TYPE (greens_mat_des), POINTER :: GDES_MAT=>NULL()
    REAL(q) :: FRNL(3,T_INFO%NIONS)
    REAL(q)    TOTEGM, DTRGHF, FLWTOT_OLD
    LOGICAL, EXTERNAL :: CALCULATE_RPA_FORCES
    COMPLEX(q):: TR_GAMMA(WDES%ISPIN)
    REAL(q) ::  FERDER( W%WDES%NB_TOT, W%WDES%NKPTS, W%WDES%ISPIN )
    ! correlated wavefunctions for CRPA calculations
    TYPE (wavespin)                  :: WCORR
    TYPE (wavespin)                  :: WCORRHF 
    ! OEP related 
    TYPE( oep_handle ), POINTER      :: OEP => NULL()
    LOGICAL :: LKINETIC=.FALSE. ! same as in chi_base.F for OEP
!    LOGICAL :: LKINETIC=.TRUE.
! not used, but required as dummy argument for polarizability
    TYPE(supercell),POINTER :: SUPER=>NULL() 
#ifdef MPI
    IF (W%WDES%COMM_KINTER%NCPU.NE.1 .AND. (.NOT. LMP2LT ) .AND. (.NOT. LSMP2LT) ) THEN
       CALL vtutor%error("CALCULATE_XI_REAL_ACFDTRK: KPAR>1 not implemented, sorry.")
    END IF
#endif
    PROFILING_START('calculate_xi_real_gwrk')
    ! hand brake if user tries to use  non-collinear version 
    IF ( WDES%LNONCOLLINEAR ) THEN
       CALL vtutor%error("non collinear calculations for low-&
          &scaling GW-RPA algorithms are not implemented yet")
    ENDIF

! finite temperature NU_DOWN not implemented yet
    IF (INFO%NUP_DOWN>0 .AND. LFINITE_TEMPERATURE ) THEN
       CALL vtutor%error("T>0 GW algorithms for NUP_DOWN > 0 not implmeneted, sorry")
    ENDIF
!pull breaks here      
!finite temperature RPA and GW work with fermi smearing only
    IF ( LFINITE_TEMPERATURE .AND.  KPOINTS%ISMEAR /= -1 ) THEN
       CALL vtutor%write(isError, FermiFiniteTemp)
    ENDIF
!
! safety check done 
!
    IF( IO%IU0>=0 .AND. LDMP1) WRITE(*,*)'Linear term in ACFDT is being calculated'

    IF (W%WDES%NKPTS == 1 .AND. SUM(W%WDES%VKPT(:,1)*W%WDES%VKPT(:,1))<G2ZERO)  THEN
       LGAMMA = .TRUE.
    ELSE
       LGAMMA = .FALSE.
    ENDIF

    CALL SET_NBANDSGW( W)

    !fermi energy update might be required    
    CALL DENSTA_SPIN( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
       INFO%NUP_DOWN,  E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
       NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)

    ! set W%AUXTOT to the derivatives of the fermi weights
    ! use FERDER instead, since AUXTOT is used in chi_base.F functions
    ! breaks compuation of HEAD and wings of CHI
    ! also W%EFERMI must be set 
    CALL DENSTA_DER( W, KPOINTS, FERDER)
      
    IF ( LFINITE_TEMPERATURE ) THEN
!       all set up properly
    ELSE
       CALL MEAN_CBM_VBM_SPIN( W, EFERMI, INFO%NUP_DOWN)

       IF (IO%IU6>=0) THEN
          WRITE(IO%IU6,"(/,' The Fermi energy was updated, please check that it is located mid-gap',/ &
          & ' values below the HOMO (VB) or above the LUMO (CB) will cause erroneous energies',/ &
          & ' E-fermi : ', 2F8.4,/)" ) W%EFERMI
          WRITE(IO%IU0,"(/,' The Fermi energy was updated, please check that it is located mid-gap',/ &
          & ' values below the HOMO (VB) or above the LUMO (CB) will cause erroneous energies',/ &
          & ' E-fermi : ', 2F8.4,/)" ) W%EFERMI
       ENDIF
    ENDIF

    IF( IO%IU0>=0 .AND. IO%NWRITE > 2 ) WRITE(*,*)'EFERMI=',W%EFERMI

    IF (IO%LOPTICS .AND. .NOT. LR_OPTIC_DONE()) THEN
       CALL START_TIMING("G")
       CALL PEAD_RESETUP_WDES(WDES, GRID, KPOINTS, LATT_CUR, LATT_CUR, IO)
       CALL LR_OPTIC( &
            P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
            T_INFO,INFO,IO,KPOINTS,SYMM,GRID,GRID_SOFT, &
            GRIDC,GRIDUS,C_TO_US,SOFT_TO_C, &
            CHTOT,DENCOR,CVTOT,CSTRF, &
            CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM, &
            CHDEN,SV,LMDIM,IRDMAX,EFERMI_IN,NEDOS, & 
            LSTORE=.TRUE., LPOT=.FALSE.)
       CALL STOP_TIMING("G",IO%IU6,'OPTICS')
    ELSE IF (.NOT. LRPAFORCE .AND. .NOT. LR_OPTIC_DONE()) THEN
       CALL READ_CDER_BETWEEN_STATES(WDES, IO%IU0, 55) ! read first derivative of wavefuntions
    ENDIF

! autoset available memory, incase MAXMEM if not set by user in INCAR 
    IF ( LAUTOSET_MAXMEM ) THEN
       CALL AUTOSET_AVAILABLE_MEMORY( MAXMEM, IO )
    ENDIF 
!
! sync orbitals and eigenvalues (just to be sure everything is ok)
    CALL KPAR_SYNC_ALL(WDES,W)
!=======================================================================
! switch off symmetry
!=======================================================================
    CALL NO_SYMMETRY(.FALSE., LGAMMA, SYMM, W, WDES, LATT_CUR, LATT_INI, &
     T_INFO, INFO, KPOINTS, GRID, NONL_S, P, WHF, IO )

!=======================================================================
! constrained RPA setup, after symmetry switched off 
!=======================================================================
    IF ( LCRPA ) THEN 
       CALL CRPA_SETUP( WDES, W, WCORR, LMDIM, T_INFO, P, CQIJ,&
          LATT_CUR, GRID, INFO, KPOINTS, S2E%NUMBER_OF_NQ, IO )  

       CALL TAG_SELF_CORRELATION_STATES( WDES, W, IO )
        
       !for disentangle and standard CRPA calculations an additional wavefunction is needed 
       IF ( LCHIC .AND. .NOT. LWEIGHTED) THEN
          WCORRHF=WCORR
          WCORRHF%WDES => WDES_FOCK
       ENDIF
   ENDIF

   CALL DUMP_ALLOCATE(IO%IU6)
!=======================================================================
! generate descriptor for response function
!=======================================================================
    ALLOCATE(WGW, GRIDWGW)
    WGW=WDES_FOCK

    WGW%NKPTS=KPOINTS_FULL%NKPTS
    WGW%NKDIM=KPOINTS_FULL%NKPTS
    WGW%NKPTS_FOR_GEN_LAYOUT=KPOINTS_FULL%NKPTS
! KPOINTS_FULL structure might be reallocated better to allocate and copy data
    ALLOCATE(WGW%VKPT(1:3,SIZE(KPOINTS_FULL%VKPT,2)),WGW%WTKPT(SIZE(KPOINTS_FULL%WTKPT,1)))
    WGW%VKPT =KPOINTS_FULL%VKPT
    WGW%WTKPT=KPOINTS_FULL%WTKPT
    WGW%ENMAX=ENCUTGW
    IF (ENCUTLF==-1) ENCUTLF=WGW%ENMAX

    ! GRIDWGW is identical to GRID_FOCK, except for GRIDWGW%FFTSCA
    GRIDWGW=GRID_FOCK
    IF (IO%IU6>=0) THEN 
       WRITE(IO%IU6,*) 'Basis sets for responsefunctions:'
       WRITE(IO%IU6,*) '================================='
    ENDIF
    CALL GEN_LAYOUT(GRIDWGW, WGW, LATT_CUR%B, LATT_CUR%B, IO%IU6,.TRUE.)
    IF (WGW%LGAMMA) THEN
       ! gamma only data layout with response functions stored as real in real space
       GRIDWGW%LREAL=.TRUE.
    ENDIF

    CALL GEN_INDEX (GRIDWGW, WGW, LATT_CUR%B, LATT_CUR%B,IO%IU6,-1, .TRUE.)
    !  init FFT (required if real to complex FFT is used)
    CALL FFTINI(WGW%NINDPW(1,1), WGW%NGVECTOR(1), WGW%NKPTS, WGW%NGDIM, GRIDWGW)

    IF(IO%IU6>=0) WRITE(IO%IU0,'(A,I8)') ' responsefunction array rank=',WGW%NGDIM
    !
    ! some safety measures for OEP
    !
    IF ( LOEP ) THEN
       IF ( (GRID_SOFT%NGX /= GRIDWGW%NGX .OR. GRID_SOFT%NGY /= GRIDWGW%NGY .OR. &
           GRID_SOFT%NGZ /= GRIDWGW%NGZ ) ) THEN
           CALL vtutor%error("currently the OEP method is supported only for PRECFOCK = Normal")
       ENDIF
       ! POT file must be present for OEP 
       IF (.NOT. INFO%LPOTOK) THEN
          CALL vtutor%write(isError, LPOTOK)
       ENDIF
    ENDIF

    ! replace kinetic energy with the Coulomb potential in the WGW%DATAKE array
    IF (ENCUTGW /= ENCUTGWSOFT .AND. ENCUTGW > 0 .AND. ENCUTGWSOFT >0 ) THEN
       CALL XI_DATAKE( WGW, LATT_CUR, ENCUTGW, ENCUTGWSOFT)
    ELSE
       CALL XI_DATAKE( WGW, LATT_CUR)
    ENDIF

    ! determine the q-points in IrBZ that are used for the response function
    CALL SETUP_IRZ_MAP(S2E, WGW, IO%IU0, IO%IU6)
    ! set the number of k-points in the original symmetry reduced irreducible BZ
    ! ? most likely the total number of k-points or?
    NKPTS_IRZ=KPOINTS_ORIG%NKPTS
! change S2E to cover entire BZ, usefull for some tests
!    S2E%NUMBER_OF_NQ=WDES%NKPTS
!    DEALLOCATE(S2E%NQ); ALLOCATE(S2E%NQ(S2E%NUMBER_OF_NQ)); S2E%NQ=(/ ( N1,N1=1, S2E%NUMBER_OF_NQ) /)

    ! transformation between IrBZ and BZ
    IF (SYMM%ISYM>=0) THEN
       CALL GENERATE_KPOINTS_TRANS(WGW%GRID, NKPTS_IRZ, WGW, KPOINTS_FULL_ORIG, RES_KPTS_TRANS)
    ENDIF
    ! relink W the WHF so that the fft grids of the Hartree Fock routine are used
    WHF=W
    WHF%WDES => WDES_FOCK
    CALL CHECK_FULL_KPOINTS

!=======================================================================
! determine frequency grid and complex shift (literal copy from chi.F)
!=======================================================================
    IMAG_GRIDS%T%NPOINTS=NOMEGA
    ! set OMEGATL
    IF (OMEGATL<0 .AND. NODES_IN_DIELECTRIC_FUNCTION>=0) THEN
       OMEGATL=NODES_IN_DIELECTRIC_FUNCTION*10
       IF ( (LACFDT .OR. LG0W0 .OR. LCRPA) .AND. (OMEGAGRID>= 100)) THEN
          OMEGATL=MAX_ENERGY_UNOCCUPIED(WDES,W)
       ELSE IF (LACFDT .OR. LCRPA) THEN
          OMEGATL=MAX(OMEGATL, MAX_ENERGY_UNOCCUPIED(WDES,W)*2.0)
       ENDIF
    ELSEIF (OMEGATL<0) THEN
       IF ( (LACFDT.OR.LG0W0 .OR. LCRPA ) .AND. (OMEGAGRID>= 100)) THEN
          OMEGATL=MAX_ENERGY_UNOCCUPIED(WDES,W)
       ELSE IF (LACFDT .OR.LG0W0 .OR. LCRPA) THEN
          OMEGATL=MAX(OMEGATL, MAX_ENERGY_UNOCCUPIED(WDES,W)*2.0)
       ENDIF
    ENDIF
    IF (OMEGATL<0) THEN
       OMEGATL=ABS(OMEGATL)
    ENDIF
    IF ( .NOT. (LACFDT .OR. OMEGAGRID<100 .OR. LCRPA) ) THEN
    ! final check on OMEGATL: the frequency grid MUST go at least to 
    ! the maximum transition energy (safeguarded by a factor 1.1)
       OMEGATL=MAX(OMEGATL, MAX_ENERGY_UNOCCUPIED(WDES,W)*1.1)
    ENDIF

    ! position of HOMO
    N1=MIN(MAX(LAST_FILLED_OPTICS_NO_MOD(W),LAST_FILLED_OPTICS_NO_MOD(W)),WDES%NB_TOT)

    ! set OMEGAMIN (usually the band gap)
    IF (OMEGAMIN<0) THEN
       ! set internal CRPA variables, this might change minimization interval
       IF (INFO%LCRPA ) &
          CALL INIT_CRPA_STATES( WDES%NB_TOT, IO%IU0, IO%IU6 )
       CALL DETERMINE_BAND_GAP_SPIN(WDES, W, WDES%NB_TOT, OMEGAMIN, OMEGATL, NOMEGA)
    ENDIF

    ! set OMEGAMAX
    IF (OMEGAMAX<0 .AND. NOMEGA <= 1) THEN
       OMEGAMAX=0
    ELSE IF (OMEGAMAX<0 .AND. NODES_IN_DIELECTRIC_FUNCTION>=0) THEN
       OMEGAMAX=NODES_IN_DIELECTRIC_FUNCTION/1.3
       ! bottom of valence to top of valence (HOMO) = valence band width
       OMEGAMAX=MAX(OMEGAMAX, MAX_ENERGY_OCC_UNOCCUPIED(WDES,W, N1)/1.3)
    ENDIF
    IF (OMEGAMAX<0) THEN
       ! bottom of valence to top of valence (HOMO) = valence band width
       OMEGAMAX=ABS(OMEGAMAX)
       OMEGAMAX=MAX(OMEGAMAX, MAX_ENERGY_OCC_UNOCCUPIED(WDES,W, N1)/1.3)
    ENDIF

    IF (IO%IU6>=0) THEN
       WRITE(IO%IU6,100) NOMEGA, OMEGAMAX, OMEGAMIN, OMEGATL
    ENDIF

100 FORMAT(' Response functions by GG contraction:'  / &
         ' ==============================================='  / &
         '   NOMEGA  =',I6,  '    number of frequencies'/ &
         '   OMEGAMAX=',F7.2,'   maximum frequency'/ &
         '   OMEGAMIN=',F7.2,'   minimum frequency'/ &
         '   OMEGATL =',F7.2,'   maximum frequency of tail'/ &
         )
      
    !we need either minimax or frequency grid + fourier matrix
    IF ( OMEGAGRID< 140 ) THEN
        CALL vtutor%write(isError, OmegaGrid140)
    ENDIF

    IF (LHFCALC_GG) THEN
       ! note that this is not needed for Finite-temperature 
       ! since 0 point is already contained in grid
       NOMEGA_GRID=NOMEGA-1
    ELSE
       NOMEGA_GRID=NOMEGA
    ENDIF

!=======================================================================
! we split up the original communicator into groups
! the frequencies are distributed between the groups, each group 
! receives IMAG_GRIDS%B%NPOINTS_IN_GROUP
! frequencies (ideally one) and works only on those
! make some sanity check before determining the grids 
!=======================================================================
    ! safety measure to prevent silly settings of user
    IF ( NOMEGAPAR>WDES%COMM_KIN%NCPU .OR. NOMEGAPAR>NOMEGA ) THEN
       !if user sets NOMEGAPAR too large, the number of omega groups is set to 
       !total number of CPUs or NOMEGAPAR
       NOMEGAPAR = MIN( NOMEGAPAR, WDES%COMM_KIN%NCPU )    
       NOMEGAPAR = MIN( NOMEGAPAR, NOMEGA )
       CALL vtutor%write(isAlert, OmegaPar, argument( ival=[NOMEGAPAR]) )
    ENDIF
    IF (NTAUPAR>WDES%COMM_KIN%NCPU .OR. NTAUPAR>NOMEGA ) THEN
       !if user sets NTAUPAR too large, the number of tau groups is set to 
       !total number of CPUs or NTAUPAR
       NTAUPAR = MIN( NTAUPAR,WDES%COMM_KIN%NCPU)
       NTAUPAR = MIN( NTAUPAR, NOMEGA )
       CALL vtutor%write(isAlert, TauPar, argument( ival=[NTAUPAR]) )
    ENDIF
    ! possible autoset of NOMEGAPAR and NTAUPAR is required, this is based on
    ! available memory per rank given by MAXMEM 
    IF (NOMEGAPAR < 0 .OR. NTAUPAR < 0) THEN
       CALL SET_NTAUPAR_NOMEGAPAR( WDES_FOCK, WGW, S2E, T_INFO, KPOINTS, LATT_CUR, LMDIM, P, IO)
    ENDIF

    IF ( MOD(WDES%COMM_KIN%NCPU,NOMEGAPAR)/=0 .OR. MOD(WDES%COMM_KIN%NCPU,NTAUPAR)/=0 ) THEN
       !NOMEGAPAR and NTAUPAR defines a 2D processor grid, 
       !so the total number of CPUs has to be dividable by these numbers
       CALL vtutor%write(isError, TauParCPU, &
         argument( ival=[WDES%COMM_KIN%NCPU,NTAUPAR, NOMEGAPAR ] ) )
    ENDIF
#ifdef MPI
    ! one last check for NTAUPAR
    CALL CHECK_GDES_MAT_SIZE( WDES%COMM_KIN, NTAUPAR, WDES%NB_TOT, IO ) 
#endif 

    ! estimate grid
    IF (NOMEGA_GRID==1) THEN
       CALL SET_IMAG_GRID_HANDLE(OMEGAMIN, OMEGATL, NOMEGA, OMEGAGRID, &
           KPOINTS%SIGMA, IMAG_GRIDS, WDES%COMM, IO)
       CALL ALLOCATE_IMAG_GRID( IMAG_GRIDS, NOMEGA )
       IMAG_GRIDS%TAU=EPSILON(1._q)
       IMAG_GRIDS%TAU_WEIGHT=1.0
       IMAG_GRIDS%BOS_RE_WEIGHT=1.0
       IF ( IO%IU0>=0 ) WRITE(IO%IU0, '(A,E10.3)' ) ' used tau point:', IMAG_GRIDS%TAU(1)
       IF ( IO%IU6>=0 ) WRITE(IO%IU6, '(A,E10.3)' ) ' used tau point:', IMAG_GRIDS%TAU(1)
    ELSE IF (NOMEGA_GRID<2) THEN
       CALL vtutor%write(isError, chiGG_NOmega)
    ELSE
       ! determine grid
       CALL SET_IMAG_GRID_HANDLE(OMEGAMIN, OMEGATL, NOMEGA, OMEGAGRID, &
           KPOINTS%SIGMA, IMAG_GRIDS, WDES%COMM, IO)

       !RPA needs cos + exp grid only
       IF ( LACFDT ) THEN
          ! for laplace transformed mp2
          ! it will determine the direct MP2 correlation energy, but yield ynonsence for RPA correlation energy
          IF ( LLTDMP2 ) THEN
             CALL ALLOCATE_IMAG_GRID( IMAG_GRIDS, NOMEGA )
             IMAG_GRIDS%BOS_RE=IMAG_GRIDS%TAU
             IMAG_GRIDS%BOS_RE_WEIGHT=IMAG_GRIDS%TAU_WEIGHT*2
             IMAG_GRIDS%TO_BOS_RE=0
             DO I = 1, IMAG_GRIDS%NOMEGA
                IMAG_GRIDS%TO_BOS_RE(I,I)=1
             ENDDO
          ENDIF 
          
       ENDIF 

       IF (IMAG_GRIDS%BOS_RE(NOMEGA_GRID)<OMEGATL/1.1 .AND. .NOT. LLTDMP2) THEN
          CALL vtutor%write(isAlert, OmegaGridOmegaTL)
       ENDIF
    ENDIF

    IF (LHFCALC_GG) THEN
       IF (IO%IU0>0) WRITE(IO%IU0,*) 'internal WARNING in VASP: adding time point, t=0 to frequency grid'
       ! note that this is not needed for Finite-temperature 
       ! since 0 point is already contained in grid
       CALL ADD_POINT_IMAG_GRID_HANDLE(IMAG_GRIDS)
    ENDIF

    IF (IO%LOPEN) CALL WFORCE(IO%IU6)
    CALLMPI( M_barrier(WDES%COMM))

    IF (IO%LOPEN) CALL WFORCE(IO%IU6)

!=======================================================================
! distribute grid points based on NOMEGAPAR and NTAUPAR 
!=======================================================================
    CALL SET_LOOP_DES( WDES%COMM_KIN, IMAG_GRIDS, 'boso', IO )
    CALL SET_LOOP_DES( WDES%COMM_KIN, IMAG_GRIDS, 'ferm', IO )
    CALL SET_LOOP_DES( WDES%COMM_KIN, IMAG_GRIDS, 'time', IO )

!=======================================================================
!  somewhat of a hack call (stochastic) Laplace transformed MP2 here
!=======================================================================
#ifdef MPI
    IF (LMP2LT) THEN
       CALL LAPLACE_TRANSFORMED_MP2( LATT_CUR, W, LMDIM, CQIJ, NONLR_S, NONL_S, P, IMAG_GRIDS%TAU, &
                                     IMAG_GRIDS%TAU_WEIGHT, WGW, ENCUTGW, ENCUTGWSOFT, IO%IU0, IO%IU6)
    ELSEIF (LSMP2LT) THEN
       CALL RND_ORBITALS_MP2( LATT_CUR, W, LMDIM, CQIJ, NONLR_S, NONL_S, P, IMAG_GRIDS%TAU, &
                              IMAG_GRIDS%TAU_WEIGHT, WGW, ENCUTGW, ENCUTGWSOFT, IO%IU0, IO%IU5, IO%IU6)
    ENDIF
#endif
!=======================================================================
! set two green function descriptors, one for time and one for frequency domain
! G(tau) is distributed among all nodes column wise. The same distribution is chosen
! for the response function in the imaginary time domain. 
! G(nu) is distributed among all nodes in groups, each group treating one freq. point
! separately. The same distribution is used for CHI in the frequency domain.
!=======================================================================
   CALL GDES_SET(GDES_TAU, WHF%WDES, WGW, T_INFO, P, IMAG_GRIDS%T%COMM_IN_GROUP)
   CALL GDES_SET(GDES, WHF%WDES, WGW, T_INFO, P, IMAG_GRIDS%B%COMM_IN_GROUP)

!  write info about required memory
   CALL DUMP_STORAGE_REQUIREMENTS( GDES_TAU, GDES, S2E, IMAG_GRIDS, WDES_FOCK, WGW, SUPER, IO )

   ALLOCATE(CHI_QPT(S2E%NUMBER_OF_NQ) ) 
   !allocate memory for CHI(nu) already here.
   DO NQ_COUNTER=1,S2E%NUMBER_OF_NQ  !  loop over q-points in the IRZ
      NQ=S2E%NQ(NQ_COUNTER)
      CALL SET_RESPONSE_KPOINT(CHI_QPT(NQ_COUNTER), WDES%VKPT(:,NQ), NQ)
      CHI=>CHI_QPT(NQ_COUNTER)
      CALL ALLOCATE_RESPONSEFUN_DISTRI( CHI, GDES%RES_NRPLWV_ROW, GDES%RES_NRPLWV_COL, &
           WGW%LGAMMA, WGW%LGAMMA, IMAG_GRIDS%B%NPOINTS_IN_GROUP)
      ! set inverse temperature 
      CHI%BETA = IMAG_GRIDS%BETA        
      CALL CLEAR_RESPONSE(CHI)
      NULLIFY(CHI)
   ENDDO

   ALLOCATE(CHI_TAU(1))
   CALL ALLOCATE_RESPONSEFUN_DISTRI( CHI_TAU(1), GDES_TAU%RES_NRPLWV_ROW, GDES_TAU%RES_NRPLWV_COL, &
        WGW%LGAMMA, WGW%LGAMMA, 1)
   ! set inverse temperature 
   CHI_TAU(1)%BETA = IMAG_GRIDS%BETA        

   ! OEP requires to invert chi(w=0) at Gamma point, 
   ! this object is stored in CHI0
   IF( LOEP ) THEN
      ALLOCATE(CHI0(1))
      CALL ALLOCATE_RESPONSEFUN_DISTRI( CHI0(1), GDES%RES_NRPLWV_ROW, GDES%RES_NRPLWV_COL, &
           WGW%LGAMMA, WGW%LGAMMA, 1)
      CHI0(1)%BETA = IMAG_GRIDS%BETA        
      CALL CLEAR_RESPONSE(CHI0(1))
   ENDIF

   CALL START_TIMING("GWLOOP")  ! super timer for e.g. one loop
   CALL START_TIMING("GREENS")  ! cummulative timer for individual steps

   IF (IO%IU6>=0) THEN
      WRITE(IO%IU6,'(A)') ' all response functions allocated'
   ENDIF
   CALL DUMP_ALLOCATE_TAG(IO%IU6, 'Response functions allocated')

   IU6_MEM=IO%IU6

   !store Fermi energy in wave function
   WHF%EFERMI=W%EFERMI
   IF (LALLOCATED_WCORR) WCORR%EFERMI=W%EFERMI

   CALL QP_PADE_READER(W,IO%IU5,IO%IU0,IO%IU6,NBANDSGW,NKPTS_IRZ)

!====================================================================
! Determine CHI_QPT on imaginary frequency axis 
! and add the q->0 limit 
!====================================================================
   CALL DETERMINE_CHI_Q_OMEGA( CHI_QPT, CHI_TAU, &
      GO_MAT, GU_MAT, GDES_MAT, &
      SUPER, S2E, W, WDES, WHF, WCORR, WCORRHF, WGW, WGWQ, &
      GDES_TAU, GDES, KPOINTS, NKPTS_IRZ, LATT_CUR, &
      NELM, IMAG_GRIDS, IO, IU6_MEM, &
      .TRUE., .TRUE., .TRUE., CHI0)

#ifdef debug
     CALL DUMPX_FREQ(1000, CHI0(1), GDES_TAU, IMAG_GRIDS%T%COMM_BETWEEN_GROUPS, &
        IMAG_GRIDS%T%COMM_IN_GROUP, WDES%COMM_KIN%NODE_ME, &
        IMAG_GRIDS%B%POINTS_LOCAL)
#endif

   ! get rid of correlated wavefunction
   IF (LALLOCATED_WCORR) CALL DEALLOCW(WCORR)

!=======================================================================
!   RPA calculations
!=======================================================================
    ! allocate correlation energy arrays
    CALL XI_ACFDT_SETUP_GG( COR, ENCUTGW, ENCUTGWSOFT)

    ! copy descriptor to current k-point
    DO NQ_COUNTER=1,S2E%NUMBER_OF_NQ  !this is over q-point (symmetry taken into account)
       NQ=S2E%NQ(NQ_COUNTER)
       CHI=>CHI_QPT(NQ_COUNTER)
       IF (NQ>NKPTS_IRZ) THEN
          ! this condition is never met, however, might be helpfull for 
          ! testing the code if the loop above is extended to all k-points
          ! equivalent k-point in IRZ
          NQ_IRZ=KPOINTS_FULL_ORIG%NEQUIV(KPOINT_IN_FULL_GRID(KPOINTS_FULL%VKPT(:,NQ),KPOINTS_FULL_ORIG))

          ! depending on whether we work in time or frequency, GDES_TAU or GDES
          DO I=1,IMAG_GRIDS%B%NPOINTS_IN_GROUP
#ifdef gammareal
             CHI%RESPONSER(:,:,I)=CHI_QPT(NQ_IRZ)%RESPONSER(:,:,I)
             CALL ROTATE_RES(CHI%RESPONSER(:,:,I), GDES, WGW%NGVECTOR(NQ), & 
               RES_KPTS_TRANS%CPHASE(:,NQ), RES_KPTS_TRANS%NINDPW(:,NQ),  &
               RES_KPTS_TRANS%LINV(NQ), RES_KPTS_TRANS%LSHIFT(NQ))
#else
             CHI%RESPONSEFUN(:,:,I)=CHI_QPT(NQ_IRZ)%RESPONSEFUN(:,:,I)
             CALL ROTATE_RES(CHI%RESPONSEFUN(:,:,I), GDES, WGW%NGVECTOR(NQ), &
               RES_KPTS_TRANS%CPHASE(:,NQ), RES_KPTS_TRANS%NINDPW(:,NQ),  &
               RES_KPTS_TRANS%LINV(NQ), RES_KPTS_TRANS%LSHIFT(NQ))
#endif
          ENDDO
       ENDIF
       CALL SETWDES(WGW, WGWQ, NQ )
#ifdef debug
! uncomment this if you want to compare to chi.F
     CALL DUMPX_FREQ(100, CHI, GDES_TAU, IMAG_GRIDS%T%COMM_BETWEEN_GROUPS, &
        IMAG_GRIDS%T%COMM_IN_GROUP, WDES%COMM_KIN%NODE_ME, &
        IMAG_GRIDS%B%POINTS_LOCAL)
#endif 
      ! evaluate the RPA correlation energy integral for a set of different
      ! energy cutoffs
      CALL CALCULATE_RPA_CORRELATION_ENERGY( WDES, CHI, GDES%RES_NRPLWV_ROW_DATA_POINTS, & 
      GDES%RES_NRPLWV_COL_DATA_POINTS, WGWQ, IMAG_GRIDS, LATT_CUR, COR, NQ, IDIR_MAX, IO)

       CALL STOP_TIMING("GREENS",IO%IU6,"ACFDT")

       IF (IO%IU6>=0) WRITE(IO%IU6,'(" q-point correlation energy ",2F14.6)') COR%CORRELATION_K(1)
       IF (IO%IU6>=0) WRITE(IO%IU6,'(" Hartree contr. to MP2      ",2F14.6)') COR%CORRMP2DIR_K(1)
    ENDDO

    ! linear regression 
    CALL LIN_REG_GG(COR, E, IO%IU6)

    ! deallocate correlation energy array
    CALL XI_ACFDT_DEALLOCATE_GG( COR )

!=======================================================================
!   if RPA calculations, finish now
!=======================================================================
  rpa_or_gw: IF ( LACFDT ) THEN

    CALL ON_SYMMETRY(LGAMMA, SYMM, W, WDES, LATT_CUR, LATT_INI, &
     T_INFO, DYN, INFO, KPOINTS, GRID, NONL_S, P, WHF, IO )

    CALL POTENTIAL_AND_CHARGE(1, &
       W, WDES, E, INFO, &
       GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
       LATT_CUR, T_INFO, P, SYMM, MIX, &
       CHTOT, CHTOTL, DENCOR, CVTOT, CSTRF, IRDMAX, &
       LMDIM, CQIJ, CDIJ, CRHODE, N_MIX_PAW, RHOLM, RHOLM_LAST, CHDEN, SV, IO)
    IF (LSINGLES) THEN
       CALL EDDIAG_RSE(HAMILTONIAN, &
            GRID, LATT_CUR, NONLR_S, NONL_S, W, WDES, SYMM, &
            LMDIM, CDIJ, CQIJ, SV, T_INFO, P, NKPTS_IRZ, E , CSINGLES, MP2SUM)
    ELSE
       CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
            LMDIM,CDIJ,CQIJ, 0,SV,T_INFO,P,IO%IU0,E%EXHF,&
            NBANDS_MAX=LAST_FILLED_OPTICS_THRESH( W, 1.E-8_q )/W%WDES%NB_PAR,EXHF_ACFDT=E%EXHF_ACFDT)
       E%EBANDSTR=BANDSTRUCTURE_ENERGY(WDES, W)
    ENDIF
    CALL STOP_TIMING("GREENS",IO%IU6,"EDDIAG")

    CALL DUMP_DIAGRAMMATIC_ENERGIES_RPA( "Energies using frozen KS orbitals", E, TOTEN, INFO, IO) 

    IF (LSINGLES) THEN
       IF (IO%IU6>=0) WRITE(IO%IU6,'(/," HF single shot energy change ",F18.8)') REAL(MP2SUM  ,q)
       IF (IO%IU6>=0) WRITE(IO%IU6,'(  " renormalized HF singles      ",F18.8)') REAL(CSINGLES,q)

       IF (IO%IU0>=0) WRITE(IO%IU0,'(" HF single shot energy change ",F18.8)') REAL(MP2SUM  ,q)
       IF (IO%IU0>=0) WRITE(IO%IU0,'(" renormalized HF singles      ",F18.8)') REAL(CSINGLES,q)
    ENDIF
  ELSE IF ( LCRPA ) THEN rpa_or_gw
!=======================================================================
! CRPA or NLRPA calculation
!=======================================================================
    !
    ! screened Coulomb potential W 
    !
    CALL SCREENED_POTENTIAL(S2E, CHI_QPT, WGW, WGWQ, &
     WDES, GDES, LATT_CUR, FSG0, IMAG_GRIDS, .TRUE., IO ) 
    IF(IO%IU0>=0)WRITE(IO%IU0,'(" Calculation of U in reciprocal space done")')
    IF(IO%IU0>=0)WRITE(17,'(" Calculation of U in reciprocal space done")')
#ifdef VASP2WANNIER90
    ! calculate coulomb matrix elements in localized basis 
    CALL COULOMB_MATRIX_ELEMENTS( &
    &   WDES,W,KPOINTS,GRID,T_INFO,INFO,P, &
    &   NONL_S,NONLR_S,SYMM,LATT_CUR,CQIJ,LMDIM, &
    &   CHI_QPT, IMAG_GRIDS, S2E, FSG0, IO)
#endif
  ELSE rpa_or_gw
!=======================================================================
!
! GW type calculations
!
!=======================================================================
    WMEAN=W ; NULLIFY(WMEAN%CELTOT, WMEAN%FERTOT,  WMEAN%AUXTOT)
    CALL ALLOCW_NOPLANEWAVE(WDES,WMEAN) ; WMEAN%FERTOT=W%FERTOT ; WMEAN%FERTOT = W%FERTOT ; WMEAN%AUXTOT=W%AUXTOT
    WHF_IN_KS=W ; NULLIFY(WHF_IN_KS%CELTOT, WHF_IN_KS%FERTOT,  WHF_IN_KS%AUXTOT)
    CALL ALLOCW_NOPLANEWAVE(WDES,WHF_IN_KS) ; WHF_IN_KS%FERTOT=W%FERTOT ; WHF_IN_KS%FERTOT = W%FERTOT ; WHF_IN_KS%AUXTOT=W%AUXTOT

    CALL GREENS_MAT_DES_INIT( GDES_MAT, IMAG_GRIDS%T, WDES%NB_TOT)

    NULLIFY(SIGMAW_MAT, GU_MAT, GO_MAT)
    IF (LGW .AND. ((.NOT. LG0W0 ) .AND.  (.NOT. LscQPGW))) THEN
       ! allocate space for self-energy in frequency domain, scaLAPACK distributed
       CALL ALLOCATE_GREENS_MAT(GDES_MAT, SIGMAW_MAT, NKPTS_IRZ, WDES%ISPIN,&
       IMAG_GRIDS%T%NPOINTS_IN_ROOT_GROUP)
       SIGMAW_MAT=0._q
       ! allocate space for Green function
       ! stores G_a,b in + and - imag time, required in the tau_gw loop if NELM>1
       IF (NELMGW>1) THEN
          ! avoid excessive storage, if not required
          CALL ALLOCATE_GREENS_MAT_GDEF(GDES_MAT, GU_MAT, NKPTS_IRZ, WDES%ISPIN,&
          IMAG_GRIDS%T%NPOINTS_IN_ROOT_GROUP)
          CALL ALLOCATE_GREENS_MAT_GDEF(GDES_MAT, GO_MAT, NKPTS_IRZ, WDES%ISPIN,&
          IMAG_GRIDS%T%NPOINTS_IN_ROOT_GROUP)
       ENDIF
    ENDIF
    !
    ! screened Coulomb potential W 
    !
    CALL SCREENED_POTENTIAL(S2E, CHI_QPT, WGW, WGWQ, &
     WDES, GDES, LATT_CUR, FSG0, IMAG_GRIDS, .TRUE., IO ) 
    IF(IO%IU0>=0)WRITE(IO%IU0,'(" Calculation of W done, Calculation of self-energy started for")')
    IF(IO%IU0>=0)WRITE(17,'(" Calculation of W done, Calculation of self-energy started for")')

! TODO: dump warning for forcing NELM = 1 if  ALGO = G0W0R ; NELM>1 is set
    IF ( NELMGW > 1 .AND. LG0W0 ) NELMGW = 1 
    FLWTOT_OLD=0.0_q
    EFERMI_SAVE=0.0_q
scf: DO NELM=1,NELMGW
    IF (NELM>1) THEN
       CALL NO_SYMMETRY( .FALSE., LGAMMA, SYMM, W, WDES, LATT_CUR, LATT_INI, &
          T_INFO, INFO, KPOINTS, GRID, NONL_S, P, WHF, IO )
    ENDIF

    CALLMPI( M_barrier( WDES%COMM_KIN ) ) 

    !-------------------------------------------------------------------
    ! recalculate screened Coulomb potential W
    !-------------------------------------------------------------------
    CALL RECALCULATE_W_AND_PHI_FUNCTIONAL( CHI_QPT, CHI_TAU, &
       GO_MAT, GU_MAT, GDES_MAT, &
       SUPER, S2E, W, WDES, WHF, WCORR, WCORRHF, WGW, WGWQ, &
       GDES_TAU, GDES, KPOINTS, NKPTS_IRZ, LATT_CUR, &
       RES_KPTS_TRANS, E, COR, FSG0, &
       NELM, IMAG_GRIDS, IO, IU6_MEM )

    ! allocate storage for self-energy in frequency space evaluated  w.r.t. Bloch states
    ! diagonal elements only
    ! these are deallocated after calculation of rotation matrix
    ALLOCATE(SIGMA_COS(NOMEGA, W%WDES%NB_TOT,NKPTS_IRZ,W%WDES%ISPIN))
    ALLOCATE(SIGMA_SIN(NOMEGA, W%WDES%NB_TOT,NKPTS_IRZ,W%WDES%ISPIN))
    SIGMA_COS=0
    SIGMA_SIN=0

    !allocate storage for self-energies in Bloch domain on time axis
    ALLOCATE(SIGMAO_TAU(IMAG_GRIDS%NOMEGA,W%WDES%NB_TOT,NKPTS_IRZ,W%WDES%ISPIN))
    ALLOCATE(SIGMAU_TAU(IMAG_GRIDS%NOMEGA,W%WDES%NB_TOT,NKPTS_IRZ,W%WDES%ISPIN))
    SIGMAO_TAU=zero
    SIGMAU_TAU=zero

    ! density matrix change from correlation contribution
    CALL ALLOCATE_GREENS_MAT_GDEF(GDES_MAT, CORR_MAT, NKPTS_IRZ, WDES%ISPIN, 1)
    CORR_MAT=0
    IF (CALCULATE_RPA_FORCES().AND. WDES%LOVERL) THEN
       CALL ALLOCATE_GREENS_MAT_GDEF(GDES_MAT, SIGMA_MAT, NKPTS_IRZ, WDES%ISPIN, 1)
       SIGMA_MAT=0
    ENDIF
    FRNL=0
    MP2SUM=0 ; MP2SUM2=0
    E_MP2=0
    IU6_MEM=IO%IU6
  !----------------------------------------------------------------------- 
  ! start loop over time poins for calculation of SIGMA = G W    
  tau_points_gw: DO NTAU_ROOT=1, IMAG_GRIDS%T%NPOINTS_IN_ROOT_GROUP
  !----------------------------------------------------------------------- 

    CALL SETUP_TAU_INDICES(NTAU_ROOT,IMAG_GRIDS%T)

    DO ISP=1,W%WDES%ISPIN
    !------------------------------------------------------------------- 
    ! calculate the occupied and unoccupied Green function
    !------------------------------------------------------------------- 

       ALLOCATE(SIGMAU(NKPTS_IRZ), SIGMAO(NKPTS_IRZ ))
       CALL ALLOCATE_G_RECIPROCAL( GDES_TAU, SIGMAO, NKPTS_IRZ)
       CALL ALLOCATE_G_RECIPROCAL( GDES_TAU, SIGMAU, NKPTS_IRZ)

       !allocate storage for self-energies in Bloch domain
       ALLOCATE(SIGMAO_DIAG(W%WDES%NB_TOT))
       ALLOCATE(SIGMAU_DIAG(W%WDES%NB_TOT))
       SIGMAO_DIAG=0
       SIGMAU_DIAG=0

       NULLIFY(SIGMAU_MAT)
       NULLIFY(SIGMAO_MAT)
     ! if commented out the full matrix is not calculated
       ! SIGMAU_MAT and SIGMAO_MAT are the self-energies in Bloch basis 
       ! for (un)occupied states in time domain
       ! deallocated at the end of tau_points_gw loop 
       ! since at that point the self-energy is transformed to frequency domain,
       ! and stored in SIGMAW_MAT
       CALL ALLOCATE_GREENS_MAT_GDEF(GDES_MAT, SIGMAU_MAT, 1, 1, 1)
       CALL ALLOCATE_GREENS_MAT_GDEF(GDES_MAT, SIGMAO_MAT, 1, 1, 1)

       !compute occupied and unoccupied Green function in tau domain
       CALL CALCULATE_G_POSSIBLY_SC( W, WDES, WHF, WCORR, WCORRHF,&
          GDES_TAU, GDES_MAT, NKPTS_IRZ, ISP, IMAG_GRIDS%T, NELM==1 .OR.  LscQPGW, NTAU_ROOT, &
          GU, GO, GU_MAT, GO_MAT, IU6_MEM, IO, .FALSE. )

       CALL STOP_TIMING("GREENS",IO%IU6,"GREENS")
    !------------------------------------------------------------------- 
    ! compute Sigma_k=G W in time domain 
    !-------------------------------------------------------------------
       DO NQ_COUNTER=1,S2E%NUMBER_OF_NQ
          IF (NQ_COUNTER/=S2E%NQ(NQ_COUNTER)) THEN
             CALL vtutor%error("error in VASP: this routine does not support NKRED or shifted grids yet")
          ENDIF
       ENDDO

       CALL START_TIMING("G")
       IU6_TMP=IU6_MEM
       DO NQ=1,WDES%NKPTS  ! loop over all difference vectors in the full BZ
          IF (IO%IU0>=0) WRITE(IO%IU0,'("NQ=",I4,3F10.4," TAU=",1F10.4", ")')&
          NQ,CHI_TAU(1)%VKPT,IMAG_GRIDS%T%POINT_CURRENT
          IF (IO%IU0>=0 .AND. NELM==1 ) WRITE(17,'("NQ=",I4,3F10.4," TAU=",1F10.4", ")')&
          NQ,CHI_TAU(1)%VKPT,IMAG_GRIDS%T%POINT_CURRENT

       ! k-point outside IRZ rotate potential from IRZ to desired k-point
          IF (NQ>NKPTS_IRZ) THEN
#ifndef gammareal
             ! equivalent k-point in IRZ (somewhat akward, but correct)
             NQ_IRZ=KPOINTS_FULL_ORIG%NEQUIV(KPOINT_IN_FULL_GRID(KPOINTS_FULL%VKPT(:,NQ),KPOINTS_FULL_ORIG))

             !CHI should point to current q-point
             CHI=>CHI_QPT(NQ_IRZ)

             IF ( IMAG_GRIDS%T%LDO_POINT_LOCAL ) CALL CLEAR_RESPONSE(CHI_TAU(1)) 
             ! we take CHI_QPT from IrBZ and transform it to CHI_TAU in BZ
             ! set VKPT to k-point in IRZ
             CALL SET_RESPONSE_KPOINT(CHI_TAU(1), WDES%VKPT(:,NQ_IRZ), NQ_IRZ)

             ! inverse cosine transformation of screened Coulomb potential W for current tau-point
             ! W is stored in CHI_TAU
             ! transform W from nu to tau; backward transformation set by LFORWARD=.FALSE. 
             CALL TRANS_TIME_FREQUENCY(GDES, GDES_TAU, CHI, CHI_TAU(1), &
             IMAG_GRIDS, NTAU_ROOT, .FALSE. , IO%IU0)
             ! note: depending on whether we work in time or frequency, GDES_TAU or GDES
             ! now transform to NQ
             CALL ROTATE_RES(CHI_TAU(1)%RESPONSEFUN(:,:,1), GDES_TAU, WGW%NGVECTOR(NQ), &
                  RES_KPTS_TRANS%CPHASE(:,NQ), RES_KPTS_TRANS%NINDPW(:,NQ),  &
                  RES_KPTS_TRANS%LINV(NQ), RES_KPTS_TRANS%LSHIFT(NQ))

             ! now set k-points to NQ
             CALL SET_RESPONSE_KPOINT(CHI_TAU(1), WDES%VKPT(:,NQ), NQ)
             CALL STOP_TIMING("G",IU6_TMP,"*transw")
#endif
          ELSE
             !CHI should point to current q-point 
             CHI=>CHI_QPT(NQ)

             IF ( IMAG_GRIDS%T%LDO_POINT_LOCAL ) CALL CLEAR_RESPONSE(CHI_TAU(1)) !one q-point at a time 
             CALL SET_RESPONSE_KPOINT(CHI_TAU(1), WDES%VKPT(:,NQ), NQ)

             ! inverse cosine transformation of screened Coulomb potential W for current tau-point
             ! W is stored in CHI_TAU, for a given q-point
             ! transform W from nu to tau for a given tau point; backward transformation set by LFORWARD=.FALSE. 
             CALL TRANS_TIME_FREQUENCY(GDES, GDES_TAU, CHI, CHI_TAU(1), &
             IMAG_GRIDS, NTAU_ROOT, .FALSE. , IO%IU0 )
             CALL STOP_TIMING("G",IU6_TMP,"*transw")
          ENDIF

          DO NK1=1,NKPTS_IRZ  ! loop restricted to IRZ
             IF (MOD(NK1-1,WDES%COMM_KINTER%NCPU)/=WDES%COMM_KINTER%NODE_ME-1) CYCLE
             I = DETERMINE_NTAU_GLOBAL(WDES%COMM_KIN%NODE_ME, IMAG_GRIDS%T%NPOINTSC,& 
             IMAG_GRIDS)
             TAU_TMP = IMAG_GRIDS%TAU_WEIGHT( I ) 

             ! Sigma is at NK1, W is at NQ, and GU at NK2=NK1-NQ and GO at NQ
             ! k2=q-k1
             NK2=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,NK1)-CHI_TAU(1)%VKPT(:),KPOINTS_FULL)

             ! compute self-energy for occupied orbital  SIGMAO = G( -tau) W(tau) = GO W
             CALL CALCULATE_SIGMA_TAU( LMDIM, WHF, WGW, GO(NK2), GU(NK1), SIGMAO(NK1), GDES_TAU, CHI_TAU(1),&
                  NQ, NK2, NK1, LATT_CUR, IU6_TMP, TAU_TMP, FRNL, CALCULATE_RPA_FORCES())

             CALL STOP_TIMING("G",IU6_TMP,"*sigmaocc")

             ! compute self-energy for unoccupied orbitals  SIGMAO = G( tau) W(tau) = GU W
             CALL CALCULATE_SIGMA_TAU( LMDIM, WHF, WGW, GU(NK2), GO(NK1), SIGMAU(NK1), GDES_TAU, CHI_TAU(1),&
                  NQ, NK2, NK1, LATT_CUR, IU6_TMP, TAU_TMP, FRNL, CALCULATE_RPA_FORCES())
             CALL STOP_TIMING("G",IU6_TMP,"*sigmaun")

             IU6_TMP=-1
          ENDDO ! k-points
       ENDDO ! q-point
       !now we have CHI_TAU(1) for all k-vectors and in tau domain
       CALL STOP_TIMING("GREENS",IO%IU6,"SIGMA")


       CALL START_TIMING("G")

       IU6_TMP=IU6_MEM
       ! now determine the matrix elements <i|Sigma(tau)|a>
       DO NK1=1,NKPTS_IRZ
          IF (MOD(NK1-1,WDES%COMM_KINTER%NCPU)/=WDES%COMM_KINTER%NODE_ME-1) CYCLE

          I = DETERMINE_NTAU_GLOBAL(WDES%COMM_KIN%NODE_ME, IMAG_GRIDS%T%NPOINTSC,&
          IMAG_GRIDS)
          TAU_TMP = KPOINTS_ORIG%WTKPT(NK1)*IMAG_GRIDS%TAU_WEIGHT( I ) 
          IF (L2ORDER) THEN
             ! for second order (L2ORDER): G(tau) G( -tau) W(tau) is the MP2 energy
             CALL CONTRACT_SIGMA_G_TAU( WHF, GU(NK1), SIGMAO(NK1), GDES_TAU, NK1, &
                  IO%IU6, TAU_TMP, IMAG_GRIDS%T%LDO_POINT_LOCAL,E_MP2 )
          ENDIF
          CALL CONTRACT_SIGMA_G_TAU_DER( W%WDES, NK1, TAU_TMP, &
               GU(NK1), SIGMAO(NK1), GDES_TAU, NONL_S, NONLR_S, LATT_CUR, FRNL,&
               IMAG_GRIDS%T%LDO_POINT_LOCAL)

          ! determine <i| G( -tau) W(tau) | a> in the orbital basis
          ! W(tau) symmetric (?)
          ! input : SIGMAO(NK1)%(g,g')
          ! output: SIGMAO_MAT(ab')
          CALL CALCULATE_SIGMA_ORBITAL(W, SIGMAO(NK1), SIGMAO_DIAG, SIGMAO_MAT, &
               GDES_MAT, GDES_TAU, IMAG_GRIDS, -1, NK1, ISP, .TRUE., MP2SUM, TAU_TMP )
          CALL STOP_TIMING("G",IU6_TMP,"*s_orbital")

          IF (L2ORDER) THEN
             ! for second order (L2ORDER): G(tau) G( -tau) W(tau) is the MP2 energy
             CALL CONTRACT_SIGMA_G_TAU( WHF, GO(NK1), SIGMAU(NK1), GDES_TAU,NK1, & 
                  IO%IU6, TAU_TMP, IMAG_GRIDS%T%LDO_POINT_LOCAL,E_MP2 )
          ENDIF
          CALL CONTRACT_SIGMA_G_TAU_DER( W%WDES, NK1, TAU_TMP, &
               GO(NK1), SIGMAU(NK1), GDES_TAU, NONL_S, NONLR_S, LATT_CUR, FRNL,&
               IMAG_GRIDS%T%LDO_POINT_LOCAL)

          ! determine <i| G( tau) W(tau) | a> in the orbital basis
          CALL CALCULATE_SIGMA_ORBITAL(W, SIGMAU(NK1), SIGMAU_DIAG, SIGMAU_MAT, &
               GDES_MAT, GDES_TAU, IMAG_GRIDS, -1, NK1, ISP, .FALSE., MP2SUM2, TAU_TMP )
          CALL STOP_TIMING("G",IU6_TMP,"*s_orbital")

          ! Fourier transform Sigma to imaginary frequency (on those cores we have information)
          IF (IMAG_GRIDS%T%LDO_POINT_LOCAL) THEN
             ! Fourier transform the diagonal part of Sigma to imaginary frequency 
             ! (on those cores we have information)
             ! SIGMA0_DIAG and SIGMAU_DIAG are in tau, SIGMA_COS and SIGMA_SIN are in nu 
             CALL FT_G_OR_SIGMA(W, SIGMAO_DIAG, SIGMAU_DIAG, NTAU_ROOT, &
             IMAG_GRIDS, SIGMA_COS, SIGMA_SIN, NK1, ISP )
             IF (ASSOCIATED(SIGMAO_MAT).AND. LscQPGW) THEN
                ! calculate G SIGMA and add it to CORR_MAT storing the effective (non-Hermitian) potential
                !  potential =\int tau G0(tau) Sigma(tau) d tau
                CALL ADD_G_SIGMA( W, SIGMAO_MAT(:, 1, 1, 1),&
                SIGMAU_MAT(:, 1, 1, 1), CORR_MAT(:, NK1, ISP, 1), & 
                NTAU_ROOT, IMAG_GRIDS, NK1, ISP, GDES_MAT%DESC )
             ELSE IF (ASSOCIATED(SIGMAO_MAT).AND. LG0W0 ) THEN
               ! obtain the first order change of density matrix induced by the self-energy Sigma
               ! G1 = G0 + G0 Sigma G0, 
               ! and \Delta \gamma=\int w G0(w) Sigma(w) G0(w) is the change of density matrix
               ! G0 are taken as mean field G with eigenvectors and eigenvalues from WAVECAR
               ! \Delta \gamma is stored in CORR_MAT, the complete change is accumulated during the tau loop 
               ! This should not be done if LGW0=.TRUE. since the correlation contribution
               ! to gamma is calculated below in GAMMA_FROM_G_W routine
               ! CORR_MAT is properly set to 0 there so this call has no effect anyway
               CALL ADD_G_SIGMA_G_new( W, SIGMAO_MAT(:, 1, 1, 1), SIGMAU_MAT(:, 1, 1, 1), CORR_MAT(:, NK1, ISP, 1), & 
                  NTAU_ROOT, IMAG_GRIDS, NK1, ISP, GDES_MAT%DESC )
               IF (ASSOCIATED(SIGMA_MAT)) &
               CALL ADD_G_SIGMA_Gxw_new( W, SIGMAO_MAT(:, 1, 1, 1), SIGMAU_MAT(:, 1, 1, 1), SIGMA_MAT(:, NK1, ISP, 1), & 
                  NTAU_ROOT, IMAG_GRIDS, NK1, ISP, GDES_MAT%DESC )
            ENDIF
            ! store diagonal self-energy elements on time axis 
            SIGMAO_TAU(DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME,&
            NTAU_ROOT, IMAG_GRIDS ),:,NK1,ISP) = SIGMAO_DIAG(:)
            SIGMAU_TAU(DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME,&
            NTAU_ROOT, IMAG_GRIDS ),:,NK1,ISP) = SIGMAU_DIAG(:)
          ENDIF
          IF (LGW .AND. ((.NOT. LG0W0 ) .AND.  (.NOT. LscQPGW))) THEN
             ! transform the selfenergy at the current tau point (SIGMAO_MAT and SIGMAU_MAT)
             ! to imaginary frequency and store it in SIGMAW_MAT
             ! the complete self-energy is accumulated during the tau loop
             CALL FT_SIGMA_MAT( W, SIGMAO_MAT(:, 1, 1, 1), SIGMAU_MAT(:, 1, 1, 1), & 
                  SIGMAW_MAT(:, NK1, ISP,:), NTAU_ROOT, IMAG_GRIDS, GDES_MAT )
          ENDIF

          IU6_TMP=-1
       ENDDO ! k-point
       CALL STOP_TIMING("GREENS",IO%IU6,"<i|S|n>")

       ! deallocate Green function arrays (frees some temporary storage)
       CALL DEALLOCATE_G( SIGMAO)
       CALL DEALLOCATE_G( SIGMAU)
       CALL DEALLOCATE_G( GU)
       CALL DEALLOCATE_G( GO)

       DEALLOCATE( SIGMAO_DIAG , SIGMAU_DIAG)
       DEALLOCATE( SIGMAO , SIGMAU)  !!!added
       DEALLOCATE(GU, GO)
       CALL DEALLOCATE_GREENS_MAT_GDEF(GDES_MAT, SIGMAU_MAT)
       CALL DEALLOCATE_GREENS_MAT_GDEF(GDES_MAT, SIGMAO_MAT)
  !------------------------------------------------------------------- 
       IU6_MEM=-1  ! no further detailed timing of memory dumps
    ENDDO ! spin
  ENDDO tau_points_gw
  !-------------------------------------------------------------------
  ! at this point we have
  ! CORR_MAT has the first order change of density matrix induced by SIGMA
  !   (to get improved mean field or natural orbitals)
  !   CORR_MAT is overwritten few lines below if LGW0=.TRUE.
  ! SIGMAW_MAT stores self-energy in imaginary frequency nu (positive only)
  !   (to obtain interacting G via Dyson equation), only stored if LGW0=.TRUE.
  ! SIGMA_COS and SIGMA_SIN are diagonal elements of self-energy in nu
  !   (hopefully to get QP corrections via analytic continuation)

    CALLMPI( M_sum_g(GDES_MAT%COMM_INTER, CORR_MAT, SIZE(CORR_MAT)))

  ! diagonal part of self-energy on imaginary time axis 
    CALLMPI( M_sum_g(IMAG_GRIDS%T%COMM_BETWEEN_GROUPS, SIGMAO_TAU, SIZE(SIGMAO_TAU)))
    CALLMPI( M_sum_g(IMAG_GRIDS%T%COMM_BETWEEN_GROUPS, SIGMAU_TAU, SIZE(SIGMAU_TAU)))
  ! compute G*Sigma in frequency domain 
    CALL GALITSKII_MIGDAL_DIAG( W, SIGMAO_TAU, SIGMAU_TAU, IMAG_GRIDS, IO)

    IF (ASSOCIATED(SIGMA_MAT)) & 
        CALLMPI( M_sum_g(GDES_MAT%COMM_INTER, SIGMA_MAT, SIZE(SIGMA_MAT)))
    CALLMPI( M_sum_z(GDES_MAT%COMM_INTER, MP2SUM, 1))
    CALLMPI( M_sum_z(GDES_MAT%COMM_INTER, MP2SUM2, 1))
    CALLMPI( M_sum_z(GDES_MAT%COMM_INTER, SIGMA_SIN, SIZE(SIGMA_SIN)))
    CALLMPI( M_sum_z(GDES_MAT%COMM_INTER, SIGMA_COS, SIZE(SIGMA_COS)))
    CALLMPI( M_sum_d(GDES_MAT%COMM_INTER, FRNL, SIZE(FRNL)))

    IF (IO%IU6>=0) THEN
       WRITE(IO%IU6,'(" correlation energy Sigma(-tau) G(tau) ",2F20.10)') MP2SUM
       WRITE(IO%IU6,'(" correlation energy Sigma(tau) G(-tau) ",2F20.10)') MP2SUM2
    ENDIF
    IF (IO%IU0>=0) THEN
       WRITE(IO%IU0,'(" correlation energy Sigma(-tau) G(tau) ",2F20.10)') MP2SUM
       WRITE(IO%IU0,'(" correlation energy Sigma(tau) G(-tau) ",2F20.10)') MP2SUM2
    ENDIF

  !-------------------------------------------------------------------
  ! now determine Hartree-Fock Hamiltonian
  !-------------------------------------------------------------------
    CALL ON_SYMMETRY(LGAMMA, SYMM, W, WDES, LATT_CUR, LATT_INI, &
     T_INFO, DYN, INFO, KPOINTS, GRID, NONL_S, P, WHF, IO )

    CALL ALLOCATE_GREENS_MAT_GDEF(GDES_MAT, CHAM_MAT, NKPTS_IRZ, WDES%ISPIN, 1)

  iterhf: DO NELM_HF=1,NELMHF
    ! Hartree term from Green function:
    ! determine Hatree and kinetic term in natural orbital basis 
    ! This is U_0 in PRA73,012511 Eqs. (3) and (B5)
    IF (NELMHF==1) THEN
       CALL POTENTIAL_AND_CHARGE(NELM, &
          W, WDES, E, INFO, &
          GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
          LATT_CUR, T_INFO, P, SYMM, MIX, &
          CHTOT, CHTOTL, DENCOR, CVTOT, CSTRF, IRDMAX, &
          LMDIM, CQIJ, CDIJ, CRHODE, N_MIX_PAW, RHOLM, RHOLM_LAST, CHDEN, SV, &
          IO, OEP)
    ELSE
       CALL POTENTIAL_AND_CHARGE(NELM_HF, &
          W, WDES, E, INFO, &
          GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
          LATT_CUR, T_INFO, P, SYMM, MIX, &
          CHTOT, CHTOTL, DENCOR, CVTOT, CSTRF, IRDMAX, &
          LMDIM, CQIJ, CDIJ, CRHODE, N_MIX_PAW, RHOLM, RHOLM_LAST, CHDEN, SV,&
          IO, OEP)
    ENDIF
    ! determine <i|H|a>= <i| T + V_Hartree + V_exchange | a> store in CHAM_MAT(i,a,nk,ispin)
    ! for LGW0 CHAM_MAT will be set to unitary matrix diagonalizing <i|H|a> U_aj = epsilon_j U_jj
    !  and WMEAN%CELTOT is set to corresponding eigenvalues
    ! for G0W0  WMEAN%CELTOT = <i | T + V_Hartree + V_exchange | i>
    E%ETRGHF = 0
    DTRGHF = 0
    CALL EDDIAG_SIMPLE(HAMILTONIAN, &
         GRID, LATT_CUR, NONLR_S, NONL_S, W, WMEAN, WHF_IN_KS, WDES, LGW .AND. ((.NOT. LG0W0 ) .AND.  (.NOT. LscQPGW)), SYMM, &
         LMDIM, CDIJ, CQIJ, SV, T_INFO, P, NKPTS_IRZ, CHAM_MAT, GDES_MAT, E, IO%IU6 )
    CALL STOP_TIMING("GREENS",IO%IU6,"EDDIAG")

    IF ((LG0W0.OR.LscQPGW) .OR. NELM==1) THEN
       IF (IO%IU6>=0) WRITE(IO%IU6,'(1X,A)', ADVANCE="NO") "QP shifts evaluated in KS or natural orbital/ Bruckner basis"
       CALL WRITE_SELF_ENERGY(WMEAN, SIGMA_COS, SIGMA_SIN, NKPTS_IRZ, IMAG_GRIDS, NBANDSGW, IO%IU6)
       !pL Analytic continuation of Sigma(iw) to real axis and calculate the QP poles
       CALL QP_PADE_FIT_FROM_TAU(W,WHF_IN_KS,IMAG_GRIDS,SIGMAU_TAU,SIGMAO_TAU,NKPTS_IRZ,NBANDSGW,IO)
    ENDIF

    TOTEN=E%EBANDSTR+E%DENC+E%XCENC+E%TEWEN+E%PSCENC+E%PAWPS+E%PAWAE+INFO%EALLAT+E%EXHF
    ! Galitskii-Migdal formula represent internal energy (see Fetter-Walecka Eq. 23.15)
    ! so exclude entropy term here 
    E%ETRGHF = TOTEN 
    ! compute change in bandstructure energy between iterations  
    DTRGHF = E%EBANDSTR
    
    CALL DUMP_DIAGRAMMATIC_ENERGIES_RPA("Energies using frozen KS orbitals", E, TOTEN, INFO, IO ) 

    ! following holds true if ALGO = GWRK | GW0RK
    !
    ! the electronic HF energy is calculated as 
    !
    ! E_HF = \sum_i f_i e_i - E_dc, where E_dc = 1/2( V_h - V_x )
    !
    ! -E_dc = E%DENC + E%EXHF + ( E%XCENC + E%PAWPS + E%PAWAE ) 
    ! 
    ! these terms essentially are the U_0+Tr( G_ref * Sigma_x ) term in PRA 73, 012511 Eq. (B5)
    IF (LGW .AND. ((.NOT. LG0W0 ) .AND.  (.NOT. LscQPGW))) THEN
       ! GW0 or GW: use HF Green's function as reference Green's function
       CALL DENSTA_SPIN( IO%IU0, IO%IU6, WDES, WMEAN, KPOINTS, INFO%NELECT, &
            INFO%NUP_DOWN,  E%EENTROPY, EFERMI_MEAN, KPOINTS%SIGMA, .FALSE., &
            NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
       CALL WRITE_EIGENVAL_NBANDS(  WMEAN%WDES, WMEAN, IO%IU6, LAST_FILLED_OPTICS(W)*2)

       ! reset Fermi-energy and one-electron occupancies to previous value
       ! we use the mean field part only to get a convergent Matsubara summation (integration)
       WMEAN%FERTOT(:,1:NKPTS_IRZ,:)=W%FERTOT(:,1:NKPTS_IRZ,:)  ! copy occupancies
       EFERMI_MEAN=W%EFERMI

       E%EBANDSTR=BANDSTRUCTURE_ENERGY(WDES, WMEAN)
       CALL DUMP_DIAGRAMMATIC_ENERGIES_RPA( "Energies after diagonalization of HF Hamiltonian (single shot)", E, TOTEN, INFO, IO) 
       IF (IO%IU6>=0) THEN
          WRITE(IO%IU6,"(/,' The Hartree-Fock Fermi energy was updated to be between HOMO (VB) and LUMO (CB) ',/ &
          & ' E-fermi : ',2F12.6,/)" ) EFERMI_MEAN
       ENDIF
       IF (IO%IU0>=0.AND. IO%NWRITE>2 ) THEN
          WRITE(IO%IU0,"(/,' The Hartree-Fock Fermi energy was updated to be between HOMO (VB) and LUMO (CB) ',/ &
          & ' E-fermi : ',2F12.6,/)" ) EFERMI_MEAN
       ENDIF
    ELSE
       ! G0W0: use DFT Green's function as reference Green's function
       EFERMI_MEAN=W%EFERMI     ! copy Fermi-energy
       WMEAN%CELTOT(:,1:NKPTS_IRZ,:)=W%CELTOT(:,1:NKPTS_IRZ,:)  ! copy eigenvalues
       WMEAN%FERTOT(:,1:NKPTS_IRZ,:)=W%FERTOT(:,1:NKPTS_IRZ,:)  ! copy occupancies
    ENDIF
    IF (LFERMIGW) THEN
       CALL ITERATE_EFERMI_FROM_SIGMA_DIAG(WMEAN, NKPTS_IRZ, INFO%NELECT, INFO%NUP_DOWN, EFERMI_MEAN, EFERMI_SLOPE, &
            IMAG_GRIDS, SIGMA_COS, SIGMA_SIN, GDES_MAT, IO )
       WMEAN%EFERMI=EFERMI_MEAN
    ENDIF
    ! in last iteration release diagonal entries of sigma 
    IF ( NELM_HF == NELMHF ) THEN        
       DEALLOCATE( SIGMA_COS , SIGMA_SIN)
       DEALLOCATE( SIGMAO_TAU , SIGMAU_TAU)
    ENDIF  
  !----------------------------------------------------------------------- 
  ! calculate correlated density matrix
  !----------------------------------------------------------------------- 
    IF (( LGW .AND. ((.NOT. LG0W0 ) .AND.  (.NOT. LscQPGW))) .AND. NELM_HF==1 ) THEN
       ! if diag_cham is set, the self-energy is always stored in the NO basis
       ! if diag_cham is not set, the self-energy is now transformed from the NO basis to the mean field basis
       !   the second option is much more efficient 
#ifndef diag_cham
       ! now rotate selfenergy into Hartree-Fock basis
       CALL ROTATE_GREEN_SIGMA(WDES, CHAM_MAT, SIGMAW_MAT, IMAG_GRIDS, GDES_MAT)
       ! Pade fit of self-energy in Hartree-Fock basis
       CALL PADE_FIT_GW0( W, WMEAN, IMAG_GRIDS, GDES_MAT, SIGMAW_MAT, NKPTS_IRZ, NBANDSGW, IO )
       ! and rotate KS orbitals to diagonalize Hartree-Fock Hamiltonian
       CALL ROTATE_KS_ORBITALS(W, WDES, CHAM_MAT, GDES_MAT)
       CALL STOP_TIMING("GREENS",IO%IU6,"ROTSIGMA")
#endif

       IF (LFERMIGW) THEN
          CALL ITERATE_EFERMI_G_FROM_SIGMA( WMEAN, NKPTS_IRZ, INFO%NELECT, &
          INFO%NUP_DOWN, EFERMI_MEAN, EFERMI_SLOPE, SIGMAW_MAT, IMAG_GRIDS,&
          GDES_MAT, CHAM_MAT, IO )
! mix the Fermi-level with previous one, and store it as the reference value for the next iteration
! mixing might accelerate convergence
!          W%EFERMI=EFERMI_MEAN*0.8+W%EFERMI*0.2
           W%EFERMI=EFERMI_MEAN
           WMEAN%EFERMI=EFERMI_MEAN
           EFERMI_VAR=SQRT(DOT_PRODUCT(EFERMI_SAVE-EFERMI_MEAN,EFERMI_SAVE-EFERMI_MEAN))
           IF ( IO%IU0>=0 .AND.  IO%NWRITE>2 ) THEN
              WRITE(IO%IU0,"(' Updated Fermi energy (new,old): ',3F14.8)" ) EFERMI_MEAN, &
              EFERMI_VAR
           ENDIF
           EFERMI_SAVE = EFERMI_MEAN
       ENDIF
       CALL STOP_TIMING("GREENS",IO%IU6,"EFERMI")
       IF (IO%IU6>=0) THEN
          WRITE(IO%IU6,"(/,' The GW-Fermi energy was updated to conserve the number of electrons',/ &
          & ' E-fermi : ',2F8.4,/)" ) EFERMI_MEAN
       ENDIF

       ! calculate interacting Greens function G = (iw - Sigma(w) - H_mean)^-1 + (iw - H_mean)^-1       
       ! H_mean is the Hartree-Fock Hamiltonian
       ! use in place version storing the Green's function in SIGMAW_MAT
       ! this routine also determines the Galitskii-Migdal Energy:
       ! E_GM = int dw Tr( G(iw) Sigma( iw ))

       ! initialize diagrammatic approximations to energy 
       E%ELOGG0 = 0   
       E%ELOG1G = 0
       E%ECGWGM = 0
       E%EKLLOG = 0
       E%EKLGG0 = 0
       CALL CALCULATE_G_FROM_SIGMA( WMEAN, NKPTS_IRZ, EFERMI_MEAN, SIGMAW_MAT, &
            IMAG_GRIDS, GDES_MAT, CHAM_MAT, E, IO )

       ! if for terms depending of O( 1/omega ) at T=0
       IF ( .NOT. LFINITE_TEMPERATURE ) THEN
          CALL SUBTRACT_GHF_FROM_G( WMEAN, NKPTS_IRZ, EFERMI_MEAN, SIGMAW_MAT, &
               IMAG_GRIDS, GDES_MAT, CHAM_MAT )
       ENDIF 

       CALL STOP_TIMING("GREENS",IO%IU6,"SIGMA2G")

       IF (ASSOCIATED(GO_MAT)) THEN
          ! TODO: move loop into INV_FT_SIGMA_MAT
          DO ISP=1,WDES%ISPIN
          DO NK1=1,NKPTS_IRZ
             DO NTAU_ROOT=1, IMAG_GRIDS%T%NPOINTS_IN_ROOT_GROUP

                CALL SETUP_TAU_INDICES( NTAU_ROOT, IMAG_GRIDS%T )
                ! transform the correlated part of Greens function to imaginary time
                ! occ and unocc components SIGMAW_MAT -> GU_MAT and GO_MAT
                CALL INV_FT_SIGMA_MAT(W, W%WDES%NB_TOT, GO_MAT(:, NK1, ISP, NTAU_ROOT), GU_MAT(:, NK1, ISP, NTAU_ROOT),&
                     SIGMAW_MAT(:, NK1, ISP,:), NTAU_ROOT, IMAG_GRIDS, GDES_MAT )
             ENDDO
          ENDDO
          ENDDO
       ENDIF

       ! write the diagonal part of SIGMAW_MAT
!      CALL WRITE_SELFENERGY_W(W, NKPTS_IRZ, SIGMAW_MAT, IMAG_GRIDS, GDES_MAT, IO%IU0)

       ! calculate contribution to density matrix from correlation contribution (overwrite CORR_MAT)
       CALL GAMMA_FROM_G_W( WMEAN, NKPTS_IRZ, SIGMAW_MAT, CORR_MAT, &
            IMAG_GRIDS, GDES_MAT, CHAM_MAT, LUNOCCUPIED=LHOLEGF  )

       ! at this point SIGMAW_MAT and CORR_MAT contain the interacting G-G_mean
       ! and the corresponding density matrix
       CALLMPI( M_sum_g(GDES_MAT%COMM_INTER, CORR_MAT, SIZE(CORR_MAT)))
       SIGMAW_MAT=0._q

       CALL STOP_TIMING("GREENS",IO%IU6,"GAMMA")
    ENDIF
    !
    ! set OEP shift 
    ! 
    IF ( LOEP ) THEN
       ! set shift to QP-correction of valence band maximum
       OEP%ESHIFT=0
       DO ISP = 1, WDES%NCDIJ
         ! below this variable is set to 0, 
         ! so slightly abuse this variable to store last maximum
         E%EBANDSTR  =-1000
         DO NK1= 1, WDES%NKPTS
            DO I= 1, WDES%NB_TOTK(NK1,ISP)
            IF (W%FERTOT(I,NK1,ISP)>0.4_q .AND.  REAL(W%CELTOT(I,NK1,ISP),q)>E%EBANDSTR) THEN
               E%EBANDSTR=W%CELTOT(I,NK1,ISP)
               IF ( I<=NBANDSGW .AND. ALLOCATED( QP_ZEROS ) ) THEN
                  OEP%ESHIFT(ISP)=(QP_ZEROS(1,I,NK1,ISP)-W%CELTOT(I,NK1,ISP))
               ELSE
                  OEP%ESHIFT(ISP)=-W%CELTOT(I,NK1,ISP)
               ENDIF
            ENDIF
            ENDDO
         ENDDO
       ENDDO

       ! switch off correlation for exact exchange OEP 
       IF ( LEXX ) THEN
          CORR_MAT = 0
          ! also set V_eff(R=inf) = 0 for EXX only 
          OEP%ESHIFT(:) = 0 
       ENDIF
       
       IF (.NOT. ALLOCATED( QP_ZEROS ) ) THEN
          IF (IO%IU0 >0) THEN
             WRITE(IO%IU0,"(' WARNING: OEP shift set to VBM',F14.5)") OEP%ESHIFT
          ENDIF
          IF (IO%IU6>=0) THEN
             WRITE(IO%IU6,"(' WARNING: OEP shift set to VBM',F14.5)") OEP%ESHIFT
          ENDIF
       ELSE   
          IF (IO%IU0 >0) THEN
             WRITE(IO%IU0,"(' OEP shift',F14.5)") OEP%ESHIFT
          ENDIF
          IF (IO%IU6>=0) THEN
             WRITE(IO%IU6,"(' OEP shift',F14.5)") OEP%ESHIFT
          ENDIF
       ENDIF
    ENDIF
  !----------------------------------------------------------------------- 
  ! update Green function and orbitals
  !----------------------------------------------------------------------- 
    MP2SUM =0       ! change of number of electrons
    TR_GAMMA = 0    ! Trace of density matrix 
    E%EBANDSTR = 0
    E%EENTROPY = 0

    CALL ALLOCATE_GREENS_MAT_GDEF(GDES_MAT, GAMMA, 1, 1, 1)

    DO ISP=1,WDES%ISPIN
    DO NK1=1,NKPTS_IRZ
       ! calculate change in Trace of density matrix from correlation effects
       ! CORR_MAT contains (in orbital basis) either G_0 Sigma G_0 or G_0 Sigma G (if LGW0=.TRUE.)
       ! integrated over omega, ie the change of the density matrix
       ! R stores the diagonal elements
       CALL DETERMINE_DIAGONALE_GDEF(WDES%NB_TOTK(NK1,ISP), CORR_MAT(:,NK1,ISP,1), R, GDES_MAT%DESC)
       CALLMPI( M_sum_g(GDES_MAT%COMM_INTRA, R, WDES%NB_TOTK(NK1, ISP) ))
       IF (.NOT.LscQPGW) MP2SUM=MP2SUM+SUM(R(1:WDES%NB_TOTK(NK1,ISP)))*KPOINTS_ORIG%WTKPT(NK1)*WDES%RSPIN

!#define verbose
#ifdef verbose
       CALL DUMP_GREENS_HAM_GDEF( "HF exchange", CHAM_MAT(:,NK1,ISP,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
#endif
       ! for GW0:  GAMMA is set to the HF density matrix in the current basis of eigenvectors (DFT in iteration 1)
       ! for G0W0: GAMMA is set to H_HF / abs(epsilon_a - epsilon_i)
       !           this is only usefull for first iteration
       ! for QPGW: GAMMA is simply set to HF matrix
       GAMMA(:,1,1,1)=CHAM_MAT(:,NK1,ISP,1)
       IF (LGW .AND. ((.NOT. LG0W0 ) .AND.  (.NOT. LscQPGW))) THEN
          ! set GAMMA to Hartree-Fock density matrix, and add mean field Green's function to GO_MAT
          IF (NELM_HF==NELMHF) THEN
             ! dont add mean field part, because so subtraction is done       
             CALL HF_DENSITY_MATRIX( WMEAN, NK1, ISP, GAMMA(:,1,1,1), GDES_MAT, IMAG_GRIDS, EFERMI_MEAN(ISP), &
                  GU_MAT, GO_MAT, LUNOCCUPIED=LHOLEGF)
          ELSE
             ! use density matrix from solution of dyson equation
             ! this is done below, so just clean array here
             CALL HF_DENSITY_MATRIX( WMEAN, NK1, ISP, GAMMA(:,1,1,1), GDES_MAT,&
                  IMAG_GRIDS, EFERMI_MEAN(ISP), LUNOCCUPIED = LHOLEGF )
         ENDIF
       ELSE IF (.NOT. LscQPGW) THEN
          IF (ASSOCIATED(SIGMA_MAT)) THEN 
             CALL DIVIDE_BY_EIGENVALUEDIFFxEIGENVAL( W%WDES%NB_TOTK(NK1,ISP),&
                GAMMA(:,1,1,1), SIGMA_MAT(:,NK1,ISP,1), W%CELTOT(:,NK1,ISP), &
                W%FERTOT(:,NK1,ISP), GDES_MAT%DESC, &
                LFINITE_TEMPERATURE, FERDER(:,NK1,ISP) )
          ENDIF 

          IF( LOEP ) THEN
!CALL DUMP_GREENS_HAM_GDEF( "H_HF", GAMMA(:,NK1,ISP,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
             ! 
             ! form H_KS, trivial in KS basis 
             RR = -W%CELTOT(:,NK1,ISP) 
             IF ( LKINETIC ) RR = 0
!IF(IO%IU0>=0 ) WRITE(*,'("-diag(H_KS)", 16F10.4)')RR(10:25)
             CALL  ADD_TO_DIAGONALE_REAL( WDES%NB_TOTK(NK1,ISP), GAMMA(:,1,1,1), RR, GDES_MAT%DESC)
!CALL DUMP_GREENS_HAM_GDEF( "H_HF-H_KS", GAMMA(:,NK1,ISP,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0, WDES%COMM%NODE_ME+90)
             ! gamma = ( H_mean - H_KS)  (f_i - f_a) / (epsilon_a - epsilon_i)
             CALL DIVIDE_BY_EIGENVALUEDIFF( W%WDES%NB_TOTK(NK1,ISP), &
                GAMMA(:,1,1,1), W%CELTOT(:,NK1,ISP)-W%EFERMI(ISP), W%FERTOT(:,NK1,ISP),&
                GDES_MAT%DESC, LFINITE_TEMPERATURE, FERDER(:,NK1,ISP) )
          ELSE
             ! add in term from entropy - delta T S =  sum_n (mu - epsilon_n) delta f_n
             ! alternatively include term in RPA_ADD_TO_POT in ilinear_response.F
             IF (LRPAFORCE) THEN
                RR = W%EFERMI(ISP)-W%CELTOT(:,NK1,ISP)
                CALL  ADD_TO_DIAGONALE_REAL( WDES%NB_TOTK(NK1,ISP), GAMMA(:,1,1,1), RR, GDES_MAT%DESC)
                ! add  (mu- e_n) e_n df/de from the diagonal of SIGMA_MAT
                IF (ASSOCIATED(SIGMA_MAT)) THEN
                   RR = RR(:)*REAL(W%CELTOT(:,NK1,ISP),q)*FERDER(:,NK1,ISP)
                   CALL  ADD_TO_DIAGONALE_REAL( WDES%NB_TOTK(NK1,ISP), SIGMA_MAT(:,NK1,ISP,1), RR, GDES_MAT%DESC)
                ENDIF
             ENDIF
             ! gamma = G_DFT + G_DFT (H_mean- H_DFT) G_DFT 
             ! in first iteration gamma = H_mean (f_a - f_i) / (epsilon_a - epsilon_i)
             CALL DIVIDE_BY_EIGENVALUEDIFF( W%WDES%NB_TOTK(NK1,ISP), &
                GAMMA(:,1,1,1), W%CELTOT(:,NK1,ISP), W%FERTOT(:,NK1,ISP), &
                GDES_MAT%DESC, LFINITE_TEMPERATURE, FERDER(:,NK1,ISP) )
             ! add unperturbed density matrix  1/ pi int_0^infty   G^W_0 d w
             CALL  ADD_TO_DIAGONALE_REAL( WDES%NB_TOTK(NK1,ISP), GAMMA(:,1,1,1), W%FERTOT(:,NK1,ISP), GDES_MAT%DESC)
          ENDIF
          IF ( LFINITE_TEMPERATURE .AND. LRPAFORCE ) CORR_MAT(:,NK1,ISP,1)=CORR_MAT(:,NK1,ISP,1) + GAMMA(:,1,1,1)
       ENDIF

#ifdef verbose
       CALL DUMP_GREENS_HAM_GDEF( "exchange density matrix", GAMMA(:,1,1,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
       CALL DUMP_GREENS_HAM_GDEF( "correlation density matrix", CORR_MAT(:,NK1,ISP,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
       IF (ASSOCIATED(GO_MAT)) &
       CALL DUMP_GREENS_HAM_GDEF( "green fun GO_MAT", GO_MAT(:,NK1,ISP,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
#endif

       ! add mean field (Hartree-Fock) and correlation contributions of the density matrix
       IF ( LOEP .OR. .NOT. LFINITE_TEMPERATURE ) THEN
          GAMMA(:,1,1,1)=GAMMA(:,1,1,1)+CORR_MAT(:,NK1,ISP,1)
       ELSE
          !finite temperature already includes all contributions
          GAMMA(:,1,1,1)=CORR_MAT(:,NK1,ISP,1)
       ENDIF
       ! overwrite the diagonale elements by the diagonale of the self-energy (Kotani-Schilfgaarde)
       IF (LscQPGW) THEN
          CALL  SET_DIAGONALE_REAL( WDES%NB_TOTK(NK1,ISP), GAMMA(:,1,1,1), REAL(SIGMA(:,NK1,ISP),q), GDES_MAT%DESC)
       ENDIF
  !----------------------------------------------------------------------- 
  ! construct rotation matrix
  ! NATURALO:   0   standard diagonalization, update FERTOT to gamma 
  !             1   rotation within unoccupied, no update of FERTOT
  !            -N    negative value also preserve ABS(NATURALO) unoccupied
  !                 orbitals
  !             2   standard diagonalization, not update FERTOT
  !             3   same as 1, store eigenvalues in WAVECAR file
  ! if 10 is added to NATURALO an  iterative diagonalization is used
  !-----------------------------------------------------------------------
    IF (.NOT. LRPAFORCE) THEN
       IF (LscQPGW) THEN
#ifdef verbose
          CALL DUMP_GREENS_HAM_GDEF( "exchange + correlation matrix", GAMMA(:,1,1,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
          CALL DETERMINE_DIAGONALE_GDEF(WDES%NB_TOTK(NK1,ISP), GAMMA(:,1,1,1), R, GDES_MAT%DESC)
          CALLMPI( M_sum_g(GDES_MAT%COMM_INTRA, R, WDES%NB_TOTK(NK1, ISP) ))
!          IF (IO%IU0>=0) THEN
!             WRITE(*,'(8F8.4)') REAL(R(1:24),q)
!          ENDIF
#endif

          CALL PDSYEV_ZHEEVD_DESC( GAMMA(:,1,1,1), RR(1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT%DESC, GDES_MAT%COMM_INTRA)
          CALLMPI( M_bcast_d_from(GDES_MAT%COMM_INTER , RR(1), SIZE(RR), 1))
          ! store eigenvalues in W%CELTOT
          W%CELTOT(1:W%WDES%NB_TOTK(NK1,ISP),NK1,ISP)=RR(1:W%WDES%NB_TOTK(NK1,ISP))
          ! TODO maybe call DENSTA to update occupancies
       ELSE
        IF (NATURALO==1 .OR. NATURALO==3 .OR. NATURALO==11 .OR. NATURALO==13 .OR. NATURALO < 0 ) THEN
          ! last "occupied" orbital
          N2= LAST_FILLED_XI_NOMOD( W, NK1, ISP)
          IF (NATURALO < 0) N2=N2+ABS(NATURALO)
          ! remove the occ-unocc blocks
          CALL CLEAR_OCCUNOCC_UNOCCOCC(W%WDES%NB_TOTK(NK1,ISP), GAMMA(:,1,1,1), N2, GDES_MAT%DESC)
          ! and set the occ-occ block to identity matrix to get similar set-up as MP2 NOs procedure
          CALL SET_OCCOCC_INCREASING(W%WDES%NB_TOTK(NK1,ISP), GAMMA(:,1,1,1), N2, GDES_MAT%DESC)
        ENDIF
        IF (NATURALO>=10) THEN
        ! use iterative diagonalization keeping original orbital ordering
        ! for insulators one iteration should be identical to original version documented out below
          CALL ITERATIVE_DIAG( W%WDES%NB_TOTK(NK1,ISP), NK1, ISP, GAMMA(:,1,1,1), R, GDES_MAT, ITER=5)
          RR=R
          CALLMPI( M_bcast_d_from(GDES_MAT%COMM_INTER , RR(1), SIZE(RR), 1))
        ELSE
          ! change sign of GAMMA(:,1,1,1) (eigenvalue ordering of diagonalization routine)
          GAMMA(:,1,1,1)=-GAMMA(:,1,1,1)
          CALL PDSYEV_ZHEEVD_DESC( GAMMA(:,1,1,1), RR(1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT%DESC, GDES_MAT%COMM_INTRA)
          CALLMPI( M_bcast_d_from(GDES_MAT%COMM_INTER , RR(1), SIZE(RR), 1))
          RR=-RR   ! undo previous sign change
          GAMMA(:,1,1,1)=-GAMMA(:,1,1,1)
        ENDIF
        ! store diagonal of density matrix in FERTOT for later use
        IF (NATURALO==0 .OR. NATURALO==10) W%FERTOT(1:W%WDES%NB_TOTK(NK1,ISP),NK1,ISP)=MAX(RR(1:W%WDES%NB_TOTK(NK1,ISP)),0.0_q)
        IF (NATURALO==3 .OR. NATURALO==13) THEN
            W%FERTOT(1:W%WDES%NB_TOTK(NK1,ISP),NK1,ISP)=MAX(RR(1:W%WDES%NB_TOTK(NK1,ISP)),0.0_q)
            W%FERTOT(1:N2,NK1,ISP)=1
        ENDIF
        TR_GAMMA(ISP)=TR_GAMMA(ISP)+SUM(RR(1:WDES%NB_TOTK(NK1,ISP)))*KPOINTS_ORIG%WTKPT(NK1)*WDES%RSPIN !get total density
        ! note, for LFINITE_TEMPERATURE = .TRUE. CORR_MAT is the full density matrix 
        ! containing also the HF part! build Tr( Gamma * ( T + V_ext + V_h + V_x) ) in HF basis ij
        DO I = 1, WDES%NB_TOTK(NK1,ISP) 
           DTRGHF = DTRGHF - MAX(RR(I),0.0_q)*WMEAN%CELTOT(I, NK1, ISP)*KPOINTS_ORIG%WTKPT(NK1)*WDES%RSPIN
        ENDDO
        ! determine entropy term of electrons and holes -TS =  Tr( GAMMA * Ln( GAMMA ) + (1-GAMMA)*Ln( 1-GAMMA ) ) 
        DO I = 1 , W%WDES%NB_TOTK(NK1,ISP)
           IF ( RR(I) > 0 .AND. RR(I) < 1) E%EENTROPY = E%EENTROPY&
              +WDES%RSPIN*KPOINTS_ORIG%WTKPT(NK1)*KPOINTS%SIGMA*&
              ( RR(I)*LOG(RR(I)) + (1-RR(I))*LOG(1-RR(I)) )
        ENDDO
       ENDIF

       ! bcast from root group to all other groups to avoid accumulation of rounding errors
       CALLMPI( M_bcast_g_from(GDES_MAT%COMM_INTER, GAMMA(1,1,1,1), SIZE(GAMMA,1), 1))
       IF( LFCIINPUT ) CALL DUMP_NATURALO_ROT_MATRIX(1000, W, NK1, ISP, GAMMA(:,1,1,1), GDES_MAT, IO%IU0)
#ifdef verbose
       CALL DUMP_GREENS_HAM_GDEF( "rotation matrix", GAMMA(:,1,1,1),W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
#endif
       IF ( .NOT. LOEP ) THEN
          ! restore original mean field HF matrix
          IF (LGW .AND. ((.NOT. LG0W0 ) .AND.  (.NOT. LscQPGW))) CALL HF_MATRIX( WMEAN, NK1, ISP, CHAM_MAT(:,NK1,ISP,1), GDES_MAT )
          ! calculate change of HF energy in rotated basis 
          CALL HF_SINGLES(  W, NK1, ISP, GAMMA(:,1,1,1), CHAM_MAT(:,NK1,ISP,1), GDES_MAT, E%EBANDSTR)
       ENDIF
    ENDIF
       ! store the rotation matrix for later use
       CHAM_MAT(:,NK1,ISP,1)=GAMMA(:,1,1,1)
#ifdef verbose
       IF (IO%IU0>=0) WRITE(IO%IU0,'("occup ", 10F10.6)') RR(1:MIN(10,W%WDES%NB_TOTK(NK1,ISP)))
       IF (IO%IU0>=0) WRITE(IO%IU0,'("occup ", 10F10.6)') RR(W%WDES%NB_TOTK(NK1,ISP)-9:W%WDES%NB_TOTK(NK1,ISP))
#endif
    ENDDO  !k-points
    ENDDO  !spin
    CALL DEALLOCATE_GREENS_MAT_GDEF(GDES_MAT, GAMMA)
    IF ( LFCIINPUT ) CLOSE( 1000 )   

    IF (IO%IU0>=0.AND. (NELMGW==1 .OR. .NOT.LFINITE_TEMPERATURE)) THEN
       IF (.NOT. LRPAFORCE) THEN
         IF ( LHOLEGF ) THEN
            WRITE(IO%IU0,'(" Trace of density matrix (holes)       ",2F20.10)') SUM(TR_GAMMA(:))
         ELSE     
            WRITE(IO%IU0,'(" Trace of density matrix (electrons)   ",2F20.10)') SUM(TR_GAMMA(:))

         ENDIF    
       ENDIF
       IF ( .NOT. LFINITE_TEMPERATURE ) &
          WRITE(IO%IU0,'(" correlated contrib. to density matrix ",2F20.10)') MP2SUM
    ENDIF
    !
    ! add correlation term to Galitskii-Migdal internal energy 
    !
    TOTEGM = E%ETRGHF + E%ECGWGM/2
    !
    ! zero-temperature case computes only Galitskii-Migdal energy 
    !
    IF ( .NOT. LFINITE_TEMPERATURE ) THEN
       TOTEN = 0 
       E%EENTROPY = 0 
    ENDIF 

    IF (IO%IU6>=0) THEN
    IF (.NOT. LRPAFORCE) THEN
       IF ( .NOT. (LG0W0 .OR. LACFDT .OR. .NOT. LFINITE_TEMPERATURE )) THEN
          CALL DUMP_DIAGRAMMATIC_ENERGIES( "Energies after update of density matrix ", &
          E, TOTEN, INFO, IO, ECHEMN=-DOT_PRODUCT(CMPLX(EFERMI_MEAN,0._q),TR_GAMMA)) 
       ELSE
          CALL DUMP_DIAGRAMMATIC_ENERGIES_RPA( "Energies after update of density matrix ", &
          E, TOTEN, INFO, IO) 
       ENDIF
       IF ( LHOLEGF ) THEN
          WRITE(IO%IU6,'("  Trace of density matrix (holes)       ",2F20.10)') SUM(TR_GAMMA(:)) 
       ELSE     
          WRITE(IO%IU6,'("  Trace of density matrix (electrons)   ",2F20.10)') SUM(TR_GAMMA(:)) 
       ENDIF    
       IF ( .NOT. LFINITE_TEMPERATURE ) &
          WRITE(IO%IU6,'("  correlated contrib. to density matrix ",2F20.10)') MP2SUM
    ELSE
       IF ( .NOT. LFINITE_TEMPERATURE ) &
          WRITE(IO%IU0,'("  correlated contrib. to density matrix ",2F20.10)') MP2SUM
    ENDIF
    WRITE(IO%IU6,'(/)')
    ENDIF

    CALL STOP_TIMING("GREENS",IO%IU6,"MAKEROT")
    CALL STOP_TIMING("GWLOOP",IO%IU6)

    CALL DUMP_ALLOCATE(IO%IU6); IF (IO%IU6>=0) CALL WFORCE(IO%IU6)
  !----------------------------------------------------------------------- 
  ! finally rotate orbitals
  !-----------------------------------------------------------------------
    ! update the one-electron orbitals, and rotate Greens function & density matrix
    IF (LRPAFORCE) THEN
       CALL STORE_FORCE_RPA_PULAY(T_INFO,FRNL)
       IF (ASSOCIATED(SIGMA_MAT)) THEN
          CALL STORE_RPA_DENSITY_MATRIX( W,  WDES, CHAM_MAT, GDES_MAT%COMM_INTRA,  GDES_MAT%DESC , SIGMA_MAT)
       ELSE
          CALL STORE_RPA_DENSITY_MATRIX( W,  WDES, CHAM_MAT, GDES_MAT%COMM_INTRA,  GDES_MAT%DESC )
       ENDIF
    ELSE
       !
       ! form natural orbitals and store them in W
       !
       CALL ROTATE_ORBITALS( W, WDES, CHAM_MAT, GDES_MAT, CORR_MAT, GU_MAT, GO_MAT)
       !
       ! determine OEP potential
       !
       CALL OEP_POTENTIAL( CHI0, GDES, S2E, WGW, IMAG_GRIDS, &
          W, WDES, INFO, E, &
          GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
          LATT_CUR, SYMM, T_INFO, P, &
          LMDIM, CDIJ, CQIJ, SV, &
          CHDEN, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX, &
          DENCOR, CVTOT, CSTRF, &
          MIX%AMIX, OEP, IO )

       ! update wavefunction with current OEP potential 
       CALL UPDATE_W_OEP(W, SV, CDIJ)
    ENDIF
    CALL STOP_TIMING("GREENS",IO%IU6,"ROTATE")
!    WDES%NB_TOTK=96
    ! the re-orthogonalization is not required, since all transformations
    ! are unitary, documenting in the call should yield exactly identical results
!    CALL ORTHCH(WDES,W, WDES%LOVERL, LMDIM, CQIJ)

#ifdef debug
    IF (ASSOCIATED(GO_MAT)) THEN
       DO ISP=1,WDES%ISPIN
          DO NK1=1,WDES%NKPTS
             CALL DUMP_GREENS_HAM_GDEF( "green fun GO after rot", GO_MAT(:,NK1,ISP,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
             CALL DUMP_GREENS_HAM_GDEF( "green fun GU after rot", GU_MAT(:,NK1,ISP,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
          ENDDO
       ENDDO
    ENDIF
#endif
  ENDDO iterhf

  CALL DEALLOCATE_GREENS_MAT_GDEF(GDES_MAT, CHAM_MAT)
  CALL DEALLOCATE_GREENS_MAT_GDEF(GDES_MAT, CORR_MAT)
  IF (ASSOCIATED(SIGMA_MAT)) CALL DEALLOCATE_GREENS_MAT_GDEF(GDES_MAT, SIGMA_MAT)

  IF (( LGW .AND. ((.NOT. LG0W0 ) .AND.  (.NOT. LscQPGW))) ) THEN
    IF (IO%IU0>=0) THEN
       FLUSH( 17 )    
       FLUSH( IO%IU0 )    
       IF ( NELM == 1 ) WRITE(17,303)
       IF ( NELM == 1 ) WRITE(IO%IU0,303)
       WRITE(17, 304)     NELM, TOTEN, TOTEN-FLWTOT_OLD, DTRGHF, TOTEN-TOTEGM, ABS(SUM(TR_GAMMA(:))),EFERMI_VAR  
       WRITE(IO%IU0, 304) NELM, TOTEN, TOTEN-FLWTOT_OLD, DTRGHF, TOTEN-TOTEGM, ABS(SUM(TR_GAMMA(:))),EFERMI_VAR  
    ENDIF

    IF ( ABS(SUM(TR_GAMMA(:))-INFO%NELECT)>ABS(SUM(TR_GAMMA(:))-INFO%NELECT+1).OR.&
         ABS(SUM(TR_GAMMA(:))-INFO%NELECT)>ABS(SUM(TR_GAMMA(:))-INFO%NELECT-1) ) THEN
       IF (IO%IU0>=0) THEN
          WRITE(IO%IU0,*) ' WARNING: number of electrons is not conserved ', INFO%NELECT, ABS(SUM(TR_GAMMA(:)) )
       ENDIF
    ENDIF
  ENDIF
  303 FORMAT('          N  F(LW)              dF            dTr[G*HF]     F(LW)-E(GM)      Tr( G )    d mu')
  304 FORMAT(' Dyson:',I4,1X,E16.9,2X,E12.5,2X,E12.5,2X,E12.5,2X,F12.6,2X,E12.5)
  ! 
  ! break condition
  ! 
  IF ( INFO%EDIFF > ABS( FLWTOT_OLD - TOTEN))  THEN
     EXIT scf
  ENDIF
  FLWTOT_OLD = TOTEN

  ENDDO scf

  IF (LGW .AND. ((.NOT. LG0W0 ) .AND.  (.NOT. LscQPGW))) THEN
     CALL DEALLOCATE_GREENS_MAT(GDES_MAT, SIGMAW_MAT)
     IF (ASSOCIATED(GO_MAT)) THEN
       CALL DEALLOCATE_GREENS_MAT_GDEF(GDES_MAT, GO_MAT)
       CALL DEALLOCATE_GREENS_MAT_GDEF(GDES_MAT, GU_MAT)
     ENDIF
  ENDIF

  CALL DEALLOCW_NOPLANEWAVE(WMEAN)
  CALL DEALLOCW_NOPLANEWAVE(WHF_IN_KS)
!=======================================================================
  ENDIF rpa_or_gw
!=======================================================================
   
    ! response function stores also W, so needs to be allocated during GW0
    ! self-consistent run
    DO NQ_COUNTER=1,SIZE(CHI_QPT)
       CALL DEALLOCATE_RESPONSEFUN(CHI_QPT(NQ_COUNTER))
    ENDDO
    DEALLOCATE(CHI_QPT)
    CALL DEALLOCATE_RESPONSEFUN(CHI_TAU(1))
    DEALLOCATE(CHI_TAU)

    ! destroy  sclapack descriptor
    IF (ASSOCIATED(GDES_MAT)) THEN
       CALL DESTROY_GDES_MAT( GDES_MAT )
    ENDIF
 
    IF( LOEP ) THEN
       CALL DEALLOCATE_RESPONSEFUN(CHI0(1))
       DEALLOCATE(CHI0)
    ENDIF

    IF (IO%LWAVE .AND. .NOT. ( LACFDT .OR. LCRPA ) ) THEN
       CALL OUTWAV(IO, WDES, W, LATT_CUR, W%EFERMI(1) )
    ENDIF

    IF (IO%LWAVEH5 .AND. .NOT. ( LACFDT .OR. LCRPA )) THEN
#ifdef VASP_HDF5
       IH5ERR = VH5_FILE_DELETE('vaspwave.h5', 112)
       IF (IH5ERR /= 0) THEN
         WRITE(*,*) "ERROR deleting old wave file, ierr = ", IH5ERR
       ENDIF
       IH5ERR = VH5_FILE_CREATE_OR_OVERWRITE('vaspwave.h5', IH5WAVEFILEID)
       IF (IH5ERR == 0) THEN
         CALL VH5_WRITE_WAVEFUNCTIONS(IH5WAVEFILEID, IO, WDES, W, LATT_INI, W%EFERMI(1))
         IH5ERR = VH5_FILE_CLOSE_WRITING(IH5WAVEFILEID)
       ELSE
         WRITE(IO%IU0,*) "Error: could not create vaspwave.h5, errorcode: ",IH5ERR
       ENDIF
#endif
    ENDIF
!-----------------------------------------------------------------------
! some deallocation
!-----------------------------------------------------------------------
    IF (SYMM%ISYM>=0) THEN
       CALL DEALLOCATE_KPOINTS_TRANS(RES_KPTS_TRANS)
    ENDIF

    CALL DEALLOCATE_IRZ_MAP(S2E)

    !destroy imaginary grid handle
    CALL DEALLOCATE_IMAG_GRID_HANDLE(IMAG_GRIDS)
    DEALLOCATE(WGW%VKPT, WGW%WTKPT)
    DEALLOCATE(WGW, GRIDWGW)
    CALL DEALLOCATE_OEP_HANDLE(OEP)

    PROFILING_STOP('calculate_xi_real_gwrk')

    CONTAINS

   SUBROUTINE UPDATE_W_OEP( W, SV, CDIJ )
      USE setexm
      TYPE( wavespin ) W
      RGRID       SV(DIMREAL(WDES%GRID%MPLWV),WDES%NCDIJ)
      OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      !local
      REAL(q) DESUM1, RMS, EF
      INTEGER :: ICOUEV, NSIM 
   
      ICOUEV = 0 
      IF( .NOT. LOEP ) RETURN 
      ! Davidson without HF contribution
      NSIM=WDES%NSIM
#ifdef MPI
      NSIM=((WDES%NSIM+WDES%COMM_INTER%NCPU-1)/WDES%COMM_INTER%NCPU)*WDES%COMM_INTER%NCPU
#endif
       
      CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
           LMDIM,CDIJ,CQIJ, 3,SV,T_INFO,P,IO%IU0,E%EXHF)
      IF (IO%IU0>0) WRITE(IO%IU0,*) 'calling exact diagonalization EDDIAG_EXACT'
      IF (IO%IU0>=0) WRITE(17, 200)  NELM,0._q,0._q,DESUM1,ICOUEV,RMS
      IF (IO%IU0>=0) &
       WRITE(IO%IU0, 200)  NELM,0._q,0._q,DESUM1,ICOUEV,RMS

200   FORMAT('DIA: ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5, &
        &       I6,'  ',E10.3)

      ! update occupancies (well this is done below as well, but we need it for eigenvalues)
      !keep fermi energy as it is for finite temperature calculations 
      CALL DENSTA_SPIN( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
        INFO%NUP_DOWN,  E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
        NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)

      IF (IO%IU6>=0) WRITE(IO%IU6,2202) EFERMI
      CALL WRITE_EIGENVAL_NBANDS( WDES, W, IO%IU6, NBANDSGW)

2202  FORMAT(' E-fermi : ', F8.4)
   
      IF (IO%LOPTICS) THEN
         CALL START_TIMING("G")
        ! VASP onboard optics
         CALL PEAD_RESETUP_WDES(WDES,GRID,KPOINTS,LATT_CUR,LATT_INI,IO)
         EF=SUM(EFERMI)/SIZE(EFERMI)
         CALL LR_OPTIC( &
            P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
            T_INFO,INFO,IO,KPOINTS,SYMM,GRID,GRID_SOFT, &
            GRIDC,GRIDUS,C_TO_US,SOFT_TO_C, &
            CHTOT,DENCOR,CVTOT,CSTRF, &
            CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM, &
            CHDEN,SV,LMDIM,IRDMAX,EF,NEDOS, & 
            LSTORE=.TRUE., LPOT=.FALSE.)
         CALL STOP_TIMING("G",IO%IU6,'OPTICS')
      ENDIF
   END SUBROUTINE UPDATE_W_OEP

  END SUBROUTINE CALCULATE_XI_REAL_GWRK

#endif

END MODULE chi_GG
