#include "symbol.inc"

#if defined(gammareal) && defined(RPAgamma)
#undef RPAgamma
#endif
!*********************************************************************
!
! this module implements the calculation of the response
! function from the Green's function
!   G(r',r, tau ) G* + (r',r, tau)
! where G is the Green's function for the unoccupied states
!   G (r',r, tau) = sum_a <r'|a> <a| r> exp ( -(e_a- e_F) tau )
! and G+ is the Green's function for the occupied states
!   G+(r',r, tau) = sum_a <r'|i> <i| r> exp (  (e_i- e_F) tau )
! it is used by chi_GG
! and allows to calculate the response function with linear scaling
! in the number of k-points
!
! note that the resulting response functions are conjugated
! compared to the definition used in chi_base.F
!
!*********************************************************************

MODULE chi_super
  USE prec
  USE greens_real_space
  USE GG_base
  USE greens_orbital
  USE rpa_force
  USE chi_glb 
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: CALCULATE_XI_REAL


CONTAINS

!***********************************************************************
!
! main calling routine
!
!***********************************************************************

  SUBROUTINE CALCULATE_XI_REAL( & 
          HAMILTONIAN, P, WDES, NONLR_S, NONL_S, W, LATT_CUR, LATT_INI, &
          GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, E, & 
          CHTOT, CHTOTL, DENCOR, CVTOT, CSTRF, IRDMAX, &
          T_INFO, DYN, INFO, IO, KPOINTS, SYMM, MIX, & 
          LMDIM, CQIJ, CDIJ, CRHODE, N_MIX_PAW, RHOLM, RHOLM_LAST, CHDEN, SV, & 
          EFERMI_IN, TOTEN, NEDOS, DOS, DOSI )
    USE base
    USE pseudo_struct_def
    USE mpimy
    USE constant
    USE wave
    USE wave_high
    USE ini
    USE fileio, ONLY: OUTWAV
    USE hamil, ONLY: BANDSTRUCTURE_ENERGY
    USE hamil_struct_def
    USE subrot, ONLY: EDDIAG
    USE acfdt_GG
    USE chi_base
    USE minimax
    USE pade_fit
    USE dmft
    USE pead, ONLY : PEAD_RESETUP_WDES
    USE mymath, ONLY : GREATEST_COMMON_DIVISOR
    USE tutor, ONLY: vtutor, isError, isAlert, ACFDTRnotSupported, &
        FermiFiniteTemp, OmegaGrid140, TauParCPU, OmegaPar, &
        TauPar, chiGG_NOmega, OmegaGridOmegaTL, argument, &
        ! OEP related
        LPOTOK  
    USE mlr_optic, ONLY: LR_OPTIC_DONE
    ! required for OEP
    USE us, ONLY: SET_CHARGE
    USE density_of_states, ONLY: DENSTA_DER
    use string, only: str
#ifdef VASP_HDF5
    USE vhdf5
#endif
    USE c2f_interface, ONLY: TIMING
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
    TYPE(supercell),POINTER :: SUPER=>NULL()
#ifdef VASP_HDF5
    INTEGER IH5ERR
#endif

#ifndef scaLAPACK
    CALL vtutor%write(isError, ACFDTRnotSupported)
  END SUBROUTINE CALCULATE_XI_REAL
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

    TYPE (correlation), POINTER :: COR=>NULL()
    LOGICAL LGAMMA
    TYPE (responsefunction), POINTER :: CHI_TAU(:)=>NULL()  ! response function in tau domain
    TYPE (responsefunction), POINTER :: CHI=>NULL()         ! response function in nu domain
    TYPE (responsefunction), POINTER :: CHI_QPT(:)=>NULL()  ! response function in nu domain
    TYPE (responsefunction), POINTER :: CHI0(:)             ! response function at w=0
    !!!need to figure out which one we need in the WORK state
    TYPE (greensf), POINTER :: GU(:)=>NULL()           ! Green function for unoccupied orbitals in imag. time  G( tau)
    TYPE (greensf), POINTER :: GO(:)=>NULL()           ! Green function for occupied orbitals in imag. time G(-tau)
    TYPE (greensf), POINTER :: SIGMAU(:)=>NULL()       ! Self-energy G( tau) W(tau)  (invokes only GU)
    TYPE (greensf), POINTER :: SIGMAO(:)=>NULL()       ! Self-energy G(-tau) W(tau)  (invokes only GO)

    TYPE (wavedes), POINTER :: WGW=>NULL()
    TYPE (wavedes1)         :: WGWQ
    INTEGER                 :: NQ, NQ_COUNTER, NK1, NK2, ISP
    TYPE (grid_3d), POINTER :: GRIDWGW=>NULL()
    TYPE (screened_2e_handle) :: S2E

    TYPE (imag_grid_handle) :: IMAG_GRIDS      ! imaginary grids handle
    TYPE (wavespin) WHF
    INTEGER :: N1, N2, NCPU, ierror
    REAL(q)    PAR(1,1,1,1,WDES%NCDIJ),DOSPAR(1,1,1,WDES%NCDIJ)
    GDEF, POINTER       :: SIGMAO_DIAG(:)=>NULL()      ! diagonal of <i| Sigma_o(tau) | i> of self-energy
    GDEF, POINTER       :: SIGMAU_DIAG(:)=>NULL()      ! diagonal of <i| Sigma_u(tau) | i> of self-energy
    GDEF, POINTER, CONTIGUOUS :: SIGMAO_MAT(:,:,:,:)=>NULL() ! matrix <i| Sigma_o(tau) | a> stored distributed
    GDEF, POINTER, CONTIGUOUS :: SIGMAU_MAT(:,:,:,:)=>NULL() ! matrix <i| Sigma_u(tau) | a> stored distributed
    GDEF, POINTER, CONTIGUOUS :: CHAM_MAT(:,:,:,:)=>NULL()   ! Hamilton matrix T + V_local + V_x or density matrix
    GDEF, POINTER, CONTIGUOUS :: GAMMA(:,:,:,:)=>NULL()      ! density matrix for one k-point
    GDEF, POINTER, CONTIGUOUS :: CORR_MAT(:,:,:,:)=>NULL()   ! correlation contribution to density matrix
    GDEF, POINTER, CONTIGUOUS :: SIGMA_MAT(:,:,:,:)=>NULL()  ! int G_0(w) (\Sigma(w)+H_HF) G_0(w) w dw
    GDEF, POINTER, CONTIGUOUS :: GO_MAT(:,:,:,:)=>NULL()     ! matrix <i| G_o(tau) | a> stored distributed
    GDEF, POINTER, CONTIGUOUS :: GU_MAT(:,:,:,:)=>NULL()     ! matrix <i| G_u(tau) | a> stored distributed
!#define MatsubaraIntegration
#ifdef MatsubaraIntegration
    GDEF, POINTER, CONTIGUOUS :: GMATU_STORE(:,:,:,:)=>NULL() ! matrix <i| Sigma_o(tau) | a> stored distributed
    GDEF, POINTER, CONTIGUOUS :: GMATO_STORE(:,:,:,:)=>NULL() ! matrix <i| Sigma_u(tau) | a> stored distributed
#endif
    COMPLEX(q), POINTER, CONTIGUOUS :: SIGMAW_MAT(:,:,:,:)=>NULL() ! frequency dependent self-energy 
    COMPLEX(q), POINTER :: SIGMA_COS(:,:,:,:)=>NULL()  ! cosine transformed self-energy
    COMPLEX(q), POINTER :: SIGMA_SIN(:,:,:,:)=>NULL()  ! sin transformed self-energy
    GDEF, POINTER       :: SIGMAO_TAU(:,:,:,:)=>NULL() ! diagonal of <i| Sigma_o(tau) | i> of self-energy
    GDEF, POINTER       :: SIGMAU_TAU(:,:,:,:)=>NULL() ! diagonal of <i| Sigma_u(tau) | i> of self-energy
    INTEGER         :: NTAU_ROOT
    INTEGER :: NOMEGA_GRID                     ! number of frequency points in grids
    INTEGER :: I, J
    REAL(q) :: TAU_TMP                         ! auxillary: temporary storage
                                               ! for local time point
    INTEGER :: IU6_TMP, IU6_MEM
    TYPE (skpoints_trans) ::  RES_KPTS_TRANS   ! handle to generate response function at k-points outside the IRZ
    TYPE (skpoints_trans) ::  KPTS_TRANS       ! handle to generate function at k-points outside the IRZ
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
    INTEGER :: IRDMAA             ! temporary
    INTEGER :: NELM               ! loop counter in GW part
    INTEGER :: NELM_HF            ! loop counter in HF part
    LOGICAL :: LUSESUPER          ! use the supercell routines
    GDEF    :: R(WDES%NB_TOT)
    REAL(q) :: RR(WDES%NB_TOT)
    TYPE (greens_mat_des), POINTER :: GDES_MAT=>NULL()
    REAL(q) :: FRNL(3,T_INFO%NIONS)
    REAL(q)    TOTEGM, DTRGHF, FLWTOT_OLD
    LOGICAL, EXTERNAL :: CALCULATE_RPA_FORCES
    COMPLEX(q):: TR_GAMMA(WDES%ISPIN)
    ! correlated wavefunctions for CRPA calculations
    TYPE (wavespin)                  :: WCORR
    TYPE (wavespin)                  :: WCORRHF 
    ! progress bar 
    INTEGER                          :: IPROGRESS, IGOAL
    INTEGER                          :: IPROGRESS_GW, IGOAL_GW
    ! OEP related 
    TYPE( oep_handle ), POINTER      :: OEP => NULL()
    TYPE( wavespin ) :: WOEP 
    REAL(q) ::  FERDER( W%WDES%NB_TOT, W%WDES%NKPTS, W%WDES%ISPIN ) 

    LOGICAL :: LKINETIC=.FALSE. ! same as in chi_base.F for OEP
!    LOGICAL :: LKINETIC=.TRUE.
    REAL(q) :: PSS ! proportional set size 
#ifdef MPI
    IF (W%WDES%COMM_KINTER%NCPU.NE.1) THEN
       CALL vtutor%error("CALCULATE_XI_REAL: KPAR>1 not implemented, sorry.")
    END IF
#endif
    PROFILING_START('calculate_xi_real')

    FLWTOT_OLD = 0.0_q

!#define mem_profiling
#ifdef mem_profiling
   CALL GET_PSS( PSS, IO%IU0, IO%IU6, "START:",0 )
#endif

    ! hand brake if user tries to use  non-collinear version 
    IF ( WDES%LNONCOLLINEAR ) THEN
       CALL vtutor%error("non collinear calculations for low-&
          &scaling GW-RPA algorithms are not implemented yet")
    ENDIF

    EFERMI_SAVE = EFERMI_IN
    EFERMI_VAR = 0 

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

! decide whether to use gamma only mode for response functions

    IF (W%WDES%NKPTS == 1 .AND. SUM(W%WDES%VKPT(:,1)*W%WDES%VKPT(:,1))<G2ZERO)  THEN
       LGAMMA = .TRUE.
    ELSE
       LGAMMA = .FALSE.
    ENDIF

! decide whether to use supercell routines for GG and Sigma G or not
#ifdef gammareal
    LUSESUPER=.FALSE.
#else
    LUSESUPER=.TRUE.
#endif
    CALL SET_NBANDSGW( W)

    ! just in case, perform fermi weight update 
    CALL DENSTA_SPIN( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
       INFO%NUP_DOWN,  E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
       NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)

    ! set W%AUXTOT to the derivatives of the fermi weights
    ! mK: AUXTOT cannot be used, because FIRST_XI_EMPTY_LOCAL depends on it
    ! introduce new temporary array 
    ! W%EFERMI must be set 
    CALL DENSTA_DER( W, KPOINTS, FERDER)

    IF ( LFINITE_TEMPERATURE ) THEN
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

! sync orbitals and eigenvalues (just to be sure everything is ok)
    CALL KPAR_SYNC_ALL(WDES,W)
!=======================================================================
! switch off symmetry
!=======================================================================
    CALL NO_SYMMETRY(.TRUE., LGAMMA, SYMM, W, WDES, LATT_CUR, LATT_INI, &
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
#ifdef mem_profiling
    CALL GET_PSS( PSS, IO%IU0, IO%IU6, "AFTER SYMMETRY",0 )
#endif

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
    ! test if GRID_SOFT and response function grid are commensurate for OEP
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

    ! relink W the WHF so that the fft grids of the Hartree Fock routine are used
    WHF=W
    WHF%WDES => WDES_FOCK
    IF (SYMM%ISYM>=0) THEN
       CALL GENERATE_KPOINTS_TRANS(WDES_FOCK%GRID, NKPTS_IRZ, WDES_FOCK, KPOINTS_FULL_ORIG, KPTS_TRANS)

    ENDIF
#ifdef mem_profiling
    CALL GET_PSS( PSS, IO%IU0, IO%IU6, "AFTER TRANS",0 )
#endif
    ! transformation between IrBZ and BZ
    IF (SYMM%ISYM>=0) THEN
       ! this is tricky: reduce the k-points to WDES%NKPTS (for RPAgamma this is smaller than KPOINTS_FULL%NKPTS)
       ! the code never goes beyond WDES%NKPTS when indexing RES_KPTS_TRANS
       ! in principle, this should be fixed around WGW%NKPTS=KPOINTS_FULL%NKPTS (-> WGW%NKPTS=WDES%NKPTS)
       ! however, this would break SETUP_IRZ_MAP; so that routine requires carefull prior revision
       WGW%NKPTS=W%WDES%NKPTS
       CALL GENERATE_KPOINTS_TRANS(WGW%GRID, NKPTS_IRZ, WGW, KPOINTS_FULL_ORIG, RES_KPTS_TRANS)
    ENDIF
    CALL CHECK_FULL_KPOINTS

!=======================================================================
! super cell initialization
!=======================================================================
    CALL INIT_SUPERCELL( SUPER, T_INFO, LATT_CUR, & 
         KPOINTS%NKPX, KPOINTS%NKPY,  KPOINTS%NKPZ, KPOINTS%B , IO%IU6)
#ifdef mem_profiling
    CALL GET_PSS( PSS, IO%IU0, IO%IU6, "INIT_SUPER",0 )
#endif

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

    IF (LHFCALC_GG ) THEN
       NOMEGA_GRID=NOMEGA-1
    ELSE
       NOMEGA_GRID=NOMEGA
    ENDIF

#ifdef mem_profiling
    CALL GET_PSS( PSS, IO%IU0, IO%IU6, "AFTER MINIMAX",0 )
#endif
!=======================================================================
! we split up the original communicator into groups
! the frequencies are distributed between the groups, each group 
! receives  B%NPOINTS_IN_GROUP (~NOMEGA/NOMEGAPAR) points
! frequencies (ideally one) and works only on those
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
          ! for laplace transformed mp2 energy computation:
          ! it will determine the direct MP2 correlation energy, but nonsence for RPA correlation energy
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

    IF (LHFCALC_GG ) THEN
       IF (IO%IU0>0) WRITE(IO%IU0,*) 'internal WARNING in VASP: adding tau=0- as grid point'
       ! note that this is not needed for Finite-temperature 
       ! since 0 point is already contained in grid
       CALL ADD_POINT_IMAG_GRID_HANDLE(IMAG_GRIDS)
    ENDIF

    IF (IO%LOPEN) CALL WFORCE(IO%IU6)
    CALLMPI( M_barrier(WDES%COMM))

    IF (IO%LOPEN) CALL WFORCE(IO%IU6)
#ifdef mem_profiling
    CALL GET_PSS( PSS, IO%IU0, IO%IU6, "DISTRIBUTING IMAG POINTS",0 )
#endif

!=======================================================================
! distribute grid points based on NOMEGAPAR and NTAUPAR 
!=======================================================================
    CALL SET_LOOP_DES( WDES%COMM_KIN, IMAG_GRIDS, 'boso', IO )
    CALL SET_LOOP_DES( WDES%COMM_KIN, IMAG_GRIDS, 'ferm', IO )
    CALL SET_LOOP_DES( WDES%COMM_KIN, IMAG_GRIDS, 'time', IO )
#ifdef mem_profiling
    CALL GET_PSS( PSS, IO%IU0, IO%IU6, "SETTING GDES",0 )
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
 
    CALL INIT_SUPERCELL_FFT( SUPER, P, LMDIM,  WHF%WDES,  WGW, GDES_TAU, LATT_CUR, ENCUTGW, IO%IU6)
#ifdef mem_profiling
   CALL GET_PSS( PSS, IO%IU0, IO%IU6, "SUPERCELL_FFT DONE",0 )
#endif

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

   ! allocate memory for CHI_TAU as well
   IF ( LACFDT) THEN
      ALLOCATE(CHI_TAU(S2E%NUMBER_OF_NQ))
      DO NQ_COUNTER=1,SIZE(CHI_TAU)  ! loop over  q-points in the IRZ
         NQ=S2E%NQ(NQ_COUNTER)
         CALL SET_RESPONSE_KPOINT(CHI_TAU(NQ_COUNTER), WDES%VKPT(:,NQ), NQ)
         CALL ALLOCATE_RESPONSEFUN_DISTRI( CHI_TAU(NQ_COUNTER), &
              GDES_TAU%RES_NRPLWV_ROW, GDES_TAU%RES_NRPLWV_COL, WGW%LGAMMA, WGW%LGAMMA, 1)
         ! set inverse temperature 
         CHI_TAU(NQ_COUNTER)%BETA = IMAG_GRIDS%BETA        
      ENDDO
   ELSE
   ! CHI_TAU needs to be allocated for all k-points, since calculation of GW in supercell requires this
   ! TODO: gK one could save memory if k-point rotation is done inside GW routine
      ALLOCATE(CHI_TAU(WDES%NKPTS))
      DO NQ=1,SIZE(CHI_TAU)
         CALL SET_RESPONSE_KPOINT(CHI_TAU(NQ), WDES%VKPT(:,NQ), NQ)
         CALL ALLOCATE_RESPONSEFUN_DISTRI( CHI_TAU(NQ), &
              GDES_TAU%RES_NRPLWV_ROW, GDES_TAU%RES_NRPLWV_COL, WGW%LGAMMA, WGW%LGAMMA, 1)
         ! set inverse temperature 
         CHI_TAU(NQ)%BETA = IMAG_GRIDS%BETA        
      ENDDO
   ENDIF

   ! OEP requires to invert chi(w=0) at Gamma point, 
   ! this object is stored in CHI0
   IF ( LOEP ) THEN
      ALLOCATE(CHI0(1))
      CALL ALLOCATE_RESPONSEFUN_DISTRI( CHI0(1), &
         GDES%RES_NRPLWV_ROW, GDES%RES_NRPLWV_COL, WGW%LGAMMA, WGW%LGAMMA, 1)
      CHI0(1)%BETA = IMAG_GRIDS%BETA        
      CALL CLEAR_RESPONSE(CHI0(1))
   ENDIF

   CALL START_TIMING("GWLOOP")  ! super timer for e.g. one loop
   CALL START_TIMING("GREENS")  ! cummulative timer for individual steps

   IF (IO%IU6>=0) THEN
      WRITE(IO%IU6,'(A)') ' all response functions allocated'
   ENDIF
   CALL DUMP_ALLOCATE_TAG(IO%IU6, 'Response functions allocated')
#ifdef mem_profiling
   CALL GET_PSS( PSS, IO%IU0, IO%IU6, "RESPONSE ALLOCATION DONE",0 )
#endif

   !store Fermi energy in wave function
   WHF%EFERMI=W%EFERMI
   IF (LALLOCATED_WCORR) WCORR%EFERMI=W%EFERMI

   CALL QP_PADE_READER(W,IO%IU5,IO%IU0,IO%IU6,NBANDSGW,NKPTS_IRZ)

!   CALL FT_TEST_GF(W, IMAG_GRIDS, IO )
!====================================================================
! Determine CHI_QPT on imaginary frequency axis 
! and add the q->0 limit 
! optionally, store chi(w=0) at Gamma point in CHI0 at 
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
            LMDIM,CDIJ,CQIJ, 0, SV,T_INFO,P,IO%IU0,E%EXHF,&
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
    CALL POTENTIAL_AND_CHARGE(1, &
       W, WDES, E, INFO, &
       GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
       LATT_CUR, T_INFO, P, SYMM, MIX, &
       CHTOT, CHTOTL, DENCOR, CVTOT, CSTRF, IRDMAX, &
       LMDIM, CQIJ, CDIJ, CRHODE, N_MIX_PAW, RHOLM, RHOLM_LAST, CHDEN, SV, IO)
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
    WMEAN%EFERMI=W%EFERMI  
    WHF_IN_KS=W ; NULLIFY(WHF_IN_KS%CELTOT, WHF_IN_KS%FERTOT,  WHF_IN_KS%AUXTOT)
    CALL ALLOCW_NOPLANEWAVE(WDES,WHF_IN_KS) ; WHF_IN_KS%FERTOT=W%FERTOT ; WHF_IN_KS%CELTOT = W%CELTOT ; WHF_IN_KS%AUXTOT=W%AUXTOT

    CALL GREENS_MAT_DES_INIT( GDES_MAT, IMAG_GRIDS%T, WDES%NB_TOT)

    NULLIFY(SIGMAW_MAT, GU_MAT, GO_MAT)
#ifdef MatsubaraIntegration
    NULLIFY(GMATU_STORE, GMATO_STORE)
    CALL ALLOCATE_GREENS_MAT_GDEF(GDES_MAT, GMATU_STORE, NKPTS_IRZ, WDES%ISPIN,& 
    IMAG_GRIDS%T%NPOINTS)
    CALL ALLOCATE_GREENS_MAT_GDEF(GDES_MAT, GMATO_STORE, NKPTS_IRZ, WDES%ISPIN,& 
    IMAG_GRIDS%T%NPOINTS)
    GMATU_STORE=zero
    GMATO_STORE=zero
#endif
    IF (LGW .AND. ((.NOT. LG0W0 ) .AND.  (.NOT. LscQPGW))) THEN
       ! allocate space for self-energy in frequency domain, scaLAPACK distributed
       CALL ALLOCATE_GREENS_MAT(GDES_MAT, SIGMAW_MAT, NKPTS_IRZ, WDES%ISPIN, &
       IMAG_GRIDS%T%NPOINTS_IN_ROOT_GROUP)
       SIGMAW_MAT=(0._q,0._q)
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
    IF(IO%IU0>=0)&
       WRITE(IO%IU0,'(" Calculation of W done, Calculation of self-energy started for")')
    IF(IO%IU0>=0)WRITE(17,'(" Calculation of W done, Calculation of self-energy started for")')

    IGOAL_GW = IMAG_GRIDS%T%NPOINTS_IN_ROOT_GROUP*W%WDES%ISPIN*NKPTS_IRZ*NELMGW
    IPROGRESS_GW =0  

! TODO: dump warning for forcing NELM = 1 if  ALGO = G0W0R ; NELM>1 is set
    IF ( NELMGW > 1 .AND. LG0W0 ) NELMGW = 1 
scf: DO NELM=1,NELMGW
    IF (NELM>1) THEN
       CALL NO_SYMMETRY(.TRUE., LGAMMA, SYMM, W, WDES, LATT_CUR, LATT_INI, &
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
    SIGMA_COS=zero
    SIGMA_SIN=zero

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
    FRNL=0._q
    MP2SUM=(0._q,0._q) 
    MP2SUM2=(0._q,0._q)
    IU6_MEM=IO%IU6

    E_MP2 = 0 
  !----------------------------------------------------------------------- 
  ! start loop over time poins for calculation of SIGMA = G W    
  tau_points_gw: DO NTAU_ROOT=1, IMAG_GRIDS%T%NPOINTS_IN_ROOT_GROUP
  !----------------------------------------------------------------------- 
    CALL SETUP_TAU_INDICES(NTAU_ROOT,IMAG_GRIDS%T)
    DO ISP=1,W%WDES%ISPIN
    !------------------------------------------------------------------- 
    ! calculate the occupied and unoccupied Green function
    !------------------------------------------------------------------- 

       !allocate storage for self-energies
       ! TODO: most likely some allocations can be moved outside of the tau loop
       ALLOCATE(SIGMAU(NKPTS_IRZ), SIGMAO(NKPTS_IRZ ))
       CALL NULLIFY_G_RECIPROCAL( GDES_TAU, SIGMAO, NKPTS_IRZ)
       CALL NULLIFY_G_RECIPROCAL( GDES_TAU, SIGMAU, NKPTS_IRZ)
       IF (.NOT. LUSESUPER) THEN
          CALL ALLOCATE_G_RECIPROCAL( GDES_TAU, SIGMAO, NKPTS_IRZ)
          CALL ALLOCATE_G_RECIPROCAL( GDES_TAU, SIGMAU, NKPTS_IRZ)
       ENDIF

       !allocate storage for self-energies in Bloch domain
       ALLOCATE(SIGMAO_DIAG(W%WDES%NB_TOT))
       ALLOCATE(SIGMAU_DIAG(W%WDES%NB_TOT))
       SIGMAO_DIAG=zero
       SIGMAU_DIAG=zero

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

       !we have G in tau and W in nu, need to loop over q-points and get W_q(nu)
       !for all q or rather k from 1st BZ
       !the screened potential W is stored in the response array CHI_QPT
       !then CALCULATE_SIGMA_SUPER can be called to get everything 

       !first transform the W_q(nu) to W_q(tau)
       CALL START_TIMING("G")
       IU6_TMP=IU6_MEM
       DO NQ=1,WDES%NKPTS  ! loop over all difference vectors in the full BZ
       ! k-point outside IRZ rotate potential from IRZ to desired k-point
          IF (NQ>NKPTS_IRZ) THEN
#ifndef gammareal
             ! equivalent k-point in IRZ (somewhat akward, but correct)
             NQ_IRZ=KPOINTS_FULL_ORIG%NEQUIV(KPOINT_IN_FULL_GRID(KPOINTS_FULL%VKPT(:,NQ),KPOINTS_FULL_ORIG))

             !CHI should point to current q-point
             CHI=>CHI_QPT(NQ_IRZ)

             IF ( IMAG_GRIDS%T%LDO_POINT_LOCAL ) CALL CLEAR_RESPONSE(CHI_TAU(NQ)) 
             ! we take CHI_QPT from IrBZ and transform it to CHI_TAU in BZ
             ! set VKPT to k-point in IRZ
             CALL SET_RESPONSE_KPOINT(CHI_TAU(NQ), WDES%VKPT(:,NQ_IRZ), NQ_IRZ)

             ! inverse cosine transformation of screened Coulomb potential W for current tau-point
             ! W is stored in CHI_TAU
             ! transform W from nu to tau; backward transformation set by LFORWARD=.FALSE. 
             CALL TRANS_TIME_FREQUENCY(GDES, GDES_TAU, CHI, CHI_TAU(NQ), &
             IMAG_GRIDS, NTAU_ROOT, .FALSE. , IO%IU0 )

             ! note: depending on whether we work in time or frequency, GDES_TAU or GDES
             ! now transform to NQ
             CALL ROTATE_RES(CHI_TAU(NQ)%RESPONSEFUN(:,:,1), GDES_TAU, WGW%NGVECTOR(NQ), &
                  RES_KPTS_TRANS%CPHASE(:,NQ), RES_KPTS_TRANS%NINDPW(:,NQ),  &
                  RES_KPTS_TRANS%LINV(NQ), RES_KPTS_TRANS%LSHIFT(NQ))

             ! now set k-points to NQ
             CALL SET_RESPONSE_KPOINT(CHI_TAU(NQ), WDES%VKPT(:,NQ), NQ)
#endif
          ELSE
             !CHI should point to current q-point 
             CHI=>CHI_QPT(NQ)

             IF ( IMAG_GRIDS%T%LDO_POINT_LOCAL ) CALL CLEAR_RESPONSE(CHI_TAU(NQ)) !one q-point at a time
             CALL SET_RESPONSE_KPOINT(CHI_TAU(NQ), WDES%VKPT(:,NQ), NQ)

             ! inverse cosine transformation of screened Coulomb potential W for current tau-point
             ! W is stored in CHI_TAU, for a given q-point
             ! transform W from nu to tau for a given tau point; backward transformation set by LFORWARD=.FALSE. 
             CALL TRANS_TIME_FREQUENCY(GDES, GDES_TAU, CHI, CHI_TAU(NQ), &
             IMAG_GRIDS, NTAU_ROOT, .FALSE. , IO%IU0 )
#ifdef debug
         CALL DUMPX_FREQ(100, CHI_TAU(NQ), GDES_TAU, IMAG_GRIDS%T%COMM_BETWEEN_GROUPS, IMAG_GRIDS%T%COMM_IN_GROUP, &
              WDES%COMM_KIN%NODE_ME, IMAG_GRIDS%T%POINTS_LOCAL)
#endif 

          ENDIF
          IF (.NOT. LUSESUPER) THEN
          DO NK1=1,NKPTS_IRZ  ! loop restricted to IRZ
             IF (MOD(NK1-1,WDES%COMM_KINTER%NCPU)/=WDES%COMM_KINTER%NODE_ME-1) CYCLE

             ! Sigma is at NK1, W is at NQ, and GU at NK2=NK1-NQ and GO at NQ
             ! k2=q-k1
             NK2=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,NK1)-CHI_TAU(NQ)%VKPT(:),KPOINTS_FULL)

             I = DETERMINE_NTAU_GLOBAL(WDES%COMM_KIN%NODE_ME, IMAG_GRIDS%T%NPOINTSC,&
             IMAG_GRIDS)
             TAU_TMP = KPOINTS_ORIG%WTKPT(NK1)*IMAG_GRIDS%TAU_WEIGHT( I ) 

             ! compute self-energy for occupied orbital  SIGMAO = G( -tau) W(tau) = GO W
             CALL CALCULATE_SIGMA_TAU( LMDIM, WHF, WGW, GO(NK2), GU(NK1), SIGMAO(NK1), GDES_TAU, CHI_TAU(NQ),&
                  NQ, NK2, NK1, LATT_CUR, IU6_TMP, TAU_TMP, FRNL, CALCULATE_RPA_FORCES())

             CALL STOP_TIMING("G",IU6_TMP,"*sigmaocc")

             ! compute self-energy for unoccupied orbitals  SIGMAU = G( tau) W(tau) = GU W = -Sigma(+tau)
             CALL CALCULATE_SIGMA_TAU( LMDIM, WHF, WGW, GU(NK2), GO(NK1), SIGMAU(NK1), GDES_TAU, CHI_TAU(NQ),&
                  NQ, NK2, NK1, LATT_CUR, IU6_TMP, TAU_TMP, FRNL, CALCULATE_RPA_FORCES())
             CALL STOP_TIMING("G",IU6_TMP,"*sigmaun")

             IU6_TMP=-1

#ifdef debug
             WRITE(100+WDES%COMM_KIN%NODE_ME,*)CHI_TAU(NQ)%RESPONSEFUN
#endif 
          ENDDO ! k-points
          ENDIF
       ENDDO ! q-point
       !now we have CHI_TAU for all k-vectors and in tau domain
       CALL STOP_TIMING("GREENS",IO%IU6,"SIGMA")

       IF (LUSESUPER) THEN
          ! now calculate the self-energy as G(r,r',tau) W(r,r',tau)
          ! for a given tau point treated currently
          I = DETERMINE_NTAU_GLOBAL(WDES%COMM_KIN%NODE_ME, IMAG_GRIDS%T%NPOINTSC,&
          IMAG_GRIDS)
          TAU_TMP = IMAG_GRIDS%TAU_WEIGHT( I ) 
          CALL CALCULATE_SIGMA_SUPER( SUPER, WHF, WGW, GO, GU, GDES_TAU, &
               CHI_TAU, SIGMAO, S2E%NUMBER_OF_NQ, S2E%NQ, LMDIM, LATT_CUR, IU6_MEM, & 
               TAU_TMP, FRNL, CALCULATE_RPA_FORCES(),IO)
          CALL STOP_TIMING("GREENS",IO%IU6,"SIGMAOC")
          
          !calculate self-energy for negative tau (unoccupied states or G)
          CALL CALCULATE_SIGMA_SUPER( SUPER, WHF, WGW, GU, GO, GDES_TAU, &
               CHI_TAU, SIGMAU, S2E%NUMBER_OF_NQ, S2E%NQ, LMDIM, LATT_CUR, IU6_MEM, &
               TAU_TMP, FRNL, CALCULATE_RPA_FORCES(),IO)
          CALL STOP_TIMING("GREENS",IO%IU6,"SIGMAUN")
       ENDIF

       CALL START_TIMING("G")

       IU6_TMP=IU6_MEM
       ! now determine the matrix elements <i|Sigma(tau)|a>
       DO NK1=1,NKPTS_IRZ
          IF (MOD(NK1-1,WDES%COMM_KINTER%NCPU)/=WDES%COMM_KINTER%NODE_ME-1) CYCLE

          I = DETERMINE_NTAU_GLOBAL(WDES%COMM_KIN%NODE_ME, IMAG_GRIDS%T%NPOINTSC,&
          IMAG_GRIDS)
          TAU_TMP = KPOINTS_ORIG%WTKPT(NK1)*IMAG_GRIDS%TAU_WEIGHT( I ) 
          IF (L2ORDER .OR. LRPAFORCE) THEN
#ifdef RPAgamma
             ! SIGMAO und SIGMAU is in uncompressed storage mode, so need to bring GU and GO also
             ! to uncompressed mode
             IF ( ICHIREAL == 2 ) CALL UNCOMPRESS_G_RECIPROCAL ( GU(NK1), GDES_TAU, NK1)
             IF ( ICHIREAL == 2 ) CALL UNCOMPRESS_G_RECIPROCAL ( GO(NK1), GDES_TAU, NK1)
#endif
             IF (L2ORDER) THEN
                ! for second order (L2ORDER): G(tau) G( -tau) W(tau) is the MP2 energy
                CALL CONTRACT_SIGMA_G_TAU( WHF, GU(NK1), SIGMAO(NK1), GDES_TAU, NK1, &
                     IO%IU6, TAU_TMP, IMAG_GRIDS%T%LDO_POINT_LOCAL,E_MP2 )
             ENDIF
             CALL CONTRACT_SIGMA_G_TAU_DER( W%WDES, NK1, TAU_TMP,&
                  GU(NK1), SIGMAO(NK1), GDES_TAU, NONL_S, NONLR_S, LATT_CUR, FRNL,&
                  IMAG_GRIDS%T%LDO_POINT_LOCAL)
          ENDIF

          ! determine <i| G( -tau) W(tau) | a> in the orbital basis
          ! W(tau) symmetric (?)
          ! input : SIGMAO(NK1)%(g,g')
          ! output: SIGMAO_MAT(ab')
          CALL CALCULATE_SIGMA_ORBITAL(W, SIGMAO(NK1), SIGMAO_DIAG, SIGMAO_MAT, &
               GDES_MAT, GDES_TAU, IMAG_GRIDS, -1, NK1, ISP, .TRUE., MP2SUM, TAU_TMP)
#ifdef debug
          IF(NTAU_ROOT>1) WRITE(*,*)'orbital done  ', WDES%COMM_KIN%NODE_ME,NTAU_ROOT, NK1,NKPTS_IRZ
          CALL STOP_TIMING("G",IU6_TMP,"*s_orbital")
#endif
          IF (L2ORDER .OR. LRPAFORCE) THEN
             IF (L2ORDER) THEN
                ! for second order (L2ORDER): G(tau) G( -tau) W(tau) is the MP2 energy
                CALL CONTRACT_SIGMA_G_TAU( WHF, GO(NK1), SIGMAU(NK1), GDES_TAU,NK1, & 
                     IO%IU6, TAU_TMP, IMAG_GRIDS%T%LDO_POINT_LOCAL,E_MP2 )
             ENDIF
             CALL CONTRACT_SIGMA_G_TAU_DER( W%WDES, NK1, TAU_TMP, &
                  GO(NK1), SIGMAU(NK1), GDES_TAU, NONL_S, NONLR_S, LATT_CUR, FRNL, IMAG_GRIDS%T%LDO_POINT_LOCAL)
          ENDIF
          ! determine <i| G( tau) W(tau) | a> in the orbital basis
          CALL CALCULATE_SIGMA_ORBITAL(W, SIGMAU(NK1), SIGMAU_DIAG, SIGMAU_MAT, &
               GDES_MAT, GDES_TAU, IMAG_GRIDS, -1, NK1, ISP, .FALSE., MP2SUM2, TAU_TMP)

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
                CALL ADD_G_SIGMA( W, SIGMAO_MAT(:, 1, 1, 1), SIGMAU_MAT(:, 1, 1, 1), CORR_MAT(:, NK1, ISP, 1), & 
                  NTAU_ROOT, IMAG_GRIDS, NK1, ISP, GDES_MAT%DESC )

                IF (IO%NWRITE>3) CALL DUMP_GREENS_HAM_GDEF( "corr mat", &
                   CORR_MAT(:,NK1,ISP,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)

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
          IF ( LGW .AND. ((.NOT. LG0W0 ) .AND.  (.NOT. LscQPGW)) ) THEN
             ! transform the selfenergy at the current tau point (SIGMAO_MAT and SIGMAU_MAT)
             ! to imaginary frequency and store it in SIGMAW_MAT
             ! the complete self-energy is accumulated during the tau loop
             CALL FT_SIGMA_MAT( W, SIGMAO_MAT(:, 1, 1, 1), SIGMAU_MAT(:, 1, 1, 1), & 
                  SIGMAW_MAT(:, NK1, ISP,:), NTAU_ROOT, IMAG_GRIDS, GDES_MAT )
          ENDIF

          IU6_TMP=-1
          ! dump progress
          IPROGRESS_GW = IPROGRESS_GW + 1 
          CALL DUMP_PROGRESS( IPROGRESS_GW, IGOAL_GW, 'Sigma', IO )
!store SIGMA matrix in time domain
#ifdef MatsubaraIntegration
          ! store diagonal self-energy elements on time axis 
          GMATU_STORE(:,NK1,ISP,DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME,&
             NTAU_ROOT, IMAG_GRIDS )) = SIGMAU_MAT(:,NK1,ISP,1)
          GMATO_STORE(:,NK1,ISP,DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME,&
             NTAU_ROOT, IMAG_GRIDS )) = SIGMAO_MAT(:,NK1,ISP,1)
#endif
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

!#define verbose
#ifdef verbose
CALL DUMP_GREENS_HAM_GDEF( "first order change", CORR_MAT(:,NKPTS_IRZ,1,1), W%WDES%NB_TOTK(NKPTS_IRZ,1), GDES_MAT, IO%IU0)
#endif

  ! diagonal part of self-energy on imaginary time axis 
    CALLMPI( M_sum_g(IMAG_GRIDS%T%COMM_BETWEEN_GROUPS, SIGMAO_TAU, SIZE(SIGMAO_TAU)))
    CALLMPI( M_sum_g(IMAG_GRIDS%T%COMM_BETWEEN_GROUPS, SIGMAU_TAU, SIZE(SIGMAU_TAU)))
  ! compute G*Sigma in frequency domain 
    CALL GALITSKII_MIGDAL_DIAG( W, SIGMAO_TAU, SIGMAU_TAU, IMAG_GRIDS, IO)

#ifdef MatsubaraIntegration
    CALLMPI( M_sum_g(IMAG_GRIDS%T%COMM_BETWEEN_GROUPS, GMATO_STORE, SIZE(GMATO_STORE)))
    CALLMPI( M_sum_g(IMAG_GRIDS%T%COMM_BETWEEN_GROUPS, GMATU_STORE, SIZE(GMATU_STORE)))
    CALL FT_G_TEST( W, GDES_MAT, GMATU_STORE, GMATO_STORE, IMAG_GRIDS, IO)
    DEALLOCATE( GMATO_STORE, GMATU_STORE )
#endif

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
    IF (IO%IU0>=0.AND. (NELMGW==1 .OR. .NOT.LFINITE_TEMPERATURE)) THEN
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
          LMDIM, CQIJ, CDIJ, CRHODE, N_MIX_PAW, RHOLM, RHOLM_LAST, CHDEN, SV,&
          IO, OEP)
    ELSE
       CALL POTENTIAL_AND_CHARGE(NELM_HF, &
          W, WDES, E, INFO, &
          GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
          LATT_CUR, T_INFO, P, SYMM, MIX, &
          CHTOT, CHTOTL, DENCOR, CVTOT, CSTRF, IRDMAX, &
          LMDIM, CQIJ, CDIJ, CRHODE, N_MIX_PAW, RHOLM, RHOLM_LAST, CHDEN, SV, &
          IO, OEP) 
    ENDIF

  !----------------------------------------------------------------------- 
  ! Determine reference Greens function
  ! CHAM_MAT is determined and contains the HF Hamiltonian matrix
  ! update fermi energy and occupancies 
  ! also, calculate correlated density matrix G(t=0-) stored in CORR_MAT
  !----------------------------------------------------------------------- 
    E%ETRGHF = 0
    DTRGHF = 0
    CALL REF_GF_HFHAM_CORR_DENS_MAT( CHAM_MAT, CORR_MAT, SIGMAW_MAT ) 
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
    MP2SUM = 0      ! change of number of electrons
    TR_GAMMA = 0    ! Trace of density matrix 
    E%EBANDSTR = 0 
    E%EENTROPY = 0

    CALL ALLOCATE_GREENS_MAT_GDEF(GDES_MAT, GAMMA, 1, 1, 1)

    DO ISP=1,WDES%ISPIN
    DO NK1=1,NKPTS_IRZ

       ! calculate change in Trace of density matrix from correlation effects
       ! CORR_MAT contains (in orbital basis) either G_0 Sigma G_0 or G_0 Sigma G (if LGW0=.TRUE.)
       ! CORR_MAT contains G_0 + G_0 Sigma G if LFINITE_TEMPERATURE = .TRUE. and LGW0 = .TRUE.
       ! CORR_MAT contains G_0 Sigma G_0 if LFINITE_TEMPERATURE = .FALSE. and LG0W0 = .TRUE.
       ! integrated over omega, ie the change of the density matrix for tau->0-
       ! R stores the diagonal elements
       CALL DETERMINE_DIAGONALE_GDEF(WDES%NB_TOTK(NK1,ISP), CORR_MAT(:,NK1,ISP,1), R, GDES_MAT%DESC)
       CALLMPI( M_sum_g(GDES_MAT%COMM_INTRA, R, WDES%NB_TOTK(NK1, ISP) ))
       IF (.NOT.LscQPGW) MP2SUM=MP2SUM+SUM(R(1:WDES%NB_TOTK(NK1,ISP)))*KPOINTS_ORIG%WTKPT(NK1)*WDES%RSPIN

!#ifdef verbose
!       CALL DUMP_GREENS_HAM_GDEF( "HF exchange", CHAM_MAT(:,NK1,ISP,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
!#endif
!#define RPA_HF_TEST
#ifdef RPA_HF_TEST
          SIGMA_MAT(:,NK1,ISP,1) = 0
          CORR_MAT(:,NK1,ISP,1) = 0
          FRNL=0
#endif
       ! for GW0:  GAMMA is set to the HF density matrix in the current basis of eigenvectors (DFT in iteration 1)
       ! for G0W0: GAMMA is set to H_HF / abs(epsilon_a - epsilon_i)
       !           this is only usefull for first iteration
       ! for QPGW: GAMMA is simply set to HF matrix
       GAMMA(:,1,1,1)=CHAM_MAT(:,NK1,ISP,1)
       IF (LGW .AND. ((.NOT. LG0W0 ) .AND.  (.NOT. LscQPGW))) THEN
          IF (NELM_HF==NELMHF) THEN
          ! set GAMMA to Hartree-Fock density matrix, and add mean field Green's function to GO_MAT
             CALL HF_DENSITY_MATRIX( WMEAN, NK1, ISP, GAMMA(:,1,1,1), GDES_MAT, IMAG_GRIDS, EFERMI_MEAN(ISP), &
                  GU_MAT, GO_MAT, LUNOCCUPIED=LHOLEGF)
          ELSE
             ! use density matrix from solution of dyson equation
             ! this is done below, so just clean array here
             CALL HF_DENSITY_MATRIX( WMEAN, NK1, ISP, GAMMA(:,1,1,1), GDES_MAT,&
                  IMAG_GRIDS, EFERMI_MEAN(ISP), LUNOCCUPIED = LHOLEGF )
          ENDIF
       ELSE IF (.NOT. LscQPGW) THEN
!#ifdef verbose
!       CALL DUMP_GREENS_HAM_GDEF( " B (sigma matrix)", SIGMA_MAT(:,NK1,ISP,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
!       CALL DUMP_GREENS_HAM_GDEF( " A (gamma)", GAMMA(:,1,1,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
!#endif
          ! in the RPA forces paper this is gamma^(1)_HF 
          ! gam_ij = F_ij * (fi ei - fj ej)/(ei-ej) 
          ! note that the variation of the G_KS = ( (iw + mu) S - H ) ^-1  =: A^-1
          ! follows from the simple formulat d A^-1 = - A^-1 (dA) A^-1  and is 
          ! d G_KS = - G_KS [ (iw+mu) dS - dH ] G_KS 
          ! this is why only CELTOT is passed and not CELTOT-EFERMI
          IF (ASSOCIATED(SIGMA_MAT)) THEN
             CALL DIVIDE_BY_EIGENVALUEDIFFxEIGENVAL( W%WDES%NB_TOTK(NK1,ISP),&
                GAMMA(:,1,1,1), SIGMA_MAT(:,NK1,ISP,1), W%CELTOT(:,NK1,ISP), &
                W%FERTOT(:,NK1,ISP), GDES_MAT%DESC, &
                LFINITE_TEMPERATURE, FERDER(:,NK1,ISP) )
          ENDIF 
!#ifdef verbose
!       CALL DUMP_GREENS_HAM_GDEF( "B + F_mn* A", SIGMA_MAT(:,NK1,ISP,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
!#endif
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
             !
             ! include entropy term for RPA forces, where the entropy is given by 
             ! 
             ! E_ent = sum_i { f_i ln(f_i) + (1-f_i) ln(1-f_i) }
             !  
             ! 
             ! d(E_ent)/d R = sum_i{  T ln{ f_i/(1-f_i) } * {d f_i / d x_i } * 1/T ( d e_i/d R )}
             !                        \_________________/   \____________________/
             !                             -(e_i-mu)             f'_i                 
             ! 
             ! the last factor can be deduced from Perturbation theory: e_i + de_i
             !
             ! <n|H+dH|n+dn> = <n|(e_n+de_n)(S + dS)|n+dn>
             !
             ! now use the fact that <n+dn|S+dS|n+dn> = 1 is required 
             ! which implies <n|S|n> + <dn|S|n> + <n|S|dn> + <n|dS|n> = 1 
             !               \_____/
             !                  =1
             !                             \
             !  => <dn|S|n> + <n|S|dn> = 0  |   => <dn|S|n> = 0 = <n|S|dn>
             !                \______/      |
             !                <dn|S|n>     / 
             ! 
             ! collecting terms proportional to d^0 yields the generalised
             ! eigenvalue problem H | n >  = e_n S |n >. But collecting terms
             ! proportional do d yields an expression for the perturbation de_n of
             ! the eigenvalue e_n
             ! 
             !  de_j = d_ij * <i| dH - e_j dS |j>
             ! 
             ! So the entropic correction  is given by 
             ! 
             ! [rho^1]_ij = -d_ij* ( e_j - mu )*f'_j 
             ! 
             ! [gam^1]_ij = -d_ij* ( e_j - mu )*f'_j * e_j
             !
             IF ( LRPAFORCE ) THEN
                !
                ! note that f'_j is multiplied below by DIVIDE_BY_EIGENVALUEDIFF
                !
                RR = W%EFERMI(ISP)-W%CELTOT(:,NK1,ISP)
                CALL  ADD_TO_DIAGONALE_REAL( WDES%NB_TOTK(NK1,ISP), GAMMA(:,1,1,1), RR, GDES_MAT%DESC)
                !
                ! for gam^1 this must be done manually here
                !
                IF ( ASSOCIATED( SIGMA_MAT ) )  THEN
                   RR = RR*( W%CELTOT(:,NK1,ISP)*FERDER(:,NK1,ISP) )
                   CALL  ADD_TO_DIAGONALE_REAL( WDES%NB_TOTK(NK1,ISP), SIGMA_MAT(:,NK1,ISP,1), RR, GDES_MAT%DESC)
                ENDIF
             ENDIF
          ! gamma = G_DFT + G_DFT (H_mean- H_DFT) G_DFT 
          ! in first iteration gamma = H_mean (f_a - f_i) / (epsilon_a - epsilon_i)
          ! in the paper this is rho^(1)_HF 
             CALL DIVIDE_BY_EIGENVALUEDIFF( W%WDES%NB_TOTK(NK1,ISP), &
                GAMMA(:,1,1,1), W%CELTOT(:,NK1,ISP)-W%EFERMI(ISP), W%FERTOT(:,NK1,ISP), &
                GDES_MAT%DESC, LFINITE_TEMPERATURE, FERDER(:,NK1,ISP) )

             ! add unperturbed density matrix  1/ pi int_0^infty   G^W_0 d w
             CALL  ADD_TO_DIAGONALE_REAL( WDES%NB_TOTK(NK1,ISP), GAMMA(:,1,1,1), W%FERTOT(:,NK1,ISP), GDES_MAT%DESC)
          ENDIF
          ! gK this is twisted CORR_MAT is set here to exchange (GAMMA) and correlation (CORR_MAT) contribution
          ! gK LOEP .AND. LFINITE_TEMPERATURE .AND. LRPAFORCE
          !
          ! note that G0 is subtracted in STORE_RPA_DENSITY_MATRIX, thats why it is added here
          !
          IF ( LFINITE_TEMPERATURE .AND. LRPAFORCE ) THEN
             CORR_MAT(:,NK1,ISP,1)=CORR_MAT(:,NK1,ISP,1) + GAMMA(:,1,1,1)
          ENDIF
       ENDIF

#ifdef verbose
       CALL DUMP_GREENS_HAM_GDEF( "exchange density matrix", GAMMA(:,1,1,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
       CALL DUMP_GREENS_HAM_GDEF( "correlation density matrix", CORR_MAT(:,NK1,ISP,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
       IF (ASSOCIATED(GO_MAT)) &
       CALL DUMP_GREENS_HAM_GDEF( "green fun GO_MAT", GO_MAT(:,NK1,ISP,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
#endif
       ! add mean field (Hartree-Fock) and correlation contributions of the density matrix
       IF ( LOEP .OR. .NOT. LFINITE_TEMPERATURE ) THEN
          ! gK: not using finite temperature; add exhange (GAMMA) and correlation (COTT_MAT) now
          GAMMA(:,1,1,1)=GAMMA(:,1,1,1)+CORR_MAT(:,NK1,ISP,1)
!CALL DUMP_GREENS_HAM_GDEF( "r.h.s.", GAMMA(:,NK1,ISP,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
       ELSE
          ! gK now CORR_MAT is copied over to GAMMA
          !in finite temperature formalism GAMMA  = lim tau->0- of interacting G
          ! G = G_0 + G_0 SIGMA G solved for G -> GAMMA = G(tau=0-)  
          GAMMA(:,1,1,1)=CORR_MAT(:,NK1,ISP,1)
       ENDIF
       ! overwrite the diagonale elements by the diagonale of the self-energy (Kotani-Schilfgaarde)
       IF (LscQPGW) THEN
          CALL SET_DIAGONALE_REAL( WDES%NB_TOTK(NK1,ISP), GAMMA(:,1,1,1), REAL(SIGMA(:,NK1,ISP),q), GDES_MAT%DESC)
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
!#ifdef verbose
!          CALL DUMP_GREENS_HAM_GDEF( "exchange + correlation matrix", GAMMA(:,1,1,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
!          CALL DETERMINE_DIAGONALE_GDEF(WDES%NB_TOTK(NK1,ISP), GAMMA(:,1,1,1), R, GDES_MAT%DESC)
!          CALLMPI( M_sum_g(GDES_MAT%COMM_INTRA, R, WDES%NB_TOTK(NK1, ISP) ))
!#endif

          CALL PDSYEV_ZHEEVD_DESC( GAMMA(:,1,1,1), RR(1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT%DESC, GDES_MAT%COMM_INTRA)
          CALLMPI( M_bcast_d_from(GDES_MAT%COMM_INTER , RR(1), SIZE(RR), 1))
!#ifdef verbose
!          IF (IO%IU0>=0) WRITE(IO%IU0,'("old   ", 10F12.6)') REAL(W%CELTOT(1:MIN(10,W%WDES%NB_TOTK(NK1,ISP)),NK1,ISP))
!          IF (IO%IU0>=0) WRITE(IO%IU0,'("occup ", 10F12.6)') RR(1:MIN(10,W%WDES%NB_TOTK(NK1,ISP)))
!#endif
          ! store eigenvalues in W%CELTOT
          W%CELTOT(1:W%WDES%NB_TOTK(NK1,ISP),NK1,ISP)=RR(1:W%WDES%NB_TOTK(NK1,ISP))
       ELSE
        !
        ! clean some blocks in Natural orbial matrix
        ! 
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
          ! undo previous sign change
          RR=-RR   
          GAMMA(:,1,1,1)=-GAMMA(:,1,1,1)
#ifdef verbose
          IF (IO%IU0>=0) WRITE(IO%IU0,'("eigvalues 1-10 ", 10F12.6)') RR(1:MIN(10,W%WDES%NB_TOTK(NK1,ISP)))
          IF (IO%IU0>=0) WRITE(IO%IU0,'("eigvalues 11-20", 10F12.6)') RR(MIN(11,W%WDES%NB_TOTK(NK1,ISP)):MIN(20,W%WDES%NB_TOTK(NK1,ISP)))
#endif
        ENDIF
        ! store diagonal of density matrix in FERTOT for later use
        IF (NATURALO==0 .OR. NATURALO==10 .OR. NATURALO==3 .OR. NATURALO==13) THEN
           IF ( LOEP ) THEN
              W%FERTOT(1:W%WDES%NB_TOTK(NK1,ISP),NK1,ISP)=RR(1:W%WDES%NB_TOTK(NK1,ISP))
           ELSE
              W%FERTOT(1:W%WDES%NB_TOTK(NK1,ISP),NK1,ISP)=MAX(RR(1:W%WDES%NB_TOTK(NK1,ISP)),0.0_q)
              IF (NATURALO==3 .OR. NATURALO == 13) W%FERTOT(1:N2,NK1,ISP)=1
           ENDIF
        ENDIF
        !
        ! Luttinger-Ward quantities
        ! 
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
!#ifdef verbose
!       CALL DUMP_GREENS_HAM_GDEF( "rotation matrix", GAMMA(:,1,1,1),W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
!#endif
       IF ( .NOT. LOEP ) THEN
         ! restore original mean field HF matrix
         IF (LGW .AND. ((.NOT. LG0W0 ) .AND.  (.NOT. LscQPGW))) CALL HF_MATRIX( WMEAN, NK1, ISP, CHAM_MAT(:,NK1,ISP,1), GDES_MAT )
         ! calculate change of HF energy in rotated basis 
         CALL HF_SINGLES(  W, NK1, ISP, GAMMA(:,1,1,1), CHAM_MAT(:,NK1,ISP,1), GDES_MAT, E%EBANDSTR)
       ENDIF
    ENDIF

!#ifdef verbose
!    CALL DUMP_GREENS_HAM_GDEF( "Gamma", GAMMA(:,1,1,1), W%WDES%NB_TOTK(NK1,ISP), GDES_MAT, IO%IU0)
!#endif
    ! store the rotation matrix for later use
    CHAM_MAT(:,NK1,ISP,1)=GAMMA(:,1,1,1)

    ENDDO  !k-points
    ENDDO  !spin
    CALL DEALLOCATE_GREENS_MAT_GDEF(GDES_MAT, GAMMA)

    ! close unit for rotation matrix    
    IF ( LFCIINPUT ) CLOSE( 1000 )   
 
    IF (IO%IU0>=0.AND.(NELMGW==1 .OR. .NOT. LFINITE_TEMPERATURE)) THEN
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
       TOTEN = TOTEGM
       E%EENTROPY = 0
    ENDIF 

    IF (.NOT. LRPAFORCE) THEN
       IF ( .NOT. LG0W0 .AND. LFINITE_TEMPERATURE ) THEN
          CALL DUMP_DIAGRAMMATIC_ENERGIES( " Energies after update of density matrix ", &
               E, TOTEN, INFO, IO, -DOT_PRODUCT(CMPLX(EFERMI_MEAN,0._q),TR_GAMMA))
       ELSE
          CALL DUMP_DIAGRAMMATIC_ENERGIES_RPA( " Energies after update of density matrix ", &
               E, TOTEN, INFO, IO)
       ENDIF
       IF ( LHOLEGF ) THEN
          IF (IO%IU6>=0) WRITE(IO%IU6,'("  Trace of density matrix (holes)       ",2F20.10)') SUM(TR_GAMMA(:)) 
       ELSE     
          IF (IO%IU6>=0) WRITE(IO%IU6,'("  Trace of density matrix (electrons)   ",2F20.10)') SUM(TR_GAMMA(:)) 
       ENDIF
    ENDIF
    IF ( .NOT. LFINITE_TEMPERATURE .AND. IO%IU6>=0) &
         WRITE(IO%IU6,'("  correlated contrib. to density matrix ",2F20.10)') MP2SUM
    IF (IO%IU6>=0) WRITE(IO%IU6,'(/)')

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
       ! form natural orbitals and store them in W or WOEP
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

       ! update wavefunction and charge with current OEP potential 
       CALL UPDATE_W_OEP(W, SV, CDIJ)

    ENDIF
    CALL STOP_TIMING("GREENS",IO%IU6,"ROTATE")
    ! the re-orthogonalization is not required, since all transformations
    ! are unitary, documenting in the call should yield exactly identical results
!    CALL ORTHCH(WDES,W, WDES%LOVERL, LMDIM, CQIJ)
  ENDDO iterhf

  CALL DEALLOCATE_GREENS_MAT_GDEF(GDES_MAT, CHAM_MAT)
  CALL DEALLOCATE_GREENS_MAT_GDEF(GDES_MAT, CORR_MAT)
  IF (ASSOCIATED(SIGMA_MAT)) CALL DEALLOCATE_GREENS_MAT_GDEF(GDES_MAT, SIGMA_MAT)

  IF (( LGW .AND. ((.NOT. LG0W0 ) .AND.  (.NOT. LscQPGW))) .AND.  LFINITE_TEMPERATURE  ) THEN
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
  ! finite-temperature case the breaking condition is when LW energy converges 
  IF ( LFINITE_TEMPERATURE ) THEN
     IF ( INFO%EDIFF > ABS( FLWTOT_OLD - TOTEN))  THEN
        EXIT scf
     ENDIF
     FLWTOT_OLD = TOTEN
  ENDIF

  ENDDO scf

  IF ( LGW .AND. ((.NOT. LG0W0 ) .AND.  (.NOT. LscQPGW)) ) THEN
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
    DO NQ_COUNTER=1,SIZE(CHI_TAU)
       CALL DEALLOCATE_RESPONSEFUN(CHI_TAU(NQ_COUNTER))
    ENDDO
    DEALLOCATE(CHI_QPT, CHI_TAU)

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
       CALL DEALLOCATE_KPOINTS_TRANS(KPTS_TRANS)
    ENDIF

    CALL DEALLOCATE_IRZ_MAP(S2E)

    !destroy imaginary grid handle
    CALL DEALLOCATE_IMAG_GRID_HANDLE(IMAG_GRIDS)
    DEALLOCATE(WGW%VKPT, WGW%WTKPT)
    DEALLOCATE(WGW, GRIDWGW)
    CALL DEALLOCATE_SUPERCELL(SUPER)
    CALL DEALLOCATE_OEP_HANDLE(OEP)

    PROFILING_STOP('calculate_xi_real')

  CONTAINS

!***********************************************************************
!
! Determine reference Greens function (exchange gf )
! update fermi energy and occupancies 
! also, calculate correlated density matrix G(t=0-)
!
! SIGMAW_MAT is set to 0 on exit
!
!***********************************************************************
    SUBROUTINE  REF_GF_HFHAM_CORR_DENS_MAT( CHAM_MAT, CORR_MAT, SIGMAW_MAT ) 
       GDEF, POINTER, INTENT(INOUT) :: CHAM_MAT(:,:,:,:)
       GDEF, POINTER, INTENT(INOUT) :: CORR_MAT(:,:,:,:)
       COMPLEX(q),POINTER, INTENT(INOUT) :: SIGMAW_MAT(:,:,:,:)
       INTEGER :: ISP, NK1, N
       REAL(q) F
       COMPLEX(q):: MP2SUM, MP2SUM2
       REAL(q) HF_DIAG(WDES%NB_TOT)
       LOGICAL :: LGWMODE 
!tester, evaluates Galitskii-Migdal correlation energy in time domain
!works, currently only for std and ncl version
!#define GMInTau
#ifdef GMInTau
       COMPLEX(q),ALLOCATABLE  :: GS(:)! auxillary 
       COMPLEX(q),ALLOCATABLE  :: CDIAG(:)! auxillary 
       INTEGER  :: NOMEGA_GLOBAL, J
       REAL(q) :: E_GMTR, E_GMTAU
#endif

    ! determine <i|H|a>= <i| T + V_Hartree + V_exchange | a> store in CHAM_MAT(i,a,nk,ispin)
    ! for LGW0 CHAM_MAT will be set to unitary matrix diagonalizing <i|H|a> U_aj = epsilon_j U_jj
    ! and WMEAN%CELTOT is set to corresponding eigenvalues
    ! for G0W0  WMEAN%CELTOT = <i | T + V_Hartree + V_exchange | i>

    CALL EDDIAG_SIMPLE(HAMILTONIAN, &
         GRID, LATT_CUR, NONLR_S, NONL_S, W, WMEAN, WHF_IN_KS, WDES, &
         (LGW.AND.((.NOT. LG0W0 ).AND.(.NOT. LscQPGW))).OR.(LFINITE_TEMPERATURE.AND.LGW0), SYMM, &
         LMDIM, CDIJ, CQIJ, SV, T_INFO, P, NKPTS_IRZ, CHAM_MAT, GDES_MAT, E, IO%IU6 )
    CALL STOP_TIMING("GREENS",IO%IU6,"EDDIAG")
    IF ((LG0W0.OR.LscQPGW) .OR. NELM==1) THEN
       IF (IO%IU6>=0) WRITE(IO%IU6,'(1X,A)', ADVANCE="NO") "QP shifts evaluated in KS or natural orbital/ Bruckner basis"
       !WHF_IN_KS stores diagonal elements of HF Hamiltonian matrix in KS orbitals
       !spline fit of SIGMA_SIN_SIN to cos grids and replace 2/3 of the SIGMA_SIN
       !CALL WRITE_SELF_ENERGY(WHF_IN_KS, SIGMA_COS, SIGMA_SIN, NKPTS_IRZ, IMAG_GRIDS, NBANDSGW, IO%IU6, NELM, 'KS/NA basis')
       CALL WRITE_SELF_ENERGY(WMEAN, SIGMA_COS, SIGMA_SIN, NKPTS_IRZ, IMAG_GRIDS, NBANDSGW, IO%IU6)
       !pL Analytic continuation of Sigma(iw) to real axis and calculate the QP poles
       CALL QP_PADE_FIT_FROM_TAU(W,WHF_IN_KS,IMAG_GRIDS,SIGMAU_TAU,SIGMAO_TAU,NKPTS_IRZ,NBANDSGW,IO, .TRUE.)
    ENDIF
   
    TOTEN=E%EBANDSTR+E%DENC+E%XCENC+E%TEWEN+E%PSCENC+E%PAWPS+E%PAWAE+INFO%EALLAT+E%EXHF
    ! Galitskii-Migdal formula represent internal energy (see Fetter-Walecka Eq. 23.15)
    ! so exclude entropy term here 
    E%ETRGHF = TOTEN 
    ! compute change in bandstructure energy between iterations  
    DTRGHF = E%EBANDSTR

    CALL DUMP_DIAGRAMMATIC_ENERGIES_RPA("Energies using frozen KS orbitals", E, TOTEN, INFO, IO ) 

    ! following holds true if ALGO = GWR | GW0R
    !
    ! the electronic HF energy is calculated as 
    !
    ! E_HF = \sum_i f_i e_i - E_dc, where E_dc = 1/2( V_h - V_x )
    !
    ! -E_dc = E%DENC + E%EXHF + ( E%XCENC + E%PAWPS + E%PAWAE ) 
    ! 
    ! these terms essentially are the U_0+Tr( G_ref * Sigma_x ) term in PRA 73, 012511 Eq. (B5)
    IF ( LGW .AND. ((.NOT. LG0W0 ) .AND. (.NOT. LscQPGW)) ) THEN

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
       CALL DUMP_DIAGRAMMATIC_ENERGIES_RPA( "Energies after diagonalization of HF Hamiltonian (single shot)", E, TOTEN,INFO, IO) 
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
            IMAG_GRIDS, SIGMA_COS, SIGMA_SIN, GDES_MAT, IO)
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
    IF ( ( LGW .AND. ((.NOT. LG0W0 ) .AND. (.NOT. LscQPGW))) .AND. NELM_HF==1 ) THEN
       ! if diag_cham is set, the self-energy is always stored in the NO basis
       ! if diag_cham is not set, the self-energy is now transformed from the NO basis to the mean field basis
       !   the second option is much more efficient 
#ifndef diag_cham
       ! now rotate correlated selfenergy into Hartree-Fock basis
       CALL ROTATE_GREEN_SIGMA(WDES, CHAM_MAT, SIGMAW_MAT, IMAG_GRIDS, GDES_MAT)

       ! Pade fit of self-energy in Hartree-Fock basis
       CALL PADE_FIT_GW0( W, WMEAN, IMAG_GRIDS, GDES_MAT, SIGMAW_MAT, NKPTS_IRZ, NBANDSGW, IO )
       ! and rotate KS orbitals to diagonalize Hartree-Fock Hamiltonian
       CALL ROTATE_KS_ORBITALS(W, WDES, CHAM_MAT, GDES_MAT)
       CALL STOP_TIMING("GREENS",IO%IU6,"ROTSIGMA")
#endif

       IF (LFERMIGW) THEN
          CALL ITERATE_EFERMI_G_FROM_SIGMA( WMEAN, NKPTS_IRZ, INFO%NELECT, INFO%NUP_DOWN, EFERMI_MEAN, EFERMI_SLOPE, & 
               SIGMAW_MAT, IMAG_GRIDS, GDES_MAT, CHAM_MAT, IO )
! mix the Fermi-level with previous one, and store it as the reference value for the next iteration
! mixing might accelerate convergence
           ! TODO: check if this is correct!
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
       
       ! this computes the GM energy in imaginary time 
#ifdef GMInTau
       E_GMTAU=0
       IF (ASSOCIATED(GO_MAT)) THEN
          ALLOCATE(GS(SIZE(SIGMAW_MAT,1)))
          GS=0

          DO ISP=1,WDES%ISPIN
          DO NK1=1,NKPTS_IRZ
             ! additional storage is needed to store Green's functions in time
             ! since GO_MAT and GU_MAT are overwritten
             NULLIFY(SIGMAU_MAT)
             NULLIFY(SIGMAO_MAT)
             CALL ALLOCATE_GREENS_MAT_GDEF(GDES_MAT, SIGMAU_MAT, 1, 1,IMAG_GRIDS%T%NPOINTS_IN_ROOT_GROUP)
             CALL ALLOCATE_GREENS_MAT_GDEF(GDES_MAT, SIGMAO_MAT, 1, 1,IMAG_GRIDS%T%NPOINTS_IN_ROOT_GROUP)
             ! size of diagonal
             N=  W%WDES%NB_TOTK(NK1, ISP)
             ALLOCATE(CDIAG(N))
             E_GMTR=0 

             DO NTAU_ROOT=1, IMAG_GRIDS%T%NPOINTS_IN_ROOT_GROUP
                CALL SETUP_TAU_INDICES( NTAU_ROOT, IMAG_GRIDS%T ) 

                CALL INV_FT_SIGMA_MAT(W, W%WDES%NB_TOT, SIGMAO_MAT(:, 1, 1, NTAU_ROOT), SIGMAU_MAT(:, 1, 1, NTAU_ROOT),&
                     SIGMAW_MAT(:, NK1, ISP,:), NTAU_ROOT, IMAG_GRIDS, GDES_MAT )

                ! contract Galtiskii-Migdal energy 
                CALL PZGEMM( 'N', 'N', N, N, N, (1._q,0.0_q),  &
                     SIGMAU_MAT(1,1,1,NTAU_ROOT), 1, 1, GDES_MAT%DESC, &
                     GO_MAT(1,NK1,ISP,NTAU_ROOT), 1, 1, GDES_MAT%DESC, &
                     (0.0_q,0.0_q),  &
                     GS(1), 1, 1, GDES_MAT%DESC )
                CALL PZGEMM( 'N', 'N', N, N, N, (1._q,0.0_q),  &
                     SIGMAO_MAT(1,1,1,NTAU_ROOT), 1, 1, GDES_MAT%DESC, &
                     GU_MAT(1,NK1,ISP,NTAU_ROOT), 1, 1, GDES_MAT%DESC, &
                     (1.0_q,0.0_q),  &
                     GS(1), 1, 1, GDES_MAT%DESC )

                !obtain diagonal 
                CDIAG = 0 
                CALL DETERMINE_DIAGONALE(N, GS(:), CDIAG(:), GDES_MAT%DESC)
                IF( .NOT.IMAG_GRIDS%T%LDO_POINT_LOCAL ) CDIAG = 0 
                !this is the global tau point of the node
                NOMEGA_GLOBAL= DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, NTAU_ROOT,&
                IMAG_GRIDS)
                !IF (NOMEGA_GLOBAL> IMAG_GRIDS%NOMEGA) CYCLE
                ! build trace
                DO J = 1 , N
!                   IF( ABS( AIMAG( CDIAG(J) )) > 1.E-6_q) WRITE(*,300)AIMAG( CDIAG(J) ) 
                   E_GMTR = E_GMTR - REAL(CDIAG(J),q)*&
                      IMAG_GRIDS%TAU_WEIGHT(NOMEGA_GLOBAL)
                ENDDO
             ENDDO
300    FORMAT(' WARNING Galiskii-Migdal energy calculated in tau-domain is complex: ',F20.10)
             ! accumulate contributions from trace 
             CALLMPI( M_sum_z(GDES_MAT%COMM_INTER, E_GMTR, 1))
             E_GMTAU = E_GMTAU + REAL(E_GMTR,q)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN
             ! get grid of auxillary arrays
             CALL DEALLOCATE_GREENS_MAT_GDEF(GDES_MAT, SIGMAU_MAT)
             CALL DEALLOCATE_GREENS_MAT_GDEF(GDES_MAT, SIGMAO_MAT)
             DEALLOCATE(CDIAG)

          ENDDO
          ENDDO
          ! auxillary arrays are not needed 
          DEALLOCATE(GS)
       ENDIF
#endif   
       ! solve Dyson equation: G = 1/( iw +mu - T -V_ext -V_h -V_x -Sigma_c(iw))
       ! SIGMAW_MAT contains G in HF basis if diag_cham is not set
       ! otherwise SIGMAW_MAT contains G in the NO basis
       CALL CALCULATE_G_FROM_SIGMA( WMEAN, NKPTS_IRZ, EFERMI_MEAN, SIGMAW_MAT, &
            IMAG_GRIDS, GDES_MAT, CHAM_MAT, E, IO )
 
#ifdef GMInTau
       IF( IO%IU0>=0) WRITE(*,'("GM in tau and omega,diff:",2F20.10,E13.6)' ) E_GMTAU , E%ECGWGM, E_GMTAU - E%ECGWGM
       ! at zero-temperature GM is hard to obtain in frequency, 
       ! so we calculate the term in time, however this deteriorates scf convergence 
       ! because we determine G_(i-1) * Sigma_i, where i is the scf loop iteration
       ! however, this is a good test to check if the code works at T=0 as well 
       E%ECGWGM = E_GMTAU
#endif

     !----------------------------------------------------------------------- 
     ! check if self-energy on real axis can be determined  
     !----------------------------------------------------------------------- 
       ! do this for self-consistent G 
!       CALL PERFORM_INTEGRATION_TEST( W, SIGMAW_MAT, IMAG_GRIDS, GDES_MAT, &
!            NKPTS_IRZ, EFERMI_MEAN, E, IO )  

       ! if for terms depending of O( 1/omega ) at T=0
       IF ( .NOT. LFINITE_TEMPERATURE ) THEN
          CALL SUBTRACT_GHF_FROM_G( WMEAN, NKPTS_IRZ, EFERMI_MEAN, SIGMAW_MAT, &
               IMAG_GRIDS, GDES_MAT, CHAM_MAT )
       ENDIF 

       CALL STOP_TIMING("GREENS",IO%IU6,"SIGMA2G")

       ! write the diagonal part of SIGMAW_MAT
!       CALL WRITE_SELFENERGY_W(W, NKPTS_IRZ, SIGMAW_MAT, IMAG_GRIDS, GDES_MAT, IO%IU0)
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


       ! calculate contribution to density matrix from correlation part (overwrite CORR_MAT)
       CALL GAMMA_FROM_G_W( WMEAN, NKPTS_IRZ, SIGMAW_MAT, CORR_MAT, &
            IMAG_GRIDS, GDES_MAT, CHAM_MAT, LUNOCCUPIED=LHOLEGF  )

       ! at this point SIGMAW_MAT and CORR_MAT contain the interacting G-G_mean
       ! and the corresponding density matrix
       CALLMPI( M_sum_g(GDES_MAT%COMM_INTER, CORR_MAT, SIZE(CORR_MAT)))
       ! G(w) is not neded anymore
       SIGMAW_MAT=0

!#define verbose
!#ifdef verbose
!CALL DUMP_GREENS_HAM_GDEF( "after GAMMA_FROM_G_W", CORR_MAT(:,NKPTS_IRZ,1,1), W%WDES%NB_TOTK(NKPTS_IRZ,1), GDES_MAT, IO%IU0)
!#endif
       CALL STOP_TIMING("GREENS",IO%IU6,"GAMMA")
    ENDIF
    END SUBROUTINE  REF_GF_HFHAM_CORR_DENS_MAT
!*********************************************************************
!
! subroutine to test frequency to time FT
!
!*********************************************************************
   SUBROUTINE TEST_FT
     INTEGER :: NO, NT
     REAL(q) :: INVFTCOS(IMAG_GRIDS%NOMEGA,IMAG_GRIDS%NOMEGA)
     REAL(q) :: INVFTSIN(IMAG_GRIDS%NOMEGA,IMAG_GRIDS%NOMEGA)
     COMPLEX(q) :: FT_PLUS (IMAG_GRIDS%NOMEGA)
     COMPLEX(q) :: FT_MINUS(IMAG_GRIDS%NOMEGA)
     COMPLEX(q) :: FW
     REAL(q)   :: EPSILON=0.4, BETA

     INVFTCOS=IMAG_GRIDS%TO_BOS_RE
     INVFTSIN=IMAG_GRIDS%TO_BOS_IM
     CALL INVERT_REAL_MATRIX(INVFTCOS,-1)
     CALL INVERT_REAL_MATRIX(INVFTSIN,-1)

     ! loop over all time points
     ! first calculate results for negative times
     ! F(t+) = 1/2( F(w) + F(w)*)*FTCOS^-1  +
     !         1/2( F(w)*- F(w) )*FTSIN^-1
     !
     ! F(t-) =-1/2( F(w) + F(w)*)*FTCOS^-1  +
     !         1/2( F(w)*- F(w) )*FTSIN^-1
     FT_MINUS=0
     FT_PLUS =0

     DO NT=1,IMAG_GRIDS%NOMEGA
        DO NO=1,IMAG_GRIDS%NOMEGA
           
           ! ( iw - epsilon)^-1
           FW=-1/(CMPLX(0._q,IMAG_GRIDS%FER_RE(NO),q)-EPSILON)

           ! F(w) + F(w)*)*FTCOS^-1 terms
           FT_PLUS(NT)=FT_PLUS(NT)+ & 
                (FW+CONJG(FW))*INVFTCOS(NT,NO)*0.25_q
           FT_MINUS(NT)=FT_MINUS(NT)- & 
                (FW+CONJG(FW))*INVFTCOS(NT,NO)*0.25_q
           
           ! (F(w)*- F(w))*FTSIN^-1 terms
           FT_PLUS(NT)=FT_PLUS(NT)- & 
                (FW-CONJG(FW))*INVFTSIN(NT,NO)*(0.0_q,0.25_q)
           FT_MINUS(NT)=FT_MINUS(NT)- & 
                (FW-CONJG(FW))*INVFTSIN(NT,NO)*(0.0_q,0.25_q)
        ENDDO
     ENDDO

     BETA=IMAG_GRIDS%TAU(IMAG_GRIDS%NOMEGA)*2
     DO NT=IMAG_GRIDS%NOMEGA,1,-1
        IF (EPSILON>0) THEN
           WRITE(*,'(4F14.7)') -IMAG_GRIDS%TAU(NT), REAL(FT_MINUS(NT),q), 0.0_q, & 
                EXP(EPSILON*IMAG_GRIDS%TAU(NT))/(1+EXP(BETA*EPSILON))
        ELSE
           WRITE(*,'(4F14.7)') -IMAG_GRIDS%TAU(NT), REAL(FT_MINUS(NT),q), EXP(-ABS(IMAG_GRIDS%TAU(NT)*EPSILON)), & 
                EXP(EPSILON*IMAG_GRIDS%TAU(NT))/(1+EXP(BETA*EPSILON))
        ENDIF
     ENDDO

     DO NT=1,IMAG_GRIDS%NOMEGA
        IF (EPSILON>0) THEN
           WRITE(*,'(4F14.7)')  IMAG_GRIDS%TAU(NT), REAL(FT_PLUS(NT),q), EXP(-ABS(IMAG_GRIDS%TAU(NT)*EPSILON)), &
                EXP(EPSILON*(BETA-IMAG_GRIDS%TAU(NT)))/(1+EXP(BETA*EPSILON))
        ELSE
           WRITE(*,'(4F14.7)')  IMAG_GRIDS%TAU(NT), REAL(FT_PLUS(NT),q),  0.0_q, &
                EXP(EPSILON*(BETA-IMAG_GRIDS%TAU(NT)))/(1+EXP(BETA*EPSILON))
        ENDIF
     ENDDO
     CALL vtutor%stopCode()
     
   END SUBROUTINE TEST_FT

 
   SUBROUTINE UPDATE_W_OEP( W, SV, CDIJ )
      USE setexm
      TYPE( wavespin ) W
      RGRID       SV(DIMREAL(WDES%GRID%MPLWV),WDES%NCDIJ)
      OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      !local
      REAL(q) DESUM1, RMS,EF
      INTEGER :: ICOUEV, NSIM 
   
      ICOUEV = 0 
      DESUM1 = 0.0_q
      RMS = 0.0_q
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

      ! wavefunction grid (GRID) -> charge density grid (GRIDC)
      CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
           GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
           LATT_CUR, P, SYMM, T_INFO, &
           CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

   END SUBROUTINE UPDATE_W_OEP
  END SUBROUTINE CALCULATE_XI_REAL



#endif

!here we should put following internal routines of both monsters:
!RECALCULATEW
!REF_GF_HFHAM_CORR_DENS_MAT

END MODULE chi_super
