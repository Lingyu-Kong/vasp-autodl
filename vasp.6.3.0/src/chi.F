!#define debug
#include "symbol.inc"

!*********************************************************************
!
! This module implements the top level GW routines
! mostly written by gK with some contributions from Maxim Shishkin
! notably the analytical Kramers Kronig transformations using finite
! element basis in the frequency domain
!
! TODO: block inversion in XI_INVERT
!       block inversion in XI_LOCAL_FIELD
! 
!*********************************************************************

MODULE xi
  USE chi_base
  USE local_field
  USE bse
  USE wpot
  USE crpa 
#ifdef PROFILING
  USE profiling
#endif
  USE acfdt
  USE chi_glb 
  IMPLICIT NONE

CONTAINS


   
  SUBROUTINE RESPONSE_SET_ENCUT(ENMAX)
    USE tutor, ONLY: vtutor, isAlert, RPAENCUTGWSOFT, argument
    REAL(q) :: ENMAX

    ! nothing in the INCAR file, set ENCUTGW to default (2/3 ENCUT)
    IF (ENCUTGW==-2) THEN
       IF (LCHI .OR. LMP2 .OR. LMP2KPAR .OR. LMP2NO .OR. LFCIDUMP .OR. LRPAX .OR. LCCSD .OR. LBRACKETST) THEN
          ! default for MP2, GW etc.
          ENCUTGW=ENMAX/3*2
       ELSE
          ! default for cases (typically OEP/LHF) methods
          ! typically use all PW components in a sphere that fits into the FFT-box
          ! sphere size is 3/4*2 if default precision is used
          IF (EXXOEP==1) THEN
             ENCUTGW=ENMAX*(3._q/4._q*2)**2
             IF (ENCUTGWSOFT<0) ENCUTGWSOFT=ENCUTGW
          ELSE IF (EXXOEP>=2) THEN
             ENCUTGW=ENMAX
             IF (ENCUTGWSOFT<0) ENCUTGWSOFT=ENCUTGW
          ENDIF
       ENDIF
       ! default for ENCUTGWSOFT
       ! ENCUTGWSOFT now consistenly defaults to ENCUTGW*0.8
       ! this increases compatabilty between different calculations
       ! also RPA force calculations converge slower than they need to when ENCUTGWSOFT = ENCUTGW
       IF ((LACFDT .OR. LMP2 .OR. LMP2KPAR  .OR. LMP2NO .OR. LFCIDUMP .OR. LRPAX &
      & .OR. LCCSD .OR. LBRACKETST) .AND. ENCUTGWSOFT_OLD ==-2 ) THEN
          ENCUTGWSOFT_OLD=ENCUTGW*0.8
       ELSE IF ( ENCUTGWSOFT_OLD<0 ) THEN
          ENCUTGWSOFT_OLD=ENCUTGW
       ENDIF

       IF (ENCUTGWSOFT < 0 ) THEN
          ENCUTGWSOFT=ENCUTGW*0.8
       ENDIF
       ! ok to be entirely save
       IF (ENCUTGWSOFT > ENCUTGW) THEN
         ENCUTGWSOFT = ENCUTGW
       ENDIF
    ENDIF
    IF ( ABS( ENCUTGWSOFT - ENCUTGWSOFT_OLD ) > 0.01_q) CALL vtutor%write(isAlert, RPAENCUTGWSOFT, argument(rval = [ENCUTGWSOFT_OLD]))
  END SUBROUTINE RESPONSE_SET_ENCUT

!**********************************************************************
!
! write the response function related variables to a file
!
!**********************************************************************

  SUBROUTINE WRITE_RESPONSE(IU6)
    USE fock
    IMPLICIT NONE
    INTEGER IU6
    
    IF (IU6>=0 .AND. (LCHI.OR. LTIME_EVOLUTION .OR.LMP2.OR.LMP2KPAR.OR.LMP2NO .OR. LFCIDUMP  .OR. LRPAX .OR.  LCCSD .OR. LBRACKETST)) THEN
       WRITE(IU6,10 )  ALGO_FROM_GW(), LFINITE_TEMPERATURE, &
            LFXC, LADDER, LRPAFORCE, LHARTREE,IBSE,KPOINT_BSE, &
            LTRIPLET,LRSRPA,LTCTC,LTCTE,LTETE,LFXCEPS,LFXHEG,  &
            NATURALO, LHOLEGF, L2ORDER, LDMP1, LMP2LT,LSMP2LT, LGWLF, &
            ENCUTGW, ENCUTGWSOFT, ENCUTLF, LESF_SPLINES, LMAXMP2, &
            SCISSOR, NOMEGA, NOMEGAR, &
            NBANDSGW, NBANDSGWLOW,NBANDSO,NBANDSV,NELMGW,NELMHF, DIMGW, IEPSILON, ANTIRES,OMEGAMAX, OMEGAMIN, OMEGATL, OMEGAGRID, & 
            SHIFT,LSELFENERGY,LSPECTRAL,LSPECTRALGW,LSINGLES,LFERMIGW,ODDONLYGW,EVENONLYGW, NKREDLFX, NKREDLFY, NKREDLFZ, TELESCOPE, &
            NTAUPAR,NOMEGAPAR,DAMP_NEWTON
       IF (LAMBDA/=1.0_q) THEN
          WRITE(IU6,"('   LAMBDA  =',F6.3,'    coupling constant (scaling of Coulomb kernel)'/)") LAMBDA
       ENDIF
   ENDIF

10  FORMAT(' Response functions by sum over occupied states:'  / &
         '   ALGO                = ',A6,  '   selected many-body perturbation algorithm'/&
         '   LFINITE_TEMPERATURE =',L6,  '    use finite temperature formalism'/&
         '   LFXC                =',L6,  '    include DFT xc kernel in TD-HF and W'/&
         '   LADDER              =',L6,  '    particle hole ladder diagrams for BSE/TD-HF and W (Nano-quanta)'/&
         '   LRPAFORCE           =',L6,  '    calculate RPA forces'/&
         '   LHARTREE            =',L6,  '    include Hartree terms (bubbles-RPA) in BSE/TD-HF'/&
         '   IBSE                =',I6,  '    BSE modus: 0 exact, 1 time-evolution (store matrix), 10 time-evolution (implicit)'/&
         '   KPOINT              =',I6,3I4,'    k-point index at which BSE equation is solved (and G at which response is eval.)'/&
         '   LTRIPLET            =',L6,  '    triplet instead of singlet excitation (BSE only)'/&
         '   LRSRPA              =',L6,  '    range separated RPA (total - short) in ACFDT'/&
         '   LTCTC               =',L6,  '    use test-charge test-charge for self-energy (with ALGO=QPBSE)'/&
         '   LTCTE               =',L6,  '    use test-charge test-electron dielectric function'/&
         '   LTETE               =',L6,  '    use test-electron test-electron dielectric function'/&
         '   LFXCEPS             =',L6,  '    determine f_xc=epsilon-1/Xi Sharma and Gross kernel'/&
         '   LFXHEG              =',L6,  '    model exchange using f_x(G,k_F) from free electron gas'/&
         '   NATURALO            =',I6,  '    0 natural orbitals, 1 only virtual, 2 only occupied-virtual rotation'/&
         '   LHOLEGF             =',L6,  '    use holes in self-consistent GW '/&
         '   L2ORDER             =',L6,  '    2nd order terms only'/&
         '   LDMP1               =',L6,  '    first order term in ACFDT routine is calculated'/&
         '   LMP2LT              =',L6,  '    use LT MP2'/&
         '   LSMP2LT             =',L6,  '    use stochastic LT MP2'/&
         '   LGWLF               =',L6,  '    use W potentials for exchange like diagrams'/&
         '   ENCUTGW             =',F6.1,'    cutoff for response function in eV'/ &
         '   ENCUTGWSOFT         =',F6.1,'    soft cutoff for Coulomb kernel in GW response function in eV'/ &
         '   ENCUTLF             =',F6.1,'    cutoff for local field effects'/ &
         '   LESF_SPLINES        =',L6,  '    tricubic splines for electronic structure factor S(q+G)'/&
         '   LMAXMP2             =',I6,  '    maximum L for one center terms'/ &
         '   SCISSOR             =',F6.1,'    scissor correction'/ &
         '   NOMEGA              =',I6,  '    number of frequencies'/ & 
         '   NOMEGAR             =',I6,  '    number of frequencies along real axis'/ & 
         '   NBANDSGW            =',I6,  '    number of bands for which selfenergy shift is calculated'/ & 
         '   NBANDSGWLOW         =',I6,  '    lowest band included in GW (to exclude core correlation)'/ & 
         '   NBANDSO             =',I6,  '    number of bands for electron-hole treatment (occupied)'/ & 
         '   NBANDSV             =',I6,  '    number of bands for electron-hole treatment (virtual)'/ & 
         '   NELM                =',I6,  '    number of iterations in solving QP/BSE equation'/ & 
         '   NELMHF              =',I6,  '    number of iterations in the inner HF iteration'/ & 
         '   DIM                 =',I6,  '    dimensionality of system (0=0D molecules, 3=3D)'/ & 
         '   IEPSILON            =',I6,  '    cartesion direction for which dielectric function is calculated (4=all)'/ & 
         '   ANTIRES             =',I6,  '    antiresonant part (0) no TDA (1) w=0 exact (2) accurate'/ &
         '   OMEGAMAX            =',F6.1,'    maximum frequency'/ & 
         '   OMEGAMIN            =',F6.1,'    minimum frequency'/ & 
         '   OMEGATL             =',F6.1,'    maximum frequency of tail'/ & 
         '   OMEGAGRID           =',I6,  '    grid type (0 default)'/ & 
         '   CSHIFT              =',F6.1,'    complex shift used in evaluation of response functions'/&
         '   LSELFENERGY         =',L6,  '    calculate selfenergy instead of QP shifts'/&
         '   LSPECTRAL           =',L6,  '    use spectral functions                   '/&
         '   LSPECTRALGW         =',L6,  '    use spectral functions to calculate int dz G(w-z) W(z) (more accurate)'/&
         '   LSINGLES            =',L6,  '    calculate the singles contribution to the correlation energy'/&
         '   LFERMIGW            =',L6,  '    update Fermi-energy in the Greens function'//&
         '   ODDONLYGW           =',L6,  '    avoid all even points for polarizability '/&
         '   EVENONLYGW          =',L6,  '    avoid all odd  points for polarizability '/&
         '   NKREDLFX            =',I6,  '    NKREDLFY  =',I4,'   NKREDLFZ  =',I4,//&
         '   TELESCOPE           =',I6,  '    use telescope k-point sampling'/&
         '   NTAUPAR             =',I6,  '    number of TAU groups used in imaginary time calculation'/&
         '   NOMEGAPAR           =',I6,  '    number of OMEGA groups used in imaginary frequency calculation'/&
         '   DAMP_NEWTON         =',F6.4,  '    damping for Newton step for position update when RPA forces are used'/&
         )

    IF (IU6>=0.AND.((LRHFCALC.AND.LRSCOR).OR.LRSRPA)) THEN
       WRITE(IU6,11 ) HFSCREEN_ORIG
    ENDIF

11  FORMAT( &
       & '   HFSCREEN=',F6.3,'    HFSCREEN for range separated RPA'/)

  END SUBROUTINE WRITE_RESPONSE

  SUBROUTINE XML_WRITE_RESPONSE
    USE vaspxml
    IMPLICIT NONE
    ! local
    INTEGER IDUM, N
    REAL(q) RDUM
    COMPLEX(q) CDUM
    LOGICAL LDUM
    CHARACTER (1) :: CHARAC

    N=1
    CALL XML_TAG("separator","response functions")
    CALL XML_INCAR('ALGO','A',IDUM,RDUM,CDUM,LDUM,ALGO_FROM_GW(),N)
    CALL XML_INCAR('LFINITE_TEMPERATURE','L',IDUM,RDUM,CDUM,LFINITE_TEMPERATURE,CHARAC,N)
    CALL XML_INCAR('LADDER','L',IDUM,RDUM,CDUM,LADDER,CHARAC,N)
    CALL XML_INCAR('LRPAFORCE','L',IDUM,RDUM,CDUM,LRPAFORCE,CHARAC,N)
    CALL XML_INCAR('LFXC','L',IDUM,RDUM,CDUM,LFXC,CHARAC,N)
    CALL XML_INCAR('LHARTREE','L',IDUM,RDUM,CDUM,LHARTREE,CHARAC,N)
    CALL XML_INCAR('IBSE','I',IBSE,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR_V('KPOINT','I',KPOINT_BSE,RDUM,CDUM,LDUM,CHARAC,4)
    CALL XML_INCAR('LTCTC','L',IDUM,RDUM,CDUM,LTCTC,CHARAC,N)
    CALL XML_INCAR('LTCTE','L',IDUM,RDUM,CDUM,LTCTE,CHARAC,N)
    CALL XML_INCAR('LTETE','L',IDUM,RDUM,CDUM,LTETE,CHARAC,N)
    CALL XML_INCAR('LTRIPLET','L',IDUM,RDUM,CDUM,LTRIPLET,CHARAC,N)
    CALL XML_INCAR('LFXCEPS','L',IDUM,RDUM,CDUM,LFXCEPS,CHARAC,N)
    CALL XML_INCAR('LFXHEG','L',IDUM,RDUM,CDUM,LFXHEG,CHARAC,N)
    CALL XML_INCAR('NATURALO','I',NATURALO,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('LHOLEGF','L',IDUM,RDUM,CDUM,LHOLEGF,CHARAC,N)
    CALL XML_INCAR('L2ORDER','L',IDUM,RDUM,CDUM,L2ORDER,CHARAC,N)
    CALL XML_INCAR('LDMP1','L',IDUM,RDUM,CDUM,LDMP1,CHARAC,N)
    CALL XML_INCAR('LMP2LT','L',IDUM,RDUM,CDUM,LMP2LT,CHARAC,N)
    CALL XML_INCAR('LSMP2LT','L',IDUM,RDUM,CDUM,LSMP2LT,CHARAC,N)
    CALL XML_INCAR('LGWLF','L',IDUM,RDUM,CDUM,LGWLF,CHARAC,N)
    CALL XML_INCAR('ENCUTGW','F',IDUM,ENCUTGW,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('ENCUTGWSOFT','F',IDUM,ENCUTGWSOFT,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('ENCUTLF','F',IDUM,ENCUTLF,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('LESF_SPLINES','L',IDUM,RDUM,CDUM,LESF_SPLINES,CHARAC,N)
    CALL XML_INCAR('LMAXMP2','I',LMAXMP2,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('SCISSOR','F',IDUM,SCISSOR,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('NOMEGA','I',NOMEGA,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('NOMEGAR','I',NOMEGAR,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('NBANDSGW','I',NBANDSGW,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('NBANDSO','I',NBANDSO,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('NBANDSV','I',NBANDSV,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('NELMGW','I',NELMGW,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('NELMHF','I',NELMHF,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('DIM','I',DIMGW,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('IESPILON','I',IEPSILON,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('ANTIRES','I',ANTIRES,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('OMEGAMAX','F',IDUM,OMEGAMAX,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('OMEGAMIN','F',IDUM,OMEGAMIN,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('OMEGATL','F',IDUM,OMEGATL,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('OMEGAGRID','I',OMEGAGRID,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('CSHIFT','F',IDUM,SHIFT,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('LSELFENERGY','L',IDUM,RDUM,CDUM,LSELFENERGY,CHARAC,N)
    CALL XML_INCAR('LSPECTRAL','L',IDUM,RDUM,CDUM,LSPECTRAL  ,CHARAC,N)
    CALL XML_INCAR('LSPECTRALGW','L',IDUM,RDUM,CDUM,LSPECTRALGW  ,CHARAC,N)
    CALL XML_INCAR('LSINGLES','L',IDUM,RDUM,CDUM,LSINGLES, CHARAC,N)
    CALL XML_INCAR('LFERMIGW','L',IDUM,RDUM,CDUM,LFERMIGW, CHARAC,N)
    CALL XML_INCAR('ODDONLYGW','L',IDUM,RDUM,CDUM,ODDONLYGW,CHARAC,N)
    CALL XML_INCAR('EVENONLYGW','L',IDUM,RDUM,CDUM,EVENONLYGW,CHARAC,N)
    CALL XML_INCAR('NKREDLFX','I',NKREDLFX,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('NKREDLFY','I',NKREDLFY,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('NKREDLFZ','I',NKREDLFZ,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('MAXMEM','I',MAXMEM,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('TELESCOPE','I',TELESCOPE,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('NTAUPAR','I',NTAUPAR,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('NOMEGAPAR','I',NOMEGAPAR,RDUM,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('DAMP_NEWTON','F',IDUM,DAMP_NEWTON,CDUM,LDUM,CHARAC,N)
    CALL XML_INCAR('LAMBDA','F',IDUM,LAMBDA,CDUM,LDUM,CHARAC,N) 

    CALL XML_CLOSE_TAG

  END SUBROUTINE XML_WRITE_RESPONSE

!*********************************************************************
!
! main module
! for GW calculations
! it allows for G0W0, GW0, GW and QP-scGW QP-scGW0 calculations
! aditionally it is the main scheduler for the BSE calculations
! in the module bse.F
! local field effects can be included in the DFT using 
! various approximation
! ) in W only (i.e. TC-TC)
! ) in W and G (i.e. TC-TE)
! ) using a BSE derived f_xc(r,r') (see local_field.F)
! 
! The frequency dependency is usually evaluated using a fast
! spectral method, allowing for very fast GW calculations
! for bulk materials
!
!> @details @ref openmp :
!> All call to chi_base::clean_responsefunction_cache and
!> chi_base::clean_responsefunction_int are placed within parallel
!> regions.
!
!*********************************************************************

  SUBROUTINE CALCULATE_XI( &
          KINEDEN, HAMILTONIAN, P, WDES, NONLR_S, NONL_S, W, LATT_CUR, LATT_INI, &
          GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C,E, & 
          CHTOT, CHTOTL, DENCOR, CVTOT, CSTRF, IRDMAX, &
          T_INFO, DYN, INFO, IO, KPOINTS, SYMM, MIX, & 
          LMDIM, CQIJ, CDIJ, CRHODE, N_MIX_PAW, RHOLM, RHOLM_LAST, CHDEN, SV, & 
          EFERMI, NEDOS, DOS, DOSI )
    USE hamil
    USE meta
    USE base
    USE lattice
    USE pseudo
    USE lattice
    USE nonl_high
    USE msymmetry
    USE mpimy
    USE mkpoints
    USE constant
    USE poscar
    USE wave
    USE us
    USE pot
    USE pawm
    USE wave_high
    USE kpoints_change
    USE full_kpoints
    USE mlr_optic
    USE dfast
    USE ini
    USE fileio
    USE wave_cacher
    USE setexm
    USE sym_grad
    USE choleski
    USE hamil_struct_def
    USE david
    USE subrot
    USE mlwf
    USE minimax_struct, ONLY : imag_grid_handle
    USE minimax, ONLY : SET_IMAG_GRID_HANDLE
    USE pead, ONLY : PEAD_RESETUP_WDES,LPEAD_SYM_RED
    USE tutor, ONLY: vtutor, isError, isAlert, FermiFiniteTemp, NUPDOWN, NOmega32, &
        OmegaGridOmegaTL, LSPECTRALinefficient, LPOTOK
#ifdef VASP_HDF5
    USE vhdf5
#endif
    IMPLICIT NONE
! structures
    TYPE (tau_handle)  KINEDEN
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
    REAL(q)     EFERMI
    INTEGER     NEDOS
    REAL(q)     DOS(NEDOS,WDES%ISPIN),DOSI(NEDOS,WDES%ISPIN)
! local
    TYPE (responsefunction) CHI, CHIR, WPOT, TBSE, TBSEA, TVXC, CHI0
    TYPE (screened_2e_handle) S2E
    REAL(q) :: XCSIF(3,3)
    TYPE (energy)      E
    REAL(q), ALLOCATABLE ::  OMEGA(:)      ! real or imaginary part of  frequencies at which polarizability is calculated
    REAL(q), ALLOCATABLE :: OMEGAWEIGHT(:) ! weights for Gauss integration
    COMPLEX(q), ALLOCATABLE :: COMEGA(:)   ! complex frequencies
    TYPE(imag_grid_handle) :: IMAG_GRIDS   ! grid handle for minimax 
    LOGICAL :: LCHIR                       ! calculate resonant part of response function
    LOGICAL :: LCHIREALLOCATE              ! reallocate chi in each electronic step
    INTEGER :: NELM                        ! number of electronic steps
    INTEGER :: IRDMAA
    INTEGER :: NSTRIP, I, N1, N1_LAST, K1, K2, NSTRIP1_ACT, ISP, NOPER, NOPER_, NQ, NQ_COUNTER, NK
    INTEGER :: NOMEGA_CHI, NOMEGA_WPOT, NOMEGA_INDEX, NOMEGA_INDEXI
    REAL(q) :: NFLOAT, NFLOAT_, FSG0, RMST
    TYPE (wavedes), POINTER :: WGW
    TYPE (grid_3d), POINTER :: GRIDWGW
    TYPE (wavedes1) :: WGWQ
    INTEGER :: NCPU
    COMPLEX(qs), ALLOCATABLE :: SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:,:,:)
    LOGICAL, EXTERNAL  :: CALCULATE_RESPONSE_FUNCTIONS
    REAL(q)    PAR(1,1,1,1,WDES%NCDIJ),DOSPAR(1,1,1,WDES%NCDIJ)
    REAL(q), ALLOCATABLE :: FERTOT_INITIAL(:,:,:)
! tables for Kramers Kronig transformation
    COMPLEX(q),POINTER :: TABLE(:,:)            ! table for Kramers Kronig of repsonse
    COMPLEX(q),POINTER :: TABLE_RES(:,:)        ! table for Kramers Kronig of repsonse resonant part
    COMPLEX(q),POINTER :: TABLE_POT_PLUS(:,:)   ! table for Kramers Kronig of potential post. shift
    COMPLEX(q),POINTER :: TABLE_POT_MIN(:,:)    ! table for Kramers Kronig of potential neg. shift
! ACFDT related handle
    TYPE (correlation), POINTER :: COR
! one center terms
    TYPE (one_center_handle), POINTER :: H
!  original/ or QP eigenvalues
    TYPE (wavespin)     :: WACC1, WACC2, W_W
    TYPE (wavespin)     :: W_CORR
! variables related to OEP method
    REAL(q) :: DLM_EXX(N_MIX_PAW,WDES%NCDIJ)
    REAL(q) :: RMS, DESUM1
    REAL(q), ALLOCATABLE :: ERRORS(:)
    INTEGER :: ICOUEV, NSIM
! cache structure to store wavefunctions for non diagonal HF
    TYPE (wave_cache), POINTER ::  WCACHE
    LOGICAL :: LWFROMFILE
    LOGICAL :: SAVE_CACHE_MEMORY=.TRUE.
    LOGICAL :: LGAMMA
    TYPE (skpoints_trans) :: KPOINTS_TRANS_WPOT
    INTEGER :: ISTAT

!$  INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
#ifdef VASP_HDF5
    INTEGER :: IH5ERR
#endif

    PROFILING_START('calculate_xi')

    ICOUEV = 0 
    RMS = 0.0_q
    DESUM1 = 0.0_q
    IF (.NOT. CALCULATE_RESPONSE_FUNCTIONS()) THEN
       PROFILING_STOP('calculate_xi')
       RETURN
    ENDIF
! finite temperature NU_DOWN not implemented yet
    IF (INFO%NUP_DOWN>0 .AND. LFINITE_TEMPERATURE ) THEN
       CALL vtutor%error("T>0 GW algorithms for NUP_DOWN > 0 not implmeneted, sorry")
    ENDIF

! autoset available memory, incase MAXMEM if not set by user in INCAR 
    IF ( LAUTOSET_MAXMEM ) THEN
       CALL AUTOSET_AVAILABLE_MEMORY( MAXMEM, IO )
    ENDIF 

    IF( IO%IU0>=0 .AND. IO%NWRITE > 2 ) WRITE(*,*)'EFERMI=',EFERMI
    W%EFERMI(:)=EFERMI
!pull breaks here      
!finite temperature RPA and GW requires fermi smearing
    IF ( LFINITE_TEMPERATURE .AND. .NOT. ( KPOINTS%ISMEAR == -1 .OR.  KPOINTS%ISMEAR == 0 )  ) THEN
        CALL vtutor%write(isError, FermiFiniteTemp)
    ENDIF
    CALL START_TIMING("GWLOOP")

    IF (W%WDES%NKPTS == 1 .AND. W%WDES%VKPT(1,1)**2+W%WDES%VKPT(2,1)**2+W%WDES%VKPT(3,1)**2<1E-10_q) THEN
       LGAMMA = .TRUE.
    ELSE
       LGAMMA = .FALSE.
    ENDIF

    CALL SET_NBANDSGW( W)
!
! scGW can either use mixing for Hartree potential
! or a damped MD on the orbitals; both work well, although I prefer
! the damped MD
    IF (MIX%IMIX==0 .AND. LscQPGW ) THEN
       IF (IO%IU0>=0) WRITE(IO%IU0,*) 'WARNING: call to DENSTA to obtain consistent occupancies and Fermi-level'
       IF (INFO%NUP_DOWN>=0) THEN
          CALL vtutor%write(isAlert, NUPDOWN)
       INFO%NUP_DOWN=-1
       ENDIF
            
       !keep fermi energy as it is for finite temperature calculations 
       IF ( .NOT. LFINITE_TEMPERATURE ) THEN     
          CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
          INFO%NUP_DOWN,  E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
          NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)

          IF (IO%IU6>=0) &
          WRITE(IO%IU6,"(/,' The Fermi energy was updated, please check that it is located mid-gap',/ &
                       & ' values below the HOMO (VB) or above the LUMO (CB) will cause erroneous energies',/ &
                       & ' E-fermi : ', F8.4,/)" ) EFERMI
       ENDIF  
       ! update charge
       CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
            GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
            LATT_CUR, P, SYMM, T_INFO, &
            CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

    ENDIF
    IF( IO%IU0>=0 .AND. LDMP1) WRITE(*,*)'Linear term in ACFDT is being calculated'

    !no update for finite temperature calculations   
    IF ( .NOT. (LCRPA .OR. LFINITE_TEMPERATURE ) ) THEN
       CALL MEAN_CBM_VBM( W, EFERMI)

       IF (IO%IU6>=0) THEN
          WRITE(IO%IU6,"(/,' The Fermi energy was updated, please check that it is located mid-gap',/ &
          & ' values below the HOMO (VB) or above the LUMO (CB) will cause erroneous energies',/ &
          & ' E-fermi : ', F8.4,/)" ) EFERMI
          WRITE(IO%IU0,"(/,' The Fermi energy was updated, please check that it is located mid-gap',/ &
          & ' values below the HOMO (VB) or above the LUMO (CB) will cause erroneous energies',/ &
          & ' E-fermi : ', F8.4,/)" ) EFERMI
       ENDIF
    ENDIF

!
! slot in calculation of Fermi- wavevecor by fitting the
! exact exchange energy using |rho(G)|^^2 F_local_field(G) 
    IF (LFXHEG .AND. LACFDT) THEN
       ! determine optimized Fermi-vector for this system
       CALL  DETERMINE_HEX( &
            P,NONLR_S,NONL_S,W,LATT_CUR, &
            T_INFO,IO,SYMM, &
            LMDIM,CQIJ)
    ENDIF

!
! sync orbitals and eigenvalues (just to be sure everything is ok)
    CALL KPAR_SYNC_ALL(WDES,W)
    IF (SYMM%ISYM>0) &
      CALL CLEANUP_CELEN(W) ! set degenerated or near degenerated eigenvalues to a single value

! precalculate HF eigenvalues (not required for OEP and BSE calculations)

    IF (.NOT. LOEP .AND. (.NOT. LBSE .OR. LQPBSE)) THEN
       CALL SET_EIGENVALUE_HARTREE_KINETIC( &
         HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR, &
         T_INFO,INFO,IO,KPOINTS,GRID,GRID_SOFT, &
         GRIDC,GRIDUS,C_TO_US,SOFT_TO_C,SYMM, &
         CHTOT,DENCOR,CVTOT,CSTRF, &
         CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM, &
         CHDEN,SV,LMDIM,IRDMAX)
    ENDIF

    IF (SYMM%ISYM>0 .AND. ASSOCIATED(CELTOT_X)) THEN
       CALL CLEANUP_CELEN_HELPER(W, NBANDSGW, CELTOT_X )
       CALL CLEANUP_CELEN_HELPER(W, NBANDSGW, CELTOT_HARTREE_KINETIC )
    ENDIF

    ! TODO gK: I am not sure, whether this should be commented out
    ! scale the exchange by the coupling constant
    !IF (ASSOCIATED(CELTOT_X)) THEN
    !   CELTOT_X=CELTOT_X* LAMBDA
    !ENDIF
    ! calculate (1-lambda) (T+V_H+V_ion+ V_xc^DFT) + lambda (T+V_H+V_ion)
    !IF (LAMBDA/=1.0_q .AND. ASSOCIATED(CELTOT_HARTREE_KINETIC)) THEN
    !   CELTOT_HARTREE_KINETIC=W%CELTOT*(1.0_q-LAMBDA)+CELTOT_HARTREE_KINETIC*LAMBDA
    !ENDIF

    ! nullify required on some compilers
    NULLIFY(CHI%RESPONSEFUN, CHIR%RESPONSEFUN, WPOT%RESPONSEFUN, &
            TBSE%RESPONSEFUN, TBSEA%RESPONSEFUN, TVXC%RESPONSEFUN)

    ! switch off all model dielectric functions
    CALL CHECK_FULL_KPOINTS ! all set up properly ?

!=======================================================================
! read second WAVECAR file
!=======================================================================
    LWFROMFILE=.FALSE.
    ! read wavefunctions W_W
    IF (LGW0) THEN
       CALL INWAV_ALTERNATIVE( IO, WDES, W_W, GRID, LATT_CUR, LWFROMFILE, 'chi')
    ENDIF
!
! calculate linear response (only for pead applicable) or try to read WAVEDER.chi file 
!
    IF (LWFROMFILE) THEN
       CALL PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,W_W)
       CALL ORTHCH(WDES,W_W, INFO%LOVERL, LMDIM,CQIJ)
       CALL REDIS_PW_OVER_BANDS(WDES, W_W)
       CALL KPAR_SYNC_ALL(WDES,W_W)
       CALL CLEANUP_CELEN(W_W)
       IF (IO%LOPTICS) THEN
          CALL START_TIMING("G")
          CALL PEAD_RESETUP_WDES(WDES, GRID, KPOINTS, LATT_CUR, LATT_CUR, IO)
          CALL LR_OPTIC( &
               P,WDES,NONLR_S,NONL_S,W_W,LATT_CUR,LATT_INI, &
               T_INFO,INFO,IO,KPOINTS,SYMM,GRID,GRID_SOFT, &
               GRIDC,GRIDUS,C_TO_US,SOFT_TO_C, &
               CHTOT,DENCOR,CVTOT,CSTRF, &
               CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM, &
               CHDEN,SV,LMDIM,IRDMAX,EFERMI,NEDOS, & 
               LSTORE=.TRUE., LPOT=.FALSE.)
          CALL STOP_TIMING("G",IO%IU6,'OPTICS')
       ELSE
          CALL READ_CDER_BETWEEN_STATES(WDES, IO%IU0, 55, 'chi')
       ENDIF
    ELSE
!
! calculate linear response (only for pead applicable) or try to read WAVEDER file
!
       IF (IO%LOPTICS .AND. .NOT. LR_OPTIC_DONE()) THEN
          CALL START_TIMING("G")
          CALL PEAD_RESETUP_WDES(WDES, GRID, KPOINTS, LATT_CUR, LATT_CUR, IO)
          CALL LR_OPTIC( &
               P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
               T_INFO,INFO,IO,KPOINTS,SYMM,GRID,GRID_SOFT, &
               GRIDC,GRIDUS,C_TO_US,SOFT_TO_C, &
               CHTOT,DENCOR,CVTOT,CSTRF, &
               CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM, &
               CHDEN,SV,LMDIM,IRDMAX,EFERMI,NEDOS, & 
               LSTORE=.TRUE., LPOT=.FALSE.)
          CALL STOP_TIMING("G",IO%IU6,'OPTICS')
       ELSE IF (.NOT. LR_OPTIC_DONE()) THEN
          CALL READ_CDER_BETWEEN_STATES(WDES, IO%IU0, 55)  ! read first derivative of wavefuntions
       ENDIF
       W_W=W
    ENDIF

    CALL APPLY_SCISSOR(W, SCISSOR)
    IF (LWFROMFILE) CALL APPLY_SCISSOR(W_W, SCISSOR)
!=======================================================================
! determine frequency grid and complex shift
!=======================================================================
    ! set OMEGATL
    IF (OMEGATL<0 .AND. NODES_IN_DIELECTRIC_FUNCTION>=0) THEN
       OMEGATL=NODES_IN_DIELECTRIC_FUNCTION*10
       IF (LACFDT .AND. (OMEGAGRID>= 100)) THEN
          OMEGATL=MAX_ENERGY_UNOCCUPIED(WDES,W_W)
       ELSE IF (LACFDT) THEN
          OMEGATL=MAX(OMEGATL, MAX_ENERGY_UNOCCUPIED(WDES,W_W)*2.0)
       ENDIF
    ELSEIF (OMEGATL<0) THEN
       IF (LACFDT .AND. (OMEGAGRID>= 100)) THEN
          OMEGATL=MAX_ENERGY_UNOCCUPIED(WDES,W_W)
       ELSE IF (LACFDT) THEN
          OMEGATL=MAX(OMEGATL, MAX_ENERGY_UNOCCUPIED(WDES,W_W)*2.0)
       ENDIF
    ENDIF
    IF (OMEGATL<0) THEN
       OMEGATL=ABS(OMEGATL)
    ENDIF
    IF ( .NOT. LACFDT .OR. OMEGAGRID<100) THEN  
    ! final check on OMEGATL: the frequency grid MUST go at least to 
    ! the maximum transition energy (safeguarded by a factor 1.1)
       OMEGATL=MAX(OMEGATL, MAX_ENERGY_UNOCCUPIED(WDES,W)*1.1)
       OMEGATL=MAX(OMEGATL, MAX_ENERGY_UNOCCUPIED(WDES,W_W)*1.1)
    ENDIF

    ! position of HOMO
    N1=MIN(MAX(LAST_FILLED_OPTICS_NO_MOD(W),LAST_FILLED_OPTICS_NO_MOD(W_W)),WDES%NB_TOT)

    ! set OMEGAMIN (usually the band gap)
    IF (OMEGAMIN<0 .AND. .NOT. LCRPA ) THEN
       CALL DETERMINE_BAND_GAP(WDES, W, WDES%NB_TOT, OMEGAMIN, OMEGATL, EFERMI, NOMEGA)
    ENDIF

    ! set OMEGAMAX
    IF (LQPBSE .AND. OMEGAMAX<=0) THEN
       ! maximum excitation in G from BSE is 2 eps_max - eps_min, eps_max - 2 eps_min 
       ! so this default should be fine
       OMEGAMAX=MAX_ENERGY_OCC_UNOCCUPIED(WDES,W,WDES%NB_TOT)*2.0
    ELSE IF (OMEGAMAX<0 .AND. NOMEGA <= 1) THEN
       OMEGAMAX=0
    ELSE IF (OMEGAMAX<0 .AND. NODES_IN_DIELECTRIC_FUNCTION>=0) THEN
       OMEGAMAX=NODES_IN_DIELECTRIC_FUNCTION/1.3
       ! bottom of valence to top of valence (HOMO) = valence band width
       OMEGAMAX=MAX(OMEGAMAX, MAX_ENERGY_OCC_UNOCCUPIED(WDES,W, N1)/1.3)
       OMEGAMAX=MAX(OMEGAMAX, MAX_ENERGY_OCC_UNOCCUPIED(WDES,W_W, N1)/1.3)
    ENDIF
    IF (OMEGAMAX<0) THEN
       ! bottom of valence to top of valence (HOMO) = valence band width
       OMEGAMAX=ABS(OMEGAMAX)
       OMEGAMAX=MAX(OMEGAMAX, MAX_ENERGY_OCC_UNOCCUPIED(WDES,W, N1)/1.3)
       OMEGAMAX=MAX(OMEGAMAX, MAX_ENERGY_OCC_UNOCCUPIED(WDES,W_W, N1)/1.3)
    ENDIF
    
    ! now set SHIFT only complex frequencies, shift makes no sense whatsoever
    IF (NOMEGAR==0 .AND. NOMEGA/=0 .AND. SHIFT <0) THEN
       SHIFT=0
    ELSE IF (LQPBSE .AND. NOMEGABSE>1 .AND. SHIFT <0) THEN
       SHIFT=OMEGAMAX/NOMEGABSE ! half the grid spacing is usually okay for Lorenzian
    ELSE IF (NOMEGA>2 .AND. SHIFT<0) THEN
       SHIFT=OMEGAMAX*1.3/NOMEGA
       IF (NOMEGA<40) THEN
          SHIFT=OMEGAMAX*1.3/40
       ENDIF
    ENDIF
    IF (SHIFT<0) THEN
       SHIFT=ABS(SHIFT)
    ENDIF

    IF (IO%IU6>=0) THEN
       WRITE(IO%IU6,100) NOMEGA, NOMEGAR, OMEGAMAX, OMEGAMIN, OMEGATL, SHIFT
    ENDIF
    
100  FORMAT(' Response functions by sum over occupied states:'  / &
            ' ==============================================='  / &
         '   NOMEGA  =',I6,  '    number of frequencies'/ & 
         '   NOMEGAR =',I6,  '    number of frequencies along real axis'/ & 
         '   OMEGAMAX=',F7.2,'   maximum frequency'/ & 
         '   OMEGAMIN=',F7.2,'   minimum frequency'/ & 
         '   OMEGATL =',F7.2,'   maximum frequency of tail'/ & 
         '   CSHIFT  =',F7.3,'   complex shift used in evaluation of response functions'/&
         )
!=======================================================================
! set frequencies
! the responsefunction chi(w)
!  chi(w) prop  int 1/ (w-w'+ i delta) + 1/ (-w-w'+ i delta) dw'
!        prop  int (w'-i delta) / ((w'-i delta)^2 + w^2) dw'
! this implies that the broadening is stronger at larger frequencies
!=======================================================================
    ALLOCATE(OMEGA(NOMEGA), OMEGAWEIGHT(NOMEGA), COMEGA(NOMEGA))

    IF (LACFDT .AND. LSPECTRAL )  THEN
       IF (IO%IU0>=0) WRITE(IO%IU0,*) 'forcing LSPECTRAL = .FALSE. for AC-FDT'
       LSPECTRAL = .FALSE.
    ENDIF

    IMAG_GRIDS%BETA= 1/KPOINTS%SIGMA
    IMAG_GRIDS%T%BETA= 1/KPOINTS%SIGMA
    IMAG_GRIDS%B%BETA= 1/KPOINTS%SIGMA
    IMAG_GRIDS%F%BETA= 1/KPOINTS%SIGMA
    IF (NOMEGA==1 .AND. LACFDT ) THEN
       OMEGA(1) =OMEGAMAX
       COMEGA(1)=OMEGA(1)*(0.0_q,1.0_q)
       OMEGAWEIGHT(1)=1
    ELSE IF (NOMEGA==1) THEN
       IF ( NOMEGAR == 0 ) THEN
          OMEGA(1) =OMEGAMAX
          COMEGA(1)=OMEGA(1)*(0._q, 1._q)
       ELSE
         OMEGA(1) =OMEGAMAX
         COMEGA(1)=OMEGA(1)
       ENDIF           
       OMEGAWEIGHT(1) =1
    ELSE
       IF (OMEGAGRID/=0) THEN
          CALL SET_OMEGA_GRID(OMEGAMAX, OMEGAMIN, OMEGATL, KPOINTS%SIGMA, OMEGA(1:NOMEGAR), IO%IU6, VERSION=OMEGAGRID)
       ELSE
          CALL SET_OMEGA_GRID(OMEGAMAX, OMEGAMIN, OMEGATL, KPOINTS%SIGMA, OMEGA(1:NOMEGAR), IO%IU6)
       ENDIF
       
       COMEGA(1:NOMEGAR)=OMEGA(1:NOMEGAR)

       ! determine weights along imaginary axis 
       IF (NOMEGAR<NOMEGA) THEN
          IF (OMEGAGRID>=100 .AND. OMEGAGRID < 170 ) THEN
             IF (NOMEGA>32 .AND. OMEGAGRID < 160  ) THEN
                CALL vtutor%write(isError, NOmega32)
             ENDIF
             !call minimax reader
             CALL SET_IMAG_GRID_HANDLE(OMEGAMIN, OMEGATL, NOMEGA, OMEGAGRID, &
                KPOINTS%SIGMA, IMAG_GRIDS, WDES%COMM, IO)

             IF ( LLTDMP2 ) THEN
                ! copy time grids into omega points for laplace transofrmed direct MP2 energy 
                OMEGA = IMAG_GRIDS%TAU
                ! make up for 2pi factor in RPA routine
                OMEGAWEIGHT = IMAG_GRIDS%TAU_WEIGHT*2
             ELSE
                IF ( LGW ) THEN   
                   OMEGA=IMAG_GRIDS%FER_RE
                   OMEGAWEIGHT=IMAG_GRIDS%FER_RE_WEIGHT
                ELSE 
                   OMEGA=IMAG_GRIDS%BOS_RE
                   ! factor 1/4 is due to definit of object function function
                   OMEGAWEIGHT=IMAG_GRIDS%BOS_RE_WEIGHT
                ENDIF     
             ENDIF

             IF (OMEGA(NOMEGA)<OMEGATL/1.1.AND. .NOT.LLTDMP2) THEN
                CALL vtutor%write(isAlert, OmegaGridOmegaTL)
             ENDIF
          ELSE IF (OMEGAGRID/=0) THEN
             CALL SET_OMEGA_GRID(OMEGAMAX, OMEGAMIN, OMEGATL, KPOINTS%SIGMA, OMEGA(NOMEGAR+1:NOMEGA), IO%IU6, & 
                  OMEGAWEIGHT(NOMEGAR+1:NOMEGA), VERSION=OMEGAGRID)
          ELSE
             CALL SET_OMEGA_GRID(OMEGAMAX, OMEGAMIN, OMEGATL, KPOINTS%SIGMA, OMEGA(NOMEGAR+1:NOMEGA), IO%IU6, & 
                  OMEGAWEIGHT(NOMEGAR+1:NOMEGA), VERSION=40)
          ENDIF
          COMEGA(NOMEGAR+1:NOMEGA)=OMEGA(NOMEGAR+1:NOMEGA)*(0.0_q,1.0_q)
       ENDIF
    ENDIF


    IF ( .NOT. LLTDMP2 ) THEN 
       IF (IO%IU0>=0) WRITE(IO%IU0,'("energies w= ",/,(8F9.3))') COMEGA(1:NOMEGA)
       IF (IO%IU0>=0) WRITE(17,'("energies w= ",/,(8F9.3))') COMEGA(1:NOMEGA)
       IF (IO%IU6>=0) THEN
          WRITE(IO%IU6,'(/," Energies iw= ",/," ============",/,(8F9.3))') COMEGA(1:NOMEGA)
          WRITE(IO%IU6,*)
       ENDIF
    ELSE
       IF (IO%IU0>=0) WRITE(IO%IU0,'("times tau= ",/,(8F9.3))') COMEGA(1:NOMEGA)
       IF (IO%IU0>=0) WRITE(17,'("times tau= ",/,(8F9.3))') COMEGA(1:NOMEGA)
       IF (IO%IU6>=0) THEN
          WRITE(IO%IU6,'(/," Times i tau= ",/," ============",/,(8F9.3))') COMEGA(1:NOMEGA)
          WRITE(IO%IU6,*)
       ENDIF
    ENDIF

    IF (LQPBSE .AND. NBANDSO <= 0) THEN
       NBANDSO=MAX(NBANDSGW,LAST_FILLED_OPTICS(W))
    ELSEIF ((LADDER .OR. LBSE) .AND. NBANDSO <= 0) THEN
       NBANDSO=LAST_FILLED_OPTICS(W)
    ENDIF
    IF (LQPBSE .AND. NBANDSV <= 0) THEN
       NBANDSV=WDES%NB_TOT   ! use all virtual (unoccupied) orbitals
    ELSEIF ((LADDER .OR. LBSE) .AND. NBANDSV <= 0) THEN
       NBANDSV=LAST_FILLED_OPTICS(W)
    ENDIF

    IF (LHFCORRECT) THEN
      CALL FOCK_ACFDT(GRID, LATT_CUR, W, LMDIM, NONLR_S, NONL_S, P , OMEGAWEIGHT, OMEGA, IO%IU6)
      CALL vtutor%stopCode()
    ENDIF
!=======================================================================
! generate wavefunctions at all k-points
! this simplifies the calculation of Xi_q(r,r',w) since
! q can be selected to lie in the irreducible wedge of the BZ
! and one can simply calculate all k-point pairs
! k_a and k_i=k_a+q without further symmetry considerations.
! certainly it takes more memory, but compared to the responsefunction
! this is almost negligible
!=======================================================================
    IF (SYMM%ISYM>=0 .AND. .NOT. LGAMMA) THEN
    ! switch off symmetry
      CALL NOSYMM(LATT_CUR%A,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS,SYMM%PTRANS, &
           SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%MAGROT,INFO%ISPIN,IO%IU6)
    ! reread k-points with LINVERSION=.FALSE.  to generate full mesh
      CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR,.FALSE., &
           T_INFO%NIONS,SYMM%ROTMAP, SYMM%MAGROT, SYMM%ISYM, IO%IU6,IO%IU0)
      CALL KPAR_SYNC_ALL(WDES,W)
      CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
! TODO: I really fail to understand this, why is this required
      CALL PEAD_RESETUP_WDES(WDES, GRID, KPOINTS, LATT_CUR, LATT_CUR, IO) 
      CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
      CALL RESETUP_FOCK( WDES, LATT_CUR) 
      IF (LWFROMFILE) &
           CALL REALLOCATE_WAVE( W_W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
    ENDIF

    CALL KINETIC_ENERGY(W)
    IF (LWFROMFILE) CALL KINETIC_ENERGY(W_W)

    IF (LACFDT) THEN
! apply soft band cutoff (does not seem to improve AC-FDT results
!  hence commented out for the time being)
       IF (TELESCOPE>0) THEN
          CALL VIRTUAL_BAND_CUTOFF(W, WDES%ENMAX, TELESCOPE)
          IF (LWFROMFILE)  CALL VIRTUAL_BAND_CUTOFF(W_W, WDES%ENMAX, TELESCOPE)
       ELSE
          CALL VIRTUAL_BAND_CUTOFF(W)
          IF (LWFROMFILE) CALL VIRTUAL_BAND_CUTOFF(W_W)
       ENDIF
       IF (ENCUTGW<=0) THEN
          CALL XI_ACFDT_SETUP( COR, WDES%ENMAX , ENCUTGWSOFT) 
       ELSE
          CALL XI_ACFDT_SETUP( COR, ENCUTGW, ENCUTGWSOFT) 
       ENDIF
    ELSE
       CALL VIRTUAL_BAND_CUTOFF(W)
       IF (LWFROMFILE) CALL VIRTUAL_BAND_CUTOFF(W_W)
    ENDIF

    IF (LscQPGW .AND. LGW0 .AND. .NOT. LWFROMFILE) THEN
       CALL ALLOCW( W%WDES, W_W)
       W_W%CPTWFP=W%CPTWFP
       W_W%CPROJ =W%CPROJ 
       W_W%CELTOT=W%CELTOT
       W_W%FERTOT=W%FERTOT
       W_W%AUXTOT=W%AUXTOT
       W_W%EFERMI=W%EFERMI
    ELSE IF (.NOT. LWFROMFILE) THEN
       W_W=W
       IF (LGW0 .AND. .NOT. LscQPGW) THEN
          ALLOCATE(W_W%CELTOT(W_W%WDES%NB_TOT,W_W%WDES%NKPTS,W_W%WDES%ISPIN), &
               W_W%FERTOT(W_W%WDES%NB_TOT,W_W%WDES%NKPTS,W_W%WDES%ISPIN))
          W_W%FERWE => W_W%FERTOT(W_W%WDES%NB_LOW:W_W%WDES%NB_TOT:W_W%WDES%NB_PAR,:,:)
          W_W%CELEN => W_W%CELTOT(W_W%WDES%NB_LOW:W_W%WDES%NB_TOT:W_W%WDES%NB_PAR,:,:)
          W_W%CELTOT=W%CELTOT
          W_W%FERTOT=W%FERTOT
          W_W%EFERMI=W%EFERMI
       ENDIF
    ENDIF
    
    ALLOCATE (FERTOT_INITIAL(W_W%WDES%NB_TOT,W_W%WDES%NKPTS,W_W%WDES%ISPIN))
    FERTOT_INITIAL=W%FERTOT

    IF(IO%IU6>=0) WRITE(IO%IU6,'(A)', ADVANCE="No") ' files read and symmetry switched off, memory is now:'
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
       ! gamma only data layout with densities stored as real in real space
       GRIDWGW%LREAL=.TRUE.
    ENDIF

    CALL GEN_INDEX (GRIDWGW, WGW, LATT_CUR%B, LATT_CUR%B,IO%IU6,-1, .TRUE.)
    !  init FFT (required if real to complex FFT is used)
    CALL FFTINI(WGW%NINDPW(1,1), WGW%NGVECTOR(1), WGW%NKPTS, WGW%NGDIM, GRIDWGW)

    IF(IO%IU6>=0) WRITE(IO%IU0,'(A,I8)') ' responsefunction array rank=',WGW%NGDIM

!
! TODO: better always pop the XC, but care must be taken, that the 
! LDA part of the functional is not included if LGWLF is set
!
    IF (.NOT. LGWLF ) THEN
       CALL POP_XC_TYPE
       IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
          CALL SETUP_LDA_XC(2,IO%IU6,IO%IU0,IO%IDIOT)
       ELSE
          CALL SETUP_LDA_XC(1,IO%IU6,IO%IU0,IO%IDIOT)
       ENDIF
    ENDIF

! soft cutoff on kinetic energy
    IF (ENCUTGW /= ENCUTGWSOFT .AND. ENCUTGW > 0 .AND. ENCUTGWSOFT >0 ) THEN
       CALL XI_DATAKE( WGW, LATT_CUR, ENCUTGW, ENCUTGWSOFT)
    ELSE
       CALL XI_DATAKE( WGW, LATT_CUR)
    ENDIF

    ! determined the q-points that are used for the response function
    CALL SETUP_IRZ_MAP(S2E, WGW, IO%IU0, IO%IU6)

    CALL SETUP_KPOINTS_ESF( LATT_CUR, SYMM, WDES, S2E, T_INFO, DYN, IO )

!=======================================================================
! the contributions to Xi are calculated considering NSTRIP
! valence bands at the same time (on each core)
!=======================================================================
    NSTRIP=1
    DO ISP=1,WDES%ISPIN
       DO K1=1,WDES%NKPTS
          NSTRIP=MAX(NSTRIP,(LAST_FILLED_XI(W,K1,ISP,IMAG_GRIDS%T%BETA)+WDES%NB_PAR-1)/WDES%NB_PAR)
       ENDDO
    ENDDO
    NSTRIP=MAX(MIN(NSTRIP_STANDARD,NSTRIP),1)

    ! memory estamete for job
    CALL DUMP_REQUIRED_STORAGE( WDES, WGW, S2E, IO )

!    NSTRIP=1
!=======================================================================
! read screened two electron integrals or prepare for their calculation
!=======================================================================

    IF (LGW) THEN
       IF (LSPECTRAL) THEN
          IF (LSELFENERGY ) THEN
             ALLOCATE(SCREENED_TWO_ELECTRON_INTEGRALS(NBANDSGW, KPOINTS_ORIG%NKPTS, WDES%NBANDS, MAXVAL(S2E%K2_STORE_INDEX),  &
                  WDES%ISPIN, NOMEGA*2), STAT = ISTAT )
             ISP = 16._q*NBANDSGW*KPOINTS_ORIG%NKPTS*WDES%NBANDS*MAXVAL(S2E%K2_STORE_INDEX)*WDES%ISPIN*NOMEGA*2
          ELSE
             ALLOCATE(SCREENED_TWO_ELECTRON_INTEGRALS(NBANDSGW, KPOINTS_ORIG%NKPTS, WDES%NBANDS, MAXVAL(S2E%K2_STORE_INDEX),  &
                  WDES%ISPIN, 2), STAT = ISTAT )
             ISP = 16._q*NBANDSGW*KPOINTS_ORIG%NKPTS*WDES%NBANDS*MAXVAL(S2E%K2_STORE_INDEX)*WDES%ISPIN*2
          ENDIF
       ELSE
          ALLOCATE(SCREENED_TWO_ELECTRON_INTEGRALS(NBANDSGW, KPOINTS_ORIG%NKPTS, WDES%NBANDS, MAXVAL(S2E%K2_STORE_INDEX),  &
               WDES%ISPIN, NOMEGA), STAT = ISTAT )
          ISP = 16._q*NBANDSGW*KPOINTS_ORIG%NKPTS*WDES%NBANDS*MAXVAL(S2E%K2_STORE_INDEX)*WDES%ISPIN*NOMEGA
       ENDIF
       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "CALCULATE_XI (SCREENED_TWO_ELECTRON_INTEGRALS) is not able to allocate "//&
             str(ISP/1024) //" kB of data on MPI rank 0." )
       ENDIF
       CALL REGISTER_ALLOCATE(8._q*SIZE(SCREENED_TWO_ELECTRON_INTEGRALS), "2eintegral")
    ENDIF

    ! determine in which parts the spectral method is used
    IF (LSPECTRAL) THEN
       LSPECTRALCHI=.TRUE.       ! response function using spectral method
    ELSE
       LSPECTRALCHI=.FALSE.
    ENDIF

!=======================================================================
! allocation
! determine how many frequencies can be done in one go
! MAXMEM Mbyte is coded as the upper memory limit
!=======================================================================
    NELM=0

    IF (LSPECTRAL .AND. S2E%NUMBER_OF_NQ==1 .AND. LGW0) THEN
       LCHIREALLOCATE=.FALSE.  ! CHI is calculated once and kept fixed during electronic iterations
    ELSE
       LCHIREALLOCATE=.TRUE.
    ENDIF

2000 NELM=NELM+1

    IF (LGW) SCREENED_TWO_ELECTRON_INTEGRALS=0

    ! allocate memory for sc GW
    IF (LscQPGW) THEN
       CALL ALLOCW(WDES,WACC1)
       CALL ALLOCW(WDES,WACC2)
       IF (SAVE_CACHE_MEMORY) THEN
          ! save memory by storing the action of the selfenergy on a maximum 
          ! of NSTRIP times number of cores orbitals
          CALL ALLOCATE_CACHER(WDES, WCACHE, LMDIM, NSTRIP*WDES%NB_PAR)
       ELSE
          CALL ALLOCATE_CACHER(WDES, WCACHE, LMDIM, NBANDSGW)
       ENDIF
    ELSE
       NULLIFY(WCACHE)
    ENDIF

    ! set number of frequency points for non-cached response function
    CALL SET_NOMEGA_CHI( WGW, NCPU, NOMEGA_CHI ) 

    CHI %SHIFT =SHIFT
    CHIR%SHIFT =SHIFT
    CHI0%SHIFT =SHIFT
    WPOT%SHIFT =SHIFT
    TVXC%SHIFT =SHIFT
    TBSE%SHIFT =SHIFT
    TBSEA%SHIFT=SHIFT

    ! also set inverse temperature  
    CHI %BETA = IMAG_GRIDS%BETA
    CHIR%BETA = IMAG_GRIDS%BETA
    CHI0%BETA = IMAG_GRIDS%BETA
    WPOT%BETA = IMAG_GRIDS%BETA
    TVXC%BETA = IMAG_GRIDS%BETA
    TBSE%BETA = IMAG_GRIDS%BETA
    TBSEA%BETA= IMAG_GRIDS%BETA

    IF ( SIZE( OMEGA ) > 0 ) THEN
       ! and zero frequency point
       CHI  %OMEGA0 = OMEGA(1)
       CHIR %OMEGA0 = OMEGA(1)
       CHI0 %OMEGA0 = OMEGA(1)
       WPOT %OMEGA0 = OMEGA(1)
       TVXC %OMEGA0 = OMEGA(1)
       TBSE %OMEGA0 = OMEGA(1)
       TBSEA%OMEGA0 = OMEGA(1)
    ENDIF

    IF (LCHIREALLOCATE .OR. NELM==1) THEN
       IF(IO%IU0>=0) WRITE(IO%IU0,'(A,I4,A,I6)') ' allocating', NOMEGA_CHI,' responsefunctions rank=',WGW%NGDIM
       CALL ALLOCATE_RESPONSEFUN(CHI, WGW%NGDIM, WGW%LGAMMA, LACFDT .AND. WGW%LGAMMA, NOMEGA_CHI)
    ENDIF

    IF (LSPECTRAL .AND. LADDER .AND. NBANDSO>0) THEN
       LCHIR=.TRUE.
       IF (IO%IU0>=0) WRITE(IO%IU0,*) 'resonant part calculated as well'
       IF(IO%IU0>=0) WRITE(IO%IU0,'(A,I4,A,I6)') ' allocating', NOMEGA_CHI,' responsefunctions rank=',WGW%NGDIM
       CALL ALLOCATE_RESPONSEFUN(CHIR, WGW%NGDIM, WGW%LGAMMA, LACFDT .AND. WGW%LGAMMA, NOMEGA_CHI)
       CALL CHI_KRAM_KRON_RES_TABLE( TABLE_RES, OMEGA, SHIFT)
    ELSE
       LCHIR=.FALSE.
    ENDIF

    ! second response function array is required if spectral method is used 
    ! or in parallel version
    IF ((LSPECTRAL .OR. NCPU>1) .AND. NOMEGA>0) THEN
       ! set number of frequency points for cached response function
       CALL DETERMINE_NOMEGA_WPOT( WGW, NCPU, NOMEGA_CHI, NOMEGA_WPOT ) 

#ifdef use_shmem
       IF(IO%IU0>=0) WRITE(IO%IU0,'(A,I4,A,I6)') ' shmem allocating', NOMEGA_WPOT,' responsefunctions rank=',WGW%NGDIM
       IF(IO%IU0>=0) WRITE(IO%IU0,'(A,I4,A,I6)') ' response function shared by NCSHMEM nodes ',WGW%COMM_SHMEM%NCPU
       IF(IO%IU6>=0) WRITE(IO%IU6,'(A,I4,A,I6)') ' shmem allocating', NOMEGA_WPOT,' responsefunctions rank=',WGW%NGDIM
       IF(IO%IU6>=0) WRITE(IO%IU6,'(A,I4,A,I6)') ' response function shared by NCSHMEM nodes ',WGW%COMM_SHMEM%NCPU

       CALL ALLOCATE_RESPONSEFUN_SHMEM(WPOT, WGW%NGDIM, WGW%LGAMMA, LACFDT .AND. WGW%LGAMMA, NOMEGA_WPOT, & 
            WGW%COMM_SHMEM, IO%IU0, IO%IU6 )

       CALL COMM_INTER_SHMEM( WPOT, WGW, INNER_LOOP_PARALLEL)
#else
       IF(IO%IU0>=0) WRITE(IO%IU0,'(A,I4,A,I6)') ' allocating', NOMEGA_WPOT,' responsefunctions rank=',WGW%NGDIM
       IF(IO%IU6>=0) WRITE(IO%IU6,'(A,I4,A,I6)') ' allocating', NOMEGA_WPOT,' responsefunctions rank=',WGW%NGDIM
       CALL ALLOCATE_RESPONSEFUN(WPOT, WGW%NGDIM, WGW%LGAMMA, LACFDT .AND. WGW%LGAMMA, NOMEGA_WPOT )
#endif
       ! allocate cache structure
    ELSE
       NOMEGA_WPOT=0
       WPOT=CHI
    ENDIF
    ! need cache structure for spectral method as well as for scGW
    IF (LSPECTRAL .OR. LscQPGW) THEN
       CALL ALLOCATE_RESPONSEFUN_CACHE( WPOT, MAXCACHE)
    ELSE
       NULLIFY(WPOT%NCACHE)
    ENDIF

! need cache structure to evaluate selfenergy between states
    IF (.NOT.ASSOCIATED(WPOT%NCACHE) .AND. LGW) THEN
       CALL ALLOCATE_RESPONSEFUN_CACHE( WPOT, MAXCACHE)
    ENDIF

    IF (LSPECTRAL) THEN
       CALL CHI_KRAM_KRON_TABLE( TABLE, OMEGA, SHIFT)

       IF (LSPECTRALGW) THEN
          CALL POT_HILBERT_TABLE_SPECTRAL( TABLE_POT_PLUS, OMEGA, SHIFT, 1)
          CALL POT_HILBERT_TABLE_SPECTRAL( TABLE_POT_MIN , OMEGA, SHIFT, -1)
       ELSE
          CALL POT_HILBERT_TABLE( TABLE_POT_PLUS, OMEGA, SHIFT, 1)
          CALL POT_HILBERT_TABLE( TABLE_POT_MIN , OMEGA, SHIFT, -1)
       ENDIF
    ENDIF

    IF (ANTIRES <0) THEN
       IF (IO%IU0>=0) WRITE(IO%IU0,*) 'only resonant part in response function'
       CALL CHI_KRAM_KRON_RES_TABLE( TABLE, OMEGA, SHIFT)
    ENDIF


    IF (LOEP .OR. LFXCEPS) THEN
       ! allocate response function for w=0 (required for OEP method)
       CALL ALLOCATE_RESPONSEFUN(CHI0, WGW%NGDIM, WGW%LGAMMA, LACFDT .AND. WGW%LGAMMA, 1)
       CHI0%COMEGA=-1
    ENDIF

    IF (IO%IU0>=0) THEN
       WRITE(IO%IU0,*) 'Doing ',NOMEGA_CHI,' frequencies on each core in blocks of ',NOMEGA_WPOT
    ENDIF

    IF (IO%LOPEN) CALL WFORCE(IO%IU6)

!=======================================================================
! calculate correlation energy from Bethe Salpeter equation
!=======================================================================
    IF (LBSE .AND. LCORBSE) THEN
       CALL CALCULATE_COR_BSE()
       CALL DUMP_FINAL_TIMING(IO%IU6)
       CALL XML_CLOSE_TAG('calculation')
       CALL STOP_XML
#ifdef PROFILING
       PROFILING_STOP('calculate_xi')
       ! started in INIT_PROFILING
       PROFILING_STOP('total_time')
       IF (IO%IU0>=0) THEN
          CALL PRINT_FULL_PROFILE(IO%IU6,PFLAT=.FALSE.)
       ENDIF
#endif
       CALL vtutor%stopCode()
    ENDIF
!=======================================================================
! calculate Green's function and QP energies from Bethe Salpeter calculation 
!=======================================================================
    IF (LBSE .AND. LQPBSE) THEN
       CALL CALCULATE_QPBSE()
       GOTO 1000
    ENDIF
!=======================================================================
! Bethe Salpeter this will also terminate VASP
!=======================================================================
    IF (LBSE) THEN
       ! sets DELTA_COND in local_field.F
       CALL CALCULATE_LOCAL_FIELD_PREPARE( 1, W, WGW, LBSE, LGWLF, & 
            .FALSE.,  & 
            IDIR_MAX, COMEGA, TBSE, TBSEA, ANTIRES, & 
            LATT_CUR, NKREDLFX, NKREDLFY, NKREDLFZ, IO%IU0, IO%IU6 )

       IF (KPOINT_BSE(1)==-1) THEN
          ! select first vector in the set of difference vectors
          ! that is usually the Gamma-point
          NQ=S2E%NQ(1)
       ELSE
          ! user selected k-point
          NQ=KPOINT_BSE(1)
       ENDIF

       ! if LGWLF is set, ALDAX and ALDAC are zero at this point
       ! hence the DFT local field effects f_xc (stored in FXCR) are exactly zero
       CALL CALCULATE_LOCAL_FIELD_DFT( &
            GRIDC, WDES, LATT_CUR, CHTOT, DENCOR, &
            WGW, NQ, 1 )
       
       IF (WDES%ISPIN >1) THEN
          ! take into account all spins
          CALL CALCULATE_BSE( &
               IBSE,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
               T_INFO,DYN,INFO,IO,KPOINTS,SYMM,GRID,LMDIM,CQIJ, &
               WGW, SHIFT, CELTOT_HARTREE_KINETIC, CELTOT_X, & 
               NBANDSO, NBANDSV, OMEGAMAX, NQ, KPOINT_BSE(2:4), -1, &
               LHARTREE, LADDER, LTRIPLET, ANTIRES, LGWLF, LFXC, NEDOS, NELMGW, NBANDSGW, L2ORDER, &
               LCORBSE, LQPBSE, LTCTC)
       ELSE
          ! take into account only first (and only) spin component
          CALL CALCULATE_BSE( &
                IBSE,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
                T_INFO,DYN,INFO,IO,KPOINTS,SYMM,GRID,LMDIM,CQIJ, &
                WGW, SHIFT, CELTOT_HARTREE_KINETIC, CELTOT_X, & 
                NBANDSO, NBANDSV, OMEGAMAX, NQ, KPOINT_BSE(2:4), 1, &
                LHARTREE, LADDER, LTRIPLET, ANTIRES, LGWLF, LFXC, NEDOS, NELMGW, NBANDSGW, L2ORDER, &
                LCORBSE, LQPBSE, LTCTC)
       ENDIF

       CALL XML_CLOSE_TAG('calculation')
       CALL DUMP_FINAL_TIMING(IO%IU6)
       CALL STOP_XML
#ifdef PROFILING
       PROFILING_STOP('calculate_xi')
       ! started in INIT_PROFILING
       PROFILING_STOP('total_time')
       IF (IO%IU0>=0) THEN
          CALL PRINT_FULL_PROFILE(IO%IU6,PFLAT=.FALSE.)
       ENDIF
#endif
       CALL vtutor%stopCode()
    ENDIF

    IF (LADDER) &
      CALL CALCULATE_LOCAL_FIELD_PREPARE( NBANDSO, W, WGW, .FALSE., LGWLF,  & 
       .NOT. LGW .AND. .NOT. LACFDT , &  ! determines whether CB is shifted, not applied for GW or ACFDT
       IDIR_MAX, COMEGA, TBSE, TBSEA, ANTIRES, & 
       LATT_CUR, NKREDLFX, NKREDLFY, NKREDLFZ,  IO%IU0, IO%IU6  )

    IF(IO%IU6>=0) WRITE(IO%IU6,'(A)', ADVANCE="No") ' all allocation done, memory is now:'
    CALL DUMP_ALLOCATE(IO%IU6)

!======================================================================
! constrained RPA
!======================================================================
    IF ( LCRPA ) &
       CALL CRPA_SETUP( WDES, W, W_CORR, LMDIM, T_INFO, P, CQIJ,&
            LATT_CUR, GRID, INFO, KPOINTS, S2E%NUMBER_OF_NQ, IO )  

!=======================================================================
! loop over spin, q-points in the BZ and frequencies
!=======================================================================
    IF (LscQPGW) THEN 
       WACC1%CW=0
       WACC2%CW=0
    ENDIF

    NOPER_=0 ; NFLOAT_=0

qpoints: DO NQ_COUNTER=1,S2E%NUMBER_OF_NQ
#ifdef MPI   
    IF (OUTER_LOOP_PARALLEL .AND. MOD(NQ_COUNTER-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
    NQ=S2E%NQ(NQ_COUNTER)
    CALL SET_RESPONSE_KPOINT(CHI, WDES%VKPT(:,NQ), NQ)
    IF (LCHIR) CALL SET_RESPONSE_KPOINT(CHIR, WDES%VKPT(:,NQ), NQ)
    CALL SET_RESPONSE_KPOINT(WPOT, WDES%VKPT(:,NQ), NQ)
    CALL SETWDES(WGW, WGWQ, NQ )

    IF (IO%IU0>=0) WRITE(IO%IU0,'("NQ=",I4,3F10.4,", ")') NQ,CHI%VKPT
    IF (IO%IU0>=0) WRITE(17,'("NQ=",I4,3F10.4,", ")') NQ,CHI%VKPT
    
130 FORMAT (5X, //, &
         &'----------------------------------------------------', &
         &'----------------------------------------------------'//)
    IF (IO%IU6>=0) THEN
       WRITE(IO%IU6,130)
       WRITE(IO%IU6,'("NQ=",I4,3F10.4,", ")') NQ,CHI%VKPT
    ENDIF

lchi: DO NOMEGA_INDEX=1, NOMEGA, SIZE(CHI%COMEGA)*NCPU

#ifdef MPI   
       IF (INNER_LOOP_PARALLEL) THEN
          CALL SET_RESPONSEFUN_FREQ_KSUM( CHI, WGW, COMEGA, NOMEGA_INDEX)
          IF (LCHIR) CALL SET_RESPONSEFUN_FREQ_KSUM( CHIR, WGW, COMEGA, NOMEGA_INDEX)
       ELSE
#endif
          CALL SET_RESPONSEFUN_FREQ( CHI, WGW, COMEGA, NOMEGA_INDEX)
          IF (LCHIR) CALL SET_RESPONSEFUN_FREQ( CHIR, WGW, COMEGA, NOMEGA_INDEX)
#ifdef MPI   
       ENDIF
#endif

       CALL START_TIMING("GW")

       NOPER=0 ; NFLOAT=0

calc_chi: IF (LCHIREALLOCATE .OR. NELM==1) THEN
! this test allows to check whether there is any read/write to CHI%RESPONSE 
! if NOMEGA is zero (test passed for ACFDT)
!      IF (CHI%NOMEGA==0) CALL DEALLOCATE_RESPONSEFUN( CHI )
       IF (CHI%NOMEGA>0) CALL CLEAR_RESPONSE(CHI)
       IF (LCHIR .AND. CHIR%NOMEGA>0) CALL CLEAR_RESPONSE(CHIR)
!-----------------------------------------------------------------------
! determine response function
!-----------------------------------------------------------------------
       CALL START_TIMING("G")


       NOMEGA_INDEXI=0
       DO
          IF( SET_RESPONSEFUN_FREQI( WPOT, WGW, COMEGA, NOMEGA_INDEXI, & 
               NOMEGA_INDEX, NOMEGA_INDEX+ SIZE(CHI%COMEGA)*NCPU-1, LSPECTRALCHI)) EXIT

          IF (LSPECTRAL .OR. NCPU>1) CALL CLEAR_RESPONSE_SHMEM(WPOT)

          ! loop over all k-points K1 (index a)
          NK=0
          DO K1=1,WDES%NKPTS
             IF (ODDONLYGW .AND. ABS(MODULO(NINT(WDES%VKPT(1,K1)*KPOINTS%NKPX+ & 
                  WDES%VKPT(2,K1)*KPOINTS%NKPY+ & 
                  WDES%VKPT(3,K1)*KPOINTS%NKPZ),2))<=1E-6) THEN
                CYCLE
             ENDIF
             IF (EVENONLYGW   .AND. ABS(MODULO(NINT(WDES%VKPT(1,K1)*KPOINTS%NKPX+ & 
                  WDES%VKPT(2,K1)*KPOINTS%NKPY+ & 
                  WDES%VKPT(3,K1)*KPOINTS%NKPZ+1),2))<=1E-6) THEN
                CYCLE
             ENDIF

             NK=NK+1
#ifdef MPI   
             IF (INNER_LOOP_PARALLEL .AND. MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
             DO ISP=1,WDES%ISPIN
             ! loop over all bands (index a) in blocks of NSTRIP
             N1_LAST=LAST_FILLED_XI(W_W,K1,ISP,CHI%BETA)/WDES%NB_PAR
             DO N1=1,N1_LAST,NSTRIP
                CALL GWPROGRESS(IO%IU0, K1, WDES%NKPTS, N1, N1_LAST)
                ! determine upper bound (avoid to go beyond the total
                ! number of bands
                NSTRIP1_ACT=MIN(LAST_FILLED_XI(W_W,K1,ISP,CHI%BETA)/WDES%NB_PAR+1-N1,NSTRIP)

                ! add contributions to Xi
                IF (LSPECTRALCHI) THEN
                   CALL ADD_XI_SPECTRAL( LMDIM, LATT_CUR, W_W, WGWQ, &
                        H, P, ISP,  &
                        K1, N1, NSTRIP1_ACT, WPOT, OMEGA, NBANDSGWLOW, NOPER, NFLOAT )
                ELSEIF( LCRPA ) THEN
                   ! in case of disentangled CRPA calculation W_CORR is the disentangled wavefunction
                   IF ( LDISENTANGLED ) THEN
                       LSUBTRACT_CHIC=.FALSE.
                       CALL ADD_XI( LMDIM, LATT_CUR, W_CORR, WGWQ, &
                            H, P, ISP, KPOINTS%ISMEAR, KPOINTS%SIGMA,&
                            K1, N1, NSTRIP1_ACT, NSTRIP_TOTAL, WPOT, NBANDSGWLOW, &
                            NOPER, NFLOAT)
                   ! in case of standard CRPA calculations W_CORR is the correlated wavefunction
                   ELSE
                       LSUBTRACT_CHIC=.FALSE.
                       CALL ADD_XI( LMDIM, LATT_CUR, W_W, WGWQ, &
                            H, P, ISP, KPOINTS%ISMEAR, KPOINTS%SIGMA,&
                            K1, N1, NSTRIP1_ACT, NSTRIP_TOTAL, WPOT, NBANDSGWLOW, &
                            NOPER, NFLOAT)
                       IF ( LCHIC ) THEN !remove correlated part
                          LSUBTRACT_CHIC=.TRUE.
                          CALL ADD_XI( LMDIM, LATT_CUR, W_CORR, WGWQ, &
                               H, P, ISP, KPOINTS%ISMEAR, KPOINTS%SIGMA,&
                               K1, N1, NSTRIP1_ACT, NSTRIP_TOTAL, WPOT, NBANDSGWLOW, &
                               NOPER, NFLOAT) 
                       ENDIF 
                   ENDIF 
                ELSE
                   LSUBTRACT_CHIC=.FALSE.
                   CALL ADD_XI( LMDIM, LATT_CUR, W_W, WGWQ, &
                        H, P, ISP, KPOINTS%ISMEAR, KPOINTS%SIGMA, &
                        K1, N1, NSTRIP1_ACT, NSTRIP_TOTAL, WPOT, NBANDSGWLOW, NOPER, NFLOAT)
                ENDIF
             ENDDO
             ENDDO
          ENDDO
!$OMP PARALLEL
          CALL CLEAN_RESPONSEFUNCTION_CACHE(WPOT, 1, WPOT%NOMEGA)
!$OMP END PARALLEL
          IF (MODULO(WDES%NKPTS,NK)/=0) THEN
             CALL vtutor%error("The ODDONLYGW flag results in a problem since the total number of &
                &k-points is not dividable by the visited k-points " // str(WDES%NKPTS) // " " // str(I))
          ENDIF

          ! proper weighting if the k points were restricted somehow
          WPOT%RESPONSEFUN=WPOT%RESPONSEFUN*(REAL(WDES%NKPTS,q)/NK)
          WPOT%HEAD=WPOT%HEAD*(REAL(WDES%NKPTS,q)/NK)
          WPOT%WING=WPOT%WING*(REAL(WDES%NKPTS,q)/NK)
          WPOT%CWING=WPOT%CWING*(REAL(WDES%NKPTS,q)/NK)
          CALL STOP_TIMING("G",IO%IU6,"CHI")

          IF (LSPECTRALCHI) THEN
             ! determine Hermitian part of response function from anti-Hermitian part
             ! using Kramers Kronig transform
             CALL ADD_DRUDE_IMAG(WPOT, WGW, WPLASMON, COMEGA, LATT_CUR%OMEGA, INNER_LOOP_PARALLEL )
             CALL DO_CHI_SUM(WPOT, WGWQ,INNER_LOOP_PARALLEL )
             CALL DO_CHI_KRAM_KRON_TABLE(WPOT, CHI, WGW, TABLE )
             IF (LCHIR) CALL DO_CHI_KRAM_KRON_TABLE(WPOT, CHIR, WGW, TABLE_RES )
             CALL STOP_TIMING("G",IO%IU6,"KRAMKRO")
          ELSE IF (NCPU>1) THEN
             CALL SCATTER_FREQU(WPOT, CHI, WGWQ, INNER_LOOP_PARALLEL) 
             CALL STOP_TIMING("G",IO%IU6,"SCATTER")
          ENDIF
       ENDDO
       CALL GWPROGRESS(IO%IU0, WDES%NKPTS, WDES%NKPTS, N1_LAST, N1_LAST)

       IF (.NOT. LSPECTRALCHI .AND. LNLRPA ) THEN
          CALL ADD_DRUDE_REAL(CHI, WGW, WPLASMON, COMEGA, LATT_CUR%OMEGA, SHIFT, INNER_LOOP_PARALLEL)
          IF (LCHIR) CALL ADD_DRUDE_REAL(CHIR, WGW, WPLASMON, COMEGA, LATT_CUR%OMEGA, SHIFT, INNER_LOOP_PARALLEL)
       ENDIF

!-----------------------------------------------------------------------
! calculate the vertex corrections from local exchange correlation pot.
! and screened exchange potential
!-----------------------------------------------------------------------
       IF (NOMEGA_INDEX==1) THEN
          IF (LFXC) THEN
             CALL CALCULATE_LOCAL_FIELD_DFT( &
             GRIDC, WDES, LATT_CUR, CHTOT, DENCOR, &
             WGW, NQ, 1, TVXC)
          ENDIF

          IF (LADDER) THEN
             ! this line is a hook for the local_field routine
             ! in order to do some tests
             ! usually it simply returns doing nothing
             CALL DETERMINE_TBSE_FROM_FXC(CHI, TBSE, TVXC,  WGW, LATT_CUR, IO%IU0, IO%IU6)

             LINVXI=.NOT. LACFDT
             CALL CALCULATE_LOCAL_FIELD_FOCK( &
                  P,WDES,W_W,LATT_CUR,T_INFO,IO,KPOINTS, &
                  WGW, TBSE, TBSEA, CHI, ENCUTLF, ENCUTGW, ENCUTGWSOFT, COMEGA, & 
                  .NOT. LACFDT , & ! include direct contributions are included
                  LINVXI , & ! determine X-1 GG V GG X-1 (LINVXI=.FALSE. TBSE= GG V GG)
                  NBANDSO, NBANDSV, NKREDLFX, NKREDLFY, NKREDLFZ, NQ, ANTIRES, LGWLF, NELM)
          ENDIF
       ENDIF
       CALL STOP_TIMING("G",IO%IU6,"LFIELD")
!-----------------------------------------------------------------------
! calculate the screened potential
!-----------------------------------------------------------------------
       CALL DUMP_RESPONSE( CHI, WGWQ, NOMEGA, & 
                  "HEAD OF MICROSCOPIC DIELECTRIC TENSOR (INDEPENDENT PARTICLE)", IO%IU6, IO%NWRITE, INNER_LOOP_PARALLEL, .FALSE., & 
                  -EDEPS/LATT_CUR%OMEGA *LAMBDA, 1.0_q)
       IF (LCHIR) &
            CALL DUMP_RESPONSE( CHIR, WGWQ, NOMEGA, &
                  "HEAD OF MICROSCOPIC DIELECTRIC TENSOR (RESONANT)", IO%IU6, IO%NWRITE, INNER_LOOP_PARALLEL, .FALSE., & 
                  -EDEPS/LATT_CUR%OMEGA *LAMBDA, 1.0_q)
       CALL DUMP_RESPONSE_G_MIC( CHI, WGWQ, LATT_CUR, "epsilon_diag_mic")
       
       IF ((LOEP .OR. LFXCEPS) .AND. CHI%LGAMMA ) THEN 
          CALL SET_RESPONSE_KPOINT(CHI0, WDES%VKPT(:,NQ), NQ)
          IF ( CHI%NOMEGA_LOW==1) THEN
             CALL COPY_CHI( CHI, 1, CHI0, 1, 1.0_q)
             CHI0%COMEGA(1)=CHI%COMEGA(1)
          ENDIF
       ENDIF
!-----------------------------------------------------------------------
! adiabatic connection 
!-----------------------------------------------------------------------
       IF (LACFDT) THEN
#ifdef debug
   DO ISP = 1, CHI%NOMEGA
      IF ( CHI%LGAMMA )  THEN
         WRITE(100+WDES%COMM%NODE_ME,*)'HEAD',OMEGA(ISP)
         DO N1 = 1, 3
            WRITE(100+WDES%COMM%NODE_ME,'(3F12.6, "iamg:",3F12.6)')REAL(CHI%HEAD(N1,1:3,ISP),q),AIMAG(CHI%HEAD(N1,1:3,ISP)) 
         ENDDO
      ENDIF  
      WRITE(100+WDES%COMM%NODE_ME,*)'WINGS',OMEGA(ISP)
      DO N1 = 1, 16
         WRITE(100+WDES%COMM%NODE_ME,'(12F12.6)')CHI%WING(N1,1:3,ISP),CHI%CWING(N1,1:3,ISP)
      ENDDO
    
      WRITE(100+WDES%COMM%NODE_ME,*)'R',OMEGA(ISP)
      DO N1 = 1, 16
         WRITE(100+WDES%COMM%NODE_ME,'(16F12.6)')REAL(CHI%RESPONSEFUN(N1,1:16,ISP),q)
      ENDDO
      WRITE(100+WDES%COMM%NODE_ME,*)'I',OMEGA(ISP)
      DO N1 = 1, 16
         WRITE(100+WDES%COMM%NODE_ME,'(16F12.6)')AIMAG(CHI%RESPONSEFUN(N1,1:16,ISP))
      ENDDO
   ENDDO
#endif
          ! new fast memory conserving version
          CALL XI_ACFDT_ALL_RPA(IO%IU0, CHI, WGWQ, LATT_CUR, OMEGAWEIGHT, COR, NQ , IDIR_MAX, & 
             LRSRPA, INNER_LOOP_PARALLEL, LADDER .OR. LFXHEG, LINVXI, TBSE )
          
          ! symmetrize electronic structure factor    
          CALL SYMMETRIZE_ELSTUFAC( SYMM, LATT_CUR, T_INFO, DYN, INFO,&
               WDES, NQ, WGW, WGWQ, COR, LGAMMA, IO)

          CALL STOP_TIMING("G",IO%IU6,"ACFDT")

          IF (IO%IU6>=0) WRITE(IO%IU6,'(" q-point correlation energy ",2F14.6)') COR%CORRELATION_K(1)
          IF (IO%IU6>=0) WRITE(IO%IU6,'(" Hartree contr. to MP2      ",2F14.6)') COR%CORRMP2DIR_K(1)
          IF (LADDER .OR. LFXHEG) THEN
          IF (IO%IU6>=0) WRITE(IO%IU6,'(" q-point  SOSEX      energy ",2F14.6)') COR%CORRSOSEX_K(1)
          IF (IO%IU6>=0) WRITE(IO%IU6,'(" exchange cont. to MP2      ",2F14.6)') COR%CORRMP2EX_K(1)
          ENDIF
!-----------------------------------------------------------------------
! W from Xi (new version)
!-----------------------------------------------------------------------
       ELSE IF (.TRUE.) THEN
          ! calculate  X_red= X_0 (1- f_xc X_0 - v X_0)^-1
          IF (.NOT. L2ORDER) THEN
             IF (LCHIR) THEN
                CALL XI_LOCAL_FIELD_SYM( IO%IU0, CHI, CHIR, TVXC, TBSE, TBSEA, WGWQ, LHARTREE=.TRUE.)
             ELSE
                CALL XI_LOCAL_FIELD( IO%IU0, CHI, TVXC, TBSE, WGWQ, LHARTREE=.TRUE.)
             ENDIF
          ENDIF

          CALL STOP_TIMING("G",IO%IU6,"XI_LOCAL")

          CALL DUMP_RESPONSE( CHI, WGWQ, NOMEGA, & 
               "1 + v P,  with REDUCIBLE POLARIZABILTY P=P_0 (1 -(v+f) P_0)^-1",& 
               IO%IU6, IO%NWRITE, INNER_LOOP_PARALLEL, .TRUE., &
               EDEPS/LATT_CUR%OMEGA *LAMBDA, 1.0_q)

          IF (.NOT. LTETE) THEN
             !  eps^-1  = (1 + (f_xc +v) X_red )
             CALL XI_RED_TO_EPS( CHI, TVXC, TBSE, WGWQ, LTCTE)
             CALL STOP_TIMING("G",IO%IU6,"XI_RED_TO_EPS")
             IF (.NOT.LADDER .AND. .NOT. LFXC) THEN
                CALL DUMP_RESPONSE( CHI, WGWQ, NOMEGA, & 
                     "INVERSE MACROSCOPIC DIELECTRIC TENSOR (including local field effects in RPA (Hartree))",& 
                     IO%IU6, IO%NWRITE, INNER_LOOP_PARALLEL, .TRUE.)
             ELSE IF (LTCTE) THEN
                CALL DUMP_RESPONSE( CHI, WGWQ, NOMEGA, & 
                     "INVERSE MACROSCOPIC DIELECTRIC TENSOR (test charge-test electron, local field effects in DFT)",& 
                     IO%IU6, IO%NWRITE, INNER_LOOP_PARALLEL, .TRUE.)
             ELSE
                CALL DUMP_RESPONSE( CHI, WGWQ, NOMEGA, & 
                     "INVERSE MACROSCOPIC DIELECTRIC TENSOR (test charge-test charge, local field effects in DFT)",& 
                     IO%IU6, IO%NWRITE, INNER_LOOP_PARALLEL, .TRUE.)
             ENDIF
             CALL DUMP_RESPONSE_G( CHI, WGWQ, LATT_CUR, "epsilon_diag")
             CALL XI_TO_W( CHI, WGWQ)
          ELSE
             CALL XI_RED_TO_TETE( CHI, TVXC, TBSE, WGWQ, LTETE)
          ENDIF

          ! set fxc = epsilon^-1 X^-1
          IF (CHI%LGAMMA .AND. LFXCEPS ) THEN 
             CALL SET_RESPONSE_KPOINT(TVXC, WDES%VKPT(:,NQ), NQ)

             IF (.NOT. ASSOCIATED(TVXC%RESPONSEFUN)) THEN
                CALL ALLOCATE_RESPONSEFUN(TVXC, WGW%NGDIM, WGW%LGAMMA, .FALSE., 1)
             ENDIF
             TVXC%RESPONSEFUN=0

             IF ( CHI%NOMEGA_LOW==1) THEN
                CALL XI_FXC_FROM_EPS( IO%IU0, CHI, CHI0, TVXC, WGWQ )
             ENDIF
          ENDIF

          CALL STOP_TIMING("G",IO%IU6,"XI_TO_W")
          CALL XI_COULOMB( CHI, WGWQ, LATT_CUR, S2E%NUMBER_OF_NQ_FULL, FSG0, LFOCK_SUBTRACT)
       !
       ! old version using symmetric dielectric matrix (not used)
       !
       ELSE
          CALL XI_LOCAL_FIELD( IO%IU0, CHI, TVXC, TBSE, WGWQ, LHARTREE=.FALSE.)
          CALL STOP_TIMING("G",IO%IU6,"XI_LOCAL")

          CALL XI_TO_EPS( CHI, WGWQ)
          CALL STOP_TIMING("G",IO%IU6,"XI_TO_EPS")
          IF (.NOT.LADDER) THEN
             CALL DUMP_RESPONSE( CHI, WGWQ, NOMEGA, & 
                  "HEAD OF MICROSCOPIC DIELECTRIC TENSOR", IO%IU6, IO%NWRITE, INNER_LOOP_PARALLEL)
          ELSE
             CALL DUMP_RESPONSE( CHI, WGWQ, NOMEGA, & 
                  "HEAD OF MICROSCOPIC DIELECTRIC TENSOR (DFT)", IO%IU6, IO%NWRITE, INNER_LOOP_PARALLEL)
          ENDIF
          CALL XI_INVERT( IO%IU0, CHI, WGWQ )
          IF (.NOT.LADDER .AND. .NOT. LFXC) THEN
             CALL DUMP_RESPONSE( CHI, WGWQ, NOMEGA, & 
                  "INVERSE MACROSCOPIC DIELECTRIC TENSOR (including local field effects in RPA (Hartree))",& 
                  IO%IU6, IO%NWRITE, INNER_LOOP_PARALLEL, .TRUE.)
          ELSE
             CALL DUMP_RESPONSE( CHI, WGWQ, NOMEGA, & 
                  "INVERSE MACROSCOPIC DIELECTRIC TENSOR (including local field effects in DFT)",& 
                  IO%IU6, IO%NWRITE, INNER_LOOP_PARALLEL, .TRUE.)
          ENDIF
          CALL DUMP_RESPONSE_G( CHI, WGWQ, LATT_CUR, "epsilon_diag")
          
          CALL XI_TO_W_SYM( CHI, WGWQ)

          CALL STOP_TIMING("G",IO%IU6,"XI_TO_W_SYM")
          CALL XI_COULOMB( CHI, WGWQ, LATT_CUR, S2E%NUMBER_OF_NQ_FULL, FSG0, LFOCK_SUBTRACT)
       ENDIF

       ! symmetrize the final screened potential 
       ! not very elegant since we switch on symmetry
       ! calculate the required operations and switch off symmetry finally
       ! should be precalculated most likely and stored in an array KPOINTS_TRANS_WPOT
       IF (SYMM%ISYM>0 .AND. .NOT. LGAMMA ) THEN
          CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
               T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
               SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
               SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,-1)
          CALL KPOINTS_TRANS_Q(WGW%GRID, WGW, KPOINTS_TRANS_WPOT, NQ, SYMM%ROTMAP, SYMM%MAGROT )

          DO K1=1,CHI%NOMEGA
             CALL SYMMETRIZE_WPOT(WGWQ, KPOINTS_TRANS_WPOT, CHI%RESPONSEFUN(:,:,K1))
          ENDDO

          CALL DEALLOCATE_KPOINTS_TRANS(KPOINTS_TRANS_WPOT)
          
          ! switch of symmetry
          CALL NOSYMM(LATT_CUR%A,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS,SYMM%PTRANS, &
               SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%MAGROT,INFO%ISPIN,-1)
          CALL STOP_TIMING("G",IO%IU6,"CHI_SYM")
       ENDIF

! difficult to interpret so leave it out
       CALL DUMP_RESPONSE( CHI, WGWQ, NOMEGA, "screened Coulomb potential", IO%IU6, IO%NWRITE, INNER_LOOP_PARALLEL)
       IF (LGW .OR. LCRPA ) THEN
          CALL WRITE_WPOT(CHI, WGWQ, FSG0, NQ, LFOCK_SUBTRACT )
          CALL WRITE_WPOT_FULL(CHI, WGWQ, FSG0, NQ, LFOCK_SUBTRACT )
       ENDIF

       IF (LSPECTRAL .AND. LSPECTRALGW) THEN
          CALL W_SPECTRAL( CHI, WGWQ)
          CALL STOP_TIMING("G",IO%IU6,"W_SPECTRAL")
       ENDIF

       ENDIF calc_chi
!-----------------------------------------------------------------------
! orbitals optimized for COH, calculation of SCREENED_TWO_ELECTRON_INTEGRAL
! can be skipped for unoccupied states
!-----------------------------------------------------------------------
       IF (LGWNO) THEN
          DO ISP=1,WDES%ISPIN
             DO K1=1,WDES%NKPTS
                W%AUXTOT(LAST_FILLED_XI(W_W,K1,ISP,CHI%BETA)+1:, K1, ISP)=0
             ENDDO
          ENDDO
       ENDIF
!-----------------------------------------------------------------------
! calculate dynamically screened two electron integrals
!-----------------------------------------------------------------------
   gw: IF (LGW) THEN
       NOMEGA_INDEXI=0
       DO
       IF( SET_RESPONSEFUN_FREQI( WPOT, WGW, COMEGA, NOMEGA_INDEXI, &
            NOMEGA_INDEX, NOMEGA_INDEX+ SIZE(CHI%COMEGA)*NCPU-1, LSPECTRALCHI)) EXIT

       IF (.NOT.LSPECTRAL) THEN
        IF (.NOT. LscQPGW .AND. .FALSE.) THEN
!-----------------------------------------------------------------------
          IF (NCPU>1) THEN
             CALL MERGE_FREQU(CHI, WPOT, WGW, INNER_LOOP_PARALLEL)
          ENDIF
          ! loop over all k-points K1 (index a)
          DO ISP=1,WDES%ISPIN
          DO K1=1,WDES%NKPTS
#ifdef MPI   
             IF (INNER_LOOP_PARALLEL .AND. MOD(K1-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
             ! loop over all bands (index a) in blocks of NSTRIP
             N1_LAST=NBANDSGW/WDES%NB_PAR
             DO N1=1,N1_LAST,NSTRIP
                CALL GWPROGRESS(IO%IU0,K1, WDES%NKPTS, N1, N1_LAST)
                ! determine upper bound (avoid to go beyond the total
                ! number of bands
                NSTRIP1_ACT=MIN(NBANDSGW/WDES%NB_PAR+1-N1,NSTRIP)
                CALL SCREENED_TWO_ELECTRON_INTEGRAL( LMDIM, LATT_CUR, W, WGWQ, NQ,  &
                     H, P, ISP, S2E,  &
                     SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:,:,WPOT%NOMEGA_LOW:WPOT%NOMEGA_LOW+WPOT%NOMEGA-1),  &
                     K1, N1, NSTRIP1_ACT, NSTRIP_TOTAL, WPOT, NOPER, NFLOAT)
             ENDDO
          ENDDO
          ENDDO
!$OMP PARALLEL
          CALL CLEAN_RESPONSEFUNCTION_INT( WGWQ, WPOT, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, WCACHE, 1, WPOT%NOMEGA )
!$OMP END PARALLEL
          CALL STOP_TIMING("G",IO%IU6,"TWOEINT")
          CALL GWPROGRESS(IO%IU0, WDES%NKPTS, WDES%NKPTS, N1_LAST, N1_LAST)
        ELSE
!-----------------------------------------------------------------------
          ! same version but using  CACHED structure
          IF (NCPU>1) THEN
             CALL MERGE_FREQU(CHI, WPOT, WGW, INNER_LOOP_PARALLEL)
          ENDIF

          ! loop over all k-points K1 (index a)
          DO ISP=1,WDES%ISPIN
          DO K1=1,WDES%NKPTS
#ifdef MPI   
             IF (INNER_LOOP_PARALLEL .AND. MOD(K1-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
             ! loop over all bands (index a) in blocks of NSTRIP
             N1_LAST=NBANDSGW/WDES%NB_PAR
             DO N1=1,N1_LAST,NSTRIP
                CALL GWPROGRESS(IO%IU0,K1, WDES%NKPTS, N1, N1_LAST)
                ! determine upper bound (avoid to go beyond the total
                ! number of bands
                NSTRIP1_ACT=MIN(NBANDSGW/WDES%NB_PAR+1-N1,NSTRIP)
                IF (SAVE_CACHE_MEMORY) CALL NBLOW_CACHER( WCACHE, N1, N1-1+NSTRIP1_ACT, WDES%NB_PAR )
                ! alternative version that uses the cache structure in WPOT
                CALL SCREENED_TWO_ELECTRON_CACHED( LMDIM, LATT_CUR, W, WGWQ, NQ, &
                     H, P, ISP, S2E, WCACHE, &
                     SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:,:,WPOT%NOMEGA_LOW:WPOT%NOMEGA_LOW+WPOT%NOMEGA-1),  &
                     K1, N1, NSTRIP1_ACT, WPOT, OMEGA, NOPER, NFLOAT, 0)
                IF (LscQPGW .AND. SAVE_CACHE_MEMORY) THEN
!$OMP PARALLEL
                   CALL CLEAN_RESPONSEFUNCTION_INT(WGWQ, WPOT, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, WCACHE, 1, WPOT%NOMEGA )
!$OMP END PARALLEL
                   CALL STORE_GW_ACC_FINAL( WCACHE, 1,  WACC1, LATT_CUR, NONLR_S, NONL_S, K1, ISP, W)
                   CALL STORE_GW_ACC_FINAL( WCACHE, 2,  WACC2, LATT_CUR, NONLR_S, NONL_S, K1, ISP, W)
                ENDIF
             ENDDO
             IF (LscQPGW .AND. .NOT. SAVE_CACHE_MEMORY) THEN
!$OMP PARALLEL
                CALL CLEAN_RESPONSEFUNCTION_INT(WGWQ, WPOT, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, WCACHE, 1, WPOT%NOMEGA )
!$OMP END PARALLEL
                CALL STORE_GW_ACC_FINAL( WCACHE, 1,  WACC1, LATT_CUR, NONLR_S, NONL_S, K1, ISP, W)
                CALL STORE_GW_ACC_FINAL( WCACHE, 2,  WACC2, LATT_CUR, NONLR_S, NONL_S, K1, ISP, W)
             ENDIF
          ENDDO
          ENDDO
!$OMP PARALLEL
          CALL CLEAN_RESPONSEFUNCTION_INT( WGWQ, WPOT, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, WCACHE, 1, WPOT%NOMEGA )
!$OMP END PARALLEL
          CALL STOP_TIMING("G",IO%IU6,"TWOEINT")
          CALL GWPROGRESS(IO%IU0, WDES%NKPTS, WDES%NKPTS, N1_LAST, N1_LAST)

        ENDIF
!-----------------------------------------------------------------------
     ELSE
          IF (LSELFENERGY) THEN
             ! same as above but now use the spectral representation of W
             ! for positive shifts i delta
             CALL DO_POT_HILBERT_TABLE(CHI, WPOT, WGW, TABLE_POT_PLUS, INNER_LOOP_PARALLEL)
             CALL STOP_TIMING("G",IO%IU6,"KRAMKRO")

             DO ISP=1,WDES%ISPIN
             DO K1=1,WDES%NKPTS
#ifdef MPI   
                IF (INNER_LOOP_PARALLEL .AND. MOD(K1-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
                N1_LAST=NBANDSGW/WDES%NB_PAR
                DO N1=1,N1_LAST,NSTRIP
                   CALL GWPROGRESS(IO%IU0, K1, WDES%NKPTS, N1, N1_LAST)
                   NSTRIP1_ACT=MIN(NBANDSGW/WDES%NB_PAR+1-N1,NSTRIP)
                   CALL SCREENED_TWO_ELECTRON_CACHED( LMDIM, LATT_CUR, W, WGWQ, NQ,  &
                        H, P, ISP, S2E, WCACHE, &
                        SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:,:,NOMEGA_INDEX:NOMEGA_INDEX+WPOT%NOMEGA-1),  &
                        K1, N1, NSTRIP1_ACT, WPOT, OMEGA, NOPER, NFLOAT, 0 )
                ENDDO
             ENDDO
             ENDDO
!$OMP PARALLEL
             CALL CLEAN_RESPONSEFUNCTION_INT( WGWQ, WPOT, S2E, & 
                  SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:,:,NOMEGA_INDEX:NOMEGA_INDEX+WPOT%NOMEGA-1), & 
                  WCACHE, 1, WPOT%NOMEGA )
!$OMP END PARALLEL
             CALL STOP_TIMING("G",IO%IU6,"TWOEINT")
             CALL GWPROGRESS(IO%IU0, WDES%NKPTS, WDES%NKPTS, N1_LAST, N1_LAST)

             ! now for negative shifts - i delta
             CALL DO_POT_HILBERT_TABLE(CHI, WPOT, WGW, TABLE_POT_MIN, INNER_LOOP_PARALLEL)
             CALL STOP_TIMING("G",IO%IU6,"KRAMKRO")

             DO ISP=1,WDES%ISPIN
             DO K1=1,WDES%NKPTS
#ifdef MPI   
                IF (INNER_LOOP_PARALLEL .AND. MOD(K1-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
                N1_LAST=NBANDSGW/WDES%NB_PAR
                DO N1=1,N1_LAST,NSTRIP
                   CALL GWPROGRESS(IO%IU0, K1, WDES%NKPTS, N1, N1_LAST)
                   NSTRIP1_ACT=MIN(NBANDSGW/WDES%NB_PAR+1-N1,NSTRIP)
                   CALL SCREENED_TWO_ELECTRON_CACHED( LMDIM, LATT_CUR, W, WGWQ, NQ,  &
                        H, P, ISP, S2E, WCACHE, &
                        SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:,:,NOMEGA_INDEX+NOMEGA:NOMEGA_INDEX+WPOT%NOMEGA-1+NOMEGA),  &
                        K1, N1, NSTRIP1_ACT, WPOT, OMEGA, NOPER, NFLOAT, 0 )
                ENDDO
             ENDDO
             ENDDO
!$OMP PARALLEL
             CALL CLEAN_RESPONSEFUNCTION_INT(WGWQ, WPOT, S2E, & 
                  SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:,:,NOMEGA_INDEX+NOMEGA:NOMEGA_INDEX+WPOT%NOMEGA-1+NOMEGA), &
                  WCACHE, 1, WPOT%NOMEGA )
!$OMP END PARALLEL
             CALL STOP_TIMING("G",IO%IU6,"TWOEINT")
             CALL GWPROGRESS(IO%IU0, WDES%NKPTS, WDES%NKPTS, N1_LAST, N1_LAST)

!-----------------------------------------------------------------------
          ELSE
             ! only small energy interval around eigenvalues
             ! this is done by supplying the sign of the supplied WPOT as the last argument
             ! to SCREENED_TWO_ELECTRON_CACHED
             CALL DO_POT_HILBERT_TABLE(CHI, WPOT, WGW, TABLE_POT_PLUS, INNER_LOOP_PARALLEL)
             CALL STOP_TIMING("G",IO%IU6,"KRAMKRO")

             DO ISP=1,WDES%ISPIN
             DO K1=1,WDES%NKPTS
#ifdef MPI   
                IF (INNER_LOOP_PARALLEL .AND. MOD(K1-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
                N1_LAST=NBANDSGW/WDES%NB_PAR
                DO N1=1,N1_LAST,NSTRIP
                   CALL GWPROGRESS(IO%IU0, K1, WDES%NKPTS, N1, N1_LAST)
                   NSTRIP1_ACT=MIN(NBANDSGW/WDES%NB_PAR+1-N1,NSTRIP)
                   IF (SAVE_CACHE_MEMORY) CALL NBLOW_CACHER( WCACHE, N1, N1-1+NSTRIP1_ACT, WDES%NB_PAR )
                   CALL SCREENED_TWO_ELECTRON_CACHED( LMDIM, LATT_CUR, W, WGWQ, NQ, &
                        H, P, ISP, S2E, WCACHE, SCREENED_TWO_ELECTRON_INTEGRALS,  &
                        K1, N1, NSTRIP1_ACT, WPOT, OMEGA, NOPER, NFLOAT, 1 )
                   IF (LscQPGW .AND. SAVE_CACHE_MEMORY) THEN
!$OMP PARALLEL
                      CALL CLEAN_RESPONSEFUNCTION_INT(WGWQ, WPOT, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, WCACHE, 1, WPOT%NOMEGA )
!$OMP END PARALLEL
                      CALL STORE_GW_ACC_FINAL( WCACHE, 1,  WACC1, LATT_CUR, NONLR_S, NONL_S, K1, ISP, W)
                      CALL STORE_GW_ACC_FINAL( WCACHE, 2,  WACC2, LATT_CUR, NONLR_S, NONL_S, K1, ISP, W)
                   ENDIF
                ENDDO
                IF (LscQPGW .AND. .NOT. SAVE_CACHE_MEMORY) THEN
!$OMP PARALLEL
                   CALL CLEAN_RESPONSEFUNCTION_INT(WGWQ, WPOT, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, WCACHE, 1, WPOT%NOMEGA )
!$OMP END PARALLEL
                   CALL STORE_GW_ACC_FINAL( WCACHE, 1,  WACC1, LATT_CUR, NONLR_S, NONL_S, K1, ISP, W)
                   CALL STORE_GW_ACC_FINAL( WCACHE, 2,  WACC2, LATT_CUR, NONLR_S, NONL_S, K1, ISP, W)
                ENDIF
             ENDDO
             ENDDO
             ! calculate all required data from the chached results
!$OMP PARALLEL
             CALL CLEAN_RESPONSEFUNCTION_INT(WGWQ, WPOT, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, WCACHE, 1, WPOT%NOMEGA )
!$OMP END PARALLEL
             CALL STOP_TIMING("G",IO%IU6,"TWOEINT")
             CALL GWPROGRESS(IO%IU0, WDES%NKPTS, WDES%NKPTS, N1_LAST, N1_LAST)

             ! now for negative shifts - i delta
             CALL DO_POT_HILBERT_TABLE(CHI, WPOT, WGW, TABLE_POT_MIN, INNER_LOOP_PARALLEL)
             CALL STOP_TIMING("G",IO%IU6,"KRAMKRO")

             DO ISP=1,WDES%ISPIN
             DO K1=1,WDES%NKPTS
#ifdef MPI   
                IF (INNER_LOOP_PARALLEL .AND. MOD(K1-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
                N1_LAST=NBANDSGW/WDES%NB_PAR
                DO N1=1,N1_LAST,NSTRIP
                   CALL GWPROGRESS(IO%IU0, K1, WDES%NKPTS, N1, N1_LAST)
                   NSTRIP1_ACT=MIN(NBANDSGW/WDES%NB_PAR+1-N1,NSTRIP)

                   IF (SAVE_CACHE_MEMORY) CALL NBLOW_CACHER( WCACHE, N1, N1-1+NSTRIP1_ACT, WDES%NB_PAR )

                   CALL SCREENED_TWO_ELECTRON_CACHED( LMDIM, LATT_CUR, W, WGWQ, NQ,  &
                        H, P, ISP, S2E, WCACHE, SCREENED_TWO_ELECTRON_INTEGRALS,  &
                        K1, N1, NSTRIP1_ACT, WPOT, OMEGA, NOPER, NFLOAT, -1 )

                   IF (LscQPGW .AND. SAVE_CACHE_MEMORY) THEN
!$OMP PARALLEL
                      CALL CLEAN_RESPONSEFUNCTION_INT(WGWQ, WPOT, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, WCACHE, 1, WPOT%NOMEGA )
!$OMP END PARALLEL
                      CALL STORE_GW_ACC_FINAL( WCACHE, 1,  WACC1, LATT_CUR, NONLR_S, NONL_S, K1, ISP, W)
                      CALL STORE_GW_ACC_FINAL( WCACHE, 2,  WACC2, LATT_CUR, NONLR_S, NONL_S, K1, ISP, W)
                   ENDIF
                ENDDO
                IF (LscQPGW .AND. .NOT. SAVE_CACHE_MEMORY) THEN
!$OMP PARALLEL
                   CALL CLEAN_RESPONSEFUNCTION_INT(WGWQ, WPOT, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, WCACHE, 1, WPOT%NOMEGA )
!$OMP END PARALLEL
                   CALL STORE_GW_ACC_FINAL( WCACHE, 1,  WACC1, LATT_CUR, NONLR_S, NONL_S, K1, ISP, W)
                   CALL STORE_GW_ACC_FINAL( WCACHE, 2,  WACC2, LATT_CUR, NONLR_S, NONL_S, K1, ISP, W)
                ENDIF
             ENDDO
             ENDDO
!$OMP PARALLEL
             CALL CLEAN_RESPONSEFUNCTION_INT(WGWQ, WPOT, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, WCACHE, 1, WPOT%NOMEGA )
!$OMP END PARALLEL
             CALL STOP_TIMING("G",IO%IU6,"TWOEINT")
             CALL GWPROGRESS(IO%IU0, WDES%NKPTS, WDES%NKPTS, N1_LAST, N1_LAST)

          ENDIF
       ENDIF

       ENDDO
       END IF gw

!       IF (LGWNO) W%AUXTOT=1

       CALL ISEMTPY_CACHER(WCACHE)

       CALLMPI( M_sum_d(WDES%COMM, NFLOAT, 1))
       CALLMPI( M_sum_i(WDES%COMM, NOPER , 1))

       NOPER_=NOPER_+NOPER/1000
       NFLOAT_=NFLOAT_+NFLOAT

       CALL STOP_TIMING("GW",IO%IU6)

       IF (IO%IU0>=0) WRITE(IO%IU0,*) 
       IF (IO%IU0>=0) WRITE(IO%IU0,10) NOPER_, 1E-9_q*NFLOAT_

       IF (IO%IU0>=0) WRITE(17,10) NOPER_, 1E-9_q*NFLOAT_

       IF (IO%IU6>=0) WRITE(IO%IU6,*) 
       IF (IO%IU6>=0) WRITE(IO%IU6,10) NOPER_, 1E-9_q*NFLOAT_
10     FORMAT(" performed ",I10,"000 updates of chi_q(r,r)",/ &
            " total number of BLAS operations",F12.2," Gflops")

    ENDDO lchi

    IF (IO%LOPEN) CALL WFORCE(IO%IU6)
 ENDDO qpoints

    IF (LACFDT .AND. OUTER_LOOP_PARALLEL) then
      CALLMPI( M_sum_z( WDES%COMM_KINTER, COR%CORRELATION, SIZE(COR%CORRELATION)))
      CALLMPI( M_sum_z( WDES%COMM_KINTER, COR%CORRMP2DIR,  SIZE(COR%CORRMP2DIR)))
      CALLMPI( M_sum_z( WDES%COMM_KINTER, COR%CORRELATION, SIZE(COR%CORRSOSEX)))
      CALLMPI( M_sum_z( WDES%COMM_KINTER, COR%CORRMP2DIR,  SIZE(COR%CORRMP2EX)))
    ENDIF

    IF (LSPECTRAL .OR. NCPU>1) THEN
       CALL DEALLOCATE_RESPONSEFUN( WPOT )
       IF (LSPECTRAL) THEN
          CALL DEALLOCATE_KRAM_KRON_TABLE( TABLE)
          CALL DEALLOCATE_KRAM_KRON_TABLE( TABLE_POT_PLUS)
          CALL DEALLOCATE_KRAM_KRON_TABLE( TABLE_POT_MIN)
       ENDIF
    ENDIF
    IF (LCHIR) THEN
       CALL DEALLOCATE_KRAM_KRON_TABLE( TABLE_RES)
       CALL DEALLOCATE_RESPONSEFUN( CHIR )
    ENDIF
    
    IF (LCHIREALLOCATE) THEN
       CALL DEALLOCATE_RESPONSEFUN( CHI )
    ENDIF
    IF (LADDER .AND. NBANDSO>0) CALL DEALLOCATE_RESPONSEFUN( TBSE )
    IF (LADDER .AND. NBANDSO>0 .AND. ANTIRES>=2)  CALL DEALLOCATE_RESPONSEFUN( TBSEA )
!-----------------------------------------------------------------------
! resolve degeneracies in screened two electron integrals
!-----------------------------------------------------------------------
    IF (LGW) THEN
    IF (IO%IU0>=0) WRITE(IO%IU0,*) 'resolving degeneracies of screened two electron integrals'
    IF (IO%IU0>=0) WRITE(17,*) 'resolving degeneracies of screened two electron integrals'
    DO NOMEGA_INDEX=1, SIZE(SCREENED_TWO_ELECTRON_INTEGRALS, 6)
       CALL CLEANUP_SCREENED_2E(W, S2E, SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:,:,NOMEGA_INDEX), SYMM%ISYM, INNER_LOOP_PARALLEL)
    ENDDO

#ifdef MPI
!    CALL WRITE_SCREENED_2E_FILE(SCREENED_TWO_ELECTRON_INTEGRALS, IO%IU0, NOMEGA, NOMEGAR, OMEGAMAX, OMEGATL, WDES%COMM_INTER%NODE_ME)
#else
!    CALL WRITE_SCREENED_2E_FILE(SCREENED_TWO_ELECTRON_INTEGRALS, IO%IU0, NOMEGA, NOMEGAR, OMEGAMAX, OMEGATL)
#endif

     ENDIF
!-----------------------------------------------------------------------
! calculate self energy and QP shifts
!-----------------------------------------------------------------------
    IF (LGW) THEN
    IF (IO%IU6>=0) WRITE(IO%IU6,130)

    IF (LSELFENERGY .AND. LSPECTRAL) THEN
       ! spectral representation of W is available
       CALL CALC_SELFENERGY_LINEAR(W, S2E, & 
            SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:,:,:), &
            SHIFT, OMEGA, CELTOT_X+CELTOT_HARTREE_KINETIC, LFOCK_SUBTRACT, IO%IU6, IO%IU0, EFERMI )
    ELSE IF (LSELFENERGY) THEN
       IF (IO%IU0>=0) WRITE(IO%IU0,*) 'calculate frequency dependent selfenergy'
       IF (IO%IU0>=0) WRITE(17,*) 'calculate frequency dependent selfenergy'

       ! analytical continuation
       IF (NOMEGAR==0) THEN
          IF (ASSOCIATED(IMAG_GRIDS%TO_FER_RE) .AND.  ASSOCIATED(IMAG_GRIDS%TO_FER_IM)) THEN
             CALL CALC_SELFENERGY_IMAG_PADE(W, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, & 
                  OMEGA, OMEGAWEIGHT, NOMEGAR, CELTOT_X+CELTOT_HARTREE_KINETIC, &
                  LFOCK_SUBTRACT, IO%IU5, IO%IU6, IO%IU0, EFERMI, IMAG_GRIDS)
          ELSE
             IF ( OMEGAGRID==46 )THEN
                CALL CALC_SELFENERGY_IMAG_GAUSS(W, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, & 
                 SHIFT, OMEGA, OMEGAWEIGHT, NOMEGAR, CELTOT_X+CELTOT_HARTREE_KINETIC, LFOCK_SUBTRACT, IO%IU6, IO%IU0,EFERMI )
              ELSE
                 CALL CALC_SELFENERGY_IMAG(W, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, & 
                 SHIFT, OMEGA, NOMEGAR, CELTOT_X+CELTOT_HARTREE_KINETIC, LFOCK_SUBTRACT, IO%IU6, IO%IU0, EFERMI )
              ENDIF
           ENDIF
       ! complex contour integrals
       ELSEIF (NOMEGAR<NOMEGA) THEN
          CALL CALC_SELFENERGY_C(W, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, & 
               SHIFT, OMEGA, OMEGAWEIGHT, NOMEGAR, CELTOT_X+CELTOT_HARTREE_KINETIC, LFOCK_SUBTRACT, IO%IU6, IO%IU0, EFERMI )
       ENDIF

       ! real axis integration
       IF (NOMEGAR>0) THEN
          CALL CALC_SELFENERGY(W, S2E, SCREENED_TWO_ELECTRON_INTEGRALS(:,:,:,:,:,1:NOMEGAR), & 
               SHIFT, OMEGA(1:NOMEGAR), CELTOT_X+CELTOT_HARTREE_KINETIC, LFOCK_SUBTRACT, IO%IU6, IO%IU0, EFERMI )
       ENDIF
    ELSE IF (ASSOCIATED(CELTOT_X) .AND. (.NOT. LscQPGW .OR. NELM==1 )) THEN
       IF (NOMEGA>2 .AND. NOMEGAR<=1) THEN
          CALL QP_SHIFT_PADE_SIGMA(W, LSPECTRAL, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, &
               SHIFT, OMEGA, OMEGAWEIGHT, NOMEGAR, CELTOT_HARTREE_KINETIC, CELTOT_X, & 
               LFOCK_SUBTRACT, NELM, LscQPGW, LG0W0, EFERMI, IMAG_GRIDS, IO%IU6, IO%IU5, IO%IU0)
       ELSE
          CALL QP_SHIFT(W, LSPECTRAL, S2E, SCREENED_TWO_ELECTRON_INTEGRALS, &
               SHIFT, OMEGA, OMEGAWEIGHT, NOMEGAR, CELTOT_HARTREE_KINETIC, CELTOT_X, & 
               LFOCK_SUBTRACT, NELM, LscQPGW, LG0W0, IO%IU6, IO%IU0)
       ENDIF

       IF (IO%IU6>=0) WRITE(IO%IU6,130)
       CALL STOP_TIMING("G",IO%IU6,"QP_SHIFT")
    ENDIF
    ELSE IF (LACFDT) THEN
       IF ( LESF_SPLINES )  THEN 
          CALL EVALUATE_ESF_FOR_RPA_ENERGY( W%WDES, COR, LATT_CUR, INNER_LOOP_PARALLEL, IO ) 
       ENDIF

       CALL LIN_REG(COR, IO%IU6)
       IF (LADDER .OR. LFXHEG) THEN
          CALL LIN_REG_EX(COR, IO%IU6)
       ENDIF
       ! in case electronic structure factor has been collected
       ! deallocate electronic structure factor
       CALL DEALLOCATE_ELSTUFAC_HANDLE( COR ) 
    ENDIF


    IF (SYMM%ISYM>0) CALL CLEANUP_CELEN(W)     ! set degenerated or near degenerated eigenvalues to a single value
!-----------------------------------------------------------------------
! SCGW aka Kotani et al.
!-----------------------------------------------------------------------
    IF (LscQPGW .AND. .NOT. LOEP) THEN
       CALL DEALLOCATE_CACHER(WCACHE )

       CALLMPI( M_sum_z(WDES%COMM_KINTER,WACC1%CW,SIZE(WACC1%CW)))
       CALLMPI( M_sum_g(WDES%COMM_KINTER,WACC1%CPROJ,SIZE(WACC1%CPROJ)))
       CALLMPI( M_sum_z(WDES%COMM_KINTER,WACC2%CW,SIZE(WACC2%CW)))
       CALLMPI( M_sum_g(WDES%COMM_KINTER,WACC2%CPROJ,SIZE(WACC2%CPROJ)))

       ! bring the contributions from k-point in the full BZ into IRZ
       CALL CONTRACT_WAVE(WACC1, GRID, WDES )
       CALL CONTRACT_WAVE(WACC2, GRID, WDES )

       ! switch to IRZ
       IF (SYMM%ISYM>=0 .AND. .NOT. LGAMMA) THEN
          IF (SYMM%ISYM>0) &
          CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
               T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
               SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
               SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,-1)
          CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
               SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
               T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,-1,IO%IU0)
          CALL KPAR_SYNC_ALL(WDES,W)
          CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
          CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
          CALL REALLOCATE_WAVE( WACC1, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
          CALL REALLOCATE_WAVE( WACC2, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
       ENDIF

       ! use HF exchange only
       IF (.NOT. LGWLF) THEN
          CALL PUSH_XC_TYPE_FOR_GW
          IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
             CALL SETUP_LDA_XC(2,-1,-1,IO%IDIOT)
          ELSE
             CALL SETUP_LDA_XC(1,-1,-1,IO%IDIOT)
          ENDIF
          CALL SET_FSG_STORE( GRIDHF, LATT_CUR, WDES)
!          DO K1=1,WDES%NKPTS
!             FSG_STORE(K1)=SET_FSG(GRIDHF, LATT_CUR, K1)
!          ENDDO
       ENDIF

       ! update charge
       CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
            GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
            LATT_CUR, P, SYMM, T_INFO, &
            CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

       ! mix charge density with the one of the previous step
       IF (MIX%IMIX/=0) INFO%TIME=1.0
       IF (NELM>=2 .AND. MIX%IMIX/=0 ) THEN
          RMST=0
          CALL MIX_SIMPLE(GRIDC,MIX,WDES%NCDIJ, CHTOT,CHTOTL, &
               N_MIX_PAW, RHOLM, RHOLM_LAST, LATT_CUR%B, LATT_CUR%OMEGA, RMST)
          IF (IO%IU6>=0) WRITE(IO%IU6,'(" charge density residual (rmsc) ",E14.5)') RMST
       ENDIF

       ! store old density
       DO ISP=1,WDES%NCDIJ
          CALL RC_ADD(CHTOT(1,ISP),1.0_q,CHTOT(1,ISP),0.0_q,CHTOTL(1,ISP),GRIDC)
       ENDDO
       RHOLM_LAST=RHOLM

       !  calculate potential
       CALL POTLOK(GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES,  &
            INFO,P,T_INFO,E,LATT_CUR, &
            CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF)

       ! add the one center augmentation related terms  
       CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
            LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

       ! finally add one center terms
       CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
            WDES%NCDIJ, LMDIM, CDIJ(1,1,1,1),  RHOLM, CRHODE(1,1,1,1), &
            E,  LMETA=.FALSE., LASPH=INFO%LASPH, LCOREL= .FALSE.  )

       IF (IO%IU6>=0) WRITE(IO%IU6,130)
       CALL EDDIAG_GW( NBANDSGW, W%WDES, LATT_CUR,  NONLR_S, NONL_S, W, WACC1, WACC2,  &
            LMDIM, CDIJ, CQIJ, SV, T_INFO, P, NELM, INFO%TIME, SYMM, .NOT. LSPECTRAL, LGWNO, &
            KPOINTS, INFO%WEIMIN, EFERMI, IO%IU0, IO%IU6 )
       IF (IO%IU6>=0) WRITE(IO%IU6,130)

       ! force orbitals to be real, at those k-points where this is possible
       ! unfortunately, this seems to screw up the code right now
       ! no idea why
       ! obviously the scGW makes the orbitals somehow complex
       ! maybe this has to do with the violation of the shell structure
       ! or it is the derivatives of the orbitals
!       IF (W%WDES%LORBITALREAL .AND. INFO%LDIAG ) THEN
!          CALL WVREAL_PRECISE(W)
!          CALL PROALL (GRID,LATT_CUR,NONLR_S,NONL_S,W)
!          CALL STOP_TIMING("G",IO%IU6,"WVREAL",XMLTAG="wvreal")
!          CALL ORTHCH(WDES,W, INFO%LOVERL, LMDIM,CQIJ)
!          CALL STOP_TIMING("G",IO%IU6,"ORTHCH",XMLTAG="orth")
!       ENDIF


       IF (SYMM%ISYM>0) CALL CLEANUP_CELEN(W) ! set degenerated or near degenerated eigenvalues to a single value

       IF (.NOT. LGWLF) THEN
          CALL POP_XC_TYPE
          IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
             CALL SETUP_LDA_XC(2,-1,-1,IO%IDIOT)
          ELSE
             CALL SETUP_LDA_XC(1,-1,-1,IO%IDIOT)
          ENDIF
       ENDIF

       ! deallocate all
       CALL DEALLOCW( WACC1)
       CALL DEALLOCW( WACC2)

       ! update occupancies (done below, but we might need EFERMI in the optics routine)
       !keep fermi energy as it is for finite temperature calculations 
       IF ( .NOT. LFINITE_TEMPERATURE ) &
       CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
         INFO%NUP_DOWN,  E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
         NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)


       IF (.NOT. LGW0 .AND.  IO%LOPTICS) THEN
          CALL START_TIMING("G")
         ! VASP onboard optics
          CALL PEAD_RESETUP_WDES(WDES,GRID,KPOINTS,LATT_CUR,LATT_INI,IO)

          CALL LR_OPTIC( &
             P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
             T_INFO,INFO,IO,KPOINTS,SYMM,GRID,GRID_SOFT, &
             GRIDC,GRIDUS,C_TO_US,SOFT_TO_C, &
             CHTOT,DENCOR,CVTOT,CSTRF, &
             CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM, &
             CHDEN,SV,LMDIM,IRDMAX,EFERMI,NEDOS, & 
             LSTORE=.TRUE., LPOT=.FALSE.)
          CALL STOP_TIMING("G",IO%IU6,'OPTICS')
        ENDIF

       ! restore full k-point grid
       IF (SYMM%ISYM>=0 .AND. .NOT. LGAMMA) THEN
          CALL NOSYMM(LATT_CUR%A,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS,SYMM%PTRANS, &
               &   SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%MAGROT,INFO%ISPIN,IO%IU6)

          CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR,.FALSE., &
               T_INFO%NIONS,SYMM%ROTMAP, SYMM%MAGROT, SYMM%ISYM,-1,IO%IU0)
          CALL KPAR_SYNC_ALL(WDES,W)
          CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
          CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
       ENDIF
       IF (.NOT. LGW0) THEN
          W_W=W
       ENDIF
       CALL STOP_TIMING("G",IO%IU6,"GW_DIAG")
    ENDIF
!-----------------------------------------------------------------------
! OEP method
!-----------------------------------------------------------------------
    IF (LOEP) THEN
       CALL DEALLOCATE_CACHER(WCACHE )

       CALLMPI( M_sum_z(WDES%COMM_KINTER,WACC1%CW,SIZE(WACC1%CW)))
       CALLMPI( M_sum_g(WDES%COMM_KINTER,WACC1%CPROJ,SIZE(WACC1%CPROJ)))
       CALLMPI( M_sum_z(WDES%COMM_KINTER,WACC2%CW,SIZE(WACC2%CW)))
       CALLMPI( M_sum_g(WDES%COMM_KINTER,WACC2%CPROJ,SIZE(WACC2%CPROJ)))

       CALL CONTRACT_WAVE(WACC1, GRID, WDES )
       CALL CONTRACT_WAVE(WACC2, GRID, WDES )

       ! switch to IRZ
       IF (SYMM%ISYM>=0 .AND. .NOT. LGAMMA) THEN
          IF (SYMM%ISYM>0) &
          CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
               T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
               SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
               SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,-1)
          CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
               SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
               T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,-1,IO%IU0)
          CALL KPAR_SYNC_ALL(WDES,W)
          CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
          CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
          CALL REALLOCATE_WAVE( WACC1, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
          CALL REALLOCATE_WAVE( WACC2, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
       ENDIF

       ! use HF exchange only
       ! correlation part of DFT-functional is set to 0 
       IF (.NOT. LGWLF) THEN
          CALL PUSH_XC_TYPE_FOR_GW
          IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
             CALL SETUP_LDA_XC(2,-1,-1,IO%IDIOT)
          ELSE
             CALL SETUP_LDA_XC(1,-1,-1,IO%IDIOT)
          ENDIF
          CALL SET_FSG_STORE( GRIDHF, LATT_CUR, WDES)
!          DO K1=1,WDES%NKPTS
!             FSG_STORE(K1)=SET_FSG(GRIDHF, LATT_CUR, K1)
!          ENDDO
       ENDIF

       ! update charge
       CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
            GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
            LATT_CUR, P, SYMM, T_INFO, &
            CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

       IF (.NOT. INFO%LPOTOK) THEN
          CALL vtutor%write(isError, LPOTOK)
          ! calculate potential here
       ENDIF

       IF (SYMM%ISYM>=0 .AND. .NOT. LGAMMA) THEN
         CALL APPLY_SMALL_SPACE_GROUP_OP( W, WACC1, NONLR_S, NONL_S, & 
            P, T_INFO%NIONS, LATT_CUR, SYMM, CQIJ, .FALSE. , IO%IU6)
     
         CALL APPLY_SMALL_SPACE_GROUP_OP( W, WACC2, NONLR_S, NONL_S, & 
            P, T_INFO%NIONS, LATT_CUR, SYMM, CQIJ, .FALSE. , -1)

          CALL STOP_TIMING("G",IO%IU6,"GWSYM")
       ENDIF

       IF (IO%IU6>=0) WRITE(IO%IU6,130)

       CALL ALLOCATE_RESPONSEFUN_CACHE( CHI0, MAXCACHE)
       CALL SETWDES(WGW, WGWQ, CHI0%NQ )

       CALL OEP_GW( NBANDSGW, W%WDES, LATT_CUR,  NONLR_S, NONL_S, W, WACC1, WACC2, LEXX, &
            LMDIM, CDIJ, CQIJ, SV, T_INFO, P, INFO%TIME, LGW0, IO%IU0, IO%IU6, &
            KPOINTS%ISMEAR, KPOINTS%SIGMA, SYMM, .NOT. LSPECTRAL , &
            INFO, WGWQ, CHI0, GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
            CHTOT, DENCOR, CVTOT, CSTRF, IRDMAX, CRHODE, N_MIX_PAW, MIX%AMIX, RHOLM )
       IF (IO%IU6>=0) WRITE(IO%IU6,130)

       CALL STOP_TIMING("G",IO%IU6,"OEPGW")

       ! Davidson without HF contribution
        NSIM=WDES%NSIM
#ifdef MPI
        NSIM=((WDES%NSIM+WDES%COMM_INTER%NCPU-1)/WDES%COMM_INTER%NCPU)*WDES%COMM_INTER%NCPU
#endif
        
        CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
             LMDIM,CDIJ,CQIJ, 3,SV,T_INFO,P,IO%IU0,E%EXHF)
        IF (IO%IU0>0) WRITE(IO%IU0,*) 'calling exact diagonalization EDDIAG_EXACT'
!        CALL EDDIAG_EXACT(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES, &
!             LMDIM,CDIJ,CQIJ, 3,SV,T_INFO,P,IO%IU0,E%EXHF)
!        CALL EDDAV(HAMILTONIAN,P, GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WDES, NSIM, &
!             LMDIM,CDIJ,CQIJ, RMS,DESUM1,ICOUEV, SV,E%EXHF, IO%IU6,IO%IU0, .FALSE., .TRUE., .FALSE.)


       ! convergence reached on eigenvalues; stop iteration by setting NELMGW= 0
!       IF (ABS(DESUM1) < ABS(INFO%EDIFF) ) NELMGW = 0

       IF (IO%IU0>=0) WRITE(17, 200)  NELM,0._q,0._q,DESUM1,ICOUEV,RMS
       IF (IO%IU0>=0) &
         WRITE(IO%IU0, 200)  NELM,0._q,0._q,DESUM1,ICOUEV,RMS

  200 FORMAT('DIA: ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5, &
          &       I6,'  ',E10.3)

       ! to be save that eigenvalues are fully converged iterated 1 more times
       ! unfortunately EDDAV breaks sometimes if many bands are calculated...
       DO K1=1,1
!          CALL EDDAV(HAMILTONIAN,P, GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WDES, NSIM, &
!            LMDIM,CDIJ,CQIJ, RMS,DESUM1,ICOUEV, SV,E%EXHF, IO%IU6,IO%IU0, .FALSE., .TRUE., .FALSE.)
       ENDDO

       ! update occupancies (well this is done below as well, but we need it for eigenvalues)
       !keep fermi energy as it is for finite temperature calculations 
       CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
         INFO%NUP_DOWN,  E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
         NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)

       IF (IO%IU6>=0) WRITE(IO%IU6,2202) EFERMI
       CALL WRITE_EIGENVAL_NBANDS( WDES, W, IO%IU6, NBANDSGW)

2202   FORMAT(' E-fermi : ', F8.4)

       IF (SYMM%ISYM>0) CALL CLEANUP_CELEN(W) ! set degenerated or near degenerated eigenvalues to a single value

       IF (.NOT. LGWLF) THEN
          CALL POP_XC_TYPE
          IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
             CALL SETUP_LDA_XC(2,-1,-1,IO%IDIOT)
          ELSE
             CALL SETUP_LDA_XC(1,-1,-1,IO%IDIOT)
          ENDIF
       ENDIF

       ! deallocate all
       CALL DEALLOCW( WACC1)
       CALL DEALLOCW( WACC2)


       IF (.NOT. LGW0 .AND. IO%LOPTICS) THEN
          CALL START_TIMING("G")
         ! VASP onboard optics
          CALL PEAD_RESETUP_WDES(WDES,GRID,KPOINTS,LATT_CUR,LATT_INI,IO)

          CALL LR_OPTIC( &
             P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
             T_INFO,INFO,IO,KPOINTS,SYMM,GRID,GRID_SOFT, &
             GRIDC,GRIDUS,C_TO_US,SOFT_TO_C, &
             CHTOT,DENCOR,CVTOT,CSTRF, &
             CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM, &
             CHDEN,SV,LMDIM,IRDMAX,EFERMI,NEDOS, & 
             LSTORE=.TRUE., LPOT=.FALSE.)
          CALL STOP_TIMING("G",IO%IU6,'OPTICS')
        ENDIF

       ! restore full k-point grid
       IF (SYMM%ISYM>=0 .AND. .NOT. LGAMMA) THEN
          CALL NOSYMM(LATT_CUR%A,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS,SYMM%PTRANS, &
               &   SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%MAGROT,INFO%ISPIN,IO%IU6)

          CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR,.FALSE., &
               T_INFO%NIONS,SYMM%ROTMAP, SYMM%MAGROT, SYMM%ISYM,-1,IO%IU0)
          CALL KPAR_SYNC_ALL(WDES,W)
          CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
          CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
       ENDIF
       IF (.NOT. LGW0) THEN
          W_W=W
       ENDIF
       CALL STOP_TIMING("G",IO%IU6,"EDDAV")

       CALL DEALLOCATE_RESPONSEFUN(CHI0)
    ENDIF

    IF (LFXCEPS) THEN
       CALL DEALLOCATE_RESPONSEFUN(CHI0)
    ENDIF
    
    ! update occupancies
    CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
         INFO%NUP_DOWN,  E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
         NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)

    CALL STOP_TIMING("GWLOOP",IO%IU6,XMLTAG='total')

!-----------------------------------------------------------------------
! ok, if the occupancies have changed
! we need to recalculate the Hartree and exchange terms
!-----------------------------------------------------------------------
    IF (NELM<NELMGW) THEN
       IF (ABS(MAXVAL(FERTOT_INITIAL-W%FERTOT))>1E-2 .AND. .NOT. LscQPGW) THEN
          FERTOT_INITIAL=W%FERTOT
          ! restore symmetry reduced k-point grid
          IF (SYMM%ISYM>=0 .AND. .NOT. LGAMMA) THEN
          IF (SYMM%ISYM>0) &
             CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
                  T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
                  SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
                  SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,-1)
             CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
                  SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
                  T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,-1,IO%IU0)
             CALL KPAR_SYNC_ALL(WDES,W)
             CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
             CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
          ENDIF
          ! use HF exchange only
          IF (.NOT. LGWLF) THEN
             CALL PUSH_XC_TYPE_FOR_GW
             IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
                CALL SETUP_LDA_XC(2,-1,-1,IO%IDIOT)
             ELSE
                CALL SETUP_LDA_XC(1,-1,-1,IO%IDIOT)
             ENDIF
             CALL SET_FSG_STORE( GRIDHF, LATT_CUR, WDES)
!             DO K1=1,WDES%NKPTS
!                FSG_STORE(K1)=SET_FSG(GRIDHF, LATT_CUR, K1)
!             ENDDO
          ENDIF
          ! update charge
          CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
               GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
               LATT_CUR, P, SYMM, T_INFO, &
               CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

          CALL SET_EIGENVALUE_HARTREE_KINETIC( &
               HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR, &
               T_INFO,INFO,IO,KPOINTS,GRID,GRID_SOFT, &
               GRIDC,GRIDUS,C_TO_US,SOFT_TO_C,SYMM, &
               CHTOT,DENCOR,CVTOT,CSTRF, &
               CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM, &
               CHDEN,SV,LMDIM,IRDMAX)
          IF (.NOT. LGWLF) THEN
             CALL POP_XC_TYPE
             IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
                CALL SETUP_LDA_XC(2,-1,-1,IO%IDIOT)
             ELSE
                CALL SETUP_LDA_XC(1,-1,-1,IO%IDIOT)
             ENDIF
          ENDIF

          ! switch off symmetry
          IF (SYMM%ISYM>=0 .AND. .NOT. LGAMMA) THEN
             CALL NOSYMM(LATT_CUR%A,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS,SYMM%PTRANS, &
                  &   SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%MAGROT,INFO%ISPIN,IO%IU6)

             CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR,.FALSE., &
                  T_INFO%NIONS,SYMM%ROTMAP, SYMM%MAGROT, SYMM%ISYM,-1,IO%IU0)
             CALL KPAR_SYNC_ALL(WDES,W)
             CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
             CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
          ENDIF
          IF (.NOT. LGW0) THEN
             W_W=W
          ELSE
             ! W_W needs to be relinked
             W_W%CPTWFP=>W%CPTWFP
             W_W%CPROJ =>W%CPROJ
             W_W%AUX   =>W%AUX
             W_W%AUXTOT=>W%AUXTOT
          ENDIF
       ENDIF
       
       ! ===============================================================
       ! next electronic step GOTO 2000, really ugly
       ! ===============================================================
       GOTO 2000
    ENDIF

    CALL DEALLOCATE_LOCAL_FIELD_FOCK

    IF (.NOT. LCHIREALLOCATE) THEN
       CALL DEALLOCATE_RESPONSEFUN( CHI )
    ENDIF
!-----------------------------------------------------------------------
! go back to original symmetry
!-----------------------------------------------------------------------
    DEALLOCATE(OMEGA, OMEGAWEIGHT, COMEGA, WGW)
    IF (LGW) THEN
       CALL DEREGISTER_ALLOCATE(8._q*SIZE(SCREENED_TWO_ELECTRON_INTEGRALS), "2eintegral")
       DEALLOCATE(SCREENED_TWO_ELECTRON_INTEGRALS)
    ENDIF

1000 CONTINUE

    IF (SYMM%ISYM>=0.AND. .NOT. LGAMMA) THEN
       IF (SYMM%ISYM>0) &
       CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
            T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
            SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
            SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)
       CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
            SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
            T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,-1,IO%IU0)
       CALL KPAR_SYNC_ALL(WDES,W)
       CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
       CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
    ENDIF
    IF (LscQPGW .AND. .NOT. LGW0 .AND. IO%LWAVE) THEN
       CALL WRT_CDER_BETWEEN_STATES(WDES, IO%IU0, 55) ! write updated nabla
    ENDIF

    !also get rid of auxilary wave functions in CRPA
    IF ( LCRPA ) THEN
       IF ( ALLOCATED( UCRPA ) ) DEALLOCATE(UCRPA)
       IF ( LALLOCATED_WCORR ) THEN
          CALL DEALLOCW( W_CORR )
          LALLOCATED_WCORR = .FALSE.
       ENDIF 
    ENDIF

    IF (IO%LWAVE .AND. .NOT. LCRPA  ) THEN
       CALL OUTWAV(IO, WDES, W, LATT_CUR, EFERMI )
    ENDIF

    IF (IO%LWAVEH5 .AND. .NOT. LCRPA ) THEN
#ifdef VASP_HDF5
       IH5ERR = VH5_FILE_DELETE('vaspwave.h5', 112)
       IF (IH5ERR /= 0) THEN
         WRITE(*,*) "ERROR deleting old wave file, ierr = ", IH5ERR
       ENDIF
       IH5ERR = VH5_FILE_CREATE_OR_OVERWRITE('vaspwave.h5', IH5WAVEFILEID)
       IF (IH5ERR == 0) THEN
         CALL VH5_WRITE_WAVEFUNCTIONS(IH5WAVEFILEID, IO, WDES, W, LATT_INI, EFERMI)
         IH5ERR = VH5_FILE_CLOSE_WRITING(IH5WAVEFILEID)
       ELSE
         WRITE(IO%IU0,*) "Error: could not create vaspwave.h5, errorcode: ",IH5ERR
       ENDIF
#endif
    ENDIF

!=======================================================================
! write the local potential
!=======================================================================
    IF (IO%LVTOT .AND. LOEP) THEN
#ifdef MPI   
      IF (WDES%COMM_KINTER%NODE_ME.EQ.1) THEN 
#endif
       IF (IO%IU0>=0) THEN
          IF (IO%LOPEN) OPEN(IO%IUVTOT,FILE='POT',STATUS='UNKNOWN')
          REWIND IO%IUVTOT
          CALL OUTPOS(IO%IUVTOT,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A,.FALSE.,DYN%POSION)
       ENDIF
       
       ! at the moment the spin up and down potential is written to the file
       CALL SET_RHO_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
            CDIJ, DLM_EXX)
       CALL OUTPOT(GRIDC, IO%IUVTOT,.TRUE.,CVTOT)
       CALL WRT_RHO_PAW(P, T_INFO, INFO%LOVERL, DLM_EXX(:,1), GRIDC%COMM, IO%IUVTOT)
       
       DO ISP=2,WDES%NCDIJ
          IF (IO%IU0>=0) THEN
             WRITE( IO%IUVTOT,'(5E20.12)') (T_INFO%ATOMOM(I),I=1,T_INFO%NIONS)
          ENDIF
          CALL OUTPOT(GRIDC, IO%IUVTOT,.TRUE.,CVTOT(1,ISP))
          
          CALL WRT_RHO_PAW(P, T_INFO, INFO%LOVERL, DLM_EXX(:,ISP), GRIDC%COMM, IO%IUVTOT )
       ENDDO
#ifdef MPI   
      ENDIF
#endif
    ENDIF
!-----------------------------------------------------------------------
! some deallocation
!-----------------------------------------------------------------------
    CALL DEALLOCATE_IRZ_MAP(S2E)

    ! set XC type back to HF as it was upon entry to the routine
    IF (.NOT. LGWLF) THEN
       CALL PUSH_XC_TYPE_FOR_GW
       IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
          CALL SETUP_LDA_XC(2,-1,-1,IO%IDIOT)
       ELSE
          CALL SETUP_LDA_XC(1,-1,-1,IO%IDIOT)
       ENDIF
    ENDIF

    PROFILING_STOP('calculate_xi')

    CONTAINS

!**********************************************************************
!
! calculate correlation energy from Bethe Salpeter equation
!
! the actual course of the calculations depends on LAMBDA
! for LAMBDA> 1 and LAMBDA an integer: a  GL coupling constant integration 
!     is performed by CALCULATE_BSE using LAMBDA data points
! for LAMBDA/=1.0: the correlation energy for a single point
!     in the coupling constant integration is evaluated by CALCULATE_BSE
! 
!**********************************************************************

  SUBROUTINE CALCULATE_COR_BSE()
    INTEGER :: NQ_COUNTER, NQ
    REAL(q) :: CORREP, CORREP_kp

       ! sets DELTA_COND in local_field.F
       CALL CALCULATE_LOCAL_FIELD_PREPARE( 1, W, WGW, LBSE, LGWLF, & 
            .FALSE.,  & 
            IDIR_MAX, COMEGA, TBSE, TBSEA, ANTIRES, & 
            LATT_CUR, NKREDLFX, NKREDLFY, NKREDLFZ, IO%IU0, IO%IU6 )


       CORREP=0.0_q
! loop over all difference vectors between k-points
       DO NQ_COUNTER=1,S2E%NUMBER_OF_NQ
          CORREP_kp=0.0_q
      ! determine corresponding point in the VASP lists
          NQ=S2E%NQ(NQ_COUNTER)  ! storage position 
     
          CALL CALCULATE_LOCAL_FIELD_DFT( &
              GRIDC, WDES, LATT_CUR, CHTOT, DENCOR, &
              WGW, NQ, 1 )
       
       ! if LGWLF is set, ALDAX and ALDAC are zero at this point
       ! hence the DFT loCal field effects f_xc (stored in FXCR) are exactly zero
         IF (WDES%ISPIN >1) THEN
          ! take into account all spins if the variable after NQ is -1 otherwise:
          ! negative number: both spin channels are treated simultaneously (only valid for ISPIN=2)
          ! ISP_IN=1 : treat only spin channel 1 (and pretend the system is non-magnetic)
          ! ISP_IN=2 : treat only spin channel 2 (and pretend the system is non-magnetic)
            CALL CALCULATE_BSE( &
               IBSE,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
               T_INFO,DYN,INFO,IO,KPOINTS,SYMM,GRID,LMDIM,CQIJ, &
               WGW, SHIFT, CELTOT_HARTREE_KINETIC, CELTOT_X, & 
               NBANDSO, NBANDSV, OMEGAMAX, NQ,  KPOINT_BSE(2:4), -1, &
               LHARTREE, LADDER, LTRIPLET, ANTIRES, LGWLF, LFXC, NEDOS, NELMGW, NBANDSGW, L2ORDER, &
               LCORBSE, LQPBSE, LTCTC, CORREP_kp, LAMBDA)
         ELSE
    ! take into account only first (and only) spin component
            CALL CALCULATE_BSE( &
               IBSE,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
               T_INFO,DYN,INFO,IO,KPOINTS,SYMM,GRID,LMDIM,CQIJ, &
               WGW, SHIFT, CELTOT_HARTREE_KINETIC, CELTOT_X,  & 
               NBANDSO, NBANDSV, OMEGAMAX, NQ,  KPOINT_BSE(2:4), 1, &
               LHARTREE, LADDER, LTRIPLET, ANTIRES, LGWLF, LFXC, NEDOS, NELMGW, NBANDSGW, L2ORDER, &
               LCORBSE, LQPBSE, LTCTC, CORREP_kp, LAMBDA)
         ENDIF
!added loop over k-points
         CORREP=CORREP+CORREP_kp
      ENDDO
!end added loop
       IF (LAMBDA > 1.0_q .AND. LAMBDA==INT(LAMBDA)) THEN
          IF (IO%IU6>=0) WRITE(IO%IU0,'(1X,A,F20.10)') 'converged correlation energy from BSE', CORREP 
          IF (IO%IU6>=0) WRITE(IO%IU6,'(1X,A,F20.10)') 'converged correlation energy from BSE', CORREP 
       ELSE 
          IF (IO%IU0>=0) WRITE(IO%IU0,'(1X,A,F20.10,A,F10.4)') 'Integrand AC path', CORREP, ' @ lambda', LAMBDA  
          IF (IO%IU6>=0) WRITE(IO%IU6,'(1X,A,F20.10,A,F10.4)') 'Integrand AC path', CORREP, ' @ lambda', LAMBDA
       ENDIF
     END SUBROUTINE CALCULATE_COR_BSE

!*********************************************************************
!
! internal subroutine  CALCULATE_QPBSE
! this subroutine calculates the QP energies by solving the BSE 
! equations and calculating G from the polarization propagator
!
!*********************************************************************

     SUBROUTINE CALCULATE_QPBSE()

       TYPE (selfenergy_from_bse) :: SEFROMBSE         ! self-energy from BSE at NQ
       REAL(q), ALLOCATABLE, TARGET :: RE_SE(:,:,:,:)  ! real part of self-energy summed over all k-points
       REAL(q), ALLOCATABLE, TARGET :: IM_SE(:,:,:,:)  ! imaginary part of self-energy summed over all k-points
       REAL(q), ALLOCATABLE :: A(:,:)                  ! local structure that allows to print one band at the time  
       INTEGER :: NELM, NELMGW_COUNT                   ! counter for self-consisten QP shift calculation
       REAL(q) :: CORREP_kp                            ! actually not needed
       INTEGER :: I

       !initialise the data structure for selfenergy_from_bse type
       SEFROMBSE%NOMEGA=NOMEGABSE
       SEFROMBSE%NBANDS=NBANDSGW
       ! this is tricky we only need the spectral function inside the irreducible wedge of the BZ
       SEFROMBSE%NKPTS =KPOINTS_ORIG%NKPTS
       SEFROMBSE%ISPIN =WDES%ISPIN

       ALLOCATE(A(3,NOMEGABSE) )
       ALLOCATE(RE_SE(NOMEGABSE, SEFROMBSE%NBANDS, SEFROMBSE%NKPTS, SEFROMBSE%ISPIN))
       ALLOCATE(IM_SE(NOMEGABSE, SEFROMBSE%NBANDS, SEFROMBSE%NKPTS, SEFROMBSE%ISPIN))

       IF (IO%IU0>=0) WRITE(IO%IU0,*) 'Self-energy computed', LQPBSE 
       IF (IO%IU0>=0) WRITE(IO%IU0,*) 'Exchange present', LADDER 
       IF (IO%IU0>=0) WRITE(IO%IU0,*) 'GW TC-TC calculation', LTCTC 
       IF (IO%IU0>=0) WRITE(IO%IU0,*) '# Omega points', NOMEGABSE

       IF (LTCTC .AND. .NOT. LADDER) THEN
          IF (IO%IU0>=0) WRITE(IO%IU0,*) 'GW TC-TC calculation only allowed with LADDER = .TRUE.' 
       ENDIF 

       ! sets DELTA_COND in local_field.F
       CALL CALCULATE_LOCAL_FIELD_PREPARE( 1, W, WGW, LBSE, LGWLF, & 
            .FALSE.,  & 
            IDIR_MAX, COMEGA, TBSE, TBSEA, ANTIRES, & 
            LATT_CUR, NKREDLFX, NKREDLFY, NKREDLFZ, IO%IU0, IO%IU6 )

       ! allocate and initialise the QP energies (for self-consistent evaluation) 
       ! if you want the BSE matrix to be included in the self-consistency cycle 
       ! point to W%CELTOT 
       SEFROMBSE%CELNEW => W%CELTOT
 
       ! in the gamma only version the loop over NELMGW is performed in the BSE
       ! module 
       IF (W%WDES%LGAMMA) THEN
          NELMGW_COUNT=1
       ELSE 
          NELMGW_COUNT=NELMGW
       ENDIF 
       ! introduce loop over NELM 
       DO NELM = 1, NELMGW_COUNT
         IF (.NOT. ASSOCIATED(SEFROMBSE%REAL_PART)) THEN
            ALLOCATE(SEFROMBSE%REAL_PART(SEFROMBSE%NOMEGA,SEFROMBSE%NBANDS,SEFROMBSE%NKPTS,SEFROMBSE%ISPIN))
            ALLOCATE(SEFROMBSE%IMAG_PART(SEFROMBSE%NOMEGA,SEFROMBSE%NBANDS,SEFROMBSE%NKPTS,SEFROMBSE%ISPIN))
            ALLOCATE(SEFROMBSE%OMEGA(SEFROMBSE%NOMEGA))
         ENDIF

         !simple linear grid from -OMEGAMAX to OMEGAMAX
         DO I=1,NOMEGABSE
            SEFROMBSE%OMEGA(I)=-OMEGAMAX+((2*OMEGAMAX)/(NOMEGABSE-1))*(I-1)
         ENDDO

         !initialisation
         A=0.0_q ; RE_SE=0.0_q ; IM_SE=0.0_q ! accumulation vectors over q-pts for selfenergy
  ! loop over all difference vectors between k-points
         DO NQ_COUNTER=1,S2E%NUMBER_OF_NQ
#ifdef MPI
            IF (MOD(NQ_COUNTER,WDES%COMM_KINTER%NCPU) /= WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
            CORREP_kp=0.0_q
            NQ=S2E%NQ(NQ_COUNTER) ! storage position 

            ! initialize the exchange correlation kernel
            ! if LGWLF is set, ALDAX and ALDAC are zero at this point
            ! hence the DFT local field effects f_xc (stored in FXCR) are exactly zero
            CALL CALCULATE_LOCAL_FIELD_DFT( &
                 GRIDC, WDES, LATT_CUR, CHTOT, DENCOR, &
                 WGW, NQ, 1 )

            IF (WDES%ISPIN >1) THEN
               ! treat up and down spin simultaneously
               CALL CALCULATE_BSE( &
                    IBSE,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
                    T_INFO,DYN,INFO,IO,KPOINTS,SYMM,GRID,LMDIM,CQIJ, &
                    WGW, SHIFT, CELTOT_HARTREE_KINETIC, CELTOT_X, & 
                    NBANDSO, NBANDSV, OMEGAMAX, NQ,  KPOINT_BSE(2:4), -1, &
                    LHARTREE, LADDER, LTRIPLET, ANTIRES, LGWLF, LFXC, NEDOS, NELMGW, NBANDSGW, L2ORDER, &
                    LCORBSE, LQPBSE, LTCTC, CORREP_kp, LAMBDA, S2E, EFERMI, SEFROMBSE)
            ELSE
               ! take into account only first (and only) spin component
               CALL CALCULATE_BSE( &
                    IBSE,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
                    T_INFO,DYN,INFO,IO,KPOINTS,SYMM,GRID,LMDIM,CQIJ, &
                    WGW, SHIFT, CELTOT_HARTREE_KINETIC, CELTOT_X, & 
                    NBANDSO, NBANDSV, OMEGAMAX, NQ,  KPOINT_BSE(2:4), 1, &
                    LHARTREE, LADDER, LTRIPLET, ANTIRES, LGWLF, LFXC, NEDOS, NELMGW, NBANDSGW, L2ORDER, &
                    LCORBSE, LQPBSE, LTCTC, CORREP_kp, LAMBDA, S2E, EFERMI, SEFROMBSE)
            ENDIF
            ! just add up self-energies from different q-points
            DO ISP=1, SEFROMBSE%ISPIN
               DO NK=1, SEFROMBSE%NKPTS
                  DO I=1, SEFROMBSE%NBANDS
                     RE_SE(:,I,NK,ISP)=RE_SE(:,I,NK,ISP)+ SEFROMBSE%REAL_PART(:,I,NK,ISP)*WDES%NKPTS
                     IM_SE(:,I,NK,ISP)=IM_SE(:,I,NK,ISP)+ SEFROMBSE%IMAG_PART(:,I,NK,ISP)*WDES%NKPTS 
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         DEALLOCATE(SEFROMBSE%REAL_PART)
         DEALLOCATE(SEFROMBSE%IMAG_PART)
  
         ! now point it to RE_SE since the Hilbert transformation routine will use SEFROMBSE
         SEFROMBSE%REAL_PART => RE_SE
         SEFROMBSE%IMAG_PART => IM_SE

         ! now perform a Hilbert transform to obtain self-energy from it's spectral representation
         CALLMPI( M_sum_d(WDES%COMM_KINTER, SEFROMBSE%REAL_PART(1,1,1,1), SIZE(SEFROMBSE%REAL_PART)))
         CALLMPI( M_sum_d(WDES%COMM_KINTER, SEFROMBSE%IMAG_PART(1,1,1,1), SIZE(SEFROMBSE%IMAG_PART)))

         ! clean up degenerate orbitals
         CALL CLEANUP_CELEN_SELFE(W, NBANDSGW, SEFROMBSE%REAL_PART)
         CALL CLEANUP_CELEN_SELFE(W, NBANDSGW, SEFROMBSE%IMAG_PART)

         CALL SE_HILBERT_TRANSF(SEFROMBSE, EFERMI, SHIFT, CELTOT_X, CELTOT_HARTREE_KINETIC)

         ! compute the quasi-particle energy for state |n> as energy where Re(Sigma(w)) and y=w intersect
         IF (IO%IU6 >=0) WRITE (IO%IU6,'(//A, I2)') ' QP shifts <psi_nk| G(iteration)W_0 |psi_nk>: iteration', NELM
         ! the QP energies are updated and overwritten in W%CELTOT
         CALL QPE_FROMSIGMA(W, SEFROMBSE, EFERMI, CELTOT_X, CELTOT_HARTREE_KINETIC, NELM, IO%IU0, IO%IU6) 
         IF (IO%IU6 >=0) WRITE (IO%IU6,130) 

         ! update the Fermi level as well as the partial occupancies
         !keep fermi energy as it is for finite temperature calculations 
         IF ( .NOT. LFINITE_TEMPERATURE ) &
         CALL DENSTA( IO%IU0, IO%IU6, WDES, W, KPOINTS, INFO%NELECT, &
              INFO%NUP_DOWN,  E%EENTROPY, EFERMI, KPOINTS%SIGMA, .FALSE., &
              NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
    !-----------------------------------------------------------------------
    ! if the occupancies have changed
    ! we need to recalculate the Hartree and exchange terms
    !-----------------------------------------------------------------------
         IF (NELM<NELMGW) THEN
            ! IF (ABS(MAXVAL(FERTOT_INITIAL-W%FERTOT))>1E-2 .AND. .NOT. LscQPGW) THEN
            IF (ABS(MAXVAL(FERTOT_INITIAL-W%FERTOT))>1E-2) THEN
               FERTOT_INITIAL=W%FERTOT
               ! restore symmetry reduced k-point grid
               IF (SYMM%ISYM>=0 .AND. .NOT. LGAMMA) THEN
               IF (SYMM%ISYM>0) &
                  CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
                       T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
                       SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
                       SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,-1)
                  CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
                       SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
                       T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,-1,IO%IU0)
                  CALL KPAR_SYNC_ALL(WDES,W)
                  CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
                  CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
               ENDIF

               ! use HF exchange only
               IF (.NOT. LGWLF) THEN
                  CALL PUSH_XC_TYPE_FOR_GW
                  IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
                     CALL SETUP_LDA_XC(2,-1,-1,IO%IDIOT)
                  ELSE
                     CALL SETUP_LDA_XC(1,-1,-1,IO%IDIOT)
                  ENDIF
                  CALL SET_FSG_STORE( GRIDHF, LATT_CUR, WDES)
!                  DO K1=1,WDES%NKPTS
!                     FSG_STORE(K1)=SET_FSG(GRIDHF, LATT_CUR, K1)
!                  ENDDO
               ENDIF
               ! update charge
               CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
                    GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
                    LATT_CUR, P, SYMM, T_INFO, &
                    CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)
     
               CALL SET_EIGENVALUE_HARTREE_KINETIC( &
                    HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR, &
                    T_INFO,INFO,IO,KPOINTS,GRID,GRID_SOFT, &
                    GRIDC,GRIDUS,C_TO_US,SOFT_TO_C,SYMM, &
                    CHTOT,DENCOR,CVTOT,CSTRF, &
                    CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM, &
                    CHDEN,SV,LMDIM,IRDMAX)
               IF (.NOT. LGWLF) THEN
                  CALL POP_XC_TYPE
                  IF (WDES%LNONCOLLINEAR .OR. INFO%ISPIN == 2) THEN
                     CALL SETUP_LDA_XC(2,-1,-1,IO%IDIOT)
                  ELSE
                     CALL SETUP_LDA_XC(1,-1,-1,IO%IDIOT)
                  ENDIF
               ENDIF

               ! switch off symmetry
               IF (SYMM%ISYM>=0 .AND. .NOT. LGAMMA) THEN
                  CALL NOSYMM(LATT_CUR%A,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS,SYMM%PTRANS, &
                       SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%MAGROT,INFO%ISPIN,IO%IU6)
                  CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR,.FALSE., &
                       T_INFO%NIONS,SYMM%ROTMAP, SYMM%MAGROT, SYMM%ISYM,-1,IO%IU0)
                  CALL KPAR_SYNC_ALL(WDES,W)
                  CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
                  CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
               ENDIF
               ! of we perform LGW0 calculations W_W 
               IF (.NOT. LGW0) THEN
                  W_W=W
               ELSE
                  ! W_W needs to be relinked
                  W_W%CPTWFP=>W%CPTWFP
                  W_W%CPROJ =>W%CPROJ
                  W_W%AUX   =>W%AUX
                  W_W%AUXTOT=>W%AUXTOT
               ENDIF
            ENDIF
         ENDIF

    ! write selfenergy from BSE calculation
         IF (IO%IU6>=0 .AND. NELM==NELMGW_COUNT) THEN
           DO ISP=1, SEFROMBSE%ISPIN
             DO NK=1, S2E%NUMBER_OF_NQ ! SEFROMBSE%NKPTS
               A = 0.0_q
               DO I=1, SEFROMBSE%NBANDS
                 A(1,:)=  SEFROMBSE%OMEGA(:)
                 A(2,:)=  SEFROMBSE%REAL_PART(:,I,NK,ISP)
                 A(3,:)=  SEFROMBSE%IMAG_PART(:,I,NK,ISP)
                 CALL XML_VECARRAY("selfenergy")
                 CALL XML_ARRAY_REAL(A)
                 CALL XML_CLOSE_TAG
               ENDDO
             ENDDO
           ENDDO
         ENDIF

    ! write spectral density from BSE calculation
    ! gK: commented this out, can be done by simple postprocessing from self-energy
    !     IF (IO%IU6>=0 .AND. NELM==NELMGW_COUNT) THEN
    !       DO ISP=1, SEFROMBSE%ISPIN
    !         DO NK=1, S2E%NUMBER_OF_NQ ! SEFROMBSE%NKPTS
    !           A = 0.0_q
    !           A(1,:)= SEFROMBSE%OMEGA(:)
    !           DO I=1, SEFROMBSE%NBANDS
    !             A(3,:)=  ((SEFROMBSE%REAL_PART(:,I,NK,ISP))**2)+((SEFROMBSE%IMAG_PART(:,I,NK,ISP))**2)
    !             A(2,:)=  A(2,:)+(ABS(SEFROMBSE%IMAG_PART(:,I,NK,ISP))/A(3,:))
    !           ENDDO
    !           CALL XML_VECARRAY("spectraldens")
    !           CALL XML_ARRAY_REAL(A)
    !           CALL XML_CLOSE_TAG
    !         ENDDO
    !       ENDDO
    !     ENDIF

         IF (ASSOCIATED(SEFROMBSE%REAL_PART)) THEN
            NULLIFY(SEFROMBSE%REAL_PART)
            NULLIFY(SEFROMBSE%IMAG_PART)
            NULLIFY(SEFROMBSE%OMEGA)
         ENDIF

       ENDDO !NELM
       DEALLOCATE( A, RE_SE, IM_SE)
       !DEALLOCATE(SEFROMBSE%CELNEW)

130 FORMAT (5X, //, &
         &'----------------------------------------------------', &
         &'----------------------------------------------------'//)

     END SUBROUTINE CALCULATE_QPBSE


  END SUBROUTINE CALCULATE_XI


!*********************************************************************
!
! set NBANDSGW
!
!*********************************************************************

  SUBROUTINE SET_NBANDSGW(W)
    TYPE (wavespin)    W
! determine the number of bands required for the HF type calculations
    IF (NBANDSGW<=0) THEN
       ! default calculate selfenergy shift for twice the number of occupied bands
       NBANDSGW=((MIN(LAST_FILLED_OPTICS(W)*2,W%WDES%NB_TOT)+W%WDES%NB_PAR-1)/W%WDES%NB_PAR)*W%WDES%NB_PAR
       IF (LGWNO) THEN
          NBANDSGW= W%WDES%NB_TOT
       ENDIF
    ELSE
       NBANDSGW=((NBANDSGW+W%WDES%NB_PAR-1)/W%WDES%NB_PAR)*W%WDES%NB_PAR
    ENDIF

    NBANDSGW=MIN(NBANDSGW,W%WDES%NB_TOT)
  END SUBROUTINE SET_NBANDSGW

!*********************************************************************
!
! main wavefunction based correlation treatment ((MP2,MP2NO,..) scheduler
! does not do a lot but to lower symmetry and call the
! calculations procedure
!
!*********************************************************************

  SUBROUTINE PSI_CORR_MAIN( &
          P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI,  &
          T_INFO, DYN, IO,KPOINTS,SYMM,GRID,INFO,AMIX,BMIX)
    USE mp2
    USE mp2kpar
    USE ump2no
    USE fcidump
!aG: to be added
#ifdef vasp6
!    USE ccsd
!    USE uccsd
!    USE bracketst
    USE rpax
!    USE grpax
!
#endif
    USE base
    USE ini
    USE lattice
    USE pseudo
    USE lattice
    USE nonl_high
    USE msymmetry
    USE poscar
    USE wave_high
    USE full_kpoints
    USE pead, ONLY : LPEAD_SYM_RED

    IMPLICIT NONE

    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    TYPE (wavedes)     WDES
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W
    TYPE (dynamics)    DYN
    TYPE (latt)        LATT_CUR, LATT_INI
    TYPE (in_struct)   IO
    TYPE (kpoints_struct) KPOINTS
    TYPE (symmetry)    SYMM
    TYPE (grid_3d)     GRID       ! grid for wavefunctions
    TYPE (info_struct) INFO
    REAL(q) :: AMIX, BMIX
   ! local
    TYPE (wavedes), POINTER :: WGW
    TYPE (grid_3d), POINTER :: GRIDWGW

    ! aG: the convergence acceleration is not always advantagous, I would prefer to
    ! neglect it.
    CALL READ_CDER_BETWEEN_STATES(WDES, IO%IU0, 55)


#ifdef MPI
    IF (WDES%COMM_KINTER%NCPU.NE.1) THEN
       CALL vtutor%error("PSI_CORR_MAIN: KPAR>1 not implemented, sorry.")
    END IF
#endif

    IF (SYMM%ISYM>=0.AND. .NOT. WDES%LGAMMA) THEN
       ! switch of symmetry
       CALL NOSYMM(LATT_CUR%A,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS,SYMM%PTRANS, &
            &   SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%MAGROT,WDES%ISPIN,IO%IU6)
       
       ! reread k-points with LINVERSION=.FALSE.  to generate full mesh
       CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR,.FALSE., &
            T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IO%IU6,IO%IU0)
       CALL KPAR_SYNC_ALL(WDES,W)
       CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_CUR,-1, IO%IU0)
       CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
    ENDIF
!=======================================================================
! generate descriptor for response function
!=======================================================================
    ALLOCATE(WGW, GRIDWGW)
    WGW=WDES_FOCK

    WGW%NKPTS=KPOINTS_FULL%NKPTS
    WGW%NKDIM=KPOINTS_FULL%NKPTS
    WGW%NKPTS_FOR_GEN_LAYOUT=KPOINTS_FULL%NKPTS
    ALLOCATE(WGW%VKPT(1:3,SIZE(KPOINTS_FULL%VKPT,2)),WGW%WTKPT(SIZE(KPOINTS_FULL%WTKPT,1)))
    WGW%VKPT =KPOINTS_FULL%VKPT
    WGW%WTKPT=KPOINTS_FULL%WTKPT

!    IF (ENCUTGW<WDES%ENMAX .AND. ENCUTGW/=-1) THEN
      WGW%ENMAX=ENCUTGW
!    ENDIF

    ! GRIDWGW is identical to GRID_FOCK, except for GRIDWGW%FFTSCA
    GRIDWGW=GRID_FOCK
    IF (IO%IU6>=0) THEN 
       WRITE(IO%IU6,*) 'Basis sets for responsefunctions:'
       WRITE(IO%IU6,*) '================================='
    ENDIF
    CALL GEN_LAYOUT(GRIDWGW, WGW, LATT_CUR%B, LATT_CUR%B, IO%IU6,.TRUE.)

    IF (WGW%LGAMMA) THEN
       ! gamma only data layout with wavefunction stored as real in real space
       GRIDWGW%LREAL=.TRUE.
    ENDIF

    CALL GEN_INDEX (GRIDWGW, WGW, LATT_CUR%B, LATT_CUR%B,IO%IU6,-1, .TRUE.)
    !  init FFT (required if real to complex FFT is used)
    CALL FFTINI(WGW%NINDPW(1,1), WGW%NGVECTOR(1), WGW%NKPTS, WGW%NGDIM, GRIDWGW)


    IF (LMP2) THEN
#ifdef scaLAPACK
       CALL CALCULATE_MP2( &
          P,WDES,W,LATT_CUR,T_INFO,IO,KPOINTS, WGW, ENCUTGW, ENCUTGWSOFT, LMAXMP2, SYMM)
#endif
    ELSEIF (LMP2KPAR) THEN
       CALL CALCULATE_MP2_KPAR( &
          P,WDES,W,LATT_CUR,T_INFO,IO,KPOINTS,WGW, ENCUTGW, ENCUTGWSOFT, LMAXMP2)
    ELSEIF (LMP2NO) THEN
#ifdef scaLAPACK
       CALL CALCULATE_FNO(P,WDES,W,LATT_INI,LATT_CUR,T_INFO,IO,KPOINTS,WGW, ENCUTGW, &
         & ENCUTGWSOFT, LMAXMP2,INFO)
#endif
    ELSEIF (LFCIDUMP) THEN
#ifdef scaLAPACK
       CALL CALC_FCIDUMP(P,WDES,W,LATT_CUR,T_INFO,IO,KPOINTS,WGW, ENCUTGW, &
         & ENCUTGWSOFT, LMAXMP2, INFO)
#endif
#ifdef vasp6
!    ELSEIF (LCCSD) THEN
! #ifdef gammareal
!       CALL CALC_UCCSD(P,WDES,W,LATT_CUR,T_INFO,IO,KPOINTS,WGW, ENCUTGW, &
!        & ENCUTGWSOFT, LMAXMP2, INFO)
! #else
!       CALL CALCULATE_CCSD(P,WDES,W,LATT_CUR,T_INFO,IO,KPOINTS,WGW,ENCUTGW, ENCUTGWSOFT, &
!         & LMAXMP2,AMIX,BMIX)
! #endif
!    ELSEIF (LBRACKETST) THEN
!       CALL CALCULATE_BRACKETST(P,WDES,W,LATT_CUR,T_INFO,IO,KPOINTS,WGW, &
!         & LMAXMP2)
    ELSEIF (LRPAX) THEN
      CALL CALCULATE_RPAX(P,WDES,W,LATT_CUR,T_INFO,IO,KPOINTS,WGW,ENCUTGW, &
        & ENCUTGWSOFT, LMAXMP2,NBANDSGW)
#endif
    ENDIF

    IF (SYMM%ISYM>=0.AND. .NOT. WDES%LGAMMA) THEN
       IF (SYMM%ISYM>0) &
       CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
            T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
            SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
            SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,-1)
       CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
            SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
            T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,-1,IO%IU0)
       CALL KPAR_SYNC_ALL(WDES,W)
       CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_CUR, IO%IU6, IO%IU0)
       CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
    ENDIF

  END SUBROUTINE PSI_CORR_MAIN
 
!*********************************************************************
!
! helper routine: apply scissor correction
!
!*********************************************************************

  SUBROUTINE APPLY_SCISSOR( W, DELTA)
    REAL(q) DELTA                   ! shift of conduction band
    TYPE (wavespin)      W          ! wavefunction array
  ! local
    INTEGER NK, N, ISP

    DO ISP=1,W%WDES%ISPIN
       DO NK=1,W%WDES%NKPTS
          DO N=1,W%WDES%NB_TOT
             IF (EMPTY_XI_ORBITAL( W%FERTOT(N, NK, ISP))) THEN
                W%CELTOT(N, NK, ISP)= W%CELTOT(N, NK, ISP)+DELTA
             ENDIF
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE APPLY_SCISSOR


!*********************************************************************
!
! set the frequencies in the CHI array that are currently
! considered 
! the data distribution in parallel mode is determined by WGW
! this is not quite a simple routine, but tested rather carefully
! the frequency points are distributed among the nodes starting
! at NOMEGA_INDEX, the maximum number of frequencies is restricted
! by the number of frequencies in CHI%OMEGA
!
!*********************************************************************

  SUBROUTINE SET_RESPONSEFUN_FREQ( CHI, WGW, COMEGA, NOMEGA_INDEX)
    USE constant
    USE wave 

    IMPLICIT NONE
    TYPE (responsefunction) CHI  
    TYPE (wavedes) WGW
    COMPLEX(q) COMEGA(:)            ! frequencies to be considered
    INTEGER NOMEGA_INDEX            ! first frequency in the array COMEGA to be considered
  ! local
    INTEGER NOMEGA
    
  ! number of frequencies, if COMEGA were distributed entirely onto all nodes
  ! starting at frequency NOMEGA_INDEX
    NOMEGA=( SIZE(COMEGA)-NOMEGA_INDEX+WGW%NB_PAR)/WGW%NB_PAR
  ! restrict by the size of CHI
    NOMEGA=MIN(NOMEGA, SIZE(CHI%OMEGA))

    CHI%NOMEGA_LOW =(WGW%NB_LOW-1)*NOMEGA+NOMEGA_INDEX
    CHI%NOMEGA_HIGH= WGW%NB_LOW   *NOMEGA+NOMEGA_INDEX-1

  ! NOMEGA_LOW restricted by  size of OMEGA
    CHI%NOMEGA_LOW=MIN(CHI%NOMEGA_LOW, SIZE(COMEGA)+1)

  ! NOMEGA_HIGH restricted by  size of COMEGA
    CHI%NOMEGA_HIGH=MIN(CHI%NOMEGA_HIGH, SIZE(COMEGA))

    CHI%NOMEGA=CHI%NOMEGA_HIGH-CHI%NOMEGA_LOW+1

    CHI%COMEGA(1:CHI%NOMEGA)=COMEGA(CHI%NOMEGA_LOW:CHI%NOMEGA_HIGH)
    CHI%OMEGA(1:CHI%NOMEGA) =ABS(COMEGA(CHI%NOMEGA_LOW:CHI%NOMEGA_HIGH))

!    WRITE(*,*) 'data distribution:',WGW%NB_LOW, CHI%NOMEGA_LOW, CHI%NOMEGA_HIGH, CHI%NOMEGA
    
  END SUBROUTINE SET_RESPONSEFUN_FREQ

!
! this version is applied when a parallelization over the inner
! k-point loop is performed
! essentially identical to the above routine, however uses entire WGW%COMM
! to determine data distribution
!

  SUBROUTINE SET_RESPONSEFUN_FREQ_KSUM( CHI, WGW, COMEGA, NOMEGA_INDEX)
    USE constant
    USE wave 

    IMPLICIT NONE
    TYPE (responsefunction) CHI  
    TYPE (wavedes) WGW
    COMPLEX(q) COMEGA(:)            ! frequencies to be considered
    INTEGER NOMEGA_INDEX            ! first frequency in the array COMEGA to be considered
  ! local
    INTEGER NOMEGA, NOMEGA_TOTAL, NCPU_STRIDE
    INTEGER NCPU, NODE_ME
#ifdef MPI    
    NCPU   =WGW%COMM%NCPU
    NODE_ME=WGW%COMM%NODE_ME
#else
    NCPU=1
    NODE_ME=1
#endif
  ! number of frequencies, if COMEGA were distributed entirely onto all nodes
  ! starting at frequency NOMEGA_INDEX
    NOMEGA=( SIZE(COMEGA)-NOMEGA_INDEX+NCPU)/NCPU

    NOMEGA_TOTAL=SIZE(COMEGA)-NOMEGA_INDEX+1
  ! restrict by the size of CHI
    NOMEGA=MIN(NOMEGA, SIZE(CHI%OMEGA))

  ! conventional round robin
    CHI%NOMEGA_LOW =(NODE_ME-1)*NOMEGA+NOMEGA_INDEX
    CHI%NOMEGA_HIGH= NODE_ME   *NOMEGA+NOMEGA_INDEX-1
  ! however, if the number of cores is larger than twice the total number
  ! of frequencies NOMEGA_TOTAL, 
  ! then only every NCPU/NOMEGA_TOTAL core receives data
  ! this results in better load/storage balancing on multi core CPUs
    IF (NOMEGA==1 .AND. NOMEGA_TOTAL*2 <= NCPU) THEN
       NCPU_STRIDE=NCPU/NOMEGA_TOTAL
       IF (MOD(NODE_ME-1,NCPU_STRIDE)==0) THEN
          CHI%NOMEGA_LOW =(NODE_ME-1)/NCPU_STRIDE+NOMEGA_INDEX
          CHI%NOMEGA_HIGH=CHI%NOMEGA_LOW
       ELSE
          ! zero size array
          CHI%NOMEGA_LOW =SIZE(COMEGA)+1
          CHI%NOMEGA_HIGH=SIZE(COMEGA)
       ENDIF
    ENDIF

  ! NOMEGA_LOW restricted by  size of OMEGA
    CHI%NOMEGA_LOW=MIN(CHI%NOMEGA_LOW, SIZE(COMEGA)+1)

  ! NOMEGA_HIGH restricted by  size of COMEGA
    CHI%NOMEGA_HIGH=MIN(CHI%NOMEGA_HIGH, SIZE(COMEGA))

    CHI%NOMEGA=CHI%NOMEGA_HIGH-CHI%NOMEGA_LOW+1

    CHI%COMEGA(1:CHI%NOMEGA)=COMEGA(CHI%NOMEGA_LOW:CHI%NOMEGA_HIGH)
    CHI%OMEGA (1:CHI%NOMEGA)=ABS(COMEGA(CHI%NOMEGA_LOW:CHI%NOMEGA_HIGH))

!    WRITE(*,*) 'data distribution:', NODE_ME, CHI%NOMEGA_LOW, CHI%NOMEGA_HIGH, CHI%NOMEGA
    
  END SUBROUTINE SET_RESPONSEFUN_FREQ_KSUM


!*********************************************************************
!
! set the frequencies in the WPOT array that should be considered
!
! the data distribution in parallel mode is determined by WGW
! the function returns true, when all frequencies have been done
! NOMEGA_INDEX is the index, that is incremented by the routine
! and this routines decides how NOMEGA_INDEX is to be interpreted
! for LSPECT:
!   the frequencies in WPOT run over all frequencies (COMEGA)
! whereas for .NOT. LSPECT
!   the frequencies in WPOT run from NOMEGA_INDEX_START to NOMEGA_INDEX_END
!
!*********************************************************************

  FUNCTION SET_RESPONSEFUN_FREQI( WPOT, WGW, COMEGA, NOMEGA_INDEX, & 
       NOMEGA_INDEX_START, NOMEGA_INDEX_END, LSPECT)
    USE constant
    USE wave 

    IMPLICIT NONE
    LOGICAL SET_RESPONSEFUN_FREQI
    TYPE (responsefunction) CHI, WPOT
    TYPE (wavedes) WGW
    COMPLEX(q) COMEGA(:)            ! frequencies to be considered
    INTEGER NOMEGA_INDEX            ! index
    INTEGER NOMEGA_INDEX_START      ! initial index COMEGA (for LSPECT)
    INTEGER NOMEGA_INDEX_END        ! final index into COMEGA (for LSPECT)
    LOGICAL LSPECT                  ! spectral representation used or not
  ! local
    INTEGER NOMEGA
    
    IF (LSPECT) THEN 
       !
       ! do all frequencies in OMEGA starting from 1 until SIZE(COMEGA)
       ! in steps of SIZE(WPOT%COMEGA)
       !
       NOMEGA=SIZE(COMEGA)
       IF (NOMEGA_INDEX+1 > NOMEGA ) THEN
          ! all frequencies done finish
          SET_RESPONSEFUN_FREQI=.TRUE.
          WPOT%NOMEGA=0
       ELSE
          ! continue loop
          SET_RESPONSEFUN_FREQI=.FALSE.

          WPOT%NOMEGA_LOW =NOMEGA_INDEX+1
          WPOT%NOMEGA_HIGH=NOMEGA_INDEX+SIZE(WPOT%COMEGA)

          ! restrict index to be below NOMEGA
          WPOT%NOMEGA_HIGH=MIN(WPOT%NOMEGA_HIGH, NOMEGA)

          ! increase NOMEGA_INDEX
          NOMEGA_INDEX=NOMEGA_INDEX+SIZE(WPOT%COMEGA)

          WPOT%NOMEGA=WPOT%NOMEGA_HIGH-WPOT%NOMEGA_LOW+1
          WPOT%COMEGA(1:WPOT%NOMEGA)=COMEGA(WPOT%NOMEGA_LOW:WPOT%NOMEGA_HIGH)
          WPOT%OMEGA(1:WPOT%NOMEGA) =ABS(COMEGA(WPOT%NOMEGA_LOW:WPOT%NOMEGA_HIGH))
       ENDIF
    ELSE
       !
       ! do all frequencies in COMEGA starting with NOMEGA_INDEX_START until NOMEGA_INDEX_END
       ! in steps of SIZE(WPOT%COMEGA)
       !
       NOMEGA=MIN(SIZE(COMEGA), NOMEGA_INDEX_END)
       IF (NOMEGA_INDEX+NOMEGA_INDEX_START > NOMEGA) THEN
          ! all frequencies done finish
          SET_RESPONSEFUN_FREQI=.TRUE.
          WPOT%NOMEGA=0
       ELSE
          ! continue loop
          SET_RESPONSEFUN_FREQI=.FALSE.

          WPOT%NOMEGA_LOW =NOMEGA_INDEX+NOMEGA_INDEX_START
          WPOT%NOMEGA_HIGH=NOMEGA_INDEX+SIZE(WPOT%COMEGA)+NOMEGA_INDEX_START-1

          ! restrict index to be below NOMEGA
          WPOT%NOMEGA_HIGH=MIN(WPOT%NOMEGA_HIGH, NOMEGA)

          ! increase NOMEGA_INDEX
          NOMEGA_INDEX=NOMEGA_INDEX+SIZE(WPOT%COMEGA)

          WPOT%NOMEGA=WPOT%NOMEGA_HIGH-WPOT%NOMEGA_LOW+1
          WPOT%COMEGA(1:WPOT%NOMEGA)=COMEGA(WPOT%NOMEGA_LOW:WPOT%NOMEGA_HIGH)
          WPOT%OMEGA(1:WPOT%NOMEGA) =ABS(COMEGA(WPOT%NOMEGA_LOW:WPOT%NOMEGA_HIGH))
       ENDIF
    ENDIF

!    WRITE(*,*) 'data distribution:',WGW%NB_LOW, WPOT%NOMEGA_LOW, WPOT%NOMEGA_HIGH, WPOT%NOMEGA
!    WRITE(*,'(A,I3,(10F14.7))') ' data distribution:',WGW%NB_LOW, REAL(WPOT%COMEGA(1:WPOT%NOMEGA),q)
    
  END FUNCTION SET_RESPONSEFUN_FREQI


!**********************************************************************
!
! set the DATAKE entry of the WGW array to the Coulomb kernel v
! instead of the kinetic energy
!
! also reallocate the DATAKE array for the GAMMA only case
! double its size, since the complex coefficients C_G are
! interpreted as cosine and sin transforms and are stored seperately
! admittedly this is really akward, but difficult to program in an
! entirely clean fashion
!
!**********************************************************************

  SUBROUTINE XI_DATAKE( WGW, LATT_CUR, ENCUT, ENCUTSOFT)
    USE fock
    USE constant
    USE acfdt, ONLY: SQUEEZED_COULOMB_KERNEL
    TYPE (wavedes) WGW
    TYPE (latt) LATT_CUR
    REAL(q), OPTIONAL :: ENCUT, ENCUTSOFT
  ! local
    TYPE (wavedes1) :: WGWQ
    INTEGER    NP, NI, NI_, NQ
    REAL(q) :: DKX, DKY, DKZ, GX, GY, GZ, GSQU, SCALE, POTFAK, E
    REAL(q) :: SFUN, DFUN
    REAL(q) :: Q1, Q2, ALPHA, BETA

    ! e^2/ volume 
    ! coupling constant is lambda is multiplied in here as well
    SCALE=EDEPS/LATT_CUR%OMEGA * LAMBDA

    IF (WGW%LGAMMA) THEN
       DEALLOCATE(WGW%DATAKE)
       ALLOCATE(WGW%DATAKE(WGW%NGDIM*2,2,WGW%NKPTS))
    ENDIF

    IF (PRESENT(ENCUT).AND. PRESENT(ENCUTSOFT)) THEN
       Q1=SQRT(ENCUTSOFT/HSQDTM)
       Q2=SQRT(ENCUT/HSQDTM)

!      use this line if the  ATTENUATE_CUTOFF_LINEAR kernel is used to determine ALPHA
!       CALL CONSERVING_KERNEL(Q1, Q2, ALPHA )
    ENDIF

    DO NQ=1,WGW%NKPTS
    CALL SETWDES(WGW, WGWQ, NQ )

    DKX=(WGW%VKPT(1,NQ))*LATT_CUR%B(1,1)+ &
        (WGW%VKPT(2,NQ))*LATT_CUR%B(1,2)+ &
        (WGW%VKPT(3,NQ))*LATT_CUR%B(1,3)
    DKY=(WGW%VKPT(1,NQ))*LATT_CUR%B(2,1)+ &
        (WGW%VKPT(2,NQ))*LATT_CUR%B(2,2)+ &
        (WGW%VKPT(3,NQ))*LATT_CUR%B(2,3)
    DKZ=(WGW%VKPT(1,NQ))*LATT_CUR%B(3,1)+ &
        (WGW%VKPT(2,NQ))*LATT_CUR%B(3,2)+ &
        (WGW%VKPT(3,NQ))*LATT_CUR%B(3,3)

    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2

    DO NI=1,NP
       NI_=NI
       IF (WGWQ%LGAMMA) NI_=(NI-1)/2+1
       
       GX=(WGWQ%IGX(NI_)*LATT_CUR%B(1,1)+WGWQ%IGY(NI_)* &
            LATT_CUR%B(1,2)+WGWQ%IGZ(NI_)*LATT_CUR%B(1,3))
       GY=(WGWQ%IGX(NI_)*LATT_CUR%B(2,1)+WGWQ%IGY(NI_)* &
            LATT_CUR%B(2,2)+WGWQ%IGZ(NI_)*LATT_CUR%B(2,3))
       GZ=(WGWQ%IGX(NI_)*LATT_CUR%B(3,1)+WGWQ%IGY(NI_)* &
            LATT_CUR%B(3,2)+WGWQ%IGZ(NI_)*LATT_CUR%B(3,3))
          
       GSQU=(DKX+GX)**2+(DKY+GY)**2+(DKZ+GZ)**2
       
       IF (ABS(GSQU)<G2ZERO) THEN
          ! head and wing
          POTFAK=SCALE
          IF (HFRCUT/=0) THEN
             ! spherical cutoff on Coloumb kernel
             ! no regular contribution since head stores q^2 like terms
             POTFAK=0
          ENDIF
       ELSE
          ! the factor 1/(2 pi)^2 is required to obtain proper reciprocal
          ! lattice vector lenght
          POTFAK=SCALE/(GSQU*TPI**2)
          IF (HFRCUT/=0) THEN
             ! spherical cutoff on Coloumb kernel
             ! see for instance C.A. Rozzi, PRB 73, 205119 (2006)
             ! POTFAK=POTFAK*(1-COS(SQRT(GSQU)*TPI*HFRCUT)*EXP(-(SQRT(GSQU)*TPI*HFRCUT)**2*HFRCUT_SMOOTH))
             ! well unfortunately this does not work at all for the correlation energy
             ! since the large G oscillations spoil the convergence
             ! a smooth cutoff however seems to work
             CALL DELSTP(3,SQRT(GSQU)*TPI*HFRCUT/10,DFUN,SFUN)
             POTFAK=POTFAK*(SFUN-0.5)*2
          ELSE IF (LRSCOR) THEN
             ! Martijn asked me to support long range HF
             POTFAK=POTFAK*EXP(-GSQU*(TPI*TPI/(4*HFSCREENC*HFSCREENC)))
          ENDIF
       ENDIF

       ! smooth cutoff function between  ENCUTSOFT and ENCUT
       E=HSQDTM*(GSQU*TPI**2)
       IF (PRESENT(ENCUT).AND. PRESENT(ENCUTSOFT)) THEN 
          IF (E>ENCUTSOFT) THEN
             IF ( .NOT. LSCK) THEN
                POTFAK=POTFAK*(1+COS((E-ENCUTSOFT)/(ENCUT-ENCUTSOFT)*PI))/2
             ELSE
                POTFAK=POTFAK*SQUEEZED_COULOMB_KERNEL(SQRT(GSQU)*TPI, Q1, Q2)
             ENDIF
          ENDIF
       ENDIF

       WGW%DATAKE(NI, 1, NQ )=POTFAK
       WGW%DATAKE(NI, 2, NQ )=POTFAK
    ENDDO
    ENDDO

  END SUBROUTINE XI_DATAKE

!********************** CONSERVING_KERNEL ****************************
!
! this small subroutine constructs a "conserving" kernel
!  int_q1^infty dq 1/q^4 = int_q1^q2  dq 1/q^4 f(q)^2 
! where 
! f(q) = sum_i alpha_i 1/2 (1+cos( (2 i+1) pi (q-q1)/(q2-q1) )
! 
! the two coefficients alpha_i are returned to the calling routine
!
!**********************************************************************

  SUBROUTINE CONSERVING_KERNEL_TEST(Q1, Q2, ALPHA, BETA)
    USE acfdt, ONLY: SQUEEZED_COULOMB_KERNEL
    IMPLICIT NONE
    REAL(q) :: Q1, Q2, ALPHA, BETA

    INTEGER, PARAMETER :: NSTEP=400
   
    ! local variables
    INTEGER :: I, N
    REAL(q) :: H, X, A, B 

    ! integrate two functions between Q1 and Q2
    H=(Q2-Q1)/NSTEP

    ! loop over ALPHA until we get approximately the right norm
    ! alas very very quick hack
    A=1/Q1**4*SQUEEZED_COULOMB_KERNEL(Q1, Q1, Q2)**2 + & 
         1/Q2**4*SQUEEZED_COULOMB_KERNEL(Q2, Q1, Q2)**2
    DO N=1,NSTEP-1
       X=Q1+H*N
       IF (MOD(N,2)==1) THEN
          A=A+1/X**4*SQUEEZED_COULOMB_KERNEL(X, Q1, Q2)**2*4
       ELSE
          A=A+1/X**4*SQUEEZED_COULOMB_KERNEL(X, Q1, Q2)**2*2
       ENDIF
    ENDDO
    A=A*H/3
    B=A
    
    A=1/Q1**6*SQUEEZED_COULOMB_KERNEL(Q1, Q1, Q2)**2 + & 
         1/Q2**6*SQUEEZED_COULOMB_KERNEL(Q2, Q1, Q2)**2
    DO N=1,NSTEP-1
       X=Q1+H*N
       IF (MOD(N,2)==1) THEN
          A=A+1/X**6*SQUEEZED_COULOMB_KERNEL(X, Q1, Q2)**2*4
       ELSE
          A=A+1/X**6*SQUEEZED_COULOMB_KERNEL(X, Q1, Q2)**2*2
       ENDIF
    ENDDO
    A=A*H/3
    WRITE(*,'(A,2F10.5,A,2E14.5,A,2E14.5)') 'ALPHA',ALPHA,BETA, ' Q4',B, 1/Q1**3/3, ' Q6',A, 1/Q1**5/5
  END SUBROUTINE CONSERVING_KERNEL_TEST

!
! this subroutine determines ALPHA for ATTENUATE_CUTOFF_LINEAR
! by stupid search
! TODO this should be done by intervall bisectioning
!
  SUBROUTINE CONSERVING_KERNEL(Q1, Q2, ALPHA)
    IMPLICIT NONE
    REAL(q) :: Q1, Q2, ALPHA

    INTEGER, PARAMETER :: NSTEP=400
   
    ! local variables
    INTEGER :: I, N
    REAL(q) :: H, X, A

    ! integrate two functions between Q1 and Q2
    H=(Q2-Q1)/NSTEP

    ! loop over ALPHA until we get approximately the right norm
    ! alas very very quick hack
    DO I=1,5000
       ALPHA=I*0.001 
       ! loop over all allowed functions
       ! Simpson to calculate integral
       A=1/Q1**4*ATTENUATE_CUTOFF_QUAD(Q1, Q1, Q2, ALPHA)**2 + & 
         1/Q2**4*ATTENUATE_CUTOFF_QUAD(Q2, Q1, Q2, ALPHA)**2
       DO N=1,NSTEP-1
          X=Q1+H*N
          IF (MOD(N,2)==1) THEN
             A=A+1/X**4*ATTENUATE_CUTOFF_QUAD(X, Q1, Q2, ALPHA)**2*4
          ELSE
             A=A+1/X**4*ATTENUATE_CUTOFF_QUAD(X, Q1, Q2, ALPHA)**2*2
          ENDIF
       ENDDO
       A=A*H/3

       IF (A>1/Q1**3/3) EXIT
    ENDDO
  END SUBROUTINE CONSERVING_KERNEL


!
! this subroutine determine ALPHA and BETA for ATTENUATE_CUTOFF_QUADRATIC
! by stupid search
! TODO this should be done by intervall bisectioning
!
    
  SUBROUTINE CONSERVING_KERNEL_QUADRATIC(Q1, Q2, ALPHA, BETA)
    IMPLICIT NONE
    REAL(q) :: Q1, Q2, ALPHA, BETA

    INTEGER, PARAMETER :: NSTEP=400

    ! local variables
    INTEGER :: I, J, N
    REAL(q) :: H, X, A, B, ERROR, EMIN, ALPHAMIN, BETAMIN

    ! integrate two functions between Q1 and Q2
    H=(Q2-Q1)/NSTEP
    EMIN = 1000

    ! loop over ALPHA until we get approximately the right norm
    ! alas very very quick hack
    DO J=-3000,-2000
       BETA=J*0.001
    DO I=-3000,-2000
       ALPHA=I*0.001 
       ! loop over all allowed functions
       ! Simpson to calculate integral
       A=1/Q1**4*ATTENUATE_CUTOFF_QUADRATIC(Q1, Q1, Q2, ALPHA,BETA)**2 + & 
         1/Q2**4*ATTENUATE_CUTOFF_QUADRATIC(Q2, Q1, Q2, ALPHA,BETA)**2
       DO N=1,NSTEP-1
          X=Q1+H*N
          IF (MOD(N,2)==1) THEN
             A=A+1/X**4*ATTENUATE_CUTOFF_QUADRATIC(X, Q1, Q2, ALPHA,BETA)**2*4
          ELSE
             A=A+1/X**4*ATTENUATE_CUTOFF_QUADRATIC(X, Q1, Q2, ALPHA,BETA)**2*2
          ENDIF
       ENDDO
       A=A*H/3
       B=A

       A=1/Q1**6*ATTENUATE_CUTOFF_QUADRATIC(Q1, Q1, Q2, ALPHA,BETA)**2 + & 
         1/Q2**6*ATTENUATE_CUTOFF_QUADRATIC(Q2, Q1, Q2, ALPHA,BETA)**2
       DO N=1,NSTEP-1
          X=Q1+H*N
          IF (MOD(N,2)==1) THEN
             A=A+1/X**6*ATTENUATE_CUTOFF_QUADRATIC(X, Q1, Q2, ALPHA,BETA)**2*4
          ELSE
             A=A+1/X**6*ATTENUATE_CUTOFF_QUADRATIC(X, Q1, Q2, ALPHA,BETA)**2*2
          ENDIF
       ENDDO
       A=A*H/3
       ERROR=((B- 1/Q1**3/3)/(1/Q1**3/3))**2*10+ ((A- 1/Q1**5/5)/(1/Q1**5/5))**2
!       WRITE(*,'(A,3F10.5,A,1E14.5,A,1E14.5)') 'KFERMI',ALPHA,BETA,ERROR, ' Q4',B- 1/Q1**3/3, ' Q6',A- 1/Q1**5/5
       IF ( ERROR <EMIN) THEN
!         WRITE(*,'(A,3F10.5,A,1E14.5,A,1E14.5)') 'KFERMI',ALPHA,BETA,ERROR, ' Q4',B- 1/Q1**3/3, ' Q6',A- 1/Q1**5/5
         EMIN=ERROR
         ALPHAMIN=ALPHA
         BETAMIN=BETA
       ENDIF
    ENDDO
    ENDDO
    WRITE(*,*) ALPHAMIN, BETAMIN, EMIN
    ALPHA=ALPHAMIN
    BETA =BETAMIN
  END SUBROUTINE CONSERVING_KERNEL_QUADRATIC

!
! increase Coloumb potential linearly starting at Q1
!
  FUNCTION ATTENUATE_CUTOFF_LINEAR( X, Q1, Q2, ALPHA)
    USE prec
    USE constant
    IMPLICIT NONE

    REAL(q) :: ATTENUATE_CUTOFF_LINEAR, X, Q1, Q2, ALPHA

    ATTENUATE_CUTOFF_LINEAR=1._q+(X-Q1)/(Q2-Q1)*ALPHA
        
  END FUNCTION ATTENUATE_CUTOFF_LINEAR

!
! increase Coloumb potential linearly starting at Q1
!
  FUNCTION ATTENUATE_CUTOFF_QUAD( X, Q1, Q2, ALPHA)
    USE prec
    USE constant
    IMPLICIT NONE

    REAL(q) :: ATTENUATE_CUTOFF_QUAD, X, Q1, Q2, ALPHA

    ATTENUATE_CUTOFF_QUAD=1._q+(X-Q1)**2/(Q2-Q1)**2*ALPHA
        
  END FUNCTION ATTENUATE_CUTOFF_QUAD


  FUNCTION ATTENUATE_CUTOFF_QUADRATIC( X, Q1, Q2, ALPHA, BETA)
    USE prec
    USE constant
    IMPLICIT NONE

    REAL(q) :: ATTENUATE_CUTOFF_QUADRATIC, X, Q1, Q2, ALPHA, BETA

    ATTENUATE_CUTOFF_QUADRATIC=1._q+ (X-Q1)/(Q2-Q1)*ALPHA+ ((X-Q1)/(Q2-Q1))**2*BETA

  END FUNCTION ATTENUATE_CUTOFF_QUADRATIC

!
! increase Coloumb potential and than bring it back to 0 smoothly
! second version by Merzuk Kaltak
! not used anymore
!
  FUNCTION ATTENUATE_CUTOFF_SMOOTH2( X, Q0, Q1, ALPHA)
    USE prec
    USE constant
    IMPLICIT NONE

    REAL(q) :: ATTENUATE_CUTOFF_SMOOTH2, X, Q0, Q1, ALPHA

    IF (ABS(Q1-X)>1E-5_q) THEN
       ATTENUATE_CUTOFF_SMOOTH2= & 
         ((-Q0 + Q1)*X*(-Q1 + X)**2*Sqrt(( ALPHA + (Q0**2 + (-2*Q0 + Q1)*X)**2/(-Q1 + X)**2)/( ALPHA + X**2)))/&
               (Q0**2 + (-2*Q0 + Q1)*X)**3
       ATTENUATE_CUTOFF_SMOOTH2=ATTENUATE_CUTOFF_SMOOTH2* X**2
    ELSE
       ATTENUATE_CUTOFF_SMOOTH2=1
    ENDIF


  END FUNCTION ATTENUATE_CUTOFF_SMOOTH2


!********************** RESTORE_HEAD *********************************
!
! restore directions in the head of a dielectric matrix 
! which are not explicitly updated
!
!**********************************************************************

  SUBROUTINE RESTORE_HEAD( CHI )
    IMPLICIT NONE
    TYPE (responsefunction) CHI
  ! local
    INTEGER NOMEGA, IDIR

    IF (CHI%LGAMMA) THEN
       DO NOMEGA=1,CHI%NOMEGA
          DO IDIR=IDIR_MAX+1,3
             CHI%HEAD(IDIR,IDIR,NOMEGA)=CHI%HEAD(IDIR_MAX,IDIR_MAX,NOMEGA)
          ENDDO
       ENDDO
    END IF

  END SUBROUTINE RESTORE_HEAD

!********************** SUBROUTINE XI_TO_EPS *************************
!
! calculate the symmetric dielectric matrix in the RPA approximation 
!        1/2      1/2
!  1 - v     Xi v 
!
! if LPLUS is present the sign is inverted i.e.
!        1/2      1/2
!  1 + v     Xi v 
!
! is calculated (useful for the reducible polarizability)
!
! what is actually coded is
!
!  1- 4 pi e^2/ (2 pi)^2 |G+q| |G'+q'| RESPONSEFUN(G, G')
!
! this follows from the fact that the indices G' have been obtained
! from 
!   1/N \sum_r' u*_k+q(r') u_k(r') e iG'r'
! therefore the polarizability (including the Bloch factor) is given by
!
!   chi(r',G) = sum_G'  RESPONSEFUN(G', G) e i(-G'-q)r'
!
! likewise the second index is related to quantities
!
!   1/N \sum_r u*_k(r) u_k+q(r) e -iGr
!
! and the polarizability in real space is hence given by
!
!   chi(G',r) = sum_G   RESPONSEFUN(G', G) e i(G+q)r
! 
!**********************************************************************

  SUBROUTINE XI_TO_EPS( CHI, WGWQ, LPLUS)
    USE prec
    USE constant
    USE mgrid
    USE wave 
    USE full_kpoints

    IMPLICIT NONE
    TYPE (responsefunction) CHI
    TYPE (wavedes1) WGWQ
    LOGICAL, OPTIONAL :: LPLUS
  ! local
    INTEGER    NI, NOMEGA, NP
    REAL(q)    POTFAK, SCALE
    
    SCALE=-1
    IF (PRESENT(LPLUS)) THEN
       SCALE=1
    ELSE
       SCALE=-1
    ENDIF
    CHI%RESPONSEFUN=CHI%RESPONSEFUN*SCALE

    IF (CHI%LGAMMA) THEN
       CHI%HEAD(:,:,1:CHI%NOMEGA) = CHI%HEAD(:,:,1:CHI%NOMEGA)*SCALE*WGWQ%DATAKE(1,1)
       CHI%WING(:,:,1:CHI%NOMEGA) = CHI%WING(:,:,1:CHI%NOMEGA)*SCALE*SQRT(WGWQ%DATAKE(1,1))
       CHI%CWING(:,:,1:CHI%NOMEGA)=CHI%CWING(:,:,1:CHI%NOMEGA)*SCALE*SQRT(WGWQ%DATAKE(1,1))
       
       DO NOMEGA=1,CHI%NOMEGA
          CHI%HEAD(1,1,NOMEGA)= 1+CHI%HEAD(1,1,NOMEGA)
          CHI%HEAD(2,2,NOMEGA)= 1+CHI%HEAD(2,2,NOMEGA)
          CHI%HEAD(3,3,NOMEGA)= 1+CHI%HEAD(3,3,NOMEGA)
       ENDDO
    ENDIF

    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2
    DO NI=1,NP
       POTFAK=SQRT(WGWQ%DATAKE(NI,1))

       DO NOMEGA=1,CHI%NOMEGA
          CHI%RESPONSEFUN(1:NP,NI,NOMEGA)=CHI%RESPONSEFUN(1:NP,NI,NOMEGA)*POTFAK
          CHI%RESPONSEFUN(NI,1:NP,NOMEGA)=CHI%RESPONSEFUN(NI,1:NP,NOMEGA)*POTFAK
          CHI%RESPONSEFUN(NI,NI,NOMEGA)=1+CHI%RESPONSEFUN(NI,NI,NOMEGA)
          IF (CHI%LGAMMA) THEN
             CHI%WING(NI,:,NOMEGA) =CHI%WING(NI,:,NOMEGA)*POTFAK
             CHI%CWING(NI,:,NOMEGA)=CHI%CWING(NI,:,NOMEGA)*POTFAK
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE XI_TO_EPS


!********************** SUBROUTINE XI_TO_W_SYM ************************
!
! calculate potential operator times inverse symmetric dielectric matrix
! minus the bare Coulomb operator            ---------
! (if LFOCK_SUBTRACT is set the bare Coulomb operator is subtracted)
!
!      1/2        -1    1/2
!    v     epsilon    v       -  v
!
!    4 pi e^2/ (2 pi)^2 |G+q| |G'+q'| RESPONSEFUN(G, G')
!  - 4 pi e^2/ (2 pi)^2 |G+q|^2
!
! 
!**********************************************************************

  SUBROUTINE XI_TO_W_SYM( CHI, WGWQ )
    USE prec
    USE constant
    USE mgrid
    USE wave 
    USE full_kpoints
    USE fock

    IMPLICIT NONE
    TYPE (responsefunction) CHI
    TYPE (wavedes1) WGWQ
    LOGICAL LFOCK_SUBTRACT
    ! local
    INTEGER    NI, NC, NQ_IN_WGW, NOMEGA, NP
    REAL(q) :: POTFAK

    IF (CHI%LGAMMA) THEN
       CHI%WING(:,:,1:CHI%NOMEGA) = CHI%WING(:,:,1:CHI%NOMEGA)*SQRT(WGWQ%DATAKE(1,1))
       CHI%CWING(:,:,1:CHI%NOMEGA)=CHI%CWING(:,:,1:CHI%NOMEGA)*SQRT(WGWQ%DATAKE(1,1))
    ENDIF

    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2

    DO NI=1,NP
       POTFAK=SQRT(WGWQ%DATAKE(NI,1))

       DO NOMEGA=1,CHI%NOMEGA
          CHI%RESPONSEFUN(1:NP,NI,NOMEGA)=CHI%RESPONSEFUN(1:NP,NI,NOMEGA)*POTFAK
          CHI%RESPONSEFUN(NI,1:NP,NOMEGA)=CHI%RESPONSEFUN(NI,1:NP,NOMEGA)*POTFAK
          IF (CHI%LGAMMA) THEN
             CHI%WING(NI,:,NOMEGA) =CHI%WING(NI,:,NOMEGA)*POTFAK
             CHI%CWING(NI,:,NOMEGA)=CHI%CWING(NI,:,NOMEGA)*POTFAK
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE XI_TO_W_SYM

!********************** SUBROUTINE XI_TO_W ****************************
!
! calculate potential operator times inverse dielectric matrix
! minus the bare Coulomb operator
! (if LFOCK_SUBTRACT is set the bare Coulomb operator is subtracted)
!
!                 -1       
!          epsilon    v
!
! the actual code calculates
!    RESPONSEFUN(G, G') 4 pi e^2/ (2 pi)^2 |G'+q'|^2
!  - 4 pi e^2/ (2 pi)^2 |G+q|^2
! 
!**********************************************************************

  SUBROUTINE XI_TO_W( CHI, WGWQ)
    USE constant
    USE wave 

    IMPLICIT NONE
    TYPE (responsefunction) CHI
    TYPE (wavedes1) WGWQ
    LOGICAL LFOCK_SUBTRACT
    ! local
    INTEGER    NOMEGA

    ! Xi v -> Xi
    DO NOMEGA=1,CHI%NOMEGA
       CALL XI_HARTREE_T( CHI, NOMEGA, WGWQ )
    ENDDO

  END SUBROUTINE XI_TO_W


!********************** SUBROUTINE W_SPECTRAL *************************
!
! calculate the spectral function of W
! this involves taking the anti-Hermitian part of the matrix W
! a each frequency point and multiplying with -i
! 
!**********************************************************************

  SUBROUTINE W_SPECTRAL( CHI, WGWQ)
    USE constant
    USE wave 

    IMPLICIT NONE
    TYPE (responsefunction) CHI
    TYPE (wavedes1) WGWQ
    LOGICAL LFOCK_SUBTRACT
    ! local
    INTEGER    NP, NI, NIP, NQ_IN_WGW, NOMEGA

    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2

    DO NOMEGA=1,CHI%NOMEGA
       DO NI=1,NP
          DO NIP=NI,NP
             CHI%RESPONSEFUN(NI,NIP,NOMEGA)=-(CHI%RESPONSEFUN(NI,NIP,NOMEGA)-CONJG(CHI%RESPONSEFUN(NIP,NI,NOMEGA)))/2*(0._q,1._q)
             CHI%RESPONSEFUN(NIP,NI,NOMEGA)=CONJG(CHI%RESPONSEFUN(NI,NIP,NOMEGA))
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE W_SPECTRAL

!********************** SUBROUTINE XI_COULOMB *************************
!
! set the convergence corrections
! limitations: presently the implementation works
! for diagonal dielectric matrices only
! ) hence it is restricted to orthorhombic systems
! ) and hexagonal systems
! 
!**********************************************************************

  SUBROUTINE XI_COULOMB( CHI, WGWQ, LATT_CUR, NQ_, FSG0, LFOCK_SUBTRACT)
    USE constant
    USE lattice
    USE wave 
    USE kpoints_change

    IMPLICIT NONE
    TYPE (responsefunction) CHI
    TYPE (wavedes1) WGWQ
    TYPE (latt) LATT_CUR
    INTEGER :: NQ_         ! total number of q-points used (differences vectors in GW)
    REAL(q) :: FSG0        ! convergence correction for bare Coulomb kernel
    LOGICAL LFOCK_SUBTRACT
  ! local
    INTEGER NOMEGA, NI, IDIR, NP
    TYPE (latt) LATT_EWALD
    REAL(q) :: B(3,3)
    REAL(q) :: SCALE1, SCALE2, SCALE3, FSG
    COMPLEX(q) :: INVEPS
    INTEGER :: NQ          ! total number of q-points used (differences vectors in GW) 

    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2

    ! subtract bare Coulomb operator
    DO NOMEGA=1,CHI%NOMEGA
       IF (LFOCK_SUBTRACT) THEN
          DO NI=1,NP
             CHI%RESPONSEFUN(NI,NI,NOMEGA)=CHI%RESPONSEFUN(NI,NI,NOMEGA)-WGWQ%DATAKE(NI,1)
          ENDDO
       ENDIF
    ENDDO

    ! gamma point
    ! get proper convergence corrections
    ! for the 1/r long range behaviour lim G-> 4 pi e^2/G^2
    IF (CHI%LGAMMA) THEN
       ! clear wings in the body of the matrix
       DO NOMEGA=1,CHI%NOMEGA
          CHI%RESPONSEFUN(1:NP,1,NOMEGA)=0
          CHI%RESPONSEFUN(1,1:NP,NOMEGA)=0
       ENDDO

       ! convergence correction for unscreened bare Coulomb kernel
       ! must be passed back to calling routine
       IF( HFRCUT==0) THEN
          LATT_EWALD=LATT_CUR
          IF (KPOINTS_FULL%NKPX /=-1 .AND. KPOINTS_FULL%NKPY /=-1 .AND. KPOINTS_FULL%NKPZ /=-1) THEN
             LATT_EWALD%A(:,1)=LATT_EWALD%A(:,1)*KPOINTS_FULL%NKPX/NKREDX
             LATT_EWALD%A(:,2)=LATT_EWALD%A(:,2)*KPOINTS_FULL%NKPY/NKREDY
             LATT_EWALD%A(:,3)=LATT_EWALD%A(:,3)*KPOINTS_FULL%NKPZ/NKREDZ
          ELSE IF (MAXVAL(KPOINTS_FULL%B)>=0) THEN
             ! determine reciprocal lattice of the generating k-point mesh
             LATT_EWALD%A=KPOINTS_FULL%B
             CALL LATTIC(LATT_EWALD)
             ! store this as direct lattice
             LATT_EWALD%A=LATT_EWALD%B
             CALL LATTIC(LATT_EWALD)
          ELSE
             CALL vtutor%error("error in XI_COULOMB: presently the GW routine only supports k-point &
                &meshes \n generated automatically")
          ENDIF
             
          IF (ODDONLY .OR. EVENONLY) THEN
             B(1,:)=(/0.5_q,0.5_q,0.5_q/)
             B(2,:)=(/-.5_q,0.5_q,0.5_q/)
             B(3,:)=(/0.5_q,-.5_q,0.5_q/)
             LATT_EWALD%A(:,:)=MATMUL(B, LATT_EWALD%A(:,:))
          ENDIF
          
          CALL LATTIC(LATT_EWALD)
          CALL EWALD_MONO(FSG0,1.0_q,LATT_EWALD)
          ! number of q-points (difference vectors included in the GW)
          NQ=NINT(LATT_EWALD%OMEGA/LATT_CUR%OMEGA)
          IF (NQ /= NQ_) THEN
             CALL vtutor%bug("internal error in VASP: XI_COULOMB new method for calculating NQ is " &
                // "not correct", __FILE__, __LINE__)
          ENDIF
          FSG0=FSG0*2*NQ
       ELSE
          FSG0=HFRCUT*HFRCUT/2*EDEPS/LATT_CUR%OMEGA
       ENDIF
       ! scale convergence correction by coupling constant strenght
       FSG0=FSG0*LAMBDA

       DO NOMEGA=1,CHI%NOMEGA
        ! isotropic average of the  dielectric matrix
        INVEPS=(CHI%HEAD(1,1,NOMEGA)+CHI%HEAD(2,2,NOMEGA)+CHI%HEAD(3,3,NOMEGA))/3
        IF( HFRCUT==0) THEN
          ! the inverse dielectric tensor defines a metric
          ! which can be used to scale the Bravais lattice
          ! the Poisson equation generally reads
          !
          !   nabla  epsilon nabla V = 4 pi e^2 rho
          !
          ! using x' = epsilon-1/2 x
          ! the usual Poisson equation is obtained (laplace' V(x') = 4 pi e^2 rho(x'))
          ! this yields the following interaction between point charges
          ! in a non uniform dielectric medium
          !   1/ sqrt( x epsilon x ) det epsilon -1/2
          ! for an isotropic medium this reduces to 
          !   1/|x| epsilon-1
          
          LATT_EWALD=LATT_CUR

          ! determine scaling relations
          SCALE1=SQRT(ABS(CHI%HEAD(1,1,NOMEGA)/INVEPS))
          SCALE2=SQRT(ABS(CHI%HEAD(2,2,NOMEGA)/INVEPS))
          SCALE3=SQRT(ABS(CHI%HEAD(3,3,NOMEGA)/INVEPS))

          LATT_EWALD=LATT_CUR

          IF (KPOINTS_FULL%NKPX /=-1 .AND. KPOINTS_FULL%NKPY /=-1 .AND. KPOINTS_FULL%NKPZ /=-1) THEN
             LATT_EWALD%A(:,1)=LATT_EWALD%A(:,1)*KPOINTS_FULL%NKPX/NKREDX
             LATT_EWALD%A(:,2)=LATT_EWALD%A(:,2)*KPOINTS_FULL%NKPY/NKREDY
             LATT_EWALD%A(:,3)=LATT_EWALD%A(:,3)*KPOINTS_FULL%NKPZ/NKREDZ
          ELSE IF (MAXVAL(KPOINTS_FULL%B)>=0) THEN
             ! determine reciprocal lattice of the generating k-point mesh
             LATT_EWALD%A=KPOINTS_FULL%B
             CALL LATTIC(LATT_EWALD)
             ! store this as direct lattice
             LATT_EWALD%A=LATT_EWALD%B
             CALL LATTIC(LATT_EWALD)
          ELSE
             CALL vtutor%error("error in XI_COULOMB: presently the GW routine only supports k-point &
                &meshes \n generated automatically")
          ENDIF

          ! multiply by number of k-points in each direction and perform scaling
          LATT_EWALD%A(1,:)=LATT_EWALD%A(1,:)*SCALE1
          LATT_EWALD%A(2,:)=LATT_EWALD%A(2,:)*SCALE2
          LATT_EWALD%A(3,:)=LATT_EWALD%A(3,:)*SCALE3

          IF (ODDONLY .OR. EVENONLY) THEN
             B(1,:)=(/0.5_q,0.5_q,0.5_q/)
             B(2,:)=(/-.5_q,0.5_q,0.5_q/)
             B(3,:)=(/0.5_q,-.5_q,0.5_q/)
             LATT_EWALD%A(:,:)=MATMUL(B, LATT_EWALD%A(:,:))
          ENDIF

          CALL LATTIC(LATT_EWALD)
          CALL EWALD_MONO(FSG,1.0_q,LATT_EWALD)

          ! k-point weight needs to be removed since the k-weight is applied
          ! later when the sum over all q=k2-k1 is performed
          ! 1/weight = number of k-points
          ! also a factor two is lacking since we correct potentials and not energies
          ! finally  det epsilon-1/2 needs to be included
          FSG=FSG*2*NQ*SCALE1*SCALE2*SCALE3
        ELSE
          FSG=HFRCUT*HFRCUT/2*EDEPS/LATT_CUR%OMEGA
          ! matter of taste, actually for molecules
          ! the dielectric function is converging to 1
          ! no correction at G=0 is in principle required
          ! if a finite cutoff is used
          INVEPS=1
        ENDIF
        ! scale convergence correction by coupling constant strenght
        FSG=FSG*LAMBDA

        ! set head of matrix to proper convergence corrections
        ! convergence correction times isotropic average of epsilon
        IF (LFOCK_SUBTRACT) THEN
           ! subtracted  bare Coulomb correction if Fock kernel is evaluated
           ! seperately
           CHI%RESPONSEFUN(1,1,NOMEGA)=FSG*INVEPS-FSG0
        ELSE
           CHI%RESPONSEFUN(1,1,NOMEGA)=FSG*INVEPS
        ENDIF

        ! finally take bare Coulomb off from HEAD 
        IF (LFOCK_SUBTRACT) THEN
           CHI%HEAD(1,1,NOMEGA)=CHI%HEAD(1,1,NOMEGA)-1
           CHI%HEAD(2,2,NOMEGA)=CHI%HEAD(2,2,NOMEGA)-1
           CHI%HEAD(3,3,NOMEGA)=CHI%HEAD(3,3,NOMEGA)-1
        ENDIF
          
        ! for gamma only the second entry is the sine transform corresponding to G=0, always zero
        IF (CHI%LREAL) CHI%RESPONSEFUN(2,:,NOMEGA)=0
        IF (CHI%LREAL) CHI%RESPONSEFUN(:,2,NOMEGA)=0
       ENDDO
    ENDIF
  END SUBROUTINE XI_COULOMB


!********************** SUBROUTINE XI_LIND    *************************
!
! calculates contribution from Gamma point for metals 
! using the Lindhard model dielectric function
! 
!**********************************************************************

  SUBROUTINE XI_LIND   ( CHI, WGWQ, LATT_CUR, NQ )
    USE constant
    USE lattice
    USE wave 
    USE kpoints_change

    IMPLICIT NONE
    TYPE (responsefunction) CHI
    TYPE (wavedes1) WGWQ
    TYPE (latt) LATT_CUR
    INTEGER :: NQ          ! total number of q-points used
    REAL(q) :: FSG0        ! convergence correction for bare Coulomb kernel
    LOGICAL LFOCK_SUBTRACT
  ! local
    INTEGER NOMEGA, N1, N2, N3
    TYPE (latt) LATT_EWALD
    REAL(q) :: B(3,3), QQ(3), F(CHI%NOMEGA)
    REAL(q) :: SCALE1, SCALE2, SCALE3, FSG
    COMPLEX(q) :: INVEPS
    INTEGER, PARAMETER  :: NINTER=2

    ! gamma point
    ! get proper convergence corrections
    ! for the 1/r long range behaviour lim G-> 4 pi e^2/G^2
    IF (CHI%LGAMMA) THEN

       ! convergence correction for unscreened bare Coulomb kernel
       ! must be passed back to calling routine
       IF( HFRCUT==0) THEN
          LATT_EWALD=LATT_CUR
          IF (KPOINTS_FULL%NKPX /=-1 .AND. KPOINTS_FULL%NKPY /=-1 .AND. KPOINTS_FULL%NKPZ /=-1) THEN
             LATT_EWALD%A(:,1)=LATT_EWALD%A(:,1)*KPOINTS_FULL%NKPX/NKREDX
             LATT_EWALD%A(:,2)=LATT_EWALD%A(:,2)*KPOINTS_FULL%NKPY/NKREDY
             LATT_EWALD%A(:,3)=LATT_EWALD%A(:,3)*KPOINTS_FULL%NKPZ/NKREDZ
          ELSE IF (MAXVAL(KPOINTS_FULL%B)>=0) THEN
             ! determine reciprocal lattice of the generating k-point mesh
             LATT_EWALD%A=KPOINTS_FULL%B
             CALL LATTIC(LATT_EWALD)
             ! store this as direct lattice
             LATT_EWALD%A=LATT_EWALD%B
             CALL LATTIC(LATT_EWALD)
          ELSE
             CALL vtutor%error("error in XI_LIND: presently the GW routine only supports k-point &
                &meshes \n generated automatically")
          ENDIF
             
          IF (ODDONLY .OR. EVENONLY) THEN
             B(1,:)=(/0.5_q,0.5_q,0.5_q/)
             B(2,:)=(/-.5_q,0.5_q,0.5_q/)
             B(3,:)=(/0.5_q,-.5_q,0.5_q/)
             LATT_EWALD%A(:,:)=MATMUL(B, LATT_EWALD%A(:,:))
          ENDIF
          
          CALL LATTIC(LATT_EWALD)
          ! at this point LATT_EWALD%B(:,:) stores the generating lattice
          ! vectors of the k-point mesh
          
          F=0
          DO N1=0,NINTER-1
             DO N2=0,NINTER-1
                DO N3=0,NINTER-1
                   QQ= REAL(2*N1-NINTER+1,q)/(2*NINTER)*LATT_EWALD%B(:,1)+ &
                       REAL(2*N2-NINTER+1,q)/(2*NINTER)*LATT_EWALD%B(:,2)+ &
                       REAL(2*N3-NINTER+1,q)/(2*NINTER)*LATT_EWALD%B(:,3) 
                   WRITE(*,'(3F14.7)') QQ
                   
                   DO NOMEGA=1,CHI%NOMEGA
                      IF (AIMAG(CHI%COMEGA(NOMEGA))/=0) THEN
                         ! F(NOMEGA)=F(NOMEGA)+LIND_IMAG(Q,CHI%OMEGA(NOMEGA))
                         ! complex version
                      ELSE
                         ! real version
                         ! F(NOMEGA)=F(NOMEGA)+LIND_REAL(Q,CHI%OMEGA(NOMEGA))
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          F=F/NINTER**3
       ELSE
          CALL vtutor%bug("internal error in XI_LIND: HFRCUT==0 not yet implemented", __FILE__, __LINE__)
       ENDIF

    ENDIF
  END SUBROUTINE XI_LIND


!********************** SUBROUTINE XI_INVERT **************************
!
! invert the dielectric matrix
! The matrix has two properties
! )  first it is Hermitian at q=0 since
!   RESPONSEFUN(G, G') = \sum_ai C_ai(G) C^*_ai(G')
!   RESPONSEFUN(G',G)* = \sum_ai C*_ai(G') C_ai(G') = RESPONSEFUN(G, G')
! ) second it posses the property 
!   RESPONSEFUN_q(G, G') = RESPONSEFUN*_-q(-G',-G)
!   this follows from the fact that the response function is real
!   in real space
! 
! after calling the routine the inverted response function
! is stored in CHI%RESPONSEFUN
! the inverse of the long range part is stored in 
!    CHI%HEAD
! for the q->0 case block inversion
! similar to the routine DETERMINE_FXC_FROM_TBSE_IDIR should be
! used
!
!**********************************************************************

  SUBROUTINE XI_INVERT( IU0, CHI, WGWQ)
    USE prec
    USE constant
    USE mgrid
    USE wave 
    USE full_kpoints

    IMPLICIT NONE
    INTEGER IU0
    TYPE (responsefunction) CHI
    TYPE (wavedes1) WGWQ
    ! local
    INTEGER    NP, I, II, NOMEGA, NQ_IN_WGW, INFO, IDIR, JDIR
    INTEGER   IPIV(CHI%NP2)
    INTEGER, PARAMETER :: NWORK=64
    COMPLEX(q) :: WORK(CHI%NP2*NWORK)
    COMPLEX(q) :: CHI_WORK(CHI%NP2, CHI%NP2)
    COMPLEX(q) :: EPS_WORK(CHI%NP2, CHI%NP2)

    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2

    DO NOMEGA=1,CHI%NOMEGA
       CALL GWPROGRESS(IU0, NOMEGA, CHI%NOMEGA, 1, 1)
       IF (CHI%LGAMMA) THEN

          ! diagonalize the matrix for three directions
          ! for a non orthorhombic material the present implementation
          ! does not work
          ! one solution
          !  ) diagonalisation of the 3x3 RPA head CHI%HEAD(:,:,NOMEGA)
          !  ) loop over all three eigenvectors a(i) i=1,3 of the head 
          !  )   set CHI_WORK(:,1) to \sum_j  a(i)_j CHI%WING(:,j,NOMEGA)
          !  )   set CHI_WORK(1,1) to \sum_jk a(i)_j CHI%HEAD(i,j,NOMEGA) a(i)_k
          !  )   diagonalize the full matrix 
          !  ) add to INV_EPSILON_MACRO \sum_jk a(i)_j CHI_WORK(1,1) a(i)_k
          ! or change to the block diagonalization
          ! of Baroni and Resta Phys. Rev. B 33, 7017 (1986).
          EPS_WORK=0
          DO IDIR=1,IDIR_MAX
             CALL SET_MAT_FROM_RESPONSE( CHI_WORK, CHI, NOMEGA, IDIR)

             INFO=0
             ! compose into uper and lower triangular matrix
             CALL ZGETRF( NP, NP, CHI_WORK(1,1), SIZE(CHI_WORK,1), IPIV, INFO )
             IF (INFO/=0) THEN
                CALL vtutor%error("error in XI_INVERT: ZGETRF returns " // str(INFO))
             ENDIF
             ! invert matrix
             CALL ZGETRI( NP, CHI_WORK(1,1), SIZE(CHI_WORK,1), IPIV, & 
                  WORK, SIZE(CHI%RESPONSEFUN,1)*NWORK, INFO )
             IF (INFO/=0) THEN
                CALL vtutor%error("error in XI_INVERT: ZGETRI returns " // str(INFO))
             ENDIF

             ! accumulate inverted matrices
             EPS_WORK=EPS_WORK+CHI_WORK
             ! set head and wing to calculated values
             CALL SET_WING_FROM_MAT( CHI_WORK, CHI, NOMEGA, IDIR)
          ENDDO
          ! clean non diagonal head elements
          CHI%HEAD(2,1, NOMEGA)=0 ; CHI%HEAD(3,1, NOMEGA)=0 ; CHI%HEAD(3,2, NOMEGA)=0 
          CHI%HEAD(1,2, NOMEGA)=0 ; CHI%HEAD(1,3, NOMEGA)=0 ; CHI%HEAD(3,2, NOMEGA)=0 

          EPS_WORK=EPS_WORK*(1.0_q/IDIR_MAX)
          ! body is simply the average
          CALL SET_RESPONSE_FROM_MAT(EPS_WORK, CHI, NOMEGA, 0)
       ELSE
          CHI_WORK=CHI%RESPONSEFUN(:,:,NOMEGA)

          INFO=0
          CALL ZGETRF( NP, NP, CHI_WORK(1,1), SIZE(CHI_WORK,1), IPIV, INFO )
          IF (INFO/=0) THEN
             CALL vtutor%error("error in XI_INVERT: ZGETRF returns " // str(INFO))
          ENDIF
          CALL ZGETRI( NP, CHI_WORK(1,1), SIZE(CHI_WORK,1), IPIV, & 
               WORK, SIZE(CHI%RESPONSEFUN,1)*NWORK, INFO )
          IF (INFO/=0) THEN
             CALL vtutor%error("error in XI_INVERT: ZGETRI returns " // str(INFO))
          ENDIF

          CALL SET_RESPONSE_FROM_MAT(CHI_WORK, CHI, NOMEGA)
       ENDIF
    ENDDO
    CALL GWPROGRESS(IU0, CHI%NOMEGA, CHI%NOMEGA, 1, 1)
    CALL RESTORE_HEAD( CHI )

  END SUBROUTINE XI_INVERT


!********************** SUBROUTINE XI_LOCAL_FIELD *********************
!                                   -1                 -1
! calculate  X_f= X_0 (1- f_xc X_0)   = ( 1-  X_0 f_xc)   X_0
! the matrix is Hermitian since it can be reformulated as
!                                    -1
!            X_0 (X_0 - X_0 f_xc X_0)   X_0
!
! if LHARTREE is .TRUE., the Hartree kernel is added i.e. the reducible
! polarizability is calculated as
!                                  -1
!  X_red= X_0 (1- f_xc X_0 - v X_0) 
!
! is calculated
!
! note that local field effects can be included in two at first sight
! dissimilar manners:
! Resta and co. use:
!    -1             -1                           -1  -1
! eps =  (1 - v X_f)   = (1 - v X_0 (1- f_xc X_0)   )
!                                           -1  -1
!     =  ((1-f_xc X_0 - v X_0) (1- f_xc X_0)   )
!                                          -1
!     =  (1-f_xc X_0) (1- f_xc X_0 - v X_0)
!                                                -1                             -1
!     =  (1-f_xc X_0-v X_0)(1- f_xc X_0 - v X_0)  + v X_0 (1- f_xc X_0 - v X_0)
!    -1                                  -1
! eps   = 1 + v X_0 (1- f_xc X_0 - v X_0)
! 
! this is the expression used by GW people
! it is sort of nicer since, since it invokes only a single
! inversion, one can easily neglect local field effects entirely
! or include them on any required level
! 
! TODO:
! local field effects should be done properly using a block inversion
! similar to the routine DETERMINE_FXC_FROM_TBSE_IDIR
!
!**********************************************************************

  SUBROUTINE XI_LOCAL_FIELD( IU0, CHI, TVXC, TBSE, WGWQ, LHARTREE)
    USE prec
    USE constant
    USE mgrid
    USE wave 
    USE full_kpoints

    IMPLICIT NONE
    TYPE (responsefunction) CHI
    TYPE (responsefunction) TVXC, TBSE
    TYPE (wavedes1) :: WGWQ
    LOGICAL    :: LHARTREE
    INTEGER    :: IU0
    ! local
    INTEGER    NP, I, II, NOMEGA, NQ_IN_WGW, INFO, IDIR, JDIR
    INTEGER   IPIV(CHI%NP2)
    INTEGER, PARAMETER :: NWORK=64
    COMPLEX(q) :: WORK(CHI%NP2*NWORK)
    COMPLEX(q) :: CHI_WORK(CHI%NP2, CHI%NP2)
    COMPLEX(q) :: EPS_WORK(CHI%NP2, CHI%NP2)
    COMPLEX(q) :: INV_EPSILON_MACRO(3,3)
      
    IF (.NOT. ASSOCIATED(TVXC%RESPONSEFUN).AND. .NOT. ASSOCIATED(TBSE%RESPONSEFUN) & 
         .AND. .NOT. LHARTREE) RETURN

    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2

    DO NOMEGA=1,CHI%NOMEGA
       CALL GWPROGRESS(IU0, NOMEGA, CHI%NOMEGA, 1, 1)
       IF (CHI%LGAMMA) THEN
          ! diagonalize the matrix for three directions
          ! this could be reworked using the block diagonalization
          ! of Baroni and Resta Phys. Rev. B 33, 7017 (1986).

          ! the stored head and wing in chi describes the 
          ! q^2 and q behaviour of Xi
          !            | head q^2  wing q|
          ! Xi_0(q)  = |                 |
          !            | wing q    body  |
          ! for TVXC the head and wing stores the 1/q^2 and 1/q behaviour
          EPS_WORK=0
          DO IDIR=1,IDIR_MAX
             CALL BODY_FROM_WING( CHI, IDIR, NOMEGA)

             IF (ASSOCIATED(TVXC%RESPONSEFUN)) CALL BODY_FROM_WING( TVXC, IDIR)
             IF (ASSOCIATED(TBSE%RESPONSEFUN)) CALL BODY_FROM_WING( TBSE, IDIR)

             ! calculate  Xi_0 f_xc and store in CHI_WORK
             IF (ASSOCIATED(TVXC%RESPONSEFUN).AND. ASSOCIATED(TBSE%RESPONSEFUN)) THEN
                CALL MATMUL_RESPONSE(CHI_WORK,CHI%RESPONSEFUN(:,:,NOMEGA),TVXC%RESPONSEFUN(:,:,1)+TBSE%RESPONSEFUN(:,:,1)*SCALE_TBSE, WGWQ)
             ELSEIF (ASSOCIATED(TVXC%RESPONSEFUN)) THEN
                CALL MATMUL_RESPONSE(CHI_WORK,CHI%RESPONSEFUN(:,:,NOMEGA),TVXC%RESPONSEFUN(:,:,1), WGWQ)
             ELSEIF (ASSOCIATED(TBSE%RESPONSEFUN)) THEN
                CALL MATMUL_RESPONSE(CHI_WORK,CHI%RESPONSEFUN(:,:,NOMEGA),TBSE%RESPONSEFUN(:,:,1)*SCALE_TBSE, WGWQ)
             ELSE
                CHI_WORK=0
             ENDIF

             IF (LHARTREE) THEN
                 CALL XI_LOCAL_FIELD_HARTREE_T( CHI_WORK, CHI, NOMEGA, WGWQ )
             ENDIF

             CHI_WORK(1:NP,1:NP)=-CHI_WORK(1:NP,1:NP)

             DO I=1,NP
                CHI_WORK(I,I)=1+CHI_WORK(I,I)
             ENDDO

             INFO=0
             ! compose into upper and lower triangular matrix
             CALL ZGETRF( NP, NP, CHI_WORK(1,1), SIZE(CHI_WORK,1), IPIV, INFO )
             IF (INFO/=0) THEN
                CALL vtutor%error("error in XI_LOCAL_FIELD: ZGETRF returns " // str(INFO))
             ENDIF
             ! invert matrix
             CALL ZGETRI( NP, CHI_WORK(1,1), SIZE(CHI_WORK,1), IPIV, & 
                  WORK, SIZE(CHI%RESPONSEFUN,1)*NWORK, INFO )
             IF (INFO/=0) THEN
                CALL vtutor%error("error in XI_LOCAL_FIELD: ZGETRI returns " // str(INFO))
             ENDIF

             ! CHI_WORK(1:NP,1:NP)=MATMUL(CHI_WORK(1:NP,1:NP),CHI%RESPONSEFUN(1:NP,1:NP,NOMEGA))
             CALL MATMUL_RIGHT(CHI_WORK, CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)
             EPS_WORK=EPS_WORK+CHI_WORK
             ! set head and wing to calculated values
             CALL SET_WING_FROM_MAT( CHI_WORK, CHI, NOMEGA, IDIR)
          ENDDO
          EPS_WORK=EPS_WORK*(1.0_q/IDIR_MAX)
          ! body is simply the average
          CALL SET_RESPONSE_FROM_MAT(EPS_WORK, CHI, NOMEGA, 0)
       ELSE
          ! calculate  Xi_0 f_xc and store in CHI_WORK
          IF (ASSOCIATED(TVXC%RESPONSEFUN).AND. ASSOCIATED(TBSE%RESPONSEFUN)) THEN
             CALL MATMUL_RESPONSE(CHI_WORK,CHI%RESPONSEFUN(:,:,NOMEGA),TVXC%RESPONSEFUN(:,:,1)+TBSE%RESPONSEFUN(:,:,1)*SCALE_TBSE, WGWQ)
          ELSEIF (ASSOCIATED(TVXC%RESPONSEFUN)) THEN
             CALL MATMUL_RESPONSE(CHI_WORK,CHI%RESPONSEFUN(:,:,NOMEGA),TVXC%RESPONSEFUN(:,:,1), WGWQ)
          ELSEIF (ASSOCIATED(TBSE%RESPONSEFUN)) THEN
             CALL MATMUL_RESPONSE(CHI_WORK,CHI%RESPONSEFUN(:,:,NOMEGA),TBSE%RESPONSEFUN(:,:,1)*SCALE_TBSE, WGWQ)
          ELSE
             CHI_WORK=0
          ENDIF

          IF (LHARTREE) CALL XI_LOCAL_FIELD_HARTREE_T( CHI_WORK, CHI, NOMEGA, WGWQ )

          CHI_WORK(1:NP,1:NP)=-CHI_WORK(1:NP,1:NP)

          DO I=1,NP
             CHI_WORK(I,I)=1+CHI_WORK(I,I)
          ENDDO

          INFO=0
          ! compose into upper and lower triangular matrix
          CALL ZGETRF( NP, NP, CHI_WORK(1,1), SIZE(CHI_WORK,1), IPIV, INFO )
          IF (INFO/=0) THEN
             CALL vtutor%error("error in XI_LOCAL_FIELD: ZGETRF returns " // str(INFO))
          ENDIF
          ! invert matrix
          CALL ZGETRI( NP, CHI_WORK(1,1), SIZE(CHI_WORK,1), IPIV, & 
               WORK, SIZE(CHI%RESPONSEFUN,1)*NWORK, INFO )
          IF (INFO/=0) THEN
             CALL vtutor%error("error in XI_LOCAL_FIELD: ZGETRI returns " // str(INFO))
          ENDIF
          
          ! CHI_WORK(1:NP,1:NP)=MATMUL(CHI_WORK(1:NP,1:NP),CHI%RESPONSEFUN(1:NP,1:NP,NOMEGA))
          CALL MATMUL_RIGHT(CHI_WORK, CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)
          CALL SET_RESPONSE_FROM_MAT(CHI_WORK, CHI, NOMEGA)
       ENDIF
    ENDDO
    ! TODO UINIT: On multiple cores only rank 0 has non-zero CHI%OMEGA. The
    ! subroutine will be executed with arguments (IU0, 0, 0, 1, 1) producing
    ! a division by zero in chi_base.F:4124 (Co_CRPA_proj). This is normally
    ! without consequences, what is the best way to avoid it?
DOESI IF (CHI%NOMEGA /= 0) &
    CALL GWPROGRESS(IU0, CHI%NOMEGA, CHI%NOMEGA, 1, 1)

    CALL RESTORE_HEAD( CHI )

  END SUBROUTINE XI_LOCAL_FIELD



!********************** SUBROUTINE XI_FXC_FROM_EPS  *******************
!
! this is the simple exchange kernel suggested in
!  S. Sharma, J.K. Dewhurst, A. Sanna, and E.K.U Gross, Phys. Rev. Lett
!
! f_xc = epsilon-1 * Xi-1
!
!**********************************************************************

  SUBROUTINE XI_FXC_FROM_EPS( IU0, EPSINV, CHI0, TVXC, WGWQ)
    USE prec
    USE constant
    USE mgrid
    USE wave 
    USE full_kpoints

    IMPLICIT NONE
    TYPE (responsefunction) EPSINV
    TYPE (responsefunction) CHI0
    TYPE (responsefunction) TVXC
    TYPE (wavedes1) :: WGWQ
    INTEGER    :: IU0
    ! local
    INTEGER    NP, I, II, NOMEGA, NQ_IN_WGW, INFO, IDIR, JDIR
    INTEGER   IPIV(EPSINV%NP2)
    INTEGER, PARAMETER :: NWORK=64
    COMPLEX(q) :: WORK(EPSINV%NP2*NWORK)
    COMPLEX(q) :: CHI_WORK(EPSINV%NP2, EPSINV%NP2)
    COMPLEX(q) :: CHI_WORK2(EPSINV%NP2, EPSINV%NP2)
    COMPLEX(q) :: EPS_WORK(EPSINV%NP2, EPSINV%NP2)
    COMPLEX(q) :: INV_EPSILON_MACRO(3,3)
!    INTEGER :: IDIR_MAX=1
      
    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2

    DO NOMEGA=1,1
       TVXC%HEAD=0
       IF (EPSINV%LGAMMA) THEN
          EPS_WORK=0
          DO IDIR=1,IDIR_MAX
!             WRITE(*,'(6F14.7)') EPSINV%HEAD(:,:,1)
!             WRITE(*,'(6F14.7)') CHI0%HEAD
             CALL BODY_FROM_WING( EPSINV,  IDIR, NOMEGA)
             CALL BODY_FROM_WING( CHI0, IDIR, NOMEGA)

             CHI_WORK(1:NP, 1:NP)=CHI0%RESPONSEFUN(1:NP, 1:NP, 1)

             INFO=0
             ! compose into upper and lower triangular matrix
             !CALL ZGETRF( NP, NP, CHI_WORK(1,1), SIZE(CHI_WORK,1), IPIV, INFO )
             !IF (INFO/=0) THEN
             !   WRITE(0,*) 'error in XI_FXC_FROM_EPS: ZGETRF returns',INFO
             !   STOP
             !ENDIF
             ! invert matrix
             !CALL ZGETRI( NP, CHI_WORK(1,1), SIZE(CHI_WORK,1), IPIV, & 
             !     WORK, SIZE(CHI0%RESPONSEFUN,1)*NWORK, INFO )
             !IF (INFO/=0) THEN
             !   WRITE(0,*) 'error in XI_FXC_FROM_EPS: ZGETRI returns',INFO
             !   STOP
             !ENDIF


             ! inversion using diagonalization and possibly removing singular contributions
          !   CALL ROTINV( CHI_WORK, NP, 6)
          !   CHI_WORK2=CHI_WORK

             ! simple model from Hardy Gross slide (confirmed private comm. Sharma)
             CHI_WORK=0
             DO I=1,NP
                CHI_WORK(I,I)=1/CHI0%HEAD(IDIR,IDIR,1)
             ENDDO

             ! f_xc = Xi_0^-1 eps^-1
             ! CHI_WORK(1:NP,1:NP)=MATMUL(CHI_WORK(1:NP,1:NP),EPSINV%RESPONSEFUN(1:NP,1:NP,NOMEGA))
             CALL MATMUL_RIGHT(CHI_WORK, EPSINV%RESPONSEFUN(:,:,NOMEGA), WGWQ)
          !   CALL MATMUL_LEFT(CHI_WORK, CHI_WORK2, WGWQ)
          !   CHI_WORK=-CHI_WORK

             EPS_WORK=EPS_WORK+CHI_WORK
             ! set head and wing to calculated values
             CALL SET_WING_FROM_MAT( CHI_WORK,  TVXC, NOMEGA, IDIR)
          ENDDO

          EPS_WORK=EPS_WORK*(1.0_q/IDIR_MAX)
!          WRITE(*,'(8F12.5)') EPS_WORK(1:4,1:4)
          ! body is simply the average
          CALL SET_RESPONSE_FROM_MAT(EPS_WORK, TVXC, NOMEGA, 0)
       ELSE
             CHI_WORK(1:NP, 1:NP)=-CHI0%RESPONSEFUN(1:NP, 1:NP, 1)

             ! inversion using diagonalization and possibly removing singular contributions
             CALL ROTINV( CHI_WORK, NP, -1)

             INFO=0
             ! compose into upper and lower triangular matrix
             !CALL ZGETRF( NP, NP, CHI_WORK(1,1), SIZE(CHI_WORK,1), IPIV, INFO )
             !IF (INFO/=0) THEN
             !   WRITE(0,*) 'error in XI_FXC_FROM_EPS: ZGETRF returns',INFO
             !   STOP
             !ENDIF
             ! invert matrix
             !CALL ZGETRI( NP, CHI_WORK(1,1), SIZE(CHI_WORK,1), IPIV, & 
             !     WORK, SIZE(CHI%RESPONSEFUN,1)*NWORK, INFO )
             !IF (INFO/=0) THEN
             !   WRITE(0,*) 'error in XI_FXC_FROM_EPS: ZGETRI returns',INFO
             !   STOP
             !ENDIF
             
             ! f_xc = Xi_0^-1 eps^-1
             ! CHI_WORK(1:NP,1:NP)=MATMUL(CHI_WORK(1:NP,1:NP),EPSINV%RESPONSEFUN(1:NP,1:NP,NOMEGA))
             CALL MATMUL_RIGHT(CHI_WORK, EPSINV%RESPONSEFUN(:,:,NOMEGA), WGWQ)
             CALL SET_RESPONSE_FROM_MAT(CHI_WORK, TVXC, NOMEGA)
       ENDIF
    ENDDO

    CALL RESTORE_HEAD( EPSINV )

  END SUBROUTINE XI_FXC_FROM_EPS


!********************** SUBROUTINE XI_LOCAL_FIELD_SYM *****************
!
! symmetric version of the previous routine
!
! calculate  the irreducible  polarizabilty from the independent 
! particle polarizability using
!                              -1
!  X_f= X_0 (X_0- X_0 f_xc X_0)   X_0
! 
! although this is algebraically equivalent to
!                           -1
!  X_f= X_0 (1 -  f_xc X_0) 
!
! the first version maintains the symmetry of f_xc, whereas the latter
! does not.
! furthermore f_xc is always left and right multiplied by X_0, making
! it numerically stable.
! finally it allows to handle elegantly resonant (anti) resonant coupling
!
!  X_f= X_0 (X_0 - X_r f^r  X_r - X_ar f^r X_ar -
!                 X_ar f^ar X_r -  X_r f^ar X_ar)^-1   X_0
!
! if LH is passed the Hartree kernel is added i.e. the reducible
! polarizability is calculated as
!                                            -1
!  X_red= X_0 (X_0- X_0 f_xc X_0 - X_0 v X_0)   X_0
!
! is calculated
!
!
!**********************************************************************

  SUBROUTINE XI_LOCAL_FIELD_SYM( IU0, CHI, CHIR, TVXC, TBSE, TBSEA, WGWQ, LHARTREE)
    USE prec
    USE constant
    USE mgrid
    USE wave 
    USE full_kpoints

    IMPLICIT NONE
    INTEGER IU0
    TYPE (responsefunction) CHI
    TYPE (responsefunction) CHIR
    TYPE (responsefunction) TVXC, TBSE, TBSEA
    TYPE (wavedes1) :: WGWQ
    LOGICAL    :: LHARTREE
    ! local
    INTEGER    NP, I, II, NOMEGA, NQ_IN_WGW, INFO, IDIR, JDIR
    INTEGER   IPIV(CHI%NP2)
    INTEGER, PARAMETER :: NWORK=64
    COMPLEX(q) :: WORK(CHI%NP2*NWORK)
    COMPLEX(q) :: CHI_WORK(CHI%NP2, CHI%NP2)
    COMPLEX(q) :: EPS_WORK(CHI%NP2, CHI%NP2)
    COMPLEX(q) :: CHI_RES(CHI%NP2, CHI%NP2)
    COMPLEX(q) :: INV_EPSILON_MACRO(3,3)
      
    IF (.NOT. ASSOCIATED(TVXC%RESPONSEFUN).AND. .NOT. ASSOCIATED(TBSE%RESPONSEFUN) & 
         .AND. .NOT. LHARTREE) RETURN

    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2

    DO NOMEGA=1,CHI%NOMEGA
       CALL GWPROGRESS(IU0, NOMEGA, CHI%NOMEGA, 1, 1)
       IF (CHI%LGAMMA) THEN

          ! diagonalize the matrix for three directions
          ! this could be reworked using the block diagonalization
          ! of Baroni and Resta Phys. Rev. B 33, 7017 (1986).

          ! the stored head and wing in chi describes the 
          ! q^2 and q behaviour of Xi
          !            | head q^2  wing q|
          ! Xi_0(q)  = |                 |
          !            | wing q    body  |
          ! for TVXC the head and wing stores the 1/q^2 and 1/q behaviour
          EPS_WORK=0
          DO IDIR=1,IDIR_MAX
             CALL BODY_FROM_WING( CHI, IDIR, NOMEGA)
             IF (ASSOCIATED(CHIR%RESPONSEFUN)) CALL BODY_FROM_WING( CHIR, IDIR, NOMEGA)
             IF (ASSOCIATED(TVXC%RESPONSEFUN)) CALL BODY_FROM_WING( TVXC, IDIR)
             IF (ASSOCIATED(TBSE%RESPONSEFUN)) CALL BODY_FROM_WING( TBSE, IDIR)
             IF (ASSOCIATED(TBSEA%RESPONSEFUN))CALL BODY_FROM_WING( TBSEA, IDIR)
             
             
             IF (ASSOCIATED(TBSE%RESPONSEFUN)) THEN
             ! if a direction dependent fxc is used (which is required if the head
             !  is taken into account at Gamma) the x, y and z are stored in 2,3,4
             ! if a speficic component is need it is copied back to 1
                IF (SIZE(TBSE%RESPONSEFUN,3)==4) THEN
                   TBSE%RESPONSEFUN(:,:,1) =TBSE%RESPONSEFUN (:,:,IDIR+1)
                ENDIF
             ENDIF

             IF (ASSOCIATED(TBSEA%RESPONSEFUN)) THEN
                IF (SIZE(TBSEA%RESPONSEFUN,3)==4) THEN
                   TBSEA%RESPONSEFUN(:,:,1)=TBSEA%RESPONSEFUN(:,:,IDIR+1)
                ENDIF
             ENDIF

             ! calculate  Xi_0 f_xc and store in CHI_WORK
             CHI_WORK=0
             IF (ASSOCIATED(TVXC%RESPONSEFUN).AND. & 
                (ASSOCIATED(TBSE%RESPONSEFUN).AND. .NOT.ASSOCIATED(CHIR%RESPONSEFUN))) THEN
                CALL MATMUL_RESPONSE(CHI_WORK,CHI%RESPONSEFUN(:,:,NOMEGA),TVXC%RESPONSEFUN(:,:,1)+TBSE%RESPONSEFUN(:,:,1)*SCALE_TBSE, WGWQ)
             ELSEIF (ASSOCIATED(TVXC%RESPONSEFUN)) THEN
                CALL MATMUL_RESPONSE(CHI_WORK,CHI%RESPONSEFUN(:,:,NOMEGA),TVXC%RESPONSEFUN(:,:,1), WGWQ)
             ELSEIF (ASSOCIATED(TBSE%RESPONSEFUN).AND. .NOT.ASSOCIATED(CHIR%RESPONSEFUN)) THEN
                CALL MATMUL_RESPONSE(CHI_WORK,CHI%RESPONSEFUN(:,:,NOMEGA),TBSE%RESPONSEFUN(:,:,1)*SCALE_TBSE, WGWQ)
             ENDIF

             IF (LHARTREE) CALL XI_LOCAL_FIELD_HARTREE_T( CHI_WORK, CHI, NOMEGA, WGWQ )
             ! apply Hartree potential but remove q=0 component 
             !IF (LHARTREE) CALL XI_HARTREEBAR_T( CHI_WORK, CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ )

             CHI_WORK(1:NP,1:NP)=-CHI_WORK(1:NP,1:NP)
             DO I=1,NP
                CHI_WORK(I,I)=1+CHI_WORK(I,I)
             ENDDO
             CALL MATMUL_RIGHT(CHI_WORK, CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)

             IF (ASSOCIATED(TBSE%RESPONSEFUN) .AND. ASSOCIATED(CHIR%RESPONSEFUN)) THEN
                ! resonant part X_r f_r X_r
                CALL MATMUL_RESPONSE(CHI_RES, CHIR%RESPONSEFUN(:,:,NOMEGA),TBSE%RESPONSEFUN(:,:,1), WGWQ)
                CALL MATMUL_RIGHT(CHI_RES, CHIR%RESPONSEFUN(:,:,NOMEGA), WGWQ)

                CHI_WORK(1:NP,1:NP)=CHI_WORK(1:NP,1:NP)-CHI_RES(1:NP,1:NP)

                ! anti-resonant part X_ar f_r X_ar
                CALL MATMUL_RESPONSE(CHI_RES,CHI%RESPONSEFUN(:,:,NOMEGA)-CHIR%RESPONSEFUN(:,:,NOMEGA),TBSE%RESPONSEFUN(:,:,1), WGWQ)
                CALL MATMUL_RIGHT(CHI_RES, CHI%RESPONSEFUN(:,:,NOMEGA)-CHIR%RESPONSEFUN(:,:,NOMEGA), WGWQ)
                CHI_WORK(1:NP,1:NP)=CHI_WORK(1:NP,1:NP)-CHI_RES(1:NP,1:NP)

                ! resonant-antiresonant part  X_r f_ar X_ar
                IF (ASSOCIATED(TBSEA%RESPONSEFUN)) THEN
                   CALL MATMUL_RESPONSE(CHI_RES,CHIR%RESPONSEFUN(:,:,NOMEGA),TBSEA%RESPONSEFUN(:,:,1), WGWQ)
                   CALL MATMUL_RIGHT(CHI_RES, CHI%RESPONSEFUN(:,:,NOMEGA)-CHIR%RESPONSEFUN(:,:,NOMEGA), WGWQ)
                   CHI_WORK(1:NP,1:NP)=CHI_WORK(1:NP,1:NP)-CHI_RES(1:NP,1:NP)

                ! antiresonant-resonant part X_ar f_ar X_r
                   CALL MATMUL_RESPONSE(CHI_RES,CHI%RESPONSEFUN(:,:,NOMEGA)-CHIR%RESPONSEFUN(:,:,NOMEGA),TBSEA%RESPONSEFUN(:,:,1), WGWQ)
                   CALL MATMUL_RIGHT(CHI_RES, CHIR%RESPONSEFUN(:,:,NOMEGA), WGWQ)
                   CHI_WORK(1:NP,1:NP)=CHI_WORK(1:NP,1:NP)-CHI_RES(1:NP,1:NP)
                ENDIF
             ENDIF

             INFO=0

             ! the second column and row are exactly zero, get rid of it
             IF (WGWQ%LGAMMA) CHI_WORK(2,2)=100

             ! now invert X_0 + X_0 (f_xc +v) X_0 this is possibly pretty unstable, unfortunately
             ! L U decomposition
             CALL ZGETRF( NP, NP, CHI_WORK(1,1), SIZE(CHI_WORK,1), IPIV, INFO )
             IF (INFO/=0) THEN
                CALL vtutor%error("error in XI_LOCAL_FIELD_SYM: ZGETRF returns " // str(INFO))
             ENDIF
             ! invert matrix
             CALL ZGETRI( NP, CHI_WORK(1,1), SIZE(CHI_WORK,1), IPIV, & 
                  WORK, SIZE(CHI%RESPONSEFUN,1)*NWORK, INFO )
             ! zero second column and row (now 1/100)
             IF (WGWQ%LGAMMA) CHI_WORK(2,2)=0

             IF (INFO/=0) THEN
                CALL vtutor%error("error in XI_LOCAL_FIELD_SYM: ZGETRI returns " // str(INFO))
             ENDIF

             ! left and right multiply by X_0
             ! CHI_WORK(1:NP,1:NP)=MATMUL(CHI_WORK(1:NP,1:NP),CHI%RESPONSEFUN(1:NP,1:NP,NOMEGA))
             CALL MATMUL_RIGHT(CHI_WORK, CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)
             CALL MATMUL_LEFT(CHI_WORK, CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)

             EPS_WORK=EPS_WORK+CHI_WORK

             CALL SET_WING_FROM_MAT( CHI_WORK, CHI, NOMEGA, IDIR)
          ENDDO
          EPS_WORK=EPS_WORK*(1.0_q/IDIR_MAX)
          ! body is simply the average
          CALL SET_RESPONSE_FROM_MAT(EPS_WORK, CHI, NOMEGA, 0)
       ELSE
             ! calculate  Xi_0 f_xc and store in CHI_WORK
             IF (ASSOCIATED(TVXC%RESPONSEFUN).AND. & 
                (ASSOCIATED(TBSE%RESPONSEFUN).AND. .NOT.ASSOCIATED(CHIR%RESPONSEFUN))) THEN
                CALL MATMUL_RESPONSE(CHI_WORK,CHI%RESPONSEFUN(:,:,NOMEGA),TVXC%RESPONSEFUN(:,:,1)+TBSE%RESPONSEFUN(:,:,1)*SCALE_TBSE, WGWQ)
             ELSEIF (ASSOCIATED(TVXC%RESPONSEFUN)) THEN
                CALL MATMUL_RESPONSE(CHI_WORK,CHI%RESPONSEFUN(:,:,NOMEGA),TVXC%RESPONSEFUN(:,:,1), WGWQ)
             ELSEIF (ASSOCIATED(TBSE%RESPONSEFUN).AND. .NOT.ASSOCIATED(CHIR%RESPONSEFUN)) THEN
                CALL MATMUL_RESPONSE(CHI_WORK,CHI%RESPONSEFUN(:,:,NOMEGA),TBSE%RESPONSEFUN(:,:,1)*SCALE_TBSE, WGWQ)
             ELSE
                CHI_WORK=0
             ENDIF

             IF (LHARTREE) CALL XI_LOCAL_FIELD_HARTREE_T( CHI_WORK, CHI, NOMEGA, WGWQ )

             CHI_WORK(1:NP,1:NP)=-CHI_WORK(1:NP,1:NP)

             DO I=1,NP
                CHI_WORK(I,I)=1+CHI_WORK(I,I)
             ENDDO
             CALL MATMUL_RIGHT(CHI_WORK, CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)

             IF (ASSOCIATED(TBSE%RESPONSEFUN) .AND. ASSOCIATED(CHIR%RESPONSEFUN)) THEN
                ! resonant part
                CALL MATMUL_RESPONSE(CHI_RES, CHIR%RESPONSEFUN(:,:,NOMEGA),TBSE%RESPONSEFUN(:,:,1), WGWQ)
                CALL MATMUL_RIGHT(CHI_RES, CHIR%RESPONSEFUN(:,:,NOMEGA), WGWQ)
                CHI_WORK(1:NP,1:NP)=CHI_WORK(1:NP,1:NP)-CHI_RES(1:NP,1:NP)

                ! anti-resonant part
                CALL MATMUL_RESPONSE(CHI_RES,CHI%RESPONSEFUN(:,:,NOMEGA)-CHIR%RESPONSEFUN(:,:,NOMEGA),TBSE%RESPONSEFUN(:,:,1), WGWQ)
                CALL MATMUL_RIGHT(CHI_RES, CHI%RESPONSEFUN(:,:,NOMEGA)-CHIR%RESPONSEFUN(:,:,NOMEGA), WGWQ)
                CHI_WORK(1:NP,1:NP)=CHI_WORK(1:NP,1:NP)-CHI_RES(1:NP,1:NP)

                ! resonant-antiresonant part
                IF (ASSOCIATED(TBSEA%RESPONSEFUN)) THEN
                   CALL MATMUL_RESPONSE(CHI_RES,CHIR%RESPONSEFUN(:,:,NOMEGA),TBSEA%RESPONSEFUN(:,:,1), WGWQ)
                   CALL MATMUL_RIGHT(CHI_RES, CHI%RESPONSEFUN(:,:,NOMEGA)-CHIR%RESPONSEFUN(:,:,NOMEGA), WGWQ)
                   CHI_WORK(1:NP,1:NP)=CHI_WORK(1:NP,1:NP)-CHI_RES(1:NP,1:NP)

                ! antiresonant-resonant part
                   CALL MATMUL_RESPONSE(CHI_RES,CHI%RESPONSEFUN(:,:,NOMEGA)-CHIR%RESPONSEFUN(:,:,NOMEGA),TBSEA%RESPONSEFUN(:,:,1), WGWQ)
                   CALL MATMUL_RIGHT(CHI_RES, CHIR%RESPONSEFUN(:,:,NOMEGA), WGWQ)
                   CHI_WORK(1:NP,1:NP)=CHI_WORK(1:NP,1:NP)-CHI_RES(1:NP,1:NP)
                ENDIF
             ENDIF

             INFO=0
             ! compose into upper and lower triangular matrix
             CALL ZGETRF( NP, NP, CHI_WORK(1,1), SIZE(CHI_WORK,1), IPIV, INFO )
             IF (INFO/=0) THEN
                CALL vtutor%error("error in XI_LOCAL_FIELD_SYM: ZGETRF returns " // str(INFO))
             ENDIF
             ! invert matrix
             CALL ZGETRI( NP, CHI_WORK(1,1), SIZE(CHI_WORK,1), IPIV, & 
                  WORK, SIZE(CHI%RESPONSEFUN,1)*NWORK, INFO )
             IF (INFO/=0) THEN
                CALL vtutor%error("error in XI_LOCAL_FIELD_SYM: ZGETRI returns " // str(INFO))
             ENDIF
             
             ! CHI_WORK(1:NP,1:NP)=MATMUL(CHI_WORK(1:NP,1:NP),CHI%RESPONSEFUN(1:NP,1:NP,NOMEGA))
             CALL MATMUL_RIGHT(CHI_WORK, CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)
             CALL MATMUL_LEFT(CHI_WORK, CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)
             CALL SET_RESPONSE_FROM_MAT(CHI_WORK, CHI, NOMEGA)
       ENDIF

    ENDDO
    CALL GWPROGRESS(IU0, CHI%NOMEGA, CHI%NOMEGA, 1, 1)

    CALL RESTORE_HEAD( CHI )

  END SUBROUTINE XI_LOCAL_FIELD_SYM


!********************** SUBROUTINE XI_RED_TO_EPS **********************
!            
! calculate the epsilon from the reducible polarizability 
! if LTCTE = LOCAL_FIELD is set (test-charge test-electron) eps-1 
! is calculated as
!     -1
!  eps   = (1 + (f_xc +v) X_red )
!
! where X_red is defined as
!                                  -1
!  X_red= X_0 (1- f_xc X_0 - v X_0) 
!
! if LOCAL_FIELD=.FALSE.
!
!     -1
!  eps   = (1 + v X_red )
!
! is evaluated
!
!**********************************************************************

  SUBROUTINE XI_RED_TO_EPS( CHI, TVXC, TBSE, WGWQ , LOCAL_FIELD)
    USE prec
    USE mgrid
    USE wave 
    USE full_kpoints

    IMPLICIT NONE
    TYPE (responsefunction) CHI
    TYPE (responsefunction) TVXC, TBSE
    TYPE (wavedes1) WGWQ
    LOGICAL :: LOCAL_FIELD
    ! local
    INTEGER    NP, I, II, NOMEGA, NQ_IN_WGW, INFO, IDIR, JDIR
    INTEGER   IPIV(CHI%NP2 )
    INTEGER, PARAMETER :: NWORK=64
    COMPLEX(q) :: WORK(CHI%NP2*NWORK)
    COMPLEX(q) :: CHI_WORK(CHI%NP2, CHI%NP2)
    COMPLEX(q) :: EPS_WORK(CHI%NP2, CHI%NP2)
      
    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2

    DO NOMEGA=1,CHI%NOMEGA
       IF (CHI%LGAMMA) THEN
          EPS_WORK=0
          DO IDIR=1,IDIR_MAX
             CALL BODY_FROM_WING( CHI, IDIR, NOMEGA)

             IF (ASSOCIATED(TVXC%RESPONSEFUN)) CALL BODY_FROM_WING( TVXC, IDIR)
             IF (ASSOCIATED(TBSE%RESPONSEFUN)) CALL BODY_FROM_WING( TBSE, IDIR)

             ! calculate  Xi_0 (f_xc+v)
             IF (ASSOCIATED(TVXC%RESPONSEFUN).AND. ASSOCIATED(TBSE%RESPONSEFUN).AND. LOCAL_FIELD) THEN
                CALL MATMUL_RESPONSE(CHI_WORK,TVXC%RESPONSEFUN(:,:,1)+TBSE%RESPONSEFUN(:,:,1)*SCALE_TBSE,CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)
             ELSEIF (ASSOCIATED(TVXC%RESPONSEFUN) .AND. LOCAL_FIELD) THEN
                CALL MATMUL_RESPONSE(CHI_WORK,TVXC%RESPONSEFUN(:,:,1),CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)
             ELSEIF (ASSOCIATED(TBSE%RESPONSEFUN) .AND. LOCAL_FIELD) THEN
                CALL MATMUL_RESPONSE(CHI_WORK,TBSE%RESPONSEFUN(:,:,1)*SCALE_TBSE,CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)
             ELSE
                CHI_WORK=0
             ENDIF

             CALL XI_LOCAL_FIELD_HARTREE( CHI_WORK, CHI, NOMEGA, WGWQ )

             DO I=1,NP
                CHI_WORK(I,I)=1+CHI_WORK(I,I)
             ENDDO
             ! second entry is the sin transform corresponding to G=0, always zero
             IF (CHI%LREAL) CHI_WORK(2,2)=0

             ! accumulate
             EPS_WORK=EPS_WORK+CHI_WORK
             ! set head and wing to calculated values
             CALL SET_WING_FROM_MAT( CHI_WORK, CHI, NOMEGA, IDIR)
          ENDDO
          EPS_WORK=EPS_WORK*(1.0_q/IDIR_MAX)
         CALL SET_RESPONSE_FROM_MAT( EPS_WORK, CHI, NOMEGA, 0)
       ELSE

             ! calculate  Xi_0 (f_xc+v)
             IF (ASSOCIATED(TVXC%RESPONSEFUN).AND. ASSOCIATED(TBSE%RESPONSEFUN).AND. LOCAL_FIELD) THEN
                CALL MATMUL_RESPONSE(CHI_WORK,TVXC%RESPONSEFUN(:,:,1)+TBSE%RESPONSEFUN(:,:,1)*SCALE_TBSE,CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)
             ELSEIF (ASSOCIATED(TVXC%RESPONSEFUN) .AND. LOCAL_FIELD) THEN
                CALL MATMUL_RESPONSE(CHI_WORK,TVXC%RESPONSEFUN(:,:,1),CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)
             ELSEIF (ASSOCIATED(TBSE%RESPONSEFUN) .AND. LOCAL_FIELD) THEN
                CALL MATMUL_RESPONSE(CHI_WORK,TBSE%RESPONSEFUN(:,:,1)*SCALE_TBSE,CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)
             ELSE
                CHI_WORK=0
             ENDIF
          
             CALL XI_LOCAL_FIELD_HARTREE( CHI_WORK, CHI, NOMEGA, WGWQ )
             
             DO I=1,NP
                CHI_WORK(I,I)=1+CHI_WORK(I,I)
             ENDDO
             
             CALL SET_RESPONSE_FROM_MAT(CHI_WORK, CHI, NOMEGA)
       ENDIF

    ENDDO
    CALL RESTORE_HEAD( CHI )

  END SUBROUTINE XI_RED_TO_EPS

!********************** SUBROUTINE XI_RED_TO_TETE **********************
!            
! calculate the test-electron test-electron W
!    
!  W  = v + (f_xc +v) X_red (f_xc +v)
!
! where X_red is defined as
!                                  -1
!  X_red= X_0 (1- f_xc X_0 - v X_0) 
!
!**********************************************************************

  SUBROUTINE XI_RED_TO_TETE( CHI, TVXC, TBSE, WGWQ , LOCAL_FIELD)
    USE prec
    USE mgrid
    USE wave 
    USE full_kpoints

    IMPLICIT NONE
    TYPE (responsefunction) CHI
    TYPE (responsefunction) TVXC, TBSE
    TYPE (wavedes1) WGWQ
    LOGICAL :: LOCAL_FIELD
    ! local
    INTEGER    NP, I, II, NOMEGA, NQ_IN_WGW, INFO, IDIR, JDIR
    INTEGER   IPIV(CHI%NP2)
    INTEGER, PARAMETER :: NWORK=64
    COMPLEX(q) :: WORK(CHI%NP2*NWORK)
    COMPLEX(q) :: CHI_WORK(CHI%NP2, CHI%NP2)
    COMPLEX(q) :: CHI_WORK2(CHI%NP2, CHI%NP2)
    REAL(q), PARAMETER :: F=0.5
      
    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2

    DO NOMEGA=1,CHI%NOMEGA
       IF (CHI%LGAMMA) THEN
          DO IDIR=1,IDIR_MAX
             CALL BODY_FROM_WING( CHI, IDIR, NOMEGA)

             IF (ASSOCIATED(TVXC%RESPONSEFUN)) CALL BODY_FROM_WING( TVXC, IDIR)
             IF (ASSOCIATED(TBSE%RESPONSEFUN)) CALL BODY_FROM_WING( TBSE, IDIR)

             ! calculate  f_xc Xi_0 
             IF (ASSOCIATED(TVXC%RESPONSEFUN).AND. ASSOCIATED(TBSE%RESPONSEFUN).AND. LOCAL_FIELD) THEN
                CALL MATMUL_RESPONSE(CHI_WORK,(TVXC%RESPONSEFUN(:,:,1)+TBSE%RESPONSEFUN(:,:,1)*SCALE_TBSE)*F,CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)
             ELSEIF (ASSOCIATED(TVXC%RESPONSEFUN) .AND. LOCAL_FIELD) THEN
                CALL MATMUL_RESPONSE(CHI_WORK,TVXC%RESPONSEFUN(:,:,1)*F,CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)
             ELSEIF (ASSOCIATED(TBSE%RESPONSEFUN) .AND. LOCAL_FIELD) THEN
                CALL MATMUL_RESPONSE(CHI_WORK,TBSE%RESPONSEFUN(:,:,1)*F*SCALE_TBSE,CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)
             ELSE
                CHI_WORK=0
             ENDIF

             ! add v Xi_0 
             CALL XI_LOCAL_FIELD_HARTREE( CHI_WORK, CHI, NOMEGA, WGWQ )

             ! calculate  (f_xc +v) Xi_0 f_xc
             IF (ASSOCIATED(TVXC%RESPONSEFUN).AND. ASSOCIATED(TBSE%RESPONSEFUN).AND. LOCAL_FIELD) THEN
                CALL MATMUL_RESPONSE(CHI_WORK2,CHI_WORK,(TVXC%RESPONSEFUN(:,:,1)+TBSE%RESPONSEFUN(:,:,1)*SCALE_TBSE)*F, WGWQ)
             ELSEIF (ASSOCIATED(TVXC%RESPONSEFUN) .AND. LOCAL_FIELD) THEN
                CALL MATMUL_RESPONSE(CHI_WORK2,CHI_WORK,TVXC%RESPONSEFUN(:,:,1)*F,WGWQ)
             ELSEIF (ASSOCIATED(TBSE%RESPONSEFUN) .AND. LOCAL_FIELD) THEN
                CALL MATMUL_RESPONSE(CHI_WORK2,CHI_WORK,TBSE%RESPONSEFUN(:,:,1)*F*SCALE_TBSE,WGWQ)
             ELSE
                CHI_WORK2=0
             ENDIF

             ! add (1+ (f_xc +v) Xi_0) v
             DO I=1,NP
                CHI_WORK(I,I)=1+CHI_WORK(I,I)
             ENDDO
             CALL XI_LOCAL_FIELD_HARTREE_T_MAT( CHI_WORK2, CHI_WORK, WGWQ )

             ! set head and wing to calculated values
             CALL SET_WING_FROM_MAT( CHI_WORK2, CHI, NOMEGA, IDIR)

          ENDDO
          CALL SET_RESPONSE_FROM_MAT( CHI_WORK2, CHI, NOMEGA, 0)

          ! slight trouble with the head: divide it by one over potential operator
          ! since the convergence corrections require epsilon instead of W
          CHI%HEAD(:,:,NOMEGA)=CHI%HEAD(:,:,NOMEGA)/WGWQ%DATAKE(1,1)

       ELSE
          ! calculate  f_xc Xi_0 
          IF (ASSOCIATED(TVXC%RESPONSEFUN).AND. ASSOCIATED(TBSE%RESPONSEFUN).AND. LOCAL_FIELD) THEN
             CALL MATMUL_RESPONSE(CHI_WORK,(TVXC%RESPONSEFUN(:,:,1)+TBSE%RESPONSEFUN(:,:,1)*SCALE_TBSE)*F,CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)
          ELSEIF (ASSOCIATED(TVXC%RESPONSEFUN) .AND. LOCAL_FIELD) THEN
             CALL MATMUL_RESPONSE(CHI_WORK,TVXC%RESPONSEFUN(:,:,1)*F,CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)
          ELSEIF (ASSOCIATED(TBSE%RESPONSEFUN) .AND. LOCAL_FIELD) THEN
             CALL MATMUL_RESPONSE(CHI_WORK,TBSE%RESPONSEFUN(:,:,1)*F*SCALE_TBSE,CHI%RESPONSEFUN(:,:,NOMEGA), WGWQ)
          ELSE
             CHI_WORK=0
          ENDIF

          ! add v Xi_0 
          CALL XI_LOCAL_FIELD_HARTREE( CHI_WORK, CHI, NOMEGA, WGWQ )

          ! calculate  (f_xc +v) Xi_0 f_xc
          IF (ASSOCIATED(TVXC%RESPONSEFUN).AND. ASSOCIATED(TBSE%RESPONSEFUN).AND. LOCAL_FIELD) THEN
             CALL MATMUL_RESPONSE(CHI_WORK2,CHI_WORK,(TVXC%RESPONSEFUN(:,:,1)+TBSE%RESPONSEFUN(:,:,1)*SCALE_TBSE)*F, WGWQ)
          ELSEIF (ASSOCIATED(TVXC%RESPONSEFUN) .AND. LOCAL_FIELD) THEN
             CALL MATMUL_RESPONSE(CHI_WORK2,CHI_WORK,TVXC%RESPONSEFUN(:,:,1)*F,WGWQ)
          ELSEIF (ASSOCIATED(TBSE%RESPONSEFUN) .AND. LOCAL_FIELD) THEN
             CALL MATMUL_RESPONSE(CHI_WORK2,CHI_WORK,TBSE%RESPONSEFUN(:,:,1)*F*SCALE_TBSE,WGWQ)
          ELSE
             CHI_WORK2=0
          ENDIF

          ! add (1+ (f_xc +v) Xi_0) v
          DO I=1,NP
             CHI_WORK(I,I)=1+CHI_WORK(I,I)
          ENDDO
          CALL XI_LOCAL_FIELD_HARTREE_T_MAT( CHI_WORK2, CHI_WORK, WGWQ )

          CALL SET_RESPONSE_FROM_MAT(CHI_WORK2, CHI, NOMEGA)
       ENDIF

    ENDDO
    CALL RESTORE_HEAD( CHI )

  END SUBROUTINE XI_RED_TO_TETE

!*********************************************************************
!
! calculate the ionic, Hartree and kinetic energy contribution to the
! eigenvalues and store it in a static array
!  CELTOT_HARTREE_KINETIC
! the Fock exchange on the plane wave grid is stored in 
!  CELTOT_X
! this routine requires that the xc type is set to Hartree-Fock
! i.e. PUSH_XC_TYPE_FOR_GW must be called otherwise the routine
! does not work properly
!
!*********************************************************************

  SUBROUTINE SET_EIGENVALUE_HARTREE_KINETIC( &
       HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR, &
       T_INFO,INFO,IO,KPOINTS,GRID,GRID_SOFT, &
       GRIDC,GRIDUS,C_TO_US,SOFT_TO_C,SYMM, &
       CHTOT,DENCOR,CVTOT,CSTRF, &
       CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM, &
       CHDEN,SV,LMDIM,IRDMAX)

    USE base
    USE lattice
    USE charge
    USE pseudo
    USE lattice
    USE nonl_high
    USE msymmetry
    USE mpimy
    USE mgrid
    USE mkpoints
    USE poscar
    USE wave
    USE pot
    USE pawm
    USE wave_high
    USE subrot
    USE fock
    USE hamil_struct_def
    IMPLICIT NONE
    TYPE (ham_handle)  HAMILTONIAN
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    TYPE (wavedes)     WDES
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W          ! wavefunction
    TYPE (latt)        LATT_CUR
    TYPE (info_struct) INFO
    TYPE (in_struct)   IO
    TYPE (kpoints_struct) KPOINTS
    TYPE (grid_3d)     GRID       ! grid for wavefunctions
    TYPE (grid_3d)     GRID_SOFT  ! grid for soft chargedensity
    TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
    TYPE (grid_3d)     GRIDUS     ! temporary grid in us.F
    TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
    TYPE (transit)     SOFT_TO_C  ! index table between GRID_SOFT and GRIDC
    TYPE (symmetry) ::   SYMM

    INTEGER LMDIM,IRDMAX,IRDMAA

    COMPLEX(q)  CHTOT(GRIDC%MPLWV,WDES%NCDIJ) ! charge-density in real / reciprocal space
    RGRID       DENCOR(GRIDC%RL%NP)           ! partial core
    COMPLEX(q)  CVTOT(GRIDC%MPLWV,WDES%NCDIJ) ! local potential
    COMPLEX(q)  CSTRF(GRIDC%MPLWV,T_INFO%NTYP)

!  augmentation related quantities
    OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
             CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
             CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
!  paw sphere charge density
    INTEGER N_MIX_PAW
    REAL(q)  RHOLM(N_MIX_PAW,WDES%NCDIJ)
!  charge-density and potential on soft grid
    COMPLEX(q)  CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
    RGRID       SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ)
!  Hamiltonian
    REAL(q)   XCSIF(3,3)
!  local
    TYPE (energy)      E
    TYPE (wavespin)    W_TMP

    IF (.NOT. LGW .AND. .NOT. LQPBSE) RETURN
    IF (IO%IU0>=0) WRITE(IO%IU0,*) 'calculate exact exchange contribution'
!=======================================================================
!  calculate the total local potential
!  ionic contribution + Hartree +
!   on-site terms (including core valence in Hartree Fock approximation)
!=======================================================================
 !  calculate potential
    CALL POTLOK(GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES,  &
         INFO,P,T_INFO,E,LATT_CUR, &
         CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF)

 ! add the one center augmentation related terms  
    CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
         LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)
 ! finally add one center terms
    CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
         WDES%NCDIJ, LMDIM, CDIJ(1,1,1,1),  RHOLM, CRHODE(1,1,1,1), &
         E,  LMETA=.FALSE., LASPH=INFO%LASPH, LCOREL= .FALSE.  )

    W_TMP=W
!=======================================================================
! CELTOT_HARTREE_KINETIC = 
!  <T + V_ion + V_H > + <T+V_H+V_ion>^1  + <V_x>^1
! <V_x>^1 includes core valence exchange
!
! the one-center exchange term is included in CELTOT_HARTREE_KINETIC
! (LDA and GGA contributions are switched off since AEXX=1.0)
!=======================================================================
    ALLOCATE(CELTOT_HARTREE_KINETIC(WDES%NB_TOT, WDES%NKPTS, WDES%ISPIN))
    CELEN_HARTREE_KINETIC => CELTOT_HARTREE_KINETIC(WDES%NB_LOW:WDES%NB_TOT:WDES%NB_PAR,:,:)

    W_TMP%CELTOT=>CELTOT_HARTREE_KINETIC
    W_TMP%CELEN =>CELEN_HARTREE_KINETIC

    ! no valence PW HF contribution
    LHFCALC=.FALSE.
    CALL EDDIAG(HAMILTONIAN, GRID, LATT_CUR, NONLR_S, NONL_S, W_TMP, WDES, SYMM, &
         LMDIM, CDIJ, CQIJ, 0, SV, T_INFO, P, IO%IU0, E%EXHF, NBANDS_MAX= NBANDSGW/WDES%NB_PAR)
    CALL MRG_CEL(WDES,W_TMP)

!    IF (IO%IU6>=0) THEN
!       WRITE(IO%IU6,'(//A)') '  <phi| T + V_ion + V_H + V_xc(core-valence)-V_xc(valence) | phi>'
!       CALL WRITE_EIGENVAL( WDES, W_TMP, IO%IU6)
!    ENDIF
!=======================================================================
!  determine CELTOT_X =   < V_x >
!  expectation value of exchange operator on plane wave grid
!=======================================================================
    LHFCALC=.TRUE.

    ALLOCATE(CELTOT_X(WDES%NB_TOT, WDES%NKPTS, WDES%ISPIN))
DOESI CELTOT_X=(0.0_q, 0.0_q)
    CELEN_X => CELTOT_X(WDES%NB_LOW:WDES%NB_TOT:WDES%NB_PAR,:,:)

    W_TMP%CELTOT=>CELTOT_X
    W_TMP%CELEN =>CELEN_X

    CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W_TMP,WDES, SYMM, &
         LMDIM,CDIJ,CQIJ, 0, SV, T_INFO, P, IO%IU0, E%EXHF, NBANDS_MAX= NBANDSGW/WDES%NB_PAR)
    CALL MRG_CEL(WDES,W_TMP)

    CELTOT_X=CELTOT_X-CELTOT_HARTREE_KINETIC

!    IF (IO%IU6>=0) THEN
!       WRITE(IO%IU6,'(//A)') ' exact exchange contribution, plane wave part only'
!    ENDIF
!    CALL WRITE_EIGENVAL( WDES, W_TMP, IO%IU6)

  END SUBROUTINE SET_EIGENVALUE_HARTREE_KINETIC

!**********************************************************************
!
!>This routine returns the band gap, more precisely the 
!!minimum transition energy beween occupied and unoccupied states
!!E1 is set to the minimal gap
!!E2 is set to the maximum transition energy (passed by calling routine)
!
!**********************************************************************

  SUBROUTINE DETERMINE_BAND_GAP(WDES, W, NB, E1, E2, EFERMI, NOMEGA)
    USE prec
    USE wave
    IMPLICIT NONE
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    REAL(q)::          E1,E2      ! E1 ... band gap
    REAL(q)::          EFERMI     ! Fermi-level
    INTEGER NB                    ! maximum band to be considered
    INTEGER             :: NOMEGA

    ! local
    INTEGER NK, N, ISP, NKP, NP
    REAL(q) :: EMAXO, EMINU, ENEW

    PROFILING_START( 'determine_band_gap' )
    EMAXO=-1000
    EMINU= 1000

    DO ISP=1,WDES%ISPIN
       DO NK=1,WDES%NKPTS
          DO N=1,NB
            IF ( ABS(W%FERTOT(N,NK,ISP))>0.5 ) THEN   
               ! largest occupied energy state
               EMAXO=MAX(EMAXO,REAL(W%CELTOT(N,NK,ISP),q))
            ELSE
               !  smallest unoccupied energy state
               EMINU=MIN(EMINU,REAL(W%CELTOT(N,NK,ISP),q))                
            ENDIF
          ENDDO
       ENDDO
    ENDDO

    E1=(EMINU-EMAXO)

    ! for RPA calculations we need the maximum and minimum of |e_nk - e_n'k'|
    ! now more carefully scan considering all allowed transitions
    E1=1000
    DO ISP=1,WDES%ISPIN
       DO NKP=1,WDES%NKPTS
       DO NK=NKP,WDES%NKPTS
          DO NP=1,NB
          DO N=1,NB
            IF ( ABS(W%FERTOT(N,NK,ISP)-W%FERTOT(NP,NKP,ISP))>0.02)  THEN
               ! this is geared towards the new chi_GG routine, where transitions
               ! are limited to occur from the Fermi-level to occupied or unoccupied states
               !
               ! for insulators the first line always applies
               ! for metals transitions between occupied-occupied or unoccupied-unoccupied
               ! occur but then those are limited by the other two lines
               ENEW=MAX(ABS(REAL(W%CELTOT(N,NK,ISP),q)-REAL(W%CELTOT(NP,NKP,ISP),q)), & 
                     MIN(ABS(REAL(W%CELTOT(N,NK,ISP),q)-EFERMI), &
                         ABS(REAL(W%CELTOT(NP,NKP,ISP),q)-EFERMI)))
               E1=MIN(E1,ENEW)
            ENDIF
          ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO

    ! GW requires to integrate even smaller energies 
    IF ( LGW ) E1 = E1/2

    PROFILING_STOP( 'determine_band_gap' )
    RETURN
  ENDSUBROUTINE DETERMINE_BAND_GAP

END MODULE xi

!**********************************************************************
!
! querry functions for GW
!
!**********************************************************************

FUNCTION CALCULATE_RESPONSE_FUNCTIONS()
  USE xi
  LOGICAL CALCULATE_RESPONSE_FUNCTIONS
  CALCULATE_RESPONSE_FUNCTIONS=LCHI
END

FUNCTION USE_OEP_IN_GW()
  USE xi
  LOGICAL USE_OEP_IN_GW
  USE_OEP_IN_GW=LOEP
END

FUNCTION ENCUTGW_IN_CHI()
  USE xi
  IMPLICIT NONE
  REAL(q) ENCUTGW_IN_CHI
  ENCUTGW_IN_CHI=ENCUTGW
END

FUNCTION FAST_FOCK()
  USE xi
  LOGICAL FAST_FOCK
  FAST_FOCK=.FALSE.
  IF (ICHIREAL>0) THEN
    FAST_FOCK=.TRUE.
  ENDIF
END

FUNCTION CALCULATE_RPA_FORCES()
  USE xi
  LOGICAL CALCULATE_RPA_FORCES
  CALCULATE_RPA_FORCES=LRPAFORCE
END

!********************** SUBROUTINE PUSH_XC_TYPE_FOR_GW*****************
!
! push all parameters read for the xc onto the stack 
! these parameters are potentially required 
! ) for the evaluation of the screened interaction in the local field effects
!  restoring them upon entry of the GW is sufficient (POP_XC)
! ) for the all-in-one new RPA code, in order to perform the DFT
!  calculations with the parameters specified in the INCAR file
!
! Further comments (gK: not entirely sure whether this still applies for vasp.6):
!  if this routine is undocumented the INCAR selected functional
!  is used inside the PAW spheres and for the core-valence interaction
!  if this routine is used, core valence interaction and interaction in
!  the PAW spheres are evaluated using HF without correlation
!
!**********************************************************************

SUBROUTINE PUSH_XC_TYPE_FOR_GW
  USE fock
  USE xi
  USE setexm

  IF (LCHI) THEN
     MODEL_GW=0
     IF (LRSCOR) THEN
     ! long range HF is really an exception, since here we keep the original
     ! correlation functional as it was supplied by the user
        CALL PUSH_XC_TYPE(LEXCH, 0.0_q, ALDAC, 0.0_q, 0.0_q, 0.0_q)
     ELSE
        CALL PUSH_XC_TYPE(LEXCH, 0.0_q, 0.0_q, 0.0_q, 0.0_q, 0.0_q)
     ENDIF
  ENDIF
      
END SUBROUTINE PUSH_XC_TYPE_FOR_GW
