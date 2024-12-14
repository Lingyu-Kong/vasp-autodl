#include "symbol.inc"

!*********************************************************************
!
! this module implements forces for the random phase approximation
! the key equations are presented in 
!  B. Ramberger, T. Schäfer, G. Kresse, G. (2017). 
!   Analytic interatomic forces in the random phase approximation. 
!  Physical review letters, 118(10), 106403.
!  https://doi.org/10.1103/PhysRevLett.118.106403
!  (in particular the supplementary)
!  Many details are also presented in the thesis of B. Ramberger
!
!*********************************************************************

MODULE rpa_force
  USE base
  USE chi_glb
  USE wave_high
  IMPLICIT NONE
  
  LOGICAL :: LINIT_RPA_FORCE = .FALSE.

!
! the force routine needs to know the number of bands used for
! the ground state KS calculations 
!

  TYPE (wavedes)     WDES_GROUNDSTATE

!
! the calculate forces this module stores a number of 
! variables that are calculated on the fly in the RPA routines
! 
! RPA density matrix calculated as 
! n = [ int  d w G_DFT(w) ( T + V_loc + V_nl + V_HF + Sigma(w)) G_DFT(w) ] (r)
! the corresponding contribtutions to the density matrix are stored in
!
  OVERLAP, ALLOCATABLE :: CRHODE_RPA(:,:,:,:)! first order change of PAW occupancies
  COMPLEX(q),ALLOCATABLE:: CHDEN_RPA(:,:)    ! charge-density in reciprocal space on coarse grid
  COMPLEX(q),ALLOCATABLE:: CHTOT_RPA(:,:)    ! charge-density in reciprocal space on dense grid
!
! RPA forces, the two arrays below will be passed to symmetrization routine
! 
  REAL(q), ALLOCATABLE :: DE_LOCAL_RPA(:,:)      ! energy change from local contribution upon displacement RPA_DISPLACEMENT
  REAL(q), ALLOCATABLE :: DE_RPA_ONEC(:,:)       ! energy change from one-center terms upon displacement RPA_DISPLACEMENT
  REAL(q), ALLOCATABLE :: RPA_DISPLACEMENT(:,:,:)! considered displacements for each ion
                                                 ! compare ND in sydmat.F and linear_response.F 
  INTEGER, ALLOCATABLE :: RPA_NDISPLACEMENTS(:)  ! number of displacement for each ion
                                                 ! compared D in sydmat.F and linear_reasponse.F
! other "standard" force components
  REAL(q), ALLOCATABLE :: FORCE_RPA_NL(:,:)
  REAL(q), ALLOCATABLE :: FORCE_GAMMA(:,:)
  REAL(q), ALLOCATABLE :: FORCE_RPA_PULAY(:,:)
!
! RPA orbitals calculated as sum over unoccupied states a 
! \sum_a |a> <a| [ int d w G_DFT(w) ( T + V_loc + V_nl + V_HF + Sigma(w)) G_DFT(w) ]|i> 
  TYPE (wavespin), ALLOCATABLE :: W_RPA      ! conjugated "RPA" orbitals (|i> \sum_a gamma_ia |a>)
  TYPE (wavespin), POINTER   :: W_S=>NULL()  ! conjugated overlap orbitals (|i> \sum_a S_ia |a>)
  
!
! arrays to store first order change of self-consistent Kohn-Sham potential
  COMPLEX(q), ALLOCATABLE :: CVTOT_RPA(:,:)! derivative of local potential (dense grid)
  OVERLAP,ALLOCATABLE :: CDIJ_RPA(:,:,:,:) ! derivative of non local strength
  REAL(q)             :: FORCE_ADD_RPA     ! additional component to force from one-centre terms

!
! approximate second derivative 
! set by calling the routine  STORE_RPA_SECOND_DERIV
  REAL(q),ALLOCATABLE :: SECOND_DERIV_RPA(:,:)
  
  CONTAINS
#ifndef scaLAPACK
!
! dummy routines if scaLAPACK is not available
!

  SUBROUTINE STORE_WDES_GROUNDSTATE(WDES)
      TYPE (wavedes)     WDES
  END SUBROUTINE STORE_WDES_GROUNDSTATE


  SUBROUTINE RPA_ZERO_POT( CVTOT1, CDIJ1)
    
    COMPLEX(q) :: CVTOT1(:,:)! derivative of local potential (dense grid)
    OVERLAP,ALLOCATABLE :: CDIJ1(:,:,:,:)! derivative of non local strength
  END SUBROUTINE RPA_ZERO_POT

  SUBROUTINE RPA_ADD_TO_POT( CVTOT1, CDIJ1, FORCE, WEIGHT)

    REAL(q)    :: WEIGHT
    COMPLEX(q) :: CVTOT1(:,:)! derivative of local potential (dense grid)
    OVERLAP,ALLOCATABLE :: CDIJ1(:,:,:,:)! derivative of non local strength
    REAL(q)    :: FORCE      ! additional component to force from mean field 

  END SUBROUTINE RPA_ADD_TO_POT


  SUBROUTINE RPA_FORCE_CHARGE( WDES, T_INFO, GRIDC, ION, IDIR)
    USE poscar
    TYPE (wavedes)     WDES    
    TYPE (type_info)   T_INFO
    TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
    INTEGER            ION        ! ion that was displaced
    INTEGER            IDIR       ! cartesian direction into which ion was displaced
  END SUBROUTINE RPA_FORCE_CHARGE


  SUBROUTINE STORE_RPA_SECOND_DERIV(SECOND_DERIV )
    REAL(q) :: SECOND_DERIV(:,:)
  END SUBROUTINE STORE_RPA_SECOND_DERIV


  SUBROUTINE WRITE_RPA_FORCE( TIFOR, T_INFO, LATT_CUR, SYMM, POSIOC, IBRION, IO )
    USE poscar
    USE lattice
    REAL(q) :: TIFOR(:,:)
    TYPE (type_info)   T_INFO
    TYPE (latt)        LATT_CUR
    TYPE (symmetry)    SYMM
    REAL(q) :: POSIOC(:,:)      ! updated ionic positions (written to CONTCAR)
    INTEGER :: IBRION           ! selected algorithm for ions
    TYPE (in_struct)   IO
  END SUBROUTINE WRITE_RPA_FORCE

  SUBROUTINE RPA_DENSITY( W, LMDIM, IRDMAX, T_INFO, LATT_CUR, P, SYMM, &
       CDIJ, CQIJ, NONLR_S, NONL_S, &
       GRID, GRID_SOFT, GRIDC, GRIDUS, SOFT_TO_C, C_TO_US) 
    USE poscar
    USE lattice
    USE pseudo
    USE charge
    USE us
    USE nonl_high

    TYPE (wavespin), TARGET ::    W
    INTEGER LMDIM 
    INTEGER IRDMAX
    TYPE (type_info)   T_INFO
    TYPE (latt)        LATT_CUR
    TYPE (potcar)      P(T_INFO%NTYP)
    TYPE (symmetry)    SYMM
    OVERLAP  CDIJ(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ), &
             CQIJ(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (grid_3d)     GRID       ! grid for wavefunctions
    TYPE (grid_3d)     GRID_SOFT  ! grid for soft chargedensity
    TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
    TYPE (grid_3d)     GRIDUS     ! temporary grid in us.F

    TYPE (transit)     SOFT_TO_C  ! index table between GRID_SOFT and GRIDC
    TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
  END SUBROUTINE RPA_DENSITY

#else
!************************ SUBROUTINE STORE_NBANDS_GROUNDSTATE **********
!
! this functions stores the number of bands used in the 
! original groundstate KS calculations
! that is usually the default value used by VASP
! 
!***********************************************************************

    SUBROUTINE STORE_WDES_GROUNDSTATE(WDES)
      TYPE (wavedes)     WDES
      LINIT_RPA_FORCE=.TRUE.
      WDES_GROUNDSTATE=WDES

    END SUBROUTINE STORE_WDES_GROUNDSTATE


!************************ SUBROUTINE STORE_RPA_PULAY_FORCES **********
!
! this functions stores the RPA Pulay forces as calculated in chi_GG
! 
!***********************************************************************

    SUBROUTINE STORE_FORCE_RPA_PULAY(T_INFO,FRNL)
      USE poscar
      TYPE (type_info)   T_INFO
      REAL(q) :: FRNL(3,T_INFO%NIONS)
      
      IF (.NOT. LINIT_RPA_FORCE) RETURN

      IF (ALLOCATED(FORCE_RPA_PULAY))  DEALLOCATE(FORCE_RPA_PULAY)
      ALLOCATE(FORCE_RPA_PULAY(3,T_INFO%NIONS))
      FORCE_RPA_PULAY=0
      FORCE_RPA_PULAY=FRNL

    END SUBROUTINE STORE_FORCE_RPA_PULAY

!***********************************************************************
! 
! according to Equ. (12) in PRL 118, 106403 (2017)
!                                _
! The forces are calculated as _/  Tr{G(F + Sigma)G ( dH^KS/dR + w.dS/dR)} dw 
! where F is the Fock operator and Sigma the correlation part of the self-energy
!
! STORE_RPA_DENSITY_MATRIX is called from chi_super.F and the calling
! routine needs to pass down the density matrix in the set of KS orbitals
!  CHAM_MAT  = gammaH_nm + delta_nm fermi_weight_n
!  SIGMA_MAT = gammaS_nm 
!
! the routine sets up two sets of orbitals W_RPA and W_S
! which are then used by the routine RPA_DENSITY to calculate 
! the first order change of the density rho^(1) and gamma^(1)
! [Equ. (13) and (15) in PRL 118, 106403 (2017)]
! as well as the corresponding non-local contributions
!
! the density matrices rho^1 and gamma^1 [Equ. (13) and (15) in PRL 118, 106403 (2017)]
! are given by
!   rho^1   = sum_nm |n> gammaH_nm <m|
!   gamma^1 = sum_nm |n> gammaS_nm <m|
! where gamma_nm is calculated by the calling routine as
!
!               < n | F | m >
! gammaH_nm=  ----------------.(f_n-f_m) + CORR_MAT_nm
!                e_n  - e_m
!             \..........  ........./
!                        \/
!       = 0 if n,m are both occ or both uocc
!
! as gammaS_nm is given by
!
!                 < n | F | m >             
! gammaS_nm= e_n.---------------- .(f_n-f_m) + f_n.f_m.< n | F | m > +
!                  e_n  - e_m    
!                  _             
!           + <n|_/  dw w.G.\Sigma.G |m> 
!
! To efficiently determine rho^(1)(r) = Tr[ |r><r| rho^(1)] 
! the orbitals in W_RPA are set to
!  |n^RPA > = \sum_m  gammaH_mn |m>
!  |n^S >   = \sum_m  gammaS_mn |m>
!
! hence the first-order change of the density can be obtained by
! summing over the *occupied* orbitals only:
!   rho^1(r)   =  sum_n <n | r><r | n^RPA> 
!   gamma^1(r) =  sum_n <n | r><r | n^S> 
!
! so |n^RPA> can be interpreted as the first-order change of the orbital |n>
!
!***********************************************************************

  SUBROUTINE STORE_RPA_DENSITY_MATRIX(W, WDES, CHAM_MAT, COMM_INTAU, DESC, SIGMA_MAT )
    USE dfast
    USE scala
    IMPLICIT NONE
    TYPE (wavespin)    W
    TYPE (wavedes)     WDES
    GDEF, POINTER :: CHAM_MAT(:,:,:,:)     ! density matrix rho^(1) in KS orbital basis
    GDEF, POINTER, OPTIONAL :: SIGMA_MAT(:,:,:,:) ! overlap density matrix gamma^(1) in KS orbital basis
    TYPE (communic) COMM_INTAU
    INTEGER,PARAMETER :: DLEN_=9
    INTEGER              :: DESC( DLEN_ )  ! distributed matrix descriptor array

  ! local
    TYPE (wavedes1)    WDES1               ! descriptor for one k-point
    TYPE (wavefuna)    WA                  ! subpointer to part of W
    INTEGER ISP, NK

    IF (.NOT. LINIT_RPA_FORCE) RETURN

    IF (ALLOCATED(W_RPA)) THEN
       CALL DEALLOCW(W_RPA)
       DEALLOCATE(W_RPA)
    ENDIF
    ALLOCATE(W_RPA)
    CALL ALLOCW(WDES,W_RPA)

    IF (SIZE(W_RPA%CW,1) /= SIZE(W%CW,1)) THEN
       CALL vtutor%bug("internal error in STORE_RPA_DENSITY_MATRIX: number of coefficients not " // &
          "correct " // str(SIZE(W_RPA%CW,1)) // " " // str(SIZE(W%CW,1)), __FILE__, __LINE__)
    ENDIF

    IF (SIZE(W_RPA%CW,2) /= WDES%NBANDS) THEN
       CALL vtutor%bug("internal error in STORE_RPA_DENSITY_MATRIX: number of bands not correct " &
          // str(SIZE(W_RPA%CW,2)) // " " // str(WDES%NBANDS), __FILE__, __LINE__)
    ENDIF

    CALL SETWDES(WDES,WDES1, 0)
    CALL NEWWAVA(WA, WDES1, WDES%NBANDS)

    DO ISP=1,WDES%ISPIN
       DO NK=1,WDES%NKPTS
          CALL SETWDES(WDES,WDES1,NK)
          ! the calling routine supplies the full density matrix, subtract unperturbed groundstate KS density matrix
          CALL  ADD_TO_DIAGONALE_REAL( WDES%NB_TOTK(NK,ISP), CHAM_MAT(:,NK,ISP,1), -W%FERTOT(:,NK,ISP), DESC)

          WA%CW   =W%CW(:,:, NK, ISP)
          WA%CPROJ=W%CPROJ(:,:, NK, ISP)

          !  redistribution over bands
          CALL REDISTRIBUTE_PROJ( WA)
          CALL REDISTRIBUTE_PW( WA)

          CALL LINCOM_DISTRI_DESC('F',WA%CW_RED(1,1),WA%CPROJ_RED(1,1), CHAM_MAT(1,NK,ISP,1), &
               WDES1%NB_TOTK(ISP), & 
               WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,WDES%NB_TOT, &
               COMM_INTAU, NBLK, DESC)

          !  back-redistribution over plane wave coefficients
          CALL REDISTRIBUTE_PROJ( WA)
          CALL REDISTRIBUTE_PW( WA)
          W_RPA%CW   (:, 1:WDES%NBANDS, NK, ISP)=WA%CW   (:, 1:WDES%NBANDS)
          W_RPA%CPROJ(:, 1:WDES%NBANDS, NK, ISP)=WA%CPROJ(:, 1:WDES%NBANDS)
       ENDDO
    ENDDO

    W_RPA%FERTOT(1:WDES%NBANDS, :, :)=0
    W_RPA%CELTOT(1:WDES%NBANDS, :, :)=W%CELTOT(1:WDES%NBANDS, :, :)
!
! now rotate the orbitals using SIGMA_MAT and store results in W_S
!
    IF ( PRESENT(SIGMA_MAT)) THEN
       IF (ASSOCIATED( W_S)) THEN
          CALL DEALLOCW(W_S)
          DEALLOCATE(W_S)
       ENDIF
       ALLOCATE(W_S)
       CALL ALLOCW(WDES,W_S)

       DO ISP=1,WDES%ISPIN
       DO NK=1,WDES%NKPTS
          CALL SETWDES(WDES,WDES1,NK)

          WA%CW   =W%CW(:,:, NK, ISP)
          WA%CPROJ=W%CPROJ(:,:, NK, ISP)

          !  redistribution over bands
          CALL REDISTRIBUTE_PROJ( WA)
          CALL REDISTRIBUTE_PW( WA)

          CALL LINCOM_DISTRI_DESC('F',WA%CW_RED(1,1),WA%CPROJ_RED(1,1), SIGMA_MAT(1,NK,ISP,1), &
               WDES1%NB_TOTK(ISP), & 
               WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,WDES%NB_TOT, &
               COMM_INTAU, NBLK, DESC)

          !  back-redistribution over plane wave coefficients
          CALL REDISTRIBUTE_PROJ( WA)
          CALL REDISTRIBUTE_PW( WA)
          W_S%CW   (:, 1:WDES%NBANDS, NK, ISP)=WA%CW   (:, 1:WDES%NBANDS)
          W_S%CPROJ(:, 1:WDES%NBANDS, NK, ISP)=WA%CPROJ(:, 1:WDES%NBANDS)
       ENDDO
       ENDDO
       W_S%FERTOT(1:WDES%NBANDS, :, :)=0
       W_S%CELTOT(1:WDES%NBANDS, :, :)=W%CELTOT(1:WDES%NBANDS, :, :)
    ENDIF

    CALL DELWAVA(WA)

  END SUBROUTINE STORE_RPA_DENSITY_MATRIX


!***********************************************************************
!
! RPA_DENSITY is called from main.F
! - it calculates the RPA charge density defined as 
!    rho^1(r)   =  sum_n <n | r><r |n^RPA> 
! - the non-local terms to the forces
!    F  =  sum_n <n | d p_i/ d R > D_ij + epsilon_n Q_ij < p_j| n^RPA> +
!          sum_n <n | p_i > D_ij + epsilon_n Q_ij < d p_j/ d R | n^RPA> +
!        + sum_n <n | d p_i/ d R > Q_ij < p_j| n^S> +
!          sum_n <n | p_i > Q_ij < d p_j/ d R | n^S> +
!  where p are the projection operators and D_ij the non-local
!  strength operators in DFT
!
! the result is stored in CHTOT_RPA, CHDEN_RPA etc.
!
!***********************************************************************

  SUBROUTINE RPA_DENSITY( W, LMDIM, IRDMAX, T_INFO, LATT_CUR, P, SYMM, &
       CDIJ, CQIJ, NONLR_S, NONL_S, &
       GRID, GRID_SOFT, GRIDC, GRIDUS, SOFT_TO_C, C_TO_US) 
    USE poscar
    USE lattice
    USE pseudo
    USE charge
    USE us
    USE nonl_high

    TYPE (wavespin), TARGET ::    W
    INTEGER LMDIM 
    INTEGER IRDMAX
    TYPE (type_info)   T_INFO
    TYPE (latt)        LATT_CUR
    TYPE (potcar)      P(T_INFO%NTYP)
    TYPE (symmetry)    SYMM
    OVERLAP  CDIJ(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ), &
             CQIJ(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (grid_3d)     GRID       ! grid for wavefunctions
    TYPE (grid_3d)     GRID_SOFT  ! grid for soft chargedensity
    TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
    TYPE (grid_3d)     GRIDUS     ! temporary grid in us.F

    TYPE (transit)     SOFT_TO_C  ! index table between GRID_SOFT and GRIDC
    TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
  ! local variables
    INTEGER :: ISP, NG
    REAL(q) DISPL(3,T_INFO%NIONS)
    TYPE (wavedes), POINTER :: WDES
    OVERLAP, ALLOCATABLE :: CRHODE0(:,:,:,:)

    IF (.NOT. LINIT_RPA_FORCE) RETURN

    ! somewhat tricky, ISYM defaults to 3 for HF type calculations
    ! in this case, orbitals are brought to all k-points in the full
    ! wedge and symmetry is not applied to density related quantities
    ! here, however, we need to symmetrize the densities and one-centre occupancies 
    ! this only happens for ISYM=2
    IF (SYMM%ISYM>0) SYMM%ISYM=2

    WDES=>W%WDES

    IF (ALLOCATED(CRHODE_RPA)) DEALLOCATE(CRHODE_RPA)
    IF (ALLOCATED(CHDEN_RPA))  DEALLOCATE(CHDEN_RPA)
    IF (ALLOCATED(CHTOT_RPA))  DEALLOCATE(CHTOT_RPA)
    IF (ALLOCATED(FORCE_RPA_NL))  DEALLOCATE(FORCE_RPA_NL)
    IF (ALLOCATED(FORCE_GAMMA))  DEALLOCATE(FORCE_GAMMA)
    IF (ALLOCATED(DE_LOCAL_RPA))  DEALLOCATE(DE_LOCAL_RPA)
    IF (ALLOCATED(DE_RPA_ONEC))  DEALLOCATE(DE_RPA_ONEC)
    IF (ALLOCATED(RPA_DISPLACEMENT))  DEALLOCATE(RPA_DISPLACEMENT)
    IF (ALLOCATED(RPA_NDISPLACEMENTS))  DEALLOCATE(RPA_NDISPLACEMENTS)
    ALLOCATE(CHTOT_RPA(GRIDC%MPLWV,WDES%NCDIJ),  &
             CHDEN_RPA(GRID_SOFT%MPLWV,WDES%NCDIJ), &
             CRHODE_RPA(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), &
             FORCE_RPA_NL(3,T_INFO%NIONS), &
             FORCE_GAMMA(3,T_INFO%NIONS), &
             DE_LOCAL_RPA(3,T_INFO%NIONS), &
             DE_RPA_ONEC(3,T_INFO%NIONS), &
             RPA_DISPLACEMENT(3,3, T_INFO%NIONS), &
             RPA_NDISPLACEMENTS( T_INFO%NIONS))

    RPA_DISPLACEMENT=0
    RPA_NDISPLACEMENTS=0

    FORCE_RPA_NL=0
    FORCE_GAMMA=0
    DE_LOCAL_RPA=0
    DE_RPA_ONEC=0

    CHDEN_RPA=0
    CHTOT_RPA=0
    ! determine first order change of one center occupancy matrix 
    CALL DEPSUM1_RPA(W, W_RPA, WDES, LMDIM, CRHODE_RPA, WDES%LOVERL)
    ! change storage convention to (total, magnetization)
    CALL US_FLIP(WDES, LMDIM, CRHODE_RPA, WDES%LOVERL, .FALSE.)

    ! symmetrize CRHODE_RPA (currently only done by DEPLE_ADD, alas...)
    ! so just call DEPLE_ADD and disregard CHTOT_RPA
    CALL DEPLE_ADD(WDES,GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
         LATT_CUR,P,T_INFO,SYMM, WDES%LOVERL, SOFT_TO_C,&
         LMDIM, CRHODE_RPA, CHTOT_RPA, CHDEN_RPA, IRDMAX)
    CHTOT_RPA=0

    ! soft pseudo charge
    CALL SOFT_CHARGE1_RPA(GRID,GRID_SOFT,W,W_RPA,WDES, CHDEN_RPA)
    ! change storage convention to (total, magnetization)
    CALL RC_FLIP(CHDEN_RPA, GRID_SOFT, WDES%NCDIJ, .FALSE.)

   ! symmetrisation of soft pseudo charge
    IF (SYMM%ISYM==2) THEN
       IF (WDES%LNONCOLLINEAR) THEN
          CALL RHOSYM(CHDEN_RPA(1,1),GRID_SOFT,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,1)
          IF (.NOT.WDES%LSPIRAL) &
               CALL SYMFIELD(CHDEN_RPA(1,2),GRID_SOFT,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,WDES%SAXIS,LATT_CUR)
       ELSE
          DO ISP=1,WDES%ISPIN
             CALL RHOSYM(CHDEN_RPA(1,ISP),GRID_SOFT,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,ISP)
          ENDDO
       ENDIF
    ENDIF
    
    ALLOCATE(CRHODE0(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ))
    ! TODO: eventually, only the contribution CHDEN_RPA should be used,
    ! and CHTOT_RPA should be removed from the code
    ! but this will be done later, and will change only little in terms of
    ! performance 

    CRHODE0=0
    ! go from coarse grid to fine grid to obtain CHTOT_RPA
    ! (augmentation charges are NOT added and treated seperately)
    CALL DEPLE_ADD(WDES,GRID_SOFT,GRIDC,GRIDUS,C_TO_US, &
         LATT_CUR,P,T_INFO,SYMM, WDES%LOVERL, SOFT_TO_C,&
         LMDIM, CRHODE0, CHTOT_RPA, CHDEN_RPA, IRDMAX)
    DEALLOCATE(CRHODE0)
    
    ! symmetrise total pseudo charge density CHTOT as well 
    IF (SYMM%ISYM==1) THEN
       IF (WDES%LNONCOLLINEAR) THEN
          CALL RHOSYM(CHTOT_RPA(1,1),GRIDC,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,1)
          IF (.NOT.WDES%LSPIRAL) &
               &   CALL SYMFIELD(CHTOT_RPA(1,2),GRIDC,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,WDES%SAXIS,LATT_CUR)
       ELSE
          DO ISP=1,WDES%ISPIN
             CALL RHOSYM(CHTOT_RPA(1,ISP),GRIDC,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,ISP)
          ENDDO
       ENDIF
    ENDIF

    DISPL=0
    IF (NONLR_S%LREAL) THEN
       IF (ASSOCIATED(W_S)) THEN
          CALL vtutor%bug("internal error in VASP: RPA_DENSITY does not yet support W_S for real space", __FILE__, __LINE__)
       ENDIF
       CALL FORNLR2(GRID,NONLR_S,P,LATT_CUR, &
            &   W, W_RPA, CDIJ, CQIJ, DISPL, FORCE_RPA_NL)
    ELSE
       CALL FORNL2(NONL_S,WDES, W, W_RPA, W_S, LATT_CUR, LMDIM,CDIJ,CQIJ, FORCE_RPA_NL, FORCE_GAMMA)
    ENDIF

! at this point we do not need W_RPA anymore, so deallocate it
    CALL DEALLOCW(W_RPA)
    DEALLOCATE(W_RPA)
    IF (ASSOCIATED(W_S)) THEN
       CALL DEALLOCW(W_S)
       DEALLOCATE(W_S)
       NULLIFY(W_S)
    ENDIF
    
  END SUBROUTINE RPA_DENSITY

!***********************************************************************
!
!  the followin subroutines allocates the necessary arrays to accumulate
!  the first order change of the local Kohn-Sham potential
!   CVTOT_RPA   potential on fine grid
!   CDIJ_RPA    non local strength parameters
!  the routine RPA_ADD_TO_POT can be called multiple times, for instance
!  when using finite differences with appropriate weights 
!
!***********************************************************************

  SUBROUTINE RPA_ZERO_POT( CVTOT1, CDIJ1)
    
    COMPLEX(q) :: CVTOT1(:,:)   ! derivative of local potential (dense grid)
    OVERLAP    :: CDIJ1(:,:,:,:)! derivative of non local strength

    IF (.NOT. LINIT_RPA_FORCE) RETURN

    IF (ALLOCATED(CVTOT_RPA)) THEN
       DEALLOCATE(CVTOT_RPA)
    ENDIF

    IF (ALLOCATED(CDIJ_RPA)) THEN
       DEALLOCATE(CDIJ_RPA)
    ENDIF

    ALLOCATE(CVTOT_RPA( SIZE(CVTOT1,1), SIZE(CVTOT1,2)))
    ALLOCATE(CDIJ_RPA(SIZE(CDIJ1,1), SIZE(CDIJ1,2), SIZE(CDIJ1,3), SIZE(CDIJ1,4)))
    
    CVTOT_RPA=0
    CDIJ_RPA=0
    FORCE_ADD_RPA=0

  END SUBROUTINE RPA_ZERO_POT

  SUBROUTINE RPA_ADD_TO_POT( CVTOT1, CDIJ1, FORCE, WEIGHT)

    REAL(q)    :: WEIGHT        ! weight with which the contribution is added
    COMPLEX(q) :: CVTOT1(:,:)   ! derivative of local potential (dense grid)
    OVERLAP    :: CDIJ1(:,:,:,:)! derivative of non local strength
    REAL(q)    :: FORCE         ! additional component to force from mean field 

    IF (.NOT. LINIT_RPA_FORCE) RETURN

    IF (.NOT. ALLOCATED(CVTOT_RPA)) THEN
       CALL vtutor%bug("internal error in VASP: RPA_ADD_TO_POT, CVTOT_RPA is not allocated", __FILE__, __LINE__)
    ENDIF
    IF (.NOT. ALLOCATED(CDIJ_RPA)) THEN
       CALL vtutor%bug("internal error in VASP: RPA_ADD_TO_POT, CDIJ_RPA is not allocated", __FILE__, __LINE__)
    ENDIF

    CVTOT_RPA=CVTOT_RPA+CVTOT1*WEIGHT
    CDIJ_RPA=CDIJ_RPA+CDIJ1*WEIGHT
    FORCE_ADD_RPA=FORCE_ADD_RPA+FORCE*WEIGHT

  END SUBROUTINE RPA_ADD_TO_POT

!****************** subroutine RPA_FORCE_CHARGE ************************
!
! calculate force contributions related to the change of the 
! Kohn-Sham potential
!
!***********************************************************************


  SUBROUTINE RPA_FORCE_CHARGE( WDES, T_INFO, GRIDC, ION, IDIR)
    USE poscar
    TYPE (wavedes)     WDES    
    TYPE (type_info)   T_INFO
    TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
    INTEGER            ION        ! ion that was displaced
    INTEGER            IDIR       ! cartesian direction into which ion was displaced
  ! local
    INTEGER ISP, NG, NC, N1, N2, N3
    COMPLEX(q) :: FOR, FOR_NL
    REAL(q)    :: FACTM
    INTEGER    :: ISPINOR, ISPINOR_, LMBASE, LMBASE_, NIS, NT, LMMAXC, NI, L, LP

    IF (.NOT. LINIT_RPA_FORCE) RETURN

  ! go from real space to reciprocal space for potential
     DO ISP=1,WDES%NCDIJ
        CALL FFT_RC_SCALE(CVTOT_RPA(1,ISP), CVTOT_RPA(1,ISP),GRIDC)
     ENDDO

  ! go to (total potential, magnetization representation)
     IF (SIZE(CVTOT_RPA,1) /=GRIDC%MPLWV) THEN
        WRITE(*,*) 'internal error in VASP: CVTOT_RPA allocated in non-compatible manner',SIZE(CVTOT_RPA,1), GRIDC%MPLWV
     ENDIF
     CALL RC_FLIP_POTENTIAL(CVTOT_RPA, GRIDC, WDES%NCDIJ, .FALSE.)

  ! now contract potential and charge density
  ! these are identical to CVTOT_RPA and CHTOT_RPA, but on a coarser grid
     FOR=0
     DO ISP=1,WDES%NCDIJ
        ! loop over reciprocal lattice vectors
         DO NG=1,GRIDC%RC%NP
           N1= MOD((NG-1),GRIDC%RC%NROW) +1
           NC= (NG-1)/GRIDC%RC%NROW+1
           N2= GRIDC%RC%I2(NC)
           N3= GRIDC%RC%I3(NC)

           SETFACT1
           SETFACT
           FOR=FOR+CVTOT_RPA(NG,ISP)*MULFACT CONJG(CHTOT_RPA(NG,ISP))
        ENDDO
     ENDDO
   ! sum all cores
     CALLMPI( M_sum_z(GRIDC%COMM, FOR, 1))

   ! now contract RPA one-center density matrix against change of one-centre strenght
     FOR_NL=0
     IF (WDES%LOVERL) THEN
       spin:   DO ISP=1,WDES%ISPIN
       spinor: DO ISPINOR =0,WDES%NRSPINORS-1
       DO ISPINOR_=0,WDES%NRSPINORS-1

         LMBASE =ISPINOR *(WDES%NPRO/2)
         LMBASE_=ISPINOR_*(WDES%NPRO/2)

         NIS   =1
         typ:  DO NT=1,WDES%NTYP
         LMMAXC=WDES%LMMAX(NT)
         IF (LMMAXC==0) GOTO 210
         
         DO NI=NIS,WDES%NITYP(NT)+NIS-1
         
!DIR$ IVDEP
!OCL NOVREC
         DO L=1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
            DO LP=1,LMMAXC
               FOR_NL=FOR_NL+ CRHODE_RPA(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)*CHGCONJG(CDIJ_RPA(LP,L,NI,ISP+ISPINOR_+2*ISPINOR))
            ENDDO
         ENDDO
         LMBASE = LMMAXC+LMBASE
         LMBASE_= LMMAXC+LMBASE_
         ENDDO

 210     NIS = NIS+WDES%NITYP(NT)
       ENDDO typ
       ENDDO
       ENDDO spinor 
       ENDDO spin
     ENDIF

   ! sum all cores
     CALLMPI( M_sum_z(WDES%COMM_INB, FOR_NL, 1))

     RPA_NDISPLACEMENTS(ION)=RPA_NDISPLACEMENTS(ION)+1
     IF (RPA_NDISPLACEMENTS(ION)>3) THEN
        CALL vtutor%bug("internal error in VASP: RPA_FORCE_CHARGE displacements larger 3 " // str(ION), __FILE__, __LINE__)
     ENDIF

     ! store displacement vector
     RPA_DISPLACEMENT(IDIR,RPA_NDISPLACEMENTS(ION),ION)=1

     !force is -dE/dR
     DE_LOCAL_RPA(RPA_NDISPLACEMENTS(ION),ION)=-REAL(FOR,q)-FORCE_ADD_RPA
     DE_RPA_ONEC(RPA_NDISPLACEMENTS(ION),ION) =-REAL(FOR_NL,q)
     DEALLOCATE(CVTOT_RPA)

   END SUBROUTINE RPA_FORCE_CHARGE

!****************** subroutine STORE_RPA_SECOND_DERIV ************************
!
! small helper routine to store approximate second derivative
! (e.g. DFT Hessian)
!
!***********************************************************************

   SUBROUTINE STORE_RPA_SECOND_DERIV(SECOND_DERIV )
     REAL(q) :: SECOND_DERIV(:,:)

     IF (.NOT. LINIT_RPA_FORCE) RETURN

     IF (ALLOCATED(SECOND_DERIV_RPA)) THEN
        DEALLOCATE(SECOND_DERIV_RPA)
     ENDIF
     ALLOCATE(SECOND_DERIV_RPA(SIZE(SECOND_DERIV,1),SIZE(SECOND_DERIV,2)))

     ! note sign inversion since the negative Hessian is calculated
     SECOND_DERIV_RPA=-SECOND_DERIV

   END SUBROUTINE STORE_RPA_SECOND_DERIV

!****************** subroutine WRITE_RPA_FORCE ************************
!
! write the final RPA forces to OUTCAR
!
!***********************************************************************

   SUBROUTINE WRITE_RPA_FORCE( TIFOR, T_INFO, LATT_CUR, SYMM, POSIOC, IBRION, IO )
     USE poscar
     USE lattice
     REAL(q) :: TIFOR(:,:)
     TYPE (type_info)   T_INFO
     TYPE (latt)        LATT_CUR
     TYPE (symmetry)    SYMM
     REAL(q) :: POSIOC(:,:)      ! updated ionic positions (written to CONTCAR)
     INTEGER :: IBRION           ! selected algorithm for ions
     TYPE (in_struct)   IO
   ! local    
     INTEGER :: NI
     REAL(q) :: VTMP(3)

     INTEGER ISYMOP, NROT, IGRPOP, NROTK, INVMAP, NPCELL
     REAL(q)  GTRANS,AP
     COMMON /SYMM/   ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
          GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL
     
   ! temporary form 
    REAL(q) :: FORCE_SYM(3,T_INFO%NIONS)
    REAL(q) :: FORCE_RPA(3,T_INFO%NIONS) 
    REAL(q) :: FORCE_RPA_ONEC(3,T_INFO%NIONS) 
    REAL(q), ALLOCATABLE :: HESSIAN(:,:)
    REAL(q) :: STEP(3*T_INFO%NIONS)
   ! test
    REAL(q) :: RPA_DISPLACEMENT_SYM(3,3,T_INFO%NIONS)

    IF (.NOT. LINIT_RPA_FORCE) RETURN

    CALL MKFOR(SYMM%ROTMAP,ISYMOP,INVMAP,NROTK,NPCELL,RPA_DISPLACEMENT,DE_LOCAL_RPA,RPA_NDISPLACEMENTS, &
         1,T_INFO%NTYP,T_INFO%NIONS,T_INFO%NITYP, & 
         LATT_CUR%A(1,1),LATT_CUR%A(1,2),LATT_CUR%A(1,3), &
         LATT_CUR%B(1,1),LATT_CUR%B(1,2),LATT_CUR%B(1,3), IO%IU0 >=0 .AND. .FALSE.)
     
    FORCE_RPA=DE_LOCAL_RPA

    IF (IO%IU6 >=0 ) THEN
        WRITE(IO%IU6,'(A45)') 'local RPA contributions before symmetrization'
        DO NI=1,T_INFO%NIONS
          WRITE(IO%IU6, '(1(3F12.7,2X))') FORCE_RPA(:,NI)
        ENDDO 
    ENDIF

    CALL MKFOR(SYMM%ROTMAP,ISYMOP,INVMAP,NROTK,NPCELL,RPA_DISPLACEMENT,DE_RPA_ONEC,RPA_NDISPLACEMENTS, &
         1,T_INFO%NTYP,T_INFO%NIONS,T_INFO%NITYP, & 
         LATT_CUR%A(1,1),LATT_CUR%A(1,2),LATT_CUR%A(1,3), &
         LATT_CUR%B(1,1),LATT_CUR%B(1,2),LATT_CUR%B(1,3), IO%IU0 >=0 .AND. .FALSE.)
     
    FORCE_RPA_ONEC=DE_RPA_ONEC

    IF (SYMM%ISYM>0) THEN
       ! non-local contributions are only evaluated in IRZ, need to be symmetrized
       CALL FORSYM(FORCE_RPA_NL,SYMM%ROTMAP(1,1,1),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIOND,SYMM%TAUROT,SYMM%WRKROT,LATT_CUR%A)
       CALL FORSYM(FORCE_GAMMA,SYMM%ROTMAP(1,1,1),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIOND,SYMM%TAUROT,SYMM%WRKROT,LATT_CUR%A)
       CALL FORSYM(FORCE_RPA_PULAY,SYMM%ROTMAP(1,1,1),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIOND,SYMM%TAUROT,SYMM%WRKROT,LATT_CUR%A)
       ! just in case, symmetrize the terms obtained by symmetry above
       CALL FORSYM(FORCE_RPA,SYMM%ROTMAP(1,1,1),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIOND,SYMM%TAUROT,SYMM%WRKROT,LATT_CUR%A)
       CALL FORSYM(FORCE_RPA_ONEC,SYMM%ROTMAP(1,1,1),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIOND,SYMM%TAUROT,SYMM%WRKROT,LATT_CUR%A)
    ENDIF

    ! write to stdout 
!      IF (IO%IU0 >=0 .AND. .FALSE. ) THEN
!	DO NI=1,T_INFO%NIONS
!           WRITE(IO%IU0, '(6(3F12.7,2X))') FORCE_RPA(:,NI), FORCE_RPA_NL(:,NI), FORCE_GAMMA(:,NI), FORCE_RPA_ONEC(:,NI), FORCE_RPA_PULAY(:,NI),TIFOR(:,NI), &
!                                          FORCE_RPA(:,NI) + FORCE_RPA_NL(:,NI) + FORCE_GAMMA(:,NI)  + FORCE_RPA_ONEC(:,NI) + FORCE_RPA_PULAY(:,NI) + TIFOR(:,NI) 
!	ENDDO
!      ENDIF

    ! now write to OUTCAR IU6
      IF (IO%IU6 >=0) THEN
        WRITE(IO%IU6,'(6(A24,7X))') 'local RPA contributions','non-local contributions','overlap Q gamma(1)','one-center contributions','Pulay contributions','Hellman-Feynman'
        DO NI=1,T_INFO%NIONS
          WRITE(IO%IU6, '(7(3F10.5,1X))') FORCE_RPA(:,NI), FORCE_RPA_NL(:,NI), FORCE_GAMMA(:,NI), FORCE_RPA_ONEC(:,NI), FORCE_RPA_PULAY(:,NI),TIFOR(:,NI) 
        ENDDO 
        WRITE(IO%IU6,*)
        WRITE(IO%IU6,72)
      ENDIF

      ! update total force
      TIFOR = FORCE_RPA + FORCE_RPA_NL + FORCE_GAMMA + FORCE_RPA_ONEC + FORCE_RPA_PULAY + TIFOR

      IF (IO%IU6 >=0) THEN
        DO NI=1,T_INFO%NIONS
           VTMP=T_INFO%POSION(1:3,NI)
           CALL  DIRKAR(1,VTMP,LATT_CUR%A)
           WRITE(IO%IU6, '(3F13.5,3X,3F14.6)') VTMP, TIFOR(:,NI)
	ENDDO
      ENDIF

 72   FORMAT( ' POSITION    ',35X,'TOTAL RPA FORCE (eV/Angst)'/ &
     &          ' ----------------------------------------------', &
     &          '-------------------------------------')

      VTMP=SUM(TIFOR,2)
      IF (IO%IU6>=0) WRITE(IO%IU6,74) VTMP

 74     FORMAT( ' ----------------------------------------------', &
     &          '-------------------------------------',/ &
     &          '    total drift:      ',20X,3F14.6)


      ! if size of second derivative matrix is fine write out 
      ! updated positions
      IF (ALLOCATED(SECOND_DERIV_RPA)) THEN
      IF (SIZE(SECOND_DERIV_RPA,1) == T_INFO%NIONS*3) THEN

         ! note only IO core has SECOND_DERIV_RPA set and enters here
         ! since only only IO core performs required calculations in linear_response.F
         ALLOCATE(HESSIAN( T_INFO%NIONS*3,  T_INFO%NIONS*3))
         HESSIAN=SECOND_DERIV_RPA
         
         ! use the Hessian matrix from DFT to estimate new optimal positions
         ! essentialy x_new = x - Hessian^-1 force_rpa
         CALL MUL_INV_SECOND_DERIV( HESSIAN, 3*T_INFO%NIONS, RESHAPE(TIFOR, (/ T_INFO%NIONS *3 /)) , STEP)
         ! apply damping using parameter DAMP_NEWTON
         STEP=STEP*DAMP_NEWTON

         IF (IO%IU6 >=0) THEN
            IF (IBRION==-1) THEN
               WRITE(IO%IU6,"(/ ' SUGGESTED UPDATED POSCAR (direct coordinates)  step (CONTCAR is also updated)'/ & 
     &          ' ----------------------------------------------', &
     &          '-------------------------------------')")
            ELSE
               WRITE(IO%IU6,"(/ ' SUGGESTED UPDATED POSCAR (direct coordinates)  step '/ & 
     &          ' ----------------------------------------------', &
     &          '-------------------------------------')")
            ENDIF
            DO NI=1,T_INFO%NIONS
               VTMP=STEP( (NI-1)*3+1: (NI-1)*3+3)
               CALL KARDIR(1,VTMP,LATT_CUR%B)
               WRITE(IO%IU6, '(3F13.8,6X,3F13.8)') VTMP+T_INFO%POSION(1:3,NI), VTMP
               IF (IBRION==-1) THEN
                  POSIOC(1:3,NI)= VTMP+T_INFO%POSION(1:3,NI)
               ENDIF
            ENDDO
            WRITE(IO%IU6,"( ' ----------------------------------------------', &
     &          '-------------------------------------')")
         ENDIF
         DEALLOCATE(HESSIAN)
      ENDIF
      ENDIF
   END SUBROUTINE


!************************ SUBROUTINE MUL_INV_SECOND_DERIV *************
!
! this subroutine calculates the eigenvalues and eigenvectors
! of the supplied Hessian and then calculates
!   CUNI^-1 F
! For each eigenmode the step size is however restricted
! to be a maximum displacement 
!
!***********************************************************************

  SUBROUTINE MUL_INV_SECOND_DERIV(CUNI, DOF , FORCE, STEP)
    INTEGER DOF
    REAL(q) CUNI(:,:)
    REAL(q) FORCE(:)
    REAL(q) STEP(:)
    ! local
    REAL(q) :: HFEIG(DOF),W(3*DOF)
    INTEGER :: N1, N2, IFAIL, NDIM
    
    REAL(q), ALLOCATABLE ::  CTMP(:,:), CEIDB(:,:)
    
    NDIM = SIZE(CUNI,1)
    ALLOCATE(CTMP(DOF,DOF),CEIDB(DOF,DOF))
!=======================================================================
! diagononalize the Hessian
!=======================================================================
    CALL DSYEV &
            ('V','U',DOF,CUNI(1,1),NDIM, &
            HFEIG,CTMP,DOF*DOF, IFAIL)

    IF (IFAIL/=0) THEN
       CALL vtutor%bug("internal ERROR in MUL_INV_SECOND_DERIV Call to routine DSYEV failed! Error " &
          // "code was " // str(IFAIL), __FILE__, __LINE__)
    ENDIF
!=======================================================================
! set up the unitary transformation CEIDB divided by eigenvalues
! skip lowest three eigenvalues which should correspond to 
! translational modes
! also skip soft modes close to zero
!=======================================================================
    CTMP=0
    ! first three eigenvalues are usually translation so skip them
    DO N2=4,DOF
       ! TODO: maybe some more tempering required
       ! the current setup damps all eigenvalues smaller than 1 by replacing them by 1
       CTMP(:,N2)=CUNI(:,N2)/MAX(HFEIG(N2),1.0_q)
    ENDDO

    CALL DGEMM( 'N',trans, DOF, DOF, DOF, one, CTMP, &
         &             DOF, CUNI(1,1), NDIM, zero, CEIDB, DOF)
   
    CUNI=0
    CUNI(1:DOF,1:DOF)=CEIDB(1:DOF,1:DOF)

    DEALLOCATE(CTMP, CEIDB)

    STEP=MATMUL(CUNI,FORCE)
    
  END SUBROUTINE MUL_INV_SECOND_DERIV

!****************** subroutine FORNLR2   *******************************
!
! this routine is not yet fully implemented; FORNL2 works though
!
! this subroutine is similar to FORNLR but allows for
! a different set of orbitals to the bra and ket states
! check comments above and SOFT_CHARGE1_RPA first
!
! the subroutine calculates the forces related to the non local
! pseudopotential in real space (factor 1/2 omitted here)
!  sum_n sum_ij  < n| d/dR p_i> D_ij - e_n Q_ij <p_j | n^RPA > + c.c.
!  sum_n sum_ij  < n| p_i> D_ij - e_n Q_ij <d/dR p_j | n^RPA > + c.c.
! approximated by finite difference
!  sum_n sum_ij  < n| p_i(+) > D_ij - e_n Q_ij  <p_j | n^RPA > + c.c.
! -sum_n sum_ij  < n| p_i(-) > D_ij - e_n Q_ij  <p_j | n^RPA > + c.c.
!
!  sum_n sum_ij  < n| p_i> D_ij - e_n Q_ij <p_j(+) | n > + c.c.
! -sum_n sum_ij  < n| p_i> D_ij - e_n Q_ij <p_j(-) | n > + c.c.
!
! orbitals n     are stored in W0  unperturbed Kohn-Sham orbitals
! orbitals n^RPA are stored in W1  first order change of each orbital
!
!***********************************************************************

  SUBROUTINE FORNLR2(GRID, NONLR_S, P, LATT_CUR, W0, W1, &
       &    CDIJ, CQIJ, DISPL0, FRNL)
    USE pseudo
    USE mpimy
    USE lattice
    USE constant
    USE nonlr

    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    TYPE (potcar)      P(NONLR_S%NTYP)
    TYPE (wavespin), TARGET ::    W1, W0
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR

    REAL(q) FRNL(3,NONLR_S%NIONS)
    REAL(q) DISPL0(3,NONLR_S%NIONS)
    OVERLAP CDIJ(:,:,:,:)
    OVERLAP CQIJ(:,:,:,:)
  ! local
    TYPE (wavespin)    WTMP1, WTMP2
    TYPE (wavedes), POINTER :: WDES
    REAL(q) ENL(NONLR_S%NIONS), EVALUE, WEIGHT
    COMPLEX(q) CE
    INTEGER IDIR, NK, ISP, N, ISPINOR, ISPINOR_, LBASE, LBASE_, NIS, NT, &
         LMMAXC, NI, L, LP, NIP
    REAL(q) DISPL1(3,NONLR_S%NIONS),DISPL2(3,NONLR_S%NIONS)
    ! allocate required work space
    GDEF,ALLOCATABLE,TARGET :: CPROW1(:,:,:,:), CPROW0(:,:,:,:)
  ! magnitude used for finite differencesq
    REAL(q) :: DIS=fd_displacement

    IF (.NOT. LINIT_RPA_FORCE) RETURN

    WDES=>W1%WDES

    ALLOCATE(CPROW1(WDES%NPROD,WDES%NBANDS,WDES%NKPTS,WDES%ISPIN), &
             CPROW0(WDES%NPROD,WDES%NBANDS,WDES%NKPTS,WDES%ISPIN))
!=======================================================================
! initialise non-local forces to zero
!=======================================================================
    FRNL=0
!=======================================================================
! calculate the contribution to the force from the nonlocal
! projection functions for displacement X using central (semianalytical)
! finite differences (about 9 digits precision)
!=======================================================================
    dir: DO IDIR=1,3
       ENL=0

       DISPL1=DISPL0
       DISPL1(IDIR,:)= DISPL0(IDIR,:)-DIS

       DISPL2=DISPL0
       DISPL2(IDIR,:)= DISPL0(IDIR,:)+DIS
       CALL RSPHER_ALL(GRID,NONLR_S,LATT_CUR,LATT_CUR,LATT_CUR, DISPL1, DISPL2, 1)

       WTMP1=W1
       WTMP1%CPROJ => CPROW1       ! relink the CPROJ array to temporary workspace

       WTMP2=W0
       WTMP2%CPROJ => CPROW0       ! relink the CPROJ array to temporary workspace

       kpoint: DO NK=1,WDES%NKPTS
#ifdef MPI
          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
          CALL PHASER(GRID,LATT_CUR,NONLR_S, NK,WDES)
          CALL RPRO(NONLR_S,WDES,WTMP1,GRID,NK)  ! this sets the content of CPROW1 to <p_k(+)-p_k(-) |n^RPAi>
          CALL RPRO(NONLR_S,WDES,WTMP2,GRID,NK)  ! this sets the content of CPROW0 to <p_k(+)-p_k(-) |n>

          spin: DO ISP=1,WDES%ISPIN

             band: DO N=1,WDES%NBANDS
                EVALUE=W0%CELEN(N,NK,ISP)
                WEIGHT=WDES%WTKPT(NK)*WDES%RSPIN*0.5

                spinor: DO ISPINOR=0,WDES%NRSPINORS-1
                   DO ISPINOR_=0,WDES%NRSPINORS-1

                      LBASE =ISPINOR *WDES%NPRO/2
                      LBASE_=ISPINOR_*WDES%NPRO/2

                      NIS=1
                      typ: DO NT=1,WDES%NTYP
                         LMMAXC=WDES%LMMAX(NT)
                         IF (LMMAXC==0) GOTO 510

                         DO NI=NIS,WDES%NITYP(NT)+NIS-1
                            CE=0
                            ! contribution <n^RPA| p_k> D(0)_kl - e(0) Q(0)_kl <p_k(+) | n >
                            CALL ECCP_NL_(SIZE(CDIJ,1),LMMAXC,CDIJ(1,1,NI,ISP+ISPINOR_+2*ISPINOR),CQIJ(1,1,NI,ISP+ISPINOR_+2*ISPINOR), &
                                 EVALUE,W1%CPROJ(LBASE_+1,N,NK,ISP),CPROW0(LBASE+1,N,NK,ISP),CE)
                            CE=0
                            ENL(NI)=ENL(NI)+CE*WEIGHT
                            
                            ! contribution <n^RPA| p_k(+) > D(0)_kl - e(0) Q(0)_kl  <p_k | n >
                            CALL ECCP_NL_(SIZE(CDIJ,1),LMMAXC,CDIJ(1,1,NI,ISP+ISPINOR_+2*ISPINOR),CQIJ(1,1,NI,ISP+ISPINOR_+2*ISPINOR), &
                                 EVALUE,CPROW1(LBASE+1,N,NK,ISP),W0%CPROJ(LBASE_+1,N,NK,ISP),CE)
                            ENL(NI)=ENL(NI)+CE*WEIGHT

                            LBASE = LMMAXC+LBASE
                            LBASE_= LMMAXC+LBASE_
                         ENDDO
510                      NIS = NIS+WDES%NITYP(NT)
                      ENDDO typ

                   ENDDO
                ENDDO spinor
             ENDDO band
          ENDDO spin
       ENDDO kpoint

       DO NI=1,WDES%NIONS
          NIP=NI_GLOBAL(NI, WDES%COMM_INB)
          FRNL(IDIR,NIP)=FRNL(IDIR,NIP)-ENL(NI)/DIS
       ENDDO

    ENDDO dir

    CALLMPI( M_sum_d(WDES%COMM, FRNL(1,1),NONLR_S%NIONS*3))

    DEALLOCATE(CPROW1, CPROW0)

    RETURN
  END SUBROUTINE FORNLR2


!************************ SUBROUTINE FORNL2  ****************************
!
! this subroutine is similar to FORNL but allows for
! a different set of orbitals for the bra and ket states
! check comments above and SOFT_CHARGE1_RPA first
!
! the subroutine calculates the forces related to the non local
! pseudopotential in reciprocal space (factor 1/2 omitted here)
! =  sum_n sum_ij <n | d/dR p_i > D_ij - e_n Q_ij < p_j| n^RPA> + c.c.
!    sum_n sum_ij <n | p_i > D_ij - e_n Q_ij  < d/dR p_j| n^RPA> + c.c.
!  + sum_n sum_ij <n | d/dR p_i > Q_ij < p_j| n^S>  + c.c.
!    sum_n sum_ij <n | p_i > Q_ij < d/dR p_j | n^S> + c.c.
!
! orbitals n     are stored in W0  unperturbed Kohn-Sham orbitals
! orbitals n^RPA are stored in W1  first order change of each orbital
!
! Algorithm:
! F(R)  = SUM(L,L')
!        (SUM(G)  QPROJ1(G ,L)  EXP( iG  R) C(G)  iG)
!          (D(L,L') + Evalue  Q(L,L'))
!        (SUM(G)  QPROJ2(G,L') EXP(iG R) C(G))^* + c.c.
!
!       +(SUM(G)  QPROJ1(G ,L)  EXP( iG  R) C(G))
!          (D(L,L') + Evalue  Q(L,L'))
!        (SUM(G)  QPROJ2(G,L') EXP( iG R) C(G) iG)^* + c.c.
!
!***********************************************************************

  SUBROUTINE FORNL2(NONL_S,WDES, W0,W1,WS, LATT_CUR, LMDIM,CDIJ,CQIJ,EINL, EINLGAMMA )
    USE nonl
    USE lattice
    USE constant

    IMPLICIT NONE

    TYPE (nonl_struct) NONL_S
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W1, W0
    TYPE (wavespin), POINTER :: WS
    TYPE (latt)        LATT_CUR
    TYPE (wavedes1)    WDES1         ! descriptor for one k-point
    INTEGER :: LMDIM
    REAL(q) :: EINL(3,NONL_S%NIONS)
    REAL(q) :: EINLGAMMA(3,NONL_S%NIONS)
    OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    OVERLAP CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    !local
    INTEGER :: ISP, ISPINOR, ISPINOR_, ISPIRAL, K, KK, LM, LMBASE, LMBASE_, LMMAXC, &
               LMP, N, NI, NIP, NIS, NK, NL,  NPL, NT
    REAL(q) :: WEIGHT
    COMPLEX(q) :: CMUL, CVAL, EVALUE
    
    ! work arrays
    REAL(q) ::  ENL(NONL_S%NIONS)
    GDEF :: CM
    GDEF,ALLOCATABLE :: CX1(:) ,CY1(:) ,CZ1(:), C1(:)
    GDEF,ALLOCATABLE :: CXL1(:),CYL1(:),CZL1(:),CL1(:)
    GDEF,ALLOCATABLE :: CX2(:) ,CY2(:) ,CZ2(:), C2(:)
    GDEF,ALLOCATABLE :: CXL2(:),CYL2(:),CZL2(:),CL2(:)

    LOGICAL :: NEW_VERSION=.TRUE.

    IF (.NOT. LINIT_RPA_FORCE) RETURN

    N =WDES%NPRO_TOT
    NL=WDES%NPRO

    ALLOCATE(CX1(N),CY1(N),CZ1(N),C1(N),CXL1(NL),CYL1(NL),CZL1(NL),CL1(NL))
    ALLOCATE(CX2(N),CY2(N),CZ2(N),C2(N),CXL2(NL),CYL2(NL),CZL2(NL),CL2(NL))

    EINL=0._q
    EINLGAMMA=0._q
    ENL =0._q
!=======================================================================
! loop over special points, and bands
!=======================================================================

    kpoint: DO NK=1,WDES%NKPTS
#ifdef MPI
       IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
       NPL= WDES%NGVECTOR(NK)
       CALL SETWDES(WDES,WDES1,NK)

       CALL PHASE(WDES,NONL_S,NK)

       spin:   DO ISP=1,WDES%ISPIN
          band: DO N=1,WDES%NBANDS

             EVALUE=W0%CELEN(N,NK,ISP)
             IF (ASSOCIATED(WS) .AND. NEW_VERSION) EVALUE=0
             WEIGHT=WDES%WTKPT(NK)*WDES%RSPIN*0.5
!=======================================================================
! first build up CX, CY, CZ for table
! this loop uses W0 and W1 (=W_RPA)
!=======================================================================
             CX1=0
             CY1=0
             CZ1=0
             C1=0

             CX2=0
             CY2=0
             CZ2=0
             C2=0

             LMBASE= 0

             ISPIRAL = 1
             spinor: DO ISPINOR=0,WDES1%NRSPINORS-1

                NIS=1
                typ:  DO NT=1,NONL_S%NTYP
                   LMMAXC=NONL_S%LMMAX(NT)
                   IF (LMMAXC==0) GOTO 100

                   ions: DO NI=NIS,NONL_S%NITYP(NT)+NIS-1
                      l_loop: DO LM=1,LMMAXC
                         CMUL=NONL_S%CQFAK(LM,NT)*CITPI
                         !DIR$ IVDEP
                         !OCL NOVREC

                         DO K=1,NPL                 
                            KK=K+NPL*ISPINOR
                            CVAL =  (NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*NONL_S%CREXP(K,NI)*W0%CPTWFP(KK,N,NK,ISP))*CMUL
                            ! new version calculate the projected wave function character on the fly
                            C2(LM+LMBASE) =C2(LM+LMBASE) +(NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*NONL_S%CREXP(K,NI)*W0%CPTWFP(KK,N,NK,ISP))*NONL_S%CQFAK(LM,NT)
#ifdef shift_der_k
                            ! new version add i k which contributes nothing to the energy derivative
                            CX2(LM+LMBASE)=CX2(LM+LMBASE)+(WDES%IGX(K,NK)+WDES%VKPT(1,NK))*CVAL
                            CY2(LM+LMBASE)=CY2(LM+LMBASE)+(WDES%IGY(K,NK)+WDES%VKPT(2,NK))*CVAL
                            CZ2(LM+LMBASE)=CZ2(LM+LMBASE)+(WDES%IGZ(K,NK)+WDES%VKPT(3,NK))*CVAL
#else
                            CX2(LM+LMBASE)=CX2(LM+LMBASE)+(WDES%IGX(K,NK))*CVAL
                            CY2(LM+LMBASE)=CY2(LM+LMBASE)+(WDES%IGY(K,NK))*CVAL
                            CZ2(LM+LMBASE)=CZ2(LM+LMBASE)+(WDES%IGZ(K,NK))*CVAL
#endif

                            CVAL =  (NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*NONL_S%CREXP(K,NI)*W1%CPTWFP(KK,N,NK,ISP))*CMUL
                            ! new version calculate the projected wave function character on the fly
                            C1(LM+LMBASE) =C1(LM+LMBASE) +(NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*NONL_S%CREXP(K,NI)*W1%CPTWFP(KK,N,NK,ISP))*NONL_S%CQFAK(LM,NT)
#ifdef shift_der_k
                            ! new version add i k which contributes nothing to the energy derivative
                            CX1(LM+LMBASE)=CX1(LM+LMBASE)+(WDES%IGX(K,NK)+WDES%VKPT(1,NK))*CVAL
                            CY1(LM+LMBASE)=CY1(LM+LMBASE)+(WDES%IGY(K,NK)+WDES%VKPT(2,NK))*CVAL
                            CZ1(LM+LMBASE)=CZ1(LM+LMBASE)+(WDES%IGZ(K,NK)+WDES%VKPT(3,NK))*CVAL
#else
                            CX1(LM+LMBASE)=CX1(LM+LMBASE)+(WDES%IGX(K,NK))*CVAL
                            CY1(LM+LMBASE)=CY1(LM+LMBASE)+(WDES%IGY(K,NK))*CVAL
                            CZ1(LM+LMBASE)=CZ1(LM+LMBASE)+(WDES%IGZ(K,NK))*CVAL
#endif

                         ENDDO

                      ENDDO l_loop
                      LMBASE= LMMAXC+LMBASE
                   ENDDO ions

100                NIS = NIS+NONL_S%NITYP(NT)
                ENDDO typ
                IF (NONL_S%LSPIRAL) ISPIRAL=2
             ENDDO spinor

             ! distribute projectors among nodes
             CALL DIS_PROJ(WDES1,C1,CL1)
             CALL DIS_PROJ(WDES1,CX1,CXL1)
             CALL DIS_PROJ(WDES1,CY1,CYL1)
             CALL DIS_PROJ(WDES1,CZ1,CZL1)

             CALL DIS_PROJ(WDES1,C2,CL2)
             CALL DIS_PROJ(WDES1,CX2,CXL2)
             CALL DIS_PROJ(WDES1,CY2,CYL2)
             CALL DIS_PROJ(WDES1,CZ2,CZL2)
!-----------------------------------------------------------------------
! sum up local contributions
! calculate SUM_LP  D(LP,L)-E Q(LP,L) * C(LP)
!-----------------------------------------------------------------------
             spinor2 : DO ISPINOR=0,WDES%NRSPINORS-1
                DO ISPINOR_=0,WDES%NRSPINORS-1

                   LMBASE =ISPINOR *WDES%NPRO/2
                   LMBASE_=ISPINOR_*WDES%NPRO/2

                   NIS   =1
                   typ2:  DO NT=1,WDES%NTYP
                      LMMAXC=WDES%LMMAX(NT)
                      IF (LMMAXC==0) GOTO 600

                      ions2: DO NI=NIS,WDES%NITYP(NT)+NIS-1
                         NIP=NI_GLOBAL(NI, WDES%COMM_INB) !  local storage index
                         l_loop2: DO LM=1,LMMAXC
                            CM=0
                            DO LMP=1,LMMAXC
                               CM=     (CDIJ(LMP,LM,NI,ISP+ISPINOR_+2*ISPINOR)-EVALUE*CQIJ(LMP,LM,NI,ISP+ISPINOR_+2*ISPINOR))* &
                                    GCONJG(CL1(LM+LMBASE))
                               ENL(NIP)=ENL(NIP) + (W0%CPROJ(LMP+LMBASE_,N,NK,ISP)*CM)*WEIGHT
                               EINL(1,NIP)=EINL(1,NIP)-(2*WEIGHT)*(CXL2(LMP+LMBASE_)*CM)
                               EINL(2,NIP)=EINL(2,NIP)-(2*WEIGHT)*(CYL2(LMP+LMBASE_)*CM)
                               EINL(3,NIP)=EINL(3,NIP)-(2*WEIGHT)*(CZL2(LMP+LMBASE_)*CM)

                               CM=     (CDIJ(LMP,LM,NI,ISP+ISPINOR_+2*ISPINOR)-EVALUE*CQIJ(LMP,LM,NI,ISP+ISPINOR_+2*ISPINOR))* &
                                    GCONJG(CL2(LM+LMBASE))
                               ENL(NIP)=ENL(NIP) + (W1%CPROJ(LMP+LMBASE_,N,NK,ISP)*CM)*WEIGHT
                               EINL(1,NIP)=EINL(1,NIP)-(2*WEIGHT)*(CXL1(LMP+LMBASE_)*CM)
                               EINL(2,NIP)=EINL(2,NIP)-(2*WEIGHT)*(CYL1(LMP+LMBASE_)*CM)
                               EINL(3,NIP)=EINL(3,NIP)-(2*WEIGHT)*(CZL1(LMP+LMBASE_)*CM)

                            ENDDO
                         ENDDO l_loop2
                         LMBASE = LMMAXC+LMBASE
                         LMBASE_= LMMAXC+LMBASE_
                      ENDDO ions2
600                   NIS = NIS+WDES%NITYP(NT)
                   ENDDO typ2
                ENDDO
             ENDDO spinor2
!=======================================================================

!=======================================================================
! now the same using WS and Q (if WS is associated)
! this loop uses W0 and WS (=W_S)
! this contribution comes from the change of the orbitals due to changes
! in the overlap operator  
!                             |\delta \phi>
!                              ----------   \delta S
!                               \delta S
!=======================================================================
             IF (ASSOCIATED(WS) .AND. NEW_VERSION) THEN

             CX1=0
             CY1=0
             CZ1=0
             C1=0

             CX2=0
             CY2=0
             CZ2=0
             C2=0

             LMBASE= 0

             ISPIRAL = 1
             spinorS: DO ISPINOR=0,WDES1%NRSPINORS-1

                NIS=1
                typS:  DO NT=1,NONL_S%NTYP
                   LMMAXC=NONL_S%LMMAX(NT)
                   IF (LMMAXC==0) GOTO 110

                   ionsS: DO NI=NIS,NONL_S%NITYP(NT)+NIS-1
                      l_loopS: DO LM=1,LMMAXC
                         CMUL=NONL_S%CQFAK(LM,NT)*CITPI
                         !DIR$ IVDEP
                         !OCL NOVREC

                         DO K=1,NPL                 
                            KK=K+NPL*ISPINOR
                            CVAL =  (NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*NONL_S%CREXP(K,NI)*W0%CPTWFP(KK,N,NK,ISP))*CMUL
                            ! new version calculate the projected wave function character on the fly
                            C2(LM+LMBASE) =C2(LM+LMBASE) +(NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*NONL_S%CREXP(K,NI)*W0%CPTWFP(KK,N,NK,ISP))*NONL_S%CQFAK(LM,NT)
#ifdef shift_der_k
                            ! new version add i k which contributes nothing to the energy derivative
                            CX2(LM+LMBASE)=CX2(LM+LMBASE)+(WDES%IGX(K,NK)+WDES%VKPT(1,NK))*CVAL
                            CY2(LM+LMBASE)=CY2(LM+LMBASE)+(WDES%IGY(K,NK)+WDES%VKPT(2,NK))*CVAL
                            CZ2(LM+LMBASE)=CZ2(LM+LMBASE)+(WDES%IGZ(K,NK)+WDES%VKPT(3,NK))*CVAL
#else
                            CX2(LM+LMBASE)=CX2(LM+LMBASE)+(WDES%IGX(K,NK))*CVAL
                            CY2(LM+LMBASE)=CY2(LM+LMBASE)+(WDES%IGY(K,NK))*CVAL
                            CZ2(LM+LMBASE)=CZ2(LM+LMBASE)+(WDES%IGZ(K,NK))*CVAL
#endif

                            CVAL =  (NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*NONL_S%CREXP(K,NI)*WS%CPTWFP(KK,N,NK,ISP))*CMUL
                            ! new version calculate the projected wave function character on the fly
                            C1(LM+LMBASE) =C1(LM+LMBASE) +(NONL_S%QPROJ(K,LM,NT,NK,ISPIRAL)*NONL_S%CREXP(K,NI)*WS%CPTWFP(KK,N,NK,ISP))*NONL_S%CQFAK(LM,NT)
#ifdef shift_der_k
                            ! new version add i k which contributes nothing to the energy derivative
                            CX1(LM+LMBASE)=CX1(LM+LMBASE)+(WDES%IGX(K,NK)+WDES%VKPT(1,NK))*CVAL
                            CY1(LM+LMBASE)=CY1(LM+LMBASE)+(WDES%IGY(K,NK)+WDES%VKPT(2,NK))*CVAL
                            CZ1(LM+LMBASE)=CZ1(LM+LMBASE)+(WDES%IGZ(K,NK)+WDES%VKPT(3,NK))*CVAL
#else
                            CX1(LM+LMBASE)=CX1(LM+LMBASE)+(WDES%IGX(K,NK))*CVAL
                            CY1(LM+LMBASE)=CY1(LM+LMBASE)+(WDES%IGY(K,NK))*CVAL
                            CZ1(LM+LMBASE)=CZ1(LM+LMBASE)+(WDES%IGZ(K,NK))*CVAL
#endif

                         ENDDO

                      ENDDO l_loopS
                      LMBASE= LMMAXC+LMBASE
                   ENDDO ionsS

110                NIS = NIS+NONL_S%NITYP(NT)
                ENDDO typS
                IF (NONL_S%LSPIRAL) ISPIRAL=2
             ENDDO spinorS

             ! distribute projectors among nodes
             CALL DIS_PROJ(WDES1,C1,CL1)
             CALL DIS_PROJ(WDES1,CX1,CXL1)
             CALL DIS_PROJ(WDES1,CY1,CYL1)
             CALL DIS_PROJ(WDES1,CZ1,CZL1)

             CALL DIS_PROJ(WDES1,C2,CL2)
             CALL DIS_PROJ(WDES1,CX2,CXL2)
             CALL DIS_PROJ(WDES1,CY2,CYL2)
             CALL DIS_PROJ(WDES1,CZ2,CZL2)
!-----------------------------------------------------------------------
! sum up local contributions
! calculate SUM_LP  Q(LP,L) * C(LP)
!-----------------------------------------------------------------------
             spinor2S : DO ISPINOR=0,WDES%NRSPINORS-1
                DO ISPINOR_=0,WDES%NRSPINORS-1

                   LMBASE =ISPINOR *WDES%NPRO/2
                   LMBASE_=ISPINOR_*WDES%NPRO/2

                   NIS   =1
                   typ2S:  DO NT=1,WDES%NTYP
                      LMMAXC=WDES%LMMAX(NT)
                      IF (LMMAXC==0) GOTO 610

                      ions2S: DO NI=NIS,WDES%NITYP(NT)+NIS-1
                         NIP=NI_GLOBAL(NI, WDES%COMM_INB) !  local storage index
                         l_loop2S: DO LM=1,LMMAXC
                            CM=0
                            DO LMP=1,LMMAXC
                               CM=     (-CQIJ(LMP,LM,NI,ISP+ISPINOR_+2*ISPINOR))* &
                                    GCONJG(CL1(LM+LMBASE))
                               ENL(NIP)=ENL(NIP) + (W0%CPROJ(LMP+LMBASE_,N,NK,ISP)*CM)*WEIGHT
                               EINLGAMMA(1,NIP)=EINLGAMMA(1,NIP)-(2*WEIGHT)*(CXL2(LMP+LMBASE_)*CM)
                               EINLGAMMA(2,NIP)=EINLGAMMA(2,NIP)-(2*WEIGHT)*(CYL2(LMP+LMBASE_)*CM)
                               EINLGAMMA(3,NIP)=EINLGAMMA(3,NIP)-(2*WEIGHT)*(CZL2(LMP+LMBASE_)*CM)

                               CM=     (-CQIJ(LMP,LM,NI,ISP+ISPINOR_+2*ISPINOR))* &
                                    GCONJG(CL2(LM+LMBASE))
                               ENL(NIP)=ENL(NIP) + (WS%CPROJ(LMP+LMBASE_,N,NK,ISP)*CM)*WEIGHT
                               EINLGAMMA(1,NIP)=EINLGAMMA(1,NIP)-(2*WEIGHT)*(CXL1(LMP+LMBASE_)*CM)
                               EINLGAMMA(2,NIP)=EINLGAMMA(2,NIP)-(2*WEIGHT)*(CYL1(LMP+LMBASE_)*CM)
                               EINLGAMMA(3,NIP)=EINLGAMMA(3,NIP)-(2*WEIGHT)*(CZL1(LMP+LMBASE_)*CM)

                            ENDDO
                         ENDDO l_loop2S
                         LMBASE = LMMAXC+LMBASE
                         LMBASE_= LMMAXC+LMBASE_
                      ENDDO ions2S
610                   NIS = NIS+WDES%NITYP(NT)
                   ENDDO typ2S
                ENDDO
             ENDDO spinor2S
!=======================================================================
             ENDIF
          ENDDO band
       ENDDO spin
    ENDDO kpoint
!=======================================================================
! forces are now in the reciprocal lattice transform it to
! cartesian coordinates
!=======================================================================
    CALLMPI( M_sum_d(WDES%COMM, EINL(1,1),NONL_S%NIONS*3))
    CALLMPI( M_sum_d(WDES%COMM, EINLGAMMA(1,1),NONL_S%NIONS*3))
    CALLMPI( M_sum_d(WDES%COMM, ENL (1)  ,NONL_S%NIONS))
    CALL  DIRKAR(NONL_S%NIONS,EINL,LATT_CUR%B)
    CALL  DIRKAR(NONL_S%NIONS,EINLGAMMA,LATT_CUR%B)

!    IF (WDES%COMM%NODE_ME==WDES%COMM%IONODE) THEN
!       WRITE(*,'(3F14.7)') EINL
!    ENDIF

    DEALLOCATE(CX1,CY1,CZ1,C1,CXL1,CYL1,CZL1,CL1)
    DEALLOCATE(CX2,CY2,CZ2,C2,CXL2,CYL2,CZL2,CL2)
    
    RETURN
  END SUBROUTINE FORNL2


!*********************** SOFT_CHARGE1_RPA *********************************
!
! subroutine SOFT_CHARGE1_RPA constructs the electronic charge density using
!
!   rho^1(r)   =  sum_n  <n | r | n^RPA>  = sum_n  <n_RPA | r | n>  = 
!              = 0.5  sum_n (n*(r) n^RPA(r) +  n(r) n^RPA*(r))
!
!  or following the usual VASP notation n -> phi0_nk ; n^RPA -> phi1_nk
!
!  rho(r) = sum_n,k (phi0_nk(r) phi1_nk*(r) + c.c. ) * 1/2
!
!***********************************************************************

      SUBROUTINE SOFT_CHARGE1_RPA(GRID,GRID_SOFT,W0,W1,WDES, CHDEN)
      USE prec
      USE mpimy
      USE mgrid
      USE wave
      USE wave_mpi
      IMPLICIT NONE 

      TYPE (wavedes)     WDES
      TYPE (wavespin)    W0
      TYPE (wavespin)    W1
      TYPE (grid_3d)     GRID,GRID_SOFT
      COMPLEX(q)   CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)

    ! work arrays
      RGRID,   ALLOCATABLE :: CDWORK(:,:)
      COMPLEX(q), ALLOCATABLE :: CPTDUM0(:), CPTDUM1(:)
      INTEGER ISPINOR, ISPINOR_, ISP, NK, N, NPL, M, MM, MM_, I
      REAL(q) :: WEIGHT

      ! MPLWV is the allocation in complex words
      ! hence if CDWORK is REAL one needs to double the allocation 
      ALLOCATE(CDWORK(GRID%MPLWV*2,WDES%NCDIJ), &
           CPTDUM0(GRID%MPLWV*WDES%NRSPINORS),CPTDUM1(GRID%MPLWV*WDES%NRSPINORS))

      CDWORK=0
!=======================================================================
! loop over k-points and bands
!=======================================================================
      spin:    DO ISP=1,WDES%ISPIN
      kpoints: DO NK=1,WDES%NKPTS
#ifdef MPI
               IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
      band:    DO N=1,WDES%NBANDS
         WEIGHT=WDES%RSPIN*WDES%WTKPT(NK)*0.5
         IF (WEIGHT==0) CYCLE

         NPL=WDES%NGVECTOR(NK)
!=======================================================================
! fourier-transformation of wave-function
! sum up  real space charge density
!=======================================================================
         DO ISPINOR=0,WDES%NRSPINORS-1
            CALL FFTWAV(NPL,WDES%NINDPW(1,NK),CPTDUM0(1+ISPINOR*GRID%MPLWV),W0%CPTWFP(1+ISPINOR*NPL,N,NK,ISP),GRID)
            CALL FFTWAV(NPL,WDES%NINDPW(1,NK),CPTDUM1(1+ISPINOR*GRID%MPLWV),W1%CPTWFP(1+ISPINOR*NPL,N,NK,ISP),GRID)
         ENDDO
         spinor: DO ISPINOR=0,WDES%NRSPINORS-1 
         DO ISPINOR_=0,WDES%NRSPINORS-1 
            DO M=1,GRID%RL%NP
              MM =M+ISPINOR *GRID%MPLWV
              MM_=M+ISPINOR_*GRID%MPLWV
              CDWORK(M,ISP+ISPINOR_+2*ISPINOR)=CDWORK(M,ISP+ISPINOR_+2*ISPINOR)+ &
                   CPTDUM0(MM)*CONJG(CPTDUM1(MM_))*WEIGHT+ &
                   CPTDUM1(MM)*CONJG(CPTDUM0(MM_))*WEIGHT
           ENDDO
        ENDDO
        ENDDO spinor
      ENDDO band
      ENDDO kpoints
      ENDDO spin

!=======================================================================
! fourier-transformation of charge-density using GRID_SOFT (see above)
!=======================================================================
      DO I=1,WDES%NCDIJ
! now merge the chargedensity from all nodes

!PK Reduce onto first KINTER nodes then broadcast
#ifdef realmode
         CALLMPI( M_sum_d(WDES%COMM_INTER, CDWORK(1,I), GRID%RL%NP))
         CALLMPI( M_sum_d(WDES%COMM_KINTER, CDWORK(1,I), GRID%RL%NP))
#else
         CALLMPI( M_sum_z(WDES%COMM_INTER, CDWORK(1,I), GRID%RL%NP))
         CALLMPI( M_sum_z(WDES%COMM_KINTER, CDWORK(1,I), GRID%RL%NP))
#endif
         CALL FFT_RC_SCALE(CDWORK(1,I),CHDEN(1,I),GRID_SOFT)
! set the charge-density of unbalanced lattic-vectors to 0
         CALL SETUNB(CHDEN(1,I),GRID_SOFT)
      ENDDO

      DEALLOCATE(CDWORK, CPTDUM0, CPTDUM1)

      RETURN
      END SUBROUTINE


!************************ SUBROUTINE DEPSUM1_RPA ****************************
!
! this subroutine calculates  the first order change of the 
! one center (on site) "occupancy" 
! matrix of each ll'mm'augmentation channel, 
! the weight of each k-point WTKPT and  the wavefunction character 
! of all bands
! result is stored in CRHODE1
!
!***********************************************************************

      SUBROUTINE DEPSUM1_RPA(W0, W1, WDES, LMDIM, CRHODE, LOVERL)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavespin)    W0
      TYPE (wavespin)    W1
      TYPE (wavedes)     WDES
      LOGICAL LOVERL
      INTEGER LMDIM
      OVERLAP CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)

      INTEGER ISP, NT, NK, N, ISPINOR, ISPINOR_, LMBASE, LMBASE_, NIS, &
           LMMAXC, NI, L, LP
      REAL(q) WEIGHT

      IF (.NOT.LOVERL) RETURN
!=======================================================================
! initialise to zero
!=======================================================================
      CRHODE=0
!=======================================================================
! loop over all bands and k-points
!=======================================================================
      spin:   DO ISP=1,WDES%ISPIN
      kpoint: DO NK=1,WDES%NKPTS
#ifdef MPI
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
      band:   DO N=1,WDES%NBANDS
      WEIGHT=WDES%RSPIN*WDES%WTKPT(NK)*0.5

      spinor: DO ISPINOR =0,WDES%NRSPINORS-1
      DO ISPINOR_=0,WDES%NRSPINORS-1

      LMBASE =ISPINOR *(WDES%NPRO/2)
      LMBASE_=ISPINOR_*(WDES%NPRO/2)

      NIS   =1
      typ:  DO NT=1,WDES%NTYP
      LMMAXC=WDES%LMMAX(NT)
      IF (LMMAXC==0) GOTO 210

      ion: DO NI=NIS,WDES%NITYP(NT)+NIS-1

!DIR$ IVDEP
!OCL NOVREC
        DO L=1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
        DO LP=1,LMMAXC
           CRHODE(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)=CRHODE(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)+ &
                WEIGHT*W1%CPROJ(L+LMBASE,N,NK,ISP)*GCONJG(W0%CPROJ(LP+LMBASE_,N,NK,ISP))+ &
                WEIGHT*W0%CPROJ(L+LMBASE,N,NK,ISP)*GCONJG(W1%CPROJ(LP+LMBASE_,N,NK,ISP))
        ENDDO
        ENDDO
      
      LMBASE = LMMAXC+LMBASE
      LMBASE_= LMMAXC+LMBASE_

      ENDDO ion

  210 NIS = NIS+WDES%NITYP(NT)
      ENDDO typ
      ENDDO
      ENDDO spinor

      ENDDO band
      ENDDO kpoint
      ENDDO spin
      ! sum over all bands
#ifdef realmode
      CALLMPI( M_sum_d(WDES%COMM_INTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
      CALLMPI( M_sum_d(WDES%COMM_KINTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ))
#else
      CALLMPI( M_sum_d(WDES%COMM_INTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ*2))
      CALLMPI( M_sum_d(WDES%COMM_KINTER,CRHODE,LMDIM*LMDIM*WDES%NIONS*WDES%NCDIJ*2))
#endif


      RETURN
      END SUBROUTINE

#endif
END MODULE rpa_force
