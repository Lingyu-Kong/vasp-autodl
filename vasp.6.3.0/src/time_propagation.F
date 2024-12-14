#include "symbol.inc"

!*********************************************************************
!
! This module implements the time-evolution code
! using a finite time propagation
! presently implemented for
!  ) DFT
!  ) HSE
! as opposed to BSE, the TDA approximation is never applied
! and the time evolution is performed without ever storing
! the BSE matrix explicitly
!
! all the parallel features of VASP are supported
! NOTE:
!  ) currently unoccupied states are required for the groundstate
!    since the orbitals are orthogonalized to the groundstate orbitals
!  ) to improved the conditioning unoccupied states
!    beyond the available occupied states are projected out from
!    the solution
!
!*********************************************************************

! unit conversions in the module are somewhat akward
! they essentially have been inherited from the time evolution BSE code
! written by Bechstedt's group
!
! ) the time step is in units of 1/eV
!   this sounds akward is however quite convenient and can be conceived as
!   dividing the real time by hbar
!     i hbar  d/dt = i   d/d (t/hbar) = i d/d (t') with t'=t/hbar
!
! ) Here is a collection of the other factors used in the routine
! routine CALC_WRITE_DF:
!      CELLVOL = LATT_CUR%OMEGA/AUTOA**3
!      HARTREE  = 2._q*REAL(RYTOEV,q)
!    4._q*REAL(PI,q)*(HARTREE**3)/(CELLVOL)
! twice there are factors
!    AUTOA/(2._q*HSQDTM) = 1 /(2* RYTOEV* AUTOA)
! once when the delta kick is done, and once when the dipole is calculated (both
! use the stored matrix elements).

! there seems to be no special physical meaning to the individual conversions,
! however, their combination is simply
!   4 PI (2*RYTOEV*AUTOA)/ LATT_CUR%OMEGA
! This conversion factor is used everywhere in VASP to scale the product of response functions
! and the Coulomb kernels to unit less dielectric functions

MODULE time_propagation
  USE prec
  USE constant
  USE bse
  IMPLICIT NONE

!global variables
!--------------------------------------------------
 COMPLEX(q),ALLOCATABLE, PRIVATE :: SCALAR(:,:)
 REAL(q),PRIVATE, PARAMETER      :: HARTREE  = 2._q*REAL(Rytoev,q) !27.2112_q
 REAL(q),PRIVATE                 :: CELLVOL
 COMPLEX(q),PRIVATE, PARAMETER   :: CZERO    = (0.0_q,0.0_q)
 REAL(q),PRIVATE, ALLOCATABLE    :: TIMESTEP_V(:,:) ! ::timestep_v(2,ntdsteps) -> Time Steps And Cumulated Time Steps For Td
 COMPLEX(q), PARAMETER, PRIVATE  :: ICI=(0.0_q,-1.0_q), CI=(0.0_q,1.0_q)
!--------------------------------------------------

  CONTAINS
  SUBROUTINE CALCULATE_TIME_PROPAGATION ( &
          KINEDEN, HAMILTONIAN, P, WDES, NONLR_S, NONL_S, W, LATT_CUR, LATT_INI, &
          GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
          CHTOT, CHTOTL, DENCOR, CVTOT, CSTRF, IRDMAX, &
          T_INFO, DYN, INFO, IO, KPOINTS, SYMM, &
          LMDIM, CQIJ, CDIJ, CRHODE, N_MIX_PAW, RHOLM, CHDEN, SV, &
          EFERMI, NEDOS, DOS, DOSI , IEPSILON, &
          NBANDSO, NBANDSV, OMEGAMAX_IN, SHIFT, NELMGW_IN, &
          LHARTREE, LADDER, LTRIPLET, LFXC, LGWLF, ANTIRES)
    USE base
    USE lattice
    USE pseudo
    USE nonl_high
    USE msymmetry
    USE constant
    USE poscar
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
    USE hamil_struct_def
    USE david
    USE subrot
    USE meta
    USE choleski
    USE pead, ONLY : LPEAD_SYM_RED
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
    INTEGER  LMDIM
    OVERLAP  CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    OVERLAP  CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
!  paw sphere charge density
    INTEGER N_MIX_PAW
    REAL(q)  RHOLM(N_MIX_PAW,WDES%NCDIJ)
!  charge-density and potential on soft grid
    COMPLEX(q)  CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
    RGRID       SV(DIMREAL(WDES%GRID%MPLWV),WDES%NCDIJ)
    REAL(q)     EFERMI
    INTEGER     NEDOS
    REAL(q)     DOS(NEDOS,WDES%ISPIN),DOSI(NEDOS,WDES%ISPIN)
    INTEGER             :: IEPSILON      ! cartesian direction (4=all)
!   flags specific to time evoluation
    INTEGER              NBANDSO    ! number of occupied states included in time evolution
    INTEGER              NBANDSV    ! number of virtual states included in time evolution
    REAL(q)              OMEGAMAX_IN! maximum considered transition energy in BSE
    REAL(q)              OMEGEMAX   ! actually maximum considered transition energy
    REAL(q)              SHIFT      ! complex frequency shift for broadening
    INTEGER              NELMGW_IN  ! number of timesteps in solving the BSE equation using time propagation

    ! I am not sure that we will can support the flags below
    ! right now we simply apply the full Hamiltonian in the time evolution
    ! maybe by updating only selected terms in the Hamiltonian this could be done ???
    LOGICAL              LHARTREE   ! include Hartree effects TODO presently not supported
    LOGICAL              LADDER     ! local field effects beyond RPA (ladder diagrams) TODO presently not supported
    LOGICAL              LTRIPLET   ! calculate triplet TODO presently not supported
    LOGICAL              LGWLF      ! use screened W from GW for exchange
    LOGICAL              LFXC       ! include DFT xc kernel
    INTEGER              ANTIRES    ! how to treat antiresonant part TODO presently not supported
    REAL(q)              OMEGAMAX   ! actual maximum frequency


! symmetry related quantities (common block)
    INTEGER  ISYMOP,NROT,IGRPOP,NROTK,INVMAP,NPCELL
    INTEGER :: NKORIG
    INTEGER :: DIR_MIN, DIR_MAX
    REAL(q)  GTRANS,AP
    COMMON /SYMM/   ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
         GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL
    REAL(q)  DESUM
    TYPE(banddesc)       BD(WDES%ISPIN)! descriptor for describing which bands are included, spin up
    INTEGER              ISP, NELMGW,NELMGW_ACTUAL
    LOGICAL              UPDATE_CHARGE ! determines whether charge and potential is updated

!=======================================================================
    INTEGER IDIR
    INTEGER IU6K
    INTEGER NCV                        ! number of conduction valence pairs
    TYPE(bse_matrix_index)  :: BSE_INDEX
    TYPE (skpoints_trans)   :: KPOINTS_TRANS
!========================================

    IF (IO%NWRITE>=3) THEN
       IU6K  = IO%IU6
    ELSE
       IU6K  = -1
    ENDIF
    NKORIG=WDES%NKPTS
    SHIFT=ABS(SHIFT)

    IF (NBANDSO <= 0) THEN
       NBANDSO=INFO%NELECT/2+0.9
    ENDIF
    IF (NBANDSV <= 0) THEN
       NBANDSV=INFO%NELECT/2+0.9
    ENDIF

    NELMGW=20000
    IF (NELMGW_IN>100) THEN
       NELMGW=NELMGW_IN
    ENDIF


    CALL READ_CDER_BETWEEN_STATES(WDES, IO%IU0, 55)
!=======================================================================
! update charge and then potential
!=======================================================================
    CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
         GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
         LATT_CUR, P, SYMM, T_INFO, &
         CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

    CALL SET_KINEDEN(GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM, &
         T_INFO%NIONS,W,WDES,KINEDEN)

    CALL STOP_TIMING("G",IO%IU6,"CHARGE")

    CALL UPDATE_POTENTIAL(  &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,LATT_CUR, &
          T_INFO,INFO,IO, &
          GRID,GRID_SOFT,GRIDC,GRIDUS,C_TO_US,SOFT_TO_C, &
          CHTOT,DENCOR,CVTOT,CSTRF, &
          LMDIM,IRDMAX,CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,CHDEN,SV)
!=======================================================================
! dump out some of the important parameters and set those that are
! not defaulted
!=======================================================================

    ! determine which bands are included in the update procedure
    DO ISP=1,WDES%ISPIN
       CALL SET_BAND_PARAMETERS(W, ISP, NBANDSO, NBANDSV, &
            BD(ISP)%VBMAX, BD(ISP)%VBMIN, BD(ISP)%CBMIN, BD(ISP)%CBMAX)
    ENDDO

    IF (.NOT. LFXC .AND. .NOT. LHARTREE) THEN
       ! no update of potential at all this is IP approximation
       ! and that is easy to implement and supported
       UPDATE_CHARGE=.FALSE.
    ELSE
       UPDATE_CHARGE=.TRUE.
       ! all other settings, update all and everything
       ! this is simply to get the correct write statements below (copied from bse.F)
       ! the code otherwise currently neglects these flags
       LFXC  = .TRUE.
       LHARTREE = .TRUE.
    ENDIF

    IF (IO%IU6>=0) THEN
       DO ISP=1,WDES%ISPIN
          WRITE(IO%IU6,'(/A,/2(A,I5,2X,A,I5/))')' Bands included in the BSE', &
            ' VB(min)=',BD(ISP)%VBMIN,'VB(max)=',BD(ISP)%VBMAX,' CB(min)=',BD(ISP)%CBMIN,'CB(max)=',BD(ISP)%CBMAX
       END DO


       IF (OMEGAMAX_IN>0) THEN
          WRITE(IO%IU6,'(A,F10.5)') 'electron-hole pairs beyond OMEGAMAX are removed OMEGAMAX=',OMEGAMAX_IN
       ENDIF


       IF (LADDER .AND. (AEXX /=0 .OR. L_MODEL_HF .OR. LGWLF) ) THEN
          IF (LGWLF) THEN
             CALL vtutor%error("W is read from the files WXXXX.tmp \n this is presently not supported")
          ELSE
             WRITE(IO%IU6,'(A,F14.7,A,F14.7)')' parameters for screened Coulomb W: AEXX=',AEXX,' HFSCREEN=',HFSCREEN
          ENDIF
          IF (.NOT. LHARTREE .OR. LTRIPLET) THEN
             CALL vtutor%error("only ladders are included (Hartree term=RPA part is switched off) \n &
                &this is presently not supported")
          ENDIF
       ELSE
          IF (.NOT. LHARTREE .OR. LTRIPLET) THEN
             WRITE(IO%IU6,'(A)') ' IP particle spectrum'
          ELSE
             WRITE(IO%IU6,'(A)')' simple RPA calculation, excitonic effects (ladders) are not included'
          ENDIF
       ENDIF
       IF (LHARTREE .AND. LFXC .AND. .NOT. LGWLF) THEN
          WRITE(IO%IU6,'(A)')' xc-kernel from DFT is included (LFXC=.TRUE.)'
       ENDIF

    ENDIF

!=======================================================================
! response to external time-dependent field
!=======================================================================
    IF (IEPSILON==1) THEN
       DIR_MIN=1
       DIR_MAX=1
    ELSE IF (IEPSILON==2) THEN
       DIR_MIN=2
       DIR_MAX=2
    ELSEIF (IEPSILON==3) THEN
       DIR_MIN=3
       DIR_MAX=3
    ELSEIF (IEPSILON==4) THEN
       DIR_MIN=1
       DIR_MAX=3
    ENDIF
    DO IDIR=DIR_MIN, DIR_MAX
       
       IF (IO%IU0>=0) THEN
          WRITE (IO%IU0,*)'Linear response to external field (no local field effect), progress :'
          WRITE (IO%IU0,'(A,I3)') &
               '  Direction: ',IDIR
          WRITE (17,*)'Linear response to external field (no local field effect), progress :'
          WRITE (17,'(A,I3)') &
               '  Direction: ',IDIR
       END IF

       IF (IO%IU6>=0) THEN
          WRITE (IO%IU6,*)'Linear response to external field (no local field effect), progress :'
          WRITE (IO%IU6,'(A,I3)') 'Direction: ',IDIR
       ENDIF

       IF (SYMM%ISYM>0) THEN
          DYN%VEL=0
          DYN%VEL(IDIR,:)=1.0
          CALL  KARDIR (T_INFO%NIONS,DYN%VEL,LATT_CUR%B)
          CALL INISYM (LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
               T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
               SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
               SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)
#ifdef oldsym
       ! we can not use TI symmetry
          CALL RE_READ_KPOINTS (KPOINTS,LATT_CUR, &
               .FALSE.
               T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IU6K,IO%IU0)
#else
       ! we can not use TI symmetry
          CALL RE_READ_KPOINTS (KPOINTS,LATT_CUR, &
               .FALSE., &
               T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IU6K,IO%IU0)
#endif
          CALL KPAR_SYNC_ALL (WDES,W)
          CALL RE_GEN_LAYOUT ( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)

          CALL REALLOCATE_WAVE ( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR, KPOINTS_TRANS)
          ! gK: added on 14.03.2014
          ! seems to be required for e.g. ZnO
          ! Loewdin perturbation theory to improve states at added k-points
          ! not yet sure whether 3 (ext diag) or 4 (Loewdin) should be used
          CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
          &    LMDIM,CDIJ,CQIJ,3,SV,T_INFO,P,IO%IU0,DESUM, NKSTART=NKORIG+1)


       ENDIF

! just in case sync all data between nodes
       CALL KPAR_SYNC_ALL(WDES,W)

       OMEGAMAX=OMEGAMAX_IN
       CALL SET_BSE_MATRIX_INDEX (1, WDES%ISPIN, WDES%NKPTS, BD, BSE_INDEX, W, OMEGAMAX, NCV, 1)

! perform electron dynamics for a field in cartesian direction IDIR
       CALL PROPAGATE ( &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W, &
          LATT_CUR, &
          T_INFO,INFO,IO,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDUS,C_TO_US,SOFT_TO_C, &
          CHTOT,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM, &
          CHDEN,SV,DOS,DOSI, &
          LMDIM,IRDMAX,NEDOS, &
          EFERMI, KPOINTS_TRANS, IDIR, &
          BD, BSE_INDEX, OMEGAMAX, SHIFT, NELMGW,NELMGW_ACTUAL, UPDATE_CHARGE, LADDER)

       CALL DEALLOCATE_BSE_MATRIX_INDEX ( BSE_INDEX )

       IF (IO%LOPEN) CALL WFORCE(IO%IU6)
    ENDDO ! IDIR

    CALL CALC_WRITE_DF(IO,OMEGAMAX, MAX(NEDOS, 1000), LATT_CUR%OMEGA, shift, W%WDES%ISPIN, IEPSILON, NELMGW_ACTUAL)

  END SUBROUTINE CALCULATE_TIME_PROPAGATION


!*********************************************************************
!
! This module implements the time-evolution for a field in direction
! IDIR  using a finite time propagation
!
!*********************************************************************

  SUBROUTINE PROPAGATE ( &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W0, &
          LATT_CUR, &
          T_INFO,INFO,IO,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDUS,C_TO_US,SOFT_TO_C, &
          CHTOT,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM, &
          CHDEN,SV,DOS,DOSI, &
          LMDIM,IRDMAX,NEDOS, &
          EFERMI, KPOINTS_TRANS, IDIR, &
          BD, BSE_INDEX, OMEGAMAX, SHIFT, NELMGW, NELMGW_ACTUAL,UPDATE_CHARGE, UPDATE_EXX)


    USE base
    USE lattice
    USE charge
    USE pseudo
    USE nonl_high
    USE msymmetry
    USE mpimy
    USE mgrid
    USE mkpoints
    USE constant
    USE pot
    USE ini
    USE hamil_lrf
    USE subrot_lr
    USE us
    USE wave_high
    USE choleski
    USE kpoints_change
    USE hamil_struct_def
    USE meta
    IMPLICIT NONE
!=======================================================================
!  structures
!=======================================================================
    TYPE (tau_handle)  KINEDEN
    TYPE (ham_handle)  HAMILTONIAN
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    TYPE (wavedes)     WDES
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W0         ! unperturbed wavefunctions
    TYPE (latt)        LATT_CUR
    TYPE (info_struct) INFO
    TYPE (in_struct)   IO
    TYPE (kpoints_struct) KPOINTS
    TYPE (symmetry)    SYMM
    TYPE (grid_3d)     GRID       ! grid for wavefunctions
    TYPE (grid_3d)     GRID_SOFT  ! grid for soft chargedensity
    TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
    TYPE (grid_3d)     GRIDUS     ! temporary grid in us.F
    TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
    TYPE (transit)     SOFT_TO_C  ! index table between GRID_SOFT and GRIDC
    TYPE (energy)      E

    INTEGER :: LMDIM,IRDMAX,NEDOS
    REAL(q) :: EFERMI

    COMPLEX(q) :: CHTOT(GRIDC%MPLWV,WDES%NCDIJ) ! charge-density in real / reciprocal space
    RGRID       :: DENCOR(GRIDC%RL%NP)           ! partial core
    COMPLEX(q)  :: CVTOT(GRIDC%MPLWV,WDES%NCDIJ) ! local potential
    COMPLEX(q)  :: CSTRF(GRIDC%MPLWV,T_INFO%NTYP)

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
!  density of states
    REAL(q)    DOS(NEDOS,WDES%ISPIN),DOSI(NEDOS,WDES%ISPIN)
! local l-projected wavefunction characters (not really used here)
    REAL(q)    PAR(1,1,1,1,WDES%NCDIJ),DOSPAR(1,1,1,WDES%NCDIJ)
    TYPE (skpoints_trans)   :: KPOINTS_TRANS
    INTEGER :: IDIR                      ! direction

    TYPE(banddesc)       BD(WDES%ISPIN)! descriptor for describing which bands are included, spin up
    TYPE(bse_matrix_index) :: BSE_INDEX
    REAL(q)              OMEGAMAX   ! maximum considered transition energy in BSE
    REAL(q)              SHIFT      ! complex frequency shift for broadening
    INTEGER              NELMGW     ! number of timesteps in solving the BSE equation using time propagation
    LOGICAL              UPDATE_CHARGE ! determines whether charge and potential is updated
    LOGICAL              UPDATE_EXX ! update the exact exchange contribution

! local variables
    TYPE (wavedes1)    WDES1
    TYPE (wavespin)    WT0          ! perturbed orbitals  at timestep t_i
    TYPE (wavespin)    WHPHI        ! action of Hamiltonian on orbitals at timestep t_i+1
    COMPLEX(q), ALLOCATABLE ::  C_TRANS_PROP(:) ! transition propabilities
    COMPLEX(q), ALLOCATABLE ::  C_OCC_UNOCC1(:) ! coefficients determining the perturbed orbitals at time step t-1
    COMPLEX(q), ALLOCATABLE ::  C_OCC_UNOCC2(:) ! coefficients determining the perturbed orbitals at time step t
    COMPLEX(q), ALLOCATABLE ::  H_OCC_UNOCC(:)  ! action of Hamiltonian on perturbed orbitals at t
    COMPLEX(q), ALLOCATABLE ::  H0_OCC_UNOCC(:) ! action of unperturbed Hamiltonian on unperturbed orbitals
    COMPLEX(q), EXTERNAL    ::  ZDOTC

! perturbation of first order orbitals
    REAL(q), PARAMETER ::  DELTA=0.5E-2_q ! maybe even 1E-1_q
    INTEGER  :: NSIM, NELM,NELMGW_ACTUAL,J,I
    INTEGER  :: ISP, K1
!    COMPLEX(q), PARAMETER :: ICI=(0.0_q,-1.0_q), CI=(0.0_q,1.0_q)
    INTEGER :: INITIAL_STEPS=10
    REAL(q) :: TIME
    TYPE (in_struct)   NOIO

    ! noio perform allows to avoid timing information in each step
    ! simply gets too much
    NOIO=IO
    NOIO%IU6=-1

 130  FORMAT (5X, //, &
     &'----------------------------------------------------', &
     &'----------------------------------------------------'//)

! first redistribute the unperturbed orbitals over plane wave coefficients
     DO ISP=1,W0%WDES%ISPIN
        DO K1=1,W0%WDES%NKPTS
           IF (W0%WDES%WTKPT(K1)==0) CYCLE
#ifdef MPI
           IF (MOD(K1-1,W0%WDES%COMM_KINTER%NCPU).NE.W0%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
           CALL SETWDES(W0%WDES,WDES1,K1)
           CALL REDISTRIBUTE_PROJ( ELEMENTS( W0, WDES1, ISP))
           CALL REDISTRIBUTE_PW( ELEMENTS( W0, WDES1, ISP))
        ENDDO
     ENDDO
!=======================================================================
! allocate all required arrays
!=======================================================================

    IF (.NOT. ALLOCATED(TIMESTEP_V)) ALLOCATE (TIMESTEP_V(2,NELMGW))

    IF (IO%IU6>=0) WRITE(IO%IU6, 130)
    CALL CALC_TIMESTEPS_TE (SHIFT,IO, TIMESTEP_V, NELMGW, NELMGW_ACTUAL, OMEGAMAX)
    IF (IO%LOPEN) CALL WFORCE(IO%IU6) ! force a write of OUTCAR

    IF (.NOT. ALLOCATED(SCALAR)) THEN
       ALLOCATE(SCALAR(NELMGW_ACTUAL,3))
       SCALAR=0
    ENDIF

    TIME=0

    CALL ALLOCW(WDES,WT0)
    CALL ALLOCW(WDES,WHPHI)

    ALLOCATE(C_TRANS_PROP(BSE_INDEX%NCV), C_OCC_UNOCC1(BSE_INDEX%NCV), &
             C_OCC_UNOCC2(BSE_INDEX%NCV), H_OCC_UNOCC(BSE_INDEX%NCV), &
             H0_OCC_UNOCC(BSE_INDEX%NCV))

    C_TRANS_PROP=(0.0_q,0.0_q)
    NSIM=WDES%NSIM*2
#ifdef MPI
    ! MPI dividable by WDES%COMM_INTER%NCPU
    NSIM=((NSIM+WDES%COMM_INTER%NCPU-1)/WDES%COMM_INTER%NCPU)*WDES%COMM_INTER%NCPU
#endif
!=======================================================================
! act with ground state Hamiltonian onto groundstate orbitals
!=======================================================================
    CALL START_TIMING("LOOP")

    C_OCC_UNOCC1=(0.0_q,0.0_q)
    CALL DETERMINE_ORBITALS (W0, WT0, C_OCC_UNOCC1, BD, BSE_INDEX, DELTA)
    CALL STOP_TIMING("G",IO%IU6,"ORBITALS")

    ! orthogonalization of this perturbed vector
    ! most likely not require
!#define use_orthch !anc
#ifdef use_orthch
    CALL ORTHCH (WT0%WDES, WT0, WT0%WDES%LOVERL, LMDIM, CQIJ)
    CALL STOP_TIMING("G",IO%IU6,"ORTHCH")
#endif

    IF (UPDATE_CHARGE) THEN
    ! update charge and then potential
       CALL SET_CHARGE (WT0, WDES, INFO%LOVERL, &
         GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
         LATT_CUR, P, SYMM, T_INFO, &
         CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

       CALL SET_KINEDEN (GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM, &
         T_INFO%NIONS,WT0,WDES,KINEDEN)

       CALL STOP_TIMING("G",IO%IU6,"CHARGE")

       CALL UPDATE_POTENTIAL (  &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,LATT_CUR, &
          T_INFO,INFO,IO, &
          GRID,GRID_SOFT,GRIDC,GRIDUS,C_TO_US,SOFT_TO_C, &
          CHTOT,DENCOR,CVTOT,CSTRF, &
          LMDIM,IRDMAX,CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,CHDEN,SV)
    ENDIF
    ! calculate H - S epsilon |i> for all states
    CALL DETERMINE_HPHI_FAST (HAMILTONIAN, &
         GRID, KPOINTS, LATT_CUR, NONLR_S, NONL_S, T_INFO, P, WT0, W0, WDES, WHPHI, &
         LMDIM, WDES%NSIM, CQIJ,CDIJ, SV, SYMM, UPDATE_EXX)
    CALL STOP_TIMING("G",IO%IU6,"HAMIL")

    CALL PROJECT_ORBITALS_GROUNDSTATE (W0, WHPHI, H0_OCC_UNOCC,  BD, BSE_INDEX, DELTA)
    CALL STOP_TIMING("G",IO%IU6,"PROJECT")

    CALL SEPERATOR_TIMING(IO%IU6)
    CALL STOP_TIMING("LOOP",IO%IU6,XMLTAG='total')

    IF (IO%LOPEN) CALL WFORCE(IO%IU6) ! force a write of OUTCAR so that we get a timing on OUTCAR

!=======================================================================
! start be setting up the approximate perturbed orbitals
!=======================================================================
    ! set |i> = |i> + \sum_a c_ia |a>, a sum over unoccupied
    ! c_ia is the transition propability betwen <a| nabla |i>/(e_a-e_i)

!==========================================================================
!   C_TRANS_PROP(NCV13) - called for each direction IDIR seperately
!==========================================================================
    CALL INIT_OSCILLATOR (W0, C_TRANS_PROP, LATT_CUR, BD, BSE_INDEX, IDIR, .FALSE.)
!===============================================================
! store oscillator strength as start vector at time step n-1
! c_trans_prop contains all coefficients for one direction IDIR
!===============================================================
    C_OCC_UNOCC1=C_TRANS_PROP
!=======================================================================
! act with Hamiltonian on C_OCC_UNOCC1 and store action in H_OCC_UNOCC
!     v(2) = v(0) + dt/(i*hbar) H v(0)
!=======================================================================
    CALL START_TIMING("LOOP")

    CALL DETERMINE_ORBITALS (W0, WT0, C_OCC_UNOCC1, BD, BSE_INDEX, DELTA)
    CALL STOP_TIMING("G",IO%IU6,"ORBITALS")

    ! orthogonalization of this perturbed vector
    ! most likely not require
#ifdef use_orthch
    CALL ORTHCH (WT0%WDES, WT0, WT0%WDES%LOVERL, LMDIM, CQIJ)
    CALL STOP_TIMING("G",IO%IU6,"ORTHCH")
#endif

    IF (UPDATE_CHARGE) THEN
    ! update charge and then potential
       CALL SET_CHARGE (WT0, WDES, INFO%LOVERL, &
         GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
         LATT_CUR, P, SYMM, T_INFO, &
         CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

       CALL SET_KINEDEN (GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM, &
         T_INFO%NIONS,WT0,WDES,KINEDEN)

       CALL STOP_TIMING("G",IO%IU6,"CHARGE")

       CALL UPDATE_POTENTIAL (  &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,LATT_CUR, &
          T_INFO,INFO,IO, &
          GRID,GRID_SOFT,GRIDC,GRIDUS,C_TO_US,SOFT_TO_C, &
          CHTOT,DENCOR,CVTOT,CSTRF, &
          LMDIM,IRDMAX,CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,CHDEN,SV)
    ENDIF
    ! calculate H - S epsilon |i> for all states
    CALL DETERMINE_HPHI_FAST (HAMILTONIAN, &
         GRID, KPOINTS, LATT_CUR, NONLR_S, NONL_S, T_INFO, P, WT0, W0, WDES, WHPHI, &
         LMDIM, WDES%NSIM, CQIJ,CDIJ, SV, SYMM, UPDATE_EXX)
    CALL STOP_TIMING("G",IO%IU6,"HAMIL")

    CALL PROJECT_ORBITALS_GROUNDSTATE (W0, WHPHI, H_OCC_UNOCC,  BD, BSE_INDEX, DELTA)
    H_OCC_UNOCC=H_OCC_UNOCC-H0_OCC_UNOCC

    CALL STOP_TIMING("G",IO%IU6,"PROJECT")

    CALL SEPERATOR_TIMING(IO%IU6)
    CALL STOP_TIMING("LOOP",IO%IU6,XMLTAG='total')

    IF (IO%LOPEN) CALL WFORCE(IO%IU6) ! force a write of OUTCAR so that we get a timing on OUTCAR

!#define second_order_start
#ifdef second_order_start
!=======================================================================
! act with Hamiltonian on C_OCC_UNOCC1 and store action in H_OCC_UNOCC
!     v(2) = v(2) + [dt/(i*hbar)]^2 H H v(0)
!=======================================================================
    ! first order start equation
    ! vector(t_n) = vector(t_n-1) + dt i  H vector(t_n-1)
    C_OCC_UNOCC2=C_OCC_UNOCC1+ICI*timestep_v(1,1)*H_OCC_UNOCC

    CALL DETERMINE_ORBITALS (W0, WT0, H_OCC_UNOCC, BD, BSE_INDEX, DELTA)
    CALL STOP_TIMING("G",IO%IU6,"ORBITALS")

! orthogonalization of this perturbed vector
#ifdef use_orthch
    CALL ORTHCH (WT0%WDES, WT0, WT0%WDES%LOVERL, LMDIM, CQIJ)
    CALL STOP_TIMING("G",IO%IU6,"ORTHCH")
#endif

    IF (UPDATE_CHARGE) THEN
! update charge and then potential
       CALL SET_CHARGE( WT0, WDES, INFO%LOVERL, &
            GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
            LATT_CUR, P, SYMM, T_INFO, &
            CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

       CALL SET_KINEDEN(GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM, &
            T_INFO%NIONS,WT0,WDES,KINEDEN)

       CALL STOP_TIMING("G",IO%IU6,"CHARGE")

       CALL UPDATE_POTENTIAL (  &
            KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,LATT_CUR, &
            T_INFO,INFO,IO, &
            GRID,GRID_SOFT,GRIDC,GRIDUS,C_TO_US,SOFT_TO_C, &
            CHTOT,DENCOR,CVTOT,CSTRF, &
            LMDIM,IRDMAX,CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,CHDEN,SV)
    ENDIF

    ! calculate H - S epsilon |i> for all states
    CALL DETERMINE_HPHI_FAST(HAMILTONIAN, &
         GRID, KPOINTS, LATT_CUR, NONLR_S, NONL_S, T_INFO, P, WT0, W0, WDES, WHPHI, &
         LMDIM, WDES%NSIM, CQIJ,CDIJ, SV, SYMM, UPDATE_EXX )
    CALL STOP_TIMING("G",IO%IU6,"HAMIL")

    CALL PROJECT_ORBITALS_GROUNDSTATE (W0, WHPHI, H_OCC_UNOCC,  BD, BSE_INDEX, DELTA)
    H_OCC_UNOCC=H_OCC_UNOCC-H0_OCC_UNOCC

    ! second order correction
    ! vector(t_n) = vector(t_n-1) + dt i  H vector(t_n-1) + (dt i  H)^2 vector(t_n-1)
    C_OCC_UNOCC2=C_OCC_UNOCC2+0.5_q*(ICI*timestep_v(1,1))**2*H_OCC_UNOCC
    CALL STOP_TIMING("G",IO%IU6,"PROJECT")

    CALL SEPERATOR_TIMING(IO%IU6)
    CALL STOP_TIMING("LOOP",IO%IU6,XMLTAG='total')
#else
!=======================================================================
! propagate for INITIAL_STEPS-1 using timstep/INITIAL_STEPS
! really this makes very little difference but maybe preferable
!=======================================================================
    ! vector(t_n) = vector(t_n-1) + dt i  H vector(t_n-1)
    C_OCC_UNOCC2=C_OCC_UNOCC1+(ICI*timestep_v(1,1)/INITIAL_STEPS)*H_OCC_UNOCC

    DO NELM=2,INITIAL_STEPS
! determine H phi at timestep n
       CALL DETERMINE_ORBITALS (W0, WT0, C_OCC_UNOCC2, BD, BSE_INDEX, DELTA)
       CALL STOP_TIMING("G",IO%IU6,"ORBITALS")

! orthogonalization of this perturbed vector
#ifdef use_orthch
       CALL ORTHCH (WT0%WDES, WT0, WT0%WDES%LOVERL, LMDIM, CQIJ)
       CALL STOP_TIMING("G",IO%IU6,"ORTHCH")
#endif

       IF (UPDATE_CHARGE) THEN
! update charge and then potential
          CALL SET_CHARGE( WT0, WDES, INFO%LOVERL, &
               GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
               LATT_CUR, P, SYMM, T_INFO, &
               CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

          CALL SET_KINEDEN(GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM, &
               T_INFO%NIONS,WT0,WDES,KINEDEN)

          CALL STOP_TIMING("G",IO%IU6,"CHARGE")

          CALL UPDATE_POTENTIAL (  &
               KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,LATT_CUR, &
               T_INFO,INFO,IO, &
               GRID,GRID_SOFT,GRIDC,GRIDUS,C_TO_US,SOFT_TO_C, &
               CHTOT,DENCOR,CVTOT,CSTRF, &
               LMDIM,IRDMAX,CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,CHDEN,SV)
       ENDIF

! calculate H |i> for all states
       CALL DETERMINE_HPHI_FAST(HAMILTONIAN, &
            GRID, KPOINTS, LATT_CUR, NONLR_S, NONL_S, T_INFO, P, WT0, W0, WDES, WHPHI, &
            LMDIM, WDES%NSIM, CQIJ,CDIJ, SV, SYMM, UPDATE_EXX )
       CALL STOP_TIMING("G",IO%IU6,"HAMIL")

       CALL PROJECT_ORBITALS_GROUNDSTATE (W0, WHPHI, H_OCC_UNOCC,  BD, BSE_INDEX, DELTA)
       H_OCC_UNOCC=H_OCC_UNOCC-H0_OCC_UNOCC

       ! determine orbital coefficients at timestep n+1
       H_OCC_UNOCC=C_OCC_UNOCC1+(2._q*ICI*timestep_v(1,1)/INITIAL_STEPS)*H_OCC_UNOCC

       ! move forward in time
       C_OCC_UNOCC1=C_OCC_UNOCC2
       C_OCC_UNOCC2=H_OCC_UNOCC

       SCALAR(NELM-1,IDIR)=ZDOTC(SIZE(C_TRANS_PROP),C_TRANS_PROP(1),1,C_OCC_UNOCC1(1),1)
       CALL STOP_TIMING("G",IO%IU6,"PROJECT")
    ENDDO
    ! C_OCC_UNOCC2 is now propagated for at total of time of timestep_v(1,1)

    CALL SEPERATOR_TIMING(IO%IU6)
    CALL STOP_TIMING("LOOP",IO%IU6,"9xLOOP",XMLTAG='total')
#endif
! ===================================================================================
! C_OCC_UNOCC2 now ready
! now we can run through time development of transition coefficients using alogrithm:
! C_OCC_UNOCC_3 = C_OCC_UNOCC_1 + 2*dt/(i*hbar) H C_OCC_UNOCC_2
! ===================================================================================
    ! just to be save copy C_TRANS_PROP again
    C_OCC_UNOCC1 = C_TRANS_PROP

    ! NOW recalculate it including k-point weights
    CALL INIT_OSCILLATOR (W0, C_TRANS_PROP, LATT_CUR, BD, BSE_INDEX, IDIR, .TRUE.)

    IF (IO%LOPEN) CALL WFORCE(IO%IU6) ! force a write of OUTCAR
!=======================================================================
! act with Hamiltonian on C_OCC_UNOCC2 and store action in H_OCC_UNOCC
!     v(2) = v(2) + [dt/(i*hbar)]^2 H H v(0)
! C_OCC_UNOCC1 stores the orbitals at time step n-1
! C_OCC_UNOCC2 stores the orbitals at time step n
!=======================================================================
    DO NELM=2,NELMGW_ACTUAL

! determine H phi at timestep n
       CALL DETERMINE_ORBITALS (W0, WT0, C_OCC_UNOCC2, BD, BSE_INDEX, DELTA)

! orthogonalization of this perturbed vector
#ifdef use_orthch
       CALL ORTHCH (WT0%WDES, WT0, WT0%WDES%LOVERL, LMDIM, CQIJ)
#endif

       IF (UPDATE_CHARGE) THEN
! update charge and then potential
          CALL SET_CHARGE( WT0, WDES, INFO%LOVERL, &
               GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
               LATT_CUR, P, SYMM, T_INFO, &
               CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

          CALL SET_KINEDEN(GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM, &
               T_INFO%NIONS,WT0,WDES,KINEDEN)

          CALL UPDATE_POTENTIAL (  &
               KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,LATT_CUR, &
               T_INFO,INFO,NOIO, &
               GRID,GRID_SOFT,GRIDC,GRIDUS,C_TO_US,SOFT_TO_C, &
               CHTOT,DENCOR,CVTOT,CSTRF, &
               LMDIM,IRDMAX,CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,CHDEN,SV)
       ENDIF

! calculate H |i> for all states
       CALL DETERMINE_HPHI_FAST(HAMILTONIAN, &
            GRID, KPOINTS, LATT_CUR, NONLR_S, NONL_S, T_INFO, P, WT0, W0, WDES, WHPHI, &
            LMDIM, WDES%NSIM, CQIJ,CDIJ, SV, SYMM , UPDATE_EXX)

       CALL PROJECT_ORBITALS_GROUNDSTATE (W0, WHPHI, H_OCC_UNOCC,  BD, BSE_INDEX, DELTA)
       H_OCC_UNOCC=H_OCC_UNOCC-H0_OCC_UNOCC

       ! determine orbital coefficients at timestep n+1
       H_OCC_UNOCC=C_OCC_UNOCC1+2._q*ICI*timestep_v(1,NELM)* H_OCC_UNOCC

       ! move forward in time
       C_OCC_UNOCC1=C_OCC_UNOCC2
       C_OCC_UNOCC2=H_OCC_UNOCC
       TIME=TIME+timestep_v(1,NELM)

       SCALAR(NELM-1,IDIR)=ZDOTC(SIZE(C_TRANS_PROP),C_TRANS_PROP(1),1,C_OCC_UNOCC1(1),1)

       IF (MOD(NELM,100)==0) THEN
          CALL SEPERATOR_TIMING(IO%IU6)
          CALL STOP_TIMING("LOOP",IO%IU6,"100xLOOP",XMLTAG='total')
          IF (IO%LOPEN) CALL WFORCE(IO%IU6) ! force a write of OUTCAR
       ENDIF
    ENDDO !enddo LOOP NELM

    ! last step's contribution
    SCALAR(NELM-1,IDIR)=ZDOTC(SIZE(C_TRANS_PROP),C_TRANS_PROP(1),1,C_OCC_UNOCC2(1),1)

    DEALLOCATE(C_TRANS_PROP, C_OCC_UNOCC1, C_OCC_UNOCC2, H_OCC_UNOCC, H0_OCC_UNOCC )
!    DEALLOCATE(timestep_v)
    CALL DEALLOCW(WT0)
    CALL DEALLOCW(WHPHI)

! now redistribute the unperturbed orbitals back over orbitals
    DO ISP=1,W0%WDES%ISPIN
       DO K1=1,W0%WDES%NKPTS
          IF (W0%WDES%WTKPT(K1)==0) CYCLE
#ifdef MPI
          IF (MOD(K1-1,W0%WDES%COMM_KINTER%NCPU).NE.W0%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
          CALL SETWDES(W0%WDES,WDES1,K1)
          CALL REDISTRIBUTE_PROJ( ELEMENTS( W0, WDES1, ISP))
          CALL REDISTRIBUTE_PW( ELEMENTS( W0, WDES1, ISP))
       ENDDO
    ENDDO

  END SUBROUTINE PROPAGATE


!************************ SUBROUTINE DETERMINE_HPHI ********************
!
! determine the action of the Hamiltonian on a set of orbitals
!
! d (H + constraint)/ d <i| =  f_i (H |i>  - sum_j S |j><j| H  |i> )
!                           - 1/2 S [H, F] |i>
!
! note that the geometrical gradient is
! d (H + constraint)/ d S <i| = S^-1 d (H + constraint)/ d <i|
!
! TODO: last term is for metals and presently not included
! will be required if we want to do metals
!
! up to second order this is equivalent to
!
! d (H + constraint)/ d <i| =  f_i (H |i>  - epsilon^0_i S |i><i|)
!
! where epsilon^0_i is the groundstate eigenvalue for orbital i
!
! the second routine (presently used)
!
! DETERMINE_HPHI_FAST
!
! the original set is stored in W
! the action of the Hamiltonian is stored in WHPHI
!
!***********************************************************************

  SUBROUTINE DETERMINE_HPHI (HAMILTONIAN, &
      GRID, KPOINTS, LATT_CUR, NONLR_S, NONL_S, T_INFO, P, W, W0, WDES, WHPHI, &
      LMDIM, NSIM, CQIJ,CDIJ, SV, SYMM, UPDATE_CHARGE )
   USE prec
   USE wave_high
   USE lattice
   USE mpimy
   USE nonl_high
   USE pseudo
   USE hamil
   USE hamil_struct_def
   USE fock
   USE lattice
   USE sym_grad
   USE pead, ONLY : PEAD_ACC_CALC_ALL,PEAD_ACC_ADD
   IMPLICIT NONE

   TYPE (ham_handle)  HAMILTONIAN
   TYPE (grid_3d)     GRID
   TYPE (latt)        LATT_CUR
   TYPE (type_info)   T_INFO
   TYPE (potcar)      P(T_INFO%NTYP)
   TYPE (nonlr_struct)NONLR_S
   TYPE (nonl_struct) NONL_S
   TYPE (wavespin)    W            ! current set of orbitals
   TYPE (wavespin)    W0           ! groundstate orbitals
   TYPE (wavespin)    WHPHI        ! derivative d H / d <i|
   TYPE (wavedes)     WDES
   TYPE (kpoints_struct) KPOINTS
   TYPE (symmetry)    SYMM
   LOGICAL   UPDATE_CHARGE         ! use W0 for evaluation of Fock term

   INTEGER NSIM                    ! simultaneously optimised bands
   INTEGER LMDIM                   ! dimension of arrays CQIJ and CDIJ

   RGRID   SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ) ! local potential
   OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
! work array
   GDEF, ALLOCATABLE ::  CHAM(:,:)
   TYPE (wavedes1)    WDES1
   TYPE (wavefun1)    W1(NSIM),WSEARCH
   GDEF,ALLOCATABLE,TARGET :: CPROW(:,:)
   COMPLEX(q),ALLOCATABLE,TARGET:: CF(:,:)    ! stores the plane wave coefficients of the gradient vector
   GDEF,ALLOCATABLE,TARGET      :: CPROF(:,:) ! stores the projected coefficients corresponding to CF
   COMPLEX(q), POINTER :: CF_RED(:,:)
   GDEF      , POINTER :: CPROW_RED(:,:),CPROF_RED(:,:)
   TYPE (wavefuna)    WA                      ! pointer to W%CW(:,:;NK,ISP)
   INTEGER :: NODE_ME, IONODE, NCPU, NSIM_LOCAL, IDUMP
   INTEGER :: I, NB_TOT, NBANDS, ISP, NK, NN, NNMAX, NBANDSK, N, NP, NBAND_LAST_FILLED
   INTEGER :: NGVECTOR, NSIM_, M, ISPINOR, NPOS_RED, NSTRIP_RED, N1, N2
   INTEGER :: NPRO_
   COMPLEX(q) :: CDCHF                        ! double counting correction from HF
   REAL(q) :: DORT
!=======================================================================
! initialise the required variables for MPI
!=======================================================================
#ifdef MPI
   NODE_ME =WDES%COMM_KIN%NODE_ME
   IONODE  =WDES%COMM_KIN%IONODE
   NCPU    =WDES%COMM_INTER%NCPU     ! number of groups (each group holds one band)
#else
   NODE_ME =1
   IONODE  =1
   NCPU    =1
#endif
!=======================================================================
! number of bands treated simultaneously this must be a multiple of NCPU
!=======================================================================
   NSIM_LOCAL=NSIM/NCPU  ! number of bands optimised on this node
   IF (NSIM_LOCAL*NCPU /= NSIM) THEN
      CALL vtutor%bug("internal ERROR in EDWAV: DETERMINE_HPHI is not correct " // str(NSIM), __FILE__, __LINE__)
   ENDIF

   IDUMP=1
#ifdef MPI
   IF (NODE_ME/=IONODE) IDUMP=0
#endif

   CALL SETWDES(WDES,WDES1,0)
   DO I=1,NSIM_LOCAL
      CALL NEWWAV(W1(I) , WDES1, .TRUE.)
   ENDDO
   CALL NEWWAV(WSEARCH, WDES1, .TRUE.)

   NB_TOT=WDES%NB_TOT
   NBANDS=WDES%NBANDS
   ALLOCATE(CPROW(WDES%NPROD,NBANDS), CF(WDES%NRPLWV,NSIM_LOCAL),CPROF(WDES%NPROD,NSIM_LOCAL), CHAM(WDES%NB_TOT,WDES%NB_TOT))

   WHPHI%CW   =0
   WHPHI%CPROJ=0

   CDCHF=0
!=======================================================================
! first loop
! calculate gradient (for wavefunctions)
!=======================================================================
   CALL PEAD_ACC_CALC_ALL(W,KPOINTS,P,CQIJ,LATT_CUR,T_INFO)
!=======================================================================
! start with HF part and store results WHPHI
!=======================================================================
   IF (LHFCALC) THEN
      WHPHI%CW=0
      DO ISP=1,WDES%ISPIN
      DO NK=1,WDES%NKPTS
#ifdef MPI
         IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
         CALL SETWDES(WDES,WDES1,NK)

         IF (NONLR_S%LREAL) THEN
            CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES)
         ELSE
            CALL PHASE(WDES,NONL_S,NK)
         ENDIF
         NBANDSK=CEILING(REAL(W%WDES%NB_TOTK(NK,ISP))/REAL(W%WDES%NB_PAR))
         ! we only need the dynamics for occupied states
         ! because of parallel version maybe a little bit more needs to be calculated
         NBANDSK=MIN(NBANDSK,LAST_FILLED_XI(W,NK,ISP)/WDES%NB_PAR)

         DO NN=1,NBANDSK,NSIM_LOCAL
            NNMAX=MIN(NN+NSIM_LOCAL-1,NBANDSK)
            CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W,   &
                    NONLR_S, NONL_S, NK, ISP, NN, NNMAX-NN+1, &
                    WHPHI%CW(:,NN:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF, LSYMGRAD=LSYMGRAD)
         ENDDO
      ENDDO
      ENDDO
      IF (LSYMGRAD) THEN
         CALL APPLY_SMALL_SPACE_GROUP_OP( W, WHPHI, NONLR_S, NONL_S,P, T_INFO%NIONS, LATT_CUR, SYMM, CQIJ, .FALSE. , -1)
      ENDIF
   ENDIF
!=======================================================================
! DFT part optimized to yield minimum number of operations
!=======================================================================
spin:   DO ISP=1,WDES%ISPIN
kpoint: DO NK=1,WDES%NKPTS
#ifdef MPI
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
      CHAM=0
      CALL SETWDES(WDES,WDES1,NK)

      WA=ELEMENTS(W, WDES1, ISP)

      IF (WDES%DO_REDIS) THEN
        CALL SET_GPOINTER(CPROW_RED, WDES1%NPROD_RED, NB_TOT, CPROW(1,1))
        CALL SET_WPOINTER(CF_RED,   WDES1%NRPLWV_RED, NSIM, CF(1,1))
        CALL SET_GPOINTER(CPROF_RED, WDES1%NPROD_RED, NSIM, CPROF(1,1))
      ELSE
        CPROW_RED => CPROW(:,:)
        CF_RED    => CF(:,:)
        CPROF_RED => CPROF(:,:)
      ENDIF

      NGVECTOR=WDES1%NGVECTOR
      CALL OVERL(WDES1, W%WDES%LOVERL,LMDIM,CQIJ(1,1,1,ISP), W%CPROJ(1,1,NK,ISP),CPROW(1,1))

!-----------------------------------------------------------------------
! calculate the phase-factor on the compressed grid (only non local)
!-----------------------------------------------------------------------
      IF (NONLR_S%LREAL) THEN
         CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES)
      ELSE
         CALL PHASE(WDES,NONL_S,NK)
      ENDIF

      NBANDSK=CEILING(REAL(W%WDES%NB_TOTK(NK,ISP))/REAL(W%WDES%NB_PAR))
      ! we only need the dynamics for occupied states
      ! because of parallel version maybe a little bit more needs to be calculated
      NBANDSK=MIN(NBANDSK,LAST_FILLED_XI(W,NK,ISP)/WDES%NB_PAR)
      NBAND_LAST_FILLED=LAST_FILLED_XI_NOMOD(W,NK,ISP)

      ! redistribute over plane wave coefficients
      IF (WDES%DO_REDIS) THEN
         IF (W%WDES%LOVERL) CALL REDIS_PROJ(WDES1, NBANDS, CPROW(1,1))
         CALL REDIS_PW  (WDES1, NBANDS, W%CW(1,1,NK,ISP))
         CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
      ENDIF
!=======================================================================
  bands: DO NN=1,NBANDSK,NSIM_LOCAL
!=======================================================================
      CF=0
      CPROF=0

      NNMAX=MIN(NN+NSIM_LOCAL-1,NBANDSK)
      NSIM_=NNMAX-NN+1

      DO N=NN,NNMAX ! each step in the loop copies NSIM_LOCAL/NSIM bands !
         NP=N-NN+1
         CF(:,NP)=WHPHI%CW(:,N,NK,ISP)
         CALL SETWDES(WDES, WDES1, NK)

         ! use as source the wavefunctions distributed over bands (if possible)
         DO M=1,WDES%NRPLWV
            W1(NP)%CW(M)=W%CW(M,N,NK,ISP)
         ENDDO
         DO M=1,WDES%NPROD
            W1(NP)%CPROJ(M)=W%CPROJ(M,N,NK,ISP)
         ENDDO
         ! redistribute over bands
         IF (WDES%DO_REDIS) THEN
            CALL REDIS_PW(WDES1, 1, W1(NP)%CW(1))
            CALL REDIS_PROJ(WDES1, 1, W1(NP)%CPROJ(1))
         ENDIF

         DO ISPINOR=0,WDES%NRSPINORS-1
            CALL FFTWAV(NGVECTOR, WDES%NINDPW(1,NK),W1(NP)%CR(1+ISPINOR*GRID%MPLWV),W1(NP)%CW(1+ISPINOR*NGVECTOR),GRID)
         ENDDO

         IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
            CALL HAMILT_LOCAL_TAU(W1(NP), SV, LATT_CUR, HAMILTONIAN%MU, ISP, CF(:,NP), LHFCALC)
         ELSE
            CALL HAMILT_LOCAL(W1(NP), SV, ISP, CF(:,NP), LHFCALC)
         ENDIF

         WDES1%NBANDS=1    ! WDES1%NBANDS is used only here to fake OVERL
         CALL OVERL(WDES1, .TRUE., LMDIM, CDIJ(1,1,1,ISP), W1(NP)%CPROJ(1), CPROF(1,NP))

         CALL PEAD_ACC_ADD(CF(:,NP),CPROF(:,NP),N,NK,ISP)

         IF ((W%WDES%NB_LOW+(N-1)*W%WDES%NB_PAR)>W%WDES%NB_TOTK(NK,ISP)) THEN
            CF(:,NP)=0
            CPROF(:,NP)=0
            W1(NP)%CW=0
            W1(NP)%CR=0
            W1(NP)%CPROJ=0
         ENDIF

         ! redistribute over plane wave coefficients
         IF (WDES%DO_REDIS) THEN
            CALL REDIS_PW(WDES1, 1, CF(1,NP))
            CALL REDIS_PROJ(WDES1, 1, CPROF(1,NP))
         ENDIF
      ENDDO

! calculate the lower triangle of the Hamilton matrix
! in the subspace spanned by the present set of wavefunctions

      NPOS_RED  =(NN-1)*NCPU+1
      NSTRIP_RED=NSIM_*NCPU
! TODO here we calculate too much, CHAM is only required up to NBAND_LAST_FILLED
      CALL ORTH1('L', &
        WA%CW_RED(1,1),CF(1,1),WA%CPROJ_RED(1,1), &
        CPROF(1,1),NB_TOT, &
        NPOS_RED, NSTRIP_RED, WDES1%NPL_RED, WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CHAM(1,1))

      CALLMPI( M_sum_g(WDES%COMM_KIN,CHAM(1,NPOS_RED),NB_TOT*NSTRIP_RED))

!     fill in upper triangle
!     in the parallel version part of the upper triangle has been destroyed
!     and needs to be refilled
      DO N1=NPOS_RED,NPOS_RED+NSTRIP_RED-1
         DO N2=N1+1,NB_TOT
            CHAM(N1,N2)=GCONJG(CHAM(N2,N1))
         ENDDO
         W%CELTOT(N1,NK,ISP)=CHAM(N1,N1)
         DO N2=1,N1-1
            CHAM(N2,N1)=GCONJG(CHAM(N1,N2))
         ENDDO
      ENDDO

      ! set all matrix elements beyond NB_TOTK to zero
      IF (W%WDES%NB_TOTK(NK,ISP)<NB_TOT) THEN
         CHAM(:,W%WDES%NB_TOTK(NK,ISP)+1:NB_TOT)=0
         CHAM(W%WDES%NB_TOTK(NK,ISP)+1:NB_TOT,:)=0
      ENDIF

!      CALL DUMP_HAM( "Hamilton matrix",W%WDES, CHAM)
!
! subtract from H phi the components parallel to the subspace
! spanned by the other occupied orbitals
      IF ( WDES1%NPL_RED/=0) &
! subtract sum_n2 | phi_n2> < phi_n2 | H_local  | phi_n1 >
      CALL GGEMM('N', 'N', m_  WDES1%NPL_RED , NSTRIP_RED, NBAND_LAST_FILLED , -one, &
     &               WA%CW_RED(1,1), m_ WDES1%NRPLWV_RED, CHAM(1,NPOS_RED), NB_TOT,  &
     &               one, CF(1,1), m_ WDES1%NRPLWV_RED)

! subtract sum_n2  Q_ij  < p_j| phi_n2 > < phi_n2 |  H_local | phi_n1 > from orbital character
      IF (W%WDES%LOVERL) THEN
         CALL GGEMM('N', 'N',  WDES1%NPRO_RED, NSTRIP_RED, NBAND_LAST_FILLED, -one, &
     &               CPROW_RED(1,1), WDES1%NPROD_RED, CHAM(1,NPOS_RED), NB_TOT,  &
     &               one, CPROF(1,1), WDES1%NPROD_RED)
      ENDIF
      ! redistribute the gradient over bands
      DO N=NN,NNMAX
         NP=N-NN+1
         IF (WDES%DO_REDIS) THEN
            CALL REDIS_PW(WDES1, 1, CF(1,NP))
            CALL REDIS_PROJ(WDES1, 1, CPROF(1,NP))
         ENDIF
      ENDDO

    innerband: DO N=NN,NNMAX
      IF ( (W%WDES%NB_LOW+(N-1)*W%WDES%NB_PAR)>W%WDES%NB_TOTK(NK,ISP)) CYCLE

      NP=N-NN+1
!-----------------------------------------------------------------------
! add the non local contribution to CF
!-----------------------------------------------------------------------
      IF (NONLR_S%LREAL) THEN
         WSEARCH%CR(:)=0
         CALL RACC0(NONLR_S,WDES1,CPROF(1,NP),WSEARCH%CR(1))
         DO ISPINOR=0,WDES1%NRSPINORS-1
            CALL FFTEXT(NGVECTOR,WDES%NINDPW(1,NK),WSEARCH%CR(1+ISPINOR*WDES1%GRID%MPLWV),CF(1+ISPINOR*NGVECTOR,NP),GRID,.TRUE.)
         ENDDO
      ELSE
        CALL VNLAC0(NONL_S,WDES1,CPROF(1,NP),CF(1,NP))
      ENDIF
!-----------------------------------------------------------------------
! calculate expected energy-change into search direction
! and store the derivative
!-----------------------------------------------------------------------
      ! gradient should be orthogonal to present wavefunction
      DORT=0
      DO M=1,WDES1%NPL
         DORT  =DORT+ CF(M,NP)*CONJG(W1(NP)%CW(M))
      ENDDO
      CALLMPI( M_sum_d(WDES%COMM_INB, DORT, 1))
      IF (ABS(DORT)>1E-4) THEN
         CALL vtutor%bug("EDWAV: internal error, the gradient is not orthogonal " // str(NK) // " " &
            // str(NP) // " " // str(DORT), __FILE__, __LINE__)
      ENDIF
!-----------------------------------------------------------------------
! move CF to WSEARCH
!-----------------------------------------------------------------------
      WSEARCH%CW=CF(:,NP)
!----------------------------------------------------------------------
!  calculate the character of the orbital
!----------------------------------------------------------------------
      DO ISPINOR=0,WDES1%NRSPINORS-1
         CALL FFTWAV(NGVECTOR,WDES1%NINDPW(1),WSEARCH%CR(1+ISPINOR*WDES1%GRID%MPLWV),WSEARCH%CW(1+ISPINOR*NGVECTOR),GRID)
      ENDDO
      IF (NONLR_S%LREAL) THEN
         CALL RPRO1(NONLR_S,WDES1,WSEARCH)
      ELSE
         CALL PROJ1(NONL_S,WDES1,WSEARCH)
      ENDIF

      DO M=1,WDES1%NPL
         WHPHI%CW(M,N,NK,ISP)=WSEARCH%CW(M)*W%FERWE(N,NK,ISP)
      ENDDO

      DO NPRO_=1,WDES%NPRO
         WHPHI%CPROJ(NPRO_,N,NK,ISP)=WSEARCH%CPROJ(NPRO_)*W%FERWE(N,NK,ISP)
      ENDDO
!=======================================================================

    ENDDO innerband
  ENDDO bands
!=======================================================================

      IF (WDES%DO_REDIS) THEN
! redistribute CPROJ over  bands
        CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
! gradient distributed over plane wave coefficients
        CALL REDIS_PW  (WDES1, NBANDS, WHPHI%CW   (1,1,NK,ISP))
      ENDIF
!=======================================================================
   ENDDO kpoint
   ENDDO spin


   DEALLOCATE(CPROW, CF, CPROF, CHAM)

   DO I=1,NSIM_LOCAL
      CALL DELWAV(W1(I), .TRUE.)
   ENDDO
   CALL DELWAV(WSEARCH, .TRUE.)

  END SUBROUTINE DETERMINE_HPHI


!***********************************************************************
!
! simpler and faster routines: determine
!
! d (H + constraint)/ d <i| =  f_i (H |i>  - epsilon^0_i S |i><i|)
!
! where epsilon^0_i is the groundstate eigenvalue for orbital i
! this is only correct to second order but much simpler
!
!***********************************************************************

  SUBROUTINE DETERMINE_HPHI_FAST (HAMILTONIAN, &
      GRID, KPOINTS, LATT_CUR, NONLR_S, NONL_S, T_INFO, P, W, W0, WDES, WHPHI, &
      LMDIM, NSIM, CQIJ,CDIJ, SV, SYMM, UPDATE_EXX )
   USE prec
   USE wave_high
   USE lattice
   USE mpimy
   USE nonl_high
   USE pseudo
   USE hamil
   USE hamil_struct_def
   USE fock
   USE lattice
   USE sym_grad
   USE pead, ONLY : PEAD_ACC_CALC_ALL
   IMPLICIT NONE

   TYPE (ham_handle)  HAMILTONIAN
   TYPE (grid_3d)     GRID
   TYPE (latt)        LATT_CUR
   TYPE (type_info)   T_INFO
   TYPE (potcar)      P(T_INFO%NTYP)
   TYPE (nonlr_struct)NONLR_S
   TYPE (nonl_struct) NONL_S
   TYPE (wavespin)    W            ! current set of orbitals
   TYPE (wavespin)    W0           ! groundstate orbitals
   TYPE (wavespin)    WHPHI        ! derivative d H / d <i|
   TYPE (wavedes)     WDES
   TYPE (kpoints_struct) KPOINTS
   TYPE (symmetry)    SYMM
   LOGICAL   UPDATE_EXX            ! use W0 for evaluation of Fock term
   INTEGER NSIM                    ! simultaneously optimised bands
   INTEGER LMDIM                   ! dimension of arrays CQIJ and CDIJ

   RGRID   SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ) ! local potential
   OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
! work array
   TYPE (wavedes1)    WDES1
   TYPE (wavefun1)    W1(NSIM)
   TYPE (wavefuna)    WH           ! stores the (H -epsilon S) |phi> except for Fock part
   REAL(q) :: EVALUE(NSIM)         ! eigenvalue used above

   INTEGER :: I, ISP, NK, NN, NNMAX, NBANDSK, N, NP
   INTEGER :: NSIM_                ! number of bands done currently
   COMPLEX(q) :: CDCHF             ! double counting correction from HF
!=======================================================================
! number of bands treated simultaneously this must be a multiple of NCPU
!=======================================================================
   CALL SETWDES(WDES,WDES1,0)
   DO I=1,NSIM
      CALL NEWWAV(W1(I) , WDES1, .TRUE.)
   ENDDO
   CALL NEWWAVA(WH, WDES1, NSIM)

   WHPHI%CW   =0
   WHPHI%CPROJ=0

   CDCHF=0
!=======================================================================
! start with HF part and store results WHPHI
!=======================================================================
! TODO: need to ask Martijn what this about
   CALL PEAD_ACC_CALC_ALL(W,KPOINTS,P,CQIJ,LATT_CUR,T_INFO)

   IF (LHFCALC) THEN
      DO ISP=1,WDES%ISPIN
      DO NK=1,WDES%NKPTS
#ifdef MPI
         IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
         CALL SETWDES(WDES,WDES1,NK)

         IF (NONLR_S%LREAL) THEN
            CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES)
         ELSE
            CALL PHASE(WDES,NONL_S,NK)
         ENDIF
         NBANDSK=CEILING(REAL(W%WDES%NB_TOTK(NK,ISP))/REAL(W%WDES%NB_PAR))
         ! we only need the dynamics for occupied states
         ! because of parallel version maybe a little bit more needs to be calculated
         NBANDSK=MIN(NBANDSK,LAST_FILLED_XI(W,NK,ISP)/WDES%NB_PAR)
         DO NN=1,NBANDSK,NSIM
            NNMAX=MIN(NN+NSIM-1,NBANDSK)
            IF (.NOT. UPDATE_EXX) THEN
               DO N=NN,NNMAX
                  NP=N-NN+1
                  CALL W1_COPY(ELEMENT(W,WDES1,N,ISP),W1(NP))
               ENDDO
               ! use W0 for evaluation of fock exchange
               CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W0,   &
                    NONLR_S, NONL_S, NK, ISP, NN, NNMAX-NN+1, &
                    WHPHI%CW(:,NN:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF, W1, LSYMGRAD=LSYMGRAD)
            ELSE
               CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W,   &
                    NONLR_S, NONL_S, NK, ISP, NN, NNMAX-NN+1, &
                    WHPHI%CW(:,NN:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF, LSYMGRAD=LSYMGRAD)
            ENDIF
         ENDDO
      ENDDO
      ENDDO
      IF (LSYMGRAD) THEN
         CALL APPLY_SMALL_SPACE_GROUP_OP( W, WHPHI, NONLR_S, NONL_S,P, T_INFO%NIONS, LATT_CUR, SYMM, CQIJ, .FALSE. , -1)
      ENDIF
   ENDIF
!=======================================================================
! DFT part
!=======================================================================
spin:   DO ISP=1,WDES%ISPIN
kpoint: DO NK=1,WDES%NKPTS
#ifdef MPI
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
      CALL SETWDES(WDES,WDES1,NK)

!-----------------------------------------------------------------------
! calculate the phase-factor on the compressed grid (only non local)
!-----------------------------------------------------------------------
      IF (NONLR_S%LREAL) THEN
         CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES)
      ELSE
         CALL PHASE(WDES,NONL_S,NK)
      ENDIF

      NBANDSK=CEILING(REAL(W%WDES%NB_TOTK(NK,ISP))/REAL(W%WDES%NB_PAR))
      ! we only need the dynamics for occupied states
      ! because of parallel version maybe a little bit more needs to be calculated
      NBANDSK=MIN(NBANDSK,LAST_FILLED_XI(W,NK,ISP)/WDES%NB_PAR)

!=======================================================================
  bands: DO NN=1,NBANDSK,NSIM
!=======================================================================
      NNMAX=MIN(NN+NSIM-1,NBANDSK)
      NSIM_=NNMAX-NN+1

      ! unmark all bands
      W1%LDO  =.FALSE.

      ! now fill relevant bands between NN and NNMAX to W1
      DO N=NN,NNMAX
         NP=N-NN+1  ! local index

         CALL W1_COPY(ELEMENT(W, WDES1, N, ISP), W1(NP))
         CALL FFTWAV_W1(W1(NP))
         W1(NP)%LDO=.TRUE.  ! mark the band if it is treated currently

         EVALUE(NP)=W%CELEN(N,NK,ISP)
      ENDDO

      !  calculate H-epsilon S | psi > and store in WH
      IF (ASSOCIATED(HAMILTONIAN%MU)) THEN
         CALL HAMILTMU_TAU(WDES1, W1, NONLR_S, NONL_S, EVALUE, &
              &     CDIJ, CQIJ, SV, LATT_CUR, HAMILTONIAN%MU, ISP, WH)
      ELSE
         CALL HAMILTMU(WDES1,W1,NONLR_S,NONL_S,EVALUE, &
              &     CDIJ,CQIJ, SV, ISP, WH)
      ENDIF

!-----------------------------------------------------------------------
! store final results back to WHPHI
! TODO: for DFT we could store the entire contribution
! already above into WHPHI
!-----------------------------------------------------------------------
      DO N=NN,NNMAX
         NP=N-NN+1  ! local index
         IF (LHFCALC) THEN
            WHPHI%CW(:,N,NK,ISP)=WH%CW(:,NP)+WHPHI%CW(:,N,NK,ISP)
         ELSE
            WHPHI%CW(:,N,NK,ISP)=WH%CW(:,NP)
         ENDIF
         ! currently the wave function characters are not used, set to 0
         WHPHI%CPROJ(:,N,NK,ISP)=0
!=======================================================================
      ENDDO
!=======================================================================
   ENDDO bands
   ENDDO kpoint
   ENDDO spin

   DO I=1,NSIM
      CALL DELWAV(W1(I), .TRUE.)
   ENDDO
   CALL DELWAVA(WH)

  END SUBROUTINE DETERMINE_HPHI_FAST


!*********************************************************************
!
! initialize the initial orbitals for the time evolution scheme
! using the precalculated transition amplitudes
!
!*********************************************************************

  SUBROUTINE INIT_OSCILLATOR (W, C_OCC_UNOCC, LATT_CUR, BD, BSE_INDEX, IDIR, LWEIGHT)
    USE wave_high
    USE lattice
    TYPE (wavespin),INTENT(IN)  :: W          ! unperturbed orbitals
    COMPLEX(q)     C_OCC_UNOCC(:)             ! coefficients determining the perturbed orbitals
    TYPE(banddesc)    :: BD(:)                ! descriptor for bands
    TYPE(bse_matrix_index) :: BSE_INDEX
    TYPE (latt) LATT_CUR                      ! basis set vectors
    INTEGER                        IDIR       ! direction of field
    LOGICAL                        LWEIGHT    ! add k-point weight
! local variables
    INTEGER ISP
    INTEGER K1
    INTEGER N1, N3
    INTEGER NB1_INTO_TOT, NB3_INTO_TOT, NCV13
    REAL(q) :: WEIGHT
    INTEGER NK1_IN_KPOINTS_FULL_ORIG
    GDEF       :: CDER_BETWEEN_STATE(3)
    TYPE (wavedes1) :: WDES1


    C_OCC_UNOCC=0

    DO ISP=1,W%WDES%ISPIN
    DO K1=1,W%WDES%NKPTS
       IF (W%WDES%WTKPT(K1)==0) CYCLE
#ifdef MPI
       IF (MOD(K1-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
       CALL SETWDES(W%WDES,WDES1,K1)

       ! determined the index of this k-point in the original full k-point grid
       NK1_IN_KPOINTS_FULL_ORIG=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,K1),KPOINTS_FULL_ORIG)
       ! empty orbital N3
       DO N3=1,BD(ISP)%CBMAX-BD(ISP)%CBMIN+1
          NB3_INTO_TOT=N3+BD(ISP)%CBMIN-1
          IF (FILLED_XI_ORBITAL(W%FERTOT(NB3_INTO_TOT,K1,ISP))) CYCLE

          ! occupied orbital N1
          DO N1=1,BD(ISP)%VBMAX-BD(ISP)%VBMIN+1
             NB1_INTO_TOT=N1+BD(ISP)%VBMIN-1
             IF (EMPTY_XI_ORBITAL(W%FERTOT(NB1_INTO_TOT,K1,ISP))) CYCLE

             WEIGHT=(W%FERTOT(NB1_INTO_TOT,K1,ISP)-W%FERTOT(NB3_INTO_TOT, K1, ISP))
             ! G-> 0 limit, consider all three directions
             CALL  CDER_BETWEEN_STATES_ROTATED(CDER_BETWEEN_STATE,LATT_CUR, NK1_IN_KPOINTS_FULL_ORIG, ISP, NB1_INTO_TOT, NB3_INTO_TOT)

             NCV13=BSE_INDEX%INDEX(NB1_INTO_TOT-BD(ISP)%VBMIN+1, NB3_INTO_TOT-BD(ISP)%CBMIN+1, K1, ISP)
             IF (NCV13>0) THEN
                ! CDER_BETWEEN_STATE =  <v_k1,n1|- i d/dq_j | c_k1+q,n3>
                ! the orbital after the pulse is (Fermi's golden rule)
                !        \sum_n3 |c_k1+q,n3> <c_k1+q,n3|- i d/dq_j | v_k1,n1>
                ! so we need to conjugate here (alas why not change CDER above,however, never change a running code)
                C_OCC_UNOCC(NCV13)=GCONJG(CDER_BETWEEN_STATE(IDIR))*AUTOA/(2._q*HSQDTM)*(0.0_q,1.0_q)
                IF (LWEIGHT) THEN
                   C_OCC_UNOCC(NCV13)=C_OCC_UNOCC(NCV13)*W%WDES%WTKPT(K1)
                ENDIF
             ENDIF

          ENDDO  ! enddo for N3-loop
       ENDDO     ! enddo for N1-loop

    ENDDO  ! enddo for K1-loop
    ENDDO  ! enddo spin loop

    CALLMPI( M_sum_z(W%WDES%COMM_KINTER, C_OCC_UNOCC, SIZE(C_OCC_UNOCC)))

  END SUBROUTINE INIT_OSCILLATOR

!*********************************************************************
!
! determine perturbed orbitals from given coefficient matrix
!
!  | i> = |i^0> + \sum_a c_ia |a^0>
!
! input the the routine are the coefficients c_ia and the
! unperturbed orbitals
!
! important note: the routine assumes that the orbitals W0
! are distributed over plane wave coefficients
!
!*********************************************************************

  SUBROUTINE DETERMINE_ORBITALS (W0, WT0, C_OCC_UNOCC, BD, BSE_INDEX, DELTA)
    USE wave_high
    USE lattice
    USE pead, ONLY : LUSEPEAD
    TYPE (wavespin),INTENT(IN)  :: W0         ! unperturbed orbitals
    TYPE (wavespin),INTENT(INOUT) :: WT0      ! perturbed orbitals  at timestep t_i
    COMPLEX(q)        :: C_OCC_UNOCC(:)       ! coefficients determining the perturbed orbitals
    TYPE(banddesc)    :: BD(:)                ! descriptor for bands
    TYPE(bse_matrix_index) :: BSE_INDEX
    REAL (q)                       DELTA      ! magnitude of perturbation
! local variables
    INTEGER ISP
    INTEGER K1
    INTEGER N1, N3
    INTEGER NB1_INTO_TOT, NB3_INTO_TOT, NCV13
    REAL(q) :: WEIGHT
    INTEGER NK1_IN_KPOINTS_FULL_ORIG
    GDEF,ALLOCATABLE,TARGET::  CHAM(:,:)      ! matrix of transition propabilities
    TYPE (wavedes1)    WDES1

    ALLOCATE(CHAM(W0%WDES%NB_TOT,W0%WDES%NB_TOT))

    DO ISP=1,W0%WDES%ISPIN
    DO K1=1,W0%WDES%NKPTS
       IF (W0%WDES%WTKPT(K1)==0) CYCLE
#ifdef MPI
       IF (MOD(K1-1,W0%WDES%COMM_KINTER%NCPU).NE.W0%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
       WT0%CW(:,:,K1,ISP)   =W0%CW(:,:,K1,ISP)
       WT0%CPROJ(:,:,K1,ISP)=W0%CPROJ(:,:,K1,ISP)
       WT0%FERTOT(:,K1,ISP) =W0%FERTOT(:,K1,ISP)
       WT0%CELTOT(:,K1,ISP) =W0%CELTOT(:,K1,ISP)

       CALL SETWDES(W0%WDES,WDES1,K1)

       ! fill elements from stored coefficient matrix into CHAM
       CHAM=0
       DO N3=1,BD(ISP)%CBMAX-BD(ISP)%CBMIN+1
          NB3_INTO_TOT=N3+BD(ISP)%CBMIN-1
          IF (FILLED_XI_ORBITAL(W0%FERTOT(NB3_INTO_TOT,K1,ISP))) CYCLE

          ! occupied orbital N1
          DO N1=1,BD(ISP)%VBMAX-BD(ISP)%VBMIN+1
             NB1_INTO_TOT=N1+BD(ISP)%VBMIN-1
             IF (EMPTY_XI_ORBITAL(W0%FERTOT(NB1_INTO_TOT,K1,ISP))) CYCLE

             NCV13=BSE_INDEX%INDEX(NB1_INTO_TOT-BD(ISP)%VBMIN+1, NB3_INTO_TOT-BD(ISP)%CBMIN+1, K1, ISP)
             IF (NCV13>0) THEN
                CHAM(NB3_INTO_TOT , NB1_INTO_TOT) = C_OCC_UNOCC(NCV13)*DELTA
                ! make matrix anti Hermitian = approximately unitary
                ! since unoccupied states are not used, not really that important I guess
                CHAM(NB1_INTO_TOT , NB3_INTO_TOT) =-GCONJG(C_OCC_UNOCC(NCV13))*DELTA
             ENDIF
          ENDDO  ! enddo for N3-loop
       ENDDO     ! enddo for N1-loop

       ! finally set diagonal component
       DO N1=1,W0%WDES%NB_TOT
          CHAM(N1,N1)=1
       ENDDO

       ! and normalize the vectors to one
       ! this will yield orthogonal vectors to second order
!       DO N1=1,W0%WDES%NB_TOT
!          CHAM(:,N1)=CHAM(:,N1)/SQRT(REAL(DOT_PRODUCT(CHAM(:,N1),GCONJG(CHAM(:,N1))),q))
!       ENDDO

       !linear combination --> results in perturbed orbitals

       ! TODO: the matrix we pass to LINCOM is a unitary matrix for small
       ! delta, however it is not exactly unitary
       ! currently we fix that a posteriori by calling ORTHCH in PROPAGATE
       ! also we need only the occupied orbitals and not the unoccupied ones
       CALL LINCOM('F',WT0%CW(:,:,K1,ISP),WT0%CPROJ(:,:,K1,ISP),CHAM(1,1), &
            W0%WDES%NB_TOTK(K1,ISP),W0%WDES%NB_TOTK(K1,ISP), &
            WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,W0%WDES%NB_TOT, &
            WT0%CW(:,:,K1,ISP),WT0%CPROJ(:,:,K1,ISP))

! now redistribute over orbitals
       CALL REDISTRIBUTE_PROJ( ELEMENTS( WT0, WDES1, ISP))
       CALL REDISTRIBUTE_PW( ELEMENTS( WT0, WDES1, ISP))

    ENDDO  ! enddo for K1-loop
    ENDDO  ! enddo spin loop

    DEALLOCATE(CHAM)
#ifdef MPI
    IF ((LHFCALC.OR.LUSEPEAD()).AND. WT0%WDES%COMM_KINTER%NCPU.GT.1) THEN
       CALL KPAR_SYNC_ALL (WT0%WDES, WT0)
    ENDIF
#endif
  END SUBROUTINE DETERMINE_ORBITALS


!*********************************************************************
!
! project the perturbed orbitals onto groundstate orbitals
!
! c_ia        =  < i^0| j >
! H_OCC_UNOCC = < W_i | WHPHI_j>
!
! the overlap operator is not included here
!
! important note: the routine assumes that the grounstate orbitals W
! are distributed over plane wave coefficients
! whereas the WHPHI are assumed to be distributed over orbitals
! upon exit WHPHI is distributed over plane wave coefficients
!
!*********************************************************************

  SUBROUTINE PROJECT_ORBITALS_GROUNDSTATE(W0, WHPHI, H_OCC_UNOCC, BD, BSE_INDEX, DELTA )

    TYPE (wavespin),INTENT(IN)  :: W0          ! unperturbed orbitals
    TYPE (wavespin)             :: WHPHI       ! H - epsilon S | phi> at current timestep
    COMPLEX(q) ::  H_OCC_UNOCC(:)              ! action of Hamiltonian on perturbed orbitals
    TYPE(banddesc)    :: BD(:)                 ! descriptor for bands
    TYPE(bse_matrix_index) :: BSE_INDEX
    REAL(q) :: DELTA
!    local
    GDEF,ALLOCATABLE,TARGET::  CHAM(:,:)      ! matrix of transition propabilities
    TYPE (wavedes1)    WDES1
    INTEGER ISP, K1
    TYPE (wavefuna)    WA             ! pointer towards one set of orbitals in W
    TYPE (wavefuna)    WHPHIA           ! pointer towards one set of orbitals in WHPHI
    INTEGER :: N3, N1, NB1_INTO_TOT, NB3_INTO_TOT, NCV13

    ALLOCATE(CHAM(W0%WDES%NB_TOT,W0%WDES%NB_TOT))
    CALL SETWDES(W0%WDES,WDES1,0)

    H_OCC_UNOCC=0

    DO ISP=1,W0%WDES%ISPIN
    DO K1=1,W0%WDES%NKPTS
       IF (W0%WDES%WTKPT(K1)==0) CYCLE
#ifdef MPI
       IF (MOD(K1-1,W0%WDES%COMM_KINTER%NCPU).NE.W0%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
       CALL SETWDES(W0%WDES,WDES1,K1)

       WA  =ELEMENTS(W0, WDES1, ISP)
       WHPHIA=ELEMENTS(WHPHI, WDES1, ISP)

       ! currently H phi is distributed over orbitals
       ! ORTH2 required distribution over plane wave coefficients
       CALL REDISTRIBUTE_PW( WHPHIA)

       CHAM=0
       ! TODO: we only want the block between occupied and unoccupied orbitals
       ! currently we calculate the entire matrix
       CALL ORTH2( &
            WA%CW(1,1), WHPHIA%CW(1,1), WA%CPROJ(1,1),WA%CPROJ(1,1), &   ! the CPROJ part is not used here
            W0%WDES%NB_TOT, 1, W0%WDES%NB_TOT, WDES1%NPL_RED, 0, WDES1%NRPLWV_RED, WDES1%NPROD_RED, CHAM(1,1))

       CALLMPI( M_sum_g(W0%WDES%COMM_KIN,CHAM(1,1),W0%WDES%NB_TOT*W0%WDES%NB_TOT))

!       IF (K1==2) CALL DUMP_HAM( "projected orbitals",W0%WDES, CHAM)

       ! store the result in compressed manner in
       ! empty orbital N3
       DO N3=1,BD(ISP)%CBMAX-BD(ISP)%CBMIN+1
          NB3_INTO_TOT=N3+BD(ISP)%CBMIN-1
          IF (FILLED_XI_ORBITAL(W0%FERTOT(NB3_INTO_TOT,K1,ISP))) CYCLE

          ! occupied orbital N1
          DO N1=1,BD(ISP)%VBMAX-BD(ISP)%VBMIN+1
             NB1_INTO_TOT=N1+BD(ISP)%VBMIN-1
             IF (EMPTY_XI_ORBITAL(W0%FERTOT(NB1_INTO_TOT,K1,ISP))) CYCLE

             NCV13=BSE_INDEX%INDEX(NB1_INTO_TOT-BD(ISP)%VBMIN+1, NB3_INTO_TOT-BD(ISP)%CBMIN+1, K1, ISP)

             IF (NCV13>0) THEN
                H_OCC_UNOCC(NCV13)=CHAM(NB3_INTO_TOT , NB1_INTO_TOT)/DELTA
             ENDIF
          ENDDO  ! enddo for N3-loop
       ENDDO     ! enddo for N1-loop

    ENDDO
    ENDDO

    DEALLOCATE(CHAM)

    CALLMPI( M_sum_z(W0%WDES%COMM_KINTER, H_OCC_UNOCC, SIZE(H_OCC_UNOCC)))

  END SUBROUTINE PROJECT_ORBITALS_GROUNDSTATE


!*********************************************************************
!
! update potential
!
!*********************************************************************


    SUBROUTINE UPDATE_POTENTIAL(  &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,LATT_CUR, &
          T_INFO,INFO,IO, &
          GRID,GRID_SOFT,GRIDC,GRIDUS,C_TO_US,SOFT_TO_C, &
          CHTOT,DENCOR,CVTOT,CSTRF, &
          LMDIM,IRDMAX,CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,CHDEN,SV)

    USE base
    USE ini
    USE morbitalmag
    USE lattice
    USE charge
    USE pseudo
    USE lattice
    USE nonl_high
    USE mgrid
    USE poscar
    USE wave
    USE pot
    USE pawm
    USE us
    USE wave_high
    USE hamil_struct_def
    USE meta
    USE LDAPLUSU_MODULE
    IMPLICIT NONE
!=======================================================================
!  structures
!=======================================================================
    TYPE (tau_handle)  KINEDEN
    TYPE (ham_handle)  HAMILTONIAN
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    TYPE (wavedes)     WDES
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (latt)        LATT_CUR
    TYPE (info_struct) INFO
    TYPE (in_struct)   IO
    TYPE (grid_3d)     GRID       ! grid for wavefunctions
    TYPE (grid_3d)     GRID_SOFT  ! grid for soft chargedensity
    TYPE (grid_3d)     GRIDC      ! grid for potentials/charge
    TYPE (grid_3d)     GRIDUS     ! temporary grid in us.F
    TYPE (transit)     C_TO_US    ! index table between GRIDC and GRIDUS
    TYPE (transit)     SOFT_TO_C  ! index table between GRID_SOFT and GRIDC
    TYPE (energy)      E
!   local
    INTEGER LMDIM,IRDMAX

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
!  local
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

    CALL STOP_TIMING("G",IO%IU6,"POTLOK")
    DWRITE0 'potlok is ok'


    CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
         LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

    CALL SETDIJ_AVEC(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
         LMDIM,CDIJ,HAMILTONIAN%AVTOT, NONLR_S, NONL_S, IRDMAX)

    CALL SET_DD_MAGATOM(WDES, T_INFO, P, LMDIM, CDIJ)

    CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
         WDES%NCDIJ, LMDIM, CDIJ(1,1,1,1),  RHOLM, CRHODE(1,1,1,1), &
         E,  LMETA=.FALSE., LASPH=INFO%LASPH, LCOREL= .FALSE.  )

    CALL UPDATE_CMBJ(GRIDC,T_INFO,LATT_CUR,IO%IU6)

    IF (USELDApU()) CALL LDAPLUSU_PRINTOCC(WDES,T_INFO%NIONS,T_INFO%ITYP,IO%IU6)

    DWRITE0 'setdij is ok'
    CALL STOP_TIMING("G",IO%IU6,"SETDIJ")

    INFO%LPOTOK=.TRUE.

    END SUBROUTINE UPDATE_POTENTIAL

!*********************************************************************
!
!     These are values for the iteration scheme. The timesteps have to
!     be smaller than the inverse max eigenvalue, otherwise it is not
!     converging. These values are a good choice, with good speed and
!     good convergence, so best is to leave them that way. However, if
!     you are encountering problems increase the values of evmax and tmax
!     if the calculation is too time-consuming, try to decrease them.
!     (safer than to decrease tmax is to increase gamma!!)
!
!     Calculate the timesteps for the time-development scheme.
!
!*********************************************************************

    SUBROUTINE CALC_TIMESTEPS_TE(SHIFT,IO, TIMESTEP_V, NELMGW, NELMGW_ACTUAL, OMEGAMAX)
      USE base

      IMPLICIT NONE
      REAL(q)   :: SHIFT             ! complex shift determines allowed broadening
      TYPE (in_struct)     IO
      REAL(q)   :: TIMESTEP_V(:,:)
      INTEGER   :: NELMGW           ! max number of timesteps
      INTEGER   :: NELMGW_ACTUAL    ! actual number of timesteps required to reach desired convergence
      REAL(q)   :: OMEGAMAX         ! maximum transition energy
! local
      INTEGER   :: STEP
      REAL(q)   :: TMAX,TIMESTEP,SUMTIME
      REAL(q), PARAMETER :: TIMESTEP_FACTOR=2 ! 1.5  ! essentially exact with 4


      IF (IO%IU6>=0) THEN
         WRITE(IO%IU6,'(1X,(A))')  'determining timestep:'
         WRITE(IO%IU6,'(A)')      '-----------------------'
         WRITE(IO%IU6,'((A),F13.6,1X,/(A),F13.6,1X,/(A),I8)') &
                                                            ' maximum frequency       OMEGAMAX=',OMEGAMAX, &
                                                            ' complex shift           CSHIFT  =',SHIFT, &
                                                            ' maximum number of steps NELM    =',NELMGW
      ENDIF
      TMAX=5._q/SHIFT
      SUMTIME=0._q

      IF (IO%IU6>=0) THEN
         WRITE(IO%IU6,'((A),F8.4,1X,(A))') &
                                                            ' maximum time (5/CSHIFT)         =',TMAX,' 1/eV'
      ENDIF

      DO STEP=1,NELMGW
!  the relative error at OMEGMAX is given by (see below)
!  SIN(1/TIMESTEP_FACTOR)*TIMESTEP_FACTOR
!  error for 2   is about 4 %
!  error for 1.5 is about 8 %
!  error for 1   is about 16 %
!  error for 0.7 is about 30 %
!  this frequency error is however removed a posteriori
!
!  Tobias TODO: unfortunately the errors in the intensity are not entirely removed
!  test for IP
         TIMESTEP=1._q/OMEGAMAX/TIMESTEP_FACTOR

         SUMTIME=SUMTIME+TIMESTEP
         TIMESTEP_V(1,STEP)=TIMESTEP
         TIMESTEP_V(2,STEP)=SUMTIME

         IF (MOD(STEP,100)==0) THEN
            IF (IO%IU6>=0) THEN
!               WRITE(IO%IU6,'(i6,2(1x,f13.6))') step,sumtime,timestep
            ENDIF
         ENDIF

         IF (SUMTIME.GE.TMAX) THEN
            NELMGW_ACTUAL=STEP
            EXIT
         ENDIF
      ENDDO

      IF (SUMTIME.LT.TMAX) THEN
         IF (IO%IU6>=0) THEN
            WRITE(IO%IU6,'(1X,(a))') ' Too many timesteps. - I bet there is something wrong.'
            WRITE(IO%IU6,'(1X,(a))') ' Calculating anyway ...'
         ENDIF
            NELMGW_ACTUAL=NELMGW
      ENDIF

      IF (IO%IU6>=0) THEN
         write(IO%IU6,*)' '
         WRITE(IO%IU6,'(1X,(A),I5)') 'Number of time steps: ',NELMGW_ACTUAL
         WRITE(IO%IU6,'(1X,(A),F13.6,1X,(A))') 'First time step: ',TIMESTEP_V(1,1),'1/eV'
         WRITE(IO%IU6,'(1X,(A),F13.6,1X,(A))') 'Last  time step: ',TIMESTEP_V(1,NELMGW_ACTUAL),'1/eV'
         WRITE(IO%IU6,*)
      ENDIF

    END SUBROUTINE CALC_TIMESTEPS_TE

!**********************************************************************
!
! calculate DF and dump it to vasprun.xml
!
!**********************************************************************

    SUBROUTINE CALC_WRITE_DF(IO, OMEGAMAX, NEDOS, OMEGA, SHIFT, ISPIN, IEPSILON, NELMGW_ACTUAL)
    USE base
    USE vaspxml

    IMPLICIT NONE
    TYPE (in_struct)     IO              ! IO-handle
    REAL(q)             :: OMEGAMAX      ! maximum frequency
    REAL(q)             :: OMEGA         ! volume
    REAL(q)             :: SHIFT         ! complex shift = broadening
    INTEGER             :: ISPIN
    INTEGER             :: IEPSILON      ! cartesian direction (4=all)
    INTEGER             :: NELMGW_ACTUAL ! timestep
! local
    COMPLEX(q), POINTER :: EPSILON(:,:)
    REAL(q)    :: NORM
    REAL(q)    :: DOMEGA                 ! frequency step
    REAL(q)    :: TIMESTEP               !
    INTEGER    :: DIR_MIN, DIR_MAX, DIR
    INTEGER    :: NEDOS
    INTEGER    :: I,J
    REAL(q)    :: W,WEIGHT
    COMPLEX(q) :: EXPFAC,SUM_WEIGHT,MEAN_TSTEP,MEAN_TSTEP_2

    IF (IEPSILON==1) THEN
       DIR_MIN=1
       DIR_MAX=1
    ELSE IF (IEPSILON==2) THEN
       DIR_MIN=2
       DIR_MAX=2
    ELSEIF (IEPSILON==3) THEN
       DIR_MIN=3
       DIR_MAX=3
    ELSEIF (IEPSILON==4) THEN
       DIR_MIN=1
       DIR_MAX=3
    ENDIF

    DOMEGA=OMEGAMAX/(NEDOS-1)

    TIMESTEP=TIMESTEP_V(1,1)

    CELLVOL=OMEGA/AUTOA**3

    ALLOCATE(EPSILON(1:NEDOS,3))

    NORM=4._q*REAL(PI,q)*(HARTREE**3)/(CELLVOL)
    IF (ISPIN==1) NORM=2._q*NORM

    EPSILON=CZERO

!      IF (node_me==node_lead) THEN
    DO DIR=DIR_MIN,DIR_MAX
       DO I=1,NEDOS
          ! rescale frequency by error function
          W=ASIN((I-1)*DOMEGA*TIMESTEP)/TIMESTEP
          ! weight is the derivative of the function above
          WEIGHT=1/SQRT(1 - TIMESTEP**2*((I-1)*DOMEGA)**2)
          DO J=1,NELMGW_ACTUAL
             EXPFAC=EXP((ICI*W-SHIFT)*TIMESTEP_V(2,J)) !expfac=exponential-factor
             EPSILON(I,DIR)=EPSILON(I,DIR)-2.0_q*AIMAG(SCALAR(J,DIR))*EXPFAC*TIMESTEP_V(1,J)*WEIGHT
          ENDDO
          !epsilon(i,dir)=cone+norm*epsilon(i,dir)
          EPSILON(I,DIR)=NORM*GCONJG(EPSILON(I,DIR))
          IF (DIR<4) EPSILON(I,DIR)=(1.0_q,0.0_q)+EPSILON(I,DIR)
       ENDDO !ENDDO omega_sum
    ENDDO !enddo dir


    DEALLOCATE(TIMESTEP_V)

!dump dielectric function to XML file
    CALL XML_EPSILON_E_DIR(DOMEGA, EPSILON, NEDOS, IEPSILON, "time-propagation")

    DEALLOCATE(epsilon)

    END SUBROUTINE CALC_WRITE_DF


END MODULE TIME_PROPAGATION
