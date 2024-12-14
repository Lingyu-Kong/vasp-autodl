#include "symbol.inc"

MODULE mlr_main
  USE prec
  USE vaspxml
  USE lr_helper
  IMPLICIT NONE

  REAL(q),ALLOCATABLE,PRIVATE :: POSION_LAST(:,:)
  REAL(q),ALLOCATABLE,PRIVATE :: BORN_CHARGES_SAVE(:,:,:)
  REAL(q),ALLOCATABLE,PRIVATE :: DIPOLE_TIME(:,:)
CONTAINS

!*********************************************************************
!
! calculate second derivatives using linear response theory
!
! implemented by gK
! the main scheduler is based on the finite difference kernel of 
! Orest Dubay (at least the skeleton is essentially identical)
! it steps through all ions (at least those that are allowed to move)
! and calls the main linear response kerner to calculate the
! linear response of the wavefunction and the second derivatives
! with respect to the ionic positions
!
!*********************************************************************

  SUBROUTINE LR_SKELETON( &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI, &
          T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,VDW_SET,DOS,DOSI,CHAM, &
          IBRION,LMDIM,IRDMAX,NEDOS, &
          TOTENIN,EFERMI,FORCE,LFAST)

    USE ini
    USE base
    USE lattice
    USE finite_differences
    USE charge
    USE pseudo
    USE lattice
    
    USE nonl_high
    USE msymmetry
    USE mpimy
    USE mgrid
    USE mkpoints
    USE constant
    USE poscar
    USE wave
    USE pot
    USE subrot
    USE pawm
    USE rmm_diis
    USE choleski
    USE david
    USE wave_high
    USE mlrf_main
    USE lri_main
    USE subrot_cluster
    USE kpoints_change
    USE full_kpoints
    USE hamil_struct_def
    USE pead, ONLY : LUSEPEAD,LPEAD_RETURN_Q_CPROJ,PEAD_DPSI_DK_ALL,LPEAD_SYM_RED,LBORN,BORN_CHARGES_PEAD
    USE fock
    USE morbitalmag
    USE meta
    USE wannier_interpolation
! related to wave interpolation
    USE wave_interpolate
    USE mlr_optic
#ifdef VASP_HDF5
    USE vhdf5
#endif
! end realted to wave interpolation
    USE crpa, ONLY : LCRPAPLOT
    USE rpa_force, ONLY : STORE_RPA_SECOND_DERIV
    USE vdwhess
    USE locproj, ONLY: WANPROJFILE
    USE mlwf, ONLY: MLWF_ROTATE_CDER, READ_WANPROJ_FILE
    USE vdwd4, ONLY: vdw_settings
    USE tutor, ONLY: vtutor, isError, isAlert, LHFLinearResponse, MetaGGAResponse, &
        LinearResponseLREAL, LinearReponsePOSCAR
! solvation__
    USE solvation
! solvation__
#ifdef PROFILING
    USE profiling
#endif
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
    TYPE (wavespin)    W          ! wavefunction
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
    TYPE (energy)      E
    TYPE (latt)        LATT_INI
    TYPE (vdw_settings), INTENT(IN) :: VDW_SET

    INTEGER IBRION
    INTEGER LMDIM,IRDMAX,IRDMAA,NEDOS
    REAL(q) TOTENIN,EFERMI

    COMPLEX(q)  CHTOT(GRIDC%MPLWV,WDES%NCDIJ) ! charge-density in real / reciprocal space
    COMPLEX(q)  CHTOTL(GRIDC%MPLWV,WDES%NCDIJ)! old charge-density
    RGRID       DENCOR(GRIDC%RL%NP)           ! partial core
    COMPLEX(q)  CVTOT(GRIDC%MPLWV,WDES%NCDIJ) ! local potential
    COMPLEX(q)  CSTRF(GRIDC%MPLWV,T_INFO%NTYP)

!  augmentation related quantities
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
!  Hamiltonian
    GDEF       CHAM(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN)
    REAL(q)   XCSIF(3,3)
    REAL(q) :: FORCE(3,T_INFO%NIONS)   ! forces in cartesian coordinates
    LOGICAL, OPTIONAL :: LFAST


! local variables related to finite difference code

    TYPE (eigenf_cluster_pointer),POINTER :: DEG_CLUSTER(:,:)
    INTEGER :: NIONS
    INTEGER :: IU6                ! OUTCAR file
    INTEGER :: IU6K               ! OUTCAR k-points
    INTEGER :: IU0                ! stdout

    TYPE (type_info)   T_INFO_0
    INTEGER :: IUDYNMAT           ! DYNMAT file
    REAL(q) :: X
    REAL(q),ALLOCATABLE      :: INITIAL_FORCE(:,:)
    REAL(q),ALLOCATABLE      :: DISPL_FORCES(:,:,:)
    REAL(q),ALLOCATABLE      :: INT_STRAIN(:,:,:), PIEZO(:,:,:)
    REAL(q),ALLOCATABLE      :: BORN_CHARGES(:,:,:), BORN_CHARGES2(:,:,:)
    REAL(q),ALLOCATABLE      :: SECOND_DERIV(:,:)
    INTEGER                  :: DOF
    INTEGER                  :: PROCESSED_DISPL
    INTEGER                  :: I,J,K,M,N, IDIR, IDISPL_ON_ION
    INTEGER                  :: II,JJ
    REAL(q),ALLOCATABLE      :: WORK(:,:)
    REAL(q),ALLOCATABLE      :: EIGENVECTORS(:,:)
    REAL(q),ALLOCATABLE      :: EIGENVALUES(:)
    INTEGER                  :: IERROR
    REAL(q) :: EPSILON(3,3),  ELASTIC(3,3,3,3),  ELASTICP(6,3,3)
! variables for the reevaluation of energy
    REAL(q) TOTEN, DESUM      ! energy and change of one electron energies
! linear response using symmetry
    INTEGER                  :: NKORIG
    REAL(q),ALLOCATABLE      :: D(:,:,:)
    INTEGER,ALLOCATABLE      :: ND(:), IDIRD(:,:)
    REAL(q)                  :: WORKD(3,3,T_INFO%NIONS)
    INTEGER                  :: IWORK(T_INFO%NIONS)
    CHARACTER(3)             :: IDIR_TEXT(3)=(/"x","y","z"/)
    LOGICAL                  :: LDO(T_INFO%NIONS)
    REAL(q)                  :: WDMAT(3,T_INFO%NIONS,3,T_INFO%NIONS),DMAT(3,3,T_INFO%NIONS,T_INFO%NIONS)
    REAL(q)                  :: ST(3,3,3,T_INFO%NIONS), AST(3,3,3,T_INFO%NIONS)
    REAL(q)                  :: FACT, TMP(3,3)
! variables to store G [H,r] phi
    COMPLEX(qs), ALLOCATABLE :: RPHI(:,:,:,:,:)
    GDEFS, ALLOCATABLE ::  RPHI_CPROJ(:,:,:,:,:)
    TYPE (skpoints_trans)   :: KPOINTS_TRANS

! symmetry related quantities (common block)
    INTEGER ISYMOP,NROT,IGRPOP,NROTK,INVMAP,NPCELL
    REAL(q) GTRANS,AP
    COMMON /SYMM/   ISYMOP(3,3,48),NROT,IGRPOP(3,3,48),NROTK, &
         GTRANS(3,48),INVMAP(48),AP(3,3),NPCELL

    ! NULLIFY to make sure a later "IF(ASSOCIATED())" works correctly.
    NULLIFY(KPOINTS_TRANS%CPHASE)

    IF (LHFCALC .AND. AEXX/=0 ) THEN
       CALL vtutor%write(isError, LHFLinearResponse)
    ENDIF

    IF (LDO_METAGGA()) THEN
       CALL vtutor%write(isError, MetaGGAResponse)
    ENDIF

    IF (INFO%LREAL) THEN
       CALL vtutor%write(isAlert, LinearResponseLREAL)
    ENDIF

    NIONS = T_INFO%NIONS
    IU6   = IO%IU6
    IU0   = IO%IU0

    IF (IO%NWRITE>=3) THEN
       IU6K  = IO%IU6
    ELSE
       IU6K  = -1
    ENDIF

    IF (IBRION>0) THEN
       IF (IBRION==8) THEN
          ALLOCATE(D(3,3,T_INFO%NIONS), ND(T_INFO%NIONS),IDIRD(3,T_INFO%NIONS))

          CALL FREDOM(SYMM%ROTMAP,ISYMOP,INVMAP,NROTK,NPCELL,D,ND,1,T_INFO%NTYP,T_INFO%NIONS, &
             T_INFO%NITYP,WORKD,IWORK, & 
             LATT_CUR%A(1,1),LATT_CUR%A(1,2),LATT_CUR%A(1,3), &
             LATT_CUR%B(1,1),LATT_CUR%B(1,2),LATT_CUR%B(1,3),IDIRD,0,[INTEGER::])
          DOF=SUM(ND)

          IF (IU6>=0) THEN
             WRITE(IU6,*)
             WRITE(IU6,'(A,I5,A)') ' Found ',DOF,' degrees of freedom:'
             WRITE(IU6,*)' ----------------------------------------------'
             WRITE(IU6,*)
             DO J=1,NIONS
                IF (ND(J).GT.0) THEN
                   WRITE(IU6,'(A,I5,A,3A2)') '     directions for atom ',J,':  ',IDIR_TEXT(IDIRD(ND(J),J))
                ENDIF
             ENDDO
          ENDIF
       ELSE IF (IBRION==7) THEN
          CALL COUNT_DOF(NIONS, T_INFO%LSFOR, T_INFO%LSDYN, DOF)
       ELSE
          DOF=0
       ENDIF
    ELSE
       DOF=0
    ENDIF

    T_INFO_0=T_INFO
    NULLIFY(T_INFO_0%POSION)
    ALLOCATE(T_INFO_0%POSION(3,NIONS))
    ALLOCATE(INITIAL_FORCE(3,NIONS))
    ALLOCATE(DISPL_FORCES(DOF,3,NIONS),BORN_CHARGES(3,3,NIONS),BORN_CHARGES2(3,3,NIONS), & 
      PIEZO(3,3,3),INT_STRAIN(DOF,3,3))

    DISPL_FORCES=0
    BORN_CHARGES=0 ; BORN_CHARGES2=0
    PIEZO=0
    INT_STRAIN=0

    T_INFO_0%POSION             = T_INFO%POSION
    INITIAL_FORCE               = FORCE
    NKORIG                      = WDES%NKPTS
!=======================================================================
! reset the potential and recalculate the ground state wavefunctions
! with very high precision for this setup
!=======================================================================
    IF (IU0>=0) THEN
       WRITE (IU0,*) 'Linear response reoptimize wavefunctions to high precision'
       WRITE ( 17,*) 'Linear response reoptimize wavefunctions'
    ENDIF

    IF (IU6>=0) THEN
       WRITE (IU6,*) 'Linear response reoptimize wavefunctions'
    ENDIF

    IF (INFO%LREAL) THEN
       CALL RSPHER(GRID,NONLR_S,LATT_CUR)
    ENDIF

    CALL POTLOK(GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES,  &
         INFO,P,T_INFO,E,LATT_CUR, &
         CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF)

!   CALL VECTORPOT(GRID, GRIDC, GRID_SOFT, SOFT_TO_C,  WDES%COMM_INTER, & 
!        LATT_CUR, T_INFO%POSION, HAMILTONIAN%AVEC, HAMILTONIAN%AVTOT)
                  
    CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
         LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

    CALL SETDIJ_AVEC(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
         LMDIM,CDIJ,HAMILTONIAN%AVTOT, NONLR_S, NONL_S, IRDMAX)
    
    CALL SET_DD_MAGATOM(WDES, T_INFO, P, LMDIM, CDIJ) 
    CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
         WDES%NCDIJ, LMDIM, CDIJ(1,1,1,1),  RHOLM, CRHODE(1,1,1,1), &
         E,  LMETA=.FALSE., LASPH=INFO%LASPH, LCOREL= .FALSE.  )

    CALL RE_OPTIMIZE_ORBITALS
    CALL KPAR_SYNC_ALL(WDES,W)

    IF (IO%IU0>=0 .AND. IO%LOPEN) CALL WFORCE(17)

!=======================================================================
! determine G [H, r] |phi> = d/ dk | phi(k)>
!=======================================================================
    NULLIFY(DEG_CLUSTER)
    CALL FIND_DEG_CLUSTERS(WDES, W, DEG_CLUSTER)

    IF (LEPSILON .OR. LMAGBLOCH ) THEN

    CALL SET_NABIJ_AUG(P,T_INFO%NTYP)

    ALLOCATE(RPHI(WDES%NRPLWV,WDES%NBANDS,WDES%NKPTS,WDES%ISPIN,3), &
             RPHI_CPROJ(WDES%NPROD,WDES%NBANDS,WDES%NKPTS,WDES%ISPIN,3))

    IF (LUSEPEAD()) THEN
       IF (LMAGBLOCH) LPEAD_RETURN_Q_CPROJ=.FALSE.

       RPHI=0
       RPHI_CPROJ=0

       CALL PEAD_DPSI_DK_ALL(W,KPOINTS,P,CQIJ,LATT_CUR,T_INFO,RPHI,RPHI_CPROJ)

       LPEAD_RETURN_Q_CPROJ=.TRUE.
    ELSE 
       RPHI=0
       RPHI_CPROJ=0

       DO IDIR=1,3
          IF (IU0>=0) THEN
             WRITE (IU0,*)'Linear response G [H, r] |phi>, progress :'
             WRITE (IU0,'(A,I3)') &
                  '  Direction: ',IDIR
             WRITE (17,*)'Linear response G [H, r] |phi>, progress :'
             WRITE (17,'(A,I3)') &
                  '  Direction: ',IDIR
          END IF

          IF (IU6>=0) THEN
             WRITE (IU6,*)'Linear response G [H, r] |phi>, progress :'
             WRITE (IU6,'(A,I3)') &
                  '  Direction: ',IDIR
          ENDIF

          CALL LRF_RPHI( &
             P,WDES,NONLR_S,NONL_S,W,LATT_CUR, &
             T_INFO,INFO,IO,GRID,GRIDC,GRIDUS,C_TO_US,IRDMAX, &
             CDIJ,CQIJ,SV,LMDIM,DEG_CLUSTER, IDIR, RPHI(:,:,:,:,IDIR), RPHI_CPROJ(:,:,:,:,IDIR), & 
             .TRUE.)

          IF (IO%LOPEN) CALL WFORCE(IO%IU6)
       ENDDO
    ENDIF

    IF (LMAGBLOCH) THEN

       CALL BLOCH_CURRENT( W, GRID_SOFT, GRIDC, GRIDUS, C_TO_US, SOFT_TO_C, P, LATT_CUR, & 
          HAMILTONIAN%AVEC, HAMILTONIAN%AVTOT, CHTOT, NONLR_S, NONL_S, & 
          RPHI, RPHI_CPROJ, CDIJ, CQIJ, SV, EFERMI, &
          T_INFO, LMDIM, CRHODE, IRDMAX, IO%IU6, IO%IU0)

       ! well this deallocates most of the quantities we have used, but not all
       DEALLOCATE(RPHI, RPHI_CPROJ)
       RETURN
    ENDIF

!=======================================================================
! response to external excluding local field effects
! i.e. in the independent particle approximation
!=======================================================================
    IF (.NOT. PRESENT (LFAST)) THEN
    DO IDIR=1,3

       IF (IU0>=0) THEN
          WRITE (IU0,*)'Linear response to external field (no local field effect), progress :'
          WRITE (IU0,'(A,I3)') &
               '  Direction: ',IDIR
          WRITE (17,*)'Linear response to external field (no local field effect), progress :'
          WRITE (17,'(A,I3)') &
               '  Direction: ',IDIR
       END IF
     
       IF (IU6>=0) THEN
          WRITE (IU6,*)'Linear response to external field (no local field effect), progress :'
          WRITE (IU6,'(A,I3)') &
               '  Direction: ',IDIR
       ENDIF

       IF (SYMM%ISYM>0) THEN
          DYN%VEL=0
          DYN%VEL(IDIR,:)=1.0
          CALL KARDIR(T_INFO%NIONS,DYN%VEL,LATT_CUR%B)
          CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
               T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
               SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
               SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)
       ENDIF

       CALL LRF_MAIN( &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,RPHI,RPHI_CPROJ, &
          LATT_CUR, &
          T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,VDW_SET,DOS,DOSI, &
          LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI, DEG_CLUSTER, KPOINTS_TRANS, EPSILON(IDIR,:), BORN_CHARGES(IDIR,:,:), PIEZO(IDIR,:,:),  & 
          IDIR, .FALSE.)

       IF (IO%LOPEN) CALL WFORCE(IO%IU6)
    ENDDO
    ! reinitialise symmetry
    IF (SYMM%ISYM>0) THEN
       DYN%VEL=0
       CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
            T_INFO%LSDYN,.FALSE.,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
            SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
            SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)

       CALL TSYM(EPSILON,ISYMOP,NROTK,LATT_CUR%A)

       IF (LPEAD_SYM_RED()) &
          CALL QSTAR_SYM_REDUCE(REAL(RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/)),KIND=q),T_INFO,SYMM,IO%IU6)
    ENDIF

    IF (IO%IU6>=0) THEN
       WRITE(IO%IU6,1110) '(INDEPENDENT PARTICLE, excluding Hartree and local field effects)',EPSILON
       WRITE(IO%IU6,130)
       CALL XML_TENSOR("dielectric_ipa",EPSILON)
#ifdef VASP_HDF5
       CALL VH5_WRITE_DIELECTRIC_STATIC(IH5OUTFILEID,"independent_particle",EPSILON)
#endif
    ENDIF
    ENDIF
!=======================================================================
! response to external field including local field effects
! if do_kpoints_chage is not set the k-point set remains unchanged
! and symmetrization of the final tensor is performed
! results are not exact but errors are sometimes small (not recommended)
!=======================================================================
#define do_kpoints_change
#ifndef do_kpoints_change
    DO IDIR=1,3

       IF (IU0>=0) THEN
          WRITE (IU0,*)'Linear response to external field, progress :'
          WRITE (IU0,'(A,I3)') &
               '  Direction: ',IDIR
          WRITE (17,*)'Linear response to external field, progress :'
          WRITE (17,'(A,I3)') &
               '  Direction: ',IDIR
       END IF
     
       IF (IU6>=0) THEN
          WRITE (IU6,*)'Linear response to external field, progress :'
          WRITE (IU6,'(A,I3)') &
               '  Direction: ',IDIR
       ENDIF

       IF (SYMM%ISYM>0) THEN
          DYN%VEL=0
          DYN%VEL(IDIR,:)=1.0
          CALL KARDIR(T_INFO%NIONS,DYN%VEL,LATT_CUR%B)
          CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
               T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
               SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
               SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)
       ENDIF

       CALL LRF_MAIN( &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,RPHI,RPHI_CPROJ, &
          LATT_CUR, &
          T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,DOS,DOSI, &
          LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI, DEG_CLUSTER, KPOINTS_TRANS, EPSILON(IDIR,:), BORN_CHARGES(IDIR,:,:), PIEZO(IDIR,:,:),  & 
          IDIR, .TRUE.)

       IF (IO%LOPEN) CALL WFORCE(IO%IU6)
    ENDDO
    ! reinitialise symmetry
    IF (SYMM%ISYM>0) THEN
       DYN%VEL=0
       CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
            T_INFO%LSDYN,.FALSE.,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
            SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
            SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)

       CALL TSYM(EPSILON,ISYMOP,NROTK,LATT_CUR%A)

       IF (LPEAD_SYM_RED()) &
          CALL QSTAR_SYM_REDUCE(REAL(RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/)),KIND=q),T_INFO,SYMM,IO%IU6)
    ENDIF

    IF (IO%IU6>=0) THEN
       IF (LRPA) THEN
          WRITE(IO%IU6,1100) '(including local field effects in RPA (Hartree))',EPSILON
          CALL XML_TENSOR("dielectric_rpa",EPSILON)
#ifdef VASP_HDF5
       CALL VH5_WRITE_DIELECTRIC_STATIC(IH5OUTFILEID,"electron",EPSILON,"rpa")
#endif
       ELSE
          WRITE(IO%IU6,1100) '(including local field effects in DFT)',EPSILON
          CALL XML_TENSOR("dielectric_dft",EPSILON)
#ifdef VASP_HDF5
       CALL VH5_WRITE_DIELECTRIC_STATIC(IH5OUTFILEID,"electron",EPSILON,"dft")
#endif
       ENDIF
       WRITE(IO%IU6,130)
       ! The tags above are more informative! We keep 'epsilon' for legacy but will discontinue it
       CALL XML_TENSOR("epsilon",EPSILON)
    ENDIF

!=======================================================================
! response to external field including local field effects
! this version changes the k-point set on the fly
! and symmetrizes the polarisation vector
!=======================================================================
#else
    DO IDIR=1,3
       IF (IU0>=0) THEN
          WRITE (IU0,*)'Linear response to external field, progress :'
          WRITE (IU0,'(A,I3)') &
               '  Direction: ',IDIR
          WRITE (17,*)'Linear response to external field, progress :'
          WRITE (17,'(A,I3)') &
               '  Direction: ',IDIR
       END IF
     
       IF (IU6>=0) THEN
          WRITE (IU6,*)'Linear response to external field, progress :'
          WRITE (IU6,'(A,I3)') &
               '  Direction: ',IDIR
       ENDIF
       !
       ! reinitialise symmetry part for field in direction IDIR
       ! presently this is done by supplying a velocity field to the ions
       !
       IF (SYMM%ISYM>0) THEN
          DYN%VEL=0
          DYN%VEL(IDIR,:)=1.0
          CALL  KARDIR(T_INFO%NIONS,DYN%VEL,LATT_CUR%B)
          CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
               T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
               SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
               SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)
#ifdef oldsym
          CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
               SYMM%ISYM>=0.AND..NOT.WDES%LSORBIT.AND..NOT.WDES%LSPIRAL, &
               T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IU6K,IO%IU0)
#else
          CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
               SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
               T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IU6K,IO%IU0)
#endif
          CALL KPAR_SYNC_ALL(WDES,W)
          CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
          CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR, KPOINTS_TRANS)

          ! Loewdin perturbation theory to improve states at added k-points
          CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
          &    LMDIM,CDIJ,CQIJ,4,SV,T_INFO,P,IO%IU0,DESUM, NKSTART=NKORIG+1)
          CALL KPAR_SYNC_ALL(WDES,W)

          CALL FIND_DEG_CLUSTERS(WDES, W, DEG_CLUSTER)
       ENDIF

       CALL LRF_MAIN( &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,RPHI,RPHI_CPROJ, &
          LATT_CUR, &
          T_INFO,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,VDW_SET,DOS,DOSI, &
          LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI, DEG_CLUSTER, KPOINTS_TRANS, EPSILON(IDIR,:), BORN_CHARGES(IDIR,:,:), PIEZO(IDIR,:,:), & 
          IDIR, .TRUE.)
       IF (IO%LOPEN) CALL WFORCE(IO%IU6)

       IF (SYMM%ISYM>0) CALL DEALLOCATE_KPOINTS_TRANS(KPOINTS_TRANS)
    ENDDO

    IF (SYMM%ISYM>0) THEN
       DYN%VEL=0
       CALL KARDIR(T_INFO%NIONS,DYN%VEL,LATT_CUR%B)
       CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
            T_INFO%LSDYN,.FALSE.,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
            SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
            SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)

       CALL TSYM(EPSILON,ISYMOP,NROTK,LATT_CUR%A)

       IF (LPEAD_SYM_RED()) &
          CALL QSTAR_SYM_REDUCE(REAL(RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/)),KIND=q),T_INFO,SYMM,IO%IU6)
#ifdef oldsym
       CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
            SYMM%ISYM>=0.AND..NOT.WDES%LSORBIT.AND..NOT.WDES%LSPIRAL, &
            T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IU6K,IO%IU0)
#else
       CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
            SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
            T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IU6K,IO%IU0)
#endif
       CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
       CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
    ENDIF

    IF (IO%IU6>=0) THEN
       IF (LRPA) THEN
          WRITE(IO%IU6,1100) '(including local field effects in RPA (Hartree))',EPSILON
          CALL XML_TENSOR("dielectric_rpa",EPSILON)
#ifdef VASP_HDF5
          CALL VH5_WRITE_DIELECTRIC_STATIC(IH5OUTFILEID,"electron",EPSILON,"rpa")
#endif
       ELSE
          WRITE(IO%IU6,1100) '(including local field effects in DFT)',EPSILON
          CALL XML_TENSOR("dielectric_dft",EPSILON)
#ifdef VASP_HDF5
          CALL VH5_WRITE_DIELECTRIC_STATIC(IH5OUTFILEID,"electron",EPSILON,"dft")
#endif
       ENDIF
       WRITE(IO%IU6,130)
       CALL XML_TENSOR("epsilon",EPSILON)
    ENDIF
#endif
1100   FORMAT(// &
            " MACROSCOPIC STATIC DIELECTRIC TENSOR ",A/, &
            " ------------------------------------------------------"/, &
            &       3(6X,3F13.6/), &
            " ------------------------------------------------------"/)

1110   FORMAT(// &
            " HEAD OF MICROSCOPIC STATIC DIELECTRIC TENSOR ",A/, &
            " ------------------------------------------------------"/, &
            &       3(6X,3F13.6/), &
            " ------------------------------------------------------"/)


    ENDIF
130 FORMAT (5X, //, &
     &'----------------------------------------------------', &
     &'----------------------------------------------------'//)
!=======================================================================
! ionic displacements
!=======================================================================

    IF (IU0>=0 .AND. DOF>0) THEN
       WRITE (IU0,*) 'Linear response DOF=',DOF
       WRITE (17,*) 'Linear response DOF=',DOF
    ENDIF

    IF (IU6>=0 .AND. DOF>0) THEN
       WRITE (IU6,*) 'Linear response:'
       WRITE (IU6,*) '  Degrees of freedom DOF   = ',DOF
    END IF

    DO PROCESSED_DISPL=1,DOF

       IF (IU0>=0) THEN
          WRITE (IU0,*)'Linear response progress:'
          WRITE (IU0,'(A,I3,A,I3)') &
               '  Degree of freedom: ',PROCESSED_DISPL,'/',DOF
          WRITE (17,*)'Linear response progress:'
          WRITE (17,'(A,I3,A,I3)') &
               '  Degree of freedom: ',PROCESSED_DISPL,'/',DOF
       END IF
     
       IF (IU6>=0) THEN
          WRITE (IU6,*)'Linear response progress:'
          WRITE (IU6,'(A,I3,A,I3)') &
               '  Degree of freedom: ',PROCESSED_DISPL,'/',DOF
       ENDIF
       T_INFO%POSION=T_INFO_0%POSION
       IF (IBRION==8) THEN
          CALL FIND_IJ_ID(.FALSE.,NIONS,PROCESSED_DISPL,ND,J,IDISPL_ON_ION)
          IDIR=IDIRD(IDISPL_ON_ION,J)
          IF (ABS(D(IDIR,IDISPL_ON_ION,J)-1.0) >1E-8) THEN
             CALL vtutor%write(isError, LinearReponsePOSCAR)
          ENDIF
       ELSE
          CALL FIND_IJ(NIONS,T_INFO%LSFOR,T_INFO%LSDYN,PROCESSED_DISPL,IDIR,J)
       ENDIF

       IF (SYMM%ISYM>0) THEN
          DYN%VEL=0
          DYN%VEL(IDIR,J)=1.0
          CALL  KARDIR(1,DYN%VEL(:,J),LATT_CUR%B)
          CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
               T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
               SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
               SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)
#ifdef oldsym
          CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
               SYMM%ISYM>=0.AND..NOT.WDES%LSORBIT.AND..NOT.WDES%LSPIRAL, &
               T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IU6K,IO%IU0)
#else
          CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
               SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
               T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IU6K,IO%IU0)
#endif
          CALL KPAR_SYNC_ALL(WDES,W)
          CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
          CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR, KPOINTS_TRANS)

          ! Loewdin perturbation theory to improve states at added k-points
          CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W,WDES,SYMM, &
               &    LMDIM,CDIJ,CQIJ,4,SV,T_INFO,P,IO%IU0,DESUM, NKSTART=NKORIG+1)
          
          ! this call causes issues in the piezoelectric tensors
          ! most likely because orbitals are rotated
          ! also does not really improve results
!          CALL RE_OPTIMIZE_ORBITALS

          CALL KPAR_SYNC_ALL(WDES,W)
          CALL FIND_DEG_CLUSTERS(WDES, W, DEG_CLUSTER)
       ENDIF

       CALL LR_MAIN( &
          KINEDEN,HAMILTONIAN,P,WDES,NONLR_S,NONL_S,W,LATT_CUR, &
          T_INFO,T_INFO_0,DYN,INFO,IO,MIX,KPOINTS,SYMM,GRID,GRID_SOFT, &
          GRIDC,GRIDB,GRIDUS,C_TO_US,B_TO_C,SOFT_TO_C,E, &
          CHTOT,CHTOTL,DENCOR,CVTOT,CSTRF, &
          CDIJ,CQIJ,CRHODE,N_MIX_PAW,RHOLM,RHOLM_LAST, &
          CHDEN,SV,VDW_SET,DOS,DOSI, &
          LMDIM,IRDMAX,NEDOS, &
          TOTEN,EFERMI,DISPL_FORCES(PROCESSED_DISPL,:,:),INT_STRAIN(PROCESSED_DISPL,:,:), & 
          IDIR, J, BORN_CHARGES2(:,IDIR,J), &
          DEG_CLUSTER, KPOINTS_TRANS, RPHI, RPHI_CPROJ)
       IF (IO%LOPEN) CALL WFORCE(IO%IU6)

       IF (SYMM%ISYM>0) CALL DEALLOCATE_KPOINTS_TRANS(KPOINTS_TRANS)
    ENDDO

    CALL FREE_DEG_CLUSTERS(WDES,DEG_CLUSTER)

    IF (LEPSILON) THEN
       DEALLOCATE(RPHI, RPHI_CPROJ)
    ENDIF

    ! restore original symmetry
    IF (SYMM%ISYM>0) THEN
       DYN%VEL=0
       CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
            T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
            SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
            SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)
#ifdef oldsym
       CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
            SYMM%ISYM>=0.AND..NOT.WDES%LSORBIT.AND..NOT.WDES%LSPIRAL, &
            T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IU6K,IO%IU0)
#else
       CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
            SYMM%ISYM>=0.AND..NOT.WDES%LNONCOLLINEAR, &
            T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IU6K,IO%IU0)
#endif
       CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
       CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
    ENDIF

!=======================================================================
!
! final processing and output 
! Born effective charges, piezoelectric tensors, and internal strain
!
!=======================================================================
    IF (SYMM%ISYM>0.AND.LPEAD_SYM_RED()) THEN
       ! switch to full symmetry, including symmetry operations that do not
       ! conserve the finite difference stencils of the PEAD routines
       DYN%VEL=0
       CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
            T_INFO%LSDYN,.FALSE.,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
            SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
            SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IU6=-1)
    ENDIF

    IF (IBRION==8) THEN
       !
       ! version for case symmetry was used to reduce the number of displacements
       !
       CALL PRINT_DYNMAT_ID(.FALSE.,NIONS,DOF,1.0_q,T_INFO%NTYP,T_INFO%NITYP,T_INFO%POMASS,DISPL_FORCES,D,ND,IU6)

       DO PROCESSED_DISPL=1,DOF
          CALL FIND_IJ_ID(.FALSE.,NIONS,PROCESSED_DISPL,ND,J,IDIR)
          DO N=1,T_INFO%NIONS
             DMAT(1:3,IDIR,N,J)=DISPL_FORCES(PROCESSED_DISPL,1:3,N)
          END DO
          ST(1:3,1:3,IDIR,J)=INT_STRAIN  (PROCESSED_DISPL,1:3,1:3)
       END DO

       CALL MKDMAT(SYMM%ROTMAP,ISYMOP,INVMAP,NROTK,NPCELL,D,DMAT,ST,ND,      &
            1,T_INFO%NTYP,T_INFO%NIONS,T_INFO%NITYP,WORKD,WDMAT,AST,IWORK,SYMM%TAU,SYMM%TAUROT,   &
            SYMM%WRKROT, &
            LATT_CUR%A(1,1),LATT_CUR%A(1,2),LATT_CUR%A(1,3), &
            LATT_CUR%B(1,1),LATT_CUR%B(1,2),LATT_CUR%B(1,3))

       LDO=.TRUE.
       CALL SYDMAT(DMAT,SYMM%ROTMAP,ISYMOP,NROTK,NPCELL,1,T_INFO%NTYP, &
            T_INFO%NIONS,T_INFO%NITYP,WDMAT,LATT_CUR%A(1,1),LATT_CUR%A(1,2),LATT_CUR%A(1,3),LDO)
       CALL STMAT(ST  ,SYMM%ROTMAP,ISYMOP,NROTK,NPCELL,1,T_INFO%NTYP, &
            T_INFO%NIONS,T_INFO%NITYP,AST  ,LATT_CUR%A(1,1),LATT_CUR%A(1,2),LATT_CUR%A(1,3),LDO)

       ! now reset DOF to NIONS*3
       ! and copy results from DMAT back to DISPL_FORCES

       CALL COUNT_DOF(NIONS, T_INFO%LSFOR, T_INFO%LSDYN, DOF)
       DEALLOCATE(DISPL_FORCES)
       ALLOCATE(DISPL_FORCES(DOF,3,NIONS))

       DO PROCESSED_DISPL=1,DOF
          CALL FIND_IJ(NIONS,T_INFO%LSFOR,T_INFO%LSDYN,PROCESSED_DISPL,IDIR,J)
          DO N=1,NIONS
             DISPL_FORCES(PROCESSED_DISPL,1:3,N)=DMAT(1:3,IDIR,N,J)
          ENDDO
       ENDDO
    ELSE
       !
       ! version for case no symmetry was used
       !
       CALL PRINT_DYNMAT(NIONS,DOF,1.0_q,T_INFO%NTYP,T_INFO%NITYP,T_INFO%POMASS,DISPL_FORCES,T_INFO%LSDYN,T_INFO%LSFOR,IU6)

       DO PROCESSED_DISPL=1,DOF
          CALL FIND_IJ(NIONS,T_INFO%LSFOR,T_INFO%LSDYN,PROCESSED_DISPL,IDIR,J)
          ST(1:3,1:3,IDIR,J)=INT_STRAIN  (PROCESSED_DISPL,1:3,1:3)
       END DO
    ENDIF

    IF (IU6>=0) THEN
       WRITE (IU6,130) 
       IF (LEPSILON) THEN
          ! for the sake of having everything at the very end print EPSILON again
          IF (LRPA) THEN
             WRITE(IO%IU6,1100) '(including local field effects in RPA (Hartree))',EPSILON
          ELSE
             WRITE(IO%IU6,1100) '(including local field effects in DFT)',EPSILON
          ENDIF

          FACT=EVTOJ*1E20_q/LATT_CUR%OMEGA
          
170       FORMAT(/ ' ',A,'  for field in x, y, z        (e  Angst)',/ &
               10X,'XX', 10X,'YY', 10X,'ZZ',10X,'XY', 10X,'YZ', 10X,'ZX'/ &
               '  ----------------------------------------------------', &
               '----------------------------')

180       FORMAT(/ ' ',A,'  for field in x, y, z        (C/m^2)',/ &
               &        10X,'XX', 10X,'YY', 10X,'ZZ', 10X,'XY', 10X,'YZ', 10X,'ZX'/ &
               &        '  ----------------------------------------------------', &
               &        '----------------------------')
          IF (.NOT.LRPA) THEN
             CALL TSYM3(PIEZO,ISYMOP,NROTK,LATT_CUR%A)

             WRITE (IU6,170) 'PIEZOELECTRIC TENSOR (including local field effects)'
             DO I =1,3
                WRITE (IU6,140) IDIR_TEXT(I),(PIEZO(I,J,J),J=1,3), & 
                     PIEZO(I,1,2),PIEZO(I,2,3),PIEZO(I,3,1)
             ENDDO

             WRITE (IU6,180) 'PIEZOELECTRIC TENSOR (including local field effects)'
             DO I =1,3
                WRITE (IU6,140) IDIR_TEXT(I),(PIEZO(I,J,J)*FACT,J=1,3), & 
                     PIEZO(I,1,2)*FACT,PIEZO(I,2,3)*FACT,PIEZO(I,3,1)*FACT
             ENDDO
140          FORMAT(2X,A1,6F12.5)
#ifdef VASP_HDF5
             CALL VH5_WRITE_PIEZOELECTRIC_TENSOR(IH5OUTFILEID, "electron", PIEZO * FACT)
#endif

             ! remove any drift from BORN charges
             IF (IO%IU0>0) WRITE(IO%IU0,*) 'Born effective charges drift removed'
             TMP=0
             DO N=1,T_INFO%NIONS
                TMP(:,:)=TMP(:,:)+BORN_CHARGES(:,:,N)/T_INFO%NIONS
             ENDDO
             DO N=1,T_INFO%NIONS
                BORN_CHARGES(:,:,N)=BORN_CHARGES(:,:,N)-TMP(:,:)
             ENDDO

             CALL TENSYM(BORN_CHARGES,SYMM%ROTMAP(1,1,1),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIOND,LATT_CUR%A)
 
             WRITE (IU6,*)
             WRITE (IU6,'(X,A)') 'BORN EFFECTIVE CHARGES (including local field effects) (in |e|, cummulative output)'
             WRITE (IU6,'(X,A)') '---------------------------------------------------------------------------------'
             
             DO N=1,T_INFO%NIONS
                WRITE (IU6,'(" ion ",I4)') N
                DO IDIR =1,3
                   WRITE (IU6,'(I5,3F12.5)') IDIR,BORN_CHARGES(IDIR,:,N)
                ENDDO
             ENDDO

             CALL XML_BORN_CHARGES(BORN_CHARGES,T_INFO%NIONS)
#ifdef VASP_HDF5
             CALL VH5_WRITE_BORN_CHARGES(IH5OUTFILEID,BORN_CHARGES)
#endif

            ! now set the Born charges in the global array defined in pead
             IF (.NOT. LBORN) THEN
                LBORN=.TRUE.
                IF (ALLOCATED(BORN_CHARGES_PEAD)) DEALLOCATE(BORN_CHARGES_PEAD)
                ALLOCATE(BORN_CHARGES_PEAD(3,3,T_INFO%NIOND))
             ENDIF
             BORN_CHARGES_PEAD=BORN_CHARGES
          ENDIF
       ENDIF
    ENDIF
    IF (DOF>0 .AND. IU6>=0) THEN
       IF (DYN%ISIF>0) THEN
160       FORMAT(/ ' INTERNAL STRAIN TENSOR FOR ION ',I4,' for displacements in x,y,z  (eV/Angst):',/ &
               10X,'X', 11X,'Y', 11X,'Z', 10X,'XY', 10X,'YZ', 10X,'ZX'/ &
               '  ----------------------------------------------------', &
               '----------------------------')
          DO N=1,T_INFO%NIONS
             WRITE(IU6,160) N
             DO I =1,3
                WRITE (IU6,140) IDIR_TEXT(I),(ST(J,J,I,N),J=1,3),ST(1,2,I,N),ST(2,3,I,N),ST(3,1,I,N)
             ENDDO
          ENDDO
       END IF
       WRITE (IU6,130)
#ifdef VASP_HDF5
       CALL VH5_WRITE_INTERNAL_STRAIN(IH5OUTFILEID,ST)
#endif
    ENDIF

!=======================================================================
! apply delta impulse in time DFIELD on ions, by adding 
! a force term to FORCE
!=======================================================================
    IF (.NOT. DFIELD_DONE) THEN
       DFIELD_DONE=.TRUE.
       DO N=1,T_INFO%NIONS
          DO IDIR =1,3 ! loop over cartesian direction of external field
             FORCE(:,N)=FORCE(:,N)+BORN_CHARGES(IDIR,:,N)*DFIELD(IDIR)/DYN%POTIM
          ENDDO
       ENDDO
       
       ALLOCATE(POSION_LAST(3,T_INFO%NIONS))
       ALLOCATE(BORN_CHARGES_SAVE(3,3,T_INFO%NIONS))

       ! store initial postions for later use (required for evolution of dipole moment)
       POSION_LAST=T_INFO_0%POSION
    ENDIF
    
    IF (ALLOCATED(BORN_CHARGES_SAVE) .AND. LEPSILON) THEN
       BORN_CHARGES_SAVE=BORN_CHARGES
    ENDIF
!=======================================================================
!
! final processing and output
! vibrational frequencies etc.
!
!=======================================================================
    IF (DOF>0 .AND. IU6>=0 ) THEN

       ALLOCATE(SECOND_DERIV(DOF,DOF))
       SECOND_DERIV=0._q      
 
       DO N=1,DOF
          CALL FIND_IJ(NIONS,T_INFO%LSFOR,T_INFO%LSDYN,N,I,J)       
          DO M=1,DOF
             SECOND_DERIV(M,N)=DISPL_FORCES(M,I,J)
             !c add analytic vdw contribution if available
             IF (LHESSEVDW) THEN
               CALL FIND_IJ(NIONS,T_INFO%LSFOR,T_INFO%LSDYN,M,II,JJ)
               SECOND_DERIV(M,N)=SECOND_DERIV(M,N)  + HESSE_VDW(3*j-3+i,3*jj-3+ii)
               !SECOND_DERIV(M,N)= HESSE_VDW(3*j-3+i,3*jj-3+ii)
             ENDIF
          END DO
       END DO
       
       IF (IU6>=0) THEN
          WRITE (IU6,*) 
          WRITE (IU6,*) 'SECOND DERIVATIVES (NOT SYMMETRIZED)'
          WRITE (IU6,*) '------------------------------------'
          CALL PRINT_SECOND_DERIV(NIONS,DOF,SECOND_DERIV,T_INFO%LSFOR,T_INFO%LSDYN,IU6)
       END IF

       ! symmetrize Hessian matrix
       DO N=1,DOF
          DO M=N+1,DOF
             X=0.5_q*(SECOND_DERIV(N,M)+SECOND_DERIV(M,N))
             SECOND_DERIV(N,M)=X
             SECOND_DERIV(M,N)=X
          END DO
       END DO

#ifdef VASP_HDF5
      CALL VH5_WRITE_FORCE_CONSTANTS(IH5OUTFILEID,SECOND_DERIV)
#endif

       ALLOCATE(WORK(DOF,32),EIGENVECTORS(DOF,DOF),EIGENVALUES(DOF))

       ! Compute Hessian (mass normalized force-constants)
       EIGENVECTORS=SECOND_DERIV
       N=1
       DO I=1,T_INFO%NTYP
          DO J=1,T_INFO%NITYP(I)
             DO K=1,3
                CALL FIND_DOF_INDEX(NIONS,T_INFO%LSFOR,T_INFO%LSDYN,K,N,M)
                IF (M>0) EIGENVECTORS(:,M)=EIGENVECTORS(:,M)/SQRT(T_INFO%POMASS(I))
                IF (M>0) EIGENVECTORS(M,:)=EIGENVECTORS(M,:)/SQRT(T_INFO%POMASS(I))
             END DO
             N=N+1
          END DO
       END DO

#ifdef VASP_HDF5
      CALL VH5_WRITE_HESSIAN(IH5OUTFILEID,EIGENVECTORS)
#endif

       CALL XML_TAG("dynmat")
       CALL XML_TAG_STRING("unit", "THz^2")

       CALL XML_VECARRAY("hessian")
       CALL XML_ARRAY_REAL(CONVERT_FREQUENCY(EIGENVECTORS, UNIT_INTERNAL, UNIT_HZ, 2._q) * 1E-24_q)
       CALL XML_CLOSE_TAG

       CALL DSYEV &
            ('V','U',DOF,EIGENVECTORS,DOF, &
            EIGENVALUES,WORK,32*DOF, IERROR)
       IF (IERROR/=0) THEN
          IF (IU6>=0) THEN
             WRITE(IU6,*) "Error while diagonalisation DSYEV INFO=",IERROR
             WRITE(IU6,*) "Some of (or all) eigenvectors and eigenvalues are not correct !"
          END IF
       END IF

       CALL XML_VEC_REAL(CONVERT_FREQUENCY(EIGENVALUES, UNIT_INTERNAL, UNIT_HZ, 2._q) * 1E-24_q, "eigenvalues", '(ES17.8E3)')
       CALL XML_VECARRAY("eigenvectors")
       CALL XML_ARRAY_REAL(EIGENVECTORS,'(ES17.8E3)')
       CALL XML_CLOSE_TAG
       CALL XML_CLOSE_TAG

       CALL PRINT_EIGENVECTORS(NIONS,DOF,T_INFO_0%POSION,LATT_CUR%A, &
            EIGENVECTORS,EIGENVALUES,       &
            T_INFO%LSFOR,T_INFO%LSDYN,IU6)

       N=1
       DO I=1,T_INFO%NTYP
          DO J=1,T_INFO%NITYP(I)
             DO K=1,3
                CALL FIND_DOF_INDEX(NIONS,T_INFO%LSFOR,T_INFO%LSDYN,K,N,M)
                IF (M>0) EIGENVECTORS(M,:)=EIGENVECTORS(M,:)/SQRT(T_INFO%POMASS(I))
             END DO
             N=N+1
          END DO
       END DO

       IF (IU6>=0 .AND. IO%NWRITE>=3) THEN
          WRITE(IU6,*) "Eigenvectors after division by SQRT(mass)"
          CALL PRINT_EIGENVECTORS(NIONS,DOF,T_INFO_0%POSION,LATT_CUR%A, &
               EIGENVECTORS,EIGENVALUES,       &
               T_INFO%LSFOR,T_INFO%LSDYN,IU6)
       ENDIF
       DEALLOCATE(WORK)

       CALL STORE_RPA_SECOND_DERIV(SECOND_DERIV)
!=======================================================================
!
! ionic contribution to dielectric, elastic and piezoelectric tensor
!
!=======================================================================
       IF (DYN%ISIF>0 .OR. LEPSILON) THEN
          WRITE(IO%IU6,130)
          ! invert the matrix of the second derivatives
          SECOND_DERIV=-SECOND_DERIV
          CALL INV_SECOND_DERIV(SECOND_DERIV, DOF, IU6 )

          ! ionic contribution to macroscopic dielectric tensor
          IF (LBORN .AND. .NOT. LRPA .AND. IO%IU6>=0 ) THEN
             CALL EPSILON_ION( T_INFO, DOF, SECOND_DERIV, BORN_CHARGES_PEAD, EPSILON )

             ! induced polariation -> field
             EPSILON=EPSILON*EDEPS/(LATT_CUR%OMEGA)
             WRITE(IO%IU6,1100) 'IONIC CONTRIBUTION',EPSILON

             CALL XML_TENSOR("epsilon_ion",EPSILON)
             ! This is to keep consistency with replacement above
             CALL XML_TENSOR("dielectric_ion",EPSILON)
#ifdef VASP_HDF5
             CALL VH5_WRITE_DIELECTRIC_STATIC(IH5OUTFILEID,"ion",EPSILON)
#endif

             CALL EPSILON_ION_OMEGA( T_INFO, DOF, EIGENVECTORS, EIGENVALUES, BORN_CHARGES_PEAD, LATT_CUR%OMEGA )
          ENDIF

          IF (DYN%ISIF>0.AND.IO%IU6>=0) THEN
100       FORMAT(/ &
            A / &
            ' Direction', &
            4X,'XX', 10X,'YY', 10X,'ZZ', 10X,'XY', 10X,'YZ', 10X,'ZX'/ &
            ' --------------------------------------------------------------------------------'/ &
            ' XX     ',6F12.4/ &
            ' YY     ',6F12.4/ &
            ' ZZ     ',6F12.4/ &
            ' XY     ',6F12.4/ &
            ' YZ     ',6F12.4/ &
            ' ZX     ',6F12.4/ &
            ' --------------------------------------------------------------------------------'/)
             CALL ELASTIC_ION( T_INFO, DOF, SECOND_DERIV, ST, ELASTIC )

             ELASTICP(1,:,:)=ELASTIC(1,1,:,:)
             ELASTICP(2,:,:)=ELASTIC(2,2,:,:)
             ELASTICP(3,:,:)=ELASTIC(3,3,:,:)
             ELASTICP(4,:,:)=ELASTIC(1,2,:,:)
             ELASTICP(5,:,:)=ELASTIC(2,3,:,:)
             ELASTICP(6,:,:)=ELASTIC(3,1,:,:)

             FACT=EVTOJ*1E22_q/LATT_CUR%OMEGA

             WRITE(IU6,100) ' ELASTIC MODULI CONTR FROM IONIC RELAXATION (kBar)', ( &
                  (ELASTICP(J,I,I)*FACT,I=1,3), &
                   ELASTICP(J,1,2)*FACT,ELASTICP(J,2,3)*FACT,ELASTICP(J,3,1)*FACT,J=1,6)
          END IF

          ! ionic contribution to piezoelectric tensor
          IF ( DYN%ISIF>0 .AND. LBORN .AND. .NOT. LRPA .AND.IO%IU6>=0 ) THEN

             CALL PIEZO_ION( T_INFO, DOF, SECOND_DERIV, ST, BORN_CHARGES_PEAD, PIEZO )
             
             WRITE (IU6,180) 'PIEZOELECTRIC TENSOR IONIC CONTR'
             FACT=EVTOJ*1E20_q/LATT_CUR%OMEGA
             DO I =1,3
                WRITE (IU6,140) IDIR_TEXT(I),(PIEZO(I,J,J)*FACT,J=1,3), & 
                     PIEZO(I,1,2)*FACT,PIEZO(I,2,3)*FACT,PIEZO(I,3,1)*FACT
             ENDDO
#ifdef VASP_HDF5
             CALL VH5_WRITE_PIEZOELECTRIC_TENSOR(IH5OUTFILEID, "ion", PIEZO * FACT)
#endif
          ENDIF
       END IF
       DEALLOCATE(SECOND_DERIV, EIGENVECTORS, EIGENVALUES)
    END IF
    IF (IO%IU6>=0) WRITE(IO%IU6,130)

    DEALLOCATE(T_INFO_0%POSION)
    DEALLOCATE(INITIAL_FORCE)
    DEALLOCATE(DISPL_FORCES,BORN_CHARGES,BORN_CHARGES2, &
      PIEZO,INT_STRAIN)

    IF (IBRION==8) THEN
      DEALLOCATE(D, ND, IDIRD)
    ENDIF

    IF (IU0>=0 .AND. DOF>0) THEN
       WRITE (IU0,*) 'Linear response finished'
       WRITE (17,*) 'Linear response finished'
    ENDIF

    IF (IO%LOPEN) CALL WFORCE(IO%IU6)

    ! and now remove the symmetry operations that do not conserve the finite
    ! difference stencils of the PEAD routines (if so desired)
    IF (SYMM%ISYM>0.AND.LPEAD_SYM_RED()) &
         CALL QSTAR_SYM_REDUCE(REAL(RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/)),KIND=q),T_INFO,SYMM,IU6=-1)
  CONTAINS

!************************ SUBROUTINE STORE_DIPOLE_AND_AUTOCORRELATE ****
!
! only if the tolerance is very tight for the occupied states
! it is possible to attain a tight tolerance for the perturbed states
! with the present DIIS algorithm
!
!***********************************************************************

  SUBROUTINE RE_OPTIMIZE_ORBITALS
    REAL(q) EDIFF            ! break condition
    INTEGER NSIM
    INTEGER ICOUEV           ! number of optimizations
    REAL(q) ::  TOTEN_, RMS  ! temporary variables

#ifdef _OPENACC_TST
    PUSH_ACC_EXEC_ON(.TRUE.)

    ACC_ASYNC_Q=ACC_ASYNC_SYNC

    CALL ACC_COPYIN_TYPED_VAR(WDES)
    CALL ACC_COPYIN_TYPED_VAR(W)
    CALL ACC_COPYIN_TYPED_VAR(NONL_S)
    CALL ACC_COPYIN_TYPED_VAR(NONLR_S)

    CALL ACC_COPYIN_TYPED_VAR(GRIDC)
    CALL ACC_COPYIN_TYPED_VAR(GRID_SOFT)
    CALL ACC_COPYIN_TYPED_VAR(P)

!$ACC ENTER DATA COPYIN(SV,CDIJ,CQIJ) __IF_ASYNC__
!$ACC ENTER DATA COPYIN(LATT_CUR) __IF_ASYNC__
!$ACC ENTER DATA COPYIN(HAMILTONIAN,HAMILTONIAN%MU) IF(ASSOCIATED(HAMILTONIAN%MU).AND.ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
!$ACC ENTER DATA COPYIN(HAMILTONIAN,HAMILTONIAN%AVEC) IF(ASSOCIATED(HAMILTONIAN%AVEC).AND.ACC_EXEC_ON) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT IF(ACC_EXEC_ON)
#endif

    EDIFF = 1E-10
    NSIM=WDES%NSIM*2
#ifdef MPI
    NSIM=((WDES%NSIM*2+WDES%COMM_INTER%NCPU-1)/WDES%COMM_INTER%NCPU)*WDES%COMM_INTER%NCPU
#endif

    X=INFO%EBREAK

    ! for some non obvious reasons 
    IF (.NOT. LMAGBLOCH) THEN 
       INFO%EBREAK=0.25*EDIFF
    ELSE
       IF (IU0>=0) WRITE(IU0,*) 'WARNING: EBREAK remains at default, EDDAV has a problem for too tight convergence criteria'
       INFO%EBREAK=EDIFF
    ENDIF

    TOTEN=TOTENIN
    DO I=1,3
       CALL EDDAV(HAMILTONIAN,P, GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WDES, NSIM, &
            LMDIM,CDIJ,CQIJ, RMS,DESUM,ICOUEV, SV, E%EXHF, IO%IU6,IO%IU0, .FALSE., .TRUE., .FALSE.)

       E%EBANDSTR=BANDSTRUCTURE_ENERGY(WDES, W)
       TOTEN_=E%EBANDSTR+E%DENC+E%XCENC+E%TEWEN+E%PSCENC+E%EENTROPY+E%PAWPS+E%PAWAE+INFO%EALLAT+E%EXHF+Ediel_sol

       IF (IO%IU0>=0) THEN
          WRITE(IO%IU0,1000) I, TOTEN_, TOTEN_-TOTEN, DESUM, ICOUEV, RMS
          WRITE(17,1000)     I, TOTEN_, TOTEN_-TOTEN, DESUM, ICOUEV, RMS
       ENDIF
1000   FORMAT('DAV: ',I3,'   ',E20.12,'   ',E12.5,'   ',E12.5,I6,'  ',E10.3)
       IF(ABS(DESUM)<EDIFF.AND.ABS(TOTEN_-TOTEN)<EDIFF) EXIT

       TOTEN=TOTEN_
    ENDDO
    INFO%EBREAK=X    ! restore the break condition

#ifdef _OPENACC_TST
!$ACC WAIT IF(ACC_EXEC_ON)
!$ACC EXIT DATA DELETE(HAMILTONIAN%MU,HAMILTONIAN) IF(ASSOCIATED(HAMILTONIAN%MU).AND.ACC_EXEC_ON)
!$ACC EXIT DATA DELETE(HAMILTONIAN%AVEC,HAMILTONIAN) IF(ASSOCIATED(HAMILTONIAN%AVEC).AND.ACC_EXEC_ON)
!$ACC EXIT DATA DELETE(SV,CDIJ,CQIJ) IF(ACC_EXEC_ON)
!$ACC EXIT DATA DELETE(LATT_CUR) IF(ACC_EXEC_ON)

    CALL ACC_DELETE_TYPED_VAR(GRIDC)
    CALL ACC_DELETE_TYPED_VAR(GRID_SOFT)
    CALL ACC_DELETE_TYPED_VAR(P)
    CALL ACC_DELETE_TYPED_VAR(NONL_S)
    CALL ACC_DELETE_TYPED_VAR(NONLR_S)

!$ACC UPDATE SELF(W%CW,W%CPROJ) IF(ACC_EXEC_ON)
    CALL ACC_DELETE_TYPED_VAR(WDES)
    CALL ACC_DELETE_TYPED_VAR(W)

    POP_ACC_EXEC_ON
#endif

  END SUBROUTINE RE_OPTIMIZE_ORBITALS
  
  END SUBROUTINE LR_SKELETON

!************************ SUBROUTINE STORE_DIPOLE_AND_AUTOCORRELATE ****
!
! Routine to store the dipole moment as a function of time
!
!***********************************************************************

  SUBROUTINE STORE_DIPOLE_AND_AUTOCORRELATE( T_INFO, LATT_CUR, POTIM, NSW, NSTEP, TEMP, IO)
    USE constant
    USE lattice 
    USE poscar
    USE base
    USE mlrf_main
    IMPLICIT NONE
    TYPE (type_info)   T_INFO     ! information on positions etc.
    TYPE (latt)        LATT_CUR   ! lattice parameters
    REAL(q)            POTIM
    INTEGER            NSW        ! total number of steps
    INTEGER            NSW2       ! number of steps to calculate correlation function
    INTEGER            NSTEP      ! current step
    REAL(q)            TEMP       ! curent simulation temperature
    TYPE (in_struct)   IO
    REAL(q),PARAMETER :: PLANK=6.626075E-34
  ! local variables
    REAL(q) :: DIPOLE(3), TMP(3), A(3,3), OMEGA_MAX, OMEGA, DAMPING, SCALE
    INTEGER :: N, I, J, K
    COMPLEX(q), ALLOCATABLE :: EPSILON(:,:), SIGMA(:,:,:)
    REAL(q) :: SIGMA0(3,3)
    REAL(q), ALLOCATABLE :: JJ(:,:,:)
    COMPLEX(q) :: EXPFAC
    
    ! first allocate if required
    IF (.NOT. ALLOCATED(DIPOLE_TIME)) THEN
       ALLOCATE(DIPOLE_TIME(3, NSW))
       DIPOLE_TIME=0
    ENDIF

    ! first calculate current dipole
    DIPOLE=0
    
    IF (SQRT(DFIELD(1)**2+DFIELD(2)**2+DFIELD(3)**2)==0) THEN
       SCALE=1.0
    ELSE
       SCALE=1.0_q/SQRT(DFIELD(1)**2+DFIELD(2)**2+DFIELD(3)**2)
    ENDIF
    DO N=1,T_INFO%NIONS
       TMP=T_INFO%POSION(:,N)-POSION_LAST(:,N)
       TMP(1)=MOD(TMP(1)+4.5_q,1.0_q)-0.5_q
       TMP(2)=MOD(TMP(2)+4.5_q,1.0_q)-0.5_q
       TMP(3)=MOD(TMP(3)+4.5_q,1.0_q)-0.5_q
       ! convert to cartesian coordinates
       CALL DIRKAR(1, TMP, LATT_CUR%A)
       ! divide by the absolute magnitude of the external field
       ! to obtain response
       DIPOLE(:)=DIPOLE(:)+MATMUL(BORN_CHARGES_SAVE(:,:,N), TMP)*SCALE
    ENDDO
    ! store current positions into last position array
    POSION_LAST=T_INFO%POSION

    IF (NSTEP==1) THEN
       ! initial dipole is assumed to be zero
       DIPOLE_TIME(:,NSTEP)=0
    ELSE
       ! accumulate the change of the dipole over time
       DIPOLE_TIME(:,NSTEP)=DIPOLE(:)+DIPOLE_TIME(:,NSTEP-1)
    ENDIF

    !
    ! final processing for the delta-kick method
    !
    IF (NSTEP==NSW .AND. SQRT(DFIELD(1)**2+DFIELD(2)**2+DFIELD(3)**2)/=0) THEN
       ALLOCATE(EPSILON(3, NSW))
       EPSILON=0

       OMEGA_MAX=PI/POTIM
       DAMPING=4.0_q/(POTIM*NSW)

       !
       ! first write the time evolution of the dipole moment
       !
       IF (IO%IU6>=0) THEN
          DO I=0,NSW-1
             WRITE(77,'(4F14.7)') POTIM*I, DIPOLE_TIME(:,I+1)
          ENDDO
       ENDIF

       ! for the delta kick method, we obtain 3 components of EPSILON
       !
       ! FFT the dipole from time to frequency
       DO I=0,NSW-1
          OMEGA=OMEGA_MAX/NSW*I
          DO J=0,NSW-1
             EXPFAC=EXP( ((0._q,1._q)*OMEGA-DAMPING )*POTIM*J) !expfac=exponential-factor
             ! weight of each contribution is given by the time step POTIM
             EPSILON(:,I+1)=EPSILON(:,I+1)-1.0_q*DIPOLE_TIME(:,J+1)*EXPFAC*POTIM
          ENDDO
       ENDDO

       EPSILON=EPSILON*EDEPS/LATT_CUR%OMEGA

       IF (IO%IU6>=0) THEN
          DO I=0,NSW-1
             WRITE(78,'(7F14.7)') OMEGA_MAX/NSW*I*1000, -REAL(EPSILON(:,I+1),q)
          ENDDO
          WRITE(78,*)
          DO I=0,NSW-1
             WRITE(78,'(7F14.7)') OMEGA_MAX/NSW*I*1000, -AIMAG(EPSILON(:,I+1))
          ENDDO
       ENDIF

       DEALLOCATE(EPSILON)
    ENDIF
    !
    ! final processing for the Kubo-Green method
    ! calculate autocorrelation function  <dot P(0) dot P(t)>
    ! and Fourier transform to frequency
    !
    IF (NSTEP==NSW .AND. SQRT(DFIELD(1)**2+DFIELD(2)**2+DFIELD(3)**2)==0) THEN
       
       ALLOCATE(JJ(3,3,NSW),SIGMA(3, 3, NSW))

       OMEGA_MAX=PI/POTIM
       DAMPING=8.0_q/(POTIM*NSW)

       !
       ! calculate current = derivative of dipoles 
       ! overwriting the dipoles
       !
       DO I=1,NSW-2
          ! central differences
          DIPOLE_TIME(:,I)=(DIPOLE_TIME(:,I+2)-DIPOLE_TIME(:,I))/2/POTIM
       ENDDO
       !
       ! calculate the autocorrelation function over half of the timesteps
       !
       NSW2=(NSW-2)/2
       JJ=0
       DO I=1,NSW2  ! loop over initial time
          DO J=1,NSW2  ! loop over time differences
             ! better use of pipeline no inner do loop
             JJ(1, 1, J)=JJ(1, 1, J)+DIPOLE_TIME(1, I)*DIPOLE_TIME(1, I+J-1)
             JJ(1, 2, J)=JJ(1, 2, J)+DIPOLE_TIME(1, I)*DIPOLE_TIME(2, I+J-1)
             JJ(1, 3, J)=JJ(1, 3, J)+DIPOLE_TIME(1, I)*DIPOLE_TIME(3, I+J-1)
             JJ(2, 1, J)=JJ(2, 1, J)+DIPOLE_TIME(2, I)*DIPOLE_TIME(1, I+J-1)
             JJ(2, 2, J)=JJ(2, 2, J)+DIPOLE_TIME(2, I)*DIPOLE_TIME(2, I+J-1)
             JJ(2, 3, J)=JJ(2, 3, J)+DIPOLE_TIME(2, I)*DIPOLE_TIME(3, I+J-1)
             JJ(3, 1, J)=JJ(3, 1, J)+DIPOLE_TIME(3, I)*DIPOLE_TIME(1, I+J-1)
             JJ(3, 2, J)=JJ(3, 2, J)+DIPOLE_TIME(3, I)*DIPOLE_TIME(2, I+J-1)
             JJ(3, 3, J)=JJ(3, 3, J)+DIPOLE_TIME(3, I)*DIPOLE_TIME(3, I+J-1)
          ENDDO
       ENDDO
       JJ=JJ*(1.0_q/NSW2)
       !
       ! write the current-current autocorrelation function
       !
       IF (IO%IU6>=0) THEN
          DO I=1,NSW2
             WRITE(77,'(10F14.7)') POTIM*(I-1), JJ(1,1,I), JJ(2,2,I), JJ(3,3,I), JJ(1,2,I), JJ(2,3,I), JJ(3,1,I)
          ENDDO
       ENDIF


       SIGMA=0
       DO I=1,NSW
          OMEGA=OMEGA_MAX/NSW*(I-1)
          ! weight is the derivative of the function above
          DO J=1,NSW2
             EXPFAC=EXP( ((0._q,1._q)*OMEGA-DAMPING )*POTIM*(J-1)) !expfac=exponential-factor
             ! damping (Lorenzian broadening) + Hann window
!             EXPFAC=EXP( ((0._q,1._q)*OMEGA-DAMPING)*POTIM*(J-1))*0.5_q*(COS(PI*(J-1)/(NSW2-1))+1.0_q)
             ! weight is POTIM
             SIGMA(:,:,I)=SIGMA(:,:,I)+JJ(:,:,J)*EXPFAC*POTIM
          ENDDO
       ENDDO

       ! convert from e^2 A^2/A^3 /fs  to eV/fs
       SIGMA=SIGMA*EDEPS/LATT_CUR%OMEGA
       ! now multiply by beta = 1/(k_B T) (in eV) so units are now fs
       SIGMA=SIGMA/(TEMP*BOLKEV)

       ! epsilon= i/w sigma
       ! to obtain polarizability divide by frequency (i/w is in units of fs)
       ! remove constant part of the optical conductivity (imaginary part of the dielectric function)
       ! essentially it forces the conductivity to be zero, as it should be (but why is there any contribution)
       ! TODO: the conductivity sigma has a constant real valued offset, which should not be there
       SIGMA0=(SIGMA(:,:,1)+SIGMA(:,:,2))/2 ! one can alternatively also use SIGMA(:,:,NSW)
       SIGMA0=0
       IF (IO%IU6>=0) THEN
          DO I=2,NSW
             OMEGA=OMEGA_MAX/NSW*(I-1)
             SIGMA(:,:,I)=(SIGMA(:,:,I)-SIGMA0)/(OMEGA)*(0.0_q,1.0_q)
          ENDDO
          ! w=0 is missing copy from lowest frequency result
          SIGMA(:,:,1)=SIGMA(:,:,2)

          ! write the real part to unit 78
          DO I=1,NSW
             OMEGA=OMEGA_MAX/NSW*(I-1)
             WRITE(78,'(7F14.7)') OMEGA*1000, REAL(SIGMA(1,1,I),q),REAL(SIGMA(2,2,I),q),REAL(SIGMA(3,3,I),q), & 
                  REAL(SIGMA(1,2,I),q),REAL(SIGMA(2,3,I),q),REAL(SIGMA(3,1,I),q)
          ENDDO
          ! write the imaginary part to unit 79
          DO I=1,NSW
             OMEGA=OMEGA_MAX/NSW*(I-1)
             WRITE(79,'(7F14.7)') OMEGA*1000, AIMAG(SIGMA(1,1,I)), AIMAG(SIGMA(2,2,I)),AIMAG(SIGMA(3,3,I)), & 
                  AIMAG(SIGMA(1,2,I)),AIMAG(SIGMA(2,3,I)),AIMAG(SIGMA(3,1,I))
          ENDDO
       ENDIF

       DEALLOCATE(JJ,SIGMA)
    ENDIF

  END SUBROUTINE STORE_DIPOLE_AND_AUTOCORRELATE

  !>
  !> @brief
  !> Interpolate the eigenvalues and eigenvectors to a dense grid,
  !> compute the band velocities optical matrix elements as well
  !> as conductivity and linear response.
  !>
  SUBROUTINE INTERPOLATE_OPTICS_DENSE(HAMILTONIAN, E, KPOINTS, GRID, LATT_CUR, LATT_INI, & 
               T_INFO, NONLR_S, NONL_S, W, NEDOS_IN, LMDIM, P, SV, CQIJ, CDIJ, SYMM, INFO, IO, &
               GRIDC, GRIDUS, C_TO_US, IRDMAX)
#ifdef _OPENACC
   USE mopenacc_struct_def
#endif
   USE prec
   USE wave_high
   USE hamil_struct_def
   USE lattice
   USE poscar
   USE mpimy
   USE mgrid
   USE nonl_high
   USE base
   USE pseudo
   USE kpoints_change
   USE constant
   USE choleski
   USE subrot
   USE subrot_cluster
   USE mlr_optic
   USE mlrf_main
   USE mkpoints
   USE ini
   USE david
#ifdef VASP_HDF5
   USE vhdf5
#endif
   USE wave_interpolate
   USE msphpro, ONLY: SPHPRO, SPHPRO_FAST, WRITE_PROCAR
   USE tutor, ONLY: vtutor, isError, isAlert, NoKPAR, GRDSHFT_NOTET, SmallCSHIFT, NoOcc
   IMPLICIT NONE

   TYPE (ham_handle)     :: HAMILTONIAN
   TYPE (energy)         :: E
   TYPE (kpoints_struct) :: KPOINTS
   TYPE (grid_3d)        :: GRID
   TYPE (latt)           :: LATT_CUR
   TYPE (latt)           :: LATT_INI
   TYPE (type_info)      :: T_INFO
   TYPE (nonlr_struct)   :: NONLR_S
   TYPE (nonl_struct)    :: NONL_S
   TYPE (potcar)         :: P(:)
   TYPE (grid_3d)        :: GRIDC 
   TYPE (grid_3d)        :: GRIDUS
   TYPE (transit)        :: C_TO_US
   TYPE (info_struct)    :: INFO
   TYPE (in_struct)      :: IO
   TYPE (symmetry)       :: SYMM
   INTEGER               :: LMDIM, IRDMAX, NEDOS_IN
   TYPE (wavespin), TARGET :: W
   TYPE (eigenf_cluster_pointer),POINTER :: DEG_CLUSTER(:,:)
   OVERLAP :: CQIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
   OVERLAP :: CDIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
   RGRID   :: SV(DIMREAL(GRID%MPLWV),W%WDES%NCDIJ)

! local variables
   TYPE (kpoints_struct)    :: KPOINTS_INTER
   TYPE (wave_interpolator) :: WAVE_INTERP
   !> Wave desriptor for interpolated wavefunction
   TYPE (wavedes)  WDES_INTER
   !> Wavefunction container for interpolated wavefunction
   TYPE (wavespin) :: W_INTER
   !> Wavefunction container for computation of optical matrix elements
   TYPE (wavespin) :: WDIFF
   TYPE (skpoints_full),POINTER :: KPOINTS_INTER_FULL  
   TYPE (skpoints_full),POINTER :: KPOINTS_FULL_TEMP  
   INTEGER :: IKPT, I, IK, ISP, IDIR, IKDENSE, IKCOARSE
   INTEGER :: NEDOS
   INTEGER :: NSTEP
   LOGICAL :: LDONE
   REAL(q) :: EFERMI, ENTROPY
   REAL(q) :: DOS(NEDOS_IN,W%WDES%ISPIN),DOSI(NEDOS_IN,W%WDES%ISPIN)
   GDEFS, ALLOCATABLE   :: CHAM_INTER(:,:,:,:,:), CHAM_TMP(:,:,:,:)
   REAL(q), ALLOCATABLE :: EDER_INTER(:,:,:,:),EDER_INTER_TMP(:,:,:,:)
#ifdef VASP_HDF5
   integer(HID_T) :: fileid, groupid, subgroupid
   integer :: ierr
#endif
   ! for linear optics
   REAL(q) :: DELTAE, EMAX, EMAX_COND
   COMPLEX(q), ALLOCATABLE  :: EPSDD(:,:,:)
   COMPLEX(q), ALLOCATABLE  :: EPSJJ(:,:,:)
   REAL(q),ALLOCATABLE :: COND_ENERGY(:,:,:,:), EDOS(:)
   REAL(q)    :: CON(3,3), BIMAG(3,3)
   COMPLEX(q) :: WPLASMA_INTER(3,3)
   ! for projected density of states
   INTEGER :: JOBPAR
   INTEGER,PARAMETER :: LTRUNC=3
   INTEGER :: NT
   INTEGER :: LDIMP, LMDIMP
   INTEGER :: LPAR, NTYPD
   REAL(q),ALLOCATABLE :: PAR_INTER(:,:,:,:,:)
   REAL(q),ALLOCATABLE :: PAR(:,:,:,:,:)
   REAL(q),ALLOCATABLE :: DOSPAR(:,:,:,:)
   COMPLEX(q),ALLOCATABLE :: PHAS_INTER(:,:,:,:,:)
   COMPLEX(q),ALLOCATABLE :: PHAS(:,:,:,:,:)
   CHARACTER(LEN=3) :: LCHAR(4)=(/'  s','  p','  d','  f'/)
   CHARACTER(LEN=5) :: LMCHAR(16)=(/'    s','   py','   pz','   px','  dxy', &
      '  dyz','  dz2','  dxz','x2-y2','fy3x2',' fxyz',' fyz2','  fz3',' fxz2',' fzx2','  fx3'/)
   CHARACTER(LEN=:), ALLOCATABLE :: SUBGROUP

#ifdef _OPENACC
   LOGICAL :: ACC_EXEC_ALLOW_SAVE
!$ACC WAIT IF(ACC_EXEC_ON)
   ACC_EXEC_ALLOW_SAVE = ACC_EXEC_ALLOW
   ACC_EXEC_ALLOW = .FALSE.
#endif

   WAVE_INTERP%HAS_WAVEDER = .FALSE.
   CALL WAVE_INTERPOLATOR_WAVEDER(WAVE_INTERP,W,KPOINTS,GRID,LATT_INI,LATT_CUR,T_INFO,&
                 NONLR_S,NONL_S,LMDIM,GRIDC,GRIDUS,C_TO_US,IRDMAX,P,CQIJ,CDIJ,SV,SYMM,INFO,IO)

   CALL WAVE_INTERPOLATOR_INIT(WAVE_INTERP, KPOINTS, GRID, &
                 LATT_CUR, LATT_INI, T_INFO, NONL_S, W, P, SYMM, IO, KPOINTS_OPT_MODE)

   ! Generate dense grid of k-points
   CALL COPY_KPOINTS(KPOINTS,KPOINTS_INTER)
#ifdef oldsym
   CALL RD_KPOINTS(KPOINTS_INTER, LATT_CUR, &
        SYMM%ISYM>=0.AND..NOT. W%WDES%LSORBIT.AND..NOT. W%WDES%LSPIRAL, &
        SYMM%ISYM<0, -1, -1, FILE_NAME="KPOINTS_OPT")
#else
   CALL RD_KPOINTS(KPOINTS_INTER, LATT_CUR, &
        SYMM%ISYM>=0.AND..NOT. W%WDES%LNONCOLLINEAR, SYMM%ISYM<0, &
        -1, -1, FILE_NAME="KPOINTS_OPT")
#endif

   ! initialize WDES_INTER to fit KPOINTS_INTER
   WDES_INTER=W%WDES
   WDES_INTER%NKPTS= KPOINTS_INTER%NKPTS
   WDES_INTER%VKPT =>KPOINTS_INTER%VKPT
   WDES_INTER%WTKPT=>KPOINTS_INTER%WTKPT

   ! Allocate a container for interpolated eigenvalues
   CALL ALLOCW_NOPLANEWAVE(WDES_INTER, W_INTER)

   ! TODO: Check if this is needed! 
   IF (LVEL.OR.IO%LOPTICS) THEN
      ALLOCATE(EDER_INTER(W%WDES%NB_TOT,KPOINTS_INTER%NKPTS,W%WDES%ISPIN,3))
      EDER_INTER = 0
   END IF
   IF (IO%LOPTICS) THEN
      ALLOCATE(CHAM_TMP(W%WDES%NB_TOT,W%WDES%NB_TOT,W%WDES%NKPTS,W%WDES%ISPIN))
      ALLOCATE(CHAM_INTER(W%WDES%NB_TOT,W%WDES%NB_TOT,KPOINTS_INTER%NKPTS,W%WDES%ISPIN,3))
      CHAM_INTER = 0
   END IF
   ! To get orbital projections
   JOBPAR=1
   DO NT=1,T_INFO%NTYP
      IF (T_INFO%RWIGS(NT)<=0._q) JOBPAR=0
   ENDDO
   IF (JOBPAR/=0 .OR. IO%LORBIT>=10) THEN
      NTYPD = SIZE(P)
      CALL LDIM_PSEUDO(IO%LORBIT, NTYPD, P, LDIMP, LMDIMP)
      IF (IO%LORBIT==1 .OR. IO%LORBIT==2 .OR. (IO%LORBIT>=11 .AND. IO%LORBIT<=14)) THEN
         LPAR=LMDIMP
      ELSE
         LPAR=LDIMP
      ENDIF

      ALLOCATE(PAR(W%WDES%NB_TOT,W%WDES%NKPTS,LPAR,T_INFO%NIONP,W%WDES%NCDIJ))
      ALLOCATE(PHAS(LMDIMP,T_INFO%NIONS,W%WDES%NKPTS,W%WDES%NB_TOT,W%WDES%ISPIN))
      ALLOCATE(PAR_INTER(W%WDES%NB_TOT,WDES_INTER%NKPTS,LPAR,T_INFO%NIONP,W%WDES%NCDIJ))
      ALLOCATE(PHAS_INTER(LMDIMP,T_INFO%NIONS,WDES_INTER%NKPTS,W%WDES%NB_TOT,W%WDES%ISPIN))
      ALLOCATE(DOSPAR(NEDOS_IN,LPAR,T_INFO%NIONP,W%WDES%NCDIJ))
      PAR_INTER=0.0_q
      PHAS_INTER=0.0_q
      DOSPAR=0.0_q
    ELSE
      ALLOCATE(PHAS_INTER(1,1,1,1,W%WDES%ISPIN))
      ALLOCATE(PAR_INTER(1,1,1,1,W%WDES%NCDIJ))
      ALLOCATE(DOSPAR(1,1,1,W%WDES%NCDIJ))
   END IF

!************************************************************************
! loop over k-points
!************************************************************************
   DO IKPT=1,KPOINTS_INTER%NKPTS
      ! If the batch is not complete add this k-point to batch
      IF (WAVE_INTERP%NKPTS_BATCH < WAVE_INTERP%MAX_NKPTS_BATCH) THEN
         ! TODO: Also Check if point is part of the initial ones and skip it
         WAVE_INTERP%NKPTS_BATCH = WAVE_INTERP%NKPTS_BATCH + 1
         WAVE_INTERP%KPOINTS_BATCH(:,WAVE_INTERP%NKPTS_BATCH) = KPOINTS_INTER%VKPT(:,IKPT)
         WAVE_INTERP%KPOINTS_INDEX(WAVE_INTERP%NKPTS_BATCH) = IKPT
      ENDIF
      ! If the batch is complete or last k-point then launch computation
      IF (WAVE_INTERP%NKPTS_BATCH == WAVE_INTERP%MAX_NKPTS_BATCH .OR. IKPT==KPOINTS_INTER%NKPTS) THEN
         IF (IO%IU0>=0) THEN
            WRITE(IO%IU0,'(A,I0,A,I0,A,I0,A)') 'k-point batch [',WAVE_INTERP%KPOINTS_INDEX(1),'-',&
                    WAVE_INTERP%KPOINTS_INDEX(WAVE_INTERP%NKPTS_BATCH),'\',KPOINTS_INTER%NKPTS,']'
         END IF
         CALL WAVE_INTERPOLATOR_KPOINTS(WAVE_INTERP, HAMILTONIAN, GRID, &
                    LATT_CUR, LATT_INI, T_INFO, NONLR_S, NONL_S, W, LMDIM, P, SV, &
                    CQIJ, CDIJ, SYMM, INFO, IO%IU6, IO%IU0)

         ! Save eigenvalues
         DO ISP=1,W%WDES%ISPIN
            DO IK=1,WAVE_INTERP%NKPTS_BATCH
               IKCOARSE = WAVE_INTERP%NKPTS_ORIG+IK
               IKDENSE = WAVE_INTERP%KPOINTS_INDEX(IK)
               W_INTER%CELTOT(:,IKDENSE,ISP) = W%CELTOT(:,IKCOARSE,ISP)
            END DO
         END DO

         ! Use the wavefunctions to compute other stuff
         ! projections
         IF (JOBPAR/=0 .OR. IO%LORBIT>=10) THEN
            NSTEP=-1
            PAR=0.0_q
            IF (IO%LORBIT>=10) THEN
                CALL SPHPRO_FAST( GRID,LATT_CUR, P,T_INFO,W, W%WDES, 71,IO%IU6,&
                                  INFO%LOVERL,LMDIM,CQIJ, LPAR, LDIMP,LMDIMP,NSTEP, IO%LORBIT,PAR, &
                                  EFERMI, KPOINTS%EMIN, KPOINTS%EMAX, SYMM, LACUMULATE=.TRUE., PHAS_OUT=PHAS)
            ELSE
                CALL SPHPRO( &
                     GRID,LATT_CUR, P,T_INFO,W, W%WDES, 71,IO%IU6,&
                     INFO%LOVERL,LMDIM,CQIJ, LPAR, LDIMP, LMDIMP, LTRUNC, IO%LORBIT,PAR)
            ENDIF

            DO IK=1,WAVE_INTERP%NKPTS_BATCH
                IKCOARSE = WAVE_INTERP%NKPTS_ORIG+IK
#ifdef MPI
                IF (MOD(IKCOARSE-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
                IKDENSE = WAVE_INTERP%KPOINTS_INDEX(IK)
                PAR_INTER(:,IKDENSE,:,:,:) = PAR(:,IKCOARSE,:,:,:)
                IF (IO%LORBIT>=12 .AND. IO%LORBIT<=14) PHAS_INTER(:,:,IKDENSE,:,:) = PHAS(:,:,IKCOARSE,:,:)
            ENDDO
         ENDIF

         ! linear optics
         IF (LVEL.OR.IO%LOPTICS) THEN
            CALL ALLOCW(W%WDES,WDIFF)
             
            ! Find degenerate energies 
            NULLIFY(DEG_CLUSTER)
            CALL FIND_DEG_CLUSTERS(W%WDES, W, DEG_CLUSTER)
   
            ! Compute optical matrix elements for each direction 
            DO IDIR=1,3
               WDIFF%CW   =0
               WDIFF%CPROJ=0

               ! We use linear response only (to allow arbitrary grids)
               ! TODO: The derivatives for the original points are always recomputed
               CALL FOCK_K_DER_ANALYT(KPOINTS, GRID, LATT_CUR, LATT_INI, &
                    T_INFO,  NONLR_S, NONL_S, W, WDIFF, LMDIM, P, CQIJ, &
                    SYMM, IDIR, LDONE, IO%IU0, IO%IU6)

               CALL LRF_RPHI0( P,NONLR_S,NONL_S,W,LATT_CUR, &
                    T_INFO,INFO,IO,GRID,GRIDC,GRIDUS,C_TO_US,IRDMAX, &
                    CDIJ,CQIJ,SV,LMDIM,DEG_CLUSTER, IDIR, WDIFF, LDONE)
    
               ! Pull out the <u_n| r | u_m > for the calculation of the dielectric function
               IF (IO%LOPTICS) THEN
                  CALL INPROD_W(WDIFF,W,CHAM_TMP(:,:,:,:),INFO%LOVERL,IO%IU0)
               ENDIF

               DO ISP=1,W%WDES%ISPIN
                  DO IK=1,WAVE_INTERP%NKPTS_BATCH
                     IKCOARSE = WAVE_INTERP%NKPTS_ORIG+IK
#ifdef MPI
                     IF (MOD(IKCOARSE-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
                     IKDENSE = WAVE_INTERP%KPOINTS_INDEX(IK)
                     ! Store energy derivatives in the dense grid
                     EDER_INTER(:,IKDENSE,ISP,IDIR) = WDIFF%CELTOT(:,IKCOARSE,ISP)
                     ! Store the <u_n| r | u_m > in the dense grid
                     IF (IO%LOPTICS) CHAM_INTER(:,:,IKDENSE,ISP,IDIR) = CHAM_TMP(:,:,IKCOARSE,ISP)
                  END DO
               END DO
            ENDDO

            CALL DEALLOCW(WDIFF)
            CALL FREE_DEG_CLUSTERS(W%WDES, DEG_CLUSTER)
         ENDIF
         
         ! Restore to compute a new batch
         WAVE_INTERP%NKPTS_BATCH = 0

      ENDIF ! end if calculate velocities
   END DO ! end loop over k-points

   IF (IO%LOPTICS) DEALLOCATE(CHAM_TMP)
   CALL WAVE_INTERPOLATOR_FREE(WAVE_INTERP, KPOINTS, GRID, &
                 LATT_CUR, LATT_INI, T_INFO, NONL_S, W, P, SYMM)

!************************************************************************
! Final computations using interpolated wavefunction and matrix elements
! and writing of the results
!************************************************************************

   ! Related to eigenvalues only
   !-----------------------------------------------

   IF (JOBPAR/=0 .OR. IO%LORBIT>=10) THEN
      CALLMPI( M_sum_d( WDES_INTER%COMM, PAR_INTER, SIZE(PAR_INTER)))
      CALLMPI( M_sum_z( WDES_INTER%COMM, PHAS_INTER, SIZE(PHAS_INTER)))
      CALL DENSTA( IO%IU0, IO%IU6, WDES_INTER, W_INTER, KPOINTS_INTER, &
           INFO%NELECT, INFO%NUP_DOWN, ENTROPY, EFERMI, KPOINTS_INTER%SIGMA, &
           .TRUE., NEDOS_IN, LDIMP, T_INFO%NIONP, DOS, DOSI, PAR_INTER, DOSPAR)
   ELSE
      CALL DENSTA( IO%IU0, IO%IU6, WDES_INTER, W_INTER, KPOINTS_INTER, &
           INFO%NELECT, INFO%NUP_DOWN, ENTROPY, EFERMI, KPOINTS_INTER%SIGMA, &
           .TRUE., NEDOS_IN, 0, 0, DOS, DOSI, PAR_INTER, DOSPAR)
   ENDIF

   CALL XML_DOS(EFERMI, KPOINTS_INTER%EMIN, KPOINTS_INTER%EMAX, .FALSE., &
        DOS, DOSI, DOSPAR, NEDOS_IN, 1, 1, WDES_INTER%NCDIJ, &
        comment='kpoints_opt')

   CALL XML_TAG("eigenvalues_kpoints_opt",comment='kpoints_opt')
   CALL XML_TAG("kpoints")
   CALL XML_KPOINTS_LIST(KPOINTS_INTER%VKPT, KPOINTS_INTER%WTKPT)
   CALL XML_CLOSE_TAG("kpoints")
   CALL XML_EIGENVALUES_EXT(REAL(W_INTER%CELTOT),1, W%WDES%NB_TOT, KPOINTS_INTER%NKPTS, W%WDES%ISPIN)
   CALL XML_CLOSE_TAG("eigenvalues_kpoints_opt")

#ifdef VASP_HDF5
   CALL VH5_WRITE_DOS(IH5OUTFILEID, WDES_INTER, KPOINTS_INTER, DOS, DOSI, DOSPAR, EFERMI, &
                      T_INFO%NIONP, -1,SUBGROUP="electron_dos_kpoints_opt")
   CALL VH5_WRITE_EIGENVAL(IH5OUTFILEID, WDES_INTER, W_INTER, KPOINTS_INTER, &
                           SUBGROUP="electron_eigenvalues_kpoints_opt")

   ! Write projections
   IF (JOBPAR/=0 .OR. IO%LORBIT>=10) THEN
     IF (IO%LORBIT==11.OR.IO%LORBIT==12) THEN
       CALL VH5_WRITE_PROJECTORS(IH5OUTFILEID, T_INFO, WDES_INTER, W_INTER, IO%LORBIT, LPAR, PAR_INTER, LMCHAR, PHAS_INTER, &
                                 SUBGROUP="projectors_kpoints_opt")
     ELSE
       CALL VH5_WRITE_PROJECTORS(IH5OUTFILEID, T_INFO, WDES_INTER, W_INTER, IO%LORBIT, LPAR, PAR_INTER, LCHAR, PHAS_INTER, &
                                 SUBGROUP="projectors_kpoints_opt")
     ENDIF
     IF (IO%IU6>=0) THEN
       CALL XML_PROCAR(PAR, W_INTER%CELTOT, W_INTER%FERTOT, WDES_INTER%NB_TOT, WDES_INTER%NKPTS, LPAR,&
                     T_INFO%NIONP, WDES_INTER%NCDIJ, TAG="projected_kpoints_opt")
       OPEN(UNIT=99,FILE='PROCAR_OPT',STATUS='UNKNOWN')
       CALL WRITE_PROCAR(99,W_INTER,WDES_INTER,IO%LORBIT,LPAR,LMDIMP,PAR_INTER,PHAS_INTER,T_INFO,KPOINTS%EMIN,KPOINTS%EMAX)
       CLOSE(99)
     ENDIF
   ENDIF
#endif

   ! Compute IP optical spectra and conductivity on dense grid
   !-----------------------------------------------------------------
   EMAX=MAX_ENERGY_UNOCCUPIED(W%WDES,W)*1.2; IF (OMEGAMAX_OPTIC/=-1) EMAX=OMEGAMAX_OPTIC
   EMAX_COND=5; IF (OMEGAMAX_OPTIC/=-1) EMAX_COND=OMEGAMAX_OPTIC

   ! smaller values are dangerous and fail more often than not
   NEDOS=MAX(NEDOS_IN,1000)

   DELTAE=EMAX/(NEDOS-1)
   ! This one is very nasty
   NBANDS_CDER = WDES_INTER%NB_TOT

   ALLOCATE(COND_ENERGY(NEDOS,3,3,W%WDES%NCDIJ))
   ALLOCATE(EDOS(NEDOS))

   DO I=1,NEDOS
      EDOS(I)=EFERMI-EMAX_COND+(EMAX_COND*2)/(NEDOS-1)*(I-1)
   ENDDO

   IF (IO%LOPTICS) THEN
#ifdef gammareal
      CALLMPI( M_sum_single(WDES_INTER%COMM_KINTER, CHAM_INTER, SIZE(CHAM_INTER)) )
#else
      CALLMPI( M_sum_single(WDES_INTER%COMM_KINTER, CHAM_INTER, 2*SIZE(CHAM_INTER)) )
#endif
      CALLMPI( M_sum_d(WDES_INTER%COMM_KINTER, EDER_INTER, SIZE(EDER_INTER)) )
      ALLOCATE(EPSDD(NEDOS,3,3),EPSJJ(NEDOS,3,3))
      CALL LRF_EPSILON_COMPUTE(WDES_INTER, W_INTER, LATT_CUR, LATT_INI, INFO, IO, KPOINTS_INTER, SYMM, &
                               EFERMI, DELTAE, EMAX, EMAX_COND, NEDOS, &
                               EDER_INTER, CHAM_INTER, &
                               WPLASMON, CON, BIMAG, WPLASMA_INTER, EPSDD, EPSJJ, COND_ENERGY )

      CALL LR_OPTIC_OUTPUT(EPSDD,EPSJJ,WPLASMA_INTER,CON,BIMAG,NEDOS,DELTAE,KPOINTS,LATT_CUR,INFO,IO)

      DO I=1,WDES_INTER%ISPIN
         IF (I==1) THEN
            CALL XML_EPSILON_COND(EDOS, COND_ENERGY(1,1,1,I), NEDOS, "spin=1" )
         ELSE
            CALL XML_EPSILON_COND(EDOS, COND_ENERGY(1,1,1,I), NEDOS, "spin=2" )
         ENDIF
      ENDDO

#ifdef VASP_HDF5
      SUBGROUP = "linear_response_kpoints_opt"
      CALL VH5_WRITE_DIELECTRIC_ENERGIES(IH5OUTFILEID, DELTAE, NEDOS, SUBGROUP)
      CALL VH5_WRITE_DIELECTRIC_DYNAMIC(IH5OUTFILEID, EPSDD, "density_density", SUBGROUP)
      CALL VH5_WRITE_DIELECTRIC_DYNAMIC(IH5OUTFILEID, EPSJJ, "current_current", SUBGROUP)
      CALL VH5_WRITE_CONDUCTIVITY_DYNAMIC(IH5OUTFILEID, NEDOS, EDOS, COND_ENERGY, EFERMI)
#endif
      DEALLOCATE(EPSDD,EPSJJ)
   ELSEIF (LVEL.AND.(.NOT.IO%LOPTICS)) THEN
      CALL CONDUCTIVITY_ENERGY_RESOLVED( WDES_INTER, W_INTER, EDER_INTER, EFERMI, &
                                        NEDOS, EMAX_COND, COND_ENERGY, INFO%NELECT, INFO%NUP_DOWN, &
                                        KPOINTS_INTER, SYMM, LATT_CUR, RTIME)
#ifdef VASP_HDF5
      CALL VH5_WRITE_CONDUCTIVITY_DYNAMIC(IH5OUTFILEID, NEDOS, EDOS, COND_ENERGY, EFERMI)
#endif
   ENDIF
   DEALLOCATE(COND_ENERGY,EDOS)

   ! Write optical matrix elements
   IF (LVEL) THEN
      CALL XML_TAG("velocities", comment="kpoints_opt")
      CALL XML_TAG("kpoints")
      CALL XML_KPOINTS_LIST(KPOINTS_INTER%VKPT, KPOINTS_INTER%WTKPT)
      CALL XML_CLOSE_TAG("kpoints")
      ! TODO: change this ugliness
      ALLOCATE(EDER_INTER_TMP(3,W%WDES%NB_TOT,KPOINTS_INTER%NKPTS,W%WDES%ISPIN))
      EDER_INTER_TMP(1,:,:,:) = EDER_INTER(:,:,:,1)
      EDER_INTER_TMP(2,:,:,:) = EDER_INTER(:,:,:,2)
      EDER_INTER_TMP(3,:,:,:) = EDER_INTER(:,:,:,3)
      CALL XML_EIGENVALUES_EXT(EDER_INTER_TMP, 3, W%WDES%NB_TOT, KPOINTS_INTER%NKPTS, W%WDES%ISPIN)
      DEALLOCATE(EDER_INTER_TMP)
      CALL XML_CLOSE_TAG("velocities")
#ifdef VASP_HDF5
      CALL VH5_WRITE_BAND_VELOCITIES(IH5OUTFILEID, EDER_INTER, SUBGROUP="electron_eigenvalues_kpoints_opt")
#endif
      DEALLOCATE(EDER_INTER)
   ENDIF

#ifdef VASP_HDF5
   IF (IO%LOPTICS) THEN
      DEALLOCATE(CHAM_INTER)
   END IF
#endif

#ifdef _OPENACC
   ACC_EXEC_ALLOW = ACC_EXEC_ALLOW_SAVE
#endif

  END SUBROUTINE INTERPOLATE_OPTICS_DENSE
    
END MODULE mlr_main
