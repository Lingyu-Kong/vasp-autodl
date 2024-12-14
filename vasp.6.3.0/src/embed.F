#include "symbol.inc"

!***********************************************************
!
! This module shall provide embedding routines
!  - construct localized orbitals
!  - set up cluster, entangled, and rest orbitals
!  - writing out WAVECARS and transformation matrices
! Tobias Schaefer
!
! The following embedding techniques are available 
!   EMBED_TYPE = 1: localization only in the occupied
!                   manifold. Excitations possible from
!                   localized orbitals to entire virtual 
!                   space.
!   EMBED_TYPE = 1: creates localized and bath orbitals.
!                   localized orbitals include also virtual
!                   manifold and bath is restricted to 
!                   occupied space projected onto the space
!                   of localized orbitals.
!
! ToDo: parallelization not yet supported
!
!************************************************************

MODULE embed
    USE base
    USE wave
    USE locproj
  IMPLICIT NONE

  PUBLIC :: IS_EMBED ! check if module is active
  
  INTEGER              :: EMBED_TYPE      ! specifies the embedding technique 
  INTEGER, ALLOCATABLE :: EMBED_SITES(:)
  INTEGER, ALLOCATABLE :: EMBED_STATES(:)
  INTEGER              :: NEMBED_STATES   ! = SIZE(EMBED_STATES) 
  LOGICAL              :: LDUMPWTRAFO
  LOGICAL              :: LEMBED          ! to activate this module

  CONTAINS

  !********************* SUBROUTINE EMBED_READER **********************
  ! reads the following flags
  !  EMBED_SITES (or preferably EMBED_STATES)
  !  EMBED_TYPE
  !********************************************************************
  SUBROUTINE EMBED_READER(IO)
    USE reader_tags
    USE tutor, ONLY: vtutor
    
    TYPE(in_struct) :: IO
    ! local
    INTEGER :: IERR, IDUM(1)
    INTEGER :: N_LOCAL_FUNC
    INTEGER :: I, J
    LOGICAL :: LOPEN
    INTEGER :: NSITES

    CALL OPEN_INCAR_IF_FOUND(IO%IU5, LOPEN)

    LEMBED=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LEMBED', LEMBED, IERR, WRITEXMLINCAR)
    IF (.NOT. LEMBED) RETURN
     
    ! read LDUMPWTRAFO
    LDUMPWTRAFO=.FALSE.
    CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LDUMPWTRAFO', LDUMPWTRAFO, IERR, WRITEXMLINCAR)

    ! read EMBED_TYPE
    CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'EMBED_TYPE', EMBED_TYPE, &
                       IERR, WRITEXMLINCAR)
    IF(EMBED_TYPE .LT. 1 .OR. EMBED_TYPE .GT. 6) THEN
      WRITE(IO%IU0,*) "ERROR: illegal value for EMBED_TYPE in INCAR", EMBED_TYPE
      STOP
    ENDIF

    ! determine the number of arguments of EMBED_STATES 
    CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'EMBED_STATES', IDUM, 1, IERR, &
                       LWRITEXML=.FALSE., FOUNDNUMBER = NEMBED_STATES)

    ! if EMBED_STATES is set, use it, otherwise we try to find EMBED_SITES
    IF(NEMBED_STATES .GE. 1) THEN
      ALLOCATE(EMBED_STATES(NEMBED_STATES))
      ! read arguments of EMBED_STATES
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'EMBED_STATES', EMBED_STATES, NEMBED_STATES, &
                         IERR, WRITEXMLINCAR)
      ! write to EMBED_STATES to OUTCAR
      IF(IO%IU6 .GE. 0) THEN
        WRITE(IO%IU6,'(/A$)') " EMBED_STATES = "
        DO I = 1, NEMBED_STATES
          WRITE(IO%IU6,'(I5$)') EMBED_STATES(I)
        ENDDO
        WRITE(IO%IU6,'(/)')
      ENDIF

    ELSE 

      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'EMBED_SITES', IDUM, 1, IERR, &
                         LWRITEXML=.FALSE., FOUNDNUMBER = NSITES)
      IF(NSITES .LT. 1) THEN
        CALL vtutor%error("ERROR: Neither NTARGET_STATES nor NTARGET_SITES is set in the INCAR file")
      ENDIF

      ALLOCATE(EMBED_SITES(NSITES))
      ! read arguments of EMBED_SITES
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'EMBED_SITES', EMBED_SITES, NSITES, &
                         IERR, WRITEXMLINCAR)

      ! get total number of local functions
      N_LOCAL_FUNC = LPRJ_GET_NUM_WAN()

      ! get number of selected states 
      NEMBED_STATES = 0
      DO I = 1, N_LOCAL_FUNC
        IF( ANY(EMBED_SITES .EQ. LPRJ_GET_SITE(I)) ) THEN
          NEMBED_STATES = NEMBED_STATES + 1
        ENDIF
      ENDDO

      ! collect selected states in EMBED_SITES
      ALLOCATE(EMBED_STATES(NEMBED_STATES))
      J = 0
      DO I = 1, N_LOCAL_FUNC
        IF( ANY(EMBED_SITES .EQ. LPRJ_GET_SITE(I)) ) THEN
          J = J + 1
          EMBED_STATES(J) = I
        ENDIF
      ENDDO

      ! write to EMBED_SITES and EMBED_STATES to OUTCAR
      ! but not if EMBED_TYPE=6
      IF((IO%IU6 .GE. 0) .AND. (EMBED_TYPE .NE. 6)) THEN
        WRITE(IO%IU6,'(/A,I5)') " #EMBED_SITES = ", SIZE(EMBED_SITES)
        WRITE(IO%IU6,'(A$)') " EMBED_SITES = "
        DO I = 1, NSITES
          WRITE(IO%IU6,'(I5$)') EMBED_SITES(I)
        ENDDO
        WRITE(IO%IU6,*)
        WRITE(IO%IU6,'(A,I5)') " #EMBED_STATES = ", SIZE(EMBED_STATES)
        WRITE(IO%IU6,'(A$)') " EMBED_STATES = "
        DO I = 1, NEMBED_STATES
          WRITE(IO%IU6,'(I5$)') EMBED_STATES(I)
        ENDDO
        WRITE(IO%IU6,'(/)')
      ENDIF

      CALL CLOSE_INCAR_IF_FOUND(IO%IU5)

    ENDIF

  END SUBROUTINE EMBED_READER



  FUNCTION IS_EMBED() RESULT(IS_LEMBED)
    LOGICAL :: IS_LEMBED
    IS_LEMBED = LEMBED
  END FUNCTION IS_EMBED



  !************** SUBROUTINE EMBED_SETUP_ORBITALS ******************
  ! decides which SETUP_ORBITAL routine should be called
  ! according to EMBED_TYPE
  !*****************************************************************
  SUBROUTINE EMBED_SETUP_ORBITALS(W, GRID, P, CQIJ, LATT_CUR, HAMILTONIAN, &
               NONLR_S, NONL_S, SYMM, CDIJ, SV, T_INFO, INFO, IO)
    USE lattice
    USE pseudo
    USE poscar
    USE fileio
    USE nonl_high
    USE hamil_struct_def
    USE subrot
    USE choleski
#ifdef PROFILING
    USE profiling
#endif

    TYPE(wavespin)      :: W
    TYPE(grid_3d)       :: GRID
    TYPE(potcar)        :: P(:)
    OVERLAP             :: CQIJ(:, :, :, :)
    TYPE(latt)          :: LATT_CUR
    TYPE (ham_handle)   :: HAMILTONIAN
    TYPE (nonlr_struct) :: NONLR_S
    TYPE (nonl_struct)  :: NONL_S
    TYPE (symmetry)     :: SYMM
    OVERLAP             :: CDIJ(:, :, :, :)
    RGRID               :: SV(DIMREAL(GRID%MPLWV),W%WDES%NCDIJ) 
    TYPE(type_info)     :: T_INFO
    TYPE(info_struct)   :: INFO
    TYPE(in_struct)     :: IO

    ! status output
    IF (IO%IU0 >= 0) THEN
       WRITE(IO%IU0,'(A)') " start embed ..." 
    ENDIF

    IF(EMBED_TYPE .EQ. 1) THEN
      CALL EMBED_SETUP_ORBITALS_M1(W, GRID, P, CQIJ, LATT_CUR, &
                                   T_INFO, INFO, IO)
    ELSEIF(EMBED_TYPE .EQ. 2)  THEN
      CALL EMBED_SETUP_ORBITALS_M2(W, GRID, P, CQIJ, LATT_CUR, & 
                                   HAMILTONIAN, NONLR_S, NONL_S, &
                                   SYMM, CDIJ, SV, T_INFO, INFO, IO)
    ELSEIF(EMBED_TYPE .EQ. 4)  THEN
      CALL EMBED_SETUP_ORBITALS_M4(W, GRID, P, CQIJ, LATT_CUR, &
                                   T_INFO, INFO, IO)
    ELSEIF(EMBED_TYPE .EQ. 5)  THEN
      CALL EMBED_DUMP_HAMIL(W, GRID, HAMILTONIAN, NONLR_S, NONL_S, &
                                   SYMM, P, CDIJ, CQIJ, SV, LATT_CUR, T_INFO, INFO, IO)
    ELSEIF(EMBED_TYPE .EQ. 6)  THEN
      CALL EMBED_REORDER_ORBITALS(W, T_INFO, LATT_CUR, IO)
    !ELSEIF(EMBED_TYPE .EQ. 3) THEN
    !  CALL EMBED_SETUP_ORBITALS_SCDM(W, GRID, P, CQIJ, LATT_CUR, &
    !                                 T_INFO, INFO, IO)
    ENDIF

    ! status output
    IF (IO%IU0 >= 0) THEN
       WRITE(IO%IU0,'(A)') " ... finished embed" 
    ENDIF

!#ifdef PROFILING
!     ! started in INIT_PROFILING
!     PROFILING_STOP('total_time')
!     IF (IO%IU0>=0) THEN
!        CALL PRINT_FULL_PROFILE(IO%IU6, PFLAT=.FALSE.)
!     ENDIF
!#endif
!     CALLMPI_C(M_exit())
!    STOP
  ENDSUBROUTINE EMBED_SETUP_ORBITALS




  !************** SUBROUTINE EMBED_SETUP_ORBITALS_M1 ******************
  ! routine for method I
  ! This routine calculates local (cluster) orbitals based
  ! on selected sites via EMBED_STATES using a projection approach.
  !********************************************************************
  SUBROUTINE EMBED_SETUP_ORBITALS_M1(W, GRID, P, CQIJ, LATT_CUR, &
                                     T_INFO, INFO, IO)
    USE lattice
    USE pseudo
    USE poscar
    USE fileio
    USE nonl_high
    USE hamil_struct_def
    USE subrot
    USE choleski
#ifdef PROFILING
    USE profiling
#endif

    TYPE(wavespin)      :: W
    TYPE(grid_3d)       :: GRID
    TYPE(potcar)        :: P(:)
    OVERLAP             :: CQIJ(:, :, :, :)
    TYPE(latt)          :: LATT_CUR
    TYPE(type_info)     :: T_INFO
    TYPE(info_struct)   :: INFO
    TYPE(in_struct)     :: IO
    ! local
    GDEF, ALLOCATABLE :: TRAFO_MATRIX(:, :, :, :)
    GDEF, ALLOCATABLE :: WORK_MATRIX(:, :, :, :)
    TYPE(wavefuna)    :: WA
    TYPE(wavedes1)    :: WDESK
    INTEGER           :: NB_TOT, NOCC, NWAN, NREST
    INTEGER           :: I, J, L, IK, ISP
    INTEGER           :: REST_IDX, LOC_IDX

#ifdef PROFILING
    CALL START_PROFILING('EMBED_SETUP_ORBITALS_M1')
#endif

    ! GAMMA ONLY and SPIN UNPOLARIZED
    IK = 1
    ISP = 1
    IF(W%WDES%ISPIN > 1 .OR. W%WDES%NKPTS > 1) THEN
      WRITE(*,*) "ERROR: At the moment, embedding is only supported for &
                 &ISPIN=1 at the GAMMA-point."
      STOP
    ENDIF
        
    NB_TOT = W%WDES%NB_TOT
    NOCC = NINT(SUM(W%FERTOT(:, IK, ISP)))
    NWAN = LPRJ_GET_NUM_WAN()
    NREST = MAX(0, NOCC - NEMBED_STATES)

    IF(NOCC .NE. NWAN) THEN
      WRITE(*,*) "ERROR: EMBED_TYPE=1 requires NOCC = NWAN, but ", NOCC, NWAN
      STOP
    ENDIF

    !! W_EMBED is allocted to store up to NB_TOT bands
    !WDES_EMBED=W%WDES
    !CALL ALLOCW(WDES_EMBED, W_EMBED)
    !W_EMBED%CELTOT=W%CELTOT
    !W_EMBED%FERTOT=W%FERTOT
      

    !---------------
    ! calc trafo to the local orbitals (cluster orbitals) by 
    ! projecting local functions onto occupied states and store 
    ! transofmration matrix in WORK_MATRIX
    !---------------
    CALL LPRJ_PROALL(W, W%WDES, P, CQIJ, LATT_CUR, LPRJ_functions, LPRJ_COVL, T_INFO, INFO, IO)
    !CALL LPRJ_CALC_WANNIER_TRAFO(W, WORK_MATRIX, IO)
    CALL LPRJ_SVD_LOCALIZE(WORK_MATRIX, W%WDES%ISPIN, [1,0], [NOCC,0], IO)

    ! reorder wannier states such that the target states are 
    ! below the fermi energy -> store in TRAFO_MATRIX
    ALLOCATE( TRAFO_MATRIX(NB_TOT, NOCC, IK, ISP) )
    TRAFO_MATRIX = 0 ! nullify
    REST_IDX = 1
    LOC_IDX = NREST + 1
    DO I = 1, NOCC
      IF( ANY(EMBED_STATES .EQ. I) ) THEN
        TRAFO_MATRIX(:,LOC_IDX,IK,ISP) = WORK_MATRIX(:,I,IK,ISP)
        LOC_IDX = LOC_IDX + 1
      ELSE
        TRAFO_MATRIX(:,REST_IDX,IK,ISP) = WORK_MATRIX(:,I,IK,ISP)
        REST_IDX = REST_IDX + 1
      ENDIF
    ENDDO

    !IF(LDUMPWTRAFO) THEN
    !  IF(IO%IU6 > 0) THEN
    !    WRITE(IO%IU6,*) 
    !    WRITE(IO%IU6,*) "Writing out TRAFO_MATRIX:"
    !    WRITE(IO%IU6,*) "# U  I  L  ISPI ISPJ"
    !    DO I = 1, NOCC
    !      DO L = 1, NOCC
    !        WRITE(IO%IU6,*) TRAFO_MATRIX(I,L,IK,ISP), I, L, IK, ISP
    !      ENDDO
    !    ENDDO
    !    WRITE(IO%IU6,*) 
    !  ENDIF
    !ENDIF
    IF(LDUMPWTRAFO) THEN
      IF(IO%IU6 > 0) THEN
        CALL WRITE_CC4S_GMATRIX(TRAFO_MATRIX(1:NOCC,:,IK,ISP), NOCC, NWAN, "trafo.dat", 1E-12_q)
      ENDIF
    ENDIF

    ! setup embedding orbitals using TRAFO_MATRIX

    CALL SETWDES(W%WDES, WDESK, IK)
    WA = ELEMENTS(W, WDESK, ISP)

    IF (W%WDES%DO_REDIS) THEN
       CALL REDISTRIBUTE_PROJ(WA)
       CALL REDISTRIBUTE_PW(WA)
    ENDIF

    CALL LINCOM('F', WA%CW_RED(:,:), WA%CPROJ_RED(:,:), TRAFO_MATRIX(:,:,IK,ISP), &
                NB_TOT, NOCC, WDESK%NPL_RED, WDESK%NPRO_RED, WDESK%NRPLWV_RED, &
                WDESK%NPROD_RED, NB_TOT, WA%CW_RED(:,:), WA%CPROJ_RED(:,:))
    
    IF (W%WDES%DO_REDIS) THEN
       CALL REDISTRIBUTE_PROJ(WA)
       CALL REDISTRIBUTE_PW(WA)
    ENDIF

    CALL OUTWAV(IO, W%WDES, W, LATT_CUR, 0.0_q, 'emb')

#ifdef PROFILING
    CALL STOP_PROFILING('EMBED_SETUP_ORBITALS_M1')
#endif

  
  ENDSUBROUTINE EMBED_SETUP_ORBITALS_M1




  !************** SUBROUTINE EMBED_SETUP_ORBITALS_M2 ******************
  ! routine for method II
  ! This routine calculates local (cluster) and entangled (bath) orbitals based
  ! on selected sites via EMBED_STATES using a projection approach.
  ! Furthermore the local+entangled orbitals are orthogonalized and 
  ! canonicalized, i.e. the Hamiltonian is re-diagonalized in this basis.
  ! The canonicalized orbitals are written into WAVECAR.emb
  ! Also a unitary matrix is written to XXX in order to transform from the 
  ! canonicalized to the local+entangled orbitals
  !
  ! ! GAMMA ONLY 
  ! ToDo: written WAVECARS are too large! OUTWAV seems to ignore NB_TOT, why?
  !****************************************************************************
  SUBROUTINE EMBED_SETUP_ORBITALS_M2(W, GRID, P, CQIJ, LATT_CUR, HAMILTONIAN, &
               NONLR_S, NONL_S, SYMM, CDIJ, SV, T_INFO, INFO, IO)
    USE lattice
    USE pseudo
    USE poscar
    USE fileio
    USE nonl_high
    USE hamil_struct_def
    USE subrot
    USE choleski
    USE tutor, ONLY: vtutor
#ifdef PROFILING
    USE profiling
#endif

    TYPE(wavespin)      :: W
    TYPE(grid_3d)       :: GRID
    TYPE(potcar)        :: P(:)
    OVERLAP             :: CQIJ(:, :, :, :)
    TYPE(latt)          :: LATT_CUR
    TYPE (ham_handle)   :: HAMILTONIAN
    TYPE (nonlr_struct) :: NONLR_S
    TYPE (nonl_struct)  :: NONL_S
    TYPE (symmetry)     :: SYMM
    OVERLAP             :: CDIJ(:, :, :, :)
    RGRID               :: SV(DIMREAL(GRID%MPLWV),W%WDES%NCDIJ) 
    TYPE(type_info)     :: T_INFO
    TYPE(info_struct)   :: INFO
    TYPE(in_struct)     :: IO
    ! local
    GDEF, ALLOCATABLE :: TRAFO_MATRIX(:, :, :, :)
    GDEF, ALLOCATABLE :: WORK_MATRIX(:, :, :, :)
    GDEF, ALLOCATABLE :: CLUSTER_PROJECTOR(:, :, :, :)
    GDEF, ALLOCATABLE :: OVERLAPS(:,:)
    TYPE(wavespin)    :: W_EMBED, W_TEMP
    TYPE(wavedes)     :: WDES_EMBED
    INTEGER           :: NB_TOT, NOCC, NREST
    INTEGER           :: I, J, L, IK, ISP
    TYPE (energy)     :: E

#ifdef PROFILING
    CALL START_PROFILING('EMBED_SETUP_ORBITALS_M2')
#endif

    ! GAMMA ONLY and SPIN UNPOLARIZED
    IK = 1
    ISP = 1
    IF(W%WDES%ISPIN > 1 .OR. W%WDES%NKPTS > 1) THEN
      CALL vtutor%error("ERROR: At the moment, embedding is only supported for ISPIN=1 at the &
         &GAMMA-point.")
    ENDIF
        
    NB_TOT = W%WDES%NB_TOT
    NOCC = NINT(SUM(W%FERTOT(:, IK, ISP)))
    NREST = MAX(0, NOCC - NEMBED_STATES)

    ! W_EMBED is allocted to store up to NB_TOT bands
    WDES_EMBED=W%WDES
    CALL ALLOCW(WDES_EMBED, W_EMBED)
    W_EMBED%CELTOT=W%CELTOT
    W_EMBED%FERTOT=W%FERTOT

    ALLOCATE( TRAFO_MATRIX(NB_TOT, NOCC+NEMBED_STATES, 1, 1) )
    TRAFO_MATRIX = 0 ! nullify


    !---------------
    ! first step: calc trafo for the local orbitals (cluster orbitals) by 
    ! projecting local functions onto occupied+few unoccupied states and
    ! store transofmration matrix in TRAFO_MATRIX
    !---------------
    CALL LPRJ_PROALL(W, W%WDES, P, CQIJ, LATT_CUR, LPRJ_functions, LPRJ_COVL, T_INFO, INFO, IO)
    CALL LPRJ_CALC_WANNIER_TRAFO(W, LPRJ_COVL, WORK_MATRIX, IO)
    CALL LPRJ_LINCOM(W, WORK_MATRIX, W_EMBED)
    CALL OUTWAV(IO, W_EMBED%WDES, W_EMBED, LATT_CUR, 0.0_q, 'wan')
    
    ! store transformation data for the selected localized orbitals 
    ! position the localized orbitals at the end of the TRAFO_MATRIX
    DO I = 1, NEMBED_STATES
      TRAFO_MATRIX(:, NOCC+I, IK, ISP) = WORK_MATRIX(:, EMBED_STATES(I), IK, ISP)
    ENDDO


    !---------------
    ! second step: calc trafo for entangled orbitals by projecting selected cluster
    ! orbitals onto all occupied states. But if NEMBED_STATES >= NOCC do nothing.
    !---------------
    
    ! the overlap between occupied and selected cluster orbitals is equal to 
    ! the entries of the transofmration matrix
    ALLOCATE( OVERLAPS(NOCC, NEMBED_STATES) )
    DO I = 1, NEMBED_STATES
      OVERLAPS(1:NOCC, I) = TRAFO_MATRIX(1:NOCC, NOCC+I, IK, ISP)
    ENDDO

    DEALLOCATE(WORK_MATRIX)
    ALLOCATE( WORK_MATRIX(NOCC,NOCC, 1, 1) )

    ! set up the transofmration for the entangled and rest orbitals 
    IF(NOCC > NEMBED_STATES) THEN
      CALL SETUP_ENTGLT_REST_TRAFO(OVERLAPS, WORK_MATRIX)
    ELSE ! nothing to do -> simply identity matrix
      WORK_MATRIX = 0 ! nullify
      DO I = 1, NOCC
        WORK_MATRIX(I,I,1,1) = 1.0_q
      ENDDO
    ENDIF

    ! copy the result to the correct potition 
    TRAFO_MATRIX(1:NOCC, 1:NOCC, IK, ISP) = WORK_MATRIX(:, :, IK, ISP)

    ! calc rest+entangled+local orbitals
    CALL LPRJ_LINCOM(W, TRAFO_MATRIX, W_EMBED)

    ! SCARY, fake the number of bands 
    ! (no guarantee that this will always work, but for OUTWAV it seems ok)
    W_EMBED%WDES%NB_TOT = NOCC + NEMBED_STATES
    CALL OUTWAV(IO, W_EMBED%WDES, W_EMBED, LATT_CUR, 0.0_q, 'ent')
    ! reset to "proper" number of bands
    W_EMBED%WDES%NB_TOT = W%WDES%NB_TOT


    !---------------
    ! third step: canonicalize the localized+entangled orbitals while
    ! freezing the rest orbitals. At the moment localized+entangled+rest are
    ! are stored in W_EMBED. Create a new structure (W_TEMP) to store the 
    ! canonicalized orbitals + rest.
    !---------------

    ! allocate the W_TEMP array and copy the entire W_EMBED to W_TEMP in
    ! order to canonicalize it.
    CALL ALLOCW(WDES_EMBED, W_TEMP)
    W_TEMP%CELTOT = W_EMBED%CELTOT
    W_TEMP%FERTOT = W_EMBED%FERTOT
    W_TEMP%CW     = W_EMBED%CW
    W_TEMP%CPROJ  = W_EMBED%CPROJ

    ! more subtle change the number of bands handled by ORTHCH and EDDIAG
    DO ISP = 1, W%WDES%ISPIN
      DO IK = 1, W%WDES%NKPTS ! KPOINTS_FULL%NKPTS
        W_EMBED%WDES%NB_TOTK(IK, ISP) = NINT(SUM(W%FERTOT(:, IK, ISP))) + NEMBED_STATES
      ENDDO
    ENDDO

    ! freeze the rest orbitals
    W_EMBED%WDES%NBANDSLOW = NREST + 1

    ! orthogonalize 
    CALL ORTHCH(W_EMBED%WDES, W_EMBED, INFO%LOVERL, SIZE(CQIJ,1), CQIJ)

    ! cononicalize localized + entangled orbitals
    CALL EDDIAG(HAMILTONIAN, GRID, LATT_CUR, NONLR_S, NONL_S, W_EMBED, &
                W_EMBED%WDES, SYMM, SIZE(CQIJ,1), CDIJ, CQIJ, 3, SV, &
                T_INFO, P, IO%IU0, E%EXHF)

    ! SCARY, fake the number of bands 
    W_EMBED%WDES%NB_TOT = NOCC + NEMBED_STATES
    CALL WRITE_EIGENVAL(W_EMBED%WDES, W_EMBED, IO%IU6)
    CALL OUTWAV(IO, W_EMBED%WDES, W_EMBED, LATT_CUR, 0.0_q, 'emb')


    !---------------
    ! fourth step: calculate the unitary transformation matrix (WORK_MATRIX)
    ! to transform between localized+entangled (W_EMBED) and the 
    ! canonicalized states (W_TEMP)
    !---------------

    DEALLOCATE(WORK_MATRIX)
    ALLOCATE( WORK_MATRIX(NOCC+NEMBED_STATES, NOCC+NEMBED_STATES, 1, 1) )
    CALL CALC_OVERLAP_MATRIX(W_TEMP, W_EMBED, CQIJ, NOCC, NEMBED_STATES, WORK_MATRIX)

    ! copy the relevant entries
    DEALLOCATE(TRAFO_MATRIX)
    ALLOCATE( TRAFO_MATRIX(NEMBED_STATES, 2*NEMBED_STATES, 1, 1) )
    TRAFO_MATRIX(:,:,1,1) = WORK_MATRIX(NOCC+1:NOCC+NEMBED_STATES, &
                                        NOCC-NEMBED_STATES+1:NOCC+NEMBED_STATES, 1, 1)

    ! calc projector onto localized (cluster) states P 
    ! TRAFO_MATRIX = TRAFO_MATRIX^H * TRAFO_MATRIX
    ALLOCATE( CLUSTER_PROJECTOR(2*NEMBED_STATES, 2*NEMBED_STATES, 1, 1) )
    CLUSTER_PROJECTOR = 0 ! nullify
    CALL GGEMM('C', 'N', 2*NEMBED_STATES, 2*NEMBED_STATES, NEMBED_STATES, one, &
               TRAFO_MATRIX, NEMBED_STATES, TRAFO_MATRIX, NEMBED_STATES, zero, &
               CLUSTER_PROJECTOR(:,:,1,1), 2*NEMBED_STATES)

    ! write out matrix
    CALL WRITE_GMATRIX(CLUSTER_PROJECTOR, 2*NEMBED_STATES, 2*NEMBED_STATES, 1, 1, "cluster_projector")


    !---------------
    ! final cleanup 
    !---------------
    CALL DEALLOCW(W_EMBED)
    CALL DEALLOCW(W_TEMP)

    DEALLOCATE(OVERLAPS)
    DEALLOCATE(CLUSTER_PROJECTOR)
    DEALLOCATE(WORK_MATRIX)
    DEALLOCATE(TRAFO_MATRIX)

#ifdef PROFILING
    CALL STOP_PROFILING('EMBED_SETUP_ORBITALS_M2')
#endif

#ifdef PROFILING
     ! started in INIT_PROFILING
     PROFILING_STOP('total_time')
     IF (IO%IU0>=0) THEN
        CALL PRINT_FULL_PROFILE(IO%IU6, PFLAT=.FALSE.)
     ENDIF
#endif
    CALL vtutor%stopCode()

    CONTAINS

      !**************** SUBROUTINE SETUP_ENTGLT_REST_TRAFO ********************
      ! this routine performs a SVD of the matrix OVERLAPS which corresponds
      ! to the overlap between occupied and localized orbitals.
      ! The result is a tranformation matrix (TRAFO_ENTGLT_REST) to calc
      ! the entangled and rest orbitals
      !**********************************************************************
      SUBROUTINE SETUP_ENTGLT_REST_TRAFO(OVERLAPS, TRAFO_ENTGLT_REST)
        GDEF :: OVERLAPS(:,:)
        GDEF :: TRAFO_ENTGLT_REST(:,:,:,:)
        ! local
        GDEF, ALLOCATABLE    :: U(:,:), VD(:,:), SIGMA(:)
        INTEGER              :: INFO, LWORK
        INTEGER              :: I
        GDEF                 :: ASKWORK
        INTEGER, ALLOCATABLE :: IWORK(:)
        GDEF, ALLOCATABLE    :: WORK(:)

        ! SVD: query optimal size of work arrays
        ALLOCATE( U(NOCC,NOCC), VD(NEMBED_STATES,NEMBED_STATES), SIGMA(NEMBED_STATES) )
        ALLOCATE( IWORK(8*MIN(NOCC,NEMBED_STATES)) )
#ifdef gammareal
        CALL DGESDD('A', NOCC, NEMBED_STATES, OVERLAPS, NOCC, SIGMA, &
                    U, NOCC, VD, NEMBED_STATES, ASKWORK, -1, IWORK, INFO)
#else
        CALL ZGESDD('A', NOCC, NEMBED_STATES, OVERLAPS, NOCC, SIGMA, &
                    U, NOCC, VD, NEMBED_STATES, ASKWORK, -1, IWORK, INFO)
#endif
        IF (INFO /= 0) THEN
           CALL vtutor%error("ERROR: GGESDD query returned: " // str(INFO))
        ENDIF

        LWORK = INT(ASKWORK)
        ALLOCATE(WORK(LWORK))

        ! SVD: perform OVERLAPS = U SIGMA VD
#ifdef gammareal
        CALL DGESDD('A', NOCC, NEMBED_STATES, OVERLAPS, NOCC, SIGMA, &
                    U, NOCC, VD, NEMBED_STATES, WORK, LWORK, IWORK, INFO)
#else
        CALL ZGESDD('A', NOCC, NEMBED_STATES, OVERLAPS, NOCC, SIGMA, &
                    U, NOCC, VD, NEMBED_STATES, WORK, LWORK, IWORK, INFO)
#endif
        IF (INFO /= 0) THEN
           CALL vtutor%error("ERROR: GGESDD query returned: " // str(INFO))
        ENDIF

        ! calc tranformation matrix for bath orbitals
        ! TRAFO_ENTGLT_REST(i,n) = sum_m^NEMBED_STATES U(i,m) * VD(m,n)
        CALL GGEMM('N', 'N', NOCC, NEMBED_STATES, NEMBED_STATES, one, &
                   U, NOCC, VD, NEMBED_STATES, zero, TRAFO_ENTGLT_REST, NOCC)

        ! reorder: shift entangled states to NOCC-NEMBED_STATES+1
        DO I = NEMBED_STATES, 1, -1 ! go backwards to not overwrite data
          TRAFO_ENTGLT_REST(:, I+NOCC-NEMBED_STATES, IK, ISP) = TRAFO_ENTGLT_REST(:, I, IK, ISP)
        ENDDO

        ! calc tranformation matrix for rest orbitals
        ! TRAFO_ENTGLT_REST(i,n-NEMBED_STATES) = U(i,n) for n = NEMBED_STATES+1,...,NOCC
        TRAFO_ENTGLT_REST(:, 1:NOCC-NEMBED_STATES, IK, ISP) = U(:, NEMBED_STATES+1:NOCC)

        DEALLOCATE(IWORK)
        DEALLOCATE(WORK)
        DEALLOCATE(SIGMA)
        DEALLOCATE(VD)
        DEALLOCATE(U)

      END SUBROUTINE SETUP_ENTGLT_REST_TRAFO


      !**************** SUBROUTINE CALC_OVERLAP_MATRIX ***********************
      ! calculates the matrix of the inner products 
      ! M(i,j,k,s) = < W1_i,k,s | S | W2_j,k,s >
      !**********************************************************************
      SUBROUTINE CALC_OVERLAP_MATRIX(W1, W2, CQIJ, NOCC, NEMBED_STATES, M)
        USE choleski, ONLY: CALC_PAW_OVERLAP
        USE wave_high

        TYPE(wavespin) :: W1
        TYPE(wavespin) :: W2
        OVERLAP        :: CQIJ(:, :, :, :)
        INTEGER        :: NOCC
        INTEGER        :: NEMBED_STATES
        GDEF           :: M(:, :, :, :)
        ! local
        TYPE(wavedes1), TARGET :: WDESK
        TYPE(wavefuna)         :: W1K, W2K
        INTEGER                :: ISPINOR, ISP, IK, I
        GDEF, ALLOCATABLE      ::  S(:,:)

        ALLOCATE(S(W1%WDES%NB_TOT, W2%WDES%NB_TOT))

        M = 0 ! nullify

        DO ISP = 1, W1%WDES%ISPIN
          DO IK = 1, W1%WDES%NKPTS 
             CALL SETWDES(W1%WDES, WDESK, IK)
             W1K = ELEMENTS(W1, WDESK, ISP)
             W2K = ELEMENTS(W2, WDESK, ISP)
             CALL CALC_PAW_OVERLAP(W1K, W2K, CQIJ, S)
             M(:, :, IK, ISP) = S(:,:)
          ENDDO
        ENDDO

        DEALLOCATE(S)

      END SUBROUTINE CALC_OVERLAP_MATRIX
  
  END SUBROUTINE EMBED_SETUP_ORBITALS_M2



  !************** SUBROUTINE EMBED_SETUP_ORBITALS_M4 ******************
  ! routine for method 4
  ! this routine simply localizes a set of bands and writes a WAVECAR
  !********************************************************************
  SUBROUTINE EMBED_SETUP_ORBITALS_M4(W, GRID, P, CQIJ, LATT_CUR, &
                                     T_INFO, INFO, IO)
    USE lattice
    USE pseudo
    USE poscar
    USE fileio
    USE nonl_high
    USE hamil_struct_def
    USE subrot
    USE choleski
#ifdef PROFILING
    USE profiling
#endif

    TYPE(wavespin)      :: W
    TYPE(grid_3d)       :: GRID
    TYPE(potcar)        :: P(:)
    OVERLAP             :: CQIJ(:, :, :, :)
    TYPE(latt)          :: LATT_CUR
    TYPE(type_info)     :: T_INFO
    TYPE(info_struct)   :: INFO
    TYPE(in_struct)     :: IO
    ! local
    GDEF, ALLOCATABLE :: TRAFO_MATRIX(:, :, :, :)
    TYPE(wavefuna)    :: WA
    TYPE(wavedes1)    :: WDESK
    INTEGER           :: NBSTART(2), NBSTOP(2), NBANDS(2) ! for spin up/down
    INTEGER           :: IK, ISP

    PROFILING_START("EMBED_SETUP_ORBITALS_M4")
    
    ! GAMMA ONLY and SPIN UNPOLARIZED
    IK = 1
    IF(W%WDES%NKPTS > 1) THEN
      CALL vtutor%error("EMBED_SETUP_ORBITALS_M4 only for one k-point")
      CALL vtutor%stopcode()
    ENDIF
        
    NBSTART = 0
    NBSTOP = 0
    NBANDS = 0
    DO ISP = 1, W%WDES%ISPIN
      NBSTART(ISP) = 1
      NBSTOP(ISP) = NINT(SUM(W%FERTOT(:, IK, ISP)))
      NBANDS(ISP) = NBSTOP(ISP) - NBSTART(ISP) + 1
    ENDDO

    ALLOCATE(TRAFO_MATRIX(MAXVAL(NBANDS),MAXVAL(NBANDS),IK,W%WDES%ISPIN))
    TRAFO_MATRIX = zero

    ! calc localization trafo via LOCPROJ and SVD
    CALL LPRJ_PROALL(W, W%WDES, P, CQIJ, LATT_CUR, LPRJ_functions, LPRJ_COVL, T_INFO, INFO, IO)
    CALL LPRJ_SVD_LOCALIZE(TRAFO_MATRIX, W%WDES%ISPIN, NBSTART, NBSTOP, IO)

    ! setup embedding orbitals using TRAFO_MATRIX

    DO ISP = 1, W%WDES%ISPIN
      CALL SETWDES(W%WDES, WDESK, IK)
      WA = ELEMENTS(W, WDESK, ISP)

      IF (W%WDES%DO_REDIS) THEN
         CALL REDISTRIBUTE_PROJ(WA)
         CALL REDISTRIBUTE_PW(WA)
      ENDIF
      
      WRITE(*,*) SIZE(TRAFO_MATRIX,1), SIZE(TRAFO_MATRIX,2),SIZE(TRAFO_MATRIX,3),SIZE(TRAFO_MATRIX,4)
      WRITE(*,*) NBANDS

      CALL LINCOM('F', &
                  WA%CW_RED(:,NBSTART(ISP):NBSTOP(ISP)), &
                  WA%CPROJ_RED(:,NBSTART(ISP):NBSTOP(ISP)), &
                  TRAFO_MATRIX(1:NBANDS(ISP),1:NBANDS(ISP),IK,ISP), &
                  NBANDS(ISP), &
                  NBANDS(ISP), &
                  WDESK%NPL_RED, &
                  WDESK%NPRO_RED, &
                  WDESK%NRPLWV_RED, &
                  WDESK%NPROD_RED, &
                  NBANDS(ISP), &
                  WA%CW_RED(:,NBSTART(ISP):NBSTOP(ISP)), &
                  WA%CPROJ_RED(:,NBSTART(ISP):NBSTOP(ISP)))
      
      IF (W%WDES%DO_REDIS) THEN
         CALL REDISTRIBUTE_PROJ(WA)
         CALL REDISTRIBUTE_PW(WA)
      ENDIF
    ENDDO

    CALL OUTWAV(IO, W%WDES, W, LATT_CUR, 0.0_q, 'emb')

    DEALLOCATE(TRAFO_MATRIX)

    PROFILING_STOP("EMBED_SETUP_ORBITALS_M4")
  
  ENDSUBROUTINE EMBED_SETUP_ORBITALS_M4



  !************** SUBROUTINE EMBED_DUMP_HAMIL ******************
  ! routine for method 5
  ! this routine calculates the matrix elements of the Hamiltonian
  ! in in the basis of the 
  !   i)  IAOs (intrinsic atomic orbitals) and
  !   ii) trial functions (LOCPROJ)
  !********************************************************************
  SUBROUTINE EMBED_DUMP_HAMIL(W, GRID, HAMILTONIAN, NONLR_S, NONL_S, &
                                     SYMM, P, CDIJ, CQIJ, SV, LATT_CUR, &
                                     T_INFO, INFO, IO)
    USE lattice
    USE pseudo
    USE poscar
    USE fileio
    USE nonl_high
    USE hamil_struct_def
    USE subrot
    USE choleski
#ifdef PROFILING
    USE profiling
#endif

    TYPE(wavespin)      :: W
    TYPE(grid_3d)       :: GRID
    TYPE (ham_handle)   :: HAMILTONIAN
    TYPE (nonlr_struct) :: NONLR_S
    TYPE (nonl_struct)  :: NONL_S
    TYPE (symmetry)     :: SYMM
    TYPE(potcar)        :: P(:)
    OVERLAP             :: CDIJ(:, :, :, :)
    OVERLAP             :: CQIJ(:, :, :, :)
    RGRID               :: SV(:,:) 
    !RGRID               :: SV(DIMREAL(GRID%MPLWV),W%WDES%NCDIJ) 
    TYPE(latt)          :: LATT_CUR
    TYPE(type_info)     :: T_INFO
    TYPE(info_struct)   :: INFO
    TYPE(in_struct)     :: IO
    ! local
    TYPE(wavespin)      :: W_BASIS
    GDEF, ALLOCATABLE   :: H(:,:,:,:) ! the hamiltonian matrix 
    INTEGER             :: LMDIM, IFLAG, I, J, RUN
    LOGICAL             :: LFIRST, LLAST
    REAL(q)             :: EXHF
    CHArACTER(LEN=3)    :: FPOSTFIX
    ! debug
    REAL(q), ALLOCATABLE :: EVLS(:)
    GDEF, ALLOCATABLE    :: GBUF(:)
    INTEGER              :: IERR, FUNIT, ISP, IK
    REAL(q)              :: POSI(3), POSJ(3), DIFF(3), DIST
    GDEF                 :: GDIM
    INTEGER              :: NMU

    PROFILING_START("EMBED_DUMP_HAMIL")

    IK = 1 !gamma only

    CALL LPRJ_PROALL(W, W%WDES, P, CQIJ, LATT_CUR, LPRJ_functions, LPRJ_COVL, T_INFO, INFO, IO)

    DO RUN = 1, 2

      SELECT CASE(RUN)
        CASE(1)
          !W_BASIS = LPRJ_GET_WIAO()
          FPOSTFIX = 'IAO'
        CASE(2)
          !W_BASIS = LPRJ_GET_WTRIALFUNC()
          WRITE(*,*) "LPRJ_GET_WTRIALFUNC DOWNT EXIST ANY MORE"
          STOP
          FPOSTFIX = 'TRL'
      END SELECT
  
      ALLOCATE( H(W_BASIS%WDES%NB_TOT, &
                  W_BASIS%WDES%NB_TOT, &
                  W_BASIS%WDES%NKPTS,  &
                  W_BASIS%WDES%ISPIN) )
      H = zero
  
      LMDIM = SIZE(CQIJ, 1)
  
      IFLAG = 1 ! only matrix elements, no diagonalization
      LFIRST = .FALSE.
      LLAST = .TRUE.
  
      ! calc matrix elements of Hamiltonian
      CALL EDDIAG(HAMILTONIAN,    &
                  GRID,           &
                  LATT_CUR,       &
                  NONLR_S,        &
                  NONL_S,         &
                  W_BASIS,        &
                  W_BASIS%WDES,   &
                  SYMM,           &
                  LMDIM,          &
                  CDIJ,           &
                  CQIJ,           &
                  IFLAG,          &
                  SV,             &
                  T_INFO,         &
                  P,              &
                  IO%IU0,         &
                  EXHF,           &
                  H(:,:,:,:),     &        
                  LFIRST,         &
                  LLAST)
  
      ! restrict to NMU since NB_TOT can be larger than then 
      ! number of actual IAOs (due to parallelization)
      NMU = LPRJ_GET_NUM_WAN()
  
      ! HACK: write out hamiltonian in IAO basis
      IF(IO%IU0 .GT. 0) THEN
        DO ISP = 1, W_BASIS%WDES%ISPIN
          FUNIT=700
          OPEN(unit = FUNIT,file = "hamil_" // FPOSTFIX // "basis_ISP" // str(ISP) // ".dat", FORM='FORMATTED', access='stream', STATUS='REPLACE')
          WRITE(FUNIT,'(I7)') NMU
          DO I = 1, NMU
            DO J = 1, NMU
              POSI = LPRJ_GET_POS(I)
              POSJ = LPRJ_GET_POS(J)
              DIFF = POSI - POSJ
              DIFF(:) = (DIFF(:) - NINT(DIFF(:))) * LATT_CUR%ANORM(:)
              DIST = SQRT(DIFF(1)**2.0+DIFF(2)**2.0+DIFF(3)**2.0)
#ifdef gammareal
              WRITE(FUNIT,'(2I7,E26.16,7F10.3)') I, J, H(I,J,IK,ISP), POSI, POSJ, DIST
#else
              WRITE(FUNIT,'(2I7,2E26.16,7F10.3)') I, J, REAL(H(I,J,IK,ISP), KIND=q), AIMAG(H(I,J,IK,ISP)), POSI, POSJ, DIST
#endif
            ENDDO
          ENDDO
          CLOSE(FUNIT)
        ENDDO
      ENDIF
  
      ! diagonalize to check eigenvalues
      ALLOCATE(EVLS(NMU))
      CALL DSYEV('V', 'L', NMU, H(1:NMU,1:NMU,1,1), NMU, EVLS(:), GDIM, -1, IERR) ! query
      IF (IERR .NE. 0) CALL vtutor%error("module LOCPROJ: LAPACK: Routine DSYEV/ZHEEV failed! INFO:" // str(IERR) )
      I = INT(REAL(GDIM, KIND=q))
      ALLOCATE(GBUF(I))
      CALL DSYEV('V', 'L', NMU, H(1:NMU,1:NMU,1,1), NMU, EVLS(:), GBUF, I, IERR) ! perform
  
      ! dump eigenvalues
      IF(IO%IU0 .GT. 0) THEN
        DO I = 1, NMU
          WRITE(IO%IU0,'(I6,3E20.10)') I, EVLS(I), REAL(W%CELTOT(I,1,1),KIND=q), REAL(EVLS(I)-W%CELTOT(I,1,1), KIND=q)
        ENDDO
        !DO I = 1, NMU
        !  IF(ABS(REAL(EVLS(I)-W%CELTOT(I,1,1), KIND=q)) .GT. 1E-4) THEN
        !    CALL vtutor%alert(" Found deviation of eigenvalues: " // str(ABS(REAL(EVLS(I)-W%CELTOT(I,1,1), KIND=q))) // " for i=" // str(I)  )
        !  ENDIF
        !ENDDO
      ENDIF
  
      DEALLOCATE(GBUF)
      DEALLOCATE(EVLS)
      DEALLOCATE(H)

    ENDDO

    PROFILING_STOP("EMBED_DUMP_HAMIL")
  
  ENDSUBROUTINE EMBED_DUMP_HAMIL




  !************** SUBROUTINE EMBED_REORDER_ORBITALS ******************
  ! routine for method 6
  ! this routine simply reorders the bands according to EMBED_SITES
  ! and the file orbitals_of_sites_ISP*
  !********************************************************************
  SUBROUTINE EMBED_REORDER_ORBITALS(W, T_INFO, LATT_CUR, IO)
    USE poscar
    USE m_mrgrnk
    USE dfast
    USE fileio
    USE wave_high
    TYPE (wavespin) :: W  
    TYPE(type_info) :: T_INFO
    TYPE(latt)      :: LATT_CUR
    TYPE(in_struct) :: IO
    ! local
    INTEGER              :: I, J, ISP, IK
    INTEGER              :: LOWEST_ONE, HIGHEST_ZERO
    INTEGER              :: IBMIN, IBMAX, ROTDIM
    INTEGER              :: FUNIT = 777
    INTEGER              :: IERR 
    INTEGER              :: VAL
    INTEGER              :: IATOM, NATOMS
    INTEGER, ALLOCATABLE :: SB(:)
    GDEF, ALLOCATABLE    :: ROTMAT(:,:)
    INTEGER, ALLOCATABLE :: NSORT(:)
    TYPE(wavefuna)       :: WA
    TYPE(wavedes1)       :: WDESK
    ! debug
    !GDEF, ALLOCATABLE :: MAT(:,:)

    PROFILING_START("EMBED_REORDER_ORBITALS")

    ! gamma-only at the moment
    IK = 1

    ! store total number of atoms (ions)
    NATOMS = T_INFO%NIONP

    DO ISP = 1, W%WDES%ISPIN

      ALLOCATE(SB(MAXVAL(W%WDES%NB_TOTK(IK,:))))
      SB = 0

      IF(IO%IU6 .GT. 0) THEN
        OPEN(unit = FUNIT, file = "orbitals_of_sites_ISP" // str(ISP), FORM='FORMATTED', access='stream', STATUS='OLD', IOSTAT=IERR)  
        IF(IERR .NE. 0) THEN
          CALL vtutor%error("internal error in routine EMBED_REORDER_ORBITALS: cannot open file orbitals_of_sites_ISP" // str(ISP) // ". Error: " // str(IERR))
          CALL vtutor%stopCode()
        ENDIF
      ENDIF

      IF(IO%IU6 .GT. 0) THEN
        IBMIN = W%WDES%NB_TOTK(IK,ISP)
        IBMAX = 1
        DO IATOM = 1, NATOMS
          DO 
            I = 0
            READ(FUNIT,'(I8)',IOSTAT=IERR,advance='no') I
            IF(I .EQ. 0) EXIT
            IF( ANY(EMBED_SITES .EQ. IATOM) ) THEN
              SB(I) = 1 ! mark band as selected
              IF(I .LT. IBMIN) IBMIN = I
              IF(I .GT. IBMAX) IBMAX = I
            ELSE
              IF(I .GT. IBMAX) IBMAX = I
            ENDIF
            IF(IERR .LT. 0) EXIT
          ENDDO
        ENDDO
      ELSE
        IBMIN = 0
        IBMAX = 0
      ENDIF
      CALLMPI(M_sum_i(W%WDES%COMM, SB, SIZE(SB(:))))
      CALLMPI(M_sum_i(W%WDES%COMM, IBMIN, 1))
      CALLMPI(M_sum_i(W%WDES%COMM, IBMAX, 1))
      !IBMIN = 1
      WRITE(*,*) "MIN, MAX", IBMIN, IBMAX

      !DO IATOM = 1, NATOMS
      !  ! check if this atom was selected with EMBED_SITES
      !  IF( ANY(EMBED_SITES .EQ. IATOM) ) THEN
      !    DO
      !      I = 0
      !      READ(FUNIT,'(I8)',IOSTAT=IERR,advance='no') I
      !      IF(I .EQ. 0) EXIT
      !      SB(I,ISP) = 1 ! mark band as selected
      !      IF(I .LT. IBMIN) IBMIN = I
      !      IF(I .GT. IBMAX) IBMAX = I
      !      IF(IERR .LT. 0) EXIT
      !    ENDDO
      !  ELSE 
      !    READ(FUNIT,*) ! dummy "read" to skip this line in file
      !  ENDIF
      !ENDDO
      !WRITE(*,*) "MIN, MAX", IBMIN, IBMAX

      !WRITE(*,*) ""
      !!DO I = IBMIN, IBMAX
      !DO I = 1, SIZE(SB(:,ISP))
      !  WRITE(*,*) SB(I,ISP)
      !ENDDO
      !WRITE(*,*) ""

      ! setup matrix to reorder bands
      ROTDIM = IBMAX-IBMIN+1
      ALLOCATE(NSORT(ROTDIM))
      CALL MRGRNK(SB(IBMIN:IBMAX), NSORT) ! sort
      ALLOCATE(ROTMAT(ROTDIM,ROTDIM))
      ROTMAT = zero
      DO I = 1, ROTDIM
        ROTMAT(NSORT(I),I) = one
      ENDDO

      !! ---------------------------------------------------
      !! setup matrix to reorder bands
      !! the strategy is to swap all zeros in SB(:,ISP) 
      !! with a lower lying one.
      !! ---------------------------------------------------
      !ROTDIM = IBMAX-IBMIN+1
      !ALLOCATE(ROTMAT(ROTDIM,ROTDIM))
      !! initialize as identity matrix
      !ROTMAT = zero
      !DO I = 1, ROTDIM
      !  ROTMAT(I,I) = one
      !ENDDO
      !!LOWEST_ONE = IBMIN   ! this is a "one" for sure
      !HIGHEST_ZERO = IBMAX ! just a guess 
      !DO I = IBMAX, IBMIN, -1
      !  IF(SB(I) .EQ. 0) THEN ! we found the highest "zero" => swap with lowest "one"
      !    HIGHEST_ZERO = I
      !    IF( (HIGHEST_ZERO-LOWEST_ONE) .LE. 1) EXIT
      !    ROTMAT(HIGHEST_ZERO-IBMIN+1,HIGHEST_ZERO-IBMIN+1) = zero
      !    ROTMAT(LOWEST_ONE  -IBMIN+1,LOWEST_ONE  -IBMIN+1) = zero
      !    ROTMAT(HIGHEST_ZERO-IBMIN+1,LOWEST_ONE  -IBMIN+1) = one
      !    ROTMAT(LOWEST_ONE  -IBMIN+1,HIGHEST_ZERO-IBMIN+1) = one
      !    DO J = LOWEST_ONE+1, IBMAX ! reset LOWEST_ONE
      !      IF(SB(J) .EQ. 1) THEN
      !        LOWEST_ONE = J
      !        EXIT
      !      ENDIF
      !    ENDDO
      !  ENDIF
      !ENDDO

      !DO I = 1, ROTDIM
      !  DO J = 1, ROTDIM
      !    WRITE(*,'(F5.1)',advance='no') ROTMAT(I,J)
      !  ENDDO
      !  WRITE(*,*) " "
      !ENDDO

      !ALLOCATE(MAT(4,ROTDIM))
      !MAT = zero
      !DO I = 1, ROTDIM
      !  MAT(:,I) = I*one
      !ENDDO
      !WRITE(*,*) " "
      !DO I = 1, 4
      !  DO J = 1, ROTDIM
      !    WRITE(*,'(F5.1)',advance='no') MAT(I,J)
      !  ENDDO
      !  WRITE(*,*) " "
      !ENDDO
      !MAT = MATMUL(MAT,ROTMAT)
      !WRITE(*,*) " "
      !DO I = 1, 4
      !  DO J = 1, ROTDIM
      !    WRITE(*,'(F5.1)',advance='no') MAT(I,J)
      !  ENDDO
      !  WRITE(*,*) " "
      !ENDDO

      !SB(IBMIN:IBMAX,ISP) = MATMUL(ROTMAT(:,:), SB(IBMIN:IBMAX,ISP))
      !WRITE(*,*) ""
      !DO I = IBMIN, IBMAX
      !  WRITE(*,*) SB(I,ISP)
      !ENDDO
      !WRITE(*,*) ""

      ! finally re-sort orbitals with LINCOM
      CALL SETWDES(W%WDES, WDESK, IK)
      WA = ELEMENTS(W, WDESK, ISP)
      IF (W%WDES%DO_REDIS) THEN
         CALL REDISTRIBUTE_PROJ(WA)
         CALL REDISTRIBUTE_PW(WA)
      ENDIF
      CALL LINCOM('F',                         &
                  WA%CW_RED(:,IBMIN:IBMAX),    &   
                  WA%CPROJ_RED(:,IBMIN:IBMAX), &
                  ROTMAT(:,:),                 & 
                  ROTDIM,                      &
                  ROTDIM,                      &
                  WDESK%NPL_RED,               &  
                  WDESK%NPRO_RED,              & 
                  WDESK%NRPLWV_RED,            & 
                  WDESK%NPROD_RED,             &
                  ROTDIM,                      &
                  WA%CW_RED(:,IBMIN:IBMAX),    &
                  WA%CPROJ_RED(:,IBMIN:IBMAX))  
      IF (W%WDES%DO_REDIS) THEN
         CALL REDISTRIBUTE_PROJ(WA)
         CALL REDISTRIBUTE_PW(WA)
      ENDIF
      DEALLOCATE(ROTMAT)
      DEALLOCATE(NSORT)

      ! write to EMBED_SITES and EMBED_STATES to OUTCAR
      IF(IO%IU6 .GE. 0) THEN
        WRITE(IO%IU6,'(/A,I2)') " for ISP = ", ISP
        WRITE(IO%IU6,'(A$)') " EMBED_SITES = "
        DO I = 1, SIZE(EMBED_SITES)
          WRITE(IO%IU6,'(I5$)') EMBED_SITES(I)
        ENDDO
        WRITE(IO%IU6,*)
        WRITE(IO%IU6,'(A,I5)') " #EMBED_SITES = ", SIZE(EMBED_SITES)
        WRITE(IO%IU6,'(A$)') " EMBED_STATES = "
        J=0
        DO I = 1, W%WDES%NB_TOTK(IK,ISP)
          IF(SB(I) .EQ. 1) THEN
            WRITE(IO%IU6,'(I5$)') I
            J=J+1
          ENDIF
        ENDDO
        WRITE(IO%IU6,*) 
        WRITE(IO%IU6,'(A,I5)') " #EMBED_STATES = ", J
      ENDIF

      DEALLOCATE(SB)

    ENDDO

    CALL OUTWAV(IO, W%WDES, W, LATT_CUR, 0.0_q, 'srt')

    PROFILING_STOP("EMBED_REORDER_ORBITALS")
                                   
  END SUBROUTINE EMBED_REORDER_ORBITALS
  SUBROUTINE EMBED_REORDER_ORBITALS_OLD(W, T_INFO, LATT_CUR, IO)
    USE poscar
    USE dfast
    USE fileio
    USE wave_high
    TYPE (wavespin) :: W  
    TYPE(type_info) :: T_INFO
    TYPE(latt)      :: LATT_CUR
    TYPE(in_struct) :: IO
    ! local
    INTEGER              :: I, J, ISP, IK
    INTEGER              :: LOWEST_ONE, HIGHEST_ZERO
    INTEGER              :: IBMIN, IBMAX, ROTDIM
    INTEGER              :: FUNIT = 777
    INTEGER              :: IERR 
    INTEGER              :: VAL
    INTEGER              :: IATOM, NATOMS
    INTEGER, ALLOCATABLE :: SB(:,:)
    GDEF, ALLOCATABLE    :: ROTMAT(:,:)
    TYPE(wavefuna)       :: WA
    TYPE(wavedes1)       :: WDESK
    ! debug
    !GDEF, ALLOCATABLE :: MAT(:,:)

    PROFILING_START("EMBED_REORDER_ORBITALS_OLD")

    ! gamma-only at the moment
    IK = 1

    ! store total number of atoms (ions)
    NATOMS = T_INFO%NIONP

    ALLOCATE(SB(MAXVAL(W%WDES%NB_TOTK(IK,:)),W%WDES%ISPIN))
    SB = 0

    DO ISP = 1, W%WDES%ISPIN

      IF(IO%IU6 .GT. 0) THEN
        OPEN(unit = FUNIT, file = "orbitals_of_sites_ISP" // str(ISP), FORM='FORMATTED', access='stream', STATUS='OLD', IOSTAT=IERR)  
        IF(IERR .NE. 0) THEN
          CALL vtutor%error("internal error in routine EMBED_REORDER_ORBITALS_OLD: cannot open file orbitals_of_sites_ISP" // str(ISP) // ". Error: " // str(IERR))
          CALL vtutor%stopCode()
        ENDIF
      ENDIF

      IF(IO%IU6 .GT. 0) THEN
        IBMIN = W%WDES%NB_TOTK(IK,ISP)
        IBMAX = 1
        DO IATOM = 1, NATOMS
          DO 
            I = 0
            READ(FUNIT,'(I8)',IOSTAT=IERR,advance='no') I
            IF(I .EQ. 0) EXIT
            IF( ANY(EMBED_SITES .EQ. IATOM) ) THEN
              SB(I,ISP) = 1 ! mark band as selected
              IF(I .LT. IBMIN) IBMIN = I
            ELSE
              IF(I .GT. IBMAX) IBMAX = I
            ENDIF
            IF(IERR .LT. 0) EXIT
          ENDDO
        ENDDO
      ELSE
        IBMIN = 0
        IBMAX = 0
      ENDIF
      CALLMPI(M_sum_i(W%WDES%COMM, SB, SIZE(SB(:,ISP))))
      CALLMPI(M_sum_i(W%WDES%COMM, IBMIN, 1))
      CALLMPI(M_sum_i(W%WDES%COMM, IBMAX, 1))
      WRITE(*,*) "MIN, MAX", IBMIN, IBMAX

      !DO IATOM = 1, NATOMS
      !  ! check if this atom was selected with EMBED_SITES
      !  IF( ANY(EMBED_SITES .EQ. IATOM) ) THEN
      !    DO
      !      I = 0
      !      READ(FUNIT,'(I8)',IOSTAT=IERR,advance='no') I
      !      IF(I .EQ. 0) EXIT
      !      SB(I,ISP) = 1 ! mark band as selected
      !      IF(I .LT. IBMIN) IBMIN = I
      !      IF(I .GT. IBMAX) IBMAX = I
      !      IF(IERR .LT. 0) EXIT
      !    ENDDO
      !  ELSE 
      !    READ(FUNIT,*) ! dummy "read" to skip this line in file
      !  ENDIF
      !ENDDO
      !WRITE(*,*) "MIN, MAX", IBMIN, IBMAX

      !WRITE(*,*) ""
      !!DO I = IBMIN, IBMAX
      !DO I = 1, SIZE(SB(:,ISP))
      !  WRITE(*,*) SB(I,ISP)
      !ENDDO
      !WRITE(*,*) ""

      ! ---------------------------------------------------
      ! setup matrix to reorder bands
      ! the strategy is to swap all zeros in SB(:,ISP) 
      ! with a lower lying one.
      ! ---------------------------------------------------
      ROTDIM = IBMAX-IBMIN+1
      ALLOCATE(ROTMAT(ROTDIM,ROTDIM))
      ! initialize as identity matrix
      ROTMAT = zero
      DO I = 1, ROTDIM
        ROTMAT(I,I) = one
      ENDDO
      LOWEST_ONE = IBMIN   ! this is a "one" for sure
      HIGHEST_ZERO = IBMAX ! just a guess 
      DO I = IBMAX, IBMIN, -1
        IF(SB(I,ISP) .EQ. 0) THEN ! we found the highest "zero" => swap with lowest "one"
          HIGHEST_ZERO = I
          IF( (HIGHEST_ZERO-LOWEST_ONE) .LE. 1) EXIT
          ROTMAT(HIGHEST_ZERO-IBMIN+1,HIGHEST_ZERO-IBMIN+1) = zero
          ROTMAT(LOWEST_ONE  -IBMIN+1,LOWEST_ONE  -IBMIN+1) = zero
          ROTMAT(HIGHEST_ZERO-IBMIN+1,LOWEST_ONE  -IBMIN+1) = one
          ROTMAT(LOWEST_ONE  -IBMIN+1,HIGHEST_ZERO-IBMIN+1) = one
          DO J = LOWEST_ONE+1, IBMAX ! reset LOWEST_ONE
            IF(SB(J,ISP) .EQ. 1) THEN
              LOWEST_ONE = J
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      !DO I = 1, ROTDIM
      !  DO J = 1, ROTDIM
      !    WRITE(*,'(F5.1)',advance='no') ROTMAT(I,J)
      !  ENDDO
      !  WRITE(*,*) " "
      !ENDDO

      !ALLOCATE(MAT(4,ROTDIM))
      !MAT = zero
      !DO I = 1, ROTDIM
      !  MAT(:,I) = I*one
      !ENDDO
      !WRITE(*,*) " "
      !DO I = 1, 4
      !  DO J = 1, ROTDIM
      !    WRITE(*,'(F5.1)',advance='no') MAT(I,J)
      !  ENDDO
      !  WRITE(*,*) " "
      !ENDDO
      !MAT = MATMUL(MAT,ROTMAT)
      !WRITE(*,*) " "
      !DO I = 1, 4
      !  DO J = 1, ROTDIM
      !    WRITE(*,'(F5.1)',advance='no') MAT(I,J)
      !  ENDDO
      !  WRITE(*,*) " "
      !ENDDO

      !SB(IBMIN:IBMAX,ISP) = MATMUL(ROTMAT(:,:), SB(IBMIN:IBMAX,ISP))
      !WRITE(*,*) ""
      !DO I = IBMIN, IBMAX
      !  WRITE(*,*) SB(I,ISP)
      !ENDDO
      !WRITE(*,*) ""

      ! finally re-sort orbitals with LINCOM
      CALL SETWDES(W%WDES, WDESK, IK)
      WA = ELEMENTS(W, WDESK, ISP)
      IF (W%WDES%DO_REDIS) THEN
         CALL REDISTRIBUTE_PROJ(WA)
         CALL REDISTRIBUTE_PW(WA)
      ENDIF
      CALL LINCOM('F',                         &
                  WA%CW_RED(:,IBMIN:IBMAX),    &   
                  WA%CPROJ_RED(:,IBMIN:IBMAX), &
                  ROTMAT(:,:),                 & 
                  ROTDIM,                      &
                  ROTDIM,                      &
                  WDESK%NPL_RED,               &  
                  WDESK%NPRO_RED,              & 
                  WDESK%NRPLWV_RED,            & 
                  WDESK%NPROD_RED,             &
                  ROTDIM,                      &
                  WA%CW_RED(:,IBMIN:IBMAX),    &
                  WA%CPROJ_RED(:,IBMIN:IBMAX))  
      IF (W%WDES%DO_REDIS) THEN
         CALL REDISTRIBUTE_PROJ(WA)
         CALL REDISTRIBUTE_PW(WA)
      ENDIF
      DEALLOCATE(ROTMAT)

      ! write to EMBED_SITES and EMBED_STATES to OUTCAR
      IF(IO%IU6 .GE. 0) THEN
        WRITE(IO%IU6,'(/A,I2)') " for ISP = ", ISP
        WRITE(IO%IU6,'(A$)') " EMBED_SITES = "
        DO I = 1, SIZE(EMBED_SITES)
          WRITE(IO%IU6,'(I5$)') EMBED_SITES(I)
        ENDDO
        WRITE(IO%IU6,*)
        WRITE(IO%IU6,'(A,I5)') " #EMBED_SITES = ", SIZE(EMBED_SITES)
        WRITE(IO%IU6,'(A$)') " EMBED_STATES = "
        J=0
        DO I = 1, W%WDES%NB_TOTK(IK,ISP)
          IF(SB(I,ISP) .EQ. 1) THEN
            WRITE(IO%IU6,'(I5$)') I
            J=J+1
          ENDIF
        ENDDO
        WRITE(IO%IU6,*) 
        WRITE(IO%IU6,'(A,I5)') " #EMBED_STATES = ", J
      ENDIF

    ENDDO

    CALL OUTWAV(IO, W%WDES, W, LATT_CUR, 0.0_q, 'srt')

    DEALLOCATE(SB)

    PROFILING_STOP("EMBED_REORDER_ORBITALS_OLD")
                                   
  END SUBROUTINE EMBED_REORDER_ORBITALS_OLD



  !***************** SUBROUTINE EMBED_BANDORDER_DISTANCE ********************
  !
  ! sorts the bands according to the distance to a point in 
  ! POSCAR (fractional coordinates).
  ! The closest band has highest band number.
  !
  !************************************************************************
  SUBROUTINE EMBED_BANDORDER_DISTANCE(W, IWF_LOW, IWF_HIGH, POS, &
                                      CENTROIDS, LATT_CUR, IO)
    USE m_mrgrnk
    USE wave_high
    USE dfast
    USE lattice
    USE fileio

    TYPE (wavespin)    :: W
    INTEGER            :: IWF_LOW
    INTEGER            :: IWF_HIGH
    REAL(q)            :: POS(3)
    REAL(q)            :: CENTROIDS(:,:)
    TYPE (latt)        :: LATT_CUR
    TYPE(in_struct)    :: IO
    ! local
    INTEGER              :: NWF
    INTEGER, ALLOCATABLE :: NSORT(:)
    GDEF, ALLOCATABLE    :: USORT(:,:) ! unitary matrix to re-sort orbitals
    REAL(q)              :: D1, D2, D3
    REAL(q), ALLOCATABLE :: DISTANCES(:)
    TYPE(wavedes1)       :: WDESK
    TYPE(wavefuna)       :: WA
    INTEGER              :: I, IK, ISP

    NWF = IWF_HIGH-IWF_LOW+1

    ALLOCATE(DISTANCES(NWF))

    ! calcualte distances
    DO I = 1, NWF
      D1 = CENTROIDS(1,I)-POS(1)
      D2 = CENTROIDS(2,I)-POS(2)
      D3 = CENTROIDS(3,I)-POS(3)

      ! minimum image convention
      D1 = (D1-NINT(D1))*LATT_CUR%ANORM(1)   
      D2 = (D2-NINT(D2))*LATT_CUR%ANORM(2) 
      D3 = (D3-NINT(D3))*LATT_CUR%ANORM(3) 

      DISTANCES(I) = SQRT( D1**2 + D2**2 + D3**2 )
    ENDDO

    ! sort the distances 
    ALLOCATE(NSORT(NWF))
    ! we take the negative distances here since 
    ! we want closest orbitals last 
    CALL MRGRNK(-DISTANCES, NSORT) 

    ! build sort unitary
    ALLOCATE(USORT(NWF,NWF))
    USORT = zero
    DO I = 1, NWF
      USORT(NSORT(I),I) = 1.0_q
    ENDDO

    ! re-sort distances for output
    DISTANCES = MATMUL(TRANSPOSE(USORT),DISTANCES)
    IF(IO%IU6 .GT. 0) THEN
      WRITE(*,*) " "
      WRITE(*,*) " sorted orbital distances to POS"
      DO I = 1, NWF
        WRITE(*,'(I6,F9.4)') I, DISTANCES(I)
      ENDDO
    ENDIF

    ! gamma only, spin unpolarized
    ISP = 1
    IK = 1

    ! re-sort orbitals via LINCOM

    CALL SETWDES(W%WDES, WDESK, IK)
    WA = ELEMENTS(W, WDESK, ISP)

    IF (W%WDES%DO_REDIS) THEN
       CALL REDISTRIBUTE_PROJ(WA)
       CALL REDISTRIBUTE_PW(WA)
    ENDIF

    CALL LINCOM('F', &
                WA%CW_RED(:,IWF_LOW:IWF_HIGH), &
                WA%CPROJ_RED(:,IWF_LOW:IWF_HIGH), &
                USORT(:,:), &
                NWF, &
                NWF, &
                WDESK%NPL_RED, &
                WDESK%NPRO_RED, &
                WDESK%NRPLWV_RED, &
                WDESK%NPROD_RED, &
                NWF, &
                WA%CW_RED(:,IWF_LOW:IWF_HIGH), &
                WA%CPROJ_RED(:,IWF_LOW:IWF_HIGH))
    
    IF (W%WDES%DO_REDIS) THEN
       CALL REDISTRIBUTE_PROJ(WA)
       CALL REDISTRIBUTE_PW(WA)
    ENDIF

    CALL OUTWAV(IO, W%WDES, W, LATT_CUR, 0.0_q, 'loc')

    DEALLOCATE(NSORT)
    DEALLOCATE(USORT)
    DEALLOCATE(DISTANCES)

  END SUBROUTINE EMBED_BANDORDER_DISTANCE



  !***************** SUBROUTINE EMBED_ORBITAL_CENTROIDS ********************
  !
  ! calculates the charge center of mass (usefull for localized orbitals)
  !
  !************************************************************************
  SUBROUTINE EMBED_ORBITAL_CENTROIDS(W, WDES, IWF_LOW, IWF_HIGH, GRID, GRIDC, &
                                     GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C,   &
                                     T_INFO, LATT_CUR, P, SYMM, INFO, LMDIM,  &
                                     N_MIX_PAW, IRDMAX, IO)
    USE poscar
    USE lattice
    USE pseudo
    USE us

    TYPE (wavespin)    :: W
    TYPE (wavedes)     :: WDES
    INTEGER            :: IWF_LOW(2)
    INTEGER            :: IWF_HIGH(2)
    TYPE (grid_3d)     :: GRID, GRIDC, GRID_SOFT, GRIDUS
    TYPE (transit)     :: C_TO_US, SOFT_TO_C 
    TYPE (type_info)   :: T_INFO
    TYPE (latt)        :: LATT_CUR
    TYPE (potcar)      :: P(T_INFO%NTYP)
    TYPE (symmetry)    :: SYMM
    TYPE (info_struct) :: INFO
    INTEGER            :: LMDIM, N_MIX_PAW, IRDMAX
    TYPE(in_struct)    :: IO
    !REAL(q)            :: CENTROIDS(:,:)
    ! local
    OVERLAP,    ALLOCATABLE :: CRHODE(:,:,:,:) ! occupancy of augmentation channels
    COMPLEX(q), ALLOCATABLE :: CHDEN(:,:)      ! soft pseudo charge
    COMPLEX(q), ALLOCATABLE :: CHTOT(:,:)      ! total charge
    REAL(q), ALLOCATABLE    :: RHOLM(:,:)   
    REAL(q), ALLOCATABLE    :: FERTOT(:,:,:) ! tmp backup of weights
    RGRID,  ALLOCATABLE     :: GRID_WORK(:)
    REAL(q),ALLOCATABLE     :: RHO_Z(:,:)
    INTEGER                 :: IK, ISP, NWF(2)
    INTEGER                 :: I, J, ITER, N, NGX, NGY, NGZ
    INTEGER                 :: NI, NX, NY, NZ
    REAL(q)                 :: DX, DY, DZ
    REAL(q)                 :: ORIGIN_SHIFT_NX, ORIGIN_SHIFT_NY, ORIGIN_SHIFT_NZ
    REAL(q)                 :: MRX, MRY, MRZ, CHANGE_CX, CHANGE_CY, CHANGE_CZ, STDDEV
    REAL(q)                 :: MAXDENS, GUESSED_MEAN_NX, GUESSED_MEAN_NY, GUESSED_MEAN_NZ 
    REAL(q), ALLOCATABLE    :: STD_DEVIATION(:)
    REAL(q), ALLOCATABLE    :: MEAN(:,:)
    LOGICAL                 :: CONVERGED
    REAL(q)                 :: RWORK
    REAL(q)                 :: POS(3)

    IF (GRIDC%NPLWV /= GRIDC%NGX*GRIDC%NGY*GRIDC%NGZ) THEN
      CALL vtutor%error("internal ERROR in EMBED_ORBITAL_CENTROIDS: NPLWV is &
                         not compatibel with NGX,NGY,NGZ \n " // &
                         str(GRIDC%NPLWV) // " " // str(GRIDC%NGX) // " " // & 
                         str(GRIDC%NGY) // " " // str(GRIDC%NGZ))
    ENDIF

    ! GAMMA ONLY, NO SPIN POLARIZATION
    IK = 1

    NGX = GRIDC%NGX
    NGY = GRIDC%NGY
    NGZ = GRIDC%NGZ

    NWF(:) = IWF_HIGH(:)-IWF_LOW(:)+1

    ALLOCATE(CRHODE(LMDIM, LMDIM, WDES%NIONS, WDES%NCDIJ))
    ALLOCATE(CHDEN(GRID_SOFT%MPLWV, WDES%NCDIJ))
    ALLOCATE(CHTOT(GRIDC%MPLWV, WDES%NCDIJ))
    ALLOCATE(RHOLM(N_MIX_PAW, WDES%NCDIJ))
    ALLOCATE(GRID_WORK(DIMREAL(GRIDC%MPLWV)))
    ALLOCATE(RHO_Z(NGX,NGY))

    !CALL SLEEP(2*W%WDES%COMM%NODE_ME)
    !WRITE(*,*) "NWF", NWF
    !WRITE(*,*) "GRID_SOFT%MPLWV", GRID_SOFT%MPLWV
    !WRITE(*,*) "GRIDC%MPLWV", GRIDC%MPLWV
    !WRITE(*,*) "GRIDC%RC%NP", GRIDC%RC%NP
    !WRITE(*,*) "WDES%NCDIJ", WDES%NCDIJ
    !WRITE(*,*) "dim", NGX*NGY*NGZ
    !WRITE(*,*) "GRIDC%RL%NCOL, GRIDC%RL%NROW", GRIDC%RL%NCOL, GRIDC%RL%NROW
    !WRITE(*,*) "LATT_CUR%A(1,:)", LATT_CUR%A(1,:)
    !WRITE(*,*) "LATT_CUR%A(2,:)", LATT_CUR%A(2,:)
    !WRITE(*,*) "LATT_CUR%A(3,:)", LATT_CUR%A(3,:)
    !WRITE(*,*) "LATT_CUR%ANORM(1)", LATT_CUR%ANORM(1)
    !WRITE(*,*) "LATT_CUR%ANORM(2)", LATT_CUR%ANORM(2)
    !WRITE(*,*) "LATT_CUR%ANORM(3)", LATT_CUR%ANORM(3)
    !WRITE(*,*) "GRIDC%RL%NFAST", GRIDC%RL%NFAST
    !WRITE(*,*) "NGX NGY NGZ", NGX, NGY, NGZ
    !STOP

    ! backup FERTOT 
    ALLOCATE(FERTOT(WDES%NB_TOT, WDES%NKPTS, WDES%ISPIN))
    FERTOT = W%FERTOT

    IK = 1

    DO ISP = 1, WDES%ISPIN
      ! output
      IF(IO%IU6 .GT. 0) THEN
        WRITE(*,'(A,I1,A,I1)') " Spin channel ", ISP, " of ", WDES%ISPIN
        WRITE(*,'(A)') " centroids in fractional coordinates (RX,RY,RZ) and standard deviation SIGMA in Angstrom"
        WRITE(*,'(A)') "  #band    RX       RY       RZ       SIGMA"
      ENDIF

      ALLOCATE(MEAN(3,NWF(ISP)))
      ALLOCATE(STD_DEVIATION(NWF(ISP)))
      MEAN = 0.0_q
      DO I = IWF_LOW(ISP), IWF_HIGH(ISP)
        W%FERTOT(:,:,:) = 0.0_q
        W%FERTOT(I,IK,ISP) = 1.0_q
        ! calcualte charge from single (local) orbital
        CALL SET_CHARGE(W, WDES, INFO%LOVERL, GRID, GRIDC, GRID_SOFT, GRIDUS, &
                        C_TO_US, SOFT_TO_C, LATT_CUR, P, SYMM, T_INFO, CHDEN, &
                        LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)
        ! transfer to GRID_WORK and FFT
        CALL RC_ADD(CHTOT, 1.0_q, CHTOT, 0.0_q, GRID_WORK, GRIDC)
        ! transform to real space representation
        CALL FFT3D(GRID_WORK, GRIDC, 1)

        ! find the position with maximum density as an estimate for the "position"
        ! of the orbital.
        MAXDENS = 0.0_q
        DO NZ = 1, NGZ ! loop over all planes 
          ! get slice of density and store in RHO_Z
          CALL MRG_GRID_RL_PLANE(GRIDC, RHO_Z, GRID_WORK, NZ) 
          DO NY = 1, NGY
            DO NX = 1, NGX
              IF(MAXDENS .LT. RHO_Z(NX,NY)) THEN
                GUESSED_MEAN_NX = 1.0_q*NX
                GUESSED_MEAN_NY = 1.0_q*NY
                GUESSED_MEAN_NZ = 1.0_q*NZ
                MAXDENS = RHO_Z(NX,NY)
              ENDIF
            ENDDO
          ENDDO
        ENDDO

        DO ITER=1, 30  ! if not conv. after 30 iter. it will never do

          ! now that we have a rough estimate for the "position" of the orbital, we
          ! shift the origin, such that this "position" is exactly in the center of
          ! the cell. We do this to avoid that sigificant parts of the orbital cut
          ! an edge of the integration volume.
          ORIGIN_SHIFT_NX = GUESSED_MEAN_NX - (1.0_q*NGX/2.0_q + 1.0_q)   
          ORIGIN_SHIFT_NY = GUESSED_MEAN_NY - (1.0_q*NGY/2.0_q + 1.0_q)
          ORIGIN_SHIFT_NZ = GUESSED_MEAN_NZ - (1.0_q*NGZ/2.0_q + 1.0_q)

          MRX = 0.0_q
          MRY = 0.0_q
          MRZ = 0.0_q
          STDDEV = 0.0_q
          CONVERGED = .FALSE.

          ! construct planes for every z
          DO NZ = 1, NGZ

            ! get slice of density and store in RHO_Z
            CALL MRG_GRID_RL_PLANE(GRIDC, RHO_Z, GRID_WORK, NZ) 
            
            ! normalize density to 1
            RHO_Z = RHO_Z / ( NGX*NGY*NGZ * ( 2.0_q - ABS(1-W%WDES%ISPIN) ) )
            
            DO NY = 1, NGY
              DO NX = 1, NGX

                ! calc distance between origin and current point (NX,NY,NZ)
                DX = 1.0_q*NX - (1.0_q + ORIGIN_SHIFT_NX)   
                DY = 1.0_q*NY - (1.0_q + ORIGIN_SHIFT_NY) 
                DZ = 1.0_q*NZ - (1.0_q + ORIGIN_SHIFT_NZ)

                IF(DX .LT. 0.0_q) DX = DX + 1.0_q*NGX
                IF(DY .LT. 0.0_q) DY = DY + 1.0_q*NGY
                IF(DZ .LT. 0.0_q) DZ = DZ + 1.0_q*NGZ
                IF(DX .GT. 1.0_q*NGX) DX = DX - 1.0_q*NGX
                IF(DY .GT. 1.0_q*NGY) DY = DY - 1.0_q*NGY
                IF(DZ .GT. 1.0_q*NGZ) DZ = DZ - 1.0_q*NGZ

                ! transform to relative distances
                DX = DX / (1.0_q*NGX)
                DY = DY / (1.0_q*NGY)
                DZ = DZ / (1.0_q*NGZ)

                ! UNNESSECARY (AND WRONG)
                !! calculate scaling factor to correctly weight the
                !! corners of the unit cell
                !F = 1.0_q
                !IF(NX .EQ. 1) F=F*2.0_q  
                !IF(NY .EQ. 1) F=F*2.0_q 
                !IF(NZ .EQ. 1) F=F*2.0_q 

                ! calc centroid <r> (mean value) 
                MRX = MRX + DX*RHO_Z(NX,NY)  
                MRY = MRY + DY*RHO_Z(NX,NY)
                MRZ = MRZ + DZ*RHO_Z(NX,NY)

                ! for calculation of variance we need to scale 
                DX = DX * LATT_CUR%ANORM(1)
                DY = DY * LATT_CUR%ANORM(2)
                DZ = DZ * LATT_CUR%ANORM(3)
                STDDEV = STDDEV + RHO_Z(NX,NY) * (DX**2 + DY**2 + DZ**2) 

              ENDDO 
            ENDDO 
          ENDDO 

          ! convert to standard deviation
          STDDEV = SQRT( STDDEV - ( &
                           (MRX*LATT_CUR%ANORM(1))**2 + &
                           (MRY*LATT_CUR%ANORM(2))**2 + &
                           (MRZ*LATT_CUR%ANORM(3))**2   &
                         )                              &
                       )
          
          ! shift back to standard origin
          MRX = MRX + ( ORIGIN_SHIFT_NX / (1.0_q*NGX) )
          MRY = MRY + ( ORIGIN_SHIFT_NY / (1.0_q*NGY) )
          MRZ = MRZ + ( ORIGIN_SHIFT_NZ / (1.0_q*NGZ) )

          ! calculate array index
          J = I+1 - IWF_LOW(ISP)

          ! calc the change of centroid (for convergence check)
          CHANGE_CX = ABS(MEAN(1,J)-MRX) 
          CHANGE_CY = ABS(MEAN(2,J)-MRY) 
          CHANGE_CZ = ABS(MEAN(3,J)-MRZ) 

          ! store
          MEAN(1,J) = MRX 
          MEAN(2,J) = MRY 
          MEAN(3,J) = MRZ 
          STD_DEVIATION(J) = STDDEV

          ! check convergence
          IF( (CHANGE_CX .LT. 0.01_q / NGX) .AND. &
              (CHANGE_CY .LT. 0.01_q / NGY) .AND. &
              (CHANGE_CZ .LT. 0.01_q / NGZ) ) THEN
            CONVERGED = .TRUE.
            EXIT
          ENDIF

          ! reset guess for centroid on (pseudo) integer grid
          GUESSED_MEAN_NX = MEAN(1,J)*NGX
          GUESSED_MEAN_NY = MEAN(2,J)*NGY
          GUESSED_MEAN_NZ = MEAN(3,J)*NGZ
        ENDDO 

        ! write out centroid and standard deviation
        IF(IO%IU6 .GT. 0) THEN
          IF(CONVERGED) THEN
            WRITE(*,'(A,I6,4F9.4)') " ", I, MEAN(1,J), MEAN(2,J), MEAN(3,J), STD_DEVIATION(J)
          ELSE
            WRITE(*,'(A,I6,4F9.4,A)') " ", I, MEAN(1,J), MEAN(2,J), MEAN(3,J), STD_DEVIATION(J), " (centroid search not converged, delocalized orbital?)"
          ENDIF
        ENDIF

      ENDDO

      IF(IO%IU6 .GT. 0) THEN
        WRITE(*,'(A,F9.4)') " average standard deviation per band: ", &
                            SUM(STD_DEVIATION(:))/(1.0_q*(NWF(ISP))) 
      ENDIF

      ! find maximal standard deviation
      RWORK = 0.0_q
      DO J = 1, NWF(ISP)
        IF(STD_DEVIATION(J) .GT. RWORK) RWORK = STD_DEVIATION(J)
      ENDDO
      IF(IO%IU6 .GT. 0) &
        WRITE(*,'(A,F9.4)') " maximal standard deviation: ", RWORK 
      ! find minimal standard deviation
      DO J = 1, NWF(ISP)
        IF(STD_DEVIATION(J) .LT. RWORK) RWORK = STD_DEVIATION(J)
      ENDDO
      IF(IO%IU6 .GT. 0) &
        WRITE(*,'(A,F9.4)') " minimal standard deviation: ", RWORK 

      DEALLOCATE(MEAN)
      DEALLOCATE(STD_DEVIATION)

    ENDDO ! spin

    !IF((SIZE(CENTROIDS,1) .NE. SIZE(MEAN,1)) &
    !    .OR. (SIZE(CENTROIDS,2) .NE. SIZE(MEAN,2))) THEN
    !  CALL vtutor%error("internal ERROR in EMBED_ORBITAL_CENTROIDS: dimension &
    !                     missmatch error.")
    !ENDIF
    !CENTROIDS = MEAN

    ! restore FERTOT
    W%FERTOT = FERTOT
    DEALLOCATE(FERTOT)

    !POS = [0.9699710620002729,  0.6490202557962372,  0.9872509517970919]
    POS = [0.206270635_q, 0.272414756_q, 0.500000000_q]
    !CALL EMBED_BANDORDER_DISTANCE(W, IWF_LOW, IWF_HIGH, POS, CENTROIDS)
    !CALL EMBED_BANDORDER_DISTANCE(W, IWF_LOW, IWF_HIGH, POS, MEAN, LATT_CUR, IO)

    DEALLOCATE(GRID_WORK, RHO_Z)
    DEALLOCATE(RHOLM, CHTOT, CHDEN, CRHODE)

  ENDSUBROUTINE EMBED_ORBITAL_CENTROIDS



  !***************** SUBROUTINE EMBED_CHARGE ********************
  ! This routine calculates local (cluster) charge density
  ! for method 1 as well as the local exchange correlation
  ! energy
  !***************************************************************
  SUBROUTINE EMBED_CHARGE(W, WDES, GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, &
                          SOFT_TO_C, T_INFO, CSTRF, DENCOR, LATT_CUR, P,    &
                          SYMM, INFO, LMDIM, N_MIX_PAW, IRDMAX, IU0, IU6)
    USE poscar
    USE lattice
    USE pseudo
    USE pawm
    USE us
    USE pot

    TYPE (wavespin)    :: W
    TYPE (wavedes)     :: WDES
    TYPE (grid_3d)     :: GRID, GRIDC, GRID_SOFT, GRIDUS
    TYPE (transit)     :: C_TO_US, SOFT_TO_C 
    TYPE (type_info)   :: T_INFO
    COMPLEX(q)         :: CSTRF(GRIDC%MPLWV, T_INFO%NTYP)
    RGRID              :: DENCOR(GRIDC%RL%NP)
    TYPE (latt)        :: LATT_CUR
    TYPE (potcar)      :: P(T_INFO%NTYP)
    TYPE (symmetry)    :: SYMM
    TYPE (info_struct) :: INFO
    INTEGER            :: LMDIM, N_MIX_PAW, IRDMAX
    INTEGER            :: IU0, IU6
    ! local
    OVERLAP,    ALLOCATABLE :: CRHODE(:,:,:,:) ! occupancy of augmentation channels
    COMPLEX(q), ALLOCATABLE :: CHDEN(:,:)      ! soft pseudo charge
    COMPLEX(q), ALLOCATABLE :: CHTOT(:,:)      ! total charge
    REAL(q), ALLOCATABLE    :: RHOLM(:,:)   
    REAL(q), ALLOCATABLE    :: FERTOT(:,:,:) ! tmp backup of weights
    COMPLEX(q), ALLOCATABLE :: CVTOT(:,:)
    RGRID, ALLOCATABLE      :: SV(:,:)
    OVERLAP, ALLOCATABLE    :: CDIJ(:,:,:,:)
    REAL(q)                 :: XCSIF(3,3)
    TYPE(energy)            :: E
    INTEGER :: NB_TOT
    INTEGER :: NOCC
    INTEGER :: NREST
    INTEGER :: IK, ISP

    ! deactivate core charges
    !DENCOR = 0.0_q

    ! GAMMA ONLY, NO SPIN POLARIZATION
    ISP = 1
    IK = 1

    NB_TOT = W%WDES%NB_TOT
    NOCC = NINT(SUM(W%FERTOT(:, IK, ISP)))
    NREST = MAX(0, NOCC - NEMBED_STATES)

    ALLOCATE(CRHODE(LMDIM, LMDIM, WDES%NIONS, WDES%NCDIJ))
    ALLOCATE(CHDEN(GRID_SOFT%MPLWV, WDES%NCDIJ))
    ALLOCATE(CHTOT(GRIDC%MPLWV, WDES%NCDIJ))
    ALLOCATE(RHOLM(N_MIX_PAW, WDES%NCDIJ))
    ALLOCATE(CVTOT(GRIDC%MPLWV, WDES%NCDIJ))
    ALLOCATE(SV(DIMREAL(GRID%MPLWV), WDES%NCDIJ))
    ALLOCATE(CDIJ(LMDIM, LMDIM, WDES%NIONS, WDES%NCDIJ))


    !!!!!!!!!!!!!!!!!!!!!! ALL !!!!!!!!!!!!!!!!!!!!!!!!!!!

    CALL SET_CHARGE(W, WDES, INFO%LOVERL, GRID, GRIDC, GRID_SOFT, GRIDUS, &
                    C_TO_US, SOFT_TO_C, LATT_CUR, P, SYMM, T_INFO, CHDEN, &
                    LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

    CALL POTLOK(GRID, GRIDC, GRID_SOFT, WDES%COMM_INTER, WDES, INFO, P,  &
                T_INFO, E, LATT_CUR, CHTOT, CSTRF, CVTOT, DENCOR, &
                SV, SOFT_TO_C, XCSIF)

    CALL SET_DD_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
                    CDIJ(1,1,1,1), RHOLM, CRHODE(1,1,1,1), E,        &
                    LMETA=.FALSE., LASPH=INFO%LASPH, LCOREL=.FALSE.)

    IF(IU6 > 0) THEN
      WRITE(IU6,*) "  "
      WRITE(IU6,*) "FOR ALL:"
      WRITE(IU6,*) "ALL: E%EXCG ", E%EXCG
      WRITE(IU6,*) "ALL: E%XCENC ", E%XCENC
      WRITE(IU6,*) "ALL: E%PAWPS ", E%PAWPS
      WRITE(IU6,*) "ALL: E%PAWAE ", E%PAWAE
      WRITE(IU6,*) "ALL: E%EXCG + E%PAWPSG + E%PAWAEG ", E%EXCG + E%PAWPSG + E%PAWAEG 
    ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    ! backup FERTOT and fake the weights such that only the  
    ! active space (i.e. local/cluster) is occupied
    ALLOCATE(FERTOT(WDES%NB_TOT, WDES%NKPTS, WDES%ISPIN))
    FERTOT = W%FERTOT


    !!!!!!!!!!!!!!!!!!!!!! BATH !!!!!!!!!!!!!!!!!!!!!!!!!!!
    W%FERTOT(NREST+1:WDES%NB_TOT,:,:) = 0

    CALL SET_CHARGE(W, WDES, INFO%LOVERL, GRID, GRIDC, GRID_SOFT, GRIDUS, &
                    C_TO_US, SOFT_TO_C, LATT_CUR, P, SYMM, T_INFO, CHDEN, &
                    LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

    CALL POTLOK(GRID, GRIDC, GRID_SOFT, WDES%COMM_INTER, WDES, INFO, P,  &
                T_INFO, E, LATT_CUR, CHTOT, CSTRF, CVTOT, DENCOR, &
                SV, SOFT_TO_C, XCSIF)

    CALL SET_DD_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
                    CDIJ(1,1,1,1), RHOLM, CRHODE(1,1,1,1), E,        &
                    LMETA=.FALSE., LASPH=INFO%LASPH, LCOREL=.FALSE.)

    IF(IU6 > 0) THEN
      WRITE(IU6,*) "  "
      WRITE(IU6,*) "FOR BATH:"
      WRITE(IU6,*) "BATH: E%EXCG ", E%EXCG
      WRITE(IU6,*) "BATH: E%XCENC ", E%XCENC
      WRITE(IU6,*) "BATH: E%PAWPS ", E%PAWPS
      WRITE(IU6,*) "BATH: E%PAWAE ", E%PAWAE
      WRITE(IU6,*) "BATH: E%EXCG + E%PAWPSG + E%PAWAEG ", E%EXCG + E%PAWPSG + E%PAWAEG 
    ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! restore FERTOT
    W%FERTOT = FERTOT


    !!!!!!!!!!!!!!!!!!!!!! ACTIVE !!!!!!!!!!!!!!!!!!!!!!!!!!!
    W%FERTOT(1:NREST,:,:) = 0

    CALL SET_CHARGE(W, WDES, INFO%LOVERL, GRID, GRIDC, GRID_SOFT, GRIDUS, &
                    C_TO_US, SOFT_TO_C, LATT_CUR, P, SYMM, T_INFO, CHDEN, &
                    LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

    CALL POTLOK(GRID, GRIDC, GRID_SOFT, WDES%COMM_INTER, WDES, INFO, P,  &
                T_INFO, E, LATT_CUR, CHTOT, CSTRF, CVTOT, DENCOR, &
                SV, SOFT_TO_C, XCSIF)

    CALL SET_DD_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
                    CDIJ(1,1,1,1), RHOLM, CRHODE(1,1,1,1), E,        &
                    LMETA=.FALSE., LASPH=INFO%LASPH, LCOREL=.FALSE.)

    IF(IU6 > 0) THEN
      WRITE(IU6,*) "  "
      WRITE(IU6,*) "FOR ACTIVE:"
      WRITE(IU6,*) "ACTIVE: E%EXCG ", E%EXCG
      WRITE(IU6,*) "ACTIVE: E%XCENC ", E%XCENC
      WRITE(IU6,*) "ACTIVE: E%PAWPS ", E%PAWPS
      WRITE(IU6,*) "ACTIVE: E%PAWAE ", E%PAWAE
      WRITE(IU6,*) "ACTIVE: E%EXCG + E%PAWPSG + E%PAWAEG ", E%EXCG + E%PAWPSG + E%PAWAEG 
    ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    ! restore FERTOT
    W%FERTOT = FERTOT
    DEALLOCATE(FERTOT)

    DEALLOCATE(CDIJ, SV, CVTOT, RHOLM, CHTOT, CHDEN, CRHODE)

  ENDSUBROUTINE EMBED_CHARGE



  !************************ SUBROUTINE WRITE_GMATRIX ***************************
  ! writes out a rectengular matrix in cc4s format
  ! GMAT(i,j)
  !****************************************************************************
  SUBROUTINE WRITE_CC4S_GMATRIX(GMAT, NROWS, NCOLS, FILENAME, THRESHOLD)
    GDEF    :: GMAT(:,:)
    INTEGER :: NROWS
    INTEGER :: NCOLS
    INTEGER :: NKPTS
    INTEGER :: ISPIN
    CHARACTER(LEN = *) :: FILENAME
    REAL(q), OPTIONAL :: THRESHOLD
    ! local
    INTEGER, PARAMETER :: FUNIT = 173
    INTEGER :: I, J, ISP, IK
    CHARACTER(LEN = 13) :: GFORMAT
    GDEF :: MATRIX_ELEMENT
    REAL(q) :: THRESH

#ifdef gammareal
    GFORMAT = '(1E26.16)'
#else
    GFORMAT = '(2E26.16)'
#endif
    
    ! if desired set threshold, otherwise set to -1.0
    IF(PRESENT(THRESHOLD)) THEN
      THRESH = THRESHOLD
    ELSE
      THRESH = -1.0
    ENDIF

    ! open file
    OPEN(UNIT=FUNIT, FILE=FILENAME, STATUS='REPLACE')

    ! write matrix dimension 
    WRITE(FUNIT,'(A,3I6)') "gmatrix", 2,  NROWS, NCOLS
    WRITE(FUNIT,'(A)') "ij "

    ! write data
    DO J = 1, NCOLS
      DO I = 1, NROWS
        MATRIX_ELEMENT = GMAT(I,J)
        ! skip if matrix element is smaller than threshold
        IF(ABS(MATRIX_ELEMENT) < THRESH) CYCLE
        WRITE(FUNIT, GFORMAT) MATRIX_ELEMENT
      ENDDO
    ENDDO

    CLOSE(FUNIT)
  END SUBROUTINE WRITE_CC4S_GMATRIX



  !************************ SUBROUTINE WRITE_GMATRIX ***************************
  ! writes out a rectengular matrix with information about spin k-points
  ! GMAT(i,j,k,s)
  !****************************************************************************
  SUBROUTINE WRITE_GMATRIX(GMAT, NROWS, NCOLS, ISPIN, NKPTS, FILENAME, THRESHOLD)
    GDEF    :: GMAT(:,:,:,:)
    INTEGER :: NROWS
    INTEGER :: NCOLS
    INTEGER :: NKPTS
    INTEGER :: ISPIN
    CHARACTER(LEN = *) :: FILENAME
    REAL(q), OPTIONAL :: THRESHOLD
    ! local
    INTEGER, PARAMETER :: FUNIT = 173
    INTEGER :: I, J, ISP, IK
    CHARACTER(LEN = 13) :: GFORMAT
    GDEF :: MATRIX_ELEMENT
    REAL(q) :: THRESH

#ifdef gammareal
    GFORMAT = '(4I6,1E26.16)'
#else
    GFORMAT = '(4I6,2E26.16)'
#endif
    
    ! if desired set threshold, otherwise set to -1.0
    IF(PRESENT(THRESHOLD)) THEN
      THRESH = THRESHOLD
    ELSE
      THRESH = -1.0
    ENDIF

    ! open file
    OPEN(UNIT=FUNIT, FILE=FILENAME, STATUS='REPLACE')

    ! write matrix dimension, spin, and k-points info
    WRITE(FUNIT,*) NROWS, NCOLS, ISPIN, NKPTS 

    ! write data
    DO ISP = 1, ISPIN
      Do IK = 1, NKPTS
        DO I = 1, NROWS
          DO J = 1, NCOLS
            MATRIX_ELEMENT = GMAT(I,J,IK,ISP)
            ! skip if matrix element is smaller than threshold
            IF(ABS(MATRIX_ELEMENT) < THRESH) CYCLE
            WRITE(FUNIT, GFORMAT) I, J , IK, ISP, MATRIX_ELEMENT
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    CLOSE(FUNIT)
  END SUBROUTINE WRITE_GMATRIX



  !************************ SUBROUTINE READ_GMATRIX ***************************
  ! reads a rectengular matrix with information about spin k-points
  ! GMAT(i,j,k,s)
  !
  ! if you want to read the matrix use
  ! CALL READ_GMATRIX(GMAT, "filename")
  !
  ! if you want to read the dimensions in the file, use
  ! CALL READ_GMATRIX(GMAT, "filename", QUERY_DIM=.TRUE., NROWS, NCOLS, ISPIN, NKPTS)
  !****************************************************************************
  SUBROUTINE READ_GMATRIX(GMAT, FILENAME, QUERY_DIM, NROWS, NCOLS, ISPIN, NKPTS)
    USE string, ONLY: str
    USE tutor, ONLY: vtutor
    GDEF    :: GMAT(:,:,:,:)
    CHARACTER(LEN = *) :: FILENAME
    LOGICAL, OPTIONAL :: QUERY_DIM
    INTEGER, OPTIONAL :: NROWS
    INTEGER, OPTIONAL :: NCOLS
    INTEGER, OPTIONAL :: NKPTS
    INTEGER, OPTIONAL :: ISPIN
    ! local
    INTEGER, PARAMETER :: FUNIT = 173
    INTEGER :: STAT
    INTEGER :: DIM1, DIM2, DIM3, DIM4
    INTEGER :: I, J, ISP, IK
    CHARACTER(LEN = 13) :: GFORMAT
    GDEF :: MATRIX_ELEMENT

#ifdef gammareal
    GFORMAT = '(4I6,1F26.16)'
#else
    GFORMAT = '(4I6,2F26.16)'
#endif

    ! open file
    OPEN(UNIT=FUNIT, FILE=FILENAME, STATUS='OLD', IOSTAT=STAT)

    IF(STAT .NE. 0) THEN
      CALL vtutor%error("ERROR IN READ_GMATRIX: cannot open file " // FILENAME)
    ENDIF
        
    READ(FUNIT,*,IOSTAT=STAT) DIM1, DIM2, DIM3, DIM4
    IF(STAT .NE. 0) THEN
      CALL vtutor%error("ERROR IN READ_GMATRIX: corrupted file " // FILENAME)
    ENDIF

    ! only return dimensions if required
    IF(PRESENT(QUERY_DIM)) THEN
      IF(QUERY_DIM) THEN
        IF(PRESENT(NROWS)) NROWS = DIM1
        IF(PRESENT(NCOLS)) NCOLS = DIM2
        IF(PRESENT(ISPIN)) ISPIN = DIM3
        IF(PRESENT(NKPTS)) NKPTS = DIM4
        RETURN
      ENDIF
    ENDIF

    ! check the dimensions to avoid segfault
    IF( (SIZE(GMAT,1) .LT. DIM1) .OR. (SIZE(GMAT,2) .LT. DIM2) .OR. &
        (SIZE(GMAT,3) .LT. DIM3) .OR. (SIZE(GMAT,4) .LT. DIM4) ) THEN
      CALL vtutor%error("ERROR IN READ_GMATRIX: missmatch of dimensions.")
    ENDIF

    GMAT = 0 ! nullify

    ! read the matrix until end-of-file (i.e. STAT < 0)
    DO
      READ(FUNIT, GFORMAT, IOSTAT=STAT) I, J , IK, ISP, MATRIX_ELEMENT
      IF(STAT .GT. 0) THEN
        CALL vtutor%error("ERROR IN READ_GMATRIX: ERROR CODE: " // str(STAT))
      ENDIF
      IF(STAT .LT. 0) EXIT
      GMAT(I, J, IK, ISP) = MATRIX_ELEMENT
    ENDDO

    CLOSE(FUNIT)
  END SUBROUTINE READ_GMATRIX






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!













!  !*****************************************************
!  ! OLD VERSION (Georg)
!  !*****************************************************
!  SUBROUTINE EMBED_LOCALIZATION_SELECT_ACTIVE(W, NTARGET_STATES, GRID, P, CQIJ, LATT_CUR, T_INFO, INFO, IO, & 
!       HAMILTONIAN,NONLR_S,NONL_S,SYMM,CDIJ,SV)
!    USE base
!    USE wave
!    USE lattice
!    USE locproj
!    USE pseudo
!    USE poscar
!    USE fileio
!    USE nonl_high
!    USE hamil_struct_def
!    USE subrot
!    USE choleski
!
!    TYPE(wavespin) :: W
!    INTEGER        :: NTARGET_STATES(:)
!    TYPE(grid_3d) :: GRID
!    TYPE(potcar) :: P(:)
!    OVERLAP :: CQIJ(:, :, :, :)
!    TYPE(latt) :: LATT_CUR
!    TYPE(type_info) :: T_INFO
!    TYPE(info_struct) :: INFO
!    TYPE(in_struct) :: IO
!  ! 
!    TYPE (ham_handle)  HAMILTONIAN
!    TYPE (nonlr_struct) NONLR_S
!    TYPE (nonl_struct) NONL_S
!    TYPE (symmetry) :: SYMM
!    OVERLAP :: CDIJ(:, :, :, :)
!    RGRID   SV(DIMREAL(GRID%MPLWV),W%WDES%NCDIJ) ! local potential
!  ! local
!    GDEF, ALLOCATABLE :: UK(:, :, :, :)
!    GDEF, ALLOCATABLE :: UK_TRAFO(:, :, :, :)
!    TYPE(wavedes)  :: WDES_WAN
!    TYPE(wavespin) :: W_WAN, W_LOCAL
!    INTEGER :: NTARGETS
!    INTEGER :: ISP, NK, NOCC
!    TYPE (energy)      E
!!-----------------------------------------------------------------------
!! first step calculate the local orbitals by projecting onto 
!! all Bloch states
!!-----------------------------------------------------------------------
!    ! W_WAN is allocted to store up to NB_TOT bands
!    WDES_WAN=W%WDES
!    CALL ALLOCW(WDES_WAN, W_WAN)
!    W_WAN%CELTOT=W%CELTOT
!    W_WAN%FERTOT=W%FERTOT
!
!    CALL LPRJ_PROALL(W, W%WDES, GRID, P, CQIJ, LATT_CUR, T_INFO, INFO, IO)
!    
!    CALL LPRJ_CALC_WANNIER_TRAFO(W, UK, IO)
!
!    CALL LPRJ_LINCOM(W, UK, W_WAN)
!
!    ! write local projected orbitals to file
!    CALL OUTWAV(IO, W_WAN%WDES, W_WAN, LATT_CUR, 0.0_q, 'wan')
!!-----------------------------------------------------------------------
!! now calculate the entangled orbitals by projection the selected
!! active local orbitals onto the occupied manyfold
!!-----------------------------------------------------------------------
!
!    ! first deallocate the LPRJ_COVL array and allocate it so that we can store the entangled orbitals
!    DEALLOCATE( LPRJ_COVL)
!    NTARGETS=SIZE(NTARGET_STATES)
!    ALLOCATE(LPRJ_COVL(W%WDES%NB_TOT, W%WDES%NKPTS, MAX(W%WDES%ISPIN,W%WDES%NRSPINORS), NTARGETS ))
!
!    ALLOCATE(UK_TRAFO(SIZE(UK,1), NTARGETS, SIZE(UK,3), SIZE(UK,4)))
!    ! save the transformation matrix from Bloch to local orbitals in UK_TRAFO
!    ! store transformation matrix from Blocj to local orbitals in LPRJ_COVL
!    CALL COPY_OVERLAP(W, W_WAN, UK, UK_TRAFO, NTARGET_STATES)
!
!    CALL LPRJ_CALC_WANNIER_TRAFO_3(W, UK, IO)
!
!    ! copy transformation matrix from UK_TRAFO back to UK (for states beyond NOCC)
!    CALL COPY_FROM_UK_TRAFO( W, UK_TRAFO, UK, NTARGET_STATES)
!
!    CALL LPRJ_LINCOM(W, UK, W_WAN)
!
!    ! scarry, fake the number of bands (no guarantee that this will always work, but for OUTWAV it seems ok)
!    W_WAN%WDES%NB_TOT=NINT(SUM(W%FERTOT(:, 1, 1)))+NTARGETS
!    CALL OUTWAV(IO, W_WAN%WDES, W_WAN, LATT_CUR, 0.0_q, 'ent')
!    ! reset to "proper" number of bands
!    W_WAN%WDES%NB_TOT=W%WDES%NB_TOT
!
!    ! now allocate the W_LOCAL array and copy the entire W_WAN to W_LOCAL
!    CALL ALLOCW(WDES_WAN, W_LOCAL)
!    W_LOCAL%CELTOT=W%CELTOT
!    W_LOCAL%FERTOT=W%FERTOT
!    W_LOCAL%CW    =W_WAN%CW
!    W_LOCAL%CPROJ =W_WAN%CPROJ
!
!    ! more subtle change the number of bands handled by ORTHCH and EDDIAG
!    DO ISP=1,W%WDES%ISPIN
!       DO NK=1,W%WDES%NKPTS ! KPOINTS_FULL%NKPTS
!          NOCC  = NINT(SUM(W%FERTOT(:, NK, ISP)))
!          W_WAN%WDES%NB_TOTK(NK,ISP)=NOCC+NTARGETS
!       ENDDO
!    ENDDO
!    W_WAN%WDES%NBANDSLOW=NOCC-NTARGETS+1
!
!    CALL ORTHCH(W_WAN%WDES,W_WAN, INFO%LOVERL, SIZE(CQIJ,1),CQIJ)
!
!    CALL EDDIAG(HAMILTONIAN,GRID,LATT_CUR,NONLR_S,NONL_S,W_WAN,W_WAN%WDES,SYMM, &
!         SIZE(CQIJ,1),CDIJ,CQIJ, 3, SV,T_INFO,P,IO%IU0,E%EXHF)
!
!    W_WAN%WDES%NB_TOT=NINT(SUM(W%FERTOT(:, 1, 1)))+NTARGETS
!    CALL WRITE_EIGENVAL( W_WAN%WDES, W_WAN, IO%IU6)
!    ! reset to "proper" number of bands
!    W_WAN%WDES%NB_TOT=W%WDES%NB_TOT
!
!    CALL  CALC_OVERLAP(W_WAN, W_LOCAL, UK, NOCC, NTARGETS)
!
!    IF (ALLOCATED(UK)) DEALLOCATE(UK)
!    IF (ALLOCATED(UK_TRAFO)) DEALLOCATE(UK_TRAFO)
!    CALL DEALLOCW(W_WAN)
!
!    CONTAINS
!
!!******************** SUBROUTINE COPY_OVERLAP **************************
!!      
!! calculate the overlap between two wave function arrays W and W_WAN
!! results is stored in the matrix UK (for each k-point and spin)
!! this routine is obviously not written for speed and essentially
!! a copy from  COPY_OVERLAP_GN in locproj.F
!! so suboptimal performance is guaranteed
!!
!!***********************************************************************
!
!      SUBROUTINE COPY_OVERLAP(W, W_WAN, UK, UK_TRAFO, NTARGET_STATES)
!        USE elphon
!        USE wave_high
!        IMPLICIT NONE
!        TYPE(wavespin) :: W
!        TYPE(wavespin) :: W_WAN
!        GDEF :: UK(:, :, :, :)
!        GDEF :: UK_TRAFO(:, :, :, :)
!        INTEGER        :: NTARGET_STATES(:)
!
!        TYPE(wavedes1), TARGET :: WDESK
!        TYPE(wavefuna) WK, W_WANK
!        INTEGER ISPINOR, ISP, NK, I, NOCC
!
!        DO ISPINOR=1,W%WDES%NRSPINORS
!        DO ISP=1,W%WDES%ISPIN
!
!        ! TODO: this works only if no symmetry is used
!        DO NK=1,W%WDES%NKPTS ! KPOINTS_FULL%NKPTS
!           ! descriptor to work on the present k-point
!           CALL SETWDES(W%WDES,WDESK,NK)
!
!           WK=ELEMENTS(W, WDESK, ISP)
!           W_WANK=ELEMENTS(W_WAN, WDESK, ISP)
!
!           NOCC  = NINT(SUM(W%FERTOT(:, NK, ISP)))
!
!           ! copy the results over to LPRJ_COVL
!           DO I=1,SIZE( NTARGET_STATES)
!              LPRJ_COVL(1:W%WDES%NB_TOT,NK,ISP+ISPINOR-1, I) = UK(1:W%WDES%NB_TOT, NTARGET_STATES(I), NK, ISP)
!
!           ENDDO
!           ! save the results over to UK_TRAFO
!           DO I=1,SIZE( NTARGET_STATES)
!              ! shortcut of UK_TRAFO is too small
!              UK_TRAFO(1:W%WDES%NB_TOT, I, NK, ISP)  = UK(1:W%WDES%NB_TOT, NTARGET_STATES(I), NK, ISP)
!           ENDDO
!        ENDDO
!        ENDDO
!        ENDDO
!
!      END SUBROUTINE COPY_OVERLAP
!
!
!      SUBROUTINE CALC_OVERLAP(W, W_LOCAL, UK, NOCC, NTARGETS)
!        USE elphon
!        USE wave_high
!        IMPLICIT NONE
!        TYPE(wavespin) :: W
!        TYPE(wavespin) :: W_LOCAL
!        GDEF :: UK(:, :, :, :)
!        INTEGER        :: NOCC,NTARGETS
!
!        TYPE(wavedes1), TARGET :: WDESK
!        TYPE(wavefuna) WK, W_LOCALK
!        INTEGER ISPINOR, ISP, NK, I
!
!      ! TODO: the overlap matrix S is identical to UK, so need to recalculate it
!      !  let's keep this for further testing for the time being
!        GDEF, ALLOCATABLE ::  S(:,:)
!
!        ALLOCATE(S(W%WDES%NB_TOT, W%WDES%NB_TOT))
!
!        DO ISPINOR=1,W%WDES%NRSPINORS
!        DO ISP=1,W%WDES%ISPIN
!
!        ! TODO: this works only if no symmetry is used
!        DO NK=1,W%WDES%NKPTS ! KPOINTS_FULL%NKPTS
!           ! descriptor to work on the present k-point
!           CALL SETWDES(W%WDES,WDESK,NK)
!
!           WK=ELEMENTS(W, WDESK, ISP)
!           W_LOCALK=ELEMENTS(W_LOCAL, WDESK, ISP)
!
!           CALL ELPH_SOVL(WK, W_LOCALK, SIZE(CQIJ,1) , CQIJ, S)
!
!           ! copy the results over to LPRJ_COVL
!           DO I=NOCC+1,NOCC+NTARGETS
!              WRITE(*,'(18F7.3)') S(NOCC-NTARGETS+1: NOCC+NTARGETS , I)
!              WRITE(*,'(1F18.14)') SUM(S(NOCC-NTARGETS+1:NOCC+NTARGETS , I )*S(NOCC-NTARGETS+1:NOCC+NTARGETS , I))
!           ENDDO
!        ENDDO
!        ENDDO
!        ENDDO
!
!        DEALLOCATE(S)
!
!      END SUBROUTINE CALC_OVERLAP
!
!
!
!      SUBROUTINE COPY_FROM_UK_TRAFO(W, UK_TRAFO, UK, NTARGET_STATES )
!        USE elphon
!        USE wave_high
!        IMPLICIT NONE
!        TYPE(wavespin) :: W
!        GDEF :: UK(:, :, :, :)
!        GDEF :: UK_TRAFO(:, :, :, :)
!        INTEGER        :: NTARGET_STATES(:)
!        INTEGER ISP, NK, I, NOCC
!
!
!        DO ISP=1,W%WDES%ISPIN
!        ! TODO: this works only if no symmetry is used
!        DO NK=1,W%WDES%NKPTS ! KPOINTS_FULL%NKPTS
!
!           NOCC  = NINT(SUM(W%FERTOT(:, NK, ISP)))
!
!           ! copy from UK_TRAFO to UK
!           DO I=1,SIZE( NTARGET_STATES)
!              IF (NOCC+I <= SIZE(UK,2)) THEN
!                 UK(1:W%WDES%NB_TOT, NOCC+I, NK, ISP)  = UK_TRAFO(1:W%WDES%NB_TOT, I, NK, ISP)
!              ENDIF
!           ENDDO
!        ENDDO
!        ENDDO
!
!      END SUBROUTINE COPY_FROM_UK_TRAFO
!  
!  END SUBROUTINE EMBED_LOCALIZATION_SELECT_ACTIVE

ENDMODULE embed

