#include "symbol.inc"

#if defined(gammareal) && defined(RPAgamma)
#undef RPAgamma
#endif
!*********************************************************************
!
! basis module for imaginary time GW implemenation 
! (chi_GG.F, chi_super.F)
! content:
!  -> forward and backward cosine transformation 
!  -> communicator routines
!  -> calculation of chi's head and wings  
!  -> screened Coulomb potential W 
!  -> contraction of GG to form CHI
!  -> contraction of GW to form SIGMA 
!
!*********************************************************************

MODULE GG_base
  USE prec
  USE poscar
  USE lattice
  USE mpimy
  USE nonl_high
  USE chi_glb
  USE chi_base
  USE xi
  USE wave_high
  USE mlr_optic
  USE ini, ONLY: START_TIMING, STOP_TIMING
  USE minimax_struct, ONLY : imag_grid_handle
  USE minimax, ONLY: LOCAL_INDEX_TO_GLOBAL
  IMPLICIT NONE
  PRIVATE
!
! size of arrays related to Greens functions
!
#ifdef gammareal
    REAL(q), PARAMETER :: wsgf=8._q
#else
    REAL(q), PARAMETER :: wsgf=16._q
#endif

  ! add contributions that grow exponentially with tau to the Green's function
  ! this can cause numerical issues, since these contributions can
  ! be huge at large tau compared to other contributions
  ! by means of this we can exactly implement what is done in chi.F
  ! namely the IP polarizability 
  !   X(iw) = f_n (1-f_m) [1/(e_n - e_m - iw) + 1/(e_n - e_m + iw)] 
  ! becomes after FT
  !   X(i tau) = f_n exp( tau  e_n) * f_m exp( tau  e_m)
  ! this allows for much faster convergence with respect to k-points
  ! assuming that the user has used values of 0.1 for sigma this yields
  ! essentially exact results (but even 0.2 works well)
  REAL(q) :: EFERMIADD=0.2_q

!
! this special flag allows to calculate the HF energy using sigma=G(r,r) v(r,r')
! this is done using one time point at tau=0 (last time point)
! adds some extra storage requirements at this point but allows better load balancing
!
  LOGICAL, PUBLIC :: LHFCALC_GG=.FALSE.

#ifdef scaLAPACK
  PUBLIC :: GDES_SET

! formaly, acont_
  PUBLIC :: greens_mat_des
  PUBLIC :: GREENS_MAT_DES_INIT
  PUBLIC :: ALLOCATE_GREENS_MAT_GDEF
  PUBLIC :: ALLOCATE_GREENS_MAT
  PUBLIC :: DUMP_GREENS_HAM_GDEF
  PUBLIC :: DUMP_GREENS_HAM
  PUBLIC :: DEALLOCATE_GREENS_MAT_GDEF
  PUBLIC :: DEALLOCATE_GREENS_MAT
  PUBLIC :: CHECK_GDES_MAT_SIZE
  PUBLIC :: DESTROY_GDES_MAT

  PUBLIC :: TRANS_TIME_FREQUENCY
  PUBLIC :: INTEGRATE_CHI_TAU
  PUBLIC :: COMPUTE_SCREENED_POTENTIAL
  PUBLIC :: DETERMINE_CURRENT_FREQ_GROUP
  PUBLIC :: DETERMINE_LOCAL_FREQ_INDEX
  PUBLIC :: DETERMINE_NTAU_GLOBAL
  PUBLIC :: CALCULATE_LONG_WAVE_LIMIT
  PUBLIC :: ALLOCATE_HEAD_HANDLE
  PUBLIC :: DEALLOCATE_HEAD_HANDLE
  PUBLIC :: SUBTRACT_SELF_CORRELATION
  PUBLIC :: DISTRIBUTE_HEAD_TO_CHI
  PUBLIC :: wsgf
  PUBLIC :: EFERMIADD
  PUBLIC :: DUMP_DIAGRAMMATIC_ENERGIES
  PUBLIC :: DUMP_DIAGRAMMATIC_ENERGIES_RPA
  PUBLIC :: SET_RESPONSE_FREQUENCIES
  PUBLIC :: DUMP_HEAD_GG
  PUBLIC :: ADD_DRUDE_TERM
  PUBLIC :: STORAGE_REQ_GG
  PUBLIC :: SET_NTAUPAR_NOMEGAPAR
  PUBLIC :: SET_LOOP_DES
  PUBLIC :: DETERMINE_BAND_GAP_SPIN
  PUBLIC :: SETUP_TAU_INDICES
  PUBLIC :: DUMP_PROGRESS
  PUBLIC :: DUMPX_FREQ
! formaly, acont
  PUBLIC :: FT_G_OR_SIGMA
  PUBLIC :: INV_FT_G_OR_SIGMA
  PUBLIC :: WRITE_SELF_ENERGY
  PUBLIC :: ADD_G_SIGMA_G
  PUBLIC :: ADD_G_SIGMA_Gxw
  PUBLIC :: ADD_G_SIGMA_G_new
  PUBLIC :: ADD_G_SIGMA_Gxw_new
  PUBLIC :: GF_WEIGHT
  PUBLIC :: ADD_G_SIGMA
  PUBLIC :: SET_G_SIGMA_DIAG
  PUBLIC :: GALITSKII_MIGDAL_DIAG
  PUBLIC :: FT_TEST_GF
 ! OEP related 
  PUBLIC :: OEP_POTENTIAL
  PUBLIC :: POTENTIAL_AND_CHARGE
  PUBLIC :: DEALLOCATE_OEP_HANDLE
! supercell 
  PUBLIC :: INIT_SUPERCELL
  PUBLIC :: INIT_SUPERCELL_FFT
  PUBLIC :: DEALLOCATE_SUPERCELL
!
! customize if required (see also scala.F, however, one might want to use
! a different blocking here)
!
  INTEGER,PRIVATE, SAVE :: NB_SCALA=16     ! blocking factor for distribution of matrices

! some BLACS definitions (see scala.F)
  INTEGER,PRIVATE, PARAMETER :: MB_=5
  INTEGER,PRIVATE, PARAMETER :: NB_=6 
  INTEGER,PRIVATE, PARAMETER :: DLEN_=9
  INTEGER,PRIVATE, PARAMETER :: M_=3
  INTEGER,PRIVATE, PARAMETER :: N_=4

  TYPE greens_mat_des
     TYPE (communic) :: COMM             ! communicator for communication between all nodes
     TYPE (communic) :: COMM_INTRA       ! communicator for nodes sharing one time/frequency
     TYPE (communic) :: COMM_INTER       ! communicator between time/frequency points
     INTEGER :: NROWS                    ! maximal matrix size (required for allocation)
     INTEGER :: MY_NROWS                 ! local row matrix size
     INTEGER :: MY_NCOLS                 ! local column matrix size 
     INTEGER :: ICTXT                    ! BLACS context defined on top of COMM_INTAU
     INTEGER :: DESC(DLEN_)              ! matrix descriptor
!
! customize if required (see also scala.F, however, one might want to use
! a different blocking here)
!
  END TYPE greens_mat_des

CONTAINS

!***********************************************************************
!
!  initialize the descriptor for the Green's function
!  the projector part is essentially the same routine as in wave.F
!   WDES_SET_NPRO
!  but since inheritance does not exist in F90, and because
!  subtle differences exist (different communicator) we need 
!  recoding
!
!***********************************************************************

  SUBROUTINE GDES_SET(GDES, WDES, WDES_RESPONSE, T_INFO, P, COMM)
    USE  wave
    USE  poscar
    USE  pseudo_struct_def

    TYPE (greensfdes)  GDES
    TYPE (wavedes)     WDES           ! descriptor for orbitals
    TYPE (wavedes)     WDES_RESPONSE  ! descriptor for response function
       ! the descriptor for the response function
       ! usually uses less data points in reciprocal space
    TYPE (type_info) :: T_INFO
    TYPE (potcar)   P(T_INFO%NTYP)
  ! local varibles
    INTEGER NALLOC,NPRO_ROW,NLM_ROW,NT,NI,NIS,NPRO, &
         LMMAXC,NIONS,LASTTYP, NI_L, NLM_COL, NODE_ME
    TYPE (communic) COMM

    GDES%COMM=COMM
    GDES%NRSPINORS=WDES%NRSPINORS
    GDES%LNONCOLLINEAR=WDES%LNONCOLLINEAR

    IF (GDES%LNONCOLLINEAR) THEN
       CALL vtutor%bug("internal ERROR in GDES_SET_NPRO: non collinear calculations not supported", __FILE__, __LINE__)
    ENDIF
!=======================================================================
!
! set up descriptor for plane wave coefficients of orbitals
!
!=======================================================================
    GDES%NRPLWV_ROW=WDES%NRPLWV
#ifndef MPI
    GDES%NRPLWV_COL=GDES%NRPLWV_ROW
    GDES%NRPLWV_POS=1
#else
    GDES%NRPLWV_COL=(GDES%NRPLWV_ROW+GDES%COMM%NCPU-1)/GDES%COMM%NCPU
    IF (GDES%NRPLWV_COL*GDES%COMM%NCPU /= GDES%NRPLWV_ROW) THEN
       CALL vtutor%bug("internal error in GDES_SET: NRPLWV not dividable by cores " // &
          str(GDES%NRPLWV_COL) // " " // str(GDES%COMM%NCPU) // " " // str(GDES%NRPLWV_ROW), __FILE__, __LINE__)
    ENDIF
    GDES%NRPLWV_POS=(GDES%COMM%NODE_ME-1)*GDES%NRPLWV_COL+1
    GDES%NRPLWV_COL=MIN(GDES%COMM%NODE_ME*GDES%NRPLWV_COL,GDES%NRPLWV_ROW)-GDES%NRPLWV_POS+1

#endif
    ! real to complex FFT is used
    GDES%NRPLWV_COL_DATA_POINTS=GDES%NRPLWV_COL
    GDES%NRPLWV_ROW_DATA_POINTS=GDES%NRPLWV_ROW
#ifdef gammareal
    GDES%NRPLWV_COL_DATA_POINTS=GDES%NRPLWV_COL*2
    GDES%NRPLWV_ROW_DATA_POINTS=GDES%NRPLWV_ROW*2  
#endif

    GDES%NRPLWV_COL_MAX_DATA_POINTS=GDES%NRPLWV_COL_DATA_POINTS
    CALLMPI( M_max_i(COMM, GDES%NRPLWV_COL_MAX_DATA_POINTS ,1))
!=======================================================================
!
! set up descriptor for plane wave coefficients of response function
!
!=======================================================================
    GDES%RES_NRPLWV_ROW=WDES_RESPONSE%NRPLWV

#ifndef MPI
    GDES%RES_NRPLWV_COL=GDES%RES_NRPLWV_ROW
    GDES%RES_NRPLWV_POS=1
#else
    GDES%RES_NRPLWV_COL=(GDES%RES_NRPLWV_ROW+GDES%COMM%NCPU-1)/GDES%COMM%NCPU
    IF (GDES%RES_NRPLWV_COL*GDES%COMM%NCPU /= GDES%RES_NRPLWV_ROW) THEN
       CALL vtutor%bug("internal error in GDES_SET: NRPLWV not dividable by cores " // &
          str(GDES%RES_NRPLWV_COL) // " " // str(GDES%COMM%NCPU) // " " // str(GDES%RES_NRPLWV_ROW), __FILE__, __LINE__)
    ENDIF
    GDES%RES_NRPLWV_POS=MIN((GDES%COMM%NODE_ME-1)*GDES%RES_NRPLWV_COL+1,GDES%RES_NRPLWV_ROW+1)
    GDES%RES_NRPLWV_COL=MIN(GDES%COMM%NODE_ME*GDES%RES_NRPLWV_COL,GDES%RES_NRPLWV_ROW)-GDES%RES_NRPLWV_POS+1
#endif
    ! real to complex FFT is used
    GDES%RES_NRPLWV_COL_DATA_POINTS=GDES%RES_NRPLWV_COL
    GDES%RES_NRPLWV_ROW_DATA_POINTS=GDES%RES_NRPLWV_ROW
#ifdef gammareal
    GDES%RES_NRPLWV_COL_DATA_POINTS=GDES%RES_NRPLWV_COL*2
    GDES%RES_NRPLWV_ROW_DATA_POINTS=GDES%RES_NRPLWV_ROW*2  
#endif

    GDES%RES_NRPLWV_COL_MAX_DATA_POINTS=GDES%RES_NRPLWV_COL_DATA_POINTS
    CALLMPI( M_max_i(COMM, GDES%RES_NRPLWV_COL_MAX_DATA_POINTS ,1))

    NALLOC=GDES%RES_NRPLWV_COL
    CALLMPI( M_sum_i(COMM, NALLOC ,1))
    IF (NALLOC/= GDES%RES_NRPLWV_ROW) THEN
       CALL vtutor%bug("internal ERROR in GDES_SET_NPRO NRPLWV_ROW: " // str(NALLOC) // " " // &
          str(GDES%RES_NRPLWV_ROW), __FILE__, __LINE__)
    ENDIF

!=======================================================================
!
!  real space presentation
!
!=======================================================================

    GDES%MPLWV_ROW=WDES%GRID%NPLWV
#ifndef MPI
    GDES%MPLWV_COL=GDES%MPLWV_ROW
    GDES%MPLWV_COL_MAX=GDES%MPLWV_ROW
    GDES%MPLWV_POS=1
#else
    GDES%MPLWV_COL=(GDES%MPLWV_ROW+GDES%COMM%NCPU-1)/GDES%COMM%NCPU
    GDES%MPLWV_POS=MIN((GDES%COMM%NODE_ME-1)*GDES%MPLWV_COL+1,GDES%MPLWV_ROW+1)
    GDES%MPLWV_COL=MIN(GDES%COMM%NODE_ME*GDES%MPLWV_COL,GDES%MPLWV_ROW)-GDES%MPLWV_POS+1

    GDES%MPLWV_COL_MAX=GDES%MPLWV_COL
    CALLMPI( M_max_i(COMM, GDES%MPLWV_COL_MAX ,1))
#endif
    NALLOC=GDES%MPLWV_COL
    CALLMPI( M_sum_i(COMM, NALLOC ,1))
    IF (NALLOC/= GDES%MPLWV_ROW) THEN
       CALL vtutor%bug("internal ERROR in GDES_SET_NPRO MPLWV_ROW: " // str(NALLOC) // " " // &
          str(GDES%MPLWV_ROW), __FILE__, __LINE__)
    ENDIF
!=======================================================================
!
!  finally most complicated bit
!  projected wave function character 
!
!=======================================================================
#ifndef MPI
    GDES%NIONS = T_INFO%NIONS
    ALLOCATE(GDES%NIONS_CORE(1))
    GDES%NIONS_CORE(1)=T_INFO%NIONS
    GDES%NTYP  = T_INFO%NTYP
    GDES%NITYP =>T_INFO%NITYP
    GDES%POSION=>T_INFO%POSION
    ALLOCATE(GDES%LMMAX(GDES%NTYP))
    ALLOCATE(GDES%NLM_LMMAX(GDES%NTYP))
    ALLOCATE(GDES%NT_GLOBAL(GDES%NTYP))
    ALLOCATE(GDES%NI_GLOBAL(T_INFO%NIONS))

!$  GDES%NT_LOCAL=>T_INFO%ITYP

    DO NT=1,T_INFO%NTYP
       GDES%LMMAX(NT)=P(NT)%LMMAX
       GDES%NLM_LMMAX(NT)=AUG_DES%LMMAX
       GDES%NT_GLOBAL(NT)=NT
    ENDDO
    DO NI=1,T_INFO%NIONS
       GDES%NI_GLOBAL(NI)=NI
    ENDDO
    GDES%NPRO_COL=SUM(GDES%LMMAX*GDES%NITYP)
    GDES%NPRO_COL_MAX=GDES%NPRO_COL
    GDES%NPRO_ROW=SUM(GDES%LMMAX*GDES%NITYP)
    NULLIFY(GDES%NPRO_POS)
    NULLIFY(GDES%NPRO_ENTRIES)

    GDES%NLM_COL=SUM(AUG_DES%LMMAX*GDES%NITYP)
    GDES%NLM_COL_MAX=SUM(AUG_DES%LMMAX*GDES%NITYP)
    GDES%NLM_ROW=SUM(AUG_DES%LMMAX*GDES%NITYP)
    NULLIFY(GDES%NLM_POS)
    NULLIFY(GDES%NLM_ENTRIES)

#else
    !-----------------------------------------------------------------------
    ! parallel version
    !-----------------------------------------------------------------------
    ! first count number of projection operators
    COMM = GDES%COMM
    NPRO_ROW=0
    NLM_ROW =0
    DO NT=1,T_INFO%NTYP
       NPRO_ROW=NPRO_ROW+P(NT)%LMMAX*T_INFO%NITYP(NT)
       NLM_ROW=NLM_ROW+ AUG_DES%LMMAX(NT)*T_INFO%NITYP(NT)
    ENDDO
    ! check implementation of NI_LOCAL and NI_GLOBAL
    DO NI=1,T_INFO%NIONS
       NI_L=NI_LOCAL(NI,COMM)
       IF (NI_L /= 0) THEN
          IF (NI /= NI_GLOBAL(NI_L,COMM)) THEN
             CALL vtutor%bug("internal error NI_GLOBAL, NI_LOCAL", __FILE__, __LINE__)
          ENDIF
       ENDIF
    ENDDO
    ! number of projection operators per node
    NALLOC=(T_INFO%NIONS+COMM%NCPU-1)/COMM%NCPU

    GDES%NIONS=0
    GDES%NTYP =0
    GDES%NPRO_COL =0
    GDES%NPRO_ROW=NPRO_ROW
    GDES%NLM_COL=0
    GDES%NLM_ROW =NLM_ROW

    ALLOCATE(GDES%NITYP(T_INFO%NTYP))
    ALLOCATE(GDES%NT_GLOBAL(T_INFO%NTYP))
    ALLOCATE(GDES%NI_GLOBAL(T_INFO%NIONS))
    ALLOCATE(GDES%LMMAX(T_INFO%NTYP))
    ALLOCATE(GDES%NPRO_POS(NALLOC,COMM%NCPU))
    ALLOCATE(GDES%NPRO_ENTRIES(NALLOC,COMM%NCPU))
    ALLOCATE(GDES%NLM_POS(NALLOC,COMM%NCPU))
    ALLOCATE(GDES%NLM_ENTRIES(NALLOC,COMM%NCPU))
    ALLOCATE(GDES%NLM_LMMAX(T_INFO%NTYP))
    ALLOCATE(GDES%POSION(3,NALLOC))
!$  ALLOCATE(GDES%NT_LOCAL(T_INFO%NIONS))

    GDES%LMMAX=0
    GDES%NITYP=0
    GDES%NPRO_POS=0
    GDES%NPRO_ENTRIES=0
    GDES%NLM_POS=0
    GDES%NLM_ENTRIES=0

    LASTTYP=0
    NIS    =1

    NPRO=0
    NLM_ROW=0

    NODE_ME=GDES%COMM%NODE_ME

    DO NT=1,T_INFO%NTYP
       LMMAXC=P(NT)%LMMAX
!       IF (LMMAXC==0) GOTO 100
       DO NI=NIS,T_INFO%NITYP(NT)+NIS-1
          ! does this element reside on local node
          IF (NI_LOCAL(NI,COMM) /=0 ) THEN
             GDES%NIONS=GDES%NIONS+1
             GDES%POSION(:,GDES%NIONS)=T_INFO%POSION(:,NI)
             GDES%NPRO_POS(GDES%NIONS,NODE_ME)=NPRO
             GDES%NPRO_ENTRIES(GDES%NIONS,NODE_ME)=LMMAXC
             GDES%NPRO_COL =GDES%NPRO_COL +LMMAXC

             GDES%NLM_POS (GDES%NIONS,NODE_ME)=NLM_ROW
             GDES%NLM_ENTRIES(GDES%NIONS,NODE_ME)=AUG_DES%LMMAX(NT)
             GDES%NLM_COL= GDES%NLM_COL+AUG_DES%LMMAX(NT)
             IF (NT /= LASTTYP) THEN
                GDES%NTYP=GDES%NTYP+1
                LASTTYP=NT
             ENDIF
             GDES%LMMAX(GDES%NTYP)=LMMAXC
             GDES%NLM_LMMAX(GDES%NTYP)=AUG_DES%LMMAX(NT)
             GDES%NITYP(GDES%NTYP)=GDES%NITYP(GDES%NTYP)+1
             GDES%NT_GLOBAL(GDES%NTYP)=NT
             GDES%NI_GLOBAL(GDES%NIONS)=NI
!$           GDES%NT_LOCAL(GDES%NIONS)=GDES%NTYP
          ENDIF
          NPRO= LMMAXC+NPRO
          NLM_ROW=AUG_DES%LMMAX(NT)+NLM_ROW
       ENDDO

100    NIS = NIS+T_INFO%NITYP(NT)
    ENDDO

    GDES%NPRO_COL_MAX= GDES%NPRO_COL
    CALLMPI( M_max_i(COMM, GDES%NPRO_COL_MAX ,1))

    GDES%NLM_COL_MAX= GDES%NLM_COL
    CALLMPI( M_max_i(COMM, GDES%NLM_COL_MAX ,1))

    ALLOCATE(GDES%NIONS_CORE(COMM%NCPU))
    GDES%NIONS_CORE=0
    GDES%NIONS_CORE(COMM%NODE_ME)=GDES%NIONS
    CALLMPI( M_sum_i(COMM, GDES%NIONS_CORE, SIZE(GDES%NIONS_CORE)) )

    CALLMPI( M_sum_i(COMM, GDES%NPRO_POS, SIZE(GDES%NPRO_POS)) )
    CALLMPI( M_sum_i(COMM, GDES%NPRO_ENTRIES, SIZE(GDES%NPRO_ENTRIES)) )

    CALLMPI( M_sum_i(COMM, GDES%NLM_POS, SIZE(GDES%NLM_POS)) )
    CALLMPI( M_sum_i(COMM, GDES%NLM_ENTRIES, SIZE(GDES%NLM_ENTRIES)) )

    ! check whether everything is right
    NPRO=SUM(GDES%LMMAX*GDES%NITYP)
    IF (NPRO/= GDES%NPRO_COL) THEN
       CALL vtutor%bug("internal ERROR 1 in GDES_SET_NPRO: " // str(NPRO) // " " // str(GDES%NPRO_COL), __FILE__, __LINE__)
    ENDIF
    NLM_ROW=SUM(GDES%NLM_LMMAX*GDES%NITYP)
    IF (NLM_ROW/= GDES%NLM_COL) THEN
       CALL vtutor%bug("internal ERROR 2 in GDES_SET_NPRO: " // str(NLM_ROW) // " " // str(GDES%NLM_COL), __FILE__, __LINE__)
    ENDIF
    
    CALLMPI( M_sum_i(COMM, NPRO ,1))
    IF (NPRO/= NPRO_ROW) THEN
       CALL vtutor%bug("internal ERROR 3 in GDES_SET_NPRO: " // str(NPRO) // " " // str(NPRO_ROW), __FILE__, __LINE__)
    ENDIF
    CALLMPI( M_sum_i(COMM, NLM_ROW ,1))
    IF (NLM_ROW/= GDES%NLM_ROW) THEN
       CALL vtutor%bug("internal ERROR 4 in GDES_SET_NPRO: " // str(NLM_ROW) // " " // str(NLM_ROW), __FILE__, __LINE__)
    ENDIF
    IF (NLM_ROW/= AUG_DES%NPRO) THEN
       CALL vtutor%bug("internal ERROR 5 in GDES_SET_NPRO: " // str(NLM_ROW) // " " // str(AUG_DES%NPRO), __FILE__, __LINE__)
    ENDIF


#endif
!
!  GDES elements related to LUSEINV are set up by INIT_SUPERCELL_FFT
!
    NULLIFY(GDES%LUSEINV, GDES%NRPLWV_ROW_DATA_POINTS_NK, GDES%NRPLWV_COL_DATA_POINTS_NK, GDES%NRPLWV_COL_MAX_DATA_POINTS_NK)
    ! TODO: think about spinors later
    GDES%NPRO_COL=GDES%NPRO_COL ! *GDES%NRSPINORS
    GDES%NPRO_ROW=GDES%NPRO_ROW ! *GDES%NRSPINORS

  END SUBROUTINE GDES_SET

!*********************************************************************
!
! the subroutine INIT_SUPERCELL initializes the supercell structure 
! at this point the simplest possible algorithm is chosen
! A supercell is created by multiplying the original cell with
! NKPX, Y and Z in the three directions
! then the replicated ions are created 
!
!*********************************************************************

  SUBROUTINE INIT_SUPERCELL( S, T_INFO, LATT_CUR, NKPX, NKPY, NKPZ, B, IU6) 
    USE lattice
    TYPE (supercell), POINTER :: S
    TYPE(type_info) :: T_INFO     ! ionic positions in original cell
    TYPE (latt)     :: LATT_CUR   ! lattice constants of original cell
    INTEGER         :: NKPX, NKPY, NKPZ  ! integer division along rec. lattice vectors for generating the mesh
    REAL(q)         :: B(3,3)     ! alternatively generating lattice for k-point mesh
    INTEGER         :: IU6        ! IO unit for output

  ! local
    INTEGER :: N, M, N1, N2, N3, NI, NI_NEW
    INTEGER         :: IA(3,3)    ! supercell in terms of original cell in integer units
    REAL(q)         :: TMP(3,3)   ! temporary vector
    REAL(q), ALLOCATABLE :: TRANS(:,:)

    IF (ASSOCIATED(S)) THEN
       CALL vtutor%bug("internal error in INIT_SUPERCELL: S is allready associated", __FILE__, __LINE__)
    ENDIF

    ALLOCATE(S)

    S%NKPX=NKPX
    S%NKPY=NKPY
    S%NKPZ=NKPZ

    IF (NKPX/=-1 .AND. NKPY/=-1 .AND. NKPZ/=-1) THEN
!
! (TODO: the supercell could be Minkovski reduced to yield near
!  orthogonal vectors)
      S%LATT_CUR%A(:,1)=LATT_CUR%A(:,1)*NKPX
      S%LATT_CUR%A(:,2)=LATT_CUR%A(:,2)*NKPY
      S%LATT_CUR%A(:,3)=LATT_CUR%A(:,3)*NKPZ

      S%NREP=NKPX*NKPY*NKPZ
      ALLOCATE(TRANS(3,S%NREP))

      CALL LATTIC(S%LATT_CUR)
      N=0
      DO N1=0,NKPX-1
         DO N2=0,NKPY-1
            DO N3=0,NKPZ-1
               N=N+1
               TRANS(1,N)=N1
               TRANS(2,N)=N2
               TRANS(3,N)=N3
            ENDDO
         ENDDO
      ENDDO
    ELSE
      ! set unit cell to generating lattice of k-point grid	
      S%LATT_CUR%A=B
      ! determine reciprocal lattice of generating lattice -> B
      CALL LATTIC(S%LATT_CUR)

      ! TMP now stores the supercell
      TMP=S%LATT_CUR%B

      CALL KARDIR(3, TMP, LATT_CUR%B)
      IA=NINT(TMP)
      IF (IU6>=0) THEN
        WRITE(IU6,*)
        WRITE(IU6,'(A)') ' Supercell in multiplies of origin unit cell:'
        WRITE(IU6,'(A)') ' ============================================'
        WRITE(IU6,'(2X,3I5)') IA
      ENDIF
      ! just in case recalculate A to reduce user rounding errors
      
      TMP=IA
      CALL DIRKAR(3, TMP, LATT_CUR%A)

      S%LATT_CUR%A=TMP
      CALL LATTIC(S%LATT_CUR)

      ! now determine generating basis in terms of original lattice vectors
      S%NREP=NINT(S%LATT_CUR%OMEGA/LATT_CUR%OMEGA)
      ALLOCATE(TRANS(3,S%NREP))
      N=0
      ! brute force search for integer Bravais points that are within the supercell 
      DO N1=0,S%NREP
         DO N2=0,S%NREP
            DO N3=0,S%NREP
               TMP(1,1)=N1
               TMP(2,1)=N2
               TMP(3,1)=N3
               ! bring to cartesian coordinates
               CALL DIRKAR(1, TMP, LATT_CUR%A)
               ! convert to direct coordinates in super cell units
               CALL KARDIR(1, TMP, S%LATT_CUR%B)
               ! bring to Wigner Seitz  super cell
               TMP(1,1)=MODULO(TMP(1,1),1.0_q)
               TMP(2,1)=MODULO(TMP(2,1),1.0_q)
               TMP(3,1)=MODULO(TMP(3,1),1.0_q)

               ! rounding issues, very close to one, make it zero 
               IF (ABS(TMP(1,1)-1.0_q)<1E-4_q) TMP(1,1)=0
               IF (ABS(TMP(2,1)-1.0_q)<1E-4_q) TMP(2,1)=0
               IF (ABS(TMP(3,1)-1.0_q)<1E-4_q) TMP(3,1)=0

               ! cartesian coordinates
               CALL DIRKAR(1, TMP, S%LATT_CUR%A)
               ! direct coordinates of original primitive cell
               CALL KARDIR(1, TMP, LATT_CUR%B)
               IF (ABS(TMP(1,1)-NINT(TMP(1,1)))>0.1  .OR. & 
                   ABS(TMP(2,1)-NINT(TMP(2,1)))>0.1  .OR. &
                   ABS(TMP(3,1)-NINT(TMP(3,1)))>0.1) THEN
                  CALL vtutor%bug("internal error in INIT_SUPERCELL: not integer multiples " // &
                     str(TMP(:,1)), __FILE__, __LINE__)
               ENDIF
               TMP(:,1)=NINT(TMP(:,1))
               
               ! compare to yet found positions
               DO M=1,N
                  IF ((TMP(1,1)-TRANS(1,M))==0 .AND. & 
                      (TMP(2,1)-TRANS(2,M))==0 .AND. &
                      (TMP(3,1)-TRANS(3,M))==0) EXIT
               ENDDO
               IF (M>N) THEN
                  N=N+1
                  IF (N>S%NREP) THEN
                     CALL vtutor%bug(str(TRANS) // "\n internal error in INIT_SUPERCELL: too many &
                        &replicated positions " // str(N) // " " // str(S%NREP) // " " // str(TMP(:,1)), &
                        __FILE__, __LINE__)
                  ENDIF
                  TRANS(1,N)=TMP(1,1)
                  TRANS(2,N)=TMP(2,1)
                  TRANS(3,N)=TMP(3,1)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
    ENDIF

    IF (N/=S%NREP) THEN
       CALL vtutor%bug(str(TRANS) // "\n internal error in INIT_SUPERCELL: not all replicated " &
          // "positions found " // str(N) // " " // str(S%NREP), __FILE__, __LINE__)
    ENDIF

    N=S%NREP

    S%T_INFO%NTYPD =T_INFO%NTYPD       ! dimension for types
    S%T_INFO%NTYP  =T_INFO%NTYP        ! number of types
    S%T_INFO%NTYPPD=T_INFO%NTYPD       ! dimension for types inc. empty spheres
    S%T_INFO%NTYPP =T_INFO%NTYP        ! number of types empty spheres
    S%T_INFO%NIOND =N*T_INFO%NIOND     ! dimension for ions
    S%T_INFO%NIONPD=N*T_INFO%NIOND     ! dimension for ions inc. empty spheres
    S%T_INFO%NIONS =N*T_INFO%NIONS     ! actual number of ions
    S%T_INFO%NIONP =N*T_INFO%NIONS     ! actual number of ions inc. empty spheres
    S%T_INFO%LSDYN =T_INFO%LSDYN       ! selective dynamics (yes/ no)
    S%T_INFO%LDIRCO=T_INFO%LDIRCO      ! positions in direct/recproc. lattice

    N=S%T_INFO%NTYPD

    ALLOCATE(S%T_INFO%NITYP(N), S%T_INFO%TYPE(N), S%T_INFO%VCA(N))
    S%T_INFO%NITYP(1:N)=T_INFO%NITYP(1:N)*S%NREP  ! number of  ions per type
    S%T_INFO%TYPE(1:N) =T_INFO%TYPE(1:N)   ! type information
    S%T_INFO%VCA(1:N)  =T_INFO%VCA(1:N)

    N=S%T_INFO%NIOND

    ALLOCATE(S%T_INFO%POSION(3,N), S%T_INFO%LSFOR(3,N), S%T_INFO%ITYP(N))
    ALLOCATE(S%POS_IN_PRIM(N), S%POSION(3,N))
    NI_NEW=0
    DO NI=1,T_INFO%NIONS
       ! maybe we should chose a minimum image convention here
       ! but results should be entirely independent of this choice
       DO M=1,S%NREP
          NI_NEW=NI_NEW+1
          ! position in fractional coordinates of original cell
          S%POSION(1,NI_NEW)=T_INFO%POSION(1,NI)+TRANS(1,M)
          S%POSION(2,NI_NEW)=T_INFO%POSION(2,NI)+TRANS(2,M)
          S%POSION(3,NI_NEW)=T_INFO%POSION(3,NI)+TRANS(3,M)
          TMP(:,1)=S%POSION(:,NI_NEW)
          ! cartesian coordinates and then fractional coordinates of supercell
          CALL DIRKAR(1, TMP, LATT_CUR%A)
          CALL KARDIR(1, TMP, S%LATT_CUR%B)

          ! bring to Wigner Seitz  super cell
!          TMP(1,1)=MODULO(TMP(1,1),1.0_q)
!          TMP(2,1)=MODULO(TMP(2,1),1.0_q)
!          TMP(3,1)=MODULO(TMP(3,1),1.0_q)
          
          ! position in fractional coordinates of super cell
          S%T_INFO%POSION(:,NI_NEW)=TMP(:,1)
          S%T_INFO%LSFOR(:,NI_NEW) =T_INFO%LSFOR(:,NI)
          S%T_INFO%ITYP(NI_NEW)    =T_INFO%ITYP(NI)
          S%POS_IN_PRIM(NI_NEW)=NI

          ! convert from supercell coordinates to unit cell coordinates
!          CALL DIRKAR(1, TMP, S%LATT_CUR%A)
!          CALL KARDIR(1, TMP, LATT_CUR%B)
!          S%POSION(:,NI_NEW)=TMP(:,1)
       ENDDO
    ENDDO
    ! pretty sure none of the remaining elements are important
    ! nullify so that we stop when something is missing
    NULLIFY( S%T_INFO%POMASS, S%T_INFO%RWIGS, S%T_INFO%ROPT, S%T_INFO%ATOMOM, &
             S%T_INFO%DARWIN_R, S%T_INFO%DARWIN_V, S%T_INFO%ZCT, S%T_INFO%RGAUS)

    DEALLOCATE(TRANS)

  END SUBROUTINE INIT_SUPERCELL

!***********************************************************************
!
! allocate inital entries in GDES that are required if inversion
! symmetry is used
!
!***********************************************************************
    
  SUBROUTINE GDES_ALLOC_K(GDES, WDES)
    USE  wave
    USE  poscar
    USE  pseudo

    TYPE (greensfdes)  GDES
    TYPE (wavedes)     WDES           ! descriptor for orbitals
    ! nasty part, we need additional data in GDES
    ALLOCATE(GDES%LUSEINV(WDES%NKPTS))
    GDES%LUSEINV=.FALSE.
    ALLOCATE(GDES%NGVECTOR_INV(WDES%NKPTS), GDES%NRPLWV_ROW_DATA_POINTS_NK(WDES%NKPTS), & 
         GDES%NRPLWV_COL_DATA_POINTS_NK(WDES%NKPTS), GDES%NRPLWV_COL_MAX_DATA_POINTS_NK(WDES%NKPTS))
    GDES%NGVECTOR_INV=0
    GDES%NRPLWV_ROW_DATA_POINTS_NK=0
    GDES%NRPLWV_COL_DATA_POINTS_NK=0
    GDES%NRPLWV_COL_MAX_DATA_POINTS_NK=0

    ALLOCATE(GDES%MAP_TO_FULL(WDES%NGDIM, 2, WDES%NKPTS))

  END SUBROUTINE GDES_ALLOC_K


!***********************************************************************
!
! set entries in GDES that are required if inversion symmetry is used
!
! input: entries in GDES:
! NGVECTOR_INV number of plane waves after applying inversion
!                    symmetry (half grid mode)
! MAP_TO_FULL(:,1,NK) index the PW G+k has in the standard mode 
!                     without inversion symmetry
! output:
! NRPLWV_ROW_DATA_POINTS_NK
! NRPLWV_COL_DATA_POINTS_NK
! NRPLWV_COL_MAX_DATA_POINTS_NK
! MAP_TO_FULL(:,2,NK) index of the PW -G-k has has in the standard mode 
!                     without inversion symmetry
!***********************************************************************

  SUBROUTINE GDES_SET_K(GDES, WDES)
    USE  wave
    USE  poscar
    USE  pseudo

    TYPE (greensfdes)  GDES
    TYPE (wavedes)     WDES           ! descriptor for orbitals
    INTEGER, ALLOCATABLE :: IGRIDIND(:,:,:)
    INTEGER :: NK, NGX, NGY, NGZ, NG1, NG2, NG3, NG1I, NG2I, NG3I, NI, NI_ORIG
    
#ifdef gammareal
    GDES%NRPLWV_ROW_DATA_POINTS_NK=GDES%NGVECTOR_INV*2
#else
    GDES%NRPLWV_ROW_DATA_POINTS_NK=GDES%NGVECTOR_INV
#endif

#ifdef MPI
    GDES%NRPLWV_COL_DATA_POINTS_NK=(GDES%NRPLWV_ROW_DATA_POINTS_NK+GDES%COMM%NCPU-1)/GDES%COMM%NCPU
    GDES%NRPLWV_ROW_DATA_POINTS_NK=GDES%NRPLWV_COL_DATA_POINTS_NK*GDES%COMM%NCPU
    GDES%NRPLWV_COL_MAX_DATA_POINTS_NK=GDES%NRPLWV_COL_DATA_POINTS_NK
    CALLMPI( M_max_i(GDES%COMM, GDES%NRPLWV_COL_MAX_DATA_POINTS_NK ,SIZE(GDES%NRPLWV_COL_MAX_DATA_POINTS_NK )))

#else
    GDES%NRPLWV_COL_DATA_POINTS_NK=GDES%NRPLWV_ROW_DATA_POINTS_NK
    GDES%NRPLWV_COL_MAXDATA_POINTS_NK=GDES%NRPLWV_ROW_DATA_POINTS_NK
#endif
    
!   now we need to set up  GDES%MAP_TO_FULL(:,2,NK)
! allocate arrays and set some values to zero
    NGX=WDES%GRID%NGX
    NGY=WDES%GRID%NGY
    NGZ=WDES%GRID%NGZ

    ALLOCATE(IGRIDIND(NGX,NGY,NGZ))

! now index the wavefunctions at the new k-point
! take G-index from old k-point, rotate G-vec, enter new index
    DO NK=1,WDES%NKPTS
       IF (GDES%LUSEINV(NK)) THEN
          ! first we set up an 3d-array, that stores the index into
          ! the compact plane wave array at the new k-point
          IGRIDIND=0
          DO NI=1,GDES%NGVECTOR_INV(NK)
             ! original index  GDES%MAP_TO_FULL(NI,1,NK)
             NI_ORIG=GDES%MAP_TO_FULL(NI,1,NK)
             NG1=WDES%IGX(NI_ORIG,NK)
             NG2=WDES%IGY(NI_ORIG,NK)
             NG3=WDES%IGZ(NI_ORIG,NK)
             IGRIDIND(MP(NG1,NGX),MP(NG2,NGY),MP(NG3,NGZ))=NI
          ENDDO
          GDES%MAP_TO_FULL(:,2,NK)=0
          ! now loop over all G vectors and search inverted ones
          ! and store their position 
          DO NI=1,WDES%NGVECTOR(NK)
             NG1I=WDES%IGX(NI,NK)
             NG2I=WDES%IGY(NI,NK)
             NG3I=WDES%IGZ(NI,NK)
             ! apply time inversion symmetry
             ! (the actual G vector is G+k -> -G-k, and then k needs to be subtracted
             NG1=-NG1I-NINT(WDES%VKPT(1,NK)*2)
             NG2=-NG2I-NINT(WDES%VKPT(2,NK)*2)
             NG3=-NG3I-NINT(WDES%VKPT(3,NK)*2)

             IF (IGRIDIND(MP(NG1,NGX),MP(NG2,NGY),MP(NG3,NGZ))>0) THEN
                IF (GDES%MAP_TO_FULL(IGRIDIND(MP(NG1,NGX),MP(NG2,NGY),MP(NG3,NGZ)),2,NK)/=0) THEN
                   CALL vtutor%bug("internal error in GDES_SET_K: found a wave vector twice " // &
                      str(NG1) // " " // str(NG2) // " " // str(NG3), __FILE__, __LINE__)
                ENDIF
                GDES%MAP_TO_FULL(IGRIDIND(MP(NG1,NGX),MP(NG2,NGY),MP(NG3,NGZ)),2,NK)=NI
                IGRIDIND(MP(NG1,NGX),MP(NG2,NGY),MP(NG3,NGZ))=0
             ENDIF
          ENDDO
          ! we must have cleared IGRIDIND completely, otherwise some problem occured
          IF (SUM(IGRIDIND)/=0) THEN
             CALL vtutor%bug("internal error in GDES_SET_K: the inverted grid does not cover all " &
                // "original wave vectors G " // str(SUM(IGRIDIND)), __FILE__, __LINE__)
          ENDIF
          ! values in GDES%MAP_TO_FULL(:,:,NK) must be all distinct, except for the Gamma-point
          ! which might be found in both maps
       ENDIF
    ENDDO
    DEALLOCATE(IGRIDIND)

    CONTAINS
!
! small helper function, produces positive indices
!

    FUNCTION MP(IND,MAXI)
      IMPLICIT NONE
      INTEGER MP,MAXI,IND
      IF (IND<=0) THEN
         MP=MAXI+IND
      ELSE
         MP=IND
      ENDIF
    END FUNCTION MP
    

  END SUBROUTINE GDES_SET_K

!*********************************************************************
!
! INIT_SUPERCELL_FFT initializes the FFT's for the supercell
!
!*********************************************************************

  SUBROUTINE INIT_SUPERCELL_FFT( S, P, LMDIM, WDES, WGW, GDES, LATT_CUR, ENCUTGW, IU6)
    USE pseudo
    USE pseudo_struct_def
    USE FOCK
    
    TYPE (supercell), POINTER :: S
    TYPE (potcar),TARGET::  P(:)
    INTEGER        :: LMDIM      ! leading dimension of arrays like CDIJ
    TYPE (wavedes) :: WDES       ! original orbital descriptor for primitive cell
    TYPE (wavedes1):: WDES1      ! original orbital descriptor for primitive cell
    TYPE (wavedes) :: WGW        ! original orbital response descriptor for primitive cell
    TYPE (greensfdes) :: GDES    ! descriptor for distributed Green's function & polarization functions
    TYPE (latt)    :: LATT_CUR   ! lattice constants of original cell
    REAL(q)        :: ENCUTGW
    INTEGER        :: IU6        ! IO unit for output
  ! local
    REAL(q) :: G1, G2, G3, F1, F2, F3, SPACING
    INTEGER :: M, NK, IG1, IG2, IG3, NG, M_INV
    INTEGER, ALLOCATABLE ::  INDEX_INTO_SUPERCELL(:,:,:)
    INTEGER, ALLOCATABLE ::  KFOUND(:,:,:)
    LOGICAL :: LINV
    REAL(q), PARAMETER :: G2ZERO=1E-12_q

    S%WDES=WDES
    ! set number of k-points to 1
    S%WDES%NKDIM=1
    S%WDES%NKPTS=1
    S%WDES%NKPTS_FOR_GEN_LAYOUT=1

    CALL WDES_SET_NPRO(S%WDES, S%T_INFO, P, WDES%LOVERL)
 
    ! set the single k-point to Gamma
    NULLIFY(S%WDES%VKPT, S%WDES%WTKPT)
    ALLOCATE(S%WDES%VKPT(3,1), S%WDES%WTKPT(1))
    S%WDES%VKPT =0
    S%WDES%WTKPT=1

    S%GRID=WDES%GRID
    IF (S%NKPX/=-1 .AND. S%NKPY/=-1 .AND. S%NKPZ/=-1) THEN
       S%GRID%NGPTAR(1)=WDES%GRID%NGPTAR(1)*S%NKPX
       S%GRID%NGPTAR(2)=WDES%GRID%NGPTAR(2)*S%NKPY
       S%GRID%NGPTAR(3)=WDES%GRID%NGPTAR(3)*S%NKPZ
    ELSE
       ! try to get same spacing as for conventional cell
       ! unfortunately the mapping is not always exact
       ! hence the FFT grids will be necessarily incompatible
       SPACING=(WDES%GRID%NGPTAR(1)/LATT_CUR%ANORM(1)+ & 
                WDES%GRID%NGPTAR(2)/LATT_CUR%ANORM(2)+ &
                WDES%GRID%NGPTAR(3)/LATT_CUR%ANORM(3))/3.0_q
       S%GRID%NGPTAR(1)=SPACING*S%LATT_CUR%ANORM(1)+0.99
       S%GRID%NGPTAR(2)=SPACING*S%LATT_CUR%ANORM(2)+0.99
       S%GRID%NGPTAR(3)=SPACING*S%LATT_CUR%ANORM(3)+0.99
    ENDIF

    CALL FFTCHK(S%GRID%NGPTAR)
    ! reinitialize the loop counters
    CALL INILGRD(S%GRID%NGPTAR(1),S%GRID%NGPTAR(2),S%GRID%NGPTAR(3),S%GRID)

    IF (IU6>=0) THEN
       WRITE(IU6,20) S%GRID%NGPTAR
20     FORMAT(/' FFT grid for supercell:   NGX =',I4,'; NGY =',I4,'; NGZ =',I4,/ & 
               ' =======================',/)
    ENDIF

    ! set S%WDES%GRID to S%GRID and generate the data layout for FFT's
#ifdef RPAgamma
    CALL GEN_LAYOUT(S%GRID,S%WDES, S%LATT_CUR%B, S%LATT_CUR%B, IU6, .TRUE., LNOGAMMA=.FALSE.)
#else
    CALL GEN_LAYOUT(S%GRID,S%WDES, S%LATT_CUR%B, S%LATT_CUR%B, IU6, .TRUE.)
#endif
    CALL GEN_INDEX (S%GRID,S%WDES, S%LATT_CUR%B, S%LATT_CUR%B, IU6, -1, .TRUE.)

    CALLMPI( MAPSET(S%GRID))  ! generate the communication maps (patterns) for FFT
    !  init FFT (required if real to complex FFT is used)
    CALL FFTINI(S%WDES%NINDPW(1,1),S%WDES%NGVECTOR(1),S%WDES%NKPTS,S%WDES%NGDIM,S%GRID)
    ! create (and destroy) the FFT plans to generate "wisdom"
    CALL FFTGRIDPLAN(S%GRID)

    ! set up the structure (FAST_AUG_FOCK) to perform the fast augmentation 
    ! of the charge density

    IF (S%WDES%LOVERL .AND. LMAX_FOCK>=0) THEN
       ! note FAST_AUG_FOCK, LMAX_FOCKAE, and NMAX_FOCKAE are global in fock.F
       ! imported via USE
       CALL FASTAUG_SUPER( S%T_INFO, P, S%WDES, S%GRID, S%LATT_CUR,  LMDIM, & 
            FAST_AUG_FOCK, S%FAST_AUG, S%TRANS_MATRIX, LMAX_FOCKAE, NMAX_FOCKAE )
       CALL RSPHER(S%GRID, S%FAST_AUG, S%LATT_CUR)
#ifdef shmem_rproj
       IF (S%FAST_AUG%COMM_shmem%NODE_ME==1) S%FAST_AUG%RPROJ= S%FAST_AUG%RPROJ*SQRT(S%LATT_CUR%OMEGA)
       CALLMPI( M_barrier(S%FAST_AUG%COMM_shmem))
#else
       S%FAST_AUG%RPROJ= S%FAST_AUG%RPROJ*SQRT(S%LATT_CUR%OMEGA)
#endif
       ! set up AUG_DES describing data layout of one center terms
       CALL SETUP_AUG_DES(S%GRID, S%WDES, S%AUG_DES, S%FAST_AUG )
    ELSE
       NULLIFY(S%TRANS_MATRIX)
       ALLOCATE(S%AUG_DES%LMMAX(SIZE(S%WDES%LMMAX)))
       ! to make live easier set LMMAX (number of channels) to zero for each type
       ! also set NPRO, NPROD etc. to zero
       S%AUG_DES%LMMAX=0
       S%AUG_DES%NPRO =0
       S%AUG_DES%NPROD=0
       S%AUG_DES%NPRO_TOT=0
    ENDIF
!
! WGW descriptor is used for FFT of overlap density related quantities
!
    S%WGW=S%WDES
    S%WGW%ENMAX=ENCUTGW
    S%GRID_RES=S%GRID
    IF (IU6>=0) THEN 
       WRITE(IU6,*)
       WRITE(IU6,*) 'Basis sets for responsefunction in the supercell:'
       WRITE(IU6,*) '================================================='
    ENDIF
#ifdef RPAgamma
    CALL GEN_LAYOUT(S%GRID_RES, S%WGW, S%LATT_CUR%B, S%LATT_CUR%B, IU6,.TRUE.,LNOGAMMA=.FALSE., LREAL=.TRUE.)
#else
    CALL GEN_LAYOUT(S%GRID_RES, S%WGW, S%LATT_CUR%B, S%LATT_CUR%B, IU6,.TRUE.)
#endif
    CALL GEN_INDEX (S%GRID_RES, S%WGW, S%LATT_CUR%B, S%LATT_CUR%B, IU6, -1, .TRUE.)

    !  init FFT (required if real to complex FFT is used)
    CALL FFTINI(S%WGW%NINDPW(1,1), S%WGW%NGVECTOR(1), S%WGW%NKPTS, S%WGW%NGDIM, S%GRID_RES)

#ifdef RPAgamma
    CALL GDES_ALLOC_K(GDES, WDES)
    S%LUSEINV     =>GDES%LUSEINV
    S%NGVECTOR_INV=>GDES%NGVECTOR_INV
#else
    NULLIFY(S%LUSEINV, S%NGVECTOR_INV)
#endif
!=====================================================================
! set up the index array to go from original orbitals in primitive
! cell to supercell indexing
!=====================================================================
    ! allocate a "huge" array to track the reciprocal lattice vectors in the supercell
    CALL SETWDES(S%WDES, S%WDES1, 1)
    NG=S%WDES1%NGVECTOR

    ALLOCATE(INDEX_INTO_SUPERCELL( & 
         MINVAL(S%WDES1%IGX(1:NG)):MAXVAL(S%WDES1%IGX(1:NG)), &
         MINVAL(S%WDES1%IGY(1:NG)):MAXVAL(S%WDES1%IGY(1:NG)), &
         MINVAL(S%WDES1%IGZ(1:NG)):MAXVAL(S%WDES1%IGZ(1:NG))))
    ALLOCATE(KFOUND( & 
         MINVAL(S%WDES1%IGX(1:NG)):MAXVAL(S%WDES1%IGX(1:NG)), &
         MINVAL(S%WDES1%IGY(1:NG)):MAXVAL(S%WDES1%IGY(1:NG)), &
         MINVAL(S%WDES1%IGZ(1:NG)):MAXVAL(S%WDES1%IGZ(1:NG))))
    ! setup the indices in that array
    INDEX_INTO_SUPERCELL=0  ! clear all entries
    KFOUND=0
    DO M=1,S%WDES1%NGVECTOR
       ! store the index M at position G_x, G_y, G_z
       INDEX_INTO_SUPERCELL(S%WDES1%IGX(M),S%WDES1%IGY(M), S%WDES1%IGZ(M))=M
    ENDDO
    ALLOCATE(S%INDEX(WDES%NGDIM, WDES%NKPTS), S%NGVECTOR(WDES%NKPTS))
    ALLOCATE(S%INDEX_INV(WDES%NGDIM, WDES%NKPTS), S%IND_IF_INV(WDES%NGDIM, WDES%NKPTS))
    S%INDEX=0
    S%INDEX_INV=0
    S%IND_IF_INV=0

    DO NK=1,WDES%NKPTS
#ifdef RPAgamma
       ! determine whether we can reduce the number of PW coefficients for this k-point
       ! this reduction will be used for G(g,g'), G(g,r'), G(r',g)
       IF (ABS(MOD(WDES%VKPT(1,NK)*2+100,1._q))> SQRT(G2ZERO) .OR. &
           ABS(MOD(WDES%VKPT(2,NK)*2+100,1._q))> SQRT(G2ZERO) .OR. &
           ABS(MOD(WDES%VKPT(3,NK)*2+100,1._q))> SQRT(G2ZERO)) THEN
          S%LUSEINV(NK)=.FALSE.
       ELSE
          S%LUSEINV(NK)=.TRUE.
       ENDIF
       M_INV=0    ! counter for the number of PW after use of inversion symmetry
#endif
       CALL SETWDES(WDES, WDES1, NK)
       S%NGVECTOR(NK)=WDES1%NGVECTOR
       DO M=1,WDES1%NGVECTOR
          ! G vector in original cell
          G1=WDES1%IGX(M)+WDES1%VKPT(1)
          G2=WDES1%IGY(M)+WDES1%VKPT(2)
          G3=WDES1%IGZ(M)+WDES1%VKPT(3)
          ! cartesian coordinates
          F1=G1*LATT_CUR%B(1,1)+G2*LATT_CUR%B(1,2)+G3*LATT_CUR%B(1,3)
          F2=G1*LATT_CUR%B(2,1)+G2*LATT_CUR%B(2,2)+G3*LATT_CUR%B(2,3)
          F3=G1*LATT_CUR%B(3,1)+G2*LATT_CUR%B(3,2)+G3*LATT_CUR%B(3,3)
          ! supercell reciprocal coordinates
          G1=F1*S%LATT_CUR%A(1,1)+F2*S%LATT_CUR%A(2,1)+F3*S%LATT_CUR%A(3,1)
          G2=F1*S%LATT_CUR%A(1,2)+F2*S%LATT_CUR%A(2,2)+F3*S%LATT_CUR%A(3,2)
          G3=F1*S%LATT_CUR%A(1,3)+F2*S%LATT_CUR%A(2,3)+F3*S%LATT_CUR%A(3,3)

          IG1=NINT(G1)
          IG2=NINT(G2)
          IG3=NINT(G3)
          IF( ABS(G1-IG1)>1E-6 .OR. ABS(G2-IG2)>1E-6 .OR. ABS(G3-IG3)>1E-6) THEN
             CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT: can not map g+k to supercell: &
                &" // str(G1) // " " // str(IG1) // " " // str(G2) // " " // str(IG2) // " " // &
                str(G3) // " " // str(IG3), __FILE__, __LINE__)
          ENDIF
#ifdef RPAgamma
          ! skip wave vectors that are not considered in half grid mode
          LINV=.FALSE.
          IF (IG1<0) THEN
             LINV=.TRUE.
          ELSE IF (IG1==0) THEN
             IF (IG2<0) THEN
                LINV=.TRUE.
             ELSE IF (IG2==0) THEN
                IF (IG3<0) THEN
                   LINV=.TRUE.
                ENDIF
             ENDIF
          ENDIF
          IF (LINV) THEN
             IF (INDEX_INTO_SUPERCELL(-IG1, -IG2, -IG3)<=0 ) THEN
                CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT: can not find inv g+k: " // &
                   str(-IG1) // " " // str(-IG2) // " " // str(-IG3) // " " // &
                   str(INDEX_INTO_SUPERCELL(-IG1,-IG2,-IG3)), __FILE__, __LINE__)
             ENDIF
             S%INDEX_INV(M, NK)=INDEX_INTO_SUPERCELL(-IG1, -IG2, -IG3)             
             KFOUND(-IG1, -IG2, -IG3)=NK ! set k-point where it was found
             CYCLE
          ELSE
             ! Green's function G_k(g,g'), G_k(r,g'), etc. are stored applying inversion symmetry
             IF (S%LUSEINV(NK)) THEN
                M_INV=M_INV+1
                S%IND_IF_INV(M_INV, NK)=INDEX_INTO_SUPERCELL(IG1, IG2, IG3)
                GDES%MAP_TO_FULL(M_INV,1,NK)=M
             ENDIF
          ENDIF
          ! test whether this G vector has already been hit does not work for RPAgamma
          ! since some vectors are found twice (once for the conjugated plane wave)
          IF (INDEX_INTO_SUPERCELL(IG1, IG2, IG3)<=0 ) THEN
#else
          IF (INDEX_INTO_SUPERCELL(IG1, IG2, IG3)<=0 .OR. KFOUND(IG1, IG2, IG3)>0 ) THEN
#endif
             CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT: can not find g+k: " // str(IG1) &
                // " " // str(IG2) // " " // str(IG3) // " " // str(INDEX_INTO_SUPERCELL(IG1,IG2,IG3)) &
                // " " // str(KFOUND(IG1,IG2,IG3)), __FILE__, __LINE__)
          ENDIF
          S%INDEX(M, NK)=INDEX_INTO_SUPERCELL(IG1, IG2, IG3)
          ! a one to one mapping must exist, so clear the one we have already found
          KFOUND(IG1, IG2, IG3)=NK ! set k-point where it was found
       ENDDO
       ! store number of required G vectors, if inversion symmetry is used
#ifdef RPAgamma
       S%NGVECTOR_INV(NK)=M_INV    ! zero if inversion symmetry does not apply at this k-points
#endif
    ENDDO

    ! check whether all G vectors in supercell have been found
    DO M=1,S%WDES1%NGVECTOR
       IF (KFOUND(S%WDES1%IGX(M),S%WDES1%IGY(M), S%WDES1%IGZ(M))==0) THEN
          CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT: can not map G to primitive cell: " &
             // str(S%WDES1%IGX(M)) // " " // str(S%WDES1%IGY(M)) // " " // str(S%WDES1%IGZ(M)), __FILE__, __LINE__)
       ENDIF
    ENDDO

#ifdef RPAgamma
    CALL GDES_SET_K(GDES, WDES)
#endif
    DEALLOCATE(INDEX_INTO_SUPERCELL, KFOUND)
!=====================================================================
! same for response function: index array to go from original orbitals in primitive
! cell to supercell indexing
!=====================================================================
    ! allocate a "huge" array to track the position in the supercell
    CALL SETWDES(S%WGW, S%WGW1, 1)
    NG=S%WGW1%NGVECTOR

    ALLOCATE(INDEX_INTO_SUPERCELL( & 
         MINVAL(S%WGW1%IGX(1:NG)):MAXVAL(S%WGW1%IGX(1:NG)), &
         MINVAL(S%WGW1%IGY(1:NG)):MAXVAL(S%WGW1%IGY(1:NG)), &
         MINVAL(S%WGW1%IGZ(1:NG)):MAXVAL(S%WGW1%IGZ(1:NG))))
    ALLOCATE(KFOUND( & 
         MINVAL(S%WGW1%IGX(1:NG)):MAXVAL(S%WGW1%IGX(1:NG)), &
         MINVAL(S%WGW1%IGY(1:NG)):MAXVAL(S%WGW1%IGY(1:NG)), &
         MINVAL(S%WGW1%IGZ(1:NG)):MAXVAL(S%WGW1%IGZ(1:NG))))
    ! setup the indices in that array
    INDEX_INTO_SUPERCELL=0  ! clear all entries
    KFOUND=0
    DO M=1,S%WGW1%NGVECTOR
       ! store the index M at position G_x, G_y, G_z
       INDEX_INTO_SUPERCELL(S%WGW1%IGX(M),S%WGW1%IGY(M), S%WGW1%IGZ(M))=M
    ENDDO

    ALLOCATE(S%INDEX_RES(WGW%NGDIM, WGW%NKPTS),S%NGVECTOR_RES(WGW%NKPTS))
    ALLOCATE(S%INDEX_RES_INV(WGW%NGDIM, WGW%NKPTS))
    S%INDEX_RES=0
    S%INDEX_RES_INV=0

    DO NK=1,WGW%NKPTS
       CALL SETWDES(WGW, WDES1, NK)
       S%NGVECTOR_RES(NK)=WDES1%NGVECTOR
       DO M=1,WDES1%NGVECTOR
          G1=WDES1%IGX(M)+WDES1%VKPT(1)
          G2=WDES1%IGY(M)+WDES1%VKPT(2)
          G3=WDES1%IGZ(M)+WDES1%VKPT(3)

          ! cartesian coordinates
          F1=G1*LATT_CUR%B(1,1)+G2*LATT_CUR%B(1,2)+G3*LATT_CUR%B(1,3)
          F2=G1*LATT_CUR%B(2,1)+G2*LATT_CUR%B(2,2)+G3*LATT_CUR%B(2,3)
          F3=G1*LATT_CUR%B(3,1)+G2*LATT_CUR%B(3,2)+G3*LATT_CUR%B(3,3)
          ! to supercell reciprocal coordinates
          G1=F1*S%LATT_CUR%A(1,1)+F2*S%LATT_CUR%A(2,1)+F3*S%LATT_CUR%A(3,1)
          G2=F1*S%LATT_CUR%A(1,2)+F2*S%LATT_CUR%A(2,2)+F3*S%LATT_CUR%A(3,2)
          G3=F1*S%LATT_CUR%A(1,3)+F2*S%LATT_CUR%A(2,3)+F3*S%LATT_CUR%A(3,3)
 
          IG1=NINT(G1)
          IG2=NINT(G2)
          IG3=NINT(G3)
          IF( ABS(G1-IG1)>1E-6 .OR. ABS(G2-IG2)>1E-6 .OR. ABS(G3-IG3)>1E-6) THEN
             CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT 2: can not map g+k to supercell: &
                &" // str(G1) // " " // str(IG1) // " " // str(G2) // " " // str(IG2) // " " // &
                str(G3) // " " // str(IG3), __FILE__, __LINE__)
          ENDIF
#ifdef RPAgamma
          ! skip wave vectors that are not considered in half grid mode
          LINV=.FALSE.
          IF (IG1<0) THEN
             LINV=.TRUE.
          ELSE IF (IG1==0) THEN
             IF (IG2<0) THEN
                LINV=.TRUE.
             ELSE IF (IG2==0) THEN
                IF (IG3<0) THEN
                   LINV=.TRUE.
                ENDIF
             ENDIF
          ENDIF
          IF (LINV) THEN
             IF ( & 
                  -IG1<LBOUND(INDEX_INTO_SUPERCELL,1).OR. -IG1>UBOUND(INDEX_INTO_SUPERCELL,1) .OR. &
                  -IG2<LBOUND(INDEX_INTO_SUPERCELL,2).OR. -IG2>UBOUND(INDEX_INTO_SUPERCELL,2) .OR. &
                  -IG3<LBOUND(INDEX_INTO_SUPERCELL,3).OR. -IG3>UBOUND(INDEX_INTO_SUPERCELL,3)) THEN
                CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT: index error " // str(IG1) // &
                   " " // str(IG2) // " " // str(IG3) // " " // str(LBOUND(INDEX_INTO_SUPERCELL)) // &
                   " " // str(UBOUND(INDEX_INTO_SUPERCELL)), __FILE__, __LINE__)
             ENDIF
             
             IF (INDEX_INTO_SUPERCELL(-IG1, -IG2, -IG3)<=0) THEN
                CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT 2: can not find inv g+k: " // &
                   str(-IG1) // " " // str(-IG2) // " " // str(-IG3) // " " // &
                   str(INDEX_INTO_SUPERCELL(-IG1,-IG2,-IG3)), __FILE__, __LINE__)
             ENDIF
             S%INDEX_RES_INV(M, NK)=INDEX_INTO_SUPERCELL(-IG1, -IG2, -IG3)             
             KFOUND(-IG1, -IG2, -IG3)=NK ! set k-point where it was found
             CYCLE
          ENDIF
          ! test whether this G vector has already been hit does not work for RPAgamma
          ! since some vectors are found twice (once for the conjugated plane wave
          IF (INDEX_INTO_SUPERCELL(IG1, IG2, IG3)<=0 ) THEN
#else
          IF (INDEX_INTO_SUPERCELL(IG1, IG2, IG3)<=0 .OR. KFOUND(IG1, IG2, IG3)>0 ) THEN
#endif
             CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT 2: can not find g+k: " // &
                str(IG1) // " " // str(IG2) // " " // str(IG3) // " " // str(INDEX_INTO_SUPERCELL(IG1,&
                IG2,IG3)) // " " // str(KFOUND(IG1,IG2,IG3)), __FILE__, __LINE__)
          ENDIF
          IF ( & 
              IG1<LBOUND(INDEX_INTO_SUPERCELL,1).OR. IG1>UBOUND(INDEX_INTO_SUPERCELL,1) .OR. &
              IG2<LBOUND(INDEX_INTO_SUPERCELL,2).OR. IG2>UBOUND(INDEX_INTO_SUPERCELL,2) .OR. &
              IG3<LBOUND(INDEX_INTO_SUPERCELL,3).OR. IG3>UBOUND(INDEX_INTO_SUPERCELL,3)) THEN
              CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT: index error " // str(IG1) // &
                 " " // str(IG2) // " " // str(IG3) // " " // str(LBOUND(INDEX_INTO_SUPERCELL)) // " " &
                // str(UBOUND(INDEX_INTO_SUPERCELL)), __FILE__, __LINE__)
          ENDIF
          S%INDEX_RES(M, NK)=INDEX_INTO_SUPERCELL(IG1, IG2, IG3)
          ! a one to one mapping must exist, so clear the one we have already found
          KFOUND(IG1, IG2, IG3)=NK ! set k-point where it was found
       ENDDO
    ENDDO
    ! check whether all G vectors in supercell have been found
    DO M=1,S%WGW1%NGVECTOR
       IF (KFOUND(S%WGW1%IGX(M),S%WGW1%IGY(M), S%WGW1%IGZ(M))==0) THEN
          CALL vtutor%bug("internal error in INIT_SUPERCELL_FFT 2: can not map G to primitive cell: " &
             // str(S%WGW1%IGX(M)) // " " // str(S%WGW1%IGY(M)) // " " // str(S%WGW1%IGZ(M)), __FILE__, __LINE__)
       ENDIF
    ENDDO

    DEALLOCATE(INDEX_INTO_SUPERCELL, KFOUND)
!=====================================================================
! finally set the weight for the projectors when they are
! collected from primitive to supercell
!=====================================================================
    ALLOCATE(S%KWEIGHT(WDES%NKPTS))
    DO NK=1,WDES%NKPTS
       S%KWEIGHT(NK)=WDES%WTKPT(NK)*S%NREP
    ENDDO
!=====================================================================
! for debugging it is conventient to get set NONL_S structure
! so that the CPROJ array can be calculated directly in the supercell
!=====================================================================
    NULLIFY(S%NONL_S)
#ifdef debug
    ALLOCATE(S%NONL_S)
    CALL NONL_ALLOC(S%NONL_S,S%T_INFO,P,S%WDES, .FALSE.)
    ! set phase once and forever
    CALL SPHER(S%GRID, S%NONL_S, P, S%WDES, S%LATT_CUR,  1)
    CALL PHASE(S%WDES,S%NONL_S,1)
#endif

  END SUBROUTINE INIT_SUPERCELL_FFT

!*********************************************************************
!
!  DEALLOCATE_SUPERCELL  deallocates all data structures
! related to supercell
!
!*********************************************************************

  SUBROUTINE DEALLOCATE_SUPERCELL( S )
    TYPE (supercell), POINTER :: S

    IF ( ASSOCIATED(S) ) THEN
       DEALLOCATE(S%T_INFO%NITYP, S%T_INFO%TYPE, S%T_INFO%VCA)
       DEALLOCATE(S%T_INFO%POSION, S%T_INFO%LSFOR, S%T_INFO%ITYP)
       DEALLOCATE(S%POS_IN_PRIM, S%POSION)

       DEALLOCATE(S%WDES%VKPT, S%WDES%WTKPT)
       DEALLOCATE(S%INDEX, S%INDEX_INV, S%NGVECTOR)
       DEALLOCATE(S%INDEX_RES, S%INDEX_RES_INV, S%NGVECTOR_RES)
       DEALLOCATE(S%KWEIGHT)
       ! S%FAST_AUG reuses some elements of FAST_AUG
       ! deallocation is therefore not save 
       ! S%AUG_DES%LMMAX is not properly deallocated either
       ! so here we have some memory leakage
       CALL DEALLOC_FASTAUG_SUPER( S%FAST_AUG, S%TRANS_MATRIX )
       CALL DEALLOCWDES(S%WDES, .TRUE.)
       CALL DEALLOC_GRD(S%GRID)
       DEALLOCATE(S)
       NULLIFY(S)
    ENDIF


  END SUBROUTINE DEALLOCATE_SUPERCELL



   SUBROUTINE GROUP_OVERLAP_ON_NODE( NTAUPAR_ON_NODE, TOVERLAP, FOVERLAP, T, F ) 
      USE main_mpi, ONLY: COMM_WORLD
      USE main_mpi, ONLY: COMM_INTRA_NODE_WORLD
      USE main_mpi, ONLY: COMM_INTER_NODE_WORLD
      INTEGER, INTENT(INOUT) :: NTAUPAR_ON_NODE
      REAL(q), INTENT(INOUT) :: TOVERLAP(1)
      REAL(q), INTENT(INOUT) :: FOVERLAP(1)
      TYPE( loop_des ) T, F
      ! local
      INTEGER :: IPTS(COMM_INTRA_NODE_WORLD%NCPU)
      INTEGER :: IPTS_(COMM_INTRA_NODE_WORLD%NCPU)
      TOVERLAP=0
      FOVERLAP=0
      NTAUPAR_ON_NODE = 1 

#ifndef MPI
      RETURN
#else
      ! figure out if a group carries more points on one node than others
      IPTS=0
      IPTS_=0
      IPTS(COMM_INTRA_NODE_WORLD%NODE_ME)=T%DISTRIBUTION(T%COMM_BETWEEN_GROUPS%NODE_ME,2)
      CALLMPI( M_sum_i( COMM_INTRA_NODE_WORLD, IPTS(1),COMM_INTRA_NODE_WORLD%NCPU ) )
      ! if so, the first and last entry in IPTS will differ
      ! signaling that there is at least one group that has more data than
      ! the others. determine the number of ranks that have more data than others
      TOVERLAP=0
      IF ( IPTS(1) /= IPTS(COMM_INTRA_NODE_WORLD%NCPU ) ) THEN
         ! determine ranks that have more data, by subtracting the data every
         ! ranks carries
         IPTS_(COMM_INTRA_NODE_WORLD%NODE_ME) = T%DISTRIBUTION(T%COMM_BETWEEN_GROUPS%NODE_ME,2)-IPTS(COMM_INTRA_NODE_WORLD%NCPU)
         ! tell the others about it
         CALLMPI( M_sum_i( COMM_INTRA_NODE_WORLD, IPTS_(1),COMM_INTRA_NODE_WORLD%NCPU ) )
         ! and determine multiplication factor of overlapping G*G
         ! contractions in time
         TOVERLAP(1) = TOVERLAP(1)+REAL(SUM(IPTS_),q)/COMM_INTRA_NODE_WORLD%NCPU
      ENDIF   
WRITE(*,'(A,4I5,F12.4)') "T",T%COMM_BETWEEN_GROUPS%NODE_ME, T%DISTRIBUTION(T%COMM_BETWEEN_GROUPS%NODE_ME,:), COMM_INTER_NODE_WORLD%NODE_ME, TOVERLAP

      ! do the same for the frequency split
      IPTS=0
      IPTS_=0
      IPTS(COMM_INTRA_NODE_WORLD%NODE_ME)=F%DISTRIBUTION(F%COMM_BETWEEN_GROUPS%NODE_ME,2)
      CALLMPI( M_sum_i( COMM_INTRA_NODE_WORLD, IPTS(1),COMM_INTRA_NODE_WORLD%NCPU ) )
      FOVERLAP(1)=0
      IF ( IPTS(1) /= IPTS(COMM_INTRA_NODE_WORLD%NCPU ) ) THEN
         IPTS_(COMM_INTRA_NODE_WORLD%NODE_ME) = F%DISTRIBUTION(F%COMM_BETWEEN_GROUPS%NODE_ME,2)-IPTS(COMM_INTRA_NODE_WORLD%NCPU)
         CALLMPI( M_sum_i( COMM_INTRA_NODE_WORLD, IPTS_(1),COMM_INTRA_NODE_WORLD%NCPU ) )
         FOVERLAP(1) = FOVERLAP(1)+REAL(SUM(IPTS_),q)/COMM_INTRA_NODE_WORLD%NCPU
      ENDIF   
WRITE(*,'(A,4I5,F12.4)') "F",F%COMM_BETWEEN_GROUPS%NODE_ME, F%DISTRIBUTION(F%COMM_BETWEEN_GROUPS%NODE_ME,:), COMM_INTER_NODE_WORLD%NODE_ME, FOVERLAP

      ! communicate to others and use the maximum on all nodes
      CALLMPI( M_max_d( COMM_WORLD, TOVERLAP(1), 1 ) )
      CALLMPI( M_max_d( COMM_WORLD, FOVERLAP(1), 1 ) )
!STOP

#endif
   END SUBROUTINE GROUP_OVERLAP_ON_NODE


!***********************************************************************
!
! sets NTAUPAR based on available memory 
! 
!***********************************************************************
   SUBROUTINE SET_NTAUPAR_NOMEGAPAR( WDES, WDES_RESPONSE, S2E, T_INFO, &
      KPOINTS, LATT_CUR, LMDIM, P, IO)
      USE tutor, ONLY: VTUTOR, ISERROR, ISALERT
      USE mpimy, ONLY: communic 
      USE main_mpi, ONLY: COMM_WORLD, COMM_INTRA_NODE_WORLD
      USE pseudo_struct_def, ONLY: potcar 
      USE poscar, ONLY: type_info 
      USE ini, ONLY: QUERRY_ALLOCATE, DUMP_ALLOCATE_TAG
      USE mymath, ONLY : GREATEST_COMMON_DIVISOR
      TYPE( wavedes )            :: WDES
      TYPE( wavedes )            :: WDES_RESPONSE
      TYPE( screened_2e_handle ) :: S2E
      TYPE (kpoints_struct)      :: KPOINTS
      TYPE (latt)                :: LATT_CUR
      TYPE (type_info)           :: T_INFO
      INTEGER  LMDIM
      TYPE (potcar)              :: P(T_INFO%NTYP)
      TYPE( in_struct )          :: IO 
      ! local           
      INTEGER               :: NCPU = 1 
      INTEGER               :: NCPU_NODE = 1 
      INTEGER               :: IPAR 
      TYPE (greensfdes)     :: GDES_TAU  
      TYPE (greensfdes)     :: GDES     
      TYPE( mem_gw_handle ) :: MEM_GG 
      REAL(q)               :: M_CHIOMEGA 
      TYPE( loop_des )   T, F              ! loop descriptor 
      INTEGER:: IPAR_NODE
      REAL(q):: TRATIO(1)
      REAL(q):: FRATIO(1)
      TYPE( supercell ), POINTER :: S => NULL()
#ifdef MPI 
      NCPU = COMM_WORLD%NCPU 
      NCPU_NODE = COMM_INTRA_NODE_WORLD%NCPU
#endif
      M_CHIOMEGA=0

      CALL INIT_MEM_GG( MEM_GG ) 
      ! first set NOMEGAPAR 
      ! auto set to gcd 
      IF ( NOMEGAPAR == 0 .OR. LMP2LT .OR. LSMP2LT) THEN
         NOMEGAPAR =GREATEST_COMMON_DIVISOR(NOMEGA,WDES%COMM_KIN%NCPU)
         IF (IO%IU6>=0) WRITE(IO%IU6,"('   NOMEGAPAR  =',I4)") NOMEGAPAR
         
      ELSE IF ( NOMEGAPAR < 0 ) THEN
         NOMEGAPAR = 1 
         ! set GDES here, 
         ! sets also communicator in GDES_TAU 
         CALL SET_THIS_GDES( GDES_TAU, T, WDES, WDES_RESPONSE, 1, T_INFO, P ) 

         ! find maximum NOMEGAPAR for NTAUPAR = 1 
         DO IPAR = MIN(NOMEGA,NCPU),  1 , -1 
            ! total number of ranks must be divisible by IPAR -> NCPU/IPAR ranks in group
            IF( MOD( NCPU, IPAR ) /= 0 ) CYCLE 
            ! number of plane waves must be divisor of number of CPUS in tau group
            IF( MOD( WDES%NRPLWV, NCPU/IPAR) /=0  ) CYCLE 

            ! set GDES here, 
            ! this also initialises a descriptor for frequency and time loops
            CALL SET_THIS_GDES( GDES, F, WDES, WDES_RESPONSE, IPAR, T_INFO, P ) 
            ! 
            ! obtain estimate for RAM to store response functions in frequency
            ! domain
            CALL STORAGE_REQ_GG( "o", MEM_GG, GDES_TAU, GDES, T, F, S2E, WDES, WDES_RESPONSE, S )

            CALL RELEASE_LOOP_DES( F )
                  
            IF ( MEM_GG%M_CHIOMEGA < MAXMEM ) THEN
               M_CHIOMEGA = MEM_GG%M_CHIOMEGA 
               NOMEGAPAR = IPAR
               EXIT 
            ENDIF
         ENDDO

         ! delete time loop descriptor 
         CALL RELEASE_LOOP_DES( T )

         IF ( IO%IU0>=0 ) WRITE(*,9999) NOMEGAPAR
         IF ( IO%IU6>=0 ) WRITE(IO%IU6,9999) NOMEGAPAR
         IF ( IO%IU6>=0 ) WRITE(IO%IU6,1001) MEM_GG%M_CHIOMEGA

      ENDIF
9999  FORMAT( ' NOMEGAPAR set to', I3,' based on MAXMEM')
1001  FORMAT( ' MEMORY estimate per rank in MB for (frequency):', F12.2)

      ! set number of plane waves in real and reciprocal space 
      CALL SET_THIS_GDES( GDES, F, WDES, WDES_RESPONSE, NOMEGAPAR, T_INFO, P ) 

      ! auto set to 1
      IF ( NTAUPAR == 0 .OR. LMP2LT .OR. LSMP2LT) THEN
         NTAUPAR=1
         IF (IO%IU6>=0) WRITE(IO%IU6,"('   NTAUPAR  =',I4)") NTAUPAR
      ! set NTAUPAR based on available memory, 
      ! start from setting that requires most memory 
      ELSE IF ( NTAUPAR < 0 ) THEN
         NTAUPAR = 1 
         IF ( IO%IU6>=0 ) WRITE(IO%IU6,1002) 
         DO IPAR = MIN(NOMEGA,NCPU),  1 , -1 
            IF( MOD( NCPU, IPAR ) /= 0 ) CYCLE 
            ! number of plane waves must be divisor of number of CPUS in tau group
            IF( MOD( WDES%NRPLWV, NCPU/IPAR ) /=0  ) CYCLE 

            CALL SET_THIS_GDES( GDES_TAU, T, WDES, WDES_RESPONSE, IPAR, T_INFO, P, S ) 

            ! determine storage for full GW job 
            CALL STORAGE_REQ_GG( "t", MEM_GG, GDES_TAU, GDES, T, F, S2E, WDES, WDES_RESPONSE, S )
!CALL GROUP_OVERLAP_ON_NODE( IPAR_NODE, TRATIO, FRATIO, T, F )
!CALL DUMP_MEM_GG( MEM_GG, IPAR, 'MEM_GG', IO) 
            IF ( IO%IU6>=0 ) WRITE(IO%IU6,1003) IPAR, MEM_GG%M_TOTAL, &
               MEM_GG%M_TOTAL/1.05_q

            CALL RELEASE_LOOP_DES( T )
            CALL DEALLOCATE_SUPERCELL(S)

            IF ( MEM_GG%M_TOTAL< MAXMEM ) THEN
               NTAUPAR = IPAR 
               EXIT 
            ENDIF
         ENDDO
      ELSE
         IPAR=NTAUPAR
         CALL SET_THIS_GDES( GDES_TAU, T, WDES, WDES_RESPONSE, IPAR, T_INFO, P, S ) 
         ! determine storage for full GW job 
         CALL STORAGE_REQ_GG( "t", MEM_GG, GDES_TAU, GDES, T, F, S2E, WDES, WDES_RESPONSE, S )
!         CALL DUMP_MEM_GG( MEM_GG, IPAR, 'MEM_GG', IO)
         IF ( IO%IU6>=0 ) WRITE(IO%IU6,1003) IPAR, MEM_GG%M_TOTAL, &
            MEM_GG%M_TOTAL/1.05_q
         CALL RELEASE_LOOP_DES( T )
         CALL DEALLOCATE_SUPERCELL(S)
      ENDIF
!CALL DUMP_ALLOCATE_TAG(IO%IU6,"ESTIMATED"); IF (IO%IU6>=0) CALL WFORCE(IO%IU6)
1002  FORMAT(/'Looking for optimal time points distribution:',/,&
             ' NTAUPAR    req. mB/rank(max)    req. mb/rank(min)')
1003  FORMAT(I8,8X,F12.1,8X,F12.1)
      CONTAINS 

      SUBROUTINE INIT_MEM_GG( MEM_GG ) 
         TYPE( mem_gw_handle )  :: MEM_GG 
         MEM_GG%M_GTAU = 0._q
         MEM_GG%M_CHITAU = 0._q
         MEM_GG%M_CHIOMEGA = 0._q
         MEM_GG%M_GTAU_G = 0._q
         MEM_GG%M_GTAU_WTAU = 0._q
         MEM_GG%M_ORB = 0._q
         MEM_GG%M_REST = 0._q
         MEM_GG%M_SCREENED2E = 0._q
         MEM_GG%M_WAVEFUN = 0._q
         MEM_GG%M_TOTAL = 0._q
      END SUBROUTINE INIT_MEM_GG
      ! sets array sizes used in GDES based on IPAR 
      SUBROUTINE DUMP_MEM_GG( MEM_GG, NTAUPAR, TAG, IO ) 
         USE main_mpi, ONLY: COMM_WORLD
         TYPE( mem_gw_handle )  :: MEM_GG 
         INTEGER                :: NTAUPAR 
         CHARACTER( LEN=* )     :: TAG
         TYPE( in_struct )      :: IO
          
!         IF ( IO%IU6 >= 0 ) THEN
            WRITE(200+COMM_WORLD%NODE_ME, "( / A, &
                ' for NTAUPAR :', I3,/&
                '---------------------------------------'/&
                'MEM_GG%M_GTAU     (mB):', F14.4/ &
                'MEM_GG%M_GTAU_G   (mB):', F14.4/ &
                'MEM_GG%M_GTAU_WTAU(mB):', F14.4/ &
                'MEM_GG%M_CHITAU   (mB):', F14.4/ &
                'MEM_GG%M_CHIOMEGA (mB):', F14.4/ &
                'MEM_GG%M_ORB      (mB):', F14.4/ &
                'MEM_GG%M_REST     (mB):', F14.4/ &
                'MEM_GG%M_WAVEFUN  (mB):', F14.4/  &
                '---------------------------------------'/&
                'MEM_GG%TOTAL      (mB):', F14.4/)" )  &
                TAG, NTAUPAR, &
                MEM_GG%M_GTAU, &
                MEM_GG%M_GTAU_G, &
                MEM_GG%M_GTAU_WTAU, &
                MEM_GG%M_CHITAU, &
                MEM_GG%M_CHIOMEGA, &
                MEM_GG%M_ORB, &
                MEM_GG%M_REST, &
                MEM_GG%M_WAVEFUN ,  &
                MEM_GG%M_TOTAL 
!         ENDIF


      END SUBROUTINE DUMP_MEM_GG

      ! sets array sizes used in GDES based on IPAR 
      SUBROUTINE SET_THIS_GDES(GDES, T, WDES,  WDES_RESPONSE, IPAR, T_INFO, P, S)
         TYPE (greensfdes)  GDES
         TYPE( loop_des )   T              ! loop descriptor 
         TYPE (wavedes)     WDES           ! descriptor for orbitals
         TYPE (wavedes)     WDES_RESPONSE  ! descriptor for response function
         TYPE (type_info)   T_INFO
         TYPE (potcar)      P(T_INFO%NTYP)
         TYPE( supercell ), POINTER,  OPTIONAL :: S 
         INTEGER, INTENT(IN) :: IPAR 
         ! local 
         REAL(q), POINTER  :: TAUP(:)

         ! split communicator according to PAR
         ALLOCATE( TAUP(NOMEGA) ) 
         NULLIFY( T%POINTS_LOCAL ) 
         ! pass same communicator as in chi_super or chi_GG
         CALL DISTRIBUTE_IMAG_GRID( WDES%COMM_KIN, T, TAUP, IPAR, 'time', -1)
         DEALLOCATE( TAUP )
         NULLIFY( TAUP)

         ! set descriptor for green's functions and response functions
         CALL GDES_SET(GDES, WDES, WDES_RESPONSE, T_INFO, P, T%COMM_IN_GROUP)

         IF ( PRESENT( S ) .AND. ICHIREAL > 1 ) THEN
            CALL INIT_SUPERCELL( S, T_INFO, LATT_CUR, KPOINTS%NKPX, KPOINTS%NKPY,  KPOINTS%NKPZ, KPOINTS%B , -1)
            CALL INIT_SUPERCELL_FFT( S, P, LMDIM,  WDES,  WDES_RESPONSE, GDES_TAU, LATT_CUR, WDES_RESPONSE%ENMAX, -1)

#ifdef RPAgamma
            ! nasty part, we need additional data in GDES
            ALLOCATE(GDES%LUSEINV(WDES%NKPTS))
            GDES%LUSEINV=.FALSE.
            ALLOCATE(GDES%NGVECTOR_INV(WDES%NKPTS), GDES%NRPLWV_ROW_DATA_POINTS_NK(WDES%NKPTS), & 
                 GDES%NRPLWV_COL_DATA_POINTS_NK(WDES%NKPTS), GDES%NRPLWV_COL_MAX_DATA_POINTS_NK(WDES%NKPTS))
            GDES%NGVECTOR_INV=0
            GDES%NRPLWV_ROW_DATA_POINTS_NK=0
            GDES%NRPLWV_COL_DATA_POINTS_NK=0
            GDES%NRPLWV_COL_MAX_DATA_POINTS_NK=0

            ALLOCATE(GDES%MAP_TO_FULL(WDES%NGDIM, 2, WDES%NKPTS))
            S%LUSEINV     =>GDES%LUSEINV
            S%NGVECTOR_INV=>GDES%NGVECTOR_INV
#else
            NULLIFY(S%LUSEINV, S%NGVECTOR_INV)
#endif
         ENDIF
      END SUBROUTINE  SET_THIS_GDES

      ! deletes loop descriptor 
      SUBROUTINE RELEASE_LOOP_DES( DES ) 
         USE minimax_struct, ONLY: loop_des
         TYPE( loop_des ) :: DES 
         IF ( ASSOCIATED( DES%POINTS_LOCAL ) ) DEALLOCATE( DES%POINTS_LOCAL )
         IF ( ASSOCIATED( DES%GROUP_OF_NODE ) ) DEALLOCATE( DES%GROUP_OF_NODE )
         IF ( ASSOCIATED( DES%DISTRIBUTION ) ) THEN
            DEALLOCATE( DES%DISTRIBUTION )
            CALLMPI( M_freec( DES%COMM_IN_GROUP ))
            CALLMPI( M_freec( DES%COMM_BETWEEN_GROUPS ))
         ENDIF
      END SUBROUTINE RELEASE_LOOP_DES

   END SUBROUTINE SET_NTAUPAR_NOMEGAPAR

!***********************************************************************
!
! distribute times over processor groups
! the distribution is as follows:
! NTAUPAR determines how many tau groups should be created. 
! T%NPOINTS_IN_GROUP gives the # of tau points in each group
! distribution is done in a round robin fashion
! for NTAUPAR=3 and NOMEGA=8 we have for instance:
!  _______________________________
! | 1 | 2 | 3 | 1 | 2 | 3 | 1 | 2 | 
!  -------------------------------
! 
!***********************************************************************

  SUBROUTINE DISTRIBUTE_IMAG_GRID(COMM, T, TAU, NTAUPAR, TEXT, IU6)
    USE minimax_struct, ONLY: loop_des
    USE mpimy, ONLY: M_divide_general
    TYPE (communic)  :: COMM
    TYPE (loop_des)  :: T
    REAL(q), TARGET  :: TAU(:)
    INTEGER          :: NTAUPAR
    CHARACTER(LEN=*) :: TEXT     ! info
    INTEGER          :: IU6     ! unit for IO (usually OUTCAR)
  ! local
    INTEGER :: NTAU_SPLIT, NCPU_INGROUP, NTAU_INGROUP
    INTEGER :: REST,NGROUPS
    INTEGER :: I
    TYPE (communic) COMM_TMP
#ifdef MPI
    INTEGER, ALLOCATABLE :: BREAK(:)
#endif
    ! total number of points are determined here
    T%NPOINTS = SIZE( TAU )
    ! points within group will be determined below, starting value is T%NPOINTS
    T%NPOINTS_IN_GROUP = NOMEGA
    T%COMM = COMM
#ifdef MPI
    !NTAUPAR is the given # of tau groups 
    
    !allocate memory for distribution matrix
    NULLIFY( T%DISTRIBUTION )
    ALLOCATE( T%DISTRIBUTION( NTAUPAR , 2 ) )
    T%DISTRIBUTION = 0         
    
    !--------------------------------------------------------------------------  
    !determine how many cores are in one tau group
    IF ( MOD( COMM%NCPU , NTAUPAR ) == 0 )  THEN
       ! if NTAUPAR is a divisor fo NCPU then
       ! there are NCPU / NTAUPAR cores in each of the NTAUPAR groups
       !save distribution as matrix form
       T%DISTRIBUTION( 1 : NTAUPAR , 1 ) = COMM%NCPU / NTAUPAR 
    ELSE
       !otherwise we need to distribute the cores in a blockwise fasion
       !this means that group 1 to MOD( NCPU , NTAUPAR ) gets one
       !core more than the remaining ones
       T%DISTRIBUTION( 1 : MOD( COMM%NCPU , NTAUPAR ) , 1 ) = ( COMM%NCPU / NTAUPAR ) + 1
       T%DISTRIBUTION( MOD ( COMM%NCPU , NTAUPAR ) + 1 : NTAUPAR , 1 ) = COMM%NCPU / NTAUPAR
    ENDIF
    !--------------------------------------------------------------------------  

    !--------------------------------------------------------------------------  
    !determine number of tau points in each group
    IF ( MOD( T%NPOINTS_IN_GROUP , NTAUPAR ) == 0 ) THEN
       !if NTAUPAR is a divisor of T%NPOINTS_IN_GROUP
       !the number of tau points in each group is T%NPOINTS_IN_GROUP / NTAUPAR
       T%DISTRIBUTION( 1 : NTAUPAR , 2 ) = T%NPOINTS_IN_GROUP / NTAUPAR 
    ELSE
       !otherwise we need to distribute the tau points in a blockwise fasion
       !this means that group 1 to MOD( T%NPOINTS_IN_GROUP , NTAUPAR ) gets one
       !point more than the remaining ones
       T%DISTRIBUTION( 1 : MOD( T%NPOINTS_IN_GROUP , NTAUPAR ) , 2 ) = ( T%NPOINTS_IN_GROUP / NTAUPAR ) + 1
       T%DISTRIBUTION( MOD ( T%NPOINTS_IN_GROUP , NTAUPAR ) + 1 : NTAUPAR , 2 ) = T%NPOINTS_IN_GROUP / NTAUPAR
    ENDIF          
    !--------------------------------------------------------------------------  

    !do consistency checks, we don't want to loose CPUs or tau points
    IF ( SUM( T%DISTRIBUTION( : , 1 ) )  /= COMM%NCPU ) THEN 
       CALL vtutor%error("ERROR in DISTRIBUTE_IMAG_GRID: We lost some CPUs " // &
          str(SUM(T%DISTRIBUTION(:,1))) // " " // str(COMM%NCPU) // " Exiting here,...")
    ENDIF
    IF ( SUM( T%DISTRIBUTION( : , 2 ) )  /=  T%NPOINTS_IN_GROUP ) THEN 
       CALL vtutor%error("ERROR in DISTRIBUTE_IMAG_GRID: We lost some taus " // &
          str(SUM(T%DISTRIBUTION(:,2))) // " " // str(T%NPOINTS_IN_GROUP) // " Exiting here,...")
    ENDIF

    ! we have NTAUPAR groups,  split the global communicator 
    CALLMPI( M_divide_general( COMM, NTAUPAR, T%COMM_BETWEEN_GROUPS, T%COMM_IN_GROUP, T%DISTRIBUTION(:,1) ) )

    !for 4cpus,4 points and NTAUPAR=3, the call M_initc changes T%COMM_BETWEEN_GROUPS%NCPU from 3 to 1???
    ALLOCATE(BREAK(COMM%NCPU)) 
    BREAK=0
    IF ( T%COMM_BETWEEN_GROUPS%NCPU /= NTAUPAR ) THEN
       BREAK(COMM%NODE_ME)=1
       WRITE(*,'(A,4I4)')'ERROR in DISTRIBUTE_IMAG_GRID #of groups and NTAUPAR inconsistent for node',&
       COMM%NODE_ME,NTAUPAR,T%COMM_BETWEEN_GROUPS%NCPU
    ELSE
       BREAK(COMM%NODE_ME)=0
    ENDIF
    CALLMPI(M_sum_i(COMM,BREAK,COMM%NCPU))
    IF ( SUM(BREAK(1:COMM%NCPU)) > 0 ) THEN
       CALL vtutor%error("I have to stop now. Try a run with #CPUs dividable by NTAUPAR")
    ENDIF
    DEALLOCATE(BREAK)
  

#ifdef debug
    WRITE(200+COMM%NODE_ME,'("NODE:",I3," #T_GR:",I3," T_GR-ID:",I3," #T_IGR:",I3," T_IGR-ID:",I3," SIM:",I3," comm",3I3," io_in:",I2," io_bet:",I2)')&
       COMM%NODE_ME,&
       T%COMM_BETWEEN_GROUPS%NCPU,&
       T%COMM_BETWEEN_GROUPS%NODE_ME,&
       T%COMM_IN_GROUP%NCPU,&
       T%COMM_IN_GROUP%NODE_ME,&
       T%DISTRIBUTION(T%COMM_BETWEEN_GROUPS%NODE_ME,2), &
       COMM%MPI_COMM ,&
       T%COMM_IN_GROUP%MPI_COMM ,&
       T%COMM_BETWEEN_GROUPS%MPI_COMM ,&
       T%COMM_IN_GROUP%IONODE,&
       T%COMM_BETWEEN_GROUPS%IONODE
#endif

    !set T%NPOINTS_IN_GROUP to the number of tau points in the group 
    T%NPOINTS_IN_GROUP = T%DISTRIBUTION( T%COMM_BETWEEN_GROUPS%NODE_ME , 2 ) 
   
    !each group I has DISTRIBUION(I,2) tau-points
    !so the T%POINTS_LOCAL array of each group I has size T%DISTRIBUTION(I,2)
    !allocate appropriate T%POINTS_LOCAL arrays for each group
    ALLOCATE( T%POINTS_LOCAL( T%NPOINTS_IN_GROUP ) )
    !and initialize it 
    T%POINTS_LOCAL=0

!better block-cyclic distribution
    !the first MOD( T%NPOINTS_IN_GROUP , NTAUPAR ) groups have T%NPOINTS_IN_GROUP/NTAUPAR + 1 points
    !the NTAUPAR - MOD( T%NPOINTS_IN_GROUP , NTAUPAR ) remaining groups have T%NPOINTS_IN_GROUP/NTAUPAR points
    DO I = 1 , T%NPOINTS_IN_GROUP 
       T%POINTS_LOCAL( I ) = TAU( T%COMM_BETWEEN_GROUPS%NODE_ME + ( I - 1 ) * NTAUPAR  )
    ENDDO

    !!show distribution
#ifdef verbose
    IF ( COMM%NODE_ME == 1 ) THEN
          WRITE( *, '( " Distributing ",I3," tau points into",I3," group(s)")')&
          T%NPOINTS_IN_GROUP, T%COMM_BETWEEN_GROUPS%NCPU   
       DO I = 1, NTAUPAR
          WRITE( *,'(" Group:",I3," has",I3," cores and", I3,'//&
               '" time point(s)")' )I , T%DISTRIBUTION( I , 1 ), T%DISTRIBUTION( I , 2 )
       ENDDO
    ENDIF
#endif
    !show distribution
    IF ( IU6 >= 0 ) THEN
          WRITE( IU6, '(/, " Distributing ",I3," ",A," points into",I3," group(s)")')&
          SUM(T%DISTRIBUTION(:,2)),TEXT, T%COMM_BETWEEN_GROUPS%NCPU
          WRITE( *, '(/, " Distributing ",I3," ",A," points into",I3," group(s)")')&
          SUM(T%DISTRIBUTION(:,2)),TEXT, T%COMM_BETWEEN_GROUPS%NCPU
       DO I = 1, NTAUPAR
          WRITE( IU6,'(" Group:",I3," has",I3," cores and", I3," ",A," point(s)")' )&
          I , T%DISTRIBUTION( I , 1 ), T%DISTRIBUTION( I , 2 ), TEXT
          WRITE( *,  '(" Group:",I3," has",I3," cores and", I3," ",A," point(s)")' )&
          I , T%DISTRIBUTION( I , 1 ), T%DISTRIBUTION( I , 2 ), TEXT
       ENDDO
       WRITE( IU6, *)
    ENDIF

    !final consistency check before we get out
    IF ( SIZE( T%POINTS_LOCAL ) /= T%NPOINTS_IN_GROUP ) THEN
       CALL vtutor%error("ERROR in DISTRIBUTE_IMAG_GRID: size of T%POINTS_LOCAL not matching with # &
          &of points " // str(COMM%NODE_ME) // " " // str(SIZE(T%POINTS_LOCAL)) // " " // &
          str(T%NPOINTS_IN_GROUP))
    ENDIF 

    !finally also determine which node is in which group 
    ALLOCATE( T%GROUP_OF_NODE( COMM%NCPU ) )
    T%GROUP_OF_NODE=0
    T%GROUP_OF_NODE( COMM%NODE_ME ) =  T%COMM_BETWEEN_GROUPS%NODE_ME 
   
    !finally communicate with each other 
    CALLMPI( M_sum_i( COMM, T%GROUP_OF_NODE, COMM%NCPU ) )  
#else
    T%COMM_BETWEEN_GROUPS=COMM
    T%COMM_INTAU=COMM
    T%POINTS_LOCAL => TAU
    ALLOCATE( T%GROUP_OF_NODE( COMM%NCPU ) )
    T%GROUP_OF_NODE( COMM%NODE_ME ) = T%COMM_BETWEEN_GROUPS%NODE_ME
    !allocate memory for distribution matrix
    ALLOCATE( T%DISTRIBUTION( 1 , 2 ) )
    T%DISTRIBUTION(1,1) = 1
    T%NPOINTS_IN_GROUP=SIZE(TAU)
    T%DISTRIBUTION(1,2) = T%NPOINTS_IN_GROUP  
#endif

    T%NPOINTS_IN_ROOT_GROUP = T%DISTRIBUTION( 1 , 2 )

  END SUBROUTINE DISTRIBUTE_IMAG_GRID

!***********************************************************************
! set loop descriptors
! 
! GTYPE = 'boso' ... sets the frequency loop_des IMAG_GRIDS%B
! GTYPE = 'ferm' ... sets the frequency loop_des IMAG_GRIDS%F
! GTYPE = 'time' ... sets the time grid loop_des IMAG_GRIDS%T
!
!***********************************************************************
   SUBROUTINE SET_LOOP_DES( COMM,  IMAG_GRIDS, GTYPE, IO )  
      USE minimax_struct, ONLY: imag_grid_handle
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      TYPE (communic)          :: COMM
      TYPE (imag_grid_handle)  :: IMAG_GRIDS
      CHARACTER( LEN=4 )       :: GTYPE 
      TYPE (in_struct)         :: IO
      ! local 
      INTEGER                  :: IU0, IU6
      
      IU0=IO%IU0
      IU6=IO%IU6

      IF ( GTYPE == 'boso' ) THEN
          ! exit here if imag_grids are not allocated 
          IF ( .NOT. ASSOCIATED( IMAG_GRIDS%BOS_RE ) ) THEN
             CALL vtutor%bug("SET_LOOP_DES called with 'freq' but BOS_RE not associated : " // GTYPE // " " // str(ASSOCIATED(IMAG_GRIDS%BOS_RE)), __FILE__, __LINE__)
          ENDIF
       
          NULLIFY(IMAG_GRIDS%B%POINTS_LOCAL)
          !distribute OMEGA points into NOMEGAPAR groups
          IF( LACFDT .AND. .NOT. LRPAFORCE ) THEN
             CALL DISTRIBUTE_IMAG_GRID(COMM, IMAG_GRIDS%B, IMAG_GRIDS%BOS_RE, &
             NOMEGAPAR, 'bosonic', IU6)
          ELSE
             CALL DISTRIBUTE_IMAG_GRID(COMM, IMAG_GRIDS%B, IMAG_GRIDS%BOS_RE, &
             NOMEGAPAR, 'bosonic/fermionic', IU6)
          ENDIF

      ELSE IF ( GTYPE == 'ferm' ) THEN
          ! exit here if imag_grids are not allocated 
          IF ( .NOT. ASSOCIATED( IMAG_GRIDS%FER_RE ) ) THEN
             CALL vtutor%bug("SET_LOOP_DES called with 'freq' but FER_RE not associated : " // GTYPE // " " // str(ASSOCIATED(IMAG_GRIDS%FER_RE)), __FILE__, __LINE__)
          ENDIF
       
          NULLIFY(IMAG_GRIDS%F%POINTS_LOCAL)
          !distribute TAU points into NTAUPAR groups
          CALL DISTRIBUTE_IMAG_GRID(COMM, IMAG_GRIDS%F, IMAG_GRIDS%FER_RE, &
          NOMEGAPAR, 'fermionic', -1)

      ELSE IF ( GTYPE == 'time' ) THEN
          ! exit here if imag_grids are not allocated 
          IF ( .NOT. ASSOCIATED( IMAG_GRIDS%TAU ) ) THEN
             CALL vtutor%bug("SET_LOOP_DES called with 'time' but TAU not associated : " // GTYPE // " " // str(ASSOCIATED(IMAG_GRIDS%TAU)), __FILE__, __LINE__)
          ENDIF

          NULLIFY(IMAG_GRIDS%T%POINTS_LOCAL)
          !distribute TAU points into NTAUPAR groups
          CALL DISTRIBUTE_IMAG_GRID(COMM, IMAG_GRIDS%T, IMAG_GRIDS%TAU, &
          NTAUPAR, 'time', IU6)
      ELSE 
         CALL vtutor%bug("SET_LOOP_DES reports unknown SYTPE: " // GTYPE, __FILE__, __LINE__)
      ENDIF


   END SUBROUTINE SET_LOOP_DES

!******************************************************************
!
! set up tau indices for tau loops 
! TODO: put this routine into the loop_des structure
!
!******************************************************************

  SUBROUTINE SETUP_TAU_INDICES(NTAU_ROOT, T)
     USE minimax_struct, ONLY: loop_des
     INTEGER         :: NTAU_ROOT 
     TYPE (loop_des) :: T 

     !if the local number of tau points in the group coincides with the
     !one of the root group, the tau loop variable is the same 
     T%NPOINTSC = NTAU_ROOT
     !other groups may have less tau points locally
     T%LDO_POINT_LOCAL = .TRUE.
     IF ( T%NPOINTSC > T%NPOINTS_IN_GROUP ) THEN
        ! set to last tau point in group to avoid going over array bounds
        T%NPOINTSC = T%NPOINTS_IN_GROUP
        ! but actually little work is done since most expensive calls are bypassed
        T%LDO_POINT_LOCAL = .FALSE.
     ENDIF    
     ! set current imaginary time tau point
     T%POINT_CURRENT = T%POINTS_LOCAL(T%NPOINTSC)

     !determine current smallest TAU value of all groups 
     !this is always the tau value of the first tau group 
     !this weird construction here is since M_min_d is not available
     T%POINTMIN_CURRENT = -T%POINT_CURRENT
     CALLMPI( M_max_d(T%COMM, T%POINTMIN_CURRENT, 1))
     T%POINTMIN_CURRENT = -T%POINTMIN_CURRENT
  END SUBROUTINE SETUP_TAU_INDICES


!***********************************************************************
! 
! performs a cos transformation of the response function from time to 
! frequency domain and its inverse if LFORWARD = .FALSE. 
!
!***********************************************************************
   SUBROUTINE TRANS_TIME_FREQUENCY(GDES_SEND, GDES_RECV, CHI_SEND, CHI_RECV,&
   IMAG_GRIDS, NTAU_ROOT, LFORWARD, INU)
      USE mathtools, ONLY : INVERT_REAL_MATRIX
      USE minimax_struct, ONLY : imag_grid_handle

      IMPLICIT NONE
      TYPE (greensfdes)       :: GDES_SEND            !Green's function descriptor of sender
      TYPE (greensfdes)       :: GDES_RECV            !Green's function descriptor of receiver
      TYPE (responsefunction) :: CHI_SEND             !response function of sender
      TYPE (responsefunction) :: CHI_RECV             !response function of receiver
      TYPE (imag_grid_handle) :: IMAG_GRIDS           !imaginary grid descriptor
      INTEGER                 :: NTAU_ROOT            !current local tau point of root group
      LOGICAL                 :: LFORWARD             !forward transformation 
      INTEGER                 :: INU                  !writing unit
      !local
      TYPE (communic)         :: GLOBALCOMM           !global communicator
      TYPE (communic)         :: COMM_BETWEENRECV     !communicator between receiver groups
      TYPE (communic)         :: COMM_INRECV          !communicator in each receiver group 
      TYPE (communic)         :: COMM_BETWEENSEND     !communicator between sender groups
      TYPE (communic)         :: COMM_INSEND          !communicator in each sender group 
      INTEGER                 :: NPT_TOTAL            !total number of grid points 
      INTEGER                 :: NPT_RECV_SIMULTANEOUS!number of points in group of receiver
      INTEGER                 :: NPT_SEND_SIMULTANEOUS!number of points in group of sender

      REAL(q)                 :: FACT                 !scaling factor for FTCOS
      REAL(q),ALLOCATABLE     :: FT(:,:)              !tranformation matrix
      INTEGER                 :: NROWS                !total number of rows in response function 
      INTEGER                 :: ICOL_RECV            !column index of receiver
      INTEGER                 :: ICOL_SEND            !column index of sender       
      INTEGER                 :: I,J
      INTEGER                 :: IPT_SEND             ! 
      INTEGER                 :: IPT_RECV             ! 
      INTEGER                 :: IPT_RECV_LOCAL       ! local point of receiver
      INTEGER                 :: IPT_SEND_LOCAL       ! local point of sender
      INTEGER                 :: IPT_RECV_GLOBAL      ! global point of receiver
      INTEGER                 :: IPT_SEND_GLOBAL      ! global point of receiver
      INTEGER                 :: NPT_SEND_TOTAL
      INTEGER                 :: NPT_RECV_TOTAL      !number of points in group of receiver
      GDEF,ALLOCATABLE        :: RECV_DATA(:)         !sending/receiving array 
      INTEGER                 :: NSENDPAR             !# of sender groups
      INTEGER                 :: NRECVPAR             !# of receiver groups
      INTEGER                 :: ierror
      INTEGER,ALLOCATABLE     :: SEND_COUNTS(:)
      INTEGER,ALLOCATABLE     :: SEND_DISPLS(:)
      INTEGER,ALLOCATABLE     :: RECV_COUNTS(:)
      INTEGER,ALLOCATABLE     :: RECV_DISPLS(:)

      PROFILING_START('trans_time_frequency')

      ! we define X = -G*G
      ! due to symmetry of the integrand, there's a factor of 2
      ! so in total the factor for the cos transform is 
      FACT = -IMAG_GRIDS%FACTOR
      ! global communicator
      GLOBALCOMM = IMAG_GRIDS%T%COMM
      ! total number of grid points is the same in T and B structures
      ! and should be equal to NOMEGA 
      NPT_TOTAL = IMAG_GRIDS%T%NPOINTS
   
      ! local FT matrix
      ALLOCATE(FT(NPT_TOTAL,NPT_TOTAL))

      IF (LFORWARD) THEN
         FT(:,:) = FACT*IMAG_GRIDS%TO_BOS_RE(:,:)
         ! set communicators based on FT type 
         ! from tau -> freq 
         COMM_BETWEENRECV = IMAG_GRIDS%B%COMM_BETWEEN_GROUPS
         COMM_INRECV = IMAG_GRIDS%B%COMM_IN_GROUP
         COMM_BETWEENSEND = IMAG_GRIDS%T%COMM_BETWEEN_GROUPS
         COMM_INSEND = IMAG_GRIDS%T%COMM_IN_GROUP
         NPT_RECV_SIMULTANEOUS = IMAG_GRIDS%B%NPOINTS_IN_GROUP
         NPT_SEND_SIMULTANEOUS = 1
      ELSE
         !compute inverse
         IF ( LHFCALC_GG .AND. .NOT. LFINITE_TEMPERATURE ) THEN
            FT(:,:) = FACT*IMAG_GRIDS%TO_BOS_RE(:,:)
            ! special case, one row and column is just zero
            CALL INVERT_REAL_MATRIX(FT(1:NPT_TOTAL-1,1:NPT_TOTAL-1),INU)
            FT(NPT_TOTAL,:)=0
            FT(:,NPT_TOTAL)=0
            FT(NPT_TOTAL,NPT_TOTAL)=1
         ELSE
            FT(:,:) = FACT*IMAG_GRIDS%TO_BOS_RE(:,:)
            CALL INVERT_REAL_MATRIX(FT,INU)
         ENDIF
         ! set communicators based on FT type 
         ! from tau <- freq 
         COMM_BETWEENRECV = IMAG_GRIDS%T%COMM_BETWEEN_GROUPS
         COMM_INRECV = IMAG_GRIDS%T%COMM_IN_GROUP
         COMM_BETWEENSEND = IMAG_GRIDS%B%COMM_BETWEEN_GROUPS
         COMM_INSEND = IMAG_GRIDS%B%COMM_IN_GROUP
         NPT_RECV_SIMULTANEOUS = 1
         NPT_SEND_SIMULTANEOUS = IMAG_GRIDS%B%NPOINTS_IN_GROUP
      ENDIF

      !determine global array size 
      !(should be the same in each group )
      NROWS = GDES_SEND%RES_NRPLWV_ROW_DATA_POINTS
      !also should NROWS coincide with # of rows in frequency domain
      IF ( NROWS /= GDES_RECV%RES_NRPLWV_ROW_DATA_POINTS ) THEN
         CALL vtutor%bug("TRANS_TIME_FREQUENCY: inconsistent # of rows in time and frequency " // str(NROWS) // " " // str(GDES_RECV%RES_NRPLWV_ROW_DATA_POINTS), __FILE__, __LINE__)
      ENDIF
      !at this stage all necessary checks are done. 

      IF ( LFORWARD ) THEN
         !number of sender groups
         NSENDPAR = NTAUPAR
         !number of receiver groups
         NRECVPAR = NOMEGAPAR
      ELSE
         !number of receiver groups
         NRECVPAR = NTAUPAR
         !number of sender groups
         NSENDPAR = NOMEGAPAR
      ENDIF

      !number of simultaneously sended points is
      NPT_SEND_TOTAL = NPT_SEND_SIMULTANEOUS
      !number of simultaneously received points is
      NPT_RECV_TOTAL = NPT_RECV_SIMULTANEOUS
      !find maximum number of simultaneously treated points 
      CALLMPI( M_max_i(GLOBALCOMM, NPT_SEND_TOTAL, 1))
      !find maximum number of simultaneously treated points 
      CALLMPI( M_max_i(GLOBALCOMM, NPT_RECV_TOTAL, 1))

#ifdef debug
IF (LFORWARD   ) THEN
 WRITE(400+GLOBALCOMM%NODE_ME,*)'complete matrix'
 DO ICOL_RECV=1, GDES_SEND%RES_NRPLWV_COL_DATA_POINTS
 WRITE(400+GLOBALCOMM%NODE_ME,'(20F10.4)')REAL(CHI_SEND%RESPONSEFUN(1:NROWS,ICOL_RECV,1))
 ENDDO
 WRITE(400+GLOBALCOMM%NODE_ME,*)'complete matrix orig'
 WRITE(400+GLOBALCOMM%NODE_ME,'(20F10.4)')REAL(CHI_SEND%RESPONSEFUN(:,:,1))
ENDIF
#endif


      !allocate storage for counts and displacement arrays for alltoallv     
      ALLOCATE(SEND_COUNTS(GLOBALCOMM%NCPU))
      ALLOCATE(SEND_DISPLS(GLOBALCOMM%NCPU))
      ALLOCATE(RECV_COUNTS(GLOBALCOMM%NCPU))
      ALLOCATE(RECV_DISPLS(GLOBALCOMM%NCPU))

      !in the worst case one node receives data from all nodes
      ALLOCATE( RECV_DATA(GLOBALCOMM%NCPU*NROWS))

      !------------------------------------------------------------------
      DO ICOL_RECV = 1 , GDES_RECV%RES_NRPLWV_COL_DATA_POINTS !loop over local columns of receiver
      !------------------------------------------------------------------

         !determine current columns of response function to be
         !exchanged between cores
         CALL SETUP_SEND_RECV_ARRAYS(GLOBALCOMM,COMM_INSEND,COMM_INRECV,ICOL_RECV,NROWS,&
            NSENDPAR,NRECVPAR,ICOL_SEND,SEND_COUNTS,SEND_DISPLS,RECV_COUNTS,RECV_DISPLS,LFORWARD)

         !at this stage the columns are exchanged, next step is the actual transformation
         DO IPT_SEND = 1, NPT_SEND_TOTAL !number of local grid points 
#ifdef gammareal
            !exchange blocks with all to all
            CALL MPI_alltoallv(&
               CHI_SEND%RESPONSER(1,1,MIN(IPT_SEND,NPT_SEND_SIMULTANEOUS)), &
               SEND_COUNTS, SEND_DISPLS, MPI_double_precision, RECV_DATA(1), &
               RECV_COUNTS, RECV_DISPLS, MPI_double_precision, GLOBALCOMM%MPI_COMM, ierror)
#else
            !exchange blocks with all to all
            CALL MPI_alltoallv(&
               CHI_SEND%RESPONSEFUN(1,1,MIN(IPT_SEND,NPT_SEND_SIMULTANEOUS)), &
               SEND_COUNTS, SEND_DISPLS, MPI_double_complex, RECV_DATA(1), &
               RECV_COUNTS, RECV_DISPLS, MPI_double_complex, GLOBALCOMM%MPI_COMM, ierror)
#endif
            IF ( ierror /= MPI_success ) &
               CALL vtutor%error('TRANS_TIME_FREQUENCY: Error in MPI_alltoallv' // str(ierror))

            !determine the global sending points:
            !first, each node has to know of which node it obtained data; RECV_COUNTS has this info
            counts: DO I = 1, GLOBALCOMM%NCPU
               !find nodes which sent data
               IF ( RECV_COUNTS(I) == 0 ) CYCLE counts
#ifdef debug
 WRITE(500+GLOBALCOMM%NODE_ME,'(6I6,"|",6I6)')SEND_COUNTS(:),SEND_DISPLS(:)
 WRITE(500+GLOBALCOMM%NODE_ME,'(6I6,"|",6I6)')RECV_COUNTS(:),RECV_DISPLS(:)
#endif
               ! I is one of this nodes, determine its global sending point 
               IF ( LFORWARD ) THEN
                  IPT_SEND_GLOBAL = DETERMINE_NTAU_GLOBAL(I,NTAU_ROOT,&
                  IMAG_GRIDS)
               ELSE
                  IPT_SEND_GLOBAL = DETERMINE_NOMEGA_GLOBAL(I,IPT_SEND,&
                  IMAG_GRIDS)
               ENDIF
#ifdef debug 
IF ( LFORWARD  )THEN
 WRITE(600+GLOBALCOMM%NODE_ME,'(4I3,":",20F10.4)')NTAU_ROOT,I,IPT_SEND_GLOBAL,IPT_SEND,&
     REAL(RECV_DATA(RECV_DISPLS(I)+1:RECV_DISPLS(I)+MIN(NROWS,10)))
ENDIF
#endif
               !loop over local point of receiver 
               DO IPT_RECV = 1, NPT_RECV_TOTAL
                  !determine the corresponding global point  
                  IF ( LFORWARD ) THEN
                     IPT_RECV_GLOBAL = DETERMINE_NOMEGA_GLOBAL(GLOBALCOMM%NODE_ME,IPT_RECV,&
                     IMAG_GRIDS)
                     IPT_RECV_LOCAL = MIN(IPT_RECV,NPT_RECV_SIMULTANEOUS)
                  ELSE
                     IPT_RECV_GLOBAL = &
                     DETERMINE_NTAU_GLOBAL(GLOBALCOMM%NODE_ME,NTAU_ROOT,&
                     IMAG_GRIDS)
                     IPT_RECV_LOCAL = MIN(NTAU_ROOT,NPT_RECV_SIMULTANEOUS)
                  ENDIF

                  !add to CHI_RECV
                  IF ( IPT_SEND_GLOBAL <= NPT_TOTAL .AND. IPT_RECV_GLOBAL <= NPT_TOTAL ) THEN
#ifdef debug 
IF ( LFORWARD  )THEN
 WRITE(700+GLOBALCOMM%NODE_ME,'(3I3,":",20F10.4)')IPT_RECV_LOCAL,IPT_RECV_GLOBAL,ICOL_RECV,&
     REAL(RECV_DATA(RECV_DISPLS(I)+1:RECV_DISPLS(I)+MIN(NROWS,10)))
ENDIF
#endif
#ifdef gammareal
                     CHI_RECV%RESPONSER(1:NROWS,ICOL_RECV,IPT_RECV_LOCAL)=&
                        CHI_RECV%RESPONSER(1:NROWS,ICOL_RECV,IPT_RECV_LOCAL)+&
#else
                     CHI_RECV%RESPONSEFUN(1:NROWS,ICOL_RECV,IPT_RECV_LOCAL)=&
                        CHI_RECV%RESPONSEFUN(1:NROWS,ICOL_RECV,IPT_RECV_LOCAL)+&
#endif            
                        FT(IPT_RECV_GLOBAL,IPT_SEND_GLOBAL)*&
                        RECV_DATA(RECV_DISPLS(I)+1:RECV_DISPLS(I)+NROWS)
                  ENDIF


               ENDDO
            ENDDO counts
        ENDDO
      !------------------------------------------------------------------
      ENDDO                         !loop over local columns of receiver
      !------------------------------------------------------------------
      DEALLOCATE(SEND_COUNTS)
      DEALLOCATE(SEND_DISPLS)
      DEALLOCATE(RECV_COUNTS)
      DEALLOCATE(RECV_DISPLS)
      DEALLOCATE(RECV_DATA)
      DEALLOCATE(FT)

      PROFILING_STOP('trans_time_frequency')

   END SUBROUTINE TRANS_TIME_FREQUENCY


!***********************************************************************
! 
! integrate the response function in imaginary time to obtain
! chi at omega=0
!
!***********************************************************************
   SUBROUTINE INTEGRATE_CHI_TAU(GDES_SEND, GDES_RECV, CHI_SEND, CHI_RECV,&
   IMAG_GRIDS, NTAU_ROOT, INU)
      USE mathtools, ONLY : INVERT_REAL_MATRIX
      USE minimax_struct, ONLY : imag_grid_handle

      IMPLICIT NONE
      TYPE (greensfdes)       :: GDES_SEND            !Green's function descriptor of sender
      TYPE (greensfdes)       :: GDES_RECV            !Green's function descriptor of receiver
      TYPE (responsefunction) :: CHI_SEND             !response function of sender
      TYPE (responsefunction) :: CHI_RECV             !response function of receiver
      TYPE (imag_grid_handle) :: IMAG_GRIDS           !imaginary grid descriptor
      INTEGER                 :: NTAU_ROOT            !current local tau point of root group
      INTEGER                 :: INU                  !writing unit
      !local
      TYPE (communic)         :: GLOBALCOMM           !global communicator
      TYPE (communic)         :: COMM_BETWEENRECV     !communicator between receiver groups
      TYPE (communic)         :: COMM_INRECV          !communicator in each receiver group 
      TYPE (communic)         :: COMM_BETWEENSEND     !communicator between sender groups
      TYPE (communic)         :: COMM_INSEND          !communicator in each sender group 
      INTEGER                 :: NPT_TOTAL            !total number of grid points 
      INTEGER                 :: NPT_RECV_SIMULTANEOUS!number of points in group of receiver
      INTEGER                 :: NPT_SEND_SIMULTANEOUS!number of points in group of sender

      REAL(q)                 :: FACT                 !scaling factor for FTCOS
      REAL(q),ALLOCATABLE     :: FT(:,:)              !tranformation matrix
      INTEGER                 :: NROWS                !total number of rows in response function 
      INTEGER                 :: ICOL_RECV            !column index of receiver
      INTEGER                 :: ICOL_SEND            !column index of sender       
      INTEGER                 :: I,J
      INTEGER                 :: IPT_SEND             ! 
      INTEGER                 :: IPT_RECV             ! 
      INTEGER                 :: IPT_RECV_LOCAL       ! local point of receiver
      INTEGER                 :: IPT_SEND_LOCAL       ! local point of sender
      INTEGER                 :: IPT_RECV_GLOBAL      ! global point of receiver
      INTEGER                 :: IPT_SEND_GLOBAL      ! global point of receiver
      INTEGER                 :: NPT_SEND_TOTAL
      INTEGER                 :: NPT_RECV_TOTAL      !number of points in group of receiver
      GDEF,ALLOCATABLE        :: RECV_DATA(:)         !sending/receiving array 
      INTEGER                 :: NSENDPAR             !# of sender groups
      INTEGER                 :: NRECVPAR             !# of receiver groups
      INTEGER                 :: ierror
      INTEGER,ALLOCATABLE     :: SEND_COUNTS(:)
      INTEGER,ALLOCATABLE     :: SEND_DISPLS(:)
      INTEGER,ALLOCATABLE     :: RECV_COUNTS(:)
      INTEGER,ALLOCATABLE     :: RECV_DISPLS(:)
      LOGICAL                 :: LFORWARD=.TRUE.      !forward transformation 

      PROFILING_START('integrate_chi_tau')

      ! we define X = -G*G
      ! due to symmetry of the integrand, there's a factor of 2
      ! so in total the factor for the cos transform is 
      FACT = -2*IMAG_GRIDS%FACTOR
      ! global communicator
      GLOBALCOMM = IMAG_GRIDS%T%COMM
      ! total number of grid points is the same in T and B structures
      ! and should be equal to NOMEGA 
      NPT_TOTAL = IMAG_GRIDS%T%NPOINTS
   
      ! local FT matrix
      ALLOCATE(FT(NPT_TOTAL,NPT_TOTAL))
      FT=0

      FT(1,:) = FACT*IMAG_GRIDS%TAU_WEIGHT(:)
      ! set communicators based on FT type 
      ! from tau -> freq 
      COMM_BETWEENRECV = IMAG_GRIDS%B%COMM_BETWEEN_GROUPS
      COMM_INRECV = IMAG_GRIDS%B%COMM_IN_GROUP
      COMM_BETWEENSEND = IMAG_GRIDS%T%COMM_BETWEEN_GROUPS
      COMM_INSEND = IMAG_GRIDS%T%COMM_IN_GROUP
      NPT_RECV_SIMULTANEOUS = 1
      NPT_SEND_SIMULTANEOUS = 1

      !determine global array size 
      !(should be the same in each group )
      NROWS = GDES_SEND%RES_NRPLWV_ROW_DATA_POINTS
      !also should NROWS coincide with # of rows in frequency domain
      IF ( NROWS /= GDES_RECV%RES_NRPLWV_ROW_DATA_POINTS ) THEN
         CALL vtutor%bug("TRANS_TIME_FREQUENCY: inconsistent # of rows in time and frequency " // str(NROWS) // " " // str(GDES_RECV%RES_NRPLWV_ROW_DATA_POINTS), __FILE__, __LINE__)
      ENDIF
      !at this stage all necessary checks are done. 

      !number of sender groups
      NSENDPAR = NTAUPAR
      !number of receiver groups
      NRECVPAR = NOMEGAPAR

      !number of simultaneously sended points is
      NPT_SEND_TOTAL = NPT_SEND_SIMULTANEOUS
      !number of simultaneously received points is
      NPT_RECV_TOTAL = NPT_RECV_SIMULTANEOUS
      !find maximum number of simultaneously treated points 
      CALLMPI( M_max_i(GLOBALCOMM, NPT_SEND_TOTAL, 1))
      !find maximum number of simultaneously treated points 
      CALLMPI( M_max_i(GLOBALCOMM, NPT_RECV_TOTAL, 1))

      !allocate storage for counts and displacement arrays for alltoallv     
      ALLOCATE(SEND_COUNTS(GLOBALCOMM%NCPU))
      ALLOCATE(SEND_DISPLS(GLOBALCOMM%NCPU))
      ALLOCATE(RECV_COUNTS(GLOBALCOMM%NCPU))
      ALLOCATE(RECV_DISPLS(GLOBALCOMM%NCPU))

      !in the worst case one node receives data from all nodes
      ALLOCATE( RECV_DATA(GLOBALCOMM%NCPU*NROWS))

      !------------------------------------------------------------------
      DO ICOL_RECV = 1 , GDES_RECV%RES_NRPLWV_COL_DATA_POINTS !loop over local columns of receiver
      !------------------------------------------------------------------

         !determine current columns of response function to be
         !exchanged between cores
         CALL SETUP_SEND_RECV_ARRAYS(GLOBALCOMM,COMM_INSEND,COMM_INRECV,ICOL_RECV,NROWS,&
            NSENDPAR,NRECVPAR,ICOL_SEND,SEND_COUNTS,SEND_DISPLS,RECV_COUNTS,RECV_DISPLS,LFORWARD)

         !at this stage the columns are exchanged, next step is the actual transformation
         DO IPT_SEND = 1, NPT_SEND_TOTAL !number of local grid points 
#ifdef gammareal
            !exchange blocks with all to all
            CALL MPI_alltoallv(&
               CHI_SEND%RESPONSER(1,1,MIN(IPT_SEND,NPT_SEND_SIMULTANEOUS)), &
               SEND_COUNTS, SEND_DISPLS, MPI_double_precision, RECV_DATA(1), &
               RECV_COUNTS, RECV_DISPLS, MPI_double_precision, GLOBALCOMM%MPI_COMM, ierror)
#else
            !exchange blocks with all to all
            CALL MPI_alltoallv(&
               CHI_SEND%RESPONSEFUN(1,1,MIN(IPT_SEND,NPT_SEND_SIMULTANEOUS)), &
               SEND_COUNTS, SEND_DISPLS, MPI_double_complex, RECV_DATA(1), &
               RECV_COUNTS, RECV_DISPLS, MPI_double_complex, GLOBALCOMM%MPI_COMM, ierror)
#endif
            IF ( ierror /= MPI_success ) &
               CALL vtutor%error('TRANS_TIME_FREQUENCY: Error in MPI_alltoallv' // str(ierror))

            !determine the global sending points:
            !first, each node has to know of which node it obtained data; RECV_COUNTS has this info
            counts: DO I = 1, GLOBALCOMM%NCPU
               !find nodes which sent data
               IF ( RECV_COUNTS(I) == 0 ) CYCLE counts
               ! I is one of this nodes, determine its global sending point 
               IPT_SEND_GLOBAL = DETERMINE_NTAU_GLOBAL(I,NTAU_ROOT,&
               IMAG_GRIDS)
               !loop over local point of receiver 
               DO IPT_RECV = 1, NPT_RECV_TOTAL
                  !determine the corresponding global point  
                  IPT_RECV_GLOBAL = DETERMINE_NOMEGA_GLOBAL(GLOBALCOMM%NODE_ME,IPT_RECV,&
                  IMAG_GRIDS)
                  IPT_RECV_LOCAL = MIN(IPT_RECV,NPT_RECV_SIMULTANEOUS)

                  !add to CHI_RECV
                  IF ( IPT_SEND_GLOBAL <= NPT_TOTAL .AND. IPT_RECV_GLOBAL <= NPT_TOTAL ) THEN
#ifdef gammareal
                     CHI_RECV%RESPONSER(1:NROWS,ICOL_RECV,IPT_RECV_LOCAL)=&
                        CHI_RECV%RESPONSER(1:NROWS,ICOL_RECV,IPT_RECV_LOCAL)+&
#else
                     CHI_RECV%RESPONSEFUN(1:NROWS,ICOL_RECV,IPT_RECV_LOCAL)=&
                        CHI_RECV%RESPONSEFUN(1:NROWS,ICOL_RECV,IPT_RECV_LOCAL)+&
#endif            
                        FT(IPT_RECV_GLOBAL,IPT_SEND_GLOBAL)*&
                        RECV_DATA(RECV_DISPLS(I)+1:RECV_DISPLS(I)+NROWS)
                  ENDIF


               ENDDO
            ENDDO counts
        ENDDO
      !------------------------------------------------------------------
      ENDDO                         !loop over local columns of receiver
      !------------------------------------------------------------------
      DEALLOCATE(SEND_COUNTS)
      DEALLOCATE(SEND_DISPLS)
      DEALLOCATE(RECV_COUNTS)
      DEALLOCATE(RECV_DISPLS)
      DEALLOCATE(RECV_DATA)
      DEALLOCATE(FT)

      PROFILING_STOP('integrate_chi_tau')

   END SUBROUTINE INTEGRATE_CHI_TAU

!***********************************************************************
! determines the global tau point, depending on the global ID NODE_ME 
! and the local tau point NTAU_LOCAL  
!***********************************************************************
  
   FUNCTION DETERMINE_NTAU_GLOBAL( NODE_SEND , NTAU_LOCAL, IMAG_GRIDS ) 
      USE minimax_struct, ONLY : imag_grid_handle
      IMPLICIT NONE
      INTEGER                  :: DETERMINE_NTAU_GLOBAL
      INTEGER                  :: NODE_SEND             !given sending node 
      INTEGER                  :: NTAU_LOCAL            !local tau point index
      TYPE( imag_grid_handle ) :: IMAG_GRIDS            !imaginary grids
      !local
      INTEGER                  :: NGROUP_ID             !group id of NODE_SEND
      INTEGER                  :: NG, NGI, NGI_1   
      INTEGER                  :: NTAU_LOCAL_MAX

! much better block-cyclic contribution
      !NG is the number of groups
      NG = SIZE( IMAG_GRIDS%T%DISTRIBUTION , 1 ) 
       
      !NGROUP_ID is the id of NODE_SEND
      NGROUP_ID = IMAG_GRIDS%T%GROUP_OF_NODE( NODE_SEND ) 

      !now NGROUP_ID is the group id of NODE_SEND
      !the global tau index is 
      DETERMINE_NTAU_GLOBAL = NGROUP_ID + ( NTAU_LOCAL - 1 ) * NG
      RETURN

      !if I'm here, the global index was not found  
      CALL vtutor%bug("DETERMINE_NTAU_GLOBAL: index not found " // str(NODE_SEND) // " " // str(NTAU_LOCAL), __FILE__, __LINE__)

   END FUNCTION DETERMINE_NTAU_GLOBAL

!***********************************************************************
! determines the global freq point, depending on the global ID NODE_ME 
! and the local freq point NOMEGA_LOCAL  
!***********************************************************************
  
   FUNCTION DETERMINE_NOMEGA_GLOBAL( NODE_SEND , NOMEGA_LOCAL, IMAG_GRIDS ) 
      USE minimax_struct, ONLY : imag_grid_handle
      IMPLICIT NONE
      INTEGER                  :: DETERMINE_NOMEGA_GLOBAL
      INTEGER                  :: NODE_SEND             !given sending node 
      INTEGER                  :: NOMEGA_LOCAL          !local tau point index
      TYPE( imag_grid_handle ) :: IMAG_GRIDS            !imaginary grids
      !local
      INTEGER                  :: NGROUP_ID             !group id of NODE_SEND
      INTEGER                  :: NG, NGI, NGI_1   
      INTEGER                  :: NOMEGA_LOCAL_MAX

      ! blockwise distribution
! much better block-cyclic contribution
      !NG is the number of groups
      NG = SIZE( IMAG_GRIDS%B%DISTRIBUTION , 1 ) 

      !NGROUP_ID is the id of NODE_SEND
      NGROUP_ID = IMAG_GRIDS%B%GROUP_OF_NODE( NODE_SEND )

      !now NGROUP_ID is the group id of NODE_SEND
      !the global tau index is 
      DETERMINE_NOMEGA_GLOBAL = NGROUP_ID + ( NOMEGA_LOCAL - 1 ) * NG
      RETURN

      !if I'm here, the global index was not found  
      CALL vtutor%bug("DETERMINE_NOMEGA_GLOBAL: index not found " // str(NODE_SEND) // " " // str(NOMEGA_LOCAL), __FILE__, __LINE__ )

   END FUNCTION DETERMINE_NOMEGA_GLOBAL

!***********************************************************************
! determines the frequency group, depending on 
! the global frquency point NOMEGA   
!***********************************************************************

   FUNCTION DETERMINE_CURRENT_FREQ_GROUP( NOMEGA_TOTAL )
      INTEGER DETERMINE_CURRENT_FREQ_GROUP , NOMEGA_TOTAL
      !block cyclic distribution 
      IF ( MOD( NOMEGA_TOTAL, NOMEGAPAR ) == 0 ) THEN
         DETERMINE_CURRENT_FREQ_GROUP = NOMEGAPAR
      ELSE
         DETERMINE_CURRENT_FREQ_GROUP = MOD( NOMEGA_TOTAL, NOMEGAPAR )
      ENDIF

   ENDFUNCTION

!***********************************************************************
! determines the local frequency index, depending on 
! the global frquency point NOMEGA   
!***********************************************************************

   FUNCTION DETERMINE_LOCAL_FREQ_INDEX( NOMEGA_TOTAL )
      INTEGER DETERMINE_LOCAL_FREQ_INDEX , NOMEGA_TOTAL
      !block cyclic distribution 
      IF ( MOD( NOMEGA_TOTAL, NOMEGAPAR ) == 0 ) THEN
         DETERMINE_LOCAL_FREQ_INDEX = NOMEGA_TOTAL/NOMEGAPAR
      ELSE
         DETERMINE_LOCAL_FREQ_INDEX = NOMEGA_TOTAL/NOMEGAPAR + 1
      ENDIF

   ENDFUNCTION



!***********************************************************************
! 
!>computes the head and the wings of chi in the imaginary time domain
!!(mostly copied from ADD_XI ) 
!
!***********************************************************************

  SUBROUTINE CALCULATE_LONG_WAVE_LIMIT( LATT_CUR, W, WDESQ, KPOINTS, &
       NSTRIP_TOTAL, CHI, HEADHAND, IMAG_GRIDS )
    USE sym_prec
    USE nonl_high
    USE wave_high
    USE lattice
    USE constant
    USE full_kpoints
    USE kpoints_change
    USE fock
    USE subrot_cluster
    USE minimax_struct, ONLY : imag_grid_handle
    USE crpa
    IMPLICIT NONE
    ! passed variables
    TYPE (latt)             LATT_CUR
    TYPE (wavespin)         W
    TYPE (kpoints_struct)   KPOINTS
    TYPE (wavedes1)         WDESQ
    INTEGER                 NSTRIP_TOTAL ! maximum number of pairs done simultaneously
    TYPE (responsefunction) CHI                 
    TYPE (head_handle)      HEADHAND     !long-wave limit handle 
    TYPE (imag_grid_handle) IMAG_GRIDS   !imaginary grids
  ! local variables
    INTEGER                 ISMEAR       ! smearing method 
    REAL(q)                 SIGMA        ! smearing parameter
    ! required to correct self correlation of orbitals
    TYPE (communic) :: COMM_INOMEGA      !communicator in tau
    TYPE (wavedes1) :: WDESK1 
    TYPE (wavefun1) :: W2
    TYPE (wavefun1),ALLOCATABLE :: W1(:)
    COMPLEX(q) :: WEIGHT
    INTEGER              ISP             ! spin index
    INTEGER              NK1             ! k-point
    INTEGER              NB1             ! first band index
    INTEGER              NSTRIPE         ! # of valence bands considered simultaneously
    INTEGER              N1_LAST
    INTEGER N, NB2, NGLB, NB1_INTO_TOT, NP, NB2_INTO_TOT
    INTEGER I, J, NK1_IN_KPOINTS_FULL_ORIG
    COMPLEX(q), ALLOCATABLE :: GCHG(:)   ! charge 
    COMPLEX(q), ALLOCATABLE :: GWORK(:)  ! charge weighted by 1/(w-e1-e2)
    COMPLEX(q), ALLOCATABLE :: CWORK(:)  ! work array for FFT
    COMPLEX(q), ALLOCATABLE :: CWORK_INV(:)! work array for FFT of conjugated charge
    GDEF, ALLOCATABLE :: CRHOLM(:)       ! augmentation occupancy matrix
    TYPE (wavespin) WHF
    GDEF :: CDER_BETWEEN_STATE(3)
    COMPLEX(q) :: COMEGA
    INTEGER :: NO
    LOGICAL :: USE_INVERSION_SYMMETRY
    REAL(q)  :: DCELTOT, EFERMI 

    PROFILING_START('calculate_long_wave_limit')

    ISMEAR = KPOINTS%ISMEAR
    SIGMA = KPOINTS%SIGMA
    ! required to correct self correlation of orbitals
    COMM_INOMEGA = IMAG_GRIDS%B%COMM_IN_GROUP

    ! clear head and wings to start with
    HEADHAND%CWING=0
    HEADHAND%WING=0
    HEADHAND%HEAD=0

    ! quick return if possible
    IF (.NOT. ALLOCATED(CDER_BETWEEN_STATES)) THEN
       PROFILING_STOP('calculate_long_wave_limit')
       RETURN
    ENDIF

    ! few sanity checks
    IF (WDESQ%LGAMMA .NEQV. CHI%LREAL) THEN
       CALL vtutor%bug("CALCULATE_LONG_WAVE_LIMIT: internal error WDESQ%LGAMMA .NEQV. CHI%LREAL", __FILE__, __LINE__)
    ENDIF
    IF (CHI%LREAL .AND. .NOT. CHI%LREALSTORE) THEN
       CALL vtutor%bug("CALCULATE_LONG_WAVE_LIMIT: CHI%LREAL and .NOT. CHI%LREALSTORE not supported",__FILE__, __LINE__)
    ENDIF

    !=======================================================================
    ! setup
    !=======================================================================

    ! use temporarily the fock WDES (cheaper FFT's for orbitals)
    WHF=W
    WHF%WDES => WDES_FOCK
    CALL CHECK_FULL_KPOINTS
    ! there might be better tests to see whether some symmetry is still used, however,
    ! this routine does not have access to e.g. SYMM
    IF (W%WDES%NKPTS /= KPOINTS_FULL_ORIG%NKPTS) THEN
       USE_INVERSION_SYMMETRY=.TRUE.
    ELSE
       USE_INVERSION_SYMMETRY=.FALSE.
    ENDIF

    ! the long wave contribution to chi is calculated considering NSTRIPE
    ! valence bands at the same time (on each core)
    NSTRIPE=1
    DO ISP=1,W%WDES%ISPIN
       DO NK1=1,W%WDES%NKPTS
          NSTRIPE=MAX(NSTRIPE,(LAST_FILLED_XI(W,NK1,ISP,IMAG_GRIDS%T%BETA)+W%WDES%NB_PAR-1)/W%WDES%NB_PAR)
       ENDDO
    ENDDO
    NSTRIPE=MAX(MIN(NSTRIP_STANDARD,NSTRIPE),1)

    ! allocate work space
    ALLOCATE(CRHOLM(AUG_DES%NPRO*WHF%WDES%NRSPINORS), & 
         GCHG( SIZE(CHI%RESPONSEFUN,1)), &
         GWORK(SIZE(CHI%RESPONSEFUN,1)), &
         CWORK(MAX(GRIDHF%MPLWV,WDESQ%GRID%MPLWV)), &
         CWORK_INV(MAX(GRIDHF%MPLWV,WDESQ%GRID%MPLWV)))

    CALL SETWDES(WHF%WDES, WDESK1, 0)

    CALL NEWWAV(W2, WDESK1, .TRUE.)

    ALLOCATE(W1(NSTRIPE*W%WDES%NB_PAR))
    DO N=1,NSTRIPE*W%WDES%NB_PAR
       CALL NEWWAV(W1(N) , WDESK1, .TRUE.)
    ENDDO

    ! set the phase factor for the Gamma point
    CALL PHASER_HF(GRIDHF, LATT_CUR, FAST_AUG_FOCK, (/0.0_q,0.0_q,0.0_q/))
    !=======================================================================
    DO ISP=1,W%WDES%ISPIN
    ! set fermi energy for current spin channel 
    EFERMI = W%EFERMI(ISP)
    !loop over all k-points 
    DO NK1 =1, W%WDES%NKPTS
    CALL SETWDES(WHF%WDES, WDESK1, NK1)
    ! determined the index of this k-point in the original full k-point grid
    NK1_IN_KPOINTS_FULL_ORIG=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,NK1),KPOINTS_FULL_ORIG)

    ! loop over all bands (index a) in blocks of NSTRIPE
    N1_LAST=LAST_FILLED_XI(W,NK1,ISP)/W%WDES%NB_PAR
    DO NB1=1,N1_LAST,NSTRIPE
    !=======================================================================    
       ! count the number of occupied bands accurately
       NGLB=0
       DO N=1,NSTRIPE*W%WDES%NB_PAR
          NB1_INTO_TOT=(NB1-1)*W%WDES%NB_PAR+N
          ! determine how many states have to be considered
          ! for T>0 those states are states propagating inverse in time
          IF ( LFINITE_TEMPERATURE ) THEN 
             IF ( .NOT.NEGATIVE_TIME_ORBITAL( WHF%FERTOT(NB1_INTO_TOT,NK1,ISP),       &
                   REAL(WHF%CELTOT(NB1_INTO_TOT,NK1,ISP)-EFERMI,q)*CHI%BETA )) THEN
                CYCLE 
             ELSE
                NGLB=MAX(NGLB,N)
             ENDIF
          ELSE
             IF (EMPTY_XI_ORBITAL(WHF%FERTOT(NB1_INTO_TOT,NK1,ISP))) THEN
                CYCLE
             ELSE
                NGLB=MAX(NGLB,N)
             ENDIF
          ENDIF
       ENDDO        
       ! collect NGLB occupied bands in current block 
       CALL W1_GATHER_N( WHF, NB1, NB1+NSTRIPE-1, ISP, W1, NGLB)

       NP=WDESQ%NGVECTOR
!-----------------------------------------------------------------------------
! calculate head and wings
!-----------------------------------------------------------------------------
       ! loop over unoccupied bands (distributed over all cores)
       mband: DO NB2=FIRST_POS_TIME_LOCAL_ORBITAL( WHF, NK1, ISP, EFERMI, IMAG_GRIDS%T%BETA), LAST_EMPTY_XI_LOCAL( WHF, NK1, ISP)

        NB2_INTO_TOT=(NB2-1)*W%WDES%NB_PAR+W%WDES%NB_LOW   

        ! skip dummy states introduced due to parallelization
        IF (WHF%AUX(NB2, NK1, ISP)==0 ) CYCLE
        ! skip non-positive time states
        IF ( LFINITE_TEMPERATURE ) THEN
          IF( .NOT. POSITIVE_TIME_ORBITAL( WHF%FERWE( NB2, NK1, ISP),&
              REAL(WHF%CELEN(NB2,NK1,ISP)-EFERMI,q)*IMAG_GRIDS%T%BETA )) CYCLE
        ! skip occupied states at T=0
        ELSE 
          IF( FILLED_XI_ORBITAL( WHF%FERWE( NB2, NK1, ISP) ) ) CYCLE 
        ENDIF
        
        CALL W1_COPY(ELEMENT(WHF, WDESK1, NB2, ISP), W2)
        CALL FFTWAV_W1(W2)
        ! loop over all occupied orbitals in the current block
        nband: DO N=1,NGLB
           !-----------------------------------------------------------------------------
           ! calculate the charge u_1(r) u_2*(r) 
           !-----------------------------------------------------------------------------
           
           ! determine index into the FERTOT (fermi-weight) array
           NB1_INTO_TOT=(NB1-1)*W%WDES%NB_PAR+N
           ! skip states excluded explicitly by user
           IF ( NB1_INTO_TOT<NBANDSGWLOW ) CYCLE
           ! skip non-negative time states at T>0
           IF ( LFINITE_TEMPERATURE ) THEN
              IF( .NOT. NEGATIVE_TIME_ORBITAL( WHF%FERTOT( NB1_INTO_TOT, NK1, ISP),&
                  REAL(WHF%CELTOT(NB1_INTO_TOT,NK1,ISP)-EFERMI,q)*IMAG_GRIDS%T%BETA )) CYCLE
           ! skip empty states at T=0
           ELSE
              IF ( EMPTY_XI_ORBITAL(WHF%FERTOT(NB1_INTO_TOT,NK1,ISP)) ) CYCLE 
           ENDIF
           
           CALL FOCK_CHARGE_NOINT( W1(N), W2, CWORK(1), CRHOLM(1), SIZE(CRHOLM))

           IF (USE_INVERSION_SYMMETRY .AND..NOT. CHI%LREAL) THEN
              CWORK_INV(1:WDESQ%GRID%RL%NP)=CONJG(CWORK(1:WDESQ%GRID%RL%NP))
           ENDIF

           ! extract data using wave function FFT (impose cutoff at the same time)
           CALL FFTEXT(WDESQ%NGVECTOR,WDESQ%NINDPW(1), &
                CWORK(1 ), GWORK(1 ),WDESQ%GRID,.FALSE.)
           ! divide by number of grid points and conjugate
           IF (CHI%LREAL) THEN
              ! GCHG = (\sum_r' u_1(r') u*_2(r') e-iG'r')
              ! interpreted as cosine and -sin transform
              GCHG(1:NP )=GWORK(1:NP )*(1.0_q/GRIDHF%NPLWV)
           ELSE
              ! GCHG = (\sum_r' u_1(r') u*_2(r') e-iG'r')*
              GCHG(1:NP )=CONJG(GWORK(1:NP  ))*(1.0_q/GRIDHF%NPLWV)
           ENDIF

           IF (USE_INVERSION_SYMMETRY .AND..NOT. CHI%LREAL) THEN
              ! determine overlap density at -k
              ! extract data using wave function FFT (impose cutoff at the same time)
              CALL FFTEXT(WDESQ%NGVECTOR,WDESQ%NINDPW(1), &
                   CWORK_INV(1 ), GWORK(1 ),WDESQ%GRID,.FALSE.)
              ! divide by number of grid points and conjugate
              ! GCHG = (\sum_r' u_1(r') u*_2(r') e-iG'r')*
              GWORK(1:NP )=CONJG(GWORK(1:NP  ))*(1.0_q/GRIDHF%NPLWV)
           ENDIF

           ! determine head and wing
           ! CDER_BETWEEN_STATE =  <u_k1,n1|- i d/dq_j | u_k1+q,n2>
           CALL  CDER_BETWEEN_STATES_ROTATED(CDER_BETWEEN_STATE,LATT_CUR, &
                NK1_IN_KPOINTS_FULL_ORIG, ISP, NB1_INTO_TOT, NB2_INTO_TOT)
           
           !loop over all frequencies
           DO NO = 1, HEADHAND%NOMEGA
              COMEGA=CMPLX( 0._q, HEADHAND%OMEGA(NO),q)
              WEIGHT=0
              DCELTOT = REAL(WHF%CELTOT(NB2_INTO_TOT, NK1, ISP)-WHF%CELTOT(NB1_INTO_TOT,NK1,ISP),q) ; 
                 
              IF ( LFINITE_TEMPERATURE ) THEN
                 ! laplace transformed direct mp2 needs polarizability on imaginary time grid
                 ! this is just an exponentially decaying factor of transition energy 
                 IF ( LLTDMP2 ) THEN
                    CALL SET_GG_TAU_WEIGHT( WEIGHT, &
                        REAL(WHF%CELTOT(NB2_INTO_TOT, NK1, ISP)-EFERMI,q),&
                        REAL(WHF%FERTOT(NB2_INTO_TOT, NK1, ISP),q),&
                        REAL(WHF%CELTOT(NB1_INTO_TOT, NK1, ISP)-EFERMI,q),&
                        REAL(WHF%FERTOT(NB1_INTO_TOT, NK1, ISP),q),&
                        IMAG_GRIDS%BETA, COMEGA )
                 ! evaluation on complex frequency axes
                 ELSE
                    CALL SET_GG_OMEGA_WEIGHT( WEIGHT, &
                        REAL(WHF%CELTOT(NB2_INTO_TOT, NK1, ISP)-EFERMI,q),&
                        REAL(WHF%FERTOT(NB2_INTO_TOT, NK1, ISP),q),&
                        REAL(WHF%CELTOT(NB1_INTO_TOT, NK1, ISP)-EFERMI,q),&
                        REAL(WHF%FERTOT(NB1_INTO_TOT, NK1, ISP),q),&
                        IMAG_GRIDS%BETA, COMEGA, IMAG_GRIDS%BOS_RE(1) )
                 ENDIF
                 WEIGHT=WHF%WDES%RSPIN* WHF%WDES%WTKPT(NK1)*WEIGHT
                 !opposite sign used in chi_base
                 WEIGHT=WEIGHT*WHF%AUX(NB2, NK1, ISP)
              ELSE
                 IF ( DCELTOT >= 0 ) THEN
                    IF ( LLTDMP2 ) THEN
                       CALL SET_XI_TAU_WEIGHT( WEIGHT, DCELTOT, COMEGA )
                    ! evaluation on complex frequency axes
                    ELSE
                       CALL SET_XI_WEIGHT( WEIGHT, &
                            REAL(WHF%CELTOT(NB1_INTO_TOT,NK1,ISP)-WHF%CELEN(NB2, NK1, ISP),q), &
                            0.0_q, COMEGA, ISMEAR, SIGMA)
                    ENDIF

                    WEIGHT=WHF%WDES%RSPIN* WHF%WDES%WTKPT(NK1)*&
                         (WHF%FERTOT(NB1_INTO_TOT,NK1,ISP)-WHF%FERWE(NB2, NK1, ISP))*WEIGHT
                    WEIGHT=WEIGHT*WHF%AUX(NB2, NK1, ISP)
                 ENDIF
              ENDIF
              
              !used by Friedrich et al in Phys. Rev. B 83, 121101
              !removes only diagonal part of correlated polarizability
              IF ( ALLOCATED( CRPA_WEIGHTS ) ) THEN
                 WEIGHT = WEIGHT*(1._q - CRPA_WEIGHTS( NB1_INTO_TOT,NK1,ISP )*&
                      CRPA_WEIGHTS( NB2_INTO_TOT,NK1,ISP ) )
              ENDIF
              !removes disentangled states completely
              IF ( ALLOCATED( LCRPA_STATE )) THEN
                IF ( LCRPA_STATE(ISP,NK1,NB1_INTO_TOT) .AND.  LCRPA_STATE(ISP,NK1,NB2_INTO_TOT) ) &
                   WEIGHT = 0 
              ENDIF 

              !remove all bands selected by NCRPA_BANDS
              IF ( ALLOCATED( LCRPA_BAND ) ) THEN
                IF ( LCRPA_BAND( NB1_INTO_TOT ) .AND. LCRPA_BAND(NB2_INTO_TOT) ) & 
                   WEIGHT = 0 
              ENDIF

              IF ( WEIGHT == 0 ) CYCLE                 
              ! addition or subtraction handled by factor attribute
              ! which is set externally 
              WEIGHT = WEIGHT*HEADHAND%FACTOR

              ! rotate WAVEDER 
              IF ( ALLOCATED( UCRPA ) .AND.                     & 
                 ( ( LPROJECTED  .AND. HEADHAND%FACTOR < 0 ).OR.&
                   ( LDISENTANGLED .AND. HEADHAND%FACTOR > 0 )  )  ) THEN 
                 CALL  CDER_CRPA(CDER_BETWEEN_STATE,LATT_CUR, NK1_IN_KPOINTS_FULL_ORIG,& 
                    ISP, NB1_INTO_TOT, NB2_INTO_TOT, UCRPA, NK1)
              ENDIF

              DO I=1,3
                 DO J=1,3
                    HEADHAND%HEAD(J,I,NO)= HEADHAND%HEAD(J,I,NO)+&
                         REAL(GCONJG(CDER_BETWEEN_STATE(J))*CDER_BETWEEN_STATE(I),q)*WEIGHT
                 ENDDO
              ENDDO
              ! the wing stores the limit G -> 0 (second index towards 0)
              ! accumulate  1/N (\sum_r' u_1(r') u*_2(r') e -iG'r')*  
              !               (\int -ir u_k1,n1(r) u*_k1-q,n2(r)) 
              !             1/N (\sum_r' u_1(r') u*_2(r') e -iG'r')* 
              !               [(\int  u*_k1,n1(r) d/d q u_k1+q,n2(r))*] *
              ! CDER_BETWEEN_STATE(I) holds <u*_k1,n1 | d/d(iq) u_k1+q,n2>
              ! multiply by i 
              IF (CHI%LREAL) THEN
                 ! real valued array with 2 NP entries = NP complex entries
                 DO I=1,3
                    HEADHAND%WING(1:NP,I,NO) =HEADHAND%WING(1:NP,I,NO) +&
                         GCHG(1:NP )*CDER_BETWEEN_STATE(I)*REAL(WEIGHT)
                    HEADHAND%CWING(1:NP,I,NO)=HEADHAND%CWING(1:NP,I,NO)+&
                         GCHG(1:NP )*CDER_BETWEEN_STATE(I)*REAL(WEIGHT)
                 ENDDO
              ELSE
                 IF (USE_INVERSION_SYMMETRY) THEN
                    WEIGHT=WEIGHT/2
                 ENDIF
                 ! complex array with NP entries
                 DO I=1,3
                    HEADHAND%WING(1:NP,I,NO) =HEADHAND%WING(1:NP,I,NO) +&
                         GCHG(1:NP )*(0.0_q,1.0_q)*GCONJG(CDER_BETWEEN_STATE(I))*WEIGHT
                    HEADHAND%CWING(1:NP,I,NO)=HEADHAND%CWING(1:NP,I,NO)+&
                         CONJG(GCHG(1:NP )*(0.0_q,1.0_q)*GCONJG(CDER_BETWEEN_STATE(I)))*WEIGHT
                    IF (USE_INVERSION_SYMMETRY) THEN
                    ! add contribution from -k as well (at -k the CDER is just conjugated)
                       HEADHAND%WING(1:NP,I,NO) =HEADHAND%WING(1:NP,I,NO) +&
                         GWORK(1:NP )*(0.0_q,1.0_q)*CDER_BETWEEN_STATE(I)*WEIGHT
                       HEADHAND%CWING(1:NP,I,NO)=HEADHAND%CWING(1:NP,I,NO)+&
                         CONJG(GWORK(1:NP )*(0.0_q,1.0_q)*CDER_BETWEEN_STATE(I))*WEIGHT
                    ENDIF
                 ENDDO
              ENDIF
           ENDDO   !omega loop
        ENDDO nband
       ENDDO mband !nb2
    !=======================================================================
    ENDDO   !band
    ENDDO   !kpoint
    ENDDO   !spin
    !=======================================================================

    DEALLOCATE(CRHOLM, GCHG, GWORK, CWORK, CWORK_INV)

    CALL DELWAV(W2,.TRUE.)
    DO N=1,NSTRIPE*W%WDES%NB_PAR
       CALL DELWAV(W1(N) ,.TRUE.)
    ENDDO
    DEALLOCATE(W1)

    PROFILING_STOP('calculate_long_wave_limit')

  END SUBROUTINE CALCULATE_LONG_WAVE_LIMIT

!***********************************************************************
!
! this subroutine is used to add the Drude contribution to the head of 
! the intra-band
! polarizablity in the imaginary frequency domain
!
!***********************************************************************

  SUBROUTINE ADD_DRUDE_TERM(HEADHAND,WPLASMON,VOLUME)
    USE constant    !EDEPS
    IMPLICIT NONE
    TYPE (head_handle)      :: HEADHAND        !long-wave limit handle 
    REAL(q)                 :: WPLASMON(3,3)   !plasma frequency square
    REAL(q)                 :: VOLUME
    !local
    INTEGER                 :: NO

   ! these lines are for metals not well tested so use with care
    IF (ABS(WPLASMON(1,1))+ABS(WPLASMON(2,2))+ABS(WPLASMON(3,3))>1.0 ) THEN
       DO NO=1,HEADHAND%NOMEGA
          IF (ABS(HEADHAND%OMEGA(NO))>1E-6 .AND. .NOT. LLTDMP2 ) THEN
             HEADHAND%HEAD(:,:,NO)=HEADHAND%HEAD(:,:,NO) &
              -WPLASMON(:,:)/(ABS(HEADHAND%OMEGA(NO))*ABS(HEADHAND%OMEGA(NO)))/(EDEPS/VOLUME)
          ELSE
             ! metallic singularity of head
             HEADHAND%HEAD(:,:,NO)=0
             HEADHAND%HEAD(1,1,NO)=-1000
             HEADHAND%HEAD(2,2,NO)=-1000
             HEADHAND%HEAD(3,3,NO)=-1000
          ENDIF
       ENDDO
    ENDIF
  END SUBROUTINE ADD_DRUDE_TERM


!***********************************************************************
!
! this routine subtracts the self-correlation contribution
! for partially filled states we get a contribution
!
!  f_n  (1-f_m)  |m><n|   |n><m|
!
! which can be result in very slow convergence of the frequency
! integration
!
!***********************************************************************

  SUBROUTINE SUBTRACT_SELF_CORRELATION( LATT_CUR, W, WDESQ, NSTRIP_TOTAL, CHI,&
    HEADHAND, IMAG_GRIDS, GDES )  
    USE sym_prec
    USE nonl_high
    USE wave_high
    USE lattice
    USE constant
    USE full_kpoints
    USE kpoints_change
    USE fock
    USE subrot_cluster
    USE crpa 
    USE minimax_struct, ONLY : imag_grid_handle
    IMPLICIT NONE
    ! passed variables
    TYPE (latt)             :: LATT_CUR
    TYPE (wavespin)         :: W
    TYPE (wavedes1)         :: WDESQ
    INTEGER                 :: NSTRIP_TOTAL         ! maximum number of pairs done
                                                    ! simultaneously
    TYPE (responsefunction) :: CHI                  ! response function        
    TYPE (head_handle)      :: HEADHAND             ! long-wave limit handle 
    TYPE (imag_grid_handle) :: IMAG_GRIDS           ! imaginary grids
    TYPE (greensfdes)       :: GDES                 ! green function in freq. domain  
    ! required to correct self correlation of orbitals
    INTEGER                 NTAU         ! number of tau points in total 
    ! local variables
    TYPE (wavedes1) :: WDESK1, WDESK2 
    TYPE (wavefun1),ALLOCATABLE :: W2(:) ! orbitals around Fermi-level at K2
    TYPE (wavefun1),ALLOCATABLE :: W1(:) ! orbitals around Fermi-level at K1
    COMPLEX(q) :: WEIGHT
    INTEGER              ISP             ! spin index
    INTEGER              NK1, NK1_FULL   ! 1st k-point
    INTEGER              NK2, NK2_FULL   ! 2nd k-point
    INTEGER N, N2, NGLB1, NGLB2, NB1_INTO_TOT, NP, NB2_INTO_TOT
    INTEGER N1_FIRST, N1_LAST, N2_FIRST, N2_LAST
    COMPLEX(q), ALLOCATABLE :: GCHG(:)   ! charge 
    COMPLEX(q), ALLOCATABLE :: GWORK(:)  ! charge weighted by 1/(w-e1-e2)
    COMPLEX(q), ALLOCATABLE :: CWORK(:)  ! work array for FFT
    COMPLEX(q) :: CPHASE(GRIDHF%MPLWV)
    GDEF, ALLOCATABLE :: CRHOLM(:)       ! augmentation occupancy matrix
    TYPE (wavespin) WHF
    LOGICAL :: LPHASE
    COMPLEX(q) :: COMEGA
    INTEGER :: NO, NTAUC, NOMEGA_LOCAL, NCURRENT_GROUP
    REAL(q) :: DECEL
    LOGICAL :: LCONJG1, LCONJG2

    PROFILING_START('subtract_self_correlation')
    NTAU  = IMAG_GRIDS%T%NPOINTS
    ! quick return if possible
    IF (EFERMIADD<=0) THEN
      PROFILING_STOP('subtract_self_correlation')
      RETURN
    ENDIF

    !=======================================================================
    ! setup
    !=======================================================================
    ! use temporarily another WDES
    WHF=W
    WHF%WDES => WDES_FOCK
    CALL CHECK_FULL_KPOINTS

    ! allocate work space
    ALLOCATE(CRHOLM(AUG_DES%NPRO*WHF%WDES%NRSPINORS), & 
         GCHG( SIZE(CHI%RESPONSEFUN,1)), &
         GWORK(SIZE(CHI%RESPONSEFUN,1)), &
         CWORK(MAX(GRIDHF%MPLWV,WDESQ%GRID%MPLWV)))


    !=======================================================================
    DO ISP=1,W%WDES%ISPIN
    !loop over all k-points 
    DO NK1_FULL =1, KPOINTS_FULL%NKPTS
    NK1=KPOINTS_FULL%NEQUIV(NK1_FULL)
    IF (NK1/=NK1_FULL) THEN
       LCONJG1=.TRUE.
       IF ( ABS(KPOINTS_FULL%VKPT(1,NK1_FULL)+W%WDES%VKPT(1,NK1))>1E-6 .OR. &
            ABS(KPOINTS_FULL%VKPT(2,NK1_FULL)+W%WDES%VKPT(2,NK1))>1E-6 .OR. &
            ABS(KPOINTS_FULL%VKPT(3,NK1_FULL)+W%WDES%VKPT(3,NK1))>1E-6 ) THEN
          CALL vtutor%bug("SUBTRACT_SELF_CORRELATION no inversion symmetry " // str(KPOINTS_FULL%VKPT(:,NK1_FULL)+W%WDES%VKPT(:,NK1)), __FILE__, __LINE__)
       ENDIF
    ELSE
       LCONJG1=.FALSE.
    ENDIF

    CALL SETWDES(WHF%WDES, WDESK1, NK1)

    ! determined the index of k2=k1-q
    NK2_FULL=KPOINT_IN_FULL_GRID(KPOINTS_FULL%VKPT(:,NK1_FULL)-CHI%VKPT(:),KPOINTS_FULL)
    NK2=KPOINTS_FULL%NEQUIV(NK2_FULL)
    IF (NK2/=NK2_FULL) THEN
       LCONJG2=.TRUE.
       IF ( ABS(KPOINTS_FULL%VKPT(1,NK2_FULL)+W%WDES%VKPT(1,NK2))>1E-6 .OR. &
            ABS(KPOINTS_FULL%VKPT(2,NK2_FULL)+W%WDES%VKPT(2,NK2))>1E-6 .OR. &
            ABS(KPOINTS_FULL%VKPT(3,NK2_FULL)+W%WDES%VKPT(3,NK2))>1E-6 ) THEN
          CALL vtutor%bug("SUBTRACT_SELF_CORRELATION no inversion symmetry" // str(KPOINTS_FULL%VKPT(:,NK2_FULL)+W%WDES%VKPT(:,NK2)), __FILE__, __LINE__)
       ENDIF
    ELSE
       LCONJG2=.FALSE.
    ENDIF

    CALL SETWDES(WHF%WDES,WDESK2,NK2)

    !=======================================================================
    ! count the number of partially occupied bands around Fermi-level within 
    ! EFERMIADD for both k-points
    !=======================================================================
    ! determine window of partially occupied states
    N1_FIRST=W%WDES%NB_TOT
    N1_LAST =1
    DO N=FIRST_EMPTY_XI_NOMOD(W,NK1,ISP), LAST_FILLED_XI_NOMOD(W,NK1,ISP)
       ! in energy window around EFERMIADD
       IF ( REAL(W%CELTOT(N, NK1, ISP),q) >W%EFERMI(ISP)-EFERMIADD  .AND. &
            REAL(W%CELTOT(N, NK1, ISP),q) <W%EFERMI(ISP)+EFERMIADD) THEN
          N1_FIRST=MIN(N1_FIRST,N)
          N1_LAST =MAX(N1_LAST ,N)
       ENDIF
    ENDDO

    N2_FIRST=W%WDES%NB_TOT
    N2_LAST =1
    DO N=FIRST_EMPTY_XI_NOMOD(W,NK2,ISP), LAST_FILLED_XI_NOMOD(W,NK2,ISP)
       ! in energy window around EFERMIADD
       IF ( REAL(W%CELTOT(N, NK2, ISP),q) >W%EFERMI(ISP)-EFERMIADD  .AND. &
            REAL(W%CELTOT(N, NK2, ISP),q) <W%EFERMI(ISP)+EFERMIADD) THEN
          N2_FIRST=MIN(N2_FIRST,N)
          N2_LAST =MAX(N2_LAST ,N)
       ENDIF
    ENDDO

    IF (N1_FIRST> N1_LAST .OR. N2_FIRST> N2_LAST) CYCLE

    NGLB1=(N1_LAST-N1_FIRST+1)
    NGLB2=(N2_LAST-N2_FIRST+1)

    ! set the phase factors to k1-k2
    CALL PHASER_HF(GRIDHF, LATT_CUR, FAST_AUG_FOCK, &
         KPOINTS_FULL%VKPT(:,NK1_FULL) - KPOINTS_FULL%VKPT(:,NK2_FULL))
    
    ! k1-k2-q might be any reciprocal lattice vector G
    ! CPHASE stores the required "shift" e+iGr to shift the result back to GAMMA
    CALL SETPHASE(KPOINTS_FULL%VKPT(:,NK1_FULL) - KPOINTS_FULL%VKPT(:,NK2_FULL)-CHI%VKPT(:), &
         GRIDHF,CPHASE(1),LPHASE)

    ALLOCATE(W1(NGLB1), W2(NGLB2))
    DO N=1,NGLB1
       CALL NEWWAV(W1(N) , WDESK1, .TRUE.)
    ENDDO
    DO N=1,NGLB2
       CALL NEWWAV(W2(N) , WDESK2, .TRUE.)
    ENDDO
    
    CALL W1_GATHER_GLB( WHF, N1_FIRST, N1_LAST, ISP, W1)
    ! if necessary conjugate orbital (apply time inversion symmetry)
    IF (LCONJG1) THEN
       DO N=1,NGLB1
          W1(N)%CR(1:WHF%WDES%GRID%RL%NP)=CONJG(W1(N)%CR(1:WHF%WDES%GRID%RL%NP))
          W1(N)%CPROJ(1:WHF%WDES%NPRO)=GCONJG(W1(N)%CPROJ(1:WHF%WDES%NPRO))
       ENDDO
    ENDIF

    CALL W1_GATHER_GLB( WHF, N2_FIRST, N2_LAST, ISP, W2)
    IF (LCONJG2) THEN
       DO N=1,NGLB2
          W2(N)%CR(1:WHF%WDES%GRID%RL%NP)=CONJG(W2(N)%CR(1:WHF%WDES%GRID%RL%NP))
          W2(N)%CPROJ(1:WHF%WDES%NPRO)=GCONJG(W2(N)%CPROJ(1:WHF%WDES%NPRO))
       ENDDO
    ENDIF

    NP=WDESQ%NGVECTOR
!-----------------------------------------------------------------------------
! correct response matrix for "self-correlation" error from
! calculating X(tau) = G+(tau) G-(tau)
! this version is currently only correct on a single core
!-----------------------------------------------------------------------------
       ! correct for self-correlation in the response function
       ! this version corrects for errors in any pair of states
       ! however the correction is limited to q=0, which is not really entirely correct
       ! but since the error decreases with the number of k-points, this should be of little relevance
       DO N2=1,NGLB2
          NB2_INTO_TOT=N2_FIRST-1+N2

          ! completely occupied or empty skip
          IF (FILLED_XI_ORBITAL(WHF%FERTOT( NB2_INTO_TOT, NK2, ISP))) CYCLE
          IF (EMPTY_XI_ORBITAL(WHF%FERTOT(NB2_INTO_TOT,NK2,ISP))) CYCLE
          
          DO N=1,NGLB1
             NB1_INTO_TOT=N1_FIRST-1+N

             ! completely occupied or empty skip
             IF (FILLED_XI_ORBITAL(WHF%FERTOT(NB1_INTO_TOT,NK1,ISP))) CYCLE
             IF (EMPTY_XI_ORBITAL(WHF%FERTOT(NB1_INTO_TOT,NK1,ISP))) CYCLE

             ! skip all terms larger than omegamin, since they don't need to be corrected
             DECEL = REAL(WHF%CELTOT(NB2_INTO_TOT,NK2,ISP)-WHF%CELTOT(NB1_INTO_TOT, NK1, ISP),q)
             ! terms with positive transition energy larger than OMEGAMIN are taken into account accurately 
             ! by the fourier transform. They dont have to be altered
             IF ( DECEL >= OMEGAMIN ) CYCLE 
            
!             WRITE(*,'(4I4,4F14.7)') NK1, NB1_INTO_TOT, NK2, NB2_INTO_TOT, REAL(W%CELTOT(NB1_INTO_TOT,NK1,ISP),q), REAL(W%CELTOT(NB2_INTO_TOT,NK2,ISP),q), &
!              MIN(REAL(W%CELTOT(NB1_INTO_TOT,NK1,ISP),q),W%EFERMI(ISP))-MAX(REAL(W%CELTOT(NB2_INTO_TOT,NK2,ISP),q),W%EFERMI(ISP))
               
             CALL FOCK_CHARGE_NOINT( W1(N), W2(N2), CWORK(1), CRHOLM(1), SIZE(CRHOLM))
             IF (LPHASE) CALL APPLY_PHASE( GRIDHF, CPHASE(1), CWORK(1), CWORK(1))
             CALL FFTEXT(WDESQ%NGVECTOR,WDESQ%NINDPW(1), &
                  CWORK(1 ), GWORK(1 ), WDESQ%GRID, .FALSE.)
             IF (CHI%LREAL) THEN
                ! GCHG = (\sum_r' u_1(r') u*_2(r') e-iG'r')
                ! interpreted as cosine and -sin transform
                GCHG(1:NP )=GWORK(1:NP )*(1.0_q/GRIDHF%NPLWV)
             ELSE
                ! GCHG = (\sum_r' u_1(r') u*_2(r') e-iG'r')*
                GCHG(1:NP )=CONJG(GWORK(1:NP  ))*(1.0_q/GRIDHF%NPLWV)
             ENDIF

             DO NO=1, HEADHAND%NOMEGA
                NOMEGA_LOCAL = DETERMINE_LOCAL_FREQ_INDEX( NO ) 
                NCURRENT_GROUP = DETERMINE_CURRENT_FREQ_GROUP(NO)
             
                IF ( HEADHAND%COMMBET%NODE_ME == NCURRENT_GROUP ) THEN
                   ! if state at NB1, NK1 is treated as in Gminor (occ) and NB2, NK2 as Gmajor (unocc)
                   ! we need to subtract the term from e( tau(e_1-u) -(e_2-u))
                   WEIGHT=0

                   ! remove all exponentially growing fourier transformed terms contained in chi_tau
                   ! including those not accurately represented by the transform, that is those with  DECEL < OMEGAMIN
                   DO NTAUC=1, NTAU
                      WEIGHT = WEIGHT - IMAG_GRIDS%TO_BOS_RE( NO, NTAUC )*EXP( -DECEL*IMAG_GRIDS%TAU(NTAUC) )*&
                         ((1-W%FERTOT(NB2_INTO_TOT,NK2,ISP))*W%FERTOT(NB1_INTO_TOT,NK1,ISP))
                   ENDDO

                   COMEGA=CMPLX( 0._q, HEADHAND%OMEGA(NO),q)
                   !  replace these terms with exact ones 
                   IF ( DECEL >= 0 ) THEN 
                      ! for LTDMP2 they form an exponentially decaying factor 
                      IF ( LLTDMP2 ) THEN
                         DO NTAUC=1, NTAU
                            IF ( ABS(DECEL) > XI_EMPTY_THRESHHOLD ) THEN
                               WEIGHT = WEIGHT + IMAG_GRIDS%TO_BOS_RE( NO, NTAUC )*EXP( -DECEL*IMAG_GRIDS%TAU(NTAUC) )*&
                                  ((1-W%FERTOT(NB2_INTO_TOT,NK2,ISP))*W%FERTOT(NB1_INTO_TOT,NK1,ISP))
                            ! terms smaller than XI_EMPTY_THRESHHOLD are considered to be 0 
                            ! is should be consistent with the compuation on the imaginary axis  
                            ELSE
                               WEIGHT = WEIGHT + IMAG_GRIDS%TO_BOS_RE( NO, NTAUC )*&
                                  ((1-W%FERTOT(NB2_INTO_TOT,NK2,ISP))*W%FERTOT(NB1_INTO_TOT,NK1,ISP))
                            ENDIF
                         ENDDO
                      ! in case computation is done on the imaginary frequency  axis
                      ELSE
                         ! at zero temperature the weight has not exponentially decaying factor 
                         IF ( ABS(DECEL) > XI_EMPTY_THRESHHOLD  ) &
                            WEIGHT=WEIGHT+(2*DECEL/(COMEGA**2+DECEL**2)*&
                               ((1-W%FERTOT(NB2_INTO_TOT,NK2,ISP))*W%FERTOT(NB1_INTO_TOT,NK1,ISP)))
                         ! note that there is a contribution at zero transition energy for w->0 which 
                         ! grows like beta at finite temperature. At T=0 this gives an infinite term
                         ! which we drop here 
                      ENDIF
                   ENDIF
                   ! remaining spin factor is the same 
                   WEIGHT=WEIGHT*KPOINTS_FULL%WTKPT(NK1)*W%WDES%RSPIN

                   IF (CHI%LREAL .AND. CHI%LREALSTORE) THEN
                      CALL SUBTRACT_RESPONSER_RANK1( CHI%RESPONSER(1,1,NOMEGA_LOCAL), SIZE(CHI%RESPONSER,1), SIZE(CHI%RESPONSER,2), &
                           GCHG(1), REAL(WEIGHT,q), GDES%RES_NRPLWV_ROW_DATA_POINTS, GDES%RES_NRPLWV_COL_DATA_POINTS, IMAG_GRIDS%B%COMM_IN_GROUP)
                   ELSE
                      CALL SUBTRACT_RESPONSE_RANK1( CHI%RESPONSEFUN(1,1,NOMEGA_LOCAL), SIZE(CHI%RESPONSEFUN,1), SIZE(CHI%RESPONSEFUN,2), &
                           GCHG(1), REAL(WEIGHT,q), GDES%RES_NRPLWV_ROW_DATA_POINTS, GDES%RES_NRPLWV_COL_DATA_POINTS, IMAG_GRIDS%B%COMM_IN_GROUP)
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    !=======================================================================
    DO N=1,NGLB1
       CALL DELWAV(W1(N) ,.TRUE.)
    ENDDO
    DO N=1,NGLB2
       CALL DELWAV(W2(N) ,.TRUE.)
    ENDDO

    DEALLOCATE(W1,W2)

    ENDDO   !kpoint
    ENDDO   !spin
    !=======================================================================

    DEALLOCATE(CRHOLM, GCHG, GWORK, CWORK)

    PROFILING_STOP('subtract_self_correlation')

   END SUBROUTINE SUBTRACT_SELF_CORRELATION

!***********************************************************************
!
! distributes the long-wave limit to the frequency groups 
! first, the global sum of HEAD, WING and CWING is done
! second, the distribution to the groups is performed.
!
!***********************************************************************
   
   SUBROUTINE DISTRIBUTE_HEAD_TO_CHI( HEADHAND, CHI ) 
      TYPE( head_handle)      :: HEADHAND
      TYPE( responsefunction) :: CHI
      !local
      INTEGER                 :: NG
      INTEGER                 :: IDIR 
      INTEGER                 :: NOMEGA_TOTAL 
      INTEGER                 :: NOMEGA_CURRENT
      INTEGER                 :: NCURRENT_GROUP
      INTEGER                 :: NOMEGA_LOCAL
      INTEGER                 :: NOMEGA_IN_GROUP
 
      !# of G-vectors 
      NG = HEADHAND%NROWS
      !# of dimensions
      IDIR = SIZE( HEADHAND%HEAD,1 )
      !# of frequency points 
      NOMEGA_TOTAL = HEADHAND%NOMEGA 

      ! do a global sum of HEAD, WING and CWING
      CALLMPI( M_sum_z(HEADHAND%COMM, HEADHAND%WING(1,1,1) , NG*IDIR*NOMEGA_TOTAL ))
      CALLMPI( M_sum_z(HEADHAND%COMM, HEADHAND%CWING(1,1,1), NG*IDIR*NOMEGA_TOTAL ))
      CALLMPI( M_sum_z(HEADHAND%COMM, HEADHAND%HEAD(1,1,1) , IDIR*IDIR*NOMEGA_TOTAL ))

      !--------------------------------------------------------------------------
      DO NOMEGA_CURRENT = 1 , NOMEGA_TOTAL   !loop over all global frequency points 
      !--------------------------------------------------------------------------
         ! determine frequency group and local frequency index
         NCURRENT_GROUP = DETERMINE_CURRENT_FREQ_GROUP(NOMEGA_CURRENT)
         IF ( NCURRENT_GROUP > HEADHAND%COMMBET%NCPU ) THEN
            CALL vtutor%bug("DISTRIBUTE_HEAD_TO_CHI: omega group not found " // str(HEADHAND%COMM%NODE_ME) // " " // str(NOMEGA_TOTAL) // " " // str(NCURRENT_GROUP), __FILE__, __LINE__)
         ENDIF 
         NOMEGA_LOCAL = DETERMINE_LOCAL_FREQ_INDEX( NOMEGA_CURRENT ) 
         IF ( NOMEGA_LOCAL > HEADHAND%DISTRIBUTION( NCURRENT_GROUP , 2 ) ) THEN
            CALL vtutor%bug("DISTRIBUTE_HEAD_TO_CHI: local omega index too large " // str(HEADHAND%COMM%NODE_ME) // " " // str(NCURRENT_GROUP) // " " // str(NOMEGA_LOCAL) // " " // str(HEADHAND%DISTRIBUTION(NCURRENT_GROUP,2)), __FILE__, __LINE__)
         ENDIF 
         
         NOMEGA_IN_GROUP = HEADHAND%DISTRIBUTION( NCURRENT_GROUP , 2 )

         ! only group members of group NCURRENT_GROUP save the result
         ! compensate for complex conjugation
         ! the chi_GG store the conjugated response function (compared to chi_base)
         IF ( HEADHAND%COMMBET%NODE_ME == NCURRENT_GROUP ) THEN
            IF (CHI%LREAL .AND. CHI%LREALSTORE) THEN
               CHI%HEAD(:,:,NOMEGA_LOCAL) = HEADHAND%HEAD( : , :, NOMEGA_CURRENT )
#ifdef gammareal
               CHI%WINGR(:,:,NOMEGA_LOCAL) = HEADHAND%WINGR( : , :, NOMEGA_CURRENT )
               CHI%CWINGR(:,:,NOMEGA_LOCAL)= HEADHAND%CWINGR( : , :, NOMEGA_CURRENT )
#else
               CHI%WING(:,:,NOMEGA_LOCAL) = HEADHAND%WING( : , :, NOMEGA_CURRENT )
               CHI%CWING(:,:,NOMEGA_LOCAL)= HEADHAND%CWING( : , :, NOMEGA_CURRENT )
#endif
            ELSE IF (CHI%LREAL) THEN
               WRITE(0,*) 'internal error in DISTRIBUTE_HEAD_TO_CHI: CHI%LREAL present not supported'
            ELSE
               CHI%HEAD(:,:,NOMEGA_LOCAL) = CONJG(HEADHAND%HEAD( : , :, NOMEGA_CURRENT ))  
#ifdef gammareal
               CHI%WINGR(:,:,NOMEGA_LOCAL) = HEADHAND%WINGR( : , :, NOMEGA_CURRENT )
               CHI%CWINGR(:,:,NOMEGA_LOCAL)= HEADHAND%CWINGR( : , :, NOMEGA_CURRENT )
#else
               CHI%WING(:,:,NOMEGA_LOCAL) = CONJG(HEADHAND%WING( : , :, NOMEGA_CURRENT ))   
               CHI%CWING(:,:,NOMEGA_LOCAL)= CONJG(HEADHAND%CWING( : , :, NOMEGA_CURRENT )) 
#endif 
            ENDIF
         ENDIF 
      !--------------------------------------------------------------------------
      ENDDO
      !--------------------------------------------------------------------------

   END SUBROUTINE DISTRIBUTE_HEAD_TO_CHI

!***********************************************************************
! allocation routine for head_handle
!***********************************************************************
 
   SUBROUTINE ALLOCATE_HEAD_HANDLE( HEADHAND, CHI, IMAG_GRIDS )
      USE minimax_struct, ONLY: imag_grid_handle
      IMPLICIT NONE 
      TYPE(head_handle)      :: HEADHAND   !handle to be initialized
      TYPE(responsefunction) :: CHI 
      TYPE(imag_grid_handle) :: IMAG_GRIDS !imaginary grids
      !local 
      TYPE( communic )       :: GLOBCOMM !global communicator
      TYPE( communic )       :: COMMIN   !communicator in frequency group
      TYPE( communic )       :: COMMBET  !communicator in frequency group
      INTEGER                :: IDIR

      GLOBCOMM = IMAG_GRIDS%B%COMM
      COMMIN = IMAG_GRIDS%B%COMM_IN_GROUP
      COMMBET = IMAG_GRIDS%B%COMM_BETWEEN_GROUPS

      !this routine works only if the array IMAG_GRIDS%B%DISTRIBUTION is allocated
      IF ( .NOT. ( ASSOCIATED( IMAG_GRIDS%B%DISTRIBUTION) ) ) THEN 
         CALL vtutor%bug("ALLOCATE_HEAD_HANDLE, B%DISTRIBUTION is not allocated", __FILE__, __LINE__)
      ENDIF 
      !small consistency check for DISTRIBUTION
      IF ( COMMBET%NCPU /= SIZE( IMAG_GRIDS%B%DISTRIBUTION, 1 ) ) THEN
         CALL vtutor%bug("ALLOCATE_HEAD_HANDLE, # of frequency groups inconsistent with DISTRIBUTION " // str(COMMBET%NCPU) // " " // str(SIZE(IMAG_GRIDS%B%DISTRIBUTION,1)), __FILE__, __LINE__ )
      ENDIF 

      !number of (complex) rows of response matrix
      HEADHAND%NROWS  = CHI%NP1
      !total number of frequency points
      HEADHAND%NOMEGA = IMAG_GRIDS%B%NPOINTS
      !set global communicator
      HEADHAND%COMM   = GLOBCOMM
      !communicator in frequency group
      HEADHAND%COMMIN = COMMIN
      !communicator between frequency groups
      HEADHAND%COMMBET= COMMBET
      
      !# of reciprocal lattice vectors
      !default, 3 dimensions
      IDIR = 3 
      !allocate head, wing and cwing
      !according to NROWS and NOMEGA
      ALLOCATE( HEADHAND%HEAD( IDIR, IDIR, IMAG_GRIDS%B%NPOINTS ) )
      ALLOCATE( HEADHAND%WING( HEADHAND%NROWS, IDIR, IMAG_GRIDS%B%NPOINTS ) )  
      ALLOCATE( HEADHAND%CWING( HEADHAND%NROWS, IDIR, IMAG_GRIDS%B%NPOINTS ) )  
      ! poke complex array to real array with twice the size in case of gamma-only
      CALL SET_RESPONSER( HEADHAND%WINGR,  2*HEADHAND%NROWS, IDIR,  IMAG_GRIDS%B%NPOINTS, HEADHAND%WING(1,1,1))
      CALL SET_RESPONSER( HEADHAND%CWINGR, 2*HEADHAND%NROWS, IDIR,  IMAG_GRIDS%B%NPOINTS, HEADHAND%CWING(1,1,1))
      !initialize to zero 
      HEADHAND%HEAD = 0 
      HEADHAND%WING = 0 
      HEADHAND%CWING = 0 

      !allocate frequency array 
      ALLOCATE( HEADHAND%OMEGA( IMAG_GRIDS%B%NPOINTS ) ) 
      !and set points to OMEGA
      HEADHAND%OMEGA(1:IMAG_GRIDS%B%NPOINTS) = IMAG_GRIDS%BOS_RE(1:IMAG_GRIDS%B%NPOINTS)

      !also allocate DISTRIBUTION
      ALLOCATE( HEADHAND%DISTRIBUTION( COMMBET%NCPU, 2 ) ) 
      !and set to DISTRIBUTION
      HEADHAND%DISTRIBUTION = IMAG_GRIDS%B%DISTRIBUTION

   END SUBROUTINE ALLOCATE_HEAD_HANDLE

!***********************************************************************
! deallocation routine for head_handle
!***********************************************************************
 
   SUBROUTINE DEALLOCATE_HEAD_HANDLE( HEADHAND ) 
      TYPE( head_handle) :: HEADHAND

      !deallocate head, wing and cwing
      DEALLOCATE( HEADHAND%HEAD )
      DEALLOCATE( HEADHAND%WING )  
      DEALLOCATE( HEADHAND%CWING)  
      DEALLOCATE( HEADHAND%OMEGA) 
      !and nullify
      NULLIFY( HEADHAND%WING, HEADHAND%CWING, HEADHAND%HEAD )
      !also nullify OMEGA and DISTRIBUTION
      NULLIFY( HEADHAND%OMEGA , HEADHAND%DISTRIBUTION ) 
      
   END SUBROUTINE DEALLOCATE_HEAD_HANDLE

!***********************************************************************
!
! computes the screened Coulomb potential W in the imaginary frequency
! domain. On entry CHI contains the response function, on exit the 
! screened potential W
!  
!***********************************************************************
SUBROUTINE COMPUTE_SCREENED_POTENTIAL( WDES, CHI, GDES, WGWQ, IMAG_GRIDS,&
    LATT_CUR, NQ , IDIRMAX, FSG0, IO)
    USE base
    USE constant
    USE mpimy
    USE ini
    USE acfdt_GG 
    USE minimax_struct, ONLY : imag_grid_handle

    IMPLICIT NONE
    TYPE (wavedes)          :: WDES
    TYPE (responsefunction) :: CHI                  
    TYPE (greensfdes)       :: GDES                 ! Green function in freq. domain
    TYPE (wavedes1)         :: WGWQ
    TYPE (imag_grid_handle) :: IMAG_GRIDS           !imaginary grids 
    TYPE (latt)             :: LATT_CUR
    INTEGER                 :: NQ                   !current q-point
    INTEGER                 :: IDIRMAX             !maximum number of directions
    REAL(q)                 :: FSG0                 !convergence correction for bare Coulomb kernel
    TYPE (in_struct)        :: IO                 
    ! local 
    INTEGER                 :: NROWS, NCOLS         !number of local rows and cols of CHI
    TYPE (communic)         :: COMM_BETWEENOMEGA    !communicator between nu groups
    TYPE (communic)         :: COMM_INOMEGA         !communicator in nu groups
    INTEGER                 :: NOMEGA_SIMULTANEOUS  !# of freq. points on one core

    INTEGER                 :: NLOCAL
    INTEGER                 :: IUNIT
    INTEGER IDIR, NOMEGA_CURRENT

    NROWS = GDES%RES_NRPLWV_ROW_DATA_POINTS
    NCOLS = GDES%RES_NRPLWV_COL_DATA_POINTS
    COMM_BETWEENOMEGA = IMAG_GRIDS%B%COMM_BETWEEN_GROUPS
    COMM_INOMEGA = IMAG_GRIDS%B%COMM_IN_GROUP
    NOMEGA_SIMULTANEOUS = IMAG_GRIDS%B%NPOINTS_IN_GROUP

    IF (CHI%LREAL .AND. .NOT. CHI%LREALSTORE) THEN
       CALL vtutor%bug("BUILD_W_FROM_INVERSE_EPSILON: for the Gamma point version CHI%LREALSTORE must be set", __FILE__, __LINE__)
    ENDIF
    IF (CHI%LREAL .AND. .NOT. WGWQ%LGAMMA) THEN
       CALL vtutor%bug("BUILD_W_FROM_INVERSE_EPSILON: for the Gamma point version WGWQ%LGAMMA must be set", __FILE__, __LINE__ )
    ENDIF

    PROFILING_START('compute_screened_potential')

    !Maybe TODO: it would be wiser to loop over frequencies outside of the routines below
    !            this would save at least two loops. 
    !            On the other side, the frequency loops are negligible compared to 
    !            matrix inversion or matrix multiplication. 

    !first compute reduced polarizability  X_red = X_0.(1-X_0.V)^-1
    !this is the most expensive step due to the inversion of (1-X_0.V) -> scaLAPACK is used
    !the processor grid is chosen in the same way as for the ACFDT routine 
    IF (.NOT. L2ORDER) THEN
       CALL COMPUTE_REDUCED_CHI( WDES, CHI, NROWS, NCOLS, WGWQ, COMM_BETWEENOMEGA, COMM_INOMEGA, &
            LATT_CUR, IMAG_GRIDS%B%NPOINTS, NOMEGA_SIMULTANEOUS, NQ , IDIRMAX, IO)
       CALL STOP_TIMING("G",IO%IU6,"XI_LOCAL_FIELD")
    ENDIF

    !dump response after this step 
    CALL DUMP_HEAD_GG( CHI, IMAG_GRIDS, &
       "1 + v P,  with REDUCIBLE POLARIZABILTY P=P_0 (1 -(v+f) P_0)^-1",& 
       IO%IU6, IO%NWRITE, .TRUE., EDEPS/LATT_CUR%OMEGA, 1.0_q)

    ! determine inverse reduced epsilon  eps^-1 = 1 + v X_red 
    ! NOTE currently   eps^-1 = v X_red 
    CALL BUILD_INVERSE_EPSILON( WDES, CHI, NROWS, NCOLS, WGWQ, COMM_BETWEENOMEGA, &
       COMM_INOMEGA, LATT_CUR, IMAG_GRIDS%B%NPOINTS, NOMEGA_SIMULTANEOUS, NQ , IDIRMAX, IO)

    !todo: dump diagonal of CHI to vasprun.xml    
    CALL DUMP_HEAD_GG( CHI, IMAG_GRIDS, & 
      "INVERSE MACROSCOPIC DIELECTRIC TENSOR (including local field effects in RPA (Hartree))",& 
      IO%IU6, IO%NWRITE, .TRUE.)

    ! determine W = eps^-1 v
    CALL BUILD_W_FROM_INVERSE_EPSILON( WDES, CHI, NROWS, NCOLS, WGWQ, COMM_BETWEENOMEGA, &
         COMM_INOMEGA, LATT_CUR, IMAG_GRIDS%B%NPOINTS, NOMEGA_SIMULTANEOUS, NQ , IDIRMAX, IO)

    ! apply convergence correction to W and subtract bare coloumb kernel
    CALL CONVERGENCE_CORRECTION_FOR_W( CHI, NROWS, NCOLS, WGWQ, WDES%COMM_KIN, COMM_INOMEGA, &
         COMM_BETWEENOMEGA, LATT_CUR, NOMEGA_SIMULTANEOUS, FSG0, &
         LFOCK_SUBTRACT, IO)

    !time used for W
    CALL STOP_TIMING("G",IO%IU6,"XI_TO_W")

    !dump response after this step 
    CALL DUMP_HEAD_GG( CHI, IMAG_GRIDS, "screened Coulomb potential",& 
      IO%IU6, IO%NWRITE)
   

#ifdef debug
!dump final W
    DO NLOCAL = 1,1 ! NOMEGA_SIMULTANEOUS
       IUNIT=COMM_BETWEENOMEGA%NODE_ME*300+ COMM_INOMEGA%NODE_ME
#ifndef gammareal
       WRITE(IUNIT,*)' This is W',NLOCAL
       DO I=1,MIN(12,NROWS)
       WRITE(IUNIT,'(12F8.4)')REAL(CHI%RESPONSEFUN(I,1:MIN(12,NCOLS),NLOCAL),q)
       ENDDO
       WRITE(IUNIT,*)' ' 
       WRITE(IUNIT,*)' ' 
       DO I=1,MIN(12,NROWS)
       WRITE(IUNIT,'(12F8.4)')IMAG(CHI%RESPONSEFUN(I,1:MIN(12,NCOLS),NLOCAL))
       ENDDO
       WRITE(IUNIT,*)' ' 
#else
       DO I=1,MIN(12,NROWS)
       WRITE(IUNIT,'(12F8.4)')REAL(CHI%RESPONSER(I,1:MIN(12,NCOLS),NLOCAL),q)
       ENDDO
       WRITE(IUNIT,*)' ' 
#endif
    ENDDO
#endif
    PROFILING_STOP('compute_screened_potential')

END SUBROUTINE COMPUTE_SCREENED_POTENTIAL

    
!***********************************************************************
!
!  computes 
!                                  -1
!  X_red= X_0 (1- f_xc X_0 - v X_0) 
!
!  with scaLAPACK
!  
!***********************************************************************

SUBROUTINE COMPUTE_REDUCED_CHI( WDES, CHI, NROWS, NCOLS, WGWQ, &
    COMM_BETWEENOMEGA, COMM_INOMEGA, LATT_CUR, NOMEGA_TOTAL, &
    NOMEGA_SIMULTANEOUS, NQ , IDIRMAX, IO)
    USE base
    USE constant
    USE mpimy
    USE ini
    USE acfdt_GG 

    IMPLICIT NONE
    TYPE (wavedes)          :: WDES
    TYPE (responsefunction) :: CHI                  
    INTEGER                 :: NROWS, NCOLS         !number of local rows and cols of CHI
    TYPE (wavedes1)         :: WGWQ
    TYPE (communic)         :: COMM_BETWEENOMEGA    !communicator between nu groups
    TYPE (communic)         :: COMM_INOMEGA         !communicator in nu groups
    TYPE (latt)             :: LATT_CUR
    INTEGER                 :: NOMEGA_TOTAL         !number of frequency points
    INTEGER                 :: NOMEGA_SIMULTANEOUS  !# of freq. points on one core
    INTEGER                 :: NQ                   !current q-point
    INTEGER                 :: IDIRMAX             !maximum number of directions
    TYPE (in_struct)        :: IO                 
  ! local
    GDEF, POINTER, DIMENSION(:,:) :: CHI_WORK=>NULL()
    GDEF, POINTER, DIMENSION(:,:) :: EPS_WORK=>NULL()
    INTEGER :: NP, I, NCURR, IDIR
    INTEGER :: NOMEGA_GLOBAL
    INTEGER :: NOMEGA_IN_ROOT_GROUP                   ! # of frequencies in root group
    INTEGER :: NOMEGA_CURRENT
    LOGICAL :: LDMODE
    
    IF (CHI%LREAL .AND. .NOT. CHI%LREALSTORE) THEN
       CALL vtutor%bug("COMPUTE_REDUCED_CHI: \n or the Gamma point version CHI%LREALSTORE must be set", __FILE__, __LINE__)
    ENDIF
    
    !each group needs to loop over the same # of frequencies internally,
    !otherwise when initializing BLACS in LNTRACE_OF_CHI the routine 
    !BLACS_PINFO in GEN_PROC_GROUP_GRIDS_BC waits for the CPUs in other group
    !to call itself. 
    !since the root group contains always the maximum # of frequencies we 
    !use NOMEGA_SIMULTANEOUS from the root node for the internal # of frequencies for all groups
    NOMEGA_IN_ROOT_GROUP = 0
    IF ( WDES%COMM_KIN%NODE_ME == 1 ) NOMEGA_IN_ROOT_GROUP = NOMEGA_SIMULTANEOUS
    CALLMPI( M_bcast_i( WDES%COMM_KIN, NOMEGA_IN_ROOT_GROUP, 1 ) ) 

    !scaLAPACK block size and processor grid is set here
    IF ( CHI%NQ == 1 ) THEN
       CALL BLOCKSIZE_AND_PROC_GRID( NROWS, NCOLS, COMM_INOMEGA%NCPU, -1)
    ENDIF 

    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2
    ALLOCATE(CHI_WORK(NROWS,NCOLS))
    !--------------------------------------------------------------------------------------
    omega_sim: DO NCURR = 1, NOMEGA_IN_ROOT_GROUP  !loop over frequencies within group
    !--------------------------------------------------------------------------------------
       IF ( NCURR <= NOMEGA_SIMULTANEOUS ) THEN
          !this is the global omega point treated currently in each group 
          NOMEGA_GLOBAL = LOCAL_INDEX_TO_GLOBAL( NCURR, COMM_BETWEENOMEGA%NODE_ME, &
              COMM_BETWEENOMEGA%NCPU, NOMEGA_TOTAL)
          NOMEGA_CURRENT = NCURR
          LDMODE = .FALSE.
       ELSE
          !group having less points than the root group use their largest local frequency
          NOMEGA_GLOBAL = LOCAL_INDEX_TO_GLOBAL( NOMEGA_SIMULTANEOUS, COMM_BETWEENOMEGA%NODE_ME, &
              COMM_BETWEENOMEGA%NCPU, NOMEGA_TOTAL)
          NOMEGA_CURRENT = NOMEGA_SIMULTANEOUS
          LDMODE = .TRUE. 
       ENDIF

       !-----------------------------------------------------------------------------------
       IF ( CHI%LGAMMA ) THEN  !gamma point is treated separately  
          ! diagonalize the matrix for three directions
          ! this could be reworked using the block diagonalization
          ! of Baroni and Resta Phys. Rev. B 33, 7017 (1986).
          ! unfortunately an additional working array needs to be allocated here
          ALLOCATE(EPS_WORK(NROWS,NCOLS)) 
          EPS_WORK = zero
          DO IDIR=1,IDIRMAX
             !for the long-wave limit we copy head and wings to RESPONSEFUN
             CALL BODY_FROM_WING_GG( CHI, IDIR, NOMEGA_CURRENT, NROWS, NCOLS, COMM_INOMEGA)
             !initialize CHI_WORK here 
             CHI_WORK = zero
  
             !build CHI_WORK <-- (1-X*V)      
             CALL BUILD_RED_CHI_OR_EPSILON(CHI_WORK, CHI, NOMEGA_CURRENT, NROWS, NCOLS, &
                 COMM_INOMEGA,WGWQ, IO%IU0, .TRUE. ) 
             !use scaLAPACK to compute the following expression
             !                                      -1
             !CHI_WORK <-- X_0 (1- f_xc X_0 - v X_0)     
             CALL INVERT_IRRED_CHI(CHI, CHI_WORK, NROWS,NCOLS, NP,COMM_INOMEGA,COMM_BETWEENOMEGA, &
                 WDES%COMM_KIN, NOMEGA_CURRENT, NOMEGA_GLOBAL, NOMEGA_SIMULTANEOUS, IO,           &
                 CHI%LREALSTORE, LDMODE , CHI%LGAMMA, IDIRMAX ) 
             !add the IDIRMAX directions to EPS_WORK   
             EPS_WORK=EPS_WORK+CHI_WORK
             ! set head and wing to calculated values 
             IF (.NOT. LDMODE) THEN
                CALL SET_WING_FROM_MAT_GG( CHI_WORK, CHI, COMM_INOMEGA, NOMEGA_CURRENT, LDMODE, IDIR)
             ENDIF
          ENDDO  !directions loop for gamma point
          EPS_WORK=EPS_WORK*(1.0_q/IDIRMAX)
          !set body, head and wings from EPS_WORK
          IF (.NOT. LDMODE) CALL SET_RESPONSE_FROM_MAT_GG( EPS_WORK, CHI, NOMEGA_CURRENT, LDMODE)
          !get rid of EPS_WORK
          DEALLOCATE(EPS_WORK)
       ELSE
          !for k-point different from gamma it's simpler
          CHI_WORK=zero
          !build CHI_WORK <-- (1-X*V)      
          CALL BUILD_RED_CHI_OR_EPSILON(CHI_WORK, CHI, NOMEGA_CURRENT, NROWS, NCOLS, &
              COMM_INOMEGA,WGWQ, IO%IU0, .TRUE. ) 
          !use scaLAPACK to compute the following expression
          !                                      -1
          !CHI_WORK <-- X_0 (1- f_xc X_0 - v X_0)     
          CALL INVERT_IRRED_CHI(CHI, CHI_WORK, NROWS, NCOLS, NP, COMM_INOMEGA, COMM_BETWEENOMEGA, &
              WDES%COMM_KIN, NOMEGA_CURRENT, NOMEGA_GLOBAL, NOMEGA_SIMULTANEOUS, IO,           &
              CHI%LREALSTORE, LDMODE , CHI%LGAMMA, IDIRMAX ) 
          !set body
          IF (.NOT. LDMODE) CALL SET_RESPONSE_FROM_MAT_GG( CHI_WORK, CHI, NOMEGA_CURRENT, LDMODE)

       !-----------------------------------------------------------------------------------
       ENDIF
    !--------------------------------------------------------------------------------------
    ENDDO omega_sim
    !--------------------------------------------------------------------------------------
    DEALLOCATE ( CHI_WORK)     

    !restore head finally 
    CALL RESTORE_HEAD(CHI) 
 RETURN
END SUBROUTINE COMPUTE_REDUCED_CHI
    
!***********************************************************************
!
! this is the adopted version of XI_RED_TO_EPS from chi.F
! on entry CHI is the reduced polarizability X_red
! on exit CHI  is the inverse epsilon 
!     -1
!  eps   = (1 + v X_red )
!
!***********************************************************************

SUBROUTINE BUILD_INVERSE_EPSILON( WDES, CHI, NROWS, NCOLS, WGWQ, &
    COMM_BETWEENOMEGA, COMM_INOMEGA, LATT_CUR, NOMEGA_TOTAL, &
    NOMEGA_SIMULTANEOUS, NQ , IDIRMAX, IO)
    USE base
    USE constant
    USE mpimy
    USE ini
    USE acfdt_GG 

    IMPLICIT NONE
    TYPE (wavedes)          :: WDES
    TYPE (responsefunction) :: CHI                  
    INTEGER                 :: NROWS, NCOLS         !number of local rows and cols of CHI
    TYPE (wavedes1)         :: WGWQ
    TYPE (communic)         :: COMM_BETWEENOMEGA    !communicator between nu groups
    TYPE (communic)         :: COMM_INOMEGA         !communicator in nu groups
    TYPE (latt)             :: LATT_CUR
    INTEGER                 :: NOMEGA_TOTAL         !number of frequency points
    INTEGER                 :: NOMEGA_SIMULTANEOUS  !# of freq. points on one core
    INTEGER                 :: NQ                   !current q-point
    INTEGER                 :: IDIRMAX             !maximum number of directions
    TYPE (in_struct)        :: IO                 
  ! local
    GDEF, POINTER, DIMENSION(:,:) :: CHI_WORK=>NULL()
    GDEF, POINTER, DIMENSION(:,:) :: EPS_WORK=>NULL()
    INTEGER :: NP, NCURR, IDIR
    INTEGER :: NOMEGA_GLOBAL
    INTEGER :: NOMEGA_IN_ROOT_GROUP                   ! # of frequencies in root group
    INTEGER :: NOMEGA_CURRENT
    LOGICAL :: LDMODE
    
    IF (CHI%LREAL .AND. .NOT. CHI%LREALSTORE) THEN
       CALL vtutor%bug("BUILD_INVERSE_EPSILON: \n or the Gamma point version CHI%LREALSTORE must be set", __FILE__, __LINE__ )
    ENDIF
    
    !each group needs to loop over the same # of frequencies internally,
    !otherwise when initializing BLACS in LNTRACE_OF_CHI the routine 
    !BLACS_PINFO in GEN_PROC_GROUP_GRIDS_BC waits for the CPUs in other group
    !to call itself. 
    !since the root group contains always the maximum # of frequencies we 
    !use NOMEGA_SIMULTANEOUS from the root node for the internal # of frequencies for all groups
    NOMEGA_IN_ROOT_GROUP = 0
    IF ( WDES%COMM_KIN%NODE_ME == 1 ) NOMEGA_IN_ROOT_GROUP = NOMEGA_SIMULTANEOUS
    CALLMPI( M_bcast_i( WDES%COMM_KIN, NOMEGA_IN_ROOT_GROUP, 1 ) ) 

    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2
    !we need a working array 
    ALLOCATE(CHI_WORK(NROWS,NCOLS)) 
 
    !--------------------------------------------------------------------------------------
    omega_sim: DO NCURR = 1, NOMEGA_IN_ROOT_GROUP  !loop over frequencies within group
    !--------------------------------------------------------------------------------------
       IF ( NCURR <= NOMEGA_SIMULTANEOUS ) THEN
          !this is the global omega point treated currently in each group 
          NOMEGA_GLOBAL = LOCAL_INDEX_TO_GLOBAL( NCURR, COMM_BETWEENOMEGA%NODE_ME, &
              COMM_BETWEENOMEGA%NCPU, NOMEGA_TOTAL)
          NOMEGA_CURRENT = NCURR
          LDMODE = .FALSE.
       ELSE
          !group having less points than the root group use their largest local frequency
          NOMEGA_GLOBAL = LOCAL_INDEX_TO_GLOBAL( NOMEGA_SIMULTANEOUS, COMM_BETWEENOMEGA%NODE_ME, &
              COMM_BETWEENOMEGA%NCPU, NOMEGA_TOTAL)
          NOMEGA_CURRENT = NOMEGA_SIMULTANEOUS
          LDMODE = .TRUE. 
       ENDIF
       
       ! if group is in dummy mode, do nothing    
       IF ( LDMODE ) CYCLE
       !-----------------------------------------------------------------------------------
       IF ( CHI%LGAMMA ) THEN  !gamma point is treated separately  
          ! unfortunately an additional working array needs to be allocated here
          ALLOCATE(EPS_WORK(NROWS,NCOLS)) 
          EPS_WORK = zero
          DO IDIR=1,IDIRMAX
             !for the long-wave limit we copy head and wings to RESPONSEFUN
             CALL BODY_FROM_WING_GG( CHI, IDIR, NOMEGA_CURRENT, NROWS, NCOLS, COMM_INOMEGA)
             !initialize CHI_WORK here 
             CHI_WORK = zero

             !build CHI_WORK <-- (1-V.*CHI)       
             CALL BUILD_RED_CHI_OR_EPSILON(CHI_WORK, CHI, NOMEGA_CURRENT, NROWS, NCOLS, &
                 COMM_INOMEGA,WGWQ, IO%IU0, .FALSE. ) 

             ! second entry is the sin transform corresponding to G=0, always zero
             IF (CHI%LREAL.AND.COMM_INOMEGA%NODE_ME==1) CHI_WORK(2,2)=0
             !add the IDIRMAX directions to EPS_WORK   
             EPS_WORK=EPS_WORK+CHI_WORK
             ! set head and wing to calculated values
             CALL SET_WING_FROM_MAT_GG( CHI_WORK, CHI, COMM_INOMEGA, NOMEGA_CURRENT,  LDMODE, IDIR)
          ENDDO  !directions loop for gamma point
          EPS_WORK=EPS_WORK*(1.0_q/IDIRMAX)
          !set body, head and wings from EPS_WORK
          CALL SET_RESPONSE_FROM_MAT_GG( EPS_WORK, CHI, NOMEGA_CURRENT, LDMODE)
          DEALLOCATE(EPS_WORK)
       ELSE
          !for k-point different from gamma it's simpler
          CHI_WORK=zero
          !build CHI_WORK <-- (1-V*CHI)      
          CALL BUILD_RED_CHI_OR_EPSILON(CHI_WORK, CHI, NOMEGA_CURRENT, NROWS, NCOLS, &
              COMM_INOMEGA,WGWQ, IO%IU0, .FALSE.) 
          !set body
          CALL SET_RESPONSE_FROM_MAT_GG( CHI_WORK, CHI, NOMEGA_CURRENT, LDMODE)
       !-----------------------------------------------------------------------------------
       ENDIF
    !--------------------------------------------------------------------------------------
    ENDDO omega_sim
    !--------------------------------------------------------------------------------------
    DEALLOCATE ( CHI_WORK)     
    
    !ok, now restore head finally 
    CALL RESTORE_HEAD(CHI) 
END SUBROUTINE BUILD_INVERSE_EPSILON

!***********************************************************************
! 
!  computes W form inverse epsilon according to
!
!         -1
!  W = eps   v
!
!***********************************************************************

SUBROUTINE BUILD_W_FROM_INVERSE_EPSILON( WDES, CHI, NROWS, NCOLS, WGWQ, &
    COMM_BETWEENOMEGA, COMM_INOMEGA, LATT_CUR, NOMEGA_TOTAL, &
    NOMEGA_SIMULTANEOUS, NQ , IDIRMAX, IO)
    USE base
    USE constant
    USE mpimy
    USE ini
    USE acfdt_GG 

    IMPLICIT NONE
    TYPE (wavedes)          :: WDES
    TYPE (responsefunction) :: CHI                  
    INTEGER                 :: NROWS, NCOLS         !number of local rows and cols of CHI
    TYPE (wavedes1)         :: WGWQ
    TYPE (communic)         :: COMM_BETWEENOMEGA    !communicator between nu groups
    TYPE (communic)         :: COMM_INOMEGA         !communicator in nu groups
    TYPE (latt)             :: LATT_CUR
    INTEGER                 :: NOMEGA_TOTAL         !number of frequency points
    INTEGER                 :: NOMEGA_SIMULTANEOUS  !# of freq. points on one core
    INTEGER                 :: NQ                   !current q-point
    INTEGER                 :: IDIRMAX             !maximum number of directions
    TYPE (in_struct)        :: IO                 
  ! local
    INTEGER :: NP, NCURR,  NI, NI_
    INTEGER :: NOMEGA_GLOBAL
    INTEGER :: NOMEGA_IN_ROOT_GROUP                   ! # of frequencies in root group
    INTEGER :: NOMEGA_CURRENT
    LOGICAL :: LDMODE
    
    !each group needs to loop over the same # of frequencies internally,
    !otherwise when initializing BLACS in LNTRACE_OF_CHI the routine 
    !BLACS_PINFO in GEN_PROC_GROUP_GRIDS_BC waits for the CPUs in other group
    !to call itself. 
    !since the root group contains always the maximum # of frequencies we 
    !use NOMEGA_SIMULTANEOUS from the root node for the internal # of frequencies for all groups
    NOMEGA_IN_ROOT_GROUP = 0
    IF ( WDES%COMM_KIN%NODE_ME == 1 ) NOMEGA_IN_ROOT_GROUP = NOMEGA_SIMULTANEOUS
    CALLMPI( M_bcast_i( WDES%COMM_KIN, NOMEGA_IN_ROOT_GROUP, 1 ) ) 

    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2
 
    !--------------------------------------------------------------------------------------
    omega_sim: DO NCURR = 1, NOMEGA_IN_ROOT_GROUP  !loop over frequencies within group
    !--------------------------------------------------------------------------------------
       IF ( NCURR <= NOMEGA_SIMULTANEOUS ) THEN
          !this is the global omega point treated currently in each group 
          NOMEGA_GLOBAL = LOCAL_INDEX_TO_GLOBAL( NCURR, COMM_BETWEENOMEGA%NODE_ME, &
              COMM_BETWEENOMEGA%NCPU, NOMEGA_TOTAL)
          NOMEGA_CURRENT = NCURR
          LDMODE = .FALSE.
       ELSE
          !group having less points than the root group use their largest local frequency
          NOMEGA_GLOBAL = LOCAL_INDEX_TO_GLOBAL( NOMEGA_SIMULTANEOUS, COMM_BETWEENOMEGA%NODE_ME, &
              COMM_BETWEENOMEGA%NCPU, NOMEGA_TOTAL)
          NOMEGA_CURRENT = NOMEGA_SIMULTANEOUS
          LDMODE = .TRUE. 
       ENDIF
       
       ! if group is in dummy mode, do nothing    
       IF ( LDMODE ) CYCLE

       !-----------------------------------------------------------------------------------
       DO NI=1,NCOLS
          ! get global index into Coulomb kernel
          NI_=(COMM_INOMEGA%NODE_ME-1)*NCOLS + NI 
          IF ( NI_ > SIZE(WGWQ%DATAKE,1) ) THEN
             WRITE(0,*)'internal error in VASP:  BUILD_W_FROM_INVERSE_EPSILON exceeding DATAKE matrix dimension', &
             COMM_INOMEGA%NODE_ME,WGWQ%COMM%NODE_ME,NI_, SIZE(WGWQ%DATAKE,1)
          ENDIF

         IF (CHI%LREALSTORE) THEN
          CHI%RESPONSER(1:NP,NI,NOMEGA_CURRENT)=CHI%RESPONSER(1:NP,NI,NOMEGA_CURRENT)*WGWQ%DATAKE(NI_,1)
          IF ( CHI%LGAMMA ) THEN
             CHI%CWINGR(NI,1,NOMEGA_CURRENT)=CHI%CWINGR(NI,1,NOMEGA_CURRENT)*WGWQ%DATAKE(NI_,1)
             CHI%CWINGR(NI,2,NOMEGA_CURRENT)=CHI%CWINGR(NI,2,NOMEGA_CURRENT)*WGWQ%DATAKE(NI_,1)
             CHI%CWINGR(NI,3,NOMEGA_CURRENT)=CHI%CWINGR(NI,3,NOMEGA_CURRENT)*WGWQ%DATAKE(NI_,1)
          ENDIF
         ELSE
          CHI%RESPONSEFUN(1:NP,NI,NOMEGA_CURRENT)=CHI%RESPONSEFUN(1:NP,NI,NOMEGA_CURRENT)*WGWQ%DATAKE(NI_,1)
          IF ( CHI%LGAMMA ) THEN
             CHI%CWING(NI,1,NOMEGA_CURRENT)=CHI%CWING(NI,1,NOMEGA_CURRENT)*WGWQ%DATAKE(NI_,1)
             CHI%CWING(NI,2,NOMEGA_CURRENT)=CHI%CWING(NI,2,NOMEGA_CURRENT)*WGWQ%DATAKE(NI_,1)
             CHI%CWING(NI,3,NOMEGA_CURRENT)=CHI%CWING(NI,3,NOMEGA_CURRENT)*WGWQ%DATAKE(NI_,1)
          ENDIF
         ENDIF
       ENDDO

       CHI%WING(:,:,NOMEGA_CURRENT) = CHI%WING(:,:,NOMEGA_CURRENT)*WGWQ%DATAKE(1,1)
    !--------------------------------------------------------------------------------------
    ENDDO omega_sim
    !--------------------------------------------------------------------------------------
 RETURN
END SUBROUTINE BUILD_W_FROM_INVERSE_EPSILON

!***********************************************************************
!
!  adds a convergence correction to the screened potential 
!  (in the same way as in XI_COUlOMB) 
!
!***********************************************************************

SUBROUTINE CONVERGENCE_CORRECTION_FOR_W(CHI, NROWS, NCOLS, WGWQ, GLOBCOMM, COMM_INOMEGA, &
    COMM_BETWEENOMEGA, LATT_CUR, NOMEGA_SIMULTANEOUS, FSG0, LFOCKSUBTRACT, IO)
    USE base
    USE constant
    USE mpimy
    USE ini
    USE acfdt_GG 

    IMPLICIT NONE
    TYPE (responsefunction) :: CHI                  
    INTEGER                 :: NROWS, NCOLS         !number of local rows and cols of CHI
    TYPE (wavedes1)         :: WGWQ
    TYPE (communic)         :: GLOBCOMM             !global communicator
    TYPE (communic)         :: COMM_INOMEGA         !communicator in nu groups
    TYPE (communic)         :: COMM_BETWEENOMEGA    !communicator between nu groups
    TYPE (latt)             :: LATT_CUR
    INTEGER                 :: NOMEGA_SIMULTANEOUS  !# of freq. points on one core
    INTEGER                 :: NQ                   !total number of q-points
    REAL(q)                 :: FSG0                 !convergence correction for bare Coulomb kernel
    LOGICAL                 :: LFOCKSUBTRACT
    TYPE (in_struct)        :: IO                 
  ! local
    INTEGER :: NP, NCURR,  NI, NI_
    INTEGER :: NOMEGA_GLOBAL
    INTEGER :: NOMEGA_IN_ROOT_GROUP                   ! # of frequencies in root group
    INTEGER :: NOMEGA_CURRENT
    LOGICAL :: LDMODE
    TYPE (latt) LATT_EWALD
    REAL(q) :: B(3,3)
    REAL(q) :: SCALE1, SCALE2, SCALE3, FSG
    COMPLEX(q) :: INVEPS
    
    IF (CHI%LREAL .AND. .NOT. CHI%LREALSTORE) THEN
       CALL vtutor%bug("CONVERGENCE_CORRECTION_FOR_W requires real valued storage of response functions", __FILE__, __LINE__)
    ENDIF
#ifdef gammareal
    IF (.NOT. CHI%LREAL ) THEN
       CALL vtutor%bug("CONVERGENCE_CORRECTION_FOR_W requires CHI%LREAL to be set", __FILE__, __LINE__ )
    ENDIF
#endif

    !size of global screened Coulomb matrix 
    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2

    IF ( NP > NROWS ) THEN
       CALL vtutor%bug("CONVERGENCE_CORRECTION_FOR_W:, NP > NROWS " // str(NP) // " " // str(NROWS) // " " // str(COMM_INOMEGA%NODE_ME), __FILE__, __LINE__)
    ENDIF


    !maybe a subtraction of the bare Coulomb potential is wanted
    IF (LFOCKSUBTRACT) THEN
       DO NOMEGA_CURRENT = 1, NOMEGA_SIMULTANEOUS
#ifdef onlyHF
          CHI%RESPONSEFUN(:,:,NOMEGA_CURRENT)=0
#endif
          DO NI = 1, NCOLS
             ! determine row index
             NI_=(COMM_INOMEGA%NODE_ME-1)*NCOLS + NI 
             IF ( NI_ > NROWS .OR. NI_ < 0) THEN
                WRITE(0,'(A,4I6)')'internal ERROR in VASP: CONVERGENCE_CORRECTION_FOR_W, NI_ > NROWS',NI_,NROWS,&
                   NP, COMM_INOMEGA%NODE_ME
             ENDIF
             IF ( NI_ > SIZE(WGWQ%DATAKE,1) ) THEN
                WRITE(0,*)'internal error in VASP: CONVERGENCE_CORRECTION_FOR_W exceeding DATAKE matrix dimension', &
                COMM_INOMEGA%NODE_ME,WGWQ%COMM%NODE_ME,NI_, SIZE(WGWQ%DATAKE,1)
             ENDIF

             IF ( CHI%LREALSTORE) THEN
                !subtract the bare Coulomb operator 
                CHI%RESPONSER(NI_,NI,NOMEGA_CURRENT) = CHI%RESPONSER(NI_,NI,NOMEGA_CURRENT) - &
                     WGWQ%DATAKE(NI_,1) 
             ELSE
                CHI%RESPONSEFUN(NI_,NI,NOMEGA_CURRENT) = CHI%RESPONSEFUN(NI_,NI,NOMEGA_CURRENT) - &
                     WGWQ%DATAKE(NI_,1) 
             ENDIF
#ifdef onlyHF
             CHI%RESPONSEFUN(NI_,NI,NOMEGA_CURRENT) = WGWQ%DATAKE(NI_,1) 
#endif
          ENDDO
       ENDDO
    ENDIF

    !gamma point treatment:
    ! get proper convergence corrections for the 1/r long range behaviour lim G-> 4 pi e^2/G^2
    IF (CHI%LGAMMA) THEN
       ! clear wings in the body of the matrix
       DO NOMEGA_CURRENT=1,NOMEGA_SIMULTANEOUS
          IF ( CHI%LREALSTORE) THEN
             IF ( COMM_INOMEGA%NODE_ME ==1 ) CHI%RESPONSER(1:NROWS,1,NOMEGA_CURRENT)=0._q
             CHI%RESPONSER(1,1:NCOLS,NOMEGA_CURRENT)=0._q
          ELSE
             IF ( COMM_INOMEGA%NODE_ME ==1 ) CHI%RESPONSEFUN(1:NROWS,1,NOMEGA_CURRENT)=(0._q,0._q)
             CHI%RESPONSEFUN(1,1:NCOLS,NOMEGA_CURRENT)=(0._q,0._q)
          ENDIF
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
             CALL vtutor%error("error in CONVERGENCE_CORRECTION_FOR_W: presently the GW routine only &
                &supports k-point meshes generated automatically")
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
          FSG0=FSG0*2*NQ
       ELSE
          FSG0=HFRCUT*HFRCUT/2*EDEPS/LATT_CUR%OMEGA
       ENDIF

       !each group needs to loop over the same # of frequencies internally,
       !since the root group contains always the maximum # of frequencies we 
       !use NOMEGA_SIMULTANEOUS from the root node for the internal # of frequencies for all groups
       NOMEGA_IN_ROOT_GROUP = 0
       IF ( GLOBCOMM%NODE_ME == 1 ) NOMEGA_IN_ROOT_GROUP = NOMEGA_SIMULTANEOUS
       CALLMPI( M_bcast_i( GLOBCOMM, NOMEGA_IN_ROOT_GROUP, 1 ) ) 

       !--------------------------------------------------------------------------------------
       omega_sim: DO NCURR = 1, NOMEGA_IN_ROOT_GROUP  !loop over frequencies within group
       !--------------------------------------------------------------------------------------
          IF ( NCURR <= NOMEGA_SIMULTANEOUS ) THEN
             !this is the global omega point treated currently in each group 
             NOMEGA_GLOBAL = LOCAL_INDEX_TO_GLOBAL( NCURR, COMM_BETWEENOMEGA%NODE_ME, &
                 COMM_BETWEENOMEGA%NCPU, NOMEGA)
             NOMEGA_CURRENT = NCURR
             LDMODE = .FALSE.
          ELSE
             !group having less points than the root group use their largest local frequency
             NOMEGA_GLOBAL = LOCAL_INDEX_TO_GLOBAL( NOMEGA_SIMULTANEOUS, COMM_BETWEENOMEGA%NODE_ME, &
                 COMM_BETWEENOMEGA%NCPU, NOMEGA)
             NOMEGA_CURRENT = NOMEGA_SIMULTANEOUS
             LDMODE = .TRUE. 
          ENDIF
           
          IF ( LDMODE ) CYCLE

          ! isotropic average of the  dielectric matrix
          INVEPS=(CHI%HEAD(1,1,NOMEGA_CURRENT)+CHI%HEAD(2,2,NOMEGA_CURRENT)+&
             CHI%HEAD(3,3,NOMEGA_CURRENT))/3

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
             SCALE1=SQRT(ABS(CHI%HEAD(1,1,NOMEGA_CURRENT)/INVEPS))
             SCALE2=SQRT(ABS(CHI%HEAD(2,2,NOMEGA_CURRENT)/INVEPS))
             SCALE3=SQRT(ABS(CHI%HEAD(3,3,NOMEGA_CURRENT)/INVEPS))

             IF (KPOINTS_FULL%NKPX /=-1 .AND.KPOINTS_FULL%NKPY /=-1.AND.KPOINTS_FULL%NKPZ/=-1) THEN
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
                CALL vtutor%error("error in CONVERGENCE_CORRECTION_FOR_W: presently the GW routine &
                   &only supports k-point meshes generated automatically")
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

          ! set head of matrix to proper convergence corrections
          ! convergence correction times isotropic average of epsilon
          IF ( COMM_INOMEGA%NODE_ME == 1 ) THEN ! only root node!!!
          IF (LFOCKSUBTRACT ) THEN
             ! subtracted  bare Coulomb correction if Fock kernel is evaluated
             ! seperately
             IF ( CHI%LREALSTORE) THEN
                CHI%RESPONSER(1,1,NOMEGA_CURRENT)=FSG*INVEPS-FSG0
             ELSE
                CHI%RESPONSEFUN(1,1,NOMEGA_CURRENT)=FSG*INVEPS-FSG0
             ENDIF
#ifdef onlyHF
             CHI%RESPONSEFUN(1,1,NOMEGA_CURRENT)=FSG0
#endif
          ELSE
             IF ( CHI%LREALSTORE) THEN
                CHI%RESPONSEFUN(1,1,NOMEGA_CURRENT)=FSG*INVEPS
             ELSE
                CHI%RESPONSER(1,1,NOMEGA_CURRENT)=FSG*INVEPS
             ENDIF
#ifdef onlyHF
             CHI%RESPONSEFUN(1,1,NOMEGA_CURRENT)=0
#endif
          ENDIF
          ENDIF 

          ! finally take bare Coulomb off from HEAD 
          IF (LFOCKSUBTRACT) THEN
             CHI%HEAD(1,1,NOMEGA_CURRENT)=CHI%HEAD(1,1,NOMEGA_CURRENT)-1
             CHI%HEAD(2,2,NOMEGA_CURRENT)=CHI%HEAD(2,2,NOMEGA_CURRENT)-1
             CHI%HEAD(3,3,NOMEGA_CURRENT)=CHI%HEAD(3,3,NOMEGA_CURRENT)-1
          ENDIF
#ifdef onlyHF
          CHI%HEAD(1,1,NOMEGA_CURRENT)=1
          CHI%HEAD(2,2,NOMEGA_CURRENT)=1
          CHI%HEAD(3,3,NOMEGA_CURRENT)=1
#endif
          ! for gamma only the second entry is the sine transform corresponding to G=0, always zero
          ! so clear that entry
          IF (CHI%LREALSTORE) CHI%RESPONSER(2,:,NOMEGA_CURRENT)=0
          !the complete column is deleted only on the root CPU
          IF (CHI%LREALSTORE .AND. COMM_INOMEGA%NODE_ME == 1 ) CHI%RESPONSER(:,2,NOMEGA_CURRENT)=0

       !--------------------------------------------------------------------------------------
       ENDDO omega_sim
       !--------------------------------------------------------------------------------------
    ENDIF
    
 RETURN
END SUBROUTINE CONVERGENCE_CORRECTION_FOR_W

!***********************************************************************
!
! Helper routine for COMPUTE_SCREENED_POTENTIAL, multiplies CHI with 
! Coulomb potential V either from right (LTRANSPOSE) or from left 
! (LTRANSPOSE=.FALSE.) and builds epsilon=1-CHI.V, respectively
!  epsilon=1+V.CHI,
!
!***********************************************************************

  SUBROUTINE BUILD_RED_CHI_OR_EPSILON( CHI_WORK, CHI, NOMEGA_TOTAL, NROWS, NCOLS, COMM_INOMEGA, &
    WGWQ, IU0, LTRANSPOSED)
    GDEF              :: CHI_WORK(:,:)
    TYPE (responsefunction) :: CHI
    INTEGER                 :: NOMEGA_TOTAL
    INTEGER                 :: NROWS, NCOLS
    TYPE(communic)          :: COMM_INOMEGA !communicator in omega group
    TYPE (wavedes1)         :: WGWQ
    INTEGER                 :: IU0
    LOGICAL                 :: LTRANSPOSED
  ! local
    INTEGER    NI, NP, NR, NI_

    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2
 
    IF( NP > NROWS ) THEN
       CALL vtutor%error("ERROR in BUILD_RED_CHI_OR_EPSILON: NP > NROWS " // str(NP) // " " // str(NROWS))
    ENDIF 
  
    !gamma point treatment, carefull CHI is blocked
    IF (CHI%LREALSTORE) THEN
       IF (LTRANSPOSED) THEN !multiplication from right
          DO NI=1,NCOLS !columns (blocked)
              NI_=(COMM_INOMEGA%NODE_ME-1)*NCOLS + NI 
              IF ( NI_ > NROWS ) THEN
                 WRITE(*,*)'ERROR NI_(1)> NROWS for node',COMM_INOMEGA%NODE_ME,WGWQ%COMM%NODE_ME
              ENDIF
              IF ( NI_ > SIZE(WGWQ%DATAKE,1) ) THEN
                 CALL vtutor%bug("exceeding DATAKE matrix dimension " // str(COMM_INOMEGA%NODE_ME) // " " // str(WGWQ%COMM%NODE_ME) // " " // str(NI_) // " " // str(SIZE(WGWQ%DATAKE,1)), __FILE__, __LINE__)
              ENDIF

              DO NR = 1, NROWS
                 IF ( NR > NP ) CYCLE
                 CHI_WORK(NR,NI)=CHI_WORK(NR,NI)+ &
                      CHI%RESPONSER(NR, NI,NOMEGA_TOTAL)*WGWQ%DATAKE(NI_,1)
              ENDDO
          ENDDO
       ELSE !mulitplication from left
          DO NR = 1,NROWS
             IF ( NR > SIZE(WGWQ%DATAKE,1) ) THEN
                CALL vtutor%bug("exceeding DATAKE matrix dimension " // str(COMM_INOMEGA%NODE_ME) // " " // str(WGWQ%COMM%NODE_ME) // " " // str(NR) // " " // str(SIZE(WGWQ%DATAKE,1)), __FILE__, __LINE__)
             ENDIF
             
             CHI_WORK(NR,1:NCOLS)=CHI_WORK(NR,1:NCOLS)+ & 
                  CHI%RESPONSER(NR,1:NCOLS,NOMEGA_TOTAL)*WGWQ%DATAKE(NR,1)
          ENDDO
       ENDIF 
    ELSE
       IF (LTRANSPOSED) THEN !multiplication from right
          DO NI=1,NCOLS !columns (blocked)
             NI_=(COMM_INOMEGA%NODE_ME-1)*NCOLS + NI 
             IF ( NI_ > NROWS ) THEN
                WRITE(*,*)'ERROR NI_(2)> NROWS for node',COMM_INOMEGA%NODE_ME,WGWQ%COMM%NODE_ME
             ENDIF  
             DO NR = 1, NROWS
                IF ( NR > NP ) CYCLE
                CHI_WORK(NR,NI)=CHI_WORK(NR,NI)+ & 
                   CHI%RESPONSEFUN(NR,NI,NOMEGA_TOTAL)*WGWQ%DATAKE(NI_,1)
             ENDDO
          ENDDO
       ELSE
          DO NR = 1,NROWS
             CHI_WORK(NR,1:NCOLS)=CHI_WORK(NR,1:NCOLS)+ & 
                 CHI%RESPONSEFUN(NR,1:NCOLS,NOMEGA_TOTAL)*WGWQ%DATAKE(NR,1)
          ENDDO
       ENDIF
    ENDIF

    !build epsilon = 1- X.V 
    IF (LTRANSPOSED ) &
    CHI_WORK(1:NROWS,1:NCOLS)=-CHI_WORK(1:NROWS,1:NCOLS)
    !else epsilon = 1 + V.X is built
    
    DO NI=1,NCOLS !columns (blocked)
       NI_=(COMM_INOMEGA%NODE_ME-1)*NCOLS + NI 
       IF ( NI_ > NROWS ) THEN
          WRITE(*,*)'ERROR NI_(3)> NROWS for node',COMM_INOMEGA%NODE_ME,WGWQ%COMM%NODE_ME
       ENDIF  
        DO NR = 1, NROWS
           IF ( NR > NP ) CYCLE
           IF (NR==NI_) THEN
              CHI_WORK(NR,NI)=1._q+CHI_WORK(NR,NI) 
           ENDIF
        ENDDO
    ENDDO

  END SUBROUTINE BUILD_RED_CHI_OR_EPSILON
  
!***********************************************************************
!
! this routine computes the following expression used in the 
! real-frequency GW code (XI_LOCAL_FIELD)
! on entry: CHI_WORK
!                                  -1
!  X_red= X_0 (1- f_xc X_0 - v X_0) 
!
!***********************************************************************

SUBROUTINE INVERT_IRRED_CHI(CHI, CHI_WORK, NROWS, NCOLS, NP, COMM_INOMEGA, COMM_BETWEENOMEGA, &
    GLOBALCOMM, NOMEGA_CURRENT, NOMEGA_GLOBAL, NOMEGA_SIMULTANEOUS, IO ,&
    LREALC , LDMODE, LGAMMA, IDIRMAX)
    USE base
    USE constant
    USE mpimy 
    USE ini
    USE acfdt_GG
    USE scala, ONLY: FERMAT_RAZOR
    IMPLICIT NONE
    TYPE (responsefunction) :: CHI                  
    GDEF, POINTER, DIMENSION(:,:)    :: CHI_WORK    !matrix to be redistributed
    INTEGER                 :: NROWS, NCOLS         !number of local rows and cols of CHI
    INTEGER                 :: NP                   !dimension of submatrix
    TYPE (communic)         :: COMM_INOMEGA         !communicator in tau
    TYPE (communic)         :: COMM_BETWEENOMEGA    !communicator between tau
    TYPE (communic)         :: GLOBALCOMM           !global communicator
    INTEGER                 :: NOMEGA_CURRENT       !current frequency point in group
    INTEGER                 :: NOMEGA_GLOBAL        !global frequency point treated by group  
    INTEGER                 :: NOMEGA_SIMULTANEOUS  !frequencies on one node
    TYPE (in_struct)        :: IO                   !IO 
    LOGICAL                 :: LREALC               !real or complex calculation
    LOGICAL                 :: LDMODE               !dummy mode?
    LOGICAL                 :: LGAMMA               !doing longwave limit?
    INTEGER, OPTIONAL       :: IDIRMAX             !dimensions for longwave limit 
    !local                
    GDEF, POINTER           :: INV_EPS(:,:)=>NULL()  !redistributed chi
    GDEF, POINTER           :: CHI_TMP(:,:)=>NULL()     !redistributed (1-chi.v)^-1
    GDEF, POINTER           :: CHI_FINAL(:,:)=>NULL()     !redistributed chi.(1-chi.v)^-1
    INTEGER                 :: MY_NROWS, MY_NCOLS !new size of redistributed chi 
    INTEGER                 :: IERROR,INFO          !error variables
    INTEGER                 :: ICTXA                !contextes for BLACS (column major)
    INTEGER                 :: ICTXB                !contextes for BLACS (block cyclic)
    INTEGER                 :: PINFO                !scaLAPACK routine info variable
    INTEGER                 :: LWORK, LIWORK        !for scaLAPACK routine 
    INTEGER                 :: LLD_A, LLD_B         !leading (row) dimensions for desciptors  
    GDEF, ALLOCATABLE       :: WORK(:)              !working array and eigenvectors 
    INTEGER,ALLOCATABLE     :: IPIV(:),IWORK(:)     !for scaLAPACK


    IF ( LGAMMA .AND. .NOT.(PRESENT(IDIRMAX)) ) THEN
       CALL vtutor%bug("INVERT_IRRED_CHI: Longwave wanted but IDIRMAX not present", __FILE__, __LINE__)
    ENDIF

    !---------------------------------------------------------------------------------------
    !initialzie 1D column-major processor grid
    ! row major ordering (but column-major ordering is identical if there is only 1 column)
    NPROW0=1
    NPCOL0=COMM_INOMEGA%NCPU
    CALL PROCMAP( COMM_INOMEGA, ICTXA, 2, NPROW0, NPCOL0)
    CALL BLACS_GRIDINFO( ICTXA, NPROW0, NPCOL0, MYROW0, MYCOL0) !get the ids for scaLAPACK

    !define corresponding array descritpor DESCA
    !since CHI_WORK is distributed blockwisely,  
    !the leading dimension of rows is simply NROWS, 
    !check this 
    LLD_A = MAX( 1, NUMROC(NROWS, NCOLS, MYROW0, 0, NPROW0) )
    IF ( LLD_A /= NROWS ) THEN 
       CALL vtutor%bug("INVERT_IRRED_CHI LLD_A has wrong size for NODE " // str(GLOBALCOMM%NODE_ME) // " " // str(LLD_A) // " " // str(NROWS), __FILE__, __LINE__)
    ENDIF 
    CALL DESCINIT(DESCA, NROWS , NROWS, NCOLS, NCOLS, 0, 0, ICTXA, &
       LLD_A, INFO )
    !Eventually take care for errors 
    IF ( INFO /= 0 ) THEN 
       CALL vtutor%bug("Error in 1D col-major initialization " // str(INFO), __FILE__, __LINE__)
    ENDIF


    ! set NPROW,NPCOL globally
    CALL FERMAT_RAZOR( COMM_INOMEGA%NCPU, NPROW, NPCOL )
    ! generate processor grid
    CALL PROCMAP( COMM_INOMEGA, ICTXB, 1, NPROW, NPCOL)
    CALL BLACS_GRIDINFO( ICTXB, NPROW, NPCOL, MYROW, MYCOL) !get the ids for scaLAPACK

    !leading dimension of local array for block cyclic distribution
    LLD_B = MAX( 1 , NUMROC(NROWS, MB, MYROW, 0, NPROW) )

    !initialize the array descritpor
    CALL DESCINIT(DESCB, NROWS , NROWS, MB, NB, 0, 0, ICTXB, &
       LLD_B , INFO )
    IF ( INFO /= 0 ) THEN 
       CALL vtutor%error("Error in 2D block-cyclic initialization " // str(INFO))
    ENDIF

    !check if CHI_WORK has proper size
    MY_NROWS = NUMROC(DESCB(M_),DESCB(MB_),MYROW,0,NPROW)
    MY_NCOLS = NUMROC(DESCB(N_),DESCB(NB_),MYCOL,0,NPCOL) 
 
    !the following is only neccessary if we are not in dummy mode
    IF( .NOT. LDMODE ) THEN

       !allocate storage for redistributed matrix
       ALLOCATE(INV_EPS(MY_NROWS,MY_NCOLS)) ; INV_EPS = 0 

       !now call the BLACS redistribution routine 
#ifdef gammareal
       CALL PDGEMR2D( NROWS, NROWS, CHI_WORK, 1, 1, DESCA, INV_EPS, 1, 1, DESCB, ICTXB)
#else
       CALL PZGEMR2D( NROWS, NROWS, CHI_WORK, 1, 1, DESCA, INV_EPS, 1, 1, DESCB, ICTXB)
#endif

       !allocate IPIV
       ALLOCATE( IPIV(MY_NROWS+MB)) 
       
       !LU decomposition 
       IF(LREALC) THEN
          !for real values only 
          CALL PDGETRF( NP, NP, INV_EPS(1,1), 1, 1, DESCB, IPIV, PINFO)
       ELSE
          !for real values only  
          CALL PZGETRF( NP, NP, INV_EPS(1,1), 1, 1, DESCB, IPIV, PINFO)
       ENDIF 

       IF ( PINFO /=0 ) THEN
          CALL vtutor%bug("INVERT_IRRED_CHI for NOMEGA_CURRENT= " // str(NOMEGA_CURRENT) // " 1st call of P?GETRF returns " // str(PINFO), __FILE__, __LINE__)
       ENDIF


       !first get the optimal size of working arrays
       ALLOCATE(WORK(1),IWORK(1))
       IF (LREALC) THEN
          CALL PDGETRI( NP, INV_EPS(1,1), 1, 1, DESCB, IPIV(1), WORK(1), -1 ,&
             IWORK(1), -1, PINFO)
       ELSE
          CALL PZGETRI( NP, INV_EPS(1,1), 1, 1, DESCB, IPIV(1), WORK(1), -1 ,&
             IWORK(1), -1, PINFO)
       ENDIF 
   
       !allocate working arrays with optimum size 
       LWORK=WORK(1)
       LIWORK=IWORK(1)
       IF ( LWORK == 0 .OR. LIWORK ==0 ) THEN
          CALL vtutor%bug("INVERT_IRRED_CHI for NOMEGA_CURRENT= " // str(NOMEGA_CURRENT) // " LWORK or LIWORK is 0 " // str(LWORK) // " " // str(LIWORK), __FILE__, __LINE__)
       ELSE
          DEALLOCATE(WORK)
          DEALLOCATE(IWORK)
          ALLOCATE(WORK(LWORK))
          ALLOCATE(IWORK(LIWORK))
       ENDIF
       
       !invert now!
       IF (LREALC) THEN
          CALL PDGETRI( NP, INV_EPS(1,1), 1, 1, DESCB, IPIV(1), WORK(1), LWORK ,&
             IWORK(1), LIWORK, PINFO)
       ELSE
          CALL PZGETRI( NP, INV_EPS(1,1), 1, 1, DESCB, IPIV(1), WORK(1), LWORK ,&
             IWORK(1), LIWORK, PINFO)
       ENDIF 
       IF ( PINFO < 0 ) THEN
          CALL vtutor%bug("INVERT_IRRED_CHI " // str(PINFO) // " 2nd call of P?GETRI returns", __FILE__, __LINE__ )
       ENDIF

!       CALL STOP_TIMING("PD",IO%IU6,"PDSYEVD")
       !deallocate scaLAPACK's  axuillary arrays
       IF ( ALLOCATED( IPIV  )) DEALLOCATE( IPIV  )
       IF ( ALLOCATED( WORK  )) DEALLOCATE( WORK  )
       IF ( ALLOCATED( IWORK )) DEALLOCATE( IWORK ) 

       !inversion is done at this stage
       !INV_EPS contains the inverse of ( 1-X.V),
       !this needs to be multiplied by X from left, 
       !since  P?GEMM is used, auxilary arrays are required again!
       !for this purpose X needs to be distributed block wisely 
       ALLOCATE(CHI_TMP(MY_NROWS,MY_NCOLS)) ; CHI_TMP = zero
       ALLOCATE(CHI_FINAL(MY_NROWS,MY_NCOLS)) ; CHI_FINAL = zero
       !now call the BLACS redistribution routine 
#ifdef gammareal
       CALL PDGEMR2D( NROWS, NROWS, CHI%RESPONSER(:,:,NOMEGA_CURRENT), 1, 1, DESCA, &
           CHI_TMP, 1, 1, DESCB, ICTXB)
#else
       CALL PZGEMR2D( NROWS, NROWS, CHI%RESPONSEFUN(:,:,NOMEGA_CURRENT), 1, 1, DESCA, &
           CHI_TMP, 1, 1, DESCB, ICTXB)
#endif

       !call P?GEMM in 
       IF(LREALC) THEN
          CALL PDGEMM( 'n' , 'n' , NP, NP, NP, 1._q, &
             INV_EPS(1,1), 1, 1, DESCB, CHI_TMP(1,1), 1, 1, DESCB, 0._q, &
             CHI_FINAL(1,1), 1, 1, DESCB )
       ELSE
          CALL PZGEMM( 'n' , 'n' , NP, NP, NP, (1._q,0._q), &
             INV_EPS(1,1), 1, 1, DESCB, CHI_TMP(1,1), 1, 1, DESCB, (0._q,0._q), &
             CHI_FINAL(1,1), 1, 1, DESCB )
       ENDIF 
       !CHI_TMP and INV_EPS is not used in the following, so let's get rid of them 
       IF (ASSOCIATED( CHI_TMP  )) DEALLOCATE( CHI_TMP  )  
       IF (ASSOCIATED( INV_EPS )) DEALLOCATE( INV_EPS )  

       !now distribute CHI_FINAL to 1D column major 
#ifdef gammareal
       CALL PDGEMR2D( NROWS, NROWS, CHI_FINAL, 1, 1, DESCB, CHI_WORK, 1, 1, DESCA, ICTXA)
#else
       CALL PZGEMR2D( NROWS, NROWS, CHI_FINAL, 1, 1, DESCB, CHI_WORK, 1, 1, DESCA, ICTXA)
#endif  
       !and get rid of it 
       IF (ASSOCIATED( CHI_FINAL )) DEALLOCATE( CHI_FINAL )  

       !exit blacs  
       CALL BLACS_GRIDEXIT( ICTXB )  
       CALL BLACS_GRIDEXIT( ICTXA )  

    ELSE
       !exit blacs  
       CALL BLACS_GRIDEXIT( ICTXB )  
       CALL BLACS_GRIDEXIT( ICTXA )  

    ENDIF
   

    RETURN
END SUBROUTINE INVERT_IRRED_CHI


!***********************************************************************
!
! This is the analogue of SET_WING_FROM_MAT in chi_base.F
! chi is set by the matrix MAT (distributed in a 1D column major fasion) 
! distribution into frequency groups is taken into account as well. 
!
!***********************************************************************

  SUBROUTINE SET_WING_FROM_MAT_GG( MAT, CHI, COMM, NOMEGA_TOTAL, LDMODE, IDIR)
    IMPLICIT NONE
    TYPE (responsefunction) :: CHI          !full response function
    GDEF                    :: MAT(:,:)     !result 
    TYPE(communic)          :: COMM         !communicator in omega group
    INTEGER                 :: NOMEGA_TOTAL !current frequency point of group
    LOGICAL                 :: LDMODE       !dummy mode?
    INTEGER                 :: IDIR
    INTEGER                 :: NODE_SEND
    INTEGER                 :: NCOL_GLOBAL
    INTEGER                 :: NROWS, NCOLS
   

    !in the dummy mode the complete frequency group returns imediately here 
    IF ( LDMODE ) RETURN 
    !the other groups continue 
    
    NROWS=SIZE(MAT,1)
    NCOLS=SIZE(MAT,2)

    !maybe there's is only one cpu in a frequency group --> easy
    IF (COMM%NCPU == 1) THEN
       IF (CHI%LREALSTORE) THEN
          CHI%WINGR(:,IDIR,NOMEGA_TOTAL)   = MAT(:,1)
          CHI%CWINGR(:,IDIR,NOMEGA_TOTAL)  = MAT(1,:)
          CHI%HEAD(IDIR,IDIR,NOMEGA_TOTAL) = MAT(1,1)
       ELSE
          CHI%WING(:,IDIR,NOMEGA_TOTAL)    = MAT(:,1)
          CHI%CWING(:,IDIR,NOMEGA_TOTAL)   = MAT(1,:)
          CHI%HEAD(IDIR,IDIR,NOMEGA_TOTAL) = MAT(1,1)
       ENDIF
    !in the other case we need to be more careful
    ELSE    
       !for the root node the Wing and HEAD is still easy
       IF ( COMM%NODE_ME == 1 ) THEN
          IF (CHI%LREALSTORE) THEN
             CHI%WINGR(1:NROWS,IDIR,NOMEGA_TOTAL)   = MAT(1:NROWS,1)
             CHI%HEAD(IDIR,IDIR,NOMEGA_TOTAL) = MAT(1,1)
          ELSE
             CHI%WING(1:NROWS,IDIR,NOMEGA_TOTAL)    = MAT(1:NROWS,1)
             CHI%HEAD(IDIR,IDIR,NOMEGA_TOTAL) = MAT(1,1)
          ENDIF
       ELSE  !other nodes clean their head and wing
          CHI%HEAD(IDIR,IDIR,NOMEGA_TOTAL) = (0._q,0._q)
          IF (CHI%LREALSTORE) THEN 
             CHI%WINGR(1:NROWS,IDIR,NOMEGA_TOTAL) =0._q
          ELSE
             CHI%WING(1:NROWS,IDIR,NOMEGA_TOTAL) =(0._q,0._q)
          ENDIF
       ENDIF
       
       !sum in group
       CALLMPI( M_sum_z(COMM, CHI%HEAD(IDIR,IDIR,NOMEGA_TOTAL) ,1 ))
       IF (CHI%LREALSTORE) THEN 
          CALLMPI( M_sum_d(COMM, CHI%WINGR(1:NROWS,IDIR,NOMEGA_TOTAL), NROWS ))
       ELSE
          CALLMPI( M_sum_z(COMM, CHI%WING(1:NROWS,IDIR,NOMEGA_TOTAL), NROWS ))
       ENDIF    

       ! however, the conjugated wing needs to be distributed
       NCOL_GLOBAL = (COMM%NODE_ME-1)*NCOLS + 1 
       IF ( NCOL_GLOBAL > NROWS ) THEN
           WRITE(*,'(A,3I6)')'internal error in VASP: SET_WING_FROM_MAT_GG: NCOL_GLOBAL > NROWS for NODE',&
           COMM%NODE_ME,NCOL_GLOBAL,NROWS
       ENDIF 
       !NCOL_GLOBAL is the column position of the CPU in the distributed matrix
       !first each node cleans the CWING? and stores its piece of MAT to CWING?
       !then a global sum is performed
       IF (CHI%LREALSTORE) THEN
          CHI%CWINGR(1:NROWS,IDIR,NOMEGA_TOTAL)=0._q
          CHI%CWINGR(NCOL_GLOBAL:NCOL_GLOBAL+NCOLS-1,IDIR,NOMEGA_TOTAL)=MAT(1,1:NCOLS)
          CALLMPI( M_sum_d(COMM, CHI%CWINGR(1:NROWS,IDIR,NOMEGA_TOTAL), NROWS ))
       ELSE
          CHI%CWING(1:NROWS,IDIR,NOMEGA_TOTAL)=(0._q,0._q)
          CHI%CWING(NCOL_GLOBAL:NCOL_GLOBAL+NCOLS-1,IDIR,NOMEGA_TOTAL)=MAT(1,1:NCOLS)
          CALLMPI( M_sum_z(COMM, CHI%CWING(1:NROWS,IDIR,NOMEGA_TOTAL), NROWS ))
       ENDIF
    ENDIF    
  END SUBROUTINE SET_WING_FROM_MAT_GG

!***********************************************************************
!
! This is the analogue of SET_RESPONSE_FROM_MAT in chi_base.F
! chi is set by the matrix MAT (distributed in a 1D column major fasion) 
! distribution into frequency groups is taken into account as well. 
!
!***********************************************************************

  SUBROUTINE SET_RESPONSE_FROM_MAT_GG( MAT, CHI, NOMEGA_TOTAL, LDMODE)
    IMPLICIT NONE
    TYPE (responsefunction) CHI         ! full response function
    GDEF :: MAT(:,:)
    INTEGER :: NOMEGA_TOTAL
    INTEGER :: NROWS, NCOLS
    LOGICAL :: LDMODE

    !in the dummy mode we return here
    IF ( LDMODE ) RETURN
    IF (CHI%LREALSTORE) THEN
       NROWS = SIZE(CHI%RESPONSER,1)
       NCOLS = SIZE(CHI%RESPONSER,2)
    ELSE
       NROWS = SIZE(CHI%RESPONSEFUN,1)
       NCOLS = SIZE(CHI%RESPONSEFUN,2)
    ENDIF

    IF ( SIZE( MAT,1) /= NROWS .OR. SIZE(MAT,2) /= NCOLS ) THEN
       CALL vtutor%bug("SET_RESPONSE_FROM_MAT_GG: " // str(SIZE(MAT,1)) // " " // str(NROWS) // " " // str(SIZE(MAT,2)) // " " // str(NCOLS), __FILE__, __LINE__ )
    ENDIF

    IF (CHI%LREALSTORE) THEN
       CHI%RESPONSER(1:NROWS,1:NCOLS,NOMEGA_TOTAL)=MAT(1:NROWS,1:NCOLS)
    ELSE
       CHI%RESPONSEFUN(1:NROWS,1:NCOLS,NOMEGA_TOTAL)=MAT(1:NROWS,1:NCOLS)
    ENDIF
  END SUBROUTINE SET_RESPONSE_FROM_MAT_GG

!***********************************************************************
!
! helper routine, sets complex frequency array in CHI
!
!***********************************************************************
  
   SUBROUTINE SET_RESPONSE_FREQUENCIES(CHI, IMAG_GRIDS)
      USE minimax_struct, ONLY : imag_grid_handle
      IMPLICIT NONE
      TYPE (responsefunction) :: CHI
      TYPE (imag_grid_handle) :: IMAG_GRIDS 
      INTEGER                 :: I

      IF ( .NOT. ASSOCIATED(CHI%COMEGA) ) THEN
         ALLOCATE(CHI%COMEGA(IMAG_GRIDS%B%NPOINTS_IN_GROUP))
      ENDIF
 
      IF( SIZE(CHI%COMEGA) < IMAG_GRIDS%B%NPOINTS_IN_GROUP ) THEN
         CALL vtutor%bug("SET_RESPONSE_FREQUENCIES, CHI%COMEGA too small "// str(SIZE(CHI%COMEGA)) // " "// str(IMAG_GRIDS%B%NPOINTS_IN_GROUP), __FILE__, __LINE__)
      ENDIF
      IF( SIZE(IMAG_GRIDS%B%POINTS_LOCAL) < IMAG_GRIDS%B%NPOINTS_IN_GROUP ) THEN
         CALL vtutor%bug("SET_RESPONSE_FREQUENCIES, B%POINTS_LOCAL too small " // str(SIZE(IMAG_GRIDS%B%POINTS_LOCAL)) // " " // str(IMAG_GRIDS%B%NPOINTS_IN_GROUP), __FILE__, __LINE__ )
      ENDIF
      DO I = 1, IMAG_GRIDS%B%NPOINTS_IN_GROUP
        CHI%COMEGA(I)=CMPLX(0._q,IMAG_GRIDS%B%POINTS_LOCAL(I),q)
      ENDDO
   
   END SUBROUTINE SET_RESPONSE_FREQUENCIES

!***********************************************************************
!
! helper routine, dumps response function 
!
!***********************************************************************
  
  SUBROUTINE DUMP_HEAD_GG( CHI, IMAG_GRIDS, STRING, &
    IU, NWRITE, LINVERS, FAKT, ADD)
    USE minimax_struct, ONLY : imag_grid_handle
    IMPLICIT NONE 
    TYPE (responsefunction) :: CHI                 !response function 
    TYPE (imag_grid_handle) :: IMAG_GRIDS          !imaginary grids
    CHARACTER (LEN=*)       :: STRING              !string to be written
    INTEGER, INTENT(IN)     :: IU                  !writing unit
    INTEGER, INTENT(IN)     :: NWRITE              !IO%NWRITE
    LOGICAL, OPTIONAL       :: LINVERS             !write inverse of head
    REAL(q), OPTIONAL       :: FAKT                !scale head with FAKT
    real(q), OPTIONAL       :: ADD                 !add term ADD to head

  ! local 
    INTEGER                 :: NOMEGA_SIMULTANEOUS !# of frequencies in one group  
    TYPE(communic)          :: COMMBET             !communicator between groups
    INTEGER    :: N, I, J, NPOS, NT
    COMPLEX(q), ALLOCATABLE :: HEAD(:,:,:)
    COMPLEX(q), ALLOCATABLE :: COMEGA(:)
    REAL(q), ALLOCATABLE    :: DOMEGA(:)
    REAL(q), ALLOCATABLE    :: EPS_REAL(:,:,:) 
    REAL(q), ALLOCATABLE    :: EPS_IMAG(:,:,:) 
    INTEGER    :: INDEX
 
    NOMEGA_SIMULTANEOUS = IMAG_GRIDS%B%NPOINTS_IN_GROUP
    COMMBET = IMAG_GRIDS%B%COMM_BETWEEN_GROUPS

    ALLOCATE( HEAD( 3, 3, IMAG_GRIDS%B%NPOINTS ) )
    HEAD = ( 0._q , 0._q )
    ALLOCATE( COMEGA( IMAG_GRIDS%B%NPOINTS ) )
    COMEGA = ( 0._q , 0._q )
    ALLOCATE( DOMEGA( IMAG_GRIDS%B%NPOINTS ) )
    DOMEGA = 0._q
    ALLOCATE( EPS_REAL( 3, 3, IMAG_GRIDS%B%NPOINTS ) )
    EPS_REAL = 0._q
    ALLOCATE( EPS_IMAG( 3, 3, IMAG_GRIDS%B%NPOINTS ) )
    EPS_IMAG = 0._q

    !if we are not calling this routine at the Gamma point return imediately
    !and degree of information flag (IO%NWRITE) is below 3  
    IF (.NOT.CHI%LGAMMA .AND. NWRITE <3) THEN
       DEALLOCATE( EPS_REAL, EPS_IMAG )
       DEALLOCATE( COMEGA, DOMEGA )
       DEALLOCATE( HEAD )
       RETURN
    ENDIF

    IF ( .NOT.ASSOCIATED(CHI%COMEGA)) THEN
       CALL vtutor%bug("DUMP_HEAD_GG, CHI%COMEGA not allocated", __FILE__, __LINE__)
    ENDIF
    IF ( SIZE(CHI%COMEGA) < NOMEGA_SIMULTANEOUS) THEN
       CALL vtutor%bug("DUMP_HEAD_GG, CHI%COMEGA too small", __FILE__, __LINE__ )
    ENDIF 
 
    HEAD=(0._q,0._q)
    COMEGA=(0._q,0._q)
    !we need to collect heads from all frequency groups first
    !only root group needs this information
    IF ( COMMBET%NODE_ME > 1 ) THEN
       NPOS=SUM(IMAG_GRIDS%B%DISTRIBUTION(1:COMMBET%NODE_ME-1,2))+1
    ELSE
       NPOS = 1 
    ENDIF
!WRITE(100+COMM%NODE_ME,'(10I3)')COMMBET%NODE_ME,NPOS,IMAG_GRIDS%B%DISTRIBUTION(COMMBET%NODE_ME,:)
    !save head of each group to HEAD
    DO N=1, NOMEGA_SIMULTANEOUS
       NT = LOCAL_INDEX_TO_GLOBAL(N,COMMBET%NODE_ME,COMMBET%NCPU,IMAG_GRIDS%B%NPOINTS)
       IF ( NT > IMAG_GRIDS%B%NPOINTS ) THEN 
          WRITE(*,*)'internal error in VASP: DUMP_HEAD_GG NT out of bounds', NT, IMAG_GRIDS%B%NPOINTS
       ENDIF 
       DO I = 1, 3
          DO J = 1, 3  
             HEAD(I,J,NT)=CHI%HEAD(I,J,N)
          ENDDO
       ENDDO
       COMEGA(NT)=CHI%COMEGA(N)
    ENDDO
    !communicate the stuff 
    CALLMPI(M_sum_z(COMMBET, HEAD, 3*3*IMAG_GRIDS%B%NPOINTS))
    CALLMPI(M_sum_z(COMMBET, COMEGA , IMAG_GRIDS%B%NPOINTS))

    IF (PRESENT(FAKT) )THEN
       HEAD=HEAD*FAKT
    ENDIF
    IF ( PRESENT(ADD)) THEN
       HEAD(1,1,:)=HEAD(1,1,:)+ADD
       HEAD(2,2,:)=HEAD(2,2,:)+ADD
       HEAD(3,3,:)=HEAD(3,3,:)+ADD
    ENDIF

    IF (IU>=0) THEN
       WRITE(IU,1100) STRING
       INDEX=0
       DO N=1,IMAG_GRIDS%B%NPOINTS
          WRITE(IU,'(" w=",2F10.3)')  COMEGA(N)
          IF (CHI%LGAMMA) THEN
             IF (PRESENT(LINVERS)) THEN
                IF (LINVERS) THEN
                HEAD(1,1,N)=1._q/HEAD(1,1,N)
                HEAD(2,2,N)=1._q/HEAD(2,2,N)
                HEAD(3,3,N)=1._q/HEAD(3,3,N)
                ENDIF
             ENDIF

             WRITE(IU,'(3(5X,3(2F10.4,2X)/))')  HEAD(:,:,N)
             WRITE(IU,'(3X,F10.3,4X,2F10.3,A)') ABS(COMEGA(N)), &
                (HEAD(1,1,N)+ HEAD(2,2,N)+ HEAD(3,3,N))/3._q, " dielectric  constant"
             INDEX=INDEX+1
             DOMEGA(INDEX)=ABS(COMEGA(N))
             EPS_REAL(:,:,INDEX)=REAL(HEAD(:,:,N),q)
             EPS_IMAG(:,:,INDEX)=AIMAG(HEAD(:,:,N))
          ENDIF
          WRITE(IU,*)
       ENDDO
       IF (CHI%LGAMMA) THEN
          CALL XML_EPSILON_E(DOMEGA, EPS_REAL, EPS_IMAG, INDEX, STRING )
       ENDIF
    ENDIF

    ! get rid of local arrays 
    IF ( ALLOCATED( EPS_REAL ) ) DEALLOCATE( EPS_REAL )
    IF ( ALLOCATED( EPS_IMAG ) ) DEALLOCATE( EPS_IMAG )
    IF ( ALLOCATED( COMEGA ) ) DEALLOCATE( COMEGA )
    IF ( ALLOCATED( DOMEGA ) ) DEALLOCATE( DOMEGA )
    IF ( ALLOCATED( HEAD) ) DEALLOCATE( HEAD )

1100 FORMAT(/" ",A/, &
             " -------------------------------------")

  END SUBROUTINE DUMP_HEAD_GG


!**********************************************************************
!
! writes the diagonal part of the 'blocked' response function to 
! vasprun.xml
! TODO: rewrite for blocked response functions
! 
!**********************************************************************

  SUBROUTINE DUMP_CHI_GG_XML( CHI, WDESQ, LATT_CUR, STRING)
    USE constant
    USE mgrid
    USE lattice
    USE wave 

    IMPLICIT NONE
    TYPE (responsefunction) CHI
    TYPE (wavedes1) WDESQ
    TYPE (latt) LATT_CUR
    ! local
    INTEGER    NP, NI, NI_, NG
    REAL(q) :: DKX, DKY, DKZ, GX, GY, GZ, GSQU, SCALE
    CHARACTER (LEN=*) :: STRING

    REAL(q), ALLOCATABLE :: A(:,:)

    IF (CHI%NOMEGA_LOW/=1) RETURN

    SCALE=TPI**2

    DKX=(CHI%VKPT(1))*LATT_CUR%B(1,1)+ &
        (CHI%VKPT(2))*LATT_CUR%B(1,2)+ &
        (CHI%VKPT(3))*LATT_CUR%B(1,3)
    DKY=(CHI%VKPT(1))*LATT_CUR%B(2,1)+ &
        (CHI%VKPT(2))*LATT_CUR%B(2,2)+ &
        (CHI%VKPT(3))*LATT_CUR%B(2,3)
    DKZ=(CHI%VKPT(1))*LATT_CUR%B(3,1)+ &
        (CHI%VKPT(2))*LATT_CUR%B(3,2)+ &
        (CHI%VKPT(3))*LATT_CUR%B(3,3)

    NP=WDESQ%NGVECTOR
    IF (WDESQ%LGAMMA) NP=NP*2

    ALLOCATE(A(2,NP))
    A=0

    DO NI=1,WDESQ%NGVECTOR
       NI_=NI
       IF (WDESQ%LGAMMA) NI_=(NI-1)/2+1

       GX=(WDESQ%IGX(NI_)*LATT_CUR%B(1,1)+WDESQ%IGY(NI_)* &
            LATT_CUR%B(1,2)+WDESQ%IGZ(NI_)*LATT_CUR%B(1,3))
       GY=(WDESQ%IGX(NI_)*LATT_CUR%B(2,1)+WDESQ%IGY(NI_)* &
            LATT_CUR%B(2,2)+WDESQ%IGZ(NI_)*LATT_CUR%B(2,3))
       GZ=(WDESQ%IGX(NI_)*LATT_CUR%B(3,1)+WDESQ%IGY(NI_)* &
            LATT_CUR%B(3,2)+WDESQ%IGZ(NI_)*LATT_CUR%B(3,3))
          
       GSQU=(DKX+GX)**2+(DKY+GY)**2+(DKZ+GZ)**2

       A(1,NI)=SQRT(GSQU*SCALE)
       IF (CHI%LREALSTORE) THEN
          A(2,NI)=CHI%RESPONSER  (NI,NI,1)
       ELSE
          A(2,NI)=CHI%RESPONSEFUN(NI,NI,1)
       ENDIF
    ENDDO

    IF (CHI%LGAMMA) THEN
       A(2,1)=(CHI%HEAD(1,1,1)+CHI%HEAD(2,2,1)+CHI%HEAD(3,3,1))/3
    ENDIF
      
    CALL XML_VECARRAY(STRING)
    CALL XML_ARRAY_REAL(A(:,1:NP))
    CALL XML_CLOSE_TAG

    DEALLOCATE(A)

  END SUBROUTINE DUMP_CHI_GG_XML

!***********************************************************************
! determines global column LOCAL_TO_GLOBAL_COL from the local point
! NLOCAL, the group id ID and the blocking factor NBLOCK
! blocked distribution is assumed
!***********************************************************************

   FUNCTION LOCAL_TO_GLOBAL_COL( NLOCAL, NBLOCK, ID,NODEME, NTOT )
      IMPLICIT NONE
      INTEGER :: LOCAL_TO_GLOBAL_COL !global point
      INTEGER :: NLOCAL              !local point in group
      INTEGER :: NBLOCK              !blocking factor 
      INTEGER :: ID                  !number of group
      INTEGER :: NODEME              !nodeme
      INTEGER, OPTIONAL :: NTOT      !size of global array

      !the global index is
      LOCAL_TO_GLOBAL_COL = (ID-1)*NBLOCK + NLOCAL

      IF ( PRESENT(NTOT) ) THEN
         IF ( LOCAL_TO_GLOBAL_COL > NTOT ) THEN
            CALL vtutor%bug("LOCAL_TO_GLOBAL_COL global index is out bounds " // str(LOCAL_TO_GLOBAL_COL) // " " // str(NTOT) // " " // str(ID) // " " // str(NODEME),  __FILE__, __LINE__ )
         ENDIF
      ENDIF
      RETURN
   END FUNCTION LOCAL_TO_GLOBAL_COL

!***********************************************************************
!
! sets up send and receive arrays for alltoallv call in cos trafo.
!
!***********************************************************************

   SUBROUTINE SETUP_SEND_RECV_ARRAYS( GLOBALCOMM, COMM_INSEND, COMM_INRECV,ICOL_RECV,NROWS,&
      NSENDPAR,NRECVPAR,ICOL_SEND,SEND_COUNTS,SEND_DISPLS,RECV_COUNTS,RECV_DISPLS,LFORWARD)
      IMPLICIT NONE
      TYPE(communic) :: GLOBALCOMM    !global communicator
      TYPE(communic) :: COMM_INSEND   !communicator in sender group 
      TYPE(communic) :: COMM_INRECV   !communicator in receiver group
      INTEGER        :: ICOL_RECV     !current local column index of receiver
      INTEGER        :: NROWS         !total number of columns
      INTEGER        :: NSENDPAR      !number of sender groups 
      INTEGER        :: NRECVPAR      !number of receiver groups
      INTEGER        :: ICOL_SEND     !local column index of sender
      INTEGER        :: SEND_COUNTS(:)!entries sender sends
      INTEGER        :: SEND_DISPLS(:)!displacement of sender
      INTEGER        :: RECV_COUNTS(:)!bits sender receives
      INTEGER        :: RECV_DISPLS(:)!displacement of receiver
      LOGICAL        :: LFORWARD
      !local variables
      INTEGER        :: NBLOCK_RECV !blocking factor of reciver 
      INTEGER        :: NBLOCK_SEND !blocking factor of sender 
      INTEGER        :: ID_SEND
      INTEGER        :: ID_RECV
      INTEGER        :: I,J,K,L
      INTEGER        :: NPOS_SENDER !lowest global column of sender 
      INTEGER        :: IRECV_COL(GLOBALCOMM%NCPU)  !global column index of receiver
      INTEGER        :: ISEND_COL(GLOBALCOMM%NCPU)  !local column index of sender
      INTEGER        :: ISEND_COL_RED(COMM_INRECV%NCPU )!red. part of loc. col. array of sender
      INTEGER        :: ARED(COMM_INRECV%NCPU,COMM_INSEND%NCPU) !irreducible adjacency matrix 
      INTEGER        :: A(GLOBALCOMM%NCPU,GLOBALCOMM%NCPU)   !full adjacency matrix 

      !determine blocking factor of receiver
      NBLOCK_RECV=NROWS/COMM_INRECV%NCPU
      !clean receiver array
      IRECV_COL=0
      !find corresponding global index of ICOL_RECV
      IRECV_COL(GLOBALCOMM%NODE_ME)=LOCAL_TO_GLOBAL_COL(ICOL_RECV,NBLOCK_RECV,&
         COMM_INRECV%NODE_ME,GLOBALCOMM%NODE_ME,NROWS)
      !communicate to others 
      CALLMPI(M_sum_i(GLOBALCOMM,IRECV_COL,GLOBALCOMM%NCPU))

      !determine blocking factor of sender
      NBLOCK_SEND=NROWS/COMM_INSEND%NCPU

      !initialize sending array
      ISEND_COL_RED=0
      !initialize adjacency matrix  
      ARED=0
      A=0

      !find out who is sending these global columns
      DO ID_SEND = 1, COMM_INSEND%NCPU
         !global starting position of ID_SEND is
         NPOS_SENDER=(ID_SEND-1)*NBLOCK_SEND+1
         !loop over all receiver
         DO ID_RECV = 1, COMM_INRECV%NCPU
            !check if the column is carried by node ID_SEND
            IF ( NPOS_SENDER <= IRECV_COL(ID_RECV) .AND. &
                 IRECV_COL(ID_RECV)<= NPOS_SENDER+NBLOCK_SEND-1 ) THEN
               !if this is the case, convert global column index to local index of sender 
               IF ( MOD(IRECV_COL(ID_RECV),NBLOCK_SEND)==0 ) THEN
                  !save this information to adjacency matrix
                  ARED(ID_RECV,ID_SEND)=NBLOCK_SEND
               ELSE
                  ARED(ID_RECV,ID_SEND)=MOD(IRECV_COL(ID_RECV),NBLOCK_SEND)
               ENDIF
            ENDIF
         ENDDO
         !if the sender has the id of ID_SEND than save this sending index
         IF ( ID_SEND == COMM_INSEND%NODE_ME) THEN
            ISEND_COL_RED(1:COMM_INRECV%NCPU)=ARED(1:COMM_INRECV%NCPU,ID_SEND)
         ENDIF
      ENDDO

      !at this stage ARED is set up in the irreducible subset of the communication graph
      !from this we construct the global communication graph 
      ISEND_COL=0

#ifdef debug
      WRITE(800+ GLOBALCOMM%NODE_ME,'(48I3)')  GLOBALCOMM%NODE_ME, ARED
#endif

      DO ID_RECV = 1, NRECVPAR
         DO ID_SEND = 1, NSENDPAR
            A(1+(ID_RECV-1)*COMM_INRECV%NCPU:ID_RECV*COMM_INRECV%NCPU,&
              1+(ID_SEND-1)*COMM_INSEND%NCPU:ID_SEND*COMM_INSEND%NCPU)=&
              ARED(1:COMM_INRECV%NCPU,1:COMM_INSEND%NCPU)
         ENDDO
         !this is the local column the sender sends
         ISEND_COL(1+(ID_RECV-1)*COMM_INRECV%NCPU:ID_RECV*COMM_INRECV%NCPU)=&
            ISEND_COL_RED(1:COMM_INRECV%NCPU)
      ENDDO

#ifdef debug
      DO ID_RECV = 1, GLOBALCOMM%NCPU
         WRITE(900+ GLOBALCOMM%NODE_ME,'(48I3)')A(ID_RECV,1:GLOBALCOMM%NCPU)
      ENDDO
#endif

      !setup send and receive arrays
      SEND_COUNTS=0
      SEND_DISPLS=0
      RECV_COUNTS=0
      RECV_DISPLS=0
      I=0
      J=0
      DO ID_SEND = 1, GLOBALCOMM%NCPU
         IF ( A(ID_SEND,GLOBALCOMM%NODE_ME)>0) THEN
            SEND_COUNTS(ID_SEND)=NROWS
            I=I+1
            !zero based displacement array (maybe more than one column needs to be sent)
            SEND_DISPLS(ID_SEND)=(A(ID_SEND,GLOBALCOMM%NODE_ME)-1)*NROWS
!            SEND_DISPLS(ID_SEND)=A(ID_SEND,GLOBALCOMM%NODE_ME)
         ENDIF
         !determine recv counts
         IF ( A(GLOBALCOMM%NODE_ME,ID_SEND) > 0  ) THEN
            J=J+1
            RECV_COUNTS(ID_SEND)=NROWS
            !zero based
            RECV_DISPLS(ID_SEND)=(J-1)*NROWS
         ENDIF
      ENDDO

#ifdef debug
      !last but not least determine local column index of sender
      !initialize to 1
      ICOL_SEND=1
      DO ID_SEND = 1, GLOBALCOMM%NODE_ME
         IF( ISEND_COL(ID_SEND)>0 ) THEN
            ICOL_SEND=ISEND_COL(ID_SEND)
            EXIT
         ENDIF
      ENDDO
#endif 
      !the adjacency matrix is now set up 
      !it describes what local columns the sender sends to the receiver
   END SUBROUTINE SETUP_SEND_RECV_ARRAYS

!***********************************************************************
!
! dump progress bar 
!
!***********************************************************************
   SUBROUTINE DUMP_PROGRESS( ITER, ITERMAX, TYP_, IO ) 
      use base
      INTEGER                   :: ITER        ! current position
      INTEGER                   :: ITERMAX     ! max iterations
      CHARACTER(LEN=*)          :: TYP_        ! what is printed
      TYPE(in_struct)           :: IO          ! who is printing
      ! local 
      LOGICAL,SAVE              :: ISTARTED=.FALSE.
      INTEGER,SAVE              :: IOLD=1
      INTEGER, PARAMETER        :: IMAX = 40
      INTEGER                   :: I,K
      REAL                      :: IP
      CHARACTER(LEN=5)          :: TYP        ! what is printed
      CHARACTER(LEN=1)          :: S
      PARAMETER( S='-')                       ! progress bar symbol
      !
      ! print only for NWRITE > 2
      !
      IF ( IO%NWRITE < 5 ) RETURN     
      !
      ! cut string to proper length
      !
      TYP=TYP_
      !
      ! progress in percent 
      !
      IP = REAL(ITER)/REAL(ITERMAX)
      !
      ! progress inside bar
      !
      I = MIN( INT(IP*IMAX), IMAX )
      !
      !print progress bar
      !
      IF ( IO%IU0 >=0 ) THEN
         !
         ! initialize progress bar
         !
         IF ( .NOT. ISTARTED )THEN
            WRITE(IO%IU0,FMT='(" ",A," :|")',ADVANCE='NO')TYP
            ISTARTED = .TRUE.
         ENDIF
         DO K = IOLD, I 
            WRITE(IO%IU0,FMT='(A)',ADVANCE='NO')S
         ENDDO
         IOLD = K
         !
         ! finish progress bar 
         !
         IF ( ITER==ITERMAX ) THEN
            IF ( IOLD < IMAX) THEN
               DO K=IOLD, IMAX
                  WRITE(IO%IU0,FMT='(A)',ADVANCE='NO')S
               ENDDO
            ENDIF
            ISTARTED = .FALSE.
            IOLD=1
            WRITE(IO%IU0,FMT='(A)')'|'
         ENDIF
      ENDIF
   END SUBROUTINE DUMP_PROGRESS

!***********************************************************************
!
! storage requirements for GW and RPA/ACFDT calculations 
!
!***********************************************************************
   SUBROUTINE STORAGE_REQ_GG(TYP, MEM_GG, GDES_TAU, GDES, T, F, S2E, WDES, WDES_RESPONSE, S )
      USE ini, ONLY: QUERRY_ALLOCATE
      USE main_mpi, ONLY: COMM_WORLD
      USE tutor, ONLY: vtutor
      TYPE( mem_gw_handle )      :: MEM_GG 
      TYPE (greensfdes)          :: GDES_TAU             ! green function in time domain  
      TYPE (greensfdes)          :: GDES                 ! green function in freq domain  
      TYPE(loop_des)             :: T, F
      TYPE( screened_2e_handle ) :: S2E
      TYPE( wavedes )            :: WDES
      TYPE( wavedes )            :: WDES_RESPONSE
      TYPE( supercell ), POINTER :: S
      CHARACTER(LEN=1)           :: TYP
      ! local 
      LOGICAL                    :: LOMEGA 
      REAL(q), PARAMETER         :: PERCENT=1.06_q
      REAL(q)                    :: PSS 

      LOMEGA = .FALSE.
      IF ( TYP == "O" .OR. TYP == "o" ) LOMEGA = .TRUE. 

      ! initialize wavefunction 
      ! account for additional wavefunction (response wavefunction )
      MEM_GG%M_WAVEFUN = 2*QUERRY_ALLOCATE()/1024._q*PERCENT
      CALL GET_PSS( PSS, -1, -1, " ", 0 )
      ! safety if virtual memory resisident set size not found
      IF ( PSS > 0 ) THEN
         MEM_GG%M_REST =  PSS
      ELSE
         MEM_GG%M_REST =  MEM_GG%M_WAVEFUN 
      ENDIF
      ! contribution from chi matrices in frequency domain
      ! and reciprocal space 
      CALL STORAGE_FOR_CHI_OMEGA( MEM_GG%M_CHIOMEGA, GDES, F, S2E%NUMBER_OF_NQ )
      CALLMPI( M_max_d( COMM_WORLD, MEM_GG%M_CHIOMEGA, 1) )
      IF ( LOMEGA ) RETURN 
      ! 
      ! storage requirement for one Green's function at one time point:
      ! corresponds to GU and GO in CALCULATE_G_POSSIBLY_SC
      ! 
      CALL STORAGE_FOR_G_TAU( MEM_GG%M_GTAU, GDES_TAU, WDES%NKPTS, T )
      CALLMPI( M_max_d( COMM_WORLD, MEM_GG%M_GTAU, 1) )

      ! contribution from chi matrices in time domain 
      CALL STORAGE_FOR_CHI_TAU( MEM_GG%M_CHITAU, GDES_TAU, S2E, T )
      CALLMPI( M_max_d( COMM_WORLD, MEM_GG%M_CHITAU, 1) )

      ! storage requirement for actual contraction routines 
      ! first consider only GG 
      CALL STORAGE_FOR_GG_TAU( MEM_GG%M_GTAU_G, GDES_TAU, WDES, S2E, T, S ) 
      CALLMPI( M_max_d( COMM_WORLD, MEM_GG%M_GTAU_G, 1) )

      ! storage for self-energy or response function
      ! second consider GW
      CALL STORAGE_FOR_GTAU_WTAU( MEM_GG%M_GTAU_WTAU, GDES_TAU, WDES, S2E, T, S, WDES_RESPONSE ) 
      CALLMPI( M_max_d( COMM_WORLD, MEM_GG%M_GTAU_WTAU, 1) )

      ! storage requirement for orbital basis 
      CALL STORAGE_FOR_SIGMA( MEM_GG%M_ORB, GDES_TAU, T, F, S2E, WDES ) 
      CALLMPI( M_max_d( COMM_WORLD, MEM_GG%M_ORB, 1) )

      ! add all contributions 
      MEM_GG%M_TOTAL = MEM_GG%M_GTAU + MAX(MEM_GG%M_GTAU_G , MEM_GG%M_GTAU_WTAU)  + &
                       MEM_GG%M_CHITAU + MEM_GG%M_CHIOMEGA + MEM_GG%M_ORB  + &
                       MEM_GG%M_WAVEFUN + MEM_GG%M_REST

#ifdef debug
      IF ( vtutor%unitOut  >= 0 ) THEN
        WRITE( vtutor%unitout, * ) "required memory composition:"
        WRITE( vtutor%unitout, '(A, F12.4, F14.4 )') 'M_TOTAL     ', MEM_GG%M_TOTAL
        WRITE( vtutor%unitout, '(A, F12.4, F14.4 )') 'M_GTAU      ', MEM_GG%M_GTAU
        WRITE( vtutor%unitout, '(A, F12.4, F14.4 )') 'M_GTAU_G    ', MEM_GG%M_GTAU_G
        WRITE( vtutor%unitout, '(A, F12.4, F14.4 )') 'M_GTAU_WTAU ', MEM_GG%M_GTAU_WTAU
        WRITE( vtutor%unitout, '(A, F12.4, F14.4 )') 'M_CHITAU    ', MEM_GG%M_CHITAU
        WRITE( vtutor%unitout, '(A, F12.4, F14.4 )') 'M_CHIOMEGA  ', MEM_GG%M_CHIOMEGA
        WRITE( vtutor%unitout, '(A, F12.4, F14.4 )') 'M_ORB       ', MEM_GG%M_ORB
        WRITE( vtutor%unitout, '(A, F12.4, F14.4 )') 'M_WAVEFUN   ', MEM_GG%M_WAVEFUN
        WRITE( vtutor%unitout, '(A, F12.4, F14.4 )') 'M_REST      ', MEM_GG%M_REST
        WRITE( vtutor%unitout, '(A, F12.4, F14.4 )') 'PSS         ', REAL(PSS,q)
      ENDIF
#endif

      CONTAINS 
     
      !************************************************************************
      ! memory requirement of single Green's function at one time point
      ! deduced from ALLOCATE_G_RECIPROCAL for G%GG, G%G_PROJ G%PROJ_PROJ
      !************************************************************************
      SUBROUTINE STORAGE_FOR_G_TAU( M_GTAU, GDES_TAU, NKPTS, T )
         REAL(q), INTENT(INOUT) :: M_GTAU
         TYPE (greensfdes)      :: GDES_TAU             ! green function in time domain  
         INTEGER, INTENT(IN)    :: NKPTS
         TYPE(loop_des)         :: T

         M_GTAU     = 0 
         IF ( ICHIREAL >= 1 ) THEN
            M_GTAU = M_GTAU                                                      &
               +(GDES_TAU%NRPLWV_ROW_DATA_POINTS*GDES_TAU%NRPLWV_COL_DATA_POINTS)&
               +(GDES_TAU%NRPLWV_ROW_DATA_POINTS*GDES_TAU%NPRO_COL)              &
               +(GDES_TAU%NPRO_ROW*GDES_TAU%NPRO_COL) 
            ! scale by number of k-points
            M_GTAU = M_GTAU * NKPTS
         ENDIF
        

         ! two greens functions are required
         ! convert to MB 
         M_GTAU = 2*M_GTAU * ndata *8._q / 1024/ 1024*PERCENT
                                                      

      END SUBROUTINE STORAGE_FOR_G_TAU

      !************************************************************************
      ! memory requirement for all response functions in time domain
      ! ALLOCATE_RESPONSEFUN_DISTRI
      !************************************************************************
      SUBROUTINE STORAGE_FOR_CHI_TAU( M_CHITAU, GDES_TAU, S2E, T )
         REAL(q), INTENT(INOUT) :: M_CHITAU
         TYPE (greensfdes)      :: GDES_TAU           
         TYPE( screened_2e_handle ) :: S2E
         TYPE(loop_des)         :: T

         M_CHITAU = 0 

         IF ( ICHIREAL == 1 ) THEN
            M_CHITAU = M_CHITAU  &
               +(GDES_TAU%RES_NRPLWV_ROW*GDES_TAU%RES_NRPLWV_COL) 
            ! deduced in CALCULATE_CHI_TAU_FROM_G for 
            ! P_G_PROJ, P_G_R, P_PROJ_G, P_R_G
            M_CHITAU = M_CHITAU                                    &
               +(GDES_TAU%RES_NRPLWV_ROW_DATA_POINTS*GDES_TAU%NLM_COL)           &
               +(GDES_TAU%RES_NRPLWV_ROW_DATA_POINTS*GDES_TAU%MPLWV_COL)         &
               +(GDES_TAU%NLM_ROW*GDES_TAU%RES_NRPLWV_COL_DATA_POINTS)           &
               +(GDES_TAU%MPLWV_ROW*GDES_TAU%RES_NRPLWV_COL_DATA_POINTS)  

         ELSE IF ( ICHIREAL == 2 ) THEN
            ! storage requirement for one Green's function at one time point:
            ! and all q-points 
            M_CHITAU = M_CHITAU  &
               +(GDES_TAU%RES_NRPLWV_ROW*GDES_TAU%RES_NRPLWV_COL*WDES%NKPTS)  

            IF ( LOEP ) THEN
               M_CHITAU = M_CHITAU  &
                  +(GDES_TAU%RES_NRPLWV_ROW*GDES_TAU%RES_NRPLWV_COL)
            ENDIF
         ENDIF
         M_CHITAU = M_CHITAU * ndata *8._q / 1024/ 1024 *PERCENT
      END SUBROUTINE STORAGE_FOR_CHI_TAU

      !************************************************************************
      ! memory requirement for all response functions in frequency domain
      ! ALLOCATE_RESPONSE_DISTRI
      !************************************************************************
      SUBROUTINE STORAGE_FOR_CHI_OMEGA( M_CHIOMEGA, GDES, F, NKPTS )
         REAL(q), INTENT(INOUT) :: M_CHIOMEGA
         TYPE (greensfdes)      :: GDES           
         TYPE(loop_des)         :: F
         INTEGER, INTENT(IN)    :: NKPTS 

         M_CHIOMEGA = 0   ! treated as peak memory usage
   
         IF ( ICHIREAL >= 1 ) THEN
            M_CHIOMEGA = M_CHIOMEGA &
                       + NKPTS*GDES%RES_NRPLWV_ROW_DATA_POINTS*GDES%RES_NRPLWV_COL_DATA_POINTS &
                       * F%NPOINTS_IN_ROOT_GROUP
            ! computation of W 
            IF ( .NOT. ( LACFDT .AND. .NOT. LRPAFORCE ) )  THEN
               M_CHIOMEGA = M_CHIOMEGA  &
                    + GDES%RES_NRPLWV_ROW_DATA_POINTS*GDES%RES_NRPLWV_COL_DATA_POINTS
            ENDIF 
         ENDIF

         ! convert to MB 
         M_CHIOMEGA = M_CHIOMEGA * ndata *8._q /1024/ 1024*PERCENT 
         ! double memory for gamma-only calculations
         IF ( WDES%NKPTS == 1 ) THEN
            M_CHIOMEGA = M_CHIOMEGA * 2 
         ENDIF
      END SUBROUTINE STORAGE_FOR_CHI_OMEGA

      !************************************************************************
      ! memory requirement for contraction of two greens functions in time 
      ! and real space and fourier transformation to reciprocal space 
      !************************************************************************

      SUBROUTINE STORAGE_FOR_GG_TAU( M_GTAU_G, GDES_TAU, WDES, S2E, T, S )
         REAL(q), INTENT(INOUT) :: M_GTAU_G
         TYPE (greensfdes)      :: GDES_TAU
         TYPE (wavedes)         :: WDES
         TYPE( screened_2e_handle ) :: S2E
         TYPE(loop_des)             :: T
         TYPE( supercell ), POINTER :: S
         INTEGER                    :: NRPLWV_COL_DATA_POINTS
         INTEGER                    :: NRPLWV_ROW_DATA_POINTS
         REAL(q)                    :: M_TOTAL 
         INTEGER                    :: NQ, NK 
         ! 
         M_GTAU_G = 0    ! treated as peak memory usage 
         M_TOTAL = 0     ! current memory usage 

         IF ( ICHIREAL == 1 ) THEN
            DO NQ=1,S2E%NUMBER_OF_NQ  ! loop over  q-points in the IRZ
            DO NK=1,WDES%NKPTS
               CALL MEM_RESPONSE_GG( M_GTAU_G, M_TOTAL, GDES_TAU)
            ENDDO
            ENDDO

         ELSE IF ( ICHIREAL == 2 ) THEN
            CALL MEM_RESPONSE_SUPER( M_GTAU_G, M_TOTAL, GDES_TAU, S, &
               WDES, S2E%NUMBER_OF_NQ)
         ENDIF

         ! convert to MB 
         M_GTAU_G = M_GTAU_G* ndata *8._q / 1024/ 1024*PERCENT                                                      
      END SUBROUTINE STORAGE_FOR_GG_TAU

      !************************************************************************
      ! memory requirement for contraction of a greens functions and a
      ! self-energy in time and real space 
      !************************************************************************

      SUBROUTINE STORAGE_FOR_GTAU_WTAU( M_GTAU_WTAU, GDES_TAU, WDES, S2E, T, S, WDES_RESPONSE )
         REAL(q), INTENT(INOUT) :: M_GTAU_WTAU
         TYPE (greensfdes)      :: GDES_TAU
         TYPE (wavedes)         :: WDES
         TYPE( screened_2e_handle ) :: S2E
         TYPE(loop_des)         :: T
         TYPE( supercell ), POINTER ::  S
         TYPE( wavedes )       ::  WDES_RESPONSE
         !local 
         INTEGER                    :: NDIR, NQ, NK
         INTEGER                    :: NRPLWV_COL_DATA_POINTS
         INTEGER                    :: NRPLWV_ROW_DATA_POINTS
         REAL(q)                    :: M_TOTAL 

         M_GTAU_WTAU = 0   ! treated as peak memory usage for CALCLUATE_SIGMA_...
         M_TOTAL = 0             ! treated as current memory usage

         ! no computation of SIGMA for ACFDT calculations
         ! perform an early return 
         IF ( LACFDT .AND. .NOT. LRPAFORCE) RETURN 

         IF ( ICHIREAL == 1 ) THEN
            DO NQ=1,WDES%NKPTS 
            DO NK=1,S2E%NUMBER_OF_NQ
               CALL MEM_SIGMA_TAU( M_GTAU_WTAU, M_TOTAL, GDES_TAU)
            ENDDO         
            ENDDO         
         ELSE IF ( ICHIREAL == 2 ) THEN
            CALL MEM_SIGMA_SUPER( M_GTAU_WTAU, M_TOTAL, GDES_TAU, S, WDES, WDES_RESPONSE, S2E%NUMBER_OF_NQ )
         ENDIF

         ! convert to MB 
         M_GTAU_WTAU = M_GTAU_WTAU* ndata *8._q / 1024/ 1024*PERCENT                                                      
      END SUBROUTINE STORAGE_FOR_GTAU_WTAU

      !************************************************************************
      ! memory requirement for self-energy 
      !************************************************************************

      SUBROUTINE STORAGE_FOR_SIGMA( M_SIGMA, GDES_TAU, T, F, S2E, WDES )
         REAL(q), INTENT(INOUT)     :: M_SIGMA
         TYPE (greensfdes)          :: GDES_TAU
         TYPE(loop_des)             :: T, F
         TYPE( screened_2e_handle ) :: S2E
         TYPE( wavedes )            :: WDES
         ! local
         REAL(q)                    :: M_TOTAL
         REAL(q)                    :: RSIZE
         INTEGER                    :: NCPU
         INTEGER                    :: SIGMAW_MAT
         INTEGER                    :: GU_MAT, GO_MAT
         INTEGER                    :: CORR_MAT
         INTEGER                    :: CHAM_MAT
         INTEGER                    :: SIGMA_MAT
         INTEGER                    :: SIGMAO_MAT
         INTEGER                    :: SIGMAU_MAT
         INTEGER                    :: GAMMA_MAT
         INTEGER                    :: SIGMAU

         M_SIGMA = 0 
         M_TOTAL = 0 
         IF ( LACFDT .AND. .NOT. LRPAFORCE) RETURN 

         NCPU = 1
#ifdef MPI
         ! number of CPUs in a time group 
         NCPU = T%COMM_IN_GROUP%NCPU
#endif
         ! least amount of matrix elements master rank has to store
         RSIZE= (WDES%NB_TOT/NCPU + MOD( WDES%NB_TOT, NCPU )) ** 2 

         ! set size in bytes 
         ! SIGMAW_MAT is double complex 
         SIGMAW_MAT = 16*RSIZE*S2E%NUMBER_OF_NQ*WDES%ISPIN*T%NPOINTS_IN_ROOT_GROUP
         ! GU_MAT and GO_MAT are of size GDEF and remainder
         GU_MAT = wsgf * RSIZE*S2E%NUMBER_OF_NQ*WDES%ISPIN*T%NPOINTS_IN_ROOT_GROUP
         GO_MAT = wsgf * RSIZE*S2E%NUMBER_OF_NQ*WDES%ISPIN*T%NPOINTS_IN_ROOT_GROUP
         CORR_MAT = wsgf * RSIZE*S2E%NUMBER_OF_NQ*WDES%ISPIN
         SIGMA_MAT = wsgf * RSIZE*S2E%NUMBER_OF_NQ*WDES%ISPIN
         SIGMAO_MAT = wsgf * RSIZE
         SIGMAU_MAT = wsgf * RSIZE
         CHAM_MAT = wsgf * RSIZE*S2E%NUMBER_OF_NQ*WDES%ISPIN
         GAMMA_MAT = wsgf * RSIZE*S2E%NUMBER_OF_NQ*WDES%ISPIN
         ! required only for gamma-only version
         SIGMAU = wsgf * S2E%NUMBER_OF_NQ * &
               ( (GDES_TAU%NRPLWV_ROW_DATA_POINTS*GDES_TAU%NRPLWV_COL_DATA_POINTS)&
               +(GDES_TAU%NRPLWV_ROW_DATA_POINTS*GDES_TAU%NPRO_COL)              &
               +(GDES_TAU%NPRO_ROW*GDES_TAU%NPRO_COL) )
         
         
         IF ( ICHIREAL > 0 ) THEN

            ! in case self-consistent GW is done, that is Dyson equation
            ! is solved, additional matrices are required 
            IF (LGW .AND. ((.NOT. LG0W0 ) .AND.  (.NOT. LscQPGW))) THEN
               ! SIGMAW_MAT for self-energy in omega space 
               CALL REG_MEM( M_SIGMA, M_TOTAL, SIGMAW_MAT )
               IF ( NELMGW>1 ) THEN
                  ! GU_MAT and GO_MAT for self-consistent solutions of Dyson eq.
                  CALL REG_MEM( M_SIGMA, M_TOTAL, GO_MAT )
                  CALL REG_MEM( M_SIGMA, M_TOTAL, GU_MAT )
               ENDIF
            ENDIF 

            ! for all GW and RPAFORCE calculation 
            CALL REG_MEM( M_SIGMA, M_TOTAL, CORR_MAT )
            IF ( LRPAFORCE .AND. WDES%LOVERL ) THEN
               CALL REG_MEM( M_SIGMA, M_TOTAL, SIGMA_MAT )
            ENDIF

#ifdef gammareal
            ! SIGMAO and SIGMAU are 
            CALL REG_MEM( M_SIGMA, M_TOTAL, SIGMAU )
            CALL REG_MEM( M_SIGMA, M_TOTAL, SIGMAU )
#endif
            ! SIGMAO_MAT And SIGMAU_MAT
            CALL REG_MEM( M_SIGMA, M_TOTAL, SIGMAO_MAT ) 
            CALL REG_MEM( M_SIGMA, M_TOTAL, SIGMAU_MAT )
         
            ! CHAM_MAT
            CALL REG_MEM( M_SIGMA, M_TOTAL, CHAM_MAT )
           
            !GAMMA 
            CALL REG_MEM( M_SIGMA, M_TOTAL, GAMMA_MAT )

         ENDIF
         ! convert to MB 
         M_SIGMA = M_SIGMA / 1024/ 1024*PERCENT
         ! double memory for gamma-only calculations
         IF ( WDES%NKPTS == 1 ) THEN
            M_SIGMA = M_SIGMA * 2 
         ENDIF
      END SUBROUTINE STORAGE_FOR_SIGMA

      !*********************************************************************** 
      ! helper routine to set total and peak memory 
      !*********************************************************************** 
      SUBROUTINE  REG_MEM( PEAK_MEM, MEM_TOTAL, M ) 
         REAL(q), INTENT(INOUT) :: PEAK_MEM
         REAL(q), INTENT(INOUT) :: MEM_TOTAL
         INTEGER, INTENT(IN)    :: M
         
         MEM_TOTAL = MEM_TOTAL + REAL(M,q)

         IF (MEM_TOTAL>PEAK_MEM) THEN
            PEAK_MEM = MEM_TOTAL 
         ENDIF

      END SUBROUTINE  REG_MEM

      !*********************************************************************** 
      ! helper routine to deregister memory
      !*********************************************************************** 
      SUBROUTINE  DEREG_MEM( MEM_TOTAL, M ) 
         REAL(q), INTENT(INOUT) :: MEM_TOTAL
         INTEGER, INTENT(IN)    :: M
         
         MEM_TOTAL = MEM_TOTAL - REAL(M,q)

      END SUBROUTINE  DEREG_MEM

      !***********************************************************************
      ! this routine presents a skeleton version of 
      ! CALCULATE_CHI_TAU_FROM_G, where all large
      ! allocate and deallocate statements are replaced
      ! with REG_MEM and DEREG_MEM to compute the memory
      ! usage of one call to CALCULATE_CHI_TAU_FROM_G
      !***********************************************************************
      SUBROUTINE MEM_RESPONSE_GG( PEAK_MEM, MEM_TOTAL, GDES )
         REAL(q), INTENT(INOUT)  :: PEAK_MEM
         REAL(q), INTENT(INOUT)  :: MEM_TOTAL
         TYPE (greensfdes) :: GDES      ! descriptor for Green function
         ! local for projectors
         INTEGER :: P_G_PROJ
         INTEGER :: P_G_R
         INTEGER :: P_R_G
         INTEGER :: P_PROJ_G
 
         ! set sizes in bits
         P_G_PROJ =GDES%RES_NRPLWV_ROW_DATA_POINTS* GDES%NLM_COL
         P_G_R    =GDES%RES_NRPLWV_ROW_DATA_POINTS* GDES%MPLWV_COL
         P_PROJ_G =GDES%NLM_ROW*  GDES%RES_NRPLWV_COL_DATA_POINTS
         P_R_G    =GDES%MPLWV_ROW*GDES%RES_NRPLWV_COL_DATA_POINTS

         ! two calls to FFT_G
         CALL MEM_FFT_G( PEAK_MEM, MEM_TOTAL, GDES)
         CALL MEM_FFT_G( PEAK_MEM, MEM_TOTAL, GDES)

         ! allocation of P_G_PROJ and P_G_R  comes next
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_G_PROJ )
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_G_R )

         ! followed by two RELEASE_FFT_G calls
         CALL MEM_RELEASE_FFT_G( MEM_TOTAL, GDES )
         CALL MEM_RELEASE_FFT_G( MEM_TOTAL, GDES )

         ! restister again two auxillary arrays
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_PROJ_G )
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_R_G )

         ! transposition of TRANSPOSE_G_PROJ_RESPONSE
         ! is equivalent to 
         CALL MEM_TRANSPOSE_G_R_RESPONSE(PEAK_MEM, MEM_TOTAL, GDES )
         CALL MEM_TRANSPOSE_G_PROJ_RESPONSE(PEAK_MEM, MEM_TOTAL, GDES )

         ! derestister 
         CALL DEREG_MEM( MEM_TOTAL, P_G_PROJ )
         CALL DEREG_MEM( MEM_TOTAL, P_G_R )
         CALL DEREG_MEM( MEM_TOTAL, P_R_G )
         CALL DEREG_MEM( MEM_TOTAL, P_PROJ_G )

      END SUBROUTINE MEM_RESPONSE_GG
      !***********************************************************************
      ! memory requirement for FFT_G
      !***********************************************************************
      SUBROUTINE MEM_FFT_G( PEAK_MEM, MEM_TOTAL, GDES)
         REAL(q), INTENT(INOUT)  :: PEAK_MEM
         REAL(q), INTENT(INOUT)  :: MEM_TOTAL
         TYPE (greensfdes) :: GDES      ! descriptor for Green function
         ! local
         INTEGER           :: R_PROJ
         INTEGER           :: PROJ_R
         INTEGER           :: R_G

         ! set size in bits 
         R_PROJ=GDES%MPLWV_ROW*GDES%NPRO_COL
         PROJ_R=GDES%NPRO_ROW*GDES%MPLWV_COL
         R_G=GDES%MPLWV_ROW*GDES%NRPLWV_COL_DATA_POINTS

         ! ALLOCATE_R_G 
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, R_G )

         !ALLOCATE_G_R
         CALL MEM_G_R(PEAK_MEM, MEM_TOTAL, GDES)

         !TRANSPOSE_R_G
         CALL MEM_TRANSPOSE_R_G( PEAK_MEM, MEM_TOTAL, GDES)

         !DEALLOCATE_R_G
         CALL DEREG_MEM( MEM_TOTAL, R_G )
       
         !ALLOCATE_R_PROJ(GDES,G)
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, R_PROJ )

         !ALLOCATE_PROJ_R
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, PROJ_R )

         !TRANSPOSE_R_PROJ
         CALL MEM_TRANSPOSE_R_PROJ( PEAK_MEM, MEM_TOTAL, GDES)
         
      END SUBROUTINE MEM_FFT_G
      !***********************************************************************
      ! memory estimate for TRANSPOSE_G_R_RESPONSE
      !***********************************************************************
      SUBROUTINE MEM_TRANSPOSE_G_R_RESPONSE(PEAK_MEM, MEM_TOTAL, GDES)
         REAL(q), INTENT(INOUT)        :: PEAK_MEM
         REAL(q), INTENT(INOUT)        :: MEM_TOTAL
         TYPE (greensfdes), INTENT(IN) :: GDES      ! descriptor for Green function
#ifdef MPI
         ! local
         INTEGER :: P_G_R_LOCAL, P_G_R_RCV 

         P_G_R_LOCAL=GDES%MPLWV_COL_MAX* NSTRIP_STANDARD* GDES%COMM%NCPU
         P_G_R_RCV  =GDES%MPLWV_COL_MAX* NSTRIP_STANDARD* GDES%COMM%NCPU

         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_G_R_LOCAL )
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_G_R_RCV )
         CALL DEREG_MEM( MEM_TOTAL, P_G_R_LOCAL )
         CALL DEREG_MEM( MEM_TOTAL, P_G_R_RCV )
#endif
      END SUBROUTINE MEM_TRANSPOSE_G_R_RESPONSE
      !***********************************************************************
      ! memory estimate for TRANSPOSE_G_PROJ_RESPONSE
      !***********************************************************************
      SUBROUTINE MEM_TRANSPOSE_G_PROJ_RESPONSE(PEAK_MEM, MEM_TOTAL, GDES)
         REAL(q), INTENT(INOUT)        :: PEAK_MEM
         REAL(q), INTENT(INOUT)        :: MEM_TOTAL
         TYPE (greensfdes), INTENT(IN) :: GDES      ! descriptor for Green function
#ifdef MPI
         ! local
         INTEGER :: P_PROJ_R_LOCAL, P_PROJ_R_RCV 

         P_PROJ_R_LOCAL=GDES%NLM_COL_MAX* NSTRIP_STANDARD* GDES%COMM%NCPU
         P_PROJ_R_RCV  =GDES%NLM_COL_MAX* NSTRIP_STANDARD* GDES%COMM%NCPU
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_PROJ_R_LOCAL )
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_PROJ_R_RCV )
         CALL DEREG_MEM( MEM_TOTAL, P_PROJ_R_LOCAL )
         CALL DEREG_MEM( MEM_TOTAL, P_PROJ_R_RCV )
#endif
      END SUBROUTINE MEM_TRANSPOSE_G_PROJ_RESPONSE

      !***********************************************************************
      ! this routine presents a skeleton version of 
      ! CALCULATE_RESPONSE_SUPER, where all large
      ! allocate and deallocate statements are replaced
      ! with REG_MEM and DEREG_MEM to compute the memory
      ! usage of one call to CALCULATE_RESPONSE_SUPER
      !***********************************************************************
      SUBROUTINE MEM_RESPONSE_SUPER( PEAK_MEM, MEM_TOTAL, GDES, S, &
           WDES, NQMAX )
         REAL(q), INTENT(INOUT)  :: PEAK_MEM
         REAL(q), INTENT(INOUT)  :: MEM_TOTAL
         TYPE (greensfdes)   :: GDES      ! descriptor for Green function
         TYPE (wavedes)      :: WDES
         TYPE( supercell )   :: S
         INTEGER, INTENT(IN) :: NQMAX   ! total number of q-points 
         ! local for projectors
         INTEGER :: P_G_PROJ
         INTEGER :: P_G_R
         INTEGER :: P_R_G
         INTEGER :: P_PROJ_G
         INTEGER :: NK, NT, LMMAXC
         ! for GREENS function
         INTEGER :: GG
         INTEGER :: G_PROJ
         INTEGER :: PROJ_PROJ
         INTEGER :: GWORK_REAL, P_TMP 
 
         ! set sizes of arrays in bits
         P_G_PROJ = GDES%RES_NRPLWV_ROW_DATA_POINTS* GDES%NLM_COL* NQMAX
         P_G_R    = GDES%RES_NRPLWV_ROW_DATA_POINTS* GDES%MPLWV_COL* NQMAX
         P_R_G    = GDES%MPLWV_ROW* GDES%RES_NRPLWV_COL_DATA_POINTS
         P_PROJ_G = GDES%NLM_ROW* GDES%RES_NRPLWV_COL_DATA_POINTS
      
         ! first part where second index involves terms such as  
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_G_PROJ ) 

         DO NT=1,GDES%NTYP 
#ifdef RPAgamma
#ifndef gammareal
            GWORK_REAL=2 * S%WGW%GRID%MPLWV 
#else
            ! factor 2 taken into account by GDEF size, becaue this array is real
            GWORK_REAL= S%WGW%GRID%MPLWV 
#endif       
#else
            GWORK_REAL= m_ S%WGW%GRID%MPLWV
#endif
            GWORK_REAL = GWORK_REAL * GDES%NLM_LMMAX(NT)
            P_TMP=m_ S%WGW1%GRID%RC%NP
            ! register GWORK
            CALL REG_MEM( PEAK_MEM, MEM_TOTAL, GWORK_REAL ) 
            ! register P_TMP
            CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_TMP ) 
            CALL REG_MEM( PEAK_MEM, MEM_TOTAL, S%AUG_DES%NPRO* GDES%NLM_LMMAX(NT) ) 

            LMMAXC=GDES%LMMAX(NT)
            IF ( LMMAXC /=0 ) THEN
               CALL REG_MEM( PEAK_MEM, MEM_TOTAL, S%WDES1%NGVECTOR *LMMAXC)! GO_GK_PROJ
               CALL REG_MEM( PEAK_MEM, MEM_TOTAL, S%WDES1%GRID%RL%NP *LMMAXC)! GO_R_PROJ
               CALL REG_MEM( PEAK_MEM, MEM_TOTAL, S%WDES1%NPRO   *LMMAXC) !GO_PROJ_PROJ
               CALL REG_MEM( PEAK_MEM, MEM_TOTAL, S%WDES1%NGVECTOR *LMMAXC) !GU_GK_PROJ
               CALL REG_MEM( PEAK_MEM, MEM_TOTAL, S%WDES1%GRID%RL%NP *LMMAXC) !GU_R_PROJ
               CALL REG_MEM( PEAK_MEM, MEM_TOTAL, S%WDES1%GRID%RL%NP)  !GO_GU_PROJ
               CALL REG_MEM( PEAK_MEM, MEM_TOTAL, S%WDES1%NPRO*  LMMAXC) ! GU_PROJ_PROJ
               ! get rid of them again
               CALL DEREG_MEM( MEM_TOTAL, S%WDES1%NGVECTOR *LMMAXC)! GO_GK_PROJ
               CALL DEREG_MEM( MEM_TOTAL, S%WDES1%GRID%RL%NP *LMMAXC)! GO_R_PROJ
               CALL DEREG_MEM( MEM_TOTAL, S%WDES1%NPRO   *LMMAXC) !GO_PROJ_PROJ
               CALL DEREG_MEM( MEM_TOTAL, S%WDES1%NGVECTOR *LMMAXC) !GU_GK_PROJ
               CALL DEREG_MEM( MEM_TOTAL, S%WDES1%GRID%RL%NP *LMMAXC) !GU_R_PROJ
               CALL DEREG_MEM( MEM_TOTAL, S%WDES1%GRID%RL%NP)  !GO_GU_PROJ
               CALL DEREG_MEM( MEM_TOTAL, S%WDES1%NPRO*  LMMAXC) ! GU_PROJ_PROJ
            ENDIF

            ! deregister
            CALL DEREG_MEM( MEM_TOTAL, GWORK_REAL ) 
            CALL DEREG_MEM( MEM_TOTAL, P_TMP ) 
            CALL DEREG_MEM( MEM_TOTAL, S%AUG_DES%NPRO* GDES%NLM_LMMAX(NT) ) 
         ENDDO
         ! FFT and transpose orbitals at all k-points to obtain G%G_R and G%PROJ_R
         DO NK=1,WDES%NKPTS
            ! FFT GU and GO to real space and set entries G_R, PROJ_R
            CALL MEM_FFT_G_SUPER(PEAK_MEM, MEM_TOTAL, GDES, WDES, NK)
            CALL MEM_FFT_G_SUPER(PEAK_MEM, MEM_TOTAL, GDES, WDES, NK)
         ENDDO

         ! at this point G%GG and G%G_PROJ can be deallocated
         ! set sizes of arrays in bits 
         ! these are parts of the Green's function type G
         GG = GDES%NRPLWV_ROW_DATA_POINTS* GDES%NRPLWV_COL_DATA_POINTS
         G_PROJ = GDES%NRPLWV_ROW_DATA_POINTS* GDES%NPRO_COL
         PROJ_PROJ = GDES%NPRO_ROW*GDES%NPRO_COL
         DO NK=1,WDES%NKPTS
            CALL DEREG_MEM( MEM_TOTAL, GG )
            CALL DEREG_MEM( MEM_TOTAL, GG )
            CALL DEREG_MEM( MEM_TOTAL, G_PROJ )
            CALL DEREG_MEM( MEM_TOTAL, G_PROJ )
            CALL DEREG_MEM( MEM_TOTAL, PROJ_PROJ )
            CALL DEREG_MEM( MEM_TOTAL, PROJ_PROJ )
         ENDDO
 
         ! now P_G_R is allocated in routine 
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_G_R )

         ! register GWORK P_TMP
#ifdef RPAgamma
#ifndef gammareal
         GWORK_REAL=2 * S%WGW%GRID%MPLWV 
#else
         ! factor 2 taken into account by GDEF size, becaue this array is real
         GWORK_REAL= S%WGW%GRID%MPLWV 
#endif    
#else
         GWORK_REAL=m_ S%WGW%GRID%MPLWV
#endif
         P_TMP=m_ S%WGW1%GRID%RC%NP
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, GWORK_REAL ) 
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_TMP ) 
         CALL DEREG_MEM(  MEM_TOTAL, GWORK_REAL ) 
         CALL DEREG_MEM(  MEM_TOTAL, P_TMP ) 

         ! free the Green's functions
         DO NK=1,WDES%NKPTS
            ! deallocate G_R, PROJ_R, and R_PROJ
            CALL MEM_RELEASE_FFT_G(MEM_TOTAL, GDES, NK )
            CALL MEM_RELEASE_FFT_G(MEM_TOTAL, GDES, NK )
         ENDDO
  
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_R_G )
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_PROJ_G )

         ! deregister remaining temporary arrays
         CALL DEREG_MEM( MEM_TOTAL, P_R_G )
         CALL DEREG_MEM( MEM_TOTAL, P_PROJ_G )

         CALL DEREG_MEM( MEM_TOTAL, P_G_R )
         CALL DEREG_MEM( MEM_TOTAL, P_G_PROJ )
    
      END SUBROUTINE MEM_RESPONSE_SUPER
      !***********************************************************************
      ! computes memory usage of FFT_G_SUPER
      !***********************************************************************
      SUBROUTINE MEM_FFT_G_SUPER( PEAK_MEM, MEM_TOTAL, GDES, WDES, NK)
         REAL(q), INTENT(INOUT)  :: PEAK_MEM
         REAL(q), INTENT(INOUT)  :: MEM_TOTAL
         TYPE (greensfdes),INTENT(IN) :: GDES         
         TYPE (wavedes),INTENT(IN) :: WDES
         INTEGER, OPTIONAL, INTENT(IN) :: NK
         ! local 
         INTEGER  :: NRPLWV_ROW_DATA_POINTS
         INTEGER  :: NRPLWV_COL_DATA_POINTS
         INTEGER  :: R_G, G_R, R_PROJ, PROJ_R 

#ifdef gammareal
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, 4*WDES%GRID%RL%NP )
#else
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, 2*WDES%GRID%RL%NP )
#endif

         NRPLWV_ROW_DATA_POINTS=GDES%NRPLWV_ROW_DATA_POINTS
         NRPLWV_COL_DATA_POINTS=GDES%NRPLWV_COL_DATA_POINTS
         IF (ASSOCIATED(GDES%LUSEINV)) THEN
            IF (GDES%LUSEINV(NK)) THEN
               NRPLWV_ROW_DATA_POINTS=GDES%NRPLWV_ROW_DATA_POINTS_NK(NK)
               NRPLWV_COL_DATA_POINTS=GDES%NRPLWV_COL_DATA_POINTS_NK(NK)
            ENDIF
         ENDIF

         R_G=GDES%MPLWV_ROW* NRPLWV_COL_DATA_POINTS
         G_R=NRPLWV_ROW_DATA_POINTS* GDES%MPLWV_COL
         R_PROJ=GDES%MPLWV_ROW* GDES%NPRO_COL
         PROJ_R=GDES%NPRO_ROW* GDES%MPLWV_COL
        
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, R_G)
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, G_R)

         CALL MEM_TRANSPOSE_R_G( PEAK_MEM, MEM_TOTAL, GDES, NK)
    
         ! we do not need R_G anymore, so destroy it
         CALL DEREG_MEM( MEM_TOTAL, R_G)
    
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, R_PROJ)
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, PROJ_R)

         CALL MEM_TRANSPOSE_R_PROJ( PEAK_MEM, MEM_TOTAL, GDES)
    
         ! we do not need R_PROJ anymore, so destroy it
         CALL DEREG_MEM( MEM_TOTAL, R_PROJ )
         CALL DEREG_MEM( MEM_TOTAL, PROJ_R) 
    
#ifdef gammareal
         CALL DEREG_MEM( MEM_TOTAL, 4*WDES%GRID%RL%NP )
#else
         CALL DEREG_MEM( MEM_TOTAL, 2*WDES%GRID%RL%NP )
#endif
      END SUBROUTINE MEM_FFT_G_SUPER
      !***********************************************************************
      ! memory for transposition of R_G
      !***********************************************************************
      SUBROUTINE MEM_TRANSPOSE_R_G ( PEAK_MEM, MEM_TOTAL, GDES, NK)
         REAL(q), INTENT(INOUT)  :: PEAK_MEM
         REAL(q), INTENT(INOUT)  :: MEM_TOTAL
         TYPE (greensfdes),INTENT(IN)  :: GDES         
         INTEGER, OPTIONAL, INTENT(IN) :: NK
#ifdef MPI
         ! local 
         INTEGER :: NRPLWV_ROW_DATA_POINTS
         INTEGER :: NRPLWV_COL_DATA_POINTS
         INTEGER :: NRPLWV_COL_MAX_DATA_POINTS
         INTEGER :: P_G_R_LOCAL , P_G_R_RCV
 
         NRPLWV_ROW_DATA_POINTS=GDES%NRPLWV_ROW_DATA_POINTS
         NRPLWV_COL_MAX_DATA_POINTS=GDES%NRPLWV_COL_MAX_DATA_POINTS
         NRPLWV_COL_DATA_POINTS=GDES%NRPLWV_COL_DATA_POINTS
         IF (PRESENT(NK) .AND. ASSOCIATED(GDES%LUSEINV)) THEN
            IF (GDES%LUSEINV(NK)) THEN
               NRPLWV_ROW_DATA_POINTS=GDES%NRPLWV_ROW_DATA_POINTS_NK(NK)
               NRPLWV_COL_MAX_DATA_POINTS=GDES%NRPLWV_COL_MAX_DATA_POINTS_NK(NK)
               NRPLWV_COL_DATA_POINTS=GDES%NRPLWV_COL_DATA_POINTS_NK(NK)
            ENDIF
         ENDIF
 
         ! communicate up to NSTRIP rows at a time
         P_G_R_LOCAL=NRPLWV_COL_MAX_DATA_POINTS* NSTRIP_STANDARD* GDES%COMM%NCPU
         P_G_R_RCV  =NRPLWV_COL_MAX_DATA_POINTS* NSTRIP_STANDARD* GDES%COMM%NCPU
 
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_G_R_LOCAL )
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_G_R_RCV )
 
         CALL DEREG_MEM( MEM_TOTAL, P_G_R_LOCAL )
         CALL DEREG_MEM( MEM_TOTAL, P_G_R_RCV )
#endif
      END SUBROUTINE MEM_TRANSPOSE_R_G
      !***********************************************************************
      ! memory for transposition of R_PROJ
      !***********************************************************************
      SUBROUTINE MEM_TRANSPOSE_R_PROJ( PEAK_MEM, MEM_TOTAL, GDES)
         REAL(q), INTENT(INOUT)  :: PEAK_MEM
         REAL(q), INTENT(INOUT)  :: MEM_TOTAL
         TYPE (greensfdes),INTENT(IN)  :: GDES         
#ifdef MPI
         ! local
         INTEGER :: P_PROJ_R_LOCAL, P_PROJ_R_RCV
   
         P_PROJ_R_LOCAL=GDES%NPRO_COL_MAX* NSTRIP_STANDARD* GDES%COMM%NCPU
         P_PROJ_R_RCV  =GDES%NPRO_COL_MAX* NSTRIP_STANDARD* GDES%COMM%NCPU
   
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_PROJ_R_LOCAL )
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_PROJ_R_RCV )

         CALL DEREG_MEM( MEM_TOTAL, P_PROJ_R_LOCAL )
         CALL DEREG_MEM( MEM_TOTAL, P_PROJ_R_RCV )
#endif
      END SUBROUTINE MEM_TRANSPOSE_R_PROJ
      !***********************************************************************
      ! helper that frees memory usage after RELEASE_FFT_G
      !***********************************************************************
      SUBROUTINE MEM_RELEASE_FFT_G( MEM_TOTAL, GDES, NK )
         REAL(q), INTENT(INOUT)        :: MEM_TOTAL 
         TYPE (greensfdes), INTENT(IN) :: GDES
         INTEGER, OPTIONAL, INTENT(IN) :: NK
         ! local
         INTEGER :: NRPLWV_ROW_DATA_POINTS
         INTEGER :: G_R, PROJ_R, R_PROJ

         ! storage for G_R 
         NRPLWV_ROW_DATA_POINTS=GDES%NRPLWV_ROW_DATA_POINTS
         IF (PRESENT(NK) .AND. ASSOCIATED(GDES%LUSEINV)) THEN
            IF (GDES%LUSEINV(NK)) THEN
               NRPLWV_ROW_DATA_POINTS=GDES%NRPLWV_ROW_DATA_POINTS_NK(NK)
            ENDIF
         ENDIF
         G_R=NRPLWV_ROW_DATA_POINTS*GDES%MPLWV_COL
         PROJ_R=GDES%NPRO_ROW*GDES%MPLWV_COL
         R_PROJ=GDES%MPLWV_ROW*GDES%NPRO_COL

         CALL DEREG_MEM( MEM_TOTAL, G_R )
         CALL DEREG_MEM( MEM_TOTAL, PROJ_R )
         CALL DEREG_MEM( MEM_TOTAL, R_PROJ )
    
      END SUBROUTINE MEM_RELEASE_FFT_G

      !***********************************************************************
      !
      ! computes memory requirement for CALCULATE_SIGMA_SUPER
      !
      !***********************************************************************
      SUBROUTINE MEM_SIGMA_SUPER( PEAK_MEM, MEM_TOTAL, GDES, S, WDES, WGW, NQMAX)
         REAL(q), INTENT(INOUT)        :: PEAK_MEM
         REAL(q), INTENT(INOUT)        :: MEM_TOTAL
         TYPE (greensfdes), INTENT(IN) :: GDES      ! descriptor for Green function
         TYPE (supercell)              :: S
         TYPE (wavedes)                :: WDES
         TYPE (wavedes)                :: WGW
         INTEGER, INTENT(IN)           :: NQMAX     ! total number of q-points 
         ! local 
         INTEGER :: NDIR, NK, NQ, IDIR,NT
         TYPE (wavedes1)                :: WDESQ
         INTEGER :: P_G_R
         INTEGER :: P_R_G
         INTEGER :: P_PROJ_G
         INTEGER :: P_G_PROJ

         INTEGER :: GU_R_PROJ
         INTEGER :: SIG_R_PROJ 
         INTEGER :: SIG_PROJ_PROJ 
         INTEGER :: GU_GK_PROJ 
         INTEGER :: GU_PROJ_PROJ 
         INTEGER :: D_G_PROJ 
         INTEGER :: D_R_PROJ 
         INTEGER :: D_PROJ_PROJ

         INTEGER :: GWORK, GWORK_REAL, P_TMP
         INTEGER :: LMMAXC, NLM_LMMAXC
         
         ! to seek further changes seek IDIR
         IF (LRPAFORCE) THEN
            NDIR=4
         ELSE
            NDIR=1
         ENDIF
         ! set size in bits
         P_G_R    = GDES%RES_NRPLWV_ROW_DATA_POINTS* GDES%MPLWV_COL* WDES%NKPTS
         P_R_G    = GDES%MPLWV_ROW* GDES%RES_NRPLWV_COL_DATA_POINTS* WDES%NKPTS
         P_PROJ_G = GDES%NLM_ROW* GDES%RES_NRPLWV_COL_DATA_POINTS * WDES%NKPTS*NDIR
         P_G_PROJ = GDES%RES_NRPLWV_ROW_DATA_POINTS* GDES%NLM_COL * WDES%NKPTS*NDIR
    
         DO NK=1,WDES%NKPTS
            ! FFT and transpose orbitals at all k-points to obtain G%G_r and G%PROJ_r
            ! FFT GU to real space and set entries G'_r, PROJ'_r
            CALL MEM_FFT_G_SUPER(PEAK_MEM, MEM_TOTAL, GDES, WDES, NK)
         ENDDO
    
         CALL REG_MEM(PEAK_MEM, MEM_TOTAL, P_G_R )
         CALL REG_MEM(PEAK_MEM, MEM_TOTAL, P_R_G )
         CALL REG_MEM(PEAK_MEM, MEM_TOTAL, P_PROJ_G )
         CALL REG_MEM(PEAK_MEM, MEM_TOTAL, P_G_PROJ )

    
         DO NK=1,WDES%NKPTS
            CALL SETWDES(WDES_RESPONSE, WDESQ, NK)
            GWORK = m_ WDESQ%GRID%MPLWV
            CALL REG_MEM(PEAK_MEM, MEM_TOTAL, GWORK )
            CALL DEREG_MEM( MEM_TOTAL, GWORK )
          
            CALL MEM_TRANSPOSE_R_G_RESPONSE(PEAK_MEM, MEM_TOTAL, GDES )
            !DO IDIR=1,NDIR
               CALL MEM_TRANSPOSE_PROJ_G_RESPONSE(PEAK_MEM, MEM_TOTAL, GDES )
            !ENDDO
         ENDDO
         
         CALL DEREG_MEM(MEM_TOTAL, P_R_G )
         CALL DEREG_MEM(MEM_TOTAL, P_PROJ_G )

         ! part 2, supercell
 
         GWORK = m_ S%WDES%GRID%MPLWV
         P_TMP = m_ S%WDES%GRID%RC%NP
         CALL REG_MEM(PEAK_MEM, MEM_TOTAL, GWORK )
         CALL REG_MEM(PEAK_MEM, MEM_TOTAL, P_TMP )

    
         DO NK=1, NQMAX ! SIGMA is allocated only in the IBZ
            CALL MEM_G_R(PEAK_MEM, MEM_TOTAL, GDES, NK)
            CALL MEM_G_PROJ(PEAK_MEM, MEM_TOTAL, GDES)
            CALL MEM_PROJ_PROJ(PEAK_MEM, MEM_TOTAL, GDES)
         ENDDO


         !contraction part over realspace grid points 
#ifdef RPAgamma
#ifndef gammareal
         GWORK_REAL=2 * S%WGW%GRID%MPLWV 
#else
         ! factor 2 taken into account by GDEF size, becaue this array is real
         GWORK_REAL= S%WGW%GRID%MPLWV 
#endif    
#else
         GWORK_REAL=m_ S%WGW%GRID%MPLWV
#endif
         CALL REG_MEM(PEAK_MEM, MEM_TOTAL, GWORK_REAL )
         CALL DEREG_MEM( MEM_TOTAL, GWORK_REAL )
         CALL DEREG_MEM(MEM_TOTAL, P_TMP )
         CALL DEREG_MEM(MEM_TOTAL, P_G_R )

         ! now the awkward part
         IF (S%WDES1%LOVERL) THEN
         DO IDIR=1,NDIR
            DO NT=1,GDES%NTYP 
               LMMAXC=GDES%LMMAX(NT)
               NLM_LMMAXC=GDES%NLM_LMMAX(NT)
               IF (LMMAXC.NE.0) THEN
                  GU_R_PROJ=S%WDES1%GRID%RL%NP* LMMAXC
                  SIG_R_PROJ=S%WDES1%GRID%RL%NP* LMMAXC
                  SIG_PROJ_PROJ=S%WDES1%NPRO* LMMAXC
                  GU_GK_PROJ=S%WDES1%NGVECTOR *LMMAXC
                  GU_PROJ_PROJ=S%WDES1%NPRO *LMMAXC
                  D_G_PROJ=S%WGW1%NGVECTOR* NLM_LMMAXC
#ifdef RPAgamma
#ifndef gammareal
                  D_R_PROJ=S%WGW%GRID%MPLWV*2* NLM_LMMAXC
#else
                  D_R_PROJ=S%WGW%GRID%MPLWV* NLM_LMMAXC
#endif
#else
                  D_R_PROJ=S%WGW%GRID%MPLWV* NLM_LMMAXC
#endif
                  D_PROJ_PROJ=S%AUG_DES%NPRO *NLM_LMMAXC
                        
                  CALL REG_MEM( PEAK_MEM, MEM_TOTAL, GU_R_PROJ  )
                  CALL REG_MEM( PEAK_MEM, MEM_TOTAL, SIG_R_PROJ ) 
                  CALL REG_MEM( PEAK_MEM, MEM_TOTAL, SIG_PROJ_PROJ )
                  CALL REG_MEM( PEAK_MEM, MEM_TOTAL, GU_GK_PROJ )
                  CALL REG_MEM( PEAK_MEM, MEM_TOTAL, GU_PROJ_PROJ )
                  CALL REG_MEM( PEAK_MEM, MEM_TOTAL, D_G_PROJ ) 
                  CALL REG_MEM( PEAK_MEM, MEM_TOTAL, D_R_PROJ )
                  CALL REG_MEM( PEAK_MEM, MEM_TOTAL, D_PROJ_PROJ )

                  CALL DEREG_MEM( MEM_TOTAL, GU_R_PROJ  )
                  CALL DEREG_MEM( MEM_TOTAL, SIG_R_PROJ ) 
                  CALL DEREG_MEM( MEM_TOTAL, SIG_PROJ_PROJ )
                  CALL DEREG_MEM( MEM_TOTAL, GU_GK_PROJ )
                  CALL DEREG_MEM( MEM_TOTAL, GU_PROJ_PROJ )
                  CALL DEREG_MEM( MEM_TOTAL, D_G_PROJ ) 
                  CALL DEREG_MEM( MEM_TOTAL, D_R_PROJ )
                  CALL DEREG_MEM( MEM_TOTAL, D_PROJ_PROJ )
               ENDIF
            ENDDO 
         ENDDO 
         ENDIF

         CALL DEREG_MEM(MEM_TOTAL, GWORK )
         CALL DEREG_MEM(MEM_TOTAL, P_G_PROJ )
    
         DO NK=1,WDES%NKPTS
            ! release G_R, PROJ_R, R_PROJ
            CALL MEM_RELEASE_FFT_G(MEM_TOTAL, GDES, NK )
         ENDDO
         DO NQ=1,NQMAX
            ! transpose sigma and fft in other direction to gain Sigma(g,g') for each k
            CALL MEM_FFT_SIGMA_SUPER ( PEAK_MEM, MEM_TOTAL, GDES, NQ )
#ifdef RPAgamma
            ! if only half of the coefficients are available in second direction
            ! uncompress using inversion symmetry
            CALL MEM_UNCOMPRESS_G_RECIPROCAL( PEAK_MEM, MEM_TOTAL, GDES, NQ)
#endif
         ENDDO !NQ loop
      END SUBROUTINE MEM_SIGMA_SUPER
      !***********************************************************************
      ! memory estimate for TRANSPOSE_R_G_RESPONSE
      !***********************************************************************
      SUBROUTINE MEM_TRANSPOSE_R_G_RESPONSE(PEAK_MEM, MEM_TOTAL, GDES)
         REAL(q), INTENT(INOUT)        :: PEAK_MEM
         REAL(q), INTENT(INOUT)        :: MEM_TOTAL
         TYPE (greensfdes), INTENT(IN) :: GDES      ! descriptor for Green function
#ifdef MPI
         ! local
         INTEGER :: P_G_R_LOCAL, P_G_R_RCV 

         P_G_R_LOCAL=GDES%RES_NRPLWV_COL_MAX_DATA_POINTS* NSTRIP_STANDARD* GDES%COMM%NCPU
         P_G_R_RCV  =GDES%RES_NRPLWV_COL_MAX_DATA_POINTS* NSTRIP_STANDARD* GDES%COMM%NCPU
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_G_R_LOCAL )
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_G_R_RCV )
         CALL DEREG_MEM( MEM_TOTAL, P_G_R_LOCAL )
         CALL DEREG_MEM( MEM_TOTAL, P_G_R_RCV )
#endif
      END SUBROUTINE MEM_TRANSPOSE_R_G_RESPONSE
      !***********************************************************************
      ! memory estimate for TRANSPOSE_PROJ_G_RESPONSE
      !***********************************************************************
      SUBROUTINE MEM_TRANSPOSE_PROJ_G_RESPONSE(PEAK_MEM, MEM_TOTAL, GDES)
         REAL(q), INTENT(INOUT)        :: PEAK_MEM
         REAL(q), INTENT(INOUT)        :: MEM_TOTAL
         TYPE (greensfdes), INTENT(IN) :: GDES      ! descriptor for Green function

         ! is equivalent to G_RPOJ_RESPONSE transpotions
         CALL MEM_TRANSPOSE_G_PROJ_RESPONSE(PEAK_MEM, MEM_TOTAL, GDES )
      END SUBROUTINE MEM_TRANSPOSE_PROJ_G_RESPONSE
      !***********************************************************************
      ! memory estimate for ALLOCATE_G_R
      !***********************************************************************
      SUBROUTINE MEM_G_R(PEAK_MEM, MEM_TOTAL, GDES, NK)
         REAL(q), INTENT(INOUT)        :: PEAK_MEM
         REAL(q), INTENT(INOUT)        :: MEM_TOTAL
         TYPE (greensfdes), INTENT(IN) :: GDES      ! descriptor for Green function
         INTEGER, OPTIONAL, INTENT(IN) :: NK
         ! local 
         INTEGER :: NRPLWV_ROW_DATA_POINTS
         INTEGER :: G_R

         NRPLWV_ROW_DATA_POINTS=GDES%NRPLWV_ROW_DATA_POINTS
         IF (PRESENT(NK) .AND. ASSOCIATED(GDES%LUSEINV)) THEN
            IF (GDES%LUSEINV(NK)) THEN
               NRPLWV_ROW_DATA_POINTS=GDES%NRPLWV_ROW_DATA_POINTS_NK(NK)
            ENDIF
         ENDIF
         G_R=NRPLWV_ROW_DATA_POINTS*GDES%MPLWV_COL
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, G_R )
      END SUBROUTINE MEM_G_R
      !***********************************************************************
      ! memory estimate for ALLOCATE_G_PROJ
      !***********************************************************************
      SUBROUTINE MEM_G_PROJ(PEAK_MEM, MEM_TOTAL, GDES)
         REAL(q), INTENT(INOUT)        :: PEAK_MEM
         REAL(q), INTENT(INOUT)        :: MEM_TOTAL
         TYPE (greensfdes), INTENT(IN) :: GDES      ! descriptor for Green function
         ! local 
         INTEGER :: G_PROJ

         G_PROJ=GDES%NRPLWV_ROW_DATA_POINTS* GDES%NPRO_COL
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, G_PROJ )
      END SUBROUTINE MEM_G_PROJ
      !***********************************************************************
      ! memory estimate for ALLOCATE_PROJ_PROJ
      !***********************************************************************
      SUBROUTINE MEM_PROJ_PROJ(PEAK_MEM, MEM_TOTAL, GDES)
         REAL(q), INTENT(INOUT)        :: PEAK_MEM
         REAL(q), INTENT(INOUT)        :: MEM_TOTAL
         TYPE (greensfdes), INTENT(IN) :: GDES      ! descriptor for Green function
         ! local 
         INTEGER :: PROJ_PROJ

         PROJ_PROJ = GDES%NPRO_ROW*GDES%NPRO_COL
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, PROJ_PROJ )
      END SUBROUTINE MEM_PROJ_PROJ
      !***********************************************************************
      ! memory estimate for FFT_SIGMA_SUPER
      !***********************************************************************
      SUBROUTINE MEM_FFT_SIGMA_SUPER( PEAK_MEM, MEM_TOTAL, GDES, NK )
         REAL(q), INTENT(INOUT)        :: PEAK_MEM
         REAL(q), INTENT(INOUT)        :: MEM_TOTAL
         TYPE (greensfdes), INTENT(IN) :: GDES      ! descriptor for Green function
         INTEGER, OPTIONAL, INTENT(IN) :: NK
         ! local 
         INTEGER                       :: NRPLWV_COL_DATA_POINTS
         INTEGER                       :: NRPLWV_ROW_DATA_POINTS
         INTEGER                       :: G_R, R_G, GG

         NRPLWV_COL_DATA_POINTS=GDES%NRPLWV_COL_DATA_POINTS
         NRPLWV_ROW_DATA_POINTS=GDES%NRPLWV_ROW_DATA_POINTS
         IF (PRESENT(NK) .AND. ASSOCIATED(GDES%LUSEINV)) THEN
            IF (GDES%LUSEINV(NK)) THEN
               NRPLWV_ROW_DATA_POINTS=GDES%NRPLWV_ROW_DATA_POINTS_NK(NK)
               NRPLWV_COL_DATA_POINTS=GDES%NRPLWV_COL_DATA_POINTS_NK(NK)
            ENDIF
         ENDIF
         ! ALLOCATE_R_G 
         R_G=GDES%MPLWV_ROW*NRPLWV_COL_DATA_POINTS
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, R_G )

         ! transpose G_R
         CALL MEM_TRANSPOSE_G_R(PEAK_MEM, MEM_TOTAL, GDES)

         !DEALLOCATE_G_R
         G_R=NRPLWV_ROW_DATA_POINTS*GDES%MPLWV_COL
         CALL DEREG_MEM( MEM_TOTAL, G_R )

         ! allocate GG
         GG=GDES%NRPLWV_ROW_DATA_POINTS*GDES%NRPLWV_COL_DATA_POINTS
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, GG )

         ! deallocate R_G
         CALL DEREG_MEM( MEM_TOTAL, R_G)
      END SUBROUTINE MEM_FFT_SIGMA_SUPER
      !***********************************************************************
      ! memory estimate for TRANSPOSE_G_R
      !***********************************************************************
      SUBROUTINE MEM_TRANSPOSE_G_R( PEAK_MEM, MEM_TOTAL, GDES)
         REAL(q), INTENT(INOUT)        :: PEAK_MEM
         REAL(q), INTENT(INOUT)        :: MEM_TOTAL
         TYPE (greensfdes), INTENT(IN) :: GDES      ! descriptor for Green function
       ! local
         INTEGER :: P_R_G_LOCAL
         INTEGER :: P_R_G_RCV
#ifdef MPI
         P_R_G_LOCAL=GDES%MPLWV_COL_MAX* NSTRIP_STANDARD* GDES%COMM%NCPU
         P_R_G_RCV  =GDES%MPLWV_COL_MAX* NSTRIP_STANDARD* GDES%COMM%NCPU
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_R_G_LOCAL )
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_R_G_RCV )
         CALL DEREG_MEM( MEM_TOTAL, P_R_G_LOCAL )
         CALL DEREG_MEM( MEM_TOTAL, P_R_G_RCV )
#endif
      END SUBROUTINE MEM_TRANSPOSE_G_R
      !***********************************************************************
      ! memory estimate for UNCOMPRESS_G_RECIPROCAL
      !***********************************************************************
      SUBROUTINE MEM_UNCOMPRESS_G_RECIPROCAL( PEAK_MEM, MEM_TOTAL, GDES, NK1)
         REAL(q), INTENT(INOUT)        :: PEAK_MEM
         REAL(q), INTENT(INOUT)        :: MEM_TOTAL
         TYPE (greensfdes), INTENT(IN) :: GDES      ! descriptor for Green function
         INTEGER, INTENT(IN)           :: NK1
         ! local 
         INTEGER                       :: GG1, GG2
         IF (ASSOCIATED(GDES%LUSEINV)) THEN
         IF (GDES%LUSEINV(NK1)) THEN
            GG1=GDES%NRPLWV_ROW_DATA_POINTS*GDES%NRPLWV_COL_DATA_POINTS
            GG2=GDES%NRPLWV_ROW_DATA_POINTS*GDES%NRPLWV_COL_DATA_POINTS
            ! register auxillary arrays 
            CALL REG_MEM( PEAK_MEM, MEM_TOTAL, GG1 )
            CALL REG_MEM( PEAK_MEM, MEM_TOTAL, GG2 )
            CALL MEM_TRANSPOSE_G_G( PEAK_MEM, MEM_TOTAL, GDES, NK1 )
            CALL MEM_TRANSPOSE_G_G( PEAK_MEM, MEM_TOTAL, GDES, NK1 )
            ! transpose again
            CALL MEM_TRANSPOSE_G_G( PEAK_MEM, MEM_TOTAL, GDES )
            CALL DEREG_MEM( MEM_TOTAL, GG1 )
            CALL DEREG_MEM( MEM_TOTAL, GG2 )
         ENDIF
         ENDIF
      END SUBROUTINE MEM_UNCOMPRESS_G_RECIPROCAL
      !***********************************************************************
      ! memory estimate for TRANSPOSE_G_G
      !***********************************************************************
      SUBROUTINE MEM_TRANSPOSE_G_G( PEAK_MEM, MEM_TOTAL, GDES, NK2 )
         REAL(q), INTENT(INOUT)        :: PEAK_MEM
         REAL(q), INTENT(INOUT)        :: MEM_TOTAL 
         TYPE (greensfdes), INTENT(IN) :: GDES
         INTEGER, OPTIONAL, INTENT(IN) :: NK2
#ifdef MPI
      ! local
        INTEGER :: P_LOCAL
        INTEGER :: P_RCV
        INTEGER :: NRPLWV_ROW_DATA_POINTS2
        INTEGER :: NRPLWV_COL_DATA_POINTS2
        INTEGER :: NRPLWV_COL_MAX_DATA_POINTS2
        
        NRPLWV_ROW_DATA_POINTS2    =GDES%NRPLWV_ROW_DATA_POINTS
        NRPLWV_COL_MAX_DATA_POINTS2=GDES%NRPLWV_COL_MAX_DATA_POINTS
        NRPLWV_COL_DATA_POINTS2    =GDES%NRPLWV_COL_DATA_POINTS
        IF (PRESENT(NK2) .AND. ASSOCIATED(GDES%LUSEINV) ) THEN
           IF (GDES%LUSEINV(NK2)) THEN
              NRPLWV_ROW_DATA_POINTS2    =GDES%NRPLWV_ROW_DATA_POINTS_NK(NK2)
              NRPLWV_COL_MAX_DATA_POINTS2=GDES%NRPLWV_COL_MAX_DATA_POINTS_NK(NK2)
              NRPLWV_COL_DATA_POINTS2    =GDES%NRPLWV_COL_DATA_POINTS_NK(NK2)
           ENDIF
        ENDIF
    
        P_LOCAL=NRPLWV_COL_MAX_DATA_POINTS2* NSTRIP_STANDARD*GDES%COMM%NCPU
        P_RCV  =NRPLWV_COL_MAX_DATA_POINTS2* NSTRIP_STANDARD*GDES%COMM%NCPU
        CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_LOCAL )
        CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_RCV )
    
        CALL DEREG_MEM( MEM_TOTAL, P_LOCAL )
        CALL DEREG_MEM( MEM_TOTAL, P_RCV )
#endif
      END SUBROUTINE MEM_TRANSPOSE_G_G
      !***********************************************************************
      !  
      ! computes memory requirement for CALCULATE_SIGMA_TAU
      ! used for ALGO = ...RK  algorithms
      !  
      !***********************************************************************
      SUBROUTINE MEM_SIGMA_TAU( PEAK_MEM, MEM_TOTAL, GDES)
         REAL(q), INTENT(INOUT)        :: PEAK_MEM
         REAL(q), INTENT(INOUT)        :: MEM_TOTAL
         TYPE (greensfdes), INTENT(IN) :: GDES      ! descriptor for Green function
         ! local 
         INTEGER :: NDIR
         INTEGER :: P_G_R
         INTEGER :: P_R_G
         INTEGER :: P_PROJ_G
         INTEGER :: P_G_PROJ
         INTEGER :: R_PROJ
         INTEGER :: G_PROJ
         INTEGER :: PROJ_PROJ
         
         ! take RPAforce calculation into account 
         IF (LRPAFORCE) THEN
            NDIR=4
         ELSE
            NDIR=1
         ENDIF
         ! set size in bits
         P_PROJ_G =GDES%NLM_ROW* GDES%RES_NRPLWV_COL_DATA_POINTS*NDIR
         P_R_G    =GDES%MPLWV_ROW* GDES%RES_NRPLWV_COL_DATA_POINTS
         P_G_PROJ =GDES%RES_NRPLWV_ROW_DATA_POINTS* GDES%NLM_COL* NDIR
         P_G_R    =GDES%RES_NRPLWV_ROW_DATA_POINTS* GDES%MPLWV_COL
 
         ! FFT GO to real space and set entries G_R, R_PROJ, PROJ_R
         CALL MEM_FFT_G( PEAK_MEM, MEM_TOTAL, GDES)

         ! allocate aux arrays 
         CALL REG_MEM(PEAK_MEM, MEM_TOTAL, P_PROJ_G )
         CALL REG_MEM(PEAK_MEM, MEM_TOTAL, P_R_G )
         CALL REG_MEM(PEAK_MEM, MEM_TOTAL, P_G_PROJ )
         CALL REG_MEM(PEAK_MEM, MEM_TOTAL, P_G_R )

         CALL MEM_TRANSPOSE_R_G_RESPONSE(PEAK_MEM, MEM_TOTAL, GDES )
         !DO IDIR=1,NDIR
            CALL MEM_TRANSPOSE_PROJ_G_RESPONSE(PEAK_MEM, MEM_TOTAL, GDES )
         !ENDDO

         CALL DEREG_MEM(MEM_TOTAL, P_R_G )
         CALL DEREG_MEM(MEM_TOTAL, P_PROJ_G )

         !ALLOCATE_G_R
         CALL MEM_G_R(PEAK_MEM, MEM_TOTAL, GDES)
         CALL DEREG_MEM(MEM_TOTAL, P_G_R )

         ! part 3
         ! ALLOCATE R_PROJ
         ! ALLOCATE PROJ_PROJ
         ! ALLOCATE G_PROJ
         R_PROJ=GDES%MPLWV_ROW*GDES%NPRO_COL
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, R_PROJ )
         CALL MEM_PROJ_PROJ(PEAK_MEM, MEM_TOTAL, GDES)
         CALL MEM_G_PROJ(PEAK_MEM, MEM_TOTAL, GDES)


         !DEALLOCATE_PROJ_PROJ
         !DEALLOCATE_G_PROJ
         !DEALLOCATE_R_PROJ
         PROJ_PROJ=GDES%NPRO_ROW*GDES%NPRO_COL
         G_PROJ=GDES%NRPLWV_ROW_DATA_POINTS* GDES%NPRO_COL
         CALL DEREG_MEM( MEM_TOTAL, PROJ_PROJ)
         CALL DEREG_MEM( MEM_TOTAL, G_PROJ )
         CALL DEREG_MEM( MEM_TOTAL, R_PROJ )

         CALL DEREG_MEM( MEM_TOTAL, P_G_PROJ )

         CALL MEM_RELEASE_FFT_G( MEM_TOTAL, GDES )
         !FFT_SIGMA contains following calls 
         ! ALLOCATE_R_G 
         !TRANSPOSE_G_R
         !DEALLOCATE_G_R
         !DEALLOCATE_R_G
         P_R_G=GDES%MPLWV_ROW*GDES%NRPLWV_COL_DATA_POINTS 
         P_G_R=GDES%NRPLWV_ROW_DATA_POINTS*GDES%MPLWV_COL
         CALL REG_MEM( PEAK_MEM, MEM_TOTAL, P_R_G)
         CALL MEM_TRANSPOSE_G_R( PEAK_MEM, MEM_TOTAL, GDES)
         CALL DEREG_MEM( MEM_TOTAL, P_G_R )
         CALL DEREG_MEM( MEM_TOTAL, P_R_G )

      END SUBROUTINE MEM_SIGMA_TAU

   END SUBROUTINE STORAGE_REQ_GG

!***********************************************************************
! 
!> this version considers only transition in one spin channel
!! which at least for RPA is what one actually wants to do
!! the spacetime GW codes in chi_GG and chi_super use this routine
! 
!***********************************************************************
   SUBROUTINE DETERMINE_BAND_GAP_SPIN(WDES, W, NB, E1, E2, NOMEGA)
      USE prec
      USE wave
      USE base 
      USE crpa 
      IMPLICIT NONE
      TYPE (wavedes)      :: WDES
      TYPE (wavespin)     :: W
      REAL(q)             :: E1,E2                ! E1 ... band gap
      INTEGER             :: NOMEGA               ! 
      INTEGER             :: NB                   ! maximum band to be considered
      !local
      INTEGER NK, N, ISP, NKP, NP
      REAL(q) :: EMAXO, EMINU, ENEW, EFERMIN
     
      PROFILING_START( 'determine_band_gap_spin' )
  
      ! first search for minimum transition energy in each spin channel
      E1=1000
      EFERMIN=1000
      DO ISP=1,WDES%ISPIN
         EMAXO=-1000
         EMINU= 1000
         DO NK=1,WDES%NKPTS
            DO N=1,NB
              ! 
              ! optionally,  neglect crpa bands 
              IF ( LCRPA .AND. ALLOCATED( LCRPA_BAND ) ) THEN
                 IF ( LCRPA_BAND( N ) ) CYCLE 
              ENDIF
              IF ( ABS(W%FERTOT(N,NK,ISP))>0.5 ) THEN   
                 ! largest occupied energy state
                 EMAXO=MAX(EMAXO,REAL(W%CELTOT(N,NK,ISP),q))
              ELSE
                 ! smallest unoccupied energy state
                 EMINU=MIN(EMINU,REAL(W%CELTOT(N,NK,ISP),q))                
              ENDIF
            ENDDO
         ENDDO
  !       E1=MIN(EMINU-EMAXO,E1)
  !       !for GW calculations we need the minimum and maximum of |e_nk-\mu| 
  !       IF ( LGW ) THEN
  !          ENEW=MIN(ABS(EMAXO-W%EFERMI(ISP)),ABS(EMINU-W%EFERMI(ISP)))
  !          EFERMIN=MIN(ENEW-0.0000001_q,EFERMIN)
  !          ! negative energies must not occur
  !          IF ( EFERMIN < 0 ) THEN
  !             IF ( WDES%COMM%NODE_ME == 1 ) WRITE(*,*) 'ERROR: negative energies found',EFERMIN
  !             STOP
  !        ENDIF
  !       ENDIF
      ENDDO
  
      ! for RPA calculations we need the maximum and minimum of |e_nk - e_n'k'|
      ! now more carefully scan considering all allowed transitions
      E1=1000
      DO ISP=1,WDES%ISPIN
         DO NKP=1,WDES%NKPTS
         DO NK=NKP,WDES%NKPTS
            DO NP=1,NB
            DO N=1,NB
              ! 
              ! optionally,  neglect crpa bands 
              IF ( LCRPA .AND. ALLOCATED( LCRPA_BAND ) ) THEN
                 IF ( LCRPA_BAND( N ) .AND. LCRPA_BAND( NP ) ) CYCLE 
              ENDIF
              IF ( ABS(W%FERTOT(N,NK,ISP)-W%FERTOT(NP,NKP,ISP))>0.02)  THEN
                 ! this is geared towards the new chi_GG routine, where transitions
                 ! are limited to occur from the Fermi-level to occupied or unoccupied states
                 !
                 ! for insulators the first line always applies
                 ! for metals transitions between occupied-occupied or unoccupied-unoccupied
                 ! occur but then those are limited by the other two lines
                 ENEW=MAX(ABS(REAL(W%CELTOT(N,NK,ISP),q)-REAL(W%CELTOT(NP,NKP,ISP),q)), & 
                       MIN(ABS(REAL(W%CELTOT(N,NK,ISP),q)-W%EFERMI(ISP)), &
                           ABS(REAL(W%CELTOT(NP,NKP,ISP),q)-W%EFERMI(ISP))))
                 E1=MIN(E1,ENEW)
              ENDIF
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO
        
      ! now take the minimum of the so determined E1 and previous EFERIMIN
      E1=MIN(E1,EFERMIN)
  
      IF ( LGW ) E1 = E1/2
  
      PROFILING_STOP( 'determine_band_gap_spin' )
   END SUBROUTINE DETERMINE_BAND_GAP_SPIN

!***********************************************************************
! 
!> this routine writes the total energy and grand potential terms
!! for self and non-self consistent Green's function methods
! 
!***********************************************************************

   SUBROUTINE DUMP_DIAGRAMMATIC_ENERGIES(STRING, E, TOTEN, INFO, IO, ECHEMN ) 
      USE base 
      USE ini 
      CHARACTER(LEN=*) :: STRING
      TYPE( energy )    :: E
      REAL(q) :: TOTEN
      TYPE (info_struct) INFO
      TYPE( in_struct ) :: IO 
      COMPLEX(q)        :: ECHEMN
      ! local 
      CHARACTER(LEN=6) :: DIGS
      REAL(q)    TOTEGM

      ! maybe more output is required
      IF ( IO%NWRITE<3 ) THEN
         DIGS='F18.8 '
      ELSE 
         DIGS='F20.10'
      ENDIF
      !
      ! note TOTEN is set in calling routine 
      !
      ! determine total energies from diagrammatic approximations to grand potential 
      ! this works only if LFINITE_TEMPERATURE = .TRUE. 

      ! add correlation term to Galitskii-Migdal internal energy 
      TOTEGM = E%ETRGHF + E%ECGWGM/2

      ! compute grand potential according to Luttinger-Ward in GW[0] approximation
      !
      ! first term is the mean-field part of the grand potential, that is
      ! expressed as 
      ! Omega_0 = { -Tr Ln{-G_0^-1} - 1/2 E_(mean-field) } 
      ! 
      ! in our case E_(mean-field) = E_Hartree + E_exchange 
      ! 
      ! the free-Energy and internal energy is then 
      ! 
      ! F_HF = Omega_0 + mu * N  ->  E_HF = F_HF + TS
      ! 
      ! ECHEMN = - mu * N 
      TOTEN  = E%PSCENC+E%TEWEN+INFO%EALLAT+E%XCENC+E%PAWPS+E%PAWAE+E%ELOGG0 - ECHEMN
      !
      ! Omega_c = -1/2 E_h[G] - E_x[G] + E_RPA_c[G] - Tr{ S_c * G + Ln[ 1-G_HF * S_c ] } 
      !
      ! E%ERPA_INF is the basis set extrapolated value, 
      ! use cutoff value instead, i.e E%ERPA_CUT
      !
      ! For GW it must hold: Phi_c^GW = -1/2 Ln(1-GGV) 
      ! because d Phi^GW / d G = G (1/(1-GGV)) V  = GW  
      ! This means, E%ERPA_CUT = Phi_c^GW in our case 
      IF ( LGREENHF ) THEN
         E%ERPA_CUT = 0 
         ! for non-interacting GF one can calculate the entropy excplicitly
         ! this can be subtracted from F(LW) to compare to GM internal energy
         TOTEN = TOTEN - E%EENTROPY
      ENDIF

      ! add double countig , GM and phi functional terms to Luttinger-Ward free energy
      TOTEN = TOTEN + E%DENC + E%EXHF + E%ERPA_CUT - (E%ECGWGM + E%ELOG1G)
      
      ! compute the entropy from the difference of the free and internal energy
      IF ( .NOT. LGREENHF ) E%EENTROPY  = TOTEN-TOTEGM

      IF ( IO%IU6 >= 0 ) THEN
         WRITE(IO%IU6,"(/, ' ' , A, / &
         &        ' Hartree-Fock free energy of the ion-electron system (eV)'/ &
         &        '  -----------------------------------------------------'/ &
         &        '  alpha Z             PSCENC = ', "//DIGS//"/ &
         &        '  Ewald energy        TEWEN  = ', "//DIGS//"/ &
         &        '  -Hartree energy     DENC   = ', "//DIGS//"/ &
         &        '  -exchange           EXHF   = ', "//DIGS//"/ &
         &        '  -V(xc)+E(xc)        XCENC  = ', "//DIGS//"/ &
         &        '  PAW double counting        = ',2"//DIGS//"/ &
         &        '  eigenvalues         EBANDS = ', "//DIGS//"/ &
         &        '  atomic energy       EATOM  = ', "//DIGS//"/ &
         &        '  -----------------------------------------------------'/ &
         &        '  HF-free energy      FHF    = ',"//DIGS//"' eV' /&
         &        '  ACFDT corr.    EXHF_ACFDT  = ',"//DIGS//"' eV (see jH, gK, PRB 81, 115126)'//&
         &        ' diagrammatic contributions to internal and free energy (eV):'/ &
         &        '  -----------------------------------------------------'/ &
         &        '  HF-Hamiltonian      ETRGHF = ', "//DIGS//"/ &
         &        '  chem. potential     ECHEMN = ',2"//DIGS//"/ &
         &        '  IP grand-potential  ELOGG0 = ', "//DIGS//"/ &
         &        '  Phi-functional GW   EPHIGW = ', "//DIGS//"/ &
         &        '  Galitskii-Migdal    ECGWGM = ', "//DIGS//"/ &
         &        '  Luttinger-Ward      ELOG1G = ', "//DIGS//"/ &
         &        '  Entropy (LW-GM)     EENTRO = ', "//DIGS//"/ &
         &        '  -----------------------------------------------------' &
         )") STRING, &
           E%PSCENC,E%TEWEN,E%DENC,E%EXHF,E%XCENC,E%PAWPS,E%PAWAE, E%EBANDSTR,INFO%EALLAT, &
           E%EBANDSTR+E%DENC+E%XCENC+E%TEWEN+E%PSCENC+E%PAWPS+E%PAWAE+INFO%EALLAT+E%EXHF,E%EXHF_ACFDT, & 
           E%ETRGHF, ECHEMN, E%ELOGG0,E%ERPA_CUT, E%ECGWGM, E%ELOG1G, E%EENTROPY
         WRITE(IO%IU6,"( '  Galitskii-Migdal internal energy TOTEGM =',"//DIGS//"' eV')")TOTEGM 
         IF( LFINITE_TEMPERATURE ) THEN
            IF ( LGREENHF ) THEN
               IF (IO%IU6>=0) WRITE(IO%IU6,"( '  Luttinger-Ward internal energy    TOTEN =',"//DIGS//"' eV')") TOTEN
            ELSE
               IF (IO%IU6>=0) WRITE(IO%IU6,"( '  Luttinger-Ward free energy        TOTEN =',"//DIGS//"' eV')") TOTEN
            ENDIF
         ENDIF
      ENDIF

   END SUBROUTINE DUMP_DIAGRAMMATIC_ENERGIES 

!***********************************************************************
! 
!> this routine writes the total energy and grand potential terms
!! for RPA performed on top of a mean field energy calculations
! 
!***********************************************************************

   SUBROUTINE DUMP_DIAGRAMMATIC_ENERGIES_RPA(STRING, E, TOTEN, INFO, IO ) 
      USE base 
      USE ini 
      CHARACTER(LEN=*) :: STRING
      TYPE( energy )    :: E
      REAL(q) :: TOTEN
      TYPE (info_struct) INFO
      TYPE( in_struct ) :: IO 
      CHARACTER(LEN=6) :: DIGS

      ! RPA energy doesn't contain ACFDT correction
      ! after studying the free-energy of Al as function of T
      ! I found that F(T) becomes smooth for dense k-point sampling
      TOTEN=E%EBANDSTR+E%DENC+E%XCENC+E%TEWEN+E%PSCENC+E%EENTROPY+E%PAWPS+E%PAWAE+INFO%EALLAT+E%EXHF+E%ERPA_CUT

      IF ( IO%NWRITE<3 ) THEN
         DIGS='F18.8 '
      ELSE 
         DIGS='F20.10'
      ENDIF
      IF ( IO%IU6 >= 0 ) THEN
         WRITE(IO%IU6,"(/, ' ' , A, / &
              ' Hartree-Fock free energy of the ion-electron system (eV)'/ &
      &        '  -----------------------------------------------------'/ &
      &        '  alpha Z             PSCENC = ', "//DIGS//"/ &
      &        '  Ewald energy        TEWEN  = ', "//DIGS//"/ &
      &        '  -Hartree energy     DENC   = ', "//DIGS//"/ &
      &        '  -exchange           EXHF   = ', "//DIGS//"/ &
      &        '  -V(xc)+E(xc)        XCENC  = ', "//DIGS//"/ &
      &        '  PAW double counting        = ',2"//DIGS//"/ &
      &        '  entropy T*S         EENTRO = ', "//DIGS//"/ &
      &        '  eigenvalues         EBANDS = ', "//DIGS//"/ &
      &        '  atomic energy       EATOM  = ', "//DIGS//"/ &
      &        '  -----------------------------------------------------'/ &
      &        '  HF-free energy      FHF    = ', "//DIGS//",' eV'/&
      &        '  HF+RPA corr. energy TOTEN  = ', "//DIGS//",' eV'/  & 
      &        '  HF+E_corr(extrapolated)    = ', "//DIGS//",' eV'/  & 
      &        '  ACFDT corr.    EXHF_ACFDT  = ', "//DIGS//",' eV (see jH, gK, PRB 81, 115126)'/ &
      )") STRING, &
              E%PSCENC,E%TEWEN,E%DENC,E%EXHF,E%XCENC,E%PAWPS,E%PAWAE, E%EENTROPY, E%EBANDSTR,INFO%EALLAT, &
              E%EBANDSTR+E%DENC+E%XCENC+E%TEWEN+E%PSCENC+E%EENTROPY+E%PAWPS+E%PAWAE+INFO%EALLAT+E%EXHF,&
              TOTEN,TOTEN+E%ERPA_INF-E%ERPA_CUT,E%EXHF_ACFDT
      ENDIF


  END SUBROUTINE DUMP_DIAGRAMMATIC_ENERGIES_RPA
      
!***********************************************************************
! dump green function w.r.t. reciprocal lattice vectors 
! used only for debugging purposes
!***********************************************************************

   SUBROUTINE DUMPX_FREQ(IU, CHI, GDES, COMMBET, COMMIN, NODE_ME, OMEGA_LOCAL)
      IMPLICIT NONE
      INTEGER                :: IU  !write unit 
      TYPE(responsefunction) :: CHI
      TYPE (greensfdes)      :: GDES
      INTEGER                :: ROW_DATA_POINTS, COL_DATA_POINTS
      TYPE(communic)         :: COMMBET,COMMIN
      INTEGER                :: NODE_ME
      REAL(q)                :: OMEGA_LOCAL(:)
      !local
      INTEGER, PARAMETER     :: NCOLMAX=16
      INTEGER, PARAMETER     :: NROWMAX=16
      INTEGER                :: NPT_SIM
      INTEGER                :: NROWS, NCOLS
      INTEGER                :: IPT
      INTEGER                :: IUNIT
      INTEGER                :: IR

      ROW_DATA_POINTS=GDES%NRPLWV_ROW_DATA_POINTS
      COL_DATA_POINTS=GDES%NRPLWV_COL_DATA_POINTS

      !determine max.  # of columns and rows
      NCOLS=MIN(NCOLMAX,COL_DATA_POINTS)
      NROWS=MIN(NROWMAX,ROW_DATA_POINTS)

      !this is the number of simultaneous points treated in CHI
#ifdef gammareal
      NPT_SIM=SIZE(CHI%RESPONSER,3)
#else
      NPT_SIM=SIZE(CHI%RESPONSEFUN,3)
#endif
      !print the response function 
      DO IPT = 1, NPT_SIM
         IUNIT= IU + 100*COMMBET%NODE_ME + 10*COMMIN%NODE_ME + IPT
!         IUNIT= IU 

!#ifdef verbose
         WRITE(IUNIT,*)'HEAD',OMEGA_LOCAL(IPT)
         DO IR = 1, 3
            WRITE(IUNIT,'(3F12.6, "iamg:",3F12.6)')REAL(CHI%HEAD(IR,1:3,IPT),q),AIMAG(CHI%HEAD(IR,1:3,IPT)) 
         ENDDO
         
         WRITE(IUNIT,*)'WINGS',OMEGA_LOCAL(IPT)
         DO IR = 1, NROWS
            WRITE(IUNIT,'(12F12.6)')-CHI%WING(IR,1:3,IPT)
         ENDDO
         WRITE(IUNIT,*)'CWINGS',OMEGA_LOCAL(IPT)
         DO IR = 1, NROWS
            WRITE(IUNIT,'(12F12.6)')CHI%CWING(IR,1:3,IPT)
         ENDDO
!#endif
!         WRITE(IUNIT,'(A,I4)')' I am the global node:',NODE_ME
         WRITE(IUNIT,'(A,F12.6)')'Freqpt:',OMEGA_LOCAL(IPT)
#ifdef gammareal
         DO IR = 1, NROWS
            WRITE(IUNIT,'(16F12.6)')-CHI%RESPONSER(IR,1:NCOLS,IPT)
         ENDDO
#else
         DO IR = 1, NROWS
            WRITE(IUNIT,'(16F12.6)')-REAL(CHI%RESPONSEFUN(IR,1:NCOLS,IPT),q)
         ENDDO
         WRITE(IUNIT,'(" Imaginary part:")')
         DO IR = 1, NROWS
            WRITE(IUNIT,'(16F12.6)')-IMAG(CHI%RESPONSEFUN(IR,1:NCOLS,IPT))
         ENDDO
#endif
      ENDDO
   END SUBROUTINE DUMPX_FREQ

! serial Fourier transformation test 
   SUBROUTINE FT_TEST_GF(W, GRIDS, IO )
      USE chi_base , ONLY :GOCC_WEIGHT_TIME, GUNO_WEIGHT_TIME
      USE base, ONLY: TOREAL
      TYPE (wavespin)          :: W
      TYPE( imag_grid_handle ) :: GRIDS
      TYPE( in_struct )        :: IO 
      ! local 
      INTEGER, PARAMETER       :: NSTATES = 10 
      GDEF,ALLOCATABLE         :: G_POS(:,:,:,:)    !diagonal positive time Green's function
      GDEF,ALLOCATABLE         :: G_NEG(:,:,:,:)    !diagonal negative time Green's function
      COMPLEX(q),ALLOCATABLE   :: G_FREQ(:,:,:,:)    !diagonal Green's function in freq
      GDEF,ALLOCATABLE         :: U(:,:,:,:)    !diagonal positive time Green's function
      GDEF,ALLOCATABLE         :: V(:,:,:,:)    !diagonal negative time Green's function

      ALLOCATE( G_POS( GRIDS%NOMEGA, W%WDES%NB_TOT, W%WDES%NKPTS, W%WDES%ISPIN ) )
      ALLOCATE( G_NEG( GRIDS%NOMEGA, W%WDES%NB_TOT, W%WDES%NKPTS, W%WDES%ISPIN ) )
      ALLOCATE( G_FREQ( GRIDS%NOMEGA, W%WDES%NB_TOT, W%WDES%NKPTS, W%WDES%ISPIN ) )
      ALLOCATE( U( GRIDS%NOMEGA, W%WDES%NB_TOT, W%WDES%NKPTS, W%WDES%ISPIN ) )
      ALLOCATE( V( GRIDS%NOMEGA, W%WDES%NB_TOT, W%WDES%NKPTS, W%WDES%ISPIN ) )
      G_POS = (0._q,0._q)
      G_NEG = (0._q,0._q)
      U = (0._q,0._q)
      V = (0._q,0._q)
      G_FREQ = (0._q,0._q)

      CALL SET_UP_G_TIME( W, GRIDS, G_POS, G_NEG, IO )

      CALL FORWARD_TRANSFORM( W, GRIDS, G_POS, G_NEG, G_FREQ, IO )

!      CALL FORWARD_TRANSFORM2( W, GRIDS, G_POS, G_NEG, U, V, IO )

      CALL BACKWARD_TRANSFORM( W, GRIDS, G_POS, G_NEG, G_FREQ, IO )

      ! perform some tests 
      CALL INT_G( W, G_POS, G_NEG, G_FREQ, GRIDS, IO )

      CONTAINS 

      SUBROUTINE SET_UP_G_TIME( W, GRIDS, G_POS, G_NEG, IO ) 
         USE chi_base , ONLY :GOCC_WEIGHT_TIME, GUNO_WEIGHT_TIME
         TYPE (wavespin)          :: W
         TYPE( imag_grid_handle ) :: GRIDS
         GDEF                     :: G_POS(:,:,:,:) 
         GDEF                     :: G_NEG(:,:,:,:) 
         TYPE( in_struct )        :: IO 
         !local  
         INTEGER                  :: ISP, NK1, NB, I
         REAL(q)                  :: EN, BETA, TAU, FE
 
         IF ( IO%IU0>=0 ) OPEN(99,FILE='G_TIME.dat',STATUS='REPLACE')
         BETA = TOREAL(GRIDS%TIME%SCALING)
         DO ISP = 1, W%WDES%ISPIN
            IF( IO%IU0>=0 ) WRITE(99, 1)ISP, BETA, W%EFERMI(ISP)
            DO NK1 = 1, W%WDES%NKPTS
               IF( IO%IU0>=0 ) WRITE(99, 2 )NK1,W%WDES%VKPT(:,NK1)
               DO I = 1, GRIDS%NOMEGA
                  TAU=GRIDS%TAU(I)
                  IF( IO%IU0>=0 ) WRITE(99, 3 )I,TAU,GRIDS%TAU_WEIGHT( I )
                  DO NB = 1, W%WDES%NB_TOT
                     EN = REAL(W%CELTOT(NB,NK1,ISP)-W%EFERMI(ISP),q)
                     FE = REAL(W%FERTOT(NB,NK1,ISP),q)
                     G_POS( I, NB, NK1, ISP ) = GUNO_WEIGHT_TIME( EN, FE, BETA, TAU )!+CMPLX(0._q,0.6_q/EN)
                     G_NEG( I, NB, NK1, ISP ) = GOCC_WEIGHT_TIME( EN, FE, BETA,-TAU )!+CMPLX(0._q,1._q/EN)
                     IF( IO%IU0>=0 .AND. NB <= NSTATES ) &
                     WRITE( 99, 4 ) NB, G_NEG( I,NB,NK1, ISP ),G_POS( I,NB,NK1, ISP )
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         IF ( IO%IU0>=0 ) CLOSE(99)
           
1        FORMAT( '# isp, beta=, e_fermi', I3, 2F12.6 )  
2        FORMAT( '# kpt, coord=', I4, 3F12.6 )  
3        FORMAT( '# i, tau, tau_weight=',I3,2F20.10 )
4        FORMAT( I5, 5F20.10 )       
      END SUBROUTINE  SET_UP_G_TIME

      SUBROUTINE FORWARD_TRANSFORM(  W, GRIDS, G_POS, G_NEG, G_FREQ, IO )
         TYPE (wavespin)          :: W
         TYPE( imag_grid_handle ) :: GRIDS
         GDEF                     :: G_POS(:,:,:,:) 
         GDEF                     :: G_NEG(:,:,:,:) 
         COMPLEX(q)               :: G_FREQ(:,:,:,:) 
         TYPE( in_struct )        :: IO 
         !local  
         INTEGER                  :: ISP, NK1, NB, I, J
         REAL(q)                  :: PHI, PSI, BETA

         G_FREQ = (0._q,0._q)
 
         IF ( IO%IU0>=0 ) OPEN(99,FILE='G_FREQ.dat',STATUS='REPLACE')
         BETA = TOREAL(GRIDS%TIME%SCALING)
         DO ISP = 1, W%WDES%ISPIN
            IF( IO%IU0>=0 ) WRITE(99, 1)ISP, BETA, W%EFERMI
            DO NK1 = 1, W%WDES%NKPTS
               IF( IO%IU0>=0 ) WRITE(99, 2 )NK1,W%WDES%VKPT(:,NK1)
               DO I = 1, GRIDS%NOMEGA
                  IF( IO%IU0>=0 ) WRITE(99, 3 )GRIDS%FER_RE(I),GRIDS%FER_RE_WEIGHT( I )
                  DO NB = 1, W%WDES%NB_TOT
                     DO J = 1, GRIDS%NOMEGA ! time loop 
                        PHI = (G_POS(J, NB, NK1, ISP )+GCONJG(G_NEG(J, NB, NK1, ISP)))/2
                        PSI = (G_POS(J, NB, NK1, ISP )-GCONJG(G_NEG(J, NB, NK1, ISP)))/2

                        G_FREQ( I, NB, NK1, ISP ) = G_FREQ( I, NB, NK1, ISP ) + &
                            PSI*GRIDS%TO_FER_RE( I, J) + & 
                           (0._q,1._q)*PHI*GRIDS%TO_FER_RE_CONJG( I, J )

                     ENDDO
                     IF( IO%IU0>=0 .AND. NB <= NSTATES ) &
                     WRITE( 99, 4 ) NB, G_FREQ( I,NB,NK1, ISP )
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         IF ( IO%IU0>=0 ) CLOSE(99)

1        FORMAT( '# isp, beta=, e_fermi', I3, 2F12.6 )  
2        FORMAT( '# kpt, coord=', I4, 3F12.6 )  
3        FORMAT( '# nu, nu_weight=',2F20.10 )
4        FORMAT( I5, 2F20.10 )       
             
      END SUBROUTINE FORWARD_TRANSFORM

      SUBROUTINE FORWARD_TRANSFORM2(  W, GRIDS, G_POS, G_NEG, U, V, IO )
         TYPE (wavespin)          :: W
         TYPE( imag_grid_handle ) :: GRIDS
         GDEF                     :: G_POS(:,:,:,:) 
         GDEF                     :: G_NEG(:,:,:,:) 
         GDEF                     :: U(:,:,:,:) 
         GDEF                     :: V(:,:,:,:) 
         TYPE( in_struct )        :: IO 
         !local  
         INTEGER                  :: ISP, NK1, NB, I, J
         REAL(q)                  :: PHI, PSI, BETA

         U = (0._q,0._q)
         V = (0._q,0._q)
 
         IF ( IO%IU0>=0 ) OPEN(99,FILE='U_V.dat',STATUS='REPLACE')
         BETA = TOREAL(GRIDS%TIME%SCALING)
         DO ISP = 1, W%WDES%ISPIN
            IF( IO%IU0>=0 ) WRITE(99, 1)ISP, BETA, W%EFERMI
            DO NK1 = 1, W%WDES%NKPTS
               IF( IO%IU0>=0 ) WRITE(99, 2 )NK1,W%WDES%VKPT(:,NK1)
               DO I = 1, GRIDS%NOMEGA
                  IF( IO%IU0>=0 ) WRITE(99, 3 )GRIDS%FER_RE(I),GRIDS%FER_RE_WEIGHT( I )
                  DO NB = 1, W%WDES%NB_TOT
                     DO J = 1, GRIDS%NOMEGA ! time loop 
                        PHI = (G_POS(J, NB, NK1, ISP )+GCONJG(G_NEG(J, NB, NK1, ISP)))/2
                        PSI = (G_POS(J, NB, NK1, ISP )-GCONJG(G_NEG(J, NB, NK1, ISP)))/2

                        U( I, NB, NK1, ISP ) = U( I, NB, NK1, ISP ) + PSI*GRIDS%TO_FER_RE( I, J)
                        V( I, NB, NK1, ISP ) = V( I, NB, NK1, ISP ) + PHI*GRIDS%TO_FER_RE_CONJG( I, J )

                     ENDDO
                     IF( IO%IU0>=0 .AND. NB <= NSTATES ) &
                     WRITE( 99, 4 ) NB, U( I,NB,NK1, ISP ), V( I,NB,NK1,ISP)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         IF ( IO%IU0>=0 ) CLOSE(99)

1        FORMAT( '# isp, beta=, e_fermi', I3, 2F12.6 )  
2        FORMAT( '# kpt, coord=', I4, 3F12.6 )  
3        FORMAT( '# nu, nu_weight=',2F20.10 )
4        FORMAT( I5, 4F20.10 )       
             
      END SUBROUTINE FORWARD_TRANSFORM2

      SUBROUTINE BACKWARD_TRANSFORM(  W, GRIDS, G_POS, G_NEG, G_FREQ, IO )
         USE mathtools, ONLY : INVERT_REAL_MATRIX
         TYPE (wavespin)          :: W
         TYPE( imag_grid_handle ) :: GRIDS
         GDEF                     :: G_POS(:,:,:,:) 
         GDEF                     :: G_NEG(:,:,:,:) 
         COMPLEX(q)               :: G_FREQ(:,:,:,:) 
         TYPE( in_struct )        :: IO 
         !local  
         INTEGER                  :: ISP, NK1, NB, I, J
         REAL(q)                  :: EN, BETA, TAU
         COMPLEX(q)               :: PHI, PSI
         REAL(q)                  :: FTCOS( GRIDS%NOMEGA, GRIDS%NOMEGA )
         REAL(q)                  :: FTSIN( GRIDS%NOMEGA, GRIDS%NOMEGA )

         G_POS = (0._q,0._q)
         G_NEG = (0._q,0._q)

         FTCOS=-GRIDS%TO_FER_RE/2
         FTSIN=-GRIDS%TO_FER_RE_CONJG/2
         CALL INVERT_REAL_MATRIX(FTCOS,-1)
         CALL INVERT_REAL_MATRIX(FTSIN,-1)
 
         IF ( IO%IU0>=0 ) OPEN(99,FILE='G_INV.dat',STATUS='REPLACE')
         BETA = TOREAL(GRIDS%TIME%SCALING)
         DO ISP = 1, W%WDES%ISPIN
            IF( IO%IU0>=0 ) WRITE(99, 1)ISP, BETA, W%EFERMI
            DO NK1 = 1, W%WDES%NKPTS
               IF( IO%IU0>=0 ) WRITE(99, 2 )NK1,W%WDES%VKPT(:,NK1)
               DO I = 1, GRIDS%NOMEGA
                  TAU=GRIDS%TAU(I)
                  IF( IO%IU0>=0 ) WRITE(99, 3 )TAU,GRIDS%TAU_WEIGHT( I )
                  DO NB = 1, W%WDES%NB_TOT
             
                     DO J = 1, GRIDS%NOMEGA ! time loop 
                        PHI = G_FREQ( J, NB, NK1, ISP ) + CONJG( G_FREQ( J, NB, NK1, ISP ))
                        PSI = G_FREQ( J, NB, NK1, ISP ) - CONJG( G_FREQ( J, NB, NK1, ISP ))

                        G_NEG( I, NB, NK1, ISP ) = G_NEG( I, NB, NK1, ISP ) &
                           +PHI*FTCOS( I, J)*0.25_q + (0._q,0.25_q)*PSI*FTSIN( I, J )

                        G_POS( I, NB, NK1, ISP ) = G_POS( I, NB, NK1, ISP ) &
                           -PHI*FTCOS( I, J)*0.25_q + (0._q, 0.25_q)*PSI*FTSIN( I, J )

                     ENDDO
                     IF( IO%IU0>=0 .AND. NB <= NSTATES ) &
                     WRITE( 99, 4 ) NB, G_NEG( I,NB,NK1, ISP ), G_POS( I, NB, NK1, ISP)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         IF ( IO%IU0>=0 ) CLOSE(99)

1        FORMAT( '# isp, beta=, e_fermi', I3, 2F12.6 )  
2        FORMAT( '# kpt, coord=', I4, 3F12.6 )  
3        FORMAT( '# tau, tau_weight=',2F20.10 )
4        FORMAT( I5, 5F20.10 )       
             
      END SUBROUTINE BACKWARD_TRANSFORM

      SUBROUTINE INT_G(  W, G_POS, G_NEG, G_FREQ, GRIDS, IO )
         TYPE (wavespin)          :: W
         GDEF                     :: G_POS(:,:,:,:) 
         GDEF                     :: G_NEG(:,:,:,:) 
         COMPLEX(q)               :: G_FREQ(:,:,:,:) 
         TYPE( imag_grid_handle ) :: GRIDS
         TYPE( in_struct )        :: IO 
         !local  
         INTEGER                  :: ISP, NK1, NB, I
         GDEF                     :: E_PHI=zero
         GDEF                     :: E_PSI=zero
         COMPLEX(q)               :: E_RE=(0._q,0._q)
         COMPLEX(q)               :: E_IM=(0._q,0._q)
         GDEF                     :: E_PHI2=zero
         GDEF                     :: E_PSI2=zero
         COMPLEX(q)               :: E_RE2=(0._q,0._q)
         COMPLEX(q)               :: E_IM2=(0._q,0._q)
         GDEF                     :: PHI, PSI

         DO ISP = 1, W%WDES%ISPIN
            DO NK1 = 1, W%WDES%NKPTS
               DO I = 1, GRIDS%NOMEGA
                  DO NB = 1, W%WDES%NB_TOT

                      PHI = G_POS( I, NB, NK1, ISP ) + G_NEG( I, NB, NK1, ISP )
                      PSI = G_POS( I, NB, NK1, ISP ) - G_NEG( I, NB, NK1, ISP )

                      E_PHI = E_PHI + PHI &
                        *GRIDS%TAU_WEIGHT(I)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN

                      E_PSI = E_PSI + PSI &
                        *GRIDS%TAU_WEIGHT(I)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN

                      E_PHI2 = E_PHI2 + PHI**2 &
                        *GRIDS%TAU_WEIGHT(I)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN

                      E_PSI2 = E_PSI2 + PSI**2 &
                        *GRIDS%TAU_WEIGHT(I)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN

                      E_RE2 = E_RE2 + REAL(G_FREQ( I, NB, NK1, ISP),q)**2&
                        *GRIDS%FER_RE_WEIGHT(I)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN

                      E_IM2 = E_IM2 + AIMAG(G_FREQ( I, NB, NK1, ISP))**2&
                        *GRIDS%FER_RE_WEIGHT(I)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN

                      E_RE = E_RE + REAL(G_FREQ( I, NB, NK1, ISP),q)&
                        *GRIDS%FER_RE_WEIGHT(I)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN

                      E_IM = E_IM + AIMAG(G_FREQ( I, NB, NK1, ISP))&
                        *GRIDS%FER_RE_WEIGHT(I)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN

                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         IF( IO%IU0 >=0 ) THEN
            WRITE( *, '(A,I6)' )     '#M', GRIDS%NOMEGA
            WRITE( *, '(A,2F12.6)' ) 'PHI ', E_PHI
            WRITE( *, '(A,2F12.6)' ) 'PSI ', E_PSI
            WRITE( *, '(A,2F12.6)' ) 'U   ', E_RE
            WRITE( *, '(A,2F12.6)' ) 'V   ', E_IM
            WRITE( *, '(A,2F12.6)' ) 'PHI2', E_PHI2
            WRITE( *, '(A,2F12.6)' ) 'PSI2', E_PSI2
            WRITE( *, '(A,2F12.6)' ) 'U2  ', E_RE2
            WRITE( *, '(A,2F12.6)' ) 'V2  ', E_IM2
         ENDIF
      END SUBROUTINE INT_G

   END SUBROUTINE FT_TEST_GF

! compute Galitskii-Migdal from diagonal of SIGMAU and SIGMAO
   SUBROUTINE GALITSKII_MIGDAL_DIAG(W, GO, GU, GRIDS, IO )
      USE chi_base , ONLY :GOCC_WEIGHT_TIME, GUNO_WEIGHT_TIME
      USE base, ONLY: TOREAL
      USE mathtools, ONLY : INVERT_REAL_MATRIX
      TYPE (wavespin)          :: W
      GDEF                     :: GU(:,:,:,:)
      GDEF                     :: GO(:,:,:,:)
      TYPE( imag_grid_handle ) :: GRIDS
      TYPE( in_struct )        :: IO 
      ! local 
      INTEGER, PARAMETER       :: NSTATES = 10 
      INTEGER                  :: NKPTS_IRZ
      INTEGER                  :: ISP, NK1, NB
      INTEGER                  :: I, J
      INTEGER                  :: I_, J_
      INTEGER                  :: NTAU_ROOT, NTAU
      COMPLEX(q)               :: E_GM0=(0._q,0._q) 
      COMPLEX(q)               :: E_GM1=(0._q,0._q) 
      COMPLEX(q)               :: E_GM2=(0._q,0._q) 
      COMPLEX(q)               :: E_GM3=(0._q,0._q) 
      COMPLEX(q)               :: E_GM4=(0._q,0._q) 
      COMPLEX(q)               :: U_GS, V_GS
      COMPLEX(q)               :: U_I, V_J
      COMPLEX(q)               :: G_OMEGA
      COMPLEX(q)               :: SIGMA1,SIGMA2, SIGMA_OMEGA
      REAL(q)                  :: BETA, TAU
      COMPLEX(q),ALLOCATABLE   :: SIGMA_COS(:,:,:,:)
      COMPLEX(q),ALLOCATABLE   :: SIGMA_SIN(:,:,:,:)
      COMPLEX(q),ALLOCATABLE   :: SIGMA_COS2(:,:,:,:)
      COMPLEX(q),ALLOCATABLE   :: SIGMA_SIN2(:,:,:,:)
      REAL(q)                  :: RE_TO_TAU( GRIDS%NOMEGA, GRIDS%NOMEGA )
      REAL(q)                  :: IM_TO_TAU( GRIDS%NOMEGA, GRIDS%NOMEGA )

      NKPTS_IRZ = SIZE( GU, 3 )

      ALLOCATE( SIGMA_COS( SIZE(GU,1), SIZE(GU,2), SIZE(GU,3),SIZE(GU,4) ) )
      ALLOCATE( SIGMA_SIN( SIZE(GU,1), SIZE(GU,2), SIZE(GU,3),SIZE(GU,4) ) )
      SIGMA_COS=(0._q,0._q)
      SIGMA_SIN=(0._q,0._q)
      ALLOCATE( SIGMA_COS2( SIZE(GU,1), SIZE(GU,2), SIZE(GU,3),SIZE(GU,4) ) )
      ALLOCATE( SIGMA_SIN2( SIZE(GU,1), SIZE(GU,2), SIZE(GU,3),SIZE(GU,4) ) )
      SIGMA_COS2=(0._q,0._q)
      SIGMA_SIN2=(0._q,0._q)
      ! inverse transform to time 
      RE_TO_TAU = GRIDS%TO_FER_RE_CONJG
      CALL INVERT_REAL_MATRIX(RE_TO_TAU,-1)
      IM_TO_TAU = GRIDS%TO_FER_IM
      CALL INVERT_REAL_MATRIX(IM_TO_TAU,-1)
      

      BETA = GRIDS%F%BETA
      E_GM1=(0._q,0._q) 
      E_GM2=(0._q,0._q) 

      DO ISP = 1, W%WDES%ISPIN
         DO NK1 = 1, NKPTS_IRZ
            ! forward Transform
            DO NTAU_ROOT = 1, GRIDS%T%NPOINTS_IN_ROOT_GROUP 
               CALL SET_INDICES( NTAU_ROOT, GRIDS%T ) 
               IF (GRIDS%T%LDO_POINT_LOCAL) THEN
                  NTAU= DETERMINE_NTAU_GLOBAL( GRIDS%T%COMM%NODE_ME, NTAU_ROOT, GRIDS ) 
                  ! Fourier transform the diagonal part of Sigma to imaginary frequency 
                  ! (on those cores we have information)
                  ! SIGMA0_DIAG and SIGMAU_DIAG are in tau, SIGMA_COS and SIGMA_SIN are in nu 
                  CALL FT_G_OR_SIGMA(W, GO(NTAU,:,NK1,ISP), GU(NTAU,:,NK1,ISP),&
                    NTAU_ROOT, GRIDS, SIGMA_COS, SIGMA_SIN, NK1, ISP )
#ifdef OLDVERSION
                  CALL FT_G_OR_SIGMA2(W, GO(NTAU,:,NK1,ISP), GU(NTAU,:,NK1,ISP),&
                    NTAU_ROOT, GRIDS, SIGMA_COS2, SIGMA_SIN2, NK1, ISP )
#else
                  CALL FT_G_OR_SIGMA_TO_FER_RE_AND_IM(W, GO(NTAU,:,NK1,ISP), GU(NTAU,:,NK1,ISP),&
                    NTAU_ROOT, GRIDS, SIGMA_COS2, SIGMA_SIN2, NK1, ISP )
#endif
               ENDIF
            ENDDO

         ENDDO
      ENDDO

      ! communicate
      CALLMPI( M_sum_z( GRIDS%T%COMM_BETWEEN_GROUPS, SIGMA_COS, SIZE(SIGMA_COS) ))
      CALLMPI( M_sum_z( GRIDS%T%COMM_BETWEEN_GROUPS, SIGMA_SIN, SIZE(SIGMA_SIN) ))
      CALLMPI( M_sum_z( GRIDS%T%COMM_BETWEEN_GROUPS, SIGMA_COS2, SIZE(SIGMA_COS2) ))
      CALLMPI( M_sum_z( GRIDS%T%COMM_BETWEEN_GROUPS, SIGMA_SIN2, SIZE(SIGMA_SIN2) ))

      E_GM0 = (0._q,0._q)
      E_GM1 = (0._q,0._q)
      E_GM2 = (0._q,0._q)
      E_GM3 = (0._q,0._q)
      E_GM4 = (0._q,0._q)

      DO ISP = 1, W%WDES%ISPIN
         DO NK1 = 1, NKPTS_IRZ
            DO I = 1, GRIDS%F%NPOINTS
               DO NB = 1, SIZE( SIGMA_COS, 2 )
#ifdef OLDVERSION
!
! if G = u + i v, Sigma = U + i V , with real valued even functions u, U and odd functions v, V
! the following integrates ( G(iw) + G(-iw))*( Sigma(iw) + Sigma(iw)^+)/2 = 4*u*U 
!
                  ! self-energy in frequency domain
                  SIGMA_OMEGA = SIGMA_COS( I, NB, NK1, ISP ) + (0._q,1._q)*SIGMA_SIN( I, NB, NK1, ISP )
                  G_OMEGA = 1._q/( CMPLX(0._q, GRIDS%FER_RE(I),q) - W%CELTOT(NB,NK1,ISP)+W%EFERMI(ISP) )
                  U_I = ( G_OMEGA + CONJG(G_OMEGA) )*( SIGMA_OMEGA + CONJG( SIGMA_OMEGA ) )/4
                  E_GM0 = E_GM0 + U_I * GRIDS%FER_RE_WEIGHT(I)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN

!
! the following integrates ( G(iw) - G(-iw))*( S(iw) - S(iw)^+)/2 = -4*v*V 
! for T=0 this should be done with the sine-grid (the grid for odd basis functions)
!
                  ! self-energy in frequency domain
                  SIGMA_OMEGA = SIGMA_COS2( I, NB, NK1, ISP ) + (0._q,1._q)*SIGMA_SIN2( I, NB, NK1, ISP )
                  G_OMEGA = 1._q/( CMPLX(0._q, GRIDS%FER_IM(I),q)-W%CELTOT(NB,NK1,ISP)+W%EFERMI(ISP) )
                  U_I = ( G_OMEGA - CONJG(G_OMEGA) )*( SIGMA_OMEGA - CONJG( SIGMA_OMEGA ) )/4
                  E_GM1 = E_GM1 + U_I * GRIDS%FER_IM_WEIGHT(I)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN


! however, for self-consistent GW calculations, one has only access to G and S
! on the same grid constructed typically for even functions u 
! two ways to handle this: 
! method 1: simply use the quadrature for even functions U on this term
!           for the T=0 grid this yield terrible convergence 
! method 2: perform inverse transform to time grid and forward transform to odd
!           grid and use quadrature for odd functions
!           is also pretty bad for both T=0 and T>0 grid 
! method 3: use averge
!
!FUN FACT: T>0 grid integrates v*V term with RE weights pretty well 
!          so this is the default


                  ! self-energy in frequency domain
                  SIGMA_OMEGA = SIGMA_COS( I, NB, NK1, ISP ) + (0._q,1._q)*SIGMA_SIN( I, NB, NK1, ISP )
                  G_OMEGA = 1._q/( CMPLX(0._q, GRIDS%FER_RE(I),q)-W%CELTOT(NB,NK1,ISP)+W%EFERMI(ISP) )
                  U_I = ( G_OMEGA - CONJG(G_OMEGA) )*( SIGMA_OMEGA - CONJG( SIGMA_OMEGA ) )/4
                  E_GM2 = E_GM2 + U_I * GRIDS%FER_RE_WEIGHT(I)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN

                  ! transform from time to freq 
                  U_I = (0._q,0._q)
                  DO J_ = 1, GRIDS%F%NPOINTS
                     ! transform odd part to time domain
                     V_J = (0._q,0._q)
                     DO I_ = 1, GRIDS%F%NPOINTS
                        SIGMA_OMEGA = SIGMA_COS( I_, NB, NK1, ISP ) + (0._q,1._q)*SIGMA_SIN( I_, NB, NK1, ISP )
                        G_OMEGA = 1._q/( CMPLX(0._q, GRIDS%FER_RE(I_),q) - W%CELTOT(NB,NK1,ISP)+W%EFERMI(ISP) )
                        V_J = V_J + RE_TO_TAU(J_,I_)*( G_OMEGA - CONJG(G_OMEGA) )*( SIGMA_OMEGA - CONJG( SIGMA_OMEGA ) )/4
                     ENDDO
                     U_I = U_I + GRIDS%TO_FER_IM(I,J_)*V_J
                  ENDDO
                  E_GM3 = E_GM3 + U_I * GRIDS%FER_IM_WEIGHT(I)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN
#else
!
! evaluate as in the case of GAMMA
! perform transformation to RE for evean and IM grid for odd part
! SIGMA_COS2 is on RE grid
! SIGMA_SIN2 is on IM grid
! G*S = G*( SIGMA_COS + i SIGMA_SIN ) 
!     = G( RE ) * SIGMA_COS( RE ) + G(-RE ) * SIGMA_COS( RE )  
!     + i* ( G( IM ) * SIGMA_SIN( IM ) - G(-IM ) * SIGMA_SIN( IM ) ) 
!   
!
                  ! G(RE)* SIGMA_COS( RE ) 
                  SIGMA_OMEGA = SIGMA_COS2( I, NB, NK1, ISP )
                  G_OMEGA = 1._q/( CMPLX(0._q, GRIDS%FER_RE(I),q) - W%CELTOT(NB,NK1,ISP)+W%EFERMI(ISP) )
                  E_GM0 = E_GM0 + G_OMEGA* SIGMA_OMEGA * &
                          GRIDS%FER_RE_WEIGHT(I)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN/2
                  ! G(-RE)* SIGMA_COS( RE ) 
                  SIGMA_OMEGA = SIGMA_COS2( I, NB, NK1, ISP )
                  G_OMEGA = 1._q/( CMPLX(0._q,-GRIDS%FER_RE(I),q) - W%CELTOT(NB,NK1,ISP)+W%EFERMI(ISP) )
                  E_GM0 = E_GM0 + G_OMEGA* SIGMA_OMEGA * &
                          GRIDS%FER_RE_WEIGHT(I)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN/2

                  ! G(IM)* SIGMA_SIN( IM ) 
                  SIGMA_OMEGA = SIGMA_SIN2( I, NB, NK1, ISP )
                  G_OMEGA = (0._q,1._q)/( CMPLX(0._q, GRIDS%FER_IM(I),q) - W%CELTOT(NB,NK1,ISP)+W%EFERMI(ISP) )
                  E_GM0 = E_GM0 + G_OMEGA* SIGMA_OMEGA * &
                          GRIDS%FER_IM_WEIGHT(I)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN/2
                  ! G(-IM)* SIGMA_SIN( IM ) 
                  SIGMA_OMEGA = SIGMA_SIN2( I, NB, NK1, ISP )
                  G_OMEGA = (0._q,1._q)/( CMPLX(0._q,-GRIDS%FER_IM(I),q) - W%CELTOT(NB,NK1,ISP)+W%EFERMI(ISP) )
                  E_GM0 = E_GM0 - G_OMEGA* SIGMA_OMEGA * &
                          GRIDS%FER_IM_WEIGHT(I)*W%WDES%WTKPT(NK1)*W%WDES%RSPIN/2

                  E_GM2=0
#endif

               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
      ! compute Integration Error
      IF ( LFINITE_TEMPERATURE ) THEN
         GRIDS%GM_INT_ERROR = E_GM1-E_GM2 + MAX( GRIDS%FER_RE_ERROR, GRIDS%FER_IM_ERROR)
         ! no output for Matsubara grid 
         IF( IO%NWRITE > 4 ) THEN
            IF ( IO%IU0>=0 ) WRITE(*,10) GRIDS%GM_INT_ERROR
            IF ( IO%IU6>=0 ) WRITE(IO%IU6,21)E_GM0+E_GM2
            IF ( IO%IU0>=0 ) WRITE(*,21)E_GM0+E_GM2
         ENDIF
      ELSE
         GRIDS%GM_INT_ERROR = E_GM1-(E_GM2+E_GM3)/2 + MAX( GRIDS%FER_IM_ERROR, GRIDS%FER_IM_ERROR ) 
         IF( IO%NWRITE > 4 ) THEN
            IF ( IO%IU0>=0 ) WRITE(*,10) GRIDS%GM_INT_ERROR
            IF ( IO%IU6>=0 ) WRITE(IO%IU6,21)E_GM0+(E_GM2+E_GM3)/2
            IF ( IO%IU0>=0 ) WRITE(*,21)E_GM0+(E_GM2+E_GM3)/2
         ENDIF
      ENDIF

21    FORMAT(' DIAGONAL GALITSKII-MIGDAL correlation energy (real and imaginary part): ',2F20.10 )
10    FORMAT(' estimated integration error for Galitskii-Migdal energy: ',2F20.10 )

IF ( IO%IU0>=0 .AND. IO%NWRITE>6) WRITE(*,11)E_GM0
IF ( IO%IU0>=0 .AND. IO%NWRITE>6) WRITE(*,12)E_GM1
IF ( IO%IU0>=0 .AND. IO%NWRITE>6) WRITE(*,13)E_GM2
IF ( IO%IU0>=0 .AND. IO%NWRITE>6) WRITE(*,14)E_GM3
IF ( IO%IU0>=0 .AND. IO%NWRITE>6) WRITE(*,15)(E_GM2+E_GM3)/2
IF ( IO%IU6>=0 .AND. IO%NWRITE>5) WRITE(IO%IU6,11)E_GM0
IF ( IO%IU6>=0 .AND. IO%NWRITE>5) WRITE(IO%IU6,12)E_GM1
IF ( IO%IU6>=0 .AND. IO%NWRITE>5) WRITE(IO%IU6,13)E_GM2
IF ( IO%IU6>=0 .AND. IO%NWRITE>5) WRITE(IO%IU6,14)E_GM3
IF ( IO%IU6>=0 .AND. IO%NWRITE>5) WRITE(IO%IU6,15)(E_GM2+E_GM3)/2
11 FORMAT(' ( G + G^+)*(S + S^+) even grid (has error of fermion grid, no correction needed):',2F20.10 )
12 FORMAT(' ( G - G^+)*(S - S^+) reference for terms below (evaluated on fer_im grid):',2F20.10 )
13 FORMAT(' ( G - G^+)*(S - S^+) evaluated on fer_re grid (used for T>0):   ',2F20.10 )
14 FORMAT(' ( G - G^+)*(S - S^+) fer_re -> tau -> fer_im:                   ',2F20.10 )
15 FORMAT(' ( G - G^+)*(S - S^+) average of last two terms  above (for T=0):',2F20.10 )

      DEALLOCATE( SIGMA_COS, SIGMA_SIN )
      DEALLOCATE( SIGMA_COS2, SIGMA_SIN2 )
      CONTAINS

      SUBROUTINE SET_INDICES(NTAU_ROOT, T)
         USE minimax_struct, ONLY: loop_des
         INTEGER         :: NTAU_ROOT 
         TYPE (loop_des) :: T 

         !if the local number of tau points in the group coincides with the
         !one of the root group, the tau loop variable is the same 
         T%NPOINTSC = NTAU_ROOT
         !other groups may have less tau points locally
         T%LDO_POINT_LOCAL = .TRUE.
         IF ( T%NPOINTSC > T%NPOINTS_IN_GROUP ) THEN
            ! set to last tau point in group to avoid going over array bounds
            T%NPOINTSC = T%NPOINTS_IN_GROUP
            ! but actually little work is done since most expensive calls are bypassed
            T%LDO_POINT_LOCAL = .FALSE.
         ENDIF    
         ! set current imaginary time tau point
         T%POINT_CURRENT = T%POINTS_LOCAL(T%NPOINTSC)

         !determine current smallest TAU value of all groups 
         !this is always the tau value of the first tau group 
         !this weird construction here is since M_min_d is not available
         T%POINTMIN_CURRENT = -T%POINT_CURRENT
         CALLMPI( M_max_d(T%COMM, T%POINTMIN_CURRENT, 1))
         T%POINTMIN_CURRENT = -T%POINTMIN_CURRENT
      END SUBROUTINE SET_INDICES

   END SUBROUTINE GALITSKII_MIGDAL_DIAG


!*******************************************************************
!>subroutine to determine weight of requested Green's function
!>for T>0
!*******************************************************************
  FUNCTION GF_WEIGHT( LOCCUPIED, TAU, BETA, EFERMI, E, F ) 
     USE chi_base, ONLY: GOCC_WEIGHT_TIME, GUNO_WEIGHT_TIME
     REAL(q) :: GF_WEIGHT !function value 
     LOGICAL :: LOCCUPIED 
     REAL(q) :: TAU
     REAL(q) :: BETA
     REAL(q) :: EFERMI
     REAL(q) :: E
     REAL(q) :: F

     ! occupied green's function defined only for negative times
     IF ( LOCCUPIED ) THEN 
        GF_WEIGHT = GOCC_WEIGHT_TIME( E-EFERMI, F, BETA, -TAU )
     ! occupied green's function defined only for positive times
     ELSE
        GF_WEIGHT = GUNO_WEIGHT_TIME( E-EFERMI, F, BETA,  TAU )
     ENDIF
     
  END FUNCTION GF_WEIGHT

!*********************************************************************
! 
! Performs a fourier transform to imaginary frequency domain of G or 
! Sigma 
!
! on entry: G*_TAU    ... Green's function on imaginary time axis TAU
! on output: G*_OMEGA ... Green's function in imaginary freq axis OMEGA 
!
!*********************************************************************

   SUBROUTINE FT_G_OR_SIGMA(W, GO_TAU, GU_TAU, NTAU_ROOT, IMAG_GRIDS, &
   G_COS, G_SIN, NQ, ISP)
      IMPLICIT NONE
      TYPE (wavespin)         :: W
      GDEF                    :: GO_TAU(:)         !diagonal of occupied GF in bloch domain 
      GDEF                    :: GU_TAU(:)         !diagonal of unoccupied GF in bloch domain 
      INTEGER                 :: NTAU_ROOT         !current tau point
      TYPE(imag_grid_handle)  :: IMAG_GRIDS        !time and frequency grid handle
      COMPLEX(q)              :: G_COS(:,:,:,:)    !diagonal of occupied GF in bloch domain 
      COMPLEX(q)              :: G_SIN(:,:,:,:)    !diagonal of occupied GF in bloch domain 
      INTEGER                 :: NQ                !current q-point
      INTEGER                 :: ISP               !current spin
      !local variables 
      INTEGER                 :: I, N
      INTEGER                 :: NTAU_GLOBAL
      INTEGER                 :: NTAU_TOTAL
      GDEF                    :: G
      REAL(q)                 :: FTCOS(IMAG_GRIDS%NOMEGA,IMAG_GRIDS%NOMEGA)
      REAL(q)                 :: FTSIN(IMAG_GRIDS%NOMEGA,IMAG_GRIDS%NOMEGA)

      PROFILING_START('ft_g_or_sigma')

      !this is the global tau point of the node
      NTAU_GLOBAL = DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, NTAU_ROOT,&
      IMAG_GRIDS ) 

      FTCOS = IMAG_GRIDS%TO_FER_RE
      FTSIN = IMAG_GRIDS%TO_FER_RE_CONJG

      !perform the transformation 
      DO I=1,IMAG_GRIDS%NOMEGA
         DO N=1,SIZE(GO_TAU,1)
            G_COS(I,N,NQ,ISP)=G_COS(I,N,NQ,ISP)+IMAG_GRIDS%TO_FER_RE(I,NTAU_GLOBAL)*&
               (GU_TAU(N)-GO_TAU(N))
            G_SIN(I,N,NQ,ISP)=G_SIN(I,N,NQ,ISP)+IMAG_GRIDS%TO_FER_RE_CONJG(I,NTAU_GLOBAL)*&
               (GU_TAU(N)+GO_TAU(N))
         ENDDO
     ENDDO
     PROFILING_STOP('ft_g_or_sigma')

   END SUBROUTINE FT_G_OR_SIGMA

!*********************************************************************
! transforms to fer_im grid 
!*********************************************************************
   SUBROUTINE FT_G_OR_SIGMA2(W, GO_TAU, GU_TAU, NTAU_ROOT, IMAG_GRIDS, &
   G_COS, G_SIN, NQ, ISP)
      IMPLICIT NONE
      TYPE (wavespin)         :: W
      GDEF                    :: GO_TAU(:)         !diagonal of occupied GF in bloch domain 
      GDEF                    :: GU_TAU(:)         !diagonal of unoccupied GF in bloch domain 
      INTEGER                 :: NTAU_ROOT         !current tau point
      TYPE(imag_grid_handle)  :: IMAG_GRIDS        !time and frequency grid handle
      COMPLEX(q)              :: G_COS(:,:,:,:)    !diagonal of occupied GF in bloch domain 
      COMPLEX(q)              :: G_SIN(:,:,:,:)    !diagonal of occupied GF in bloch domain 
      INTEGER                 :: NQ                !current q-point
      INTEGER                 :: ISP               !current spin
      !local variables 
      INTEGER                 :: I, N
      INTEGER                 :: NTAU_GLOBAL
      INTEGER                 :: NTAU_TOTAL
      GDEF                    :: G
      REAL(q)                 :: FTCOS(IMAG_GRIDS%NOMEGA,IMAG_GRIDS%NOMEGA)
      REAL(q)                 :: FTSIN(IMAG_GRIDS%NOMEGA,IMAG_GRIDS%NOMEGA)

      PROFILING_START('ft_g_or_sigma2')

      !this is the global tau point of the node
      NTAU_GLOBAL = DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, NTAU_ROOT,&
      IMAG_GRIDS ) 


      FTCOS = IMAG_GRIDS%TO_FER_IM_CONJG
      FTSIN = IMAG_GRIDS%TO_FER_IM

      !perform the transformation 
      DO I=1,IMAG_GRIDS%NOMEGA
         DO N=1,SIZE(GO_TAU,1)
            G_COS(I,N,NQ,ISP)=G_COS(I,N,NQ,ISP)+FTCOS(I,NTAU_GLOBAL)*&
               (GU_TAU(N)-GO_TAU(N))
            G_SIN(I,N,NQ,ISP)=G_SIN(I,N,NQ,ISP)+FTSIN(I,NTAU_GLOBAL)*&
               (GU_TAU(N)+GO_TAU(N))
         ENDDO
     ENDDO
     PROFILING_STOP('ft_g_or_sigma2')

   END SUBROUTINE FT_G_OR_SIGMA2

   SUBROUTINE FT_G_OR_SIGMA_TO_FER_RE_AND_IM(W, GO_TAU, GU_TAU, NTAU_ROOT, IMAG_GRIDS, &
   G_COS, G_SIN, NQ, ISP)
      IMPLICIT NONE
      TYPE (wavespin)         :: W
      GDEF                    :: GO_TAU(:)         !diagonal of occupied GF in bloch domain 
      GDEF                    :: GU_TAU(:)         !diagonal of unoccupied GF in bloch domain 
      INTEGER                 :: NTAU_ROOT         !current tau point
      TYPE(imag_grid_handle)  :: IMAG_GRIDS        !time and frequency grid handle
      COMPLEX(q)              :: G_COS(:,:,:,:)    !diagonal of occupied GF in bloch domain 
      COMPLEX(q)              :: G_SIN(:,:,:,:)    !diagonal of occupied GF in bloch domain 
      INTEGER                 :: NQ                !current q-point
      INTEGER                 :: ISP               !current spin
      !local variables 
      INTEGER                 :: I, N
      INTEGER                 :: NTAU_GLOBAL
      INTEGER                 :: NTAU_TOTAL
      GDEF                    :: G
      REAL(q)                 :: FTCOS(IMAG_GRIDS%NOMEGA,IMAG_GRIDS%NOMEGA)
      REAL(q)                 :: FTSIN(IMAG_GRIDS%NOMEGA,IMAG_GRIDS%NOMEGA)

      PROFILING_START('ft_g_or_sigma_to_fer_re_and_im')

      !this is the global tau point of the node
      NTAU_GLOBAL = DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, NTAU_ROOT,&
      IMAG_GRIDS ) 

      FTCOS = IMAG_GRIDS%TO_FER_RE
      FTSIN = IMAG_GRIDS%TO_FER_IM

      !perform the transformation 
      DO I=1,IMAG_GRIDS%NOMEGA
         DO N=1,SIZE(GO_TAU,1)
            G_COS(I,N,NQ,ISP)=G_COS(I,N,NQ,ISP)+FTCOS(I,NTAU_GLOBAL)*&
               (GU_TAU(N)-GO_TAU(N))
            G_SIN(I,N,NQ,ISP)=G_SIN(I,N,NQ,ISP)+FTSIN(I,NTAU_GLOBAL)*&
               (GU_TAU(N)+GO_TAU(N))
         ENDDO
     ENDDO
     PROFILING_STOP('ft_g_or_sigma_to_fer_re_and_im')

   END SUBROUTINE FT_G_OR_SIGMA_TO_FER_RE_AND_IM

!*********************************************************************
! 
! Performs an inverse fourier transform to imaginary time domain of G or 
! Sigma 
!
! on entry : GOMEGA ... Green's function on imaginary frequency axis OMEGA
! on output: G*_TAU ... Green's function in imaginary freq axis TAU 
!
! FT SIGMA from frequency to grid is given by        
!
! occupied part:                                   -1 
! GO_TAU = 1/2( GOMEGA + GOMEGA*)*FTCOS  +
!                                                  -1  
!              i/2( GOMEGA*- GOMEGA )*FTSIN
!
! unoccupied part:                                 -1
! GU_TAU =-1/2( GOMEGA + GOMEGA*)*FTCOS  +
!                                                  -1  
!              i/2( GOMEGA*- GOMEGA )*FTSIN
!
!*********************************************************************

   SUBROUTINE INV_FT_G_OR_SIGMA(W, GOMEGA, GO_TAU, GU_TAU, IMAG_GRIDS, &
     NKPTS_IRZ)
      USE mathtools, ONLY : INVERT_REAL_MATRIX
      IMPLICIT NONE
      TYPE (wavespin)         :: W
      COMPLEX(q)              :: GOMEGA(:,:,:,:)   !diagonal of frequency dependent GF
      GDEF                    :: GO_TAU(:,:,:,:)   !diagonal of occupied GF in bloch domain 
      GDEF                    :: GU_TAU(:,:,:,:)   !diagonal of unoccupied GF in bloch domain 
      TYPE(imag_grid_handle)  :: IMAG_GRIDS        !time and frequency grid handle
      INTEGER                 :: NKPTS_IRZ
      !local variables 
      INTEGER                 :: NQ, ISP 
      INTEGER                 :: I, N, J

      REAL(q),ALLOCATABLE     :: INVFTCOS(:,:)
      REAL(q),ALLOCATABLE     :: INVFTSIN(:,:)
      COMPLEX(q)              :: SIGMA1,SIGMA2

      PROFILING_START('inv_ft_g_or_sigma')

      ALLOCATE(INVFTCOS(IMAG_GRIDS%NOMEGA,IMAG_GRIDS%NOMEGA))
      ALLOCATE(INVFTSIN(IMAG_GRIDS%NOMEGA,IMAG_GRIDS%NOMEGA))
      INVFTCOS=IMAG_GRIDS%TO_FER_RE
      INVFTSIN=IMAG_GRIDS%TO_FER_RE_CONJG
      CALL INVERT_REAL_MATRIX(INVFTCOS,-1)
      CALL INVERT_REAL_MATRIX(INVFTSIN,-1)
      GO_TAU = zero
      GU_TAU = zero

      DO ISP=1,SIZE(GO_TAU,4)
         DO NQ=1,SIZE(GO_TAU,3)
            !perform the transformation 
            DO J=1,IMAG_GRIDS%NOMEGA
               DO I=1,IMAG_GRIDS%NOMEGA
                  DO N=1,SIZE(GO_TAU,2)
                     ! real and imaginary part of selfenergy, factor 2 included in matrices
                     SIGMA1=GOMEGA( I, N, NQ, ISP)+CONJG( GOMEGA( I, N, NQ, ISP) )
                     SIGMA2=GOMEGA( I, N, NQ, ISP)-CONJG( GOMEGA( I, N, NQ, ISP) )

                     GO_TAU(J,N,NQ,ISP)=GO_TAU(J,N,NQ,ISP)+INVFTCOS(J,I)*SIGMA1*0.25_q
                     GO_TAU(J,N,NQ,ISP)=GO_TAU(J,N,NQ,ISP)+INVFTSIN(J,I)*SIGMA2*(0._q,0.25_q)

                     GU_TAU(J,N,NQ,ISP)=GU_TAU(J,N,NQ,ISP)-INVFTCOS(J,I)*SIGMA1*0.25_q
                     GU_TAU(J,N,NQ,ISP)=GU_TAU(J,N,NQ,ISP)+INVFTSIN(J,I)*SIGMA2*(0._q,0.25_q)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO


      DEALLOCATE(INVFTCOS)
      DEALLOCATE(INVFTSIN)

      PROFILING_STOP('inv_ft_g_or_sigma')

   END SUBROUTINE INV_FT_G_OR_SIGMA

!*********************************************************************
! 
! this routine determines the change of the density matrix gamma.
! Specifically, the integral
!
!  delta gamma = 1/ pi int_0^infty G_0(iw) sigma(iw) G_0(iw) dw
!
! on the cosine grid (since this is an even function)
! since the routine receives the sigma as a function
! of imaginary time, the function is first FT to iw
! then the one electron Green functions are multiplied in
! Mean field representation is used for the Green functions
! (i.e. eigenvectors and eigenvalues are needed and read from W)
!
! on input: 
!  SIGMAO_MAT and SIGMAU_MAT are the self-energies for a given tau point
!      and orbital representation distributed by scaLAPACK  
!
! on output:
!  CORR_MAT contains the contribution from the current tau point
!      to the change of the density matrix in orbital representation
! 
!*********************************************************************

   SUBROUTINE ADD_G_SIGMA_G( W, SIGMAO_MAT, SIGMAU_MAT, CORR_MAT, & 
        NTAU, IMAG_GRIDS, NQ, ISP, DESCA )
      USE scala
      USE constant
      IMPLICIT NONE
      TYPE (wavespin)        :: W
      GDEF                   :: SIGMAO_MAT(:) ! matrix <i| Sigma_o(tau) | a> stored distributed
      GDEF                   :: SIGMAU_MAT(:) ! matrix <i| Sigma_u(tau) | a> stored distributed     
      GDEF                   :: CORR_MAT(:)   ! contribution to density matrix from correlation effects
      INTEGER                :: NTAU          !current tau point
      TYPE(imag_grid_handle) :: IMAG_GRIDS    !time and frequency grid handle
      INTEGER                :: NQ            !current q-point
      INTEGER                :: ISP           !current spin
      INTEGER                :: DESCA(*)          
    !local variables 
      REAL(q)                :: EFERMI
      INTEGER                :: I
      INTEGER                :: NTAU_GLOBAL
      COMPLEX(q),ALLOCATABLE :: SIGMA_COS(:)
      COMPLEX(q)             :: CELTOTINV( SIZE(W%CELTOT,1))
      
      PROFILING_START('add_g_sigma_g')

      EFERMI = W%EFERMI(ISP)
      ALLOCATE(SIGMA_COS(SIZE(SIGMAO_MAT,1)))

      !this is the global tau point of the node
      NTAU_GLOBAL = DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, NTAU,&
      IMAG_GRIDS )

      !
      ! only the term Re( G * Sigma * G ) is relevant, so 
      !
      SIGMA_COS=(0._q,0._q)
      DO I=1,IMAG_GRIDS%NOMEGA
         ! take all four terms into account at once
         ! multiply \Sigma_cos from left and right with G and integrate only real part
         ! cosine transformation to i w
         SIGMA_COS(:)=(SIGMAO_MAT(:)-SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE(I,NTAU_GLOBAL) &
                    - (SIGMAO_MAT(:)+SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE_CONJG(I,NTAU_GLOBAL)*(0._q,1.0_q)

         ! for the Kohn-Sham Green's function the expression is almost trivial
         CELTOTINV=1._q/(CMPLX(0._q,IMAG_GRIDS%FER_RE(I),q)-(W%CELTOT(:,NQ,ISP)-EFERMI))
         !now compute Re(G.\Sigma.G)
         CALL MULTIPLY_BY_IVEIGENVALUE(W%WDES%NB_TOTK(NQ,ISP), SIGMA_COS, CELTOTINV, DESCA)
         ! and add to total first order change of density matrix multiplied by weight
         CORR_MAT(:)=CORR_MAT(:)+SIGMA_COS(:)*(IMAG_GRIDS%FER_RE_WEIGHT(I)/2._q)  ! make up for special scaling in OMEGAWEIGHT
         ! take all four terms into account at once
         ! multiply \Sigma_cos from left and right with G and integrate only real part
         ! cosine transformation to i w
         SIGMA_COS(:)=(SIGMAO_MAT(:)-SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE(I,NTAU_GLOBAL) &
                    + (SIGMAO_MAT(:)+SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE_CONJG(I,NTAU_GLOBAL)*(0._q,1.0_q)

         ! for the Kohn-Sham Green's function the expression is almost trivial
         CELTOTINV=1._q/(CMPLX(0._q,-IMAG_GRIDS%FER_IM(I),q)-(W%CELTOT(:,NQ,ISP)-EFERMI))
         !now compute Re(G.\Sigma.G)
         CALL MULTIPLY_BY_IVEIGENVALUE(W%WDES%NB_TOTK(NQ,ISP), SIGMA_COS, CELTOTINV, DESCA)
         ! and add to total first order change of density matrix multiplied by weight
         CORR_MAT(:)=CORR_MAT(:)+SIGMA_COS(:)*(IMAG_GRIDS%FER_IM_WEIGHT(I)/2._q)  ! make up for special scaling in OMEGAWEIGHT
      ENDDO

      DEALLOCATE(SIGMA_COS)

      PROFILING_STOP('add_g_sigma_g')

    END SUBROUTINE ADD_G_SIGMA_G


    SUBROUTINE ADD_G_SIGMA_G_new( W, SIGMAO_MAT, SIGMAU_MAT, CORR_MAT, & 
        NTAU, IMAG_GRIDS, NQ, ISP, DESCA )
      USE scala
      USE constant
      IMPLICIT NONE
      TYPE (wavespin)        :: W
      GDEF                   :: SIGMAO_MAT(:) ! matrix <i| Sigma_o(tau) | a> stored distributed
      GDEF                   :: SIGMAU_MAT(:) ! matrix <i| Sigma_u(tau) | a> stored distributed     
      GDEF                   :: CORR_MAT(:)   ! contribution to density matrix from correlation effects
      INTEGER                :: NTAU          !current tau point
      TYPE(imag_grid_handle) :: IMAG_GRIDS    !time and frequency grid handle
      INTEGER                :: NQ            !current q-point
      INTEGER                :: ISP           !current spin
      INTEGER                :: DESCA(*)          
    !local variables 
      REAL(q)                :: EFERMI
      INTEGER                :: I
      INTEGER                :: NTAU_GLOBAL
      COMPLEX(q),ALLOCATABLE :: SIGMA_w(:)
      COMPLEX(q)             :: CELTOTINV( SIZE(W%CELTOT,1))
      PROFILING_START('add_g_sigma_g_new')

      IF (IMAG_GRIDS%FER_IM_WEIGHT(1)==0) THEN
         ! if the weight of the first point on the sine grid is 0, call old routine
         ! which does not require sine grid
         CALL ADD_G_SIGMA_G( W, SIGMAO_MAT, SIGMAU_MAT, CORR_MAT, & 
              NTAU, IMAG_GRIDS, NQ, ISP, DESCA )
         PROFILING_STOP('add_g_sigma_g_new')
         RETURN
      ENDIF
      EFERMI = W%EFERMI(ISP)
      
      ALLOCATE(SIGMA_w(SIZE(SIGMAO_MAT,1)))

      ! additional factor 1/2: because at T=0 cos transform transforms bosons as too

      !this is the global tau point of the node
      NTAU_GLOBAL = DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, NTAU,&
      IMAG_GRIDS )
      
      !First part only calculates the contribution of SIGMA_even on COS-grid
      SIGMA_w=(0._q,0._q)
      DO I=1,IMAG_GRIDS%NOMEGA
         ! multiply \Sigma_w (even) from left and right with G and perform integration for \omega > 0:
         !------------------------------------------------------------------------
         ! cosine transformation of SIGMA_even to i w
         SIGMA_w(:)=(SIGMAO_MAT(:)-SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE(I,NTAU_GLOBAL)
         !first the even part for \omega > 0 is calculated
         ! for the Kohn-Sham Green's function the expression G(iw) is almost trivial
         CELTOTINV=1._q/(CMPLX(0._q,IMAG_GRIDS%FER_RE(I),q)-(W%CELTOT(:,NQ,ISP)-EFERMI))
         !now compute (G.\Sigma_even.G)
         CALL MULTIPLY_BY_IVEIGENVALUE(W%WDES%NB_TOTK(NQ,ISP), SIGMA_w, CELTOTINV, DESCA)
         ! and calculate 1st part of the first order change of the density matrix 
         ! due to the change of the overlap operator S, by integrating 
         ! (G.\Sigma_even.G)(iw) (i\omega+e_fermi) on the COS-grid from 0 to \infty
         CORR_MAT(:)=CORR_MAT(:)+SIGMA_w(:)*(IMAG_GRIDS%FER_RE_WEIGHT(I)/2._q)
         !-------------------------------------------------------------------------

         ! repeat for \omega < 0
         SIGMA_w(:)=(SIGMAO_MAT(:)-SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE(I,NTAU_GLOBAL)
         CELTOTINV=1._q/(CMPLX(0._q,-IMAG_GRIDS%FER_RE(I),q)-(W%CELTOT(:,NQ,ISP)-EFERMI))
         !now compute (G.\Sigma_even.G)
         CALL MULTIPLY_BY_IVEIGENVALUE(W%WDES%NB_TOTK(NQ,ISP), SIGMA_w, CELTOTINV, DESCA)
         CORR_MAT(:)=CORR_MAT(:)+SIGMA_w(:)*(IMAG_GRIDS%FER_RE_WEIGHT(I)/2._q)
      ENDDO
      !Second part only calculates the contribution of SIGMA_odd on SIN-grid
      SIGMA_w=(0._q,0._q)
      DO I=1,IMAG_GRIDS%NOMEGA
         ! multiply \Sigma_w (odd) from left and right with G and perform integration for \omega > 0:
         !------------------------------------------------------------------------
         ! sine transformation of SIGMA_odd to i w
         SIGMA_w(:)=-(SIGMAO_MAT(:)+SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_IM(I,NTAU_GLOBAL)*(0._q,1.0_q)
         ! for the Kohn-Sham Green's function the expression G(iw) is almost trivial
         CELTOTINV=1._q/(CMPLX(0._q,IMAG_GRIDS%FER_IM(I),q)-(W%CELTOT(:,NQ,ISP)-EFERMI))
         !now compute (G.\Sigma_odd.G)
         CALL MULTIPLY_BY_IVEIGENVALUE(W%WDES%NB_TOTK(NQ,ISP), SIGMA_w, CELTOTINV, DESCA)
         ! and calculate 1st part of the first order change of the density matrix 
         ! due to the change of the overlap operator S, by integrating 
         ! (G.\Sigma_odd.G)(iw) (i\omega+e_fermi) on the COS-grid from 0 to \infty
         CORR_MAT(:)=CORR_MAT(:)+SIGMA_w(:)*(IMAG_GRIDS%FER_IM_WEIGHT(I)/2._q)
         !-------------------------------------------------------------------------

         ! repeat for \omega < 0
         SIGMA_w(:)=+(SIGMAO_MAT(:)+SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_IM(I,NTAU_GLOBAL)*(0._q,1.0_q)
         CELTOTINV=1._q/(CMPLX(0._q,-IMAG_GRIDS%FER_IM(I),q)-(W%CELTOT(:,NQ,ISP)-EFERMI))
         !now compute (G.\Sigma_odd.G)
         CALL MULTIPLY_BY_IVEIGENVALUE(W%WDES%NB_TOTK(NQ,ISP), SIGMA_w, CELTOTINV, DESCA)
         CORR_MAT(:)=CORR_MAT(:)+SIGMA_w(:)*(IMAG_GRIDS%FER_IM_WEIGHT(I)/2._q)
      ENDDO
      DEALLOCATE(SIGMA_w)

      PROFILING_STOP('add_g_sigma_g_new')

    END SUBROUTINE ADD_G_SIGMA_G_new

!*********************************************************************
! 
! this routine determines the change of the density matrix gamma.
! Specifically, the integral
!
!  delta gamma = 1/ (2 pi) int_-infty^infty G_0(iw) sigma(iw) G_0(iw) iw dw
!
! on the cosine grid (since this is an even function)
! since the routine receives the sigma as a function
! of imaginary time, the function is first FT to iw
! then the one electron Green functions are multiplied in
! Mean field representation is used for the Green functions
! (i.e. eigenvectors and eigenvalues are needed and read from W)
!
! on input: 
!  SIGMAO_MAT and SIGMAU_MAT are the self-energies for a given tau point
!      and orbital representation distributed by scaLAPACK  
!
! on output:
!  S_MAT contains the contribution from the current tau point
!      to the change of the density matrix with respect to the change
!      in the overlap S in orbital representation
!
!*********************************************************************

   SUBROUTINE ADD_G_SIGMA_Gxw( W, SIGMAO_MAT, SIGMAU_MAT, S_MAT, & 
        NTAU, IMAG_GRIDS, NQ, ISP, DESCA )
      USE scala
      USE constant
      IMPLICIT NONE
      TYPE (wavespin)        :: W
      GDEF                   :: SIGMAO_MAT(:) ! matrix <i| Sigma_o(tau) | a> stored distributed
      GDEF                   :: SIGMAU_MAT(:) ! matrix <i| Sigma_u(tau) | a> stored distributed     
      GDEF                   :: S_MAT(:)      ! contribution to density matrix from correlation effects
      INTEGER                :: NTAU          !current tau point
      TYPE(imag_grid_handle) :: IMAG_GRIDS    !time and frequency grid handle
      INTEGER                :: NQ            !current q-point
      INTEGER                :: ISP           !current spin
      INTEGER                :: DESCA(*)          
    !local variables 
      REAL(q)                :: EFERMI        ! fermi energy 
      INTEGER                :: I
      INTEGER                :: NTAU_GLOBAL
      COMPLEX(q),ALLOCATABLE :: SIGMA_COS(:)
      COMPLEX(q)             :: CELTOTINV( SIZE(W%CELTOT,1))
      
      PROFILING_START('add_g_sigma_gxw')

      EFERMI = W%EFERMI(ISP)
      ALLOCATE(SIGMA_COS(SIZE(SIGMAO_MAT,1)))

      !this is the global tau point of the node
      NTAU_GLOBAL = DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, NTAU,&
      IMAG_GRIDS )
      !
      ! only the term Re( G * Sigma * G ) is relevant, so 
      !
      SIGMA_COS=(0._q,0._q)
      DO I=1,IMAG_GRIDS%NOMEGA
         ! take all four terms into account at once
         ! multiply \Sigma_cos from left and right with G and integrate only real part
         ! cosine transformation to i w
         SIGMA_COS(:)=(SIGMAO_MAT(:)-SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE(I,NTAU_GLOBAL) &
                    - (SIGMAO_MAT(:)+SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE_CONJG(I,NTAU_GLOBAL)*(0._q,1._q)

         ! for the Kohn-Sham Green's function the expression is almost trivial
         CELTOTINV=1._q/(CMPLX(0._q,IMAG_GRIDS%FER_RE(I),q)-(W%CELTOT(:,NQ,ISP)-EFERMI))
         !now compute Re(G.\Sigma.G)
         CALL MULTIPLY_BY_IVEIGENVALUE(W%WDES%NB_TOTK(NQ,ISP), SIGMA_COS, CELTOTINV, DESCA)
         ! and add to total first order change of density matrix multiplied by weight
         S_MAT(:)=S_MAT(:)+SIGMA_COS(:)*(IMAG_GRIDS%FER_RE_WEIGHT(I)/2._q*(CMPLX(0._q,IMAG_GRIDS%FER_RE(I),q)+EFERMI))
         ! make up for special scaling in OMEGAWEIGHT, and multiply with (i\omega+e_fermi)

         ! take all four terms into account at once
         ! multiply \Sigma_cos from left and right with G and integrate only real part
         ! cosine transformation to i w
         SIGMA_COS(:)=(SIGMAO_MAT(:)-SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE(I,NTAU_GLOBAL) &
                    + (SIGMAO_MAT(:)+SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE_CONJG(I,NTAU_GLOBAL)*(0._q,1.0_q)

         ! for the Kohn-Sham Green's function the expression is almost trivial
         CELTOTINV=1._q/(CMPLX(0._q,-IMAG_GRIDS%FER_RE(I),q)-(W%CELTOT(:,NQ,ISP)-EFERMI))
         !now compute Re(G.\Sigma.G)
         CALL MULTIPLY_BY_IVEIGENVALUE(W%WDES%NB_TOTK(NQ,ISP), SIGMA_COS, CELTOTINV, DESCA)
         ! and add to total first order change of density matrix multiplied by weight
         S_MAT(:)=S_MAT(:)+SIGMA_COS(:)*(IMAG_GRIDS%FER_RE_WEIGHT(I)/2._q*(CMPLX(0.0_q,-IMAG_GRIDS%FER_RE(I),q)+EFERMI))
         ! make up for special scaling in OMEGAWEIGHT, and multiply with \omega
      ENDDO

      DEALLOCATE(SIGMA_COS)

      PROFILING_STOP('add_g_sigma_gxw')

    END SUBROUTINE ADD_G_SIGMA_Gxw

!*********************************************************************
! 
! this routine determines the change of the density matrix gamma.
! Specifically, the integral
!
!  delta gamma = 1/ (2 pi) int_-infty^infty G_0(iw) sigma(iw) G_0(iw) iw dw
!
! on the cosine grid (since this is an even function)
! since the routine receives the sigma as a function
! of imaginary time, the function is first FT to iw
! then the one electron Green functions are multiplied in
! Mean field representation is used for the Green functions
! (i.e. eigenvectors and eigenvalues are needed and read from W)
!
! on input: 
!  SIGMAO_MAT and SIGMAU_MAT are the self-energies for a given tau point
!      and orbital representation distributed by scaLAPACK  
!
! on output:
!  S_MAT contains the contribution from the current tau point
!      to the change of the density matrix with respect to the change
!      in the overlap S in orbital representation
!
!*********************************************************************

   SUBROUTINE ADD_G_SIGMA_Gxw_new( W, SIGMAO_MAT, SIGMAU_MAT, S_MAT, & 
        NTAU, IMAG_GRIDS, NQ, ISP, DESCA )
      USE scala
      USE constant
      IMPLICIT NONE
      TYPE (wavespin)        :: W
      GDEF                   :: SIGMAO_MAT(:) ! matrix <i| Sigma_o(tau) | a> stored distributed
      GDEF                   :: SIGMAU_MAT(:) ! matrix <i| Sigma_u(tau) | a> stored distributed     
      GDEF                   :: S_MAT(:)      ! contribution to density matrix from correlation effects
      INTEGER                :: NTAU          !current tau point
      TYPE(imag_grid_handle) :: IMAG_GRIDS    !time and frequency grid handle
      INTEGER                :: NQ            !current q-point
      INTEGER                :: ISP           !current spin
      INTEGER                :: DESCA(*)          
    !local variables 
      REAL(q)                :: EFERMI        ! fermi energy 
      INTEGER                :: I
      INTEGER                :: NTAU_GLOBAL
      COMPLEX(q),ALLOCATABLE :: SIGMA_w(:)
      COMPLEX(q)             :: CELTOTINV( SIZE(W%CELTOT,1))

      PROFILING_START('add_g_sigma_gxw_new')

      IF (IMAG_GRIDS%FER_IM_WEIGHT(1)==0) THEN
         ! if the weight of the first point on the sin grid is 0, call old routine
         CALL ADD_G_SIGMA_Gxw( W, SIGMAO_MAT, SIGMAU_MAT, S_MAT, & 
              NTAU, IMAG_GRIDS, NQ, ISP, DESCA )
         PROFILING_STOP('add_g_sigma_gxw_new')
         RETURN
      ENDIF
      EFERMI = W%EFERMI(ISP)
      
      ALLOCATE(SIGMA_w(SIZE(SIGMAO_MAT,1)))

      !this is the global tau point of the node
      NTAU_GLOBAL = DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, NTAU,&
      IMAG_GRIDS )
      
      !First part only calculates the contribution of SIGMA_even on COS-grid
      SIGMA_w=(0._q,0._q)
      DO I=1,IMAG_GRIDS%NOMEGA
         ! multiply \Sigma_w (even) from left and right with G and perform integration for \omega > 0:
         !------------------------------------------------------------------------
         ! cosine transformation of SIGMA_even to i w
         SIGMA_w(:)=(SIGMAO_MAT(:)-SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE(I,NTAU_GLOBAL)
         !first the even part for \omega > 0 is calculated
         ! for the Kohn-Sham Green's function the expression G(iw) is almost trivial
         CELTOTINV=1._q/(CMPLX(0._q,IMAG_GRIDS%FER_RE(I),q)-(W%CELTOT(:,NQ,ISP)-EFERMI))
         !now compute (G.\Sigma_even.G)(iw)
         CALL MULTIPLY_BY_IVEIGENVALUE(W%WDES%NB_TOTK(NQ,ISP), SIGMA_w, CELTOTINV, DESCA)
         ! and calculate 1st part of the first order change of the density matrix 
         ! due to the change of the overlap operator S, by integrating 
         ! (G.\Sigma_even.G)(iw) (i\omega+e_fermi) on the COS-grid from 0 to \infty
         S_MAT(:)=S_MAT(:)+SIGMA_w(:)*(IMAG_GRIDS%FER_RE_WEIGHT(I)/2._q*(CMPLX(0._q,IMAG_GRIDS%FER_RE(I),q)+EFERMI))
         !-------------------------------------------------------------------------

         ! repeat for \omega < 0
         SIGMA_w(:)=(SIGMAO_MAT(:)-SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE(I,NTAU_GLOBAL)
         CELTOTINV=1._q/(CMPLX(0._q,-IMAG_GRIDS%FER_RE(I),q)-(W%CELTOT(:,NQ,ISP)-EFERMI))
         !now compute (G.\Sigma_even.G)(-iw)
         CALL MULTIPLY_BY_IVEIGENVALUE(W%WDES%NB_TOTK(NQ,ISP), SIGMA_w, CELTOTINV, DESCA)
         S_MAT(:)=S_MAT(:)+SIGMA_w(:)*(IMAG_GRIDS%FER_RE_WEIGHT(I)/2._q*(CMPLX(0._q,-IMAG_GRIDS%FER_RE(I),q)+EFERMI))

      ENDDO

      !Second part only calculates the contribution of SIGMA_odd on SIN-grid
      SIGMA_w=(0._q,0._q)
      DO I=1,IMAG_GRIDS%NOMEGA
         ! multiply \Sigma_w (odd) from left and right with G and perform integration for \omega > 0:
         !------------------------------------------------------------------------
         ! sine transformation of SIGMA_odd to i w
         SIGMA_w(:)=-(SIGMAO_MAT(:)+SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_IM(I,NTAU_GLOBAL)*(0._q,1._q)
         ! for the Kohn-Sham Green's function the expression G(iw) is almost trivial
         CELTOTINV=1._q/(CMPLX(0._q,IMAG_GRIDS%FER_IM(I),q)-(W%CELTOT(:,NQ,ISP)-EFERMI))
         !now compute (G.\Sigma_odd.G)(iw)
         CALL MULTIPLY_BY_IVEIGENVALUE(W%WDES%NB_TOTK(NQ,ISP), SIGMA_w, CELTOTINV, DESCA)
         ! and calculate 1st part of the first order change of the density matrix 
         ! due to the change of the overlap operator S, by integrating 
         ! (G.\Sigma_odd.G)(iw) (i\omega+e_fermi) on the COS-grid from 0 to \infty
         S_MAT(:)=S_MAT(:)+SIGMA_w(:)*(IMAG_GRIDS%FER_IM_WEIGHT(I)/2._q*(CMPLX(0._q,IMAG_GRIDS%FER_IM(I),q)+EFERMI))
         !-------------------------------------------------------------------------

         ! repeat for \omega < 0
         SIGMA_w(:)= (SIGMAO_MAT(:)+SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_IM(I,NTAU_GLOBAL)*(0._q,1._q)
         CELTOTINV=1._q/(CMPLX(0._q,-IMAG_GRIDS%FER_IM(I),q)-(W%CELTOT(:,NQ,ISP)-EFERMI))
         !now compute (G.\Sigma_even.G)(-iw)
         CALL MULTIPLY_BY_IVEIGENVALUE(W%WDES%NB_TOTK(NQ,ISP), SIGMA_w, CELTOTINV, DESCA)
         S_MAT(:)=S_MAT(:)+SIGMA_w(:)*(IMAG_GRIDS%FER_IM_WEIGHT(I)/2._q*(CMPLX(0._q,-IMAG_GRIDS%FER_IM(I),q)+EFERMI))
         
      ENDDO

      DEALLOCATE(SIGMA_w)

      PROFILING_STOP('add_g_sigma_gxw_new')

    END SUBROUTINE ADD_G_SIGMA_Gxw_new


!*********************************************************************
! 
! this routine determines the Bruckner RPA potential
! defined as
!
!  delta gamma_ai = int_0^infty  d tau
!   <a| sigma( tau) e( -tau |e_i|)  + e(-tau |e_a|) sigma( -tau) |i>
!  delta gamma_ia = int_0^infty d tau
!   <i| e( -tau |e_i|) sigma( tau)  + sigma(- tau) e(-tau |e_a|)  |a>
!
! on the cosine grid (since this is an even function)
!
! on input: 
!  SIGMAO_MAT and SIGMAU_MAT are the self-energies for a given tau point
!      and orbital representation distributed by scaLAPACK
!      SIGMAU_MAT is 
!
! on output:
!  CORR_MAT contains the contribution from the current tau point
!      to the change of the density matrix in orbital representation
!
!*********************************************************************

   SUBROUTINE ADD_G_SIGMA( W, SIGMAO_MAT, SIGMAU_MAT, CORR_MAT, & 
        NTAU, IMAG_GRIDS, NK, ISP, DESCA )
      USE scala
      USE constant
      IMPLICIT NONE
      TYPE (wavespin)        :: W
      GDEF                   :: SIGMAO_MAT(:) ! matrix <i| Sigma_o(tau) | a> stored distributed
      GDEF                   :: SIGMAU_MAT(:) ! matrix <i| Sigma_u(tau) | a> stored distributed     
      GDEF                   :: CORR_MAT(:)   ! contribution to density matrix from correlation effects
      INTEGER                :: NTAU          ! current tau point
      TYPE(imag_grid_handle) :: IMAG_GRIDS    ! time and frequency grid handle
      INTEGER                :: NK            ! current q-point
      INTEGER                :: ISP           ! current spin
      INTEGER                :: DESCA(:)          
    !local variables 
      REAL(q)                :: EFERMI        ! fermi energy 
      INTEGER                :: N, NB
      INTEGER                :: NTAU_GLOBAL
      REAL(q)                :: R( SIZE(W%CELTOT,1)), THETA(  SIZE(W%CELTOT,1))
      REAL(q)                :: TAU
      GDEF,ALLOCATABLE       :: SIGMA(:)

      PROFILING_START( 'ADD_G_SIGMA' )
      EFERMI = W%EFERMI(ISP)

      ALLOCATE(SIGMA(SIZE(SIGMAO_MAT,1)))

      !this is the global tau point of the node
      NTAU_GLOBAL = DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, NTAU,&
      IMAG_GRIDS )

      N= W%WDES%NB_TOTK(NK, ISP)
      R= REAL(W%CELTOT(:, NK, ISP),q)
      ! <a| sigma( tau) e( -tau |e_i|) |i>
      ! occupied Green function G_0(-tau) = exp( tau (e-EFERMI))
      ! calculate -SIGMA_unocc exp( tau (e-EFERMI))
      ! note that in VASP SIGMAU "lacks a minus sign" since the unoccupied G lacks the minus sign in the code
      THETA=0
      ! careful in case of finite temperature one requieres different potential
      IF ( LFINITE_TEMPERATURE ) THEN
         ! occupied greens function requiered
         DO NB=1,N
            THETA(NB)=GF_WEIGHT( .TRUE., IMAG_GRIDS%TAU(NTAU_GLOBAL), IMAG_GRIDS%T%BETA, EFERMI, &
               R(NB), REAL(W%FERTOT(NB,NK,ISP),q) )
         ENDDO
      ELSE
         DO NB=1,N
            IF (R(NB)-EFERMI<0.0_q) THEN
               THETA(NB)=EXP(IMAG_GRIDS%TAU(NTAU_GLOBAL)*(R(NB)-EFERMI))*REAL(W%FERTOT(NB,NK,ISP),q)
            ENDIF
         ENDDO
      ENDIF
       
      SIGMA=-SIGMAU_MAT
      CALL RIGHT_MULTIPLY_BY_R(N, SIGMA, THETA, DESCA )
      ! mask to unoccupied-occupied block
      ! CALL LEFT_MULTIPLY_BY_R(N, SIGMA, 1-W%FERTOT(:, NK, ISP), DESCA )
      CORR_MAT(:)=CORR_MAT(:)+IMAG_GRIDS%TAU_WEIGHT(NTAU_GLOBAL)*SIGMA(:)

      !    <i| e( -tau |e_i|) sigma( tau) |a>
      SIGMA=-SIGMAU_MAT
      CALL LEFT_MULTIPLY_BY_R(N, SIGMA, THETA, DESCA )
      ! mask to occupied-unoccupied block
      ! CALL RIGHT_MULTIPLY_BY_R(N, SIGMA, 1-W%FERTOT(:, NK, ISP), DESCA )
      CORR_MAT(:)=CORR_MAT(:)+IMAG_GRIDS%TAU_WEIGHT(NTAU_GLOBAL)*SIGMA(:)

      !  <i| sigma(- tau) e(-tau |e_a|)  |a>
      ! unoccupied Green function -exp( tau (EFERMI-e))
      ! calculate - (SIGMAO   -exp( tau (EFERMI-e)))
      THETA=0
      IF ( LFINITE_TEMPERATURE ) THEN
         ! unoccupied greens function requiered
         DO NB=1,N
            THETA(NB)=GF_WEIGHT( .FALSE., IMAG_GRIDS%TAU(NTAU_GLOBAL), IMAG_GRIDS%T%BETA, EFERMI, &
               R(NB), REAL(W%FERTOT(NB,NK,ISP),q) )
         ENDDO
      ELSE
         DO NB=1,N
            IF (R(NB)-EFERMI>=0.0_q) THEN
               THETA(NB)=EXP(IMAG_GRIDS%TAU(NTAU_GLOBAL)*(EFERMI-R(NB)))*(1-REAL(W%FERTOT(NB,NK,ISP),q))
            ENDIF
         ENDDO
      ENDIF
      SIGMA=SIGMAO_MAT
      CALL RIGHT_MULTIPLY_BY_R(N, SIGMA, THETA, DESCA )
      ! mask to occupied-unoccupied block
      ! CALL LEFT_MULTIPLY_BY_R(N, SIGMA, W%FERTOT(:, NK, ISP), DESCA )
      CORR_MAT(:)=CORR_MAT(:)+IMAG_GRIDS%TAU_WEIGHT(NTAU_GLOBAL)*SIGMA(:)

      !  <a| e(-tau |e_a|) sigma( -tau) |i>
      SIGMA=SIGMAO_MAT
      CALL LEFT_MULTIPLY_BY_R(N, SIGMA, THETA, DESCA )
      ! mask to unoccupied-occupied block
      ! CALL RIGHT_MULTIPLY_BY_R(N, SIGMA, W%FERTOT(:, NK, ISP), DESCA )
      CORR_MAT(:)=CORR_MAT(:)+IMAG_GRIDS%TAU_WEIGHT(NTAU_GLOBAL)*SIGMA(:)

      ! set occupied-occupied and unoccupied-unoccupied block from sigma(w->0)
      ! use this only if the masks are used in the previous lines (4 mask operations)
      !
      ! add sigma for occupied-occupied block sigma(w->0)
!      SIGMA=SIGMAO_MAT-SIGMAU_MAT
!      CALL RIGHT_MULTIPLY_BY_R(N, SIGMA, W%FERTOT(:, NK, ISP), DESCA )
!      CALL LEFT_MULTIPLY_BY_R(N, SIGMA,  W%FERTOT(:, NK, ISP), DESCA )
!      CORR_MAT(:)=CORR_MAT(:)+IMAG_GRIDS%TAU_WEIGHT(NTAU_GLOBAL)*SIGMA(:)

!      ! add sigma for unoccupied-unoccupied block sigma(w->0)
!      SIGMA=SIGMAO_MAT-SIGMAU_MAT
!      CALL RIGHT_MULTIPLY_BY_R(N, SIGMA, 1-W%FERTOT(:, NK, ISP), DESCA )
!      CALL LEFT_MULTIPLY_BY_R(N, SIGMA,  1-W%FERTOT(:, NK, ISP), DESCA )
!      CORR_MAT(:)=CORR_MAT(:)+IMAG_GRIDS%TAU_WEIGHT(NTAU_GLOBAL)*SIGMA(:)


!      the RPA total energy is invariant for unitary transformations within
!      the occupied (or unoccupied) block
!      here is a small test to show this
!      <k| sigma( tau) e( -tau |e_i|) - e( -tau |e_k|) sigma( tau) |i>
!      is a anti-Hermitian matrix (so unitary transformation leave the energy invariant)
!      this is analogous to density functional theory (energy invariant under unitary 
!      orbital transformations)
!      ! <k| sigma( tau) e( -tau |e_i|) |i>
!      THETA=0
!      DO NB=1,N
!         IF (R(NB)-EFERMI<0.0_q) THEN
!            THETA(NB)=EXP(IMAG_GRIDS%TAU(NTAU_GLOBAL)*(R(NB)-EFERMI))
!         ENDIF
!      ENDDO
!      SIGMA=-SIGMAU_MAT
!      CALL RIGHT_MULTIPLY_BY_R(N, SIGMA, THETA, DESCA )
!      ! mask to occupied-occupied block
!      CALL LEFT_MULTIPLY_BY_R(N, SIGMA, W%FERTOT(:, NK, ISP), DESCA )
!      CORR_MAT(:)=CORR_MAT(:)+IMAG_GRIDS%TAU_WEIGHT(NTAU_GLOBAL)*SIGMA(:)
!
!      ! - <k| e( -tau |e_k|) sigma( tau)  |i>
!      SIGMA=SIGMAU_MAT
!      CALL LEFT_MULTIPLY_BY_R(N, SIGMA, THETA, DESCA )
!      ! mask to occupied-occupied block
!      CALL RIGHT_MULTIPLY_BY_R(N, SIGMA, W%FERTOT(:, NK, ISP), DESCA )
!      CORR_MAT(:)=CORR_MAT(:)+IMAG_GRIDS%TAU_WEIGHT(NTAU_GLOBAL)*SIGMA(:)


! various alternatives
!      CALL FREQUENCYINT  ! test for alternative frequency integration (identical to above)

! int w G(w) sigma(w) G(w) w dw 
!      IF (IMAG_GRIDS%FER_IM_WEIGHT(1)==0) THEN
!         CALL ADD_G_SIGMA_Gxw_occocc
!      ELSE
!         CALL ADD_G_SIGMA_Gxw_sin_occocc
!      ENDIF

      DEALLOCATE(SIGMA)

      PROFILING_STOP( 'ADD_G_SIGMA' )

      CONTAINS

! this is a small test routine to transform Sigma(tau) to frequency and
! perform the integration in frequency regime
! obviously the relation \int d w Sigma(w) G(w)  =  \int d t Sigma(t) G(-t)
! implies that it is easier to perform the calculations in time
! this yields identical results as above if the masking
! and the occ-occ and unocc-unocc contributions are removed

    SUBROUTINE FREQUENCYINT
      INTEGER :: I
      COMPLEX(q),ALLOCATABLE :: SIGMA_COS(:), SIGMA_COS2(:)
      COMPLEX(q)         :: CELTOTINV( SIZE(W%CELTOT,1))
      
      ALLOCATE(SIGMA_COS(SIZE(SIGMAO_MAT,1)), SIGMA_COS2(SIZE(SIGMAO_MAT,1)))

      !
      ! only the term Re( G * Sigma * G ) is relevant, so 
      !
      SIGMA_COS=(0._q,0._q)
      SIGMA_COS2=(0._q,0._q)
      DO I=1,IMAG_GRIDS%NOMEGA
         ! take all four terms into account at once
         ! multiply \Sigma_cos from left and right with G and integrate only real part
         ! cosine transformation to i w
         SIGMA_COS(:)=(SIGMAO_MAT(:)-SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE(I,NTAU_GLOBAL)/2 &
                    - (SIGMAO_MAT(:)+SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE_CONJG(I,NTAU_GLOBAL)*(0._q,0.5_q)

         ! for the Kohn-Sham Green's function the expression is almost trivial
         CELTOTINV=1._q/(CMPLX(0._q,IMAG_GRIDS%FER_RE(I),q)-(W%CELTOT(:,NK,ISP)-EFERMI))
         !now compute Re(G.\Sigma.G)
         SIGMA_COS2=SIGMA_COS
         CALL RIGHT_MULTIPLY_BY_C(N, SIGMA_COS, CELTOTINV, DESCA )
         ! and add to total first order change of density matrix multiplied by weight
         CORR_MAT(:)=CORR_MAT(:)-SIGMA_COS(:)*(IMAG_GRIDS%FER_RE_WEIGHT(I)/2._q)/2

         SIGMA_COS=SIGMA_COS2
         CALL LEFT_MULTIPLY_BY_C(N, SIGMA_COS, CELTOTINV, DESCA )
         ! and add to total first order change of density matrix multiplied by weight
         CORR_MAT(:)=CORR_MAT(:)-SIGMA_COS(:)*(IMAG_GRIDS%FER_RE_WEIGHT(I)/2._q)/2
         ! take all four terms into account at once
         ! multiply \Sigma_cos from left and right with G and integrate only real part
         ! cosine transformation to i w
         SIGMA_COS(:)=(SIGMAO_MAT(:)-SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE(I,NTAU_GLOBAL)/2 &
                    + (SIGMAO_MAT(:)+SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE_CONJG(I,NTAU_GLOBAL)*(0._q,0.5_q)

         ! for the Kohn-Sham Green's function the expression is almost trivial
         CELTOTINV=1._q/(CMPLX(0._q,-IMAG_GRIDS%FER_IM(I),q)-(W%CELTOT(:,NK,ISP)-EFERMI))
         !now compute Re(G.\Sigma.G)
         SIGMA_COS2=SIGMA_COS
         CALL RIGHT_MULTIPLY_BY_C(N, SIGMA_COS, CELTOTINV, DESCA )
         ! and add to total first order change of density matrix multiplied by weight
         CORR_MAT(:)=CORR_MAT(:)-SIGMA_COS(:)*(IMAG_GRIDS%FER_IM_WEIGHT(I)/2._q)/2
         SIGMA_COS=SIGMA_COS2
         CALL LEFT_MULTIPLY_BY_C(N, SIGMA_COS, CELTOTINV, DESCA )
         ! and add to total first order change of density matrix multiplied by weight
         CORR_MAT(:)=CORR_MAT(:)-SIGMA_COS(:)*(IMAG_GRIDS%FER_IM_WEIGHT(I)/2._q)/2

      ENDDO

      DEALLOCATE(SIGMA_COS,SIGMA_COS2)

   END SUBROUTINE FREQUENCYINT

   SUBROUTINE ADD_G_SIGMA_Gxw_occocc
    !local variables 
      INTEGER            :: I
      COMPLEX(q),ALLOCATABLE :: SIGMA_w(:)
      COMPLEX(q)         :: CELTOTINV( SIZE(W%CELTOT,1))

      ALLOCATE(SIGMA_w(SIZE(SIGMAO_MAT,1)))

      !this is the global tau point of the node
      NTAU_GLOBAL = DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, NTAU,&
      IMAG_GRIDS )
      !
      ! only the term Re( G * Sigma * G ) is relevant, so 
      !
      SIGMA_w=(0._q,0._q)
      DO I=1,IMAG_GRIDS%NOMEGA
         ! take all four terms into account at once
         ! multiply \Sigma_cos from left and right with G and integrate only real part
         ! cosine transformation to i w
         SIGMA_w(:)=(SIGMAO_MAT(:)-SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE(I,NTAU_GLOBAL)/2 &
                  - (SIGMAO_MAT(:)+SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE_CONJG(I,NTAU_GLOBAL)*(0._q,0.5_q)

         ! for the Kohn-Sham Green's function the expression is almost trivial
         CELTOTINV=1._q/(CMPLX(0._q,IMAG_GRIDS%FER_RE(I),q)-(W%CELTOT(:,NK,ISP)-EFERMI))
         !now compute Re(G.\Sigma.G)
         CALL MULTIPLY_BY_IVEIGENVALUE(N, SIGMA_w, CELTOTINV, DESCA)
         ! and add to total first order change of density matrix multiplied by weight
         CALL CLEAR_OCCUNOCC_UNOCCOCC_FER(N, SIGMA_w, W%FERTOT(:,NK,ISP), DESCA)
         CORR_MAT(:)=CORR_MAT(:)+SIGMA_w(:)*(IMAG_GRIDS%FER_RE_WEIGHT(I)/2._q*CMPLX(EFERMI,IMAG_GRIDS%FER_RE(I),q))
         ! make up for special scaling in OMEGAWEIGHT, and multiply with (i\omega+e_fermi)

         ! take all four terms into account at once
         ! multiply \Sigma_cos from left and right with G and integrate only real part
         ! cosine transformation to i w
         SIGMA_w(:)=(SIGMAO_MAT(:)-SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE(I,NTAU_GLOBAL)/2 &
                  + (SIGMAO_MAT(:)+SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE_CONJG(I,NTAU_GLOBAL)*(0._q,0.5_q)

         ! for the Kohn-Sham Green's function the expression is almost trivial
         CELTOTINV=1._q/(CMPLX(0._q,-IMAG_GRIDS%FER_RE(I),q)-(W%CELTOT(:,NK,ISP)-EFERMI))
         !now compute Re(G.\Sigma.G)
         CALL MULTIPLY_BY_IVEIGENVALUE(N, SIGMA_w, CELTOTINV, DESCA)
         ! and add to total first order change of density matrix multiplied by weight
         CALL CLEAR_OCCUNOCC_UNOCCOCC_FER(N, SIGMA_w, W%FERTOT(:,NK,ISP), DESCA)
         CORR_MAT(:)=CORR_MAT(:)+SIGMA_w(:)*(IMAG_GRIDS%FER_RE_WEIGHT(I)/2._q*CMPLX(EFERMI,-IMAG_GRIDS%FER_RE(I),q))
         ! make up for special scaling in OMEGAWEIGHT, and multiply with \omega
      ENDDO

      DEALLOCATE(SIGMA_w)

   END SUBROUTINE ADD_G_SIGMA_Gxw_occocc


   SUBROUTINE ADD_G_SIGMA_Gxw_sin_occocc
    !local variables 
      INTEGER            :: I
      COMPLEX(q),ALLOCATABLE :: SIGMA_w(:)
      COMPLEX(q)         :: CELTOTINV( SIZE(W%CELTOT,1))
      
      ALLOCATE(SIGMA_w(SIZE(SIGMAO_MAT,1)))

      !this is the global tau point of the node
      NTAU_GLOBAL = DETERMINE_NTAU_GLOBAL( W%WDES%COMM_KIN%NODE_ME, NTAU,&
      IMAG_GRIDS )
      
      !First part only calculates the contribution of SIGMA_even on COS-grid
      SIGMA_w=(0._q,0._q)
      DO I=1,IMAG_GRIDS%NOMEGA
         ! multiply \Sigma_w (even) from left and right with G and perform integration for \omega > 0:
         !------------------------------------------------------------------------
         ! cosine transformation of SIGMA_even to i w
         SIGMA_w(:)=(SIGMAO_MAT(:)-SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE(I,NTAU_GLOBAL)/2
         !first the even part for \omega > 0 is calculated
         ! for the Kohn-Sham Green's function the expression G(iw) is almost trivial
         CELTOTINV=1._q/(CMPLX(0._q,IMAG_GRIDS%FER_RE(I),q)-(W%CELTOT(:,NK,ISP)-EFERMI))
         !now compute (G.\Sigma_even.G)(iw)
         CALL MULTIPLY_BY_IVEIGENVALUE(N, SIGMA_w, CELTOTINV, DESCA)
         ! and calculate 1st part of the first order change of the density matrix 
         ! due to the change of the overlap operator S, by integrating 
         ! (G.\Sigma_even.G)(iw) (i\omega+e_fermi) on the COS-grid from 0 to \infty
         CALL CLEAR_OCCUNOCC_UNOCCOCC_FER(N, SIGMA_w, W%FERTOT(:,NK,ISP), DESCA)

         CORR_MAT(:)=CORR_MAT(:)+SIGMA_w(:)*(IMAG_GRIDS%FER_RE_WEIGHT(I)/2._q*CMPLX(EFERMI,IMAG_GRIDS%FER_RE(I),q))
         !-------------------------------------------------------------------------

         ! repeat for \omega < 0
         SIGMA_w(:)=(SIGMAO_MAT(:)-SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_RE(I,NTAU_GLOBAL)/2
         CELTOTINV=1._q/(CMPLX(0._q,-IMAG_GRIDS%FER_RE(I),q)-(W%CELTOT(:,NK,ISP)-EFERMI))
         !now compute (G.\Sigma_even.G)(-iw)
         CALL MULTIPLY_BY_IVEIGENVALUE(N, SIGMA_w, CELTOTINV, DESCA)
         CALL CLEAR_OCCUNOCC_UNOCCOCC_FER(N, SIGMA_w, W%FERTOT(:,NK,ISP), DESCA)

         CORR_MAT(:)=CORR_MAT(:)+SIGMA_w(:)*(IMAG_GRIDS%FER_RE_WEIGHT(I)/2._q*CMPLX(EFERMI,-IMAG_GRIDS%FER_RE(I),q))

      ENDDO

      !Second part only calculates the contribution of SIGMA_odd on SIN-grid
      SIGMA_w=(0._q,0._q)
      DO I=1,IMAG_GRIDS%NOMEGA
         ! multiply \Sigma_w (odd) from left and right with G and perform integration for \omega > 0:
         !------------------------------------------------------------------------
         ! sine transformation of SIGMA_odd to i w
         SIGMA_w(:)=-(SIGMAO_MAT(:)+SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_IM(I,NTAU_GLOBAL)*(0._q,0.5_q)
         ! for the Kohn-Sham Green's function the expression G(iw) is almost trivial
         CELTOTINV=1._q/(CMPLX(0._q,IMAG_GRIDS%FER_IM(I),q)-(W%CELTOT(:,NK,ISP)-EFERMI))
         !now compute (G.\Sigma_odd.G)(iw)
         CALL MULTIPLY_BY_IVEIGENVALUE(N, SIGMA_w, CELTOTINV, DESCA)
         ! and calculate 1st part of the first order change of the density matrix 
         ! due to the change of the overlap operator S, by integrating 
         ! (G.\Sigma_odd.G)(iw) (i\omega+e_fermi) on the COS-grid from 0 to \infty
         CALL CLEAR_OCCUNOCC_UNOCCOCC_FER(N, SIGMA_w, W%FERTOT(:,NK,ISP), DESCA)
         CORR_MAT(:)=CORR_MAT(:)+SIGMA_w(:)*(IMAG_GRIDS%FER_IM_WEIGHT(I)/2._q*(CMPLX(EFERMI,IMAG_GRIDS%FER_IM(I),q)))
         !-------------------------------------------------------------------------

         ! repeat for \omega < 0
         SIGMA_w(:)= (SIGMAO_MAT(:)+SIGMAU_MAT(:))*IMAG_GRIDS%TO_FER_IM(I,NTAU_GLOBAL)*(0._q,0.5_q)
         CELTOTINV=1._q/(CMPLX(0._q,-IMAG_GRIDS%FER_IM(I),q)-(W%CELTOT(:,NK,ISP)-EFERMI))
         !now compute (G.\Sigma_even.G)(-iw)
         CALL MULTIPLY_BY_IVEIGENVALUE(N, SIGMA_w, CELTOTINV, DESCA)
         CALL CLEAR_OCCUNOCC_UNOCCOCC_FER(N, SIGMA_w, W%FERTOT(:,NK,ISP), DESCA)
         CORR_MAT(:)=CORR_MAT(:)+SIGMA_w(:)*(IMAG_GRIDS%FER_IM_WEIGHT(I)/2._q*(CMPLX(EFERMI,-IMAG_GRIDS%FER_IM(I),q)))
         
      ENDDO

      DEALLOCATE(SIGMA_w)

    END SUBROUTINE ADD_G_SIGMA_Gxw_sin_occocc


    END SUBROUTINE ADD_G_SIGMA

!*********************************************************************
! 
! this routine sets the diagonal of the approximate Hermitian potential
! to
!   int_0^infty G(iw) iw d w
! it yields pretty identical results as the perturbative approach
! implemented above (at least at convergence)
!
!*********************************************************************

   SUBROUTINE SET_G_SIGMA_DIAG( W, EFERMI, G_COS, G_SIN, & 
        CORR_MAT, IMAG_GRIDS, DESCA )
      USE scala
      USE constant
      IMPLICIT NONE
      TYPE (wavespin)        :: W
      REAL(q)                :: EFERMI(W%WDES%ISPIN)
      COMPLEX(q)             :: G_COS(:,:,:,:)      !diagonal of cos transformed GF in bloch domain 
      COMPLEX(q)             :: G_SIN(:,:,:,:)      !diagonal of sin transformed GF in bloch domain 
      GDEF                   :: CORR_MAT(:,:,:,:)   ! contribution to density matrix from correlation effects
      TYPE(imag_grid_handle) :: IMAG_GRIDS          ! time and frequency grid handle
      INTEGER                :: DESCA(:)
   ! local
      INTEGER                :: ISP, NK, N, NB, I
      REAL(q)                :: E
      COMPLEX(q)             :: GW, GW_INT
      REAL(q)                :: GW_DIAG(1:W%WDES%NB_TOT)


      DO ISP=1,W%WDES%ISPIN
      DO NK=1,W%WDES%NKPTS
         N= W%WDES%NB_TOTK(NK,ISP)
         ! loop over all bands
         GW_DIAG=0
         DO NB=1, N
            E=W%CELTOT(NB,NK,ISP)-EFERMI(ISP)
            GW_INT=0
            DO I=1,SIZE(G_COS,1)  ! loop over all frequency points
               ! diagonal of Green function  G(iw) = 1/(i w + EFERMI - Sigma_c - H_HF)
               GW=   1.0_q/(CMPLX(0._q,IMAG_GRIDS%FER_RE(I),q)-(G_COS(I,NB,NK,ISP)- G_SIN(I,NB,NK,ISP)*(0.0_q,1.0_q)+E))
               ! now subtract mean field Green's function, but shift poles if it is too close for 
               ! numerical integration
               GW=GW-1.0_q/(CMPLX(0._q,IMAG_GRIDS%FER_RE(I),q)-E)
               GW_INT=GW_INT+GW*(IMAG_GRIDS%FER_RE_WEIGHT(I)/2._q*CMPLX(EFERMI(ISP),IMAG_GRIDS%FER_RE(I),q))

               ! negative w
               ! diagonal of Green function  G(iw) = 1/(i w + EFERMI - Sigma_c - H_HF)
               GW=   1.0_q/(CMPLX(0._q,-IMAG_GRIDS%FER_RE(I),q)-(G_COS(I,NB,NK,ISP)+ G_SIN(I,NB,NK,ISP)*(0.0_q,1.0_q)+E))
               ! now subtract mean field Green's function, but shift poles if it is too close for 
               ! numerical integration
               GW=GW-1.0_q/(CMPLX(0._q,-IMAG_GRIDS%FER_RE(I),q)-E)
               GW_INT=GW_INT+GW*(IMAG_GRIDS%FER_RE_WEIGHT(I)/2._q*CMPLX(EFERMI(ISP),-IMAG_GRIDS%FER_RE(I),q))

            ENDDO
            ! finally add eigenvalue E to obtain centre of QP
            GW_DIAG(NB)=GW_INT ! +W%CELTOT(NB,NK,ISP)
         ENDDO
         IF (W%WDES%COMM%NODE_ME==W%WDES%COMM%IONODE) THEN
            WRITE(*,'(8F14.7)') GW_DIAG(1:24)
         ENDIF
!         CALL SET_DIAGONALE_REAL(N, CORR_MAT(:, NK, ISP, 1), GW_DIAG, DESCA )
      
      ENDDO  !k-point
      ENDDO  !spin
      
    END SUBROUTINE SET_G_SIGMA_DIAG

!*********************************************************************
!
! on entry, SIGMA_COS, SIGMA_SIN, SIGMA_SIN_SIN
! on exit, 2/3 of SIGMA_SIN is replaced by the spline fit of SIGMA_SIN_SIN
!
!*********************************************************************
   SUBROUTINE SPLINE_FIT_SIN(W, SIGMA_COS, SIGMA_SIN, SIGMA_SIN_SIN, NKPTS_IRZ, IMAG_GRIDS)
      USE wave
      USE minimax_struct
      IMPLICIT NONE
      TYPE (wavespin) :: W                       !wave function 
      COMPLEX(q)      :: SIGMA_COS(:,:,:,:)      !even part of the self-energy (cos grids) 
      COMPLEX(q)      :: SIGMA_SIN(:,:,:,:)      !odd part of the self-energy (cos grids) 
      COMPLEX(q)      :: SIGMA_SIN_SIN(:,:,:,:)  !odd part of the self-energy at sine grids
      INTEGER         :: NKPTS_IRZ
      TYPE(imag_grid_handle)  :: IMAG_GRIDS      !time and frequency grid handle
      !local variables 
      INTEGER                :: I,NB,NK,ISP
      REAL(q), ALLOCATABLE   :: P(:,:),F(:)
      REAL(q)                :: FDER

      ALLOCATE(P(2*IMAG_GRIDS%NOMEGA,5)) ; P=0.0_q
      ALLOCATE(F(IMAG_GRIDS%NOMEGA)) ; F=0.0_q

      !loop over all Bloch states 
      DO ISP=1, W%WDES%ISPIN
        DO NK=1,NKPTS_IRZ
           DO NB=1,W%WDES%NB_TOT

                !spline interpolate SIGMA_SIN_SIN with open boundary conditions
                P=0.0_q 
                DO I=1,IMAG_GRIDS%NOMEGA
                   P(I,1)=-IMAG_GRIDS%FER_IM(IMAG_GRIDS%NOMEGA-I+1)
                   P(I,2)=-REAL(SIGMA_SIN_SIN(IMAG_GRIDS%NOMEGA-I+1,NB,NK,ISP),q)
                   P(I+IMAG_GRIDS%NOMEGA,1)=IMAG_GRIDS%FER_IM(I)
                   P(I+IMAG_GRIDS%NOMEGA,2)=REAL(SIGMA_SIN_SIN(I,NB,NK,ISP),q)
                ENDDO
                !first, get the coefficients of spline fit
                CALL SPLCOF(P,2*IMAG_GRIDS%NOMEGA,2*IMAG_GRIDS%NOMEGA,10E30_q)
                !second, spline-interpolate SIGMA_SIN_SIN to cos grids
                F=0.0_q
                DO I=1,IMAG_GRIDS%NOMEGA
                   CALL SPLVAL(IMAG_GRIDS%FER_RE(I),F(I),FDER,P,2*IMAG_GRIDS%NOMEGA,2*IMAG_GRIDS%NOMEGA)
                ENDDO

                !now we overwrite the first 2/3 of SIGMA_SIN data points with interpolated data
                DO I=1,MAX(1,2*IMAG_GRIDS%NOMEGA/3)
                   SIGMA_SIN(I,NB,NK,ISP)=F(I)
                ENDDO

           ENDDO   !NB
        ENDDO   !NK
      ENDDO   !ISP

      DEALLOCATE(P,F)
    END SUBROUTINE SPLINE_FIT_SIN

!*********************************************************************
!
! write self-energy to vasprun.xml
! (using the diagonal part of the self-energy matrix)
!
!*********************************************************************

   SUBROUTINE WRITE_SELF_ENERGY(W, G_COS, G_SIN, NKPTS_IRZ, IMAG_GRIDS, NBANDSGW, IU6)
      USE constant
      USE wave
      USE vaspxml
      IMPLICIT NONE
      TYPE (wavespin) :: W               !wave function 
      COMPLEX(q)      :: G_COS(:,:,:,:)  !diagonal of cos transformed GF in bloch domain 
      COMPLEX(q)      :: G_SIN(:,:,:,:)  !diagonal of sin transformed GF in bloch domain 
      INTEGER         :: NKPTS_IRZ
      TYPE(imag_grid_handle)  :: IMAG_GRIDS        !time and frequency grid handle
      INTEGER         :: NBANDSGW        !number of bands to be written to file
      INTEGER         :: IU6             !for output 
      !local variables 
      INTEGER                :: NOMEGA_TOT,I
      INTEGER                :: NB,NK,ISP
      REAL(q), ALLOCATABLE :: A(:,:)

      !if the frequencies for the sine transformation are distinct
      !the code makes a quick return for the time being
      IF ( .NOT. LSELFENERGY ) RETURN

      PROFILING_START( 'WRITE_SELF_ENERGY' )

      ALLOCATE(A(3, IMAG_GRIDS%NOMEGA))

      IF (IU6>0) THEN
      !loop over all Bloch states 
      DO ISP=1, W%WDES%ISPIN
      DO NK=1,NKPTS_IRZ
         DO NB=1,NBANDSGW
            !store Green's function to auxillary arrays            
            DO I=1,IMAG_GRIDS%NOMEGA
               A(1,I)=IMAG_GRIDS%FER_RE(I)
               A(2,I)= G_COS(I,NB,NK,ISP)+W%CELTOT(NB,NK,ISP)
               A(3,I)=-G_SIN(I,NB,NK,ISP)
            ENDDO
            CALL XML_VECARRAY("selfenergy along imaginary axis")
            CALL XML_ARRAY_REAL(A,'(F24.16)')
            CALL XML_CLOSE_TAG
         ENDDO
      ENDDO                     !k-point
      ENDDO                     !spin
      ENDIF

      DEALLOCATE(A)
      PROFILING_STOP( 'WRITE_SELF_ENERGY' )
    END SUBROUTINE WRITE_SELF_ENERGY

!***********************************************************************
!
! update potential and charge
!
!***********************************************************************

  SUBROUTINE POTENTIAL_AND_CHARGE(NELM_HF,&
     W, WDES, E, INFO, &
     GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
     LATT_CUR, T_INFO, P, SYMM, MIX, &
     CHTOT, CHTOTL, DENCOR, CVTOT, CSTRF, IRDMAX, &
     LMDIM, CQIJ, CDIJ, CRHODE, N_MIX_PAW, RHOLM, RHOLM_LAST, CHDEN, SV, & 
     IO, OEP)
    USE us, ONLY: SET_CHARGE
    USE pot, ONLY: POTLOK
    USE pawm, ONLY: SET_DD_PAW
    USE setexm, ONLY: SETUP_LDA_XC
    USE pseudo_struct_def 
    INTEGER, INTENT(IN) :: NELM_HF
! structures
    TYPE (wavespin)     :: W
    TYPE (wavedes)      :: WDES
    TYPE (energy)       :: E
    TYPE (info_struct)  :: INFO
    TYPE (type_info)    :: T_INFO
    TYPE (potcar)       :: P(T_INFO%NTYP)
    TYPE (latt)         :: LATT_CUR
    TYPE (grid_3d)      :: GRID       ! grid for wavefunctions
    TYPE (grid_3d)      :: GRID_SOFT  ! grid for soft chargedensity
    TYPE (grid_3d)      :: GRIDC      ! grid for potentials/charge
    TYPE (grid_3d)      :: GRIDUS     ! temporary grid in us.F
    TYPE (transit)      :: C_TO_US    ! index table between GRIDC and GRIDUS
    TYPE (transit)      :: SOFT_TO_C  ! index table between GRID_SOFT and GRIDC
    TYPE (symmetry)     :: SYMM
    TYPE (mixing)       :: MIX
    TYPE (in_struct)    :: IO
! charge 
    COMPLEX(q) CHTOTL(GRIDC%MPLWV,WDES%NCDIJ)! old charge-density
    RGRID      SV(DIMREAL(GRID%MPLWV), WDES%NCDIJ)
    COMPLEX(q) CSTRF(GRIDC%MPLWV, T_INFO%NTYP), &
               CHTOT(GRIDC%MPLWV, WDES%NCDIJ), &
               CVTOT(GRIDC%MPLWV, WDES%NCDIJ)
    RGRID      DENCOR(GRIDC%RL%NP)
    OVERLAP  CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    OVERLAP  CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    INTEGER     :: IRDMAX
    INTEGER     :: LMDIM
!  paw sphere charge density
    INTEGER N_MIX_PAW
    REAL(q)  RHOLM(N_MIX_PAW,WDES%NCDIJ),RHOLM_LAST(N_MIX_PAW,WDES%NCDIJ)
!  charge-density and potential on soft grid
    COMPLEX(q)  CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)

    TYPE (oep_handle), POINTER, OPTIONAL    :: OEP
    ! local 
    REAL(q) :: RMST
    INTEGER :: IRDMAA             ! temporary
    REAL(q) :: XCSIF(3,3)         ! temporary (stress tensor)
    INTEGER :: ISP, K1, N
    REAL(q) :: EMAX 

    PROFILING_START( 'potential_and_charge' )

    IF ( LOEP ) THEN
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
    ENDIF

    ! update charge
    CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
         GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
         LATT_CUR, P, SYMM, T_INFO, &
         CHDEN, LMDIM, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX)

    ! backup original potential for OEP
    IF ( LOEP ) THEN
       CALL ALLOCATE_OEP_HANDLE( OEP, GRID, GRIDC, LMDIM, WDES )
       OEP%SV = SV
       OEP%CDIJ = CDIJ
       OEP%CVTOT = CVTOT
       ! force full fock contribution
       HFSCREEN=0.0
       AEXX    =1.0
    ENDIF

    ! mix charge density with the one of the previous step
    ! copy from chi.F, but maybe not required
    IF (MIX%IMIX/=0) INFO%TIME=1.0
    IF (NELM_HF>=2 .AND. MIX%IMIX/=0 ) THEN
       RMST=0
       CALL MIX_SIMPLE(GRIDC,MIX,WDES%NCDIJ, CHTOT, CHTOTL, &
            N_MIX_PAW, RHOLM, RHOLM_LAST, LATT_CUR%B, LATT_CUR%OMEGA, RMST)
       IF (IO%IU6>=0) WRITE(IO%IU6,'(" charge density residual (rmsc) ",E14.5)') RMST
    ENDIF

    ! store old density to CHTOTL and ROHML_LAST (for mixer in next step)
    DO ISP=1,WDES%NCDIJ
       CALL RC_ADD(CHTOT(1,ISP),1.0_q,CHTOT(1,ISP),0.0_q,CHTOTL(1,ISP),GRIDC)
    ENDDO
    RHOLM_LAST=RHOLM

    ! calculate potential (Hartree + ionic only)
    CALL POTLOK(GRID,GRIDC,GRID_SOFT, WDES%COMM_INTER, WDES,  &
         INFO,P,T_INFO,E,LATT_CUR, &
         CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF)

    IF ( LOEP ) THEN
       CDIJ=0
       ! add the one center augmentation related terms  
       CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
            LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

       ! force compuation of fock part in PAW, USEFOCK_CONTRIBUTION will return T
       LHFCALC_FORCE=.TRUE.
       CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
            WDES%NCDIJ, LMDIM, CDIJ(1,1,1,1),  RHOLM, CRHODE(1,1,1,1), &
            E,  LMETA=.FALSE., LASPH=INFO%LASPH, LCOREL= .FALSE.  )
       LHFCALC_FORCE=.FALSE.
    ELSE
       ! add the one center augmentation related terms  
       CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
            LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

       ! finally add one center terms
       CALL SET_DD_PAW(WDES, P , T_INFO, INFO%LOVERL, &
            WDES%NCDIJ, LMDIM, CDIJ(1,1,1,1),  RHOLM, CRHODE(1,1,1,1), &
            E, .FALSE., INFO%LASPH, .FALSE.  )
    ENDIF

    PROFILING_STOP( 'potential_and_charge' )
  END SUBROUTINE POTENTIAL_AND_CHARGE


!************************************************************************
!  
! DESCRIPTION:                                      
!> allocates OEP handle that stores potential terms 
!************************************************************************
  SUBROUTINE ALLOCATE_OEP_HANDLE( OEP, GRID, GRIDC, LMDIM, WDES ) 
     TYPE( oep_handle ), POINTER :: OEP
     TYPE( grid_3d )             :: GRID
     TYPE( grid_3d )             :: GRIDC
     INTEGER                     :: LMDIM
     TYPE( wavedes )             :: WDES
     ! local
     INTEGER                     :: ISTAT
     IF( ASSOCIATED( OEP ) ) RETURN 

     ALLOCATE( OEP )  

     ALLOCATE( OEP%SV(DIMREAL(WDES%GRID%MPLWV), WDES%NCDIJ), STAT=ISTAT )
     IF ( ISTAT/=0 ) THEN
        CALL VTUTOR%ERROR( "ALLOCATE_OEP_HANDLE is not able to allocate "//&
          str( 8*DIMREAL(WDES%GRID%MPLWV)*WDES%NCDIJ) //&
          " kB of data for OEP%SV." )
     ENDIF
     OEP%SV = 0 

     ALLOCATE( OEP%CVTOT(GRIDC%MPLWV, WDES%NCDIJ), STAT = ISTAT )
     IF ( ISTAT/=0 ) THEN
        CALL VTUTOR%ERROR( "ALLOCATE_OEP_HANDLE is not able to allocate "//&
          str( 16*GRIDC%MPLWV*WDES%NCDIJ) //&
          " kB of data for OEP%CVTOT." )
     ENDIF
     OEP%CVTOT = 0 
     ALLOCATE( OEP%CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ) )
     OEP%CDIJ = 0 

     ALLOCATE( OEP%ESHIFT(WDES%NCDIJ ) ) 
     OEP%ESHIFT = 0 

  END SUBROUTINE ALLOCATE_OEP_HANDLE

!************************************************************************
!  
! DESCRIPTION:                                      
!> deallocates OEP handle that stores potential terms 
!************************************************************************
  SUBROUTINE DEALLOCATE_OEP_HANDLE( OEP )
     TYPE( oep_handle ), POINTER :: OEP
     TYPE( grid_3d )             :: GRID
     TYPE( grid_3d )             :: GRIDC
     INTEGER                     :: NCDIJ
     IF( .NOT. ASSOCIATED( OEP ) ) RETURN 

     IF ( ASSOCIATED( OEP%SV ) ) DEALLOCATE( OEP%SV )
     NULLIFY( OEP%SV )
     IF ( ASSOCIATED( OEP%CVTOT ) ) DEALLOCATE( OEP%CVTOT )
     NULLIFY( OEP%CVTOT )
     IF ( ASSOCIATED( OEP%CDIJ ) ) DEALLOCATE( OEP%CDIJ )
     NULLIFY( OEP%CDIJ )
     IF ( ASSOCIATED( OEP%ESHIFT ) ) DEALLOCATE( OEP%ESHIFT )
     NULLIFY( OEP%ESHIFT )

     IF ( ASSOCIATED( OEP ) ) DEALLOCATE( OEP )
     NULLIFY( OEP ) 
     
  END SUBROUTINE DEALLOCATE_OEP_HANDLE

!************************************************************************
!  
! DESCRIPTION: 
!> Implements OEP method, specifically the change in \f$ V_{xc} \f$
!> \f$
!>   \Delta V_{xc} = V_{xc} - V_{KS}
!> \f$
!>
!> derived from linearizing following expression 
!> 
!> \f$
!>  0 = G_{KS} - G = G_{KS}( \Sigma_{xc} -V_{KS} )G|_{ G = G_{KS} }
!> \f$
!>  
!> which after addition and removal of \f$ H_{KS} = T + V_h + V_{KS} \f$ 
!> yields 
!>  
!> \f$
!> \chi \Delta V_{xc} = [ ( F_{ij} - \epsilon^{KS}_j \delta_{ij} ) 
!> \frac{ f_j - f_i }{\epsilon_i - \epsilon_j} + 
!>          + \Gamma_{ij}^c ]_{ j=i } 
!> \f$
!> 
!> @param[in] CHI0 contains the polarizability at \f$ \tau=0-\f$
!> @param[in] W contains the eigenvalues of the r.h.s. matrix
! 
!************************************************************************

   SUBROUTINE OEP_POTENTIAL( CHI0, GDES, S2E, WGW, IMAG_GRIDS, &
      W, WDES, INFO, E,&
      GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
      LATT_CUR, SYMM, T_INFO, P, &
      LMDIM, CDIJ, CQIJ, SV, &
      CHDEN, CRHODE, CHTOT, RHOLM, N_MIX_PAW, IRDMAX, &
      DENCOR, CVTOT, CSTRF, &
      AMIX, OEP, IO )

      USE tutor, ONLY: vtutor
      USE chi_base, ONLY: SET_MAT_FROM_RESPONSE
      USE scala, ONLY: PSEUDO_INVERSE_DESC
      USE ini, ONLY: REGISTER_ALLOCATE, DEREGISTER_ALLOCATE
      USE us, ONLY: SET_CHARGE
      USE pot, ONLY:  SET_SV
      USE poscar_struct_def 
      USE pseudo_struct_def 
      ! chi related
      TYPE (responsefunction), POINTER :: CHI0(:)           ! response function at w=0
      TYPE (greensfdes)         :: GDES
      TYPE (screened_2e_handle) :: S2E
      TYPE (wavedes), POINTER   :: WGW   ! response function descriptor
      TYPE (imag_grid_handle)   :: IMAG_GRIDS
      ! charge related
      TYPE (wavespin)           :: W     ! r.h.s. matrix 
      TYPE (wavedes)            :: WDES  ! wave function descriptor
      TYPE (info_struct)        :: INFO
      TYPE (energy)             :: E 

      TYPE (grid_3d)            :: GRID       ! grid for wavefunctions
      TYPE (grid_3d)            :: GRIDC      ! grid for potentials/charge
      TYPE (grid_3d)            :: GRID_SOFT  ! grid for soft chargedensity
      TYPE (grid_3d)            :: GRIDUS     ! temporary grid in us.F
      TYPE (transit)            :: C_TO_US    ! index table between GRIDC and GRIDUS
      TYPE (transit)            :: SOFT_TO_C  ! index table between GRID_SOFT and GRIDC

      TYPE (latt)               :: LATT_CUR
      TYPE (symmetry)           :: SYMM
      TYPE (type_info)          :: T_INFO
      TYPE (potcar)             :: P(T_INFO%NTYP)

      INTEGER                   :: LMDIM
      OVERLAP                   :: CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ) ! local potential
      OVERLAP                   :: CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ) ! local potential
      RGRID                     :: SV(DIMREAL(WDES%GRID%MPLWV),WDES%NCDIJ) ! local potential

      COMPLEX(q)                :: CHDEN(GRID_SOFT%MPLWV,WDES%NCDIJ)
      OVERLAP                   :: CRHODE(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q)                :: CHTOT(GRIDC%MPLWV, WDES%NCDIJ)
      REAL(q)                   :: RHOLM(N_MIX_PAW,WDES%NCDIJ)
      INTEGER                   :: N_MIX_PAW
      INTEGER                   :: IRDMAX

      RGRID                     :: DENCOR(GRIDC%RL%NP)
      COMPLEX(q)                :: CVTOT(GRIDC%MPLWV,WDES%NCDIJ)
      COMPLEX(q)                :: CSTRF(GRIDC%MPLWV,T_INFO%NTYP)

      REAL(q)                   :: AMIX
      TYPE (oep_handle), POINTER ::  OEP         ! oep handle 
      TYPE (in_struct)          :: IO
      ! local 
      GDEF, POINTER, CONTIGUOUS :: CHI_INV(:)=>NULL() ! response matrix workspace
      INTEGER                   :: I, J
      INTEGER                   :: ISP
      INTEGER                   :: NQ, NQ_COUNTER
      TYPE (wavedes1)           :: WGWQ
      TYPE (wavedes1)           :: WDESQ
      INTEGER                   :: NP
      INTEGER                   :: NP1, NP2
      INTEGER                   :: NP2_GLB, NP2_END
      INTEGER                   :: IROW, JCOL
      TYPE(greens_mat_des), POINTER :: RDES1D=>NULL() ! distributed response
      TYPE(greens_mat_des), POINTER :: RDES2D=>NULL() ! distributed response
      INTEGER                   :: ISTAT
      ! OEP charge quantities
      COMPLEX(q), ALLOCATABLE   :: CHG(:,:)
      GDEF, ALLOCATABLE         :: CHG1(:)
      GDEF, ALLOCATABLE         :: CHG2(:)
      RGRID                     :: SV1(DIMREAL(WDES%GRID%MPLWV),WDES%NCDIJ) ! local potential
      COMPLEX(q)                :: CHTOT2(GRIDC%MPLWV, WDES%NCDIJ)
      OVERLAP                   :: CRHODE2(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      COMPLEX(q)                :: CHDEN2(GRID_SOFT%MPLWV,WDES%NCDIJ)
      REAL(q)                   :: RHOLM2(N_MIX_PAW,WDES%NCDIJ)
      COMPLEX(q), ALLOCATABLE   :: CV1(:)
      ! dummy arguments needed for determination of pot 
      INTEGER                   :: IRDMAA     ! dummy arguments to determine pot
      REAL(q)                   :: XCSIF(3,3)

      INTEGER                   :: NODE_ME = 1 
!#define unit_test
#ifdef unit_test
      GDEF, ALLOCATABLE         :: CHI_TMP(:)
      GDEF, ALLOCATABLE        :: CHI2(:)
#endif
      IF ( .NOT. LOEP ) RETURN

DOESI SV1=0

      PROFILING_START( 'OEP_POTENTIAL' )
#ifdef MPI
      NODE_ME = IMAG_GRIDS%B%COMM_IN_GROUP%NODE_ME
#endif

      ! determine potential -V_OEP by calculating the difference of
      ! current local potential, ionic and Hartree contribution 

      ! set q-point to Gamma
      qloop: DO NQ_COUNTER=1,S2E%NUMBER_OF_NQ 
         NQ=S2E%NQ(NQ_COUNTER)
         CALL SETWDES(WGW, WGWQ, NQ )
         IF( SQRT( DOT_PRODUCT(WDES%VKPT(:,NQ), WDES%VKPT(:,NQ)) )<1.E-8_q ) EXIT  qloop
      ENDDO qloop 

      IF ( NQ_COUNTER /= 1) THEN
         CALL vtutor%bug("OEP_POTENTIAL: Gamma-point not found", __FILE__, __LINE__)
      ENDIF

      ! WGWQ is now the response function descriptor at Gamma
      NP=WGWQ%NGVECTOR
      IF (WGWQ%LGAMMA) NP=NP*2

      ! clean head and wings
      IF (CHI0(NQ_COUNTER)%LREALSTORE) THEN
         ! all ranks clean their distributed row wing
         CHI0(NQ_COUNTER)%RESPONSER(1:2,:,:)= 0 
         ! master rank cleans its column wing
         IF ( NODE_ME == 1 ) THEN
            CHI0(NQ_COUNTER)%RESPONSER(:,1:2,:)= 0 
            CHI0(NQ_COUNTER)%RESPONSER(1,1,:)= 1 
         ENDIF 
         NP1 = SIZE( CHI0(NQ_COUNTER)%RESPONSER,1) 
         NP2 = SIZE( CHI0(NQ_COUNTER)%RESPONSER,2) 
      ELSE
         ! all ranks clean their distributed row wing
         CHI0(NQ_COUNTER)%RESPONSEFUN(1,:,:)= 0 
         ! master rank cleans its column wing
         IF ( NODE_ME == 1 ) THEN
            CHI0(NQ_COUNTER)%RESPONSEFUN(:,1,:)= 0 
            CHI0(NQ_COUNTER)%RESPONSEFUN(1,1,:)= 1 
         ENDIF 
         NP1 = SIZE( CHI0(NQ_COUNTER)%RESPONSEFUN,1) 
         NP2 = SIZE( CHI0(NQ_COUNTER)%RESPONSEFUN,2) 
      ENDIF

      ! initialize ScaLAPACK descriptor ( column major ) 
      CALL RESPONSE_DES_INIT_COLMAJ( RDES1D, IMAG_GRIDS%B, GDES )
      ! initialize ScaLAPACK descriptor ( block cyclic )
      CALL RESPONSE_DES_INIT_BC( RDES2D, IMAG_GRIDS%B, GDES )

      !allocate storage for redistributed matrix
      ALLOCATE(CHI_INV(RDES2D%MY_NROWS*RDES2D%MY_NCOLS), STAT = ISTAT) 
      IF ( ISTAT/=0 ) THEN
         CALL VTUTOR%ERROR( "OEP_POTENTIAL is not able to allocate "//&
           str( wsgf*RDES2D%MY_NROWS*RDES2D%MY_NCOLS) //&
           " kB of data on MPI rank 0." )
      ENDIF
      CALL REGISTER_ALLOCATE(wsgf*SIZE(CHI_INV,KIND=qi8), "OEP_chi")
      CHI_INV = 0 
      !
      ! redistribute matrix for ScaLAPACK 
      !
      !now call the BLACS redistribution routine 
      IF ( CHI0(NQ_COUNTER)%LREALSTORE ) THEN
         CALL PDGEMR2D(RDES1D%NROWS, RDES1D%NROWS, &
            CHI0(NQ_COUNTER)%RESPONSER(:,:,1), 1, 1, RDES1D%DESC, &
            CHI_INV, 1, 1, RDES2D%DESC, RDES2D%ICTXT)
         ! force Chi to be hermitean, use pdgemm routine
         CALL HERMITIZE_CHI( CHI_INV, NP1, RDES2D) 
      ELSE
         CALL PZGEMR2D(RDES1D%NROWS, RDES1D%NROWS, &
            CHI0(NQ_COUNTER)%RESPONSEFUN(:,:,1), 1, 1, RDES1D%DESC, &
            CHI_INV, 1, 1, RDES2D%DESC, RDES2D%ICTXT)
         ! force Chi to be hermitean, use pzgemm routine
         CALL HERMITIZE_CHI( CHI_INV, NP1, RDES2D) 
      ENDIF
#ifdef unit_test
      !
      ! test, redistribute back to original layout
      !
      IF ( CHI0(NQ_COUNTER)%LREALSTORE ) THEN
         CALL PDGEMR2D(RDES2D%NROWS, RDES2D%NROWS, &
            CHI_INV, 1, 1, RDES2D%DESC, &
            CHI0(NQ_COUNTER)%RESPONSER(:,:,1), 1, 1, RDES1D%DESC, RDES1D%ICTXT)
      ELSE
         CALL PZGEMR2D(RDES2D%NROWS, RDES2D%NROWS, &
            CHI_INV, 1, 1, RDES2D%DESC, &
            CHI0(NQ_COUNTER)%RESPONSEFUN(1,1,1), 1, 1, RDES1D%DESC, RDES1D%ICTXT)
      ENDIF
      DO J = 1 , NP2
      DO I = 1 , NP1
        WRITE(200+WDES%COMM%NODE_ME,'(2I4,2F20.10)')I,J,CHI0(NQ_COUNTER)%RESPONSEFUN(I,J,1)
      ENDDO
      ENDDO
      
      
      !CALL DUMP_GREENS_HAM_GDEF( "CHI-RESHAPED", RESHAPE( CHI0(NQ_COUNTER)%RESPONSEFUN(:,:,1), (/GDES%RES_NRPLWV_ROW_DATA_POINTS*GDES%RES_NRPLWV_COL_DATA_POINTS/)), NP, RDES1D, IO%IU0)
      CALL DUMP_GREENS_HAM_GDEF("CHI-hermitized distributed", CHI_INV, NP, RDES2D, IO%IU0 )
      ALLOCATE( CHI_TMP( SIZE( CHI_INV) ) )
      ALLOCATE( CHI2( SIZE( CHI_INV) ) )
      CHI_TMP = CHI_INV
      CHI2 = CHI_INV
      CALL DUMP_GREENS_HAM_GDEF("X distributed", CHI_INV, NP, RDES2D, IO%IU0 )
#endif 
      !
      ! call diagonalizer for soft inversion, for finite temperature
      ! 
      CALL PSEUDO_INVERSE_DESC( CHI_INV,NP,RDES2D%DESC, 1.E-4_q, IO%IU6 )

#ifdef unit_test
      CALL DUMP_GREENS_HAM_GDEF("V+ 1/eig V", CHI_INV, NP, RDES2D, IO%IU0 )
      !
      ! test A * A ^-1 
      !
      CALL PGGEMM( 'N', 'N', NP, NP, NP, one,  &
         CHI_TMP(1), 1, 1, RDES2D%DESC,           &
         CHI_INV(1), 1, 1, RDES2D%DESC,        &
         zero, CHI2(1), 1, 1, RDES2D%DESC )
      CALL DUMP_GREENS_HAM_GDEF(" X^-1 . X", CHI2, NP, RDES2D, IO%IU0 )
      DEALLOCATE( CHI_TMP )
      DEALLOCATE( CHI2 )
#endif 
#undef unit_test
      !
      ! redistribute matrix to column major distribution
      !
      !now call the BLACS redistribution routine 
      IF ( CHI0(NQ_COUNTER)%LREALSTORE ) THEN
         CALL PDGEMR2D(RDES2D%NROWS, RDES2D%NROWS, &
            CHI_INV, 1, 1, RDES2D%DESC, &
            CHI0(NQ_COUNTER)%RESPONSER(:,:,1), 1, 1, RDES1D%DESC, RDES1D%ICTXT)
      ELSE
         CALL PZGEMR2D(RDES2D%NROWS, RDES2D%NROWS, &
            CHI_INV, 1, 1, RDES2D%DESC, &
            CHI0(NQ_COUNTER)%RESPONSEFUN(:,:,1), 1, 1, RDES1D%DESC, RDES1D%ICTXT)
      ENDIF
      CALL DEREGISTER_ALLOCATE(wsgf*SIZE(CHI_INV,KIND=qi8), "OEP_chi")
      DEALLOCATE(CHI_INV)
      NULLIFY(CHI_INV)
      ! get rid of block-cyclic ScaLAPACK descriptor 
      CALL DESTROY_GDES_MAT( RDES2D )
      !
      ! at this stage CHI0(NQ_COUNTER)%RESPONSE[R|FUN] stores inverse CHI at Gamma
      ! get rid of temporary storage for computation of inverse chi
      !
      ! 
      ! now, determine V_xc - V_xc^DFT:
      ! first step is to determine "change in charge" from W
      ! this is the r.h.s. of the third Eq. given in the description 
      ! of this subroutine and is determined from natural orbitals stored in W
      ! 
      ! wavefunction grid (GRID) -> charge density grid (GRIDC)
      CALL SET_CHARGE(W, WDES, INFO%LOVERL, &
           GRID, GRIDC, GRID_SOFT, GRIDUS, C_TO_US, SOFT_TO_C, &
           LATT_CUR, P, SYMM, T_INFO, &
           CHDEN2, LMDIM, CRHODE2, CHTOT2, RHOLM2, N_MIX_PAW, IRDMAX)
      ! at this stage CHTOT2 is the total OEP charge on reciprocal space grid
      ! described by charge density grid handle GRIDC (C ... CHARGE)
      !
      ! remark on GRID, GRID_SOFT:
      ! GRID is the grid on which orbitals are stored.
      ! arrays stored on GRID are always complex valued
      ! GRID_SOFT is the same as GRID with the exception that 
      ! arrays stored on GRID_SOFT are real-valued (like densities)
      !
      ! now downsample CHTOT2 to GRID_SOFT and store result in SVOEP
      CALL SET_SV( W%WDES%GRID, GRIDC, GRID_SOFT, W%WDES%COMM_INTER, SOFT_TO_C,&
         W%WDES%NCDIJ, SV1, CHTOT2)
      ! SV1 stores the OEP charge on GRID_SOFT represented in real-space
      ! now extract non-zero modes of SV1 in reciprocal space:
      ! this requires, first, to restore SV1 into a complex array.
      ! Second, forward transform to reciprocal space and impose cutoff.
      ! The result is stored in CHG and is a relatively small array
      ! of the same size as the number of columns of response function 
      ! stored locally. At this stage this number is determined by NP1.
      ALLOCATE( CHG(NP1, WDES%NCDIJ ) ) 
      CHG = 0 
      ! helper array 
      ALLOCATE( CV1(WGWQ%GRID%MPLWV ) )
      DO ISP = 1, WDES%NCDIJ 
         CV1 = 0 
         CALL RL_COPY_RGRID_WAVE(SV1(1,ISP), CV1(1), 1.0_q, GRID)
         ! second, perform a forward transform and impose cutoff 
         ! at the same time. This neglects high-frequency modes
         CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), &
              CV1(1), CHG(1,ISP), WGWQ%GRID, .FALSE.)
         CHG(:,ISP)=CHG(:,ISP)*(1.0_q/WGWQ%GRID%NPLWV)
      ENDDO
      DEALLOCATE( CV1 )


      ! CHG stores the OEP charge in reciprocal space  (r.h.s. of Eq.)
      ! and has the same dimension as NP
      ! also, all ranks have the same copy of CHG
      ! symmetrize this charge again
      CALL SYMMETRIZE_CHARGE( CHG, SV1, WDES, WGWQ, GRID, GRID_SOFT )
!WRITE(200+WDES%COMM%NODE_ME, '(F20.10)')CHG
      ! 
      ! now, obtain V_xc - V_xc^DFT = D rho * chi^-1
      !
      ! allocate some additional storage 
      NP2_GLB = (RDES1D%COMM_INTRA%NODE_ME-1) * GDES%RES_NRPLWV_COL_DATA_POINTS + 1 
      NP2_END =  RDES1D%COMM_INTRA%NODE_ME    * GDES%RES_NRPLWV_COL_DATA_POINTS
      ALLOCATE( CHG2( NP2 ) )

      IF (CHI0(NQ_COUNTER)%LREAL) THEN
      !
      ! CHI0 is real-valued, but CHG is a half sized complex valued array
      !
         ALLOCATE( CHG1( NP1 ) )
         DO ISP = 1, WDES%NCDIJ
            CHG1 = 0 
            CHG2 = 0 
            ! restride complex array to real array of double size
            CALL DAXPY( NP, 1.0_q, CHG(1,ISP), 1,  CHG1(1), 1)
            ! evaluate vector * matrix multiplication
            !store stripe of vector 
            CHG2(:) = CHG1(NP2_GLB:NP2_END)
            ! CHI * CHG, CHI is blocked
            CALL DGEMV('N', NP, NP2, 1._q, &
                 CHI0(NQ_COUNTER)%RESPONSER(1,1,1), NP1, &
                 CHG2(1), 1, &
                 0.0_q, CHG1(1), 1 )
            CALLMPI( M_sum_g( RDES1D%COMM_INTRA, CHG1(1), NP ) )
            ! restride real array into complex array of half size 
            CALL DAXPY( NP, 1.0_q, CHG1(1), 1,  CHG(1,ISP), 1)
         ENDDO
         DEALLOCATE( CHG1 )
      !
      ! full k-point version
      !
      ELSE
         DO ISP = 1, WDES%NCDIJ
!           ! CHG * CHI^-1 = CHG_OEP distributed among COMM_INTRA
            ! example call for ?GEMM:
            ! dgemm('n','n', m, n, k, 1.d+0, a, m, b, k, 0.d+0, c, m)
            !  ! a_mk B_kn = C_mn
            ! however, use ?GEMV
            ! matrix vector multiplication, blocked 
            ! CHI * CHG , CHI is blocked 
            !store stripe of vector 
            CHG2(:) = CHG(NP2_GLB:NP2_END,ISP)
            ! matrix vector multiplication 
            CALL ZGEMV('N', NP, NP2, (1._q, 0._q), &
                 CHI0(NQ_COUNTER)%RESPONSEFUN(1,1,1), NP1, &
                 CHG2(1), 1, &
                 (0.0_q,0._q), CHG(1,ISP), 1 )
            ! synchronize data
            CALLMPI( M_sum_z( RDES1D%COMM_INTRA, CHG(1,ISP), NP ) )
         ENDDO
      ENDIF
      DEALLOCATE( CHG2 ) 
      ! ESHIFT shifts OEP potential to specific value in real-space at infinity 
      ! this value can be chosen freely
!      OEP%ESHIFT=0
!WRITE(*,*) "WARNING SETTING ESHIFT = ", OEP%ESHIFT
      ! 
      ! set value of potential at infinity
      ! which is basically a shift of OEP eigenenergies
      CHG(1,:)=OEP%ESHIFT
!WRITE(300+WDES%COMM%NODE_ME, '(F20.10)')CHG
      ! symmetrize potential and store in SV
      CALL SYMMETRIZE_CHARGE( CHG, SV1, WDES, WGWQ, GRID, GRID_SOFT, SV)
      ! also multiply stored SV with factor AMIX to dampen response
      SV=SV*AMIX
      ! at this stage, OEP corrected xc-potential is stored in SV on soft grid
      ! non-synchronized SV is stored in SV
!WRITE(400+WDES%COMM%NODE_ME, '(F20.10)')CHG
!WRITE(500+WDES%COMM%NODE_ME, '(F20.10)')SV
      DEALLOCATE( CHG ) 
!=======================================================================
! now update OEP handle, that is add corrections to potential
!=======================================================================
      ! restore original potential: V_ext + V_h + V^{KS}_xc
      CDIJ = OEP%CDIJ
      CVTOT = OEP%CVTOT
      IF (IO%IU0>=0) WRITE(*,'(8F14.7)') CVTOT(1:8,1)

#ifdef verbose
WRITE(100+WDES%COMM%NODE_ME,'(F20.10)')SV
WRITE(200+WDES%COMM%NODE_ME,'(F20.10)')CDIJ
WRITE(300+WDES%COMM%NODE_ME,'(F20.10)')CVTOT
#endif

      ! first determine the one-center contributions
      OEP%CDIJ=CDIJ
      ! determine CDIJ from from original DFT potential CVTOT
      CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
           LMDIM,OEP%CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)
      ! subtract the plane wave part from the CDIJ that were read from the file 
      ! to obtain the OEP contribution only
      OEP%CDIJ=CDIJ-OEP%CDIJ
      !
      ! OEP%CDIJ = C_old - C_new
      ! add correction SV to CVTOT
      DO ISP=1,WDES%NCDIJ
         CALL FFT_RC_SCALE(CVTOT(1,ISP),CVTOT(1,ISP),GRIDC)
! GR  ID_SOFT -> GRIDC , SOFT_TO_C handle contains info about mapping 
         CALL ADD_GRID(GRIDC, GRID_SOFT, SOFT_TO_C, SV(1,ISP), CVTOT(1,ISP))
      ENDDO
      
      ! now set SV from CVTOT
      CALL SET_SV( GRID, GRIDC, GRID_SOFT, WDES%COMM_INTER, SOFT_TO_C, WDES%NCDIJ, SV, CVTOT)
      IF (IO%IU0>=0) WRITE(*,'(8F14.7)') CVTOT(1:8,1)

      CALL SETDIJ(WDES,GRIDC,GRIDUS,C_TO_US,LATT_CUR,P,T_INFO,INFO%LOVERL, &
          LMDIM,CDIJ,CQIJ,CVTOT,IRDMAA,IRDMAX)

      ! add one center OEP contributions
      CDIJ = CDIJ + OEP%CDIJ

      OEP%CDIJ = CDIJ
      OEP%CVTOT = CVTOT
      OEP%SV = SV 
!=======================================================================
! write the local potential
!=======================================================================
      CALL DUMP_POTENTIAL_TO_FILE( OEP, 'OEP.vasp' )
      ! no need to keep this 
      CALL DESTROY_GDES_MAT( RDES1D )
    
      PROFILING_STOP( 'OEP_POTENTIAL' )

      CONTAINS
 
      ! performs a hermitization of chi0 using scalapack 
      SUBROUTINE HERMITIZE_CHI( CHI, N, RDES ) 
         GDEF                              :: CHI(:)
         INTEGER                           :: N
         TYPE(greens_mat_des), POINTER     :: RDES
         ! local 
         GDEF                        :: CWORK1(SIZE(CHI))
         GDEF                        :: CWORK2(SIZE(CHI))
         GDEF                        :: CDIAG( N ) 

         CWORK1 = 0
         CDIAG(:) = 1
#ifdef gammareal
         CALL ADD_TO_DIAGONALE_REAL(N, CWORK1, CDIAG, RDES%DESC )
#else
         CALL ADD_TO_DIAGONALE(N, CWORK1, CDIAG, RDES%DESC )
#endif
         CWORK2 = CHI
         CALL PGGEMM( 'C', 'N', N, N, N, (0.5_q,0.0_q),  &
              CWORK2(1), 1, 1, RDES%DESC, &
              CWORK1(1), 1, 1, RDES%DESC, &
              (0.5_q,0.0_q),  &
              CHI(1), 1, 1, RDES%DESC )

      END SUBROUTINE HERMITIZE_CHI

      !
      ! symmetrize the final potential
      ! since this can not be done using CHG, go to potential on GRID_SOFT
      ! symmetrize and go back
      !
      SUBROUTINE SYMMETRIZE_CHARGE( CHG, SV, WDES, WGWQ, GRID, GRID_SOFT, SVSTORE ) 
         USE main_mpi
         COMPLEX( q )    :: CHG(:,:)
         RGRID           :: SV(:,:)
         TYPE( wavedes ) :: WDES
         TYPE( wavedes1) :: WGWQ
         TYPE( grid_3d ) :: GRID
         TYPE( grid_3d ) :: GRID_SOFT
         RGRID, OPTIONAL :: SVSTORE(:,:)
         ! local
         COMPLEX(q)      :: CWORK(WDES%GRID%MPLWV)
         INTEGER         :: ISP

         DO ISP=1,WDES%ISPIN
            ! FFT CHG to reciprocal space 
            CALL FFTWAV(WGWQ%NGVECTOR, WGWQ%NINDPW(1), &
                 CWORK(1), CHG(1,ISP), WGWQ%GRID)
            ! copy wave function based array (COMPLEX) to RGRID
            CALL RL_COPY_WAVE_RGRID( CWORK(1), SV(1,ISP), 1.0_q, GRID)
         
            CALL FFT_RC_SCALE(SV(1,ISP), CWORK(1), GRID_SOFT)
            ! CWORK is given in reciprocal space 
            CALL RC_ADD(CWORK(1),1.0_q,CWORK(1),0.0_q,SV(1,ISP),GRID_SOFT)
         ENDDO
     
         ! change storage convention to (total, magnetization)
         CALL RC_FLIP(SV,GRID,WDES%NCDIJ,.FALSE.)
     
         IF (SYMM%ISYM==2) THEN
            ! hand brake for NONCOLLINEAR version pulled in chi_super
            DO ISP=1,WDES%ISPIN
               CALL RHOSYM(SV(1,ISP),GRID_SOFT,SYMM%PTRANS,T_INFO%NIOND,SYMM%MAGROT,ISP)
            ENDDO
         ENDIF
     
         ! change storage convention back
         CALL RC_FLIP(SV,GRID,WDES%NCDIJ,.TRUE.)

         DO ISP=1,WDES%ISPIN
            IF ( PRESENT( SVSTORE ) ) THEN
               ! store SV  (update potential)
               SVSTORE(:,ISP)=SV(:,ISP)
            ENDIF
            ! go to real space
            CALL FFT3D( SV(1,ISP), GRID_SOFT, 1 )
            ! correct data on all nodes
#ifdef realmode
            CALLMPI( M_bcast_d(COMM_INTER, SV(1,ISP), GRID%RL%NP))
            CALLMPI( M_bcast_d(COMM_KINTER, SV(1,ISP), GRID%RL%NP))
#else
            CALLMPI( M_bcast_z(COMM_INTER, SV(1,ISP), GRID%RL%NP))
            CALLMPI( M_bcast_z(COMM_KINTER, SV(1,ISP), GRID%RL%NP))
#endif
            ! go back to CHG
            CALL RL_COPY_RGRID_WAVE(SV(1,ISP), CWORK(1), 1.0_q, GRID)
            ! go to reciprocal space 
            CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), &
                 CWORK(1), CHG(1,ISP), WGWQ%GRID, .FALSE.)
            CHG(:,ISP)=CHG(:,ISP)*(1.0_q/WGWQ%GRID%NPLWV)
         ENDDO
      END SUBROUTINE SYMMETRIZE_CHARGE
    
      !
      ! write local potential stored in OEP to file 
      !
      SUBROUTINE DUMP_POTENTIAL_TO_FILE( OEP, FILEN )
         USE paw, ONLY: SET_RHO_PAW
         USE pawm, ONLY: WRT_RHO_PAW
         TYPE( oep_handle ) :: OEP  
         CHARACTER( LEN=* ), OPTIONAL :: FILEN
         ! local 
         INTEGER            :: NODE_ME = 1 
         REAL(q)            :: DLM_EXX(N_MIX_PAW,WDES%NCDIJ)

#ifdef MPI   
         NODE_ME = WDES%COMM_KINTER%NODE_ME
#endif

         IF (IO%LVTOT .AND. LOEP .AND. NODE_ME == 1) THEN
            IF (IO%IU0>=0) THEN
               IF ( PRESENT( FILEN ) ) THEN
                  IF (IO%LOPEN) OPEN(IO%IUVTOT,FILE=FILEN,STATUS='UNKNOWN')
               ELSE
                  IF (IO%LOPEN) OPEN(IO%IUVTOT,FILE='POT',STATUS='UNKNOWN')
               ENDIF
               REWIND IO%IUVTOT
               CALL OUTPOS(IO%IUVTOT,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A,.FALSE.,T_INFO%POSION)
            ENDIF
            
            ! at the moment the spin up and down potential is written to the file
            CALL SET_RHO_PAW(WDES, P, T_INFO, INFO%LOVERL, WDES%NCDIJ, LMDIM, &
                 OEP%CDIJ, DLM_EXX)
            CALL OUTPOT(GRIDC, IO%IUVTOT,.TRUE.,OEP%CVTOT)
            CALL WRT_RHO_PAW(P, T_INFO, INFO%LOVERL, DLM_EXX(:,1), GRIDC%COMM, IO%IUVTOT)
            
            DO ISP=2,WDES%NCDIJ
               IF (IO%IU0>=0) THEN
                  WRITE( IO%IUVTOT,'(5E20.12)') (T_INFO%ATOMOM(I),I=1,T_INFO%NIONS)
               ENDIF
               CALL OUTPOT(GRIDC, IO%IUVTOT,.TRUE.,OEP%CVTOT(1,ISP))
               CALL WRT_RHO_PAW(P, T_INFO, INFO%LOVERL, DLM_EXX(:,ISP), GRIDC%COMM, IO%IUVTOT )
            ENDDO
            ! close unit
            IF( IO%IU0>=0 .AND. IO%LOPEN ) CLOSE( IO%IUVTOT) 
         ENDIF
      END SUBROUTINE DUMP_POTENTIAL_TO_FILE

      ! testing routine, dump charge given in reciprocal space 
      ! so CHG is complex valued
      SUBROUTINE DUMP_CHARGE_GRID( FILEN, GRID, CHG)
         USE fileio, ONLY: OUTCHG
         CHARACTER(LEN=*) :: FILEN
         TYPE (grid_3d)   :: GRID
         COMPLEX(q)       :: CHG(:,:)
         !local 
         INTEGER          :: ISP
         ! print this charge to file and compare with charge on fine grid
         IF( IO%IU0>=0 ) THEN
            OPEN(UNIT=99,FILE=FILEN,STATUS='UNKNOWN')
            CALL OUTPOS(99,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A, &
           &             .FALSE., T_INFO%POSION)
         ENDIF
         DO ISP=1,SIZE(CHG,2)
            CALL OUTCHG(GRID,99,.TRUE.,CHG(:,ISP))
         ENDDO
         IF( IO%IU0>=0 ) CLOSE(99)
      END SUBROUTINE DUMP_CHARGE_GRID

      ! helper, dumps charge density given on RGRID 
      ! so CHG is real valued
      SUBROUTINE DUMP_CHARGE_RGRID( FILEN, GRID, CHG)
         CHARACTER(LEN=*) :: FILEN
         TYPE (grid_3d)   :: GRID
         RGRID            :: CHG(:,:)
         !local 
         INTEGER          :: ISP
         ! print this charge to file and compare with charge on fine grid
         IF( IO%IU0>=0 ) THEN
            OPEN(UNIT=99,FILE=FILEN,STATUS='UNKNOWN')
            CALL OUTPOS(99,.FALSE.,INFO%SZNAM1,T_INFO,LATT_CUR%SCALE,LATT_CUR%A, &
           &             .FALSE., T_INFO%POSION)
         ENDIF
         ! dump charge, which is already given in real space representation
         DO ISP=1,SIZE(CHG,2)
            CALL OUTCHG_RGRID(GRID,99,.TRUE.,CHG(:,ISP))
         ENDDO
         IF( IO%IU0>=0 ) CLOSE(99)
      END SUBROUTINE DUMP_CHARGE_RGRID
   END SUBROUTINE OEP_POTENTIAL

   !
   ! simple copy from OUTCHG, however
   ! no FFT to real space is done, so CHTOT must be provided on 
   ! real-space grid 
   !
   SUBROUTINE OUTCHG_RGRID(GRID, IU, LLONG,SV)
      USE prec
      USE mpimy
      USE mgrid
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      TYPE (grid_3d) GRID

      RGRID  SV(DIMREAL(GRID%MPLWV))
      LOGICAL LLONG        ! long or short format
! work arrays
      REAL(q),ALLOCATABLE ::  WORK (:)
      INTEGER NALLOC,NZ, NWRITE, NWRITTEN
      CHARACTER (40) FORM
      INTEGER ISTAT, NODE_ME, IONODE, IU, N, IERR

      NALLOC=GRID%NGX*GRID%NGY

      ALLOCATE(WORK(NALLOC),STAT=ISTAT)
      IF (ISTAT>0) RETURN ! can not write the charge immediate exit

      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=GRID%COMM%NODE_ME
      IONODE =GRID%COMM%IONODE
#endif

      IF (GRID%NPLWV/= GRID%NGX*GRID%NGY*GRID%NGZ) THEN
        CALL vtutor%bug("OUTCHG_RGRID NPLWV is not compatibel with NGX,NGY,NGZ \n " // str(GRID%NPLWV) // " " // str(GRID%NGX) // " " // str(GRID%NGY) // " " // str(GRID%NGZ), __FILE__, __LINE__ )
      ENDIF
      ! original subroutine stores CHTOT (complex) to SV (real)
      ! and transforms SV to real-space 
      ! note that SV has 2 times the dimension of CHTOT
      ! so that Re and Im part of CHTOT fits into SV 
      ! this requires that RC_ADD to reside outside of a module 
      ! for the fortran compiler to understand this
!      ! transfer to SV and FFT
!      CALL RC_ADD(CHTOT,1.0_q,CHTOT,0.0_q,SV,GRID)
!      ! transform to real space representation
!      CALL FFT3D(SV,GRID,1)

      IF (LLONG) THEN
        FORM='(1(1X,E17.11))'
        NWRITE=5
      ELSE
        FORM='(1(1X,G11.5))'
        NWRITE=10
      ENDIF

      do_io WRITE(IU,'(3I5)') GRID%NGX,GRID%NGY,GRID%NGZ

      NWRITTEN=0
      DO NZ=1,GRID%NGZ
         CALL MRG_GRID_RL_PLANE(GRID, WORK, SV, NZ)
         io_begin
         DO N=1,NALLOC
            NWRITTEN=NWRITTEN+1
            IF ( MOD(NWRITTEN,NWRITE)==0 ) THEN
               WRITE(IU,FORM) WORK(N)
            ELSE
               WRITE(IU,FORM,ADVANCE='NO') WORK(N)
            ENDIF
         ENDDO
         io_end
      ENDDO
      IF ( MOD(NWRITTEN,NWRITE)/=0 ) WRITE(IU,*)' '

      DEALLOCATE(WORK)

      RETURN
   END SUBROUTINE OUTCHG_RGRID
#undef debug

!
! destroys distributed response descriptor
!
   SUBROUTINE DESTROY_GDES_MAT( RDES )
      TYPE(greens_mat_des), POINTER :: RDES 
 
      IF ( .NOT. ASSOCIATED( RDES ) ) RETURN

      CALL BLACS_GRIDEXIT( RDES%ICTXT ) 
      DEALLOCATE( RDES ) 
      NULLIFY( RDES )
   END SUBROUTINE DESTROY_GDES_MAT

!
!  initalize column major distrubuted scalapack matrix descriptor
!
   SUBROUTINE RESPONSE_DES_INIT_COLMAJ( RDES, LOOPD, GDES )
      TYPE(greens_mat_des), POINTER :: RDES ! distributed response
      TYPE(loop_des)   :: LOOPD
      TYPE(greensfdes) :: GDES
      ! local 
      INTEGER        :: NPROW0 
      INTEGER        :: NPCOL0
      INTEGER        :: MYROW0
      INTEGER        :: MYCOL0
      INTEGER        :: LLD
      INTEGER        :: INFO 
      INTEGER        :: NROWS
      INTEGER        :: NCOLS
      INTEGER, EXTERNAL :: NUMROC

      ! mark descriptor as initialized 
      ALLOCATE( RDES )

      ! set communicators
      RDES%COMM_INTRA=LOOPD%COMM_IN_GROUP
      RDES%COMM_INTER=LOOPD%COMM_BETWEEN_GROUPS
      RDES%COMM=LOOPD%COMM

      !---------------------------------------------------------------------------------------
      !initialzie 1D column-major processor grid
      ! row major ordering (but column-major ordering is identical if there is only 1 column)
      NPROW0=1
      NPCOL0=LOOPD%COMM_IN_GROUP%NCPU
      ! actual size of local matrix
      NROWS  = GDES%RES_NRPLWV_ROW_DATA_POINTS
      !number of cols ( in parallel blocking factor )
      NCOLS  = GDES%RES_NRPLWV_COL_DATA_POINTS
      ! 
      RDES%NROWS = NROWS

      CALL PROCMAP( RDES%COMM_INTRA, RDES%ICTXT, 2, NPROW0, NPCOL0)
      CALL BLACS_GRIDINFO( RDES%ICTXT, NPROW0, NPCOL0, MYROW0, MYCOL0) !get the ids for scaLAPACK
      !define corresponding array descritpor DESCA
      !since CHI_WORK is distributed blockwisely,  
      !the leading dimension of rows is simply NROWS, 
      !check this 
      LLD = MAX( 1, NUMROC(NROWS, NCOLS, MYROW0, 0, NPROW0) )
      IF ( LLD /= NROWS ) THEN 
         CALL vtutor%bug("RESPONSE_DES_INIT_COLMAJ LLD has wrong size for NODE " // str(RDES%COMM%NODE_ME) // " " // str(LLD) // " " // str(NROWS), __FILE__, __LINE__)
      ENDIF 
      CALL DESCINIT(RDES%DESC, NROWS , NROWS, NCOLS, NCOLS, 0, 0, RDES%ICTXT, &
         LLD, INFO )
      !Eventually take care for errors 
      IF ( INFO /= 0 ) THEN 
         CALL vtutor%bug("Error in 1D col-major initialization " // str(INFO), __FILE__,__LINE__)
      ENDIF

      !check if CHI_WORK has proper size
      RDES%MY_NROWS = NUMROC(RDES%DESC(M_),RDES%DESC(MB_),MYROW0,0,NPROW0)
      RDES%MY_NCOLS = NUMROC(RDES%DESC(N_),RDES%DESC(NB_),MYCOL0,0,NPCOL0) 

   END SUBROUTINE RESPONSE_DES_INIT_COLMAJ

!
!  initalize block cyclic distrubuted scalapack matrix descriptor
!
   SUBROUTINE RESPONSE_DES_INIT_BC( RDES, LOOPD, GDES )
      USE acfdt_GG, ONLY: MB, NB
      TYPE(greens_mat_des), POINTER :: RDES ! distributed response
      TYPE(loop_des)   :: LOOPD
      TYPE(greensfdes) :: GDES
      ! local 
      INTEGER        :: NPROW
      INTEGER        :: NPCOL
      INTEGER        :: MYROW
      INTEGER        :: MYCOL
      INTEGER        :: LLD
      INTEGER        :: INFO 
      INTEGER        :: NROWS
      INTEGER        :: NCOLS
      INTEGER, EXTERNAL :: NUMROC
      
      ! mark descriptor as initialized 
      ALLOCATE( RDES )

      ! set communicators
      RDES%COMM_INTRA=LOOPD%COMM_IN_GROUP
      RDES%COMM_INTER=LOOPD%COMM_BETWEEN_GROUPS
      RDES%COMM=LOOPD%COMM

      !---------------------------------------------------------------------------------------
      !initialzie block-cyclic 2D processor grid
      CALL FERMAT_RAZOR( RDES%COMM_INTRA%NCPU, NPROW, NPCOL )
      ! actual size of local matrix
      NROWS  = GDES%RES_NRPLWV_ROW_DATA_POINTS
      !number of cols ( in parallel blocking factor )
      NCOLS  = GDES%RES_NRPLWV_COL_DATA_POINTS
      RDES%NROWS = NROWS

      CALL PROCMAP( RDES%COMM_INTRA, RDES%ICTXT, 1, NPROW, NPCOL)
      CALL BLACS_GRIDINFO( RDES%ICTXT, NPROW, NPCOL, MYROW, MYCOL) !get the ids for scaLAPACK
      !define corresponding array descritpor DESCA
      !since CHI_WORK is distributed blockwisely,  
      !the leading dimension of rows is simply NROWS, 
      !check this 
      !leading dimension of local array for block cyclic distribution
      LLD = MAX( 1 , NUMROC(NROWS, MB, MYROW, 0, NPROW) )
      CALL DESCINIT(RDES%DESC, NROWS , NROWS, MB, NB, 0, 0, RDES%ICTXT, LLD, INFO )
      !Eventually take care for errors 
      IF ( INFO /= 0 ) THEN 
         CALL vtutor%error("Error in 2D block-cyclic initialization " // str(INFO))
      ENDIF

      !check if CHI_WORK has proper size
      RDES%MY_NROWS = NUMROC(RDES%DESC(M_),RDES%DESC(MB_),MYROW,0,NPROW)
      RDES%MY_NCOLS = NUMROC(RDES%DESC(N_),RDES%DESC(NB_),MYCOL,0,NPCOL) 
   END SUBROUTINE RESPONSE_DES_INIT_BC

!*********************************************************************
! matrix size checker 
! TODO: change blocking factor NB to support more NCPUs
! but QUERRY_SCALA_NB sets NB to the value in scala.F
!*********************************************************************
   SUBROUTINE CHECK_GDES_MAT_SIZE( COMM, NGROUPS, NB_TOT, IO ) 
      USE base, ONLY: in_struct
      USE scala, ONLY: FERMAT_RAZOR
      USE tutor, ONLY: vtutor, GDESMatBlockRed, GDESMatSize, isAlert, isError, argument
      TYPE( communic )  :: COMM
      INTEGER           :: NGROUPS
      INTEGER           :: NB_TOT
      TYPE( in_struct ) :: IO
      ! local
      INTEGER           :: NCPU
      INTEGER           :: NODE_ME
      INTEGER           :: NPROCS
      INTEGER           :: I,J
      INTEGER           :: NPROW, NPCOL
      TYPE(argument)    :: arg

      PROFILING_START( 'CHECK_GDES_MAT_SIZE' )

      NCPU = COMM%NCPU
      NODE_ME = COMM%NODE_ME
 
      ! estimate maximum number of CPUs among each tau group
      IF ( MOD( NCPU, NGROUPS ) ==  0 ) THEN
         NPROCS = NCPU/NGROUPS 
      ELSE
         NPROCS = NCPU/NGROUPS + 1 
      ENDIF   

      CALL FERMAT_RAZOR( NPROCS,  NPROW, NPCOL )
     
      ! check if any cpus in group will not carry data due to blocking 
      ! factor and NGROUPS chosen
      DO I = 1, NPCOL - 1 
         IF ( I*NB_SCALA > NB_TOT ) EXIT
      ENDDO

      DO J = 1, NPROW - 1 
         IF ( J*NB_SCALA > NB_TOT ) EXIT
      ENDDO
       
!      ! Blocking size is too large, decrease it and drop error for user
!      IF ( I < NPCOL .OR. J < NPROW ) THEN
!         ALLOCATE(arg%ival(5))
!         NB_SCALA = NB_SCALA/2 
!         arg%ival(1)=NPROW
!         arg%ival(2)=NPCOL
!         arg%ival(3)=NB_SCALA*2
!         arg%ival(4)=NB_TOT
!         arg%ival(5)=NGROUPS
!         ! if this works out, drop a tutor warning
!         CALL vtutor%write(isAlert, GDESMatBlockRed, arg)
!
!      ENDIF
! 
!      ! try again with smaller block size 
!
!      ! check if any cpus in group will not carry data due to blocking 
!      ! factor and NGROUPS chosen
!      DO I = 1, NPCOL - 1 
!         IF ( I*NB_SCALA > NB_TOT ) EXIT
!      ENDDO
!
!      DO J = 1, NPROW - 1 
!         IF ( J*NB_SCALA > NB_TOT ) EXIT
!      ENDDO

      ! if at least one CPUs will not have data, print error
      IF ( I < NPCOL .OR. J < NPROW  ) THEN
         ALLOCATE(arg%ival(5))
         arg%ival(1)=NPROW
         arg%ival(2)=NPCOL
         arg%ival(3)=NB_SCALA
         arg%ival(4)=NB_TOT
         arg%ival(5)=NGROUPS
         CALL vtutor%write(isError, GDESMatSize, arg)
      ENDIF

      PROFILING_STOP( 'CHECK_GDES_MAT_SIZE' )
   END SUBROUTINE CHECK_GDES_MAT_SIZE

!*********************************************************************
!
!
! initialize descriptor 
!
!
!*********************************************************************

  SUBROUTINE GREENS_MAT_DES_INIT( GDES_MAT, IGRID_DES, NB_TOT)

    USE scala
    TYPE (greens_mat_des), POINTER :: GDES_MAT   ! out: pointer to generated descriptor
    TYPE (loop_des)                :: IGRID_DES
    INTEGER                        :: NB_TOT     ! in: maximum number of orbitals
  ! local
    INTEGER :: NPROCS, NPROW, NPCOL, MYROW, MYCOL, NP, NQ
    INTEGER,EXTERNAL :: NUMROC
    INTEGER :: IERROR

    ALLOCATE(GDES_MAT)
    
    GDES_MAT%COMM=IGRID_DES%COMM
    GDES_MAT%COMM_INTRA=IGRID_DES%COMM_IN_GROUP
    GDES_MAT%COMM_INTER=IGRID_DES%COMM_BETWEEN_GROUPS

    ! few sanity tests
    IF (.NOT. LscaLAPACK) THEN
       CALL vtutor%bug("GREENS_MAT_DES_INIT requires scaLAPACK", __FILE__, __LINE__)
    ENDIF

    IF (GDES_MAT%COMM%NCPU /= GDES_MAT%COMM_INTRA%NCPU * GDES_MAT%COMM_INTER%NCPU) THEN
       CALL vtutor%bug("GREENS_MAT_DES_INIT core number incompatible " // str(GDES_MAT%COMM%NCPU) // " " // str(GDES_MAT%COMM_INTRA%NCPU) // " " // str(GDES_MAT%COMM_INTER%NCPU), __FILE__, __LINE__)
    ENDIF

    GDES_MAT%NROWS=NB_TOT
    GDES_MAT%ICTXT=0
    NPROCS=GDES_MAT%COMM_INTRA%NCPU
!#define distribute_columns
#ifdef distribute_columns    
    ! we chose to distribute the columns in a round robin fashion
    ! at least for the time being
    NPROW=1
    NPCOL=NPROCS
    
    CALL MPI_barrier( GDES_MAT%COMM_INTRA%MPI_COMM, IERROR )
    CALL PROCMAP( GDES_MAT%COMM_INTRA, GDES_MAT%ICTXT, 2, NPROW, NPCOL)

! calculate local size of matrices

    CALL BLACS_GRIDINFO( GDES_MAT%ICTXT, NPROW, NPCOL, MYROW, MYCOL )

    NP = NUMROC(NB_TOT,NB_TOT,MYROW,0,NPROW)   ! get number of rows on proc
    NQ = NUMROC(NB_TOT,1, MYCOL,0, NPCOL)      ! get number of cols on proc

    CALL DESCINIT(GDES_MAT%DESC, NB_TOT,NB_TOT, NB_TOT,1, 0,0, GDES_MAT%ICTXT,MAX( 1, NP ), IERROR ) ! setup descriptor
    IF (IERROR.NE.0) THEN
      CALL vtutor%bug("GREENS_MAT_DES_INIT: DESC, DESCINIT, INFO: " // str(IERROR), __FILE__, __LINE__)
    ENDIF
#else
! determine processor grid, according to naive guess
    CALL FERMAT_RAZOR( NPROCS, NPROW, NPCOL )

    CALL MPI_barrier( GDES_MAT%COMM_INTRA%MPI_COMM, IERROR )
    CALL PROCMAP( GDES_MAT%COMM_INTRA, GDES_MAT%ICTXT, 2, NPROW, NPCOL)

! calculate local size of matrices

    CALL BLACS_GRIDINFO( GDES_MAT%ICTXT, NPROW, NPCOL, MYROW, MYCOL )

    ! set blocking factor compatible to scala.F
    ! this allows to call INIT_scala and use simple functions
    ! defined in scala.F
    !CALL QUERRY_SCALA_NB( NB_SCALA )
    ! mK: INIT_scala is actually never used, so we allow NB_SCALA to differ from
    !     the value in scala.F

    NP = NUMROC(NB_TOT,NB_SCALA, MYROW,0,NPROW)   ! get number of rows on proc
    NQ = NUMROC(NB_TOT,NB_SCALA, MYCOL,0,NPCOL)   ! get number of cols on proc

    CALL DESCINIT(GDES_MAT%DESC, NB_TOT,NB_TOT, NB_SCALA,NB_SCALA, 0,0, GDES_MAT%ICTXT, MAX( 1, NP ), IERROR ) ! setup descriptor
    IF (IERROR.NE.0) THEN
      CALL vtutor%bug("GREENS_MAT_DES_INIT: DESC, DESCINIT, INFO: " // str(IERROR) // " " // str(NB_TOT) // " " // str(NB_SCALA) // " " // str(NP),__FILE__, __LINE__)
    ENDIF
#endif
    GDES_MAT%MY_NROWS = NP
    GDES_MAT%MY_NCOLS = NQ
  END SUBROUTINE GREENS_MAT_DES_INIT

  
!*********************************************************************
!
! allocate an array for the Green function stored in the orbital
! basis
!
!*********************************************************************

  SUBROUTINE ALLOCATE_GREENS_MAT(GDES_MAT, CHAM_MAT,  NKPTS_IRZ, ISPIN, NOMEGA)
    USE ini
    USE tutor, ONLY: VTUTOR
    TYPE (greens_mat_des), POINTER :: GDES_MAT ! in: descriptor
    COMPLEX(q), POINTER :: CHAM_MAT(:,:,:,:)   ! out: matrix in orbital basis
    INTEGER             :: NKPTS_IRZ           ! in: number of k-points
    INTEGER             :: ISPIN               ! in: number of spin degrees of freedom
    INTEGER             :: NOMEGA              ! in: number of frequency points (or 1)
  ! local
    INTEGER :: NPROW, NPCOL, MYROW, MYCOL, NP, NQ
    INTEGER,EXTERNAL :: NUMROC
    INTEGER          :: ISTAT 

! calculate local size of matrices

    CALL BLACS_GRIDINFO( GDES_MAT%ICTXT, NPROW, NPCOL, MYROW, MYCOL )

    NP = NUMROC(GDES_MAT%NROWS, GDES_MAT%DESC(MB_), MYROW, 0, NPROW)
    NQ = NUMROC(GDES_MAT%NROWS, GDES_MAT%DESC(NB_), MYCOL, 0, NPCOL)

    ALLOCATE(CHAM_MAT(NP*NQ, NKPTS_IRZ, ISPIN, NOMEGA), STAT = ISTAT ) 
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "ALLOCATE_GREENS_MAT is not able to allocate "//&
         str( 16._q*NP*NQ*NKPTS_IRZ*ISPIN*NOMEGA/1024 ) //&
         " kB of data on MPI rank 0." )
    ENDIF
    CALL REGISTER_ALLOCATE(16._q*SIZE(CHAM_MAT,KIND=qi8), "Greenorb")
    
  END SUBROUTINE ALLOCATE_GREENS_MAT


  SUBROUTINE DEALLOCATE_GREENS_MAT(GDES_MAT, CHAM_MAT)
    USE ini
    TYPE (greens_mat_des), POINTER :: GDES_MAT ! in: descriptor
    COMPLEX(q), POINTER :: CHAM_MAT(:,:,:,:)   ! out: matrix in orbital basis

    IF (ASSOCIATED(CHAM_MAT)) THEN
       CALL DEREGISTER_ALLOCATE(16._q*SIZE(CHAM_MAT,KIND=qi8), "Greenorb")
       DEALLOCATE(CHAM_MAT)
       NULLIFY(CHAM_MAT)
    ELSE
       WRITE(*,*) 'internal warning in DEALLOCATE_GREENS_MAT: matrix is not allocated'
    ENDIF
    
  END SUBROUTINE DEALLOCATE_GREENS_MAT

!
! same routines as above
! in imaginary time the matrix is real valued (and symmetric) at the Gamma point
! hence once can use GDEF instead of COMPLEX(q)
!

  SUBROUTINE ALLOCATE_GREENS_MAT_GDEF(GDES_MAT, CHAM_MAT,  NKPTS_IRZ, ISPIN, NOMEGA)
    USE ini
    USE tutor, ONLY: VTUTOR
    TYPE (greens_mat_des), POINTER :: GDES_MAT ! in: descriptor
    GDEF, POINTER       :: CHAM_MAT(:,:,:,:)   ! out: matrix in orbital basis
    INTEGER             :: NKPTS_IRZ           ! in: number of k-points
    INTEGER             :: ISPIN               ! in: number of spin degrees of freedom
    INTEGER             :: NOMEGA              ! in: number of frequency points (or 1)
  ! local
    INTEGER :: NPROW, NPCOL, MYROW, MYCOL, NP, NQ
    INTEGER,EXTERNAL :: NUMROC
    INTEGER          :: ISTAT 
    
! calculate local size of matrices

    CALL BLACS_GRIDINFO( GDES_MAT%ICTXT, NPROW, NPCOL, MYROW, MYCOL )

    NP = NUMROC(GDES_MAT%NROWS, GDES_MAT%DESC(MB_), MYROW, 0, NPROW)
    NQ = NUMROC(GDES_MAT%NROWS, GDES_MAT%DESC(NB_), MYCOL, 0, NPCOL)

    ALLOCATE(CHAM_MAT(NP*NQ, NKPTS_IRZ, ISPIN, NOMEGA), STAT = ISTAT)
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "ALLOCATE_GREENS_MAT_GDEF is not able to allocate "//&
         str( wsgf*NP*NQ*NKPTS_IRZ*ISPIN*NOMEGA/1024 ) //&
         " kB of data on MPI rank 0." )
    ENDIF
    CHAM_MAT=0
    CALL REGISTER_ALLOCATE(wsgf*SIZE(CHAM_MAT,KIND=qi8), "Greenorb")
    
  END SUBROUTINE ALLOCATE_GREENS_MAT_GDEF


  SUBROUTINE DEALLOCATE_GREENS_MAT_GDEF(GDES_MAT, CHAM_MAT)
    USE ini
    TYPE (greens_mat_des), POINTER :: GDES_MAT ! in: descriptor
    GDEF, POINTER :: CHAM_MAT(:,:,:,:)   ! out: matrix in orbital basis

    IF (ASSOCIATED(CHAM_MAT)) THEN
       CALL DEREGISTER_ALLOCATE(wsgf*SIZE(CHAM_MAT,KIND=qi8), "Greenorb")
       DEALLOCATE(CHAM_MAT)
       NULLIFY(CHAM_MAT)
    ELSE
       WRITE(*,*) 'internal warning in DEALLOCATE_GREENS_MAT: matrix is not allocated'
    ENDIF
    
  END SUBROUTINE DEALLOCATE_GREENS_MAT_GDEF

!=======================================================================
!
! small routine to dump a distributed matrix
! the descriptor in DESCA must properly describe the matrix
! since RECON_SLICE calls check 
!
!=======================================================================
    
  SUBROUTINE DUMP_GREENS_HAM_GDEF( STRING, CHAM_DISTRI, NB_TOT, GDES_MAT, IU_,IUNIT)
    USE wave
    USE scala
    IMPLICIT NONE
    CHARACTER (LEN=*) :: STRING            ! string to dump
    GDEF              :: CHAM_DISTRI(*)    ! distributed matrix
    TYPE (greens_mat_des):: GDES_MAT       ! in: descriptor
    INTEGER           :: NB_TOT
    INTEGER           :: IU_                ! unit to write to (not dump for IU<0)
    INTEGER,OPTIONAL  :: IUNIT             ! unit to write to (not dump for IU<0)
  ! local
    INTEGER :: NDUMP=16
    GDEF,ALLOCATABLE   :: CHAM(:,:)
    INTEGER IU
    INTEGER N1, N2
    INTEGER COLUMN_HIGH, COLUMN_LOW
    LOGICAL :: LGAMMA_PRINT = .TRUE.
LGAMMA_PRINT = .FALSE. 

    IF ( PRESENT( IUNIT ) ) THEN 
       IU = IUNIT
       NDUMP =  NB_TOT
    ELSE
       IU = IU_
    ENDIF
    COLUMN_LOW=1

    COLUMN_HIGH=MIN(COLUMN_LOW+NDUMP-1,NB_TOT)

    ALLOCATE(CHAM(NB_TOT, COLUMN_HIGH-COLUMN_LOW+1))

    CALL RECON_SLICE(CHAM, NB_TOT, NB_TOT, CHAM_DISTRI,  GDES_MAT%DESC, COLUMN_LOW, COLUMN_HIGH)

    CALLMPI( M_sum_g(GDES_MAT%COMM_INTRA, CHAM(1,1), SIZE(CHAM)))


    IF (IU>=0) THEN
    WRITE(IU,*) STRING
    DO N1=1,NDUMP
    IF ( PRESENT( IUNIT ) ) THEN 
       DO N2 = 1, NDUMP
          WRITE(IU,10)N1+COLUMN_LOW-1, N2+COLUMN_LOW-1, CHAM(N1+COLUMN_LOW-1,N2)
       ENDDO
    ELSE
       WRITE(IU,1)N1+COLUMN_LOW-1,(REAL( CHAM(N1+COLUMN_LOW-1,N2) ,KIND=q) ,N2=1,NDUMP)
    ENDIF
    ENDDO
    WRITE(IU,*)
#ifndef gammareal
IF ( LGAMMA_PRINT ) THEN
    DO N1=1,NDUMP
    IF ( .NOT. PRESENT( IUNIT ) ) THEN 
       WRITE(IU,2)N1+COLUMN_LOW-1,(AIMAG( CHAM(N1+COLUMN_LOW-1,N2)),N2=1,NDUMP)
    ENDIF
    ENDDO
    WRITE(IU,*)
ENDIF
#endif
    ENDIF

    DEALLOCATE(CHAM)
1   FORMAT(1I3,3X,40F10.4)
2   FORMAT(1I3,3X,40F10.4)
10   FORMAT(2I4,2F20.10)

  END SUBROUTINE DUMP_GREENS_HAM_GDEF


  SUBROUTINE DUMP_GREENS_HAM( STRING, CHAM_DISTRI, NB_TOT, GDES_MAT, IU_,IUNIT)
    USE wave
    USE scala
    IMPLICIT NONE
    CHARACTER (LEN=*) :: STRING            ! string to dump
    COMPLEX(q)        :: CHAM_DISTRI(*)    ! distributed matrix
    TYPE (greens_mat_des):: GDES_MAT       ! in: descriptor
    INTEGER           :: NB_TOT
    INTEGER           :: IU_               ! unit to write to (not dump for IU<0)
    INTEGER,OPTIONAL  :: IUNIT             ! unit to write to (not dump for IU<0)
  ! local
    INTEGER           :: IU                ! unit to write to (not dump for IU<0)
    INTEGER, PARAMETER :: NDUMP=8
    COMPLEX(q),ALLOCATABLE  :: CHAM(:,:)
    INTEGER N1, N2
    INTEGER COLUMN_HIGH, COLUMN_LOW

    IF ( PRESENT( IUNIT ) ) THEN 
       IU = IUNIT
    ELSE
       IU = IU_
    ENDIF
    COLUMN_LOW=1

    COLUMN_HIGH=MIN(COLUMN_LOW+NDUMP-1,NB_TOT)

    ALLOCATE(CHAM(NB_TOT, COLUMN_HIGH-COLUMN_LOW+1))

    CALL RECON_SLICE_C(CHAM, NB_TOT, NB_TOT, CHAM_DISTRI,  GDES_MAT%DESC, COLUMN_LOW, COLUMN_HIGH)

    CALLMPI( M_sum_z(GDES_MAT%COMM_INTRA, CHAM(1,1), SIZE(CHAM)))

    IF (IU>=0) THEN
    WRITE(IU,*) STRING
    DO N1=1,NDUMP
       WRITE(IU,1)N1+COLUMN_LOW-1,(REAL( CHAM(N1,N2+COLUMN_LOW-1) ,KIND=q) ,N2=1,NDUMP)
    ENDDO
    WRITE(IU,*)
#ifndef gammareal
    DO N1=1,NDUMP
       WRITE(IU,2)N1+COLUMN_LOW-1,(AIMAG( CHAM(N1,N2+COLUMN_LOW-1)),N2=1,NDUMP)
    ENDDO
    WRITE(IU,*)
#endif
    ENDIF

    DEALLOCATE(CHAM)
1   FORMAT(1I2,3X,40F16.8)
2   FORMAT(1I2,3X,40F16.8)

  END SUBROUTINE DUMP_GREENS_HAM


#endif
END MODULE GG_base

!*********************************************************************** 
!
! helper F77 functions to subtract density-density contributions
! from response function
!
!*********************************************************************** 


   SUBROUTINE SUBTRACT_RESPONSE_RANK1( RESPONSEFUN, ND1, ND2, GCHG, WEIGHT,  NROWS, NCOLS, COMM_IN)
     USE prec
     USE mpimy
     IMPLICIT NONE
     INTEGER    :: ND1, ND2                    ! dimensions of response function array
     COMPLEX(q) :: RESPONSEFUN(ND1, ND2)       ! response function array
     COMPLEX(q) :: GCHG(ND1)
     REAL(q)    :: WEIGHT
     INTEGER                 :: NROWS          ! number of rows in CHI%RESPONSEFUN
     INTEGER                 :: NCOLS          ! number of columns in CHI%RESPONSEFUN
     TYPE (communic)         :: COMM_IN        ! communicator in group 
   ! local
     INTEGER  :: I, J
     INTEGER  :: NODE_ME = 1
#ifdef MPI
     NODE_ME =COMM_IN%NODE_ME
#endif
     
     IF (NROWS > ND1 .OR. NCOLS> ND2 ) THEN
        WRITE(*,*) 'internal errors in SUBTRACT_RESPONSE_RANK1: wrong dimensions'
        WRITE(*,*) NROWS, ND1, NCOLS, ND2
     ENDIF
!        WRITE(*,*) NROWS, ND1, NCOLS, ND2
     ! subtract state from response function
     ! conjugated version of what is used in chi_base.F
     DO J=1,NCOLS
        DO I=1,NROWS
           RESPONSEFUN(I,J)=RESPONSEFUN(I,J)-CONJG(GCHG(I))*GCHG((NODE_ME-1)*NCOLS+J )*WEIGHT
        ENDDO
     ENDDO
     
   END SUBROUTINE SUBTRACT_RESPONSE_RANK1

   SUBROUTINE SUBTRACT_RESPONSER_RANK1( RESPONSER, ND1, ND2, GCHG, WEIGHT,  NROWS, NCOLS, COMM_IN)
     USE prec
     USE mpimy
     INTEGER    :: ND1, ND2                    ! dimensions of response function array
     REAL(q)    :: RESPONSER(ND1, ND2)       ! response function array
     REAL(q)    :: GCHG(ND1)
     REAL(q)    :: WEIGHT
     INTEGER                 :: NROWS          ! number of rows in CHI%RESPONSER
     INTEGER                 :: NCOLS          ! number of columns in CHI%RESPONSER
     TYPE (communic)         :: COMM_IN        ! communicator in group 
   ! local
     INTEGER  :: I, J
     INTEGER  :: NODE_ME = 1
#ifdef MPI
     NODE_ME =COMM_IN%NODE_ME
#endif
     
     IF (NROWS > ND1 .OR. NCOLS> ND2 ) THEN
        WRITE(*,*) 'internal errors in SUBTRACT_RESPONSE_RANK1: wrong dimensions'
        WRITE(*,*) NROWS, ND1, NCOLS, ND2
     ENDIF
     ! subtract state from response function
     ! conjugated version of what is used in chi_base.F
     DO J=1,NCOLS
        DO I=1,NROWS
           RESPONSER(I,J)=RESPONSER(I,J)-GCHG(I)*GCHG((NODE_ME-1)*NCOLS+J )*WEIGHT
        ENDDO
     ENDDO
     
   END SUBROUTINE SUBTRACT_RESPONSER_RANK1

!************************ SUBROUTINE RC_RL_ADD  ***************************
!
!  This modul contains a helper routine to perform the operation
!     C = A * SCALE1 + B * SCALE2
!  in RECIPROCAL SPACE
!  C might be equal to A
!
!***********************************************************************

  SUBROUTINE RC_RL_ADD(A,SCALE1,B,SCALE2,C,GRID)
#ifdef _OPENACC
    USE mopenacc_struct_def
#endif
    USE mgrid
    IMPLICIT REAL(q) (A-H,O-Z)

    TYPE (grid_3d) GRID
    COMPLEX(q) C(GRID%RC%NP),A(GRID%RC%NP),B(GRID%RC%NP)
    NP=GRID%RC%NP
    !
    !   C=A*SCALE1
    !
    !  )  case SCALE1=0
    IF (SCALE1==0)  THEN
!$ACC PARALLEL LOOP PRESENT(C) __IF_ASYNC__
       DO 100 N=1,NP
          C(N)=0
100    ENDDO
       !  )  case SCALE1=1
    ELSE IF (SCALE1==1)  THEN
!$ACC PARALLEL LOOP PRESENT(C,A) __IF_ASYNC__
       DO 200 N=1,NP
          C(N)=A(N)
200    ENDDO
       !  )  else
    ELSE
!$ACC PARALLEL LOOP PRESENT(C,A) __IF_ASYNC__
       DO 300 N=1,NP
          C(N)=A(N)*SCALE1
300    ENDDO
    ENDIF

    !
    !   C=C+ B*SCALE2
    !
    IF (SCALE2==0) THEN
    ELSE IF (SCALE2==1) THEN
!$ACC PARALLEL LOOP PRESENT(C,B) __IF_ASYNC__
       DO 400 N=1,NP
          C(N)=C(N)+B(N)
400    ENDDO
    ELSE
!$ACC PARALLEL LOOP PRESENT(C,B) __IF_ASYNC__
       DO 500 N=1,NP
          C(N)=C(N)+B(N)*SCALE2
500    ENDDO
    ENDIF

    RETURN
  END SUBROUTINE RC_RL_ADD



