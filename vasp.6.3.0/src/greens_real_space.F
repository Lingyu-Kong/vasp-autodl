!#define debug
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
!
! similar routines for the self-energy are included  
!
! the following routines are used in chi_super and chi_GG
! which differ mostly in the way k-points are treated in the 
! computation of the polarizability CHI and the self-energy SIGMA
! 
! chi_GG.F uses a simple convolution over k-points in the primitive 
! cell (relatively small memory footprint).
! 
! chi_super.F uses a more elaborate supercell approach that achieves
! linear scaling in the number of k-points, but has a larger memory 
! footprint. 
!
!*********************************************************************

MODULE greens_real_space
  USE chi_glb
  USE chi_base
  USE wave_high
  USE mlr_optic
  USE GG_base
  USE greens_orbital
  USE rpa_force
  USE minimax_struct
  USE ini
  IMPLICIT NONE
  PRIVATE

#ifdef scaLAPACK
  PUBLIC :: DUMP_STORAGE_REQUIREMENTS
  PUBLIC :: DETERMINE_CHI_Q_OMEGA
  PUBLIC :: SCREENED_POTENTIAL
  PUBLIC :: RECALCULATE_W_AND_PHI_FUNCTIONAL
  PUBLIC :: NULLIFY_G_RECIPROCAL
  PUBLIC :: CALCULATE_G_POSSIBLY_SC
  PUBLIC :: ROTATE_RES
  PUBLIC :: EDDIAG_RSE
  PUBLIC :: CALCULATE_SIGMA_ORBITAL
  PUBLIC :: DEALLOCATE_G
  PUBLIC :: EDDIAG_SIMPLE
  PUBLIC :: CONTRACT_SIGMA_G_TAU
  PUBLIC :: ROTATE_ORBITALS
  PUBLIC :: ROTATE_KS_ORBITALS
  PUBLIC :: NO_SYMMETRY
  PUBLIC :: ON_SYMMETRY
  ! exclusive chi_super routines 
  PUBLIC :: CALCULATE_SIGMA_SUPER
  ! exclusive chi_GG routines
  PUBLIC :: CALCULATE_SIGMA_TAU
  PUBLIC :: CONTRACT_SIGMA_G_TAU_DER
  PUBLIC :: ALLOCATE_G_RECIPROCAL
  PUBLIC :: UNCOMPRESS_G_RECIPROCAL
 CONTAINS

!***********************************************************************
!
! allocate the Greens function
! there is one version to allocate the Greens function
! in real space
! and one version to allocate in reciprocal space
! the reciprocal space version saves a lot of space
!
!***********************************************************************

  SUBROUTINE ALLOCATE_G_RECIPROCAL(GDES, G, NKPT)
    USE ini
    TYPE (greensfdes):: GDES
    TYPE (greensf) :: G(:)
    INTEGER  :: NKPT
  ! local
    INTEGER  :: IK
    INTEGER  :: ISTAT 

    CALL  NULLIFY_G_RECIPROCAL(GDES, G, NKPT)

    DO IK=1,NKPT
       ALLOCATE( G(IK)%GG    (GDES%NRPLWV_ROW_DATA_POINTS, GDES%NRPLWV_COL_DATA_POINTS), & 
          STAT = ISTAT ) 
       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "ALLOCATE_G_RECIPROCAL (G("//str(IK)//"%GG) is not able to allocate "//&
            str( ndata * 8._q*GDES%NRPLWV_ROW_DATA_POINTS*GDES%NRPLWV_COL_DATA_POINTS/1024) //&
            " kB of data on MPI rank 0." )
       ENDIF

       ALLOCATE( G(IK)%G_PROJ(GDES%NRPLWV_ROW_DATA_POINTS, GDES%NPRO_COL  ), & 
          STAT = ISTAT ) 
       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "ALLOCATE_G_RECIPROCAL (G("//str(IK)//"%G_PROJ) is not able to allocate "//&
            str( ndata * 8._q*GDES%NRPLWV_ROW_DATA_POINTS*GDES%NPRO_COL/1024) //&
            " kB of data on MPI rank 0." )
       ENDIF

       ALLOCATE( G(IK)%PROJ_PROJ(GDES%NPRO_ROW,GDES%NPRO_COL), STAT = ISTAT ) 
       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "ALLOCATE_G_RECIPROCAL (G("//str(IK)//"%G_PROJ_PROJ) is not able to allocate "//&
            str( ndata * 8._q*GDES%NPRO_ROW*GDES%NPRO_COL/1024) //&
            " kB of data on MPI rank 0." )
       ENDIF

       CALL REGISTER_ALLOCATE(wsgf*SIZE(G(IK)%GG,KIND=qi8), "Greensfun")
       CALL REGISTER_ALLOCATE(wsgf*SIZE(G(IK)%G_PROJ,KIND=qi8), "Greensfun")
       CALL REGISTER_ALLOCATE(wsgf*SIZE(G(IK)%PROJ_PROJ,KIND=qi8), "Greensfun")
       G(IK)%GG=0
       G(IK)%G_PROJ=0
       G(IK)%PROJ_PROJ=0
    ENDDO
  END SUBROUTINE ALLOCATE_G_RECIPROCAL

  SUBROUTINE NULLIFY_G_RECIPROCAL(GDES, G, NKPT)
    USE ini
    TYPE (greensfdes):: GDES
    TYPE (greensf) :: G(:)
    INTEGER  :: NKPT
  ! local
    INTEGER  ::  IK
   
    DO IK=1 , NKPT
       NULLIFY( G(IK)%RR         ,&          
                G(IK)%R_PROJ     ,&
                G(IK)%PROJ_R     ,&
                G(IK)%PROJ_PROJ  ,&
                G(IK)%GG         ,&
                G(IK)%G_PROJ     ,&
                G(IK)%G_R        ,&
                G(IK)%R_G        ,&
                G(IK)%PROJ_G )   
    ENDDO
  END SUBROUTINE NULLIFY_G_RECIPROCAL

!***********************************************************************
!
! This routine calculates all Green's functions GO and GU
! in the plane wave basis 
!
!***********************************************************************

  SUBROUTINE CALCULATE_G_POSSIBLY_SC( W, WDES, WHF, WCORR, WCORRHF, &
    GDES_TAU, GDES_MAT, NKPTS_IRZ, ISP,  T, ONE_ELECTRON_GREENS, NTAU_ROOT, &
    GU, GO, GU_MAT, GO_MAT, IU6_MEM, IO, LCORRELATED)
    USE minimax_struct, ONLY: loop_des
    USE base, ONLY: in_struct
    USE greens_orbital
    TYPE (wavespin)   :: W
    TYPE (wavespin)   :: WCORR
    TYPE (wavedes)    :: WDES
    TYPE (wavespin)   :: WHF
    TYPE (wavespin)   :: WCORRHF
    TYPE (greensfdes) :: GDES_TAU            ! Green function in time  domain
    TYPE (greens_mat_des), POINTER :: GDES_MAT
    INTEGER               :: NKPTS_IRZ
    INTEGER           :: ISP                 ! spin channel
    TYPE (loop_des)   :: T                   ! tau-loop handle 
    LOGICAL           :: ONE_ELECTRON_GREENS ! whether one electron Green's function is used
    TYPE (greensf), POINTER   :: GU(:)         ! Green function for tau>0 
    TYPE (greensf), POINTER   :: GO(:)         ! Green function for tau<0
    GDEF, POINTER, CONTIGUOUS :: GU_MAT(:,:,:,:)! matrix <i| G_u(tau) | a> stored distributed
    GDEF, POINTER, CONTIGUOUS :: GO_MAT(:,:,:,:)! matrix <i| G_u(tau) | a> stored distributed
    INTEGER                   ::  NTAU_ROOT                ! tau point 
    TYPE (in_struct)          :: IO            ! input output handle 
    INTEGER                   :: IU6_MEM       ! dump of timings or memory
    LOGICAL                   :: LCORRELATED 
  ! local 
    INTEGER IU6_TMP, NK1, NQ_IRZ
    LOGICAL           :: LCOMPUTE_CORRELATED
   
    PROFILING_START( 'calculate_g_possibly_sc' )

    ! in case correlated Green's function is calculated
    LCOMPUTE_CORRELATED = LCORRELATED 

    ALLOCATE(GU(W%WDES%NKPTS), GO(W%WDES%NKPTS ))
    ! nullify
    CALL NULLIFY_G_RECIPROCAL( GDES_TAU, GU, W%WDES%NKPTS)
    CALL NULLIFY_G_RECIPROCAL( GDES_TAU, GO, W%WDES%NKPTS)
    ! and allocate
    CALL ALLOCATE_G_RECIPROCAL( GDES_TAU, GU, W%WDES%NKPTS)
    CALL ALLOCATE_G_RECIPROCAL( GDES_TAU, GO, W%WDES%NKPTS)

    IU6_TMP=IU6_MEM

    CALL START_TIMING("G")
    !compute occupied and unoccupied Green function in tau domain
    IF (ONE_ELECTRON_GREENS) THEN
       ! use non-interacting Green function in the first step
       DO NK1=1,W%WDES%NKPTS
          IF (MOD(NK1-1,WDES%COMM_KINTER%NCPU)/=WDES%COMM_KINTER%NODE_ME-1) CYCLE
          ! computes correlated greens function of occupied states
          IF ( LCOMPUTE_CORRELATED ) THEN
             ! with original wave function
             IF ( ALLOCATED( LCRPA_BAND ) .OR. LWEIGHTED ) THEN
                CALL CALCULATE_G_RECIPROCAL( WHF, GO(NK1), GDES_TAU, T, IO,&
                    IU6_TMP, NK1, NBANDSGWLOW, LOCCUPIED=.TRUE., ISP=ISP,&
                    LCORRELATED = .TRUE.)
             ! with manipulated wave function
             ELSE
                CALL CALCULATE_G_RECIPROCAL( WCORRHF, GO(NK1), GDES_TAU, T, IO, &
                    IU6_TMP, NK1, NBANDSGWLOW, LOCCUPIED=.TRUE., ISP=ISP,    &
                    LCORRELATED = .TRUE.)
             ENDIF
          ELSE
             CALL CALCULATE_G_RECIPROCAL( WHF, GO(NK1), GDES_TAU, T, IO,&
                 IU6_TMP, NK1, NBANDSGWLOW, LOCCUPIED=.TRUE., ISP=ISP,&
                 LCORRELATED = .FALSE.)
          ENDIF

#ifdef RPAgamma
          IF ( ICHIREAL == 2 ) CALL COMPRESS_G_RECIPROCAL( GO(NK1), GDES_TAU, NK1)
#endif
          CALL STOP_TIMING("G",IU6_TMP,"*gocc")
          ! computes correlated greens function of unoccupied states
          IF ( LCOMPUTE_CORRELATED ) THEN
             ! with original wave function
             IF ( ALLOCATED( LCRPA_BAND ) .OR. LWEIGHTED ) THEN
                CALL CALCULATE_G_RECIPROCAL( WHF, GU(NK1), GDES_TAU, T, IO,&
                     IU6_TMP, NK1, NBANDSGWLOW, LOCCUPIED=LDMP1,ISP=ISP,&
                     LCORRELATED = .TRUE.)
             ! with manipulated wave function
             ELSE
                CALL CALCULATE_G_RECIPROCAL( WCORRHF, GU(NK1), GDES_TAU, T, IO,&
                     IU6_TMP, NK1, NBANDSGWLOW, LOCCUPIED=LDMP1,ISP=ISP,&
                     LCORRELATED = .TRUE.)
             ENDIF
          ELSE
             CALL CALCULATE_G_RECIPROCAL( WHF, GU(NK1), GDES_TAU, T, IO,&
                  IU6_TMP, NK1, NBANDSGWLOW, LOCCUPIED=LDMP1,ISP=ISP,&
                  LCORRELATED = .FALSE.)
          ENDIF
#ifdef RPAgamma
          IF ( ICHIREAL == 2 ) CALL COMPRESS_G_RECIPROCAL( GU(NK1), GDES_TAU, NK1)
#endif
          CALL STOP_TIMING("G",IU6_TMP,"*gunocc")
          IU6_TMP=-1
          ! dump progress
!          IPROGRESS = IPROGRESS + 1 
!          CALL DUMP_PROGRESS( IPROGRESS, IGOAL, 'G', IO )
       ENDDO
    ELSE
       IF ( LCRPA ) THEN
          CALL vtutor%error("Sorry self-consistent Greens functions for CRPA not implemented")
       ENDIF

       ! get the Green from GU_MAT and GO_MAT
       ! transform them from Bloch basis to plane-waves and projectors
       IF (LG0W0) THEN
          CALL vtutor%error("internal error in VASP: G0W0 method chosen but NELM larger than 1 \n for &
             &self-consistent GW, set ALGO=GW0R in INCAR \n will stop now")
       ENDIF
       DO NK1=1,W%WDES%NKPTS
          IF (MOD(NK1-1,WDES%COMM_KINTER%NCPU)/=WDES%COMM_KINTER%NODE_ME-1) CYCLE
          
          ! fix a) "GO_MAT_q"="GO_MAT_k", since these are expansion in the
          ! orbital basis
          IF (NK1>NKPTS_IRZ) THEN
             ! if we are outside IrBZ, use index to IrBZ for collecting
             ! GU_MAT
             NQ_IRZ=KPOINTS_FULL_ORIG%NEQUIV(KPOINT_IN_FULL_GRID(KPOINTS_FULL%VKPT(:,NK1),KPOINTS_FULL_ORIG))
          ELSE
             ! can use what we have now
             NQ_IRZ=NK1
          ENDIF
          
          ! in self-consistent version we need to transform 
          ! G(a,b',tau) to G(g,g',tau)
          ! transform from Bloch basis to plane-waves and projectors for
          ! unoccupied states
          ! the GU_MAT stores Greens function for all local TAU points, need
          ! to pass down current tau point
          ! T%NPOINTSC or NTAU_ROOT, not entirely sure
          CALL TRANSFORM_G_ORBIT_PW(W, GU(NK1), GU_MAT, GDES_MAT, GDES_TAU,IO%IU6, NK1, NQ_IRZ, ISP, &
               T%COMM_IN_GROUP, T%COMM_BETWEEN_GROUPS, T%LDO_POINT_LOCAL, NTAU_ROOT)
          CALL STOP_TIMING("G",IU6_TMP,"*gocc")
#ifdef RPAgamma
          IF ( ICHIREAL == 2 ) CALL COMPRESS_G_RECIPROCAL( GU(NK1), GDES_TAU, NK1)
#endif
          
          ! transform from Bloch basis to plane-waves and projectors for
          ! occupied states
          CALL TRANSFORM_G_ORBIT_PW(W, GO(NK1), GO_MAT, GDES_MAT, GDES_TAU,IO%IU6, NK1, NQ_IRZ, ISP, &
               T%COMM_IN_GROUP, T%COMM_BETWEEN_GROUPS, T%LDO_POINT_LOCAL, NTAU_ROOT)
#ifdef RPAgamma
          IF ( ICHIREAL == 2 ) CALL COMPRESS_G_RECIPROCAL( GO(NK1), GDES_TAU, NK1)
#endif
          CALL STOP_TIMING("G",IU6_TMP,"*gunocc")
          IU6_TMP=-1
          ! dump progress
!          CALL DUMP_PROGRESS( NK1, WDES%NKPTS, 'G_dyson', IO )
       ENDDO
    ENDIF !selection of NELM
    
    ! sync if KPAR is used (gK: commented in on 2017.07.05)
    CALL SYNC_G_RECIPROCAL(WDES, GO)
    CALL SYNC_G_RECIPROCAL(WDES, GU)
    
    PROFILING_STOP( 'calculate_g_possibly_sc' )
    CONTAINS 

    !***********************************************************************
    !
    ! allocate the Greens function
    ! there is one version to allocate the Greens function
    ! in real space
    ! and one version to allocate in reciprocal space
    ! the reciprocal space version saves a lot of space
    !
    !***********************************************************************
  
    SUBROUTINE SYNC_G_RECIPROCAL(WDES, G)
       USE ini
       TYPE (wavedes):: WDES
       TYPE (greensf) :: G(:)
     ! local
       INTEGER  ::  IK
   
       IF (WDES%COMM_KINTER%NCPU==1) RETURN
       DO IK=1,SIZE(G,1)
               CALLMPI( M_sum_z(WDES%COMM_KINTER, G(IK)%GG, SIZE(G(IK)%GG)) )
               CALLMPI( M_sum_z(WDES%COMM_KINTER, G(IK)%G_PROJ, SIZE(G(IK)%G_PROJ)) )
               CALLMPI( M_sum_z(WDES%COMM_KINTER, G(IK)%PROJ_PROJ, SIZE(G(IK)%PROJ_PROJ)) )
       ENDDO
    END SUBROUTINE SYNC_G_RECIPROCAL

    END SUBROUTINE CALCULATE_G_POSSIBLY_SC

!***********************************************************************
!
! calculate the Greens function in reciprocal space G(G,G',tau)
! in principle this the simpler, more efficient 
! and more space conserving version
! but it requires more FFT's at a later point
!
! finite-T: the appropriate finite temperature Green's function is
!  1/(iw - eps) <-> - exp(eps (beta-tau))/(exp(eps beta)+1)
!                 = - exp(-eps tau) 1/(1+ exp(-eps beta))
! this is well defined and finite in the intervall [0, beta]
!***********************************************************************

  SUBROUTINE CALCULATE_G_RECIPROCAL( W, G, GDES, T, IO, &
    IU6_, NK1, NBANDSGWLOW, LOCCUPIED, ISP, LCORRELATED, LTEST)
    USE dfast
    USE ini
    USE minimax_struct, ONLY: loop_des
    TYPE (wavespin)   :: W
    TYPE (greensf)    :: G
    TYPE (greensfdes) :: GDES  
    TYPE (loop_des)   :: T 
    TYPE (in_struct)  :: IO
    INTEGER :: IU6_            ! unit for IO (usually OUTCAR)
    INTEGER :: NK1             ! considered k-point
    INTEGER :: NBANDSGWLOW     ! lowest band to be included in Green's function
    INTEGER :: ISP             ! considered spin
    LOGICAL :: LOCCUPIED       ! Green's function for occupied or unoccupied orbitals
    LOGICAL :: LCORRELATED
    LOGICAL, OPTIONAL :: LTEST 
  ! local
    REAL(q) :: TAU
    REAL(q) :: TAUMIN          ! smallest time (tau) value of all groups
    LOGICAL :: LDO_TAU_LOCAL   ! may be GEMM calls can be skipped

    INTEGER :: NB_FIRST, NB_LAST
    GDEF, ALLOCATABLE       :: CPROJ(:,:)
    COMPLEX(q), ALLOCATABLE :: CG(:,:)
    COMPLEX(q), ALLOCATABLE :: CG_COLUMN   (:,:)
    GDEF, ALLOCATABLE       :: CPROJ_COLUMN(:,:)
    TYPE (wavefun1):: W1
    TYPE (wavedes1) :: WDES1
    INTEGER :: NB1, NB2, NSTRIP, NSTRIP_GLOBAL, N, NN, NGLB, NB1_INTO_TOT
    INTEGER :: NTAU, IERROR
    REAL(q), ALLOCATABLE :: WEIGHT(:)
    REAL(q) :: SWEIGHT,EXPO
    REAL(q) :: EFERMI          ! Fermi-energy, or chemical potential for T>0
  
    INTEGER,SAVE :: IU6=-1            ! unit for IO (usually OUTCAR)

    PROFILING_START('calculate_g_reciprocal')

    ! verbosity level changed with NWRITE in INCAR       
    IF ( IO%NWRITE> 2  ) THEN
       IF( IU6 < 0 )  IU6=IU6_
    ELSE  
       IU6=IU6_
    ENDIF 
     
    TAU = T%POINT_CURRENT
    LDO_TAU_LOCAL = T%LDO_POINT_LOCAL
    EFERMI = W%EFERMI(ISP)

    ! In case exact exchange is calculated with 
    ! LDMP1 = T; NOMEGA = 1 ; LFINITE_TEMPERATURE  = T
    ! one needs to set up G_occ( 0-) G^uno (beta - 0-)
    IF ( LFINITE_TEMPERATURE .AND. LDMP1 .AND. NOMEGA == 1 ) THEN 
       IF( LOCCUPIED ) THEN
          TAU = 0
       ELSE
          TAU = T%BETA 
       ENDIF
       T%POINTMIN_CURRENT = 0
    ENDIF 
    !having a loop over local time points in each group separately
    !has the effect that TAUMIN differs in each group.
    !consequently the number of unoccupied bands sent and received 
    !may differ from group to group, which can cause a deadlock
    !solution: take here the smallest tau value from group 1 for all nodes
    TAUMIN =T%POINTMIN_CURRENT

    IF (LOCCUPIED) THEN
       ! seek last unoccupied band to be included in Green's function at negative time
       DO NB1=W%WDES%NB_TOT, 1, -1
          ! seek last "occupied" band
          ! for T>0 those states propagate negatively in time
          IF ( LFINITE_TEMPERATURE ) THEN
             ! these are essentially fully and partially occupied states
             ! find last partially occupied state included in G_<
             IF ( GF_WEIGHT( .TRUE., TAUMIN, T%BETA, EFERMI, &
                     REAL(W%CELTOT(NB1,NK1,ISP),q), & 
                     REAL(W%FERTOT(NB1,NK1,ISP),q))  > EXP(-GPOS_THRESHHOLD) ) EXIT
          ELSE
             IF ( REAL(W%CELTOT(NB1,NK1,ISP),q)-(EFERMI+EFERMIADD)<0.0_q) EXIT
          ENDIF
       ENDDO 
       ! correlated band could be even lower
       IF ( LCORRELATED ) NB1 = MIN( NB1, NCRPAMAX )  
       ! round to next larger value modulo W%WDES%NB_PAR
       NB1=((NB1+W%WDES%NB_PAR-1)/W%WDES%NB_PAR)*W%WDES%NB_PAR
       ! convert global index into local index
       ! NB_LAST =(NB1-1)/W%WDES%NB_PAR+1  ! gK 06.03.2015  identical and used in other places
       NB_LAST =NB1/W%WDES%NB_PAR
       NB_FIRST=1

       IF (NBANDSGWLOW>0) THEN
          ! round to next smaller value
          NB1=((NBANDSGWLOW-1)/W%WDES%NB_PAR)*W%WDES%NB_PAR+1
          ! convert to local index
          NB_FIRST=MIN((NB1-1)/W%WDES%NB_PAR+1, W%WDES%NBANDS)
       ENDIF
       ! correlated band could be higher
       IF ( LCORRELATED ) THEN
          ! round to next smaller value
          NB1=((NCRPAMIN-1)/W%WDES%NB_PAR)*W%WDES%NB_PAR+1
          ! convert to local index
          NB_FIRST=MIN((NB1-1)/W%WDES%NB_PAR+1, W%WDES%NBANDS)
       ENDIF
    ELSE
       ! seek first band to be included in positive time Green's function
       DO NB1=1,W%WDES%NB_TOT
          ! for T>0 those states propagate positively in time
          IF ( LFINITE_TEMPERATURE ) THEN
             EXIT 
!orig
!             IF ( POSITIVE_TIME_ORBITAL( W%FERTOT(NB1,NK1,ISP),  &
!                   REAL(W%CELTOT(NB1,NK1,ISP)-EFERMI,q)*T%BETA ) ) EXIT
             ! first band that has essentially no contribution determines NB1
!gK
!             IF (REAL(W%CELTOT(NB1,NK1,ISP),q)-EFERMI < 0 .AND. & 
!                 GF_WEIGHT( .FALSE., TAU, T%BETA, EFERMI, &
!                      REAL(W%CELTOT(NB1,NK1,ISP),q), & 
!                      REAL(W%FERTOT(NB1,NK1,ISP),q))  < EXP(-GPOS_THRESHHOLD) ) EXIT
          ELSE
             ! CB(min)
             IF ( REAL(W%CELTOT(NB1,NK1,ISP),q)-(EFERMI-EFERMIADD)>0.0_q) EXIT
          ENDIF
       ENDDO
       ! correlated GF does not contain states below NCRPAMIN
       IF ( LCORRELATED ) NB1 = MAX( NB1, NCRPAMIN )  
       ! round to next smaller value
       NB1=((NB1-1)/W%WDES%NB_PAR)*W%WDES%NB_PAR+1
       ! convert to local index
       NB_FIRST=MIN((NB1-1)/W%WDES%NB_PAR+1, W%WDES%NBANDS)

       ! current time point might yield a negilgible contribution for unoccupied GF
       N = W%WDES%NB_TOT
       IF ( LCORRELATED ) N = NCRPAMAX
       ! seek last band that yields a non negligable contribution
       ! this contribution is easy to dermine for T=0, 
       ! namely exponent must be smaller than GPOS_THRESHHOLD 
       IF ( LFINITE_TEMPERATURE ) THEN
          ! for T>0, this is the last state that propagates positively in time
          DO NB1 = N, 1, -1
!orig
!             IF ( POSITIVE_TIME_ORBITAL( W%FERTOT(NB1,NK1,ISP),  &
!                   REAL(W%CELTOT(NB1,NK1,ISP)-EFERMI,q)*T%BETA ) ) EXIT
!new
             IF ( GF_WEIGHT( .FALSE., TAUMIN, T%BETA, EFERMI, &
                      REAL(W%CELTOT(NB1,NK1,ISP),q), & 
                      REAL(W%FERTOT(NB1,NK1,ISP),q))  > EXP(-GPOS_THRESHHOLD) ) EXIT
          ENDDO
       ELSE
          ! T=0 code 
          DO NB1 = N, 1, -1
             IF (ABS(TAUMIN*(REAL(W%CELTOT(NB1,NK1,ISP),q)-EFERMI)) < GPOS_THRESHHOLD) EXIT
          ENDDO
       ENDIF

       ! round to next larger value modulo W%WDES%NB_PAR
       NB1=((NB1+W%WDES%NB_PAR-1)/W%WDES%NB_PAR)*W%WDES%NB_PAR
       ! convert global index into local index
       ! NB_LAST =(NB1-1)/W%WDES%NB_PAR+1  ! gK 06.03.2015  identical and used in other places
       NB_LAST =NB1/W%WDES%NB_PAR
    ENDIF

    IF (IU6>=0) THEN
       IF (LOCCUPIED) THEN
          IF ( LCORRELATED ) THEN
             WRITE(IU6,'(/A,/2(A,I5,2X,A,I5/))')' Correlated GF: Bands included as occupied valence bands',&
               ' VB(min)=',(NB_FIRST-1)*W%WDES%NB_PAR+1,'VB(max)=',(NB_LAST)*W%WDES%NB_PAR
          ELSE
             WRITE(IU6,'(/A,/2(A,I5,2X,A,I5/))')' Bands included as occupied valence bands',&
               ' VB(min)=',(NB_FIRST-1)*W%WDES%NB_PAR+1,'VB(max)=',(NB_LAST)*W%WDES%NB_PAR
          ENDIF
       ELSE
          IF ( LCORRELATED ) THEN
             WRITE(IU6,'(/A,/2(A,I5,2X,A,I5/))')' Correlated GF: Bands included as un-occupied conduction bands', & 
               ' CB(min)=',(NB_FIRST-1)*W%WDES%NB_PAR+1,'CB(max)=',(NB_LAST)*W%WDES%NB_PAR
          ELSE
             WRITE(IU6,'(/A,/2(A,I5,2X,A,I5/))')' Bands included as un-occupied conduction bands', & 
               ' CB(min)=',(NB_FIRST-1)*W%WDES%NB_PAR+1,'CB(max)=',(NB_LAST)*W%WDES%NB_PAR
          ENDIF
       ENDIF
    ENDIF

    NSTRIP_GLOBAL=NSTRIP_STANDARD_GLOBAL
#ifdef MPI
    NSTRIP=(NSTRIP_GLOBAL+W%WDES%COMM_INTER%NCPU-1)/W%WDES%COMM_INTER%NCPU
    NSTRIP_GLOBAL=NSTRIP*W%WDES%COMM_INTER%NCPU
#else
    NSTRIP=NSTRIP_GLOBAL
#endif
    ALLOCATE(CG(GDES%NRPLWV_ROW,NSTRIP_GLOBAL),CPROJ(GDES%NPRO_ROW,NSTRIP_GLOBAL), &
             CG_COLUMN(GDES%NRPLWV_COL,NSTRIP_GLOBAL),CPROJ_COLUMN(GDES%NPRO_COL,NSTRIP_GLOBAL), &
             WEIGHT(NSTRIP_GLOBAL))

    CPROJ=0.0_q
    CG=0.0_q
    CPROJ_COLUMN=0.0_q
    CG_COLUMN=0.0_q


    CALL SETWDES(W%WDES,WDES1,NK1)
    CALL NEWWAV(W1, WDES1, .TRUE.)

    G%GG=0
    G%G_PROJ=0
    G%PROJ_PROJ=0
! use weight of first (Gamma) point
    SWEIGHT=SQRT(W%WDES%WTKPT(1)*W%WDES%RSPIN)
!=======================================================================
! loop over all bands 
!=======================================================================
    nblock: DO NB1=NB_FIRST, NB_LAST, NSTRIP
       NB2=MIN(NB_LAST, NB1+NSTRIP-1)
          ! gather all bands with local index NB1 up to and including NB2
       CALL W1_GATHER_ARRAY_RECIPROCAL(W, NB1, NB2, ISP, W1, CG, CPROJ)

       ! total number of bands  W1_GATHER_ARRAY has collected
       NGLB=(NB2-NB1+1)*W%WDES%NB_PAR

       ! copy relevant part over to CG_COLUMN
       CG_COLUMN(:,1:NGLB)=CG(GDES%NRPLWV_POS:GDES%NRPLWV_POS+GDES%NRPLWV_COL-1,1:NGLB)

       IF (GDES%NPRO_COL/=0) THEN
          CALL GDIS_PROJ(GDES, CPROJ(:,1:NGLB), CPROJ_COLUMN(:,1:NGLB))
       ENDIF
       
       ! TODO it would be nicer to copy the relevant part only once
       ! and multiply with the product new WEIGHT / last WEIGHT
       ! weight is determined on CPU 
       WEIGHT=0
       ngathered:DO N=1,NGLB
          NB1_INTO_TOT=(NB1-1)*W%WDES%NB_PAR+N
          ! if contribution is outside regime of included bands, skip contribution
          IF (NB1_INTO_TOT<NBANDSGWLOW) THEN
             WEIGHT(N)=0
          ! calculate contribution
          ELSE
             IF ( LFINITE_TEMPERATURE ) THEN
                WEIGHT(N) = SWEIGHT* GF_WEIGHT( LOCCUPIED, TAU, T%BETA, EFERMI, &
                   REAL(W%CELTOT(NB1_INTO_TOT,NK1,ISP),q), REAL(W%FERTOT(NB1_INTO_TOT,NK1,ISP),q) )
             ELSE   
                ! T=0 case 
                IF (ABS(TAU*(REAL(W%CELTOT(NB1_INTO_TOT,NK1,ISP),q)-EFERMI))>GPOS_THRESHHOLD) THEN
                   WEIGHT(N)=0
                ELSE 
                   IF (LOCCUPIED .AND. REAL(W%CELTOT(NB1_INTO_TOT,NK1,ISP),q)-(EFERMI+EFERMIADD)<0.0_q) THEN
                      EXPO = TAU*(REAL(W%CELTOT(NB1_INTO_TOT,NK1,ISP),q)-EFERMI) 
                      IF ( EXPO < 0.0_q ) THEN
                         ! TODO: the occupied and unoccupied Greens function has the sqrt of the k-point weight and the
                         ! the spin multiplicity multiplied in
                         ! this is not very elegant, but simple to do, ideally those weights
                         ! should be taken care of upon contraction e.g. X = G G 
                         WEIGHT(N)=EXP(  TAU*(REAL(W%CELTOT(NB1_INTO_TOT,NK1,ISP),q)-EFERMI))*SWEIGHT
                      ELSE
                         WEIGHT(N)=SWEIGHT
                      ENDIF
                      IF (EFERMIADD>0) WEIGHT(N)=WEIGHT(N)*W%FERTOT(NB1_INTO_TOT,NK1,ISP)
                   ELSE IF (.NOT. LOCCUPIED .AND.  REAL(W%CELTOT(NB1_INTO_TOT,NK1,ISP),q)-(EFERMI-EFERMIADD) >0.0_q) THEN
                      EXPO = TAU*(REAL(W%CELTOT(NB1_INTO_TOT,NK1,ISP),q)-EFERMI)
                      IF (EXPO>0.0_q) THEN
                         WEIGHT(N)=EXP(-(TAU*(REAL(W%CELTOT(NB1_INTO_TOT,NK1,ISP),q)-EFERMI)))*SWEIGHT
                      ELSE
                         WEIGHT(N)=SWEIGHT
                      ENDIF
                      IF (EFERMIADD>0) WEIGHT(N)=WEIGHT(N)*(1-W%FERTOT(NB1_INTO_TOT,NK1,ISP))
                   ENDIF
                ENDIF
             ENDIF
          ! TODO: the occupied and unoccupied Greens function has the sqrt of the k-point weight and the
          ! the spin multiplicity multiplied in
          ! this is not very elegant, but simple to do, ideally those weights
          ! should be taken care of upon contraction e.g. X = G G 
          ENDIF
#ifdef debug
          ! this line allows to include only partially occupied bands, &
          ! which are taken out in SUBTRACT_SELF_CORRELATION. 
          ! It can be very usefull for tests
          IF ( ABS(W%FERTOT(NB1_INTO_TOT,NK1,ISP))< XI_EMPTY_THRESHHOLD .OR. &
               ABS(W%FERTOT(NB1_INTO_TOT,NK1,ISP)-1)< XI_EMPTY_THRESHHOLD) WEIGHT(N)=0
#endif 

          ! in case correlated Green's function is requested
          IF ( LCORRELATED ) THEN
             !weighted CRPA computes the correlated Green's function as follows:
             IF ( ALLOCATED( CRPA_WEIGHTS) ) THEN
                WEIGHT(N)=WEIGHT(N)*CRPA_WEIGHTS(NB1_INTO_TOT,NK1,ISP)  
             ENDIF
             !correlated greens function if full bands are removed
             IF ( ALLOCATED( LCRPA_BAND ) ) THEN
                IF ( .NOT. LCRPA_BAND(NB1_INTO_TOT) ) WEIGHT(N)=0
             ENDIF
             !disentangled CRPA computes the correlated Green's function as follows:
             IF ( ALLOCATED( LCRPA_STATE)) THEN
                IF ( .NOT. LCRPA_STATE(ISP,NK1,NB1_INTO_TOT) ) WEIGHT(N)=0
             ENDIF
          ENDIF 
       ENDDO ngathered

       ! apply exponential phase factor to each entry
       DO N=1,NGLB
          CG_COLUMN(:,N)   =CG_COLUMN(:,N)*WEIGHT(N)
       ENDDO
       IF (GDES%NPRO_COL/=0) THEN
          DO N=1,NGLB
             CPROJ_COLUMN(:,N)=CPROJ_COLUMN(:,N)*WEIGHT(N)
          ENDDO
       ENDIF

          !skip GEMM calls if group has to wait for other groups 
       IF ( .NOT. LDO_TAU_LOCAL ) THEN
          CYCLE  nblock
       ENDIF
#ifdef gammareal
       ! construct all relevant Green's functions
       ! interpret coefficients as sine and cosine transforms
       ! we need all "pairs" as real final coefficients
       CALL __DGEMM__('N','T', GDES%NRPLWV_ROW_DATA_POINTS, GDES%NRPLWV_COL_DATA_POINTS, NGLB, 1.0_q, &
            CG(1,1), 2*SIZE(CG,1), CG_COLUMN(1,1), 2*SIZE(CG_COLUMN,1), &
            1.0_q, G%GG(1,1), SIZE(G%GG, 1))
       
       ! missuse DGEMM first array is complex (plane wave) result is also complex
       IF (GDES%NPRO_COL/=0) THEN
          CALL __DGEMM__('N','T', GDES%NRPLWV_ROW_DATA_POINTS, GDES%NPRO_COL, NGLB, 1.0_q, &
               CG(1,1), 2*SIZE(CG,1), CPROJ_COLUMN(1,1), SIZE(CPROJ_COLUMN,1), &
               1.0_q, G%G_PROJ(1,1), SIZE(G%G_PROJ, 1))
       ENDIF

#else
       CALL __ZGEMM__('N','C', GDES%NRPLWV_ROW, GDES%NRPLWV_COL, NGLB, one, &
            CG(1,1), SIZE(CG,1), CG_COLUMN(1,1), SIZE(CG_COLUMN,1), &
            one, G%GG(1,1), SIZE(G%GG, 1))
       

       IF (GDES%NPRO_COL/=0) THEN
          CALL __ZGEMM__('N','C', GDES%NRPLWV_ROW, GDES%NPRO_COL, NGLB, one, &
               CG(1,1), SIZE(CG,1), CPROJ_COLUMN(1,1), SIZE(CPROJ_COLUMN,1), &
               one, G%G_PROJ(1,1), SIZE(G%G_PROJ, 1))
       ENDIF
#endif
       IF (GDES%NPRO_COL/=0) THEN
          CALL __GGEMM__('N',trans, GDES%NPRO_ROW, GDES%NPRO_COL, NGLB, one, &
               CPROJ(1,1), SIZE(CPROJ,1), CPROJ_COLUMN(1,1), SIZE(CPROJ_COLUMN,1), &
               one, G%PROJ_PROJ(1,1), SIZE(G%PROJ_PROJ, 1))
       ENDIF
    ENDDO nblock
   
    CALL DELWAV(W1, .TRUE.)    
    DEALLOCATE(CG, CPROJ, CG_COLUMN, CPROJ_COLUMN, WEIGHT)

    PROFILING_STOP('calculate_g_reciprocal')

  END SUBROUTINE CALCULATE_G_RECIPROCAL

!***********************************************************************
!
! the following subroutine is similar to DIS_PROJ
! it copies the relevant part of the orbital character CPROJ
! from a globally known array CPROJ to the one treated
! locally CPROJ_LOCAL
!
!***********************************************************************

   SUBROUTINE GDIS_PROJ(GDES, CPROJ, CPROJ_LOCAL)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (greensfdes) GDES
      GDEF CPROJ(:,:)
      GDEF CPROJ_LOCAL(:,:)
#ifdef MPI
      INTEGER NIS,NPRO,NT,LMMAXC,NPRO_POS,L,NI,NPRO_ROW
!
! quick copy if possible
      IF (GDES%COMM%NCPU==1) THEN
         CPROJ_LOCAL(1:GDES%NPRO_COL,:)=CPROJ(1:GDES%NPRO_COL,:)
         RETURN
      ENDIF

      NIS =1
      NPRO=0
      DO NT=1,GDES%NTYP
        LMMAXC=GDES%LMMAX(NT)
        IF (LMMAXC/=0) THEN
           DO NI=NIS,GDES%NITYP(NT)+NIS-1
              NPRO_POS=GDES%NPRO_POS(NI,GDES%COMM%NODE_ME)
              DO L=1,LMMAXC
                 CPROJ_LOCAL(L+NPRO,:)=CPROJ(L+NPRO_POS,:)
              ENDDO
              NPRO= LMMAXC+NPRO
           ENDDO
        ENDIF
        NIS = NIS+GDES%NITYP(NT)
      ENDDO

      IF (GDES%LNONCOLLINEAR) THEN
         NIS =1
         NPRO_ROW= GDES%NPRO_ROW/2

         DO NT=1,GDES%NTYP
            LMMAXC=GDES%LMMAX(NT)
            IF (LMMAXC/=0) THEN
               DO NI=NIS,GDES%NITYP(NT)+NIS-1
                  NPRO_POS=GDES%NPRO_POS(NI,GDES%COMM%NODE_ME)
                  DO L=1,LMMAXC
                     CPROJ_LOCAL(L+NPRO,:)=CPROJ(L+NPRO_POS+NPRO_ROW,:)
                  ENDDO
                  NPRO= LMMAXC+NPRO
               ENDDO
            ENDIF
            NIS = NIS+GDES%NITYP(NT)
         ENDDO
      ENDIF
#else
      CPROJ_LOCAL(1:GDES%NPRO_COL,:)=CPROJ(1:GDES%NPRO_COL,:)
#endif
    END SUBROUTINE GDIS_PROJ


!***********************************************************************
!
! reduce the Green's function or self-energy storage
! by compressing along the second storage index
!
! the second routine restores the Green's function using
! inversion symmetry 
!
! the relevant relation is G(g+k,g'+k)=G(-g-k,-g'-k)
! therefore compression is possible along either the first
! or the second index
! currently we use compression along the second index
! albeit those are stored distributed
!
! TODO storage on temporary arrayes needs to be reduced
!
!***********************************************************************

  SUBROUTINE COMPRESS_G_RECIPROCAL( G, GDES, NK1)
    TYPE (greensf) :: G
    TYPE (greensfdes) :: GDES  
    INTEGER :: NK1             ! considered k-point
    GDEF, ALLOCATABLE :: GG(:,:)
    INTEGER :: N,M             ! row and column index

    IF (ASSOCIATED(GDES%LUSEINV)) THEN
    IF (GDES%LUSEINV(NK1)) THEN
       ALLOCATE(GG(GDES%NRPLWV_ROW_DATA_POINTS, GDES%NRPLWV_COL_DATA_POINTS))
       CALL TRANSPOSE_G_G(G%GG, GG, GDES )

       ! now compress data by using inversion symmetry on the first index
       DO M=1,GDES%NRPLWV_COL_DATA_POINTS
          DO N=1,GDES%NGVECTOR_INV(NK1) 
             ! use out of place compression to make sure no dependencies emerge
             G%GG(N,M)=GG(GDES%MAP_TO_FULL(N,1,NK1),M)
          ENDDO
       ENDDO
       ! copy back
       GG=G%GG
       ! transpose again and use compressed storage mode along first dimension
       ! as a result G%GG is compressed along second dimension
       CALL TRANSPOSE_G_G(GG, G%GG, GDES, NK1=NK1 )

       DEALLOCATE(GG)

    ENDIF
    ENDIF
    
  END SUBROUTINE COMPRESS_G_RECIPROCAL

  SUBROUTINE UNCOMPRESS_G_RECIPROCAL( G, GDES, NK1)
    TYPE (greensf) :: G
    TYPE (greensfdes) :: GDES  
    INTEGER :: NK1             ! considered k-point
    GDEF, ALLOCATABLE :: GG1(:,:)
    GDEF, ALLOCATABLE :: GG2(:,:)
    INTEGER :: N,M             ! row and column index

    IF (ASSOCIATED(GDES%LUSEINV)) THEN
    IF (GDES%LUSEINV(NK1)) THEN
       ALLOCATE(GG1(GDES%NRPLWV_ROW_DATA_POINTS, GDES%NRPLWV_COL_DATA_POINTS))
       ALLOCATE(GG2(GDES%NRPLWV_ROW_DATA_POINTS, GDES%NRPLWV_COL_DATA_POINTS))

       CALL TRANSPOSE_G_G(G%GG, GG1, GDES, NK2=NK1 )
       ! now apply inversion symmetry on first index
       GG2=G%GG
       DO M=1,GDES%NRPLWV_COL_DATA_POINTS
          DO N=1,GDES%NGVECTOR_INV(NK1) 
             ! interchange inverted and un-inverted rows
             G%GG(GDES%MAP_TO_FULL(N,1,NK1),M)=GG2(GDES%MAP_TO_FULL(N,2,NK1),M)
             G%GG(GDES%MAP_TO_FULL(N,2,NK1),M)=GG2(GDES%MAP_TO_FULL(N,1,NK1),M)
          ENDDO
       ENDDO
       CALL TRANSPOSE_G_G(G%GG, GG2, GDES, NK2=NK1 )
       ! now un compress data by using inversion symmetry on the first index
       DO M=1,GDES%NRPLWV_COL_DATA_POINTS
          DO N=1,GDES%NGVECTOR_INV(NK1) 
             ! use out of place compression to make sure no dependencies emerge
             G%GG(GDES%MAP_TO_FULL(N,2,NK1),M)=GCONJG(GG2(N,M))
          ENDDO
       ENDDO
       ! now without inversion symmetry
       DO M=1,GDES%NRPLWV_COL_DATA_POINTS
          DO N=1,GDES%NGVECTOR_INV(NK1) 
             ! use out of place compression to make sure no dependencies emerge
             G%GG(GDES%MAP_TO_FULL(N,1,NK1),M)=GG1(N,M)
          ENDDO
       ENDDO

       ! copy back
       GG1=G%GG
       ! transpose again
       CALL TRANSPOSE_G_G(GG1, G%GG, GDES )

       DEALLOCATE(GG1)
       DEALLOCATE(GG2)

    ENDIF
    ENDIF
    
  END SUBROUTINE UNCOMPRESS_G_RECIPROCAL

!***********************************************************************
!
! determine the selfenergy  sigma_k2 = W_q GO_k1
! where k2 must equal q + k1 (upon entry)
!
! subroutine allocates SIGMA(:)%PROJ_PROJ, %G_PROJ and %GG
! if LADD=.FALSE.
! for LADD=.TRUE., the corresponding arrays must be allocated prior
! to calling the routine
!
!***********************************************************************

  SUBROUTINE CALCULATE_SIGMA_TAU( LMDIM, W, WGW, GO, GU, SIGMA, GDES, & 
       CHI, NQ, NK1, NK2, LATT_CUR, IU6, TAU_WEIGHT, FRNL, LFORCE )
    USE ini
    USE dfast
    USE fock
    USE constant
    TYPE (wavespin):: W
    TYPE (wavedes) :: WGW
    TYPE (greensf) :: GO         ! Green function for occupied orbitals at NK1
    TYPE (greensf) :: GU         ! Green function for unoccupied orbitals at NK2
    TYPE (greensf) :: SIGMA      ! self energy at NK2
    TYPE (greensfdes) :: GDES   
    TYPE (responsefunction) CHI  ! response function
    INTEGER :: NQ, NK1, NK2
    TYPE (latt)        LATT_CUR
    INTEGER :: IU6
    REAL(q) :: TAU_WEIGHT
    REAL(q) :: FRNL(:,:)! forces related to rigid shift of augmentation charges
    LOGICAL :: LFORCE
  ! local
    TYPE (wavedes1) :: WDESQ
    TYPE (wavedes1) :: WDES1, WDES2
    INTEGER :: NRP, NT, NI, NIS, LMMAXC, L, LP, LM, LMBASE, NAUG, NT_GLOBAL
    INTEGER :: NK1_IN_KPOINTS_FULL_ORIG
    LOGICAL :: LPHASE
    COMPLEX(q) :: CPHASE(GRIDHF%MPLWV)

    GDEF, ALLOCATABLE :: GWORK(:,:)
    GDEF, ALLOCATABLE,TARGET :: CRHOLM(:) ! augmentation occupancy matrix
    GDEF, ALLOCATABLE,TARGET :: P_NLM(:,:)
    GDEF, ALLOCATABLE :: P_G_PROJ(:,:,:)
    GDEF, ALLOCATABLE,TARGET :: P_PROJ_G(:,:,:)
    GDEF, ALLOCATABLE :: P_G_R(:,:)
    GDEF, ALLOCATABLE :: P_R_G(:,:)
    TYPE (wavefun1) :: WTMP
    INTEGER LMDIM
    GDEF            :: CDIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    REAL(q)         :: SCALE
!force
    TYPE (greensf)  :: SIGMA_TMP  ! required to temporarily store PROJ_PROJ and G_PROJ
    INTEGER         :: NDIR, IDIR
    TYPE (nonlr_struct), ALLOCATABLE :: FAST_AUG(:)
!force
    COMPLEX(q), ALLOCATABLE :: CWORK1(:),CWORK2(:)
    GDEF,       ALLOCATABLE :: CPROJ(:)

    REAL(q) :: RINPL
    INTEGER :: ISTAT 
!$  INTEGER :: LASTTYP

    PROFILING_START('calculate_sigma_tau')

    ALLOCATE(CRHOLM(GDES%NLM_ROW))
    CALL SETWDES(W%WDES,WDES1,NK1)
    CALL SETWDES(W%WDES,WDES2,NK2)

!force
    ! to seek further changes seek IDIR
    IF (LFORCE) THEN
       NDIR=4
    ELSE
       NDIR=1
    ENDIF
    CALL SETUP_FAST_AUG_NDIR
!force

    NK1_IN_KPOINTS_FULL_ORIG=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,NK1),KPOINTS_FULL_ORIG)

    CALL PHASER_HF(GRIDHF, LATT_CUR, FAST_AUG_FOCK, -(KPOINTS_FULL%VKPT(:,NK1)-KPOINTS_FULL%VKPT(:,NK2)))

    ! gK: note we take the negative phase here compared to calculation of response function
    CALL SETPHASE(-(W%WDES%VKPT(:,NK2)-W%WDES%VKPT(:,NK1)-CHI%VKPT(:)), &
         GRIDHF,CPHASE(1),LPHASE) 

    ALLOCATE(CWORK1(WDES1%GRID%MPLWV*WDES1%NRSPINORS))
    ALLOCATE(CWORK2(WDES2%GRID%MPLWV*WDES2%NRSPINORS),CPROJ(WDES2%NPROD))
    ! this is the q-point for the screened exchange W q=k2-k1
    CALL SETWDES(WGW, WDESQ, NQ )

    IF (WDES1%GRID%RL%NP/= WDESQ%GRID%RL%NP) THEN
       WRITE(*,*) 'internal error in CALCULATE_SIGMA_TAU: grids distinct ',WDES1%GRID%RL%NP, WDESQ%GRID%RL%NP
    ENDIF
    IF (WDES1%NPRO /= WDESQ%NPRO) THEN
       WRITE(*,*) 'internal error in CALCULATE_SIGMA_TAU: NPRO distinct ',WDES1%NPRO /= WDESQ%NPRO
    ENDIF
 
     ! FFT GO to real space and set entries G_R, R_PROJ, PROJ_R
     ! passed down that way now
     CALL FFT_G(W%WDES, GO, GDES, NK1)
!=======================================================================
! FFT of W(G,G') along first direction 
! yielding W(r,G') and W(NLM, G')
! then transpose and conjugate to obtain W(G',r ) and W(G', NLM)
!=======================================================================
       PROFILING_START('calculate_sigma_tau1')

       ALLOCATE(GWORK(m_ WDESQ%GRID%MPLWV,1))
!force
       ALLOCATE(P_PROJ_G (GDES%NLM_ROW,    GDES%RES_NRPLWV_COL_DATA_POINTS, NDIR),STAT=ISTAT)
       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "CALCULATE_SIGMA_TAU (P_PROJ_G) is not able to allocate "//&
            str( ndata * 8._q*GDES%NLM_ROW*GDES%RES_NRPLWV_COL_DATA_POINTS*NDIR/1024) //&
            " kB of data on MPI rank 0." )
       ENDIF
!force
       ALLOCATE(P_R_G    (GDES%MPLWV_ROW, GDES%RES_NRPLWV_COL_DATA_POINTS),STAT=ISTAT)
       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "CALCULATE_SIGMA_TAU (P_R_G) is not able to allocate "//&
            str( ndata * 8._q*GDES%MPLWV_ROW*GDES%RES_NRPLWV_COL_DATA_POINTS/1024) //&
            " kB of data on MPI rank 0." )
       ENDIF

       CALL REGISTER_ALLOCATE(wsgf*(SIZE(P_R_G,KIND=qi8)+SIZE(P_PROJ_G,KIND=qi8)), "GGwork")

       IF (WDES1%LOVERL) AUG_DES%RINPL=1._q/WDESQ%GRID%NPLWV

     ! loop over second index G'
!!!$OMP PARALLEL DO PRIVATE(NRP,GWORK,IDIR,WTMP)
       DO NRP=1,GDES%RES_NRPLWV_COL_DATA_POINTS
          IF (CHI%LREALSTORE) THEN
             ! FFT of W(G,G') along first direction -> GWORK(1,1)
             CALL FFTWAV(WDESQ%NGVECTOR,WDESQ%NINDPW(1), &
                  GWORK(1,1), &
                  CHI%RESPONSER(1, NRP, 1),WDESQ%GRID)

             IF (LPHASE) CALL APPLY_PHASE_GDEF( GRIDHF, CPHASE(1), GWORK(1,1) , GWORK(1,1))

          ELSE
             ! FFT along first direction -> GWORK(1,1)
             CALL FFTWAV(WDESQ%NGVECTOR,WDESQ%NINDPW(1), &
                  GWORK(1,1), &
                  CHI%RESPONSEFUN(1, NRP, 1),WDESQ%GRID)

             IF (LPHASE) CALL APPLY_PHASE_GDEF( GRIDHF, CPHASE(1), GWORK(1,1) , GWORK(1,1))
          ENDIF

          ! add augmentation part  P(NLM, r) Q_NLM(r')
          IF (WDES1%LOVERL) THEN
             DO IDIR=1,NDIR
                ! int d^3 r Q_LM(r) W(r,G)
                WTMP%CPROJ=> P_PROJ_G(:, NRP, IDIR)   ! RPRO1_HF only dereferences WTMP%CPROJ
                CALL RPRO1_HF(FAST_AUG(IDIR), AUG_DES, WTMP, GWORK(1,1))
             ENDDO
          END IF

!!          P_R_G(1:WDESQ%GRID%RL%NP,NRP)=GWORK(1:WDESQ%GRID%RL%NP,1) ! or (1._q/WDESQ%GRID%NPLWV)
          CALL OMP_DCOPY(ndata*WDESQ%GRID%RL%NP,GWORK(1,1),1,P_R_G(1,NRP),1)
       ENDDO
!!!$OMP END PARALLEL DO
       ! transpose P_R_G    to P_G_R
       ! and       P_PROJ_G to P_G_PROJ
!force
       ALLOCATE(P_G_PROJ (GDES%RES_NRPLWV_ROW_DATA_POINTS, GDES%NLM_COL, NDIR),STAT=ISTAT)
       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "CALCULATE_SIGMA_TAU (P_G_PROJ) is not able to allocate "//&
            str( ndata * 8._q*GDES%RES_NRPLWV_ROW_DATA_POINTS*GDES%NLM_COL*NDIR/1024) //&
            " kB of data on MPI rank 0." )
       ENDIF
!force
       CALL REGISTER_ALLOCATE(wsgf*SIZE(P_G_PROJ,KIND=qi8), "GGwork")

       ALLOCATE(P_G_R    (GDES%RES_NRPLWV_ROW_DATA_POINTS, GDES%MPLWV_COL))
       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "CALCULATE_SIGMA_TAU (P_G_R) is not able to allocate "//&
            str( ndata * 8._q*GDES%RES_NRPLWV_ROW_DATA_POINTS*GDES%MPLWV_COL/1024) //&
            " kB of data on MPI rank 0." )
       ENDIF
       CALL REGISTER_ALLOCATE(wsgf*SIZE(P_G_R,KIND=qi8), "GGwork")

       CALL DUMP_ALLOCATE_TAG(IU6,"SIGMA_TAU"); IF (IU6>=0) CALL WFORCE(IU6) ! here we reach the first memory bottleneck

       CALL TRANSPOSE_R_G_RESPONSE   (P_R_G, P_G_R,  GDES)
       IF (WDES1%LOVERL) THEN
          DO IDIR=1,NDIR
             CALL TRANSPOSE_PROJ_G_RESPONSE(P_PROJ_G(:,:,IDIR), P_G_PROJ(:,:,IDIR),  GDES)
          ENDDO
       ENDIF
       CALL DEREGISTER_ALLOCATE(wsgf*(SIZE(P_R_G, KIND=qi8 )+SIZE(P_PROJ_G, KIND=qi8)), "GGwork")
       DEALLOCATE( P_R_G, P_PROJ_G)

       PROFILING_STOP('calculate_sigma_tau1')
!=======================================================================
! FFT of W(G', r) along second direction (now again the first index
! since the matrix has been transposed)
!  sigma(G', r) -> SIGMA%G_R
! optional (precompiler): sigma(alpha, r) -> SIGMA%PROJ_R
!=======================================================================
       PROFILING_START('calculate_sigma_tau2')

      ! allocate SIGMA 
       CALL ALLOCATE_G_R(GDES, SIGMA)
       CALL DUMP_ALLOCATE_TAG(IU6,"SIGMA_TAU2"); IF (IU6>=0) CALL WFORCE(IU6) ! here things might get tight again

  ! loop over all space points r
!!!$OMP PARALLEL DO PRIVATE(NRP,GWORK,CWORK1,CWORK2)
       DO NRP=1,GDES%MPLWV_COL
          ! FFT of W(G', r) -> W(r' ,r)
          CALL FFTWAV(WDESQ%NGVECTOR,WDESQ%NINDPW(1), &
               GWORK(1,1), &
               P_G_R(1,  NRP),WDESQ%GRID )
          IF (LPHASE) CALL APPLY_PHASE_GDEF( GRIDHF, CPHASE(1), GWORK(1,1), GWORK(1,1))
          ! FFT for Greens function G(G, r) -> G(r', r)
          CALL FFTWAV(WDES1%NGVECTOR, WDES1%NINDPW(1), &
               CWORK1(1), &
               GO%G_R(1,NRP),WDES1%GRID)
          
          ! Sigma(r',r) = W(r',r) G(r',r)  for all r
!!          CWORK2(1:WDES2%GRID%RL%NP)=CWORK1(1:WDES2%GRID%RL%NP)*GWORK(1:WDES2%GRID%RL%NP,1)*(1.0_q/WDES2%GRID%NPLWV)
          RINPL=1.0_q/WDES2%GRID%NPLWV
#ifndef gammareal
          CALL CMPLX_CMPLX_CMPLX_MUL(WDES2%GRID%RL%NP,.FALSE.,GWORK(1,1),CWORK1(1),RINPL,CWORK2(1),0._q,CWORK2(1))
#else
          CALL REAL_CMPLX_CMPLX_MUL (WDES2%GRID%RL%NP,GWORK(1,1),CWORK1(1),RINPL,CWORK2(1),0._q,CWORK2(1))
#endif
          ! FFT to reciprocal space -> Sigma(G, r)
          CALL FFTEXT(WDES2%NGVECTOR, WDES2%NINDPW(1), &
               CWORK2(1), &
               SIGMA%G_R(1,NRP),WDES2%GRID, .FALSE.)
       ENDDO
!!!$OMP END PARALLEL DO
       CALL DEREGISTER_ALLOCATE(wsgf*SIZE(P_G_R, KIND=qi8), "GGwork")
       DEALLOCATE(P_G_R)

       DEALLOCATE(GWORK)

       PROFILING_STOP('calculate_sigma_tau2')
!=======================================================================
! now the akward part dealing with P_G_PROJ =  W(G', NLM)
!=======================================================================
       IF (W%WDES%LOVERL) THEN
       PROFILING_START('calculate_sigma_tau3')

       IF (.NOT. ASSOCIATED(SIGMA%PROJ_PROJ)) THEN
          CALL vtutor%bug("internal error in CALCULATE_SIGMA_TAU: SIGMA%PROJ_PROJ is not associated", __FILE__, __LINE__)
       ENDIF

       CALL ALLOCATE_R_PROJ(GDES, SIGMA)

       CALL ALLOCATE_PROJ_PROJ(GDES, SIGMA_TMP)
       CALL ALLOCATE_G_PROJ(GDES, SIGMA_TMP)

       dir: DO IDIR=1,NDIR

    ! zero, since we add contributions
       SIGMA%R_PROJ=0
       SIGMA_TMP%G_PROJ=0
       SIGMA_TMP%PROJ_PROJ=0

       AUG_DES%RINPL=1.0_q/WDESQ%GRID%NPLWV

    ! loop over types that are stored locally (second = column index)
       LMBASE =0
       NAUG=0
#ifndef _OPENMP
       NIS =1
       DO NT=1,GDES%NTYP  ! note NT is a local type index
                          ! and might not span all atomic types
          NT_GLOBAL=GDES%NT_GLOBAL(NT)  ! global type index

          ! WDESQ%GRID%MPLWV is the number of complex words required to do in place FFT
          ALLOCATE(P_NLM(GDES%NLM_ROW,GDES%NLM_LMMAX(NT)),  & 
                   GWORK(m_ WDESQ%GRID%MPLWV,GDES%NLM_LMMAX(NT)))

          IF (AUG_DES%LMMAX(NT_GLOBAL)/= GDES%NLM_LMMAX(NT)) THEN
             CALL vtutor%bug("internal error 1 in CALCULATE_SIGMA_TAU: " // &
                str(AUG_DES%LMMAX(NT_GLOBAL)) // " " // str(GDES%NLM_LMMAX(NT)), __FILE__, __LINE__)
          ENDIF

          LMMAXC=GDES%LMMAX(NT)
          IF (LMMAXC/=0) THEN
             !
             ! main loop over local ion index N (second or column index)
             !
             P_NLM=0
             GWORK=0

             DO NI=NIS,GDES%NITYP(NT)+NIS-1
#else
       LASTTYP=0
!!!$OMP PARALLEL DO PRIVATE(NI,NT,NT_GLOBAL,P_NLM,GWORK,LMMAXC,LM,WTMP,L,LP,CRHOLM,CDIJ,CPROJ,RINPL,CWORK2) &
!!!$OMP FIRSTPRIVATE(LASTTYP,LMBASE,NAUG)
       DO NI=1,GDES%NIONS

          NT=GDES%NT_LOCAL(NI) ; NT_GLOBAL=GDES%NT_GLOBAL(NT)

          IF (NT /= LASTTYP) THEN
             IF (AUG_DES%LMMAX(NT_GLOBAL)/= GDES%NLM_LMMAX(NT)) THEN
                CALL vtutor%bug("internal error 1 in CALCULATE_SIGMA_TAU: " // &
                   str(AUG_DES%LMMAX(NT_GLOBAL)) // " " // str(GDES%NLM_LMMAX(NT)), __FILE__, __LINE__)
             ENDIF

             IF (ALLOCATED(P_NLM)) DEALLOCATE(P_NLM)
             IF (ALLOCATED(GWORK)) DEALLOCATE(GWORK)

             ALLOCATE(P_NLM(GDES%NLM_ROW,GDES%NLM_LMMAX(NT)),  & 
                      GWORK(m_ WDESQ%GRID%MPLWV,GDES%NLM_LMMAX(NT)))

             LASTTYP=NT
          ENDIF

          LMMAXC=GDES%LMMAX(NT)

          IF (LMMAXC/=0) THEN
             IF (NI>1) THEN
                LMBASE=SUM(GDES%NPRO_ENTRIES(1:NI-1,GDES%COMM%NODE_ME))
                NAUG  =SUM(GDES%NLM_ENTRIES (1:NI-1,GDES%COMM%NODE_ME))
             ENDIF

!            P_NLM=0 ; GWORK=0
#endif
                ! extract data using wave function FFT (impose cutoff at the same time)
                !  W(r', NLM) = \sum_G' W(G', NLM) e (-i G' r')
                DO LM=1,GDES%NLM_LMMAX(NT)
                   CALL FFTWAV(WDESQ%NGVECTOR,WDESQ%NINDPW(1), &
                        GWORK(1,LM), &
                        P_G_PROJ(1, NAUG+LM, IDIR),WDESQ%GRID)

                   IF (NAUG+LM >  GDES%NLM_COL) THEN
                      CALL vtutor%bug("internal error 2 in CALCULATE_SIGMA_TAU: " // str(NAUG+LM) &
                         // " " // str(GDES%NLM_COL), __FILE__, __LINE__)
                   ENDIF
                   IF (LPHASE) CALL APPLY_PHASE_GDEF( GRIDHF, CPHASE(1), GWORK(1,LM), GWORK(1,LM))
                ENDDO

                ! for all LM int d^3 r Q_N'L'M'(r') W(r',NLM) = W(N'L'M', NLM) -> P_NLM
                DO LM=1,GDES%NLM_LMMAX(NT)
                   WTMP%CPROJ=> P_NLM(:,LM)
                   CALL RPRO1_HF(FAST_AUG_FOCK, AUG_DES, WTMP, GWORK(1,LM))
                ENDDO

                DO L=1,LMMAXC   ! alpha
                DO LP=1,LMMAXC  ! beta
                   CRHOLM=0
                   ! Clebsch Gordon transformation for  *second* index
                   ! for all N'L'M' calculate D(N'L'M') = \sum_LM W(N'L'M', NLM)  T_alpha, beta^LM
                   DO LM=1,GDES%NLM_LMMAX(NT)
                      ! index build from NI and LM
                      CRHOLM(:)=CRHOLM(:)+P_NLM(:,LM)*TRANS_MATRIX_FOCK(LP,L,LM,NT_GLOBAL)
                   ENDDO

                   ! Clebsch Gordon transformation over first index ->  D(N' L'M') -> D(alpha'beta')
                   CALL CALC_DLLMM_TRANS(W%WDES, AUG_DES, TRANS_MATRIX_FOCK, CDIJ(:,:,:,:),CRHOLM)

                   ! sum_alpha',beta' D(alpha'beta') x G(alpha', alpha) 
                   CALL OVERL_FOCK(W%WDES, SIZE(CDIJ,1), CDIJ(1,1,1,1),GO%PROJ_PROJ(1,LMBASE+L),CPROJ(1), .FALSE.)
                   ! add to Sigma(beta', beta )
                   SIGMA_TMP%PROJ_PROJ(1:WDES2%NPRO,LMBASE+LP)=SIGMA_TMP%PROJ_PROJ(1:WDES2%NPRO,LMBASE+LP)+CPROJ(1:WDES2%NPRO)
                ENDDO
                ENDDO

                CALL CALC_DLLMM_VECTOR( SIGMA%R_PROJ, GO%R_PROJ, GWORK, GDES%MPLWV_ROW, TRANS_MATRIX_FOCK(:,:,:,NT_GLOBAL), LMMAXC, LMBASE, GDES%NLM_LMMAX(NT) )

                ! FFT SIGMA%R_PROJ and add result to  SIGMA%G_PROJ
                DO LP=1,LMMAXC
!!                   CWORK2(1:WDES2%GRID%RL%NP)=SIGMA%R_PROJ(1:WDES2%GRID%RL%NP, LP+LMBASE)*(1.0_q/WDES2%GRID%NPLWV)
                   RINPL=1.0_q/WDES2%GRID%NPLWV
#ifndef gammareal
                   CALL OMP_DSCPY(2*WDES2%GRID%RL%NP,RINPL,SIGMA%R_PROJ(1,LP+LMBASE),1,CWORK2(1),1)
#else
                   CALL REAL_CMPLX_SCALE(WDES2%GRID%RL%NP,SIGMA%R_PROJ(1,LP+LMBASE),RINPL,CWORK2(1))
#endif
                   CALL FFTEXT(WDES2%NGVECTOR,WDES2%NINDPW(1), &
                        CWORK2, SIGMA_TMP%G_PROJ(1,LP+LMBASE), &
                        WDES2%GRID,.TRUE.)
                ENDDO
#ifndef _OPENMP
                LMBASE = LMMAXC+LMBASE
                NAUG = NAUG+ GDES%NLM_LMMAX(NT)
             ENDDO
          ENDIF

          NIS = NIS+GDES%NITYP(NT)
          DEALLOCATE(P_NLM, GWORK)
#else
          ENDIF
#endif
       ENDDO
!!!$OMP END PARALLEL DO
!$     IF (ALLOCATED(P_NLM)) DEALLOCATE(P_NLM)
!$     IF (ALLOCATED(GWORK)) DEALLOCATE(GWORK)

       ! get non-local contributions to the energy
       IF (LFORCE) &
            CALL CONTRACT_SIGMA_G_TAU_NONL( WDES2, KPOINTS_ORIG%WTKPT(NK2)*TAU_WEIGHT, & 
            GU, SIGMA_TMP, GDES, FRNL, IDIR )
       ENDDO dir

       ! add results from final dir loop to SIGMA (corresponds to undisplaced augmentation charges)
       SIGMA%PROJ_PROJ=SIGMA%PROJ_PROJ+SIGMA_TMP%PROJ_PROJ
       SIGMA%G_PROJ   =SIGMA%G_PROJ   +SIGMA_TMP%G_PROJ

       CALL DEALLOCATE_PROJ_PROJ(SIGMA_TMP)
       CALL DEALLOCATE_G_PROJ(SIGMA_TMP)

       CALL DEALLOCATE_R_PROJ(SIGMA)

       PROFILING_STOP('calculate_sigma_tau3')
   
       ENDIF

       CALL DEREGISTER_ALLOCATE(wsgf*SIZE(P_G_PROJ,KIND=qi8), "GGwork")
       DEALLOCATE(P_G_PROJ )

    ! deallocate all objects allocated in initial part
    CALL RELEASE_FFT_G( GO)
    DEALLOCATE(CWORK1,CWORK2,CPROJ)
    DEALLOCATE(CRHOLM)
    CALL FREE_FAST_AUG_NDIR()
    CALL FFT_SIGMA ( W%WDES, SIGMA, GDES, NK2 )

    PROFILING_STOP('calculate_sigma_tau')

    RETURN

    CONTAINS
!==========================================================================
! setup first derivative of augmentation charges with respect
! to ionic positions
!==========================================================================

      SUBROUTINE SETUP_FAST_AUG_NDIR
        INTEGER :: IDIR
        REAL(q) :: DIS=fd_displacement
        
        REAL(q)    :: DISPL1(3,FAST_AUG_FOCK%NIONS),DISPL2(3,FAST_AUG_FOCK%NIONS)

        IF (.NOT.WDES1%LOVERL) RETURN

        ALLOCATE(FAST_AUG(NDIR))
        ! last entry is identical to FAST_AUG_FOCK
!       FAST_AUG(NDIR)=FAST_AUG_FOCK
        CALL NONLR_ASSIGN(FAST_AUG(NDIR),FAST_AUG_FOCK)

        ! set first derivative of in first 3 entries
        DO IDIR=1,NDIR-1
           CALL COPY_FASTAUG(FAST_AUG_FOCK, FAST_AUG(IDIR))
           DISPL1=0
           DISPL1(IDIR,:)=-DIS
           
           DISPL2=0
           DISPL2(IDIR,:)=+DIS
           CALL RSPHER_ALL(GRIDHF, FAST_AUG(IDIR),LATT_CUR,LATT_CUR,LATT_CUR, DISPL1, DISPL2,1)
           ! phase factor identical to FAST_AUG_FOCK
           IF (ASSOCIATED(FAST_AUG(IDIR)%CRREXP)) DEALLOCATE(FAST_AUG(IDIR)%CRREXP)
           FAST_AUG(IDIR)%CRREXP=>FAST_AUG_FOCK%CRREXP
#ifdef shmem_rproj
           IF (FAST_AUG(IDIR)%COMM_shmem%NODE_ME==1) FAST_AUG(IDIR)%RPROJ=FAST_AUG(IDIR)%RPROJ*SQRT(LATT_CUR%OMEGA)*(1._q/(2*DIS))
           CALLMPI( M_barrier(FAST_AUG(IDIR)%COMM_shmem))
#else
           FAST_AUG(IDIR)%RPROJ=FAST_AUG(IDIR)%RPROJ*SQRT(LATT_CUR%OMEGA)*(1._q/(2*DIS))
#endif
        ENDDO
      END SUBROUTINE SETUP_FAST_AUG_NDIR
      
      SUBROUTINE FREE_FAST_AUG_NDIR
        IF (.NOT.WDES1%LOVERL) RETURN

        DO IDIR=1,NDIR-1
           NULLIFY(FAST_AUG(IDIR)%CRREXP)
           CALL NONLR_DEALLOC(FAST_AUG(IDIR))
        ENDDO
        DEALLOCATE(FAST_AUG)
      END SUBROUTINE FREE_FAST_AUG_NDIR

  END SUBROUTINE CALCULATE_SIGMA_TAU
!***********************************************************************
!
! contract nonlocal part of selfenergy and Green's function G
! this is required for force calculations
! similar to previous routine
!
!***********************************************************************

  SUBROUTINE CONTRACT_SIGMA_G_TAU_NONL(WDES2, WEIGHT, GU, SIGMA, GDES, FRNL, IDIR )
    USE ini
    USE dfast
    USE fock
    USE constant
    TYPE (wavedes1):: WDES2
    REAL(q) :: WEIGHT
    TYPE (greensf) :: GU         ! Green function for unoccupied orbitals 
    TYPE (greensf) :: SIGMA      ! Green function for occupied orbitals
    TYPE (greensfdes) :: GDES
    REAL(q) :: FRNL(:,:)
    INTEGER  :: IDIR             ! component of force to be calculated
! local
    COMPLEX(q)      :: CSUM
    INTEGER :: ICOL, NG, NGP, NPRO
    INTEGER :: LMMAXC, NIS, NI, NIP, NT, LM
    COMPLEX(q)      :: CFOR(SIZE(FRNL,2))

! set k-point descriptor
    CFOR=0

    NIS=1
    ICOL=0
    DO NT=1,GDES%NTYP
       LMMAXC=GDES%LMMAX(NT)
       IF (LMMAXC>0) THEN
          DO NI=NIS,GDES%NITYP(NT)+NIS-1
             NIP=GDES%NI_GLOBAL(NI) ! global index of ion 
             DO LM=1, LMMAXC
                ! only use actual data points
                DO NG=1,m_ WDES2%NGVECTOR
                   CFOR(NIP)=CFOR(NIP)+SIGMA%G_PROJ(NG,ICOL+LM)*GCONJG( GU%G_PROJ(NG,ICOL+LM))
                ENDDO
                ! only use actual data points
                DO NPRO=1,WDES2%NPRO ! same as GDES%NPRO_ROW
                   CFOR(NIP)=CFOR(NIP)+SIGMA%PROJ_PROJ(NPRO,ICOL+LM)*GCONJG( GU%PROJ_PROJ(NPRO,ICOL+LM))
                ENDDO
             ENDDO
          ICOL=ICOL+LMMAXC
          ENDDO
       ENDIF
       NIS = NIS+GDES%NITYP(NT)
    ENDDO

    ! sum over current time point
    CALLMPI( M_sum_z(GDES%COMM, CFOR , SIZE(CFOR)))
    CFOR=CFOR*WEIGHT

    CSUM=SUM(CFOR)
#ifdef MPI
    IF (GDES%COMM%NODE_ME==GDES%COMM%IONODE .AND. IDIR<=SIZE(FRNL,1) ) THEN
#endif
!       WRITE(*,'(A,I2,2F20.10)') ' CONTRACT_SIGMA_G_TAU_NONL     ',IDIR,CSUM
!       WRITE(*,'(8F10.6)') CFOR
#ifdef MPI
    ENDIF
#endif
    ! add force to FRNL
    IF (IDIR<=SIZE(FRNL,1)) THEN
       FRNL(IDIR,:)=FRNL(IDIR,:)+CFOR
    ENDIF

  END SUBROUTINE CONTRACT_SIGMA_G_TAU_NONL

!***********************************************************************
!
! FFT_G fourier transform along first index from G to R
!
!    from G%G_PROJ  ->  G%PROJ_R and G%R_PROJ are calculated by FFT
!    from G%GG      ->  G%G_R     is calculated by FFT
!
! the three new entires  G%PROJ_R, G%R_PROJ and  G%G_R are allocated
!
!***********************************************************************

  SUBROUTINE FFT_G( WDES, G, GDES, NK)
    USE ini
    TYPE (wavedes) :: WDES
    TYPE (greensf) :: G
    TYPE (greensfdes) :: GDES         
  ! local
    INTEGER :: ICOL
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1
    INTEGER :: NK

    PROFILING_START('fft_g')

    IF (ASSOCIATED(G%R_G).OR.ASSOCIATED(G%R_PROJ) .OR. ASSOCIATED(G%PROJ_R) ) THEN
       CALL vtutor%bug("internal error in FFT_G: G%R_G is associated " // str(ASSOCIATED(G%R_G)) // &
          " " // str(ASSOCIATED(G%R_PROJ)) // " " // str(ASSOCIATED(G%PROJ_R)), __FILE__, __LINE__)
    ENDIF

! so that we know how much storage is taken
    CALL ALLOCATE_R_G(GDES,G)

! Fock FFT is used
    CALL SETWDES(WDES, WDES1, NK)
    CALL NEWWAV(W1, WDES1,.TRUE.)

    DO ICOL=1, GDES%NRPLWV_COL_DATA_POINTS
       ! possibly complex to real FFT
       ! nasty CR is complex valued but only the real components
       ! matter in the Gamma-point only case
       CALL FFTWAV(WDES1%NGVECTOR, WDES1%NINDPW(1), &
            W1%CR(1), &
            G%GG(1,ICOL),WDES1%GRID)

!!       G%R_G(1:W1%WDES1%GRID%RL%NP, ICOL)=W1%CR(1:W1%WDES1%GRID%RL%NP)
!!       G%R_G(W1%WDES1%GRID%RL%NP+1:GDES%MPLWV_ROW, ICOL)=0  ! pad CR with zeros (just in case)
#ifndef gammareal
       CALL OMP_DCOPY(2*W1%WDES1%GRID%RL%NP,W1%CR(1),1,G%R_G(1,ICOL),1)
#else
       CALL OMP_DCOPY(  W1%WDES1%GRID%RL%NP,W1%CR(1),2,G%R_G(1,ICOL),1)
#endif
    ENDDO

    CALL ALLOCATE_G_R(GDES,G)
    CALL TRANSPOSE_R_G( G%R_G, G%G_R, GDES)

    ! we do not need R_G anymore, so destroy it
    CALL DEALLOCATE_R_G(G)


    CALL ALLOCATE_R_PROJ(GDES,G)
    DO ICOL=1,  GDES%NPRO_COL
       CALL FFTWAV(WDES1%NGVECTOR, WDES1%NINDPW(1), &
            W1%CR(1), G%G_PROJ(1,ICOL),WDES1%GRID)
       
!!       G%R_PROJ(1:W1%WDES1%GRID%RL%NP, ICOL)=W1%CR(1:W1%WDES1%GRID%RL%NP)
!!       G%R_PROJ(W1%WDES1%GRID%RL%NP+1:GDES%MPLWV_ROW, ICOL)=0  ! pad CR with zeros (just in case)
#ifndef gammareal
       CALL OMP_DCOPY(2*W1%WDES1%GRID%RL%NP,W1%CR(1),1,G%R_PROJ(1,ICOL),1)
#else
       CALL OMP_DCOPY(  W1%WDES1%GRID%RL%NP,W1%CR(1),2,G%R_PROJ(1,ICOL),1)
#endif
    ENDDO
       
    CALL ALLOCATE_PROJ_R(GDES,G)
    CALL TRANSPOSE_R_PROJ( G%R_PROJ, G%PROJ_R, GDES)
     
    CALL DELWAV(W1, .TRUE.)

    PROFILING_STOP('fft_g')

  END SUBROUTINE FFT_G

!
! deallocate intermediate functions
!   G%PROJ_R, G%R_PROJ and  G%G_R
!
  SUBROUTINE RELEASE_FFT_G( G)
    USE ini
    TYPE (greensf) :: G

    IF (ASSOCIATED(G%G_R)) THEN
       CALL DEALLOCATE_G_R(G)
    ENDIF
    IF (ASSOCIATED(G%PROJ_R)) THEN
       CALL DEALLOCATE_PROJ_R(G)
    ENDIF
    IF (ASSOCIATED(G%R_PROJ)) THEN
       CALL DEALLOCATE_R_PROJ(G)
    ENDIF

  END SUBROUTINE RELEASE_FFT_G


!***********************************************************************
!
! FFT_SIGMA fourier transform along first index from R to G
!
!    from SIGMA%G_R     ->  SIGMA%GG      is calculated by FFT
!
! SIGMA%G_R is deallocated
!
! scaling issue: to make this the exact reciprocal of FFT_G the 
! results need to be scaled by the number of grid points
! this might not be required to FFT the selfenergy however
! 
!***********************************************************************

  SUBROUTINE FFT_SIGMA( WDES, SIGMA, GDES, NK)
    USE ini
    TYPE (wavedes) :: WDES
    TYPE (greensf) :: SIGMA
    TYPE (greensfdes) :: GDES         
    INTEGER :: NK             ! k-point presently considered.
  ! local
    INTEGER :: ICOL
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1
    REAL(q) :: RINPL

    PROFILING_START('fft_sigma')

    IF (.NOT. ASSOCIATED(SIGMA%GG).OR. .NOT. ASSOCIATED(SIGMA%G_PROJ)) THEN
       CALL vtutor%bug("internal error in FFT_SIGMA: SIGMA%GG is not associated " // &
          str(ASSOCIATED(SIGMA%GG)) // " " // str(ASSOCIATED(SIGMA%G_PROJ)), __FILE__, __LINE__)
    ENDIF

! Fock FFT is used
    CALL SETWDES(WDES, WDES1, NK)
    CALL NEWWAV(W1, WDES1,.TRUE.)

    IF (WDES%LOVERL) THEN
    ENDIF
! now handle G_R: transpose G_R -> R_G, and FFT to obtain GG
    CALL ALLOCATE_R_G(GDES, SIGMA )

    CALL TRANSPOSE_G_R(SIGMA%G_R, SIGMA%R_G, GDES)
! deallocate the SIGMA%G_R entry
    CALL DEALLOCATE_G_R(SIGMA)

! allocate SIGMA%GG
    DO ICOL=1, GDES%NRPLWV_COL_DATA_POINTS
!!       W1%CR(1:W1%WDES1%GRID%RL%NP)=SIGMA%R_G(1:W1%WDES1%GRID%RL%NP, ICOL)*(1.0_q/W1%WDES1%GRID%NPLWV)
       RINPL=1.0_q/W1%WDES1%GRID%NPLWV
#ifndef gammareal
       CALL OMP_DSCPY(2*W1%WDES1%GRID%RL%NP,RINPL,SIGMA%R_G(1,ICOL),1,W1%CR(1),1)
#else
       CALL REAL_CMPLX_SCALE(W1%WDES1%GRID%RL%NP,SIGMA%R_G(1,ICOL),RINPL,W1%CR(1))
#endif
       ! possibly complex to real FFT
       ! nasty CR is complex valued but only the real components
       ! matter in the Gamma-point only case
       CALL FFTEXT(WDES1%NGVECTOR, WDES1%NINDPW(1), &
            W1%CR(1), &
            SIGMA%GG(1,ICOL),WDES1%GRID,.TRUE.)
    ENDDO

    CALL DEALLOCATE_R_G(SIGMA)

    CALL DELWAV(W1, .TRUE.)

    PROFILING_STOP('fft_sigma')

  END SUBROUTINE FFT_SIGMA

!***********************************************************************
!
! contract selfenergy and Green's function G
! this can be used to determine the correlation energy in second order
! since
!  G(-tau) G(tau) W(tau) = E_direct_MP2
!
! is W(tau) = v X(tau) v
! e.g. the second order expansion of the screened potential
!
!***********************************************************************

  SUBROUTINE CONTRACT_SIGMA_G_TAU( W, GU, SIGMA, GDES, NK2, IU6,&
    TAUWEIGHT, LDO_TAU_LOCAL, MP2 )
    USE ini
    USE dfast
    USE fock
    USE constant
    IMPLICIT NONE
    TYPE (wavespin):: W
    TYPE (greensf) :: GU         ! Green function for unoccupied orbitals 
    TYPE (greensf) :: SIGMA      ! Green function for occupied orbitals
    TYPE (greensfdes) :: GDES   
    INTEGER :: NK2               ! k-point at which GU and SIGMA is 
    INTEGER :: IU6
    REAL(q) :: TAUWEIGHT
    LOGICAL :: LDO_TAU_LOCAL
    COMPLEX(q) :: MP2
! local
    TYPE (wavedes1) :: WDES2
    COMPLEX(q)      :: CSUM_GG 
    COMPLEX(q)      :: CSUM_G_PROJ
    COMPLEX(q)      :: CSUM_PROJ_PROJ
    COMPLEX(q)      :: CSUM
    INTEGER :: ICOL, NG, NGP, NPRO

    IF ( .NOT. LDO_TAU_LOCAL ) RETURN

    PROFILING_START('contract_sigma_g_tau')

! set k-point descriptor
    CALL SETWDES(W%WDES,WDES2,NK2)

!   start with plane wave part
    CSUM_GG=0

    DO ICOL=1, GDES%NRPLWV_COL_DATA_POINTS
       ! determine corresponding plane wave coefficient (note sine cosine transform)
       NGP=ICOL-2+ m_ GDES%NRPLWV_POS
       ! make sure that we do not index beyond last relevant storage position
       ! see notes below for m_ 
       IF (NGP > m_ WDES2%NGVECTOR) CYCLE

       ! only use actual data points
       ! in the Gamma point only mode the number of entries in the real valued
       ! Greens function is doubled (sine and cosine transformation)
       ! m_ is define to be 2 for gammareal version taking care of this issue
       DO NG=1,m_ WDES2%NGVECTOR
          CSUM_GG=CSUM_GG+SIGMA%GG(NG,ICOL)*GCONJG( GU%GG(NG,ICOL))
       ENDDO
    ENDDO


! non local parts only for PAW and US-PP
    IF (W%WDES%LOVERL) THEN

!   first mixed part
    CSUM_G_PROJ=0
    DO ICOL=1, GDES%NPRO_COL
       ! only use actual data points
       DO NG=1,m_ WDES2%NGVECTOR
          CSUM_G_PROJ=CSUM_G_PROJ+SIGMA%G_PROJ(NG,ICOL)*GCONJG( GU%G_PROJ(NG,ICOL))
       ENDDO
    ENDDO

!   projector part
    CSUM_PROJ_PROJ=0
    IF (ASSOCIATED(SIGMA%PROJ_PROJ)) THEN
    DO ICOL=1, GDES%NPRO_COL
       ! only use actual data points
       DO NPRO=1,WDES2%NPRO ! same as GDES%NPRO_ROW
          CSUM_PROJ_PROJ=CSUM_PROJ_PROJ+SIGMA%PROJ_PROJ(NPRO,ICOL)*GCONJG( GU%PROJ_PROJ(NPRO,ICOL))
       ENDDO
    ENDDO
    ENDIF

    ENDIF

    CALLMPI( M_sum_z(GDES%COMM, CSUM_GG , 1 ))
    CALLMPI( M_sum_z(GDES%COMM, CSUM_G_PROJ , 1 ))
    CALLMPI( M_sum_z(GDES%COMM, CSUM_PROJ_PROJ , 1 ))

    CSUM=CSUM_GG+2*REAL(CSUM_G_PROJ,KIND=q)+  CSUM_PROJ_PROJ

    !MP2 shoud be 1/2 Tr(G W^(2) G)
    MP2=CSUM*TAUWEIGHT/2

#ifdef MPI
    IF (W%WDES%COMM%NODE_ME==1  ) THEN
#endif
!       WRITE(*,'(A,I4,4F20.10)') ' MP2  CONTRACT_SIGMA_G_TAU     ',NK2,REAL(MP2,q), & 
!       REAL(CSUM_GG,q)*TAUWEIGHT/2,REAL(CSUM_G_PROJ,q)*TAUWEIGHT,REAL(CSUM_PROJ_PROJ,q)*TAUWEIGHT/2
#ifdef MPI
    ENDIF
#endif
    PROFILING_STOP('contract_sigma_g_tau')

  END SUBROUTINE CONTRACT_SIGMA_G_TAU



!***********************************************************************
!
! calculate from the self-energy in reciprocal space Sigma(G,G',tau)
! the self-energy in the basis of the orbitals
!
!  <i |G> Sigma(G,G',tau) <G' | a > = Sigma(i,a, tau)
!
! on input 
!   SIGMA stores the self-energy in plane-wave basis for one tau point
!
! on output 
!  SIGMA_DIAG stores the diagonal elements of the self-energy in orbital basis
!  SIGMA_MAT(NROW*NCOL, NKPTS_IRZ, ISPIN, NOMEGA=1 )  
!       (if associated) stores the complete self-energy matrix in orbital
!              basis for currently treated k-point
!
!***********************************************************************

  SUBROUTINE CALCULATE_SIGMA_ORBITAL( W, SIGMA, SIGMA_DIAG, SIGMA_MAT, GDES_MAT, & 
       GDES, IMAG_GRIDS, IU6, NK1, ISP, LOCCUPIED, MP2SUM, TAUWEIGHT )
    USE dfast
    USE ini
    USE scala
    USE greens_orbital
    USE minimax_struct, ONLY: imag_grid_handle
    TYPE (wavespin)                :: W
    TYPE (greensf)                 :: SIGMA              ! Sigma(G,G',tau)
    TYPE (greensfdes)              :: GDES               ! descriptor for Green's function
    GDEF                           :: SIGMA_DIAG(:)      ! diagonal components of self-energy <i| Sigma | i>
    GDEF, POINTER                  :: SIGMA_MAT(:,:,:,:) ! self-energy matrix stored in distributed fashion (scaLAPACK)
                                                         ! (if not a null pointer)
    TYPE (greens_mat_des), POINTER :: GDES_MAT
    TYPE (imag_grid_handle)        :: IMAG_GRIDS         ! imaginary grids 
    INTEGER                        :: IU6                ! unit for IO (usually OUTCAR)
    INTEGER                        :: NK1                ! considered k-point
    INTEGER                        :: ISP                ! considered spin
    LOGICAL                        :: LOCCUPIED          ! calculate correlation energy (only implemented for Sigma at negative tau
    COMPLEX(q)                     :: MP2SUM             ! correlation energy using Galinski-Migdal
    REAL(q)                        :: TAUWEIGHT 
  ! local
    REAL(q)                        :: TAU                ! time point to be considered
    TYPE( communic)                :: COMM_INTRA         ! communicator in the time point for Greens(G,G',tau)
    TYPE( communic)                :: COMM_INTER    ! communicator between time points  for Greens(G,G',tau)
    LOGICAL                        :: LDO_TAU_LOCAL      ! maybe GEMM calls can be skipped
    INTEGER                        :: NB_FIRST, NB_LAST
    GDEF, ALLOCATABLE              :: CPROJ(:,:), CPROJP(:,:)
    COMPLEX(q), ALLOCATABLE        :: CG(:,:), CGP(:,:)
    COMPLEX(q), ALLOCATABLE        :: CG_COLUMN   (:,:)
    GDEF, ALLOCATABLE              :: CPROJ_COLUMN(:,:)
    TYPE (wavefun1)                :: W1
    TYPE (wavedes1)                :: WDES1
    INTEGER                        :: NB1, NB2, NSTRIP, NSTRIP_GLOBAL, N, NN, NGLB, NB1_INTO_TOT, NPRO, NG
    COMPLEX(q)                     :: CSUM_GG
    GDEF                           :: CSUM_PROJ
    COMPLEX(q)                     :: CORR
    GDEF,ALLOCATABLE               :: MAT(:,:)           ! matrix to temporarily store the self-energy matrix elements
    GDEF                           :: R(W%WDES%NB_TOT)   ! temporary array to accumulate the results from this k-point  
    INTEGER                        :: NTAU
    REAL(q)                        :: SWEIGHT

    PROFILING_START('calculate_sigma_orbital')

    TAU = IMAG_GRIDS%T%POINT_CURRENT
    COMM_INTRA = IMAG_GRIDS%T%COMM_IN_GROUP
    COMM_INTER = IMAG_GRIDS%T%COMM_BETWEEN_GROUPS
    LDO_TAU_LOCAL = IMAG_GRIDS%T%LDO_POINT_LOCAL 

    IF (GDES%NPRO_ROW.GT.0) THEN
       ! if SIGMA%PROJ_G is not yet associated calculate it (force routines have possible done this)
       IF (.NOT. ASSOCIATED(SIGMA%PROJ_G)) THEN
          CALL ALLOCATE_PROJ_G(GDES,SIGMA)
          CALL TRANSPOSE_G_PROJ(SIGMA%G_PROJ, SIGMA%PROJ_G, GDES)
       ENDIF
    ENDIF

    CORR=(0._q,0._q)

    NB_FIRST=1
    NB_LAST =W%WDES%NBANDS

    NSTRIP_GLOBAL=NSTRIP_STANDARD_GLOBAL
#ifdef MPI
    NSTRIP=(NSTRIP_GLOBAL+W%WDES%COMM_KIN%NCPU-1)/W%WDES%COMM_KIN%NCPU
    NSTRIP_GLOBAL=NSTRIP*W%WDES%COMM_KIN%NCPU
#else
    NSTRIP=NSTRIP_GLOBAL
#endif
    ALLOCATE(CG(GDES%NRPLWV_ROW,NSTRIP_GLOBAL),CPROJ(GDES%NPRO_ROW,NSTRIP_GLOBAL), &
             CGP(GDES%NRPLWV_ROW,NSTRIP_GLOBAL),CPROJP(GDES%NPRO_ROW,NSTRIP_GLOBAL), &
             CG_COLUMN(GDES%NRPLWV_COL,NSTRIP_GLOBAL),CPROJ_COLUMN(GDES%NPRO_COL,NSTRIP_GLOBAL))

    CPROJ=zero
    CG=(0._q,0._q)
    CPROJ_COLUMN=zero
    CG_COLUMN=(0._q,0._q)

    CALL SETWDES(W%WDES,WDES1,NK1)

    CALL NEWWAV(W1, WDES1, .TRUE.)
    R=zero


! use weight of first (Gamma) point
    SWEIGHT=SQRT(W%WDES%WTKPT(1)*W%WDES%RSPIN)
!=======================================================================
! loop over all bands in blocks of NSTRIP
!=======================================================================
    DO NB1=NB_FIRST, NB_LAST, NSTRIP
       NB2=MIN(NB_LAST, NB1+NSTRIP-1)
       ! gather all bands with local index NB1 up to and including NB2
       ! this requires communication between all cores, also between different tau 
       CALL W1_GATHER_ARRAY_RECIPROCAL(W, NB1, NB2, ISP, W1, CG, CPROJ)

       ! total number of bands  W1_GATHER_ARRAY has collected
       NGLB=(NB2-NB1+1)*W%WDES%NB_PAR

       ! copy relevant part over to CG_COLUMN
       CG_COLUMN(:,1:NGLB)=CG(GDES%NRPLWV_POS:GDES%NRPLWV_POS+GDES%NRPLWV_COL-1,1:NGLB)
       
       ! same for CPROJ_COLUMN
       IF (GDES%NPRO_COL/=0) THEN
          CALL GDIS_PROJ(GDES, CPROJ(:,1:NGLB), CPROJ_COLUMN(:,1:NGLB))
       ENDIF

       ! loop over all time/frequencies
       !skip GEMM calls if group has no time point
       IF ( LDO_TAU_LOCAL ) THEN
#ifdef gammareal
          ! construct all relevant Green's functions
          ! interpret coefficients as sine and cosine transforms
          ! we need all "pairs" as real final coefficients
          
          ! plane-wave  plane-wave part: a'(G) = \sum_G' Sigma(G,G', tau) a(G')
          ! (gamma only SIGMA%GG(1,1) real valued with  GDES%NRPLWV_ROW_DATA_POINTSx GDES%NRPLWV_COL_DATA_POINTS)
          CALL DGEMM('N','N', GDES%NRPLWV_ROW_DATA_POINTS, NGLB, GDES%NRPLWV_COL_DATA_POINTS, one, &
               SIGMA%GG(1,1), SIZE(SIGMA%GG, 1), CG_COLUMN(1,1), 2*SIZE(CG_COLUMN,1), &
               zero, CG(1,1), 2*SIZE(CG,1))
          
          ! add plane-wave  projector part: a'(G) = \sum_alpha' Sigma(G,alpha', tau) a(alpha')
          ! must be counted twice since the conjugated part is missing
          IF (GDES%NPRO_COL/=0) THEN
             CALL DGEMM('N','N', GDES%NRPLWV_ROW_DATA_POINTS, NGLB, GDES%NPRO_COL,  one, &
                  SIGMA%G_PROJ(1,1), SIZE(SIGMA%G_PROJ, 1), CPROJ_COLUMN(1,1), SIZE(CPROJ_COLUMN,1), &
                  one, CG(1,1), 2*SIZE(CG,1))
          ENDIF
          ! projector  projector part: a'(alpha) = \sum_alpha' Sigma(alpha,alpha', tau) a(alpha')
          IF (GDES%NPRO_COL/=0) THEN
             CALL DGEMM('N','N', GDES%NPRO_ROW, NGLB, GDES%NPRO_COL, one, &
                  SIGMA%PROJ_PROJ(1,1), SIZE(SIGMA%PROJ_PROJ, 1), CPROJ_COLUMN(1,1), SIZE(CPROJ_COLUMN,1), &
                  zero, CPROJ(1,1), SIZE(CPROJ,1))
          ELSE
             CPROJ=0
          ENDIF
          ! finally add projector plane-wave part
          IF (GDES%NPRO_ROW/=0) THEN
             CALL DGEMM('N','N', GDES%NPRO_ROW, NGLB, GDES%NRPLWV_COL_DATA_POINTS, one, &
                  SIGMA%PROJ_G(1,1), SIZE(SIGMA%PROJ_G, 1), CG_COLUMN(1,1), 2*SIZE(CG_COLUMN,1), &
                  one, CPROJ(1,1), SIZE(CPROJ,1) )
          ENDIF
#else
          ! plane-wave  plane-wave part: a'(G) = \sum_G' Sigma(G,G', tau) a(G')
          CALL ZGEMM('N','N', GDES%NRPLWV_ROW, NGLB, GDES%NRPLWV_COL, one, &
               SIGMA%GG(1,1), SIZE(SIGMA%GG, 1), CG_COLUMN(1,1), SIZE(CG_COLUMN,1), &
               zero,  CG(1,1), SIZE(CG,1))
          ! add plane-wave  projector part: a'(G) = \sum_alpha' Sigma(G,alpha', tau) a(alpha')
          IF (GDES%NPRO_COL/=0) THEN
             CALL ZGEMM('N','N', GDES%NRPLWV_ROW, NGLB, GDES%NPRO_COL, one, &
                  SIGMA%G_PROJ(1,1), SIZE(SIGMA%G_PROJ, 1), CPROJ_COLUMN(1,1), SIZE(CPROJ_COLUMN,1), &
                  one, CG(1,1), SIZE(CG,1) )
          ENDIF
          ! projector  projector part: a'(alpha) = \sum_alpha' Sigma(alpha,alpha', tau) a(alpha')
          IF (GDES%NPRO_COL/=0) THEN
             CALL GGEMM('N','N', GDES%NPRO_ROW, NGLB, GDES%NPRO_COL, one, &
                  SIGMA%PROJ_PROJ(1,1), SIZE(SIGMA%PROJ_PROJ, 1), CPROJ_COLUMN(1,1), SIZE(CPROJ_COLUMN,1), &
                  zero, CPROJ(1,1), SIZE(CPROJ,1))
          ELSE
             CPROJ=0
          ENDIF
          ! finally add projector plane-wave part
          IF (GDES%NPRO_ROW/=0) THEN
             CALL ZGEMM('N','N', GDES%NPRO_ROW, NGLB, GDES%NRPLWV_COL, one, &
                  SIGMA%PROJ_G(1,1), SIZE(SIGMA%PROJ_G, 1), CG_COLUMN(1,1), SIZE(CG_COLUMN,1), &
                  one, CPROJ(1,1), SIZE(CPROJ,1) )
          ENDIF
#endif
       ELSE
          CG=(0._q,0._q)
          CPROJ=zero
       ENDIF
       ! the summation was done over distributed index G'
       ! sum all local data over all columns G' (20 % of DGEMM)
       CALLMPI( M_sum_z(COMM_INTRA, CG, SIZE(CG)))
       CALLMPI( M_sum_g(COMM_INTRA, CPROJ, SIZE(CPROJ)))
       !=======================================================================
       ! we now have all resulting coefficients a'(G) and a'(alpha)  on ALL 
       ! cores sharing one time point
       ! contract this with the orbitals stored locally on each core (WDES%NBANDS)
       ! i.e. CG(a',G) is shared on all cores, a' is a stripe of orbitals currently 
       ! treated and G are all coefficients
       !=======================================================================
      DO NTAU=1, COMM_INTER%NCPU
       CGP=CG
       CPROJP=CPROJ
       ! broadcast to nodes that hold other time points
       CALLMPI( M_bcast_z_from(COMM_INTER, CGP, SIZE(CGP), NTAU))
       CALLMPI( M_bcast_g_from(COMM_INTER, CPROJP, SIZE(CPROJP), NTAU))
       !=======================================================================
       ! the following version calculates only the diagonal part 
       ! of the self-energy matrix
       !=======================================================================
       IF (.NOT. ASSOCIATED( SIGMA_MAT)) THEN
          ! loop over all a' = sigma(tau) a 
          DO N=1,NGLB
             ! for the time being only calculate "diagonal" element <n | Sigma(tau) | n>
             NB1_INTO_TOT=(NB1-1)*W%WDES%NB_PAR+N
             ! corresponding local band (if stored locally)
             NN=NB_LOCAL(NB1_INTO_TOT,WDES1)
             IF (NN>0) THEN
                ! is always complex 
                CSUM_GG=0
                DO NG=1,m_ WDES1%NGVECTOR
                   CSUM_GG=CSUM_GG+CONJG(W%CW(NG, NN, NK1, ISP))*CGP(NG, N)
                ENDDO
                
                CSUM_PROJ=zero
                IF (W%WDES%LOVERL) THEN
                   DO NPRO=1,WDES1%NPRO ! same as GDES%NPRO_ROW
                      CSUM_PROJ=CSUM_PROJ+GCONJG( W%CPROJ(NPRO, NN, NK1, ISP))*CPROJP(NPRO, N)
                   ENDDO
                ENDIF

                ! store self-energy here 
                R(NB1_INTO_TOT)=(CSUM_GG+CSUM_PROJ)*W%WDES%WTKPT(1)
                ! stupid G has a scaling factor SWEIGHT
                ! that should rather occur in contraction G G
                R(NB1_INTO_TOT)=R(NB1_INTO_TOT)*(1.0_q/SWEIGHT)
             ENDIF
          ENDDO
       !=======================================================================
       ! determine full matrix and store it in distributed array SIGMA_MAT
       !    <n'| Sigma(tau) |n>
       !=======================================================================
       ELSE
          
          ! parallelization strategy has been chosen, with the orbitals n' distributed
          ! and Sigma(tau) | n> stored on all cores
          ALLOCATE(MAT(WDES1%NB_TOT, NGLB))
          MAT=zero
          !multiply the current CGP(a',G,tau) with orbitals stored locally in W
          ! and store in MAT(all local orbitals, stripe size)
          CALL GGEMM(trans,'N',  WDES1%NBANDS, NGLB, m_ WDES1%NPL, one, &
               W%CW(1, 1, NK1, ISP), m_ WDES1%NRPLWV, & 
               CGP(1,1), m_ SIZE(CGP,1), zero, MAT(1,1), WDES1%NB_TOT)

          IF (WDES1%NPRO>0) THEN
             CALL GGEMM(trans,'N',  WDES1%NBANDS, NGLB, WDES1%NPRO, one, &
               W%CPROJ(1, 1, NK1, ISP), WDES1%NPROD, &
               CPROJP(1,1), SIZE(CPROJP,1), one, MAT(1,1), WDES1%NB_TOT)
          ENDIF

          ! stupid G has a scaling factor SWEIGHT
          ! that should rather occur in contraction G G
          MAT(1:WDES1%NBANDS,:)=MAT(1:WDES1%NBANDS,:)*(W%WDES%WTKPT(1)*(1.0_q/SWEIGHT))
          ! set remaining elements to zero
          MAT(WDES1%NBANDS+1:,:)=zero
          ! the orbitals are stored locally one after another, as they are distributed.
          ! to be able to use MPI_SUM we need to change the position
          ! of the band in the matrix MAT from local band number to global band number
          !
          ! now place calculated matrix elements in the proper position
          DO N=WDES1%NBANDS,1,-1
             ! global storage index of that band
             NN=WDES1%NB_LOW+(N-1)*WDES1%NB_PAR
             ! move data to proper position for seriel data layout
             MAT(NN,:)=MAT(N,:)
             ! clear source 
             IF (NN/=N) THEN
                MAT(N,:)=zero
             ENDIF
          ENDDO
          ! sum over all nodes of which bands are distributed
          CALLMPI( M_sum_g(W%WDES%COMM_KIN, MAT(1,1), WDES1%NB_TOT*NGLB))

          ! now distribute the results back over time points
          IF (NTAU==GDES_MAT%COMM_INTER%NODE_ME) THEN
             CALL DISTRI_SLICE_NOINT(MAT(1,1),SIZE(MAT,1), WDES1%NB_TOTK(ISP), SIGMA_MAT(1,1,1,1), GDES_MAT%DESC, & 
                  (NB1-1)*W%WDES%NB_PAR+1,NB2*W%WDES%NB_PAR)
          ELSE IF (NTAU>GDES_MAT%COMM_INTER%NCPU) THEN
             CALL vtutor%bug("internal error in VASP: NOMEGAPAR must be larger than NTAUPAR \n there " &
                // "are not enough arrays to store the result", __FILE__, __LINE__)
          ENDIF
          DEALLOCATE(MAT)
       ENDIF
      ENDDO
    ENDDO
    CALL DEALLOCATE_PROJ_G(SIGMA)

    IF (ASSOCIATED( SIGMA_MAT)) THEN
       CALL DETERMINE_DIAGONALE_GDEF(WDES1%NB_TOTK(ISP), SIGMA_MAT(:,1,1,1), R, GDES_MAT%DESC)
    ENDIF
    ! global sum to get all elements on all cores
    CALLMPI( M_sum_g(GDES_MAT%COMM_INTRA, R, WDES1%NB_TOTK(ISP) ))

    ! store in SIGMA_DIAG
    SIGMA_DIAG(1:WDES1%NB_TOTK(ISP))=R(1:WDES1%NB_TOTK(ISP))

    CORR=(0._q,0._q)
    IF ( LFINITE_TEMPERATURE ) THEN
       DO NB1_INTO_TOT=1,WDES1%NB_TOTK(ISP)
          ! for T>0 compute the negative GF if self-energy for positive times is
          ! passed to this routine and vice versa
          CORR = CORR + GF_WEIGHT( .NOT.LOCCUPIED, TAU, IMAG_GRIDS%T%BETA, W%EFERMI(ISP), &
             REAL(W%CELTOT(NB1_INTO_TOT,NK1,ISP),q), W%FERTOT(NB1_INTO_TOT,NK1,ISP) )* &
             R(NB1_INTO_TOT)
       ENDDO
    ELSE
       ! in case of positive time propagation
       IF (.NOT.LOCCUPIED) THEN
          DO NB1_INTO_TOT=1,WDES1%NB_TOTK(ISP)
             ! occupied orbital add to total correlation energy with proper weight
             ! note the weight comes in a second time because of double loop over k-points
             IF (REAL(W%CELTOT(NB1_INTO_TOT,NK1,ISP),q)-W%EFERMI(ISP)<0.0_q) THEN
                CORR=CORR+EXP( TAU*(REAL(W%CELTOT(NB1_INTO_TOT,NK1,ISP),q)-W%EFERMI(ISP)))&
                   *R(NB1_INTO_TOT)*W%FERTOT( NB1_INTO_TOT,NK1,ISP)
             ENDIF
          ENDDO
       ! in case of negative time propagation
       ELSE
          DO NB1_INTO_TOT=1,WDES1%NB_TOTK(ISP)
             ! unoccupied orbital add to total correlation energy with proper weight
             IF (REAL(W%CELTOT(NB1_INTO_TOT,NK1,ISP),q)-W%EFERMI(ISP)>=0.0_q) THEN
                CORR=CORR+EXP(-TAU*(REAL(W%CELTOT(NB1_INTO_TOT,NK1,ISP),q)-W%EFERMI(ISP)))&
                   *R(NB1_INTO_TOT)*(1-W%FERTOT(NB1_INTO_TOT,NK1,ISP))
             ENDIF
          ENDDO
       ENDIF
    ENDIF

    MP2SUM=MP2SUM+CORR*TAUWEIGHT*W%WDES%RSPIN/2
  
    IF (IU6>=0) THEN
       WRITE(*,'(L1,A,8F20.10)') LOCCUPIED,' energy of <n|G|n><n|sigma|n> ', CORR*TAUWEIGHT
    ENDIF

    ! TODO: routine is quite slow
    ! maybe connected with following fact: 
    ! observed NaN for pgi and gnu_mkl for jellium 
    IF ( MP2SUM /= MP2SUM ) THEN
       CALL vtutor%error("internal ERROR in VASP: CALCULATE_SIGMA_ORBITAL &
                         &produced NaN. Try using less MPI ranks for this job.")
    ENDIF


    CALL DELWAV(W1, .TRUE.)    
    DEALLOCATE(CG, CPROJ, CGP, CPROJP, CG_COLUMN, CPROJ_COLUMN)

    PROFILING_STOP('calculate_sigma_orbital')

  END SUBROUTINE CALCULATE_SIGMA_ORBITAL


!***********************************************************************
!
! calculate from the Green's function in orbital representation G(a,b',tau)
! the Green's function in reciprocal space and projector components
!
!  G(G,G',tau)  = <G |i>G(i,a, tau)<a|G'>
!
!
!***********************************************************************

  SUBROUTINE TRANSFORM_G_ORBIT_PW( W, G, G_MAT, GDES_MAT, & 
       GDES, IU6, NK1, NQ_IRZ, ISP, & 
       COMM_INTRA, COMM_INTER, LDO_TAU_LOCAL, NTAU_ROOT)
    USE dfast
    USE ini
    USE scala
    USE greens_orbital
    TYPE (wavespin):: W
    TYPE (greensf) :: G           ! G(G,G',tau)
    TYPE (greensfdes) :: GDES     ! descriptor for Green's function
    GDEF, POINTER :: G_MAT(:,:,:,:) ! G matrix stored in distributed fashion (scaLAPACK)
                                  ! (if not a null pointer)
    TYPE (greens_mat_des), POINTER :: GDES_MAT
    INTEGER :: IU6                ! unit for IO (usually OUTCAR)
    INTEGER :: NK1                ! considered k-point, index into G
    INTEGER :: NQ_IRZ             ! index into G_MAT
    INTEGER :: ISP                ! considered spin
    TYPE( communic) :: COMM_INTRA ! communicator in the time point for Greens(G,G',tau)
    TYPE( communic) :: COMM_INTER ! communicator between time points  for Greens(G,G',tau)
    LOGICAL :: LDO_TAU_LOCAL      ! maybe GEMM calls can be skipped
    INTEGER :: NTAU_ROOT
  ! local
    INTEGER :: NB_FIRST, NB_LAST
    GDEF, ALLOCATABLE       :: CPROJ(:,:), CPROJP(:,:)
    COMPLEX(q), ALLOCATABLE :: CG(:,:), CGP(:,:)
    COMPLEX(q), ALLOCATABLE :: CG_COLUMN   (:,:)
    GDEF, ALLOCATABLE       :: CPROJ_COLUMN(:,:)
    TYPE (wavefun1):: W1
    TYPE (wavedes1) :: WDES1
    INTEGER :: NB1, NB2, NSTRIP, NSTRIP_GLOBAL, N, NN, NGLB, NB1_INTO_TOT, NPRO, NG, NB1_TOT, NB2_TOT
    GDEF, ALLOCATABLE :: MAT(:,:)       ! matrix to temporarily store the self-energy matrix elements
    GDEF, ALLOCATABLE :: MATS(:,:)      ! matrix to temporarily store the self-energy matrix elements
    COMPLEX(q) :: R(W%WDES%NB_TOT)      ! temporary array to accumulate the results from this k-point  
    INTEGER :: NTAU
    REAL(q) :: SWEIGHT

    PROFILING_START('transform_g_orbit_pw')

! use weight of first (Gamma) point
    SWEIGHT=SQRT(W%WDES%WTKPT(1)*W%WDES%RSPIN)

    IF (W%WDES%COMM%NODE_ME.eq.0) THEN
!    WRITE (*,*) 'g 1'
!    WRITE (*,'(10f12.8)') G%GG(1:10,1:10)
!    WRITE (*,*) 'g g proj', SIZE(G%G_PROJ,1), SIZE(G%G_PROJ,2)
!    WRITE (*,'(8f12.8)') G%G_PROJ(1:8,1:8)
!    WRITE (*,*) 'g proj proj', SIZE(G%PROJ_PROJ,1), SIZE(G%PROJ_PROJ,2)
!    WRITE (*,'(8f12.8)') G%PROJ_PROJ(1:8,1:8)
    !WRITE (*,'(10f12.8)') SUM(G%GG,2)
    ENDIF
 
    !!!!!!!!!!!!!
    G%GG=0._q
    G%G_PROJ=0._q
    G%PROJ_PROJ=0._q
    !!!!!!!!!!!!! 
    

    NB_FIRST=1
    NB_LAST =W%WDES%NBANDS

    NSTRIP_GLOBAL=NSTRIP_STANDARD_GLOBAL
#ifdef MPI
    NSTRIP=(NSTRIP_GLOBAL+W%WDES%COMM_KIN%NCPU-1)/W%WDES%COMM_KIN%NCPU
    NSTRIP_GLOBAL=NSTRIP*W%WDES%COMM_KIN%NCPU
#else
    NSTRIP=NSTRIP_GLOBAL
#endif
    ALLOCATE(CG(GDES%NRPLWV_ROW,NSTRIP_GLOBAL),CPROJ(GDES%NPRO_ROW,NSTRIP_GLOBAL), &
             CGP(GDES%NRPLWV_ROW,NSTRIP_GLOBAL),CPROJP(GDES%NPRO_ROW,NSTRIP_GLOBAL), &
             CG_COLUMN(GDES%NRPLWV_COL,NSTRIP_GLOBAL),CPROJ_COLUMN(GDES%NPRO_COL,NSTRIP_GLOBAL))

    CPROJ=0.0_q
    CG=0.0_q
    CPROJ_COLUMN=0.0_q
    CG_COLUMN=0.0_q

    CALL SETWDES(W%WDES,WDES1,NK1)

    CALL NEWWAV(W1, WDES1, .TRUE.)

    !first step:
    ! transform G(a,b',tau) to G(G,b',tau)
    !WRITE (*,*) 'node ', W%WDES%COMM%NODE_ME, GDES_MAT%COMM_INTER%NODE_ME,  COMM_INTER%NCPU,W%WDES%COMM_KIN%NODE_ME, GDES_MAT%COMM_INTRA%NODE_ME, COMM_INTRA%NODE_ME
    !loop over stripes of b' (all orbitals in total)
    DO NB1=NB_FIRST, NB_LAST, NSTRIP
       NB2=MIN(NB_LAST, NB1+NSTRIP-1)  !?? needs to be properly checked
       NB1_TOT=(NB1-1)*W%WDES%NB_PAR+1
       NB2_TOT=NB2*W%WDES%NB_PAR
       NGLB=NB2_TOT-NB1_TOT+1
       ALLOCATE(MAT(WDES1%NB_TOT, NGLB))
       ALLOCATE(MATS(WDES1%NBANDS, NGLB))

       !loop over all tau points currently treated
       DO NTAU=1, COMM_INTER%NCPU
          MAT=0._q  !just to be sure
          MATS=0._q !can be removed
          
          !collect G(a, b':b'+N, tau) and send out to all nodes
          IF (NTAU==GDES_MAT%COMM_INTER%NODE_ME) THEN
             !fills the array MAT(a, b':b'+N) with locally available data 
             !WRITE (*,*) 'getting slice ', NTAU, (NB1-1)*W%WDES%NB_PAR+1, NB2*W%WDES%NB_PAR, WDES1%NB_TOT, NGLB
             CALL RECON_SLICE(MAT(1,1), SIZE(MAT,1), WDES1%NB_TOTK(ISP), G_MAT(:,NQ_IRZ,ISP,NTAU_ROOT), &
                   GDES_MAT%DESC, (NB1-1)*W%WDES%NB_PAR+1,NB2*W%WDES%NB_PAR)
             CALLMPI( M_sum_g(GDES_MAT%COMM_INTRA, MAT(1,1), WDES1%NB_TOT*NGLB)) 
          ENDIF

          !broadcast to nodes that hold other tau points
          CALLMPI(M_bcast_g_from(GDES_MAT%COMM_INTER, MAT, SIZE(MAT), NTAU))

       IF (W%WDES%COMM%NODE_ME.eq.0) THEN
       !this looks correct as well
!       WRITE (101,*) 'MAT in TRANSFORM_G_GW'
!       !WRITE (*,'(10f12.8)') MAT(1:10,1:10)
!       WRITE (101,'(16f12.8)') REAL(MAT(:,:))
#ifndef gammareal
!       WRITE (101,*) ' '
!       WRITE (101,'(16f12.8)') IMAG(MAT(:,:))
#endif
       ENDIF

          !change from global to local idx in MAT (also could reduce the size of the array)
          DO N=1,WDES1%NBANDS
             !global storage index 
             NN=WDES1%NB_LOW+(N-1)*WDES1%NB_PAR
             !move data to local storage idx
             MATS(N,:)=MAT(NN,:) *SWEIGHT
          ENDDO
          !WRITE (*,*) 'MAT: ', SUM(MATS), SUM(MAT)
          !MATS=0.0
          !DO N=1,4
          !   MATS(N,N)=1.0
          !ENDDO
          !clear the other bands
          !MAT(WDES1%NBANDS+1:,:)=0

          !multiply with the stored W coefficients to get G(G,b',tau)
          ! we want basically an opposite to
          !MAT(NB_TOT,NGLB)=CGP(PW_row,NGLB)*W%CW(PW,NGLB)
          ! so CGP(PW_row,NGLB)=W%CW(PW,NBANDS)*MATS(NBANDS,NGLB)
          ! m_ WDES1%NPL    : current number of PWs
          ! m_ WDES1%NRPLWV : size of the wvfc array 
          CGP=0._q
          CPROJP=0._q

          !IF (W%WDES%COMM%NODE_ME.eq.0) THEN
          !   WRITE(*,*) 'multiplying W%CW * MATS = CGP '
          !   WRITE (*,*) GDES%NRPLWV_ROW, NGLB, WDES1%NBANDS
          !   WRITE (*,*) SIZE(W%CW,1), SIZE(W%CW,2), GDES%NRPLWV_ROW, WDES1%NBANDS, m_ WDES1%NRPLWV
          !   WRITE (*,*) SIZE(MATS,1), SIZE(MATS,2), NGLB, WDES1%NBANDS, WDES1%NBANDS
          !   WRITE (*,*) SIZE(CGP,1), SIZE(CGP,2), GDES%NRPLWV_ROW, NGLB, m_ SIZE(CGP,1)
          !   WRITE (*,*) 'CW ', W%WDES%COMM%NODE_ME, SUM(W%CW(:, 1: WDES1%NBANDS, NK1, ISP))
          !   WRITE (*,*) 'MATS ', W%WDES%COMM%NODE_ME, SUM(MATS)
          !ENDIF
          CALL GGEMM('N','N', m_ GDES%NRPLWV_ROW, NGLB, WDES1%NBANDS, one, &
                 W%CW(1, 1, NK1, ISP), m_ WDES1%NRPLWV, &
                 MATS(1,1), WDES1%NBANDS, zero, &
                 CGP(1,1), m_ SIZE(CGP,1))

          IF (WDES1%NPRO>0) THEN
             IF (WDES1%LOVERL) THEN
                ! W%CPROJ (NPROD, NBANDS,NKPTS,ISPIN),  NPROD =WDES%NPROD, NBANDS=WDES%NBANDS
                ! projections for locally stored bands, so OK
                CALL GGEMM('N','N', GDES%NPRO_ROW, NGLB, WDES1%NBANDS, one, &
                       W%CPROJ(1, 1, NK1, ISP), SIZE(W%CPROJ(:, :, NK1, ISP),1), &
                       MATS(1,1), WDES1%NBANDS, zero, &
                       CPROJP(1,1), SIZE(CPROJP,1))
                !not sure about the size of W%CPROJ again, so use the usual trick
             ELSE 
                CPROJP=0._q
             ENDIF
          ENDIF
          IF (W%WDES%COMM%NODE_ME.eq.0) THEN
!             WRITE (*,*) 'CGP 1'
!             WRITE (*,'(8f12.8)') CGP(1:8,1:8)
!            !WRITE (*,'(10f12.8)') SUM(G%GG,2)
!             WRITE (*,*) 'sum ', W%WDES%COMM%NODE_ME, SUM(CGP)
          ENDIF

          !locally the CGP array contains only contribution to G(G,b',tau) from local bands
          !need to be globally summed to obtain the complete matrix for all a bands
          ! sum over all nodes of which bands are distributed
          CALLMPI( M_sum_z(W%WDES%COMM_KIN, CGP(1,1), SIZE(CGP,1)*NGLB))
          CALLMPI( M_sum_g(W%WDES%COMM_KIN, CPROJP(1,1), SIZE(CPROJP,1)*NGLB))
          !the complete G(G,b',tau) for one tau on all cores 
          IF (W%WDES%COMM%NODE_ME.eq.0) THEN
!             WRITE (*,*) 'CGP 2'
!             WRITE (*,'(8f12.8)') CGP(1:8,1:8)
!             !WRITE (*,'(10f12.8)') SUM(G%GG,2)
          ENDIF
    
          !store on appropriate node with the proper tau point
          IF (NTAU==GDES_MAT%COMM_INTER%NODE_ME) THEN
             !store the complete matrix locally
             CG=CGP
             CPROJ=CPROJP
          ENDIF
       ENDDO
       !end tau loop



       CGP=0._q     ! not required if everything was done properly
       CPROJP=0._q  ! or so hopefully
       !collect stripe of bands b'(G') on all cores
       CALL W1_GATHER_ARRAY_RECIPROCAL(W, NB1, NB2, ISP, W1, CGP, CPROJP)
       ! copy relevant part over to CG_COLUMN
       CG_COLUMN=0._q
       CG_COLUMN(:,1:NGLB)=CGP(GDES%NRPLWV_POS:GDES%NRPLWV_POS+GDES%NRPLWV_COL-1,1:NGLB)

       ! same for CPROJ_COLUMN
       IF (GDES%NPRO_COL/=0) THEN
          CALL GDIS_PROJ(GDES, CPROJP(:,1:NGLB), CPROJ_COLUMN(:,1:NGLB))
       ENDIF

       !add G(G,b'->b'+N,tau)<b'|G'> to local G(G,G',tau)
       !we want opposite to b'(G) = \sum_G' Sigma(G,G', tau) b(G')   
       ! i.e.  G(G,G',tau)=\sum b' G(G,b',tau) b'(G')
       !  G%GG(GDES%NRPLWV_ROW,GDES%NRPLWV_COL)= \sum NGLB  
       !                CG(GDES%NRPLWV_ROW,NGLB) CG_COLUMN(GDES%NRPLWV_COL ,NGLB)
       ! TODO most likely wrong if NGLB .ne. NSTRIP_GLOBAL, CG_COLUMN ( , NSTRIP_GLOBAL)
       IF (W%WDES%COMM%NODE_ME.eq.0) THEN
!       WRITE (*,*) 'multiplying CG*CG_COLUMN=G%GG'
!       WRITE (*,*) GDES%NRPLWV_ROW, GDES%NRPLWV_COL, NGLB
!       WRITE (*,*) SIZE(CG,1), SIZE(CG,2), GDES%NRPLWV_ROW,  NGLB, SIZE(CG,1)
!       WRITE (*,*) SIZE(CG_COLUMN,1), SIZE(CG_COLUMN,2), GDES%NRPLWV_COL, NGLB, SIZE(CG_COLUMN,1)
!       WRITE (*,*) SIZE(G%GG,1), SIZE(G%GG,2), GDES%NRPLWV_ROW, GDES%NRPLWV_COL, SIZE(G%GG, 1)
       ENDIF

#ifdef gammareal
       CALL DGEMM('N','T', GDES%NRPLWV_ROW_DATA_POINTS, GDES%NRPLWV_COL_DATA_POINTS, NGLB, 1.0_q, &
            CG(1,1), 2*SIZE(CG,1), CG_COLUMN(1,1), 2*SIZE(CG_COLUMN,1), &
            1.0_q, G%GG(1,1), SIZE(G%GG, 1))
#else
       CALL ZGEMM('N', 'C', GDES%NRPLWV_ROW, GDES%NRPLWV_COL, NGLB, one, &
            CG(1,1), SIZE(CG,1), CG_COLUMN(1,1), SIZE(CG_COLUMN,1), &
            one, G%GG(1,1), SIZE(G%GG, 1))
#endif

       IF (GDES%NPRO_COL/=0) THEN 
          IF (WDES1%LOVERL) THEN
#ifdef gammareal
          CALL DGEMM('N','T', GDES%NRPLWV_ROW_DATA_POINTS, GDES%NPRO_COL, NGLB, 1.0_q, &
               CG(1,1), 2*SIZE(CG,1), CPROJ_COLUMN(1,1), SIZE(CPROJ_COLUMN,1), &
               1.0_q, G%G_PROJ(1,1), SIZE(G%G_PROJ, 1))
#else
          CALL ZGEMM('N', 'C', GDES%NRPLWV_ROW, GDES%NPRO_COL, NGLB, one, &
               CG(1,1), SIZE(CG,1), CPROJ_COLUMN(1,1), SIZE(CPROJ_COLUMN,1), &
               one, G%G_PROJ(1,1), SIZE(G%G_PROJ, 1))
#endif
          CALL GGEMM('N', trans, GDES%NPRO_ROW, GDES%NPRO_COL, NGLB, one, &
               CPROJ(1,1), SIZE(CPROJ,1), CPROJ_COLUMN(1,1), SIZE(CPROJ_COLUMN,1), &
               one, G%PROJ_PROJ(1,1), SIZE(G%PROJ_PROJ, 1))
          ELSE
             G%G_PROJ=0._q
             G%PROJ_PROJ=0._q
          ENDIF
       ENDIF
       DEALLOCATE(MAT,MATS)
    ENDDO
!    WRITE(*,*) 'leaving the subroutine'
    !end NB1 loop
    IF (W%WDES%COMM%NODE_ME.eq.0) THEN
!    WRITE (*,*) 'G 2'
!    WRITE (*,'(10f12.8)') G%GG(1:10,1:10)
!    WRITE (*,*) 'G g proj', SIZE(G%G_PROJ,1), SIZE(G%G_PROJ,2)
!    WRITE (*,'(8f12.8)') G%G_PROJ(1:8,1:8)
!    WRITE (*,*) 'G proj proj', SIZE(G%PROJ_PROJ,1), SIZE(G%PROJ_PROJ,2)
!    WRITE (*,'(8f12.8)') G%PROJ_PROJ(1:8,1:8)
!    !WRITE (*,'(10f12.8)') SUM(G%GG,2)
    ENDIF

    CALL DELWAV(W1, .TRUE.)
    DEALLOCATE(CG, CPROJ, CGP, CPROJP, CG_COLUMN, CPROJ_COLUMN)

    PROFILING_STOP('transform_g_orbit_pw')

    END SUBROUTINE TRANSFORM_G_ORBIT_PW


!***********************************************************************
!
! individual allocation routines for the Green's function parts 
!
!***********************************************************************

  SUBROUTINE ALLOCATE_R_PROJ(GDES,G)
    USE ini
    TYPE (greensfdes):: GDES
    TYPE (greensf) :: G
    INTEGER        :: ISTAT 

    ALLOCATE(G%R_PROJ(GDES%MPLWV_ROW, GDES%NPRO_COL),STAT=ISTAT )
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "ALLOCATE_R_PROJ is not able to allocate "//&
         str( ndata * 8._q*GDES%MPLWV_ROW*GDES%NPRO_COL/1024) //&
         " kB of data on MPI rank 0." )
    ENDIF
    CALL REGISTER_ALLOCATE(wsgf*SIZE(G%R_PROJ,KIND=qi8), "Greensfun")

  END SUBROUTINE ALLOCATE_R_PROJ

  SUBROUTINE ALLOCATE_R_G(GDES, G, NK)
    USE ini
    TYPE (greensfdes):: GDES
    TYPE (greensf) :: G
    INTEGER, OPTIONAL:: NK
  ! local
    INTEGER :: NRPLWV_COL_DATA_POINTS
    INTEGER :: ISTAT 

    NRPLWV_COL_DATA_POINTS=GDES%NRPLWV_COL_DATA_POINTS
    IF (PRESENT(NK) .AND. ASSOCIATED(GDES%LUSEINV)) THEN
       IF (GDES%LUSEINV(NK)) THEN
          NRPLWV_COL_DATA_POINTS=GDES%NRPLWV_COL_DATA_POINTS_NK(NK)
       ENDIF
    ENDIF

    ALLOCATE(G%R_G(GDES%MPLWV_ROW, NRPLWV_COL_DATA_POINTS),STAT=ISTAT)
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "ALLOCATE_R_G is not able to allocate "//&
         str( ndata * 8._q*GDES%MPLWV_ROW*NRPLWV_COL_DATA_POINTS/1024) //&
         " kB of data on MPI rank 0." )
    ENDIF
    CALL REGISTER_ALLOCATE(wsgf*SIZE(G%R_G,KIND=qi8), "Greensfun")

  END SUBROUTINE ALLOCATE_R_G

  SUBROUTINE ALLOCATE_G_R(GDES, G, NK)
    USE ini
    TYPE (greensfdes):: GDES
    TYPE (greensf)   :: G
    INTEGER, OPTIONAL:: NK
  ! local
    INTEGER :: NRPLWV_ROW_DATA_POINTS
    INTEGER :: ISTAT 

    NRPLWV_ROW_DATA_POINTS=GDES%NRPLWV_ROW_DATA_POINTS
    IF (PRESENT(NK) .AND. ASSOCIATED(GDES%LUSEINV)) THEN
       IF (GDES%LUSEINV(NK)) THEN
          NRPLWV_ROW_DATA_POINTS=GDES%NRPLWV_ROW_DATA_POINTS_NK(NK)
       ENDIF
    ENDIF

    ALLOCATE(G%G_R(NRPLWV_ROW_DATA_POINTS, GDES%MPLWV_COL), STAT=ISTAT)
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "ALLOCATE_G_R is not able to allocate "//&
         str( ndata * 8._q*NRPLWV_ROW_DATA_POINTS*GDES%MPLWV_COL/1024) //&
         " kB of data on MPI rank 0." )
    ENDIF
    CALL REGISTER_ALLOCATE(wsgf*SIZE(G%G_R,KIND=qi8), "Greensfun")

  END SUBROUTINE ALLOCATE_G_R

  SUBROUTINE ALLOCATE_PROJ_R(GDES,G)
    USE ini
    TYPE (greensfdes):: GDES
    TYPE (greensf) :: G
    INTEGER :: ISTAT 

    ALLOCATE(G%PROJ_R(GDES%NPRO_ROW,  GDES%MPLWV_COL),STAT=ISTAT)
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "ALLOCATE_PROJ_R is not able to allocate "//&
         str( ndata * 8._q*GDES%NPRO_ROW*GDES%MPLWV_COL/1024) //&
         " kB of data on MPI rank 0." )
    ENDIF
    CALL REGISTER_ALLOCATE(wsgf*SIZE(G%PROJ_R,KIND=qi8), "Greensfun")

  END SUBROUTINE ALLOCATE_PROJ_R

  SUBROUTINE ALLOCATE_GG(GDES,G)
    USE ini
    TYPE (greensfdes):: GDES
    TYPE (greensf) :: G
    INTEGER :: ISTAT 
  
    ALLOCATE(G%GG(GDES%NRPLWV_ROW_DATA_POINTS, GDES%NRPLWV_COL_DATA_POINTS),STAT=ISTAT)
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "ALLOCATE_GG is not able to allocate "//&
         str( ndata * 8._q*GDES%NRPLWV_ROW_DATA_POINTS*GDES%NRPLWV_COL_DATA_POINTS/1024) //&
         " kB of data on MPI rank 0." )
    ENDIF
    CALL REGISTER_ALLOCATE(wsgf*SIZE(G%GG,KIND=qi8), "Greensfun")

  END SUBROUTINE ALLOCATE_GG

  SUBROUTINE ALLOCATE_G_PROJ(GDES,G)
    USE ini
    TYPE (greensfdes):: GDES
    TYPE (greensf) :: G
    INTEGER :: ISTAT 

    ALLOCATE(G%G_PROJ(GDES%NRPLWV_ROW_DATA_POINTS, GDES%NPRO_COL  ),STAT=ISTAT)
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "ALLOCATE_G_PROJ is not able to allocate "//&
         str( ndata * 8._q*GDES%NRPLWV_ROW_DATA_POINTS*GDES%NPRO_COL/1024) //&
         " kB of data on MPI rank 0." )
    ENDIF
    CALL REGISTER_ALLOCATE(wsgf*SIZE(G%G_PROJ,KIND=qi8), "Greensfun")

  END SUBROUTINE ALLOCATE_G_PROJ

  SUBROUTINE ALLOCATE_PROJ_G(GDES,G)
    USE ini
    TYPE (greensfdes):: GDES
    TYPE (greensf) :: G
    INTEGER :: ISTAT 

    ALLOCATE(G%PROJ_G(GDES%NPRO_ROW, GDES%NRPLWV_COL_DATA_POINTS ),STAT=ISTAT)
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "ALLOCATE_PROJ_G is not able to allocate "//&
         str( ndata * 8._q*GDES%NPRO_ROW*GDES%NRPLWV_COL_DATA_POINTS/1024) //&
         " kB of data on MPI rank 0." )
    ENDIF
    CALL REGISTER_ALLOCATE(wsgf*SIZE(G%PROJ_G,KIND=qi8), "Greensfun")

  END SUBROUTINE ALLOCATE_PROJ_G

  SUBROUTINE ALLOCATE_PROJ_PROJ(GDES,G)
    USE ini
    TYPE (greensfdes):: GDES
    TYPE (greensf) :: G
    INTEGER :: ISTAT 

    ALLOCATE(G%PROJ_PROJ(GDES%NPRO_ROW,GDES%NPRO_COL),STAT=ISTAT)
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "ALLOCATE_PROJ_PROJ is not able to allocate "//&
         str( ndata * 8._q*GDES%NPRO_ROW*GDES%NPRO_COL/1024) //&
         " kB of data on MPI rank 0." )
    ENDIF
    CALL REGISTER_ALLOCATE(wsgf*SIZE(G%PROJ_PROJ,KIND=qi8), "Greensfun")

  END SUBROUTINE ALLOCATE_PROJ_PROJ


!***********************************************************************
!
! individual deallocation routines for the Green's function parts 
!
!***********************************************************************

  SUBROUTINE DEALLOCATE_G(G)
    USE ini
    TYPE (greensf) :: G(:)
  ! local
    INTEGER ::  IK

    DO IK=1, SIZE(G,1)
       IF (ASSOCIATED(G(IK)%RR)) THEN
          CALL DEREGISTER_ALLOCATE(wsgf*SIZE(G(IK)%RR,KIND=qi8), "Greensfun")
          DEALLOCATE(G(IK)%RR)
          NULLIFY(G(IK)%RR)
       ENDIF
       IF (ASSOCIATED(G(IK)%R_PROJ)) THEN
          CALL DEREGISTER_ALLOCATE(wsgf*SIZE(G(IK)%R_PROJ,KIND=qi8), "Greensfun")
          DEALLOCATE(G(IK)%R_PROJ)
          NULLIFY(G(IK)%R_PROJ)
       ENDIF
       IF (ASSOCIATED(G(IK)%PROJ_R)) THEN
          CALL DEREGISTER_ALLOCATE(wsgf*SIZE(G(IK)%PROJ_R,KIND=qi8), "Greensfun")
          DEALLOCATE(G(IK)%PROJ_R)
          NULLIFY(G(IK)%PROJ_R)
       ENDIF
       IF (ASSOCIATED(G(IK)%GG)) THEN
          CALL DEREGISTER_ALLOCATE(wsgf*SIZE(G(IK)%GG,KIND=qi8), "Greensfun")
          DEALLOCATE(G(IK)%GG)
          NULLIFY(G(IK)%GG)
       ENDIF
       IF (ASSOCIATED(G(IK)%G_PROJ)) THEN
          CALL DEREGISTER_ALLOCATE(wsgf*SIZE(G(IK)%G_PROJ,KIND=qi8), "Greensfun")
          DEALLOCATE(G(IK)%G_PROJ)
          NULLIFY(G(IK)%G_PROJ)
       ENDIF
       IF (ASSOCIATED(G(IK)%G_R)) THEN
          CALL DEREGISTER_ALLOCATE(wsgf*SIZE(G(IK)%G_R,KIND=qi8), "Greensfun")
          DEALLOCATE(G(IK)%G_R)
          NULLIFY(G(IK)%G_R)
       ENDIF
       IF (ASSOCIATED(G(IK)%R_G)) THEN
          CALL DEREGISTER_ALLOCATE(wsgf*SIZE(G(IK)%R_G,KIND=qi8), "Greensfun")
          DEALLOCATE(G(IK)%R_G)
          NULLIFY(G(IK)%R_G)
       ENDIF
       IF (ASSOCIATED(G(IK)%PROJ_PROJ)) THEN
          CALL DEREGISTER_ALLOCATE(wsgf*SIZE(G(IK)%PROJ_PROJ,KIND=qi8), "Greensfun")
          DEALLOCATE(G(IK)%PROJ_PROJ)
          NULLIFY(G(IK)%PROJ_PROJ)
       ENDIF
       IF (ASSOCIATED(G(IK)%PROJ_G)) THEN
          CALL DEREGISTER_ALLOCATE(wsgf*SIZE(G(IK)%PROJ_G,KIND=qi8), "Greensfun")
          DEALLOCATE(G(IK)%PROJ_G)
          NULLIFY(G(IK)%PROJ_G)
       ENDIF
    ENDDO

  END SUBROUTINE DEALLOCATE_G

!
! individual deallocation routines
!

  SUBROUTINE DEALLOCATE_R_PROJ(G)
    USE ini
    TYPE (greensf) :: G

    IF (ASSOCIATED(G%R_PROJ)) THEN
       CALL DEREGISTER_ALLOCATE(wsgf*SIZE(G%R_PROJ,KIND=qi8), "Greensfun")
       DEALLOCATE(G%R_PROJ)
       NULLIFY(G%R_PROJ)
    ENDIF
  END SUBROUTINE DEALLOCATE_R_PROJ

  SUBROUTINE DEALLOCATE_R_G(G)
    USE ini
    TYPE (greensf) :: G

    IF (ASSOCIATED(G%R_G)) THEN
       CALL DEREGISTER_ALLOCATE(wsgf*SIZE(G%R_G,KIND=qi8), "Greensfun")
       DEALLOCATE(G%R_G)
       NULLIFY(G%R_G)
    ENDIF
  END SUBROUTINE DEALLOCATE_R_G

  SUBROUTINE DEALLOCATE_G_R(G)
    USE ini
    TYPE (greensf) :: G

    IF (ASSOCIATED(G%G_R)) THEN
       CALL DEREGISTER_ALLOCATE(wsgf*SIZE(G%G_R,KIND=qi8), "Greensfun")
       DEALLOCATE(G%G_R)
       NULLIFY(G%G_R)
    ENDIF
  END SUBROUTINE DEALLOCATE_G_R

  SUBROUTINE DEALLOCATE_PROJ_R(G)
    USE ini
    TYPE (greensf) :: G

    IF (ASSOCIATED(G%PROJ_R)) THEN
       CALL DEREGISTER_ALLOCATE(wsgf*SIZE(G%PROJ_R,KIND=qi8), "Greensfun")
       DEALLOCATE(G%PROJ_R)
       NULLIFY(G%PROJ_R)
    ENDIF
  END SUBROUTINE DEALLOCATE_PROJ_R

  SUBROUTINE DEALLOCATE_GG(G)
    USE ini
    TYPE (greensf) :: G
  
    IF (ASSOCIATED(G%GG)) THEN
       CALL DEREGISTER_ALLOCATE(wsgf*SIZE(G%GG,KIND=qi8), "Greensfun")
       DEALLOCATE(G%GG)
       NULLIFY(G%GG)
    ENDIF
  END SUBROUTINE DEALLOCATE_GG

  SUBROUTINE DEALLOCATE_G_PROJ(G)
    USE ini
    TYPE (greensf) :: G

    IF (ASSOCIATED(G%G_PROJ)) THEN
       CALL DEREGISTER_ALLOCATE(wsgf*SIZE(G%G_PROJ,KIND=qi8), "Greensfun")
       DEALLOCATE(G%G_PROJ)
       NULLIFY(G%G_PROJ)
    ENDIF
  END SUBROUTINE DEALLOCATE_G_PROJ

  SUBROUTINE DEALLOCATE_PROJ_G(G)
    USE ini
    TYPE (greensf) :: G

    IF (ASSOCIATED(G%PROJ_G)) THEN
       CALL DEREGISTER_ALLOCATE(wsgf*SIZE(G%PROJ_G,KIND=qi8), "Greensfun")
       DEALLOCATE(G%PROJ_G)
       NULLIFY(G%PROJ_G)
    ENDIF
  END SUBROUTINE DEALLOCATE_PROJ_G


  SUBROUTINE DEALLOCATE_PROJ_PROJ(G)
    USE ini
    TYPE (greensf) :: G
    IF (ASSOCIATED(G%PROJ_PROJ)) THEN
       CALL DEREGISTER_ALLOCATE(wsgf*SIZE(G%PROJ_PROJ,KIND=qi8), "Greensfun")
       DEALLOCATE(G%PROJ_PROJ)
       NULLIFY(G%PROJ_PROJ)
    ENDIF

  END SUBROUTINE DEALLOCATE_PROJ_PROJ

!***********************************************************************
!
! rotate one electron orbitals given a unitary 
! rotation matrix
! the transformed orbitals are given by
!   <G|a'> = <G|a> U_a a'   (LINCOM_DISTRI)
!
! matrices (such as the Green's function)  in the transformed basis  
! are than given by
!
!   <b'| M |a'>  =  U+_b' b   <b| M |a>  U_a a'
!
! SIGMAW_MAT storing the Green's function in imaginary time 
! are accordingly transformed
!
!***********************************************************************


  SUBROUTINE ROTATE_ORBITALS(W, WDES, CHAM_MAT, GDES_MAT, CORR_MAT, GU_MAT, GO_MAT)
    USE greens_orbital
    IMPLICIT NONE
    TYPE (wavespin)    W
    TYPE (wavedes)     WDES
    GDEF, POINTER :: CHAM_MAT(:,:,:,:)                  ! unitary transformation
    TYPE (greens_mat_des), POINTER :: GDES_MAT
    GDEF, POINTER, OPTIONAL        :: CORR_MAT(:,:,:,:) ! density matrix
    GDEF, POINTER, OPTIONAL        :: GU_MAT(:,:,:,:)   ! time dependent Green function
    GDEF, POINTER, OPTIONAL        :: GO_MAT(:,:,:,:)   ! time dependent Green function
  ! local
    TYPE (wavedes1)    WDES1               ! descriptor for one k-point
    TYPE (wavefuna)    WA                  ! subpointer to part of W
    INTEGER ISP, NK, I, N
    GDEF, ALLOCATABLE :: G1(:)
    COMPLEX(q) :: R(WDES%NB_TOT)

    PROFILING_START('rotate_orbitals')

    IF (PRESENT(CORR_MAT) .OR. PRESENT(GU_MAT) .OR. PRESENT(GO_MAT)) THEN
       ALLOCATE(G1(SIZE(CHAM_MAT,1)))
    ENDIF

    DO ISP=1,WDES%ISPIN
       DO NK=1,WDES%NKPTS
          CALL SETWDES(WDES,WDES1,NK)
          WA=ELEMENTS(W, WDES1, ISP)
          !  redistribution over bands
          CALL REDISTRIBUTE_PROJ( WA)
          CALL REDISTRIBUTE_PW( WA)

          CALL LINCOM_DISTRI_DESC('F',WA%CW_RED(1,1),WA%CPROJ_RED(1,1), CHAM_MAT(1,NK,ISP,1), &
               WDES1%NB_TOTK(ISP), & 
               WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,WDES%NB_TOT, &
               GDES_MAT%COMM_INTRA, NBLK, GDES_MAT%DESC)

          !  back-redistribution over plane wave coefficients
          CALL REDISTRIBUTE_PROJ( WA)
          CALL REDISTRIBUTE_PW( WA)

          ! rotate density matrix
          IF (PRESENT(CORR_MAT)) THEN
          IF (ASSOCIATED(CORR_MAT)) THEN
             DO I=1,SIZE(CORR_MAT,4)
                N=  W%WDES%NB_TOTK(NK, ISP)
                ! G1= U+ G
                CALL PGGEMM( trans, 'N', N, N, N, one,  &
                     CHAM_MAT(1,NK,ISP,1), 1, 1, GDES_MAT%DESC, &
                     CORR_MAT(1, NK, ISP, I), 1, 1, GDES_MAT%DESC, &
                     zero,  &
                     G1(1), 1, 1, GDES_MAT%DESC )
                ! G = G1 U
                CALL PGGEMM( 'N', 'N', N, N, N, one,  &
                     G1(1), 1, 1, GDES_MAT%DESC, &
                     CHAM_MAT(1,NK,ISP,1), 1, 1, GDES_MAT%DESC, &
                     zero,  &
                     CORR_MAT(1, NK, ISP, I), 1, 1, GDES_MAT%DESC )
             ENDDO
          ENDIF
          ENDIF

          ! rotate Green's function in orbital basis
          IF (PRESENT(GU_MAT)) THEN
          IF (ASSOCIATED(GU_MAT)) THEN
             DO I=1,SIZE(GU_MAT,4)
                N=  W%WDES%NB_TOTK(NK, ISP)
                ! G1= U+ G
                CALL PGGEMM( trans, 'N', N, N, N, one,  &
                     CHAM_MAT(1,NK,ISP,1), 1, 1, GDES_MAT%DESC, &
                     GU_MAT(1, NK, ISP, I), 1, 1, GDES_MAT%DESC, &
                     zero,  &
                     G1(1), 1, 1, GDES_MAT%DESC )
                ! G = G1 U
                CALL PGGEMM( 'N', 'N', N, N, N, one,  &
                     G1(1), 1, 1, GDES_MAT%DESC, &
                     CHAM_MAT(1,NK,ISP,1), 1, 1, GDES_MAT%DESC, &
                     zero,  &
                     GU_MAT(1, NK, ISP, I), 1, 1, GDES_MAT%DESC )
             ENDDO
          ENDIF
          ENDIF

         ! rotate Green's function in orbital basis
          IF (PRESENT(GO_MAT)) THEN
          IF (ASSOCIATED(GO_MAT)) THEN
             DO I=1,SIZE(GO_MAT,4)
                N=  W%WDES%NB_TOTK(NK, ISP)
                ! G1= U+ G
                CALL PGGEMM( trans, 'N', N, N, N, one,  &
                     CHAM_MAT(1,NK,ISP,1), 1, 1, GDES_MAT%DESC, &
                     GO_MAT(1, NK, ISP, I), 1, 1, GDES_MAT%DESC, &
                     zero,  &
                     G1(1), 1, 1, GDES_MAT%DESC )
                ! G = G1 U
                CALL PGGEMM( 'N', 'N', N, N, N, one,  &
                     G1(1), 1, 1, GDES_MAT%DESC, &
                     CHAM_MAT(1,NK,ISP,1), 1, 1, GDES_MAT%DESC, &
                     zero,  &
                     GO_MAT(1, NK, ISP, I), 1, 1, GDES_MAT%DESC )
             ENDDO

!             CALL DETERMINE_DIAGONALE_GDEF(WDES%NB_TOTK(NK,ISP), GO_MAT(:,NK,ISP,1), R, GDES_MAT%DESC)
!             WRITE(*,'(10F14.7)') R(1:10)

          ENDIF
          ENDIF

       ENDDO
    ENDDO

    IF (ALLOCATED(G1)) DEALLOCATE(G1)

    PROFILING_STOP('rotate_orbitals')

  END SUBROUTINE ROTATE_ORBITALS

!
! rotate KS orbitals into HF basis
! 

  SUBROUTINE ROTATE_KS_ORBITALS(W, WDES, CHAM_MAT, GDES_MAT)
    USE greens_orbital
    IMPLICIT NONE
    TYPE (wavespin)    W
    TYPE (wavedes)     WDES
    TYPE (greens_mat_des), POINTER :: GDES_MAT
    GDEF, POINTER :: CHAM_MAT(:,:,:,:)                  ! unitary transformation
  ! local
    TYPE (wavedes1)    WDES1               ! descriptor for one k-point
    TYPE (wavefuna)    WA                  ! subpointer to part of W
    INTEGER ISP, NK

    PROFILING_START('rotate_ks_orbitals')

    DO ISP=1,WDES%ISPIN
       DO NK=1,WDES%NKPTS
          CALL SETWDES(WDES,WDES1,NK)
          WA=ELEMENTS(W, WDES1, ISP)
          !  redistribution over bands
          CALL REDISTRIBUTE_PROJ( WA)
          CALL REDISTRIBUTE_PW( WA)

          CALL LINCOM_DISTRI_DESC('F',WA%CW_RED(1,1),WA%CPROJ_RED(1,1), CHAM_MAT(1,NK,ISP,1), &
               WDES1%NB_TOTK(ISP), & 
               WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,WDES%NB_TOT, &
               GDES_MAT%COMM_INTRA, NBLK, GDES_MAT%DESC)

          !  back-redistribution over plane wave coefficients
          CALL REDISTRIBUTE_PROJ( WA)
          CALL REDISTRIBUTE_PW( WA)
       ENDDO
    ENDDO

    PROFILING_STOP('rotate_ks_orbitals')

  END SUBROUTINE ROTATE_KS_ORBITALS

!***********************************************************************
!  
! determine polarizability Chi_q = GO_k1 GU^*_k2
!
!***********************************************************************

  SUBROUTINE CALCULATE_CHI_TAU_FROM_G( W, WGW, GO, GU, GDES, CHI, NQ, NK1, NK2, LATT_CUR, IU6)
    USE ini
    USE dfast
    USE fock
    TYPE (wavespin):: W
    TYPE (wavedes) :: WGW
    TYPE (greensf) :: GU         ! Green function for unoccupied orbitals 
    TYPE (greensf) :: GO         ! Green function for occupied orbitals
    TYPE (greensfdes) :: GDES   
    TYPE (responsefunction) CHI  ! response function
    INTEGER :: NQ, NK1, NK2
    TYPE (latt)        LATT_CUR
    INTEGER :: IU6
  ! local
    TYPE (wavedes1) :: WDESQ
    TYPE (wavedes1) :: WDES1, WDES2
    INTEGER :: NR, NRP, NT, NI, NIS, LMMAXC, L, LP, LM, LMBASE, NAUG, NT_GLOBAL
    INTEGER :: NK1_IN_KPOINTS_FULL_ORIG
    LOGICAL :: LPHASE
    COMPLEX(q) :: CPHASE(GRIDHF%MPLWV)

    GDEF, ALLOCATABLE :: GWORK(:,:)
    GDEF, ALLOCATABLE :: CRHOLM(:) ! augmentation occupancy matrix
    GDEF, ALLOCATABLE :: P_NLM(:,:)
    GDEF, ALLOCATABLE :: P_G_PROJ(:,:)
    GDEF, ALLOCATABLE :: P_PROJ_G(:,:)
    GDEF, ALLOCATABLE :: GO_GU_PROJ(:)
    GDEF, ALLOCATABLE :: P_G_R(:,:)
    GDEF, ALLOCATABLE :: P_R_G(:,:)
    GDEF, ALLOCATABLE :: CHWRK(:)
    TYPE (wavefun1) :: W1, W2

    REAL(q) :: RINPL
    INTEGER :: I
    INTEGER :: ISTAT 
!$  INTEGER :: LASTTYP

    PROFILING_START('calculate_chi_tau_from_g')

    IF (CHI%LREALSTORE) THEN
#ifndef gammareal
       CALL vtutor%bug("internal error in CALCULATE_CHI_TAU_FROM_G: CHWRK is not real", __FILE__, __LINE__)
#endif
    ELSE
#ifdef gammareal
       CALL vtutor%bug("internal error in CALCULATE_CHI_TAU_FROM_G: CHWRK is not", __FILE__, __LINE__)
#endif
    ENDIF
    ALLOCATE(CRHOLM(GDES%NLM_ROW))

    CALL SETWDES(W%WDES,WDES1,NK1)
    CALL SETWDES(W%WDES,WDES2,NK2)

    NK1_IN_KPOINTS_FULL_ORIG=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,NK1),KPOINTS_FULL_ORIG)

    CALL PHASER_HF(GRIDHF, LATT_CUR, FAST_AUG_FOCK, KPOINTS_FULL%VKPT(:,NK1)-KPOINTS_FULL%VKPT(:,NK2))

    CALL SETPHASE(W%WDES%VKPT(:,NK1)-W%WDES%VKPT(:,NK2)-CHI%VKPT(:), &
         GRIDHF,CPHASE(1),LPHASE)

    CALL NEWWAV(W1, WDES1,.TRUE.)
    CALL NEWWAV(W2, WDES2,.TRUE.)
    CALL SETWDES(WGW, WDESQ, NQ )

    IF (WDES1%GRID%RL%NP/= WDESQ%GRID%RL%NP) THEN
       WRITE(*,*) 'internal error in CALCULATE_CHI_TAU_FROM_G: grids distinct ',WDES1%GRID%RL%NP, WDESQ%GRID%RL%NP
    ENDIF
    IF (WDES1%NPRO /= WDESQ%NPRO) THEN
       WRITE(*,*) 'internal error in CALCULATE_CHI_TAU_FROM_G: NPRO distinct ',WDES1%NPRO /= WDESQ%NPRO
    ENDIF
 
     ! FFT GU and GO to real space and set entries G_R, R_PROJ, PROJ_R
     CALL FFT_G(W%WDES, GO, GDES, NK1)
     CALL FFT_G(W%WDES, GU, GDES, NK2)


       ! TODO consider spinors later
       ! DO ISPINOR =0,0 ! GDES%NRSPINORS-1
!=======================================================================
! contract over first direction
!=======================================================================
       ALLOCATE(P_G_PROJ (GDES%RES_NRPLWV_ROW_DATA_POINTS, GDES%NLM_COL),STAT=ISTAT)
       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "CALCULATE_CHI_TAU_FROM_G (P_G_RPROJ) is not able to allocate "//&
            str( ndata * 8._q*GDES%RES_NRPLWV_ROW_DATA_POINTS*GDES%NLM_COL/1024) //&
            " kB of data on MPI rank 0." )
       ENDIF
       ALLOCATE(P_G_R    (GDES%RES_NRPLWV_ROW_DATA_POINTS, GDES%MPLWV_COL),STAT=ISTAT)
       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "CALCULATE_CHI_TAU_FROM_G (P_G_R) is not able to allocate "//&
            str( ndata * 8._q*GDES%RES_NRPLWV_ROW_DATA_POINTS*GDES%MPLWV_COL/1024) //&
            " kB of data on MPI rank 0." )
       ENDIF

       CALL REGISTER_ALLOCATE(wsgf*(SIZE(P_G_R,KIND=qi8)+SIZE(P_G_PROJ,KIND=qi8)), "GGwork")
       CALL DUMP_ALLOCATE_TAG(IU6,"CHI_TAU_FROM_G"); IF (IU6>=0) CALL WFORCE(IU6)
       
       IF (W%WDES%LOVERL) THEN

       AUG_DES%RINPL=1._q

    ! contract over beta, beta' to N'L'M'
       LMBASE =0
       NAUG=0
#ifndef _OPENMP
       NIS =1
       DO NT=1,GDES%NTYP  ! note NT is a local type index
                          ! and might not span all atomic types
          NT_GLOBAL=GDES%NT_GLOBAL(NT)  ! global type index

          ! WDESQ%GRID%MPLWV is the number of complex words required to do in place FFT
          ALLOCATE(P_NLM(GDES%NLM_ROW,GDES%NLM_LMMAX(NT)),  & 
                   GWORK(m_ WDESQ%GRID%MPLWV,GDES%NLM_LMMAX(NT)), GO_GU_PROJ(m_ WDESQ%GRID%MPLWV))

          IF (AUG_DES%LMMAX(NT_GLOBAL)/= GDES%NLM_LMMAX(NT)) THEN
             CALL vtutor%bug("internal error 1 in CALCULATE_CHI_TAU_FROM_G: " // &
                str(AUG_DES%LMMAX(NT_GLOBAL)) // " " // str(GDES%NLM_LMMAX(NT)), __FILE__, __LINE__)
          ENDIF

          LMMAXC=GDES%LMMAX(NT)
          IF (LMMAXC/=0) THEN
             !
             ! main loop over local ion index N (second = column index)
             !
             DO NI=NIS,GDES%NITYP(NT)+NIS-1
#else
       ALLOCATE(GO_GU_PROJ(m_ WDESQ%GRID%MPLWV))
       LASTTYP=0
!!!$OMP PARALLEL DO PRIVATE(NI,NT,NT_GLOBAL,P_NLM,GWORK,GO_GU_PROJ,LMMAXC,LM,L,LP,CRHOLM,RINPL,I) &
!!!$OMP FIRSTPRIVATE(LASTTYP,LMBASE,NAUG)
       DO NI=1,GDES%NIONS

          NT=GDES%NT_LOCAL(NI) ; NT_GLOBAL=GDES%NT_GLOBAL(NT)

          IF (NT /= LASTTYP) THEN
             IF (AUG_DES%LMMAX(NT_GLOBAL)/= GDES%NLM_LMMAX(NT)) THEN
                CALL vtutor%bug("internal error 1 in CALCULATE_SIGMA_TAU: " // &
                   str(AUG_DES%LMMAX(NT_GLOBAL)) // " " // str(GDES%NLM_LMMAX(NT)), __FILE__, __LINE__)
             ENDIF

             IF (ALLOCATED(P_NLM)) DEALLOCATE(P_NLM)
             IF (ALLOCATED(GWORK)) DEALLOCATE(GWORK)

             ALLOCATE(P_NLM(GDES%NLM_ROW,GDES%NLM_LMMAX(NT)),  & 
                      GWORK(m_ WDESQ%GRID%MPLWV,GDES%NLM_LMMAX(NT)))

             LASTTYP=NT
          ENDIF

          LMMAXC=GDES%LMMAX(NT)

          IF (LMMAXC/=0) THEN
             IF (NI>1) THEN
                LMBASE=SUM(GDES%NPRO_ENTRIES(1:NI-1,GDES%COMM%NODE_ME))
                NAUG  =SUM(GDES%NLM_ENTRIES (1:NI-1,GDES%COMM%NODE_ME))
             ENDIF
#endif
                P_NLM=0
                GWORK=0

                ! T_beta, beta'^NLM' G(N' beta, N alpha) G(N' beta', N alpha') T_alpha, alpha'^NLM
                DO L=1,LMMAXC   ! alpha
                DO LP=1,LMMAXC  ! alpha'
                   ! contraction T_beta, beta'^NLM' G(N' beta, N alpha) G(N' beta', N alpha')
                   CALL DEPSUM_TWO_BANDS_RHOLM_TRACE( &
                     GO%PROJ_PROJ(:,LMBASE+L),GU%PROJ_PROJ(:,LMBASE+LP), & 
                     WDES1, AUG_DES, &
                     TRANS_MATRIX_FOCK, CRHOLM, 1._q, WDES1%LOVERL)

                   ! contraction over T_alpha, alpha'^NLM
                   DO LM=1,GDES%NLM_LMMAX(NT)
                      ! index build from NI and LM
                      P_NLM(:,LM)=P_NLM(:,LM)+CRHOLM(:)*TRANS_MATRIX_FOCK(LP,L,LM,NT_GLOBAL)
                   ENDDO
                ENDDO
                ENDDO

                ! calculate G(r', N alpha) G*+(r', N alpha') T_alpha alpha'^NLM
                CALL DEPSUM_VECTOR( GO%R_PROJ, GU%R_PROJ, GWORK, GDES%MPLWV_ROW, TRANS_MATRIX_FOCK(:,:,:,NT_GLOBAL), LMMAXC, LMBASE, GDES%NLM_LMMAX(NT) )

                ! for each atom the matrix P_NLM must be symmetric
                ! at leat at the atom under consideration
                 
                ! now add P(N'L'M', NLM) Q_N'L'M'(r) to density
                IF (WDES1%LOVERL ) THEN
                   DO LM=1,GDES%NLM_LMMAX(NT)
!                     AUG_DES%RINPL=1._q ! multiplicator used by RACC0
                      CALL RACC0_HF(FAST_AUG_FOCK, AUG_DES, P_NLM(1,LM), GWORK(1,LM))
                   ENDDO
                ENDIF
                ! extract data using wave function FFT (impose cutoff at the same time)
                DO LM=1,GDES%NLM_LMMAX(NT)
                   IF (NAUG+LM >  GDES%NLM_COL) THEN
                      CALL vtutor%bug("internal error 2 in CALCULATE_CHI_TAU_FROM_G: " // &
                         str(NAUG+LM) // " " // str(GDES%NLM_COL), __FILE__, __LINE__)
                   ENDIF
                   IF (LPHASE) CALL APPLY_PHASE_GDEF( GRIDHF, CPHASE(1), GWORK(1,LM), GWORK(1,LM))

                   ! collect into P_G_PROJ = P(r, NLM) e (-i G r)
                   ! in the gamma point only mode the entries are 
                   ! interpreted as coefficients of sin and cosine transforms

                   CALL FFTEXT(WDESQ%NGVECTOR,WDESQ%NINDPW(1), &
                        GWORK(1,LM), &
                        P_G_PROJ(1, NAUG+LM),WDESQ%GRID,.FALSE.)

!!                   P_G_PROJ(1:GDES%RES_NRPLWV_ROW_DATA_POINTS, NAUG+LM)=P_G_PROJ(1:GDES%RES_NRPLWV_ROW_DATA_POINTS, NAUG+LM)*(1.0_q/GRIDHF%NPLWV)
                   RINPL=1.0_q/GRIDHF%NPLWV
                   CALL OMP_DSCAL(ndata*GDES%RES_NRPLWV_ROW_DATA_POINTS,RINPL,P_G_PROJ(1,NAUG+LM),1)
                ENDDO
#ifndef _OPENMP
                LMBASE = LMMAXC+LMBASE
                NAUG = NAUG+ GDES%NLM_LMMAX(NT)
             ENDDO
          ENDIF

          NIS = NIS+GDES%NITYP(NT)
          DEALLOCATE(P_NLM, GWORK, GO_GU_PROJ)
#else
          ENDIF
#endif
       ENDDO
!!!$OMP END PARALLEL DO
!$     IF (ALLOCATED(GO_GU_PROJ)) DEALLOCATE(GO_GU_PROJ)
!$     IF (ALLOCATED(P_NLM)) DEALLOCATE(P_NLM)
!$     IF (ALLOCATED(GWORK)) DEALLOCATE(GWORK)
       ENDIF

       ALLOCATE(GWORK(m_ WDESQ%GRID%MPLWV,1))

  ! contract over real space grid points
       DO NRP=1,GDES%MPLWV_COL

          CALL FFTWAV(WDES1%NGVECTOR, WDES1%NINDPW(1), &
               W1%CR(1), &
               GO%G_R(1,NRP),WDES1%GRID)
          CALL FFTWAV(WDES2%NGVECTOR, WDES2%NINDPW(1), &
               W2%CR(1), &
               GU%G_R(1,NRP),WDES2%GRID)
          
          CALL PW_CHARGE_TRACE(WDES1, GWORK(1,1), &
               W1%CR(1), W2%CR(1))

          IF (WDES1%LOVERL ) THEN
             CALL DEPSUM_TWO_BANDS_RHOLM_TRACE( &
                  GO%PROJ_R(:,NRP),GU%PROJ_R(:,NRP), WDES1, AUG_DES, &
                  TRANS_MATRIX_FOCK, CRHOLM, 1._q, WDES1%LOVERL)
             AUG_DES%RINPL=1._q ! multiplicator used by RACC0
             CALL RACC0_HF(FAST_AUG_FOCK, AUG_DES, CRHOLM(1), GWORK(1,1))
          END IF
          IF (LPHASE) CALL APPLY_PHASE_GDEF( GRIDHF, CPHASE(1), GWORK(1,1), GWORK(1,1))
          ! FFT to  reciprocal space and 
          CALL FFTEXT(WDESQ%NGVECTOR,WDESQ%NINDPW(1), &
               GWORK(1,1), &
               P_G_R(1,  NRP),WDESQ%GRID,.FALSE.)
!!          P_G_R(1:GDES%RES_NRPLWV_ROW_DATA_POINTS, NRP)=P_G_R(1:GDES%RES_NRPLWV_ROW_DATA_POINTS, NRP)*(1.0_q/WDESQ%GRID%NPLWV)
          RINPL=1.0_q/WDESQ%GRID%NPLWV
          CALL OMP_DSCAL(ndata*GDES%RES_NRPLWV_ROW_DATA_POINTS,RINPL,P_G_R(1,NRP),1)
       ENDDO

       CALL  RELEASE_FFT_G( GU)
       CALL  RELEASE_FFT_G( GO)
!=======================================================================
! contract over second direction, which
! now involves mostly an FFT in the second direction
!=======================================================================
       ! start with a transpose  (and conjugation for complex orbitals)
       ALLOCATE(P_PROJ_G (GDES%NLM_ROW,    GDES%RES_NRPLWV_COL_DATA_POINTS ), STAT=ISTAT)
       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "CALCULATE_CHI_TAU_FROM_G (P_PROJ_G) is not able to allocate "//&
            str( ndata * 8._q*GDES%NLM_ROW*GDES%RES_NRPLWV_COL_DATA_POINTS/1024) //&
            " kB of data on MPI rank 0." )
       ENDIF
       ALLOCATE( P_R_G    (GDES%MPLWV_ROW,  GDES%RES_NRPLWV_COL_DATA_POINTS), STAT=ISTAT)
       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "CALCULATE_CHI_TAU_FROM_G (P_R_G) is not able to allocate "//&
            str( ndata * 8._q*GDES%MPLWV_ROW*GDES%RES_NRPLWV_COL_DATA_POINTS/1024) //&
            " kB of data on MPI rank 0." )
       ENDIF

       CALL REGISTER_ALLOCATE(wsgf*(SIZE(P_R_G,KIND=qi8)+SIZE(P_PROJ_G,KIND=qi8)), "GGwork")
       CALL DUMP_ALLOCATE_TAG(IU6,"CHI_TAU_FROM_G"); IF (IU6>=0) CALL WFORCE(IU6)

       CALL TRANSPOSE_G_R_RESPONSE   (P_G_R   , P_R_G, GDES)
       CALL TRANSPOSE_G_PROJ_RESPONSE(P_G_PROJ, P_PROJ_G, GDES)

       CALL DEREGISTER_ALLOCATE(wsgf*(SIZE(P_G_R,KIND=qi8)+SIZE(P_G_PROJ,KIND=qi8)), "GGwork")
       DEALLOCATE(P_G_PROJ, P_G_R)


       ALLOCATE(CHWRK(GDES%RES_NRPLWV_ROW_DATA_POINTS))
DOESI  CHWRK=0
       DO NRP=1,GDES%RES_NRPLWV_COL_DATA_POINTS
          ! copy data into GWORK
!!          GWORK=0
!!          GWORK(1:GDES%MPLWV_ROW,1)=P_R_G(1:GDES%MPLWV_ROW,NRP)
          CALL OMP_DCOPY(ndata*GDES%MPLWV_ROW,P_R_G(1,NRP),1,GWORK(1,1),1)
!         IF (GDES%MPLWV_ROW < m_ WDESQ%GRID%MPLWV) GWORK(GDES%MPLWV_ROW+1:m_ WDESQ%GRID%MPLWV,1)=0

          ! add augmentation part  P(NLM, r) Q_NLM(r')
          IF (WDES1%LOVERL) THEN
             AUG_DES%RINPL=1._q ! multiplicator used by RACC0
             CALL RACC0_HF(FAST_AUG_FOCK, AUG_DES, P_PROJ_G(1, NRP), GWORK(1,1))
          END IF
          ! store final results in RESPONSER
          IF (CHI%LREALSTORE) THEN
             IF (LPHASE) CALL APPLY_PHASE_GDEF( GRIDHF, CPHASE(1), GWORK(1,1) , GWORK(1,1))
             CALL FFTEXT(WDESQ%NGVECTOR,WDESQ%NINDPW(1), &
                  GWORK(1,1), &
                  CHWRK(1),WDESQ%GRID,.FALSE.)
!!             CHI%RESPONSER(1:GDES%RES_NRPLWV_ROW_DATA_POINTS, NRP, 1)= & 
!!                  CHI%RESPONSER(1:GDES%RES_NRPLWV_ROW_DATA_POINTS, NRP, 1)+ & 
!!                  CHWRK(1:GDES%RES_NRPLWV_ROW_DATA_POINTS)*(1.0_q/WDESQ%GRID%NPLWV)
             RINPL=1.0_q/WDESQ%GRID%NPLWV
             CALL OMP_DAXPY(GDES%RES_NRPLWV_ROW_DATA_POINTS,RINPL,CHWRK(1),ndata,CHI%RESPONSER(1,NRP,1),1)
          ELSE
             IF (LPHASE) CALL APPLY_PHASE_GDEF( GRIDHF, CPHASE(1), GWORK(1,1) , GWORK(1,1))
             CALL FFTEXT(WDESQ%NGVECTOR,WDESQ%NINDPW(1), &
                  GWORK(1,1), &
                  CHWRK(1),WDESQ%GRID,.FALSE.)
!!             CHI%RESPONSEFUN(1:GDES%RES_NRPLWV_ROW_DATA_POINTS, NRP, 1)= &
!!                  CHI%RESPONSEFUN(1:GDES%RES_NRPLWV_ROW_DATA_POINTS, NRP, 1)+ & 
!!                  CHWRK(1:GDES%RES_NRPLWV_ROW_DATA_POINTS)*(1.0_q/WDESQ%GRID%NPLWV)
             RINPL=1.0_q/WDESQ%GRID%NPLWV
#ifndef gammareal
             CALL OMP_DAXPY(ndata*GDES%RES_NRPLWV_ROW_DATA_POINTS,RINPL,CHWRK(1),1,CHI%RESPONSEFUN(1,NRP,1),1)
#else
             CALL OMP_DAXPY(GDES%RES_NRPLWV_ROW_DATA_POINTS,RINPL,CHWRK(1),1,CHI%RESPONSEFUN(1,NRP,1),ndata)
#endif
          ENDIF
       ENDDO
       DEALLOCATE(CHWRK)

       CALL DEREGISTER_ALLOCATE(wsgf*(SIZE(P_R_G,KIND=qi8)+SIZE(P_PROJ_G,KIND=qi8)), "GGwork")
       DEALLOCATE( P_R_G, P_PROJ_G)

       DEALLOCATE(GWORK)


    DEALLOCATE(CRHOLM)
    CALL DELWAV(W1, .TRUE.)
    CALL DELWAV(W2, .TRUE.)

    PROFILING_STOP('calculate_chi_tau_from_g')

  END SUBROUTINE CALCULATE_CHI_TAU_FROM_G


!***********************************************************************
!
! FFT_G_SUPER fourier transform the Green function along first index from 
!  g' to r'
!
!  G(r',g)  =  exp(i k r') \sum_g' exp(i g'r') G(g', g)     -> G%G_R
!  G(r',N a)=  exp(i k r') \sum_g' exp(i g'r') G(g', N a)   -> G%PROJ_R
! 
! a phase factor exp(i k r') is also multiplied in
! data are distributed over g and N respectively
! afterwards the arrays are transposed
!
!  G(g, r') =  G*(r', g) 
!  G(N a, r') =  G*(r', N a)
!
! data are now distributed over r'
!
! except for
! ) the extra phase factor the routine is equivalent to FFT_G
! ) the routine deallocates G%R_PROJ
!
!***********************************************************************

  SUBROUTINE FFT_G_SUPER( WDES, G, GDES, NK)
    USE ini
    USE twoelectron4o
    TYPE (wavedes),INTENT(IN)    :: WDES
    TYPE (greensf),INTENT(INOUT) :: G
    TYPE (greensfdes),INTENT(IN) :: GDES         
  ! local
    INTEGER :: ICOL
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1
    INTEGER :: NK
    LOGICAL LPHASE           ! determines whether multiplication with phase factor is required
    COMPLEX(q) :: CPHASE(WDES%GRID%RL%NP)
    INTEGER :: NRPLWV_COL_DATA_POINTS, NRPLWV_ROW_DATA_POINTS
    INTEGER :: ISTAT

    PROFILING_START('fft_g_super')

    IF (ASSOCIATED(G%R_G).OR.ASSOCIATED(G%R_PROJ) .OR. ASSOCIATED(G%PROJ_R) ) THEN
       CALL vtutor%bug("internal error in FFT_G: G%R_G is associated " // str(ASSOCIATED(G%R_G)) // &
          " " // str(ASSOCIATED(G%R_PROJ)) // " " // str(ASSOCIATED(G%PROJ_R)), __FILE__, __LINE__)
    ENDIF
    NRPLWV_ROW_DATA_POINTS=GDES%NRPLWV_ROW_DATA_POINTS
    NRPLWV_COL_DATA_POINTS=GDES%NRPLWV_COL_DATA_POINTS
    IF (ASSOCIATED(GDES%LUSEINV)) THEN
       IF (GDES%LUSEINV(NK)) THEN
          NRPLWV_ROW_DATA_POINTS=GDES%NRPLWV_ROW_DATA_POINTS_NK(NK)
          NRPLWV_COL_DATA_POINTS=GDES%NRPLWV_COL_DATA_POINTS_NK(NK)
       ENDIF
    ENDIF

    ALLOCATE(G%R_G   (GDES%MPLWV_ROW, NRPLWV_COL_DATA_POINTS),STAT=ISTAT)
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "FFT_G_SUPER (G_R) is not able to allocate "//&
         str( ndata * 8._q*GDES%MPLWV_ROW*NRPLWV_COL_DATA_POINTS/1024)  // &
         " kB of data on MPI rank 0." )
    ENDIF

! Fock FFT is used
    CALL SETWDES(WDES, WDES1, NK)
! set exp(i k r')
    CALL SETPHASE_NOCHK(WDES%VKPT(1:3,NK), WDES1%GRID, CPHASE, LPHASE)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(W1,ICOL)
    CALL NEWWAV(W1, WDES1,.TRUE.)
!$OMP DO
    DO ICOL=1, NRPLWV_COL_DATA_POINTS
       ! possibly complex to real FFT
       ! nasty CR is complex valued but only the real components
       ! matter in the Gamma-point only case
       CALL FFTWAV(WDES1%NGVECTOR, WDES1%NINDPW(1), &
            W1%CR(1), &
            G%GG(1,ICOL),WDES1%GRID)

       IF (LPHASE) CALL APPLY_PHASE_INPLACE( WDES1%GRID, CPHASE, W1%CR(1))

!      G%R_G(1:W1%WDES1%GRID%RL%NP, ICOL)=W1%CR(1:W1%WDES1%GRID%RL%NP)
!      G%R_G(W1%WDES1%GRID%RL%NP+1:GDES%MPLWV_ROW, ICOL)=0  ! pad with zeros (just in case)
#ifndef gammareal
       CALL OMP_DCOPY(2*W1%WDES1%GRID%RL%NP,W1%CR(1),1,G%R_G(1,ICOL),1)
#else
       CALL OMP_DCOPY(  W1%WDES1%GRID%RL%NP,W1%CR(1),2,G%R_G(1,ICOL),1)
#endif
    ENDDO
!$OMP END DO
!$  CALL DELWAV(W1, .TRUE.)    
!$OMP END PARALLEL
    ALLOCATE(G%G_R(NRPLWV_ROW_DATA_POINTS, GDES%MPLWV_COL), STAT=ISTAT)
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "FFT_G_SUPER (G_R) is not able to allocate "//&
         str( ndata * 8._q*NRPLWV_ROW_DATA_POINTS*GDES%MPLWV_COL/1024)  // " kB of data on MPI rank 0." )
    ENDIF
    CALL REGISTER_ALLOCATE(wsgf*(SIZE(G%G_R,KIND=qi8)), "Greensfun")

    CALL TRANSPOSE_R_G( G%R_G, G%G_R, GDES, NK)

    ! we do not need R_G anymore, so destroy it
    DEALLOCATE(G%R_G)
    NULLIFY(G%R_G)

    ALLOCATE(G%R_PROJ(GDES%MPLWV_ROW, GDES%NPRO_COL),STAT=ISTAT)
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "FFT_G_SUPER (G_R_PROJ) is not able to allocate "//&
         str( ndata * 8._q*GDES%MPLWV_ROW*GDES%NPRO_COL/1024)  // " kB of data on MPI rank 0." )
    ENDIF

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(W1,ICOL)
!$  CALL NEWWAV(W1, WDES1,.TRUE.)
!$OMP DO
    DO ICOL=1,  GDES%NPRO_COL
       CALL FFTWAV(WDES1%NGVECTOR, WDES1%NINDPW(1), &
            W1%CR(1), G%G_PROJ(1,ICOL),WDES1%GRID)

       IF (LPHASE) CALL APPLY_PHASE_INPLACE( WDES1%GRID, CPHASE, W1%CR(1))

!      G%R_PROJ(1:W1%WDES1%GRID%RL%NP, ICOL)=W1%CR(1:W1%WDES1%GRID%RL%NP)
!      G%R_PROJ(W1%WDES1%GRID%RL%NP+1:GDES%MPLWV_ROW, ICOL)=0  ! pad with zeros (just in case)
#ifndef gammareal
       CALL OMP_DCOPY(2*W1%WDES1%GRID%RL%NP,W1%CR(1),1,G%R_PROJ(1,ICOL),1)
#else
       CALL OMP_DCOPY(  W1%WDES1%GRID%RL%NP,W1%CR(1),2,G%R_PROJ(1,ICOL),1)
#endif
    ENDDO
!$OMP END DO
    CALL DELWAV(W1, .TRUE.)    
!$OMP END PARALLEL
       
    ALLOCATE(G%PROJ_R(GDES%NPRO_ROW,  GDES%MPLWV_COL), STAT=ISTAT)
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "FFT_G_SUPER (PROJ_R) is not able to allocate "//&
         str( ndata * 8._q*GDES%NPRO_ROW*GDES%MPLWV_COL/1024)// " kB of data on MPI rank 0." )
    ENDIF
    CALL REGISTER_ALLOCATE(wsgf*(SIZE(G%PROJ_R,KIND=qi8)), "Greensfun")

    CALL TRANSPOSE_R_PROJ( G%R_PROJ, G%PROJ_R, GDES)

    ! we do not need R_PROJ anymore, so destroy it
    DEALLOCATE(G%R_PROJ)
    NULLIFY(G%R_PROJ)

#ifdef debug
    WRITE(77+GDES%COMM%NODE_ME,*) GDES%NPRO_ROW, GDES%NPRO_COL
    WRITE(77+GDES%COMM%NODE_ME,*)
    DO ICOL=1,16
       CALL FFTWAV(WDES1%NGVECTOR, WDES1%NINDPW(1), &
            W1%CR(1),G%G_R(1,ICOL),WDES1%GRID)
       WRITE(77+GDES%COMM%NODE_ME,'(16F7.3)') REAL(W1%CR(1:16))   
    ENDDO
    WRITE(77+GDES%COMM%NODE_ME,*)
    WRITE(77+GDES%COMM%NODE_ME,'(16F7.3)') G%PROJ_PROJ(1:16,1:GDES%NPRO_COL)
    WRITE(77+GDES%COMM%NODE_ME,*)
    WRITE(77+GDES%COMM%NODE_ME,'(16F7.3)') G%PROJ_R(1:16,1:16)
#endif

    PROFILING_STOP('fft_g_super')

  END SUBROUTINE FFT_G_SUPER


!***********************************************************************
!
! build supercell "orbitals" from the primitive cell entries
! 
!  G(g+k, X)     =  G_k(g, X) 
!  G(N' alpha, X) =  sum_k exp^(i k R_N' ) G_k(N alpha, X)
!
! uses the previously set up G%G_R and G%PROJ_R of FFT_G_SUPER
!
! as input the routines expect the currently considered position
! index X=r' in G_R and PROJ_R
! (maybe we can subindex)
!
!***********************************************************************

  SUBROUTINE  COLLECT_ORBITALS_SUPER( S, WDES, GO_GK_R, GO_PROJ_R, GO, NRP)
    USE constant

    TYPE (supercell), POINTER, INTENT(IN) :: S   ! supercell descriptor
    TYPE (wavedes)                        :: WDES
    COMPLEX(q), INTENT(OUT) :: GO_GK_R(:) ! G(g+k, X) 
    GDEF, INTENT(OUT) :: GO_PROJ_R(:)     ! G(N' alpha, X) 
    TYPE (greensf), INTENT(IN) :: GO(:)   ! Green function for each k-point
    INTEGER :: NRP                        ! place from which to collect the Greens function
  ! local 
    INTEGER :: NK, N, NT, NI, NIS, NIS_DEST, LMBASE, LMMAXC, LMBASE_DEST, NP, NGVECTOR
    REAL(q) :: SCALE
    LOGICAL :: LUSEINV

    PROFILING_START('collect_orbitals_super')

#ifndef RPAgamma
    IF (SIZE(GO) /= S%NREP) THEN
       CALL vtutor%bug("internal error in COLLECT_ORBITALS_SUPER: mismatched G and k-points " // &
          str(SIZE(GO)) // " " // str(S%NREP), __FILE__, __LINE__)
    ENDIF
#endif

    GO_GK_R=0
    GO_PROJ_R=0
    DO NK=1,SIZE(GO)
       NGVECTOR=S%NGVECTOR(NK)
       LUSEINV=.FALSE.

       IF (ASSOCIATED(S%LUSEINV)) THEN
          IF (S%LUSEINV(NK)) THEN
             NGVECTOR=S%NGVECTOR_INV(NK)
             LUSEINV=.TRUE.
          ENDIF
       ENDIF
             
       IF (.NOT. LUSEINV .AND. WDES%NGVECTOR(NK)/=NGVECTOR) THEN
          CALL vtutor%bug("internal error in COLLECT_ORBITALS_SUPER: NGVECTOR not correct " // &
             str(NK) // " " // str(WDES%NGVECTOR(NK)) // " " // str(NGVECTOR), __FILE__, __LINE__)
       ENDIF

       IF (NGVECTOR> SIZE(GO(NK)%G_R,1)) THEN
          WRITE(0,*) 'internal error in COLLECT_ORBITALS_SUPER: index in G_R exceeded:', & 
               NGVECTOR, SIZE(GO(NK)%G_R,1)
       ENDIF

#ifdef RPAgamma
       IF (.NOT. LUSEINV) THEN
!DIR$ IVDEP
!OCL NOVREC
         DO N=1,NGVECTOR
          IF (S%INDEX(N,NK)/=0) THEN
             ! in half grid mode the PW coefficients on reduced grid need to be multiplied with FFTSCA(N,1)
             GO_GK_R(S%INDEX(N,NK))=GO(NK)%G_R(N, NRP)*S%GRID%FFTSCA(S%INDEX(N,NK),1)
          ELSE IF (S%INDEX_INV(N,NK)/=0) THEN
             ! if INDEX is missing use inversion symmetry
             GO_GK_R(S%INDEX_INV(N,NK))=GCONJG(GO(NK)%G_R(N, NRP))*S%GRID%FFTSCA(S%INDEX_INV(N,NK),1)
          ENDIF
         ENDDO
       ELSE
!DIR$ IVDEP
!OCL NOVREC
        ! only half the coefficients are available
        DO N=1,NGVECTOR
          IF (S%IND_IF_INV(N,NK)/=0) THEN
             ! in half grid mode the PW coefficients on reduced grid need to be multiplied with FFTSCA(N,1)
             GO_GK_R(S%IND_IF_INV(N,NK))=GO(NK)%G_R(N, NRP)*S%GRID%FFTSCA(S%IND_IF_INV(N,NK),1)
          ENDIF
       ENDDO
     ENDIF
#else
!DIR$ IVDEP
!OCL NOVREC
      DO N=1,NGVECTOR
          GO_GK_R(S%INDEX(N,NK))=GO(NK)%G_R(N, NRP)
      ENDDO
#endif
    ENDDO

    NIS=1         ! loop over all ions
    LMBASE=0      ! index into PROJ_R
    NIS_DEST=1    ! ion in the supercell
    LMBASE_DEST=0 ! index into GO_PROJ_R

    ! scaling of stored CPROJ is somewhat odd and involves a 1/sqrt(volume) term
    SCALE= SQRT(1.0_q/S%NREP)

    DO NT=1,WDES%NTYP
       LMMAXC=WDES%LMMAX(NT)
       DO NI=1,WDES%NITYP(NT)
          DO NP=1,S%NREP  ! loop over replicated ions
             ! loop over all k-points in the IRZ
             ! VKPT is in reciprocal coordinates of original cell, POSION in supercell coordinates
             ! this puts in the proper phase factor e^(i k R_N) corresponding
             ! to the "orbitals" constructed above (tested carefully)
             DO NK=1,SIZE(GO)
                GO_PROJ_R(LMBASE_DEST+1:LMBASE_DEST+LMMAXC)=GO_PROJ_R(LMBASE_DEST+1:LMBASE_DEST+LMMAXC)+ &
                  ! position of atom in fractional coordinates of original cell
                  ! is required here
#ifdef RPAgamma
                  ! coefficients are real valued, and k-points are reduced
                  REAL(EXP(CITPI*(WDES%VKPT(1,NK)*S%POSION(1,NIS_DEST)+  &
                             WDES%VKPT(2,NK)*S%POSION(2,NIS_DEST)+  &
                             WDES%VKPT(3,NK)*S%POSION(3,NIS_DEST))) & 
                   *GO(NK)%PROJ_R(LMBASE+1:LMBASE+LMMAXC, NRP),q)*SCALE*S%KWEIGHT(NK)
#else
                  EXP(CITPI*(WDES%VKPT(1,NK)*S%POSION(1,NIS_DEST)+  &
                             WDES%VKPT(2,NK)*S%POSION(2,NIS_DEST)+  &
                             WDES%VKPT(3,NK)*S%POSION(3,NIS_DEST))) & 
                   *GO(NK)%PROJ_R(LMBASE+1:LMBASE+LMMAXC, NRP)*SCALE
#endif
             ENDDO
             LMBASE_DEST=LMBASE_DEST+LMMAXC
             NIS_DEST=NIS_DEST+1
          ENDDO
          LMBASE =LMBASE+LMMAXC
          NIS=NIS+1
       ENDDO
    ENDDO

    PROFILING_STOP('collect_orbitals_super')

  END SUBROUTINE COLLECT_ORBITALS_SUPER


!***********************************************************************
!
!  take the response function in supercell and distribute it to k-point
!  response function arrays 
! 
!  G(g+k, X)     =  G_k(g, X) 
!  G(N' alpha, X) =  sum_k exp^(i k R_N' ) G_k(N alpha, X)
!       ! for green's function
!       CALL COLLECT_ORBITALS_SUPER(S, W%WDES, W1%CW, W1%CPROJ, GU(T%NPOINTS,:), NRP)
!       ! for screening, use chi descriptor
!       CALL COLLECT_RESPONSE_SUPER(S, WGW, WQ%CW, WQ%CPROJ, P_G_R, NRP)
!
!***********************************************************************

  SUBROUTINE  COLLECT_RESPONSE_SUPER( S, WGW, P_GK_R, P_TMP, NRP)
    USE constant
    TYPE (supercell), POINTER, INTENT(IN) :: S   ! supercell descriptor
    TYPE (wavedes)                        :: WGW
    COMPLEX(q), INTENT(OUT) :: P_GK_R(:) ! W(g+k, X) 
    GDEF, INTENT(IN)   :: P_TMP(:,:,:)   ! W for each k-point: g', r, nq 
    INTEGER :: NRP                       ! index for second dimension of P_TMP
  ! local 
    INTEGER :: NK, N, NT, NI, NIS, NIS_DEST, LMBASE, LMMAXC, LMBASE_DEST, NP
    REAL(q) :: SCALE

#ifndef RPAgamma
    IF (SIZE(P_TMP,3) /= S%NREP) THEN
       CALL vtutor%bug("internal error in COLLECT_RESPONSE_SUPER: mismatched G and k-points " // &
          str(SIZE(P_TMP,3)) // " " // str(S%NREP), __FILE__, __LINE__)
    ENDIF 
#endif
          
    P_GK_R=0
    DO NK=1,SIZE(P_TMP,3)
       IF (WGW%NGVECTOR(NK)/=S%NGVECTOR_RES(NK)) THEN
       !%IF (WGW%NGVECTOR(NK)/=S%NGVECTOR(NK)) THEN
          CALL vtutor%bug("internal error in COLLECT_RESPONSE_SUPER: NGVECTOR_RES not correct " // &
             str(NK) // " " // str(WGW%NGVECTOR(NK)) // " " // str(S%NGVECTOR_RES(NK)), __FILE__, __LINE__)
       ENDIF
          
       IF (S%NGVECTOR_RES(NK)> SIZE(P_TMP,1)) THEN
       !%IF (S%NGVECTOR(NK)> SIZE(P_TMP,1)) THEN
          WRITE(0,*) 'internal error in COLLECT_RESPONSE_SUPER: index in G_R exceeded:', &
               S%NGVECTOR_RES(NK), SIZE(P_TMP,1)
       ENDIF

       DO N=1,S%NGVECTOR_RES(NK)
#ifdef RPAgamma
          ! in half grid mode the PW coefficients on reduced grid need to be multiplied with FFTSCA(N,1)
          IF (S%INDEX_RES(N,NK)/=0) THEN
             P_GK_R(S%INDEX_RES(N,NK))=P_TMP(N,NRP,NK)*S%GRID_RES%FFTSCA(S%INDEX_RES(N,NK),1)
          ELSE IF (S%INDEX_RES_INV(N,NK)/=0) THEN
             P_GK_R(S%INDEX_RES_INV(N,NK))=GCONJG(P_TMP(N,NRP,NK))*S%GRID_RES%FFTSCA(S%INDEX_RES_INV(N,NK),1)
          ENDIF
#else
          P_GK_R(S%INDEX_RES(N,NK))=P_TMP(N,NRP,NK)
#endif
       ENDDO
    ENDDO
  END SUBROUTINE COLLECT_RESPONSE_SUPER


!***********************************************************************
!
! build supercell "orbital" characters from the primitive cell entries
! 
!  G(g+k, N beta)     =  G_k(g, N beta) 
!  G(N'' alpha, N beta) =  sum exp^(i k R_N'' ) G_k(N' alpha, N beta)
!
! uses the previously set up G%G_PROJ and G%PROJ_PROJ
! as determined by the routine CALCULATE_G_RECIPROCAL
!
! as input the routines expect the currently considered storage positions
! LMBASE_+1: LMBASE_+ LMMAXC_ and the corresponding position POSION 
! 
! important note: to calculate the supercell characters a phase
! factor e^-ik r_N needs to be added in for the ion at position N,
! because VASP stores the orbital characters without this phase factor
! (see nonl.F). When going to the supercell quantities  (at Gamma)
! this Bloch phase factor is missing.
!
!***********************************************************************

  SUBROUTINE COLLECT_PROJECTORS_SUPER(S, WDES, GO_GK_PROJ, GO_PROJ_PROJ, GO, LMBASE_, LMMAXC_, POSION)
    USE constant
    TYPE (supercell), POINTER, INTENT(IN) :: S
    TYPE (wavedes)                        :: WDES
    GDEF, INTENT(OUT) :: GO_GK_PROJ(:,:)  
    GDEF, INTENT(OUT) :: GO_PROJ_PROJ(:,:) 
    TYPE (greensf), INTENT(IN) :: GO(:)    ! Green function for each k-points
    INTEGER :: LMBASE_                     !
    INTEGER :: LMMAXC_
    REAL(q) :: POSION(3)
  ! local 
    INTEGER :: NK, N, NT, NI, NIS, NIS_DEST, LMBASE, LMMAXC, LMBASE_DEST, NP
    REAL(q) :: SCALE
    COMPLEX(q) :: CPHASE

    PROFILING_START('collect_projectors_super')

#ifndef RPAgamma
    IF (SIZE(GO) /= S%NREP) THEN
       CALL vtutor%bug("internal error in COLLECT_PROJECTORS_SUPER: mismatched G and k-points " // &
          str(SIZE(GO)) // " " // str(S%NREP), __FILE__, __LINE__)
    ENDIF
#endif

    GO_GK_PROJ=0
    GO_PROJ_PROJ=0

    DO NK=1,SIZE(GO)
       IF (WDES%NGVECTOR(NK)/=S%NGVECTOR(NK)) THEN
          CALL vtutor%bug("internal error in COLLECT_PROJECTORS_SUPER: NGVECTOR not correct " // &
             str(NK) // " " // str(WDES%NGVECTOR(NK)) // " " // str(S%NGVECTOR(NK)), __FILE__, __LINE__)
       ENDIF

       IF (S%NGVECTOR(NK)> SIZE(GO(NK)%G_PROJ,1)) THEN
          WRITE(0,*) 'internal error in COLLECT_PROJECTORS_SUPER: index in G_R exceeded:', & 
               S%NGVECTOR(NK), SIZE(GO(NK)%G_PROJ,1)
       ENDIF
       ! in VASP the phase factor e^i k r_N is missing (see also nonl.F)
       ! here we need to include it to get the real Bloch Green's function
       ! the second index in G is conjugated
! CHECKED AGAINST APPLY_PHASE_NLM  SAME RESULTS
       CPHASE=EXP(CITPI*(-WDES%VKPT(1,NK)*POSION(1) &
                         -WDES%VKPT(2,NK)*POSION(2) &
                         -WDES%VKPT(3,NK)*POSION(3)))
       DO N=1,S%NGVECTOR(NK)
#ifdef RPAgamma
          IF (S%INDEX(N,NK)/=0) THEN
             ! in half grid mode the PW coefficients on reduced grid need to be multiplied with FFTSCA(N,1)
             ! kweight because NK loop is reduced
             GO_GK_PROJ(S%INDEX(N,NK),1:LMMAXC_)=GO(NK)%G_PROJ(N, LMBASE_+1:LMBASE_+LMMAXC_)*CPHASE*S%GRID%FFTSCA(S%INDEX(N,NK),1)
          ELSE IF (S%INDEX_INV(N,NK)/=0) THEN
             GO_GK_PROJ(S%INDEX_INV(N,NK),1:LMMAXC_)=CONJG(GO(NK)%G_PROJ(N, LMBASE_+1:LMBASE_+LMMAXC_)*CPHASE)*S%GRID%FFTSCA(S%INDEX_INV(N,NK),1)
          ENDIF
#else
          GO_GK_PROJ(S%INDEX(N,NK),1:LMMAXC_)=GO(NK)%G_PROJ(N, LMBASE_+1:LMBASE_+LMMAXC_)*CPHASE
#endif
       ENDDO
    ENDDO

    NIS=1         ! loop over all ions
    LMBASE=0      ! index into PROJ_R
    NIS_DEST=1    ! ion in the supercell
    LMBASE_DEST=0 ! index into GO_PROJ_R

    ! scaling of stored CPROJ is somewhat odd and involves a 1/sqrt(volume) term
    SCALE= SQRT(1.0_q/S%NREP)

    DO NT=1,WDES%NTYP
       LMMAXC=WDES%LMMAX(NT)
       DO NI=1,WDES%NITYP(NT)
          DO NP=1,S%NREP  ! loop over replicated ions
             ! loop over all k-points in the IRZ
             ! VKPT is in reciprocal coordinates of original cell, POSION in supercell coordinates
             DO NK=1,SIZE(GO)
#ifdef RPAgamma
                ! coefficients are real valued, and kweight because k-points are reduced
                GO_PROJ_PROJ(LMBASE_DEST+1:LMBASE_DEST+LMMAXC, 1:LMMAXC_)= &
                GO_PROJ_PROJ(LMBASE_DEST+1:LMBASE_DEST+LMMAXC, 1:LMMAXC_)+ &
                     REAL(EXP(CITPI*(WDES%VKPT(1,NK)*S%POSION(1,NIS_DEST) &
                            +WDES%VKPT(2,NK)*S%POSION(2,NIS_DEST) &
                            +WDES%VKPT(3,NK)*S%POSION(3,NIS_DEST) &
                            -WDES%VKPT(1,NK)*POSION(1) &
                            -WDES%VKPT(2,NK)*POSION(2) &
                            -WDES%VKPT(3,NK)*POSION(3))) &
                    *GO(NK)%PROJ_PROJ(LMBASE+1:LMBASE+LMMAXC, LMBASE_+1:LMBASE_+LMMAXC_),q)*SCALE*S%KWEIGHT(NK)
#else
                GO_PROJ_PROJ(LMBASE_DEST+1:LMBASE_DEST+LMMAXC, 1:LMMAXC_)= &
                GO_PROJ_PROJ(LMBASE_DEST+1:LMBASE_DEST+LMMAXC, 1:LMMAXC_)+ &
                     EXP(CITPI*(WDES%VKPT(1,NK)*S%POSION(1,NIS_DEST) &
                            +WDES%VKPT(2,NK)*S%POSION(2,NIS_DEST) &
                            +WDES%VKPT(3,NK)*S%POSION(3,NIS_DEST) &
                            -WDES%VKPT(1,NK)*POSION(1) &
                            -WDES%VKPT(2,NK)*POSION(2) &
                            -WDES%VKPT(3,NK)*POSION(3))) &
                    *GO(NK)%PROJ_PROJ(LMBASE+1:LMBASE+LMMAXC, LMBASE_+1:LMBASE_+LMMAXC_)*SCALE
#endif
             ENDDO
             LMBASE_DEST=LMBASE_DEST+LMMAXC
             NIS_DEST=NIS_DEST+1
          ENDDO
          LMBASE =LMBASE+LMMAXC
          NIS=NIS+1
       ENDDO
    ENDDO

    PROFILING_STOP('collect_projectors_super')

  END SUBROUTINE COLLECT_PROJECTORS_SUPER

!***********************************************************************
!
! transform projector part in supercell by summing over all 
! contributions with an appropriate phase factor
!
! as input the routines expect the currently considered storage positions
! LMBASE_+1: LMBASE_+ LMMAXC_ and the corresponding position POSION 
! 
! important note: to calculate the supercell characters a phase
! factor e^-ik r_N needs to be added in for the ion at position N,
! because VASP stores the orbital characters without this phase factor
! (see nonl.F). When going to the supercell quantities  (at Gamma)
! this Bloch phase factor is missing.
!
!***********************************************************************

  SUBROUTINE COLLECT_PROJ_R_SUPER(S, WDES, GU_PROJ_R, GU, NRP)
    USE constant
    TYPE (supercell), POINTER, INTENT(IN) :: S
    TYPE (wavedes)                        :: WDES
    GDEF, INTENT(OUT) :: GU_PROJ_R(:,:) 
    TYPE (greensf), INTENT(IN) :: GU(:)    ! Green function for each k-point
    INTEGER :: NRP
  ! local 
    INTEGER :: NK, N, NT, NI, NIS, NIS_DEST, LMBASE, LMMAXC, LMBASE_DEST, NP
    REAL(q) :: SCALE

#ifndef RPAgamma
    IF (SIZE(GU) /= S%NREP) THEN
       CALL vtutor%bug("internal error in COLLECT_PROJ_R_SUPER: mismatched G and k-points " // &
          str(SIZE(GU)) // " " // str(S%NREP), __FILE__, __LINE__)
    ENDIF
#endif

    GU_PROJ_R=0

    NIS=1         ! loop over all ions
    LMBASE=0      ! index into PROJ_R
    NIS_DEST=1    ! ion in the supercell
    LMBASE_DEST=0 ! index into GO_PROJ_R

    ! scaling of stored CPROJ is somewhat odd and involves a 1/sqrt(volume) term
    SCALE= SQRT(1.0_q/S%NREP)

    DO NT=1,WDES%NTYP
       LMMAXC=WDES%LMMAX(NT)
       DO NI=1,WDES%NITYP(NT)
          DO NP=1,S%NREP  ! loop over replicated ions
             ! loop over all k-points in the IRZ
             ! VKPT is in reciprocal coordinates of original cell, POSION in supercell coordinates
             DO NK=1,SIZE(GU)
#ifndef RPAgamma
                ! coefficients are real valued, and k-points are reduced
                GU_PROJ_R(LMBASE_DEST+1:LMBASE_DEST+LMMAXC, 1)= &
                GU_PROJ_R(LMBASE_DEST+1:LMBASE_DEST+LMMAXC, 1)+ &
                     REAL(EXP(CITPI*(WDES%VKPT(1,NK)*S%POSION(1,NIS_DEST) &
                            +WDES%VKPT(2,NK)*S%POSION(2,NIS_DEST) &
                            +WDES%VKPT(3,NK)*S%POSION(3,NIS_DEST))) &
                    *GU(NK)%PROJ_R(LMBASE+1:LMBASE+LMMAXC, NRP),q)*SCALE*S%KWEIGHT(NK)
#else
                GU_PROJ_R(LMBASE_DEST+1:LMBASE_DEST+LMMAXC, 1)= &
                GU_PROJ_R(LMBASE_DEST+1:LMBASE_DEST+LMMAXC, 1)+ &
                     EXP(CITPI*(WDES%VKPT(1,NK)*S%POSION(1,NIS_DEST) &
                            +WDES%VKPT(2,NK)*S%POSION(2,NIS_DEST) &
                            +WDES%VKPT(3,NK)*S%POSION(3,NIS_DEST))) &
                    *GU(NK)%PROJ_R(LMBASE+1:LMBASE+LMMAXC, NRP)*SCALE
#endif
             ENDDO
             LMBASE_DEST=LMBASE_DEST+LMMAXC
             NIS_DEST=NIS_DEST+1
          ENDDO
          LMBASE =LMBASE+LMMAXC
          NIS=NIS+1
       ENDDO
    ENDDO

  END SUBROUTINE COLLECT_PROJ_R_SUPER

!***********************************************************************
!
! build supercell "orbital" characters from the primitive cell entries
! 
!  G(g+k, N beta)     =  G_k(g, N beta) 
!  G(N'' alpha, N beta) =  sum exp^(i k R_N'' ) G_k(N' alpha, N beta)
!
! uses the previously set up G%G_PROJ and G%PROJ_PROJ
! as determined by the routine CALCULATE_G_RECIPROCAL
!
! as input the routines expect the currently considered storage positions
! LMBASE_+1: LMBASE_+ LMMAXC_ and the corresponding position POSION 
! 
! important note: to calculate the supercell characters a phase
! factor e^-ik r_N needs to be added in for the ion at position N,
! because VASP stores the orbital characters without this phase factor
! (see nonl.F). When going to the supercell quantities  (at Gamma)
! this Bloch phase factor is missing.
!
!***********************************************************************

  SUBROUTINE COLLECT_D_PROJECTORS_SUPER(S, WGW, D_GK_PROJ, P_G_PROJ, LMBASE_, LMMAXC_, POSION, LPHASE)
    USE constant
    TYPE (supercell), POINTER, INTENT(IN) :: S
    TYPE (wavedes)                        :: WGW
    GDEF, INTENT(OUT) :: D_GK_PROJ(:,:)  
    !GDEF, INTENT(OUT) :: D_PROJ_PROJ(:,:) 
    GDEF, INTENT(INOUT) :: P_G_PROJ(:,:,:)    ! Green function for each k-points
    INTEGER :: LMBASE_                     !
    INTEGER :: LMMAXC_
    REAL(q) :: POSION(3)
    LOGICAL, OPTIONAL :: LPHASE
  ! local 
    INTEGER :: NK, N, NT, NI, NIS, NIS_DEST, LMBASE, LMMAXC, LMBASE_DEST, NP
    COMPLEX(q) :: CPHASE
    LOGICAL :: LCHECK(S%WDES1%NGVECTOR)

    LCHECK=.FALSE.

#ifndef RPAgamma
    IF (SIZE(P_G_PROJ, 3) /= S%NREP) THEN
       CALL vtutor%bug("internal error in COLLECT_D_PROJECTORS_SUPER: mismatched G and k-points " &
          // str(SIZE(P_G_PROJ,3)) // " " // str(S%NREP), __FILE__, __LINE__)
    ENDIF
#endif

    D_GK_PROJ=0
    !GO_PROJ_PROJ=0

    DO NK=1,SIZE(P_G_PROJ, 3)
       IF (WGW%NGVECTOR(NK)/=S%NGVECTOR_RES(NK)) THEN
          CALL vtutor%bug("internal error in COLLECT_D_PROJECTORS_SUPER: NGVECTOR not correct " // &
             str(NK) // " " // str(WGW%NGVECTOR(NK)) // " " // str(S%NGVECTOR_RES(NK)), __FILE__, __LINE__)
       ENDIF

       IF (S%NGVECTOR_RES(NK)> SIZE(P_G_PROJ(:,:,NK),1)) THEN
          WRITE(0,*) 'internal error in COLLECT_D_PROJECTORS_SUPER: index in G_PROJ exceeded:', & 
               S%NGVECTOR_RES(NK), SIZE(P_G_PROJ(:,:,NK),1)
       ENDIF
       ! in VASP the phase factor e^i k r_N is missing (see also nonl.F)
       ! here we need to include it to get the real Bloch Green's function
       ! the second index in G is conjugated
! CHECKED AGAINST APPLY_PHASE_NLM  SAME RESULTS
       CPHASE=EXP(CITPI*(-WGW%VKPT(1,NK)*POSION(1) &
                         -WGW%VKPT(2,NK)*POSION(2) &
                         -WGW%VKPT(3,NK)*POSION(3)))
       IF (PRESENT(LPHASE)) THEN
         IF (.NOT.(LPHASE)) CPHASE=1._q
       ENDIF
       !WRITE (*,*) 'phase in COLLECT_D_PROJECTORS_SUPER'
       !WRITE (*,*) SUM(POSION), NK, CPHASE

       DO N=1,S%NGVECTOR_RES(NK)
#ifdef RPAgamma
          ! in half grid mode the coefficients on reduced grid need to be multiplied with FFTSCA(N,1)
          IF (S%INDEX_RES(N,NK)/=0) THEN
             D_GK_PROJ(S%INDEX_RES(N,NK),1:LMMAXC_)=P_G_PROJ(N, LMBASE_+1:LMBASE_+LMMAXC_,NK)*CPHASE*S%GRID_RES%FFTSCA(S%INDEX_RES(N,NK),1)
          ELSE IF (S%INDEX_RES_INV(N,NK)/=0) THEN
             D_GK_PROJ(S%INDEX_RES_INV(N,NK),1:LMMAXC_)=CONJG(P_G_PROJ(N, LMBASE_+1:LMBASE_+LMMAXC_,NK)*CPHASE)*S%GRID_RES%FFTSCA(S%INDEX_RES_INV(N,NK),1)
          ENDIF
#else
          D_GK_PROJ(S%INDEX_RES(N,NK),1:LMMAXC_)=P_G_PROJ(N, LMBASE_+1:LMBASE_+LMMAXC_,NK)*CPHASE
#endif

       ENDDO
       !WRITE (*,*) SUM(D_GK_PROJ), SUM(P_G_PROJ(:,LMBASE_+1:LMBASE_+LMMAXC_,1:NK))
    ENDDO

    NIS=1         ! loop over all ions
    LMBASE=0      ! index into PROJ_R
    NIS_DEST=1    ! ion in the supercell
    LMBASE_DEST=0 ! index into GO_PROJ_R

  END SUBROUTINE COLLECT_D_PROJECTORS_SUPER

!***********************************************************************
!
! transform projector part in supercell by summing over all 
! contributions with an appropriate phase factor
!
! as input the routines expect the currently considered storage positions
! LMBASE_+1: LMBASE_+ LMMAXC_ and the corresponding position POSION 
! 
! important note: to calculate the supercell characters a phase
! factor e^-ik r_N needs to be added in for the ion at position N,
! because VASP stores the orbital characters without this phase factor
! (see nonl.F). When going to the supercell quantities  (at Gamma)
! this Bloch phase factor is missing.
!
!***********************************************************************

  SUBROUTINE DISTRIBUTE_PROJ_R_SUPER(S, WDES, W, SIGMA, NRP, NQMAX, NQ_INDEX)
    USE constant
    TYPE (supercell), POINTER, INTENT(IN) :: S
    TYPE (wavedes)                        :: WDES
    TYPE (wavefun1), INTENT(IN) :: W    ! Green function for each k-point
    TYPE (greensf), INTENT(OUT) :: SIGMA(:) 
    INTEGER :: NRP
    INTEGER :: NQMAX
    INTEGER :: NQ_INDEX(:)
  ! local 
    INTEGER :: N, NT, NI, NIS, NIS_DEST, LMBASE, LMMAXC, LMBASE_DEST, NP, NQ
    REAL(q) :: SCALE
    COMPLEX(q) :: CPHASE

    IF (SIZE(SIGMA) /= NQMAX) THEN
       CALL vtutor%bug("internal error in DISTRIBUTE_PROJ_R_SUPER: mismatched G and k-points " // &
          str(SIZE(SIGMA)) // " " // str(NQMAX) // " " // str(S%NREP), __FILE__, __LINE__)
    ENDIF
   
    !this is the array that we want to fill
    !the index for sigma is for k-points
    DO NQ=1, NQMAX !S%NREP
       SIGMA(NQ)%PROJ_R(:,NRP)=0._q
    ENDDO

    !we have W%CPROJ(1:NPRO) to distribute into the array

    NIS=1         ! loop over all ions
    LMBASE=0      ! index into PROJ_R
    NIS_DEST=1    ! ion in the supercell
    LMBASE_DEST=0 ! index into GO_PROJ_R

    ! scaling of stored CPROJ is somewhat odd and involves a 1/sqrt(volume) term
    SCALE= SQRT(1.0_q/S%NREP)

    DO NT=1,WDES%NTYP
       LMMAXC=WDES%LMMAX(NT)
       DO NI=1,WDES%NITYP(NT)
          DO NP=1,S%NREP   ! loop over replicated ions
             DO NQ=1,NQMAX ! loop over k-points in the IRZ
                ! VKPT is in reciprocal coordinates of original cell, POSION in supercell coordinates
                !added a minus sign since backward transformation
                !LMBASE <-> LMBASE_DEST
                !NIS or NIS_DEST for the ionic position?
                !NIS_DEST since that corresponds to replicated ions (sum over displacements)
                SIGMA(NQ)%PROJ_R(LMBASE+1:LMBASE+LMMAXC, NRP)= &
                SIGMA(NQ)%PROJ_R(LMBASE+1:LMBASE+LMMAXC, NRP)+ &
!gKnew TODO: maybe force W%CPROJ to be real valued
                     EXP(-CITPI*(WDES%VKPT(1,NQ)*S%POSION(1,NIS_DEST) &
                            +WDES%VKPT(2,NQ)*S%POSION(2,NIS_DEST) &
                            +WDES%VKPT(3,NQ)*S%POSION(3,NIS_DEST))) &
                    *W%CPROJ(LMBASE_DEST+1:LMBASE_DEST+LMMAXC)*SCALE

             ENDDO
             LMBASE_DEST=LMBASE_DEST+LMMAXC
             NIS_DEST=NIS_DEST+1
          ENDDO
          LMBASE =LMBASE+LMMAXC
          NIS=NIS+1
       ENDDO
    ENDDO

  END SUBROUTINE DISTRIBUTE_PROJ_R_SUPER

!***********************************************************************
!
! this subroutine takes out the phase factor e^ik R_N from the
! supercell orbital characters to yield the original
! VASP convention (very similar to APPLY_PHASE_ONE_CENTER in fast_aug.F
!
!***********************************************************************

  SUBROUTINE APPLY_PHASE_NLM(AUG_DES, FAST_AUG, CPROJ, VKPT, NQ)
     USE constant
     TYPE (wavedes1)  :: AUG_DES
     TYPE (nonlr_struct) FAST_AUG
     GDEF             :: CPROJ(:)
     REAL(q)          :: VKPT(3)
     INTEGER          :: NQ
  ! local
     INTEGER LMMAXC, NIS, NI, NIP, NT, NPRO
     COMPLEX(q) :: CPHASE
    
     PROFILING_START( 'apply_phase_nlm' )

     NIS=1
     NPRO=0
     DO NT=1,AUG_DES%NTYP
        LMMAXC=AUG_DES%LMMAX(NT)
        IF (LMMAXC>0) THEN
           DO NI=NIS,AUG_DES%NITYP(NT)+NIS-1
              ! stupid need the global ion index to properly index FAST_AUG%POSION
              NIP=NI_GLOBAL(NI,AUG_DES%COMM_INB)  ! global index of ion 
                                                  ! since  COMM_INB contains only a single core NIP=NI
              CPHASE=EXP(CITPI*(-VKPT(1)*FAST_AUG%POSION(1,NIP) &
                                -VKPT(2)*FAST_AUG%POSION(2,NIP) &
                                -VKPT(3)*FAST_AUG%POSION(3,NIP)))
              CPROJ(NPRO+1:NPRO+LMMAXC)=CPROJ(NPRO+1:NPRO+LMMAXC)*CPHASE
              NPRO= LMMAXC+NPRO
           ENDDO
        ENDIF
        NIS = NIS+AUG_DES%NITYP(NT)
     ENDDO

     PROFILING_STOP( 'apply_phase_nlm' )

   END SUBROUTINE APPLY_PHASE_NLM

!***********************************************************************
!
! this subroutine takes out the phase factor e^ik R_N from the
! supercell orbital characters to yield the original
! VASP convention (very similar to APPLY_PHASE_ONE_CENTER in fast_aug.F
! G_PROJ can be either G_PROJ or PROJ_PROJ
!
!***********************************************************************

  SUBROUTINE APPLY_PHASE_ALPHA(WDES, GDES, G_PROJ, VKPT, NQ)
     USE constant
     TYPE (wavedes)  :: WDES
     TYPE (greensfdes):: GDES
     TYPE (nonlr_struct) FAST_AUG
     GDEF             :: G_PROJ(:,:)
     REAL(q)          :: VKPT(:)
     INTEGER          :: NQ
  ! local
     INTEGER LMMAXC, NIS, NI, NIP, NT, NPRO
     COMPLEX(q) :: CPHASE

     NIS=1
     NPRO=0
     DO NT=1,GDES%NTYP
        LMMAXC=GDES%LMMAX(NT)
        IF (LMMAXC>0) THEN
           DO NI=NIS,GDES%NITYP(NT)+NIS-1
              NIP=NI
              CPHASE=EXP(CITPI*(VKPT(1)*GDES%POSION(1,NIP)+ &
                                VKPT(2)*GDES%POSION(2,NIP)+ &
                                VKPT(3)*GDES%POSION(3,NIP)))
              !WRITE (*,*) 'NT, NI, NIP, CPHASE, POSION', NT, NI, NIP, CPHASE, GDES%POSION(:,:)
              !WRITE (*,*) 'phase in APPLY_PHASE_ALPHA', NPRO, WDES%LMMAX(NT)
              !WRITE (*,*) SUM(GDES%POSION(:,NIP)), NQ, CPHASE
              G_PROJ(:,NPRO+1:NPRO+LMMAXC)=G_PROJ(:,NPRO+1:NPRO+LMMAXC)*CPHASE
              NPRO= LMMAXC+NPRO
           ENDDO
        ENDIF
        NIS = NIS+GDES%NITYP(NT)
     ENDDO
   END SUBROUTINE APPLY_PHASE_ALPHA

!***********************************************************************
!
! extract primitive cell polarizability from supercell polarizability
! 
!  P_k(g, X)     =  P(g+k, X)
!
! index X=r  or X=N alpha
!
! first version: index X=r
!
!***********************************************************************

  SUBROUTINE DISTRIBUTE_RESPONSE_SUPER(S, WGW, P_TMP, P_G_R, NRP, NQMAX, NQ_INDEX, SCALE )

    TYPE (supercell), POINTER, INTENT(IN) :: S
    TYPE (wavedes) :: WGW          ! descriptor for response function in IRZ
    INTEGER :: NQMAX               ! total number of q-points that are considered
    INTEGER :: NQ_INDEX(:)         ! index of the corresponding q-point in the IRZ (used when NKRED is set)
    INTEGER :: NRP                 ! where to place data in P_G_R
    GDEF :: P_G_R(:,:,:)           ! destination array (response function)
    GDEF :: P_TMP(:)               ! polarizability in supercell
    REAL(q) :: SCALE               ! scaling factor
  ! local
    INTEGER :: NQ_COUNTER, NQ, N 

    P_G_R(:,NRP,:)=0

    DO NQ_COUNTER=1,NQMAX
       NQ=NQ_INDEX(NQ_COUNTER)
       IF (WGW%NGVECTOR(NQ) /=S%NGVECTOR_RES(NQ)) THEN
          CALL vtutor%bug("internal error in DISTRIBUTE_RESPONSE_SUPER: NGVECTOR not correct " // &
             str(NQ) // " " // str(WGW%NGVECTOR(NQ)) // " " // str(S%NGVECTOR_RES(NQ)), __FILE__, __LINE__)
       ENDIF
!DIR$ IVDEP
!OCL NOVREC
       DO N=1,S%NGVECTOR_RES(NQ)
#ifdef RPAgamma
          ! in half grid mode the coefficients on reduced grid are multiplied with FFTSCA(N,1)
          IF (S%INDEX_RES(N,NQ)/=0) THEN
             ! undo this factor
             P_G_R(N, NRP, NQ_COUNTER)=P_TMP(S%INDEX_RES(N,NQ))*SCALE*S%GRID_RES%FFTSCA(S%INDEX_RES(N,NQ),2)
          ENDIF
          IF (S%INDEX_RES_INV(N,NQ)/=0) THEN
             ! undo this factor
             P_G_R(N, NRP, NQ_COUNTER)=GCONJG(P_TMP(S%INDEX_RES_INV(N,NQ)))*SCALE*S%GRID_RES%FFTSCA(S%INDEX_RES_INV(N,NQ),2)
          ENDIF
#else
          P_G_R(N, NRP, NQ_COUNTER)=P_TMP(S%INDEX_RES(N,NQ))*SCALE
#endif
       ENDDO
    ENDDO
    
  END SUBROUTINE DISTRIBUTE_RESPONSE_SUPER

!
! second version: X=N alpha
!

  SUBROUTINE DISTRIBUTE_RESPONSE_PROJ_SUPER(S, WGW, P_TMP, P_G_PROJ, NAUG , NQMAX, NQ_INDEX, SCALE )

    TYPE (supercell), POINTER, INTENT(IN) :: S    
    TYPE (wavedes) :: WGW          ! descriptor for response function in IRZ
    INTEGER :: NQMAX               ! total number of q-points that are considered
    INTEGER :: NQ_INDEX(:)         ! index of the corresponding q-point in the IRZ (used when NKRED is set)
    INTEGER :: NAUG                ! where to place the data in P_G_PROJ
    GDEF :: P_G_PROJ(:,:,:)        ! destination array (response function)
    GDEF :: P_TMP(:)               ! polarizability in supercell
    REAL(q) :: SCALE               ! scaling factor
  ! local
    INTEGER :: NQ_COUNTER, NQ, N

    PROFILING_START('distribute_projectors_super')

    P_G_PROJ(:,NAUG,:)=0

    DO NQ_COUNTER=1,NQMAX
       NQ=NQ_INDEX(NQ_COUNTER)
       IF (WGW%NGVECTOR(NQ) /=S%NGVECTOR_RES(NQ)) THEN
          CALL vtutor%bug("internal error in DISTRIBUTE_RESPONSE_PROJ_SUPER: NGVECTOR not correct " &
             // str(NQ) // " " // str(WGW%NGVECTOR(NQ)) // " " // str(S%NGVECTOR_RES(NQ)), __FILE__, __LINE__)
       ENDIF
!DIR$ IVDEP
!OCL NOVREC
       DO N=1,S%NGVECTOR_RES(NQ)
#ifdef RPAgamma
          ! in half grid mode the coefficients on reduced grid are multiplied with FFTSCA(N,1)
          IF (S%INDEX_RES(N,NQ)/=0) THEN
             ! undo this factor
             P_G_PROJ(N, NAUG, NQ_COUNTER)=P_TMP(S%INDEX_RES(N,NQ))*SCALE*S%GRID_RES%FFTSCA(S%INDEX_RES(N,NQ),2)
          ENDIF
          IF (S%INDEX_RES_INV(N,NQ)/=0) THEN
             ! undo this factor
             P_G_PROJ(N, NAUG, NQ_COUNTER)=GCONJG(P_TMP(S%INDEX_RES_INV(N,NQ)))*SCALE*S%GRID_RES%FFTSCA(S%INDEX_RES_INV(N,NQ),2)
          ENDIF
#else
          P_G_PROJ(N, NAUG, NQ_COUNTER)=P_TMP(S%INDEX_RES(N,NQ))*SCALE
#endif
       ENDDO
    ENDDO

    PROFILING_STOP('distribute_projectors_super')

  END SUBROUTINE DISTRIBUTE_RESPONSE_PROJ_SUPER


!***********************************************************************
!
! extract primitive cell self-energy from supercell self-energy
! 
!  Sigma_k(g, X)     =  Sigma(g+k, X)
!
! index X=r  or X=N alpha
!
! first version: index X=r
!
!***********************************************************************

  SUBROUTINE DISTRIBUTE_SIGMA_SUPER(S, WDES, P_TMP, P_G_R, NRP, NQMAX, SCALE )

    TYPE (supercell), POINTER, INTENT(IN) :: S
    TYPE (wavedes) :: WDES         ! descriptor for orbitals function in IRZ
    INTEGER :: NQMAX               ! total number of q-points that are considered
    INTEGER :: NRP                 ! where to place data in P_G_R
    TYPE (greensf) :: P_G_R(:)     ! destination array (self-energy)
    GDEF :: P_TMP(:)               ! self-energy in supercell
    REAL(q) :: SCALE               ! scaling factor
  ! local
    INTEGER :: NQ, N, NGVECTOR
    LOGICAL :: LUSEINV

    DO NQ=1,NQMAX
       NGVECTOR=S%NGVECTOR(NQ)
       LUSEINV=.FALSE.

       IF (ASSOCIATED(S%LUSEINV)) THEN
          IF (S%LUSEINV(NQ)) THEN
             NGVECTOR=S%NGVECTOR_INV(NQ)
             LUSEINV=.TRUE.
          ENDIF
       ENDIF

       IF (.NOT. LUSEINV .AND. WDES%NGVECTOR(NQ) /=S%NGVECTOR(NQ)) THEN
          CALL vtutor%bug("internal error in DISTRIBUTE_SIGMA_SUPER: NGVECTOR not correct " // &
             str(NQ) // " " // str(WDES%NGVECTOR(NQ)) // " " // str(S%NGVECTOR(NQ)), __FILE__, __LINE__)
       ENDIF

       P_G_R(NQ)%G_R(:,NRP)=zero

#ifdef RPAgamma
       IF (.NOT. LUSEINV) THEN
!DIR$ IVDEP
!OCL NOVREC
          ! coefficients are only available for half the G+k vectors (after real to complex FFT)
          ! now all coefficients are restored applying inversion symmetry
          DO N=1,NGVECTOR
             ! in half grid mode the coefficients on reduced grid are multiplied with FFTSCA(N,1)
             IF (S%INDEX(N,NQ)/=0) THEN
                ! undo this factor
                P_G_R(NQ)%G_R(N, NRP)=P_TMP(S%INDEX(N,NQ))*SCALE*S%GRID%FFTSCA(S%INDEX(N,NQ),2)
             ENDIF
             IF (S%INDEX_INV(N,NQ)/=0) THEN
                ! undo this factor
                P_G_R(NQ)%G_R(N, NRP)=GCONJG(P_TMP(S%INDEX_INV(N,NQ)))*SCALE*S%GRID%FFTSCA(S%INDEX_INV(N,NQ),2)
             ENDIF
          ENDDO
       ELSE
!DIR$ IVDEP
!OCL NOVREC
          ! store only half of the coefficients (applying inversion symmetry in stored G_R)
          DO N=1,NGVECTOR
             ! in half grid mode the coefficients on reduced grid are multiplied with FFTSCA(N,1)
             IF (S%IND_IF_INV(N,NQ)/=0) THEN
                ! undo this factor
                P_G_R(NQ)%G_R(N, NRP)=P_TMP(S%IND_IF_INV(N,NQ))*SCALE*S%GRID%FFTSCA(S%IND_IF_INV(N,NQ),2)
             ENDIF
          ENDDO
       ENDIF
#else
!DIR$ IVDEP
!OCL NOVREC
       DO N=1,S%NGVECTOR(NQ)
          P_G_R(NQ)%G_R(N, NRP)=P_TMP(S%INDEX(N,NQ))*SCALE
       ENDDO
#endif
       ENDDO

  END SUBROUTINE DISTRIBUTE_SIGMA_SUPER

!
! second version: X=N alpha
!

  SUBROUTINE DISTRIBUTE_SIGMA_PROJ_SUPER(S, WDES, SIG_WORK, SIG_PROJ_PROJ, SIGMA, NAUG, NQMAX, SCALE )

    USE constant
    TYPE (supercell), POINTER, INTENT(IN) :: S    
    TYPE (wavedes) :: WDES         ! descriptor for GF in IRZ
    INTEGER :: NQMAX               ! total number of q-points that are considered
    INTEGER :: NAUG                ! where to place the data in P_G_PROJ
    TYPE(greensf), INTENT(INOUT) :: SIGMA(:)        ! destination array (response function)
    GDEF, INTENT(IN) :: SIG_WORK(:)               ! 
    GDEF, INTENT(IN) :: SIG_PROJ_PROJ(:)
    REAL(q) :: SCALE               ! scaling factor
  ! local
    INTEGER :: NQ, N
    INTEGER :: NIS, LMBASE, NIS_DEST, LMBASE_DEST, LMMAXC, NT, NI, NP

    PROFILING_START('distribute_sigma_projectors_super')


    DO NQ=1,NQMAX
       !WRITE (*,*) 'distribute ', NQMAX, NQ
       SIGMA(NQ)%G_PROJ(:,NAUG)=zero
    ENDDO

    DO NQ=1,NQMAX
       IF (WDES%NGVECTOR(NQ) /=S%NGVECTOR(NQ)) THEN
          CALL vtutor%bug("internal error in DISTRIBUTE_SIGMA_PROJ_SUPER: NGVECTOR not correct " // &
             str(NQ) // " " // str(WDES%NGVECTOR(NQ)) // " " // str(S%NGVECTOR(NQ)), __FILE__, __LINE__)
       ENDIF
!DIR$ IVDEP
!OCL NOVREC
       DO N=1,S%NGVECTOR(NQ)
#ifdef RPAgamma
          ! in half grid mode the coefficients on reduced grid are multiplied with FFTSCA(N,1)
          IF (S%INDEX(N,NQ)/=0) THEN
             ! undo this factor
             SIGMA(NQ)%G_PROJ(N,NAUG)=SIG_WORK(S%INDEX(N,NQ))*SCALE*S%GRID%FFTSCA(S%INDEX(N,NQ),2)
          ENDIF
          IF (S%INDEX_INV(N,NQ)/=0) THEN
             ! undo this factor
             SIGMA(NQ)%G_PROJ(N,NAUG)=GCONJG(SIG_WORK(S%INDEX_INV(N,NQ)))*SCALE*S%GRID%FFTSCA(S%INDEX_INV(N,NQ),2)
          ENDIF
#else
          ! we need here to use the index for orbitals (and Sigma), not for response
          SIGMA(NQ)%G_PROJ(N,NAUG)=SIG_WORK(S%INDEX(N,NQ))*SCALE
#endif
       ENDDO
    ENDDO

    !this is the array that we want to fill
    !the index for sigma is for k-points
    DO NQ=1,NQMAX
       SIGMA(NQ)%PROJ_PROJ(:,NAUG)=zero
    ENDDO

    !we have W%CPROJ(1:NPRO) to distribute into the array

    NIS=1         ! loop over all ions
    LMBASE=0      ! index into PROJ_R
    NIS_DEST=1    ! ion in the supercell
    LMBASE_DEST=0 ! index into GO_PROJ_R

    ! scaling of stored CPROJ is somewhat odd and involves a 1/sqrt(volume) term
    SCALE= SQRT(1.0_q/S%NREP)

    DO NT=1,WDES%NTYP
       LMMAXC=WDES%LMMAX(NT)
       DO NI=1,WDES%NITYP(NT)
          DO NP=1,S%NREP  ! loop over replicated ions
             ! loop over all k-points in the IRZ
             ! VKPT is in reciprocal coordinates of original cell, POSION in supercell coordinates
             DO NQ=1,NQMAX
                !added a minus sign since backward transformation
                !LMBASE <-> LMBASE_DEST
                !NIS or NIS_DEST for the ionic position?
                !NIS_DEST since that corresponds to replicated ions (sum over displacements)
!gKnew TODO: maybe force SIG_PROJ_PROJ to be real valued
                SIGMA(NQ)%PROJ_PROJ(LMBASE+1:LMBASE+LMMAXC, NAUG)= &
                SIGMA(NQ)%PROJ_PROJ(LMBASE+1:LMBASE+LMMAXC, NAUG)+ &
                     EXP(-CITPI*(WDES%VKPT(1,NQ)*S%POSION(1,NIS_DEST) &
                            +WDES%VKPT(2,NQ)*S%POSION(2,NIS_DEST) &
                            +WDES%VKPT(3,NQ)*S%POSION(3,NIS_DEST))) &
                    *SIG_PROJ_PROJ(LMBASE_DEST+1:LMBASE_DEST+LMMAXC)*SCALE

             ENDDO
             LMBASE_DEST=LMBASE_DEST+LMMAXC
             NIS_DEST=NIS_DEST+1
          ENDDO
          LMBASE =LMBASE+LMMAXC
          NIS=NIS+1
       ENDDO
    ENDDO

    PROFILING_STOP('distribute_sigma_projectors_super')

  END SUBROUTINE DISTRIBUTE_SIGMA_PROJ_SUPER


!***********************************************************************
!
! FFT_SIGMA fourier transform along first index from R to G
!
!    from SIGMA%G_R     ->  SIGMA%GG      is calculated by FFT
! the entriy SIGMA%GG are allocated
!
! SIGMA%G_R is deallocated
!
! scaling issue: to make this the exact reciprocal of FFT_G the 
! results need to be scaled by the number of grid points
! this might not be required to FFT the selfenergy however
! 
!***********************************************************************

  SUBROUTINE FFT_SIGMA_SUPER( WDES, SIGMA, GDES, NK, CPHASE, LPHASE )
    USE ini
    TYPE (wavedes) :: WDES
    TYPE (greensf) :: SIGMA
    TYPE (greensfdes) :: GDES
    COMPLEX(q) :: CPHASE(:)
    LOGICAL :: LPHASE
  ! local
    INTEGER :: ICOL
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1
    INTEGER :: NK
    INTEGER :: NRPLWV_COL_DATA_POINTS

    IF (ASSOCIATED(SIGMA%GG)) THEN
       CALL vtutor%bug("internal error in FFT_SIGMA_SUPER: SIGMA%GG is associated " // &
          str(ASSOCIATED(SIGMA%GG)), __FILE__, __LINE__)
    ENDIF

! Fock FFT is used
    CALL SETWDES(WDES, WDES1, NK)
    CALL NEWWAV(W1, WDES1,.TRUE.)
! handle G_R: transpose G_R -> R_G, and FFT to obtain GG
    CALL ALLOCATE_R_G(GDES, SIGMA, NK) ! allocate the SIGMA%R_G entry
    CALL TRANSPOSE_G_R(SIGMA%G_R, SIGMA%R_G, GDES, NK)
    CALL DEALLOCATE_G_R(SIGMA)  ! deallocate the SIGMA%G_R entry

    NRPLWV_COL_DATA_POINTS=GDES%NRPLWV_COL_DATA_POINTS
    IF (ASSOCIATED(GDES%LUSEINV)) THEN
       IF (GDES%LUSEINV(NK)) THEN
          NRPLWV_COL_DATA_POINTS=GDES%NRPLWV_COL_DATA_POINTS_NK(NK)
       ENDIF
    ENDIF

! allocate SIGMA%GG
    CALL ALLOCATE_GG(GDES, SIGMA)
DOESI SIGMA%GG=0

! second index is possibly only available for half the G+k vectors
    DO ICOL=1, NRPLWV_COL_DATA_POINTS
       IF (LPHASE) THEN
          W1%CR(1:W1%WDES1%GRID%RL%NP)=SIGMA%R_G(1:W1%WDES1%GRID%RL%NP, ICOL)*(1.0_q/W1%WDES1%GRID%NPLWV)*CPHASE(1:W1%WDES1%GRID%RL%NP)
       ELSE
          W1%CR(1:W1%WDES1%GRID%RL%NP)=SIGMA%R_G(1:W1%WDES1%GRID%RL%NP, ICOL)*(1.0_q/W1%WDES1%GRID%NPLWV)
       ENDIF
       ! possibly complex to real FFT
       ! nasty CR is complex valued but only the real components
       ! matter in the Gamma-point only case
       CALL FFTEXT(WDES1%NGVECTOR, WDES1%NINDPW(1), &
            W1%CR(1), &
            SIGMA%GG(1,ICOL),WDES1%GRID, .FALSE.)
    ENDDO

    CALL DEALLOCATE_R_G(SIGMA)
    CALL DELWAV(W1, .TRUE.)
    !WRITE (*,*) 'FFT_SIGMA_SUPER done'

  END SUBROUTINE FFT_SIGMA_SUPER


!***********************************************************************
!
! another routine required for first derivate (e.g. Pulay force) contributions
! related to the derivative of the non-local projectors with 
! respect to e.g. the positions
! this contribution is given by 
!
! Sum Sigma^*(alpha,g') <d p_alpha/ d R|g> G(g,g') + c.c.
! Sum Sigma^*(alpha,alpha') <d p_alpha/ d R|g> G(g,alpha') + c.c.
! 
! the contributions are calculated at individual tau points and summed
! and G is either the major (occupied) Green's function and Sigma
! the corresponding Green's function or the minor Green's function
!
! g are PW, alpha non-local projection operator indices
!
! the strategy is to calculate and store
! kappa(alpha,g) = Sum_g' Sigma(alpha,g') G*(g,g') 
!                + Sum_alpha' Sigma(alpha,alpha') G*(g,alpha')
! and contract  <d p_alpha |g> kappa(alpha,g)
! 
! TODO: 
! - right now, no optimization (e.g. BLAS) is used this
!   obviously GEMM could be utilized
! - should we take care of LDO_TAU_LOCAL
!
!***********************************************************************

  SUBROUTINE CONTRACT_SIGMA_G_TAU_DER( WDES, NK, WEIGHT, GU, SIGMA, GDES, & 
       NONL_S, NONLR_S, LATT_CUR, FRNL, LDO_TAU_LOCAL )
    USE ini
    USE dfast
    USE fock
    USE constant
    USE lattice
    USE nonl
    IMPLICIT NONE
    TYPE (wavedes) :: WDES       ! full wave descriptor
    INTEGER :: NK
    REAL(q) :: WEIGHT
    TYPE (greensf) :: GU         ! Green function for unoccupied orbitals 
    TYPE (greensf) :: SIGMA      ! Green function for occupied orbitals
    TYPE (greensfdes) :: GDES    
    TYPE (nonl_struct) NONL_S    ! non-local pseudopotential
    TYPE (nonlr_struct) NONLR_S
    TYPE (latt)    :: LATT_CUR  ! lattice vectors
    REAL (q)       :: FRNL(:,:) ! non-local contribution to forces
    LOGICAL :: LDO_TAU_LOCAL
! local
    TYPE (wavedes1):: WDES1      ! wave descriptor for current k-point
    INTEGER :: ICOL, NG, NGP, NPRO, NPROP
    GDEF, ALLOCATABLE :: PROJ_G_DER(:,:)
    COMPLEX(q), ALLOCATABLE :: CW(:,:)
    GDEF, ALLOCATABLE ::       CPROJXYZ(:,:,:)
! use blocking over 64 band blocks
    INTEGER, PARAMETER :: NBLOCK=64
    TYPE (wavedes) :: WDES_ORB
    LOGICAL, EXTERNAL :: CALCULATE_RPA_FORCES

    IF (.NOT. WDES%LOVERL .OR. .NOT. CALCULATE_RPA_FORCES()) RETURN

    CALL SETWDES(WDES, WDES1, NK)

    IF (GDES%NPRO_ROW > 0) THEN
       ! if SIGMA%PROJ_G is not yet associated allocate and calculate
       IF (.NOT. ASSOCIATED(SIGMA%PROJ_G)) THEN
          CALL ALLOCATE_PROJ_G(GDES,SIGMA)
          CALL TRANSPOSE_G_PROJ(SIGMA%G_PROJ, SIGMA%PROJ_G, GDES)
       ENDIF
    ENDIF

    ! copy wave descriptor and set number of orbitals to NPRO
    WDES_ORB=WDES
    WDES_ORB%NBANDS=((WDES1%NPRO+WDES%NB_PAR-1)/WDES%NB_PAR)
    WDES_ORB%NB_TOT=WDES_ORB%NBANDS*WDES%NB_PAR


    ! calculate in small blocks
    ALLOCATE(CW(WDES_ORB%NRPLWV, WDES_ORB%NBANDS))
    CW=(0._q,0._q)

    ALLOCATE(PROJ_G_DER(m_ WDES%NRPLWV , NBLOCK))


!   outer loop over blocks of projectors
 blk: DO NPRO=1, WDES1%NPRO, NBLOCK
    PROJ_G_DER=zero
    IF (LDO_TAU_LOCAL) THEN
!=======================================================================
! calculate kappa(alpha,g)
!=======================================================================
       !   loop over g'  Sum_g' Sigma*(alpha,g') G(g,g')
       DO ICOL=1,GDES%NRPLWV_COL_DATA_POINTS ! data distributed over g'   
          ! determine corresponding plane wave coefficient (note sine cosine transform)
          NGP=ICOL-2+ m_ GDES%NRPLWV_POS
          ! make sure that we do not index beyond last relevant storage position
          ! see notes below for m_ 
          IF (NGP > m_ WDES1%NGVECTOR) CYCLE
          ! only use actual data points
          ! in the Gamma point only mode the number of entries in the real valued
          ! Greens function is doubled (sine and cosine transformation)
          ! m_ is define to be 2 for gammareal version taking care of this issue
          DO NPROP=NPRO,MIN(NPRO+NBLOCK-1,WDES1%NPRO)
             DO NG=1,m_ WDES1%NGVECTOR
                PROJ_G_DER(NG, NPROP-NPRO+1)=PROJ_G_DER(NG, NPROP-NPRO+1)+GCONJG(SIGMA%PROJ_G(NPROP,ICOL))*GU%GG(NG,ICOL)
             ENDDO
          ENDDO
       ENDDO
       
       !   projector part Sum_alpha' Sigma(alpha,alpha') G*(g,alpha')
       IF (ASSOCIATED(SIGMA%PROJ_PROJ)) THEN  ! should be always allocated for LOVERL I guess
          DO ICOL=1, GDES%NPRO_COL               ! data distributed over alpha'
             ! only use actual data points
             DO NPROP=NPRO,MIN(NPRO+NBLOCK-1,WDES1%NPRO)
                DO NG=1,m_ WDES1%NGVECTOR
                   PROJ_G_DER(NG, NPROP-NPRO+1)=PROJ_G_DER(NG, NPROP-NPRO+1)+GCONJG(SIGMA%PROJ_PROJ(NPROP,ICOL))*GU%G_PROJ(NG,ICOL)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       ! multiply by time and k-point weight
       PROJ_G_DER=PROJ_G_DER*WEIGHT
    ENDIF

    ! sum over all nodes working at this time point AND all time points presently considered
#ifdef gammareal
    CALLMPI( M_sum_d(WDES%COMM, PROJ_G_DER, SIZE(PROJ_G_DER)))
#else
    CALLMPI( M_sum_z(WDES%COMM, PROJ_G_DER, SIZE(PROJ_G_DER)))
#endif
    ! now distribute round robin onto cores
    DO NPROP=NPRO,MIN(NPRO+NBLOCK-1,WDES1%NPRO)
       ! local then store
       IF (MOD(NPROP-1, WDES%NB_PAR)+1==WDES%NB_LOW) THEN
#ifdef gammareal
          ! go back to complex data representation for  kappa(alpha,g)
          DO  NG=1,WDES1%NGVECTOR
             CW(NG, (NPROP-1)/WDES%NB_PAR+1)=CMPLX(PROJ_G_DER((NG-1)*2+1, NPROP-NPRO+1),PROJ_G_DER(NG*2, NPROP-NPRO+1),KIND=q)
          ENDDO
#else
          CW(:, (NPROP-1)/WDES%NB_PAR+1)=PROJ_G_DER(:, NPROP-NPRO+1)
#endif
       ENDIF
    ENDDO

  ENDDO blk
    
    DEALLOCATE(PROJ_G_DER)
    ALLOCATE( CPROJXYZ( WDES_ORB%NPROD, WDES_ORB%NBANDS, 3))
!=======================================================================
! contract kappa(alpha,g)  <d p_alpha |g>
! right now this calculates all elements, specifically
!    <d p_beta |g> kappa(alpha,g)
!=======================================================================
    IF (NONLR_S%LREAL) THEN
       CALL vtutor%bug("internal error in VASP: CONTRACT_SIGMA_G_TAU_DER does not support real space " &
          // "projection", __FILE__, __LINE__)
!       CALL PHASER(WDES%GRID,LATT_CUR,NONLR_S,NK,WDES)
!       CALL RPROXYZ(WDES%GRID, NONLR_S, P, LATT_CUR, W, W%WDES, ISP, NK, CPROJXYZ)
    ELSE
       CALL PHASE(WDES,NONL_S,NK)
       CALL PROJXYZ_WA(NONL_S, WDES_ORB, CW, LATT_CUR, NK, CPROJXYZ)
    ENDIF

    CALL STORE_FORCE

    DEALLOCATE(CW, CPROJXYZ)


    RETURN

  CONTAINS
!=======================================================================
! internal subroutine to extract "diagonal" components from 
! CPROJXYZ <d p_beta |g> kappa(alpha,g)
!=======================================================================

    SUBROUTINE STORE_FORCE
      INTEGER :: LMBASE, LMMAXC, NT, NI, NIS, NIP, LM, NPROP
      COMPLEX(q) :: CFOR(3, 1:SIZE(FRNL,2))

      CFOR=0
! now collect diagonal elements and store as contribution to forces
      LMBASE =0
      NIS   =1
      DO NT=1,WDES%NTYP
         LMMAXC=WDES%LMMAX(NT)
         IF (LMMAXC==0) GOTO 600
         
         DO NI=NIS,WDES%NITYP(NT)+NIS-1
            NIP=NI_GLOBAL(NI, WDES%COMM_INB) !  local storage index
            DO LM=1,LMMAXC
               NPROP=LM+LMBASE
               ! global index is NPROP, determine local index
               IF (MOD(NPROP-1, WDES%NB_PAR)+1==WDES%NB_LOW) THEN
                  CFOR(1,NIP)=CFOR(1,NIP)-CPROJXYZ( NPROP, (NPROP-1)/WDES%NB_PAR+1,  1)
                  CFOR(2,NIP)=CFOR(2,NIP)-CPROJXYZ( NPROP, (NPROP-1)/WDES%NB_PAR+1,  2)
                  CFOR(3,NIP)=CFOR(3,NIP)-CPROJXYZ( NPROP, (NPROP-1)/WDES%NB_PAR+1,  3)
               ENDIF
            ENDDO
            LMBASE = LMMAXC+LMBASE
         ENDDO
600      NIS = NIS+WDES%NITYP(NT)
      ENDDO
      ! FRNL is summed over GDES_MAT%COMM_INTER at a later point
      ! hence we only use GDES%COMM=COMM_INTRA here
      CALLMPI(  M_sum_z(GDES%COMM, CFOR , SIZE(CFOR)))
      ! add the FRNL
      FRNL=FRNL-CFOR
    END SUBROUTINE STORE_FORCE

  END SUBROUTINE CONTRACT_SIGMA_G_TAU_DER

!***********************************************************************
!
! the following subroutines redistribute that data and 
! usually at the same time conjugate the data
!  P_G_R(g,r) = P_G_R(r,g)*
! the second routine performs the oposite operation
!  P_G_R(r,g) = P_G_R(g,r)*
!
! if an additional k-point argument is supplied the inversion symmetry
! reduced data are used if LUSEINV(NK) is .TRUE.
!   
!***********************************************************************

  SUBROUTINE TRANSPOSE_R_G (P_R_G, P_G_R, GDES, NK)
    USE dfast
    USE mpimy

    GDEF :: P_R_G(:,:)
    GDEF :: P_G_R(:,:)
   TYPE (greensfdes) :: GDES
    INTEGER, OPTIONAL :: NK
  ! local
    INTEGER :: NG, NR, NG_MAX, NG_MAX_BLOCK, NGP, NRPLWV_START, NRPLWV_MAX, NCPU
    INTEGER :: ierror
    INTEGER :: NSTRIP
    GDEF, ALLOCATABLE :: P_G_R_LOCAL(:,:,:),  P_G_R_RCV(:,:,:)
    INTEGER :: NRPLWV_ROW_DATA_POINTS, NRPLWV_COL_MAX_DATA_POINTS, NRPLWV_COL_DATA_POINTS

    PROFILING_START('transpose_r_g')

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

#ifndef MPI
    IF (SIZE(P_G_R) .NE. SIZE( P_R_G)) THEN
       WRITE(*,*) 'internal error 1 in  TRANSPOSE_R_G ',SIZE(P_G_R),SIZE( P_R_G)
    ENDIF
    

    DO NG=1,NRPLWV_ROW_DATA_POINTS
       DO NR=1, GDES%MPLWV_ROW
          P_G_R(NG, NR)=GCONJG(P_R_G(NR, NG))
       ENDDO
    ENDDO
#else
    NSTRIP=NSTRIP_STANDARD

    ! communicate up to NSTRIP rows at a time
    ALLOCATE(P_G_R_LOCAL(NRPLWV_COL_MAX_DATA_POINTS, NSTRIP, GDES%COMM%NCPU), &
             P_G_R_RCV  (NRPLWV_COL_MAX_DATA_POINTS, NSTRIP, GDES%COMM%NCPU))

    ! to avoid sending NaN
    P_G_R_LOCAL=0
    
    DO NG=1, GDES%MPLWV_COL_MAX, NSTRIP
       ! the sender transposes and collects the data
       ! loop over all receiver blocks
       NG_MAX=MIN(NG+NSTRIP-1, GDES%MPLWV_COL_MAX)

!$OMP PARALLEL DO PRIVATE(NCPU,NG_MAX_BLOCK,NGP) SHARED(NG,NG_MAX,GDES,P_G_R_LOCAL,P_R_G)
       DO NCPU=0,GDES%COMM%NCPU-1
          NG_MAX_BLOCK=MIN(NG_MAX, GDES%MPLWV_ROW-NCPU* GDES%MPLWV_COL_MAX)

          ! allocation: R_G(GDES%MPLWV_ROW, NRPLWV_COL_DATA_POINTS)
          DO NGP=NG, NG_MAX_BLOCK
             P_G_R_LOCAL(1:NRPLWV_COL_DATA_POINTS, NGP-NG+1, NCPU+1)=GCONJG(P_R_G(NGP+NCPU* GDES%MPLWV_COL_MAX,1:NRPLWV_COL_DATA_POINTS))
          ENDDO
          ! zero remaining entries, just in case
!DIR$ IVDEP
!OCL NOVREC
          DO NGP=MAX(NG, NG_MAX_BLOCK+1), NG_MAX
             P_G_R_LOCAL(1:NRPLWV_COL_DATA_POINTS, NGP-NG+1, NCPU+1)=0
          ENDDO

       ENDDO
!$OMP END PARALLEL DO

       IF (NG_MAX>=NG) THEN
#ifdef gammareal
          CALL MPI_alltoall( &
            P_G_R_LOCAL(1,1,1), NRPLWV_COL_MAX_DATA_POINTS* NSTRIP , MPI_double_precision, &
            P_G_R_RCV(1,1,1),   NRPLWV_COL_MAX_DATA_POINTS* NSTRIP , MPI_double_precision, &
            GDES%COMM%MPI_COMM, ierror )
#else
          CALL MPI_alltoall( &
            P_G_R_LOCAL(1,1,1), NRPLWV_COL_MAX_DATA_POINTS* NSTRIP, MPI_double_complex, &
            P_G_R_RCV(1,1,1),   NRPLWV_COL_MAX_DATA_POINTS* NSTRIP, MPI_double_complex, &
            GDES%COMM%MPI_COMM, ierror )
#endif
          IF ( ierror /= MPI_success ) THEN
             CALL vtutor%bug("internal error in TRANSPOSE_R_G: MPI_alltoall finished with " // &
                str(ierror), __FILE__, __LINE__)
          ENDIF
       ENDIF

       ! put data received from each CPU to proper position
!$OMP PARALLEL DO PRIVATE(NCPU,NRPLWV_START,NRPLWV_MAX,NGP) SHARED(GDES,NG,NG_MAX,P_G_R,P_G_R_RCV)
       DO NCPU=0,GDES%COMM%NCPU-1
          ! allocation: P_G_R(NRPLWV_ROW_DATA_POINTS, GDES%MPLWV_COL)

          NRPLWV_START=NRPLWV_COL_MAX_DATA_POINTS*NCPU+1
          NRPLWV_MAX  =MIN(NRPLWV_COL_MAX_DATA_POINTS*(NCPU+1),NRPLWV_ROW_DATA_POINTS)
          IF ( NRPLWV_MAX >= NRPLWV_START) THEN
             DO NGP=NG, MIN(NG_MAX, GDES%MPLWV_COL)
                P_G_R(NRPLWV_START:NRPLWV_MAX, NGP )=P_G_R_RCV(1:NRPLWV_MAX-NRPLWV_START+1, NGP-NG+1, NCPU+1)
             ENDDO
          END IF
       ENDDO
!$OMP END PARALLEL DO
    ENDDO

    DEALLOCATE(P_G_R_LOCAL, P_G_R_RCV)

#endif
    PROFILING_STOP('transpose_r_g')

  END SUBROUTINE TRANSPOSE_R_G

!
! "inverse" transposition
!

  SUBROUTINE TRANSPOSE_G_R( P_G_R, P_R_G, GDES, NK)
    USE dfast
    USE mpimy

    GDEF :: P_G_R(:,:)
    GDEF :: P_R_G(:,:)
    TYPE (greensfdes) :: GDES
    INTEGER, OPTIONAL :: NK
  ! local
    INTEGER :: NG, NR, NG_MAX, NGP, MPLWV_START, MPLWV_MAX, NCPU
    INTEGER :: ierror
    INTEGER :: NSTRIP
    GDEF, ALLOCATABLE :: P_R_G_LOCAL(:,:,:),  P_R_G_RCV(:,:,:)
    INTEGER :: NRPLWV_ROW_DATA_POINTS, NRPLWV_COL_MAX_DATA_POINTS, NRPLWV_COL_DATA_POINTS

    PROFILING_START('transpose_g_r')

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

#ifndef MPI
    IF (SIZE(P_G_R) .NE. SIZE( P_R_G)) THEN
       WRITE(*,*) 'internal error 1 in  TRANSPOSE_G_R ',SIZE(P_G_R),SIZE( P_R_G)
    ENDIF

    DO NR=1, GDES%MPLWV_COL
       DO NG=1,NRPLWV_ROW_DATA_POINTS
          P_R_G(NR, NG) = GCONJG(P_G_R(NG, NR))
       ENDDO
    ENDDO
#else
    NSTRIP=NSTRIP_STANDARD

    ! communicate up to NSTRIP rows at a time
    ALLOCATE(P_R_G_LOCAL(GDES%MPLWV_COL_MAX, NSTRIP, GDES%COMM%NCPU), &
             P_R_G_RCV  (GDES%MPLWV_COL_MAX, NSTRIP, GDES%COMM%NCPU))

    DO NG=1, NRPLWV_COL_DATA_POINTS, NSTRIP
       ! the sender transposes and collects the data
       ! loop over all receiver blocks
       NG_MAX=MIN(NG+NSTRIP-1, NRPLWV_COL_DATA_POINTS)

       ! to avoid sending NaN
       P_R_G_LOCAL=0
       P_R_G_RCV=0

       DO NCPU=0,GDES%COMM%NCPU-1
          IF (NG_MAX+NCPU* NRPLWV_COL_DATA_POINTS > NRPLWV_ROW_DATA_POINTS) THEN
             WRITE(*,*) 'internal error in  TRANSPOSE_G_R: exceed allowed index ',NCPU,NG_MAX, NCPU* NRPLWV_COL_DATA_POINTS, NRPLWV_ROW_DATA_POINTS
          ENDIF

          ! allocation: P_G_R    (NRPLWV_ROW_DATA_POINTS, GDES%MPLWV_COL)
          DO NGP=NG, NG_MAX
             P_R_G_LOCAL(1:GDES%MPLWV_COL,NGP-NG+1, NCPU+1)=GCONJG(P_G_R(NGP+NCPU* NRPLWV_COL_DATA_POINTS , 1:GDES%MPLWV_COL))
          ENDDO
       ENDDO
       IF (NG_MAX>=NG) THEN
#ifdef gammareal
          CALL MPI_alltoall( &
            P_R_G_LOCAL(1,1,1), GDES%MPLWV_COL_MAX* NSTRIP , MPI_double_precision, &
            P_R_G_RCV(1,1,1),   GDES%MPLWV_COL_MAX* NSTRIP , MPI_double_precision, &
            GDES%COMM%MPI_COMM, ierror )
#else
          CALL MPI_alltoall( &
            P_R_G_LOCAL(1,1,1), GDES%MPLWV_COL_MAX* NSTRIP , MPI_double_complex, &
            P_R_G_RCV(1,1,1),   GDES%MPLWV_COL_MAX* NSTRIP , MPI_double_complex, &
            GDES%COMM%MPI_COMM, ierror )
#endif
          IF ( ierror /= MPI_success ) THEN
             CALL vtutor%bug("internal error in TRANSPOSE_G_R: MPI_alltoall finished with " // &
                str(ierror), __FILE__, __LINE__)
          ENDIF
       ENDIF

       ! put data received from each CPU to proper position
       DO NCPU=0,GDES%COMM%NCPU-1
          ! allocation: P_R_G    (GDES%MPLWV_ROW,  NRPLWV_COL_DATA_POINTS)
          MPLWV_START=GDES%MPLWV_COL_MAX*NCPU+1
          MPLWV_MAX  =MIN(GDES%MPLWV_COL_MAX*(NCPU+1),GDES%MPLWV_ROW)
          IF (MPLWV_MAX>=MPLWV_START) THEN
             DO NGP=NG, NG_MAX
                P_R_G(MPLWV_START:MPLWV_MAX, NGP )=P_R_G_RCV(1:MPLWV_MAX-MPLWV_START+1, NGP-NG+1, NCPU+1)
             ENDDO
          END IF
       ENDDO
    ENDDO

    DEALLOCATE(P_R_G_LOCAL, P_R_G_RCV)

#endif
    PROFILING_STOP('transpose_g_r')

  END SUBROUTINE TRANSPOSE_G_R

!=======================================================================
!
!  TRANSPOSE_R_PROJ
!   array is distributed over PROJ (orbital character) inititally, 
!   R is the position index
!   after the transpose PROJ is the first and fast index, and
!   distribution is now over the real space positions
!   since data distribution in PROJ is round robin (when distributed)
!   this involves pretty akward rearrangment to get back to 
!   the correct ordering
!
!=======================================================================

  SUBROUTINE TRANSPOSE_R_PROJ   (P_R_PROJ, P_PROJ_R, GDES)
    GDEF :: P_R_PROJ(:,:)
    GDEF :: P_PROJ_R(:,:)
    TYPE (greensfdes) :: GDES       
  ! local
    INTEGER :: NG, NR, NG_MAX, NG_MAX_BLOCK, NGP, NCPU, NI, LMBASE
    INTEGER :: ierror
    INTEGER :: NSTRIP
    GDEF, ALLOCATABLE :: P_PROJ_R_LOCAL(:,:,:),  P_PROJ_R_RCV(:,:,:)

#ifndef MPI
    IF (SIZE(P_R_PROJ) .NE. SIZE( P_PROJ_R)) THEN
       WRITE(*,*) 'internal error 1 in  TRANSPOSE_G_PORJ ',SIZE(P_R_PROJ),SIZE( P_PROJ_R)
    ENDIF
    DO NG=1, GDES%NPRO_ROW
       DO NR=1,GDES%MPLWV_ROW
          P_PROJ_R(NG, NR) = GCONJG(P_R_PROJ(NR, NG))
       ENDDO
    ENDDO
#else
    NSTRIP=NSTRIP_STANDARD

    ! communicate up to NSTRIP rows at a time
    ALLOCATE(P_PROJ_R_LOCAL(GDES%NPRO_COL_MAX, NSTRIP, GDES%COMM%NCPU), &
             P_PROJ_R_RCV  (GDES%NPRO_COL_MAX, NSTRIP, GDES%COMM%NCPU))

    ! to avoid sending NaN
    P_PROJ_R_LOCAL=0
    
    DO NG=1, GDES%MPLWV_COL_MAX, NSTRIP
       ! the sender transposes and collects the data
       ! loop over all receiver blocks
       NG_MAX=MIN(NG+NSTRIP-1, GDES%MPLWV_COL_MAX)

!$OMP PARALLEL DO PRIVATE(NCPU,NG_MAX_BLOCK,NGP) SHARED(NG,NG_MAX,GDES,P_PROJ_R_LOCAL,P_R_PROJ)
       DO NCPU=0,GDES%COMM%NCPU-1
          NG_MAX_BLOCK=MIN(NG_MAX, GDES%MPLWV_ROW-NCPU* GDES%MPLWV_COL_MAX)
          
          ! allocation: P_R_PROJ    (GDES%MPLWV_ROW, GDES%NPRO_COL)
          DO NGP=NG, NG_MAX_BLOCK
             P_PROJ_R_LOCAL(1:GDES%NPRO_COL, NGP-NG+1, NCPU+1)=GCONJG(P_R_PROJ(NGP+NCPU* GDES%MPLWV_COL_MAX , 1:GDES%NPRO_COL))
          ENDDO
          ! zero remaining entries, just in case
          DO NGP=MAX(NG, NG_MAX_BLOCK+1), NG_MAX
             P_PROJ_R_LOCAL(1:GDES%NPRO_COL, NGP-NG+1, NCPU+1)=0
          ENDDO

       ENDDO
!$OMP END PARALLEL DO

       IF (NG_MAX>=NG) THEN
#ifdef gammareal
          CALL MPI_alltoall( &
            !P_PROJ_R_LOCAL(1,1,1), GDES%NPRO_COL_MAX* (NG_MAX-NG+1) , MPI_double_precision, &
            !P_PROJ_R_RCV(1,1,1),   GDES%NPRO_COL_MAX* (NG_MAX-NG+1) , MPI_double_precision, &
            P_PROJ_R_LOCAL(1,1,1), GDES%NPRO_COL_MAX* NSTRIP , MPI_double_precision, &
            P_PROJ_R_RCV(1,1,1),   GDES%NPRO_COL_MAX* NSTRIP , MPI_double_precision, &
            GDES%COMM%MPI_COMM, ierror )
#else
          CALL MPI_alltoall( &
            P_PROJ_R_LOCAL(1,1,1), GDES%NPRO_COL_MAX* NSTRIP, MPI_double_complex, &
            P_PROJ_R_RCV(1,1,1),   GDES%NPRO_COL_MAX* NSTRIP, MPI_double_complex, &
            GDES%COMM%MPI_COMM, ierror )
#endif
          IF ( ierror /= MPI_success ) THEN
             CALL vtutor%bug("internal error in TRANSPOSE_R_PROJ: MPI_alltoall finished with " // &
                str(ierror), __FILE__, __LINE__)
          ENDIF
       ENDIF

       ! put data received from each CPU to proper position
!$OMP PARALLEL DO PRIVATE(NCPU,LMBASE,NI,NGP) SHARED(GDES,NG,NG_MAX,P_PROJ_R,P_PROJ_R_RCV)
       DO NCPU=1,GDES%COMM%NCPU
          LMBASE=1
          DO NI=1,GDES%NIONS_CORE(NCPU)
!          WRITE(*,*) GDES%COMM%NODE_ME, NCPU,NI, SIZE(P_PROJ_R,1),GDES%NPRO_POS(NI,NCPU)+1,GDES%NPRO_POS(NI,NCPU)+GDES%NPRO_ENTRIES(NI,NCPU),GDES%NPRO_ENTRIES(NI,NCPU)
          ! allocation: P_PROJ_R    (GDES%NLM_ROW,    GDES%MPLWV_COL )
             DO NGP=NG, MIN(NG_MAX, GDES%MPLWV_COL)
                P_PROJ_R(GDES%NPRO_POS(NI,NCPU)+1:GDES%NPRO_POS(NI,NCPU)+GDES%NPRO_ENTRIES(NI,NCPU), NGP )= & 
                     P_PROJ_R_RCV(LMBASE:LMBASE+GDES%NPRO_ENTRIES(NI,NCPU)-1, NGP-NG+1, NCPU)
             ENDDO
             LMBASE=LMBASE+GDES%NPRO_ENTRIES(NI,NCPU)
          ENDDO
       ENDDO
!$OMP END PARALLEL DO
    ENDDO
 
    DEALLOCATE(P_PROJ_R_LOCAL, P_PROJ_R_RCV)

#endif
    
  END SUBROUTINE TRANSPOSE_R_PROJ

!
! "inverse" routine
!

  SUBROUTINE TRANSPOSE_PROJ_R   (P_PROJ_R, P_R_PROJ, GDES)
    GDEF :: P_PROJ_R(:,:)
    GDEF :: P_R_PROJ(:,:)
    TYPE (greensfdes) :: GDES       
  ! local
    INTEGER :: NG, NR, NG_MAX, NG_MAX_BLOCK, NGP, NCPU, NI, LMBASE
    INTEGER :: ierror
    INTEGER :: NSTRIP
    GDEF, ALLOCATABLE :: P_PROJ_R_LOCAL(:,:,:),  P_PROJ_R_RCV(:,:,:)

#ifndef MPI
    IF (SIZE(P_R_PROJ) .NE. SIZE( P_PROJ_R)) THEN
       WRITE(*,*) 'internal error 1 in  TRANSPOSE_G_PORJ ',SIZE(P_R_PROJ),SIZE( P_PROJ_R)
    ENDIF
    DO NG=1, GDES%NPRO_ROW
       DO NR=1,GDES%MPLWV_ROW
          P_R_PROJ(NR, NG)= GCONJG(P_PROJ_R(NG, NR))
       ENDDO
    ENDDO
#else
    NSTRIP=NSTRIP_STANDARD

    ! communicate up to NSTRIP rows at a time
    ALLOCATE(P_PROJ_R_LOCAL(GDES%NPRO_COL_MAX, NSTRIP, GDES%COMM%NCPU), &
             P_PROJ_R_RCV  (GDES%NPRO_COL_MAX, NSTRIP, GDES%COMM%NCPU))

    ! to avoid sending NaN
    P_PROJ_R_LOCAL=0
    
    DO NG=1, GDES%MPLWV_COL_MAX, NSTRIP
       NG_MAX=MIN(NG+NSTRIP-1, GDES%MPLWV_COL_MAX)

       ! collect data
       DO NCPU=1,GDES%COMM%NCPU
          LMBASE=1
          DO NI=1,GDES%NIONS_CORE(NCPU)
!          WRITE(*,*) GDES%COMM%NODE_ME, NCPU,NI, SIZE(P_PROJ_R,1),GDES%NPRO_POS(NI,NCPU)+1,GDES%NPRO_POS(NI,NCPU)+GDES%NPRO_ENTRIES(NI,NCPU),GDES%NPRO_ENTRIES(NI,NCPU)
          ! allocation: P_PROJ_R    (GDES%NLM_ROW,    GDES%MPLWV_COL )
             DO NGP=NG, MIN(NG_MAX, GDES%MPLWV_COL)
                P_PROJ_R_RCV(LMBASE:LMBASE+GDES%NPRO_ENTRIES(NI,NCPU)-1, NGP-NG+1, NCPU)=P_PROJ_R(GDES%NPRO_POS(NI,NCPU)+1:GDES%NPRO_POS(NI,NCPU)+GDES%NPRO_ENTRIES(NI,NCPU), NGP )
             ENDDO
             LMBASE=LMBASE+GDES%NPRO_ENTRIES(NI,NCPU)
          ENDDO
       ENDDO

       IF (NG_MAX>=NG) THEN
#ifdef gammareal
          CALL MPI_alltoall( &
            P_PROJ_R_RCV(1,1,1), GDES%NPRO_COL_MAX* NSTRIP , MPI_double_precision, &
            P_PROJ_R_LOCAL(1,1,1), GDES%NPRO_COL_MAX* NSTRIP , MPI_double_precision, &
            GDES%COMM%MPI_COMM, ierror )
#else
          CALL MPI_alltoall( &
            P_PROJ_R_RCV(1,1,1), GDES%NPRO_COL_MAX* NSTRIP, MPI_double_complex, &
            P_PROJ_R_LOCAL(1,1,1), GDES%NPRO_COL_MAX* NSTRIP, MPI_double_complex, &
            GDES%COMM%MPI_COMM, ierror )
#endif
          IF ( ierror /= MPI_success ) THEN
             CALL vtutor%bug("internal error in TRANSPOSE_R_PROJ: MPI_alltoall finished with " // &
                str(ierror), __FILE__, __LINE__)
          ENDIF
       ENDIF

       ! the receiver  transposes and collects the data
       ! loop over all receiver blocks
       DO NCPU=0,GDES%COMM%NCPU-1
          NG_MAX_BLOCK=MIN(NG_MAX, GDES%MPLWV_ROW-NCPU* GDES%MPLWV_COL_MAX)
          
          ! allocation: P_R_PROJ    (GDES%MPLWV_ROW, GDES%NPRO_COL)
          DO NGP=NG, NG_MAX_BLOCK
             P_R_PROJ(NGP+NCPU* GDES%MPLWV_COL_MAX , 1:GDES%NPRO_COL)=GCONJG(P_PROJ_R_LOCAL(1:GDES%NPRO_COL, NGP-NG+1, NCPU+1))
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(P_PROJ_R_LOCAL, P_PROJ_R_RCV)
#endif
    
  END SUBROUTINE TRANSPOSE_PROJ_R


!=======================================================================
!
!  TRANSPOSE_G_R_RESPONSE
!   array is distributed over R inititally, G is the plane
!   wave index = row index (fast index in Fortran)
!   after the transpose R is the first and fast index, and
!   distribution is now over the plane wave coefficients
!
! this routine is equivalent to TRANSPOSE_G_R
! however, since the number of plane waves in the response function
! can be different the following "replacements" have been done
! GDES%NRPLWV_ -> GDES%RES_NRPLWV_
!
!=======================================================================

  SUBROUTINE TRANSPOSE_G_R_RESPONSE   (P_G_R   , P_R_G, GDES)
    USE dfast
    USE mpimy

    GDEF :: P_G_R(:,:)
    GDEF :: P_R_G(:,:)
    TYPE (greensfdes) :: GDES        
  ! local
    INTEGER :: NG, NR, NG_MAX, NGP, MPLWV_START, MPLWV_MAX, NCPU
    INTEGER :: ierror
    INTEGER :: NSTRIP
    GDEF, ALLOCATABLE :: P_R_G_LOCAL(:,:,:),  P_R_G_RCV(:,:,:)

    PROFILING_START('transpose_g_r_response')

#ifndef MPI
    IF (SIZE(P_G_R) .NE. SIZE( P_R_G)) THEN
       WRITE(*,*) 'internal error 1 in  TRANSPOSE_G_R ',SIZE(P_G_R),SIZE( P_R_G)
    ENDIF

    DO NR=1, GDES%MPLWV_COL
       DO NG=1,GDES%RES_NRPLWV_ROW_DATA_POINTS
          P_R_G(NR, NG) = GCONJG(P_G_R(NG, NR))
       ENDDO
    ENDDO
#else
    NSTRIP=NSTRIP_STANDARD

    ! communicate up to NSTRIP rows at a time
    ALLOCATE(P_R_G_LOCAL(GDES%MPLWV_COL_MAX, NSTRIP, GDES%COMM%NCPU), &
             P_R_G_RCV  (GDES%MPLWV_COL_MAX, NSTRIP, GDES%COMM%NCPU))

    ! to avoid sending NaN
    P_R_G_LOCAL=0

    DO NG=1, GDES%RES_NRPLWV_COL_DATA_POINTS, NSTRIP
       ! the sender transposes and collects the data
       ! loop over all receiver blocks
       NG_MAX=MIN(NG+NSTRIP-1, GDES%RES_NRPLWV_COL_DATA_POINTS)

       P_R_G_LOCAL=0
       P_R_G_RCV=0

!$OMP PARALLEL DO PRIVATE(NCPU,NGP) SHARED(NG,NG_MAX,GDES,P_R_G_LOCAL,P_G_R)
       DO NCPU=0,GDES%COMM%NCPU-1
          IF (NG_MAX+NCPU* GDES%RES_NRPLWV_COL_DATA_POINTS > GDES%RES_NRPLWV_ROW_DATA_POINTS) THEN
!$OMP CRITICAL (omp_wrt_stdout)
             WRITE(*,*) 'internal error in  TRANSPOSE_G_R: exceed allowed index ',NCPU,NG_MAX, NCPU* GDES%RES_NRPLWV_COL_DATA_POINTS, GDES%RES_NRPLWV_ROW_DATA_POINTS
!$OMP END CRITICAL (omp_wrt_stdout)
          ENDIF

          ! allocation: P_G_R    (GDES%RES_NRPLWV_ROW_DATA_POINTS, GDES%MPLWV_COL)
          DO NGP=NG, NG_MAX
             P_R_G_LOCAL(1:GDES%MPLWV_COL,NGP-NG+1, NCPU+1)=GCONJG(P_G_R(NGP+NCPU* GDES%RES_NRPLWV_COL_DATA_POINTS , 1:GDES%MPLWV_COL))
          ENDDO
       ENDDO
!$OMP END PARALLEL DO
       IF (NG_MAX>=NG) THEN
#ifdef gammareal
          CALL MPI_alltoall( &
!            P_R_G_LOCAL(1,1,1), GDES%MPLWV_COL_MAX* (NG_MAX-NG+1) , MPI_double_precision, &
!            P_R_G_RCV(1,1,1),   GDES%MPLWV_COL_MAX* (NG_MAX-NG+1) , MPI_double_precision, &
            P_R_G_LOCAL(1,1,1), GDES%MPLWV_COL_MAX* NSTRIP , MPI_double_precision, &
            P_R_G_RCV(1,1,1),   GDES%MPLWV_COL_MAX* NSTRIP , MPI_double_precision, &
            GDES%COMM%MPI_COMM, ierror )
#else
          CALL MPI_alltoall( &
         !IMHO bug but doesn't seem to change anything in the old version...(weird)
         !   P_R_G_LOCAL(1,1,1), GDES%MPLWV_COL_MAX* (NG_MAX-NG+1), MPI_double_complex, &
         !   P_R_G_RCV(1,1,1),   GDES%MPLWV_COL_MAX* (NG_MAX-NG+1), MPI_double_complex, &
            P_R_G_LOCAL(1,1,1), GDES%MPLWV_COL_MAX* NSTRIP , MPI_double_complex, &
            P_R_G_RCV(1,1,1),   GDES%MPLWV_COL_MAX* NSTRIP , MPI_double_complex, &
            GDES%COMM%MPI_COMM, ierror )
#endif
          IF ( ierror /= MPI_success ) THEN
             CALL vtutor%bug("internal error in TRANSPOSE_G_R: MPI_alltoall finished with " // &
                str(ierror), __FILE__, __LINE__)
          ENDIF
       ENDIF

       ! put data received from each CPU to proper position
!$OMP PARALLEL DO PRIVATE(NCPU,MPLWV_START,MPLWV_MAX,NGP) SHARED(GDES,NG,NG_MAX,P_R_G,P_R_G_RCV)
       DO NCPU=0,GDES%COMM%NCPU-1
          ! allocation: P_R_G    (GDES%MPLWV_ROW,  GDES%RES_NRPLWV_COL_DATA_POINTS)
          MPLWV_START=GDES%MPLWV_COL_MAX*NCPU+1
          MPLWV_MAX  =MIN(GDES%MPLWV_COL_MAX*(NCPU+1),GDES%MPLWV_ROW)
          IF (MPLWV_MAX>=MPLWV_START) THEN
             DO NGP=NG, NG_MAX
                P_R_G(MPLWV_START:MPLWV_MAX, NGP )=P_R_G_RCV(1:MPLWV_MAX-MPLWV_START+1, NGP-NG+1, NCPU+1)
             ENDDO
          END IF
       ENDDO
!$OMP END PARALLEL DO
    ENDDO

    DEALLOCATE(P_R_G_LOCAL, P_R_G_RCV)

#endif
    PROFILING_STOP('transpose_g_r_response')

  END SUBROUTINE TRANSPOSE_G_R_RESPONSE

!
!  "inverse" transposition
!

  SUBROUTINE TRANSPOSE_R_G_RESPONSE   (P_R_G   , P_G_R, GDES, LCONJG)
    USE dfast
    USE mpimy

    GDEF :: P_R_G(:,:)
    GDEF :: P_G_R(:,:)
   TYPE (greensfdes) :: GDES       
    LOGICAL, OPTIONAL :: LCONJG
  ! local
    INTEGER :: NG, NR, NG_MAX, NG_MAX_BLOCK, NGP, NRPLWV_START, NRPLWV_MAX, NCPU
    INTEGER :: ierror
    INTEGER :: NSTRIP
    GDEF, ALLOCATABLE :: P_G_R_LOCAL(:,:,:),  P_G_R_RCV(:,:,:)

    PROFILING_START('transpose_r_g_response')

#ifndef MPI
    IF (SIZE(P_G_R) .NE. SIZE( P_R_G)) THEN
       WRITE(*,*) 'internal error 1 in  TRANSPOSE_R_G ',SIZE(P_G_R),SIZE( P_R_G)
    ENDIF

    DO NG=1,GDES%RES_NRPLWV_ROW_DATA_POINTS
       DO NR=1, GDES%MPLWV_ROW
          IF (PRESENT(LCONJG)) THEN
             IF (LCONJG) THEN
                P_G_R(NG, NR)=GCONJG(P_R_G(NR, NG))
             ELSE
                P_G_R(NG, NR)=P_R_G(NR, NG)
             ENDIF
          ELSE           
             P_G_R(NG, NR)=GCONJG(P_R_G(NR, NG))
          ENDIF
       ENDDO
    ENDDO
#else
    NSTRIP=NSTRIP_STANDARD

    ! communicate up to NSTRIP rows at a time
    ALLOCATE(P_G_R_LOCAL(GDES%RES_NRPLWV_COL_MAX_DATA_POINTS, NSTRIP, GDES%COMM%NCPU), &
             P_G_R_RCV  (GDES%RES_NRPLWV_COL_MAX_DATA_POINTS, NSTRIP, GDES%COMM%NCPU))

    ! to avoid sending NaN
    P_G_R_LOCAL=0
    
    DO NG=1, GDES%MPLWV_COL_MAX, NSTRIP
       ! the sender transposes and collects the data
       ! loop over all receiver blocks
       NG_MAX=MIN(NG+NSTRIP-1, GDES%MPLWV_COL_MAX)

       DO NCPU=0,GDES%COMM%NCPU-1
          NG_MAX_BLOCK=MIN(NG_MAX, GDES%MPLWV_ROW-NCPU* GDES%MPLWV_COL_MAX)

          ! allocation: R_G(GDES%MPLWV_ROW, GDES%RES_NRPLWV_COL_DATA_POINTS)
          DO NGP=NG, NG_MAX_BLOCK
             IF (PRESENT(LCONJG)) THEN
                IF (LCONJG) THEN
                   P_G_R_LOCAL(1:GDES%RES_NRPLWV_COL_DATA_POINTS, NGP-NG+1, NCPU+1)=GCONJG(P_R_G(NGP+NCPU* GDES%MPLWV_COL_MAX, 1:GDES%RES_NRPLWV_COL_DATA_POINTS))
                ELSE
                   WRITE (*,*) 'TRANSPOSE_R_G_RESPONSE not conjugating'
                   P_G_R_LOCAL(1:GDES%RES_NRPLWV_COL_DATA_POINTS, NGP-NG+1, NCPU+1)=(P_R_G(NGP+NCPU* GDES%MPLWV_COL_MAX, 1:GDES%RES_NRPLWV_COL_DATA_POINTS))
                ENDIF
             ELSE
                P_G_R_LOCAL(1:GDES%RES_NRPLWV_COL_DATA_POINTS, NGP-NG+1, NCPU+1)=GCONJG(P_R_G(NGP+NCPU* GDES%MPLWV_COL_MAX, 1:GDES%RES_NRPLWV_COL_DATA_POINTS))
             ENDIF
          ENDDO
          ! zero remaining entries, just in case
          DO NGP=MAX(NG, NG_MAX_BLOCK+1), NG_MAX
             P_G_R_LOCAL(1:GDES%RES_NRPLWV_COL_DATA_POINTS, NGP-NG+1, NCPU+1)=0
          ENDDO

       ENDDO
       IF (NG_MAX>=NG) THEN
#ifdef gammareal
          CALL MPI_alltoall( &
    !        P_G_R_LOCAL(1,1,1), GDES%RES_NRPLWV_COL_MAX_DATA_POINTS* (NG_MAX-NG+1) , MPI_double_precision, &
    !        P_G_R_RCV(1,1,1),   GDES%RES_NRPLWV_COL_MAX_DATA_POINTS* (NG_MAX-NG+1) , MPI_double_precision, &
            P_G_R_LOCAL(1,1,1), GDES%RES_NRPLWV_COL_MAX_DATA_POINTS* NSTRIP , MPI_double_precision, &
            P_G_R_RCV(1,1,1),   GDES%RES_NRPLWV_COL_MAX_DATA_POINTS* NSTRIP , MPI_double_precision, &
            GDES%COMM%MPI_COMM, ierror )
#else
          CALL MPI_alltoall( &
    !        P_G_R_LOCAL(1,1,1), GDES%RES_NRPLWV_COL_MAX_DATA_POINTS* (NG_MAX-NG+1), MPI_double_complex, &
    !        P_G_R_RCV(1,1,1),   GDES%RES_NRPLWV_COL_MAX_DATA_POINTS* (NG_MAX-NG+1), MPI_double_complex, &
            P_G_R_LOCAL(1,1,1), GDES%RES_NRPLWV_COL_MAX_DATA_POINTS* NSTRIP, MPI_double_complex, &
            P_G_R_RCV(1,1,1),   GDES%RES_NRPLWV_COL_MAX_DATA_POINTS* NSTRIP, MPI_double_complex, &
            GDES%COMM%MPI_COMM, ierror )
#endif
          IF ( ierror /= MPI_success ) THEN
             CALL vtutor%bug("internal error in TRANSPOSE_R_G: MPI_alltoall finished with " // &
                str(ierror), __FILE__, __LINE__)
          ENDIF
       ENDIF

       ! put data received from each CPU to proper position
       DO NCPU=0,GDES%COMM%NCPU-1
          ! allocation: P_G_R(GDES%RES_NRPLWV_ROW_DATA_POINTS, GDES%MPLWV_COL)

          NRPLWV_START=GDES%RES_NRPLWV_COL_MAX_DATA_POINTS*NCPU+1
          NRPLWV_MAX  =MIN(GDES%RES_NRPLWV_COL_MAX_DATA_POINTS*(NCPU+1),GDES%RES_NRPLWV_ROW_DATA_POINTS)
          IF ( NRPLWV_MAX >= NRPLWV_START) THEN
             DO NGP=NG, MIN(NG_MAX, GDES%MPLWV_COL)
                P_G_R(NRPLWV_START:NRPLWV_MAX, NGP )=P_G_R_RCV(1:NRPLWV_MAX-NRPLWV_START+1, NGP-NG+1, NCPU+1)
             ENDDO
          END IF
       ENDDO
    ENDDO

    DEALLOCATE(P_G_R_LOCAL, P_G_R_RCV)

#endif
    PROFILING_STOP('transpose_r_g_response')

  END SUBROUTINE TRANSPOSE_R_G_RESPONSE

!=======================================================================
!
!  TRANSPOSE_G_G_RESPONSE
! transpose an array stored in reciprocal space
! this version is for response like quantities
!
!=======================================================================


  SUBROUTINE TRANSPOSE_G_G_RESPONSE   (P_GG_1, P_GG_2, GDES )
    USE dfast
    USE mpimy
    GDEF :: P_GG_1(:,:)
    GDEF :: P_GG_2(:,:)
    TYPE (greensfdes) :: GDES       
  ! local
    INTEGER :: NG, NR, NG_MAX, NG_MAX_BLOCK, NGP, NRPLWV_START, NRPLWV_MAX, NCPU
    INTEGER :: ierror
    INTEGER :: NSTRIP
    GDEF, ALLOCATABLE :: P_LOCAL(:,:,:),  P_RCV(:,:,:)

    PROFILING_START('transpose_g_g_response')

!    WRITE(*,*) 'enter'
#ifndef MPI
    IF (SIZE(P_GG_2) .NE. SIZE( P_GG_1)) THEN
       WRITE(*,*) 'internal error 1 in  TRANSPOSE_G_G_RESPONSE ',SIZE(P_GG_2),SIZE( P_GG_1)
    ENDIF

    DO NG=1,GDES%RES_NRPLWV_ROW_DATA_POINTS
       DO NR=1,GDES%RES_NRPLWV_COL_DATA_POINTS
          P_GG_2(NG, NR)=GCONJG(P_GG_1(NR, NG))
       ENDDO
    ENDDO
#else
    NSTRIP=NSTRIP_STANDARD

    ! communicate up to NSTRIP rows at a time
    ALLOCATE(P_LOCAL(GDES%RES_NRPLWV_COL_MAX_DATA_POINTS, NSTRIP, GDES%COMM%NCPU), &
             P_RCV  (GDES%RES_NRPLWV_COL_MAX_DATA_POINTS, NSTRIP, GDES%COMM%NCPU))

    ! to avoid sending NaN
    P_LOCAL=0
    
    DO NG=1, GDES%RES_NRPLWV_COL_MAX_DATA_POINTS, NSTRIP
       ! the sender transposes and collects the data
       ! loop over all receiver blocks
       NG_MAX=MIN(NG+NSTRIP-1, GDES%RES_NRPLWV_COL_MAX_DATA_POINTS)
!       WRITE(*,*) NG, NG_MAX

       DO NCPU=0,GDES%COMM%NCPU-1
          NG_MAX_BLOCK=MIN(NG_MAX, GDES%RES_NRPLWV_ROW_DATA_POINTS-NCPU* GDES%RES_NRPLWV_COL_MAX_DATA_POINTS)

          ! allocation: GG(GDES%RES_NRPLWV_ROW, GDES%RES_NRPLWV_COL_DATA_POINTS)
          DO NGP=NG, NG_MAX_BLOCK
             P_LOCAL(1:GDES%RES_NRPLWV_COL_DATA_POINTS, NGP-NG+1, NCPU+1)=GCONJG(P_GG_1(NGP+NCPU* GDES%RES_NRPLWV_COL_MAX_DATA_POINTS, 1:GDES%RES_NRPLWV_COL_DATA_POINTS))
          ENDDO
          ! zero remaining entries, just in case
          DO NGP=MAX(NG, NG_MAX_BLOCK+1), NG_MAX
             P_LOCAL(1:GDES%RES_NRPLWV_COL_DATA_POINTS, NGP-NG+1, NCPU+1)=0
          ENDDO

       ENDDO

!       WRITE(*,*) 'send'
       IF (NG_MAX>=NG) THEN
#ifdef gammareal
          CALL MPI_alltoall( &
    !        P_LOCAL(1,1,1), GDES%RES_NRPLWV_COL_MAX_DATA_POINTS* (NG_MAX-NG+1) , MPI_double_precision, &
    !        P_RCV(1,1,1),   GDES%RES_NRPLWV_COL_MAX_DATA_POINTS* (NG_MAX-NG+1) , MPI_double_precision, &
            P_LOCAL(1,1,1), GDES%RES_NRPLWV_COL_MAX_DATA_POINTS* NSTRIP , MPI_double_precision, &
            P_RCV(1,1,1),   GDES%RES_NRPLWV_COL_MAX_DATA_POINTS* NSTRIP , MPI_double_precision, &
            GDES%COMM%MPI_COMM, ierror )
#else
          CALL MPI_alltoall( &
    !        P_LOCAL(1,1,1), GDES%RES_NRPLWV_COL_MAX_DATA_POINTS* (NG_MAX-NG+1), MPI_double_complex, &
    !        P_RCV(1,1,1),   GDES%RES_NRPLWV_COL_MAX_DATA_POINTS* (NG_MAX-NG+1), MPI_double_complex, &
            P_LOCAL(1,1,1), GDES%RES_NRPLWV_COL_MAX_DATA_POINTS* NSTRIP, MPI_double_complex, &
            P_RCV(1,1,1),   GDES%RES_NRPLWV_COL_MAX_DATA_POINTS* NSTRIP, MPI_double_complex, &
            GDES%COMM%MPI_COMM, ierror )
#endif
          IF ( ierror /= MPI_success ) THEN
             CALL vtutor%bug("internal error in TRANSPOSE_G_G_RESPONSE: MPI_alltoall finished with " &
                // str(ierror), __FILE__, __LINE__)
          ENDIF
       ENDIF

!       WRITE(*,*) 'store'
       ! put data received from each CPU to proper position
       DO NCPU=0,GDES%COMM%NCPU-1
          ! allocation: P_GG_2(GDES%RES_NRPLWV_ROW_DATA_POINTS, GDES%RES_NRPLWV_COL_DATA_POINTS)

          NRPLWV_START=GDES%RES_NRPLWV_COL_MAX_DATA_POINTS*NCPU+1
          NRPLWV_MAX  =MIN(GDES%RES_NRPLWV_COL_MAX_DATA_POINTS*(NCPU+1),GDES%RES_NRPLWV_ROW_DATA_POINTS)
          IF ( NRPLWV_MAX >= NRPLWV_START) THEN
             DO NGP=NG, MIN(NG_MAX, GDES%RES_NRPLWV_COL_DATA_POINTS)
                P_GG_2(NRPLWV_START:NRPLWV_MAX, NGP )=P_RCV(1:NRPLWV_MAX-NRPLWV_START+1, NGP-NG+1, NCPU+1)
             ENDDO
          END IF
!       WRITE(*,*) 'done'
       ENDDO
    ENDDO

    DEALLOCATE(P_LOCAL, P_RCV)

#endif
    PROFILING_STOP('transpose_g_g_response')

  END SUBROUTINE TRANSPOSE_G_G_RESPONSE


!=======================================================================
!
!  TRANSPOSE_G_G
! transpose an array stored in reciprocal space
! this version is for Green's function related data
! if an additional k-point index is passed to the routine
! the compressed storage mode (applying inversion symmetry) is used
! if NK1 is passed, 
!    it is assumed that the input array applies inversion
!    symmetry along first dimension, and thus the output array applies
!    inversion symmetry along the second dimension
! if NK2 is passed, 
!    it is assumed that the input array applies inversion
!    symmetry along second dimension, and thus the output array applies
!    inversion symmetry along the first dimension
!
!=======================================================================


  SUBROUTINE TRANSPOSE_G_G   (P_GG_1, P_GG_2, GDES, NK1, NK2 )
    USE dfast
    USE mpimy
    GDEF :: P_GG_1(:,:)
    GDEF :: P_GG_2(:,:)
    TYPE (greensfdes) :: GDES
    INTEGER, OPTIONAL :: NK1, NK2
  ! local
    INTEGER :: NG, NR, NG_MAX, NG_MAX_BLOCK, NGP, NRPLWV_START, NRPLWV_MAX, NCPU
    INTEGER :: ierror
    INTEGER :: NSTRIP
    GDEF, ALLOCATABLE :: P_LOCAL(:,:,:),  P_RCV(:,:,:)
    INTEGER :: NRPLWV_ROW_DATA_POINTS1, NRPLWV_COL_MAX_DATA_POINTS1, NRPLWV_COL_DATA_POINTS1
    INTEGER :: NRPLWV_ROW_DATA_POINTS2, NRPLWV_COL_MAX_DATA_POINTS2, NRPLWV_COL_DATA_POINTS2


    PROFILING_START('transpose_g_g_response')
    
    NRPLWV_ROW_DATA_POINTS1    =GDES%NRPLWV_ROW_DATA_POINTS
    NRPLWV_COL_MAX_DATA_POINTS1=GDES%NRPLWV_COL_MAX_DATA_POINTS
    NRPLWV_COL_DATA_POINTS1    =GDES%NRPLWV_COL_DATA_POINTS
    IF (PRESENT(NK1) .AND. ASSOCIATED(GDES%LUSEINV) ) THEN
       IF (GDES%LUSEINV(NK1)) THEN
          NRPLWV_ROW_DATA_POINTS1    =GDES%NRPLWV_ROW_DATA_POINTS_NK(NK1)
          NRPLWV_COL_MAX_DATA_POINTS1=GDES%NRPLWV_COL_MAX_DATA_POINTS_NK(NK1)
          NRPLWV_COL_DATA_POINTS1    =GDES%NRPLWV_COL_DATA_POINTS_NK(NK1)
       ENDIF
    ENDIF

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


!    WRITE(*,*) 'enter'
#ifndef MPI
    IF (SIZE(P_GG_2) .NE. SIZE( P_GG_1)) THEN
       WRITE(*,*) 'internal error 1 in  TRANSPOSE_G_G ',SIZE(P_GG_2),SIZE( P_GG_1)
    ENDIF

    DO NG=1,NRPLWV_ROW_DATA_POINTS2
       DO NR=1,NRPLWV_COL_DATA_POINTS1
          P_GG_2(NG, NR)=GCONJG(P_GG_1(NR, NG))
       ENDDO
    ENDDO
#else
    NSTRIP=NSTRIP_STANDARD

    ! communicate up to NSTRIP rows at a time
    ALLOCATE(P_LOCAL(NRPLWV_COL_MAX_DATA_POINTS2, NSTRIP, GDES%COMM%NCPU), &
             P_RCV  (NRPLWV_COL_MAX_DATA_POINTS2, NSTRIP, GDES%COMM%NCPU))

    ! to avoid sending NaN
    P_LOCAL=0
    
    DO NG=1, NRPLWV_COL_MAX_DATA_POINTS1, NSTRIP
       ! the sender transposes and collects the data
       ! loop over all receiver blocks
       NG_MAX=MIN(NG+NSTRIP-1, NRPLWV_COL_MAX_DATA_POINTS1)
!       WRITE(*,*) NG, NG_MAX

       DO NCPU=0,GDES%COMM%NCPU-1
          NG_MAX_BLOCK=MIN(NG_MAX, NRPLWV_ROW_DATA_POINTS1-NCPU* NRPLWV_COL_MAX_DATA_POINTS1)

          ! allocation: GG(NRPLWV_ROW1, NRPLWV_COL_DATA_POINTS2)
          DO NGP=NG, NG_MAX_BLOCK
             P_LOCAL(1:NRPLWV_COL_DATA_POINTS2, NGP-NG+1, NCPU+1)=GCONJG(P_GG_1(NGP+NCPU* NRPLWV_COL_MAX_DATA_POINTS1,1:NRPLWV_COL_DATA_POINTS2))
          ENDDO
          ! zero remaining entries, just in case
          DO NGP=MAX(NG, NG_MAX_BLOCK+1), NG_MAX
             P_LOCAL(1:NRPLWV_COL_DATA_POINTS2, NGP-NG+1, NCPU+1)=0
          ENDDO
       ENDDO

!       WRITE(*,*) 'send'
       IF (NG_MAX>=NG) THEN
#ifdef gammareal
          CALL MPI_alltoall( &
            P_LOCAL(1,1,1), NRPLWV_COL_MAX_DATA_POINTS2* NSTRIP , MPI_double_precision, &
            P_RCV(1,1,1),   NRPLWV_COL_MAX_DATA_POINTS2* NSTRIP , MPI_double_precision, &
            GDES%COMM%MPI_COMM, ierror )
#else
          CALL MPI_alltoall( &
            P_LOCAL(1,1,1), NRPLWV_COL_MAX_DATA_POINTS2* NSTRIP, MPI_double_complex, &
            P_RCV(1,1,1),   NRPLWV_COL_MAX_DATA_POINTS2* NSTRIP, MPI_double_complex, &
            GDES%COMM%MPI_COMM, ierror )
#endif
          IF ( ierror /= MPI_success ) THEN
             CALL vtutor%bug("internal error in TRANSPOSE_G_G: MPI_alltoall finished with " // &
                str(ierror), __FILE__, __LINE__)
          ENDIF
       ENDIF

!       WRITE(*,*) 'store'
       ! put data received from each CPU to proper position
       DO NCPU=0,GDES%COMM%NCPU-1
          ! allocation: P_GG_2(NRPLWV_ROW_DATA_POINTS2, NRPLWV_COL_DATA_POINTS1)

          NRPLWV_START=NRPLWV_COL_MAX_DATA_POINTS2*NCPU+1
          NRPLWV_MAX  =MIN(NRPLWV_COL_MAX_DATA_POINTS2*(NCPU+1),NRPLWV_ROW_DATA_POINTS2)
          IF ( NRPLWV_MAX >= NRPLWV_START) THEN
             DO NGP=NG, MIN(NG_MAX, NRPLWV_COL_DATA_POINTS1)
                P_GG_2(NRPLWV_START:NRPLWV_MAX, NGP )=P_RCV(1:NRPLWV_MAX-NRPLWV_START+1, NGP-NG+1, NCPU+1)
             ENDDO
          END IF
!       WRITE(*,*) 'done'
       ENDDO
    ENDDO

    DEALLOCATE(P_LOCAL, P_RCV)

#endif
    PROFILING_STOP('transpose_g_g_response')

  END SUBROUTINE TRANSPOSE_G_G



!=======================================================================
!
!  TRANSPOSE_G_PROJ_RESPONSE
!   array is distributed over PROJ (NLM) inititally, G is the plane
!   wave index = row index (fast index in Fortran)
!   after the transpose PROJ (NLM) is the first and fast index, and
!   distribution is now over the plane wave coefficients
!   since data distribution in PROJ is round robin (when distributed)
!   this involves pretty akward rearrangment to get back to 
!   the correct ordering
!   TODO: these are really routines for RESPONSE; and should
!   hence carry a _RESPONSE at the end
!
!=======================================================================

  SUBROUTINE TRANSPOSE_G_PROJ_RESPONSE(P_G_PROJ, P_PROJ_G, GDES)
    GDEF :: P_G_PROJ(:,:)
    GDEF :: P_PROJ_G(:,:)
    TYPE (greensfdes) :: GDES   
  ! local
    INTEGER :: NG, NR, NG_MAX, NGP, NCPU, NI, LMBASE
    INTEGER :: ierror
    INTEGER :: NSTRIP
    GDEF, ALLOCATABLE :: P_PROJ_G_LOCAL(:,:,:),  P_PROJ_G_RCV(:,:,:)
   
    !IF (GDES%NLM_COL_MAX.eq.0) RETURN

    PROFILING_START( 'transpose_g_proj_response' )

#ifndef MPI
    IF (SIZE(P_G_PROJ) .NE. SIZE( P_PROJ_G)) THEN
       WRITE(*,*) 'internal error 1 in  TRANSPOSE_G_PROJ_RESPONSE',SIZE(P_G_PROJ),SIZE( P_PROJ_G)
    ENDIF
    DO NR=1, GDES%NLM_COL
       DO NG=1,GDES%RES_NRPLWV_ROW_DATA_POINTS
          P_PROJ_G(NR, NG) = GCONJG(P_G_PROJ(NG, NR))
       ENDDO
    ENDDO
#else
    NSTRIP=NSTRIP_STANDARD

    ! communicate up to NSTRIP rows at a time
    ALLOCATE(P_PROJ_G_LOCAL(GDES%NLM_COL_MAX, NSTRIP, GDES%COMM%NCPU), &
             P_PROJ_G_RCV  (GDES%NLM_COL_MAX, NSTRIP, GDES%COMM%NCPU))

    ! to avoid sending NaN
    P_PROJ_G_LOCAL=0
    
    DO NG=1, GDES%RES_NRPLWV_COL_DATA_POINTS, NSTRIP
       ! the sender transposes and collects the data
       ! loop over all receiver blocks
       NG_MAX=MIN(NG+NSTRIP-1, GDES%RES_NRPLWV_COL_DATA_POINTS)

!$OMP PARALLEL DO PRIVATE(NCPU,NGP) SHARED(NG,NG_MAX,GDES,P_PROJ_G_LOCAL,P_G_PROJ)
       DO NCPU=0,GDES%COMM%NCPU-1
          IF (NG_MAX+NCPU* GDES%RES_NRPLWV_COL_DATA_POINTS > GDES%RES_NRPLWV_ROW_DATA_POINTS) THEN
!$OMP CRITICAL (omp_wrt_stdout)
             WRITE(*,*) 'internal error in  TRANSPOSE_G_PROJ_RESPONSE: exceed allowed index ',NCPU,NG_MAX, NCPU* GDES%RES_NRPLWV_COL_DATA_POINTS, GDES%RES_NRPLWV_ROW_DATA_POINTS
!$OMP END CRITICAL (omp_wrt_stdout)
          ENDIF

          ! allocation: P_G_PROJ    (GDES%RES_NRPLWV_ROW_DATA_POINTS, GDES%NLM_COL)
          DO NGP=NG, NG_MAX
             P_PROJ_G_LOCAL(1:GDES%NLM_COL, NGP-NG+1, NCPU+1)=GCONJG(P_G_PROJ(NGP+NCPU* GDES%RES_NRPLWV_COL_DATA_POINTS , 1:GDES%NLM_COL))
          ENDDO
       ENDDO
!$OMP END PARALLEL DO
       IF (NG_MAX>=NG) THEN
#ifdef gammareal
          CALL MPI_alltoall( &
            P_PROJ_G_LOCAL(1,1,1), GDES%NLM_COL_MAX* NSTRIP , MPI_double_precision, &
            P_PROJ_G_RCV(1,1,1),   GDES%NLM_COL_MAX* NSTRIP , MPI_double_precision, &
            GDES%COMM%MPI_COMM, ierror )
#else
          CALL MPI_alltoall( &
            P_PROJ_G_LOCAL(1,1,1), GDES%NLM_COL_MAX* NSTRIP, MPI_double_complex, &
            P_PROJ_G_RCV(1,1,1),   GDES%NLM_COL_MAX* NSTRIP, MPI_double_complex, &
            GDES%COMM%MPI_COMM, ierror )
#endif
          IF ( ierror /= MPI_success ) THEN
             CALL vtutor%bug("internal error in TRANSPOSE_G_PROJ_RESPONSE: MPI_alltoall finished " &
                // "with " // str(ierror), __FILE__, __LINE__)
          ENDIF
       ENDIF

       ! put data received from each CPU to proper position
!$OMP PARALLEL DO PRIVATE(NCPU,LMBASE,NI,NGP) SHARED(GDES,NG,NG_MAX,P_PROJ_G,P_PROJ_G_RCV)
       DO NCPU=1,GDES%COMM%NCPU
          LMBASE=1
          DO NI=1,GDES%NIONS_CORE(NCPU)
!          WRITE(*,*) GDES%COMM%NODE_ME, NCPU,NI, SIZE(P_PROJ_G,1),GDES%NLM_POS(NI,NCPU)+1,GDES%NLM_POS(NI,NCPU)+GDES%NLM_ENTRIES(NI,NCPU),GDES%NLM_ENTRIES(NI,NCPU)
          ! allocation: P_PROJ_G    (GDES%NLM_ROW,    GDES%RES_NRPLWV_COL_DATA_POINTS )
             DO NGP=NG, NG_MAX
                P_PROJ_G(GDES%NLM_POS(NI,NCPU)+1:GDES%NLM_POS(NI,NCPU)+GDES%NLM_ENTRIES(NI,NCPU), NGP )= & 
                     P_PROJ_G_RCV(LMBASE:LMBASE+GDES%NLM_ENTRIES(NI,NCPU)-1, NGP-NG+1, NCPU)
             ENDDO
             LMBASE=LMBASE+GDES%NLM_ENTRIES(NI,NCPU)
          ENDDO
       ENDDO
!$OMP END PARALLEL DO
    ENDDO

    DEALLOCATE(P_PROJ_G_LOCAL, P_PROJ_G_RCV)

#endif

    PROFILING_STOP( 'transpose_g_proj_response' )
    
  END SUBROUTINE TRANSPOSE_G_PROJ_RESPONSE

!
! "inverse" transposition
!

  SUBROUTINE TRANSPOSE_PROJ_G_RESPONSE(P_PROJ_G, P_G_PROJ, GDES)
    GDEF :: P_PROJ_G(:,:)
    GDEF :: P_G_PROJ(:,:)
    TYPE (greensfdes) :: GDES   
  ! local
    INTEGER :: NG, NR, NG_MAX, NGP, NCPU, NI, LMBASE
    INTEGER :: ierror
    INTEGER :: NSTRIP
    GDEF, ALLOCATABLE :: P_PROJ_G_LOCAL(:,:,:),  P_PROJ_G_RCV(:,:,:)

    PROFILING_START( 'transpose_proj_g_response' )

#ifndef MPI
    IF (SIZE(P_G_PROJ) .NE. SIZE( P_PROJ_G)) THEN
       WRITE(*,*) 'internal error 1 in  TRANSPOSE_G_PROJ_RESPONSE ',SIZE(P_G_PROJ),SIZE( P_PROJ_G)
    ENDIF
    DO NR=1, GDES%NLM_COL
       DO NG=1,GDES%RES_NRPLWV_ROW_DATA_POINTS
          P_G_PROJ(NG, NR)= GCONJG(P_PROJ_G(NR, NG))
       ENDDO
    ENDDO
#else
    NSTRIP=NSTRIP_STANDARD

    ! communicate up to NSTRIP rows at a time
    ALLOCATE(P_PROJ_G_LOCAL(GDES%NLM_COL_MAX, NSTRIP, GDES%COMM%NCPU), &
             P_PROJ_G_RCV  (GDES%NLM_COL_MAX, NSTRIP, GDES%COMM%NCPU))

    ! to avoid sending NaN
    P_PROJ_G_RCV=0
    
    DO NG=1, GDES%RES_NRPLWV_COL_DATA_POINTS, NSTRIP
       ! the sender transposes and collects the data
       ! loop over all receiver blocks
       NG_MAX=MIN(NG+NSTRIP-1, GDES%RES_NRPLWV_COL_DATA_POINTS)

       ! put data received from each CPU to proper position
       DO NCPU=1,GDES%COMM%NCPU
          LMBASE=1
          DO NI=1,GDES%NIONS_CORE(NCPU)
!          WRITE(*,*) GDES%COMM%NODE_ME, NCPU,NI, SIZE(P_PROJ_G,1),GDES%NLM_POS(NI,NCPU)+1,GDES%NLM_POS(NI,NCPU)+GDES%NLM_ENTRIES(NI,NCPU),GDES%NLM_ENTRIES(NI,NCPU)
          ! allocation: P_PROJ_G    (GDES%NLM_ROW,    GDES%RES_NRPLWV_COL_DATA_POINTS )
             DO NGP=NG, NG_MAX
                P_PROJ_G_RCV(LMBASE:LMBASE+GDES%NLM_ENTRIES(NI,NCPU)-1, NGP-NG+1, NCPU)=P_PROJ_G(GDES%NLM_POS(NI,NCPU)+1:GDES%NLM_POS(NI,NCPU)+GDES%NLM_ENTRIES(NI,NCPU), NGP )                     
             ENDDO
             LMBASE=LMBASE+GDES%NLM_ENTRIES(NI,NCPU)
          ENDDO
       ENDDO
       IF (NG_MAX>=NG) THEN
#ifdef gammareal
          CALL MPI_alltoall( &
            P_PROJ_G_RCV(1,1,1), GDES%NLM_COL_MAX* NSTRIP , MPI_double_precision, &
            P_PROJ_G_LOCAL(1,1,1),   GDES%NLM_COL_MAX* NSTRIP , MPI_double_precision, &
            GDES%COMM%MPI_COMM, ierror )
#else
          CALL MPI_alltoall( &
            P_PROJ_G_RCV(1,1,1), GDES%NLM_COL_MAX* NSTRIP, MPI_double_complex, &
            P_PROJ_G_LOCAL(1,1,1),   GDES%NLM_COL_MAX* NSTRIP, MPI_double_complex, &
            GDES%COMM%MPI_COMM, ierror )
#endif
          IF ( ierror /= MPI_success ) THEN
             CALL vtutor%bug("internal error in TRANSPOSE_G_PROJ_RESPONSE: MPI_alltoall finished " &
                // "with " // str(ierror), __FILE__, __LINE__)
          ENDIF
       ENDIF

       DO NCPU=0,GDES%COMM%NCPU-1
          IF (NG_MAX+NCPU* GDES%RES_NRPLWV_COL_DATA_POINTS > GDES%RES_NRPLWV_ROW_DATA_POINTS) THEN
             WRITE(*,*) 'internal error in  TRANSPOSE_G_PROJ_RESPONSE: exceed allowed index ',NCPU,NG_MAX, NCPU* GDES%RES_NRPLWV_COL_DATA_POINTS, GDES%RES_NRPLWV_ROW_DATA_POINTS
          ENDIF

          ! allocation: P_G_PROJ    (GDES%RES_NRPLWV_ROW_DATA_POINTS, GDES%NLM_COL)
          DO NGP=NG, NG_MAX
             P_G_PROJ(NGP+NCPU* GDES%RES_NRPLWV_COL_DATA_POINTS , 1:GDES%NLM_COL)=GCONJG(P_PROJ_G_LOCAL(1:GDES%NLM_COL, NGP-NG+1, NCPU+1))
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(P_PROJ_G_LOCAL, P_PROJ_G_RCV)

#endif
    
    PROFILING_STOP( 'transpose_proj_g_response' )

  END SUBROUTINE TRANSPOSE_PROJ_G_RESPONSE

!
! similar to routine above, however, using
! wave function character instead of augmentation channels
!

  SUBROUTINE TRANSPOSE_G_PROJ(P_G_PROJ, P_PROJ_G, GDES)
    GDEF :: P_G_PROJ(:,:)
    GDEF :: P_PROJ_G(:,:)
    TYPE (greensfdes) :: GDES   
  ! local
    INTEGER :: NG, NR, NG_MAX, NGP, NCPU, NI, LMBASE
    INTEGER :: ierror
    INTEGER :: NSTRIP
    GDEF, ALLOCATABLE :: P_PROJ_G_LOCAL(:,:,:),  P_PROJ_G_RCV(:,:,:)

#ifndef MPI
    IF (SIZE(P_G_PROJ) .NE. SIZE( P_PROJ_G)) THEN
       WRITE(*,*) 'internal error 1 in  TRANSPOSE_G_PROJ ',SIZE(P_G_PROJ),SIZE( P_PROJ_G)
    ENDIF
    DO NR=1, GDES%NPRO_COL
       DO NG=1,GDES%NRPLWV_ROW_DATA_POINTS
          P_PROJ_G(NR, NG) = GCONJG(P_G_PROJ(NG, NR))
       ENDDO
    ENDDO
#else
    NSTRIP=NSTRIP_STANDARD

    ! communicate up to NSTRIP rows at a time
    ALLOCATE(P_PROJ_G_LOCAL(GDES%NPRO_COL_MAX, NSTRIP, GDES%COMM%NCPU), &
             P_PROJ_G_RCV  (GDES%NPRO_COL_MAX, NSTRIP, GDES%COMM%NCPU))

    ! to avoid sending NaN
    P_PROJ_G_LOCAL=0
    
    DO NG=1, GDES%NRPLWV_COL_DATA_POINTS, NSTRIP
       ! the sender transposes and collects the data
       ! loop over all receiver blocks
       NG_MAX=MIN(NG+NSTRIP-1, GDES%NRPLWV_COL_DATA_POINTS)

       DO NCPU=0,GDES%COMM%NCPU-1
          IF (NG_MAX+NCPU* GDES%NRPLWV_COL_DATA_POINTS > GDES%NRPLWV_ROW_DATA_POINTS) THEN
             WRITE(*,*) 'internal error in  TRANSPOSE_G_PROJ: exceed allowed index ',NCPU,NG_MAX, NCPU* GDES%NRPLWV_COL_DATA_POINTS, GDES%NRPLWV_ROW_DATA_POINTS
          ENDIF

          ! allocation: P_G_PROJ    (GDES%NRPLWV_ROW_DATA_POINTS, GDES%NPRO_COL)
          DO NGP=NG, NG_MAX
             P_PROJ_G_LOCAL(1:GDES%NPRO_COL, NGP-NG+1, NCPU+1)=GCONJG(P_G_PROJ(NGP+NCPU* GDES%NRPLWV_COL_DATA_POINTS , 1:GDES%NPRO_COL))
          ENDDO
       ENDDO
       IF (NG_MAX>=NG) THEN
#ifdef gammareal
          CALL MPI_alltoall( &
            P_PROJ_G_LOCAL(1,1,1), GDES%NPRO_COL_MAX* NSTRIP , MPI_double_precision, &
            P_PROJ_G_RCV(1,1,1),   GDES%NPRO_COL_MAX* NSTRIP , MPI_double_precision, &
            GDES%COMM%MPI_COMM, ierror )
#else
          CALL MPI_alltoall( &
            P_PROJ_G_LOCAL(1,1,1), GDES%NPRO_COL_MAX* NSTRIP, MPI_double_complex, &
            P_PROJ_G_RCV(1,1,1),   GDES%NPRO_COL_MAX* NSTRIP, MPI_double_complex, &
            GDES%COMM%MPI_COMM, ierror )
#endif
          IF ( ierror /= MPI_success ) THEN
             CALL vtutor%bug("internal error in TRANSPOSE_G_PROJ: MPI_alltoall finished with " // &
                str(ierror), __FILE__, __LINE__)
          ENDIF
       ENDIF

       ! put data received from each CPU to proper position
       DO NCPU=1,GDES%COMM%NCPU
          LMBASE=1
          DO NI=1,GDES%NIONS_CORE(NCPU)
!          WRITE(*,*) GDES%COMM%NODE_ME, NCPU,NI, SIZE(P_PROJ_G,1),GDES%NPRO_POS(NI,NCPU)+1,GDES%NPRO_POS(NI,NCPU)+GDES%NPRO_ENTRIES(NI,NCPU),GDES%NPRO_ENTRIES(NI,NCPU)
          ! allocation: P_PROJ_G    (GDES%NPRO_ROW,    GDES%NRPLWV_COL_DATA_POINTS )
             DO NGP=NG, NG_MAX
                P_PROJ_G(GDES%NPRO_POS(NI,NCPU)+1:GDES%NPRO_POS(NI,NCPU)+GDES%NPRO_ENTRIES(NI,NCPU), NGP )= & 
                     P_PROJ_G_RCV(LMBASE:LMBASE+GDES%NPRO_ENTRIES(NI,NCPU)-1, NGP-NG+1, NCPU)
             ENDDO
             LMBASE=LMBASE+GDES%NPRO_ENTRIES(NI,NCPU)
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(P_PROJ_G_LOCAL, P_PROJ_G_RCV)

#endif
    
  END SUBROUTINE TRANSPOSE_G_PROJ

!***********************************************************************
! 
! rotate the response function
! from one q-point to another
! according the VASP convention the second index transforms like G+q
! whereas the the first one transforms -(G+q) 
! [compare ADD_XI comments]
! this routine closely follows ROTATE_WPOT but with the added
! complexity of a distributed response function
! furthermore response functions are transposed compared
! to the convention in chi.F 
! this implies that the phase factors are conjugated as compared
! to the version in mkpoints_full.F
!
!***********************************************************************

  SUBROUTINE ROTATE_RES( P_GG, GDES, NGVECTOR, &
       CPHASE, NINDPW, LINV, LSHIFT)
    USE prec
    USE mpimy
    USE mgrid
    IMPLICIT NONE
! this subroutine only works and applies to the many k-point version
    GDEF :: P_GG(:,:)
    TYPE (greensfdes) :: GDES
    INTEGER    :: NGVECTOR    ! actual number of plane wave coefficients
    COMPLEX(q) :: CPHASE(:)   ! complex phase factor
    INTEGER       NINDPW(:)   ! plane wave index upon rotation
    LOGICAL :: LINV           ! conjugate at destination k-point
    LOGICAL :: LSHIFT         ! apply complex phase factor
#ifndef gammareal
    ! local
    INTEGER M
    GDEF,ALLOCATABLE :: P_TMP(:,:)

! allocation would be most likely ok for Gamma only, but everthing else is just non sense
    ALLOCATE(P_TMP(GDES%RES_NRPLWV_ROW_DATA_POINTS, GDES%RES_NRPLWV_COL_DATA_POINTS))
    IF (SIZE(P_TMP) /= SIZE(P_GG))  THEN
       CALL vtutor%bug("internal error in VASP: ROTATE_RES size mismatch " // str(SHAPE(P_TMP)) // &
          " " // str(SHAPE(P_GG)), __FILE__, __LINE__)
    ENDIF
DOESI P_TMP=0

    IF (NGVECTOR > GDES%RES_NRPLWV_ROW) THEN
       CALL vtutor%bug("internal error in VASP: ROTATE_RES dimension too small " // str(NGVECTOR) &
          // " " // str(GDES%RES_NRPLWV_ROW), __FILE__, __LINE__)
    ENDIF

! first build the complex conjugated if required
    IF (LINV) THEN
       P_GG=CONJG(P_GG)
    ENDIF
! second apply phase factors
    IF (LSHIFT) THEN
       DO M=1,NGVECTOR
          P_GG(M,1:GDES%RES_NRPLWV_COL )=P_GG(M,1:GDES%RES_NRPLWV_COL  )*CPHASE(M)
       ENDDO
       DO M=1,GDES%RES_NRPLWV_COL
          ! cycle if actual number of plane wave coefficients is exhausted
          IF ( M+GDES%RES_NRPLWV_POS-1 > NGVECTOR) CYCLE
          P_GG(1:GDES%RES_NRPLWV_ROW,M)=P_GG(1:GDES%RES_NRPLWV_ROW,M)*CONJG(CPHASE(M+GDES%RES_NRPLWV_POS-1))
       ENDDO
    ENDIF

! third step, re-index rows
    DO M=1,NGVECTOR
       P_TMP(NINDPW(M),1:GDES%RES_NRPLWV_COL)=P_GG(M,1:GDES%RES_NRPLWV_COL)
    ENDDO
! interchange rows and columns
    CALL TRANSPOSE_G_G_RESPONSE   (P_TMP, P_GG, GDES )
! re-index columns (which are rows after the transpose)
    DO M=1,NGVECTOR
       P_TMP(NINDPW(M),1:GDES%RES_NRPLWV_COL)=P_GG(M,1:GDES%RES_NRPLWV_COL)
    ENDDO
! interchange row and columns again
    CALL TRANSPOSE_G_G_RESPONSE   (P_TMP, P_GG, GDES )

    DEALLOCATE(P_TMP)
#endif    

  END SUBROUTINE ROTATE_RES


!***********************************************************************
!
! dump storage requirements for RPAR[K] calculations 
!
!***********************************************************************
   SUBROUTINE DUMP_STORAGE_REQUIREMENTS( GDES_TAU, GDES, S2E, IGRIDS, WDES, WGW, S, IO )
      USE tutor, ONLY: argument, VTUTOR, ISERROR, ISALERT, MAXMEM_LOW
      USE minimax_struct, ONLY : imag_grid_handle
      USE main_mpi, ONLY: COMM_INTRA_NODE_WORLD
      TYPE (greensfdes)          :: GDES_TAU             ! green function in time domain  
      TYPE (greensfdes)          :: GDES                 ! green function in freq domain  
      TYPE( screened_2e_handle ) :: S2E
      TYPE( imag_grid_handle )   :: IGRIDS
      TYPE( wavedes )            :: WDES
      TYPE( wavedes )            :: WGW
      TYPE( supercell ), POINTER :: S
      TYPE( in_struct )          :: IO 
      ! local           
      TYPE(argument) :: ARG
      INTEGER        :: IVERBOSITY=3
      LOGICAL        :: LADD_TEMP_ARRAYS = .TRUE. 
      INTEGER        :: NCPU = 1
      TYPE( mem_gw_handle ) :: MEM_GG 

#ifdef MPI 
      NCPU = COMM_INTRA_NODE_WORLD%NCPU
#endif 
     
      ! initialize 
      IF ( ICHIREAL == 0 )THEN
         IF( IO%IU0>=0 ) &
            WRITE(IO%IU0,*) ' storage requirement for chi.F not implemented yet'

      ELSE

         ! note that NTAUPAR is a multiplicative factor in GDES_TAU%RES_NRPLWV_COL_DATA_POINTS 
         ! and one may write GDES_TAU%RES_NRPLWV_COL_DATA_POINTS = GDES_TAU%RES_NRPLWV_COL * ( NCPU / NTAUPAR ) 

#ifdef debug
         IF ( IO%IU6>=0 ) WRITE(IO%IU6, "( &
            &'GDES_TAU for NTAUPAR', I4, /&
             'GDES_TAU%NRPLWV_ROW_DATA_POINTS     :', I8,/&
             'GDES_TAU%NRPLWV_COL_DATA_POINTS     :', I8,/&
             'GDES_TAU%NPRO_ROW                   :', I8,/&
             'GDES_TAU%NPRO_COL                   :', I8,/&  
             'GDES_TAU%RES_NRPLWV_ROW             :', I8,/&
             'GDES_TAU%RES_NRPLWV_COL             :', I8,/&
             'GDES_TAU%MPLWV_ROW                  :', I8,/&
             'GDES_TAU%MPLWV_COL                  :', I8,/&
             'GDES_TAU%RES_NRPLWV_ROW_DATA_POINTS :', I8,/&
             'GDES_TAU%RES_NRPLWV_COL_DATA_POINTS :', I8,/&
             'GDES_TAU%NLM_ROW                    :', I8,/&
             'GDES_TAU%NLM_COL                    :', I8,/&
             'GDES    %NRPLWV_ROW_DATA_POINTS     :', I8,/&
             'GDES    %NRPLWV_COL_DATA_POINTS     :', I8,/&
             'GDES    %NPRO_ROW                   :', I8,/&
             'GDES    %NPRO_COL                   :', I8,/&  
             'GDES    %RES_NRPLWV_ROW             :', I8,/&
             'GDES    %RES_NRPLWV_COL             :', I8,/&
             'GDES    %MPLWV_ROW                  :', I8,/&
             'GDES    %MPLWV_COL                  :', I8,/&
             'GDES    %RES_NRPLWV_ROW_DATA_POINTS :', I8,/&
             'GDES    %RES_NRPLWV_COL_DATA_POINTS :', I8,/&
             'GDES    %NLM_ROW                    :', I8,/&
             'GDES    %NLM_COL                    :', I8,/&
            )") NTAUPAR , &
              GDES_TAU%NRPLWV_ROW_DATA_POINTS    ,&
              GDES_TAU%NRPLWV_COL_DATA_POINTS    ,&
              GDES_TAU%NPRO_ROW                  ,&
              GDES_TAU%NPRO_COL                  ,&  
              GDES_TAU%RES_NRPLWV_ROW            ,&
              GDES_TAU%RES_NRPLWV_COL            ,&
              GDES_TAU%MPLWV_ROW                 ,&
              GDES_TAU%MPLWV_COL                 ,&
              GDES_TAU%RES_NRPLWV_ROW_DATA_POINTS,&
              GDES_TAU%RES_NRPLWV_COL_DATA_POINTS,&
              GDES_TAU%NLM_ROW                   ,&
              GDES_TAU%NLM_COL                   ,&
              GDES    %NRPLWV_ROW_DATA_POINTS    ,&
              GDES    %NRPLWV_COL_DATA_POINTS    ,&
              GDES    %NPRO_ROW                  ,&
              GDES    %NPRO_COL                  ,&  
              GDES    %RES_NRPLWV_ROW            ,&
              GDES    %RES_NRPLWV_COL            ,&
              GDES    %MPLWV_ROW                 ,&
              GDES    %MPLWV_COL                 ,&
              GDES    %RES_NRPLWV_ROW_DATA_POINTS,&
              GDES    %RES_NRPLWV_COL_DATA_POINTS,&
              GDES    %NLM_ROW                   ,&
              GDES    %NLM_COL                   

#endif
         
         CALL STORAGE_REQ_GG( "t", MEM_GG, GDES_TAU, GDES, IGRIDS%T, IGRIDS%B, S2E, WDES, WGW,  S )

         ! drop a warning if MAXMEM is too low 
         IF( MAXMEM < MEM_GG%M_TOTAL ) THEN
            ALLOCATE(ARG%IVAL(2))
            ARG%IVAL(1)=MAXMEM
            ARG%IVAL(2)=INT(MEM_GG%M_TOTAL)
            CALL VTUTOR%WRITE(ISALERT, MAXMEM_LOW, ARG)
            DEALLOCATE( ARG%IVAL )
         ENDIF 

                ! rank         group       
         IF ( IO%IU0>=0 ) WRITE(IO%IU0, 990 ) MEM_GG%M_TOTAL, MEM_GG%M_TOTAL*NCPU 
         IF ( IO%IU6>=0 ) WRITE(IO%IU6, 990 ) MEM_GG%M_TOTAL, MEM_GG%M_TOTAL*NCPU 

      ENDIF 

990   FORMAT( /," estimated memory requirement per rank ", &
           F8.1," MB, per node ", F8.1, " MB",/ )

   END SUBROUTINE DUMP_STORAGE_REQUIREMENTS


   ! dumps diagonal matrix elements of self-energy

   SUBROUTINE DUMP_DIAG( SIGMA_COS, SIGMA_SIN, IMAG_GRIDS, IO )
      USE minimax_struct, ONLY: imag_grid_handle
      COMPLEX(q)               :: SIGMA_COS(:,:,:,:)
      COMPLEX(q)               :: SIGMA_SIN(:,:,:,:)
      TYPE( imag_grid_handle ) :: IMAG_GRIDS
      TYPE( in_struct )        :: IO
      ! local 
      INTEGER                  :: I 
      INTEGER                  :: NB, NK,ISP
      REAL(q), ALLOCATABLE     :: SIGMAO(:,:,:) 
      REAL(q), ALLOCATABLE     :: SIGMAU(:,:,:) 
      
      IF ( IO%IU0>=0 ) THEN
      OPEN( 100, FILE='sigmaw',STATUS='REPLACE')
      DO ISP = 1, SIZE( SIGMA_COS,4)
         !DO NB = 1, SIZE( SIGMA_COS,2)
         DO NB = 2, 2
            DO NK = 1, SIZE( SIGMA_COS,3)
            WRITE(100,*)ISP,NB,NK
            DO I = 1, IMAG_GRIDS%NOMEGA 
               WRITE(100,'(5F25.16)') IMAG_GRIDS%FER_RE( I ), SIGMA_COS(I,NB,NK, ISP ), SIGMA_SIN( I, NB,NK, ISP )
            ENDDO
            ENDDO
         ENDDO
      ENDDO
      CLOSE(100)
      ENDIF
   END SUBROUTINE DUMP_DIAG

!************************ SUBROUTINE EDDIAG_SIMPLE *********************
!
! this subroutine calculates the matrix
!
!  <i | T + V_local + V_exact_exchange | a> = 
!  <i | T + V_ext + V_hartree + V_exact_exchange | a>
!
! and stores the matrix in a scaLAPACK compatible fashion
! it is essentially a copy of EDDIAG but cleaned up to do only
! what it needs to do for self consistent GW calculations 
! 
! if LGW0 is set, the routines also diagonalizes the Hamiltonian
! and sets WMEAN%CELTOT accordingly
! 
!***********************************************************************

  SUBROUTINE EDDIAG_SIMPLE(HAMILTONIAN, &
       GRID, LATT_CUR, NONLR_S, NONL_S, W, WMEAN, WHF_IN_KS, WDES, LGW0_, SYMM, &
       LMDIM, CDIJ, CQIJ, SV, T_INFO, P, NKPTS_IRZ, CHAM_MAT, GDES_MAT, E, IU6 )
    USE prec
    USE wave_high
    USE lattice
    USE nonl_high
    USE hamil
    USE hamil_struct_def
    USE scala
    USE main_mpi
    USE fock
    USE pseudo
    USE ini
    USE sym_grad
    IMPLICIT NONE
    TYPE (ham_handle)  HAMILTONIAN
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W
    TYPE (wavespin)    WMEAN         ! identical to W, except for CELTOT and FERTOT which are used to store mean field
                                     ! part of Hamiltonian
    TYPE (wavespin)    WHF_IN_KS     ! identical to W, except for CELTOT and FERTOT (CELTOT is set to diagonals of HF Hamiltonian)
    TYPE (wavedes)     WDES
    LOGICAL            LGW0_         ! diagonalize mean field Hamiltonian
    TYPE (symmetry) :: SYMM
    INTEGER LMDIM
    OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ),CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    RGRID   SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ)   ! local potential
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    GDEF, POINTER   :: CHAM_MAT(:,:,:,:)         ! Hamilton matrix (time independent part) T+ V_x + V_local
    TYPE (greens_mat_des), POINTER :: GDES_MAT
    INTEGER         :: NKPTS_IRZ
    TYPE (energy)   :: E
    INTEGER         :: IU6    
    
! local variables
    TYPE (wavedes1)    WDES1          ! descriptor for one k-point
    TYPE (wavefun1)    W1             ! structure to hold one orbital
    TYPE (wavefuna)    WA             ! subpointer to part of W
    TYPE (wavespin)    WFOCK          ! array to store the Fock (exchange contribution)
    TYPE (wavefuna)    WNONL          ! array to hold non local part D * wave function character
    TYPE (wavefuna)    WHAM           ! array to store accelerations for a selected block
    COMPLEX(q)         CDCHF          ! HF double counting energy
    REAL(q)         :: RR(WDES%NB_TOT)
    GDEF            :: R(WDES%NB_TOT)
    INTEGER         :: NK, ISP
    INTEGER :: NBANDS, NBANDSK, NSTRIP, N, NP, NPOS, NSTRIP_ACT, &
         NPOS_RED, NSTRIP_RED, NSIM_LOCAL

    PROFILING_START('eddiag_simple')

    CDCHF=0

    NBANDS=WDES%NBANDS
    NSTRIP=NSTRIP_STANDARD

    CALL SETWDES(WDES,WDES1,0)
 
    ! allocate work space
    ALLOCATE(W1%CR(GRID%MPLWV*WDES%NRSPINORS))
    CALL NEWWAVA(WHAM, WDES1, NSTRIP)
    CALL NEWWAVA_PROJ(WNONL, WDES1)

    ! force full fock contribution
    HFSCREEN=0
    AEXX    =1._q
    E%EBANDSTR=0
    E%EXHF_ACFDT=0

!=======================================================================
! start with HF part and local part and store results WFOCK
! since the orbitals are later redistributed we need
! to do these calculations before hand
!=======================================================================
    CALL ALLOCW(WDES,WFOCK)
    WFOCK%CW=0
    WFOCK%CPROJ=0

    NSIM_LOCAL=(W%WDES%NSIM*2+W%WDES%NB_PAR-1)/W%WDES%NB_PAR

    DO ISP=1,WDES%ISPIN
    DO NK=1,NKPTS_IRZ
       CHAM_MAT(:,NK,ISP,1)=0
       ! reduce bands to NB_TOTK (last usefull orbital)
       NBANDSK=MIN((WDES%NB_TOTK(NK,ISP)+WDES%NB_PAR-1)/WDES%NB_PAR, NBANDS)

       IF (MOD(NK-1,WDES%COMM_KINTER%NCPU)/=WDES%COMM_KINTER%NODE_ME-1) CYCLE

       CALL SETWDES(WDES,WDES1,NK)

       IF (NONLR_S%LREAL) THEN
          CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES)
       ELSE
          CALL PHASE(WDES,NONL_S,NK)
       ENDIF

       DO NPOS=1,NBANDSK,NSIM_LOCAL
          NSTRIP_ACT=MIN(NBANDSK+1-NPOS,NSIM_LOCAL)
          CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W,  &
               NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIP_ACT, &
               WFOCK%CW(:,NPOS:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF, EXHF_ACFDT=E%EXHF_ACFDT )
          ! add local part
          DO N=NPOS,NPOS+NSTRIP_ACT-1
             NP=N-NPOS+1
             CALL SETWAV(W, W1, WDES1, N, ISP)
             CALL FFTWAV_W1(W1)
! this flag allows to calculate the exchange interaction only
! can be also set in GG_base.F to get half the exchange energy from the correlation part
#ifndef onlyHF
             CALL HAMILT_LOCAL(W1, SV, ISP, WFOCK%CW(:,N,NK,ISP), .TRUE.)
#endif
          ENDDO
       ENDDO

!=======================================================================

       CALL SETWDES(WDES,WDES1,NK)

!=======================================================================
!  IFLAG /= 0 calculate Hamiltonian CHAM_MAT
!=======================================================================
       WA=ELEMENTS(W, WDES1, ISP)

       CALL OVERL(WDES1, .TRUE.,LMDIM,CDIJ(1,1,1,ISP), WA%CPROJ(1,1),WNONL%CPROJ(1,1))
#ifdef onlyHF
       WNONL%CPROJ=0
#endif

       ! redistribute the orbitals and orbital characters 
       CALL REDISTRIBUTE_PW( WA)
       CALL REDISTRIBUTE_PROJ(WA)
       CALL REDISTRIBUTE_PROJ(WNONL)

       stripe: DO NPOS=1,NBANDSK,NSTRIP
          NSTRIP_ACT=MIN(NBANDSK+1-NPOS,NSTRIP)
          
          !  calculate V_{local} |phi> + T | phi >
          !  for a block containing NSTRIP orbitals
          
          ! set Fock contribution
          WHAM%CW(:,1:NSTRIP_ACT)=WFOCK%CW(:,NPOS:NPOS+NSTRIP_ACT-1,NK,ISP)
          
          ! redistribute H |orbital>
          CALL REDISTRIBUTE_PW( ELEMENTS( WHAM, 1, NSTRIP_ACT))
          
          NPOS_RED  =(NPOS-1)*WDES%NB_PAR+1
          NSTRIP_RED=NSTRIP_ACT*WDES%NB_PAR

          CALL ORTH1_DISTRI_HERM_DESC('U', &
               WA%CW_RED(1,1),WHAM%CW(1,1),WA%CPROJ_RED(1,1), &
               WNONL%CPROJ_RED(1,NPOS_RED),WDES%NB_TOT, &
               NPOS_RED, NSTRIP_RED, WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CHAM_MAT(1,NK,ISP,1), & 
               WDES%COMM_KIN, WDES%NB_TOTK(NK,ISP), GDES_MAT%DESC)
       ENDDO stripe

       !  back redistribution over bands
       CALL REDISTRIBUTE_PW( WA)
       CALL REDISTRIBUTE_PROJ( WA)

       ! for testing compatibility with Loewdin perturbation theory in subrot.F
       R=zero
       CALL  DETERMINE_DIAGONALE_GDEF( WDES1%NB_TOTK(ISP), CHAM_MAT(:,NK,ISP,1),R(:), GDES_MAT%DESC)
       CALLMPI( M_sum_g(GDES_MAT%COMM_INTRA, R(1), WDES1%NB_TOTK(ISP)))
       WMEAN%CELTOT(:, NK, ISP) = R(:)
       WHF_IN_KS%CELTOT=WMEAN%CELTOT

       ! evaulate E_HF^GM = Tr { GAMMA ( h + F ) }, h+F is the HF Hamiltonian
       ! note that GAMMA is the density matrix of electrons (tau->0-)
       ! store the result in E%BANDSTR
       DO N=1,WDES%NB_TOTK(NK,ISP)
          ! calculate sum of diagonal components of HF Hamiltonian times occupancies
          ! this is fully compatible with ALGO = Eigenval
          E%EBANDSTR=E%EBANDSTR+WDES%RSPIN* REAL( WMEAN%CELTOT(N, NK, ISP) ,KIND=q) *WDES%WTKPT(NK)*W%FERTOT(N,NK,ISP)
       ENDDO
       
       !CALL DUMP_GREENS_HAM( "mean HF orig", CHAM_MAT(:,NK,ISP,1), WDES%NB_TOTK(NK,ISP), GDES_MAT, 6)
!       CALL DUMP_GREENS_HAM( "mean HF orig", CHAM_MAT(:,NK,ISP,1), WDES%NB_TOTK(NK,ISP), GDES_MAT, IU6)
       IF (LGW0_) THEN
          ! diagonalize Hartree-Fock Hamiltonian and store unitary matrix in CHAM_MAT
          CALL PDSYEVX_ZHEEVX_DESC(CHAM_MAT(:,NK,ISP,1), RR,  W%WDES%NB_TOTK(NK, ISP), GDES_MAT%DESC, GDES_MAT%COMM_INTRA)
          ! bcast from root group to all other groups to avoid issues with degenerate orbitals
          CALLMPI( M_bcast_g_from(GDES_MAT%COMM_INTER, CHAM_MAT(1,NK,ISP,1), SIZE(CHAM_MAT,1), 1))
          CALLMPI( M_bcast_d_from(GDES_MAT%COMM_INTER, RR, SIZE(RR), 1))

          ! store eigenvalues
          WMEAN%CELTOT(:, NK, ISP)=RR
       ENDIF
       WMEAN%CELTOT( W%WDES%NB_TOTK(NK, ISP)+1:, NK, ISP) = 10000.0
!       CALL DUMP_GREENS_HAM( "mean HF, maybe diagonalized", CHAM_MAT(:,NK,ISP,1), WDES%NB_TOTK(NK,ISP), GDES_MAT, IU6)
    END DO
    END DO
    CALLMPI( M_sum_d(WDES%COMM_KINTER, E%EBANDSTR, 1))

    CALLMPI( M_sum_z(WDES%COMM_KIN,CDCHF,1))
    CALLMPI( M_sum_z(WDES%COMM_KINTER, CDCHF,1))

    CALLMPI( M_sum_d(WDES%COMM_KIN,E%EXHF_ACFDT,1))
    CALLMPI( M_sum_d(WDES%COMM_KINTER,E%EXHF_ACFDT,1))

    ! global sum, if KPAR is set
    IF (WDES%COMM_KINTER%NCPU>1) THEN
       CALLMPI( M_sum_g(WDES%COMM_KINTER, CHAM_MAT, SIZE(CHAM_MAT)))
       CALL KPAR_SYNC_CELTOT(WDES,WMEAN)  ! just in case CELTOT has been updated
    ENDIF
       
    E%EXHF=CDCHF

    ! deallocation ...
    DEALLOCATE(W1%CR)
    CALL DELWAVA(WHAM)
    CALL DELWAVA_PROJ(WNONL)
    CALL DEALLOCW(WFOCK)

    PROFILING_STOP('eddiag_simple')

  END SUBROUTINE EDDIAG_SIMPLE

!************************ SUBROUTINE EDDIAG_RSE ************************
!
! this subroutine calculates the "renormalized" singles contribution
! to the correlation energy in GW calculations, specifically
!
!  <i | T + V_local + V_exact_exchange | a> - \delta_ai epsilon_a
!
! the routine is mostly a copy of EDDIAG_SIMPLE with some modifications
! to avoid storing the Hamiltonian matrix for each k-points
!
!***********************************************************************

  SUBROUTINE EDDIAG_RSE(HAMILTONIAN, &
       GRID, LATT_CUR, NONLR_S, NONL_S, W, WDES, SYMM, &
       LMDIM, CDIJ, CQIJ, SV, T_INFO, P, NKPTS_IRZ, E, CSINGLES, CSINGLE_SHOT)
    USE prec
    USE wave_high
    USE lattice
    USE nonl_high
    USE hamil_struct_def
    USE hamil
    USE scala
    USE main_mpi
    USE fock
    USE pseudo
    USE ini
    USE sym_grad
    IMPLICIT NONE
    TYPE (ham_handle)  HAMILTONIAN
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W
    TYPE (wavedes)     WDES
    TYPE (symmetry) :: SYMM
    INTEGER LMDIM
    OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ),CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    RGRID   SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ)   ! local potential
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    GDEF, DIMENSION(:), POINTER   :: CHAM =>NULL()                  ! Hamilton matrix (time independent part) T+ V_x + V_local
    GDEF, DIMENSION(:), POINTER   :: CUNI =>NULL()                  ! unitary matrix diagonalizing the Hamiltonian
    GDEF, DIMENSION(:), POINTER   :: CTMP =>NULL()                  ! temporary matrix
    INTEGER         :: NKPTS_IRZ
    TYPE (energy)   :: E
    
! local variables
    TYPE (wavedes1)    WDES1          ! descriptor for one k-point
    TYPE (wavefun1)    W1             ! structure to hold one orbital
    TYPE (wavefuna)    WA             ! subpointer to part of W
    TYPE (wavespin)    WFOCK          ! array to store the Fock (exchange contribution)
    TYPE (wavefuna)    WNONL          ! array to hold non local part D * wave function character
    TYPE (wavefuna)    WHAM           ! array to store accelerations for a selected block
    COMPLEX(q)         CDCHF          ! HF double counting energy
    COMPLEX(q)     :: CSINGLES        ! singles contribution
    COMPLEX(q)     :: CSINGLE_SHOT    ! sum of HF eigenvalue- original values in diagonal of Hamiltonian
    COMPLEX(q)     :: CSUM            ! temporary
    REAL(q)  :: R(WDES%NB_TOT)        ! eigenvalues of Hamiltonian
    GDEF            :: GDIAG(WDES%NB_TOT)  ! diagonal of Hamiltonian
    INTEGER         :: NK, ISP
    INTEGER :: NBANDS, NBANDSK, NSTRIP, N, NP, NPOS, NSTRIP_ACT, &
         NPOS_RED, NSTRIP_RED, NSIM_LOCAL

    PROFILING_START('eddiag_rse')

    CDCHF=0

    NBANDS=WDES%NBANDS
    NSTRIP=NSTRIP_STANDARD

    CALL SETWDES(WDES,WDES1,0)
 
    ! use "standard" VASP scalapack descriptors
    CALL INIT_scala(WDES%COMM_KIN, WDES%NB_TOT)
    ALLOCATE(CHAM(SCALA_NP()*SCALA_NQ()), CUNI(SCALA_NP()*SCALA_NQ()))
    ALLOCATE(CTMP(SCALA_NP()*SCALA_NQ()))
    ! allocate work space
    ALLOCATE(W1%CR(GRID%MPLWV*WDES%NRSPINORS))

    CALL NEWWAVA(WHAM, WDES1, NSTRIP)
    CALL NEWWAVA_PROJ(WNONL, WDES1)

    ! force full fock contribution
    HFSCREEN=0.0
    AEXX    =1.0
    E%EBANDSTR=0
    CSINGLES=0
    CSINGLE_SHOT=0
!=======================================================================
! start with HF part and local part and store results WFOCK
! since the orbitals are later redistributed we need
! to these calculations before hand
!=======================================================================
    CALL ALLOCW(WDES,WFOCK)
    NSIM_LOCAL=(W%WDES%NSIM*2+W%WDES%NB_PAR-1)/W%WDES%NB_PAR

    DO ISP=1,WDES%ISPIN
    DO NK=1,NKPTS_IRZ
       ! set present number of bands
       CALL INIT_scala(WDES%COMM_KIN, WDES%NB_TOTK(NK,ISP))
       CHAM=0
       ! reduce bands to NB_TOTK (last usefull orbital)
       NBANDSK=MIN((WDES%NB_TOTK(NK,ISP)+WDES%NB_PAR-1)/WDES%NB_PAR, NBANDS)

       IF (MOD(NK-1,WDES%COMM_KINTER%NCPU)/=WDES%COMM_KINTER%NODE_ME-1) CYCLE

       CALL SETWDES(WDES,WDES1,NK)

       IF (NONLR_S%LREAL) THEN
          CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES)
       ELSE
          CALL PHASE(WDES,NONL_S,NK)
       ENDIF
       
       DO NPOS=1,NBANDSK,NSIM_LOCAL
          NSTRIP_ACT=MIN(NBANDSK+1-NPOS,NSIM_LOCAL)
          CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W,  &
               NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIP_ACT, &
               WFOCK%CW(:,NPOS:,NK,ISP), P, CQIJ(1,1,1,1), CDCHF, EXHF_ACFDT=E%EXHF_ACFDT )
          ! add local part
          DO N=NPOS,NPOS+NSTRIP_ACT-1
             NP=N-NPOS+1
             CALL SETWAV(W, W1, WDES1, N, ISP)
             CALL FFTWAV_W1(W1)
             CALL HAMILT_LOCAL(W1, SV, ISP, WFOCK%CW(:,N,NK,ISP), .TRUE.)
          ENDDO
       ENDDO
!=======================================================================

       CALL SETWDES(WDES,WDES1,NK)

!=======================================================================
!  IFLAG /= 0 calculate Hamiltonian CHAM
!=======================================================================
       WA=ELEMENTS(W, WDES1, ISP)

       CALL OVERL(WDES1, .TRUE.,LMDIM,CDIJ(1,1,1,ISP), WA%CPROJ(1,1),WNONL%CPROJ(1,1))

       ! redistribute the orbitals and orbital characters 
       CALL REDISTRIBUTE_PW( WA)
       CALL REDISTRIBUTE_PROJ(WA)
       CALL REDISTRIBUTE_PROJ(WNONL)

       strip: DO NPOS=1,NBANDSK,NSTRIP
          NSTRIP_ACT=MIN(NBANDSK+1-NPOS,NSTRIP)
          
          !  calculate V_{local} |phi> + T | phi >
          !  for a block containing NSTRIP orbitals
          
          ! set Fock contribution
          WHAM%CW(:,1:NSTRIP_ACT)=WFOCK%CW(:,NPOS:NPOS+NSTRIP_ACT-1,NK,ISP)
          
          ! redistribute H |orbital>
          CALL REDISTRIBUTE_PW( ELEMENTS( WHAM, 1, NSTRIP_ACT))
          
          NPOS_RED  =(NPOS-1)*WDES%NB_PAR+1
          NSTRIP_RED=NSTRIP_ACT*WDES%NB_PAR
          
          CALL ORTH1_DISTRI_HERM_DESC('L', &
               WA%CW_RED(1,1),WHAM%CW(1,1),WA%CPROJ_RED(1,1), &
               WNONL%CPROJ_RED(1,NPOS_RED),WDES%NB_TOT, &
               NPOS_RED, NSTRIP_RED, WDES1%NPL_RED,WDES1%NPRO_RED,WDES1%NRPLWV_RED,WDES1%NPROD_RED,CHAM(1), & 
               WDES%COMM_KIN, WDES%NB_TOTK(NK,ISP), DESCSTD)
       ENDDO strip
       
       !  back redistribution over bands
       CALL REDISTRIBUTE_PW( WA)
       CALL REDISTRIBUTE_PROJ( WA)

       NBANDSK=  W%WDES%NB_TOTK(NK, ISP)

       NPOS=0 ! set NPOS to last occupied orbital
       DO N=1,NBANDSK
          IF ( W%FERTOT(N, NK, ISP)>0.5) THEN
             NPOS=N
          ENDIF
       ENDDO

       CALL DETERMINE_DIAGONALE_GDEF( NBANDSK, CHAM(:), GDIAG(:), DESCSTD)
       CALLMPI( M_sum_g(WDES%COMM_KIN, GDIAG(1), SIZE(GDIAG,1)))

       ! determine band structure energy
       DO N=1,NBANDSK
          E%EBANDSTR=E%EBANDSTR+WDES%RSPIN* GDIAG(N)*WDES%WTKPT(NK)*W%FERTOT(N,NK,ISP)
       ENDDO

       ! get HF eigenvalues by diagonalization (not sure about performance
       ! maybe faster to call PDSYEVX_ZHEEVX_DESC or PDSYEV_ZHEEVD_DESC
       ! however those crash for 'N')
       CUNI=CHAM
       CALL PDSYEV_ZHEEV_DESC(CUNI(:), R, NBANDSK, DESCSTD, WDES%COMM_KIN, 'N')

       ! determine single shot band structure energy minus diagonals summed above
       DO N=1,NBANDSK
          CSINGLE_SHOT=CSINGLE_SHOT+WDES%RSPIN*(R(N)-GDIAG(N))*WDES%WTKPT(NK)*W%FERTOT(N,NK,ISP)
       ENDDO

       ! "renormalized" perturbation theory ala Ren et al.
       ! remove the occ-unocc blocks (no obvious generalization to partial occupancies
       CUNI=CHAM
       CALL CLEAR_OCCUNOCC_UNOCCOCC(NBANDSK, CUNI(:), NPOS, DESCSTD)
       ! determine rotation matrix in occupied-occupied and unoccupied-unoccupied block
       CALL PDSYEVX_ZHEEVX_DESC(CUNI(:), R, NBANDSK, DESCSTD, WDES%COMM_KIN)

       ! subtract original (usually DFT eigenvalues) from HF-Hamiltonian
       CALL ADD_TO_DIAGONALE_REAL( NBANDSK, CHAM(:), -REAL(W%CELTOT(:, NK, ISP)), DESCSTD)
       ! now rotate the <a |H |i> - epsilon_DFT delta_ai using previously determined unitary matrix
       ! TMP= U+ H
       ! mK: use following order, for some reason gnu+mkl openmpi-4 breaks after the
       ! first PGGEM call
       ! TMP+= H+ U
       CTMP = 0
       CALL PGGEMM(  trans, 'N',  NBANDSK, NBANDSK, NBANDSK, one, &
            CHAM(1), 1, 1, DESCSTD, &
            CUNI(1), 1, 1, DESCSTD, &
            zero,  &
            CTMP(1), 1, 1, DESCSTD )
       ! H = TMP+ U                                         
       CALL PGGEMM( trans, 'N', NBANDSK, NBANDSK, NBANDSK, one, &
            CTMP(1), 1, 1, DESCSTD, &
            CUNI(1), 1, 1, DESCSTD, &
            zero,  &
            CHAM(1), 1, 1, DESCSTD )

!       CALL DUMP_HAM_DISTRI( 'U+H U', WDES, CHAM, NBANDSK, DESCSTD, 100+WDES%COMM%NODE_ME)
       ! perturbation theory using the "updated" renormalized HF eigenvalues
       CALL SUM_SECOND_ORDER(NBANDSK, CHAM, R, W%FERTOT(:, NK, ISP), CSUM,&
          DESCSTD, LFINITE_TEMPERATURE )
       ! factor 2 because lower and upper square are included
       CSINGLES=CSINGLES+CSUM*WDES%RSPIN*WDES%WTKPT(NK)/2
    END DO
    END DO

    CALLMPI( M_sum_d(WDES%COMM_KINTER, E%EBANDSTR, 1))

    CALLMPI( M_sum_z(WDES%COMM_KIN,CDCHF,1))
    CALLMPI( M_sum_z(WDES%COMM_KINTER, CDCHF,1))

    CALLMPI( M_sum_d(WDES%COMM_KIN,E%EXHF_ACFDT,1))
    CALLMPI( M_sum_d(WDES%COMM_KINTER,E%EXHF_ACFDT,1))

    CALLMPI( M_sum_z(WDES%COMM_KIN, CSINGLES, 1))    ! sum within k-point (within scaLAPACK universe)
    CALLMPI( M_sum_z(WDES%COMM_KINTER, CSINGLES, 1)) ! sum between k-points
    CALLMPI( M_sum_z(WDES%COMM_KINTER, CSINGLE_SHOT, 1)) ! sum between k-points
       
    E%EXHF=CDCHF

    ! deallocation ...
    DEALLOCATE(CHAM)
    DEALLOCATE(CUNI)
    DEALLOCATE(CTMP)
    DEALLOCATE(W1%CR)
    CALL DELWAVA(WHAM)
    CALL DELWAVA_PROJ(WNONL)
    CALL DEALLOCW(WFOCK)

    PROFILING_STOP('eddiag_rse')

  END SUBROUTINE EDDIAG_RSE

!***********************************************************************
!
! calculate the response function 
! from the Green's function using the supercell approach
! this version is fairly simple, since it follows in virtually
! all areas the 
!
!***********************************************************************

  SUBROUTINE CALCULATE_RESPONSE_SUPER( S, W, WGW, GO, GU, GDES, CHI, & 
       NQMAX, NQ_INDEX, LATT_CUR, IU6)
    USE ini
    USE dfast
    USE fock
    USE string, ONLY: str
    TYPE (supercell), POINTER :: S
    TYPE (wavespin):: W            ! one electron orbitals
    TYPE (wavedes) :: WGW          ! descriptor for response function in IRZ
    TYPE (greensf) :: GU(:)        ! Green function for unoccupied orbitals 
    TYPE (greensf) :: GO(:)        ! Green function for occupied orbitals
    TYPE (greensfdes) :: GDES      ! descriptor for Green function
    TYPE (responsefunction) CHI(:) ! response function at each q point in the IRZ
    INTEGER :: NQMAX               ! total number of q-points that are considered
    INTEGER :: NQ_INDEX(:)         ! index of the corresponding q-point in the IRZ (used when NKRED is set)
    TYPE (latt)        LATT_CUR
    INTEGER :: IU6
  ! local
    TYPE (wavedes1) :: WDESQ
    INTEGER :: NK, NR, NRP, NT, NI, NIS, LMMAXC, L, LP, LM, LMBASE, NAUG, NT_GLOBAL, NQ_COUNTER, NQ

    GDEF, ALLOCATABLE :: GWORK(:,:)
#ifdef RPAgamma
    REAL(q), ALLOCATABLE :: GWORK_REAL(:,:)
#else
    GDEF, ALLOCATABLE :: GWORK_REAL(:,:)
#endif
    GDEF, ALLOCATABLE :: CRHOLM(: __omp_arg(:)) ! augmentation occupancy matrix
    GDEF, ALLOCATABLE :: P_NLM(:,:)
    GDEF, ALLOCATABLE :: P_G_PROJ(:,:,:)
    GDEF, ALLOCATABLE :: P_G_R(:,:,:)
    GDEF, ALLOCATABLE :: P_PROJ_G(:,:)
    GDEF, ALLOCATABLE :: P_R_G(:,:)
    
    GDEF, ALLOCATABLE :: GO_GK_PROJ(:,:)
    GDEF, ALLOCATABLE :: GO_R_PROJ(:,:)
    GDEF, ALLOCATABLE :: GO_PROJ_PROJ(:,:)
    GDEF, ALLOCATABLE :: GU_GK_PROJ(:,:)
    GDEF, ALLOCATABLE :: GO_GU_PROJ(:)
    GDEF, ALLOCATABLE :: GU_R_PROJ(:,:)
    GDEF, ALLOCATABLE :: GU_PROJ_PROJ(:,:)

    GDEF, ALLOCATABLE :: P_TMP(:)

    TYPE (wavefun1), POINTER, CONTIGUOUS :: W1(:), W2(:)
    LOGICAL LPHASE
    COMPLEX(q) :: CPHASE(WGW%GRID%RL%NP)
    INTEGER :: ISTAT
    INTEGER(qi8) :: NTMAX 
    REAL(q) :: PSS 
!$  INTEGER ::  __omp_id
!$  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS

    PROFILING_START('calculate_response_super')
!#define mem_profiling
#ifdef mem_profiling
    CALL GET_PSS( PSS, IU6, IU6, "SUPER START", 0 )
#endif

    ALLOCATE(CRHOLM(S%AUG_DES%NPRO __omp_get_max_threads))  ! CRHOLM total number of channels in supercell

    CALL NEWWAV_OMP(W1, S%WDES1,.TRUE.)
    CALL NEWWAV_OMP(W2, S%WDES1,.TRUE.)
    
       ! TODO consider spinors later
       ! DO ISPINOR =0,0 ! GDES%NRSPINORS-1
!=======================================================================
! contract over first direction
!=======================================================================
!=======================================================================
! first part where second index involves terms such as  
!    G(X,N alpha) G+(X,N alpha') T^NLM N,alpha alpha'
!=======================================================================
       IF (S%WDES1%LOVERL) THEN
       ALLOCATE(P_G_PROJ (GDES%RES_NRPLWV_ROW_DATA_POINTS, GDES%NLM_COL, NQMAX), STAT=ISTAT)
       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "CALCULATE_RESPONSE_SUPER (P_G_PROJ) is not able to allocate "//&
          str( ndata * 8._q*GDES%RES_NRPLWV_ROW_DATA_POINTS*GDES%NLM_COL*NQMAX/1024)  // &
          " kB of data on MPI rank 0.")
       ENDIF
       P_G_PROJ=0
    ! contract over beta, beta' to N'L'M'
       LMBASE =0
       NAUG=0
       NIS =1

       NTMAX = 0 

#ifdef mem_profiling
       CALL GET_PSS( PSS, IU6, IU6, "P_G_PROJ:",0 )
#endif
       DO NT=1,GDES%NTYP  ! note NT is a local type index in the primitive cell
                          ! and might not span all atomic types
          NT_GLOBAL=GDES%NT_GLOBAL(NT)  ! global type index

          ALLOCATE(P_NLM(S%AUG_DES%NPRO, GDES%NLM_LMMAX(NT)),  & 
                   P_TMP(m_ S%WGW1%GRID%RC%NP))  ! double allocation since declared as GDEF
          ! S%WGW%GRID%MPLWV is the number of complex words required to do in place FFT
#ifdef RPAgamma
          ALLOCATE(GWORK_REAL(2* S%WGW%GRID%MPLWV,GDES%NLM_LMMAX(NT)), STAT=ISTAT)
          IF ( ISTAT/=0 ) THEN
             CALL VTUTOR%ERROR( "CALCULATE_RESPONSE_SUPER (GWORK_REAL) is not able to allocate "//&
             str( 16._q* S%WGW%GRID%MPLWV*GDES%NLM_LMMAX(NT)/1024)  // &
             " kB of data on MPI rank 0.")
          ENDIF
          NTMAX = MAX( NTMAX , 16*SIZE( GWORK_REAL ,KIND=qi8) )
#else
          ALLOCATE(GWORK_REAL(m_ S%WGW%GRID%MPLWV,GDES%NLM_LMMAX(NT)), STAT=ISTAT)
          IF ( ISTAT/=0 ) THEN
             CALL VTUTOR%ERROR( "CALCULATE_RESPONSE_SUPER (GWORK_REAL) is not able to allocate "//&
             str( m_ 8._q* S%WGW%GRID%MPLWV*GDES%NLM_LMMAX(NT)/1024)  // &
             " kB of data on MPI rank 0.")
          ENDIF
          NTMAX = MAX( NTMAX , m_ 8* SIZE( GWORK_REAL ,KIND=qi8) )
#endif
          CALL REGISTER_ALLOCATE(REAL(NTMAX,q), "GGwork")

          ! TODO: remove this consistency check at a later point
          IF (AUG_DES%LMMAX(NT_GLOBAL)/= GDES%NLM_LMMAX(NT)) THEN
             CALL vtutor%bug("internal error 1 in CALCULATE_RESPONSE_FROM_G: " // &
                str(AUG_DES%LMMAX(NT_GLOBAL)) // " " // str(GDES%NLM_LMMAX(NT)), __FILE__, __LINE__)
          ENDIF

          LMMAXC=GDES%LMMAX(NT)
          IF (LMMAXC/=0) THEN
             ALLOCATE(GO_GK_PROJ(S%WDES1%NGVECTOR ,LMMAXC), & ! stores orbital(g+k, alpha) see below
                  GO_R_PROJ(S%WDES1%GRID%RL%NP ,LMMAXC), & ! stores orbital(r, alpha) see below
                  GO_PROJ_PROJ(S%WDES1%NPRO   ,LMMAXC), &  ! stores projectors(N'' alpha')
                  GU_GK_PROJ(S%WDES1%NGVECTOR ,LMMAXC), &  ! stores orbital(g+k, alpha) see below
                  GU_R_PROJ(S%WDES1%GRID%RL%NP ,LMMAXC), & ! stores orbital(r, alpha) see below
                  GO_GU_PROJ(S%WDES1%GRID%RL%NP), &        ! temporary
                  GU_PROJ_PROJ(S%WDES1%NPRO   ,LMMAXC))    ! stores projectors(N'' alpha')
             !
             ! main loop over local ion index N (second = column index)
             ! this index is in the primitive cell
             !
             DO NI=NIS,GDES%NITYP(NT)+NIS-1
                ! merge orbital(g+k, alpha)         =  G_k(g, N alpha) from G_PROJ
                ! merge projectors(N'' alpha') = sum_k exp^ik(R_N''-R_N') G_k(N' alpha', N alpha) from PROJ_PROJ
                ! phase factors are added for the projector parts 
                CALL COLLECT_PROJECTORS_SUPER(S, W%WDES, GO_GK_PROJ, GO_PROJ_PROJ, GO(:), & 
                     LMBASE, LMMAXC, GDES%POSION(:,NI))
                CALL COLLECT_PROJECTORS_SUPER(S, W%WDES, GU_GK_PROJ, GU_PROJ_PROJ, GU(:), & 
                     LMBASE, LMMAXC, GDES%POSION(:,NI))

                ! supercell FFT orbital(g+k, alpha') -> orbital(r, alpha')
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(__omp_id,L)
!$              __omp_id=OMP_GET_THREAD_NUM()+1
!$OMP DO
                DO L=1, LMMAXC
                   W1(__omp_id)%CW(1:S%WDES1%NGVECTOR)=GO_GK_PROJ(1:S%WDES1%NGVECTOR, L)
                   CALL FFTWAV(S%WDES1%NGVECTOR, S%WDES1%NINDPW(1), &
                        W1(__omp_id)%CR(1),W1(__omp_id)%CW(1),S%WDES1%GRID) ! GO_GK_PROJ(1, L),S%WDES1%GRID)
                   GO_R_PROJ(1:S%WDES1%GRID%RL%NP, L)=W1(__omp_id)%CR(1:S%WDES1%GRID%RL%NP)

                   W2(__omp_id)%CW(1:S%WDES1%NGVECTOR)=GU_GK_PROJ(1:S%WDES1%NGVECTOR, L)
                   CALL FFTWAV(S%WDES1%NGVECTOR, S%WDES1%NINDPW(1), &
                        W2(__omp_id)%CR(1),W2(__omp_id)%CW(1),S%WDES1%GRID)  ! GU_GK_PROJ(1, L),S%WDES1%GRID)
                   GU_R_PROJ(1:S%WDES1%GRID%RL%NP, L)=W2(__omp_id)%CR(1:S%WDES1%GRID%RL%NP)
                   
                   ! test consistency of GO_PROJ_PROJ with GO_GK_PROJ
                   IF (ASSOCIATED(S%NONL_S)) THEN
                      CALL PROJ1(S%NONL_S,S%WDES1,W1(__omp_id))
                      GO_PROJ_PROJ(1:S%WDES1%NPRO,L)=W1(__omp_id)%CPROJ(1:S%WDES1%NPRO)
                      CALL PROJ1(S%NONL_S,S%WDES1,W2(__omp_id))
                      GU_PROJ_PROJ(1:S%WDES1%NPRO,L)=W2(__omp_id)%CPROJ(1:S%WDES1%NPRO)
                   ENDIF
                ENDDO
!$OMP END DO
!$OMP END PARALLEL

                P_NLM=0
                GWORK_REAL=0

                ! T_beta, beta'^NLM' G(N' beta, N alpha) G(N' beta', N alpha') T_alpha, alpha'^NLM
                DO L=1,LMMAXC   ! alpha
                DO LP=1,LMMAXC  ! alpha'
                   ! contraction T_beta, beta'^NLM' G(N' beta, N alpha) G(N' beta', N alpha')
                   CALL DEPSUM_TWO_BANDS_RHOLM_TRACE( &
                     GO_PROJ_PROJ(:,L),GU_PROJ_PROJ(:,LP), S%WDES1, S%AUG_DES, &
                     S%TRANS_MATRIX, CRHOLM(: __omp_arg(1)), 1._q, S%WDES1%LOVERL)

                   ! contraction over T_alpha, alpha'^NLM
                   DO LM=1,GDES%NLM_LMMAX(NT)
                      ! index build from NI and LM
                      P_NLM(:,LM)=P_NLM(:,LM)+CRHOLM(: __omp_arg(1))*TRANS_MATRIX_FOCK(LP,L,LM,NT_GLOBAL)
                   ENDDO
                ENDDO
                ENDDO
                ! calculate G(r', N alpha) G*+(r', N alpha') T_alpha alpha'^NLM
#ifdef RPAgamma
                CALL DEPSUM_VECTOR_REAL( GO_R_PROJ, GU_R_PROJ, GWORK_REAL, S%GRID%RL%NP, TRANS_MATRIX_FOCK(:,:,:,NT_GLOBAL), LMMAXC, 0, GDES%NLM_LMMAX(NT) )
#else
                CALL DEPSUM_VECTOR( GO_R_PROJ, GU_R_PROJ, GWORK_REAL, S%GRID%RL%NP, TRANS_MATRIX_FOCK(:,:,:,NT_GLOBAL), LMMAXC, 0, GDES%NLM_LMMAX(NT) )
#endif

                ! for each atom the matrix P_NLM must be symmetric
                ! at least at the atom under consideration
!                WRITE(*,*) GDES%NLM_LMMAX(NT), AUG_DES%LMMAX(NT)
!                WRITE(*,'(9E10.2)') P_NLM(1:9, 1:9)
!                WRITE(*,*)
!                WRITE(*,'(9E10.2)') P_NLM(10:18,1:9)
!                STOP
                 
                ! now add P(N'L'M', NLM) Q_N'L'M'(r) to density in supercell direction
                IF (S%WDES1%LOVERL ) THEN
                   DO LM=1,GDES%NLM_LMMAX(NT)
                      S%AUG_DES%RINPL=1._q  ! actual multiplicator used by RACC0
#ifdef RPAgamma
                      CALL RACC0_REAL(S%FAST_AUG, S%AUG_DES, P_NLM(1,LM), GWORK_REAL(1,LM))
#else
                      CALL RACC0_HF(S%FAST_AUG, S%AUG_DES, P_NLM(1,LM), GWORK_REAL(1,LM))
#endif
                   ENDDO
                ENDIF
                ! extract data using wave function FFT (impose cutoff at the same time)
                DO LM=1,GDES%NLM_LMMAX(NT)
                   IF (NAUG+LM >  GDES%NLM_COL) THEN
                      CALL vtutor%bug("internal error 2 in CALCULATE_RESPONSE_FROM_G: " // &
                         str(NAUG+LM) // " " // str(GDES%NLM_COL), __FILE__, __LINE__)
                   ENDIF
                   ! FFT to  reciprocal space using supercell FFT: P(g+k)= sum_g+k e -i(k+g) r P(r)
                   CALL FFTEXT(S%WGW1%NGVECTOR,S%WGW1%NINDPW(1), &
                        GWORK_REAL(1,LM), &
                        P_TMP(1),S%WGW1%GRID,.FALSE.)
                   ! extract data points after FFT into appropriate 
                   ! P_k(g, r') = P(g+k )  [P_G_PROJ]
                   CALL DISTRIBUTE_RESPONSE_PROJ_SUPER(S, WGW, P_TMP, P_G_PROJ, NAUG+LM , NQMAX, NQ_INDEX, (1.0_q/S%GRID%NPLWV) )
                ENDDO

                LMBASE = LMMAXC+LMBASE
                NAUG = NAUG+ GDES%NLM_LMMAX(NT)
             ENDDO
             DEALLOCATE(GO_GK_PROJ, &
                  GO_R_PROJ, &
                  GO_PROJ_PROJ, &
                  GU_GK_PROJ, &
                  GO_GU_PROJ, &
                  GU_R_PROJ, &
                  GU_PROJ_PROJ)
          ENDIF

          NIS = NIS+GDES%NITYP(NT)
          DEALLOCATE(P_NLM, P_TMP, GWORK_REAL)
          CALL DEREGISTER_ALLOCATE(REAL(NTMAX,q), "GGwork")
       ENDDO
       ENDIF
!=======================================================================
! second part, where second index is a position index G(X,r) G+(X,r)
!=======================================================================
#ifdef mem_profiling
       CALL GET_PSS( PSS, IU6, IU6, "BEFORE_FFT",0 )
#endif


     ! FFT and transpose orbitals at all k-points to obtain G%G_R and G%PROJ_R
       DO NK=1,W%WDES%NKPTS
          ! FFT GU and GO to real space and set entries G_R, PROJ_R
          CALL FFT_G_SUPER(W%WDES, GO(NK), GDES, NK)
          CALL FFT_G_SUPER(W%WDES, GU(NK), GDES, NK)
       ENDDO
#ifdef mem_profiling
       CALL GET_PSS( PSS, IU6, IU6, "AFTER_FFT",0 )
!       CALL DUMP_FINAL_TIMING( IU6) 
!       CALL vtutor%stopCode()
#endif
       ! at this point G%G_G and G%G_PROJ can be deallocated
       DO NK=1,W%WDES%NKPTS
          CALL DEALLOCATE_GG(GU(NK))
          CALL DEALLOCATE_GG(GO(NK))
          CALL DEALLOCATE_G_PROJ(GU(NK))
          CALL DEALLOCATE_G_PROJ(GO(NK))
          CALL DEALLOCATE_PROJ_PROJ(GU(NK))
          CALL DEALLOCATE_PROJ_PROJ(GO(NK))
       ENDDO
#ifdef mem_profiling
       CALL GET_PSS( PSS, IU6, IU6, "AFTER_FFT_DEALLOC",0 )
!       CALL DUMP_FINAL_TIMING( IU6) 
!       CALL vtutor%stopCode()
#endif

       ALLOCATE(P_G_R    (GDES%RES_NRPLWV_ROW_DATA_POINTS, GDES%MPLWV_COL, NQMAX),STAT=ISTAT)
       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "CALCULATE_RESPONSE_SUPER (P_G_R) is not able to allocate "//&
          str( ndata * 8._q*GDES%RES_NRPLWV_ROW_DATA_POINTS*GDES%MPLWV_COL*NQMAX/1024)  // &
          " kB of data on MPI rank 0.")
       ENDIF
       CALL REGISTER_ALLOCATE(wsgf*SIZE(P_G_R,KIND=qi8), "GGwork")
       CALL DUMP_ALLOCATE_TAG(IU6,"RESPONSE_SUPER"); IF (IU6>=0) CALL WFORCE(IU6)

       P_G_R   =0

  ! contract over real space grid points
       S%AUG_DES%RINPL=1._q ! multiplicator used by RACC0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(GWORK_REAL,P_TMP,__omp_id,NRP)
!$     __omp_id=OMP_GET_THREAD_NUM()+1
#ifdef RPAgamma
       ALLOCATE(GWORK_REAL(2* S%WGW%GRID%MPLWV,1), P_TMP(m_ S%WGW%GRID%RC%NP))  ! double allocation, declared as REAL(q)
#else
       ALLOCATE(GWORK_REAL(m_ S%WGW%GRID%MPLWV,1), P_TMP(m_ S%WGW%GRID%RC%NP))  ! double allocation, declared as GDEF
#endif
       GWORK_REAL=0
       P_TMP=0
#ifdef mem_profiling
       CALL GET_PSS( PSS, IU6, IU6, "LOOP",0 )
!       CALL DUMP_FINAL_TIMING( IU6) 
!       CALL vtutor%stopCode()
#endif
!$OMP DO
       ! loop over all real space grid points in primitive cell r' (for r' stored locally)
       DO NRP=1,GDES%MPLWV_COL
          ! merge orbital(g+k)         =  G_k(g, r') from G_R
          ! merge projectors(N' alpha) = sum_k exp^ik(R_N'-R_N) G_k(N alpha, r') from PROJ_R
          CALL COLLECT_ORBITALS_SUPER(S, W%WDES, W1(__omp_id)%CW, W1(__omp_id)%CPROJ, GO(:), NRP)
          CALL COLLECT_ORBITALS_SUPER(S, W%WDES, W2(__omp_id)%CW, W2(__omp_id)%CPROJ, GU(:), NRP)
          ! these statements allow to calculate WX%CPROJ from the orbitals in the supercell
          ! instead of from the CPROJ in the primitive cell (usually not used since S%NONLR_S is not assoc
          IF (ASSOCIATED(S%NONL_S)) THEN
             CALL PROJ1(S%NONL_S,S%WDES1,W1(__omp_id))
             CALL PROJ1(S%NONL_S,S%WDES1,W2(__omp_id))
          ENDIF
          !WRITE (*,*) 'WPROJ1 ', SUM(W1(__omp_id)%CPROJ)
          !WRITE (*,*) 'WPROJ2 ', SUM(W2(__omp_id)%CPROJ)

          CALL FFTWAV(S%WDES1%NGVECTOR, S%WDES1%NINDPW(1), &
               W1(__omp_id)%CR(1), W1(__omp_id)%CW(1), S%WDES1%GRID)
          CALL FFTWAV(S%WDES1%NGVECTOR, S%WDES1%NINDPW(1), &
               W2(__omp_id)%CR(1), W2(__omp_id)%CW(1), S%WDES1%GRID)
#ifdef RPAgamma
          CALL PW_CHARGE_TRACE_REAL(S%WDES1, GWORK_REAL(1,1), W1(__omp_id)%CR(1), W2(__omp_id)%CR(1))
#else          
          CALL PW_CHARGE_TRACE(S%WDES1, GWORK_REAL(1,1), W1(__omp_id)%CR(1), W2(__omp_id)%CR(1))
#endif
          ! add augmentation charges
          IF (S%WDES1%LOVERL ) THEN
             CALL DEPSUM_TWO_BANDS_RHOLM_TRACE( &
                  W1(__omp_id)%CPROJ(:), W2(__omp_id)%CPROJ(:), S%WDES1, S%AUG_DES, &
                  S%TRANS_MATRIX, CRHOLM(: __omp_arg_id), 1._q, S%WDES1%LOVERL)
             S%AUG_DES%RINPL=1._q ! multiplicator used by RACC0
#ifdef RPAgamma
             CALL RACC0_REAL(S%FAST_AUG, S%AUG_DES, CRHOLM(1 __omp_arg_id), GWORK_REAL(1,1))
#else
             CALL RACC0_HF(S%FAST_AUG, S%AUG_DES, CRHOLM(1 __omp_arg_id), GWORK_REAL(1,1))
#endif
          END IF
          
          ! FFT to  reciprocal space using supercell FFT: P(g+k)= sum_g+k e -i(k+g) r P(r)
          CALL FFTEXT(S%WGW1%NGVECTOR,S%WGW1%NINDPW(1), &
               GWORK_REAL(1,1), &
               P_TMP(1),S%WGW1%GRID,.FALSE.)

          ! extract data points after FFT to correct P_G_R in primitive cell
          ! P_k(g, r') = P(g+k, r')
          CALL DISTRIBUTE_RESPONSE_SUPER(S, WGW, P_TMP, P_G_R, NRP, NQMAX, NQ_INDEX, (1.0_q/S%GRID%NPLWV) )

       ENDDO
!$OMP END DO
       DEALLOCATE(GWORK_REAL, P_TMP)
!$OMP END PARALLEL

#ifdef mem_profiling
       CALL GET_PSS( PSS, IU6, IU6, "BEFOR RELEASE_FFT_G",0 )
!       CALL DUMP_FINAL_TIMING( IU6) 
!       CALL vtutor%stopCode()
#endif

     ! free the Green's functions
       DO NK=1,W%WDES%NKPTS
          ! deallocate G_R, PROJ_R, and R_PROJ
          CALL  RELEASE_FFT_G( GU(NK))
          CALL  RELEASE_FFT_G( GO(NK))
       ENDDO
!=======================================================================
! contract over second direction, which
! now involves mostly an FFT in the second direction
!=======================================================================
       ALLOCATE(GWORK(m_ WGW%GRID%MPLWV,1))


       ALLOCATE(P_R_G    (GDES%MPLWV_ROW, GDES%RES_NRPLWV_COL_DATA_POINTS),STAT=ISTAT)
       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "CALCULATE_RESPONSE_SUPER (P_R_G) is not able to allocate "//&
          str( ndata * 8._q*GDES%MPLWV_ROW*GDES%RES_NRPLWV_COL_DATA_POINTS/1024)  // &
          " kB of data on MPI rank 0.")
       ENDIF
       IF (S%WDES1%LOVERL) THEN
          ALLOCATE(P_PROJ_G (GDES%NLM_ROW, GDES%RES_NRPLWV_COL_DATA_POINTS), STAT=ISTAT)
          IF ( ISTAT/=0 ) THEN
             CALL VTUTOR%ERROR( "CALCULATE_RESPONSE_SUPER (P_PROJ_G) is not able to allocate "//&
             str( ndata * 8._q*GDES%NLM_ROW*GDES%RES_NRPLWV_COL_DATA_POINTS/1024)  // &
             " kB of data on MPI rank 0.")
          ENDIF
       ENDIF

       CALL REGISTER_ALLOCATE(wsgf*SIZE(P_R_G,KIND=qi8), "GGwork")
       CALL DUMP_ALLOCATE_TAG(IU6,"RESPONSE_SUPER"); IF (IU6>=0) CALL WFORCE(IU6)

       DO NQ_COUNTER=1,NQMAX
          NQ=NQ_INDEX(NQ_COUNTER)
          CALL SETWDES(WGW, WDESQ, NQ )

          ! start with a transpose  (and conjugation for complex orbitals)
          CALL TRANSPOSE_G_R_RESPONSE   (P_G_R(:,:,NQ_COUNTER)   , P_R_G, GDES)
          IF (S%WDES1%LOVERL) THEN
             CALL TRANSPOSE_G_PROJ_RESPONSE(P_G_PROJ(:,:,NQ_COUNTER), P_PROJ_G, GDES)
          ENDIF

          ! set exp(-i q r)
          CALL SETPHASE_NOCHK(-WDESQ%VKPT(1:3), WDESQ%GRID, CPHASE, LPHASE)

          ! note RACC0_HF below applies the conjugated phase factor
          CALL PHASER_HF(GRIDHF, LATT_CUR, FAST_AUG_FOCK, WDESQ%VKPT(:))

!$OMP PARALLEL DEFAULT(SHARED) FIRSTPRIVATE(GWORK)
!$OMP DO 
          DO NRP=1,GDES%RES_NRPLWV_COL_DATA_POINTS
             ! copy data into GWORK
             GWORK=0
             GWORK(1:GDES%MPLWV_ROW,1)=P_R_G(1:GDES%MPLWV_ROW,NRP)
             
             ! apply phase factor to obtain cell periodic part
             ! note that this will crash at Gamma but LPHASE is .FALSE. in this case
             IF (LPHASE) THEN
                CALL APPLY_PHASE_INPLACE( WDESQ%GRID, CPHASE, GWORK(1,1))
                ! undo phase factor in augmentation part as well (e^-i k R_N)
                ! to get back "usual" VASP phase factors that lacks e^i k R_N
                IF (S%WDES1%LOVERL) THEN
                   CALL APPLY_PHASE_NLM( AUG_DES, FAST_AUG_FOCK, P_PROJ_G(:, NRP), WDESQ%VKPT(1:3), NQ)
                ENDIF
             ENDIF
             
             ! add augmentation part  P(NLM, r) Q_NLM(r') to cell periodic part
             IF (S%WDES1%LOVERL) THEN
                AUG_DES%RINPL=1._q ! multiplicator used by RACC0
                CALL RACC0_HF(FAST_AUG_FOCK, AUG_DES, P_PROJ_G(1, NRP), GWORK(1,1))
             END IF
             ! store final results in RESPONSER
             IF (CHI(NQ_COUNTER)%LREALSTORE) THEN
                CALL FFTEXT(WDESQ%NGVECTOR,WDESQ%NINDPW(1), &
                     GWORK(1,1), &
                     CHI(NQ_COUNTER)%RESPONSER(1, NRP, 1),WDESQ%GRID,.FALSE.)
                CHI(NQ_COUNTER)%RESPONSER(1:GDES%RES_NRPLWV_ROW_DATA_POINTS, NRP, 1)= & 
                     CHI(NQ_COUNTER)%RESPONSER(1:GDES%RES_NRPLWV_ROW_DATA_POINTS, NRP, 1)*(1.0_q/GRIDHF%NPLWV)
             ELSE
                CALL FFTEXT(WDESQ%NGVECTOR,WDESQ%NINDPW(1), &
                     GWORK(1,1), &
                     CHI(NQ_COUNTER)%RESPONSEFUN(1, NRP, 1),WDESQ%GRID,.FALSE.)
                CHI(NQ_COUNTER)%RESPONSEFUN(1:GDES%RES_NRPLWV_ROW_DATA_POINTS, NRP, 1)= &
                     CHI(NQ_COUNTER)%RESPONSEFUN(1:GDES%RES_NRPLWV_ROW_DATA_POINTS, NRP, 1)*(1.0_q/GRIDHF%NPLWV)
             ENDIF
          ENDDO
!$OMP END DO
!$OMP END PARALLEL
       ENDDO

       CALL DEREGISTER_ALLOCATE(wsgf*SIZE(P_R_G,KIND=qi8), "GGwork")
       DEALLOCATE(P_R_G)
       IF (S%WDES1%LOVERL) THEN
          DEALLOCATE(P_PROJ_G)
       ENDIF

    DEALLOCATE(GWORK)

    CALL DEREGISTER_ALLOCATE(wsgf*SIZE(P_G_R,KIND=qi8), "GGwork")
    DEALLOCATE(P_G_R)

    IF ( S%WDES1%LOVERL ) THEN
       DEALLOCATE(P_G_PROJ)
    ENDIF

    DEALLOCATE(CRHOLM)
    CALL DELWAV_OMP(W1, .TRUE.)
    CALL DELWAV_OMP(W2, .TRUE.)

#ifdef mem_profiling
    CALL GET_PSS( PSS, IU6, IU6, "SUPER DONE", 5 )
#endif
    PROFILING_STOP('calculate_response_super')

  END SUBROUTINE CALCULATE_RESPONSE_SUPER

!***********************************************************************
!            
! calculate the self energy from screened interaction  
! and the Green's function using the supercell approach
!   
!***********************************************************************

  SUBROUTINE CALCULATE_SIGMA_SUPER( S, W, WGW, GU, GO, GDES, &
                CHI, SIGMA, NQMAX, NQ_INDEX, LMDIM, LATT_CUR, IU6, TAU_WEIGHT, FRNL, LFORCE ,IO)
    USE ini
    USE dfast
    USE fock
    USE base , ONLY: in_struct
    IMPLICIT NONE
    TYPE (supercell), POINTER :: S
    TYPE (wavespin):: W            ! one electron orbitals
    TYPE (wavedes) :: WGW          ! descriptor for response function in IRZ
    TYPE (greensf) :: GU(:)        ! Green function for positive/negative time in orbital basis
    TYPE (greensf) :: GO(:)        ! Green function for negative/positive time in orbital basis
                                   ! quantities are contracted against GO to yield total energies
    TYPE (greensfdes) :: GDES      ! descriptor for Green function
    TYPE (responsefunction) :: CHI(:) ! response function at each q point in the IRZ
    TYPE (greensf) :: SIGMA(:)     ! Green function for occupied orbitals, first index tau (not used), second kpt
    INTEGER :: NQMAX               ! total number of q-points that are considered
    INTEGER :: NQ_INDEX(:)         ! index of the corresponding q-point in the IRZ (used when NKRED is set)
    INTEGER :: LMDIM
    TYPE (latt)        LATT_CUR
    INTEGER :: IU6
    REAL(q) :: TAU_WEIGHT
    REAL(q) :: FRNL(:,:)! forces related to rigid shift of augmentation charges
    LOGICAL :: LFORCE
    TYPE(in_struct) :: IO
  ! local variables
    TYPE (wavedes1) :: WDESQ
    INTEGER :: NK, NRP, NT, NI, NIS, LMMAXC, L, LP, LM, LMBASE, NAUG, NT_GLOBAL, NQ_COUNTER, NQ
    INTEGER :: NLM_LMMAXC
    INTEGER :: ISTAT, IERROR
    GDEF, ALLOCATABLE :: GWORK(:)
#ifdef RPAgamma
    REAL(q), ALLOCATABLE :: GWORK_REAL(:,:)
#else
    GDEF, ALLOCATABLE :: GWORK_REAL(:,:)
#endif
    GDEF, ALLOCATABLE, TARGET :: CRHOLM(:) ! augmentation occupancy matrix
    GDEF, ALLOCATABLE, TARGET :: P_NLM(:,:)
    GDEF, ALLOCATABLE :: P_G_PROJ(:,:,:,:)
    GDEF, ALLOCATABLE, TARGET :: P_PROJ_G(:,:,:,:)
    GDEF, ALLOCATABLE :: P_G_R(:,:,:)
    GDEF, ALLOCATABLE :: P_R_G(:,:,:)
    GDEF, ALLOCATABLE :: D_G_PROJ(:,:)
#ifdef RPAgamma
    REAL(q), ALLOCATABLE :: D_R_PROJ(:,:)
#else
    GDEF, ALLOCATABLE :: D_R_PROJ(:,:)
#endif
    GDEF, ALLOCATABLE, TARGET :: D_PROJ_PROJ(:,:)

    GDEF, ALLOCATABLE :: GO_GK_PROJ(:,:)
    GDEF, ALLOCATABLE :: SIG_R_PROJ(:,:)
    GDEF, ALLOCATABLE :: SIG_PROJ_PROJ(:,:)
    GDEF, ALLOCATABLE :: GO_PROJ_PROJ(:,:)
    GDEF, ALLOCATABLE :: GU_GK_PROJ(:,:)
    GDEF, ALLOCATABLE :: GU_R_PROJ(:,:)
    GDEF, ALLOCATABLE :: GU_PROJ_PROJ(:,:)
    GDEF, ALLOCATABLE :: GU_PROJ_R(:,:)

    GDEF, ALLOCATABLE :: P_TMP(:), P_TMPA(:)

    TYPE (wavefun1) :: W1, W2, WQ, WTMP
    LOGICAL LPHASE
    COMPLEX(q)      :: CPHASE(MAX(WGW%GRID%RL%NP,W%WDES%GRID%RL%NP))
    GDEF            :: CDIJ (LMDIM,LMDIM,S%WDES%NIONS,S%WDES%NCDIJ)
!force
    INTEGER         :: NDIR, IDIR
    TYPE (nonlr_struct), ALLOCATABLE :: FAST_AUG(:)
!force

    PROFILING_START('calculate_sigma_super')

    !issues: we don't have real space structures for the response function,
    ! we rely on P_G_R or GWORK or so
    ! this makes FFT_G_SUPER awkward
    ! also COLLECT_ORBITALS_SUPER needs to be changed so that it makes sense
    ALLOCATE(CRHOLM(S%AUG_DES%NPRO))  ! CRHOLM total number of channels in supercell
    CRHOLM=zero
!force
    ! to seek further changes seek IDIR
    IF (LFORCE) THEN
       NDIR=4
    ELSE
       NDIR=1
    ENDIF
    CALL SETUP_FAST_AUG_NDIR

!force
    CALL NEWWAV(W1, S%WDES1,.TRUE.)
    CALL NEWWAV(W2, S%WDES1,.TRUE.)
    CALL NEWWAV(WQ, S%WGW1,.TRUE.)   ! we need to collect the vectors in the supercell

    DO NK=1,W%WDES%NKPTS
       ! FFT and transpose orbitals at all k-points to obtain G%G_r and G%PROJ_r
       ! FFT GU to real space and set entries G'_r, PROJ'_r
       CALL FFT_G_SUPER(W%WDES, GU(NK), GDES, NK)
    ENDDO

    ALLOCATE(P_G_R    (GDES%RES_NRPLWV_ROW_DATA_POINTS, GDES%MPLWV_COL, W%WDES%NKPTS),STAT=ISTAT)  !for W
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "CALCUALTE_SIGMA_SUPER (P_G_R) is not able to allocate "//&
       str( ndata * 8._q*GDES%RES_NRPLWV_ROW_DATA_POINTS*GDES%MPLWV_COL*W%WDES%NKPTS/1024)  // &
       " kB of data on MPI rank 0.")
    ENDIF
    ALLOCATE(P_R_G    (GDES%MPLWV_ROW, GDES%RES_NRPLWV_COL_DATA_POINTS, W%WDES%NKPTS),STAT=ISTAT)  !for W, response
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "CALCUALTE_SIGMA_SUPER (P_R_G) is not able to allocate "//&
       str( ndata * 8._q*GDES%MPLWV_ROW*GDES%RES_NRPLWV_COL_DATA_POINTS*W%WDES%NKPTS/1024)  // &
       " kB of data on MPI rank 0.")
    ENDIF
    CALL REGISTER_ALLOCATE(wsgf*SIZE(P_R_G,KIND=qi8), "GGwork")
    CALL REGISTER_ALLOCATE(wsgf*SIZE(P_G_R,KIND=qi8), "GGwork")
!force
    ALLOCATE(P_PROJ_G    (GDES%NLM_ROW, GDES%RES_NRPLWV_COL_DATA_POINTS, W%WDES%NKPTS,NDIR),STAT=ISTAT) !for D_k(alpha beta, g)
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "CALCUALTE_SIGMA_SUPER (P_PROJ_G) is not able to allocate "//&
       str( ndata * 8._q*GDES%NLM_ROW*GDES%RES_NRPLWV_COL_DATA_POINTS*W%WDES%NKPTS*NDIR/1024)  // &
       " kB of data on MPI rank 0.")
    ENDIF
    ALLOCATE(P_G_PROJ    (GDES%RES_NRPLWV_ROW_DATA_POINTS, GDES%NLM_COL, W%WDES%NKPTS,NDIR),STAT=ISTAT) !transpose
    IF ( ISTAT/=0 ) THEN
       CALL VTUTOR%ERROR( "CALCUALTE_SIGMA_SUPER (P_G_PROJ) is not able to allocate "//&
       str( ndata * 8._q*GDES%RES_NRPLWV_ROW_DATA_POINTS*GDES%NLM_COL*W%WDES%NKPTS*NDIR/1024)  // &
       " kB of data on MPI rank 0.")
    ENDIF
!force
    CALL REGISTER_ALLOCATE(wsgf*SIZE(P_PROJ_G,KIND=qi8), "GGwork")
    CALL REGISTER_ALLOCATE(wsgf*SIZE(P_G_PROJ,KIND=qi8), "GGwork")
    CALL DUMP_ALLOCATE_TAG(IU6,"SIGMA_SUPER"); IF (IU6>=0) CALL WFORCE(IU6)

    P_G_R=zero
    P_R_G=zero
    P_PROJ_G=zero
    P_G_PROJ=zero
!=======================================================================
! FFT of W(G,G') along first direction 
! yielding W(r,G') and W(NLM, G')
! then transpose and conjugate to obtain W(G',r ) and W(G', NLM)
!=======================================================================
    PROFILING_START('calculate_sigma_super1')



    DO NK=1,W%WDES%NKPTS
       CALL SETWDES(WGW, WDESQ, NK)
       !phase e^ikr to add Bloch factor
       CALL SETPHASE_NOCHK(W%WDES%VKPT(1:3,NK), WDESQ%GRID, CPHASE, LPHASE)
       ! sets CREEXP, same k-point, 
       CALL PHASER_HF(GRIDHF, LATT_CUR, FAST_AUG_FOCK, W%WDES%VKPT(1:3,NK))

       ALLOCATE(GWORK(m_ WDESQ%GRID%MPLWV))
       GWORK=zero
    
       !this is a fourier transform for every g' to go g->r, 
       ! response sized quantity
       DO NRP=1,GDES%RES_NRPLWV_COL_DATA_POINTS
          ! FFT of W_k(g,g') for all k-points and set W_k(r,g')
          !need to use P_R_G(NK) array for this and standard FFT
          CALL FFTWAV(WDESQ%NGVECTOR,WDESQ%NINDPW(1), &
                GWORK(1), &
                CHI(NK)%RESPONSEFUN(1, NRP, 1), WDESQ%GRID )
          ! phase factor added afterwards with APPLY_PHASE_ALPHA
          IF (S%WDES1%LOVERL) THEN
             AUG_DES%RINPL=1._q/WDESQ%GRID%NPLWV
             DO IDIR=1,NDIR
                ! int d^3 r Q_LM(r) W(r,G)
                WTMP%CPROJ=> P_PROJ_G(:, NRP, NK, IDIR)   ! RPRO1_HF only dereferences WTMP%CPROJ
                CALL RPRO1_HF(FAST_AUG(IDIR), AUG_DES, WTMP, GWORK(1))
             ENDDO
          ENDIF

          !add the Bloch factor for real space supercell quantity
          IF (LPHASE) CALL APPLY_PHASE_INPLACE( WDESQ%GRID, CPHASE, GWORK(1))
          P_R_G(1:WDESQ%GRID%RL%NP,NRP,NK)=GWORK(1:WDESQ%GRID%RL%NP) ! or (1._q/WDESQ%GRID%NPLWV)
       ENDDO
       DEALLOCATE(GWORK)

       CALL TRANSPOSE_R_G_RESPONSE(P_R_G(:,:,NK), P_G_R(:,:,NK), GDES)
       ! the routine below uses RES_ data descriptors
       ! and NLM array sizes, that's what we need
       IF (S%WDES1%LOVERL) THEN
          DO IDIR=1,NDIR
             CALL TRANSPOSE_PROJ_G_RESPONSE(P_PROJ_G(:,:,NK,IDIR), P_G_PROJ(:,:,NK,IDIR), GDES)
          ENDDO
       ENDIF

       !the phase factor e^ikR_n is added later in COLLECT_D_PROJECTORS_SUPER call
    ENDDO
    CALL DEREGISTER_ALLOCATE(wsgf*SIZE(P_R_G,KIND=qi8), "GGwork")
    CALL DEREGISTER_ALLOCATE(wsgf*SIZE(P_PROJ_G,KIND=qi8), "GGwork")
    DEALLOCATE(P_R_G)
    DEALLOCATE(P_PROJ_G)

    PROFILING_STOP('calculate_sigma_super1')
!=======================================================================
! memory allocation
!=======================================================================
    PROFILING_START('calculate_sigma_super2')

    ALLOCATE(GWORK(m_ S%WDES%GRID%MPLWV), P_TMP(m_ S%WDES%GRID%RC%NP))  ! double allocation, declared as GDEF
    GWORK=zero
    CALL REGISTER_ALLOCATE(wsgf*SIZE(GWORK,KIND=qi8), "GGwork")
    CALL REGISTER_ALLOCATE(wsgf*SIZE(P_TMP,KIND=qi8), "GGwork")

    DO NK=1, NQMAX ! SIGMA is allocated only in the IBZ
       CALL ALLOCATE_G_R(GDES, SIGMA(NK), NK)
       IF (ASSOCIATED(SIGMA(NK)%PROJ_PROJ)) THEN
          CALL vtutor%bug("internal error in CALCULATE_SIGMA_SUPER: SIGMA%PROJ_PROJ is associated", __FILE__, __LINE__)
       ENDIF
       IF (ASSOCIATED(SIGMA(NK)%G_PROJ)) THEN
          CALL vtutor%bug("internal error in CALCULATE_SIGMA_SUPER: SIGMA%G_PROJ is associated", __FILE__, __LINE__)
       ENDIF
       CALL ALLOCATE_G_PROJ(GDES, SIGMA(NK))
       CALL ALLOCATE_PROJ_PROJ(GDES, SIGMA( NK))
       SIGMA( NK)%G_PROJ=0
       SIGMA( NK)%PROJ_PROJ=0       
    ENDDO
    ALLOCATE(GU_PROJ_R(S%WDES1%NPRO,1))

    CALL DUMP_ALLOCATE_TAG(IU6,"SIGMA_SUPER"); IF (IU6>=0) CALL WFORCE(IU6)
!=======================================================================
! contract over real space grid points
! use P_G_R (NK) first
! i) expand it to super cell and contract with G to get SIGMA_R_R
! ii) project to obtain P_SPROJ_R (in super cell)  and contract with G_SPROJ_R to get SIGMA_SPROJ_R
! I think we can save some memory by doing FFT and redistribution after this (avoiding FFT_SIGMA)
!=======================================================================
#ifdef RPAgamma
    ALLOCATE(GWORK_REAL(2* S%WGW%GRID%MPLWV,1))
    GWORK_REAL=0._q
#else
    ALLOCATE(GWORK_REAL(m_ S%WGW%GRID%MPLWV,1))
    GWORK_REAL=zero
#endif
    CALL REGISTER_ALLOCATE(wsgf*SIZE(GWORK_REAL,KIND=qi8), "GGwork")

    ! THIS PART WORKS OK, both (R,r) and (PROJ,r) contraction
    ! loop over all local columns in real space
    DO NRP=1,GDES%MPLWV_COL
       ! transform from storage over k-points (g'+k) to supercell reciprocal space (G')
       ! for green's function
       CALL COLLECT_ORBITALS_SUPER(S, W%WDES, W1%CW, W1%CPROJ, GU(:), NRP)
       ! for screening, use chi descriptor
       CALL COLLECT_RESPONSE_SUPER(S, WGW, WQ%CW, P_G_R, NRP)

       ! FFT supercell index to real space R', Green's function
       CALL FFTWAV(S%WDES1%NGVECTOR, S%WDES1%NINDPW(1), &
            W1%CR(1), W1%CW(1), S%WDES1%GRID)
       ! FFT of screening to real space R' (result is a real valued vector for RPAgamma)
       CALL FFTWAV(S%WGW1%NGVECTOR, S%WGW1%NINDPW(1), &
            GWORK_REAL(1,1), WQ%CW(1), S%WGW1%GRID)  ! gKnew: was previously WQ%CR(1)
!gKnew:       GWORK(1:S%WDES1%GRID%RL%NP)=W1%CR(1:S%WDES1%GRID%RL%NP)*WQ%CR(1:S%WDES1%GRID%RL%NP)
!      complex = complex * complex  |   complex=complex * real, only real part matters (for RPAgamma)
#ifndef RPAgamma
       GWORK(1:S%WDES1%GRID%RL%NP)=W1%CR(1:S%WDES1%GRID%RL%NP)*GWORK_REAL(1:S%WDES1%GRID%RL%NP,1)
#else
       GWORK(1:S%WDES1%GRID%RL%NP)=REAL(W1%CR(1:S%WDES1%GRID%RL%NP),q)*GWORK_REAL(1:S%WDES1%GRID%RL%NP,1)
#endif
       ! FFT to  reciprocal space using supercell FFT: P(g+k)= sum_g+k e -i(k+g) r P(r)
       ! GWORK contains R' data, P_TMP contains G' data
       ! use orbital descriptor since \Sigma has the structure of G
       CALL FFTEXT(S%WDES1%NGVECTOR,S%WDES1%NINDPW(1), &
            GWORK(1), &
            P_TMP(1),S%WDES1%GRID,.FALSE.)
            
       ! distribute the coefficients to k-point structure P_G_R, from G' to g'+k
       ! P_G_R can be SIGMA?
       CALL DISTRIBUTE_SIGMA_SUPER(S, W%WDES, P_TMP, SIGMA(:), NRP, NQMAX, (1.0_q/S%GRID%NPLWV) )
    ENDDO

    CALL DEREGISTER_ALLOCATE(wsgf*SIZE(GWORK_REAL,KIND=qi8), "GGwork")
    DEALLOCATE(GWORK_REAL)
    DEALLOCATE(P_TMP, GU_PROJ_R)
    CALL DEREGISTER_ALLOCATE(wsgf*SIZE(P_G_R,KIND=qi8), "GGwork")
    CALL DEREGISTER_ALLOCATE(wsgf*SIZE(P_TMP,KIND=qi8), "GGwork")
    DEALLOCATE(P_G_R)

    PROFILING_STOP('calculate_sigma_super2')
    !at this point we have SIGMA(1,:)%G'_r and SIGMA(1,:)%PROJ'_r
    !transpose back and do fft for r is done in FFT_SIGMA
!=======================================================================
! now the akward part dealing with P_G_PROJ =  W(G', NLM)
!
! i)  loop over atoms and projectors and expand in super cell to P_SR_PROJ    [SR like Supercell R]
! ii) project to get P_SPROJ_PROJ
!     contract with G, distribute, and backtransform
! not sure which order will work the best
!=======================================================================

    !start with the SR_PROJ part first
    ! we have P_G_PROJ(:,:,NK), where the second index is in the NLM basis
    ! we also need to prepare GWORK(m_ S%WGW1%GRID%MPLWV) for collecting data and doing FFTs?
    ! and GU_R_PROJ(S%WDES1%GRID%RL%NP, LMMAXC)  for the collected Green's function
    ! where we will store G (R', \alpha)
    ! both of them then need to be contracted to \Sigma(R',\beta)
    IF (S%WDES1%LOVERL) THEN
    PROFILING_START('calculate_sigma_super3')

    dir: DO IDIR=1,NDIR
    ! contract over beta, beta' to N'L'M'
    LMBASE =0
    NAUG=0
    NIS =1
    DO NT=1,GDES%NTYP  ! note NT is a local type index in the primitive cell
                       ! and might not span all atomic types
       NT_GLOBAL=GDES%NT_GLOBAL(NT)  ! global type index

       LMMAXC=GDES%LMMAX(NT)
       NLM_LMMAXC=GDES%NLM_LMMAX(NT)
       IF (LMMAXC.NE.0) THEN
          ! WGW1%GRID%MPLWV is the number of complex words required to do in place FFT for \chi
          ! sized routines
          ALLOCATE(GU_R_PROJ(S%WDES1%GRID%RL%NP, LMMAXC) )   ! stores orbital (r', alpha)
          ALLOCATE(SIG_R_PROJ(S%WDES1%GRID%RL%NP, LMMAXC))   ! stores sigma (r', beta) for one atom
          ALLOCATE(SIG_PROJ_PROJ(S%WDES1%NPRO, LMMAXC))      ! stores sigma (beta, beta) for one atom
          ALLOCATE(GU_GK_PROJ(S%WDES1%NGVECTOR ,LMMAXC))     ! stores orbital(g+k, alpha)
          ALLOCATE(GU_PROJ_PROJ(S%WDES1%NPRO ,LMMAXC))       ! stores orbital(g+k, alpha)
          ALLOCATE(D_G_PROJ(S%WGW1%NGVECTOR, NLM_LMMAXC))    ! stores W (g'+k, alpha beta)
#ifdef RPAgamma
          ALLOCATE(D_R_PROJ(S%WGW%GRID%MPLWV*2, NLM_LMMAXC)) ! stores W (r', alpha beta) & in place FFT
#else
          ALLOCATE(D_R_PROJ(S%WGW%GRID%MPLWV, NLM_LMMAXC))   ! stores W (r', alpha beta) & in place FFT
#endif
          ALLOCATE(D_PROJ_PROJ(S%AUG_DES%NPRO ,NLM_LMMAXC))  ! stores W(alpha beta, alpha beta)

       
          ! main loop over local ion index N (second = column index)
          ! this index is in the primitive cell
          ! all contributions on one ion are treated at the same time
          ! this loop goes only over the unit cell as it is the first index (r)
          ! the other index is the supercell index 
          DO NI=NIS,GDES%NITYP(NT)+NIS-1
             SIG_R_PROJ=zero
             SIG_PROJ_PROJ=zero
             GU_GK_PROJ=zero
             GU_R_PROJ=zero
             GU_PROJ_PROJ=zero
             D_G_PROJ=zero
             D_R_PROJ=zero
             D_PROJ_PROJ=zero

             !transform to supercell in the first index of G(R',lm)
             !phase factor applied
             CALL COLLECT_PROJECTORS_SUPER(S, W%WDES, GU_GK_PROJ, GU_PROJ_PROJ, GU(:), &
                     LMBASE, LMMAXC, GDES%POSION(:,NI))
             DO  L=1, LMMAXC
                W1%CW(1:S%WDES1%NGVECTOR)=GU_GK_PROJ(1:S%WDES1%NGVECTOR, L)
                CALL FFTWAV(S%WDES1%NGVECTOR, S%WDES1%NINDPW(1), &
                        W1%CR(1),W1%CW(1),S%WDES1%GRID) ! GO_GK_PROJ(1, L),S%WDES1%GRID)
                GU_R_PROJ(1:S%WDES1%GRID%RL%NP, L)=W1%CR(1:S%WDES1%GRID%RL%NP)
             ENDDO !  L=1, LMMAXC
             !transform to supercell for P_G_PROJ(:,:,NK), similar, but with LM instead of lm
             !from P_G_PROJ(~NGVECTOR, ~NLM_NPRO, NK) to DWORK2(m_ S%WGW1%GRID%RC%NP, NLM_LMMAX)
             !all can be done in one loop or separated
             CALL COLLECT_D_PROJECTORS_SUPER(S, WGW, D_G_PROJ, P_G_PROJ(:,:,:,IDIR), &
                     NAUG, NLM_LMMAXC, GDES%POSION(:,NI),.TRUE.)

             !and transform to real space in supercell
             DO LM=1, NLM_LMMAXC
                CALL FFTWAV(S%WGW1%NGVECTOR,S%WGW1%NINDPW(1), &
                        D_R_PROJ(1,LM), &
                        D_G_PROJ(1,LM),S%WGW1%GRID)
                !project 
                S%AUG_DES%RINPL=1.0_q/S%WGW1%GRID%NPLWV   !AUG_DES or S%AUG_DES ???
                WTMP%CPROJ=>D_PROJ_PROJ(:,LM)
#ifdef RPAgamma
                CALL RPRO1_REAL(S%FAST_AUG, S%AUG_DES, WTMP, D_R_PROJ(:,LM))
#else
                CALL RPRO1_HF(S%FAST_AUG, S%AUG_DES, WTMP, D_R_PROJ(:,LM))
#endif
             ENDDO ! LM=1, NLM_LMMAXC
             !loop over index for GU
             DO L=1,LMMAXC
             DO LP=1,LMMAXC
                CRHOLM=zero
                DO LM=1,NLM_LMMAXC  
                   !transform second index  D_PROJ_PROJ(LM,L'M')->CRHOLM(LM) for alpha' beta'
                   CRHOLM(:)=CRHOLM(:)+D_PROJ_PROJ(:,LM)*TRANS_MATRIX_FOCK(LP,L,LM,NT_GLOBAL)
                ENDDO
                !transform over first index CRHOLM(NLM)->CDIJ(alpha beta)
                CALL CALC_DLLMM_TRANS(S%WDES,S%AUG_DES,S%TRANS_MATRIX, CDIJ(:,:,:,:), CRHOLM)

                !sum over alpha,
                CALL OVERL_FOCK(S%WDES, SIZE(CDIJ,1), CDIJ(1,1,1,1), GU_PROJ_PROJ(:,L), SIG_PROJ_PROJ(:,LP), .TRUE.)
                !W2%CPROJ now has contributions to Sigma%PROJ_PROJ(SPROJ,ion)
                !this needs to be redistributed in the first index to 
                !SIGMA(NK)%PROJ_PROJ(:,LMBASE+LP)
             ENDDO !LP
             ENDDO !L

#ifdef RPAgamma
             CALL CALC_DLLMM_VECTOR_REAL( SIG_R_PROJ, GU_R_PROJ, D_R_PROJ, S%WDES1%GRID%RL%NP, TRANS_MATRIX_FOCK(:,:,:,NT_GLOBAL), LMMAXC, 0,  NLM_LMMAXC)
#else
             CALL CALC_DLLMM_VECTOR( SIG_R_PROJ, GU_R_PROJ, D_R_PROJ, S%WDES1%GRID%RL%NP, TRANS_MATRIX_FOCK(:,:,:,NT_GLOBAL), LMMAXC, 0,  NLM_LMMAXC)
#endif

             !distribute SIGMA 
             DO LP=1,LMMAXC
                CALL FFTEXT(S%WDES1%NGVECTOR,S%WDES1%NINDPW(1), &
                         SIG_R_PROJ(1,LP), GWORK(1), &
                         S%WDES1%GRID,.FALSE.)
                ! set SIGMA%PROJ_PROJ and SIGMA%G_PROJ = GWORK
                CALL DISTRIBUTE_SIGMA_PROJ_SUPER(S, W%WDES, GWORK, SIG_PROJ_PROJ(:,LP), SIGMA(:), &
                                 LMBASE+LP, NQMAX, (1.0_q/S%GRID%NPLWV) )
             ENDDO !LP
             LMBASE = LMMAXC+LMBASE
             NAUG = NAUG + GDES%NLM_LMMAX(NT)

          ENDDO !NI=NIS,GDES%NITYP(NT)+NIS-1

          DEALLOCATE(GU_R_PROJ) ! stores orbital (r', alpha)
          DEALLOCATE(GU_GK_PROJ,GU_PROJ_PROJ)   ! stores orbital(g+k, alpha)
          DEALLOCATE(D_G_PROJ, SIG_R_PROJ, D_R_PROJ, D_PROJ_PROJ, SIG_PROJ_PROJ)

       ENDIF !LMMAXC.NE.O
       NIS = NIS+GDES%NITYP(NT)
    ENDDO !NT=1,GDES%NTYP

    DO NQ=1,NQMAX
       CALL SETWDES(W%WDES, WDESQ, NQ )
       CALL APPLY_PHASE_ALPHA(W%WDES, GDES, SIGMA( NQ)%G_PROJ, +WDESQ%VKPT(1:3), NQ)  
       CALL APPLY_PHASE_ALPHA(W%WDES, GDES, SIGMA( NQ)%PROJ_PROJ, +WDESQ%VKPT(1:3), NQ)

       ! get non-local contributions to the energy
       IF (LFORCE) THEN
          CALL CONTRACT_SIGMA_G_TAU_NONL( WDESQ, KPOINTS_ORIG%WTKPT(NQ)*TAU_WEIGHT, & 
               GO(NQ), SIGMA(NQ), GDES, FRNL, IDIR )
       ENDIF
    ENDDO


    ENDDO dir
    PROFILING_STOP('calculate_sigma_super3')

   
    ENDIF !S%WDES1%LOVERL

    DEALLOCATE (GWORK)
    CALL DEREGISTER_ALLOCATE(wsgf*SIZE(GWORK,KIND=qi8), "GGwork")
    DEALLOCATE(P_G_PROJ)
    CALL DEREGISTER_ALLOCATE(wsgf*SIZE(P_G_PROJ,KIND=qi8), "GGwork")

    DO NK=1,W%WDES%NKPTS
       ! release G_R, PROJ_R, R_PROJ
       CALL RELEASE_FFT_G( GU(NK))
    ENDDO

    PROFILING_START('calculate_sigma_super4')
    DO NQ=1,NQMAX
       CALL SETWDES(W%WDES, WDESQ, NQ )
       ! set exp(-i q r)
       !minus works with pluses above
       CALL SETPHASE_NOCHK(-WDESQ%VKPT(1:3), WDESQ%GRID, CPHASE, LPHASE)
       ! transpose sigma and fft in other direction to gain Sigma(g,g') for each k
       CALL FFT_SIGMA_SUPER ( W%WDES, SIGMA( NQ), GDES, NQ, CPHASE, LPHASE )

! if only half of the coefficients are available in second direction
! uncompress using inversion symmetry
#ifdef RPAgamma
       CALL UNCOMPRESS_G_RECIPROCAL( SIGMA(NQ), GDES, NQ)
#endif

    ENDDO !NK loop
    PROFILING_STOP('calculate_sigma_super4')

    CALL DELWAV(W1, .TRUE.)
    CALL DELWAV(W2, .TRUE.)
    CALL DELWAV(WQ, .TRUE.)

    CALL FREE_FAST_AUG_NDIR()

    PROFILING_STOP('calculate_sigma_super')

    CONTAINS

!==========================================================================
! setup first derivative of augmentation charges with respect
! to ionic positions
!==========================================================================

      SUBROUTINE SETUP_FAST_AUG_NDIR
        INTEGER :: IDIR
        REAL(q) :: DIS=fd_displacement
        
        REAL(q)    :: DISPL1(3,FAST_AUG_FOCK%NIONS),DISPL2(3,FAST_AUG_FOCK%NIONS)

        IF (.NOT.S%WDES1%LOVERL) RETURN

        ALLOCATE(FAST_AUG(NDIR))
        ! last entry is identical to FAST_AUG_FOCK
!       FAST_AUG(NDIR)=FAST_AUG_FOCK
        CALL NONLR_ASSIGN(FAST_AUG(NDIR),FAST_AUG_FOCK)

        ! set first derivative of in first 3 entries
        DO IDIR=1,NDIR-1
           CALL COPY_FASTAUG(FAST_AUG_FOCK, FAST_AUG(IDIR))
           DISPL1=0
           DISPL1(IDIR,:)=-DIS
           
           DISPL2=0
           DISPL2(IDIR,:)=+DIS
           CALL RSPHER_ALL(GRIDHF, FAST_AUG(IDIR),LATT_CUR,LATT_CUR,LATT_CUR, DISPL1, DISPL2,1)
           ! phase factor identical to FAST_AUG_FOCK
           IF (ASSOCIATED(FAST_AUG(IDIR)%CRREXP)) DEALLOCATE(FAST_AUG(IDIR)%CRREXP)
           FAST_AUG(IDIR)%CRREXP=>FAST_AUG_FOCK%CRREXP
#ifdef shmem_rproj
           IF (FAST_AUG(IDIR)%COMM_shmem%NODE_ME==1) FAST_AUG(IDIR)%RPROJ=FAST_AUG(IDIR)%RPROJ*SQRT(LATT_CUR%OMEGA)*(1._q/(2*DIS))
           CALLMPI( M_barrier(FAST_AUG(IDIR)%COMM_shmem))
#else
           FAST_AUG(IDIR)%RPROJ=FAST_AUG(IDIR)%RPROJ*SQRT(LATT_CUR%OMEGA)*(1._q/(2*DIS))
#endif
        ENDDO
      END SUBROUTINE SETUP_FAST_AUG_NDIR
      
      SUBROUTINE FREE_FAST_AUG_NDIR
        IF (.NOT.S%WDES1%LOVERL) RETURN

        DO IDIR=1,NDIR-1
           NULLIFY(FAST_AUG(IDIR)%CRREXP)
           CALL NONLR_DEALLOC(FAST_AUG(IDIR))
        ENDDO
        DEALLOCATE(FAST_AUG)
      END SUBROUTINE FREE_FAST_AUG_NDIR

  END SUBROUTINE CALCULATE_SIGMA_SUPER


!=======================================================================
! medium level routines used in cubic GW schedulers
!=======================================================================

!***********************************************************************
!
! switch off symmetry
!
!***********************************************************************

  SUBROUTINE NO_SYMMETRY(LSUPER, LGAMMA, SYMM, W, WDES, LATT_CUR, LATT_INI, &
     T_INFO, INFO, KPOINTS, GRID, NONL_S, P, WHF, IO )
    USE msymmetry, ONLY: NOSYMM
    USE pseudo_struct_def 
    LOGICAL            LSUPER
    LOGICAL            LGAMMA
    TYPE (symmetry)    SYMM
    TYPE (wavespin)    W
    TYPE (wavedes)     WDES
    TYPE (latt)        LATT_CUR, LATT_INI
    TYPE (type_info)   T_INFO
    TYPE (info_struct) INFO
    TYPE (kpoints_struct) KPOINTS
    TYPE (grid_3d)     GRID       ! grid for wavefunctions
    TYPE (nonl_struct) NONL_S
    TYPE (potcar)      P(T_INFO%NTYP)
    TYPE (wavespin)    WHF
    TYPE (in_struct)   IO
    !local
    LOGICAL            LRE_READ

    LRE_READ=.FALSE.
#ifdef RPAgamma
    IF ( LSUPER ) THEN
       ! used in supercell approach
       LRE_READ=SYMM%ISYM>=0.AND..NOT.( WDES%LNONCOLLINEAR .OR. LCRPA)
    ENDIF
#endif
    IF (SYMM%ISYM>=0 .AND. .NOT. LGAMMA) THEN
       ! switch of symmetry
       CALL NOSYMM(LATT_CUR%A,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS,SYMM%PTRANS, &
            SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%MAGROT,INFO%ISPIN,IO%IU6)
       CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
            LRE_READ,&
            T_INFO%NIONS,SYMM%ROTMAP, SYMM%MAGROT, SYMM%ISYM, -1,IO%IU0)
       CALL KPAR_SYNC_ALL(WDES,W)
       CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
       CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
       
       WHF=W
       WHF%WDES => WDES_FOCK
       ! presumably W from file
       !IF (LWFROMFILE) &
       !     CALL REALLOCATE_WAVE( W_W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)

       IF (IO%IU6>=0) WRITE(IO%IU6,'(A)') ' orbital symmetry switched off'

    ENDIF

  END SUBROUTINE NO_SYMMETRY


!***********************************************************************
!
! switch on symmetry
!
!***********************************************************************
  SUBROUTINE ON_SYMMETRY( LGAMMA, SYMM, W, WDES, LATT_CUR, LATT_INI, &
     T_INFO, DYN, INFO, KPOINTS, GRID, NONL_S, P, WHF, IO )
    USE pead, ONLY : LPEAD_SYM_RED
    USE msymmetry, ONLY : INISYM
    USE pseudo_struct_def 
    LOGICAL            LGAMMA
    TYPE (symmetry)    SYMM
    TYPE (wavespin)    W
    TYPE (wavedes)     WDES
    TYPE (latt)        LATT_CUR, LATT_INI
    TYPE (type_info)   T_INFO
    TYPE (dynamics)    DYN
    TYPE (info_struct) INFO
    TYPE (kpoints_struct) KPOINTS
    TYPE (grid_3d)     GRID       ! grid for wavefunctions
    TYPE (nonl_struct) NONL_S
    TYPE (potcar)      P(T_INFO%NTYP)
    TYPE (wavespin)    WHF
    TYPE (in_struct)   IO

    IF (SYMM%ISYM>=0.AND. .NOT. LGAMMA) THEN
       IF (SYMM%ISYM>0) &
       CALL INISYM(LATT_CUR%A,T_INFO%POSION,DYN%VEL,T_INFO%LSFOR, &
            T_INFO%LSDYN,LPEAD_SYM_RED(),T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS, &
            SYMM%PTRANS,SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%TAU,SYMM%TAUROT,SYMM%WRKROT, &
            SYMM%INDROT,T_INFO%ATOMOM,WDES%SAXIS,SYMM%MAGROT,WDES%NCDIJ,IO%IU6)
       CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR, &
            SYMM%ISYM>=0.AND..NOT.( WDES%LNONCOLLINEAR .OR. LCRPA ), &
            T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,-1,IO%IU0)
       CALL KPAR_SYNC_ALL(WDES,W)
       CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_INI, IO%IU6, IO%IU0)
       CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
       IF (IO%IU6>=0) WRITE(IO%IU6,'(A)') ' orbital symmetry switched on'

       WHF=W
       WHF%WDES => WDES_FOCK

    ENDIF
  END SUBROUTINE ON_SYMMETRY

!*******************************************************************
!
! screened Coulomb potential W
! the potential is stored distributed with B%NPOINTS_IN_GROUP
! frequencies per group
! the frequency in the array OMEGA can be determined as
! OMEGA(B%COMM_BETWEEN_GROUPS%NODE_ME+ (NO%POINTS_LOCAL-1)* B%COMM_BETWEEN_GROUPS%NCPU)
!
!*******************************************************************
   SUBROUTINE SCREENED_POTENTIAL(S2E,CHI_QPT, WGW, WGWQ, &
      WDES, GDES, LATT_CUR, FSG0, IMAG_GRIDS, LDUMP_QPT, IO )
      USE wpot, ONLY: WRITE_WPOT_FULL
      TYPE (screened_2e_handle) :: S2E
      TYPE (responsefunction), POINTER :: CHI_QPT(:)
      TYPE (wavedes), POINTER :: WGW
      TYPE (wavedes1)         :: WGWQ
      TYPE (wavedes)          :: WDES
      TYPE (greensfdes)       :: GDES     
      TYPE (latt)             :: LATT_CUR
      REAL(q)                 :: FSG0
      TYPE (imag_grid_handle) :: IMAG_GRIDS  
      LOGICAL                 :: LDUMP_QPT 
      TYPE (in_struct)        :: IO
      ! local
      TYPE (responsefunction) :: WPOT_DUMP  
      INTEGER                 :: NQ_COUNTER, NQ
      INTEGER                 :: NOMEGA_LOCAL 
      INTEGER                 :: NOMEGAPAR_ID
      TYPE (responsefunction), POINTER :: CHI
      CHARACTER(LEN=16)       :: TEXT

      DO NQ_COUNTER=1,S2E%NUMBER_OF_NQ  !this is over q-point (symmetry taken into account) 
         NQ=S2E%NQ(NQ_COUNTER)
         ! this is the response function in frequency space
         CHI=>CHI_QPT(NQ_COUNTER)
         CALL SETWDES(WGW, WGWQ, NQ )
  
         CALL COMPUTE_SCREENED_POTENTIAL( WDES, CHI, GDES, WGWQ, IMAG_GRIDS, &
         LATT_CUR, NQ, IDIR_MAX, FSG0, IO)
         !dump timing for each q-point separately 
         WRITE(TEXT,*) NQ
         IF ( LDUMP_QPT ) &
             CALL STOP_TIMING("GREENS",IO%IU6,'W_q='//TRIM(ADJUSTL(TEXT)))
         ! 
         ! no need to dump W to WFULL????.tmp by default 
         ! 
         IF ( NOMEGA_DUMP < 1 ) THEN
            CYCLE
         ENDIF
         ! NOMEGA_DUMP is the global frequency point,
         ! get local frequency point and group holding this point 
         NOMEGAPAR_ID = DETERMINE_CURRENT_FREQ_GROUP( NOMEGA_DUMP ) 
         !
         ! in case user selects NOMEGAPAR group with NOMEGA_DUMP and 
         ! and local frequency point of the group with NOMEGAPAR_ID
         ! WFULL????.tmp files are writen containing the corresponding 
         ! screened potential 
         ! 
         ! WPOT_DUMP is the screened potential, which will be written to 
         ! first it has to be gathered from the ranks in the NOMEGAPAR group
         ! holding the correct frequency point
         ! allocate complete response matrix with only one frequency point
         CALL ALLOCATE_RESPONSEFUN(WPOT_DUMP, GDES%RES_NRPLWV_ROW_DATA_POINTS, WGWQ%LGAMMA, WGW%LGAMMA,1)
         ! set inverse temperature 
         WPOT_DUMP%BETA = IMAG_GRIDS%BETA        
         !clear the response array
         CALL CLEAR_RESPONSE(WPOT_DUMP)
         !set current q-point of reponse function
         CALL SET_RESPONSE_KPOINT(WPOT_DUMP, WDES%VKPT(:,NQ), NQ)
         ! pick frequency group with ID NOMEGAPAR_ID
         NOMEGA_LOCAL = 0 
         IF ( NOMEGAPAR_ID == IMAG_GRIDS%B%COMM_BETWEEN_GROUPS%NODE_ME ) THEN
            ! get the local frequency index from NOMEGA_DUMP
            NOMEGA_LOCAL = DETERMINE_LOCAL_FREQ_INDEX( NOMEGA_DUMP ) 
            !collect data 
            IF ( CHI%LREALSTORE ) THEN
               WPOT_DUMP%RESPONSER(1:GDES%RES_NRPLWV_ROW_DATA_POINTS, &
                  GDES%RES_NRPLWV_POS:GDES%RES_NRPLWV_POS+GDES%RES_NRPLWV_COL_DATA_POINTS-1, 1 )= &
                  CHI%RESPONSER(1:GDES%RES_NRPLWV_ROW_DATA_POINTS,1:GDES%RES_NRPLWV_COL_DATA_POINTS,NOMEGA_LOCAL)
               CALLMPI( M_sum_d(IMAG_GRIDS%B%COMM_IN_GROUP,WPOT_DUMP%RESPONSER(1,1,1),GDES%RES_NRPLWV_ROW_DATA_POINTS**2))
            ELSE
               WPOT_DUMP%RESPONSEFUN(1:GDES%RES_NRPLWV_ROW_DATA_POINTS, &
                  GDES%RES_NRPLWV_POS:GDES%RES_NRPLWV_POS+GDES%RES_NRPLWV_COL_DATA_POINTS-1, 1 )= &
                  CHI%RESPONSEFUN(1:GDES%RES_NRPLWV_ROW_DATA_POINTS,1:GDES%RES_NRPLWV_COL_DATA_POINTS,NOMEGA_LOCAL)
               CALLMPI( M_sum_z(IMAG_GRIDS%B%COMM_IN_GROUP,WPOT_DUMP%RESPONSEFUN(1,1,1),GDES%RES_NRPLWV_ROW_DATA_POINTS**2))
            ENDIF
            ! is is set to indicate which rank is going to write the files 
            ! in the WRITE_WPOT_FULL routine 
            IF ( IMAG_GRIDS%B%COMM_IN_GROUP%NODE_ME /= 1 ) NOMEGA_LOCAL = 0
         ENDIF

         !root node in root frequency group is going to write the file
         IF ( NOMEGA_LOCAL > 0 ) THEN
            !set NOMEGA_LOW to 1, otherwise WRITE_WPOT does nothing
            WPOT_DUMP%NOMEGA_LOW = 1  
            IF ( NQ_COUNTER == S2E%NUMBER_OF_NQ ) &
              WRITE(*,'(" WFULL written for iw=",F14.8)')IMAG_GRIDS%B%POINTS_LOCAL( NOMEGA_LOCAL )
         ELSE
            WPOT_DUMP%NOMEGA_LOW = -1
         ENDIF

         !write potential
         CALL WRITE_WPOT_FULL(WPOT_DUMP, WGWQ, FSG0, NQ, LFOCK_SUBTRACT )
 
         !get rid of CHI_TAU
         CALL DEALLOCATE_RESPONSEFUN(WPOT_DUMP)
      ENDDO
   END SUBROUTINE SCREENED_POTENTIAL

!******************************************************************
! 
!  Determine CHI_QPT on imaginary freq. axis
!  includes the limit q->0 correction if requested with ADD_Q0_LIMIT
!  both, supercell and primitive cell approach are used 
!  based on value of ICHIREAL
!
!    basic form of routine is as follows 
!
!    tau loop
!        occupied and unoccupied Green's functions
!    
!        chi_q(i tau) = sum_{k} G_q(it) G_q-k(-it)
!        Contraction+FT: CHI_q(iw) = int dt chi_q(-it) cos(wt)
!    end tau loop
!
!    +long wave limit q->0
!
! on ENTRY: CHI_TAU, CHI_QPT arrays allocated with proper sizes
! 
! on EXIT:  CHI_QPT contains chi(omega)
! 
!*******************************************************************
   SUBROUTINE DETERMINE_CHI_Q_OMEGA( CHI_QPT, CHI_TAU, &
      GO_MAT, GU_MAT, GDES_MAT, &
      SUPER, S2E, W, WDES, WHF, WCORR, WCORRHF, WGW, WGWQ, &
      GDES_TAU, GDES, KPOINTS, NKPTS_IRZ, LATT_CUR, &
      NELM, IMAG_GRIDS, IO, IU6_MEM, &
      ADD_Q0_LIMIT, LDO_CRPA, ONE_ELECTRON_GREENS, CHI0  )
      USE constant
USE c2f_interface, ONLY: TIMING
      TYPE (responsefunction), POINTER :: CHI_QPT(:)
      TYPE (responsefunction), POINTER :: CHI_TAU(:)
      GDEF, POINTER, CONTIGUOUS        :: GO_MAT(:,:,:,:)
      GDEF, POINTER, CONTIGUOUS        :: GU_MAT(:,:,:,:)
      TYPE (greens_mat_des), POINTER   :: GDES_MAT  ! scGF in time
      TYPE(supercell), POINTER         :: SUPER
      TYPE (screened_2e_handle)        :: S2E
      TYPE (wavespin)                  :: W
      TYPE (wavedes)                   :: WDES
      TYPE (wavespin)                  :: WHF
      TYPE (wavespin)                  :: WCORR
      TYPE (wavespin)                  :: WCORRHF 
      TYPE (wavedes), POINTER          :: WGW
      TYPE (wavedes1)                  :: WGWQ
      TYPE (greensfdes)                :: GDES_TAU  ! GF in time
      TYPE (greensfdes)                :: GDES      ! GF in freq. 
      TYPE (kpoints_struct)            :: KPOINTS
      INTEGER                          :: NKPTS_IRZ ! tot # of k-points in the IRZ
      TYPE (latt)                      :: LATT_CUR
      INTEGER                          :: NELM
      TYPE (imag_grid_handle)          :: IMAG_GRIDS
      TYPE (in_struct)                 :: IO
      INTEGER                          :: IU6_MEM
      LOGICAL :: ADD_Q0_LIMIT          ! add q->0 limit?
      LOGICAL :: LDO_CRPA              ! compute CRPA polarizability?
      LOGICAL :: ONE_ELECTRON_GREENS   ! use G_0 or self-consistent G?
      TYPE (responsefunction), POINTER, OPTIONAL :: CHI0(:)
      ! local temporary arrays
      TYPE (responsefunction), POINTER :: CHI=>NULL() ! temporary response function in nu domain
      TYPE( head_handle)               :: HEADHAND    ! handle for the long-wave limit of chi
      TYPE (greensf), POINTER          :: GU(:)=>NULL() ! temp. G for postive times
      TYPE (greensf), POINTER          :: GO(:)=>NULL() ! temp. G for negative times
      ! some loop variables
      INTEGER         :: NTAU_ROOT
      INTEGER         :: ISP
      INTEGER         :: NQ, NQ_COUNTER
      LOGICAL         :: LUSESUPER
      INTEGER         :: I,J
      REAL(q)         :: PSS
#ifdef gammareal
      LUSESUPER=.FALSE.
#else
      LUSESUPER=.TRUE.
#endif

#ifdef mem_profiling
      CALL GET_PSS( PSS, IO%IU0, IO%IU6, "CHI START", 5 )
#endif
      ! chi_GG uses primitive cell apporach
      IF( ICHIREAL == 1 ) LUSESUPER = .FALSE.

      IU6_MEM=IO%IU6

      tau_points_gg: DO NTAU_ROOT=1, IMAG_GRIDS%T%NPOINTS_IN_ROOT_GROUP
         ! set up individual tau indices for each tau group
         CALL SETUP_TAU_INDICES(NTAU_ROOT,IMAG_GRIDS%T)
         DO ISP=1,W%WDES%ISPIN
            ! compute the occupied and unoccupied Green functions at the
            ! current local tau point of the group 
            CALL CALCULATE_G_POSSIBLY_SC( W, WDES, WHF, WCORR, WCORRHF, &
               GDES_TAU, GDES_MAT, NKPTS_IRZ, ISP, IMAG_GRIDS%T, ONE_ELECTRON_GREENS, NTAU_ROOT, & 
               GU, GO, GU_MAT, GO_MAT, IU6_MEM, IO, .FALSE. )
            CALL STOP_TIMING("GREENS",IO%IU6,"GREENS")

#ifdef mem_profiling
            CALL GET_PSS( PSS, IO%IU0, IO%IU6, "G DONE",0 )
#endif
   
            ! calculate X = G G and fourier transform frequency domain
            CALL CHI_QPT_IN_FREQ_DOMAIN()

            !CRPA, subtract correlated polarizability 
            IF ( LDO_CRPA .AND. LCHIC ) THEN
               ! compute the occupied and unoccupied Green functions at the
               ! current local tau point of the group 
               IF( NTAU_ROOT == 1 .AND. ISP == 1 )IU6_MEM = IO%IU6
               CALL CALCULATE_G_POSSIBLY_SC( W, WDES, WHF, WCORR, WCORRHF, &
                  GDES_TAU, GDES_MAT, NKPTS_IRZ, ISP, IMAG_GRIDS%T, ONE_ELECTRON_GREENS, NTAU_ROOT, &
                  GU, GO, GU_MAT, GO_MAT, IU6_MEM, IO , .TRUE. )
               CALL STOP_TIMING("GREENS",IO%IU6,"GREENS")
   
               ! calculate X =-G_c G_c and fourier transform frequency domain
               IMAG_GRIDS%FACTOR = -IMAG_GRIDS%FACTOR
               CALL CHI_QPT_IN_FREQ_DOMAIN()
               IMAG_GRIDS%FACTOR = -IMAG_GRIDS%FACTOR
            ENDIF

         ENDDO ! spin
      ENDDO tau_points_gg

      ! sync the response over KPAR nodes
      CALL SYNC_CHI(WDES, CHI_QPT)


      ! determine long-wave limit of CHI in frequency domain
      ! TODO gK: think about this for scf GW 
      IF ( .NOT. ADD_Q0_LIMIT ) THEN
         RETURN
      ENDIF
      ! calculates q->0 limit of CHI for a set of imaginary frequencies
      ! and adds the result to the RESPONSEFUN structure
      PROFILING_START( 'add_head_and_wings_to_chi_qpt' )

      !set up head and wings of chi in imaginary frequency
      !this needs to be done only for q=0, so first find the corresponding q-point
      DO NQ_COUNTER=1,S2E%NUMBER_OF_NQ
         NQ=S2E%NQ(NQ_COUNTER)
         CHI=>CHI_QPT(NQ_COUNTER)
         CALL SETWDES(WGW, WGWQ, NQ )
         !set local complex frequency arrays in CHI
         CALL SET_RESPONSE_FREQUENCIES(CHI,IMAG_GRIDS)
         IF ( ABS(SUM(CHI%VKPT*CHI%VKPT))<G2ZERO ) EXIT
      ENDDO
  
      IF (ALLOCATED(CDER_BETWEEN_STATES) .OR. EFERMIADD>0 .OR. SUM(ABS(WPLASMON))>0) THEN
         !allocate long wave limit handle  
         CALL ALLOCATE_HEAD_HANDLE( HEADHAND, CHI, IMAG_GRIDS ) 
  
         !calculate long wave limit and save it to HEADHAND
         IF ( LCRPA .AND. LDISENTANGLED ) THEN
            CALL CALCULATE_LONG_WAVE_LIMIT( LATT_CUR, WCORR, WGWQ, KPOINTS, &
            NSTRIP_TOTAL, CHI, HEADHAND, IMAG_GRIDS)
         ELSE
            CALL CALCULATE_LONG_WAVE_LIMIT( LATT_CUR, W, WGWQ, KPOINTS, &
            NSTRIP_TOTAL, CHI, HEADHAND, IMAG_GRIDS)
         ENDIF
         
         !remove correlated head, for projected CRPA only 
         IF ( LCRPA .AND. LPROJECTED ) THEN 
            !calculate long wave limit and save it to HEADHAND
            HEADHAND%FACTOR = -HEADHAND%FACTOR 
            CALL CALCULATE_LONG_WAVE_LIMIT( LATT_CUR, WCORR, WGWQ, KPOINTS, &
            NSTRIP_TOTAL, CHI, HEADHAND, IMAG_GRIDS)
            HEADHAND%FACTOR = -HEADHAND%FACTOR 
         ENDIF
 
         IF ( .NOT. LFINITE_TEMPERATURE ) THEN
            !pL: add Drude contribution to the intra-band transitions
            IF ( .NOT. (LCRPA )  ) THEN
              CALL ADD_DRUDE_TERM(HEADHAND,WPLASMON,LATT_CUR%OMEGA)
            ENDIF 
        
            DO NQ_COUNTER=1,S2E%NUMBER_OF_NQ
               NQ=S2E%NQ(NQ_COUNTER)
               CALL SETWDES(WGW, WGWQ, NQ )
               ! for the time being skip this for CRPA
               ! TODO: CRPA correction should be done inside the routine below
               ! subtract self correlation from (1-f) f terms
               IF ( (.NOT. LCRPA) ) &
               CALL SUBTRACT_SELF_CORRELATION( LATT_CUR, W, WGWQ, NSTRIP_TOTAL, &
                  CHI_QPT(NQ_COUNTER), HEADHAND, IMAG_GRIDS, GDES )
            ENDDO
         ENDIF
         !distribute HEADHAND to CHI
         CALL DISTRIBUTE_HEAD_TO_CHI( HEADHAND, CHI )
  
         !destroy handle 
         CALL DEALLOCATE_HEAD_HANDLE( HEADHAND )
  
         CALL STOP_TIMING("GREENS",IO%IU6,"HEAD")
         !dump progress to stdout and OSZICAR
         IF (IO%IU0>=0.AND.NELM==1) WRITE(IO%IU0,'(" Calculation of long-wave limit done ")')
      ELSE
         CHI%HEAD  = 0
         CHI%WING  = 0
         CHI%CWING = 0
      ENDIF
  
      !dump HEAD
      !at this stage CHI still points to polarizability at Gamma
      CALL DUMP_HEAD_GG( CHI, IMAG_GRIDS, &
         'HEAD OF MICROSCOPIC DIELECTRIC TENSOR (INDEPENDENT PARTICLE)',&
         IO%IU6, IO%NWRITE, .FALSE., -EDEPS/LATT_CUR%OMEGA, 1._q)

      PROFILING_STOP( 'add_head_and_wings_to_chi_qpt' )
      CONTAINS

      !*******************************************************************
      !! helper: 
      !> calculates CHI in frequency domain inside a tau loop 
      !
      !*******************************************************************
      SUBROUTINE CHI_QPT_IN_FREQ_DOMAIN()
         !local 
         INTEGER           :: NQ_COUNTER, NQ
         INTEGER           :: NK1, NK2
         INTEGER           :: IU6_TMP
         PROFILING_START( 'chi_qpt_in_freq_domain' )

         IU6_TMP=IU6_MEM
         ! 
         ! use the supercell version to determine the response CHI_TAU(NQ_COUNTER)
         ! 
         IF ( LUSESUPER ) THEN

            ! clean CHI_TAU
            DO NQ_COUNTER=1,S2E%NUMBER_OF_NQ  ! loop over  q-points in the IRZ
               NQ=S2E%NQ(NQ_COUNTER)
               ! only cleared if chi(tau) is really computed locally
               IF ( IMAG_GRIDS%T%LDO_POINT_LOCAL ) CALL CLEAR_RESPONSE(CHI_TAU(NQ_COUNTER))
            ENDDO 
            ! determine chi(tau) if node carries data
            IF ( IMAG_GRIDS%T%LDO_POINT_LOCAL ) THEN
               ! use the supercell version to determine the response CHI_TAU(NQ_COUNTER)
               CALL CALCULATE_RESPONSE_SUPER( SUPER, WHF, WGW, GO, GU, GDES_TAU, &
                    CHI_TAU, S2E%NUMBER_OF_NQ, S2E%NQ, LATT_CUR, IU6_MEM)
               IU6_MEM=-1
            ENDIF
            CALL STOP_TIMING("GREENS",IO%IU6,"*chi_super")

            IU6_TMP=IU6_MEM
            DO NQ_COUNTER=1,S2E%NUMBER_OF_NQ
               NQ=S2E%NQ(NQ_COUNTER)
               CHI=>CHI_QPT(NQ_COUNTER)
               ! at this point we have CHI_TAU at some tau and q we need to distribute it
               ! cos transform of chi from imaginary time to imaginary frequency
               CALL TRANS_TIME_FREQUENCY(GDES_TAU, GDES, CHI_TAU(NQ_COUNTER), CHI, &
                IMAG_GRIDS, NTAU_ROOT, .TRUE. , IO%IU0 )
               CALL STOP_TIMING("GREENS",IU6_TMP,"COS_TRANS")
               ! 
               ! store chi_GG'(q=0,w=0) to CHI0 for OEP at Gamma point
               ! 
               IF( ( PRESENT(CHI0).AND.LOEP ) .AND. SQRT(DOT_PRODUCT(WDES%VKPT(:,NQ),WDES%VKPT(:,NQ))) < 1E-8_q ) THEN
                  CALL INTEGRATE_CHI_TAU(GDES_TAU, GDES, CHI_TAU(NQ_COUNTER), CHI0(1), &
                     IMAG_GRIDS, NTAU_ROOT, IO%IU0 )
               ENDIF
               IU6_TMP=-1
            ENDDO 
         ELSE
         ! 
         ! use the primitive cell version to determine the response CHI_TAU(NQ_COUNTER)
         ! 
            IU6_TMP=IU6_MEM
            DO NQ_COUNTER=1,S2E%NUMBER_OF_NQ  ! loop over  q-points in the IRZ
               NQ=S2E%NQ(NQ_COUNTER)
               CHI=>CHI_QPT(NQ_COUNTER)
               IF(IO%IU0>=0) WRITE(IO%IU0,'("NQ=",I4,3F10.4," TAU=",1F10.4", ")')&
                  NQ_COUNTER,CHI%VKPT(1:3),IMAG_GRIDS%T%POINT_CURRENT
               IF(IO%IU0>=0) WRITE(17,'("NQ=",I4,3F10.4," TAU=",1F10.4", ")')&
                  NQ_COUNTER,CHI%VKPT(1:3),IMAG_GRIDS%T%POINT_CURRENT

               ! this is where super and primitive cell approach differ
               IF ( IMAG_GRIDS%T%LDO_POINT_LOCAL ) CALL CLEAR_RESPONSE(CHI_TAU(1)) !one q-point at a time
               CALL SET_RESPONSE_KPOINT(CHI_TAU(1), WDES%VKPT(:,NQ), NQ)
    
               ! use the primitive cell version to determine the response CHI_TAU(NQ_COUNTER)
               IF ( IMAG_GRIDS%T%LDO_POINT_LOCAL ) THEN 
                  DO NK1=1,WDES%NKPTS
                    IF (MOD(NK1-1,WDES%COMM_KINTER%NCPU)/=WDES%COMM_KINTER%NODE_ME-1) CYCLE
   
                    NK2=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,NK1)-CHI_TAU(1)%VKPT(:),KPOINTS_FULL)
                    !contract GO and GU into CHI_TAU
                    CALL CALCULATE_CHI_TAU_FROM_G( WHF, WGW, GO(NK1), GU(NK2), GDES_TAU, &
                         CHI_TAU(1), NQ, NK1, NK2, LATT_CUR, IU6_MEM)
                    IU6_MEM=-1
                  ENDDO
               ENDIF
               CALL STOP_TIMING("GREENS",IU6_TMP,"*chi_tau")

               ! at this point we have CHI_TAU at some tau and q we need to distribute it
               ! cos transform of chi from imaginary time to imaginary frequency
               CALL TRANS_TIME_FREQUENCY(GDES_TAU, GDES, CHI_TAU(1), CHI, &
                IMAG_GRIDS, NTAU_ROOT, .TRUE. , IO%IU0 )
               ! 
               ! store chi_GG'(q=0,w=0) to CHI0 for OEP at Gamma point
               ! 
               IF( ( PRESENT(CHI0).AND.LOEP ) .AND. SQRT(DOT_PRODUCT(WDES%VKPT(:,NQ),WDES%VKPT(:,NQ))) < 1E-8_q ) THEN
                  CALL INTEGRATE_CHI_TAU(GDES_TAU, GDES, CHI_TAU(1), CHI0(1), &
                     IMAG_GRIDS, NTAU_ROOT, IO%IU0 )
               ENDIF
               CALL STOP_TIMING("G",IU6_TMP,"*cos_trans")
               IU6_TMP=-1
            ENDDO
         ENDIF

         CALL STOP_TIMING("GREENS",IO%IU6,"CHI")
 
         ! deallocate Green function arrays (frees some temporary storage)
         ! gK: TODO should be done inside CALCULATE_RESPONSE_SUPER to give us some memory
         CALL DEALLOCATE_G( GU)
         CALL DEALLOCATE_G( GO)
         DEALLOCATE(GU, GO)

         PROFILING_STOP( 'chi_qpt_in_freq_domain' )
      END SUBROUTINE CHI_QPT_IN_FREQ_DOMAIN


      !*******************************************************************
      ! synchronize 
      !*******************************************************************
        SUBROUTINE SYNC_CHI(WDES, CHI)
          USE ini
          TYPE (wavedes):: WDES
          TYPE (responsefunction) :: CHI(:)  ! response function
        ! local
          INTEGER :: I
      
          IF (WDES%COMM_KINTER%NCPU==1) RETURN
      
          DO I=1,SIZE(CHI)
             IF (CHI(I)%LREALSTORE) THEN
                CALLMPI( M_sum_d(WDES%COMM_KINTER, CHI(I)%RESPONSER, SIZE(CHI(I)%RESPONSER)) )
             ELSE
                CALLMPI( M_sum_z(WDES%COMM_KINTER, CHI(I)%RESPONSEFUN, SIZE(CHI(I)%RESPONSEFUN)) )
             ENDIF
          ENDDO
        END SUBROUTINE SYNC_CHI
   END SUBROUTINE DETERMINE_CHI_Q_OMEGA

!*************************************************************************
!
! subroutine to recalculate W in the self-consistent case
!
! This routine also determines the Phi functional of the grand 
! canonical potential. That is  the
! correlation part of the Luttinger Ward term in the GW approximation
!
! Phi^GW_c = 1/2 * Tr [ Log( 1 - GGV ) + GGV ]
!
! note, G is the interacting Green's function
!
!*************************************************************************

   SUBROUTINE RECALCULATE_W_AND_PHI_FUNCTIONAL( CHI_QPT, CHI_TAU, &
      GO_MAT, GU_MAT, GDES_MAT, &
      SUPER, S2E, W, WDES, WHF, WCORR, WCORRHF, WGW, WGWQ, &
      GDES_TAU, GDES, KPOINTS, NKPTS_IRZ, LATT_CUR, &
      RES_KPTS_TRANS, E, COR, FSG0, &
      NELM, IMAG_GRIDS, IO, IU6_MEM )

      USE acfdt_GG

      TYPE (responsefunction), POINTER :: CHI_QPT(:)
      TYPE (responsefunction), POINTER :: CHI_TAU(:)
      GDEF, POINTER, CONTIGUOUS        :: GO_MAT(:,:,:,:)
      GDEF, POINTER, CONTIGUOUS        :: GU_MAT(:,:,:,:)
      TYPE (greens_mat_des), POINTER   :: GDES_MAT  ! scGF in time
      TYPE (supercell), POINTER        :: SUPER
      TYPE (screened_2e_handle)        :: S2E
      TYPE (wavespin)                  :: W
      TYPE (wavedes)                   :: WDES
      TYPE (wavespin)                  :: WHF
      TYPE (wavespin)                  :: WCORR
      TYPE (wavespin)                  :: WCORRHF 
      TYPE (wavedes), POINTER          :: WGW
      TYPE (wavedes1)                  :: WGWQ
      TYPE (greensfdes)                :: GDES_TAU  ! GF in time
      TYPE (greensfdes)                :: GDES      ! GF in freq. 
      TYPE (kpoints_struct)            :: KPOINTS
      INTEGER                          :: NKPTS_IRZ ! tot # of k-points in the IRZ
      TYPE (latt)                      :: LATT_CUR
      TYPE (skpoints_trans)            :: RES_KPTS_TRANS  ! handle to generate response function at k-points outside the IRZ
      TYPE (energy)                    :: E 
      TYPE(correlation), POINTER       :: COR
      REAL(q)                          :: FSG0
      INTEGER                          :: NELM
      TYPE (imag_grid_handle)          :: IMAG_GRIDS
      TYPE (in_struct)                 :: IO
      INTEGER                          :: IU6_MEM
      ! local
      TYPE (responsefunction), POINTER :: CHI=>NULL() ! temporary response function in nu domain
      INTEGER                          :: ISP, I
      INTEGER                          :: NQ, NQ_COUNTER
      INTEGER                          :: NQ_IRZ

      IF ( NELM < 2 ) RETURN
      !               ALGO = GWR   .OR.           ALGO=GW0R at T>0
      IF ( (LGW.AND.((.NOT. LG0W0 ).AND.(.NOT. LscQPGW)) ) .OR. &
         (LFINITE_TEMPERATURE.AND.LGW0) ) THEN

         IF (IO%IU0>=0 .AND. IO%NWRITE>3) WRITE (IO%IU0,*) 'Iteration ', NELM, ' for W started'
         IF (IO%IU6>=0 .AND. IO%NWRITE>3) WRITE (IO%IU6,*) 'Iteration ', NELM, ' for W started'
  
         DO NQ_COUNTER=1,S2E%NUMBER_OF_NQ  !  loop over q-points in the IRZ
            NQ=S2E%NQ(NQ_COUNTER)
            CALL SET_RESPONSE_KPOINT(CHI_QPT(NQ_COUNTER), WDES%VKPT(:,NQ), NQ)
            CHI=>CHI_QPT(NQ_COUNTER)
            CALL CLEAR_RESPONSE(CHI)
            NULLIFY(CHI)
         ENDDO

         !====================================================================
         ! Determine CHI_QPT on imaginary frequency axis 
         ! no q->0 limit is added though, for scf GW runs this limit is
         ! interpolated
         ! leave out CRPA routines here 
         !====================================================================
         CALL DETERMINE_CHI_Q_OMEGA( CHI_QPT, CHI_TAU, &
            GO_MAT, GU_MAT, GDES_MAT, &
            SUPER, S2E, W, WDES, WHF, WCORR, WCORRHF, WGW, WGWQ, &
            GDES_TAU, GDES, KPOINTS, NKPTS_IRZ, LATT_CUR, &
            NELM, IMAG_GRIDS, IO, IU6_MEM, &
            .FALSE., .FALSE.,NELM.EQ.1 .OR. LscQPGW  ) 
   
         ! LFINITE_TEMPERATURE determines also Luttinger-Ward functional
         ! determine contribution of phi functional in RPA 
         IF ( LFINITE_TEMPERATURE ) THEN
            CALL SET_PHI_RPA_FUNCTIONAL()
         ENDIF
  
         ! screened Coulomb potential W
         CALL SCREENED_POTENTIAL(S2E, CHI_QPT, WGW, WGWQ, &
          WDES, GDES, LATT_CUR, FSG0, IMAG_GRIDS, .FALSE., IO ) 

      ENDIF

      CONTAINS
      ! computes correlation part of Luttinger Ward term in GW approximation
      ! Phi_GW^c = 1/2 * Tr [ Log( 1 - GGV ) + GGV ]
      SUBROUTINE SET_PHI_RPA_FUNCTIONAL()
         E%ERPA_CUT = 0
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
      
         ! reset CHI to G_0.G_0 
         IF ( LGW0 ) THEN
            DO NQ_COUNTER=1,S2E%NUMBER_OF_NQ  !  loop over q-points in the IRZ
               NQ=S2E%NQ(NQ_COUNTER)
               CALL SET_RESPONSE_KPOINT(CHI_QPT(NQ_COUNTER), WDES%VKPT(:,NQ), NQ)
               CHI=>CHI_QPT(NQ_COUNTER)
               CALL CLEAR_RESPONSE(CHI)
               NULLIFY(CHI)
            ENDDO
            CALL DETERMINE_CHI_Q_OMEGA( CHI_QPT, CHI_TAU, &
               GO_MAT, GU_MAT, GDES_MAT, &
               SUPER, S2E, W, WDES, WHF, WCORR, WCORRHF, WGW, WGWQ, &
               GDES_TAU, GDES, KPOINTS, NKPTS_IRZ, LATT_CUR, &
               NELM, IMAG_GRIDS, IO, IU6_MEM, &
               .TRUE., .FALSE., .TRUE.  )
         ENDIF
      END SUBROUTINE SET_PHI_RPA_FUNCTIONAL
  
   END SUBROUTINE RECALCULATE_W_AND_PHI_FUNCTIONAL


#endif

END MODULE greens_real_space
