#include "symbol.inc"

! use time propagation code of Friedhelm Bechstedts group (added by C. Roedl and jF)
!#define use_bse_te
!furth Currently, bse_te.F code cannot be compiled with -DwNGZhalf because this
!      part is fully complex, so we need to switch off use_bse_te in this case
#ifdef gammareal
#undef use_bse_te
#endif
#ifndef MPI
#undef use_bse_te
#endif
   
MODULE bse
  USE chi_glb
  USE local_field
  USE scala
  USE wave_high
#ifdef use_bse_te
  USE bse_te
#endif
  IMPLICIT NONE

!**********************************************************************
!
! this subroutine solves the BSE equation 
! presently in the Tamm Dancoff approximation
!
! it is based on the local_field subroutine and relies on the 
! routines implemented therein
! the required integrals correspond to the publication 
! Sander, Maggio, Kresse, Phys. Rev. B 92, 045209
! 
!  k1,n1=  i
!  k2,n2=  j
!  k3,n3=  a
!  k4,n4=  b
!
! The  direct contribution (A = <e,h |W| e',h'>)
! (two electrons at the same spatial coordinate, and two holes)
!
!     int_d3r d3r' c*_k1+q,n3(r') c'_k2+q,n4(r') v'*_k2,n2(r)   v_k1,n1(r)  W(r',r)
!
! and a second  exchange like contribution (B= <e,h | h',e'>)
! (lacking in the Tamm Dancoff approximation)
!    int_d3r d3r' c_k1-q,n3(r')  v'*_k2,n2(r') c'_k2+q,n4(r)  v*_k1,n1(r)  W(r',r)
!
! and a term stemming from the variation of the Hartree potential
!
!    int_d3r d3r' v'*_k2,n2(r') c'_k2+q,n4(r') c*_k1+q,n3(r) v_k1,n1(r) v(r',r)
!
!
! Note that the notation A and B refers to the quantum chemists notation,
! where the full TDFT-TDHF response function is often written as
!
!        K_ai,bj    K_ai,jb         A_ai,bj   B_ai,bj
!                                =  
!        K_ia,bj    K_ia,jb         B*_ai,bj  A*_ai,bj
!
! v states are restricted to valence band states 
! whereas c are conduction band states
! the sums are restricted to NBANDSLF bands below and above the Fermi-level

! loops over n1/n2/n3 and n4 are done in blocks of NSTRIP/NLOC
!  where NSTRIP is the number of local bands, and NLOC the total number
!
! outline of the structure of the routine:
!
!   allocate array W1 for merged block n1
!   allocate array W2 for merged block n2
!
!   allocate array W3 for merged block n3
!   allocate array W4 for merged block n4
!   allocate charge array (WKAPPA) of size merged block n1 x NBANDSLF
!     (this costs possibly a lot of memory and restricts n1)
!
!   allocate charge array of size merged WA block n4 x merged block n1
!
!   loop over all k1-points
!     loop overal all v1-states in blocks (preferably all n1)
!        gather wavefunction in the current block
!        loop overall k2-points
!           loop overal all v2-states in blocks (preferably all n2)
!              gather wavefunction in the current block
!
!              calculate F(n1,n3,n2,n4) = c*_k1+q,n3  c'_k2+q,n4(r') W(r',r) v'*_k2,n2(r), v_k1,n1(r)
!
! 
! ANTIRES = 0 corresponds to the Tamm Dancoff approximation, only direct terms
! ANTIRES = 1 applies a fancy approximation that includes
!             resonant-antiresonant coupling with w=0 (zero frequency) exact
! ANTIRES = 2 full version
!
!
!              
!**********************************************************************

! the routine has plenty of flags to determine the precise
! behaviour
!
! the flag single_prec_bse allows to select single precision storage of BSE 
! matrix and thus saves a factor 2 in storage (the default is double_prec_bse)
!
#if ! defined(double_prec_bse) && ! defined(single_prec_bse)
#define double_prec_bse
#endif

!
! the flag determines whether ZHEEVX is used
! ZHEEVX is about twice to three times faster than ZHEEV, but requires twice as much memory 
!
#define USE_ZHEEVX
!#define USE_ZHEEVD
#ifdef USE_ZHEEVD
#undef USE_ZHEEVX
#endif

! this is already defined in scala.F, but private there
! it is possible to undef scaLAPACK here, then MPI is still used but scaLAPACK is avoided
#ifdef scaLAPACK
      INTEGER,PARAMETER,PRIVATE :: BLOCK_CYCLIC_2D=1, DTYPE_=1, &
                           CSRC_ =8, CTXT_=2, DLEN_=9, LLD_=9, &
                           MB_=5, M_=3, NB_=6, N_=4, RSRC_=7
#else
! we use some descriptors DESC so set at least DLEN_
      INTEGER,PARAMETER :: DLEN_=9
#endif

  TYPE banddesc
    INTEGER :: VBMIN, VBMAX      ! band index for the VB minimum and VB maximum
    INTEGER :: CBMIN, CBMAX      ! band index for the CB minimum and CB maximum
    INTEGER :: CBMIN4, CBMAX4    ! band index for the CB minimum and CB maximum for n4
                                 ! this might be reduced for parallelization over bands
    INTEGER :: NGLB              ! block size for conduction bands = total number of conduction bands
    INTEGER :: NGLB4             ! block size for conduction bands for fourth index (N4 usually)
                                 ! [ N4 is distributed over nodes when band parallelization is used
                                 !   in this case NGLB4 is equal NGLB/ number of nodes]
  END TYPE banddesc

  TYPE bse_matrix_index
    INTEGER :: NCV
    INTEGER, ALLOCATABLE ::  N1(:)  ! valence band index
    INTEGER, ALLOCATABLE ::  N3(:)  ! conduction band index
    INTEGER, ALLOCATABLE ::  NK(:)  ! k-point index
    INTEGER, ALLOCATABLE ::  ISP(:) ! spin index
    INTEGER, ALLOCATABLE ::  INDEX(:,:,:,:) ! yield index into Hamilton matrix for given valence, conduction, k-point and spin
  END TYPE bse_matrix_index


  TYPE selfenergy_from_bse 
     INTEGER :: NOMEGA
     INTEGER :: NBANDS
     INTEGER :: NKPTS
     INTEGER :: ISPIN
     REAL(q), POINTER :: OMEGA(:)  => NULL()
     REAL(q), POINTER :: REAL_PART(:,:,:,:) => NULL() !1st index over frequencies, 2nd index over bands
     REAL(q), POINTER :: IMAG_PART(:,:,:,:) => NULL()!3rd index over kpoints, 4th index over spin 
     COMPLEX(q), POINTER :: CELNEW(:,:,:) => NULL()  !bands, kpoints, spin 
  END TYPE selfenergy_from_bse


CONTAINS

  SUBROUTINE CALCULATE_BSE( &
       IBSE,P,WDES,NONLR_S,NONL_S,W,LATT_CUR,LATT_INI,  &
       T_INFO,DYN,INFO,IO,KPOINTS,SYMM,GRID,LMDIM,CQIJ, &
       WGW, SHIFT,CELTOT_HK, CELTOT_X, & 
       NBANDSO,NBANDSV, OMEGAMAX_IN, NQPOINT, KPOINT_BSE, ISP_IN, & 
       LHARTREE, LADDER, LTRIPLET, ANTIRES, LGWLF, LFXC, NEDOS_IN, NELMGW, NBANDSGW, L2ORDER, &
       LCORBSE, LQP, LTCTC, CORREP, LAMBDA, S2E, EFERMI, SEFROMBSE, FIRST_EMPTYSTATE, LAST_OCCSTATE)

    USE base
    USE pseudo
    USE nonl_high
    USE msymmetry
    USE mpimy
    USE mgrid
    USE mkpoints
    USE constant
    USE poscar
    USE pot
    USE pawm
    USE kpoints_change
    USE full_kpoints
    USE mlr_optic
    USE dfast
    USE choleski
    USE ini
    USE gauss_quad 
    USE mymath
    USE tutor, ONLY: vtutor, isAlert, isError, BSEantires, BSEncv
    IMPLICIT NONE
    ! structures
    INTEGER              IBSE       ! 0 matrix diagonalization, 1 time evolution, 2 iterative (not supported)
    TYPE (type_info)     T_INFO
    TYPE (potcar)        P(T_INFO%NTYP)
    TYPE (wavedes)       WDES
    TYPE (nonlr_struct)  NONLR_S
    TYPE (nonl_struct)   NONL_S
    TYPE (wavespin)      W
    TYPE (latt)          LATT_CUR, LATT_INI
    TYPE (dynamics)      DYN
    TYPE (info_struct)   INFO
    TYPE (in_struct)     IO
    TYPE (kpoints_struct) KPOINTS
    TYPE (symmetry)      SYMM
    TYPE (grid_3d)       GRID       ! grid for wavefunctions
    INTEGER  LMDIM
    OVERLAP  CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    TYPE (wavedes)       WGW        ! descriptor for basis set of response function
    REAL(q)              SHIFT      ! complex frequency shift for broadening
    COMPLEX(q) :: CELTOT_HK(:,:,:)  ! this contains the Hartree + kinetic term
    COMPLEX(q) :: CELTOT_X(:,:,:)   ! this contains the exchange term
    INTEGER              NBANDSO    ! number of occupied states included in BSE
    INTEGER              NBANDSV    ! number of virtual states included in BSE
    REAL(q)              OMEGAMAX_IN! maximum considered transition energy in BSE
    INTEGER              NQPOINT    ! q-point for which excitations are calculated (usually optical e.g. q=0)
    INTEGER              KPOINT_BSE(3) ! reciprocal lattice vector G (g=G+q) at which response is evaluated
    INTEGER              ISP_IN     ! in: spin index to be considered:
                                    !  ISP_IN<0: both spin channels are treated simultaneously (only valid for ISPIN=2)
                                    !  ISP_IN=1 : treat only spin channel 1 (and pretend the system is non-magnetic)
                                    !  ISP_IN=2 : treat only spin channel 2 (and pretend the system is non-magnetic)
    LOGICAL              LHARTREE   ! include Hartree effects
    LOGICAL              LTRIPLET   ! calculate triplet
    LOGICAL              LADDER     ! local field effects beyond RPA (ladder diagrams)
    INTEGER              ANTIRES    ! how to treat antiresonant part
    ! presently implemented 0 = Tamm-Dankoff
    ! presently implemented 1 -> correct w=0 properties 
    LOGICAL              LGWLF      ! use screened W from GW for exchange interactions
    LOGICAL              LFXC       ! include DFT xc kernel
    INTEGER              NEDOS_IN   ! number grid points in DOS (read from INCAR)
    INTEGER              NELMGW     ! number of timesteps in solving the BSE equation using time propagation
    INTEGER              NELM       ! NELMGW is also the number of iterations in
                                    ! the GW0Gamma0 QP energies evaluations
                                    ! reqquired here for the gammareal version
    INTEGER              NBANDSGW   ! number of bands for QP calculation
    LOGICAL              L2ORDER    ! only second order contribution
    LOGICAL              LCORBSE    ! determine correlation energy from BSE
    LOGICAL              LQP        ! determine self-energy
    LOGICAL              LTCTC      ! if true computes the test charge-test charge selfenergy
    REAL(q), OPTIONAL :: LAMBDA     ! coupling constant value can be used to scale the Coulumb integrals
    TYPE (screened_2e_handle), OPTIONAL :: S2E !screened two electron integrals
    REAL(q), OPTIONAL :: EFERMI
    TYPE(selfenergy_from_bse), OPTIONAL :: SEFROMBSE
    INTEGER, OPTIONAL :: FIRST_EMPTYSTATE(:,:,:), LAST_OCCSTATE(:,:,:)
! local
    REAL(q)              OMEGAMAX   ! maximum considered transition energy in BSE
    INTEGER              NEDOS      ! number grid points in DOS
    TYPE(banddesc) :: BD(WDES%ISPIN)! descriptor for describing which bands are included, spin up
    TYPE(banddesc) :: B1,B2         ! current descriptor for the spin loops
    TYPE(bse_matrix_index) :: BSE_INDEX
    TYPE(bse_matrix_index) :: POT_INDEX
    ! the following matrices store are sub range of the full matrices AMAT, BMAT and VMAT
    GDEF, ALLOCATABLE :: TWOELECTRON3O(:,:,:,:)   ! part of the A matrix resonant-resonant coupling
    GDEF, ALLOCATABLE :: BTWOELECTRON3O(:,:,:,:)  ! part of the B matrix resonant-antiresonant coupling 
    GDEF, ALLOCATABLE :: TWOE3O_HARTREE(:,:,:,:)  ! part of Coulomb kernel in four orbital representation (Hartree only no exchange)
    TYPE (wavespin) WHF
    TYPE (wavefun1),ALLOCATABLE :: W1(:), W2(:), W3(:), W4(:)
    TYPE (wavedes1), TARGET :: WDESK1, WDESK2, WDESK3, WDESK4
    INTEGER :: NSTRIP               ! block size used internally in TWOELECTRON4O routines
    INTEGER :: NSTRIPV              ! block size for valence bands (valence bands are not distributed over nodes)
    INTEGER :: N 
    INTEGER :: NPOS1, NSTRIP1       ! base index and width of the n1 block
    INTEGER :: NPOS2, NSTRIP2       ! base index and width of the n1 block
    INTEGER :: K1, K2, K3, K4, K2_LOCAL, K2_COLLECT, K2_DONE, K4_LOCAL
    INTEGER :: ALLOC_VALENCE_BANDS, ALLOC_CONDUCTION_BANDS, ALLOC_KPOINTS
    INTEGER :: NCBD, NCBD4          ! max # of conduction bands over different spin orientations
    REAL(q) :: NFFTW
    REAL(q) :: NFLOAT4O, NFFT4O     ! number of BLAS3 operations in 4 orbital routines
    LOGICAL :: LEX_INTERPOLATED     ! use interpolation of the exchange like contribution (B)
    LOGICAL :: LscaLAPACKaware      ! distribute data using BLACS
    INTEGER :: NCV                  ! number of conduction valence band pairs 
    INTEGER :: NAB                  ! NBANDSGW x (conduction bands + valence bands)
    INTEGER :: IFAIL               
    INTEGER :: ISP, ISP2, ISP_LOW, ISP2_LOW, ISP_HIGH
    LOGICAL :: W1EQUALW2, LKPOINT_PARALLEL
    REAL(q), PARAMETER :: G2ZERO=1E-12_q
    REAL(q) :: AHARTREE
    LOGICAL :: LFULL
    LOGICAL :: LTRIANGLE            ! only calculate lower triangle of A and B matrix and add remainder using conjugation
    REAL(q) :: POTFAK(GRIDHF%MPLWV)       ! 1/(G+dk)**2 (G)
#ifdef single_prec_bse
    GDEFS, ALLOCATABLE :: AMAT(:,:), AMAT_SCALA(:)  ! A matrix used in BSE
    GDEFS, ALLOCATABLE :: BMAT(:,:), BMAT_SCALA(:)  ! B matrix used in BSE
    GDEFS, ALLOCATABLE :: VMAT(:,:), VMAT_SCALA(:)  ! Coulomb kernel in four orbital representation (Hartree only no exchange)
    GDEFS, ALLOCATABLE :: AVpW(:,:), AVpW_SCALA(:)  ! A matrix, including elements between occup-occup and unoccup-unoccup
    GDEFS, ALLOCATABLE :: BVpW(:,:), BVpW_SCALA(:)  ! B matrix, including elements between occup-occup and unoccup-unoccup
#endif

#ifdef double_prec_bse
    GDEF, ALLOCATABLE :: AMAT(:,:), AMAT_SCALA(:)
    GDEF, ALLOCATABLE :: BMAT(:,:), BMAT_SCALA(:)
    GDEF, ALLOCATABLE :: VMAT(:,:), VMAT_SCALA(:)
    GDEF, ALLOCATABLE :: AVpW(:,:), AVpW_SCALA(:)
    GDEF, ALLOCATABLE :: BVpW(:,:), BVpW_SCALA(:)
#endif
    COMPLEX(q) :: CSUM_A    ! sum of diagonal components of matrix A
    COMPLEX(q) :: SUM_ALPHA ! accumulation variable for the correlation energy
                            ! along the AC path
    COMPLEX(q) :: TR_ALPHA  ! =tr{V*P_alpha} with P the polarisability
    REAL(q), OPTIONAL :: CORREP       ! correlation energy from plasmon formula @ given k-point 
    LOGICAL :: LALPHA       ! true if the coupling constant value LAMBDA has being passed by calling routine
    REAL(q) :: ALPHAV(10), W_GL(10) ! variables for the AC integration, max 10 
    INTEGER :: NALPHA       ! number of points considered in coupling constant integration
    INTEGER :: IALPHA       ! current point considered in coupling constant integration
    REAL(q), PARAMETER :: SCALE_PERT=1._q ! this parameter can be used to scale the Coulomb interaction

    INTEGER :: NSIZE        ! counter for allocation
    
    GDEF                 :: TMP
    REAL(q), ALLOCATABLE :: R(:), RB(:)
! scaLAPACK
    INTEGER :: NP,NQ,NPROW,NPCOL,MYROW,MYCOL
    INTEGER, EXTERNAL    ::     NUMROC
    INTEGER :: nd
#ifdef scaLAPACK
    TYPE(scalapack_des) :: BSE_GSD, POT_GSD
#endif
    INTEGER :: BSE_DESC(DLEN_), POT_DESC(DLEN_)
    INTEGER :: ISTAT  ! allocation error code 


#ifdef MPI
    IF (WDES%COMM_KINTER%NCPU /= 1 .AND. .NOT. LQP) THEN
       CALL vtutor%error("CALCULATE_BSE: KPAR>1 not implemented, sorry.")
    END IF
#endif
    MULTIPLY_FRACTIONAL_OCC=.TRUE.
    LTRIANGLE=.TRUE.

    IF (LQP) THEN
       ! usually local_field.F restricts calculations to particle-hole pairs (BSE)
       ! unsetting FAST_CYCLE and MULTIPLY_FRACTIONAL_OCC allows to calculate all
       ! which is required to obtain the entire Coulomb kernel (including e.g. particle-particle)
       ! this causes issues for partially filled states in the BSE,
       ! hence avoid to perform quasi-particle calculations for metals
       ! TODO for LQP:
       !   - try to get rid of LHEAD_DIRECT tag
       !   - full scaLAPACK version of SELFEN
       !   - around 5047 see comment gKtest
       !   - it would be better to centre the self-energy at the mean-field energies 
       !     (or at least Fermi-energy) to allow a slight reduction of default OMEGAMAX
       !   - C3 with 4x4x4 k-points crashed but maybe NOMEGA was missing
       FAST_CYCLE=.FALSE.
       MULTIPLY_FRACTIONAL_OCC=.FALSE. 
       LCORBSE=.FALSE.
       LTRIANGLE=.FALSE.
       LHEAD_DIRECT=.TRUE.
    ENDIF
!=======================================================================
! somewhat akward: the routine can solve the BSE equation
! at various Coulomb couplings in the list ALPHAV(1:NALPHA)
! this allows to perform a coupling constant integration to obtain
! the correlation energy
!=======================================================================
    LALPHA=.FALSE.              ! true if the coupling constant value LAMBDA has being passed by calling routine
    IF ( LCORBSE ) THEN
       NALPHA=6                 ! use six point Gauss-Legendre integration as default
       IF (PRESENT(LAMBDA)) THEN! if LAMBDA is an integer value use that as NALPHA
          IF (LAMBDA==INT(LAMBDA) .AND. LAMBDA>=2.0) NALPHA=MIN(INT(LAMBDA),10)
       ENDIF
       IF (IO%IU0>0) WRITE(IO%IU0,*) 'Gauss-Legendre quadrature with ',NALPHA, 'points'
       IF (IO%IU6>0) WRITE(IO%IU6,*) 'Gauss-Legendre quadrature with ',NALPHA, 'points'
       
       CALL GAUSS_LEGENDRE(0.0_q, 1.0_q, ALPHAV, W_GL, NALPHA)
    ELSE
       NALPHA=1
       W_GL=1.0; ALPHAV=1.0     ! weight 1 and unscaled Coulomb kernel
    ENDIF
    IF ( PRESENT(LAMBDA) ) THEN
       IF (LAMBDA /=1.0_q .AND. LAMBDA/=INT(LAMBDA)) THEN
          LALPHA=.TRUE.
          NALPHA=1
          W_GL=1.0; ALPHAV=LAMBDA ! weight 1, but scaled Coulomb kernel
       ENDIF
    ENDIF

    ! smaller values are dangerous and are usually not what the user desires
    NEDOS=MAX(NEDOS_IN,1000)

    IF (ISP_IN<0 .AND. WDES%ISPIN==1) THEN
       CALL vtutor%bug("internal error in CALCULATE_BSE: ISP is negative but ISPIN=1", __FILE__, __LINE__)
    ENDIF

    IF (WDES%LGAMMA .AND. (.NOT. WGW%GRID%REAL2CPLX .OR. .NOT. WGW%GRID%LREAL)) THEN
       CALL vtutor%bug("internal error in CALCULATE_BSE: WGW must apply real to complex FFT", __FILE__, __LINE__)
    ENDIF
!=======================================================================
! preparation
!=======================================================================
    ! read in WPOT files
    ! safer here, since it can happen that not all nodes call GET_WPOT simultaneously
    ! (e.g. k-point parallelization)
    ! this might cause problems with SHMEM calls
    IF (LADDER .AND. LGWLF .AND. .NOT. ASSOCIATED(WPOTH)) THEN
       CALL INIT_WPOT_HANDLE( WPOTH, WGW, KPOINTS_FULL%NKPTS, IO%IU6, IO%IU0, 1, 1, 1 )
       DO K1=1,KPOINTS_FULL%NKPTS
          IF (W%WDES%WTKPT(K1)==0) CYCLE
! NKREDLF not presently used by bse.F
!          IF (SKIP_THIS_KPOINT_IN_LF(W%WDES%VKPT(:,K1), NKREDLFX, NKREDLFY, NKREDLFZ)) CYCLE
          K2=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,K1)-W%WDES%VKPT(:,1),KPOINTS_FULL)
          IF (ASSOCIATED(WPOTH)) CALL GET_WPOT( WPOTH, K2, POTFAK, LFULL)
       ENDDO
    ENDIF

    WHF=W      ! use temporarily another WDES
    WHF%WDES => WDES_FOCK

!  switch off NKRED
    NKREDX=1
    NKREDY=1
    NKREDZ=1

! determine wether the k1+k2 is found in the full k-point grid or not
    DO K1=1,KPOINTS_FULL%NKPTS
       IF (LIDENTICAL_KPOINT(WHF%WDES%VKPT(:,1)+WHF%WDES%VKPT(:,MIN(KPOINTS_FULL%NKPTS,2)),KPOINTS_FULL%VKPT(:,K1))) EXIT
    ENDDO
   
    IF (K1==KPOINTS_FULL%NKPTS+1 .AND. ANTIRES>=1) THEN
       LEX_INTERPOLATED=.TRUE.      ! not found use LEX_INTERPOLATED (slower)
    ELSE
       LEX_INTERPOLATED=.FALSE.     ! found, LEX_INTERPOLATED not required
    ENDIF

! set the wavefunction descriptors
    CALL SETWDES(WHF%WDES,WDESK1,0)
    CALL SETWDES(WHF%WDES,WDESK2,0)
    CALL SETWDES(WHF%WDES,WDESK3,0)
    CALL SETWDES(WHF%WDES,WDESK4,0)

    IF (IO%IU6>=0) THEN
      WRITE(IO%IU6,'("NQ=",I4,3F10.4,", ")') NQPOINT,WHF%WDES%VKPT(:,NQPOINT)
    ENDIF

    W1EQUALW2=.FALSE.

    DO ISP=1,WDES%ISPIN
       IF (LQP) THEN
             CALL SET_BAND_PARAMETERS(W, ISP, NBANDSO, NBANDSV,  & 
                  BD(ISP)%VBMAX, BD(ISP)%VBMIN, BD(ISP)%CBMIN, BD(ISP)%CBMAX, & 
                  NSTRIPV, NSTRIP, BD(ISP)%NGLB, BD(ISP)%CBMIN4, BD(ISP)%CBMAX4, BD(ISP)%NGLB4, W1EQUALW2, LKPOINT_PARALLEL, &
                  1, NBANDSGW)
          ELSE
             CALL SET_BAND_PARAMETERS(W, ISP, NBANDSO, NBANDSV,  & 
                  BD(ISP)%VBMAX, BD(ISP)%VBMIN, BD(ISP)%CBMIN, BD(ISP)%CBMAX, & 
                  NSTRIPV, NSTRIP, BD(ISP)%NGLB, BD(ISP)%CBMIN4, BD(ISP)%CBMAX4, BD(ISP)%NGLB4, W1EQUALW2, LKPOINT_PARALLEL)
          ENDIF
    ENDDO

    ! if ISP_IN is set, only the channel ISP_IN is considered
    IF (ISP_IN>0) THEN
       IF (LQP) THEN
          CALL SET_BAND_PARAMETERS(W, ISP_IN, NBANDSO, NBANDSV, & 
               BD(ISP_IN)%VBMAX, BD(ISP_IN)%VBMIN, BD(ISP_IN)%CBMIN, BD(ISP_IN)%CBMAX, & 
               NSTRIPV, NSTRIP, BD(ISP_IN)%NGLB, & 
               BD(ISP_IN)%CBMIN4, BD(ISP_IN)%CBMAX4, BD(ISP_IN)%NGLB4, W1EQUALW2, LKPOINT_PARALLEL, &
               1, NBANDSGW)
       ELSE
          CALL SET_BAND_PARAMETERS(W, ISP_IN, NBANDSO, NBANDSV, & 
               BD(ISP_IN)%VBMAX, BD(ISP_IN)%VBMIN, BD(ISP_IN)%CBMIN, BD(ISP_IN)%CBMAX, & 
               NSTRIPV, NSTRIP, BD(ISP_IN)%NGLB, & 
               BD(ISP_IN)%CBMIN4, BD(ISP_IN)%CBMAX4, BD(ISP_IN)%NGLB4, W1EQUALW2, LKPOINT_PARALLEL)
       ENDIF
       ALLOC_VALENCE_BANDS=BD(ISP_IN)%VBMAX-BD(ISP_IN)%VBMIN+1
       ALLOC_CONDUCTION_BANDS=BD(ISP_IN)%NGLB
       ISP_LOW =ISP_IN
       ISP_HIGH=ISP_IN
    ELSE 
    ! spin polarized calculations
       ALLOC_VALENCE_BANDS=MAX(BD(1)%VBMAX-BD(1)%VBMIN+1,BD(2)%VBMAX-BD(2)%VBMIN+1)
       ALLOC_CONDUCTION_BANDS=MAX(BD(1)%NGLB,BD(2)%NGLB)
       ISP_LOW =1
       ISP_HIGH=2
       ! TODO this requires a little bit more consideration
       ! for spin polarized calculations at the gamma point this increases the memory demand
       W1EQUALW2=.FALSE.
    ENDIF

    IF (ABS(WHF%WDES%VKPT(1,NQPOINT))>1E-6_q .OR. ABS(WHF%WDES%VKPT(2,NQPOINT))>1E-6_q .OR. ABS(WHF%WDES%VKPT(3,NQPOINT))>1E-6_q) THEN
       ! allocate space for all k-points in the full wedge
       CALL  ALLOCATE_CDER_BETWEEN_STATES_Q( WHF%WDES%NB_TOT, MAX(BD(ISP_LOW)%VBMAX, BD(ISP_HIGH)%VBMAX), & 
            WHF%WDES%NKPTS , WHF%WDES%ISPIN)
    ENDIF

    IF (LKPOINT_PARALLEL) THEN
       IF (IO%IU6>=0) WRITE(IO%IU6,*) 'parallelization over k-points'
    ELSE
       IF (IO%IU6>=0) WRITE(IO%IU6,*) 'parallelization over bands'
    ENDIF

    ! allocate the orbitals: these contain the states in the conduction or valence bands 
    ALLOCATE(W1 (ALLOC_VALENCE_BANDS))
    ALLOCATE(W2 (ALLOC_VALENCE_BANDS))
    ALLOCATE(W3 (ALLOC_CONDUCTION_BANDS))
    ALLOCATE(W4 (ALLOC_CONDUCTION_BANDS))

    DO N=1,ALLOC_VALENCE_BANDS
       CALL NEWWAV(W1(N) , WDESK1,.TRUE.)
       IF (.NOT. W1EQUALW2) CALL NEWWAV(W2(N) , WDESK2,.TRUE.)
    ENDDO

    DO N=1,ALLOC_CONDUCTION_BANDS
       CALL NEWWAV(W3(N) , WDESK3,.TRUE.)
       IF (.NOT. W1EQUALW2) CALL NEWWAV(W4(N) , WDESK4,.TRUE.)
    ENDDO

    ! set variable nd equal to ndata or twice ndata
    IF (ANTIRES >=2) THEN
       nd=2*ndata
    ELSE
       nd=ndata
    ENDIF

    OMEGAMAX=OMEGAMAX_IN
    IF (IO%IU6>=0) THEN
       DO ISP=ISP_LOW,ISP_HIGH
          WRITE(IO%IU6,'(/A,I2,/2(A,I5,2X,A,I5/))')' Bands included in the BSE spin=',ISP, & 
               ' VB(min)=',BD(ISP)%VBMIN,'VB(max)=',BD(ISP)%VBMAX,' CB(min)=',BD(ISP)%CBMIN,'CB(max)=',BD(ISP)%CBMAX
       ENDDO

       IF (OMEGAMAX>0 .AND. .NOT. LQP ) THEN
          WRITE(IO%IU6,'(A,F10.5)') ' electron-hole pairs beyond OMEGAMAX are removed OMEGAMAX=',OMEGAMAX
       ENDIF

       IF (LADDER .AND. (AEXX /=0 .OR. L_MODEL_HF .OR. LGWLF) ) THEN
          IF (LGWLF) THEN
             WRITE(IO%IU6,'(A,3F14.7)')' W is read from the files WXXXX.tmp or WFULLXXXX.tmp if present'
          ELSE
             WRITE(IO%IU6,'(A,F14.7,A,F14.7)')' parameters for screened Coulomb W: AEXX=',AEXX,' HFSCREEN=',HFSCREEN
          ENDIF
          IF (.NOT. LHARTREE .OR. LTRIPLET) THEN
             WRITE(IO%IU6,'(A)') ' only ladders are included (Hartree term=RPA part is switched off)'
          ENDIF
       ELSE
          IF (.NOT. LHARTREE .OR. LTRIPLET) THEN
             WRITE(IO%IU6,'(A)') ' IP particle spectrum'
          ELSE
             WRITE(IO%IU6,'(A)')' simple RPA calculation, excitonic effects (ladders) are not included'
          ENDIF
       ENDIF
       IF (LHARTREE .AND. LFXC .AND. .NOT. LGWLF) THEN
          WRITE(IO%IU6,'(A)')' xc-kernel from LOCAL part of DFT is included (LFXC=.TRUE.)'
       ENDIF
    ENDIF

!==========================================================================
! few BSE specific things
!==========================================================================
    ALLOC_KPOINTS=0
    DO K1=1,WDES%NKPTS
       IF (WHF%WDES%WTKPT(K1)==0) CYCLE
       ALLOC_KPOINTS=ALLOC_KPOINTS+1
    ENDDO
    IF (LQP) THEN
       CALL SET_BSE_MATRIX_INDEX( ISP_LOW, ISP_HIGH, ALLOC_KPOINTS, BD, BSE_INDEX, WHF, -OMEGAMAX, NCV, NQPOINT )
       CALL SET_POT_MATRIX_INDEX( ISP_LOW, ISP_HIGH, ALLOC_KPOINTS, BD, POT_INDEX, WHF, NAB, NQPOINT )
    ELSE
       CALL SET_BSE_MATRIX_INDEX( ISP_LOW, ISP_HIGH, ALLOC_KPOINTS, BD, BSE_INDEX, WHF, OMEGAMAX, NCV, NQPOINT )
    ENDIF
    IF (NCV ==0) THEN
       CALL vtutor%write(isError, BSEncv)
    ENDIF

    LscaLAPACKaware=.FALSE.
#ifdef scaLAPACK
    IF (LscaLAPACK) THEN
       LscaLAPACKaware=.TRUE.
    ENDIF
    ! time evolution algorithm does not handle distributed matrices, so switch it off
    IF (IBSE>0) LscaLAPACKaware=.FALSE.
    ! QP calculation does not yet support scaLAPACK
    IF (LQP)    LscaLAPACKaware=.FALSE.
#endif
!==========================================================================
! intialize scaLAPACK with appropriate matrix dimensions
!==========================================================================
    IF ( IBSE==0) THEN
! if scaLAPACK is available initialize BSE_DESC
#ifdef scaLAPACK
    CALL INIT_scala_DESC( WDES%COMM_INTER, NCV, BSE_DESC, BSE_GSD )
    CALL BLACS_GRIDINFO(BSE_DESC(CTXT_), NPROW, NPCOL,MYROW,MYCOL)
    NP = NUMROC(NCV,BSE_DESC(MB_),MYROW,0,NPROW)
    NQ = NUMROC(NCV,BSE_DESC(NB_),MYCOL,0,NPCOL)
#endif
    IF ( LscaLAPACKaware) THEN
#ifdef scaLAPACK
       ALLOCATE( AMAT_SCALA( MAX(NP*NQ,1)), STAT=ISTAT)
       CALL CHECK_ALLOCATION_ERROR( "AMAT_SCALA", MAX(NP*NQ,1), ISTAT )
       ALLOCATE( AMAT( 1, 1), R(NCV))
       NSIZE=SIZE(AMAT_SCALA)

       IF (ANTIRES>=2) THEN
          ALLOCATE(BMAT_SCALA(MAX(NP*NQ,1)), STAT=ISTAT)
          CALL CHECK_ALLOCATION_ERROR( "BMAT_SCALA", MAX(NP*NQ,1), ISTAT )
          NSIZE=NSIZE+SIZE(BMAT_SCALA)
          ALLOCATE(BMAT(1,1), RB(NCV))
       ENDIF
       IF ((.NOT. LQP) .AND. (NALPHA > 1 .OR. LALPHA))  THEN
          ALLOCATE( VMAT_SCALA(MAX(NP*NQ,1)), STAT = ISTAT )
          CALL CHECK_ALLOCATION_ERROR( "VMAT_SCALA", MAX(NP*NQ,1), ISTAT )
          NSIZE=NSIZE+SIZE(VMAT_SCALA)
          ALLOCATE(VMAT(1,1))
       ENDIF
#ifdef single_prec_bse
       IF (IO%IU0>=0) WRITE(IO%IU0,'(A,F8.3,A,I7)') ' BSE (scaLAPACK) single prec attempting allocation of',1.0_q/1E9*NP*NQ*4*nd ,' Gbyte  rank=',NCV
       IF (IO%IU6>=0) WRITE(IO%IU6,'(A,F8.3,A,I7)') ' BSE (scaLAPACK) single prec attempting allocation of',1.0_q/1E9*NP*NQ*4*nd ,' Gbyte  rank=',NCV
       CALL REGISTER_ALLOCATE(4._q*ndata* NSIZE, "bse")
#endif

#ifdef double_prec_bse
       IF (IO%IU0>=0) WRITE(IO%IU0,'(A,F8.3,A,I7)') ' BSE (scaLAPACK) attempting allocation of',1.0_q/1E9*NP*NQ*8*nd ,' Gbyte  rank=',NCV
       IF (IO%IU6>=0) WRITE(IO%IU6,'(A,F8.3,A,I7)') ' BSE (scaLAPACK) attempting allocation of',1.0_q/1E9*NP*NQ*8*nd ,' Gbyte  rank=',NCV
       CALL REGISTER_ALLOCATE(8._q*ndata* NSIZE, "bse")
#endif

       IF (LQP) THEN
          CALL INIT_scala_DESC( WDES%COMM_INTER, NAB, POT_DESC, POT_GSD )
          CALL BLACS_GRIDINFO(POT_DESC(CTXT_), NPROW, NPCOL,MYROW,MYCOL)
          NP = NUMROC(NAB,POT_DESC(MB_),MYROW,0,NPROW)
          NQ = NUMROC(NAB,POT_DESC(NB_),MYCOL,0,NPCOL)
          ALLOCATE(AVpW_SCALA( MAX(NP*NQ,1)), STAT = ISTAT )
          CALL CHECK_ALLOCATION_ERROR( "AVpW_SCALA", MAX(NP*NQ,1), ISTAT )

          ALLOCATE( VMAT_SCALA(MAX(NP*NQ,1)), STAT = ISTAT )
          CALL CHECK_ALLOCATION_ERROR( "VMAT_SCALA", MAX(NP*NQ,1), ISTAT )

          NSIZE=SIZE(AVpW_SCALA)+SIZE(VMAT_SCALA)
          ALLOCATE( AVpW(1,1), VMAT(1,1) ) 
          IF (ANTIRES>=2) THEN
             ALLOCATE(BVpW_SCALA(MAX(NP*NQ,1)), STAT = ISTAT)
             CALL CHECK_ALLOCATION_ERROR( "BVpW_SCALA", MAX(NP*NQ,1), ISTAT )
             ALLOCATE(BVpW(1,1))
             NSIZE=NSIZE+SIZE(BVpW_SCALA)
          ENDIF
#ifdef single_prec_bse
          IF (IO%IU0>=0) WRITE(IO%IU0,'(A,F8.3,A,I7)') ' BSE attempting extended potential matrix allocation of ',1.0_q/1E9*NP*NQ*4*(nd+ndata) ,' Gbyte  rank=',NAB
          IF (IO%IU6>=0) WRITE(IO%IU6,'(A,F8.3,A,I7)') ' BSE attempting allocation of',1.0_q/1E9*NP*NQ*4*(nd+ndata) ,' Gbyte  rank=',NAB
          CALL REGISTER_ALLOCATE(4._q*ndata* NSIZE , "pot")
#endif
#ifdef double_prec_bse
          IF (IO%IU0>=0) WRITE(IO%IU0,'(A,F8.3,A,I7)') ' BSE attempting extended potential matrix allocation of',1.0_q/1E9*NP*NQ*8*(nd+ndata) ,' Gbyte  rank=',NAB
          IF (IO%IU6>=0) WRITE(IO%IU6,'(A,F8.3,A,I7)') ' BSE attempting extended potential matrix allocation of',1.0_q/1E9*NP*NQ*8*(nd+ndata) ,' Gbyte  rank=',NAB
          CALL REGISTER_ALLOCATE(8._q*ndata* NSIZE , "pot")
#endif
       ENDIF
#endif
!==========================================================================
    ELSE ! initialisation for the non-scaLAPACK-aware case
!==========================================================================
       ALLOCATE(AMAT(NCV, NCV), STAT = ISTAT )
       CALL CHECK_ALLOCATION_ERROR( "AMAT", NCV*NCV, ISTAT )
       ALLOCATE(R(NCV),  AMAT_SCALA(1) )
       NSIZE=SIZE(AMAT)
       IF (ANTIRES>=2) THEN
          ALLOCATE(BMAT(NCV, NCV), STAT = ISTAT )
          CALL CHECK_ALLOCATION_ERROR( "BMAT", NCV*NCV, ISTAT )

          ALLOCATE( BMAT_SCALA(1), RB(NCV))
          NSIZE=NSIZE+SIZE(BMAT)
       ENDIF
       IF ((.NOT. LQP) .AND. (NALPHA > 1 .OR. LALPHA))  THEN
          ALLOCATE(VMAT(NCV, NCV), STAT = ISTAT )
          CALL CHECK_ALLOCATION_ERROR( "VMAT", NCV*NCV, ISTAT )

          ALLOCATE(VMAT_SCALA(1))
          NSIZE=NSIZE+SIZE(VMAT)
       ENDIF
#ifdef single_prec_bse
       IF (IO%IU0>=0) WRITE(IO%IU0,'(A,F8.3,A,I7)') ' BSE single prec attempting allocation of',1.0_q/1E9*NCV*NCV*4*nd ,' Gbyte  rank=',NCV
       IF (IO%IU6>=0) WRITE(IO%IU6,'(A,F8.3,A,I7)') ' BSE single prec attempting allocation of',1.0_q/1E9*NCV*NCV*4*nd ,' Gbyte  rank=',NCV
       CALL REGISTER_ALLOCATE(4._q*ndata* NSIZE , "bse")
#endif

#ifdef double_prec_bse
       IF (IO%IU0>=0) WRITE(IO%IU0,'(A,F8.3,A,I7)') ' BSE attempting matrix allocation of',1.0_q/1E9*NCV*NCV*8*nd,' Gbyte  rank=',NCV
       IF (IO%IU6>=0) WRITE(IO%IU6,'(A,F8.3,A,I7)') ' BSE attempting matrix allocation of',1.0_q/1E9*NCV*NCV*8*nd,' Gbyte  rank=',NCV
       CALL REGISTER_ALLOCATE(8._q*ndata* NSIZE , "bse")
#endif
       
       IF (LQP) THEN
          ALLOCATE(AVpW(NAB,NCV), STAT = ISTAT )
          CALL CHECK_ALLOCATION_ERROR( "AVpW", NAB*NCV, ISTAT )

          ALLOCATE(VMAT(NAB,NCV), STAT = ISTAT )
          CALL CHECK_ALLOCATION_ERROR( "VMAT", NAB*NCV, ISTAT )

          NSIZE=SIZE(AVpW)+SIZE(VMAT)

          ALLOCATE( AVpW_SCALA(1), VMAT_SCALA(1))
          IF (ANTIRES>=2) THEN
             ALLOCATE(BVpW(NAB,NCV), STAT = ISTAT )
             CALL CHECK_ALLOCATION_ERROR( "BVpW", NAB*NCV, ISTAT )

             ALLOCATE(BVpW_SCALA(1))
             NSIZE=NSIZE+SIZE(BVpW)
          ENDIF
#ifdef single_prec_bse
          IF (IO%IU0>=0) WRITE(IO%IU0,'(A,F8.3,A,I7)') ' BSE single prec attempting allocation of',1.0_q/1E9*NAB*NCV*4*(nd) ,' Gbyte  rank=',NAB
          IF (IO%IU6>=0) WRITE(IO%IU6,'(A,F8.3,A,I7)') ' BSE single prec attempting allocation of',1.0_q/1E9*NAB*NCV*4*(nd) ,' Gbyte  rank=',NAB
          CALL REGISTER_ALLOCATE(4._q*ndata* NSIZE , "pot")
#endif
#ifdef double_prec_bse
          IF (IO%IU0>=0) WRITE(IO%IU0,'(A,F8.3,A,I7)') ' BSE attempting extended potential matrix allocation of',1.0_q/1E9*NAB*NCV*8*(nd),' Gbyte  rank=',NAB
          IF (IO%IU6>=0) WRITE(IO%IU6,'(A,F8.3,A,I7)') ' BSE attempting extended potential matrix allocation of',1.0_q/1E9*NAB*NCV*8*(nd),' Gbyte  rank=',NAB
          CALL REGISTER_ALLOCATE(8._q*ndata* NSIZE , "pot")
#endif
       ENDIF
    ENDIF
    ENDIF
#ifdef use_bse_te
!==========================================================================
! allocation for IBSE>0, timeevolution code
!==========================================================================
    IF (IBSE>0) THEN
       ALLOCATE( AMAT(1,1), AMAT_SCALA(1), R(NCV))
       IF (ANTIRES>=2) ALLOCATE(BMAT(1,1), BMAT_SCALA(1), RB(NCV))
       CALL TW4O_STORE_STRIP_PREPARE(WHF, IO, NCV,  LKPOINT_PARALLEL, WDES%NKPTS) 
    ENDIF
#endif
!==========================================================================

    NCBD=MAX( BD(ISP_LOW)%NGLB,BD(ISP_HIGH)%NGLB) ; NCBD4=MAX( BD(ISP_LOW)%NGLB4,BD(ISP_HIGH)%NGLB4)
    IF (IO%IU0>=0) WRITE(IO%IU0,'(A,4I5)') ' allocating two-electron 4 orbital integral table',NSTRIPV, NCBD, NSTRIPV, NCBD4

    ALLOCATE(TWOELECTRON3O( NSTRIPV, NCBD, NSTRIPV, NCBD4), STAT = ISTAT)
    CALL CHECK_ALLOCATION_ERROR_Q( "TWOELECTRON3O",NSTRIPV*NCBD*NSTRIPV*NCBD4 , ISTAT )

    CALL REGISTER_ALLOCATE(8._q*ndata* SIZE(TWOELECTRON3O) , "bse")
    IF (ANTIRES>=2) THEN
       ALLOCATE(BTWOELECTRON3O( NSTRIPV, NCBD, NSTRIPV, NCBD4), STAT = ISTAT )
       CALL CHECK_ALLOCATION_ERROR_Q( "BTWOELECTRON3O",NSTRIPV*NCBD*NSTRIPV*NCBD4 , ISTAT )

       CALL REGISTER_ALLOCATE(8._q*ndata* SIZE(BTWOELECTRON3O) , "bse")
    ENDIF
    IF (ALLOCATED(VMAT).OR. ALLOCATED(VMAT_SCALA)) THEN
       ALLOCATE(TWOE3O_HARTREE( NSTRIPV, NCBD, NSTRIPV, NCBD4), STAT = ISTAT )
       CALL CHECK_ALLOCATION_ERROR_Q( "TWOE3O_HARTREE",NSTRIPV*NCBD*NSTRIPV*NCBD4 , ISTAT )

       CALL REGISTER_ALLOCATE(8._q*ndata* SIZE(TWOE3O_HARTREE) , "bse")
    ENDIF
!==========================================================================

    CALL DUMP_ALLOCATE(IO%IU6)
    IF (IO%LOPEN) CALL WFORCE(IO%IU6)
!==========================================================================
! start cycle over values of the interaction strength in order 
! to evaluate the correlation energy within the adiabatic-connection
! formalism.
!==========================================================================
    IF (ALLOCATED(VMAT)) VMAT=0
    IF (ALLOCATED(VMAT_SCALA)) VMAT_SCALA=0

    IF (LQP) THEN
       AVpW=0; AVpW_SCALA=0
       BVpW=0; BVpW_SCALA=0
    ENDIF
!==========================================================================
! loop over coupling constant integration points alpha
    SUM_ALPHA=0.0
ccintegration: DO IALPHA=1,NALPHA
!==========================================================================
    TR_ALPHA=0.0
    NFFTW = 0
    NFLOAT4O=0 ; NFFT4O=0


    IF ( LscaLAPACKaware) THEN
       AMAT_SCALA=0
       IF (ANTIRES>=2) BMAT_SCALA=0
    ELSE
       AMAT=0
       IF (ANTIRES>=2) BMAT=0
    ENDIF

    IF (IO%IU0>=0) WRITE(IO%IU0,*) 'BSE setting up matrix'
    CALL START_TIMING("G")

    ! Tamm-Dankoff with correct properties at large w

    IF (ISP_IN<0 .OR. WHF%WDES%LNONCOLLINEAR ) THEN
       ! spin polarized calculation weight for Hartree term is one
       AHARTREE=1.0_q
    ELSE
       ! non spin polarized: weight for Hartree is two
       AHARTREE=2.0_q
    ENDIF

    ! correct properties at w=0
    ! at w=0  resonant-resonant coupling = antires-antires=resonant-antires
    ! so Hartree needs to be counted twice (this however spoils high frequency behaviour)
    IF (ANTIRES==1) AHARTREE=AHARTREE*2.0_q
!==========================================================================
! outer loop body for first set of pair states
!==========================================================================
    DO ISP=ISP_LOW,ISP_HIGH
    B1=BD(ISP)
    DO K1=1,WDES%NKPTS
       ! generate the proper descriptor for W1 wavefunctions
       CALL SETWDES(WHF%WDES,WDESK1,K1)
       IF (WHF%WDES%WTKPT(K1)==0) CYCLE

       K3 =KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K1)+WHF%WDES%VKPT(:,NQPOINT),KPOINTS_FULL)

       ! skip if there are no occupied bands
       IF (B1%VBMAX==0) CYCLE

       CALL W1_GATHER_GLB( WHF, B1%VBMIN, B1%VBMAX, ISP, W1)
       NFFTW=NFFTW+B1%VBMAX-B1%VBMIN+1

       ! generate the proper descriptor for W3 wavefunctions at k-point K3
       CALL SETWDES(WHF%WDES,WDESK3 ,K3)

       CALL W1_GATHER_GLB( WHF, B1%CBMIN, B1%CBMAX, ISP, W3)
       NFFTW=NFFTW+B1%CBMAX-B1%CBMIN+1

       
       IF (ABS(WHF%WDES%VKPT(1,NQPOINT))>1E-6_q .OR. ABS(WHF%WDES%VKPT(2,NQPOINT))>1E-6_q .OR. ABS(WHF%WDES%VKPT(3,NQPOINT))>1E-6_q) THEN
          ! set transition matrix elements
          CALL ONEELECTRON_TRANS(WHF, P, LATT_CUR, ISP, WGW, & 
               W1, K1, B1%VBMIN, B1%VBMAX, & 
               W3, K3, B1%CBMIN, B1%CBMAX, KPOINT_BSE, NFFTW )
       ENDIF


       IF (NSTRIPV == 0) NSTRIPV = 1 

       DO NPOS1=B1%VBMIN, B1%VBMAX, NSTRIPV
          CALL GWPROGRESS(IO%IU0, K1, WDES%NKPTS, NPOS1-B1%VBMIN+1, B1%VBMAX-B1%VBMIN+1 )
          NSTRIP1=MIN(B1%VBMAX+1-NPOS1,NSTRIPV)

!==========================================================================
! inner loop body for second set of pair states
!==========================================================================
          IF (LTRIANGLE) THEN
             ISP2_LOW=ISP
          ELSE
             ISP2_LOW=ISP_LOW
          ENDIF
          DO ISP2=ISP2_LOW,ISP_HIGH
          B2=BD(ISP2)
          IF (LTRIANGLE) THEN
             K2=K1-1
          ELSE
             K2=0
          ENDIF
          IF (ISP2/=ISP) K2=0
          DO
             K2=K2+1
             ! distribute the bands over k-points in a round robin fashion
             K2_DONE =0                   ! counts the number of k-points that have been collected
             K2_LOCAL=-1                  ! determine which k-point treated locally
             DO K2_COLLECT=K2, WDES%NKPTS ! loop from present K2 up to number of k-points in wavefunction descriptor
                ! generate the proper descriptor for W2 wavefunctions
                CALL SETWDES(WHF%WDES,WDESK2,K2_COLLECT)

                ! skip bands that are empty
                IF (WHF%WDES%WTKPT(K2_COLLECT)==0) CYCLE

                K2_DONE=K2_DONE+1  ! new k-point to be included

                K4=KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K2_COLLECT)+WHF%WDES%VKPT(:,NQPOINT),KPOINTS_FULL)

                IF (LKPOINT_PARALLEL) THEN
                   IF ( B2%VBMAX-B2%VBMIN >=0) &
                        CALL W1_GATHER_KSEL( WHF, B2%VBMIN, B2%VBMAX, ISP2, W2, K2_DONE)
                   NFFTW=NFFTW+(B2%CBMAX-B2%CBMIN+1)/WDES%NB_PAR
                ELSE
                   IF (W1EQUALW2) THEN
                      W2=W1
                   ELSE
                      IF ( B2%VBMAX-B2%VBMIN >=0) &
                           CALL W1_GATHER_GLB( WHF, B2%VBMIN, B2%VBMAX, ISP2, W2)
                      NFFTW=NFFTW+B2%VBMAX-B2%VBMIN+1
                   ENDIF
                ENDIF

                CALL SETWDES(WHF%WDES,WDESK4,K4)

                IF (LKPOINT_PARALLEL) THEN
                   IF ( B2%CBMAX-B2%CBMIN >=0) &
                        CALL W1_GATHER_KSEL( WHF, B2%CBMIN, B2%CBMAX, ISP2, W4, K2_DONE)
                   NFFTW=NFFTW+(B2%CBMAX-B2%CBMIN+1)/WDES%NB_PAR
                ELSE
                   IF (W1EQUALW2) THEN
                      W4=W3
                   ELSE
                      IF ( B2%CBMAX-B2%CBMIN >=0) &
                           CALL W1_GATHER_GLB( WHF, B2%CBMIN, B2%CBMAX, ISP2, W4)
                      NFFTW=NFFTW+B2%CBMAX-B2%CBMIN+1
                   ENDIF
                ENDIF

                IF (K2_DONE==WDES%NB_LOW .OR. .NOT. LKPOINT_PARALLEL) THEN
                   K2_LOCAL=K2_COLLECT
                   K4_LOCAL=K4
                ENDIF
                IF (K2_DONE==WDES%NB_PAR .OR. .NOT. LKPOINT_PARALLEL) EXIT
             ENDDO
!==========================================================================
!  set up BSE matrix as well as four point Coulomb and exchange integrals
!  for the subset of bands and k-points considered
!==========================================================================
             ! at this point each CPU holds the k-points corresponding to K2_LOCAL 
             ! in W2 (and corresponding to K2-NQ in W4)
             ! set K2 and K4 again
             ! (Mind the condition K2>=1, some nodes might skip this loop altogheter)
             K2=K2_LOCAL
             K4=K4_LOCAL
             IF (K2>=1) THEN
                CALL SETWDES(WHF%WDES,WDESK2,K2)
                CALL SETWDES(WHF%WDES,WDESK4,K4)
             ENDIF
             DO NPOS2=B2%VBMIN, B2%VBMAX, NSTRIPV
                NSTRIP2=MIN(B2%VBMAX+1-NPOS2,NSTRIPV)

                IF (K2>=1) THEN

                   TWOELECTRON3O=0
                   IF (ISP2==ISP) THEN

                     IF (LHARTREE .AND. .NOT. LTRIPLET) THEN
                       CALL TWOELECTRON4O_ACC_HARTREE(WHF, P, LATT_CUR, ISP, ISP2, WGW, & 
                         LFXC .AND. .NOT. LGWLF, .FALSE. , AHARTREE,  &
                         W1, K1, NPOS1, NSTRIP1, W2, K2, NPOS2, NSTRIP2, W3, K3, W4, K4, &
                         TWOELECTRON3O, B1%CBMIN, B1%CBMAX, B1%CBMIN4, B1%CBMAX4, B1%VBMIN, NFLOAT4O, NFFT4O, NSTRIP)

                       ! store Coulomb integrals in TWOE3O_HARTREE (contributions from variation of Hartree term)
                       IF (ALLOCATED(TWOE3O_HARTREE)) TWOE3O_HARTREE=TWOELECTRON3O
 
                       ! include head correction in the polarisation
                       CALL TWOELECTRON4O_ACC_HEAD(WHF, P, LATT_CUR, ISP, ISP2, WGW, & 
                         AHARTREE,  &
                         W1, K1, NPOS1, NSTRIP1, W2, K2, NPOS2, NSTRIP2, W3, K3, W4, K4, &
                         TWOELECTRON3O, B1%CBMIN, B1%CBMAX, B1%CBMIN4, B1%CBMAX4, NSTRIP)
                     ENDIF

                     ! copy Coulomb integrals to resonant-antiresonant part
                     IF (ANTIRES>=2) BTWOELECTRON3O=TWOELECTRON3O
                     
                     IF (LADDER .AND. (AEXX /=0 .OR. L_MODEL_HF .OR. LGWLF) ) THEN
                        ! calculate direct term from exchange (contribution to A matrix)
                        CALL TWOELECTRON4O_ACC_DIRECT(WHF, P, LATT_CUR, ISP, WGW, DELTA_COND,  &
                             W1, K1, NPOS1, NSTRIP1, W2, K2, NPOS2, NSTRIP2, W3, K3, W4, K4, &
                             TWOELECTRON3O, B1%CBMIN, B1%CBMAX, B1%CBMIN4, B1%CBMAX4, B1%VBMIN, NFLOAT4O, NFFT4O, NSTRIP)
                 
                        ! calculate resonant-antiresonant term from exchange and add it to the A matrix
                        IF (ANTIRES==1 .AND. .NOT. LEX_INTERPOLATED) &
                           CALL TWOELECTRON4O_ACC_EX(WHF, P, LATT_CUR, ISP, WGW, 0.0_q, &
                                W1, K1, NPOS1, NSTRIP1, W2, K2, NPOS2, NSTRIP2, W3, K3, W4, K4, &
                                TWOELECTRON3O, B1%CBMIN, B1%CBMAX, B1%CBMIN4, B1%CBMAX4, B1%VBMIN, NFLOAT4O, NFFT4O, NSTRIP)
                        ! calculate resonant-antiresonant term from exchange and add it to the B matrix
                        IF (ANTIRES>=2 .AND. .NOT. LEX_INTERPOLATED) &
                           CALL TWOELECTRON4O_ACC_EX(WHF, P, LATT_CUR, ISP, WGW, 0.0_q, &
                                W1, K1, NPOS1, NSTRIP1, W2, K2, NPOS2, NSTRIP2, W3, K3, W4, K4, &
                                BTWOELECTRON3O, B1%CBMIN, B1%CBMAX, B1%CBMIN4, B1%CBMAX4, B1%VBMIN, NFLOAT4O, NFFT4O, NSTRIP)
                 
                        ! calculate resonant-antiresonant term from exchange and add it to the A matrix (interpolated)
                        IF (ANTIRES==1 .AND. LEX_INTERPOLATED) &
                           CALL TWOELECTRON4O_ACC_EX_INTER(WHF, P, LATT_CUR, ISP, WGW, 0.0_q, &
                                W1, K1, NPOS1, NSTRIP1, W2, K2, NPOS2, NSTRIP2, W3, K3, W4, K4, &
                                TWOELECTRON3O, B1%CBMIN, B1%CBMAX, B1%CBMIN4, B1%CBMAX4, B1%VBMIN, & 
                                NFLOAT4O, NFFT4O, NSTRIP, LATT_CUR%B)
                        ! calculate resonant-antiresonant term from exchange and add it to the B matrix (interpolated)
                        IF (ANTIRES>=2 .AND. LEX_INTERPOLATED) &
                           CALL TWOELECTRON4O_ACC_EX_INTER(WHF, P, LATT_CUR, ISP, WGW, 0.0_q, &
                                W1, K1, NPOS1, NSTRIP1, W2, K2, NPOS2, NSTRIP2, W3, K3, W4, K4, &
                                BTWOELECTRON3O, B1%CBMIN, B1%CBMAX, B1%CBMIN4, B1%CBMAX4, B1%VBMIN, & 
                                NFLOAT4O, NFFT4O, NSTRIP, LATT_CUR%B)
                     ENDIF ! condition LADDER=.TRUE. .AND. (AEXX /=0 .OR. L_MODEL_HF .OR. LGWLF)
       
                   ELSE ! condition ISP2/=ISP
                       ! spin-up spin-down coupling only via Hartree like Coulomb integrals
                       IF (LHARTREE .AND. .NOT. LTRIPLET) &
                       CALL TWOELECTRON4O_ACC_HARTREE_SFLIP(WHF, P, LATT_CUR, ISP, ISP2, WGW, & 
                          LFXC .AND. .NOT. LGWLF, .FALSE. ,AHARTREE,  &
                          W1, K1, NPOS1, NSTRIP1, W2, K2, NPOS2, NSTRIP2, W3, K3, W4, K4, &
                          TWOELECTRON3O, B1%CBMIN, B1%CBMAX, B1%VBMIN, & 
                          B2%CBMIN, B2%CBMIN4, B2%CBMAX4, B2%VBMIN, NFLOAT4O, NFFT4O, NSTRIP)

                       IF (ALLOCATED(TWOE3O_HARTREE)) TWOE3O_HARTREE=TWOELECTRON3O

                       ! include head correction only in the polarisation
                       ! propagator, not in the Coulomb integrals
                       CALL TWOELECTRON4O_ACC_HEAD(WHF, P, LATT_CUR, ISP, ISP2, WGW, & 
                         AHARTREE,  &
                         W1, K1, NPOS1, NSTRIP1, W2, K2, NPOS2, NSTRIP2, W3, K3, W4, K4, &
                         TWOELECTRON3O, B1%CBMIN, B1%CBMAX, B2%CBMIN4, B2%CBMAX4, NSTRIP)

                       ! copy Hartree-term to resonant-antiresonant part
                       IF (ANTIRES>=2) BTWOELECTRON3O=TWOELECTRON3O

                   ENDIF ! condition ISP2==ISP 
                ENDIF ! condition K2>=1

                IF (IBSE==0) THEN
!==========================================================================
! now store the calculated contributions in the final arrays
!               IBSE=0, scaLAPACK-aware storage
!==========================================================================
                IF (LscaLAPACKaware) THEN
#ifdef scaLAPACK
                   CALL TWOELECTRON4O_STORE_SCALA(WHF, ISP, ISP2, ISP_LOW, &
                        K1, NPOS1, NSTRIP1, K2, NPOS2, NSTRIP2, K3, K4, &
                        TWOELECTRON3O, B1, B2, AMAT_SCALA, NCV, BSE_INDEX,1, ALPHAV(IALPHA),SCALE_PERT, BSE_DESC)
                   IF (IALPHA==1) THEN
                     IF (LQP) THEN
                        ! 
                        CALL TWOELECTRON4O_STORE_SCALA(WHF, ISP, ISP2, ISP_LOW, &
                             K1, NPOS1, NSTRIP1, K2, NPOS2, NSTRIP2, K3, K4, &
                             TWOE3O_HARTREE, B1, B2, VMAT_SCALA, NAB, POT_INDEX,0,1.0_q,SCALE_PERT,POT_DESC)
                        CALL TWOELECTRON4O_STORE_SCALA(WHF, ISP, ISP2, ISP_LOW, &
                             K1, NPOS1, NSTRIP1, K2, NPOS2, NSTRIP2, K3, K4, &
                             TWOELECTRON3O, B1, B2, AVpW_SCALA, NAB, POT_INDEX,0, 1.0_q,SCALE_PERT,POT_DESC)
                     ELSE
                        IF (ALLOCATED(VMAT_SCALA)) &
                          CALL TWOELECTRON4O_STORE_SCALA(WHF, ISP, ISP2, ISP_LOW, &
                             K1, NPOS1, NSTRIP1, K2, NPOS2, NSTRIP2, K3, K4, &
                             TWOE3O_HARTREE, B1, B2, VMAT_SCALA, NCV, BSE_INDEX,0,1.0_q,SCALE_PERT,BSE_DESC)
                     ENDIF
                   ENDIF

                   IF (ANTIRES>=2) THEN 
                     CALL TWOELECTRON4O_STORE_SCALA(WHF, ISP, ISP2, ISP_LOW, &
                          K1, NPOS1, NSTRIP1, K2, NPOS2, NSTRIP2, K3, K4, &
                          BTWOELECTRON3O, B1, B2, BMAT_SCALA, NCV, BSE_INDEX,0, ALPHAV(IALPHA),SCALE_PERT,BSE_DESC)
                     IF (IALPHA==1 .AND. LQP) &
                       CALL TWOELECTRON4O_STORE_SCALA(WHF, ISP, ISP2, ISP_LOW, &
                            K1, NPOS1, NSTRIP1, K2, NPOS2, NSTRIP2, K3, K4, &
                            BTWOELECTRON3O, B1, B2, BVpW_SCALA, NAB, POT_INDEX,0, 1.0_q, SCALE_PERT,POT_DESC)
                   ENDIF
                
#endif
!==========================================================================
!               IBSE=0, non scaLAPACK-aware storage
!==========================================================================
                ELSE IF (K2>=1) THEN
                    CALL TWOELECTRON4O_STORE(WHF, ISP, ISP2, ISP_LOW,  &
                         K1, NPOS1, NSTRIP1, K2, NPOS2, NSTRIP2, K3, K4, &
                         TWOELECTRON3O, B1, B2, AMAT, BSE_INDEX,1, ALPHAV(IALPHA), SCALE_PERT)

                    IF (IALPHA==1) THEN
                      IF (LQP) THEN
                        IF (ALLOCATED(VMAT)) & 
                         CALL TWOELECTRON4O_STORE_POT(WHF, ISP, ISP2, ISP_LOW,  &
                              K1, NPOS1, NSTRIP1, K2, NPOS2, NSTRIP2, K3, K4, &
                              TWOE3O_HARTREE, B1, B2, VMAT, BSE_INDEX, POT_INDEX, 1.0_q, LTRIANGLE)
                         CALL TWOELECTRON4O_STORE_POT(WHF, ISP, ISP2, ISP_LOW,  &
                              K1, NPOS1, NSTRIP1, K2, NPOS2, NSTRIP2, K3, K4, &
                              TWOELECTRON3O, B1, B2, AVpW, BSE_INDEX, POT_INDEX, 1.0_q, LTRIANGLE)
                      ELSE
                        IF (ALLOCATED(VMAT)) & 
                          CALL TWOELECTRON4O_STORE(WHF, ISP, ISP2, ISP_LOW,  &
                              K1, NPOS1, NSTRIP1, K2, NPOS2, NSTRIP2, K3, K4, &
                              TWOE3O_HARTREE, B1, B2, VMAT, BSE_INDEX,0, 1.0_q, SCALE_PERT)
                      ENDIF
                    ENDIF 

                    IF (ANTIRES>=2) THEN
                      CALL TWOELECTRON4O_STORE(WHF, ISP, ISP2, ISP_LOW,  &
                           K1, NPOS1, NSTRIP1, K2, NPOS2, NSTRIP2, K3, K4, &
                           BTWOELECTRON3O, B1, B2, BMAT, BSE_INDEX,0, ALPHAV(IALPHA), SCALE_PERT)
                      IF (IALPHA==1 .AND. LQP) &
                        CALL TWOELECTRON4O_STORE_POT(WHF, ISP, ISP2, ISP_LOW,  &
                             K1, NPOS1, NSTRIP1, K2, NPOS2, NSTRIP2, K3, K4, &
                             BTWOELECTRON3O, B1, B2, BVpW, BSE_INDEX, POT_INDEX, 1.0_q, LTRIANGLE)
                    ENDIF
                ENDIF ! (LscaLAPACKaware)
                ENDIF ! (IBSE==0)
#ifdef use_bse_te
!==========================================================================
!               IBSE>0, time-evolution code
!==========================================================================
                IF (IBSE>0) THEN
                   CALL TW4O_STORE_STRIP(WHF, ISP, ISP2, ISP_LOW, &
                        K1, NPOS1, NSTRIP1, K2, NPOS2, NSTRIP2, K3, K4, &
                        TWOELECTRON3O, B1, B2, BSE_INDEX )
                ENDIF
#endif
             ENDDO
!==========================================================================
!  end of loop body
!==========================================================================
             K2=K2_COLLECT
             IF (K2>=WDES%NKPTS) EXIT
          ENDDO
          ENDDO
       ENDDO
    ENDDO
    ENDDO

    CALLMPI( M_sum_d(WGW%COMM_INTER, NFLOAT4O, 1))
    CALLMPI( M_sum_d(WGW%COMM_INTER, NFFT4O  , 1))

10  FORMAT(" BLAS level 3 operations / number of FFT's:"/ &
         " number of FFTs for wave wavefunctions          ",F10.0," fft"/ & 
         " number of operations in four-orbital integrals ",F10.2," Gflops, ",F10.0," fft")

    IF (IO%IU6>0) THEN
       WRITE(IO%IU6,10) NFFTW,NFLOAT4O/1E9, NFFT4O
    ENDIF

    IF (IALPHA==NALPHA) THEN 
       !deallocate wavefunctions only if you are evaluating the last value for 
       !the coupling constant integration
       ! ok finished get rid of WPOTH
       IF (ASSOCIATED(WPOTH)) CALL DESTROY_WPOT_HANDLE(WPOTH)
!==========================================================================
! deallocation
!==========================================================================
       DO N=1,ALLOC_CONDUCTION_BANDS
          CALL DELWAV(W3(N) ,.TRUE.)
          IF (.NOT. W1EQUALW2) CALL DELWAV(W4(N) ,.TRUE.)
       ENDDO

       DO N=1,ALLOC_VALENCE_BANDS
          CALL DELWAV(W1(N) ,.TRUE.)
          IF (.NOT. W1EQUALW2) CALL DELWAV(W2(N) ,.TRUE.)
       ENDDO
       
       CALL DEREGISTER_ALLOCATE(8._q*ndata* SIZE(TWOELECTRON3O) , "bse")
       DEALLOCATE(TWOELECTRON3O)
       IF (ANTIRES>=2) THEN
          CALL DEREGISTER_ALLOCATE(8._q*ndata* SIZE(BTWOELECTRON3O) , "bse")
          DEALLOCATE(BTWOELECTRON3O)
       ENDIF
       IF (ALLOCATED(TWOE3O_HARTREE)) THEN
          CALL DEREGISTER_ALLOCATE(8._q*ndata* SIZE(TWOE3O_HARTREE) , "bse")
          DEALLOCATE(TWOE3O_HARTREE)
       ENDIF
       DEALLOCATE(W1,W2,W3,W4)

! deallocate orbitals
! this is usually allowed but obviously an issue if we desire correlation
! energies, since the routine is called repeatedly
! deallocation frees precious memory for the scaLAPACK calls
       IF (.NOT. LCORBSE .AND. (.NOT. LQP)) THEN
          CALL DEALLOCW_CW(WHF)
       ENDIF
    ENDIF
!==========================================================================
! BSE  redistribution of matrix elements
!==========================================================================
    CALL STOP_TIMING("G",IO%IU6,'BSESET')

    IF (IO%IU0>=0) WRITE(IO%IU0,*)
    IF (IO%IU0>=0) WRITE(IO%IU0,*) 'BSE redistributing all elements'

    IF ((IBSE==0).AND.(.NOT. LscaLAPACKaware)) THEN
#ifdef single_prec_bse
#ifdef gammareal
       CALLMPI( M_sum_single(WGW%COMM_INTER, AMAT  , SIZE(AMAT)))
       IF (ANTIRES>=2) THEN
          CALLMPI( M_sum_single(WGW%COMM_INTER, BMAT  , SIZE(BMAT)))
       ENDIF
       IF (ALLOCATED(VMAT)) CALLMPI( M_sum_single(WGW%COMM_INTER, VMAT, SIZE(VMAT)))
       IF (ALLOCATED(AVpW)) CALLMPI( M_sum_single(WGW%COMM_INTER, AVpW, SIZE(AVpW)))
       IF (ALLOCATED(BVpW)) CALLMPI( M_sum_single(WGW%COMM_INTER, BVpW, SIZE(BVpW)))
#else 
       CALLMPI( M_sum_single(WGW%COMM_INTER, AMAT  , 2*SIZE(AMAT)))
       IF (ANTIRES>=2) THEN
           CALLMPI( M_sum_single(WGW%COMM_INTER, BMAT  , 2*SIZE(BMAT)))
       ENDIF 
       IF (ALLOCATED(VMAT)) CALLMPI( M_sum_single(WGW%COMM_INTER, VMAT, 2*SIZE(VMAT)))
       IF (ALLOCATED(AVpW)) CALLMPI( M_sum_single(WGW%COMM_INTER, AVpW, 2*SIZE(AVpW)))
       IF (ALLOCATED(BVpW)) CALLMPI( M_sum_single(WGW%COMM_INTER, BVpW, 2*SIZE(BVpW)))
#endif
#endif

#ifdef double_prec_bse
       CALLMPI( M_sum_g(WGW%COMM_INTER, AMAT  , SIZE(AMAT)))
       IF (ANTIRES>=2) THEN
          CALLMPI( M_sum_g(WGW%COMM_INTER, BMAT  , SIZE(BMAT)))
       ENDIF
       IF (ALLOCATED(VMAT)) THEN
          CALLMPI( M_sum_g(WGW%COMM_INTER, VMAT, SIZE(VMAT)))
       ENDIF
       IF (ALLOCATED(AVpW)) THEN
          CALLMPI( M_sum_g(WGW%COMM_INTER, AVpW, SIZE(AVpW)))
       ENDIF
       IF (ALLOCATED(BVpW)) THEN
          CALLMPI( M_sum_g(WGW%COMM_INTER, BVpW, SIZE(BVpW)))
       ENDIF
#endif
    ENDIF

    CALL STOP_TIMING("G",IO%IU6,'BSESUM')

    IF (IO%LOPEN) CALL WFORCE(IO%IU6)
!==========================================================================
! at this point the BSE matrix has been set up
! for the plasmon formulate of the correlation energy we need
! the sum of the diagonal components of the A matrix
!==========================================================================
    IF (IBSE==0) THEN
    CSUM_A=(0.0_q,0.0_q)
    IF (.NOT. LADDER) THEN 
       IF (LscaLAPACKaware) THEN
#ifdef scaLAPACK
          CALL TWOELECTRON4O_GET_DIAG_SCALA( AMAT_SCALA, NCV, CSUM_A, WHF%WDES%COMM_INTER, BSE_DESC)
#endif
       ELSE
          CALL TWOELECTRON4O_GET_DIAG(AMAT, NCV, CSUM_A)
       ENDIF
       IF (REAL(CSUM_A) /= REAL(CSUM_A) ) WRITE (*,*) 'WARNING: diag(A) contains NaN elements'
    ENDIF

    ! store A-B in A, and A+B in B
    IF (ANTIRES>=2) THEN
!==========================================================================
!  full beyond TDA case   (store A+B in A, and A-B in B)
!==========================================================================
       IF (LscaLAPACKaware) THEN
          DO K1=1,SIZE(AMAT_SCALA)
             TMP=AMAT_SCALA(K1)
             AMAT_SCALA(K1)=TMP+BMAT_SCALA(K1)
             BMAT_SCALA(K1)=TMP-BMAT_SCALA(K1)
          ENDDO
       ELSE
          DO K1=1,NCV
             DO K2=K1,NCV
                TMP=AMAT(K1,K2)
                AMAT(K1,K2)=TMP+BMAT(K1,K2)
                BMAT(K1,K2)=TMP-BMAT(K1,K2)
                ! for AMAT the lower triangle needs to be properly set up
                ! since below we multiply the matrix from the left and right
                AMAT(K2,K1)=GCONJG(AMAT(K1,K2))
                ! strictly speaking it suffices to set the upper triangle of B
                ! just to be sure, set the lower triangle as well
                BMAT(K2,K1)=GCONJG(BMAT(K1,K2))
             ENDDO
          ENDDO
       ENDIF

       ! diagonalize A-B 
       IF (.NOT. L2ORDER) CALL DIAG_BSE_MATRIX( NCV, BMAT, BMAT_SCALA, RB, LscaLAPACKaware, WDES%COMM_INTER, IO%IU0, BSE_DESC)
       
       CALL STOP_TIMING("G",IO%IU6,'BSEDIAG1')
       IF (IO%IU6>=0) CALL WFORCE(IO%IU6)
       ! calculate  (A-B)^1/2 (A+B) (A-B)^1/2
       IF ( LscaLAPACKaware) THEN
#ifdef scaLAPACK
          CALL TDA_BSE_SCA( NCV, AMAT_SCALA, BMAT_SCALA, RB , IO%IU0, BSE_DESC)
#endif
       ELSE
          CALL TDA_BSE( NCV, AMAT, BMAT, RB , IO%IU0)
       ENDIF
       CALL STOP_TIMING("G",IO%IU6,'SQUARE')

       ! diagonalize (A-B)^1/2 (A+B) (A-B)^1/2
       IF (.NOT. L2ORDER) CALL DIAG_BSE_MATRIX( NCV, AMAT, AMAT_SCALA, R, LscaLAPACKaware, WDES%COMM_INTER, IO%IU0, BSE_DESC)

       CALL STOP_TIMING("G",IO%IU6,'BSEDIAG2')
       IF (IO%IU6>=0) CALL WFORCE(IO%IU6)
       ! eigenvalues of original problem
       CALL SQRT_EIGENVALUE_INPLACE(R, IO%IU0)

       IF (IO%IU6>=0) CALL WFORCE(IO%IU6)

       ! compute the correlation energy using the plasmon formula only when
       ! LAMBDA=1.0 is set in the INCAR file (and there are no exchange terms)git 
       IF (ALPHAV(IALPHA)==1.0_q .AND. NALPHA==1 .AND. (.NOT. LADDER) .AND. (.NOT. LQP) ) THEN
           IF (IO%IU0>=0) WRITE(IO%IU0,'(A,1F20.10,A,I6)') ' plasmon correlation energy ', &
                              0.5*(SUM(R)-REAL(CSUM_A,q))*KPOINTS_ORIG%WTKPT(NQPOINT)/SCALE_PERT**2, ' @ NQ ', NQPOINT
           IF (IO%IU6>=0) WRITE(IO%IU6,'(A,1F20.10,A,I6)') ' plasmon correlation energy ', &
                              0.5*(SUM(R)-REAL(CSUM_A,q))*KPOINTS_ORIG%WTKPT(NQPOINT)/SCALE_PERT**2, ' @ NQ ', NQPOINT
       ENDIF

       ! left multiply by (A-B)^1/2 to obtain eigenvectors 
       IF ( LscaLAPACKaware) THEN
#ifdef scaLAPACK
          IF (LQP .OR. LCORBSE) THEN
             IF (.NOT. L2ORDER) CALL BSE_EIGENVECTORS_FULL_SCA ( NCV, BMAT_SCALA, RB, AMAT_SCALA, R, IO%IU0, BSE_DESC )
          ELSE
             ! to calculate optical properties we need only eigenvectors of squared problem
             IF (.NOT. L2ORDER) CALL TDA_BSE_EIGENVECTORS_SCA( NCV, BMAT_SCALA, RB, AMAT_SCALA, R, IO%IU0, BSE_DESC )
          ENDIF
#endif
       ELSE
          IF (LQP .OR. LCORBSE) THEN
             IF (.NOT. L2ORDER) CALL BSE_EIGENVECTORS_FULL( NCV, BMAT, RB, AMAT, R, IO%IU0)
          ELSE
             ! to calculate optical properties we need only eigenvectors of squared problem
             IF (.NOT. L2ORDER) CALL TDA_BSE_EIGENVECTORS( NCV, BMAT, RB, AMAT, R, IO%IU0)
          ENDIF
       ENDIF

       CALL STOP_TIMING("G",IO%IU6,'EIGVEC')
       IF (IO%IU6>=0) CALL WFORCE(IO%IU6)
!========================================================================
! now call subroutine to evaluate the density matrix P_alpha = |X,Y><X,Y|
! which is then contracted with the Hartree term. This is implemented in 
! AC_CORRE(_SCA) and needs to call BSE_EIGENVECTORS_FULL routines.
! ACQ_CORRE(_SCA) instead uses the quantity Q_alpha=(X_a+Y_a)(X_a+Y_a)T-E
! where X_a==X_alpha T is transposition E the identity matrix
!========================================================================
       ! compute AC-correlation energy if there is more than one integration point 
       IF(NALPHA > 1 ) THEN
         IF ( LscaLAPACKaware) THEN
#ifdef scaLAPACK
           CALL AC_CORRE_SCA(NCV, AMAT_SCALA, BMAT_SCALA, VMAT_SCALA, TR_ALPHA, WHF%WDES%COMM_INTER, BSE_DESC) 
#endif
         ELSE 
           CALL AC_CORRE (NCV, AMAT, BMAT, VMAT, TR_ALPHA)
         ENDIF
         SUM_ALPHA=SUM_ALPHA+(TR_ALPHA*W_GL(IALPHA))
       ENDIF

       CALL STOP_TIMING("G",IO%IU6,'CORREN')
       IF (IO%IU6>=0) CALL WFORCE(IO%IU6)
 
    ELSE

!==========================================================================
! simple TDA standard algorithms
!==========================================================================

       CALL DIAG_BSE_MATRIX( NCV, AMAT, AMAT_SCALA, R, LscaLAPACKaware, WDES%COMM_INTER, IO%IU0, BSE_DESC )
       CALL STOP_TIMING("G",IO%IU6,'BSEDIAG')
    ENDIF ! ANTIRES>=2

    IF (IO%IU6>=0) CALL WFORCE(IO%IU6)
!================================================================================
! self-energy evaluation at each q-point: for the gammareal case this includes a
! loop over the number of self-consistent  iterations, for the complex case this
! happens in chi.F because different q-points have to be summed over.
!================================================================================
    IF (LQP) THEN
      IF(LscaLAPACKaware) THEN
#ifdef scaLAPACK
         IF (L2ORDER) THEN 
            CALL SELFEN_2ndO_SCALA(BSE_INDEX, POT_INDEX, BD, S2E, WHF, EFERMI, VMAT_SCALA, &
                 NQPOINT, CELTOT_HK+CELTOT_X , AHARTREE, IO%IU6, IO%IU0, SEFROMBSE, POT_DESC)
         ELSE
            CALL vtutor%bug("internal error in BSE: self-energy calculation not supported in " &
               // "scaLAPACK version", __FILE__, __LINE__)
         ENDIF
#endif
      ELSE !non-scaLAPACK case
         IF (L2ORDER) THEN 
            CALL SELFEN_2ndO(BSE_INDEX, POT_INDEX, BD, S2E, WHF, VMAT, AVpW, BVpW, NQPOINT, &
                 LATT_CUR, AHARTREE, IO%IU6, IO%IU0, SEFROMBSE, LADDER)
         ELSE
! iterateG allows to iterate the Green's function inside right here (only works for Gamma only)
!#define iterateG
#if defined(gammareal)  && defined(iterateG)
            DO NELM =1, NELMGW
               IF(ANTIRES>0) CALL SELFEN_FULL(BSE_INDEX, POT_INDEX, BD, S2E, WHF, AMAT, BMAT, VMAT, &
                    AVpW, BVpW, R, NQPOINT, LATT_CUR, AHARTREE, IO%IU6, IO%IU0, &
                    SEFROMBSE, LADDER, LTCTC, .TRUE.)

               ! clean up degenerate orbitals
               CALL CLEANUP_CELEN_SELFE(W, NBANDSGW, SEFROMBSE%REAL_PART)
               CALL CLEANUP_CELEN_SELFE(W, NBANDSGW, SEFROMBSE%IMAG_PART)
               
               ! now perform a Hilbert transform to obtain self-energy from it's spectral representation
               CALL SE_HILBERT_TRANSF(SEFROMBSE, EFERMI, SHIFT, CELTOT_X, CELTOT_HK)
               
               ! compute the quasi-particle energy for state |n> as energy where
               ! Re(Sigma(w)) and y=w intersect
               IF (IO%IU6 >=0) WRITE (IO%IU6,'(//A, I2)') ' QP shifts <psi_nk| G(iteration)W_0 |psi_nk>: iteration', NELM
               CALL QPE_FROMSIGMA(W, SEFROMBSE, EFERMI, CELTOT_X, CELTOT_HK, NELM, IO%IU0, IO%IU6) 
               IF (IO%IU6 >=0) WRITE (IO%IU6,130) 
            ENDDO
            ! finally call SELFEN_FULL once more so that top level routine can still continue
#endif
            IF(ANTIRES>0) CALL SELFEN_FULL(BSE_INDEX, POT_INDEX, BD, S2E, WHF, AMAT, BMAT, VMAT, &
                 AVpW, BVpW, R, NQPOINT, LATT_CUR, AHARTREE, IO%IU6, IO%IU0, &
                 SEFROMBSE, LADDER, LTCTC)
         ENDIF
      ENDIF  !LscaLAPACKaware
    ENDIF ! LQP
    ENDIF ! IBSE == 0
!================================================================================
  ENDDO ccintegration ! over values of the coupling constant (IALPHA) 
!================================================================================
 
130 FORMAT (5X, //, &
         &'----------------------------------------------------', &
         &'----------------------------------------------------'//)

    CALL STOP_TIMING("G",IO%IU6,'SELFEN')
    IF (IO%IU6>=0) CALL WFORCE(IO%IU6)
!================================================================================
! AC-correlation energy from coupling constant integration
!================================================================================
    IF (IBSE==0) THEN
    IF (NALPHA > 1)  THEN
       CORREP=0.5*(REAL(SUM_ALPHA,q))*KPOINTS_ORIG%WTKPT(NQPOINT)/SCALE_PERT**2 

      IF (IO%IU0>=0) THEN
         WRITE(IO%IU0,'(1X,A,I4,F20.10)') 'AC-correlation energy @ NQ', NQPOINT, CORREP
         WRITE(IO%IU6,'(1X,A,I4,F20.10)') 'AC-correlation energy @ NQ', NQPOINT, CORREP
      ENDIF
    ENDIF
    IF ( LALPHA ) THEN
       CORREP=0.5*(REAL(TR_ALPHA,q))*KPOINTS_ORIG%WTKPT(NQPOINT)/SCALE_PERT**2 

       IF (IO%IU0>=0) THEN
          WRITE(IO%IU0,'(1X,A,F20.10,A,I4,A,F10.4)') 'AC-energy integrand',CORREP, 'NQ', NQPOINT,'Lambda', LAMBDA
          WRITE(IO%IU0,'(1X,A,F20.10)') 'Tr_lambda',TR_ALPHA 
          
          WRITE(IO%IU6,'(1X,A,F20.10,A,I4,A,F10.4)') 'AC-energy integrand',CORREP, 'NQ', NQPOINT,'Lambda', LAMBDA
          WRITE(IO%IU6,'(1X,A,F20.10)') 'Tr_lambda',TR_ALPHA 
       ENDIF
    ENDIF ! LALPHA
    ENDIF ! IBSE == 0

#ifdef use_bse_te
!==========================================================================
! use time evolution algorithm 
!==========================================================================
    IF (IBSE>0) THEN
       CALL TW4O_CHECK_STRIP(WHF)
       CALL INIT_BSE_OSZI_STRENGTH(WHF, LATT_CUR, ISP_LOW, ISP_HIGH, BD, BSE_INDEX)
       CALL CALCULATE_BSE_TIME_EVOLUTION(SHIFT,LATT_CUR%OMEGA, &
            WDES%ISPIN, WDES%NKPTS, NCV, WDES%WTKPT, OMEGAMAX, LADDER .OR. (LHARTREE .AND. .NOT. LTRIPLET), ANTIRES==2, NELMGW, IO)
    ENDIF
#endif

!==========================================================================
! BSEFATBANDS
!==========================================================================
    IF (IBSE==0) THEN
    IF (IO%IU0>=0) &
       WRITE(IO%IU0,"(A17,I6,A13)") ' BSE writing out ',NBSEEIG,' eigenvectors'

    IF ( LscaLAPACKaware) THEN
#ifdef scaLAPACK
       CALL WRITE_BSE_SCALA(WHF, ISP_LOW, ISP_HIGH, &
            BD, BSE_INDEX, R, AMAT_SCALA, NCV, NBSEEIG, IO%IU6, BSE_DESC)
#endif
    ELSE
       IF (IO%IU6>=0) THEN
         CALL WRITE_BSE_NOTSCALA(WHF, NQPOINT, KPOINT_BSE, ISP_LOW, ISP_HIGH,  &
         BD, BSE_INDEX, R, AMAT, NBSEEIG, IO%IU6)
       ENDIF
    ENDIF
    ENDIF

!==========================================================================
! calculate oscillator strength
!==========================================================================

    IF (IBSE==0 .AND. (.NOT. LCORBSE) .AND. (.NOT. LQP)) THEN
    ! the optical properties calculation is no longer compatible with the rest
    ! of the code when QP energies are computed because SET_BAND_PARAMETERS
    ! includes a few unoccupied states in the VB.
      ISP=1
      IF (IO%IU0>=0) WRITE(IO%IU0,*) 'BSE calculating oscillator strength'
      IF ( LscaLAPACKaware) THEN
#ifdef scaLAPACK
         CALL CALCULATE_BSE_OSZI_SCALA(WHF, NQPOINT, KPOINT_BSE, LATT_CUR, ISP_LOW, ISP_HIGH, & 
              SHIFT, BD, BSE_INDEX, R, AMAT_SCALA, NCV, ANTIRES, NEDOS, OMEGAMAX, IO%IU6, BSE_DESC)
#endif
      ELSE
         IF (IO%IU0>=0) WRITE(IO%IU0,*) 'calling CALCULATE_BSE_OSZI_STRENGTH'
         CALL CALCULATE_BSE_OSZI_STRENGTH(WHF, NQPOINT, KPOINT_BSE, LATT_CUR, ISP_LOW, ISP_HIGH, & 
              SHIFT, BD, BSE_INDEX, R, AMAT, ANTIRES, NEDOS, OMEGAMAX, IO%IU6)
      ENDIF
      CALL STOP_TIMING("G",IO%IU6,'BSEOSZI')
      IF (IO%IU6>=0) CALL WFORCE(IO%IU6)
    ENDIF
!==========================================================================
! deallocation 
!==========================================================================
    CALL DEALLOCATE_CDER_BETWEEN_STATES_Q
    CALL DEALLOCATE_BSE_MATRIX_INDEX( BSE_INDEX )
    IF ( LscaLAPACKaware) THEN
#ifdef scaLAPACK
       NSIZE=SIZE(AMAT_SCALA)
       DEALLOCATE(AMAT, AMAT_SCALA, R)
       IF (ANTIRES>=2) THEN
          NSIZE=NSIZE+SIZE(BMAT_SCALA)
          DEALLOCATE(BMAT, BMAT_SCALA, RB)
       ENDIF
       IF ((.NOT. LQP) .AND. (NALPHA > 1 .OR. LALPHA)) THEN
          NSIZE=NSIZE+SIZE(VMAT_SCALA)
          DEALLOCATE(VMAT_SCALA, VMAT)
       ENDIF
#ifdef single_prec_bse
       CALL DEREGISTER_ALLOCATE(4._q*ndata* NSIZE , "bse")
#endif
#ifdef double_prec_bse
       CALL DEREGISTER_ALLOCATE(8._q*ndata* NSIZE , "bse")
#endif


       IF (LQP) THEN
          NSIZE=SIZE(AVpW_SCALA)+SIZE(VMAT_SCALA)+SIZE(BVpW_SCALA)
          DEALLOCATE(VMAT_SCALA, VMAT, AVpW_SCALA, AVpW, BVpW_SCALA, BVpW)
#ifdef single_prec_bse
          CALL DEREGISTER_ALLOCATE(4._q*ndata* NSIZE , "pot")
#endif
#ifdef double_prec_bse
          CALL DEREGISTER_ALLOCATE(8._q*ndata* NSIZE , "pot")
#endif

          CALL DEALLOCATE_BSE_MATRIX_INDEX( POT_INDEX )
       ENDIF
#endif
    ELSE
       NSIZE=SIZE(AMAT)
       DEALLOCATE(AMAT, AMAT_SCALA, R)

       IF (ANTIRES>=2) THEN
          NSIZE=NSIZE+SIZE(BMAT)
          DEALLOCATE(BMAT, BMAT_SCALA, RB)
       ENDIF
       IF ((.NOT. LQP) .AND. (NALPHA > 1 .OR. LALPHA)) THEN
          NSIZE=NSIZE+SIZE(VMAT)
          DEALLOCATE(VMAT_SCALA, VMAT)
       ENDIF
#ifdef single_prec_bse
       CALL DEREGISTER_ALLOCATE(4._q*ndata* NSIZE , "bse")
#endif
#ifdef double_prec_bse
       CALL DEREGISTER_ALLOCATE(8._q*ndata* NSIZE , "bse")
#endif

       IF (LQP) THEN
          NSIZE=SIZE(AVpW)+SIZE(VMAT)+SIZE(BVpW)
          DEALLOCATE(VMAT_SCALA, VMAT, AVpW_SCALA, AVpW, BVpW_SCALA, BVpW)
#ifdef single_prec_bse
          CALL DEREGISTER_ALLOCATE(4._q*ndata* NSIZE , "pot")
#endif
#ifdef double_prec_bse
          CALL DEREGISTER_ALLOCATE(8._q*ndata* NSIZE , "pot")
#endif
          CALL DEALLOCATE_BSE_MATRIX_INDEX( POT_INDEX )
       ENDIF
    ENDIF

    CALL DUMP_ALLOCATE(IO%IU6)

    RETURN

    CONTAINS

 
    ! ***********************************************************
    ! checks for allocation error and stops code with message
    ! ***********************************************************
    SUBROUTINE CHECK_ALLOCATION_ERROR( STRING, ISIZE, ISTAT ) 
       CHARACTER( LEN=* )  :: STRING
       INTEGER, INTENT(IN) :: ISIZE
       INTEGER, INTENT(IN) :: ISTAT

       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "CALCULATE_BSE ("//STRING//") is not able to allocate "//&
#ifdef single_prec_bse
            str( 4._q*ISIZE/1024) //" kB of data on MPI rank 0." )
#endif
#ifdef double_prec_bse
            str( 8._q*ISIZE/1024) //" kB of data on MPI rank 0." )
#endif
       ENDIF
       
    END SUBROUTINE CHECK_ALLOCATION_ERROR

    ! ***********************************************************
    ! same as above, but prints only double precision sizes
    ! ***********************************************************
    SUBROUTINE CHECK_ALLOCATION_ERROR_Q( STRING, ISIZE, ISTAT ) 
       CHARACTER( LEN=* )  :: STRING
       INTEGER, INTENT(IN) :: ISIZE
       INTEGER, INTENT(IN) :: ISTAT

       IF ( ISTAT/=0 ) THEN
          CALL VTUTOR%ERROR( "CALCULATE_BSE ("//STRING//") is not able to allocate "//&
            str( 8._q*ISIZE/1024) //" kB of data on MPI rank 0." )
       ENDIF
       
    END SUBROUTINE CHECK_ALLOCATION_ERROR_Q


  END SUBROUTINE CALCULATE_BSE



!**********************************************************************
!
! some preprocessor statements to allow for simpler coding
! essentially map e.g
! GGEMMbse ->   SGEMM   CGEMM      DGEMM   ZGEMM
! one      ->   1.      (1.0,0.0)  1._q    (1.0_q,0.0_q)
!
!**********************************************************************

#ifdef single_prec_bse

#ifdef gammareal
#define GGEMMbse SGEMM
#define GETRIbse SGETRI
#define GETRFbse SGETRF
#define onebse   1._qs
#define zerobse  0._qs
#else
#define GGEMMbse CGEMM
#define GETRIbse CGETRI
#define GETRFbse CGETRF
#define onebse   (1._qs,0._qs)
#define zerobse  (0._qs,0._qs)
#endif

#else

#ifdef gammareal
#define GGEMMbse DGEMM
#define GETRIbse DGETRI
#define GETRFbse DGETRF
#define onebse   1._q
#define zerobse  0._q
#else
#define GGEMMbse ZGEMM
#define GETRIbse ZGETRI
#define GETRFbse ZGETRF
#define onebse   (1._q,0._q)
#define zerobse  (0._q,0._q)
#endif

#endif

#ifdef scaLAPACK

! define some shorthands for scaLAPACK
#ifdef single_prec_bse

#ifdef gammareal
#define GPGEMMbse PSGEMM
#define GPTRANbse PSTRAN
#else
#define GPGEMMbse PCGEMM
#define GPTRANbse PCTRANU
#endif

#else

#ifdef gammareal
#define GPGEMMbse PDGEMM
#define GPTRANbse PDTRAN
#else
#define GPGEMMbse PZGEMM
#define GPTRANbse PZTRANU
#endif

#endif

#endif


!****************** SUBROUTINE  TWOELECTRON4O_STORE *******************
!
! this subroutine stores the two electron 4 orbital integrals
! in the BSE / TDFT Hamilton matrix AMAT
! this version of the routine is not scaLAPACK "aware" and
! requires that the full BSE matrix is stored on every node
! 
!**********************************************************************

  SUBROUTINE TWOELECTRON4O_STORE(WHF, ISP, ISP2, ISP_LOW, &
       K1, NPOS1, NSTRIP1, K2, NPOS2, NSTRIP2, K3, K4, &
       TWOELECTRON3O, B1, B2, AMAT, BSE_INDEX, IDIAG, ALPHAPT, SCALE_PERT)

    USE constant
    IMPLICIT NONE

    ! passed variables
    TYPE (wavespin) WHF
    INTEGER ISP, ISP2, ISP_LOW
    INTEGER K1, K2, K3, K4, NPOS1, NSTRIP1, NPOS2, NSTRIP2
    TYPE(banddesc) :: B1,B2         ! current descriptor for the spin loops
    GDEF ::  TWOELECTRON3O(:,:,:,:)
    REAL(q) :: ALPHAPT ! coupling constant value for AC integration
    REAL(q):: SCALE_PERT
! actual dim TWOELECTRON3O( NSTRIPV, NGLB, NSTRIPV, NGLB4))
#ifdef single_prec_bse
    GDEFS :: AMAT(:,:)
#endif

#ifdef double_prec_bse
    GDEF :: AMAT(:,:)
#endif
  ! local variables
    INTEGER N1, N2, N3, N4, NB1_INTO_TOT, NB2_INTO_TOT, NB3_INTO_TOT, NB4_INTO_TOT
    INTEGER NB3, NB4
    INTEGER NCV13, NCV24
    TYPE(bse_matrix_index) :: BSE_INDEX
    INTEGER :: IDIAG     ! add eigenvalues

    DO N1=1,NSTRIP1
       DO N2=1,NSTRIP2
          NB1_INTO_TOT=NPOS1-1+N1
          NB2_INTO_TOT=NPOS2-1+N2
          DO N4=1,B2%CBMAX4-B2%CBMIN4+1
             NB4_INTO_TOT=B2%CBMIN4-1+N4
             NB4         =N4
#ifdef check_index
             IF (N1 > SIZE(TWOELECTRON3O,1)) THEN
                CALL vtutor%bug("internal error in TWOELECTRON4O_STORE: out of bounds 1 " // str(N1), __FILE__, __LINE__)
             ENDIF

             IF (N2 > SIZE(TWOELECTRON3O,3)) THEN
                CALL vtutor%bug("internal error in TWOELECTRON4O_STORE: out of bounds 3 " // str(N2), __FILE__, __LINE__)
             ENDIF

             IF (NB4 > SIZE(TWOELECTRON3O,4)) THEN
                CALL vtutor%bug("internal error in TWOELECTRON4O_STORE: out of bounds 4 " // str(NB4), __FILE__, __LINE__)
             ENDIF
#endif
           ! NCV24=NB2_INTO_TOT-B1%VBMIN+1 +(B1%VBMAX-B1%VBMIN+1)*(NB4_INTO_TOT-B1%CBMIN+(B1%CBMAX-B1%CBMIN+1)*(K2-1))
             NCV24=BSE_INDEX%INDEX(NB2_INTO_TOT-B2%VBMIN+1, NB4_INTO_TOT-B2%CBMIN+1, K2, ISP2-ISP_LOW+1)
             IF (NCV24==0) CYCLE ! cycle if this pair is not included
#ifdef check_index
             IF ( NCV24 > SIZE(AMAT,2) ) THEN
                CALL vtutor%bug(str(NB2_INTO_TOT-B2%VBMIN+1) // " " // str(B2%VBMAX-B2%VBMIN+1) // &
                   " " // str(NB4_INTO_TOT-B2%CBMIN+1) // " " // str(B2%CBMAX-B2%CBMIN+1) // " " // &
                   str(K2-1) // "\n internal error in TWOELECTRON4O_STORE: NCV exceeds AMAT &
                   &size-boundaries " // str(NCV13) // " " // str(NCV24) // " " // str(SIZE(AMAT,1)) &
                   // " " // str(SIZE(AMAT,2)), __FILE__, __LINE__)
             ENDIF
#endif
             DO N3=1,B1%CBMAX-B1%CBMIN+1
                NB3_INTO_TOT=B1%CBMIN-1+N3
                NB3         =N3
#ifdef check_index
                IF ( NB3> SIZE(TWOELECTRON3O,2)) THEN
                   CALL vtutor%bug("internal error in TWOELECTRON4O_STORE: out of bounds 2 " // &
                      str(NB3_INTO_TOT) // " " // str(NB3), __FILE__, __LINE__)
                ENDIF
#endif
                ! build index NCV13 and NCV24
                ! NCV13=NB1_INTO_TOT-B1%VBMIN+1 +(B1%VBMAX-B1%VBMIN+1)*(NB3_INTO_TOT-B1%CBMIN+(B1%CBMAX-B1%CBMIN+1)*(K1-1))
                NCV13=BSE_INDEX%INDEX(NB1_INTO_TOT-B1%VBMIN+1, NB3_INTO_TOT-B1%CBMIN+1, K1, ISP -ISP_LOW+1)
                IF (NCV13==0) CYCLE ! cycle if this pair is not included
                ! check whether NCV13 is smaller or equal than SIZE(AMAT,1)
#ifdef check_index
                IF ( NCV13 > SIZE(AMAT,1) ) THEN
                   CALL vtutor%bug(str(NB2_INTO_TOT-B2%VBMIN+1) // " " // str(B2%VBMAX-B2%VBMIN+1) &
                      // " " // str(NB4_INTO_TOT-B2%CBMIN+1) // " " // str(B2%CBMAX-B2%CBMIN+1) // "&
                      & " // str(K2-1) // "\n internal error in TWOELECTRON4O_STORE: NCV exceeds AMAT &
                      &size-boundaries " // str(NCV13) // " " // str(NCV24) // " " // str(SIZE(AMAT,&
                      1)) // " " // str(SIZE(AMAT,2)), __FILE__, __LINE__)
                ENDIF
                IF (AMAT(NCV13,NCV24) /= 0) THEN
                   CALL vtutor%bug("internal error in TWOELECTRON4O_STORE: AMAT was already written &
                      &" // str(NCV13) // " " // str(NCV24) // "\n attempt of overwriting points &
                      &towards inconsistent indexing \n correct the formulas for NCV13 and NCV24 in bse.F", &
                      __FILE__, __LINE__)
                ENDIF
#endif
                AMAT(NCV13,NCV24)=TWOELECTRON3O( N1, NB3 , N2, NB4)*WHF%WDES%WTKPT(1)*SCALE_PERT*ALPHAPT

                ! diagonal components 
                ! subtract eigenvalue difference for the A matrix (IDIAG is 0.0 when BMAT is passed on)
                IF (K1==K2 .AND. NB1_INTO_TOT==NB2_INTO_TOT .AND. NB3_INTO_TOT==NB4_INTO_TOT .AND.  ISP== ISP2) THEN
                   IF (NCV13/= NCV24) THEN
                      CALL vtutor%bug("internal error in TWOELECTRON4O_STORE: not a diagonal element " &
                         // str(NCV13) // " " // str(NCV24) // "\n " // str(N1) // " " // str(N2) &
                         // " " // str(N3) // " " // str(N4) // " " // str(K1) // " " // str(K2), __FILE__, __LINE__)
                   ENDIF
                   IF(IDIAG>0) THEN
                      !! catch exception if the aufbau principle is violated (i.e.: E_k1,n1 > E_k3,n3)
                      IF((REAL(WHF%CELTOT(NB1_INTO_TOT,K1, ISP))>REAL(WHF%CELTOT(NB3_INTO_TOT,K3, ISP)))) THEN
                         WRITE(*,*) 'internal WARNING: Aufbau principle is violated, might need reexamination'
                      ELSE
                          AMAT(NCV13, NCV24) = AMAT(NCV13, NCV24)-IDIAG*(WHF%CELTOT(NB1_INTO_TOT,K1, ISP)-WHF%CELTOT(NB3_INTO_TOT,K3, ISP))
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE TWOELECTRON4O_STORE


!****************** SUBROUTINE  TWOELECTRON4O_STORE_POT ***************
!
! this subroutine stores the two electron 4 orbital integrals
! in the BSE / TDFT Hamilton matrix AMAT
! this version of the routine is not scaLAPACK "aware" and
! requires that the full BSE matrix is stored on every node
! 
!**********************************************************************

  SUBROUTINE TWOELECTRON4O_STORE_POT(WHF, ISP, ISP2, ISP_LOW, &
       K1, NPOS1, NSTRIP1, K2, NPOS2, NSTRIP2, K3, K4, &
       TWOELECTRON3O, B1, B2, AMAT, BSE_INDEX, POT_INDEX, ALPHAPT, LTRIANGLE)

    USE constant
    IMPLICIT NONE

    ! passed variables
    TYPE (wavespin) WHF
    INTEGER ISP, ISP2, ISP_LOW
    INTEGER K1, K2, K3, K4, NPOS1, NSTRIP1, NPOS2, NSTRIP2
    TYPE(banddesc) :: B1,B2         ! current descriptor for the spin loops
    GDEF ::  TWOELECTRON3O(:,:,:,:)
    REAL(q) :: ALPHAPT              ! coupling constant value for AC integration
    LOGICAL :: LTRIANGLE            ! true if upper and lower triangle are calculated anyway
! actual dim TWOELECTRON3O( NSTRIPV, NGLB, NSTRIPV, NGLB4))
#ifdef single_prec_bse
    GDEFS :: AMAT(:,:)
#endif

#ifdef double_prec_bse
    GDEF :: AMAT(:,:)
#endif
  ! local variables
    INTEGER N1, N2, N3, N4, NB1_INTO_TOT, NB2_INTO_TOT, NB3_INTO_TOT, NB4_INTO_TOT
    INTEGER NB3, NB4
    INTEGER NCV13, NCV24, NCV13P, NCV24P
    TYPE(bse_matrix_index) :: BSE_INDEX, POT_INDEX

    DO N1=1,NSTRIP1
       DO N2=1,NSTRIP2
          NB1_INTO_TOT=NPOS1-1+N1
          NB2_INTO_TOT=NPOS2-1+N2
          DO N4=1,B2%CBMAX4-B2%CBMIN4+1
             NB4_INTO_TOT=B2%CBMIN4-1+N4
             NB4         =N4
             NCV24 =BSE_INDEX%INDEX(NB2_INTO_TOT-B2%VBMIN+1, NB4_INTO_TOT-B2%CBMIN+1, K2, ISP2-ISP_LOW+1)
             NCV13P=POT_INDEX%INDEX(NB2_INTO_TOT-B2%VBMIN+1, NB4_INTO_TOT-B2%CBMIN+1, K2, ISP2-ISP_LOW+1)

             DO N3=1,B1%CBMAX-B1%CBMIN+1
                NB3_INTO_TOT=B1%CBMIN-1+N3
                NB3         =N3
                ! build index NCV13 and NCV24
                NCV13 =POT_INDEX%INDEX(NB1_INTO_TOT-B1%VBMIN+1, NB3_INTO_TOT-B1%CBMIN+1, K1, ISP -ISP_LOW+1)
                NCV24P=BSE_INDEX%INDEX(NB1_INTO_TOT-B1%VBMIN+1, NB3_INTO_TOT-B1%CBMIN+1, K1, ISP -ISP_LOW+1)

                IF (NCV13/=0 .AND. NCV24/=0) THEN
                   IF (AMAT(NCV13,NCV24) /= 0) THEN
                      CALL vtutor%bug("internal error in TWOELECTRON4O_STORE_POT: AMAT was already &
                         &written " // str(NCV13) // " " // str(NCV24) // "\n attempt of overwriting &
                         &points towards inconsistent indexing \n correct the formulas for NCV13 and " &
                         // "NCV24 in bse.F", __FILE__, __LINE__)
                   ENDIF
                   AMAT(NCV13,NCV24)=TWOELECTRON3O( N1, NB3 , N2, NB4)*WHF%WDES%WTKPT(1)*ALPHAPT
                ENDIF
                ! add conjugated elements: causes issues on many cores
                ! TODO: it is unclear why this was ever present
                ! since a similar thing misses in TWOELECTRON4O_STORE
                IF (NCV13P/=0 .AND. NCV24P/=0 .AND. LTRIANGLE) THEN
                   AMAT(NCV13P,NCV24P)=GCONJG(TWOELECTRON3O( N1, NB3 , N2, NB4))*WHF%WDES%WTKPT(1)*ALPHAPT
                ENDIF                
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE TWOELECTRON4O_STORE_POT

#ifdef scaLAPACK
!****************** SUBROUTINE  TWOELECTRON4O_STORE_SCALA  ************
!
! this subroutine stores the two electron 4 orbital integrals
! in the BSE / TDFT Hamilton matrix AMAT
! this version is scaLAPACK "aware" and matter of fact only
! available if scalapack is used 
! 
!**********************************************************************

  SUBROUTINE TWOELECTRON4O_STORE_SCALA(WHF, ISP, ISP2, ISP_LOW, &
       K1, NPOS1, NSTRIP1, K2_, NPOS2, NSTRIP2, K3, K4_, &
       TWOELECTRON3O, B1, B2, AMAT, NCV, BSE_INDEX, IDIAG, ALPHAPT, SCALE_PERT, DESC)

    USE constant
    IMPLICIT NONE

    ! passed variables
    TYPE (wavespin) WHF
    INTEGER ISP, ISP2, ISP_LOW
    INTEGER K1, K2_, K3, K4_, NPOS1, NSTRIP1, NPOS2, NSTRIP2
    TYPE(banddesc) :: B1,B2         ! current descriptor for the spin loops
    GDEF ::  TWOELECTRON3O(:,:,:,:)
! actual dim TWOELECTRON3O( NSTRIPV, NGLB, NSTRIPV, NGLB4))
    TYPE(bse_matrix_index) :: BSE_INDEX
    INTEGER :: IDIAG     ! add eigenvalues
    REAL(q) :: ALPHAPT !value of the coupling constant for AC integration 
    REAL(q) :: SCALE_PERT 
    INTEGER :: DESC(DLEN_)

#ifdef single_prec_bse
    GDEFS :: AMAT(:)
#endif
#ifdef double_prec_bse
    GDEF :: AMAT(:)
#endif
    INTEGER :: NCV
  ! local variables
    GDEF, ALLOCATABLE ::  TWOELECTRON3O_LOCAL(:,:,:,:)
    INTEGER N1, N2, NB1_INTO_TOT, NB2_INTO_TOT, NB3_INTO_TOT, NB4_INTO_TOT
    INTEGER NB3, NB4
    INTEGER NCV13, NCV24
    INTEGER K2, K4, CBMIN4, CBMAX4
    INTEGER N
    INTEGER ITWOELECTRON3O(4)
! BLACS variables
    INTEGER, EXTERNAL ::     NUMROC
    INTEGER MYROW, MYCOL, NPROW, NPCOL, NP, NQ
    INTEGER I1RES, J1RES, IROW, JCOL
    INTEGER I1, I2, J1, J2
    INTEGER IFAIL

! allocate work array into which data is subsequently merged from all nodes
    ITWOELECTRON3O=SHAPE(TWOELECTRON3O)
    CALL M_max_i(WHF%WDES%COMM_INTER, ITWOELECTRON3O, 4 )

    ALLOCATE(TWOELECTRON3O_LOCAL(ITWOELECTRON3O(1),ITWOELECTRON3O(2), & 
                                 ITWOELECTRON3O(3),ITWOELECTRON3O(4)))

!==========================================================================
! loop over all nodes and broadcast their TWOELECTRO3O array
! to other nodes also broadcast the storage position (B2%CBMIN4, B2%CBMAX4)
! and the local k-point indices
!==========================================================================
cp: DO N=1,WHF%WDES%COMM_INTER%NCPU

! broadcast K2, K4 and all the other integers to all nodes
    IF (N==WHF%WDES%COMM_INTER%NODE_ME) THEN
       K2=K2_
       K4=K4_
       CBMIN4=B2%CBMIN4
       CBMAX4=B2%CBMAX4
    ENDIF
    CALL M_bcast_i_from(WHF%WDES%COMM_INTER, K2, 1, n)
    CALL M_bcast_i_from(WHF%WDES%COMM_INTER, K4, 1, n)
    CALL M_bcast_i_from(WHF%WDES%COMM_INTER, CBMIN4, 1, n)
    CALL M_bcast_i_from(WHF%WDES%COMM_INTER, CBMAX4, 1, n)

! if K2 is properly set distribute the matrix elements
    IF (K2>=1) THEN

    TWOELECTRON3O_LOCAL=0
    IF (N==WHF%WDES%COMM_INTER%NODE_ME) THEN
       TWOELECTRON3O_LOCAL(1:SIZE(TWOELECTRON3O,1),1:SIZE(TWOELECTRON3O,2),1:SIZE(TWOELECTRON3O,3),1:SIZE(TWOELECTRON3O,4))=TWOELECTRON3O
    ENDIF
#ifdef gammareal 
    CALL M_bcast_d_from(WHF%WDES%COMM_INTER, TWOELECTRON3O_LOCAL, SIZE(TWOELECTRON3O_LOCAL), N)
#else
    CALL M_bcast_z_from(WHF%WDES%COMM_INTER, TWOELECTRON3O_LOCAL, SIZE(TWOELECTRON3O_LOCAL), N)
#endif
!==========================================================================
! now distribute the matrix elements among nodes
!==========================================================================
    CALL BLACS_GRIDINFO(DESC(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

    NP = NUMROC(NCV,DESC(MB_),MYROW,0,NPROW)
    NQ = NUMROC(NCV,DESC(NB_),MYCOL,0,NPCOL)

    ! loop over all columns of the global matrix
    JCOL=0
    DO J1=1,NQ,DESC(NB_)
       J1RES=MIN(DESC(NB_),NQ-J1+1)
       DO J2=1,J1RES
          JCOL=JCOL+1
        ! global column index NCV24
          NCV24=DESC(NB_)*MYCOL+NPCOL*(J1-1)+J2
        ! spin test
          IF (ISP2 /= BSE_INDEX%ISP(NCV24)) CYCLE
        ! k-point test, CYCLE if not identical
        ! IF (K2 /= 1+(NCV24-1)/((B1%VBMAX-B1%VBMIN+1)*(B1%CBMAX-B1%CBMIN+1))) CYCLE
          IF (K2 /= BSE_INDEX%NK(NCV24)) CYCLE

        ! test on NB2
        ! NB2_INTO_TOT=MOD(NCV24-1,B1%VBMAX-B1%VBMIN+1)+B1%VBMIN
          NB2_INTO_TOT=BSE_INDEX%N1(NCV24)
          N2=NB2_INTO_TOT-NPOS2+1
          IF (N2<1 .OR. N2 > NSTRIP2) CYCLE

        ! test on NB4_INTO_TOT
        ! NB4_INTO_TOT=MOD((NCV24-1)/(B1%VBMAX-B1%VBMIN+1),(B1%CBMAX-B1%CBMIN+1))+B1%CBMIN
          NB4_INTO_TOT=BSE_INDEX%N3(NCV24)
          IF (NB4_INTO_TOT<CBMIN4 .OR. NB4_INTO_TOT >CBMAX4) CYCLE

          NB4=NB4_INTO_TOT-CBMIN4+1

#ifdef check_index
          ! internal consistency test (could be removed)
        ! IF (NCV24/=NB2_INTO_TOT-B2%VBMIN+1 +(B2%VBMAX-B2%VBMIN+1)*(NB4_INTO_TOT-B2%CBMIN+(B2%CBMAX-B2%CBMIN+1)*(K2-1))) THEN
          IF (NCV24/=BSE_INDEX%INDEX(NB2_INTO_TOT-B2%VBMIN+1, NB4_INTO_TOT-B2%CBMIN+1, K2, ISP2-ISP_LOW+1)) THEN
             CALL vtutor%bug("internal error 1 in TWOELECTRON4O_STORE_SCALA: " // str(NCV24) // " &
                &" // str(NB2_INTO_TOT-B2%VBMIN+1) // " " // str(NB4_INTO_TOT-B2%CBMIN) // " " // &
                str(K2) // " " // str(ISP2) // " " // str(BSE_INDEX%INDEX(NB2_INTO_TOT-B2%VBMIN+1,&
                NB4_INTO_TOT-B2%CBMIN+1,K2,ISP2-ISP_LOW+1)) // " " // str(BSE_INDEX%ISP(NCV24)), __FILE__, __LINE__)
          ENDIF
#endif
          IROW=0
          DO I1=1,NP,DESC(MB_)
             I1RES=MIN(DESC(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW=IROW+1
              ! global row index NCV13
                NCV13=DESC(MB_)*MYROW+NPROW*(I1-1)+I2
              ! spin test
                IF (ISP /= BSE_INDEX%ISP(NCV13)) CYCLE
              ! k-point test
              ! IF (K1 /= 1+(NCV13-1)/((B1%VBMAX-B1%VBMIN+1)*(B1%CBMAX-B1%CBMIN+1))) CYCLE
                IF (K1 /= BSE_INDEX%NK(NCV13)) CYCLE

              ! test on NB1
              ! NB1_INTO_TOT=MOD(NCV13-1,B1%VBMAX-B1%VBMIN+1)+B1%VBMIN
                NB1_INTO_TOT=BSE_INDEX%N1(NCV13)
                N1=NB1_INTO_TOT-NPOS1+1
                IF (N1<1 .OR. N1 > NSTRIP1) CYCLE

              ! test on NB3
              ! NB3_INTO_TOT=MOD((NCV13-1)/(B1%VBMAX-B1%VBMIN+1),(B1%CBMAX-B1%CBMIN+1))+B1%CBMIN
                NB3_INTO_TOT=BSE_INDEX%N3(NCV13)
                NB3=NB3_INTO_TOT-B1%CBMIN+1

#ifdef check_index
                ! internal consistency test (could be removed)
              ! IF (NCV13/=NB1_INTO_TOT-B1%VBMIN+1 +(B1%VBMAX-B1%VBMIN+1)*(NB3_INTO_TOT-B1%CBMIN+(B1%CBMAX-B1%CBMIN+1)*(K1-1))) THEN
                IF (NCV13/=BSE_INDEX%INDEX(NB1_INTO_TOT-B1%VBMIN+1, NB3_INTO_TOT-B1%CBMIN+1, K1, ISP -ISP_LOW+1)) THEN
                   CALL vtutor%bug("internal error 2 in TWOELECTRON4O_STORE_SCALA: " // str(NCV13) &
                      // " " // str(NB1_INTO_TOT-B1%VBMIN+1) // " " // str(NB3_INTO_TOT) // " " // &
                      str(K1) // " " // str(ISP) // " " // str(BSE_INDEX%INDEX(NB1_INTO_TOT-B1%VBMIN+1,&
                      NB3_INTO_TOT-B1%CBMIN+1,K1,ISP -ISP_LOW+1)) // " " // str(BSE_INDEX%ISP(NCV13)), __FILE__, __LINE__)
                ENDIF

                IF ( N1>SIZE(TWOELECTRON3O_LOCAL,1) .OR. NB3>SIZE(TWOELECTRON3O_LOCAL,2) & 
                 .OR.N2>SIZE(TWOELECTRON3O_LOCAL,3) .OR. NB4>SIZE(TWOELECTRON3O_LOCAL,4)) THEN
                   CALL vtutor%bug("internal error 3 in TWOELECTRON4O_STORE_SCALA: " // str(N1) // &
                      " " // str(SIZE(TWOELECTRON3O_LOCAL,1)) // " " // str(NB3) // " " // &
                      str(SIZE(TWOELECTRON3O_LOCAL,2)) // " " // str(N2) // " " // &
                      str(SIZE(TWOELECTRON3O_LOCAL,3)) // " " // str(NB4) // " " // &
                      str(SIZE(TWOELECTRON3O_LOCAL,4)) // " " // str(CBMAX4), __FILE__, __LINE__)
                ENDIF
                IF ( IROW+(JCOL-1)*DESC(LLD_)>SIZE(AMAT) ) THEN
                   CALL vtutor%bug("internal error 4 in TWOELECTRON4O_STORE_SCALA: " // &
                      str(IROW+(JCOL-1))(LLD_) // " " // str(SIZE(AMAT)) // " " // str(IROW) // " " &
                      // str(JCOL) // " " // str(DESC(LLD_)), __FILE__, __LINE__)
                ENDIF
#endif
                ! finally if all test have a green light store the element in AMAT
                AMAT(IROW+(JCOL-1)*DESC(LLD_))= &
                     TWOELECTRON3O_LOCAL( N1, NB3 , N2, NB4)*WHF%WDES%WTKPT(1)
                AMAT(IROW+(JCOL-1)*DESC(LLD_))= AMAT(IROW+(JCOL-1)*DESC(LLD_))*SCALE_PERT*ALPHAPT

                ! add diagonal element
                IF (NCV13==NCV24 .AND. IDIAG>0) THEN

                  !! catch exception if the aufbau principle is violated (ie: E_k1,n1 > E_k3,n3)
                  IF((REAL(WHF%CELTOT(NB1_INTO_TOT,K1, ISP)) > REAL(WHF%CELTOT(NB3_INTO_TOT,K3, ISP)))) THEN
                      WRITE(*,*) 'internal WARNING: aufbau principle is violated'
                      !AMAT(IROW+(JCOL-1)*DESCSTD(LLD_)) = AMAT(IROW+(JCOL-1)*DESCSTD(LLD_))
                  ELSE
                      AMAT(IROW+(JCOL-1)*DESC(LLD_)) = AMAT(IROW+(JCOL-1)*DESC(LLD_)) & 
                             -IDIAG*(WHF%CELTOT(NB1_INTO_TOT,K1, ISP)-WHF%CELTOT(NB3_INTO_TOT,K3, ISP))
                  ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !==========================================================================
    ! now fill in conjugated elements in lower triangle of matrix
    !==========================================================================

#define lower_triangle
#ifdef lower_triangle
    JCOL=0
    DO J1=1,NQ,DESC(NB_)
       J1RES=MIN(DESC(NB_),NQ-J1+1)
       DO J2=1,J1RES
          JCOL=JCOL+1
          ! global column index NCV24 species now the row in the conjugated matrix
          NCV24=DESC(NB_)*MYCOL+NPCOL*(J1-1)+J2

          ! spin test
          IF (ISP /= BSE_INDEX%ISP(NCV24)) CYCLE
          ! k-point test
          IF (K1 /= BSE_INDEX%NK(NCV24)) CYCLE
          
          ! test on NB1
          NB1_INTO_TOT=BSE_INDEX%N1(NCV24)
          N1=NB1_INTO_TOT-NPOS1+1
          IF (N1<1 .OR. N1 > NSTRIP1) CYCLE
          
          ! test on NB3
          NB3_INTO_TOT=BSE_INDEX%N3(NCV24)

          NB3=NB3_INTO_TOT-B1%CBMIN+1

          IROW=0
          DO I1=1,NP,DESC(MB_)
             I1RES=MIN(DESC(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW=IROW+1
              ! global row index NCV13 species the column position in the conjugated matrix
                NCV13=DESC(MB_)*MYROW+NPROW*(I1-1)+I2

                ! spin test
                IF (ISP2 /= BSE_INDEX%ISP(NCV13)) CYCLE
                ! k-point test, CYCLE if not identical
                IF (K2 /= BSE_INDEX%NK(NCV13)) CYCLE

                ! test on NB2
                NB2_INTO_TOT=BSE_INDEX%N1(NCV13)
                N2=NB2_INTO_TOT-NPOS2+1
                IF (N2<1 .OR. N2 > NSTRIP2) CYCLE
                
                ! test on NB4_INTO_TOT
                NB4_INTO_TOT=BSE_INDEX%N3(NCV13)
                IF (NB4_INTO_TOT<CBMIN4 .OR. NB4_INTO_TOT >CBMAX4) CYCLE

                NB4=NB4_INTO_TOT-CBMIN4+1

                IF (NCV13 > NCV24) THEN
                   AMAT(IROW+(JCOL-1)*DESC(LLD_))= &
                        GCONJG(TWOELECTRON3O_LOCAL(N1, NB3, N2, NB4 ))*WHF%WDES%WTKPT(1)*SCALE_PERT*ALPHAPT
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
#endif

    ENDIF
    ENDDO cp

    DEALLOCATE(TWOELECTRON3O_LOCAL)
  END SUBROUTINE TWOELECTRON4O_STORE_SCALA
#endif

!****************** TWOELECTRON4O_GET_DIAG ****************************
!
! this subroutine returns the sum of the diagonal elements of matrix A
! when this is not stored in a scaLAPACK aware way
! 
!**********************************************************************

  SUBROUTINE TWOELECTRON4O_GET_DIAG(AMAT, NCV, CSUM )

    USE constant
    IMPLICIT NONE
    INTEGER NCV

#ifdef single_prec_bse
    GDEFS :: AMAT(:,:)
#endif

#ifdef double_prec_bse
    GDEF :: AMAT(:,:)
#endif
    COMPLEX(q) :: CSUM
    INTEGER I,J

    CSUM=0.0_q
    DO I=1,NCV
       CSUM=CSUM+AMAT(I,I)
    ENDDO


  END SUBROUTINE TWOELECTRON4O_GET_DIAG

#ifdef scaLAPACK
!****************** TWOELECTRON4O_GET_DIAG_SCALA  *********************
!
! this subroutine determines the sum of the diagonal elements 
! of a matrix A
! 
!**********************************************************************

  SUBROUTINE TWOELECTRON4O_GET_DIAG_SCALA( AMAT, NCV, CSUM, MY_COMM, DESC, CTXT )

    USE constant
    IMPLICIT NONE

#ifdef single_prec_bse
    GDEFS :: AMAT(:)          ! input matrix 
#endif
#ifdef double_prec_bse
    GDEF :: AMAT(:)           ! input matrix 
#endif
    INTEGER :: NCV            ! dimension of global matrix
    COMPLEX(q) :: CSUM        ! sum of diagonal components
    TYPE(communic) :: MY_COMM ! communicator to globally sum the diagonal components
    INTEGER :: DESC(DLEN_)
    INTEGER, OPTIONAL :: CTXT

  ! local variables
    INTEGER NCV13, NCV24
  ! BLACS variables
    INTEGER, EXTERNAL ::     NUMROC
    INTEGER MYROW, MYCOL, NPROW, NPCOL, NP, NQ
    INTEGER IROW, JCOL
    INTEGER I1, I2, J1, J2, J1RES, I1RES

    CALL BLACS_GRIDINFO(DESC(CTXT_), NPROW, NPCOL, MYROW, MYCOL)

    NP = NUMROC(NCV,DESC(MB_),MYROW,0,NPROW)
    NQ = NUMROC(NCV,DESC(NB_),MYCOL,0,NPCOL)
    ! loop over all columns of the global matrix
    JCOL=0
    ! loop over block along column index
    DO J1=1,NQ,DESC(NB_)
       J1RES=MIN(DESC(NB_),NQ-J1+1)
       ! loop over elements in this block
       DO J2=1,J1RES
          JCOL=JCOL+1

          ! original global column index (if matrix where stored sequentially)
          NCV24=DESC(NB_)*MYCOL+NPCOL*(J1-1)+J2

           IROW=0
         ! loop over blocks along row index
          DO I1=1,NP,DESC(MB_)
             I1RES=MIN(DESC(MB_),NP-I1+1)
            ! loop over elements in this block
             DO I2=1,I1RES
                IROW=IROW+1
                ! original global row index (if matrix where stored sequentially)
                NCV13=DESC(MB_)*MYROW+NPROW*(I1-1)+I2

                IF (NCV13 == NCV24) THEN
                   CSUM= CSUM + AMAT(IROW+(JCOL-1)*DESC(LLD_))
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    CALLMPI( M_sum_z(MY_COMM, CSUM, 1))

  END SUBROUTINE TWOELECTRON4O_GET_DIAG_SCALA
#endif


#ifdef use_bse_te
!****************** SUBROUTINE TW4O_STORE_STRIP ***********************
!
! store matrix elements as stripes on target nodes
! certainly not the most clever routine since data are sent
! everywhere and than the target nodes pick their required data
!
!**********************************************************************

  SUBROUTINE TW4O_STORE_STRIP(WHF, ISP, ISP2, ISP_LOW, &
       K1, NPOS1, NSTRIP1, K2_, NPOS2, NSTRIP2, K3, K4_, &
       TWOELECTRON3O, B1, B2, BSE_INDEX )

    USE constant
    IMPLICIT NONE

    ! passed variables
    TYPE (wavespin) WHF
    INTEGER ISP, ISP2, ISP_LOW
    INTEGER K1, K2_, K3, K4_, NPOS1, NSTRIP1, NPOS2, NSTRIP2
    GDEF ::  TWOELECTRON3O(:,:,:,:)
! actual dim TWOELECTRON3O( NSTRIPV, NGLB, NSTRIPV, NGLB4))
    TYPE(banddesc) :: B1,B2         ! current descriptor for the spin loops
    TYPE(bse_matrix_index) :: BSE_INDEX
! local variables
    INTEGER N1, N2, N3, N4, NB1_INTO_TOT, NB2_INTO_TOT, NB3_INTO_TOT, NB4_INTO_TOT
    INTEGER NB3, NB4
    INTEGER NCV13, NCV24
    INTEGER NCV13_LOCAL, NCV24_LOCAL
    INTEGER K2, K4, CBMIN4, CBMAX4
    INTEGER ITWOELECTRON3O(4)
    INTEGER N

    GDEF, ALLOCATABLE ::  TWOELECTRON3O_LOCAL(:,:,:,:)


! allocate work array into which data is subsequently merged from all nodes
    ITWOELECTRON3O=SHAPE(TWOELECTRON3O)
    CALL M_max_i(WHF%WDES%COMM_INTER, ITWOELECTRON3O, 4 )

    ALLOCATE(TWOELECTRON3O_LOCAL(ITWOELECTRON3O(1),ITWOELECTRON3O(2), & 
                                 ITWOELECTRON3O(3),ITWOELECTRON3O(4)))
!==========================================================================
! loop over all nodes and broadcast their TWOELECTRO3O array
! to other nodes also broadcast the storage position (B2%CBMIN4, B2%CBMAX4)
! and the local k-point indices
!==========================================================================
cp: DO N=1,WHF%WDES%COMM_INTER%NCPU

! broadcast K2, K4 and all the other integers to all nodes
    IF (N==WHF%WDES%COMM_INTER%NODE_ME) THEN
       K2=K2_
       K4=K4_
       CBMIN4=B2%CBMIN4
       CBMAX4=B2%CBMAX4
    ENDIF
    CALL M_bcast_i_from(WHF%WDES%COMM_INTER, K2, 1, n)
    CALL M_bcast_i_from(WHF%WDES%COMM_INTER, K4, 1, n)
    CALL M_bcast_i_from(WHF%WDES%COMM_INTER, CBMIN4, 1, n)
    CALL M_bcast_i_from(WHF%WDES%COMM_INTER, CBMAX4, 1, n)

! if K2 is properly set distribute the matrix elements
    IF (K2>=1) THEN

    TWOELECTRON3O_LOCAL=0
    IF (N==WHF%WDES%COMM_INTER%NODE_ME) THEN
       TWOELECTRON3O_LOCAL(1:SIZE(TWOELECTRON3O,1),1:SIZE(TWOELECTRON3O,2),1:SIZE(TWOELECTRON3O,3),1:SIZE(TWOELECTRON3O,4))=TWOELECTRON3O
    ENDIF
#ifdef gammareal 
    CALL M_bcast_d_from(WHF%WDES%COMM_INTER, TWOELECTRON3O_LOCAL, SIZE(TWOELECTRON3O_LOCAL), n)
#else
    CALL M_bcast_z_from(WHF%WDES%COMM_INTER, TWOELECTRON3O_LOCAL, SIZE(TWOELECTRON3O_LOCAL), n)
#endif
!==========================================================================
! now store the matrix elements on the local nodes
! simply walk through all received data and decide whether they
! are stored locally (stupid but fairly simple)
!==========================================================================
    DO N1=1,NSTRIP1
       DO N2=1,NSTRIP2
          NB1_INTO_TOT=NPOS1-1+N1
          NB2_INTO_TOT=NPOS2-1+N2
          DO N4=1,CBMAX4-CBMIN4+1
             NB4_INTO_TOT=CBMIN4-1+N4
             NB4         =N4
#ifdef check_index
             IF (NB4 > SIZE(TWOELECTRON3O_LOCAL,4)) THEN
                CALL vtutor%bug("internal error in TW4O_STORE_STRIP: out of bounds 4 " // str(NB4), __FILE__, __LINE__)
             ENDIF
#endif
             NCV24=BSE_INDEX%INDEX(NB2_INTO_TOT-B2%VBMIN+1, NB4_INTO_TOT-B2%CBMIN+1, K2, ISP2-ISP_LOW+1)
             IF (NCV24==0) CYCLE ! cycle if this pair is not included
             NCV24_LOCAL=NCV24-FIRST_ROW_INDEX+1
#ifdef check_index
                ! check whether NCV24 is smaller or equal than SIZE(BSEMATRIX,2)
                IF ( NCV24 > SIZE(BSEMATRIX,2)) THEN
                   CALL vtutor%bug("internal error in TW4O_STORE_STRIP: NCV exceeds BSEMATRIX &
                      &size-boundaries " // str(NCV13) // " " // str(NCV24) // " " // &
                      str(SIZE(BSEMATRIX,1)) // " " // str(SIZE(BSEMATRIX,2)), __FILE__, __LINE__)
                ENDIF
#endif
             DO N3=1,B1%CBMAX-B1%CBMIN+1
                NB3_INTO_TOT=B1%CBMIN-1+N3
                NB3         =N3
#ifdef check_index
                IF ( NB3> SIZE(TWOELECTRON3O_LOCAL,2)) THEN
                   CALL vtutor%bug("internal error in TW4O_STORE_STRIP: out of bounds 2 " // &
                      str(NB3_INTO_TOT) // " " // str(NB3), __FILE__, __LINE__)
                ENDIF
#endif
                ! build index NCV13 and NCV24
                NCV13=BSE_INDEX%INDEX(NB1_INTO_TOT-B1%VBMIN+1, NB3_INTO_TOT-B1%CBMIN+1, K1, ISP -ISP_LOW+1)
                IF (NCV13==0) CYCLE ! cycle if this pair is not included
                NCV13_LOCAL=NCV13-FIRST_ROW_INDEX+1
#ifdef check_index
                ! check whether NCV24 or  NCV13 is smaller or equal than SIZE(BSEMATRIX,2)
                IF ( NCV13 > SIZE(BSEMATRIX,2) ) THEN
                   CALL vtutor%bug("internal error in TW4O_STORE_STRIP: NCV exceeds BSEMATRIX &
                      &size-boundaries " // str(NCV13) // " " // str(NCV24) // " " // &
                      str(SIZE(BSEMATRIX,1)) // " " // str(SIZE(BSEMATRIX,2)), __FILE__, __LINE__)
                ENDIF
#endif
                IF ( NCV24_LOCAL>0 .AND. NCV24_LOCAL<=SIZE(BSEMATRIX,1)) THEN
                   BSEMATRIX(NCV24_LOCAL,NCV13)=TWOELECTRON3O_LOCAL( N1, NB3 , N2, NB4)*WHF%WDES%WTKPT(1)
                   ! diagonal components subtract eigenvalue difference
                   IF (NCV13== NCV24) THEN
                      BSEMATRIX( NCV24_LOCAL, NCV13 ) = BSEMATRIX( NCV24_LOCAL, NCV13 )-(WHF%CELTOT(NB1_INTO_TOT,K1, ISP)-WHF%CELTOT(NB3_INTO_TOT,K3, ISP))
                   ENDIF
                ENDIF
                IF ( NCV13_LOCAL>0 .AND. NCV13_LOCAL<=SIZE(BSEMATRIX,1)) THEN
                   BSEMATRIX(NCV13_LOCAL,NCV24)=CONJG(TWOELECTRON3O_LOCAL( N1, NB3 , N2, NB4)*WHF%WDES%WTKPT(1))
                   ! diagonal components subtract eigenvalue difference
                   IF (NCV13== NCV24) THEN
                      BSEMATRIX( NCV13_LOCAL, NCV24 ) = BSEMATRIX( NCV13_LOCAL, NCV24 )-(WHF%CELTOT(NB1_INTO_TOT,K1, ISP)-WHF%CELTOT(NB3_INTO_TOT,K3, ISP))
! store matrix in transposed form, this allows to change the DGEMM calls in bse_te.F GGEMM('N', -> GGEMM('T',
! also upon allocation of matrix the dimensions must be exchanged!
!                      BSEMATRIX( NCV24, NCV13_LOCAL ) = BSEMATRIX( NCV24, NCV13_LOCAL )-(WHF%CELTOT(NB1_INTO_TOT,K1, ISP)-WHF%CELTOT(NB3_INTO_TOT,K3, ISP))
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ENDIF
  ENDDO cp
  DEALLOCATE(TWOELECTRON3O_LOCAL)

  END SUBROUTINE TW4O_STORE_STRIP
#endif

!**********************************************************************
!
! calculate the BSE oscillator strength
! AMAT stores the eigenvectors of the BSE Hamiltonian
! the BSE Hamiltonian is calculated as
! Hamiltonian(I, J) = < c_k3,n3 v'_k2,n2 | v+W | v_k1,n1  c'_k4,n4 >
! with I = (k1,n1, k3, n3)   J= (k2,n2, k4, n4)
!
! after diagonalization AMAT stores U(I, lambda) with
! Hamiltonian(I, J)= U Lambda U+ =
!        \sum_lambda U(I, lambda) lambda U(J, lambda)*
! (this has been carefully checked by dumping matrix elements)
!
! the oscillator strength is given by
! \sum_lambda \sum_I <v_k1,n1|O|c_k3,n3> U(I, lambda)  1/ (w- lambda )
!             [sum_J U( J, lambda) <v_k2,n2|O|c_k4,n4>]*
! where O is an arbitrary operator
! To obtain the polarizabilty at a wavevector q=k3-k1 
!  < c_k3,n3| e^i(k3-k1 ) | v_k1,n1 >
!
!**********************************************************************

  SUBROUTINE CALCULATE_BSE_OSZI_STRENGTH (WHF, NQPOINT, KPOINT_BSE, LATT_CUR, ISP_LOW, ISP_HIGH,  &
       SHIFT, BD, BSE_INDEX, R, AMAT, ANTIRES, NOMEGA_DIM, OMEGAMAX, IU6)

    USE constant
    USE full_kpoints
    USE kpoints_change
    USE pseudo
    USE mymath
    IMPLICIT NONE

    ! passed variables
    TYPE (wavespin) WHF
    INTEGER :: NQPOINT       ! q-point considered
    INTEGER :: KPOINT_BSE(3) ! G-point considered
    TYPE (latt) LATT_CUR
    INTEGER ISP_LOW, ISP_HIGH
    REAL(q) SHIFT   ! complex shift for broadening
    TYPE(banddesc) :: BD(:)
    TYPE(bse_matrix_index) :: BSE_INDEX
    REAL(q) :: R(:)
#ifdef single_prec_bse
    GDEFS :: AMAT(:,:)
#else
    GDEF  :: AMAT(:,:)
#endif
    INTEGER :: ANTIRES
    INTEGER :: NOMEGA_DIM
    REAL(q) :: OMEGAMAX
    INTEGER :: IU6
    ! local 
    INTEGER :: NCV13
    INTEGER K1, K3
    INTEGER N1, N3, ISP
    INTEGER NB1_INTO_TOT, NB3_INTO_TOT
    COMPLEX(q) :: WEIGHT
    INTEGER I, J, NK1_IN_KPOINTS_FULL_ORIG, INFO
    INTEGER     :: NCV
    INTEGER     :: LAMBDA, NOMEGA
    GDEF        :: CDER_BETWEEN_STATE(3)
    COMPLEX(q)  :: CTRANS(3)
    COMPLEX(q)  :: CTRANS_SQUARE(3,3)
    COMPLEX(q)  :: CTRANS_SQUARE_ANTIRES(3,3)
    REAL(q)     :: DKX, DKY, DKZ
    COMPLEX(q)  :: EPS(3,3,NOMEGA_DIM)
    REAL(q)     :: OMEGA(NOMEGA_DIM),DOMEGA
    REAL(q)     :: EPS_REAL(NOMEGA_DIM, 3,3), EPS_IMAG(NOMEGA_DIM, 3, 3)
    REAL(q)     :: R_AND_INTENSITY(2,SIZE(R))
    REAL(q)     :: RSPIN

!  three cases need to be considered
!  non-collinear RSPIN=1.0
    IF (WHF%WDES%LNONCOLLINEAR) THEN
       RSPIN=1.0_q
    ELSE
!  collinear non spin polarized: RSPIN=2.0
!  collinear spin polarized:     RSPIN=1.0 for  ISP_LOW=1 and ISP_HIGH=2
!  I believe the case ISP_HIGH=ISP_LOW for collinear spin polarized makes little sense
!  but it is not used right now anyway
       RSPIN=2._q/(ISP_HIGH-ISP_LOW+1)
    ENDIF

    NCV=SIZE(AMAT,1)

    DOMEGA=R(NCV)/(NOMEGA_DIM-1)
    IF (OMEGAMAX>0) THEN
       DOMEGA=OMEGAMAX/(NOMEGA_DIM-1)
    ENDIF

    EPS = 0
    lambda_bse: DO LAMBDA=1,NCV
       CTRANS=0

       DO ISP=ISP_LOW,ISP_HIGH
       DO K1=1,WHF%WDES%NKPTS
          IF (WHF%WDES%WTKPT(K1)==0) CYCLE
          K3 =KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K1)+WHF%WDES%VKPT(:,NQPOINT),KPOINTS_FULL)

          ! determined the index of this k-point in the original full k-point grid
          NK1_IN_KPOINTS_FULL_ORIG=KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K1),KPOINTS_FULL_ORIG)

          DO N3=1,BD(ISP)%CBMAX-BD(ISP)%CBMIN+1
             NB3_INTO_TOT=N3+BD(ISP)%CBMIN-1
             IF (FILLED_XI_ORBITAL(WHF%FERTOT(NB3_INTO_TOT,K3,ISP))) CYCLE
             DO N1=1,BD(ISP)%VBMAX-BD(ISP)%VBMIN+1
                NB1_INTO_TOT=N1+BD(ISP)%VBMIN-1
                IF (EMPTY_XI_ORBITAL(WHF%FERTOT(NB1_INTO_TOT,K1,ISP))) CYCLE

                WEIGHT=(WHF%FERTOT(NB1_INTO_TOT,K1,ISP)-WHF%FERTOT(NB3_INTO_TOT, K3, ISP))

                ! CDER_BETWEEN_STATE =  <v_k1,n1|- i d/dq_j c_k1+q,n3> =  - <v_k1,n1| r | c_k1,n3>  for q->0
                ! or   <v_k1,n1| e-iqr  |c_k1+q,n3>  at finite q
                ! right now we have to used  <c_k1+q,n3| e iqr  |v_k1,n1> which is very odd
                ! see ONEELECTRON_TRANS in local_field.F
                CALL  CDER_BETWEEN_STATES_ROTATED(CDER_BETWEEN_STATE,LATT_CUR, NK1_IN_KPOINTS_FULL_ORIG, ISP, NB1_INTO_TOT, NB3_INTO_TOT)
                ! sort  fort.78 | uniq >sort.78.sort
                !WRITE(78,'(4I4,2F14.7)')  K1, K3, NB1_INTO_TOT, NB3_INTO_TOT, CDER_BETWEEN_STATE(1)

              !  NCV13=NB1_INTO_TOT-BD(ISP)%VBMIN+1 +(BD(ISP)%VBMAX-BD(ISP)%VBMIN+1)*(NB3_INTO_TOT-BD(ISP)%CBMIN+(BD(ISP)%CBMAX-BD(ISP)%CBMIN+1)*(K1-1))

                NCV13=BSE_INDEX%INDEX(NB1_INTO_TOT-BD(ISP)%VBMIN+1, NB3_INTO_TOT-BD(ISP)%CBMIN+1, K1, ISP -ISP_LOW+1)
                IF (NCV13==0) CYCLE ! cycle if this pair is not included                
                CTRANS=CTRANS+CDER_BETWEEN_STATE*AMAT(NCV13,LAMBDA)*WEIGHT
             ENDDO  ! enddo for N3-loop
          ENDDO     ! enddo for N1-loop
       ENDDO  ! enddo K1 loop
       ENDDO  ! enddo ISP loop
       DO I=1,3
          DO J=1,3
             CTRANS_SQUARE(I,J)=CTRANS(I)*CONJG(CTRANS(J))
          ENDDO
       ENDDO

       R_AND_INTENSITY(1,LAMBDA)= R(LAMBDA)
       R_AND_INTENSITY(2,LAMBDA)= (REAL(CTRANS_SQUARE(1,1),q)+REAL(CTRANS_SQUARE(2,2),q)+REAL(CTRANS_SQUARE(3,3),q))

       ! add to reducible polarizability X(w)
       IF (ANTIRES<0) THEN
          ! only antiresonant part (helpful for debugging)
          DO NOMEGA=1,NOMEGA_DIM
             EPS(:,:,NOMEGA)=EPS(:,:,NOMEGA)+RSPIN*WHF%WDES%WTKPT(1)*CTRANS_SQUARE* & 
                  (1/((NOMEGA-1)*DOMEGA-R(LAMBDA)+CMPLX(0,SHIFT,q)))
          ENDDO
       ELSE
          DO NOMEGA=1,NOMEGA_DIM
             EPS(:,:,NOMEGA)=EPS(:,:,NOMEGA)+RSPIN*WHF%WDES%WTKPT(1)*CTRANS_SQUARE* & 
                  (1/( (NOMEGA-1)*DOMEGA-R(LAMBDA)+CMPLX(0,SHIFT,q)) &
                  +1/(-(NOMEGA-1)*DOMEGA-R(LAMBDA)-CMPLX(0,SHIFT,q)))
          ENDDO
       ENDIF
    ENDDO lambda_bse

    CALL XML_VECARRAY("opticaltransitions")
    CALL XML_ARRAY_REAL(R_AND_INTENSITY,"(F10.4,' ')")
    CALL XML_CLOSE_TAG("varray")

    ! multiply by Coloumb kernel 4 pi e^2 / q^2 (in a.u.)
    ! to obtain v X^red 
    ! for q->0 the 1/q^2 is accounted for in the CDER_BETWEEN_STATES 
    EPS=EPS*EDEPS/LATT_CUR%OMEGA
    
    IF (ABS(WHF%WDES%VKPT(1,NQPOINT)+KPOINT_BSE(1))>1E-6_q .OR. & 
        ABS(WHF%WDES%VKPT(2,NQPOINT)+KPOINT_BSE(2))>1E-6_q .OR. & 
        ABS(WHF%WDES%VKPT(3,NQPOINT)+KPOINT_BSE(3))>1E-6_q) THEN       
       ! divide by q^2 at finite q
       DKX=(WHF%WDES%VKPT(1,NQPOINT)+KPOINT_BSE(1))*LATT_CUR%B(1,1)+ &
           (WHF%WDES%VKPT(2,NQPOINT)+KPOINT_BSE(2))*LATT_CUR%B(1,2)+ &
           (WHF%WDES%VKPT(3,NQPOINT)+KPOINT_BSE(3))*LATT_CUR%B(1,3)
       DKY=(WHF%WDES%VKPT(1,NQPOINT)+KPOINT_BSE(1))*LATT_CUR%B(2,1)+ &
           (WHF%WDES%VKPT(2,NQPOINT)+KPOINT_BSE(2))*LATT_CUR%B(2,2)+ &
           (WHF%WDES%VKPT(3,NQPOINT)+KPOINT_BSE(3))*LATT_CUR%B(2,3)
       DKZ=(WHF%WDES%VKPT(1,NQPOINT)+KPOINT_BSE(1))*LATT_CUR%B(3,1)+ &
           (WHF%WDES%VKPT(2,NQPOINT)+KPOINT_BSE(2))*LATT_CUR%B(3,2)+ &
           (WHF%WDES%VKPT(3,NQPOINT)+KPOINT_BSE(3))*LATT_CUR%B(3,3)
       
       EPS=EPS/((DKX**2+DKY**2+DKZ**2)*TPI**2)
       ! eps_mac = ( 1 + v X^red(w))^-1
       DO NOMEGA=1,NOMEGA_DIM
          EPS(1,1,NOMEGA)=1/(1+EPS(1,1,NOMEGA))
          EPS(2,2,NOMEGA)=1/(1+EPS(2,2,NOMEGA))
          EPS(3,3,NOMEGA)=1/(1+EPS(3,3,NOMEGA))
!          EPS(1,1,NOMEGA)=(1+EPS(1,1,NOMEGA))
!          EPS(2,2,NOMEGA)=(1+EPS(2,2,NOMEGA))
!          EPS(3,3,NOMEGA)=(1+EPS(3,3,NOMEGA))
!          CALL INVERSE_Z(EPS(1,1,NOMEGA),3)
       ENDDO
    ELSE
       IF (LHEAD_DIRECT) THEN
          WRITE(*,*) 'doing head direct'
          ! at q=0 have included the singular contribution
          ! eps_mac = ( 1 + v X^red(w))^-1
          DO NOMEGA=1,NOMEGA_DIM
             EPS(1,1,NOMEGA)=1/(1+EPS(1,1,NOMEGA))
             EPS(2,2,NOMEGA)=EPS(1,1,NOMEGA)
             EPS(3,3,NOMEGA)=EPS(1,1,NOMEGA)
          ENDDO
       ELSE
          ! at q=0, we have calculated EPS= v X'^red, where X'^red is the
          ! the reducible polarizability with the amputed Coulomb kernel
          ! (component G->0 removed)
          ! eps_mac = 1 - v X'^red (Onida et al. Rev. Mod. Phys. 74, 601)
          ! so, first create " - v X'^red" (do the multiplication with -1) ...
          EPS(:,:,1:NOMEGA_DIM)=-1*EPS(:,:,1:NOMEGA_DIM)
          ! ... and since we have already "absorbed" a minus sign above we have
          ! now to takes "1 + [ - v X'^red ]" on the diagonal
          DO NOMEGA=1,NOMEGA_DIM
             EPS(1,1,NOMEGA)=1+EPS(1,1,NOMEGA)
             EPS(2,2,NOMEGA)=1+EPS(2,2,NOMEGA)
             EPS(3,3,NOMEGA)=1+EPS(3,3,NOMEGA)
          ENDDO
       ENDIF
    ENDIF

    DO NOMEGA=1,NOMEGA_DIM
       OMEGA(NOMEGA)=(NOMEGA-1)*DOMEGA
       EPS_REAL(NOMEGA,:,:)=REAL(EPS(:,:,NOMEGA),q)
       EPS_IMAG(NOMEGA,:,:)=AIMAG(EPS(:,:,NOMEGA))
    ENDDO

    CALL XML_EPSILON_W(DOMEGA, EPS_REAL, EPS_IMAG, NOMEGA_DIM )

  END SUBROUTINE CALCULATE_BSE_OSZI_STRENGTH

#ifdef scaLAPACK
!**********************************************************************
!
! calculate the BSE oszillator strength
! scaLAPACK aware version
! see comments above
!
!**********************************************************************

  SUBROUTINE CALCULATE_BSE_OSZI_SCALA (WHF, NQPOINT, KPOINT_BSE, LATT_CUR, ISP_LOW, ISP_HIGH,  &
       SHIFT, BD, BSE_INDEX, R, AMAT, NCV, ANTIRES, NOMEGA_DIM, OMEGAMAX, IU6, DESC)

    USE constant
    USE full_kpoints
    USE kpoints_change
    USE pseudo
    USE mymath
    IMPLICIT NONE

    ! passed variables
    TYPE (wavespin) WHF
    INTEGER :: NQPOINT        ! q-point considered
    INTEGER :: KPOINT_BSE(3)  ! G-vector considered
    TYPE (latt) LATT_CUR
    INTEGER ISP_LOW, ISP_HIGH
    REAL(q) SHIFT   ! complex shift for broadening
    TYPE(banddesc) :: BD(:)
    TYPE(bse_matrix_index) :: BSE_INDEX
    REAL(q) :: R(:)
#ifdef single_prec_bse
    GDEFS :: AMAT(:)
#endif

#ifdef double_prec_bse
    GDEF  :: AMAT(:)
#endif
    INTEGER :: NCV
    INTEGER :: ANTIRES
    INTEGER :: NOMEGA_DIM
    REAL(q) :: OMEGAMAX
    INTEGER :: IU6
    INTEGER, OPTIONAL :: DESC(DLEN_)
 ! local 
    INTEGER :: NCV13
    INTEGER K1, K1_LOOKED_UP, K3, ISP
    INTEGER N1, N3
    INTEGER NB1_INTO_TOT, NB3_INTO_TOT
    COMPLEX(q) :: WEIGHT
    INTEGER I, J, NK1_IN_KPOINTS_FULL_ORIG
    INTEGER     :: LAMBDA, NOMEGA
    GDEF        :: CDER_BETWEEN_STATE(3)
    COMPLEX(q)  :: CTRANS(3)
    COMPLEX(q)  :: CTRANS_SQUARE(3,3)
    REAL(q)     :: DKX, DKY, DKZ
    COMPLEX(q)  :: EPS(3,3,NOMEGA_DIM)
    REAL(q)     :: EPS_REAL(NOMEGA_DIM, 3,3), EPS_IMAG(NOMEGA_DIM, 3, 3)
    REAL(q)     :: OMEGA(NOMEGA_DIM), DOMEGA
    REAL(q)     :: R_AND_INTENSITY(2,SIZE(R))
! BLACS variables
    INTEGER, EXTERNAL ::     NUMROC
    INTEGER MYROW, MYCOL, NPROW, NPCOL, NP,NQ
    INTEGER I1RES, J1RES, IROW, JCOL
    INTEGER I1, I2, J1, J2, NLAMBDA
    REAL(q)     :: RSPIN

!  three cases need to be considered
!  non-collinear RSPIN=1.0
    IF (WHF%WDES%LNONCOLLINEAR) THEN
       RSPIN=1.0_q
    ELSE
!  collinear non spin polarized: RSPIN=2.0
!  collinear spin polarized:     RSPIN=1.0 for  ISP_LOW=1 and ISP_HIGH=2
!  I believe the case ISP_HIGH=ISP_LOW for collinear spin polarized makes little sense
!  but it is not used right now anyway
       RSPIN=2._q/(ISP_HIGH-ISP_LOW+1)
    ENDIF

    CALL BLACS_GRIDINFO(DESC(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
    NP = NUMROC(NCV,DESC(MB_),MYROW,0,NPROW)
    NQ = NUMROC(NCV,DESC(NB_),MYCOL,0,NPCOL)

    EPS = 0
    K1_LOOKED_UP=-1
    NLAMBDA=0

    DOMEGA=R(NCV)/(NOMEGA_DIM-1)
    IF (OMEGAMAX>0) THEN
       DOMEGA=OMEGAMAX/(NOMEGA_DIM-1)
    ENDIF

!   loop over column indices (eigenvectors)
    DO LAMBDA=1, NCV
       ! local storage index in block (J1) and local index of block (J2)
       J1=   ((LAMBDA-DESC(NB_)*MYCOL-1)/(NPCOL*DESC(NB_)))*DESC(NB_)+1
       J2=MOD((LAMBDA-DESC(NB_)*MYCOL-1) ,NPCOL*DESC(NB_))+1

       JCOL=J1+J2-1
       IF (J2 <= 0 .OR. J2 > DESC(NB_)) THEN
          JCOL  =-1
       ELSE
          NLAMBDA=NLAMBDA+1
          IF (LAMBDA /= DESC(NB_)*MYCOL+NPCOL*(J1-1)+J2) THEN
             CALL vtutor%bug("internal error in CALCULATE_BSE_OSZI_SCALA: " // str(LAMBDA) // " " &
                // str(DESC(NB_)) // " " // str(NPCOL*(J1-1)) // " " // str(J2), __FILE__, __LINE__)
          ENDIF
       ENDIF

       IF (JCOL==-1) THEN
          ! no data on local node for eigenvector/eigenvalue pair LAMBDA_DONE
          ! simply clear CTRANS
          CTRANS=0
       ELSE
          ! local data then calculate contribution to CTRANS
          CTRANS=0
          IROW=0
          DO I1=1,NP,DESC(MB_)
             I1RES=MIN(DESC(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW=IROW+1
                
                NCV13=DESC(MB_)*MYROW+NPROW*(I1-1)+I2
                ISP= BSE_INDEX%ISP(NCV13)
              ! K1 = 1+(NCV13-1)/((BD(ISP)%VBMAX-BD(ISP)%VBMIN+1)*(BD(ISP)%CBMAX-BD(ISP)%CBMIN+1))
                K1 = BSE_INDEX%NK(NCV13)
                IF (WHF%WDES%WTKPT(K1)==0) CYCLE
                
                IF (K1/= K1_LOOKED_UP) THEN
                   NK1_IN_KPOINTS_FULL_ORIG=KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K1),KPOINTS_FULL_ORIG)
                   K1_LOOKED_UP=K1
                   K3 =KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K1)+WHF%WDES%VKPT(:,NQPOINT),KPOINTS_FULL)
                ENDIF
                      
              ! NB1_INTO_TOT=MOD(NCV13-1,BD(ISP)%VBMAX-BD(ISP)%VBMIN+1)+BD(ISP)%VBMIN
                NB1_INTO_TOT=BSE_INDEX%N1(NCV13)
                N1=NB1_INTO_TOT-BD(ISP)%VBMIN+1
                
              ! NB3_INTO_TOT=MOD((NCV13-1)/(BD(ISP)%VBMAX-BD(ISP)%VBMIN+1),(BD(ISP)%CBMAX-BD(ISP)%CBMIN+1))+BD(ISP)%CBMIN
                NB3_INTO_TOT=BSE_INDEX%N3(NCV13)
                N3=NB3_INTO_TOT-BD(ISP)%CBMIN+1
                      
              !  IF (NCV13/=NB1_INTO_TOT-BD(ISP)%VBMIN+1 +(BD(ISP)%VBMAX-BD(ISP)%VBMIN+1)*(NB3_INTO_TOT-BD(ISP)%CBMIN+(BD(ISP)%CBMAX-BD(ISP)%CBMIN+1)*(K1-1))) THEN
                IF (NCV13/=BSE_INDEX%INDEX(NB1_INTO_TOT-BD(ISP)%VBMIN+1, NB3_INTO_TOT-BD(ISP)%CBMIN+1, K1, ISP -ISP_LOW+1)) THEN
                   CALL vtutor%bug("internal error 1 in CALCULATE_BSE_OSZI_SCALA: " // str(NCV13) &
                      // " " // str(NB1_INTO_TOT-BD(ISP)%VBMIN+1) // " " // str(NB3_INTO_TOT) // " " // &
                      str(K1) // " " // str(ISP) // " " // &
                      str(BSE_INDEX%INDEX(NB1_INTO_TOT-BD(ISP)%VBMIN+1,NB3_INTO_TOT-BD(ISP)%CBMIN+1,&
                      K1,ISP -ISP_LOW+1)), __FILE__, __LINE__)
                ENDIF
                IF (FILLED_XI_ORBITAL(WHF%FERTOT(NB3_INTO_TOT,K3,ISP))) CYCLE
                IF (EMPTY_XI_ORBITAL(WHF%FERTOT(NB1_INTO_TOT,K1,ISP)))  CYCLE
                
                WEIGHT=(WHF%FERTOT(NB1_INTO_TOT,K1,ISP)-WHF%FERTOT(NB3_INTO_TOT, K3, ISP))

                ! CDER_BETWEEN_STATE =  <v_k1,n1|- i d/dq_j c_k1+q,n3> =  - <v_k1,n1| r | c_k1,n3>  for q->0
                ! or   <v_k1,n1| e-iqr  |c_k1+q,n3>  at finite q
                CALL  CDER_BETWEEN_STATES_ROTATED(CDER_BETWEEN_STATE,LATT_CUR, NK1_IN_KPOINTS_FULL_ORIG, ISP, NB1_INTO_TOT, NB3_INTO_TOT)
                CTRANS=CTRANS+CDER_BETWEEN_STATE*AMAT(IROW+(JCOL-1)*DESC(LLD_))*WEIGHT
             ENDDO ! I2
          ENDDO  ! I1
       ENDIF
       
       ! now sum CTRANS over all nodes for the current LAMBDA
       CALLMPI( M_sum_z(WHF%WDES%COMM_INTER, CTRANS, 3))
       
       DO I=1,3
          DO J=1,3
             CTRANS_SQUARE(I,J)=CTRANS(I)*CONJG(CTRANS(J))
          ENDDO
       ENDDO

       R_AND_INTENSITY(1,LAMBDA)= R(LAMBDA)
       R_AND_INTENSITY(2,LAMBDA)= (REAL(CTRANS_SQUARE(1,1),q)+REAL(CTRANS_SQUARE(2,2),q)+REAL(CTRANS_SQUARE(3,3),q))

       ! add to reducible polarizability X(w)
       IF (ANTIRES<0) THEN
          ! only antiresonant part (helpful for debugging)
          DO NOMEGA=1,NOMEGA_DIM
             EPS(:,:,NOMEGA)=EPS(:,:,NOMEGA)+RSPIN*WHF%WDES%WTKPT(1)*CTRANS_SQUARE* & 
                  (1/((NOMEGA-1)*DOMEGA-R(LAMBDA)+CMPLX(0,SHIFT,q)))
          ENDDO
       ELSE
          DO NOMEGA=1,NOMEGA_DIM
             EPS(:,:,NOMEGA)=EPS(:,:,NOMEGA)+RSPIN*WHF%WDES%WTKPT(1)*CTRANS_SQUARE* & 
                  (1/( (NOMEGA-1)*DOMEGA-R(LAMBDA)+CMPLX(0,SHIFT,q)) &
                  +1/(-(NOMEGA-1)*DOMEGA-R(LAMBDA)-CMPLX(0,SHIFT,q)))
          ENDDO
       ENDIF
    ENDDO

    CALL XML_VECARRAY("opticaltransitions")
    CALL XML_ARRAY_REAL(R_AND_INTENSITY,"(F10.4,' ')")
    CALL XML_CLOSE_TAG("varray")

    ! multiply by Coloumb kernel 4 pi e^2 / q^2 (in a.u.)
    ! to obtain v X^red 
    ! for q->0 the 1/q^2 is accounted for in the CDER_BETWEEN_STATES at q=0
    EPS=EPS*EDEPS/LATT_CUR%OMEGA

    IF (ABS(WHF%WDES%VKPT(1,NQPOINT)+KPOINT_BSE(1))>1E-6_q .OR. & 
        ABS(WHF%WDES%VKPT(2,NQPOINT)+KPOINT_BSE(2))>1E-6_q .OR. & 
        ABS(WHF%WDES%VKPT(3,NQPOINT)+KPOINT_BSE(3))>1E-6_q) THEN       
       ! divide by q^2 at finite q
       DKX=(WHF%WDES%VKPT(1,NQPOINT)+KPOINT_BSE(1))*LATT_CUR%B(1,1)+ &
           (WHF%WDES%VKPT(2,NQPOINT)+KPOINT_BSE(2))*LATT_CUR%B(1,2)+ &
           (WHF%WDES%VKPT(3,NQPOINT)+KPOINT_BSE(3))*LATT_CUR%B(1,3)
       DKY=(WHF%WDES%VKPT(1,NQPOINT)+KPOINT_BSE(1))*LATT_CUR%B(2,1)+ &
           (WHF%WDES%VKPT(2,NQPOINT)+KPOINT_BSE(2))*LATT_CUR%B(2,2)+ &
           (WHF%WDES%VKPT(3,NQPOINT)+KPOINT_BSE(3))*LATT_CUR%B(2,3)
       DKZ=(WHF%WDES%VKPT(1,NQPOINT)+KPOINT_BSE(1))*LATT_CUR%B(3,1)+ &
           (WHF%WDES%VKPT(2,NQPOINT)+KPOINT_BSE(2))*LATT_CUR%B(3,2)+ &
           (WHF%WDES%VKPT(3,NQPOINT)+KPOINT_BSE(3))*LATT_CUR%B(3,3)
       
       EPS=EPS/((DKX**2+DKY**2+DKZ**2)*TPI**2)
       ! eps_mac = ( 1 + v X^red(w))^-1
       DO NOMEGA=1,NOMEGA_DIM
          EPS(1,1,NOMEGA)=1/(1+EPS(1,1,NOMEGA))
          EPS(2,2,NOMEGA)=1/(1+EPS(2,2,NOMEGA))
          EPS(3,3,NOMEGA)=1/(1+EPS(3,3,NOMEGA))
!          EPS(1,1,NOMEGA)=(1+EPS(1,1,NOMEGA))
!          EPS(2,2,NOMEGA)=(1+EPS(2,2,NOMEGA))
!          EPS(3,3,NOMEGA)=(1+EPS(3,3,NOMEGA))
!          CALL INVERSE_Z(EPS(1,1,NOMEGA),3)
       ENDDO
    ELSE
       IF (LHEAD_DIRECT) THEN
          IF (WHF%WDES%COMM%NODE_ME==WHF%WDES%COMM%IONODE) WRITE(*,*) 'doing head direct'
          DO NOMEGA=1,NOMEGA_DIM
             EPS(1,1,NOMEGA)=1/(1+EPS(1,1,NOMEGA))
             EPS(2,2,NOMEGA)=EPS(1,1,NOMEGA)
             EPS(3,3,NOMEGA)=EPS(1,1,NOMEGA)
          ENDDO
       ELSE
          ! at q=0, we have calculated EPS= v X'^red, where X'^red is the
          ! the reducible polarizability with the amputed Coulomb kernel
          ! (component G->0 removed)
          ! epsilon_mac = 1 - v X'^red (Onida et al. Rev. Mod. Phys. 74, 601)
          ! so, first create " - v X'^red" (do the multiplication with -1) ...
          EPS(:,:,1:NOMEGA_DIM)=-1*EPS(:,:,1:NOMEGA_DIM)
          ! ... and since we have already "absorbed" a minus sign above we have
          ! now to takes "1 + [ - v X'^red ]" on the diagonal
          DO NOMEGA=1,NOMEGA_DIM
             EPS(1,1,NOMEGA)=1+EPS(1,1,NOMEGA)
             EPS(2,2,NOMEGA)=1+EPS(2,2,NOMEGA)
             EPS(3,3,NOMEGA)=1+EPS(3,3,NOMEGA)
          ENDDO
       ENDIF
    ENDIF

    DO NOMEGA=1,NOMEGA_DIM
       OMEGA(NOMEGA)=(NOMEGA-1)*DOMEGA
       EPS_REAL(NOMEGA,:,:)=REAL(EPS(:,:,NOMEGA),q)
       EPS_IMAG(NOMEGA,:,:)=AIMAG(EPS(:,:,NOMEGA))
    ENDDO

    CALL XML_EPSILON_W(DOMEGA, EPS_REAL, EPS_IMAG, NOMEGA_DIM )

  END SUBROUTINE CALCULATE_BSE_OSZI_SCALA
#endif


#ifdef use_bse_te
!****************** SUBROUTINE INIT_BSE_OSZI_STRENGTH *****************
!
! initialize the BSE oscillator strenght
! this is done quite analogously to the calculation of the 
! determination of the CALCULATE_BSE_OSZI_STRENGTH routine
!
!**********************************************************************

  SUBROUTINE INIT_BSE_OSZI_STRENGTH(WHF, LATT_CUR, ISP_LOW, ISP_HIGH, BD, BSE_INDEX )
    USE constant
    USE full_kpoints
    USE kpoints_change
    USE pseudo
    IMPLICIT NONE

    ! passed variables
    TYPE (wavespin) WHF
    TYPE (latt) LATT_CUR
    INTEGER :: ISP_LOW, ISP_HIGH
    TYPE(banddesc) :: BD(:)
    TYPE(bse_matrix_index) :: BSE_INDEX
    ! local
    INTEGER ISP
    INTEGER :: NCV13
    INTEGER K1
    INTEGER N1, N3
    INTEGER NB1_INTO_TOT, NB3_INTO_TOT
    REAL(q) :: WEIGHT
    INTEGER I, J, NK1_IN_KPOINTS_FULL_ORIG
    GDEF       :: CDER_BETWEEN_STATE(3)

    DO ISP=ISP_LOW,ISP_HIGH
    DO K1=1,WHF%WDES%NKPTS
       IF (WHF%WDES%WTKPT(K1)==0) CYCLE

       ! determined the index of this k-point in the original full k-point grid
       NK1_IN_KPOINTS_FULL_ORIG=KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K1),KPOINTS_FULL_ORIG)
       
       DO N3=1,BD(ISP)%CBMAX-BD(ISP)%CBMIN+1
          NB3_INTO_TOT=N3+BD(ISP)%CBMIN-1
          IF (FILLED_XI_ORBITAL(WHF%FERTOT(NB3_INTO_TOT,K1,ISP))) CYCLE
          DO N1=1,BD(ISP)%VBMAX-BD(ISP)%VBMIN+1
             NB1_INTO_TOT=N1+BD(ISP)%VBMIN-1
             IF (EMPTY_XI_ORBITAL(WHF%FERTOT(NB1_INTO_TOT,K1,ISP))) CYCLE
             
             WEIGHT=(WHF%FERTOT(NB1_INTO_TOT,K1,ISP)-WHF%FERTOT(NB3_INTO_TOT, K1, ISP))
             ! G-> 0 limit, consider all three directions
             CALL  CDER_BETWEEN_STATES_ROTATED(CDER_BETWEEN_STATE,LATT_CUR, NK1_IN_KPOINTS_FULL_ORIG, ISP, NB1_INTO_TOT, NB3_INTO_TOT)
             
             ! CDER_BETWEEN_STATE =  <v_k1,n1|- i d/dq_j | c_k1+q,n3>
             ! we need \int r v*_k1,n1(r) c_k1+q,n3(r) = i CDER_BETWEEN_STATE
           ! NCV13=NB1_INTO_TOT-BD(ISP)%VBMIN+1 +(BD(ISP)%VBMAX-BD(ISP)%VBMIN+1)*(NB3_INTO_TOT-BD(ISP)%CBMIN+(BD(ISP)%CBMAX-BD(ISP)%CBMIN+1)*(K1-1))
             NCV13=BSE_INDEX%INDEX(NB1_INTO_TOT-BD(ISP)%VBMIN+1, NB3_INTO_TOT-BD(ISP)%CBMIN+1, K1, ISP -ISP_LOW+1)
             IF (NCV13==0) CYCLE ! cycle if this pair is not included
!      OPTMAT (on Transmatrix)
!      was first of all multiplied with EC-EV in order to obtain the true matrix
!      element (for consistency with the old momentum-operator implementation)
!      but was then divided again by this factor upon reading in get_exciton.F;
!      Only a factor AUTOA/(2*HSQDTM) for conversion into atomic units should be
!      left here! Furthermore the k-point weight is also misplaced here!
             CDER_BETWEEN_STATE=CDER_BETWEEN_STATE*WEIGHT*AUTOA/(2._q*HSQDTM)
             OPTMAT(NCV13,:)=CDER_BETWEEN_STATE
             LDAGW(NCV13)=WHF%CELTOT(NB3_INTO_TOT,K1,ISP)-WHF%CELTOT(NB1_INTO_TOT, K1, ISP)
          ENDDO  ! enddo for N3-loop
       ENDDO     ! enddo for N1-loop
    ENDDO  ! enddo for K1-loop
    ENDDO  ! enddo spin loop

  END SUBROUTINE INIT_BSE_OSZI_STRENGTH
#endif

!**********************************************************************
!
! small subroutine that determine
! which band window is condsidered in the BSE
! In fact the window is usually different for up and
! down spin making things really difficult
!
!**********************************************************************

  SUBROUTINE SET_BAND_PARAMETERS(W, ISP, NBANDSO, NBANDSV, & 
       VBMAX, VBMIN, CBMIN, CBMAX, & 
       NSTRIPV, NSTRIP, NGLB, CBMIN4, CBMAX4, NGLB4, W1EQUALW2,  LKPOINT_PARALLEL, FIRST_EMPTYSTATE, LAST_OCCSTATE ) 
    IMPLICIT NONE
    TYPE (wavespin) W
    INTEGER :: ISP, NBANDSO, NBANDSV
    INTEGER :: VBMAX, VBMIN, CBMIN, CBMAX
    INTEGER,OPTIONAL :: NSTRIPV, NSTRIP, NGLB, CBMIN4, CBMAX4, NGLB4
    INTEGER,OPTIONAL, INTENT(IN) :: FIRST_EMPTYSTATE, LAST_OCCSTATE
    LOGICAL,OPTIONAL :: W1EQUALW2,  LKPOINT_PARALLEL
  ! local
    INTEGER :: N, NK
    
    INTEGER LAST_FILLED, FIRST_EMPTY
 
    LAST_FILLED=0
    FIRST_EMPTY=W%WDES%NB_TOT
    DO NK=1,W%WDES%NKPTS
       IF (W%WDES%WTKPT(NK)==0) CYCLE   
       LAST_FILLED=MAX(LAST_FILLED, LAST_FILLED_XI_NOMOD(W,NK,ISP))
       FIRST_EMPTY=MIN(FIRST_EMPTY, FIRST_EMPTY_XI_NOMOD(W,NK,ISP))
    ENDDO
    IF (PRESENT (FIRST_EMPTYSTATE)) THEN
        FIRST_EMPTY=FIRST_EMPTYSTATE
    ENDIF
    IF (PRESENT(LAST_OCCSTATE)) THEN
       LAST_FILLED=LAST_OCCSTATE
    ENDIF

    ! determine VBMIN and VBMAX
    VBMAX=LAST_FILLED
    VBMIN=MAX(LAST_FILLED-NBANDSO+1,1)
    IF (VBMIN > VBMAX) VBMIN=VBMAX

    CBMIN=FIRST_EMPTY
    CBMAX=MIN(FIRST_EMPTY+NBANDSV-1,W%WDES%NB_TOT)

    ! return if optional arguments are missing
    IF (.NOT. PRESENT(NSTRIPV)) RETURN

    ! NSTRIPV determines the blocking for the valence band states
    ! it can be made smaller than NBANDSO
    ! seems to work in the cases I have tested
    NSTRIPV=VBMAX-VBMIN+1
    NSTRIPV=MIN(NSTRIPV,64)

! more k-points than cores, always parallelization over k-points
! or number of k-points larger than number of conduction bands/4
#ifdef MPI
    IF (W%WDES%NKPTS >= W%WDES%NB_PAR*2 .OR. W%WDES%NKPTS*4 >= (CBMAX-CBMIN)) THEN
       LKPOINT_PARALLEL=.TRUE.
       W1EQUALW2=.FALSE.
    ELSE
       LKPOINT_PARALLEL=.FALSE.
       W1EQUALW2=.FALSE.
       IF (W%WDES%NKPTS==1) W1EQUALW2=.TRUE.
    ENDIF
#else
    LKPOINT_PARALLEL=.FALSE.
    W1EQUALW2=.FALSE.
    IF (W%WDES%NKPTS==1) W1EQUALW2=.TRUE.
#endif      
       
    IF (LKPOINT_PARALLEL) THEN
       CBMIN4=CBMIN
       CBMAX4=CBMAX
    ELSE
       ! loops over the fourth index are distributed over nodes
       ! determine data distribution
       NGLB4=(CBMAX-CBMIN)/W%WDES%NB_PAR+1
       CBMIN4=CBMIN
       DO N=1,W%WDES%NB_PAR
          CBMAX4=MIN(CBMIN4+NGLB4-1,CBMAX)
          IF (N==W%WDES%NB_LOW) EXIT
          CBMIN4=CBMAX4+1
       ENDDO
    ENDIF

    NGLB =CBMAX -CBMIN +1
    NGLB4=CBMAX4-CBMIN4+1

    ! NSTRIP determines the blocking for the conduction band states
    ! applied only internally in the TWOELECTRON4O_ACC
    ! it limits the matrix size in BLAS level three, but values around 32-64 should
    ! be fine for maximum performance
    NSTRIP= MIN(NGLB,64)

  END SUBROUTINE SET_BAND_PARAMETERS

!**********************************************************************
!
! subroutine which determines the index arrays that allow
! to determine the indices N1(valence),N2(conduction),NK,ISP
! from the global index into the BSE matrix
! also an index array is set up to determine the global index
! from the band indices, k-point and spin
!
!**********************************************************************

  SUBROUTINE SET_BSE_MATRIX_INDEX( ISP_LOW, ISP_HIGH, ALLOC_KPOINTS, BD, BSE_INDEX, W, OMEGAMAX_IN, NCV, NQPOINT)
       
    INTEGER, INTENT(IN) :: ISP_LOW, ISP_HIGH, ALLOC_KPOINTS
    TYPE (banddesc), INTENT(IN) :: BD(:)
    TYPE (bse_matrix_index), INTENT(OUT) :: BSE_INDEX
    TYPE (wavespin), INTENT(IN) :: W
    ! local
    INTEGER :: ALLOC_VALENCE_BANDS, ALLOC_CONDUCTION_BANDS, ALLOC_SPIN
    INTEGER :: N1, N3, K1, K3
    INTEGER :: ISP, INDEX
    REAL(q) :: OMEGAMAX_IN    ! maximum energy difference between conduction and valence band to be considered in BSE
    REAL(q) :: OMEGAMAX_ACT   ! actual maximum energy difference between conduction and valence band
    INTEGER :: NCV            ! number of valence conduction band pairs
    INTEGER :: NQPOINT        ! q point at which exciton is calculated 

    REAL(q) :: OMEGAMAX

    IF (OMEGAMAX_IN>0) THEN
       OMEGAMAX=OMEGAMAX_IN
    ELSE
       OMEGAMAX=1E10
    ENDIF

    OMEGAMAX_ACT=0

    ALLOC_SPIN=ISP_HIGH-ISP_LOW+1
    ALLOC_VALENCE_BANDS=0
    ALLOC_CONDUCTION_BANDS=0

    DO ISP=ISP_LOW,ISP_HIGH
       ALLOC_VALENCE_BANDS   =MAX(ALLOC_VALENCE_BANDS   ,BD(ISP)%VBMAX-BD(ISP)%VBMIN+1)
       ALLOC_CONDUCTION_BANDS=MAX(ALLOC_CONDUCTION_BANDS,BD(ISP)%CBMAX-BD(ISP)%CBMIN+1)
    ENDDO
    
    ALLOCATE(BSE_INDEX%INDEX(ALLOC_VALENCE_BANDS,ALLOC_CONDUCTION_BANDS,ALLOC_KPOINTS,ALLOC_SPIN))
    ALLOCATE(BSE_INDEX%N1(ALLOC_VALENCE_BANDS*ALLOC_CONDUCTION_BANDS*ALLOC_KPOINTS*ALLOC_SPIN))
    ALLOCATE(BSE_INDEX%N3(ALLOC_VALENCE_BANDS*ALLOC_CONDUCTION_BANDS*ALLOC_KPOINTS*ALLOC_SPIN))
    ALLOCATE(BSE_INDEX%NK(ALLOC_VALENCE_BANDS*ALLOC_CONDUCTION_BANDS*ALLOC_KPOINTS*ALLOC_SPIN))
    ALLOCATE(BSE_INDEX%ISP(ALLOC_VALENCE_BANDS*ALLOC_CONDUCTION_BANDS*ALLOC_KPOINTS*ALLOC_SPIN))
    
    BSE_INDEX%INDEX=0
    BSE_INDEX%N1=0
    BSE_INDEX%N3=0
    BSE_INDEX%NK=0
    BSE_INDEX%ISP=0
! first spin component always U
    INDEX=0
    DO ISP=ISP_LOW,ISP_HIGH
       ! K1 is the k-point index for the valence band
       DO K1=1, ALLOC_KPOINTS
          ! K3 is the k-point index for the conduction band state (unoccupied)
          K3 =KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,K1)+W%WDES%VKPT(:,NQPOINT),KPOINTS_FULL)

          DO N3=BD(ISP)%CBMIN, BD(ISP)%CBMAX
             ! safe guard, it might be that the number of bands (NB_TOT) exceeds the
             ! number of plane wave coefficients stored in NB_TOTK(NK,ISP)
             ! skip those useless bands
             IF (N3> W%WDES%NB_TOTK(K3,ISP)) CYCLE

             ! if occupancy is larger than 0.5 the band is occupied and we do not treat it as unoccupied state
             IF (FILLED_XI_ORBITAL(W%FERTOT(N3,K3,ISP))) CYCLE

             DO N1=BD(ISP)%VBMIN, BD(ISP)%VBMAX
                ! if occupancy of the valence band is smaller than 0.5 we skip the orbital as well
                IF (EMPTY_XI_ORBITAL(W%FERTOT(N1,K1,ISP))) CYCLE
                !! another safeguard: if the transition from k1 to k3 has energy difference smaller or eq to 0.0 then cycle
                IF (REAL(W%CELTOT(N3,K3,ISP)-W%CELTOT(N1,K1,ISP),q) <= 0.0_q) CYCLE
 
                OMEGAMAX_ACT=MAX(OMEGAMAX_ACT, REAL(W%CELTOT(N3,K3,ISP)-W%CELTOT(N1,K1,ISP),q))
                IF (REAL(W%CELTOT(N3,K3,ISP)-W%CELTOT(N1,K1,ISP),q) < OMEGAMAX ) THEN
                   INDEX=INDEX+1
                   BSE_INDEX%INDEX(N1-BD(ISP)%VBMIN+1,N3-BD(ISP)%CBMIN+1,K1,ISP-ISP_LOW+1)=INDEX
                   BSE_INDEX%N1(INDEX) =  N1
                   BSE_INDEX%N3(INDEX) =  N3
                   BSE_INDEX%NK(INDEX) =  K1
                   BSE_INDEX%ISP(INDEX)= ISP

                ELSE
                   BSE_INDEX%INDEX(N1-BD(ISP)%VBMIN+1,N3-BD(ISP)%CBMIN+1,K1,ISP-ISP_LOW+1)=0
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

! number of pairs equals final index
    NCV=INDEX
    BSE_INDEX%NCV=NCV

    IF (OMEGAMAX_IN<=0) THEN
       ! be shure that the highest transition is never thrown away
       ! add a small safety measure
       OMEGAMAX_IN=OMEGAMAX_ACT*(1._q+1E-5_q)
    ENDIF

  END SUBROUTINE SET_BSE_MATRIX_INDEX

!**********************************************************************
!
! subroutine which determines the index arrays that allow
! to determine the indices N1(valence),N2(conduction),NK,ISP
! from the global index into the matrix for the coulomb potential
! also an index array is set up to determine the global index
! from the band indices, k-point and spin
!
! this is different from the previous routine, since we allow
! for filled-filled as well as unoccupied-unoccupied combinations
! (which are missing in the previous routine)
!
!**********************************************************************

  SUBROUTINE SET_POT_MATRIX_INDEX( ISP_LOW, ISP_HIGH, ALLOC_KPOINTS, BD, POT_INDEX, W, NAB, NQPOINT)
       
    INTEGER, INTENT(IN) :: ISP_LOW, ISP_HIGH, ALLOC_KPOINTS
    TYPE (banddesc), INTENT(IN) :: BD(:)
    TYPE (bse_matrix_index), INTENT(OUT) :: POT_INDEX
    TYPE (wavespin), INTENT(IN) :: W
    INTEGER,INTENT(OUT)  :: NAB      ! number of valence conduction band pairs
    INTEGER,INTENT(IN) :: NQPOINT    ! q point at which exciton is calculated 
  ! local 
    INTEGER :: ALLOC_VALENCE_BANDS, ALLOC_CONDUCTION_BANDS, ALLOC_SPIN
    INTEGER :: N1, N3, K1, K3
    INTEGER :: ISP, INDEX


    ALLOC_SPIN=ISP_HIGH-ISP_LOW+1
    ALLOC_VALENCE_BANDS=0
    ALLOC_CONDUCTION_BANDS=0

    DO ISP=ISP_LOW,ISP_HIGH
       ALLOC_VALENCE_BANDS   =MAX(ALLOC_VALENCE_BANDS   ,BD(ISP)%VBMAX-BD(ISP)%VBMIN+1)
       ALLOC_CONDUCTION_BANDS=MAX(ALLOC_CONDUCTION_BANDS,BD(ISP)%CBMAX-BD(ISP)%CBMIN+1)
    ENDDO
    
    ALLOCATE(POT_INDEX%INDEX(ALLOC_VALENCE_BANDS,ALLOC_CONDUCTION_BANDS,ALLOC_KPOINTS,ALLOC_SPIN))
    ALLOCATE(POT_INDEX%N1(ALLOC_VALENCE_BANDS*ALLOC_CONDUCTION_BANDS*ALLOC_KPOINTS*ALLOC_SPIN))
    ALLOCATE(POT_INDEX%N3(ALLOC_VALENCE_BANDS*ALLOC_CONDUCTION_BANDS*ALLOC_KPOINTS*ALLOC_SPIN))
    ALLOCATE(POT_INDEX%NK(ALLOC_VALENCE_BANDS*ALLOC_CONDUCTION_BANDS*ALLOC_KPOINTS*ALLOC_SPIN))
    ALLOCATE(POT_INDEX%ISP(ALLOC_VALENCE_BANDS*ALLOC_CONDUCTION_BANDS*ALLOC_KPOINTS*ALLOC_SPIN))
    
    POT_INDEX%INDEX=0
    POT_INDEX%N1=0
    POT_INDEX%N3=0
    POT_INDEX%NK=0
    POT_INDEX%ISP=0
! first spin component always U
    INDEX=0
    DO ISP=ISP_LOW,ISP_HIGH
       ! K1 is the k-point index for the valence band
       DO K1=1, ALLOC_KPOINTS
          ! K3 is the k-point index for the conduction band state (unoccupied)
          K3 =KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,K1)+W%WDES%VKPT(:,NQPOINT),KPOINTS_FULL)

          DO N3=BD(ISP)%CBMIN, BD(ISP)%CBMAX
             ! safe guard, it might be that the number of bands (NB_TOT) exceeds the
             ! number of plane wave coefficients stored in NB_TOTK(NK,ISP)
             ! skip those useless bands
             IF (N3> W%WDES%NB_TOTK(K3,ISP)) CYCLE

             DO N1=BD(ISP)%VBMIN, BD(ISP)%VBMAX
                IF (N1> W%WDES%NB_TOTK(K3,ISP)) CYCLE
                INDEX=INDEX+1
                POT_INDEX%INDEX(N1-BD(ISP)%VBMIN+1,N3-BD(ISP)%CBMIN+1,K1,ISP-ISP_LOW+1)=INDEX
                POT_INDEX%N1(INDEX) =  N1
                POT_INDEX%N3(INDEX) =  N3
                POT_INDEX%NK(INDEX) =  K1
                POT_INDEX%ISP(INDEX)= ISP
             ENDDO
          ENDDO
       ENDDO
    ENDDO

! number of pairs equals final index
    NAB=INDEX
    POT_INDEX%NCV=NAB
  END SUBROUTINE SET_POT_MATRIX_INDEX


!**********************************************************************
!
! deallocate the BSE_INDEX array
!
!**********************************************************************

  SUBROUTINE DEALLOCATE_BSE_MATRIX_INDEX( BSE_INDEX )
    TYPE (bse_matrix_index) :: BSE_INDEX

    DEALLOCATE(BSE_INDEX%INDEX)
    DEALLOCATE(BSE_INDEX%N1)
    DEALLOCATE(BSE_INDEX%N3)
    DEALLOCATE(BSE_INDEX%NK)
    DEALLOCATE(BSE_INDEX%ISP)
    
  END SUBROUTINE DEALLOCATE_BSE_MATRIX_INDEX

!**********************************************************************
!
! diagonalize BSE matrix using one of the many possible calls
!
!**********************************************************************

  SUBROUTINE DIAG_BSE_MATRIX( NCV, AMAT, AMAT_SCALA, R, LscaLAPACKaware, MY_COMM, IU0, DESC )
    USE wave_high
    IMPLICIT NONE
    INTEGER, INTENT(IN OUT):: NCV           ! number of pairs (dimension of problem)
#ifdef single_prec_bse
    GDEFS, INTENT(IN OUT) :: AMAT(:,:)      ! in: matrix to be diagonalized
                                            ! out: eigenvectors
    GDEFS, INTENT(IN OUT), ALLOCATABLE :: AMAT_SCALA(:)  ! in: matrix distributed according to DESC
#endif

#ifdef double_prec_bse
    GDEF, INTENT(IN OUT)  :: AMAT(:,:)
    GDEF, INTENT(IN OUT), ALLOCATABLE :: AMAT_SCALA(:)
#endif
    REAL(q), INTENT(OUT)::  R(:)            ! out: eigenvalues of matrix
    LOGICAL, INTENT(IN) :: LscaLAPACKaware  ! distribute data using BLACS
    TYPE(communic), INTENT(IN) :: MY_COMM   ! communicator to be passed to the subroutine that computes the trace 
    INTEGER, INTENT(IN) :: IU0              ! VASP unit for debug messages
    INTEGER, INTENT(IN) :: DESC(DLEN_)      ! SCALAPACK matrix descriptor
! local
    REAL(q) ::  ABSTOL
    REAL(q) :: VL, VU

    INTEGER, PARAMETER :: LWORK=32
#ifdef single_prec_bse
    REAL(qs), ALLOCATABLE ::  RWORK(:)
    GDEFS, ALLOCATABLE    ::  CWRK(:)
    REAL(qs) :: R_TMP(SIZE(R))
#endif

#ifdef double_prec_bse
    REAL(q), ALLOCATABLE ::  RWORK(:)
    GDEF, ALLOCATABLE    ::  CWRK(:)
#endif
#if defined (USE_ZHEEVX) || defined (USE_ZHEEVD)
    INTEGER, ALLOCATABLE :: IWORK(:), INFOZ(:)
    INTEGER :: IL, IU, NB_CALC
#endif
#ifdef single_prec_bse
    GDEFS, ALLOCATABLE :: COVL(:,:)
#endif
#ifdef double_prec_bse
    GDEF, ALLOCATABLE :: COVL(:,:)
#endif
    INTEGER :: IFAIL, NPROW, NPCOL, MYROW, MYCOL, NP, NQ
    INTEGER, EXTERNAL    ::     NUMROC
! some [sca]LAPACK routines require parameter "ABSTOL" which controls the
! final accuracy; it is most important for single precision routines, for
! double precision routines just a "small enough number" is sufficient
! (e.g. 1.E-8). Now the "optimum parameter" is enforced for the case
! "single_prec_bse" while for "double_prec_bse" 1.E-8 would be used (the
! default) but alternatively also the optimal value for double precision
! can be enforced for "double_prec_bse" if also "high_prec_abstol" is set
!#define high_prec_abstol
#if defined (single_prec_bse)
    REAL(qs) :: SLAMCH
#ifdef scaLAPACK
    REAL(qs) :: PSLAMCH
#endif
#endif
#if defined (double_prec_bse) && defined (high_prec_abstol)
    REAL(q) :: DLAMCH
#ifdef scaLAPACK
    REAL(q) :: PDLAMCH
#endif
#endif

! the default for the double precision case ...
    ABSTOL=1.E-8_q
! this is the recommended value for maximum accuracy in the single precision
! case (always enforced for maximum safety/accuracy for "single_prec_bse"!)
#if defined (single_prec_bse)
#ifdef scaLAPACK
    ABSTOL=2._q*PSLAMCH(DESC(CTXT_),'S')
#else
    ABSTOL=2._q*SLAMCH('S')
#endif
#endif
!     this is then the corresponding recommended value for double precision
!     (only set if "high_prec_abstol" is defined, otherwise 1.E-8 is used)
#if defined (double_prec_bse) && defined (high_prec_abstol)
#ifdef scaLAPACK
    ABSTOL=2._q*PDLAMCH(DESC(CTXT_),'S')
#else
    ABSTOL=2._q*DLAMCH('S')
#endif
#endif
!==========================================================================
! version using scaLAPACK 
! matrix is either already distributed
! or redistributed using DISTRI and RECON routines in scala.F
!==========================================================================
#ifdef scaLAPACK
       IF ( .NOT. LscaLAPACKaware) THEN
          IF (IU0>0) WRITE(IU0,*) 'BSE distributing matrix using scaLAPACK distribution'
          ! distribute the global matrix AMAT to AMAT_SCALA
          CALL BLACS_GRIDINFO(DESC(CTXT_), NPROW, NPCOL,MYROW,MYCOL)
          NP = NUMROC(NCV,DESC(MB_),MYROW,0,NPROW)
          NQ = NUMROC(NCV,DESC(NB_),MYCOL,0,NPCOL)
          DEALLOCATE( AMAT_SCALA)
          ALLOCATE( AMAT_SCALA( MAX(NP*NQ,1)))
#ifdef single_prec_bse
          CALL DISTRI_SINGLE_SINGLE(AMAT,SIZE(AMAT,1),NCV,AMAT_SCALA,DESC)
#else
          CALL DISTRI(AMAT,SIZE(AMAT,1),NCV,AMAT_SCALA,DESC)
#endif

       ENDIF
#ifdef single_prec_bse
#ifdef USE_ZHEEVX
       IF (IU0>0) WRITE(IU0,*) 'BSE diagonalizing single prec. matrix (PSSYEVX_CHEEVX)'
       CALL PSSYEVX_CHEEVX_DESC(AMAT_SCALA, R_TMP, NCV, DESC, MY_COMM, ABSTOL_= ABSTOL )
#elif defined( USE_ZHEEVD ) 
       IF (IU0>0) WRITE(IU0,*) 'BSE diagonalizing double prec. matrix (PSSYEVD_CHEEVD)'
       CALL PSSYEV_CHEEVD_DESC(AMAT_SCALA, R_TMP, NCV, DESC, MY_COMM )
#else
       IF (IU0>0) WRITE(IU0,*) 'BSE diagonalizing single prec. matrix (PSSYEV_CHEEV)'
       CALL PSSYEV_CHEEV_DESC(AMAT_SCALA, R_TMP, NCV, DESC, MY_COMM )
#endif
       R=R_TMP
#endif

#ifdef double_prec_bse
#ifdef USE_ZHEEVX
       IF (IU0>0) WRITE(IU0,*) 'BSE diagonalizing double prec. matrix (PDSYEVX_ZHEEVX)'
       CALL PDSYEVX_ZHEEVX_DESC(AMAT_SCALA, R, NCV, DESC, MY_COMM, ABSTOL_= ABSTOL)
#elif defined( USE_ZHEEVD ) 
       IF (IU0>0) WRITE(IU0,*) 'BSE diagonalizing double prec. matrix (PDSYEVD_ZHEEVD)'
       CALL PDSYEV_ZHEEVD_DESC(AMAT_SCALA, R, NCV, DESC, MY_COMM)
#else
       IF (IU0>0) WRITE(IU0,*) 'BSE diagonalizing double prec. matrix (PDSYEV_ZHEEV)'
       CALL PDSYEV_ZHEEV_DESC(AMAT_SCALA, R, NCV, DESC, MY_COMM )
#endif
#endif
       IF ( .NOT. LscaLAPACKaware) THEN
          ! recombine AMAT_SCALA into AMAT and sum over all nodes
#ifdef single_prec_bse
          CALL RECON_SINGLE_SINGLE(AMAT,SIZE(AMAT,1),NCV,AMAT_SCALA,DESC)
#ifdef gammareal
          CALLMPI( M_sum_single(MY_COMM, AMAT  , SIZE(AMAT)))
#else
          CALLMPI( M_sum_single(MY_COMM, AMAT  , 2*SIZE(AMAT)))
#endif
#else
          CALL RECON(AMAT,SIZE(AMAT,1),NCV,AMAT_SCALA,DESC)
          CALLMPI( M_sum_g(MY_COMM, AMAT  , SIZE(AMAT)))
#endif

          DEALLOCATE(AMAT_SCALA)
          ALLOCATE(AMAT_SCALA(1))
       ENDIF
#else


#ifdef USE_ZHEEVX 
!==========================================================================
! version using LAPACK ZHEEVX
!==========================================================================
       VL=0 ; VU=0 ; IL=0 ; IU=0
       ALLOCATE(COVL(NCV,NCV), CWRK(NCV*LWORK), RWORK(7*NCV), IWORK(5*NCV), INFOZ(NCV))
#ifdef  gammareal
#ifdef single_prec_bse
       IF (IU0>=0) WRITE(IU0,*) 'BSE diagonalizing single prec. matrix (SSYEVX)'
       CALL SSYEVX( 'V', 'A', 'U', NCV, AMAT(1,1) , NCV, VL, VU, IL, IU, &
            ABSTOL , NB_CALC , R_TMP, COVL(1,1), NCV, CWRK, &
            LWORK*NCV, RWORK, IWORK, INFOZ, IFAIL )
       R=R_TMP
#else
       IF (IU0>=0) WRITE(IU0,*) 'BSE diagonalizing double prec. matrix (DSYEVX)'
       CALL DSYEVX( 'V', 'A', 'U', NCV, AMAT(1,1) , NCV, VL, VU, IL, IU, &
            ABSTOL , NB_CALC , R, COVL(1,1), NCV, CWRK, &
            LWORK*NCV, RWORK, IWORK, INFOZ, IFAIL )         
#endif
#else
#ifdef single_prec_bse
       IF (IU0>=0) WRITE(IU0,*) 'BSE diagonalizing single prec. matrix (CHEEVX)'
       CALL CHEEVX( 'V', 'A', 'U', NCV, AMAT(1,1) , NCV, VL, VU, IL, IU, &
            ABSTOL , NB_CALC , R_TMP, COVL(1,1), NCV, CWRK, &
            LWORK*NCV, RWORK, IWORK, INFOZ, IFAIL )
       R=R_TMP
#else
       IF (IU0>=0) WRITE(IU0,*) 'BSE diagonalizing double prec. matrix (ZHEEVX)'
       CALL ZHEEVX( 'V', 'A', 'U', NCV, AMAT(1,1) , NCV, VL, VU, IL, IU, &
            ABSTOL , NB_CALC , R, COVL(1,1), NCV, CWRK, &
            LWORK*NCV, RWORK, IWORK, INFOZ, IFAIL )
#endif
#endif
       AMAT=COVL   !overwriting eigenvectors to the original matrix
       DEALLOCATE(COVL, CWRK, RWORK, IWORK, INFOZ)
#else
!==========================================================================
! simplest version LAPACK ZHEEV
!==========================================================================
       ALLOCATE( CWRK(NCV*LWORK), RWORK(3*NCV))
#ifdef  gammareal
#ifdef  single_prec_bse
       IF (IU0>=0) WRITE(IU0,*) 'BSE diagonalizing single prec. matrix (SSYEV)'

       CALL SSYEV &
            ('V','U',NCV, AMAT(1,1), NCV, &
            R_TMP, CWRK, LWORK*NCV, IFAIL)
       R=R_TMP
#else
       IF (IU0>=0) WRITE(IU0,*) 'BSE diagonalizing double prec. matrix (DSYEV)'
       CALL DSYEV &
            ('V','U',NCV, AMAT(1,1), NCV, &
            R, CWRK, LWORK*NCV, IFAIL)
#endif
#else
#ifdef  single_prec_bse
       IF (IU0>=0) WRITE(IU0,*) 'BSE diagonalizing single prec. matrix (CHEEV)'

       CALL CHEEV &
            ('V','U',NCV, AMAT(1,1), NCV, &
            R_TMP, CWRK, LWORK*NCV, RWORK,  IFAIL)
       R=R_TMP
#else
       IF (IU0>0) WRITE(IU0,*) 'BSE diagonalizing double prec. matrix (ZHEEV)'
    
       CALL ZHEEV &
            ('V','U',NCV, AMAT(1,1), NCV, &
            R, CWRK, LWORK*NCV, RWORK,  IFAIL)
#endif
#endif
       DEALLOCATE( CWRK, RWORK )
#endif

#endif
  END SUBROUTINE DIAG_BSE_MATRIX

!**********************************************************************
!
! left and right multiply the matrix A+B by (A-B)^(1/2)
! compare Sander, Maggio, Kresse, PRB 045209 (2015)  Equ. (40)
! 
! the matrix A+B is supplied in A
! the matrix A-B is supplied as U R U+
!
! the transformed  matrix is therefore given by
! C = U R U+ A U R U+
! to save storage four matrix matrix products are calculated
! only one temporary matrix is used
!
!**********************************************************************

  SUBROUTINE TDA_BSE( NCV, A, U, RB ,IU0)
    USE wave_high
    INTEGER, INTENT(IN)   :: NCV            ! number of pairs (dimension of problem)
    REAL(q), INTENT(IN)   :: RB(:)          ! eigenvalues of (A-B)
#ifdef single_prec_bse
    GDEFS, INTENT(IN OUT) :: A(:,:)         ! in: matrix A
                                            ! out:  U R U+ A U R U+
    GDEFS, INTENT(IN)     :: U(:,:)         ! eigenvectors of (A-B)
    GDEFS,ALLOCATABLE :: T(:,:)
#else
    GDEF, INTENT(IN OUT)  :: A(:,:)
    GDEF, INTENT(IN)      :: U(:,:)
    GDEF,ALLOCATABLE :: T(:,:)
#endif
    INTEGER :: I, J
    REAL(q) :: SR(SIZE(RB))
    INTEGER :: IU0

    ALLOCATE(T(NCV,NCV))

    !  U+ (A+B) -> T
    CALL GGEMMbse( trans,'N', NCV, NCV, NCV, onebse, & 
         U(1,1), NCV, A(1,1), NCV, zerobse, T(1,1), NCV)
    !  [U+ (A+B)] U -> A
    CALL GGEMMbse( 'N','N', NCV, NCV, NCV, onebse, & 
         T(1,1), NCV, U(1,1), NCV, zerobse, A(1,1), NCV)
    ! multiply with sqrt of eigenvalues from left and right

    !SR=SQRT(RB)
    CALL SQRT_EIGENVALUE(RB ,SR, IU0 )

    ! TODO gK: we need one subroutine to handle negative eigenvalues in RB
    ! the way it is handled here is not even remotely gracefull
    DO J=1,NCV
       DO I=1,NCV
          A(I,J)=A(I,J)*SR(I)*SR(J)
       ENDDO
    ENDDO
    !  U A -> T
    CALL GGEMMbse( 'N', 'N', NCV, NCV, NCV, onebse, & 
         U(1,1), NCV, A(1,1), NCV, zerobse, T(1,1), NCV)
    !  [U A ] U+ -> A
    CALL GGEMMbse( 'N',trans, NCV, NCV, NCV, onebse, & 
         T(1,1), NCV, U(1,1), NCV, zerobse, A(1,1), NCV)

    DEALLOCATE(T)

  END SUBROUTINE TDA_BSE

#ifdef scaLAPACK

!**********************************************************************
!
! same is TDA_BSE but using SCALAPACK (see documentation there)
!
!**********************************************************************

 SUBROUTINE TDA_BSE_SCA( NCV, A, U, RB , IU0, DESC)
    USE wave_high
    INTEGER, INTENT(IN) :: NCV
    REAL(q), INTENT(IN) :: RB(:)
#ifdef single_prec_bse
    GDEFS, INTENT(IN OUT) :: A(:)                       
    GDEFS, INTENT(IN) :: U(:)
    GDEFS,ALLOCATABLE :: T(:)
#else
    GDEF, INTENT(IN OUT)  :: A(:)
    GDEF, INTENT(IN)  :: U(:)
    GDEF,ALLOCATABLE :: T(:)
#endif
    INTEGER :: I, J
    REAL(q) :: SR(SIZE(RB))
    INTEGER :: IU0
! BLACS variables
    INTEGER, EXTERNAL ::     NUMROC
    INTEGER MYROW, MYCOL, NPROW, NPCOL, NP, NQ
    INTEGER I1RES,J1RES,IROW,JCOL
    INTEGER I1,I2,J1,J2
    INTEGER :: DESC(DLEN_) ! descriptor set up for matrix distribution
    CALL BLACS_GRIDINFO(DESC(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

    NP = NUMROC(NCV,DESC(MB_),MYROW,0,NPROW)
    NQ = NUMROC(NCV,DESC(NB_),MYCOL,0,NPCOL)

    ALLOCATE(T(SIZE(A)))

    !  U+ (A+B) -> T
    CALL GPGEMMbse( trans,'N', NCV, NCV, NCV, onebse, & 
         U(1), 1, 1, DESC, A(1), 1, 1, DESC, zerobse, T(1), 1, 1, DESC)
    !  [U+ (A+B)] U -> A
    CALL GPGEMMbse( 'N','N', NCV, NCV, NCV, onebse, & 
         T(1), 1, 1, DESC, U(1), 1, 1, DESC, zerobse, A(1), 1, 1, DESC)

    ! multiply with sqrt of eigenvalues from left and right

    !SR=SQRT(RB)
    CALL SQRT_EIGENVALUE(RB ,SR, IU0 )

    JCOL=0
    DO J1=1,NQ,DESC(NB_)
       J1RES=MIN(DESC(NB_),NQ-J1+1)
       DO J2=1,J1RES
          JCOL=JCOL+1
          J=DESC(NB_)*MYCOL+NPCOL*(J1-1)+J2

          IROW=0
          DO I1=1,NP,DESC(MB_)
             I1RES=MIN(DESC(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW=IROW+1
                I=DESC(MB_)*MYROW+NPROW*(I1-1)+I2
          ! A(I,J)=A(I,J)*SR(I)*SR(J)
                A(IROW+(JCOL-1)*DESC(LLD_))=A(IROW+(JCOL-1)*DESC(LLD_))*SR(I)*SR(J)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !  U A -> T
    CALL GPGEMMbse( 'N', 'N', NCV, NCV, NCV, onebse, & 
         U(1), 1, 1, DESC, A(1), 1, 1, DESC, zerobse, T(1), 1, 1, DESC)
    !  [U A ] U+ -> A
    CALL GPGEMMbse( 'N',trans, NCV, NCV, NCV, onebse, & 
         T(1), 1, 1, DESC, U(1), 1, 1, DESC, zerobse, A(1), 1, 1, DESC)

    DEALLOCATE(T)

  END SUBROUTINE TDA_BSE_SCA
#endif

!**********************************************************************
!
! small helper routine to calculate the sqrt(R)
! and prompt a warning if R holds negative values
!
!**********************************************************************

  SUBROUTINE SQRT_EIGENVALUE(R ,SR, IU0 )
    REAL(q), INTENT(IN) :: R(:)           ! eigenvalues
    REAL(q), INTENT(OUT):: SR(:)          ! sqrt of eigenvalues
    INTEGER, INTENT(IN) :: IU0            ! unit for error logging
  ! local
    INTEGER :: I, NSR

    NSR=0

    !SR =SQRT(R)
    DO I=1,SIZE(R)
       IF (R(I)<=0.0) THEN
          SR(I)=0.0_q
          NSR=NSR+1
       ELSE
          SR(I) =SQRT(R(I))
       ENDIF
    ENDDO

    IF (NSR>0 .AND. IU0>=0 ) THEN 
       WRITE(IU0,*) 'WARNING in SQRT_EIGENVALUE: some eigenvalues in R are negative',NSR
       WRITE(IU0,*) ' this indicates instabilities'
    ENDIF
   END SUBROUTINE SQRT_EIGENVALUE


!**********************************************************************
!
! small helper routine to calculate the sqrt(R)
! and prompt a warning if R holds negative values
!
! this routine is used to determine the eigenvalues
! of the original problem Omega_lambda from the eigenvalues of the
! squared problem (Omega_lambda^2)
! compare Sander, Maggio, Kresse, PRB 045209 (2015)  Equ. (42)
!
!**********************************************************************

  SUBROUTINE SQRT_EIGENVALUE_INPLACE(R , IU0 )
    REAL(q) :: R(:)
    REAL(q), ALLOCATABLE :: SR(:)
    INTEGER :: I, NSR
    INTEGER :: IU0

    ALLOCATE(SR(SIZE(R)))
    NSR=0

    !SR =SQRT(R)
    DO I=1,SIZE(R)
       IF (R(I)<=0.0) THEN
          SR(I)=0.0_q
          NSR=NSR+1
       ELSE
          SR(I) =SQRT(R(I))
       ENDIF
    ENDDO

    IF (NSR>0 .AND. IU0>=0 ) THEN 
       WRITE(IU0,*) 'WARNING in SQRT_EIGENVALUE_INPLACE: some eigenvalues of the squared problem are negative',NSR
       WRITE(IU0,*) ' this indicates instabilities'
    ENDIF
    R(:)=SR(:)

    DEALLOCATE(SR)

   END SUBROUTINE SQRT_EIGENVALUE_INPLACE



!**********************************************************************
!
! determine the eigenvectors of the original problem e.g.
!  A  B    X     w  0   X
!             =          
!  B  A    Y     0  -w  Y
!
! it is fairly easy to see that
!  (X+Y) = (A-B)^( 1/2) Z R^(-1/2)
!  (X-Y) = (A-B)^(-1/2) Z R^1/2
! compare Sander, Maggio, Kresse, PRB 045209 (2015)  Eqs. (44-45)
!
! where Z are the eigenvectors of the 
! squared problem [(A-B)^(1/2) (A+B) (A-B)^(1/2) Z = Z R]
! and R the eigenvalues omega_lambda of the original problem
! furthermore one can show the
!  (X+Y) (X-Y)+  = 1
! and 
!  X Y+  - Y X+  = 0
! and 
!  X X+ - Y Y+   = 1
!
! input:
! Z=z (eigenvectors of squared problem) and R
! A-B is passed as A-B  = U RB U+
!
! the first routine TDA_BSE_EIGENVECTORS returns only
! Z = X+Y, which is required to calculate the dielectric function
! compare Sander, Maggio, Kresse, PRB 045209 (2015)  Eqs. (58)
!
! the second routine BSE_EIGENVECTORS_FULL 
! returns X and Y, but requires additional work space
!
!**********************************************************************
    
  SUBROUTINE TDA_BSE_EIGENVECTORS( NCV, U, RB, Z,  R , IU0)
    USE wave_high
    INTEGER, INTENT(IN) :: NCV            ! number of pairs (dimension of problem)
    REAL(q), INTENT(IN) :: RB(:)          ! eigenvalues of B and squared problem
    REAL(q), INTENT(IN) :: R(:)           ! sqrt of eigenvalues of squared problem 
    INTEGER, INTENT(IN) :: IU0            ! unit for error logging
#ifdef single_prec_bse
    GDEFS, INTENT(IN OUT) :: Z(:,:)       ! in:  eigenvectors of squared problem
                                          ! out: X+Y
    GDEFS, INTENT(IN OUT) :: U(:,:)       ! in: eigenvectors of A-B 
                                          ! out: 0 
    GDEFS,ALLOCATABLE :: T(:,:)           ! auxilary
#else
    GDEF, INTENT(IN OUT)   :: Z(:,:)
    GDEF, INTENT(IN OUT)   :: U(:,:)      
    GDEF,ALLOCATABLE :: T(:,:)
#endif
! local
    INTEGER :: I, J
    REAL(q) :: SRB(SIZE(R)),SR(SIZE(R))
    INTEGER :: NSR, NSRB

    ALLOCATE(T(NCV,NCV))

    !  U+ Z -> T2
    CALL GGEMMbse( trans,'N', NCV, NCV, NCV, onebse, & 
         U(1,1), NCV, Z(1,1), NCV, zerobse, T(1,1), NCV)
    NSR=0
    !SR =1/SQRT(R)
    DO I=1,SIZE(R)

       IF (R(I)<=0.0_q ) THEN
          SR(I)=0.0_q
          NSR=NSR+1
       ELSE
          SR(I) =1/SQRT(R(I))
       ENDIF
    ENDDO

    NSRB=0
    DO I=1,SIZE(R)
       IF (RB(I)<=0.0_q ) THEN
          SRB(I)=0.0_q
          NSRB=NSRB+1
       ELSE
          SRB(I) =SQRT(RB(I))
       ENDIF
    ENDDO

    ! T1= SQRT(RB) U+ Z 1/SQRT(R)
    DO J=1,NCV
       DO I=1,NCV
          T(I,J)=SRB(I)*T(I,J)*SR(J)
       ENDDO
    ENDDO

    !  Z =(X+Y)= U T1 -> A
    CALL GGEMMbse( 'N', 'N', NCV, NCV, NCV, onebse, & 
         U(1,1), NCV, T(1,1), NCV, zerobse, Z(1,1), NCV)
    U=0

    DEALLOCATE(T)

    IF (NSR>0 .AND. IU0>=0) WRITE(*,'(1X,A,I5)') 'WARNING in TDA_BSE_EIGENVECTORS: some eigenvalues in R are negative',NSR
    IF (NSRB>0.AND. IU0>=0) WRITE(*,'(1X,A,I5)') 'WARNING in TDA_BSE_EIGENVECTORS: some eigenvalues in RB are negative',NSRB
      
  END SUBROUTINE TDA_BSE_EIGENVECTORS


  SUBROUTINE BSE_EIGENVECTORS_FULL( NCV, U, RB, Z,  R , IU0 )
    USE wave_high
    INTEGER, INTENT(IN) :: NCV            ! number of pairs (dimension of problem)
    REAL(q), INTENT(IN) :: RB(:), R(:)    ! eigenvalues of B and squared problem
    INTEGER, INTENT(IN) :: IU0            ! unit for error logging
#ifdef single_prec_bse
    GDEFS, INTENT(IN OUT) :: Z(:,:)       ! in:  eigenvectors of squared problem
                                          ! out: X
    GDEFS, INTENT(IN OUT) :: U(:,:)       ! eigenvectors of A-B 
                                          ! out: Y
    GDEFS,ALLOCATABLE :: T1(:,:),T2(:,:)  ! auxilary
#else
    GDEF, INTENT(IN OUT)   :: Z(:,:)
    GDEF, INTENT(IN OUT)   :: U(:,:)      
    GDEF,ALLOCATABLE :: T1(:,:),T2(:,:)
#endif
! local
    INTEGER :: I, J
    REAL(q) :: SRB(SIZE(R)),SR(SIZE(R))
    INTEGER :: NSR, NSRB

    ALLOCATE(T1(NCV,NCV),T2(NCV,NCV))

    !  U+ Z -> T2
    CALL GGEMMbse( trans,'N', NCV, NCV, NCV, onebse, & 
         U(1,1), NCV, Z(1,1), NCV, zerobse, T2(1,1), NCV)


    ! T1= SQRT(RB) U+ Z 1/SQRT(R)

    NSR=0
    !SR =1/SQRT(R)
    DO I=1,SIZE(R)
       IF (R(I)<=0.0) THEN
          SR(I)=0.0_q
          NSR=NSR+1
       ELSE
          SR(I) =1/SQRT(R(I))
       ENDIF
    ENDDO

    NSRB=0
    DO I=1,SIZE(R)
       IF (RB(I)<=0.0) THEN
          SRB(I)=0.0_q
          NSRB=NSRB+1
       ELSE
          SRB(I) =SQRT(RB(I))
       ENDIF
    ENDDO

    DO I=1,NCV
       DO J=1,NCV
          T1(I,J)=SRB(I)*T2(I,J)*SR(J)
       ENDDO
    ENDDO

    ! T2= 1/SQRT(RB) U+ Z SQRT(R)
    DO I=1,SIZE(R)
       IF (R(I)<=0.0) THEN
          SR(I)=0.0_q
       ELSE
          SR(I) =SQRT(R(I))
       ENDIF
    ENDDO

    !SRB=1/SQRT(RB)
    DO I=1,SIZE(R)
       IF (RB(I)<=0.0) THEN
          SRB(I)=0.0_q
       ELSE
          SRB(I) =1/SQRT(RB(I))
       ENDIF
    ENDDO

    DO J=1,NCV
       DO I=1,NCV
          T2(I,J)=SRB(I)*T2(I,J)*SR(J)
       ENDDO
    ENDDO

    !  Z =(X+Y)= U T1 -> A
    CALL GGEMMbse( 'N', 'N', NCV, NCV, NCV, onebse, & 
         U(1,1), NCV, T1(1,1), NCV, zerobse, Z(1,1), NCV)

    !  T1=(X-Y)= U T2 -> A
    CALL GGEMMbse( 'N', 'N', NCV, NCV, NCV, onebse, & 
         U(1,1), NCV, T2(1,1), NCV, zerobse, T1(1,1), NCV)

    !  test multiply Z x U should be equal unity matrix
    ! U=T1
    ! CALL GGEMMbse( 'N', trans, NCV, NCV, NCV, onebse, & 
    !     Z(1,1), NCV, T1(1,1), NCV, zerobse, T2(1,1), NCV)
    !DO I=1,NCV
    !   WRITE(77,'(10E20.10)') T2(:,I)
    !ENDDO

    ! U = Y
    U=(Z-T1)/2
    ! Z = X
    Z=(Z+T1)/2

    !WRITE(*,*) 'X:'
    !WRITE(*,*) Z(:,:)
    !WRITE(*,*) 'Y:'
    !WRITE(*,*) U(:,:)
    ! test calculate X X+ - Y Y+  should be equal unity matrix
    ! CALL GGEMMbse( 'N', trans, NCV, NCV, NCV, onebse, & 
    !   Z(1,1), NCV, Z(1,1), NCV, zerobse, T1(1,1), NCV)
    ! CALL GGEMMbse( 'N', trans, NCV, NCV, NCV, -onebse, & 
    !   U(1,1), NCV, U(1,1), NCV, onebse, T1(1,1), NCV)
    ! DO I=1,NCV
    !    WRITE(77,'(10E20.10)') T1(:,I)
    ! ENDDO

    DEALLOCATE(T1,T2)

    IF (NSR>0  .AND. IU0>=0) WRITE(*,'(1X,A,I5)') 'WARNING in BSE_EIGENVECTORS_FULL: some eigenvalues in R are negative',NSR
    IF (NSRB>0 .AND. IU0>=0) WRITE(*,'(1X,A,I5)') 'WARNING in BSE_EIGENVECTORS_FULL: some eigenvalues in RB are negative',NSRB
      
  END SUBROUTINE BSE_EIGENVECTORS_FULL

#ifdef scaLAPACK
!**********************************************************************
!
! same as TDA_BSE_EIGENVECTORS but using SCALAPACK (see documentation there)
!
! input:
! Z=z (eigenvectors of squared problem) and R
! A-B is passed as A-B  = U RB U+
!
! the first routine TDA_BSE_EIGENVECTORS_SCA returns only
! Z = X+Y, which is required to calculate the dielectric function
!
! the second routine BSE_EIGENVECTORS_FULL_SCA 
! returns X and Y, but requires additional work space
!
!**********************************************************************


  SUBROUTINE TDA_BSE_EIGENVECTORS_SCA( NCV, U, RB, Z,  R , IU0, DESC)
    USE wave_high
    INTEGER, INTENT(IN) :: NCV            ! number of pairs (dimension of problem)
    REAL(q), INTENT(IN) :: RB(:)          ! eigenvalues of B
    REAL(q), INTENT(IN) :: R(:)           ! sqrt of eigenvalues of squared problem 
    INTEGER, INTENT(IN) :: IU0            ! unit for error logging
    INTEGER, INTENT(IN) :: DESC(DLEN_)    ! SCALAPACK matrix descriptor
#ifdef single_prec_bse
    GDEFS, INTENT(IN OUT) :: Z(:)         ! in: eigenvectors of squared problem
                                          ! out: X+Y
    GDEFS, INTENT(IN OUT) :: U(:)         ! in: eigenvectors of A-B 
                                          ! out: 0 
    GDEFS,ALLOCATABLE :: T(:)             ! auxilary
#else
    GDEF, INTENT(IN OUT)   :: Z(:)
    GDEF, INTENT(IN OUT)   :: U(:)      
    GDEF,ALLOCATABLE :: T(:)
#endif
! local
    INTEGER :: I, J
    REAL(q) :: SRB(SIZE(R)),SR(SIZE(R))
    INTEGER :: NSR, NSRB
! BLACS variables
    INTEGER, EXTERNAL ::     NUMROC
    INTEGER MYROW, MYCOL, NPROW, NPCOL, NP, NQ
    INTEGER I1RES,J1RES,IROW,JCOL
    INTEGER I1,I2,J1,J2
    CALL BLACS_GRIDINFO(DESC(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

    NP = NUMROC(NCV,DESC(MB_),MYROW,0,NPROW)
    NQ = NUMROC(NCV,DESC(NB_),MYCOL,0,NPCOL)

    ALLOCATE(T(SIZE(Z)))

    !  U+ Z -> T2
    CALL GPGEMMbse( trans,'N', NCV, NCV, NCV, onebse, & 
         U(1), 1, 1, DESC, Z(1), 1, 1, DESC, zerobse, T(1), 1, 1, DESC)
    NSR=0
    !SR =1/SQRT(R)
    DO I=1,SIZE(R)

       IF (R(I)<=0.0_q ) THEN
          SR(I)=0.0_q
          NSR=NSR+1
       ELSE
          SR(I) =1/SQRT(R(I))
       ENDIF
    ENDDO

    NSRB=0
    DO I=1,SIZE(R)
       IF (RB(I)<=0.0_q ) THEN
          SRB(I)=0.0_q
          NSRB=NSRB+1
       ELSE
          SRB(I) =SQRT(RB(I))
       ENDIF
    ENDDO

    JCOL=0
    DO J1=1,NQ,DESC(NB_)
       J1RES=MIN(DESC(NB_),NQ-J1+1)
       DO J2=1,J1RES
          JCOL=JCOL+1
          J=DESC(NB_)*MYCOL+NPCOL*(J1-1)+J2

          IROW=0
          DO I1=1,NP,DESC(MB_)
             I1RES=MIN(DESC(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW=IROW+1
                I=DESC(MB_)*MYROW+NPROW*(I1-1)+I2
          ! T(I,J)=SRB(I)*T(I,J)*SR(J)
                T(IROW+(JCOL-1)*DESC(LLD_))=T(IROW+(JCOL-1)*DESC(LLD_))*SRB(I)*SR(J)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !  Z =(X+Y)= U T1 -> A
    CALL GPGEMMbse( 'N', 'N', NCV, NCV, NCV, onebse, & 
         U(1), 1, 1, DESC, T(1), 1, 1, DESC, zerobse, Z(1), 1, 1, DESC)
    U=0

    DEALLOCATE(T)

    IF (NSR>0 .AND. IU0>=0) WRITE(*,'(1X,A,I5)') 'WARNING in TDA_BSE_EIGENVECTORS: some eigenvalues in R are negative',NSR
    IF (NSRB>0.AND. IU0>=0) WRITE(*,'(1X,A,I5)') 'WARNING in TDA_BSE_EIGENVECTORS: some eigenvalues in RB are negative',NSRB

  END SUBROUTINE TDA_BSE_EIGENVECTORS_SCA


  SUBROUTINE BSE_EIGENVECTORS_FULL_SCA( NCV, U, RB, Z,  R , IU0, DESC)
! input:
! Z=z (eigenvectors of squared problem) and R eigenvalues of the original problem
! A-B is passed as A-B  = U RB U+
    USE wave_high
    INTEGER, INTENT(IN) :: NCV            ! number of pairs (dimension of problem)
    REAL(q), INTENT(IN) :: RB(:), R(:)    ! eigenvalues of B and squared problem
    INTEGER, INTENT(IN) :: IU0            ! unit for error logging
    INTEGER, INTENT(IN) :: DESC(DLEN_)    ! SCALAPACK matrix descriptor
#ifdef single_prec_bse
    GDEFS, INTENT(IN OUT) :: Z(:)         ! in: eigenvectors of squared problem
                                          ! out: X
    GDEFS, INTENT(IN OUT) :: U(:)         ! in: eigenvectors of A-B 
                                          ! out: Y 
    GDEFS,ALLOCATABLE :: T1(:), T2(:)     ! auxilary
#else
    GDEF, INTENT(IN OUT)   :: Z(:)
    GDEF, INTENT(IN OUT)   :: U(:)      
    GDEF,ALLOCATABLE :: T1(:), T2(:)
#endif
! local
    INTEGER :: I, J
    REAL(q) :: SRB(SIZE(R)),SR(SIZE(R))
    INTEGER :: NSR, NSRB
! BLACS variables
    INTEGER, EXTERNAL ::     NUMROC
    INTEGER MYROW, MYCOL, NPROW, NPCOL, NP, NQ
    INTEGER I1RES,J1RES,IROW,JCOL
    INTEGER I1,I2,J1,J2
    CALL BLACS_GRIDINFO(DESC(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

    NP = NUMROC(NCV,DESC(MB_),MYROW,0,NPROW)
    NQ = NUMROC(NCV,DESC(NB_),MYCOL,0,NPCOL)


    ALLOCATE(T1(SIZE(Z)),T2(SIZE(Z)))

    !  U+ Z -> T2
    CALL GPGEMMbse( trans,'N', NCV, NCV, NCV, onebse, & 
         U(1), 1, 1, DESC, Z(1), 1, 1, DESC, zerobse, T2(1), 1, 1, DESC)

    NSR=0
    !SR =1/SQRT(R)
    DO I=1,SIZE(R)
       IF (R(I)<=0.0) THEN
          SR(I)=0.0_q
          NSR=NSR+1
       ELSE
          SR(I) =1/SQRT(R(I))
       ENDIF
    ENDDO

    NSRB=0
    DO I=1,SIZE(R)
       IF (RB(I)<=0.0) THEN
          SRB(I)=0.0_q
          NSRB=NSRB+1
       ELSE
          SRB(I) =SQRT(RB(I))
       ENDIF
    ENDDO

    JCOL=0
    DO J1=1,NQ,DESC(NB_)
       J1RES=MIN(DESC(NB_),NQ-J1+1)
       DO J2=1,J1RES
          JCOL=JCOL+1
          J=DESC(NB_)*MYCOL+NPCOL*(J1-1)+J2

          IROW=0
          DO I1=1,NP,DESC(MB_)
             I1RES=MIN(DESC(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW=IROW+1
                I=DESC(MB_)*MYROW+NPROW*(I1-1)+I2
          ! T1(I,J)=SRB(I)*T2(I,J)*SR(J)
                T1(IROW+(JCOL-1)*DESC(LLD_))=T2(IROW+(JCOL-1)*DESC(LLD_))*(SRB(I)*SR(J))
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! T2= 1/SQRT(RB) U+ Z SQRT(R)
    DO I=1,SIZE(R)
       IF (R(I)<=0.0) THEN
          SR(I)=0.0_q
       ELSE
          SR(I) =SQRT(R(I))
       ENDIF
    ENDDO

    !SRB=1/SQRT(RB)
    DO I=1,SIZE(R)
       IF (RB(I)<=0.0) THEN
          SRB(I)=0.0_q
       ELSE
          SRB(I) =1/SQRT(RB(I))
       ENDIF
    ENDDO

    JCOL=0
    DO J1=1,NQ,DESC(NB_)
       J1RES=MIN(DESC(NB_),NQ-J1+1)
       DO J2=1,J1RES
          JCOL=JCOL+1
          J=DESC(NB_)*MYCOL+NPCOL*(J1-1)+J2

          IROW=0
          DO I1=1,NP,DESC(MB_)
             I1RES=MIN(DESC(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW=IROW+1
                I=DESC(MB_)*MYROW+NPROW*(I1-1)+I2
                T2(IROW+(JCOL-1)*DESC(LLD_))=T2(IROW+(JCOL-1)*DESC(LLD_))*(SRB(I)*SR(J))
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !  Z =(X+Y)= U T1 -> A
    CALL GPGEMMbse( 'N', 'N', NCV, NCV, NCV, onebse, & 
         U(1), 1, 1, DESC, T1(1), 1, 1, DESC, zerobse, Z(1), 1, 1, DESC)

    !  T1=(X-Y)= U T2 -> A
    CALL GPGEMMbse( 'N', 'N', NCV, NCV, NCV, onebse, & 
         U(1), 1, 1, DESC, T2(1), 1, 1, DESC, zerobse, T1(1), 1, 1, DESC)

    !  test multiply Z x U should be equal unity matrix
    ! U=T1
    ! CALL GPGEMMbse( 'N', trans, NCV, NCV, NCV, onebse, & 
    !     Z(1), 1, 1, DESC, T1(1), 1, 1, DESC, zerobse, T2(1), 1, 1, DESC)
    !DO I=1,NCV
    !   WRITE(77,'(10E20.10)') T2(:,I)
    !ENDDO

    ! U = Y
    U=(Z-T1)/2
    ! Z = X
    Z=(Z+T1)/2

    !WRITE(*,*) 'X:'
    !WRITE(*,*) Z(:)
    !WRITE(*,*) 'Y:'
    !WRITE(*,*) U(:)
    ! test calculate X X+ - Y Y+  should be equal unity matrix
    ! CALL GPGEMMbse( 'N', trans, NCV, NCV, NCV, onebse, & 
    !   Z(1), 1, 1, DESC, Z(1), 1, 1, DESC, zerobse, T1(1), 1, 1, DESC)
    ! CALL GPGEMMbse( 'N', trans, NCV, NCV, NCV, -onebse, & 
    !   U(1), 1, 1, DESC, U(1), 1, 1, DESC, onebse, T1(1), 1, 1, DESC)
    ! DO I=1,NCV*NCV
    !    WRITE(77,'(10E20.10)') T1(I)
    ! ENDDO

    DEALLOCATE(T1,T2)

    IF (NSR>0  .AND. IU0>=0) WRITE(*,'(1X,A,I5)') 'WARNING in BSE_EIGENVECTORS_FULL: some eigenvalues in R are negative',NSR
    IF (NSRB>0 .AND. IU0>=0) WRITE(*,'(1X,A,I5)') 'WARNING in BSE_EIGENVECTORS_FULL: some eigenvalues in RB are negative',NSRB

      
  END SUBROUTINE BSE_EIGENVECTORS_FULL_SCA  
#endif


!****************** SUBROUTINE  AC_CORRE *****************************
!
! in this subroutine the density matrix P_alpha = |X,Y><X,Y| is 
! constructed and some orthogonality properties are tested
! The matrix is then contracted with the Coulomb superoperator
! to yield the correlation energy using the plasmon equation
!
!*********************************************************************

  SUBROUTINE AC_CORRE(NCV, X, Y, H, SUM_ALPHA)

    IMPLICIT NONE
    !passed variables
    INTEGER :: NCV
#ifdef single_prec_bse
    GDEFS :: X(:,:), Y(:,:), H(:,:)
    GDEFS,ALLOCATABLE :: T1(:,:),T2(:,:),T3(:,:),TMP(:,:),EYE(:,:)
#else
    GDEF  :: X(:,:), Y(:,:), H(:,:)
    GDEF,ALLOCATABLE :: T1(:,:),T2(:,:),T3(:,:),TMP(:,:),EYE(:,:)
#endif

    !local variables
    INTEGER :: I, J
    COMPLEX(q) :: SUM_ALPHA
 
    ALLOCATE(T1(NCV,NCV),EYE(NCV,NCV),T2(NCV,NCV),T3(NCV,NCV),TMP(NCV,NCV))
 
    EYE(:,:)=0.0_q
    DO I=1,NCV
       EYE(I,I) = 1.0_q
    ENDDO

    ! test calculate X X+ - Y Y+  should be equal unity matrix
    !  CALL GGEMMbse( 'N', trans, NCV, NCV, NCV, onebse, & 
    !    X(1,1), NCV, X(1,1), NCV, zerobse, T2(1,1), NCV)
    !  CALL GGEMMbse( 'N', trans, NCV, NCV, NCV, -onebse, & 
    !    Y(1,1), NCV, Y(1,1), NCV, onebse, T2(1,1), NCV)
    !  DO I=1,NCV
    !     WRITE(77,*) T2(:,I)
    !  ENDDO
    !  WRITE(77,*) '    '
    !  T2(:,:)=0.0_q
    !end test
   
    ! XX+ --> T1
      CALL GGEMMbse( 'N', trans, NCV, NCV, NCV, onebse, & 
        X(1,1), NCV, X(1,1), NCV, zerobse, T1(1,1), NCV)

    ! YY+ --> T2
      CALL GGEMMbse( 'N', trans, NCV, NCV, NCV, onebse, & 
        Y(1,1), NCV, Y(1,1), NCV, zerobse, T2(1,1), NCV)

    ! XY+ --> T3
      CALL GGEMMbse( 'N', trans, NCV, NCV, NCV, onebse, & 
        X(1,1), NCV, Y(1,1), NCV, zerobse, T3(1,1), NCV)

    ! test: according to Scuseria et al. (J. Chem. Phys. 2008, 129, 231101)
    ! Y+X = X+Y (where + indicates adjoint); check that XY+ = YX+
    !  CALL GGEMMbse( 'N', trans, NCV, NCV, NCV, onebse, & 
    !    Y(1,1), NCV, X(1,1), NCV, -onebse, T3(1,1), NCV)
    !  DO I=1,NCV
    !     WRITE(78,*) T3(:,I)
    !  ENDDO
    !  WRITE(78,*) '    '
    ! end test 

    ! T1 + T2 - EYE --> TMP; 
    ! T3 --> 2*T3 ;
    DO I=1,NCV
      DO J=1,NCV
        TMP(I,J) = T1(I,J) + T2(I,J) - EYE(I,J)
        T3(I,J) = 2.0_q*T3(I,J)
      ENDDO
    ENDDO

    T1(:,:)=0.0 ; T2(:,:)=0.0  

    ! H*TMP --> T1 
    ! H*T3 --> T2
    ! now contract T1 and T3 with the Hartree kernel:
     CALL GGEMMbse( 'N', 'N', NCV, NCV, NCV, onebse, & 
       H(1,1), NCV, TMP(1,1), NCV, zerobse, T1(1,1), NCV)
     CALL GGEMMbse( 'N', 'N', NCV, NCV, NCV, onebse, & 
       H(1,1), NCV, T3(1,1), NCV, zerobse, T2(1,1), NCV)

    ! and take the trace
     SUM_ALPHA=0.0_q
     DO I=1,NCV
        SUM_ALPHA=SUM_ALPHA+T1(I,I)+T2(I,I)
     ENDDO

    DEALLOCATE(T1,T2,T3,TMP,EYE)

  END SUBROUTINE AC_CORRE

#ifdef scaLAPACK

!**********************************************************************
!
! same as above but for scaLAPACK
!
!**********************************************************************

  SUBROUTINE AC_CORRE_SCA(NCV, X, Y, H, SUM_ALPHA, MY_COMM, DESC )

  USE wave_high

    INTEGER :: NCV
    COMPLEX(q) :: SUM_ALPHA, SUM_TMP
    TYPE(communic) :: MY_COMM ! communicator to be passed to the subroutine that computes the trace 
    
#ifdef single_prec_bse
    GDEFS :: X(:), Y(:), H(:)
    GDEFS,ALLOCATABLE :: T1(:),T2(:),T3(:),TMP(:)
#else
    GDEF  :: X(:), Y(:), H(:)
    GDEF,ALLOCATABLE :: T1(:),T2(:),T3(:),TMP(:)
#endif
    INTEGER :: I, J
    INTEGER :: DESC(DLEN_)
! BLACS variables
    INTEGER, EXTERNAL ::     NUMROC
    INTEGER MYROW, MYCOL, NPROW, NPCOL, NP, NQ
    INTEGER I1RES,J1RES,IROW,JCOL
    INTEGER I1,I2,J1,J2
    CALL BLACS_GRIDINFO(DESC(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

    NP = NUMROC(NCV,DESC(MB_),MYROW,0,NPROW)
    NQ = NUMROC(NCV,DESC(NB_),MYCOL,0,NPCOL)
 
    ALLOCATE(T1(SIZE(X)),T2(SIZE(Y)),T3(SIZE(X)),TMP(SIZE(Y)))!,IDIAG(SIZE(X)))


    ! XX+ --> T1
      CALL GPGEMMbse( 'N',trans , NCV, NCV, NCV, onebse, & 
         X(1), 1, 1, DESC, X(1), 1, 1, DESC, zerobse, T1(1), 1, 1, DESC)

    ! YY+ --> T2
      CALL GPGEMMbse( 'N',trans , NCV, NCV, NCV, onebse, & 
         Y(1), 1, 1, DESC, Y(1), 1, 1, DESC, zerobse, T2(1), 1, 1, DESC)

    ! test calculate X X+ - Y Y+ should be equal unity matrix
    ! CALL GPGEMMbse( 'N', trans, NCV, NCV, NCV, onebse, & 
    !     X(1), 1, 1, DESC, X(1), 1, 1, DESC, zerobse, T1(1), 1, 1, DESC)
    ! CALL GPGEMMbse( 'N', trans, NCV, NCV, NCV, -onebse, & 
    !     Y(1), 1, 1, DESC, Y(1), 1, 1, DESC, onebse, T1(1), 1, 1, DESC)
    !  JCOL=0
    !  DO J1=1,NQ,DESC(NB_)
    !     J1RES=MIN(DESC(NB_),NQ-J1+1)
    !     DO J2=1,J1RES
    !        JCOL=JCOL+1
    !        J=DESC(NB_)*MYCOL+NPCOL*(J1-1)+J2
    !        IROW=0
    !        DO I1=1,NP,DESC(MB_)
    !           I1RES=MIN(DESC(MB_),NP-I1+1)
    !           DO I2=1,I1RES
    !              IROW=IROW+1
    !              I=DESC(MB_)*MYROW+NPROW*(I1-1)+I2
    !              IF (I==J) THEN
    !               WRITE(77,*) T1(IROW+(JCOL-1)*DESC(LLD_))
    !              ENDIF
    !              WRITE(78,*) T1(IROW+(JCOL-1)*DESC(LLD_))
    !           ENDDO
    !        ENDDO
    !     ENDDO
    !  ENDDO
    ! end test (passed for LiH @ \Gamma on 04/11/14)

    ! XY+ --> T3
      CALL GPGEMMbse( 'N',trans , NCV, NCV, NCV, onebse, & 
         X(1), 1, 1, DESC, Y(1), 1, 1, DESC, zerobse, T3(1), 1, 1, DESC)

    ! Now construct the matrices:
    ! T1 + T2 - EYE --> TMP ; T3 --> 2*T3 

      JCOL=0
      DO J1=1,NQ,DESC(NB_)
         J1RES=MIN(DESC(NB_),NQ-J1+1)
         DO J2=1,J1RES
            JCOL=JCOL+1
            J=DESC(NB_)*MYCOL+NPCOL*(J1-1)+J2
            IROW=0
            DO I1=1,NP,DESC(MB_)
               I1RES=MIN(DESC(MB_),NP-I1+1)
               DO I2=1,I1RES
                  IROW=IROW+1
                  I=DESC(MB_)*MYROW+NPROW*(I1-1)+I2
                  TMP(IROW+(JCOL-1)*DESC(LLD_))=T1(IROW+(JCOL-1)*DESC(LLD_))+T2(IROW+(JCOL-1)*DESC(LLD_))
                  T3(IROW+(JCOL-1)*DESC(LLD_))=(T3(IROW+(JCOL-1)*DESC(LLD_)))*2.0_q
                  IF (I==J) THEN
                     TMP(IROW+(JCOL-1)*DESC(LLD_))=TMP(IROW+(JCOL-1)*DESC(LLD_))-1.0_q
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO


    ! V*TMP --> T1 ; W*T3 --> T2
      CALL GPGEMMbse( 'N','N' , NCV, NCV, NCV, onebse, & 
         H(1), 1, 1, DESC, TMP(1), 1, 1, DESC, zerobse, T1(1), 1, 1, DESC)

      CALL GPGEMMbse( 'N','N', NCV, NCV, NCV, onebse, & 
         H(1), 1, 1, DESC, T3(1), 1, 1, DESC, zerobse, T2(1), 1, 1, DESC)

    ! Now take the trace tr{T1}+tr{T2}
      SUM_TMP=0.0_q

      CALL TWOELECTRON4O_GET_DIAG_SCALA( T1, NCV, SUM_TMP, MY_COMM, DESC )
      SUM_ALPHA=SUM_TMP
      SUM_TMP=0.0_q
      CALL TWOELECTRON4O_GET_DIAG_SCALA( T2, NCV, SUM_TMP, MY_COMM, DESC )
      SUM_ALPHA=SUM_ALPHA+SUM_TMP

    DEALLOCATE(T1,T2,T3,TMP)     

  END SUBROUTINE AC_CORRE_SCA
#endif


!****************** SUBROUTINE SELFEN_2ndO *******************************
!
! in this subroutine the self-energy is evaluated according to the 
! expression provided by van Leeuwen in PRB (2006) 74, 195105
! also consider the diagram analysis in Radi Jishi "Feynman diagram techniques 
! in Cond Matt Phys" Cambridge (ch 9, page 199 & ff.)
! the idea is to implement something like
! \Sigma_c = \V GGG \V
! where \V is the bare Hartree contribution, 
! the propagators in the expression above are bare propagators and
! the energy denominators are obtained using the single particle 
! energies of the respective orbitals (no reference to the BSE eigenvalues)
! This is a Gamma point only version that does not support the data structure
! for accumulation over q-points
! 
!*********************************************************************

  SUBROUTINE SELFEN_2ndO(BSE_INDEX, POT_INDEX, BD, S2E, W, V, IRRVA, IRRVB, NQ, &
       LATT_CUR, AHARTREE, IU6, IU0, SEFROMBSE, LADDER)
  ! passed variables 
    USE constant
    TYPE (bse_matrix_index) :: BSE_INDEX, POT_INDEX
    TYPE (banddesc) :: BD(:)
    TYPE (screened_2e_handle) :: S2E
    TYPE (selfenergy_from_bse) :: SEFROMBSE
    TYPE (wavespin) :: W
    
#ifdef single_prec_bse
    GDEFS :: V(:,:), IRRVA(:,:), IRRVB(:,:)  
#else
    GDEF  :: V(:,:), IRRVA(:,:), IRRVB(:,:) 
#endif
    INTEGER :: NQ ! q-point 
    !COMPLEX(q) :: CELTOT(:,:,:) ! this contains the Hartree + kinetic + exchange term
    TYPE (latt) LATT_CUR
    REAL(q) :: AHARTREE
    INTEGER :: IU6 !output unit
    INTEGER :: IU0
    LOGICAL :: LADDER
    
    ! local variables
#ifdef single_prec_bse
    GDEFS, ALLOCATABLE :: TMP(:,:)
#else
    GDEF, ALLOCATABLE  :: TMP(:,:)
#endif
    COMPLEX(q),ALLOCATABLE :: SE_SPECTR(:,:,:,:)   ! and its spectral representation
    REAL(q) :: OMEGAMAX, OMEGAMIN, DOMEGA, DWEIGHT
    REAL(q) :: E1, E2, E3, E4           ! internal energies in the GF denominators
    REAL(q) :: ERES, ERES1, ERES2       ! and their combination 
    REAL(q) :: P1G, P2G, P3G, P4G       ! (1 - Fermi function weight for the state)
    REAL(q) :: P1L, P2L, P3L, P4L       ! Fermi function weight for the state 
    REAL(q) :: P_PLUS, P_MINUS     ! combination of weights consistent with whether we are summing over 
    ! 2 particle or 2 hole states
    INTEGER :: VBMAX, CBMIN, CBMAX
    INTEGER :: NOMEGA
    INTEGER :: I, J, NPIN ! counters of the global BSE matrix
    INTEGER :: NCV        ! number of e-h pairs in the BSE matrix
    INTEGER :: K1, K2, K3, K4 ! momentum counters
    INTEGER :: N1, N2, N3, N4 ! band counters
    INTEGER :: K1_IN_IRZ, K3_IN_IRZ ! momenta in the irreducible Brillouin zone
    INTEGER :: IV      ! counter over the potential matrix which includes also particle-particle integrals
    INTEGER :: ISP, ISP2, ISP_LOW, ISP_HIGH  !spin variables
    GDEF    :: CDER_BETWEEN_STATE(3), AHEAD
    INTEGER :: K2_IN_FULL_ORIG

    ! 0) allocation 
    NCV   = BSE_INDEX%NCV
    VBMAX = BD(1)%VBMAX ! BD(2) seems to be not set

    NOMEGA   = SEFROMBSE%NOMEGA
    OMEGAMAX = SEFROMBSE%OMEGA(NOMEGA)
    OMEGAMIN = SEFROMBSE%OMEGA(1)
    DOMEGA   =(OMEGAMAX-OMEGAMIN)/(NOMEGA-1)

    ! for the second order with exchange only (SOSEX) remove the Hartree contribution from IRRVB and IRRVA
    !IF (LADDER) THEN
    !   ALLOCATE(TMP(POT_INDEX%NCV,BSE_INDEX%NCV))
    !   TMP=0._q
    !   TMP=IRRVB-V
    !   IRRVB=TMP
    !   TMP=0._q
    !   TMP=IRRVA-V
    !   IRRVA=TMP
    !ELSE
    !   ALLOCATE(TMP(1,1))
    !ENDIF
    !DEALLOCATE(TMP)

    ALLOCATE( SE_SPECTR(NOMEGA, VBMAX, SEFROMBSE%NKPTS, SEFROMBSE%ISPIN))

    ! initialisation 
    SE_SPECTR=0._q

    IF (IU0>=0) THEN
       WRITE(IU0,'(A,2F7.2)') ' calculating selfenergy between w=',OMEGAMIN,OMEGAMAX
    ENDIF
!================================================================================
! Regular contribution to the spectral representation of the self-energy  
! this is made up by two contributions: resonant and antiresonant and the
! frequencies are pinned to the values given by the sum (resonant) or difference
! (antiresonant) between the plasmon's energy and the scattered state's (E3).
! For the resonant process the scattered state is a particle (unoccupied state)
! whereas for the antiresonant componentthe scattered state is a hole (occupied
! state). Since the incoming state |N1,K1,ISP> can either be occupied or
! unoccupied all possible self-energy diagrams are constructed 
!================================================================================
  ! 2) loop over first counter I --> |n1, k1, ISP >
      DO I = 1, POT_INDEX%NCV
         N1  = POT_INDEX%N1(I)
         K1  = POT_INDEX%NK(I)
         ISP = POT_INDEX%ISP(I)
         !E1  = W%CELTOT(N1, K1, ISP)
         E1  = SEFROMBSE%CELNEW(N1, K1, ISP)
         P1G = -(1.0_q-W%FERTOT(N1,K1,ISP))
         P1L = W%FERTOT(N1,K1,ISP)
  ! 2a) knowing NQ obtain |n3, k3, ISP >
         N3  = POT_INDEX%N3(I)
         K3  = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,K1)+W%WDES%VKPT(:,NQ),KPOINTS_FULL)
         
         ! bring K1 into IRZ, and apply the symmetry operation that brings K1 to K1_IN_IRZ 
         ! to K3 yielding K3_IN_IRZ 
         K1_IN_IRZ=S2E%K1_IN_IRZ(NQ, K1)
         IF (K1_IN_IRZ > SIZE(SE_SPECTR,3)) THEN
            CALL vtutor%bug("internal error in SELFEN_2ndO: K1_IN_IRZ exceeds SE_SPECTR " // &
               str(K1_IN_IRZ) // " " // str(SIZE(SE_SPECTR,3)), __FILE__, __LINE__)
         ENDIF

         K3_IN_IRZ=S2E%K2_IN_IRZ(NQ, K1)

         IV = POT_INDEX%INDEX(N1,N3,K1,ISP)
         ! if there is no potential matrix element calculated, CYCLE
         IF( IV==0) CYCLE
  
  ! 2b) obtain energy of |n3, k3, ISP > and shift wrt EFERMI
         !E3  = W%CELTOT(N3, K3, ISP)
         E3  = SEFROMBSE%CELNEW(N3, K3, ISP)
  ! 2c) construct appropriate weight for the spectral function P3=P(n3, k3, ISP)
         P3G = -(1.0_q-W%FERTOT(N3,K3,ISP))
         P3L = W%FERTOT(N3,K3,ISP)
         IF (LADDER) THEN
            DO J = 1, NCV
               N2  = BSE_INDEX%N1(J)
               K2  = BSE_INDEX%NK(J)
               ISP2= BSE_INDEX%ISP(J)
     ! 3a) construct spectral weight P2; shift energy wrt EFERMI
               !E2  = W%CELTOT(N2, K2, ISP2)
               E2  = SEFROMBSE%CELNEW(N2, K2, ISP2)
               P2G = -(1.0_q-W%FERTOT(N2,K2,ISP2))
               P2L = W%FERTOT(N2,K2,ISP2)
     ! 3b) knowing NQ obtain |n4, k4, ISP2 >
               N4 = BSE_INDEX%N3(J)
               K4 = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,K2)+W%WDES%VKPT(:,NQ),KPOINTS_FULL)
     ! 3c) construct spectral weight P4; shift energy wrt EFERMI
               !E4  = W%CELTOT(N4, K4, ISP2)
               E4  = SEFROMBSE%CELNEW(N4, K4, ISP2)
               P4G = -(1.0_q-W%FERTOT(N4,K4,ISP2))
               P4L = W%FERTOT(N4,K4,ISP2)
     ! 3d) construct the energy denominator and the occupancy product
               ! if state 3 is empty same sign as 4
               P_PLUS = P3G*P2L*P4G
               ERES1 =  E3-E2+E4
              ! if state 3 is occupied, same sign as 2
               P_MINUS =P3L*P2L*P4G
               ERES2 = -E3-E2+E4
   
               NPIN=INT((ERES1-OMEGAMIN)/((OMEGAMAX-OMEGAMIN)/REAL(NOMEGA-1))+1.5)
               NPIN=MIN(MAX(1,NPIN),NOMEGA)
               SE_SPECTR(NPIN,N1,K1_IN_IRZ,ISP)=SE_SPECTR(NPIN,N1,K1_IN_IRZ,ISP)- V(IV,J)*GCONJG(IRRVB(IV,J))*P_PLUS*S2E%WTKPT(K1_IN_IRZ, K3_IN_IRZ)
               NPIN=INT((-ERES2-OMEGAMIN)/((OMEGAMAX-OMEGAMIN)/REAL(NOMEGA-1))+1.5)
               NPIN=MIN(MAX(1,NPIN),NOMEGA)
               SE_SPECTR(NPIN,N1,K1_IN_IRZ,ISP)=SE_SPECTR(NPIN,N1,K1_IN_IRZ,ISP)+ V(IV,J)*GCONJG(IRRVA(IV,J))*P_MINUS*S2E%WTKPT(K1_IN_IRZ, K3_IN_IRZ)
            ENDDO
         ELSE
            DO J = 1, NCV
               N2  = BSE_INDEX%N1(J)
               K2  = BSE_INDEX%NK(J)
               ISP2= BSE_INDEX%ISP(J)
     ! 3a) construct spectral weight P2; shift energy wrt EFERMI
               !E2  = W%CELTOT(N2, K2, ISP2)
               E2  = SEFROMBSE%CELNEW(N2, K2, ISP2)
               P2G = -(1.0_q-W%FERTOT(N2,K2,ISP2))
               P2L = W%FERTOT(N2,K2,ISP2)
     ! 3b) knowing NQ obtain |n4, k4, ISP2 >
               N4 = BSE_INDEX%N3(J)
               K4 = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,K2)+W%WDES%VKPT(:,NQ),KPOINTS_FULL)
     ! 3c) construct spectral weight P4; shift energy wrt EFERMI
               !E4  = W%CELTOT(N4, K4, ISP2)
               E4  = SEFROMBSE%CELNEW(N4, K4, ISP2)
               P4G = -(1.0_q-W%FERTOT(N4,K4,ISP2))
               P4L = W%FERTOT(N4,K4,ISP2)
     ! 3d) construct the energy denominator and the occupancy product
               ! if state 3 is empty same sign as 4
               P_PLUS = P3G*P2L*P4G
               ERES1 =  E3-E2+E4
              ! if state 3 is occupied, same sign as 2
               P_MINUS =P3L*P2L*P4G
               ERES2 = -E3-E2+E4
   
               NPIN=INT((ERES1-OMEGAMIN)/((OMEGAMAX-OMEGAMIN)/REAL(NOMEGA-1))+1.5)
               NPIN=MIN(MAX(1,NPIN),NOMEGA)
               SE_SPECTR(NPIN,N1,K1_IN_IRZ,ISP)=SE_SPECTR(NPIN,N1,K1_IN_IRZ,ISP)- V(IV,J)*GCONJG(V(IV,J))*P_PLUS*S2E%WTKPT(K1_IN_IRZ, K3_IN_IRZ)
               NPIN=INT((-ERES2-OMEGAMIN)/((OMEGAMAX-OMEGAMIN)/REAL(NOMEGA-1))+1.5)
               NPIN=MIN(MAX(1,NPIN),NOMEGA)
               SE_SPECTR(NPIN,N1,K1_IN_IRZ,ISP)=SE_SPECTR(NPIN,N1,K1_IN_IRZ,ISP)+ V(IV,J)*GCONJG(V(IV,J))*P_MINUS*S2E%WTKPT(K1_IN_IRZ, K3_IN_IRZ)
            ENDDO ! over J 
         ENDIF
      ENDDO ! over I

!==========================================================================
! contribution from head
!==========================================================================
  ! 2) loop over first counter I --> |n1, k1, ISP >
      DO I = 1, POT_INDEX%NCV
         N1  = POT_INDEX%N1(I)
         K1  = POT_INDEX%NK(I)
         ISP = POT_INDEX%ISP(I)
         !E1  = W%CELTOT(N1, K1, ISP)
         E1  = SEFROMBSE%CELNEW(N1, K1, ISP)
         P1G = -(1.0_q-W%FERTOT(N1,K1,ISP))
         P1L = W%FERTOT(N1,K1,ISP)
  ! 2a) knowing NQ obtain |n3, k3, ISP >
         N3  = POT_INDEX%N3(I)
         K3  = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,K1)+W%WDES%VKPT(:,NQ),KPOINTS_FULL)
         
         IF (N1 /= N3 .OR. K1 /= K3) CYCLE
         ! we need contributions only inside IRZ, so short-cut to IRZ here
         IF (K1 > KPOINTS_ORIG%NKPTS) CYCLE

  ! 2b) obtain energy of |n3, k3, ISP > and shift wrt EFERMI
         !E3  = W%CELTOT(N3, K3, ISP)
         E3  = SEFROMBSE%CELNEW(N3, K3, ISP)
  ! 2c) construct appropriate weight for the spectral function P3=P(n3, k3, ISP)
         P3G = -(1.0_q-W%FERTOT(N3,K3,ISP))
         P3L = W%FERTOT(N3,K3,ISP)

  ! 3) loop over second global counter J --> |n2, k2, IPS2 > 
         DO J = 1, NCV
            N2  = BSE_INDEX%N1(J)
            K2  = BSE_INDEX%NK(J)
            ISP2= BSE_INDEX%ISP(J)
  ! 3a) construct spectral weight P2; shift energy wrt EFERMI
            !E2  = W%CELTOT(N2, K2, ISP2)
            E2  = SEFROMBSE%CELNEW(N2, K2, ISP2)
            P2G = -(1.0_q-W%FERTOT(N2,K2,ISP2))
            P2L = W%FERTOT(N2,K2,ISP2)
  ! 3b) knowing NQ obtain |n4, k4, ISP2 >
            N4 = BSE_INDEX%N3(J)
            K4 = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,K2)+W%WDES%VKPT(:,NQ),KPOINTS_FULL)

            IF (K2 /= K4) CYCLE
            K2_IN_FULL_ORIG=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,K2),KPOINTS_FULL_ORIG)
            CALL  CDER_BETWEEN_STATES_ROTATED_NOERR( &
                 CDER_BETWEEN_STATE,LATT_CUR, K2_IN_FULL_ORIG, ISP, N2 , N4 )

            AHEAD=(CDER_BETWEEN_STATE(1)*GCONJG(CDER_BETWEEN_STATE(1)))*(EDEPS/LATT_CUR%OMEGA)*W%WDES%RSPIN*W%WDES%RSPIN*W%WDES%WTKPT(K2)
            
  ! 3c) construct spectral weight P4; shift energy wrt EFERMI
            !E4  = W%CELTOT(N4, K4, ISP2)
            E4  = SEFROMBSE%CELNEW(N4, K4, ISP2)
            P4G = -(1.0_q-W%FERTOT(N4,K4,ISP2))
            P4L = W%FERTOT(N4,K4,ISP2)
  ! 3d) construct the energy denominator and the occupancy product
            ! if state 3 is empty same sign as 4
            P_PLUS = P3G*P2L*P4G
            ERES1 =  E3-E2+E4

           ! if state 3 is occupied, same sign as 2
            P_MINUS =P3L*P2L*P4G
            ERES2 = -E3-E2+E4

            ! right now the reason for W%WDES%NKPTS is not entirely obvious (most likely related to FSG)
            NPIN=INT((ERES1-OMEGAMIN)/((OMEGAMAX-OMEGAMIN)/REAL(NOMEGA-1))+1.5)
            NPIN=MIN(MAX(1,NPIN),NOMEGA)
            SE_SPECTR(NPIN,N1,K1,ISP)=SE_SPECTR(NPIN,N1,K1,ISP)-AHEAD*P_PLUS*(FSG_STORE(1)/W%WDES%NKPTS)
  
            NPIN=INT((-ERES2-OMEGAMIN)/((OMEGAMAX-OMEGAMIN)/REAL(NOMEGA-1))+1.5)
            NPIN=MIN(MAX(1,NPIN),NOMEGA)
            SE_SPECTR(NPIN,N1,K1,ISP)=SE_SPECTR(NPIN,N1,K1,ISP)+AHEAD*P_MINUS*(FSG_STORE(1)/W%WDES%NKPTS)
         ENDDO ! over J 
      ENDDO ! over I


  ! 10) print self-energy for each state |n1, k1, ISP >
  ! the factor 1/2 originates from the Hartree terms being multiplied by 2 in the stored
  ! two electron four orbital matrix elements
  ! this means we obtain 4 times the 2nd order direct diagram G (V P_0 V) (instead of twice)
  ! and twice the 2nd order direct diagram G (I P_0 V) (instead of once)
    IF (VBMAX>SEFROMBSE%NBANDS .AND. IU0>=0) THEN
       WRITE(IU0,'(A, 2I3)') 'SELFEN_2ndO: you are throwing away computed bands!', VBMAX, SEFROMBSE%NBANDS
    ENDIF 

    DO ISP=1,SEFROMBSE%ISPIN
       DO K1=1,SEFROMBSE%NKPTS
          DO N1=1,SEFROMBSE%NBANDS
             SEFROMBSE%REAL_PART(:,N1,K1,ISP) = REAL(SE_SPECTR(:,N1,K1,ISP),q)*(1.0_q/AHARTREE)
             SEFROMBSE%IMAG_PART(:,N1,K1,ISP) = 0._q
          ENDDO
       ENDDO
    ENDDO
    
    DEALLOCATE(SE_SPECTR)

 END SUBROUTINE SELFEN_2ndO

#ifdef scaLAPACK
!****************** SUBROUTINE SELFEN_2ndO_SCALA *******************************
!
! same as previous routine but for scaLAPACK
! 
!******************************************************************************

  SUBROUTINE SELFEN_2ndO_SCALA(BSE_INDEX, POT_INDEX, BD, S2E, W, EFERMI, VMAT, NQPOINT, &
                               CELTOT, AHARTREE, IU6, IU0, SE, DESC)

  ! passed variables 
    USE constant
    TYPE (bse_matrix_index) :: BSE_INDEX, POT_INDEX
    TYPE (banddesc) :: BD(:)
    TYPE (screened_2e_handle) :: S2E
    TYPE (selfenergy_from_bse) :: SE
    TYPE (wavespin) :: W
    
#ifdef single_prec_bse
    GDEFS :: VMAT(:) 
#else
    GDEF  :: VMAT(:)
#endif
    REAL(q) :: EFERMI
    INTEGER :: NQPOINT ! q-point 
    COMPLEX(q) :: CELTOT(:,:,:) ! this contains the Hartree + kinetic + exchange term
    REAL(q) :: AHARTREE
    INTEGER :: IU6 !output unit
    INTEGER :: IU0
    
    ! local variables
    COMPLEX(q),ALLOCATABLE :: SIGMA_SPECTRAL(:,:,:,:) ! self-energy for occupied states 
    COMPLEX(q),ALLOCATABLE :: CTMP(:) 
    COMPLEX(q) :: OMEGACPLX 
    REAL(q), ALLOCATABLE :: F(:,:)      ! output function with frequencies and self-energy components 
    REAL(q) :: OMEGAMAX, OMEGAMIN, DOMEGA, DWEIGHT
    REAL(q) :: E1, E2, E3, E4           ! internal energies in the GF denominators
    REAL(q) :: ERES1, ERES2             ! and their combination 
    REAL(q) :: P1G, P2G, P3G, P4G       ! (1 - Fermi function weight for the state)
    REAL(q) :: P1L, P2L, P3L, P4L       ! Fermi function weight for the state 
    REAL(q) :: P_PLUS, P_MINUS     ! combination of weights consistent with whether we are summing over 
    ! 2 particle or 2 hole states
    INTEGER :: VBMAX, CBMIN, CBMAX
    INTEGER :: I, J         ! counters of the global BSE matrix
    INTEGER :: IB           ! auxiliary counter
    INTEGER :: NCV, NAB     ! number of e-h pairs in the BSE matrix and POT matrix respectively
    INTEGER :: NCV13, NCV24 ! index of the electron-hole pair (or particle-particle pair)
    INTEGER :: K1, K2, K3, K4 ! momentum counters
    INTEGER :: N1, N2, N3, N4 ! band counters
    INTEGER :: NB1_INTO_TOT, NB2_INTO_TOT, NB3_INTO_TOT, NB4_INTO_TOT
    INTEGER :: IV, JV  ! counters over the potential matrix which includes also particle-particle integrals
    INTEGER :: ISP, ISP2   !spin variables
    INTEGER :: NOMEGA
    INTEGER :: OMG ! counter over frequencies 
    INTEGER :: NPIN_RES, NPIN_ANT ! resonant or antiresonant frequencies where selfenergy is computed

    ! BLACS variables
    INTEGER, EXTERNAL ::     NUMROC
    INTEGER MYROW, MYCOL, NPROW, NPCOL, NP, NQ
    INTEGER I1RES, J1RES, IROW, JCOL
    INTEGER I1, I2, J1, J2
    INTEGER IFAIL
    INTEGER :: DESC(DLEN_) ! descriptor set up for potential matrix distribution
    
  ! 0) allocation 
    NCV = BSE_INDEX%NCV
    NAB = POT_INDEX%NCV
    VBMAX = BD(1)%VBMAX ! BD(2) seems to be not set
    CBMIN = BD(1)%CBMIN
    CBMAX = BD(1)%CBMAX

    NOMEGA   = SE%NOMEGA
    OMEGAMAX = SE%OMEGA(NOMEGA)
    OMEGAMIN = SE%OMEGA(1)
    DOMEGA   =(OMEGAMAX-OMEGAMIN)/(NOMEGA-1)

    ALLOCATE(CTMP(NOMEGA))
    ALLOCATE(SIGMA_SPECTRAL(NOMEGA, VBMAX, SE%NKPTS, SE%ISPIN))

    ! initialisation 
    SIGMA_SPECTRAL=0._q

    IF (IU0>=0) THEN
       WRITE(IU0,'(A,2F7.2)') ' calculating selfenergy between w=',OMEGAMIN,OMEGAMAX
    ENDIF

    CALL BLACS_GRIDINFO(DESC(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
    NP = NUMROC(NAB,DESC(MB_),MYROW,0,NPROW)
    NQ = NUMROC(NAB,DESC(NB_),MYCOL,0,NPCOL)
    
  ! 2) loop over first counter J --> |n1, k1, ISP >
    DO J = 1,NAB
       J1=   ((J-DESC(NB_)*MYCOL-1)/(NPCOL*DESC(NB_)))*DESC(NB_)+1
       J2=MOD((J-DESC(NB_)*MYCOL-1) ,NPCOL*DESC(NB_))+1

       JCOL=J1+J2-1
       IF (J2 <= 0 .OR. J2 > DESC(NB_)) THEN
          JCOL  =-1
       ELSE
          IF (J /= DESC(NB_)*MYCOL+NPCOL*(J1-1)+J2) THEN
             CALL vtutor%bug("internal error in SELFEN_2ndO_SCALA: " // str(J) // " " // &
                str(DESC(NB_)) // " " // str(NPCOL*(J1-1)) // " " // str(J2), __FILE__, __LINE__)
          ENDIF
       ENDIF
       IF (JCOL<0) THEN
          CTMP=0._q
       ELSE
          
          NCV13=DESC(NB_)*MYCOL+NPCOL*(J1-1)+J2
          N1  = POT_INDEX%N1(NCV13)
          K1  = POT_INDEX%NK(NCV13)
          ISP = POT_INDEX%ISP(NCV13)
          E1  = W%CELTOT(N1, K1, ISP)
          P1G = -(1.0_q-W%FERTOT(N1,K1,ISP))
          P1L = W%FERTOT(N1,K1,ISP)
  ! 2a) knowing NQ obtain |n3, k3, ISP >
          N3  = POT_INDEX%N3(NCV13)
          K3  = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,K1)+W%WDES%VKPT(:,NQPOINT),KPOINTS_FULL)
       
          JV = POT_INDEX%INDEX(N1,N3,K1,ISP)
          IF(JV/=NCV13) THEN
            CALL vtutor%bug("internal error1 in SE2ndORD_SCALA: " // str(NCV13) // " " // str(JV) &
               // " " // str(DESC(NB_)) // " " // str(NPCOL*(J1-1)), __FILE__, __LINE__)
          ENDIF
       ! if there is no potential matrix element calculated, CYCLE
          IF( JV==0) CYCLE

  ! 2b) obtain energy of |n3, k3, ISP > and shift wrt EFERMI
          E3  = W%CELTOT(N3, K3, ISP)
  ! 2c) construct appropriate weight for the spectral function P3=P(n3, k3, ISP)
          P3G = -(1.0_q-W%FERTOT(N3,K3,ISP))
          P3L = W%FERTOT(N3,K3,ISP)
        
  ! 3) loop over second global counter I --> |n2, k2, IPS2 > 
          IROW=0
          CTMP=0._q
          !DO I1=1,NP,DESCSTD(MB_)
          !   I1RES=MIN(DESCSTD(MB_),NP-I1+1)
          DO I1=1,NP,DESC(MB_)
             I1RES=MIN(DESC(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW=IROW+1 !this way to access the column index might be incorrect
          !DO I = 1,NAB
          !   I1=   ((I-DESCSTD(MB_)*MYROW-1)/(NPROW*DESCSTD(MB_)))*DESCSTD(MB_)+1
          !   I2=MOD((I-DESCSTD(MB_)*MYROW-1) ,NPROW*DESCSTD(MB_))+1
          !   IROW=I1+I2-1

             !IF (I2 <= 0 .OR. I2 > DESCSTD(MB_)) THEN
             !   IROW  =-1
             !ELSE
             !   IF (I /= DESCSTD(MB_)*MYROW+NPROW*(I1-1)+I2) THEN
             !      WRITE(*,*) 'internal error2 in SELFEN_2ndO_SCALA:',I,DESCSTD(MB_)*MYROW, NPROW*(I1-1), I2
             !      STOP
             !   ENDIF
             !ENDIF
             !IF (IROW<0) THEN
             !   !CTMP=0._q
             !ELSE
              ! global column index NCV24
                NCV24=DESC(MB_)*MYROW+NPROW*(I1-1)+I2

                N2  = POT_INDEX%N1(NCV24)
                K2  = POT_INDEX%NK(NCV24)
                ISP2= POT_INDEX%ISP(NCV24)
  ! 3a) construct spectral weight P2; shift energy wrt EFERMI
                E2  = W%CELTOT(N2, K2, ISP2)
                P2G = -(1.0_q-W%FERTOT(N2,K2,ISP2))
                P2L = W%FERTOT(N2,K2,ISP2)
  ! 3b) knowing NQ obtain |n4, k4, ISP2 >
                N4 = POT_INDEX%N3(NCV24)
                K4 = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,K2)+W%WDES%VKPT(:,NQPOINT),KPOINTS_FULL)
                
                ! if the index does not identify a particle-hole pair then  CYCLE
                IV = POT_INDEX%INDEX(N2,N4,K2,ISP2)
                IB = BSE_INDEX%INDEX(N2,N4,K2,ISP2)
                IF(IV/=NCV24) THEN
                !WRITE(*,*) 'internal error3 in SE2ndORD_SCALA:',NCV24,IV,DESCSTD(MB_)*MYROW,NPROW*(I1-1)+I2
                CALL vtutor%bug("internal error3 in SE2ndORD_SCALA: " // str(NCV24) // " " // &
                   str(IV) // " " // str(DESC(MB_)) // " " // str(NPROW*(I1-1)) // "\n internal error3 &
                   &in SE2ndORD_SCALA: " // str(NCV24) // " " // str(IV) // " " // str(N2) // " " // &
                   str(N4) // " " // str(K2) // " " // str(ISP2) // " " // str(IB), __FILE__, __LINE__)
                ENDIF
                ! if the index does not identify a particle-hole pair then  CYCLE
                IF(NCV24==0) CYCLE 
                IF( IB==0) CYCLE

  ! 3c) construct spectral weight P4; shift energy wrt EFERMI
                E4  = W%CELTOT(N4, K4, ISP2)
                P4G = -(1.0_q-W%FERTOT(N4,K4,ISP2))
                P4L = W%FERTOT(N4,K4,ISP2)
  ! 3d) construct the energy denominator and the occupancy product
                ! if state 3 is empty same sign as 4
                P_PLUS = P3G*P2L*P4G
                ERES1 =  E3-E2+E4

                NPIN_RES=INT((ERES1-OMEGAMIN)/((OMEGAMAX-OMEGAMIN)/REAL(NOMEGA-1))+1.5)
                NPIN_RES=MIN(MAX(1,NPIN_RES),NOMEGA)
                CTMP(NPIN_RES)=CTMP(NPIN_RES)-VMAT(IROW+(JCOL-1)*DESC(LLD_))*GCONJG(VMAT(IROW+(JCOL-1)*DESC(LLD_)))*P_PLUS

                IF((REAL(CTMP(NPIN_RES)) /= REAL(CTMP(NPIN_RES)) ) .OR. ( AIMAG(CTMP(NPIN_RES)) /= AIMAG(CTMP(NPIN_RES))) ) THEN 
                 !WRITE(*,*) 'NaN @ RES',NPIN_RES, N1,N3,N2,N4,NCV24,IB, DESCSTD(MB_)*MYROW, NPROW*(I1-1),I2
                  CALL vtutor%error("NaN @ RES " // str(NPIN_RES) // " " // str(N1) // " " // str(N3) &
                     // " " // str(N2) // " " // str(N4) // " " // str(NCV24) // " " // str(IB) // "&
                     & " // str(DESC(MB_)) // " " // str(NPROW*(I1-1)) // " " // str(I2))
                ENDIF
                !IF(I==1) WRITE(400+NQPOINT,'(2I3,2F8.3,2I3,2F8.3)') N2, K2, P2G, P2L, N4, K4, P4G, P4L 
                ! if state 3 is occupied, same sign as 2
                P_MINUS =P3L*P2L*P4G
                ERES2 = -E3-E2+E4

                NPIN_ANT=INT((-ERES2-OMEGAMIN)/((OMEGAMAX-OMEGAMIN)/REAL(NOMEGA-1))+1.5)
                NPIN_ANT=MIN(MAX(1,NPIN_ANT),NOMEGA)
                CTMP(NPIN_ANT)=CTMP(NPIN_ANT)+VMAT(IROW+(JCOL-1)*DESC(LLD_))*GCONJG(VMAT(IROW+(JCOL-1)*DESC(LLD_)))*P_MINUS

                IF(REAL(CTMP(NPIN_ANT)) /= REAL(CTMP(NPIN_ANT))) THEN 
                 !WRITE(*,*) 'NaN @ ANT',NPIN_ANT, N1,N3,N2,N4,NCV24,IB, DESCSTD(MB_)*MYROW, NPROW*(I1-1),I2
                  CALL vtutor%error("NaN @ ANT " // str(NPIN_ANT) // " " // str(N1) // " " // str(N3) &
                     // " " // str(N2) // " " // str(N4) // " " // str(NCV24) // " " // str(IB) // "&
                     & " // str(DESC(MB_)) // " " // str(NPROW*(I1-1)) // " " // str(I2))
                ENDIF

             ENDDO ! over I2
             !ENDIF ! IROW<0
          ENDDO !over I1 
          !ENDDO !over I 
       ! now sum CTMP over all nodes for the current J
       CALLMPI( M_sum_z(W%WDES%COMM_INTER, CTMP, NOMEGA))
       SIGMA_SPECTRAL(:,N1,K1,ISP)=SIGMA_SPECTRAL(:,N1,K1,ISP)+CTMP(:)
       ENDIF ! JCOL not present on node
    ENDDO ! over J
  ! 10) print self-energy for each state |n1, k1, ISP >
    IF (VBMAX>SE%NBANDS .AND. IU0>=0) THEN
       WRITE(IU0,'(A, 2I3)') 'SELFEN_2ndO_SCALA: you are throwing away computed bands!', VBMAX, SE%NBANDS
    ENDIF
  ! the factor 1/2 originates from the Hartree terms being multiplied by 2 in the stored
  ! two electron four orbital matrix elements
  ! this means we obtain 4 times the 2nd order direct diagram G (V P_0 V) (instead of twice)
  ! and twice the 2nd order direct diagram G (I P_0 V) (instead of once)
    DO ISP=1,SE%ISPIN
       DO K1=1,SE%NKPTS
          DO N1=1,SE%NBANDS
             SE%REAL_PART(:,N1,K1,ISP)=(REAL(SIGMA_SPECTRAL(:,N1,K1,ISP),q))*(1.0_q/AHARTREE)
             SE%IMAG_PART(:,N1,K1,ISP)= 0.0_q
           ENDDO
       ENDDO
    ENDDO
    
    DEALLOCATE(CTMP)
    DEALLOCATE(SIGMA_SPECTRAL)
     
 END SUBROUTINE SELFEN_2ndO_SCALA
#endif

!****************** SUBROUTINE SELFEN_FULL ************************************
!
! in this subroutine the self-energy is evaluated using an expression
! given first by Chelikowsky in PRB (2006) 73, 205334 with the difference
! that therein the eigenvectors used were from the squared RPA EVP
!
! Details can be found in
! Maggio and Kresse J. Chem. Theory Comput. 13, 4765-4778 (2017)
! Equ. (22) and erratum
! 
!*******************************************************************************

 SUBROUTINE SELFEN_FULL(BSE_INDEX, POT_INDEX, BD, S2E, W, X, Y, V, IRRVA, IRRVB, R, NQ, &
                            LATT_CUR, AHARTREE, IU6, IU0, SEFROMBSE, LADDER, LTCTC, LRESTORE )

  ! passed variables 
    USE constant
    USE dfast
    USE ini
    IMPLICIT NONE
    TYPE (bse_matrix_index) :: BSE_INDEX, POT_INDEX
    TYPE (banddesc) :: BD(:)
    TYPE (screened_2e_handle) :: S2E
    TYPE (selfenergy_from_bse) :: SEFROMBSE
    TYPE (wavespin) :: W
  
#ifdef single_prec_bse
    GDEFS :: X(:,:)      ! eigenvectors of BSE contributions X
    GDEFS :: Y(:,:)      ! eigenvectors of BSE contributions Y
    GDEFS :: V(:,:)      ! Coulomb kernel spans also particle-particle manyfold
    GDEFS :: IRRVA(:,:)  ! A matrix of BSE spans also particle-particle manyfold
    GDEFS :: IRRVB(:,:)  ! B matrix of BSE spans also particle-particle manyfold
#else
    GDEF :: X(:,:)
    GDEF :: Y(:,:)
    GDEF :: V(:,:)
    GDEF :: IRRVA(:,:)
    GDEF :: IRRVB(:,:)
#endif
    REAL(q) :: R(:)      ! eigenvalues of the BSE equation (for the spectral representation of L)
    INTEGER :: NQ        ! q-point 
    TYPE (latt) LATT_CUR 
    INTEGER :: IU6  ! output unit
    INTEGER :: IU0
    REAL(q) :: AHARTREE
    LOGICAL :: LADDER, LTCTC
    LOGICAL, OPTIONAL :: LRESTORE ! restore V, IRRVA, IRRVB upon exit

  ! local variables
    REAL(q),ALLOCATABLE :: SE_SPECTR(:,:,:,:)   ! and its spectral representation
    REAL(q) :: OMEGAMAX, OMEGAMIN, DOMEGA, DWEIGHT
    REAL(q) :: E3                 ! single particle energy in the GF denominator 
    REAL(q) :: ERES, EANT         ! and its combination with the plasmon energy 
    REAL(q) :: P3G                ! (1 - Fermi function weight for the state)
    REAL(q) :: P3L                ! Fermi function weight for the state 
    REAL(q) :: P_RES, P_ANT       ! combination of weights consistent with whether we are summing over 
                                  ! 2 particle or 2 hole states
    INTEGER :: VBMAX, NBANDS
   
    INTEGER :: NOMEGA
    INTEGER :: I, J     ! counters of the global matrix
    INTEGER :: K        ! auxiliary counter for the spectral representation of the polarisation propagator
    INTEGER :: NPINRES, NPINANT    ! counter for the spectral density
    COMPLEX(q) :: CRES1, CRES2 , CANT1, CANT2  ! tmp work storage
    INTEGER :: IV, JV   ! counters over the potential matrix which includes also p-p integrals
    INTEGER :: NCV      ! number of e-h pairs in the BSE matrix
    INTEGER :: K1, K2, K3, K4  ! momentum counters
    INTEGER :: K1_IN_IRZ, K3_IN_IRZ ! momenta in the irred Brillouin zone
    INTEGER :: N1, N2, N3, N4  ! band counters
    INTEGER :: ISP, ISP2       ! spin variables
    INTEGER :: ILEN, IBLOCK    ! for strip-wise matrix multiplication
    GDEF, ALLOCATABLE :: AHEAD(:), TMPHEAD(:)
    GDEF    :: CDER_BETWEEN_STATE(3)
    INTEGER :: K2_IN_FULL_ORIG
!$  INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM

    ! pre-contraction work arrays
#ifdef single_prec_bse
    GDEFS, ALLOCATABLE :: S(:,:)
#else
    GDEF, ALLOCATABLE :: S(:,:)
#endif
    ! pre-contraction result arrays
#ifdef single_prec_bse
    GDEFS, ALLOCATABLE :: VR(:,:), VN(:,:)
#else
    GDEF, ALLOCATABLE ::  VR(:,:), VN(:,:)
#endif

    ! pre-contraction save arrays, needed only for the gamma-only version where
    ! the kernel cannot be overwritten 
#ifdef single_prec_bse
    GDEFS, ALLOCATABLE :: V_SAVE(:,:), IRRVA_SAVE(:,:), IRRVB_SAVE(:,:)
#else
    GDEF, ALLOCATABLE :: V_SAVE(:,:), IRRVA_SAVE(:,:), IRRVB_SAVE(:,:)
#endif

    NCV   = BSE_INDEX%NCV
    VBMAX = BD(1)%VBMAX ! BD(2) seems to be not set

    NOMEGA   = SEFROMBSE%NOMEGA
    OMEGAMAX = SEFROMBSE%OMEGA(NOMEGA)
    OMEGAMIN = SEFROMBSE%OMEGA(1)
    DOMEGA   =(OMEGAMAX-OMEGAMIN)/(NOMEGA-1)

    IF (IU0>=0) THEN
       WRITE(IU0,'(A,2F7.2)') ' calculating selfenergy between w=',OMEGAMIN,OMEGAMAX
    ENDIF

    ! preliminary sanity tests
    IF (VBMAX>SEFROMBSE%NBANDS .AND. IU0>=0) THEN
       CALL vtutor%error("internal warning in SELFEN_FULL: you are throwing away computed bands! " // &
          str(VBMAX) // " " // str(SEFROMBSE%NBANDS))
    ENDIF

    IF ( .NOT. ASSOCIATED(SEFROMBSE%REAL_PART)) THEN
       CALL vtutor%bug("internal error in SELF_FULL: SELFENERGY data structure not allocated", __FILE__, __LINE__)
    ENDIF

    ! 0) allocation 
    ALLOCATE(SE_SPECTR(NOMEGA, VBMAX, SEFROMBSE%NKPTS, SEFROMBSE%ISPIN))

!============================================================================
! Pre-contraction of the BSE eigenvectors X and Y
! against the Coulomb kernel or the full irreducible kernel (+ xc).
! Since the kernel here spans also the particle/particle manifold beyond the
! usual e-/h manifold, to efficiently contract the kernel against the
! polarisation propagator the former is down-sized back to the original e-/h
! manifold along one dimension.  
!=============================================================================

! the kernel is needed if the self-consistency is performed inside the bse
! module, make a copy here before it gets overwritten
    IF (PRESENT(LRESTORE)) THEN
       ALLOCATE(V_SAVE(POT_INDEX%NCV, NCV),IRRVA_SAVE(POT_INDEX%NCV, NCV),IRRVB_SAVE(POT_INDEX%NCV, NCV))
       V_SAVE(:,:)     = V(:,:)
       IRRVA_SAVE(:,:) = IRRVA(:,:)
       IRRVB_SAVE(:,:) = IRRVB(:,:)
    ENDIF

#ifdef unblocked_version
    ! this is a simpler version that does not use blocking (required more workspace)
    ALLOCATE(S(NCV,NCV)) ! work arrays
    ALLOCATE(VR(POT_INDEX%NCV, NCV))
    IF (LADDER) ALLOCATE(VN(POT_INDEX%NCV, NCV))
    S=X+Y
    CALL GGEMMbse( 'N', 'N', POT_INDEX%NCV, NCV, NCV, onebse, &
       V(1,1), POT_INDEX%NCV, S(1,1), NCV, zerobse, VR(1,1), POT_INDEX%NCV)

    !overwrite the Hartree potential V with the product V (X+Y)
    V=VR

    IF (LADDER .AND. .NOT. LTCTC) THEN
       ! VR = IRRVA*X + IRRVB*Y --> IRRVA
       ! VN = IRRVA*X + IRRVB*Y --> IRRVB
       CALL GGEMMbse( 'N', 'N', POT_INDEX%NCV, NCV, NCV, onebse, &
          IRRVA(1,1), POT_INDEX%NCV, X(1,1), NCV, zerobse, VR(1,1), POT_INDEX%NCV)
       CALL GGEMMbse( 'N', 'N', POT_INDEX%NCV, NCV, NCV, onebse, &
          IRRVA(1,1), POT_INDEX%NCV, Y(1,1), NCV, zerobse, VN(1,1), POT_INDEX%NCV)
       CALL GGEMMbse( 'N', 'N', POT_INDEX%NCV, NCV, NCV, onebse, &
          IRRVB(1,1), POT_INDEX%NCV, X(1,1), NCV, onebse, VN(1,1), POT_INDEX%NCV)
       CALL GGEMMbse( 'N', 'N', POT_INDEX%NCV, NCV, NCV, onebse, &
          IRRVB(1,1), POT_INDEX%NCV, Y(1,1), NCV, onebse, VR(1,1), POT_INDEX%NCV)
    !overwrite the Irreducible vertex 
       IRRVA=VR
       IRRVB=VN
    ENDIF

    DEALLOCATE(S)
    DEALLOCATE(VR)
    IF (LADDER) DEALLOCATE(VN)
#else
    ALLOCATE(VR(NBLK, NCV))
    IF (LADDER) ALLOCATE(VN(NBLK, NCV))

    ! to save storage work on stripes of V and other matrices
    DO IBLOCK=0,POT_INDEX%NCV-1, NBLK 
      ILEN=MIN(NBLK, POT_INDEX%NCV-IBLOCK)

      ! set V_lambda = V (X_lambda+Y_lambda)
      CALL GGEMMbse( 'N', 'N', ILEN, NCV, NCV, onebse, &
         V(IBLOCK+1,1), POT_INDEX%NCV, X(1,1), NCV, zerobse, VR(1,1), ILEN )
     
      CALL GGEMMbse( 'N', 'N', ILEN, NCV, NCV, onebse, &
         V(IBLOCK+1,1), POT_INDEX%NCV, Y(1,1), NCV, onebse, VR(1,1), ILEN )
     
      V(IBLOCK+1:IBLOCK+ILEN,:)=VR(1:ILEN,:)
     
      IF (LADDER .AND. .NOT. LTCTC) THEN
         ! VR = IRRVA*X + IRRVB*Y --> IRRVA
         ! VN = IRRVA*X + IRRVB*Y --> IRRVB
         CALL GGEMMbse( 'N', 'N', ILEN, NCV, NCV, onebse, &
            IRRVA(IBLOCK+1,1), POT_INDEX%NCV, X(1,1), NCV, zerobse, VR(1,1), ILEN )
     
         CALL GGEMMbse( 'N', 'N', ILEN, NCV, NCV, onebse, &
            IRRVA(IBLOCK+1,1), POT_INDEX%NCV, Y(1,1), NCV, zerobse, VN(1,1), ILEN )
     
         CALL GGEMMbse( 'N', 'N', ILEN, NCV, NCV, onebse, &
            IRRVB(IBLOCK+1,1), POT_INDEX%NCV, X(1,1), NCV, onebse, VN(1,1), ILEN )
     
         CALL GGEMMbse( 'N', 'N', ILEN, NCV, NCV, onebse, &
            IRRVB(IBLOCK+1,1), POT_INDEX%NCV, Y(1,1), NCV, onebse, VR(1,1), ILEN )
      !overwrite the irreducible vertex 
         IRRVA(IBLOCK+1:IBLOCK+ILEN,:)=VR(1:ILEN,:)
         IRRVB(IBLOCK+1:IBLOCK+ILEN,:)=VN(1:ILEN,:)
      ENDIF
    ENDDO

    DEALLOCATE(VR)
    IF (LADDER) DEALLOCATE(VN)
#endif
!==========================================================================
! setup optical matrix elements, when the difference vector q is 0
!==========================================================================
    IF (  W%WDES%VKPT(1,NQ)**2+W%WDES%VKPT(2,NQ)**2+W%WDES%VKPT(3,NQ)**2<1E-10_q) THEN
       ALLOCATE(AHEAD(NCV),TMPHEAD(NCV))
       AHEAD=0
       DO J = 1, NCV
          N2  = BSE_INDEX%N1(J)
          K2  = BSE_INDEX%NK(J)
          ISP2= BSE_INDEX%ISP(J)
          N4 = BSE_INDEX%N3(J)
          K4 = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,K2)+W%WDES%VKPT(:,NQ),KPOINTS_FULL)

          K2_IN_FULL_ORIG=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,K2),KPOINTS_FULL_ORIG)
          CALL  CDER_BETWEEN_STATES_ROTATED_NOERR( &
               CDER_BETWEEN_STATE,LATT_CUR, K2_IN_FULL_ORIG, ISP2, N2 , N4 )
          ! get  <u_m|  -i d/dk_j u_n> = - <u_m| r_j |u_n>
          AHEAD(J)=CDER_BETWEEN_STATE(1)*SQRT((EDEPS/LATT_CUR%OMEGA)*W%WDES%RSPIN*W%WDES%RSPIN*W%WDES%WTKPT(K2))
       ENDDO
       ! transform transition matrix elements by X_lambda+Y_lambda
       ! AHEAD_lambda =  AHEAD (X_lambda + Y_lambda)
       TMPHEAD=0
       DO K = 1, NCV
          DO J = 1, NCV
             TMPHEAD(K)=TMPHEAD(K)+AHEAD(J)*(X(J,K)+Y(J,K))
          ENDDO
       ENDDO
       AHEAD=TMPHEAD
       DEALLOCATE(TMPHEAD)
    ENDIF

    ! initialisation 
    SE_SPECTR=0._q
!================================================================================
! Regular contribution to the spectral representation of the self-energy  
! this is made up by two contributions: resonant and antiresonant and the
! frequencies are pinned to the values given by the sum (resonant) or difference
! (antiresonant) between the plasmon's energy and the scattered state's (E3).
! For the resonant process the scattered state is a particle (unoccupied state)
! whereas for the antiresonant component the scattered state is a hole (occupied
! state). Since the incoming state |N1,K1,ISP> can either be occupied or
! unoccupied all possible self-energy diagrams are constructed 
!================================================================================
  ! 2) loop over first counter I --> |n1, k1, ISP >

!$OMP PARALLEL DO ORDERED PRIVATE(I,N1,K1,ISP,N3,K3,E3,P3G,P3L,P_RES,P_ANT,K1_IN_IRZ,K3_IN_IRZ,IV,K,ERES,NPINRES,EANT,NPINANT,CRES1, CRES2, CANT1, CANT2 )
    DO I = 1, POT_INDEX%NCV
       N1  = POT_INDEX%N1(I)
!gKtest: change N1-1 to I-1, should result in better load balancing 
!!$    IF (MOD(I-1,OMP_GET_NUM_THREADS())/=OMP_GET_THREAD_NUM()) CYCLE
#ifdef MPI
       IF (MOD(I-1,W%WDES%COMM_INTER%NCPU)+1/=W%WDES%COMM_INTER%NODE_ME) CYCLE
#endif
       K1  = POT_INDEX%NK(I)
       ISP = POT_INDEX%ISP(I)
       ! skip those that will have zero weights below (WTKPT)
       IF (.NOT. REQUIRED_SCREENED_2E(S2E, NQ, K1)) CYCLE

  ! 2a) knowing NQ obtain |n3, k3, ISP >
       N3  = POT_INDEX%N3(I)
       K3  = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,K1)+W%WDES%VKPT(:,NQ),KPOINTS_FULL)
  ! 2b) obtain energy of |n3, k3, ISP > and shift wrt EFERMI
       !E3  = W%CELTOT(N3, K3, ISP)
       E3  = SEFROMBSE%CELNEW(N3, K3, ISP)
  ! 2c) construct appropriate weight for the spectral function P3=P(n3, k3, ISP)
       P3G = -(1.0_q-W%FERTOT(N3,K3,ISP))
       P3L = W%FERTOT(N3,K3,ISP)
       !  hole in occupied orbital, scatters into another hole emitting a  plasmon
       P_RES =P3G
       !  hole in occupied orbital annihilates against electron adsorbing a plasmon
       P_ANT =P3L

       ! bring K1 into IRZ, and apply the symmetry operation that brings K1 to K1_IN_IRZ 
       ! to K3 yielding K3_IN_IRZ 
       K1_IN_IRZ=S2E%K1_IN_IRZ(NQ, K1)
       IF (K1_IN_IRZ > SIZE(SE_SPECTR,3)) THEN
          CALL vtutor%bug("internal error in SELFEN_FULL: K1_IN_IRZ exceeds SE_SPECTR " // &
             str(K1_IN_IRZ) // " " // str(SIZE(SE_SPECTR,3)), __FILE__, __LINE__)
       ENDIF
       K3_IN_IRZ=S2E%K2_IN_IRZ(NQ, K1)

       IV = POT_INDEX%INDEX(N1,N3,K1,ISP)
       ! if there is no potential matrix element calculated, CYCLE
       IF( IV/=I )  THEN
         CALL vtutor%bug("internal error in SELFEN_FULL: potential matrix element index do not match " &
            // str(I) // " " // str(IV), __FILE__, __LINE__)
       ENDIF

  ! 3)  loop over eigenvalue/ eigenvector pairs of BSE matrix 
!$OMP ORDERED
       DO K= 1, NCV
  ! 3a) construct the energy denominator and the occupancy product

          CRES1=0 ; CANT1=0 ; CRES2=0 ; CANT2=0

          IF (LTCTC) THEN
             CRES1=CRES1+V(IV,K)
             CRES2=CRES2+V(IV,K)
             ! now anti-resonant part (negative frequency components in L)
             CANT1=CANT1+V(IV,K)
             CANT2=CANT2+V(IV,K)
          ELSEIF (LADDER) THEN
       ! IRRVA*A + IRRVB*B --> IRRVA
       ! IRRVA*B + IRRVB*A --> IRRVB
             CRES1=CRES1+V(IV,K)
             CRES2=CRES2+IRRVB(IV,K)
             ! now anti-resonant part (negative frequency components in L)
             CANT1=CANT1+V(IV,K)
             CANT2=CANT2+IRRVA(IV,K)
          ELSE ! no Exchange in kernel
             CRES1=V(IV,K)
             CRES2=V(IV,K)
             CANT1=V(IV,K)
             CANT2=V(IV,K)
          ENDIF

          ERES   = E3+ R(K)
          NPINRES= INT((ERES-OMEGAMIN)/DOMEGA)+1
          NPINRES= MIN(MAX(1,NPINRES),NOMEGA-1)
          DWEIGHT= (ERES-OMEGAMIN-(NPINRES-1)*DOMEGA)/DOMEGA

          SE_SPECTR(NPINRES,N1,K1_IN_IRZ,ISP)=SE_SPECTR(NPINRES,N1,K1_IN_IRZ,ISP)+&
                (CRES1*CONJG(CRES2)*P_RES*S2E%WTKPT(K1_IN_IRZ, K3_IN_IRZ))*(1-DWEIGHT)
          SE_SPECTR(NPINRES+1,N1,K1_IN_IRZ,ISP)=SE_SPECTR(NPINRES+1,N1,K1_IN_IRZ,ISP)+&
                (CRES1*CONJG(CRES2)*P_RES*S2E%WTKPT(K1_IN_IRZ, K3_IN_IRZ))*DWEIGHT

          EANT   = -E3+ R(K)
          NPINANT= INT((-EANT-OMEGAMIN)/DOMEGA)+1
          NPINANT= MIN(MAX(1,NPINANT),NOMEGA-1)
          DWEIGHT= (-EANT-OMEGAMIN-(NPINANT-1)*DOMEGA)/DOMEGA

          SE_SPECTR(NPINANT,N1,K1_IN_IRZ,ISP)=SE_SPECTR(NPINANT,N1,K1_IN_IRZ,ISP)-&
                CANT1*CONJG(CANT2)*P_ANT*S2E%WTKPT(K1_IN_IRZ, K3_IN_IRZ)*(1-DWEIGHT)
          SE_SPECTR(NPINANT+1,N1,K1_IN_IRZ,ISP)=SE_SPECTR(NPINANT+1,N1,K1_IN_IRZ,ISP)-&
                CANT1*CONJG(CANT2)*P_ANT*S2E%WTKPT(K1_IN_IRZ, K3_IN_IRZ)*DWEIGHT
       ENDDO ! over K
!$OMP END ORDERED
    ENDDO ! over I
!$OMP END PARALLEL DO

! once that the modified kernels have been used, overwrite them with the
! original matrices
    IF (PRESENT(LRESTORE)) THEN
       V(:,:)     = V_SAVE(:,:)      
       IRRVA(:,:) = IRRVA_SAVE(:,:)  
       IRRVB(:,:) = IRRVB_SAVE(:,:)  
       DEALLOCATE(V_SAVE,IRRVA_SAVE,IRRVB_SAVE)
    ENDIF
!==========================================================================
! contribution from head
! head contribution is idential for each band but shifted by the eigenenergy
! (essentially proportional to FSG x frequency dependent dielectric function)
!==========================================================================
  ! 2) loop over first counter I --> |n1, k1, ISP >
    IF (ALLOCATED(AHEAD)) THEN
      DO I = 1, POT_INDEX%NCV
#ifdef MPI
       IF (MOD(I-1,W%WDES%COMM_INTER%NCPU)+1/=W%WDES%COMM_INTER%NODE_ME) CYCLE
#endif
         N1  = POT_INDEX%N1(I)
         K1  = POT_INDEX%NK(I)
         ISP = POT_INDEX%ISP(I)
         ! skip those that will have zero weights below (WTKPT)
         IF (.NOT. REQUIRED_SCREENED_2E(S2E, NQ, K1)) CYCLE
  
    ! 2a) knowing NQ obtain |n3, k3, ISP >
         N3  = POT_INDEX%N3(I)
         K3  = KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,K1)+W%WDES%VKPT(:,NQ),KPOINTS_FULL)
         ! intermediate state is the same as initial state, this yields a singular 
         ! contribution in Coulomb integral that needs to be considered
         IF (N1 /= N3 .OR. K1 /= K3) CYCLE
         ! we need contributions only inside IRZ, so short-cut to IRZ here
         IF (K1 > KPOINTS_ORIG%NKPTS) CYCLE
  
    ! 2b) obtain energy of |n3, k3, ISP > 
         !E3  = W%CELTOT(N3, K3, ISP)
         E3  = SEFROMBSE%CELNEW(N3, K3, ISP)
    ! 2c) construct appropriate weight for the spectral function P3=P(n3, k3, ISP)
         P3G = -(1.0_q-W%FERTOT(N3,K3,ISP))
         P3L = W%FERTOT(N3,K3,ISP)
         !  hole in occupied orbital, scatters into another hole emiting a  plasmon
         P_RES =P3G
         !  hole in occupied orbital annihilates against electron adsorbing a plasmon
         P_ANT =P3L
  
    ! 3)  loop over eigenvalue/ eigenvector pairs of BSE matrix 
         DO K= 1, NCV
    ! 3a) construct the energy denominator and the occupancy product
            CRES1=AHEAD(K)

           ! if state 3 is empty same sign as 4
            ERES   = E3+ R(K)
            NPINRES= INT((ERES-OMEGAMIN)/DOMEGA)+1
            NPINRES= MIN(MAX(1,NPINRES),NOMEGA-1)
            DWEIGHT= (ERES-OMEGAMIN-(NPINRES-1)*DOMEGA)/DOMEGA

            SE_SPECTR(NPINRES,N1,K1,ISP)=SE_SPECTR(NPINRES,N1,K1,ISP)+&
                  CRES1*CONJG(CRES1)*P_RES*(FSG_STORE(1)/W%WDES%NKPTS)*(1-DWEIGHT)
            SE_SPECTR(NPINRES+1,N1,K1,ISP)=SE_SPECTR(NPINRES+1,N1,K1,ISP)+&
                  CRES1*CONJG(CRES1)*P_RES*(FSG_STORE(1)/W%WDES%NKPTS)*DWEIGHT

            ! if state 3 is occupied, same sign as 2
            EANT   = -E3+ R(K)
            NPINANT= INT((-EANT-OMEGAMIN)/DOMEGA)+1
            NPINANT= MIN(MAX(1,NPINANT),NOMEGA-1)
            DWEIGHT= (-EANT-OMEGAMIN-(NPINANT-1)*DOMEGA)/DOMEGA

            SE_SPECTR(NPINANT,N1,K1,ISP)=SE_SPECTR(NPINANT,N1,K1,ISP)-&
                  CRES1*CONJG(CRES1)*P_ANT*(FSG_STORE(1)/W%WDES%NKPTS)*(1-DWEIGHT)
            SE_SPECTR(NPINANT+1,N1,K1,ISP)=SE_SPECTR(NPINANT+1,N1,K1,ISP)-&
                  CRES1*CONJG(CRES1)*P_ANT*(FSG_STORE(1)/W%WDES%NKPTS)*DWEIGHT

         ENDDO ! over K
      ENDDO ! over I
    ENDIF ! include head contribution

    CALLMPI( M_sum_d(W%WDES%COMM_INTER, SE_SPECTR, SIZE(SE_SPECTR)))

  ! 5) store self-energy for each state |n1, k1, ISP >
  ! the factor 1/2 originates from the Hartree terms being multiplied by 2 in the stored
  ! two electron four orbital matrix elements
  ! this means we obtain 4 times the 2nd order direct diagram G (V P_0 V) (instead of twice)
  ! and twice the 2nd order direct diagram G (I P_0 V) (instead of once)
    DO ISP=1,SEFROMBSE%ISPIN
       DO K1=1,SEFROMBSE%NKPTS
          DO N1=1,SEFROMBSE%NBANDS
             SEFROMBSE%REAL_PART(:,N1,K1,ISP) = REAL(SE_SPECTR(:,N1,K1,ISP),q)*(1.0_q/AHARTREE)
             ! the imaginary part should be zero, not quite obvious for beyond RPA
             SEFROMBSE%IMAG_PART(:,N1,K1,ISP) = 0
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(SE_SPECTR)
    IF (ALLOCATED(AHEAD)) DEALLOCATE(AHEAD)

 END SUBROUTINE SELFEN_FULL

 
!**************************************************************************
!
! Hilbert transformation of the imaginary part of the selfenergy
! specifically, the calling routine calculates the spectral function
! of the selfenergy
! the selfenergy is then obtained as
!
!  int_-infty^efermi d w' Sigma^s(w')/(w-w'+i delta)+
!  int_efermi^infty  d w' Sigma^s(w')/(w-w'-i delta)
!
! where Sigma^s(w') is the spectral function of the selfenergy 
! Sigma^s(w') = -1/pi Imag Sigma(w') sgn(w'-efermi) 
!
! Above the spectral function is defined such that a single delta peak
! yields a "norm" of int d w' delta(w') = 1 -> sum_i s(w_i) = 1 
! thus integrals are just replaced by sum_i
!
!**************************************************************************
 
 SUBROUTINE SE_HILBERT_TRANSF(SEFROMBSE, EFERMI, SHIFT, CELTOT_X, CELTOT_HARTREE_KINETIC)
  
    ! passed variables
    TYPE (selfenergy_from_bse) :: SEFROMBSE
           ! in: SEFROMBSE%REAL_PART  holds the spectral representation of the selfenergy
           ! out: SEFROMBSE%REAL_PART real part of selfenergy
           ! out: SEFROMBSE%IMAG_PART imaginary part of selfenergy
    REAL(q), INTENT(IN)    :: EFERMI                 ! Fermi energy
    REAL(q), INTENT(IN)    :: SHIFT                  ! complex shift for Lorenzian broadening
    COMPLEX(q), INTENT(IN) :: CELTOT_X(:,:,:)        ! mean field exchange contribtion
    COMPLEX(q), INTENT(IN) :: CELTOT_HARTREE_KINETIC(:,:,:) ! mean field kinetic and Hatree contribution
    ! local variables
    INTEGER :: I,J
    INTEGER :: K      ! momentum counter
    INTEGER :: N      ! band counter
    INTEGER :: ISP    ! spin variable
    REAL(q) :: RSIGN
    REAL(q) :: OMEGAMAX
    REAL(q), ALLOCATABLE :: OMEGA(:)                ! energies at which the self-energy is required
    COMPLEX(q) :: OMEGACPLX                         ! complex energy variable
    COMPLEX(q),ALLOCATABLE :: ZSE(:,:,:,:)          ! self-energy for occupied states 
    COMPLEX(q),ALLOCATABLE :: ZSE_SPECTR(:,:,:,:)   ! and its spectral representation

    !allocation
    ALLOCATE( OMEGA(SEFROMBSE%NOMEGA))
    ALLOCATE( ZSE(SEFROMBSE%NOMEGA, SEFROMBSE%NBANDS, SEFROMBSE%NKPTS, SEFROMBSE%ISPIN), &
             ZSE_SPECTR(SEFROMBSE%NOMEGA, SEFROMBSE%NBANDS, SEFROMBSE%NKPTS, SEFROMBSE%ISPIN))
    ! reconstruct the complex integrand from the separate Re and Im of the spectral representation

    ZSE=0._q 
    ZSE_SPECTR(:,:,:,:)=CMPLX(SEFROMBSE%REAL_PART(:,:,:,:),SEFROMBSE%IMAG_PART(:,:,:,:),q)

 ! frequency variable parameters
    OMEGA(:) = SEFROMBSE%OMEGA(:)
    OMEGAMAX = MAXVAL(OMEGA)

    DO ISP=1,SEFROMBSE%ISPIN
       DO K=1,SEFROMBSE%NKPTS
          DO N=1,SEFROMBSE%NBANDS
             DO J=1,SEFROMBSE%NOMEGA 
!!DIR$ IVDEP
!!OCL NOVREC
                DO I=1,SEFROMBSE%NOMEGA
                   IF(OMEGA(I)< EFERMI) THEN
                      RSIGN= 1.0_q
                   ELSE
                      RSIGN=-1.0_q
                   ENDIF
                   ZSE(J,N,K,ISP)=ZSE(J,N,K,ISP)+&
                        ZSE_SPECTR(I,N,K,ISP)/(-OMEGA(J)+CMPLX(OMEGA(I),RSIGN*SHIFT,q))
                ENDDO
             ENDDO

! incorrect version
!             DO J=1,SEFROMBSE%NOMEGA 
!                IF(OMEGA(J)< EFERMI) THEN
!                  RSIGN=-1.0_q
!                ELSE
!                  RSIGN=1.0_q
!                ENDIF
!                OMEGACPLX = CMPLX(OMEGA(J),RSIGN*SHIFT,q) 
!!DIR$ IVDEP
!!OCL NOVREC
!                DO I=1,SEFROMBSE%NOMEGA
!                   ZSE(J,N,K,ISP)=ZSE(J,N,K,ISP)+&
!                        ZSE_SPECTR(I,N,K,ISP)/(-OMEGACPLX+OMEGA(I))*0.5_q
!                ENDDO
!             ENDDO
          ENDDO
       ENDDO
    ENDDO

  ! copy back the self-energy for each state |w; n, k, ISP > into storage

    DO ISP=1,SEFROMBSE%ISPIN
       DO K=1,SEFROMBSE%NKPTS
          DO N=1,SEFROMBSE%NBANDS
             SEFROMBSE%REAL_PART(:,N,K,ISP) = REAL(ZSE(:,N,K,ISP),q)+CELTOT_X(N,K,ISP)+CELTOT_HARTREE_KINETIC(N,K,ISP)
             SEFROMBSE%IMAG_PART(:,N,K,ISP) = AIMAG(ZSE(:,N,K,ISP))
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(OMEGA, ZSE, ZSE_SPECTR)

 END SUBROUTINE SE_HILBERT_TRANSF


!**************************************************************************
! Find the quasi-particle energies at the intersection between the 
! real part of the selfenergy with the straight line y=w
!**************************************************************************

 SUBROUTINE QPE_FROMSIGMA(W, SEFROMBSE, EFERMI, CELTOT_X, CELTOT_HARTREE_KINETIC, NELM, IU0, IU6)
    IMPLICIT NONE
    ! passed variables
    TYPE (wavespin) :: W
    TYPE (selfenergy_from_bse) :: SEFROMBSE
    COMPLEX(q) :: CELTOT_X(:,:,:), CELTOT_HARTREE_KINETIC(:,:,:)
    REAL(q) :: EFERMI
    INTEGER :: NELM
    INTEGER :: IU0, IU6

    ! local variables
    REAL(q) :: EVAL        ! independent particle energy
    COMPLEX(q) :: QPE      ! quasi-particle energy
    REAL(q) :: Z           ! renormalisation factor
    REAL(q) :: QPE_LINEAR  ! QPE for linearized version
    REAL(q) :: SIGMA       ! selfenergy at independent particle energy
    REAL(q) :: Z_LINEAR    ! renormalisation factor for linearized version
    ! loop variables
    INTEGER :: ISP, K, N ! spin, k-point, band number

    DO ISP = 1, SEFROMBSE%ISPIN
       IF (W%WDES%ISPIN==2 .AND. IU6>=0 ) WRITE(IU6,'(/A,I1)') ' spin component ',ISP
       DO K = 1,SEFROMBSE%NKPTS
         IF (K > KPOINTS_ORIG%NKPTS) CYCLE
         IF (IU6>=0) WRITE(IU6,100) K, W%WDES%VKPT(:,K)

         DO N = 1,SEFROMBSE%NBANDS
            EVAL = SEFROMBSE%CELNEW(N,K,ISP)

            CALL BSE_LINEARISE_SE( SEFROMBSE%OMEGA, SEFROMBSE%REAL_PART(:, N, K, ISP), & 
                 EVAL, SIGMA, QPE_LINEAR, Z_LINEAR)
            CALL BSE_SEARCH_ROOT ( SEFROMBSE%OMEGA, SEFROMBSE%REAL_PART(:, N, K, ISP), SEFROMBSE%IMAG_PART(:, N, K, ISP),  &
                 EVAL, QPE, Z)

            SEFROMBSE%CELNEW(N,K,ISP) = REAL(QPE,q)
            IF ( IU6>=0 ) WRITE(IU6,200) &
                 N, EVAL, QPE_LINEAR, SIGMA, Z_LINEAR, & 
                 REAL(QPE,q), Z, W%FERTOT(N,K,ISP)*W%WDES%RSPIN, AIMAG(QPE)
         ENDDO
       ENDDO
    ENDDO
        
    IF (IU6 >= 0) WRITE(IU6, '(//)')

100 FORMAT(/' k-point ',I3,' :',3X,3F10.4/  &
   & " band No.  KS-energies   QP-e(linear)    sigma(KS)    Z         QP-e(zeros)     Z        occupation    Imag(E_QP)"/)
200 FORMAT((3X,I4,3X,9(F10.4,3X)))

 END SUBROUTINE QPE_FROMSIGMA


!**********************************************************************
!
! find the roots of the linearlized QP equation
!
!**********************************************************************

   SUBROUTINE BSE_LINEARISE_SE(OMEGA, SELFE, EVAL, SIGMA, QPE, Z)
     REAL(q), INTENT(IN) :: SELFE(:)  ! real part of selfenergy on grid
     REAL(q), INTENT(IN) :: OMEGA(:)  ! frequency grid
     REAL(q), INTENT(IN) :: EVAL      ! mean field eigenvalue
     REAL(q), INTENT(OUT) ::SIGMA     ! self-energy at mean field
     REAL(q), INTENT(OUT) :: QPE      ! QP energy
     REAL(q), INTENT(OUT) :: Z        ! renormalization factor
     ! local
     INTEGER :: NOMEGA, I, IEVAL, ID
     REAL(q) :: DOMEGA                ! distance to closest point
     REAL(q) :: SLOPE                 ! slope of selfenergy @ intersection 
     REAL(q) :: W                     ! frequency at which SELFE(w)=w
     
     NOMEGA=SIZE(OMEGA)
     DOMEGA=1E10
     ! search nearest frequency point on grid
     DO I = 1, NOMEGA
        IF(ABS(OMEGA(I)-EVAL) <= DOMEGA) THEN
           IEVAL =I
           DOMEGA=ABS(OMEGA(I)-EVAL)
        END IF
     ENDDO
     ! saveguard against under or overflow
     IEVAL=MAX(2,MIN(IEVAL,NOMEGA-1))
     ! determine slope (use two nearest points)
     ! seek a reasonable increment to get sensible Z
     DO ID=1,3
        SLOPE = (SELFE(IEVAL+ID)-SELFE(IEVAL-ID))/(OMEGA(IEVAL+ID)-OMEGA(IEVAL-ID))
        Z = 1.0_q/(1.0_q-SLOPE)
        IF (Z>0.4) CYCLE
     ENDDO
     ! estimate self-energy at EVAL from chosen points
     SIGMA = (SELFE(IEVAL-ID)*(OMEGA(IEVAL+ID)-EVAL)+SELFE(IEVAL+ID)*(EVAL-OMEGA(IEVAL-ID)))/(OMEGA(IEVAL+ID)-OMEGA(IEVAL-ID))
     ! now calculate linear estimate of intersection point
     QPE = EVAL + Z*(SIGMA-EVAL)
     
   END SUBROUTINE BSE_LINEARISE_SE

!**********************************************************************
!
! find the exact QP roots assuming picewise linearity in the self-energy
! i.e. determining the roots of  SELFE(w) -w = 0
! the routine first seeks the node closest to the
! mean field value EVAL and then seeks the root in the linearized
! function
! it is also possible to set EVAL first using BSE_LINEARISED_SE
! the best strategy would be to determine Z at each node and
! return the one with the largest Z factor
!
!**********************************************************************

   SUBROUTINE BSE_SEARCH_ROOT( OMEGA, SELFE, SELFE_IMAG, EVAL, QPE, Z)
     REAL(q), INTENT(IN) :: SELFE(:)  ! real part of selfenergy on grid
     REAL(q), INTENT(IN) :: SELFE_IMAG(:)  ! imaginary part of selfenergy on grid
     REAL(q), INTENT(IN) :: OMEGA(:)  ! frequency grid
     REAL(q), INTENT(IN) :: EVAL      ! mean field eigenvalue
     COMPLEX(q), INTENT(OUT) :: QPE   ! QP energy
     REAL(q), INTENT(OUT) :: Z        ! renormalization factor
     ! local<s
     INTEGER :: NOMEGA, I, IEVAL    
     REAL(q) :: DOMEGA                ! distance to closest point
     REAL(q) :: S                     ! store sign of SELFENERGY-OMEGA
     REAL(q) :: SLOPE                 ! slope of selfenergy @ intersection
     REAL(q) :: W                     ! frequency at which SELFE(w)=w
     
     NOMEGA=SIZE(OMEGA)
     ! search all nodes and select one closest to EVAL
     S=SIGN(1.0_q, SELFE(1)-OMEGA(1))
     IEVAL=1
     DOMEGA=1E10
     DO I = 2, NOMEGA
        IF (S*SIGN(1.0_q, SELFE(I)-OMEGA(I)) < 0) THEN
           S=SIGN(1.0_q, SELFE(I)-OMEGA(I)) ! store sign
           ! change of sign found a node in SELFE(w)-w
           IF (ABS(OMEGA(I)-EVAL)<DOMEGA) THEN
              DOMEGA=ABS(OMEGA(I)-EVAL)
              IEVAL=I
           ENDIF
        ENDIF
     ENDDO
     ! saveguard against under or overflow and underflow
     IEVAL=MAX(2,MIN(IEVAL,NOMEGA))
     ! determine slope of self-energy using neigh. points  
     SLOPE = (SELFE(IEVAL)-SELFE(IEVAL-1))/(OMEGA(IEVAL)-OMEGA(IEVAL-1))
     ! Z = 1/(w-self(w))'= 1/(1-self(w)')
     Z = 1.0_q/(1.0_q-SLOPE)
     ! estimate node of real part of SELFE(w)-w
     W   = OMEGA(IEVAL-1)+(SELFE(IEVAL-1)-OMEGA(IEVAL-1))*Z
     ! test whether correct W ! = interpolated selfenergy
     IF (ABS(W-(SELFE(IEVAL-1)+(SELFE(IEVAL)-SELFE(IEVAL-1)) & 
           *(W-OMEGA(IEVAL-1))/(OMEGA(IEVAL)-OMEGA(IEVAL-1))))>1E-10) THEN
        CALL vtutor%bug("internal error in BSE_SEARCH_ROOT: root not properly determined: W " // &
           str(SELFE(IEVAL-1)) // " " // str(SELFE(IEVAL)) // " " // str(OMEGA(IEVAL)) // " " // &
           str(OMEGA(IEVAL)), __FILE__, __LINE__)
     ENDIF
     
     QPE = CMPLX(W, (SELFE_IMAG(IEVAL-1)+(SELFE_IMAG(IEVAL)-SELFE_IMAG(IEVAL-1)) & 
          *(W-OMEGA(IEVAL-1))/(OMEGA(IEVAL)-OMEGA(IEVAL-1))),KIND=q)

   END SUBROUTINE BSE_SEARCH_ROOT

!**********************************************************************
!
! Write out BSE eigenfunctions and eigenvalues
! 			M.Bokdam Nov 2014
!**********************************************************************

  SUBROUTINE WRITE_BSE(WHF,NCV,BSE_INDEX,EIGVECT,R,NBSEEIG)

    USE constant
    USE full_kpoints
    USE kpoints_change
    USE pseudo
    IMPLICIT NONE

    ! passed variables
    TYPE (wavespin) WHF
    TYPE (latt) LATT_CUR
    TYPE(bse_matrix_index) :: BSE_INDEX
    INTEGER :: NCV
    INTEGER :: NBSEEIG
    COMPLEX(q) :: EIGVECT(:)
    REAL(q) :: R(:)
 ! local 
    INTEGER K1,K2,K3,K4,NPOS1,NPOS2,MINL 
    INTEGER I, J
    INTEGER     :: LAMBDA, ISP
    REAL(q) :: RIP(NCV),MINV 

      K1=1
      DO NPOS2=1,SIZE(R)
        K2=BSE_INDEX%NK(NPOS2)
        K3=BSE_INDEX%N1(NPOS2)
        K4=BSE_INDEX%N3(NPOS2)
        ISP=BSE_INDEX%ISP(NPOS2)
        RIP(NPOS2)=REAL(WHF%CELTOT(K4,K2,ISP)-WHF%CELTOT(K3,K2,ISP))
        K1=K1+1
      ENDDO
!MB very inefficient SORTING routine, but ok for now
      DO NPOS1=1,NBSEEIG
        MINV=MINVAL(RIP(NPOS1:NCV)) !locate
        MINL=INT(MINLOC(RIP(NPOS1:NCV),1))+NPOS1-1
        RIP(MINL)=RIP(NPOS1)        !swap
        RIP(NPOS1)=MINV
      ENDDO

      OPEN(unit=78, file='BSEFATBAND', status="UNKNOWN")
      WRITE(78,'(2I18)') SIZE(R),NBSEEIG
      K1=1
      DO NPOS1=1,NBSEEIG
        WRITE(78,'(I6,A14,F14.8,A20,F14.8)') NPOS1,"BSE eigenvalue:",R(NPOS1), "IP-eigenvalue:",RIP(NPOS1)
        DO NPOS2=1,SIZE(R)
          K2=BSE_INDEX%NK(NPOS2)
          K3=BSE_INDEX%N1(NPOS2)
          K4=BSE_INDEX%N3(NPOS2)
          ISP=BSE_INDEX%ISP(NPOS2)
          WRITE(78,'(3F9.5,3F14.7,2I6,F14.6,A5,F14.6)') WHF%WDES%VKPT(:,K2),REAL(WHF%CELTOT(K3,K2,ISP)), &
          & REAL(WHF%CELTOT(K4,K2,ISP)),ABS(EIGVECT(K1))/WHF%WDES%WTKPT(K2),K3,K4,REAL(EIGVECT(K1)), &
          & " +i* ",AIMAG(EIGVECT(K1))
          K1=K1+1
        ENDDO
      ENDDO
      CLOSE(78)
      
  END SUBROUTINE WRITE_BSE


!**********************************************************************
!
! Write out BSE eigenfunctions and eigenvalues
! .NOT. scaLAPACKaware version
! Based on CALCULATE_BSE_OSZI_STRENGHT
! 			M.Bokdam August 2018
!**********************************************************************

  SUBROUTINE WRITE_BSE_NOTSCALA(WHF, NQPOINT, KPOINT_BSE, ISP_LOW, ISP_HIGH,  &
       BD, BSE_INDEX, R, AMAT, NBSEEIG, IU6)

    USE constant
    USE full_kpoints
    USE kpoints_change
    USE pseudo
    USE mymath
    IMPLICIT NONE

    ! passed variables
    TYPE (wavespin) WHF
    INTEGER :: NQPOINT       ! q-point considered
    INTEGER :: KPOINT_BSE(3) ! G-point considered
!    TYPE (latt) LATT_CUR
    INTEGER ISP_LOW, ISP_HIGH
    TYPE(banddesc) :: BD(:)
    TYPE(bse_matrix_index) :: BSE_INDEX
    REAL(q) :: R(:)
    INTEGER :: NBSEEIG !number of BSE eigenvectors to be written
#ifdef single_prec_bse
    GDEFS :: AMAT(:,:)
#else
    GDEF  :: AMAT(:,:)
#endif
    INTEGER :: ANTIRES
    INTEGER :: NOMEGA_DIM
    REAL(q) :: OMEGAMAX
    INTEGER :: IU6
    ! local 
    INTEGER :: NCV13
    INTEGER K1, K3
    INTEGER N1, N3, ISP
    INTEGER NB1_INTO_TOT, NB3_INTO_TOT
    INTEGER I, J, NK1_IN_KPOINTS_FULL_ORIG, INFO
    INTEGER :: NCV,LAMBDA


      NCV=SIZE(AMAT,1)
      OPEN(unit=78, file='BSEFATBAND', status="UNKNOWN")
      WRITE(78,'(2I18)') SIZE(R),NBSEEIG

    lambda_bse: DO LAMBDA=1,NBSEEIG
        WRITE(78,'(I6,A16,F14.8)') LAMBDA," BSE eigenvalue:",R(LAMBDA)

        DO ISP=ISP_LOW,ISP_HIGH
          DO K1=1,WHF%WDES%NKPTS
            IF (WHF%WDES%WTKPT(K1)==0) CYCLE
            K3 =KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K1)+WHF%WDES%VKPT(:,NQPOINT),KPOINTS_FULL)
            ! determined the index of this k-point in the original full k-point grid
            NK1_IN_KPOINTS_FULL_ORIG=KPOINT_IN_FULL_GRID(WHF%WDES%VKPT(:,K1),KPOINTS_FULL_ORIG)

            DO N3=1,BD(ISP)%CBMAX-BD(ISP)%CBMIN+1
              NB3_INTO_TOT=N3+BD(ISP)%CBMIN-1
              IF (FILLED_XI_ORBITAL(WHF%FERTOT(NB3_INTO_TOT,K3,ISP))) CYCLE
              DO N1=1,BD(ISP)%VBMAX-BD(ISP)%VBMIN+1
                NB1_INTO_TOT=N1+BD(ISP)%VBMIN-1
                IF (EMPTY_XI_ORBITAL(WHF%FERTOT(NB1_INTO_TOT,K1,ISP))) CYCLE
                NCV13=BSE_INDEX%INDEX(NB1_INTO_TOT-BD(ISP)%VBMIN+1, NB3_INTO_TOT-BD(ISP)%CBMIN+1, K1, ISP -ISP_LOW+1)
                IF (NCV13==0) CYCLE ! cycle if this pair is not included   
#ifdef  gammareal             
                WRITE(78,'(3F9.5,3F14.7,2I6,F14.6)') WHF%WDES%VKPT(:,K1),REAL(WHF%CELTOT(NB1_INTO_TOT, K1, ISP)), &
             &  REAL(WHF%CELTOT(NB3_INTO_TOT,K3,ISP)),ABS(AMAT(NCV13,LAMBDA))/WHF%WDES%WTKPT(K1),NB1_INTO_TOT,NB3_INTO_TOT, &
             &  AMAT(NCV13,LAMBDA)
#else
                WRITE(78,'(3F9.5,3F14.7,2I6,F14.6,F14.6)') WHF%WDES%VKPT(:,K1),REAL(WHF%CELTOT(NB1_INTO_TOT, K1, ISP)), &
             &  REAL(WHF%CELTOT(NB3_INTO_TOT,K3,ISP)),ABS(AMAT(NCV13,LAMBDA))/WHF%WDES%WTKPT(K1),NB1_INTO_TOT,NB3_INTO_TOT, &
             &  REAL(AMAT(NCV13,LAMBDA)),AIMAG(AMAT(NCV13,LAMBDA))
#endif
              ENDDO  ! enddo for N3-loop
            ENDDO     ! enddo for N1-loop
          ENDDO  ! enddo K1 loop
        ENDDO  ! enddo ISP loop

      ENDDO lambda_bse
      CLOSE(78)

  END SUBROUTINE WRITE_BSE_NOTSCALA


#ifdef scaLAPACK
!**********************************************************************
!
! Write out BSE eigenfunctions and eigenvalues
! scaLAPACK aware version
! Based on CALCULATE_BSE_OSZI_SCALA
! 			M.Bokdam April  2016, Aug 2018
!**********************************************************************

  SUBROUTINE WRITE_BSE_SCALA(WHF, ISP_LOW, ISP_HIGH,  &
       BD, BSE_INDEX, R, AMAT, NCV, NBSEEIG, IU6, DESC)

    USE constant
    USE full_kpoints
    USE kpoints_change
    USE pseudo
    IMPLICIT NONE

    ! passed variables
    TYPE (wavespin) WHF
    INTEGER ISP_LOW, ISP_HIGH
    TYPE(banddesc) :: BD(:)
    TYPE(bse_matrix_index) :: BSE_INDEX
    REAL(q) :: R(:)
#ifdef single_prec_bse
    GDEFS :: AMAT(:)
#else
    GDEF  :: AMAT(:)
#endif
    INTEGER :: NCV
    INTEGER :: ANTIRES
    INTEGER :: IU6
 ! local 
    INTEGER :: NCV13
    INTEGER K1, ISP
    INTEGER N1, N3
    INTEGER K2,K3,K4,NPOS1,NPOS2 !MB
    INTEGER NB1_INTO_TOT, NB3_INTO_TOT
    INTEGER I, J
    INTEGER     :: LAMBDA
    INTEGER, OPTIONAL :: DESC(DLEN_)
!MB
    COMPLEX(q), ALLOCATABLE :: EIGVECT(:)
    INTEGER :: NBSEEIG

! BLACS variables
    INTEGER, EXTERNAL ::     NUMROC
    INTEGER MYROW, MYCOL, NPROW, NPCOL, NP,NQ
    INTEGER I1RES, J1RES, IROW, JCOL
    INTEGER I1, I2, J1, J2

    IF (NBSEEIG.GE.NCV) THEN
       IF (IU6>=0) THEN
          WRITE(*,*) "Warning: NBSEEIG cannot be larger than the number of BSE eigenstates."
          WRITE(*,*) "NBSEEIG is set to the total number of BSE eigenstates."
       ENDIF
       NBSEEIG=NCV
    ENDIF

!Large NBSEEIG value can overrun the memory of a local node.    
    ALLOCATE(EIGVECT(NBSEEIG*NCV)) 
    EIGVECT=0.0_q

    CALL BLACS_GRIDINFO(DESC(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
    NP = NUMROC(NCV,DESC(MB_),MYROW,0,NPROW)
    NQ = NUMROC(NCV,DESC(NB_),MYCOL,0,NPCOL)

!   loop over column indices (eigenvectors)
    DO LAMBDA=1, NBSEEIG
       ! local storage index in block (J1) and local index of block (J2)
       J1=   ((LAMBDA-DESC(NB_)*MYCOL-1)/(NPCOL*DESC(NB_)))*DESC(NB_)+1
       J2=MOD((LAMBDA-DESC(NB_)*MYCOL-1) ,NPCOL*DESC(NB_))+1

       JCOL=J1+J2-1
       IF (J2 <= 0 .OR. J2 > DESC(NB_)) THEN
          JCOL  =-1
       ELSE
          IF (LAMBDA /= DESC(NB_)*MYCOL+NPCOL*(J1-1)+J2) THEN
             CALL vtutor%bug("internal error in WRITE_BSE_SCALA: " // str(LAMBDA) // " " // &
                str(DESC(NB_)) // " " // str(NPCOL*(J1-1)) // " " // str(J2), __FILE__, __LINE__)
          ENDIF
       ENDIF

       IF (JCOL==-1) THEN
          ! no data on local node for eigenvector/eigenvalue pair LAMBDA_DONE
       ELSE
          IROW=0
          DO I1=1,NP,DESC(MB_)
             I1RES=MIN(DESC(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW=IROW+1
                
                NCV13=DESC(MB_)*MYROW+NPROW*(I1-1)+I2
                ISP= BSE_INDEX%ISP(NCV13)
                K1 = BSE_INDEX%NK(NCV13)

                IF (WHF%WDES%WTKPT(K1)==0) CYCLE
                
                NB1_INTO_TOT=BSE_INDEX%N1(NCV13)
                N1=NB1_INTO_TOT-BD(ISP)%VBMIN+1
                
                NB3_INTO_TOT=BSE_INDEX%N3(NCV13)
                N3=NB3_INTO_TOT-BD(ISP)%CBMIN+1
                      
                IF (NCV13/=BSE_INDEX%INDEX(NB1_INTO_TOT-BD(ISP)%VBMIN+1, NB3_INTO_TOT-BD(ISP)%CBMIN+1, K1, ISP -ISP_LOW+1)) THEN
                   CALL vtutor%bug("internal error in WRITE_BSE_SCALA: " // str(NCV13) // " " // &
                      str(NB1_INTO_TOT-BD(ISP)%VBMIN+1) // " " // str(NB3_INTO_TOT) // " " // str(K1) // " "&
                       // str(ISP) // " " // str(BSE_INDEX%INDEX(NB1_INTO_TOT-BD(ISP)%VBMIN+1,&
                      NB3_INTO_TOT-BD(ISP)%CBMIN+1,K1,ISP -ISP_LOW+1)), __FILE__, __LINE__)
                ENDIF

                IF (FILLED_XI_ORBITAL(WHF%FERTOT(NB3_INTO_TOT,K1,ISP))) CYCLE
                IF (EMPTY_XI_ORBITAL(WHF%FERTOT(NB1_INTO_TOT,K1,ISP)))  CYCLE
                
                EIGVECT((LAMBDA-1)*NCV+NCV13)=AMAT(IROW+(JCOL-1)*DESC(LLD_))
             ENDDO
          ENDDO
       ENDIF
    ENDDO

!Gather all parts of the eigenvectors stored on the different nodes
!    CALLMPI( M_sum_z(WHF%WDES%COMM, EIGVECT, NBSEEIG*NCV))
    CALLMPI( M_sum_z(WHF%WDES%COMM_INTER, EIGVECT, NBSEEIG*NCV))
    IF (IU6>=0) CALL WRITE_BSE(WHF,NCV,BSE_INDEX,EIGVECT,R,NBSEEIG)

    DEALLOCATE(EIGVECT)
 
  END SUBROUTINE WRITE_BSE_SCALA
#endif



      
END MODULE bse
