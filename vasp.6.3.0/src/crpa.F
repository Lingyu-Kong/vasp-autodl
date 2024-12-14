#include "symbol.inc"


!*********************************************************************
!
!  this module contains all necessary routines for the constrained RPA
!  mK
!
! TODO: implement LOCPROJ basis set
!*********************************************************************

MODULE crpa
   USE prec 
   USE mlwf
   IMPLICIT NONE
   !
   ! these INCAR flags handle CRPA flavors: 
   !
   ! this defines the target states in localized basis
   INTEGER, ALLOCATABLE, SAVE   :: NTARGET_STATES(:)
   !target CRPA bands, to be removed completely
   INTEGER, ALLOCATABLE, SAVE   :: NCRPA_BANDS(:)
   ! non-local or unscreened RPA 
   LOGICAL, SAVE                :: LNLRPA = .FALSE.
   ! projected CRPA method
   LOGICAL, SAVE                :: LPROJECTED = .FALSE.
   ! disentanglement of Miyake
   LOGICAL, SAVE                :: LDISENTANGLED = .FALSE.  
   ! weighted method of Friedrich and Shih
   LOGICAL, SAVE                :: LWEIGHTED= .FALSE.  
   !include Drude term in chi ?
   LOGICAL, SAVE                :: LDRUDE=.FALSE.
   !plot CRPA banstructure    
   LOGICAL, SAVE                :: LCRPAPLOT = .FALSE. 
   !regularize
   LOGICAL, PRIVATE, SAVE       :: LREGULARIZE = .TRUE.  
   !number of interations in Jacobi diagonalizer
   INTEGER, SAVE                :: MAXLIE = 150
   ! plot disentangled bands 
   LOGICAL, PRIVATE, SAVE       :: LPLOTDIS = .FALSE. 
   !drump FCI input
   LOGICAL, SAVE                :: LFCIINPUT = .FALSE. 
   !
   ! following flags are set internally set by VASP 
   !
   ! used basis sets, set by DMFT_BASIS tag in INCAR 
   !
   !use bloch basis  
   LOGICAL, SAVE                :: LUSEBLOCH = .FALSE.   
   !use MLWF basis (this is the default)
   LOGICAL, SAVE                :: LUSEMLWF = .FALSE.
   !use LOCPROJ basis 
   LOGICAL, SAVE                :: LUSELOCPROJ = .FALSE.
   !use orthogonal LOCPROJs
   LOGICAL, SAVE                :: LUSEORTH_LOCPROJS=.FALSE.
   !optionally the Wannier projections are read from UWAN
   LOGICAL, SAVE                :: LUFROMFILE=.FALSE.
   !correlated polariazbility needs to be computed?
   LOGICAL, SAVE                :: LCHIC=.FALSE. 
   !correlated polariazbility needs to be subtracted?
   LOGICAL, SAVE                :: LSUBTRACT_CHIC=.FALSE. 
   !correlated wavefunction is allocated?
   LOGICAL, SAVE                :: LALLOCATED_WCORR=.FALSE. 

   !relevant bands for correlated Green's function
   INTEGER, SAVE                :: NCRPAMIN = -1
   INTEGER, SAVE                :: NCRPAMAX = 1E8
 
   ! used to mark metallic states 
   LOGICAL,ALLOCATABLE :: LSELFCORRELATION_STATE(:,:,:)

   !Projection of wavefunction onto correlated sub subspace
   GDEF, ALLOCATABLE, SAVE      :: UCRPA(:,:,:,:)
   !Projection of derivative of wavefunction onto correlated sub subspace
   GDEF, ALLOCATABLE, SAVE      :: UCDER(:,:,:,:)
   !this maskes the target states in localized basis (like a bit mask)
   INTEGER, ALLOCATABLE, SAVE   :: NTARGET_MASK(:)
   !target specific Bloch states for disentangled CRPA
   LOGICAL, ALLOCATABLE, SAVE   :: LCRPA_STATE(:,:,:)
   !target whole Bloch bands for "hard" CRPA
   LOGICAL, ALLOCATABLE, SAVE   :: LCRPA_BAND(:)
   !weights for weighted CRPA
   REAL(q), ALLOCATABLE, SAVE   :: CRPA_WEIGHTS(:,:,:)
   !Wannier projections, used for projection and plotting 
   GDEF,ALLOCATABLE, SAVE       :: UWAN(:,:,:,:)
   !number of states
   INTEGER, SAVE                :: CRPA_NUM_WAN = -1 

   !some constants
   REAL(q), PRIVATE, PARAMETER  :: TINY=1E-11_q
   REAL(q), PRIVATE, PARAMETER  :: TINYLARGE=1E-6_q
 
CONTAINS

!*********************************************************************
!
! reader of important flags for CRPA 
!
!*********************************************************************
   SUBROUTINE CRPA_READER( INFO, IO )
      USE vaspxml
      USE full_kpoints
      USE mkpoints
      USE base
      USE mlwf, ONLY: LWANNIER90,LWANNIER90_RUN
      USE pead, ONLY: PEAD_REQUEST
      USE main_mpi, ONLY: NCSHMEM
      USE reader_tags
      USE tutor, ONLY: vtutor, isError, SHMEM_CRPA

      IMPLICIT NONE
      TYPE( in_struct ) IO
      TYPE( info_struct ) INFO
      ! local
      INTEGER IU5, IU6, IU0
      LOGICAL LCRPA, L2E4W, L2E4W_ALL
      INTEGER I, IERR, IDUM(1)
      INTEGER N, NDAT
      LOGICAL LOPEN
      CHARACTER (40)   INPLIN

      PROFILING_START('crpa_reader' )

      ! input output units
      IU0 = IO%IU0
      IU5 = IO%IU5
      IU6 = IO%IU6
      ! algorithm tags
      LCRPA = INFO%LCRPA
      L2E4W = INFO%L2E4W
      L2E4W_ALL = INFO%L2E4W_ALL

      !security break, SHMEM not supported for CRPA calculations
      IF ( (LCRPA .OR. L2E4W .OR. L2E4W_ALL) .AND. NCSHMEM > 1 )  THEN
         CALL vtutor%write(isError, SHMEM_CRPA)
      ENDIF 

      CALL OPEN_INCAR_IF_FOUND(IU5, LOPEN)

      !disentanglement approach is done
      LPROJECTED=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LPROJECTED', LPROJECTED, IERR, WRITEXMLINCAR)
   
      !disentanglement approach is done
      LDISENTANGLED=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LDISENTANGLED', LDISENTANGLED, IERR, WRITEXMLINCAR)
  
      !weighted approach is done
      LWEIGHTED=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LWEIGHTED', LWEIGHTED, IERR, WRITEXMLINCAR)
  
      ! non-local RPA is used
      LNLRPA=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LNLRPA', LNLRPA, IERR, WRITEXMLINCAR)

      !NCRPA_BANDS selects those Bloch bands, which are removed in the polarizability
      !first determine the number of arguments of NCRPA_BANDS
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NCRPA_BANDS', IDUM, 1, IERR, LWRITEXML=.FALSE., FOUNDNUMBER=NDAT )
      N=MAX(1,NDAT)
      ALLOCATE(NCRPA_BANDS(N))
      ! read arguments of NCRPA_BANDS
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NCRPA_BANDS', NCRPA_BANDS, N, IERR, WRITEXMLINCAR )
      ! in case tag not found, initialize to -1 
      IF ( NDAT == 0 ) THEN
         NCRPA_BANDS = -1
         ! set projector method as a default
         IF ( (.NOT. LDISENTANGLED ) .AND. (.NOT. LWEIGHTED) .AND. (.NOT. LNLRPA ) ) THEN
            LPROJECTED = .TRUE.
         ENDIF
      ENDIF

      ! plot bands used for correlated CHI in CRPA 
      LCRPAPLOT=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LCRPAPLOT', LCRPAPLOT, IERR, WRITEXMLINCAR)

      ! disentangled bands are plot 
      LPLOTDIS=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LPLOTDIS', LPLOTDIS, IERR, WRITEXMLINCAR)
      !also force downsampling of WAVECAR if CRPA bands are plotted 
      IF ( LCRPAPLOT .AND. .NOT. IO%LDOWNSAMPLE ) THEN
         IF ( IU0>=0 ) WRITE(IU0,*)'LDOWNSAMPLE forced by LCRPAPLOT '
         IO%LDOWNSAMPLE = LCRPAPLOT 
      ENDIF 

      ! disentangled bands are plot 
      LFCIINPUT=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LFCIINPUT', LFCIINPUT, IERR, WRITEXMLINCAR)
      IF ( LFCIINPUT .AND. .NOT. IO%LDOWNSAMPLE ) THEN
         IF ( IU0>=0 ) WRITE(IU0,*)'LDOWNSAMPLE forced by LFCIINPUT '
         IO%LDOWNSAMPLE = LFCIINPUT 
      ENDIF 

      !don't use orthogonal LOCPROJs by default
      LUSEORTH_LOCPROJS=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LUSEORTH_LOCPROJS', LUSEORTH_LOCPROJS, IERR, WRITEXMLINCAR)
   
      !NTARGET_STATES decides which Wannier states belong to target space in the constrained RPA
      ! first determine the number of arguments of NTARGET_STATES 
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NTARGET_STATES', IDUM, 1, IERR, LWRITEXML=.FALSE., FOUNDNUMBER=NDAT )
      N=MAX(1,NDAT)
      ALLOCATE(NTARGET_STATES(N))
      ! read arguments of NTARGET_STATES
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NTARGET_STATES', NTARGET_STATES, N, IERR, WRITEXMLINCAR )
      ! in case tag not found, initialize to -1 
      IF ( NDAT <= 0 ) THEN
         NTARGET_STATES = -1
      ENDIF

      !find out which basis is desired 
      INPLIN="--"
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LOCALIZED_BASIS', INPLIN, 40, IERR, WRITEXMLINCAR, FOUNDNUMBER=NDAT)
      !set correct Flags for basis choice
      CALL BASIS_FROM_INCAR(INPLIN, LCRPA, L2E4W, NDAT, IU0, IU6)
    
  
      ! regularization of correlated projectors
      LREGULARIZE=.TRUE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LREGULARIZE', LREGULARIZE, IERR, WRITEXMLINCAR)

      ! include Drude term (by default)
      LDRUDE=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LDRUDE', LDRUDE, IERR, WRITEXMLINCAR)

      IF ( LCRPA .OR. L2E4W .OR. L2E4W_ALL ) THEN
         CALL PEAD_REQUEST

         !force downsampling of WAVECAR 
         IF ( .NOT. IO%LDOWNSAMPLE ) THEN
            IF ( IU0>=0) WRITE(IU0,*)'Downsampling of WAVECAR forced by CRPA_READER'
            IF ( IU6>=0) WRITE(IU6,*)'Downsampling of WAVECAR forced by CRPA_READER'
            IO%LDOWNSAMPLE=.TRUE.
         ENDIF 
         IF( LUSEMLWF ) THEN
            !switch on wannier90 ( default)
            LWANNIER90=.TRUE.
            LWANNIER90_RUN=.TRUE.
#ifndef VASP2WANNIER90
            CALL vtutor%error("CRPA_READER: ERROR: VASP was compiled without wannier90 library, &
               &exiting now ...")
#endif

         ENDIF 
      ENDIF

      !maximum iterations in diagonalization routines
      MAXLIE = 150
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MAXLIE', MAXLIE, IERR, WRITEXMLINCAR)
  
      CALL CLOSE_INCAR_IF_FOUND(IU5)

      PROFILING_STOP('crpa_reader' )
   END SUBROUTINE CRPA_READER

!**********************************************************************
! initializes NTARGET_MASK, after NBANDS is known 
! also sets CRPA flavor related flags 
!
!**********************************************************************
   SUBROUTINE INIT_CRPA_STATES( NBANDS, IU0, IU6 ) 
      USE tutor, ONLY: vtutor, isError, isAlert, TargetNotSet, TargetOut, NCRPA_BANDSout
      INTEGER NBANDS 
      INTEGER IU0, IU6
      INTEGER I     
      INTEGER, ALLOCATABLE :: TMP(:)
      INTEGER              :: NWAN=0
      PROFILING_START( 'init_crpa_states' )

      ! mask NTARGET_MASK based on selection of NTARGET_STATES 
      IF ( ALLOCATED( NTARGET_MASK ) ) DEALLOCATE( NTARGET_MASK )
      ALLOCATE( NTARGET_MASK( NBANDS ) ) 
      NTARGET_MASK = 0 
      
      ! in case NTARGET_STATES was not set by the user 
      IF ( SIZE( NTARGET_STATES ) == 1 .AND. NTARGET_STATES(1) == -1  ) THEN
         DEALLOCATE( NTARGET_STATES ) 
         ALLOCATE( NTARGET_STATES( NBANDS ) )
         CALL vtutor%write(isAlert, TargetNotSet)
         DO I = 1, NBANDS
            NTARGET_STATES(I) = I 
         ENDDO
      ! in case it was set, 
      ELSE IF ( SIZE( NTARGET_STATES ) > NBANDS ) THEN

         CALL vtutor%error("ERROR: More NTARGET_STATES given than NBANDS available")
      ELSE
         ALLOCATE(TMP( SIZE(NTARGET_STATES) ) )
         TMP = NTARGET_STATES 
         DEALLOCATE( NTARGET_STATES ) 
         ALLOCATE( NTARGET_STATES( NBANDS ) ) 
         NTARGET_STATES = -1 
         NTARGET_STATES( 1:SIZE( TMP )  ) = TMP( 1: SIZE(TMP) ) 
         DEALLOCATE( TMP ) 
      ENDIF  

      DO I = 1, NBANDS
         IF( NTARGET_STATES( I ) <= NBANDS)  THEN
            IF ( NTARGET_STATES( I ) > 0 ) NTARGET_MASK( NTARGET_STATES(I) ) = 1   
         ELSE
            CALL vtutor%write(isError, TargetOut)
         ENDIF
      ENDDO

      ! 
      ! set LCRPA_BAND based on NCRPA_BANDS
      ! 
      IF ( ALLOCATED( LCRPA_BAND ) ) DEALLOCATE( LCRPA_BAND )
      ALLOCATE( LCRPA_BAND( NBANDS ) ) 
      LCRPA_BAND = .FALSE. 
      
      ! in case NCRPA_BANDS was not set by the user
      IF ( SIZE( NCRPA_BANDS ) == 1 .AND. NCRPA_BANDS(1) == -1  ) THEN
         ! in case bloch basis is used, identify LCRPA_BAND with NTARGET_STATES
         IF( LUSEBLOCH .AND. .NOT. LNLRPA ) THEN
            DO I = 1 , NBANDS 
               IF( NTARGET_STATES( I ) <= NBANDS)  THEN
                  IF ( NTARGET_STATES( I ) > 0 ) LCRPA_BAND( NTARGET_STATES(I) ) = .TRUE.
               ELSE
                  CALL vtutor%write(isError, TargetOut)
               ENDIF
            ENDDO
         ENDIF
      ! in case it was not set properly
      ELSE IF ( SIZE( NCRPA_BANDS ) > NBANDS ) THEN
         CALL vtutor%error("ERROR: More NCRPA_BANDS given than NBANDS available")
      ! in case it was set proberly 
      ELSE
         DO I = 1, SIZE( NCRPA_BANDS )
            IF( NCRPA_BANDS( I ) <= NBANDS)  THEN
               ! mark the band as CRPA band
               IF ( NCRPA_BANDS( I ) > 0 ) LCRPA_BAND( NCRPA_BANDS(I) ) = .TRUE.
            ELSE
               CALL vtutor%write(isError, NCRPA_BANDSout)
            ENDIF
         ENDDO
      ENDIF  
      PROFILING_STOP( 'init_crpa_states' )
   END SUBROUTINE INIT_CRPA_STATES

!**********************************************************************
!
! small helper routine, set the basis for Coulomb integrals 
! and DMFT interaction matrix elements   
! 
!**********************************************************************

   SUBROUTINE BASIS_FROM_INCAR( TEXT, LCRPA, L2E4W, NDAT, IU0, IU6 )
      USE fock
      USE locproj, ONLY: LPRJ_GET_NUM_WAN, WANPROJFILE
      USE tutor, ONLY: vtutor, isError, isAlert, LOCPROJnotSet, BlochBasisUsed
      CHARACTER (LEN=*) :: TEXT
      LOGICAL           :: LCRPA
      LOGICAL           :: L2E4W
      INTEGER           :: NDAT, IU0, IU6
      INTEGER, EXTERNAL :: LENGTH
      LOGICAL           :: LWANPROJ_FOUND = .FALSE.
      INTEGER           :: NUM_WAN
      INTEGER           :: N
 
      CALL STRIP(TEXT,N,'L')
      CALL LOWER(TEXT)
      N=LENGTH(TEXT)
      LUSELOCPROJ = .FALSE.
      LUSEBLOCH = .FALSE.
      LUSEMLWF = .FALSE.
 
      IF ( TEXT(1:N)=='mlwf' ) THEN
         LUSEMLWF = .TRUE.
#ifndef VASP2WANNIER90
         ! pull breaks here, if VASP was not compiled with wannier90 support 
         CALL vtutor%error("VASP was not compiled with VASP2WANNIER90!")
#endif
         IF ( LCRPA .AND. IU0 >= 0 ) &
            WRITE(IU0,*)'Coulomb integrals are evaluated in MLWF basis'
      ELSEIF( TEXT(1:N)=='lprj' )THEN
         LUSELOCPROJ = .TRUE.
         IF ( LCRPA .AND. IU0 >= 0 ) &
            WRITE(IU0,*)'Coulomb integrals are evaluated in basis defined by LOCPROJ'
         IF( LUSEORTH_LOCPROJS ) THEN
           
            IF ( LCRPA .AND. IU0 >= 0 ) &
               WRITE(IU0,*)'orthogonalized LOCPROJ functions used'
         ENDIF 
         NUM_WAN = LPRJ_GET_NUM_WAN()
         ! print error message if locproj is chosen but not set in INCAR
         IF ( LCRPA .AND. NUM_WAN == 0 ) THEN
            CALL vtutor%write(isError, LOCPROJnotSet)
         ENDIF
      ! bloch basis is chosen 
      ELSEIF( TEXT(1:N)=='bloch' )THEN
         IF ( LCRPA .AND. IU0 >= 0 ) &
            WRITE(IU0,*)'Coulomb integrals are evaluated in Bloch basis'
         LUSEBLOCH = .TRUE.

      ! otherwise, check if wannier90 or LOCPROJ is set
      ELSE
         IF ( LCRPA .OR. L2E4W .OR. (NDAT == 0 .AND. SUM( NTARGET_STATES(:)) >0) ) THEN
            ! check out if wannier90.win or WANPROJ file is present 
            INQUIRE(FILE=WANPROJFILE, EXIST = LUFROMFILE )
            ! locproj settings in INCAR are favored
            ! if at least one entry is found, select LOCPROJ basis
            NUM_WAN = LPRJ_GET_NUM_WAN()
            IF ( NUM_WAN > 0 ) THEN
               LUSELOCPROJ = .TRUE. 
               IF ( IU0 >= 0 ) &
                  WRITE(IU0,*)'Coulomb integrals are evaluated in basis defined by LOCPROJ'
            ! if no LOCPROJ entry is found, look for a WANPROJ file
            ELSE IF ( LUFROMFILE ) THEN
!! TODO: USE file written by LOCPROJ and get rid of WANPROJFILE 
               IF ( IU0 >= 0 ) &
                  WRITE(IU0,*)'found '//WANPROJFILE//' file, will use it for CRPA'
               
            ELSE IF ( LEN_TRIM(WANNIER90_WIN) > 0 ) THEN
#ifndef VASP2WANNIER90
               ! pull breaks here, if VASP was not compiled with wannier90 support 
               CALL vtutor%error("VASP was not compiled with VASP2WANNIER90!")
#endif
               LUSEMLWF = .TRUE. 
               IF ( IU0 >= 0 ) &
                  WRITE(IU0,*)'found wannier90 input file, will use it for CRPA'
            ELSE
               LUSEBLOCH = .TRUE. 
               CALL vtutor%write(isAlert, BlochBasisUsed)
            ENDIF
         ENDIF
      ENDIF 

   END SUBROUTINE BASIS_FROM_INCAR 

!**********************************************************************
!
! set flags for calculation of effective interaction parameter U,U',J
!
!**********************************************************************

   SUBROUTINE CRPA_SETUP( WDES, W, W_CORR, LMDIM, T_INFO, P, CQIJ,&
      LATT_CUR, GRID, INFO, KPOINTS, NQ_IBZ, IO )  
      USE base
      USE pseudo
      USE wave_high
      USE mkpoints
      USE lattice
      USE ini
      USE constant
      USE poscar
      USE mlwf, ONLY: READ_WANPROJ_FILE
      USE locproj
      USE full_kpoints
      USE fileio
      !needed for interpolation
      USE hamil_struct_def
      USE mgrid
      USE nonl_high
      USE vaspxml
      TYPE (wavedes)     WDES
      TYPE (wavespin)    W
      TYPE (wavespin)    W_CORR
      INTEGER            LMDIM
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P(T_INFO%NTYP)
      OVERLAP  CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      TYPE (latt)        LATT_CUR
      TYPE (grid_3d)     GRID
      TYPE (info_struct) INFO
      TYPE (kpoints_struct) KPOINTS
      INTEGER            NQ_IBZ
      TYPE (in_struct)   IO
      !local
      INTEGER, ALLOCATABLE :: WTARGET(:)
      INTEGER              :: I, COUNTER, NK,ISP
      PROFILING_START ('crpa_setup')

      IF ( .NOT. INFO%LCRPA ) THEN
         PROFILING_STOP( 'crpa_setup' )
         RETURN 
      ENDIF
    
      CALL INIT_CRPA_STATES( WDES%NB_TOT, IO%IU0, IO%IU6 )

      ! only in the following cases correlated polarizability has to be calculated 
      !
      ! For INFO%LCRPA = .TRUE. the following table shows a summary whether LCHIC
      ! is calculated or not
      !
      ! LRPOJECTED LDISENTANGLED LWEIGHTED LNLRPA | LCHIC
      !     F            F           F       F        F
      !     T            F           F       F        T
      !     F            T           F       F        F
      !     F            F           T       F        F
      !     F            F           F       T        F
      !     F            T           F       T        F
      !     T            F           F       T        F
      !set default CRPA flavor
      LCHIC = LPROJECTED 
      ! correlated polarizability is always needed in case spacetime GW is executed
      IF ( INFO%ICHIREAL > 1 ) LCHIC = .TRUE. 
      ! except for NLRPA computaitons 
      IF ( LNLRPA ) LCHIC = .FALSE.
      !---------------------------------------------------------------------------------------
      IF ( LUSEBLOCH ) THEN  !use Bloch basis 
      !---------------------------------------------------------------------------------------
         IF ( INFO%ICHIREAL == 0  ) LCHIC = .FALSE.
         ! if the interaction parameter are computed in the Bloch basis 
         ! we don't need Wannier functions
         ! number of basis elements is then the total # of bands 
         CRPA_NUM_WAN = WDES%NB_TOT
        
         !set targetbands from NTARGET_STATES if not full RPA is done
         IF ( .NOT. LNLRPA ) THEN
            CALL PRINT_INFO_ABOUT_BANDS
         ELSE
            IF ( IO%IU0>=0 ) WRITE(IO%IU0,'(" Target bands: ")',ADVANCE="NO")
            IF ( IO%IU6>=0 ) WRITE(IO%IU6,'(" Target bands: ")',ADVANCE="NO")

            CALL SET_RELEVANT_BANDS_FOR_G( .TRUE. )
         ENDIF
         IF ( IO%IU0>=0 ) WRITE(IO%IU0,*)
         IF ( IO%IU6>=0 ) WRITE(IO%IU6,*)
         
      !---------------------------------------------------------------------------------------
      ELSE      ! use Wannier or LOCPROJ basis 
      !---------------------------------------------------------------------------------------
         !otherwise we need Wannier functions 
         IF( LUSEMLWF ) THEN
             !first, look for an existing WANPROJ file
             CALL READ_WANPROJ_FILE(WANPROJFILE,WDES,IO,LUFROMFILE,CRPA_NUM_WAN,UWAN)
             ! project only if file not present
             IF ( .NOT. LUFROMFILE ) THEN
                CALL MLWF_MAIN( WDES, W, P, KPOINTS, CQIJ, T_INFO, LATT_CUR, INFO, IO, MLWF_GLOBAL )
                !set CRPA_NUM_WAN to CRPA_NUM_WAN
                CRPA_NUM_WAN = MLWF_GLOBAL%NUM_WANN
             ENDIF 
         ELSEIF( LUSELOCPROJ) THEN
             !first, look for an existing WANPROJ file
             CALL READ_WANPROJ_FILE(WANPROJFILE,WDES,IO,LUFROMFILE,CRPA_NUM_WAN,UWAN)
             ! project only if file not present
             IF ( .NOT. LUFROMFILE ) THEN
                CALL LPRJ_PROALL(W,WDES,P,CQIJ,LATT_CUR,LPRJ_functions,LPRJ_COVL,T_INFO,INFO,IO)
                CRPA_NUM_WAN = LPRJ_GET_NUM_WAN()
             ENDIF
         ELSE
             ! read wanproj file and set CRPA_NUM_WAN
             CALL READ_WANPROJ_FILE(WANPROJFILE,WDES,IO,LUFROMFILE,CRPA_NUM_WAN,UWAN)
         ENDIF 
       
         IF ( IO%IU0>=0 .AND. LUFROMFILE ) WRITE(IO%IU0,*)'using projections found in '//WANPROJFILE

         !maybe some bands are removed completely (highest priority)
         CALL PRINT_INFO_ABOUT_BANDS
         CALL SET_RELEVANT_BANDS_FOR_G( .FALSE. )

         IF ( LWEIGHTED ) THEN
            !this computes Friedrich and Shih's CRPA weights
            CALL DETERMINE_CRPA_WEIGHTS( WDES, W, NQ_IBZ, IO ) 
            IF ( IO%IU0>=0 ) WRITE(IO%IU0,'(A)')' polarizability determined '//&
               'using CRPA weights'
            IF ( IO%IU6>=0 ) WRITE(IO%IU6,'(A)')' polarizability determined '//&
               'using CRPA weights' 
            LCHIC = .FALSE.
         ELSE IF ( LDISENTANGLED .OR. LPROJECTED ) THEN
            !otherwise disentanglement or projector method is done.
            !additional wavefunction is necessary in these cases
            CALL ALLOCW( W%WDES, W_CORR )
            W_CORR%WDES  =W%WDES
            ! copy wavefunction to W_LOC, W_LOC is going to be the projected wave function
            W_CORR%CPTWFP=W%CPTWFP
            W_CORR%CPROJ =W%CPROJ 
            W_CORR%FERWE =W%FERWE
            W_CORR%AUX   =W%AUX
            W_CORR%CELEN =W%CELEN
            W_CORR%FERTOT=W%FERTOT
            W_CORR%AUXTOT=W%AUXTOT
            W_CORR%CELTOT=W%CELTOT
            W_CORR%OVER_BAND=W%OVER_BAND
            ! mark that allocated wavefunction is allocated 
            LALLOCATED_WCORR = .TRUE. 

            !now decide whether to disentangle the correlated wavefunction or
            !to use the projector approach to subdivide the hilbert space
            IF ( LDISENTANGLED ) THEN 
               !disentangle the wavefunction using Miyake's disentanglement approach
               CALL DISENTANGLE( WDES, W_CORR, LMDIM, CQIJ, LATT_CUR, & 
                  INFO, KPOINTS, NQ_IBZ, IO )  
               IF ( IO%IU0>=0 ) WRITE(IO%IU0,'(A)')' polarizability constrained '//&
                  'using disentangled wavefunctions'
               IF ( IO%IU6>=0 ) WRITE(IO%IU6,'(A)')' polarizability constrained '//&
                  'using disentangled wavefunctions'
            ELSEIF ( LPROJECTED ) THEN
               !remove correlated polarizability from full polarizability
#ifdef old_version
               !replace Bloch states by projected Wannier states
               CALL SEPARATE_STATES( WDES, W_ALL , W_CORR, LMDIM, CQIJ, LATT_CUR, & 
                  INFO, KPOINTS, NQ_IBZ, IO )  
#else
               !use numerically more stable version. 
               CALL SEPARATE_CORRELATED_STATES( WDES, W , W_CORR, LMDIM, CQIJ, LATT_CUR, &
                  INFO, KPOINTS, NQ_IBZ, IO )  
#endif

               IF ( IO%IU0>=0 ) WRITE(IO%IU0,'(A)')' polarizability constrained '//&
                  'using projected wavefunctions'
               IF ( IO%IU6>=0 ) WRITE(IO%IU6,'(A)')' polarizability constrained '//&
                  'using projected wavefunctions'
            ENDIF 
         ENDIF
         !print information about correlated subspace
         CALL PRINT_CSPACE_INFO(IO%IU0,IO%IU6)
      !---------------------------------------------------------------------------------------
      ENDIF   ! basis choice 
      !---------------------------------------------------------------------------------------
 
      !maybe RPA is done
      IF( LNLRPA ) THEN 
         !maybe CRPA_WEIGHTS are allocated -> get rid of it 
         IF ( ALLOCATED(CRPA_WEIGHTS) ) DEALLOCATE(CRPA_WEIGHTS)
         IF ( ALLOCATED(LCRPA_BAND) ) DEALLOCATE(LCRPA_BAND)
         !in case of disentangled RPA we need this
         IF ( ALLOCATED(LCRPA_STATE) ) LCRPA_STATE=.FALSE.
         IF ( ALLOCATED(UCRPA) .AND. .NOT. (LDISENTANGLED .OR. LNLRPA) ) DEALLOCATE(UCRPA) 
         IF ( IO%IU0>=0 ) WRITE(IO%IU0,'(A)')' RPA polarizability is used '
         IF ( IO%IU6>=0 ) WRITE(IO%IU6,'(A)')' RPA polarizability is used '
      ENDIF 
      
      !stop here if Wannier projection matrices are only written
      IF ( LCRPAPLOT ) THEN
         CALL STOP_XML
         PROFILING_STOP ('crpa_setup')
         IF ( IO%IU0>=0 ) WRITE(IO%IU0,'(A)')'U.WAN written, proceed with plotting bands'
         IF ( IO%IU6>=0 ) WRITE(IO%IU6,'(A)')'U.WAN written, proceed with plotting bands'
         CALL vtutor%stopCode()
      ENDIF 

      PROFILING_STOP ('crpa_setup')
      CONTAINS

      SUBROUTINE PRINT_INFO_ABOUT_BANDS
         INTEGER :: I, J , K
         ! count them
         K=0
         DO I = 1, WDES%NB_TOT
            IF( LCRPA_BAND( I ) ) K=K+1
         ENDDO

         DO I = 1, WDES%NB_TOT
            IF( LCRPA_BAND( I ) ) EXIT
         ENDDO
          
         ! this signals that there is at least one CRPA band
         IF ( I <= WDES%NB_TOT ) THEN
            IF ( IO%IU0>=0 ) THEN
               WRITE(IO%IU0,'(A)')' Bands neglected within CRPA:'
               J = 0 
               DO I = 1, WDES%NB_TOT
                  IF( LCRPA_BAND(I) ) THEN
                     J = J + 1 
                     IF ( J == K ) THEN
                        WRITE(IO%IU0,'(I5)') I
                        J=0                      
                     ELSE
                        WRITE(IO%IU0,'(I5)',ADVANCE="No") I
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
            IF ( IO%IU6>=0 ) THEN
               WRITE(IO%IU6,'(A)')' Bands neglected within CRPA:'
               J = 0 
               DO I = 1, WDES%NB_TOT
                  IF( LCRPA_BAND(I) ) THEN
                     J = J + 1 
                     IF ( J == K ) THEN
                        WRITE(IO%IU6,'(I5)') I
                        J=0                      
                     ELSE
                        WRITE(IO%IU6,'(I5)',ADVANCE="No") I
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      END SUBROUTINE PRINT_INFO_ABOUT_BANDS

      SUBROUTINE SET_RELEVANT_BANDS_FOR_G( LBLOCHBASIS )
         LOGICAL :: LBLOCHBASIS
         !relevant bands for correlated Green's function
         NCRPAMIN=WDES%NB_TOT
         NCRPAMAX= 1 
         IF ( LBLOCHBASIS ) THEN
            DO I = 1, CRPA_NUM_WAN
               IF ( NTARGET_MASK (I) == 1 ) THEN
                  IF ( IO%IU0>=0 ) WRITE(IO%IU0,'(I4)',ADVANCE="NO")I
                  IF ( IO%IU6>=0 ) WRITE(IO%IU6,'(I4)',ADVANCE="NO")I
                  NCRPAMIN=MIN(NCRPAMIN,I)
                  NCRPAMAX=MAX(NCRPAMAX,I)
               ENDIF
            ENDDO
            IF ( IO%IU0>=0 ) WRITE(IO%IU0,*)
            IF ( IO%IU6>=0 ) WRITE(IO%IU6,*)
         ELSE
            DO I = 1, WDES%NB_TOT
               IF ( LCRPA_BAND(I) ) THEN
                  NCRPAMIN=MIN(NCRPAMIN,I)
                  NCRPAMAX=MAX(NCRPAMAX,I)
               ENDIF
            ENDDO
         ENDIF
      END SUBROUTINE SET_RELEVANT_BANDS_FOR_G

   END SUBROUTINE CRPA_SETUP

!**********************************************************************
!
! disentanglement of W into correlated and non-correlated part 
! using Miyake's approch from PRB. 80, 155134 (2009) 
! 
!**********************************************************************

   SUBROUTINE DISENTANGLE( WDES, W , LMDIM, CQIJ, LATT_CUR, & 
      INFO, KPOINTS, NQ_IBZ, IO )  
      USE base
      USE wave_high
      USE full_kpoints
      USE lattice
      USE constant
      USE mlwf!, ONLY : U_matrix,U_matrix_opt,lexclude_band,lwindow,kpt_latt,num_kpts
      USE choleski
      USE m_mrgrnk 
      USE mkpoints 
      USE ini
      USE tutor, ONLY: vtutor, isError, argument, CRPADisentangle
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      INTEGER LMDIM
      OVERLAP  CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      TYPE (latt) LATT_CUR
      TYPE (info_struct) INFO
      TYPE (kpoints_struct) KPOINTS
      INTEGER NQ_IBZ
      TYPE(in_struct) IO
      ! local variables
      TYPE(wavedes1) WDES1
      TYPE(wavefuna) WA
      INTEGER I,J,K,IP
      INTEGER ISP,NK
      COMPLEX(q), ALLOCATABLE :: NEWEIG(:,:,:)
      COMPLEX(q), ALLOCATABLE :: NEWEIGC(:,:,:)
      REAL(q), ALLOCATABLE ::  COREIG(:)  !correlated eigenvalues
      REAL(q), ALLOCATABLE :: NCOREIG(:)  !noncorrelated
      REAL(q), ALLOCATABLE :: ALLEIG(:)  !disentangled eigenvalues 
      GDEF, ALLOCATABLE :: U(:,:)
      GDEF, ALLOCATABLE :: UNCORR(:,:)
      GDEF, ALLOCATABLE :: UPLOT(:,:)
      INTEGER,ALLOCATABLE:: NSORT(:)
      INTEGER :: NROWMIN1, NROWMAX1
      INTEGER :: NCOLMIN1, NCOLMAX1
      INTEGER :: NMIN, NMAX
      TYPE(argument) arg
     
      CALL CHECK_FULL_KPOINTS
  
      !consistency check: Disentanglement needs a minimial basis 
      !that is: the Wannier basis should not be larger than the correlated sub space
      J = 0
      DO I = 1, CRPA_NUM_WAN
         IF ( NTARGET_MASK(I) == 1 ) J = J+1
      ENDDO 
      IF ( J /= CRPA_NUM_WAN ) THEN
         ALLOCATE(arg%ival(2))
         arg%ival(1)=CRPA_NUM_WAN
         arg%ival(2)=J
         CALL vtutor%write(isError, CRPADisentangle, arg)
      ENDIF 

      !allocate wannier projection matrices
      ALLOCATE(U(WDES%NB_TOT,WDES%NB_TOT))
      ALLOCATE(UNCORR(WDES%NB_TOT,WDES%NB_TOT))
      ALLOCATE(UPLOT(WDES%NB_TOT,WDES%NB_TOT))
  
      !nullify CRPA projector and reallocate it
      IF ( ALLOCATED ( UCRPA ) ) DEALLOCATE(UCRPA) 
      ALLOCATE(UCRPA(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS, WDES%ISPIN) )
      UCRPA = 0
  
      !new eigenvalues
      ALLOCATE(NEWEIG(WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN))
      NEWEIG=zero
      ALLOCATE(NEWEIGC(WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN))
      NEWEIGC=zero
  
      !entangled and disentangled eigenvalues 
      ALLOCATE(ALLEIG(WDES%NB_TOT))    
      ALLOCATE(COREIG(WDES%NB_TOT))    
      ALLOCATE(NCOREIG(WDES%NB_TOT))    
  
      !sorting array 
      ALLOCATE( NSORT( WDES%NB_TOT ) ) 
      IF (ALLOCATED(LCRPA_STATE)) DEALLOCATE(LCRPA_STATE)
      ALLOCATE( LCRPA_STATE( WDES%ISPIN, WDES%NKPTS, WDES%NB_TOT  ) )
      LCRPA_STATE=.FALSE.
  
      !find largest set of bloch bands considered 
      NROWMIN1 = WDES%NB_TOT
      NROWMAX1 = 1 
      NCOLMIN1 = WDES%NB_TOT
      NCOLMAX1 = 1 
      NMIN=WDES%NB_TOT
      NMAX=1
      ! find sub block of Wannier functions
      DO ISP=1,WDES%ISPIN
         ksearch: DO NK=1,WDES%NKPTS
            !determine local projector U
            CALL DETERMINE_WANNIER_PROJECTOR(WDES,NK,ISP,UCRPA(:,:,NK,ISP),IO%IU0)
            !----------------------------------------------------------------------  
            !construct overlap
            U=zero
            DO I=1,WDES%NB_TOT
               DO J=1,WDES%NB_TOT
                  DO K = 1, WDES%NB_TOT
                     U(I,J) = U(I,J) + UCRPA(I,K,NK,ISP)*GCONJG(UCRPA(J,K,NK,ISP))
                  ENDDO 
               ENDDO
            ENDDO
            !extract nonzero CRPA_NUM_WAN x CRPA_NUM_WAN block
            !find largest set of bands considered 
            DO I = 1, WDES%NB_TOT
               DO  J = 1, WDES%NB_TOT
                  IF ( ABS( U(I,J)) > TINY ) THEN
                      NROWMIN1 = MIN( I, NROWMIN1 )
                      NROWMAX1 = MAX( I, NROWMAX1 )
                  ENDIF
               ENDDO
            ENDDO 
            NMIN=MIN(NMIN,NROWMIN1)
            NMAX=MAX(NMAX,NROWMAX1)
         ENDDO ksearch
      ENDDO
      !block not found if NROWMAX1 - NROWMIN1 + 1 /= CRPA_NUM_WAN
      IF (  NROWMAX1-NROWMIN1+1< CRPA_NUM_WAN ) THEN
         CALL vtutor%bug("internal error in VASP from DISENTANGLE: sub block not found", __FILE__, __LINE__)
      ENDIF 

      !relevant bands for correlated Green's function
      NCRPAMIN=NROWMIN1 
      NCRPAMAX=NROWMAX1 

      !----------------------------------------------------------------------------------------
      ! loop over spin and all k-points
      spin: DO ISP=1,WDES%ISPIN
      !----------------------------------------------------------------------------------------
         kpts: DO NK=1,WDES%NKPTS
         !-------------------------------------------------------------------------------------
            !determine local projector U
            CALL DETERMINE_WANNIER_PROJECTOR(WDES,NK,ISP,UCRPA(:,:,NK,ISP),IO%IU0)
            !to plot correlated states 
            IF ( LCRPAPLOT ) UPLOT(:,:) = UCRPA(:,:,NK,ISP)
         
            !orthogonalize correlated sub space
            CALL DIAGONALIZE_CORR_SYSTEM(W,UCRPA(:,:,NK,ISP),CRPA_NUM_WAN,COREIG,NK,ISP,.TRUE.)
            ! diagonalizing the Hamiltonian for each k-point 
            ! destroys the maximally localized property of the Wannier functions
            ! also it changes the symmetry of the band structure. 
            U(:,:)=UCRPA(:,:,NK,ISP)
   
            !------------------------------------------------------------------------------
            !determine non-correlated eigenvalues 
            !UNCORR is the projector onto non-correlated system
            UNCORR=zero 
            DO I = 1, WDES%NB_TOT
               DO K = 1, WDES%NB_TOT
                  DO J = 1, CRPA_NUM_WAN
                     UNCORR(K,I)=UNCORR(K,I)-GCONJG(U(I,J))*U(K,J)
                  ENDDO         
               ENDDO              
            ENDDO          
            DO I=1,WDES%NB_TOT
               UNCORR(I,I)= UNCORR(I,I) + 1._q
            ENDDO

            CALL DIAGONALIZE_NONCORR_SYSTEM(W,UNCORR,WDES%NB_TOT-CRPA_NUM_WAN,NCOREIG,NK,NROWMIN1,NROWMAX1,ISP)
            !backup UNCORR for chi_base
            IF ( LPLOTDIS ) UPLOT(:,:) = UNCORR(:,:)

            !rearrange 
            !------|---------|----------------|
            ! nc   1   cor   2     nc         3
            ! nc  - non-correlated 
            ! cor - correlated 
            ALLEIG(1:NROWMIN1-1)=NCOREIG(1:NROWMIN1-1)
            ALLEIG(NROWMIN1:NROWMIN1+CRPA_NUM_WAN-1)=COREIG(1:CRPA_NUM_WAN)
            ALLEIG(NROWMIN1+CRPA_NUM_WAN:WDES%NB_TOT)=NCOREIG(NROWMIN1:WDES%NB_TOT-CRPA_NUM_WAN)
            !also projectors have to be set up accordingly
            !put the correlated block to 1-2
            U=UNCORR
            DO I = 1, WDES%NB_TOT
               U(I,NROWMIN1:NROWMIN1+CRPA_NUM_WAN-1)=UCRPA(I,1:CRPA_NUM_WAN,NK,ISP)
            ENDDO
            !shift the correlated block to 1-2
            UNCORR(:,:)=UCRPA(:,:,NK,ISP)
            UCRPA(:,:,NK,ISP)=zero
            DO I = 1, WDES%NB_TOT
               UCRPA(I,NROWMIN1:NROWMIN1+CRPA_NUM_WAN-1,NK,ISP)=UNCORR(I,1:CRPA_NUM_WAN)
            ENDDO
            UNCORR=U
             
            !UNCORR ->  corr + non-corr
            !UCRPA ->  corr only 
            !----------------------------------------------------------------------
            !sort eigenvalues 
            CALL MRGRNK(ALLEIG,NSORT) 
            !determine sorting matrix
            !U gives sorted eigenvalues 
            U=zero
            DO I=1,WDES%NB_TOT
               IP=NSORT(I)
               U(IP,I) = 1._q
               !identify correlated states for later purpose
               IF ( NROWMIN1 <= IP .AND. IP <= NROWMIN1+CRPA_NUM_WAN-1 ) THEN
                  LCRPA_STATE(ISP,NK,I)=.TRUE.
               ENDIF
            ENDDO 
            !all new disentangled eigenvalues are sorted and stored in NEWEIG,
            !correlated ones are in NEWEIGC
            NEWEIG(:,NK,ISP)=zero
            NEWEIGC(:,NK,ISP)=zero
            DO I = 1,WDES%NB_TOT
               DO J=1,WDES%NB_TOT
                  NEWEIG(I,NK,ISP)=NEWEIG(I,NK,ISP)+U(J,I)*ALLEIG(J)
               IF ( J >= NROWMIN1 .AND. J <= NROWMIN1+CRPA_NUM_WAN-1 ) THEN 
                     NEWEIGC(I,NK,ISP)=NEWEIGC(I,NK,ISP)+U(J,I)*ALLEIG(J)
                  ENDIF 
               ENDDO
            ENDDO
            !for MCRPA we need to rotate the HEAD appropriately
            UCRPA(:,:,NK,ISP)=zero
            DO I = 1,WDES%NB_TOT
               DO J=1,WDES%NB_TOT
                  DO K = 1, WDES%NB_TOT 
                     UCRPA(I,J, NK,ISP) = UCRPA(I,J,NK,ISP)+UNCORR(I,K)*U(K,J)
                  ENDDO
              ENDDO
            ENDDO
#ifdef debug
WRITE(500,*)'all NK=',NK
WRITE(500,'(16F12.6)')REAL(NEWEIG(:,NK,ISP),q)
DO I = 1, 16
WRITE(500,'(16F12.6)')REAL(UCRPA(I,1:16,NK,ISP),q)
ENDDO
#endif

            !first rotate correlated + noncorrelated function
            CALL SETWDES(WDES,WDES1,NK)
            WA=ELEMENTS(W,WDES1,ISP)
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF
            ! build linear combinations of the wave functions in
            ! accordance with the effective rotation matrix UCRPA
            CALL LINCOM('F',WA%CW_RED,WA%CPROJ_RED,UCRPA(1,1,NK,ISP), &
             &   WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), &
             &   WDES1%NPL_RED,WDES1%NPRO_RED,WDES%NRPLWV_RED,WDES%NPROD_RED,WDES%NB_TOT, &
             &   WA%CW_RED,WA%CPROJ_RED)
  
            ! redistribute back over bands
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF

#ifdef debug
            !first rotate correlated + noncorrelated function
            CALL SETWDES(WDES,WDES1,NK)
            WA=ELEMENTS(W,WDES1,ISP)
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF
            ! build linear combinations of the wave functions in
            ! accordance with the effective rotation matrix UCRPA
            CALL LINCOM('F',WA%CW_RED,WA%CPROJ_RED,UNCORR(1,1), &
             &   WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), &
             &   WDES1%NPL_RED,WDES1%NPRO_RED,WDES%NRPLWV_RED,WDES%NPROD_RED,WDES%NB_TOT, &
             &   WA%CW_RED,WA%CPROJ_RED)
  
            ! redistribute back over bands
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF
            !sort
            CALL SETWDES(WDES,WDES1,NK)
            WA=ELEMENTS(W,WDES1,ISP)
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF
            CALL LINCOM('F',WA%CW_RED,WA%CPROJ_RED,U(1,1), &
             &   WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), &
             &   WDES1%NPL_RED,WDES1%NPRO_RED,WDES%NRPLWV_RED,WDES%NPROD_RED,WDES%NB_TOT, &
             &   WA%CW_RED,WA%CPROJ_RED)
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF
  
            !for MCRPA we need to rotate the HEAD appropriately
            UCRPA(:,:,NK,ISP)=zero
            DO I = 1,WDES%NB_TOT
               DO J=1,WDES%NB_TOT
                  DO K = 1, WDES%NB_TOT 
                     UCRPA(I,J, NK,ISP) = UCRPA(I,J,NK,ISP)+UNCORR(I,K)*U(K,J)
                  ENDDO
              ENDDO
            ENDDO
#endif
            !maybe we want to plot the non-correlated disentangled bands
            IF (LPLOTDIS ) UCRPA(:,:,NK,ISP)=UPLOT(:,:)
!!TEST  non-correlated states should determine head
!            !for MCRPA we need to rotate the HEAD appropriately
!            UCRPA(:,:,NK,ISP)=zero
!            DO I = 1,WDES%NB_TOT
!               DO J=1,WDES%NB_TOT
!                  DO K = 1, WDES%NB_TOT 
!                     UCRPA(I,J, NK,ISP) = UCRPA(I,J,NK,ISP)+UPLOT(I,K)*U(K,J)
!                  ENDDO
!              ENDDO
!            ENDDO
!!TEST
!WRITE(600,*)'all NK=',NK
#ifdef debug
WRITE(600,'(16F12.6)')REAL(NEWEIG(:,NK,ISP),q)
DO I = 1, 16
WRITE(600,'(16F12.6)')REAL(UCRPA(I,1:16,NK,ISP),q)
ENDDO
#endif
        !--------------------------------------------------------------------------------------
         ENDDO kpts
      !----------------------------------------------------------------------------------------
      ENDDO spin
      !----------------------------------------------------------------------------------------
      DEALLOCATE(U)
      DEALLOCATE(UNCORR,ALLEIG,COREIG,NCOREIG)

      !backup old eigenvalues  
      NEWEIGC(1:WDES%NB_TOT,1:WDES%NKPTS,1:WDES%ISPIN) = &
         W%CELTOT(1:WDES%NB_TOT,1:WDES%NKPTS,1:WDES%ISPIN)
      
      !last but not least save eigenvalues 
      W%CELTOT(1:WDES%NB_TOT,1:WDES%NKPTS,1:WDES%ISPIN)= &
         NEWEIG(1:WDES%NB_TOT,1:WDES%NKPTS,1:WDES%ISPIN)
  
#ifdef debug
      !check S-orthogonality 
      CALL DUMP_S( WDES, W, INFO%LOVERL, LMDIM, CQIJ ) 
#endif 
      DO ISP = 1,WDES%ISPIN
         IF (WDES%ISPIN==2) WRITE(IO%IU6,'(/A,I1)') ' spin component ',ISP
         DO NK = 1,WDES%NKPTS
             IF ( IO%NWRITE < 3 .AND. NK > NQ_IBZ ) CYCLE
             IF( IO%IU6>=0 ) THEN
                 WRITE(IO%IU6,2202)NK,WDES%VKPT(:,NK), &
                 (I,REAL(NEWEIGC(I,NK,ISP),q),REAL( W%CELTOT(I,NK,ISP),KIND=q), &
                 W%FERTOT(I,NK,ISP)*WDES%RSPIN, &
                 LCRPA_STATE(ISP,NK,I), I=NMIN,NMAX )
             ENDIF
         ENDDO
      ENDDO
  
      DEALLOCATE(NEWEIG)
      DEALLOCATE(NEWEIGC)
  
      !dump wannier matrices for ploting later on
      IF ( LCRPAPLOT ) THEN
         IF(IO%IU0>=0)THEN
            ! write rotation matrices
            OPEN(UNIT=199,FILE='U.WAN',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
            WRITE(199) WDES%NB_TOT
            WRITE(199) UCRPA
            CLOSE(199)       
         ENDIF  
      ENDIF 
  
      ! in case RPA with disentangled wave function is selected
      ! get rid of LCRPA_STATE 
      IF ( LNLRPA ) THEN
         IF ( IO%IU0>=0) WRITE(IO%IU0,'(A)')' Disentangled fully screened interaction chosen'
         DEALLOCATE(LCRPA_STATE)
      ENDIF  

2202  FORMAT(/' k-point ',I3,' :',3X,3F10.4/ &
         &    '  band No.  old energies   dis. energies     occupation   cor. state?'/ &
         &    (3X,I4,3X,F10.4,5X,F10.4,7X,F10.5,3X,L4,4X))
  
      RETURN
   END SUBROUTINE DISENTANGLE

!**********************************************************************
! 
! This Routine separates correlated states from noncorrelated states. 
! First, the original Bloch states (W) are replaced by the eigensystem of 
! the Wannier Hamiltonian. This is necessary, since the projection is 
! not perfect for surjective maps from Bloch to the Wannier space.
! Second, the correlated states (WC) are determined using the regularized
! correlated projector. 
! 
!**********************************************************************

   SUBROUTINE SEPARATE_STATES( WDES, W , WC, LMDIM, CQIJ, LATT_CUR, & 
      INFO, KPOINTS, NQ_IBZ, IO )  
      USE base
      USE wave_high
      USE full_kpoints
      USE lattice
      USE constant
      USE mlwf!, ONLY : U_matrix,U_matrix_opt,lexclude_band,lwindow,kpt_latt,num_kpts
      USE choleski
      USE m_mrgrnk 
      USE mkpoints 
      USE ini
!      USE lie 
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      TYPE(wavespin) WC
      INTEGER LMDIM
      OVERLAP  CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      TYPE (latt) LATT_CUR
      TYPE (info_struct) INFO
      TYPE (kpoints_struct) KPOINTS
      INTEGER NQ_IBZ
      TYPE(in_struct) IO
      ! local variables
      TYPE(wavedes1) WDES1
      TYPE(wavefuna) WA
      INTEGER I,J,K,IP
      INTEGER ISP,NK
      COMPLEX(q), ALLOCATABLE :: NEWEIG(:,:,:)
      REAL(q), ALLOCATABLE ::  COREIG(:)  !correlated eigenvalues
      GDEF, ALLOCATABLE :: WINDOWFUN(:,:,:)
      GDEF, ALLOCATABLE :: U(:,:)
      GDEF, ALLOCATABLE :: UPLOT(:,:)
      INTEGER :: NMIN, NMAX
      INTEGER :: NMIN1, NMAX1
      GDEF, ALLOCATABLE :: D(:,:,:)
      INTEGER :: NK_IN_KPOINTS_FULL_ORIG
     
      CALL CHECK_FULL_KPOINTS
  
      !allocate wannier projection matrices
      ALLOCATE(U(WDES%NB_TOT,WDES%NB_TOT))
      ALLOCATE(UPLOT(WDES%NB_TOT,WDES%NB_TOT))
  
      !nullify CRPA projector and reallocate it
      IF ( ALLOCATED ( UCRPA ) ) DEALLOCATE(UCRPA) 
      ALLOCATE(UCRPA(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS, WDES%ISPIN) )
      UCRPA = zero
      IF ( ALLOCATED ( UCDER ) ) DEALLOCATE(UCDER) 
      ALLOCATE(UCDER(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS, WDES%ISPIN) )
      UCDER = zero
      !nullify weight and allocate the appropriate memory for the weighting factor
      IF (ALLOCATED(CRPA_WEIGHTS)) DEALLOCATE(CRPA_WEIGHTS)
      ALLOCATE( CRPA_WEIGHTS(WDES%NB_TOT, WDES%NKPTS, WDES%ISPIN ) )
  
      !new eigenvalues
      ALLOCATE(NEWEIG(WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN))
      NEWEIG=zero
    
      !window function will be necessary
      ALLOCATE(WINDOWFUN(WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN))
      WINDOWFUN=zero

      !for jacobi diagonalization
      ALLOCATE(D(WDES%NB_TOT,WDES%NB_TOT,WDES%ISPIN))
  
      NMIN=NCRPAMIN
      NMAX=NMIN+CRPA_NUM_WAN-1

      NMIN1=WDES%NB_TOT
      NMAX1=1
      ALLOCATE(COREIG(WDES%NB_TOT))    
  
      !----------------------------------------------------------------------------------------
      ! loop over spin and all k-points
      spin: DO ISP=1,WDES%ISPIN
      !----------------------------------------------------------------------------------------
         kpts: DO NK=1,WDES%NKPTS
         !-------------------------------------------------------------------------------------
            !initialize the identity matrix
            DO I = 1, WDES%NB_TOT
                UCRPA(I,I,NK,ISP)=one
            ENDDO
            !determine local projector U
            CALL DETERMINE_WANNIER_PROJECTOR(WDES,NK,ISP,U(:,:),IO%IU0)
            !replace the bands represented by the localized basis set
            UCRPA(1:WDES%NB_TOT,NMIN:NMAX,NK,ISP)=U(1:WDES%NB_TOT,1:CRPA_NUM_WAN)
            IF ( LCRPAPLOT ) THEN
               UPLOT = zero
               UPLOT(:,:) = UCRPA(:,:,NK,ISP)
            ENDIF

            !orthogonalize correlated sub space
            CALL ORTHOGONALIZE_SYSTEM(W,UCRPA(:,:,NK,ISP),WDES%NB_TOT,COREIG,NK,ISP)
            !backup rotation matrix 
            U(:,:)=UCRPA(:,:,NK,ISP)
            !new eigenvalues
            NEWEIG(1:WDES%NB_TOT,NK,ISP)=COREIG(1:WDES%NB_TOT)

            !backup for CDER rotation
            !the roation of the CDER needs to be done on the full mesh unfortunately  
            UCDER(:,:,NK,ISP)=UCRPA(:,:,NK,ISP)

            !first rotate all states using U
            CALL SETWDES(WDES,WDES1,NK)
            WA=ELEMENTS(W,WDES1,ISP)
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF
            ! build linear combinations of the wave functions in
            ! accordance with the effective rotation matrix U
            CALL LINCOM('F',WA%CW_RED,WA%CPROJ_RED,U(1,1), &
             &   WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), &
             &   WDES1%NPL_RED,WDES1%NPRO_RED,WDES%NRPLWV_RED,WDES%NPROD_RED,WDES%NB_TOT, &
             &   WA%CW_RED,WA%CPROJ_RED)
  
            ! redistribute back over bands
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF

            !determine Wannier projector
            CALL DETERMINE_WANNIER_PROJECTOR( WDES, NK, ISP, U, IO%IU0)
            !corresponding correlated projector
            CALL DETERMINE_CORRELATED_PROJECTOR(WDES%NB_TOT,U,UCRPA(:,:,NK,ISP),NK,ISP,IO%IU0)
            !save the correlated probabilities
            DO I=1, WDES%NB_TOT
               CRPA_WEIGHTS(I,NK,ISP)=REAL(UCRPA(I,I,NK,ISP),q)
               IF ( ABS(CRPA_WEIGHTS(I,NK,ISP))>TINYLARGE) THEN
                  NMIN1=MIN(NMIN1,I)
                  NMAX1=MAX(NMAX1,I)
               ENDIF
            ENDDO
            !diagonalize the Hamiltonian using Givens rotations (Jacobi rotation) 
            CALL JOINT_DIAGONALIZER(D,UCRPA,INFO%EDIFF,WDES%NB_TOT,NK,WDES%ISPIN,&
               IO%IU0,IO%IU6,MAXLIE,.TRUE.)
            !the diagonalized projector will contain only 1 and 0 entries 
            !This defines the Window function
            DO I=1,WDES%NB_TOT
               WINDOWFUN(I,NK,ISP)=UCRPA(I,I,NK,ISP)
               CRPA_WEIGHTS(I,NK,ISP)=WINDOWFUN(I,NK,ISP)
            ENDDO
            CALL DETERMINE_CORRELATED_PROJECTOR(WDES%NB_TOT,U,UCRPA(:,:,NK,ISP),NK,ISP,IO%IU0,&
               WINDOWFUN)

            !rotate correlated function
            CALL SETWDES(WDES,WDES1,NK)
            WA=ELEMENTS(WC,WDES1,ISP)
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, WC%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, WC%CW   (1,1,NK,ISP))
            ENDIF
            ! build linear combinations of the wave functions in
            ! accordance with the effective rotation matrix U
            CALL LINCOM('F',WA%CW_RED,WA%CPROJ_RED,UCRPA(1,1,NK,ISP), &
             &   WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), &
             &   WDES1%NPL_RED,WDES1%NPRO_RED,WDES%NRPLWV_RED,WDES%NPROD_RED,WDES%NB_TOT, &
             &   WA%CW_RED,WA%CPROJ_RED)
  
            ! redistribute back over bands
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, WC%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, WC%CW   (1,1,NK,ISP))
            ENDIF
  
        !--------------------------------------------------------------------------------------
         ENDDO kpts
      !----------------------------------------------------------------------------------------
      ENDDO spin
      !----------------------------------------------------------------------------------------
      DEALLOCATE(U)
 
      !correlated part has updated eigenvalues
      WC%CELTOT(1:WDES%NB_TOT,1:WDES%NKPTS,1:WDES%ISPIN)= &
         NEWEIG(1:WDES%NB_TOT,1:WDES%NKPTS,1:WDES%ISPIN)
      !backup old eigenvalues  
      NEWEIG(1:WDES%NB_TOT,1:WDES%NKPTS,1:WDES%ISPIN) = &
         W%CELTOT(1:WDES%NB_TOT,1:WDES%NKPTS,1:WDES%ISPIN)
      !and update old eigenvalues
      W%CELTOT(1:WDES%NB_TOT,1:WDES%NKPTS,1:WDES%ISPIN)= &
         WC%CELTOT(1:WDES%NB_TOT,1:WDES%NKPTS,1:WDES%ISPIN)
 
#ifdef debug
      !check S-orthogonality 
      CALL DUMP_S( WDES, W, INFO%LOVERL, LMDIM, CQIJ ) 
#endif 
      !update occupancies of correlated bands
      IF ( IO%IU6>=0 ) WRITE(IO%IU6,'(/A)')'Original and rotated energies:'
      DO ISP = 1,WDES%ISPIN
         IF (WDES%ISPIN==2) WRITE(IO%IU6,'(/A,I1)') ' spin component ',ISP
         DO NK = 1,WDES%NKPTS
             IF ( IO%NWRITE < 3 .AND. NK > NQ_IBZ ) CYCLE
             IF( IO%IU6>=0 ) THEN
                 WRITE(IO%IU6,2202)NK,WDES%VKPT(:,NK), &
                 (I,REAL(NEWEIG(I,NK,ISP),q),REAL( W%CELTOT(I,NK,ISP),KIND=q), &
                 ABS(W%CELTOT(I,NK,ISP)-NEWEIG(I,NK,ISP)),CRPA_WEIGHTS(I,NK,ISP),&
                 REAL(WINDOWFUN(I,NK,ISP),q),I=NMIN1,NMAX1)
             ENDIF
         ENDDO
      ENDDO
  
      DEALLOCATE(NEWEIG)
      DEALLOCATE(D,WINDOWFUN)
 
      IF (ALLOCATED(CRPA_WEIGHTS)) DEALLOCATE(CRPA_WEIGHTS)

      !only for LNLRPA UCDER is kept
      IF (.NOT. LNLRPA ) DEALLOCATE(UCDER)

      ! maybe the disentangled full RPA is wanted, 
      ! in this case we kill LCRPA_STATE,
      ! that all states are included in screending 
      IF ( LNLRPA ) THEN
         IF ( IO%IU0>=0) WRITE(IO%IU0,'(A)')' rotated fully screened interaction chosen'
      ENDIF  
  
2202  FORMAT(/' k-point ',I3,' :',3X,3F10.4/ &
         &    '  band No.  original   rotated      error  weights  window'/ &
         &    (I10,F10.4,F10.4,F11.4,F9.3,F8.3))
      RETURN
   END SUBROUTINE SEPARATE_STATES 

!**********************************************************************
!
! separates correlated states from non-correlated states by using
! regularized correlated projectors at each k-point.
! the original wavefunction is untouched. 
! this algorithm is quite stable if the projection of the Bloch to 
! Wannier state is surjective.
! 
!**********************************************************************

   SUBROUTINE SEPARATE_CORRELATED_STATES( WDES, W , WC, LMDIM, CQIJ, LATT_CUR, & 
       INFO, KPOINTS, NQ_IBZ, IO )  
      USE base
      USE wave_high
      USE full_kpoints
      USE lattice
      USE constant
      USE mlwf!, ONLY : U_matrix,U_matrix_opt,lexclude_band,lwindow,kpt_latt,num_kpts
      USE choleski
      USE m_mrgrnk 
      USE mkpoints 
!      USE lie
      USE ini
      IMPLICIT NONE
      TYPE(wavedes) WDES
      TYPE(wavespin) W
      TYPE(wavespin) WC
      INTEGER LMDIM
      OVERLAP  CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      TYPE (latt) LATT_CUR
      TYPE (info_struct) INFO
      TYPE (kpoints_struct) KPOINTS
      INTEGER         NQ_IBZ
      TYPE(in_struct) IO
      ! local variables
      TYPE(wavedes1) WDES1
      TYPE(wavefuna) WA
      INTEGER I 
      INTEGER ISP,NK
      GDEF, ALLOCATABLE :: U(:,:)
      GDEF, ALLOCATABLE :: WINDOWFUN(:,:,:)
      GDEF, ALLOCATABLE :: D(:,:,:)
      INTEGER :: NMIN, NMAX
     
      PROFILING_START( 'separate_correlated_states' )

      CALL CHECK_FULL_KPOINTS
  
      !allocate wannier projection matrices
      ALLOCATE(U(WDES%NB_TOT,WDES%NB_TOT))
      ALLOCATE(WINDOWFUN(WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN))
      WINDOWFUN=zero
  
      !nullify CRPA projector and reallocate it
      IF ( ALLOCATED ( UCRPA ) ) DEALLOCATE(UCRPA) 
      ALLOCATE(UCRPA(WDES%NB_TOT,WDES%NB_TOT,WDES%NKPTS, WDES%ISPIN) )
      UCRPA = 0
      !nullify weight and allocate the appropriate memory for the weighting factor
      IF (ALLOCATED(CRPA_WEIGHTS)) DEALLOCATE(CRPA_WEIGHTS)
      ALLOCATE( CRPA_WEIGHTS(WDES%NB_TOT, WDES%NKPTS, WDES%ISPIN ) )
  
#ifdef debug
      CALL FIND_ON_LIE_MANIFOLD(UCRPA(:,:,1,1),4,'MIN',NTARGET_MASK,INFO%EDIFF,MAXLIE,1,IO,&
         .FALSE.,.TRUE.) 
      CALL vtutor%stopCode()
#endif 
  
      ALLOCATE(D(WDES%NB_TOT,WDES%NB_TOT,WDES%ISPIN))
  
      NMIN=WDES%NB_TOT
      NMAX=1
      DO ISP=1,WDES%ISPIN
         DO NK=1, WDES%NKPTS
            !determine Wannier projector
            CALL DETERMINE_WANNIER_PROJECTOR( WDES, NK, ISP, U, IO%IU0)
            !corresponding correlated projector
            CALL DETERMINE_CORRELATED_PROJECTOR(WDES%NB_TOT,U,UCRPA(:,:,NK,ISP),NK,ISP,IO%IU0)
            !save the correlated probabilities
            DO I=1, WDES%NB_TOT
               CRPA_WEIGHTS(I,NK,ISP)=REAL(UCRPA(I,I,NK,ISP),q)
            ENDDO
            !diagonalize the Hamiltonian using Givens rotations (Jacobi rotation) 
            CALL JOINT_DIAGONALIZER(D,UCRPA,INFO%EDIFF,WDES%NB_TOT,NK,WDES%ISPIN,&
               IO%IU0,IO%IU6,MAXLIE,.TRUE.)
            !the diagonalized projector will contain only 1 and 0 entries 
            !This defines the Window function
            DO I=1,WDES%NB_TOT
               WINDOWFUN(I,NK,ISP)=UCRPA(I,I,NK,ISP)
               CRPA_WEIGHTS(I,NK,ISP)=WINDOWFUN(I,NK,ISP)
               IF ( ABS(WINDOWFUN(I,NK,ISP)) > 0 )THEN
                  NMIN=MIN(NMIN,I)
                  NMAX=MAX(NMAX,I)
               ENDIF
            ENDDO
            !determine again the correlated projector, this time multiplied with the Window function
            CALL DETERMINE_CORRELATED_PROJECTOR(WDES%NB_TOT,U,UCRPA(:,:,NK,ISP),NK,ISP,IO%IU0,&
               WINDOWFUN)
  
            !rotate the Wavefunction 
            CALL SETWDES(WDES,WDES1,NK)
            WA=ELEMENTS(WC,WDES1,ISP)
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, WC%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, WC%CW   (1,1,NK,ISP))
            ENDIF
            ! build linear combinations of the wave functions in
            ! accordance with the effective rotation matrix U
            CALL LINCOM('F',WA%CW_RED,WA%CPROJ_RED,UCRPA(1,1,NK,ISP), &
             &   WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), &
             &   WDES1%NPL_RED,WDES1%NPRO_RED,WDES%NRPLWV_RED,WDES%NPROD_RED,WDES%NB_TOT, &
             &   WA%CW_RED,WA%CPROJ_RED)
            ! redistribute back over bands
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, WC%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, WC%CW   (1,1,NK,ISP))
            ENDIF
         ENDDO
      ENDDO
      !get rid of auxillary matrix 
      DEALLOCATE(U)

      !relevant bands for correlated Green's function
      NCRPAMIN=NMIN 
      NCRPAMAX=NMAX 

      !dump some information 
      DO ISP = 1,WDES%ISPIN
         IF (WDES%ISPIN==2 .AND. IO%IU6>= 0 ) WRITE(IO%IU6,'(/A,I1)') ' spin component ',ISP
         DO NK = 1,WDES%NKPTS
             !per default print only irreducible wedge of Brillouin zone
             IF ( IO%NWRITE < 3 .AND. NK > NQ_IBZ ) CYCLE
             IF( IO%IU6>=0 ) THEN
                 WRITE(IO%IU6,2202)NK,WDES%VKPT(:,NK), &
                 (I,REAL(W%CELTOT(I,NK,ISP),q),REAL( WC%CELTOT(I,NK,ISP),KIND=q), &
                 REAL(WINDOWFUN(I,NK,ISP),q), CRPA_WEIGHTS(I,NK,ISP), &
                  I=NMIN,NMAX )
             ENDIF
         ENDDO
      ENDDO
  
      DEALLOCATE(CRPA_WEIGHTS)
  
2202  FORMAT(/' k-point ',I3,' :',3X,3F10.4/ &
         &    '  band No.  old energies   rot. energies     windowfun    weights'/ &
         &    (3X,I4,3X,F10.4,5X,F10.4,7X,F10.5,3X,F10.5,4X))
  
      PROFILING_STOP( 'separate_correlated_states' )
      RETURN
   END SUBROUTINE SEPARATE_CORRELATED_STATES

!**********************************************************************
!
! determines CRPA_WEIGHTS for weighted approach of cRPA 
!
!**********************************************************************

   SUBROUTINE DETERMINE_CRPA_WEIGHTS( WDES, W, NQ_IBZ, IO ) 
      USE base 
      USE wave_struct_def
      USE full_kpoints
      IMPLICIT NONE
      TYPE(wavedes)   WDES
      TYPE(wavespin)  W
      INTEGER         NQ_IBZ
      TYPE(in_struct) IO
      ! local variables
      INTEGER IU0
      INTEGER IU6
      INTEGER I 
      INTEGER ISP,NK
      INTEGER NMIN,NMAX
      GDEF, ALLOCATABLE    :: U(:,:)
      GDEF, ALLOCATABLE    :: UP(:,:)
      REAL(q) :: RN 
 
      PROFILING_START( 'determine_crpa_weights' )
 
      CALL CHECK_FULL_KPOINTS
      IU0=IO%IU0
      IU6=IO%IU6
 
      !allocate memory for effective rotation matrices
      ALLOCATE(U(WDES%NB_TOT,WDES%NB_TOT))
      ALLOCATE(UP(WDES%NB_TOT,WDES%NB_TOT))
    
      !FCRPA 
      !nullify weight and allocate the appropriate memory for the weighting factor
      IF (ALLOCATED(CRPA_WEIGHTS)) DEALLOCATE(CRPA_WEIGHTS)
      ALLOCATE( CRPA_WEIGHTS(WDES%NB_TOT, WDES%NKPTS, WDES%ISPIN ) )
 
      NMIN=WDES%NB_TOT
      NMAX=1
      !----------------------------------------------------------------------------------------
      spin: DO ISP=1,WDES%ISPIN  !spin
      !----------------------------------------------------------------------------------------
         kpoints: DO NK=1,WDES%NKPTS  !k-points
         !-------------------------------------------------------------------------------------
            !get the CRPA weights 
            CALL DETERMINE_WANNIER_PROJECTOR( WDES, NK, ISP, U, IU0)
            CALL DETERMINE_CORRELATED_PROJECTOR(WDES%NB_TOT,U,UP,NK,ISP,IU0)
  
            !now U is the rotation matrix at kpoint NK with spin ISP
            !finally calculate the weight for CRPA
            DO I = 1, WDES%NB_TOT
               CRPA_WEIGHTS(I,NK,ISP) = REAL(UP(I,I),q)
               IF ( REAL(UP(I,I),q)>TINYLARGE ) THEN
                  NMIN = MIN(NMIN,I)
                  NMAX = MAX(NMAX,I)
               ENDIF 
            ENDDO              
         !-------------------------------------------------------------------------------------
         ENDDO kpoints
      !----------------------------------------------------------------------------------------
      ENDDO spin
      !----------------------------------------------------------------------------------------
      DEALLOCATE(U)
      DEALLOCATE(UP)

      !set relevant bands for correlated Green's function
      NCRPAMIN=NMIN
      NCRPAMAX=NMAX

      !dump weights
      DO ISP=1,WDES%ISPIN  !spin
         IF (WDES%ISPIN==2) THEN
            IF ( IU6>=0 ) WRITE(IU6,'(/A,I1)') ' spin component ',ISP
         ENDIF 
         DO NK=1,WDES%NKPTS  !k-points
            !default output is the irreducible wedge of the BZ
            IF ( IO%NWRITE < 3 .AND. NK > NQ_IBZ ) CYCLE
            IF( IU6>=0  ) THEN
               WRITE(IU6,1000)NK,WDES%VKPT(:,NK), &
               (I,REAL( W%CELTOT(I,NK,ISP),q), CRPA_WEIGHTS(I,NK,ISP),&
               W%FERTOT(I,NK,ISP)*WDES%RSPIN, I=NMIN,NMAX)
            ENDIF
            !also # of states at each k-point sound be number of states J
            RN = SUM(CRPA_WEIGHTS(1:WDES%NB_TOT,NK,ISP))
            IF ( ABS(RN) < TINYLARGE ) THEN
               IF (IU0>=0) WRITE(IU0,'(A,E13.6,A,I3,A,2I4)')' WARNING: # of target states is ',RN,&
                  ' and not ',SUM(NTARGET_MASK(1:CRPA_NUM_WAN)),' for k-point,ISP',NK,ISP
            ENDIF
         ENDDO
      ENDDO 
 
      PROFILING_STOP( 'determine_crpa_weights' )
      RETURN
1000  FORMAT(/' k-point ',I3,' :',3X,3F10.4/ &
      &         '  band No.     energies   CRPA-weights   occupation'/ &
      &           (3X,I4,3X,F10.4,3X,F12.4,3X,F10.5))
   END SUBROUTINE DETERMINE_CRPA_WEIGHTS
 
!**********************************************************************
!
! Determine Wannier transformation matrix
!
!**********************************************************************

   SUBROUTINE  DETERMINE_WANNIER_PROJECTOR(WDES,NK,ISP,U,IU0)
      USE prec
      USE wave
      USE mlwf
      USE locproj
      USE full_kpoints
      IMPLICIT NONE
      TYPE(wavedes) WDES
      INTEGER       NK
      INTEGER       ISP
      GDEF          U(:,:) 
      INTEGER       IU0
      !local
      INTEGER I,J,K,L,IP,IPP
      INTEGER       NKP 
      INTEGER NW,NB,NBEXCL,NBwin
      INTEGER ISPINOR
      GDEF CTMP

      PROFILING_START('determine_wannier_projector') 

      ! number of wannier function must be defined
      IF ( CRPA_NUM_WAN < 1 ) THEN
         CALL vtutor%error("DETERMINE_WANNIER_PROJECTOR reports CRPA_NUM_WAN not set")
      ENDIF

      ! if projector is read from file, do this and get out
      IF ( LUFROMFILE ) THEN
         !read projector from file 
         IF ( .NOT. ALLOCATED(UWAN) )THEN
            CALL vtutor%error("DETERMINE_WANNIER_PROJECTOR ERROR: UWAN not allocated")
         ENDIF
         ! check band consistency 
         IF ( SIZE(UWAN,1) < WDES%NB_TOT .OR. SIZE(UWAN,2) < CRPA_NUM_WAN )THEN
            CALL vtutor%error("DETERMINE_WANNIER_PROJECTOR reports inconsistent bands in UWAN " // &
               str(SIZE(UWAN,1)) // "<" // str(WDES%NB_TOT) // " or " // str(SIZE(UWAN,2)) // "<" // str(WDES%NB_TOT))
         ENDIF
         ! check k-point and spin consistency 
         IF ( SIZE(UWAN,3) < NK .OR. SIZE(UWAN,4) < ISP )THEN
            CALL vtutor%error("DETERMINE_WANNIER_PROJECTOR reports inconsistent k,spin in UWAN " // &
               str(SIZE(UWAN,3)) // " " // str(SIZE(UWAN,4)) // " " // str(NK) // " " // str(ISP))
         ENDIF
         U = 0
         U(1:WDES%NB_TOT,1:CRPA_NUM_WAN)=UWAN(1:WDES%NB_TOT,1:CRPA_NUM_WAN,NK,ISP)
      ENDIF

      !depending on the basis choice the transformation matrix
      !needs to be set up appropriately
      !----------------------------------------------------------------
      IF( LUSEMLWF ) THEN   ! maximally localized Wannier functions
      !----------------------------------------------------------------
         ! find the corresponding entry in KPOINTS_FULL, that contains
         ! the set of k-points used in the computation of the rotation matrices
         DO NKP=1,KPOINTS_FULL%NKPTS
            IF (LIDENTICAL_KPOINT(WDES%VKPT(:,NK),KPOINTS_FULL%VKPT(:,NKP))) exit
         ENDDO
         IF (NKP>KPOINTS_FULL%NKPTS) THEN
            CALL vtutor%error("DETERMINE_WANNIER_PROJECTOR ERROR: no matching k-point found in &
               &kpt_latt " // str(NK))
         ENDIF

         ! setup the effective rotation matrix U at the present k-point
         CALL MLWF_GET_U_IK_ISP(MLWF_GLOBAL,NKP,ISP,U)
      !----------------------------------------------------------------
      ELSEIF ( LUSELOCPROJ ) THEN
      !----------------------------------------------------------------
         !the transformation matrix must be allocated for this basis choice
         IF( .NOT.ALLOCATED(LPRJ_COVL) ) THEN
            CALL vtutor%error("DETERMINE_WANNIER_PROJECTOR ERROR: LPRJ_COVL not allocated!")
         ENDIF 

         ! sizes must be consistent 
         IF ( SIZE( LPRJ_COVL , 1 ) < WDES%NB_TOT  ) THEN
            CALL vtutor%error("DETERMINE_WANNIER_PROJECTOR ERROR: LPRJ_COVL has wrong number of bands &
               &" // str(WDES%NB_TOT) // " " // str(SIZE(LPRJ_COVL,1)))
         ENDIF

         ! sizes must be consistent 
         IF ( SIZE( LPRJ_COVL , 4 ) < CRPA_NUM_WAN  ) THEN
            CALL vtutor%error("DETERMINE_WANNIER_PROJECTOR ERROR: LPRJ_COVL has wrong number of &
               &Wannier functions " // str(CRPA_NUM_WAN) // " " // str(SIZE(LPRJ_COVL,4)))
         ENDIF
   
         ! find the corresponding entry in KPOINTS_FULL, that contains
         ! the set of k-points used in the computation of the rotation matrices
         DO NKP=1,KPOINTS_FULL%NKPTS
            IF (LIDENTICAL_KPOINT(WDES%VKPT(:,NK),KPOINTS_FULL%VKPT(:,NKP))) exit
         ENDDO
         IF (NKP>KPOINTS_FULL%NKPTS) THEN
            CALL vtutor%error("DETERMINE_WANNIER_PROJECTOR ERROR: no matching k-point found in &
               &KPOINTS_FULL " // str(NK))
         ENDIF

         U=zero
         DO ISPINOR = 1, WDES%NRSPINORS  
            DO I = 1, WDES%NB_TOT        
               DO J = 1, CRPA_NUM_WAN
!               IF (.NOT.WNPR_use_band(I,NKP,ISP)) CYCLE
               U(I,J)=LPRJ_COVL( I, NKP, ISP+ISPINOR-1,J)
               ENDDO
            ENDDO
         ENDDO

      !----------------------------------------------------------------
      ELSEIF ( LUSEBLOCH ) THEN
      !----------------------------------------------------------------
         !maybe to most simple tranformation :)
         U=zero
         DO I=1,WDES%NB_TOT
            U(I,I) = 1._q
         ENDDO
      !----------------------------------------------------------------
      ENDIF 
      !----------------------------------------------------------------
      PROFILING_STOP( 'determine_wannier_projector') 
   END SUBROUTINE DETERMINE_WANNIER_PROJECTOR

!**********************************************************************
!
! Determine correlated projector P
!
!**********************************************************************

   SUBROUTINE  DETERMINE_CORRELATED_PROJECTOR(NB_TOT,U,P,NK,ISP,IU0,W)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INTEGER       NB_TOT
      GDEF          U(:,:) 
      GDEF          P(:,:) 
      INTEGER       IU0,NK,ISP
      GDEF,OPTIONAL :: W(:,:,:) !window function
      !local
      INTEGER I,J,K
 
      ! test will be necessary 
      IF ( SIZE( U,1) /= SIZE(P,1) .AND. SIZE(P,1) < NB_TOT ) THEN
         CALL vtutor%bug("internal error in VASP from DETERMINE_CORRELATED_PROJECTOR: Size 1 of U &
            &and P inconsistent " // str(SIZE(U,1)) // " " // str(SIZE(P,1)) // " " // str(NB_TOT), &
            __FILE__, __LINE__)
      ENDIF 
      IF ( SIZE( U,2) /= SIZE(P,2) .AND. SIZE(P,2) < NB_TOT ) THEN
         CALL vtutor%bug("internal error in VASP from DETERMINE_CORRELATED_PROJECTOR: Size 2 of U &
            &and P inconsistent " // str(SIZE(U,2)) // " " // str(SIZE(P,2)) // " " // str(NB_TOT), &
            __FILE__, __LINE__)
      ENDIF 
 
      P(:,:)=zero 
      !determine correlated projector U 
      DO I = 1, NB_TOT
         DO K = 1, NB_TOT
            !only correlated states contribute to projector 
            DO J = 1, CRPA_NUM_WAN
               IF ( NTARGET_MASK(J) == 1) THEN
                  P(K,I)=P(K,I)+GCONJG(U(I,J))*U(K,J)
               ENDIF
            ENDDO         
         ENDDO              
      ENDDO    
 
      IF ( ALLOCATED(CRPA_WEIGHTS))THEN
        DO I = 1, NB_TOT
           CRPA_WEIGHTS(I,NK,ISP)=REAL(P(I,I),q)
        ENDDO
      ENDIF
      !if Window function is present we multiply the projector with it
      IF (PRESENT(W) .AND. LREGULARIZE) THEN
         !check if sizes are correct
         IF ( SIZE(W,1) /= NB_TOT ) THEN
            CALL vtutor%bug("internal error in VASP from DETERMINE_CORRELATED_PROJECTOR: Size 1 of " &
               // "W and U inconsistent " // str(SIZE(W,1)) // " " // str(NB_TOT), __FILE__, __LINE__)
         ENDIF 
         !check if sizes are correct
         IF ( SIZE(W,2) < NK ) THEN
            CALL vtutor%bug("internal error in VASP from DETERMINE_CORRELATED_PROJECTOR: Size 2 of " &
               // "W and U inconsistent " // str(SIZE(W,2)) // " " // str(NK), __FILE__, __LINE__)
         ENDIF 
         !check if sizes are correct
         IF ( SIZE(W,3) < ISP) THEN
            CALL vtutor%bug("internal error in VASP from DETERMINE_CORRELATED_PROJECTOR: Size 3 of " &
               // "W and U inconsistent " // str(SIZE(W,3)) // " " // str(ISP), __FILE__, __LINE__)
         ENDIF 
 
         DO I = 1, NB_TOT
            P(:,I)=P(:,I)*W(I,NK,ISP)
         ENDDO              
#ifdef verbose 
         IF ( IU0>=0 ) THEN
            WRITE(*,'(A,I4)')'******************* WINDOW *********** NK=',NK
            CALL DUMP_PROJ(P,NB_TOT,3, 17, 33,.TRUE.)
         ENDIF 
      ELSE 
         IF ( IU0>=0 ) THEN
            WRITE(*,'(A,I4)')'++++++++++++++++++ original ++++++++++ NK=',NK
            CALL DUMP_PROJ(P,NB_TOT,3, 17, 33, .TRUE.)
          ENDIF 
#endif
      ENDIF
 
 
   END SUBROUTINE DETERMINE_CORRELATED_PROJECTOR

!**********************************************************************
!
! orthogonalizes non-correlated sub system
! 
!**********************************************************************

  SUBROUTINE DIAGONALIZE_NONCORR_SYSTEM(W, U, NCORR, EIG , NK,NCOLSTART,NCOLSTOP, ISP)
      USE base
      USE wave_high
      USE constant
      IMPLICIT NONE
      TYPE(wavespin) :: W
      GDEF           :: U(:,:)
      INTEGER        :: NCORR
      REAL(q)        :: EIG(:)
      INTEGER        :: NK , ISP
      INTEGER        :: NCOLSTART,NCOLSTOP
      !local
      INTEGER                :: NB_TOT
      INTEGER                :: I, J, K
      GDEF, ALLOCATABLE      :: S(:,:)
      GDEF, ALLOCATABLE      :: H(:,:), HP(:,:)
      REAL(q), ALLOCATABLE   :: SIGMA(:)
      INTEGER                :: NC
      REAL(q)                :: RTMP
      LOGICAL, ALLOCATABLE   :: LTAG(:)
     
      !the total number of bands is NBANDS
      NB_TOT = W%WDES%NB_TOT
 
      !check the size of EIG needs to be OK
      IF ( SIZE(EIG)/=NB_TOT ) THEN
         CALL vtutor%error("ERROR: EIG not large enough " // str(SIZE(EIG)) // " " // str(NB_TOT))
      ENDIF 
      !also check if U is allocated and has proper size 
      IF ( SIZE(U,1)/=NB_TOT .OR. SIZE(U,2)/=NB_TOT ) THEN
         CALL vtutor%error("ERROR: U not large enough " // str(SIZE(U,1)) // " " // str(SIZE(U,2)) // &
            " " // str(NB_TOT))
      ENDIF 

      !allocate auxillary matrices 
      ALLOCATE(S(NB_TOT,NB_TOT))
      ALLOCATE(H(NB_TOT,NB_TOT))
      ALLOCATE(HP(NB_TOT,NB_TOT))
      ALLOCATE(SIGMA(NB_TOT))
      S(:,:)=U(:,:)
      H=zero
      SIGMA=0
 
      ! build correlated hamiltonian
      H=zero
      DO I = 1, NB_TOT
         DO J = 1, NB_TOT
             DO K = 1, NB_TOT
                H(I,J) = H(I,J) + U(K, J)*REAL(W%CELTOT(K,NK,ISP),q)*GCONJG(U(K,I))
             ENDDO
         ENDDO
      ENDDO

      !now diagonalize Hamiltonian 
      CALL DIAGH( H, NB_TOT, EIG)
      ! check result
      ! find number of non-zero eigenvalues (should be NCORR eigenvalues )
      NC = 0
      SIGMA=0
      SIGMA=EIG
      EIG=0  
      DO I = 1, NB_TOT
 
         IF ( ABS(SIGMA(I))> TINYLARGE) THEN
           NC = NC + 1 
           EIG(NC)=SIGMA(I)
         ENDIF
 
      ENDDO
      IF ( NC /= NCORR ) THEN
         CALL vtutor%bug("internal error in VASP from DIAGONALIZE_NONCORR_SYSTEM # of corr. states &
            &inconsistent " // str(NC) // " " // str(NCORR) // " " // str(NK), &
            __FILE__, __LINE__)
      ENDIF 
 
      !compute effective rotation matrix
      U=zero
#ifdef gammareal
      CALL DGEMM( 'N', 'N', NB_TOT, NB_TOT, NB_TOT, 1._q, S(1,1), NB_TOT,&
         H(1,1), NB_TOT, 0._q, U(1,1), NB_TOT)
#else
      CALL ZGEMM( 'N', 'N', NB_TOT, NB_TOT, NB_TOT, 1._q, S(1,1), NB_TOT,&
         H(1,1), NB_TOT, 0._q, U(1,1), NB_TOT)
#endif 

      !resulting matrix has a block of mostly zero entries between NCOLSTART and NCOLSTOP
      !these band indices are determined by the energy window in the Wannier projection
      !there are only NCOLSTOP-NCOLSTART+1-DMFT_num_states non-zero states in this block 
      !these columns have to be identified and shifted to the right in this block 
      H=zero
      H=U
      ALLOCATE(LTAG(NB_TOT))
      LTAG(:)=.FALSE. 
      !scan relevant block for non-zero elements 
!WRITE(*,*)NK,NB_TOT
      DO J = NCOLSTART, NCOLSTOP 
         RTMP=0
         DO I = 1, NB_TOT
            RTMP = RTMP + ABS(H(I,J))
            !RTMP = RTMP + REAL(H(I,J),q)
         ENDDO 
         !tag column for non-zero entries 
         IF ( RTMP > 1.E-12_q ) LTAG(J)=.TRUE. 
! WRITE(*,'(I4,L4,F12.6," for ",16F10.4)')J,LTAG(J),RTMP,REAL(H(:,J))
      ENDDO
     
      S=H
      S(:,NCOLSTART:NCOLSTOP)=0
      J=0
      DO I = NCOLSTART, NCOLSTOP
         IF ( LTAG(I)) THEN
            J=J+1
            S(:,NCOLSTOP-J+1)=H(:,I)
         ENDIF
      ENDDO

#ifdef debug
WRITE(100,*)'non-corr NK=',NK,NCOLSTART,NCOLSTOP
WRITE(100,'(16L12)')LTAG(1:16)
DO I = 1, NB_TOT
WRITE(100,'(16F12.6)')REAL(U(I,1:16),q)
ENDDO
WRITE(100,*)'IMAG'
DO I = 1, NB_TOT
WRITE(100,'(16F12.6)')IMAG(U(I,1:16))
ENDDO
      
WRITE(100,*)'rearranged for NK=',NK,NCOLSTART,NCOLSTOP
WRITE(100,'(16L12)')LTAG(1:16)
WRITE(100,'(16F12.6)')EIG(1:16)
DO I = 1, NB_TOT
WRITE(100,'(16F12.6)')REAL(S(I,1:16),q)
ENDDO
WRITE(100,*)'IMAG'
DO I = 1, NB_TOT
WRITE(100,'(16F12.6)')IMAG(S(I,1:16))
ENDDO
#endif

      !restore 
      U=S
      !get rid of auxilary arrays 
      DEALLOCATE( S ) 
      DEALLOCATE( H ) 
      DEALLOCATE( HP ) 
      DEALLOCATE( LTAG )

      RETURN
   END SUBROUTINE DIAGONALIZE_NONCORR_SYSTEM 

  SUBROUTINE ORTHOGONALIZE_SYSTEM(W, U, NCORR, EIG , NK, ISP)
      USE base
      USE wave_high
      USE constant
      IMPLICIT NONE
      TYPE(wavespin) :: W
      GDEF           :: U(:,:)
      INTEGER        :: NCORR
      REAL(q)        :: EIG(:)
      INTEGER        :: NK , ISP
      !local
      INTEGER                :: NB_TOT
      INTEGER                :: I, J, K
      GDEF, ALLOCATABLE      :: S(:,:)
      GDEF, ALLOCATABLE      :: H(:,:), HP(:,:)
      REAL(q), ALLOCATABLE   :: SIGMA(:)
      INTEGER                :: NC
     
      !the total number of bands is NBANDS
      NB_TOT = W%WDES%NB_TOT
 
      !check the size of EIG needs to be OK
      IF ( SIZE(EIG)/=NB_TOT ) THEN
         CALL vtutor%error("ERROR: EIG not large enough " // str(SIZE(EIG)) // " " // str(NB_TOT))
      ENDIF 
      !also check if U is allocated and has proper size 
      IF ( SIZE(U,1)/=NB_TOT .OR. SIZE(U,2)/=NB_TOT ) THEN
         CALL vtutor%error("ERROR: U not large enough " // str(SIZE(U,1)) // " " // str(SIZE(U,2)) // &
            " " // str(NB_TOT))
      ENDIF 
 
      !allocate auxillary matrices 
      ALLOCATE(S(NB_TOT,NB_TOT))
      ALLOCATE(H(NB_TOT,NB_TOT))
      ALLOCATE(HP(NB_TOT,NB_TOT))
      ALLOCATE(SIGMA(NB_TOT))
      S(:,:)=U(:,:)
      H=zero
      SIGMA=0
 
      ! build correlated hamiltonian
      H=zero
      DO I = 1, NB_TOT
         DO J = 1, NB_TOT
             DO K = 1, NB_TOT
                H(I,J) = H(I,J) + U(K, J)*W%CELTOT(K,NK,ISP)*GCONJG(U(K,I))
             ENDDO
         ENDDO
      ENDDO
 
      !now diagonalize Hamiltonian 
      CALL DIAGH( H, NB_TOT, EIG)
      ! check result
      ! find number of non-zero eigenvalues (should be NCORR eigenvalues )
      NC = 0
      SIGMA=0
      SIGMA=EIG
      EIG=0  
      DO I = 1, NB_TOT
         EIG(I)=SIGMA(I)
      ENDDO
 
      !compute effective rotation matrix
      U=zero
#ifdef gammareal
      CALL DGEMM( 'N', 'N', NB_TOT, NB_TOT, NB_TOT, 1._q, S(1,1), NB_TOT,&
         H(1,1), NB_TOT, 0._q, U(1,1), NB_TOT)
#else
      CALL ZGEMM( 'N', 'N', NB_TOT, NB_TOT, NB_TOT, 1._q, S(1,1), NB_TOT,&
         H(1,1), NB_TOT, 0._q, U(1,1), NB_TOT)
#endif 
 
      !now we can get rid of the auxilary matrices S and H
      DEALLOCATE( S ) 
      DEALLOCATE( H ) 
      DEALLOCATE( HP ) 
  
      RETURN
   END SUBROUTINE ORTHOGONALIZE_SYSTEM 

!**********************************************************************
!
! orthogonalizes correlated sub system
! 
!**********************************************************************

   SUBROUTINE DIAGONALIZE_CORR_SYSTEM(W, U, NCORR, EIG , NK, ISP, LORTHOGONALIZE)
      USE base
      USE wave_high
      USE constant
      IMPLICIT NONE
      TYPE(wavespin) :: W
      GDEF           :: U(:,:)
      INTEGER        :: NCORR
      REAL(q)        :: EIG(:)
      INTEGER        :: NK , ISP
      LOGICAL        :: LORTHOGONALIZE
      ! local
      INTEGER                :: NB_TOT, NB_SUB
      INTEGER                :: I, J, K
      GDEF, ALLOCATABLE      :: S(:,:)
      GDEF, ALLOCATABLE      :: H(:,:), HP(:,:)
      REAL(q), ALLOCATABLE   :: SIGMA(:)
      INTEGER                :: NC
     
      !the total number of bands is NBANDS
      NB_TOT = W%WDES%NB_TOT
      !number of states in the subspace
      NB_SUB = NCORR
 
      !check the size of EIG needs to be OK
      IF ( SIZE(EIG)/=NB_TOT ) THEN
         CALL vtutor%error("ERROR: EIG not large enough " // str(SIZE(EIG)) // " " // str(NB_TOT))
      ENDIF 
      !also check if U is allocated and has proper size 
      IF ( SIZE(U,1)/=NB_TOT .OR. SIZE(U,2)/=NB_TOT ) THEN
         CALL vtutor%error("ERROR: U not large enough " // str(SIZE(U,1)) // " " // str(SIZE(U,2)) // &
            " " // str(NB_TOT))
      ENDIF 
 
      !allocate auxillary matrices 
      ALLOCATE(S(NB_TOT,NB_TOT))
      ALLOCATE(H(NB_TOT,NB_TOT))
      ALLOCATE(HP(NB_TOT,NB_TOT))
      ALLOCATE(SIGMA(NB_TOT))
      S(:,:)=U(:,:)
      H=zero
      SIGMA=0
 
      ! build correlated hamiltonian
      H=zero
      DO I = 1, NB_TOT
         DO J = 1, NB_TOT
             DO K = 1, NB_TOT
                H(I,J) = H(I,J) + U(K, J)*REAL(W%CELTOT(K,NK,ISP),q)*GCONJG(U(K,I))
             ENDDO
         ENDDO
      ENDDO
       
      !now diagonalize Hamiltonian 
      CALL DIAGH( H, NB_SUB, EIG)
      ! check result
      ! find number of non-zero eigenvalues (should be NCORR eigenvalues )
      NC = 0
      SIGMA=0
      SIGMA=EIG
      EIG=0  
      DO I = 1, NB_SUB
         IF ( ABS(SIGMA(I))>TINYLARGE ) THEN
           NC = NC + 1 
           EIG(NC)=SIGMA(I)
         ENDIF
      ENDDO
      IF ( NC /= NCORR ) THEN
         CALL vtutor%bug("internal error in VASP from DIAGONALIZE_CORR_SYSTEM # of corr. states &
            &inconsistent " // str(NC) // " " // str(NCORR) // " " // str(NK), & 
            __FILE__, __LINE__)
      ENDIF 
      !compute effective rotation matrix
      U=zero
#ifdef gammareal
      CALL DGEMM( 'N', 'N', NB_TOT, NB_TOT, NB_TOT, 1._q, S(1,1), NB_TOT,&
         H(1,1), NB_TOT, 0._q, U(1,1), NB_TOT)
#else
      CALL ZGEMM( 'N', 'N', NB_TOT, NB_TOT, NB_TOT, 1._q, S(1,1), NB_TOT,&
         H(1,1), NB_TOT, 0._q, U(1,1), NB_TOT)
#endif 
 

      IF ( LORTHOGONALIZE ) THEN
         S=U
         H=U
         HP=zero
         !build overlap
#ifdef gammareal
         CALL DGEMM( 'T', 'N', NB_TOT, NB_TOT, NB_TOT, 1._q, S(1,1), NB_TOT,&
            H(1,1), NB_TOT, 0._q, HP(1,1), NB_TOT)
#else
         CALL ZGEMM( 'C', 'N', NB_TOT, NB_TOT, NB_TOT, 1._q, S(1,1), NB_TOT,&
            H(1,1), NB_TOT, 0._q, HP(1,1), NB_TOT)
#endif 

#ifdef debug
WRITE(100,*)'overlap NK=',NK
DO I = 1, 16
WRITE(100,'(16F12.6)')REAL(HP(I,1:16),q)
ENDDO
WRITE(100,*)'IMAG'
DO I = 1, 16
WRITE(100,'(16F12.6)')IMAG(HP(I,1:16))
ENDDO
#endif

      ENDIF

      !now we can get rid of the auxilary matrices S and H
      DEALLOCATE( S ) 
      DEALLOCATE( H ) 
      DEALLOCATE( HP ) 

      RETURN
   END SUBROUTINE DIAGONALIZE_CORR_SYSTEM 

!**********************************************************************
!
! Diagonalizes a matrix CUNI and returns its eigenvalues in HFEIG
!
!**********************************************************************

   SUBROUTINE DIAGH(CUNI, NBANDS, HFEIG )
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      GDEF    :: CUNI(:,:)
      INTEGER    :: NBANDS
      REAL(q) :: HFEIG(NBANDS)
      INTEGER :: NEIG
    ! local
      INTEGER, PARAMETER :: LWORK=32
      GDEF    ::CWRK(LWORK*NBANDS)
      REAL(q) :: W(3*NBANDS)
      REAL(q):: EV(NBANDS,NBANDS)
      INTEGER :: IFAIL, NDIM, N
      NDIM = SIZE(CUNI,1)

#ifdef gammareal
         CALL DSYEV &
           ('V','U',NBANDS,CUNI(1,1),NDIM,HFEIG,CWRK,LWORK*NBANDS, &
           W,  IFAIL)
      IF (IFAIL/=0) THEN
         CALL vtutor%bug("internal error in VASP from DIAGH: Call to routine DSYEV failed! Error " &
            // "code was " // str(IFAIL), __FILE__, __LINE__)
#else
         CALL ZHEEV &
           ('V','U',NBANDS,CUNI(1,1),NDIM,HFEIG,CWRK,LWORK*NBANDS, &
           W,  IFAIL)
      IF (IFAIL/=0) THEN
         CALL vtutor%bug("internal error in VASP from DIAGH: Call to routine ZHEEV failed! Error " &
            // "code was " // str(IFAIL), __FILE__, __LINE__)
#endif 
      ENDIF
   END SUBROUTINE DIAGH

!**********************************************************************
! helper: prints information about correlated subspace
!**********************************************************************

   SUBROUTINE PRINT_CSPACE_INFO(IU0,IU6)
      IMPLICIT NONE
      INTEGER IU0, IU6
      !local 
      INTEGER, ALLOCATABLE :: WTARGET(:)
      INTEGER              :: COUNTER 
      INTEGER              :: I
 
      ALLOCATE(WTARGET(CRPA_NUM_WAN))
      WTARGET=0
      COUNTER=0
      DO I = 1, CRPA_NUM_WAN
         IF ( NTARGET_MASK (I) == 1 ) THEN
            COUNTER = COUNTER + 1
            WTARGET(COUNTER) = I
         ENDIF
      ENDDO
      IF ( COUNTER >0 )THEN
         IF ( IU0>=0 ) WRITE(IU0,'(" Correlated subspace contains the Wannier'//&
            ' states #: ",/10I4)')WTARGET(1:COUNTER)
         IF ( IU6>=0 ) WRITE(IU6,'(" Correlated subspace contains the Wannier'//&
            ' states #: ",/10I4)')WTARGET(1:COUNTER)
      ELSE
         IF ( IU0>=0 ) WRITE(IU0,'(" No states in correalted subspace present")')
         IF ( IU6>=0 ) WRITE(IU6,'(" No states in correalted subspace present")')
      ENDIF
      DEALLOCATE(WTARGET)
     
   END SUBROUTINE PRINT_CSPACE_INFO 

!*************************SUBROUTINE INWAV_DOWNFOLD*********************
!
!   read wavefunctions from file
!
!   the routine can read WAVECAR files generated for different
!   lattice vectors and cutoff (both in the serial and parallel case)
!
!   it is also possible to restart a spin polarised calculation
!   from a non spin polarised wavefunction file or
!
!***********************************************************************

      SUBROUTINE INWAV_DOWNFOLD(IO, WDES, W, SYMM, KPOINTS, GRID, LATT_CUR, LATT_INI, ISTART, EFERMI)
      USE prec
      USE base
      USE wave_high
      USE mgrid
      USE lattice
      USE mkpoints
      USE fileio 
      USE mlr_optic , ONLY:IKPT_SELECT
 
      TYPE (in_struct)   IO
      TYPE (wavedes)  WDES
      TYPE (wavedes1) WDES1
      TYPE (wavespin) W
      TYPE (symmetry)    SYMM
      TYPE (kpoints_struct) :: KPOINTS
      TYPE (grid_3d)  GRID
      TYPE (latt)     LATT_CUR,LATT_INI
      INTEGER ISTART 
      REAL (q) :: EFERMI

! local work arrays
      REAL(q) VKPT(3)
      REAL(q) VKPT2(3)
      INTEGER,ALLOCATABLE    :: IND(:),INDI(:)
      COMPLEX(q),ALLOCATABLE :: CW1(:),CW2(:)
      LOGICAL ::  SINGLE_PREC, SPAN_RECORDS
      COMPLEX(qs),  ALLOCATABLE :: CRD(:)
      LOGICAL, POINTER :: LCONJG(:)=>null()
      INTEGER, POINTER :: INDEX(:)=>null()
      REAL(q), ALLOCATABLE :: INBUF(:)
!wavecar read in
      COMPLEX(q),ALLOCATABLE:: CELTOT(:,:)   ! global array for eigenvalues
      REAL(q)   ,ALLOCATABLE:: FERTOT(:,:)   ! global array for fermi-weights
      INTEGER    :: NK
      INTEGER    :: NKPT_COUNT
 
      INTEGER  IU0, ISP, J, K, I , IS
      INTEGER  IFAIL , IRECL 
      INTEGER  IONODE,NALLO, ISPINREAD
      INTEGER  NODE_ME, NBANDF, MALLOC
      INTEGER   NCOMM, NALLOC
      REAL(q) RTAG, RDUM, RISPIN,  RKPTSF,RBANDF 
      REAL(q) ENMAXF, RNPL
      INTEGER IREC, NKPTSF, NPLREAD, NPL
      CHARACTER(LEN=:), ALLOCATABLE :: MSG

      IU0=IO%IU0

      NALLOC = MAXVAL(WDES%NPLWKP_TOT)

      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE
#endif
!
! parse the header
!
        RTAG=0
        READ(12,REC=1,ERR=200) RDUM,RISPIN,RTAG
        IF (RTAG==53300) THEN
           IF (IU0>=0) WRITE(IU0,*) "VASP.5.3 WAVECAR encountered"
           SINGLE_PREC=.TRUE.
           SPAN_RECORDS=.TRUE.
        ELSE IF (RTAG==53310) THEN
           IF (IU0>=0) WRITE(IU0,*) "VASP.5.3 double precision WAVECAR encountered"
           SINGLE_PREC=.FALSE.
           SPAN_RECORDS=.TRUE.
        ELSE IF (RTAG==45200) THEN
           SINGLE_PREC=.TRUE.
           SPAN_RECORDS=.FALSE.
        ELSE
           IF (IU0>=0) WRITE(IU0,*) "double precision WAVECAR encountered, converting it"
           SINGLE_PREC=.FALSE.
           SPAN_RECORDS=.FALSE.
        ENDIF

        ISPINREAD=NINT(RISPIN)
        READ(12,REC=2,ERR=200) RKPTSF,RBANDF,ENMAXF, &
     &                         ((LATT_INI%A(I,J),I=1,3),J=1,3),EFERMI

        IREC=2
        NKPTSF=NINT(RKPTSF)
        NBANDF=NINT(RBANDF)

        ALLOCATE(CELTOT(NBANDF,ISPINREAD))
        CELTOT=0
        ALLOCATE(FERTOT(NBANDF,ISPINREAD))
        FERTOT=0

        !is every k-point taken into account?
        ALLOCATE(IKPT_SELECT(NKPTSF))
        DO K = 1, NKPTSF
           IKPT_SELECT(K)=K
        ENDDO
      
        CALL LATTIC(LATT_INI)
!
! if number of k-points in WAVECAR > then current number of k-points, we need to select
! appropriate k-points 
!        
       IF ( NKPTSF > WDES%NKPTS ) THEN
         IF ( IU0>=0 )WRITE(IU0,'(A,I4,A)',ADVANCE='NO')&
            ' WAVECAR contains more than',WDES%NKPTS,&
            ' k-points, downsampling...'
         DO K = 1, NKPTSF
            IKPT_SELECT(K)=-1
         ENDDO
          CALL FIND_KPOINTS_IN_WAVECAR(WDES,LATT_CUR,SYMM,KPOINTS,NKPTSF,IO%IU0)
       ENDIF 
       
!=======================================================================
! read WAVECAR file, number of bands agree
!=======================================================================
      IF (.NOT.(WDES%NRSPINORS==2 .AND. NBANDF*WDES%NRSPINORS == WDES%NB_TOT)) THEN
        spin:    DO ISP=1, MIN(WDES%ISPIN, ISPINREAD)
        NKPT_COUNT=0
        kpts: DO K=1,NKPTSF
 
        IF ( IKPT_SELECT(K)> NKPTSF ) THEN
           IREC=IRECL
        ELSE
           IRECL=IREC
        ENDIF

        IF ( IKPT_SELECT(K) > 0 ) THEN
           NKPT_COUNT = NKPT_COUNT + 1
           CALL SETWDES(WDES,WDES1,IKPT_SELECT(K))
#ifdef gammareal
           CALL GEN_INDEX_GAMMA(WDES1, LATT_CUR%B, WDES%ENMAX, INDEX, LCONJG)
#endif
        io_begin
        IFAIL=0
        IF (.NOT. SPAN_RECORDS) THEN
           IREC=IREC+1
           !read energy and occupancy for kpoint of file
           READ(12,REC=IREC,ERR=230) RNPL,VKPT, &
                (W%CELTOT(J,IKPT_SELECT(K),ISP),W%FERTOT(J,IKPT_SELECT(K),ISP),J=1,MIN(WDES%NB_TOT,NBANDF))
        ELSE
           ALLOCATE(INBUF(4+NBANDF*3))
           CALL READ_TO_BUF(12, IREC, IO%IRECLW/IO%ICMPLX*2, INBUF, IFAIL)
           RNPL   =INBUF(1)
           VKPT(1)=INBUF(2)
           VKPT(2)=INBUF(3)
           VKPT(3)=INBUF(4)
           DO J=1,MIN(WDES%NB_TOT,NBANDF)
              W%CELTOT(J,IKPT_SELECT(K),ISP)=CMPLX(INBUF(5+(J-1)*3),INBUF(6+(J-1)*3),q)
              W%FERTOT(J,IKPT_SELECT(K),ISP)=INBUF(7+(J-1)*3)
           ENDDO
           DEALLOCATE(INBUF)

        ENDIF
        NPLREAD=NINT(RNPL)
        io_end

        CALLMPI( M_bcast_i( WDES%COMM, IFAIL, 1))
        IF (IFAIL /=0) GOTO 230

        CALLMPI( M_bcast_i( WDES%COMM, IREC, 1))
        CALLMPI( M_bcast_i( WDES%COMM, NPLREAD, 1))
        NPL=WDES%NPLWKP_TOT(IKPT_SELECT(K))

        MALLOC=MAX(NPL, NPLREAD)
        ALLOCATE(CW1(MALLOC),CW2(MALLOC),CRD(MALLOC),IND(MALLOC),INDI(MALLOC))

        io_begin
        ! create an index to allow for change of cutoff or cell size
        CALL REPAD_INDEX_ARRAY(GRID, WDES%VKPT(:,IKPT_SELECT(K)), VKPT, LATT_CUR%B,  LATT_INI%B, &
                WDES%ENMAX, ENMAXF, NPL/WDES%NRSPINORS, NPLREAD/WDES%NRSPINORS, IND, INDI, MALLOC, IFAIL )
        io_end

        CALLMPI( M_bcast_i( WDES%COMM, IFAIL, 1))
        IF (IFAIL /=0) GOTO 220

        band: DO J=1,NBANDF
          IREC=IREC+1
          IF (J>WDES%NB_TOT) CYCLE
          io_begin
            IF (SINGLE_PREC) THEN
               READ(12,REC=IREC,ERR=240) (CRD(I),I=1,NPLREAD)
               CW2(1:NPLREAD)=CRD(1:NPLREAD)
            ELSE
               READ(12,REC=IREC,ERR=240) (CW2(I),I=1,NPLREAD)
            ENDIF

            CW1=0
            DO IS=1,WDES%NRSPINORS
               ! store the wave function coefficients according to new cutoff
               CALL REPAD_WITH_INDEX_ARRAY( MALLOC, IND, INDI, &
                  CW1((IS-1)*NPL/WDES%NRSPINORS+1), CW2((IS-1)*NPLREAD/WDES%NRSPINORS+1))
            ENDDO

          io_end
#ifdef MPI
          IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
             CALLMPI( M_bcast_z(WDES%COMM_KINTER,CW1,SIZE(CW1)) )
          END IF
#endif
#ifdef gammareal
          ! gamma only go from seriel layout to parallel layout if required
          CALL DIS_PW_BAND_GAMMA(WDES1, J, CW1, W%CW(1,1,IKPT_SELECT(K),ISP), INDEX(1), LCONJG(1))
#else
          CALL DIS_PW_BAND(WDES1, J, CW1, W%CW(1,1,IKPT_SELECT(K),ISP))
#endif
        ENDDO band

#ifdef gammareal
        CALL FREE_INDEX_GAMMA(INDEX, LCONJG)
#endif
        DEALLOCATE(CW1,CW2,CRD,IND,INDI)
        ELSE    !in this case we don't store but proceed to read the WAVECAR
           io_begin
              IF(IO%NWRITE>2) WRITE(IU0,'(A,I5,A)')' kpoint,',K,' skipped.'
              IFAIL=0
              IF (.NOT. SPAN_RECORDS) THEN
                 IREC=IREC+1
                 !read energy and occupancy for kpoint of file
                 READ(12,REC=IREC,ERR=230) RNPL,VKPT, &
                      (CELTOT(J,ISP),FERTOT(J,ISP),J=1,MIN(WDES%NB_TOT,NBANDF))
              ELSE
                 ALLOCATE(INBUF(4+NBANDF*3))
                 CALL READ_TO_BUF(12, IREC, IO%IRECLW/IO%ICMPLX*2, INBUF, IFAIL)
                 RNPL   =INBUF(1)
                 VKPT(1)=INBUF(2)
                 VKPT(2)=INBUF(3)
                 VKPT(3)=INBUF(4)
                 DO J=1,MIN(WDES%NB_TOT,NBANDF)
                    CELTOT(J,ISP)=CMPLX(INBUF(5+(J-1)*3),INBUF(6+(J-1)*3),q)
                    FERTOT(J,ISP)=INBUF(7+(J-1)*3)
                 ENDDO
                 DEALLOCATE(INBUF)
   
              ENDIF
              NPLREAD=NINT(RNPL)
           io_end
   
           CALLMPI( M_bcast_i( WDES%COMM, IFAIL, 1))
           IF (IFAIL /=0) GOTO 230
   
           CALLMPI( M_bcast_i( WDES%COMM, IREC, 1))
           CALLMPI( M_bcast_i( WDES%COMM, NPLREAD, 1))
           NPL=1
   
           MALLOC=MAX(NPL, NPLREAD)
           ALLOCATE(CW1(MALLOC),CW2(MALLOC),CRD(MALLOC),IND(MALLOC),INDI(MALLOC))
   
           band2: DO J=1,NBANDF
             IREC=IREC+1
             IF (J>WDES%NB_TOT) CYCLE
             io_begin
               IF (SINGLE_PREC) THEN
                  READ(12,REC=IREC,ERR=240) (CRD(I),I=1,NPLREAD)
                  CW2(1:NPLREAD)=CRD(1:NPLREAD)
               ELSE
                  READ(12,REC=IREC,ERR=240) (CW2(I),I=1,NPLREAD)
               ENDIF
             io_end
           ENDDO band2
   
#ifdef gammareal
           CALL FREE_INDEX_GAMMA(INDEX, LCONJG)
#endif
           DEALLOCATE(CW1,CW2,CRD,IND,INDI)
   
           ENDIF   !decide if something is stored
        ENDDO kpts

        IF (NKPTSF > WDES%NKPTS .AND. ISP == MIN(WDES%ISPIN, ISPINREAD) ) THEN
           IF ( NKPT_COUNT == WDES%NKPTS ) THEN
              IF(IU0>=0)WRITE(IU0,*)'Success'
           ELSE
              IF(IU0>=0)WRITE(IU0,*)' FAILED!'
           ENDIF 
           IREC=IREC+(NKPTSF-WDES%NKPTS)*(WDES%NB_TOT+1)
        ENDIF

        ! copy eigenvalues and weights to all nodes
        NCOMM=WDES%NB_TOT*WDES%NKPTS
        CALLMPI( M_bcast_z(WDES%COMM, W%CELTOT(1,1,ISP),NCOMM ))
        CALLMPI( M_bcast_d(WDES%COMM, W%FERTOT(1,1,ISP),NCOMM ))

        ENDDO spin

        IF (ISPINREAD > WDES%ISPIN .AND. IU0>=0 ) THEN
           WRITE(IU0,*) 'down-spin wavefunctions not read'
        ENDIF

        IF (NBANDF<WDES%NB_TOT) THEN
           IF (IU0>=0) WRITE(IU0,*) 'random initialization beyond band ',NBANDF
           CALL WFINIT(WDES, W, 1E10_q, NBANDF+1) ! ENINI=1E10 not cutoff restriction
        ENDIF

        IF (IU0 >= 0) WRITE(IU0,*) 'the WAVECAR file was read successfully'

        IF (WDES%ISPIN<=ISPINREAD) THEN
           IF ( NKPTSF==WDES%NKPTS) DEALLOCATE(IKPT_SELECT)
           RETURN
        ENDIF
!
!  spin down is missing
!
        IF (IU0>=0) &
        WRITE(IU0,*) 'No down-spin wavefunctions found', &
     &             ' --> setting down-spin equal up-spin ...'

        DO K=1,WDES%NKPTS
          W%CELTOT(1:WDES%NB_TOT,K,2)=W%CELTOT(1:WDES%NB_TOT,K,1)
          W%FERTOT(1:WDES%NB_TOT,K,2)=W%FERTOT(1:WDES%NB_TOT,K,1)
          NPL=WDES%NPLWKP(K)
          W%CW(1:NPL,1:WDES%NBANDS,K,2)=W%CW(1:NPL,1:WDES%NBANDS,K,1)
        ENDDO

        IF ( NKPTSF==WDES%NKPTS) DEALLOCATE(IKPT_SELECT)
        RETURN
!=======================================================================
! read collinear WAVECAR file for a non collinear run
!=======================================================================
      ELSE
        CALL vtutor%error("Sorry: downsampling for non-collinear WAVECARs currently not supported")
      ENDIF

      MSG = "ERROR: while reading WAVECAR, file is incompatible"
      IF (WDES%ENMAX /= ENMAXF) &
         MSG = MSG // "the energy cutoff has changed (new,old) " // str(WDES%ENMAX) // " " // str(ENMAXF)
      IF  (NBANDF /= WDES%NB_TOT) &
         MSG = MSG // "the number of bands has changed (new,old) " // str(WDES%NB_TOT) // " " // str(NBANDF)
      IF (NKPTSF /= WDES%NKPTS) &
         MSG = MSG // "the number of k-points has changed (new,old) " // str(WDES%NKPTS) // " " // str(NKPTSF)
      CALL vtutor%error(MSG)
!=======================================================================
! can not do anything with WAVECAR 
! hard stop, pull all breaks
!=======================================================================
  200 CONTINUE

      CALL vtutor%error("ERROR: while reading WAVECAR, header is corrupt")

  220 CONTINUE

      CALL vtutor%error("ERROR: while reading WAVECAR, plane wave coefficients changed " // str(NPL) &
         // " " // str(NPLREAD))

  230 CONTINUE

      CALL vtutor%error("ERROR: while reading eigenvalues from WAVECAR " // str(K) // " " // str(ISP))

  240 CONTINUE

      CALL vtutor%error("ERROR: while reading plane wave coef. from WAVECAR " // str(K) // " " // &
         str(ISP) // " " // str(J))
      
      END SUBROUTINE INWAV_DOWNFOLD


!**********************************************************************
! finds correct k-points in WAVECAR
!**********************************************************************

   SUBROUTINE FIND_KPOINTS_IN_WAVECAR(WDES,LATT_CUR,SYMM,KPOINTS,&
      IRZ_WAVECAR,IU0)
      USE base
      USE wave_high
      USE mkpoints
      USE constant
      USE full_kpoints
      USE lattice
      USE mlrf_main ! contains KINTER
      USE mlr_optic , ONLY: IKPT_SELECT

      IMPLICIT NONE           
      !kpoints of WAVECAR
      TYPE (kpoints_struct) :: KPOINTS
      TYPE (kpoints_struct) :: KPOINTS_WAVECAR
      TYPE (wavedes)        :: WDES
      TYPE (latt)           :: LATT_CUR
      TYPE (symmetry)       :: SYMM
      INTEGER               :: IRZ_WAVECAR
      INTEGER               :: IU0
      !local 
      INTEGER, PARAMETER    :: KMULT_MAX=36
      INTEGER               :: KMULT
      INTEGER               :: NK,NKP
      INTEGER, DIMENSION(3) :: MULTIPL
      
      !find correct k-point grid
      !create a denser k-point grid 
      KPOINTS_WAVECAR=KPOINTS
      nkpts_grow:DO KMULT = 1, KMULT_MAX
         MULTIPL(:) = KMULT
         CALL RD_KPOINTS(KPOINTS_WAVECAR, LATT_CUR, &
            SYMM%ISYM>=0, SYMM%ISYM<0, -1, -1, MULTIPL = MULTIPL)
         IF ( KPOINTS_WAVECAR%NKPTS == IRZ_WAVECAR) EXIT nkpts_grow
      ENDDO nkpts_grow

      !if k-point grid found  KMULT is below KMULT_MAX
      IF  ( KMULT > KMULT_MAX ) THEN
         CALL vtutor%error("Error FIND_KPOINTS_IN_WAVECAR: \n cannot find k-point grid corresponding &
            &to WAVECAR")
      ENDIF 

      !now we have the correct grid, 
      !find out which kpoints of original grid corresponds to current grid
      KMULT=0
      DO NK = 1, KPOINTS_WAVECAR%NKPTS
         !check if current k-point of list corresponds to k-point of file
         DO NKP=1, WDES%NKPTS
            IF (LIDENTICAL_KPOINT(WDES%VKPT(:,NKP),KPOINTS_WAVECAR%VKPT(:,NK))) THEN
               KMULT=KMULT+1
               IKPT_SELECT(NK)=NKP
            ENDIF
         ENDDO
#ifdef debug
      IF ( IU0>=0 )WRITE(*,*)NK,IKPT_SELECT(NK)
#endif 
      ENDDO

      IF ( KMULT/=KPOINTS%NKPTS) THEN
         CALL vtutor%error("Error FIND_KPOINTS_IN_WAVECAR: not correct points found " // str(KMULT) &
            // " " // str(KPOINTS%NKPTS))
      ENDIF     
 
   END SUBROUTINE FIND_KPOINTS_IN_WAVECAR


!**********************************************************************
!
! small dumping routine for P, P^2, P^3, ....
!
!**********************************************************************
  SUBROUTINE DUMP_PROJ( P, NB_TOT, NORDER, NDUMP, NSTOP, LPRINT)
    GDEF              :: P (NB_TOT,NB_TOT)
    INTEGER           :: NB_TOT, NORDER
    INTEGER,OPTIONAL  :: NDUMP
    INTEGER,OPTIONAL  :: NSTOP
    LOGICAL,OPTIONAL  :: LPRINT
    !local 
    INTEGER           :: N1, N2, NPL2
    INTEGER           :: I 
    GDEF              :: AUX(NB_TOT,NB_TOT)
    GDEF              :: BUX(NB_TOT,NB_TOT)
    GDEF              :: CUX(NB_TOT,NB_TOT)

    AUX=P
    CUX=P
    BUX=P
 
DO   I = 1, MIN(NORDER,4)  

    WRITE(*,'(A,I2)')' P^',I

    IF ( PRESENT( NDUMP ) .AND. PRESENT(NSTOP) ) THEN
       DO N1=NDUMP,NSTOP
          WRITE(*,1)N1,(REAL( AUX(N1,N2) ,KIND=q) ,N2=NDUMP,NSTOP)
       ENDDO
       WRITE(*,*)
#ifndef gammareal
       IF ( PRESENT (LPRINT ) ) THEN
       DO N1=NDUMP,NSTOP
          WRITE(6,2)N1,(AIMAG( AUX(N1,N2) ) ,N2=NDUMP,NSTOP)
       ENDDO
       WRITE(*,*)
       ENDIF
#endif
    ELSE IF ( PRESENT( NDUMP ) ) THEN
       NPL2=MIN(NDUMP,NB_TOT)
       DO N1=1,NPL2
          WRITE(*,1)N1,(REAL( AUX(N1,N2) ,KIND=q) ,N2=1,NPL2)
       ENDDO
       WRITE(*,*)

#ifndef gammareal
       IF ( PRESENT (LPRINT ) ) THEN
       DO N1=1,NPL2
          WRITE(6,2)N1,(AIMAG( AUX(N1,N2)),N2=1,NPL2)
       ENDDO
       WRITE(*,*)
       ENDIF 
#endif
    ELSE
       NPL2=MIN(12,NB_TOT)
       DO N1=1,NPL2
          WRITE(*,1)N1,(REAL( AUX(N1,N2) ,KIND=q) ,N2=1,NPL2)
       ENDDO
       WRITE(*,*)
#ifndef gammareal
       IF ( PRESENT (LPRINT ) ) THEN
       DO N1=1,NPL2
          WRITE(6,2)N1,(AIMAG( AUX(N1,N2)),N2=1,NPL2)
       ENDDO
       WRITE(*,*)
       ENDIF 
#endif
    ENDIF 

       CUX=TRANSPOSE(GCONJG(AUX))
       BUX=MATMUL(CUX,AUX)
       AUX=BUX
ENDDO 

1   FORMAT(1I3,3X,40F9.5)
2   FORMAT(1I3,3X,40F9.5)

  END SUBROUTINE DUMP_PROJ


!the following routines are for joint diagonalization 

!**********************************************************************
!
! diagonalizes a set of matrices as good as possible 
! if LSINGLE is true only matrix at k=NKTOT is diagonalized with 
! Jacobi algorithm. The matlab code is written by 
!  author       = "Jean-Fran\c{c}ois Cardoso and Antoine Souloumiac",
!  journal      = "{SIAM} J. Mat. Anal. Appl.",
!  title        = "Jacobi angles for simultaneous diagonalization",
!  pages        = "161--164",
!  volume       = "17",
!  number       = "1",
!  month        = jan,
!  year         = {1996}
!
!**********************************************************************

   SUBROUTINE JOINT_DIAGONALIZER(D,MAT,THRESH,M,NKPTS,ISPIN,IU0,IU6,MAXITER,LSINGLE)
      USE prec
      USE constant
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      GDEF :: D(:,:,:)     !matrix which diagonalizes all input matrices 
      GDEF :: MAT(:,:,:,:) !set of matrices
      REAL(q)    :: THRESH       !convergence criterium 
      INTEGER    :: M            !dimension of one matrix
      INTEGER    :: NKPTS        !number of matrices to be diagonalized
      INTEGER    :: ISPIN        !spin number
      INTEGER    :: IU0, IU6     !for writing
      INTEGER    :: MAXITER      !number of allowed iterations
      LOGICAL,OPTIONAL :: LSINGLE
      !local 
      INTEGER                :: I, J, K, ISP, ITER
      INTEGER                :: ALGOINFO 
      INTEGER                :: NM
      COMPLEX(q),ALLOCATABLE :: A(:,:,:)     !working array 
      COMPLEX(q),ALLOCATABLE :: AT1(:,:)     !working array 
      COMPLEX(q)             :: B(3,3)       !rotation matrix
      COMPLEX(q)             :: BT(3,3)      !the conjugate transpose of B
      INTEGER,ALLOCATABLE    :: IP(:)
      INTEGER,ALLOCATABLE    :: IQ(:)
      COMPLEX(q),ALLOCATABLE :: G(:,:)       !Givens rotation matrix
      COMPLEX(q),ALLOCATABLE :: GU(:,:)      !Givens rotation matrix
      COMPLEX(q)             :: VCP(3,3)
      COMPLEX(q)             :: DIA(3,3)
      INTEGER                :: KU(3)        !sort array 
      REAL(q)                :: KMAT(3,3)
      REAL(q)                :: ANGLES(3)
      INTEGER                :: PAIR(2)
      REAL(q)                :: C
      COMPLEX(q)             :: S
      COMPLEX(q),ALLOCATABLE :: V(:,:)       !diagonalizer matrix
      COMPLEX(q),ALLOCATABLE :: V1(:,:)      !rotated matrix
      COMPLEX(q),ALLOCATABLE :: V2(:,:)      !rotated matrix
      COMPLEX(q),ALLOCATABLE :: A1(:,:)      !rotated matrix
      COMPLEX(q),ALLOCATABLE :: A2(:,:)      !rotated matrix
      REAL(q)                :: LA(3)
      INTEGER                :: P,QI
      REAL(q)                :: RAUX(3,3)
      REAL(q)                :: SBAK
      LOGICAL                :: LCONTINUE 
      INTEGER                :: NKTOT
      
      IF ( PRESENT(LSINGLE) .AND. LSINGLE ) THEN
         NKTOT=1
      ELSE
         NKTOT=NKPTS
      ENDIF 

      
      !some standard test
      IF ( SIZE(MAT,1) /= SIZE(MAT,2) ) THEN
         CALL vtutor%bug("internal error in VASP from JOINT_DIAGONALIZER: input matrices are not " &
            // "square! " // str(SIZE(MAT,1)) // " " // str(SIZE(MAT,2)), __FILE__, __LINE__)
      ENDIF 
      IF ( SIZE(D,1) /= SIZE(MAT,1) ) THEN
         CALL vtutor%bug("internal error in VASP from JOINT_DIAGONALIZER: sizes of matrices " &
            // "inconsistent! " // str(SIZE(D,1)) // " " // str(SIZE(MAT,1)), __FILE__, __LINE__)
      ENDIF 
      IF ( SIZE(D,1) /= SIZE(D,2) ) THEN
         CALL vtutor%bug("internal error in VASP from JOINT_DIAGONALIZER: diagonalizing matrix must " &
            // "be square! " // str(SIZE(D,1)) // " " // str(SIZE(D,2)), __FILE__, __LINE__)
      ENDIF 
      IF ( SIZE(D,3) < ISPIN ) THEN
         CALL vtutor%bug("internal error in VASP from JOINT_DIAGONALIZER: diagonalizing matrix must " &
            // "contain enough spin dimensions! " // str(SIZE(D,3)) // " " // str(ISPIN), __FILE__, __LINE__)
      ENDIF 
      IF ( SIZE(MAT,2) < M ) THEN
         CALL vtutor%bug("internal error in VASP from JOINT_DIAGONALIZER: dimension of input " &
            // "matrices too small! " // str(SIZE(MAT,2)) // " " // str(M), __FILE__, __LINE__)
      ENDIF 
      IF ( SIZE(MAT,3) < NKTOT ) THEN
         CALL vtutor%bug("internal error in VASP from JOINT_DIAGONALIZER: too few input matrices " &
            // "given! " // str(SIZE(MAT,3)) // " " // str(NKTOT), __FILE__, __LINE__)
      ENDIF 
      IF ( SIZE(MAT,4) < ISPIN ) THEN
         CALL vtutor%bug("internal error in VASP from JOINT_DIAGONALIZER: too few input matrices " &
            // "given (ISPIN)! " // str(SIZE(MAT,3)) // " " // str(ISPIN), __FILE__, __LINE__)
      ENDIF 
    
      !A will be the working array 
      NM=NKTOT*M
      ALLOCATE(A(M,M,NKTOT))
      ALLOCATE(AT1(M,2*NKTOT))
      !initialize loop arrays
      ALLOCATE(IP(NKTOT),IQ(NKTOT))
      !Givens matrices
      ALLOCATE(G(3,NKTOT),GU(2,2))
      ALLOCATE(V(M,M)) 
      ALLOCATE(V1(M,2),V2(M,2)) 
      ALLOCATE(A1(2,NM),A2(2,NM)) 

      !----------------------------------------------------------------------
      spin : DO ISP=1, ISPIN 
      !----------------------------------------------------------------------
      
         A=(0._q,0._q)
         !resave input matrices 
         IF ( PRESENT(LSINGLE) .AND. LSINGLE ) THEN
            DO I=1, M
               DO J=1, M
#ifdef gammareal
                  A(I,J,NKTOT)=CMPLX(MAT(I,J,NKPTS,ISP),0._q,q)
#else
                  A(I,J,NKTOT)=MAT(I,J,NKPTS,ISP)
#endif 
               ENDDO
            ENDDO
         ELSE
            DO I=1, M
               DO J=1, M
                  DO K=1, NKTOT
#ifdef gammareal
                     A(I,J,K)=CMPLX(MAT(I,J,K,ISP),0._q,q)
#else
                     A(I,J,K)=MAT(I,J,K,ISP)
#endif
                  ENDDO
               ENDDO
            ENDDO
         ENDIF 
         !set B
         B(1,1)=(1._q,0._q)
         B(1,2)=(0._q,0._q)
         B(1,3)=(0._q,0._q)
         B(2,1)=(0._q,0._q)
         B(2,2)=(1._q,0._q)
         B(2,3)=(1._q,0._q)
         B(3,1)=(0._q,0._q)
         B(3,2)=(0._q,-1._q)
         B(3,3)=(0._q,1._q)
         BT=CONJG(TRANSPOSE(B))

         IP=0
         IQ=0
   
         G=(0._q,0._q)
         GU=(0._q,0._q)

         !initialize VCP and DIA
         VCP=(0._q,0._q)
         DIA=VCP
         LA=0
         KMAT=0._q
         ANGLES=0._q
         PAIR=0
         C=0._q
         S=(0._q,0._q)
          
         !initialize V to identity 
         V=(0._q,0._q)
         DO I=1,M
            V(I,I)=(1._q,0._q)
         ENDDO
      
         !----------------------------------------------------------------------  
         minloop: DO ITER = 1, MAXITER
         !----------------------------------------------------------------------  
            LCONTINUE=.FALSE.  !sould we continue 
            !p is the row index, q is the column index
            DO P=1,M-1
               DO QI=P+1,M

                  !we consider 3x3 submatrices 
                  DO I = 1, NKTOT
                     IP(I)=P+(I-1)*M
                     IQ(I)=QI+(I-1)*M
                  ENDDO
                  !compute the Givens angles
                  DO I=1,NKTOT
                     G(1,I)=A(P,IP(I),1)-A(QI,IQ(I),1) 
                     G(2,I)=A(P,IQ(I),1) 
                     G(3,I)=A(QI,IP(I),1) 
                  ENDDO
                  RAUX=REAL(MATMUL(B,MATMUL(MATMUL(G,CONJG(TRANSPOSE(G))),BT)),q)

                  !get the eigensystem
                  CALL RSYS_EIGENSYSTEM(RAUX,VCP,LA,IU0)
                  !maybe sort eigenvalues (LAPACK should do this anyhow)
                  KU(1)=1; KU(2)=2; KU(3)=3 
 
                  !define angles
                  ANGLES(1:3)=REAL(VCP(1:3,KU(3)),q)
                  
                  !maybe the first angle is negative, then multiply by (-1)
                  IF ( ANGLES(1) < 0 ) ANGLES = - ANGLES
                  
                  C= SQRT(0.5_q+ANGLES(1)/2._q)
                  S=0.5_q*( ANGLES(2)- (0._q,1._q)*ANGLES(3))/C              
                  !------------------------------------------------------------
                  !convergence test 
                  IF ( ABS(S) > THRESH) THEN
                     LCONTINUE = .TRUE.  !continue iterations in this case 
                     !if test not passed -> update A and V by a Givens rotation
                     PAIR(1) = P
                     PAIR(2) = QI
                     !Givens rotation matrix 
                     GU(1,1) = CMPLX(C,0._q,q)
                     GU(1,2) = -CONJG(S) 
                     GU(2,1) = S 
                     GU(2,2) = CMPLX(C,0._q,q)
 
                     !rotate relevant part of V 
                     DO I =1 ,2
                        V1(1:M,I)=V(1:M,PAIR(I))
                     ENDDO
                     V2 = MATMUL(V1,GU)
                     DO I =1 ,2
                        V(1:M,PAIR(I))=V2(1:M,I)
                     ENDDO

                     !rotate A as well
                     DO I =1 ,2
                        A1(I,1:NM)=A(PAIR(I),1:NM,1)
                     ENDDO
                     A2 = MATMUL(CONJG(TRANSPOSE(GU)),A1)
                     DO I =1 ,2
                        A(PAIR(I),1:NM,1)=A2(I,1:NM)
                     ENDDO
                    
                     !save relevant part of A temporarily to ATMP
                     AT1=(0._q,0._q)
                     DO I = 1, M
                        DO J = 1, 2*NKTOT
                           IF ( J <= NKTOT ) THEN
                              AT1(I,J)=C*A(I,IP(J),1)+S*A(I,IQ(J),1)
                           ELSE
                              AT1(I,J)=-CONJG(S)*A(I,IP(J-NKTOT),1)+C*A(I,IQ(J-NKTOT),1)
                           ENDIF 
                        ENDDO
                     ENDDO
                     !overwrite these matrix elements of A by AT1
                     DO I = 1, M
                        DO J = 1, 2*NKTOT
                           IF ( J <= NKTOT ) THEN
                              A(I,IP(J),1) = AT1(I,J)  
                           ELSE
                              A(I,IQ(J-NKTOT),1) = AT1(I,J)  
                           ENDIF 
                        ENDDO
                     ENDDO
                  ENDIF   !update necessary?
                  !----------------------------------------------------------
               ENDDO  !columns
            ENDDO  !rows

            !alternative breaking contition for single diagonalizations
            IF ( PRESENT(LSINGLE) .AND. LSINGLE ) THEN
               IF ( ITER == 1 ) THEN
                  SBAK=ABS(S) 
               ELSE
                  SBAK=ABS(ABS(S)-SBAK) 
               ENDIF 
               IF ( SBAK < THRESH ) THEN
                  LCONTINUE=.FALSE.              
               ELSE
                  SBAK=ABS(S)
               ENDIF
            ENDIF

            !breaking condition    
            IF ( ITER < MAXITER ) THEN
               IF ( PRESENT(LSINGLE) .AND. LSINGLE ) THEN
               ELSE
                  !dump convergence info 
                  IF ( IU0 >= 0 ) THEN
                     WRITE(IU0,'(A,I4,A,E13.6)')' iteration: ', ITER, ' error ',ABS(S)
                  ENDIF 
               ENDIF 
               IF ( .NOT. LCONTINUE ) THEN
                  ALGOINFO=1
                  EXIT minloop
               ENDIF 
            ELSE
               IF ( IU0>=0 ) WRITE(IU0,'(A,I5,A,E13.4)') 'WARNING from JOINT_DIAGONALIZER:'//&
                  'not converged after ',MAXITER,' iterations. ERROR=',ABS(S)
                  ALGOINFO=-1
            ENDIF 
            

         !-------------------------------------------------------------------
         ENDDO minloop  
         !-------------------------------------------------------------------

         !save the result
         !resave input matrices 
         IF ( PRESENT(LSINGLE) .AND. LSINGLE ) THEN
            DO I=1, M
               DO J=1, M
#ifdef gammareal
                  MAT(I,J,NKPTS,ISP)=REAL(A(I,J,NKTOT),q)
#else
                  MAT(I,J,NKPTS,ISP)=A(I,J,NKTOT)
#endif
               ENDDO
            ENDDO
            !dump convergence info 
            !IF ( IU0 >= 0 ) THEN
            !   WRITE(IU0,'(A,I4,A,E13.6)')'K-Point ',NKPTS,' done E=',ABS(S)
            !ENDIF 
         ELSE
            DO I=1, M
               DO J=1, M
                  DO K=1, NKTOT
#ifdef gammareal
                     MAT(I,J,K,ISP)=A(I,J,K)
#else
                     MAT(I,J,K,ISP)=REAL(A(I,J,K),q)
#endif
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         !V is the diagonalizer 
         DO I=1, M
            DO J=1, M
#ifdef gammareal
               D(I,J,ISP)=REAL(V(I,J),q)
#else
               D(I,J,ISP)= V(I,J)
#endif
            ENDDO
         ENDDO

#ifdef debug
WRITE(*,*)'Diagonalizer for spin ',ISP
DO I = 1, M
   WRITE(*,'(12F14.4)') REAL(D(I,1:M,ISP),q)
ENDDO
WRITE(*,*)' Imaginary part'
DO I = 1, M
   WRITE(*,'(12F14.4)') REAL((0._q,-1._q)*D(I,1:M,ISP),q)
ENDDO
K=3
WRITE(*,*)'NK=',K
DO I = 1, M
   WRITE(*,'(12F20.16)') REAL(MAT(I,1:M,K,ISP),q)
ENDDO
WRITE(*,*)' Imaginary part'
DO I = 1, M
   WRITE(*,'(12F20.16)') REAL((0._q,-1._q)*MAT(I,1:M,K,ISP),q)
ENDDO
#endif 

      !----------------------------------------------------------------------
      ENDDO spin 
      !----------------------------------------------------------------------
      RETURN
   END SUBROUTINE JOINT_DIAGONALIZER


!**********************************************************************
!
! gets the eigensystem of real symmetric matrix 
!
!**********************************************************************
    
   SUBROUTINE RSYS_EIGENSYSTEM(C,EV,E,IU0)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      REAL(q)                :: C(:,:)
      COMPLEX(q)             :: EV(:,:)
      REAL(q)                :: E(:)
      INTEGER                :: IU0
      !local
      INTEGER                :: M, N
      INTEGER                :: LDC
      INTEGER                :: I, J
      INTEGER,PARAMETER      :: LWMAX=10000
      INTEGER                :: INFO, LWORK
      REAL(q),ALLOCATABLE    :: RWORK(:) 
      REAL(q),ALLOCATABLE    :: CC(:,:)            
      REAL(q),ALLOCATABLE    :: W(:)            
      REAL(q),ALLOCATABLE    :: WORK(:)            

      !get the size of the matrix and the LDs
      M=SIZE(C,2)
      N=SIZE(C,1)
      LDC = M
    
      IF ( SIZE(EV,1) /= N ) THEN
         CALL vtutor%bug("internal error in VASP from RSYS_EIGENSYSTEM: EV not proper size " // &
            str(SIZE(EV,1)) // " " // str(N), __FILE__, __LINE__)
      ENDIF 
      IF ( SIZE(E) /= N ) THEN
         CALL vtutor%bug("internal error in VASP from RSYS_EIGENSYSTEM: E not proper size " // &
            str(SIZE(E,1)) // " " // str(N), __FILE__, __LINE__)
      ENDIF 
      IF ( M /= N ) THEN
         CALL vtutor%bug("internal error in VASP from RSYS_EIGENSYSTEM: matrix not square " // &
            str(M) // " " // str(N), __FILE__, __LINE__)
      ENDIF 

      !we work in CC
      ALLOCATE(CC(N,M))
      CC=C

      ALLOCATE(RWORK(2*LDC))
      ALLOCATE(W(LDC), WORK(LWMAX))

      ! Query the optimal workspace.
      LWORK = -1
      CALL DSYEV( 'Vectors', 'Upper', N, CC, LDC, W, WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

      ! solve the problem 
      CALL DSYEV( 'Vectors', 'Upper', N, CC, LDC, W, WORK, LWORK, INFO )

      ! Check for convergence.
      IF( INFO /= 0 ) THEN
         CALL vtutor%bug("internal error in VASP from RSYS_EIGENSYSTEM: DSYEV returns " // str(INFO), __FILE__, __LINE__)
      ENDIF
       
      !save eigenvectors
      DO I = 1, M
         DO J = 1, M
            EV(I,J)=CMPLX(CC(I,J),0._q,q)
         ENDDO
         E(I)=W(I)
      ENDDO
       
      !get rid of auxilary arrays
      DEALLOCATE(CC)
      DEALLOCATE(RWORK,W,WORK)
      RETURN
   END SUBROUTINE  RSYS_EIGENSYSTEM

!**********************************************************************
!
! tags joint-density of states with vanishing energy difference. 
! these states appear only in metals and are neglected in chi_base.F
!
!**********************************************************************
   SUBROUTINE TAG_SELF_CORRELATION_STATES(WDES,W , IO)
      USE base
      USE wave_high
      USE fileio
      USE subrot_cluster
      TYPE (wavedes)    :: WDES
      TYPE (wavespin)   :: W
      TYPE (in_struct)  :: IO
      !local 
      REAL(q), PARAMETER  :: EMPTY_THRESHHOLD=0.00001 !the same as in chi_base.F
      INTEGER             :: N1, N2, NK1, NK2, ISP
   
      ALLOCATE(LSELFCORRELATION_STATE(WDES%NB_TOT,WDES%NKPTS,WDES%ISPIN))
      LSELFCORRELATION_STATE=.FALSE.   
     
      DO ISP = 1, WDES%ISPIN
         DO NK1 = 1, WDES%NKPTS
         DO NK2 = 1, WDES%NKPTS
            DO N1 = 1, WDES%NB_TOT    !occupied
            IF ( ABS(W%FERTOT(N1,NK1,ISP)) < EMPTY_THRESHHOLD ) CYCLE !N1,NK1 -> occupied
            DO N2 = WDES%NB_TOT,1,-1  !unoccupied
            IF ( ABS(W%FERTOT(N2,NK2,ISP)-1) < EMPTY_THRESHHOLD ) CYCLE !N2,NK2 -> unoccupied

               IF ( ABS(W%CELTOT(N2, NK2, ISP)-W%CELTOT(N1,NK1,ISP)) < DEGENERACY_THRESHOLD )THEN 
                  !only tag states with indices n1, nk1, and isp
                  LSELFCORRELATION_STATE(N1,NK1,ISP)=.TRUE.
               ENDIF 
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDDO

#ifdef debug
IF ( WDES%COMM%NODE_ME == 1 ) THEN
      WRITE(*,*)' These states are self-correlation states '
      DO ISP = 1, WDES%ISPIN
         DO NK1 = 1, WDES%NKPTS
            DO N1 = 1, WDES%NB_TOT
              IF ( LSELFCORRELATION_STATE(N1,NK1,ISP))&
               WRITE(*,'(3I4,3F12.6)')N1,NK1,ISP,&
                 REAL(W%CELTOT(N1,NK1,ISP),q),REAL(W%FERTOT(N1,NK1,ISP),q),DEGENERACY_THRESHOLD
            ENDDO
         ENDDO
      ENDDO
ENDIF
CALL vtutor%stopCode()
#endif
   ENDSUBROUTINE TAG_SELF_CORRELATION_STATES

!******************** ROTATE_ACTIVE_OCCUPIED_ORBITALS ******************
!! rotates occupied orbitals such that states defined in IDXACTIVE are 
!! highest occupied states
!***********************************************************************
      SUBROUTINE ROTATE_ACTIVE_OCCUPIED_ORBITALS( WDES, W, IDXACTIVE, &
         SYMM,LATT_CUR,LATT_INI,T_INFO,KPOINTS,GRID,NONL_S, P, DYN, IO)
      USE dfast
      USE wave_high
      USE mkpoints
      USE full_kpoints
      USE kpoints_change
      USE base
      USE lattice
      USE poscar
      USE nonl_struct_def 
      USE pseudo_struct_def
      USE msymmetry
      USE chi_glb, ONLY: XI_EMPTY_THRESHHOLD
      USE pead, ONLY: PEAD_RESETUP_WDES,LPEAD_SYM_RED
      USE mlwf, ONLY: READ_WANPROJ_FILE, LWRITE_WANPROJ
      USE locproj, ONLY: WANPROJFILE
      IMPLICIT NONE
      TYPE(wavedes)         :: WDES
      TYPE(wavespin)        :: W
      INTEGER,INTENT(IN)    :: IDXACTIVE(:)
      TYPE (symmetry)       :: SYMM
      TYPE (latt)           :: LATT_CUR
      TYPE (latt)           :: LATT_INI
      TYPE (type_info)      :: T_INFO
      TYPE (kpoints_struct) :: KPOINTS
      TYPE (grid_3d)        :: GRID
      TYPE (nonl_struct)    :: NONL_S
      TYPE (potcar)         :: P(T_INFO%NTYP)
      TYPE (dynamics)       :: DYN
      TYPE(in_struct)       :: IO
      ! local variables
      TYPE(wavedes1) WDES1 
      TYPE(wavefuna) WA
      INTEGER I,J,K
      INTEGER ISP,NK
      INTEGER IHOMO
      GDEF, ALLOCATABLE :: U(:,:)
      GDEF, ALLOCATABLE :: UEFF(:,:)
      GDEF, ALLOCATABLE :: S(:,:)
      INTEGER, ALLOCATABLE :: ISORT(:)
      INTEGER, ALLOCATABLE :: NSORT(:)

      PROFILING_START('rotate_active_occupied_orbitals')

      ! nothing to do if IDXACTIVE is not set 
      IF( SUM( IDXACTIVE(:) )<1 ) THEN
         PROFILING_STOP('rotate_active_occupied_orbitals')
         RETURN
      ENDIF

      IF ( LUSEBLOCH ) THEN
         CRPA_NUM_WAN = WDES%NB_TOT
      ELSE 
         IF ( LUSEMLWF ) THEN
            CRPA_NUM_WAN = MLWF_GLOBAL%num_wann
         ELSE
            ! read wanproj file and set CRPA_NUM_WAN
            CALL READ_WANPROJ_FILE(WANPROJFILE,WDES,IO,LUFROMFILE,CRPA_NUM_WAN,UWAN)
         ENDIF
      ENDIF

      ! sanity check
      DO I = 1, SIZE( IDXACTIVE ) 
         IF ( IDXACTIVE( I ) > CRPA_NUM_WAN .OR. IDXACTIVE( I ) < 0 )  THEN
            CALL vtutor%error("IDXACTIVE seems corrupt " // str(IDXACTIVE(I)) // " " // &
               str(CRPA_NUM_WAN) // " " // str(WDES%NB_TOT))
         ENDIF
      ENDDO

      !switch off symmetry if UWAN files are written
      IF ((SYMM%ISYM>=0.AND. .NOT.WDES%LGAMMA) .AND. LWRITE_WANPROJ) THEN
         ! switch of symmetry
         CALL NOSYMM(LATT_CUR%A,T_INFO%NTYP,T_INFO%NITYP,T_INFO%NIONS,SYMM%PTRANS, &
              SYMM%NROT,SYMM%NPTRANS,SYMM%ROTMAP,SYMM%MAGROT,WDES%ISPIN,IO%IU6)
         ! reread k-points with LINVERSION=.FALSE. to generate full mesh
         CALL RE_READ_KPOINTS(KPOINTS,LATT_CUR,.FALSE., &
              T_INFO%NIONS,SYMM%ROTMAP,SYMM%MAGROT,SYMM%ISYM,IO%IU6,IO%IU0)
         CALL KPAR_SYNC_ALL(WDES,W)
         CALL RE_GEN_LAYOUT( GRID, WDES, KPOINTS, LATT_CUR, LATT_CUR,-1, IO%IU0)
         CALL PEAD_RESETUP_WDES(WDES, GRID, KPOINTS, LATT_CUR, LATT_CUR, IO) 
         CALL REALLOCATE_WAVE( W, GRID, WDES, NONL_S, T_INFO, P, LATT_CUR)
      ENDIF
 
      CALL CHECK_FULL_KPOINTS

      ! first thing is to setup permutation matrix 
      ALLOCATE(U(WDES%NB_TOT,WDES%NB_TOT))
      U = zero
     
      ! obtain highest occupied state 
      IHOMO = 1 
      DO ISP = 1, WDES%ISPIN 
         DO NK=1,KPOINTS_FULL%NKPTS
            DO I = 1, WDES%NB_TOT
               IF ( ABS( W%FERTOT(I,NK,ISP))> XI_EMPTY_THRESHHOLD ) THEN
                  IHOMO = MAX( IHOMO, I ) 
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      ! further sanity check, this will probably never happen
      DO I = 1, SIZE(IDXACTIVE ) 
         IF ( IDXACTIVE( I ) > IHOMO ) THEN
            CALL vtutor%error("ERROR: NTARGET_STATES contains entry above Fermi level " // str(IHOMO) &
               // " " // str(IDXACTIVE(I)))
         ENDIF 
      ENDDO
      
      ! set up permutation matrix 
      ALLOCATE(S(WDES%NB_TOT,WDES%NB_TOT))
      S = zero
      DO ISP = IHOMO+1, WDES%NB_TOT 
         S(ISP,ISP)=1
      ENDDO
      ! meddling is done only in the occupied-occupied block
      ! based on IDXACTIVE , the size of IDXACTIVE determines the number 
      ! of active states.  
      ALLOCATE( ISORT( IHOMO - SIZE( IDXACTIVE ) ) )
      ISORT = 0 
      ! look for states not included in IDXACTIVE
      K=0
      DO J = 1, IHOMO
         ! look up if IDXACTIVE contains this state 
         DO I = 1, SIZE(IDXACTIVE) 
            IF ( J == IDXACTIVE( I ) ) EXIT
         ENDDO
         ! check if index is not contained in IDXACTIVE 
         IF ( I > SIZE( IDXACTIVE ) ) THEN
            K = K + 1  
            IF ( K > SIZE( ISORT ) ) THEN
               CALL vtutor%error("ERROR, something went wrong with the non-active space " // str(K) &
                  // " " // str(IHOMO) // " " // str(SIZE(ISORT)))
            ENDIF
            ISORT( K ) = J 
         ENDIF
      ENDDO
      ! one more sanity check 
      IF ( K /= SIZE( ISORT ) ) THEN
         CALL vtutor%error("ERROR, some states have been lost " // str(K) // " " // str(SIZE(ISORT)))
      ENDIF

      ALLOCATE( NSORT( CRPA_NUM_WAN ) )
      ! merge indices
      NSORT( 1:SIZE(ISORT) ) = ISORT(:)
      NSORT( 1+SIZE(ISORT) : CRPA_NUM_WAN ) = IDXACTIVE(:)
      ! first entries are non-active states
      DO I = 1, CRPA_NUM_WAN 
         S( NSORT( I ) , I ) = 1 
      ENDDO

      ALLOCATE(UEFF(WDES%NB_TOT,WDES%NB_TOT))
      UEFF = zero

#ifdef debug
IF ( IO%IU0>=0 ) THEN
WRITE(*,*) REAL(W%CW(10,1:5,1,1))
ENDIF
W%CW(15,2,1,1) = W%CW(15,1,1,1)+0.1
W%CW(15,3,1,1) = W%CW(15,1,1,1)+0.3
W%CW(15,4,1,1) = W%CW(15,1,1,1)+0.4
#endif
      DO ISP=1,WDES%ISPIN
         DO NK=1,KPOINTS_FULL%NKPTS

            CALL DETERMINE_WANNIER_PROJECTOR(WDES,NK,ISP,U(:,:),IO%IU0)

            ! replace remainder in unoccupied space with identity matrix
            DO I = 1+IHOMO, WDES%NB_TOT
               U(I,I) = 1
            ENDDO

            ! sort after projection that is UEFF = S.U 
            CALL GGEMM( 'N', 'N', WDES%NB_TOT, WDES%NB_TOT, WDES%NB_TOT, &
               1._q, U(1,1), WDES%NB_TOT, S(1,1), WDES%NB_TOT, &
               0._q, UEFF(1,1), WDES%NB_TOT)
#ifdef debug
IF ( IO%IU0>=0 ) THEN
WRITE(*,*)NSORT
 WRITE(*,'(I2,35I4)')0,(I, I=1,35)
 WRITE(*,'(A)')'  -----------------------------------------------------------------'
DO I = 1, 48
 WRITE(*,'(I2,"|",48F4.1)')I,UEFF(I,1:48)
ENDDO
ENDIF
#endif
            CALL SETWDES(WDES,WDES1,NK)
            WA=ELEMENTS(W,WDES1,ISP)

            ! redistribute over plane wave coefficients
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF

            ! build linear combinations of the wave functions in
            ! accordance with the effective rotation matrix U
            CALL LINCOM('F',WA%CW_RED,WA%CPROJ_RED,UEFF(1,1), &
           &     WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), &
           &     WDES1%NPL_RED,WDES1%NPRO_RED,WDES%NRPLWV_RED,WDES%NPROD_RED,WDES%NB_TOT, &
           &     WA%CW_RED,WA%CPROJ_RED)

            ! redistribute back over bands
            IF (WDES%DO_REDIS) THEN
               CALL REDIS_PROJ(WDES1, WDES%NBANDS, W%CPROJ(1,1,NK,ISP))
               CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
            ENDIF
         ENDDO !kpoints
      ENDDO !spin

      DEALLOCATE(U)
      DEALLOCATE(S)
      DEALLOCATE(ISORT)
      DEALLOCATE(NSORT)

      IF ((SYMM%ISYM>=0.AND. .NOT.WDES%LGAMMA) .AND. LWRITE_WANPROJ) THEN
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

#ifdef debug
IF ( IO%IU0>=0 ) THEN
!WRITE(*,'(16E10.3)') REAL(W%CW(15,1:16,1,1))*1000
!WRITE(*,'(16E10.3)') REAL(W%CW(15,17:32,1,1))*1000
!WRITE(*,'(16E10.3)') REAL(W%CW(15,33:48,1,1))*1000
WRITE(*,'(5F15.6)') REAL(W%CW(15,1:5,1,1))*1000
ENDIF
#endif
      PROFILING_STOP('rotate_active_occupied_orbitals')
      RETURN
   END SUBROUTINE ROTATE_ACTIVE_OCCUPIED_ORBITALS
 
END MODULE crpa
