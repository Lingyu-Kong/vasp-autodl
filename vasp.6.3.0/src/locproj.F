#include "symbol.inc"
      MODULE locproj
      USE prec
      USE c2f_interface
      USE wave_high

      IMPLICIT NONE

      PUBLIC LPRJ_READER,LPRJ_PROALL,&
             LPRJ_COVL,LPRJ_COVL_IAO,&
             LPRJ_WRITE, LPRJ_WRITE_FUNCTIONS,&
             LPRJ_DEALLOC,LPRJ_LCAO,LPRJ_LDApU,&
             LPRJ_CALC_WANNIER_TRAFO, &
             LPRJ_CALC_WANNIER_TRAFO_QAO ,&
             LPRJ_SVD_LOCALIZE ,&
             LPRJ_ASSIGN_ORBITALS_TO_SITES,&
             !LPRJ_CALC_NONITERATIVE_PM, &
             LPRJ_CALC_INTRINSIC_ATOMIC_ORBITALS, &
             LPRJ_LINCOM, &
             LPRJ_GET_POS, LPRJ_GET_NUM_WAN, &
             LPRJ_GET_SITE, LPRJ_IS_INITIALIZED, LPRJ_IS_CALLED_PROALL, &
             CALC_OVERLAP_GN, LPRJ_ORBITAL_DEFINITIONS, &
             RGRID_PARAMETERS, SETRGRID, &
             RADIAL_HYDROGEN_FUNCTION, RADIAL_GAUSSIAN_FUNCTION, &
             RADIAL_SLATER_FUNCTION, RADIAL_CALC_MOMENT, RADIAL_NORM2, &
             WRITE_WANPROJ_FILE, LPRJ_function, &
             LPRJ_NUM_ORBITALS_ON_ATOM, &
             LPRJ_functions, LPRJ_CALC_SVD, LPRJ_FROM_PSEUDO, LOCPROJ_AUTO 

      PRIVATE
!
!> Data structure that stores the definitions of the local functions
!> on which we will project our Bloch orbitals.
!
      TYPE LPRJ_function
         ! radial part
         INTEGER radial_type    !< 1="optimal" PAW projector, 2=PAW partial wave (PS), 3=Hydrogen-like (=wannier90)
         INTEGER species        !< index of the species at which the local orbital is located
         INTEGER n              !< principal qunatum number
         REAL(q) za             !< real number determining how large the radius of the hydrogenic orbitals should be
         LOGICAL WANNIER90_ORBITAL_DEFINITIONS !< logical governing wether to use wannier90 or LOCPROJ orbital definitions
         ! spherical part
         INTEGER l              !< l angular quantum number
         INTEGER m              !< m magnetic quantum number
         LOGICAL LROTYLM        !< logical to rotate the radial part of the orbitals
         REAL(q) proj_x(3)      !< x axis to use for the rotated projections
         REAL(q) proj_z(3)      !< z axis to use for the rotated projections
         ! site part
         INTEGER poscar_site    !< index of the position in the poscar file
         REAL(q) R(3)           !< position of the localized orbital in reduced coordinates
         ! spin part
         INTEGER spinor         !< spinor index 1 for spin up and 2 for spin down
         REAL(q) spin_qaxis(3)  !< spin quantization axis
         ! related to PS or PR projectors
         INTEGER iproj          !< index of projector or partial wave to be projected on
         INTEGER ibase          !< first pseudopotential channel with angular moment L=l
         INTEGER nproj          !< number of projectors with same L

         !> contains the optimized projectors computed by SPHPRO or SPHPRO_OPTPROJ
         !> and stored in the pseudo datastructure
         REAL(q), POINTER :: optproj(:)

         !> spline representation of the radial part
         REAL(q), ALLOCATABLE :: SPLINE(:,:,:)
      END TYPE LPRJ_function

      !> array containing the information of the localized orbitals specified in the LOCPROJ tag
      TYPE (LPRJ_function), ALLOCATABLE, SAVE :: LPRJ_functions(:)

      INTEGER, PARAMETER :: LMAX=3
      INTEGER, PARAMETER :: NMAX=1000

      !> Contains the projections of the Kohn-Sham orbitals onto i
      !> localized orbitals < \tilde \psi_{nb,k} | S | RYlm_i >.
      !> The dimensions are (NBANDS,NK,ISP,NPROJ)
      GDEF, ALLOCATABLE, SAVE :: LPRJ_COVL(:,:,:,:)
      !> Contains the projections of the Kohn-Sham orbitals onto IAOs
      !> The dimensions are (NBANDS,NK,ISP,NPROJ)
      GDEF, ALLOCATABLE, SAVE :: LPRJ_COVL_IAO(:,:,:,:)

      !> switch on calculaton of IAO?
      LOGICAL, PUBLIC :: LOCPROJ_IAO = .FALSE.

      !> orhtogonalize the IAOs
      LOGICAL        :: LOCPROJ_ORTH_IAO = .TRUE.
      !> lower bound for the Bloch states to consider when building IAOs
      INTEGER        :: LOCPROJ_IAO_NBLOCH_LOW(2) = -1
      !> upper bound for the Bloch states to consider when building IAOs
      INTEGER        :: LOCPROJ_IAO_NBLOCH_HIGH(2) = -1

      !> filename of local projection matrices
      CHARACTER(LEN=7),PUBLIC,PARAMETER :: WANPROJFILE='WANPROJ'

!
!> locproj::locproj_needs_ps  will be set to .TRUE. by locproj::lprj_reader
!! when one or more local functions are of locproj::lprj_function::radial_type=2
!! (pseudo partial waves). This tells VASP to call lcao::lcao_init.
!
      LOGICAL :: LOCPROJ_NEEDS_PS=.FALSE.

      LOGICAL :: LOCPROJ_AUTO = .FALSE.

      !> logical controling wether the augmentation part of the projections should be computed
      LOGICAL :: LOCPROJ_AUG = .TRUE.

      LOGICAL, SAVE  :: IS_INITIALIZED = .FALSE.
      LOGICAL, SAVE  :: IS_CALLED_PROALL = .FALSE.

      CONTAINS


!******************** SUBROUTINE LPRJ_READER ***************************
!
!> Parse the INCAR file looking for occurences of the LOCPROJ or LOCPROJ_AUTO tag.
!> If LOCPROJ is found the LPRJ_FROM_LOCPROJ subroutine is called to parse the string.
!> If LOCPROJ_AUTO is found the LPRJ_FROM_PSEUDO subroutine is called to initialize
!> the functions from the pseudopotential information
!
!***********************************************************************

      SUBROUTINE LPRJ_READER(T_INFO,P,NONCOLLINEAR,IO)
      USE base
      USE poscar
      USE pseudo
      USE full_kpoints
      USE pead, ONLY : PEAD_REQUEST
      USE vaspxml
      USE string, ONLY: str
      USE reader_tags, ONLY: PROCESS_INCAR

      TYPE (type_info) T_INFO
      TYPE (potcar) P(T_INFO%NTYP)
      LOGICAL,INTENT(IN) :: NONCOLLINEAR !< Whether this is a noncollinear calculation or not
      TYPE (in_struct) :: IO

      INTEGER N
      LOGICAL LOPEN
      INTEGER IERR
      INTEGER,PARAMETER :: MAX_LOCPROJ_INPUT=4096
      CHARACTER(LEN=MAX_LOCPROJ_INPUT) :: LOCPROJ_INPUT
      INTEGER :: NUMB, LOCPROJ_SIZE

      REAL(q), PARAMETER :: TINYE=1.E-3_q
      REAL(q), PARAMETER :: TINYO=1.E-1_q
      REAL(q), PARAMETER :: TINYR=1.E-4_q

      ! First parse the file using specialized parser
      LOCPROJ_INPUT=""
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LOCPROJ', LOCPROJ_INPUT, NUMB, IERR, FOUNDNUMBER=LOCPROJ_SIZE)
      ! A newline is strictly required after every command by the flex/bison grammar
      ! to avoid shift/reduce warnings in the syntax.
      ! The user might not explicitely write it, so I add it here
      LOCPROJ_INPUT(LEN_TRIM(LOCPROJ_INPUT)+1:LEN_TRIM(LOCPROJ_INPUT)+1) = NEW_LINE('A')
      IF (IERR==0) CALL LPRJ_FROM_LOCPROJ(T_INFO,P,NONCOLLINEAR,LOCPROJ_INPUT,LPRJ_functions)

      ! Find if automatic projections are requested
      LOCPROJ_AUTO = .FALSE.
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LOCPROJ_AUTO', LOCPROJ_AUTO, IERR)
      IF (LOCPROJ_AUTO) CALL LPRJ_FROM_PSEUDO(T_INFO,P,NONCOLLINEAR,LPRJ_functions)

      ! IAO variables
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LOCPROJ_IAO', LOCPROJ_IAO, IERR, .FALSE.)
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LOCPROJ_ORTH_IAO', LOCPROJ_ORTH_IAO, IERR, .TRUE.)

      ! read lower bound for bloch states
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LOCPROJ_IAO_NBLOCH_LOW', N, IERR, FOUNDNUMBER=NUMB)
      IF (NUMB==0) THEN
        LOCPROJ_IAO_NBLOCH_LOW = -1
      ELSE
        CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LOCPROJ_IAO_NBLOCH_LOW', LOCPROJ_IAO_NBLOCH_LOW, MIN(NUMB,2), IERR)
        IF (NUMB==1) LOCPROJ_IAO_NBLOCH_LOW(2)=LOCPROJ_IAO_NBLOCH_LOW(1)
      ENDIF

      ! read upper bound for bloch states
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LOCPROJ_IAO_NBLOCH_HIGH', N, IERR, FOUNDNUMBER=NUMB)
      IF (NUMB==0) THEN
        LOCPROJ_IAO_NBLOCH_HIGH = -1
      ELSE
        CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LOCPROJ_IAO_NBLOCH_HIGH', LOCPROJ_IAO_NBLOCH_HIGH, MIN(NUMB,2), IERR)
        IF (NUMB==1) LOCPROJ_IAO_NBLOCH_HIGH(2)=LOCPROJ_IAO_NBLOCH_HIGH(1)
      ENDIF

      ! read whether to compute augmentation charges
      LOCPROJ_AUG = .TRUE.
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LOCPROJ_AUG', LOCPROJ_AUG, IERR)

      ! Prepare all that will be needed
      IF (.NOT.ALLOCATED(LPRJ_functions)) RETURN
      IF (LOCPROJ_AUTO) LOCPROJ_NEEDS_PS=.TRUE.
      CALL USE_FULL_KPOINTS
#ifndef gammareal
      CALL PEAD_REQUEST
#endif

      IS_INITIALIZED = .TRUE.
      RETURN
      END SUBROUTINE LPRJ_READER

!****************** SUBROUTINE LPRJ_FROM_LOCPROJ *************************
!
!> @brief Set a list of projections from LOCPORJ tag
!
!! locproj::lprj_reader calls the parser and allocated and fills the
!! locproj::lprj_functions data structure.
!!
!! At the end the following members of LPRJ_functions(i) will have been
!! set (for all functions i in the list):
!!
!! locproj::lprj_function::radial_type:
!!   specifies the type of radial funtion. This can be
!!   radial_type=1  "optimal" PAW projector
!!                  (set by sphpro::sphpro_fast for LORBIT=14),
!!   radial_type=2  PAW pseudo partial wave, or
!!   radial_type=3  Hydrogen-like (Slater type).
!!
!! locproj::lprj_function::poscar_site = n,
!!   where n is the number of the site in the POSCAR file, on which the local
!!   function is centered (if specified in this manner).
!!   Alternatively the position of the local function may have been
!!   specified in terms of direct coordinates of the real space lattice.
!!   In the latter case n=-1.
!!
!! locproj::lprj_function::R: 
!!   specifies the site on which the local function is centered in
!!   direct coordinates of the real space lattice.
!!
!! locproj::lprj_function::l:
!!   l-quantum number of this local function:
!!    l= 0: s,  l= 1: p,   l= 2: d,   l= 3: f,
!!    l=-1: sp, l=-2: sp2, l=-3: sp3, l=-4: sp3d, l=-5: sp3d2
!!   for radial_type/=3 only pure spherical harmonics are used
!!   (i.e., l>=0).
!!
!! locproj::lprj_function::m:
!!   orbital magnetic quantum number of this local function.
!!
!! locproj::lprj_function::species = n
!!   signifies that projectors or pseudo partial waves derived from
!!   the n-th PAW dataset in the POTCAR are to be used to define the
!!   radial part of a particular local function (n=-1 in case
!!   radial_type=3).
!!
!! locproj::lprj_function::n:
!!   in case radial_type/=3: use the "n"-th projector or pseudo partial
!!   wave of angular momentum l.
!!   For radial_type=3: use Slater-type function with n nodes.
!!   Default: locproj::lprj_function::n=1
!!
!! locproj::lprj_function::za:
!!   defines diffusivity of Slater type function (default za=1.0). 
!!   Only relevant for hydrogen-like local functions.
!!
!! locproj::lprj_function::iproj:
!!   the channel number of the PAW projector or pseudo partial wave
!!   to be used as radial part of the local function (only set for
!!   radial_type=1 or 2, otherwise locproj::lprj_function::iproj=-1).
!!   For radial_type=1 or 2 the projector or pseudo partial wave
!!   of the highest lying bound atomic state with angular moment l
!!   will be used. Strictly this is only guaranteed for radial_type=2:
!!   for radial_type=1 VASP will use linear combinations of the
!!   projectors read in from POTCAR in accordance with potcar::optproj
!!   (set in pawm::set_opt_proj in a way that tries to recover the
!!   original order, i.e., the one in potcar::atomic_e).
!!   For LORBIT=14 potcar::optproj is set in sphpro::sphpro_fast, and
!!   only a single "optimal" projector will be available. In that case
!!   the value of iproj is ignored in locproj::lprj_proall below.
!!
!! locproj::lprj_function::ibase:
!!   the channel number of the first PAW projector or pseudo partial wave
!!   with angular moment l (only set for radial_type=1 or 2).
!!
!! locproj::lprj_function::nproj:
!!   number PAW channels with angular moment l
!!   (only set for radial_type=1 or 2).
!
!*************************************************************************
      SUBROUTINE LPRJ_FROM_LOCPROJ(T_INFO,P,NONCOLLINEAR,LOCPROJ_INPUT,LPRJ_functions)
      USE base
      USE poscar
      USE pseudo
      USE vaspxml
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      TYPE (type_info) T_INFO
      TYPE (potcar) P(T_INFO%NTYP)
      LOGICAL,INTENT(IN) :: NONCOLLINEAR !< Whether this is a noncollinear calculation or not
      CHARACTER(LEN=*),INTENT(IN) :: LOCPROJ_INPUT !< String containing the input from LOCPROJ
      TYPE(LPRJ_function), ALLOCATABLE :: LPRJ_functions(:) !< List of functions for project onto

      ! local variables
      TYPE (parser_function), ALLOCATABLE :: functions(:)
      REAL(q) DE,DEOLD,E
      INTEGER N,I,J,IFOUND1,IFOUND2,IFOUND3,ITMP(3)
      LOGICAL LFOUND
      INTEGER ISPINOR, NRSPINOR, IS

      REAL(q), PARAMETER :: TINYE=1.E-3_q
      REAL(q), PARAMETER :: TINYO=1.E-1_q
      REAL(q), PARAMETER :: TINYR=1.E-4_q

      CALL parse_string(N,LEN_TRIM(LOCPROJ_INPUT)+1,LOCPROJ_INPUT)

      IF (N==0) THEN
         CALL free_parser(); RETURN
      ENDIF

      ALLOCATE(functions(N))

      DO I=1,N
         CALL fill_basis_info(functions(I),I)
      ENDDO

      CALL free_parser()

      ! Now we initialize the LPRJ_functions datastructure
      NRSPINOR = 1; IF (NONCOLLINEAR.eqv..TRUE.) NRSPINOR = 2

      IF (ALLOCATED(LPRJ_functions)) DEALLOCATE(LPRJ_functions)
      ALLOCATE(LPRJ_functions(N*NRSPINOR))

      LPRJ_functions(:)%iproj=-1

      CALL XML_TAG("locproj")

      ! validate and copy input to LPRJ_functions
      DO I=1,N
         DO ISPINOR=1,NRSPINOR
         IF (functions(I)%poscar_site/=-1) THEN
            ! a poscar site index has to be in range [1,NIONS]
            IF (functions(I)%poscar_site<1.OR.functions(I)%poscar_site>T_INFO%NIONS) THEN
               CALL vtutor%error("LPRJ_READER: ERROR: site index of function: " // str(I) // " out of &
                  &range: " // str(functions(I)%poscar_site))
            ENDIF
            functions(I)%posx=T_INFO%POSION(1,functions(I)%poscar_site)
            functions(I)%posy=T_INFO%POSION(2,functions(I)%poscar_site)
            functions(I)%posz=T_INFO%POSION(3,functions(I)%poscar_site)
!        ELSE
!           IF (functions(I)%radial_type==1) THEN
!              IF (IU0>=0) WRITE(*,'(A,I4,A)') 'LPRJ_READER: ERROR: function:',I, &
!             &   ': no POSCAR site specified'
!              STOP
!           ENDIF
         ENDIF
         IS = ((I-1)*NRSPINOR+ISPINOR)
         ! copy site information
         LPRJ_functions(IS)%poscar_site=functions(I)%poscar_site
         LPRJ_functions(IS)%R(1)=MOD(functions(I)%posx+100.5_q,1._q)-0.5_q
         LPRJ_functions(IS)%R(2)=MOD(functions(I)%posy+100.5_q,1._q)-0.5_q
         LPRJ_functions(IS)%R(3)=MOD(functions(I)%posz+100.5_q,1._q)-0.5_q
         LPRJ_functions(IS)%nproj=-1

         IF (functions(I)%radial_type==1.OR.functions(I)%radial_type==2) THEN
            IF (functions(I)%l<0) THEN
               CALL vtutor%error("LPRJ_READER ERROR: function: " // str(I) //&
                    " is not a simple spherical harmonic L= " // str(functions(I)%l))
            ENDIF
            IF (functions(I)%species==-1) THEN
               IF (functions(I)%poscar_site==-1) THEN
                  CALL vtutor%error("LPRJ_READER ERROR: both POSCAR as well as POTCAR index of &
                     &function: " // str(I) // " not specified")
               ENDIF
               functions(I)%species=T_INFO%ITYP(functions(I)%poscar_site)  
            ENDIF
            IF (functions(I)%species<1.OR.functions(I)%species>T_INFO%NTYP) THEN
               CALL vtutor%error("LPRJ_READER: POTCAR index of function: " // str(I) //&
                    " out of range: " // str(functions(I)%species))
            ENDIF

            IFOUND1=0; IFOUND2=0
            DO J=1,P(functions(I)%species)%LMAX
               IF (P(functions(I)%species)%LPS(J)==functions(I)%l .AND. IFOUND1==0) IFOUND1=J
               IF (P(functions(I)%species)%LPS(J)==functions(I)%l) IFOUND2=IFOUND2+1
            ENDDO
           
            IF (IFOUND2==0) THEN
               CALL vtutor%error("LPRJ_READER: function: " // str(I) // " no state with L= "//&
                    str(functions(I)%l) // " on species " // TRIM(P(functions(I)%species)%ELEMENT))
            ENDIF
            
            LPRJ_functions(IS)%ibase=IFOUND1          ! first index of projectors
            LPRJ_functions(IS)%nproj=IFOUND2          ! total number of projectors

            ALLOCATE(LPRJ_functions(IS)%optproj(IFOUND2)); LPRJ_functions(IS)%optproj=0

            IF (functions(I)%n/=-1) THEN
               ! Take the n-th projector with angular moment l (if it exists, of course).
               IFOUND1=0; LFOUND=.FALSE.
               DO J=1,P(functions(I)%species)%LMAX
                  IF (P(functions(I)%species)%LPS(J)==functions(I)%l) IFOUND1=IFOUND1+1
                  IF (IFOUND1==functions(I)%n) THEN
                     LFOUND=.TRUE.; IFOUND2=J; EXIT
                  ENDIF
               ENDDO
               IF (.NOT.LFOUND) THEN
                  CALL vtutor%error("LPRJ_READER: function: " // str(I) // " no " // &
                     str(functions(I)%n) // "-th projector, with L=" // str(functions(I)%l) // &
                     " on species " // TRIM(P(functions(I)%species)%ELEMENT))
               ELSE
                  LPRJ_functions(IS)%iproj=IFOUND2
               ENDIF
            ELSE
               ! Take the projector that corresponds to the highest lying bound state with
               ! angular momentum l.
               ! First identify the corresponding atomic state.
               LFOUND=.FALSE.; IFOUND1=0; E=-1.E+6_q
               DO J=1,SIZE(P(functions(I)%species)%ATOMIC_N)
                  IF (P(functions(I)%species)%ATOMIC_L(J)==functions(I)%l &
                 &   .AND.P(functions(I)%species)%ATOMIC_OCC(J)>TINYO &
                 &   .AND.P(functions(I)%species)%ATOMIC_E(J)>E) THEN
                     LFOUND=.TRUE.; IFOUND1=J; E=P(functions(I)%species)%ATOMIC_E(J)
                  ENDIF
               ENDDO
               IF (.NOT.LFOUND) THEN
                  CALL vtutor%error("LPRJ_READER: function: " // str(I) // " no state with L=" //&
                        str(functions(I)%l) // " on species " // TRIM(P(functions(I)%species)%ELEMENT))
               ENDIF
               ! find the highest valence projector
               IFOUND2=0; IFOUND3=0; DEOLD=1.E+6_q
               DO J=1,P(functions(I)%species)%LMAX
                  IF (P(functions(I)%species)%LPS(J)==functions(I)%l) IFOUND3=IFOUND3+1
                  DE=ABS(P(functions(I)%species)%ATOMIC_E(IFOUND1)-P(functions(I)%species)%E(J))
                  IF (DE<DEOLD) THEN
                     IFOUND2=J; DEOLD=DE
                     functions(I)%n=IFOUND3
                  ENDIF
               ENDDO            
               IF (DEOLD>TINYE) THEN
                  CALL vtutor%error("LPRJ_READER: function: " // str(I) // " no projector for &
                     &highest valence state with L=" // str(functions(I)%l) // " on species " // &
                     str(functions(I)%poscar_site))
               ELSE
                  LPRJ_functions(IS)%iproj=IFOUND2
               ENDIF
            ENDIF
         ENDIF

         ! this tells VASP to call LCAO_INIT
         IF (functions(I)%radial_type==2) LOCPROJ_NEEDS_PS=.TRUE.

         LPRJ_functions(IS)%radial_type=functions(I)%radial_type
         LPRJ_functions(IS)%WANNIER90_ORBITAL_DEFINITIONS=.FALSE.

         LPRJ_functions(IS)%species=functions(I)%species

         LPRJ_functions(IS)%n=functions(I)%n
         LPRJ_functions(IS)%l=functions(I)%l
         LPRJ_functions(IS)%m=functions(I)%m
         LPRJ_functions(IS)%spinor=ISPINOR
         ! TODO: this should be read from the input file, should be done when we modify the parser
         LPRJ_functions(IS)%spin_qaxis = [0.0_q,0.0_q,1.0_q]
         LPRJ_functions(IS)%LROTYLM = .FALSE.

         LPRJ_functions(IS)%za=functions(I)%sigma


         ITMP(1)=LPRJ_functions(IS)%n; ITMP(2)=LPRJ_functions(IS)%l; ITMP(3)=LPRJ_functions(IS)%m

         CALL XML_TAG("function")
         CALL XML_TAG_INT("radial_type",LPRJ_functions(IS)%radial_type)
         CALL XML_TAG_INT("species",LPRJ_functions(IS)%species)
         CALL XML_VEC_INT(ITMP,"nlm")
         CALL XML_TAG_REAL("za",LPRJ_functions(IS)%za)
         CALL XML_TAG_INT("poscar_site",LPRJ_functions(IS)%poscar_site)
         CALL XML_VEC_REAL(LPRJ_functions(IS)%R,"position")
         CALL XML_CLOSE_TAG("function")


!        IF (IU0>=0) WRITE(*,'(I4,X,"n=",I2,X,"l=",I2,X,"m=",I2,X,"rtype=",I2,X,"species=",I2,X,"za=",F14.7,X,"site=",I4,X,"R=",3F14.7,X,"ipr=",I2)') &
!       & I,LPRJ_functions(IS)%n,LPRJ_functions(IS)%l,LPRJ_functions(IS)%m,LPRJ_functions(IS)%radial_type, &
!       &   LPRJ_functions(IS)%species,LPRJ_functions(IS)%za,LPRJ_functions(IS)%poscar_site,LPRJ_functions(IS)%R, &
!       &   LPRJ_functions(Is)%iproj
         ENDDO
      ENDDO

      CALL XML_CLOSE_TAG("locproj")

      DEALLOCATE(functions)

      END SUBROUTINE LPRJ_FROM_LOCPROJ


!******************** SUBROUTINE LPRJ_FROM_PSEUDO ***************************
!
!> @brief Set a list of projections from the pseudos
!
!***********************************************************************
      SUBROUTINE LPRJ_FROM_PSEUDO(T_INFO,P,NONCOLLINEAR,LPRJ_functions,NUM_WANN)
      USE base
      USE poscar
      USE pseudo
      USE lattice
      USE core_rel, ONLY:IS_BOUND

      TYPE (type_info) T_INFO
      TYPE (potcar) P(T_INFO%NTYP)
      LOGICAL,INTENT(IN) :: NONCOLLINEAR !< Whether this is a noncollinear calculation or not
      TYPE(LPRJ_function), ALLOCATABLE :: LPRJ_functions(:) !< List of functions for project onto
      INTEGER,OPTIONAL :: NUM_WANN !< If NUM_WANN is present don't use LPRJ_functions, only return its size in NUM_WANN

      ! Local variables
      TYPE(LPRJ_function), ALLOCATABLE :: LPRJ_functions_tmp(:)
      REAL(q), PARAMETER :: LOWBOUND_E=-80._q
      REAL(q), PARAMETER :: TINYO=1.E-1_q
      REAL(q), PARAMETER :: TINYR=1.E-4_q
      INTEGER :: I,NP,ITYP,L,ISPINOR,M,NRSPINORS,STEP,CHANNEL
      INTEGER :: NCHANNELS_REF(0:3), NCHANNELS(0:3)

      NRSPINORS = 1; IF (NONCOLLINEAR) NRSPINORS = 2

      STEP = NRSPINORS*T_INFO%NIONS*3 !initial guess for number of projectors
      IF (.NOT.PRESENT(NUM_WANN)) THEN
         IF (ALLOCATED(LPRJ_functions)) DEALLOCATE(LPRJ_functions)
         ALLOCATE(LPRJ_functions(STEP))
      ENDIF

      ! Determine projections from the pseudos
      NP = 0
      DO I=1,T_INFO%NIONS
         ITYP = T_INFO%ITYP(I)
         ! Get all available projectors for this type
         NCHANNELS_REF = 0
         !WRITE(*,*) 'all'
         DO CHANNEL=1,P(ITYP)%LMAX
            IF (P(ITYP)%QATO(CHANNEL,CHANNEL) == 0._q )  CYCLE ! only bounds states
            L = P(ITYP)%LPS(CHANNEL)
            NCHANNELS_REF(L) = NCHANNELS_REF(L)+1
            !WRITE(*,*) CHANNEL, L, NCHANNELS_REF(L)
         ENDDO
         ! Use the projectors with largest N
         NCHANNELS = 0
         !WRITE(*,*) 'selected'
         DO CHANNEL=1,P(ITYP)%LMAX
            IF (P(ITYP)%QATO(CHANNEL,CHANNEL) == 0._q )  CYCLE ! only bound states
            L = P(ITYP)%LPS(CHANNEL)
            NCHANNELS(L) = NCHANNELS(L)+1
            IF (NCHANNELS(L) > 1) CYCLE ! only take smallest N
            !WRITE(*,*) CHANNEL, L, NCHANNELS(L)
            DO M=1,L*2+1
               DO ISPINOR=1,NRSPINORS
                  NP = NP+1
                  IF (PRESENT(NUM_WANN)) CYCLE
                  ! if size is not enough, resize
                  IF (NP>SIZE(LPRJ_functions)) THEN
                     ALLOCATE(LPRJ_functions_tmp(SIZE(LPRJ_functions)+STEP))
                     LPRJ_functions_tmp(:SIZE(LPRJ_functions)) = LPRJ_functions
                     CALL MOVE_ALLOC(LPRJ_functions_tmp,LPRJ_functions)
                  ENDIF
                  LPRJ_functions(NP)%R = T_INFO%POSION(:,I)
                  LPRJ_functions(NP)%poscar_site = I
                  LPRJ_functions(NP)%nproj = -1
                  LPRJ_functions(NP)%n = NCHANNELS(L)
                  LPRJ_functions(NP)%l = L
                  LPRJ_functions(NP)%m = M
                  LPRJ_functions(NP)%iproj = CHANNEL
                  LPRJ_functions(NP)%za = 1 ! this is hardcoded in the parser as well!
                  LPRJ_functions(NP)%LROTYLM = .FALSE.
                  LPRJ_functions(NP)%species = ITYP
                  LPRJ_functions(NP)%spinor = ISPINOR
                  LPRJ_functions(NP)%spin_qaxis = [0.0_q,0.0_q,1.0_q]
                  LPRJ_functions(NP)%radial_type = 2 ! 1="optimal" PAW projector, 2=PAW partial wave (PS), 3=Hydrogen-like (=wannier90)
                  LPRJ_functions(NP)%WANNIER90_ORBITAL_DEFINITIONS=.FALSE.
               END DO
            END DO
         END DO
      END DO

      IF (PRESENT(NUM_WANN)) THEN
         NUM_WANN = NP
         RETURN
      ENDIF

      ! remove any free space
      ALLOCATE(LPRJ_functions_tmp(NP))
      LPRJ_functions_tmp = LPRJ_functions(:NP)
      CALL MOVE_ALLOC(LPRJ_functions_tmp,LPRJ_functions)

      !DO N=1,SIZE(LPRJ_functions)
      !   WRITE(*,*) N,LPRJ_functions(N)%R,LPRJ_functions(N)%L,LPRJ_functions(N)%M
      !END DO

      IF (ANY(LPRJ_functions(:)%radial_type == 2)) LOCPROJ_NEEDS_PS=.TRUE.

      END SUBROUTINE LPRJ_FROM_PSEUDO

!******************** SUBROUTINE LPRJ_PROALL ***************************
!
!> Compute the projections of the bloch states onto a list of
!> local functions specified in the LOCPROJ tag.
!> The list of functions is stored in the locproj::lprj_functions
!> datastructure.
!>
!> This routine uses a different strategy than the previous implementation
!> Instead of computing a trial orbital and projecting all the bloch states
!> (every spin, band and kpoint) onto it, here we compute a wavefuna object containing
!> the trial WFs and the compute the overlap.
!> This allows to distribute the generation of trial orbitals over bands
!>
!> @details @ref openacc :
!> On entering this subroutine W%CW and W%CPROJ are updated on the host
!> (if present on the device) and then OpenACC execution is switched off.
!> On exit the OpenACC execution mode reverts to its original status.
!
!***********************************************************************
      SUBROUTINE LPRJ_PROALL(W,WDES,P,CQIJ,LATT_CUR,LPRJ_functions,LPRJ_COVL,T_INFO,INFO,IO,VKPT,W_TRIALFUNC)

#ifdef _OPENACC
      USE mopenacc
#endif
      USE base
      USE pseudo
      USE mgrid
      USE wave_high
      USE full_kpoints
      USE constant
      USE lattice
      USE poscar
      USE nonl_high
      USE pead, ONLY : LUSEPEAD,WDES_FULL_PEAD,PEAD_WA_ROTATE
      USE choleski, ONLY: CALC_PAW_OVERLAP

      TYPE (wavespin) W
      TYPE (wavedes) WDES
      TYPE (potcar) P(:)
      TYPE (latt) LATT_CUR
      TYPE (LPRJ_function),ALLOCATABLE :: LPRJ_functions(:)
      GDEF, ALLOCATABLE :: LPRJ_COVL(:,:,:,:)
      TYPE (type_info) T_INFO
      TYPE (info_struct) INFO
      OVERLAP CQIJ(:,:,:,:)
      TYPE (in_struct) :: IO
      REAL(q),INTENT(IN),OPTIONAL,TARGET :: VKPT(:,:) !< list of k-points in which to compute orbital projections
      !> Contains the localized orbitals as constructed by
      !> the LPRJ_GEN_TRIALORBITALS routine according to the LOCPROJ tag.
      !> These are needed for the CALC_INTRINSIC_ATOMIC_ORBITALS routine.
      TYPE(wavespin),OPTIONAL :: W_TRIALFUNC
      ! Local variables
      TYPE(wavespin) :: WP
      TYPE(wavedes), POINTER :: WDES_TRIALFUNC
      TYPE(wavedes1), TARGET :: WDESK
      TYPE(wavedes1), TARGET :: WDES1_TRIALFUNC
      TYPE(wavefuna) :: WA_TRIALFUNC, WK
      TYPE(wavefuna) :: WA_TRIALFUNC_SAVE
      TYPE(nonl_struct)  :: NONL_S !< temporary structure to compute the non-local part
      TYPE(nonlr_struct) :: NONLR_S!< temperary structure to compute the non-local part (real space)
      INTEGER, ALLOCATABLE :: LPRJ_functions_index(:)
      REAL(q), POINTER :: P_VKPT(:,:)
      INTEGER :: ISP, IKPT
      INTEGER :: NKPTS, NK
      INTEGER :: LMDIM
      INTEGER :: I, IFUNC, IFUNC_SAVE, IFUNC_LOCAL
      INTEGER :: NBANDS, NPAR
      INTEGER :: NLPRJ_functions, NLPRJ_functions_batch
      LOGICAL :: LREALLOCATE
      REAL(q) :: KPT(3)
      GDEF,ALLOCATABLE :: LPRJ_COVL_LOCAL(:,:)

      PROFILING_START('lprj_proall')

!$ACC UPDATE SELF(W%CW,W%CPROJ) IF_PRESENT IF(ACC_EXEC_ON)
      PUSH_ACC_EXEC_ON(.FALSE.)

      WP=W
      ! Use the PEAD WF descriptor
      IF (LUSEPEAD()) WP%WDES=>WDES_FULL_PEAD

      NLPRJ_functions = SIZE(LPRJ_functions)

      CALL CHECK_FULL_KPOINTS

      NKPTS = KPOINTS_FULL%NKPTS
      P_VKPT => KPOINTS_FULL%VKPT
      IF (PRESENT(VKPT)) THEN
         NKPTS = SIZE(VKPT,2)
         P_VKPT => VKPT
      ENDIF

      NONLR_S%LREAL = INFO%LREAL
      IF (INFO%LREAL) THEN
#ifndef shmem_rproj
         CALL NONLR_SETUP(NONLR_S,T_INFO,P,INFO%LREAL,WDES%LSPIRAL)
#else
         CALL NONLR_SETUP(NONLR_S,T_INFO,P,INFO%LREAL,WDES%LSPIRAL,WDES%COMM_intra_node)
#endif
         CALL REAL_OPTLAY(WDES%GRID,LATT_CUR,NONLR_S,.TRUE.,LREALLOCATE, -1,-1)
         CALL NONLR_ALLOC(NONLR_S)
      ELSE
         CALL NONL_ALLOC(NONL_S, T_INFO, P, WP%WDES, .FALSE.)
      ENDIF

      ! allocate wavespin to store LPRJ_functions if desired
      IF (PRESENT(W_TRIALFUNC)) THEN
         ALLOCATE(WDES_TRIALFUNC)
         WDES_TRIALFUNC = WP%WDES
         NPAR = WDES%NB_PAR
         NBANDS = NLPRJ_functions
         NBANDS = ((NBANDS+NPAR-1)/NPAR)*NPAR
         WDES_TRIALFUNC%NB_TOT = NBANDS
         WDES_TRIALFUNC%NBANDS=WDES_TRIALFUNC%NB_TOT/NPAR
         CALL ALLOCW(WDES_TRIALFUNC, W_TRIALFUNC)
         W_TRIALFUNC%CW = cmplx(0.0_q,0.0_q,q)
         W_TRIALFUNC%CPROJ = zero
      ENDIF

      ! Allocations
      IF (ALLOCATED(LPRJ_COVL)) DEALLOCATE(LPRJ_COVL)
      ALLOCATE(LPRJ_COVL(WDES%NB_TOT,NKPTS,WDES%ISPIN,NLPRJ_functions))
      LPRJ_COVL=zero
      ALLOCATE(LPRJ_functions_index(WDES%NB_TOT))
      ALLOCATE(LPRJ_COVL_LOCAL(WDES%NB_TOT,WDES%NB_TOT))

      LMDIM = SIZE(CQIJ,1)
      DO ISP=1,WDES%ISPIN
         DO IKPT=1,NKPTS
#ifdef MPI
            IF (MOD(IKPT-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
            KPT = P_VKPT(:,IKPT)

            ! search for kpoint k in BZ
            NK=KPOINT_IN_FULL_GRID(KPT,KPOINTS_FULL)
            CALL SETWDES(WP%WDES,WDESK,NK)
            ! rotate WF to k-point
            IF (NK==KPOINTS_FULL%NEQUIV(NK)) THEN
               ! k is a kpoint in the IBZ
               WK=ELEMENTS(WP,WDESK,ISP)
            ELSE
               ! k is not a kpoint in the IBZ
               CALL NEWWAVA(WK,WDESK,WDESK%NBANDS)
               CALL PEAD_WA_ROTATE(WP,P,LATT_CUR,ISP,WK)
            ENDIF

            IF (INFO%LREAL) THEN
               CALL RSPHER(WDES%GRID,NONLR_S,LATT_CUR)
               CALL PHASER(WDES%GRID,LATT_CUR,NONLR_S,NK,WP%WDES)
            ELSE
               CALL SPHER(WDES%GRID, NONL_S, P, WP%WDES, LATT_CUR, 1, NK)
               CALL PHASE(WP%WDES, NONL_S, NK)
            ENDIF

            ! Allocate WA_TRIALFUNC
            CALL NEWWAVA(WA_TRIALFUNC,WDESK,WDESK%NBANDS)
            WA_TRIALFUNC%CW = cmplx(0.0_q,0.0_q,q)
            WA_TRIALFUNC%CPROJ = zero

            ! Compute the overlap of blocks of local projections with the wavefunctions
            NLPRJ_functions_batch=0
            IFUNC_SAVE=0
            DO IFUNC=1,NLPRJ_functions
               ! if the batch of functions is not complete then add this one
               IF (NLPRJ_functions_batch < WDESK%NB_TOT) THEN
                  NLPRJ_functions_batch = NLPRJ_functions_batch + 1
                  LPRJ_functions_index(NLPRJ_functions_batch) = IFUNC
               ENDIF
               ! if the batch of functions is complete or last k-point then compute
               IF (NLPRJ_functions_batch==WDESK%NB_TOT .OR. IFUNC==NLPRJ_functions) THEN
                  ! generate trial orbitals at k-point
                  CALL LPRJ_GEN_TRIALORBITALS(WA_TRIALFUNC,P,CQIJ,LATT_CUR,&
                               NLPRJ_functions_batch,LPRJ_functions_index,LPRJ_functions,&
                               NONL_S,NONLR_S,INFO,IO)

                  ! compute overlap
                  CALL CALC_PAW_OVERLAP(WK,WA_TRIALFUNC,CQIJ,LPRJ_COVL_LOCAL)
                  DO I=1,NLPRJ_functions_batch
                     IFUNC_LOCAL = LPRJ_functions_index(I)
                     LPRJ_COVL(:,IKPT,ISP,IFUNC_LOCAL) = LPRJ_COVL_LOCAL(:,I)
                  ENDDO

                  ! save orbitals from IAOs
                  IF (PRESENT(W_TRIALFUNC)) THEN
                     ! get WA from W_TRIALFUNC where the WF will be saved
                     CALL SETWDES(W_TRIALFUNC%WDES,WDES1_TRIALFUNC,IKPT)
                     WA_TRIALFUNC_SAVE = ELEMENTS(W_TRIALFUNC,WDES1_TRIALFUNC,ISP)
                     DO I=1,WA_TRIALFUNC%WDES1%NBANDS
                        IFUNC_SAVE=IFUNC_SAVE+1
                        IF (IFUNC_SAVE>WA_TRIALFUNC_SAVE%WDES1%NBANDS) CYCLE
                        CALL W1_COPY(ELEMENT(WA_TRIALFUNC,I),&
                                     ELEMENT(WA_TRIALFUNC_SAVE,IFUNC_SAVE))
                     ENDDO
                  ENDIF

                  !Reset the batch size
                  NLPRJ_functions_batch = 0
               ENDIF
            ENDDO ! loop over projections

            CALL DELWAVA(WA_TRIALFUNC)
         ENDDO ! loop over bands
      ENDDO ! loop over spins
      CALLMPI( M_sum_g( WDES%COMM_KINTER, LPRJ_COVL, SIZE(LPRJ_COVL)))

      ! free splines
      DO IFUNC=1,NLPRJ_functions
         IF(ALLOCATED(LPRJ_functions(IFUNC)%SPLINE)) THEN
            DEALLOCATE(LPRJ_functions(IFUNC)%SPLINE)
         ENDIF
      ENDDO

      IS_CALLED_PROALL = .TRUE.
      IF (INFO%LREAL) THEN
         CALL NONLR_DEALLOC(NONLR_S)
      ELSE
         CALL NONL_DEALLOC(NONL_S)
      ENDIF
      DEALLOCATE(LPRJ_functions_index)
      DEALLOCATE(LPRJ_COVL_LOCAL)

      POP_ACC_EXEC_ON

      PROFILING_STOP('lprj_proall')

      END SUBROUTINE LPRJ_PROALL

!******************** SUBROUTINE LPRJ_GEN_TRIALORBITALS *************
!> Given a k-point compute the localized orbitals of the WFs
!> At the end of this routine the SPLINE varaible of each LPRJ_functions
!> will be allocated. This is done so that consecutive calls don't
!> require multiple evaluations of the SPLINE.
!> At the end the SPLINES should be deallocated for each LPRJ_functions
!********************************************************************

      SUBROUTINE LPRJ_GEN_TRIALORBITALS(WA_TRIALFUNC,P,CQIJ,LATT_CUR,NFUNCS,LPRJ_functions_index,LPRJ_functions,NONL_S,NONLR_S,INFO,IO)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE base
      USE pseudo
      USE mgrid
      USE fileio
      USE wave_high
      USE full_kpoints
      USE constant
      USE lattice
      USE poscar
      USE nonl_struct_def, ONLY: nonl_struct
      USE nonlr_struct_def, ONLY: nonlr_struct
      USE lcao, ONLY : ATOM_LCAO

      TYPE (wavefuna) :: WA_TRIALFUNC
      TYPE (potcar) P(:)
      TYPE (latt) LATT_CUR
      INTEGER :: NFUNCS
      INTEGER :: LPRJ_functions_index(NFUNCS)
      TYPE (LPRJ_function) :: LPRJ_functions(:)
      TYPE (info_struct) :: INFO
      TYPE (nonl_struct) :: NONL_S
      TYPE (nonlr_struct) :: NONLR_S
      OVERLAP CQIJ(:,:,:,:)
      TYPE (in_struct) :: IO

      ! local variables
      INTEGER IFNC,IFNC_TMP,ISPINOR,I,D,IFNC_LOCAL
      INTEGER LSTART,LSTOP
      INTEGER N,L,M,IPROJ,IBASE,ITYPE,ISPEC,NPROJ
      INTEGER IU0,IU6
      INTEGER LTMP(4)
      REAL(q) R(3),ZA
      REAL(q) PROJ_X(3),PROJ_Z(3)
      REAL(q) SPIN_QAXIS(3)
      REAL(q) HYBRID_ORBITAL((LMAX+1)**2)
      INTEGER, ALLOCATABLE :: POSSIBLE_ENTRIES(:)
      LOGICAL LNORM,LAUG,LROTYLM
      LOGICAL LWANNIER90_ORBITAL_DEFINITIONS
      TYPE(wavefun1) :: W1
      REAL(q), ALLOCATABLE :: PSPNL_SPLINE(:,:,:,:)

      PROFILING_START('lprj_gen')

      IU0 = IO%IU0
      IU6 = IO%IU6

      !IF (IU0>=0) WRITE(*,*) 'Computing projections with LOCPROJ'
      funcs: DO IFNC_TMP=1,NFUNCS
         IFNC = LPRJ_functions_index(IFNC_TMP)
         N=LPRJ_functions(IFNC)%n
         L=LPRJ_functions(IFNC)%l
         M=LPRJ_functions(IFNC)%m
         R=LPRJ_functions(IFNC)%R

         IBASE=LPRJ_functions(IFNC)%ibase
         IPROJ=LPRJ_functions(IFNC)%iproj
         ISPEC=LPRJ_functions(IFNC)%species
         ITYPE=LPRJ_functions(IFNC)%radial_type
         NPROJ=LPRJ_functions(IFNC)%nproj
         ISPINOR=LPRJ_functions(IFNC)%spinor
         SPIN_QAXIS=LPRJ_functions(IFNC)%spin_qaxis
         LROTYLM = LPRJ_functions(IFNC)%LROTYLM
         PROJ_X = LPRJ_functions(IFNC)%PROJ_X
         PROJ_Z = LPRJ_functions(IFNC)%PROJ_Z
         LWANNIER90_ORBITAL_DEFINITIONS = LPRJ_functions(IFNC)%WANNIER90_ORBITAL_DEFINITIONS

         ZA=LPRJ_functions(IFNC)%za

         ! angular types
         IF (LWANNIER90_ORBITAL_DEFINITIONS) THEN
            CALL WANNIER90_ORBITAL_DEFINITIONS(L,M,HYBRID_ORBITAL)
         ELSE
            CALL LPRJ_ORBITAL_DEFINITIONS(L,M,HYBRID_ORBITAL)
         ENDIF

         ! radial types
         SELECT CASE(ITYPE)

         ! project onto PAW projector or "optimal" PAW projector (LORBIT=14)
         CASE(1)

            ! slightly paranoid tests
            IF (IPROJ-IBASE<0) THEN
               CALL vtutor%bug("internal error in: LPRJ_GEN_TRIALORBITALS IBASE and IPROJ incompatible " // &
                  str(IBASE-IPROJ), __FILE__, __LINE__)
            ENDIF
            IF (IPROJ-IBASE>=NPROJ) THEN
               CALL vtutor%bug("internal error in: LPRJ_GEN_TRIALORBITALS IPROJ too large " // str(IBASE-IPROJ) &
                  // " " // str(NPROJ), __FILE__, __LINE__)
            ENDIF
            ! this condition can only occur if SPHPRO_FAST has reduced the number of available
            ! vectors to 1 (LORBIT=14); then we pick the only available projector
            ! (unfortunately this information can not be made available to LPRJ_READER)
            IF (IPROJ-IBASE>=SIZE(P(ISPEC)%OPTPROJ,2)) THEN
               LPRJ_functions(IFNC)%iproj=IBASE
               LPRJ_functions(IFNC)%n=1
               IPROJ=IBASE
            ENDIF
            IF (.NOT.ALLOCATED(LPRJ_functions(IFNC)%SPLINE)) THEN
               ! calculate non-local projectors (rec. space spline fits)
               CALL CALC_PSPNL_SPLINE(P,PSPNL_SPLINE)
               ! copy relevant one to SPLINE
               ALLOCATE(LPRJ_functions(IFNC)%SPLINE(1:NPSNL,5,LMAX+1))
               ! copy positions from first projector
               LPRJ_functions(IFNC)%SPLINE(1:NPSNL,1,L+1)=PSPNL_SPLINE(1:NPSNL,1,IBASE,ISPEC)
               ! initialize rest
               LPRJ_functions(IFNC)%SPLINE(1:NPSNL,2:,L+1)=0
               !IF (IU0>=0) WRITE(*,'("using NT=",I4," L=",I4," IBASE=",I4," IPROJ=",I4,4F14.7)') &
               !     ISPEC, L, IBASE, IPROJ-IBASE+1, P(ISPEC)%OPTPROJ(1:NPROJ,IPROJ-IBASE+1,L+1)
               ! now linear combine the projectors
               DO I=1,NPROJ
                  LPRJ_functions(IFNC)%SPLINE(1:NPSNL,2:,L+1)=&
                  LPRJ_functions(IFNC)%SPLINE(1:NPSNL,2:,L+1)+&
                  PSPNL_SPLINE(1:NPSNL,2:,IBASE+I-1,ISPEC)*P(ISPEC)%OPTPROJ(I,IPROJ-IBASE+1,L+1)
               ENDDO
               DEALLOCATE(PSPNL_SPLINE)
            ENDIF
            LPRJ_functions(IFNC)%optproj(1:NPROJ)=P(ISPEC)%OPTPROJ(1:NPROJ,IPROJ-IBASE+1,L+1)
            ! not S orthogonal
            LSTART=L; LSTOP=LSTART; LNORM=.FALSE.; LAUG=.FALSE.
            ! make them S orthogonal
!           LSTART=L; LSTOP=LSTART; LNORM=.TRUE.; LAUG=.TRUE.

         ! project onto PAW partial wave (pseudized)
         CASE(2)

            IF (.NOT.ALLOCATED(LPRJ_functions(IFNC)%SPLINE)) THEN
               CALL FIND_ATOM_LCAO_ENTRY(ISPEC,L,IPROJ,POSSIBLE_ENTRIES,I)
               IF (I==-1) THEN
                CALL vtutor%error('Could not find projector '//str(IPROJ)//&
                  ' with angular momentum l='//str(L)//' for species-number '//str(ISPEC)//&
                  ' with symbol '//TRIM(TRIM(P(ISPEC)%ELEMENT))//' possible options are '//&
                  str(POSSIBLE_ENTRIES(:)))
               ENDIF
               ! copy relevant entry of ATOM_LCAO%WFCT_G to SPLINE
               ALLOCATE(LPRJ_functions(IFNC)%SPLINE(SIZE(ATOM_LCAO(ISPEC)%WFCT_G,1),5,LMAX+1))
               LPRJ_functions(IFNC)%SPLINE(:,:,L+1)=ATOM_LCAO(ISPEC)%WFCT_G(:,:,I)
            ENDIF
            LSTART=L; LSTOP=LSTART; LNORM=.TRUE.; LAUG=.TRUE.

         ! project onto Hydrogen-like functions (like in the VASP2WANNIER90 interface)
         CASE(3)

            IF (.NOT.ALLOCATED(LPRJ_functions(IFNC)%SPLINE)) THEN
               ALLOCATE(LPRJ_functions(IFNC)%SPLINE(NMAX,5,LMAX+1))
               ! calculate hydrogen like functions (rec. space spline fits)
               CALL CALC_HYDRO_SPLINE(N,ZA,LMAX,INFO%ENMAX,LPRJ_functions(IFNC)%SPLINE,IO)
            ENDIF
            ! based on which HYBRID_ORBITAL are selected we determine LSTART and LSTOP
            LTMP = 0
            IF(ANY(HYBRID_ORBITAL(1:1)/=0))   LTMP(1)=0
            IF(ANY(HYBRID_ORBITAL(2:4)/=0))   LTMP(2)=1
            IF(ANY(HYBRID_ORBITAL(5:9)/=0))   LTMP(3)=2
            IF(ANY(HYBRID_ORBITAL(10:16)/=0)) LTMP(4)=3
            LSTART=MINVAL(LTMP); LSTOP=MAXVAL(LTMP); LNORM=.TRUE.; LAUG=.TRUE.
         CASE DEFAULT
            CALL vtutor%error("LPRJ_GEN_TRIALORBITALS: ERROR: function: " // str(IFNC) // "is of an unknown &
               &radial type: " // str(ITYPE))
         END SELECT

         ! check if IFNC means work for this node? (parallelization)
         D = IFNC_TMP-WA_TRIALFUNC%WDES1%NB_LOW
         IF( (D .GE. 0) .AND. (MOD(D,WA_TRIALFUNC%WDES1%NB_PAR) .EQ. 0)) THEN
           IFNC_LOCAL = 1 + D/WA_TRIALFUNC%WDES1%NB_PAR
           ! create wavefun1 and store in WA_TRIALFUNC as band number IFNC
           W1 = ELEMENT(WA_TRIALFUNC, IFNC_LOCAL)
           CALL W1_FROM_LPRJ_function(LSTART, LSTOP, LPRJ_functions(IFNC)%SPLINE, R, ISPINOR,   &
                                      SPIN_QAXIS, LROTYLM, PROJ_X, PROJ_Z, CQIJ, LATT_CUR, NONL_S, NONLR_S, &
                                      LNORM, LAUG, HYBRID_ORBITAL, W1)
         ENDIF

      ENDDO funcs

      PROFILING_STOP('lprj_gen')

      END SUBROUTINE LPRJ_GEN_TRIALORBITALS

!******************** SUBROUTINE LPRJ_CALC_WANNIER_TRAFO **********************
!
! Calculates the unitary matrix U which transforms from the Bloch basis into
! the orthonormalized Wannier basis for all k-points. This is done using a
! singular-value decomposition of the Overlap matrix < Psi_m | g_n > between
! Bloch state Psi_m and localized 'initial-guess' functions g_n.
!  ______________________________________________
! |                                              |
! |   | w_n > = \sum_m | Psi_m > U_mn            |
! |                                              |
! |   U_mn = (V   DELTA^T   SVD_W^dagger)_mn     |
! |                                              |
! |   < Psi_m | g_n > = ( V SIGMA W^dagger)_mn   |
! |______________________________________________|
!
! In general (DELTA_ij) = Kronecker_delta_ij is rectangular.
!
! < Psi_m | g_n > is stored in LPRJ_COVL, which means that LPRJ_PROALL must
! have been called first.
! The end result is returned to the calling environment.
! U_mn is exactly what LINCOM uses as CTRANS.
!
! _test_mani
! I double checked this. Should be good!
! Only works for GAMMA at the moment.
!******************************************************************************

      SUBROUTINE LPRJ_CALC_WANNIER_TRAFO(W, LPRJ_COVL, UK, IO)
         ! USE new_tutor
         USE mkpoints
         USE wave
         USE base
         USE full_kpoints
         USE string, ONLY: str
         USE tutor, ONLY: vtutor
         IMPLICIT none

         TYPE(wavespin), INTENT(IN)  :: W
         !> Projection matrices (contents of the A_{mn}^k file)
         GDEF, INTENT(IN)            :: LPRJ_COVL(:,:,:,:)
         !> U_{nm}^k tranformation matrix between Bloch and Wannier functions
         !> (band, wannier idx, kpoint, spin)
         GDEF, ALLOCATABLE, INTENT(OUT) :: UK(:, :, :, :)
         TYPE(in_struct), INTENT(IN)    :: IO

         ! locals
         GDEF, ALLOCATABLE    :: A(:, :)
         INTEGER              :: IK, ISP
         INTEGER              :: NB_TOT, NKPTS, NSPIN, NWAN
         INTEGER              :: WAN_IDX, BAND_IDX, IB
         GDEF, ALLOCATABLE    :: SVD_V(:, :), SVD_WD(:, :)
         REAL(q), ALLOCATABLE :: SVD_SIG(:)
         REAL(q), ALLOCATABLE :: WEIGHTS(:, :, :)

         INTEGER              :: INFO, LWORK
         INTEGER              :: IKIBZ
         INTEGER              :: MN
         GDEF                 :: ASKWORK
         INTEGER, ALLOCATABLE :: IWORK(:)
         REAL(q), ALLOCATABLE :: RWORK(:)
         GDEF, ALLOCATABLE    :: WORK(:)

         CALL CHECK_FULL_KPOINTS

         NB_TOT = SIZE(LPRJ_COVL, 1)
         NKPTS = SIZE(LPRJ_COVL, 2)
         NSPIN = SIZE(LPRJ_COVL, 3)
         NWAN = SIZE(LPRJ_COVL, 4)

         IF (NB_TOT < NWAN) THEN
            CALL vtutor%bug("internal error in VASP: LPRJ_CALC_WANNIER_TRAFO: NB_TOT should be " &
               // "larger or equal to NWAN", __FILE__, __LINE__)
         ENDIF
         ALLOCATE(A(NB_TOT, NWAN))

         IF (ALLOCATED(UK)) DEALLOCATE(UK)
         ALLOCATE(UK(NB_TOT, NWAN, NKPTS, NSPIN))

         ALLOCATE(SVD_SIG(MIN(NB_TOT, NWAN)), SVD_V(NB_TOT, NB_TOT), SVD_WD(NWAN, NWAN))
         ALLOCATE(IWORK(8 * NWAN))

         ! If there are more bands than Wannier functions available,
         ! we apply a smooth cutoff using a Fermi function.
         CALL CALC_FERMI_SMEARING(W, NWAN, 1._q, WEIGHTS) ! _test_mani

         ! Loop over all spin and k components.
         spin: DO ISP = 1, NSPIN
         kpoints: DO IK = 1, NKPTS
            IKIBZ = KPOINTS_FULL%NEQUIV(IK)

            ! We require contiguous memory.
            A = LPRJ_COVL(:, IK, ISP, :)
            DO WAN_IDX = 1, NWAN
               A(:, WAN_IDX) = A(:, WAN_IDX) * WEIGHTS(:, IKIBZ, ISP)
            ENDDO

#ifdef gammareal
            ! Query optimal size of work arrays.
            CALL DGESDD(   'A', NB_TOT, NWAN, A, NB_TOT, SVD_SIG, &
                           SVD_V, NB_TOT, SVD_WD, NWAN, ASKWORK, -1, IWORK, INFO)
#else
            MN = MIN(NB_TOT,NWAN)
            ALLOCATE( RWORK(5*MN*MN + 5*MN) )
            CALL ZGESDD(   'A', NB_TOT, NWAN, A, NB_TOT, SVD_SIG, &
                           SVD_V, NB_TOT, SVD_WD, NWAN, ASKWORK, -1, RWORK, IWORK, INFO)
#endif

            IF (INFO /= 0) CALL vtutor%error("LPRJ_CALC_WANNIER_TRAFO: DGESDD query returned: " // str(INFO))

            LWORK = INT(ASKWORK)
            ALLOCATE(WORK(LWORK))

#ifdef gammareal
            ! Perform singular value decomposition of A
            CALL DGESDD(   'A', NB_TOT, NWAN, A, NB_TOT, SVD_SIG, &
                           SVD_V, NB_TOT, SVD_WD, NWAN, WORK, LWORK, IWORK, INFO)
#else
            CALL ZGESDD(   'A', NB_TOT, NWAN, A, NB_TOT, SVD_SIG, &
                           SVD_V, NB_TOT, SVD_WD, NWAN, WORK, LWORK, RWORK, IWORK, INFO)
            DEALLOCATE( RWORK )
#endif

            IF (INFO /= 0) CALL vtutor%error("LPRJ_CALC_WANNIER_TRAFO: DGESDD returned: " // str(INFO))

            ! Construct transformation matrix U_k
            UK(:, :, IK, ISP) = 0
            DO WAN_IDX = 1, NWAN
              DO IB = 1, MIN(NWAN, NB_TOT)
                DO BAND_IDX = 1, NB_TOT
                  UK(BAND_IDX, WAN_IDX, IK, ISP) = &
                    UK(BAND_IDX, WAN_IDX, IK, ISP) + &
                    SVD_V(BAND_IDX, IB) * SVD_WD(IB, WAN_IDX)
                ENDDO
              ENDDO
            ENDDO
            DEALLOCATE(WORK)

         ENDDO kpoints
         ENDDO spin

      ENDSUBROUTINE LPRJ_CALC_WANNIER_TRAFO





!***************** SUBROUTINE LPRJ_SVD_LOCALIZE *******************
!
! same as LPRJ_CALC_WANNIER_TRAFO but an arbitrary set of bands can 
! be localized (not necessarily occupieds). Furthermore assumes:
! #bands_to_localize = #trial_functions
!
!******************************************************************************

      SUBROUTINE LPRJ_SVD_LOCALIZE(UK, ISPIN, NBSTART, NBSTOP, IO)
         ! USE new_tutor
         USE mkpoints
         USE wave
         USE base
         IMPLICIT none

         GDEF,              INTENT(OUT) :: UK(:,:,:,:) ! (band, wannier idx, kpoint, spin)
         INTEGER,           INTENT(IN)  :: ISPIN
         INTEGER,           INTENT(IN)  :: NBSTART(2), NBSTOP(2) ! global band indices
         TYPE(in_struct),   INTENT(IN)  :: IO
         ! local
         GDEF, ALLOCATABLE    :: A(:, :)
         INTEGER              :: IK, ISP
         INTEGER              :: NB_TOT, NKPTS, NWAN
         INTEGER              :: WAN_IDX, BAND_IDX, IB
         GDEF, ALLOCATABLE    :: SVD_SIG(:), SVD_V(:, :), SVD_WD(:, :)
         INTEGER              :: INFO, LWORK
         INTEGER              :: NBANDS(2)
         GDEF                 :: ASKWORK
         INTEGER, ALLOCATABLE :: IWORK(:)
         GDEF, ALLOCATABLE    :: WORK(:)

#ifndef gammareal
           CALL vtutor%bug("oops: the SVD in LPRJ_SVD_LOCALIZE is not ready for complex orbitals",__FILE__,__LINE__)
#endif

         IF (.NOT. ALLOCATED(LPRJ_COVL)) THEN
            CALL vtutor%bug("LPRJ_COVL not allocated", __FILE__, __LINE__)
         ENDIF

         ! GAMMA ONLY and SPIN UNPOLARIZED
         IK = 1
         ISP = 1

         NB_TOT = SIZE(LPRJ_COVL, 1)
         NKPTS = SIZE(LPRJ_COVL, 2)
         NWAN = SIZE(LPRJ_COVL, 4)
         NBANDS = NBSTOP - NBSTART + 1

         IF((NBSTOP(1) .GT. NB_TOT) .OR. (NBSTOP(2) .GT. NB_TOT)) THEN
           CALL vtutor%error("LPRJ_SVD_LOCALIZE: NBSTOP must not be larger than NB_TOT")
           CALL vtutor%stopCode()
         ENDIF
         IF(NBANDS(1) .NE. NWAN) THEN
           CALL vtutor%error("LPRJ_SVD_LOCALIZE: #bands_to_localize is not equal to #trial_functions")
           CALL vtutor%stopCode()
         ENDIF
         IF((NBANDS(2) .GT. 1) .AND. (NBANDS(2) .NE. NWAN)) THEN
           CALL vtutor%error("LPRJ_SVD_LOCALIZE: #bands_to_localize is not equal to #trial_functions")
           CALL vtutor%stopCode()
         ENDIF

         ! Loop over all spin and k components.
         spin: DO ISP = 1, ISPIN
         kpoints: DO IK = 1, NKPTS

           ALLOCATE(A(NBANDS(ISP), NWAN))
           ALLOCATE(SVD_SIG(NWAN), SVD_V(NBANDS(ISP),NBANDS(ISP)), SVD_WD(NWAN, NWAN))
           ALLOCATE(IWORK(8*NWAN))

           ! We require contiguous memory.
           A(1:NBANDS(ISP), :) = LPRJ_COVL(NBSTART(ISP):NBSTOP(ISP), IK, ISP, :)

           ! Query optimal size of work arrays.
           CALL DGESDD('A', NBANDS(ISP), NWAN, A, NBANDS(ISP), SVD_SIG, &
                       SVD_V, NBANDS(ISP), SVD_WD, NWAN, ASKWORK, -1, IWORK, INFO)

           IF (INFO .NE. 0) THEN
             CALL vtutor%error("LPRJ_SVD_LOCALIZE: DGESDD query returned: " // str(INFO))
             CALL vtutor%stopCode()
           ENDIF

           LWORK = INT(ASKWORK)
           ALLOCATE(WORK(LWORK))

           ! Perform singular value decomposition of A
           CALL DGESDD('A', NBANDS(ISP), NWAN, A, NBANDS(ISP), SVD_SIG, &
                       SVD_V, NBANDS(ISP), SVD_WD, NWAN, WORK, LWORK, IWORK, INFO)

           IF (INFO .NE. 0) THEN
             CALL vtutor%error("LPRJ_SVD_LOCALIZE: DGESDD returned: " // str(INFO))
             CALL vtutor%stopCode()
           ENDIF
           DEALLOCATE(WORK)

           WRITE(*,*) " sum of singular values / NWAN) ", SUM(SVD_SIG)/NWAN

           ! Construct transformation matrix U_k
           UK(:,:,:,ISP) = zero
           DO WAN_IDX = 1, NWAN
             DO IB = 1, NWAN
               DO BAND_IDX = 1, NBANDS(ISP)
                 UK(BAND_IDX, WAN_IDX, IK, ISP) = &
                   UK(BAND_IDX, WAN_IDX, IK, ISP) + &
                   SVD_V(BAND_IDX, IB) * SVD_WD(IB, WAN_IDX)
               ENDDO
             ENDDO
           ENDDO

           DEALLOCATE(IWORK)
           DEALLOCATE(SVD_SIG, SVD_V, SVD_WD)
           DEALLOCATE(A)

         ENDDO kpoints
         ENDDO spin


      ENDSUBROUTINE LPRJ_SVD_LOCALIZE



!******************** SUBROUTINE LPRJ_CALC_WANNIER_TRAFO_QAO **********************
!
! Calculates the transformation matrix T(a,n) which can be used to transform 
! the Bloch orbitals |phi_a> into Wannier states |W_n> = sum_a T(a.n) |phi_a>
!
! QAO = quasi atomic orbitals
! 
! Notation:
!   NOCC = #occupied orbitals
!   NVIRT = #virtual orbitals
!   NB_TOT = NOCC + NVIRT
!   NWAN = #wannier functions = #trial functions to project on
!
! The transformation matrix T is obtained in two steps.
!
! (1) We ask for the optimal NWAN-NOCC dimensional subspace of the virtual space 
!     in order to aign the Wannier states with the target functions |alpha_n>
!     => M(a,n) = < phi_(NOCC+a) | alpha_n > , a=1...NVIRT , n=1...NWAN
!        M = U SIGMA VD (singular value decomposition)
!        optimal subspace is now spanned by states
!        | X_n > = sum_a^NVIRT sum_k^NWAN U(a,k) VD(k,n) | phi_a > , n=1...(NWAN-NOCC)
!
! (2) Define overlap matrix
! 
!                /
!               /   < phi_n | alpha_m >      , if 1 <= n <= NOCC
!     M(m,n) = I
!               \   < X_(n-NOCC) | alpha_m > , if NOCC < n <= NWAN
!                \
!    
!     M = U SIGMA VD (again singular value decompose)
!     and define transformation matrix T such that we find the Wannier states
! 
!                                                       /
!                                                      /  | phi_m > , 1 <= m <= NOCC 
!     | W_n > = sum_m^NWAN sum_k^NWAN U(m,k) VD(k,n)  |
!                                                      \  | X_m >   , NOCC+1 < m < NWAN
!                                                       \
!             = sum_a^NB_TOT T(a,n) | phi_a >
!     
!******************************************************************************

      SUBROUTINE LPRJ_CALC_WANNIER_TRAFO_QAO(W, LPRJ_COVL, T, IO)
        USE mkpoints
        USE wave
        USE base
        USE full_kpoints
        USE string, ONLY: str
        USE tutor, ONLY: vtutor

        IMPLICIT none
        TYPE(wavespin),    INTENT(IN)  :: W
        COMPLEX(q), ALLOCATABLE,  INTENT(IN) :: LPRJ_COVL(:,:,:,:) !< Projection matrices (contents of the A_{mn}^k file)
        !> U_{nm}^k tranformation matrix between Bloch and Wannier functions (band, wannier idx, kpoint, spin)
        GDEF, ALLOCATABLE, INTENT(OUT) :: T(:, :, :, :)
        TYPE(in_struct),   INTENT(IN)  :: IO
        ! local
        GDEF, ALLOCATABLE :: M(:, :), U(:,:), VD(:,:), T1(:,:), T2(:,:)
        REAL(q), ALLOCATABLE :: SIGMA(:)
        INTEGER :: NB_TOT, NKPTS, NSPIN, NWAN, NVIRT, NOCC, MN
        INTEGER :: ISP, IK, IKIBZ, I, A
        GDEF :: ASKWORK(1)
        INTEGER, ALLOCATABLE :: IWORK(:)
        GDEF, ALLOCATABLE :: WORK(:)
        REAL(q), ALLOCATABLE :: RWORK(:)
        INTEGER :: LWORK, INFO

        CALL CHECK_FULL_KPOINTS

        IF (.NOT. ALLOCATED(LPRJ_COVL)) THEN
          CALL vtutor%error("LPRJ_CALC_WANNIER_TRAFO_QAO: LPRJ_COVL is not allocated")
          CALL vtutor%stopCode()
        ENDIF

        NB_TOT = SIZE(LPRJ_COVL, 1)
        NKPTS  = SIZE(LPRJ_COVL, 2)
        NSPIN  = SIZE(LPRJ_COVL, 3)
        NWAN   = SIZE(LPRJ_COVL, 4)

        IF (NB_TOT < NWAN) THEN
         CALL vtutor%bug("internal error in VASP: LPRJ_CALC_WANNIER_TRAFO_QAO: NB_TOT should be " &
            // "larger than NWAN", __FILE__, __LINE__)
         CALL vtutor%stopCode()
        ENDIF

        IF (ALLOCATED(T)) DEALLOCATE(T)
        ALLOCATE(T(NB_TOT, NWAN, NKPTS, NSPIN))
        T = 0 ! nullify
        
        ALLOCATE( IWORK(8*NWAN) ) ! auxilary array for SVD

        ! Loop over all spin and k components.
        spin: &
        DO ISP = 1, NSPIN
          kpoints: &
          DO IK = 1, NKPTS
            IKIBZ = KPOINTS_FULL%NEQUIV(IK)

            ! get number of occupied and unoccupied (virtual) bands
            NOCC  = NINT(SUM(W%FERTOT(:, IKIBZ, ISP)))
            NVIRT = NB_TOT - NOCC

            ALLOCATE( M(NVIRT,NWAN) )

            !!!!!!!!!!! first we determine the optimal unoccupied space for localization !!!!!!!!!!

            ! fetch the overlaps between virtuals and local functions
            DO A = NOCC+1, NB_TOT
               M(A-NOCC, :) =  LPRJ_COVL(A, IK, ISP, :) 
            ENDDO

            ! SVD: query optimal size of work arrays
            ALLOCATE( U(NVIRT,NVIRT) , VD(NWAN,NWAN),  SIGMA(NWAN) )
#ifdef gammareal
            CALL DGESDD('A', NVIRT, NWAN, M, NVIRT, SIGMA, &
                        U, NVIRT, VD, NWAN, ASKWORK, -1, IWORK, INFO)
#else
            MN = MIN(NVIRT,NWAN)
            ALLOCATE( RWORK(5*MN*MN + 5*MN) )
            CALL ZGESDD('A', NVIRT, NWAN, M, NVIRT, SIGMA, &
                        U, NVIRT, VD, NWAN, ASKWORK, -1, RWORK, IWORK, INFO)
#endif

            IF (INFO /= 0) THEN
              CALL vtutor%error("ERROR: GGESDD query returned: " // str(INFO))
              CALL vtutor%stopCode()
            ENDIF

            LWORK = ASKWORK(1)
            ALLOCATE(WORK(LWORK))

            ! SVD: perform M = U SIGMA VD
#ifdef gammareal
            CALL DGESDD('A', NVIRT, NWAN, M, NVIRT, SIGMA, &
                        U, NVIRT, VD, NWAN, WORK, LWORK, IWORK, INFO)
#else
            CALL ZGESDD('A', NVIRT, NWAN, M, NVIRT, SIGMA, &
                        U, NVIRT, VD, NWAN, WORK, LWORK, RWORK, IWORK, INFO)
            DEALLOCATE( RWORK )
#endif

            ! construct T1(a,m) = sum_n^NWAN U(a,n) VD(n,m) where m is from [1,NWAN-NOCC]
            ALLOCATE( T1(NVIRT,NWAN-NOCC) )
            T1 = 0
            CALL GGEMM('N', 'N', NVIRT, NWAN-NOCC, NWAN, one, &
                       U, NVIRT, VD, NWAN, zero, T1, NVIRT)

            DEALLOCATE(WORK)
            DEALLOCATE(SIGMA)
            DEALLOCATE(VD)
            DEALLOCATE(U)
            DEALLOCATE(M)

            !!!!!!!!!! now we are able to construct the overlap matrix with equal rank as !!!!!!!!!!
            !!!!!!!!!! the localized starting orbitals (e.g atomic orbitals) = NWAN       !!!!!!!!!!

            ALLOCATE( M(NWAN,NWAN) )

            ! first, fetch the overlaps with occupied space
            M = 0
            DO I = 1, NOCC
              M(I, :) =  LPRJ_COVL(I, IK, ISP, :)
            ENDDO

            ! second, fetch the overlaps with the optimized unoccupied space 
            ! for this we missuse U and VD as temp/auxilary matrices
            ALLOCATE( U(NVIRT,NWAN), VD(NWAN-NOCC,NWAN) )
            U(:,:) = LPRJ_COVL(NOCC+1:NB_TOT, IK, ISP, :) ! store overlaps in U

            ! VD = T1 * U
            CALL GGEMM('C', 'N', NWAN-NOCC, NWAN, NVIRT, one, &
                       T1, NVIRT, U, NVIRT, zero, VD, NWAN-NOCC)

            M(NOCC+1:NWAN,:) = VD(:,:) ! copy back to correct index positions

            DEALLOCATE(U)
            DEALLOCATE(VD)

            ! SVD: query optimal size of work arrays
            ALLOCATE( U(NWAN,NWAN), VD(NWAN,NWAN), SIGMA(NWAN) )
#ifdef gammareal
            CALL DGESDD('A', NWAN, NWAN, M, NWAN, SIGMA, &
                        U, NWAN, VD, NWAN, ASKWORK, -1, IWORK, INFO)
#else
            ALLOCATE( RWORK(5*NWAN*NWAN + 5*NWAN) )
            CALL ZGESDD('A', NWAN, NWAN, M, NWAN, SIGMA, &
                        U, NWAN, VD, NWAN, ASKWORK, -1, RWORK, IWORK, INFO)
#endif
            IF (INFO /= 0) THEN
               CALL vtutor%error("ERROR: GGESDD query returned: " // str(INFO))
               CALL vtutor%stopCode()
            ENDIF

            LWORK = ASKWORK(1)
            ALLOCATE(WORK(LWORK))

            ! SVD: perform M = U SIGMA VD
#ifdef gammareal
            CALL DGESDD('A', NWAN, NWAN, M, NWAN, SIGMA, &
                        U, NWAN, VD, NWAN, WORK, LWORK, IWORK, INFO)
#else
            CALL ZGESDD('A', NWAN, NWAN, M, NWAN, SIGMA, &
                        U, NWAN, VD, NWAN, WORK, LWORK, RWORK, IWORK, INFO)
            DEALLOCATE( RWORK )
#endif
            IF (INFO /= 0) THEN
               CALL vtutor%error("ERROR: GGESDD query returned: " // str(INFO))
               CALL vtutor%stopCode()
            ENDIF

            ! construct T2 = U VD
            ALLOCATE( T2(NWAN,NWAN) )
            CALL GGEMM('N', 'N', NWAN, NWAN, NWAN, one, &
                       U, NWAN, VD, NWAN, zero, T2, NWAN)

            DEALLOCATE(WORK)
            DEALLOCATE(SIGMA)
            DEALLOCATE(VD)
            DEALLOCATE(U)
            DEALLOCATE(M)

            !!!!!!!!!! construct the final transformation matrix !!!!!!!!!!

            ! for the occupied indices, T is equal to T2
            T(1:NOCC,:, IK, ISP) = T2(1:NOCC,:)

            ! missuse U and VD as temp/auxilary matrices
            ALLOCATE( U(NWAN-NOCC,NWAN), VD(NVIRT,NWAN) )
            U(:,:) = T2(NOCC+1:NWAN,:)
            ! VD = T1 * U
            CALL GGEMM('N', 'N', NVIRT, NWAN, NWAN-NOCC, one, &
                       T1, NVIRT, U, NWAN-NOCC, zero, VD, NVIRT)
            
            T(NOCC+1:NB_TOT,:,IK,ISP) = VD(:,:) ! copy to correct index positions

            DEALLOCATE(U)
            DEALLOCATE(VD)

            DEALLOCATE(T2)
            DEALLOCATE(T1)
          ENDDO kpoints
        ENDDO spin

        DEALLOCATE(IWORK)
      END SUBROUTINE LPRJ_CALC_WANNIER_TRAFO_QAO

      SUBROUTINE CALC_FERMI_SMEARING(W, NB_EDGE, SIGMA, WEIGHTS)
         USE wave
         USE mkpoints
         IMPLICIT NONE

         TYPE(wavespin),               INTENT(IN)  :: W
         INTEGER,                      INTENT(IN)  :: NB_EDGE
         REAL(q),                      INTENT(IN)  :: SIGMA
         REAL(q), ALLOCATABLE, TARGET, INTENT(OUT) :: WEIGHTS(:, :, :)

         INTEGER, PARAMETER   :: NEDOS=256
         REAL(q)              :: ENTROPY, EFERMI, OCC
         REAL(q)              :: DOS(NEDOS,W%WDES%ISPIN),DOSI(NEDOS,W%WDES%ISPIN)
         REAL(q)              :: PAR(1,1,1,1,W%WDES%NCDIJ),DOSPAR(1,1,1,W%WDES%NCDIJ)
         TYPE (wavespin)      :: W_TMP
         
         IF (ALLOCATED(WEIGHTS)) DEALLOCATE(WEIGHTS)
         ALLOCATE(WEIGHTS(W%WDES%NB_TOT, W%WDES%NKPTS, W%WDES%ISPIN))

         ! Only perform smearing if the edge is inside the range of bands.
         IF (W%WDES%NB_TOT > NB_EDGE) THEN
            ! Some trickery is required.
            W_TMP=W
            W_TMP%FERTOT => WEIGHTS
            W_TMP%FERWE => W_TMP%FERTOT(W%WDES%NB_LOW : W%WDES%NB_TOT : W%WDES%NB_PAR, :, :)

            ! Calculate the partial occupancies using smearing.
            OCC = 2._q; IF (W%WDES%ISPIN==2.OR.W%WDES%NRSPINORS==2) OCC = 1._q
            CALL DENSTA(-1, -1, W%WDES, W_TMP, KPOINTS_, NB_EDGE * OCC, &
               0._q, ENTROPY, EFERMI, SIGMA, .FALSE., &
               NEDOS, 0, 0, DOS, DOSI, PAR, DOSPAR)
         ELSE
            WEIGHTS = 1
         ENDIF
      END SUBROUTINE CALC_FERMI_SMEARING

!******************** SUBROUTINE W1_FROM_LPRJ_function ***********************
!
!> @brief Calculate \f$ | RYlm_i >\f$ using CONSTRUCT_RYLM for the
!! current LPRJ_function as given by FG. Here, i runs
!! over a set of local functions at site R=POS with
!! l = [LSTART,LSTOP], and m = [1,2*l+1]
!! Note: this routine is essentially a copy of the first part of
!! CALC_OVERLAP_GN
!
!*****************************************************************************
    SUBROUTINE W1_FROM_LPRJ_function(LSTART, LSTOP, FG, POS, ISPINOR,       &
                                     SPIN_QAXIS, LROTYLM, PROJ_X, PROJ_Z, CQIJ, &
                                     LATT_CUR, NONL_S, NONLR_S, &
                                     LNORM, LAUG, HYBRID_ORBITAL, W1)
      USE poscar
      USE pseudo
      USE dfast
      USE lattice
      USE full_kpoints
      USE wave_high
      USE nonl_high
      USE base, ONLY : info_struct
      INTEGER         :: LSTART !< L index at which to start generating orbitals
      INTEGER         :: LSTOP !< L index at which to stop generating orbitals
      REAL(q)         :: FG(:,:,:) !< information about the radial part (number of points,coefficients,nylm)
      REAL(q)         :: POS(3) !< position in fractional coordinates where to center the orbital
      INTEGER         :: ISPINOR !< index of teh spinor for this orbital can be 1 (up) or 2 (down)
      REAL(q)         :: SPIN_QAXIS(3) !< spin quantization axis
      LOGICAL         :: LROTYLM !< whether a rotation of the angular part is neede or not
      REAL(q)         :: PROJ_X(3) !< rotate angular part along the X direction
      REAL(q)         :: PROJ_Z(3) !< rotate angular part along the Y direction
      OVERLAP         :: CQIJ(:,:,:,:)
      TYPE(nonl_struct)  :: NONL_S !< structure to compute the non-local part
      TYPE(nonlr_struct) :: NONLR_S !< structure to compute the non-local part (real space)
      TYPE(latt)      :: LATT_CUR
      LOGICAL         :: LNORM !< normalize the function
      LOGICAL         :: LAUG !< compute augmentation charges
      !> specify coefficients for linear combination of orbitals
      REAL(q)         :: HYBRID_ORBITAL((LMAX+1)**2)
      TYPE(wavefun1)  :: W1 !< output of this routine
      ! local
      GDEF, ALLOCATABLE      :: ROTYLM(:,:)
      GDEF, ALLOCATABLE      :: GHYBRID_ORBITAL(:)
      TYPE(wavedes1), POINTER:: WDESK
      TYPE(wavefuna)         :: WRYLM
      REAL(q)                :: WSCAL
      INTEGER                :: NYLM
      INTEGER                :: LMBASE

      PROFILING_START("w1_from_lprj_function")

      LMBASE = LSTART**2
      NYLM=(LSTOP+1)**2-LSTART**2

      WDESK=>W1%WDES1
      CALL NEWWAVA(WRYLM, WDESK, NYLM)

      ! Build R*YLM functions |RYLM>
      CALL CONSTRUCT_RYLM(LSTART, LSTOP, FG, LATT_CUR, POS, NONLR_S, NONL_S, &
                          ISPINOR, SPIN_QAXIS, WRYLM, LAUG)

      ! if needed rotate the ylm
      IF (LROTYLM) THEN
          ALLOCATE(ROTYLM((LMAX+1)**2,(LMAX+1)**2))
          CALL SETROTYLM(PROJ_X,PROJ_Z,LMAX,ROTYLM)

          ! |RYLM> = sum_LM' ROTYLM(LM,LM') |RYLM'>
          CALL LINCOM('F', &
                      WRYLM%CW(:,:), &
                      WRYLM%CPROJ(:,:), &
                      ROTYLM(LMBASE+1:LMBASE+NYLM,LMBASE+1:LMBASE+NYLM), &
                      NYLM, &
                      NYLM, &
                      WDESK%NPL, &
                      WDESK%NPRO, &
                      WDESK%NRPLWV, &
                      WDESK%NPROD, &
                      NYLM, &
                      WRYLM%CW(:,:), &
                      WRYLM%CPROJ(:,:))
          DEALLOCATE(ROTYLM)
      ENDIF

      ! build hybrid orbital that is a linear combination of |RYLM>
      ALLOCATE(GHYBRID_ORBITAL(NYLM))
#ifdef gammareal
      GHYBRID_ORBITAL = 0.0_q
      GHYBRID_ORBITAL(:) = HYBRID_ORBITAL(LMBASE+1:LMBASE+NYLM)
#else
      GHYBRID_ORBITAL = (0.0_q, 0.0_q)
      GHYBRID_ORBITAL(:) = &
            (1.0_q, 0.0_q) * HYBRID_ORBITAL(LMBASE+1:LMBASE+NYLM)
#endif

      ! the final result here is given by
      ! |RPHI> = sum_LM  GHYBRID_ORBITAL(LM') |RYLM'>
      CALL LINCOM('F', &
                  WRYLM%CW(:,:), &
                  WRYLM%CPROJ(:,:), &
                  GHYBRID_ORBITAL(:), &
                  NYLM, &
                  1, &
                  WDESK%NPL, &
                  WDESK%NPRO, &
                  WDESK%NRPLWV, &
                  WDESK%NPROD, &
                  NYLM, &
                  WRYLM%CW(:,:), &
                  WRYLM%CPROJ(:,:))

      ! save the first element
      ! because it contains the linear combination of the orbitals
      CALL W1_COPY(ELEMENT(WRYLM, 1), W1)

      ! normalize the final result
      IF (LNORM) THEN
        CALL CNORMN(W1, CQIJ, 1, WSCAL) ! ISP = 1
      ENDIF

      CALL DELWAVA(WRYLM)

      PROFILING_STOP("w1_from_lprj_function")

    END SUBROUTINE W1_FROM_LPRJ_function



!***************** SUBROUTINE LPRJ_CALC_INTRINSIC_ATOMIC_ORBITALS **************
!
! calculates the intrinsic atomic orbitals according to Knizia (DOI: 10.1021/ct400687b)
! actually we use here the simpler version as documented in Schaefer et.al. (DOI: 10.1063/5.0074936)
! 
! The orbitals defined by LOCPROJ in INCAR are assumed to define an atom
! centered minimal basis. (basis B2 in Knizias paper)
! These functions are stored in W_TRIALFUNC.
!
! these orbitals can then be used to define IBOs (intrinsic bonding orbitals)
!
!******************************************************************************
    SUBROUTINE LPRJ_CALC_INTRINSIC_ATOMIC_ORBITALS(W, W_TRIALFUNC, LPRJ_COVL, LPRJ_COVL_IAO, W_IAO, LATT_CUR, IO, INFO, CQIJ)
      USE base
      USE dfast
      USE tutor
      USE choleski
      USE lattice
      USE fileio
      USE wave_high,  ONLY: W1_GATHER_GLB_NOFFT
      ! debug
      USE mathtools,  ONLY: SV_DECOMPOSER
      !> Kohn-sham orbitals from which to build the IAOS
      TYPE(wavespin) :: W
      TYPE(wavespin), INTENT(IN) :: W_TRIALFUNC
      GDEF, ALLOCATABLE, INTENT(IN) :: LPRJ_COVL(:,:,:,:)
      !> Contains the projections of the Kohn-Sham orbitals onto IAOs
      !> The dimensions are (NBANDS,NK,ISP,NPROJ)
      GDEF, ALLOCATABLE, INTENT(OUT) :: LPRJ_COVL_IAO(:,:,:,:)
      !> Container for the intrinsic atomic orbitals
      TYPE(wavespin) :: W_IAO
      TYPE(latt)        :: LATT_CUR
      TYPE (in_struct)  :: IO
      TYPE(info_struct) :: INFO
      OVERLAP           :: CQIJ(:,:,:,:)
      ! local
      INTEGER                     :: NB2    ! size of basis B2
      INTEGER                     :: I,J, I_GLOBAL, J_GLOBAL
      INTEGER                     :: MU, NU, MU_GLOBAL, NU_GLOBAL
      INTEGER                     :: IK, ISP, ISP_MU, MAX_NTARGET
      INTEGER                     :: NTARGET(2)
      INTEGER                     :: NTARGET_LOCAL_MPIROUNDED(2)
      INTEGER                     :: NTARGET_TOT_MPIROUNDED(2)
      INTEGER                     :: NBLOCH_HIGH_MPIROUNDED(2)
      INTEGER                     :: NBLOCH_LOW(2), NBLOCH_HIGH(2)
      INTEGER                     :: NB_COPY_LOCAL_LOW(2), NB_COPY_LOCAL_HIGH(2)
      GDEF, ALLOCATABLE           :: A(:,:)     ! overlap <mu|nu> of non-orthog. basis B2
      GDEF, ALLOCATABLE           :: S(:,:)     ! inverse overlap of non-orthog. basis B2
      GDEF, ALLOCATABLE           :: OVL_MUI(:,:,:)  ! overlap <mu|i>
      GDEF, ALLOCATABLE           :: OVL_JTMU(:,:,:) ! overlap <j^tilde|mu>
      GDEF, ALLOCATABLE           :: OVL_IAO(:,:,:) ! overlap <IAO|IAO>
      GDEF, ALLOCATABLE           :: GWORK(:,:)  
      GDEF, ALLOCATABLE           :: GWORK1(:)  
      REAL(q), ALLOCATABLE        :: RWORK(:)  
      COMPLEX(q), ALLOCATABLE     :: CW_RED(:,:)        ! tmp wave function array
      COMPLEX(q), ALLOCATABLE     :: CW_RED_SUM(:,:)    ! tmp wave function array
      GDEF, ALLOCATABLE           :: CPROJ_RED(:,:)     ! tmp wave function array
      GDEF, ALLOCATABLE           :: CPROJ_RED_SUM(:,:) ! tmp wave function array
      TYPE(wavefuna)              :: WA
      TYPE(wavefuna)              :: WA_TARGET
      TYPE(wavefuna)              :: WA_TILDE
      TYPE(wavefuna)              :: WA_IAO
      TYPE(wavefuna)              :: WA_TRIALFUNC
      TYPE(wavefun1), ALLOCATABLE :: W1GATHER(:) ! for blocking
      TYPE(wavespin)              :: W_TARGET       
      TYPE(wavespin)              :: W_TILDE
      TYPE(wavedes),POINTER       :: WDES_MU ! wavedes for trial functions and similar functions (IAOs,...)
      TYPE(wavedes)               :: WDESTARGET
      TYPE(wavedes),POINTER       :: WDES_IAO
      TYPE(wavedes1)              :: WDESK
      INTEGER                     :: IERR
      REAL(q)                     :: WSCAL, TRACE
      TYPE(SV_DECOMPOSER)         :: SVD

#ifndef gammareal
      CALL vtutor%error('Calculation of intrinsic atomic orbitals is only implemented for the vasp_gam version.')
#endif

      PROFILING_START('lprj_calc_intrinsic_atomic_orbitals')
        
      ! status output
      IF(IO%IU0 .GT. 0) THEN
        WRITE(IO%IU0,*) " constructing IAOs"
      ENDIF

      NB2 = SIZE(LPRJ_functions)
      ALLOCATE(A(NB2,NB2))
      ALLOCATE(S(NB2,NB2))
      A=zero
      S=zero

      ! no spin polarization for trial functions (always 1)
      ISP_MU = 1

      ! Gamma only 
      IK = 1

      ! gather all function in W_TRIALFUNC to this node
      ! and store in W1GATHER
      WDES_MU => W_TRIALFUNC%WDES
      CALL SETWDES(WDES_MU, WDESK, 0)
      ALLOCATE(W1GATHER(NB2))
      DO J = 1, NB2
        CALL NEWWAV(W1GATHER(J), WDESK, .FALSE.)
      ENDDO
      CALL SETWDES(WDES_MU, WDESK, IK)
      CALL W1_GATHER_GLB_NOFFT(W_TRIALFUNC, 1, NB2, ISP_MU, W1GATHER)

      CALL SETWDES(WDES_MU, WDESK, IK)
      WA_TRIALFUNC = ELEMENTS(W_TRIALFUNC, WDESK, ISP_MU) 

      ! status output
      IF(IO%IU0 .GT. 0) THEN
        WRITE(IO%IU0,*) "   calculating S=<mu|nu> and S^-1"
      ENDIF

      ! calc inverse of overlap matrix <mu|nu> of non-orthogonal basis B2
      ! first construct A=<mu|nu> (symmetric/hermitian)
      DO J_GLOBAL = 1, NB2
        DO I = 1, WDES_MU%NBANDS ! only local I to parallelize work
          I_GLOBAL = (I-1) * WDES_MU%NB_PAR + WDES_MU%NB_LOW
          IF(I_GLOBAL .GT. NB2) CYCLE
          A(I_GLOBAL,J_GLOBAL) = W1_DOT(ELEMENT(WA_TRIALFUNC,I), &
                                        W1GATHER(J_GLOBAL), & ! W1GATHER = TRIALFUNC
                                        CQIJ)
        ENDDO
      ENDDO
      ! gather results from all cores
      CALLMPI(M_sum_g(WDES_MU%COMM_INTER, A(1,1), SIZE(A)))

      ! invert overlap matrix A, S=A^-1 (using SVD and pseudo-inverse)
      CALL SVD%INIT('A', NB2, NB2, IERR) 
      CALL SVD%CALC(A, IERR) 
      S = GCONJG(TRANSPOSE(SVD%V_DAGGER_MAT))
      DO MU = 1, NB2
        IF(SVD%SINGULAR_VALS(MU) .GT. 1E-10_q) THEN
          S(:,MU) = S(:,MU) / SVD%SINGULAR_VALS(MU)
        ELSE
          S(:,MU) = zero
          CALL vtutor%alert("IAOs: overlap matrix of local trial functions not full rank (SVD found singular value " // &
               str(SVD%SINGULAR_VALS(MU)) // ") => only pseudo-inverse possible. This probably leads to deficient IAOs."  )
        ENDIF
      ENDDO
      S = MATMUL(S,GCONJG(TRANSPOSE(SVD%U_MAT)))
      CALL SVD%DEALLOC()

      ! set the target to the occupied orbitals
      NBLOCH_LOW(:) = 1
      NBLOCH_HIGH(:) = 0
      DO ISP = 1, W%WDES%ISPIN
        NBLOCH_HIGH(ISP) = NINT(SUM(W%FERTOT(:, IK, ISP)))
      ENDDO
      ! overwrite with INCAR flags
      IF(LOCPROJ_IAO_NBLOCH_LOW(1) .GT. 0)  NBLOCH_LOW(1)  = LOCPROJ_IAO_NBLOCH_LOW(1)
      IF(LOCPROJ_IAO_NBLOCH_LOW(2) .GT. 0)  NBLOCH_LOW(2)  = LOCPROJ_IAO_NBLOCH_LOW(2)
      IF(LOCPROJ_IAO_NBLOCH_HIGH(1) .GT. 0) NBLOCH_HIGH(1) = LOCPROJ_IAO_NBLOCH_HIGH(1)
      IF(LOCPROJ_IAO_NBLOCH_HIGH(2) .GT. 0) NBLOCH_HIGH(2) = LOCPROJ_IAO_NBLOCH_HIGH(2)

      NTARGET(:) = NBLOCH_HIGH(:) - NBLOCH_LOW(:) + 1
      MAX_NTARGET = MAX(NTARGET(1), NTARGET(2))
      NTARGET_LOCAL_MPIROUNDED(:) = (NTARGET(:) + W%WDES%NB_PAR - 1) / W%WDES%NB_PAR
      NTARGET_TOT_MPIROUNDED(:) = NTARGET_LOCAL_MPIROUNDED(:) * W%WDES%NB_PAR
      NBLOCH_HIGH_MPIROUNDED(:) = NBLOCH_LOW(:) + NTARGET_TOT_MPIROUNDED(:) - 1


      !DO ISP = 1, W%WDES%ISPIN
      !  ALLOCATE(GWORK(NTARGET(ISP),NTARGET(ISP)))
      !  GWORK(:,:) = MATMUL(                                                                       &
      !                 MATMUL(                                                                     &
      !                   LPRJ_COVL(NBLOCH_LOW(ISP):NBLOCH_HIGH(ISP), IK, ISP, :),                  &
      !                   S(:,:)                                                                    &
      !                 ),                                                                          &
      !                 TRANSPOSE(GCONJG(LPRJ_COVL(NBLOCH_LOW(ISP):NBLOCH_HIGH(ISP), IK, ISP, :)))  &
      !               )
      !  TRACE = 0.0_q
      !  DO I = 1, NTARGET(ISP)
      !    TRACE = TRACE + ABS(GWORK(I,I))
      !  ENDDO
      !  IF(IO%IU0 .GT. 0) WRITE(IO%IU0,'(A,E14.6)') "trace of proj: ", TRACE
      !  IF(IO%IU6 .GT. 0) WRITE(IO%IU6,'(A,E14.6)') "trace of proj: ", TRACE
      !  !STOP
      !  TRACE = (1.0_q - TRACE / (1.0_q*NTARGET(ISP))) * 100.0_q
      !  !IF(IO%IU0 .GT. 0) WRITE(IO%IU0,'(A,F9.3)') "    Density spillage in %: ", TRACE
      !  !IF(IO%IU6 .GT. 0) WRITE(IO%IU6,'(A,F9.3)') "    Density spillage in %: ", TRACE
      !
      !  DEALLOCATE(GWORK)
      !ENDDO


      ! find lowest local index of considered bloch bands to copy into WTARGET
      DO ISP = 1, W%WDES%ISPIN
        DO I = 1, W%WDES%NBANDS
          I_GLOBAL = (I-1) * W%WDES%NB_PAR + W%WDES%NB_LOW
          IF(I_GLOBAL .GE. NBLOCH_LOW(ISP)) THEN
            NB_COPY_LOCAL_LOW(ISP) = I
            EXIT
          ENDIF
        ENDDO
      ENDDO
      NB_COPY_LOCAL_HIGH(:) = NB_COPY_LOCAL_LOW(:) + NTARGET_LOCAL_MPIROUNDED(:) - 1

      ! round WDESTARGET%NBANDS to next larger value dividable by WDESTARGET%NB_PAR
      WDESTARGET = W%WDES
      WDESTARGET%NB_TOT = MAX(NTARGET_TOT_MPIROUNDED(1), NTARGET_TOT_MPIROUNDED(2))
      WDESTARGET%NBANDS = MAX(NTARGET_LOCAL_MPIROUNDED(1), NTARGET_LOCAL_MPIROUNDED(2))
      WDESTARGET%NB_TOTK(IK,1) = NTARGET(1)
      IF(W%WDES%ISPIN .EQ. 2) WDESTARGET%NB_TOTK(IK,2) = NTARGET(2)
      CALL ALLOCW(WDESTARGET, W_TARGET) 
      ! copy original bloch orbitals into this "target" wave function array
      DO ISP = 1, W%WDES%ISPIN
        W_TARGET%CW(:,1:NTARGET_LOCAL_MPIROUNDED(ISP),:,ISP)    = W%CW(:,NB_COPY_LOCAL_LOW(ISP):NB_COPY_LOCAL_HIGH(ISP),:,ISP)
        W_TARGET%CPROJ(:,1:NTARGET_LOCAL_MPIROUNDED(ISP),:,ISP) = W%CPROJ(:,NB_COPY_LOCAL_LOW(ISP):NB_COPY_LOCAL_HIGH(ISP),:,ISP)
        W_TARGET%CELTOT(1:NTARGET_TOT_MPIROUNDED(ISP),:,ISP)  = W%CELTOT(NBLOCH_LOW(ISP):NBLOCH_HIGH_MPIROUNDED(ISP),:,ISP)
        W_TARGET%FERTOT(1:NTARGET_TOT_MPIROUNDED(ISP),:,ISP)  = W%FERTOT(NBLOCH_LOW(ISP):NBLOCH_HIGH_MPIROUNDED(ISP),:,ISP)
        W_TARGET%AUXTOT(1:NTARGET_TOT_MPIROUNDED(ISP),:,ISP)  = W%AUXTOT(NBLOCH_LOW(ISP):NBLOCH_HIGH_MPIROUNDED(ISP),:,ISP)
      ENDDO

      ! status output
      IF(IO%IU0 .GT. 0) THEN
        WRITE(IO%IU0,*) "   constructing |i tilde> "
      ENDIF

      CALL SETWDES(WDESTARGET, WDESK, IK)

      ! create wavespin to store |i^tilde>
      CALL ALLOCW(WDESTARGET, W_TILDE)

      ALLOCATE(OVL_MUI(NB2,MAX_NTARGET,WDESTARGET%ISPIN))
      OVL_MUI=zero

      CALL SETWDES(WDES_MU, WDESK, IK)
      ALLOCATE(CW_RED(WDESK%NRPLWV_RED, MAX_NTARGET))
      ALLOCATE(CPROJ_RED(WDESK%NPROD_RED, MAX_NTARGET))

      ALLOCATE(GWORK(NB2,MAX_NTARGET))

      DO ISP = 1, WDESTARGET%ISPIN
        
        CW_RED = (0.0_q, 0.0_q)
        CPROJ_RED = zero

        ! calc overlap <mu|i> 
        CALL SETWDES(WDESTARGET, WDESK, IK)
        WA_TARGET = ELEMENTS(W_TARGET, WDESK, ISP)
        DO I = 1, WDESTARGET%NBANDS
          I_GLOBAL = (I-1) * WDESTARGET%NB_PAR + WDESTARGET%NB_LOW
          IF(I_GLOBAL .GT. NTARGET(ISP)) CYCLE
          DO MU = 1, NB2
            OVL_MUI(MU,I_GLOBAL,ISP) = W1_DOT(W1GATHER(MU), &   ! W1GATHER = TRIALFUNC
                                              ELEMENT(WA_TARGET,I), &
                                              CQIJ)
          ENDDO
        ENDDO
        ! gather results from all cores
        CALLMPI(M_sum_g(WDESTARGET%COMM_INTER, OVL_MUI(:,:,ISP), SIZE(OVL_MUI(:,:,ISP))))

        ! check if trial functions are a good choice (have overlap with targets)
        DO I = 1, WDESTARGET%NBANDS
          I_GLOBAL = (I-1) * WDESTARGET%NB_PAR + WDESTARGET%NB_LOW
          IF(I_GLOBAL .GT. NTARGET(ISP)) CYCLE
          TRACE = 0.0_q
          DO MU = 1, NB2
            TRACE = TRACE + ABS(OVL_MUI(MU,I_GLOBAL,ISP))**2.0_q
          ENDDO
          IF(TRACE .LT. 0.1_q) THEN
            CALL vtutor%alert("IAO: The target band i=" // str(NBLOCH_LOW(ISP)+I_GLOBAL-1) // " is orthogonal &
                             &to your chosen trial functions: sum_mu <i|mu><mu|i> = " //str(TRACE) // ". &
                             &The construction of the IAOs will probably fail.")
          ENDIF
        ENDDO

        GWORK = zero
        GWORK(:,1:NTARGET(ISP)) = MATMUL(S, OVL_MUI(:,1:NTARGET(ISP),ISP)) 

        CALL SETWDES(WDES_MU, WDESK, IK)
        WA_TRIALFUNC = ELEMENTS(W_TRIALFUNC, WDESK, ISP_MU)

        ! build |i^tilde> in tmp arrays CW_RED & CPROJ_RED
        IF (WA_TRIALFUNC%WDES1%DO_REDIS) THEN
           CALL REDISTRIBUTE_PROJ(WA_TRIALFUNC)
           CALL REDISTRIBUTE_PW(WA_TRIALFUNC)
        ENDIF

        CALL GGEMM('N',                      & ! PW coefficients
                   'N',                      &
                   m_ WDESK%NPL_RED,         &
                   MAX_NTARGET,              &
                   NB2,                      &
                   one,                      &
                   WA_TRIALFUNC%CW_RED(1,1), &
                   m_ WDESK%NRPLWV_RED,      &
                   GWORK(1,1),               &
                   NB2,                      &
                   zero,                     &
                   CW_RED(1,1),              &
                   m_ WDESK%NRPLWV_RED)
        CALL GGEMM('N',                         & ! projectors
                   'N',                         &
                   WDESK%NPRO_RED,              &
                   MAX_NTARGET,                 &
                   NB2,                         &
                   one,                         &
                   WA_TRIALFUNC%CPROJ_RED(1,1), &
                   WDESK%NPROD_RED,             &
                   GWORK(1,1),                  &
                   NB2,                         &
                   zero,                        &
                   CPROJ_RED(1,1),              &
                   WDESK%NPROD_RED)
        IF (WA_TRIALFUNC%WDES1%DO_REDIS) THEN
           CALL REDISTRIBUTE_PROJ(WA_TRIALFUNC)
           CALL REDISTRIBUTE_PW(WA_TRIALFUNC)
        ENDIF

        CALL SETWDES(W_TILDE%WDES, WDESK, IK)
        WA_TILDE = ELEMENTS(W_TILDE, WDESK, ISP)
        IF (WA_TILDE%WDES1%DO_REDIS) THEN
           CALL REDISTRIBUTE_PROJ(WA_TILDE)
           CALL REDISTRIBUTE_PW(WA_TILDE)
        ENDIF
        WA_TILDE%CW_RED(:, 1:MAX_NTARGET)    = CW_RED(:, 1:MAX_NTARGET)
        WA_TILDE%CPROJ_RED(:, 1:MAX_NTARGET) = CPROJ_RED(:, 1:MAX_NTARGET)
        IF (WA_TILDE%WDES1%DO_REDIS) THEN
           CALL REDISTRIBUTE_PROJ(WA_TILDE)
           CALL REDISTRIBUTE_PW(WA_TILDE)
        ENDIF

      ENDDO ! ISP

      DEALLOCATE(GWORK)
      DEALLOCATE(CW_RED)
      DEALLOCATE(CPROJ_RED)

      ! symmetrically orthogonalize W_TILDE
      CALL LPRJ_LOEWDIN_ORTHONORMAL(W_TILDE, WDESTARGET, W_TILDE%WDES%ISPIN, [1,1], NTARGET(:), CQIJ)

      !CALL OUTWAV(IO, W_TILDE%WDES, W_TILDE, LATT_CUR, 0.0_q, 'tld')

      !
      ! construct IAOs
      ! 
      ! now that we have <G|mu>, <G|i^tilde>, and <G|i> 
      ! we can finally build  <G|IAO_mu> via
      !
      !  <G|IAO_mu> = <G|f1_mu> - <G|f2_mu> + <G|mu>
      !  where
      !  <G|f1_mu> = sum_i <G|i>  <i|mu>
      !  <G|f2_mu> = sum_i <G|i^tilde> <i^tilde|mu>
      !
      !  then renormalize (and orthogonalize) <G|IAO_mu>
      !

      ! status output
      IF(IO%IU0 .GT. 0) THEN
        WRITE(IO%IU0,*) "   calculating coefficients of IAOs "
      ENDIF

      ! create wavespin to store IAOs 
      ALLOCATE(WDES_IAO)
      WDES_IAO = WDES_MU
      WDES_IAO%ISPIN = WDESTARGET%ISPIN
      CALL ALLOCW(WDES_IAO, W_IAO)

      ! calc <i^tilde|mu>
      ALLOCATE(OVL_JTMU(MAX_NTARGET,NB2,WDESTARGET%ISPIN))
      OVL_JTMU = zero
      DO ISP = 1, WDESTARGET%ISPIN
        CALL SETWDES(W_TILDE%WDES, WDESK, IK)
        WA_TILDE = ELEMENTS(W_TILDE, WDESK, ISP)
        DO MU = 1, NB2
          DO I = 1, WDESTARGET%NBANDS
            I_GLOBAL = (I-1) * WDESTARGET%NB_PAR + WDESTARGET%NB_LOW
            IF(I_GLOBAL .GT. NTARGET(ISP)) CYCLE
            OVL_JTMU(I_GLOBAL,MU,ISP) = W1_DOT(ELEMENT(WA_TILDE,I), &
                                               W1GATHER(MU), &  ! W1GATHER = TRIALFUNC
                                               CQIJ)
          ENDDO
        ENDDO
      ENDDO
      ! gather results from all cores
      CALLMPI(M_sum_g(WDESTARGET%COMM_INTER, OVL_JTMU(1,1,1), SIZE(OVL_JTMU)))

      CALL SETWDES(WDES_MU, WDESK, IK)
      ALLOCATE(CW_RED(WDESK%NRPLWV_RED, NB2))
      ALLOCATE(CW_RED_SUM(WDESK%NRPLWV_RED, NB2))
      ALLOCATE(CPROJ_RED(WDESK%NPROD_RED, NB2))
      ALLOCATE(CPROJ_RED_SUM(WDESK%NPROD_RED, NB2))

      ALLOCATE(GWORK(MAX_NTARGET,NB2))

      DO ISP = 1, WDESTARGET%ISPIN

        CW_RED=(0.0_q, 0.0_q)
        CW_RED_SUM=(0.0_q, 0.0_q)
        CPROJ_RED=zero
        CPROJ_RED_SUM=zero

        GWORK = zero

        ! build <G|f1_mu> and store in CW_RED_SUM and CPROJ_RED_SUM
        GWORK(1:NTARGET(ISP),:) = GCONJG(TRANSPOSE(OVL_MUI(:,1:NTARGET(ISP),ISP)))
        CALL SETWDES(WDESTARGET, WDESK, IK)
        WA_TARGET = ELEMENTS(W_TARGET, WDESK, ISP) 
        IF (WA_TARGET%WDES1%DO_REDIS) THEN
           CALL REDISTRIBUTE_PROJ(WA_TARGET)
           CALL REDISTRIBUTE_PW(WA_TARGET)
        ENDIF
        CALL GGEMM('N',                      & ! PW coefficients
                   'N',                      &
                   m_ WDESK%NPL_RED,         &
                   NB2,                      &
                   MAX_NTARGET,              &
                   one,                      &
                   WA_TARGET%CW_RED(1,1),       &
                   m_ WDESK%NRPLWV_RED,      &
                   GWORK(1,1),               &
                   MAX_NTARGET,              &
                   zero,                     &
                   CW_RED(1,1),              &
                   m_ WDESK%NRPLWV_RED)
        CALL GGEMM('N',                         & ! projectors
                   'N',                         &
                   WDESK%NPRO_RED,              &
                   NB2,                         &
                   MAX_NTARGET,                 &
                   one,                         &
                   WA_TARGET%CPROJ_RED(1,1),       &
                   WDESK%NPROD_RED,             &
                   GWORK(1,1),                  &
                   MAX_NTARGET,                 &
                   zero,                        &
                   CPROJ_RED(1,1),              &
                   WDESK%NPROD_RED)
        IF (WA_TARGET%WDES1%DO_REDIS) THEN
           CALL REDISTRIBUTE_PROJ(WA_TARGET)
           CALL REDISTRIBUTE_PW(WA_TARGET)
        ENDIF
        CW_RED_SUM = CW_RED
        CPROJ_RED_SUM = CPROJ_RED

        ! build <G|f2_mu> and add to CW_RED_SUM and CPROJ_RED_SUM
        CALL SETWDES(WDES_MU, WDESK, IK)
        WA_TILDE  = ELEMENTS(W_TILDE, WDESK, ISP)
        IF (WA_TILDE%WDES1%DO_REDIS) THEN
           CALL REDISTRIBUTE_PROJ(WA_TILDE)
           CALL REDISTRIBUTE_PW(WA_TILDE)
        ENDIF
        CALL GGEMM('N',                      & ! PW coefficients
                   'N',                      &
                   m_ WDESK%NPL_RED,         &
                   NB2,                      &
                   MAX_NTARGET,              &
                   one,                      &
                   WA_TILDE%CW_RED(1,1),     &
                   m_ WDESK%NRPLWV_RED,      &
                   OVL_JTMU(1,1,ISP),        &
                   MAX_NTARGET,              &
                   zero,                     &
                   CW_RED(1,1),              &
                   m_ WDESK%NRPLWV_RED)
        CALL GGEMM('N',                         & ! projectors
                   'N',                         &
                   WDESK%NPRO_RED,              &
                   NB2,                         &
                   MAX_NTARGET,                 &
                   one,                         &
                   WA_TILDE%CPROJ_RED(1,1),     &
                   WDESK%NPROD_RED,             &
                   OVL_JTMU(1,1,ISP),           &
                   MAX_NTARGET,                 &
                   zero,                        &
                   CPROJ_RED(1,1),              &
                   WDESK%NPROD_RED)
        IF (WA_TILDE%WDES1%DO_REDIS) THEN
           CALL REDISTRIBUTE_PROJ(WA_TILDE)
           CALL REDISTRIBUTE_PW(WA_TILDE)
        ENDIF
        CW_RED_SUM = CW_RED_SUM - CW_RED
        CPROJ_RED_SUM = CPROJ_RED_SUM - CPROJ_RED

        ! finally add <G|mu>
        CALL SETWDES(WDES_MU, WDESK, IK)
        WA_TRIALFUNC = ELEMENTS(W_TRIALFUNC, WDESK, ISP_MU)
        IF (WA_TRIALFUNC%WDES1%DO_REDIS) THEN
           CALL REDISTRIBUTE_PROJ(WA_TRIALFUNC)
           CALL REDISTRIBUTE_PW(WA_TRIALFUNC)
        ENDIF
        CW_RED_SUM = CW_RED_SUM + WA_TRIALFUNC%CW_RED
        CPROJ_RED_SUM = CPROJ_RED_SUM + WA_TRIALFUNC%CPROJ_RED
        IF (WA_TRIALFUNC%WDES1%DO_REDIS) THEN
           CALL REDISTRIBUTE_PROJ(WA_TRIALFUNC)
           CALL REDISTRIBUTE_PW(WA_TRIALFUNC)
        ENDIF

        CALL SETWDES(WDES_IAO, WDESK, IK)
        WA_IAO = ELEMENTS(W_IAO, WDESK, ISP)

        IF (WA_IAO%WDES1%DO_REDIS) THEN
           CALL REDISTRIBUTE_PROJ(WA_IAO)
           CALL REDISTRIBUTE_PW(WA_IAO)
        ENDIF
        WA_IAO%CW_RED = CW_RED_SUM
        WA_IAO%CPROJ_RED = CPROJ_RED_SUM
        IF (WA_IAO%WDES1%DO_REDIS) THEN
           CALL REDISTRIBUTE_PROJ(WA_IAO)
           CALL REDISTRIBUTE_PW(WA_IAO)
        ENDIF

        ! re-normalize 
        DO MU = 1, WDES_IAO%NBANDS
          MU_GLOBAL = (MU-1) * WDES_IAO%NB_PAR + WDES_IAO%NB_LOW
          IF(MU_GLOBAL .GT. NB2) CYCLE
          CALL CNORMN(ELEMENT(WA_IAO,MU), CQIJ, ISP, WSCAL)
        ENDDO
         
      ENDDO ! ISP

      DEALLOCATE(GWORK)

      ! symmetrically  orthogonalize IAOs (Loewdin)
      IF(LOCPROJ_ORTH_IAO) THEN
        CALL LPRJ_LOEWDIN_ORTHONORMAL(W_IAO, WDES_IAO, W_IAO%WDES%ISPIN, [1,1], [NB2,NB2], CQIJ, OVL_IAO)
      ENDIF

      !
      ! finished construction of IAOs
      ! 

      ! status output
      IF(IO%IU0 .GT. 0) THEN
        WRITE(IO%IU0,*) "   calculating overlaps <i|IAO> "
      ENDIF

      ! calc overlap <i|IAO> and store in LPRJ_COVL_IAO
      CALL SETWDES(WDES_IAO, WDESK, IK)
      ALLOCATE(LPRJ_COVL_IAO(W%WDES%NB_TOT,1,W%WDES%ISPIN,NB2))
      LPRJ_COVL_IAO = zero
      DO ISP = 1, W%WDES%ISPIN
        CALL W1_GATHER_GLB_NOFFT(W_IAO, IK, NB2, ISP, W1GATHER)
        CALL SETWDES(W%WDES, WDESK, IK)
        WA = ELEMENTS(W, WDESK, ISP) 
        DO MU = 1, NB2
          DO I = 1, W%WDES%NBANDS
            I_GLOBAL = (I-1) * W%WDES%NB_PAR + W%WDES%NB_LOW
            LPRJ_COVL_IAO(I_GLOBAL,IK,ISP,MU) = W1_DOT(ELEMENT(WA,I), &
                                                       W1GATHER(MU), & ! W1GATHER = IAOs
                                                       CQIJ)
          ENDDO
        ENDDO
      ENDDO
      ! gather results from all cores
      CALLMPI(M_sum_g(W%WDES%COMM_INTER, LPRJ_COVL_IAO(1,1,1,1), SIZE(LPRJ_COVL_IAO)))

      ! status output
      IF(IO%IU0 .GT. 0) THEN
        WRITE(IO%IU0,*) "   calculating overlaps <IAO|IAO> "
      ENDIF

      ! calc and dump overlap matrix of IAOs
      CALL SETWDES(WDES_IAO, WDESK, IK)
      ALLOCATE(OVL_IAO(NB2,NB2,W%WDES%ISPIN))
      OVL_IAO = zero
      DO ISP = 1, W%WDES%ISPIN
        CALL W1_GATHER_GLB_NOFFT(W_IAO, IK, NB2, ISP, W1GATHER)
        CALL SETWDES(W%WDES, WDESK, IK)
        WA_IAO = ELEMENTS(W_IAO, WDESK, ISP) 
        DO MU = 1, NB2
          DO NU = 1, W_IAO%WDES%NBANDS
            NU_GLOBAL = (NU-1) * W_IAO%WDES%NB_PAR + W_IAO%WDES%NB_LOW
            IF(NU_GLOBAL .GT. NB2) CYCLE
            OVL_IAO(NU_GLOBAL,MU,ISP) = W1_DOT(ELEMENT(WA_IAO,NU), &
                                               W1GATHER(MU), & ! W1GATHER = IAOs
                                               CQIJ)
          ENDDO
        ENDDO
      ENDDO
      ! gather results from all cores
      CALLMPI(M_sum_g(W%WDES%COMM_INTER, OVL_IAO(1,1,1), SIZE(OVL_IAO)))

#ifdef tschaef_verbose
      ! dump
      IF(IO%IU0 .GT. 0) THEN
        DO ISP = 1, W_IAO%WDES%ISPIN
          FUNIT=700
          OPEN(unit = FUNIT,file = "SIAO_ISP" // str(ISP) // ".dat", FORM='FORMATTED', access='stream', STATUS='REPLACE')
          WRITE(FUNIT,'(I7)') NB2
          DO MU = 1, NB2
            DO NU = 1, NB2
#ifdef gammareal
              WRITE(FUNIT,'(2I7,E26.16)')  MU, NU, OVL_IAO(MU,NU,ISP)
#else
              WRITE(FUNIT,'(2I7,2E26.16)') MU, NU, REAL(OVL_IAO(MU,NU,ISP), KIND=q), AIMAG(OVL_IAO(MU,NU,ISP))
#endif
            ENDDO
          ENDDO
          CLOSE(FUNIT)
        ENDDO
      ENDIF
#endif

      ! check if IAOS are a valid orthonormal basis for target orbitals 
      ! by calcualting sum_{i,j,mu} <i|mu><mu|j> while i,j is restricted to target
      DO ISP = 1, W%WDES%ISPIN
        ALLOCATE(GWORK1(NB2))
        GWORK1 = zero
        DO MU = 1, NB2
          GWORK1(MU) = SUM(LPRJ_COVL_IAO(NBLOCH_LOW(ISP):NBLOCH_HIGH(ISP),IK,ISP,MU))
          GWORK1(MU) = GWORK1(MU) * SUM(GCONJG(LPRJ_COVL_IAO(NBLOCH_LOW(ISP):NBLOCH_HIGH(ISP),IK,ISP,MU)))
        ENDDO
        TRACE = REAL(SUM(GWORK1(:)), KIND=q)
        !IF( ABS(TRACE-1.0_q*NTARGET(ISP)) .GT. 1E-10_q) THEN
        !  CALL vtutor%alert('The constructed IAOs do not form a valid orthonormal basis for your target orbitals. &
        !                    &If you set LOCPROJ_ORTH_IAO=.FALSE. this was to be expected. In any other &
        !                    &case, you probably chose inappropriate local trial functions.')
        !  IF(IO%IU0 .GT. 0) THEN
        !    WRITE(IO%IU0,'(A,F22.12,I10)') "  TRACE vs NTARGET ", TRACE, NTARGET(ISP)
        !  ENDIF
        !ENDIF
        !! DEBUG always dump
        !IF(IO%IU0 .GT. 0) THEN
        !  WRITE(IO%IU0,'(A,F22.12,I10)') "  TRACE vs NTARGET ", TRACE, NTARGET(ISP)
        !ENDIF
        DEALLOCATE(GWORK1)
      ENDDO

      ! status output
      IF(IO%IU0 .GT. 0) THEN
        WRITE(IO%IU0,*) " finished IAOs "
      ENDIF
       
      DEALLOCATE(CW_RED)
      DEALLOCATE(CW_RED_SUM)
      DEALLOCATE(CPROJ_RED)
      DEALLOCATE(CPROJ_RED_SUM)

      DEALLOCATE(OVL_IAO)
      DEALLOCATE(OVL_JTMU)

      CALL DEALLOCW(W_TARGET) 
      CALL DEALLOCW(W_TILDE)

      DO J = 1, NB2
        CALL DELWAV(W1GATHER(J), .FALSE.)
      ENDDO
      DEALLOCATE(W1GATHER)

      DEALLOCATE(OVL_MUI)
      DEALLOCATE(S)
      DEALLOCATE(A)

      PROFILING_STOP('lprj_calc_intrinsic_atomic_orbitals')

    END SUBROUTINE LPRJ_CALC_INTRINSIC_ATOMIC_ORBITALS



!**************** SUBROUTINE LPRJ_LOEWDIN_ORTHONORMAL *******************
!
! performs a symmetric Loewdin orthogonalization of the spin
! orbitals in W for bands NBSTART(ISP) to NBSTOP(ISP), ISP=1,ISPIN,
! then normalizes all orbitals
!
!***********************************************************************
    SUBROUTINE LPRJ_LOEWDIN_ORTHONORMAL(W, WDES, ISPIN, NBSTART, NBSTOP, CQIJ, SOUT)
      USE dfast 
      USE wave_high, ONLY: W1_GATHER_GLB_NOFFT
      USE mathtools, ONLY: SV_DECOMPOSER
      TYPE(wavespin) :: W
      TYPE(wavedes)  :: WDES
      INTEGER        :: ISPIN
      INTEGER        :: NBSTART(2)  ! for up/down
      INTEGER        :: NBSTOP(2)
      OVERLAP        :: CQIJ(:,:,:,:)
      GDEF, OPTIONAL :: SOUT(:,:,:)   ! overlap matrix as optional output
      ! local
      INTEGER                     :: IK, ISP
      INTEGER                     :: I, I_GLOBAL, J, J_GLOBAL
      INTEGER                     :: NBANDS(2)
      TYPE(wavedes1)              :: WDESK
      TYPE(wavefun1), ALLOCATABLE :: W1GATHER(:) ! for blocking
      TYPE(wavefuna)              :: WA
      GDEF, ALLOCATABLE           :: A(:,:)     ! overlap matrix
      GDEF, ALLOCATABLE           :: S(:,:)     ! inverse overlap 
      INTEGER                     :: IERR
      REAL(q)                     :: WSCAL
      TYPE(SV_DECOMPOSER)         :: SVD

      PROFILING_START('lprj_calc_loewdin_orthonorm')

      ! gamma-only so far
      IK = 1

      NBANDS(:) = NBSTOP(:) - NBSTART(:) + 1

      ! gather all Ws into W1GATHER
      DO ISP = 1, ISPIN
        CALL SETWDES(WDES, WDESK, 0)
        ALLOCATE(W1GATHER(NBANDS(ISP)))
        DO J = 1, NBANDS(ISP)
          CALL NEWWAV(W1GATHER(J), WDESK, .FALSE.)
        ENDDO
        CALL SETWDES(WDES, WDESK, IK)
        CALL W1_GATHER_GLB_NOFFT(W, NBSTART(ISP), NBSTOP(ISP), ISP, W1GATHER)
        WA = ELEMENTS(W, WDESK, ISP)

        ALLOCATE(A(NBANDS(ISP),NBANDS(ISP)))
        ALLOCATE(S(NBANDS(ISP),NBANDS(ISP)))
        A=zero
        S=zero
        DO I_GLOBAL = NBSTART(ISP), NBSTOP(ISP)
          DO J = 1, WDES%NBANDS 
            J_GLOBAL = (J-1) * WDES%NB_PAR + WDES%NB_LOW
            IF(J_GLOBAL .LT. NBSTART(ISP)) CYCLE
            IF(J_GLOBAL .GT. NBSTOP(ISP)) CYCLE
            A(I_GLOBAL,J_GLOBAL) = W1_DOT(W1GATHER(I_GLOBAL), & 
                                          ELEMENT(WA,J), & 
                                          CQIJ)
          ENDDO
        ENDDO
        ! gather results from all cores
        CALLMPI(M_sum_g(WDES%COMM_INTER, A(1,1), SIZE(A)))

        ! copy overlap matrix as output
        IF(PRESENT(SOUT)) SOUT(:,:,ISP) = A(:,:)

        ! invert overlap matrix A, S=A^-1 (using SVD and pseudo-inverse)
        CALL SVD%INIT('A', NBANDS(ISP), NBANDS(ISP), IERR) 
        CALL SVD%CALC(A, IERR) 
        S = GCONJG(TRANSPOSE(SVD%V_DAGGER_MAT))
        DO I = 1, NBANDS(ISP)
          IF(SVD%SINGULAR_VALS(I) .GT. 1E-10_q) THEN
            S(:,I) = S(:,I) / SQRT(SVD%SINGULAR_VALS(I))
          ELSE
            S(:,I) = zero
            CALL vtutor%alert("LPRJ_LOEWDIN_ORTHONORMAL: overlap matrix not full rank (SVD found singular value " // str(SVD%SINGULAR_VALS(I)) // ")!"  )
          ENDIF
        ENDDO
        S = MATMUL(S,GCONJG(TRANSPOSE(SVD%U_MAT)))
        CALL SVD%DEALLOC()

        ! apply S to symmetrically orthogonalize
        IF (WA%WDES1%DO_REDIS) THEN
           CALL REDISTRIBUTE_PROJ(WA)
           CALL REDISTRIBUTE_PW(WA)
        ENDIF
        CALL LINCOM('F', &
                    WA%CW_RED(:,NBSTART(ISP):NBSTOP(ISP)), &
                    WA%CPROJ_RED(:,NBSTART(ISP):NBSTOP(ISP)), &
                    S(:,:), &
                    NBANDS(ISP), &
                    NBANDS(ISP), &
                    WDESK%NPL_RED, &
                    WDESK%NPRO_RED, &
                    WDESK%NRPLWV_RED, &
                    WDESK%NPROD_RED, &
                    NBANDS(ISP), &
                    WA%CW_RED(:,NBSTART(ISP):NBSTOP(ISP)), &
                    WA%CPROJ_RED(:,NBSTART(ISP):NBSTOP(ISP)))
        IF (WA%WDES1%DO_REDIS) THEN
           CALL REDISTRIBUTE_PROJ(WA)
           CALL REDISTRIBUTE_PW(WA)
        ENDIF
  
        ! re-normalize
        DO I = 1, WDES%NBANDS
          I_GLOBAL = (I-1) * WDES%NB_PAR + WDES%NB_LOW
          IF(I_GLOBAL .LT. NBSTART(ISP)) CYCLE
          IF(I_GLOBAL .GT. NBSTOP(ISP)) CYCLE
          CALL CNORMN(ELEMENT(WA,I), CQIJ, ISP, WSCAL)
        ENDDO

        DEALLOCATE(A)
        DEALLOCATE(S)
        DO J = 1, NBANDS(ISP)
          CALL DELWAV(W1GATHER(J), .FALSE.)
        ENDDO
        DEALLOCATE(W1GATHER)

      ENDDO

      PROFILING_STOP('lprj_calc_loewdin_orthonorm')

    END SUBROUTINE LPRJ_LOEWDIN_ORTHONORMAL




!************************ SUBROUTINE LPRJ_CALC_SVD **************************
!
! Thin wrapper around DGESDD/ZGESDD for convenience to calculate a
! singular-value decomposition.
!
!***********************************************************************
      SUBROUTINE LPRJ_CALC_SVD(MODE, A, U, VD, SIGMA, IO)
         USE base
         USE string, ONLY  : str
         USE tutor,  ONLY  : vtutor

         CHARACTER(1),           INTENT(IN)     :: MODE
         GDEF,                   INTENT(INOUT)  :: A(:, :) !! Content is destroyed on exit
         GDEF, ALLOCATABLE,      INTENT(OUT)    :: U(:, :)
         GDEF, ALLOCATABLE,      INTENT(OUT)    :: VD(:, :)
         REAL(q), ALLOCATABLE,   INTENT(OUT)    :: SIGMA(:)
         TYPE(in_struct),        INTENT(IN)     :: IO

         INTEGER  :: DIM_1, DIM_2, MIN_DIM

         INTEGER              :: INFO, LWORK
         GDEF                 :: WORK_DIM
         INTEGER, ALLOCATABLE :: IWORK(:)
         GDEF, ALLOCATABLE    :: WORK(:)
         REAL(q), ALLOCATABLE :: RWORK(:)
         
         DIM_1 = SIZE(A, 1)
         DIM_2 = SIZE(A, 2)
         MIN_DIM = MIN(DIM_1, DIM_2)

         ! Make sure everything is deallocated first.
         IF (ALLOCATED(U)) DEALLOCATE(U)
         IF (ALLOCATED(VD)) DEALLOCATE(VD)
         IF (ALLOCATED(SIGMA)) DEALLOCATE(SIGMA)

         SELECT CASE (MODE)
            CASE('A')
               ALLOCATE(U(DIM_1, DIM_1), VD(DIM_2, DIM_2))
            CASE('S')
               ALLOCATE(U(DIM_1, MIN_DIM), VD(MIN_DIM, DIM_2))
            CASE DEFAULT
               CALL vtutor%bug("LPRJ_CALC_SVD: Mode '" // MODE // "' not implemented!", __FILE__, __LINE__)
         END SELECT

         ALLOCATE(SIGMA(MIN_DIM))
         ALLOCATE(IWORK(8 * MIN_DIM))

         ! Query optimal size of work arrays.
#ifdef gammareal
         CALL DGESDD(MODE, DIM_1, DIM_2, A, DIM_1, SIGMA, U, SIZE(U, 1), VD, SIZE(VD, 1), WORK_DIM, -1, IWORK, INFO)
#else
         ALLOCATE(RWORK(MIN_DIM * MAX(5 * (MIN_DIM + 1), 2 * (MAX(DIM_1, DIM_2) + MIN_DIM) + 1)))
         CALL ZGESDD(MODE, DIM_1, DIM_2, A, DIM_1, SIGMA, U, SIZE(U, 1), VD, SIZE(VD, 1), WORK_DIM, -1, RWORK, IWORK, INFO)
#endif
         IF (INFO /= 0) CALL vtutor%error("LPRJ_CALC_SVD: XGESDD query returned: " // str(INFO))

         LWORK = INT(WORK_DIM)
         ALLOCATE(WORK(LWORK))

         ! Perform singular value decomposition of A.
#ifdef gammareal
         CALL DGESDD(MODE, DIM_1, DIM_2, A, DIM_1, SIGMA, U, SIZE(U, 1), VD, SIZE(VD, 1), WORK, LWORK, IWORK, INFO)
#else
         !ALLOCATE(RWORK(LWORK))
         CALL ZGESDD(MODE, DIM_1, DIM_2, A, DIM_1, SIGMA, U, SIZE(U, 1), VD, SIZE(VD, 1), WORK, LWORK, RWORK, IWORK, INFO)
         DEALLOCATE(RWORK)
#endif
         IF (INFO /= 0) CALL vtutor%error("LPRJ_CALC_SVD: XGESDD returned: " // str(INFO))
      END SUBROUTINE LPRJ_CALC_SVD



!************************ SUBROUTINE LPRJ_ASSIGN_SITES_TO_ORBITALS **************************
!
! For each site a list of orbitals is generated based on the partial atomic charge.
! This allows to assign sites to localized orbitals for sitautions where the orbitals are not atom
! centered.
! Note that an orbital can appear in the list of multiple sites. In other words: 
! dont confuse "assign sites to orbitals" with "assign orbitals to sites".
!
!********************************************************************************************
      SUBROUTINE LPRJ_ASSIGN_ORBITALS_TO_SITES(LPRJ_COVL_IAO, IBLOW, IBHIGH, T_INFO, IO, UTRAFO_UP, UTRAFO_DOWN)
        USE poscar_struct_def
        USE base
        GDEF              :: LPRJ_COVL_IAO(:,:,:,:)
        INTEGER           :: IBLOW(2)  ! for spin up&down
        INTEGER           :: IBHIGH(2) ! for spin up&down
        TYPE(type_info)   :: T_INFO
        TYPE(in_struct)   :: IO
        GDEF, OPTIONAL    :: UTRAFO_UP(:,:)
        GDEF, OPTIONAL    :: UTRAFO_DOWN(:,:)
        ! local
        INTEGER              :: NAO_TOT, NB_TOT
        INTEGER              :: NATOMS, IATOM
        INTEGER              :: NB(2)
        INTEGER              :: I, J, MU
        INTEGER              :: IK, ISP, ISPIN
        GDEF, ALLOCATABLE    :: GOVL(:,:)
        REAL(q), ALLOCATABLE :: PBCHG(:,:) ! band:atom
        REAL(q)              :: THRESH, VAL
        INTEGER              :: FUNIT = 777

        ! Gamma-only at the moment
        IK = 1

        IF( (.NOT. LPRJ_IS_INITIALIZED()) .OR. (.NOT. LPRJ_IS_CALLED_PROALL()) ) THEN
          CALL vtutor%error("internal error in VASP: Module LOCPROJ is not &
                            &initialized but necessary for LPRJ_ASSIGN_ORBITALS_TO_SITES.")
          CALL vtutor%stopCode()
        ENDIF

        NB(:) = IBHIGH(:)-IBLOW(:)+1
        IF(PRESENT(UTRAFO_UP)) THEN
          IF(NB(1) .NE. SIZE(UTRAFO_UP(:,:),1)) THEN
            CALL vtutor%error("dimension missmatch in LPRJ_ASSIGN_ORBITALS_TO_SITES")
            CALL vtutor%stopCode()
          ENDIF
        ENDIF
        IF(PRESENT(UTRAFO_DOWN)) THEN
          IF(NB(2) .NE. SIZE(UTRAFO_DOWN(:,:),1)) THEN
            CALL vtutor%error("dimension missmatch in LPRJ_ASSIGN_ORBITALS_TO_SITES")
            CALL vtutor%stopCode()
          ENDIF
        ENDIF

        ! get number of atomic functions from LOCPROJ
        NAO_TOT = LPRJ_GET_NUM_WAN() 

        ! store total number of atoms (ions)
        NATOMS = T_INFO%NIONP

        ! get NB_TOT and ISPIN via LPRJ_COVL_*
        IF(LOCPROJ_IAO) THEN
          ISPIN = SIZE(LPRJ_COVL_IAO,3)
          NB_TOT = SIZE(LPRJ_COVL_IAO,1)
        ELSE 
          ISPIN = SIZE(LPRJ_COVL,3)
          NB_TOT = SIZE(LPRJ_COVL,1)
        ENDIF

        ! allocate projector
        ALLOCATE(GOVL(NB_TOT, NAO_TOT))
        GOVL = zero 

        DO ISP = 1, ISPIN

          ! get overlap matrix and apply trafo
          IF(LOCPROJ_IAO) THEN
            GOVL = LPRJ_COVL_IAO(:,IK,ISP,:)
          ELSE
            GOVL = LPRJ_COVL(:,IK,ISP,:)
          ENDIF
          IF(ISP .EQ. 1) THEN
            IF(PRESENT(UTRAFO_UP)) THEN
              GOVL(IBLOW(1):IBHIGH(1),:) = MATMUL(GCONJG(TRANSPOSE(UTRAFO_UP(:,:))), GOVL(IBLOW(1):IBHIGH(1),:))
            ENDIF
          ELSE
            IF(PRESENT(UTRAFO_DOWN)) THEN
              GOVL(IBLOW(2):IBHIGH(2),:) = MATMUL(GCONJG(TRANSPOSE(UTRAFO_DOWN(:,:))), GOVL(IBLOW(2):IBHIGH(2),:))
            ENDIF
          ENDIF

          ALLOCATE(PBCHG(NB(ISP),NATOMS))
          PBCHG = 0.0_q
          DO MU = 1, NAO_TOT
            IATOM = LPRJ_GET_SITE(MU)
            DO I = 1, NB(ISP)
              J = IBLOW(ISP) + I - 1
              PBCHG(I,IATOM) = PBCHG(I,IATOM) + ABS(GOVL(J,MU))**2.0_q
            ENDDO
          ENDDO

          ! find optimal threshold: the smallest of the "largest overlaps per local orbital"
          THRESH = 1.0_q
          DO I = 1, NB(ISP)
            VAL = MAXVAL(PBCHG(I,:))
            IF(VAL .LT. THRESH) THRESH = VAL
          ENDDO
          THRESH = THRESH * 0.9999_q ! make it a little smaller to avoid problems with IF >=

          IF(IO%IU6 .GT. 0) THEN
            OPEN(unit = FUNIT,file = "orbitals_of_sites_ISP" // str(ISP), FORM='FORMATTED', access='stream', STATUS='REPLACE')
            DO IATOM = 1, NATOMS
              DO I = 1, NB(ISP)
                IF(PBCHG(I,IATOM) .GT. THRESH) THEN
                  WRITE(FUNIT,'(I8)',advance='no') IBLOW(ISP) + I - 1
                ENDIF
              ENDDO
              WRITE(FUNIT,'(A)') " "
            ENDDO
            CLOSE(FUNIT)
          ENDIF

          DEALLOCATE(PBCHG)


        ENDDO

        DEALLOCATE(GOVL)

      END SUBROUTINE LPRJ_ASSIGN_ORBITALS_TO_SITES


!******************** SUBROUTINE LPRJ_LINCOM ***************************
! _test_mani more
! Apply the wannier rotation to the Bloch wave functions.
!***********************************************************************

      SUBROUTINE LPRJ_LINCOM(W, UK, W_WANNIER)
         USE dfast
         USE wave
         USE wave_high
         USE subrot
         USE tutor,  ONLY  : vtutor
         USE string, ONLY  : str
         IMPLICIT none

         TYPE(wavespin),   INTENT(IN)     :: W
         GDEF,             INTENT(IN)     :: UK(:, :, :, :)
         TYPE(wavespin),   INTENT(INOUT)  :: W_WANNIER

         ! LOCALS
         INTEGER IK, ISP
         INTEGER NWAN
         INTEGER NB_TOT

         TYPE(WAVEDES1) WDESK
         TYPE(WAVEFUNA) WA, WA_WAN

         NB_TOT = SIZE(UK, 1)
         NWAN = SIZE(UK, 2)

         ! This should never happen.
         IF (NB_TOT /= W%WDES%NB_TOT) CALL vtutor%bug("Internal error in LPRJ_LINCOM: " // &
            "Number of bands inconsistent: " // str(NB_TOT) // ", " // str(W%WDES%NB_TOT) // ".", &
            __FILE__, __LINE__)

         DO ISP = 1, W%WDES%ISPIN
            DO IK = 1, W%WDES%NKPTS
               CALL SETWDES(W%WDES, WDESK, IK)
               WA = ELEMENTS(W, WDESK, ISP)
               WA_WAN = ELEMENTS(W_WANNIER, WDESK, ISP)

               CALL WA_COPY(WA, WA_WAN)

               IF (W%WDES%DO_REDIS) THEN
                  CALL REDISTRIBUTE_PROJ(WA_WAN)
                  CALL REDISTRIBUTE_PW(WA_WAN)
               ENDIF

               CALL LINCOM('F', WA_WAN%CW_RED(:, :), WA_WAN%CPROJ_RED(:, :), &
                  UK(:, :, IK, ISP), NB_TOT, NWAN, WDESK%NPL_RED, WDESK%NPRO_RED, &
                  WDESK%NRPLWV_RED, WDESK%NPROD_RED, NB_TOT, &
                  WA_WAN%CW_RED(:, :), WA_WAN%CPROJ_RED(:, :))

               IF (W%WDES%DO_REDIS) THEN
                  CALL REDISTRIBUTE_PROJ(WA_WAN)
                  CALL REDISTRIBUTE_PW(WA_WAN)
               ENDIF

            ENDDO
         ENDDO
      END SUBROUTINE LPRJ_LINCOM


!******************** SUBROUTINE LPRJ_WRITE_FUNCTIONS ******************
!> Write LOCPROJ functions to OUTCAR
!***********************************************************************
      SUBROUTINE LPRJ_WRITE_FUNCTIONS(LPRJ_functions,LNONCOLLINEAR,IO)
         use base, ONLY: in_struct
         TYPE(LPRJ_function),INTENT(IN) :: LPRJ_functions(:)
         LOGICAL,INTENT(IN) :: LNONCOLLINEAR
         TYPE(in_struct),INTENT(IN) :: IO
         !local
         INTEGER :: IS

         IF (IO%IU6>=0) THEN
            WRITE(IO%IU6,*)
            WRITE(IO%IU6,*) 'LOCPROJ orbitals'
            IF (LNONCOLLINEAR) THEN
                WRITE(IO%IU6,*) '-----------------------------------------------------------------------------------------------------------------------------'
                WRITE(IO%IU6,*) ' n  l  m   za        pos                       proj_x                    proj_z                   spin  spin_qaxis'
                WRITE(IO%IU6,*) '-----------------------------------------------------------------------------------------------------------------------------'
                DO IS=1,SIZE(LPRJ_functions)
                   WRITE(IO%IU6,'(I3,I3,I3,1F8.3,A,3F8.3,A,3F8.3,A,3F8.3,A,I3,A,3F8.3)') LPRJ_functions(IS)%n, &
                                   LPRJ_functions(IS)%l, LPRJ_functions(IS)%m, &
                                   LPRJ_functions(IS)%za,'  ', LPRJ_functions(IS)%R,'  ',&
                                   LPRJ_functions(IS)%PROJ_X,'  ', LPRJ_functions(IS)%PROJ_Z,'  ',&
                                   LPRJ_functions(IS)%spinor, '  ',LPRJ_functions(IS)%spin_qaxis
                ENDDO
            ELSE
                WRITE(IO%IU6,*) '----------------------------------------------------------------------------------------------'
                WRITE(IO%IU6,*) ' n  l  m   za        pos                       proj_x                    proj_z'
                WRITE(IO%IU6,*) '----------------------------------------------------------------------------------------------'
                DO IS=1,SIZE(LPRJ_functions)
                   WRITE(IO%IU6,'(I3,I3,I3,1F8.3,A,3F8.3,A,3F8.3,A,3F8.3)') LPRJ_functions(IS)%n, &
                                   LPRJ_functions(IS)%l, LPRJ_functions(IS)%m, &
                                   LPRJ_functions(IS)%za,'  ', LPRJ_functions(IS)%R,'  ',&
                                   LPRJ_functions(IS)%PROJ_X,'  ', LPRJ_functions(IS)%PROJ_Z
                ENDDO
            ENDIF
            WRITE(IO%IU6,*)
         ENDIF

      END SUBROUTINE LPRJ_WRITE_FUNCTIONS

!******************** SUBROUTINE LPRJ_WRITE ****************************
!
!> Writes the projection of the Bloch orbitals onto the set of local
!! functions to the PROJCAR, LOCPROJ, and vasprun.xml files.
!!
!! The information written to these files is essentially the same,
!! but only PROJCAR is explicitly constructed to be human readable.
!!
!! The first line of the LOCPROJ file lists: NS, NK, NB, and NF
!! (# of spin channels, # of k-points, # of bands, and # of local functions)
!! followed by a block of NF lines describing the local functions.
!! After this block the projections are written:
!!
!! ( IFNC,LPRJ_COVL(IB,IK,ISP,IFNC), IFNC=1,NF, IB=1,NB, IK=1,NK, ISP=1,NS )
!!
!! i.e., IFNC is the fastest index, and ISP the slowest.
!
!***********************************************************************

      SUBROUTINE LPRJ_WRITE(LPRJ_functions,LPRJ_COVL,IU6,IU0,W)
      USE vaspxml
      USE wave
      USE full_kpoints
      TYPE (LPRJ_function), ALLOCATABLE :: LPRJ_functions(:)
      GDEF, ALLOCATABLE :: LPRJ_COVL(:,:,:,:)
      INTEGER IU6,IU0
      TYPE(wavespin) :: W
      ! local variables
      INTEGER NB,NK,NS,NF
      INTEGER IB,IK,IK_IRZ,ISP,IFNC,IBASE,ITMP(3)
      REAL(q), PARAMETER :: TINYR=1.E-4_q
      REAL(q) RTMP(2), DR
      LOGICAL, ALLOCATABLE :: LDONE(:), LDONOW(:)

     ! CHARACTER(LEN=10) :: LABEL(1:7,-5:3)= &
     !&(/"  sp3d2-1 ","  sp3d2-2 ","  sp3d2-3 ","  sp3d2-4 ","  sp3d2-5 ","  sp3d2-6 ","          ", &
     !&  "   sp3d-1 ","   sp3d-2 ","   sp3d-3 ","   sp3d-4 ","   sp3d-5 ","          ","          ", &
     !&  "   sp3-1  ","   sp3-2  ","   sp3-3  ","   sp3-4  ","          ","          ","          ", &
     !&  "   sp2-1  ","   sp2-2  ","   sp2-3  ","          ","          ","          ","          ", &
     !&  "    sp-1  ","    sp-2  ","          ","          ","          ","          ","          ", &
     !&  "      s   ","          ","          ","          ","          ","          ","          ", &
     !&  "     py   ","     pz   ","     px   ","          ","          ","          ","          ", &
     !&  "    dxy   ","    dyz   ","    dz2   ","    dxz   ","   dx2-y2 ","          ","          ", &
     !&  "fy(3x2-y2)","    fxyz  ","    fyz2  ","    fz3   ","    fxz2  "," fz(x2-y2)","fx(x2-3y2)" /)
      CHARACTER(LEN=10) :: LABEL(1:7,-5:3)
      DATA LABEL(1:7,-5) &
     &/"  sp3d2-1 ","  sp3d2-2 ","  sp3d2-3 ","  sp3d2-4 ","  sp3d2-5 ","  sp3d2-6 ","          "/
      DATA LABEL(1:7,-4) &
     &/"   sp3d-1 ","   sp3d-2 ","   sp3d-3 ","   sp3d-4 ","   sp3d-5 ","          ","          "/
      DATA LABEL(1:7,-3) &
     &/"   sp3-1  ","   sp3-2  ","   sp3-3  ","   sp3-4  ","          ","          ","          "/
      DATA LABEL(1:7,-2) &
     &/"   sp2-1  ","   sp2-2  ","   sp2-3  ","          ","          ","          ","          "/
      DATA LABEL(1:7,-1) &
     &/"    sp-1  ","    sp-2  ","          ","          ","          ","          ","          "/
      DATA LABEL(1:7, 0) &
     &/"      s   ","          ","          ","          ","          ","          ","          "/
      DATA LABEL(1:7, 1) &
     &/"     py   ","     pz   ","     px   ","          ","          ","          ","          "/
      DATA LABEL(1:7, 2) &
     &/"    dxy   ","    dyz   ","    dz2   ","    dxz   ","   dx2-y2 ","          ","          "/
      DATA LABEL(1:7, 3) &
     &/"fy(3x2-y2)","    fxyz  ","    fyz2  ","    fz3   ","    fxz2  "," fz(x2-y2)","fx(x2-3y2)"/

      CHARACTER(LEN=15) :: PRKIND(3)=(/"PAW projector  ","PS partial wave","Hydrogen-like  "/)

      CHARACTER (LEN=40) :: strcounter
 
      IF (.NOT.ALLOCATED(LPRJ_COVL)) RETURN

      NB=SIZE(LPRJ_COVL,1)
      NK=SIZE(LPRJ_COVL,2)
      NS=SIZE(LPRJ_COVL,3)
      NF=SIZE(LPRJ_COVL,4)


      IF (IU0>=0) THEN
!=======================================================================
         OPEN(UNIT=99,FILE='PROJCAR',STATUS='REPLACE')

         ALLOCATE(LDONE(NF), LDONOW(NF))
         LDONE=.FALSE.

         funcs1: DO IBASE=1,NF

            ! skip entries already printed
            samesite: IF (.NOT.LDONE(IBASE)) THEN
            
               ! print site and type
               WRITE(99,'(3X,A,I6,2X,A,3F14.7,2X,A,2X,A,A)',ADVANCE='No') &
              &  'ISITE:',LPRJ_functions(IBASE)%poscar_site,'  R= ',LPRJ_functions(IBASE)%R, &
              &   PRKIND(LPRJ_functions(IBASE)%radial_type),':'
       
               IF (LPRJ_functions(IBASE)%radial_type==3) THEN
                  WRITE(99,'(2X,A,I2,2X,A,F8.4)') &
                 &   'n=',LPRJ_functions(IBASE)%n,'za=',LPRJ_functions(IBASE)%za
               ELSE
                  WRITE(99,'(2X,A,I2)') &
                 &   'species  =',LPRJ_functions(IBASE)%species
               ENDIF
       
               ! now collect all functions on this site with a similar "recipy":
               ! for radial_type==3 this means all functions at the same site of
               ! the same type and with the same n and za. In case radial_type/=3
               ! all functions on the same site of the same type and the same POTCAR
               ! species are grouped together.
               LDONOW=.FALSE.
               DO IFNC=IBASE,NF
                  ! same site as before
                  DR=ABS(LPRJ_functions(IBASE)%R(1)-LPRJ_functions(IFNC)%R(1))+ &
                       &   ABS(LPRJ_functions(IBASE)%R(2)-LPRJ_functions(IFNC)%R(2))+ &
                       &   ABS(LPRJ_functions(IBASE)%R(3)-LPRJ_functions(IFNC)%R(3))
                  
                  IF (.NOT. LDONE(IFNC) .AND.  &
                       & DR<TINYR.AND. &
                       & LPRJ_functions(IBASE)%radial_type==LPRJ_functions(IFNC)%radial_type.AND. &
                       & LPRJ_functions(IBASE)%species==LPRJ_functions(IFNC)%species.AND. &
                       & LPRJ_functions(IBASE)%n==LPRJ_functions(IFNC)%n.AND. &
                       & LPRJ_functions(IBASE)%za==LPRJ_functions(IFNC)%za) THEN
                     LDONOW(IFNC)=.TRUE.
                  ENDIF
               ENDDO
       
               spin1: DO ISP=1,NS
               kpoints1: DO IK=1,NK
                  WRITE(99,'(/X,A,I6,2X,A,I4/)') 'k-point:',IK,'spin:',ISP
       
                  WRITE(99,'(3X,A)',ADVANCE='No') 'band'
                  DO IFNC=IBASE,NF
                     IF (LDONOW(IFNC)) &
                     WRITE(99,'(4X,A,X)',ADVANCE='No') LABEL(LPRJ_functions(IFNC)%m,LPRJ_functions(IFNC)%l)
                  ENDDO
                  WRITE(99,*)
       
                  DO IB=1,NB
                     WRITE(99,'(X,I6,X)',ADVANCE='No') IB
                     DO IFNC=IBASE,NF
                        IF (LDONOW(IFNC)) &
                        WRITE(99,'(2X,F6.3,X,F6.3)',ADVANCE='No') GCONJG(LPRJ_COVL(IB,IK,ISP,IFNC))
                     ENDDO
                     WRITE(99,*)
                  ENDDO
       
               ENDDO kpoints1
               ENDDO spin1
       
               LDONE=LDONE.OR.LDONOW

               WRITE(99,*)
            ENDIF samesite
         ENDDO funcs1

         DEALLOCATE(LDONE, LDONOW)
         CLOSE(99)
!=======================================================================
         OPEN(UNIT=99,FILE='LOCPROJ',STATUS='REPLACE')
         WRITE(99,'(4I6,"  # of spin, # of k-points, # of bands, # of proj" )') NS,NK,NB,NF
         DO IFNC=1,NF
            WRITE(99,'(3X,A,I6,2X,A,3F14.7,2X,A,2X,A,A)',ADVANCE='No') &
           &  'ISITE:',LPRJ_functions(IFNC)%poscar_site,'  R= ',LPRJ_functions(IFNC)%R, &
           &   PRKIND(LPRJ_functions(IFNC)%radial_type),':'
            WRITE(99,'(A)') LABEL(LPRJ_functions(IFNC)%m,LPRJ_functions(IFNC)%l)
         ENDDO
         WRITE(99,*)
         spin2: DO ISP=1,NS
            kpoints2: DO IK=1,NK
               IK_IRZ=KPOINTS_FULL%NEQUIV(IK)
               bands2: DO IB=1,NB
                  WRITE(99,'("orbital",3I6,2F20.10)') ISP,IK,IB,REAL(W%CELTOT(IB,IK_IRZ,ISP),q),W%FERTOT(IB,IK_IRZ,ISP)
                  funcs2: DO IFNC=1,NF
                     WRITE(99,'(1I6,2F20.10)') IFNC,GCONJG(LPRJ_COVL(IB,IK,ISP,IFNC))
                  ENDDO funcs2
                  WRITE(99,*)
               ENDDO bands2
            ENDDO kpoints2
         ENDDO spin2
         CLOSE(99)
      ENDIF

!=======================================================================
      CALL XML_TAG("locprojected")
      funcs3: DO IFNC=1,NF
         WRITE(strcounter,"(I6)") IFNC
         CALL XML_TAG("set", comment="function "//TRIM(ADJUSTL(strcounter)))

         ITMP(1)=LPRJ_functions(IFNC)%n; ITMP(2)=LPRJ_functions(IFNC)%l; ITMP(3)=LPRJ_functions(IFNC)%m

         CALL XML_TAG_INT("radial_type",LPRJ_functions(IFNC)%radial_type)
         CALL XML_TAG_INT("species",LPRJ_functions(IFNC)%species)
         CALL XML_VEC_INT(ITMP,"nlm")
         CALL XML_TAG_REAL("za",LPRJ_functions(IFNC)%za)
         CALL XML_TAG_INT("poscar_site",LPRJ_functions(IFNC)%poscar_site)
         CALL XML_VEC_REAL(LPRJ_functions(IFNC)%R,"position")
         CALL XML_TAG_STRING("label",TRIM(ADJUSTL(LABEL(LPRJ_functions(IFNC)%m,LPRJ_functions(IFNC)%l))))

         IF (LPRJ_functions(IFNC)%radial_type==1) &
            CALL XML_VEC_REAL(LPRJ_functions(IFNC)%optproj,"optproj")

         spin3: DO ISP=1,NS
            WRITE(strcounter,"(I2)") ISP
            CALL XML_TAG("set", comment="spin "//TRIM(ADJUSTL(strcounter)))
      
            kpoints3: DO IK=1,NK
               WRITE(strcounter,"(I6)") IK
               CALL XML_TAG("set", comment="kpoint "//TRIM(ADJUSTL(strcounter)))
      
               DO IB=1,NB
#ifdef gammareal
                  RTMP(1)=LPRJ_COVL(IB,IK,ISP,IFNC)
#else
                  RTMP(1)=AIMAG(LPRJ_COVL(IB,IK,ISP,IFNC))
                  RTMP(2)=REAL(LPRJ_COVL(IB,IK,ISP,IFNC),q)
#endif
                  CALL XML_ROW_DATA(RTMP, FORM='(F14.7,1X)')
               ENDDO
               CALL XML_CLOSE_TAG("set")
            ENDDO kpoints3
            CALL XML_CLOSE_TAG("set")
         ENDDO spin3
         CALL XML_CLOSE_TAG("set")
      ENDDO funcs3
      CALL XML_CLOSE_TAG("locprojected")

      END SUBROUTINE LPRJ_WRITE

!***********************************************************************
!> @brief writes U_{nm}^k tranformation matrix between
!  Bloch and Wannier functions to file
!***********************************************************************
      SUBROUTINE WRITE_WANPROJ_FILE(WANPROJFILE,NB_TOT,ISPIN,NKPTS,NW,VKPT,U)
         USE wave_struct_def
         USE full_kpoints
         USE string, ONLY: str
         USE tutor, ONLY: vtutor
         IMPLICIT NONE
         CHARACTER(LEN=*),INTENT(IN) :: WANPROJFILE !< Name of the file in which to write the projections
         INTEGER,INTENT(IN) :: NB_TOT !< Number of bands
         INTEGER,INTENT(IN) :: ISPIN !< Number of spins
         INTEGER,INTENT(IN) :: NKPTS !< Number of k-points
         INTEGER,INTENT(IN) :: NW !< Number of Wannier orbitals
         REAL(q),INTENT(IN) :: VKPT(3,NKPTS) !< Coordinates of the k-points
         !> U_(nm)^K rotation matrix (WDES%NB_TOT,NW,WDES%NKPTS,WDES%ISPIN)
         GDEF,ALLOCATABLE,INTENT(IN) :: U(:,:,:,:)
         ! local variables
         INTEGER I,J
         INTEGER ISP,NK

         IF ( NB_TOT/=SIZE(U,1) ) THEN
            CALL vtutor%error("WRITE_WANPROJ_FILE: Incompatible number of bands "//STR(NB_TOT)//" != "//STR(SIZE(U,1)))
         ENDIF
         IF ( NW/=SIZE(U,2) ) THEN
            CALL vtutor%error("WRITE_WANPROJ_FILE: Incompatible number of projections "//STR(NW)//" != "//STR(SIZE(U,2)))
         ENDIF
         IF ( NKPTS/=SIZE(U,3) ) THEN
            CALL vtutor%error("WRITE_WANPROJ_FILE: Incompatible number of k-points "//STR(NKPTS)//" != "//STR(SIZE(U,3)))
         ENDIF
         IF ( ISPIN/=SIZE(U,4) ) THEN
            CALL vtutor%error("WRITE_WANPROJ_FILE: Incompatible number of spins "//STR(ISPIN)//" != "//STR(SIZE(U,4)))
         ENDIF

         OPEN(UNIT=99, FILE=WANPROJFILE, STATUS='REPLACE')

         WRITE(99,*)'# This file was created by VASP, do not edit it!'
         WRITE(99,*)ISPIN,NKPTS,NB_TOT,NW
         DO NK=1,NKPTS
            WRITE(99,'(I6, 3F26.16)') NK,VKPT(:,NK)
         ENDDO

         spin: DO ISP=1,ISPIN
            kpts: DO NK=1,NKPTS
               WRITE(99,'(2I6,3F26.16)')ISP,NB_TOT,VKPT(:,NK)

               ! setup the effective rotation matrix U at the present k-point
               DO I=1,NB_TOT
                  DO J=1,NW
                     WRITE(99,'(2I6,2F26.16)')I,J,U(I,J,NK,ISP)
                  ENDDO
               ENDDO

            ENDDO kpts
         ENDDO spin
         CLOSE(99)

      END SUBROUTINE WRITE_WANPROJ_FILE


!******************** SUBROUTINE LPRJ_LDApU ****************************
!
! this is a kind of hacked in version of LDA+U
! instead of doing LDA+U directly, this bit of code writes
! a density matrix file GAMMA that can be used to correct
! the density matrix (using ICHARG=5)
!
!***********************************************************************

      SUBROUTINE LPRJ_LDApU(IO, W)
      USE base
      USE wave
      USE full_kpoints
      USE fileio
      TYPE (in_struct)   IO
      TYPE (wavespin) W
    ! local
      INTEGER ISP, NS, NK, IK, IK_IRZ, IB, IB2, NB, IFNC, IFNC2, I, NF
      REAL(q) DR
      REAL(q), PARAMETER :: TINYR=1.E-4_q
      GDEF, ALLOCATABLE :: DENSITY_MATRIX(:,:,:), U(:,:,:), GAMMA(:,:)
!#define testgamma_matrix
#ifdef testgamma_matrix
  ! LAPACK    
      REAL(q)    R(W%WDES%NB_TOT)
      INTEGER :: IFAIL
      INTEGER, PARAMETER :: LWORK=32
      GDEF       CWRK(LWORK*W%WDES%NB_TOT)
      REAL(q)    RWORK(3*W%WDES%NB_TOT)
#endif
      
      IF (.NOT.ALLOCATED(LPRJ_COVL)) RETURN
      
      NB=SIZE(LPRJ_COVL,1)
      NK=SIZE(LPRJ_COVL,2)
      NS=SIZE(LPRJ_COVL,3)
      NF=SIZE(LPRJ_COVL,4)

      ALLOCATE(DENSITY_MATRIX(NF,NF, NS), U(NF,NF, NS))
      DENSITY_MATRIX=0

      ! the first step is to determine the occupancy matrix in the
      ! space of the projected functions
      DO ISP=1,NS
         DO IK=1,NK
            IK_IRZ=KPOINTS_FULL%NEQUIV(IK)
            DO IB=1,NB
               DO IFNC=1,NF
               DO IFNC2=1,NF
                  DENSITY_MATRIX(IFNC,IFNC2,ISP)=DENSITY_MATRIX(IFNC,IFNC2,ISP)+ & 
                       LPRJ_COVL(IB,IK,ISP,IFNC)*GCONJG(LPRJ_COVL(IB,IK,ISP,IFNC2))*W%FERTOT(IB,IK_IRZ,ISP)
               ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DENSITY_MATRIX=DENSITY_MATRIX/KPOINTS_FULL%NKPTS
      ! now work on one site only
      U=0
      IFNC=1
      DO
         DO I=IFNC,NF
            ! same site as before
            DR=ABS(LPRJ_functions(IFNC)%R(1)-LPRJ_functions(I)%R(1))+ &
                 &   ABS(LPRJ_functions(IFNC)%R(2)-LPRJ_functions(I)%R(2))+ &
                 &   ABS(LPRJ_functions(IFNC)%R(3)-LPRJ_functions(I)%R(3))
            IF (DR<TINYR) THEN
               IFNC2=I
            ENDIF
         ENDDO
         ! for each site and spin calculate (1/2-n)
         DO ISP=1,NS
            U(IFNC:IFNC2,IFNC:IFNC2,ISP)=-DENSITY_MATRIX(IFNC:IFNC2,IFNC:IFNC2,ISP)
            DO I=IFNC,IFNC2
               U(I,I,ISP)=0.5_q+U(I,I,ISP)
            ENDDO
         ENDDO
         IFNC=IFNC2+1
         IF (IFNC2>=NF) EXIT
      ENDDO

    ! hard coded U
      U=U*7

#ifdef testgamma_matrix
      DO ISP=1,NS
         WRITE(*,*) 'density matrix'
         WRITE(*,'(10F8.3)') REAL(DENSITY_MATRIX(:,:,ISP),q)
      ENDDO

      DO ISP=1,NS
         WRITE(*,*) 'one-center U matrix'
         WRITE(*,'(10F8.3)') REAL(U(:,:,ISP),q)
      ENDDO
#endif
    ! finally construct the correction to the Hamiltonian and
    ! write it to the file GAMMA
      CALL OPENGAMMA
      CALL WRITEGAMMA_HEAD_HAM( W%WDES%NKPTS, W%WDES%NB_TOT, IO)

      ALLOCATE(GAMMA(W%WDES%NB_TOT,W%WDES%NB_TOT))
      
      DO ISP=1,NS
         DO IK=1,W%WDES%NKPTS
            GAMMA=0
            DO IB=1,NB
            DO IB2=1,NB
               DO IFNC=1,NF
               DO IFNC2=1,NF
                  GAMMA(IB,IB2)=GAMMA(IB,IB2)+GCONJG(LPRJ_COVL(IB,IK,ISP,IFNC))*LPRJ_COVL(IB2,IK,ISP,IFNC2)* & 
                       U(IFNC,IFNC2,ISP)
               ENDDO
               ENDDO
            ENDDO
            ENDDO
            CALL WRITEGAMMA(IK, W%WDES%NB_TOT, W%WDES%NB_TOT, GAMMA, .FALSE.)
            DO IB=1,W%WDES%NB_TOT
               GAMMA(IB,IB)=GAMMA(IB,IB)+REAL(W%CELTOT(IB,IK,ISP),q)
            ENDDO
#ifdef testgamma_matrix
            IF (IK==1 .AND. ISP==1) THEN
               WRITE(*,'(8F14.7)') (REAL(GAMMA(I,I),q),I=1, W%WDES%NB_TOT)
            ENDIF
#ifdef  gammareal
            CALL DSYEV &
              ('V','U',NB,GAMMA(1,1),W%WDES%NB_TOT, &
              R,CWRK,LWORK*WDES%NB_TOT, IFAIL)
#else
            CALL ZHEEV &
              ('V','U',NB,GAMMA(1,1),W%WDES%NB_TOT, &
              R,CWRK,LWORK*W%WDES%NB_TOT, RWORK,  IFAIL)
#endif
            IF (IK==1 .AND. ISP==1) THEN
               WRITE(*,'(8F14.7)') (R(I),I=1, W%WDES%NB_TOT)
               WRITE(*,*)
            ENDIF
#endif
         ENDDO
      ENDDO
      CALL CLOSEGAMMA

      DEALLOCATE(GAMMA)
      DEALLOCATE(DENSITY_MATRIX,U)

      END SUBROUTINE LPRJ_LDApU


!******************** SUBROUTINE LPRJ_DEALLOC *********************
!> Deallocate variables that care created by calls to LPRJ_PROALL
!***********************************************************************
      SUBROUTINE LPRJ_DEALLOC()
      IF (ALLOCATED(LPRJ_COVL)) DEALLOCATE(LPRJ_COVL)
      RETURN
      END SUBROUTINE LPRJ_DEALLOC


!******************** LPRJ_LCAO ****************************************
!
!> returns the value of locproj::locproj_needs_ps
!! In case locproj::lprj_lcao=.TRUE. VASP will call lcao::lcao_init
!! (from main.F).
!
!***********************************************************************

      FUNCTION LPRJ_LCAO()
      IMPLICIT NONE
      LOGICAL LPRJ_LCAO
      LPRJ_LCAO=LOCPROJ_NEEDS_PS
      END FUNCTION LPRJ_LCAO


!******************** SUBROUTINE CALC_OVERLAP_GN ***********************
!      
!> @brief Calculate \f$S(nb,k,i) = < \tilde \psi_{nb,k} | S | RYlm_i >\f$
!! where i runs over a set of local functions at site R=POS with
!! l = [LSTART,LSTOP], and m = [1,2*l+1], and nb runs over all
!! orbitals at k-point k.
!!
!! For LNORM=.TRUE. the local functions themselves are S-normalized.
!! For LAUG=.FALSE. the augmentation part in the contraction of S(nb,i)
!! set to zero, the overlap operator is taken to be the identity
!! operator, S=1.
!
!***********************************************************************
      SUBROUTINE CALC_OVERLAP_GN( LSTART,LSTOP,FG,POS,W,K,ISP,ISPINOR,SPIN_QAXIS,P,CQIJ,LATT_CUR,T_INFO,INFO,LNORM,LAUG,S )
      USE pead, ONLY : WDES_FULL_PEAD,PEAD_WA_ROTATE
      USE base, ONLY : info_struct
      USE poscar
      USE pseudo
      USE lattice
      USE full_kpoints
      USE wave_high
      USE nonl_high
      IMPLICIT NONE
      TYPE(wavespin) W
      TYPE(potcar) P(:)
      TYPE(latt) LATT_CUR
      TYPE(type_info) T_INFO
      TYPE(info_struct) INFO
      INTEGER LSTART,LSTOP
      INTEGER ISP
      INTEGER ISPINOR
      REAL(q) K(3)
      REAL(q) FG(:,:,:)
      REAL(q) POS(3)
      OVERLAP CQIJ(:,:,:,:)
      GDEF S(:,:)
      LOGICAL LNORM,LAUG
      ! local variables
      TYPE(wavespin) WP
      TYPE(wavefuna) WK,WRYLM
      TYPE(wavedes1), TARGET :: WDESK
      TYPE(nonl_struct) NONL_S
      TYPE(nonlr_struct) NONLR_S
      LOGICAL :: LREALLOCATE

      TYPE(rotation_handle), POINTER :: ROT_HANDLE

      GDEF C
      REAL(q) WSCAL
      REAL(q) SPIN_QAXIS(3)
      INTEGER NK,NB,N,NYLM

      PROFILING_START('calc_overlap_gn')
      
      WP=W
      WP%WDES=>WDES_FULL_PEAD
            
      CALL CHECK_FULL_KPOINTS

      NULLIFY(ROT_HANDLE)      

      NYLM=(LSTOP+1)**2-LSTART**2

      ! search for kpoint k in BZ
      NK=KPOINT_IN_FULL_GRID(K,KPOINTS_FULL)
      CALL SETWDES(WP%WDES,WDESK,NK)
      IF (NK==KPOINTS_FULL%NEQUIV(NK)) THEN
         ! k is a kpoint in the IBZ
         WK=ELEMENTS(WP,WDESK,ISP)
      ELSE
         ! k is not a kpoint in the IBZ
         CALL NEWWAVA(WK,WDESK,WDESK%NBANDS)
         CALL PEAD_WA_ROTATE(WP,P,LATT_CUR,ISP,WK)
      ENDIF

      NONLR_S%LREAL = INFO%LREAL
      IF (LAUG) THEN
         IF (INFO%LREAL) THEN
            CALL NONLR_SETUP(NONLR_S,T_INFO,P,INFO%LREAL,WP%WDES%LSPIRAL)
#ifndef shmem_rproj
            CALL NONLR_SETUP(NONLR_S,T_INFO,P,INFO%LREAL,WP%WDES%LSPIRAL)
#else
            CALL NONLR_SETUP(NONLR_S,T_INFO,P,INFO%LREAL,WP%WDES%LSPIRAL,WP%WDES%COMM_intra_node)
#endif
            CALL REAL_OPTLAY(WP%WDES%GRID,LATT_CUR,NONLR_S,.TRUE.,LREALLOCATE, -1,-1)
            CALL NONLR_ALLOC(NONLR_S)
            CALL RSPHER(WP%WDES%GRID,NONLR_S,LATT_CUR)
            CALL PHASER(WP%WDES%GRID,LATT_CUR,NONLR_S,NK,WP%WDES)
         ELSE
            CALL NONL_ALLOC(NONL_S,T_INFO,P,WP%WDES,.FALSE.)
            CALL SPHER(WP%WDES%GRID,NONL_S,P,WP%WDES,LATT_CUR,1,NK)
            CALL PHASE(WP%WDES,NONL_S,NK)
         ENDIF
      ENDIF
      CALL NEWWAVA(WRYLM,WDESK,NYLM); WRYLM%CW=0

      CALL CONSTRUCT_RYLM(LSTART,LSTOP,FG,LATT_CUR,POS,NONLR_S,NONL_S,ISPINOR,SPIN_QAXIS,WRYLM,LAUG)
      ! and normalize the functions WRYLM
      IF (LNORM) THEN
         DO N=1,NYLM
            CALL CNORMN(ELEMENT(WRYLM,N),CQIJ,1,WSCAL)
         ENDDO
      ENDIF

      IF (.NOT.LAUG) WRYLM%CPROJ=0

      ! calculate overlap between Wk and WRYLM: < w_{m,k1} | S | RYlm >
      S=0    
      DO NB=1,WP%WDES%NBANDS
         DO N=1,NYLM
            C=W1_DOT(ELEMENT(WK,NB),ELEMENT(WRYLM,N),CQIJ)
            S(WDESK%NB_LOW+WDESK%NB_PAR*(NB-1),N)=C
!           WRITE(*,*) WDESK%NB_LOW+WDESK%NB_PAR*(NB-1),N,C
         ENDDO
      ENDDO

      CALLMPI( M_sum_g(WDESK%COMM_INTER,S(1,1),WDESK%NB_TOT*NYLM) )
      
      ! some deallocation to be done      
      CALL DEALLOCATE_ROT_HANDLE(ROT_HANDLE)

      CALL DELWAVA(WRYLM)
      IF (NK/=KPOINTS_FULL%NEQUIV(NK)) CALL DELWAVA(WK)      

      IF (LAUG) THEN
         IF (INFO%LREAL) THEN
            CALL NONLR_DEALLOC(NONLR_S)
         ELSE
            CALL NONL_DEALLOC(NONL_S)
         ENDIF
      ENDIF

      PROFILING_STOP('calc_overlap_gn')

      RETURN
      END SUBROUTINE CALC_OVERLAP_GN


!******************** SUBROUTINE CONSTRUCT_RYLM ***************
!
!> @brief Construct the local functions |RYlm> at site R=POS in reciprocal space:
!
! \f$< k+G | RYlm > = 1/sqrt(Omega) (-i)^l F(|k+G|) Y_lm(k+G) e^{-i(k+G)R}\f$
!
! (see Sec. 10.5.1 of the thesis of gK). The plane wave part of these
! local functions are stored in WRYLM%CW, and the projections onto
! the PAW projectors in WRYLM%CPROJ.
!
! This is done for a set of local functions with l = [LSTART,LSTOP],
! and m = [1, 2*l+1] (i.e., -l,..,l )
!
!***********************************************************************
      SUBROUTINE CONSTRUCT_RYLM( LSTART,LSTOP,F,LATT_CUR,POS,NONLR_S,NONL_S,ISPINOR,SPIN_QAXIS,WRYLM,LAUG )
      USE ini
      USE asa
      USE lattice
      USE constant
      USE wave_high
      USE nonl_high
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      USE spinsym, only: sigmadotn
      IMPLICIT NONE
      TYPE(latt) LATT_CUR
      TYPE(wavefuna) WRYLM
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      REAL(q) F(:,:,:) !< Spline coefficients of the radial part (NMAX,5,NYLM)
      REAL(q) POS(3)
      REAL(q) SPIN_QAXIS(3) !< Spin quantization axis
      INTEGER LSTART,LSTOP
      INTEGER ISPINOR !< spinor index for this orbital can be 1 (up) or 2 (down)
      LOGICAL LAUG
      ! local variables
      INTEGER N1,N2,N3,IND,NPL,LMMAX,L,M,LM,LMBASE,IG
      REAL(q) G1,G2,G3,GKX,GKY,GKZ,FACTM,FAKT,FDER,GMAX
      REAL(q), ALLOCATABLE :: XS(:),YS(:),ZS(:),YLM(:,:)
      REAL(q), ALLOCATABLE :: G(:),FG(:)
      COMPLEX(q) CSET,CGDR
      COMPLEX(q), ALLOCATABLE :: CFAKTX(:)
      TYPE(wavefun1) :: WS(1)
      ! spin
      COMPLEX(q) :: PV_MATRIX(2,2), SIG(2,2,3)
      COMPLEX(q) :: ZW(2), ZVL(2,2), ZVR(2,2), ZWORK(4)
      REAL(q) :: ALPHA, BETA
      REAL(q) :: RWORK(4)
      INTEGER :: ISPINOR_, INFO, PERM(2)

      PROFILING_START('construct_rylm')
      LMMAX=(LMAX+1)**2

      NPL=WRYLM%WDES1%NGVECTOR      

      IF (ISPINOR/=1.AND.(.NOT.WRYLM%WDES1%LNONCOLLINEAR)) THEN
         CALL vtutor%error("CONSTRUCT_FUNCTION_RYlm: ERROR: ISPINOR= " // str(ISPINOR) // " but &
            &LNONCOLLINEAR=.FALSE.")
      ENDIF

      ALLOCATE(G(NPL),FG(NPL))
      ALLOCATE(XS(NPL),YS(NPL),ZS(NPL),CFAKTX(NPL))

      ! loop over all G-vectors in the basis at this k-point
      DO IND=1,WRYLM%WDES1%NGVECTOR
         N1=MOD(WRYLM%WDES1%IGX(IND)+WRYLM%WDES1%GRID%NGX,WRYLM%WDES1%GRID%NGX)+1
         N2=MOD(WRYLM%WDES1%IGY(IND)+WRYLM%WDES1%GRID%NGY,WRYLM%WDES1%GRID%NGY)+1 
         N3=MOD(WRYLM%WDES1%IGZ(IND)+WRYLM%WDES1%GRID%NGZ,WRYLM%WDES1%GRID%NGZ)+1

         G1=(WRYLM%WDES1%GRID%LPCTX(N1)+WRYLM%WDES1%VKPT(1))
         G2=(WRYLM%WDES1%GRID%LPCTY(N2)+WRYLM%WDES1%VKPT(2))
         G3=(WRYLM%WDES1%GRID%LPCTZ(N3)+WRYLM%WDES1%VKPT(3))

         FACTM=1._q
         IF (WRYLM%WDES1%LGAMMA .AND. (N1/=1 .OR. N2/=1 .OR. N3/=1)) FACTM=SQRT(2._q)

         GKX=(G1*LATT_CUR%B(1,1)+G2*LATT_CUR%B(1,2)+G3*LATT_CUR%B(1,3))*TPI
         GKY=(G1*LATT_CUR%B(2,1)+G2*LATT_CUR%B(2,2)+G3*LATT_CUR%B(2,3))*TPI
         GKZ=(G1*LATT_CUR%B(3,1)+G2*LATT_CUR%B(3,2)+G3*LATT_CUR%B(3,3))*TPI

         G(IND)=MAX(SQRT(GKX*GKX+GKY*GKY+GKZ*GKZ),1E-10_q)

         ! phase factor e^{-i(k+G)R} where R is the origin
         ! of the localized function
         CGDR=CITPI*(G1*POS(1)+G2*POS(2)+G3*POS(3))
! test_
!        CGDR=CITPI*(WRYLM%WDES1%GRID%LPCTX(N1)*POS(1)+WRYLM%WDES1%GRID%LPCTY(N2)*POS(2)+WRYLM%WDES1%GRID%LPCTZ(N3)*POS(3))
!        CGDR=0
! test_
         CFAKTX(IND)=FACTM*EXP(-CGDR)


         XS(IND)  =GKX/G(IND)
         YS(IND)  =GKY/G(IND)
         ZS(IND)  =GKZ/G(IND)
      ENDDO

      ALLOCATE(YLM(NPL,LMMAX))
      ! get me all the Y_lm up to and including l=LMAX
      CALL SETYLM(LMAX,NPL,YLM,XS,YS,ZS)

      ! Setup the plane wave part of the desired function
      FAKT= 1/SQRT(LATT_CUR%OMEGA)
      CSET=CMPLX(0._q,-1._q,q)

#if 0
      ! Diagonalize the Pauli vector times quantization axis vector dot(SPIN_QAXIS,PAULI_VECTOR)
      ! to obtain the eigenvectors and eigenvalues.
      ! The PAULI_VECTOR is (Sx,Sy,Sx) where Sx,Sy,Sz are the Pauli matrices
      SPIN_QAXIS = SPIN_QAXIS/SQRT(SUM(SPIN_QAXIS**2))
      PV_MATRIX = SIGMADOTN(SPIN_QAXIS)
      CALL ZGEEV( 'N', 'V', 2, PV_MATRIX, 2, ZW, ZVL, 1, &
           ZVR, 2, ZWORK, 4, RWORK, INFO)
      ! debug: Verify that these are eigenvectors
      !PV_MATRIX = SIGMADOTN(SPIN_QAXIS)
      !WRITE(*,*) DOT_PRODUCT(ZVR(:,1),MATMUL(PV_MATRIX,ZVR(:,1)))
      !WRITE(*,*) DOT_PRODUCT(ZVR(:,2),MATMUL(PV_MATRIX,ZVR(:,2)))
#else
      ! Use spin rotation matrix to rotate quantization axis.
      ! The rotation matrix is obtained from two spin rotations.
      ! A spin rotation is defined by an angle and an axis
      ! R_s(theta,axis) = exp(-I*theta/2*dot(axis,pauli_vector))
      !
      ! The two rotations are the same as the EULER routine:
      ! 1. rotation around y-axis of beta
      ! R_s(beta,y) = [ cos(beta/2)   -sin(beta/2) ]
      !               [ sin(beta/2     cos(beta/2) ]
      !
      ! 2. rotation around z-axis of alpha
      ! R_s(alpha,z) = [ exp(-I*alpha/2)              0 ]
      !                [               0 exp(I*alpha/2) ]
      !
      ! The spin rotation is given by
      ! R = R_s(alpha,z)*R_s(beta,y)
      !
      ! Rotating the spin up and down states of Sz yields:
      ! ZVR(:,1) = R [1,0]
      ! ZVR(:,2) = R [0,1]
      !
      ! This is equivalent to the case above because:
      ! R Sz Dagger(R) = dot(SPIN_QAXIS,PAULI_VECTOR)
      ! i.e. rotating the Sz operator we obtain the same operator that we diagonalize above
      ! This is equivalent to rotating the eigenvectors of the Sz operator (stored in ZVR)
      CALL EULER(SPIN_QAXIS,ALPHA,BETA)
      ZVR(:,1) = [ COS(BETA/2)*EXP(-CMPLX(0.0_q,ALPHA/2,q)),SIN(BETA/2)*EXP(CMPLX(0.0_q,ALPHA/2,q))]
      ZVR(:,2) = [-SIN(BETA/2)*EXP(-CMPLX(0.0_q,ALPHA/2,q)),COS(BETA/2)*EXP(CMPLX(0.0_q,ALPHA/2,q))]
      ! debug: Verify that these are eigenvectors
      !PV_MATRIX = SIGMADOTN(SPIN_QAXIS)
      !WRITE(*,*) DOT_PRODUCT(ZVR(:,1),MATMUL(PV_MATRIX,ZVR(:,1)))
      !WRITE(*,*) DOT_PRODUCT(ZVR(:,2),MATMUL(PV_MATRIX,ZVR(:,2)))
#endif

      WRYLM%CW=0
      LMBASE=LSTART**2+1; LM=1
      DO L=LSTART,LSTOP
         ! get me the Bessel transform of the radial function
         GMAX=F(SIZE(F,1),1,L+1); FG=0
         DO IG=1,SIZE(G)
            IF (G(IG)>=GMAX) CYCLE
            CALL SPLVAL(G(IG),FG(IG),FDER,F(:,:,L+1),SIZE(F,1),SIZE(F,1))
         ENDDO
         DO M=1,2*L+1
            DO ISPINOR_=1,WRYLM%WDES1%NRSPINORS
               DO IND=1,NPL
                  WRYLM%CW(IND+(ISPINOR_-1)*NPL,LM+M-1)=ZVR(ISPINOR_,ISPINOR)*FAKT*(CSET**L)*CFAKTX(IND)*FG(IND)*YLM(IND,LMBASE+M-1)
               ENDDO ! plane-waves
            ENDDO ! spinor
         ENDDO ! m magnetic number
         LM=LM+2*L+1
         LMBASE=LMBASE+2*L+1
      ENDDO

      ! and get the projections of RYlm onto the PAW projectors
      WRYLM%CPROJ=0
      IF (LAUG) THEN
         DO LM=1,SIZE(WRYLM%CPROJ,2)
            WS(1)=ELEMENT(WRYLM,LM)
            IF ( NONLR_S%LREAL ) THEN
               ALLOCATE(WS(1)%CR(WRYLM%WDES1%GRID%MPLWV*WRYLM%WDES1%NRSPINORS))
               CALL FFTWAV_W1(WS(1))
            ENDIF
            CALL W1_PROJALL(WRYLM%WDES1,WS,NONLR_S,NONL_S)
            IF ( NONLR_S%LREAL ) DEALLOCATE(WS(1)%CR)
         ENDDO
      ENDIF

      DEALLOCATE(G,FG,XS,YS,ZS,CFAKTX,YLM)

      PROFILING_STOP('construct_rylm')

      RETURN
      END SUBROUTINE CONSTRUCT_RYLM


!******************** SUBROUTINE CALC_PSPNL_SPLINE *********************
!      
!***********************************************************************
      SUBROUTINE CALC_PSPNL_SPLINE(P,PSPNL_SPLINE)
      USE pseudo

      TYPE (potcar) P(:)
      REAL(q), ALLOCATABLE, INTENT(OUT) :: PSPNL_SPLINE(:,:,:,:)

      ! local variables
      INTEGER ITYP,IC,I

      ALLOCATE(PSPNL_SPLINE(0:NPSNL,5,MAXVAL(P(:)%LDIM),SIZE(P)))

      DO ITYP=1,SIZE(P)
         DO IC=1,P(ITYP)%LMAX
            DO I=0,NPSNL
               PSPNL_SPLINE(I,1,IC,ITYP)=P(ITYP)%PSMAXN/NPSNL*(I-1)
            ENDDO
            PSPNL_SPLINE(:,2,IC,ITYP)=P(ITYP)%PSPNL(:,IC)

            CALL SPLCOF(PSPNL_SPLINE(0,1,IC,ITYP),NPSNL+1,NPSNL+1,1.E+30_q)
         ENDDO
      ENDDO

      END SUBROUTINE CALC_PSPNL_SPLINE


!******************** FUNCTION FIND_ATOM_LCAO_ENTRY ********************
!> Choose LCAO entry from which to build localized orbitals
!> Currently only occupied states can be selected
!***********************************************************************
      SUBROUTINE FIND_ATOM_LCAO_ENTRY(ITYP,L,ICHANNEL,LCAO_POSSIBLE_ENTRIES,LCAO_ENTRY)
      USE lcao, ONLY : ATOM_LCAO
      INTEGER :: ITYP !< type of the atom for which to find projectors
      INTEGER :: L !< Quantum number
      INTEGER :: ICHANNEL !< projector index to find
      INTEGER, ALLOCATABLE, INTENT(OUT) :: LCAO_POSSIBLE_ENTRIES(:)  !< returns angular momentum and iproj
      INTEGER,INTENT(OUT) :: LCAO_ENTRY !< returns the entry if found or -1 if not found
      ! local variables
      INTEGER :: NOCCUPIED_STATES
      INTEGER :: NOCCUPIED_STATES_WITHL
      INTEGER :: I
      LCAO_ENTRY=-1
      NOCCUPIED_STATES=ATOM_LCAO(ITYP)%OCCUPIED_VALENCE_STATES
      DO I=1,NOCCUPIED_STATES
         IF (ATOM_LCAO(ITYP)%L(I,4)==ICHANNEL) THEN
            LCAO_ENTRY=I; EXIT
         ENDIF
      ENDDO
      ! we did not find the desired projector
      ! need to warn the user and recommend what to do
      IF (LCAO_ENTRY==-1) THEN
         NOCCUPIED_STATES_WITHL = COUNT(ATOM_LCAO(ITYP)%L(:,1)==L)
         ALLOCATE(LCAO_POSSIBLE_ENTRIES(NOCCUPIED_STATES_WITHL))
         NOCCUPIED_STATES_WITHL=0
         DO I=1,NOCCUPIED_STATES
           IF (ATOM_LCAO(ITYP)%L(I,1)/=L) CYCLE
           ! store the index of the projector
           NOCCUPIED_STATES_WITHL=NOCCUPIED_STATES_WITHL+1
           LCAO_POSSIBLE_ENTRIES(NOCCUPIED_STATES_WITHL) = ATOM_LCAO(ITYP)%L(I,4)
         ENDDO
      ENDIF
      END SUBROUTINE FIND_ATOM_LCAO_ENTRY

!******************** SUBROUTINE CALC_HYDRO_SPLINE *********************
!> Build hydrogenic type orbitals, bessel transformation and spline
!***********************************************************************
      SUBROUTINE CALC_HYDRO_SPLINE(N,ZA,LMAX,ENMAX,HYDRO_SPLINE,IO)
      USE constant
      USE radial
      USE base

      INTEGER N,LMAX
      REAL(q) ZA,ENMAX
      !> coefficients of spline function describing the radial part (NPOINTS,5,LMAX+1)
      REAL(q), ALLOCATABLE, INTENT(OUT) :: HYDRO_SPLINE(:,:,:)
      TYPE(in_struct) :: IO

      ! local variables
      INTEGER L
      REAL(q), ALLOCATABLE :: FTMP(:)
      REAL(q) :: FAC

      TYPE (rgrid) R
      REAL(q) :: RSTART,REND,H,R1,R2

      ALLOCATE(HYDRO_SPLINE(NMAX,5,LMAX+1))

      ! set up a logarithmic grid
      CALL RGRID_PARAMETERS(RSTART,REND,H)
      CALL SETRGRID(RSTART,REND,H,R)
      ALLOCATE(FTMP(R%NMAX))
      ! construct hydrogenic radial functions
      IF (1<=N.AND.N<=3) THEN
        FAC = 1.0_q
        CALL RADIAL_HYDROGEN_FUNCTION(N,FAC,R,ZA,FTMP)
        CALL RADIAL_CALC_MOMENT(FTMP,R,R1,R2)
        !CALL RADIAL_WRITE_INFO("radial hydrogen functions",N,ZA,R1,R2,IO)
      ! These are hidden options for the moment
      ELSE IF (100<=N.AND.N<200) THEN
        ! construct slater functions: r^*(NSTAR-1) * EXP(-ZA*r/NSTAR)
        CALL RADIAL_SLATER_FUNCTION(N-100,R,ZA,FTMP)
        CALL RADIAL_CALC_MOMENT(FTMP,R,R1,R2)
        !CALL RADIAL_WRITE_INFO("radial slater-type functions",N,ZA,R1,R2,IO)
      ELSE IF (200<=N) THEN
        ! contruct gaussians: r^(N) * exp(-ZA*r^2)
        CALL RADIAL_GAUSSIAN_FUNCTION(N-200,R,ZA,FTMP)
        CALL RADIAL_CALC_MOMENT(FTMP,R,R1,R2)
        !CALL RADIAL_WRITE_INFO("radial gaussian functions",N,ZA,R1,R2,IO)
      ELSE
        CALL vtutor%error('LOCPROJ: Unknown radial function type '//str(N))
      ENDIF

      ! compute bessel transformation of the radial function
      DO L=0,LMAX
         CALL BESSEL_TRANSFORM_RADIAL_FUNCTION(L,R,FTMP,REAL(SQRT(2._q*ENMAX/HSQDTM)/NMAX,KIND=q),HYDRO_SPLINE(:,:,L+1))
      ENDDO
 
      DEALLOCATE(FTMP)

      DEALLOCATE(R%R,R%SI); NULLIFY(R%R,R%SI)
      RETURN
      END SUBROUTINE CALC_HYDRO_SPLINE


!******************** SUBROUTINE RADIAL_FUNCTION ***********************
!> Build hydrogenic type radial functions
!***********************************************************************
      SUBROUTINE RADIAL_HYDROGEN_FUNCTION(ITYP,FAC,R,ZA,FR)
      USE radial, ONLY: rgrid
      INTEGER :: ITYP !< principal quantum number of the hydrogen orbital
      REAL(q) :: FAC !< conversion factor of the radial grid
      TYPE(rgrid) :: R !< logaritmic radial grid
      REAL(q) :: ZA !< alpha parameter
      REAL(q),INTENT(OUT) :: FR(:) !< value of the function f(r)
      ! debug
      REAL(q) :: R1, R2

      ! FAC changes the units of the radial grid see R(:)/FAC
      ! This was used to convert from Angstroem to bohr when FAC=AUTOA
      ! The units of ZA are inverse Angstroem so converting the units of the
      ! radial grid to Bohr here is strange.

      SELECT CASE(ITYP)

       CASE (1)
         FR(:)=2._q*ZA**(3._q/2._q)*EXP(-ZA*R%R(:)/FAC)
       CASE (2)
         FR(:)=1._q/(2._q*SQRT(2._q))*ZA**(3._q/2._q)*(2._q-ZA*R%R(:)/FAC)*EXP(-ZA*R%R(:)/FAC/2._q)
       CASE (3)
         FR(:)=SQRT(4._q/27._q)*ZA**(3._q/2._q)* &
        &   (1._q-2._q/3._q*ZA*R%R(:)/FAC+2._q/27._q*ZA*ZA*R%R(:)*R%R(:)/FAC/FAC)*EXP(-ZA*R%R(:)/FAC/3._q)
       CASE DEFAULT
         call vtutor%bug('RADIAL_HYDROGEN_FUNCTION: TYPE '//str(ITYP)//'does not exist',__FILE__,__LINE__)
      END SELECT

      ! Normalize
      FR=FR/SQRT(RADIAL_NORM2(FR,R))

      END SUBROUTINE RADIAL_HYDROGEN_FUNCTION

!******************** SUBROUTINE RADIAL_SLATER_FUNCTION ***********************
!> Build Slater type radial functions
!******************************************************************************
      SUBROUTINE RADIAL_SLATER_FUNCTION(ITYP,R,ZA,FR)
      USE radial, ONLY: rgrid
      TYPE(rgrid) R
      INTEGER ITYP
      REAL(q) ZA
      REAL(q) FR(:)
      ! local
      REAL(q) NSTAR

      ! Slater-type (according to Slater's rule)
      NSTAR = 0.0_q
      SELECT CASE(ITYP)
         CASE(1)
           NSTAR = 1.0_q
         CASE(2)
           NSTAR = 2.0_q
         CASE(3)
           NSTAR = 3.0_q
         CASE(4)
           NSTAR = 3.7_q
         CASE(5)
           NSTAR = 4.0_q
         CASE(6)
           NSTAR = 4.2_q
         CASE DEFAULT
           call vtutor%bug('RADIAL_SLATER_FUNCTION: TYPE '//str(ITYP)//'does not exist',__FILE__,__LINE__)
      END SELECT
      FR(:)=(R%R(:))**(NSTAR-1.0_q) * EXP(-ZA*R%R(:)/NSTAR)

      ! Normalize
      FR=FR/SQRT(RADIAL_NORM2(FR,R))

      END SUBROUTINE RADIAL_SLATER_FUNCTION


!******************** SUBROUTINE RADIAL_GAUSSIAN_FUNCTION ***********************
!> Compute gaussian radial functions 
!********************************************************************************
      SUBROUTINE RADIAL_GAUSSIAN_FUNCTION(N,R,ZA,FR)
      USE radial
      USE constant
      USE base
      IMPLICIT NONE
      TYPE(rgrid) R
      INTEGER,INTENT(IN) :: N
      REAL(q) ZA
      REAL(q) FR(:)
      ! local
      INTEGER I

      DO I = 1, SIZE(FR)
        FR(I) = (R%R(I))**(1.0_q*N) * EXP(-ZA*(R%R(I))**2.0_q )
      ENDDO

      ! Normalize
      FR=FR/SQRT(RADIAL_NORM2(FR,R))

      END SUBROUTINE RADIAL_GAUSSIAN_FUNCTION

!******************** FUNCTION RADIAL_NORM **************************************
!> Compute norm^2 of a function defined on a logarithmic radial grid
!> int dr r^N * (r*f(r))**2
!> when N=0 the norm is obtained
!> when N=1 the first moment
!> when N=2 the second moment and so on...
!********************************************************************************
      FUNCTION RADIAL_NORM2(FR,R,N) RESULT(RES)
       use radial, only: rgrid
       REAL(q) :: FR(:) !< values of the function f(r) at the radial points
       TYPE(rgrid) :: R !< radial grid
       INTEGER,OPTIONAL :: N !< R power
       REAL(q) :: RES !< result
       ! local
       INTEGER :: I
       INTEGER :: MY_N
       MY_N=0; IF (PRESENT(N)) MY_N=N
       ! compute norm
       RES = 0.0_q
       DO I = 1, SIZE(FR)
         RES = RES + R%SI(I) * R%R(I)**MY_N * ABS(R%R(I)*FR(I))**2.0_q
       ENDDO
      END FUNCTION RADIAL_NORM2

!******************** FUNCTION RADIAL_CALC_MOMENT *******************************
!> Compute 1st 2nd and 3rd moment if the radial function
!********************************************************************************
      SUBROUTINE RADIAL_CALC_MOMENT(FR,R,R1,R2)
       USE radial, ONLY: rgrid
       USE constant, ONLY: AUTOA
       REAL(q) :: FR(:)
       TYPE(rgrid) :: R
       REAL(q) :: R1
       REAL(q) :: R2

       ! calc <r> (not to be confused with <vec r> which is zero)
       R1 = RADIAL_NORM2(FR,R,1)

       ! calc <r^2>
       R2 = RADIAL_NORM2(FR,R,2)

      END SUBROUTINE RADIAL_CALC_MOMENT

!******************** SUBROUTINE RADIAL_WRITE_INFO ******************************
!> Write information about radial function
!********************************************************************************
      SUBROUTINE RADIAL_WRITE_INFO(RADIAL_TYPE,N,ZA,R1,R2,IO)
       USE base, ONLY: in_struct
       CHARACTER(LEN=*) :: RADIAL_TYPE
       INTEGER,INTENT(IN) :: N
       REAL(q),INTENT(IN) :: ZA
       REAL(q),INTENT(IN) :: R1
       REAL(q),INTENT(IN) :: R2
       TYPE(in_struct) :: IO
       IF(IO%IU6 >= 0) THEN
         WRITE(IO%IU6,'(A,I2,A,F8.4)') " "//TRIM(RADIAL_TYPE)//" n=", N, " alpha=", ZA
         WRITE(IO%IU6,'(A,F8.4)') "   <r> =", R1               ! corresponds to average electron radius of the orbital
         WRITE(IO%IU6,'(A,F8.4)') "   sqrt(<r^2>) =", SQRT(R2) ! corresponds to radial spread of final orbital
         WRITE(IO%IU6,'(A,F8.4)') "   sqrt(<r^2>-<r>^2) =", SQRT(R2-R1**2.0_q)
         WRITE(IO%IU6,'(A)') " "
       ENDIF
      END SUBROUTINE RADIAL_WRITE_INFO


!******************** SUBROUTINE RGRID_PARAMETERS **********************
!> Get parameters for the logarithmic radial grid
!***********************************************************************
      SUBROUTINE RGRID_PARAMETERS(RSTART,REND,H)
       REAL(q) :: RSTART !< starting radius (should avoid zero)
       REAL(q) :: REND !< ending radius
       REAL(q) :: H !< step
       ! This used to be the values here but to describe
       ! n=2 and n=3 a larger radius might be required
       ! the value of ZA also plays a role
       !RSTART=0.0025_q;REND=10._q;H=0.025_q
       RSTART=0.00002_q;REND=200._q;H=0.025_q
      END SUBROUTINE RGRID_PARAMETERS

!******************** SUBROUTINE SETRGRID ******************************
!> Initializer for the radial grid structure
!***********************************************************************
      SUBROUTINE SETRGRID( &
     & RSTART,REND,H,R &
     &)
      USE radial
      TYPE(rgrid) R
      REAL(q) RSTART,REND,H
      ! local variables
      INTEGER I
      
      I=0
      DO 
        I=I+1
        IF (RSTART*EXP(H*(I-1))>=REND) EXIT
      ENDDO
      ! make sure I is uneven
      I=I+(MODULO(I,2)+1)
      R%NMAX=I
      R%RSTART=RSTART; R%H=H; R%D=EXP(H)
      IF (ASSOCIATED(R%R)) THEN
         DEALLOCATE(R%R); NULLIFY(R%R)
      ENDIF
      ALLOCATE(R%R(I))
      DO I=1,R%NMAX
         R%R(I)=R%RSTART*EXP(H*(I-1))
      ENDDO      
      R%REND=R%R(R%NMAX)

      IF (ASSOCIATED(R%SI)) THEN
         DEALLOCATE(R%SI); NULLIFY(R%SI)
      ENDIF
      CALL SET_SIMP(R)

      RETURN
      END SUBROUTINE SETRGRID


!******************** SUBROUTINE LPRJ_ORBITAL_DEFINITIONS **************
!
! VASP sets up the YLM functions in the following order
! (see for instance SETYLM in asa.F)
! YLM(:,1)     -> s
! YLM(:,2:4)   -> p:= y, z, x
! YLM(:,5:9)   -> d:= xy, yz, z2, xz, x2
! YLM(:,10:16) -> f:= y(3x2-y2), xyz, yz2, z3, xz2, z(x2-y2), x(x2-3y2)           
!
! This routine provides the aforementioned and adds the
! sp, sp2, sp3, sp3d and sp3d2 orbital definitions used in wannier90
! 
!***********************************************************************
      SUBROUTINE LPRJ_ORBITAL_DEFINITIONS( &
     &   L,M,HYBRID_ORBITAL &
     &)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INTEGER L,M
      REAL(q) HYBRID_ORBITAL(:)
      ! local variables
      
      ! we should be able to deal with anything up to and including L=3
      IF (SIZE(HYBRID_ORBITAL)<16) THEN
         CALL vtutor%error("LPRJ_ORBITAL_DEFINITIONS: ERROR: HYBRID_ORBITAL array too small " // &
            str(SIZE(HYBRID_ORBITAL)))
      ENDIF
      
      HYBRID_ORBITAL=0
      
      SELECT CASE(L)
         CASE(0)
            ! s-function
            IF (M==1) THEN
               HYBRID_ORBITAL(1)=1._q
            ELSE
               CALL vtutor%error("LPRJ_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // str(M) &
                  // " not implemented")
            ENDIF
         CASE(1)
            ! p-functions
            IF (M==1) THEN
               HYBRID_ORBITAL(2)=1._q
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(3)=1._q
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(4)=1._q
            ELSE
               CALL vtutor%error("LPRJ_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // str(M) &
                  // " not implemented")
            ENDIF
         CASE(2)
            ! d-functions
            IF (M==1) THEN
               HYBRID_ORBITAL(5)=1._q
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(6)=1._q
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(7)=1._q
            ELSEIF (M==4) THEN
               HYBRID_ORBITAL(8)=1._q
            ELSEIF (M==5) THEN
               HYBRID_ORBITAL(9)=1._q
            ELSE
               CALL vtutor%error("LPRJ_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // str(M) &
                  // " not implemented")
            ENDIF
         CASE(3)
            ! f-functions
            IF (M==1) THEN
               HYBRID_ORBITAL(10)=1._q
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(11)=1._q
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(12)=1._q
            ELSEIF (M==4) THEN
               HYBRID_ORBITAL(13)=1._q
            ELSEIF (M==5) THEN
               HYBRID_ORBITAL(14)=1._q
            ELSEIF (M==6) THEN
               HYBRID_ORBITAL(15)=1._q
            ELSEIF (M==7) THEN
               HYBRID_ORBITAL(16)=1._q
            ELSE
               CALL vtutor%error("LPRJ_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // str(M) &
                  // " not implemented")
            ENDIF
         CASE(4:)
            ! nothing beyond L=3 yet
            CALL vtutor%error("LPRJ_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // str(M) // &
               " not implemented")
         CASE(-1)
            ! sp-hybrids
            IF (M==1) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(2._q)
               HYBRID_ORBITAL(4)= 1._q/SQRT(2._q)
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(2._q)
               HYBRID_ORBITAL(4)=-1._q/SQRT(2._q)
            ELSE
               CALL vtutor%error("LPRJ_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // str(M) &
                  // " not implemented")
            ENDIF
         CASE(-2)
            ! sp2-hybrids
            IF (M==1) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(3._q)
               HYBRID_ORBITAL(4)=-1._q/SQRT(6._q)
               HYBRID_ORBITAL(2)= 1._q/SQRT(2._q)
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(3._q)
               HYBRID_ORBITAL(4)=-1._q/SQRT(6._q)
               HYBRID_ORBITAL(2)=-1._q/SQRT(2._q)
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(3._q)
               HYBRID_ORBITAL(4)= 2._q/SQRT(6._q)
            ELSE
               CALL vtutor%error("LPRJ_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // str(M) &
                  // " not implemented")
            ENDIF
         CASE(-3)
            ! sp3-hybrids
            IF (M==1) THEN
               HYBRID_ORBITAL(1)= 1._q/2._q
               HYBRID_ORBITAL(4)= 1._q/2._q
               HYBRID_ORBITAL(2)= 1._q/2._q
               HYBRID_ORBITAL(3)= 1._q/2._q
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(1)= 1._q/2._q
               HYBRID_ORBITAL(4)= 1._q/2._q
               HYBRID_ORBITAL(2)=-1._q/2._q
               HYBRID_ORBITAL(3)=-1._q/2._q
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(1)= 1._q/2._q
               HYBRID_ORBITAL(4)=-1._q/2._q
               HYBRID_ORBITAL(2)= 1._q/2._q
               HYBRID_ORBITAL(3)=-1._q/2._q
            ELSEIF (M==4) THEN
               HYBRID_ORBITAL(1)= 1._q/2._q
               HYBRID_ORBITAL(4)=-1._q/2._q
               HYBRID_ORBITAL(2)=-1._q/2._q
               HYBRID_ORBITAL(3)= 1._q/2._q
            ELSE
               CALL vtutor%error("LPRJ_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // str(M) &
                  // " not implemented")
            ENDIF
         CASE(-4)
            ! spd3d-hybrids
            IF (M==1) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(3._q)
               HYBRID_ORBITAL(4)=-1._q/SQRT(6._q)
               HYBRID_ORBITAL(2)= 1._q/SQRT(2._q)
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(3._q)
               HYBRID_ORBITAL(4)=-1._q/SQRT(6._q)
               HYBRID_ORBITAL(2)=-1._q/SQRT(2._q)
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(3._q)
               HYBRID_ORBITAL(4)= 2._q/SQRT(6._q)
            ELSEIF (M==4) THEN
               HYBRID_ORBITAL(3)= 1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)= 1._q/SQRT(2._q)
            ELSEIF (M==5) THEN
               HYBRID_ORBITAL(3)=-1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)= 1._q/SQRT(2._q)
            ELSE
               CALL vtutor%error("LPRJ_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // str(M) &
                  // " not implemented")
            ENDIF
         CASE(-5)
            ! sp3d2-hybrids
            IF (M==1) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(6._q)
               HYBRID_ORBITAL(4)=-1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)=-1._q/SQRT(12._q)
               HYBRID_ORBITAL(9)= 1._q/2._q
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(6._q)
               HYBRID_ORBITAL(4)= 1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)=-1._q/SQRT(12._q)
               HYBRID_ORBITAL(9)= 1._q/2._q
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(6._q)
               HYBRID_ORBITAL(2)=-1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)=-1._q/SQRT(12._q)
               HYBRID_ORBITAL(9)=-1._q/2._q
            ELSEIF (M==4) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(6._q)
               HYBRID_ORBITAL(2)= 1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)=-1._q/SQRT(12._q)
               HYBRID_ORBITAL(9)=-1._q/2._q
            ELSEIF (M==5) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(6._q)
               HYBRID_ORBITAL(3)=-1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)= 1._q/SQRT(3._q)
            ELSEIF (M==6) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(6._q)
               HYBRID_ORBITAL(3)= 1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)= 1._q/SQRT(3._q)
            ELSE
               CALL vtutor%error("LPRJ_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // str(M) &
                  // " not implemented")
            ENDIF
         CASE(:-6)
            CALL vtutor%error("LPRJ_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // str(M) // &
               " not implemented")
         
      END SELECT
      
      RETURN
      END SUBROUTINE LPRJ_ORBITAL_DEFINITIONS

!******************** SUBROUTINE WANNIER90_ORBITAL_DEFINITIONS *********
!
! VASP sets up the YLM functions in the following order
! (see for instance SETYLM in asa.F)
! YLM(:,1)     -> s
! YLM(:,2:4)   -> p:= y, z, x
! YLM(:,5:9)   -> d:= xy, yz, z2, xz, x2
! YLM(:,10:16) -> f:= y(3x2-y2), xyz, yz2, z3, xz2, z(x2-y2), x(x2-3y2)
!
! This routine provides a translation between the aforementioned and
! the orbital definitions used in wannier90
!
!***********************************************************************
      SUBROUTINE WANNIER90_ORBITAL_DEFINITIONS( &
     &   L,M,HYBRID_ORBITAL &
     &)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      INTEGER L,M
      REAL(q) HYBRID_ORBITAL(:)
      ! local variables

      ! we should be able to deal with anything up to and including L=3
      IF (SIZE(HYBRID_ORBITAL)<16) THEN
         CALL vtutor%error("WANNIER90_ORBITAL_DEFINITIONS: ERROR: HYBRID_ORBITAL array too small " // &
            str(SIZE(HYBRID_ORBITAL)))
      ENDIF

      HYBRID_ORBITAL=0

      SELECT CASE(L)
         CASE(0)
            ! s-function
            IF (M==1) THEN
               HYBRID_ORBITAL(1)=1._q
            ELSE
               CALL vtutor%error("WANNIER90_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // &
                  str(M) // " not implemented")
            ENDIF
         CASE(1)
            ! p-functions
            IF (M==1) THEN
               HYBRID_ORBITAL(3)=1._q
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(4)=1._q
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(2)=1._q
            ELSE
               CALL vtutor%error("WANNIER90_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // &
                  str(M) // " not implemented")
            ENDIF
         CASE(2)
            ! d-functions
            IF (M==1) THEN
               HYBRID_ORBITAL(7)=1._q
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(8)=1._q
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(6)=1._q
            ELSEIF (M==4) THEN
               HYBRID_ORBITAL(9)=1._q
            ELSEIF (M==5) THEN
               HYBRID_ORBITAL(5)=1._q
            ELSE
               CALL vtutor%error("WANNIER90_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // &
                  str(M) // " not implemented")
            ENDIF
         CASE(3)
            ! f-functions
            IF (M==1) THEN
               HYBRID_ORBITAL(13)=1._q
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(14)=1._q
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(12)=1._q
            ELSEIF (M==4) THEN
               HYBRID_ORBITAL(15)=1._q
            ELSEIF (M==5) THEN
               HYBRID_ORBITAL(11)=1._q
            ELSEIF (M==6) THEN
               HYBRID_ORBITAL(16)=1._q
            ELSEIF (M==7) THEN
               HYBRID_ORBITAL(10)=1._q
            ELSE
               CALL vtutor%error("WANNIER90_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // &
                  str(M) // " not implemented")
            ENDIF
         CASE(4:)
            ! nothing beyond L=3 yet
            CALL vtutor%error("WANNIER90_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // &
               str(M) // " not implemented")
         CASE(-1)
            ! sp-hybrids
            IF (M==1) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(2._q)
               HYBRID_ORBITAL(4)= 1._q/SQRT(2._q)
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(2._q)
               HYBRID_ORBITAL(4)=-1._q/SQRT(2._q)
            ELSE
               CALL vtutor%error("WANNIER90_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // &
                  str(M) // " not implemented")
            ENDIF
         CASE(-2)
            ! sp2-hybrids
            IF (M==1) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(3._q)
               HYBRID_ORBITAL(4)=-1._q/SQRT(6._q)
               HYBRID_ORBITAL(2)= 1._q/SQRT(2._q)
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(3._q)
               HYBRID_ORBITAL(4)=-1._q/SQRT(6._q)
               HYBRID_ORBITAL(2)=-1._q/SQRT(2._q)
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(3._q)
               HYBRID_ORBITAL(4)= 2._q/SQRT(6._q)
            ELSE
               CALL vtutor%error("WANNIER90_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // &
                  str(M) // " not implemented")
            ENDIF
         CASE(-3)
            ! sp3-hybrids
            IF (M==1) THEN
               HYBRID_ORBITAL(1)= 1._q/2._q
               HYBRID_ORBITAL(4)= 1._q/2._q
               HYBRID_ORBITAL(2)= 1._q/2._q
               HYBRID_ORBITAL(3)= 1._q/2._q
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(1)= 1._q/2._q
               HYBRID_ORBITAL(4)= 1._q/2._q
               HYBRID_ORBITAL(2)=-1._q/2._q
               HYBRID_ORBITAL(3)=-1._q/2._q
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(1)= 1._q/2._q
               HYBRID_ORBITAL(4)=-1._q/2._q
               HYBRID_ORBITAL(2)= 1._q/2._q
               HYBRID_ORBITAL(3)=-1._q/2._q
            ELSEIF (M==4) THEN
               HYBRID_ORBITAL(1)= 1._q/2._q
               HYBRID_ORBITAL(4)=-1._q/2._q
               HYBRID_ORBITAL(2)=-1._q/2._q
               HYBRID_ORBITAL(3)= 1._q/2._q
            ELSE
               CALL vtutor%error("WANNIER90_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // &
                  str(M) // " not implemented")
            ENDIF
         CASE(-4)
            ! spd3d-hybrids
            IF (M==1) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(3._q)
               HYBRID_ORBITAL(4)=-1._q/SQRT(6._q)
               HYBRID_ORBITAL(2)= 1._q/SQRT(2._q)
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(3._q)
               HYBRID_ORBITAL(4)=-1._q/SQRT(6._q)
               HYBRID_ORBITAL(2)=-1._q/SQRT(2._q)
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(3._q)
               HYBRID_ORBITAL(4)= 2._q/SQRT(6._q)
            ELSEIF (M==4) THEN
               HYBRID_ORBITAL(3)= 1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)= 1._q/SQRT(2._q)
            ELSEIF (M==5) THEN
               HYBRID_ORBITAL(3)=-1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)= 1._q/SQRT(2._q)
            ELSE
               CALL vtutor%error("WANNIER90_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // &
                  str(M) // " not implemented")
            ENDIF
         CASE(-5)
            ! sp3d2-hybrids
            IF (M==1) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(6._q)
               HYBRID_ORBITAL(4)=-1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)=-1._q/SQRT(12._q)
               HYBRID_ORBITAL(9)= 1._q/2._q
            ELSEIF (M==2) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(6._q)
               HYBRID_ORBITAL(4)= 1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)=-1._q/SQRT(12._q)
               HYBRID_ORBITAL(9)= 1._q/2._q
            ELSEIF (M==3) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(6._q)
               HYBRID_ORBITAL(2)=-1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)=-1._q/SQRT(12._q)
               HYBRID_ORBITAL(9)=-1._q/2._q
            ELSEIF (M==4) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(6._q)
               HYBRID_ORBITAL(2)= 1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)=-1._q/SQRT(12._q)
               HYBRID_ORBITAL(9)=-1._q/2._q
            ELSEIF (M==5) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(6._q)
               HYBRID_ORBITAL(3)=-1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)= 1._q/SQRT(3._q)
            ELSEIF (M==6) THEN
               HYBRID_ORBITAL(1)= 1._q/SQRT(6._q)
               HYBRID_ORBITAL(3)= 1._q/SQRT(2._q)
               HYBRID_ORBITAL(7)= 1._q/SQRT(3._q)
            ELSE
               CALL vtutor%error("WANNIER90_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // &
                  str(M) // " not implemented")
            ENDIF
         CASE(:-6)
            CALL vtutor%error("WANNIER90_ORBITAL_DEFINITIONS: ERROR: L= " // str(L) // " M= " // &
               str(M) // " not implemented")

      END SELECT

      RETURN
      END SUBROUTINE WANNIER90_ORBITAL_DEFINITIONS

      PURE FUNCTION LPRJ_IS_INITIALIZED()
         IMPLICIT NONE
         LOGICAL  :: LPRJ_IS_INITIALIZED
         LPRJ_IS_INITIALIZED = IS_INITIALIZED
      ENDFUNCTION LPRJ_IS_INITIALIZED

      PURE FUNCTION LPRJ_IS_CALLED_PROALL()
         IMPLICIT NONE
         LOGICAL  :: LPRJ_IS_CALLED_PROALL
         LPRJ_IS_CALLED_PROALL = IS_CALLED_PROALL
      ENDFUNCTION LPRJ_IS_CALLED_PROALL

      PURE FUNCTION LPRJ_GET_NUM_WAN() RESULT(NWAN)
         IMPLICIT NONE
         INTEGER  :: NWAN
         IF ( ALLOCATED( LPRJ_functions ) ) THEN
            NWAN = SIZE(LPRJ_functions, 1)
         ELSE
            NWAN = 0 
         ENDIF
      ENDFUNCTION LPRJ_GET_NUM_WAN

      PURE FUNCTION LPRJ_NUM_ORBITALS_ON_ATOM(ATOM_IDX) RESULT(NUM_ORBS)
         IMPLICIT NONE
         INTEGER  :: NUM_ORBS
         INTEGER, INTENT(IN)  :: ATOM_IDX

         IF (ALLOCATED(LPRJ_functions)) THEN
            NUM_ORBS = COUNT(LPRJ_functions(:)%poscar_site == ATOM_IDX)
         ELSE
            NUM_ORBS = 0
         ENDIF
      END FUNCTION LPRJ_NUM_ORBITALS_ON_ATOM

      PURE FUNCTION LPRJ_GET_POS(WAN_IDX) RESULT(POS)
         IMPLICIT NONE

         INTEGER, INTENT(IN)  :: WAN_IDX
         REAL(q)              :: POS(3)

         POS(:) = LPRJ_functions(WAN_IDX)%R(:)
      ENDFUNCTION LPRJ_GET_POS

      PURE FUNCTION LPRJ_GET_SITE(WAN_IDX) RESULT(NSITE)
         IMPLICIT NONE

         INTEGER, INTENT(IN)  :: WAN_IDX
         INTEGER              :: NSITE

         NSITE = LPRJ_functions(WAN_IDX)%poscar_site
      ENDFUNCTION LPRJ_GET_SITE

!******************** SUBROUTINE SETROTYLM *****************************
!> @brief Setup the Ylm rotation matrix in accordance with XIN and ZIN axis
!***********************************************************************
      SUBROUTINE SETROTYLM( XIN,ZIN,LMAX,ROTYLM )
         USE asa
         USE constant
         USE string, ONLY: str
         USE tutor, ONLY: vtutor
         INTEGER,INTENT(IN) :: LMAX !< Maximum value of l chanel
         REAL(q),INTENT(IN) :: XIN(3) !< xaxis defining the rotation
         REAL(q),INTENT(IN) :: ZIN(3) !< yaxis defining the rotation
         GDEF,INTENT(OUT) :: ROTYLM(:,:) !< Ylm rotation matrix
         ! local variables
         INTEGER LMMAX
         REAL(q) XN,ZN,XDOTZ
         REAL(q) X(3),Y(3),Z(3),U(3,3)
         REAL(q), PARAMETER :: TINY=1.0E-6_q

         INTEGER I,J,IFAIL,LM,LMP
         INTEGER PHPTS,THPTS,NPTS,NP
         REAL(q) SCALE,DELTAPHI,SIM_FAKT
         REAL(q), ALLOCATABLE :: RADPTS(:,:),XYZPTS(:,:),UXYZPTS(:,:)
         REAL(q), ALLOCATABLE :: YLM(:,:),YLMP(:,:)
         REAL(q), ALLOCATABLE :: WEIGHT(:),ABSCIS(:)
         EXTERNAL GAUSSI2

         X=XIN; Z=ZIN

         LMMAX=(LMAX+1)**2

         ! check size of ROTYLM
         IF (SIZE(ROTYLM,1)<LMMAX.OR.SIZE(ROTYLM,2)<LMMAX) THEN
            CALL vtutor%error("SETROTYLM: ERROR: ROTYLM too small: " // str(LMMAX) // " " // &
               str(SIZE(ROTYLM,1)) // " " // str(SIZE(ROTYLM,2)))
         ENDIF
         ! check size of X and Z
         XN=SQRT(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
         ZN=SQRT(Z(1)*Z(1)+Z(2)*Z(2)+Z(3)*Z(3))
         IF (XN<TINY.OR.ZN<TINY) THEN
            CALL vtutor%error("SETROTYLM: ERROR: |X| or |Z| very small: " // str(XN) // " " // str(ZN))
         ENDIF
         ! and normalize to 1
         X=X/XN; Z=Z/ZN
         ! check orthogonality of X and Z
         XDOTZ=X(1)*Z(1)+X(2)*Z(2)+X(3)*Z(3)
         IF (ABS(XDOTZ)>TINY) THEN
            CALL vtutor%error("SETROTYLM: ERROR: X and Z are not orthogonal (enough): " // str(ABS(XDOTZ)))
         ENDIF
         ! y=Z \times X
         Y(1)=(Z(2)*X(3)-X(2)*Z(3))
         Y(2)=(Z(3)*X(1)-X(3)*Z(1))
         Y(3)=(Z(1)*X(2)-X(1)*Z(2))
         ! transformation matrix
         U(1,:)=X(:); U(2,:)=Y(:); U(3,:)=Z(:)
         !test
         !  WRITE(*,*) 'transformation matrix'
         !  WRITE(*,'(3F10.5)') U(1,:)
         !  WRITE(*,'(3F10.5)') U(2,:)
         !  WRITE(*,'(3F10.5)') U(3,:)
         !test
         SCALE=2*SQRT(PI) ! 1/Y00

         !========================================================================
         ! number of theta and phi pivot points to perform angular integration
         ! since Exc=f(a*Yllmax,m) we need more pivot points than theoretically
         ! needed to integrate Yllmax,m.
         ! the factor 2 is the minium, 3 is more accurate
         !========================================================================
         PHPTS=3*(LMAX+1)
         THPTS=3*FLOOR(REAL(LMAX/2+1,KIND=q))
         NPTS=PHPTS*THPTS
         DELTAPHI=REAL(2_q*PI/PHPTS,KIND=q)
         ! allocate arrays
         ALLOCATE(XYZPTS(NPTS,3),UXYZPTS(NPTS,3),RADPTS(NPTS,2),WEIGHT(THPTS),ABSCIS(THPTS))
         ALLOCATE(YLM(NPTS,LMMAX),YLMP(NPTS,LMMAX))

         ! set phi positions, equally spaced
         RADPTS=0; WEIGHT=0; ABSCIS=0
         DO I=1,PHPTS
            DO J=1,THPTS
               RADPTS((J-1)*PHPTS+I,2)=(I-1)*DELTAPHI
            ENDDO
         ENDDO
         ! get theta positions (actually get cos(theta)) (Gauss integration)
         CALL GAUSSI(GAUSSI2,-1._q,1._q,0,THPTS,WEIGHT,ABSCIS,IFAIL)
         DO I=1,THPTS
            RADPTS((I-1)*PHPTS+1:I*PHPTS,1)=ABSCIS(I)
         ENDDO
         ! convert radial to cartesian coordinates
         DO I=1,NPTS
            XYZPTS(I,1)=COS(RADPTS(I,2))*SQRT(1_q-RADPTS(I,1)**2_q) ! x
            XYZPTS(I,2)=SIN(RADPTS(I,2))*SQRT(1_q-RADPTS(I,1)**2_q) ! y
            XYZPTS(I,3)=RADPTS(I,1)                                 ! z
         ENDDO

         YLM=0
         CALL SETYLM(LMAX,NPTS,YLM,XYZPTS(:,1),XYZPTS(:,2),XYZPTS(:,3))

         ! rotate XYZPTS, using accordance with U
         DO I=1,NPTS
            UXYZPTS(I,1)=U(1,1)*XYZPTS(I,1)+U(1,2)*XYZPTS(I,2)+U(1,3)*XYZPTS(I,3)
            UXYZPTS(I,2)=U(2,1)*XYZPTS(I,1)+U(2,2)*XYZPTS(I,2)+U(2,3)*XYZPTS(I,3)
            UXYZPTS(I,3)=U(3,1)*XYZPTS(I,1)+U(3,2)*XYZPTS(I,2)+U(3,3)*XYZPTS(I,3)
         ENDDO

         YLMP=0
         CALL SETYLM(LMAX,NPTS,YLMP,UXYZPTS(:,1),UXYZPTS(:,2),UXYZPTS(:,3))

         ROTYLM=0
         ! loop over all points in the angular grid
         points: DO NP=1,NPTS

            ! weight of this points
            SIM_FAKT=DELTAPHI*WEIGHT((INT((NP-1)/PHPTS)+1))

            DO LM=1,LMMAX
            DO LMP=1,LMMAX
               ROTYLM(LMP,LM)=ROTYLM(LMP,LM)+YLMP(NP,LMP)*YLM(NP,LM)*SIM_FAKT
            ENDDO
            ENDDO

         ENDDO points

         DEALLOCATE(XYZPTS,UXYZPTS,RADPTS,WEIGHT,ABSCIS,YLM,YLMP)
         RETURN
      END SUBROUTINE SETROTYLM

      END MODULE locproj
