#include "symbol.inc"
       MODULE dynconstr
        USE random_seeded
        USE prec
        USE constant
        USE poscar
        USE base
        USE lattice
        USE ini
        USE chain
        USE vaspxml
        USE mymath
        USE internals
        USE wave
        USE main_mpi
        USE npt_dynamics
        IMPLICIT NONE 
        !INCLUDE "hills.inc"

        TYPE hills_io
          INTEGER :: STRUCTINPUT = 0
          INTEGER :: PENALTY = 0
        END TYPE hills_io

        TYPE penalty_data
          INTEGER         :: number = 0
          TYPE(gauss_peak),POINTER :: gauss(:) => NULL()
          REAL(q),POINTER :: force(:) => NULL()
          REAL(q),POINTER :: wall(:,:) => NULL()
        END TYPE penalty_data

        TYPE gauss_peak
          REAL(q),POINTER :: position(:) => NULL()
          REAL(q)         :: high = 0.0_q
          REAL(q)         :: width = 0.0_q
        END TYPE gauss_peak

        TYPE hills_data
          TYPE(gauss_peak),POINTER :: gauss(:) => NULL()
          REAL(q),POINTER :: velocity(:) => NULL()
          REAL(q),POINTER :: position(:) => NULL()
          REAL(q),POINTER :: mass(:) => NULL()
          REAL(q),POINTER :: force_constant(:) => NULL()
          REAL(q),POINTER :: force(:) => NULL()
          REAL(q)         :: SNOSE(4) = 0.0_q
          REAL(q)         :: SQQ = 0.0_q
          REAL(q)         :: potential = 0.0_q
          REAL(q)         :: stride = 0.0_q
          REAL(q)         :: andersen_prob = 0.0_q
          REAL(q)         :: temperature = 0.0_q !temperature of fict. particles
          INTEGER         :: number = 0
          INTEGER         :: bin = 0
          INTEGER         :: maxstride = 0
          LOGICAL         :: variable_width = .FALSE.
        END TYPE

        !c initialization done or not
        LOGICAL :: LINIT_GL = .FALSE.

        !c parameters read in from           INCAR
        INTEGER :: MDALGO_GL=0               ! MDALGO
        REAL(q) :: ANDERSEN_PROB_GL=0._q     ! ANDERSEN_PROB 
        INTEGER :: SCALING_GL=0              ! SCALING
        INTEGER :: NSUBSYS_GL(3)             ! NSUBSYS
        REAL(q) :: TSUBSYS_GL(3)             ! TSUBSYS
        REAL(q) :: PSUBSYS_GL(3)             ! PSUBSYS
        REAL(q),ALLOCATABLE :: GAMMA_GL(:)   ! LANGEVIN_GAMMA
        REAL(q) :: GAMMA_L_GL=0._q           ! LANGEVIN_GAMMA_L
        REAL(q) :: AMASS_GL=1000._q          ! PMASS
        LOGICAL :: LATTICE_CONSTRAINTS_GL(3)=(/.TRUE.,.TRUE.,.TRUE./)  ! LATTICE_CONSTRAINTS
        REAL(q) :: CNEXP_GL(2)=(/9.,14./)    ! CNEXP
        REAL(q) :: TILAMBDA_GL=-1.           ! TILAMBDA
        REAL(q) :: TIDOT_GL=0._q             ! TIDOT
        REAL(q),ALLOCATABLE :: valueA_GL(:)  ! VALUE_MIN
        REAL(q),ALLOCATABLE :: valueB_GL(:)  ! VALUE_MAX 
        REAL(q),ALLOCATABLE :: INCREM_GL(:)  ! INCREM
        INTEGER :: EQUI_REGIME_GL=0          ! EQUI_REGIME_GL
        LOGICAL :: LBLUEOUT_GL=.FALSE.       ! LBLUEOUT 
        REAL(q) :: BMTOL_GL=1e-5             ! SHAKETOL
        REAL(q) :: BMTOLSOFT_GL=1e-5         ! SHAKETOLSOFT
        REAL(q) :: BMSCA_GL=2._q             ! SHAKESCA
        INTEGER :: BMITER_GL=1000            ! SHAKEITER
                                             ! HILLS_BIN: hills_GL%bin
                                             ! HILLS_W: hills_GL%gauss(1)%width
                                             ! HILLS_H: hills_GL%gauss(1)%high


        !c module's internal variables
        LOGICAL :: Lhspring_V0_GL=.FALSE.
        LOGICAL :: LTIDOT_GL=.FALSE. !c is the TI speed set up by a user?
        REAL(q), ALLOCATABLE :: hspring_R0_GL(:),hspring_V0_GL(:),hspring_K_GL(:),hspring_force_GL(:)
        REAL(q), ALLOCATABLE :: fbias_R0_GL(:),fbias_D_GL(:),fbias_A_GL(:)
        !REAL(q), ALLOCATABLE :: spring_accel(:,:)
        REAL(q), ALLOCATABLE :: HARMONIC_HESSIAN_GL(:,:) ,POSITIONS0_GL(:,:)
        REAL(q) ::  HARMONIC_ENERGY0_GL
        LOGICAL :: LMIN_GL=.TRUE.
        LOGICAL :: LMAX_GL=.TRUE.
        TYPE(coordstructure) :: ICOORDINATES_GL
        TYPE(hills_data) :: hills_GL
        TYPE(penalty_data) :: penalty_GL
        INTEGER :: LATTICE_DOF_GL=0


         
        CONTAINS

        SUBROUTINE DYNCONSTR_INIT(DYN,T_INFO,INFO,LATT_CUR, IO,g_io)
        !c initialization for all MD methods of this module
        !c the input parameters read in from INCAR are currently writtent 
        !c into the REPORT file within the following readers: 
        !c DYN_READER, READER_ICONST0, READER_ICONST7, HILLS_READER
        !c PENALTY_READER, FBIAS_READER, HARMONIC_SPRING_READER 
          TYPE(dynamics) :: DYN
          TYPE(type_info) :: T_INFO
          TYPE (info_struct) :: INFO
          TYPE(latt) :: LATT_CUR
          TYPE (in_struct) :: IO
          TYPE(gadget_io) :: g_io
          INTEGER :: i

          IF (LINIT_GL) RETURN

          NSUBSYS_GL=T_INFO%NIONS
          TSUBSYS_GL=0._q
          PSUBSYS_GL=0._q
          ALLOCATE(GAMMA_GL(T_INFO%NTYP))
          GAMMA_GL=0._q
          CALL DYN_READER(IO,g_io,MDALGO_GL,ANDERSEN_PROB_GL,SCALING_GL,NSUBSYS_GL,TSUBSYS_GL,PSUBSYS_GL,GAMMA_GL, &
          &    GAMMA_L_GL,T_INFO%NIONS,T_INFO%NTYP,DYN%ISIF,CNEXP_GL,&
          &    TILAMBDA_GL,TIDOT_GL,LTIDOT_GL,AMASS_GL,LATTICE_CONSTRAINTS_GL,LATTICE_DOF_GL)

          !c read in harmonic forcefield for a TI calculation
          IF (TILAMBDA_GL .GE. 0.) THEN
            ALLOCATE(HARMONIC_HESSIAN_GL(3*T_INFO%NIONS,3*T_INFO%NIONS),POSITIONS0_GL(3,T_INFO%NIONS))
            CALL READ_HARMONIC_FORCEFIELD(IO,g_io,DYN,T_INFO,HARMONIC_HESSIAN_GL,&
            &  HARMONIC_ENERGY0_GL,POSITIONS0_GL)
          ENDIF

          nullify(ICOORDINATES_GL%COORDSTRUCT)
          ALLOCATE(ICOORDINATES_GL%COORDSTRUCT(10))
          ICOORDINATES_GL%NUMINTERNALS=0
          ICOORDINATES_GL%CNEXP=CNEXP_GL

          CALL constraints_reader(T_INFO,g_io,DYN%POSION,TRANSPOSE(LATT_CUR%A),LATT_CUR%B,ICOORDINATES_GL,IO)

          !c transformation rate for the STATUS==0 coordinates
          IF (ICOORDINATES_GL%iconst0>0) THEN
            ALLOCATE(INCREM_GL(ICOORDINATES_GL%iconst0))
            INCREM_GL=0._q
            CALL READER_ICONST0(IO,g_io,ICOORDINATES_GL%iconst0,BMTOL_GL,BMTOLSOFT_GL,BMSCA_GL, &
            &    BMITER_GL,LBLUEOUT_GL,EQUI_REGIME_GL,INCREM_GL)           
          END IF
          

          !c limits for the monitored coordinates
          IF (ICOORDINATES_GL%iconst7>0) THEN
            ALLOCATE(valueA_GL(ICOORDINATES_GL%iconst7),valueB_GL(ICOORDINATES_GL%iconst7))
            valueA_GL=-1000._q
            valueB_GL= 1000._q
            CALL READER_ICONST7(IO,g_io,ICOORDINATES_GL,LMIN_GL,LMAX_GL,valueA_GL,valueB_GL)
          ENDIF


          !c reference positions, their velocities, and force-constants for spring connections, internal coordinates
          IF (ICOORDINATES_GL%iconst4>0) THEN
            ALLOCATE(fbias_R0_GL(ICOORDINATES_GL%iconst4),fbias_D_GL(ICOORDINATES_GL%iconst4),fbias_A_GL(ICOORDINATES_GL%iconst4))
            fbias_R0_GL=0._q;fbias_D_GL=0._q;fbias_A_GL=0._q
            CALL FBIAS_READER(IO,ICOORDINATES_GL%iconst4,g_io,fbias_R0_GL,fbias_D_GL,fbias_A_GL)
          ELSE IF(ICOORDINATES_GL%iconst4.EQ.0) THEN
            ALLOCATE(fbias_R0_GL(1:1), fbias_D_GL(1:1), fbias_A_GL(1:1)) 
          ENDIF
              
          !c reference positions, their velocities, and force-constants for spring connections, internal coordinates
          IF (ICOORDINATES_GL%iconst8>0 ) THEN
            ALLOCATE(hspring_R0_GL(ICOORDINATES_GL%iconst8),hspring_V0_GL(ICOORDINATES_GL%iconst8), &
            &  hspring_K_GL(ICOORDINATES_GL%iconst8),hspring_force_GL(ICOORDINATES_GL%iconst8))
            hspring_R0_GL=0._q;hspring_K_GL=0._q;hspring_force_GL=0._q;hspring_V0_GL=0._q
            CALL HARMONIC_SPRING_READER(IO,ICOORDINATES_GL%iconst8,g_io,hspring_R0_GL,hspring_V0_GL,hspring_K_GL,Lhspring_V0_GL)
          ELSE IF(ICOORDINATES_GL%iconst8.EQ.0) THEN
            ALLOCATE(hspring_R0_GL(1:1),hspring_V0_GL(1:1), &
            &  hspring_K_GL(1:1),hspring_force_GL(1:1))
            hspring_R0_GL=0._q;hspring_K_GL=0._q;hspring_force_GL=0._q
          ENDIF
                                                   
          !c mMD stuff
          hills_GL%number=ICOORDINATES_GL%iconst5+ICOORDINATES_GL%iconst6
   
          IF (hills_GL%number>0) THEN
            IF (MDALGO_GL==10 .OR. MDALGO_GL==1 .OR. MDALGO_GL==20 .OR. MDALGO_GL==2 .OR. MDALGO_GL==3) THEN 
              IF (hills_GL%number<1 .AND. (MDALGO_GL==10 .OR. MDALGO_GL==20)) THEN
                CALL vtutor%error("Error meta-MD: at least one coord. with STATUS=5 must be defined &
                   &if MDALGO=10|20")
              ENDIF
              ALLOCATE(hills_GL%mass(hills_GL%number),hills_GL%force_constant(hills_GL%number))
              ALLOCATE(hills_GL%velocity(hills_GL%number),hills_GL%position(hills_GL%number))
              ALLOCATE(hills_GL%force(hills_GL%number))
              hills_GL%mass=0._q; hills_GL%force_constant=0._q
              hills_GL%velocity=0._q; hills_GL%position=0._q
              hills_GL%potential=0._q; hills_GL%force=0._q

              CALL HILLS_READER(IO,DYN,hills_GL,g_io,ICOORDINATES_GL,MDALGO_GL)
              CALL PENALTY_READER(hills_GL,penalty_GL,IO,g_io)
            ENDIF
          ENDIF 


          LINIT_GL=.TRUE.

        END SUBROUTINE DYNCONSTR_INIT

        SUBROUTINE STEP_tb(DYN,T_INFO,INFO,LATT_CUR, &
          &      EKIN,EKIN_LAT,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,IO,EPOT,g_io,WDES,SIF,ISCALE,TEIN,TIFOR)
        !c main driver for different MD algorithms
          TYPE(dynamics) :: DYN
          TYPE(type_info) :: T_INFO
          TYPE (info_struct) :: INFO
          TYPE(latt) :: LATT_CUR
          TYPE (in_struct) :: IO
          TYPE(gadget_io) :: g_io          
          TYPE (wavedes)  ::   WDES   
          REAL(q) :: EKIN,EKIN_Lat,EPS,ES,DISMAX,EPOT,TEMPER,TEIN
          INTEGER :: i,j,k,ii
          INTEGER :: NDEGREES_OF_FREEDOM,cDOF
          REAL(q),SAVE,ALLOCATABLE :: REF_C(:,:)
          INTEGER,SAVE ::NDEGREES_OF_FREEDOM_
          INTEGER,SAVE :: counter=0
          REAL(q) :: SIF(3,3)
          INTEGER :: ISCALE
          INTEGER IDUM
          LOGICAL ::LDUM
          REAL(q) :: RDUM
          COMPLEX(q) :: CDUM
          CHARACTER*1 :: CHARAC
          REAL(q) :: TIFOR(3,T_INFO%NIONS)
          INTEGER :: Ltxyz(3)
          LOGICAL,SAVE :: LSCALE=.TRUE. !c rescale initial velocities so as to have exact Tinst
          LOGICAL, SAVE :: LCMASS=.TRUE. !c remove the center-of-mass motion?
          REAL(q) :: CMASS(3,1)
          REAL(q) :: ESPRING
          INTEGER :: ios
          REAL(q) :: HARMONIC_ENERGY
        

          counter=counter+1
          !CALL XML_VEL(T_INFO%NIONS,DYN%VEL)

          !DYN%POSIOC=DYN%POSION

          TEIN=0._q;ESPRING=0._q

          IF (counter==1) THEN
            !c make sure the initialization has been done  
            CALL DYNCONSTR_INIT(DYN,T_INFO,INFO,LATT_CUR, IO,g_io)
            
            IF (MDALGO_GL==3 ) LCMASS=.FALSE.

            !c  georg's fix 
            IF (MDALGO_GL>0) THEN
            ! copy current positions considered in VASP to POSIOC
               DYN%POSIOC=DYN%POSION
            ENDIF
                            
            NDEGREES_OF_FREEDOM_=NDEGREES_OF_FREEDOM

            !c Langevin dynamics does not fix total momentum
            !c hence NDEGREES_OF_FREEDOM=3*T_INFO%NIONS, not 3*T_INFO%NIONS-3
            IF (MDALGO_GL==3) THEN
              Ltxyz=1
              DO i=1,T_INFO%NIONS
                DO j=1,3
                  IF (.NOT. T_INFO%LSFOR(j,i)) THEN
                    Ltxyz(j)=0
                  ENDIF
                ENDDO
              ENDDO
              NDEGREES_OF_FREEDOM_=NDEGREES_OF_FREEDOM+SUM(Ltxyz)
            ENDIF

            !c write out some outpout
            IF (MDALGO_GL>0) THEN
          
              !c mMD stuff
              IF (MDALGO_GL==10 .OR. MDALGO_GL==1 .OR. MDALGO_GL==20 .OR. MDALGO_GL==2 .OR. MDALGO_GL==3) THEN 
                OPEN(UNIT=g_io%STRUCTINPUT,FILE='HILLSPOT',STATUS='UNKNOWN')
                IF (g_io%STRUCTINPUT>=0) THEN
                  DO i=1,penalty_GL%number
                    DO j=1,hills_GL%number
                      write(g_io%STRUCTINPUT,ADVANCE='NO',FMT='(X,F9.5)') penalty_GL%gauss(i)%position(j)
                    ENDDO
                    write(g_io%STRUCTINPUT,ADVANCE='NO',FMT='(X,F9.5)') penalty_GL%gauss(i)%high
                    write(g_io%STRUCTINPUT,ADVANCE='YES',FMT='(X,F9.5)') penalty_GL%gauss(i)%width
                  ENDDO
                ENDIF
                CLOSE(g_io%STRUCTINPUT)
              ENDIF

              cDOF=0
              IF (DYN%ISIF==3) THEN
                CALL CONSTRAINED_DOF(DYN,T_INFO,3*T_INFO%NIONS+9,LATT_CUR,ICOORDINATES_GL,cDOF)
              ELSE
                CALL CONSTRAINED_DOF(DYN,T_INFO,3*T_INFO%NIONS,LATT_CUR,ICOORDINATES_GL,cDOF)
              ENDIF
            
              NDEGREES_OF_FREEDOM=NDEGREES_OF_FREEDOM_-cDOF
          
              IF (IO%IU6>0) THEN
                write(g_io%REPORT,FMT='(/,3X,A35,X,I5)')    'original number of atomic DOF:     ', NDEGREES_OF_FREEDOM_
                write(g_io%REPORT,FMT='(3X,A35,X,I5)')    'number of independent constraints: ', cDOF
                write(g_io%REPORT,FMT='(3X,A35,X,I5)')    'number of active DOF:              ', NDEGREES_OF_FREEDOM
              ENDIF
            END IF
            !c write summaray of flags in vasprun.xml
            !CALL DYN_XML_WRITE(MDALGO_GL,ANDERSEN_PROB_GL,SCALING_GL,NSUBSYS_GL,TSUBSYS_GL,PSUBSYS_GL)
          ELSE 
            !c Georg's fix
            IF (MDALGO_GL>0) THEN
             ! copy current positions considered in VASP to POSIOC
                DYN%POSIOC=DYN%POSION
            ENDIF

          !CALL XML_TAG("stepmd")
          !CALL XML_INCAR('stepnum','I',counter,RDUM,CDUM,LDUM,CHARAC,1)

          ENDIF


          IF (MDALGO_GL>0) THEN
            IF (IO%IU6>0) THEN
              write(g_io%REPORT,FMT='(/,A40)') '========================================'
              write(g_io%REPORT,FMT='(A20,X,I7)') 'MD step No.',counter
              write(g_io%REPORT,FMT='(A40)') '========================================'
            ENDIF

            IF (DYN%INIT==-1) THEN
            !IF (DYN%INIT<1) THEN
              LSCALE=.TRUE.
              !IF (MDALGO_GL==3) LCMASS=.FALSE.
              CALL init_velocities(DYN%VEL,T_INFO,DYN,LATT_CUR,NDEGREES_OF_FREEDOM_,LSCALE,LCMASS)                
              IF (IO%IU6>0) write(g_io%REPORT,FMT='(/,2X,A40)') "Atomic velocities initialized by STEP_tb"
              DYN%INIT=0
            ENDIF
            
            !c define the reference coordinates needed
            !c to avoid problems with definition of internal parameters
            IF (counter==1) THEN
              ALLOCATE(REF_C(3,T_INFO%NIONS))
              REF_C=0._q
            ELSE
              CALL minimize_difference(DYN%POSIOC,REF_C,T_INFO%NIONS)
            ENDIF

            !c update position of minimum for harmonic spring
            IF ((ICOORDINATES_GL%iconst8>0) .AND. Lhspring_V0_GL) THEN
              hspring_R0_GL(:)=hspring_R0_GL(:)+DYN%POTIM*hspring_V0_GL(:)
              IF (IO%IU6>0) THEN
                write(g_io%REPORT,ADVANCE='NO',FMT='(3X,A22)') '          SPRING_R0 = '
                DO i=1,(ICOORDINATES_GL%iconst8)
                  write(g_io%REPORT,ADVANCE='NO',FMT='(X,F10.5)') hspring_R0_GL(i)
                ENDDO
              ENDIF
            ENDIF

            REF_C=DYN%POSIOC        
            CALL minimize_difference(DYN%POSION,DYN%POSIOC,T_INFO%NIONS)

            !c some atom might have escaped from the original cell so we 
            !c need to update WHERE in definition of ICOORDINATES
            CALL update_where(T_INFO,DYN%POSIOC,TRANSPOSE(LATT_CUR%A),ICOORDINATES_GL)
              
              
            IF (TILAMBDA_GL .GE. 0.) THEN
              !c update TILAMBDA if the user defined a transformation speed TIDOT
              IF (LTIDOT_GL) THEN
                TILAMBDA_GL=TILAMBDA_GL+TIDOT_GL
                !c make sure TILAMBDA stays in interval <0,1>
                IF (TILAMBDA_GL .GT. 1._q) TILAMBDA_GL=1._q
                IF (TILAMBDA_GL .LT. 0._q) TILAMBDA_GL=0._q
                IF (IO%IU6>0) THEN
                  write(IO%IU6,FMT='(3X,A22,X,F10.8)') '           TILAMBDA = ',TILAMBDA_GL
                  write(g_io%REPORT,FMT='(3X,A22,X,F10.8)') '           TILAMBDA = ',TILAMBDA_GL
                ENDIF
              ENDIF 

              CALL minimize_difference(POSITIONS0_GL,DYN%POSIOC,T_INFO%NIONS)
              CALL HARMONIC_FORCEFIELD_COMPUTE(IO,g_io,LATT_CUR,DYN,T_INFO,TILAMBDA_GL,POSITIONS0_GL,ICOORDINATES_GL,&
              &  HARMONIC_HESSIAN_GL,EPOT,HARMONIC_ENERGY0_GL,HARMONIC_ENERGY)
            ENDIF

            !c update internal coordinates
            CALL DEAL_XYZ(T_INFO,DYN%POSIOC,REF_C,LATT_CUR%A,ICOORDINATES_GL)

            !c monitor coordinates with STATUS==7, stop if we get out of 
            !c limits
            IF (ICOORDINATES_GL%iconst7>0) CALL tis_stop(DYN,T_INFO,INFO,LATT_CUR,IO,g_io,ICOORDINATES_GL,LMIN_GL,LMAX_GL,VALUEA_GL,VALUEB_GL)

            !c monitor values of coordinates affected by harmonic potential given by HESSEMAT 
            IF ((ICOORDINATES_GL%iconst3>0) .AND. (IO%IU6>0)) THEN
              write(g_io%REPORT,FMT='(/,2X,A15)') '>Harmonic_coord'
              DO i=1,ICOORDINATES_GL%NUMINTERNALS
                IF (ICOORDINATES_GL%COORDSTRUCT(i)%STATUS==3) THEN
                  !write(g_io%REPORT,FMT='(3X,A3,X,A2,X,F15.5)') &
                  !&  'hc>',ICOORDINATES_GL%COORDSTRUCT(i)%TAG,ICOORDINATES_GL%COORDSTRUCT(i)%VALUE
                  write(g_io%REPORT,FMT='(3X,A3,X,A2,X,F15.5,X,F15.5)') &
                  &  'hc>',ICOORDINATES_GL%COORDSTRUCT(i)%TAG,ICOORDINATES_GL%COORDSTRUCT(i)%VALUE, ICOORDINATES_GL%COORDSTRUCT(i)%DVALUE
                ENDIF
              ENDDO
            ENDIF
            
            !c monitor values of coordinates affected by harmonic bias potential 
            IF ((ICOORDINATES_GL%iconst8>0) .AND. (IO%IU6>0)) THEN
              write(g_io%REPORT,FMT='(/,2X,A12)') '>Spring_coord'
              DO i=1,ICOORDINATES_GL%NUMINTERNALS
                IF (ICOORDINATES_GL%COORDSTRUCT(i)%STATUS==8) THEN
                  write(g_io%REPORT,FMT='(3X,A3,X,A2,X,F15.5)') &
                  &  'sc>',ICOORDINATES_GL%COORDSTRUCT(i)%TAG,ICOORDINATES_GL%COORDSTRUCT(i)%VALUE
                ENDIF
              ENDDO
            ENDIF

            !c monitor values of coordinates affected by Fermi-shape bias potential 
            IF ((ICOORDINATES_GL%iconst4>0) .AND. (IO%IU6>0)) THEN
              write(g_io%REPORT,FMT='(/,2X,A12)') '>Fbias_coord'
              DO i=1,ICOORDINATES_GL%NUMINTERNALS
                IF (ICOORDINATES_GL%COORDSTRUCT(i)%STATUS==4) THEN
                  write(g_io%REPORT,FMT='(3X,A3,X,A2,X,F15.5)') &
                  &  'fc>',ICOORDINATES_GL%COORDSTRUCT(i)%TAG,ICOORDINATES_GL%COORDSTRUCT(i)%VALUE
                ENDIF
              ENDDO
            ENDIF
          ENDIF

          !c rescale ionic velocities
          IF (DYN%ISIF==2 .AND. (MDALGO_GL>0 .AND. MDALGO_GL/=13) .AND. SCALING_GL >0 ) THEN
            IF (MOD(counter,SCALING_GL)==0) THEN
              CALL KINETIC_E(EKIN,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,&
                       & DYN%POTIM,LATT_CUR%A,DYN%VEL)
              TEMPER=2*EKIN/BOLKEV/NDEGREES_OF_FREEDOM
              DYN%VEL=DYN%VEL*(DYN%TEMP/TEMPER)**0.5

              IF (IO%IU6>0) THEN
                  write(g_io%REPORT,FMT='(/,2X,A29)') 'Velocities have been rescaled'
              ENDIF
            ENDIF
          ENDIF


          !c take care of the cetenter-of-mass motion
          IF (.NOT. T_INFO%LSDYN) THEN
            !c remove drift from velocities
            IF (LCMASS) CALL SYMVEL_WARNING(T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS, &
            &            DYN%VEL,IO%IU6,IO%IU0)
            !c remove drift from forces
            CALL SYMVEL(T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS, &
            &            DYN%POSION,DYN%D2C,LATT_CUR%A,LATT_CUR%B)
          ELSE
            !c remove velocities for frozen directions
            CALL SET_SELECTED_VEL_ZERO(T_INFO, DYN%VEL,LATT_CUR)
          ENDIF


          SELECT CASE(MDALGO_GL)
            !c Andersen thermostat (direct metaDynamics)
            CASE(1)
              CALL HILLS_METHOD_DIRECT(DYN,ICOORDINATES_GL,hills_GL,penalty_GL,T_INFO,INFO,LATT_CUR, &
              &      EKIN,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,IO,EPOT,ANDERSEN_PROB_GL,&
              &      g_io,ICOORDINATES_GL%iconst0,ICOORDINATES_GL%iconst2,ICOORDINATES_GL%iconst4,&
              &      ICOORDINATES_GL%iconst5,ICOORDINATES_GL%iconst8,hspring_R0_GL, & 
              &      hspring_K_GL,hspring_force_GL,fbias_R0_GL,fbias_D_GL,fbias_A_GL,TEIN)
              
            !c Nose-Hoover thermostat (direct metaDynamics)
            CASE(2)
              CALL HILLS_METHOD_DIRECT_nose(DYN,ICOORDINATES_GL,hills_GL,penalty_GL,T_INFO,INFO,LATT_CUR, &
              &      EKIN,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,IO,EPOT,ANDERSEN_PROB_GL,& 
              &      g_io,ICOORDINATES_GL%iconst0,ICOORDINATES_GL%iconst2,ICOORDINATES_GL%iconst4,& 
              &      ICOORDINATES_GL%iconst5,ICOORDINATES_GL%iconst8,hspring_R0_GL,& 
              &      hspring_K_GL,hspring_force_GL,fbias_R0_GL,fbias_D_GL,fbias_A_GL,TEIN)
              IF (IO%IU6>0) CALL NOSE_OUT(g_io,DYN)
          
            !c Langevin dynamics, velocity verlet (to be replaced by leap-frog)
            CASE(3)
              IF (DYN%ISIF==3) THEN
                CALL STEP_LANGEVIN_ISIF3(DYN,ICOORDINATES_GL,hills_GL,penalty_GL,T_INFO,INFO,LATT_CUR, &
                &      EKIN,EKIN_Lat,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,NDEGREES_OF_FREEDOM_,IO,&
                &      EPOT,GAMMA_GL,GAMMA_L_GL,g_io,SIF,ICOORDINATES_GL%iconst0,ICOORDINATES_GL%iconst2,& 
                &      ICOORDINATES_GL%iconst4,ICOORDINATES_GL%iconst5, &
                &      ICOORDINATES_GL%iconst8,hspring_R0_GL,hspring_K_GL,hspring_force_GL,fbias_R0_GL,&
                &      fbias_D_GL,fbias_A_GL,TEIN,AMASS_GL,LATTICE_CONSTRAINTS_GL,LATTICE_DOF_GL)
              ELSE 
                CALL STEP_LANGEVIN(DYN,ICOORDINATES_GL,hills_GL,penalty_GL,T_INFO,INFO,LATT_CUR,EKIN,EPS,ES,DISMAX,&
                &      NDEGREES_OF_FREEDOM, &
                &      NDEGREES_OF_FREEDOM_,IO,EPOT,GAMMA_GL,g_io,ICOORDINATES_GL%iconst0,ICOORDINATES_GL%iconst2,&
                &      ICOORDINATES_GL%iconst4,ICOORDINATES_GL%iconst5,0,&
                &      ICOORDINATES_GL%iconst8,hspring_R0_GL,hspring_K_GL,hspring_force_GL,fbias_R0_GL,fbias_D_GL,fbias_A_GL,TEIN)
              ENDIF

            !c Andersen thermostat (metaDynamics with fictituous variables - this is presently a "hidden" feature)
            CASE(10)
              CALL HILLS_METHOD(DYN,ICOORDINATES_GL,hills_GL,penalty_GL,T_INFO,INFO,LATT_CUR, &
              &      EKIN,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,IO,EPOT,ANDERSEN_PROB_GL,&
              &      g_io,WDES,ICOORDINATES_GL%iconst0,ICOORDINATES_GL%iconst2,TEIN)

            
            !c Nose-Hoover thermostat (metaDynamics with fictituous variables - this is presently a "hidden" feature)
            CASE(20)
              CALL HILLS_METHOD_nose(DYN,ICOORDINATES_GL,hills_GL,penalty_GL,T_INFO,INFO,LATT_CUR, &
              &      EKIN,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,IO,EPOT,ANDERSEN_PROB_GL,& 
              &      g_io,WDES,ICOORDINATES_GL%iconst0,ICOORDINATES_GL%iconst2,TEIN)
              IF (IO%IU6>0) CALL NOSE_OUT(g_io,DYN)

            !c dynamics with multiple Andersen thermostats
            CASE(13)
              CALL STEP_NANDERSEN(DYN,ICOORDINATES_GL,T_INFO,INFO,LATT_CUR,EKIN,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,&
              &              NDEGREES_OF_FREEDOM_,IO,EPOT,g_io,NSUBSYS_GL,TSUBSYS_GL,PSUBSYS_GL,ICOORDINATES_GL%iconst0,ICOORDINATES_GL%iconst2,TEIN)

            !c use original MD routines from vasp before version 5.4.4    
            CASE DEFAULT
              IF ( DYN%ISIF==3 .OR. DYN%ISIF==8 .OR. DYN%ISIF==9 .OR. DYN%ISIF==10 ) then 
                CALL step_NPT( DYN%INIT, DYN%ISIF, NDEGREES_OF_FREEDOM+3, T_INFO%NIONS, T_INFO%NITYP,  &
                &       T_INFO%NTYP, LATT_CUR,DYN%PSTRESS, &
                &       EPOT, &
                &       T_INFO%POMASS, DYN%POSION, &
                &       DYN%POSIOC, DYN%POTIM, DYN%SMASS, DYN%SNOSE, DYN%TEMP, TIFOR,SIF, DYN%VEL, &
                &       DYN%D2C, DYN%D2, IO%IU6,IO%IU5,IO%IU0)
                CALL KINETIC_E(EKIN,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,&
                        & DYN%POTIM,LATT_CUR%A,DYN%D2)             
                TEIN = 2*EKIN/BOLKEV/(NDEGREES_OF_FREEDOM+3)
              ELSE
                !c standard MD module
                CALL STEP(DYN%INIT,ISCALE,T_INFO%NIONS,LATT_CUR%A,LATT_CUR%ANORM,DYN%D2C,DYN%SMASS,DYN%POSION,DYN%POSIOC, &
                 DYN%POTIM,T_INFO%POMASS,T_INFO%NTYP,T_INFO%ITYP,DYN%TEMP,DYN%VEL,DYN%D2,DYN%D3,DYN%SNOSE, &
                 EKIN,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM, IO%IU6)
                TEIN = 2*EKIN/BOLKEV/NDEGREES_OF_FREEDOM
              ENDIF
          END SELECT   
          
#ifdef MPI
          CALL M_bcast_d(WDES%COMM, DYN%POSION , T_INFO%NIONS*3)
          CALL M_bcast_d(WDES%COMM, DYN%VEL , T_INFO%NIONS*3)
          CALL M_bcast_d(WDES%COMM, LATT_CUR%A, 9)
          CALL M_bcast_d(WDES%COMM, LATT_CUR%Avel, 9) 
          CALL LATTIC(LATT_CUR)
#endif

          IF (MDALGO_GL>0) THEN
            DYN%INIT=1
            IF (DYN%ISIF>=3) THEN
              IF (LATT_CUR%INITlatv==0) LATT_CUR%INITlatv=1
            ENDIF
            CALL put_in_box(T_INFO%NIONS,DYN%POSION)
            !c on output, POSIOC contains the previous step
            !c and POSION the current one
            !DYN%POSIOC=DYN%POSION
            IF (IO%IU6>=0) THEN
              CALL WFORCE(g_io%REPORT)
            ENDIF
          ENDIF

          !c stepmd
          !CALL XML_CLOSE_TAG
        END SUBROUTINE STEP_tb

        SUBROUTINE DYN_XML_WRITE(MDALGO,ANDERSEN_PROB,SCALING,NSUBSYS,TSUBSYS,PSUBSYS)
          INTEGER IDUM,SCALING,NIONS,MDALGO
          LOGICAL ::LDUM
          REAL(q) :: RDUM,ANDERSEN_PROB
          COMPLEX(q) :: CDUM
          CHARACTER*1 :: CHARAC
          INTEGER :: NSUBSYS(3)
          REAL(q) :: TSUBSYS(3),PSUBSYS(3)

          CALL XML_TAG("parameters_embedded")
          CALL XML_TAG("separator","mdynamics")
          CALL XML_INCAR('MDALGO','I',MDALGO,RDUM,CDUM,LDUM,CHARAC,1)
          IF (MDALGO==1 .OR. MDALGO==10) THEN
            CALL XML_INCAR('ANDERSEN_PROB','F',IDUM,ANDERSEN_PROB,CDUM,LDUM,CHARAC,1)
         
                !write(g_io%REPORT,FMT='(3X,A22,X,I6)')    '            SCALING = ' , SCALING
          ELSE IF (MDALGO==13) THEN   
            CALL XML_INCAR_V('NSUBSYS','I',NSUBSYS,RDUM,CDUM,LDUM,CHARAC,3)
            CALL XML_INCAR_V('TSUBSYS','F',IDUM,TSUBSYS,CDUM,LDUM,CHARAC,3)
            CALL XML_INCAR_V('PSUBSYS','F',IDUM,PSUBSYS,CDUM,LDUM,CHARAC,3)
          ENDIF  
          CALL XML_CLOSE_TAG
          CALL XML_CLOSE_TAG
       
        END SUBROUTINE DYN_XML_WRITE

        SUBROUTINE DYN_READER(IO,g_io,MDALGO,ANDERSEN_PROB,SCALING,NSUBSYS,TSUBSYS,PSUBSYS,GAMM,&
        &  GAMMA_L,NIONS,NTYP,ISIF,CNEXP,TILAMBDA,TIDOT,LTIDOT,AMASS,LATTICE_CONSTRAINTS,LATTICE_DOF)
        !c read general input parameters for MD
          USE reader_tags

          TYPE (in_struct) :: IO
          TYPE(gadget_io) :: g_io
          INTEGER IDUM,THERMOSTAT,SCALING,NIONS,NTYP,ISIF,MDALGO
          INTEGER :: K,N,IERR,NSUB,i
          LOGICAL :: LOPEN,LDUM
          REAL(q) :: RDUM,ANDERSEN_PROB
          COMPLEX(q) :: CDUM
          CHARACTER*1 :: CHARAC
          INTEGER :: NSUBSYS(3)
          REAL(q) :: TSUBSYS(3),PSUBSYS(3)
          REAL(q) :: GAMM(NTYP) !c friction coeficient in ps^(-1)
          REAL(q) :: GAMMA_L !c friction coeficient in ps^(-1) for lattice parameters
          REAL(q) :: CNEXP(2)
          REAL(q) :: TILAMBDA,TIDOT
          LOGICAL :: LTIDOT
          LOGICAL :: LATTICE_CONSTRAINTS(3) ! specifies which of the lattice parameters (a_1,a_2,a_3) are free (true) or fixed (false)
          INTEGER :: LATTICE_DOF ! the total number of the lattice degrees of freedom according to the applied lattice constra
          REAL(q) :: AMASS


          CALL OPEN_INCAR_IF_FOUND(IO%IU5, LOPEN)

          !c THERMOSTAT - obsolete, replaced by MDALGO
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'THERMOSTAT', THERMOSTAT, IERR, WRITEXMLINCAR)
          IF (IERR==3) THERMOSTAT=0
          IF (THERMOSTAT<0) THERMOSTAT=0

          !c MDALGO
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'MDALGO', MDALGO, IERR, WRITEXMLINCAR)
          IF (IERR==3) MDALGO=THERMOSTAT
          IF (MDALGO<0) MDALGO=0
          !c replace the obsolete values
          IF (MDALGO==11) THEN
            MDALGO=1
            CALL vtutor%warning("the obsolete value of MDALGO=11 reset to 1")
          ELSE IF (MDALGO==21) THEN
            MDALGO=2
            CALL vtutor%warning("the obsolete value of MDALGO=21 reset to 2")
          ELSE IF (MDALGO==31) THEN
            MDALGO=3
            CALL vtutor%warning("the obsolete value of MDALGO=31 reset to 3")
          ENDIF

          !c basic sanity test of MDALGO setting
          IF (.NOT. ((MDALGO==0) .OR. (MDALGO==1) .OR. (MDALGO==10) .OR. (MDALGO==13) .OR. &
          &   (MDALGO==2) .OR. (MDALGO==20) .OR. (MDALGO==3))) THEN
            CALL vtutor%error("invalid value of MDALGO:")
          ENDIF
 

          !c multiple thermostat stuff
          IF (MDALGO==13) THEN
            !c NSUBSYS - subsystems with different temperature
            CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'NSUBSYS', NSUBSYS, 3, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)
            DO i=1,NIONS
              IF (NSUBSYS(i) .GT. NIONS) THEN
                CALL vtutor%error("NSUBSYS(i) must not be greater than the number of atoms defined in &
                   &POSCAR")
              END IF
            END DO
            IF (N==2 .AND. NSUBSYS(2)==NIONS) N=1
            IF (N==3 .AND. NSUBSYS(3)==NIONS) N=2
            IF ((N==1 .OR. N==2) .AND. (NSUBSYS(1)<NSUBSYS(2)) .AND. (NSUBSYS(2) .LE. NSUBSYS(3))) THEN
              NSUB=N+1
              CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'PSUBSYS', PSUBSYS, 3, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)
              IF (N == NSUB) THEN
                CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'TSUBSYS', TSUBSYS, 3, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)
                IF (N /= NSUB) THEN
                  CALL vtutor%error("Error reading item 'TSUBSYS' from file INCAR")
                ENDIF
                DO i=1,NSUB
                  IF (TSUBSYS(i) .LE. 0) THEN
                    CALL vtutor%error("zero in item 'TSUBSYS' from file INCAR, all TSUBSYS(i) must be &
                     &greater than 0")
                  END IF
                ENDDO
              ELSE
                CALL vtutor%error("Error reading item 'PSUBSYS' from file INCAR")
              END IF
            ELSE
              CALL vtutor%error("Error reading item 'NSUBSYS' from file INCAR, two or three &
                 &sub-systems must be defined")
            END IF           
          END IF

          !c ANDERSEN_PROB - probability of collision
          IF (MDALGO==1 .OR. MDALGO==10 .OR. MDALGO==13) THEN
            CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'ANDERSEN_PROB', ANDERSEN_PROB, IERR, WRITEXMLINCAR)
            IF (IERR==3) ANDERSEN_PROB=0._q 
            IF (ANDERSEN_PROB <0._q) ANDERSEN_PROB=0._q
          ENDIF

          IF (MDALGO==3) THEN
            !c friction coeficients for atomic positions
            CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LANGEVIN_GAMMA', GAMM, NTYP, IERR, WRITEXMLINCAR)
            IF (IERR==3) GAMM=0._q 
            !c terminate if negative value of GAMM is defined 
            DO i=1,NTYP
              IF (GAMM(i) <0._q) THEN
                CALL vtutor%error("Error reading item 'LANGEVIN_GAMMA' from file INCAR: all &
                   &LANGEVIN_GAMMA must be positive numbers")
              ENDIF
            ENDDO
           
            !c friction coeficient for lattice parameters
            GAMMA_L=0._q
            CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LANGEVIN_GAMMA_L', GAMMA_L, IERR, WRITEXMLINCAR)
            IF (IERR==3) GAMMA_L=0._q
            IF (GAMMA_L <0._q) THEN
              CALL vtutor%error("Error reading item 'LANGEVIN_GAMMA_L' from file INCAR: &
                   &LANGEVIN_GAMMAL must be a positive number")
            ENDIF 
          ENDIF


          !c SCALING - in how many steps should the velocities be rescaled?
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'SCALING', SCALING, IERR, WRITEXMLINCAR)
          IF (IERR==3) SCALING=0
          IF (SCALING<0) SCALING=0

          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'CNEXP', CNEXP, 2, IERR, WRITEXMLINCAR)
          IF (IERR==3) THEN 
            CNEXP(1)=9.;CNEXP(2)=14.
          ENDIF
          
          !c terminate if negative values of CNEXP are defined 
          DO i=1,2
            IF (CNEXP(i) <0._q) THEN
              CALL vtutor%error("Error reading item 'CNEXP' from file INCAR: all CNEXP must be &
                 &positive numbers")
            ENDIF 
          ENDDO

          TILAMBDA=-1.
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'TILAMBDA', TILAMBDA, IERR, WRITEXMLINCAR)
          IF (IERR==3) THEN 
            TILAMBDA=-1.
          ENDIF
          IF ((IERR==0) .AND. ((TILAMBDA .LT. 0._q) .OR. (TILAMBDA .GT. 1._q))) THEN
            CALL vtutor%error("Error reading item ''TILAMBDA'' from file INCAR: TILAMBDA must be from interval <0,1>")
          ENDIF

          IF (TILAMBDA>=0.) THEN
            CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'TIDOT', TIDOT, IERR, WRITEXMLINCAR)
            IF (IERR==3) THEN
              TIDOT=0._q
              LTIDOT=.FALSE.
            !ELSE IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
            !  CALL vtutor%error("Error reading item ''TIDOT'' from file INCAR")
            ELSE IF ((ABS(TIDOT) .GT. 0.1_q)) THEN
              CALL vtutor%error("Error reading item ''TIDOT'' from file INCAR: ABS(TIDOT) must be less than 0.1")
            ELSE IF ((TILAMBDA .GE. 1._q) .AND. (TIDOT .GT. 0._q)) THEN
              CALL vtutor%error("Error reading item ''TIDOT'' from file INCAR: TIDOT must be negative for TILAMBDA==1")
            ELSE IF ((TILAMBDA .LE. 0._q) .AND. (TIDOT .LT. 0._q)) THEN
              CALL vtutor%error("Error reading item ''TIDOT'' from file INCAR: TIDOT must be positive for TILAMBDA==0")
            ELSE
              LTIDOT=.TRUE.
            ENDIF
          ENDIF


          IF (ISIF==3) THEN
            CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'PMASS', AMASS, IERR, WRITEXMLINCAR)
            IF (IERR==3) Amass=1000._q

            LATTICE_CONSTRAINTS = .TRUE.
            CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LATTICE_CONSTRAINTS', LATTICE_CONSTRAINTS, 3, IERR, WRITEXMLINCAR)
            IF (IERR==3) LATTICE_CONSTRAINTS = .TRUE.

            ! count the lattice degrees of freedom
            LATTICE_DOF = 0
            DO I = 1, 3
              IF (LATTICE_CONSTRAINTS(I)) LATTICE_DOF = LATTICE_DOF + 1
            END DO
            LATTICE_DOF = LATTICE_DOF * LATTICE_DOF

         ENDIF


          !c write out some outpout
          IF (MDALGO>0) THEN
            IF (IO%IU6>0) THEN
              !c write the input parameter summary into OUTCAR
              write(IO%IU6,FMT='(/,2X,A22)') 'MD-specific parameters'

              IF ((TILAMBDA.GE. 0.) .AND. (.NOT. LTIDOT)) THEN
                write(IO%IU6,FMT='(3X,A22,X,F10.8)') '           TILAMBDA = ',TILAMBDA
              ELSE IF ((TILAMBDA_GL .GE. 0.) .AND. (LTIDOT)) THEN
                write(IO%IU6,FMT='(3X,A22,X,F10.8)') '           TIDOT= ',TIDOT
              ENDIF


              write(IO%IU6,FMT='(3X,A22,X,I2)')    '             MDALGO = ', MDALGO
              IF (MDALGO==1 .OR. MDALGO==10) THEN
                write(IO%IU6,FMT='(3X,A22,X,F10.8)')         '      ANDERSEN_PROB = ',ANDERSEN_PROB
              ENDIF
              IF (MDALGO==3) THEN
                write(IO%IU6,ADVANCE='NO',FMT='(3X,A22)')  '     LANGEVIN_GAMMA = '
                DO i=1,NTYP
                  write(IO%IU6,ADVANCE='NO',FMT='(X,F8.3)') GAMM(i)
                ENDDO
                write(IO%IU6,ADVANCE='NO',FMT='(/)')
                IF (ISIF==3) THEN
                  write(IO%IU6,FMT='(3X,A22,X,F8.3)')  '   LANGEVIN_GAMMA_L = ',GAMMA_L
                ENDIF
              ENDIF

              !write(g_io%REPORT,FMT='(3X,A22,X,I6)')    '            SCALING = ' , SCALING
              IF (MDALGO==13) THEN
                write(IO%IU6,FMT='(3X,A22,X,I5,X,I5,X,I5)')    '           NSUBSYS = ', NSUBSYS(1), NSUBSYS(2),NSUBSYS(3)
                write(IO%IU6,FMT='(3X,A22,X,F8.3,X,F8.3,X,F8.3)')    '           TSUBSYS = ', TSUBSYS(1), TSUBSYS(2),TSUBSYS(3)
                write(IO%IU6,FMT='(3X,A22,X,F8.3,X,F8.3,X,F8.3)')    '           PSUBSYS = ', PSUBSYS(1), PSUBSYS(2),PSUBSYS(3)
              ENDIF

              !c write the input parameter summary into REPORT
              IF ((TILAMBDA.GE. 0.) .AND. (.NOT. LTIDOT)) THEN
                write(g_io%REPORT,FMT='(3X,A22,X,F10.8)') '           TILAMBDA = ',TILAMBDA
              ELSE IF ((TILAMBDA.GE. 0.) .AND. (LTIDOT)) THEN
                write(g_io%REPORT,FMT='(3X,A22,X,F10.8)') '           TIDOT= ',TIDOT
              ENDIF

              write(g_io%REPORT,FMT='(3X,A22,X,I2)')    '             MDALGO = ', MDALGO
              IF (MDALGO==1 .OR. MDALGO==10) THEN
                write(g_io%REPORT,FMT='(3X,A22,X,F10.8)') '      ANDERSEN_PROB = ',ANDERSEN_PROB
              ENDIF

              IF (MDALGO==3) THEN
                write(g_io%REPORT,ADVANCE='NO',FMT='(3X,A22)')  '     LANGEVIN_GAMMA = '
                DO i=1,NTYP
                  write(g_io%REPORT,ADVANCE='NO',FMT='(X,F8.3)') GAMM(i)
                ENDDO
                write(g_io%REPORT,ADVANCE='NO',FMT='(/)')
                IF (ISIF==3) THEN
                  write(g_io%REPORT,FMT='(3X,A22,X,F8.3)')  '   LANGEVIN_GAMMA_L = ',GAMMA_L
                ENDIF
              ENDIF

              write(g_io%REPORT,FMT='(3X,A22,X,I6)')    '            SCALING = ' , SCALING
              IF (MDALGO==13) THEN
                write(g_io%REPORT,FMT='(3X,A22,X,I5,X,I5,X,I5)')    '           NSUBSYS = ', NSUBSYS(1), NSUBSYS(2),NSUBSYS(3)
                write(g_io%REPORT,FMT='(3X,A22,X,F8.3,X,F8.3,X,F8.3)')    '           TSUBSYS = ', TSUBSYS(1), TSUBSYS(2),TSUBSYS(3)
                write(g_io%REPORT,FMT='(3X,A22,X,F8.3,X,F8.3,X,F8.3)')    '           PSUBSYS = ', PSUBSYS(1), PSUBSYS(2),PSUBSYS(3)
              ENDIF

              write(IO%IU6,FMT='(3X,A22,X,F8.3,X,F8.3)')   '              CNEXP = ',CNEXP(1),CNEXP(2)
              write(g_io%REPORT,FMT='(3X,A22,X,F8.3,X,F8.3)')      '              CNEXP = ',CNEXP(1),CNEXP(2)
            ENDIF
          ENDIF

          IF ((MDALGO==3) .AND. (ISIF==3)) THEN
            IF (IO%IU6>0) THEN
              write(g_io%REPORT,FMT='(3X,A22,X,F10.5)')    '              PMASS = ', AMASS  
              write(g_io%REPORT,FMT='(3X,A22,X,3L,A,I1,A)') 'LATTICE_CONSTRAINTS = ', (LATTICE_CONSTRAINTS(I), I=1,3), &
              & ' (X,Y,Z) = ', LATTICE_DOF, ' degree(s) of freedom'
            ENDIF
         ENDIF


          CALL CLOSE_INCAR_IF_FOUND(IO%IU5)

        END SUBROUTINE DYN_READER

        SUBROUTINE READER_ICONST0(IO,g_io,iconst0,BMTOL,BMTOLSOFT,BMSCA,BMITER,LBLUEOUT,EQUI_REGIME,INCREM)
        !c read input parameters specific for constrained dynamics
          USE reader_tags

          TYPE (in_struct) :: IO
          TYPE(gadget_io) :: g_io
! fixme[vasp]: shouldnt ii be initialized to 1?
          INTEGER IDUM, K,N,IERR,i,ii,iconst0
          LOGICAL :: LOPEN,LDUM
          REAL(q) :: RDUM
          COMPLEX(q) :: CDUM
          CHARACTER*1 :: CHARAC
          REAL(q) :: BMTOL,BMTOLSOFT,BMSCA,INCREM(iconst0)
          INTEGER :: BMITER,EQUI_REGIME
          LOGICAL :: LBLUEOUT

          if (iconst0<1) RETURN

          CALL OPEN_INCAR_IF_FOUND(IO%IU5, LOPEN)
            
          !c tolerance for the shake algorithm
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'SHAKETOL', BMTOL, IERR, WRITEXMLINCAR)
          IF (IERR==3) BMTOL=1e-5
          IF (BMTOL<0.) BMTOL=1e-5

          !c soft tolerance for the shake algorithm
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'SHAKETOLSOFT', BMTOLSOFT, IERR, WRITEXMLINCAR)
          IF (IERR==3) BMTOLSOFT=BMTOL
          IF (BMTOLSOFT<0.) BMTOLSOFT=BMTOL

          !c scaling of step taken in one iteration in the SHAKE algorithm
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'SHAKESCA', BMSCA, IERR, WRITEXMLINCAR)
          IF (IERR==3) BMSCA=2._q
          IF (BMSCA<0.) BMSCA=2._q

          !c maximal number of iterations in the SHAKE algorithm
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'SHAKEMAXITER', BMITER, IERR, WRITEXMLINCAR)
          IF (IERR==3) BMITER=1000
          IF (BMITER<10) BMITER=1000

          !c compute and write down the output for blue-moon?
          !c might be time-consuming, useful only for the free-energy calculations
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LBLUEOUT', LBLUEOUT, IERR, WRITEXMLINCAR)
          IF (IERR==3) LBLUEOUT=.FALSE.

          !c EQUI_REGIME - length of the equilibration preriod (don't constrain
          !c anything during this period)
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'EQUI_REGIME', EQUI_REGIME, IERR, WRITEXMLINCAR)
          IF (IERR==3) EQUI_REGIME=0
          IF (EQUI_REGIME<0) EQUI_REGIME=0

          !c INCREM - the values of coordinates with STATUS==0 defined in ICONST
          !c are increased by this ammount 
          IF (iconst0>0) THEN 
            INCREM=0._q
            CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'INCREM', INCREM, ICONST0, IERR, WRITEXMLINCAR)
            IF (IERR==3) INCREM=0._q
          END IF

          IF (IO%IU6>0) THEN
            !write(g_io%REPORT,FMT='(/,2X,A31)')    '>parameters for SHAKE algorithm'
            write(g_io%REPORT,FMT='(3X,A22,X,I5)')    '   EQUI_REGIME = ',EQUI_REGIME_GL
            write(g_io%REPORT,FMT='(3X,A22,X,I8)')         'SHAKEMAXITER = ',BMITER_GL
            write(g_io%REPORT,FMT='(3X,A22,X,E13.6E2)')    '    SHAKETOL = ',BMTOL_GL
            write(g_io%REPORT,FMT='(3X,A22,X,E13.6E2)')    'SHAKETOLSOFT = ',BMTOLSOFT_GL
            write(g_io%REPORT,FMT='(3X,A22,X,F12.8)')      '    SHAKESCA = ',BMSCA_GL
          END IF 

          CALL CLOSE_INCAR_IF_FOUND(IO%IU5)
        END SUBROUTINE READER_ICONST0

        SUBROUTINE READER_ICONST7(IO,g_io,ICOORDINATES,LMIN,LMAX,VALUEA,VALUEB)
        !c read input parameters specific coordinates with STATUS=7
          USE reader_tags

          TYPE (in_struct) :: IO
          TYPE(gadget_io) :: g_io
          TYPE(coordstructure) :: ICOORDINATES
          INTEGER IDUM, K,N,IERR,i,ii
          LOGICAL :: LOPEN,LDUM
          REAL(q) :: RDUM
          COMPLEX(q) :: CDUM
          CHARACTER*1 :: CHARAC
          LOGICAL :: LMIN,LMAX
          REAL(q) :: VALUEA(ICOORDINATES%iconst7),VALUEB(ICOORDINATES%iconst7)

          IF (ICOORDINATES%iconst7<1)  RETURN

          CALL OPEN_INCAR_IF_FOUND(IO%IU5, LOPEN)

          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'VALUE_MIN', VALUEA, ICOORDINATES%ICONST7, IERR, WRITEXMLINCAR)
          IF (IERR/=0) THEN
            LMIN=.FALSE.
          ENDIF

          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'VALUE_MAX', VALUEB, ICOORDINATES%ICONST7, IERR, WRITEXMLINCAR)
          IF (IERR/=0) THEN
            LMAX=.FALSE.
          ENDIF
              
          CALL CLOSE_INCAR_IF_FOUND(IO%IU5)

          DO i=1,ICOORDINATES%iconst7
            IF (VALUEB(i)<VALUEA(i)) THEN
              LMIN=.FALSE.
              LMAX=.FALSE.
            END IF
          ENDDO

          IF (LMIN) THEN
            IF (IO%IU6>0) THEN
              write(g_io%REPORT,FMT='(3X,A22)',ADVANCE='NO')  '          VALUE_MIN = '
              ii=0
              DO i=1,ICOORDINATES%NUMINTERNALS
                IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==7) THEN
                  ii=ii+1
                  write(g_io%REPORT,FMT='(X,F9.5)',ADVANCE='NO') VALUEA(ii)
                ENDIF
              ENDDO
              write(g_io%REPORT,FMT='(X)',ADVANCE='YES')
            ENDIF
          ENDIF

          IF (LMAX) THEN
            IF (IO%IU6>0) THEN
              write(g_io%REPORT,FMT='(3X,A22)',ADVANCE='NO')  '          VALUE_MAX = '
              ii=0
              DO i=1,ICOORDINATES%NUMINTERNALS
                IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==7) THEN
                  ii=ii+1
                  write(g_io%REPORT,FMT='(X,F9.5)',ADVANCE='NO') VALUEB(ii)
                ENDIF
              ENDDO
              write(g_io%REPORT,FMT='(X)',ADVANCE='YES')
            ENDIF
          ENDIF
           
        END SUBROUTINE READER_ICONST7

        SUBROUTINE READ_HARMONIC_FORCEFIELD(IO,g_io,DYN,T_INFO,HARMONIC_HESSIAN,&
        &  HARMONIC_ENERGY0,POSITIONS0)
          TYPE (in_struct) :: IO
          TYPE(gadget_io) :: g_io
          TYPE(dynamics) :: DYN
          TYPE(type_info) :: T_INFO
          REAL(q) :: HARMONIC_HESSIAN(3*T_INFO%NIONS,3*T_INFO%NIONS)
          REAL(q) :: POSITIONS0(3,T_INFO%NIONS),d(3*T_INFO%NIONS)
          REAL(q) :: HARMONIC_ENERGY0
          REAL(q) :: EIGENMODES(3*T_INFO%NIONS,3*T_INFO%NIONS)
          INTEGER :: ios,i,j,iblock,imax,counter
          
          OPEN(UNIT=g_io%HESSEMAT,FILE='HESSEMAT',STATUS='UNKNOWN')
          iblock=T_INFO%NIONS+1
          imax=(3*T_INFO%NIONS)*iblock
          
          !c read in energy and positions
          POSITIONS0=0._q
          READ(g_io%HESSEMAT,FMT=*,IOSTAT=ios) HARMONIC_ENERGY0
          DO i=1,T_INFO%NIONS
            READ(g_io%HESSEMAT,FMT=*,IOSTAT=ios)  POSITIONS0(1,i),POSITIONS0(2,i),POSITIONS0(3,i)
          ENDDO  
          
          !c read in force constants and eigenvectors
          HARMONIC_HESSIAN=0._q;EIGENMODES=0._q
          counter=0;j=0
          DO i=1,imax
            IF (MOD(i,iblock)==1) THEN   
              counter=counter+1
              j=0
              READ(g_io%HESSEMAT,FMT=*,IOSTAT=ios) HARMONIC_HESSIAN(counter,counter)
            ELSE
              READ(g_io%HESSEMAT,FMT=*,IOSTAT=ios) EIGENMODES(j+1,counter),EIGENMODES(j+2,counter),EIGENMODES(j+3,counter)
              j=j+3
            ENDIF
          ENDDO
          
          CLOSE(g_io%HESSEMAT)
          
          !HARMONIC_HESSIAN=MATMUL(HARMONIC_HESSIAN,EIGENMODES)
          !HARMONIC_HESSIAN=MATMUL(TRANSPOSE(EIGENMODES),HARMONIC_HESSIAN)

          HARMONIC_HESSIAN=MATMUL(HARMONIC_HESSIAN,TRANSPOSE(EIGENMODES))
          HARMONIC_HESSIAN=MATMUL((EIGENMODES),HARMONIC_HESSIAN)
          
          !DO i=1,3*T_INFO%NIONS
          !  do j=1,3*T_INFO%NIONS
          !    IF (IO%IU6>0)  write(*,*) "hm1",i,j,HARMONIC_HESSIAN(i,j)
          !  enddo
          !ENDDO

          CALL EIGVAL_GENERAL(HARMONIC_HESSIAN,3*T_INFO%NIONS,d)
          
          !IF (IO%IU6>0) THEN
          !DO i=1,3*T_INFO%NIONS
          !  write(*,*) 'aa',i,d(i)
          !ENDDO
          !ENDIF
        END SUBROUTINE READ_HARMONIC_FORCEFIELD
        
        SUBROUTINE HARMONIC_FORCEFIELD_COMPUTE(IO,g_io,LATT_CUR,DYN,T_INFO,TILAMBDA,POSITIONS0,ICOORDINATES,&
        &  HARMONIC_HESSIAN,EPOT,HARMONIC_ENERGY0,HARMONIC_ENERGY)
          TYPE (in_struct) :: IO
          TYPE(gadget_io) :: g_io
          TYPE(latt) :: LATT_CUR
          TYPE(dynamics) :: DYN
          TYPE(type_info) :: T_INFO
          TYPE(coordstructure) :: ICOORDINATES 
          REAL(q) :: TILAMBDA
          REAL(q) :: HARMONIC_HESSIAN(3*T_INFO%NIONS,3*T_INFO%NIONS)
          REAL(q) :: POSITIONS0(3,T_INFO%NIONS),dummy(3,T_INFO%NIONS)
          REAL(q) :: HARMONIC_ENERGY0,EPOT
          REAL(q) :: HARMONIC_ENERGY
          REAL(q) :: BMAT(ICOORDINATES%NUMINTERNALS,3*T_INFO%NIONS)
          REAL(q) :: HARMONIC_BMAT(ICOORDINATES%iconst3,3*T_INFO%NIONS)
          REAL(q) :: GMAT(ICOORDINATES%iconst3,ICOORDINATES%iconst3)
          REAL(q) :: DMAT(3*T_INFO%NIONS,ICOORDINATES%iconst3)
          REAL(q) :: dq(ICOORDINATES%iconst3),df(ICOORDINATES%iconst3),dfc(3*T_INFO%NIONS)
          REAL(q),SAVE,ALLOCATABLE :: AMAT(:,:)
          LOGICAL,SAVE :: LFIRST=.TRUE.
          INTEGER :: i,j,k,idum,NI,NT
          REAL(q) :: FACT
           REAL(q) :: d(3*T_INFO%NIONS)
          
          !c this routine requires coordinates with status=1
          !c exit it they are not available
          IF (ICOORDINATES%iconst3==0) THEN
            RETURN
          ENDIF
          
          dq=0.;df=0.
          !write(*,*) "ICOORDINATES%iconst3", ICOORDINATES%iconst3 
          IF (LFIRST .EQV. .TRUE.) THEN
            ALLOCATE(AMAT(3*T_INFO%NIONS,ICOORDINATES%iconst3))
            AMAT=0.
            CALL DEAL_XYZ(T_INFO,POSITIONS0,POSITIONS0,LATT_CUR%A,ICOORDINATES)
            CALL BMATRIX(T_INFO,POSITIONS0,POSITIONS0,TRANSPOSE(LATT_CUR%A),ICOORDINATES,BMAT,.TRUE.)
            HARMONIC_BMAT=0.
            j=0
            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==3) THEN
                j=j+1
                HARMONIC_BMAT(j,:)=BMAT(i,:)
                ICOORDINATES%COORDSTRUCT(i)%DVALUE=ICOORDINATES%COORDSTRUCT(i)%VALUE  !xx
              ENDIF
            ENDDO
            GMAT=MATMUL(HARMONIC_BMAT,TRANSPOSE(HARMONIC_BMAT))
            CALL SVDINVERSE(GMAT,ICOORDINATES%iconst3,idum)
            AMAT=MATMUL(TRANSPOSE(HARMONIC_BMAT),GMAT)
            LFIRST=.FALSE.
          ENDIF

          BMAT=0.
          CALL DEAL_XYZ(T_INFO,DYN%POSIOC,DYN%POSIOC,LATT_CUR%A,ICOORDINATES)
          CALL BMATRIX(T_INFO,DYN%POSIOC,DYN%POSIOC,TRANSPOSE(LATT_CUR%A),ICOORDINATES,BMAT,.TRUE.)
          !!CALL BMATRIX(T_INFO,DYN%POSIOC,DYN%POSIOC,TRANSPOSE(LATT_CUR%A),ICOORDINATES,BMAT,.FALSE.)
          HARMONIC_BMAT=0.
          j=0
          DO i=1,ICOORDINATES%NUMINTERNALS
            IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==3) THEN
              j=j+1
              HARMONIC_BMAT(j,:)=BMAT(i,:)
              dq(j)=ICOORDINATES%COORDSTRUCT(i)%VALUE-ICOORDINATES%COORDSTRUCT(i)%DVALUE

              !c avoid problems with straight angles: simply switch off the problematic coordinate
              IF (ICOORDINATES%COORDSTRUCT(i)%TAG=='A ') THEN
                IF (ABS(ICOORDINATES%COORDSTRUCT(i)%VALUE) .GT. 3.1) THEN
                  dq(j)=0.
                  HARMONIC_BMAT(j,:)=0.
                  IF (IO%IU6>0) write(g_io%REPORT,FMT='(3X,A33,I5,A26)') 'Warning: invalid value of coord. ', i,', temporarily switched off'
                ENDIF
              ENDIF

              !c take care of PBC
              IF ((ICOORDINATES%COORDSTRUCT(i)%TAG=='T ') .OR. (ICOORDINATES%COORDSTRUCT(i)%TAG=='ST')) THEN  
                DO
                  IF (dq(j)>pi) THEN
                    !dq(j)=2*pi-dq(j)
                    dq(j)=dq(j)-2*pi
                  ELSE
                    EXIT
                  ENDIF
                ENDDO
                DO
                  IF (dq(j)<=-pi) THEN
                    dq(j)=2*pi+dq(j)
                  ELSE
                    EXIT
                  ENDIF
                ENDDO
              ENDIF
              
              ICOORDINATES%COORDSTRUCT(i)%VALUE=ICOORDINATES%COORDSTRUCT(i)%DVALUE+dq(j)

              !IF (IO%IU6>0) write(*,*) 'lq',ICOORDINATES%COORDSTRUCT(i)%VALUE,ICOORDINATES%COORDSTRUCT(i)%DVALUE
              !IF (IO%IU6>0) write(*,*) 'qq',dq(j)

              !IF ((IO%IU6>0)) THEN
                !write(g_io%REPORT,FMT='(/,2X,A15)') '>Harmonic_coord'
              !  write(g_io%REPORT,FMT='(3X,A3,X,A2,X,F15.5,X,F15.5,X,F15.5)') &
              !  &  'hc>',ICOORDINATES%COORDSTRUCT(i)%TAG,ICOORDINATES%COORDSTRUCT(i)%VALUE, ICOORDINATES%COORDSTRUCT(i)%DVALUE, dq(j)
              !ENDIF

            ENDIF
          ENDDO

          !DO i=1, ICOORDINATES%iconst1
          !  DO j=1, ICOORDINATES%iconst1
          !    IF (IO%IU6>0)  write(*,*) "gm",GMAT(i,j)
          !  ENDDO
          !ENDDO
           
          !DO i=1, ICOORDINATES%iconst3
          !  do j=1,T_INFO%NIONS
          !    IF (IO%IU6>0)  write(*,*) "bm",HARMONIC_BMAT(i,3*j-2),HARMONIC_BMAT(i,3*j-1),HARMONIC_BMAT(i,3*j)
          !  enddo
          !ENDDO
          !DO i=1, ICOORDINATES%iconst1
          !  do j=1,T_INFO%NIONS
          !    IF (IO%IU6>0)  write(*,*) "am",AMAT(3*j-2,i),AMAT(3*j-1,i),AMAT(3*j,i)
          !  enddo
          !ENDDO
           
          !GMAT=MATMUL(HARMONIC_BMAT,TRANSPOSE(HARMONIC_BMAT))
          !CALL SVDINVERSE(GMAT,ICOORDINATES%iconst1,idum)
          !AMAT=MATMUL(TRANSPOSE(HARMONIC_BMAT),GMAT)
            
          DMAT=MATMUL(HARMONIC_HESSIAN,AMAT)
          GMAT=MATMUL(TRANSPOSE(AMAT),DMAT)

          !IF (IO%IU6>0) THEN
          !      write(*,*) 'GMAT',GMAT
          !ENDIF
          
          !c force in internal coords
          df=MATMUL(GMAT,(dq))
          !df=0.
          !CALL DGEMV('N',ICOORDINATES%iconst1 ,ICOORDINATES%iconst1,1.,GMAT,3*T_INFO%NIONS,dq, 1, 0., df, 1)
          !c update energy
          HARMONIC_ENERGY=0.5*SUM(dq*df)+HARMONIC_ENERGY0
          
          !DMAT=MATMUL(TRANSPOSE(HARMONIC_BMAT),GMAT)
          
          !df2=MATMUL(DMAT,(dq))
          !HARMONIC_ENERGY=SUM(dq*df)+HARMONIC_ENERGY0
          
          IF (IO%IU6>0) THEN
            write(g_io%REPORT,FMT='(/,2X,A26)') '>Thermodynamic Integration'
            write(g_io%REPORT,FMT='(7X,3(X,A17))')      'E(lambda=0)','E(lambda=1)','E(1)-E(0)'
            write(g_io%REPORT,FMT='(3X,A5,3(X,E17.8E2))') 'e_ti>', HARMONIC_ENERGY,EPOT,EPOT-HARMONIC_ENERGY
           ! write(*,*) 'df',dq,df
          ENDIF
          
          !c force in Cartesian coordinates
          !dfc=-MATMUL(TRANSPOSE(HARMONIC_BMAT),df)
          dfc=-(MATMUL(df,HARMONIC_BMAT))
          !dfc=(MATMUL(df,HARMONIC_BMAT))
          CALL ONETOTHREE(T_INFO%NIONS,dummy,dfc)
          !dummy=MATMUL(TRANSPOSE(LATT_CUR%B),dummy)
          CALL KARDIR(T_INFO%NIONS,dummy,LATT_CUR%B)
          
          FACT=(DYN%POTIM**2)*EVTOJ/AMTOKG *1E-10_q
          NI=1
          DO NT=1,T_INFO%NTYP
            DO NI=NI,T_INFO%NITYP(NT)+NI-1
              IF (T_INFO%LSFOR(1,NI)) DYN%D2C(1,NI)=TILAMBDA*DYN%D2C(1,NI)+(1.-TILAMBDA)*dummy(1,NI)*FACT/2/T_INFO%POMASS(NT)
              IF (T_INFO%LSFOR(2,NI)) DYN%D2C(2,NI)=TILAMBDA*DYN%D2C(2,NI)+(1.-TILAMBDA)*dummy(2,NI)*FACT/2/T_INFO%POMASS(NT)
              IF (T_INFO%LSFOR(3,NI)) DYN%D2C(3,NI)=TILAMBDA*DYN%D2C(3,NI)+(1.-TILAMBDA)*dummy(3,NI)*FACT/2/T_INFO%POMASS(NT)
              !IF (IO%IU6>0) write(*,*) "hg",dummy(1,NI),dummy(2,NI),dummy(3,NI)
            ENDDO 
          ENDDO
          !!DYN%VEL=DYN%VEL+2*(DYN%D2C-hills_accel-penalty_accel-hspring_accel)
          !!spring_accel=MATMUL(LATT_CUR%B,dummy) !!!!!!!!!
          !!dummy=MATMUL(TRANSPOSE(LATT_CUR%B),dummy)
          
        END SUBROUTINE HARMONIC_FORCEFIELD_COMPUTE
        
        SUBROUTINE ENERGY_OUT(g_io,DYN,T_INFO,EKIN,EPOT,ECONST,ESPRING,EPS,ES,counter)
        !c write some output from MD
          TYPE(gadget_io) :: g_io
          TYPE(dynamics) :: DYN
          TYPE(type_info) :: T_INFO
          INTEGER :: counter
          REAL(q) :: ETOTAL,EKIN,EPOT,EPS,ES,ECONST,ESPRING

          ETOTAL=EKIN+EPOT+ECONST+EPS+ES+ESPRING
            

            write(g_io%REPORT,FMT='(/,2X,A9)') '>Energies'
            !write(g_io%REPORT,FMT='(7X,7(X,A17))')      'E_tot','E_pot','E_kin','E_spring','E_const','EPS','ES'
            !write(g_io%REPORT,FMT='(3X,A4,6(X,E17.8E2))') 'e_b>', ETOTAL,EPOT,EKIN,ESPRING,ECONST,EPS,ES
            write(g_io%REPORT,FMT='(7X,5(X,A17))')      'E_tot','E_pot','E_kin','EPS','ES'
            write(g_io%REPORT,FMT='(3X,A4,4(X,E17.8E2))') 'e_b>', ETOTAL,EPOT,EKIN,EPS,ES

            !write(g_io%REPORT,FMT='(/,2X,A11)') '>Positions:'
            !DO i=1,T_INFO%NIONS
            !  write(g_io%REPORT,FMT='(3X,A4, I5,3(X,F9.5))') 'x_b>', i,DYN%POSION(:,i)
            !ENDDO

            !write(g_io%REPORT,FMT='(/,2X,A12)') '>Velocities:'
            !DO i=1,T_INFO%NIONS
            !  write(g_io%REPORT,FMT='(3X,A4, I5,3(X,F9.5))') 'v_b>', i,DYN%VEL(:,i)
            !ENDDO

            !write(g_io%REPORT,FMT='(/,2X,A14)') '>Acceleration:'
            !DO i=1,T_INFO%NIONS
            !  write(g_io%REPORT,FMT='(3X,A4, I5,3(X,F9.5))') 'a_b>', i,DYN%D2C(:,i)
            !ENDDO

            !write(g_io%REPORT,FMT='(/,2X,A20)') '>Acceleration_const:'
            !DO i=1,T_INFO%NIONS
            !  write(g_io%REPORT,FMT='(3X,A4, I5,3(X,F9.5))') 'a_c>', i,C_FORCE(:,i)
            !ENDDO

        END SUBROUTINE ENERGY_OUT

        SUBROUTINE NOSE_OUT(g_io,DYN)
        !c write some output from MD
          TYPE(gadget_io) :: g_io
          TYPE(dynamics) :: DYN

            write(g_io%REPORT,FMT='(/,2X,A12)') '>Thermostat:'
            write(g_io%REPORT,FMT='(7X,5(X,A9))') &
                     "s(n+1)","sdot(n+1)","s(n)","sdot(n)","Qs"
            write(g_io%REPORT,FMT='(3X,A4,5(X,F9.5))') 's_b>',&
                     DYN%SNOSE(1),DYN%SNOSE(2),DYN%SNOSE(4),DYN%SNOSE(3),DYN%SMASS
        END SUBROUTINE NOSE_OUT

        SUBROUTINE TEMPERATURE_OUT(g_io,DYN,T_new)
        !c write some output from MD
          TYPE(gadget_io) :: g_io
          TYPE(dynamics) :: DYN
           REAL(q) :: T_new

          write(g_io%REPORT,FMT='(/,2X,A12)') '>Temperature'
          write(g_io%REPORT,FMT='(9X,2(X,A10))') &
                & "T_sim","T_inst"    
          write(g_io%REPORT,FMT='(3X,A6,2(X,F10.3))') &
                & 'tmprt>', DYN%TEMP,T_new    
        END SUBROUTINE TEMPERATURE_OUT


        SUBROUTINE CONSTRAINED_DOF(DYN,T_INFO,bdim,LATT_CUR,ICOORDINATES,cDOF)
         !c identify number of constrained (non-redundant) degrees of freedom
          TYPE(dynamics) :: DYN
          TYPE(type_info) :: T_INFO
          TYPE(latt) :: LATT_CUR
          TYPE(coordstructure) :: ICOORDINATES
          INTEGER :: i,bdim,cDOF
          REAL(q) :: BMAT(ICOORDINATES%NUMINTERNALS,bdim)
          REAL(q) :: norm

          cDOF=0
                  
          IF (ICOORDINATES%NUMINTERNALS>=1) THEN          
            BMAT=0._q
            CALL BMATRIX(T_INFO,DYN%POSION,DYN%POSION,TRANSPOSE(LATT_CUR%A),ICOORDINATES,BMAT,.TRUE.)
            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS/=0) THEN
                BMAT(i,:)=0._q
              ENDIF
            ENDDO
            IF (ICOORDINATES%NUMINTERNALS>1) CALL ORTHO_NORMALIZE(BMAT)
            DO i=1,ICOORDINATES%NUMINTERNALS
              norm=VECTORSIZE(bdim,BMAT(i,:))
              IF (norm>1e-5) THEN
                cDOF=cDOF+1 
              END IF
            ENDDO
          ENDIF 
        END SUBROUTINE CONSTRAINED_DOF

        SUBROUTINE STEP_ANDERSEN(DYN,ICOORDINATES,T_INFO,INFO,LATT_CUR, &
          &      EKIN,ESPRING,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,NDEGREES_OF_FREEDOM_,IO,&
          &      EPOT,ANDERSEN_PROB,g_io,iconst0,iconst2,TEIN)
        !c Andersen thermostat
          TYPE(coordstructure) :: ICOORDINATES     !c working and reference coords.
          TYPE(dynamics) :: DYN
          TYPE(type_info) :: T_INFO
          TYPE (info_struct) INFO
          TYPE(latt) :: LATT_CUR
          TYPE (in_struct) :: IO
          TYPE(gadget_io) :: g_io
          REAL(q) :: EKIN,EPS,ES,DISMAX,TEMPER,EPOT,ECONST,TEIN,ESPRING
          INTEGER :: i,j,ii,jj,NDEGREES_OF_FREEDOM,NDEGREES_OF_FREEDOM_
          INTEGER,SAVE :: counter=0
          INTEGER :: iconst0,iconst2
          REAL(q) :: CMASS(3,1)
          REAL(q) :: VEL_tmp(3,T_INFO%NIONS)
          REAL(q) :: C_FORCE(3,T_INFO%NIONS)
          REAL(q) :: dummy
          REAL(q) :: ANDERSEN_PROB
          INTEGER :: random_counter
          REAL(q) :: AMASS,C_FORCE_L(3,3)
          LOGICAL :: LSCALE,LCMASS

          LSCALE=.FALSE.
          LCMASS=.FALSE.
          !c add this step to total count
          counter=counter+1

          AMASS=0._q

          VEL_tmp=0._q
          CALL init_velocities(VEL_tmp,T_INFO,DYN,LATT_CUR,NDEGREES_OF_FREEDOM_,LSCALE,LCMASS)
          CMASS=0._q

          DYN%VEL=DYN%VEL+DYN%D2C
          DYN%VEL=DYN%VEL+DYN%D2C
          random_counter=0
          IF (DYN%INIT==1) THEN
            DO i=1,T_INFO%NIONS
              CALL RANE_ION(dummy)
              IF (dummy<ANDERSEN_PROB) THEN
                random_counter=random_counter+1
                DYN%VEL(1:3,i)=VEL_tmp(1:3,i)
              ENDIF
            ENDDO
          ENDIF

          CALL GIVE_CMASS(T_INFO,DYN%VEL,CMASS)
          DO i=1,T_INFO%NIONS
            DO j=1,3
              IF (T_INFO%LSFOR(j,i)) DYN%VEL(j,i)=DYN%VEL(j,i)-CMASS(j,1)
            ENDDO
          ENDDO
          !CALL GIVE_CMASS(T_INFO,DYN%VEL,CMASS)

          C_FORCE=0._q
          DYN%POSION=DYN%POSIOC+DYN%VEL

          ECONST=0._q
          IF (ICOORDINATES%NUMINTERNALS>0) THEN
            CALL Lproj_out_step(T_INFO,INFO,DYN,LATT_CUR,IO,C_FORCE,C_FORCE_L,ICOORDINATES,g_io,3*T_INFO%NIONS,AMASS,LATT_CUR%A,ECONST,iconst0,iconst2)
            DYN%POSION=DYN%POSION+C_FORCE
            DYN%VEL=DYN%VEL+C_FORCE
          ENDIF

          CALL KINETIC_E(EKIN,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,&
                       & DYN%POTIM,LATT_CUR%A,DYN%VEL)
          TEIN=2*EKIN/BOLKEV/NDEGREES_OF_FREEDOM

          IF (IO%IU6>0) THEN
            CALL ENERGY_OUT(g_io,DYN,T_INFO,EKIN,EPOT,ECONST,ESPRING,0._q,0._q,counter)
            CALL TEMPERATURE_OUT(g_io,DYN,TEIN)
            write(g_io%REPORT,FMT='(/,2X,A32,X,I16)') '>Thermostat, num. of collisions:',random_counter            
          ENDIF


          DISMAX=0
          DO I=1,T_INFO%NIONS
             DISMAX=MAX(DISMAX, &
             &(DYN%VEL(1,I))**2+ &
             &(DYN%VEL(2,I))**2+ &
             &(DYN%VEL(3,I))**2)
          ENDDO
          DISMAX=SQRT(DISMAX)
          EPS=0.0_q
          ES =0.0_q
        END SUBROUTINE STEP_ANDERSEN

        SUBROUTINE STEP_ANDERSEN_STOCH(DYN,ICOORDINATES,T_INFO,INFO,LATT_CUR, &
          &      EKIN,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,NDEGREES_OF_FREEDOM_,IO,&
          &      EPOT,ANDERSEN_PROB,g_io,iconst0,iconst2,TEIN)
        ! this version removes the center of mass drift from the final velocities
        ! in each step
        ! before randomizing the individual velocities
        ! a random center of mass drift is added 
        ! see below
        ! this routine can be used as an alternative to STEP_ANDERSEN
        ! the main advantage is that the routine is "stateless"

          TYPE(coordstructure) :: ICOORDINATES     !c working and reference coords.
          TYPE(dynamics) :: DYN
          TYPE(type_info) :: T_INFO
          TYPE (info_struct) INFO
          TYPE(latt) :: LATT_CUR
          TYPE (in_struct) :: IO
          TYPE(gadget_io) :: g_io
          REAL(q) :: EKIN,EPS,ES,DISMAX,TEMPER,EPOT,ECONST,TEIN
          INTEGER :: i,j,ii,jj,NDEGREES_OF_FREEDOM,NDEGREES_OF_FREEDOM_
          INTEGER,SAVE :: counter=0
          INTEGER :: iconst0,iconst2
          REAL(q) :: CMASS(3,1)
          REAL(q) :: VEL_tmp(3,T_INFO%NIONS)
          REAL(q) :: C_FORCE(3,T_INFO%NIONS)
          REAL(q) :: dummy
          REAL(q) :: ANDERSEN_PROB
          INTEGER :: random_counter
          REAL(q) :: AMASS,C_FORCE_L(3,3)
          REAL(q) :: TOTAL_MASS
          INTEGER,SAVE :: Ltxyz(3)

          !c add this step to total count
          counter=counter+1

          !c determine if coordinates for some direction
          !c are frozen
          IF (counter==1) THEN
            Ltxyz(:)=1
            DO j=1,3
              DO i=1,T_INFO%NIONS
                IF (.NOT. T_INFO%LSFOR(j,i)) Ltxyz(j)=0
              ENDDO
            ENDDO

            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0) THEN
                IF (ICOORDINATES%COORDSTRUCT(i)%TAG=="X ") Ltxyz(1)=0
                IF (ICOORDINATES%COORDSTRUCT(i)%TAG=="Y ") Ltxyz(2)=0
                IF (ICOORDINATES%COORDSTRUCT(i)%TAG=="Z ") Ltxyz(3)=0
              ENDIF
            ENDDO
          ENDIF

          AMASS=0._q

          VEL_tmp=0._q

          ! Add a random velicity that is exactly
          ! distributed like the centre of mass velocity
          TOTAL_MASS = 0.0_q
          DO I = 1, T_INFO%NIONS
            TOTAL_MASS = TOTAL_MASS + T_INFO%POMASS(T_INFO%ITYP(I))
          END DO
          CMASS(:,1) = BOLTZMANN_VELOCITY(DYN, LATT_CUR, TOTAL_MASS)

          !c don't add translation for frozen direction
          DO i=1,3
            IF (Ltxyz(i)==0) CMASS(i,1)=0._q
          ENDDO

          DO i=1,T_INFO%NIONS
            DO j=1,3
              DYN%VEL(j,i)=DYN%VEL(j,i)+CMASS(j,1)
            ENDDO
          ENDDO

          CALL boltzmann_velocities(T_INFO,DYN,LATT_CUR,VEL_tmp)

          DYN%VEL=DYN%VEL+DYN%D2C
          DYN%VEL=DYN%VEL+DYN%D2C
          random_counter=0
          IF (DYN%INIT==1) THEN
            DO i=1,T_INFO%NIONS
              CALL RANE_ION(dummy)
              IF (dummy<ANDERSEN_PROB) THEN
                random_counter=random_counter+1
                DYN%VEL(1:3,i)=VEL_tmp(1:3,i)
              ENDIF
            ENDDO
          ENDIF

          CMASS=0._q
          !c GIVE_CMASS takes care of frozen directions (corresponding
          !c component set to zero)
          CALL GIVE_CMASS(T_INFO,DYN%VEL,CMASS)

          ! store the centre of mass drift
          ! so that we can add it again before swaping numbers
          DO i=1,T_INFO%NIONS
            DO j=1,3
              DYN%VEL(j,i)=DYN%VEL(j,i)-CMASS(j,1)
            ENDDO
         ENDDO
          !CALL GIVE_CMASS(T_INFO,DYN%VEL,CMASS)

         !c velocity for the step t-1/2dt
          VEL_tmp=DYN%VEL

          !c velocity for the step t+1/2dt
          !DYN%VEL=DYN%VEL+DYN%D2C
          !DYN%VEL=DYN%VEL+DYN%D2C
          
          DYN%POSION=DYN%POSIOC+DYN%VEL

          ECONST=0._q
          IF (iconst0>0 .OR. iconst2>0) THEN
            C_FORCE=0._q
            CALL Lproj_out_step(T_INFO,INFO,DYN,LATT_CUR,IO,C_FORCE,C_FORCE_L,ICOORDINATES,g_io,3*T_INFO%NIONS,AMASS,LATT_CUR%A,ECONST,iconst0,iconst2)
            DYN%POSION=DYN%POSION+C_FORCE
            DYN%VEL=DYN%VEL+C_FORCE
          ENDIF

          !!c velocity for the step t
          !!VEL_tmp=(VEL_tmp+DYN%VEL)/2
          !!CALL KINETIC_E(EKIN,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,&
          !!             & DYN%POTIM,LATT_CUR%A,VEL_tmp)

          !c due to compatibility reasons we still compute
          !c kinetic energy and T using velocities from the step t+1/2
          !c but in the future we should use VEL_tmp (velocity for step t) instead
          CALL KINETIC_E(EKIN,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,&
                       & DYN%POTIM,LATT_CUR%A,DYN%VEL)
          TEIN=2*EKIN/BOLKEV/NDEGREES_OF_FREEDOM

          IF (IO%IU6>0) THEN
            CALL ENERGY_OUT(g_io,DYN,T_INFO,EKIN,EPOT,ECONST,0._q,0._q,0._q,counter)
            CALL TEMPERATURE_OUT(g_io,DYN,TEIN)
            write(g_io%REPORT,FMT='(/,2X,A32,X,I16)') '>Thermostat, num. of collisions:',random_counter
          ENDIF


          DISMAX=0
          DO I=1,T_INFO%NIONS
             DISMAX=MAX(DISMAX, &
             &(DYN%VEL(1,I))**2+ &
             &(DYN%VEL(2,I))**2+ &
             &(DYN%VEL(3,I))**2)
          ENDDO
          DISMAX=SQRT(DISMAX)
          EPS=0.0_q
          ES =0.0_q
        END SUBROUTINE STEP_ANDERSEN_STOCH


        SUBROUTINE STEP_LANGEVIN(DYN,ICOORDINATES,hills,penalty,T_INFO,INFO,LATT_CUR, &
          &      EKIN,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,NDEGREES_OF_FREEDOM_,IO,&
          &      EPOT,GAMM,g_io,iconst0,iconst2,iconst4,iconst5,iconst6,iconst8,& 
          &      hspring_R0,hspring_K,hspring_force,fbias_R0,fbias_D,fbias_A,TEIN)
        !c Langevin thermostat, velocity verlet algorithm!!!
          TYPE(coordstructure) :: ICOORDINATES     !c working and reference coords.
          TYPE(dynamics) :: DYN
          TYPE(type_info) :: T_INFO
          TYPE (info_struct) INFO
          TYPE(latt) :: LATT_CUR
          TYPE (in_struct) :: IO
          TYPE(gadget_io) :: g_io
          TYPE(penalty_data) :: penalty
          TYPE(hills_data) :: hills
          REAL(q) :: EKIN,EPS,ES,DISMAX,TEMPER,EPOT,ECONST
          INTEGER :: i,j,ii,jj,NDEGREES_OF_FREEDOM,NDEGREES_OF_FREEDOM_
          INTEGER,SAVE :: counter=0
          INTEGER :: iconst0,iconst2,iconst4,iconst5,iconst6,iconst8
          REAL(q) :: CMASS1(3,1),CMASS2(3,1)
          REAL(q) :: VEL_tmp(3,T_INFO%NIONS)
          REAL(q) :: V_FORCE(3,T_INFO%NIONS),C_FORCE(3,T_INFO%NIONS)
          REAL(q) :: dummy
          INTEGER :: random_counter
          REAL(q) :: AMASS,TEIN
          REAL(q) :: C_FORCE_L(3,3),V_FORCE_L(3,3)
          REAL(q) :: GAMM(T_INFO%NTYP) !c friction coeficient (in ps^(-1)) for each atomic type
          REAL(q) :: D2C_stoch(3,T_INFO%NIONS),D2C_momentum(3,T_INFO%NIONS)
          REAL(q) :: ekin1,ekin2,tein1,tein2
          REAL(q),PARAMETER :: tol=1e-6
          INTEGER,PARAMETER :: MAXITER=1000
          REAL(q) :: err2
          REAL(q) :: vel_(3,T_INFO%NIONS), vel_old(3,T_INFO%NIONS)
          REAL(q) :: hills_accel(3,T_INFO%NIONS),penalty_accel(3,T_INFO%NIONS),hspring_accel(3,T_INFO%NIONS)
          REAL(q) :: hills_accel_L(3,3),penalty_accel_L(3,3),hspring_accel_L(3,3)
          INTEGER :: NI, NT
          REAL(q) :: FACT,VTMP(3)
          REAL(q) :: maxForce
          LOGICAL :: LSTOP2
          REAL(q),SAVE :: EPOT_old
          REAL(q) :: hspring_R0(iconst8),hspring_K(iconst8),hspring_force(iconst8)
          REAL(q) :: fbias_R0(iconst4),fbias_D(iconst4),fbias_A(iconst4),fbias_force(iconst4)
          REAL(q) :: fbias_accel(3,T_INFO%NIONS),fbias_accel_L(3,3)
          
         
 72     FORMAT( ' POSITION    ',35X,'REST-FORCE (eV/Angst)'/ &
     &          ' ----------------------------------------------', &
     &          '-------------------------------------')
 76     FORMAT((3F13.5,3X,3F14.6))


          !c add this step to total count
          counter=counter+1

          !c make sure that energy convergence criterion is not fullfiled 
          !c in the first relaxation step (irelevant for MD)
          IF (counter==1) EPOT_old=EPOT+10*ABS(DYN%EDIFFG)

          V_FORCE_L=0._q
          C_FORCE_L=0._q
          V_FORCE=0._q
          C_FORCE=0._q
          hills_accel_L=0._q;penalty_accel_L=0._q;hspring_accel_L=0._q
          hills_accel=0._q;penalty_accel=0._q;hspring_accel=0._q
          fbias_accel=0._q;fbias_accel_L=0._q

      !c metadynamics stuff 
          !c update collective variables
          IF (iconst5 >0) THEN
            j=0
            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==5) THEN
                j=j+1
                hills%position(j)=ICOORDINATES%COORDSTRUCT(i)%VALUE 
              ENDIF
            ENDDO
          
            !c add new hill
            IF (counter .GE. hills%bin) THEN
              CALL hills_bias_potential(hills,penalty,ICOORDINATES,DYN,counter,g_io,IO)
            ENDIF 

            !c calculate the penalty potential
            CALL penalty_potential(penalty,hills,ICOORDINATES)
 
            !c compute contributions from bias potentials
            CALL hills_bias_direct(T_INFO,LATT_CUR,DYN,hills%force,ICOORDINATES,DYN%POSIOC,hills_accel,hills_accel_L,3*T_INFO%NIONS,0._q)
            CALL hills_bias_direct(T_INFO,LATT_CUR,DYN,penalty%force,ICOORDINATES,DYN%POSIOC,penalty_accel,penalty_accel_L,3*T_INFO%NIONS,0._q)
            !hills%force=hills%force+penalty%force
            !CALL hills_bias_direct(T_INFO,LATT_CUR,DYN,hills%force,hspring_force,ICOORDINATES,DYN%POSIOC,hills_accel,hills_accel_L,3*T_INFO%NIONS,0._q)
          ENDIF
          
          !c calculate forces due to Fermi-shaped bias potential
          IF (iconst4>0) THEN
            fbias_force=0.
            CALL fermi_bias_force(iconst4,fbias_r0,fbias_d,fbias_a,fbias_force,ICOORDINATES)
            CALL harmonic_spring_accel(4,T_INFO,LATT_CUR,DYN,fbias_force,ICOORDINATES,DYN%POSIOC,fbias_accel,fbias_accel_L,3*T_INFO%NIONS,0._q)
          ENDIF

          !c calculate forces due to harmonic springs
          IF (iconst8>0) THEN
            CALL harmonic_spring_force(iconst8,hspring_R0,hspring_K,hspring_force,ICOORDINATES)
            CALL harmonic_spring_accel(8,T_INFO,LATT_CUR,DYN,hspring_force,ICOORDINATES,DYN%POSIOC,hspring_accel,hspring_accel_L,3*T_INFO%NIONS,0._q)
          ENDIF

          IF ((iconst4>0) .OR. (iconst8>0)) hspring_accel=hspring_accel+fbias_accel

          !c compute accel. due to friction term
          CALL friction_Forces_stoch(D2C_stoch,T_INFO,DYN,LATT_CUR,GAMM)

          !c v(t)dt=v´(t)dt+f(t)/(2m)*dt^2
          DYN%VEL=DYN%VEL+DYN%D2C+D2C_stoch

          !IF (iconst5 .GT. 0) DYN%VEL=DYN%VEL-(hills_accel+penalty_accel)
          IF ((iconst5 .GT. 0) .OR. (iconst8 .GT. 0) .OR. (iconst4 .GT. 0)) DYN%VEL=DYN%VEL-(hills_accel+penalty_accel+hspring_accel)
          !c factor of 1/2 was missing - CHECK!!!
          !IF ((iconst5 .GT. 0) .OR. (iconst8 .GT. 0)) DYN%VEL=DYN%VEL-(hills_accel+penalty_accel+hspring_accel)/2

   
          vel_=DYN%VEL

          !c velocities for the time slice (t) depend on themselves
          !c iterative solution is needed
          i=0
          DO
            i=i+1
            IF (i .GT. MAXITER) THEN
              CALL vtutor%error("Error: convergence problem!")
            ENDIF 
            vel_old=DYN%VEL
  
            !c compute accel. due to friction term 
            CALL friction_Forces_momentum(D2C_momentum,T_INFO,DYN,GAMM)
   
            DYN%VEL=vel_+D2C_momentum 

            IF (ICOORDINATES%NUMINTERNALS>0) THEN                            
               CALL Rattle(counter,T_INFO,INFO,DYN,LATT_CUR,IO,C_FORCE,C_FORCE_L,V_FORCE,V_FORCE_L,ICOORDINATES,g_io,3*T_INFO%NIONS,&
               &   AMASS,LATT_CUR%A,ECONST,iconst0,iconst2,1)
              !CALL Rattle_vel(T_INFO,INFO,DYN,LATT_CUR,IO,C_FORCE,C_FORCE_L,V_FORCE,V_FORCE_L,ICOORDINATES,g_io,3*T_INFO%NIONS,AMASS,LATT_CUR%A,ECONST,iconst0,iconst2)
              DYN%VEL=DYN%VEL+V_FORCE
            ENDIF
                            
            vel_old=ABS(vel_old-DYN%VEL)
            err2=MAXVAL(vel_old)
            IF (err2<tol) EXIT
          ENDDO

          !DYN%VEL=DYN%VEL+DYN%D2C+D2C_momentum+D2C_stoch
          CALL langevin_temperature(EKIN1,EKIN2,TEIN1,TEIN2,T_INFO,DYN,LATT_CUR%A,GAMM)
          
          CALL KINETIC_E(EKIN,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,&
                       & DYN%POTIM,LATT_CUR%A,DYN%VEL)
          TEIN=2*EKIN/BOLKEV/NDEGREES_OF_FREEDOM

!           IF (IO%IU6>0) THEN
!             CALL ENERGY_OUT(g_io,DYN,T_INFO,EKIN,EPOT,ECONST,0._q,0._q,counter)
!             CALL TEMPERATURE_OUT(g_io,DYN,TEIN) 
!             IF (IO%NWRITE==3) THEN
!               write(g_io%REPORT,FMT='(9X,2(X,A10))') &
!                  & "T(Lang.)","T(nLang.)"    
!               write(g_io%REPORT,FMT='(3X,A6,2(X,F10.3))') &
!                  & 'tmpLg>',  TEIN1,TEIN2   
!             ENDIF      
!           ENDIF

          !c v'(t+dt)*dt=v(t)*dt+f(t)/(2m)*dt^2
          DYN%VEL=DYN%VEL+DYN%D2C+D2C_momentum+D2C_stoch

          !c metydynamics or biased MD
          !IF (iconst5 .GT. 0 ) DYN%VEL=DYN%VEL-(hills_accel+penalty_accel)
          IF ((iconst5 .GT. 0) .OR. (iconst8 .GT. 0) .OR. (iconst4 .GT. 0)) DYN%VEL=DYN%VEL-(hills_accel+penalty_accel+hspring_accel)
          !c factor of 1/2 was missing - CHECK!!!
          !IF ((iconst5 .GT. 0) .OR. (iconst8 .GT. 0)) DYN%VEL=DYN%VEL-(hills_accel+penalty_accel+hspring_accel)/2


!           IF (ICOORDINATES%NUMINTERNALS>0) THEN           
!            DYN%VEL=DYN%VEL+V_FORCE/2._q
!           ENDIF

          !c compute position of c. of mass for the slice t:
          CALL GIVE_CMASS(T_INFO,DYN%POSIOC,CMASS1)

          !c r(t+dt)=r(t)+v'(t+dt)*dt
          DYN%POSION=DYN%POSIOC+DYN%VEL

          
          ECONST=0._q
          IF (ICOORDINATES%NUMINTERNALS>0) THEN
            !c and now RATTLE
            CALL Rattle(counter,T_INFO,INFO,DYN,LATT_CUR,IO,C_FORCE,C_FORCE_L,V_FORCE,V_FORCE_L,ICOORDINATES,g_io,3*T_INFO%NIONS,&
            &   AMASS,LATT_CUR%A,ECONST,iconst0,iconst2,0)
            DYN%POSION=DYN%POSION+C_FORCE
            DYN%VEL=DYN%VEL+C_FORCE
          ENDIF

          !c if this MD subroutine is used for relaxation
          IF (DYN%TEBEG<0.1) THEN
            !c compute the total restoring force
            IF (iconst0>0) THEN
              C_FORCE=C_FORCE+V_FORCE
              V_FORCE=0._q

              V_FORCE=DYN%D2C+C_FORCE

              FACT=(DYN%POTIM**2)*EVTOJ/AMTOKG *1E-10_q
              NI=1
              DO NT=1,T_INFO%NTYP
                DO NI=NI,T_INFO%NITYP(NT)+NI-1
                  V_FORCE(1,NI)=V_FORCE(1,NI)/FACT*2*T_INFO%POMASS(NT)
                  V_FORCE(2,NI)=V_FORCE(2,NI)/FACT*2*T_INFO%POMASS(NT)
                  V_FORCE(3,NI)=V_FORCE(3,NI)/FACT*2*T_INFO%POMASS(NT)
                ENDDO
              ENDDO

              CALL DIRKAR(T_INFO%NIONS,V_FORCE,LATT_CUR%A)

              !c write ionic forces after removing the restoring forces
              IF (IO%IU0>0) THEN
                WRITE(IO%IU6,72)
                DO J=1,T_INFO%NIONS
                  VTMP=T_INFO%POSION(1:3,J)
                  CALL  DIRKAR(1,VTMP,LATT_CUR%A)
                  WRITE(IO%IU6,76) VTMP,(V_FORCE(I,J),I=1,3)
                ENDDO
              ENDIF 
            ENDIF
            
            !c test the convergence criterion
            maxForce=0._q
            DO i=1,T_INFO%NIONS
              dummy=V_FORCE(1,i)*V_FORCE(1,i)+V_FORCE(2,i)*V_FORCE(2,i)+V_FORCE(3,i)*V_FORCE(3,i)
              IF (dummy>maxForce) &
                & maxForce=dummy
            END DO
            maxForce=maxForce**0.5
            IF (IO%IU0>0) THEN
              write(*,*) 'maxForce',maxForce,DYN%EDIFFG
            ENDIF 

            IF (DYN%EDIFFG<0.) THEN
              LSTOP2=.TRUE.

              IF (maxForce>ABS(DYN%EDIFFG)) LSTOP2=.FALSE.

              IF (LSTOP2) THEN
                INFO%LSTOP=.TRUE.
                !c take care of this POSION/POSIOC issue
                DYN%POSIOC=DYN%POSION
                RETURN
              END IF
            ELSE IF (counter>1) THEN
              IF (IO%IU0>0) THEN
              write(*,*) 'dEPOT',ABS(EPOT-EPOT_old),DYN%EDIFFG
            ENDIF
              INFO%LSTOP=(ABS(EPOT-EPOT_old)<DYN%EDIFFG)
            ENDIF
          ENDIF

          !!c compute position of c. of mass for the slice t+dt:
          !CALL GIVE_CMASS(T_INFO,DYN%POSION,CMASS2)

          !remove shift of c. of mass:
!           DO i=1,T_INFO%NIONS
!             DO j=1,3
!               IF (T_INFO%LSFOR(j,i)) DYN%POSION(j,i)=DYN%POSION(j,i)-(CMASS2(j)-CMASS1(j))
!             ENDDO
!           ENDDO
 
          IF (IO%IU6>0) THEN
            CALL ENERGY_OUT(g_io,DYN,T_INFO,EKIN,EPOT,ECONST,0._q,0._q,0._q,counter)
            CALL TEMPERATURE_OUT(g_io,DYN,TEIN) 
            IF (IO%NWRITE==3) THEN
              write(g_io%REPORT,FMT='(9X,2(X,A10))') &
                 & "T(Lang.)","T(nLang.)"    
              write(g_io%REPORT,FMT='(3X,A6,2(X,F10.3))') &
                 & 'tmpLg>',  TEIN1,TEIN2   
            ENDIF   

            !c output for metadynamics
            IF (iconst5>0) THEN
              write(g_io%REPORT,ADVANCE='YES',FMT='(/,2X,A13)') '>Metadynamics'
              DO i=1,ICOORDINATES%NUMINTERNALS
                IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==5) THEN
                  write(g_io%REPORT,ADVANCE='YES',FMT='(3X,A6,X,F9.5)') &
                  & 'fic_p>',ICOORDINATES%COORDSTRUCT(i)%VALUE
                ENDIF
              ENDDO   
            ENDIF    
          ENDIF

          DISMAX=0
          DO I=1,T_INFO%NIONS
             DISMAX=MAX(DISMAX, &
             &(DYN%VEL(1,I))**2+ &
             &(DYN%VEL(2,I))**2+ &
             &(DYN%VEL(3,I))**2)
          ENDDO
          DISMAX=SQRT(DISMAX)
          EPS=0.0_q
          ES =0.0_q

          EPOT_old=EPOT
        END SUBROUTINE STEP_LANGEVIN

        SUBROUTINE STEP_LANGEVIN_ISIF3(DYN,ICOORDINATES,hills,penalty,T_INFO,INFO,LATT_CUR, &
          &      EKIN,EKIN_Lat,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,NDEGREES_OF_FREEDOM_,IO,&
          &      EPOT,GAMM,GAMMA_L,g_io,SIF,iconst0,iconst2,iconst4,iconst5,iconst8,& 
          &      hspring_R0,hspring_K,hspring_force,fbias_R0,fbias_D,fbias_A,TEIN,AMASS,LATTICE_CONSTRAINTS,LATTICE_DOF)
        !c Langevin thermostat, velocity verlet algorithm
          USE reader_tags

          TYPE(coordstructure) :: ICOORDINATES     !c working and reference coords.
          TYPE(dynamics) :: DYN
          TYPE(type_info) :: T_INFO
          TYPE (info_struct) INFO
          TYPE(latt) :: LATT_CUR
          TYPE (in_struct) :: IO
          TYPE(gadget_io) :: g_io
          TYPE(penalty_data) :: penalty
          TYPE(hills_data) :: hills
          REAL(q) :: EKIN,EKIN_Lat,EPS,ES,DISMAX,TEMPER,EPOT,ECONST
          INTEGER :: i,j,ii,jj,NDEGREES_OF_FREEDOM,NDEGREES_OF_FREEDOM_
          INTEGER,SAVE :: counter=0
          INTEGER :: iconst0,iconst2,iconst4,iconst5,iconst8
          REAL(q) :: CMASS(3,1)
          REAL(q) :: V_FORCE(3,T_INFO%NIONS),C_FORCE(3,T_INFO%NIONS)
          REAL(q) :: GAMM(T_INFO%NTYP) !c friction coeficient (in ps^(-1)) for each atomic type
          REAL(q) :: GAMMA_L !c friction coeficient (in ps^(-1)) for lattice components
          REAL(q) :: D2C_stoch(3,T_INFO%NIONS),D2C_momentum(3,T_INFO%NIONS)
          REAL(q) :: AC(3,3),Aaccel(3,3),Avel_tmp(3,3),Ginv(3,3),Gdot(3,3)
          REAL(q) :: AMASS !,Avel(3,3)
          REAL(q) :: dumM(9),FACT
          REAL(q) :: SIF(3,3) 
          REAL(q) :: TEIN,TEIN_at,TEIN_lat,EKIN_at
          REAL(q) :: V_FORCE_L(3,3),C_FORCE_L(3,3)
          REAL(q),SAVE :: ROTMAT1(3,3)
          REAL(q) :: ROTMAT2(3,3)
          REAL(q) :: kinetic_stress(3,3),SIF_momentum(3,3),SIF_stoch(3,3)
          INTEGER :: IDUM,IERR,N
          LOGICAL :: LDUM,LOPEN
          COMPLEX(q) :: CDUM
          REAL(q) :: RDUM
          CHARACTER*1 :: CHARAC
          REAL(q),PARAMETER :: tol=1e-6
          INTEGER,PARAMETER :: MAXITER=100
          REAL(q) :: err1,err2
          REAL(q) :: Avel_(3,3),Avel_old(3,3)
          REAL(q) :: vel_(3,T_INFO%NIONS), vel_old(3,T_INFO%NIONS)
          LOGICAL :: LATTICE_CONSTRAINTS(3) ! specifies which of the lattice parameters (a_1,a_2,a_3) are free (true) or fixed (false)
          INTEGER :: LATTICE_DOF ! the total number of the lattice degrees of freedom according to the applied lattice constraints
          REAL(q) :: hills_accel(3,T_INFO%NIONS),penalty_accel(3,T_INFO%NIONS),hspring_accel(3,T_INFO%NIONS)
          REAL(q) :: hills_accel_L(3,3),penalty_accel_L(3,3),hspring_accel_L(3,3)
          INTEGER :: NI, NT
          REAL(q) :: VTMP(3)
          REAL(q) :: dummy, maxForce
          LOGICAL :: LSTOP2
          REAL(q),SAVE :: EPOT_old
          REAL(q) :: hspring_R0(iconst8),hspring_K(iconst8),hspring_force(iconst8)
          REAL(q) :: fbias_R0(iconst4),fbias_D(iconst4),fbias_A(iconst4),fbias_force(iconst4)
          REAL(q) :: fbias_accel(3,T_INFO%NIONS),fbias_accel_L(3,3)

 72     FORMAT( ' POSITION    ',35X,'REST-FORCE (eV/Angst)'/ &
     &          ' ----------------------------------------------', &
     &          '-------------------------------------')
 76     FORMAT((3F13.5,3X,3F14.6))

          !c add this step to total count
          counter=counter+1

          V_FORCE_L=0._q
          C_FORCE_L=0._q
          V_FORCE=0._q
          C_FORCE=0._q
          hills_accel_L=0._q;penalty_accel_L=0._q;hspring_accel_L=0._q
          hills_accel=0._q;penalty_accel=0._q;hspring_accel=0._q
          fbias_accel=0._q;fbias_accel_L=0._q

          !c first rotate cell with cell vectors to triangular matrix
          CALL ROTATE_CELL(IO%IU6,LATT_CUR,ROTMAT1)

          !c transform the triangular matrix with lattice vectors back to the original form 
          CALL ROTATE_BACK(LATT_CUR,ROTMAT1)

          IF (counter==1) THEN
            !c make sure that energy convergence criterion is not fullfiled 
            !c in the first relaxation step (irelevant for MD)
            IF (counter==1) EPOT_old=EPOT+10*ABS(DYN%EDIFFG)

            ! read-in mass for lattice parameters
            !CALL OPEN_INCAR_IF_FOUND(IO%IU5, LOPEN)
           
            !CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'PMASS', AMASS, IERR, WRITEXMLINCAR)
            !IF (IERR==3) Amass=1000._q

            !LATTICE_CONSTRAINTS = .TRUE.
            !CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LATTICE_CONSTRAINTS', LATTICE_CONSTRAINTS, 3, IERR, WRITEXMLINCAR)
            !IF (IERR==3) LATTICE_CONSTRAINTS = .TRUE.

            !CALL CLOSE_INCAR_IF_FOUND(IO%IU5)

            ! count the lattice degrees of freedom
            !LATTICE_DOF = 0
            !DO I = 1, 3
            !  IF (LATTICE_CONSTRAINTS(I)) LATTICE_DOF = LATTICE_DOF + 1
            !END DO
            !LATTICE_DOF = LATTICE_DOF * LATTICE_DOF

            !IF (IO%IU6>0) THEN
            !  write(g_io%REPORT,FMT='(/,2X,A35)') '>Damping parameter for lattice DOF:'
            !  write(g_io%REPORT,FMT='(3X,A22,X,F10.5)')    '              PMASS = ', AMASS  
            !  write(g_io%REPORT,FMT='(3X,A22,X,3L,A,I1,A)') 'LATTICE_CONSTRAINTS = ', (LATTICE_CONSTRAINTS(I), I=1,3), &
            !  & ' (X,Y,Z) = ', LATTICE_DOF, ' degree(s) of freedom'
            !  write(g_io%REPORT,FMT='(/)')
            !ENDIF 

            !c initialize the lattice velocities if needed 
            IF (LATT_CUR%INITlatv<0) THEN
              IF (IO%IU6>0) write(g_io%REPORT,FMT='(/,2X,A41)')  "Lattice velocities initialized by STEP_tb"
              CALL init_Avelocities(LATT_CUR%Avel,DYN,Amass)
              CALL APPLY_STRESS_CONSTRAINT(LATT_CUR%AVEL, Aaccel, LATTICE_CONSTRAINTS)
              LATT_CUR%INITlatv=1
            ENDIF
          ENDIF


       !c metadynamics stuff 
           !c update history of fictioous particles and of collective variables
           IF (iconst5 >0) THEN
             j=0
             DO i=1,ICOORDINATES%NUMINTERNALS
               IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==5) THEN
                 j=j+1
                 hills%position(j)=ICOORDINATES%COORDSTRUCT(i)%VALUE 
               ENDIF
             ENDDO
          
             !c add new hill
             IF (counter .GE. hills%bin) THEN
               CALL hills_bias_potential(hills,penalty,ICOORDINATES,DYN,counter,g_io,IO)
             ENDIF 

             !c calculate the penalty potential
             CALL penalty_potential(penalty,hills,ICOORDINATES)
 
             !c compute contributions from bias potentials
             CALL hills_bias_direct(T_INFO,LATT_CUR,DYN,hills%force,ICOORDINATES,DYN%POSIOC,hills_accel,hills_accel_L,3*T_INFO%NIONS+9,AMASS)
             CALL hills_bias_direct(T_INFO,LATT_CUR,DYN,penalty%force,ICOORDINATES,DYN%POSIOC,penalty_accel,penalty_accel_L,3*T_INFO%NIONS+9,AMASS)
             !hills%force=hills%force+penalty%force
             !CALL hills_bias_direct(T_INFO,LATT_CUR,DYN,hills%force,hspring_force,ICOORDINATES,DYN%POSIOC,hills_accel,hills_accel_L,3*T_INFO%NIONS+9,AMASS)
           ENDIF

          !c calculate forces due to Fermi-shaped bias potential
          IF (iconst4>0) THEN
            fbias_force=0.
            CALL fermi_bias_force(iconst4,fbias_r0,fbias_d,fbias_a,fbias_force,ICOORDINATES)
            CALL harmonic_spring_accel(4,T_INFO,LATT_CUR,DYN,fbias_force,ICOORDINATES,DYN%POSIOC,fbias_accel,fbias_accel_L,3*T_INFO%NIONS+9,AMASS)
          ENDIF  

          !c calculate forces due to harmonic springs
          IF (iconst8>0) THEN
            CALL harmonic_spring_force(iconst8,hspring_R0,hspring_K,hspring_force,ICOORDINATES)
            CALL harmonic_spring_accel(8,T_INFO,LATT_CUR,DYN,hspring_force,ICOORDINATES,DYN%POSIOC,hspring_accel,hspring_accel_L,3*T_INFO%NIONS+9,AMASS)
          ENDIF
 
          IF ((iconst4>0) .OR. (iconst8>0)) THEN
            hspring_accel=hspring_accel+fbias_accel
            hspring_accel_L=hspring_accel_L+fbias_accel_L
          ENDIF

          !!!DYN%POSIOC=DYN%POSION

          !c compute internal stress for the current time slice
          Aaccel=SIF !!+kinetic_stress

          !c take care of external pressure
          DO i=1,3
            Aaccel(i,i)=Aaccel(i,i)-DYN%PSTRESS/(EVTOJ*1E22_q)*LATT_CUR%OMEGA
          ENDDO
          
          !c compute forces
          Aaccel=MATMUL(Aaccel,(LATT_CUR%B))
         
          !c compute (one half of)acceleration
          FACT=(DYN%POTIM**2)*EVTOJ/AMTOKG *1E-10_q
          Aaccel=Aaccel*FACT/2/Amass 
       
          !c compute stochastic forces
          CALL friction_Forces_Lstoch(SIF_stoch,DYN,GAMMA_L,AMASS)
          CALL friction_Forces_stoch(D2C_stoch,T_INFO,DYN,LATT_CUR,GAMM)
          Aaccel=Aaccel+SIF_stoch


          ! update partly the velocities
          LATT_CUR%Avel=LATT_CUR%Avel+Aaccel
          ! set Avel here
          CALL APPLY_STRESS_CONSTRAINT(LATT_CUR%AVEL, Aaccel, LATTICE_CONSTRAINTS)
          DYN%VEL=DYN%VEL+DYN%D2C+D2C_stoch

!           IF (iconst5 .GT. 0) THEN
!             LATT_CUR%Avel=LATT_CUR%Avel-(hills_accel_L+penalty_accel_L)
!             DYN%VEL=DYN%VEL-(hills_accel+penalty_accel)
!           ENDIF
          IF ((iconst5 .GT. 0) .OR. (iconst8 .GT. 0) .OR. (iconst4 .GT. 0)) THEN
            LATT_CUR%Avel=LATT_CUR%Avel-(hills_accel_L+penalty_accel_L+hspring_accel_L)
            DYN%VEL=DYN%VEL-(hills_accel+penalty_accel+hspring_accel)
          ENDIF

          Avel_=LATT_CUR%Avel
          vel_=DYN%VEL
          
          !c the line below was incorrect!!
          !!Ginv=MATMUL(LATT_CUR%B,TRANSPOSE(LATT_CUR%B))
          Ginv=MATMUL(TRANSPOSE(LATT_CUR%B),LATT_CUR%B)

          !c velocities for time slice (t) depend on v(t)
          !c - iterative solution is needed
          i=0
          DO
            i=i+1
            IF (i .GT. MAXITER) THEN
              CALL vtutor%error("Error: convergence problem in MD!")
            ENDIF 
            Avel_old=LATT_CUR%Avel
            vel_old=DYN%VEL
            !c the line below was incorrect!!
            !!Gdot=MATMUL(LATT_CUR%Avel,TRANSPOSE(LATT_CUR%A))
            Gdot=MATMUL(TRANSPOSE(LATT_CUR%Avel),(LATT_CUR%A))
            Gdot=0.5*(Gdot+TRANSPOSE(Gdot))
            Gdot=MATMUL(Ginv,Gdot)

            !c compute kinetic stress
            CALL COMPUTE_KINETIC_STRESS(T_INFO,LATT_CUR,DYN,kinetic_stress,SIF,-1)
            kinetic_stress=MATMUL(kinetic_stress,(LATT_CUR%B))
            FACT=(DYN%POTIM**2)*EVTOJ/AMTOKG *1E-10_q
            kinetic_stress=kinetic_stress*FACT/2/Amass 

            !c compute accel. due to friction term 
            CALL friction_Forces_Lmomentum(SIF_momentum,LATT_CUR%Avel,DYN,GAMMA_L)
            CALL friction_Forces_momentum(D2C_momentum,T_INFO,DYN,GAMM)

            LATT_CUR%Avel=Avel_+SIF_momentum+kinetic_stress

            ! set Avel here
            CALL APPLY_STRESS_CONSTRAINT(LATT_CUR%AVEL, Aaccel, LATTICE_CONSTRAINTS)

            DYN%VEL=vel_+D2C_momentum-MATMUL(Gdot,DYN%VEL)

            !c geometric constraints
            IF (ICOORDINATES%NUMINTERNALS>0) THEN                            
               CALL Rattle(counter,T_INFO,INFO,DYN,LATT_CUR,IO,C_FORCE,C_FORCE_L,V_FORCE,V_FORCE_L,ICOORDINATES,g_io,3*T_INFO%NIONS+9,&
               &   AMASS,LATT_CUR%A,ECONST,iconst0,iconst2,1)
              !CALL Rattle_vel(T_INFO,INFO,DYN,LATT_CUR,IO,C_FORCE,C_FORCE_L,V_FORCE,V_FORCE_L,ICOORDINATES,g_io,3*T_INFO%NIONS,AMASS,LATT_CUR%A,ECONST,iconst0,iconst2)
              DYN%VEL=DYN%VEL+V_FORCE
              LATT_CUR%Avel=LATT_CUR%Avel+V_FORCE_L
            ENDIF

            Avel_old=ABS(Avel_old-LATT_CUR%Avel)
            vel_old=ABS(vel_old-DYN%VEL)
            err1=MAXVAL(Avel_old)
            err2=MAXVAL(vel_old)
            IF ((err1<tol) .AND. (err2<tol)) EXIT
          ENDDO

          !update acceleration for lattice components for the time slice (t)
          Aaccel=Aaccel+SIF_momentum+kinetic_stress

          CALL APPLY_STRESS_CONSTRAINT(LATT_CUR%AVEL, Aaccel, LATTICE_CONSTRAINTS)

          !c compute temperature for the slice t
          CALL KINETIC_E(EKIN_at,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,&
                       & DYN%POTIM,LATT_CUR%A,DYN%VEL)

          TEIN_at = 2*EKIN_at/BOLKEV/NDEGREES_OF_FREEDOM_

          CALL KINETIC_E_Lat(EKIN_Lat,AMASS,DYN%POTIM,LATT_CUR%AVEL)
   
          !!c upper triangle only
          ECONST=0._q
          !!TEIN_lat = 2*EKIN_Lat/BOLKEV/6
          !!TEIN=2*(EKIN+EKIN_Lat) /BOLKEV/(NDEGREES_OF_FREEDOM_+6) 
          TEIN_lat = 2*EKIN_Lat/BOLKEV/LATTICE_DOF
          EKIN=EKIN_at+EKIN_lat
          TEIN=2*EKIN/BOLKEV/(NDEGREES_OF_FREEDOM_+LATTICE_DOF) 
          

          !c v(t+dt/2)=v(t)dt+f(t)/(2m)*dt^2
          LATT_CUR%Avel=LATT_CUR%Avel+Aaccel 

          ! set Avel here
          CALL APPLY_STRESS_CONSTRAINT(LATT_CUR%AVEL, Aaccel, LATTICE_CONSTRAINTS)

          DYN%VEL=DYN%VEL+DYN%D2C+D2C_momentum+D2C_stoch-MATMUL(Gdot,DYN%VEL)

          !c add contribution from bias potential
!           IF (iconst5 .GT. 0) THEN
!             LATT_CUR%Avel=LATT_CUR%Avel-(hills_accel_L+penalty_accel_L)
!             DYN%VEL=DYN%VEL-(hills_accel+penalty_accel)
!           ENDIF

          IF ((iconst5 .GT. 0) .OR. (iconst8 .GT. 0) .OR. (iconst4 .GT. 0) ) THEN
            LATT_CUR%Avel=LATT_CUR%Avel-(hills_accel_L+penalty_accel_L+hspring_accel_L)
            DYN%VEL=DYN%VEL-(hills_accel+penalty_accel+hspring_accel)
          ENDIF

          !c r(t+dt)=r(t)+v'(t+dt)*dt
          AC=LATT_CUR%A+LATT_CUR%Avel
          DYN%POSION=DYN%POSIOC+DYN%VEL

          !c geometric constraints (RATTLE)
          ECONST=0._q
          IF (ICOORDINATES%NUMINTERNALS>0) THEN
            CALL Rattle(counter,T_INFO,INFO,DYN,LATT_CUR,IO,C_FORCE,C_FORCE_L,V_FORCE,V_FORCE_L,ICOORDINATES,g_io,3*T_INFO%NIONS+9,&
            &   AMASS,AC,ECONST,iconst0,iconst2,0)
            DYN%POSION=DYN%POSION+C_FORCE
            AC=AC+C_FORCE_L
            DYN%VEL=DYN%VEL+C_FORCE
            LATT_CUR%Avel=LATT_CUR%Avel+C_FORCE_L
          ENDIF

          IF (IO%IU6>0) THEN     
            write(g_io%REPORT,FMT='(/,3X,A4,4(X,F10.3))') 't_b>', DYN%TEMP,TEIN_at, TEIN_lat,TEIN
            CALL ENERGY_OUT(g_io,DYN,T_INFO,EKIN,EPOT,ECONST,0._q,0._q,0._q,counter)

            !c output for metadynamics
            IF (iconst5>0) THEN
              write(g_io%REPORT,ADVANCE='YES',FMT='(/,2X,A13)') '>Metadynamics'
              DO i=1,ICOORDINATES%NUMINTERNALS
                IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==5) THEN
                  write(g_io%REPORT,ADVANCE='YES',FMT='(3X,A6,X,F9.5)') &
                  & 'fic_p>',ICOORDINATES%COORDSTRUCT(i)%VALUE
                ENDIF
              ENDDO   
            ENDIF 
          ENDIF

          !c if this MD subroutine is used for relaxation
          IF (DYN%TEBEG<0.1) THEN

            !c take care of constraints
            IF (iconst0>0) THEN
              !c compute the total restoring force
              C_FORCE=C_FORCE+V_FORCE
              V_FORCE=0._q

              V_FORCE=DYN%D2C+C_FORCE

              FACT=(DYN%POTIM**2)*EVTOJ/AMTOKG *1E-10_q
              NI=1
              DO NT=1,T_INFO%NTYP
                DO NI=NI,T_INFO%NITYP(NT)+NI-1
                  V_FORCE(1,NI)=V_FORCE(1,NI)/FACT*2*T_INFO%POMASS(NT)
                  V_FORCE(2,NI)=V_FORCE(2,NI)/FACT*2*T_INFO%POMASS(NT)
                  V_FORCE(3,NI)=V_FORCE(3,NI)/FACT*2*T_INFO%POMASS(NT)
                ENDDO
              ENDDO

              CALL DIRKAR(T_INFO%NIONS,V_FORCE,LATT_CUR%A)

             !c write ionic forces after removing the restoring forces
              IF (IO%IU0>0) THEN
                WRITE(IO%IU6,72)
                DO J=1,T_INFO%NIONS
                  VTMP=T_INFO%POSION(1:3,J)
                  CALL  DIRKAR(1,VTMP,LATT_CUR%A)
                  WRITE(IO%IU6,76) VTMP,(V_FORCE(I,J),I=1,3)
                ENDDO
              ENDIF 
            ENDIF

            !c test the convergence criteria
            maxForce=0._q
            DO i=1,T_INFO%NIONS
              dummy=V_FORCE(1,i)*V_FORCE(1,i)+V_FORCE(2,i)*V_FORCE(2,i)+V_FORCE(3,i)*V_FORCE(3,i)
              IF (dummy>maxForce) &
                & maxForce=dummy
            END DO
            maxForce=maxForce**0.5
            IF (IO%IU0>0) THEN
              write(*,*) 'maxForce',maxForce,DYN%EDIFFG
            ENDIF 

            IF (DYN%EDIFFG<0.) THEN
              LSTOP2=.TRUE.

              IF (maxForce>ABS(DYN%EDIFFG)) LSTOP2=.FALSE.

              IF (LSTOP2) THEN
                INFO%LSTOP=.TRUE.
                !DYN%POSIOC=DYN%POSION
                DYN%POSION=DYN%POSIOC
                RETURN
              END IF
            
            ELSE IF (counter>1) THEN
            !IF (counter>1) THEN
              IF (IO%IU0>0) THEN
                write(*,*) 'dEPOT',ABS(EPOT-EPOT_old),DYN%EDIFFG
              ENDIF
              INFO%LSTOP=(ABS(EPOT-EPOT_old)<DYN%EDIFFG)
            ENDIF
          ENDIF


          !c update direct and reciprocal lattice:
          LATT_CUR%A=AC
          CALL LATTIC(LATT_CUR)

         !c first rotate cell with cell vectors to triangular matrix
         CALL ROTATE_CELL(IO%IU6,LATT_CUR,ROTMAT2)

         !c transform the triangular matrix with lattice vectors back to the original form 
         CALL ROTATE_BACK(LATT_CUR,ROTMAT1) 

         ROTMAT1=ROTMAT2
!           !c transform the triangular matrix with lattice vectors back to the original form 
!           CALL ROTATE_BACK(LATT_CUR,ROTMAT)
 
          DISMAX=0
          DO I=1,T_INFO%NIONS
             DISMAX=MAX(DISMAX, &
             &(DYN%VEL(1,I))**2+ &
             &(DYN%VEL(2,I))**2+ &
             &(DYN%VEL(3,I))**2)
          ENDDO
          DISMAX=SQRT(DISMAX)
          EPS=0.0_q
          ES =0.0_q
        END SUBROUTINE STEP_LANGEVIN_ISIF3

        SUBROUTINE isotrop( m)
          REAL(q) :: m(3,3)
          REAL(q) :: mean

          mean=(m(1,1)+m(2,2)+m(3,3))/3
          m=0
          m(1,1)=mean
          m(2,2)=mean
          m(3,3)=mean
        END SUBROUTINE isotrop

        SUBROUTINE APPLY_STRESS_CONSTRAINT(AVEL, STRESS, FREE_AXES)
          REAL(q), INTENT(INOUT) :: AVEL(3,3), STRESS(3,3)
          LOGICAL, INTENT(IN) :: FREE_AXES(3)
          INTEGER :: I, J
!          RETURN
          DO I = 1,3
            IF (.NOT. FREE_AXES(I)) THEN
              AVEL(I,:) = 0.0_q
              AVEL(:,I) = 0.0_q
!              STRESS(I,:) = 0.0_q
!              STRESS(:,I) = 0.0_q
            END IF
          END DO
        END SUBROUTINE APPLY_STRESS_CONSTRAINT

        SUBROUTINE ROTATE_CELL(IU,LATT_CUR,ROTMAT)
          TYPE(latt) :: LATT_CUR
          REAL(q) :: newA(3,3)
          REAL(q) :: oldA(3,3)
          REAL(q) :: atmp(3,3) 
          REAL(q) :: vtmp(3),vtmpNorm
          REAL(q) :: ROTMAT(3,3),ROTMAT_inv(3,3)
          REAL(q) :: r1,r2,r3
          REAL(q) :: cos23,cosX,cosY,sin23,sinX,sinY
          INTEGER :: ninfo,IU


          ROTMAT=0._q

          oldA=LATT_CUR%A
          atmp=MATMUL(transpose(LATT_CUR%A),LATT_CUR%A) 

          newA=0._q
          newA(3,3) = sqrt( atmp(3,3) )
          newA(3,2) = atmp(3,2) / sqrt( atmp(3,3) )
          newA(2,2) = sqrt( atmp(2,2) -   &
          &   atmp(3,2) * atmp(3,2) / atmp(3,3) )
          newA(3,1) = atmp(3,1) / sqrt( atmp(3,3) )
          newA(2,1) = ( atmp(3,3) * atmp(2,1) - atmp(3,2) * atmp(3,1) ) / &
          &      sqrt( atmp(3,3) * atmp(3,3) * atmp(2,2) - &
          &      atmp(3,3) * atmp(3,2) * atmp(3,2) )
          newA(1,1) = sqrt( atmp(1,1) -  &
          &    newA(3,1) * newA(3,1) - newA(2,1) * newA(2,1) )


!           r1=SQRT(atmp(1,1))
!           r2=SQRT(atmp(2,2))
!           r3=SQRT(atmp(3,3))
!           cos23=atmp(2,3)/r2/r3
!           sin23=SQRT(1._q-cos23**2)
! 
!           newA=0._q
!           newA(3,3)=r3
!           newA(3,2)=r2*cos23
!           newA(2,2)=r2*sin23
! 
!           !c vector orthogonal to the (b,c) plane
!           vtmp=CROSSPROD(3,oldA(1:3,2),oldA(1:3,3))
!           vtmpNorm=SQRT(SUM(vtmp**2))
!           !vtmp=vtmp/vtmpNorm
!           cosX=SUM(oldA(1:3,1)*vtmp)/r1/vtmpNorm
!           !sinX=SQRT(1._q-cosX**2)
!           newA(1,1)=r1*ABS(cosX)
! 
!           !newA(1,1)=r1*sinX
! 
!           !c the part of vector a that is orthogonal to the (b,c) plane
!           vtmp=oldA(1:3,1)-r1*cosX*vtmp/vtmpNorm
!           !vtmp=oldA(1:3,1)-SUM(oldA(1:3,1)*vtmp)*vtmp
!           vtmpNorm=SQRT(SUM(vtmp**2))
!           cosY=SUM(vtmp*oldA(1:3,3))/r3/vtmpNorm
!           sinY=SQRT(1-cosY**2)
!           newA(3,1)=vtmpNorm*cosY
!           newA(2,1)=vtmpNorm*sinY

         

!           IF (IU>=0) THEN
!             write(*,*) "oldA"
!             write(*,*) oldA(1,1),oldA(1,2),oldA(1,3)
!             write(*,*) oldA(2,1),oldA(2,2),oldA(2,3)
!             write(*,*) oldA(3,1),oldA(3,2),oldA(3,3)
! 
!             write(*,*) "newA"
!             write(*,*) newA(1,1),newA(1,2),newA(1,3)
!             write(*,*) newA(2,1),newA(2,2),newA(2,3)
!             write(*,*) newA(3,1),newA(3,2),newA(3,3)
!           ENDIF

          !c compute the matrix that transforms the triangular matrix
          !c to the original matrix
          ROTMAT=newA
          CALL SVDINVERSE(ROTMAT,3,ninfo)
          ROTMAT=MATMUL(oldA,ROTMAT)
          ROTMAT_inv=ROTMAT
          CALL SVDINVERSE(ROTMAT_inv,3,ninfo)

          atmp=MATMUL(ROTMAT,newA)   
!           IF (IU>=0) THEN
!             write(*,*) "atmp"
!             write(*,*) atmp(1,1),atmp(1,2),atmp(1,3)
!             write(*,*) atmp(2,1),atmp(2,2),atmp(2,3)
!             write(*,*) atmp(3,1),atmp(3,2),atmp(3,3)
!           ENDIF

          !c update LATT_CUR
          LATT_CUR%A=newA
          CALL LATTIC(LATT_CUR)
          LATT_CUR%AVEL=MATMUL(ROTMAT_inv,LATT_CUR%AVEL)

        END SUBROUTINE ROTATE_CELL

        SUBROUTINE ROTATE_BACK(LATT_CUR,ROTMAT)
        !c rotate the triangular to the original matrix
          TYPE(latt) :: LATT_CUR
          REAL(q) :: newA(3,3)
          REAL(q) :: oldA(3,3)
          REAL(q) :: ROTMAT(3,3)

          oldA=LATT_CUR%A
          newA=MATMUL(ROTMAT,oldA)          
          LATT_CUR%A=newA
          CALL LATTIC(LATT_CUR)
          LATT_CUR%AVEL=MATMUL(ROTMAT,LATT_CUR%AVEL)
        END SUBROUTINE ROTATE_BACK

        SUBROUTINE COMPUTE_KINETIC_STRESS(T_INFO,LATT_CUR,DYN,kinetic_stress,TSIF,IU)
          TYPE(type_info) :: T_INFO
          TYPE(latt) :: LATT_CUR
          TYPE(dynamics) :: DYN 
          REAL(q) :: Vtmp(3,T_INFO%NIONS)
          REAL(q) :: tmp(3,3),kinetic_stress(3,3),TSIF(3,3)
          INTEGER ::NI,NT,IU
          REAL(q) :: fac,fakt

          kinetic_stress=0._q

          Vtmp=DYN%VEL/DYN%POTIM
          Vtmp=MATMUL(LATT_CUR%A,Vtmp)
          
          NI=1
          DO NT=1, T_INFO%NTYP
            DO NI=NI,T_INFO%NITYP(NT)+NI-1
              tmp=OUTERPROD(3,Vtmp(:,NI),Vtmp(:,NI))
              kinetic_stress=kinetic_stress+tmp*T_INFO%POMASS(nt)
            ENDDO
          ENDDO

          kinetic_stress=(kinetic_stress+transpose(kinetic_stress))/2
          kinetic_stress=kinetic_stress/LATT_CUR%OMEGA

          fakt = evtoj*1e22_q/latt_cur%omega
          fac  = amtokg *  &
          &     1e5_q  *  &
          &     1e5_q  *  &
          &     1e30_q *  &
          &     1e-8_q

          IF (IU>=0) THEN
          write(iu,'(''  ideal gas correction = '',f9.2,'' kB'')') &
     &     fac*(kinetic_stress(1,1)+kinetic_stress(2,2)+kinetic_stress(3,3))/3
      write(iu,'(''  Total+kin._'',6f12.3)') &
     &     tsif(1,1)*fakt + fac*kinetic_stress(1,1),  &
     &     tsif(2,2)*fakt + fac*kinetic_stress(2,2),  &
     &     tsif(3,3)*fakt + fac*kinetic_stress(3,3),  &
     &     tsif(1,2)*fakt + fac*kinetic_stress(1,2),  &
     &     tsif(2,3)*fakt + fac*kinetic_stress(2,3),  &
     &     tsif(3,1)*fakt + fac*kinetic_stress(3,1)
          ENDIF
       
          kinetic_stress=kinetic_stress*fac/fakt

        END SUBROUTINE COMPUTE_KINETIC_STRESS


        SUBROUTINE STEP_NANDERSEN(DYN,ICOORDINATES,T_INFO,INFO,LATT_CUR, &
          &      EKIN,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,NDEGREES_OF_FREEDOM_,IO,EPOT,&
          &      g_io,NSUBSYS,TSUBSYS,PSUBSYS,iconst0,iconst2,TEIN)
        !c Andersen thermostat
          TYPE(coordstructure) :: ICOORDINATES     !c working and reference coords.
          TYPE(dynamics) :: DYN
          TYPE(type_info) :: T_INFO
          TYPE (info_struct) INFO
          TYPE(latt) :: LATT_CUR
          TYPE (in_struct) :: IO
          TYPE(gadget_io) :: g_io
          REAL(q) :: EKIN,EPS,ES,DISMAX,TEMPER,EPOT,ECONST
          INTEGER :: i,j,ii,jj,NDEGREES_OF_FREEDOM,NDEGREES_OF_FREEDOM_
          INTEGER,SAVE :: counter=0
          INTEGER :: iconst0,iconst2
          REAL(q) :: CMASS(3,1)
          REAL(q) :: VEL_tmp1(3,T_INFO%NIONS),VEL_tmp2(3,T_INFO%NIONS),VEL_tmp3(3,T_INFO%NIONS)
          REAL(q) :: VEL_tmp(3,T_INFO%NIONS)
          REAL(q) :: C_FORCE(3,T_INFO%NIONS)
          REAL(q) :: dummy
          INTEGER :: random_counter(3)
          INTEGER :: NSUBSYS(3)
          REAL(q) :: TSUBSYS(3),PSUBSYS(3)
          REAL(q) :: EKIN1,EKIN2,EKIN3,TEIN1,TEIN2,TEIN3,TEIN
          INTEGER,SAVE :: NDOF(3)
          REAL(q) :: AMASS,C_FORCE_L(3,3)
          LOGICAL :: LSCALE,LCMASS

          LSCALE=.FALSE.
          LCMASS=.TRUE.

          !c add this step to total count
          counter=counter+1
          !!!DYN%POSIOC=DYN%POSION

          AMASS=0._q

          !c generate Anderesen particles with desired temperatures (3 subsystems assumed)
          VEL_tmp1=0._q;VEL_tmp2=0._q;VEL_tmp3=0._q
          DYN%TEMP=TSUBSYS(1)
          CALL init_velocities(VEL_tmp1,T_INFO,DYN,LATT_CUR,NDEGREES_OF_FREEDOM_,LSCALE,LCMASS)
          IF (NSUBSYS(2)>NSUBSYS(1)) THEN
            DYN%TEMP=TSUBSYS(2)
            CALL init_velocities(VEL_tmp2,T_INFO,DYN,LATT_CUR,NDEGREES_OF_FREEDOM_,LSCALE,LCMASS)
            IF (NSUBSYS(3)>NSUBSYS(2)) THEN
              DYN%TEMP=TSUBSYS(3)
              CALL init_velocities(VEL_tmp3,T_INFO,DYN,LATT_CUR,NDEGREES_OF_FREEDOM_,LSCALE,LCMASS)
            ENDIF
          ENDIF

          VEL_tmp=VEL_tmp1
          
          IF (NSUBSYS(1)<NSUBSYS(2)) THEN
            DO i=NSUBSYS(1)+1,NSUBSYS(2)
              VEL_tmp(1:3,i)=VEL_tmp2(1:3,i)
            ENDDO
            IF (NSUBSYS(2)<NSUBSYS(3)) THEN
              DO i=NSUBSYS(2)+1,NSUBSYS(3)
                VEL_tmp(1:3,i)=VEL_tmp3(1:3,i)
              ENDDO
            ENDIF
          ENDIF

          !c remove translational motion
          !CMASS=0._q
          !CALL GIVE_CMASS(T_INFO,VEL_tmp,CMASS)
          !DO i=1,T_INFO%NIONS
          !  DO j=1,3
          !    IF (T_INFO%LSFOR(j,i)) VEL_tmp(j,i)=VEL_tmp(j,i)-CMASS(j,1)
          !  ENDDO
          !ENDDO
          !CMASS=0._q        

          IF (counter==1) THEN
            !c velocity initialization
            IF (DYN%INIT==0) DYN%VEL=VEL_tmp

            !c determine number of degrees of freedom for each subsystem
            NDOF=0
            DO i=1,NSUBSYS(1)
              DO j=1,3
                IF (T_INFO%LSFOR(j,i)) NDOF(1)=NDOF(1)+1
              ENDDO
            ENDDO
            DO i=NSUBSYS(1)+1,NSUBSYS(2)
              DO j=1,3
                IF (T_INFO%LSFOR(j,i)) NDOF(2)=NDOF(2)+1
              ENDDO
            ENDDO
            DO i=NSUBSYS(2)+1,NSUBSYS(3)
              DO j=1,3
                IF (T_INFO%LSFOR(j,i)) NDOF(3)=NDOF(3)+1
              ENDDO
            ENDDO
          END IF

          DYN%VEL=DYN%VEL+2*DYN%D2C
          random_counter=0

          IF (DYN%INIT==1) THEN
            DO i=1,NSUBSYS(1)
              CALL RANE_ION(dummy)
              IF (dummy<PSUBSYS(1)) THEN
                random_counter(1)=random_counter(1)+1
                DYN%VEL(1:3,i)=VEL_tmp(1:3,i)
              ENDIF
            ENDDO
            DO i=NSUBSYS(1)+1,NSUBSYS(2)
              CALL RANE_ION(dummy)
              IF (dummy<PSUBSYS(2)) THEN
                random_counter(2)=random_counter(2)+1
                DYN%VEL(1:3,i)=VEL_tmp(1:3,i)
              ENDIF
            ENDDO
            DO i=NSUBSYS(2)+1,NSUBSYS(3)
              CALL RANE_ION(dummy)
              IF (dummy<PSUBSYS(3)) THEN
                random_counter(3)=random_counter(3)+1
                DYN%VEL(1:3,i)=VEL_tmp(1:3,i)
              ENDIF
            ENDDO
          ENDIF

          CALL GIVE_CMASS(T_INFO,DYN%VEL,CMASS)
          DO i=1,T_INFO%NIONS
            DO j=1,3
              IF (T_INFO%LSFOR(j,i)) DYN%VEL(j,i)=DYN%VEL(j,i)-CMASS(j,1)
            ENDDO
          ENDDO
          !CALL GIVE_CMASS(T_INFO,DYN%VEL,CMASS)

          C_FORCE=0._q
          DYN%POSION=DYN%POSIOC+DYN%VEL

          ECONST=0._q
          IF (ICOORDINATES%NUMINTERNALS>0) THEN
            CALL Lproj_out_step(T_INFO,INFO,DYN,LATT_CUR,IO,C_FORCE,C_FORCE_L,ICOORDINATES, &
            &   g_io,3*T_INFO%NIONS,AMASS,LATT_CUR%A,ECONST,iconst0,iconst2)
            DYN%POSION=DYN%POSION+C_FORCE
            DYN%VEL=DYN%VEL+C_FORCE
          ENDIF
          
          CALL KINETIC_E(EKIN,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,&
                       & DYN%POTIM,LATT_CUR%A,DYN%VEL)

          TEIN = 2*EKIN/BOLKEV/NDEGREES_OF_FREEDOM_

          IF (IO%IU6>0) CALL ENERGY_OUT(g_io,DYN,T_INFO,EKIN,EPOT,ECONST,0._q,0._q,0._q,counter)

          !c compute temperatures for all subsystems
          VEL_tmp=0._q
          VEL_tmp(1:3,1:NSUBSYS(1))=DYN%VEL(1:3,1:NSUBSYS(1))

          CALL KINETIC_E(EKIN1,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,&
                       & DYN%POTIM,LATT_CUR%A,VEL_tmp)
          TEIN1 = 2*EKIN1/BOLKEV/NDOF(1)*SUM(NDOF)/NDEGREES_OF_FREEDOM_ 


          IF (NSUBSYS(2)>NSUBSYS(1)) THEN
            VEL_tmp=0._q
            VEL_tmp(1:3,NSUBSYS(1)+1:NSUBSYS(2))=DYN%VEL(1:3,NSUBSYS(1)+1:NSUBSYS(2))
            CALL KINETIC_E(EKIN2,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,&
                        & DYN%POTIM,LATT_CUR%A,VEL_tmp)
            TEIN2 = 2*EKIN2/BOLKEV/NDOF(2)*SUM(NDOF)/NDEGREES_OF_FREEDOM_
            IF (NSUBSYS(3)>NSUBSYS(2)) THEN
              VEL_tmp=0._q
              VEL_tmp(1:3,NSUBSYS(2)+1:)=DYN%VEL(1:3,NSUBSYS(2)+1:)
              CALL KINETIC_E(EKIN3,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,&
                         & DYN%POTIM,LATT_CUR%A,VEL_tmp)                
              TEIN3 = 2*EKIN3/BOLKEV/NDOF(3)*SUM(NDOF)/NDEGREES_OF_FREEDOM_              
            END IF
          END IF

           
          IF (IO%IU6>0) THEN
            IF (NSUBSYS(2)>NSUBSYS(1)) THEN
              IF (NSUBSYS(3)>NSUBSYS(2)) THEN            
                write(g_io%REPORT,FMT='(/,2X,A32,X,I5,X,I5,X,I5)') '>Thermostat, num. of collisions:',random_counter(1),random_counter(2),random_counter(3)
                write(g_io%REPORT,FMT='(/,2X,A5,X,F10.3,X,F10.3,X,F10.3)') 'Tsub:',TEIN1,TEIN2,TEIN3
              ELSE
                write(g_io%REPORT,FMT='(/,2X,A32,X,I5,X,I5)') '>Thermostat, num. of collisions:',random_counter(1),random_counter(2)
                write(g_io%REPORT,FMT='(/,2X,A5,X,F10.3,X,F10.3)') 'Tsub:',TEIN1,TEIN2 
              ENDIF
            END IF 
          ENDIF

          DISMAX=0
          DO I=1,T_INFO%NIONS
             DISMAX=MAX(DISMAX, &
             &(DYN%VEL(1,I))**2+ &
             &(DYN%VEL(2,I))**2+ &
             &(DYN%VEL(3,I))**2)
          ENDDO
          DISMAX=SQRT(DISMAX)
          EPS=0.0_q
          ES =0.0_q
        END SUBROUTINE STEP_NANDERSEN


        SUBROUTINE HILLS_READER(IO,DYN,hills,g_io,ICOORDINATES,MDALGO)
        !c metadynamics - read input parameters
          USE reader_tags

          TYPE(coordstructure) :: ICOORDINATES
          TYPE (in_struct) :: IO
          TYPE(dynamics) :: DYN
          TYPE(hills_data) :: hills
          TYPE(gadget_io) :: g_io
          INTEGER IDUM,i,j,K,N,IERR
          INTEGER :: MDALGO
          LOGICAL :: LOPEN,LDUM
          REAL(q) :: RDUM
          COMPLEX(q) :: CDUM
          CHARACTER*1 :: CHARAC
          REAL(q) :: width,high,FACT,BMP,ti

          CALL OPEN_INCAR_IF_FOUND(IO%IU5, LOPEN)

          !c size of the bin to update the bias potential
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'HILLS_BIN', HILLS%BIN, IERR, WRITEXMLINCAR)
          IF (IERR==3) hills%bin=DYN%NSW
          IF (hills%bin<1) hills%bin=DYN%NSW
          IF (hills%bin>DYN%NSW) hills%bin=DYN%NSW

          ALLOCATE(hills%gauss(DYN%NSW/hills%bin))
          
           !c stride 
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'HILLS_MAXSTRIDE', HILLS%MAXSTRIDE, IERR, WRITEXMLINCAR)
          IF (IERR==3) hills%maxstride=1
          IF (hills%maxstride<1) hills%maxstride=1

          !c 
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'HILLS_STRIDE', HILLS%STRIDE, IERR, WRITEXMLINCAR)
          IF (IERR==3) hills%stride=1.5
          IF ( hills%stride .LE. 0._q) hills%stride=1.5

          !c hight of the gauss peak in the bias potential (eV)
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'HILLS_H', HIGH, IERR, WRITEXMLINCAR)
          IF (IERR==3) high=0.001
          IF ( high<0._q) high=0.001

          !c width of the gauss peak in the bias potential 
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'HILLS_W', WIDTH, IERR, WRITEXMLINCAR)
          IF (IERR==3) width=0.001
          IF ( width .LE. 0._q) width=0.001

          DO i=1,DYN%NSW/hills%bin
            ALLOCATE(hills%gauss(i)%position(hills%number))
            hills%gauss(i)%position(:)=0._q
            hills%gauss(i)%high=high
            hills%gauss(i)%width=width
          ENDDO

          !c is variable width desired?
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'HILLS_VARIABLE_W', HILLS%VARIABLE_WIDTH, IERR, WRITEXMLINCAR)
          IF (IERR==3) hills%variable_width=.FALSE.

          !c only for mMD with fictituous particles
          IF ((MDALGO==10 .OR. MDALGO==20) .AND. (hills%number>0)) THEN
            !c force constant to couple the fictious and real particles 
            K=SIZE(hills%force_constant)
            CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'HILLS_K', HILLS%FORCE_CONSTANT, K, IERR, WRITEXMLINCAR)
            IF (IERR==3) hills%force_constant=0.5
            DO i=1,K
              IF ( hills%force_constant(i)<0._q) hills%force_constant(i)=0.5
            ENDDO

            !c force constant to couple the fictious and real particles 
            K=SIZE(hills%mass)
            CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'HILLS_M', HILLS%MASS, K, IERR, WRITEXMLINCAR)
            IF (IERR==3) hills%mass=10.0
            DO i=1,K
              IF ( hills%mass(i) .LE. 0._q) hills%mass(i)=10.0
            ENDDO

            !c temperature of the fictitious particles
            CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'HILLS_TEMPERATURE', HILLS%TEMPERATURE, IERR, WRITEXMLINCAR)
            IF (IERR==3) hills%temperature=DYN%TEBEG
            IF ( hills%temperature<=0._q) hills%temperature=DYN%TEBEG

            !c metapositions
            CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'HILLS_POSITION', HILLS%POSITION, HILLS%NUMBER, IERR, WRITEXMLINCAR)
            IF (IERR==3) THEN
              j=0
              DO i =1,ICOORDINATES%NUMINTERNALS
                IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==5 .OR. ICOORDINATES%COORDSTRUCT(i)%STATUS==6) THEN
                  j=j+1
                  hills%position(j)=ICOORDINATES%COORDSTRUCT(i)%VALUE
                ENDIF
              ENDDO
            ENDIF

            !c metavelocities
            CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'HILLS_VELOCITY', HILLS%VELOCITY, HILLS%NUMBER, IERR, WRITEXMLINCAR)
            IF (IERR==3) THEN
              FACT= (AMTOKG/EVTOJ)*(1E-10_q/(DYN%POTIM*1E-15_q))**2
              DO i=1, hills%number
                BMP=SQRT(hills%temperature*BOLKEV/(hills%mass(i)*FACT))
                hills%velocity(i)=boltzmann_distribution(0._q,BMP)
              ENDDO
              ti=0._q
              DO i=1, hills%number
                ti=ti+(hills%velocity(i)/DYN%POTIM)**2*hills%mass(i)
              ENDDO
              ti=ti*AMTOKG/BOLKEV/EVTOJ*1E10_q/hills%number
              hills%velocity(:)= hills%velocity(:)*SQRT(hills%temperature/ti)
            ENDIF

            !c probability of collision for the andersen thermostat
            !c for the fictitious particles
            CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'HILLS_ANDERSEN_PROB', HILLS%ANDERSEN_PROB, IERR, WRITEXMLINCAR)
            IF (IERR==3) hills%andersen_prob=0._q
            IF ( hills%andersen_prob<0._q) hills%andersen_prob=0._q

            !c mass for nose-hover for fictitious particles
            CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'HILLS_SQQ', HILLS%SQQ, IERR, WRITEXMLINCAR)
            IF (IERR==3) hills%sqq=100._q
            IF ( hills%sqq<0._q) hills%sqq=100._q
          ENDIF

          !write (*,*) 'hills%mass',hills%mass
          IF (IO%IU6>0) THEN
            !c write the input parameters summary into OUTCAR
            write(IO%IU6,FMT='(3X,A22,X,I10)')   '          HILLS_BIN = ',  hills%bin
            write(IO%IU6,FMT='(3X,A22,X,F10.8)') '            HILLS_H = ',  high
            write(IO%IU6,FMT='(3X,A22,X,F10.8)') '            HILLS_W = ' , width
            write(IO%IU6,FMT='(3X,A22,X,I5)')    '    HILLS_MAXSTRIDE = ',  hills%maxstride            
            write(IO%IU6,FMT='(3X,A22,X,F10.8)') '       HILLS_STRIDE = ',  hills%stride
            write(IO%IU6,FMT='(3X,A22,X,L10)')   '   HILLS_VARIABLE_W = ' , hills%variable_width

            !c write the input parameters summary into REPORT 
            write(g_io%REPORT,FMT='(3X,A22,X,I10)')   '          HILLS_BIN = ',  hills%bin
            write(g_io%REPORT,FMT='(3X,A22,X,F10.8)') '            HILLS_H = ',  high
            write(g_io%REPORT,FMT='(3X,A22,X,F10.8)') '            HILLS_W = ' , width
            write(g_io%REPORT,FMT='(3X,A22,X,I5)')    '    HILLS_MAXSTRIDE = ',  hills%maxstride            
            write(g_io%REPORT,FMT='(3X,A22,X,F10.8)') '       HILLS_STRIDE = ',  hills%stride
            write(g_io%REPORT,FMT='(3X,A22,X,L10)')   '   HILLS_VARIABLE_W = ' , hills%variable_width

            !c some parameters for the mMD with fictitious variables
            IF ((MDALGO==10 .OR. MDALGO==20) .AND. (hills%number>0) ) THEN
              !c write the input parameters summary into OUTCAR
              write(IO%IU6,FMT='(3X,A22,X,F10.3)') '  HILLS_TEMPERATURE = ',  hills%temperature
              write(IO%IU6,FMT='(3X,A22,X,F10.8)') 'HILLS_ANDERSEN_PROB = ',  hills%andersen_prob

              write(IO%IU6,ADVANCE='NO',FMT='(3X,A22)') '            HILLS_K = '
              DO i=1,k
                write(IO%IU6,ADVANCE='NO',FMT='(X,F10.5)') hills%force_constant(i)
              ENDDO
              write(IO%IU6,ADVANCE='NO',FMT='(/)')

              write(IO%IU6,ADVANCE='NO',FMT='(3X,A22)') '            HILLS_M = '
              DO i=1,k
                write(IO%IU6,ADVANCE='NO',FMT='(X,F10.6)') hills%mass(i)
              ENDDO
              write(IO%IU6,ADVANCE='NO',FMT='(/)')

              write(IO%IU6,ADVANCE='NO',FMT='(3X,A22)') '     HILLS_POSITION = '
              DO i=1,hills%number
                write(IO%IU6,ADVANCE='NO',FMT='(X,F10.6)') hills%position(i)
              ENDDO
              write(IO%IU6,ADVANCE='NO',FMT='(/)')

              write(IO%IU6,ADVANCE='NO',FMT='(3X,A22)') '     HILLS_VELOCITY = '
              DO i=1,hills%number
                write(IO%IU6,ADVANCE='NO',FMT='(X,F10.6)') hills%velocity(i)
              ENDDO
              write(IO%IU6,ADVANCE='NO',FMT='(/)')

              write(IO%IU6,FMT='(3X,A22,X,F10.3)')   '          HILLS_SQQ = ',  hills%sqq

              !c write the input parameters summary into REPORT 
              write(g_io%REPORT,FMT='(3X,A22,X,F10.3)') '  HILLS_TEMPERATURE = ',  hills%temperature
              write(g_io%REPORT,FMT='(3X,A22,X,F10.8)') 'HILLS_ANDERSEN_PROB = ',  hills%andersen_prob

              write(g_io%REPORT,ADVANCE='NO',FMT='(3X,A22)') '            HILLS_K = '
              DO i=1,k
                write(g_io%REPORT,ADVANCE='NO',FMT='(X,F10.5)') hills%force_constant(i)
              ENDDO
              write(g_io%REPORT,ADVANCE='NO',FMT='(/)')

              write(g_io%REPORT,ADVANCE='NO',FMT='(3X,A22)') '            HILLS_M = '
              DO i=1,k
                write(g_io%REPORT,ADVANCE='NO',FMT='(X,F10.6)') hills%mass(i)
              ENDDO
              write(g_io%REPORT,ADVANCE='NO',FMT='(/)')

              write(g_io%REPORT,ADVANCE='NO',FMT='(3X,A22)') '     HILLS_POSITION = '
              DO i=1,hills%number
                write(g_io%REPORT,ADVANCE='NO',FMT='(X,F10.6)') hills%position(i)
              ENDDO
              write(g_io%REPORT,ADVANCE='NO',FMT='(/)')

              write(g_io%REPORT,ADVANCE='NO',FMT='(3X,A22)') '     HILLS_VELOCITY = '
              DO i=1,hills%number
                write(g_io%REPORT,ADVANCE='NO',FMT='(X,F10.6)') hills%velocity(i)
              ENDDO
              write(g_io%REPORT,ADVANCE='NO',FMT='(/)')

              write(g_io%REPORT,FMT='(3X,A22,X,F10.3)')   '          HILLS_SQQ = ',  hills%sqq
            ENDIF
          ENDIF
          
          CALL CLOSE_INCAR_IF_FOUND(IO%IU5)

        END SUBROUTINE HILLS_READER

        SUBROUTINE FBIAS_READER(IO,iconst4,g_io,fbias_R0,fbias_D,fbias_A)
        !c Fermi-function-shaped bias potential (status = 4) - read in the input parameters
          USE reader_tags

          TYPE (in_struct) :: IO
          TYPE(gadget_io) :: g_io
          INTEGER IDUM,i,j,K,N,IERR
          LOGICAL :: LOPEN,LDUM
          REAL(q) :: RDUM
          COMPLEX(q) :: CDUM
          CHARACTER*1 :: CHARAC
          INTEGER :: iconst4  
          REAL(q) :: fbias_R0(iconst4),fbias_D(iconst4),fbias_A(iconst4) 

          IF (iconst4 .LE. 0) RETURN
     
          fbias_R0=0._q;fbias_D=0._q;fbias_A=0._q

          CALL OPEN_INCAR_IF_FOUND(IO%IU5, LOPEN)

          !c position 
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'FBIAS_R0', fbias_R0, ICONST4, IERR, WRITEXMLINCAR)
          IF (IERR==3) THEN
            CALL vtutor%error("Error: coordinates for fbias defined but the parameter 'FBIAS_R0' not &
               &provided in INCAR")
          ENDIF
          !c sanity test
          IF (IERR==0) THEN             
            DO i=1,iconst4
              IF (fbias_R0(i)==0) THEN
                CALL vtutor%error("Error 'FBIAS_R0(i)=0',i= " // str(i))
              ENDIF
            ENDDO
          ENDIF

          !c smearing parameter 
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'FBIAS_D', fbias_D, ICONST4, IERR, WRITEXMLINCAR)
          IF (IERR==3) THEN
            CALL vtutor%error("Error: coordinates for fbias defined but the parameter 'FBIAS_D' not &
               &provided in INCAR")
          ENDIF

          !c scaling parameter
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'FBIAS_A', fbias_A, ICONST4, IERR, WRITEXMLINCAR)
          IF (IERR==3) THEN
            CALL vtutor%error("Error: coordinates for fbias defined but the parameter 'FBIAS_A' not &
               &provided in INCAR")
          ENDIF

          IF (IO%IU6>0) THEN

            write(g_io%REPORT,ADVANCE='NO',FMT='(3X,A22)') '           FBIAS_R0 = '
            DO i=1,iconst4
                write(g_io%REPORT,ADVANCE='NO',FMT='(X,F10.5)') FBIAS_R0(i)
            ENDDO
            write(g_io%REPORT,ADVANCE='NO',FMT='(/)')
            write(g_io%REPORT,ADVANCE='NO',FMT='(3X,A22)') '            FBIAS_D = '
            DO i=1,iconst4
                write(g_io%REPORT,ADVANCE='NO',FMT='(X,F10.5)') FBIAS_D(i)
            ENDDO
            write(g_io%REPORT,ADVANCE='NO',FMT='(/)')
            write(g_io%REPORT,ADVANCE='NO',FMT='(3X,A22)') '            FBIAS_A = '
            DO i=1,iconst4
                write(g_io%REPORT,ADVANCE='NO',FMT='(X,F10.5)') FBIAS_A(i)
            ENDDO 
            write(g_io%REPORT,ADVANCE='NO',FMT='(/)')


            !c write the input parameters summary into OUTCAR
            write(IO%IU6,ADVANCE='NO',FMT='(3X,A22)') '           FBIAS_R0 = '
            DO i=1,iconst4
              write(IO%IU6,ADVANCE='NO',FMT='(X,F10.5)') FBIAS_R0(i)
            ENDDO
            write(IO%IU6,ADVANCE='NO',FMT='(/)')

            write(IO%IU6,ADVANCE='NO',FMT='(3X,A22)') '           FBIAS_D = '
            DO i=1,iconst4
              write(IO%IU6,ADVANCE='NO',FMT='(X,F10.5)') FBIAS_D(i)
            ENDDO
            write(IO%IU6,ADVANCE='NO',FMT='(/)')

            write(IO%IU6,ADVANCE='NO',FMT='(3X,A22)') '           FBIAS_A = '
            DO i=1,iconst4
              write(IO%IU6,ADVANCE='NO',FMT='(X,F10.5)') FBIAS_A(i)
            ENDDO
            write(IO%IU6,ADVANCE='NO',FMT='(/)')
          ENDIF

          CALL CLOSE_INCAR_IF_FOUND(IO%IU5)
        END SUBROUTINE FBIAS_READER
        
        SUBROUTINE HARMONIC_SPRING_READER(IO,iconst8,g_io,spring_R0,spring_V0,spring_K,Lspring_V0)
        !c harmonic springs acting on internal coordinates (status = 8) - read in the input parameters
          USE reader_tags

          TYPE (in_struct) :: IO
          TYPE(gadget_io) :: g_io
          INTEGER IDUM,i,j,K,N,IERR
          LOGICAL :: LOPEN,LDUM
          REAL(q) :: RDUM
          COMPLEX(q) :: CDUM
          CHARACTER*1 :: CHARAC
          INTEGER :: iconst8
          REAL(q) :: spring_R0(iconst8),spring_V0(iconst8),spring_K(iconst8)
          LOGICAL :: Lspring_V0
           
          IF (iconst8 .LE. 0) RETURN
        
          SPRING_R0=0._q
          SPRING_K=0._q
          spring_V0=0._q
        
          CALL OPEN_INCAR_IF_FOUND(IO%IU5, LOPEN)
         
          !c force constant (in a.u.) 
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'SPRING_K', SPRING_K, ICONST8, IERR, WRITEXMLINCAR)
          IF (IERR==3) THEN
            CALL vtutor%error("Error: spring coordinates defined but the parameter 'SPRING_K' not &
               &provided in INCAR")
          ENDIF

          !c reference value (in A, rad., etc...)  
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'SPRING_R0', SPRING_R0, ICONST8, IERR, WRITEXMLINCAR)
          IF (IERR==3) THEN
            CALL vtutor%error("Error: spring coordinates defined but the parameter 'SPRING_R0' not &
               &provided in INCAR")
          ENDIF

          !c velocity for R0 (in A/fs, rad./fs, etc...)  
          !c optional, only for harmonic springs
          SPRING_V0=0._q
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'SPRING_V0', SPRING_V0, ICONST8, IERR, WRITEXMLINCAR)
          IF (IERR==3) THEN
            WRITE(IO%IU0,*)'Warning: spring coordinates defined but the parameter ''SPRING_V0'' not provided in INCAR'
          ELSE
            Lspring_V0=.TRUE.
          ENDIF

          IF (IO%IU6>0) THEN
                    
            write(g_io%REPORT,ADVANCE='NO',FMT='(3X,A22)') '           SPRING_K = '
            DO i=1,iconst8
                write(g_io%REPORT,ADVANCE='NO',FMT='(X,F10.5)') SPRING_K(i)
            ENDDO
            write(g_io%REPORT,ADVANCE='NO',FMT='(/)')
            write(g_io%REPORT,ADVANCE='NO',FMT='(3X,A22)') '          SPRING_R0 = '
            DO i=1,iconst8
                write(g_io%REPORT,ADVANCE='NO',FMT='(X,F10.5)') SPRING_R0(i)
            ENDDO
            write(g_io%REPORT,ADVANCE='NO',FMT='(/)')
            
            IF (Lspring_V0) THEN
              write(g_io%REPORT,ADVANCE='NO',FMT='(3X,A22)') '          SPRING_V0 = '
              DO i=1,iconst8
                write(g_io%REPORT,ADVANCE='NO',FMT='(X,F10.5)') SPRING_V0(i)
              ENDDO
              write(g_io%REPORT,ADVANCE='NO',FMT='(/)')
            ENDIF 

            !c write the input parameters summary into OUTCAR
            write(IO%IU6,ADVANCE='NO',FMT='(3X,A22)') '           SPRING_K = '
            DO i=1,iconst8
              write(IO%IU6,ADVANCE='NO',FMT='(X,F10.5)') SPRING_K(i)
            ENDDO
            write(IO%IU6,ADVANCE='NO',FMT='(/)')
            
            write(IO%IU6,ADVANCE='NO',FMT='(3X,A22)') '          SPRING_R0 = '
            DO i=1,iconst8
              write(IO%IU6,ADVANCE='NO',FMT='(X,F10.5)') SPRING_R0(i)
            ENDDO
            write(IO%IU6,ADVANCE='NO',FMT='(/)')

            IF (Lspring_V0) THEN
              write(IO%IU6,ADVANCE='NO',FMT='(3X,A22)') '          SPRING_V0 = '
              DO i=1,iconst8
                write(IO%IU6,ADVANCE='NO',FMT='(X,F10.5)') SPRING_V0(i)
              ENDDO
              write(IO%IU6,ADVANCE='NO',FMT='(/)')
            ENDIF
          ENDIF

          CALL CLOSE_INCAR_IF_FOUND(IO%IU5)
          
        END SUBROUTINE HARMONIC_SPRING_READER

        SUBROUTINE PENALTY_READER(hills,penalty,IO,g_io)
        !c metadynamics - read existing bias potential and rigid walls 
          USE reader_tags

          TYPE(penalty_data) :: penalty
          TYPE(hills_data) :: hills
          TYPE(gadget_io) :: g_io
          TYPE(gauss_peak),POINTER :: reallocate_(:)
          TYPE (in_struct) :: IO
          INTEGER :: counter,ios,i,mold
          INTEGER IDUM,K,N,IERR
          LOGICAL :: LOPEN,LDUM
          REAL(q) :: RDUM
          COMPLEX(q) :: CDUM
          CHARACTER*1 :: CHARAC
          CHARACTER(LEN=2000) :: line
          
          ALLOCATE(penalty%wall(hills%number,2))
          penalty%wall(:,1)=-1000.0
          penalty%wall(:,2)=1000.0

          CALL OPEN_INCAR_IF_FOUND(IO%IU5, LOPEN)

          !c 
          K=hills%number
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'HILLS_WALL_LOWER', PENALTY%WALL(:,1), K, IERR, WRITEXMLINCAR)

          IF (IERR==3) penalty%wall(:,1)=-1000.0
         
          CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'HILLS_WALL_UPPER', PENALTY%WALL(:,2), K, IERR, WRITEXMLINCAR)
          IF (IERR==3) penalty%wall(:,2)=1000.0

          CALL CLOSE_INCAR_IF_FOUND(IO%IU5)

          ALLOCATE(penalty%force(hills%number))
          penalty%force=0._q

          !OPEN(UNIT=g_io%PENALTY,FILE=DIR_APP(1:DIR_LEN)//'PENALTYPOT',STATUS='UNKNOWN')
          OPEN(UNIT=g_io%PENALTY,FILE='PENALTYPOT',STATUS='UNKNOWN')
          
          penalty%number=10
          ALLOCATE(penalty%gauss(penalty%number))
          DO i=1,penalty%number
            ALLOCATE(penalty%gauss(i)%position(hills%number))
            penalty%gauss(i)%position(:)=0._q
          ENDDO

          counter=1
          DO
            IF (counter .GT. penalty%number) THEN
              mold=penalty%number
              penalty%number=penalty%number+10
              ALLOCATE(reallocate_(penalty%number))
              IF (.NOT. ASSOCIATED(penalty%gauss)) RETURN
              reallocate_(1:mold)=penalty%gauss(1:mold)
              DEALLOCATE(penalty%gauss)
              penalty%gauss=>reallocate_
              DO i=mold+1,penalty%number
                ALLOCATE(penalty%gauss(i)%position(hills%number))
                penalty%gauss(i)%position(:)=0._q
              ENDDO
              penalty%gauss(mold+1:)%high=0._q;penalty%gauss(mold+1:)%width=0._q 
            ENDIF

            READ(g_io%PENALTY,FMT='(A2000)',IOSTAT=ios) line
            IF (ios/=0) THEN 
              counter=counter-1
              EXIT
            ENDIF
            IF (LEN_TRIM(line)>0) THEN
              line=TRIM(ADJUSTL(line))
              READ(line,FMT=*,IOSTAT=ios) penalty%gauss(counter)%position(:),&
              & penalty%gauss(counter)%high,penalty%gauss(counter)%width
              IF (ios/=0) THEN 
                CALL vtutor%error("Error reading PENALTYPOT (item " // str(counter) // ")")
              ELSE
                IF (penalty%gauss(counter)%width <= 0._q) THEN
                  CALL vtutor%error("Error reading PENALTYPOT (item " // str(counter) // "): W must &
                     &be a positive number")
                ENDIF
                counter=counter+1
              ENDIF
            ENDIF
          ENDDO

          penalty%number=counter
          ALLOCATE(reallocate_(counter))
          IF (.NOT. ASSOCIATED(penalty%gauss)) RETURN
          reallocate_(1:counter)=penalty%gauss(1:counter)
          DEALLOCATE(penalty%gauss)
          penalty%gauss=>reallocate_
              
          CLOSE(g_io%PENALTY)
        END SUBROUTINE PENALTY_READER
 

        SUBROUTINE HILLS_METHOD(DYN,ICOORDINATES,hills,penalty,T_INFO,INFO,LATT_CUR, &
          &      EKIN,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,IO,EPOT,ANDERSEN_PROB,g_io,& 
          &      WDES,iconst0,iconst2,TEIN)
        !c NVT metadynamics (Andersen thermostat) with fictitious coordinates 
          TYPE(coordstructure) :: ICOORDINATES     !c working and reference coords.
          TYPE(dynamics) :: DYN
          TYPE(type_info) :: T_INFO
          TYPE (info_struct) INFO
          TYPE(latt) :: LATT_CUR
          TYPE (in_struct) :: IO
          TYPE(gadget_io) :: g_io
          TYPE(hills_data) :: hills
          TYPE (wavedes)  ::   WDES
          TYPE(penalty_data) :: penalty
          REAL(q) :: C_FORCE(3,T_INFO%NIONS)
          REAL(q) :: EKIN,EPS,ES,DISMAX,TEMPER,EPOT,ti,ti_old,T_old,T_new,ECONST,TEIN
          INTEGER :: i,j,ii,jj,NDEGREES_OF_FREEDOM,random_counter
          INTEGER,SAVE :: counter=0        
          INTEGER :: iconst0,iconst2
          REAL(q) :: dummy,FACT,BMP,ANDERSEN_PROB
          REAL(q) :: CMASS(3,1)
          REAL(q) :: VEL_tmp(3,T_INFO%NIONS)
          REAL(q) :: spring_accel(3,T_INFO%NIONS)
          REAL(q),SAVE,ALLOCATABLE :: spring_forceS(:)
          REAL(q) :: AMASS,C_FORCE_L(3,3)
          !LOGICAL :: LSCALE,LCMASS
          REAL(q) :: TOTAL_MASS
          INTEGER,SAVE :: Ltxyz(3)

          !LSCALE=.FALSE.
          !LCMASS=.TRUE.


          !c add this step to total count
          counter=counter+1

          !c determine if coordinates for some direction
          !c are frozen
          IF (counter==1) THEN
            Ltxyz(:)=1
            DO j=1,3
              DO i=1,T_INFO%NIONS
                IF (.NOT. T_INFO%LSFOR(j,i)) Ltxyz(j)=0
              ENDDO
            ENDDO

            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0) THEN
                IF (ICOORDINATES%COORDSTRUCT(i)%TAG=="X ") Ltxyz(1)=0
                IF (ICOORDINATES%COORDSTRUCT(i)%TAG=="Y ") Ltxyz(2)=0
                IF (ICOORDINATES%COORDSTRUCT(i)%TAG=="Z ") Ltxyz(3)=0
              ENDIF
            ENDDO
          ENDIF
         
          AMASS=0._q
          ti_old=0.

          IF (counter==1) THEN
            ALLOCATE(spring_forceS(hills%number))            
            spring_forceS=0._q
            j=0
            DO i=1,ICOORDINATES%NUMINTERNALS
              !IF (penalty%number==0) THEN
                IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==5 .OR. ICOORDINATES%COORDSTRUCT(i)%STATUS==6) THEN
                  j=j+1
                 !!! hills%position(j)=ICOORDINATES%COORDSTRUCT(i)%VALUE
                  IF (hills%position(j)<=(penalty%wall(j,1)+hills%gauss(1)%width)) THEN
                    hills%position(j)=hills%position(j)+2*hills%gauss(1)%width
                  ENDIF
                  IF (hills%position(j)>=(penalty%wall(j,2)-hills%gauss(1)%width)) THEN
                    hills%position(j)=hills%position(j)-2*hills%gauss(1)%width
                  ENDIF
                ENDIF
            ENDDO
          ENDIF 

          !c calculate additional potentaials (spring, bias, and penalty)
          CALL hills_spring_force(T_INFO,LATT_CUR,DYN,hills,ICOORDINATES, DYN%POSIOC,spring_forceS,spring_accel)

          IF (counter .GE. hills%bin) THEN
            CALL hills_bias_potential(hills,penalty,ICOORDINATES,DYN,counter,g_io,IO)
          ENDIF 

          CALL penalty_potential(penalty,hills,ICOORDINATES)

          !c real variables +  andersen thermostat
          VEL_tmp=0._q
          !CALL init_velocities(VEL_tmp,T_INFO,DYN,LATT_CUR,NDEGREES_OF_FREEDOM,LSCALE,LCMASS)
            
          ! Add a random velicity that is exactly
          ! distributed like the centre of mass velocity
          TOTAL_MASS = 0.0_q
          DO I = 1, T_INFO%NIONS
            TOTAL_MASS = TOTAL_MASS + T_INFO%POMASS(T_INFO%ITYP(I))
          END DO
          CMASS(:,1) = BOLTZMANN_VELOCITY(DYN, LATT_CUR, TOTAL_MASS)

          !c don't add translation for frozen direction
          DO i=1,3
            IF (Ltxyz(i)==0) CMASS(i,1)=0._q
          ENDDO

          DO i=1,T_INFO%NIONS
            DO j=1,3
              DYN%VEL(j,i)=DYN%VEL(j,i)+CMASS(j,1)
            ENDDO
          ENDDO

!          CALL init_velocities(VEL_tmp,T_INFO,DYN,LATT_CUR,NDEGREES_OF_FREEDOM_)
          CALL boltzmann_velocities(T_INFO,DYN,LATT_CUR,VEL_tmp)

          DYN%VEL=DYN%VEL+2*(DYN%D2C-spring_accel)
          CALL KINETIC_E(EKIN,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,&
                       & DYN%POTIM,LATT_CUR%A,DYN%VEL)
          T_old = 2*EKIN/BOLKEV/NDEGREES_OF_FREEDOM

          CMASS=0._q
          random_counter=0
          IF (DYN%INIT==1) THEN
            DO i=1,T_INFO%NIONS
              CALL RANE_ION(dummy)
              IF (dummy<ANDERSEN_PROB) THEN
                random_counter=random_counter+1
                DYN%VEL(:,i)=VEL_tmp(:,i)
              ENDIF
            ENDDO
          ENDIF

          CMASS=0._q
          CALL GIVE_CMASS(T_INFO,DYN%VEL,CMASS)
          DO i=1,T_INFO%NIONS
            DO j=1,3
              DYN%VEL(j,i)=DYN%VEL(j,i)-CMASS(j,1)
            ENDDO
          ENDDO
         
          !c velocity for the step t-1/2dt
          VEL_tmp=DYN%VEL


          DYN%POSION=DYN%POSIOC+DYN%VEL
          ECONST=0._q
          IF (iconst0>0 .OR. iconst2>0) THEN
            C_FORCE=0._q
            CALL Lproj_out_step(T_INFO,INFO,DYN,LATT_CUR,IO,C_FORCE,C_FORCE_L,ICOORDINATES, &
            &  g_io,3*T_INFO%NIONS,AMASS,LATT_CUR%A,ECONST,iconst0,iconst2)
            DYN%POSION=DYN%POSION+C_FORCE
            DYN%VEL=DYN%VEL+C_FORCE
          ENDIF

          !!c velocity for the step t
          !!VEL_tmp=(VEL_tmp+DYN%VEL)/2
          !!CALL KINETIC_E(EKIN,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,&
          !!             & DYN%POTIM,LATT_CUR%A,VEL_tmp)

          !c due to compatibility reasons we still compute
          !c kinetic energy and T using velocities from the step t+1/2
          !c but in the future we should use VEL_tmp (velocity for step t) instead
          CALL KINETIC_E(EKIN,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,&
                       & DYN%POTIM,LATT_CUR%A,DYN%VEL)
          TEIN=2*EKIN/BOLKEV/NDEGREES_OF_FREEDOM

          IF (IO%IU6>0) CALL ENERGY_OUT(g_io,DYN,T_INFO,EKIN,EPOT,ECONST,0._q,0._q,0._q,counter)


          !c fictitious variables
          FACT=(DYN%POTIM**2)*EVTOJ/AMTOKG *1E-10_q

          IF (hills%number>0) THEN
            hills%velocity=hills%velocity+FACT*(spring_forceS-hills%force-penalty%force)/hills%mass

            !c current temperature of fict. particles
            ti=0._q
            DO i=1, hills%number
              ti=ti+(hills%velocity(i)/DYN%POTIM)**2*hills%mass(i)
            ENDDO
            ti=ti*AMTOKG/BOLKEV/EVTOJ*1E10_q/hills%number
            ti_old=ti

            !c apply Andersen thermostat for fictitious particles
            FACT= (AMTOKG/EVTOJ)*(1E-10_q/(DYN%POTIM*1E-15_q))**2
            DO i=1, hills%number
              CALL RANE_ION(dummy)
              IF (dummy<hills%andersen_prob) THEN
                BMP=SQRT(hills%temperature*BOLKEV/(hills%mass(i)*FACT))
                hills%velocity(i)=boltzmann_distribution(0._q,BMP)
              ENDIF
            ENDDO

            !c new temperature of fict. particles
            ti=0._q
            DO i=1, hills%number
              ti=ti+(hills%velocity(i)/DYN%POTIM)**2*hills%mass(i)
            ENDDO
            ti=ti*AMTOKG/BOLKEV/EVTOJ*1E10_q/hills%number

            hills%position=hills%position+hills%velocity

            !c reflecting walls
            DO i=1,hills%number
              IF (hills%position(i)<=penalty%wall(i,1)) THEN
                hills%position(i)=penalty%wall(i,1)+(penalty%wall(i,1)-hills%position(i))
                hills%velocity(i)=-hills%velocity(i)
              ENDIF
              IF (hills%position(i)>=penalty%wall(i,2)) THEN
                hills%position(i)=penalty%wall(i,2)+(penalty%wall(i,2)-hills%position(i))
                hills%velocity(i)=-hills%velocity(i)
              ENDIF
            ENDDO

#ifdef MPI
            CALL  M_bcast_d(WDES%COMM, hills%position, hills%number)
#endif
            IF (IO%IU6>0) THEN
              j=0
              write(g_io%REPORT,ADVANCE='YES',FMT='(/,2X,A13)') '>Metadynamics'
              write(g_io%REPORT,ADVANCE='YES',FMT='(10X,A9,X,A9)') &
                  & "S(q)","s"
              DO i=1,ICOORDINATES%NUMINTERNALS
                IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==5 .OR. ICOORDINATES%COORDSTRUCT(i)%STATUS==6) THEN
                  j=j+1                
                  write(g_io%REPORT,ADVANCE='YES',FMT='(3X,A6,X,F9.5,X,F9.5)') &
                  & 'fic_p>',ICOORDINATES%COORDSTRUCT(i)%VALUE,hills%position(j)
                ENDIF
              ENDDO     
            ENDIF
          ENDIF

          IF (IO%IU6>0) THEN
            write(g_io%REPORT,ADVANCE='YES',FMT='(/,2X,A12)') '>Temperature'
            write(g_io%REPORT,ADVANCE='YES',FMT='(10X,A12,X,A12,X,A12,X,A12)') &
                & "T_sim","T_inst","t_sim","t_inst"
            write(g_io%REPORT,ADVANCE='YES',FMT='(3X,A6,X,F12.5,X,F12.5,X,F12.5,X,F12.5)') &
                & 'tmprt>', DYN%TEMP,T_old,hills%temperature,ti_old
            !write(g_io%REPORT,ADVANCE='YES',FMT='(3X,A6,X,F12.5,X,F12.5,X,F12.5,X,F12.5)') &
            !    & 'tmprt>',T_old,T_new,ti_old,ti
             
          ENDIF

          IF (IO%IU6>0) THEN
            write(g_io%REPORT,FMT='(/,2X,A32,X,I16)') '>Thermostat, num. of collisions:',random_counter
          ENDIF

          DISMAX=0
          DO I=1,T_INFO%NIONS
             DISMAX=MAX(DISMAX, &
             &(DYN%VEL(1,I))**2+ &
             &(DYN%VEL(2,I))**2+ &
             &(DYN%VEL(3,I))**2)
          ENDDO
          DISMAX=SQRT(DISMAX)
          EPS=0.0_q
          ES =0.0_q
        END SUBROUTINE HILLS_METHOD

        SUBROUTINE HILLS_METHOD_nose(DYN,ICOORDINATES,hills,penalty,T_INFO,INFO,LATT_CUR, &
          &      EKIN,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,IO,EPOT,ANDERSEN_PROB,g_io, & 
          &      WDES,iconst0,iconst2,TEIN)
        !c NVT metadynamics (Nose-Hover thermostat) with fictitious coordinates 
          TYPE(coordstructure) :: ICOORDINATES     !c working and reference coords.
          TYPE(dynamics) :: DYN
          TYPE(type_info) :: T_INFO
          TYPE (info_struct) INFO
          TYPE(latt) :: LATT_CUR
          TYPE (in_struct) :: IO
          TYPE(gadget_io) :: g_io
          TYPE(hills_data) :: hills
          TYPE (wavedes)  ::   WDES
          TYPE(penalty_data) :: penalty
          REAL(q) :: C_FORCE(3,T_INFO%NIONS)
          REAL(q) :: EKIN,EPS,ES,DISMAX,TEMPER,EPOT,ti,ti_old,T_old,T_new,ECONST,TEIN
          INTEGER :: i,j,ii,jj,NDEGREES_OF_FREEDOM,random_counter
          INTEGER,SAVE :: counter=0
          INTEGER :: iconst0,iconst2        
          REAL(q) :: dummy,FACT,BMP,ANDERSEN_PROB
          REAL(q) :: CMASS(3,1)
          REAL(q) :: VEL_tmp(3,T_INFO%NIONS)
          REAL(q) :: spring_accel(3,T_INFO%NIONS)
          REAL(q),SAVE,ALLOCATABLE :: spring_forceS(:)
          REAL    :: SQQ,hEKIN,UT,UL
          REAL(q) :: AMASS,C_FORCE_L(3,3)

          counter=counter+1
          UL =1E-10_q*LATT_CUR%ANORM(1)
          UT =DYN%POTIM*1E-15_q
          SQQ=DYN%SMASS*(AMTOKG/EVTOJ)*(UL/UT)**2

          AMASS=0._q
          ti_old=0.

          IF (DYN%INIT==0) THEN
             DYN%SNOSE=0._q
          ENDIF

          IF (counter==1) THEN
            hills%SNOSE=0.0
            ALLOCATE(spring_forceS(hills%number))
            spring_forceS=0._q
 
            j=0
            DO i=1,ICOORDINATES%NUMINTERNALS
                IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==5 .OR. ICOORDINATES%COORDSTRUCT(i)%STATUS==6) THEN
                  j=j+1
                  !!!hills%position(j)=ICOORDINATES%COORDSTRUCT(i)%VALUE
                  IF (hills%position(j)<=(penalty%wall(j,1)+hills%gauss(1)%width)) THEN
                    hills%position(j)=hills%position(j)+2*hills%gauss(1)%width
                  ENDIF
                  IF (hills%position(j)>=(penalty%wall(j,2)-hills%gauss(1)%width)) THEN
                    hills%position(j)=hills%position(j)-2*hills%gauss(1)%width
                  ENDIF
                ENDIF
            ENDDO
          ENDIF 
          
          !c calculate additional potentaials (spring, bias, and penalty)
          CALL hills_spring_force(T_INFO,LATT_CUR,DYN,hills,ICOORDINATES,DYN%POSIOC,spring_forceS,spring_accel)
          IF (counter .GE. hills%bin) THEN
            CALL hills_bias_potential(hills,penalty,ICOORDINATES,DYN,counter,g_io,IO)
          ENDIF 
          CALL penalty_potential(penalty,hills,ICOORDINATES)

          !c real variables +  nose-hoover thermostat
          CALL KINETIC_E(EKIN,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,&
                       & DYN%POTIM,LATT_CUR%A,DYN%VEL)
          T_old = 2*EKIN/BOLKEV/(NDEGREES_OF_FREEDOM)
          
          IF (counter>1) DYN%VEL=(1._q-DYN%SNOSE(2)/2._q)*DYN%VEL+DYN%D2C-spring_accel
          EPS=0.5_q*(DYN%SNOSE(2)**2)*SQQ
          ES =NDEGREES_OF_FREEDOM*BOLKEV*DYN%TEMP*(DYN%SNOSE(1))
          DYN%VEL=(DYN%VEL+DYN%D2C-spring_accel)/(1._q+DYN%SNOSE(2)/2._q)
          CALL KINETIC_E(EKIN,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,&
                       & DYN%POTIM,LATT_CUR%A,DYN%VEL)
          T_new = 2*EKIN/BOLKEV/(NDEGREES_OF_FREEDOM)
          TEIN=T_new
          DYN%SNOSE(3)=DYN%SNOSE(2)
          DYN%SNOSE(2)=DYN%SNOSE(2)+(2*EKIN-NDEGREES_OF_FREEDOM*BOLKEV*DYN%TEMP)/SQQ
          DYN%SNOSE(4)=DYN%SNOSE(1)
          DYN%SNOSE(1)=DYN%SNOSE(1)+(DYN%SNOSE(2)+DYN%SNOSE(3))/2._q

          !write(*,*) 'spring_accel',spring_accel
          DYN%POSION=DYN%POSIOC+DYN%VEL
          ECONST=0._q
          IF (ICOORDINATES%NUMINTERNALS>0) THEN
            C_FORCE=0._q
            CALL Lproj_out_step(T_INFO,INFO,DYN,LATT_CUR,IO,C_FORCE,C_FORCE_L,ICOORDINATES,&
            &  g_io,3*T_INFO%NIONS,AMASS,LATT_CUR%A,ECONST,iconst0,iconst2)
            DYN%POSION=DYN%POSION+C_FORCE
            DYN%VEL=DYN%VEL+C_FORCE
          ENDIF

          IF (IO%IU6>0) CALL ENERGY_OUT(g_io,DYN,T_INFO,EKIN,EPOT,ECONST,0._q,EPS,ES,counter)
         
          IF (hills%number>0) THEN
            !c fictitious variables
            !c current temperature of fict. particles
            ti=0._q
            DO i=1, hills%number
              ti=ti+(hills%velocity(i))**2*hills%mass(i)
            ENDDO
            ti=ti*AMTOKG/BOLKEV/EVTOJ*1E10_q/hills%number/(DYN%POTIM**2)
            ti_old=ti

            UL=1E-10_q
            UT=DYN%POTIM*1E-15_q
            IF (hills%andersen_prob > 1E-5_q) THEN
              !c apply Andersen thermostat for fictitious particles
              FACT=(DYN%POTIM**2)*EVTOJ/AMTOKG *1E-10_q
              hills%velocity=hills%velocity+FACT*(spring_forceS-hills%force-penalty%force)/hills%mass
              FACT= (AMTOKG/EVTOJ)*(1E-10_q/(DYN%POTIM*1E-15_q))**2
              DO i=1, hills%number
                CALL RANE_ION(dummy)
                IF (dummy<hills%andersen_prob) THEN
                  BMP=SQRT(hills%temperature*BOLKEV/(hills%mass(i)*FACT))
                  hills%velocity(i)=boltzmann_distribution(0._q,BMP)
                ENDIF
              ENDDO
              hills%position=hills%position+hills%velocity
            ELSE
              !c apply Nose-Hover thermostat for fictitious particles
              FACT=(DYN%POTIM**2)*EVTOJ/AMTOKG *1E-10_q
              IF (counter>1)  hills%velocity=(1._q-hills%SNOSE(2)/2._q)*hills%velocity+&
                                         FACT*(spring_forceS-hills%force-penalty%force)/hills%mass/2.0

              hills%velocity=(hills%velocity+FACT*(spring_forceS-hills%force-penalty%force)/hills%mass/2.0)/&
                             (1._q+hills%SNOSE(2)/2._q)
              hills%position=hills%position+hills%velocity

              hEKIN=0._q
              DO i=1, hills%number
                hEKIN=hEKIN+hills%velocity(i)**2*hills%mass(i)
              ENDDO
          
              !!!!!!!!!!!!!:check:!!!!!!!!1 hEKIN=hEKIN*AMTOKG/EVTOJ*(UL/UT)**2
              hEKIN=0.5*hEKIN*AMTOKG/EVTOJ*(UL/UT)**2

              FACT=(DYN%POTIM**2)*EVTOJ/AMTOKG *1E-10_q
              !FACT=EVTOJ/AMTOKG *UT**2/UL
              hills%SNOSE(3)=hills%SNOSE(2)
              hills%SNOSE(2)=hills%SNOSE(2)+FACT*(2*hEKIN-hills%number*BOLKEV*hills%temperature)/(hills%SQQ)
              !!!hills%SNOSE(2)=hills%SNOSE(2)+(2*hEKIN-hills%number*BOLKEV*hills%temperature)/(hills%SQQ)
              hills%SNOSE(4)=hills%SNOSE(1)
              hills%SNOSE(1)=hills%SNOSE(1)+(hills%SNOSE(2)+hills%SNOSE(3))/2._q 
            ENDIF

            !c new temperature of fict. particles
            ti=0._q
            DO i=1, hills%number
              ti=ti+(hills%velocity(i))**2*hills%mass(i)
            ENDDO
            ti=ti*AMTOKG/BOLKEV/EVTOJ*1E10_q/hills%number/(DYN%POTIM**2)
          

            !c reflecting walls
            DO i=1,hills%number
              IF (hills%position(i)<=penalty%wall(i,1)) THEN
                hills%position(i)=penalty%wall(i,1)+(penalty%wall(i,1)-hills%position(i))
                hills%velocity(i)=-hills%velocity(i)
              ENDIF
              IF (hills%position(i)>=penalty%wall(i,2)) THEN
                hills%position(i)=penalty%wall(i,2)+(penalty%wall(i,2)-hills%position(i))
                hills%velocity(i)=-hills%velocity(i)
              ENDIF
            ENDDO

#ifdef MPI
            CALL  M_bcast_d(WDES%COMM, hills%position, hills%number)
#endif

            IF (IO%IU6>0) THEN
              j=0
              write(g_io%REPORT,ADVANCE='YES',FMT='(/,2X,A13)') '>Metadynamics'
              write(g_io%REPORT,ADVANCE='YES',FMT='(10X,A9,X,A9)') &
                  & "S(q)","s"
              DO i=1,ICOORDINATES%NUMINTERNALS
                IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==5 .OR. ICOORDINATES%COORDSTRUCT(i)%STATUS==6) THEN
                  j=j+1                
                  write(g_io%REPORT,ADVANCE='YES',FMT='(3X,A6,X,F9.5,X,F9.5)') &
                  & 'fic_p>',ICOORDINATES%COORDSTRUCT(i)%VALUE,hills%position(j)
                ENDIF
              ENDDO     
            ENDIF
          ENDIF

          IF (IO%IU6>0) THEN
            write(g_io%REPORT,ADVANCE='YES',FMT='(/,2X,A12)') '>Temperature'
            write(g_io%REPORT,ADVANCE='YES',FMT='(10X,A12,X,A12,X,A12,X,A12)') &
                & "T_sim","T_inst","t_sim","t_inst"
            write(g_io%REPORT,ADVANCE='YES',FMT='(3X,A6,X,F12.5,X,F12.5,X,F12.5,X,F12.5)') &
                & 'tmprt>', DYN%TEMP,T_old,hills%temperature,ti_old
            !write(g_io%REPORT,ADVANCE='YES',FMT='(3X,A6,X,F12.5,X,F12.5,X,F12.5,X,F12.5)') &
            !    & 'tmprt>',T_old,T_new,ti_old,ti
             
          ENDIF

          DISMAX=0
          DO I=1,T_INFO%NIONS
             DISMAX=MAX(DISMAX, &
             &(DYN%VEL(1,I))**2+ &
             &(DYN%VEL(2,I))**2+ &
             &(DYN%VEL(3,I))**2)
          ENDDO
          DISMAX=SQRT(DISMAX)
        END SUBROUTINE HILLS_METHOD_nose      
        
        SUBROUTINE hills_spring_force(T_INFO,LATT_CUR,DYN,hills,ICOORDINATES,REF_C,spring_forceS,spring_accel)
          !c metadynamics with fictitious coordinates - calculates the spring force necessary to 
          !c keep values of chosen geometric parameters
          !c close to those for collective variables
          TYPE(type_info) :: T_INFO
          TYPE(latt) :: LATT_CUR
          TYPE(dynamics) :: DYN 
          TYPE(hills_data) :: hills
          TYPE(coordstructure) :: ICOORDINATES
          INTEGER :: i,j,NI,NT
          REAL(q) :: REF_C(3,T_INFO%NIONS)
          REAL(q) :: FACT,deltaS
          REAL(q) :: s_force(3*T_INFO%NIONS)
          REAL(q) :: spring_force(3,T_INFO%NIONS),spring_accel(3,T_INFO%NIONS)
          REAL(q) :: spring_forceS(:)
          REAL(q) :: dummy(3,T_INFO%NIONS)
          REAL(q) :: BMAT(ICOORDINATES%NUMINTERNALS,3*T_INFO%NIONS)
          !REAL(q),ALLOCATABLE :: BMAT(:,:)


          IF (hills%number<1) RETURN

          BMAT=0._q
          CALL BMATRIX(T_INFO,REF_C,REF_C,TRANSPOSE(LATT_CUR%A),ICOORDINATES,BMAT,.FALSE.)
          !CALL ICONST_BMAT(BMAT,ICOORDINATES,T_INFO)
          spring_force=0._q
          spring_accel=0._q
          
          j=0
          DO i=1,ICOORDINATES%NUMINTERNALS
            IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==5 .OR. ICOORDINATES%COORDSTRUCT(i)%STATUS==6) THEN
              j=j+1
              deltaS=ICOORDINATES%COORDSTRUCT(i)%VALUE-hills%position(j)
              IF (ICOORDINATES%COORDSTRUCT(i)%TAG=='A ' .OR. &
                ICOORDINATES%COORDSTRUCT(i)%TAG=='T ') THEN
                DO
                  IF (deltaS>pi) THEN
                    !!!deltaS=2*pi-deltaS
                    deltaS=deltaS-2*pi
                  ELSE
                    EXIT
                  ENDIF
                ENDDO
                DO
                  IF (deltaS<=-pi) THEN
                    deltaS=2*pi+deltaS
                  ELSE
                    EXIT
                  ENDIF
                ENDDO
              ENDIF
              spring_forceS(j)=hills%force_constant(j)*deltaS
              s_force=spring_forceS(j)*BMAT(i,:)
              dummy=0._q
              CALL ONETOTHREE(T_INFO%NIONS,dummy,s_force)
              dummy=MATMUL(LATT_CUR%B,dummy) !!!!!!!!!
              spring_force=spring_force+dummy
            ENDIF
          ENDDO
          FACT=(DYN%POTIM**2)*EVTOJ/AMTOKG *1E-10_q
          NI=1
          DO NT=1,T_INFO%NTYP
            DO NI=NI,T_INFO%NITYP(NT)+NI-1
              IF (T_INFO%LSFOR(1,NI)) spring_accel(1,NI)=spring_force(1,NI)*FACT/2/T_INFO%POMASS(NT)
              IF (T_INFO%LSFOR(2,NI)) spring_accel(2,NI)=spring_force(2,NI)*FACT/2/T_INFO%POMASS(NT)
              IF (T_INFO%LSFOR(3,NI)) spring_accel(3,NI)=spring_force(3,NI)*FACT/2/T_INFO%POMASS(NT)
            ENDDO 
          ENDDO

          !!spring_accel=MATMUL(LATT_CUR%B,spring_accel) !!!!!!!!!
          spring_accel=MATMUL(TRANSPOSE(LATT_CUR%B),spring_accel)
          !DEALLOCATE(BMAT)
          !write(*,*) 'spring_forceS',spring_forceS
          !write(*,*) 'spring_accel',spring_accel
        END SUBROUTINE hills_spring_force

        SUBROUTINE hills_bias_potential(hills,penalty,ICOORDINATES,DYN,step,g_io,IO)
          !c here we calculate the history-dependent bias potential and forces
          TYPE(hills_data) :: hills
          TYPE(penalty_data) :: penalty
          TYPE(gadget_io) :: g_io
          TYPE(coordstructure) :: ICOORDINATES
          TYPE(dynamics) :: DYN 
          TYPE (in_struct) :: IO
          REAL(q) :: snorm,w1,w2
          INTEGER,SAVE :: counter=0            !c meta-time
          INTEGER  ::step
          INTEGER :: i,j,k,dmn1,dmn2
          INTEGER,SAVE :: last_upgrade=0
          REAL(q) :: dummy,ti

          IF (step<hills%bin .OR. hills%number<1) RETURN
          !OPEN(UNIT=g_io%STRUCTINPUT,FILE=DIR_APP(1:DIR_LEN)//'HILLSPOT',STATUS='UNKNOWN',POSITION='APPEND')
          OPEN(UNIT=g_io%STRUCTINPUT,FILE='HILLSPOT',STATUS='UNKNOWN',POSITION='APPEND')
          
          IF (MOD(step,hills%bin)==0) THEN

            IF (counter==0) THEN
              !c add another hill
              last_upgrade=step
              counter=counter+1
              hills%gauss(counter)%position(:)=hills%position(:)
              DO j=1,hills%number
                IF (g_io%STRUCTINPUT>=0 .AND. (IO%IU0>0)) THEN
                  write(g_io%STRUCTINPUT,ADVANCE='NO',FMT='(X,F9.5)') hills%gauss(counter)%position(j)
                ENDIF
              ENDDO

              IF ((g_io%STRUCTINPUT>=0) .AND. (IO%IU0>0)) THEN
                write(g_io%STRUCTINPUT,ADVANCE='NO',FMT='(X,F9.5)') hills%gauss(counter)%high
                write(g_io%STRUCTINPUT,ADVANCE='YES',FMT='(X,F9.5)') hills%gauss(counter)%width
              ENDIF

            ELSE IF ((SQRT(SUM((hills%gauss(counter)%position(:)-hills%position(:))**2))>& 
            &    hills%stride*hills%gauss(counter)%width) & 
            &    .OR. (MOD((step-last_upgrade),hills%maxstride*hills%bin)==0)) THEN
              !c add another hill
              last_upgrade=step
              counter=counter+1
              hills%gauss(counter)%position(:)=hills%position(:)
              DO j=1,hills%number
                IF (g_io%STRUCTINPUT>=0 .AND. (IO%IU0>0)) THEN
                  write(g_io%STRUCTINPUT,ADVANCE='NO',FMT='(X,F9.5)') hills%gauss(counter)%position(j)
                ENDIF
              ENDDO
              
              !c if variable width is desired  - similarly as in the original metadynamics
              !c implementation
              IF (hills%variable_width) THEN
                w1=(hills%gauss(counter)%width)**2
                w2=SUM((hills%gauss(counter)%position(:)-hills%gauss(counter-1)%position(:))**2)
                hills%gauss(counter)%width=SQRT(w1*w2/(w1+w2))
              ENDIF

              IF ((g_io%STRUCTINPUT>=0) .AND. (IO%IU0>0)) THEN
                write(g_io%STRUCTINPUT,ADVANCE='NO',FMT='(X,F9.5)') hills%gauss(counter)%high
                write(g_io%STRUCTINPUT,ADVANCE='YES',FMT='(X,F9.5)') hills%gauss(counter)%width
              ENDIF
            
            ENDIF 
          ENDIF
          
          hills%potential=0._q
          hills%force=0._q
          DO i=1,counter   !c sum over all bins
            snorm=0._q
            !c sum up the potential for the current bin
            k=0
            DO j=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(j)%STATUS==5 .OR. ICOORDINATES%COORDSTRUCT(j)%STATUS==6) THEN
                k=k+1
                IF (ICOORDINATES%COORDSTRUCT(j)%STATUS==5) snorm=snorm+(hills%position(k)-hills%gauss(i)%position(k))**2
              ENDIF
            ENDDO

            !c sum up the forces for the current bin 
            k=0
            DO j=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(j)%STATUS==5 .OR. ICOORDINATES%COORDSTRUCT(j)%STATUS==6) THEN
                k=k+1
                IF (ICOORDINATES%COORDSTRUCT(j)%STATUS==5) THEN
                  hills%force(k)=hills%force(k) &
                  & -(hills%position(k)-hills%gauss(i)%position(k))*hills%gauss(i)%high/(hills%gauss(i)%width**2)* &
                  & exp(-snorm/2._q/hills%gauss(i)%width**2)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
          CLOSE(g_io%STRUCTINPUT)
        END SUBROUTINE hills_bias_potential

        SUBROUTINE hills_bias_direct(T_INFO,LATT_CUR,DYN,internal_force,ICOORDINATES,REF_C,penalty_accel,penalty_accel_L,IBDIM,AMASS)
          !c bias potential for metadynamics without fictitious particles
          TYPE(type_info) :: T_INFO
          TYPE(latt) :: LATT_CUR
          TYPE(dynamics) :: DYN 
          TYPE(coordstructure) :: ICOORDINATES
          INTEGER :: i,j,k,NI,NT
          INTEGER :: IBDIM
          REAL(q) :: REF_C(3,T_INFO%NIONS)
          REAL(q) :: FACT
          REAL(q) :: internal_force(:)
          !REAL(q) :: hspring_force(:)
          REAL(q) :: p_force(3*T_INFO%NIONS)
          REAL(q) :: p_force_L(9)
          REAL(q) :: penalty_force(3,T_INFO%NIONS),penalty_accel(3,T_INFO%NIONS)
          REAL(q) :: penalty_force_L(3,3),penalty_accel_L(3,3)
          REAL(q) :: dummy(3,T_INFO%NIONS),dummy_L(3,3)
          REAL(q) :: BMAT(ICOORDINATES%NUMINTERNALS,IBDIM)
          REAL(q) :: AMASS !c mass for lattice DOFs

          BMAT=0._q
          CALL BMATRIX(T_INFO,REF_C,REF_C,TRANSPOSE(LATT_CUR%A),ICOORDINATES,BMAT,.FALSE.)
          !CALL ICONST_BMAT(BMAT,ICOORDINATES,T_INFO)

          penalty_force=0._q
          penalty_accel=0._q

          penalty_force_L=0._q
          penalty_accel_L=0._q

          j=0;k=0
          DO i=1,ICOORDINATES%NUMINTERNALS
            IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==5) THEN
              j=j+1
              p_force=internal_force(j)*BMAT(i,1:3*T_INFO%NIONS)
              dummy=0._q
              CALL ONETOTHREE(T_INFO%NIONS,dummy,p_force)
              dummy=MATMUL(LATT_CUR%B,dummy) !!!!!!!!!
              penalty_force=penalty_force+dummy

              !c lattice dynamics
              IF (IBDIM==3*T_INFO%NIONS+9) THEN
                p_force_L=internal_force(j)*BMAT(i,3*T_INFO%NIONS+1:3*T_INFO%NIONS+9)
                dummy_L=0._q
                CALL ONETOTHREE(3,dummy_L,p_force_L)
                !!!dummy_L=TRANSPOSE(dummy_L)
                penalty_force_L=penalty_force_L+dummy_L
              ENDIF

            ENDIF
           
!             !c harmonic spring contribution 
!             IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==8) THEN
!               k=k+1
!               p_force=hspring_force(k)*BMAT(i,1:3*T_INFO%NIONS)
!               dummy=0._q
!               CALL ONETOTHREE(T_INFO%NIONS,dummy,p_force)
!               dummy=MATMUL(LATT_CUR%B,dummy) !!!!!!!!!
!               penalty_force=penalty_force+dummy
! 
!               !c lattice dynamics
!               IF (IBDIM==3*T_INFO%NIONS+9) THEN
!                 p_force_L=hspring_force(k)*BMAT(i,3*T_INFO%NIONS+1:3*T_INFO%NIONS+9)
!                 dummy_L=0._q
!                 CALL ONETOTHREE(3,dummy_L,p_force_L)
!                 dummy_L=TRANSPOSE(dummy_L)
!                 penalty_force_L=penalty_force_L+dummy_L
!               ENDIF
!             ENDIF

          ENDDO


          FACT=(DYN%POTIM**2)*EVTOJ/AMTOKG *1E-10_q
          NI=1
          DO NT=1,T_INFO%NTYP
            DO NI=NI,T_INFO%NITYP(NT)+NI-1
              IF (T_INFO%LSFOR(1,NI)) penalty_accel(1,NI)=penalty_force(1,NI)*FACT/2/T_INFO%POMASS(NT)
              IF (T_INFO%LSFOR(2,NI)) penalty_accel(2,NI)=penalty_force(2,NI)*FACT/2/T_INFO%POMASS(NT)
              IF (T_INFO%LSFOR(3,NI)) penalty_accel(3,NI)=penalty_force(3,NI)*FACT/2/T_INFO%POMASS(NT)
            ENDDO 
          ENDDO
          
          penalty_accel=MATMUL(TRANSPOSE(LATT_CUR%B),penalty_accel) 
          !penalty_accel=MATMUL(LATT_CUR%B,penalty_accel) 

          IF (IBDIM==3*T_INFO%NIONS+9) THEN
            penalty_accel_L=penalty_force_L*FACT/2/AMASS    
          ENDIF

        END SUBROUTINE hills_bias_direct

        SUBROUTINE harmonic_spring_accel(istatus,T_INFO,LATT_CUR,DYN,hspring_force,ICOORDINATES,REF_C,penalty_accel,penalty_accel_L,IBDIM,AMASS)
          !c bias potential for metadynamics without fictitious particles
          TYPE(type_info) :: T_INFO
          TYPE(latt) :: LATT_CUR
          TYPE(dynamics) :: DYN 
          TYPE(coordstructure) :: ICOORDINATES
          INTEGER :: istatus !c 4,8
          INTEGER :: i,j,k,NI,NT
          INTEGER :: IBDIM
          REAL(q) :: REF_C(3,T_INFO%NIONS)
          REAL(q) :: FACT
          REAL(q) :: hspring_force(:)
          REAL(q) :: p_force(3*T_INFO%NIONS)
          REAL(q) :: p_force_L(9)
          REAL(q) :: penalty_force(3,T_INFO%NIONS),penalty_accel(3,T_INFO%NIONS)
          REAL(q) :: penalty_force_L(3,3),penalty_accel_L(3,3)
          REAL(q) :: dummy(3,T_INFO%NIONS),dummy_L(3,3)
          REAL(q) :: BMAT(ICOORDINATES%NUMINTERNALS,IBDIM)
          REAL(q) :: AMASS !c mass for lattice DOFs

          BMAT=0._q
          CALL BMATRIX(T_INFO,REF_C,REF_C,TRANSPOSE(LATT_CUR%A),ICOORDINATES,BMAT,.FALSE.)
          !CALL ICONST_BMAT(BMAT,ICOORDINATES,T_INFO)

          penalty_force=0._q
          penalty_accel=0._q

          penalty_force_L=0._q
          penalty_accel_L=0._q

          j=0;k=0
          DO i=1,ICOORDINATES%NUMINTERNALS           
            !c harmonic spring contribution 
            IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==istatus) THEN
              k=k+1
              p_force=hspring_force(k)*BMAT(i,1:3*T_INFO%NIONS)
              dummy=0._q
              CALL ONETOTHREE(T_INFO%NIONS,dummy,p_force)
              dummy=MATMUL(LATT_CUR%B,dummy) !!!!!!!!!
              penalty_force=penalty_force+dummy

              !c lattice dynamics
              IF (IBDIM==3*T_INFO%NIONS+9) THEN
                p_force_L=hspring_force(k)*BMAT(i,3*T_INFO%NIONS+1:3*T_INFO%NIONS+9)
                dummy_L=0._q
                CALL ONETOTHREE(3,dummy_L,p_force_L)
                !!!dummy_L=TRANSPOSE(dummy_L)
                penalty_force_L=penalty_force_L+dummy_L
              ENDIF
            ENDIF

          ENDDO


          FACT=(DYN%POTIM**2)*EVTOJ/AMTOKG *1E-10_q
          NI=1
          DO NT=1,T_INFO%NTYP
            DO NI=NI,T_INFO%NITYP(NT)+NI-1
              IF (T_INFO%LSFOR(1,NI)) penalty_accel(1,NI)=penalty_force(1,NI)*FACT/2/T_INFO%POMASS(NT)
              IF (T_INFO%LSFOR(2,NI)) penalty_accel(2,NI)=penalty_force(2,NI)*FACT/2/T_INFO%POMASS(NT)
              IF (T_INFO%LSFOR(3,NI)) penalty_accel(3,NI)=penalty_force(3,NI)*FACT/2/T_INFO%POMASS(NT)
            ENDDO 
          ENDDO
          
          penalty_accel=MATMUL(TRANSPOSE(LATT_CUR%B),penalty_accel) 
          !penalty_accel=MATMUL(LATT_CUR%B,penalty_accel) 

          IF (IBDIM==3*T_INFO%NIONS+9) THEN
            penalty_accel_L=penalty_force_L*FACT/2/AMASS    
          ENDIF

        END SUBROUTINE harmonic_spring_accel

        SUBROUTINE penalty_potential(penalty,hills,ICOORDINATES)
          !c user provided bias potential
          TYPE(hills_data) :: hills
          TYPE(penalty_data) :: penalty
          TYPE(coordstructure) :: ICOORDINATES
          REAL(q) :: pnorm
          INTEGER :: i,j,k
 
          IF (hills%number<1) RETURN

          penalty%force=0._q
          DO i=1,penalty%number  
            pnorm=0._q !!!
            k=0
            DO j=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(j)%STATUS==5 .OR. ICOORDINATES%COORDSTRUCT(j)%STATUS==6) THEN
                k=k+1
                 IF (ICOORDINATES%COORDSTRUCT(j)%STATUS==5) &
                 &  pnorm=pnorm+(hills%position(k)-penalty%gauss(i)%position(k))**2
              ENDIF
            ENDDO

            !c sum up the forces 
            k=0
            DO j=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(j)%STATUS==5 .OR. ICOORDINATES%COORDSTRUCT(j)%STATUS==6) THEN
                k=k+1
                IF (ICOORDINATES%COORDSTRUCT(j)%STATUS==5) THEN
                  penalty%force(k)=penalty%force(k)&   !!!!
                  & -(hills%position(k)-penalty%gauss(i)%position(k))* &
                  & penalty%gauss(i)%high/(penalty%gauss(i)%width**2)* &
                  & exp(-pnorm/2._q/penalty%gauss(i)%width**2)   !!!!
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        END SUBROUTINE penalty_potential

        SUBROUTINE fermi_bias_force(iconst4,fbias_r0,fbias_d,fbias_a,fbias_force,ICOORDINATES)
          !c forces due to harmonic springs attached to internal coordinates
          TYPE(coordstructure) :: ICOORDINATES
          INTEGER :: i,j,k
          INTEGER :: iconst4
          REAL(q) :: fbias_r0(iconst4),fbias_d(iconst4),fbias_a(iconst4),fbias_force(iconst4)
          REAL(q) :: deltaS

          fbias_force=0._q
          k=0
          DO i=1,ICOORDINATES%NUMINTERNALS
            IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==4) THEN
              k=k+1
              deltaS=EXP(-fbias_D(k)*(ICOORDINATES%COORDSTRUCT(i)%VALUE/fbias_R0(k)-1.))
              fbias_force(k)=fbias_A(k)*fbias_D(k)/fbias_R0(k)*deltaS/((1+deltaS)*(1+deltaS))
            ENDIF
          ENDDO
        END SUBROUTINE fermi_bias_force

        SUBROUTINE harmonic_spring_force(iconst8,hspring_R0,hspring_K,hspring_force,ICOORDINATES)
          !c forces due to harmonic springs attached to internal coordinates
          TYPE(coordstructure) :: ICOORDINATES
          INTEGER :: i,j,k
          INTEGER :: iconst8
          REAL(q) :: hspring_R0(iconst8),hspring_K(iconst8),hspring_force(iconst8)
          REAL(q) :: deltaS
 
          hspring_force=0._q
          k=0
          DO i=1,ICOORDINATES%NUMINTERNALS
            IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==8) THEN
              k=k+1
              deltaS=(ICOORDINATES%COORDSTRUCT(i)%VALUE-hspring_R0(k))
              !c take care of periodicity in definition of angle and torsion angle
              IF (ICOORDINATES%COORDSTRUCT(i)%TAG=='A ' .OR. &
              &   ICOORDINATES%COORDSTRUCT(i)%TAG=='T ') THEN
                DO
                  IF (deltaS>pi) THEN
                    !!!deltaS=2*pi-deltaS
                    deltaS=deltaS-2*pi
                  ELSE
                    EXIT
                  ENDIF
                ENDDO
                DO
                  IF (deltaS<=-pi) THEN
                    deltaS=2*pi+deltaS
                  ELSE
                    EXIT
                  ENDIF
                ENDDO
              ENDIF
              hspring_force(k)=hspring_force(k)+hspring_K(k)*deltaS
              !hspring_force(k)=hspring_K(k)*deltaS
            ENDIF
          ENDDO
   
        END SUBROUTINE harmonic_spring_force

        SUBROUTINE HILLS_METHOD_DIRECT(DYN,ICOORDINATES,hills,penalty,T_INFO,INFO,LATT_CUR, &
          &      EKIN,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,IO,EPOT,ANDERSEN_PROB,&
          &      g_io,iconst0,iconst2,iconst4,iconst5,iconst8,hspring_R0, & 
          &      hspring_K,hspring_force,fbias_R0,fbias_D,fbias_A,TEIN)
        !c NVT metadynamics (Andersen thermostat) without fictitious particles 
          TYPE(coordstructure) :: ICOORDINATES     !c working and reference coords.
          TYPE(dynamics) :: DYN
          TYPE(type_info) :: T_INFO
          TYPE (info_struct) INFO
          TYPE(latt) :: LATT_CUR
          TYPE (in_struct) :: IO
          TYPE(gadget_io) :: g_io
          TYPE(hills_data) :: hills
          TYPE(penalty_data) :: penalty
          REAL(q) :: C_FORCE(3,T_INFO%NIONS)
          REAL(q) :: EKIN,EPS,ES,DISMAX,TEMPER,EPOT,T_old,T_new,ECONST,TEIN
          INTEGER :: i,j,ii,jj,NDEGREES_OF_FREEDOM,random_counter
          INTEGER,SAVE :: counter=0
          INTEGER :: iconst0,iconst2,iconst4,iconst5,iconst8        
          REAL(q) :: dummy,BMP,ANDERSEN_PROB
          REAL(q) :: CMASS(3,1)
          REAL(q) :: VEL_tmp(3,T_INFO%NIONS),VEL_init(3,T_INFO%NIONS)
          REAL(q) :: hills_accel(3,T_INFO%NIONS),penalty_accel(3,T_INFO%NIONS),hspring_accel(3,T_INFO%NIONS)
          REAL(q) :: hills_accel_L(3,3),penalty_accel_L(3,3),hspring_accel_L(3,3)
          REAL(q) :: fbias_accel(3,T_INFO%NIONS),fbias_accel_L(3,3)
          REAL(q) :: AMASS,C_FORCE_L(3,3)
          !LOGICAL :: LSCALE,LCMASS
          REAL(q) :: hspring_R0(iconst8),hspring_K(iconst8),hspring_force(iconst8)
          REAL(q) :: fbias_R0(iconst4),fbias_D(iconst4),fbias_A(iconst4),fbias_force(iconst4)
          REAL(q) :: TOTAL_MASS
          INTEGER,SAVE :: Ltxyz(3)

          !c add this step to total count
          counter=counter+1

          !c determine if coordinates for some direction
          !c are frozen
          IF (counter==1) THEN
            Ltxyz(:)=1
            DO j=1,3
              DO i=1,T_INFO%NIONS
                IF (.NOT. T_INFO%LSFOR(j,i)) Ltxyz(j)=0
              ENDDO
            ENDDO

            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0) THEN
                IF (ICOORDINATES%COORDSTRUCT(i)%TAG=="X ") Ltxyz(1)=0
                IF (ICOORDINATES%COORDSTRUCT(i)%TAG=="Y ") Ltxyz(2)=0
                IF (ICOORDINATES%COORDSTRUCT(i)%TAG=="Z ") Ltxyz(3)=0
              ENDIF
            ENDDO
          ENDIF

          AMASS=0._q 
          hills_accel_L=0._q;penalty_accel_L=0._q;hspring_accel_L=0._q
          hills_accel=0._q;penalty_accel=0._q;hspring_accel=0._q
          fbias_accel=0._q;fbias_accel_L=0._q
           
          IF (iconst5>0) THEN
            !c update history of fictioous particles and of collerctive variables
            j=0
            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==5) THEN
                j=j+1
                hills%position(j)=ICOORDINATES%COORDSTRUCT(i)%VALUE 
              ENDIF
            ENDDO
          
            !c add new hill
            IF (counter .GE. hills%bin) THEN
              CALL hills_bias_potential(hills,penalty,ICOORDINATES,DYN,counter,g_io,IO)
            ENDIF 

            !c calculate the penalty potential
            CALL penalty_potential(penalty,hills,ICOORDINATES)

            !c compute contributions from bias potentials
            CALL hills_bias_direct(T_INFO,LATT_CUR,DYN,hills%force,ICOORDINATES,DYN%POSIOC,hills_accel,hills_accel_L,3*T_INFO%NIONS,0._q)
            CALL hills_bias_direct(T_INFO,LATT_CUR,DYN,penalty%force,ICOORDINATES,DYN%POSIOC,penalty_accel,penalty_accel_L,3*T_INFO%NIONS,0._q)
            !hills%force=hills%force+penalty%force
            !CALL hills_bias_direct(T_INFO,LATT_CUR,DYN,hills%force,hspring_force,ICOORDINATES,DYN%POSIOC,hills_accel,hills_accel_L,3*T_INFO%NIONS,0._q)
          ENDIF

          !c calculate forces due to Fermi-shaped bias potential
          IF (iconst4>0) THEN
            fbias_force=0.
            CALL fermi_bias_force(iconst4,fbias_r0,fbias_d,fbias_a,fbias_force,ICOORDINATES)
            CALL harmonic_spring_accel(4,T_INFO,LATT_CUR,DYN,fbias_force,ICOORDINATES,DYN%POSIOC,fbias_accel,fbias_accel_L,3*T_INFO%NIONS,0._q)
          ENDIF

          !c calculate forces due to harmonic springs
          IF (iconst8>0) THEN
            CALL harmonic_spring_force(iconst8,hspring_R0,hspring_K,hspring_force,ICOORDINATES)
            CALL harmonic_spring_accel(8,T_INFO,LATT_CUR,DYN,hspring_force,ICOORDINATES,DYN%POSIOC,hspring_accel,hspring_accel_L,3*T_INFO%NIONS,0._q)
          ENDIF

          IF ((iconst4>0) .OR. (iconst8>0)) hspring_accel=hspring_accel+fbias_accel

          !c store velocity for step t-1/2
          VEL_init=DYN%VEL

          !c special treatment of velocities needed only if the likelihood of 
          !c collisions with virtual medium is nonvanishing
          IF (ANDERSEN_PROB>0.) THEN
            !c real variables +  andersen thermostat
            VEL_tmp=0._q
            !!CALL init_velocities(VEL_tmp,T_INFO,DYN,LATT_CUR,NDEGREES_OF_FREEDOM,LSCALE,LCMASS)
        
            ! Add a random velicity that is exactly
            ! distributed like the centre of mass velocity
            TOTAL_MASS = 0.0_q
            DO I = 1, T_INFO%NIONS
              TOTAL_MASS = TOTAL_MASS + T_INFO%POMASS(T_INFO%ITYP(I))
            END DO
            CMASS(:,1) = BOLTZMANN_VELOCITY(DYN, LATT_CUR, TOTAL_MASS)

            !c don't add translation for frozen direction
            DO i=1,3
              IF (Ltxyz(i)==0) CMASS(i,1)=0._q
            ENDDO

            !c add add randomly generated motion of center of mass
            DO i=1,T_INFO%NIONS
               DO j=1,3
                  DYN%VEL(j,i)=DYN%VEL(j,i)+CMASS(j,1)
               ENDDO
            ENDDO
          ENDIF

          !c boltzmann_velocities takes care of frozen atoms (velocities
          !c set to zero)
          CALL boltzmann_velocities(T_INFO,DYN,LATT_CUR,VEL_tmp)

          DYN%VEL=DYN%VEL+2*(DYN%D2C-hills_accel-penalty_accel-hspring_accel)
      
          !c collisions with the virtual medium
          random_counter=0
          IF (ANDERSEN_PROB>0.) THEN
            IF (DYN%INIT==1) THEN
              DO i=1,T_INFO%NIONS
                CALL RANE_ION(dummy)
                IF (dummy<ANDERSEN_PROB) THEN
                  random_counter=random_counter+1
                  DYN%VEL(:,i)=VEL_tmp(:,i)
                ENDIF
              ENDDO
            ENDIF

            CMASS=0._q
            !c GIVE_CMASS takes care of frozen directions (corresponding
            !c component set to zero)
            CALL GIVE_CMASS(T_INFO,DYN%VEL,CMASS)

            !c remove randomly generated velocity of center of mass
            DO i=1,T_INFO%NIONS
              DO j=1,3
                DYN%VEL(j,i)=DYN%VEL(j,i)-CMASS(j,1)
              ENDDO
            ENDDO
          ENDIF

          DYN%POSION=DYN%POSIOC+DYN%VEL
          ECONST=0._q
          IF (iconst0>0 .OR. iconst2>0) THEN
            C_FORCE=0._q
            CALL Lproj_out_step(T_INFO,INFO,DYN,LATT_CUR,IO,C_FORCE,C_FORCE_L,ICOORDINATES,& 
            &   g_io,3*T_INFO%NIONS,AMASS,LATT_CUR%A,ECONST,iconst0,iconst2)
            DYN%POSION=DYN%POSION+C_FORCE
            DYN%VEL=DYN%VEL+C_FORCE
          ENDIF

          !!c estimate velocity for the step t and compute T
          VEL_tmp=0.5*(VEL_init+DYN%VEL)
          CALL KINETIC_E(EKIN,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,&
                       & DYN%POTIM,LATT_CUR%A,VEL_tmp)

          TEIN=2*EKIN/BOLKEV/NDEGREES_OF_FREEDOM

          IF (IO%IU6>0) CALL ENERGY_OUT(g_io,DYN,T_INFO,EKIN,EPOT,ECONST,0._q,0._q,0._q,counter)

          !IF (g_io%REPORT>=0) THEN
          !c output for metadynamics
          IF (iconst5>0) THEN
            IF (IO%IU6>0) THEN
              j=0
              write(g_io%REPORT,ADVANCE='YES',FMT='(/,2X,A13)') '>Metadynamics'
              DO i=1,ICOORDINATES%NUMINTERNALS
                IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==5) THEN
                  j=j+1
                  write(g_io%REPORT,ADVANCE='YES',FMT='(3X,A6,X,F9.5)') &
                  & 'fic_p>',ICOORDINATES%COORDSTRUCT(i)%VALUE
                  !write(*,*) 'test_val',ICOORDINATES%COORDSTRUCT(i)%VALUE,hills%position(j)
                ENDIF
              ENDDO     
            ENDIF
          ENDIF

          IF (IO%IU6>0) THEN
            CALL TEMPERATURE_OUT(g_io,DYN,TEIN)            
          ENDIF

          !IF (g_io%REPORT>=0) THEN
          IF (IO%IU6>0) THEN
            write(g_io%REPORT,FMT='(/,2X,A32,X,I16)') '>Thermostat, num. of collisions:',random_counter
          ENDIF

          DISMAX=0
          DO I=1,T_INFO%NIONS
             DISMAX=MAX(DISMAX, &
             &(DYN%VEL(1,I))**2+ &
             &(DYN%VEL(2,I))**2+ &
             &(DYN%VEL(3,I))**2)
          ENDDO
          DISMAX=SQRT(DISMAX)
          EPS=0.0_q
          ES =0.0_q
        END SUBROUTINE HILLS_METHOD_DIRECT

        !SUBROUTINE RIGID_WALLS_direct()
        !  TYPE(type_info) :: T_INFO
        !  TYPE(latt) :: LATT_CUR
        !  TYPE(dynamics) :: DYN 
        ! ! TYPE(hills_data) :: hills
        !  TYPE(coordstructure) :: ICOORDINATES
        !  INTEGER :: i,j,NI,NT
        !  REAL(q) :: REF_C(3,T_INFO%NIONS)
        !  REAL(q) :: FACT
        !  REAL(q) :: internal_force(:)
        !  REAL(q) :: p_force(3*T_INFO%NIONS)
        !  REAL(q) :: penalty_force(3,T_INFO%NIONS),penalty_accel(3,T_INFO%NIONS)
        !  REAL(q) :: dummy(3,T_INFO%NIONS)
        !  REAL(q),ALLOCATABLE :: BMAT(:,:)
        !
        !  ALLOCATE(BMAT(ICOORDINATES%NUMINTERNALS,3*T_INFO%NIONS))
        !  BMAT=0._q
        !  CALL BMATRIX(T_INFO,REF_C,REF_C,TRANSPOSE(LATT_CUR%A),ICOORDINATES,BMAT,.FALSE.)
        !  CALL ICONST_BMAT(BMAT,ICOORDINATES,T_INFO)
        !
        !  CALL DEAL_XYZ(T_INFO,REF_C,REF_C,LATT_CUR%A,ICOORDINATES)
        !  j=0
        !  DO i=1,ICOORDINATES%NUMINTERNALS
        !    IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==5) THEN
        !      j=j+1
        !      IF (ICOORDINATES%COORDSTRUCT(i)%VALUE<=(penalty%wall(j,1)+hills%gauss(1)%width)) THEN
        !        DYN%POSION=DYN%POSION   SUM(BMAT(i,:)*DYN
        !        hills%velocity(i)=-hills%velocity(i)
        !        hills%position(j)=hills%position(j)+2*hills%gauss(1)%width
        !      ENDIF
        !      IF (hills%position(j)>=(penalty%wall(j,2)-hills%gauss(1)%width)) THEN
        !        hills%position(j)=hills%position(j)-2*hills%gauss(1)%width
        !      ENDIF
        !    ENDIF
        ! ENDDO
        !
        !END SUBROUTINE RIGID_WALLS_direc
       

        SUBROUTINE HILLS_METHOD_DIRECT_nose(DYN,ICOORDINATES,hills,penalty,T_INFO,INFO,LATT_CUR, &
          &      EKIN,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,IO,EPOT,ANDERSEN_PROB,& 
          &      g_io,iconst0,iconst2,iconst4,iconst5,iconst8,hspring_R0,& 
          &      hspring_K,hspring_force,fbias_R0,fbias_D,fbias_A,TEIN)
        !c NVT metadynamics (Nose-Hover thermostat) without fictitious particles 
          TYPE(coordstructure) :: ICOORDINATES     !c working and reference coords.
          TYPE(dynamics) :: DYN
          TYPE(type_info) :: T_INFO
          TYPE (info_struct) INFO
          TYPE(latt) :: LATT_CUR
          TYPE (in_struct) :: IO
          TYPE(gadget_io) :: g_io
          TYPE(hills_data) :: hills
          TYPE(penalty_data) :: penalty
          REAL(q) :: C_FORCE(3,T_INFO%NIONS)
          REAL(q) :: EKIN,EPS,ES,DISMAX,TEMPER,EPOT,ECONST,TEIN
          INTEGER :: i,j,ii,jj,NDEGREES_OF_FREEDOM,random_counter
          INTEGER,SAVE :: counter=0        
          INTEGER :: iconst0,iconst2,iconst4,iconst5,iconst8
          REAL(q) :: dummy,BMP,ANDERSEN_PROB
          REAL    :: SQQ,hEKIN,UT,UL
          REAL(q) :: hills_accel(3,T_INFO%NIONS),penalty_accel(3,T_INFO%NIONS),hspring_accel(3,T_INFO%NIONS)
          REAL(q) :: hills_accel_L(3,3),penalty_accel_L(3,3),hspring_accel_L(3,3)
          REAL(q) :: AMASS,C_FORCE_L(3,3), VEL_tmp(3,T_INFO%NIONS)
          REAL(q) :: hspring_R0(iconst8),hspring_K(iconst8),hspring_force(iconst8)
          REAL(q) :: fbias_R0(iconst4),fbias_D(iconst4),fbias_A(iconst4),fbias_force(iconst4)
          REAL(q) :: fbias_accel(3,T_INFO%NIONS),fbias_accel_L(3,3)

          counter=counter+1
          !DYN%POSIOC=DYN%POSION
          UL =1E-10_q*LATT_CUR%ANORM(1)
          UT =DYN%POTIM*1E-15_q
          IF (DYN%SMASS>0.) THEN 
            SQQ=DYN%SMASS*(AMTOKG/EVTOJ)*(UL/UT)**2
          ELSE
            SQQ=0.
          ENDIF

          AMASS=0._q
          hills_accel_L=0._q;penalty_accel_L=0._q;hspring_accel_L=0._q
          hills_accel=0._q;penalty_accel=0._q;hspring_accel=0._q
          fbias_accel=0._q;fbias_accel_L=0._q

          IF (DYN%INIT==0) THEN
            DYN%SNOSE=0._q
          ENDIF

          !c store the initial velocity
          VEL_tmp=DYN%VEL 

          IF (iconst5>0) THEN 
            !c update history of fictioous particles and of collerctive variables
            j=0
            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==5) THEN
                j=j+1
                hills%position(j)=ICOORDINATES%COORDSTRUCT(i)%VALUE
              ENDIF
            ENDDO
          
            !c calculate additional potentaials (bias, and penalty)
            IF (counter .GE. hills%bin) THEN
              CALL hills_bias_potential(hills,penalty,ICOORDINATES,DYN,counter,g_io,IO)
            ENDIF 
            CALL penalty_potential(penalty,hills,ICOORDINATES)

            !c compute contributions from bias potentials
            CALL hills_bias_direct(T_INFO,LATT_CUR,DYN,hills%force,ICOORDINATES,DYN%POSIOC,hills_accel,hills_accel_L,3*T_INFO%NIONS,0._q)
            CALL hills_bias_direct(T_INFO,LATT_CUR,DYN,penalty%force,ICOORDINATES,DYN%POSIOC,penalty_accel,penalty_accel_L,3*T_INFO%NIONS,0._q)
            !hills%force=hills%force+penalty%force
            !CALL hills_bias_direct(T_INFO,LATT_CUR,DYN,hills%force,hspring_force,ICOORDINATES,DYN%POSIOC,hills_accel,hills_accel_L,3*T_INFO%NIONS,0._q)
          ENDIF

          !c calculate forces due to Fermi-shaped bias potential
          IF (iconst4>0) THEN
            fbias_force=0.
            CALL fermi_bias_force(iconst4,fbias_r0,fbias_d,fbias_a,fbias_force,ICOORDINATES)
            CALL harmonic_spring_accel(4,T_INFO,LATT_CUR,DYN,fbias_force,ICOORDINATES,DYN%POSIOC,fbias_accel,fbias_accel_L,3*T_INFO%NIONS,0._q)
          ENDIF

          !c calculate forces due to harmonic springs
          IF (iconst8>0) THEN
            CALL harmonic_spring_force(iconst8,hspring_R0,hspring_K,hspring_force,ICOORDINATES)
            CALL harmonic_spring_accel(8,T_INFO,LATT_CUR,DYN,hspring_force,ICOORDINATES,DYN%POSIOC,hspring_accel,hspring_accel_L,3*T_INFO%NIONS,0._q)
          ENDIF

          IF ((iconst4>0) .OR. (iconst8>0)) hspring_accel=hspring_accel+fbias_accel


          IF (SQQ>0) THEN
            !c real variables +  nose-hoover thermostat
            IF (DYN%INIT==1) DYN%VEL=(1._q-DYN%SNOSE(2)/2._q)*DYN%VEL+DYN%D2C-hills_accel-penalty_accel
            EPS=0.5_q*(DYN%SNOSE(2)**2)*SQQ
            ES =NDEGREES_OF_FREEDOM*BOLKEV*DYN%TEMP*(DYN%SNOSE(1))
            DYN%VEL=(DYN%VEL+DYN%D2C-hills_accel-penalty_accel-hspring_accel)/(1._q+DYN%SNOSE(2)/2._q)       
          ELSE 
            !c real variables +  nose-hoover thermostat
            !IF (DYN%INIT==1) DYN%VEL=DYN%VEL+DYN%D2C-hills_accel-penalty_accel
            EPS=0.
            ES =0.
            DYN%VEL=DYN%VEL+2*(DYN%D2C-hills_accel-penalty_accel-hspring_accel)
          ENDIF 
          
          DYN%POSION=DYN%POSIOC+DYN%VEL
          ECONST=0._q
          IF (ICOORDINATES%NUMINTERNALS>0) THEN
            C_FORCE=0._q
            CALL Lproj_out_step(T_INFO,INFO,DYN,LATT_CUR,IO,C_FORCE,C_FORCE_L,ICOORDINATES,&
            &  g_io,3*T_INFO%NIONS,AMASS,LATT_CUR%A,ECONST,iconst0,iconst2)
            DYN%POSION=DYN%POSION+C_FORCE
            DYN%VEL=DYN%VEL+C_FORCE
          ENDIF

          CALL KINETIC_E(EKIN,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,&
                       & DYN%POTIM,LATT_CUR%A,DYN%VEL)

          IF (SQQ>0) THEN
            DYN%SNOSE(3)=DYN%SNOSE(2)
            DYN%SNOSE(2)=DYN%SNOSE(2)+(2*EKIN-NDEGREES_OF_FREEDOM*BOLKEV*DYN%TEMP)/SQQ
            DYN%SNOSE(4)=DYN%SNOSE(1)
            DYN%SNOSE(1)=DYN%SNOSE(1)+(DYN%SNOSE(2)+DYN%SNOSE(3))/2._q
          ENDIF

          !c estimate velocity for the step t 
          VEL_tmp=0.5*(DYN%VEL+VEL_tmp)
          CALL KINETIC_E(EKIN,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS,&
                       & DYN%POTIM,LATT_CUR%A,VEL_tmp)
          TEIN= 2*EKIN/BOLKEV/(NDEGREES_OF_FREEDOM)


          IF (IO%IU6>0) CALL ENERGY_OUT(g_io,DYN,T_INFO,EKIN,EPOT,ECONST,0._q,EPS,ES,counter)

          !c output for metadynamics
          IF (iconst5>0) THEN
            IF (IO%IU6>0) THEN
              j=0
              write(g_io%REPORT,ADVANCE='YES',FMT='(/,2X,A13)') '>Metadynamics'
              DO i=1,ICOORDINATES%NUMINTERNALS
                IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==5) THEN
                  j=j+1
                  write(g_io%REPORT,ADVANCE='YES',FMT='(3X,A6,X,F9.5)') &
                  & 'fic_p>',ICOORDINATES%COORDSTRUCT(i)%VALUE
                  !write(*,*) 'test_val',ICOORDINATES%COORDSTRUCT(i)%VALUE,hills%position(j)
                ENDIF
              ENDDO     
            ENDIF
          ENDIF

          IF (IO%IU6>0) THEN
            CALL TEMPERATURE_OUT(g_io,DYN,TEIN)
          ENDIF


          DISMAX=0
          DO I=1,T_INFO%NIONS
             DISMAX=MAX(DISMAX, &
             &(DYN%VEL(1,I))**2+ &
             &(DYN%VEL(2,I))**2+ &
             &(DYN%VEL(3,I))**2)
          ENDDO
          DISMAX=SQRT(DISMAX)
        END SUBROUTINE HILLS_METHOD_DIRECT_nose

        SUBROUTINE tis_stop(DYN,T_INFO,INFO,LATT_CUR,IO,g_io,ICOORDINATES,LMIN,LMAX,valA,valB)
          !c monitor coordinates with STATUS=7
          !c terminate if the value doesn't fall into interval (VALUE_MIN,VALUE_MAX)
          !c this is useful mainly for the transition interface sampling calculations:
          USE reader_tags
          TYPE(coordstructure) :: ICOORDINATES
          TYPE(type_info) :: T_INFO
          TYPE (info_struct) INFO
          TYPE(dynamics) :: DYN      
          TYPE (in_struct) :: IO
          TYPE(latt) :: LATT_CUR   
          TYPE(gadget_io) :: g_io
          REAL(q) :: valA(ICOORDINATES%iconst7),valB(ICOORDINATES%iconst7)
          INTEGER,SAVE :: counter=0
          INTEGER :: IDUM,IERR,N,i,ii
          LOGICAL :: LMIN
          LOGICAL :: LMAX
          LOGICAL :: LOPEN,LDUM
          REAL(q) :: RDUM
          COMPLEX(q) :: CDUM
          CHARACTER*1 :: CHARAC

          !c dont continue if no geometric parameters with STATUS==7
          IF (ICOORDINATES%iconst7==0) THEN
            RETURN
          ENDIF

          IF (IO%IU6>0) write(g_io%REPORT,FMT='(/,2X,A12)') '>Monit_coord'
          DO i=1,ICOORDINATES%NUMINTERNALS
            IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==7) THEN
              IF (IO%IU6>0) write(g_io%REPORT,FMT='(3X,A3,X,A2,X,F15.5)') &
                'mc>',ICOORDINATES%COORDSTRUCT(i)%TAG,ICOORDINATES%COORDSTRUCT(i)%VALUE
            ENDIF
          ENDDO
   
          IF (LMIN) THEN
            INFO%LSTOP=.FALSE.
            ii=0
            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==7) THEN
                ii=ii+1
                IF (ICOORDINATES%COORDSTRUCT(i)%VALUE .LE. valA(ii)) THEN
                  INFO%LSTOP=.TRUE.
                  EXIT
                ENDIF
              ENDIF
            ENDDO
            IF (INFO%LSTOP) THEN
              CALL vtutor%warning("step_tb: lower limit for STATUS=7 coord. reached")
              !IF (IO%IU0>=0) write(*,*) 'REACHED LOWER LIMIT FOR STATUS=7 COORD.'
              IF (INFO%LSTOP) RETURN
            ENDIF
          ENDIF

          IF (LMAX) THEN
            INFO%LSTOP=.FALSE.
            ii=0
            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==7) THEN
                ii=ii+1
                IF (ICOORDINATES%COORDSTRUCT(i)%VALUE .GE. valB(ii)) THEN
                  INFO%LSTOP=.TRUE.
                  EXIT
                ENDIF
              ENDIF
            ENDDO
            IF (INFO%LSTOP) THEN
              CALL vtutor%warning("step_tb: upper limit for STATUS=7 coord. reached")
              !IF (IO%IU0>=0) write(*,*) 'REACHED UPPER LIMIT FOR STATUS=7 COORD.'
              IF (INFO%LSTOP) RETURN
            ENDIF
          ENDIF
        END SUBROUTINE tis_stop

        SUBROUTINE Lproj_out_step(T_INFO,INFO,DYN,LATT_CUR, &
                IO,C_FORCE,C_FORCE_L,ICOORDINATES,g_io,IBDIM,AMASS,NewA,ECONST,iconst0,iconst2)
          !c this is basically the SHAKE algorithm...
          TYPE(type_info) :: T_INFO
          TYPE(dynamics) :: DYN      
          TYPE (in_struct) :: IO
          TYPE (info_struct) INFO
          TYPE(latt) :: LATT_CUR     
          TYPE(coordstructure) :: ICOORDINATES     !c working and reference coords.
          TYPE(gadget_io) :: g_io
          REAL(q) ::xpos(3*T_INFO%NIONS)  ! positions in row format
          INTEGER :: i,j,k  ,ii ,jj,idummy                 ! indeces
          INTEGER :: ninfo,IBDIM
          INTEGER,SAVE :: counter=0
          REAL(q) :: BMAT_0(ICOORDINATES%NUMINTERNALS,IBDIM)
          REAL(q) :: BMAT(ICOORDINATES%NUMINTERNALS,IBDIM)
          REAL(q) :: GMAT(ICOORDINATES%NUMINTERNALS) !c B matrix for working and reference coords.
          REAL(q) :: PRIMS1(ICOORDINATES%NUMINTERNALS),PRIMS2(ICOORDINATES%NUMINTERNALS)
          REAL(q) :: gamm(ICOORDINATES%NUMINTERNALS),gamm0(ICOORDINATES%NUMINTERNALS)
          REAL(q) :: ZET(ICOORDINATES%NUMINTERNALS),ZET_DET
          REAL(q),ALLOCATABLE,SAVE :: MASSES(:)
          REAL(q) :: C_FORCE(3,T_INFO%NIONS),C_FORCE_(3,T_INFO%NIONS)
          REAL(q) :: POSION_(3,T_INFO%NIONS),POSIOC_(3,T_INFO%NIONS)
          REAL(q) :: err,err1
          REAL(q) :: UL,UT,FACT,ECONST
          INTEGER IDUM,IERR,N
          INTEGER :: iconst0,iconst2
          LOGICAL :: LOPEN,LDUM
          REAL(q) :: RDUM
          COMPLEX(q) :: CDUM
          CHARACTER*1 :: CHARAC
          !REAL(q),SAVE,ALLOCATABLE :: INCREM(:)
          !INTEGER,SAVE :: EQUI_REGIME=0
          !LOGICAL,SAVE :: LBLUEOUT=.FALSE.
          !REAL(q),SAVE :: BMTOL,BMTOLSOFT,BMSCA 
          !INTEGER,SAVE :: BMITER
          REAL(q) :: AMASS,NewA(3,3),NewA_(3,3),NewB(3,3),C_FORCE_L(3,3),C_FORCE_L_(3,3),xposL(9)
          REAL(q) :: BMAT_square(ICOORDINATES%NUMINTERNALS,ICOORDINATES%NUMINTERNALS),BMAT_inv(IBDIM,ICOORDINATES%NUMINTERNALS)

          !c don't continue if no geometric constraints are defined
          IF (iconst0==0 .AND. iconst2==0) THEN
            RETURN
          ENDIF

          counter=counter+1
          UL=1E-10_q
          UT=DYN%POTIM*1E-15_q
          C_FORCE=0._q
          C_FORCE_L=0._q

          POSIOC_=DYN%POSIOC
          POSION_=DYN%POSION
          NewA_=NewA
          NewB=TRANSPOSE(NewA)
          CALL SVDINVERSE(NewB,3,idummy)
          
          IF (counter==1) THEN
            !ICOORDINATES%COORDSTRUCT(:)%DVALUE=ICOORDINATES%COORDSTRUCT(:)%VALUE  !xx
            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0) THEN
                ICOORDINATES%COORDSTRUCT(i)%DVALUE=ICOORDINATES%COORDSTRUCT(i)%VALUE  !xx
              ENDIF
            ENDDO
            !c masses in amu
            ALLOCATE(MASSES(IBDIM))
            MASSES=0._q
            DO i=1,T_INFO%NIONS
              MASSES(3*i-2:3*i)=(T_INFO%POMASS(T_INFO%ITYP(i)))
            ENDDO
            IF (IBDIM==3*T_INFO%NIONS+9) THEN
              DO i=1,3
                DO j=1,3
                  MASSES(3*T_INFO%NIONS+3*(i-1)+j)=AMASS
                ENDDO
              ENDDO
            ENDIF

            !ALLOCATE(INCREM(iconst0))
            !INCREM=0._q
            !c read in some method-specific input parameters
            !CALL READER_ICONST0(IO,g_io,iconst0,BMTOL,BMTOLSOFT,BMSCA,BMITER,LBLUEOUT,EQUI_REGIME,INCREM)
          END IF

          !c dont fix coordinates in the equilibration period
          IF (EQUI_REGIME_GL .GE. counter) RETURN

          !c everything for old coordinates
          !c read constraints
          j=0
          DO i=1,ICOORDINATES%NUMINTERNALS
            IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0) THEN
              j=j+1
              PRIMS1(i)=ICOORDINATES%COORDSTRUCT(i)%DVALUE+counter*INCREM_GL(j) 
            ELSE IF  (ICOORDINATES%COORDSTRUCT(i)%STATUS==2) THEN
              j=j+1
              !PRIMS1(i)=ICOORDINATES%COORDSTRUCT(i)%DVALUE
              PRIMS1(i)=ICOORDINATES%COORDSTRUCT(i)%VALUE
            ENDIF
          ENDDO

          BMAT_0=0._q
          CALL BMATRIX(T_INFO,DYN%POSIOC,DYN%POSIOC,TRANSPOSE(LATT_CUR%A),ICOORDINATES,BMAT_0,.TRUE.)
          !!CALL BMATRIX(T_INFO,DYN%POSION,DYN%POSION,TRANSPOSE(LATT_CUR%A),ICOORDINATES,BMAT_0,.FALSE.)
          !c upper triangle inly
          IF (IBDIM==3*T_INFO%NIONS+9) THEN
            BMAT_0(:,3*T_INFO%NIONS+4)=0._q
            BMAT_0(:,3*T_INFO%NIONS+7)=0._q
            BMAT_0(:,3*T_INFO%NIONS+8)=0._q
          ENDIF

          !CALL ICONST_BMAT(BMAT_0,ICOORDINATES,T_INFO)
          BMAT=0._q

          err=0._q
          gamm=0._q
          gamm0=0._q
          ii=0

          DO i=1,ICOORDINATES%NUMINTERNALS      
            IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0 .OR. ICOORDINATES%COORDSTRUCT(i)%STATUS==2) THEN
              err=err+1
            ENDIF
          ENDDO
        
          DO
            !c exit if no constraints are defined
            IF (ii==0 .AND. err<BMTOL_GL) THEN
              EXIT
            ENDIF
            ii=ii+1

            !c terminate with error message if SHAKE doesn't converge in BMITER steps
            IF (ii>BMITER_GL) THEN  
              IF (IO%IU6>0) write(*,*) 'Error: SHAKE algorithm did not converge! err=',err
              IF (abs(err).LE. BMTOLSOFT_GL) THEN
                IF (IO%IU6>0) write(*,*) "I'll try to recover this problem in the next step"
                EXIT
              ELSE
                CALL vtutor%error("Error too large, I have to terminate this calculation!")
              ENDIF
              !RETURN
            END IF
           
            !c values for constraint parameters in the current step 
            CALL DEAL_XYZ(T_INFO,DYN%POSION,DYN%POSION,NewA,ICOORDINATES)
            PRIMS2=0._q
            PRIMS2=ICOORDINATES%COORDSTRUCT(:)%VALUE
          
            CALL DECYCLE(PRIMS1,PRIMS2,ICOORDINATES)
            ICOORDINATES%COORDSTRUCT(:)%VALUE=PRIMS2(:)
            IF (err<BMTOL_GL) EXIT
            BMAT=0._q
       
            !c Jacobi matrix
            CALL BMATRIX(T_INFO,DYN%POSION,DYN%POSION,TRANSPOSE(NewA),ICOORDINATES,BMAT,.TRUE.)
            !c upper triangle only
            IF (IBDIM==3*T_INFO%NIONS+9) THEN
              BMAT(:,3*T_INFO%NIONS+4)=0._q
              BMAT(:,3*T_INFO%NIONS+7)=0._q
              BMAT(:,3*T_INFO%NIONS+8)=0._q
            ENDIF

            jj=0
            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0 .OR. ICOORDINATES%COORDSTRUCT(i)%STATUS==2) THEN
                jj=jj+1
                IF (jj==1) err=ABS(PRIMS2(i)-PRIMS1(i))
                err1=ABS(PRIMS2(i)-PRIMS1(i))
                IF (err1 .GT. err) err=err1
                gamm(i)=(PRIMS2(i)-PRIMS1(i))/SUM(BMAT_0(i,:)/MASSES(:)*BMAT(i,:))/BMSCA_GL     !!!
                gamm0(i)=gamm0(i)+(PRIMS2(i)-PRIMS1(i))/SUM(BMAT_0(i,:)/MASSES(:)*BMAT(i,:))/BMSCA_GL   !!!

                !c lattice dynamics
                IF (IBDIM==3*T_INFO%NIONS+9) THEN
                  C_FORCE_L_=0._q
                  xposL=-gamm(i)*BMAT_0(i,3*T_INFO%NIONS+1:3*T_INFO%NIONS+9)/MASSES(3*T_INFO%NIONS+1:3*T_INFO%NIONS+9)
                  CALL ONETOTHREE(3,C_FORCE_L_,xposL) 
                  C_FORCE_L_=TRANSPOSE(C_FORCE_L_) !!
                  NewA(:,:)=NewA(:,:)+C_FORCE_L_
                  NewB=TRANSPOSE(NewA)
                  CALL SVDINVERSE(NewB,3,idummy)
                  IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0) THEN
                    C_FORCE_L=C_FORCE_L+C_FORCE_L_
                  ELSE IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==2) THEN
                    C_FORCE_L=C_FORCE_L+2*C_FORCE_L_
                  ENDIF
                ENDIF

                xpos=-gamm(i)*BMAT_0(i,1:3*T_INFO%NIONS)/MASSES(1:3*T_INFO%NIONS)
                !!!xpos=-gamm(i)*BMAT(i,1:3*T_INFO%NIONS)/MASSES(1:3*T_INFO%NIONS)
                C_FORCE_=0._q
                CALL ONETOTHREE(T_INFO%NIONS,C_FORCE_,xpos)
                !C_FORCE_=MATMUL(TRANSPOSE(LATT_CUR%B),C_FORCE_)
                C_FORCE_=MATMUL(TRANSPOSE(NewB),C_FORCE_)
                !beg 
                DO j=1,T_INFO%NIONS
                  IF (.NOT. T_INFO%LSFOR(1,j)) C_FORCE_(1,j)=0._q
                  IF (.NOT. T_INFO%LSFOR(2,j)) C_FORCE_(2,j)=0._q
                  IF (.NOT. T_INFO%LSFOR(3,j)) C_FORCE_(3,j)=0._q
                ENDDO
                !end
                DYN%POSION(:,:)=DYN%POSION(:,:)+C_FORCE_(:,:)
                IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0) THEN
                   C_FORCE=C_FORCE+C_FORCE_
                ELSE IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==2) THEN
                   C_FORCE=C_FORCE+2*C_FORCE_
                ENDIF

              END IF
            ENDDO                     
            !!BMAT_0=BMAT
            !IF (IO%IU0>=0) write(IO%IU0,*) 'err',err
          ENDDO

          gamm=-gamm0

          !CALL XML_INCAR('shakeiter','I',ii,RDUM,CDUM,LDUM,CHARAC,1)
          IF (IO%IU6>0) write(g_io%REPORT,FMT='(/,2X,A19,I4,A6)') 'SHAKE converged in ',ii,' steps'

          !CALL XML_TAG("separator","gconstraints")
          IF (IO%IU6>0) write(g_io%REPORT,FMT='(/,2X,A12)') '>Const_coord'
          DO i=1,ICOORDINATES%NUMINTERNALS
            IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0 .OR. ICOORDINATES%COORDSTRUCT(i)%STATUS==2) THEN
              
              err=PRIMS2(i)-PRIMS1(i)
              IF (IO%IU6>0) write(g_io%REPORT,FMT='(3X,A3,X,A2,2(X,F9.5),X,E13.6E2)') 'cc>',ICOORDINATES%COORDSTRUCT(i)%TAG,&
                  PRIMS1(i),PRIMS2(i),err
            ENDIF
          ENDDO
          !CALL XML_CLOSE_TAG

!!!!!!!!!!!!!!!!!!!!!!!!!
          !c write down the componenets needed to compute the free energy gradients
          !c this can be time-consuming if too many constraints are defined
          IF (LBLUEOUT_GL) THEN                
            CALL make_Z(T_INFO,IBDIM,DYN%POSIOC,LATT_CUR%A,ICOORDINATES,MASSES,ZET) 
            ZET_DET=z_det(ZET,ICOORDINATES,IO) 
            

            CALL make_G(T_INFO,DYN,LATT_CUR,ICOORDINATES,MASSES,IBDIM,ZET,ZET_DET,GMAT,IO)
            GMAT=GMAT*BOLKEV*DYN%TEMP      
         
            FACT=(DYN%POTIM**2)*EVTOJ/(AMTOKG) *1E-10_q
            gamm=gamm/FACT

            LDUM=.TRUE.
            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0 .AND. IO%IU6>0) THEN 
                IF (LDUM) THEN
                  write(g_io%REPORT,FMT='(/,2X,A10)') '>Blue_moon'
                  write(g_io%REPORT,FMT='(8X,A6,8X,A10,4X,A3,11X,A23)') 'lambda','|z|^(-1/2)','GkT','|z|^(-1/2)*(lambda+GkT)'
                  LDUM=.FALSE.
                ENDIF
                write(g_io%REPORT,FMT='(3X,A4,(X,E13.6E2),(X,E13.6E2),(X,E13.6E2),X,E13.6E2)') &
                    'b_m>', gamm(i),1._q/SQRT(ZET_DET),GMAT(i),1._q/SQRT(ZET_DET)*(gamm(i)+GMAT(i))
              ENDIF
            ENDDO
          ENDIF

          !calculate the energy contributions due to constraints
          ECONST=0._q
          CALL DEAL_XYZ(T_INFO,POSION_,DYN%POSION,LATT_CUR%A,ICOORDINATES)
          PRIMS2=0._q
          PRIMS2=ICOORDINATES%COORDSTRUCT(:)%VALUE
          CALL DECYCLE(PRIMS1,PRIMS2,ICOORDINATES)
          DO i=1,ICOORDINATES%NUMINTERNALS
            IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0) THEN
              ECONST=ECONST+gamm(i)*(PRIMS2(i)-PRIMS1(i))
            ENDIF
          ENDDO

          DYN%POSIOC=POSIOC_
          DYN%POSION=POSION_
          NewA=NewA_
        END SUBROUTINE Lproj_out_step

        SUBROUTINE Rattle(counter,T_INFO,INFO,DYN,LATT_CUR, &
                IO,C_FORCE,C_FORCE_L,V_FORCE,V_FORCE_L,ICOORDINATES,g_io,IBDIM,AMASS,NewA,ECONST,iconst0,iconst2,IMODE)
          TYPE(type_info) :: T_INFO
          TYPE(dynamics) :: DYN      
          TYPE (in_struct) :: IO
          TYPE (info_struct) INFO
          TYPE(latt) :: LATT_CUR     
          TYPE(coordstructure) :: ICOORDINATES     !c working and reference coords.
          TYPE(gadget_io) :: g_io
          REAL(q) ::xpos(3*T_INFO%NIONS)  ! positions in row format
          INTEGER :: i,j,k  ,ii ,jj,idummy                 ! indeces
          INTEGER :: ninfo,IBDIM
          INTEGER :: IMODE
          INTEGER :: counter
          !INTEGER,SAVE :: counter=0
          REAL(q) :: BMAT_0(ICOORDINATES%NUMINTERNALS,IBDIM)
          REAL(q) :: BMAT(ICOORDINATES%NUMINTERNALS,IBDIM)
          REAL(q) :: GMAT(ICOORDINATES%NUMINTERNALS) !c B matrix for working and reference coords.
          REAL(q) :: PRIMS1(ICOORDINATES%NUMINTERNALS),PRIMS2(ICOORDINATES%NUMINTERNALS)
          REAL(q) :: gamm(ICOORDINATES%NUMINTERNALS),gamm0(ICOORDINATES%NUMINTERNALS),theta(ICOORDINATES%NUMINTERNALS)
          REAL(q) :: ZET(ICOORDINATES%NUMINTERNALS),ZET_DET
          REAL(q),ALLOCATABLE,SAVE :: MASSES(:)
          REAL(q) :: C_FORCE(3,T_INFO%NIONS),C_FORCE_(3,T_INFO%NIONS)
          REAL(q) :: V_FORCE(3,T_INFO%NIONS),V_FORCE_L(3,3)
          REAL(q) :: POSION_(3,T_INFO%NIONS),POSIOC_(3,T_INFO%NIONS), VEL_(3,T_INFO%NIONS)
          REAL(q) :: err,err1
          REAL(q) :: UL,UT,FACT,ECONST
          INTEGER IDUM,IERR,N
          INTEGER :: iconst0,iconst2
          LOGICAL :: LOPEN,LDUM
          LOGICAL,SAVE :: LFIRST=.TRUE.
          REAL(q) :: RDUM
          COMPLEX(q) :: CDUM
          CHARACTER*1 :: CHARAC
          !REAL(q),SAVE,ALLOCATABLE :: INCREM(:)
          !INTEGER,SAVE :: EQUI_REGIME=0
          !LOGICAL,SAVE :: LBLUEOUT=.FALSE.
          !REAL(q),SAVE :: BMTOL,BMTOLSOFT,BMSCA 
          !INTEGER,SAVE :: BMITER
          REAL(q) :: AMASS,NewA(3,3),NewA_(3,3),NewB(3,3),C_FORCE_L(3,3),C_FORCE_L_(3,3),xposL(9)
          REAL(q) :: AVEL_(3,3)
          REAL(q) :: BMAT_square(ICOORDINATES%NUMINTERNALS,ICOORDINATES%NUMINTERNALS),BMAT_inv(IBDIM,ICOORDINATES%NUMINTERNALS)

          V_FORCE_L=0._q
          C_FORCE_L=0._q
          V_FORCE=0._q
          C_FORCE=0._q

          !c don't continue if no geometric constraints are defined
          IF (iconst0==0 .AND. iconst2==0) THEN
            RETURN
          ENDIF

          !counter=counter+1

          UL=1E-10_q
          UT=DYN%POTIM*1E-15_q
          C_FORCE=0._q
          C_FORCE_L=0._q

          POSIOC_=DYN%POSIOC
          POSION_=DYN%POSION
          NewA_=NewA
          NewB=TRANSPOSE(NewA)
          CALL SVDINVERSE(NewB,3,idummy)
          
          IF (LFIRST) THEN
            LFIRST=.FALSE.
            !ICOORDINATES%COORDSTRUCT(:)%DVALUE=ICOORDINATES%COORDSTRUCT(:)%VALUE  !xx
            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0) THEN
                ICOORDINATES%COORDSTRUCT(i)%DVALUE=ICOORDINATES%COORDSTRUCT(i)%VALUE 
              ENDIF
            ENDDO
            
            !c masses in amu
            ALLOCATE(MASSES(IBDIM))
            MASSES=0._q
            DO i=1,T_INFO%NIONS
              MASSES(3*i-2:3*i)=(T_INFO%POMASS(T_INFO%ITYP(i)))
            ENDDO
            IF (IBDIM==3*T_INFO%NIONS+9) THEN
              DO i=1,3
                DO j=1,3
                  MASSES(3*T_INFO%NIONS+3*(i-1)+j)=AMASS
                ENDDO
              ENDDO
            ENDIF

            !ALLOCATE(INCREM(iconst0))
            !INCREM=0._q
            !c read in some method-specific input parameters
            !CALL READER_ICONST0(IO,g_io,iconst0,BMTOL,BMTOLSOFT,BMSCA,BMITER,LBLUEOUT,EQUI_REGIME,INCREM)
          END IF

          !c dont fix coordinates in the equilibration period
          IF (EQUI_REGIME_GL .GE. counter) RETURN

          !c everything for old coordinates
          !c read constraints
          j=0
          DO i=1,ICOORDINATES%NUMINTERNALS
            IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0) THEN
              j=j+1
              PRIMS1(i)=ICOORDINATES%COORDSTRUCT(i)%DVALUE+counter*INCREM_GL(j) 
            ELSE IF  (ICOORDINATES%COORDSTRUCT(i)%STATUS==2) THEN
              j=j+1
              !PRIMS1(i)=ICOORDINATES%COORDSTRUCT(i)%DVALUE
              PRIMS1(i)=ICOORDINATES%COORDSTRUCT(i)%VALUE
            ENDIF
          ENDDO

          SELECT CASE(IMODE)
          CASE(0)
            BMAT_0=0._q
            CALL BMATRIX(T_INFO,DYN%POSIOC,DYN%POSIOC,TRANSPOSE(LATT_CUR%A),ICOORDINATES,BMAT_0,.TRUE.)
 
            BMAT=0._q

            err=0._q
            gamm=0._q
            gamm0=0._q
            ii=0

            DO i=1,ICOORDINATES%NUMINTERNALS      
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0 .OR. ICOORDINATES%COORDSTRUCT(i)%STATUS==2) THEN
                err=err+1
              ENDIF
            ENDDO
        
            DO
              !c exit if no constraints are defined
              IF (ii==0 .AND. err<BMTOL_GL) THEN
                EXIT
              ENDIF
              ii=ii+1

              !c terminate with error message if SHAKE doesn't converge in BMITER steps
              IF (ii>BMITER_GL) THEN  
                !INFO%LSTOP=.TRUE.
                IF (IO%IU6>0) write(*,*) 'Error: RATTLE algorithm did not converge! err=',err
                IF (abs(err).LE. BMTOLSOFT_GL) THEN
                  IF (IO%IU6>0) write(*,*) "I'll try to recover this problem in the next step"
                  EXIT
                ELSE
                  CALL vtutor%error("Error too large, I have to terminate this calculation!")
                ENDIF
                !RETURN
              END IF

           
              !c values for constraint parameters in the current step 
              CALL DEAL_XYZ(T_INFO,DYN%POSION,DYN%POSION,NewA,ICOORDINATES)
              PRIMS2=0._q
              PRIMS2=ICOORDINATES%COORDSTRUCT(:)%VALUE
          
              CALL DECYCLE(PRIMS1,PRIMS2,ICOORDINATES)
              ICOORDINATES%COORDSTRUCT(:)%VALUE=PRIMS2(:)
              IF (err<BMTOL_GL) EXIT
              BMAT=0._q
       
              !c Jacobi matrix
              CALL BMATRIX(T_INFO,DYN%POSION,DYN%POSION,TRANSPOSE(NewA),ICOORDINATES,BMAT,.TRUE.)
 
              jj=0
              DO i=1,ICOORDINATES%NUMINTERNALS
                IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0 .OR. ICOORDINATES%COORDSTRUCT(i)%STATUS==2) THEN
                  jj=jj+1
                  IF (jj==1) err=ABS(PRIMS2(i)-PRIMS1(i))
                  err1=ABS(PRIMS2(i)-PRIMS1(i))
                  IF (err1 .GT. err) err=err1
                  gamm(i)=(PRIMS2(i)-PRIMS1(i))/SUM(BMAT_0(i,:)/MASSES(:)*BMAT(i,:))/BMSCA_GL     !!!
                  gamm0(i)=gamm0(i)+(PRIMS2(i)-PRIMS1(i))/SUM(BMAT_0(i,:)/MASSES(:)*BMAT(i,:))/BMSCA_GL   !!!

                  !c lattice dynamics
                  IF (IBDIM==3*T_INFO%NIONS+9) THEN
                    C_FORCE_L_=0._q
                    xposL=-gamm(i)*BMAT_0(i,3*T_INFO%NIONS+1:3*T_INFO%NIONS+9)/MASSES(3*T_INFO%NIONS+1:3*T_INFO%NIONS+9)
                    CALL ONETOTHREE(3,C_FORCE_L_,xposL) 
                    !!C_FORCE_L_=TRANSPOSE(C_FORCE_L_) !!
                    NewA(:,:)=NewA(:,:)+C_FORCE_L_
                    NewB=TRANSPOSE(NewA)
                    CALL SVDINVERSE(NewB,3,idummy)
                    IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0) THEN
                      C_FORCE_L=C_FORCE_L+C_FORCE_L_
!                     ELSE IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==2) THEN
!                       C_FORCE_L=C_FORCE_L+2*C_FORCE_L_
                    ENDIF
                  ENDIF

                  xpos=-gamm(i)*BMAT_0(i,1:3*T_INFO%NIONS)/MASSES(1:3*T_INFO%NIONS)
                  C_FORCE_=0._q
                  CALL ONETOTHREE(T_INFO%NIONS,C_FORCE_,xpos)
                  !C_FORCE_=MATMUL(TRANSPOSE(LATT_CUR%B),C_FORCE_)
                  C_FORCE_=MATMUL(TRANSPOSE(NewB),C_FORCE_)

                  !c remove force for frozen fractional coords 
                  DO j=1,T_INFO%NIONS
                    IF (.NOT. T_INFO%LSFOR(1,j)) C_FORCE_(1,j)=0._q
                    IF (.NOT. T_INFO%LSFOR(2,j)) C_FORCE_(2,j)=0._q
                    IF (.NOT. T_INFO%LSFOR(3,j)) C_FORCE_(3,j)=0._q
                  ENDDO

                  DYN%POSION(:,:)=DYN%POSION(:,:)+C_FORCE_(:,:)
                  IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0) THEN
                     C_FORCE=C_FORCE+C_FORCE_
!                   ELSE IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==2) THEN
!                      C_FORCE=C_FORCE+2*C_FORCE_
                  ENDIF

                END IF
              ENDDO                     
            !!BMAT_0=BMAT
            !IF (IO%IU0>=0) write(IO%IU0,*) 'err',err
            ENDDO

            gamm=-gamm0

            !CALL XML_INCAR('shakeiter','I',ii,RDUM,CDUM,LDUM,CHARAC,1)
            IF (IO%IU6>0) write(g_io%REPORT,FMT='(/,2X,A24,I4,A6)') 'RATTLE_pos converged in ',ii,' steps'

            !CALL XML_TAG("separator","gconstraints")
            IF (IO%IU6>0) write(g_io%REPORT,FMT='(/,2X,A12)') '>Const_coord'
            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0 .OR. ICOORDINATES%COORDSTRUCT(i)%STATUS==2) THEN              
                err=PRIMS2(i)-PRIMS1(i)
                IF (IO%IU6>0) write(g_io%REPORT,FMT='(3X,A3,X,A2,2(X,F9.5),X,E13.6E2)') 'cc>',ICOORDINATES%COORDSTRUCT(i)%TAG,&
                &   PRIMS1(i),PRIMS2(i),err
              ENDIF
            ENDDO
            !CALL XML_CLOSE_TAG

            !calculate the energy contributions due to constraints
            ECONST=0._q
            CALL DEAL_XYZ(T_INFO,POSION_,DYN%POSION,LATT_CUR%A,ICOORDINATES)
            PRIMS2=0._q
            PRIMS2=ICOORDINATES%COORDSTRUCT(:)%VALUE
            CALL DECYCLE(PRIMS1,PRIMS2,ICOORDINATES)
            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0) THEN
                ECONST=ECONST+gamm(i)*(PRIMS2(i)-PRIMS1(i))
              ENDIF
            ENDDO

            DYN%POSIOC=POSIOC_
            DYN%POSION=POSION_
            NewA=NewA_


            !c write down the componenets needed to compute the free energy gradients
            !c this can be time-consuming if too many constraints are defined
            IF (LBLUEOUT_GL) THEN                
              !c presently available only for NVT
              IF (IBDIM==3*T_INFO%NIONS) THEN
                CALL make_Z(T_INFO,IBDIM,DYN%POSIOC,LATT_CUR%A,ICOORDINATES,MASSES,ZET) 
                ZET_DET=z_det(ZET,ICOORDINATES,IO) 

                CALL make_G(T_INFO,DYN,LATT_CUR,ICOORDINATES,MASSES,IBDIM,ZET,ZET_DET,GMAT,IO)
                GMAT=GMAT*BOLKEV*DYN%TEMP      
         
                FACT=(DYN%POTIM**2)*EVTOJ/(AMTOKG) *1E-10_q
                gamm=gamm/FACT

                LDUM=.TRUE.
                DO i=1,ICOORDINATES%NUMINTERNALS
                  IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0 .AND. IO%IU6>0) THEN 
                    IF (LDUM) THEN
                      write(g_io%REPORT,FMT='(/,2X,A10)') '>Blue_moon'
                      write(g_io%REPORT,FMT='(8X,A6,8X,A10,4X,A3,11X,A23)') 'lambda','|z|^(-1/2)','GkT','|z|^(-1/2)*(lambda+GkT)'
                      LDUM=.FALSE.
                    ENDIF
                    write(g_io%REPORT,FMT='(3X,A4,(X,E13.6E2),(X,E13.6E2),(X,E13.6E2),X,E13.6E2)') &
                        'b_m>', gamm(i),1._q/SQRT(ZET_DET),GMAT(i),1._q/SQRT(ZET_DET)*(gamm(i)+GMAT(i))
                  ENDIF
                ENDDO
              ENDIF
            ENDIF

          CASE(1)
            BMAT=0._q
       
            !c Jacobi matrix
            CALL BMATRIX(T_INFO,DYN%POSIOC,DYN%POSIOC,TRANSPOSE(NewA),ICOORDINATES,BMAT,.TRUE.)

            CALL ICONST_BMAT(BMAT,ICOORDINATES,T_INFO)

            VEL_=DYN%VEL
            V_FORCE=0._q

            AVEL_=LATT_CUR%Avel
            V_FORCE_L=0._q

            err=0._q

            DO i=1,ICOORDINATES%NUMINTERNALS      
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0 .OR. ICOORDINATES%COORDSTRUCT(i)%STATUS==2) THEN
                err=err+1
              ENDIF
            ENDDO

            !c start iterations to fullfill constraints for velocities
            DO ii=1,BMITER_GL+1

              !c exit if no constraints are defined
              IF (ii==1 .AND. err<BMTOL_GL) THEN
                !V_FORCE=VEL_-DYN%VEL
                EXIT
              ENDIF
           
              !c terminate with error message if SHAKE doesn't converge in BMITER steps
              IF (ii>BMITER_GL) THEN  
                IF (IO%IU6>0) write(*,*) 'Error: RATTLE_vel algorithm did not converge! err=',err
                IF (abs(err).LE. BMTOLSOFT_GL) THEN
                  IF (IO%IU6>0) write(*,*) "I'll try to recover this problem in the next step"
                  EXIT
                ELSE
                  CALL vtutor%error("Error too large, I have to terminate this calculation!")
                ENDIF
                !RETURN
              END IF
           
              IF (err<BMTOL_GL) THEN 
                EXIT
              ENDIF

              PRIMS2=0._q
              VEL_=MATMUL(LATT_CUR%A,VEL_)
              CALL vel_bmat_product(T_INFO%NIONS,ICOORDINATES%NUMINTERNALS,IBDIM,VEL_,AVEL_,BMAT,PRIMS2)
              VEL_=MATMUL(TRANSPOSE(LATT_CUR%B),VEL_)

             ! IF (IO%IU0>=0) write(*,*) "errVEL:",ii,MAXVAL(ABS(PRIMS2))

              jj=0
              DO i=1,ICOORDINATES%NUMINTERNALS
                IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0 .OR. ICOORDINATES%COORDSTRUCT(i)%STATUS==2) THEN
                  jj=jj+1
                  IF (jj==1) err=ABS(PRIMS2(i)-INCREM_GL(jj))
                  err1=ABS(PRIMS2(i)-INCREM_GL(jj))
                  IF (err1 .GT. err) err=err1
                  RDUM=SUM(BMAT(i,:)/MASSES(:)*BMAT(i,:))
                  theta(i)=0.
                  IF (ABS(RDUM) .GT. 1e-6) theta(i)=(PRIMS2(i)-INCREM_GL(jj))/RDUM/BMSCA_GL     !!!
                
                  !c lattice dynamics
                  IF (IBDIM==3*T_INFO%NIONS+9) THEN
                    C_FORCE_L_=0._q
                    xposL=-theta(i)*BMAT(i,3*T_INFO%NIONS+1:3*T_INFO%NIONS+9)/MASSES(3*T_INFO%NIONS+1:3*T_INFO%NIONS+9)
                    CALL ONETOTHREE(3,C_FORCE_L_,xposL) 
                    !!C_FORCE_L_=TRANSPOSE(C_FORCE_L_) !!
                    AVEL_=AVEL_+2*C_FORCE_L_
                  ENDIF
             
                  xpos=-theta(i)*BMAT(i,1:3*T_INFO%NIONS)/MASSES(1:3*T_INFO%NIONS)
                  C_FORCE_=0._q
                  CALL ONETOTHREE(T_INFO%NIONS,C_FORCE_,xpos)

                  !C_FORCE_=MATMUL(TRANSPOSE(LATT_CUR%B),C_FORCE_)
                  C_FORCE_=MATMUL(TRANSPOSE(NewB),C_FORCE_)

                  !c remove force for frozen fractional coords 
                  DO j=1,T_INFO%NIONS
                    IF (.NOT. T_INFO%LSFOR(1,j)) C_FORCE_(1,j)=0._q
                    IF (.NOT. T_INFO%LSFOR(2,j)) C_FORCE_(2,j)=0._q
                    IF (.NOT. T_INFO%LSFOR(3,j)) C_FORCE_(3,j)=0._q
                  ENDDO

                  VEL_=VEL_+2*C_FORCE_(:,:)
               END IF
              ENDDO                     
            ENDDO
            V_FORCE=VEL_-DYN%VEL
            IF (IBDIM==3*T_INFO%NIONS+9) V_FORCE_L=AVEL_-LATT_CUR%AVEL

            !CALL XML_INCAR('shakeiter','I',ii,RDUM,CDUM,LDUM,CHARAC,1)
            IF (IO%IU6>0) write(g_io%REPORT,FMT='(/,2X,A24,I4,A6)') 'RATTLE_vel converged in ',ii,' steps'

            !CALL XML_TAG("separator","gconstraints")
            IF (IO%IU6>0) write(g_io%REPORT,FMT='(/,2X,A10)') '>Const_vel'
            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0 .OR. ICOORDINATES%COORDSTRUCT(i)%STATUS==2) THEN
                err=PRIMS2(i)
                IF (IO%IU6>0) write(g_io%REPORT,FMT='(3X,A3,X,A2,X,E13.6E2)') 'cv>',ICOORDINATES%COORDSTRUCT(i)%TAG,&
                &    err
              ENDIF
            ENDDO
          END SELECT

        END SUBROUTINE Rattle


        SUBROUTINE vel_bmat_product(n,m,ibdim,vel,avel,bmat,prod)
         INTEGER :: n,m,ibdim
         INTEGER :: i,j,indx
         REAL(q) :: vel(3,N),avel(3,3)
         REAL(q) :: bmat(M,ibdim)  
         REAL(q) :: prod(M)

         prod=0._q
         DO i=1,n 
           DO j=1,m
             indx=3*i-2
             prod(j)=prod(j)+vel(1,i)*bmat(j,indx)+vel(2,i)*bmat(j,indx+1)+vel(3,i)*bmat(j,indx+2)
           ENDDO
         ENDDO

         IF (ibdim==(3*N+9)) THEN
           DO j=1,m
             prod(j)=prod(j)+avel(1,1)*bmat(j,3*N+1)+avel(2,1)*bmat(j,3*N+2)+avel(3,1)*bmat(j,3*N+3)
             prod(j)=prod(j)+avel(1,2)*bmat(j,3*N+4)+avel(2,2)*bmat(j,3*N+5)+avel(3,2)*bmat(j,3*N+6)
             prod(j)=prod(j)+avel(1,3)*bmat(j,3*N+7)+avel(2,3)*bmat(j,3*N+8)+avel(3,3)*bmat(j,3*N+9)
           ENDDO
         ENDIF

        END SUBROUTINE vel_bmat_product

        SUBROUTINE make_Z(T_INFO,IBDIM,x,A,ICOORDINATES,MASSES,ZET)        
        !c compute the mass metric tensor (Z) for Bluemoon calculations
        !c due to ortogonalization of BMAT, Z is a diagonal matrix
        !c elements for redundant coordinates are zero
          TYPE(type_info) :: T_INFO
          TYPE(coordstructure) :: ICOORDINATES 
          INTEGER :: IBDIM,i
          REAL(q) :: x(3,T_INFO%NIONS)
          REAL(q) :: A(3,3)
          REAL(q) :: BMAT(ICOORDINATES%NUMINTERNALS,IBDIM)
          REAL(q) :: MASSES(IBDIM)
          REAL(q) :: ZET(ICOORDINATES%NUMINTERNALS)

          ZET=0._q
          
          IF (ICOORDINATES%NUMINTERNALS>0) THEN   
            BMAT=0._q      
            CALL BMATRIX(T_INFO,x,x,TRANSPOSE(A),ICOORDINATES,BMAT,.TRUE.) 
            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS/=0) THEN
                BMAT(i,:)=0._q
              ELSE 
                BMAT(i,:)=BMAT(i,:)/SQRT(MASSES(:))
              ENDIF
            ENDDO
            CALL ICONST_BMAT(BMAT,ICOORDINATES,T_INFO)
          
            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0) THEN
                ZET(i)=SUM(BMAT(i,:)*BMAT(i,:))
                !ZET(i)=SUM(BMAT(i,:)*BMAT(i,:)/MASSES(:))
              END IF
            ENDDO
          ENDIF 

        END SUBROUTINE make_Z

       SUBROUTINE make_Z2(T_INFO,IBDIM,x,A,ICOORDINATES,MASSES,ZET)        
        !c compute the mass metric tensor (Z) for Bluemoon calculations
        !c elements for redundant coordinates are zero
          TYPE(type_info) :: T_INFO
          TYPE(coordstructure) :: ICOORDINATES 
          INTEGER :: IBDIM,i,j
          REAL(q) :: x(3,T_INFO%NIONS)
          REAL(q) :: A(3,3)
          REAL(q) :: BMAT(ICOORDINATES%NUMINTERNALS,IBDIM)
          REAL(q) :: MASSES(IBDIM)
          REAL(q) :: ZET(ICOORDINATES%NUMINTERNALS,ICOORDINATES%NUMINTERNALS)

          ZET=0._q
          
          IF (ICOORDINATES%NUMINTERNALS>0) THEN   
            BMAT=0._q      
            CALL BMATRIX(T_INFO,x,x,TRANSPOSE(A),ICOORDINATES,BMAT,.TRUE.) 
            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS/=0) THEN
                BMAT(i,:)=0._q
              ELSE 
                BMAT(i,:)=BMAT(i,:)/SQRT(MASSES(:))
              ENDIF
            ENDDO
            
          
            DO i=1,ICOORDINATES%NUMINTERNALS
              DO j=1,ICOORDINATES%NUMINTERNALS
                IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0 .AND. ICOORDINATES%COORDSTRUCT(j)%STATUS==0) THEN
                  ZET(i,j)=SUM(BMAT(i,:)*BMAT(j,:))
                  !ZET(i)=SUM(BMAT(i,:)*BMAT(i,:)/MASSES(:))
                END IF
              ENDDO
            ENDDO
          ENDIF 
        END SUBROUTINE make_Z2

        FUNCTION z_det(ZET,ICOORDINATES,IO) RESULT(zdet)
        !c compute det(Z) for a diagonal matrix Z
          INTEGER :: i
          TYPE(coordstructure) :: ICOORDINATES   
          REAL(q) :: ZET(ICOORDINATES%NUMINTERNALS)
          REAL(q) :: zdet
          TYPE (in_struct) :: IO

          zdet=1._q
          DO i=1,ICOORDINATES%NUMINTERNALS
            IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0) THEN
              IF (ABS(ZET(i))<1e-7) THEN
                CALL vtutor%error("Error LBLUEOUT: lineary dependent constraints,|z|=0")
              ENDIF
              !IF (ABS(ZET(i))>1e-5) THEN
              zdet=zdet*ZET(i)
              !ENDIF
            END IF
          ENDDO
        END FUNCTION z_det


        SUBROUTINE make_G(T_INFO,DYN,LATT_CUR,ICOORDINATES,MASSES,IBDIM,ZET_0,ZET_DET_0,GMAT,IO)
          TYPE(type_info) :: T_INFO
          TYPE(dynamics) :: DYN      
          TYPE(latt) :: LATT_CUR 
          TYPE(coordstructure) :: ICOORDINATES  
          REAL(q) :: BMAT(ICOORDINATES%NUMINTERNALS,IBDIM),BMAT_m(ICOORDINATES%NUMINTERNALS,IBDIM)
          REAL(q) :: ZET_0(ICOORDINATES%NUMINTERNALS),ZET_0_inv(ICOORDINATES%NUMINTERNALS)
          REAL(q) :: ZET(ICOORDINATES%NUMINTERNALS)
          REAL(q) :: ZET2(ICOORDINATES%NUMINTERNALS,ICOORDINATES%NUMINTERNALS),ZET2_inv(ICOORDINATES%NUMINTERNALS,ICOORDINATES%NUMINTERNALS)
          REAL(q) :: GMAT(ICOORDINATES%NUMINTERNALS)
          REAL(q) :: MASSES(IBDIM)
          REAL(q) :: POSIOC_b(3,T_INFO%NIONS),POSIOC_cart(3,T_INFO%NIONS)
          REAL(q) :: STEP
          REAL(q) :: Zpd,Zmd,ZET_DET_0
          INTEGER :: IBDIM
          INTEGER :: i,j,k,l,kk,ii
          REAL(q) :: dZx(IBDIM)
          TYPE (in_struct) :: IO



          STEP=1.E-4
          CALL DEAL_XYZ(T_INFO,DYN%POSIOC,DYN%POSIOC,LATT_CUR%A,ICOORDINATES)
          CALL make_Z2(T_INFO,IBDIM,DYN%POSIOC,LATT_CUR%A,ICOORDINATES,MASSES,ZET2)

          BMAT=0._q
          IF (ICOORDINATES%NUMINTERNALS>=1) THEN          
            CALL BMATRIX(T_INFO,DYN%POSIOC,DYN%POSIOC,TRANSPOSE(LATT_CUR%A),ICOORDINATES,BMAT,.TRUE.)  
            DO i=1,ICOORDINATES%NUMINTERNALS
              IF (ICOORDINATES%COORDSTRUCT(i)%STATUS/=0) THEN
                BMAT(i,:)=0._q
              ENDIF
            ENDDO
            !CALL ICONST_BMAT(BMAT,ICOORDINATES,T_INFO)
          ENDIF 

          BMAT_m=0._q
          DO i=1,ICOORDINATES%NUMINTERNALS
            BMAT_m(i,:)=BMAT(i,:)/MASSES(:)
          ENDDO
          
          POSIOC_b=DYN%POSIOC
          POSIOC_cart=MATMUL(LATT_CUR%A,POSIOC_b)
          GMAT=0._q
          
         
          DO ii=1,3*T_INFO%NIONS
            j=MOD(ii,3)
            if (j==0) j=3
            i=(ii-j)/3+1
           
            !c take a forward step
            DYN%POSIOC=POSIOC_cart
            DYN%POSIOC(j,i)=POSIOC_cart(j,i)+STEP
            DYN%POSIOC=MATMUL(TRANSPOSE(LATT_CUR%B),DYN%POSIOC)
            CALL DEAL_XYZ(T_INFO,DYN%POSIOC,DYN%POSIOC,LATT_CUR%A,ICOORDINATES)
            CALL make_Z(T_INFO,IBDIM,DYN%POSIOC,LATT_CUR%A,ICOORDINATES,MASSES,ZET)   
            Zpd=z_det(ZET,ICOORDINATES,IO)    
            
            !c take a backward step
            DYN%POSIOC=POSIOC_cart
            DYN%POSIOC(j,i)=POSIOC_cart(j,i)-STEP
            DYN%POSIOC=MATMUL(TRANSPOSE(LATT_CUR%B),DYN%POSIOC)
            CALL DEAL_XYZ(T_INFO,DYN%POSIOC,DYN%POSIOC,LATT_CUR%A,ICOORDINATES)
            CALL make_Z(T_INFO,IBDIM,DYN%POSIOC,LATT_CUR%A,ICOORDINATES,MASSES,ZET)
            Zmd=z_det(ZET,ICOORDINATES,IO)    
            dZx(ii)=(Zpd-Zmd)/2/STEP
          ENDDO

          GMAT=MATMUL(BMAT_m,dZx)

          ZET2_inv=ZET2
          CALL SVDINVERSE(ZET2_inv,ICOORDINATES%NUMINTERNALS,i)
          GMAT=MATMUL(ZET2_inv,GMAT)

          GMAT=GMAT/2/ZET_DET_0

          !c restore the original values 
          DYN%POSIOC=POSIOC_b
          CALL DEAL_XYZ(T_INFO,DYN%POSIOC,DYN%POSIOC,LATT_CUR%A,ICOORDINATES)
        END SUBROUTINE make_G
 
      END MODULE dynconstr
