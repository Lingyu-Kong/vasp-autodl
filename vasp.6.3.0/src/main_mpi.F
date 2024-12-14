#include "symbol.inc"
      MODULE main_mpi
      USE prec
      USE base
      USE mpimy
      IMPLICIT NONE

!***********************************************************************
! RCS:  $Id: main_mpi.F,v 1.3 2002/04/16 07:28:45 kresse Exp $
!
! This module initializes the communication univeres used in VASP
!
! we have five communicator (see below)
!
!***********************************************************************

      TYPE (communic),TARGET :: COMM_WORLD ! our world wide communicator
      TYPE (communic),TARGET :: COMM_CHAIN ! communication between images
      TYPE (communic),TARGET :: COMM       ! one image communicator
      TYPE (communic),TARGET :: COMM_KINTER ! between k-points communicator
      TYPE (communic),TARGET :: COMM_KIN   ! in k-point communicator
      TYPE (communic),TARGET :: COMM_INTER ! between-band communicator
      TYPE (communic),TARGET :: COMM_INB   ! in-band communicator
#ifdef use_shmem
      TYPE (communic),TARGET :: COMM_SHMEM
      TYPE (communic),TARGET :: COMM_intra_node,COMM_inter_node
#endif
      TYPE (communic),TARGET :: COMM_intra_node_world,COMM_inter_node_world
      INTEGER :: IMAGES=0                  ! total number of images
      INTEGER :: KIMAGES=0                 ! distribution of k-points
      INTEGER :: NCORE_IN_IMAGE1=0         ! if IMAGES is 2 it is possible to specify the number of cores
                                           ! that are in the first image
      CHARACTER(LEN=10) ::  DIR_APP
      INTEGER           ::  DIR_LEN=0
      INTEGER           ::  KPAR,NCORE
#ifdef use_shmem
      INTEGER :: NCSHMEM=1
#else
      INTEGER, PARAMETER :: NCSHMEM=1
#endif

!***********************************************************************
!
! initialize MPI and read entry IMAGES and NPAR from INCAR
! and sub divide communicators
! because redistribution of nodes
! might be done one must do this as early as possible
!
!> @details @ref openmp :
!> forces main_mpi::ncshmem=1 and main_mpi::NCORE=1, and gives
!> additional output.
!
!***********************************************************************

      CONTAINS
#if defined(MPI) || defined(MPI_CHAIN)
      SUBROUTINE INIT_MPI(NPAR,IO)

      USE command_line
      USE reader_tags
      USE incar_reader, ONLY: INCAR_FROM_FILE
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
#ifdef VASP_HDF5
      USE vhdf5_base
#endif
#ifdef _OPENACC
      USE openacc, ONLY : ACC_GET_DEVICE_TYPE,ACC_DEVICE_HOST
#endif
#ifdef USENCCL
      USE mopenacc_struct_def, ONLY : LUSENCCL
#endif
      TYPE (in_struct) :: IO
      INTEGER          :: NPAR
! local
      INTEGER :: IERR, IDUM, UNIT_
      LOGICAL :: LVCAIMAGES, INCAR_SUBDIR
      REAL(q) :: VCAIMAGES
#ifdef VASP_HDF5
      INTEGER :: IH5ERR
#else
      LOGICAL :: INCAR_FOUND
#endif
!$    INTEGER, EXTERNAL :: OMP_GET_MAX_THREADS

#ifdef debug
      CHARACTER*(MPI_MAX_PROCESSOR_NAME) myname
      INTEGER :: I,resultlen,ierror
#endif

      CALL M_init(COMM_WORLD)

      CALL PARSE_COMMAND_LINE(COMM_WORLD%NODE_ME == COMM_WORLD%IONODE)

! in case a demo is compiled
      CALL VASP_DEMO( COMM_WORLD, IO%IU0 )

#ifdef VASP_HDF5
      IH5ERR = VH5_START()
      !
      ! incar, vaspin
      !
      INQUIRE(FILE=INCAR,EXIST=INCAR_FOUND)
      INQUIRE(FILE='vaspin.h5',EXIST=HDF5_FOUND)
      IF (.NOT.HDF5_FOUND) THEN
       IF (.NOT.INCAR_FOUND) CALL vtutor%error("No INCAR or vaspin.h5 found, STOPPING")
      ELSE
        IH5ERR = VH5_FILE_OPEN_READ('vaspin.h5', IH5INFILEID)
        IH5ERR = VH5_GROUP_OPEN(IH5INFILEID, GRP_INPUT, IH5ININPUTGROUP_ID)
        INCAR_F%FROM_HDF5 = .TRUE.
        INCAR_F%ERROR = ""
        ALLOCATE(INCAR_F%TAGS(0))
      ENDIF
#else
      INQUIRE(FILE=INCAR,EXIST=INCAR_FOUND)
      IF (.NOT.INCAR_FOUND) THEN
         CALL vtutor%error("No INCAR found, STOPPING")
      END IF
#endif
      IF (INCAR_FOUND) INCAR_F = INCAR_FROM_FILE(INCAR)
#ifdef VASP_HDF5
      INCAR_F%TO_HDF5 = .TRUE.
#endif

      VCAIMAGES=-1
      CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5,'VCAIMAGES',VCAIMAGES,IERR)

      IF (VCAIMAGES==-1) THEN
         LVCAIMAGES=.FALSE.
      ELSE
         LVCAIMAGES=.TRUE.
      ENDIF

      NCORE_IN_IMAGE1=0
      IF (LVCAIMAGES) THEN
         IMAGES=2
         CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5,'NCORE_IN_IMAGE1',NCORE_IN_IMAGE1,IERR)
      ELSE
         IMAGES=0
         CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5,'IMAGES',IMAGES,IERR)
      ENDIF

      KIMAGES=0
      CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5, 'KIMAGES',KIMAGES,IERR)
      IF (KIMAGES>0) THEN
         IDUM=0
         CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5, 'FOURORBIT',IDUM,IERR)
         IF (IDUM/=1) THEN
            CALL vtutor%error("Distribution of k-points over KIMAGES only works \n in combination with &
               &FOURORBIT=1, sorry, stopping ...")
         ENDIF
      ENDIF
! KPAR division of kpoints,  default to unity in case only 1 k-point
      KPAR=1
      CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5, 'KPAR',KPAR,IERR)
      IF (KPAR>1.AND.KIMAGES>0) THEN
         CALL vtutor%error("Untested combination of FOURORBIT with KPAR \n (k-point parallelization),&
            & sorry, stopping...")
      END IF

#ifdef USENCCL
      CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5, 'LUSENCCL',LUSENCCL,IERR)
#endif

!----------------------------------------------------------------------
! M_divide: creates a 2 dimensional cartesian topology
!  for seperate images or work groups
!  each work group (image) will run VASP independently in one
!  sub directory (01-99) of the current directory
! this mode is required to support either independent calculations
! on parallel machines or the nudged elastic band method
!----------------------------------------------------------------------
      IF (ABS(IMAGES)>0) THEN
         IF (IMAGES==2 .AND. NCORE_IN_IMAGE1/=0) THEN
            CALL M_initc( COMM_WORLD)
            ! M_divide2 calls M_initc for COMM and COMM_CHAIN
            CALL M_divide2( COMM_WORLD, NCORE_IN_IMAGE1, COMM_CHAIN, COMM)
         ELSE
            CALL M_divide( COMM_WORLD, ABS(IMAGES), COMM_CHAIN, COMM, .TRUE. )
            CALL M_initc( COMM_WORLD)
            CALL M_initc( COMM)
            CALL M_initc( COMM_CHAIN)
         ENDIF
      ELSEIF (KIMAGES>0) THEN
         CALL M_divide( COMM_WORLD, KIMAGES, COMM_CHAIN, COMM, .TRUE. )
         CALL M_initc( COMM_WORLD)
         CALL M_initc( COMM)
         CALL M_initc( COMM_CHAIN)
      ELSE
         CALL M_initc( COMM_WORLD)
         COMM=COMM_WORLD
         COMM_CHAIN=COMM
      ENDIF
#ifdef _OPENMP
!$    IF ( COMM_WORLD%NODE_ME == 1 ) THEN
#ifndef VASPDemo
!$      WRITE(IO%IU0,*)'----------------------------------------------------'
!$      WRITE(IO%IU0,*)'   OOO  PPPP  EEEEE N   N M   M PPPP'
!$      WRITE(IO%IU0,*)'  O   O P   P E     NN  N MM MM P   P'
!$      WRITE(IO%IU0,*)'  O   O PPPP  EEEEE N N N M M M PPPP   -- VERSION'
!$      WRITE(IO%IU0,*)'  O   O P     E     N  NN M   M P'
!$      WRITE(IO%IU0,*)'   OOO  P     EEEEE N   N M   M P'
!$      WRITE(IO%IU0,*)'----------------------------------------------------'
#endif
!$      WRITE(IO%IU0,'(" running ",I4," mpi-ranks, with ",I4," threads/rank")') COMM_WORLD%NCPU,OMP_GET_MAX_THREADS()
!$    ENDIF
#else
      IF ( COMM_WORLD%NODE_ME == 1 ) &
        WRITE(IO%IU0,'(" running on ",I4," total cores")') COMM_WORLD%NCPU
#endif
      IF ( COMM_WORLD%NODE_ME == 1 .AND. (ABS(IMAGES)>0.OR.KIMAGES>0) ) &
        WRITE(IO%IU0,'(" each image running on ",I4," cores")') COMM%NCPU

      IF (KPAR>=1) THEN
!----------------------------------------------------------------------
! M_divide: creates a 2 dimensional cartesian topology within one
!  work group (image)
!  this is required for simultaneous distribution over k-points and bands
!----------------------------------------------------------------------
         CALL M_divide( COMM, KPAR, COMM_KINTER, COMM_KIN, .FALSE.)
         CALL M_initc( COMM)   ! probably not required but who knows
         CALL M_initc( COMM_KINTER)
         CALL M_initc( COMM_KIN)
         IF ( COMM_WORLD%NODE_ME == 1 ) &
         WRITE(IO%IU0,'(" distrk:  each k-point on ",I4," cores, ",I4," groups")') &
                        COMM_KIN%NCPU,COMM_KINTER%NCPU
      ELSE
         COMM_KINTER = COMM
         COMM_KIN    = COMM
      ENDIF

! NCORE species onto how many cores a band is distributed
! often this values can be now set to the number of cores per node
! this is more handy than NPAR in most cases
      NCORE=1
      CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5, 'NCORES_PER_BAND',NCORE,IERR)
      CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5, 'NCORE',NCORE,IERR)
      NCORE = MAX(MIN(COMM_KIN%NCPU, NCORE), 1)
      IF (MOD(COMM_KIN%NCPU, NCORE) /= 0) NCORE = 1

! NPAR number of bands distributed over processors, defaults to
! COMM%NCPU/NCORE
      NPAR = 0
      CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5, 'NPAR',NPAR,IERR)
      IF ((NPAR < 1).OR.(NPAR > COMM_KIN%NCPU).OR.(MOD(COMM_KIN%NCPU, NPAR) /= 0)) THEN
         NPAR = MAX(COMM_KIN%NCPU / NCORE, 1)
      ELSE
         NCORE = COMM_KIN%NCPU / NPAR
      END IF

#ifdef _OPENMP
! When using more than one OpenMP thread we only support NCORE=1 for now
      IF (OMP_GET_MAX_THREADS()>1) THEN
         NCORE=1; NPAR=MAX(COMM_KIN%NCPU/NCORE,1)
      ENDIF
#endif
#ifdef _OPENACC
! When using a GPU accelerator one should really only use NCORE=1
      IF (ACC_GET_DEVICE_TYPE()/=ACC_DEVICE_HOST) THEN
         NCORE=1; NPAR=MAX(COMM_KIN%NCPU/NCORE,1)
      ENDIF
#endif

      IF (NPAR>=1) THEN
!----------------------------------------------------------------------
! M_divide: creates a 2 dimensional cartesian topology within one
!  work group (image)
!  this is required for simultaneous distribution over bands and plane
!  wave coefficients
!  communicators are created for inter-band in intra-band communication
!  NPAR is the number of bands over which is parallelized
! the resulting layout will be the following (NPAR=4, NCPU=8):
! (MPI uses allways Row-major layout)
! wave1           0 (0,0)        1 (0,1)
! wave2           2 (1,0)        3 (1,1)
! wave3           4 (2,0)        5 (2,1)
! wave4           6 (3,0)        7 (3,1)
! wave5           0 (0,0)        1 (0,1)
! etc.
!
! the sub-communicators COMM_INB are one dimensional communicators
! which allow communication within one row i.e. nodes are grouped to
!  0-1     2-3        4-5      6-7
! for shortness these groups will we called in-band-groups
! there communicators are called in-band-communicator
!
! the sub-communicators COMM_INTER are one dimensional communicators
! which allow communication within one column i.e nodes are grouped
!  0-2-4-6     and    1-3-5-7
! these groups will we called inter-band-groups
!
! the most complicated thing is the FFT of soft chargedensities
! the following algorithm is used:
! the soft chargedensity is calculated in real space,
! charge from all bands is merged to processor 0 and 1 using
! the COMM_INTER
! than a FFT involving all processors is done (on processors 2-7
! no components exist in real space)
! the final result in reciprocal space is defined on all processors
! (see SET_RL_GRID for more information)
!----------------------------------------------------------------------
         CALL M_divide( COMM_KIN, NPAR, COMM_INTER, COMM_INB, .FALSE.)
         CALL M_initc( COMM_KIN)   ! propably not required but who knows
         CALL M_initc( COMM_INTER)
         CALL M_initc( COMM_INB)
         IF ( COMM_WORLD%NODE_ME == 1 ) &
         WRITE(IO%IU0,'(" distr:  one band on ",I4," cores, ",I4," groups")') &
                        COMM_INB%NCPU,COMM_INTER%NCPU
      ELSE
         COMM_INTER = COMM
         COMM_INB   = COMM
      ENDIF


#ifdef use_shmem
#ifndef _OPENMP
      CALL PROCESS_INCAR(IO%LOPEN,IO%IU0,IO%IU5, 'NCSHMEM',NCSHMEM,IERR)
#endif
      CALL M_divide_shmem(COMM_WORLD,COMM_INTER,NCSHMEM,COMM_SHMEM)
      CALL M_divide_intra_inter_node(COMM_WORLD,COMM_INTER,COMM_intra_node,COMM_inter_node)
#endif

      CALL M_divide_intra_inter_node(COMM_WORLD,COMM_WORLD,COMM_intra_node_world,COMM_inter_node_world)

      IF (COMM%NODE_ME /= COMM%IONODE) IO%IU6 = -1
      IF (COMM%NODE_ME /= COMM%IONODE) IO%IU0 = -1

      IF (KIMAGES>0.AND.COMM_WORLD%NODE_ME/=COMM_WORLD%IONODE) IO%IU6 = -1
      IF (KIMAGES>0.AND.COMM_WORLD%NODE_ME/=COMM_WORLD%IONODE) IO%IU0 = -1

      IF (IMAGES<0.AND.COMM_WORLD%NODE_ME/=COMM_WORLD%IONODE) IO%IU6 = -1
      IF (IMAGES<0.AND.COMM_WORLD%NODE_ME/=COMM_WORLD%IONODE) IO%IU0 = -1

      IF (KPAR>1.AND.COMM%NODE_ME/=COMM%IONODE) IO%IU6 = -1
      IF (KPAR>1.AND.COMM%NODE_ME/=COMM%IONODE) IO%IU0 = -1

      IF ( IO%IU0/=-1 .AND. IO%IU6 == -1) THEN
         CALL vtutor%bug("internal ERROR: io-unit problem " // str(IO%IU0) // " " // str(IO%IU6), __FILE__, __LINE__)
      ENDIF

      CALL MAKE_DIR_APP(COMM_CHAIN%NODE_ME)
      IMAGES = ABS(IMAGES)

! if all nodes should write (giving complete mess) do not use the
! following line
      IF (COMM_WORLD%NODE_ME /= COMM_WORLD%IONODE .AND. IO%IU0>0) THEN
         OPEN(UNIT=IO%IU0,FILE=DIR_APP(1:DIR_LEN)//'stdout',STATUS='UNKNOWN')
      ENDIF
!----------------------------------------------------------------------
! try to go to subdir INCAR's for the rest
!----------------------------------------------------------------------
      ! TODO fix images for HDF5
      IF (DIR_LEN > 0) THEN
         INQUIRE(FILE=DIR_APP(1:DIR_LEN)//INCAR, EXIST=INCAR_SUBDIR)
         IF (INCAR_SUBDIR) THEN
            INCAR=DIR_APP(1:DIR_LEN)//INCAR
            IF (IO%IU0>=0) WRITE(IO%IU0,*) 'using from now: ',INCAR
            INCAR_F = INCAR_FROM_FILE(INCAR)
         ELSE
            INCAR_F%IMAGE = COMM_CHAIN%NODE_ME
         ENDIF
      ENDIF

#ifdef VASP_HDF5
      !
      ! if no INCAR, try to close vaspin.h5
      !
      IF (.NOT. INCAR_FOUND .AND. HDF5_FOUND) THEN
        IH5ERR = VH5_GROUP_CLOSE(IH5ININPUTGROUP_ID)
        IH5ERR = VH5_FILE_CLOSE(IH5INFILEID)
        IH5ERR = VH5_END()
      ENDIF
#endif
#ifdef debug
      DO I=1,COMM_WORLD%NCPU
         IF (COMM_WORLD%NODE_ME==I) THEN
            CALL MPI_Get_processor_name(myname,resultlen,ierror)
            WRITE(*,'(5(A,I4,X),A,A)') 'world id:',COMM_WORLD%NODE_ME,'inb:',COMM_INB%NODE_ME,'int:',COMM_INTER%NODE_ME, &
                                       'intra:',COMM_intra_node%NODE_ME,'inter:',COMM_inter_node%NODE_ME,'cpu:',myname(1:resultlen)
         ENDIF
         CALL MPI_barrier(COMM_WORLD%MPI_COMM,ierror)
      ENDDO
#endif
      RETURN
      END SUBROUTINE
#endif
!***********************************************************************
!
! make the  directory entry which is used for
! fileio
!
!***********************************************************************

      SUBROUTINE MAKE_DIR_APP(node)
      INTEGER node

! in principle one can chose here any string one wants to use
! only DIR_LEN must be adjusted
      WRITE (DIR_APP  , "(I1,I1,'/')") MOD(node/10,10),MOD(node,10)
      IF (IMAGES<=0.OR.KIMAGES>0) THEN
         DIR_LEN=0
      ELSE
         DIR_LEN=3
      ENDIF

      END SUBROUTINE MAKE_DIR_APP

!***********************************************************************
!
! once unit 6 is open write number of nodes and
! all other parameters
!
!> @details @ref openmp :
!> gives additional output.
!
!***********************************************************************
      SUBROUTINE WRT_DISTR(IU6)
      USE openmp, ONLY : omp_nthreads
      INTEGER IU6
#ifdef MPI
      IF (IU6>=0) THEN
#ifdef _OPENMP
        WRITE(IU6,'(" running ",I4," mpi-ranks, with ",I4," threads/rank")') COMM_WORLD%NCPU,omp_nthreads
#else
        WRITE(IU6,'(" running on ",I4," total cores")') COMM_WORLD%NCPU
#endif
        IF (IMAGES>0 ) &
        WRITE(IU6,'(" each image running on ",I4," cores")') COMM%NCPU
        WRITE(IU6,'(" distrk:  each k-point on ",I4," cores, ",I4," groups")') &
                        COMM_KIN%NCPU,COMM_KINTER%NCPU
        WRITE(IU6,'(" distr:  one band on NCORE=",I4," cores, ",I4," groups")') &
                        COMM_INB%NCPU,COMM_INTER%NCPU

      ENDIF
#else
      WRITE(IU6,*)'serial version'
#endif
      END SUBROUTINE WRT_DISTR


#if defined(MPI) || defined(MPI_CHAIN)
!***********************************************************************
! VASP DEMO
!***********************************************************************
      SUBROUTINE VASP_DEMO( COMM_WORLD, IU0 )
      USE tutor, ONLY: vtutor
      TYPE( communic )  :: COMM_WORLD
      INTEGER           :: IU0
      ! local
      INTEGER, PARAMETER :: NMAX = 4
!$    INTEGER, EXTERNAL :: OMP_GET_MAX_THREADS

      ! immeadiate return if not DEMO
#ifndef VASPDemo
      RETURN
#endif
      ! write a banner to indicate that demo is executed
      IF ( COMM_WORLD%NODE_ME == 1  ) THEN
         WRITE(IU0,*)'-----------------------------------------------------------------------------'
#ifdef _OPENMP
         IF( COMM_WORLD%NCPU > NMAX .OR. OMP_GET_MAX_THREADS() > NMAX) THEN
#else
         IF( COMM_WORLD%NCPU > NMAX ) THEN
#endif
         WRITE(IU0,*)'         This                                                                 '
         ENDIF
         WRITE(IU0,*)'                                                                              '
         WRITE(IU0,*)'          V    V   AA    SSSS  PPPPP    DDDDD  EEEEEE  MM MM   OOOO           '
         WRITE(IU0,*)'          V    V  A  A  S    S P    P   D    D E      M  M  M O    O          '
         WRITE(IU0,*)'          V    V A    A  SSS   PPPP     D    D EEEEE  M  M  M O    O          '
         WRITE(IU0,*)'          V    V AAAAAA      S P        D    D E      M  M  M O    O          '
         WRITE(IU0,*)'           V  V  A    A S    S P        D    D E      M     M O    O          '
         WRITE(IU0,*)'            VV   A    A  SSSS  P        DDDDD  EEEEEE M     M  OOOO           '
         WRITE(IU0,*)'                                                                              '
      ENDIF

      ! if requirements are not fullfiled, stop here
#ifdef _OPENMP
      IF( COMM_WORLD%NCPU > NMAX .OR. OMP_GET_MAX_THREADS() > NMAX) THEN
#else
      IF( COMM_WORLD%NCPU > NMAX ) THEN
#endif
         IF ( COMM_WORLD%NODE_ME == 1  ) THEN
            WRITE(IU0,'(A,I2,A,I2,A)') '          is restricted to maximally',NMAX,&
                                              ' MPI Ranks and',NMAX,' OpenMP threads!'
            WRITE(IU0,*)'-----------------------------------------------------------------------------'
         ENDIF
         CALL vtutor%stopCode()
      ENDIF
      IF ( COMM_WORLD%NODE_ME == 1  ) &
         WRITE(IU0,*)'-----------------------------------------------------------------------------'
      END SUBROUTINE
#endif

      END MODULE
