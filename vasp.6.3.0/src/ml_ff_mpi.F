#if defined(MPI) || defined(MPI_CHAIN)
#include "symbol.inc"
!****************************************************************************************************
! Module for MPI
!****************************************************************************************************

      MODULE MPI_DATA

        USE MPI_HELP
        USE ML_FF_CONSTANT
        USE ML_FF_PREC
        IMPLICIT NONE

        CONTAINS

        SUBROUTINE M_initc (COMM)
           USE ML_FF_STRUCT, ONLY: ML_MPI_PAR
           IMPLICIT NONE
           TYPE(ML_MPI_PAR) :: COMM
           INTEGER          :: id_in_group
           INTEGER          :: ierror
           CALL MPI_comm_rank( COMM%MPI_COMM, id_in_group, ierror )
           CALL MPI_comm_size( COMM%MPI_COMM, COMM%NCPU, ierror )
           COMM%NODE_ME = id_in_group !+ 1
           CALL MPI_barrier( COMM%MPI_COMM, ierror )
        END SUBROUTINE M_initc

        SUBROUTINE M_divide_intra_inter_node(COMM_WORLD,COMM,COMM_intra,COMM_inter)
           USE ml_ff_c2f_interface
           USE ML_FF_STRUCT, ONLY: ML_MPI_PAR
           USE PRINT2TERMINAL
           IMPLICIT NONE
           TYPE (ML_MPI_PAR) :: COMM_WORLD,COMM,COMM_intra,COMM_inter,COMM_test
           ! local variables
           CHARACTER*(MPI_MAX_PROCESSOR_NAME) myname
           CHARACTER*(MPI_MAX_PROCESSOR_NAME) pname
           INTEGER :: COMM_group
           INTEGER :: COMM_intra_group
           INTEGER :: COMM_inter_group
           INTEGER :: group(0:COMM%NCPU-1)
           INTEGER :: resultlen
           INTEGER :: myid
           INTEGER :: ierror
           INTEGER :: I
           INTEGER :: IERR
           INTEGER :: IGRP
           ! shmem variables for sanity check
           INTEGER(c_int)    :: shmid
           TYPE(c_ptr)       :: address
           INTEGER(c_size_t) :: k

           INTEGER*4, POINTER :: IDSHMEM(:)
           INTEGER :: ID
           ! attempt automatic division of COMM
           CALL MPI_comm_split_type(COMM%MPI_COMM,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,COMM_intra%MPI_COMM,ierror)
           CALL MPI_comm_group(COMM_intra%MPI_COMM,COMM_intra_group,ierror)
           CALL MPI_comm_group(COMM%MPI_COMM,COMM_group,ierror)
           CALL M_initc(COMM_intra)
#ifdef use_shmem
           ! sanity check
           k=1
!           IF (COMM_intra%NODE_ME==1) CALL getshmem(INT(8*k,KIND=c_size_t),shmid)
           IF (COMM_intra%NODE_ME==0) CALL getshmem(INT(8*k,KIND=c_size_t),shmid)
           CALL M_BCAST(COMM_intra,shmid,0)
           CALL attachshmem(shmid,address)
           CALL c_f_pointer(address,IDSHMEM,[k])
           
!           IF (COMM_intra%NODE_ME==1) THEN
           IF (COMM_intra%NODE_ME==0) THEN
               ID=COMM_WORLD%NODE_ME; IDSHMEM=ID
           ENDIF
           CALL M_BCAST(COMM_intra,ID,0)

           ierror=0
           IF (ID/=IDSHMEM(1)) ierror=1
           CALL MPI_ALLREDUCE(MPI_IN_PLACE,ierror,1,MPI_INTEGER,MPI_SUM,COMM_WORLD%MPI_COMM,IERR)

           CALL detachshmem(address)
!           IF (COMM_intra%NODE_ME==1) CALL destroyshmem(shmid)
           IF (COMM_intra%NODE_ME==0) CALL destroyshmem(shmid)

           IF (ierror>0) THEN
              CALL ERROR_PRINT2TERMINAL("M_divide_intra_inter_node: &
              not all procs in COMM_intra seem to be on the same node.")
           ENDIF
#endif

           IGRP=0
           CALL MPI_group_rank(COMM_intra_group,myid,ierror)
           DO I=1,COMM%NCPU
              IF (I==(COMM%NODE_ME+1)) ID=myid
              CALL M_BCAST(COMM,ID,I-1)
              IF (ID==myid) THEN
                 group(IGRP)=I-1
                 IGRP=IGRP+1
              ENDIF
           ENDDO

           CALL MPI_group_incl(COMM_group,IGRP,group(0),COMM_inter_group,ierror)
           CALL MPI_comm_create(COMM%MPI_COMM,COMM_inter_group,COMM_inter%MPI_COMM,ierror)

           CALL M_initc(COMM_inter)

        END SUBROUTINE M_divide_intra_inter_node

!******************************************************************************************
! SUBROUTINE CALC_COLOR: 
! This subroutine calculates the coloer for subcommunicators from an existing communicator
!******************************************************************************************
        SUBROUTINE CALC_COLOR(NPROC, N_LAST, N_FIRST, MYPROC, COLOR)
          IMPLICIT NONE
          INTEGER, INTENT(IN)  :: NPROC                ! number of processers in old communicator
          INTEGER, INTENT(IN)  :: N_LAST               ! last element for calculation of color
          INTEGER, INTENT(IN)  :: N_FIRST              ! first element for calculation of color
          INTEGER, INTENT(IN)  :: MYPROC               ! process rank in old communicator
          INTEGER, INTENT(OUT) :: COLOR                ! the color for the new communicator
          INTEGER              :: I
          INTEGER              :: J
          INTEGER              :: NP_PARTIAL           ! helper variable (number of processors for a given color in new communicator)
          NP_PARTIAL = NPROC/(N_LAST-N_FIRST+1)
          IF (NP_PARTIAL.GT.0) THEN
             DO I = N_FIRST, N_LAST
                DO J = 1, NP_PARTIAL
                   IF (MYPROC.EQ.(I-N_FIRST)*NP_PARTIAL+J-1) THEN
                      COLOR = I - N_FIRST
                   ENDIF
                ENDDO
             ENDDO
             DO I = (N_LAST-N_FIRST)*NP_PARTIAL+NP_PARTIAL, NPROC
                IF (I.EQ.MYPROC) THEN
                   COLOR = N_LAST-(I-((N_LAST-N_FIRST)*NP_PARTIAL+NP_PARTIAL))-1
                ENDIF
             ENDDO
          ELSE
             COLOR = MYPROC
          ENDIF
        END SUBROUTINE CALC_COLOR

!******************************************************************************************
! SUBROUTINE MPI_MAKE_NEW_COMMUNICATOR: 
! This subroutine splits an existing communicator into subcommunicators. 
! The color identity for each subcommunicator has to be calculated previously and is 
! provided as an input argument to the subroutine.
!******************************************************************************************
        SUBROUTINE MPI_MAKE_NEW_COMMUNICATOR(OLD_COMM, OLD_MYID, COLOR, NEW_COMM, NEW_NPE, NEW_MYID)
          IMPLICIT NONE
          INTEGER, INTENT(IN)  :: OLD_COMM           ! old communicator
          INTEGER, INTENT(IN)  :: OLD_MYID           ! process rank in old communicator
          INTEGER, INTENT(IN)  :: COLOR              ! the previously calculated color for the process
          INTEGER, INTENT(OUT) :: NEW_COMM           ! new communicator
          INTEGER, INTENT(OUT) :: NEW_NPE            ! number of processes in new communicator
          INTEGER, INTENT(OUT) :: NEW_MYID           ! process rank in new communicator
          INTEGER              :: IERR               ! error handling variable
          CALL MPI_COMM_RANK(OLD_COMM, OLD_MYID, IERR)                   ! determine rank of old communicator
          CALL MPI_COMM_SPLIT(OLD_COMM, COLOR, OLD_MYID, NEW_COMM, IERR) ! make new communicator
          CALL MPI_COMM_SIZE(NEW_COMM, NEW_NPE, IERR)                    ! determine numb of processes in new communicator
          CALL MPI_COMM_RANK(NEW_COMM, NEW_MYID, IERR)                   ! determine process rank in new communictor
        END SUBROUTINE MPI_MAKE_NEW_COMMUNICATOR

!******************************************************************************************
! SUBROUTINE MPI_DESTROY_COMMUNICATOR
! This subroutine deallocates given communicator.
!******************************************************************************************
        SUBROUTINE MPI_DESTROY_COMMUNICATOR(MY_COMM)
           IMPLICIT NONE
           INTEGER :: MY_COMM
           INTEGER :: IERR
           CALL MPI_COMM_FREE(MY_COMM,IERR)
        END SUBROUTINE MPI_DESTROY_COMMUNICATOR


#ifdef scaLAPACK
!******************************************************************************************
! SUBROUTINE SCALAPACK_INIT: 
! This subroutine sets up a process grid necessary for SCALAPACK.
!******************************************************************************************

        SUBROUTINE SCALAPACK_INIT(NPE,COMMUNIC,NDIM_IN,NDIM_OUT,ICTXT,NPROW,NPCOL,MYROW,MYCOL,NPROCS)
          USE ML_FF_STRUCT, ONLY : ML_MPI_PAR
          IMPLICIT NONE
          TYPE(ML_MPI_PAR), INTENT(IN) :: COMMUNIC
          INTEGER, INTENT(IN)  :: NPE                ! Desired number of
                                                     ! processors in Scalapack grid
          INTEGER, INTENT(IN)  :: NDIM_IN            ! Input dimension of scaLAPACK grid
          INTEGER, INTENT(OUT) :: NDIM_OUT           ! Input dimension of scaLAPACK grid
          INTEGER, INTENT(OUT) :: ICTXT              ! SCALAPACK context handle 
          INTEGER, INTENT(OUT) :: NPROW              ! number of rows in SCALAPACK process grid
          INTEGER, INTENT(OUT) :: NPCOL              ! number of columns in SCALAPACK process grid
          INTEGER, INTENT(OUT) :: MYROW              ! current rank of row in SCALAPACK process grid
          INTEGER, INTENT(OUT) :: MYCOL              ! current rank of column in SCALAPACK process grid
          INTEGER, INTENT(OUT) :: NPROCS             ! number of processors in context after making 
                                                     ! scaLAPACK grid. Should be equal to NPE!
! Local variables
          INTEGER :: I
          INTEGER :: J
          INTEGER :: K
          LOGICAL :: LCHECK
          INTEGER :: IERR
          INTEGER :: NPROW_TMP
          INTEGER :: NPCOL_TMP
          INTEGER :: MYROW_TMP
          INTEGER :: MYCOL_TMP
          INTEGER :: TMPCONTXT
          INTEGER :: MY_BLACS_ID
          INTEGER, ALLOCATABLE :: IMAP(:,:)
          INTEGER, ALLOCATABLE :: BLACS_IDS(:)

! External functions and subroutines
          INTEGER,EXTERNAL :: BLACS_PNUM
          EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDMAP

! Control whether 1D or 2D grid is used
          NDIM_OUT=NDIM_IN

! Set number of rowas and columns (NPROW and NPCOL)
          LCHECK = .FALSE.
          IF(NDIM_IN.EQ.1) THEN
             NPROW=1
             NPCOL=NPE
          ELSE
             DO I=1,1000
                DO J=1,I
                   IF ( I*J.EQ.NPE) THEN
                      NPCOL = J
                      NPROW = I
                      LCHECK = .TRUE.
                      EXIT
                   ENDIF
                ENDDO
                IF (LCHECK) THEN
                   EXIT
                ENDIF
             ENDDO
          ENDIF

! Allocate neccessary helping arrays IMAP and BLACS_IDS
          ALLOCATE(IMAP(NPROW,NPCOL))
          ALLOCATE(BLACS_IDS(NPE))

! Temporarily map all processes into 1 x NPROCS grid
          TMPCONTXT = COMMUNIC%MPI_COMM
          CALL MPI_BARRIER(COMMUNIC%MPI_COMM, IERR)
          CALL BLACS_GRIDINIT(TMPCONTXT, 'Row', 1, NPE)
          CALL MPI_BARRIER(COMMUNIC%MPI_COMM, IERR)

! Obtain information on one dimensional grid
          CALL BLACS_GRIDINFO( TMPCONTXT, NPROW_TMP, NPCOL_TMP, MYROW_TMP, MYCOL_TMP)
          CALL MPI_BARRIER(COMMUNIC%MPI_COMM, IERR)
          CALL MPI_BARRIER(COMMUNIC%MPI_COMM, IERR)

! Free temporary context
          CALL MPI_barrier(COMMUNIC%MPI_COMM, IERR)
 
          IF (NPE.NE.1) THEN
             MY_BLACS_ID = BLACS_PNUM(TMPCONTXT, MYROW_TMP, MYCOL_TMP)
             CALL BLACS_GRIDEXIT(TMPCONTXT)

! Distribute information on Grid 
             CALL MPI_ALLGATHER(MY_BLACS_ID,1,MPI_INTEGER,BLACS_IDS(1),1, &
                                MPI_INTEGER,COMMUNIC%MPI_COMM,IERR)
             IF (IERR.NE.0) THEN
                IF (COMMUNIC%NODE_ME.EQ.0) THEN
                   WRITE(*,666)
                   WRITE(*,*) "MLFF: Error in set up of communicator."
                   WRITE(*,*) "Exiting..."
                ENDIF
                CALL MPI_BARRIER(COMMUNIC%MPI_COMM,IERR)
                CALL MPI_ABORT(COMMUNIC%MPI_COMM,IERR)
             ENDIF
             K=1
             DO I = 1, NPROW
                DO J = 1, NPCOL
                   IMAP(I, J) = BLACS_IDS( K )
                   K = K + 1
                END DO
             END DO
          ELSE
             IF (MYROW_TMP.NE.-1.AND.MYCOL_TMP.NE.-1) THEN
                MY_BLACS_ID = BLACS_PNUM(TMPCONTXT, MYROW_TMP,MYCOL_TMP)
                CALL BLACS_GRIDEXIT(TMPCONTXT)
             ELSE
                MY_BLACS_ID = 0
             ENDIF
             CALL MPI_ALLREDUCE(MPI_IN_PLACE,MY_BLACS_ID,1,MPI_INTEGER, &
                                MPI_SUM,COMMUNIC%MPI_COMM,IERR)
             IMAP(1,1)=MY_BLACS_ID
          ENDIF
          
          CALL MPI_barrier(COMMUNIC%MPI_COMM, IERR)
          CALL MPI_barrier(COMMUNIC%MPI_COMM, IERR)
! Make individual process grid
          ICTXT = COMMUNIC%MPI_COMM
          CALL BLACS_GRIDMAP(ICTXT, IMAP, NPROW, NPROW, NPCOL )
          MYROW=-1 ! preliminary initialization of rows to -1
          MYCOL=-1 ! same for columns
          CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)   ! determine basic information of SCALAPACK grid
          NPROCS = NPROW*NPCOL
          CALL MPI_barrier(COMMUNIC%MPI_COMM, IERR)
          CALL MPI_barrier(COMMUNIC%MPI_COMM, IERR)

! Deallocate neccessary helping arrays IMAP and BLACS_IDS
          DEALLOCATE(IMAP)
          DEALLOCATE(BLACS_IDS)

! Format.
666       FORMAT(/' ------------------------------------------------', &
                  '----------------------------- '/, &
                  '|                                                ', &
                  '                             |'/, &
                  '|     EEEEEEE  RRRRRR   RRRRRR   OOOOOOO  RRRRRR ', &
                  '     ###     ###     ###     |'/, &
                  '|     E        R     R  R     R  O     O  R     R', &
                  '     ###     ###     ###     |'/, &
                  '|     E        R     R  R     R  O     O  R     R', &
                  '     ###     ###     ###     |'/, &
                  '|     EEEEE    RRRRRR   RRRRRR   O     O  RRRRRR ', &
                  '      #       #       #      |'/, &
                  '|     E        R   R    R   R    O     O  R   R  ', &
                  '                             |'/, &
                  '|     E        R    R   R    R   O     O  R    R ', &
                  '     ###     ###     ###     |'/, &
                  '|     EEEEEEE  R     R  R     R  OOOOOOO  R     R', &
                  '     ###     ###     ###     |'/, &
                  '|                                                ', &
                  '                             |')

        END SUBROUTINE SCALAPACK_INIT

!******************************************************************************************
! SUBROUTINE SCALAPACK_SET_DESCRIPTOR: 
! This subroutine sets all the necessary variables for one descriptor used by SCALAPACK.
!******************************************************************************************
        SUBROUTINE SCALAPACK_SET_DESCRIPTOR(DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT, MYROW, MYCOL, NROW, NCOL, LLD, LOC_M, LOC_N)
           IMPLICIT NONE
           INTEGER, INTENT(IN)  ::  M            ! number of rows in global array
           INTEGER, INTENT(IN)  ::  N            ! number of columns in global array
           INTEGER, INTENT(IN)  ::  MB           ! blocking factor for rows
           INTEGER, INTENT(IN)  ::  NB           ! blocking factor for columns
           INTEGER, INTENT(IN)  ::  IRSRC        ! process row of first row of matrix distribution
           INTEGER, INTENT(IN)  ::  ICSRC        ! process column of first column of matrix distribution
           INTEGER, INTENT(IN)  ::  ICTXT        ! BLACS context handle
           INTEGER, INTENT(IN)  ::  MYROW        ! row process index in SCALAPACK process grid
           INTEGER, INTENT(IN)  ::  MYCOL        ! column process index in SCALAPACK process grid
           INTEGER, INTENT(IN)  ::  NROW         ! number of processes in row of SCALAPACK process grid
           INTEGER, INTENT(IN)  ::  NCOL         ! number of processes in column of SCALAPACK process grid
           INTEGER, INTENT(OUT) ::  DESC(:) !(9)      ! descriptor array for SCALAPACK
           INTEGER, INTENT(OUT) ::  LLD          ! Leading dimension of local array(MAX(1, LOCr(M))
           INTEGER, INTENT(OUT) ::  LOC_M        ! local size of M specific to one process
           INTEGER, INTENT(OUT) ::  LOC_N        ! local size of N specific to one process
! Local variables
           INTEGER              ::  INFO
           CALL SCALAPACK_CALCULATE_LOCAL_SIZE(M, MB, MYROW, IRSRC, NROW, LOC_M)  ! calculate local size for row
           CALL SCALAPACK_CALCULATE_LOCAL_SIZE(N, NB, MYCOL, ICSRC, NCOL, LOC_N)  ! calculated local size for column
           LLD = MAX(1, LOC_M)
           CALL DESCINIT(DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT, LOC_M, INFO)    ! set SCALAPACK descriptor
        END SUBROUTINE SCALAPACK_SET_DESCRIPTOR

!******************************************************************************************
! SUBROUTINE SCALAPACK_CALCULATE_LOCAL_SIZE: 
! This subroutine calculates local size in one dimension for a block-cyclic distributed array used in SCALAPACK
!******************************************************************************************
        SUBROUTINE SCALAPACK_CALCULATE_LOCAL_SIZE(SIZE_GLOBAL, BLOCKSIZE, MYCOORD, FIRSTCOORD, NP, SIZE_LOCAL)
          IMPLICIT NONE
          INTEGER, INTENT(IN)  :: SIZE_GLOBAL    ! global size of a given dimension (row or colun) in array
          INTEGER, INTENT(IN)  :: BLOCKSIZE      ! blocksize in given dimension
          INTEGER, INTENT(IN)  :: MYCOORD        ! process rank in given dimension
          INTEGER, INTENT(IN)  :: FIRSTCOORD     ! coordinate of process that posses first row or column fo distributed matrix
          INTEGER, INTENT(IN)  :: NP             ! number of processes in given dimension
          INTEGER, INTENT(OUT) :: SIZE_LOCAL     ! local size of given dimension in the array
          INTEGER              :: NUMROC
          EXTERNAL NUMROC
          SIZE_LOCAL = MAX(1, NUMROC(SIZE_GLOBAL, BLOCKSIZE, MYCOORD, FIRSTCOORD, NP))
        END SUBROUTINE SCALAPACK_CALCULATE_LOCAL_SIZE

!******************************************************************************************
! SUBROUTINE GLOBAL2LOCAL: 
! This function calculates the local coordinate index in one dimension of a global index
! of an array that is blockcyclic distributed
!******************************************************************************************
        INTEGER FUNCTION GLOBAL2LOCAL(IGLOBAL,NP,BLOCKSIZE)
           IMPLICIT NONE
           INTEGER, INTENT(IN) ::  IGLOBAL    ! global coordinate index in given dimension
           INTEGER, INTENT(IN) ::  NP         ! number of processes in given dimension
           INTEGER, INTENT(IN) ::  BLOCKSIZE  ! blocking factor in given dimension
           INTEGER ILOCAL
           ILOCAL=(IGLOBAL-1)/(NP*BLOCKSIZE)
           ILOCAL=ILOCAL*BLOCKSIZE + MOD(IGLOBAL-1,BLOCKSIZE) + 1
           GLOBAL2LOCAL=ILOCAL
        END FUNCTION GLOBAL2LOCAL

!******************************************************************************************
! SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDSYEV_INIT
! This subroutine sets LWORK as -1 and allocate WORK for calculating the minimum size of 
! work arrays
!******************************************************************************************
        SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDSYEV_INIT(SCAWORK)
           USE ML_FF_STRUCT, ONLY : SCALA_WORK 
           IMPLICIT NONE
           TYPE(SCALA_WORK) :: SCAWORK
           SCAWORK%LWORK=-1
           IF (ALLOCATED(SCAWORK%WORK)) THEN
             DEALLOCATE(SCAWORK%WORK)
           ENDIF
           ALLOCATE(SCAWORK%WORK(1))
        END SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDSYEV_INIT

!******************************************************************************************
! SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDSYEV_FIN
! This subroutine deallocate work array that was used to determine the necessary size of 
! work arrays.
!******************************************************************************************
        SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDSYEV_FIN(SCAWORK)
           USE ML_FF_STRUCT, ONLY : SCALA_WORK
           IMPLICIT NONE
           TYPE(SCALA_WORK) :: SCAWORK
           SCAWORK%LWORK=MAX(1,INT(SCAWORK%WORK(1)))
           IF (ALLOCATED(SCAWORK%WORK)) THEN
             DEALLOCATE(SCAWORK%WORK)
           ENDIF
        END SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDSYEV_FIN

!******************************************************************************************
! SUBROUTINE SCALAPACK_WORK_ARRAY_PDSYEV_INIT
! This subroutine sets LWORK as -1 and allocate WORK for calculating the minimum size of 
! work arrays
!******************************************************************************************
        SUBROUTINE SCALAPACK_WORK_ARRAY_PDSYEV_INIT(SCAWORK)
           USE ML_FF_STRUCT, ONLY : SCALA_WORK
           IMPLICIT NONE
           TYPE(SCALA_WORK) :: SCAWORK
           IF (ALLOCATED(SCAWORK%WORK)) THEN
             DEALLOCATE(SCAWORK%WORK)
           ENDIF
           ALLOCATE(SCAWORK%WORK(SCAWORK%LWORK))
        END SUBROUTINE SCALAPACK_WORK_ARRAY_PDSYEV_INIT

!******************************************************************************************
! SUBROUTINE SCALAPACK_WORK_ARRAY_PDSYEV_FIN
! This subroutine deallocate work array that was used to determine the necessary size of 
! work arrays.
!******************************************************************************************
        SUBROUTINE SCALAPACK_WORK_ARRAY_PDSYEV_FIN(SCAWORK)
           USE ML_FF_STRUCT, ONLY : SCALA_WORK
           IMPLICIT NONE
           TYPE(SCALA_WORK) :: SCAWORK
           IF (ALLOCATED(SCAWORK%WORK)) THEN
             DEALLOCATE(SCAWORK%WORK)
           ENDIF
        END SUBROUTINE SCALAPACK_WORK_ARRAY_PDSYEV_FIN

!******************************************************************************************
! SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDGESVD_INIT
! This subroutine sets LWORK as -1 and allocate WORK for calculating the minimum size of 
! work arrays
!******************************************************************************************
        SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDGESVD_INIT(SCAWORK)
           USE ML_FF_STRUCT, ONLY : SCALA_WORK
           IMPLICIT NONE
           TYPE(SCALA_WORK) :: SCAWORK
           SCAWORK%LWORK=-1
           IF (ALLOCATED(SCAWORK%WORK)) THEN
             DEALLOCATE(SCAWORK%WORK)
           ENDIF
           ALLOCATE(SCAWORK%WORK(1))
        END SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDGESVD_INIT

!******************************************************************************************
! SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDGESVD_FIN
! This subroutine deallocate work array that was used to determine the necessary size of 
! work arrays.
!******************************************************************************************
        SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDGESVD_FIN(SCAWORK)
           USE ML_FF_STRUCT, ONLY : SCALA_WORK
           IMPLICIT NONE
           TYPE(SCALA_WORK) :: SCAWORK
           SCAWORK%LWORK=MAX(1,INT(SCAWORK%WORK(1)))
           IF (ALLOCATED(SCAWORK%WORK)) THEN
             DEALLOCATE(SCAWORK%WORK)
           ENDIF
        END SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDGESVD_FIN

!******************************************************************************************
! SUBROUTINE SCALAPACK_WORK_ARRAY_PDGESVD_INIT
! This subroutine sets LWORK as -1 and allocate WORK for calculating the minimum size of 
! work arrays
!******************************************************************************************
        SUBROUTINE SCALAPACK_WORK_ARRAY_PDGESVD_INIT(SCAWORK)
           USE ML_FF_STRUCT, ONLY : SCALA_WORK
           IMPLICIT NONE
           TYPE(SCALA_WORK) :: SCAWORK
           IF (ALLOCATED(SCAWORK%WORK)) THEN
             DEALLOCATE(SCAWORK%WORK)
           ENDIF
           ALLOCATE(SCAWORK%WORK(SCAWORK%LWORK))
        END SUBROUTINE SCALAPACK_WORK_ARRAY_PDGESVD_INIT

!******************************************************************************************
! SUBROUTINE SCALAPACK_WORK_ARRAY_PDSYEV_FIN
! This subroutine deallocate work array that was used to determine the necessary size of 
! work arrays.
!******************************************************************************************
        SUBROUTINE SCALAPACK_WORK_ARRAY_PDGESVD_FIN(SCAWORK)
           USE ML_FF_STRUCT, ONLY : SCALA_WORK
           IMPLICIT NONE
           TYPE(SCALA_WORK) :: SCAWORK
           IF (ALLOCATED(SCAWORK%WORK)) THEN
             DEALLOCATE(SCAWORK%WORK)
           ENDIF
        END SUBROUTINE SCALAPACK_WORK_ARRAY_PDGESVD_FIN

!******************************************************************************************
! SUBROUTINE SCALAPACK_WORK_ARRAY_PDGETRF_INIT
! This subroutine determine the size of IPIV necessary for PDGETRF
!******************************************************************************************
        SUBROUTINE SCALAPACK_WORK_ARRAY_PDGETRF_INIT (SCAWORK,MB,MROW,ME_ROW,NP_ROW)
           USE ML_FF_STRUCT, ONLY : SCALA_WORK
           IMPLICIT NONE
! Input variable
           TYPE(SCALA_WORK) :: SCAWORK
           INTEGER, INTENT(IN) :: MB
           INTEGER, INTENT(IN) :: ME_ROW
           INTEGER, INTENT(IN) :: MROW
           INTEGER, INTENT(IN) :: NP_ROW
           INTEGER              :: NUMROC
           EXTERNAL NUMROC
!          SCAWORK%LOC_M=MAX(1,NUMROC(MROW,MB,ME_ROW,0,NP_ROW))+MB
           SCAWORK%LOC_M=(((MROW+MB-1)/MB+NP_ROW-1)/NP_ROW)*MB+MB
           IF (ALLOCATED(SCAWORK%IPIV)) THEN
             DEALLOCATE(SCAWORK%IPIV)
           ENDIF
           ALLOCATE(SCAWORK%IPIV(1:SCAWORK%LOC_M))
        END SUBROUTINE SCALAPACK_WORK_ARRAY_PDGETRF_INIT

!******************************************************************************************
! SUBROUTINE SCALAPACK_WORK_ARRAY_PDGETRF_FIN
! This subroutine deallocates IPIV, IWORK and WORK that were used for determining
! necessary sizes of these arrays.
!******************************************************************************************
        SUBROUTINE SCALAPACK_WORK_ARRAY_PDGETRF_FIN(SCAWORK)
           USE ML_FF_STRUCT, ONLY : SCALA_WORK
           IMPLICIT NONE
           TYPE(SCALA_WORK) :: SCAWORK
           IF(ALLOCATED(SCAWORK%IPIV)) THEN
              DEALLOCATE(SCAWORK%IPIV)
           ENDIF
        END SUBROUTINE SCALAPACK_WORK_ARRAY_PDGETRF_FIN


!******************************************************************************************
! SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDGETRI_INIT
! This subroutine sets IPIV for calculations of inverse matrix
!******************************************************************************************
        SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDGETRI_INIT (SCAWORK,MB,MROW,ME_ROW,NP_ROW)
           USE ML_FF_STRUCT, ONLY : SCALA_WORK
           IMPLICIT NONE
! Input variable
           TYPE(SCALA_WORK) :: SCAWORK
           INTEGER, INTENT(IN) :: MB
           INTEGER, INTENT(IN) :: ME_ROW
           INTEGER, INTENT(IN) :: MROW
           INTEGER, INTENT(IN) :: NP_ROW
           SCAWORK%LIWORK=-1
           SCAWORK%LWORK=-1
           IF(ALLOCATED(SCAWORK%IWORK)) THEN
              DEALLOCATE(SCAWORK%IWORK)
           ENDIF
           ALLOCATE(SCAWORK%IWORK(1:1))
           IF(ALLOCATED(SCAWORK%WORK)) THEN
              DEALLOCATE(SCAWORK%WORK)
           ENDIF
           ALLOCATE(SCAWORK%WORK(1:1))
        END SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDGETRI_INIT

!******************************************************************************************
! SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDGETRI_FIN
! This subroutine deallocates IPIV, IWORK and WORK that were used for determining
! necessary sizes of these arrays.
!******************************************************************************************
        SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDGETRI_FIN(SCAWORK)
           USE ML_FF_STRUCT, ONLY : SCALA_WORK
           IMPLICIT NONE
           TYPE(SCALA_WORK) :: SCAWORK
           SCAWORK%LWORK=MAX(1,INT(SCAWORK%WORK(1)))
           SCAWORK%LIWORK=MAX(1,SCAWORK%IWORK(1))
           IF(ALLOCATED(SCAWORK%IWORK)) THEN
              DEALLOCATE(SCAWORK%IWORK)
           ENDIF
           IF(ALLOCATED(SCAWORK%WORK)) THEN
              DEALLOCATE(SCAWORK%WORK)
           ENDIF
        END SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDGETRI_FIN

!******************************************************************************************
! SUBROUTINE SCALAPACK_WORK_ARRAY_PDGETRI_INIT
! This subroutine allocates IPIV, IWORK and WORK for calculations of inverse matrix
!******************************************************************************************
        SUBROUTINE SCALAPACK_WORK_ARRAY_PDGETRI_INIT(SCAWORK)
           USE ML_FF_STRUCT, ONLY : SCALA_WORK
           IMPLICIT NONE
           TYPE(SCALA_WORK) :: SCAWORK
           IF(ALLOCATED(SCAWORK%IWORK)) THEN
              DEALLOCATE(SCAWORK%IWORK)
           ENDIF
           ALLOCATE(SCAWORK%IWORK(1:SCAWORK%LIWORK))
           IF(ALLOCATED(SCAWORK%WORK)) THEN
              DEALLOCATE(SCAWORK%WORK)
           ENDIF
           ALLOCATE(SCAWORK%WORK(1:SCAWORK%LWORK))
        END SUBROUTINE SCALAPACK_WORK_ARRAY_PDGETRI_INIT

!******************************************************************************************
! SUBROUTINE SCALAPACK_WORK_ARRAY_PDGETRI_FIN
! This subroutine deallocates SCAWORK%IPIV that was used in LU factorization
!******************************************************************************************
        SUBROUTINE SCALAPACK_WORK_ARRAY_PDGETRI_FIN(SCAWORK)
           USE ML_FF_STRUCT, ONLY : SCALA_WORK
           IMPLICIT NONE
           TYPE(SCALA_WORK) :: SCAWORK
           IF(ALLOCATED(SCAWORK%IPIV)) THEN
              DEALLOCATE(SCAWORK%IPIV)
           ENDIF
           IF(ALLOCATED(SCAWORK%IWORK)) THEN
              DEALLOCATE(SCAWORK%IWORK)
           ENDIF
           IF(ALLOCATED(SCAWORK%WORK)) THEN
              DEALLOCATE(SCAWORK%WORK)
           ENDIF
        END SUBROUTINE SCALAPACK_WORK_ARRAY_PDGETRI_FIN

!******************************************************************************************
! SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDGELS_INIT
! This subroutine sets LWORK as -1 and allocate WORK for calculating the minimum size of 
! work arrays
!******************************************************************************************
        SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDGELS_INIT(SCAWORK)
           USE ML_FF_STRUCT, ONLY : SCALA_WORK
           IMPLICIT NONE
           TYPE(SCALA_WORK) :: SCAWORK
           SCAWORK%LWORK=-1
           IF (ALLOCATED(SCAWORK%WORK)) THEN
             DEALLOCATE(SCAWORK%WORK)
           ENDIF
           ALLOCATE(SCAWORK%WORK(1))
        END SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDGELS_INIT

!******************************************************************************************
! SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDGELS_FIN
! This subroutine deallocate work array that was used to determine the necessary size of 
! work arrays.
!******************************************************************************************
        SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDGELS_FIN(SCAWORK)
           USE ML_FF_STRUCT, ONLY : SCALA_WORK
           IMPLICIT NONE
           TYPE(SCALA_WORK) :: SCAWORK
           SCAWORK%LWORK=MAX(1,INT(SCAWORK%WORK(1)))
           IF (ALLOCATED(SCAWORK%WORK)) THEN
             DEALLOCATE(SCAWORK%WORK)
           ENDIF
        END SUBROUTINE SCALAPACK_WORK_ARRAY_SIZE_PDGELS_FIN

!******************************************************************************************
! SUBROUTINE SCALAPACK_WORK_ARRAY_PDGELS_INIT
! This subroutine sets LWORK as -1 and allocate WORK for calculating the minimum size of 
! work arrays
!******************************************************************************************
        SUBROUTINE SCALAPACK_WORK_ARRAY_PDGELS_INIT(SCAWORK)
           USE ML_FF_STRUCT, ONLY : SCALA_WORK
           IMPLICIT NONE
           TYPE(SCALA_WORK) :: SCAWORK
           IF (ALLOCATED(SCAWORK%WORK)) THEN
             DEALLOCATE(SCAWORK%WORK)
           ENDIF
           ALLOCATE(SCAWORK%WORK(SCAWORK%LWORK))
        END SUBROUTINE SCALAPACK_WORK_ARRAY_PDGELS_INIT

!******************************************************************************************
! SUBROUTINE SCALAPACK_WORK_ARRAY_PDGELS_FIN
! This subroutine deallocate work array that was used to determine the
! necessary size of 
! work arrays.
!******************************************************************************************
        SUBROUTINE SCALAPACK_WORK_ARRAY_PDGELS_FIN(SCAWORK)
           USE ML_FF_STRUCT, ONLY : SCALA_WORK
           IMPLICIT NONE
           TYPE(SCALA_WORK) :: SCAWORK
           IF (ALLOCATED(SCAWORK%WORK)) THEN
             DEALLOCATE(SCAWORK%WORK)
           ENDIF
        END SUBROUTINE SCALAPACK_WORK_ARRAY_PDGELS_FIN

!******************************************************************************************
! SUBROUTINE SCALAPACK_REDUCE_TO_ONE_CORE
! This subroutine allocates the global array GENERAL_SINGLE and copies a SCALAPACK 
! distributed array CURRENT_ARRAY to it. GENERAL_SINGLE is only available on 1 core.
! The main purpose of this subroutine is for checking of the correctness of arrays
! distributed in parallel. 
!******************************************************************************************
        SUBROUTINE SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE,DESCRIPTOR,CURRENT_ARRAY,CONTEXT_CURRENT)
           USE ML_FF_STRUCT, ONLY : PARALLEL_SUPER, PARALLEL_DESCRIPTOR, SCALAPACK_GRID, ML_IO_WRITE
           IMPLICIT NONE
           TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
           TYPE (PARALLEL_DESCRIPTOR) :: DESCRIPTOR
           TYPE (SCALAPACK_GRID)      :: CONTEXT_CURRENT
           REAL(q)                    :: CURRENT_ARRAY(:,:) !(DESCRIPTOR%LOC_M,DESCRIPTOR%LOC_N)
           ! Local variables
           INTEGER                    :: I
           INTEGER                    :: IDIVISION
           INTEGER                    :: IERR
           INTEGER                    :: HELP_M
           INTEGER                    :: HELP_N
           INTEGER                    :: MSTART
           INTEGER                    :: NSTART
           INTEGER(KIND=8)            :: LIMIT_INTEGER
           INTEGER(KIND=8)            :: HELP_PRODUCT_DIMENSION
           PAR_SUP_HANDLE%DESC_GENERAL_SINGLE=DESCRIPTOR
           PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%MB=1
           PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%NB=1
           ! check here whether helping array was allocated from previous call,
           ! if yes it first needs to be deallocated
           IF (ALLOCATED(PAR_SUP_HANDLE%GENERAL_SINGLE)) THEN
              WRITE(*,*) "MLFF: Error in allocating helping array."
              WRITE(*,*) "Please deallocate array before, exiting."
              CALL MPI_BARRIER(PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
              CALL MPI_ABORT(PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
           ENDIF
           ! do allocation here on core 1, other nodes need allocation of (1,1)
           IF (PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_ROW.EQ.0.AND.PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_COL.EQ.0) THEN
              CALL SCALAPACK_SET_DESCRIPTOR (PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%DESC,PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%M, &
                PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%N,PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%MB,PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%NB, &
                PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%IRSRC,PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%ICSRC, &
                PAR_SUP_HANDLE%CONTEXT_SINGLE%ICTXT,PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_ROW,PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_COL,  &
                PAR_SUP_HANDLE%CONTEXT_SINGLE%NP_ROW,PAR_SUP_HANDLE%CONTEXT_SINGLE%NP_COL,PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%LLD, &
                PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%LOC_M,PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%LOC_N)
                ALLOCATE(PAR_SUP_HANDLE%GENERAL_SINGLE(PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%M,PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%N))
           ELSE
              PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%DESC=-1
              ALLOCATE(PAR_SUP_HANDLE%GENERAL_SINGLE(1,1))
           ENDIF
           ! We need to block over the number of elements if the product of both times 8 (for real type)
           ! is larger than the maximum allowed numberfor integer(kind=4), 
           ! since scalapack arrays use the smaller integer type
           LIMIT_INTEGER = INT(INT(2,KIND=8)**31,KIND=8)-1
           IDIVISION=1
           HELP_M=DESCRIPTOR%M
           HELP_N=DESCRIPTOR%N
           ! Determine the size of the blocking that needs to be done for PDGEMR2D
           DO 
              ! We need a type casting here to long integer
              HELP_PRODUCT_DIMENSION=8*INT(HELP_M,KIND=8)*INT(HELP_N,KIND=8)
              IF (HELP_PRODUCT_DIMENSION.GT.LIMIT_INTEGER) THEN
                 IDIVISION=IDIVISION+1
                 IF (DESCRIPTOR%M.GE.DESCRIPTOR%N) THEN
                    HELP_M=HELP_M/IDIVISION
                 ELSE
                    HELP_N=HELP_N/IDIVISION
                 ENDIF
              ELSE
                 EXIT
              ENDIF
           ENDDO
           ! Loop over the blocks
           DO I=1, IDIVISION
              IF (DESCRIPTOR%M.GE.DESCRIPTOR%N) THEN
                 MSTART=1+(I-1)*HELP_M
                 NSTART=1
              ELSE
                 MSTART=1
                 NSTART=1+(I-1)*HELP_N
              ENDIF
              CALL PDGEMR2D(HELP_M,HELP_N,CURRENT_ARRAY,MSTART,NSTART,DESCRIPTOR%DESC, &
                            PAR_SUP_HANDLE%GENERAL_SINGLE,MSTART,NSTART, &
                            PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%DESC,CONTEXT_CURRENT%ICTXT)
           ENDDO
           ! Do the remaining block size if there are any
           IF ((DESCRIPTOR%M-IDIVISION*HELP_M).GT.0 .OR. (DESCRIPTOR%N-IDIVISION*HELP_N).GT.0) THEN
              IF (DESCRIPTOR%M.GE.DESCRIPTOR%N) THEN
                 MSTART=IDIVISION*HELP_M+1
                 NSTART=1
                 HELP_M=DESCRIPTOR%M-IDIVISION*HELP_M
              ELSE
                 MSTART=1
                 NSTART=IDIVISION*HELP_N+1
                 HELP_N=DESCRIPTOR%N-IDIVISION*HELP_N
              ENDIF
              CALL PDGEMR2D(HELP_M,HELP_N,CURRENT_ARRAY,MSTART,NSTART,DESCRIPTOR%DESC, &
                            PAR_SUP_HANDLE%GENERAL_SINGLE,MSTART,NSTART, &
                            PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%DESC,CONTEXT_CURRENT%ICTXT)
           ENDIF
        END SUBROUTINE SCALAPACK_REDUCE_TO_ONE_CORE

!******************************************************************************************
! SUBROUTINE FINIT_SCALAPACK_REDUCE_TO_ONE_CORE
! This subroutine simply deallocates GENERAL_SINGLE. 
!******************************************************************************************
        SUBROUTINE FINIT_SCALAPACK_REDUCE_TO_ONE_CORE(PAR_SUP_HANDLE)
           USE ML_FF_STRUCT, ONLY : PARALLEL_SUPER
           IMPLICIT NONE
           TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
           DEALLOCATE(PAR_SUP_HANDLE%GENERAL_SINGLE)
        END SUBROUTINE FINIT_SCALAPACK_REDUCE_TO_ONE_CORE

!******************************************************************************************
! SUBROUTINE SCALAPACK_DISTRIBUTE_TO_ALL_CORES_ALLOCATE
! This subroutine allocates the global array GENERAL_SINGLE that will store
! all non-distributed data
!******************************************************************************************

        SUBROUTINE SCALAPACK_DISTRIBUTE_TO_ALL_CORES_ALLOCATE (PAR_SUP_HANDLE,DESCRIPTOR)
           USE ML_FF_STRUCT, ONLY : PARALLEL_SUPER, PARALLEL_DESCRIPTOR
           IMPLICIT NONE
           TYPE (PARALLEL_SUPER)      :: PAR_SUP_HANDLE
           TYPE (PARALLEL_DESCRIPTOR) :: DESCRIPTOR
           PAR_SUP_HANDLE%DESC_GENERAL_SINGLE=DESCRIPTOR
           PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%MB=1
           PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%NB=1
!           WRITE(*,*) CONTEXT_SINGLE%ME_ROW,CONTEXT_SINGLE%ME_COL,DESC_GENERAL_SINGLE%M,DESC_GENERAL_SINGLE%N
           IF (PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_ROW.EQ.0.AND.PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_COL.EQ.0) THEN
              CALL SCALAPACK_SET_DESCRIPTOR(PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%DESC,PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%M, &
                PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%N,PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%MB,PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%NB, &
                PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%IRSRC,PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%ICSRC, &
                PAR_SUP_HANDLE%CONTEXT_SINGLE%ICTXT,PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_ROW,PAR_SUP_HANDLE%CONTEXT_SINGLE%ME_COL,  &
                PAR_SUP_HANDLE%CONTEXT_SINGLE%NP_ROW,PAR_SUP_HANDLE%CONTEXT_SINGLE%NP_COL,PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%LLD, &
                PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%LOC_M,PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%LOC_N)
              IF (.NOT.ALLOCATED(PAR_SUP_HANDLE%GENERAL_SINGLE)) THEN
                 ALLOCATE(PAR_SUP_HANDLE%GENERAL_SINGLE(PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%M,PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%N))
              ENDIF
           ELSE
              PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%DESC=-1
              IF (.NOT.ALLOCATED(PAR_SUP_HANDLE%GENERAL_SINGLE)) THEN
                 ALLOCATE(PAR_SUP_HANDLE%GENERAL_SINGLE(1,1))
              ENDIF
           ENDIF
        END SUBROUTINE SCALAPACK_DISTRIBUTE_TO_ALL_CORES_ALLOCATE

!******************************************************************************************
! SUBROUTINE SCALAPACK_DISTRIBUTE_TO_ALL_CORES
! This subroutine distributes the global array GENERAL_SINGLE to all cores
!******************************************************************************************

        SUBROUTINE SCALAPACK_DISTRIBUTE_TO_ALL_CORES (PAR_SUP_HANDLE,DESCRIPTOR,CURRENT_ARRAY,CONTEXT_CURRENT)
           USE ML_FF_STRUCT, ONLY : PARALLEL_SUPER, PARALLEL_DESCRIPTOR, SCALAPACK_GRID
           IMPLICIT NONE
           TYPE (PARALLEL_SUPER)      :: PAR_SUP_HANDLE
           TYPE (PARALLEL_DESCRIPTOR) :: DESCRIPTOR
           TYPE (SCALAPACK_GRID)      :: CONTEXT_CURRENT
           REAL(q)                    :: CURRENT_ARRAY(:,:) !(DESCRIPTOR%LOC_M,DESCRIPTOR%LOC_N)
           ! Local variables
           INTEGER                    :: I
           INTEGER                    :: IDIVISION
           INTEGER                    :: HELP_M
           INTEGER                    :: HELP_N
           INTEGER                    :: MSTART
           INTEGER                    :: NSTART
           INTEGER(KIND=8)            :: LIMIT_INTEGER
           INTEGER(KIND=8)            :: HELP_PRODUCT_DIMENSION
           ! We need to block over the number of elements if the product of both times 8 (for real type)
           ! is larger than the maximum allowed numberfor integer(kind=4), 
           ! since scalapack arrays use the smaller integer type
           LIMIT_INTEGER = INT(INT(2,KIND=8)**31,KIND=8)-1
           IDIVISION=1
           HELP_M=PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%M
           HELP_N=PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%N
           ! Determine the size of the blocking that needs to be done for PDGEMR2D
           DO 
              ! We need a type casting here to long integer
              HELP_PRODUCT_DIMENSION=8*INT(HELP_M,KIND=8)*INT(HELP_N,KIND=8)
              IF (HELP_PRODUCT_DIMENSION.GT.LIMIT_INTEGER) THEN
                 IDIVISION=IDIVISION+1
                 IF (PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%M.GE.PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%N) THEN
                    HELP_M=HELP_M/IDIVISION
                 ELSE
                    HELP_N=HELP_N/IDIVISION
                 ENDIF
              ELSE
                 EXIT
              ENDIF
           ENDDO
           ! Loop over the blocks
           DO I=1, IDIVISION
              IF (PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%M.GE.PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%N) THEN
                 MSTART=1+(I-1)*HELP_M
                 NSTART=1
              ELSE
                 MSTART=1
                 NSTART=1+(I-1)*HELP_N
              ENDIF
              CALL PDGEMR2D(HELP_M,HELP_N, &
                   PAR_SUP_HANDLE%GENERAL_SINGLE,MSTART,NSTART,PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%DESC, &
                   CURRENT_ARRAY,MSTART,NSTART,DESCRIPTOR%DESC,CONTEXT_CURRENT%ICTXT)
           ENDDO
           ! Do the remaining block size if there are any
           IF ((PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%M-IDIVISION*HELP_M).GT.0 .OR. &
               (PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%N-IDIVISION*HELP_N).GT.0) THEN
              IF (PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%M.GE.PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%N) THEN
                 MSTART=IDIVISION*HELP_M+1
                 NSTART=1
                 HELP_M=DESCRIPTOR%M-IDIVISION*HELP_M
              ELSE
                 MSTART=1
                 NSTART=IDIVISION*HELP_N+1
                 HELP_N=DESCRIPTOR%N-IDIVISION*HELP_N
              ENDIF
              CALL PDGEMR2D(HELP_M,HELP_N, &
                   PAR_SUP_HANDLE%GENERAL_SINGLE,MSTART,NSTART,PAR_SUP_HANDLE%DESC_GENERAL_SINGLE%DESC, &
                   CURRENT_ARRAY,MSTART,NSTART,DESCRIPTOR%DESC,CONTEXT_CURRENT%ICTXT)
           ENDIF
        END SUBROUTINE SCALAPACK_DISTRIBUTE_TO_ALL_CORES

!******************************************************************************************
! SUBROUTINE FINIT_SCALAPACK_DISTRIBUTE_TO_ALL_CORES
! This subroutine simply deallocates GENERAL_SINGLE. 
!******************************************************************************************
        SUBROUTINE FINIT_SCALAPACK_DISTRIBUTE_TO_ALL_CORES(PAR_SUP_HANDLE)
           USE ML_FF_STRUCT, ONLY : PARALLEL_SUPER
           IMPLICIT NONE
           TYPE (PARALLEL_SUPER) :: PAR_SUP_HANDLE
           DEALLOCATE(PAR_SUP_HANDLE%GENERAL_SINGLE)
        END SUBROUTINE FINIT_SCALAPACK_DISTRIBUTE_TO_ALL_CORES
#endif

      END MODULE MPI_DATA
#endif
