#include "symbol.inc"
      MODULE ini
      USE prec
!**********************************************************************
!
!  this module implements a small timer utility to time
!  subroutine
!  (see START_TIMING)
!**********************************************************************

      ! this allows to set a maxmimum of 10 internested timers
      ! that should suffice for VASP
      INTEGER, PARAMETER, PRIVATE :: maxtimer=10

      INTEGER,SAVE :: used_timers=0
      CHARACTER (LEN=6), PRIVATE  :: timer_tag(maxtimer)
      REAL(q)      :: timer_vpu(maxtimer),timer_cpu(maxtimer)

      ! this allows to set a maxmimum of registered allocates
      INTEGER, PARAMETER, PRIVATE :: maxalloc=20
      INTEGER,SAVE                :: used_allocs=0
      CHARACTER (LEN=10), PRIVATE :: alloc_tag(maxalloc)
      REAL(q),PRIVATE             :: alloc_event(maxalloc)=0
      REAL(q),PRIVATE             :: alloc_total=0

      INTEGER, PRIVATE            :: used_allocs_hi
      REAL(q), PRIVATE            :: alloc_total_hi=0
      REAL(q), PRIVATE            :: alloc_event_hi(maxalloc)
      CHARACTER (LEN=10), PRIVATE :: alloc_tag_hi(maxalloc)

      INTEGER, PRIVATE ::  MINPGF,MAJPGF,ISWPS,IOOPS,IVCSW
      REAL(q), PRIVATE ::  UTIME,STIME,ETIME,RSIZM,AVSIZ,DAYTIM

      CONTAINS

!***********************************************************************
!
! timing routines
! START_TIMING(TAG) 
! registers a new timing routine with a specific name and 
! initialises the timer
!
!***********************************************************************

      SUBROUTINE START_TIMING(TAG)
      USE c2f_interface, ONLY : VTIME
      IMPLICIT NONE
      INTEGER ENTRY
      CHARACTER (LEN=*) :: TAG
      REAL(q) TV,TC

      CALL SEARCH_TIMING(TAG,ENTRY)
      IF (ENTRY==0) RETURN
    
      CALL VTIME(TV,TC)
      timer_vpu(ENTRY)=TV
      timer_cpu(ENTRY)=TC
      timer_tag(ENTRY)=TAG 

      END SUBROUTINE


      SUBROUTINE STOP_TIMING(TAG,IU,NAME,XMLTAG)
      USE vaspxml
      USE c2f_interface, ONLY : VTIME
      IMPLICIT NONE
      INTEGER ENTRY
      CHARACTER (LEN=*) :: TAG
      INTEGER :: IU
      CHARACTER (LEN=*),OPTIONAL :: XMLTAG
      CHARACTER (LEN=*),OPTIONAL :: NAME
      REAL(q) TV,TC

      CALL SEARCH_TIMING(TAG,ENTRY)
      IF (ENTRY==0) RETURN
    
      CALL VTIME(TV,TC)
      IF (PRESENT(NAME)) THEN
        IF (IU>=0) WRITE(IU,100) NAME,TV-timer_vpu(ENTRY),TC-timer_cpu(ENTRY)
      ELSE
        IF (IU>=0) WRITE(IU,100) TAG,TV-timer_vpu(ENTRY),TC-timer_cpu(ENTRY)
      ENDIF
      IF (PRESENT(XMLTAG)) &
           CALL XML_TIMING(TV-timer_vpu(ENTRY),TC-timer_cpu(ENTRY),name=XMLTAG)  


100   FORMAT(2X,A8,':  cpu time',F12.4,': real time',F12.4)

      
      timer_vpu(ENTRY)=TV
      timer_cpu(ENTRY)=TC
      timer_tag(ENTRY)=TAG 
      
      END SUBROUTINE

      SUBROUTINE SEPERATOR_TIMING(IU)

      IF (IU>0) WRITE(IU,100)
100   FORMAT(2X,'  --------------------------------------------')
      END SUBROUTINE SEPERATOR_TIMING


      SUBROUTINE SEARCH_TIMING(TAG,ENTRY)
      IMPLICIT NONE
      INTEGER ENTRY
      CHARACTER (LEN=*) :: TAG
      
      ! search for entry
      DO ENTRY=1,used_timers
         IF (timer_tag(ENTRY)==TAG) THEN
            RETURN
         ENDIF
      END DO

      IF (ENTRY>maxtimer) THEN
         ! no more entry available
         ENTRY=0
         WRITE(0,*) 'internal ERROR in SEARCH_TIMING: no more timing slot available'
      ELSE
         used_timers=used_timers+1
      ENDIF

      
      END SUBROUTINE SEARCH_TIMING

!***********************************************************************
!
! the following routines can be used to keep track of allocate
! and deallocate commands
! registered allocate calles also supply a tag
!
!***********************************************************************

      SUBROUTINE REGISTER_ALLOCATE(NALLOC, TAG)
      IMPLICIT NONE
      REAL(q) NALLOC
      CHARACTER (LEN=*), OPTIONAL :: TAG
      INTEGER ENTRY

      IF (PRESENT(TAG)) THEN
         CALL SEARCH_ALLOC(TAG, ENTRY)
         alloc_tag(ENTRY)  =TAG
         alloc_event(ENTRY)=alloc_event(ENTRY)+AINT(NALLOC/1000)
      END IF

      alloc_total=alloc_total+AINT(NALLOC/1000)

      IF (alloc_total>alloc_total_hi) THEN
         used_allocs_hi=used_allocs
         alloc_total_hi=alloc_total
         alloc_event_hi(1:used_allocs)=alloc_event(1:used_allocs)
         alloc_tag_hi(1:used_allocs)=alloc_tag(1:used_allocs)
      ENDIF

      END SUBROUTINE

      SUBROUTINE DEREGISTER_ALLOCATE(NALLOC, TAG)
      IMPLICIT NONE
      REAL(q) NALLOC
      CHARACTER (LEN=*), OPTIONAL :: TAG
      INTEGER ENTRY

      IF (PRESENT(TAG)) THEN
         CALL SEARCH_ALLOC(TAG, ENTRY)
         alloc_event(ENTRY)=alloc_event(ENTRY)-AINT(NALLOC/1000)
      END IF

      alloc_total=alloc_total-AINT(NALLOC/1000)

      END SUBROUTINE

      FUNCTION QUERRY_ALLOCATE()
      IMPLICIT NONE
      INTEGER QUERRY_ALLOCATE

      QUERRY_ALLOCATE=alloc_total
      END FUNCTION QUERRY_ALLOCATE


      SUBROUTINE DUMP_ALLOCATE(IU)
      IMPLICIT NONE
      INTEGER IU
      INTEGER ENTRY

      IF (IU>=0) THEN
      WRITE(IU,'(/1X,A,F10.0,A/A/)') 'total amount of memory used by VASP MPI-rank0',alloc_total,' kBytes', &
                                '======================================================================='

      DO ENTRY=1,used_allocs
         WRITE(IU,'(3X,A,A,F10.0,A)') alloc_tag(ENTRY),':  ',alloc_event(ENTRY),' kBytes'
      ENDDO
      WRITE(IU,*)
      ENDIF
      END SUBROUTINE DUMP_ALLOCATE


      SUBROUTINE DUMP_ALLOCATE_HI(IU)
      IMPLICIT NONE
      INTEGER IU
      INTEGER ENTRY

      IF (IU>=0) THEN
      WRITE(IU,'(/1X,A,F10.0,A/A/)') 'maximum amount of memory used by VASP MPI-rank0',alloc_total_hi,' kBytes', &
                                '======================================================================='

      DO ENTRY=1,used_allocs_hi
         WRITE(IU,'(3X,A,A,F10.0,A)') alloc_tag_hi(ENTRY),':  ',alloc_event_hi(ENTRY),' kBytes'
      ENDDO
      WRITE(IU,*)
      ENDIF
      END SUBROUTINE DUMP_ALLOCATE_HI


      SUBROUTINE DUMP_ALLOCATE_TAG(IU,TAG)
      IMPLICIT NONE
      INTEGER IU
      INTEGER ENTRY
      CHARACTER (LEN=*) :: TAG

      IF (IU>=0) THEN
      WRITE(IU,'(/1X,A,A,F10.0,A/A/)') 'memory high mark on MPI-rank0 inside ',TAG, alloc_total,' kBytes', &
                                '========================================================================'

      DO ENTRY=1,used_allocs
         WRITE(IU,'(3X,A,A,F10.0,A)') alloc_tag(ENTRY),':  ',alloc_event(ENTRY),' kBytes'
      ENDDO
      WRITE(IU,*)
      ENDIF
      END SUBROUTINE DUMP_ALLOCATE_TAG


      SUBROUTINE SEARCH_ALLOC(TAG,ENTRY)
      IMPLICIT NONE
      INTEGER ENTRY
      CHARACTER (LEN=*) :: TAG
      
      ! search for entry
      DO ENTRY=1,used_allocs
         IF (alloc_tag(ENTRY)==TAG) THEN
            RETURN
         ENDIF
      END DO

      IF (ENTRY>maxalloc) THEN
         ! no more entry available
         ENTRY=0
         WRITE(0,*) 'internal ERROR in SEARCH_ALLOC: no more registered allocation slots available'
      ELSE
         used_allocs=used_allocs+1
      ENDIF

      
      END SUBROUTINE SEARCH_ALLOC

      FUNCTION SEARCH_ALLOC_MEMORY(TAG)
      IMPLICIT NONE
      REAL(q)  SEARCH_ALLOC_MEMORY
      CHARACTER (LEN=*) :: TAG
      INTEGER ENTRY
      
      ! search for entry
      DO ENTRY=1,used_allocs
         IF (alloc_tag(ENTRY)==TAG) THEN
            SEARCH_ALLOC_MEMORY=alloc_event(ENTRY)
            RETURN
         ENDIF
      END DO

      SEARCH_ALLOC_MEMORY=0
      END FUNCTION SEARCH_ALLOC_MEMORY


!***********************************************************************
!
! dump some information on paging memory etc.
!
!***********************************************************************

      SUBROUTINE INIT_FINAL_TIMING()
      USE c2f_interface, ONLY : TIMING
      INTEGER IERR

      CALL TIMING(0,UTIME,STIME,ETIME,MINPGF,MAJPGF, &
     &            RSIZM,AVSIZ,ISWPS,IOOPS,IVCSW,IERR)
      IF (IERR/=0) ETIME=0._q

      END SUBROUTINE INIT_FINAL_TIMING

      SUBROUTINE DUMP_FINAL_TIMING(TIU6)
      USE c2f_interface, ONLY : TIMING
      INTEGER TIU6
      ! local
      INTEGER IERR
      INTEGER NODE_ME, IONODE

      CALL TIMING(0,UTIME,STIME,DAYTIM,MINPGF,MAJPGF, &
     &            RSIZM,AVSIZ,ISWPS,IOOPS,IVCSW,IERR)

      IF (TIU6>=0) THEN

      ETIME=DAYTIM-ETIME

      TOTTIM=UTIME+STIME
      WRITE(TIU6,*) ' '
      WRITE(TIU6,*) ' '
      WRITE(TIU6,'(A)') &
     &   ' General timing and accounting informations for this job:'
      WRITE(TIU6,'(A)') &
     &   ' ========================================================'
      WRITE(TIU6,*) ' '
      WRITE(TIU6,'(17X,A,F12.3)') ' Total CPU time used (sec): ',TOTTIM
      WRITE(TIU6,'(17X,A,F12.3)') '           User time (sec): ',UTIME
      WRITE(TIU6,'(17X,A,F12.3)') '         System time (sec): ',STIME
      WRITE(TIU6,'(17X,A,F12.3)') '        Elapsed time (sec): ',ETIME
      WRITE(TIU6,*) ' '
      WRITE(TIU6,'(17X,A,F12.0)') '  Maximum memory used (kb): ',RSIZM
      IF( AVSIZ == 0 ) THEN
         WRITE(TIU6,'(17X,A)') '  Average memory used (kb):          N/A'
      ELSE
         WRITE(TIU6,'(17X,A,F12.0)') '  Average memory used (kb): ',AVSIZ
      ENDIF
      WRITE(TIU6,*) ' '
      WRITE(TIU6,'(17X,A,I12)')   '         Minor page faults: ',MINPGF
      WRITE(TIU6,'(17X,A,I12)')   '         Major page faults: ',MAJPGF
      WRITE(TIU6,'(17X,A,I12)')   'Voluntary context switches: ',IVCSW
      ENDIF

      END SUBROUTINE



      END MODULE

!**************** SUBROUTINE SPLCOF, SPLCOF_N0 *************************
! RCS:  $Id: ini.F,v 1.3 2002/08/14 13:59:39 kresse Exp $
!
!  Subroutine for calculating spline-coefficients
!  using the routines of the book 'numerical  recipes'
!  on input P(1,N) must contain x-values
!           P(2,N) must contain function-values
!  YP is the first derivatives at the first point
!  if >= 10^30 natural boundary-contitions (y''=0) are used
!
!  for point N always natural boundary-conditions are used in
!  SPLCOF, whereas SPLCOF_N0 assume 0 derivative at N
!  SPLCOF_NDER allows to specify a boundary condition
!  at both end points
!
!***********************************************************************

      SUBROUTINE SPLCOF(P,N,NDIM,Y1P)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION P(NDIM,5)
!
!     determination of spline coefficients
!     ------------------------------------
!     f = ((d*dx+c)*dx+b)*dx+a
!         between adjacent x - values
!
!     result
!     P-ARRAY
!     P(I,1) = X(I)
!     P(I,2) = A(I) = F(I)
!     P(I,3) = B(I)
!     P(I,4) = C(I)
!     P(I,5) = D(I)
!
      IF (Y1P> .99E30_q) THEN
        P(1,4)=0.0_q
        P(1,3)=0.0_q
      ELSE
        P(1,4)=-.5_q
        P(1,3)=(3._q/(P(2,1)-P(1,1)))*((P(2,2)-P(1,2))/ &
     &             (P(2,1)-P(1,1))-Y1P)
      ENDIF

      DO 20 I=2,N-1
        S=(P(I,1)-P(I-1,1))/(P(I+1,1)-P(I-1,1))
        R=S*P(I-1,4)+2._q
        P(I,4)=(S-1._q)/R
        P(I,3)=(6*((P(I+1,2)-P(I,2))/(P(I+1,1)-P(I,1))- &
     &          (P(I,2)-P(I-1,2))/(P(I,1)-P(I-1,1)))/ &
     &          (P(I+1,1)-P(I-1,1))-S*P(I-1,3))/R
   20 ENDDO

      P(N,4)=0.0_q
      P(N,3)=0.0_q
!
      DO 30 I=N-1,1,-1
        P(I,4)=P(I,4)*P(I+1,4)+P(I,3)
  30  ENDDO
!
      DO 50 I=1,N-1
        S= P(I+1,1)-P(I,1)
        R=(P(I+1,4)-P(I,4))/6
        P(I,5)=R/S
        P(I,4)=P(I,4)/2.0_q
        P(I,3)=(P(I+1,2)-P(I,2))/S-(P(I,4)+R)*S
   50 ENDDO
      RETURN
      END SUBROUTINE



      SUBROUTINE SPLCOF_N0(P,N,NDIM,Y1P)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION P(NDIM,5)
!
!     determination of spline coefficients
!     ------------------------------------
!     f = ((d*dx+c)*dx+b)*dx+a
!         between adjacent x - values
!
!     result
!     P-ARRAY
!     P(I,1) = X(I)
!     P(I,2) = A(I) = F(I)
!     P(I,3) = B(I)
!     P(I,4) = C(I)
!     P(I,5) = D(I)
!
      IF (Y1P> .99E30_q) THEN
        P(1,4)=0.0_q
        P(1,3)=0.0_q
      ELSE
        P(1,4)=-.5_q
        P(1,3)=(3._q/(P(2,1)-P(1,1)))*((P(2,2)-P(1,2))/ &
     &             (P(2,1)-P(1,1))-Y1P)
      ENDIF

      DO 20 I=2,N-1
        S=(P(I,1)-P(I-1,1))/(P(I+1,1)-P(I-1,1))
        R=S*P(I-1,4)+2._q
        P(I,4)=(S-1._q)/R
        P(I,3)=(6*((P(I+1,2)-P(I,2))/(P(I+1,1)-P(I,1))- &
     &          (P(I,2)-P(I-1,2))/(P(I,1)-P(I-1,1)))/ &
     &          (P(I+1,1)-P(I-1,1))-S*P(I-1,3))/R
   20 ENDDO
      YNP=0
      IF (YNP> .99E30_q) THEN
        QN=0
        UN=0
      ELSE
        QN=0.5_q
        UN=(3._q/(P(N,1)-P(N-1,1)))*(YNP-(P(N,2)-P(N-1,2))/ &
     &             (P(N,1)-P(N-1,1)))
      ENDIF
      P(N,4)=(UN-QN*P(N-1,3))/(QN*P(N-1,4)+1.)
      P(N,3)=0  ! never used
!
      DO 30 I=N-1,1,-1
        P(I,4)=P(I,4)*P(I+1,4)+P(I,3)
  30  ENDDO
!
      DO 50 I=1,N-1
        S= P(I+1,1)-P(I,1)
        R=(P(I+1,4)-P(I,4))/6
        P(I,5)=R/S
        P(I,4)=P(I,4)/2.0_q
        P(I,3)=(P(I+1,2)-P(I,2))/S-(P(I,4)+R)*S
   50 ENDDO
      RETURN
      END SUBROUTINE


      SUBROUTINE SPLCOF_NDER(P,N,NDIM,Y1P,YNP)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION P(NDIM,5)
!
!     determination of spline coefficients
!     ------------------------------------
!     f = ((d*dx+c)*dx+b)*dx+a
!         between adjacent x - values
!
!     result
!     P-ARRAY
!     P(I,1) = X(I)
!     P(I,2) = A(I) = F(I)
!     P(I,3) = B(I)
!     P(I,4) = C(I)
!     P(I,5) = D(I)
!
      IF (Y1P> .99E30_q) THEN
        P(1,4)=0.0_q
        P(1,3)=0.0_q
      ELSE
        P(1,4)=-.5_q
        P(1,3)=(3._q/(P(2,1)-P(1,1)))*((P(2,2)-P(1,2))/ &
     &             (P(2,1)-P(1,1))-Y1P)
      ENDIF

      DO 20 I=2,N-1
        S=(P(I,1)-P(I-1,1))/(P(I+1,1)-P(I-1,1))
        R=S*P(I-1,4)+2._q
        P(I,4)=(S-1._q)/R
        P(I,3)=(6*((P(I+1,2)-P(I,2))/(P(I+1,1)-P(I,1))- &
     &          (P(I,2)-P(I-1,2))/(P(I,1)-P(I-1,1)))/ &
     &          (P(I+1,1)-P(I-1,1))-S*P(I-1,3))/R
   20 ENDDO
      IF (YNP> .99E30_q) THEN
        QN=0
        UN=0
      ELSE
        QN=0.5_q
        UN=(3._q/(P(N,1)-P(N-1,1)))*(YNP-(P(N,2)-P(N-1,2))/ &
     &             (P(N,1)-P(N-1,1)))
      ENDIF
      P(N,4)=(UN-QN*P(N-1,3))/(QN*P(N-1,4)+1.)
      P(N,3)=0  ! never used
!
      DO 30 I=N-1,1,-1
        P(I,4)=P(I,4)*P(I+1,4)+P(I,3)
  30  ENDDO
!
      DO 50 I=1,N-1
        S= P(I+1,1)-P(I,1)
        R=(P(I+1,4)-P(I,4))/6
        P(I,5)=R/S
        P(I,4)=P(I,4)/2.0_q
        P(I,3)=(P(I+1,2)-P(I,2))/S-(P(I,4)+R)*S
   50 ENDDO
      RETURN
      END SUBROUTINE

!
!  helper routine, which copies X and Y arrays to P
!  and than performes the fit on the array Y
!
      SUBROUTINE SPLCPY(X,Y,P,NAC,NDIM,Y1P)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION P(NDIM,5)
      DIMENSION X(NAC)
      DIMENSION Y(NAC)
      DO 100 N=1,NAC
        P(N,1)=X(N)
        P(N,2)=Y(N)
  100 CONTINUE
      CALL SPLCOF(P,NAC,NDIM,Y1P)
      RETURN
      END SUBROUTINE
!
!  helper routine, which evaluates the spline fit at a specific
!  position
!
      SUBROUTINE SPLVAL(X,F,FDER,P,NAC,NDIM)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION P(NDIM,5)
!  interval bisectioning
      I=1
      J=NAC
      IF (X   <P(I,1)) GO TO 60
      IF (X   <P(J,1)) GO TO 70
      K=J-1
      GOTO 90
   60 K=1
      GOTO 90
   70 K=(I+J)/2
      IF(I==K) GOTO 90
      IF (X   <P(K,1)) GO TO 80
      I=K
      GOTO 70
   80 J=K
      GOTO 70
!
   90 DX=X   -P(K,1)
      F   =((P(K,5)*DX+P(K,4))*DX+P(K,3))*DX+P(K,2)
      FDER=(3.0_q*P(K,5)*DX+2.0_q*P(K,4))*DX+P(K,3)
      END SUBROUTINE


!***********************************************************************
!
! system name date and time
!
!***********************************************************************

      SUBROUTINE MY_DATE_AND_TIME(IU6)
      USE prec
      USE vaspxml
      IMPLICIT NONE
      INTEGER IU6
      CHARACTER (8)  DATE
      CHARACTER (10) TIME
#ifndef HOST
#define HOST "unknown"
#endif

      CALL DATE_AND_TIME( DATE,TIME)
      IF (IU6>=0) &
      WRITE(IU6,"(' executed on ',A20,' date ',A4,'.',A2,'.',A2,'  ',A2,':',A2,':',A2 )") & 
           HOST,DATE(1:4),DATE(5:6),DATE(7:8),TIME(1:2),TIME(3:4),TIME(5:6)

      CALL XML_TAG_STRING("platform" , HOST)
      CALL XML_TAG_STRING("date" , DATE(1:4)//" "//DATE(5:6)//" "//DATE(7:8))
      CALL XML_TAG_STRING("time" , TIME(1:2)//":"//TIME(3:4)//":"//TIME(5:6))

      END SUBROUTINE


!***********************************************************************
!
! write out current memory requirements
!
!***********************************************************************

      SUBROUTINE MEMORY_CHECK(LOOP,STR)
      USE prec
      USE c2f_interface, ONLY : TIMING
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      CHARACTER (LEN=*) STR
      REAL(q) SUM

      CALL TIMING(0,UTIME,STIME,DAYTIM,MINPGF,MAJPGF, &
     &            RSIZM,AVSIZ,ISWPS,IOOPS,IVCSW,IERR)
      IF (IERR/=0) WRITE(*,*) 'WARNING main: call to TIMING failed.'
      IF (IERR/=0) WRITE(*,*) 'WARNING main: call to TIMING failed.'
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,'(A)') &
     &   ' General timing and accounting informations for this job:'
      WRITE(*,'(A)') &
     &   ' ========================================================'
      WRITE(*,*) ' '
      WRITE(*,'(17X,A,F12.0)') '  Maximum memory used (kb): ',RSIZM
      WRITE(*,'(17X,A,F12.0)') '  Average memory used (kb): ',AVSIZ
      WRITE(*,*) ' '
      WRITE(*,'(17X,A,I12)')   '         Minor page faults: ',MINPGF
      WRITE(*,'(17X,A,I12)')   '         Major page faults: ',MAJPGF
      WRITE(*,'(17X,A,I12)')   'Voluntary context switches: ',IVCSW

      END SUBROUTINE


!***********************************************************************
!
! DESCRIPTION:
!>This routine sets MAXMEM to a specific percentage of the available 
!>memory found on the node. The available memory is extracted from 
!>the file "/proc/meminfo" by searching for the line "FreeMem:". 
!>If this file is not present, nothing is done  
!
!***********************************************************************
      SUBROUTINE AUTOSET_AVAILABLE_MEMORY( MAXMEM, IO )
      USE tutor, ONLY: argument, VTUTOR, ISALERT, AUTOSET_MEM
      USE base
      USE main_mpi, only: COMM_WORLD, COMM_inter_node_world, COMM_intra_node_world
      INTEGER, INTENT(INOUT) :: MAXMEM 
      TYPE( in_struct )      :: IO
      !local 
      CHARACTER( LEN=64 )    :: MYNAME
      INTEGER                :: RESULTLEN, IERROR
      TYPE( argument )       :: ARG
      LOGICAL                :: LEXIST
      INTEGER                :: MEMSIZE
      INTEGER                :: NCPU = 1 
      INTEGER                :: NCPU_NODE = 1 
      INTEGER                :: NODE_ME = 1 
      INTEGER                :: IUNIT = 101
      REAL(q), PARAMETER     :: PERCENT=0.90_q ! how much percent of available
                                              ! memory should be used actually
#ifdef MPI
      NCPU = COMM_INTRA_NODE_WORLD%NCPU
      NCPU_NODE = COMM_INTRA_NODE_WORLD%NCPU
      NODE_ME = COMM_INTRA_NODE_WORLD%NODE_ME
      IUNIT=100 + COMM_WORLD%NODE_ME
#endif

      IERROR = 0 
      MEMSIZE = 0 
      ! only master rank on node reads /proc/meminfo
      IF ( NODE_ME == 1 ) THEN
         INQUIRE( FILE='/proc/meminfo', EXIST=LEXIST ) 
         IF ( .NOT. LEXIST ) THEN
            IERROR = 1
            GOTO 1000
         ENDIF
         
         ! open file that contains information about memory 
         OPEN( UNIT=IUNIT, FILE='/proc/meminfo', STATUS='OLD', IOSTAT=IERROR )
         IF ( IERROR /=0 ) THEN
            GOTO 1000 
         ENDIF
  
         ! find correct entry in file and proceed 
         CALL FIND_FREE_MEMORY_ENTRY( IUNIT, IERROR, MEMSIZE ) 

         ! 
         ! determine available RAM per rank on node
         ! 
         MEMSIZE = MEMSIZE / NCPU_NODE
      ENDIF 

                   
1000  CONTINUE 
      CALLMPI(M_sum_i( COMM_WORLD, IERROR , 1) )
      IF ( IERROR /=0 ) THEN
         ALLOCATE(ARG%IVAL(1))
         ARG%IVAL(1)=MAXMEM
         CALL VTUTOR%WRITE(ISALERT, AUTOSET_MEM, ARG)
         DEALLOCATE( ARG%IVAL )
         RETURN
      ENDIF
  
      ! communicate inside node 
      CALLMPI( M_sum_i( COMM_INTRA_NODE_WORLD, MEMSIZE, 1 ) ) 

      ! communicate only minimum storage available on all nodes 
      MEMSIZE = -MEMSIZE 
      CALLMPI( M_max_i( COMM_WORLD, MEMSIZE, 1 ) ) 
      MEMSIZE = -MEMSIZE 

      ! and use only a percentage of available memory by default
      MAXMEM = PERCENT*MEMSIZE/1024._q
  
      IF ( IO%IU0>=0 ) THEN
         WRITE(IO%IU0,100) PERCENT*MEMSIZE/1024._q/1024._q,  MAXMEM
      ENDIF 
      IF ( IO%IU6>=0 ) THEN
         WRITE(IO%IU6,100) PERCENT*MEMSIZE/1024._q/1024._q,  MAXMEM
      ENDIF 
  
100   FORMAT( ' available memory per node: ', F7.2, ' GB, setting MAXMEM to ',I7 )
  
      CONTAINS  
         
      !**************************************************************
      ! This routines looks in unit IUNIT for 
      ! specific entry that specifices the available memory in kB
      !**************************************************************
      SUBROUTINE FIND_FREE_MEMORY_ENTRY( IUNIT, IERROR, MEMSIZE ) 
         INTEGER, INTENT(IN)    :: IUNIT 
         INTEGER, INTENT(INOUT) :: IERROR ! breaking condition code
         INTEGER, INTENT(INOUT) :: MEMSIZE 
         ! local 
         CHARACTER( LEN=13 ), PARAMETER :: FINDME='MemAvailable:'
         CHARACTER( LEN=64)     :: LINE
         INTEGER                :: I_kB, I_MB
         CHARACTER( LEN=32)     :: CMEMSIZE
         INTEGER                :: I, J 
         INTEGER                :: IU, IMAX
         
         ! /proc/meminfo contains following entry
         ! 
         ! MemAvailable: An estimate of how much memory is available for starting new
         !     applications, without swapping. Calculated from MemFree,
         !     SReclaimable, the size of the file LRU lists, and the low
         !     watermarks in each zone.
         !     The estimate takes into account that the system needs some
         !     page cache to function well, and that not all reclaimable
         !     slab will be reclaimable, due to items being in use. The
         !     impact of those factors will vary from system to system.
  
         ! looking for a line with following format:
         !MemAvailable:   62159032 kB
         reading: DO 
            READ( IUNIT, '(A)', IOSTAT=IERROR) LINE
            IF ( IERROR /=0 ) THEN
               CLOSE( IUNIT )              
               EXIT reading 
            ENDIF
  
            ! line with "FreeAvailable:" has been found 
            IF ( INDEX( LINE, FINDME )>0 ) THEN
               
               ! find units of memory, take upper and lowercase into account
               I_kB = INDEX( LINE, 'kB')                   
               IF ( I_kB == 0 ) I_kB = INDEX( LINE, 'KB')                   
               I_MB = INDEX( LINE, 'mB')                   
               IF ( I_MB == 0 ) I_MB = INDEX( LINE, 'MB')                   
               
               ! calculations are done in kB, 
               ! if RAM is less than ~ 2 TB a 32 bit integer can store all digits
               IF ( I_kB == 0 .AND. I_MB == 0 ) THEN
                  IERROR = 4 
                  CLOSE( IUNIT )              
                  EXIT reading 
               ELSE IF ( I_kB /= 0 ) THEN
                  IMAX = I_kB
                  IU = 1 
               ELSE IF ( I_MB /= 0 ) THEN
                  IMAX = I_MB
                  IU = 1024
               ENDIF
  
               ! first collect digits
               CMEMSIZE='                                '
               J=0
               DO I = INDEX( LINE, ':' ) , IMAX 
                   ! if it is a numeric character add to 
                   IF ( &
                   LINE(I:I) == '1' .OR.  LINE(I:I) == '2' .OR. &
                   LINE(I:I) == '3' .OR.  LINE(I:I) == '4' .OR. & 
                   LINE(I:I) == '5' .OR.  LINE(I:I) == '6' .OR. & 
                   LINE(I:I) == '7' .OR.  LINE(I:I) == '8' .OR. & 
                   LINE(I:I) == '9' .OR.  LINE(I:I) == '0' ) THEN
                   J=J+1
                   CMEMSIZE(J:J) = LINE(I:I) 
                   ENDIF 
               ENDDO                  
               ! convert to integer  
               READ(CMEMSIZE(1:J), * ) MEMSIZE 
               MEMSIZE = MEMSIZE * IU 
  
               CLOSE( IUNIT )              
               EXIT reading 
            ENDIF
         ENDDO reading

      END SUBROUTINE FIND_FREE_MEMORY_ENTRY
  
      END SUBROUTINE AUTOSET_AVAILABLE_MEMORY

!***********************************************************************
!
! DESCRIPTION:
!>This routine extracts the currently used PSS of VASP 
!>The used memory is extracted from 
!>the file "/proc/self/status" by searching for the line "Pss:". 
!>If this file is not present, 0 is retured
!> Pss is the proportional set size, that is the number of pages in memrory
!> divided by the number of process that share this page
!> this is the proper way to estimate the correct memory usage of a process
!
!***********************************************************************
      SUBROUTINE GET_PSS( PSS, IU0, IU6, FSTRING, T )
      USE tutor, ONLY: argument, VTUTOR, ISALERT
      USE base
      USE main_mpi, only: COMM_WORLD, COMM_inter_node_world, COMM_intra_node_world
      USE c2f_interface, ONLY: getpid, getppid
      USE STRING, ONLY: STR
      REAL(q), INTENT(INOUT)       :: PSS
      INTEGER, INTENT(IN)          :: IU0
      INTEGER, INTENT(IN)          :: IU6
      CHARACTER(LEN=*), INTENT(IN) :: FSTRING
      INTEGER, INTENT(IN)          :: T
      !local 
      INTEGER                :: IERROR
      LOGICAL                :: LEXIST
      INTEGER                :: MEMSIZE
      INTEGER                :: NCPU = 1 
      INTEGER                :: NCPU_NODE = 1 
      INTEGER                :: NODE_ME = 1 
      INTEGER                :: IUNIT = 1001
      REAL(q)                :: RSIZM

      IF( T > 0  ) THEN
         CALL SLEEP( T ) 
      ENDIF

#ifdef MPI
      NCPU = COMM_INTRA_NODE_WORLD%NCPU
      NCPU_NODE = COMM_INTRA_NODE_WORLD%NCPU
      NODE_ME = COMM_INTRA_NODE_WORLD%NODE_ME
      IUNIT=1000 + COMM_WORLD%NODE_ME
#endif

      ! obtain peak memory usage 
      CALL PEAK_MEM( RSIZM ) 

      IERROR = 0 
      MEMSIZE = 0 

      ! only master rank on node reads /proc/self/smaps
      !INQUIRE( FILE='/proc/'//str(getpid())//'/status', EXIST=LEXIST ) 
      INQUIRE( FILE='/proc/self/smaps', EXIST=LEXIST ) 
      IF ( .NOT. LEXIST ) THEN
         IERROR = 1
         GOTO 1000
      ENDIF
      
      ! open file that contains information about memory 
      ! specifically, the proportional set size (PSS) is required 
      ! this is the amount of memory used by the process + shared memory / number of processes sharing it 
      !OPEN( NEWUNIT=IUNIT, FILE='/proc/'//str(getpid())//'/smaps', STATUS='OLD', IOSTAT=IERROR )
      OPEN( NEWUNIT=IUNIT, FILE='/proc/self/smaps', STATUS='OLD', IOSTAT=IERROR )
      IF ( IERROR /=0 ) THEN
         GOTO 1000 
      ENDIF
      ! find correct entry in file and proceed 
      CALL FIND_PSS_ENTRY( IUNIT, IERROR, MEMSIZE ) 
      ! 
      ! determine available RAM per rank on node
      ! 
      MEMSIZE = MEMSIZE / NCPU_NODE

1000  CONTINUE 
      CALLMPI(M_sum_i( COMM_WORLD, IERROR , 1) )
      IF ( IERROR /=0 ) THEN
         PSS=0
         IF ( IU0 >=0 ) WRITE( *, 1001 ) FSTRING, PSS, RSIZM
         IF ( IU6 >=0 ) WRITE( IU6, 1001 ) FSTRING, PSS, RSIZM 
         RETURN
      ENDIF
  
      ! communicate inside node 
      CALLMPI( M_sum_i( COMM_INTRA_NODE_WORLD, MEMSIZE, 1 ) ) 

      ! communicate maximum storage available on all nodes 
      CALLMPI( M_max_i( COMM_WORLD, MEMSIZE, 1 ) ) 
      
      ! scale to Mb 
      PSS = MEMSIZE /1024._q
      IF ( IU0 >=0 ) WRITE( *, 1001 ) FSTRING, PSS, RSIZM
      IF ( IU6 >=0 ) WRITE( IU6, 1001 ) FSTRING, PSS, RSIZM
1001  FORMAT( 1X, "**** Pss: ", A, " **** ", F14.2, F14.2 )
      CONTAINS  
         
      !**************************************************************
      ! This routines looks in unit IUNIT for 
      ! specific entry that specifices the used memory in Kb
      ! this is the number of pages divided by the number of processes
      ! sharing it
      !**************************************************************
      SUBROUTINE FIND_PSS_ENTRY( IUNIT, IERROR, MEMSIZE ) 
         INTEGER, INTENT(IN)    :: IUNIT 
         INTEGER, INTENT(INOUT) :: IERROR ! breaking condition code
         INTEGER, INTENT(INOUT) :: MEMSIZE 
         ! local 
         CHARACTER( LEN=6 ), PARAMETER :: FINDME='Pss:'
         CHARACTER( LEN=64)     :: LINE
         INTEGER                :: I_kB, I_MB
         CHARACTER( LEN=32)     :: CMEMSIZE
         INTEGER                :: MEMSIZE_ 
         INTEGER                :: I, J 
         INTEGER                :: IU, IMAX
         
         ! /proc/self/smaps
         ! 
         ! PSS: available in kernel newer than 2.6.14
         !      number of pages in memory divided by the number of 
         !      processes using it
         !      
         ! looking for a line with following format:
         !Pss:      59032 kB
         MEMSIZE = 0 
         reading: DO 
            READ( IUNIT, '(A)', IOSTAT=IERROR, END=999 ) LINE
            IF ( IERROR /=0 ) THEN
               CLOSE( IUNIT )              
               EXIT reading 
            ENDIF

            ! line with "Pss:" has been found 
            IF ( INDEX( LINE, FINDME )>0 ) THEN
               
               ! find units of memory, take upper and lowercase into account
               I_kB = INDEX( LINE, 'kB')                   
               IF ( I_kB == 0 ) I_kB = INDEX( LINE, 'KB')                   
               I_MB = INDEX( LINE, 'mB')                   
               IF ( I_MB == 0 ) I_MB = INDEX( LINE, 'MB')                   
               
               ! calculations are done in kB, 
               ! if RAM is less than ~ 2 TB a 32 bit integer can store all digits
               IF ( I_kB == 0 .AND. I_MB == 0 ) THEN
                  IERROR = 4 
                  CLOSE( IUNIT )              
                  EXIT reading 
               ELSE IF ( I_kB /= 0 ) THEN
                  IMAX = I_kB
                  IU = 1 
               ELSE IF ( I_MB /= 0 ) THEN
                  IMAX = I_MB
                  IU = 1024
               ENDIF
  
               ! first collect digits
               CMEMSIZE='                                '
               J=0
               DO I = INDEX( LINE, ':' ) , IMAX 
                   ! if it is a numeric character add to 
                   IF ( &
                   LINE(I:I) == '1' .OR.  LINE(I:I) == '2' .OR. &
                   LINE(I:I) == '3' .OR.  LINE(I:I) == '4' .OR. & 
                   LINE(I:I) == '5' .OR.  LINE(I:I) == '6' .OR. & 
                   LINE(I:I) == '7' .OR.  LINE(I:I) == '8' .OR. & 
                   LINE(I:I) == '9' .OR.  LINE(I:I) == '0' ) THEN
                   J=J+1
                   CMEMSIZE(J:J) = LINE(I:I) 
                   ENDIF 
               ENDDO                  
               ! convert to integer  
               READ(CMEMSIZE(1:J), * ) MEMSIZE_ 
               MEMSIZE = MEMSIZE + MEMSIZE_ * IU 
  
            ENDIF
         ENDDO reading

         RETURN

         ! special treatment for end of file
999      CONTINUE 
         IERROR = 0 
         CLOSE( IUNIT )              
      END SUBROUTINE FIND_PSS_ENTRY
  
      END SUBROUTINE GET_PSS

      !**************************************************************
      ! This routines obtains the peak memory usage  
      !**************************************************************
      SUBROUTINE PEAK_MEM( MYRSIZM )  
         USE main_mpi, only: COMM_WORLD
         USE base
         USE c2f_interface, ONLY: TIMING
         REAL(q), INTENT(INOUT)  :: MYRSIZM
         ! local
         INTEGER ::  MYMINPGF,MYMAJPGF,MYISWPS,MYIOOPS,MYIVCSW, MYIERR
         REAL(q) ::  MYMUTIME,MYSTIME,MYETIME,MYAVSIZ,MYDAYTIM
                  
          CALL TIMING(0,MYMUTIME,MYSTIME,MYETIME,MYMINPGF,MYMAJPGF, &
         &            MYRSIZM,MYAVSIZ,MYISWPS,MYIOOPS,MYIVCSW,MYIERR)
         MYRSIZM = MYRSIZM/1024._q

         CALLMPI( M_max_d( COMM_WORLD, MYRSIZM, 1 ) ) 

      END SUBROUTINE PEAK_MEM 
