
      MODULE internals
        USE prec
        USE mymath
        USE poscar
        USE base
        IMPLICIT NONE 
        CONTAINS
        
        SUBROUTINE update_where(T_INFO,DPOS,lattmat,ICOORDINATES)
          TYPE(type_info) :: T_INFO
          TYPE(coordstructure) :: ICOORDINATES
          TYPE(gadget_io) :: g_io
          REAL(q) :: DPOS(3,T_INFO%NIONS)        !c positions in direct
          REAL(q) :: lattmat(3,3)  
          INTEGER :: i,j,k
          
          DO i=1,ICOORDINATES%NUMINTERNALS
            CALL where_shortest(DPOS,lattmat,ICOORDINATES%COORDSTRUCT(i))
          ENDDO

        END SUBROUTINE update_where

        SUBROUTINE constraints_reader(T_INFO,g_io,DPOS,lattmat,lattinv,CCOORDINATES,IO)
          USE string, ONLY: str
          USE tutor, ONLY: vtutor
          TYPE(type_info) :: T_INFO
          TYPE(coordstructure) :: CCOORDINATES
          TYPE(gadget_io) :: g_io
          TYPE (in_struct) :: IO
          REAL(q) :: DPOS(3,T_INFO%NIONS)        !c positions in direct
          REAL(q) :: lattmat(3,3),lattinv(3,3)          !c direct and reciprocal lattice
          CHARACTER (LEN=2) :: tag
          INTEGER :: counter,primcounter,what(4),ios,info,dummy,stat
          INTEGER :: complexcounter
          REAL(q) :: dummyp,dummyq,dummyr,dummys,dummyt
          REAL(q) :: auxcoef(4)
          INTEGER :: i,j,k,ii,jj,kk
          CHARACTER (LEN=20000) :: line 
          LOGICAL :: LIS=.FALSE.,LIZ=.FALSE.
          INTEGER :: IIS=0,IIZ=0
          REAL(q) :: criteria(3)
          INTEGER :: crit(3)
          REAL(q), PARAMETER :: RMAX=30. !c cutoff radius for PBC adapted coordination numbers (WP)

          counter=0
          primcounter=0;complexcounter=0
          OPEN(UNIT=g_io%CONSTRAINTS,FILE='ICONST',STATUS='UNKNOWN')

          CCOORDINATES%iconst0=0
          CCOORDINATES%iconst1=0
          CCOORDINATES%iconst2=0
          CCOORDINATES%iconst3=0
          CCOORDINATES%iconst4=0
          CCOORDINATES%iconst5=0
          CCOORDINATES%iconst6=0
          CCOORDINATES%iconst7=0
          CCOORDINATES%iconst8=0
          CCOORDINATES%iconst9=0
          CCOORDINATES%LIRC=.FALSE.
          CCOORDINATES%IRCDIM=0
          !CCOORDINATES%CNEXP(1)=9.;CCOORDINATES%CNEXP(2)=14.
          !CCOORDINATES%CNEXP(1)=8.;CCOORDINATES%CNEXP(2)=20.

          !c determine the the number of images of simulation cell
          !c in each direction needed for a given
          !c cutoff radius
          criteria(1)=RMAX*SUM(lattinv(:,1)**2)**0.5
          criteria(2)=RMAX*SUM(lattinv(:,2)**2)**0.5
          criteria(3)=RMAX*SUM(lattinv(:,3)**2)**0.5
          !crit=int(criteria)+1
          crit=int(criteria)
   
          DO
            what=0;auxcoef=0._q
            READ(g_io%CONSTRAINTS,FMT='(A20000)',IOSTAT=ios) line
            IF (ios<0) EXIT
            !write(*,*) LEN_TRIM(line),'LEN_TRIM(line)'
            IF (LEN_TRIM(line)>0) THEN
            line=TRIM(ADJUSTL(line))
            tag=line(1:2)
            line=line(3:)

            !READ(g_io%CONSTRAINTS,FMT='(A2)',IOSTAT=ios,ADVANCE='NO') tag
            !IF (ios<0) EXIT
            counter=counter+1

            IF (counter>SIZE(CCOORDINATES%COORDSTRUCT)) THEN
              CCOORDINATES%COORDSTRUCT=>coordreallocate(CCOORDINATES%COORDSTRUCT,counter+10)
            ENDIF
            nullify(CCOORDINATES%COORDSTRUCT(counter)%COEFS)
            
            stat=100
            SELECT CASE(tag)
              CASE('X ')  !c fractional 1
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) what(1),stat
                READ(line,FMT=*,IOSTAT=ios) what(1),stat
                CALL sanity_test(T_INFO%NIONS,what,1,IO,counter)
                CCOORDINATES%COORDSTRUCT(counter)%TAG='X '
                primcounter=primcounter+1
              CASE('Y ')  !c fractional 2
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) what(1),stat
                READ(line,FMT=*,IOSTAT=ios) what(1),stat
                CALL sanity_test(T_INFO%NIONS,what,1,IO,counter)
                CCOORDINATES%COORDSTRUCT(counter)%TAG='Y '
                primcounter=primcounter+1
              CASE('Z ')  !c fractional 3
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) what(1),stat
                READ(line,FMT=*,IOSTAT=ios) what(1),stat
                CALL sanity_test(T_INFO%NIONS,what,1,IO,counter)
                CCOORDINATES%COORDSTRUCT(counter)%TAG='Z '
                primcounter=primcounter+1
              CASE('cX')  !c cartesian x
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) what(1),stat
                READ(line,FMT=*,IOSTAT=ios) what(1),stat
                CALL sanity_test(T_INFO%NIONS,what,1,IO,counter)
                CCOORDINATES%COORDSTRUCT(counter)%TAG='cX'
                primcounter=primcounter+1
              CASE('cY')  !c cartesian y
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) what(1),stat
                READ(line,FMT=*,IOSTAT=ios) what(1),stat
                CALL sanity_test(T_INFO%NIONS,what,1,IO,counter)
                CCOORDINATES%COORDSTRUCT(counter)%TAG='cY'
                primcounter=primcounter+1
              CASE('cZ')  !c cartesian z
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) what(1),stat
                READ(line,FMT=*,IOSTAT=ios) what(1),stat
                CALL sanity_test(T_INFO%NIONS,what,1,IO,counter)
                CCOORDINATES%COORDSTRUCT(counter)%TAG='cZ'
                primcounter=primcounter+1  
              CASE('R ')  !c interatomic distance
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) what(1),what(2),stat
                READ(line,FMT=*,IOSTAT=ios) what(1),what(2),stat
                CALL sanity_test(T_INFO%NIONS,what,2,IO,counter)
                CCOORDINATES%COORDSTRUCT(counter)%TAG='R '
                primcounter=primcounter+1
              CASE('Q ')  !c interatomic distance^2
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) what(1),what(2),stat
                READ(line,FMT=*,IOSTAT=ios) what(1),what(2),stat
                CALL sanity_test(T_INFO%NIONS,what,2,IO,counter)
                CCOORDINATES%COORDSTRUCT(counter)%TAG='Q '
                primcounter=primcounter+1
              CASE('A ')  !c bond angle
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) what(1),what(2),what(3),stat
                READ(line,FMT=*,IOSTAT=ios) what(1),what(2),what(3),stat
                CALL sanity_test(T_INFO%NIONS,what,3,IO,counter)
                CCOORDINATES%COORDSTRUCT(counter)%TAG='A '
                primcounter=primcounter+1
              CASE('T ')  !c torsion
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) what(1),what(2),what(3),what(4),stat
                READ(line,FMT=*,IOSTAT=ios) what(1),what(2),what(3),what(4),stat
                CALL sanity_test(T_INFO%NIONS,what,4,IO,counter)
                CCOORDINATES%COORDSTRUCT(counter)%TAG='T '
                primcounter=primcounter+1
              CASE('M ')  !c atom-midpoint distance
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) what(1),what(2),what(3),stat
                READ(line,FMT=*,IOSTAT=ios) what(1),what(2),what(3),stat
                CALL sanity_test(T_INFO%NIONS,what,3,IO,counter)
                CCOORDINATES%COORDSTRUCT(counter)%TAG='M '
                primcounter=primcounter+1
              CASE('B ')  !c midpoint-midpoint distance
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) what(1),what(2),what(3),what(4),stat
                READ(line,FMT=*,IOSTAT=ios) what(1),what(2),what(3),what(4),stat
                CALL sanity_test(T_INFO%NIONS,what,4,IO,counter)
                CCOORDINATES%COORDSTRUCT(counter)%TAG='B '
                primcounter=primcounter+1
              CASE('P ')  !c ratio between two distances
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) what(1),what(2),what(3),what(4),stat
                READ(line,FMT=*,IOSTAT=ios) what(1),what(2),what(3),what(4),stat
                CALL sanity_test(T_INFO%NIONS,what,4,IO,counter)
                CCOORDINATES%COORDSTRUCT(counter)%TAG='P '
                primcounter=primcounter+1
              CASE('W ')  !c a single contribution to coordination number (D)
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) what(1),what(2),what(3),what(4),stat
                READ(line,FMT=*,IOSTAT=ios) what(1),what(2),auxcoef(1),auxcoef(2),auxcoef(3),stat
                !c make sure that R0 in (R/R0)**m is greater than 0
                IF (auxcoef(1).LE.0) THEN
                  CALL vtutor%error("Error reading ICONST (item " // str(primcounter+1) // "): invalid &
                     &definition of coordinate")
                ENDIF
                !c avoid division by 0 in [1-(R/R0)**m]/[1-(R/R0)**n] 
                IF (auxcoef(3).EQ.0) THEN
                  CALL vtutor%error("Error reading ICONST (item " // str(primcounter+1) // "): invalid &
                     &definition of coordinate")
                ENDIF
                CALL sanity_test(T_INFO%NIONS,what,2,IO,counter)
                CCOORDINATES%COORDSTRUCT(counter)%TAG='W '
                primcounter=primcounter+1  
              CASE('WP')  !c a single contribution to coordination number (D), PBC adapted
                READ(line,FMT=*,IOSTAT=ios) what(1),what(2),auxcoef(1),auxcoef(2),auxcoef(3),stat
                !c make sure that R0 in (R/R0)**m is greater than 0
                IF (auxcoef(1).LE.0) THEN
                  IF (IO%IU0>=0) write(IO%IU0,*) 'Error reading ICONST (item',primcounter+1,'): invalid definition of coordinate'
                  STOP
                ENDIF
                !c avoid division by 0 in [1-(R/R0)**m]/[1-(R/R0)**n]
                IF (auxcoef(3).EQ.0) THEN
                  IF (IO%IU0>=0) write(IO%IU0,*) 'Error reading ICONST (item',primcounter+1,'): invalid definition of coordinate'
                  STOP
                ENDIF
                CALL sanity_test(T_INFO%NIONS,what,2,IO,counter)
                CCOORDINATES%COORDSTRUCT(counter)%TAG='WP'
                CCOORDINATES%COORDSTRUCT(counter)%TRANS=crit
                primcounter=primcounter+1
              CASE('C ')  !c norm of vector of coordinates
                CCOORDINATES%COORDSTRUCT(counter)%TAG='C '
                ALLOCATE(CCOORDINATES%COORDSTRUCT(counter)%COEFS(primcounter))
                CCOORDINATES%COORDSTRUCT(counter)%COEFS=0
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) CCOORDINATES%COORDSTRUCT(counter)%COEFS(:),stat
                READ(line,FMT=*,IOSTAT=ios) CCOORDINATES%COORDSTRUCT(counter)%COEFS(:),stat
                CCOORDINATES%COORDSTRUCT(counter)%STATUS=stat
                complexcounter=complexcounter+1
              CASE('S ')  !c sum of coordinates
                CCOORDINATES%COORDSTRUCT(counter)%TAG='S '
                ALLOCATE(CCOORDINATES%COORDSTRUCT(counter)%COEFS(primcounter))
                CCOORDINATES%COORDSTRUCT(counter)%COEFS=0
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) CCOORDINATES%COORDSTRUCT(counter)%COEFS(:),stat
                READ(line,FMT=*,IOSTAT=ios) CCOORDINATES%COORDSTRUCT(counter)%COEFS(:),stat
                CCOORDINATES%COORDSTRUCT(counter)%STATUS=stat
                complexcounter=complexcounter+1
              CASE('ST')  !c sum of coordinates
                CCOORDINATES%COORDSTRUCT(counter)%TAG='ST'
                ALLOCATE(CCOORDINATES%COORDSTRUCT(counter)%COEFS(primcounter))
                CCOORDINATES%COORDSTRUCT(counter)%COEFS=0
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) CCOORDINATES%COORDSTRUCT(counter)%COEFS(:),stat
                READ(line,FMT=*,IOSTAT=ios) CCOORDINATES%COORDSTRUCT(counter)%COEFS(:),stat
                CCOORDINATES%COORDSTRUCT(counter)%STATUS=stat
                complexcounter=complexcounter+1
              CASE('D ')  !c coordination number
                CCOORDINATES%COORDSTRUCT(counter)%TAG='D '
                ALLOCATE(CCOORDINATES%COORDSTRUCT(counter)%COEFS(primcounter))
                CCOORDINATES%COORDSTRUCT(counter)%COEFS=0
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) CCOORDINATES%COORDSTRUCT(counter)%COEFS(:),stat
                READ(line,FMT=*,IOSTAT=ios) CCOORDINATES%COORDSTRUCT(counter)%COEFS(:),stat
                CCOORDINATES%COORDSTRUCT(counter)%STATUS=stat
                complexcounter=complexcounter+1
              CASE('IS')  !c position on IRC (defined in file IRCPOINTS)
                !!IF (.NOT. LIS) THEN
                  CCOORDINATES%COORDSTRUCT(counter)%TAG='IS'
                  ALLOCATE(CCOORDINATES%COORDSTRUCT(counter)%COEFS(primcounter))
                  CCOORDINATES%COORDSTRUCT(counter)%COEFS=0
                  READ(line,FMT=*,IOSTAT=ios) CCOORDINATES%COORDSTRUCT(counter)%COEFS(:),stat
                  CCOORDINATES%COORDSTRUCT(counter)%STATUS=stat
                  !c make sure that the coefficients for IS and IZ are synchronized 
                  IF (LIZ) CCOORDINATES%COORDSTRUCT(counter)%COEFS(:)=CCOORDINATES%COORDSTRUCT(iiz)%COEFS(:)
                  IIS=counter
                  LIS=.TRUE.
                  CCOORDINATES%LIRC=.TRUE.
                  IF (.NOT. LIZ) THEN
                    DO j=1,primcounter 
                      IF  (CCOORDINATES%COORDSTRUCT(counter)%COEFS(j)>0.) CCOORDINATES%IRCDIM(1)=CCOORDINATES%IRCDIM(1)+1 
                    ENDDO
                  ENDIF
                !!ELSE
                !!  CALL vtutor%error("Error ICONST (item " // str(counter) // "): multiple IS are not &
                !!  &   supported")
                !!ENDIF
              CASE('IZ')  !c distance from IRC (defined in file IRCPOINTS)
                !!IF (.NOT. LIZ) THEN
                  CCOORDINATES%COORDSTRUCT(counter)%TAG='IZ'
                  ALLOCATE(CCOORDINATES%COORDSTRUCT(counter)%COEFS(primcounter))
                  CCOORDINATES%COORDSTRUCT(counter)%COEFS=0
                  READ(line,FMT=*,IOSTAT=ios) CCOORDINATES%COORDSTRUCT(counter)%COEFS(:),stat
                  CCOORDINATES%COORDSTRUCT(counter)%STATUS=stat
                  !c make sure that the coefficients for IS and IZ are synchronized 
                  IF (LIS) CCOORDINATES%COORDSTRUCT(counter)%COEFS(:)=CCOORDINATES%COORDSTRUCT(iis)%COEFS(:)
                  IIZ=counter
                  LIZ=.TRUE.
                  CCOORDINATES%LIRC=.TRUE.
                  IF (.NOT. LIS) THEN
                    DO j=1,primcounter 
                      IF  (CCOORDINATES%COORDSTRUCT(counter)%COEFS(j)>0.) CCOORDINATES%IRCDIM(1)=CCOORDINATES%IRCDIM(1)+1 
                    ENDDO
                  ENDIF
                !!ELSE
                !!  CALL vtutor%error("Error ICONST (item " // str(counter) // "): multiple IZ are not &
                !!  &  supported")
                !!ENDIF
              CASE('LR')  !c interatomic distance
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) what(1),stat
                READ(line,FMT=*,IOSTAT=ios) what(1),stat
                !CALL sanity_test(T_INFO%NIONS,what,1,IO,counter)
                CCOORDINATES%COORDSTRUCT(counter)%TAG='LR'
                primcounter=primcounter+1
              CASE('LA')  !c interatomic distance
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) what(1),what(2),stat
                READ(line,FMT=*,IOSTAT=ios) what(1),what(2),stat
                !CALL sanity_test(T_INFO%NIONS,what,2,IO,counter)
                CCOORDINATES%COORDSTRUCT(counter)%TAG='LA'
                primcounter=primcounter+1
              CASE('LV')  !c interatomic distance
                !READ(g_io%CONSTRAINTS,FMT=*,IOSTAT=ios) stat
                READ(line,FMT=*,IOSTAT=ios) stat
                CCOORDINATES%COORDSTRUCT(counter)%TAG='LV'
                primcounter=primcounter+1
              CASE DEFAULT
                CALL vtutor%error("Error reading ICONST (item " // str(counter) // "): unsupported &
                   &coordinate type")
            END SELECT

            !c terminate if STATUS was not defined for a coordinate in ICONST
            IF ( .NOT.  (stat==0 .OR. stat==1 .OR. stat==3 .OR. stat==4 .OR. stat==5 .OR. stat==7 .OR. stat==8)) THEN
              CALL vtutor%error("Error reading ICONST (item " // str(counter) // "): missing or &
                  &invalid STATUS")
            ENDIF

            complexcounter=counter-primcounter
            IF (counter==primcounter) THEN
              CCOORDINATES%COORDSTRUCT(counter)%WHAT=what
              CCOORDINATES%COORDSTRUCT(counter)%AUXCOEF=auxcoef
              CALL where_shortest(DPOS,lattmat,CCOORDINATES%COORDSTRUCT(counter))
              CCOORDINATES%COORDSTRUCT(counter)%STATUS=stat
              SELECT CASE(tag)
                CASE('X ')
                  CCOORDINATES%COORDSTRUCT(counter)%VALUE=DPOS(1,CCOORDINATES%COORDSTRUCT(counter)%WHAT(1))
                CASE('Y ')
                  CCOORDINATES%COORDSTRUCT(counter)%VALUE=DPOS(2,CCOORDINATES%COORDSTRUCT(counter)%WHAT(1))
                CASE('Z ')
                  CCOORDINATES%COORDSTRUCT(counter)%VALUE=DPOS(3,CCOORDINATES%COORDSTRUCT(counter)%WHAT(1))
                CASE('cX')
                  CCOORDINATES%COORDSTRUCT(counter)%VALUE=&
                  & DPOS(1,CCOORDINATES%COORDSTRUCT(counter)%WHAT(1))*lattmat(1,1)+&
                  & DPOS(2,CCOORDINATES%COORDSTRUCT(counter)%WHAT(1))*lattmat(2,1)+&
                  & DPOS(3,CCOORDINATES%COORDSTRUCT(counter)%WHAT(1))*lattmat(3,1)
                CASE('cY')
                  CCOORDINATES%COORDSTRUCT(counter)%VALUE=&
                  & DPOS(1,CCOORDINATES%COORDSTRUCT(counter)%WHAT(1))*lattmat(1,2)+&
                  & DPOS(2,CCOORDINATES%COORDSTRUCT(counter)%WHAT(1))*lattmat(2,2)+&
                  & DPOS(3,CCOORDINATES%COORDSTRUCT(counter)%WHAT(1))*lattmat(3,2)
                CASE('cZ')
                  CCOORDINATES%COORDSTRUCT(counter)%VALUE=&
                  & DPOS(1,CCOORDINATES%COORDSTRUCT(counter)%WHAT(1))*lattmat(1,3)+&
                  & DPOS(2,CCOORDINATES%COORDSTRUCT(counter)%WHAT(1))*lattmat(2,3)+&
                  & DPOS(3,CCOORDINATES%COORDSTRUCT(counter)%WHAT(1))*lattmat(3,3)
                CASE('R ')
                  CCOORDINATES%COORDSTRUCT(counter)%VALUE=CAL_DISTANCE(T_INFO%NIONS,TRANSPOSE(lattmat),&
                    & DPOS,what(1),what(2),CCOORDINATES%COORDSTRUCT(counter)%WHERE(1,:),&
                  & CCOORDINATES%COORDSTRUCT(counter)%WHERE(2,:))
                CASE('Q ')
                  CCOORDINATES%COORDSTRUCT(counter)%VALUE=CAL_DISTANCE(T_INFO%NIONS,TRANSPOSE(lattmat),&
                    & DPOS,what(1),what(2),CCOORDINATES%COORDSTRUCT(counter)%WHERE(1,:),&
                  & CCOORDINATES%COORDSTRUCT(counter)%WHERE(2,:))
                  CCOORDINATES%COORDSTRUCT(counter)%VALUE=CCOORDINATES%COORDSTRUCT(counter)%VALUE**2
                CASE('A ')
                  CCOORDINATES%COORDSTRUCT(counter)%VALUE=CAL_ANGLE(T_INFO%NIONS,transpose(lattmat),&
                  & DPOS,what(1),what(2),what(3),CCOORDINATES%COORDSTRUCT(counter)%WHERE(1,:),&
                   &   CCOORDINATES%COORDSTRUCT(counter)%WHERE(2,:),CCOORDINATES%COORDSTRUCT(counter)%WHERE(3,:),info)
                CASE('T ')
                  CCOORDINATES%COORDSTRUCT(counter)%VALUE=CAL_TORSION(T_INFO%NIONS,transpose(lattmat),&
                  & DPOS,what(1),what(2),what(3),what(4),CCOORDINATES%COORDSTRUCT(counter)%WHERE(1,:),&
                  & CCOORDINATES%COORDSTRUCT(counter)%WHERE(2,:),CCOORDINATES%COORDSTRUCT(counter)%WHERE(3,:),&
                  & CCOORDINATES%COORDSTRUCT(counter)%WHERE(4,:),info)
                CASE('M ')
                  CCOORDINATES%COORDSTRUCT(counter)%VALUE=CAL_MIDDLE(T_INFO%NIONS,transpose(lattmat),&
                  & DPOS,what(1),what(2),what(3),CCOORDINATES%COORDSTRUCT(counter)%WHERE(1,:),&
                   &   CCOORDINATES%COORDSTRUCT(counter)%WHERE(2,:),CCOORDINATES%COORDSTRUCT(counter)%WHERE(3,:))
                CASE('B ')
                  CCOORDINATES%COORDSTRUCT(counter)%VALUE=CAL_MIDDLE2(T_INFO%NIONS,transpose(lattmat),&
                  & DPOS,what(1),what(2),what(3),what(4),CCOORDINATES%COORDSTRUCT(counter)%WHERE(1,:),&
                   &   CCOORDINATES%COORDSTRUCT(counter)%WHERE(2,:),CCOORDINATES%COORDSTRUCT(counter)%WHERE(3,:),&
                  & CCOORDINATES%COORDSTRUCT(counter)%WHERE(4,:))
                CASE('P ')
                  dummyp=CAL_DISTANCE(T_INFO%NIONS,TRANSPOSE(lattmat),&
                    & DPOS,what(1),what(2),CCOORDINATES%COORDSTRUCT(counter)%WHERE(1,:),&
                  & CCOORDINATES%COORDSTRUCT(counter)%WHERE(2,:))
                  dummyq=CAL_DISTANCE(T_INFO%NIONS,TRANSPOSE(lattmat),&
                    & DPOS,what(3),what(4),CCOORDINATES%COORDSTRUCT(counter)%WHERE(3,:),&
                  & CCOORDINATES%COORDSTRUCT(counter)%WHERE(4,:))
                  IF (dummyq>0.) THEN
                    CCOORDINATES%COORDSTRUCT(counter)%VALUE=dummyp/dummyq
                  ELSE
                    CALL vtutor%error("Error ICONST (item " // str(counter) // "): incorrect value of P")
                  ENDIF
                !CASE('E ')
                !  CCOORDINATES%COORDSTRUCT(counter)%VALUE=0.
                !CASE('F ')
                !  CCOORDINATES%COORDSTRUCT(counter)%VALUE=0.
                !CASE('G ')
                !  CCOORDINATES%COORDSTRUCT(counter)%VALUE=0.
                CASE('W ')
                  dummyp=CAL_DISTANCE(T_INFO%NIONS,TRANSPOSE(lattmat),&
                  & DPOS,what(1),what(2),CCOORDINATES%COORDSTRUCT(counter)%WHERE(1,:),&
                  & CCOORDINATES%COORDSTRUCT(counter)%WHERE(2,:))
                  dummyp=dummyp/CCOORDINATES%COORDSTRUCT(counter)%AUXCOEF(1)
                  dummyq=1.-dummyp**CCOORDINATES%COORDSTRUCT(counter)%AUXCOEF(2)
                  CCOORDINATES%COORDSTRUCT(counter)%VALUE=dummyq/(1.-dummyp**CCOORDINATES%COORDSTRUCT(counter)%AUXCOEF(3))
                CASE('WP')
                  CCOORDINATES%COORDSTRUCT(counter)%VALUE=0.
                  DO ii=-CCOORDINATES%COORDSTRUCT(counter)%TRANS(1),CCOORDINATES%COORDSTRUCT(counter)%TRANS(1)
                    DO jj=-CCOORDINATES%COORDSTRUCT(counter)%TRANS(2),CCOORDINATES%COORDSTRUCT(counter)%TRANS(2)
                      DO kk=-CCOORDINATES%COORDSTRUCT(counter)%TRANS(3),CCOORDINATES%COORDSTRUCT(counter)%TRANS(3)
                        dummyp=CAL_DISTANCE_trans(T_INFO%NIONS,TRANSPOSE(lattmat),&
                        & DPOS,what(1),what(2),CCOORDINATES%COORDSTRUCT(counter)%WHERE(1,:),&
                        & CCOORDINATES%COORDSTRUCT(counter)%WHERE(2,:), &
                        &  (/ii,jj,kk/))
                         dummyp=dummyp/CCOORDINATES%COORDSTRUCT(counter)%AUXCOEF(1)
                         dummyq=1.-dummyp**CCOORDINATES%COORDSTRUCT(counter)%AUXCOEF(2)
                         CCOORDINATES%COORDSTRUCT(counter)%VALUE= CCOORDINATES%COORDSTRUCT(counter)%VALUE + &
                         &              dummyq/(1.-dummyp**CCOORDINATES%COORDSTRUCT(counter)%AUXCOEF(3))
                      ENDDO
                    ENDDO
                  ENDDO
                CASE('LR')
                  CCOORDINATES%COORDSTRUCT(counter)%VALUE=CAL_LR(TRANSPOSE(lattmat),what(1))
                CASE('LA')
                  CCOORDINATES%COORDSTRUCT(counter)%VALUE=CAL_LA(TRANSPOSE(lattmat),what(1),what(2))
                CASE('LV')
                  CCOORDINATES%COORDSTRUCT(counter)%VALUE=CAL_LV(TRANSPOSE(lattmat))
              END SELECT
            ENDIF
            ENDIF
          ENDDO

          CLOSE(g_io%CONSTRAINTS)
          CCOORDINATES%NUMINTERNALS=counter
          CCOORDINATES%NUMPRIMITIVES=primcounter
          CCOORDINATES%COORDSTRUCT=>coordreallocate(CCOORDINATES%COORDSTRUCT,counter)

          !c read in IRC points if needed (i.e. if IS or/and IZ coordinate types are defined in ICONST)
          IF (CCOORDINATES%LIRC) THEN
            OPEN(UNIT=g_io%IRCPOINTS,FILE='IRCCAR',STATUS='UNKNOWN')
            READ(g_io%IRCPOINTS,FMT=*,IOSTAT=ios) CCOORDINATES%IRCDIM(2)
            !write(*,*) "CCOORDINATES%IRCDIM(2)",CCOORDINATES%IRCDIM(1),CCOORDINATES%IRCDIM(2)
            IF (ios<0) THEN
              CALL vtutor%error("constraints_reader: Error reading IRCCAR")
            ENDIF
            ALLOCATE(CCOORDINATES%IRCPOINTS(CCOORDINATES%IRCDIM(1),CCOORDINATES%IRCDIM(2)))
                
            DO i=1,CCOORDINATES%IRCDIM(2)
              READ(g_io%IRCPOINTS,FMT=*,IOSTAT=ios) CCOORDINATES%IRCPOINTS(:,i)
              IF (ios<0) THEN
                CALL vtutor%error("constraints_reader: Error reading IRCCAR, missing line")
              ENDIF
            ENDDO

            !DO i=1,CCOORDINATES%IRCDIM(2)
            !     DO j=1,CCOORDINATES%IRCDIM(1)
            !       IF (IO%IU0>0) write(*,*) "irc",i,j,CCOORDINATES%IRCPOINTS(j,i)
            !     ENDDO
            !ENDDO
            CLOSE(g_io%IRCPOINTS)
          ENDIF

          IF (counter/=primcounter) THEN
            DO i=1,primcounter
              CCOORDINATES%COORDSTRUCT(i)%STATUS=1
            ENDDO
            DO i=primcounter+1,counter
              CCOORDINATES%COORDSTRUCT(i)%VALUE=0._q
              SELECT CASE(CCOORDINATES%COORDSTRUCT(i)%TAG)
                CASE('C ')
                  DO j=1,SIZE(CCOORDINATES%COORDSTRUCT(i)%COEFS)
                    CCOORDINATES%COORDSTRUCT(i)%VALUE=CCOORDINATES%COORDSTRUCT(i)%VALUE+&
                    (CCOORDINATES%COORDSTRUCT(j)%VALUE*CCOORDINATES%COORDSTRUCT(i)%COEFS(j))**2
                  ENDDO
                  CCOORDINATES%COORDSTRUCT(i)%VALUE=CCOORDINATES%COORDSTRUCT(i)%VALUE**0.5
                  !CCOORDINATES%COORDSTRUCT(i)%STATUS=0
                CASE('S ')
                  dummyp=0.
                  dummyq=0.
                  dummyr=0.
                  dummys=0.
                  dummyt=0.
                  DO j=1,SIZE(CCOORDINATES%COORDSTRUCT(i)%COEFS)
                    !c torsions require a special care
                    IF (CCOORDINATES%COORDSTRUCT(j)%TAG=='T ') THEN
                      IF (CCOORDINATES%COORDSTRUCT(i)%COEFS(j)>0.) THEN
                        CALL Min_Image_Cyclic(dummyp,CCOORDINATES%COORDSTRUCT(j)%VALUE)
                        dummyq=dummyq+CCOORDINATES%COORDSTRUCT(j)%VALUE*CCOORDINATES%COORDSTRUCT(i)%COEFS(j)
                        dummyp=CCOORDINATES%COORDSTRUCT(j)%VALUE
                        dummys=dummys+CCOORDINATES%COORDSTRUCT(i)%COEFS(j)
                      ELSE IF (CCOORDINATES%COORDSTRUCT(i)%COEFS(j)<0.) THEN
                        CALL Min_Image_Cyclic(dummyp,CCOORDINATES%COORDSTRUCT(j)%VALUE)
                        dummyr=dummyr+CCOORDINATES%COORDSTRUCT(j)%VALUE*CCOORDINATES%COORDSTRUCT(i)%COEFS(j)
                        dummyp=CCOORDINATES%COORDSTRUCT(j)%VALUE
                        dummyt=dummyt+(CCOORDINATES%COORDSTRUCT(i)%COEFS(j))
                      ENDIF
                    ELSE
                      CCOORDINATES%COORDSTRUCT(i)%VALUE=CCOORDINATES%COORDSTRUCT(i)%VALUE+&
                      CCOORDINATES%COORDSTRUCT(j)%VALUE*CCOORDINATES%COORDSTRUCT(i)%COEFS(j)
                    ENDIF
                  ENDDO
                  !c apply use minimum image along with CBC for torsions
                  IF (dummys > 0.) THEN
                    dummyq=dummyq/dummys
                    CALL Min_Image_Cyclic(0._q,dummyq)
                  ENDIF
                  IF (dummyt < 0.) THEN
                    dummyr=dummyr/dummyt
                    CALL Min_Image_Cyclic(dummyq,dummyr)
                  ENDIF
                  dummyq=dummyq*dummys
                  dummyr=dummyr*dummyt
                  CCOORDINATES%COORDSTRUCT(i)%VALUE=CCOORDINATES%COORDSTRUCT(i)%VALUE+dummyq+dummyr
                  !CCOORDINATES%COORDSTRUCT(i)%STATUS=0
               CASE('ST')
                  dummyp=0.
                  dummyq=0.
                  dummyr=0.
                  dummys=0.
                  dummyt=0.
                  DO j=1,SIZE(CCOORDINATES%COORDSTRUCT(i)%COEFS)
                    !c torsions require a special care
                    IF (CCOORDINATES%COORDSTRUCT(j)%TAG=='T ') THEN
                      IF (CCOORDINATES%COORDSTRUCT(i)%COEFS(j)>0.) THEN
                        CALL Min_Image_Cyclic(dummyp,CCOORDINATES%COORDSTRUCT(j)%VALUE)
                        dummyq=dummyq+CCOORDINATES%COORDSTRUCT(j)%VALUE*CCOORDINATES%COORDSTRUCT(i)%COEFS(j)
                        dummyp=CCOORDINATES%COORDSTRUCT(j)%VALUE
                        dummys=dummys+CCOORDINATES%COORDSTRUCT(i)%COEFS(j)
                      ELSE IF (CCOORDINATES%COORDSTRUCT(i)%COEFS(j)<0.) THEN
                        CALL Min_Image_Cyclic(dummyp,CCOORDINATES%COORDSTRUCT(j)%VALUE)
                        dummyr=dummyr+CCOORDINATES%COORDSTRUCT(j)%VALUE*CCOORDINATES%COORDSTRUCT(i)%COEFS(j)
                        dummyp=CCOORDINATES%COORDSTRUCT(j)%VALUE
                        dummyt=dummyt+(CCOORDINATES%COORDSTRUCT(i)%COEFS(j))
                      ENDIF
                    ENDIF
                  ENDDO
                  !c apply use minimum image along with CBC for torsions
                  IF (dummys > 0.) THEN
                    dummyq=dummyq/dummys
                    CALL Min_Image_Cyclic(0._q,dummyq)
                  ENDIF
                  IF (dummyt < 0.) THEN
                    dummyr=dummyr/dummyt
                    CALL Min_Image_Cyclic(dummyq,dummyr)
                  ENDIF
                  dummyq=dummyq*dummys
                  dummyr=dummyr*dummyt
                  CCOORDINATES%COORDSTRUCT(i)%VALUE=dummyq+dummyr
                  !CCOORDINATES%COORDSTRUCT(i)%STATUS=0
                CASE('IS')
                  dummyp=0.
                  dummyq=0.
                  DO k=1,CCOORDINATES%IRCDIM(2)
                    dummyr=0.
                    DO j=1,CCOORDINATES%NUMPRIMITIVES
                      IF (CCOORDINATES%COORDSTRUCT(i)%COEFS(j)>0.) THEN
                        dummyr=dummyr+CCOORDINATES%COORDSTRUCT(i)%COEFS(j) * &
                        &   (CCOORDINATES%COORDSTRUCT(j)%VALUE-CCOORDINATES%IRCPOINTS(j,k))**2
                      ENDIF
                    ENDDO
                    dummyr=exp(-dummyr)
                    dummyp=dummyp+(k-1)*dummyr
                    dummyq=dummyq+dummyr
                  ENDDO  
                  IF (dummyq>0.) CCOORDINATES%COORDSTRUCT(i)%VALUE=dummyp/(dummyq*(CCOORDINATES%IRCDIM(2)-1))
                CASE('IZ')
                  dummyq=0.
                  dummyp=0.
                  DO k=1,CCOORDINATES%IRCDIM(2)
                    dummyr=0.
                    DO j=1,CCOORDINATES%NUMPRIMITIVES
                      IF (CCOORDINATES%COORDSTRUCT(i)%COEFS(j)>0.) THEN
                        dummyp=CCOORDINATES%COORDSTRUCT(i)%COEFS(j)
                        dummyr=dummyr+CCOORDINATES%COORDSTRUCT(i)%COEFS(j) * &
                        &   (CCOORDINATES%COORDSTRUCT(j)%VALUE-CCOORDINATES%IRCPOINTS(j,k))**2
                      ENDIF
                    ENDDO
                    dummyr=exp(-dummyr)
                    dummyq=dummyq+dummyr
                  ENDDO  
                  IF (dummyq>0.) CCOORDINATES%COORDSTRUCT(i)%VALUE=-LOG(dummyq)/dummyp 
                CASE('D ')
                  DO j=1,SIZE(CCOORDINATES%COORDSTRUCT(i)%COEFS)
                    IF (ABS(CCOORDINATES%COORDSTRUCT(i)%COEFS(j))>1e-4) THEN
                      dummyq=CCOORDINATES%COORDSTRUCT(j)%VALUE/ABS(CCOORDINATES%COORDSTRUCT(i)%COEFS(j))
                      IF (ABS(dummyq-1._q)<1e-4) THEN
                        dummyq=1._q+1e-4
                      ENDIF
                      IF (CCOORDINATES%COORDSTRUCT(i)%COEFS(j)>0._q) THEN
                        CCOORDINATES%COORDSTRUCT(i)%VALUE=CCOORDINATES%COORDSTRUCT(i)%VALUE+&
                        & (1.0-(dummyq)**CCOORDINATES%CNEXP(1))/(1.0-(dummyq)**CCOORDINATES%CNEXP(2))
                      ELSE
                        CCOORDINATES%COORDSTRUCT(i)%VALUE=CCOORDINATES%COORDSTRUCT(i)%VALUE-&
                        & (1.0-(dummyq)**CCOORDINATES%CNEXP(1))/(1.0-(dummyq)**CCOORDINATES%CNEXP(2))
                      ENDIF
                    ENDIF
                  ENDDO
                  
              END SELECT
            ENDDO
          ENDIF

          CCOORDINATES%COORDSTRUCT(:)%DVALUE=0._q

          !c identify type of coordinates defined in ICONST
          DO i=1,CCOORDINATES%NUMINTERNALS
            SELECT CASE(CCOORDINATES%COORDSTRUCT(i)%STATUS)
              CASE(0) !c number of hard constraints
                CCOORDINATES%iconst0=CCOORDINATES%iconst0+1
              CASE(1) !c basis for complex coordinates
                CCOORDINATES%iconst1=CCOORDINATES%iconst1+1
              CASE(2)
                CCOORDINATES%iconst2=CCOORDINATES%iconst2+1
              CASE(3) !c harmonmic forcefiled HESSEMAT
                CCOORDINATES%iconst3=CCOORDINATES%iconst3+1
              CASE(4) !c fermi-shape potential
                CCOORDINATES%iconst4=CCOORDINATES%iconst4+1
              CASE(5) !c number of collective variables for mMD
                CCOORDINATES%iconst5=CCOORDINATES%iconst5+1
              CASE(6)
                CCOORDINATES%iconst6=CCOORDINATES%iconst6+1
              CASE(7) !c number of monitored coordinates
                CCOORDINATES%iconst7=CCOORDINATES%iconst7+1
              CASE(8) !c harmonic springs will act on these internal coordinates
                CCOORDINATES%iconst8=CCOORDINATES%iconst8+1
              !CASE(9) !c number of dummy atoms conections
              !  CCOORDINATES%iconst9=CCOORDINATES%iconst9+1 
            END SELECT  
          ENDDO

        END SUBROUTINE constraints_reader

        SUBROUTINE sanity_test(NIONS,WHAT,N,IO,num)
         !c check that the numbers used to define a constraint are in the interval <1,NIONS>
         !c quit with error message otherwise
          USE string, ONLY: str
          USE tutor, ONLY: vtutor
          TYPE (in_struct) :: IO
          INTEGER :: NIONS,N
          INTEGER :: WHAT(4)
          INTEGER :: i,num
          LOGICAL :: ldum

          ldum=.FALSE.
          DO i=1,N
            IF ((WHAT(i)>NIONS) .OR. (WHAT(i)<1)) THEN
              ldum=.TRUE.
            END IF
          ENDDO
          
          IF (ldum) THEN
            CALL vtutor%error("Error reading ICONST (item " // str(num) // "): invalid definition of &
               &coordinate")
          ENDIF
        END SUBROUTINE sanity_test

        !FUNCTION CAL_BOND(A,POINT1,POINT2) RESULT(BOND)
        !  REAL(q) :: POINT1(3),POINT2(3),VECT(3)
        !  REAL(q) :: A(3,3)
        !  REAL(q) :: BOND
        !  INTEGER :: i
 
        !  VECT=POINT1-POINT2
        !  CALL min_imageV(3,VECT)
        !  VECT=MATMUL(A,VECT)
        !  BOND=(SUM(VECT*VECT))**0.5
        !END FUNCTION

          FUNCTION CAL_DISTANCE(NIONS,A,DX,p1,p2,w1,w2) RESULT(dist)
            INTEGER :: NIONS
            REAL(q) :: A(3,3)
            REAL(q) :: DX(3,NIONS)
            REAL(q) :: vect(3)
            INTEGER :: p1,p2
            INTEGER :: w1(3),w2(3)
            REAL(q) :: dist

            vect=DX(:,p2)+w2 - (DX(:,p1)+w1)
            CALL min_imageV(3,vect)
            vect=MATMUL(A,vect)
            dist=VECTORSIZE(3,vect)
          END FUNCTION

          FUNCTION CAL_DISTANCE_trans(NIONS,A,DX,p1,p2,w1,w2,trans) RESULT(dist)
            INTEGER :: NIONS
            REAL(q) :: A(3,3)
            REAL(q) :: DX(3,NIONS)
            REAL(q) :: vect(3)
            INTEGER :: p1,p2
            INTEGER :: w1(3),w2(3),trans(3)
            REAL(q) :: dist

            vect=DX(:,p2)+w2 - (DX(:,p1)+w1)
            CALL min_imageV(3,VECT)
            vect=MATMUL(A,vect+trans)
            dist=VECTORSIZE(3,vect)
          END FUNCTION

          FUNCTION CAL_LR(A,p1) RESULT(dist)
            REAL(q) :: A(3,3)
            REAL(q) :: vect(3)
            INTEGER :: p1
            REAL(q) :: dist

            vect=A(:,p1)
            dist=VECTORSIZE(3,vect)
          END FUNCTION

         FUNCTION CAL_LA(A,p1,p2) RESULT(angle)
            REAL(q) :: A(3,3)
            REAL(q) :: vect1(3),vect2(3)
            INTEGER :: p1,p2
            REAL(q) :: angle

            vect1=A(:,p1)
            vect2=A(:,p2)
            angle=sum(vect1*vect2)/VECTORSIZE(3,vect1)/VECTORSIZE(3,vect2)
            angle=ACOS(angle)
          END FUNCTION

          FUNCTION CAL_LV(A) RESULT(dist)
            REAL(q) :: A(3,3)
            REAL(q) :: vect(3)
            REAL(q) :: dist

            vect=CROSSPROD(3,A(:,1),A(:,2))
            dist=sum(vect*A(:,3))
          END FUNCTION

          FUNCTION CAL_ANGLE(NIONS,A,DX,p1,p2,p3,w1,w2,w3,info) RESULT(angle)
            INTEGER :: NIONS
            REAL(q) :: A(3,3)
            REAL(q) :: DX(3,NIONS)
            REAL(q) :: vect1(3),vect2(3)
            INTEGER :: p1,p2,p3
            INTEGER :: w1(3),w2(3),w3(3)
            INTEGER :: info
            REAL(q) :: angle

            info=0
            vect1=DX(:,p2)+w2 - (DX(:,p1)+w1)
            CALL min_imageV(3,vect1)
            vect2=DX(:,p2)+w2 - (DX(:,p3)+w3)
            CALL min_imageV(3,vect2)
            vect1=MATMUL(A,vect1)
            vect2=MATMUL(A,vect2)
            angle=DOT_PRODUCT(vect1,vect2)
            angle=angle/VECTORSIZE(3,vect1)/VECTORSIZE(3,vect2)
            IF (angle>1) angle=1._q
            IF (angle<-1) angle=-1._q
            angle=ACOS(angle)
            IF (SIN(angle)>0.3) info=1
          END FUNCTION

          FUNCTION CAL_TORSION(NIONS,A,DX,p1,p2,p3,p4,w1,w2,w3,w4,info) RESULT(torsion)
            INTEGER :: NIONS,p1,p2,p3,p4
            INTEGER :: w1(3),w2(3),w3(3),w4(3)
            INTEGER :: info
            REAL(q) :: A(3,3)
            REAL(q) :: DX(3,NIONS)
            REAL(q) :: vect1(3),vect2(3),vect3(3),vect4(3)
            REAL(q) :: cross1(3),cross2(3)
            REAL(q) :: torsion,tst
           
            info=0
            vect1=DX(:,p1)+w1 - (DX(:,p2)+w2)
            CALL min_imageV(3,vect1)
            vect2=DX(:,p2)+w2 - (DX(:,p3)+w3)
            CALL min_imageV(3,vect2)
            vect3=DX(:,p3)+w3 - (DX(:,p4)+w4)
            CALL min_imageV(3,vect3)
            vect1=MATMUL(A,vect1)
            vect2=MATMUL(A,vect2)
            vect3=MATMUL(A,vect3)
            tst=DOT_PRODUCT(vect1,vect2)/VECTORSIZE(3,vect1)/VECTORSIZE(3,vect2)
            tst=(1-tst**2)**0.5
            IF (tst>0.15) info=1
            tst=DOT_PRODUCT(vect2,vect3)/VECTORSIZE(3,vect2)/VECTORSIZE(3,vect3)
            tst=(1-tst**2)**0.5
            IF (tst<0.15) info=0
            cross1=CROSSPROD(3,vect1,vect2)          
            cross1=cross1/VECTORSIZE(3,cross1)
            cross2=CROSSPROD(3,vect2,vect3)           
            cross2=cross2/VECTORSIZE(3,cross2)
            torsion=DOT_PRODUCT(cross1,cross2)
            IF (torsion>1) torsion=1_q
            IF (torsion<-1) torsion=-1_q
            torsion=ACOS(torsion)
            IF (DOT_PRODUCT(cross1,vect3)>0) torsion=-torsion         
          END FUNCTION

          FUNCTION CAL_MIDDLE(NIONS,A,DX,p1,p2,p3,w1,w2,w3)  RESULT(dist)
            INTEGER :: NIONS
            REAL(q) :: A(3,3)
            REAL(q) :: DX(3,NIONS)
            REAL(q) :: vect(3)
            INTEGER :: p1,p2,p3
            INTEGER :: w1(3),w2(3),w3(3)
            REAL(q) :: dist

            vect=((DX(:,p2)+w2) + (DX(:,p3)+w3))/2-(DX(:,p1)+w1)
            CALL min_imageV(3,vect)
            vect=MATMUL(A,vect)
            dist=VECTORSIZE(3,vect)
          END FUNCTION CAL_MIDDLE

          FUNCTION CAL_MIDDLE2(NIONS,A,DX,p1,p2,p3,p4,w1,w2,w3,w4)  RESULT(dist)
            INTEGER :: NIONS
            REAL(q) :: A(3,3)
            REAL(q) :: DX(3,NIONS)
            REAL(q) :: vect(3)
            INTEGER :: p1,p2,p3,p4
            INTEGER :: w1(3),w2(3),w3(3),w4(3)
            REAL(q) :: dist

            vect=((DX(:,p3)+w3) + (DX(:,p4)+w4))/2-((DX(:,p1)+w1) + (DX(:,p2)+w2))/2
            CALL min_imageV(3,vect)
            vect=MATMUL(A,vect)
            dist=VECTORSIZE(3,vect)
          END FUNCTION CAL_MIDDLE2

          SUBROUTINE ROT_MOM(T_INFO,A,POS,POS_old,RM)
            TYPE(type_info) :: T_INFO
            REAL(q) :: POS(3,T_INFO%NIONS),POS_old(3,T_INFO%NIONS)
            REAL(q) :: cPOS(3,T_INFO%NIONS),cPOS_old(3,T_INFO%NIONS)
            REAL(q) :: A(3,3)
            REAL(q) :: RM(3,1)
            INTEGER :: i

            RM=0._q
            cPOS=POS;cPOS_old=POS_old
            CALL minimize_difference(cPOS,cPOS_old,T_INFO%NIONS)
            cPOS=MATMUL(A,cPOS);cPOS_old=MATMUL(A,cPOS_old)
            DO i=1,T_INFO%NIONS
              RM(1,1)=RM(1,1)+T_INFO%POMASS(T_INFO%ITYP(i))*(cPOS(2,i)*cPOS_old(3,i)-cPOS(3,i)*cPOS_old(2,i))
              RM(2,1)=RM(2,1)+T_INFO%POMASS(T_INFO%ITYP(i))*(cPOS(3,i)*cPOS_old(1,i)-cPOS(1,i)*cPOS_old(3,i))
              RM(3,1)=RM(3,1)+T_INFO%POMASS(T_INFO%ITYP(i))*(cPOS(1,i)*cPOS_old(2,i)-cPOS(2,i)*cPOS_old(1,i))
            ENDDO

            !DO i=1,T_INFO%NIONS
            !  RM(1,1)=RM(1,1)+T_INFO%POMASS(T_INFO%ITYP(i))*((POS(2,i)-CM(2,1))*CM(3,1)-(POS(3,i)-CM(3,1))*CM(2,1))
            !  RM(2,1)=RM(2,1)+T_INFO%POMASS(T_INFO%ITYP(i))*((POS(3,i)-CM(3,1))*CM(1,1)-(POS(1,i)-CM(1,1))*CM(3,1))
            !  RM(3,1)=RM(3,1)+T_INFO%POMASS(T_INFO%ITYP(i))*((POS(1,i)-CM(1,1))*CM(2,1)-(POS(2,i)-CM(2,1))*CM(1,1))
            !  write(*,*) rm,'rm2'
            !ENDDO
            
          END SUBROUTINE ROT_MOM

          SUBROUTINE BMATRIX(T_INFO,DX,CX,lattmat,ICOORDINATES,BMAT,LCART)
          !c Jacobi matrix for transdformation from fractional (or cartesain) to internal coordinates
            TYPE(type_info) :: T_INFO
            TYPE(coordstructure) :: ICOORDINATES
            INTEGER :: i,j,k
            REAL(q) :: lattmat(3,3)
            REAL(q) :: DX(3,T_INFO%NIONS),CX(3,T_INFO%NIONS)
            REAL(q),DIMENSION(:,:) :: BMAT
            INTEGER :: dimn(2),dummy
            REAL(q) :: complexcoord,dummyA,dummyB,dummyC,dummyD
            REAL(q) :: dummyp,dummyq,dummyr
            LOGICAL :: LCART
            REAL(q),ALLOCATABLE :: brow(:),brow2(:),brow3(:)  
            INTEGER :: ii,jj,kk

            dimn=SHAPE(BMAT)
            BMAT=0._q
            ALLOCATE(brow(dimn(2)),brow2(dimn(2)),brow3(dimn(2)))
            brow=0._q;brow2=0._q;brow3=0._q
            

            DO i=1,ICOORDINATES%NUMINTERNALS 
              SELECT CASE (ICOORDINATES%COORDSTRUCT(i)%TAG)
                CASE ('R ')
                  BMAT(i,:)=PREPARE_L(T_INFO%NIONS,DX,lattmat,ICOORDINATES%COORDSTRUCT(i),dimn(2),LCART)
                CASE ('Q ')
                  BMAT(i,:)=PREPARE_Q(T_INFO%NIONS,DX,lattmat,ICOORDINATES%COORDSTRUCT(i),dimn(2),LCART)
                CASE ('A ')
                  BMAT(i,:)=PREPARE_A(T_INFO%NIONS,DX,lattmat,ICOORDINATES%COORDSTRUCT(i),dimn(2),LCART)
                CASE ('T ')
                  BMAT(i,:)=PREPARE_T(T_INFO%NIONS,DX,lattmat,ICOORDINATES%COORDSTRUCT(i),dimn(2),LCART)
                CASE ('X ','Y ','Z ')
                  BMAT(i,:)=PREPARE_XYZ(DX,lattmat,ICOORDINATES%COORDSTRUCT(i),dimn(2),LCART)
                CASE ('cX','cY','cZ')
                  BMAT(i,:)=PREPARE_cXYZ(DX,lattmat,ICOORDINATES%COORDSTRUCT(i),dimn(2),LCART)  
                CASE ('M ')
                  BMAT(i,:)=PREPARE_M(T_INFO%NIONS,DX,lattmat,ICOORDINATES%COORDSTRUCT(i),dimn(2),LCART)
                CASE ('B ')
                  BMAT(i,:)=PREPARE_B(T_INFO%NIONS,DX,lattmat,ICOORDINATES%COORDSTRUCT(i),dimn(2),LCART)
                CASE ('P ')
                  BMAT(i,:)=PREPARE_P(T_INFO%NIONS,DX,lattmat,ICOORDINATES%COORDSTRUCT(i),dimn(2),LCART)
                CASE ('W ')
                  brow=PREPARE_L(T_INFO%NIONS,DX,lattmat,ICOORDINATES%COORDSTRUCT(i),dimn(2),LCART)
                  
                  dummyB=CAL_DISTANCE(T_INFO%NIONS,TRANSPOSE(lattmat),&
                  &  DX,ICOORDINATES%COORDSTRUCT(i)%what(1),ICOORDINATES%COORDSTRUCT(i)%what(2), &
                  &  ICOORDINATES%COORDSTRUCT(i)%WHERE(1,:),ICOORDINATES%COORDSTRUCT(i)%WHERE(2,:))
                  
                  dummyA=dummyB/ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(1)
                  dummyC=1.0-(dummyA)**ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(2)
                  dummyD=1.0-(dummyA)**ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(3)
                  BMAT(i,:)=-ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(2)*brow*&
                  &  (dummyA**ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(2)/dummyB)/dummyD
                  BMAT(i,:)=BMAT(i,:)+ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(3)*dummyC*brow*&
                  &  (dummyA**ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(3)/dummyB)/dummyD**2
                CASE('WP')
                  DO ii=-ICOORDINATES%COORDSTRUCT(i)%TRANS(1),ICOORDINATES%COORDSTRUCT(i)%TRANS(1)
                    DO jj=-ICOORDINATES%COORDSTRUCT(i)%TRANS(2),ICOORDINATES%COORDSTRUCT(i)%TRANS(2)
                      DO kk=-ICOORDINATES%COORDSTRUCT(i)%TRANS(3),ICOORDINATES%COORDSTRUCT(i)%TRANS(3)
                        brow=PREPARE_L_trans(T_INFO%NIONS,DX,lattmat,ICOORDINATES%COORDSTRUCT(i), &
                        & (/ii,jj,kk/),dimn(2),LCART)

                        dummyB=CAL_DISTANCE_trans(T_INFO%NIONS,TRANSPOSE(lattmat),&
                        &  DX,ICOORDINATES%COORDSTRUCT(i)%what(1),ICOORDINATES%COORDSTRUCT(i)%what(2), &
                        &  ICOORDINATES%COORDSTRUCT(i)%WHERE(1,:),ICOORDINATES%COORDSTRUCT(i)%WHERE(2,:),&
                        &  (/ii,jj,kk/))

                        dummyA=dummyB/ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(1)
                        dummyC=1.0-(dummyA)**ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(2)
                        dummyD=1.0-(dummyA)**ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(3)
                        BMAT(i,:)=BMAT(i,:)-ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(2)*brow*&
                        &  (dummyA**ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(2)/dummyB)/dummyD
                        BMAT(i,:)=BMAT(i,:)+ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(3)*dummyC*brow*&
                        &  (dummyA**ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(3)/dummyB)/dummyD**2
                      ENDDO
                    ENDDO
                  ENDDO
                !CASE('E ')
                !  BMAT(i,:)=PREPARE_ROT(T_INFO,CX,lattmat,dimn(2),1,LCART)
                !CASE('F ')
                !  BMAT(i,:)=PREPARE_ROT(T_INFO,CX,lattmat,dimn(2),2,LCART)
                !CASE('G ')
                !  BMAT(i,:)=PREPARE_ROT(T_INFO,CX,lattmat,dimn(2),3,LCART)
                CASE ('LR')
                  BMAT(i,:)=PREPARE_LR(T_INFO%NIONS,lattmat,ICOORDINATES%COORDSTRUCT(i),dimn(2))
                CASE ('LA')
                  BMAT(i,:)=PREPARE_LA(T_INFO%NIONS,lattmat,ICOORDINATES%COORDSTRUCT(i),dimn(2))
                CASE ('LV')
                  BMAT(i,:)=PREPARE_LV(T_INFO%NIONS,lattmat,dimn(2))
                CASE ('C ')
                  complexcoord=0
                  DO j=1,SIZE(ICOORDINATES%COORDSTRUCT(i)%COEFS)
                    IF (ABS(ICOORDINATES%COORDSTRUCT(i)%COEFS(j))>0.) THEN
                      BMAT(i,:)=BMAT(i,:)+BMAT(j,:)*&
                      ICOORDINATES%COORDSTRUCT(i)%COEFS(j)**2*ICOORDINATES%COORDSTRUCT(j)%VALUE
                      complexcoord=complexcoord+((ICOORDINATES%COORDSTRUCT(i)%COEFS(j)**2)*&
                      ICOORDINATES%COORDSTRUCT(j)%VALUE**2)
                    ENDIF
                  ENDDO
                  BMAT(i,:)=BMAT(i,:)/complexcoord**(0.5)
                CASE ('S ')
                  DO j=1,SIZE(ICOORDINATES%COORDSTRUCT(i)%COEFS)
                    IF (ABS(ICOORDINATES%COORDSTRUCT(i)%COEFS(j))>0.) THEN
                      BMAT(i,:)=BMAT(i,:)+BMAT(j,:)*ICOORDINATES%COORDSTRUCT(i)%COEFS(j)
                    ENDIF
                  ENDDO
                CASE ('ST')
                  DO j=1,SIZE(ICOORDINATES%COORDSTRUCT(i)%COEFS)
                    IF ((ABS(ICOORDINATES%COORDSTRUCT(i)%COEFS(j))>0.) .AND. (ICOORDINATES%COORDSTRUCT(j)%TAG=='T ')) THEN
                      BMAT(i,:)=BMAT(i,:)+BMAT(j,:)*ICOORDINATES%COORDSTRUCT(i)%COEFS(j)
                    ENDIF
                  ENDDO
                CASE('IS')
                  dummyp=0.
                  dummyq=0.
                  brow2=0.
                  brow3=0.
                  DO k=1,ICOORDINATES%IRCDIM(2)
                    dummyr=0.
                    brow=0.
                    DO j=1,ICOORDINATES%NUMPRIMITIVES
                      IF (ICOORDINATES%COORDSTRUCT(i)%COEFS(j)>0.) THEN
                        !dummyr=dummyr+ICOORDINATES%COORDSTRUCT(i)%COEFS(j)**2 * &
                        !&   (ICOORDINATES%COORDSTRUCT(j)%VALUE-ICOORDINATES%IRCPOINTS(j,k))**2
                        !brow=brow+ICOORDINATES%COORDSTRUCT(i)%COEFS(j)**2 * &
                        !&   (ICOORDINATES%COORDSTRUCT(j)%VALUE-ICOORDINATES%IRCPOINTS(j,k))*BMAT(j,:)
                        dummyr=dummyr+ICOORDINATES%COORDSTRUCT(i)%COEFS(j) * &
                        &   (ICOORDINATES%COORDSTRUCT(j)%VALUE-ICOORDINATES%IRCPOINTS(j,k))**2
                        brow=brow+ICOORDINATES%COORDSTRUCT(i)%COEFS(j) * &
                        &   (ICOORDINATES%COORDSTRUCT(j)%VALUE-ICOORDINATES%IRCPOINTS(j,k))*BMAT(j,:)
                      ENDIF
                    ENDDO

                    dummyr=exp(-dummyr)
                    brow=2*brow*dummyr
                    
                    dummyp=dummyp+(k-1)*dummyr
                    dummyq=dummyq+dummyr
                    brow2=brow2+(k-1)*brow
                    brow3=brow3+brow
                  ENDDO  
                  IF (dummyp>0. .AND. dummyq>0.) THEN
                    dummyD=dummyp/(dummyq*(ICOORDINATES%IRCDIM(2)-1))
                    BMAT(i,:)=dummyD*brow3/dummyq - &
                    &  dummyD*brow2/dummyp
                  ENDIF  
                CASE('IZ')
                  dummyq=0.
                  brow3=0.
                  DO k=1,ICOORDINATES%IRCDIM(2)
                    dummyr=0.
                    brow=0.
                    DO j=1,ICOORDINATES%NUMPRIMITIVES
                      IF (ICOORDINATES%COORDSTRUCT(i)%COEFS(j)>0.) THEN
                        dummyr=dummyr+ICOORDINATES%COORDSTRUCT(i)%COEFS(j) * &
                        &   (ICOORDINATES%COORDSTRUCT(j)%VALUE-ICOORDINATES%IRCPOINTS(j,k))**2
                        brow=brow+(ICOORDINATES%COORDSTRUCT(j)%VALUE-ICOORDINATES%IRCPOINTS(j,k))*BMAT(j,:)
                      ENDIF
                    ENDDO
                
                    dummyr=exp(-dummyr)
                    brow=2*brow*dummyr
                    
                    dummyq=dummyq+dummyr
                    brow3=brow3+brow
                  ENDDO  
                  IF (dummyq>0.) THEN
                    BMAT(i,:)=brow3/dummyq
                  ENDIF
                  
                CASE ('D ')
                  DO j=1,SIZE(ICOORDINATES%COORDSTRUCT(i)%COEFS)
                    IF (ICOORDINATES%COORDSTRUCT(i)%COEFS(j)>1e-4) THEN 
                      dummyA=ICOORDINATES%COORDSTRUCT(j)%VALUE/ABS(ICOORDINATES%COORDSTRUCT(i)%COEFS(j))
                      IF (ABS(dummyA-1._q)<1e-4) dummyA=1._q+1e-4
                      dummyC=1.0-(dummyA)**ICOORDINATES%CNEXP(1)
                      dummyD=1.0-(dummyA)**ICOORDINATES%CNEXP(2)
                      BMAT(i,:)=BMAT(i,:)-ICOORDINATES%CNEXP(1)*BMAT(j,:)*(dummyA**ICOORDINATES%CNEXP(1)/ICOORDINATES%COORDSTRUCT(j)%VALUE)/dummyD
                      BMAT(i,:)=BMAT(i,:)+ICOORDINATES%CNEXP(2)*dummyC*BMAT(j,:)*(dummyA**ICOORDINATES%CNEXP(2)/ICOORDINATES%COORDSTRUCT(j)%VALUE)&
                                /dummyD**2
                    ELSE IF (ICOORDINATES%COORDSTRUCT(i)%COEFS(j)<-1e-4) THEN
                      dummyA=ICOORDINATES%COORDSTRUCT(j)%VALUE/ABS(ICOORDINATES%COORDSTRUCT(i)%COEFS(j))
                      IF (ABS(dummyA-1._q)<1e-4) dummyA=1._q+1e-4
                      dummyC=1.0-(dummyA)**ICOORDINATES%CNEXP(1)
                      dummyD=1.0-(dummyA)**ICOORDINATES%CNEXP(2)
                      BMAT(i,:)=BMAT(i,:)+ICOORDINATES%CNEXP(1)*BMAT(j,:)*(dummyA**ICOORDINATES%CNEXP(1)/ICOORDINATES%COORDSTRUCT(j)%VALUE)/dummyD
                      BMAT(i,:)=BMAT(i,:)-ICOORDINATES%CNEXP(2)*dummyC*BMAT(j,:)*(dummyA**ICOORDINATES%CNEXP(2)/ICOORDINATES%COORDSTRUCT(j)%VALUE)&
                                /dummyD**2

                    ENDIF

                  ENDDO
              END SELECT
            ENDDO
  
            DEALLOCATE(brow3,brow2,brow)
            
            CONTAINS
            FUNCTION PREPARE_XYZ(DX,lattmat, COORDSTRUCT,din,LCART) RESULT(ROW)
            !c B-matrix components for fractional coordinates
              TYPE(coordinate) :: COORDSTRUCT
              INTEGER :: din
              REAL(q),DIMENSION(:,:) :: DX
              REAL(q) :: ROW(din)
              REAL(q) :: lattmat(3,3),lattinv(3,3)
              INTEGER :: dummy,ninfo
              LOGICAL :: LCART

              ROW=0_q
              SELECT CASE (COORDSTRUCT%TAG)
                CASE ('X ')
                  dummy=2
                CASE ('Y ')
                  dummy=1
                CASE ('Z ')
                  dummy=0
              END SELECT
              ROW(3*COORDSTRUCT%WHAT(1)-dummy)=1._q

              !c transform to Cartesian coordinates if desired
              IF (LCART) THEN
                 lattinv=lattmat
                 CALL SVDINVERSE(lattinv,3,ninfo)
                 ROW(3*COORDSTRUCT%WHAT(1)-2:3*COORDSTRUCT%WHAT(1))=MATMUL(ROW(3*COORDSTRUCT%WHAT(1)-2:3*COORDSTRUCT%WHAT(1)),TRANSPOSE(lattinv))
              ENDIF 
            END FUNCTION

            FUNCTION PREPARE_cXYZ(DX,lattmat, COORDSTRUCT,din,LCART) RESULT(ROW)
            !c B-matrix components for Cartesian coordinates
              TYPE(coordinate) :: COORDSTRUCT
              INTEGER :: din
              REAL(q),DIMENSION(:,:) :: DX
              REAL(q) :: ROW(din)
              REAL(q) :: lattmat(3,3),lattinv(3,3)
              INTEGER :: dummy,ninfo
              LOGICAL :: LCART

              ROW=0_q
              SELECT CASE (COORDSTRUCT%TAG)
                CASE ('cX')
                  dummy=2
                CASE ('cY')
                  dummy=1
                CASE ('cZ')
                  dummy=0
              END SELECT
              ROW(3*COORDSTRUCT%WHAT(1)-dummy)=1._q

              !c transform to fractional coordinates if desired
              IF (.NOT. LCART) THEN
                 ROW(3*COORDSTRUCT%WHAT(1)-2:3*COORDSTRUCT%WHAT(1))=MATMUL(ROW(3*COORDSTRUCT%WHAT(1)-2:3*COORDSTRUCT%WHAT(1)),TRANSPOSE(lattmat))                
              ENDIF 
            END FUNCTION
            
            FUNCTION PREPARE_L(NIONS,DX,lattmat,COORDSTRUCT,din,LCART) RESULT(ROW)
              TYPE(coordinate) :: COORDSTRUCT
              INTEGER :: NIONS,din
              REAL(q),DIMENSION(:,:) :: DX
              REAL(q) :: lattmat(3,3)
              REAL(q) :: pa(3),pb(3),vect(3)
              REAL(q) :: ROW(din)
              REAL(q) :: LATDER(3,3)
              LOGICAL :: LCART

              ROW=0_q
              pa=DX(:,COORDSTRUCT%WHAT(1))+COORDSTRUCT%WHERE(1,:)
              pb=DX(:,COORDSTRUCT%WHAT(2))+COORDSTRUCT%WHERE(2,:)
              vect=(pa-pb)
              CALL min_imageV(3,vect)
              vect=MATMUL(vect,lattmat)
              !vect=vect/COORDSTRUCT%VALUE
              vect=vect/(vect(1)*vect(1)+vect(2)*vect(2)+vect(3)*vect(3))**0.5

              pa=vect
              pb=-vect

              !c take care of frozen fractional coordinates
              !IF (.NOT. T_INFO%LSFOR(1,COORDSTRUCT%WHAT(1))) pa(1)=0._q
              !IF (.NOT. T_INFO%LSFOR(2,COORDSTRUCT%WHAT(1))) pa(2)=0._q
              !IF (.NOT. T_INFO%LSFOR(3,COORDSTRUCT%WHAT(1))) pa(3)=0._q
              !IF (.NOT. T_INFO%LSFOR(1,COORDSTRUCT%WHAT(2))) pb(1)=0._q
              !IF (.NOT. T_INFO%LSFOR(2,COORDSTRUCT%WHAT(2))) pb(2)=0._q
              !IF (.NOT. T_INFO%LSFOR(3,COORDSTRUCT%WHAT(2))) pb(3)=0._q

              !c convert to fractional coordinates if desired
              IF (.NOT. LCART) THEN
                pa=MATMUL(pa,TRANSPOSE(lattmat))
                pb=MATMUL(pb,TRANSPOSE(lattmat))
              ENDIF

              ROW(3*COORDSTRUCT%WHAT(1)-3+1:3*COORDSTRUCT%WHAT(1)-3+3)=pa
              ROW(3*COORDSTRUCT%WHAT(2)-3+1:3*COORDSTRUCT%WHAT(2)-3+3)=pb
              IF (din==3*NIONS+9) THEN
                LATDER=0._q
                LATDER=OUTERPROD2(3,COORDSTRUCT%WHERE(1,:),pa)+&
                         OUTERPROD2(3,COORDSTRUCT%WHERE(2,:),pb)
                ROW(3*NIONS+1:3*NIONS+3)=LATDER(1,:)
                ROW(3*NIONS+4:3*NIONS+6)=LATDER(2,:)
                ROW(3*NIONS+7:3*NIONS+9)=LATDER(3,:)
              ENDIF
            END FUNCTION PREPARE_L

            FUNCTION PREPARE_L_trans(NIONS,DX,lattmat,COORDSTRUCT,trans,din,LCART) RESULT(ROW)
              TYPE(coordinate) :: COORDSTRUCT
              INTEGER :: NIONS,din
              REAL(q),DIMENSION(:,:) :: DX
              REAL(q) :: lattmat(3,3)
              REAL(q) :: pa(3),pb(3),vect(3)
              REAL(q) :: ROW(din)
              REAL(q) :: LATDER(3,3)
              INTEGER :: trans(3)
              LOGICAL :: LCART

              ROW=0_q
              pa=DX(:,COORDSTRUCT%WHAT(1))+COORDSTRUCT%WHERE(1,:)
              pb=DX(:,COORDSTRUCT%WHAT(2))+COORDSTRUCT%WHERE(2,:)
              vect=(pa-pb)
              CALL min_imageV(3,vect)
              vect=MATMUL(vect-trans,lattmat)
              !vect=vect/COORDSTRUCT%VALUE
              vect=vect/(vect(1)*vect(1)+vect(2)*vect(2)+vect(3)*vect(3))**0.5

              pa=vect
              pb=-vect

              !c convert to fractional coordinates if desired
              IF (.NOT. LCART) THEN
                pa=MATMUL(pa,TRANSPOSE(lattmat))
                pb=MATMUL(pb,TRANSPOSE(lattmat))
              ENDIF

              ROW(3*COORDSTRUCT%WHAT(1)-3+1:3*COORDSTRUCT%WHAT(1)-3+3)=pa
              ROW(3*COORDSTRUCT%WHAT(2)-3+1:3*COORDSTRUCT%WHAT(2)-3+3)=pb
              IF (din==3*NIONS+9) THEN
                LATDER=0._q
                LATDER=OUTERPROD2(3,COORDSTRUCT%WHERE(1,:),pa)+&
                         OUTERPROD2(3,COORDSTRUCT%WHERE(2,:),pb)
                ROW(3*NIONS+1:3*NIONS+3)=LATDER(1,:)
                ROW(3*NIONS+4:3*NIONS+6)=LATDER(2,:)
                ROW(3*NIONS+7:3*NIONS+9)=LATDER(3,:)
              ENDIF
            END FUNCTION PREPARE_L_trans

            FUNCTION PREPARE_LR(NIONS,lattmat,COORDSTRUCT,din) RESULT(ROW)
              TYPE(coordinate) :: COORDSTRUCT
              INTEGER :: NIONS,din
              REAL(q) :: lattmat(3,3)
              REAL(q) :: ROW(din)
              REAL(q) :: dist

              ROW=0_q
              IF (din==3*NIONS+9) THEN
                dist=CAL_LR(transpose(lattmat),COORDSTRUCT%WHAT(1))
                ROW(3*NIONS+3*(COORDSTRUCT%WHAT(1)-1)+1)=lattmat(COORDSTRUCT%WHAT(1),1)/dist 
                ROW(3*NIONS+3*(COORDSTRUCT%WHAT(1)-1)+2)=lattmat(COORDSTRUCT%WHAT(1),2)/dist
                ROW(3*NIONS+3*(COORDSTRUCT%WHAT(1)-1)+3)=lattmat(COORDSTRUCT%WHAT(1),3)/dist
              ENDIF
            END FUNCTION PREPARE_LR

            FUNCTION PREPARE_LA(NIONS,lattmat,COORDSTRUCT,din) RESULT(ROW)
              TYPE(coordinate) :: COORDSTRUCT
              INTEGER :: NIONS,din
              REAL(q) :: lattmat(3,3)
              REAL(q) :: ROW(din)
              REAL(q) :: diffav(3),diffbv(3),dala(3),dalc(3)
              REAL(q) :: d1,d2,alpha

              ROW=0_q
              IF (din==3*NIONS+9) THEN
                diffav=lattmat(COORDSTRUCT%WHAT(1),:)
                diffbv=lattmat(COORDSTRUCT%WHAT(2),:)
                d1=CAL_LR(transpose(lattmat),COORDSTRUCT%WHAT(1))
                d2=CAL_LR(transpose(lattmat),COORDSTRUCT%WHAT(2))
                alpha=CAL_LA(transpose(lattmat),COORDSTRUCT%WHAT(1),COORDSTRUCT%WHAT(2))
                dala=-(diffbv/(d1*d2)-diffav*COS(alpha)/(d1**2))/SIN(alpha)
                dalc=-(diffav/(d1*d2)-diffbv*COS(alpha)/(d2**2))/SIN(alpha)
                ROW(3*NIONS+3*(COORDSTRUCT%WHAT(1)-1)+1)=dala(1)
                ROW(3*NIONS+3*(COORDSTRUCT%WHAT(1)-1)+2)=dala(2)
                ROW(3*NIONS+3*(COORDSTRUCT%WHAT(1)-1)+3)=dala(3) 
                ROW(3*NIONS+3*(COORDSTRUCT%WHAT(2)-1)+1)=dalc(1)
                ROW(3*NIONS+3*(COORDSTRUCT%WHAT(2)-1)+2)=dalc(2)
                ROW(3*NIONS+3*(COORDSTRUCT%WHAT(2)-1)+3)=dalc(3) 
              ENDIF
            END FUNCTION PREPARE_LA
            

            FUNCTION PREPARE_LV(NIONS,lattmat,din) RESULT(ROW)
              INTEGER :: NIONS,din
              REAL(q) :: lattmat(3,3)
              REAL(q) :: ROW(din)
              REAL(q) :: dl1(3),dl2(3),dl3(3)

              ROW=0_q
              IF (din==3*NIONS+9) THEN
                dl1=CROSSPROD(3,lattmat(2,:),lattmat(3,:))
                dl2=CROSSPROD(3,lattmat(3,:),lattmat(1,:))
                dl3=CROSSPROD(3,lattmat(1,:),lattmat(2,:))          
                ROW(3*NIONS+1:3*NIONS+3)=dl1 
                ROW(3*NIONS+4:3*NIONS+6)=dl2
                ROW(3*NIONS+7:3*NIONS+9)=dl3
              ENDIF
            END FUNCTION PREPARE_LV         


            FUNCTION PREPARE_Q(NIONS,DX,lattmat,COORDSTRUCT,din,LCART) RESULT(ROW)
              TYPE(coordinate) :: COORDSTRUCT
              INTEGER :: NIONS,din
              REAL(q),DIMENSION(:,:) :: DX
              REAL(q) :: lattmat(3,3)
              REAL(q) :: pa(3),pb(3),vect(3)
              REAL(q) :: ROW(din)
              REAL(q) :: LATDER(3,3)
              LOGICAL :: LCART

              ROW=0_q
              pa=DX(:,COORDSTRUCT%WHAT(1))+COORDSTRUCT%WHERE(1,:)
              pb=DX(:,COORDSTRUCT%WHAT(2))+COORDSTRUCT%WHERE(2,:)
              vect=(pa-pb)
              CALL min_imageV(3,vect)
              vect=MATMUL(vect,lattmat)
              vect=2*vect

              pa=vect
              pb=-vect

              !c take care of frozen fractional coordinates
              !IF (.NOT. T_INFO%LSFOR(1,COORDSTRUCT%WHAT(1))) pa(1)=0._q
              !IF (.NOT. T_INFO%LSFOR(2,COORDSTRUCT%WHAT(1))) pa(2)=0._q
              !IF (.NOT. T_INFO%LSFOR(3,COORDSTRUCT%WHAT(1))) pa(3)=0._q
              !IF (.NOT. T_INFO%LSFOR(1,COORDSTRUCT%WHAT(2))) pb(1)=0._q
              !IF (.NOT. T_INFO%LSFOR(2,COORDSTRUCT%WHAT(2))) pb(2)=0._q
              !IF (.NOT. T_INFO%LSFOR(3,COORDSTRUCT%WHAT(2))) pb(3)=0._q

              !c convert to fractional coordinates if desired
              IF (.NOT. LCART) THEN
                pa=MATMUL(pa,TRANSPOSE(lattmat))
                pb=MATMUL(pb,TRANSPOSE(lattmat))
              END IF

              ROW(3*COORDSTRUCT%WHAT(1)-3+1:3*COORDSTRUCT%WHAT(1)-3+3)=pa
              ROW(3*COORDSTRUCT%WHAT(2)-3+1:3*COORDSTRUCT%WHAT(2)-3+3)=pb
              IF (din==3*NIONS+9) THEN
                LATDER=0._q
                LATDER=OUTERPROD2(3,COORDSTRUCT%WHERE(1,:),pa)+&
                         OUTERPROD2(3,COORDSTRUCT%WHERE(2,:),pb)
                ROW(3*NIONS+1:3*NIONS+3)=LATDER(1,:)
                ROW(3*NIONS+4:3*NIONS+6)=LATDER(2,:)
                ROW(3*NIONS+7:3*NIONS+9)=LATDER(3,:)
              ENDIF
            END FUNCTION PREPARE_Q


            FUNCTION PREPARE_A(NIONS,DX,lattmat,COORDSTRUCT,din,LCART) RESULT(ROW)
              TYPE(coordinate) :: COORDSTRUCT
              INTEGER :: NIONS,din
              REAL(q),DIMENSION(:,:) :: DX
              REAL(q) :: lattmat(3,3)
              REAL(q) :: pa(3),pb(3),pc(3),vect1(3),vect2(3)
              REAL(q) :: ROW(din)
              REAL(q) :: LATDER(3,3)
              LOGICAL :: LCART

              ROW=0_q
              pa=DX(:,COORDSTRUCT%WHAT(1))+COORDSTRUCT%WHERE(1,:)
              pb=DX(:,COORDSTRUCT%WHAT(2))+COORDSTRUCT%WHERE(2,:)
              pc=DX(:,COORDSTRUCT%WHAT(3))+COORDSTRUCT%WHERE(3,:)
              vect1=(pa-pb)
              CALL min_imageV(3,vect1)
              vect1=MATMUL(vect1,lattmat)
              vect2=(pc-pb)
              CALL min_imageV(3,vect2)
              vect2=MATMUL(vect2,lattmat)

              pa=-(vect2/(VECTORSIZE(3,vect1)*VECTORSIZE(3,vect2))- &
                   vect1*COS(COORDSTRUCT%VALUE)/(VECTORSIZE(3,vect1))**2)/SIN(COORDSTRUCT%VALUE)
              pc=-(vect1/(VECTORSIZE(3,vect1)*VECTORSIZE(3,vect2))- &
                 vect2*COS(COORDSTRUCT%VALUE)/(VECTORSIZE(3,vect2))**2)/SIN(COORDSTRUCT%VALUE)
              pb=-pa-pc

              !c take care of frozen fractional coordinates
              !IF (.NOT. T_INFO%LSFOR(1,COORDSTRUCT%WHAT(1))) pa(1)=0._q
              !IF (.NOT. T_INFO%LSFOR(2,COORDSTRUCT%WHAT(1))) pa(2)=0._q
              !IF (.NOT. T_INFO%LSFOR(3,COORDSTRUCT%WHAT(1))) pa(3)=0._q
              !IF (.NOT. T_INFO%LSFOR(1,COORDSTRUCT%WHAT(2))) pb(1)=0._q
              !IF (.NOT. T_INFO%LSFOR(2,COORDSTRUCT%WHAT(2))) pb(2)=0._q
              !IF (.NOT. T_INFO%LSFOR(3,COORDSTRUCT%WHAT(2))) pb(3)=0._q
              !IF (.NOT. T_INFO%LSFOR(1,COORDSTRUCT%WHAT(3))) pc(1)=0._q
              !IF (.NOT. T_INFO%LSFOR(2,COORDSTRUCT%WHAT(3))) pc(2)=0._q
              !IF (.NOT. T_INFO%LSFOR(3,COORDSTRUCT%WHAT(3))) pc(3)=0._q

            

              !c transform to fractional coordinates if desired
              IF (.NOT. LCART) THEN
                pa=MATMUL(pa,TRANSPOSE(lattmat))
                pb=MATMUL(pb,TRANSPOSE(lattmat))
                pc=MATMUL(pc,TRANSPOSE(lattmat))
              END IF 

              IF (COORDSTRUCT%WHAT(1)==COORDSTRUCT%WHAT(3)) pa=pa+pc

          
              ROW(3*COORDSTRUCT%WHAT(1)-3+1:3*COORDSTRUCT%WHAT(1)-3+3)=pa
              ROW(3*COORDSTRUCT%WHAT(2)-3+1:3*COORDSTRUCT%WHAT(2)-3+3)=pb
              IF (COORDSTRUCT%WHAT(1) /= COORDSTRUCT%WHAT(3)) THEN
                ROW(3*COORDSTRUCT%WHAT(3)-3+1:3*COORDSTRUCT%WHAT(3)-3+3)=pc
              ENDIF
              IF (din==3*NIONS+9) THEN
                LATDER=0._q
                LATDER=OUTERPROD2(3,COORDSTRUCT%WHERE(1,:),pa)+ &
                       OUTERPROD2(3,COORDSTRUCT%WHERE(2,:),pb)+ &
                       OUTERPROD2(3,COORDSTRUCT%WHERE(3,:),pc)
                ROW(3*NIONS+1:3*NIONS+3)=LATDER(1,:)
                ROW(3*NIONS+4:3*NIONS+6)=LATDER(2,:)
                ROW(3*NIONS+7:3*NIONS+9)=LATDER(3,:) 
              ENDIF          
            END FUNCTION

            FUNCTION PREPARE_M(NIONS,DX,lattmat,COORDSTRUCT,din,LCART) RESULT(ROW)
              TYPE(coordinate) :: COORDSTRUCT
              INTEGER :: NIONS,din
              REAL(q),DIMENSION(:,:) :: DX
              REAL(q) :: lattmat(3,3)
              REAL(q) :: pa(3),pb(3),pc(3),vect(3)
              REAL(q) :: ROW(din)
              REAL(q) :: LATDER(3,3)
              LOGICAL :: LCART

              ROW=0_q
              pa=DX(:,COORDSTRUCT%WHAT(1))+COORDSTRUCT%WHERE(1,:)
              pb=DX(:,COORDSTRUCT%WHAT(2))+COORDSTRUCT%WHERE(2,:)
              pc=DX(:,COORDSTRUCT%WHAT(3))+COORDSTRUCT%WHERE(3,:)
              vect=(pb+pc)-2*pa
              CALL min_imageV(3,vect)
              vect=MATMUL(vect,lattmat)
              vect=vect/COORDSTRUCT%VALUE

              !pa=MATMUL(-vect/2,TRANSPOSE(lattmat))
              !pb=MATMUL(vect/4,TRANSPOSE(lattmat))
              !pc=MATMUL(vect/4,TRANSPOSE(lattmat))

              pa=-vect/2
              pb=vect/4
              pc=vect/4

              !c take care of frozen fractional coordinates
              !IF (.NOT. T_INFO%LSFOR(1,COORDSTRUCT%WHAT(1))) pa(1)=0._q
              !IF (.NOT. T_INFO%LSFOR(2,COORDSTRUCT%WHAT(1))) pa(2)=0._q
              !IF (.NOT. T_INFO%LSFOR(3,COORDSTRUCT%WHAT(1))) pa(3)=0._q
              !IF (.NOT. T_INFO%LSFOR(1,COORDSTRUCT%WHAT(2))) pb(1)=0._q
              !IF (.NOT. T_INFO%LSFOR(2,COORDSTRUCT%WHAT(2))) pb(2)=0._q
              !IF (.NOT. T_INFO%LSFOR(3,COORDSTRUCT%WHAT(2))) pb(3)=0._q
              !IF (.NOT. T_INFO%LSFOR(1,COORDSTRUCT%WHAT(3))) pc(1)=0._q
              !IF (.NOT. T_INFO%LSFOR(2,COORDSTRUCT%WHAT(3))) pc(2)=0._q
              !IF (.NOT. T_INFO%LSFOR(3,COORDSTRUCT%WHAT(3))) pc(3)=0._q

              !c transform to fractional coordinates if desired
              IF (.NOT. LCART) THEN
                pa=MATMUL(pa,TRANSPOSE(lattmat))
                pb=MATMUL(pb,TRANSPOSE(lattmat))
                pc=MATMUL(pc,TRANSPOSE(lattmat))
              END IF
             
              ROW(3*COORDSTRUCT%WHAT(1)-3+1:3*COORDSTRUCT%WHAT(1)-3+3)=pa
              ROW(3*COORDSTRUCT%WHAT(2)-3+1:3*COORDSTRUCT%WHAT(2)-3+3)=pb
              ROW(3*COORDSTRUCT%WHAT(3)-3+1:3*COORDSTRUCT%WHAT(3)-3+3)=pc

              IF (din==3*NIONS+9) THEN
                LATDER=0._q
                LATDER=OUTERPROD2(3,COORDSTRUCT%WHERE(1,:),pa)+ &
                       OUTERPROD2(3,COORDSTRUCT%WHERE(2,:),pb)+ &
                       OUTERPROD2(3,COORDSTRUCT%WHERE(3,:),pc)
                ROW(3*NIONS+1:3*NIONS+3)=LATDER(1,:)
                ROW(3*NIONS+4:3*NIONS+6)=LATDER(2,:)
                ROW(3*NIONS+7:3*NIONS+9)=LATDER(3,:) 
              ENDIF          
            END FUNCTION PREPARE_M

            FUNCTION PREPARE_B(NIONS,DX,lattmat,COORDSTRUCT,din,LCART) RESULT(ROW)
              TYPE(coordinate) :: COORDSTRUCT
              INTEGER :: NIONS,din
              REAL(q),DIMENSION(:,:) :: DX
              REAL(q) :: lattmat(3,3)
              REAL(q) :: pa(3),pb(3),pc(3),pd(3),vect(3)
              REAL(q) :: ROW(din)
              REAL(q) :: LATDER(3,3)
              LOGICAL :: LCART

              ROW=0_q
              pa=DX(:,COORDSTRUCT%WHAT(1))+COORDSTRUCT%WHERE(1,:)
              pb=DX(:,COORDSTRUCT%WHAT(2))+COORDSTRUCT%WHERE(2,:)
              pc=DX(:,COORDSTRUCT%WHAT(3))+COORDSTRUCT%WHERE(3,:)
              pd=DX(:,COORDSTRUCT%WHAT(3))+COORDSTRUCT%WHERE(3,:)
              vect=(pa+pb)-(pc+pd)
              CALL min_imageV(3,vect)
              vect=MATMUL(vect,lattmat)
              vect=vect/COORDSTRUCT%VALUE

              pa=vect/4
              pb=vect/4
              pc=-vect/4
              pd=-vect/4

              !c take care of frozen fractional coordinates
              !IF (.NOT. T_INFO%LSFOR(1,COORDSTRUCT%WHAT(1))) pa(1)=0._q
              !IF (.NOT. T_INFO%LSFOR(2,COORDSTRUCT%WHAT(1))) pa(2)=0._q
              !IF (.NOT. T_INFO%LSFOR(3,COORDSTRUCT%WHAT(1))) pa(3)=0._q
              !IF (.NOT. T_INFO%LSFOR(1,COORDSTRUCT%WHAT(2))) pb(1)=0._q
              !IF (.NOT. T_INFO%LSFOR(2,COORDSTRUCT%WHAT(2))) pb(2)=0._q
              !IF (.NOT. T_INFO%LSFOR(3,COORDSTRUCT%WHAT(2))) pb(3)=0._q
              !IF (.NOT. T_INFO%LSFOR(1,COORDSTRUCT%WHAT(3))) pc(1)=0._q
              !IF (.NOT. T_INFO%LSFOR(2,COORDSTRUCT%WHAT(3))) pc(2)=0._q
              !IF (.NOT. T_INFO%LSFOR(3,COORDSTRUCT%WHAT(3))) pc(3)=0._q
              !IF (.NOT. T_INFO%LSFOR(1,COORDSTRUCT%WHAT(3))) pd(1)=0._q
              !IF (.NOT. T_INFO%LSFOR(2,COORDSTRUCT%WHAT(3))) pd(2)=0._q
              !IF (.NOT. T_INFO%LSFOR(3,COORDSTRUCT%WHAT(3))) pd(3)=0._q


              !c transform to fractional coordinates if desired
              IF (.NOT. LCART) THEN
                pa=MATMUL(pa,TRANSPOSE(lattmat))
                pb=MATMUL(pb,TRANSPOSE(lattmat))
                pc=MATMUL(pc,TRANSPOSE(lattmat))
                pd=MATMUL(pd,TRANSPOSE(lattmat))
              END IF


              ROW(3*COORDSTRUCT%WHAT(1)-3+1:3*COORDSTRUCT%WHAT(1)-3+3)=pa
              ROW(3*COORDSTRUCT%WHAT(2)-3+1:3*COORDSTRUCT%WHAT(2)-3+3)=pb
              ROW(3*COORDSTRUCT%WHAT(3)-3+1:3*COORDSTRUCT%WHAT(3)-3+3)=pc
              ROW(3*COORDSTRUCT%WHAT(4)-3+1:3*COORDSTRUCT%WHAT(4)-3+3)=pd

              IF (din==3*NIONS+9) THEN
                LATDER=0._q
                LATDER=OUTERPROD2(3,COORDSTRUCT%WHERE(1,:),pa)+ &
                       OUTERPROD2(3,COORDSTRUCT%WHERE(2,:),pb)+ &
                       OUTERPROD2(3,COORDSTRUCT%WHERE(3,:),pc)+ &
                       OUTERPROD2(3,COORDSTRUCT%WHERE(4,:),pd)
                ROW(3*NIONS+1:3*NIONS+3)=LATDER(1,:)
                ROW(3*NIONS+4:3*NIONS+6)=LATDER(2,:)
                ROW(3*NIONS+7:3*NIONS+9)=LATDER(3,:) 
              ENDIF          
            END FUNCTION PREPARE_B 

            FUNCTION PREPARE_P(NIONS,DX,lattmat,COORDSTRUCT,din,LCART) RESULT(ROW)
              TYPE(coordinate) :: COORDSTRUCT
              INTEGER :: NIONS,din
              REAL(q),DIMENSION(:,:) :: DX
              REAL(q) :: lattmat(3,3)
              REAL(q) :: pa(3),pb(3),pc(3),pd(3),vect(3)
              REAL(q) :: dummyp,dummyq
              REAL(q) :: ROW(din)
              REAL(q) :: LATDER(3,3)
              LOGICAL :: LCART

              ROW=0_q
              pa=DX(:,COORDSTRUCT%WHAT(1))+COORDSTRUCT%WHERE(1,:)
              pb=DX(:,COORDSTRUCT%WHAT(2))+COORDSTRUCT%WHERE(2,:)
              pc=DX(:,COORDSTRUCT%WHAT(3))+COORDSTRUCT%WHERE(3,:)
              pd=DX(:,COORDSTRUCT%WHAT(4))+COORDSTRUCT%WHERE(4,:)
              vect=(pa-pb)
              CALL min_imageV(3,vect)
              vect=MATMUL(vect,lattmat)
              dummyp=(vect(1)*vect(1)+vect(2)*vect(2)+vect(3)*vect(3))**0.5
              vect=vect/dummyp
              pa=vect
              pb=-vect

              vect=(pc-pd)
              CALL min_imageV(3,vect)
              vect=MATMUL(vect,lattmat)
              dummyq=(vect(1)*vect(1)+vect(2)*vect(2)+vect(3)*vect(3))**0.5
              vect=vect/dummyq

              pa=pa/dummyq
              pb=pb/dummyq
              pc=-dummyp/(dummyq*dummyq)*vect
              pd=dummyp/(dummyq*dummyq)*vect

              !c convert to fractional coordinates if desired
              IF (.NOT. LCART) THEN
                pa=MATMUL(pa,TRANSPOSE(lattmat))
                pb=MATMUL(pb,TRANSPOSE(lattmat))
                pc=MATMUL(pc,TRANSPOSE(lattmat))
                pd=MATMUL(pd,TRANSPOSE(lattmat))
              ENDIF

              ROW(3*COORDSTRUCT%WHAT(1)-3+1:3*COORDSTRUCT%WHAT(1)-3+3)=pa
              ROW(3*COORDSTRUCT%WHAT(2)-3+1:3*COORDSTRUCT%WHAT(2)-3+3)=pb
              ROW(3*COORDSTRUCT%WHAT(3)-3+1:3*COORDSTRUCT%WHAT(3)-3+3)=ROW(3*COORDSTRUCT%WHAT(3)-3+1:3*COORDSTRUCT%WHAT(3)-3+3)+pc
              ROW(3*COORDSTRUCT%WHAT(4)-3+1:3*COORDSTRUCT%WHAT(4)-3+3)=ROW(3*COORDSTRUCT%WHAT(4)-3+1:3*COORDSTRUCT%WHAT(4)-3+3)+pd
              IF (din==3*NIONS+9) THEN
                LATDER=0._q
                LATDER=OUTERPROD2(3,COORDSTRUCT%WHERE(1,:),pa)+ &
                       OUTERPROD2(3,COORDSTRUCT%WHERE(2,:),pb)+ &
                       OUTERPROD2(3,COORDSTRUCT%WHERE(3,:),pc)+ &
                       OUTERPROD2(3,COORDSTRUCT%WHERE(4,:),pd)
                ROW(3*NIONS+1:3*NIONS+3)=LATDER(1,:)
                ROW(3*NIONS+4:3*NIONS+6)=LATDER(2,:)
                ROW(3*NIONS+7:3*NIONS+9)=LATDER(3,:)
              ENDIF
            END FUNCTION PREPARE_P

            FUNCTION PREPARE_T(NIONS,DX,lattmat,COORDSTRUCT,din,LCART) RESULT(ROW)
              TYPE(coordinate) :: COORDSTRUCT
              INTEGER :: NIONS,din
              REAL(q),DIMENSION(:,:) :: DX
              REAL(q) :: lattmat(3,3)
              REAL(q) :: pa(3),pb(3),pc(3),pd(3),vect1(3),vect2(3),vect3(3)
              REAL(q) :: dr1,dr2,dr3,cospsi2,cospsi3,sinpsi2,sinpsi3
              REAL(q) :: ROW(din)
              REAL(q) :: LATDER(3,3)
              LOGICAL :: LCART

              ROW=0_q
              pa=DX(:,COORDSTRUCT%WHAT(1))+COORDSTRUCT%WHERE(1,:)
              pb=DX(:,COORDSTRUCT%WHAT(2))+COORDSTRUCT%WHERE(2,:)
              pc=DX(:,COORDSTRUCT%WHAT(3))+COORDSTRUCT%WHERE(3,:)
              pd=DX(:,COORDSTRUCT%WHAT(4))+COORDSTRUCT%WHERE(4,:)
              vect1=(pa-pb)
              CALL min_imageV(3,vect1)
              vect1=MATMUL(vect1,lattmat)
              vect2=(pb-pc)
              CALL min_imageV(3,vect2)
              vect2=MATMUL(vect2,lattmat)
              vect3=(pc-pd)
              CALL min_imageV(3,vect3)
              vect3=MATMUL(vect3,lattmat)

              dr1=VECTORSIZE(3,vect1)
              dr2=VECTORSIZE(3,vect2)
              dr3=VECTORSIZE(3,vect3)
              vect1=vect1/dr1
              vect2=vect2/dr2
              vect3=vect3/dr3
              cospsi2=DOT_PRODUCT(vect1,vect2)
              sinpsi2=(1-cospsi2**2)**0.5
              cospsi3=DOT_PRODUCT(vect2,vect3)
              sinpsi3=(1-cospsi3**2)**0.5
              
              pa=-CROSSPROD(3,vect1,vect2)/(dr1*sinpsi2**2)
              pb= (dr2+dr1*cospsi2)/(dr1*dr2*sinpsi2)* &
                  CROSSPROD(3,vect1,vect2)/sinpsi2+    &
                  cospsi3/(dr2*sinpsi3)*               &
                  CROSSPROD(3,vect2,vect3)/sinpsi3
              pc= -(dr2+dr3*cospsi3)/(dr2*dr3*sinpsi3)* &
                   CROSSPROD(3,vect2,vect3)/sinpsi3-    &
                   cospsi2/(dr2*sinpsi2)*               &
                  CROSSPROD(3,vect1,vect2)/sinpsi2
              pd= CROSSPROD(3,vect2,vect3)/(dr3*sinpsi3**2)             

              !c take care of frozen fractional coordinates
              !IF (.NOT. T_INFO%LSFOR(1,COORDSTRUCT%WHAT(1))) pa(1)=0._q
              !IF (.NOT. T_INFO%LSFOR(2,COORDSTRUCT%WHAT(1))) pa(2)=0._q
              !IF (.NOT. T_INFO%LSFOR(3,COORDSTRUCT%WHAT(1))) pa(3)=0._q
              !IF (.NOT. T_INFO%LSFOR(1,COORDSTRUCT%WHAT(2))) pb(1)=0._q
              !IF (.NOT. T_INFO%LSFOR(2,COORDSTRUCT%WHAT(2))) pb(2)=0._q
              !IF (.NOT. T_INFO%LSFOR(3,COORDSTRUCT%WHAT(2))) pb(3)=0._q
              !IF (.NOT. T_INFO%LSFOR(1,COORDSTRUCT%WHAT(3))) pc(1)=0._q
              !IF (.NOT. T_INFO%LSFOR(2,COORDSTRUCT%WHAT(3))) pc(2)=0._q
              !IF (.NOT. T_INFO%LSFOR(3,COORDSTRUCT%WHAT(3))) pc(3)=0._q
              !IF (.NOT. T_INFO%LSFOR(1,COORDSTRUCT%WHAT(4))) pd(1)=0._q
              !IF (.NOT. T_INFO%LSFOR(2,COORDSTRUCT%WHAT(4))) pd(2)=0._q
              !IF (.NOT. T_INFO%LSFOR(3,COORDSTRUCT%WHAT(4))) pd(3)=0._q


              !c transform to fractional coordinates if desired
              IF (.NOT. LCART) THEN
                pa=MATMUL(pa,TRANSPOSE(lattmat))
                pb=MATMUL(pb,TRANSPOSE(lattmat))
                pc=MATMUL(pc,TRANSPOSE(lattmat))
                pd=MATMUL(pd,TRANSPOSE(lattmat))
              END IF 

              ROW(3*COORDSTRUCT%WHAT(1)-3+1:3*COORDSTRUCT%WHAT(1)-3+3)=pa
              ROW(3*COORDSTRUCT%WHAT(2)-3+1:3*COORDSTRUCT%WHAT(2)-3+3)=pb
              ROW(3*COORDSTRUCT%WHAT(3)-3+1:3*COORDSTRUCT%WHAT(3)-3+3)=pc
              ROW(3*COORDSTRUCT%WHAT(4)-3+1:3*COORDSTRUCT%WHAT(4)-3+3)=pd
              IF (din==3*NIONS+9) THEN 
                LATDER=0._q
                LATDER=OUTERPROD2(3,COORDSTRUCT%WHERE(1,:),pa)+ &
                       OUTERPROD2(3,COORDSTRUCT%WHERE(2,:),pb)+ &
                       OUTERPROD2(3,COORDSTRUCT%WHERE(3,:),pc)+ &
                       OUTERPROD2(3,COORDSTRUCT%WHERE(4,:),pd)
                ROW(3*NIONS+1:3*NIONS+3)=LATDER(1,:)
                ROW(3*NIONS+4:3*NIONS+6)=LATDER(2,:)
                ROW(3*NIONS+7:3*NIONS+9)=LATDER(3,:)  
              ENDIF            
            END FUNCTION
            
            FUNCTION PREPARE_ROT(T_INFO,CX,lattmat,din,tag,LCART) RESULT(ROW)
              TYPE(type_info) :: T_INFO
              INTEGER :: din,i
              REAL(q) :: CX(3,T_INFO%NIONS),CX_(3,T_INFO%NIONS)
              REAL(q) :: lattmat(3,3)
              REAL(q) :: ROW(din)
              INTEGER :: tag
              LOGICAL :: LCART

              ROW=0_q
              CX_=MATMUL(TRANSPOSE(lattmat),CX)
              DO i=1,T_INFO%NIONS
                SELECT CASE(tag)
                  CASE(1)
                    ROW(3*i-3+2)= T_INFO%POMASS(T_INFO%ITYP(i))*CX_(3,i)
                    ROW(3*i-3+3)=-T_INFO%POMASS(T_INFO%ITYP(i))*CX_(2,i)
                  CASE(2)
                    ROW(3*i-3+1)=-T_INFO%POMASS(T_INFO%ITYP(i))*CX_(3,i)
                    ROW(3*i-3+3)= T_INFO%POMASS(T_INFO%ITYP(i))*CX_(1,i)
                  CASE(3)
                    ROW(3*i-3+1)= T_INFO%POMASS(T_INFO%ITYP(i))*CX_(2,i)
                    ROW(3*i-3+2)=-T_INFO%POMASS(T_INFO%ITYP(i))*CX_(1,i)
                END SELECT

                IF (.NOT. LCART) THEN
                  ROW(3*i-3+1:3*i-3+3)=MATMUL(ROW(3*i-3+1:3*i-3+3),TRANSPOSE(lattmat))
                ENDIF
              ENDDO              
            END FUNCTION PREPARE_ROT

          END SUBROUTINE


          SUBROUTINE DEAL_XYZ(T_INFO,DX,CX,A,ICOORDINATES)
          !c compute values of internal coordinates
            TYPE(type_info) :: T_INFO
            TYPE(coordstructure) :: ICOORDINATES
            INTEGER :: i,j,k,info,dummy 
            REAL(q) :: A(3,3)
            REAL(q) :: DX(3,T_INFO%NIONS),CX(3,T_INFO%NIONS)
            REAL(q) :: CMASS(3,1),RM(3,1)
            REAL(q) :: pa(3),pb(3)
            REAL(q) :: complexcoord,dummyp,dummyq,dummyr,dummys,dummyt
            INTEGER :: ii,jj,kk

            DO i=1,SIZE(ICOORDINATES%COORDSTRUCT)
              SELECT CASE (ICOORDINATES%COORDSTRUCT(i)%TAG)
                CASE('X ')
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=DX(1,ICOORDINATES%COORDSTRUCT(i)%WHAT(1))
                CASE('Y ')
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=DX(2,ICOORDINATES%COORDSTRUCT(i)%WHAT(1))
                CASE('Z ')
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=DX(3,ICOORDINATES%COORDSTRUCT(i)%WHAT(1))
                CASE('cX')
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=&
                  & DX(1,ICOORDINATES%COORDSTRUCT(i)%WHAT(1))*A(1,1)+&
                  & DX(2,ICOORDINATES%COORDSTRUCT(i)%WHAT(1))*A(1,2)+&
                  & DX(3,ICOORDINATES%COORDSTRUCT(i)%WHAT(1))*A(1,3)
                CASE('cY')
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=&
                  & DX(1,ICOORDINATES%COORDSTRUCT(i)%WHAT(1))*A(2,1)+&
                  & DX(2,ICOORDINATES%COORDSTRUCT(i)%WHAT(1))*A(2,2)+&
                  & DX(3,ICOORDINATES%COORDSTRUCT(i)%WHAT(1))*A(2,3)
                CASE('cZ')
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=&
                  & DX(1,ICOORDINATES%COORDSTRUCT(i)%WHAT(1))*A(3,1)+&
                  & DX(2,ICOORDINATES%COORDSTRUCT(i)%WHAT(1))*A(3,2)+&
                  & DX(3,ICOORDINATES%COORDSTRUCT(i)%WHAT(1))*A(3,3)
                CASE('R ')
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=CAL_DISTANCE(T_INFO%NIONS,A,&
                  & DX,ICOORDINATES%COORDSTRUCT(i)%what(1),ICOORDINATES%COORDSTRUCT(i)%what(2),&
                  & ICOORDINATES%COORDSTRUCT(i)%WHERE(1,:),&
                  & ICOORDINATES%COORDSTRUCT(i)%WHERE(2,:))
                CASE('Q ')
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=CAL_DISTANCE(T_INFO%NIONS,A,&
                  & DX,ICOORDINATES%COORDSTRUCT(i)%what(1),ICOORDINATES%COORDSTRUCT(i)%what(2),&
                  & ICOORDINATES%COORDSTRUCT(i)%WHERE(1,:),&
                  & ICOORDINATES%COORDSTRUCT(i)%WHERE(2,:))
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=ICOORDINATES%COORDSTRUCT(i)%VALUE**2
                CASE('A ')
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=CAL_ANGLE(T_INFO%NIONS,A,DX,&
                  ICOORDINATES%COORDSTRUCT(i)%WHAT(1),&
                  ICOORDINATES%COORDSTRUCT(i)%WHAT(2),&
                  ICOORDINATES%COORDSTRUCT(i)%WHAT(3), &
                  ICOORDINATES%COORDSTRUCT(i)%WHERE(1,:),&
                  ICOORDINATES%COORDSTRUCT(i)%WHERE(2,:), &
                  ICOORDINATES%COORDSTRUCT(i)%WHERE(3,:),info)
                CASE('T ')
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=CAL_TORSION(T_INFO%NIONS,A,DX,&
                  ICOORDINATES%COORDSTRUCT(i)%WHAT(1),&
                  ICOORDINATES%COORDSTRUCT(i)%WHAT(2), &
                  ICOORDINATES%COORDSTRUCT(i)%WHAT(3),&
                  ICOORDINATES%COORDSTRUCT(i)%WHAT(4), &
                  ICOORDINATES%COORDSTRUCT(i)%WHERE(1,:),&
                  ICOORDINATES%COORDSTRUCT(i)%WHERE(2,:), &
                  ICOORDINATES%COORDSTRUCT(i)%WHERE(3,:),&
                  ICOORDINATES%COORDSTRUCT(i)%WHERE(4,:),info)
                CASE('M')
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=CAL_MIDDLE(T_INFO%NIONS,A,DX,&
                  ICOORDINATES%COORDSTRUCT(i)%WHAT(1),&
                  ICOORDINATES%COORDSTRUCT(i)%WHAT(2),&
                  ICOORDINATES%COORDSTRUCT(i)%WHAT(3), &
                  ICOORDINATES%COORDSTRUCT(i)%WHERE(1,:),&
                  ICOORDINATES%COORDSTRUCT(i)%WHERE(2,:), &
                  ICOORDINATES%COORDSTRUCT(i)%WHERE(3,:))
                CASE('B')
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=CAL_MIDDLE2(T_INFO%NIONS,A,DX,&
                  ICOORDINATES%COORDSTRUCT(i)%WHAT(1),&
                  ICOORDINATES%COORDSTRUCT(i)%WHAT(2),&
                  ICOORDINATES%COORDSTRUCT(i)%WHAT(3), &
                  ICOORDINATES%COORDSTRUCT(i)%WHAT(4), &
                  ICOORDINATES%COORDSTRUCT(i)%WHERE(1,:),&
                  ICOORDINATES%COORDSTRUCT(i)%WHERE(2,:),&
                  ICOORDINATES%COORDSTRUCT(i)%WHERE(3,:),&
                  ICOORDINATES%COORDSTRUCT(i)%WHERE(4,:))
                CASE('P ')
                  dummyp=CAL_DISTANCE(T_INFO%NIONS,A,&
                  & DX,ICOORDINATES%COORDSTRUCT(i)%what(1),ICOORDINATES%COORDSTRUCT(i)%what(2),&
                  & ICOORDINATES%COORDSTRUCT(i)%WHERE(1,:),&
                  & ICOORDINATES%COORDSTRUCT(i)%WHERE(2,:))
                  dummyq=CAL_DISTANCE(T_INFO%NIONS,A,&
                  & DX,ICOORDINATES%COORDSTRUCT(i)%what(3),ICOORDINATES%COORDSTRUCT(i)%what(4),&
                  & ICOORDINATES%COORDSTRUCT(i)%WHERE(3,:),&
                  & ICOORDINATES%COORDSTRUCT(i)%WHERE(4,:))
                  
                  IF (dummyq>0.) THEN
                    ICOORDINATES%COORDSTRUCT(i)%VALUE=dummyp/dummyq
                  ENDIF
                !CASE('E')
                !  ICOORDINATES%COORDSTRUCT(i)%VALUE=RM(1,1)
                !CASE('F')
                !  ICOORDINATES%COORDSTRUCT(i)%VALUE=RM(2,1)
                !CASE('G')
                !  ICOORDINATES%COORDSTRUCT(i)%VALUE=RM(3,1)
                CASE('W ')
                  dummyp=CAL_DISTANCE(T_INFO%NIONS,A,&
                  & DX,ICOORDINATES%COORDSTRUCT(i)%what(1),ICOORDINATES%COORDSTRUCT(i)%what(2),&
                  & ICOORDINATES%COORDSTRUCT(i)%WHERE(1,:),&
                  & ICOORDINATES%COORDSTRUCT(i)%WHERE(2,:))
                  dummyp=dummyp/ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(1)
                  dummyq=1.-dummyp**ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(2)
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=dummyq/(1.-dummyp**ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(3))
                CASE('WP')
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=0.
                  DO ii=-ICOORDINATES%COORDSTRUCT(i)%TRANS(1),ICOORDINATES%COORDSTRUCT(i)%TRANS(1)
                    DO jj=-ICOORDINATES%COORDSTRUCT(i)%TRANS(2),ICOORDINATES%COORDSTRUCT(i)%TRANS(2)
                      DO kk=-ICOORDINATES%COORDSTRUCT(i)%TRANS(3),ICOORDINATES%COORDSTRUCT(i)%TRANS(3)
                        dummyp=CAL_DISTANCE_trans(T_INFO%NIONS,A,&
                        &  DX,ICOORDINATES%COORDSTRUCT(i)%what(1),ICOORDINATES%COORDSTRUCT(i)%what(2), &
                        &  ICOORDINATES%COORDSTRUCT(i)%WHERE(1,:),ICOORDINATES%COORDSTRUCT(i)%WHERE(2,:),&
                        &  (/ii,jj,kk/))
                         dummyp=dummyp/ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(1)
                         dummyq=1.-dummyp**ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(2)
                         ICOORDINATES%COORDSTRUCT(i)%VALUE= ICOORDINATES%COORDSTRUCT(i)%VALUE + &
                         &              dummyq/(1.-dummyp**ICOORDINATES%COORDSTRUCT(i)%AUXCOEF(3))
                      ENDDO
                    ENDDO
                  ENDDO
                CASE('LR')
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=CAL_LR(A,ICOORDINATES%COORDSTRUCT(i)%WHAT(1))
                CASE('LA')
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=CAL_LA(A,ICOORDINATES%COORDSTRUCT(i)%WHAT(1),ICOORDINATES%COORDSTRUCT(i)%WHAT(2))
                CASE('LV')
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=CAL_LV(A)
                CASE('C ')
                  complexcoord=0
                  DO j=1,SIZE(ICOORDINATES%COORDSTRUCT(i)%COEFS)
                    IF (ABS(ICOORDINATES%COORDSTRUCT(i)%COEFS(j))>0.) THEN
                      complexcoord=complexcoord+(ICOORDINATES%COORDSTRUCT(i)%COEFS(j)**2*&
                      & ICOORDINATES%COORDSTRUCT(j)%VALUE**2)
                    ENDIF 
                  ENDDO
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=complexcoord**0.5
                CASE('S')
                  dummyp=0.
                  dummyq=0.
                  dummyr=0.
                  dummys=0.
                  dummyt=0.
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=0._q
                  DO j=1,SIZE(ICOORDINATES%COORDSTRUCT(i)%COEFS)
                    !c torsions require a special care
                    IF (ICOORDINATES%COORDSTRUCT(j)%TAG=='T ') THEN
                      IF (ICOORDINATES%COORDSTRUCT(i)%COEFS(j)>0.) THEN
                        CALL Min_Image_Cyclic(dummyp,ICOORDINATES%COORDSTRUCT(j)%VALUE)
                        dummyq=dummyq+ICOORDINATES%COORDSTRUCT(j)%VALUE*ICOORDINATES%COORDSTRUCT(i)%COEFS(j)
                        dummyp=ICOORDINATES%COORDSTRUCT(j)%VALUE
                        dummys=dummys+ICOORDINATES%COORDSTRUCT(i)%COEFS(j)
                      ELSE IF (ICOORDINATES%COORDSTRUCT(i)%COEFS(j)<0.) THEN
                        CALL Min_Image_Cyclic(dummyp,ICOORDINATES%COORDSTRUCT(j)%VALUE)
                        dummyr=dummyr+ICOORDINATES%COORDSTRUCT(j)%VALUE*ICOORDINATES%COORDSTRUCT(i)%COEFS(j)
                        dummyp=ICOORDINATES%COORDSTRUCT(j)%VALUE
                        dummyt=dummyt+(ICOORDINATES%COORDSTRUCT(i)%COEFS(j))
                      ENDIF
                    ELSE
                      IF (ABS(ICOORDINATES%COORDSTRUCT(i)%COEFS(j))>0.) THEN
                        ICOORDINATES%COORDSTRUCT(i)%VALUE=ICOORDINATES%COORDSTRUCT(i)%VALUE+&
                        &  ICOORDINATES%COORDSTRUCT(j)%VALUE*ICOORDINATES%COORDSTRUCT(i)%COEFS(j)
                      ENDIF
                    ENDIF
                  ENDDO
                  !c apply use minimum image along with CBC for torsions
                  IF (dummys > 0.) THEN
                    dummyq=dummyq/dummys
                    CALL Min_Image_Cyclic(0._q,dummyq)
                  ENDIF
                  IF (dummyt < 0.) THEN
                    dummyr=dummyr/dummyt
                    CALL Min_Image_Cyclic(dummyq,dummyr)
                  ENDIF
                  dummyq=dummyq*dummys
                  dummyr=dummyr*dummyt
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=ICOORDINATES%COORDSTRUCT(i)%VALUE+dummyq+dummyr
                CASE('ST')
                  dummyp=0.
                  dummyq=0.
                  dummyr=0.
                  dummys=0.
                  dummyt=0.
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=0._q
                  DO j=1,SIZE(ICOORDINATES%COORDSTRUCT(i)%COEFS)
                    !c torsions require a special care
                    IF (ICOORDINATES%COORDSTRUCT(j)%TAG=='T ') THEN
                      IF (ICOORDINATES%COORDSTRUCT(i)%COEFS(j)>0.) THEN
                        CALL Min_Image_Cyclic(dummyp,ICOORDINATES%COORDSTRUCT(j)%VALUE)
                        dummyq=dummyq+ICOORDINATES%COORDSTRUCT(j)%VALUE*ICOORDINATES%COORDSTRUCT(i)%COEFS(j)
                        dummyp=ICOORDINATES%COORDSTRUCT(j)%VALUE
                        dummys=dummys+ICOORDINATES%COORDSTRUCT(i)%COEFS(j)
                      ELSE IF (ICOORDINATES%COORDSTRUCT(i)%COEFS(j)<0.) THEN
                        CALL Min_Image_Cyclic(dummyp,ICOORDINATES%COORDSTRUCT(j)%VALUE)
                        dummyr=dummyr+ICOORDINATES%COORDSTRUCT(j)%VALUE*ICOORDINATES%COORDSTRUCT(i)%COEFS(j)
                        dummyp=ICOORDINATES%COORDSTRUCT(j)%VALUE
                        dummyt=dummyt+(ICOORDINATES%COORDSTRUCT(i)%COEFS(j))
                      ENDIF
                    ENDIF
                  ENDDO
                  !c apply use minimum image along with CBC for torsions
                  IF (dummys > 0.) THEN
                    dummyq=dummyq/dummys
                    CALL Min_Image_Cyclic(0._q,dummyq)
                  ENDIF
                  IF (dummyt < 0.) THEN
                    dummyr=dummyr/dummyt
                    CALL Min_Image_Cyclic(dummyq,dummyr)
                  ENDIF
                  dummyq=dummyq*dummys
                  dummyr=dummyr*dummyt
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=dummyq+dummyr
                CASE('IS')
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=0._q
                  dummyp=0.
                  dummyq=0.
                  DO k=1,ICOORDINATES%IRCDIM(2)
                    dummyr=0.
                    DO j=1,ICOORDINATES%NUMPRIMITIVES
                      IF (ICOORDINATES%COORDSTRUCT(i)%COEFS(j)>0.) THEN
                        !dummyr=dummyr+ICOORDINATES%COORDSTRUCT(i)%COEFS(j)**2 * &
                        !&   (ICOORDINATES%COORDSTRUCT(j)%VALUE-ICOORDINATES%IRCPOINTS(j,k))**2
                        dummyr=dummyr+ICOORDINATES%COORDSTRUCT(i)%COEFS(j) * &
                        &   (ICOORDINATES%COORDSTRUCT(j)%VALUE-ICOORDINATES%IRCPOINTS(j,k))**2
                      ENDIF
                    ENDDO
                    !dummyr=dummyr**0.5
                    dummyr=exp(-dummyr)
                    dummyp=dummyp+(k-1)*dummyr
                    dummyq=dummyq+dummyr
                  ENDDO  
                  IF (dummyq>0.) ICOORDINATES%COORDSTRUCT(i)%VALUE=dummyp/(dummyq*(ICOORDINATES%IRCDIM(2)-1))
                CASE('IZ')
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=0._q
                  dummyq=0.
                  dummyp=0.
                  DO k=1,ICOORDINATES%IRCDIM(2)
                    dummyr=0.
                    DO j=1,ICOORDINATES%NUMPRIMITIVES
                      IF (ICOORDINATES%COORDSTRUCT(i)%COEFS(j)>0.) THEN
                        dummyp=ICOORDINATES%COORDSTRUCT(i)%COEFS(j)
                        dummyr=dummyr+ICOORDINATES%COORDSTRUCT(i)%COEFS(j) * &
                        &   (ICOORDINATES%COORDSTRUCT(j)%VALUE-ICOORDINATES%IRCPOINTS(j,k))**2
                      ENDIF
                    ENDDO
                    dummyr=exp(-dummyr)
                    dummyq=dummyq+dummyr
                  ENDDO  
                  IF (dummyq>0.) ICOORDINATES%COORDSTRUCT(i)%VALUE=-LOG(dummyq)/dummyp  
                CASE('D')
                  ICOORDINATES%COORDSTRUCT(i)%VALUE=0._q
                  DO j=1,SIZE(ICOORDINATES%COORDSTRUCT(i)%COEFS)
                    IF (ABS(ICOORDINATES%COORDSTRUCT(i)%COEFS(j))>1e-4) THEN
                      dummyq=ICOORDINATES%COORDSTRUCT(j)%VALUE/ABS(ICOORDINATES%COORDSTRUCT(i)%COEFS(j))
                      IF (ABS(dummyq-1._q)<1e-4) THEN
                        dummyq=1._q+1e-4
                      ENDIF
                      IF (ICOORDINATES%COORDSTRUCT(i)%COEFS(j)>0._q) THEN
                        ICOORDINATES%COORDSTRUCT(i)%VALUE=ICOORDINATES%COORDSTRUCT(i)%VALUE+&
                        & (1.0-(dummyq)**ICOORDINATES%CNEXP(1))/(1.0-(dummyq)**ICOORDINATES%CNEXP(2))
                      ELSE
                        ICOORDINATES%COORDSTRUCT(i)%VALUE=ICOORDINATES%COORDSTRUCT(i)%VALUE-&
                        & (1.0-(dummyq)**ICOORDINATES%CNEXP(1))/(1.0-(dummyq)**ICOORDINATES%CNEXP(2))
                      ENDIF
                    ENDIF
                  ENDDO
              END SELECT
            ENDDO                       
          END SUBROUTINE DEAL_XYZ

          SUBROUTINE MAKE_UMAT(BMAT,UTRANS)
            !c prepare U matrix which defines
            !c delocalised intern
            REAL(q),DIMENSION(:,:) :: BMAT,UTRANS
            REAL(q),ALLOCATABLE :: GMAT(:,:)
            INTEGER :: shp(2)
            INTEGER :: i,j            
            REAL(q),ALLOCATABLE :: d(:)         !c eigenvalues of the matrix        
            REAL(q),ALLOCATABLE :: u(:,:)     !c eigenvectors of the matrix
            REAL(q),ALLOCATABLE :: vt(:,:)

            shp=SHAPE(BMAT)
            ALLOCATE(GMAT(shp(2),shp(2)))
            GMAT=MATMUL(TRANSPOSE(BMAT),BMAT)
            ALLOCATE(d(shp(2)),u(shp(2),shp(2)),vt(shp(2),shp(2)))
            CALL SVDVALVEC(GMAT,shp(2),d,u,vt) !c calculates evectors and evalues
            DEALLOCATE(GMAT)
            ALLOCATE(GMAT(SIZE(UTRANS(:,1)),shp(2)))
            GMAT=0._q
            j=0
            DO i=1,shp(2)
              IF (ABS(d(i))>1e-06) THEN 
                j=j+1 
                GMAT(j,:)=u(:,i)/(ABS(d(i)))**0.5
              ENDIF
            ENDDO
            UTRANS=MATMUL(GMAT,TRANSPOSE(BMAT))
            DEALLOCATE(vt,u,d,GMAT)
          END SUBROUTINE
      END MODULE
