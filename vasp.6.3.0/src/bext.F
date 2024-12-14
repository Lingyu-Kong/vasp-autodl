#include "symbol.inc"
      MODULE bexternal
      USE prec
      IMPLICIT NONE

      PUBLIC BEXT_READER,LBEXTERNAL,BEXT,BEXT_ADDV

      PRIVATE

      REAL(q), SAVE :: BEXT(3)=0

      LOGICAL, SAVE :: LBEXT=.FALSE.

      CONTAINS

!***********************************************************************
!******************** PUBLIC PROCEDURES ********************************
!***********************************************************************

!******************** SUBROUTINE BEXT_READER ***************************
!
! Reads BEXT from the INCAR file
!
!***********************************************************************

      SUBROUTINE BEXT_READER(IU0,IU5)
      USE base
      USE vaspxml
      USE reader_tags

      INTEGER       :: IU5,IU6,IU0
      ! local variables
      INTEGER       :: IDUM, N, IERR
      REAL(q)       :: RDUM
      COMPLEX(q)    :: CDUM
      LOGICAL       :: LOPEN,LDUM
      CHARACTER (1) :: CHARAC

      REAL(q)       :: BNORM2
      INTEGER       :: ISPIN,NREQ
      LOGICAL       :: LSORBIT,LNONCOLLINEAR

      CALL OPEN_INCAR_IF_FOUND(IU5, LOPEN)
! read in flag ISPIN
      ISPIN=1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ISPIN', ISPIN, IERR, WRITEXMLINCAR)
! read in flag LSORBIT
      LSORBIT=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSORBIT', LSORBIT, IERR, WRITEXMLINCAR)
! read in flag LNONCOLLINEAR
      LNONCOLLINEAR=LSORBIT
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LNONCOLLINEAR', LNONCOLLINEAR, IERR, WRITEXMLINCAR)
      IF (LSORBIT) LNONCOLLINEAR=LSORBIT

      NREQ=0
      IF (LNONCOLLINEAR) THEN
         NREQ=3
      ELSEIF (ISPIN==2) THEN
         NREQ=1
      ENDIF

      IF (NREQ/=0) THEN
! read in flag BEXT
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'BEXT', BEXT(1:NREQ), NREQ, IERR, WRITEXMLINCAR)
         BNORM2=BEXT(1)*BEXT(1)+BEXT(2)*BEXT(2)+BEXT(3)*BEXT(3)
         IF (BNORM2>1.E-6_q) THEN
            LBEXT=.TRUE.
         ENDIF
      ENDIF

      CALL CLOSE_INCAR_IF_FOUND(IU5)
      RETURN
      END SUBROUTINE BEXT_READER


!******************** FUNCTION LBEXTERNAL ******************************
!
!***********************************************************************

      FUNCTION LBEXTERNAL()
      LOGICAL LBEXTERNAL
      LBEXTERNAL=LBEXT
      END FUNCTION LBEXTERNAL


!******************** SUBROUTINE BEXT_ADDV *****************************
!
!***********************************************************************

      SUBROUTINE BEXT_ADDV(CVTOT,GRIDC,NCDIJ)
      USE mgrid
      
      TYPE (grid_3d) GRIDC

      COMPLEX(q)  :: CVTOT(GRIDC%MPLWV,NCDIJ)
      INTEGER     :: NCDIJ
      ! local variables
      INTEGER     :: I

      DO I=2,NCDIJ
         CALL ADD2VG0(CVTOT(1,I),GRIDC,CMPLX(BEXT(I-1),0._q,KIND=q))
      ENDDO

      RETURN
      END SUBROUTINE BEXT_ADDV

      END MODULE bexternal


!******************** SUBROUTINE ADD2VG0 *******************************
!
!***********************************************************************

      SUBROUTINE ADD2VG0(CVTOT,GRIDC,C)
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE mgrid
      TYPE (grid_3d) GRIDC
      COMPLEX(q)  :: CVTOT(GRIDC%RC%NROW,GRIDC%RC%NCOL)
      COMPLEX(q)  :: C
      ! local variables
      INTEGER N1,N2,N3,NC

      N1=1; N2=1; N3=1

!$ACC PARALLEL LOOP PRESENT(CVTOT) __IF_ASYNC__
      DO NC=1,GRIDC%RC%NCOL
         IF (GRIDC%RC%I2(NC)==N2 .AND. GRIDC%RC%I3(NC)==N3) THEN
            CVTOT(NC,N1)=CVTOT(NC,N1)+C
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE ADD2VG0
