#if defined(MPI) || defined(MPI_CHAIN)
#include "symbol.inc"
!****************************************************************************************************
! Module containing subroutines and functions necessary for output information on display.
!****************************************************************************************************

      MODULE LOGFILE
        USE ML_FF_CONSTANT
        USE MPI_DATA
        USE ML_FF_PREC
        USE ML_FF_STRUCT, ONLY: ML_IO_WRITE

        IMPLICIT NONE

        !> Width of first column keyword.
        INTEGER, PARAMETER :: C1W = 15

        PRIVATE :: STRLJ

        PRIVATE :: PRINT_TAG, PRINT_TAG_REAL, PRINT_TAG_INTEGER, PRINT_TAG_LOGICAL

        INTERFACE PRINT_TAG
           PROCEDURE PRINT_TAG_REAL
           PROCEDURE PRINT_TAG_INTEGER
           PROCEDURE PRINT_TAG_LOGICAL
        END INTERFACE PRINT_TAG

        CONTAINS

         SUBROUTINE PRINT_TAG_REAL(IU, TAG_LIST, TAG, VAL)
            USE ML_FF_TAGLIST
            INTEGER,          INTENT(IN) :: IU
            TYPE(TAGLIST),    INTENT(IN) :: TAG_LIST
            CHARACTER(LEN=*), INTENT(IN) :: TAG
            REAL(q),          INTENT(IN) :: VAL

            WRITE(IU, '(A," : ",ES13.5,1X,A,1X,A)') &
               TAGLIST_DESCRIPTION(TAG),            &
               VAL,                                 &
               TAG_LIST%SYMBOL(TAG),                &
               TAG

         END SUBROUTINE PRINT_TAG_REAL

         SUBROUTINE PRINT_TAG_INTEGER(IU, TAG_LIST, TAG, VAL)
            USE ML_FF_TAGLIST
            INTEGER,          INTENT(IN) :: IU
            TYPE(TAGLIST),    INTENT(IN) :: TAG_LIST
            CHARACTER(LEN=*), INTENT(IN) :: TAG
            INTEGER,          INTENT(IN) :: VAL

            WRITE(IU, '(A," : ",I13,1X,A,1X,A)') &
               TAGLIST_DESCRIPTION(TAG),         &
               VAL,                              &
               TAG_LIST%SYMBOL(TAG),             &
               TAG

         END SUBROUTINE PRINT_TAG_INTEGER

         SUBROUTINE PRINT_TAG_LOGICAL(IU, TAG_LIST, TAG, VAL)
            USE ML_FF_TAGLIST
            INTEGER,          INTENT(IN) :: IU
            TYPE(TAGLIST),    INTENT(IN) :: TAG_LIST
            CHARACTER(LEN=*), INTENT(IN) :: TAG
            LOGICAL,          INTENT(IN) :: VAL

            WRITE(IU, '(A," : ",L13,1X,A,1X,A)') &
               TAGLIST_DESCRIPTION(TAG),         &
               VAL,                              &
               TAG_LIST%SYMBOL(TAG),             &
               TAG

         END SUBROUTINE PRINT_TAG_LOGICAL

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "STATUS" line (+ preceding separator line).
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_STATUS(FF, LOGF, LDO_AB_INITIO)

            USE ML_FF_STRUCT, ONLY: FF_PAR, ML_MPI_PAR, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),     INTENT(IN) :: FF !< FF setup.
            TYPE(LOG_PAR),    INTENT(IN) :: LOGF !< Log file setup.
            LOGICAL,          INTENT(IN) :: LDO_AB_INITIO !< Ab initio T/F?

            ENUM, BIND(C)
               ENUMERATOR :: UNKNOWN    = -1
               ENUMERATOR :: STATE_ENUM = 0
               ENUMERATOR :: ACCURATE   = 1
               ENUMERATOR :: THRESHOLD  = 2
               ENUMERATOR :: LEARNING   = 3
               ENUMERATOR :: CRITICAL   = 4
               ENUMERATOR :: PREDICT    = 5
            END ENUM

            CHARACTER(LEN=*), PARAMETER :: STATE_STRINGS(UNKNOWN:PREDICT) = [&
               "unknown  ", &
               "?????????", &
               "accurate ", &
               "threshold", &
               "learning ", &
               "critical ", &
               "predict  " ]

            INTEGER(KIND(STATE_ENUM)) :: STATE !< State descriptor.

            STATE = STATE_ENUM

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            ! Separator line.
            WRITE(LOGF%IU, '(80("-"))')

            ! Decide current state.
            IF (FF%LMLONLY) THEN
               STATE = PREDICT
            ELSE
               IF (LDO_AB_INITIO) THEN
                  IF (FF%LDOUB) THEN
                     STATE = CRITICAL
                  ELSE
                     IF (FF%LGENFF) THEN
                        STATE = LEARNING
                     ELSE
                        IF (FF%LSAMPLE) THEN
                           STATE = THRESHOLD
                        ELSE
                           STATE = UNKNOWN
                        END IF
                     END IF
                  END IF
               ELSE
                  STATE = ACCURATE
               END IF
            END IF

            ! STATUS line.
            WRITE(LOGF%IU,                                                    &
                  '(A,1X,I9,1X,A,1X,I2,1X,L6,1X,L6,1X,I9,1X,I9)') &
               STRLJ("STATUS", C1W),                                          &
               FF%NSTEP,                                                      &
               STATE_STRINGS(STATE),                                          &
               STATE,                                                         &
               LDO_AB_INITIO,                                                 &
               FF%IFF,                                                        &
               FF%NSTEP_FROM_SAMPLING,                                        &
               FF%NSTEP_FROM_GENFF

         END SUBROUTINE LOGLINE_STATUS

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "LCONF" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_LCONF(FF, FFM, LOGF)

            USE ML_FF_STRUCT, ONLY: FF_PAR, FFM_PAR, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),  INTENT(IN) :: FF !< FF setup.
            TYPE(FFM_PAR), INTENT(IN) :: FFM !< FFM setup.
            TYPE(LOG_PAR), INTENT(IN) :: LOGF !< Log file setup.

            INTEGER :: I

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            ! LCONF line.
            WRITE(LOGF%IU, '(A,1X,I9)', ADVANCE='NO') &
               STRLJ("LCONF", C1W),                   &
               FF%NSTEP
            DO I = 1, FF%MTYP
               WRITE(LOGF%IU, '(1X,A2,2(1X,I9))', ADVANCE='NO') &
                  FF%TYPE(FF%LADD(I)),                                   &
                  FFM%NB(FF%LADD(I)),                                    &
                  FFM%NB(FF%LADD(I)) + FFM%NB_NEW(FF%LADD(I))
            END DO
            WRITE(LOGF%IU, *)

         END SUBROUTINE LOGLINE_LCONF

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "STDAB" and "ERR" lines.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_STDAB_ERR(FF, LOGF)

            USE ML_FF_STRUCT, ONLY: FF_PAR, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),     INTENT(IN) :: FF !< FF setup.
            TYPE(LOG_PAR),    INTENT(IN) :: LOGF !< Log file setup.

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            ! STDAB line.
            ASSOCIATE(S => FF%SIG)
            WRITE(LOGF%IU, '(A,1X,I9,3(1X,ES16.8))')                &
               STRLJ("STDAB", C1W),                                 &
               FF%NSTEP,                                            &
               S(1) * EUNIT,                                        &
               SQRT((S(2)**2 + S(3)**2 + S(4)**2) / 3.0_q) * FUNIT, &
               SQRT((S(5)**2 + S(6)**2 + S(7)**2 +                  &
                     S(8)**2 + S(9)**2 + S(10)**2) / 6.0_q) * SUNIT
            END ASSOCIATE
            WRITE(LOGF%IU, '(A,1X,I9,3(1X,ES16.8))') &
               STRLJ("ERR", C1W),                    &
               FF%NSTEP,                             &
               FF%STOTEN * EUNIT,                    &
               FF%STIFOR * FUNIT,                    &
               FF%STSIF * SUNIT

         END SUBROUTINE LOGLINE_STDAB_ERR

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "LASTE"lines.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_LAST_ERR(FF, LOGF)
            USE ML_FF_STRUCT, ONLY: FF_PAR, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),     INTENT(IN) :: FF !< FF setup.
            TYPE(LOG_PAR),    INTENT(IN) :: LOGF !< Log file setup.

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            WRITE(LOGF%IU, '(A,1X,I9,3(1X,ES16.8))') &
               STRLJ("LASTE", C1W),                    &
               FF%NSTEP,                             &
               FF%STOTEN_LAST * EUNIT,                    &
               FF%STIFOR_LAST * FUNIT,                    &
               FF%STSIF_LAST  * SUNIT
           
         END SUBROUTINE LOGLINE_LAST_ERR

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "BEE" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_BEE(FF, LOGF, ABN, FINAL_CALL)

            USE ML_FF_STRUCT, ONLY: FF_PAR, LOG_PAR, ABINITIONEW
            IMPLICIT NONE

            TYPE(FF_PAR),      INTENT(IN) :: FF !< FF setup.
            TYPE(LOG_PAR),     INTENT(IN) :: LOGF !< Log file setup.
            TYPE(ABINITIONEW), INTENT(IN) :: ABN !< New ab initio data.
            !> Whether this is the final call in this time step.
            LOGICAL,           INTENT(IN) :: FINAL_CALL

            CHARACTER(LEN=:), ALLOCATABLE :: KEYWORD

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            KEYWORD = "BEE"
            IF (FINAL_CALL) KEYWORD = "BEEF"

            WRITE(LOGF%IU, '(A,1X,I9,6(1X,ES16.8))') &
               STRLJ(KEYWORD, C1W),                  &
               FF%NSTEP,                             &
               ABN%BEMAX(1, ABN%NCONF) * EUNIT,      &
               ABN%BEMAX(2, ABN%NCONF) * FUNIT,      &
               ABN%BEAV(2, ABN%NCONF) * FUNIT,       &
               FF%CTIFOR * FUNIT,                    &
               ABN%BEMAX(3, ABN%NCONF) * SUNIT,      &
               ABN%BEAV(3, ABN%NCONF) * SUNIT

         END SUBROUTINE LOGLINE_BEE

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "CFE" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_CFE(LOGF, NSTEP, FRMSE, FMAE, FMAX)

            USE ML_FF_STRUCT, ONLY: LOG_PAR
            IMPLICIT NONE

            TYPE(LOG_PAR),    INTENT(IN) :: LOGF !< Log file setup.
            INTEGER,          INTENT(IN) :: NSTEP !< Current time step.
            REAL(q),          INTENT(IN) :: FRMSE !< RMSE of forces.
            REAL(q),          INTENT(IN) :: FMAE !< MAE of forces.
            REAL(q),          INTENT(IN) :: FMAX !< Maximum error of forces.

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            ! CFE line.
            WRITE(LOGF%IU, '(A,1X,I9,3(1X,ES16.8))') &
               STRLJ("CFE", C1W),                    &
               NSTEP,                                &
               FRMSE * FUNIT,                        &
               FMAE * FUNIT,                         &
               FMAX * FUNIT

         END SUBROUTINE LOGLINE_CFE

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "THRHIST" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_THRHIST(FF, LOGF)

            USE ML_FF_STRUCT, ONLY: FF_PAR, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),     INTENT(IN) :: FF !< FF setup.
            TYPE(LOG_PAR),    INTENT(IN) :: LOGF !< Log file setup.

            INTEGER :: I

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            ! THRHIST line.
            DO I = 1, FF%NHIS
               WRITE(LOGF%IU, '(A,1X,I9,1X,I4,1X,ES16.8)') &
                  STRLJ("THRHIST", C1W),                   &
                  FF%NSTEP,                                &
                  I,                                       &
                  FF%BEHIS(I) * FUNIT
            END DO

         END SUBROUTINE LOGLINE_THRHIST

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "THRUPD" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_THRUPD(FF, LOGF)

            USE ML_FF_STRUCT, ONLY: FF_PAR, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),     INTENT(IN) :: FF !< FF setup.
            TYPE(LOG_PAR),    INTENT(IN) :: LOGF !< Log file setup.

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            ! THRUPD line.
            WRITE(LOGF%IU, '(A,1X,I9,4(1X,ES16.8))') &
               STRLJ("THRUPD", C1W),                 &
               FF%NSTEP,                             &
               FF%CTIFOR_PREV * FUNIT,               &
               FF%CTIFOR * FUNIT,                    &
               FF%CSIG_LOG,                          &
               FF%CSLOPE_LOG

         END SUBROUTINE LOGLINE_THRUPD

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "SPRSC" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_SPRSC(FF, FFM, LOGF, NB0, NCONF0, NCONF)

            USE ML_FF_STRUCT, ONLY: FF_PAR, FFM_PAR, LOG_PAR
            IMPLICIT NONE

            TYPE(FF_PAR),  INTENT(IN) :: FF !< FF setup.
            TYPE(FFM_PAR), INTENT(IN) :: FFM !< FFM setup.
            TYPE(LOG_PAR), INTENT(IN) :: LOGF !< Log file setup.
            INTEGER,       INTENT(IN) :: NB0(:) !< LRC before spars.
            INTEGER,       INTENT(IN) :: NCONF0 !< Conf. before spars.
            INTEGER,       INTENT(IN) :: NCONF !< Conf. after spars.

            INTEGER :: I

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            ! SPRSC line.
            WRITE(LOGF%IU, '(A,1X,I9,1x,I9,1x,I9)', ADVANCE='NO') &
               STRLJ("SPRSC", C1W),                   &
               FF%NSTEP,                              &
               NCONF0,                                &
               NCONF
            DO I = 1, FF%MTYP
               WRITE(LOGF%IU, '(1X,A2,2(1X,I9))', ADVANCE='NO') &
                  FF%TYPE(I),                                   &
                  NB0(I),                                       &
                  FFM%NB(I)
            END DO
            WRITE(LOGF%IU, *)

         END SUBROUTINE LOGLINE_SPRSC

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "REGR" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_REGR(LOGF,     &
                                 NSTEP,    &
                                 LIN_REG,  &
                                 ITE,      &
                                 SIGW2,    &
                                 SIGV2,    &
                                 REGRATIO, &
                                 WRATIO,   &
                                 WMAX,     &
                                 ZTR,      &
                                 LEXIT)

            USE ML_FF_STRUCT, ONLY: LOG_PAR
            IMPLICIT NONE

            TYPE(LOG_PAR), INTENT(IN) :: LOGF !< Log file setup.
            INTEGER,       INTENT(IN) :: NSTEP !< Current time step.
            INTEGER,       INTENT(IN) :: LIN_REG !< BLR, QR or SVD
            INTEGER,       INTENT(IN) :: ITE !< Regression loop index.
            REAL(q),       INTENT(IN) :: SIGW2 !< sigma_W^2.
            REAL(q),       INTENT(IN) :: SIGV2 !< sigma_V^2.
            REAL(q),       INTENT(IN) :: REGRATIO !< (SIGW2 / SIGV2) / WMAX.
            REAL(q),       INTENT(IN) :: ZTR !< Error of force field
            REAL(q),       INTENT(IN) :: WRATIO !< WMAX / WMIN.
            REAL(q),       INTENT(IN) :: WMAX !< highest eigenvalue/singular value
            LOGICAL,       INTENT(IN) :: LEXIT !< If last loop iteration.

            CHARACTER(LEN=C1W) :: COL1

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            ! Last iteration uses slightly different keyword ("F" for final).
            COL1 = "REGR"
            IF (LEXIT) COL1 = "REGRF"

            ! REGR line.
            WRITE(LOGF%IU, '(A,1X,I9,1x,I4,1X,I4,5(1X,ES16.8))', ADVANCE='NO') &
               STRLJ(COL1, C1W),                                         &
               NSTEP,                                                    &
               LIN_REG,                                                  &
               ITE,                                                      &
               SIGW2,                                                    &
               SIGV2,                                                    &
               REGRATIO,                                                 &
               ZTR

            ! Write WRATIO and WMAX only in last line because it does not change over
            ! iterations.
            IF (LEXIT) WRITE(LOGF%IU, '(1X,ES16.8)', ADVANCE='NO') WRATIO
            IF (LEXIT) WRITE(LOGF%IU, '(1X,ES16.8)', ADVANCE='NO') WMAX

            ! Newline.
            WRITE(LOGF%IU, '()')

         END SUBROUTINE LOGLINE_REGR

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file output: "DCOUPLE" line.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGLINE_DCOUPLE(LOGF, NSTEP, DHDLAMBDA)

            USE ML_FF_STRUCT, ONLY: LOG_PAR
            IMPLICIT NONE

            TYPE(LOG_PAR),    INTENT(IN) :: LOGF !< Log file setup.
            INTEGER,          INTENT(IN) :: NSTEP !< Current time step.
            REAL(q),          INTENT(IN) :: DHDLAMBDA !< dH/dlambda.

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            ! DCOUPLE line.
            WRITE(LOGF%IU, '(A,1X,I9,1X,ES16.8)') &
               STRLJ("DCOUPLE", C1W),             &
               NSTEP,                             &
               DHDLAMBDA * EUNIT

         END SUBROUTINE LOGLINE_DCOUPLE

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file section about input parameters.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGSECTION_ML_FF_PARAMETERS(LOGF,FF,FFM)
            USE ML_FF_CONSTANT
            USE ML_FF_STRUCT, ONLY: LOG_PAR, FF_PAR, FFM_PAR
            IMPLICIT NONE

            ! Input variables
            TYPE(LOG_PAR), INTENT(IN)    :: LOGF !< Log file setup.
            TYPE(FF_PAR),  INTENT(IN)    :: FF !< FF setup.
            TYPE(FFM_PAR), INTENT(IN)    :: FFM !< FFM setup.
            ! Local variables
            INTEGER, PARAMETER           :: NWRITE = 3
            INTEGER                      :: ICOL
            INTEGER                      :: IROW
            INTEGER                      :: NCOL
            INTEGER                      :: NMOD
            CHARACTER(LEN=12)            :: STRING_HELP(1:3)

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            NCOL=FF%MTYP/NWRITE
            NMOD=MOD(FF%MTYP,NWRITE)

            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') "* PARAMETERS FROM ML_FF FILE **************&
               &**************************************************************&
               &***********************************"
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU,'(A)') "ML_FF format typ"
            IF (FF%ML_FF_BINARY) THEN
               WRITE(LOGF%IU,'(1x,A)') "binary"
            ELSE
               WRITE(LOGF%IU,'(1x,A)') "ASCII"
            ENDIF
            WRITE(LOGF%IU,'(A)') "Version number of ML_FF file"
            WRITE(STRING_HELP(1),'(I0)') FF%ML_FF_VERSION_READ(1)
            WRITE(STRING_HELP(2),'(I0)') FF%ML_FF_VERSION_READ(2)
            WRITE(STRING_HELP(3),'(I0)') FF%ML_FF_VERSION_READ(3)
            WrITE(LOGF%IU,'(1x,A)') ADJUSTL(TRIM(STRING_HELP(1))) // "." // &
                                 ADJUSTL(TRIM(STRING_HELP(2))) // "." // &
                                 ADJUSTL(TRIM(STRING_HELP(3))) 
            WRITE(LOGF%IU,'(A)') "Number of element types:"
            WRITE(LOGF%IU,'(1x,I6)') FF%MTYP
            WRITE(LOGF%IU,'(A)') "Element types: "
            DO ICOL=1, NCOL, 1
               WRITE(LOGF%IU,'(1x,A2,1x,A2,1x,A2)') (FF%TYPE(NWRITE*(ICOL-1)+IROW),IROW=1, NWRITE)
            ENDDO
            IF(NMOD .NE. 0) THEN
               WRITE(LOGF%IU,'(1x,A2,1x,A2,1x,A2)') (FF%TYPE(NWRITE*NCOL+IROW),IROW=1, NMOD)
            ENDIF
            WRITE(LOGF%IU,'(A)') "Reference atomic energies:"
            DO ICOL=1, NCOL, 1
               WRITE(LOGF%IU,'(1x,ES16.8,1x,ES16.8,1x,ES16.8)') (FF%EATOM_REF(NWRITE*(ICOL-1)+IROW)*EUNIT,IROW=1, NWRITE)
            ENDDO
            IF(NMOD .NE. 0) THEN
               WRITE(LOGF%IU,'(1x,ES16.8,1x,ES16.8,1x,ES16.8)') (FF%EATOM_REF(NWRITE*NCOL+IROW)*EUNIT,IROW=1, NMOD)
            ENDIF
            WRITE(LOGF%IU,'(A)') "Mean square root error in energy (eV atom^-1)"
            WRITE(LOGF%IU,'(1x,ES16.8)') FF%STOTEN*EUNIT
            WRITE(LOGF%IU,'(A)') "Mean square root error in forces (eV Angst^-1)"
            WRITE(LOGF%IU,'(1x,ES16.8)') FF%STIFOR/FUNIT
            WRITE(LOGF%IU,'(A)') "Mean square root error in stress (kB)"
            WRITE(LOGF%IU,'(1x,ES16.8)') FF%STSIF/SUNIT
            WRITE(LOGF%IU,'(A)') "The number of local reference configurations"
            DO ICOL=1, FF%MTYP
               WRITE(LOGF%IU,'(1x,A2,1x,I9)') FF%TYPE(ICOL),FFM%NB(ICOL)
            ENDDO
            IF (FFM%W2.GT.0.0_q) THEN
               WRITE(LOGF%IU,'(A)') "Number of descriptors"
               DO ICOL=1, FF%MTYP
                  WRITE(LOGF%IU,'(1x,A2,1x,I9)') FF%TYPE(ICOL),FFM%NNVAR2(ICOL)
               ENDDO
               WRITE(LOGF%IU,'(A)') "Number of descriptors for SIC"
               DO ICOL=1, FF%MTYP
                  WRITE(LOGF%IU,'(1x,A2,1x,I9)') FF%TYPE(ICOL),FFM%NNVAR2_SIC(ICOL)
               ENDDO
            ENDIF
         END SUBROUTINE LOGSECTION_ML_FF_PARAMETERS

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file section about input parameters.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGSECTION_PARAMETERS(TAG_LIST, LOGF, FF, FFM)

            USE ML_FF_CONSTANT
            USE ML_FF_TAGLIST
            USE ML_FF_STRUCT, ONLY: LOG_PAR, FF_PAR, FFM_PAR
            IMPLICIT NONE

            TYPE(TAGLIST), INTENT(INOUT) :: TAG_LIST !< List with tag occurrences.
            TYPE(LOG_PAR), INTENT(IN)    :: LOGF !< Log file setup.
            TYPE(FF_PAR),  INTENT(IN)    :: FF !< FF setup.
            TYPE(FFM_PAR), INTENT(IN)    :: FFM !< FFM setup.

            CHARACTER(LEN=42)            :: CURRENT_TAG
            INTEGER                      :: I
            INTEGER                      :: NUM_TAGS

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') "* MACHINE LEARNING SETTINGS ***************&
               &**************************************************************&
               &***********************************"
            WRITE(LOGF%IU, *)

            WRITE(LOGF%IU, '(A)') 'This section lists the available &
               &machine-learning related settings with a short description, &
               &their'
            WRITE(LOGF%IU, '(A)') 'selected values and the INCAR tags. The &
               &column between the value and the INCAR tag may contain a'
            WRITE(LOGF%IU, '(A)') '"state indicator" highlighting the origin &
               &of the value. Here is a list of possible indicators:'
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(" * ",A3," : ",A)')                          &
               TAGLIST_STATE_SYMBOL(TL_DEFAULT),                          &
               "(empty) Tag was not provided in the INCAR file, a default &
               &value was chosen automatically."
            WRITE(LOGF%IU, '(" * ",A3," : ",A)') &
               TAGLIST_STATE_SYMBOL(TL_INCAR),   &
               "Value was provided in the INCAR file."
            WRITE(LOGF%IU, '(" * ",A3," : ",A)')   &
               TAGLIST_STATE_SYMBOL(TL_INCAR_ALT), &
               "Value was provided in the INCAR file, deprecated tag."
            WRITE(LOGF%IU, '(" * ",A3," : ",A)')                       &
               TAGLIST_STATE_SYMBOL(TL_OVERRIDE),                      &
               "A value found in the INCAR file was overwritten by the &
               &contents of the ML_FF or ML_AB file."
            WRITE(LOGF%IU, '(" * ",A3," : ",A)')                            &
               TAGLIST_STATE_SYMBOL(TL_UNSET),                              &
               "The value for this tag was never set (please report this to &
               &the VASP developers)."
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') 'Tag values with associated units are given &
               &here in Angstrom/eV, if not specified otherwise.'
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') 'Please refer to the VASP online manual for &
               &a detailed description of available INCAR tags.'

            CALL SECTION_HEADER("General settings")
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_ISTART               ", FF%ISTART)

            CALL SECTION_HEADER("Descriptor settings")
            WRITE(LOGF%IU, '(A)') "Radial descriptors:"
            WRITE(LOGF%IU, '(A)') "-------------------"
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_RCUT1                ", FFM%RCUT1 * AUTOA)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_SION1                ", FFM%SION1 * AUTOA)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_MRB1                 ", FFM%MRB1)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') "Angular descriptors:"
            WRITE(LOGF%IU, '(A)') "--------------------"
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_RCUT2                ", FFM%RCUT2 * AUTOA)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_SION2                ", FFM%SION2 * AUTOA)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_MRB2                 ", FFM%MRB2)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LMAX2                ", FFM%LMAX2)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LAFILT2              ", FFM%LAFILT2)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_AFILT2               ", FFM%AFILT2)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_IAFILT2              ", FFM%IAFILT2)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LSPARSDES            ", FFM%LSPARSDES)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_NRANK_SPARSDES       ", FFM%NRANK_SPARSDES)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_RDES_SPARSDES        ", FFM%RDES_SPARSDES)

            CALL SECTION_HEADER("Kernel settings")
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_W1                   ", FFM%W1)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_NHYP                 ", FFM%NHYP1)

            CALL SECTION_HEADER("Bayesian error estimation")
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_ICRITERIA            ", FF%ICRITERIA)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_IUPDATE_CRITERIA     ", FF%IUPDATE_CRITERIA)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LCTIFOR_USE_FROM_FILE", FF%LCTIFOR_USE_FROM_FILE)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_CTIFOR               ", FF%CTIFOR * FUNIT)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_SCLC_CTIFOR          ", FF%SCALE_CTIFOR)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_CSIG                 ", FF%CSIG)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_CSLOPE               ", FF%CSLOPE)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_CX                   ", FF%CX)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_CDOUB                ", FF%CDOUB)

            CALL SECTION_HEADER("Sparsification and regression")
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_EPS_LOW              ", FF%ML_EPS_LOW)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_EPS_REG              ", FF%ML_EPS_REG)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_IALGO_LINREG         ", FF%IALGO_LINREG)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_ISVD                 ", FF%ISVD)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_IREG                 ", FFM%IREG)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_SIGV0                ", FFM%SIGV0)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_SIGW0                ", FFM%SIGW0)

            CALL SECTION_HEADER("Weighting and sampling")
            ! We need to take care of the different units if IWEIGHT=1
            IF (FF%IWEIGHT.EQ.1) THEN
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_WTOTEN               ", FF%WTOTEN * EUNIT)
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_WTIFOR               ", FF%WTIFOR * FUNIT)
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_WTSIF                ", FF%WTSIF * SUNIT)
            ELSE
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_WTOTEN               ", FF%WTOTEN)
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_WTIFOR               ", FF%WTIFOR)
               CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_WTSIF                ", FF%WTSIF)
            ENDIF
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_MHIS                 ", FF%MHIS)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_NMDINT               ", FF%NMDINT)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_IWEIGHT              ", FF%IWEIGHT)

            CALL SECTION_HEADER("Static array sizes")
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_MB                   ", FFM%MB)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_MCONF                ", FF%MCONF)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_MCONF_NEW            ", FF%MCONF_NEW)

            CALL SECTION_HEADER("Special features")
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LHEAT                ", FFM%LHEAT)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LCOUPLE              ", FFM%LCOUPLE)
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_NATOM_COUPLED        ", FFM%NATOM_COUPLED)
            WRITE(LOGF%IU, '(A," : ",I13,1X,A,1X,A)') &
               TAGLIST_DESCRIPTION("ML_ICOUPLE"),     &
               SIZE(FFM%ICOUPLE),                     &
               TAG_LIST%SYMBOL("ML_ICOUPLE"),         &
               "SIZE( ML_ICOUPLE )"
            BLOCK
            CHARACTER(LEN=TAGLIST_GET_DESCRIPTION_WIDTH()) :: TMP_STRING
            DO I = 1, SIZE(FFM%ICOUPLE)
               TMP_STRING = "Ion index"
               WRITE(LOGF%IU, '(A," : ",I13,1X,A,1X,A,I0,A)') &
                  ADJUSTR(TMP_STRING),                        &
                  FFM%ICOUPLE(I),                             &
                  TAG_LIST%SYMBOL("ML_ICOUPLE"),              &
                  "ML_ICOUPLE(",                              &
                  I,                                          &
                  ")" 
            ENDDO
            END BLOCK
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_RCOUPLE              ", FFM%RCOUPLE)

            CALL SECTION_HEADER("Reference energies")
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_ISCALE_TOTEN         ", FFM%ISCALE_TOTEN)
            WRITE(LOGF%IU, '(A," : ",I13,1X,A,1X,A)') &
               TAGLIST_DESCRIPTION("ML_EATOM_REF"),   &
               SIZE(FF%EATOM_REF),                    &
               TAG_LIST%SYMBOL("ML_EATOM_REF"),       &
               "SIZE( ML_EATOM_REF )"
            BLOCK
            CHARACTER(LEN=TAGLIST_GET_DESCRIPTION_WIDTH() - 3) :: TMP_STRING
            DO I = 1, SIZE(FF%EATOM_REF)
               TMP_STRING = "Element"
               WRITE(LOGF%IU, '(A," ",A," : ",ES13.5,1X,A,1X,A,I0,A)') &
                  ADJUSTR(TMP_STRING),                                 &
                  FF%TYPE(I),                                          &
                  FF%EATOM_REF(I) * EUNIT,                             &
                  TAG_LIST%SYMBOL("ML_EATOM_REF"),                     &
                  "ML_EATOM_REF(",                                     &
                  I,                                                   &
                  ")" 
            ENDDO
            END BLOCK

            CALL SECTION_HEADER("Output options")
            CALL PRINT_TAG(LOGF%IU, TAG_LIST, "ML_LEATOM               ", FFM%LEATOM)

            ! Write input tag deprecation warnings.
            BLOCK
            !CHARACTER(LEN=:), ALLOCATABLE :: DEPRECATED(:,:)
            CHARACTER(LEN=TAGLIST_GET_TAG_WIDTH()), &
               ALLOCATABLE :: DEPRECATED(:,:)
            CALL TAG_LIST%GET_DEPRECATED(DEPRECATED)
            IF (SIZE(DEPRECATED, 1) > 0) THEN
               WRITE(LOGF%IU, *)
               WRITE(LOGF%IU, '(A)') 'WARNING: Deprecated input tags were &
                  &found, please consider the replacements:'
               DO I = 1, SIZE(DEPRECATED, 1)
                  WRITE(LOGF%IU, '(("* "),A,(" => "),A)') &
                     ADJUSTL(TRIM(DEPRECATED(I, 1))),     &
                     ADJUSTL(TRIM(DEPRECATED(I, 2)))
               END DO
            END IF
            END BLOCK

            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(140("*"))')

            CONTAINS

            SUBROUTINE SECTION_HEADER(HEADER)

               CHARACTER(LEN=*), INTENT(IN) :: HEADER

               WRITE(LOGF%IU, *)
               WRITE(LOGF%IU, *)
               WRITE(LOGF%IU, '(A)') HEADER
               WRITE(LOGF%IU, '(140("-"))')

            END SUBROUTINE SECTION_HEADER

         END SUBROUTINE LOGSECTION_PARAMETERS

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file section about memory consumption.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGSECTION_MEMORY(LOGF, LTYPE, MEM)

            USE ML_FF_STRUCT, ONLY: LOG_PAR, MEMORY
            IMPLICIT NONE

            TYPE(LOG_PAR),    INTENT(IN) :: LOGF !< Log file setup.
            !> Whether this output contains estimated (0) or actual (1) memory
            !> consumption.
            INTEGER,          INTENT(IN) :: LTYPE
            !> Structure containing memory consumption information.
            TYPE(MEMORY),     INTENT(IN) :: MEM

            INTEGER                     :: I
            ! Width of memory string, match with MF string below!
            INTEGER, PARAMETER          :: W = 45
            ! Memory line output format.
            CHARACTER(LEN=*), PARAMETER :: MF = '(A45," : ",F8.1)'
            REAL(q)                     :: RMEM_DYN

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            ! Maximum dynamical memory.
            RMEM_DYN = MEM%RMEM(5)
            DO I = 6, 11
               RMEM_DYN = MAX(RMEM_DYN, MEM%RMEM(I))
            END DO

            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') "* MEMORY INFORMATION **********************&
               &**************************************************************&
               &***********************************"
            WRITE(LOGF%IU, *)

            IF (LTYPE == 0) THEN
               WRITE(LOGF%IU, '(A)') "Estimated memory consumption for ML &
                                     &force field generation (MB):"
            ELSE IF (LTYPE == 1) THEN
               WRITE(LOGF%IU, '(A)') "Actual memory consumption for ML &
                                     &force field generation (MB):"
            END IF
            WRITE(LOGF%IU, *)
            !WRITE(LOGF%IU, MF) LJ("Static  Coulomb variables"), MEM%RMEM(1)
            !WRITE(LOGF%IU, MF) LJ("Coulomb  descriptors"), MEM%RMEM(6)
            !WRITE(LOGF%IU, MF) LJ("Coulomb  regression"), MEM%RMEM(8)
            !WRITE(LOGF%IU, MF) LJ("Coulomb  prediction"), MEM%RMEM(10)
            WRITE(LOGF%IU, MF) LJ("Persistent allocations for force field"), &
                               MEM%RMEM(2) + MEM%RMEM_CMAT_ALL +             &
                               MEM%RMEM_PS_ALL + MEM%RMEM_SOAP_ALL +         &
                               MEM%RMEM_WMAT_ALL
            WRITE(LOGF%IU, '(A)') "|"
            WRITE(LOGF%IU, MF) LJ("|-- CMAT for basis"), &
                               MEM%RMEM_CMAT + MEM%RMEM_CMAT_ALL
            WRITE(LOGF%IU, MF) LJ("|-- FMAT for basis"), MEM%RMEM_FMAT
            WRITE(LOGF%IU, MF) LJ("|-- DESC for basis"), &
                               MEM%RMEM_PS + MEM%RMEM_PS_ALL
            WRITE(LOGF%IU, MF) LJ("|-- DESC product matrix"), &
                               MEM%RMEM_SOAP + MEM%RMEM_SOAP_ALL
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, MF) LJ("Persistent allocations for ab initio &
                                  &data"), &
                               MEM%RMEM(3) + MEM%RMEM(4)
            WRITE(LOGF%IU, '(A)') "|"
            WRITE(LOGF%IU, MF) LJ("|-- Ab initio data"), MEM%RMEM(3)
            WRITE(LOGF%IU, MF) LJ("|-- Ab initio data (new)"), MEM%RMEM(4)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, MF) LJ("Temporary allocations for &
                                  &sparsification"), &
                               MEM%RMEM(5)
            WRITE(LOGF%IU, '(A)') "|"
            WRITE(LOGF%IU, MF) LJ("|-- SVD matrices"), MAX(MEM%RMEM_SVD,MEM%RMEM_RST)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, MF) LJ("Other temporary allocations"), &
                               MEM%RMEM(7) + MEM%RMEM(9) + MEM%RMEM(11)
            WRITE(LOGF%IU, '(A)') "|"
            WRITE(LOGF%IU, MF) LJ("|-- Descriptors"), MEM%RMEM(7)
            WRITE(LOGF%IU, MF) LJ("|-- Regression"), MEM%RMEM(9)
            WRITE(LOGF%IU, MF) LJ("|-- Prediction"), MEM%RMEM(11)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, MF) LJ("Total memory consumption"), &
                               MEM%RMEM(2) + MEM%RMEM_CMAT_ALL +             &
                               MEM%RMEM_PS_ALL + MEM%RMEM_SOAP_ALL +         &
                               MEM%RMEM_WMAT_ALL +                           &
                               MEM%RMEM(3) + MEM%RMEM(4) +                   &
                               MEM%RMEM(5) +                                 &
                               MEM%RMEM(7) + MEM%RMEM(9) + MEM%RMEM(11) 
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(140("*"))')

            CONTAINS

            FUNCTION LJ(STRING)
               CHARACTER(LEN=*), INTENT(IN) :: STRING !< Input string.
               CHARACTER(LEN=W)             :: LJ !< Output string.
               LJ = STRLJ(STRING, W)
            END FUNCTION LJ

         END SUBROUTINE LOGSECTION_MEMORY

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Log file section about memory consumption.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         SUBROUTINE LOGSECTION_LOOP_HEADER(FF, LOGF)

            USE ML_FF_STRUCT, ONLY: LOG_PAR, FF_PAR
            IMPLICIT NONE

            TYPE(LOG_PAR), INTENT(IN) :: LOGF !< Log file setup.
            TYPE(FF_PAR),  INTENT(IN) :: FF !< FF setup.

            INTEGER :: I
            !> An automatically allocatable string for first column keyword.
            CHARACTER(LEN=:), ALLOCATABLE :: TMPK
            !> An automatically allocatable string for format specifiers.
            CHARACTER(LEN=:), ALLOCATABLE :: TMPF
            !> An automatically allocatable string for separator lines.
            CHARACTER(LEN=:), ALLOCATABLE :: TMPS
            !> A very long string to capture the column list.
            CHARACTER(LEN=1000) :: TMPL
            !> Another very long string to capture the column list.
            CHARACTER(LEN=1000) :: TMPL2

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') "* MAIN LOOP *******************************&
               &**************************************************************&
               &***********************************"
            WRITE(LOGF%IU, *)

            TMPK = "STATUS"
            TMPF = '(A,1X,A9,1X,A9,1X,A2,1X,A6,1X,A6,1X,A9,1X,A9)'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "state",                  &
               "is",                     &
               "doabin",                 &
               "iff",                    &
               "nsample",                &
               "ngenff"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line describes the overall status of each step.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ..... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'state ..... One-word description of step action'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '            - "accurate"  (1) : Errors are low, force field is used'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '            - "threshold" (2) : Errors exceeded threshold, structure is sampled from ab initio'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '            - "learning"  (3) : Stored configurations are used for training force field'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '            - "critical"  (4) : Errors are high, ab initio sampling and learning is enforced'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // '            - "predict"   (5) : Force field is used in prediction mode only, no error checking'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'is ........ Integer representation of above one-word description (integer in parenthesis)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'doabin .... Perform ab initio calculation (T/F)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'iff ....... Force field available (T/F, False after startup hints to possible convergence problems)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nsample ... Number of steps since last reference structure collection (sample = T)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'ngenff .... Number of steps since last force field generation (genff = T)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 8)
            WRITE(LOGF%IU, *)

            TMPK = "STDAB"
            TMPF = '(A,1X,A9,3(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "std_energy",             &
               "std_force",              &
               "std_stress"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line contains the standard deviation of the collected ab initio reference data.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ........ MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'std_energy ... Standard deviation in energy (eV atom^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'std_force .... Standard deviation in forces (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'std_stress ... Standard deviation in stress (kB)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 5)
            WRITE(LOGF%IU, *)

            TMPK = "ERR"
            TMPF = '(A,1X,A9,3(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "rmse_energy",            &
               "rmse_force",             &
               "rmse_stress"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line contains the RMSEs of the predictions with respect to ab initio results for the training data.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ......... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'rmse_energy ... RMSE of energies (eV atom^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'rmse_force .... RMSE of forces (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'rmse_stress ... RMSE of stress (kB)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 5)
            WRITE(LOGF%IU, *)

            TMPK = "LASTE"
            TMPF = '(A,1X,A9,3(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "last_energy",            &
               "last_force",             &
               "last_stress"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line contains the RMSEs of the predictions with respect to ab initio results for the last configuration in the training data.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ......... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'last_energy ... RMSE of energy for last configuration (eV atom^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'last_force .... RMSE of forces for last configuration (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'last_stress ... RMSE of stress for last configuration (kB)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 5)
            WRITE(LOGF%IU, *)

            TMPK = "CFE"
            TMPF = '(A,1X,A9,3(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "cur_for_mae",            &
               "cur_for_rmse",           &
               "cur_for_maxe"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line contains errors on forces with respect to ab initio results for the current structure.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep .......... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'cur_for_mae .... Mean absolute error of forces (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'cur_for_rmse ... Root mean square error of forces (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'cur_for_maxe ... Maximum error of forces (eV Angst^-1)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 5)
            WRITE(LOGF%IU, *)

            TMPK = "THRHIST"
            TMPF = '(A,1X,A9,1X,A4,1X,A16)'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "nhis",                   &
               "behis"                  
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line contains the last ML_MHIS entries of the Bayesian errors of the forces for training structures.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'Entries are listed from newest to oldest, i.e. the most recent addition is listed first.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nhis .... Number of the Bayesian error of the force in the history'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'behis ... Bayesian error of the force (eV Angst^-1)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 4)
            WRITE(LOGF%IU, *)

            TMPK = "THRUPD"
            TMPF = '(A,1X,A9,4(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "ctifor_prev",            &
               "ctifor_new",             &
               "std_sig",                &
               "slope_sig"                 
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line contains the new and old threshold for the maximum Bayesian error of forces.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ......... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'ctifor_prev ... Previous threshold for the maximum Bayesian error of forces (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'ctifor_new .... New threshold for the maximum Bayesian error of forces (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'std_sig ....... Standard deviation of the collected Bayesion errors of forces (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'slope_sig ..... Slope of the collected Bayesian errors of forces'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 6)
            WRITE(LOGF%IU, *)

            TMPK = "REGR"
            TMPF = '(A,1X,A9,1X,A4,1X,A4,4(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "linreg",                 &
               "ite",                    &
               "sigma_w_sqr",            &
               "sigma_v_sqr",            &
               "regratio",               &
               "error"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line contains results for each step of linear regression before the last step.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ......... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'linreg ........ Linear regression method (ML_IALGO_LINREG)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'ite ........... Number of iteration'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sigma_w_sqr ... Square of sigma_w (precision parameter) in evidence approximation'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sigma_v_sqr ... Square of sigma_v (noise parameter) in evidence approximation'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'regratio ...... Ratio between (sigma_v_sqr/sigma_w_sqr) and the highest eigenvalue (wmax)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'error.......... 0.5 times squared norm of vector containing combined errors of force field'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 8) 
            WRITE(LOGF%IU, *)



            TMPK = "REGRF"
            TMPF = '(A,1X,A9,1X,A4,1x,A4,6(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "linreg",                 &
               "ite",                    &
               "sigma_w_sqr",            &
               "sigma_v_sqr",            &
               "regratio",               &
               "error",                  &
               "wratio",                 &
               "wmax"                   
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line contains results of the last step of linear regression.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ......... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'linreg ........ Linear regression method (ML_IALGO_LINREG)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'ite ........... Number of iteration'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sigma_w_sqr ... Square of sigma_w (precision parameter) in evidence approximation'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'sigma_v_sqr ... Square of sigma_v (noise parameter) in evidence approximation'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'regratio ...... Ratio between (sigma_v_sqr/sigma_w_sqr) and the highest eigenvalue (wmax)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'error.......... 0.5 times squared norm of vector containing combined errors of force field'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'wratio ........ Ratio between the highest (wmax) and lowest (wmin) eigenvalue/singular value (wmax/wmin)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'wmax .......... Highest eigenvalue/singular value (wmax)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 10)
            WRITE(LOGF%IU, *)


            TMPK = "BEE"
            TMPF = '(A,1X,A9,6(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "bee_energy",             &
               "bee_max_force",          &
               "bee_ave_force",          &
               "threshold",              &
               "bee_max_stress",         &
               "bee_ave_stress"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows the Bayesian error estimations and the current threshold criterion,'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'immediately before re-training the force field.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ............ MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_energy ....... BEE of energy per atom (eV atom^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_max_force .... BEE of forces (max) (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_ave_force .... BEE of forces (average) (kB)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'threshold ........ Current value of threshold criterion (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_max_stress ... BEE of stresses (max) (kB)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_ave_stress ... BEE of stresses (average) (kB)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 8)
            WRITE(LOGF%IU, *)

            TMPK = "BEEF"
            TMPF = '(A,1X,A9,6(1X,A16))'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "bee_energy",             &
               "bee_max_force",          &
               "bee_ave_force",          &
               "threshold",              &
               "bee_max_stress",         &
               "bee_ave_stress"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows the Bayesian error estimations and the current threshold criterion.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ............ MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_energy ....... BEE of energy per atom (eV atom^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_max_force .... BEE of forces (max) (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_ave_force .... BEE of forces (average) (kB)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'threshold ........ Current value of threshold criterion (eV Angst^-1)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_max_stress ... BEE of stresses (max) (kB)'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'bee_ave_stress ... BEE of stresses (average) (kB)'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 8)
            WRITE(LOGF%IU, *)

            TMPK = "DCOUPLE"
            TMPF = '(A,1X,A9,1X,A16)'
            WRITE(TMPL, TMPF)            &
               STRLJ("# " // TMPK, C1W), &
               "nstep",                  &
               "der_H_lambda"
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows the derivative of the Hamiltonian with respect to coupling constant (dH/dlambda).'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep .......... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'der_H_lambda ... dH/dlambda'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 3)
            WRITE(LOGF%IU, *)

            TMPK = "LCONF"
            TMPF = '(A,1X,A9'
            DO I = 1, FF%MTYP
               TMPF = TMPF // ',1X,A2,2(1X,A9)'
            END DO
            TMPF = TMPF // ')'
            WRITE(TMPL, '(A,1X,A9)') STRLJ("# " // TMPK, C1W), "nstep"
            DO I = 1, FF%MTYP
               TMPL2 = TMPL
               WRITE(TMPL, '(A,1X,A2,2(1X,A9))') TRIM(TMPL2), &
               "el",                                          &
               "nlrc_old",                                    &
               "nlrc_new"
            END DO
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows the number of local configurations'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'which were sampled from ab initio reference calculations.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ...... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'el ......... Element symbol'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nlrc_old ... Previous number of local reference configurations for this element'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nlrc_new ... Current number of local reference configurations for this element'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 2 + 3 * FF%MTYP)
            WRITE(LOGF%IU, *)

            TMPK = "SPRSC"
            TMPF = '(A,1X,A9'
            DO I = 1, FF%MTYP
               TMPF = TMPF // ',1X,A2,4(1X,A9)'
            END DO
            TMPF = TMPF // ')'
            WRITE(TMPL, '(A,1X,A9)') STRLJ("# " // TMPK, C1W), "nstep"
            DO I = 1, FF%MTYP
               TMPL2 = TMPL
               WRITE(TMPL, '(A,2(1X,A9),1x,A2,2(1x,A9))') TRIM(TMPL2), &
               "nstr_prev",                                   &
               "nstr_spar",                                   &
               "el",                                          &
               "nlrc_prev",                                   &
               "nlrc_spar"
            END DO
            TMPS = CREATE_SEPARATOR(TMPK, TMPL)
            WRITE(LOGF%IU, '(A)') TMPS
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'This line shows the results of sparsification regarding the number'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'of reference structures and local reference configurations.'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " "
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstep ....... MD time step or input structure counter'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstr_prev ... Number of reference structures before sparsification'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nstr_spar ... Number of reference structures after sparsification'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'el .......... Element symbol'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nlrc_prev ... Number of local reference configurations before sparsification for this element'
            WRITE(LOGF%IU, '(A)') "# " // TMPK // " " // 'nlrc_spar ... Number of local reference configurations after sparsification for this element'
            CALL CREATE_HEADER(TMPK, TMPL, TMPS, TMPF, 2 + 5 * FF%MTYP)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(140("*"))')

            CONTAINS

            FUNCTION CREATE_SEPARATOR(KEYWORD, LINE) RESULT(SEPARATOR)

               CHARACTER(LEN=*), INTENT(IN) :: KEYWORD
               CHARACTER(LEN=*), INTENT(IN) :: LINE
               CHARACTER(LEN=:), ALLOCATABLE :: SEPARATOR

               SEPARATOR = "# " // KEYWORD // " " // &
                           REPEAT("#", LEN(TRIM(LINE)) - LEN(KEYWORD) - 3)

            END FUNCTION CREATE_SEPARATOR

            SUBROUTINE CREATE_HEADER(KEYWORD, LINE, SEPARATOR, FRMT, NCOL)

               CHARACTER(LEN=*), INTENT(IN) :: KEYWORD
               CHARACTER(LEN=*), INTENT(IN) :: LINE
               CHARACTER(LEN=*), INTENT(IN) :: SEPARATOR
               CHARACTER(LEN=*), INTENT(IN) :: FRMT
               INTEGER,          INTENT(IN) :: NCOL

               ! String representation of column numbers (up to 99, 2 chars).
               CHARACTER(LEN=2) :: COL_NUMBERS(1:99)

               DO I = 1, SIZE(COL_NUMBERS)
                  WRITE(COL_NUMBERS(I), '(I2)') I
               END DO

               WRITE(LOGF%IU, '(A)') SEPARATOR
               WRITE(LOGF%IU, '(A)') TRIM(LINE)
               IF (NCOL <= SIZE(COL_NUMBERS)) THEN
                  WRITE(LOGF%IU, FRMT)            &
                     STRLJ("# " // KEYWORD, C1W), &
                     COL_NUMBERS(2:NCOL)
               END IF
               WRITE(LOGF%IU, '(A)') SEPARATOR

            END SUBROUTINE CREATE_HEADER

         END SUBROUTINE LOGSECTION_LOOP_HEADER

         SUBROUTINE LOGSECTION_AB_INITIO_DATA(LOGF, AB, ABN)

            USE ML_FF_STRUCT, ONLY: LOG_PAR, ABINITIO, ABINITIONEW
            IMPLICIT NONE

            TYPE(LOG_PAR),     INTENT(IN) :: LOGF !< Log file setup.
            TYPE(ABINITIO),    INTENT(IN) :: AB
            TYPE(ABINITIONEW), INTENT(IN) :: ABN

            INTEGER :: I
            INTEGER :: J
            INTEGER :: TMP_NITYP

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') "* AVAILABLE AB INITIO DATA ****************&
               &**************************************************************&
               &***********************************"
            WRITE(LOGF%IU, *)


            WRITE(LOGF%IU, '(A,I9,A,I9,A)')                         &
               "Number of stored (maximum) ab initio structures: ", &
               AB%NCONF,                                            &
               " (",                                                &
               AB%MCONF - ABN%MCONF, ")"
            DO I = 1, AB%NSYS
               WRITE(LOGF%IU, '(A,I3,A,I9,A,A,A)') &
                  " * System ",                    &
                  I,                               &
                  " : ",                           &
                  AB%NCONF_SYS(I),                 &
                  ' , name: "',                    &
                  ADJUSTL(TRIM(AB%SZNAM2_ADD(I))), &
                  '"'
            END DO
            WRITE(LOGF%IU, '(A)') "Maximum number of atoms per element:"
            DO I = 1, AB%MTYP
               TMP_NITYP = 0
               DO J = 1, AB%NCONF
                  TMP_NITYP = MAX(AB%NITYP(I, J), TMP_NITYP)
               END DO
               WRITE(LOGF%IU, '(A,A2,A,I9)') &
                  " * Element ",             &
                  AB%TYPE(I),                &
                  " : ",                     &
                  TMP_NITYP
            END DO

            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(140("*"))')

         END SUBROUTINE LOGSECTION_AB_INITIO_DATA

         SUBROUTINE LOGSECTION_SPARS_DES(LOGF, FF, FFM, NNVAR2_OLD)

            USE ML_FF_STRUCT, ONLY: LOG_PAR, FF_PAR, FFM_PAR
            IMPLICIT NONE

            TYPE(LOG_PAR), INTENT(IN) :: LOGF !< Log file setup.
            TYPE(FF_PAR),  INTENT(IN) :: FF !< FF setup.
            TYPE(FFM_PAR), INTENT(IN) :: FFM !< FFM setup.
            INTEGER,       INTENT(IN) :: NNVAR2_OLD(:)

            INTEGER :: I

            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') "* INITIAL SPARSIFICATION OF DESCRIPTORS ***&
               &**************************************************************&
               &***********************************"
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(13X,A)') "Number of descriptors"
            WRITE(LOGF%IU, '(A9,1X,A9,1X,A9,1X,A)') &
               "Element",                           &
               "before",                            &
               "after",                             &
               "sparsification"
            DO I = 1, FF%MTYP
               IF (FFM%NB(I) > 0) THEN
                  WRITE(LOGF%IU, '(A9,1X,I9,1X,I9)') &
                     FF%TYPE(I),                     &
                     NNVAR2_OLD(I),                  &
                     FFM%NNVAR2(I)
               END IF
            END DO
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(140("*"))')

         END SUBROUTINE LOGSECTION_SPARS_DES

         SUBROUTINE LOGSECTION_TIME(LOGF, CPUTIM, VPUTIM)

            USE ML_FF_STRUCT, ONLY: LOG_PAR
            IMPLICIT NONE

            TYPE(LOG_PAR), INTENT(IN)    :: LOGF !< Log file setup.
            REAL(q),       INTENT(INOUT) :: CPUTIM(:) ! (1:12)
            REAL(q),       INTENT(INOUT) :: VPUTIM(:) ! (1:12)

            INTEGER :: I
            CHARACTER(LEN=:), ALLOCATABLE :: TFMT
            CHARACTER(LEN=50) :: STRTIM(12)

            STRTIM( 1) = "Setup (file I/O, parameters,...)"
            STRTIM( 2) = "Coulomb descriptor and design matrix"
            STRTIM( 3) = "pairwise descriptor and design matrix"
            STRTIM( 4) = "Descriptor and design matrix" ! formerly: "Many-body"
            STRTIM( 5) = "Sparsification of configurations"
            STRTIM( 6) = "Regression for Coulomb interactions"
            STRTIM( 7) = "Regression for pairwise interactions"
            STRTIM( 8) = "Regression" ! formerly: "for many-body interactions"
            STRTIM( 9) = "Prediction for Coulomb interactions"
            STRTIM(10) = "Prediction for pairwise interactions"
            STRTIM(11) = "Prediction" ! formerly: "for many-body interactions
            STRTIM(12) = "TOTAL"

            CPUTIM(12) = 0.0_q
            VPUTIM(12) = 0.0_q
            DO I = 1, 11
               CPUTIM(12) = CPUTIM(12) + CPUTIM(I)
               VPUTIM(12) = VPUTIM(12) + VPUTIM(I)
            ENDDO
            ! Only certains MPI ranks write output, others return early.
            IF (.NOT. ML_IO_WRITE) RETURN

            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A)') "* TIMING INFORMATION **********************&
               &**************************************************************&
               &***********************************"
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(A50,3X,A18,3X,A18)') &
               STRLJ("Program part", 50),         &
               "system clock (sec)",              &
               "cpu time (sec)"
            WRITE(LOGF%IU, '(A50,("-|-"),A18,("-|-"),A18)') &
               REPEAT("-", 50),                             &
               REPEAT("-", 18),                             &
               REPEAT("-", 18)
            TFMT = '(A50,(" | "),F18.3,(" | "),F18.3)'
            WRITE(LOGF%IU, TFMT) STRTIM( 1), VPUTIM( 1), CPUTIM( 1)
            ! WRITE(LOGF%IU, TFMT) STRTIM( 2), VPUTIM( 2), CPUTIM( 2)
            ! WRITE(LOGF%IU, TFMT) STRTIM( 3), VPUTIM( 3), CPUTIM( 3)
            WRITE(LOGF%IU, TFMT) STRTIM( 4), VPUTIM( 4), CPUTIM( 4)
            WRITE(LOGF%IU, TFMT) STRTIM( 5), VPUTIM( 5), CPUTIM( 5)
            ! WRITE(LOGF%IU, TFMT) STRTIM( 6), VPUTIM( 6), CPUTIM( 6)
            ! WRITE(LOGF%IU, TFMT) STRTIM( 7), VPUTIM( 7), CPUTIM( 7)
            WRITE(LOGF%IU, TFMT) STRTIM( 8), VPUTIM( 8), CPUTIM( 8)
            ! WRITE(LOGF%IU, TFMT) STRTIM( 9), VPUTIM( 9), CPUTIM( 9)
            ! WRITE(LOGF%IU, TFMT) STRTIM(10), VPUTIM(10), CPUTIM(10)
            WRITE(LOGF%IU, TFMT) STRTIM(11), VPUTIM(11), CPUTIM(11)
            WRITE(LOGF%IU, '(A50,("-|-"),A18,("-|-"),A18)') &
               REPEAT("-", 50),                             &
               REPEAT("-", 18),                             &
               REPEAT("-", 18)
            WRITE(LOGF%IU, TFMT) STRTIM(12), VPUTIM(12), CPUTIM(12)
            WRITE(LOGF%IU, *)
            WRITE(LOGF%IU, '(140("*"))')

         END SUBROUTINE LOGSECTION_TIME

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Output sparsification result
!****************************************************************************************************

        SUBROUTINE OUT_SPARS_DES_CUR (COMM_WORLD,LOGF,NNVAR2_NEW,NNVAR2_OLD,NRANK,TYPE,WRATIO_NEW,WRATIO_OLD)
          USE ML_FF_STRUCT, ONLY: ML_MPI_PAR, LOG_PAR
          IMPLICIT NONE
          TYPE (ML_MPI_PAR)               :: COMM_WORLD
          TYPE (LOG_PAR)               :: LOGF
          INTEGER         , INTENT(IN) :: NNVAR2_NEW
          INTEGER         , INTENT(IN) :: NNVAR2_OLD
          INTEGER         , INTENT(IN) :: NRANK
          CHARACTER(LEN=2), INTENT(IN) :: TYPE
          REAL(q)         , INTENT(IN) :: WRATIO_NEW
          REAL(q)         , INTENT(IN) :: WRATIO_OLD
          IF(COMM_WORLD%NODE_ME.EQ.0) THEN
             WRITE(LOGF%IU,1)
             WRITE(LOGF%IU,2) TYPE,NNVAR2_OLD,NNVAR2_NEW
          ENDIF
1         FORMAT (100("/")/5X,'Spars. des. results')
2         FORMAT (100("/")/5X,'Spars. results on',1X,A2/100("-")/ &
                          10X,'Descriptors from',1X,I5,1X,'to',1X,I5)
        END SUBROUTINE OUT_SPARS_DES_CUR
#endif

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !> Left-justify a given string with desired output width.
         !>
         !> If width > len(string) the output string is truncated (cut off at
         !> the right side).
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         PURE FUNCTION STRLJ(STRING, WIDTH)
            CHARACTER(LEN=*), INTENT(IN) :: STRING !< Input string.
            INTEGER, INTENT(IN)          :: WIDTH !< Desired width of output.
            CHARACTER(LEN=WIDTH)         :: STRLJ !< Output string.
            IF (WIDTH < LEN(STRING)) THEN
               STRLJ = STRING(1:WIDTH)
            ELSE
               STRLJ = ADJUSTL(STRING)
            END IF
         END FUNCTION STRLJ

      END MODULE LOGFILE
#endif
