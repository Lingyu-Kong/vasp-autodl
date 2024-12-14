#if defined(MPI) || defined(MPI_CHAIN)

!****************************************************************************************************
! This module contains all subroutines used to print some information to terminal.
!****************************************************************************************************

      MODULE PRINT2TERMINAL
        USE MPI_HELP
        IMPLICIT NONE

        CONTAINS

!****************************************************************************************************
! This subroutine checks for paramters that don't work together and aborts the code.
!****************************************************************************************************

        SUBROUTINE CHECK_PARAMTER_COMPATIBILITY(FF)
           USE ML_FF_STRUCT, ONLY : FF_PAR
           IMPLICIT NONE
           ! Input variables
           TYPE(FF_PAR),  INTENT(INOUT)     :: FF
           
           IF (.NOT.FF%LMLMB) THEN
              CALL ERROR_PRINT2TERMINAL ("ERROR: LMLMB must be .TRUE.!")
           ENDIF
           IF (FF%ISTART == 3 .AND. FF%NMDINT > 1) THEN
              CALL ERROR_PRINT2TERMINAL ("ERROR: ML_NMDINT must be 1 if &
                 &ML_ISTART=3! ")
           ENDIF
           IF (FF%ISTART.EQ.3.AND.FF%NSW.GT.1) THEN
              CALL ERROR_PRINT2TERMINAL("ERROR: ML_ISTART=3 &
                 &but NSW>1. Please set NSW to 1. Exiting...")
           ENDIF
           IF (FF%ISTART.EQ.3.AND.FF%IALGO_LINREG.GT.1) THEN
              CALL ERROR_PRINT2TERMINAL("ERROR: ML_ISTART=3 &
                 &and ML_IALGO_LINREG>1. Only ML_IALGO_LINREG=1 &
                 &allowed for ML_ISTART=3.")
           ENDIF
            
        END SUBROUTINE CHECK_PARAMTER_COMPATIBILITY

!****************************************************************************************************
! Subroutine to print error message on terminal.
!****************************************************************************************************

        SUBROUTINE ERROR_PRINT2TERMINAL (MESSAGE)
          USE ML_FF_STRUCT, ONLY : ML_MPI_PAR, ML_IO_WRITE
          IMPLICIT NONE
! Input variables
          CHARACTER(LEN=*), INTENT(IN) :: MESSAGE

          IF(ML_IO_WRITE) THEN
             WRITE(*,1)
             WRITE(*,*) MESSAGE
          ENDIF
! Stop the operations.
          CALL stopMpi(1)
          stop
! Format.
1         FORMAT(/' ----------------------------------------------------------------------------- '/, &
                  '|                                                                             |'/, &
                  '|     EEEEEEE  RRRRRR   RRRRRR   OOOOOOO  RRRRRR      ###     ###     ###     |'/, &
                  '|     E        R     R  R     R  O     O  R     R     ###     ###     ###     |'/, &
                  '|     E        R     R  R     R  O     O  R     R     ###     ###     ###     |'/, &
                  '|     EEEEE    RRRRRR   RRRRRR   O     O  RRRRRR       #       #       #      |'/, &
                  '|     E        R   R    R   R    O     O  R   R                               |'/, &
                  '|     E        R    R   R    R   O     O  R    R      ###     ###     ###     |'/, &
                  '|     EEEEEEE  R     R  R     R  OOOOOOO  R     R     ###     ###     ###     |'/, &
                  '|                                                                             |'/, &
                  ' ----------------------------------------------------------------------------- ')
        END SUBROUTINE ERROR_PRINT2TERMINAL

        ! Borrowed from tutor.F.
        subroutine stopMpi(errCode)
            integer, intent(in), optional :: errCode
            integer errCode_, request, ierr
            logical complete
            errCode_ = 0
            if (present(errCode)) errCode_ = errCode
            if (errCode_ /= 0) then
                call MPI_ibarrier(MPI_comm_world, request, ierr)
                complete = tryWaitingForOtherCPU(request, ierr)
                if (ierr /= MPI_success.or..not.complete) &
                    call MPI_abort(MPI_comm_world, 1, ierr)
            end if
            call MPI_finalize(ierr)
            if (ierr /= MPI_success) &
                call MPI_abort(MPI_comm_world, 1, ierr)
        end subroutine stopMpi

        ! Borrowed from tutor.F.
        logical function tryWaitingForOtherCPU(request, ierr) result(complete)
            integer, intent(in) :: request
            integer, intent(inout) :: ierr
            integer try
            if (ierr /= MPI_success) return
            do try = 1, 30
                call MPI_test(request, complete, MPI_STATUS_IGNORE, ierr)
                if (ierr /= MPI_success.or.complete) exit
                call sleep(1)
            end do
        end function tryWaitingForOtherCPU

      END MODULE

#endif
