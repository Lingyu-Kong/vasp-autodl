#include "symbol.inc"

module electron_common

    use base, only: in_struct, info_struct, q
    use mpimy, only: communic

    implicit none

    private
    public checkAbort, testBreakCondition, testNumberOfStep

contains

    subroutine testBreakCondition(energyDiff, evalDiff, numStep, noMixing, info, abortWithoutConv)
        real(q), intent(in) :: energyDiff, evalDiff
        integer, intent(in) :: numStep
        logical, intent(in) :: noMixing
        type(info_struct), intent(inout) :: info
        logical, intent(out) :: abortWithoutConv
        ! eigenvalues and energy must be converged
        info%lAbort = max(abs(energyDiff), abs(evalDiff)) < info%eDiff
        ! charge-density not constant and in last cycle no change of charge
        if (.not.info%lMix .and. .not.info%lChCon .and. .not.noMixing) info%lAbort = .false.
        call testNumberOfStep(numStep, info%nElm, info, abortWithoutConv)
    end subroutine testBreakCondition

    subroutine testNumberOfStep(numStep, maxStep, info, abortWithoutConv)
        integer, intent(in) :: numStep, maxStep
        type(info_struct), intent(inout) :: info
        logical, intent(out) :: abortWithoutConv
        abortWithoutConv = .false.
        ! do not stop during the non-selfconsistent startup phase
        if (numStep <= abs(info%nElmDl)) info%lAbort = .false.
        ! do not stop before minimum number of iterations is reached
        if (numStep < min(abs(info%nElMin), maxStep)) info%lAbort = .false.
        ! but stop after INFO%NELM steps no matter where we are now
        if (numStep >= maxStep) then
             ! if abort was not set the code did not converged except if we
             ! performed a single shot calculation
             abortWithoutConv = .not.info%lAbort .and. (maxStep > 1)
             info%lAbort = .true.
        end if
    end subroutine testNumberOfStep

    subroutine checkAbort(abortWithoutConv, io, comm, info)
        logical, intent(in) :: abortWithoutConv
        type(in_struct), intent(in) :: io
        type(communic), intent(in) :: comm
        type(info_struct), intent(inout) :: info
        if (info%lAbort) then
            call writeAbortReason(abortWithoutConv, io, comm)
            return
        end if
        info%lSoft = readStopcar(io, comm)
        if (info%lSoft) call writeHardStop(io, comm)
    end subroutine checkAbort

    subroutine writeAbortReason(abortWithoutConv, io, comm)
        logical, intent(in) :: abortWithoutConv
        type(in_struct), intent(in) :: io
        type(communic), intent(in) :: comm
        integer node_me, ionode
#ifdef MPI
        node_me = comm%node_me
        ionode = comm%ionode
#endif
        io_begin
        write(io%iu6,'(a)') new_line('n')
        if (.not.abortWithoutConv) then
            write(io%iu6, '(a)') '------------------------ aborting loop because&
                & EDIFF is reached ----------------------------------------'
        else
            write(io%iu6, '(a)') '------------------------ aborting loop EDIFF &
                &was not reached (unconverged)  ----------------------------'
        end if
        write(io%iu6,'(a)') new_line('n')
        io_end
    end subroutine writeAbortReason

    logical function readStopcar(io, comm) result (lSoft)
        type(in_struct), intent(in) :: io
        type(communic), intent(in) :: comm
        integer idum, itmp, ierr, ncount
        real(q) rdum
        complex(q) cdum
        character charac
        lSoft = .false.
#ifndef noSTOPCAR
#ifndef F90_T3D
        call RDATAB(io%lOpen, 'STOPCAR', 99, 'LABORT', '=', '#', ';', 'L', &
            idum, rdum, cdum, lSoft, charac, ncount, 1, ierr)
        itmp = 0
        if (lSoft) itmp = 1
        CALLMPI(M_sum_i(comm, itmp, 1))
        lSoft = itmp > 0
#endif
#endif
    end function readStopcar

    subroutine writeHardStop(io, comm)
        type(in_struct), intent(in) :: io
        type(communic), intent(in) :: comm
        integer node_me, ionode
#ifdef MPI
        node_me = comm%node_me
        ionode = comm%ionode
#endif
        io_begin
        if (io%iu0 >= 0) write(io%iu0,*) 'hard stop encountered!  aborting job ...'
        write(io%iu6,'(a)') new_line('n')
        write(io%iu6,'(a)') '------------------------ aborting loop because hard &
            &stop was set ---------------------------------------'
        write(io%iu6,'(a)') new_line('n')
        io_end
    end subroutine writeHardStop

end module electron_common
