module command_line

    use tutor, only: vtutor, isAlert

    implicit none

    integer, parameter :: default_status = -1

    type parsed_argument
        character(len=:), allocatable :: description
        logical :: should_stop = .false.
        integer :: status_ = default_status
    end type parsed_argument

contains

    subroutine parse_command_line(should_write)
        logical, intent(in) :: should_write
        character(len=256) text
        type(parsed_argument) argument
        integer index_
        do index_ = 1, command_argument_count()
            call get_command_argument(index_, text)
            argument = parse_argument(trim(text))
            if (should_write) call log_description(argument%description, argument%status_)
            if (argument%should_stop) call vtutor%stopCode()
        end do
    end subroutine parse_command_line

    pure type(parsed_argument) function parse_argument(text) result (argument)
        use version
        character(len=*), intent(in) :: text
        select case(text)
        case ('--version', '-v')
            argument%description = vasp()
            argument%should_stop = .true.
        case default
            argument%description = "Command line argument '" // trim(text) // "' was not understood."
            argument%status_ = isAlert
        end select
    end function parse_argument

    subroutine log_description(description, status_)
        use iso_fortran_env, only: output_unit
        character(len=*), intent(in) :: description
        integer, intent(in) :: status_
        select case (status_)
        case (default_status)
            write(output_unit, '(a)') description
        case (isAlert)
            call vtutor%alert(description)
        end select
    end subroutine log_description

end module command_line
