    logical, intent(in) :: lopen ! ignored
    integer, intent(in) :: iu0, iu5 ! ignored
    character(len=*), intent(in)  :: flag_name
    logical, optional :: lwritexml ! ignored
    logical, optional :: lcontinue
    integer, optional :: foundnumber
    logical ignore_errors
    integer ierr, num_elements
#ifdef READ_ARRAY
    integer, intent(in) :: numb
    TYPE_VAL, pointer :: flag_value(:)
    flag_value(1:numb) => flag_value_
#   undef READ_ARRAY
#   undef TYPE_VAL
#endif
    !
    num_elements = count_elements(incar_f, flag_name)
    if (present(foundnumber)) foundnumber = num_elements
    if (num_elements == 0) then
        ierr = not_found_error
        return
    end if
    !
    ierr = 0
    ignore_errors = .false.
    if (present(lcontinue)) ignore_errors = lcontinue
    if (ignore_errors) then
        call process_incar(incar_f, flag_name, flag_value, ierr)
    else
        call process_incar(incar_f, flag_name, flag_value)
    end if
