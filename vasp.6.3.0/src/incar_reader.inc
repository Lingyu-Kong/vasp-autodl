#if RANK == 0
#   define RANK_SUFFIX
#elif RANK == 1
#   define RANK_SUFFIX (:)
#elif RANK == 2
#   define RANK_SUFFIX (:,:)
#elif RANK == 3
#   define RANK_SUFFIX (:,:,:)
#endif

#if defined(TYPE_INTEGER)
#   define TYPE_VAL integer
#   define TYPE_CHECK "I"
#elif defined(TYPE_LOGICAL)
#   define TYPE_VAL logical
#   define TYPE_CHECK "L"
#elif defined(TYPE_REAL)
#   define TYPE_VAL real(q)
#   define TYPE_CHECK "F"
#elif defined(TYPE_CHARACTER)
#   define TYPE_VAL character(len=:), allocatable
#   define TYPE_CHECK "S"
#endif

#ifdef INCLUDE_PROCESS_INCAR_BODY
#   ifdef VASP_HDF5
        use vhdf5_base, only: subgrp_incar
#   endif
    type(incar_file), intent(inout) :: incar
    character(len=*), intent(in) :: key
    TYPE_VAL, intent(inout) :: val RANK_SUFFIX
    integer, intent(out), optional :: ierr
    logical first_usage
    integer ierr_
    character(len=:), allocatable :: group_key
#   ifdef VASP_HDF5
        group_key = subgrp_incar // "/" // key
#   else
        group_key = key
#   endif
    !
    if (incar%from_hdf5) then
        call read_from_hdf5(group_key, val, ierr_)
        first_usage = .true.
    else
        first_usage = read_from_memory(incar, key, val, ierr_)
    end if
    !
    call incar_error_handling(key, ierr_, ierr)
    if (ierr_ /= 0) return
    !
    if (incar%to_hdf5) call write_to_hdf5(group_key, val)
    if (incar%to_xml .and. first_usage) call write_to_xml(key, val)
#endif

#ifdef INCLUDE_READ_INCAR_BODY
    type(incar_file), intent(inout) :: incar
    character(len=*), intent(in) :: key
    TYPE_VAL, intent(inout) :: val RANK_SUFFIX
    integer, intent(out) :: ierr
    integer indx
    !
    indx = find_key(incar, key)
    if (indx == 0) then
        first_usage = .false.
        ierr = not_found_error
    else
#       ifdef TYPE_CHARACTER
            val = incar%tags(indx)%val; ierr = 0
#       else
            if (count_elements_in_value(incar%tags(indx)%val, TYPE_CHECK) > 0) then
                read(incar%tags(indx)%val, *, iostat=ierr) val
            else
                ierr = read_type_error
            end if
#       endif
        if (ierr == 0) then
            first_usage = .not.incar%tags(indx)%used
            incar%tags(indx)%used = .true.
        else
            first_usage = .false.
            ierr = read_type_error
        end if
    end if
#endif

#ifdef INCLUDE_READ_HDF5_BODY
#   ifdef VASP_HDF5
        use vhdf5_base, only: ih5ininputgroup_id, vh5_read, vh5_error
#   else
        use tutor, only: vtutor
#   endif
    character(len=*), intent(in) :: group_key
    TYPE_VAL, intent(inout) :: val RANK_SUFFIX
    integer, intent(out) :: ierr
    !
    ierr = 0
#   ifdef VASP_HDF5
        ierr = vh5_read(ih5ininputgroup_id, group_key, val)
        if (ierr /= not_found_error) then; VH5_CHECK(ierr); end if
#   else
        call vtutor%bug("Flag from_hdf5 set, but code compiled without HDF5 support", __FILE__, __LINE__)
#   endif
#endif

#ifdef INCLUDE_WRITE_HDF5_BODY
#   ifdef VASP_HDF5
        use vhdf5_base, only: ih5outinputgroup_id, vh5_write, vh5_error
#   else
        use tutor, only: vtutor
#   endif
    character(len=*), intent(in) :: group_key
    TYPE_VAL, intent(in) :: val RANK_SUFFIX
    integer ierr
#   ifdef VASP_HDF5
        ierr = vh5_write(ih5outinputgroup_id, group_key, val, skip_present=.true.); VH5_CHECK(ierr)
#   else
        call vtutor%bug("Flag to_hdf5 set, but code compiled without HDF5 support", __FILE__, __LINE__)
#   endif
#endif

#ifdef INCLUDE_WRITE_XML_BODY
    character(len=*), intent(in) :: key
    TYPE_VAL, intent(in) :: val RANK_SUFFIX
    integer idum
    real(q) rdum
    complex(q) cdum
    logical ldum
    character sdum
#   if RANK == 0
#       if defined(TYPE_INTEGER)
            call xml_incar(key, TYPE_CHECK, val, rdum, cdum, ldum, sdum, 1)
#       elif defined(TYPE_LOGICAL)
            call xml_incar(key, TYPE_CHECK, idum, rdum, cdum, val, sdum, 1)
#       elif defined(TYPE_REAL)
            call xml_incar(key, TYPE_CHECK, idum, val, cdum, ldum, sdum, 1)
#       elif defined(TYPE_CHARACTER)
            call xml_incar(key, TYPE_CHECK, idum, rdum, cdum, ldum, val, len(val))
#       endif
#   else
#       if defined(TYPE_INTEGER)
            call xml_incar_v(key, TYPE_CHECK, val, rdum, cdum, ldum, sdum, size(val))
#       elif defined(TYPE_LOGICAL)
            call xml_incar_v(key, TYPE_CHECK, idum, rdum, cdum, val, sdum, size(val))
#       elif defined(TYPE_REAL)
            call xml_incar_v(key, TYPE_CHECK, idum, val, cdum, ldum, sdum, size(val))
#       endif
#   endif
#endif

#undef TYPE_VAL
#undef TYPE_CHECK
#undef RANK_SUFFIX
#undef RANK
