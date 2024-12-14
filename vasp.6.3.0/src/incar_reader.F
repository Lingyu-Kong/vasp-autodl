#include "symbol.inc"
module incar_reader
    !! Utilities to read input data from file.
    !!
    !! This module deals with input information of VASP. This input can be given
    !! in two different forms: Either as a text file where tags are matched to
    !! flags or as an HDF5 file containing the tags as elements.
    !!
    !! In case of the textual file, the file is read once line-by-line and each
    !! line is parsed to generate a map of tag -> corresponding assignment in the
    !! file. Then in a second step, when this tag is used in the code, it is
    !! processed via the default Fortran read to generate the integer, float, or
    !! whatever data is required. For the HDF5 input the first step is skipped
    !! and the data is just read as it is needed.

    use prec, only: q
    use string, only: indexer, any_newline

    implicit none

    !! Error issued when a requested tag is not found
    integer, parameter, public :: not_found_error = 3
    !! Error issued when a tag cannot be parsed to the desired format
    integer, parameter, public :: read_type_error = 5

    !! Stores the map from key used in the incar file and associated value
    !! as well as some data used for error handling and output
    type, public :: incar_tag
        character(len=:), allocatable :: key, val
        integer :: line_number = 1
        logical :: used = .false.
    end type incar_tag

    !! Contains the parsed information from the text input file.
    !!
    !! If the mode from HDF5 is used the tags are ignored and not read.
    !! The error variable contains a text with errors occuring during
    !! the reading.
    type, public :: incar_file
        type(incar_tag), allocatable :: tags(:)
        logical :: from_hdf5 = .false., to_hdf5 = .false., to_xml = .false.
        character(len=:), allocatable :: error
        integer :: image = 0
    end type incar_file

    !! Write a representation of the data in the incar_file type to a string.
    !!
    !! The string is generated in such a way that it can be subsequently read
    !! by the INCAR reader again. There is an optional argument to adjust
    !! the indent for every new line.
    interface str
        module procedure incar_to_string
    end interface str

    !! Read data from memory or the HDF5 file to a variable
    !!
    !! In case of reading from memory the string is converted using the Fortran
    !! read logic. If the tag is not found, the input value is kept. If you
    !! want a nested group separate the group name and the key by a slash. You
    !! can provide an optional ierr that will return error code
    !! 3 - if the tag is not found
    !! 5 - if converting the string to the desired type fails
    interface process_incar
        module procedure process_incar_character0, process_incar_integer0, &
            process_incar_integer1, process_incar_integer2, process_incar_integer3, &
            process_incar_logical0, process_incar_logical1, process_incar_logical2, &
            process_incar_logical3, process_incar_real0, process_incar_real1, &
            process_incar_real2, process_incar_real3
    end interface process_incar

    !
    ! Internal types and procedures
    !

    type, private :: key_comparison
        integer, allocatable :: new_groups(:), old_groups(:)
        integer :: first_different_group, first_character = 1
    end type key_comparison

    type, private :: line_content
        logical :: is_comment = .false., has_tag = .false., valid_group = .true.
        integer :: first = 1
        character(len=:), allocatable :: termination, group, error
        type(incar_tag) tag
        type(indexer) index_
    end type line_content

    type, private :: parse_status
        character(len=:), allocatable :: text
        integer, allocatable :: lines(:), line_numbers(:), tags(:)
        integer :: num_tags = 0, first = 1, last, indx = 0
        logical continue_
    end type parse_status

    integer, parameter :: indent_width = 4
    character, parameter :: separator = "/", assignment_ = "=", &
        delimiter_open = "{", delimiter_close = "}", &
        comment_bash = "#", comment_fortran = "!", &
        end_line = new_line('n'), end_inline = ";", &
        escape = "\", ignore_between = '"'
    character(len=*), parameter :: valid_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789" // separator
    character(len=*), parameter :: image_prefix = "IMAGE_"

    interface read_from_hdf5
        module procedure read_from_hdf5_character0, read_from_hdf5_integer0, &
            read_from_hdf5_integer1, read_from_hdf5_integer2, read_from_hdf5_integer3, &
            read_from_hdf5_logical0, read_from_hdf5_logical1, read_from_hdf5_logical2, &
            read_from_hdf5_logical3, read_from_hdf5_real0, read_from_hdf5_real1, &
            read_from_hdf5_real2, read_from_hdf5_real3
    end interface

    interface write_to_hdf5
        module procedure write_to_hdf5_character0, write_to_hdf5_integer0, &
            write_to_hdf5_integer1, write_to_hdf5_integer2, write_to_hdf5_integer3, &
            write_to_hdf5_logical0, write_to_hdf5_logical1, write_to_hdf5_logical2, &
            write_to_hdf5_logical3, write_to_hdf5_real0, write_to_hdf5_real1, &
            write_to_hdf5_real2, write_to_hdf5_real3
    end interface write_to_hdf5

    interface read_from_memory
        module procedure read_from_memory_character0, read_from_memory_integer0, &
            read_from_memory_integer1, read_from_memory_integer2, read_from_memory_integer3, &
            read_from_memory_logical0, read_from_memory_logical1, read_from_memory_logical2, &
            read_from_memory_logical3, read_from_memory_real0, read_from_memory_real1, &
            read_from_memory_real2, read_from_memory_real3
    end interface read_from_memory

    interface write_to_xml
        module procedure write_to_xml_character0, write_to_xml_integer0, &
            write_to_xml_integer1, write_to_xml_integer2, write_to_xml_integer3, &
            write_to_xml_logical0, write_to_xml_logical1, write_to_xml_logical2, &
            write_to_xml_logical3, write_to_xml_real0, write_to_xml_real1, &
            write_to_xml_real2, write_to_xml_real3
    end interface write_to_xml

    interface append_error
        module procedure append_error_content, append_error_incar
    end interface append_error

    private
    public str, incar_from_string, incar_from_unit, incar_from_file, process_incar, &
        count_elements, unused_incar_tags

contains

    !! Parse the contents of a string into memory of the form key -> val
    !!
    !! The string is read line by line, looking for special characters that
    !! indicate the assignment of a tag, a nesting in a group, or a multiline
    !! string. When these special characters are escaped by a \ they will not
    !! be interpreted.
    pure type(incar_file) function incar_from_string(text) result (incar)
        character(len=*), intent(in) :: text
        character(len=:), allocatable :: line
        type(parse_status) stat
        type(line_content) content
        integer line_number
        !
        incar%error = ""
        content = setup_content(group=separator, termination="")
        stat = setup_parser(content%index_, text)
        allocate(incar%tags(size(stat%tags)))
        do while (stat%continue_)
            call next_line(stat, line, line_number)
            content = parse_line(line, line_number, content)
            call update_status(stat, content)
            if (content%has_tag) incar%tags(stat%num_tags) = content%tag
            call append_error(incar, content%error)
        end do
        call reduce_to_actual_tags(incar, stat%num_tags)
        if (content%group /= separator) &
            call append_error(incar, "Group " // content%group // " never closed.")
    end function incar_from_string

    type(incar_file) function incar_from_file(filename) result(incar)
        use string, only: string_from_file
        use tutor, only: vtutor
        character(len=*), intent(in) :: filename
        character(len=:), allocatable :: file_content
        integer ierr
        file_content = string_from_file(filename, ierr)
        if (ierr /= 0) call vtutor%error("Error reading the file " // filename // &
            "please make sure the file exists and has the necessary read permissions.")
        incar = incar_from_string(file_content)
    end function incar_from_file

    type(incar_file) function incar_from_unit(iunit) result(incar)
        use string, only: string_from_unit
        use tutor, only: vtutor
        integer, intent(in) :: iunit
        character(len=:), allocatable :: file_content
        integer ierr
        file_content = string_from_unit(iunit, ierr)
        if (ierr /= 0) call vtutor%bug("Reading the incar file from the given unit" // &
            " failed. Check if it is opened for stream access.", __FILE__, __LINE__)
        incar = incar_from_string(file_content)
    end function incar_from_unit

    !! This routine counts the number of available elements within a certain tag.
    !! If the tag is not present 0 is returned.
    integer function count_elements(incar, key)
#ifdef VASP_HDF5
        use vhdf5_base, only: ih5ininputgroup_id, subgrp_incar, vh5_number_elements, vh5_error
#else
        use tutor, only: vtutor
#endif
        type(incar_file), intent(in) :: incar
        character(len=*), intent(in) :: key
        integer indx, ierr
#ifdef VASP_HDF5
        character(len=:), allocatable :: group_key
        group_key = subgrp_incar // "/" // key
#endif
        if (incar%from_hdf5) then
#ifdef VASP_HDF5
            ierr = vh5_number_elements(ih5ininputgroup_id, group_key, count_elements)
            if (ierr == 3) ierr = 0; VH5_CHECK(ierr)
#else
            call vtutor%bug("Flag from_hdf5 set, but code compiled without HDF5 support", __FILE__, __LINE__)
#endif
        else
            count_elements = 0
            indx = find_key(incar, key)
            if (indx == 0) return
            count_elements = count_elements_in_value(incar%tags(indx)%val)
        end if
    end function count_elements

    integer function count_elements_in_value(val, type_) result (num_elements)
        use string, only: uppercase
        character(len=*), intent(in) :: val
        character, intent(in), optional :: type_
        logical, parameter :: expand = .false.
        character type_loc
        character(len=32767) work
        integer, external :: nitems
        type_loc = "X" ! don't check type
        if (present(type_)) type_loc = type_
        num_elements = nitems(uppercase(val), work, expand, type_loc)
    end function count_elements_in_value

    !! Return a list of all tags in the INCAR file not used by any reader.
    pure function unused_incar_tags(incar) result (res)
        type(incar_file), intent(in) :: incar
        character(len=:), allocatable :: res
        type(line_content) content
        integer tag
        content%error = ""
        do tag = 1, size(incar%tags)
            content%tag = incar%tags(tag)
            if (unused_on_image(incar, content%tag%key)) cycle
            if (.not.content%tag%used) &
                call append_error(content, content%tag%key)
        end do
        res = content%error
    end function unused_incar_tags

    pure logical function unused_on_image(incar, key)
        use string, only: str
        type(incar_file), intent(in) :: incar
        character(len=*), intent(in) :: key
        character(len=:), allocatable :: current_image
        integer image_overwrite
        !
        current_image = image_prefix // str(incar%image) // separator
        if (index(key, image_prefix) == 0) then
            image_overwrite = find_key_inner(incar, current_image // key)
            unused_on_image = image_overwrite /= 0
        else
            unused_on_image = index(key, current_image) == 0
        end if
    end function unused_on_image

    function incar_to_string(incar, indentation) result (res)
        type(incar_file), intent(in) :: incar
        integer, intent(in), optional :: indentation
        character(len=:), allocatable :: res
        integer i, indent
        res = ""
        indent = 0
        if (present(indentation)) indent = indentation
        if (.not.allocated(incar%tags)) return
        if (size(incar%tags) == 0) return
        res = res // open_close_groups(indent, "", incar%tags(1)%key)
        do i = 1, size(incar%tags)
            if (i > 1) res = res // open_close_groups(indent, incar%tags(i-1)%key, incar%tags(i)%key)
            res = res // write_tag(indent, incar%tags(i))
        end do
        if (i > 1) res = res // open_close_groups(indent, incar%tags(i-1)%key, "")
    end function incar_to_string

    !
    ! The implementation of the process_incar routines are almost identical so
    ! we use preprocessor flags for rank and type to realize the different
    ! realizations. The implementation is in the incar_reader.inc file.
    !
#define INCLUDE_PROCESS_INCAR_BODY !--------------------------------------------

#define TYPE_CHARACTER
    subroutine process_incar_character0(incar, key, val, ierr)
#       define RANK 0
#       include "incar_reader.inc"
    end subroutine process_incar_character0
#undef TYPE_CHARACTER

#define TYPE_INTEGER
    subroutine process_incar_integer0(incar, key, val, ierr)
#       define RANK 0
#       include "incar_reader.inc"
    end subroutine process_incar_integer0
    !
    subroutine process_incar_integer1(incar, key, val, ierr)
#       define RANK 1
#       include "incar_reader.inc"
    end subroutine process_incar_integer1
    !
    subroutine process_incar_integer2(incar, key, val, ierr)
#       define RANK 2
#       include "incar_reader.inc"
    end subroutine process_incar_integer2
    !
    subroutine process_incar_integer3(incar, key, val, ierr)
#       define RANK 3
#       include "incar_reader.inc"
    end subroutine process_incar_integer3
#undef TYPE_INTEGER

#define TYPE_LOGICAL
    subroutine process_incar_logical0(incar, key, val, ierr)
#       define RANK 0
#       include "incar_reader.inc"
    end subroutine process_incar_logical0
    !
    subroutine process_incar_logical1(incar, key, val, ierr)
#       define RANK 1
#       include "incar_reader.inc"
    end subroutine process_incar_logical1
    !
    subroutine process_incar_logical2(incar, key, val, ierr)
#       define RANK 2
#       include "incar_reader.inc"
    end subroutine process_incar_logical2
    !
    subroutine process_incar_logical3(incar, key, val, ierr)
#       define RANK 3
#       include "incar_reader.inc"
    end subroutine process_incar_logical3
#undef TYPE_LOGICAL

#define TYPE_REAL
    subroutine process_incar_real0(incar, key, val, ierr)
#       define RANK 0
#       include "incar_reader.inc"
    end subroutine process_incar_real0
    !
    subroutine process_incar_real1(incar, key, val, ierr)
#       define RANK 1
#       include "incar_reader.inc"
    end subroutine process_incar_real1
    !
    subroutine process_incar_real2(incar, key, val, ierr)
#       define RANK 2
#       include "incar_reader.inc"
    end subroutine process_incar_real2
    !
    subroutine process_incar_real3(incar, key, val, ierr)
#       define RANK 3
#       include "incar_reader.inc"
    end subroutine process_incar_real3
#undef TYPE_REAL

#undef INCLUDE_PROCESS_INCAR_BODY
#define INCLUDE_READ_HDF5_BODY !------------------------------------------------

#define TYPE_CHARACTER
    subroutine read_from_hdf5_character0(group_key, val, ierr)
#       define RANK 0
#       include "incar_reader.inc"
    end subroutine read_from_hdf5_character0
#undef TYPE_CHARACTER

#define TYPE_INTEGER
    subroutine read_from_hdf5_integer0(group_key, val, ierr)
#       define RANK 0
#       include "incar_reader.inc"
    end subroutine read_from_hdf5_integer0
    !
    subroutine read_from_hdf5_integer1(group_key, val, ierr)
#       define RANK 1
#       include "incar_reader.inc"
    end subroutine read_from_hdf5_integer1
    !
    subroutine read_from_hdf5_integer2(group_key, val, ierr)
#       define RANK 2
#       include "incar_reader.inc"
    end subroutine read_from_hdf5_integer2
    !
    subroutine read_from_hdf5_integer3(group_key, val, ierr)
#       define RANK 3
#       include "incar_reader.inc"
    end subroutine read_from_hdf5_integer3
#undef TYPE_INTEGER

#define TYPE_LOGICAL
    subroutine read_from_hdf5_logical0(group_key, val, ierr)
#       define RANK 0
#       include "incar_reader.inc"
    end subroutine read_from_hdf5_logical0
    !
    subroutine read_from_hdf5_logical1(group_key, val, ierr)
#       define RANK 1
#       include "incar_reader.inc"
    end subroutine read_from_hdf5_logical1
    !
    subroutine read_from_hdf5_logical2(group_key, val, ierr)
#       define RANK 2
#       include "incar_reader.inc"
    end subroutine read_from_hdf5_logical2
    !
    subroutine read_from_hdf5_logical3(group_key, val, ierr)
#       define RANK 3
#       include "incar_reader.inc"
    end subroutine read_from_hdf5_logical3
#undef TYPE_LOGICAL

#define TYPE_REAL
    subroutine read_from_hdf5_real0(group_key, val, ierr)
#       define RANK 0
#       include "incar_reader.inc"
    end subroutine read_from_hdf5_real0
    !
    subroutine read_from_hdf5_real1(group_key, val, ierr)
#       define RANK 1
#       include "incar_reader.inc"
    end subroutine read_from_hdf5_real1
    !
    subroutine read_from_hdf5_real2(group_key, val, ierr)
#       define RANK 2
#       include "incar_reader.inc"
    end subroutine read_from_hdf5_real2
    !
    subroutine read_from_hdf5_real3(group_key, val, ierr)
#       define RANK 3
#       include "incar_reader.inc"
    end subroutine read_from_hdf5_real3
#undef TYPE_REAL

#undef INCLUDE_READ_HDF5_BODY
#define INCLUDE_WRITE_HDF5_BODY !-----------------------------------------------

#define TYPE_CHARACTER
    subroutine write_to_hdf5_character0(group_key, val)
#       define RANK 0
#       include "incar_reader.inc"
    end subroutine write_to_hdf5_character0
#undef TYPE_CHARACTER

#define TYPE_INTEGER
    subroutine write_to_hdf5_integer0(group_key, val)
#       define RANK 0
#       include "incar_reader.inc"
    end subroutine write_to_hdf5_integer0
    !
    subroutine write_to_hdf5_integer1(group_key, val)
#       define RANK 1
#       include "incar_reader.inc"
    end subroutine write_to_hdf5_integer1
    !
    subroutine write_to_hdf5_integer2(group_key, val)
#       define RANK 2
#       include "incar_reader.inc"
    end subroutine write_to_hdf5_integer2
    !
    subroutine write_to_hdf5_integer3(group_key, val)
#       define RANK 3
#       include "incar_reader.inc"
    end subroutine write_to_hdf5_integer3
#undef TYPE_INTEGER

#define TYPE_LOGICAL
    subroutine write_to_hdf5_logical0(group_key, val)
#       define RANK 0
#       include "incar_reader.inc"
    end subroutine write_to_hdf5_logical0
    !
    subroutine write_to_hdf5_logical1(group_key, val)
#       define RANK 1
#       include "incar_reader.inc"
    end subroutine write_to_hdf5_logical1
    !
    subroutine write_to_hdf5_logical2(group_key, val)
#       define RANK 2
#       include "incar_reader.inc"
    end subroutine write_to_hdf5_logical2
    !
    subroutine write_to_hdf5_logical3(group_key, val)
#       define RANK 3
#       include "incar_reader.inc"
    end subroutine write_to_hdf5_logical3
#undef TYPE_LOGICAL

#define TYPE_REAL
    subroutine write_to_hdf5_real0(group_key, val)
#       define RANK 0
#       include "incar_reader.inc"
    end subroutine write_to_hdf5_real0
    !
    subroutine write_to_hdf5_real1(group_key, val)
#       define RANK 1
#       include "incar_reader.inc"
    end subroutine write_to_hdf5_real1
    !
    subroutine write_to_hdf5_real2(group_key, val)
#       define RANK 2
#       include "incar_reader.inc"
    end subroutine write_to_hdf5_real2
    !
    subroutine write_to_hdf5_real3(group_key, val)
#       define RANK 3
#       include "incar_reader.inc"
    end subroutine write_to_hdf5_real3
#undef TYPE_REAL

#undef INCLUDE_WRITE_HDF5_BODY
#define INCLUDE_READ_INCAR_BODY !-----------------------------------------------

#define TYPE_CHARACTER
    logical function read_from_memory_character0(incar, key, val, ierr) result (first_usage)
#       define RANK 0
#       include "incar_reader.inc"
    end function read_from_memory_character0
#undef TYPE_CHARACTER

#define TYPE_INTEGER
    logical function read_from_memory_integer0(incar, key, val, ierr) result (first_usage)
#       define RANK 0
#       include "incar_reader.inc"
    end function read_from_memory_integer0
    !
    logical function read_from_memory_integer1(incar, key, val, ierr) result (first_usage)
#       define RANK 1
#       include "incar_reader.inc"
    end function read_from_memory_integer1
    !
    logical function read_from_memory_integer2(incar, key, val, ierr) result (first_usage)
#       define RANK 2
#       include "incar_reader.inc"
    end function read_from_memory_integer2
    !
    logical function read_from_memory_integer3(incar, key, val, ierr) result (first_usage)
#       define RANK 3
#       include "incar_reader.inc"
    end function read_from_memory_integer3
#undef TYPE_INTEGER

#define TYPE_LOGICAL
    logical function read_from_memory_logical0(incar, key, val, ierr) result (first_usage)
#       define RANK 0
#       include "incar_reader.inc"
    end function read_from_memory_logical0
    !
    logical function read_from_memory_logical1(incar, key, val, ierr) result (first_usage)
#       define RANK 1
#       include "incar_reader.inc"
    end function read_from_memory_logical1
    !
    logical function read_from_memory_logical2(incar, key, val, ierr) result (first_usage)
#       define RANK 2
#       include "incar_reader.inc"
    end function read_from_memory_logical2
    !
    logical function read_from_memory_logical3(incar, key, val, ierr) result (first_usage)
#       define RANK 3
#       include "incar_reader.inc"
    end function read_from_memory_logical3
#undef TYPE_LOGICAL

#define TYPE_REAL
    logical function read_from_memory_real0(incar, key, val, ierr) result (first_usage)
#       define RANK 0
#       include "incar_reader.inc"
    end function read_from_memory_real0
    !
    logical function read_from_memory_real1(incar, key, val, ierr) result (first_usage)
#       define RANK 1
#       include "incar_reader.inc"
    end function read_from_memory_real1
    !
    logical function read_from_memory_real2(incar, key, val, ierr) result (first_usage)
#       define RANK 2
#       include "incar_reader.inc"
    end function read_from_memory_real2
    !
    logical function read_from_memory_real3(incar, key, val, ierr) result (first_usage)
#       define RANK 3
#       include "incar_reader.inc"
    end function read_from_memory_real3
#undef TYPE_REAL

#undef INCLUDE_READ_INCAR_BODY
#define INCLUDE_WRITE_XML_BODY !------------------------------------------------

#define TYPE_CHARACTER
    subroutine write_to_xml_character0(key, val)
#       define RANK 0
#       include "incar_reader.inc"
    end subroutine write_to_xml_character0
#undef TYPE_CHARACTER

#define TYPE_INTEGER
    subroutine write_to_xml_integer0(key, val)
#       define RANK 0
#       include "incar_reader.inc"
    end subroutine write_to_xml_integer0
    !
    subroutine write_to_xml_integer1(key, val)
#       define RANK 1
#       include "incar_reader.inc"
    end subroutine write_to_xml_integer1
    !
    subroutine write_to_xml_integer2(key, val)
#       define RANK 2
#       include "incar_reader.inc"
    end subroutine write_to_xml_integer2
    !
    subroutine write_to_xml_integer3(key, val)
#       define RANK 3
#       include "incar_reader.inc"
    end subroutine write_to_xml_integer3
#undef TYPE_INTEGER

#define TYPE_LOGICAL
    subroutine write_to_xml_logical0(key, val)
#       define RANK 0
#       include "incar_reader.inc"
    end subroutine write_to_xml_logical0
    !
    subroutine write_to_xml_logical1(key, val)
#       define RANK 1
#       include "incar_reader.inc"
    end subroutine write_to_xml_logical1
    !
    subroutine write_to_xml_logical2(key, val)
#       define RANK 2
#       include "incar_reader.inc"
    end subroutine write_to_xml_logical2
    !
    subroutine write_to_xml_logical3(key, val)
#       define RANK 3
#       include "incar_reader.inc"
    end subroutine write_to_xml_logical3
#undef TYPE_LOGICAL

#define TYPE_REAL
    subroutine write_to_xml_real0(key, val)
#       define RANK 0
#       include "incar_reader.inc"
    end subroutine write_to_xml_real0
    !
    subroutine write_to_xml_real1(key, val)
#       define RANK 1
#       include "incar_reader.inc"
    end subroutine write_to_xml_real1
    !
    subroutine write_to_xml_real2(key, val)
#       define RANK 2
#       include "incar_reader.inc"
    end subroutine write_to_xml_real2
    !
    subroutine write_to_xml_real3(key, val)
#       define RANK 3
#       include "incar_reader.inc"
    end subroutine write_to_xml_real3
#undef TYPE_REAL

#undef INCLUDE_WRITE_XML_BODY

    !
    ! private procedures
    !

    subroutine incar_error_handling(key, ierr_, ierr)
        use string, only: str
        use tutor, only: vtutor
        character(len=*), intent(in) :: key
        integer, intent(in) :: ierr_
        integer, intent(out), optional :: ierr
        logical return_
        if (present(ierr)) ierr = ierr_
        return_ = ierr_ == 0 .or. ierr_ == not_found_error &
            .or. (ierr_ == read_type_error .and. present(ierr))
        if (return_) return
        call vtutor%error("Error reading item " // key // " from file INCAR.\n&
            &Error code was IERR= " // str(ierr_) // " ... .")
    end subroutine incar_error_handling

    pure subroutine append_error_incar(incar, error)
        type(incar_file), intent(inout) :: incar
        character(len=*), intent(in) :: error
        if (len(error) == 0) return
        if (len(incar%error) == 0) then
            incar%error = error
        else
            incar%error = incar%error // end_line // error
        end if
    end subroutine append_error_incar

    pure subroutine append_error_content(content, error)
        use string, only: str
        type(line_content), intent(inout) :: content
        character(len=*), intent(in) :: error
        character(len=:), allocatable :: line_error
        if (len(error) == 0) return
        line_error = "In line " // str(content%tag%line_number) // ": " // error
        if (len(content%error) == 0) then
            content%error = line_error
        else
            content%error = content%error // end_line // line_error
        end if
    end subroutine append_error_content

    pure integer function find_key(incar, key) result (indx)
        use string, only: str, uppercase
        type(incar_file), intent(in) :: incar
        character(len=*), intent(in) :: key
        character(len=:), allocatable :: key_common, key_image
        key_common = uppercase(key)
        key_image = image_prefix // str(incar%image) // separator // key_common
        indx = find_key_inner(incar, key_image)
        if (indx /= 0) return
        indx = find_key_inner(incar, key_common)
    end function find_key

    pure integer function find_key_inner(incar, key) result (indx)
        type(incar_file), intent(in) :: incar
        character(len=*), intent(in) :: key
        do indx = 1, size(incar%tags)
            if (key == incar%tags(indx)%key) return
        end do
        indx = 0
    end function find_key_inner

    pure type(line_content) function setup_content(group, termination) result (content)
        character(len=*), intent(in) :: group, termination
        content%group = group
        content%termination = termination
        content%error = ""
        content%index_ = indexer(escape, ignore_between)
    end function setup_content

    pure type(parse_status) function setup_parser(index_, text) result (stat)
        type(indexer), intent(in) :: index_
        character(len=*), intent(in) :: text
        stat%text = text
        call index_%find_in(text, [any_newline, end_inline], stat%lines)
        call index_%find_in(text, assignment_, stat%tags)
        call setup_line_numbers(text, stat%lines, stat%line_numbers)
        stat%continue_ = size(stat%tags) > 0
    end function setup_parser

    pure subroutine setup_line_numbers(text, lines, line_numbers)
        character(len=*), intent(in) :: text
        integer, intent(in) :: lines(:)
        integer, allocatable, intent(out) :: line_numbers(:)
        integer, allocatable :: line_breaks(:)
        type(indexer) index_
        integer line
        call index_%find_in(text, end_line, line_breaks)
        allocate(line_numbers(size(lines) + 1))
        line_numbers(1) = 1
        do line = 1, size(lines)
           line_numbers(line + 1) = count(line_breaks <= lines(line)) + 1
        end do
    end subroutine setup_line_numbers

    pure subroutine next_line(stat, line, line_number)
        type(parse_status), intent(inout) :: stat
        character(len=:), allocatable, intent(out) :: line
        integer, intent(out) :: line_number
        stat%indx = stat%indx + 1
        if (stat%indx > size(stat%lines)) then
            stat%last = len(stat%text)
        else
            stat%last = stat%lines(stat%indx)
        end if
        line = trim(adjustl(stat%text(stat%first:stat%last)))
        line_number = stat%line_numbers(stat%indx)
    end subroutine next_line

    pure subroutine update_status(stat, content)
        type(parse_status), intent(inout) :: stat
        type(line_content), intent(in) :: content
        stat%first = stat%last + len(content%termination)
        if (content%has_tag) stat%num_tags = stat%num_tags + 1
        stat%continue_ = stat%last < len(stat%text)
    end subroutine update_status

    pure subroutine reduce_to_actual_tags(incar, num_tags)
        type(incar_file), intent(inout) :: incar
        integer, intent(in) :: num_tags
        type(incar_tag), allocatable :: copy(:)
        if (num_tags < size(incar%tags)) then
            call move_alloc(incar%tags, copy)
            allocate(incar%tags(num_tags))
            incar%tags = copy(:num_tags)
        end if
    end subroutine reduce_to_actual_tags

    pure type(line_content) function parse_line(line, line_number, previous_line) result (res)
        character(len=*), intent(in) :: line
        integer, intent(in) :: line_number
        type(line_content), intent(in) :: previous_line
        character(len=:), allocatable :: trim_line
        !
        res = setup_content(previous_line%group, termination(previous_line%index_, line))
        res%tag%line_number = line_number
        res%is_comment = is_comment(line, previous_line)
        if (res%is_comment) return
        trim_line = trim_trailing_comment(res, line)
        call parse_open_group_if_present(res, trim_line)
        call parse_tag_if_present(res, trim_line)
        call parse_close_group_if_present(res, trim_line)
    end function parse_line

    pure function termination(index_, line)
        type(indexer), intent(in) :: index_
        character(len=*), intent(in) :: line
        character(len=:), allocatable :: termination
        integer, allocatable :: indices(:)
        call index_%find_in(line, [any_newline, end_inline], indices)
        if (size(indices) == 0) then
            termination = ""
        else
            termination = line(indices(1):)
        end if
    end function termination

    elemental logical function is_comment(line, previous_line)
        character(len=*), intent(in) :: line
        type(line_content), intent(in) :: previous_line
        is_comment = &
            (previous_line%is_comment.and.previous_line%termination == end_inline) &
            .or.(line(1:1) == comment_bash).or.(line(1:1) == comment_fortran)
    end function is_comment

    pure function trim_trailing_comment(content, line) result (res)
        type(line_content), intent(in) :: content
        character(len=*), intent(in) :: line
        character(len=:), allocatable :: res
        integer, allocatable :: indices(:)
        call content%index_%find_in(line, [comment_bash, comment_fortran], indices)
        if (size(indices) > 0) then
            res = line(:indices(1) - 1) // content%termination
        else
            res = line
        end if
    end function trim_trailing_comment

    pure subroutine parse_open_group_if_present(content, line)
        use string, only: uppercase
        type(line_content), intent(inout) :: content
        character(len=*), intent(in) :: line
        character(len=:), allocatable :: group
        integer, allocatable :: indices(:)
        integer indx
        !
        call content%index_%find_in(line, delimiter_open, indices)
        content%first = 1
        do indx = 1, size(indices)
            group = uppercase(trim(adjustl(line(content%first:indices(indx)-1))))
            content%first = indices(indx) + len(delimiter_open)
            content%group = content%group // group // separator
            if (len(group) == 0) &
                call append_error(content, "Missing group name before " // delimiter_open // ".")
            if (verify(group, valid_chars) /= 0) &
                call append_error(content, "Groupname " // content%group // " contains invalid characters.")
        end do
        content%valid_group = (index(content%group, repeat(separator, 2)) == 0) &
           .and. (verify(content%group, valid_chars) == 0)
    end subroutine parse_open_group_if_present

    pure subroutine parse_tag_if_present(content, line)
        type(line_content), intent(inout) :: content
        character(len=*), intent(in) :: line
        call split_at_assignment(content, line)
        if (.not.content%has_tag) return
        content%tag%key = format_key(content)
        content%tag%val = format_val(content%tag%val)
        call check_key_val_correct(content)
    end subroutine parse_tag_if_present

    pure subroutine split_at_assignment(content, line)
        type(line_content), intent(inout) :: content
        character(len=*), intent(in) :: line
        integer indx, last
        indx = index(line, assignment_)
        last = end_of_line_or_group(content, line)
        if (indx > last) &
            call append_error(content, "Found " // assignment_ // " after closing " // delimiter_close // ".")
        content%has_tag = 0 < indx .and. indx <= last .and. content%valid_group
        if (content%has_tag) then
            content%tag%key = trim(adjustl(line(content%first:indx-1)))
            content%tag%val = trim(adjustl(line(indx+1:last)))
        end if
    end subroutine split_at_assignment

    pure integer function end_of_line_or_group(content, line) result (last)
        type(line_content), intent(in) :: content
        character(len=*), intent(in) :: line
        integer, allocatable :: indices(:)
        call content%index_%find_in(line, delimiter_close, indices)
        if (size(indices) > 0) then
            last = indices(1) - 1
        else
            last = len(line) - len_trim(content%termination)
        end if
    end function end_of_line_or_group

    pure function format_key(content) result (key)
        use string, only: uppercase
        type(line_content), intent(in) :: content
        character(len=:), allocatable :: key
#if GNU_BEFORE_7_5
        key = content%group(2:len(content%group)) // uppercase(content%tag%key)
#else
        key = content%group(2:) // uppercase(content%tag%key)
#endif
    end function format_key

    pure function format_val(val_in) result (val)
        character(len=*), intent(in) :: val_in
        character(len=:), allocatable :: val, first_chars
        val = remove_escape_characters(val_in)
        first_chars = val(:min(len(ignore_between), len(val)))
        if (first_chars /= ignore_between) val = replace_newline_with_space(val)
        val = adjustl(trim(val))
    end function format_val

    pure function replace_newline_with_space(text) result (res)
        character(len=*), intent(in) :: text
        character(len=:), allocatable :: res
        integer, allocatable :: indices(:)
        type(indexer) index_
        integer indx, pos
        !
        call index_%find_in(text, any_newline, indices)
        res = ""
        pos = 1
        do indx = 1, size(indices)
            res = res // text(pos:indices(indx)-1) // " "
            pos = indices(indx) + 1
        end do
        res = res // text(pos:)
    end function replace_newline_with_space

    pure function remove_escape_characters(text) result (res)
        character(len=*), intent(in) :: text
        character(len=:), allocatable :: res
        integer, allocatable :: indices(:)
        integer indx, next_char, pos
        type(indexer) index_
        !
        call index_%find_in(text, escape, indices)
        res = ""
        pos = 1
        do indx = 1, size(indices)
            next_char = indices(indx) + 1
            if (next_char > len(text)) exit
            if (indx > 1 .and. next_char - pos == 1) cycle
            if (should_remove_escape_character(text(next_char:))) then
               res = res // text(pos:indices(indx) - 1)
               pos = next_char
            end if
        end do
        res = res // text(pos:)
    end function remove_escape_characters

    pure logical function should_remove_escape_character(text) result (res)
        ! skip if next character is special or if only whitespace in front of newline
        character(len=*), intent(in) :: text
        integer next_char
        next_char = next_escaped_character(text)
        if (next_char == index(text, end_line)) then
            res = next_char == verify(text, " ")
        else
            res = next_char == 1
        end if
    end function should_remove_escape_character

    pure integer function next_escaped_character(text) result (res)
        character(len=*), intent(in) :: text
        integer indices(9)
        logical mask(9)
        indices = index(text, [end_line, end_inline, delimiter_open, delimiter_close, &
            assignment_, escape, comment_bash, comment_fortran, ignore_between])
        mask = indices /= 0
        if (any(mask)) then
            res = minval(indices, mask)
        else
            res = 0
        end if
    end function next_escaped_character

    pure subroutine check_key_val_correct(content)
        type(line_content), intent(inout) :: content
        if (verify(content%tag%key, valid_chars) /= 0) then
            call append_error(content, "Tagname " // content%tag%key // " contains invalid characters.")
            content%has_tag = .false.
        end if
        if (len(content%tag%val) == 0) then
            call append_error(content, "Tag " // content%tag%key // " is not assigned.")
            content%has_tag = .false.
        end if
    end subroutine check_key_val_correct

    pure subroutine parse_close_group_if_present(content, line)
        use string, only: str
        type(line_content), intent(inout) :: content
        character(len=*), intent(in) :: line
        integer indx, pos, num_elements
        !
        call content%index_%count_in(line, delimiter_close, num_elements)
        do indx = 1, num_elements
            pos = index(content%group(:(len(content%group)-1)), separator, back=.true.)
            if (pos /= 0) then
                content%group = content%group(:pos)
            else
                call append_error(content, "Found unmatched closing " // delimiter_close // ".")
                content%group = separator
                exit
            end if
        end do
    end subroutine parse_close_group_if_present

    function open_close_groups(indent, old_key, new_key) result (res)
        integer, intent(inout) :: indent
        character(len=*), intent(in) :: old_key, new_key
        character(len=:), allocatable :: res
        type(key_comparison) stat
        !
        stat = compare_keys(old_key, new_key)
        res = close_groups(indent, stat)
        res = res // open_groups(indent, stat, new_key)
    end function open_close_groups

    pure type(key_comparison) function compare_keys(old_key, new_key) result (stat)
        character(len=*), intent(in) :: old_key, new_key
        integer group, first, last
        type(indexer) index_
        !
        call index_%find_in(new_key, separator, stat%new_groups)
        call index_%find_in(old_key, separator, stat%old_groups)
        first = 1
        do group = 1, min(size(stat%old_groups), size(stat%new_groups))
            last = stat%old_groups(group) - 1
            if (stat%new_groups(group) /= stat%old_groups(group)) exit
            if (new_key(first:last) /= old_key(first:last)) exit
        end do
        stat%first_different_group = group
        if (group > 1) stat%first_character = stat%new_groups(group - 1) + len(separator)
    end function compare_keys

    function close_groups(indent, stat) result (res)
        integer, intent(inout) :: indent
        type(key_comparison), intent(in) :: stat
        character(len=:), allocatable :: res
        integer group
        !
        res = ""
        do group = stat%first_different_group, size(stat%old_groups)
            indent = indent - indent_width
            res = res // write_indent(indent) // delimiter_close // end_line
        end do
    end function close_groups

    function open_groups(indent, stat, new_key) result (res)
        integer, intent(inout) :: indent
        type(key_comparison), intent(in) :: stat
        character(len=*), intent(in) :: new_key
        character(len=:), allocatable :: res
        integer group, first, last
        !
        res = ""
        first = stat%first_character
        do group = stat%first_different_group, size(stat%new_groups)
            last = stat%new_groups(group) - 1
            res = res // write_indent(indent) // new_key(first:last) // " " &
                // delimiter_open // end_line
            first = last + len(separator) + 1
            indent = indent + indent_width
        end do
    end function open_groups

    pure function write_tag(indent, tag) result (res)
        integer, intent(in) :: indent
        type(incar_tag), intent(in) :: tag
        character(len=:), allocatable :: res
        integer first
        first = index(tag%key, separator, back=.true.) + 1
#if GNU_BEFORE_7_5
        res = write_indent(indent) // tag%key(first:len(tag%key)) // " " // assignment_ &
            // " " // write_escaped(tag%val) // end_line
#else
        res = write_indent(indent) // tag%key(first:) // " " // assignment_ &
            // " " // write_escaped(tag%val) // end_line
#endif
    end function write_tag

    pure function write_indent(indent) result (res)
        integer, intent(in) :: indent
        character(len=:), allocatable :: res
        res = repeat(" ", indent)
    end function write_indent

    pure function write_escaped(text) result (res)
        character(len=*), intent(in) :: text
        character(len=:), allocatable :: res
        integer first, last
        res = ""
        first = 1
        last = next_escaped_character(text)
        do while (first <= last)
            res = res // text(first:last-1) // escape // text(last:last)
            first = last + 1
            last = last + next_escaped_character(text(first:))
        end do
        res = res // text(first:)
    end function write_escaped

end module incar_reader
