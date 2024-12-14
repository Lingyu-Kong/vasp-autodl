module reader_tags
#include "symbol.inc"
  use prec
#ifdef VASP_HDF5
  use vhdf5_base, only: incar_found
#endif
  use incar_reader, only: incar_file, process_incar, count_elements, not_found_error

  implicit none

  type(incar_file) incar_f

  logical, save :: writexmlincar

  interface process_incar
     module procedure read_incar_logical_scalar
     module procedure read_incar_logical_array
     module procedure read_incar_integer_scalar
     module procedure read_incar_integer_array
     module procedure read_incar_double_scalar
     module procedure read_incar_double_array
     module procedure read_incar_double_array_2d
     module procedure read_incar_double_array_3d
     module procedure read_incar_string
  end interface process_incar

contains

  subroutine open_incar_if_found(iu5, lopen)
    integer :: iu5
    logical, optional :: lopen
! deprecated routine for historic compatibility
  end subroutine

  subroutine close_incar_if_found(iu5)
    integer :: iu5
! deprecated routine for historic compatibility
  end subroutine

  subroutine check_error_and_number(flag_name, iu0, ierr, n, numb, lcont)
    use string, only: str
    use tutor, only: vtutor
    character(len=*), intent(in)  :: flag_name
    integer, intent(in) :: iu0, n, numb
    integer, intent(inout) :: ierr
    logical, intent(inout) :: lcont

    ! add ierr=7 to denote n<numb
    if ((ierr==0).and.(n<numb)) ierr=7

    ! if an error was signaled and lcont=.TRUE. we return in error
    if ((ierr/=0).and.lcont) then
       lcont=.false. ; return
    endif

    ! stop if lcont=.FALSE. and the error is not "tag-not-present" (ierr/=3)
    if ((ierr/=0) .and.(ierr/=3)) then
       CALL vtutor%error("Error reading item " // flag_name // " from file INCAR. \n Error code &
          &was IERR= " // str(ierr) // " ... . Found N= " // str(n) // " data.")
    endif

    ! if the error is "tag-not-present" (ierr=3) we return in error
    if (ierr/=0) then
       lcont=.false. ; return
    endif

    ! if no error was signaled
    lcont=.true.

    return
  end subroutine check_error_and_number

  !
  ! logical scalar
  !
  subroutine read_incar_logical_scalar(lopen, iu0, iu5, flag_name, flag_value, ierr, lwritexml, lcontinue, foundnumber)
    logical, intent(inout) :: flag_value
#   include "reader_base.inc"
  end subroutine read_incar_logical_scalar
  !
  ! logical array 1d
  !
  subroutine read_incar_logical_array(lopen, iu0, iu5, flag_name, flag_value_, numb, ierr, lwritexml, lcontinue, foundnumber)
    logical, dimension(:), target, contiguous, intent(inout) :: flag_value_
#   define READ_ARRAY
#   define TYPE_VAL logical
#   include "reader_base.inc"
  end subroutine read_incar_logical_array
  !
  ! integer scalar
  !
  subroutine read_incar_integer_scalar(lopen, iu0, iu5, flag_name, flag_value, ierr, lwritexml, lcontinue, foundnumber)
    integer, intent(inout) :: flag_value
#   include "reader_base.inc"
  end subroutine read_incar_integer_scalar
  !
  ! integer array 1d
  !
  subroutine read_incar_integer_array(lopen, iu0, iu5, flag_name, flag_value_, numb, ierr, lwritexml, lcontinue, foundnumber)
    integer, dimension(:), target, contiguous, intent(inout) :: flag_value_
#   define READ_ARRAY
#   define TYPE_VAL integer
#   include "reader_base.inc"
  end subroutine read_incar_integer_array
  !
  ! real scalar
  !
  subroutine read_incar_double_scalar(lopen, iu0, iu5, flag_name, flag_value, ierr, lwritexml, lcontinue, foundnumber)
    real(q), intent(inout) :: flag_value
#   include "reader_base.inc"
  end subroutine read_incar_double_scalar
  !
  ! real array 1d
  !
  subroutine read_incar_double_array(lopen, iu0, iu5, flag_name, flag_value_, numb, ierr, lwritexml, lcontinue, foundnumber)
    real(q), dimension(:), target, contiguous, intent(inout) :: flag_value_
#   define READ_ARRAY
#   define TYPE_VAL real(q)
#   include "reader_base.inc"
  end subroutine read_incar_double_array

  subroutine read_incar_double_array_2d(lopen, iu0, iu5, flag_name, flag_value_, numb, ierr, lwritexml, lcontinue, foundnumber)
    real(q), dimension(:,:), target, contiguous, intent(inout) :: flag_value_
#   define READ_ARRAY
#   define TYPE_VAL real(q)
#   include "reader_base.inc"
  end subroutine read_incar_double_array_2d

  subroutine read_incar_double_array_3d(lopen, iu0, iu5, flag_name, flag_value_, numb, ierr, lwritexml, lcontinue, foundnumber)
    real(q), dimension(:,:,:), target, contiguous, intent(inout) :: flag_value_
#   define READ_ARRAY
#   define TYPE_VAL real(q)
#   include "reader_base.inc"
  end subroutine read_incar_double_array_3d

  subroutine read_incar_string(lopen, iu0, iu5, flag_name, flag_value_, numb, ierr, lwritexml, lcontinue, foundnumber)
    character(len=*), intent(inout) :: flag_value_
    character(len=:), allocatable :: flag_value
    integer, intent(in) :: numb
#   include "reader_base.inc"
    if (allocated(flag_value)) flag_value_ = flag_value
  end subroutine read_incar_string

end module reader_tags
