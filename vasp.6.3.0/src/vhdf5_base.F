#ifdef VASP_HDF5
#include "symbol.inc"
#define BAIL(e)  if (e/=0) return

module vhdf5_base
  use hdf5

  implicit none

  character(LEN=*), parameter :: GRP_INFOHEADER = "infoheader"
  character(LEN=*), parameter :: GRP_GRID = "grid"
  character(LEN=*), parameter :: GRP_EIGENVAL = "electron_eigenvalues"
  character(LEN=*), parameter :: GRP_PHONONS = "phonons"
  character(LEN=*), parameter :: GRP_DOS = "electron_dos"
  character(LEN=*), parameter :: GRP_POSITIONS = "positions"
  character(LEN=*), parameter :: GRP_HISTORY = "ion_dynamics"
  character(LEN=*), parameter :: GRP_CHARGE = "charge"
  character(LEN=*), parameter :: GRP_CONDUCTIVITY = "conductivity"
  character(LEN=*), parameter :: GRP_LINEAR_RESPONSE = "linear_response"
  character(LEN=*), parameter :: GRP_DIELECTRIC = "dielectric"
  character(LEN=*), parameter :: GRP_PROJECTORS = "projectors"
  character(LEN=*), parameter :: GRP_LCAO = "lcao"
  character(LEN=*), parameter :: GRP_RHFATOM = "rhfatom"
  character(LEN=*), parameter :: GRP_WAVE = "wave"
  character(LEN=*), parameter :: GRP_INTERMEDIATE = "intermediate"
  character(LEN=*), parameter :: GRP_RESULTS = "results"
  character(LEN=*), parameter :: GRP_INPUT = "input"
  character(LEN=*), parameter :: GRP_ORIGINAL = "original"
  character(LEN=*), parameter :: GRP_VERSION = "version"
  character(LEN=*), parameter :: SUBGRP_KPOINTS = "kpoints"
  character(LEN=*), parameter :: SUBGRP_INCAR = "incar"
  character(LEN=*), parameter :: SUBGRP_POSCAR = "poscar"
  character(LEN=*), parameter :: SUBGRP_POTCAR = "potcar"

  integer(HID_T):: ih5intermediategroup_id
  integer(HID_T):: ih5infileid
  integer(HID_T):: ih5outfileid
  integer(HID_T):: ih5wavefileid
  integer(HID_T):: ih5ininputgroup_id
  integer(HID_T):: ih5outinputgroup_id
  integer(HID_T):: ih5outoriginalgroup_id

  logical :: incar_found
  logical :: hdf5_found
  logical :: kpoints_found
  logical :: kpoints_hdf5_found
  logical :: poscar_found
  logical :: potcar_found       ! POTCAR found
  logical :: potcar_h5_found  ! POTCAR.h5 found
  logical :: vaspinh5_potcar_found  ! POTCAR dataset in vaspin.h5 found
  logical :: wavecar_found
  logical :: vaspwave_found

  logical :: lwriteh5 =.FALSE.

  integer(HID_T) :: create_inter_prop

  interface vh5_read
     module procedure vh5_read_logical_scalar
     module procedure vh5_read_logical_array_1d
     module procedure vh5_read_logical_array_2d
     module procedure vh5_read_logical_array_3d
     module procedure vh5_read_integer_scalar
     module procedure vh5_read_integer_array_1d
     module procedure vh5_read_integer_array_2d
     module procedure vh5_read_integer_array_3d
     module procedure vh5_read_real_scalar
     module procedure vh5_read_real_array_1d
     module procedure vh5_read_real_array_2d
     module procedure vh5_read_real_array_3d
     module procedure vh5_read_real_array_4d
     module procedure vh5_read_double_scalar
     module procedure vh5_read_double_array_1d
     module procedure vh5_read_double_array_2d
     module procedure vh5_read_double_array_3d
     module procedure vh5_read_double_array_4d
     module procedure vh5_read_double_array_5d
     module procedure vh5_read_double_array_6d
     module procedure vh5_read_double_complex_array_1d
     module procedure vh5_read_double_complex_array_2d
     module procedure vh5_read_double_complex_array_3d
     module procedure vh5_read_double_complex_array_5d
     module procedure vh5_read_string_scalar
     module procedure vh5_read_string_array_1d
  end interface vh5_read

  interface vh5_write
     module procedure vh5_write_logical_scalar
     module procedure vh5_write_logical_array_1d
     module procedure vh5_write_logical_array_2d
     module procedure vh5_write_logical_array_3d
     module procedure vh5_write_integer_scalar
     module procedure vh5_write_integer_array_1d
     module procedure vh5_write_integer_array_2d
     module procedure vh5_write_integer_array_3d
     module procedure vh5_write_real_scalar
     module procedure vh5_write_real_array_1d
     module procedure vh5_write_real_array_2d
     module procedure vh5_write_real_array_3d
     module procedure vh5_write_real_array_4d
     module procedure vh5_write_real_array_5d
     module procedure vh5_write_double_scalar
     module procedure vh5_write_double_array_1d
     module procedure vh5_write_double_array_2d
     module procedure vh5_write_double_array_3d
     module procedure vh5_write_double_array_4d
     module procedure vh5_write_double_array_5d
     module procedure vh5_write_double_array_6d
     module procedure vh5_write_double_complex_array_1d
     module procedure vh5_write_double_complex_array_2d
     module procedure vh5_write_double_complex_array_3d
     module procedure vh5_write_double_complex_array_4d
     module procedure vh5_write_double_complex_array_5d
     module procedure vh5_write_complex_array_1d
     module procedure vh5_write_complex_array_2d
     module procedure vh5_write_complex_array_3d
     module procedure vh5_write_complex_array_4d
     module procedure vh5_write_complex_array_5d
     module procedure vh5_write_string_scalar
     module procedure vh5_write_string_array_1d
  end interface vh5_write

  public :: vh5_dump_original_to_h5
  public :: read_potcar_from_potcarh5
  public :: read_potcar_from_vaspin

contains

  !======================================================================================
  ! non-hdf5 utility functions
  !======================================================================================
  elemental function logical_to_integer(logical_value) result(logasint)
    logical, intent(in) :: logical_value
    integer:: logasint
    !--------------------------------------------
    if (logical_value) then
       logasint = 1
    else
       logasint = 0
    end if
  end function logical_to_integer

  elemental function integer_to_logical(integer_value) result(intaslog)
    integer, intent(in) :: integer_value
    logical :: intaslog
    !--------------------------------------------
    if (integer_value == 1) then
       intaslog = .true.
    else
       intaslog = .false.
    end if
  end function integer_to_logical

  integer function check_rank(expected_rank, actual_rank) result(ierr)
    integer, intent(in) :: expected_rank, actual_rank
    if (expected_rank /= actual_rank) then
      ierr = actual_rank
    else
      ierr = 0
    end if
  end function check_rank

  !======================================================================================
  ! Initializing HDF5
  !======================================================================================

  !> @brief Initialize hdf5 library
  !>
  !> Needs to be called at the start of the application, before using any hdf5
  !> functionality
  !>
  !> @returns hdf5 lib error code
  function vh5_start() result(ierr)
    integer:: ierr
    !-------------------------
    call h5open_f(ierr); BAIL(ierr)
    call h5pcreate_f(H5P_LINK_CREATE_F, create_inter_prop, ierr); BAIL(ierr)
    call h5pset_create_inter_group_f(create_inter_prop, 1, ierr)
  end function vh5_start

  !> @brief close hdf5 library
  !>
  !> Needs to be called at the end of the application, to free global data structures
  !> allocated by vh5start
  !> @returns hdf5 lib error code
  function vh5_end() result(ierr)
    integer :: ierr
    !--------------------------
    call h5close_f(ierr)
  end function vh5_end

  !======================================================================================
  ! HDF5 file interface: open, close, create and delete h5 files
  !======================================================================================

  !> @brief open vasp hdf5 file for reading
  !>
  !> @param filename name of file to open
  !> @param fileid hdf5 library handle (HID_T) to file
  !> @returns hdf5 lib error code
  function vh5_file_open_read(filename, fileid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(out) :: fileid
    character(len=*), intent(in) :: filename
    !------------------------------------------------------
    ! open file
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, fileid, ierr)
  end function vh5_file_open_read

  !> @brief open vasp hdf5 file for read and write access
  !>
  !> @param filename name of file to open
  !> @param fileid hdf5 library handle (HID_T) to file
  !> @returns hdf5 lib error code
  function vh5_file_open_readwrite(filename, fileid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(out) :: fileid
    character(len=*), intent(in) :: filename
    !------------------------------------------------------
    ! open file
    ierr = 0; if (.not.lwriteh5) return
    call h5fopen_f(filename, H5F_ACC_RDWR_F, fileid, ierr)
  end function vh5_file_open_readwrite

  !> @brief close hdf5 file
  !>
  !> @param fileid hdf5 library handle (HID_T) to file
  !> @returns hdf5 lib error code
  function vh5_file_close(fileid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(out) :: fileid
    !------------------------------------------------------
!    call h5fflush_f(fileid, ierr)
    call h5fclose_f(fileid, ierr)
  end function vh5_file_close

  !> @brief create vasp hdf5 file
  !>
  !> create hdf5 file, fail if file already exists
  !>
  !> @param filename name of file to open
  !> @param fileid hdf5 library handle (HID_T) to file
  !> @returns hdf5 lib error code
  function vh5_file_create(filename, fileid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(out) :: fileid
    character(len=*), intent(in) :: filename
    !------------------------------------------------------
    ! create file
    ierr = 0; if (.not.lwriteh5) return
    CALL h5fcreate_f(filename, H5F_ACC_EXCL_F, fileid, ierr)
  end function vh5_file_create

  !> @brief create vasp hdf5 file
  !>
  !> create hdf5 file, overwrite if file already exists
  !> @param filename name of file to open
  !> @param fileid hdf5 library handle (HID_T) to file
  !> @returns hdf5 lib error code
  function vh5_file_create_or_overwrite(filename, fileid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(out) :: fileid
    character(len=*), intent(in) :: filename
    !------------------------------------------------------
    ! create file
    ierr = 0; if (.not.lwriteh5) return
    CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, fileid, ierr)
  end function vh5_file_create_or_overwrite


  function vh5_file_close_writing(fileid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(out) :: fileid
    !------------------------------------------------------
!    call h5fflush_f(fileid, ierr)
    ierr = 0; if (.not.lwriteh5) return
    call h5fclose_f(fileid, ierr)
  end function vh5_file_close_writing

  !> @brief delete file
  !>
  !> deletes file (any file, not just hdf5). works by opening the
  !> file with status='old', then closing it with status='delete'.
  !> I.e.: you need to have the proper access rights to open the file,
  !> and you need to pass a unit used for opening.
  !>
  !> @param filename name of file to delete
  !> @param unit io unit to use for deleting the file
  !> @returns error code. 0 if the file could be deleted, otherwise
  !>     the error status of the last open or close command
  function vh5_file_delete(filename, iounit) result(ierr)
    integer:: ierr
    integer:: iounit
    logical:: file_exists
    character(len=*), intent(in) :: filename
    !----------------------------------------
    integer:: istat
    !----------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    inquire(file=filename, exist=file_exists)
    if (file_exists) then
      open(unit=iounit, iostat=istat, file=filename, status='old')
      if (istat == 0) then
        close(iounit, iostat=ierr, status='delete')
      else
        ierr = istat
      endif
    else
      ierr = 0
    endif
  end function vh5_file_delete

  function vh5_file_flush(fileid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(in) :: fileid
    !------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    call h5fflush_f(fileid, H5F_SCOPE_LOCAL_F, ierr)
  end function vh5_file_flush

  subroutine vh5_error(ierr, srcfile, srcline)
    use string
    use tutor, only: vtutor
    integer,intent(in) :: ierr, srcline
    character(len=*), intent(in) :: srcfile
    if (ierr==0) return
    call vtutor%bug('HDF5 call in '//srcfile//':'//str(srcline)//' produced error: ' &
       // str(ierr), srcfile, srcline)
  end subroutine vh5_error

  !======================================================================================
  ! HDF5 group interface: create, open, close,... groups in H5 files
  !======================================================================================

  !> @brief open group, fail if it does not exist
  !>
  !> @close locid hdf5 handle (HID_T) to file or group
  !> @param groupname name of group to create
  !> @param groupid handle of the created group
  !> @returns hdf5 lib error code
  function vh5_group_open(locid, groupname, groupid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(in) ::   locid
    character(len=*), intent(in) :: groupname
    integer(HID_T), intent(out) ::  groupid
    !----------------------------------------
    call h5gopen_f(locid, groupname, groupid, ierr)
  end function vh5_group_open

  !> @brief close group, fail if it does not exist
  !>
  !> @param groupid hdf5 handle (HID_T) to group
  !> @returns hdf5 lib error code
  function vh5_group_close(groupid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(in) ::   groupid
    !----------------------------------------
    call h5gclose_f(groupid, ierr)
  end function vh5_group_close

  !> @brief create group in file, fail if group exists
  !>
  !> @param locid hdf5 handle (HID_T) to file or group
  !> @param groupname name of group to create
  !> @param groupid handle of the created group
  !> @returns hdf5 lib error code
  function vh5_group_create(locid, groupname, groupid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(in) ::   locid
    character(len=*), intent(in) :: groupname
    integer(HID_T), intent(out) ::  groupid
    !----------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Create a property list and set it to allow creation of intermediate groups
    call h5gcreate_f(locid, groupname, groupid, ierr, lcpl_id = create_inter_prop)
  end function vh5_group_create

  !> @brief create group in file, fail if group exists
  !>
  !> @param locid hdf5 handle (HID_T) to file or group
  !> @param groupname name of group to create
  !> @param groupid handle of the created group
  !> @returns hdf5 lib error code
  function vh5_group_open_or_create(locid, groupname, groupid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(in) ::   locid
    character(len=*), intent(in) :: groupname
    integer(HID_T), intent(out) ::  groupid
    !----------------------------------------
    logical :: link_exists
    !----------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    call h5lexists_f(locid, groupname, link_exists, ierr); BAIL(ierr)
    if (link_exists) then
       ! group already exists, open it
       ierr = vh5_group_open(locid, groupname, groupid); BAIL(ierr);
    else
       ! group does not exist, create
       ierr = vh5_group_create(locid, groupname, groupid); BAIL(ierr);
    end if
  end function vh5_group_open_or_create

  !> @brief close group if it was opened for writing, fail if it does not exist
  !>
  !> @param groupid hdf5 handle (HID_T) to group
  !> @returns hdf5 lib error code
  function vh5_group_close_writing(groupid) result(ierr)
    integer:: ierr
    integer(HID_T), intent(in) ::   groupid
    !----------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    call h5gclose_f(groupid, ierr)
  end function vh5_group_close_writing

  !======================================================================================
  ! HDF5 links interface
  !======================================================================================

  !> @brief check for the existence of a link
  !>
  !> @param locid hdf5 handle (HID_T) to group, dataset or attribute
  !> @param link_name name of the link
  !> @returns true if the link is found, false otherwise, in case of internal error stop the code
  recursive function vh5_exists(locid, link_name) result(exists)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: link_name
    logical :: exists
    integer :: ierr
    integer :: sep_loc
    ! Search for '/' from the right
    sep_loc = scan(link_name, '/', .true.)
    if (sep_loc <= 1) then
       ! Reached root/top-level name
       call h5lexists_f(locid, link_name, exists, ierr); VH5_CHECK(ierr)
    else
       ! Recursively call to get parent name
       exists = vh5_exists(locid, link_name(1 : sep_loc - 1))
       if (exists) then
          call h5lexists_f(locid, link_name, exists, ierr); VH5_CHECK(ierr)
       endif
    endif
  end function vh5_exists

  !======================================================================================
  ! HDF5 attributes interface
  !======================================================================================
  !> @brief Add a simple attribute to a dataset
  !>
  function vh5_add_attribute_dataset(groupid, dataset_name, attr_name, attr_val, skip_present) result(ierr)
    integer(HID_T), intent(in) :: groupid
    character(LEN=*), intent(in) :: dataset_name
    character(LEN=*), intent(in) :: attr_name
    character(LEN=*), intent(in) :: attr_val
    logical, optional, intent(in) :: skip_present
    !---------------------------------------------------------
    integer(HID_T) :: dset_id, attr_id, aspace, string_type
    integer(HSIZE_T) :: dimensions(1)
    integer(HSIZE_T) :: strlen
    logical :: my_skip_present, exists
    integer :: ierr
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    ! Create string datatype
    strlen = len(attr_val)
    if (strlen == 0) ierr = 1
    BAIL(ierr);
    call h5tcopy_f(H5T_FORTRAN_S1, string_type, ierr); BAIL(ierr)
    call h5tset_size_f(string_type, strlen, ierr); BAIL(ierr)
    ! Open dataspace
    call h5dopen_f(groupid, dataset_name, dset_id, ierr); BAIL(ierr)
    call h5screate_f(H5S_SCALAR_F, aspace, ierr); BAIL(ierr)
    call h5aexists_f(dset_id, attr_name, exists, ierr); BAIL(ierr)
    if (exists) then
       if (my_skip_present) return
       call h5aopen_f(dset_id, attr_name, attr_id, ierr); BAIL(ierr)
    else
       call h5acreate_f(dset_id, attr_name, string_type, aspace, attr_id, ierr); BAIL(ierr)
    endif
    ! Write the data
    call h5awrite_f(attr_id, string_type, attr_val, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5tclose_f(string_type, ierr); BAIL(ierr)
    call h5aclose_f(attr_id, ierr); BAIL(ierr)
    call h5sclose_f(aspace, ierr); BAIL(ierr)
    call h5dclose_f(dset_id, ierr)
  end function vh5_add_attribute_dataset

  !======================================================================================
  ! HDF5 dataset interface: read, write datasets
  !======================================================================================

  !> @brief if link to dataset exists open it
  !>
  !> @param dataset_name is the path to the dataset to open
  !> @param if the link is not found ierr=3
  !> @returns hdf5 handle (HID_T) to dataset
  function vh5_open_dataset_ifexists(locid, dataset_name, dset) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer(HID_T), intent(out) :: dset
    !----------------------------------------
    integer :: ierr
    !----------------------------------------
    if (vh5_exists(locid, dataset_name)) then
        call h5dopen_f(locid, dataset_name, dset, ierr)
    else
        ierr = 3
    endif
  end function vh5_open_dataset_ifexists

  !> @brief if link to dataset exists open it otherwise create it
  function vh5_create_open_dataset(locid, dataset_name, dset) result(ierr)
    integer(HID_T), intent(in) :: locid          !> hdf5 handle (HID_T) to root group
    character(LEN=*), intent(in) :: dataset_name !> string with path to the dataset to open
    integer(HID_T), intent(out) :: dset          !> hdf5 handle (HID_T) to dataset
    !------------------------------------------
    integer :: ierr
    integer(HID_T) :: dspace
    !------------------------------------------
    if (vh5_exists(locid, dataset_name)) then
       ! Open dataset
       call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    else
       ! Create dataspace
       call h5screate_f(H5S_SCALAR_F, dspace, ierr); BAIL(ierr)
       ! Create the dataset
       call h5dcreate_f(locid, dataset_name, H5T_IEEE_F32LE, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
       ! Close dataspace
       call h5sclose_f(dspace, ierr); BAIL(ierr)
    endif
  end function vh5_create_open_dataset

  !---------------------------------------------------------
  ! single precision
  !---------------------------------------------------------
  ! read
  !---------------------------------------------------------
  function vh5_read_real_scalar(locid, dataset_name, x) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(out) :: x
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dset
    integer(HSIZE_T) :: dimensions(1) ! dummy, not used
    !---------------------------------------------------------
    !
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    ! read the data
    call h5dread_f(dset, H5T_NATIVE_REAL, x, dimensions, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_read_real_scalar

  function vh5_read_real_subarray_nd(locid, dataset_name, rank, start, count, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: start(rank), count(rank)
    real(kind(0.)), intent(out) :: array(rank)
    !---------------------------------------------------------
    integer :: ierr
    integer(HID_T) :: dspace, memspace, dset
    !---------------------------------------------------------
    !
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    ! Get a slice of the array
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, start-1, count, ierr); BAIL(ierr)
    call h5screate_simple_f(rank, count, memspace, ierr); BAIL(ierr)
    ! Read the data from the dataset
    call h5dread_f(dset, H5T_NATIVE_REAL, array, count, ierr, memspace, dspace); BAIL(ierr)
    ! Close and release resources
    call h5sclose_f(memspace, ierr); BAIL(ierr)
    call h5sclose_f(dspace, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_read_real_subarray_nd

  function vh5_read_real_array_nd(locid, dataset_name, rank, dimensions, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(out) :: dimensions(*)
    real(kind(0.)), intent(out) :: array(*)
    !---------------------------------------------------------
    integer :: ierr, rank_or_error
    integer(HID_T) :: dspace, dset
    integer(HSIZE_T) :: maxdims(rank)
    !---------------------------------------------------------
    !
    dimensions(:rank) = 0
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    call h5sget_simple_extent_dims_f(dspace, dimensions, maxdims, rank_or_error)
    ierr = check_rank(rank, rank_or_error); BAIL(ierr)
    ! Read the data from the dataset
    call h5dread_f(dset, H5T_NATIVE_REAL, array, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5sclose_f(dspace, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_read_real_array_nd

  function vh5_read_real_array_1d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(out) :: array(:)
    integer, optional, intent(out) :: dimensions
    !---------------------------------------------------------
    integer :: ierr
    integer(HSIZE_T) :: my_dimensions(1)
    !---------------------------------------------------------
    ierr = vh5_read_real_array_nd(locid, dataset_name, 1, my_dimensions, array)
    dimensions = my_dimensions(1)
  end function vh5_read_real_array_1d

  function vh5_read_real_array_2d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(out) :: array(:,:)
    integer, optional, intent(out) :: dimensions(2)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(2)
    !---------------------------------------------------------
    ierr = vh5_read_real_array_nd(locid, dataset_name, 2, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_real_array_2d

  function vh5_read_real_array_3d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(out) :: array(:,:,:)
    integer, optional, intent(out) :: dimensions(3)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(3)
    !---------------------------------------------------------
    ierr = vh5_read_real_array_nd(locid, dataset_name, 3, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_real_array_3d

  function vh5_read_real_array_4d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(out) :: array(:,:,:,:)
    integer, optional, intent(out) :: dimensions(4)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(4)
    !---------------------------------------------------------
    ierr = vh5_read_real_array_nd(locid, dataset_name, 4, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_real_array_4d
  !
  ! write
  !
  function vh5_write_real_scalar(locid, dataset_name, x, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(in) :: x
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset
    integer(HSIZE_T) :: dimensions(1) ! dummy, not used
    logical :: my_skip_present
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    if (vh5_exists(locid,dataset_name)) then
      if(my_skip_present) return
      ! Open dataset
      call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    else
      ! Create dataspace and dataset
      call h5screate_f(H5S_SCALAR_F, dspace, ierr); BAIL(ierr)
      call h5dcreate_f(locid, dataset_name, H5T_IEEE_F32LE, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
      call h5sclose_f(dspace, ierr); BAIL(ierr)
    endif
    ! Write the data to the dataset
    call h5dwrite_f(dset, H5T_NATIVE_REAL, x, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5dclose_f(dset, ierr); BAIL(ierr)
    call h5sclose_f(dspace, ierr); BAIL(ierr)
  end function vh5_write_real_scalar

  function vh5_create_real_array_nd(locid, dataset_name, rank, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: dimensions(rank)
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    ! Create dataspace
    call h5screate_simple_f(rank, dimensions, dspace, ierr); BAIL(ierr)
    ! Create the dataset
    call h5dcreate_f(locid, dataset_name, H5T_IEEE_F32LE, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
    ! Close and release resources
    call h5dclose_f(dset, ierr); BAIL(ierr)
    call h5sclose_f(dspace, ierr); BAIL(ierr)
  end function vh5_create_real_array_nd

  function vh5_write_real_subarray_nd(locid, dataset_name, rank, start, count, array, extend) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: start(rank), count(rank)
    real(kind(0.)), intent(in) :: array(*)
    integer(HSIZE_T), optional, intent(in) :: extend(rank)
    integer :: ierr
    !---------------------------------------------------------
    integer :: i, rank_or_error
    integer(HID_T) :: dspace, memspace, dset
    integer(HSIZE_T) :: dimensions(rank), maxdims(rank), new_dimensions(rank)
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    !
    if (present(extend)) then
        ! Get dimensions
        call h5sget_simple_extent_dims_f(dspace, dimensions, maxdims, rank_or_error)
        ierr = check_rank(rank, rank_or_error); BAIL(ierr)
        new_dimensions = dimensions
        do i=1,rank
            if ((extend(i)>0).and.(start(i)>dimensions(i))) new_dimensions(i) = new_dimensions(i) + max(extend(i),start(i)-dimensions(i))
        end do
        if (any(new_dimensions /= dimensions)) then
            call h5dset_extent_f(dset, new_dimensions, ierr); BAIL(ierr)
            call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
        endif
    endif
    !
    ! Get a slice of the array
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, start-1, count, ierr); BAIL(ierr)
    call h5screate_simple_f(rank, count, memspace, ierr); BAIL(ierr)
    ! Write the data to the dataset
    call h5dwrite_f(dset, H5T_NATIVE_REAL, array, count, ierr, memspace, dspace); BAIL(ierr)
    !
    ! Close and release resources
    call h5sclose_f(memspace, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
    call h5sclose_f(dspace, ierr); BAIL(ierr)
  end function vh5_write_real_subarray_nd

  function vh5_write_real_array_nd(locid, dataset_name, rank, dimensions, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: dimensions(*)
    real(kind(0.)), intent(in) :: array(*)
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    ! Create dataspace
    call h5screate_simple_f(rank, dimensions, dspace, ierr); BAIL(ierr)
    ! Create the dataset
    call h5dcreate_f(locid, dataset_name, H5T_IEEE_F32LE, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
    ! Write the data to the dataset
    call h5dwrite_f(dset, H5T_NATIVE_REAL, array, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5dclose_f(dset , ierr); BAIL(ierr)
    call h5sclose_f(dspace, ierr); BAIL(ierr)
  end function vh5_write_real_array_nd

  function vh5_write_real_array_1d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(in) :: array(:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(1)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    if (vh5_exists(locid,dataset_name).and.my_skip_present) return
    dimensions(1) = size(array)
    ierr = vh5_write_real_array_nd(locid, dataset_name, 1, dimensions, array)
  end function vh5_write_real_array_1d

  function vh5_write_real_array_2d(locid, dataset_name, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(in) :: array(:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(2)
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    dimensions = shape(array)
    ierr = vh5_write_real_array_nd(locid, dataset_name, 2, dimensions, array)
  end function vh5_write_real_array_2d

  function vh5_write_real_array_3d(locid, dataset_name, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(in) :: array(:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(3)
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    dimensions = shape(array)
    ierr = vh5_write_real_array_nd(locid, dataset_name, 3, dimensions, array)
  end function vh5_write_real_array_3d

  function vh5_write_real_array_4d(locid, dataset_name, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(in) :: array(:,:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(4)
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    dimensions = shape(array)
    ierr = vh5_write_real_array_nd(locid, dataset_name, 4, dimensions, array)
  end function vh5_write_real_array_4d

  function vh5_write_real_array_5d(locid, dataset_name, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(in) :: array(:,:,:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(5)
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    dimensions = shape(array)
    ierr = vh5_write_real_array_nd(locid, dataset_name, 5, dimensions, array)
  end function vh5_write_real_array_5d

  function vh5_write_real_array_6d(locid, dataset_name, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.)), intent(in) :: array(:,:,:,:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(6)
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    dimensions = shape(array)
    ierr = vh5_write_real_array_nd(locid, dataset_name, 6, dimensions, array)
  end function vh5_write_real_array_6d

  !---------------------------------------------------------
  ! double precision
  !---------------------------------------------------------
  ! write
  !---------------------------------------------------------
  function vh5_write_double_scalar(locid, dataset_name, x, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(in) :: x
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset
    integer(HSIZE_T) :: dimensions(1) ! dummy, not used
    logical :: my_skip_present
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    if (vh5_exists(locid,dataset_name)) then
      if (my_skip_present) return
      ! Open dataset
      call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    else
      ! Create dataspace and dataset
      call h5screate_f(H5S_SCALAR_F, dspace, ierr); BAIL(ierr)
      call h5dcreate_f(locid, dataset_name, H5T_IEEE_F64LE, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
      call h5sclose_f(dspace, ierr); BAIL(ierr)
    endif
    ! Write the data to the dataset
    call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, x, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5dclose_f(dset , ierr); BAIL(ierr)
  end function vh5_write_double_scalar

  function vh5_create_double_array_nd(locid, dataset_name, rank, dimensions, max_dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: dimensions(rank)
    integer(HSIZE_T), optional, intent(in):: max_dimensions(rank)
    integer :: ierr
    !---------------------------------------------------------
    integer :: i
    integer(HSIZE_T) :: chunk_size(rank)
    integer(HID_T) :: dspace, dset, crp_list
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    ! Create dataspace
    if (present(max_dimensions)) then
        ! Modify dataset creation properties, i.e. enable chunking
        call h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, ierr); BAIL(ierr)
        chunk_size = 1
        do i=1,rank
            if (max_dimensions(i) /= H5S_UNLIMITED_F) cycle
            chunk_size(i) = 100
        end do
        call h5pset_chunk_f(crp_list, rank, chunk_size, ierr); BAIL(ierr)
        !create extensible dataset
        call h5screate_simple_f(rank, dimensions, dspace, ierr, max_dimensions); BAIL(ierr)
    else
        call h5screate_simple_f(rank, dimensions, dspace, ierr); BAIL(ierr)
    endif
    ! Create the dataset
    if (present(max_dimensions)) then
        call h5dcreate_f(locid, dataset_name, H5T_IEEE_F64LE, dspace, dset, ierr, crp_list, lcpl_id = create_inter_prop); BAIL(ierr)
        call h5pclose_f(crp_list, ierr); BAIL(ierr)
    else
        call h5dcreate_f(locid, dataset_name, H5T_IEEE_F64LE, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
    endif
    ! Close and release resources
    call h5dclose_f(dset , ierr); BAIL(ierr)
    call h5sclose_f(dspace, ierr); BAIL(ierr)
  end function vh5_create_double_array_nd

  function vh5_write_double_subarray_nd(locid, dataset_name, rank, start, count, array, extend) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: start(rank), count(rank)
    real(kind(0.d0)), intent(in) :: array(*)
    integer(HSIZE_T), optional, intent(in) :: extend(rank)
    integer :: ierr
    !---------------------------------------------------------
    integer :: i, rank_or_error
    integer(HID_T) :: dspace, memspace, dset
    integer(HSIZE_T) :: dimensions(rank), maxdims(rank), new_dimensions(rank)
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    !
    if (present(extend)) then
        ! Get dimensions
        call h5sget_simple_extent_dims_f(dspace, dimensions, maxdims, rank_or_error)
        ierr = check_rank(rank, rank_or_error); BAIL(ierr)
        new_dimensions = dimensions
        do i=1,rank
            if ((extend(i)>0).and.(start(i)>dimensions(i))) new_dimensions(i) = new_dimensions(i) + max(extend(i),start(i)-dimensions(i))
        end do
        if (any(new_dimensions /= dimensions)) then
            call h5dset_extent_f(dset, new_dimensions, ierr); BAIL(ierr)
            call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
        endif
    endif
    !
    ! Get a slice of the array
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, start-1, count, ierr); BAIL(ierr)
    call h5screate_simple_f(rank, count, memspace, ierr); BAIL(ierr)
    ! Write the data to the dataset
    call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, array, count, ierr, memspace, dspace); BAIL(ierr)
    !
    ! Close and release resources
    call h5sclose_f(memspace, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
    call h5sclose_f(dspace, ierr); BAIL(ierr)
  end function vh5_write_double_subarray_nd

  function vh5_write_double_array_nd(locid, dataset_name, rank, dimensions, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: dimensions(*)
    real(kind(0.d0)), intent(in) :: array(*)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset
    logical :: my_skip_present
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    if (vh5_exists(locid,dataset_name)) then
      if(my_skip_present) return
      ! Open dataset
      call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    else
      ! Create dataspace and dataset
      call h5screate_simple_f(rank, dimensions, dspace, ierr); BAIL(ierr)
      call h5dcreate_f(locid, dataset_name, H5T_IEEE_F64LE, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
      call h5sclose_f(dspace, ierr); BAIL(ierr)
    endif
    ! Write the data to the dataset
    call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, array, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5dclose_f(dset , ierr); BAIL(ierr)
  end function vh5_write_double_array_nd

  function vh5_write_double_array_1d(locid, dataset_name, array, numb, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(in) :: array(:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    integer, optional :: numb
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(1)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    dimensions(1) = size(array); if (present(numb)) dimensions(1) = numb
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    ierr = vh5_write_double_array_nd(locid, dataset_name, 1, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_double_array_1d

  function vh5_write_double_array_2d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(in) :: array(:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(2)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    dimensions = shape(array)
    ierr = vh5_write_double_array_nd(locid, dataset_name, 2, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_double_array_2d

  function vh5_write_double_array_3d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(in) :: array(:,:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(3)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    dimensions = shape(array)
    ierr = vh5_write_double_array_nd(locid, dataset_name, 3, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_double_array_3d

  function vh5_write_double_array_4d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(in) :: array(:,:,:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(4)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    dimensions = shape(array)
    ierr = vh5_write_double_array_nd(locid, dataset_name, 4, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_double_array_4d

  function vh5_write_double_array_5d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(in) :: array(:,:,:,:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(5)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    dimensions = shape(array)
    ierr = vh5_write_double_array_nd(locid, dataset_name, 5, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_double_array_5d

  function vh5_write_double_array_6d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(in) :: array(:,:,:,:,:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(6)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    dimensions = shape(array)
    ierr = vh5_write_double_array_nd(locid, dataset_name, 6, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_double_array_6d

  !---------------------------------------------------------
  ! read
  !---------------------------------------------------------
  function vh5_read_double_scalar(locid, dataset_name, x) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(out) :: x
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dset
    integer(HSIZE_T) :: dimensions(1) ! dummy, not used
    !---------------------------------------------------------
    !
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    ! read the data to the dataset
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, x, dimensions, ierr); BAIL(ierr)
    call h5dclose_f(dset , ierr); BAIL(ierr)
  end function vh5_read_double_scalar

  function vh5_read_double_subarray_nd(locid, dataset_name, rank, start, count, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: start(rank), count(rank)
    real(kind(0.d0)), intent(out) :: array(rank)
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, memspace, dset
    !---------------------------------------------------------
    !
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    ! Get a slice of the array
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, start-1, count, ierr); BAIL(ierr)
    call h5screate_simple_f(rank, count, memspace, ierr); BAIL(ierr)
    ! Read the data from the dataset
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, array, count, ierr, memspace, dspace); BAIL(ierr)
    ! Close and release resources
    call h5sclose_f(memspace, ierr); BAIL(ierr)
    call h5sclose_f(dspace, ierr); BAIL(ierr)
    call h5dclose_f(dset , ierr); BAIL(ierr)
  end function vh5_read_double_subarray_nd

  function vh5_read_double_array_nd(locid, dataset_name, rank, dimensions, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(out) :: dimensions(*)
    real(kind(0.d0)), intent(out) :: array(*)
    integer :: ierr, rank_or_error
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset
    integer(HSIZE_T) :: maxdims(rank)
    !---------------------------------------------------------
    !
    dimensions(:rank)=0
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    call h5sget_simple_extent_dims_f(dspace, dimensions, maxdims, rank_or_error)
    ierr = check_rank(rank, rank_or_error); BAIL(ierr)
    ! Read the data from the dataset
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, array, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5sclose_f(dspace, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_read_double_array_nd

  function vh5_read_double_array_1d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(out) :: array(:)
    integer, optional, intent(out) :: dimensions
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(1)
    !---------------------------------------------------------
    ierr = vh5_read_double_array_nd(locid, dataset_name, 1, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions(1)
  end function vh5_read_double_array_1d

  function vh5_read_double_array_2d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(out) :: array(:,:)
    integer, optional, intent(out) :: dimensions(2)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(2)
    !---------------------------------------------------------
    ierr = vh5_read_double_array_nd(locid, dataset_name, 2, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_double_array_2d

  function vh5_read_double_array_3d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(out) :: array(:,:,:)
    integer, optional, intent(out) :: dimensions(3)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(3)
    !---------------------------------------------------------
    ierr = vh5_read_double_array_nd(locid, dataset_name, 3, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_double_array_3d

  function vh5_read_double_array_4d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(out) :: array(:,:,:,:)
    integer, optional, intent(out) :: dimensions(4)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(4)
    !---------------------------------------------------------
    ierr = vh5_read_double_array_nd(locid, dataset_name, 4, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_double_array_4d

  function vh5_read_double_array_5d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(out) :: array(:,:,:,:,:)
    integer, optional, intent(out) :: dimensions(5)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(5)
    !---------------------------------------------------------
    ierr = vh5_read_double_array_nd(locid, dataset_name, 5, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_double_array_5d

  function vh5_read_double_array_6d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    real(kind(0.d0)), intent(out) :: array(:,:,:,:,:,:)
    integer, optional, intent(out) :: dimensions(6)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(6)
    !---------------------------------------------------------
    ierr = vh5_read_double_array_nd(locid, dataset_name, 6, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_double_array_6d

  !---------------------------------------------------------
  ! complex (double precision)
  !---------------------------------------------------------
  ! write
  !---------------------------------------------------------
  function vh5_create_double_complex_array_nd(locid, dataset_name, rank, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: dimensions(rank)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(rank+1) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:rank+1)= dimensions(1:rank)
    ierr = vh5_create_double_array_nd(locid, dataset_name, rank+1, rdimensions); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_create_double_complex_array_nd

  function vh5_write_double_complex_subarray_nd(locid, dataset_name, rank, start, count, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: start(rank), count(rank)
    complex(kind(0.d0)), intent(in), target :: array(*)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.d0)), pointer :: rarray(:)
    integer(HSIZE_T) :: rstart(rank+1), rcount(rank+1)
    !---------------------------------------------------------
    !
    rstart(1) = 1
    rstart(2:rank+1) = start(1:rank)
    rcount(1) = 2
    rcount(2:rank+1) = count(1:rank)
    call c_f_pointer(c_loc(array),rarray,[product(rcount)])
    ierr = vh5_write_double_subarray_nd(locid, dataset_name, rank+1, rstart, rcount, rarray)
  end function vh5_write_double_complex_subarray_nd

 !> @brief  write a complex array of ndimensions
 !>
 !> The strategy here is to cast the nD array to a 1D real array
 !> while still keeping the information about the dimensions
 !>
 function vh5_write_double_complex_array_nd(locid, dataset_name, rank, dimensions, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(out) :: dimensions(*)
    complex(kind(0.d0)), intent(in), target :: array(*)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.d0)), pointer :: rarray(:)
    integer(HSIZE_T) :: rdimensions(rank+1) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:rank+1)= dimensions(1:rank)
    call c_f_pointer(c_loc(array),rarray,[product(rdimensions)])
    ierr = vh5_write_double_array_nd(locid, dataset_name, rank+1, rdimensions, rarray); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_double_complex_array_nd

  function vh5_write_double_complex_array_1d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(in), target :: array(:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.d0)), pointer :: rarray(:,:)
    integer(HSIZE_T) :: rdimensions(2) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_double_array_2d(locid, dataset_name, rarray); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_double_complex_array_1d

  function vh5_write_double_complex_array_2d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(in), target :: array(:,:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.d0)), pointer :: rarray(:,:,:)
    integer(HSIZE_T) :: rdimensions(3) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_double_array_3d(locid, dataset_name, rarray); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_double_complex_array_2d

  function vh5_write_double_complex_array_3d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(in), target :: array(:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.d0)), pointer :: rarray(:,:,:,:)
    integer(HSIZE_T) :: rdimensions(4) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_double_array_4d(locid, dataset_name, rarray); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_double_complex_array_3d

  function vh5_write_double_complex_array_4d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(in), target :: array(:,:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.d0)), pointer :: rarray(:,:,:,:,:)
    integer(HSIZE_T) :: rdimensions(5) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_double_array_5d(locid, dataset_name, rarray)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_double_complex_array_4d

  function vh5_write_double_complex_array_5d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(in), target :: array(:,:,:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.d0)), pointer :: rarray(:,:,:,:,:,:)
    integer(HSIZE_T) :: rdimensions(6) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_double_array_6d(locid, dataset_name, rarray); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_double_complex_array_5d

  !---------------------------------------------------------
  ! read
  !---------------------------------------------------------
  function vh5_read_double_complex_subarray_nd(locid, dataset_name, rank, start, count, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: start(rank), count(rank)
    complex(kind(0.d0)), intent(in), target :: array(*)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.d0)), pointer :: rarray(:)
    integer(HSIZE_T) :: rstart(rank+1), rcount(rank+1)
    !---------------------------------------------------------
    !
    rstart(1) = 1
    rstart(2:rank+1) = start(1:rank)
    rcount(1) = 2
    rcount(2:rank+1) = count(1:rank)
    call c_f_pointer(c_loc(array),rarray,[product(rcount)])
    ierr = vh5_read_double_subarray_nd(locid, dataset_name, rank+1, rstart, rcount, rarray)
  end function vh5_read_double_complex_subarray_nd

  function vh5_read_double_complex_array_1d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(out), target :: array(:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(2) !> dimensions of real array
    real(kind(0.d0)), pointer :: rarray(:,:)
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_read_double_array_2d(locid, dataset_name, rarray)
  end function vh5_read_double_complex_array_1d

  function vh5_read_double_complex_array_2d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(out), target :: array(:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(3) !> dimensions of real array
    real(kind(0.d0)), pointer :: rarray(:,:,:)
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_read_double_array_3d(locid, dataset_name, rarray)
  end function vh5_read_double_complex_array_2d

  function vh5_read_double_complex_array_3d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(out), target :: array(:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(4) !> dimensions of real array
    real(kind(0.d0)), pointer :: rarray(:,:,:,:)
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_read_double_array_4d(locid, dataset_name, rarray)
  end function vh5_read_double_complex_array_3d

  function vh5_read_double_complex_array_5d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.d0)), intent(out), target :: array(:,:,:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(6) !> dimensions of real array
    real(kind(0.d0)), pointer :: rarray(:,:,:,:,:,:)
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_read_double_array_6d(locid, dataset_name, rarray)
  end function vh5_read_double_complex_array_5d

  !---------------------------------------------------------
  ! complex
  !---------------------------------------------------------
  ! write
  !---------------------------------------------------------
   function vh5_create_complex_array_nd(locid, dataset_name, rank, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: dimensions(rank)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(rank+1) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:rank+1)= dimensions(1:rank)
    ierr = vh5_create_real_array_nd(locid, dataset_name, rank+1, rdimensions); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_create_complex_array_nd

  function vh5_write_complex_subarray_nd(locid, dataset_name, rank, start, count, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: start(rank), count(rank)
    complex(kind(0.)), intent(in), target :: array(*)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.)), pointer :: rarray(:)
    integer(HSIZE_T) :: rstart(rank+1), rcount(rank+1)
    !---------------------------------------------------------
    !
    rstart(1) = 1
    rstart(2:rank+1) = start(1:rank)
    rcount(1) = 2
    rcount(2:rank+1) = count(1:rank)
    call c_f_pointer(c_loc(array),rarray,[product(rcount)])
    ierr = vh5_write_real_subarray_nd(locid, dataset_name, rank+1, rstart, rcount, rarray)
  end function vh5_write_complex_subarray_nd

 !> @brief  write a complex array of ndimensions
 !>
 !> The strategy here is to cast the nD array to a 1D real array
 !> while still keeping the information about the dimensions
 !>
 function vh5_write_complex_array_nd(locid, dataset_name, rank, dimensions, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(out) :: dimensions(*)
    complex(kind(0.0)), intent(in), target :: array(*)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.0)), pointer :: rarray(:)
    integer(HSIZE_T) :: rdimensions(rank+1) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:rank+1)= dimensions(1:rank)
    call c_f_pointer(c_loc(array),rarray,[product(rdimensions)])
    ierr = vh5_write_real_array_nd(locid, dataset_name, rank+1, rdimensions, rarray)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_complex_array_nd

  function vh5_write_complex_array_1d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.)), intent(in), target :: array(:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.)), pointer :: rarray(:,:)
    integer(HSIZE_T) :: rdimensions(2) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_real_array_2d(locid, dataset_name, rarray); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_complex_array_1d

  function vh5_write_complex_array_2d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.)), intent(in), target :: array(:,:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.)), pointer :: rarray(:,:,:)
    integer(HSIZE_T) :: rdimensions(3) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_real_array_3d(locid, dataset_name, rarray); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_complex_array_2d

  function vh5_write_complex_array_3d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.)), intent(in), target :: array(:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.)), pointer :: rarray(:,:,:,:)
    integer(HSIZE_T) :: rdimensions(4) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_real_array_4d(locid, dataset_name, rarray); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_complex_array_3d

  function vh5_write_complex_array_4d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.)), intent(in), target :: array(:,:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.)), pointer :: rarray(:,:,:,:,:)
    integer(HSIZE_T) :: rdimensions(5) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_real_array_5d(locid, dataset_name, rarray)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_complex_array_4d

  function vh5_write_complex_array_5d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.)), intent(in), target :: array(:,:,:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.)), pointer :: rarray(:,:,:,:,:,:)
    integer(HSIZE_T) :: rdimensions(6) !> dimensions of real array
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_write_real_array_6d(locid, dataset_name, rarray)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "complex")
  end function vh5_write_complex_array_5d


  !---------------------------------------------------------
  ! read
  !---------------------------------------------------------
  function vh5_read_complex_subarray_nd(locid, dataset_name, rank, start, count, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: start(rank), count(rank)
    complex(kind(0.)), intent(in), target :: array(*)
    integer :: ierr
    !---------------------------------------------------------
    real(kind(0.)), pointer :: rarray(:)
    integer(HSIZE_T) :: rstart(rank+1), rcount(rank+1)
    !---------------------------------------------------------
    !
    rstart(1) = 1
    rstart(2:rank+1) = start(1:rank)
    rcount(1) = 2
    rcount(2:rank+1) = count(1:rank)
    call c_f_pointer(c_loc(array),rarray,[product(rcount)])
    ierr = vh5_read_real_subarray_nd(locid, dataset_name, rank+1, rstart, rcount, rarray)
  end function vh5_read_complex_subarray_nd

  function vh5_read_complex_array_1d(locid, dataset_name, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.)), intent(out), target :: array(:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(2) !> dimensions of real array
    real(kind(0.)), pointer :: rarray(:,:)
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_read_real_array_2d(locid, dataset_name, rarray)
  end function vh5_read_complex_array_1d

  function vh5_read_complex_array_2d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.)), intent(out), target :: array(:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(3) !> dimensions of real array
    real(kind(0.)), pointer :: rarray(:,:,:)
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_read_real_array_3d(locid, dataset_name, rarray)
  end function vh5_read_complex_array_2d

  function vh5_read_complex_array_3d(locid, dataset_name, array) result(ierr)
    use iso_c_binding
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    complex(kind(0.)), intent(out), target :: array(:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: rdimensions(4) !> dimensions of real array
    real(kind(0.)), pointer :: rarray(:,:,:,:)
    !---------------------------------------------------------
    !
    rdimensions(1) = 2
    rdimensions(2:) = shape(array)
    call c_f_pointer(c_loc(array),rarray,rdimensions)
    ierr = vh5_read_real_array_4d(locid, dataset_name, rarray)
  end function vh5_read_complex_array_3d

  !---------------------------------------------------------
  ! integer
  !---------------------------------------------------------
  ! write
  !---------------------------------------------------------
  function vh5_write_integer_scalar(locid, dataset_name, x, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: x
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset
    integer(HSIZE_T) :: dimensions(1) ! dummy, not used
    logical :: my_skip_present
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    if (vh5_exists(locid,dataset_name)) then
      if(my_skip_present) return
      ! Open dataset
      call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    else
      ! Create dataspace and dataset
      call h5screate_f(H5S_SCALAR_F, dspace, ierr); BAIL(ierr)
      call h5dcreate_f(locid, dataset_name, H5T_STD_I32LE, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
      call h5sclose_f(dspace, ierr); BAIL(ierr)
    endif
    ! Write the data to the dataset
    call h5dwrite_f(dset, H5T_NATIVE_INTEGER, x, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_write_integer_scalar

  function vh5_write_integer_array_nd(locid, dataset_name, rank, dimensions, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: dimensions(*)
    integer, intent(in) :: array(*)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset
    logical :: my_skip_present
    !---------------------------------------------------------
    !
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    if (vh5_exists(locid,dataset_name)) then
      if(my_skip_present) return
      ! Open dataset
      call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    else
      ! Create dataspace and dataset
      call h5screate_simple_f(rank, dimensions, dspace, ierr); BAIL(ierr)
      call h5dcreate_f(locid, dataset_name, H5T_STD_I32LE, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
      call h5sclose_f(dspace, ierr); BAIL(ierr)
    endif
    ! Write the data to the dataset
    call h5dwrite_f(dset, H5T_NATIVE_INTEGER, array, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_write_integer_array_nd

  function vh5_write_integer_array_1d(locid, dataset_name, array, numb, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: array(:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    integer, optional :: numb
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(1)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    ! Check if explicit dimensions are passed
    dimensions(1) = size(array); if (present(numb)) dimensions(1) = numb
    ! Finally write
    ierr = vh5_write_integer_array_nd(locid, dataset_name, 1, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_integer_array_1d

  function vh5_write_integer_array_2d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: array(:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(2)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    dimensions = shape(array)
    ierr = vh5_write_integer_array_nd(locid, dataset_name, 2, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_integer_array_2d

  function vh5_write_integer_array_3d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: array(:,:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(3)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    dimensions = shape(array)
    ierr = vh5_write_integer_array_nd(locid, dataset_name, 3, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_integer_array_3d

  !---------------------------------------------------------
  ! read
  !---------------------------------------------------------
  function vh5_read_integer_scalar(locid, dataset_name, x) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(out) :: x
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dset
    integer(HSIZE_T) :: dimensions(1) ! dummy, not used
    !---------------------------------------------------------
    !
    ! Open the dataset
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    ! Read the data
    call h5dread_f(dset, H5T_NATIVE_INTEGER, x, dimensions, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_read_integer_scalar

  function vh5_read_integer_array_nd(locid, dataset_name, rank, dimensions, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(out) :: dimensions(rank)
    integer, intent(out) :: array(*)
    integer :: ierr, rank_or_error
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset
    integer(HSIZE_T) :: maxdims(rank) ! Array to store max dimension sizes
    !---------------------------------------------------------
    !
    dimensions = 0
    ! Open the dataset
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    call h5sget_simple_extent_dims_f(dspace, dimensions, maxdims, rank_or_error)
    ierr = check_rank(rank, rank_or_error); BAIL(ierr)
    ! Read the data
    call h5dread_f(dset, H5T_NATIVE_INTEGER, array, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5sclose_f(dspace, ierr); BAIL(ierr)
    call h5dclose_f(dset , ierr); BAIL(ierr)
  end function vh5_read_integer_array_nd

  function vh5_read_integer_array_1d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(out) :: array(:)
    integer, optional, intent(out) :: dimensions
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(1)
    !---------------------------------------------------------
    ierr = vh5_read_integer_array_nd(locid, dataset_name, 1, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions(1)
  end function vh5_read_integer_array_1d

  function vh5_read_integer_array_2d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(out) :: array(:,:)
    integer, optional, intent(out) :: dimensions(2)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(2)
    !---------------------------------------------------------
    ierr = vh5_read_integer_array_nd(locid, dataset_name, 2, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_integer_array_2d

  function vh5_read_integer_array_3d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(out) :: array(:,:,:)
    integer, optional, intent(out) :: dimensions(3)
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(3)
    !---------------------------------------------------------
    ierr = vh5_read_integer_array_nd(locid, dataset_name, 3, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions
  end function vh5_read_integer_array_3d

  !---------------------------------------------------------
  ! logical
  !
  ! Note that we can cast integer to logical as the fortran standard
  ! implies that the size of logical is equal to the size of integer
  !---------------------------------------------------------
  ! write
  !---------------------------------------------------------
  function vh5_write_logical_scalar(locid, dataset_name, x, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    logical, intent(in) :: x
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !-------------------------------------------------------
    logical :: my_skip_present
    !-------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    ! write stuff
    ierr = vh5_write_integer_scalar(locid, dataset_name, logical_to_integer(x), skip_present=my_skip_present); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "bool", skip_present=my_skip_present)
  end function vh5_write_logical_scalar

  function vh5_write_logical_array_1d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    logical, target, intent(in) :: array(:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer,pointer :: iarray(:)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    ! Create pointer
    call c_f_pointer(c_loc(array),iarray,shape(array))
    ! Write stuff
    ierr = vh5_write_integer_array_1d(locid, dataset_name, iarray, skip_present=my_skip_present); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "bool", skip_present=my_skip_present)
  end function vh5_write_logical_array_1d

  function vh5_write_logical_array_2d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    logical, target, intent(in) :: array(:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer,pointer :: iarray(:,:)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    ! Create pointer
    call c_f_pointer(c_loc(array),iarray,shape(array))
    ! Write stuff
    ierr = vh5_write_integer_array_2d(locid, dataset_name, iarray, skip_present=my_skip_present); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "bool", skip_present=my_skip_present)
  end function vh5_write_logical_array_2d

  function vh5_write_logical_array_3d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    logical, target, intent(in) :: array(:,:,:)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer,pointer :: iarray(:,:,:)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    ! Create pointer
    call c_f_pointer(c_loc(array),iarray,shape(array))
    ! Write stuff
    ierr = vh5_write_integer_array_3d(locid, dataset_name, iarray, skip_present=my_skip_present); BAIL(ierr)
    ierr = vh5_add_attribute_dataset(locid, dataset_name, "dtype", "bool", skip_present=my_skip_present)
  end function vh5_write_logical_array_3d

  !---------------------------------------------------------
  ! read
  !---------------------------------------------------------
  function vh5_read_logical_scalar(locid, dataset_name, x) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    logical, intent(out) :: x
    integer :: y, ierr
    !-------------------------------------------------------
    ierr = vh5_read_integer_scalar(locid, dataset_name, y)
    x = integer_to_logical(y)
  end function vh5_read_logical_scalar

  function vh5_read_logical_array_1d(locid, dataset_name, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    logical, target, intent(out) :: array(:)
    integer :: ierr
    !---------------------------------------------------------
    integer,pointer :: iarray(:)
    !---------------------------------------------------------
    call c_f_pointer(c_loc(array),iarray,shape(array))
    ierr = vh5_read_integer_array_1d(locid, dataset_name, iarray)
  end function vh5_read_logical_array_1d

  function vh5_read_logical_array_2d(locid, dataset_name, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    logical, target, intent(out) :: array(:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer,pointer :: iarray(:,:)
    !---------------------------------------------------------
    call c_f_pointer(c_loc(array),iarray,shape(array))
    ierr = vh5_read_integer_array_2d(locid, dataset_name, iarray)
  end function vh5_read_logical_array_2d

  function vh5_read_logical_array_3d(locid, dataset_name, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    logical, target, intent(out) :: array(:,:,:)
    integer :: ierr
    !---------------------------------------------------------
    integer,pointer :: iarray(:,:,:)
    !---------------------------------------------------------
    call c_f_pointer(c_loc(array),iarray,shape(array))
    ierr = vh5_read_integer_array_3d(locid, dataset_name, iarray)
  end function vh5_read_logical_array_3d

  !---------------------------------------------------------
  ! string (aka character array)
  !
  ! Note: I didn't implement writing string arrays, because without
  ! jagged arrays, I don't think they are really used in Fortran
  !---------------------------------------------------------
  ! write
  !---------------------------------------------------------
  function vh5_write_string_scalar(locid, dataset_name, str, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    character(LEN=*), intent(in) :: str
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset, string_type
    integer(HSIZE_T) :: dimensions(1) ! dummy, not used
    integer(HSIZE_T) :: strlen
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    strlen = len(str)
    if (strlen == 0) ierr = 1
    BAIL(ierr);
    ! Create datatype
    call h5tcopy_f(H5T_FORTRAN_S1, string_type, ierr); BAIL(ierr)
    call h5tset_size_f(string_type, strlen, ierr); BAIL(ierr)
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    if (vh5_exists(locid,dataset_name)) then
      if (my_skip_present) return
      ! Open dataset
      call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    else
      ! Create dataspace and dataset
      call h5screate_f(H5S_SCALAR_F, dspace, ierr); BAIL(ierr)
      call h5dcreate_f(locid, dataset_name, string_type, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
      call h5sclose_f(dspace, ierr); BAIL(ierr)
    endif
    ! Write the data to the dataset
    call h5dwrite_f(dset, string_type, str, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5tclose_f(string_type, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_write_string_scalar

  function vh5_write_string_array_nd(locid, dataset_name, rank, dimensions, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(in) :: dimensions(*)
    character(LEN=*), intent(in) :: array(*)
    logical, optional, intent(in) :: skip_present
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset, string_type
    integer(HSIZE_T) :: strlen
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    strlen = len(array(1))
    if (strlen == 0) ierr = 1
    BAIL(ierr);
    call h5tcopy_f(H5T_FORTRAN_S1, string_type, ierr); BAIL(ierr)
    call h5tset_size_f(string_type, strlen, ierr); BAIL(ierr)
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    if (vh5_exists(locid,dataset_name)) then
      if (my_skip_present) return
      ! Open dataset
      call h5dopen_f(locid, dataset_name, dset, ierr); BAIL(ierr)
    else
      ! Create dataspace dataset
      call h5screate_simple_f(rank, dimensions, dspace, ierr); BAIL(ierr)
      call h5dcreate_f(locid, dataset_name, string_type, dspace, dset, ierr, lcpl_id = create_inter_prop); BAIL(ierr)
      call h5sclose_f(dspace, ierr); BAIL(ierr)
    endif
    ! Write the data to the dataset
    call h5dwrite_f(dset, string_type, array, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5tclose_f(string_type, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_write_string_array_nd

  function vh5_write_string_array_1d(locid, dataset_name, array, skip_present) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    character(LEN=*), intent(in) :: array(:)
    integer :: ierr
    logical, optional, intent(in) :: skip_present
    !---------------------------------------------------------
    integer(HSIZE_T) :: dimensions(1)
    logical :: my_skip_present
    !---------------------------------------------------------
    ierr = 0; if (.not.lwriteh5) return
    ! Check if present
    my_skip_present = .false.; if (present(skip_present)) my_skip_present = skip_present
    if (vh5_exists(locid,dataset_name).and.my_skip_present) return
    ! Get dimensions and write data
    dimensions(1) = size(array)
    ierr = vh5_write_string_array_nd(locid, dataset_name, 1, dimensions, array, skip_present=my_skip_present)
  end function vh5_write_string_array_1d

  !---------------------------------------------------------
  ! read
  !---------------------------------------------------------
  function vh5_read_string_scalar(locid, dataset_name, str) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    character(LEN=:), allocatable, intent(out) :: str
    integer :: ierr
    !---------------------------------------------------------
    integer(HID_T) :: dset, file_type, string_type
    integer :: data_class, string_class
    integer(HSIZE_T) :: dimensions(1) ! dummy, not used
    integer(HSIZE_T) :: strlen
    !---------------------------------------------------------
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    ! Read the data to the dataset
    call H5dget_type_f(dset, file_type, ierr); BAIL(ierr)
    call h5tget_class_f(file_type, data_class, ierr); BAIL(ierr)
    call h5tget_class_f(H5T_FORTRAN_S1, string_class, ierr); BAIL(ierr)
    if (data_class /= string_class) then
      ierr = 5
    else
      call H5tget_size_f(file_type, strlen, ierr); BAIL(ierr)
      allocate(character(len=strlen) :: str)
      call h5tcopy_f(H5T_FORTRAN_S1, string_type, ierr); BAIL(ierr)
      call h5tset_size_f(string_type, strlen, ierr); BAIL(ierr)
      !
      call h5dread_f(dset, string_type, str, dimensions, ierr); BAIL(ierr)
    endif
    ! Close and release resources
    call h5dclose_f(dset, ierr); BAIL(ierr)
    call h5tclose_f(string_type, ierr); BAIL(ierr)
  end function vh5_read_string_scalar

  function vh5_read_string_array_nd(locid, dataset_name, rank, dimensions, array) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    integer, intent(in) :: rank
    integer(HSIZE_T), intent(out) :: dimensions(*)
    character(LEN=*), intent(out) :: array(*)
    integer :: ierr, rank_or_error
    !---------------------------------------------------------
    integer(HID_T) :: dspace, dset, string_type
    integer(HSIZE_T) :: strlen
    integer(HSIZE_T) :: maxdims(rank)
    !---------------------------------------------------------
    dimensions(:rank) = 0
    strlen = len(array(1))
    call h5tcopy_f(H5T_FORTRAN_S1, string_type, ierr); BAIL(ierr)
    call h5tset_size_f(string_type, strlen, ierr); BAIL(ierr)
    !
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    call h5sget_simple_extent_dims_f(dspace, dimensions, maxdims, rank_or_error)
    ierr = check_rank(rank, rank_or_error); BAIL(ierr)
    ! Read the data
    call h5dread_f(dset, string_type, array, dimensions, ierr); BAIL(ierr)
    ! Close and release resources
    call h5sclose_f(dspace, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
    call h5tclose_f(string_type, ierr); BAIL(ierr)
  end function vh5_read_string_array_nd

  function vh5_read_string_array_1d(locid, dataset_name, array, dimensions) result(ierr)
    integer(HID_T), intent(in) :: locid
    character(LEN=*), intent(in) :: dataset_name
    character(LEN=*), intent(out) :: array(:)
    integer, optional, intent(out) :: dimensions
    integer :: ierr
    !---------------------------------------------------------
    integer(HSIZE_T) :: my_dimensions(1)
    !---------------------------------------------------------
    ierr = vh5_read_string_array_nd(locid, dataset_name, 1, my_dimensions, array)
    if (present(dimensions)) dimensions = my_dimensions(1)
  end function vh5_read_string_array_1d

  !================================================
  !  Write file contents in HDF5
  !================================================

  !> @brief dump a file in text format to hdf5
  function vh5_dump_original_to_h5(groupid, subgroup, filename) result(ierr)
    use string, only: string_from_file
    integer(HID_T), intent(in)       :: groupid
    character(LEN=*), intent(in)    :: subgroup
    character(LEN=*), intent(in)    :: filename
    character(len=:), allocatable :: content
    integer(HID_T) :: subgroupid
    integer :: ierr
    content = string_from_file(filename, ierr); BAIL(ierr)
    VH5_CHECK(vh5_group_open_or_create(groupid, subgroup, subgroupid))
    VH5_CHECK(vh5_write(subgroupid, 'content', content))
    VH5_CHECK(vh5_group_close_writing(subgroupid))
    deallocate(content)
  end function vh5_dump_original_to_h5

  ! @brief function for retrieving potcar from h5
  function read_potcar_from_potcarh5(filenameh5, potentialtype, potcarname, sha256, content) result(ierr)

    character(len=*), intent(in) :: potentialtype
    character(len=*), intent(in) :: potcarname
    character(len=*), intent(inout) :: sha256
    character(len=:), allocatable :: content

    character(len=*) :: filenameh5
    ! needs to be larger by 1 --> 65 otherwise last character is truncated
    character(len=65) :: name_buffer ! buffer to hold object's name
    integer(HID_T) :: fileid, groupid, potid, contentid, datatype, string_type
    integer(SIZE_T) :: datasize, buf_size
    integer(HSIZE_T), dimension(1) :: dimensions ! dummy, not used
    integer :: n, i
    integer :: ih5err
    integer :: ierr
    integer :: nmembers, npotcars_found
    integer :: type
    logical :: exists
    integer :: length_user_sha256
    logical :: found_sha256

    BAIL(vh5_file_open_read(filenameh5, fileid))
    BAIL(vh5_group_open(fileid, potentialtype, groupid))
!    write(*,*) 'Search for '//potcarname
    call h5lexists_f(groupid, potcarname, exists, ih5err); VH5_CHECK(ih5err)
    if (.not.exists) then
      VH5_CHECK(vh5_group_close(groupid))
      VH5_CHECK(vh5_file_close(fileid))
      ierr = -1
      return
    end if
    call h5gopen_f(groupid, potcarname, potid, ih5err); VH5_CHECK(ih5err)
    ! check for SHA256
    !write(*,*)'CREATE_POTCAR, sha256: ', sha256
    if (sha256=='') then
      ! get latest attribute
      call h5lexists_f(potid, 'latest', exists, ih5err); VH5_CHECK(ih5err)
      if (.not.exists) then
        VH5_CHECK(vh5_group_close(potid))
        VH5_CHECK(vh5_group_close(groupid))
        VH5_CHECK(vh5_file_close(fileid))
        ierr = -1
        return
      end if
      call h5dopen_f(potid, 'latest', contentid, ih5err); VH5_CHECK(ih5err)
      buf_size = len(sha256)
      call h5gget_linkval_f(potid, 'latest', buf_size, sha256, ih5err); VH5_CHECK(ih5err)
    else
      found_sha256=.false.
      call h5gn_members_f(groupid, potcarname, nmembers, ih5err); VH5_CHECK(ih5err)
      npotcars_found = 0
      do i=0,nmembers-1
        call h5gget_obj_info_idx_f(potid, '.' , i, name_buffer, type, ih5err); VH5_CHECK(ih5err)
        length_user_sha256 = index(sha256, ' ') - 1
        if (sha256(1:length_user_sha256)==name_buffer(1:length_user_sha256)) then
          call h5dopen_f(potid, name_buffer, contentid, ih5err); VH5_CHECK(ih5err)
          npotcars_found = npotcars_found + 1
          found_sha256 = .true.
          sha256 = name_buffer(1:len(sha256)-1)
        end if
      enddo
      if (npotcars_found > 1) then
        VH5_CHECK(vh5_group_close(groupid))
        VH5_CHECK(vh5_file_close(fileid))
        ierr = -2
        return
      end if
      if (.not.found_sha256) then
        VH5_CHECK(vh5_group_close(groupid))
        VH5_CHECK(vh5_file_close(fileid))
        ierr = -1
        return
      end if
    end if

    call h5dget_type_f(contentid, datatype, ih5err); VH5_CHECK(ih5err)
    call h5tget_size_f(datatype, datasize, ih5err); VH5_CHECK(ih5err)
    if (allocated(content)) deallocate(content)
    allocate(character(datasize) :: content)
    call h5tcopy_f(H5T_FORTRAN_S1, string_type, ih5err); VH5_CHECK(ih5err)
    call h5tset_size_f(string_type, datasize, ih5err); VH5_CHECK(ih5err)
    call h5dread_f(contentid, string_type, content, dimensions, ih5err); VH5_CHECK(ih5err)
    call h5dclose_f(contentid, ih5err); VH5_CHECK(ih5err)
    ! close groups and file
    call h5gclose_f(potid, ih5err); VH5_CHECK(ih5err)
    call h5gclose_f(groupid, ih5err); VH5_CHECK(ih5err)
    ierr = vh5_file_close_writing(fileid)
  end function read_potcar_from_potcarh5

  function read_potcar_from_vaspin(content) result(ierr)
    character(len=:), allocatable :: content
    integer(HID_T) fileid, groupid, potid, contentid, datatype, string_type, dspace
    integer(SIZE_T) datasize
    integer(HSIZE_T), dimension(1) :: dimensions ! dummy, not used
    integer :: n, i
    integer :: ih5err
    integer :: ierr
    integer :: nmembers, npotcars_found
    integer :: type
    logical :: exists
    integer :: length_user_sha256
    logical :: found_sha256

    VH5_CHECK(vh5_group_open(ih5ininputgroup_id, subgrp_potcar, groupid))
    call h5dopen_f(groupid, 'content', contentid, ih5err); VH5_CHECK(ih5err)
    call h5dget_type_f(contentid, datatype, ih5err); VH5_CHECK(ih5err)
    call h5tget_size_f(datatype, datasize, ih5err); VH5_CHECK(ih5err)
    if (allocated(content)) deallocate(content)
    allocate(character(datasize) :: content)
    call h5tcopy_f(H5T_FORTRAN_S1, string_type, ih5err); VH5_CHECK(ih5err)
    call h5tset_size_f(string_type, datasize, ih5err); VH5_CHECK(ih5err)
    call h5dread_f(contentid, string_type, content, dimensions, ih5err); VH5_CHECK(ih5err)
    call h5dclose_f(contentid, ih5err); VH5_CHECK(ih5err)
    call h5gclose_f(groupid, ih5err); VH5_CHECK(ih5err)

  end function read_potcar_from_vaspin

  integer function vh5_number_elements(locid, dataset_name, number_elements) result (ierr)
    integer(HID_T), intent(in) :: locid
    character(len=*), intent(in) :: dataset_name
    integer, intent(out) :: number_elements
    integer(HSIZE_T) npoints
    integer(HID_T) dspace, dset
    number_elements = 0
    ierr = vh5_open_dataset_ifexists(locid, dataset_name, dset); BAIL(ierr)
    call h5dget_space_f(dset, dspace, ierr); BAIL(ierr)
    call h5sget_simple_extent_npoints_f(dspace, npoints, ierr); BAIL(ierr)
    number_elements = npoints
    call h5sclose_f(dspace, ierr); BAIL(ierr)
    call h5dclose_f(dset, ierr); BAIL(ierr)
  end function vh5_number_elements

end module vhdf5_base
#else
!--------------------------------------------------------------------
! Dummy module
!--------------------------------------------------------------------
module vhdf5_base
  contains
  subroutine vhdf5_base_dummy
  write(*,*)'Im a DEC compiler so I need this line'
  end subroutine
end module vhdf5_base
#endif
