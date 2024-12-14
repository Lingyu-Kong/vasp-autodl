#include "symbol.inc"
#ifdef VASP_HDF5

module vhdf5
  use prec
  use poscar, only : type_info, dynamics
  use base, only : info_struct
  use wave, only: wavedes, wavespin
  use lattice, only : latt, dirkar
  use mkpoints, only : kpoints_struct
  use fileio, only : closewav, write_from_buf
  use vhdf5_base

  implicit none

  integer,parameter :: MAX_LEN_GROUP=50

  public :: start_h5
  public :: stop_h5
  public :: vh5_write_dosheader
  public :: vh5_write_eigenvalheader
  public :: vh5_write_eigenval
  public :: vh5_write_dos
  public :: vh5_write_lattice_ion
  public :: vh5_write_lattice_ion_hist_start
  public :: vh5_write_lattice_ion_hist_step
  public :: vh5_write_energies_hist_start
  public :: vh5_write_energies_hist_step
  public :: vh5_write_stress_forces_hist_start
  public :: vh5_write_stress_forces_hist_step
  public :: vh5_write_positions_trailer
  public :: vh5_write_charge
  public :: vh5_write_paw_occupancies
  public :: vh5_write_atommom
  public :: vh5_write_wavefunctions
  public :: vh5_read_wavefunction_header
  public :: vh5_read_wavefunctions
  public :: vh5_read_force_constants
  public :: vh5_write_projectors
  public :: vh5_dump_original_to_h5


contains

  subroutine start_h5(subdir, vaspin, vaspout)
    character(len=*), intent(in), optional :: subdir, vaspin, vaspout
    character(len=:), allocatable :: subdir_, vaspin_, vaspout_
    logical :: lh5

    VH5_CHECK(vh5_start())
    !
    ! incar, vaspin
    !
    if (hdf5_found) then
      if (present(vaspin)) then
        vaspin_ = vaspin
      else
        vaspin_ = 'vaspin.h5'
      end if
      VH5_CHECK(vh5_file_open_read(vaspin_, ih5infileid))
      VH5_CHECK(vh5_group_open(ih5infileid, grp_input, ih5ininputgroup_id))
    endif
    !
    ! write vaspout.h5
    !
    ! vh5_start() initializes the hdf5 library and needs to be called before using any
    ! hdf5 related functionality
    if (present(vaspout)) then
      vaspout_ = vaspout
    else
      vaspout_ = 'vaspout.h5'
    end if
    if (present(subdir)) then
      subdir_ = subdir
    else
      subdir_ = ''
    end if
    VH5_CHECK(vh5_file_create_or_overwrite(subdir_ // vaspout_, ih5outfileid))
    VH5_CHECK(vh5_group_open_or_create(ih5outfileid, grp_input, ih5outinputgroup_id))
    VH5_CHECK(vh5_group_open_or_create(ih5outfileid, grp_intermediate, ih5intermediategroup_id))
    call vh5_write_version(ih5outfileid)
  end subroutine

  subroutine stop_h5(subdir)
    character(len=*), intent(in), optional :: subdir
    character(len=:), allocatable :: subdir_
    logical image_incar
    !
    ! if no INCAR, try to close vaspin.h5
    !
    if (.not. incar_found .and. hdf5_found) then
      VH5_CHECK(vh5_group_close(ih5ininputgroup_id))
      VH5_CHECK(vh5_file_close(ih5infileid))
    endif
    !
    ! write original  input data to vaspout.h5
    !
    if (present(subdir)) then
      subdir_ = subdir
    else
      subdir_ = ''
    end if
    if (poscar_found.or.incar_found.or.kpoints_found) then
      VH5_CHECK(vh5_group_open_or_create(ih5outfileid, grp_original, ih5outoriginalgroup_id))
    endif
    if (poscar_found) then
      VH5_CHECK(vh5_dump_original_to_h5(ih5outoriginalgroup_id, subgrp_poscar, subdir_ // "POSCAR"))
    endif
    if (kpoints_found) then
      VH5_CHECK(vh5_dump_original_to_h5(ih5outoriginalgroup_id, subgrp_kpoints, subdir_ // "KPOINTS"))
    endif
    if (incar_found) then
      inquire(file=subdir_ // 'INCAR', exist=image_incar)
      if (image_incar) then
        VH5_CHECK(vh5_dump_original_to_h5(ih5outoriginalgroup_id, subgrp_incar, subdir_ // "INCAR"))
      else
        VH5_CHECK(vh5_dump_original_to_h5(ih5outoriginalgroup_id, subgrp_incar, "INCAR"))
      end if
    endif
    if (poscar_found.or.incar_found.or.kpoints_found) then
      VH5_CHECK(vh5_group_close_writing(ih5outoriginalgroup_id))
    endif
    !
    ! close vaspout.h5
    !
    VH5_CHECK(vh5_group_close_writing(ih5outinputgroup_id))
    VH5_CHECK(vh5_group_close_writing(ih5intermediategroup_id))
    VH5_CHECK(vh5_file_close_writing(ih5outfileid))
    ! vh5_end() frees memory allocated by vh5_start() and must be called at the end
    ! of any program using hdf5
    VH5_CHECK(vh5_end())
  end subroutine stop_h5

  !> @brief write Vasp version into hdf5 file
  subroutine vh5_write_version(fileid)
    use version, only: major, minor, patch
    integer(HID_T), intent(in) :: fileid
    integer(HID_T) :: groupid
    VH5_CHECK(vh5_group_open_or_create(fileid, grp_version, groupid))
    VH5_CHECK(vh5_write(groupid, "major", major))
    VH5_CHECK(vh5_write(groupid, "minor", minor))
    VH5_CHECK(vh5_write(groupid, "patch", patch))
    VH5_CHECK(vh5_group_close_writing(groupid))
  end subroutine vh5_write_version

  !> @brief write header of doscar files.
  subroutine vh5_write_dosheader(fileid, typeinfo, jobpar, mdynamics, wave, lattice, info, kpoints, group, subgroup)
    integer, intent(in)              :: jobpar !< if io%lorbit > 10 then jobpar is 1
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    type(type_info), intent(in)      :: typeinfo !< typeinfo data structure (atom types, number of ions, etc.)
    type(dynamics), intent(in)       :: mdynamics !< MD data structure (atom positions, velocities, etc.)
    type(wavedes), intent(in)        :: wave !< wave function description data
    type(latt), intent(in)           :: lattice !< lattice description
    type(info_struct), intent(in)    :: info !< info structure
    type(kpoints_struct), intent(in) :: kpoints !< kpoints data
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_DOS)
    !-------------------------------------------------------------------------------------
    real(q) :: aomega
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_DOS;      if (present(subgroup)) my_subgroup = trim(subgroup)
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    aomega = lattice%omega/typeinfo%nions
    VH5_CHECK(vh5_write(subgroupid, "nions", typeinfo%nions))
    VH5_CHECK(vh5_write(subgroupid, "jobpar", jobpar))
    VH5_CHECK(vh5_write(subgroupid, "aomega", aomega))
    VH5_CHECK(vh5_write(subgroupid, "anorm", lattice%anorm*1E-10_q))
    VH5_CHECK(vh5_write(subgroupid, "potim", mdynamics%potim*1e-15_q))
    VH5_CHECK(vh5_write(subgroupid, "temperature", mdynamics%temp))
    VH5_CHECK(vh5_write(subgroupid, "system", info%sznam1))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
  end subroutine vh5_write_dosheader

  !> @brief write header of eigenval files.
  subroutine vh5_write_eigenvalheader(fileid, typeinfo, mdynamics, wave, lattice, info, kpoints, group, subgroup)
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    type(type_info), intent(in)      :: typeinfo !< typeinfo data structure (atom types, number of ions, etc.)
    type(dynamics), intent(in)       :: mdynamics !< MD data structure (atom positions, velocities, etc.)
    type(wavedes), intent(in)        :: wave !< wave function description data
    type(latt), intent(in)           :: lattice !< lattice description
    type(info_struct), intent(in)    :: info !< info structure
    type(kpoints_struct), intent(in) :: kpoints !< kpoints data
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_EIGENVAL)
    !-------------------------------------------------------------------------------------
    real(q) :: aomega
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------

    my_group = GRP_RESULTS;     if (present(group))    my_group = group
    my_subgroup = GRP_EIGENVAL; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    aomega = lattice%omega/typeinfo%nions
    VH5_CHECK(vh5_write(subgroupid, "nions", typeinfo%nions))
    VH5_CHECK(vh5_write(subgroupid, "nblocks", mdynamics%nblock * mdynamics%kblock))
    VH5_CHECK(vh5_write(subgroupid, "ispin", wave%ispin))
    VH5_CHECK(vh5_write(subgroupid, "aomega", aomega))
    VH5_CHECK(vh5_write(subgroupid, "anorm", lattice%anorm*1E-10_q))
    VH5_CHECK(vh5_write(subgroupid, "potim", mdynamics%potim*1e-15_q))
    VH5_CHECK(vh5_write(subgroupid, "temperature", mdynamics%temp))
    VH5_CHECK(vh5_write(subgroupid, "system", info%sznam1))
    VH5_CHECK(vh5_write(subgroupid, "nelectrons", info%nelect))
    VH5_CHECK(vh5_write(subgroupid, "kpoints", kpoints%nkpts))
    VH5_CHECK(vh5_write(subgroupid, "nb_tot", wave%nb_tot))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
  end subroutine vh5_write_eigenvalheader

  !> @brief write electronic energies
  subroutine vh5_write_eigenval(fileid, wave, wavesp, kpoints, group, subgroup)
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    type(wavedes), intent(in)        :: wave !< wave function description data
    type(wavespin), intent(in)       :: wavesp !< wavefunctions including band index and spin
    type(kpoints_struct), intent(in) :: kpoints !< kpoints data
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_EIGENVAL)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_EIGENVAL; if (present(subgroup)) my_subgroup = trim(subgroup)
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "kpoint_coords", kpoints%vkpt))
    VH5_CHECK(vh5_write(subgroupid, "kpoints_symmetry_weight", kpoints%wtkpt))
    VH5_CHECK(vh5_write(subgroupid, "eigenvalues", real(wavesp%celtot)))
    VH5_CHECK(vh5_write(subgroupid, "fermiweights", wavesp%fertot))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_eigenval

  !> @brief write electronic energies
  subroutine vh5_write_band_velocities(fileid, velocities, group, subgroup)
    integer(HID_T), intent(in) :: fileid !< hdf5 file handle
    real(q), intent(in) :: velocities(:,:,:,:) !< band velocities
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_EIGENVAL)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_EIGENVAL; if (present(subgroup)) my_subgroup = trim(subgroup)
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "velocities", velocities))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_band_velocities

  !> @brief write phonon frequencies and eigenvectors
  subroutine vh5_write_phonons(fileid, q_mesh, frequ, evecs, group, subgroup)
    integer(HID_T), intent(in)        :: fileid !< hdf5 file handle
    type (kpoints_struct), intent(in) :: q_mesh
    real(q), intent(in) :: frequ(:, :)
    complex(q), intent(in) :: evecs(:, :, :)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_EIGENVAL)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    integer :: nmodes
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_PHONONS; if (present(subgroup)) my_subgroup = trim(subgroup)
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "kpoint_coords", q_mesh%vkpt))
    VH5_CHECK(vh5_write(subgroupid, "kpoints_symmetry_weight", q_mesh%wtkpt))

    nmodes = size(frequ,1)
    VH5_CHECK(vh5_write(subgroupid, "nions", nmodes/3))
    VH5_CHECK(vh5_write(subgroupid, "frequencies", frequ))
    VH5_CHECK(vh5_write(subgroupid, "eigenvectors", evecs))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_phonons

  !> @brief write electronic density of states
  subroutine vh5_write_dos(fileid, wave, kpoints, dos, dosi, dospar, efermi, nionp, lpar, group, subgroup, scale)
    integer(HID_T), intent(in)       :: fileid !< hdf5 file handle
    type(wavedes), intent(in)        :: wave !< wave function description data
    type(kpoints_struct), intent(in) :: kpoints !< kpoints data
    real(q), intent(in) :: dos(:,:) !< density of states
    real(q), intent(in) :: dosi(:,:) !< integrated dos
    real(q), intent(in) :: dospar(:,:,:,:) !< partial density of states
    real(q), intent(in) :: efermi !< fermi energy
    integer, intent(in) :: nionp !< number of ions
    integer, intent(in) :: lpar
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_DOS)
    real(q), intent(in),optional :: scale !< scale dos and idos before writing
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    integer :: i, nedos
    real(q) :: my_scale, deltae
    real(q) :: energies(size(dos,1))
    character(len=MAX_LEN_GROUP) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------

    my_group    = GRP_RESULTS;  if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_DOS;      if (present(subgroup)) my_subgroup = trim(subgroup)
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    my_scale    = 1.d0; if (present(scale)) my_scale = scale
    nedos = size(dos,1)

    VH5_CHECK(vh5_write(subgroupid, "nionp", nionp))
    VH5_CHECK(vh5_write(subgroupid, "ncdij", wave%ncdij))
    VH5_CHECK(vh5_write(subgroupid, "lpar", lpar))

    VH5_CHECK(vh5_write(subgroupid, "emax", kpoints%emax))
    VH5_CHECK(vh5_write(subgroupid, "emin", kpoints%emin))
    VH5_CHECK(vh5_write(subgroupid, "nedos", nedos))
    VH5_CHECK(vh5_write(subgroupid, "efermi", efermi))

    deltae=(kpoints%emax-kpoints%emin)/(nedos-1)
    do i=1,nedos
        energies(i)=kpoints%emin+deltae*(i-1)
    end do
    VH5_CHECK(vh5_write(subgroupid, "energies", energies))
    VH5_CHECK(vh5_write(subgroupid, "dos", dos*my_scale))
    VH5_CHECK(vh5_write(subgroupid, "dosi", dosi*my_scale))
    VH5_CHECK(vh5_write(subgroupid, "dospar", dospar))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_dos

  !> @brief write lattice parameters and positions
  subroutine vh5_write_lattice_ion(fileid, sznam, typeinfo, scale, amat, selective_dynamics, posions, group, subgroup)
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    character(LEN=*), intent(in)  :: sznam !< string with system name read from poscar
    type(type_info), intent(in)   :: typeinfo !< typeinfo data structure (atom types, number of ions, etc.)
    real(q), intent(in) :: scale !< scale for the lattice vectors amat
    real(q), intent(in) :: amat(3,3) !< lattice vectors
    logical :: selective_dynamics !< logical determining whether we are running selective dynamics
    real(q), dimension(3,typeinfo%nions), intent(in) :: posions !< ion positions
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_POSITIONS)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: posgroupid
    integer(HID_T) :: unitgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;      if (present(group))    my_group = group
    my_subgroup = GRP_POSITIONS; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), unitgroupid))
    VH5_CHECK(vh5_group_open_or_create(unitgroupid, trim(my_subgroup), posgroupid))

    ! lattice and system information
    VH5_CHECK(vh5_write(posgroupid, "system", sznam))
    VH5_CHECK(vh5_write(posgroupid, "scale", scale))
    VH5_CHECK(vh5_write(posgroupid, "lattice_vectors", amat/scale))
    VH5_CHECK(vh5_write(posgroupid, "direct_coordinates", 1))
    VH5_CHECK(vh5_write(posgroupid, "ion_types", typeinfo%type(1:typeinfo%ntyp)))
    VH5_CHECK(vh5_write(posgroupid, "ion_sha256", typeinfo%sha256(1:typeinfo%ntyp)))
    VH5_CHECK(vh5_write(posgroupid, "number_ion_types", typeinfo%nityp(1:typeinfo%ntyp)))

    !positions of the ions
    VH5_CHECK(vh5_write(posgroupid, "position_ions", posions))

    !selective dynamics
    VH5_CHECK(vh5_write(posgroupid, "selective_dynamics", selective_dynamics))
    if (selective_dynamics) then
       VH5_CHECK(vh5_write(posgroupid, "selective_dynamics_ions", typeinfo%lsfor))
    end if

    VH5_CHECK(vh5_group_close_writing(posgroupid))
    VH5_CHECK(vh5_group_close_writing(unitgroupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_lattice_ion

  !> @brief start writting of lattice vectors and ion positions history to hdf5 file
  subroutine vh5_write_lattice_ion_hist_start(fileid, subgroup, nions)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer,intent(in) :: nions
    !-------------------------------------------------------------
    integer(HID_T) :: groupid
    integer(SIZE_T) :: dims(3), maxdims(3)
    !-------------------------------------------------------------
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), groupid))

    dims = [3,nions,1]; maxdims = [integer(SIZE_T)::3,nions,H5S_UNLIMITED_F]
    VH5_CHECK(vh5_create_double_array_nd(groupid, "position_ions", 3, dims, maxdims))

    dims = [3,3,1]; maxdims = [integer(SIZE_T)::3,3,H5S_UNLIMITED_F]
    VH5_CHECK(vh5_create_double_array_nd(groupid, "lattice_vectors", 3, dims, maxdims))

    VH5_CHECK(vh5_group_close_writing(groupid))
  end subroutine vh5_write_lattice_ion_hist_start

  !> @brief write the lattice vectors and ion positions for one ionic step to hdf5 file
  subroutine vh5_write_lattice_ion_hist_step(fileid, subgroup, nstep, amat, posions)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in)          :: nstep
    real(q), intent(in) :: amat(3,3)
    real(q), intent(in) :: posions(:,:)
    !-------------------------------------------------------------------------------------
    integer :: nions
    integer(HID_T) :: groupid
    integer(SIZE_T) :: start(3), count(3), extend(3)
    !-------------------------------------------------------------------------------------
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), groupid))

    nions = size(posions,2)
    start = [1,1,nstep]
    count = [3,nions,1]
    extend = [0,0,1]
    VH5_CHECK(vh5_write_double_subarray_nd(groupid, "position_ions", 3, start, count, posions, extend))

    count = [3,3,1]
    VH5_CHECK(vh5_write_double_subarray_nd(groupid, "lattice_vectors", 3, start, count, amat, extend))

    VH5_CHECK(vh5_group_close_writing(groupid))
  end subroutine vh5_write_lattice_ion_hist_step

  !> @brief start writing of enegies history to hdf5 file
  subroutine vh5_write_energies_hist_start(fileid, subgroup, ibrion)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in) :: ibrion
    !-------------------------------------------------------------
    integer :: extent
    integer(HID_T) :: groupid
    integer(SIZE_T) :: dims(2), maxdims(2)
    character(LEN=24) :: energies_tags(7)
    !-------------------------------------------------------------
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), groupid))

    if (ibrion == 0) then
        energies_tags = [ 'ion-electron   TOTEN   ',&
                          'kinetic energy EKIN    ',&
                          'kin. lattice   EKIN_LAT',&
                          'temperature    TEIN    ',&
                          'nose potential ES      ',&
                          'nose kinetic   EPS     ',&
                          'total energy   ETOTAL  ' ]
        VH5_CHECK(vh5_write(groupid, "energies_tags", energies_tags))
        extent = 7
    else
        energies_tags(1:3) = [ 'free energy    TOTEN   ',&
                               'energy without entropy ',&
                               'energy(sigma->0)       ']
        VH5_CHECK(vh5_write(groupid, "energies_tags", energies_tags(1:3)))
        extent = 3
    endif

    ! Create an extensible array that will contain all the energies
    dims = [extent,1]; maxdims = [integer(SIZE_T)::extent,H5S_UNLIMITED_F]
    VH5_CHECK(vh5_create_double_array_nd(groupid, "energies", 2, dims, maxdims))

    VH5_CHECK(vh5_group_close_writing(groupid))
  end subroutine vh5_write_energies_hist_start

  !> @brief write the energies for one ionic step to hdf5 file
  subroutine vh5_write_energies_hist_step(fileid, subgroup, nstep, energies)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in)          :: nstep
    real(q), intent(in) :: energies(:)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid
    integer(SIZE_T) :: start(2), count(2), extend(2)
    !-------------------------------------------------------------------------------------
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), groupid))

    start = [1,nstep]
    count = [size(energies),1]
    extend = [0,1]
    VH5_CHECK(vh5_write_double_subarray_nd(groupid, "energies", 2, start, count, energies, extend))

    VH5_CHECK(vh5_group_close_writing(groupid))
  end subroutine vh5_write_energies_hist_step

  !> @brief start writting of stress and forces history to hdf5 file
  subroutine vh5_write_stress_forces_hist_start(fileid, subgroup, nions)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer,intent(in) :: nions
    !-------------------------------------------------------------
    integer(HID_T) :: posgroupid
    integer(SIZE_T) :: dims(3), maxdims(3)
    !-------------------------------------------------------------
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), posgroupid))

    dims = [3,nions,1]; maxdims = [integer(SIZE_T)::3,nions,H5S_UNLIMITED_F]
    VH5_CHECK(vh5_create_double_array_nd(posgroupid, "forces", 3, dims, maxdims))

    dims = [3,3,1]; maxdims = [integer(SIZE_T)::3,3,H5S_UNLIMITED_F]
    VH5_CHECK(vh5_create_double_array_nd(posgroupid, "stress", 3, dims, maxdims))

    VH5_CHECK(vh5_group_close_writing(posgroupid))
  end subroutine vh5_write_stress_forces_hist_start

  !> @brief write the stress and the forces for one ionic step to hdf5 file
  subroutine vh5_write_stress_forces_hist_step(fileid, subgroup, nstep, stress, forces)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in)          :: nstep
    real(q), intent(in) :: stress(3,3)
    real(q), intent(in) :: forces(:,:)
    !-------------------------------------------------------------------------------------
    integer :: nions
    integer(HID_T) :: posgroupid
    integer(SIZE_T) :: start(3), count(3), extend(3)
    !-------------------------------------------------------------------------------------
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), posgroupid))

    nions = size(forces,2)
    start = [1,1,nstep]
    count = [3,nions,1]
    extend = [0,0,1]
    VH5_CHECK(vh5_write_double_subarray_nd(posgroupid, "forces", 3, start, count, forces, extend))

    count = [3,3,1]
    VH5_CHECK(vh5_write_double_subarray_nd(posgroupid, "stress", 3, start, count, stress, extend))

    VH5_CHECK(vh5_group_close_writing(posgroupid))
  end subroutine vh5_write_stress_forces_hist_step

  !> @brief start writing of charges and magnetic moments history to hdf5 file
  subroutine vh5_write_magnetism_hist_start(fileid, subgroup, shape_)
    integer(HID_T), intent(in)   :: fileid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in) :: shape_(3)
    !-------------------------------------------------------------
    integer(HID_T) :: groupid
    integer(SIZE_T) :: dims(4), maxdims(4)
    character(len=17), parameter :: component_tags(4) = ['total charge     ', &
        'magnetization (x)', 'magnetization (y)', 'magnetization (z)']
    character(len=1), parameter :: orbital_tags(4) = ['s', 'p', 'd', 'f']
    !-------------------------------------------------------------
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), groupid))
    VH5_CHECK(vh5_write(groupid, "magnetism/component_tags", component_tags(:shape_(3))))
    VH5_CHECK(vh5_write(groupid, "magnetism/orbital_tags", orbital_tags(:shape_(1))))

    ! Create an extensible array that will contain all the charges and magnetic moments
    dims = [shape_, 1]; maxdims = [int(shape_, kind=SIZE_T), H5S_UNLIMITED_F]
    VH5_CHECK(vh5_create_double_array_nd(groupid, "magnetism/moments", size(dims), dims, maxdims))

    VH5_CHECK(vh5_group_close_writing(groupid))
  end subroutine vh5_write_magnetism_hist_start

  !> @brief write the charges and magnetic moments for one ionic step to hdf5 file
  subroutine vh5_write_magnetism_hist_step(locid, subgroup, nstep, magnetic_moments)
    integer(HID_T), intent(in)   :: locid
    character(LEN=*), intent(in) :: subgroup
    integer, intent(in)          :: nstep
    real(q), intent(in) :: magnetic_moments(:,:,:)
    !-------------------------------------------------------------------------------------
    integer ierr
    integer(HID_T) :: groupid
    integer(SIZE_T) :: start(4), count_(4), extend(4)
    !-------------------------------------------------------------------------------------
    VH5_CHECK(vh5_group_open_or_create(locid, trim(subgroup), groupid))

    start = [1,1,1,nstep]
    count_ = [shape(magnetic_moments),1]
    extend = [0,0,0,1]
    ierr = vh5_write_double_subarray_nd(groupid, "magnetism/moments", size(start), &
        start, count_, magnetic_moments, extend); VH5_CHECK(ierr)

    VH5_CHECK(vh5_group_close_writing(groupid))
  end subroutine vh5_write_magnetism_hist_step

  !> @brief write contcar trailer
  !>
  !> writes the informations contained in the contcar trailer in the "traditional"
  !> VASP output files to the hdf5 file. Writes lattice velocities and velocities
  !>
  !> @param fileid hdf5 file handle
  !> @param lattice lattice description (def. in lattice.F)
  !> @param typeinfo typeinfo data structure (atom types, number of ions, etc. defined in
  !>        poscar.f)
  !> @param mdynamics MD data structure (atom positions, velocities, etc. defined in poscar.f)
  subroutine vh5_write_positions_trailer(fileid, lattice, typeinfo, mdynamics, group, subgroup)
    integer(HID_T), intent(in) :: fileid
    type(latt), intent(in):: lattice
    type(type_info), intent(in):: typeinfo
    type(dynamics), intent(in):: mdynamics
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_POSITIONS)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    integer:: ni
    real(q), dimension(3,typeinfo%nions) :: velocities_out
    real(q), dimension(3) :: tmp
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;      if (present(group))    my_group = group
    my_subgroup = GRP_POSITIONS; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    if (mdynamics%ibrion==0 .and. mdynamics%isif==3) then
       VH5_CHECK(vh5_write(subgroupid, "lattice_velocities", lattice%avel/mdynamics%potim))
    end if

    if (typeinfo%nionp > typeinfo%nions) then
       VH5_CHECK(vh5_write(subgroupid, "number_empty_sphere_types", typeinfo%nityp(typeinfo%ntyp+1:typeinfo%ntypp)))
       VH5_CHECK(vh5_write(subgroupid, "position_empty_spheres", mdynamics%posion(:,typeinfo%nions+1:typeinfo%nionp)))
    endif

    if (mdynamics%ibrion/=44 .and. mdynamics%ibrion/=40) then
       do ni = 1, typeinfo%nions
          tmp = mdynamics%vel(:,ni)/mdynamics%potim
          call dirkar(1,tmp,lattice%a)
          velocities_out(:,ni) = tmp
       end do
    else
       velocities_out = mdynamics%vel
    end if
    VH5_CHECK(vh5_write(subgroupid, "ion_velocities", velocities_out))
    VH5_CHECK(vh5_write(subgroupid, "direct_coordinates_velocities", 0))
    if (mdynamics%init==1) then
       VH5_CHECK(vh5_write(subgroupid, "dyn_init", mdynamics%init))
       VH5_CHECK(vh5_write(subgroupid, "potim", mdynamics%potim))
       VH5_CHECK(vh5_write(subgroupid, "nose_thermostat", mdynamics%snose))
       VH5_CHECK(vh5_write(subgroupid, "predictor_coordinates", mdynamics%posion(:,1:typeinfo%nions)))
       VH5_CHECK(vh5_write(subgroupid, "predictor_coordinates_2", mdynamics%d2))
       VH5_CHECK(vh5_write(subgroupid, "predictor_coordinates_3", mdynamics%d3))
    end if

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_positions_trailer

  !> @brief write force constants
  subroutine vh5_write_force_constants(fileid, force_constants, group, subgroup)
    integer(HID_T), intent(in)       :: fileid
    real(q) :: force_constants(:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;           if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "force_constants", force_constants))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_force_constants

  !> @brief write hessian matrix
  subroutine vh5_write_hessian(fileid, hessian, group, subgroup)
    integer(HID_T), intent(in)       :: fileid
    real(q) :: hessian(:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "hessian", hessian))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_hessian

  !> @brief write born effective charges
  subroutine vh5_write_elastic_modulus(fileid, prefix, elastic_modulus, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    character(len=*), intent(in) :: prefix
    real(q),intent(in) :: elastic_modulus(:,:,:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, prefix // "_elastic_modulus", elastic_modulus))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
  end subroutine vh5_write_elastic_modulus

  !> @brief write born effective charges
  subroutine vh5_write_born_charges(fileid, born_charges, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    real(q),intent(in) :: born_charges(:,:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "born_charges", born_charges))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
  end subroutine vh5_write_born_charges

  !> @brief write internal strain tensor
  subroutine vh5_write_internal_strain(fileid, internal_strain, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    real(q),intent(in) :: internal_strain(:,:,:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "internal_strain", internal_strain))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
  end subroutine vh5_write_internal_strain

  !> @brief write piezoelectric tensor
  subroutine vh5_write_piezoelectric_tensor(fileid, prefix, piezoelectric_tensor, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    character(len=*), intent(in) :: prefix
    real(q),intent(in) :: piezoelectric_tensor(:,:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, prefix // "_piezoelectric_tensor", piezoelectric_tensor))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
  end subroutine vh5_write_piezoelectric_tensor

  !> @brief write dipole moment
  subroutine vh5_write_dipole_moment(fileid, prefix, dipole_moment, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    character(len=*), intent(in) :: prefix
    real(q),intent(in) :: dipole_moment(:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, prefix // "_dipole_moment", dipole_moment))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
  end subroutine vh5_write_dipole_moment

  !> @brief write spin resolved dipole moment
  subroutine vh5_write_spin_resolved_dipole_moment(fileid, dipole_moment, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    real(q),intent(in) :: dipole_moment(:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, "spin_resolved_dipole_moment", dipole_moment))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
  end subroutine vh5_write_spin_resolved_dipole_moment

  !> @brief write dielectric tensor
  subroutine vh5_write_dielectric_static(fileid, dielectric_tag, dielectric_tensor, method, group, subgroup)
    integer(HID_T), intent(in)  :: fileid
    character(LEN=*),intent(in) :: dielectric_tag
    real(q),intent(in) :: dielectric_tensor(:,:)
    character(LEN=*), intent(in), optional :: method
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_DIELECTRIC)
    !-------------------------------------------------------------------------------------
    integer(HID_T):: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, dielectric_tag // "_dielectric_tensor", dielectric_tensor))
    if (present(method)) &
       VH5_CHECK(vh5_write(subgroupid, "method_dielectric_tensor", method))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_dielectric_static

  subroutine vh5_write_dielectric_energies(fileid, energy_step, num_energies, group, subgroup)
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    real(q), intent(in) :: energy_step
    integer, intent(in) :: num_energies
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    real(q) :: energies(num_energies)
    integer :: ii
    !-------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    energies = [(energy_step * (ii - 1), ii = 1, num_energies)]
    VH5_CHECK(vh5_write(subgroupid, "energies_dielectric_function", energies))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_dielectric_energies

  !> @brief write energy dependent dielectric function
  subroutine vh5_write_dielectric_dynamic(fileid, dielectric_function, prefix, group, subgroup)
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    complex(q), intent(in) :: dielectric_function(:,:,:)
    character(len=*), intent(in) :: prefix
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_LINEAR_RESPONSE)
    !-------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------
    my_group = GRP_RESULTS;            if (present(group))    my_group = group
    my_subgroup = GRP_LINEAR_RESPONSE; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    VH5_CHECK(vh5_write(subgroupid, prefix // "_dielectric_function", dielectric_function))

    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_dielectric_dynamic

  !> @brief write energy dependent conductivity
  subroutine vh5_write_conductivity_dynamic(fileid, nedos, edos, conductivity, efermi, group, subgroup)
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    real(q),intent(in) :: efermi
    real(q),intent(in) :: edos(nedos), conductivity(nedos,3,3)
    integer,intent(in) :: nedos
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_CONDUCTIVITY)
    !-------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    real(q) :: energies(nedos)
    integer :: i
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------
    my_group = GRP_RESULTS;         if (present(group))    my_group = group
    my_subgroup = GRP_CONDUCTIVITY; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, my_group, groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, my_subgroup, subgroupid))

    ! write stuff
    VH5_CHECK(vh5_write(subgroupid, "energies", edos))
    VH5_CHECK(vh5_write(subgroupid, "conductivity_dynamic", conductivity))
    VH5_CHECK(vh5_write(subgroupid, "efermi", efermi))

    ! close/deallocate
    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
  end subroutine vh5_write_conductivity_dynamic

  !> @brief write contents of waveder file
  subroutine vh5_write_waveder(fileid, cder_between_states, kpoints, group, subgroup)
    integer(HID_T), intent(in)    :: fileid !< hdf5 file handle
    GDEF,intent(in) :: cder_between_states(:,:,:,:,:)
    real(q),intent(in) :: kpoints(:,:)
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_CONDUCTIVITY)
    !-------------------------------------------------------------
    real(qs) :: test
    integer(HID_T) :: groupid, subgroupid
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------
    my_group = GRP_RESULTS;         if (present(group))    my_group = group
    my_subgroup = GRP_CONDUCTIVITY; if (present(subgroup)) my_subgroup = subgroup
    VH5_CHECK(vh5_group_open_or_create(fileid, my_group, groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, my_subgroup, subgroupid))

    ! write stuff
    VH5_CHECK(vh5_write(subgroupid, "cder_between_states", cder_between_states))
    VH5_CHECK(vh5_write(subgroupid, "kpoint_coords", kpoints))

    ! close/deallocate
    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
  end subroutine

  !> @brief write charge (single spin or all spins)
  !>
  !> @param fileid hdf5 file handle
  !> @param subgroup subgroup to store the position data in
  !> @param grid grid_3d grid description type (from mgrid.F)
  !> @param chargetot charge array. shape is (ispin, nodes)
  subroutine vh5_write_charge(fileid, subgroup, grid, chargetot)
    USE prec
    USE string
    USE mpimy
    use mgrid, only : grid_3d
    use radial, only : rgrid
    use tutor, only: vtutor
    integer(HID_T), intent(in)       :: fileid
    character(LEN=*), intent(in)    :: subgroup
    type(grid_3d) ::  grid
    complex(q) :: chargetot(:,:)
    !-------------------------------------------------------------------------------------
    integer(HID_T) :: groupid
    integer(HSIZE_T) :: start(4), count(4), dimensions(4)
    type(rgrid), allocatable :: cwork(:)
    real(q), allocatable :: work(:)
    integer :: nalloc, istat, nx, ny, nz, nspins, is
    ! integer :: NODE_ME, IONODE
    !-------------------------------------------------------------------------------------

    call vh5_write_version(fileid)

    ! check for consistency of grid data
    ! taken from OUTCHG (file fileio.F around line 1051)
    if(grid%nplwv /= grid%ngx * grid%ngy * grid%ngz) then
       call vtutor%bug('internal ERROR: vh5_write_charge NPLWV = ' // str(grid%nplwv) // &
         'is not compatible with NGX = ' // str(grid%ngx) // ', NGY = ' // str(grid%ngy) &
         // ', NGZ = ' // str(grid%ngz), __FILE__, __LINE__)
    end if
    nspins = size(chargetot,2)

    !
    ! write grid data
    !
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), groupid))
    VH5_CHECK(vh5_write(groupid, "grid", [grid%ngx, grid%ngy, grid%ngz]))
    !
    ! allocate work arrays
    !
    nalloc = grid%ngx*grid%ngy
    allocate(cwork(DIMREAL(grid%mplwv)), work(nalloc), stat=istat)
    if (istat > 0) then
       write (*,*) "WARNING: insufficient memory to write charges to hdf5 file"
       VH5_CHECK(vh5_group_close_writing(groupid))
       return
    end if

    !
    ! write charge densities
    !
    dimensions = [grid%ngx, grid%ngy, grid%ngz, nspins]
    VH5_CHECK(vh5_create_double_array_nd(groupid, "charge", 4, dimensions))
    do is = 1, nspins
       !
       ! FFT
       !
       call rc_add(chargetot(1,is), 1.0_q, chargetot(1,is), 0.0_q, cwork, grid)
       call FFT3D(cwork, grid, 1)

       !
       ! prepare charge output
       !
       do nz = 1, grid%ngz
          call mrg_grid_rl_plane(grid, work, cwork, nz)
          start = [1,1,nz,is]
          count = [grid%ngx,grid%ngy,1,1]
          VH5_CHECK(vh5_write_double_subarray_nd(groupid, "charge", 4, start, count, work))
       end do
    end do
    !
    ! deallocate, close groups
    !
    deallocate(cwork, work)
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_charge
  !*******************************************************************
  !
  ! WRT_RHO_PAW write the PAW occupancies to a file specified by
  ! IU
  !
  !*******************************************************************
  subroutine vh5_write_paw_occupancies(fileid, subgroup, p, t_info, loverl, rholm_store, comm, ispin)
    use pseudo
    use poscar
    use wave
    use constant
    use paw

    integer(HID_T), intent(in)       :: fileid
    character(LEN=*), intent(in)    :: subgroup
    type (type_info)   t_info
    type (potcar),target::  p(t_info%ntyp)
    logical  loverl          ! overlap matrix used ?
    real(q) rholm_store(:)   ! storage for the channel occupancies
    type(communic) :: comm
    integer :: ispin

  ! local variables
    integer(HID_T) :: unitgroupid
    type (potcar),pointer :: pp
    integer nt, ni, i
    integer ibase, iadd, lymax, lmmax
    integer, allocatable::  nelements(:)
    integer, external :: maxl_aug
!    integer node_me, ionode
    character(LEN=40) :: datasetname

    real(q), allocatable::  buffer(:,:)

!     node_me=0
!     ionode=0
! #ifdef MPI
!     node_me=comm%node_me
!     ionode =comm%ionode
! #endif

    allocate(nelements(t_info%nions))

    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), unitgroupid))
!=======================================================================
! quick return if possible
!=======================================================================
    if (.not.loverl .or. mimic_us ) return

    lymax =maxl_aug(t_info%ntyp,p)
    lmmax=(lymax+1)**2

    allocate( buffer(lmmax*lmmax, t_info%nions))

!=======================================================================
! cycle all ions and write the required elements
!=======================================================================
    ibase=1

    ion: do ni=1,t_info%nions
       buffer(:,ni)=0
       iadd  =0

       nt=t_info%ityp(ni)
       pp=> p(nt)

       nelements(ni) = 0
       if (do_local(ni)) then
          call retrieve_rholm( buffer(:,ni), rholm_store(ibase:), &
               metric(ibase:), iadd, pp, .true., nelements(ni))

          if (nelements(ni) > lmmax*lmmax) then
             call vtutor%bug('internal error: wrt_rho_paw running out of buffer', __FILE__, __LINE__)
          endif
          ibase=ibase+iadd
       endif

       CALLMPI( m_sum_i(comm, nelements(ni), 1))
       CALLMPI( m_sum_d(comm, buffer(:,ni), nelements(ni)))

      !  io_begin
      !  write(iu,'("augmentation occupancies",2i4)') ni, nelements
      !  write(iu,'(5e15.7)') (buffer(i),i=1,nelements)
      !  io_end

    enddo ion
    write(datasetname,'(A,I0.2)') "aug_occupancies_ions_s",ispin
    VH5_CHECK(vh5_write(unitgroupid, datasetname, nelements))
    write(datasetname,'(A,I0.2)') "aug_occupancies_s",ispin
    VH5_CHECK(vh5_write(unitgroupid, datasetname, buffer))

#ifdef rholm_complex
    ibase=1

    again: do ni=1,t_info%nions
       buffer(:,ni)=0
       iadd  =0

       nt=t_info%ityp(ni)
       pp=> p(nt)

       nelements =0
       if (do_local(ni)) then
          call retrieve_rholm( buffer(:,ni), rholm_store(ibase+ubound(rholm_store,1)/2:), &
               metric(ibase+ubound(rholm_store,1)/2:), iadd, pp, .true., nelements(ni))

          if (nelements(ni) > lmmax*lmmax) then
             call vtutor%bug('internal error: wrt_rho_paw running out of buffer', __FILE__, __LINE__)
          endif

          ibase=ibase+iadd
       endif

       CALLMPI( m_sum_i(comm, nelements(ni), 1))
       CALLMPI( m_sum_d(comm, buffer(:,ni), nelements(ni)))

      !  io_begin
      !  write(iu,'("augmentation occupancies (imaginary part)",2i4)') ni, nelements
      !  write(iu,'(5e15.7)') (buffer(i),i=1,nelements)
      !  io_end

    enddo again
    write(datasetname,'(A,I0.2)') "imag_aug_occupancies_ions_s",ispin
    VH5_CHECK(vh5_write(unitgroupid, datasetname, nelements))
    write(datasetname,'(A,I0.2)') "imag_aug_occupancies_ions_s",ispin
    VH5_CHECK(vh5_write(unitgroupid, datasetname, buffer))

#endif
    deallocate(nelements)
    deallocate(buffer)
    VH5_CHECK(vh5_group_close_writing(unitgroupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_paw_occupancies

  subroutine vh5_write_atommom(fileid, subgroup, t_info)
    integer(HID_T), intent(in)       :: fileid
    character(LEN=*), intent(in)    :: subgroup
    type (type_info)   t_info

    integer(HID_T) :: unitgroupid

    VH5_CHECK(vh5_group_open_or_create(fileid, trim(subgroup), unitgroupid))
    VH5_CHECK(vh5_write(unitgroupid, "atomic_moments", t_info%atomom))
    VH5_CHECK(vh5_group_close_writing(unitgroupid))
    VH5_CHECK(vh5_file_flush(fileid))

  end subroutine vh5_write_atommom


  !> @brief write wavefunctions
  !>
  !> @param fileid hdf5 file handle
  !> @param io io handle
  !> @param wave_des wave function descriptor (see wave.f)
  !> @param wave_spin wave functions (see wave.f)
  !> @param latt_ini lattice definition
  !> @param efermi fermi energy
  subroutine vh5_write_wavefunctions(fileid, io, wave_des, wave_spin, latt_ini, efermi)
    USE base
    USE string
    USE wave
    USE lattice
    USE main_mpi

    integer(HID_T), intent(in)       :: fileid
    type(in_struct) :: io
    type(wavedes)   :: wave_des
    type(wavespin)  :: wave_spin
    type(latt)      :: latt_ini
    real (q)        :: efermi
    !----------------------------------------------------------------------
    ! local variables
    integer(HID_T) :: groupid, spin_groupid, k_groupid
    integer :: ierr, i, j, k, npl, irec, ireclw_old, isp
    integer(HSIZE_T) :: dimensions(2), start(2), count(2)
!    integer :: node_me, ionode
    integer npl_tot, iu0
    real(q) :: rdum, rispin, rnb_tot, rnkpts, rnpl
    logical, pointer :: lconjg(:)=>null()
    integer, pointer :: index(:)=>null()
    character(LEN=256) :: spingroupname, kgroupname
    !----------------------------------------------------------------------
    ! local work arrays
    type(wavedes1) wdes1
    ! FIXME CW is replaced by CPTFWP via macro in symbol.inc. still need to find out
    ! why, and what the exact effects are
    complex(q), allocatable :: CW(:), eig(:)
    !----------------------------------------------------------------------
    !
!     node_me=0
!     ionode=0
! #ifdef MPI
!     node_me=wave_des%comm%node_me
!     ionode =wave_des%comm%ionode
! #endif
    call vh5_write_version(fileid)
    !
    ! report to  stdout
    !

    if (io%iu6>=0) write(io%iu6,*)'writing wavefunctions to vaspwave.h5'
    if (io%iu0>=0) write(io%iu0,*)'writing wavefunctions to vaspwave.h5'
    !
    ! set up constants
    !
    npl_tot = maxval(wave_des%nplwkp_tot)
#ifdef WAVECAR_double
    io%ireclw=max(npl_tot,7)*io%icmplx
    ireclw_old=max(max(npl_tot,6),((wave_des%nb_tot*3+1)/2+2))*io%icmplx
#else
    io%ireclw=max((npl_tot+1)/2,7)*io%icmplx
    ireclw_old=max(max((npl_tot+1)/2,6),((wave_des%nb_tot*3+1)/2+2))*io%icmplx
#endif
    !
    ! allocate work arrays
    !
    allocate(CW(npl_tot),eig(wave_des%nb_tot))
    rispin=wave_des%ispin
    rdum  =io%ireclw
    !
    ! write constants
    !
    ! io_begin
    VH5_CHECK(vh5_group_open_or_create(fileid, GRP_WAVE, groupid))

    VH5_CHECK(vh5_write(groupid, "rdum", rdum))
    VH5_CHECK(vh5_write(groupid, "rispin", rispin))
    ! io_end

    irec=2
    ! io_begin
    ! in order to increase exchangeability of WAVECAR files across IEEE platforms
    ! avoid INTEGERS on output, write REAL(q) items instead (same below with RNPL)
    rnkpts =wave_des%nkpts
    rnb_tot=wave_des%nb_tot
    !
    ! write some more parameters
    !
    VH5_CHECK(vh5_write(groupid, "rnkpts", rnkpts))
    VH5_CHECK(vh5_write(groupid, "rnb_tot", rnb_tot))
    VH5_CHECK(vh5_write(groupid, "efermi", efermi))
    VH5_CHECK(vh5_write(groupid, "enmax", wave_des%enmax))
    VH5_CHECK(vh5_write(groupid, "amat", latt_ini%a))
    ! io_end
    !
    ! loop over spin, kpoints, bands
    !
    spin: DO ISP=1,wave_des%ISPIN
       write(spingroupname, '(2A)') "spin_"//trim(str(isp))
       VH5_CHECK(vh5_group_open_or_create(groupid, trim(spingroupname), spin_groupid))
       kpoints: DO K=1,wave_des%NKPTS
          write(kgroupname, '(2A)') "kpoint_"//trim(str(k))
          VH5_CHECK(vh5_group_open_or_create(spin_groupid, trim(kgroupname), k_groupid))
#ifdef MPI
          IF (wave_des%COMM_KINTER%NCPU.GT.1) THEN
             IF (wave_des%COMM_KINTER%NODE_ME.NE.1) THEN
                IF (MOD(K-1,wave_des%COMM_KINTER%NCPU).NE.wave_des%COMM_KINTER%NODE_ME-1) THEN
                   CYCLE
                ELSE
                   CALLMPI( M_send_z(wave_des%COMM_KINTER, 1, wave_spin%CW(1,1,K,ISP), SIZE(wave_spin%CW,1)*SIZE(wave_spin%CW,2)) )
                END IF
             ELSE
                IF (MOD(K-1,wave_des%COMM_KINTER%NCPU).NE.wave_des%COMM_KINTER%NODE_ME-1) THEN
                   CALLMPI( M_recv_z(wave_des%COMM_KINTER, MOD(K-1,wave_des%COMM_KINTER%NCPU)+1, wave_spin%CW(1,1,K,ISP), SIZE(wave_spin%CW,1)*SIZE(wave_spin%CW,2)) )
                ENDIF
             ENDIF
          END IF
#endif
          CALL SETWDES(wave_des,WDES1,K)
#ifdef gammareal
          CALL GEN_INDEX_GAMMA(WDES1, latt_ini%B, wave_des%ENMAX, INDEX, LCONJG)
#endif
          npl=wave_des%nplwkp_tot(k)
          rnpl=npl
          ! write number of plane waves, k-point coordinates and all eigenvalues and
          ! occupation numbers for current k-point K
          ! write eigenvalues in real format
          ! io_begin
          VH5_CHECK(vh5_write(k_groupid, "num_planewaves", npl))
          VH5_CHECK(vh5_write(k_groupid, "vkpt", wave_des%vkpt(:,k)))
          VH5_CHECK(vh5_write(k_groupid, "celtot", wave_spin%celtot(:,k,isp)))
          VH5_CHECK(vh5_write(k_groupid, "fertot", wave_spin%fertot(:,k,isp)))
          ! io_end

          dimensions = [npl, wave_des%nb_tot]
#ifdef WAVECAR_double
          VH5_CHECK(vh5_create_double_complex_array_nd(k_groupid, "wave", 2, dimensions))
#else
          VH5_CHECK(vh5_create_complex_array_nd(k_groupid, "wave", 2, dimensions))
#endif

          DO J=1,wave_des%NB_TOT
#ifdef gammareal
             CALL MRG_PW_BAND_GAMMA(WDES1, J, CW, wave_spin%CW(1,1,K,ISP), INDEX(1), LCONJG(1))
#else
             CALL MRG_PW_BAND(WDES1, J, CW, wave_spin%CW(1,1,K,ISP))
#endif
             irec=irec+1

             start = [1,J]
             count = [npl,1]
             VH5_CHECK(vh5_write_double_complex_subarray_nd(k_groupid, "wave", 2, start, count, CW))
          enddo

#ifdef gammareal
          call free_index_gamma(index, lconjg)
#endif
          VH5_CHECK(vh5_group_close_writing(k_groupid))
       enddo kpoints
       VH5_CHECK(vh5_group_close_writing(spin_groupid))
    enddo spin
    ! io_begin
    ! FIXME: this is an example for demonstrating the possibilities of the file format
    ! the same information is also writen per k-point and spin in the above loop. unlike
    ! CW, these arrays do not depend on the kpoint (while the size of CW does), so
    ! we could also write them as arrays, not as subgroups. in the end, we'll need to
    ! decide on one of the possibilities
    VH5_CHECK(vh5_write(groupid, "nplwkp_tot", wave_des%NPLWKP_TOT))
    VH5_CHECK(vh5_write(groupid, "vkpt", wave_des%VKPT))
    VH5_CHECK(vh5_write(groupid, "celtot", wave_spin%CELTOT))
    VH5_CHECK(vh5_write(groupid, "fertot", wave_spin%FERTOT))
    ! io_end
    !
    ! deallocate/close
    !
    deallocate(CW,eig)
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_wavefunctions

  !> read force constants from hdf5 file
  subroutine vh5_read_force_constants(fileid,force_constants)
    use tutor, only: vtutor
    use string, only: str
    integer(HID_T),intent(in) :: fileid
    real(q),allocatable,intent(inout) :: force_constants(:,:)
    ! local variables
    integer :: dims(2)
    real(q) :: fc(1,1)
    integer(HID_T) :: ih5inresultsgroup_id
    integer :: ierr

    ierr = vh5_read(fileid, '/results/linear_response/force_constants', force_constants, dimensions=dims)
    if (ierr/=0) then
        call vtutor%error("Could not find '/results/linear_response/force_constants' in vaspin.h5")
    endif
    ! check the dimensions
    if (any(shape(force_constants)/=dims)) then
        call vtutor%error("The dimensions of the force_constants array in the HDF5 file "//str(dims(1)/3)//&
                          " are incompatible with the current VASP run "//str(size(force_constants,1)))
    endif
    VH5_CHECK(vh5_read(fileid, '/results/linear_response/force_constants', force_constants))
  end subroutine vh5_read_force_constants

  !> @brief read wave function header
  !>
  !> @param fileid hdf5 file handle
  !> @param wave_des wave function descriptor (see wave.f)
  !> @param wave_spin wave functions (see wave.f)
  !> @param latt_ini lattice definition
  !> @param latt_cur lattice definition
  !> @param enmaxi
  !> @param istart
  !> @param iu0
  function vh5_read_wavefunction_header(fileid, wave_des, latt_ini, latt_cur, enmaxi, istart, iu0 ) result(ierr)
    USE base
    USE wave
    USE lattice
    USE main_mpi

    integer(HID_T), intent(in)       :: fileid
    type(wavedes)   :: wave_des
    type(wavespin)  :: wave_spin
    type(latt)      :: latt_ini
    type(latt)      :: latt_cur
    real (q)        :: enmaxi
    integer         :: istart
    integer         :: iu0
    !-------------------------------------------------------------
    ! local variables
    integer(HID_T) :: groupid
    integer :: ierr, i, j, k, n, ionode, irec, ireclw_old, isp
    integer :: node_me, npl
    integer :: npl_tot, nbandf, nkptsf
    real(q) :: rdum, rispin, rnb_tot, rnkpts, rnpl!, enmax
    logical :: ldiff
    !-------------------------------------------------------------

    node_me=0
    ionode=0
#ifdef MPI
    node_me=wave_des%comm%node_me
    ionode =wave_des%comm%ionode
#endif

    ! open group
    ierr = vh5_group_open(fileid, GRP_WAVE, groupid)
    if(ierr .EQ. 0) then
       if (iu0 >= 0) then
          do_io write(iu0,*) 'found vaspwave.h5, reading the header'
       endif

       VH5_CHECK(vh5_read(groupid, "rnkpts", rnkpts))
       VH5_CHECK(vh5_read(groupid, "rnb_tot", rnb_tot))
       VH5_CHECK(vh5_read(groupid, "enmax", enmaxi))
       VH5_CHECK(vh5_read(groupid, "amat", latt_ini%a))


       NKPTSF=NINT(rnkpts)
       NBANDF=NINT(rnb_tot)

       CALL LATTIC(LATT_INI)

       IF (ISTART==2 .AND. ENMAXI /= wave_des%ENMAX) THEN
          CALL vtutor%error('ERROR: ENMAX changed please set ISTART to 1')
       ENDIF

       IF (.NOT. (NBANDF==wave_des%NB_TOT .OR. (wave_des%NRSPINORS==2 .AND.  NBANDF*wave_des%NRSPINORS==wave_des%NB_TOT))) THEN
          IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
               'number of bands has changed, file:',NBANDF,' present:',wave_des%NB_TOT
          IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
               'trying to continue reading WAVECAR, but it might fail'
       ENDIF
       IF (NKPTSF/=wave_des%NKPTS) THEN
          IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
               'number of k-points has changed, file:',NKPTSF,' present:',wave_des%NKPTS
          IF (IU0 >= 0) WRITE(IU0,'(2X,A,I6,A,I6)') &
               'trying to continue reading WAVECAR, but it might fail'
       ENDIF
       IF (ISTART ==1) THEN

          LDIFF=.FALSE.
          DO I=1,3
             DO J=1,3
                IF (ABS(LATT_INI%A(I,J)-LATT_CUR%A(I,J)) > 1E-4) LDIFF=.TRUE.
             ENDDO
          ENDDO
          IF (ENMAXI /= wave_des%ENMAX) LDIFF=.TRUE.
          IF (LDIFF) THEN
             IF (IU0>=0) &
                  WRITE(IU0,*)'WAVECAR: different cutoff or change in lattice found'
          ENDIF

       ENDIF

       RETURN
    else
       latt_ini%a = latt_cur%a
       call lattic(latt_ini)
       if (istart==3) then
          if (iu0 >= 0) then
             write(iu0,*) "ERROR: can't restart (ISTART=3) with wavefunctions on file"
          endif
       else
          istart = 0
       endif
    endif

  end function vh5_read_wavefunction_header


  !> @brief read wave functions
  !>
  !> @param fileid hdf5 file handle
  !> @param io io handle
  !> @param wave_des wave function descriptor (see wave.f)
  !> @param wave_spin wave functions (see wave.f)
  !> @param latt_cur lattice definition
  !> @param latt_ini lattice definition
  !> @param istart
  !> @param efermi fermi energy
  function vh5_read_wavefunctions(fileid, io, wave_des, wave_spin, grid, latt_cur, latt_ini, istart, efermi ) result(ierr)
    USE base
    USE string
    USE wave_high
    use mgrid
    USE lattice
    USE main_mpi

    integer(HID_T), intent(in)       :: fileid
    type(in_struct) :: io
    type(wavedes)   :: wave_des
    type(wavespin)  :: wave_spin
    type(grid_3d)   :: grid
    type(latt)      :: latt_ini
    type(latt)      :: latt_cur
    real (q)        :: efermi
    integer         :: istart
    !-------------------------------------------------------------
    ! local variables
    integer(HID_T) :: groupid, spin_groupid, k_groupid
    integer :: ierr, i, j, k, n, ionode, irec, ireclw_old, isp, kpl
    integer :: node_me, npl, nmaxf, ifail, irecl, is, ispinread, iu0
    integer :: npl_tot, nbandf, nkptsf, malloc, ncomm, nplread
    real(q) :: rdum, rispin, rnb_tot, rnkpts, rnpl, enmax, rbandf, rkptsf, rtag, enmaxf
    logical :: ldiff
    logical, pointer :: lconjg(:)=>null()
    integer, pointer :: index(:)=>null()
    character(LEN=256) :: spingroupname, kgroupname
    !----------------------------------------------------------------------
    ! local work arrays
    type(wavedes1) wave_des1
    ! FIXME CW is replaced by CPTFWP via macro in symbol.inc. still need to find out
    ! why, and what the exact effects are
    complex(q), allocatable :: CW(:), eig(:)
    REAL(q) VKPT(3)
    INTEGER,ALLOCATABLE    :: IND(:),INDI(:)
    COMPLEX(q),ALLOCATABLE :: CW1(:),CW2(:)
    CHARACTER(:),ALLOCATABLE :: MSG
    LOGICAL ::  SINGLE_PREC, SPAN_RECORDS
    REAL(q), ALLOCATABLE :: INBUF(:)
    integer(HSIZE_T) :: start(2), count(2)
    !----------------------------------------------------------------------

    iu0 = io%iu0

    node_me=0
    ionode=0
#ifdef MPI
    node_me=wave_des%comm%node_me
    ionode =wave_des%comm%ionode
#endif

    ierr = vh5_group_open(fileid, GRP_WAVE, groupid)
    if(ierr .EQ. 0) then
       if (iu0 >= 0) then
          do_io write(iu0,*) 'found vaspwave.h5, reading file, group = ', GRP_WAVE
       endif

       RTAG=0
       VH5_CHECK(vh5_read_double_scalar(groupid, "rispin", rispin))

       ISPINREAD=NINT(RISPIN)
       VH5_CHECK(vh5_read(groupid, "rnkpts", rkptsf))
       VH5_CHECK(vh5_read(groupid, "rnb_tot", rbandf))
       VH5_CHECK(vh5_read(groupid, "enmax", enmaxf))
       VH5_CHECK(vh5_read(groupid, "amat", latt_ini%a))
       VH5_CHECK(vh5_read(groupid, "efermi", efermi))

       IREC=2
       NKPTSF=NINT(RKPTSF)
       NBANDF=NINT(RBANDF)

       CALL LATTIC(LATT_INI)
       !=======================================================================
       ! read WAVECAR file, number of bands agree
       !=======================================================================
       IF (.NOT.(wave_des%NRSPINORS==2 .AND. NBANDF*wave_des%NRSPINORS == wave_des%NB_TOT)) THEN

          spin:    DO ISP=1, MIN(wave_des%ISPIN, ISPINREAD)
             write(spingroupname, '(2A)') "spin_"//trim(str(isp))
             VH5_CHECK(vh5_group_open(groupid, trim(spingroupname), spin_groupid))
             kpoints: DO K=1,wave_des%NKPTS
                write(kgroupname, '(2A)') "kpoint_"//trim(str(k))
                VH5_CHECK(vh5_group_open(spin_groupid, trim(kgroupname), k_groupid))
                IF ( K> NKPTSF ) THEN
                   IREC=IRECL
                ELSE
                   IRECL=IREC
                ENDIF

                CALL SETWDES(wave_des,wave_des1,K)
#ifdef gammareal
                CALL GEN_INDEX_GAMMA(wave_des1, LATT_CUR%B, wave_des%ENMAX, INDEX, LCONJG)
#endif

                io_begin
                ifail = 0
                VH5_CHECK(vh5_read(k_groupid, "num_planewaves", nplread))
                VH5_CHECK(vh5_read(k_groupid, "vkpt", vkpt))
                VH5_CHECK(vh5_read(k_groupid, "celtot", wave_spin%CELTOT(:,K,ISP)))
                VH5_CHECK(vh5_read(k_groupid, "fertot", wave_spin%FERTOT(:,K,ISP)))
                io_end

                CALLMPI( M_bcast_i( wave_des%COMM, IFAIL, 1))
                IF (IFAIL /=0) GOTO 230

                CALLMPI( M_bcast_i( wave_des%COMM, IREC, 1))
                CALLMPI( M_bcast_i( wave_des%COMM, NPLREAD, 1))
                NPL=wave_des%NPLWKP_TOT(K)

                MALLOC=MAX(NPL, NPLREAD)
                ALLOCATE(CW1(MALLOC),CW2(MALLOC),IND(MALLOC),INDI(MALLOC))

                io_begin
                ! create an index to allow for change of cutoff or cell size
                IF (ISTART==2) THEN
                   ! gK: fix ISTART=2 in INWAV_FAST
                   CALL REPAD_INDEX_ARRAY(GRID, wave_des%VKPT(:,K), VKPT, LATT_INI%B,  LATT_INI%B, &
                        wave_des%ENMAX, ENMAXF, NPL/wave_des%NRSPINORS, NPLREAD/wave_des%NRSPINORS, IND, INDI, MALLOC, IFAIL )
                ELSE
                   CALL REPAD_INDEX_ARRAY(GRID, wave_des%VKPT(:,K), VKPT, LATT_CUR%B,  LATT_INI%B, &
                        wave_des%ENMAX, ENMAXF, NPL/wave_des%NRSPINORS, NPLREAD/wave_des%NRSPINORS, IND, INDI, MALLOC, IFAIL )
                ENDIF
                io_end

                CALLMPI( M_bcast_i( wave_des%COMM, IFAIL, 1))
                IF (IFAIL /=0) GOTO 220

                band: DO J=1,NBANDF
                   IREC=IREC+1
                   IF (J>wave_des%NB_TOT) CYCLE
                   start = [1,J]
                   count = [nplread,1]

                   io_begin
                   VH5_CHECK(vh5_read_double_complex_subarray_nd(k_groupid, "wave", 2, start, count, cw2))

                   CW1=0
                   DO IS=1,wave_des%NRSPINORS
                      ! store the wave function coefficients according to new cutoff
                      CALL REPAD_WITH_INDEX_ARRAY( MALLOC, IND, INDI, &
                           CW1((IS-1)*NPL/wave_des%NRSPINORS+1), CW2((IS-1)*NPLREAD/wave_des%NRSPINORS+1))
                   ENDDO

                   io_end
#ifdef MPI
                   IF (wave_des%COMM_KINTER%NCPU.GT.1) THEN
                      CALLMPI( M_bcast_z(wave_des%COMM_KINTER,CW1,SIZE(CW1)) )
                   END IF
#endif
#ifdef gammareal
                   ! gamma only go from seriel layout to parallel layout if required
                   CALL DIS_PW_BAND_GAMMA(wave_des1, J, CW1, wave_spin%CW(1,1,K,ISP), INDEX(1), LCONJG(1))
#else
                   CALL DIS_PW_BAND(wave_des1, J, CW1, wave_spin%CW(1,1,K,ISP))
#endif
                ENDDO band

#ifdef gammareal
                CALL FREE_INDEX_GAMMA(INDEX, LCONJG)
#endif
                DEALLOCATE(CW1,CW2,IND,INDI)

                VH5_CHECK(vh5_group_close(k_groupid))
             ENDDO kpoints


             IF (NKPTSF > wave_des%NKPTS) THEN
                IREC=IREC+(NKPTSF-wave_des%NKPTS)*(wave_des%NB_TOT+1)
             ENDIF

             ! copy eigenvalues and weights to all nodes
             NCOMM=wave_des%NB_TOT*wave_des%NKPTS
             CALLMPI( M_bcast_z(wave_des%COMM, wave_spin%CELTOT(1,1,ISP),NCOMM ))
             CALLMPI( M_bcast_d(wave_des%COMM, wave_spin%FERTOT(1,1,ISP),NCOMM ))

             VH5_CHECK(vh5_group_close(spin_groupid))
          ENDDO spin

          IF (ISPINREAD > wave_des%ISPIN .AND. IU0>=0 ) THEN
             WRITE(IU0,*) 'down-spin wavefunctions not read'
          ENDIF

          IF (NBANDF<wave_des%NB_TOT) THEN
             IF (IU0>=0) WRITE(IU0,*) 'random initialization beyond band ',NBANDF
             CALL WFINIT(wave_des, wave_spin, 1E10_q, NBANDF+1) ! ENINI=1E10 not cutoff restriction
          ENDIF

          IF (IU0 >= 0) WRITE(IU0,*) 'the vaspwave.h5 file was read successfully'

          IF (wave_des%ISPIN<=ISPINREAD) THEN
             RETURN
          ENDIF
          !
          !  spin down is missing
          !
          IF (IU0>=0) &
               WRITE(IU0,*) 'No down-spin wavefunctions found', &
               &             ' --> setting down-spin equal up-spin ...'

          DO K=1,wave_des%NKPTS
             wave_spin%CELTOT(1:wave_des%NB_TOT,K,2)=wave_spin%CELTOT(1:wave_des%NB_TOT,K,1)
             wave_spin%FERTOT(1:wave_des%NB_TOT,K,2)=wave_spin%FERTOT(1:wave_des%NB_TOT,K,1)
             NPL=wave_des%NPLWKP(K)
             wave_spin%CW(1:NPL,1:wave_des%NBANDS,K,2)=wave_spin%CW(1:NPL,1:wave_des%NBANDS,K,1)
          ENDDO


          RETURN
          !=======================================================================
          ! read collinear WAVECAR file for a non collinear run
          !=======================================================================
       ELSE
          wave_spin%CW=0

          spin2:    DO ISP=1,ISPINREAD
             IF (IU0>=0.AND. ISP==1) &
                  WRITE(IU0,*) 'reading wavefunctions of collinear run, up'
             IF (IU0>=0.AND. ISP==2) &
                  WRITE(IU0,*) 'reading wavefunctions of collinear run, down'
             write(spingroupname, '(A,I0.2)') "spin_", isp
             VH5_CHECK(vh5_group_open(groupid, trim(spingroupname), spin_groupid))
             kpoints2: DO K=1,NKPTSF

                write(kgroupname, '(A,I0.4)') "kpoint_", k
                VH5_CHECK(vh5_group_open(spin_groupid, trim(kgroupname), k_groupid))

                IF ( K> NKPTSF ) THEN
                   IREC=IRECL
                ELSE
                   IRECL=IREC
                ENDIF

                CALL SETWDES(wave_des,wave_des1,K)


                io_begin
                ifail = 0
                VH5_CHECK(vh5_read_integer_scalar(k_groupid, "num_planewaves", nplread))
                VH5_CHECK(vh5_read_double_array_1d(k_groupid, "vkpt", vkpt))
                VH5_CHECK(vh5_read_double_complex_array_1d(k_groupid, "celtot", wave_spin%CELTOT(:,K,ISP)))
                VH5_CHECK(vh5_read_double_array_1d(k_groupid, "fertot", wave_spin%FERTOT(:,K,ISP)))
                io_end

                CALLMPI( M_bcast_i( wave_des%COMM, IFAIL, 1))
                IF (IFAIL /=0) GOTO 230

                CALLMPI( M_bcast_i( wave_des%COMM, IREC, 1))
                CALLMPI( M_bcast_i( wave_des%COMM, NPLREAD, 1))
                NPL=wave_des%NPLWKP_TOT(K)/2

                MALLOC=MAX(NPL, NPLREAD)
                ALLOCATE(CW1(2*MALLOC),CW2(2*MALLOC),IND(MALLOC),INDI(MALLOC))

                io_begin
                IF (ISTART==2) THEN
                   ! gK: fix ISTART=2 in INWAV_FAST
                   CALL REPAD_INDEX_ARRAY(GRID, wave_des%VKPT(:,K), VKPT, LATT_INI%B,  LATT_INI%B, &
                        wave_des%ENMAX, ENMAXF, NPL, NPLREAD, IND, INDI, MALLOC, IFAIL )
                ELSE
                   CALL REPAD_INDEX_ARRAY(GRID, wave_des%VKPT(:,K), VKPT, LATT_CUR%B,  LATT_INI%B, &
                        wave_des%ENMAX, ENMAXF, NPL, NPLREAD, IND, INDI, MALLOC, IFAIL )
                ENDIF ! istart == 2
                io_end

                CALLMPI( M_bcast_i( wave_des%COMM, IFAIL, 1))
                IF (IFAIL /=0) GOTO 220

                IF (ISP==1) THEN
                   DO J=1,NBANDF
                      IREC=IREC+1
                      start = [1,J]
                      count = [nplread,1]

                      io_begin
                      VH5_CHECK(vh5_read_double_complex_subarray_nd(k_groupid, "wave", 2, start, count, cw2))

                      CW1=0
                      CALL REPAD_WITH_INDEX_ARRAY( MALLOC, IND, INDI, CW1, CW2)
                      io_end

#ifdef MPI
                      IF (wave_des%COMM_KINTER%NCPU.GT.1) THEN
                         CALLMPI( M_bcast_z(wave_des%COMM_KINTER,CW1,SIZE(CW1)) )
                      END IF ! wave_des
#endif
                      CALL DIS_PW_BAND(wave_des1, J, CW1, wave_spin%CW(1,1,K,1))

                      ! copy immediately to second panel
                      wave_spin%CELTOT(J+NBANDF,K,1)=wave_spin%CELTOT(J,K,1)
                      wave_spin%FERTOT(J+NBANDF,K,1)=wave_spin%FERTOT(J,K,1)

                      CW1(NPL+1:2*NPL)=CW1(1:NPL)
                      CW1(1:NPL)=0

                      CALL DIS_PW_BAND(wave_des1, J+NBANDF, CW1, wave_spin%CW(1,1,K,1))
                   ENDDO
                ELSE ! ISP==1
                   DO J=1,NBANDF
                      IREC=IREC+1
                      start = [1,J]
                      count = [nplread,1]

                      io_begin
                      VH5_CHECK(vh5_read_double_complex_subarray_nd(k_groupid, "wave", 2, start, count, cw2))

                      CW1=0
                      CALL REPAD_WITH_INDEX_ARRAY( MALLOC, IND, INDI, CW1(NPL+1), CW2)
                      io_end
#ifdef MPI
                      IF (wave_des%COMM_KINTER%NCPU.GT.1) THEN
                         CALLMPI( M_bcast_z(wave_des%COMM_KINTER,CW1,SIZE(CW1)) )
                      END IF
#endif
                      CALL DIS_PW_BAND(wave_des1, J+NBANDF, CW1, wave_spin%CW(1,1,K,1))
                   ENDDO

                ENDIF  ! ISP == 1

                DEALLOCATE(CW1,CW2,IND,INDI)

                VH5_CHECK(vh5_group_close(k_groupid))
             ENDDO kpoints2

             IF (NKPTSF > wave_des%NKPTS) THEN
                IREC=IREC+(NKPTSF-wave_des%NKPTS)*(NBANDF+1)
             ENDIF

             VH5_CHECK(vh5_group_close(spin_groupid))
          ENDDO spin2

          ! copy eigenvalues and weights to all nodes
          NCOMM=wave_des%NB_TOT*wave_des%NKPTS
          CALLMPI( M_bcast_z(wave_des%COMM, wave_spin%CELTOT(1,1,1),NCOMM ))
          CALLMPI( M_bcast_d(wave_des%COMM, wave_spin%FERTOT(1,1,1),NCOMM ))

          IF (IU0 >= 0) WRITE(IU0,*) 'the vaspwave.h5 file was read successfully'
          RETURN
       ENDIF ! num bands agree

       MSG = 'ERROR: while reading vaspwave.h5, file is incompatible'
       IF (wave_des%ENMAX /= ENMAXF) &
          MSG = MSG // ' the energy cutoff has changed (new,old) ' // &
             str(wave_des%ENMAX) // " " // str(ENMAXF)
       IF (NBANDF /= wave_des%NB_TOT) &
          MSG = MSG // ' the number of bands has changed (new,old) ' // &
             str(wave_des%NB_TOT) // " " // str(NBANDF)
       IF (NKPTSF /= wave_des%NKPTS) &
          MSG = MSG // ' the number of k-points has changed (new,old) ' // &
             str(wave_des%NKPTS) // " " // str(NKPTSF)

       CALL vtutor%error(MSG)

    else
       ! could not open the group WAVE_GRP in file, handle error

       CALL vtutor%error('Error reading hdf5 wave function file: group does not exist')

    endif

    ! can not do anything with WAVECAR
    ! hard stop, pull all breaks
    !=======================================================================
200 CONTINUE

    CALL vtutor%error('ERROR: while reading vaspwave.h5, header is corrupt')

220 CONTINUE

    CALL vtutor%error('ERROR: while reading vaspwave.h5, plane wave coefficients changed ' &
       // str(NPL) // ' ' // str(NPLREAD))

230 CONTINUE

    CALL vtutor%error('ERROR: while reading eigenvalues from vaspwave.h5 ' &
       // str(K) // ' ' // str(ISP))

240 CONTINUE

    CALL vtutor%error('ERROR: while reading plane wave coef. from vaspwave.h5 ' &
       // str(K) // ' ' // str(ISP) // ' ' // str(J))

  end function vh5_read_wavefunctions

  !> @brief write projectors
  subroutine vh5_write_projectors(fileid, typeinfo, wave_des, wave_spin, lorbit, &
       lpar, par, lchar, phase, group, subgroup)
    integer(HID_T), intent(in)   :: fileid !< hdf5 file handle
    type(type_info), intent(in)  :: typeinfo
    type(wavedes)   :: wave_des
    type(wavespin)  :: wave_spin
    integer :: lorbit, lpar
    real(q), dimension(wave_des%nb_tot, wave_des%nkpts, lpar, typeinfo%nionp, wave_des%ncdij) :: par
    character(len=*), dimension(:) :: lchar
    complex(q), dimension(:,:,:,:,:) :: phase
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_PROJECTORS)
    !-------------------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    character(len=256) :: lcharout
    integer :: nl
    integer(HSIZE_T), dimension(5) :: par_shape;
    integer(HSIZE_T), dimension(5) :: phase_shape;
    character(LEN=40) :: my_group, my_subgroup
    !-------------------------------------------------------------
    my_group    = GRP_RESULTS;    if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_PROJECTORS; if (present(subgroup)) my_subgroup = trim(subgroup)
    VH5_CHECK(vh5_group_open_or_create(fileid, trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    ! write procar type
    if (lorbit == 11) then
       VH5_CHECK(vh5_write(subgroupid, "procar_type", 'PROCAR lm decomposed'))
    else if (lorbit == 12) then
       VH5_CHECK(vh5_write(subgroupid, "procar_type", 'PROCAR lm decomposed + phase'))
    else
       VH5_CHECK(vh5_write(subgroupid, "procar_type", 'PROCAR new format'))
    end if
    !
    ! write parameters
    ! FIXME these parameters are also written elsewhere. we need to discuss how to structure
    ! a file where these things are only written once. ATM, I'm writing them out because
    ! I'm not sure whether they correspond to the same state of the program as the other
    ! instances
    VH5_CHECK(vh5_write(subgroupid, "nb_tot", wave_des%nb_tot))
    VH5_CHECK(vh5_write(subgroupid, "nionp", typeinfo%nionp))
    VH5_CHECK(vh5_write(subgroupid, "nkpoints", wave_des%nkpts))
    VH5_CHECK(vh5_write(subgroupid, "kpoints", wave_des%vkpt))
    VH5_CHECK(vh5_write(subgroupid, "celtot", wave_spin%celtot))
    !
    ! write lchars
    !
    VH5_CHECK(vh5_write(subgroupid, "lchar", lchar))
    !
    ! write par
    !
    par_shape = shape(par)
    VH5_CHECK(vh5_write_double_array_nd(subgroupid, "par", 5, par_shape, par))
    !
    ! write phase
    !
    if (lorbit == 12) then
       phase_shape = shape(phase)
       VH5_CHECK(vh5_write_double_complex_array_nd(subgroupid, "phase", 5, phase_shape, phase))
    end if
    !
    ! calculate and write parsum and sumion
    !

    !
    ! close/deallocate
    !
    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_projectors

end module vhdf5
#else
!--------------------------------------------------------------------
! Dummy module
!--------------------------------------------------------------------
module vhdf5
contains
  subroutine vasphdf5_dummy
    write(*,*)'Im a DEC compiler so I need this line'
  end subroutine vasphdf5_dummy
end module vhdf5
#endif
