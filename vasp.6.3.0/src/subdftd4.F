module vdwd4

   use prec, only: q
   use base, only: in_struct
#ifdef DFTD4
   use class_set, only: dftd_options
   use class_param, only: dftd_parameter
#endif

   implicit none

   integer, parameter :: UnknownFunctional = 1, UnknownAtom = -1, VdwDFTD4 = 13

#ifdef DFTD4
   logical, parameter :: dftd4Available = .true.
#else
   logical, parameter :: dftd4Available = .false.
   type dftd_options
   end type dftd_options
   type dftd_parameter
   end type dftd_parameter
#endif

   type vdw_settings
      integer iVdw
      type(in_struct) io
      type(dftd_options) options
      type(dftd_parameter) parameters
   end type vdw_settings

   type vdw_structure
      real(q), pointer :: cell(:,:) => null()
      integer, pointer :: types(:) => null()
      character(len=2), pointer :: elements(:) => null()
      real(q), pointer :: positions(:,:) => null()
   end type vdw_structure

   type vdw_interaction
      real(q) energy, stress(3,3)
      real(q), allocatable :: forces(:,:)
   end type vdw_interaction

   interface vdw_structure
      module procedure makeStructure
   end interface vdw_structure

   private
   public vdw_settings, defaultSettings, vdw_forces_D4, UnknownFunctional, &
      vdw_interaction, vdw_structure, vdw_reader, VdwDFTD4, dftd4Available

contains

   type(vdw_interaction) function vdw_forces_D4(settings, structure) result (res)
#ifdef DFTD4
      use class_molecule, only: molecule
      use class_results, only: dftd_results
      use mctc_environment, only: mctc_logger
      use dispersion_calculator, only: d4_calculation
#else
      use tutor, only: vtutor
#endif
      type(vdw_settings), intent(in) :: settings
      type(vdw_structure), intent(in) :: structure
#ifdef DFTD4
      type(molecule) mol
      type(mctc_logger) env
      type(dftd_results) dresults
      mol = convertStructureToMolecule(structure)
      call d4_calculation(settings%io%iu6, env, settings%options, mol, settings%parameters, dresults)
      res = convertUnitsOfResults(dresults, structure)
#else
      call vtutor%error("The code was not compiled with DFTD4 support.")
#endif
   end function vdw_forces_D4

   type(vdw_structure) function makeStructure(latt_cur, dyn, t_info, elem) result (structure)
      use lattice, only: latt
      use poscar, only: dynamics, type_info
      type(latt), target, intent(in) :: latt_cur
      type(type_info), target, intent(in) :: t_info
      type(dynamics), target, intent(in) :: dyn
      character(len=2), target, intent(in) :: elem(:)
      structure%cell => latt_cur%a
      structure%types => t_info%iTyp
      structure%elements => elem
      structure%positions => dyn%posIon
   end function makeStructure

#ifdef DFTD4
   type (molecule) function convertStructureToMolecule(structure) result (mol)
      use constant, only: AuToA
      use class_molecule, only: molecule
      use pbc_tools, only: dlat_to_cell, dlat_to_dvol, dlat_to_rlat
      use tutor, only: vtutor
      type(vdw_structure), intent(in) :: structure
      interface
         subroutine generate_wsc(mol, wsc)
            use class_molecule
            use class_wsc
            type(molecule), intent(inout) :: mol
            type(ws_cell), intent(inout) :: wsc
         end subroutine generate_wsc
      end interface
      call mol%allocate(size(structure%types), .false.)
      mol%at = atomicNumber(structure%elements(structure%types))
      if (any(mol%at == UnknownAtom)) then
         call vtutor%bug("Could not convert the elements to atomic numbers", __FILE__, __LINE__)
      end if
      mol%lattice = structure%cell / AuToA
      mol%xyz = matmul(mol%lattice, structure%positions)
      mol%chrg = 0.0_q
      mol%npbc = 3
      mol%pbc = .true.
      mol%volume = dlat_to_dvol(mol%lattice)
      call dlat_to_cell(mol%lattice, mol%cellpar)
      call dlat_to_rlat(mol%lattice, mol%rec_lat)
      call mol%wrap_back
      call mol%calculate_distances
      call generate_wsc(mol, mol%wsc)
   end function convertStructureToMolecule

   type (vdw_interaction) function convertUnitsOfResults(dresults, structure) result(res)
      use constant, only: RyToEv, AuToA
      use class_results, only: dftd_results
      type(dftd_results), intent(in) :: dresults
      type(vdw_structure), intent(in) :: structure
      real(q), parameter :: HtrToEv = 2 * RyToEv
      allocate(res%forces(3, size(structure%types)))
      res%energy = dresults%energy * HtrToEv
      res%forces = -dresults%gradient * HtrToEv / AuToA
      res%stress = -dresults%lattice_gradient * HtrToEv / AuToA
      res%stress = matmul(res%stress, transpose(structure%cell))
   end function convertUnitsOfResults

   elemental integer function atomicNumber(element)
      use string, only: lowercase
      character(len=2), intent(in) :: element
      character(len=2), parameter :: elements(94) = ['h ', 'he', &
         'li', 'be', 'b ', 'c ', 'n ', 'o ', 'f ', 'ne', &
         'na', 'mg', 'al', 'si', 'p ', 's ', 'cl', 'ar', &
         'k ', 'ca', 'sc', 'ti', 'v ', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', &
            'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr', &
         'rb', 'sr', 'y ', 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', &
            'cd', 'in', 'sn', 'sb', 'te', 'i ', 'xe', &
         'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy', &
            'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w ', 're', 'os', 'ir', 'pt', &
            'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn', &
         'fr', 'ra', 'ac', 'th', 'pa', 'u ', 'np', 'pu']
      character(len=2) element_
      integer i
      element_ = lowercase(element)
      atomicNumber = UnknownAtom
      do i = 1, size(elements)
         if (element_ == elements(i)) then
            atomicNumber = i
            exit
         end if
      end do
   end function atomicNumber
#endif

   subroutine vdw_reader(io, xc_functional, settings)
      use base, only: in_struct
      use reader_tags, only: writexmlincar, process_incar, open_incar_if_found, close_incar_if_found
      use tutor, only: vtutor
      type(in_struct), intent(in) :: io
      character(len=*), intent(in) :: xc_functional
      type(vdw_settings), intent(out) :: settings
      integer ierr, ivdw
      logical lopen
      !
      ivdw = 0
      call open_incar_if_found(io%iu5, lopen)
      call process_incar(lopen, io%iu0, io%iu5, 'IVDW', ivdw, ierr, writexmlincar)
      if (ivdw == 13) call read_dftd4_tags(settings)
      call close_incar_if_found(io%iu5)
      settings%ivdw = ivdw
      settings%io = io
      !
   contains
      !
      subroutine read_dftd4_tags(settings)
         type(vdw_settings), intent(out) :: settings
#ifdef DFTD4
         character(len=:), allocatable :: userDefined
         real(q) rdum
         logical unknown
         !
         settings = defaultSettings(xc_functional, ierr)
         unknown = ierr == UnknownFunctional
         userDefined = ''
         !
         associate(param => settings%parameters)
            call process_incar(lopen, io%iu0, io%iu5, 'VDW_S6', rdum, ierr, writexmlincar)
            if (ierr == 0) call overwriteIfSane('VDW_S6', [0.0_q, 2.0_q], rdum, param%s6, userDefined)
            call process_incar(lopen, io%iu6, io%iu5, 'VDW_S8', rdum, ierr, writexmlincar)
            if (ierr == 0) call overwriteIfSane('VDW_S8', [0.0_q, 10.0_q], rdum, param%s8, userDefined)
            call process_incar(lopen, io%iu6, io%iu5, 'VDW_A1', rdum, ierr, writexmlincar)
            if (ierr == 0) call overwriteIfSane('VDW_A1', [0.0_q, 10.0_q], rdum, param%a1, userDefined)
            call process_incar(lopen, io%iu6, io%iu5, 'VDW_A2', rdum, ierr, writexmlincar)
            if (ierr == 0) call overwriteIfSane('VDW_A2', [0.0_q, 10.0_q], rdum, param%a2, userDefined)
         end associate
         !
         if (unknown) then
            if (len_trim(userDefined) > 0) then
               call vtutor%warning("Using user supplied parameters" // userDefined // &
                  "\nbecause the xc functional is not supported. Required parameters not &
                  &defined in the INCAR file will be replaced with PBE defaults ones.")
            else
               call vtutor%error('Unsupported xc functional, please specify VDW_A1, &
                  &VDW_A2, VDW_S6, and VDW_S8 in the INCAR file.')
            end if
         else if (len_trim(userDefined) > 0) then
            call vtutor%alert("Specifying the vdW parameters" // userDefined // &
               "\nin the INCAR file will overwrite the defaults for " // trim(xc_functional) // &
               ". Please make sure that is what you intended.")
         end if
#else
         call vtutor%error("The code was not compiled with DFTD4 support. Please &
            recompile the code adding -DDFTD4 to the CPP_OPTIONS and add the &
            necessary include and link parameters.")
#endif
      end subroutine read_dftd4_tags
      !
      subroutine overwriteIfSane(tag, saneInterval, newValue, useValue, userDefined)
         use string, only: str
         character(len=*), intent(in) :: tag
         real(q), intent(in) :: saneInterval(2), newValue
         real(q), intent(inout) :: useValue
         character(len=:), allocatable, intent(inout) :: userDefined
         if ((minval(saneInterval) <= newValue).and.(newValue <= maxval(saneInterval))) then
            userDefined = userDefined // '\n' // tag // ' = ' // str(newValue)
            useValue = newValue
         else
            call vtutor%alert(tag // ' is not reasonable. Taking default.')
         end if
      end subroutine overwriteIfSane
      !
   end subroutine vdw_reader

   type(vdw_settings) function defaultSettings(xc_functional, ierr) result (settings)
#ifdef DFTD4
      use dfuncpar, only: d4par
      use mctc_environment, only: mctc_logger
      type(mctc_logger) env
#endif
      character(len=*), intent(in) :: xc_functional
      integer, intent(out) :: ierr
      ierr = 0
#ifdef DFTD4
      settings%options = defaultOptions()
      call d4par(xc_functional, settings%parameters, settings%options%lmbd, env)
      if (.not.env%sane) then
         ierr = UnknownFunctional
         associate(param => settings%parameters)
            ! use PBE defaults
            param%s6 = 1.0_q
            param%s8 = 0.95948085_q
            param%a1 = 0.38574991_q
            param%a2 = 4.80688534_q
         end associate
      end if
#endif
   end function defaultSettings

#ifdef DFTD4
   type(dftd_options) function defaultOptions() result (options)
      options%lmbd        = 3          ! Axilrod-Teller-Muto term
      options%refq        = 5          ! EEQ charge model
      options%wf          = 6.0_q      ! Gaussian weighting factor
      options%g_a         = 3.0_q      ! charge scaling factor height
      options%g_c         = 2.0_q      ! charge scaling factor steepness
      options%lmolpol     = .false.    ! skip molecular polarizability calculation
      options%lenergy     = .true.     ! calculate dispersion energy
      options%lgradient   = .true.     ! calculate dispersion energy derivatives
      options%lhessian    = .false.    ! skip hessian calculation
      options%print_level = 2          ! dftd4 print level
   end function defaultOptions
#endif

end module vdwd4
