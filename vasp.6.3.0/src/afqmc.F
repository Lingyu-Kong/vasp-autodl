#include "symbol.inc"
module afqmc

   use afqmc_struct
   use base, only: q, in_struct, symmetry, info_struct
   use hamil_struct_def, only: ham_handle
   use lattice, only: latt
   use meta, only: tau_handle
   use mgrid, only: grid_3d, transit
   use mkpoints, only: kpoints_struct
   use nonl, only: nonl_struct
   use nonlr, only: nonlr_struct
   use poscar, only: type_info
   use pseudo, only: potcar
   use wave, only: wavedes, wavespin

   implicit none

   private
   public afqmc_settings, afqmc_reader, afqmc_propagation

contains

   subroutine afqmc_reader(incar, num_images, settings)
      use incar_reader, only: incar_file, process_incar
      use string, only: lowercase
      use tutor, only: vtutor
      type(incar_file), intent(inout) :: incar
      integer, intent(in) :: num_images
      type(afqmc_settings), intent(out) :: settings
      character(len=:), allocatable :: algo
      !
      algo = ''
      call process_incar(incar, 'ALGO', algo)
      if (lowercase(trim(adjustl(algo))) == 'afqmc') then
          call vtutor%error("AFQMC is not implemented yet.")
      end if
      !
   end subroutine afqmc_reader

   subroutine afqmc_propagation(settings, kineden, hamiltonian, p, w, wdes, nonlr_s, nonl_s, &
         latt_cur, t_info, info, io, grid, grid_soft, gridc, gridus, c_to_us, &
         soft_to_c, symm, chtot, cvtot, cstrf, chden, dencor, sv, cdij, cqij, &
         crhode, rholm, lmdim, irdmax, n_mix_paw, kpoints)
      ! parameters
      type(afqmc_settings), intent(in) :: settings
      type(tau_handle) kineden
      type(ham_handle) hamiltonian
      type(potcar) p(:)
      type(wavespin) w
      type(wavedes) wdes
      type(nonlr_struct) nonlr_s
      type(nonl_struct) nonl_s
      type(latt) latt_cur
      type(type_info) t_info
      type(info_struct) info
      type(in_struct) io
      type(grid_3d) grid, grid_soft, gridc, gridus
      type(transit) c_to_us, soft_to_c
      type(symmetry) symm
      type(kpoints_struct) kpoints
      complex(q) chtot(:,:), cvtot(:,:), cstrf(:,:), chden(:,:)
      RGRID dencor(:), sv(:,:)
      OVERLAP cdij(:,:,:,:), cqij(:,:,:,:), crhode(:,:,:,:)
      real(q) rholm(:,:)
      integer lmdim, irdmax, n_mix_paw
      !
   end subroutine afqmc_propagation

end module afqmc
