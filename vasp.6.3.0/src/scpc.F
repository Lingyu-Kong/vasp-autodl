#include "symbol.inc"

#if defined(SCPC) && !defined(MPI)
#error SCPC uses MPI commands
#endif

module scpc

   use prec,    only: q
   use mgrid,   only: grid_3d
   use lattice, only: latt
   use pseudo,  only: potcar
   use poscar,  only: type_info
#ifdef MPI
   use mpimy
#endif

   implicit none

   real(q), parameter :: ryd2ev = 13.605698066_q

   ! contains the settings required to determine how the potential correction
   ! is conducted
   type scpc_settings
      logical :: active = .false.
      integer :: unitout = 7777
      integer :: nstep = 0
      integer :: ntarget = 1
      character(7) :: fout = "SCPCOUT"
      character(80) :: fpot = "REFPOT", fchg = "REFCHG"
      character(80) :: release = 'release: vasp(6.3.0) scpc(rev7)'
      real(q) :: qtot  = 0.0_q
      real(q) :: zlow  = 0.0_q, zhig  = 1.0_q
      real(q) :: diel  = 1.0_q
      real(q) :: broad = 0.4_q
      real(q) :: dnx = 0.0_q, dny = 0.0_q, dnz = 0.0_q
      integer :: ngx = 0, ngy = 0, ngz = 0, nxyz = 0
      integer :: mgx = -1, mgy = -1, mgz = -1, mxyz = 0
      real(q) :: dmx = 0.0_q, dmy = 0.0_q, dmz = 0.0_q
      real(q) :: omega = 0.0_q
      real(q) :: alat = 0.0_q, blat = 0.0_q, clat = 0.0_q
      real(q) :: rx0 = 0.0_q, ry0 = 0.0_q, rz0 = 0.0_q
      real(q), allocatable :: rchg(:,:,:), mchg(:,:,:), rhotot(:,:,:)
      real(q), allocatable :: rpot(:,:,:), mpot(:,:,:), vhar_loc(:,:,:)
      real(q), allocatable :: meps(:,:,:), mdeps(:,:,:)
      real(q), allocatable :: vadd(:,:,:)
      complex(q), allocatable :: cvhar(:)
      logical :: prtx = .false., prty = .false., prtz = .true.
      real(q) :: ecor1 = 0.0_q, ecor2 = 0.0_q, ecor_tot = 0.0_q
      real(q) :: mixrho = 1.00_q
      logical :: isbulk = .false.
      real(q) :: rxcut = -1.0_q, rycut = -1.0_q, rzcut = -1.0_q
   end type scpc_settings

   ! container for VASP
   type(scpc_settings) scpc_set

   private
   public scpc_settings, scpc_set, scpc_apply, scpc_reader

contains

   subroutine scpc_reader(incar, settings)
      use incar_reader, only: incar_file, process_incar
      use tutor, only: vtutor
      type(incar_file), intent(inout) :: incar
      type(scpc_settings), intent(out) :: settings
      call process_incar(incar, 'SCPC/USE', settings%active)
      if (settings%active) then
#ifndef SCPC
         call vtutor%error("The code was not compiled with -DSCPC option. &
            &Please recompile the code with it if you want to use the SCPC &
            &package. Keep in mind you need to link to the DL_MG Poisson &
            &solver and the PSPFFT isolate potential solver.")
#endif
         call process_incar(incar, 'SCPC/PRTX', settings%prtx)
         call process_incar(incar, 'SCPC/PRTY', settings%prty)
         call process_incar(incar, 'SCPC/PRTZ', settings%prtz)
         call process_incar(incar, 'SCPC/IN', settings%ntarget)
         call process_incar(incar, 'SCPC/QTOT', settings%qtot)
         call process_incar(incar, 'SCPC/ZLOW', settings%zlow)
         call process_incar(incar, 'SCPC/ZHIG', settings%zhig)
         call process_incar(incar, 'SCPC/DIEL', settings%diel)
         call process_incar(incar, 'SCPC/BROAD', settings%broad)
         call process_incar(incar, 'SCPC/RX0', settings%rx0)
         call process_incar(incar, 'SCPC/RY0', settings%ry0)
         call process_incar(incar, 'SCPC/RZ0', settings%rz0)
         call process_incar(incar, 'SCPC/MGX', settings%mgx)
         call process_incar(incar, 'SCPC/MGY', settings%mgy)
         call process_incar(incar, 'SCPC/MGZ', settings%mgz)
         call process_incar(incar, 'SCPC/RXCUT', settings%rxcut)
         call process_incar(incar, 'SCPC/RYCUT', settings%rycut)
         call process_incar(incar, 'SCPC/RZCUT', settings%rzcut)
      end if
   end subroutine scpc_reader

   subroutine scpc_apply(gridc, latt_cur, chtot, p, t_info, cstrf, &
      cvtot, escpc)
      type(grid_3d), intent(in) :: gridc
      type(latt), intent(in) :: latt_cur
      type(potcar), intent(in) :: p(:)
      type(type_info), intent(in) :: t_info
      complex(q), intent(in) :: chtot(:,:), cstrf(:,:)
      complex(q), intent(inout) :: cvtot(:,:)
      real(q), intent(out) :: escpc
      integer :: node_me, ionode, iout
      integer, save :: ncount = 0

      escpc = 0.0_q
      if (.not.scpc_set%active) return

#ifdef MPI
      node_me = gridc%comm%node_me
      ionode  = gridc%comm%ionode
#endif
      CALLMPI(M_bcast_i(gridc%comm, scpc_set%ntarget, 1))
      CALLMPI(M_bcast_l(gridc%comm, scpc_set%active, 1))
      CALLMPI(M_barrier(gridc%comm))

      ncount = ncount + 1
      if (ncount == 1) then
         io_begin
            call print_header()
         io_end
         call scpc_setup(gridc, latt_cur)
         call scpc_alloc(gridc)
         call scpc_get_reference(gridc)
         call scpc_build_diel(gridc)
      else if (ncount == scpc_set%ntarget) then
         io_begin
            write(*,'(/,2x,a,/)') "<<<  Self-consistent Potential Correction is ON  >>>"
         io_end
      end if
      CALLMPI(M_barrier(gridc%comm))

      if (ncount >= scpc_set%ntarget) then

         iout = scpc_set%unitout
         io_begin
            open(iout, file=scpc_set%fout, position='append')
            write(iout,'(//,4x,a,i5,a,/)') "<<<< SCPC Cycle:", ncount," >>>>"
            close(iout)
         io_end

         call scpc_build_model(gridc, p, latt_cur, t_info, chtot, cstrf)
         call scpc_build_correction(gridc)
         call scpc_addpot(gridc, cvtot)
         call scpc_prtall(gridc)

      end if

      escpc = scpc_set%ecor_tot
      CALLMPI(M_barrier(gridc%comm))

   end subroutine scpc_apply

   subroutine print_header()
      use tutor, only: vtutor
      integer iout
      iout = scpc_set%unitout
      open(iout, file=scpc_set%fout, position='append')
      write(iout,*)
      write(iout,*) "  <<<< SELF-CONSISTENT POTENTIAL CORRECTION - SCPC >>>>  "
      write(iout,*)
      write(iout,*) trim(adjustl(scpc_set%release))
      write(iout,*)
      write(iout,'(a)') "  Please cite:"
      write(iout,'(a)') "  M. Chagas da Silva, M. Lorke, B. Aradi, M. Farzalipour Tabriz,"
      write(iout,'(a)') "  T. Frauenheim, A. Rubio, D. Rocca, and P. Deak,"
      write(iout,'(a)') '  "Self-consistent potential correction for charged periodic systems",'
      write(iout,'(a)') "  Phys. Rev. Lett. 126, 076401 (2021)."
      write(iout,*)
      write(iout,'(a,i3)') "  SCPC will start at cycle :", scpc_set%ntarget
      write(iout,*)
      write(iout,'(a,f4.1)') "  Total model charge:", scpc_set%qtot
      write(iout,*)
      write(iout,'(a,f6.2)') "  Interface coord. low bound:", scpc_set%zlow
      write(iout,'(a,f6.2)') "  Interface cood. high bound:", scpc_set%zhig
      write(iout,'(a,f6.2)') "  Interface broadening:", scpc_set%broad
      write(iout,'(a,f6.2)') "  Interface avg. dielectric:", scpc_set%diel
      close(iout)
   end subroutine print_header

   subroutine scpc_setup(gridc, latt_cur)
      type(grid_3d), intent(in) :: gridc
      type(latt), intent(in) :: latt_cur
      integer :: iout, node_me, ionode

      scpc_set%qtot = -scpc_set%qtot

      scpc_set%alat = latt_cur%a(1,1)
      scpc_set%blat = latt_cur%a(2,2)
      scpc_set%clat = latt_cur%a(3,3)

      scpc_set%omega = latt_cur%omega

      ! external dense grid
      scpc_set%ngx = gridc%ngx
      scpc_set%ngy = gridc%ngy
      scpc_set%ngz = gridc%ngz
      scpc_set%nxyz = (scpc_set%ngx*scpc_set%ngy*scpc_set%ngz)

      scpc_set%dnx = scpc_set%alat / scpc_set%ngx
      scpc_set%dny = scpc_set%blat / scpc_set%ngy
      scpc_set%dnz = scpc_set%clat / scpc_set%ngz

      if (scpc_set%mgx == -1) then
         scpc_set%mgx = int(0.60 * gridc%ngx)
         call getgrid(scpc_set%mgx)
      end if

      if (scpc_set%mgy == -1) then
         scpc_set%mgy = int(0.60 * gridc%ngy)
         call getgrid(scpc_set%mgy)
      end if

      if (scpc_set%mgz == -1) then
         scpc_set%mgz = int(0.60 * gridc%ngz)
         call getgrid(scpc_set%mgz)
      end if

      scpc_set%mxyz = scpc_set%mgx * scpc_set%mgy * scpc_set%mgz
      scpc_set%dmx = scpc_set%alat / scpc_set%mgx
      scpc_set%dmy = scpc_set%blat / scpc_set%mgy
      scpc_set%dmz = scpc_set%clat / scpc_set%mgz

#ifdef MPI
      node_me = gridc%comm%node_me
      ionode  = gridc%comm%ionode
#endif

      io_begin
         iout = scpc_set%unitout
         open(iout, file=scpc_set%fout, position='append')
         write(iout,*)
         write(iout,'(a,e20.10)') "  Cell volume:", scpc_set%omega
         write(iout,'(a,3f12.6)') "  Cell lattice:", scpc_set%alat, scpc_set%blat, scpc_set%clat
         write(iout,*)
         write(iout,'(a,3i5)') "  Model grid:", scpc_set%mgx, scpc_set%mgy, scpc_set%mgz
         write(iout,'(a,3f12.6)') "  Model step size grid:", scpc_set%dmx, scpc_set%dmy, scpc_set%dmz
         write(iout,*)
         close(iout)
      io_end

   end subroutine scpc_setup

   subroutine getgrid(mdim)
      integer, intent(inout) :: mdim
      integer :: i, j, ilow, jlow
      real(q) :: x0, x1, x2, rlow
      x0 = dble(mdim)
      rlow = huge(1.0_q)
      do i = 1, 20
         do j = 1, 20
            x1 = 2.0_q**i * dble(j)
            x2 = abs(x0 - x1)
            if (x2 <= rlow) then
               rlow = x2
               ilow = i
               jlow = j
            end if
         end do
      end do
      mdim = int(2.0_q**ilow * dble(jlow)) + 1
   end subroutine getgrid

   subroutine scpc_alloc(gridc)
      type(grid_3d), intent(in) :: gridc
      integer :: nx, ny, nz
      integer :: node_me, ionode

#ifdef MPI
      node_me = gridc%comm%node_me
      ionode  = gridc%comm%ionode
#endif

      nx = scpc_set%mgx
      ny = scpc_set%mgy
      nz = scpc_set%mgz

      io_begin
         if (.not.allocated(scpc_set%rchg)) allocate(scpc_set%rchg(nx, ny, nz))
         if (.not.allocated(scpc_set%mchg)) allocate(scpc_set%mchg(nx, ny, nz))
         if (.not.allocated(scpc_set%rhotot)) allocate(scpc_set%rhotot(nx, ny, nz))
         if (.not.allocated(scpc_set%rpot)) allocate(scpc_set%rpot(nx, ny, nz))
         if (.not.allocated(scpc_set%mpot)) allocate(scpc_set%mpot(nx, ny, nz))
         if (.not.allocated(scpc_set%vhar_loc)) allocate(scpc_set%vhar_loc(nx, ny, nz))
         if (.not.allocated(scpc_set%meps)) allocate(scpc_set%meps(nx, ny, nz))
         if (.not.allocated(scpc_set%mdeps)) allocate(scpc_set%mdeps(nx, ny, nz))
         if (.not.allocated(scpc_set%vadd)) allocate(scpc_set%vadd(nx, ny, nz))
         scpc_set%rchg = 0.0_q
         scpc_set%mchg = 0.0_q
         scpc_set%rhotot = 0.0_q
         scpc_set%rpot = 0.0_q
         scpc_set%mpot = 0.0_q
         scpc_set%vhar_loc = 0.0_q
         scpc_set%meps = 1.0_q
         scpc_set%mdeps = 0.0_q
         scpc_set%vadd = 0.0_q
      io_end
      if (.not.allocated(scpc_set%cvhar)) allocate(scpc_set%cvhar(gridc%mplwv))
      scpc_set%cvhar = 0.0_q

      CALLMPI(M_barrier(gridc%comm))

   end subroutine scpc_alloc

   subroutine scpc_dealloc(gridc)
      type(grid_3d), intent(in) :: gridc
      if (allocated(scpc_set%rchg)) deallocate(scpc_set%rchg)
      if (allocated(scpc_set%mchg)) deallocate(scpc_set%mchg)
      if (allocated(scpc_set%rhotot)) deallocate(scpc_set%rhotot)
      if (allocated(scpc_set%rpot)) deallocate(scpc_set%rpot)
      if (allocated(scpc_set%mpot)) deallocate(scpc_set%mpot)
      if (allocated(scpc_set%meps)) deallocate(scpc_set%meps)
      if (allocated(scpc_set%mdeps)) deallocate(scpc_set%mdeps)
      if (allocated(scpc_set%cvhar)) deallocate(scpc_set%cvhar)
      CALLMPI(M_barrier(gridc%comm))
   end subroutine scpc_dealloc

   subroutine scpc_get_reference(gridc)
      use tutor, only: vtutor
      type(grid_3d), intent(in) :: gridc
      logical :: ltest
      integer :: node_me, ionode
      integer :: nx, ny, nz, mx, my, mz
      integer :: ix, iy, iz, j
      real(q) :: w1(scpc_set%nxyz), w2(scpc_set%mxyz)
      real(q) :: r3d(scpc_set%mgx, scpc_set%mgy, scpc_set%mgz)
      real(q) :: dnx, dny, dnz, dmx, dmy, dmz
      real(q) :: tmp, r1, r2
      integer :: iout, ioer

#ifdef MPI
      node_me = gridc%comm%node_me
      ionode = gridc%comm%ionode
#endif

      io_begin

         ! check for the external reference files first
         inquire(file=scpc_set%fchg, exist=ltest)
         if (.not.ltest) &
            call vtutor%error("scpc error: missing reference file for charge density")
         inquire(file=scpc_set%fpot, exist=ltest)
         if (.not.ltest) &
            call vtutor%error("scpc error: missing reference file for hartree potential")

         ! dense grid
         nx = scpc_set%ngx
         ny = scpc_set%ngy
         nz = scpc_set%ngz
         dnx = scpc_set%dnx
         dny = scpc_set%dny
         dnz = scpc_set%dnz
         ! reduced grid
         mx = scpc_set%mgx
         my = scpc_set%mgy
         mz = scpc_set%mgz
         dmx = scpc_set%dmx
         dmy = scpc_set%dmy
         dmz = scpc_set%dmz

         iout = scpc_set%unitout

         ! reference charge density must be divided by the unit cell volume
         open(213, file=scpc_set%fchg)

         do j=1,8
            read(213,*)
         end do

         do
            read(213,'(3i5)', iostat=ioer) ix, iy, iz
            if (ioer == 0) exit
         end do

         read(213,*) ix, iy, iz

         !if((ix/=nx).or.(iy/=ny).or.(iz/=nz)) then
         !  write(*,*) "scpc error: reference charge grid is not the same in actual system"
         !  call mpi_abort(mpi_comm_world,1,ier)
         !end if
         ! loop order z , y , x(fast inner loop)
         w1 = 0.0_q
         do j = 1, scpc_set%nxyz
            if ((mod(j,5) == 0).or.(j == scpc_set%nxyz)) then
               read(213, 104, advance='yes') tmp
            else
               read(213, 104, advance='no') tmp
            end if
            w1(j) = tmp / scpc_set%omega
104         format(1(1x,e17.11))
         end do
         close(213)
         ! reduce grid
         w2 = 0.0_q
         call scpc_mltscale(w1, nx, ny, nz, dnx, dny, dnz, &
                            w2, mx, my, mz, dmx, dmy, dmz)
         ! renormalize the charge density to the new reduced grid
         r1 = sum(w1) / dble(scpc_set%nxyz)
         r2 = sum(w2) / dble(scpc_set%mxyz)
         w2 = w2 * (r1 / r2)

         r3d = 0.0_q
         j = 1
         do iz = 1, mz
            do iy = 1, my
               do ix = 1,mx
                  r3d(ix,iy,iz) = w2(j)
                  j = j + 1
               end do
            end do
         end do
         scpc_set%rchg = r3d

         open(iout, file=scpc_set%fout, position='append')
         write(iout,*)
         tmp = scpc_set%omega / scpc_set%mxyz
         tmp = sum(scpc_set%rchg) * tmp
         write(iout,'(a,e20.10)') "  Integrated reference rho:", tmp
         close(iout)

         ! reference potential
         open(213, file=scpc_set%fpot)

         do j=1,8
            read(213,*)
         end do

         do
            read(213,'(3i5)', iostat=ioer) ix, iy, iz
            if (ioer == 0) exit
         end do

         read(213,*) ix, iy, iz

         !if((ix/=nx).or.(iy/=ny).or.(iz/=nz)) then
         !  write(*,*) "scpc error: reference potential grid is not the same in actual system"
         !  call mpi_abort(mpi_comm_world,1,ier)
         !end if
         ! loop order z , y , x(fast inner loop)
         w1 = 0.0_q
         do j = 1, scpc_set%nxyz
            if ((mod(j,5) == 0).or.(j == scpc_set%nxyz)) then
               read(213, 104, advance='yes') tmp
            else
               read(213, 104, advance='no') tmp
            end if
            w1(j) = tmp
         end do
         close(213)
         ! reduce grid
         w2 = 0.0_q
         call scpc_mltscale(w1, nx, ny, nz, dnx, dny, dnz, &
                            w2, mx, my, mz, dmx, dmy, dmz)
         ! setting the average potential to zero
         r2 = sum(w2) / dble(scpc_set%mxyz)
         w2 = w2 - r2

         r3d = 0.0_q
         j = 1
         do iz = 1, mz
            do iy = 1, my
               do ix = 1, mx
                  r3d(ix,iy,iz) = w2(j)
                  j = j + 1
               end do
            end do
         end do
         scpc_set%rpot = r3d

         open(iout, file=scpc_set%fout, position='append')
         write(iout,*)
         tmp = scpc_set%omega / scpc_set%mxyz
         tmp = sum(scpc_set%rpot) * tmp
         write(iout,'(a,e20.10)') "  Integrated reference vhar:", tmp
         close(iout)

      io_end

      CALLMPI(M_barrier(gridc%comm))

   end subroutine scpc_get_reference

   subroutine scpc_mltscale(w1, nx, ny, nz, dnx, dny, dnz, &
      w2, mx, my, mz, dmx, dmy, dmz)
      ! 3DGRID INTERPOLATION WITH pbc CONDITIONS
      integer, intent(in) :: nx, ny, nz, mx, my, mz
      real(q), intent(in) :: dnx, dny, dnz, dmx, dmy, dmz
      real(q), intent(in) :: w1(nx * ny * nz)
      real(q), intent(out) :: w2(mx * my * mz)
      real(q) :: xo, yo, zo
      real(q) :: xi, yi, zi
      integer :: ix, iy, iz
      integer :: jx, jy, jz
      integer :: k
      !
      w2 = dzero
      xi = dzero
      do ix = 1, mx
         call validx(xi, dnx, nx, xo, jx)
         yi = dzero
         do iy = 1, my
            call validx(yi, dny, ny, yo, jy)
            zi = dzero
            do iz = 1, mz
               k = idx(ix, iy, iz, mx, my, mz)
               call validx(zi, dnz, nz, zo, jz)
               call interpol(w2(k), w1, jx, jy, jz, xo, yo, zo, nx, ny, nz)
               zi = zi + dmz
            end do
            yi = yi + dmy
         end do
         xi = xi + dmx
      end do
      !
   end subroutine

   integer function idx(ix, iy, iz, nx, ny, nz)
      integer, intent(in) :: ix, iy, iz
      integer, intent(in) :: nx, ny, nz
      idx = modulo(ix - 1, nx) + 1
      idx = idx + modulo(iy - 1, ny) * nx
      idx = idx + modulo(iz - 1, nz) * ny * nx
   end function idx

   subroutine validx(u, du, nu, res, idu)
      real(q), intent(in) :: u, du
      integer, intent(in) :: nu
      real(q), intent(out) :: res
      integer, intent(out) :: idu
      res = u / du
      idu = min(int(res) + 1, nu - 1)
      res = res - dble(idu - 1)
   end subroutine validx

   subroutine interpol(res, vec, ir1, ir2, ir3, tr1, tr2, tr3, nr1, nr2, nr3)
      real(q), intent(inout) :: res
      real(q), intent(in) :: vec(:)
      real(q), intent(in) :: tr1, tr2, tr3
      integer, intent(in) :: ir1, ir2, ir3, nr1, nr2, nr3
      integer :: l1, l2, l3, il1, il2, il3
      real(q) :: f0, f1, f2, f3
      real(q) :: p11, p12, p13, p21, p22, p23
      real(q) :: v00, v10, v11, v20, v21, v30, v31
      real(q), parameter :: half = 0.5d0
      do l1 = 0, 1
         il1 = ir1 + l1
         p11 = pol1(tr1, l1)
         p21 = pol2(tr1, l1)
         do l2 = 0, 1
            il2 = ir2 + l2
            p12 = pol1(tr2, l2)
            p22 = pol2(tr2, l2)
            do l3 = 0, 1
               il3 = ir3 + l3
               p13 = pol1(tr3, l3)
               p23 = pol2(tr3, l3)
               v00 = vec(idx(il1  , il2  , il3  , nr1, nr2, nr3))
               v10 = vec(idx(il1+1, il2  , il3  , nr1, nr2, nr3))
               v11 = vec(idx(il1-1, il2  , il3  , nr1, nr2, nr3))
               v20 = vec(idx(il1  , il2+1, il3  , nr1, nr2, nr3))
               v21 = vec(idx(il1  , il2-1, il3  , nr1, nr2, nr3))
               v30 = vec(idx(il1  , il2  , il3+1, nr1, nr2, nr3))
               v31 = vec(idx(il1  , il2  , il3-1, nr1, nr2, nr3))
               f0  = v00 * p11 * p12 * p13
               f1  = half * (v10 - v11) * p21 * p12 * p13
               f2  = half * (v20 - v21) * p11 * p22 * p13
               f3  = half * (v30 - v31) * p11 * p12 * p23
               res = res + f0 + f1 + f2 + f3
            end do
         end do
      end do
   end subroutine interpol

   real(q) function pol1(x, fac)
      real(q), intent(in) :: x
      integer, intent(in) :: fac
      real(q) :: x2, x3
      real(q), parameter :: two = 2.0_q, three = 3.0_q
      x2 = x * x
      x3 = x2 * x
      pol1 = three * x2 - two * x3
      if (fac == 0) pol1 = done - pol1
   end function pol1

   real(q) function pol2(x, fac)
      real(q), intent(in) :: x
      integer, intent(in) :: fac
      real(q) :: x2, x3
      x2 = x * x
      x3 = x2 * x
      pol2 = -x2 + x3
      if (fac == 0) pol2 = x + pol2 - x2
   end function pol2

   subroutine scpc_build_diel(gridc)
      type(grid_3d), intent(in) :: gridc
      real(q), parameter :: c0 = 0.25_q
      real(q) :: t0, t1, t2
      real(q) :: r0, r1, r2
      real(q) :: z, beta
      integer :: ix, iy, iz
      integer :: node_me, ionode

#ifdef MPI
      node_me = gridc%comm%node_me
      ionode = gridc%comm%ionode
#endif

      io_begin

         if ((scpc_set%zlow == dzero).and.(scpc_set%zhig == done)) then
            scpc_set%meps = scpc_set%diel
            scpc_set%mdeps = dzero
            scpc_set%isbulk = .true.
         else
            beta = done / scpc_set%broad
            t0 = scpc_set%zlow*scpc_set%clat
            t1 = scpc_set%zhig*scpc_set%clat
            t2 = done / scpc_set%diel
            z  = 0.0_q
            do iz = 1, scpc_set%mgz
               r0 = done + erf((z - t0) * beta)
               r1 = done + erf((t1 - z) * beta)
               r2 = c0 * (t2 - done) * r0 * r1
               r2 = done + r2
               scpc_set%meps(:,:,iz) = done / r2
               z = z + scpc_set%dmz
            end do
            call scpc_numderv('T', 2, scpc_set%dmz, scpc_set%mgz, &
               scpc_set%meps(1,1,:), scpc_set%mdeps(1,1,:))
            do ix=1,scpc_set%mgx
               do iy=1,scpc_set%mgy
                  do iz=1,scpc_set%mgz
                     scpc_set%mdeps(ix,iy,iz) = scpc_set%mdeps(1,1,iz)
                  end do
               end do
            end do
         end if

         if ((scpc_set%prtz).and.(scpc_set%zlow /= dzero).and.(scpc_set%zhig /= done)) then
            call scpc_averages(gridc,'Z','z-diel.dat', scpc_set%meps, scpc_set%dmz, 1.0_q)
         end if

      io_end

      CALLMPI(M_barrier(gridc%comm))

   end subroutine scpc_build_diel

   subroutine scpc_numderv(pbc, ityp, h, n, v, dfv)
      character(len=1), intent(in) :: pbc
      integer, intent(in) :: ityp
      real(q), intent(in) :: h
      integer, intent(in) :: n
      real(q), intent(in) :: v(n)
      real(q), intent(out) :: dfv(n)
      real(q) :: fac
      real(q), parameter :: k12 = 12.0_q, k08 = 8.0_q
      integer :: i
      if (ityp == 2) then
         fac = 1.0_q / (2.0_q * h)
         do i = 2, n - 1
            dfv(i) = fac * (v(i + 1) - v(i - 1))
         end do
         if (pbc == 'T') then
            dfv(1) = fac * (v(2) - v(n))
            dfv(n) = fac * (v(1) - v(n - 1))
         end if
         dfv(1) = dfv(2)
         dfv(n) = dfv(n - 1)
      else
         fac = 1.0_q / (k12 * h)
         do i = 3, n - 2
            dfv(i) = fac * (v(i - 2) - k08 * v(i - 1) + k08 * v(i + 1) - v(i + 2))
         end do
         if (pbc == 'T') then
           dfv(1)   = fac * (v(n - 1) - k08*v(n)     + k08 * v(2) - v(3))
           dfv(2)   = fac * (v(n)     - k08*v(1)     + k08 * v(3) - v(4))
           dfv(n)   = fac * (v(n - 2) - k08*v(n - 1) + k08 * v(1) - v(2))
           dfv(n-1) = fac * (v(n - 3) - k08*v(n - 2) + k08 * v(n) - v(1))
         end if
      end if
   end subroutine scpc_numderv

   subroutine scpc_build_model(gridc, p, latt_cur, t_info, chtot, cstrf)
      type(grid_3d), intent(in) :: gridc
      type(potcar), intent(in) :: p(:)
      type(latt), intent(in) :: latt_cur
      type(type_info), intent(in) :: t_info
      complex(q), intent(in) :: chtot(:,:), cstrf(:,:)
      real(q), allocatable, save :: rhosave(:,:,:)
      real(q) :: r1, mixing
      integer :: node_me, ionode
      integer :: ix, iy, iz
      integer, save :: imix = 0

#ifdef MPI
      node_me = gridc%comm%node_me
      ionode  = gridc%comm%ionode
#endif

      io_begin
         if (.not.allocated(rhosave)) then
            allocate(rhosave(scpc_set%mgx,scpc_set%mgy,scpc_set%mgz))
            rhosave = 0.0_q
            mixing = 1.0_q
         else
            mixing = scpc_set%mixrho
         end if
      io_end
      CALLMPI(M_barrier(gridc%comm))

      call scpc_get_vhar(gridc, p, latt_cur, t_info, chtot, cstrf, scpc_set%cvhar)
      call scpc_merge(scpc_set%omega, gridc, chtot, scpc_set%rhotot, 0, 'RHO  :')
      CALLMPI(M_barrier(gridc%comm))

      io_begin
         scpc_set%mchg = scpc_set%rhotot - scpc_set%rchg
         r1 = sum(scpc_set%mchg) * scpc_set%omega / dble(scpc_set%mxyz)
         scpc_set%mchg = scpc_set%mchg * scpc_set%qtot / r1
         !if(scpc_set%isbulk.and.(imix.eq.1) then
         if (imix == 1) then
           ! INITIAL TEST WAS DONE USING rxcut=rycut=rzcut=0.15
           if (scpc_set%rxcut /= -1.0_q) call myrhomixer(mixing, rhosave, scpc_set%mchg, scpc_set%rxcut, 'X')
           if (scpc_set%rycut /= -1.0_q) call myrhomixer(mixing, rhosave, scpc_set%mchg, scpc_set%rycut, 'Y')
           if (scpc_set%isbulk.and.(scpc_set%rzcut /= -1.0_q)) &
              call myrhomixer(mixing, rhosave, scpc_set%mchg, scpc_set%rzcut, 'Z')
           !scpc_set%mchg = mixing*scpc_set%mchg + (1.0_q-mixing)*rhosave
           r1 = sum(scpc_set%mchg) * scpc_set%omega / dble(scpc_set%mxyz)
           scpc_set%mchg = scpc_set%mchg * scpc_set%qtot / r1
         end if
         rhosave = scpc_set%mchg
      io_end
      CALLMPI(M_barrier(gridc%comm))

      imix = 1
      call scpc_merge(1.0_q, gridc, scpc_set%cvhar, scpc_set%mpot, 1, 'VHAR :')

      io_begin
         scpc_set%vhar_loc = scpc_set%mpot
         scpc_set%mpot = scpc_set%mpot - scpc_set%rpot
         r1 = sum(scpc_set%mpot) * scpc_set%omega / dble(scpc_set%mxyz)
         scpc_set%mpot = scpc_set%mpot - r1
      io_end
      CALLMPI(M_barrier(gridc%comm))

   end subroutine scpc_build_model

   subroutine scpc_get_vhar(gridc, p, latt_cur, t_info, chtot, cstrf, cvhar)
      type(grid_3d), intent(in) :: gridc
      type(latt), intent(in) :: latt_cur
      type(potcar), intent(in) :: p(:)
      type(type_info), intent(in) :: t_info
      complex(q), intent(in) :: chtot(:,:)
      complex(q), intent(out) :: cvhar(gridc%mplwv)
      complex(q), intent(in) :: cstrf(:,:)
      complex(q) :: cwork(gridc%mplwv), cwork1(gridc%mplwv)
      real(q) :: rdum
      integer :: i

      ! cvhar holds the electrostatic potential per mpi process
      cvhar = 0.0_q

      ! hartree potential
      cwork = 0.0_q
      call pothar(gridc, latt_cur, chtot, cwork, rdum)
      cvhar = cwork

      ! local ionic potential
      cwork = 0.0_q
      call potion(gridc, p, latt_cur, t_info, cwork, cwork1, cstrf, rdum)
      cvhar = cvhar + cwork

      CALLMPI(M_barrier(gridc%comm))

   end subroutine scpc_get_vhar

   subroutine scpc_merge(fac, gridc, cinp, rout, inorm, strg)
      real(q), intent(in) :: fac
      type(grid_3d), intent(in) :: gridc
      complex(q), intent(in) :: cinp(gridc%mplwv)
      real(q), intent(out) :: rout(:,:,:)
      integer, intent(in) :: inorm
      character(*), intent(in) :: strg
      real(q), allocatable :: cwrk(:), wrk(:)
      real(q), allocatable :: nwrk(:), mwrk(:)
      integer :: ix, iy, iz, j, k
      integer :: nx, ny, nz, nxy, mx, my, mz
      real(q) :: r1, r2, tmp, dnx, dny, dnz, dmx, dmy, dmz
      integer :: node_me, ionode, iout
!
#ifdef MPI
      node_me = gridc%comm%node_me
      ionode  = gridc%comm%ionode
#endif

      ! dense grid
      nx = gridc%ngx
      ny = gridc%ngy
      nz = gridc%ngz
      nxy = nx * ny
      dnx = scpc_set%dnx
      dny = scpc_set%dny
      dnz = scpc_set%dnz
      ! reduced grid
      mx = scpc_set%mgx
      my = scpc_set%mgy
      mz = scpc_set%mgz
      dmx = scpc_set%dmx
      dmy = scpc_set%dmy
      dmz = scpc_set%dmz

      iout = scpc_set%unitout

      ! dimreal(gridc%mplwv) == gridc%mplwv*2
      allocate(cwrk(DIMREAL(gridc%mplwv)))
      allocate(nwrk(nxy * nz))
      allocate(wrk(nxy))

      call rc_add(cinp, 1.0_q, cinp, 0.0_q, cwrk, gridc)
      call FFT3D(cwrk, gridc, 1)

      k = 1
      do iz = 1, nz
         call mrg_grid_rl_plane(gridc, wrk, cwrk, iz)
         io_begin
            do j = 1, nxy
               nwrk(k) = wrk(j) / fac
               k = k + 1
            end do
         io_end
      end do
      if(allocated(cwrk)) deallocate(cwrk)
      if(allocated(wrk))  deallocate(wrk)
!
      io_begin
         allocate(mwrk(scpc_set%mxyz))
         mwrk = 0.0_q
         call scpc_mltscale(nwrk, nx, ny, nz, dnx, dny, dnz, &
                            mwrk, mx, my, mz, dmx, dmy, dmz)
         ! fixing some numerical integration
         if (inorm == 0) then
            r1 = sum(nwrk) / dble(scpc_set%nxyz)
            r2 = sum(mwrk) / dble(scpc_set%mxyz)
            tmp = r1 / r2
            mwrk = mwrk * tmp
         else
            r1 = sum(mwrk) / dble(scpc_set%mxyz)
            mwrk = mwrk - r1
         end if

         rout = 0.0_q
         k = 1
         do iz = 1, mz
            do iy = 1, my
               do ix = 1, mx
                  rout(ix,iy,iz) = mwrk(k)
                  k = k + 1
               end do
            end do
         end do
         if(allocated(mwrk)) deallocate(mwrk)

         tmp = scpc_set%omega / scpc_set%mxyz
         tmp = sum(rout) * tmp

         open(iout, file=scpc_set%fout, position='append')
         write(iout,'(a,a,e20.10)') "  Integrated System ", trim(adjustl(strg)), tmp
         close(iout)
      io_end

      if(allocated(nwrk)) deallocate(nwrk)
      CALLMPI(M_barrier(gridc%comm))

   end subroutine scpc_merge

   subroutine myrhomixer(fac, mat0, mat1, rcut, axis)
      real(q), intent(in) :: fac
      real(q), intent(in) :: mat0(:,:,:)
      real(q), intent(inout) :: mat1(:,:,:)
      real(q), intent(in) :: rcut
      character(len=1), intent(in) :: axis
      integer :: ix, iy, iz, mx, my, mz
      integer :: ilow, ihig
      real(q) :: y1,y2

      mx = scpc_set%mgx
      my = scpc_set%mgy
      mz = scpc_set%mgz

      select case (axis)
      case ('x', 'X')
         ilow = int(floor(rcut * mx))
         ihig = mx - ilow
         do iz = 1, mz
            do iy = 1, my
               do ix = 1, mx
                  y1 = fac!1.0_q - dble(ix) / dble(ilow)
                  y2 = fac!dble(ix) / dble(ihig)
                  if (ix <= ilow) then
                     mat1(ix,iy,iz) = y1 * mat0(ix,iy,iz) + (1 - y1) * mat1(ix,iy,iz)
                  else if (ix >= ihig) then
                     mat1(ix,iy,iz) = y2 * mat0(ix,iy,iz) + (1 - y2) * mat1(ix,iy,iz)
                  end if
               end do
            end do
         end do
      case ('y', 'Y')
         ilow = int(floor(rcut * my))
         ihig = my - ilow
         do ix = 1, mx
            do iz = 1, mz
               do iy = 1, my
                  y1 = fac!1.0_q - dble(iy) / dble(ilow)
                  y2 = fac!dble(iz) / dble(ihig)
                  if (iy <= ilow) then
                     mat1(ix,iy,iz) = y1 * mat0(ix,iy,iz) + (1 - y1) * mat1(ix,iy,iz)
                  else if (iy >= ihig) then
                     mat1(ix,iy,iz) = y2 * mat0(ix,iy,iz) + (1 - y2) * mat1(ix,iy,iz)
                  end if
               end do
            end do
         end do
      case ('z', 'Z')
         ilow = int(floor(rcut * mz))
         ihig = mz - ilow
         do ix = 1, mx
            do iy = 1, my
               do iz = 1, mz
                  y1 = fac!1.0_q - dble(iz) / dble(ilow)
                  y2 = fac!dble(iz) / dble(ihig)
                  if (iz <= ilow) then
                     mat1(ix,iy,iz) = y1 * mat0(ix,iy,iz) + (1 - y1) * mat1(ix,iy,iz)
                  else if (iz >= ihig) then
                     mat1(ix,iy,iz) = y2 * mat0(ix,iy,iz) + (1 - y2) * mat1(ix,iy,iz)
                  end if
               end do
            end do
         end do
      end select
   end subroutine myrhomixer

   subroutine scpc_build_correction(gridc)
      type(grid_3d), intent(in) :: gridc
      real(q), allocatable :: viso(:,:,:), vper(:,:,:), vcor(:,:,:), raux(:,:,:)
!     real(q), allocatable, save :: vold(:,:,:)
      integer :: iout
      real(q) :: ravg, fac, ealig
      integer :: node_me, ionode
!     real :: mixing

#ifdef MPI
      node_me = gridc%comm%node_me
      ionode  = gridc%comm%ionode
#endif

      io_begin

         !if(.not.allocated(vold)) then
         !  allocate(vold(scpc_set%mgx,scpc_set%mgy,scpc_set%mgz))
         !  vold = 0.0_q
         !  mixing = 1.0_q
         !else
         !  mixing = scpc_set%mixpot
         !end if

         allocate(viso(scpc_set%mgx, scpc_set%mgy, scpc_set%mgz))
         allocate(vper(scpc_set%mgx, scpc_set%mgy, scpc_set%mgz))
         allocate(vcor(scpc_set%mgx, scpc_set%mgy, scpc_set%mgz))
         allocate(raux(scpc_set%mgx, scpc_set%mgy, scpc_set%mgz))

         viso = 0.0_q
         call scpc_build_isopot(gridc, viso)
         viso = viso * ryd2ev

         ravg = sum(scpc_set%mchg) / dble(scpc_set%mxyz)
         vper = scpc_set%mpot / ryd2ev ! Put in Rydberg
         call scpc_mgsolver(gridc, scpc_set%mchg - ravg, vper, scpc_set%meps, "PER")
         vper = vper * ryd2ev
         fac = sum(vper) / dble(scpc_set%mxyz)
         vper = vper - fac

         vcor = (viso - vper) / ryd2ev
         raux = sum(scpc_set%mchg) / dble(scpc_set%mxyz)

         !call scpc_averages(gridc,'z','Z-DIF-vcor.DAT',vcor*ryd2ev,scpc_set%dmz,1.0_q)

         call scpc_mgsolver(gridc, raux, vcor, scpc_set%meps, "COR")
         vcor = vcor * ryd2ev

         !call scpc_averages(gridc,'z','Z-RESOLV-vcor.DAT',vcor,scpc_set%dmz,1.0_q)

!         vold = mixing*scpc + (1.0_q-mixing)*vold

         !vcor_set%vadd = vcor !vold
         call vref_alignment(gridc, vcor, scpc_set%vadd, ealig)

         ealig = 0.50_q * ealig * scpc_set%qtot
         fac = 0.50_q * scpc_set%omega / dble(scpc_set%mxyz)
         scpc_set%ecor1 = sum(scpc_set%mchg * scpc_set%vadd) * fac
         scpc_set%ecor2 = 2.0_q * sum(scpc_set%rhotot * scpc_set%vadd) * fac
         scpc_set%ecor_tot = scpc_set%ecor1 - scpc_set%ecor2 - ealig

         iout = scpc_set%unitout
         open(iout, file=scpc_set%fout, position='append')
         write(iout,'(/,2X,a,f8.4)') "Mixing rho:", scpc_set%mixrho
         write(iout,'(/,2X,a,f8.4)') "Mixing rxcut:", scpc_set%rxcut
         write(iout,'(/,2X,a,f8.4)') "Mixing rycut:", scpc_set%rycut
         write(iout,'(/,2X,a,f8.4)') "Mixing rzcut:", scpc_set%rzcut
         !write(iout,'(/,2X,a,f8.4)') "Mixing pot:", scpc_set%mixpot
         !write(iout,'(/,2X,a,f8.4,a)') "Energy correction :", scpc_set%ecor1, " Ev"
         !write(iout,'(/,2X,a,f8.4,a)') "Double counting :", scpc_set%ecor2, " Ev"
         !write(iout,'(/,2X,a,f8.4,a)') "Alignment term", ealig, " Ev"
         write(iout,'(/,2X,a,f8.4,a)') "Energy Correction :", scpc_set%ecor1-ealig, " Ev"
         close(iout)

         call scpc_potalignment(gridc, scpc_set%mpot, vper)

         if (scpc_set%prtx) then
            call scpc_averages(gridc, 'X', 'x-vcor.dat', vcor, scpc_set%dmx, 1.0_q)
         end if
         if (scpc_set%prty) then
            call scpc_averages(gridc, 'Y', 'y-vcor.dat', vcor, scpc_set%dmy, 1.0_q)
         end if
         if (scpc_set%prtz) then
            call scpc_averages(gridc, 'Z', 'z-vcor.dat', vcor, scpc_set%dmz, 1.0_q)
         end if

         deallocate(viso)
         deallocate(vper)
         deallocate(vcor)
         deallocate(raux)

      io_end

      CALLMPI(M_barrier(gridc%comm))
      CALLMPI(M_bcast_d(gridc%comm, scpc_set%ecor_tot, 1))

   end subroutine scpc_build_correction

   subroutine scpc_build_isopot(gridc, mypot)
      type(grid_3d), intent(in) :: gridc
      real(q), intent(inout) :: mypot(scpc_set%mgx, scpc_set%mgy, scpc_set%mgz)
      real(q), parameter :: fac1 = -25.1327412287183_q * 0.50_q
      real(q), parameter :: upfac = 0.50d0
      integer, parameter :: nmax = 300
      real(q), dimension(:,:,:), allocatable :: raux, raux2, vaux
      real(q), dimension(:,:,:), allocatable :: dzraux2, dzvaux
      real(q), dimension(:,:,:), allocatable :: rho_tot, prefac
      integer :: ix, iy, iz, iter, mpime, mpi_err
      real(q) :: tmp1,tmp2, e_old, e_new, delta_e
      real(q) :: k0, k1
      integer :: node_me, ionode

#ifdef MPI
      node_me = gridc%comm%node_me
      ionode  = gridc%comm%ionode
#endif

      io_begin

         allocate(vaux(scpc_set%mgx, scpc_set%mgy, scpc_set%mgz))
         allocate(raux(scpc_set%mgx, scpc_set%mgy, scpc_set%mgz))
         raux = scpc_set%mchg / scpc_set%diel
         raux = raux * fac1

         if ((scpc_set%zlow == 0.0_q).and.(scpc_set%zhig == 1.0_q)) then
            vaux = 0.0_q
            call scpc_solve_isopot(raux, vaux)

         else

            allocate(raux2(scpc_set%mgx,scpc_set%mgy,scpc_set%mgz))
            allocate(dZraux2(scpc_set%mgx,scpc_set%mgy,scpc_set%mgz))
            allocate(dZvaux(scpc_set%mgx,scpc_set%mgy,scpc_set%mgz))
            allocate( rho_tot(scpc_set%mgx,scpc_set%mgy,scpc_set%mgz) )
            allocate( prefac(scpc_set%mgx,scpc_set%mgy,scpc_set%mgz) )

            k0 = scpc_set%omega * 0.50_q * ryd2ev / dble( scpc_set%mxyz )
            k1 = fac1 / scpc_set%diel

            do ix=1,scpc_set%mgx
              do iy=1,scpc_set%mgy
                do iz=1,scpc_set%mgz
                  prefac(ix,iy,iz) = fac1 / scpc_set%meps(ix,iy,iz)
                end do
              end do
            end do

            raux  = scpc_set%mchg
            raux2 = 0.0_q
            e_old = 0.0_q
            e_new = 0.0_q

            do iter = 1, nmax
               e_old = e_new
               vaux = 0.0_q
               do ix=1,scpc_set%mgx
                 do iy=1,scpc_set%mgy
                   do iz=1,scpc_set%mgz
                     rho_tot(ix,iy,iz) = prefac(ix,iy,iz) * ( raux(ix,iy,iz) + raux2(ix,iy,iz) )
                   end do
                 end do
               end do
               call scpc_solve_isopot(rho_tot,vaux)
               do ix = 1, scpc_set%mgx
                 do iy = 1, scpc_set%mgy
                   call scpc_numderv('t',2,scpc_set%dmz,scpc_set%mgz, &
     &             vaux(ix,iy,1:scpc_set%mgz),dzvaux(ix,iy,1:scpc_set%mgz))
                   do iz = 1, scpc_set%mgz
                     tmp1 = scpc_set%mdeps(ix,iy,iz)*dzvaux(ix,iy,iz) / dabs(fac1)
                     tmp2 = (1.0_q-upfac)*raux2(ix,iy,iz)
                     raux2(ix,iy,iz) = upfac*tmp1 + tmp2
                   end do
                 end do
               end do
               e_new = sum( scpc_set%mchg * vaux ) * k0
               delta_e = abs( e_new - e_old )
               if ( abs(delta_e) <= 1.0d-4 ) exit
            end do

            if (allocated(raux2  )) deallocate(raux2  )
            if (allocated(dZraux2)) deallocate(dZraux2)
            if (allocated(dZvaux )) deallocate(dZvaux )
            if ( allocated( rho_tot ) ) deallocate( rho_tot )
            if ( allocated( prefac ) ) deallocate( prefac )

         end if

         mypot = vaux

         if (allocated(raux)) deallocate(raux)
         if (allocated(vaux)) deallocate(vaux)

      io_end

   end subroutine scpc_build_isopot

   subroutine scpc_solve_isopot(myrho, mypot)
#ifdef SCPC
      use pspfft
#endif
      real(q), intent(in) :: myrho(:,:,:)
      real(q), intent(out) :: mypot(:,:,:)
#ifdef SCPC
      type(arrayReal_3d_base) :: ar(3)
      character(ll) :: variablename(3)
      type(variableArrayReal_3d_groupBase) :: var(1)
      integer :: ntotalCells(3)
      real(kr) :: cellWidth(3)
      type(pspfft_form), pointer :: ps

      allocate(ar(1)%data(scpc_set%mgx, scpc_set%mgy, scpc_set%mgz))
      allocate(ar(2)%data(scpc_set%mgx, scpc_set%mgy, scpc_set%mgz))
      allocate(ar(3)%data(scpc_set%mgx, scpc_set%mgy, scpc_set%mgz))

      ar(1)%data = myrho
      ar(2)%data = 0.0_q
      ar(3)%data = 0.0_q
      cellWidth = [scpc_set%dmx, scpc_set%dmy, scpc_set%dmz]
      ntotalCells = [scpc_set%mgx, scpc_set%mgy, scpc_set%mgz]

      variablename(1) = 'Source'
      variablename(2) = 'AnalyticalSolution'
      variablename(3) = 'NumericalSolution'

      call initialize(var(1), ar, variableOption = variablename)

      var(1)%data(3)%data = var(1)%data(1)%data

      call create(ps, cellWidth, ntotalCells, MPI_comm_self)

      call solve(ps, var(1)%data(3:3))
      call destroy(ps)

      mypot = var(1)%data(3)%data

      call finalize(var(1))

      if (allocated(ar(1)%data)) deallocate(ar(1)%data)
      if (allocated(ar(2)%data)) deallocate(ar(2)%data)
      if (allocated(ar(3)%data)) deallocate(ar(3)%data)
#endif

   end subroutine scpc_solve_isopot

   subroutine scpc_mgsolver(gridc, myrho, mypot, mydiel, loc_name)
#ifdef SCPC
      use DL_mg
#endif
      type(grid_3d), intent(in) :: gridc
      real(q), intent(in) :: myrho(:,:,:)
      real(q), intent(inout) :: mypot(:,:,:)
      real(q), intent(in) :: mydiel(:,:,:)
      character(*), intent(in) :: loc_name
      integer, parameter :: mgnmax = 50
      real(q), parameter :: mgrtol = 1.0e-6_q
      real(q), parameter :: cte0 = -25.1327412287183_q * 0.50_q ! -fpi*two
      integer, parameter :: pgrid(3) = 1
      logical, parameter :: pperiods(3) = .false.
      real(q), allocatable :: eps_half(:,:,:,:)
!
      integer :: fd_order = 2
      integer :: mg_comm = 75498
      integer, save :: io_mg = 100
      integer :: my_out
!
      integer :: mpier = 0
      integer :: idx_s(3), idx_e(3)
      integer :: boundary_condition(3)
      character(len=80) :: fout
!
      integer :: viter(2) = 500
!
      integer :: ix, iy, iz
      real(q) :: tmp
!
#ifdef SCPC
      my_out = io_mg + gridc%comm%node_me
      write(fout,*) my_out
      fout = "MGSolver_" // trim(adjustl(loc_name)) // trim(adjustl(fout))

      call MPI_Cart_create(MPI_comm_self, 3, pgrid, pperiods, .false., mg_comm, mpier)

      idx_s = 1
      idx_e = [scpc_set%mgx, scpc_set%mgy, scpc_set%mgz]
      boundary_condition = dl_mg_bc_dirichlet

      call DL_mg_init(scpc_set%mgx, scpc_set%mgy, scpc_set%mgz, scpc_set%dmx, scpc_set%dmy, scpc_set%dmz, &
         boundary_condition, idx_s, idx_e, mg_comm, my_out, fout, errors_return=.false., ierror=mpier)

      allocate(eps_half(scpc_set%mgx, scpc_set%mgy, scpc_set%mgz, 3))
      eps_half = 1.0_q
      do ix = 1, scpc_set%mgx - 1
         do iy = 1, scpc_set%mgy - 1
            do iz = 1, scpc_set%mgz - 1
               tmp = mydiel(ix,iy,iz)
               eps_half(ix,iy,iz,1) = 0.50_q * (mydiel(ix+1,iy,iz) + tmp)
               eps_half(ix,iy,iz,2) = 0.50_q * (mydiel(ix,iy+1,iz) + tmp)
               eps_half(ix,iy,iz,3) = 0.50_q * (mydiel(ix,iy,iz+1) + tmp)
            end do
         end do
      end do

      call DL_mg_solver(mydiel, eps_half, cte0, myrho, mypot, fd_order, &
         tol_res_rel=mgrtol, tol_res_abs=mgrtol, tol_pot_rel=mgrtol, &
         tol_pot_abs=mgrtol, tol_vcycle_res_rel=mgrtol, use_damping=.false., &
         ierror=mpier,max_iters_defco=mgnmax, max_iters_vcycle=50, &
         max_iters_level_1=5000, v_iterations=viter)

      call DL_mg_free()

      io_mg = io_mg + 1

      if (allocated(eps_half)) deallocate(eps_half)
#endif
   end subroutine scpc_mgsolver

   subroutine vref_alignment(gridc, vinp, vout, emvout)
      type(grid_3d), intent(in) :: gridc
      real(q), intent(in) :: vinp(scpc_set%mgx, scpc_set%mgy, scpc_set%mgz)
      real(q), intent(out) :: vout(scpc_set%mgx, scpc_set%mgy, scpc_set%mgz), emvout
      real(q) :: dv(scpc_set%mgx,scpc_set%mgy,scpc_set%mgz)
      integer :: node_me, ionode, iout
      real(q) :: mdv
#ifdef MPI
      node_me = gridc%comm%node_me
      ionode  = gridc%comm%ionode
#endif
      io_begin
         emvout = 0.0_q
         dv = scpc_set%rpot - (scpc_set%vhar_loc + vinp)
         !
         iout = scpc_set%unitout
         open(iout, file=scpc_set%fout, position='append')
         mdv = sum(dv(1,:,:)) / dble(scpc_set%mgy * scpc_set%mgz)
         emvout = emvout + mdv
         !write(iout,*) "Pot alig ref averaged x:", mdv
         mdv = sum(dv(:,1,:)) / dble(scpc_set%mgx * scpc_set%mgz)
         emvout = emvout + mdv
         !write(iout,*) "Pot alig ref averaged y:", mdv
         mdv = sum(dv(:,:,1)) / dble(scpc_set%mgx * scpc_set%mgy)
         emvout = emvout + mdv
         !write(iout,*) "Pot alig ref averaged z:", mdv
         emvout = emvout / 3.0_q
         !write(iout,*) "Pot alig ref mean:", emvout
         close(iout)
         !
         vout = vinp + emvout
      io_end
   end subroutine vref_alignment

   subroutine scpc_potalignment(gridc, pot1, pot2)
      type(grid_3d), intent(in) :: gridc
      real(q), intent(in) :: pot1(:,:,:), pot2(:,:,:)
      integer :: nx, ny, nz, nxy, nyz, nxz, iout
      real(q) :: vmx, vmy, vmz, fac
      integer :: node_me, ionode

#ifdef MPI
      node_me = gridc%comm%node_me
      ionode  = gridc%comm%ionode
#endif

      nx = scpc_set%mgx
      ny = scpc_set%mgy
      nz = scpc_set%mgz

      nxy = nx * ny
      nxz = nx * nz
      nyz = ny * nz

      vmx = sum(pot1(1,:,:) ) - sum(pot2(1,:,:) )
      fac = sum(pot1(nx,:,:)) - sum(pot2(nx,:,:))
      vmx = (vmx + fac) * 0.5_q / dble(nyz)

      vmy = sum(pot1(:,1,:) ) - sum(pot2(:,1,:) )
      fac = sum(pot1(:,ny,:)) - sum(pot2(:,ny,:))
      vmy = (vmy + fac) * 0.5_q / dble(nxz)

      vmz = sum(pot1(:,:,1) ) - sum(pot2(:,:,1) )
      fac = sum(pot1(:,:,nz)) - sum(pot2(:,:,nz))
      vmz = (vmz + fac) * 0.5_q / dble(nxy)

      iout = scpc_set%unitout
      io_begin
         open(iout, file=scpc_set%fout, position='append')
         write(iout,'(/,2x,a,3f12.6,x,a)') "Potential Alignment (x,y,z):", vmx, vmy, vmz, " eV "
         close(iout)
      io_end

   end subroutine scpc_potalignment

   subroutine scpc_addpot(gridc, cvout)
      type(grid_3d), intent(in) :: gridc
      complex(q), intent(inout) :: cvout(gridc%mplwv)
      integer :: ix, iy, iz, k
      integer :: node_me, ionode
      integer :: nx, ny, nz, mx, my, mz
      real(q) :: dnx, dny, dnz, dmx, dmy, dmz
      real(q) :: rinpl
      real(q), allocatable :: nwork(:), mwork(:), rwork1(:), rwork2(:)
      real(q), allocatable :: rmat(:,:,:)
      complex(q), allocatable :: vadd(:)

#ifdef MPI
      node_me = gridc%comm%node_me
      ionode  = gridc%comm%ionode
#endif

      ! dense grid
      nx  = gridc%ngx
      ny  = gridc%ngy
      nz  = gridc%ngz
      dnx = scpc_set%dnx
      dny = scpc_set%dny
      dnz = scpc_set%dnz
      ! reduced grid
      mx  = scpc_set%mgx
      my  = scpc_set%mgy
      mz  = scpc_set%mgz
      dmx = scpc_set%dmx
      dmy = scpc_set%dmy
      dmz = scpc_set%dmz

      allocate(rwork1(DIMREAL(gridc%mplwv)))
      allocate(rwork2(gridc%ngx * gridc%ngy))
      allocate(vadd(gridc%mplwv))

      io_begin

         allocate(nwork(gridc%ngx * gridc%ngy * gridc%ngz))
         allocate(mwork(scpc_set%mgx * scpc_set%mgy * scpc_set%mgz))
         allocate(rmat(gridc%ngx, gridc%ngy, gridc%ngz))

         k = 1
         do iz = 1, mz
            do iy = 1, my
               do ix = 1, mx
                  mwork(k) = scpc_set%vadd(ix,iy,iz)
                  k = k + 1
                end do
            end do
         end do
         nwork = 0.0_q
         call scpc_mltscale(mwork, mx, my, mz, dmx, dmy, dmz, &
                            nwork, nx, ny, nz, dnx, dny, dnz)
         deallocate(mwork)

         k = 1
         do iz = 1, nz
            do iy = 1, ny
               do ix = 1, nx
                  rmat(ix,iy,iz) = nwork(k)
                  k = k + 1
               end do
            end do
         end do
         deallocate(nwork)

      io_end
      CALLMPI(M_barrier(gridc%comm))

      do iz=1,nz
         io_begin
            k = 1
            do iy=1,ny
               do ix=1,nx
                  rwork2(k) = rmat(ix,iy,iz)
                  k = k + 1
               end do
            end do
         io_end
         call dis_grid_rl_plane(gridc, rwork2, rwork1, .true., iz)
      end do
      CALLMPI(M_barrier(gridc%comm))

      rinpl = 1.0_q / dble(gridc%nplwv)
      call rl_add(rwork1, rinpl, rwork1, 0.0_q, vadd, gridc)
      call FFT3D(vadd, gridc, -1)

      do k = 1, gridc%mplwv
         cvout(k) = cvout(k) + vadd(k)
      end do

      deallocate(rwork1)
      deallocate(rwork2)
      deallocate(vadd)
      io_begin
         deallocate(rmat)
      io_end
      CALLMPI(M_barrier(gridc%comm))

   end subroutine scpc_addpot

   subroutine scpc_prtall(gridc)
      type(grid_3d), intent(in) :: gridc
      if(scpc_set%prtx) then
         call scpc_averages(gridc, 'X', 'x-mrho.dat', scpc_set%mchg, scpc_set%dmx, 1.0_q)
         call scpc_averages(gridc, 'X', 'x-mpot.dat', scpc_set%mpot, scpc_set%dmx, 1.0_q)
      end if
      if(scpc_set%prty) then
         call scpc_averages(gridc, 'Y', 'y-mrho.dat', scpc_set%mchg, scpc_set%dmy, 1.0_q)
         call scpc_averages(gridc, 'Y', 'y-mpot.dat', scpc_set%mpot, scpc_set%dmy, 1.0_q)
      end if
      if(scpc_set%prtz) then
         call scpc_averages(gridc, 'Z', 'z-mrho.dat', scpc_set%mchg, scpc_set%dmz, 1.0_q)
         call scpc_averages(gridc, 'Z', 'z-mpot.dat', scpc_set%mpot, scpc_set%dmz, 1.0_q)
      end if
      CALLMPI(M_barrier(gridc%comm))
   end subroutine scpc_prtall

   subroutine scpc_averages(gridc, axis, fio, mat, ds, fac)
      type(grid_3d), intent(in) :: gridc
      character(len=1), intent(in) :: axis
      character(*), intent(in) :: fio
      real(q), intent(in) :: mat(:,:,:)
      real(q), intent(in) :: ds, fac
      integer, parameter :: avgio = 5656
      integer :: ii, ionode, node_me
      real(q) :: s, c0, tmp
      integer :: nx, ny, nz

#ifdef MPI
      ionode  = gridc%comm%ionode
      node_me = gridc%comm%node_me
#endif

      io_begin
         nx = size(mat,1)
         ny = size(mat,2)
         nz = size(mat,3)
         open(avgio, file=fio, position='append')
         write(avgio,*) "#", fio
         s = 0.0_q
         select case (axis)
         case ('x','X')
            c0 = fac / dble(ny * nz)
            do ii = 1, nx
               tmp = sum(mat(ii,:,:)) * c0
               write(avgio,'(2e20.10)') s, tmp
               s = s + ds
            end do
         case ('y','Y')
            c0 = fac / dble(nx * nz)
            do ii = 1, ny
               tmp = sum(mat(:,ii,:)) * c0
               write(avgio,'(2e20.10)') s, tmp
               s = s + ds
            end do
         case ('z','Z')
            c0 = fac / dble(nx * ny)
            do ii = 1, nz
               tmp = sum(mat(:,:,ii)) * c0
               write(avgio,'(2e20.10)') s, tmp
               s = s + ds
            end do
         end select
         write(avgio,*) ""
         close(avgio)
      io_end

   end subroutine scpc_averages

end module scpc
