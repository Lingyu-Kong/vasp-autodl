#include "symbol.inc"
!****************** module use external potential  *********************
! This is a dummy module for the Yu/Libisch/Carter PAW DFET embedding
!***********************************************************************

  module mextpot_glb
    use prec
    implicit none
    
    type extpot
      logical :: lextpot, dextpot
      integer :: nx, ny, nz ! should be equal to the dimension of potential grid
      real(q) :: FCORR
    endtype extpot

    type (extpot) EXTPT
    RGRID, allocatable :: EXTPT_GRID(:,:,:) ! the external potential grid in normal layout
    RGRID, allocatable :: VLM(:,:,:)        ! external potential for each channel of each ion
    RGRID, allocatable :: SCX(:), SCY(:), SCZ(:)
    RGRID, allocatable :: SGRID(:,:,:)      ! the B-spline interpolation coefficients
    ! parameters for real space VLM projection
    integer, parameter :: KMAX=20, di_smear=5, NGRID=50
    real(q), parameter :: Delta = 0.1_q
  end module mextpot_glb

  module mextpot
   use mextpot_glb
   implicit none

   private
   public :: extpt_reader,extpt_lextpot,extpt_dextpot

   contains

   ! read INCAR FILE for external potential information
   subroutine EXTPT_READER(NGXF,NGYF,NGZF,IU0,IU5,IU6)
      use base
      use vaspxml
      implicit none

      integer :: NGXF, NGYF, NGZF, IU0, IU5, IU6

      EXTPT%lextpot=.FALSE.
      EXTPT%dextpot=.FALSE.
   end subroutine EXTPT_READER

   function EXTPT_LEXTPOT()
      logical :: EXTPT_LEXTPOT
      EXTPT_LEXTPOT=EXTPT%lextpot
   end function EXTPT_LEXTPOT

   function EXTPT_DEXTPOT()
      logical :: EXTPT_DEXTPOT
      EXTPT_DEXTPOT=EXTPT%dextpot
   end function EXTPT_DEXTPOT
   
  end module mextpot
