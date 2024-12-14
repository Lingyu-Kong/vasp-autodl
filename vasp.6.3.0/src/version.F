#include "symbol.inc"
module version

   implicit none

   integer, parameter :: major = 6, minor = 3, patch = 0
   character(len=*), parameter :: date = '20Jan22'

contains

   pure character(len=80) function vasp()
      use string, only: str
      character(len=:), allocatable :: label, build, version_
#ifdef VASP_VERSION
      label = VASP_VERSION
#else
      label = 'vasp.' // str(major) // '.' // str(minor) // '.' // str(patch) // ' ' // date
#endif
      build = '(build ' // __DATE__ // ' ' // __TIME__ // ')'
#ifdef gammareal
      version_ = 'gamma-only'
#else
      version_ = 'complex'
#endif
      vasp = label // ' ' // build // ' ' // version_
   end function vasp

end module version
