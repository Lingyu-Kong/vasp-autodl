#include "symbol.inc"
!************************************************************************
!
!  this module contains some control data structures for
!  machine-learning force field.
!
!***********************************************************************
      MODULE ML_FF_CONSTANT
      USE ML_FF_PREC
! Constants used for machine learning force field

!  Some important Parameters, to convert to a.u.
!  - AUTOA  = 1. a.u. in Angstroem
!  - RYTOEV = 1 Ry in Ev
!  - EVTOJ  = 1 eV in Joule
!  - AMTOKG = 1 atomic mass unit ("proton mass") in kg

      REAL(q), PARAMETER :: AUTOA=0.529177249_q,RYTOEV=13.605826_q
      REAL(q), PARAMETER :: EVTOJ=1.60217733E-19_q,AMTOKG=1.6605402E-27_q
      REAL(q),PARAMETER  :: PI =3.141592653589793238_q

! EUNIT = 1 Hartree in ev units
! MUNIT = Electron mass in atomic units
! RUNIT = Bohr Radius in Angstrom1
! FUNIT = Unit of force in Hartree/Bohr
! SUNIT = Unit of stress tensor
! TUNIT = Unit of time in atomic unit

      REAL(q), PARAMETER :: EUNIT = 2.0_q*RYTOEV, &
                            MUNIT = 1.054571726E-034_q**2/2.0_q/(RYTOEV*EVTOJ*AUTOA**2*1.0E-20_q)/AMTOKG, &
                            FUNIT = 2.0_q*RYTOEV/AUTOA, &
                            SUNIT = 2.0_q*RYTOEV/(AUTOA*1.0E-10_q)**3*EVTOJ/1.0E8_q, &
                            TUNIT = 1.054571726E-034_q/(2.0_q*RYTOEV*EVTOJ)*1.0E15_q

! Parameters for Normalization (to avoid it if norm becomes zero)

      REAL(q),PARAMETER :: TOLERANCE_NORM = 1.0E-08_q

      INTEGER, PARAMETER :: NBLOCK_SCALAPACK = 32

      END MODULE
