#include "symbol.inc"
!************************************************************************
! RCS:  $Id: base.F,v 1.2 2001/02/20 14:44:56 kresse Exp $
!
! this module contains some control data structures for VASP
!
!***********************************************************************
      MODULE ML_FF_PREC
      INTEGER, PARAMETER :: q =SELECTED_REAL_KIND(10)
#ifndef noQuadPrecision
      INTEGER, PARAMETER :: qd=SELECTED_REAL_KIND(30)
#else
      INTEGER, PARAMETER :: qd=SELECTED_REAL_KIND(10)
#endif
      INTEGER, PARAMETER :: qs=SELECTED_REAL_KIND(5)
!
! this parameter controls the step width in for numerical
! differentiation
! in some VASP routines
! 1E-5 is very reliable yielding at least 7 digits in all contributions
! to the forces
! 1E-4, however, is better suited for second derivatives
! for reasons of consistentcy with previous versions 1E-5
      REAL(q), PARAMETER :: fd_displacement=1E-5
      END MODULE ML_FF_PREC
