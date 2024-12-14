#include "symbol.inc"

!***********************************************************************
! 
! MODULE: minimax_functions
!
!> @author
!> Merzuk Kaltak, University of Vienna
! 
! DESCRIPTION:
!> contains all fitting and error functions used in minimax module
!>
!***********************************************************************

MODULE minimax_functions2D
   USE prec
   USE minimax_functions1D    
#ifdef qd_emulate
   USE qdmodule
#endif 

   CONTAINS 

!***********************************************************************
! DESCRIPTION:
!>calculates basis function:
!> \f$ 
!>     \phi(x,\lambda) = e^{-z \lambda}
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION EXPF( Z, L )
      QDPREAL               :: EXPF
      QDPREAL , INTENT(IN)  :: Z ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      PROFILING_START('expf')
      
DOESI IF (-Z*L > 100.0_q) THEN; EXPF = 1.0E10_q; PROFILING_STOP('expf'); RETURN; END IF
      EXPF = EXP( -Z*L )
      PROFILING_STOP('expf')
   END FUNCTION EXPF

!***********************************************************************
! DESCRIPTION:
!>calculates derivative of basis function w.r.t. first argument
!> \f$ 
!>     \phi(z,\lambda) = e^{-z \lambda}
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION D_EXPF_DZ( Z, L )
      QDPREAL               :: D_EXPF_DZ
      QDPREAL , INTENT(IN)  :: Z ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      PROFILING_START('d_expf_dz')
      
      D_EXPF_DZ = -L*EXPF( Z,L )

      PROFILING_STOP('d_expf_dz')
   END FUNCTION D_EXPF_DZ

!***********************************************************************
! DESCRIPTION:
!>calculates first derivative of basis function w.r.t. second argument
!> \f$ 
!>     \phi(z,\lambda) = e^{-z \lambda}
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION D_EXPF_DL( Z, L )
      QDPREAL               :: D_EXPF_DL
      QDPREAL , INTENT(IN)  :: Z ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      PROFILING_START('d_expf_dl')
      
      D_EXPF_DL = -Z*EXPF( Z,L )

      PROFILING_STOP('d_expf_dl')
   END FUNCTION D_EXPF_DL

!***********************************************************************
! DESCRIPTION:
!>calculates second derivative of basis function w.r.t. first argument
!> \f$ 
!>     \phi(z,\lambda) = e^{-z \lambda}
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION D2_EXPF_DZ2( Z, L )
      QDPREAL               :: D2_EXPF_DZ2
      QDPREAL , INTENT(IN)  :: Z ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      PROFILING_START('d2_expf_dz2')
      
      D2_EXPF_DZ2 = -L*D_EXPF_DZ( Z,L )

      PROFILING_STOP('d2_expf_dz2')
   END FUNCTION D2_EXPF_DZ2

!***********************************************************************
! DESCRIPTION:
!>calculates basis function:
!> \f$ 
!>     \phi(z,\lambda) = \frac{ z }{ z^2 + \lambda^2 }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION UFREQ( Z, L )
      QDPREAL              :: UFREQ
      QDPREAL, INTENT(IN)  :: Z ! sampling value
      QDPREAL, INTENT(IN)  :: L ! coefficient defining function
      PROFILING_START('ufreq')
      
      UFREQ = Z/( Z*Z+L*L )

      PROFILING_STOP('ufreq')
   END FUNCTION UFREQ 

!***********************************************************************
! DESCRIPTION:
!>calculates derivative of basis function w.r.t. first argument
!> \f$ 
!>     \phi(z,\lambda) = \frac{ z }{ z^2 + \lambda^2 }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION D_UFREQ_DZ( Z, L )
      QDPREAL               :: D_UFREQ_DZ
      QDPREAL , INTENT(IN)  :: Z ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      PROFILING_START('d_ufreq_dz')
      
      D_UFREQ_DZ = VFREQ(Z,L)**2 - UFREQ(Z,L)**2

      PROFILING_STOP('d_ufreq_dz')
   END FUNCTION D_UFREQ_DZ

!***********************************************************************
! DESCRIPTION:
!>calculates derivative of basis function w.r.t. second argument
!> \f$ 
!>     \phi(z,\lambda) = \frac{ z }{ z^2 + \lambda^2 }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION D_UFREQ_DL( Z, L )
      QDPREAL               :: D_UFREQ_DL
      QDPREAL , INTENT(IN)  :: Z ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      PROFILING_START('d_ufreq_dl')
      
      D_UFREQ_DL = -QDPP(2)*VFREQ(Z,L)*UFREQ(Z,L)

      PROFILING_STOP('d_ufreq_dl')
   END FUNCTION D_UFREQ_DL

!***********************************************************************
! DESCRIPTION:
!>calculates the second two derivatives of this function:
!> \f$ 
!>     \frac{ x }{ x^2 + \lambda^2 }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D2_UFREQ_DZ2( X, L )
      QDPREAL               :: D2_UFREQ_DZ2
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      
      PROFILING_START('d2_ufreq_dz2')

      D2_UFREQ_DZ2 = QDPP(2)*UFREQ(X,L)*&
         (UFREQ(X,L)**2 -QDPP(3)*VFREQ(X,L)**2 )

      PROFILING_STOP('d2_ufreq_dz2')
   END FUNCTION D2_UFREQ_DZ2

!***********************************************************************
! DESCRIPTION:
!>calculates basis function:
!> \f$ 
!>     \phi(x,\lambda)=\frac{ \lambda }{ x^2 + \lambda^2 }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION VFREQ( X, L )
      QDPREAL               :: VFREQ
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      
      PROFILING_START( 'vfreq' )

      VFREQ = L/( X*X+L*L )

      PROFILING_STOP( 'vfreq' )
   END FUNCTION VFREQ 

!***********************************************************************
! DESCRIPTION:
!>calculates the first derivative of this function (w.r.t. to x)
!> \f$ 
!>     \frac{ \lambda }{ x^2 + \lambda^2 }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D_VFREQ_DZ( X, L )
      QDPREAL               :: D_VFREQ_DZ
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd_vfreq_dz' )
      
      D_VFREQ_DZ = -QDPP(2)*UFREQ(X,L)*VFREQ(X,L)
      
      PROFILING_STOP( 'd_vfreq_dz' )

   END FUNCTION D_VFREQ_DZ

!***********************************************************************
! DESCRIPTION:
!>calculates the first derivative of this function (w.r.t. to L)
!> \f$ 
!>     \frac{ \lambda }{ x^2 + \lambda^2 }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D_VFREQ_DL( Z, L )
      QDPREAL               :: D_VFREQ_DL
      QDPREAL , INTENT(IN)  :: Z ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd_vfreq_dl' )
      
      D_VFREQ_DL = UFREQ(Z,L)**2 - VFREQ(Z,L)**2
      
      PROFILING_STOP( 'd_vfreq_dl' )

   END FUNCTION D_VFREQ_DL

!***********************************************************************
! DESCRIPTION:
!>calculates the secondderivative of this function (w.r.t. to x)
!> \f$ 
!>     \frac{ \lambda }{ x^2 + \lambda^2 }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D2_VFREQ_DZ2( X, L )
      QDPREAL               :: D2_VFREQ_DZ2
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd2_vfreq_dz2' )
      
      D2_VFREQ_DZ2 = -QDPP(2)*( D_UFREQ_DZ(X,L)*VFREQ(X,L) &
                              + D_VFREQ_DZ(X,L)*UFREQ(X,L) )
      
      PROFILING_STOP( 'd2_vfreq_dz2' )

   END FUNCTION D2_VFREQ_DZ2

!***********************************************************************
! DESCRIPTION:
!>calculates basis function for time domain :
!> \f$ 
!>     \frac{ cosh(x/2(1-2*L)) }{ cosh(x/2) }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION UTIME( X, L )
      QDPREAL               :: UTIME
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      PROFILING_START( 'utime' )
      
      UTIME = COSHCOSH( X/QDPP(2), QDPP(1)-QDPP(2)*L )

      PROFILING_STOP( 'utime' )
   END FUNCTION UTIME

!***********************************************************************
! DESCRIPTION:
!>calculates odd function for time domain :
!> \f$ 
!>     \frac{ sinh(x/2(1-2*L)) }{ cosh(x/2) }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION VTIME( X, L )
      QDPREAL               :: VTIME
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      PROFILING_START( 'vtime' )
      
      VTIME = SINHCOSH( X/QDPP(2), QDPP(1)-QDPP(2)*L )

      PROFILING_STOP( 'vtime' )
   END FUNCTION VTIME

!***********************************************************************
! DESCRIPTION:
!>calculates first derivative of follwing function (w.r.t. to x ) 
!> \f$ 
!>     \frac{ cosh(z/2(1-2*L)) }{ cosh(z/2) }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION D_UTIME_DZ( Z, L )
      QDPREAL               :: D_UTIME_DZ
      QDPREAL , INTENT(IN)  :: Z ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      PROFILING_START( 'd_utime_dz' )
      
      D_UTIME_DZ = DCOSHCOSHX( Z/QDPP(2), QDPP(1)-QDPP(2)*L )/QDPP(2)

      PROFILING_STOP( 'd_utime_dz' )
   END FUNCTION D_UTIME_DZ

!***********************************************************************
! DESCRIPTION:
!>calculates first derivative of follwing function (w.r.t. to L ) 
!> \f$ 
!>     \frac{ cosh(z/2(1-2*L)) }{ cosh(z/2) }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION D_UTIME_DL( Z, L )
      QDPREAL               :: D_UTIME_DL
      QDPREAL , INTENT(IN)  :: Z ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      PROFILING_START( 'd_utime_dl' )
      
      D_UTIME_DL = -QDPP(2)*DCOSHCOSHY( Z/QDPP(2), QDPP(1)-QDPP(2)*L )

      PROFILING_STOP( 'd_utime_dl' )
   END FUNCTION D_UTIME_DL

!***********************************************************************
! DESCRIPTION:
!>calculates second derivative of follwing function (w.r.t. to x ) 
!> \f$ 
!>     \frac{ cosh(z/2(1-2*L)) }{ cosh(z/2) }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  Z  (argument of error function)
!
!***********************************************************************

   FUNCTION D2_UTIME_DZ2( Z, L )
      QDPREAL               :: D2_UTIME_DZ2
      QDPREAL , INTENT(IN)  :: Z ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      PROFILING_START( 'd2_utime_dz2' )
      
      D2_UTIME_DZ2 = DCOSHCOSHX2( Z/QDPP(2), QDPP(1)-QDPP(2)*L )/QDPP(4)

      PROFILING_STOP( 'd2_utime_dz2' )
   END FUNCTION D2_UTIME_DZ2

!***********************************************************************
! DESCRIPTION:
!>calculates basis function for bosonic Matsubara freuencies:
!> \f$ 
!>     \frac{ x tanh(x/2) }{ x^2+L^2 }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION UTANH( X, L )
      QDPREAL               :: UTANH
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      PROFILING_START( 'utanh' )
      
      IF ( L == 0 ) THEN
         UTANH = F_THZ_OVER_Z( X ) 
      ELSE
         ! taylor expansion for small arguments
         IF ( ABS( X ) < QDPP( 1E-8 )  ) THEN
            UTANH = QDPP(0.5)*(X/L)**2                 & 
                          - (QDPP(12)+L**2)/QDPP(24)*(X/L)**4  &
                          + (QDPP(120)+QDPP(10)*L**2+L**4)/QDPP(240)*(X/L)**4
         ELSE
            UTANH = TH( X/QDPP(2) ) * UFREQ( X, L ) 
         ENDIF
      ENDIF
      PROFILING_STOP( 'utanh' )
   END FUNCTION UTANH 

!***********************************************************************
! DESCRIPTION:
!>calculates first derivative of basis function w.r.t. x :
!> \f$ 
!>     \frac{ x tanh(x/2) }{ x^2+L^2 }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D_UTANH_DZ( X, L )
      QDPREAL               :: D_UTANH_DZ
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      ! local 
      QDPREAL               :: T,C,CP
      PROFILING_START('d_utanh_dz')
    
      T = TH( X/2 )
      C  = UFREQ(X,L)
      CP = D_UFREQ_DZ(X,L)

      D_UTANH_DZ = (QDPP(1)-T**2)*C/QDPP(2) + T*CP

      PROFILING_STOP('d_utanh_dz')
   END FUNCTION D_UTANH_DZ 

!***********************************************************************
! DESCRIPTION:
!>calculates first derivative of basis function w.r.t. L :
!> \f$ 
!>     \frac{ x tanh(x/2) }{ x^2+L^2 }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D_UTANH_DL( X, L )
      QDPREAL               :: D_UTANH_DL
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      ! local 
      PROFILING_START('d_utanh_dl')
    
      D_UTANH_DL = D_UFREQ_DL(X,L)*TH(X/2)

      PROFILING_STOP('d_utanh_dl')
   END FUNCTION D_UTANH_DL 

!***********************************************************************
! DESCRIPTION:
!>calculates second derivative of basis function w.r.t. x :
!> \f$ 
!>     \frac{ x tanh(x/2) }{ x^2+L^2 }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D2_UTANH_DZ2( X, L )
      QDPREAL               :: D2_UTANH_DZ2
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      ! local 
      QDPREAL               :: T,C,CP(2)
      PROFILING_START('d2_utanh_dz2')
    
      T = TH( X/2 )
      C  = UFREQ(X,L)
      CP(1) = D_UFREQ_DZ(X,L)
      CP(2) = D2_UFREQ_DZ2(X,L)

      D2_UTANH_DZ2 = (QDPP(1)-T**2)*CP(1) + T*CP(2) &
                           - T*(QDPP(1)-T**2)*C/QDPP(2)

      PROFILING_STOP('d2_utanh_dz2')
   END FUNCTION D2_UTANH_DZ2 

!***********************************************************************
! DESCRIPTION:
!>calculates basis function for bosonic Matsubara freuencies:
!> \f$ 
!>     \frac{ L \tanh x/2}{ x^2+L^2 }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION VTANH( X, L )
      QDPREAL               :: VTANH
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      
      VTANH = VFREQ( X, L )*TH(X/2) 

   END FUNCTION VTANH 

!***********************************************************************
! DESCRIPTION:
!>calculates frist derivative of this basis function w.r.t. to X
!> \f$ 
!>     \frac{ L \tanh x/2}{ x^2+L^2 }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D_VTANH_DZ( X, L )
      QDPREAL               :: D_VTANH_DZ
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      QDPREAL               :: T
      PROFILING_START( 'd_vtanh_dz' )
      T = TH(X/2)

      D_VTANH_DZ = D_VFREQ_DZ( X, L )*T + VFREQ(X,L)*(QDPP(1)-T**2)/QDPP(2)

      PROFILING_STOP( 'd_vtanh_dz' )
   END FUNCTION D_VTANH_DZ

!***********************************************************************
! DESCRIPTION:
!>calculates frist derivative of this basis function w.r.t. to L
!> \f$ 
!>     \frac{ L \tanh x/2}{ x^2+L^2 }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D_VTANH_DL( X, L )
      QDPREAL               :: D_VTANH_DL
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      PROFILING_START( 'd_vtanh_dl' )

      D_VTANH_DL = D_VFREQ_DL( X, L )*TH(X/2)

      PROFILING_STOP( 'd_vtanh_dl' )
   END FUNCTION D_VTANH_DL

!***********************************************************************
! DESCRIPTION:
!>calculates frist derivative of this basis function w.r.t. to X
!> \f$ 
!>     \frac{ L \tanh x/2}{ x^2+L^2 }
!> \f$ 
!
!> @param[in]  L  (parameter \f$\lambda\f$)
!> @param[in]  X  (argument of error function)
!
!***********************************************************************

   FUNCTION D2_VTANH_DZ2( X, L )
      QDPREAL               :: D2_VTANH_DZ2
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function
      QDPREAL               :: T,DT
      PROFILING_START( 'd2_vtanh_dz2' )
      T = TH(X/2)
      DT = (QDPP(1)-T**2)

      D2_VTANH_DZ2 = T*D2_VFREQ_DZ2(X,L) &
                   + DT*( D_VFREQ_DZ(X,L)-T*VFREQ(X,L)/QDPP(2) )

      PROFILING_STOP( 'd2_vtanh_dz2' )
   END FUNCTION D2_VTANH_DZ2

!****************************************************************************
!
! DESCRIPTION:
!>calculates scalar product of two vectors in quad-precision:
!  
!> @param[in]  A  (first vector)
!> @param[in]  B  (second vector)
! 
!****************************************************************************

FUNCTION DOT_PRODUCT_MPR( A , B ) 
   USE prec
#ifdef qd_emulate
   USE qdmodule
#endif
   QDPREAL :: A(:), B(:)
   QDPREAL :: DOT_PRODUCT_MPR, P
   INTEGER :: I, N

   PROFILING_START('dot_product_mpr')

   P=0
   N = SIZE( A ) 
   IF ( N == SIZE( B ) ) THEN 
#ifdef qd_emulate
      DO I = 1, N
          P = P + ( A(I)*B(I) )
      ENDDO
      DOT_PRODUCT_MPR = P 
#else
      DOT_PRODUCT_MPR = DOT_PRODUCT( a, b ) 
#endif
   ELSE
      WRITE( *, * )&
      ' INTERNAL ERROR IN VASP: DOT_PRODUCT_MPR REPORTS INCONSISTENT SIZES',&
      N, SIZE( B ) 
   ENDIF 

   PROFILING_STOP('dot_product_mpr')

ENDFUNCTION DOT_PRODUCT_MPR 

!****************************************************************************
! DESCRIPTION:
!>calculates following function
!> \f$ 
!>   \frac{cosh( x y)}{cosh(x)}
!> \f$ 
!  
!> @param[in]  X  (argument)
!> @param[in]  Y  (argument)
!****************************************************************************

   FUNCTION COSHCOSH( X, Y ) 
      USE prec
#ifdef qd_emulate
      USE qdmodule
#endif
      QDPREAL :: X, Y
      QDPREAL :: COSHCOSH
      ! local 
      QDPREAL :: AX, AXY, D
   
      PROFILING_START( 'coshcosh' )
    
      ! This function is a beast!
      ! For finite temperature quadratures a naive evaluation will
      ! certainly produce an over or underflow or even a NaN 
      AX = ABS( X )
      AXY = ABS( X*Y )

      ! if both arguments are larger than 1000 : 
      IF ( AX > QDPD(MAXE)/QDPD(2) .AND. AXY > QDPD(MAXE)/QDPD(2) ) THEN
         ! here the difference of arguments is important
         D = X*Y - X 
         COSHCOSH = EXP( SIGN_QD(D)*MIN( ABS(D), QDPD(MAXE)/QDPD(2) ) ) 
      ELSE
         ! in case denominator is larger than nummerator, 
         ! prevent overflow in nummerator
         IF ( AX - AXY > QDPD(MAXE) ) THEN
            COSHCOSH = EXP( -MIN( ABS( AX - AXY ) , QDPD(MAXE) ) )
         ELSE IF ( AXY - AX > QDPD(MAXE) ) THEN
            COSHCOSH = EXP( MIN( ABS( AX - AXY ) , QDPD(MAXE) ) )
         ELSE
            COSHCOSH = COSH(AXY)/COSH(AX)
         ENDIF
#ifdef debug
         ! make sure program can deal with this
         D = EXPONENT( COSHCOSH )
         IF ( ABS(D) > QDPD(MAXE) ) THEN
            WRITE(*,*)'over or underflow in COSHCOSH',AXY,AX,COSHCOSH,EXP(SIGN_QD(D)*QDPD(MAXE))
            COSHCOSH = EXP( SIGN_QD(D)*QDPD(MAXE) )
         ENDIF
#endif
      ENDIF

      PROFILING_STOP( 'coshcosh' )
   END FUNCTION COSHCOSH

!****************************************************************************
! DESCRIPTION:
!>calculates following function
!> \f$ 
!>   \frac{sinh( x y)}{cosh(x)}
!> \f$ 
!  
!> @param[in]  X  (argument)
!> @param[in]  Y  (argument)
!****************************************************************************

   FUNCTION SINHCOSH( X, Y ) 
      USE prec
#ifdef qd_emulate
      USE qdmodule
#endif
      QDPREAL :: X, Y
      QDPREAL :: SINHCOSH
      ! local 
      QDPREAL :: AX, AXY, D
   
      PROFILING_START( 'sinhcosh' )
    
      ! This function is probably also a beast!
      AX = ABS( X )
      AXY = ABS( X*Y )

      ! if both arguments are larger than 1000 : 
      IF ( AX > QDPD(MAXE)/QDPD(2) .AND. AXY > QDPD(MAXE)/QDPP(2) ) THEN
         ! here the difference of arguments is important
         D = X*Y - X 
         SINHCOSH = SIGN_QD(X*Y)*EXP( SIGN_QD(D)*MIN( ABS(D), QDPD(MAXE)/QDPP(2) ) ) 
      ELSE
         ! in case denominator is larger than nummerator, 
         ! prevent overflow in nummerator
         IF ( AX - AXY > QDPD(MAXE) ) THEN
            SINHCOSH = SIGN_QD(X*Y)*EXP( -MIN( ABS( AX - AXY ) , QDPD(MAXE) ) )
         ELSE IF ( AXY - AX > QDPD(MAXE) ) THEN
            SINHCOSH = SIGN_QD(X*Y)*EXP( MIN( ABS( AX - AXY ) , QDPD(MAXE) ) )
         ELSE
            SINHCOSH = SINH(AXY)/COSH(AX)
         ENDIF
#ifdef debug
         ! make sure program can deal with this
         D = EXPONENT( SINHCOSH )*LOG(2._q)
         IF ( ABS(D) > QDPD(MAXE) ) THEN
            WRITE(*,*)'over or underflow in COSHCOSH',AXY,AX,SINHCOSH,EXP(SIGN_QD(D)*QDPD(MAXE))
            WRITE(*,*)AXY, EXPONENT( CXY), AX, EXPONENT( CX)
            SINHCOSH = EXP( SIGN_QD(D)*QDPD(MAXE) )
         ENDIF
#endif
      ENDIF
   
      PROFILING_STOP( 'sinhcosh' )
   END FUNCTION SINHCOSH

!****************************************************************************
! DESCRIPTION:
!>calculates following derivative
!> \f$ 
!>   \frac{\partial}{\partial x} \frac{cosh( x y)}{cosh(x)}
!>  =y\frac{sinh( x y )}{cosh( x )} - \frac{cosh( x y)}{cosh(x)} tanh(x)
!> \f$ 
!  
!> @param[in]  X  (argument)
!> @param[in]  Y  (argument)
!****************************************************************************

   FUNCTION DCOSHCOSHX( X, Y ) 
      USE prec
#ifdef qd_emulate
      USE qdmodule
#endif
      QDPREAL :: X, Y 
      QDPREAL :: DCOSHCOSHX
      PROFILING_START( 'dcoshcosx' )
   
      DCOSHCOSHX = Y * SINHCOSH( X, Y ) - COSHCOSH( X, Y )*TH( X )
   
      PROFILING_STOP( 'dcoshcosx' )
   END FUNCTION DCOSHCOSHX

!****************************************************************************
! DESCRIPTION:
!>calculates following derivative
!> \f$ 
!>   \frac{\partial^2}{\partial x^2} \frac{cosh( x y)}{cosh(x)}
!> \f$ 
!  
!> @param[in]  X  (argument)
!> @param[in]  Y  (argument)
!****************************************************************************

   FUNCTION DCOSHCOSHX2( X, Y ) 
      USE prec
#ifdef qd_emulate
      USE qdmodule
#endif
      QDPREAL :: X, Y 
      QDPREAL :: DCOSHCOSHX2
   
      PROFILING_START( 'dcoshcosx2' )
      DCOSHCOSHX2 = COSHCOSH( X, Y )*( 2*TH( X )**2 + Y**2 - QDPP(1) ) &
                  - 2*Y*SINHCOSH( X, Y )*TH( X )
   
      PROFILING_STOP( 'dcoshcosx2' )
   END FUNCTION DCOSHCOSHX2

!****************************************************************************
! DESCRIPTION:
!>calculates following function
!> \f$ 
!>   \frac{cosh( x y)^2}{cosh(x)^2}
!> \f$ 
!  
!> @param[in]  X  (argument)
!> @param[in]  Y  (argument)
!****************************************************************************

   FUNCTION COSH2COSH( X, Y ) 
      QDPREAL :: X, Y
      QDPREAL :: COSH2COSH
   
      COSH2COSH = COSHCOSH(X,Y)**2

   END FUNCTION COSH2COSH

!****************************************************************************
! DESCRIPTION:
!>calculates following derivative
!> \f$ 
!>   \frac{\partial}{\partial x} \frac{cosh( x y)^2}{cosh(x)^2}
!> \f$ 
!  
!> @param[in]  X  (argument)
!> @param[in]  Y  (argument)
!****************************************************************************

   FUNCTION DCOSH2COSHX( X, Y ) 
      USE prec
#ifdef qd_emulate
      USE qdmodule
#endif
      QDPREAL :: X, Y 
      QDPREAL :: DCOSH2COSHX
   
      DCOSH2COSHX = QDPP(2)*COSH2COSH(X,Y)*( Y * TH( X*Y ) - TH( X ) ) 
   
   END FUNCTION DCOSH2COSHX

!****************************************************************************
! DESCRIPTION:
!>calculates following derivative
!> \f$ 
!>   \frac{\partial^2}{\partial x^2} \frac{cosh( x y)^2}{cosh(x)^2}
!> \f$ 
!  
!> @param[in]  X  (argument)
!> @param[in]  Y  (argument)
!****************************************************************************

   FUNCTION D2COSH2COSHX2( X, Y ) 
      USE prec
#ifdef qd_emulate
      USE qdmodule
#endif
      QDPREAL :: X, Y 
      QDPREAL :: D2COSH2COSHX2
   
      D2COSH2COSHX2 = QDPP(2)*COSH2COSH(X,Y)*( &
                      Y**2 * TH( X*Y )**2      &
                    - QDPP(4)*Y*TH( X )*TH( X*Y )  &
                    + QDPP(3)*TH( X )**2        &
                    + Y**2-QDPP(1)  ) 
   
   END FUNCTION D2COSH2COSHX2

!****************************************************************************
! DESCRIPTION:
!>calculates following derivative
!> \f$ 
!>   \frac{\partial}{\partial y} \frac{cosh( x y)^2}{cosh(x)^2}
!> \f$ 
!  
!> @param[in]  X  (argument)
!> @param[in]  Y  (argument)
!****************************************************************************

   FUNCTION DCOSH2COSHY( X, Y ) 
      USE prec
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
#ifdef qd_emulate
      USE qdmodule
#endif
      QDPREAL :: X, Y 
      QDPREAL :: DCOSH2COSHY
   
      ! need to carefully evaluate this function 
      DCOSH2COSHY = QDPP( 2 ) * X * COSH2COSH( X, Y ) * TH( X*Y )  
   
   END FUNCTION DCOSH2COSHY

!****************************************************************************
! DESCRIPTION:
!>calculates following derivative
!> \f$ 
!>   \frac{\partial}{\partial y} \frac{cosh( x y)}{cosh(x)}
!> \f$ 
!  
!> @param[in]  X  (argument)
!> @param[in]  Y  (argument)
!****************************************************************************

   FUNCTION DCOSHCOSHY( X, Y ) 
      USE prec
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
#ifdef qd_emulate
      USE qdmodule
#endif
      QDPREAL :: X, Y 
      QDPREAL :: DCOSHCOSHY
   
      !
      DCOSHCOSHY = X * SINHCOSH( X, Y ) 
   
   END FUNCTION DCOSHCOSHY


!=========================================================================
! Basis definitions, squared functions basically
!=========================================================================
   FUNCTION EXPF2( X, L )
      QDPREAL               :: EXPF2
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'expf2' )
      
      EXPF2 = EXPF(2*X,L)
      
      PROFILING_STOP( 'expf2' )

   END FUNCTION EXPF2

   FUNCTION D_EXPF2_DZ( X, L )
      QDPREAL               :: D_EXPF2_DZ
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd_expf2_dz' )
      
      D_EXPF2_DZ = QDPP(2)*D_EXPF_DZ(QDPP(2)*X,L)
      
      PROFILING_STOP( 'd_expf2_dz' )

   END FUNCTION D_EXPF2_DZ

   FUNCTION D2_EXPF2_DZ2( X, L )
      QDPREAL               :: D2_EXPF2_DZ2
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd2_expf2_dz2' )
      
      D2_EXPF2_DZ2 = QDPP(4)*D2_EXPF_DZ2(QDPP(2)*X,L)
      
      PROFILING_STOP( 'd2_expf2_dz2' )

   END FUNCTION D2_EXPF2_DZ2

   FUNCTION D_EXPF2_DL( X, L )
      QDPREAL               :: D_EXPF2_DL
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd_expf2_dl' )
      
      D_EXPF2_DL = D_EXPF_DL(QDPP(2)*X,L)
      
      PROFILING_STOP( 'd_expf2_dl' )

   END FUNCTION D_EXPF2_DL

   FUNCTION UFREQ2( X, L )
      QDPREAL               :: UFREQ2
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'ufreq2' )
      
      UFREQ2 = UFREQ(X,L)**2
      
      PROFILING_STOP( 'ufreq2' )

   END FUNCTION UFREQ2

   FUNCTION D_UFREQ2_DZ( X, L )
      QDPREAL               :: D_UFREQ2_DZ
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd_ufreq2_dz' )
      
      D_UFREQ2_DZ = QDPP(2)*UFREQ(X,L)*D_UFREQ_DZ(X,L)
      
      PROFILING_STOP( 'd_ufreq2_dz' )

   END FUNCTION D_UFREQ2_DZ

   FUNCTION D_UFREQ2_DL( X, L )
      QDPREAL               :: D_UFREQ2_DL
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd_ufreq2_dl' )
      
      D_UFREQ2_DL = QDPP(2)*UFREQ(X,L)*D_UFREQ_DL(X,L)
      
      PROFILING_STOP( 'd_ufreq2_dl' )

   END FUNCTION D_UFREQ2_DL

   FUNCTION D2_UFREQ2_DZ2( X, L )
      QDPREAL               :: D2_UFREQ2_DZ2
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd2_ufreq2_dz2' )
      
      D2_UFREQ2_DZ2 = QDPP(2)*( D_UFREQ_DZ(X,L)**2 &
                                   + UFREQ(X,L)*D2_UFREQ_DZ2(X,L) )
      
      PROFILING_STOP( 'd2_ufreq2_dz2' )

   END FUNCTION D2_UFREQ2_DZ2

   FUNCTION VFREQ2( X, L )
      QDPREAL               :: VFREQ2
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'vfreq2' )
      
      VFREQ2 = VFREQ(X,L)**2
      
      PROFILING_STOP( 'vfreq2' )

   END FUNCTION VFREQ2

   FUNCTION D_VFREQ2_DZ( X, L )
      QDPREAL               :: D_VFREQ2_DZ
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd_vfreq2_dz' )
      
      D_VFREQ2_DZ = QDPP(2)*VFREQ(X,L)*D_VFREQ_DZ(X,L)
      
      PROFILING_STOP( 'd_vfreq2_dz' )

   END FUNCTION D_VFREQ2_DZ

   FUNCTION D_VFREQ2_DL( X, L )
      QDPREAL               :: D_VFREQ2_DL
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd_vfreq2_dl' )
      
      D_VFREQ2_DL = QDPP(2)*VFREQ(X,L)*D_VFREQ_DL(X,L)
      
      PROFILING_STOP( 'd_vfreq2_dl' )

   END FUNCTION D_VFREQ2_DL

   FUNCTION D2_VFREQ2_DZ2( X, L )
      QDPREAL               :: D2_VFREQ2_DZ2
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd2_vfreq2_dz2' )
      
      D2_VFREQ2_DZ2 = QDPP(2)*( D_VFREQ_DZ(X,L)**2 &
                                   + VFREQ(X,L)*D2_VFREQ_DZ2(X,L) )
      
      PROFILING_STOP( 'd2_vfreq2_dz2' )

   END FUNCTION D2_VFREQ2_DZ2

! T>0 basis functions 
   FUNCTION UTIME2( X, L )
      QDPREAL               :: UTIME2
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'utime2' )
      
      UTIME2 = UTIME(X,L)**2
      
      PROFILING_STOP( 'utime2' )

   END FUNCTION UTIME2

   FUNCTION D_UTIME2_DZ( X, L )
      QDPREAL               :: D_UTIME2_DZ
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd_utime2_dz' )
      
      D_UTIME2_DZ = QDPP(2)*UTIME(X,L)*D_UTIME_DZ(X,L)
      
      PROFILING_STOP( 'd_utime2_dz' )

   END FUNCTION D_UTIME2_DZ

   FUNCTION D2_UTIME2_DZ2( X, L )
      QDPREAL               :: D2_UTIME2_DZ2
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd2_utime2_dz2' )
      
      D2_UTIME2_DZ2 = QDPP(2)*( D_UTIME_DZ(X,L)**2 &
                                   + UTIME(X,L)*D2_UTIME_DZ2(X,L) )
      
      PROFILING_STOP( 'd2_utime2_dz2' )

   END FUNCTION D2_UTIME2_DZ2

   FUNCTION D_UTIME2_DL( X, L )
      QDPREAL               :: D_UTIME2_DL
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd_utime2_dl' )
      
      D_UTIME2_DL = QDPP(2)*UTIME(X,L)*D_UTIME_DL(X,L)
      
      PROFILING_STOP( 'd_utime2_dl' )

   END FUNCTION D_UTIME2_DL

! T>0 basis functions 
   FUNCTION UTANH2( X, L )
      QDPREAL               :: UTANH2
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'utanh2' )
      
      UTANH2 = UTANH(X,L)**2
      
      PROFILING_STOP( 'utanh2' )

   END FUNCTION UTANH2

   FUNCTION D_UTANH2_DZ( X, L )
      QDPREAL               :: D_UTANH2_DZ
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd_utanh2_dz' )
      
      D_UTANH2_DZ = QDPP(2)*UTANH(X,L)*D_UTANH_DZ(X,L)
      
      PROFILING_STOP( 'd_utanh2_dz' )

   END FUNCTION D_UTANH2_DZ

   FUNCTION D2_UTANH2_DZ2( X, L )
      QDPREAL               :: D2_UTANH2_DZ2
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd2_utanh2_dz2' )
      
      D2_UTANH2_DZ2 = QDPP(2)*( D_UTANH_DZ(X,L)**2 &
                                   + UTANH(X,L)*D2_UTANH_DZ2(X,L) )
      
      PROFILING_STOP( 'd2_utanh2_dz2' )

   END FUNCTION D2_UTANH2_DZ2

   FUNCTION D_UTANH2_DL( X, L )
      QDPREAL               :: D_UTANH2_DL
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd_utanh2_dl' )
      
      D_UTANH2_DL = QDPP(2)*UTANH(X,L)*D_UTANH_DL(X,L)
      
      PROFILING_STOP( 'd_utanh2_dl' )

   END FUNCTION D_UTANH2_DL

! T>0 basis functions 
   FUNCTION VTANH2( X, L )
      QDPREAL               :: VTANH2
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'vtanh2' )
      
      VTANH2 = VTANH(X,L)**2
      
      PROFILING_STOP( 'vtanh2' )

   END FUNCTION VTANH2

   FUNCTION D_VTANH2_DZ( X, L )
      QDPREAL               :: D_VTANH2_DZ
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd_vtanh2_dz' )
      
      D_VTANH2_DZ = QDPP(2)*VTANH(X,L)*D_VTANH_DZ(X,L)
      
      PROFILING_STOP( 'd_vtanh2_dz' )

   END FUNCTION D_VTANH2_DZ

   FUNCTION D2_VTANH2_DZ2( X, L )
      QDPREAL               :: D2_VTANH2_DZ2
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd2_vtanh2_dz2' )
      
      D2_VTANH2_DZ2 = QDPP(2)*( D_VTANH_DZ(X,L)**2 &
                                   + VTANH(X,L)*D2_VTANH_DZ2(X,L) )
      
      PROFILING_STOP( 'd2_vtanh2_dz2' )

   END FUNCTION D2_VTANH2_DZ2

   FUNCTION D_VTANH2_DL( X, L )
      QDPREAL               :: D_VTANH2_DL
      QDPREAL , INTENT(IN)  :: X ! sampling value
      QDPREAL , INTENT(IN)  :: L ! coefficient defining function

      PROFILING_START( 'd_vtanh2_dl' )
      
      D_VTANH2_DL = QDPP(2)*VTANH(X,L)*D_VTANH_DL(X,L)
      
      PROFILING_STOP( 'd_vtanh2_dl' )

   END FUNCTION D_VTANH2_DL


END MODULE minimax_functions2D
