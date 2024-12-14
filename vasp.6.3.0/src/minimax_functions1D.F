#include "symbol.inc"
!***********************************************************************
! 
! MODULE: minimax_functions
!
!> @author
!> Merzuk Kaltak, VASP GmbH
! 
! DESCRIPTION:
!> contains all one-dimensional functions 
!>
!***********************************************************************
MODULE minimax_functions1D
   USE prec
#ifdef qd_emulate
   USE qdmodule
#define ACC3     qdreal('1.E-8')
#define PIQ      qdpi()
   IMPLICIT NONE 
#else
   IMPLICIT NONE 
   !truncation threshold for Taylor expansion for finite temperature error functions
   QDPREAL, PARAMETER :: ACC3=1.E-8_qd   
   QDPREAL ,PARAMETER :: PIQ=3.1415926535897932384626433832795028841971693993751_qd
#endif
   ! maximum exponent to base 2 is obtained from intrinsic function 
   ! MAXEXPONENT, to base e the maximum exponent is MAXEXPONENT(X) * LN 2
   ! for QD MAXEXPONENT return usually 16384* LN 2 = 11356 
   ! we use a 1/2.5 part of this (functions are being squared) ->4500
   INTEGER, PARAMETER :: MAXE=100

   CONTAINS

!***********************************************************************
! DESCRIPTION:
!>calculates following function:
!> \f$ 
!>     \frac{1}{z}
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
  FUNCTION INVERSE_Z(Z)
      QDPREAL, INTENT(IN) :: Z
      QDPREAL             :: INVERSE_Z
       
      IF (ABS(Z)<QDPP(1.E-9)) THEN
         INVERSE_Z=0
      ELSE
         INVERSE_Z=QDPD(1)/Z
      ENDIF
   END FUNCTION INVERSE_Z

!***********************************************************************
! DESCRIPTION:
!>calculates following function:
!> \f$ 
!>     -\frac{1}{z^2}
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
  FUNCTION D1_INVERSE_Z(Z)
      QDPREAL, INTENT(IN) :: Z
      QDPREAL             :: D1_INVERSE_Z
       
      IF (ABS(Z)<QDPP(1.E-9)) THEN
         D1_INVERSE_Z=0
      ELSE
         D1_INVERSE_Z=-QDPD(1)/(Z*Z)
      ENDIF
   END FUNCTION D1_INVERSE_Z

!***********************************************************************
! DESCRIPTION:
!>calculates following function:
!> \f$ 
!>     \frac{2}{z^3}
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
  FUNCTION D2_INVERSE_Z(Z)
      QDPREAL, INTENT(IN) :: Z
      QDPREAL             :: D2_INVERSE_Z
       
      IF (ABS(Z)<QDPP(1.E-9)) THEN
         D2_INVERSE_Z=0
      ELSE
         D2_INVERSE_Z=QDPD(2)/(Z**3)
      ENDIF
   END FUNCTION D2_INVERSE_Z

!****************************************************************************
!
! DESCRIPTION:
!>calculates the hyperbolic tangent tanh
!  
!> @param[in]  Z  (argument)
! 
!****************************************************************************

   FUNCTION TH( Z ) 
      USE prec
#ifdef qd_emulate
      USE qdmodule
#endif
      QDPREAL :: Z
      QDPREAL :: TH
   
      PROFILING_START( 'th' )
      ! error is smaller than 7.4*10^{-44}
      IF ( ABS( Z ) >  50) THEN
         TH = SIGN_QD(Z)
      ELSE IF ( ABS( Z) < QDPP(1E-8) )THEN
         TH = Z**2/QDPP(2)-Z**3/QDPP(24)+Z**5/QDPP(240)
      ELSE
         TH = TANH( Z ) 
      ENDIF
      PROFILING_STOP( 'th' )
   END FUNCTION TH

!****************************************************************************
!
! DESCRIPTION:
!>calculates the derivative of hyperbolic tangent tanh
!  
!> @param[in]  Z  (argument)
! 
!****************************************************************************

   FUNCTION DTH( Z ) 
      USE prec
#ifdef qd_emulate
      USE qdmodule
#endif
      QDPREAL :: Z
      QDPREAL :: DTH
   
      PROFILING_START( 'dth' )
      ! error is smaller than 7.4*10^{-44}
      IF ( ABS( Z ) >  50) THEN
         DTH = 0
      ELSE
         DTH = QDPP(1)-TH(Z)**2
      ENDIF
      PROFILING_STOP( 'dth' )
   END FUNCTION DTH

!****************************************************************************
!
! DESCRIPTION:
!>calculates: 
!>  tanh(z/2)/2
!  
!> @param[in]  Z  (argument)
! 
!****************************************************************************

   FUNCTION F_TANH( Z ) 
      USE prec
      QDPREAL, INTENT(IN) :: Z
      QDPREAL :: F_TANH
      PROFILING_START( 'f_tanh' )
      F_TANH = TH( Z/QDPP(2) ) /QDPP(2)
      PROFILING_STOP( 'f_tanh' )
   END FUNCTION F_TANH

!****************************************************************************
!
! DESCRIPTION:
!>calculates first derivative of: 
!>  tanh(z/2)/2
!  
!> @param[in]  Z  (argument)
! 
!****************************************************************************

   FUNCTION D1_F_TANH( Z ) 
      USE prec
      QDPREAL, INTENT(IN) :: Z
      QDPREAL :: D1_F_TANH
      PROFILING_START( 'd1_f_tanh' )
      D1_F_TANH = (QDPP(1)-QDPP(4)*F_TANH(Z)**2)/QDPP(4)
      PROFILING_STOP( 'd1_f_tanh' )
   END FUNCTION D1_F_TANH

!****************************************************************************
!
! DESCRIPTION:
!>calculates secon derivative of: 
!>  tanh(z/2)/2
!  
!> @param[in]  Z  (argument)
! 
!****************************************************************************

   FUNCTION D2_F_TANH( Z ) 
      USE prec
      QDPREAL, INTENT(IN) :: Z
      QDPREAL :: T
      QDPREAL :: D2_F_TANH
      PROFILING_START( 'd2_f_tanh' )
      T = F_TANH(Z)
      D2_F_TANH = -QDPP(2)*T*(QDPP(1)-QDPP(4)*T**2)/QDPP(4)
      PROFILING_STOP( 'd2_f_tanh' )
   END FUNCTION D2_F_TANH

!***********************************************************************
! DESCRIPTION:
!>calculates following function:
!> \f$ 
!>     \frac{\tanh\frac{z}2}{z}
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION F_THZ_OVER_Z( Z ) 
      QDPREAL               :: F_THZ_OVER_Z
      QDPREAL , INTENT(IN)  :: Z 
      PROFILING_START( 'f_thz_over_z' )

      ! use a taylor expansion for small arguments 
      IF ( ABS( Z ) < QDPP(1E-8) ) THEN
         F_THZ_OVER_Z=QDPP(0.5)-Z**2/QDPP(24)+Z**4/QDPP(240)&
                    -QDPP(17)*Z**6/QDPP(40320)
      ELSE
         F_THZ_OVER_Z =TH( Z/QDPP(2) )/Z
      ENDIF
      PROFILING_STOP( 'f_thz_over_z' )
   END FUNCTION F_THZ_OVER_Z

!***********************************************************************
! DESCRIPTION:
!>calculates first derivative of following function:
!> \f$ 
!>     \frac{\tanh\frac{z}2}{z}
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D1_THZ_OVER_Z( Z ) 
      QDPREAL               :: D1_THZ_OVER_Z
      QDPREAL , INTENT(IN)  :: Z 
      QDPREAL               :: THZ
      PROFILING_START( 'd1_thz_over_z' )

      ! use a taylor expansion for small arguments 
      IF ( ABS( Z ) < QDPP(1E-8) ) THEN
         D1_THZ_OVER_Z= -Z/QDPP(12)+Z**3/QDPP(60)-QDPP(17)*Z**5/QDPP(6720)
      ELSE
         THZ = TH(Z/QDPP(2))
         D1_THZ_OVER_Z= (QDPP(1)-THZ**2)/(QDPP(2)*Z) - THZ/Z**2
      ENDIF
      PROFILING_STOP( 'd1_thz_over_z' )
   END FUNCTION D1_THZ_OVER_Z

!***********************************************************************
! DESCRIPTION:
!>calculates second derivative of following function:
!> \f$ 
!>     \frac{\tanh\frac{z}2}{z}
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D2_THZ_OVER_Z( Z ) 
      QDPREAL               :: D2_THZ_OVER_Z
      QDPREAL , INTENT(IN)  :: Z 
      QDPREAL               :: THZ
      PROFILING_START( 'd2_thz_over_z' )

      ! use a taylor expansion for small arguments 
      IF ( ABS( Z ) < QDPP(1E-8) ) THEN
         D2_THZ_OVER_Z= -QDPP(1)/QDPP(12)                    &
                        +Z**2/QDPP(20)                       &
                        -QDPP(17)*Z**4/QDPP(1344)            &
                        +QDPP(31)*Z**6/QDPP(12960)
      ELSE
         THZ = TH(Z/QDPP(2))
         D2_THZ_OVER_Z= QDPP(2)*THZ/Z**3                   &
                      - (QDPP(1)-THZ**2)/Z**2              &
                      - THZ*(QDPP(1)-THZ**2)/(QDPP(2)*Z)   
      ENDIF
      PROFILING_STOP( 'd2_thz_over_z' )
   END FUNCTION D2_THZ_OVER_Z

!***********************************************************************
! DESCRIPTION:
!>calculates following function:
!> \f$ 
!>     \frac{1}{8}(1-\tanh^2\frac{z}{2}+\frac{2}{z}\tanh\frac{z}{2}})
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION F_BETA_LRS( Z ) 
      QDPREAL               :: F_BETA_LRS
      QDPREAL , INTENT(IN)  :: Z 
      QDPREAL               :: ARG 
      PROFILING_START( 'f_beta_lrs' )

      ARG = Z/QDPP(2)
      ! use a taylor expansion for small arguments 
      IF ( ABS( Z ) < QDPP(1E-8) ) THEN
         F_BETA_LRS=QDPP(0.25)-Z**2/QDPP(24)+Z**4/QDPP(160)&
                       -QDPP(17)*Z**6/QDPP(20160)
      ELSE
         F_BETA_LRS =( QDPP(0.125)-QDPP(0.125)*TH( ARG )**2 )&
                        +( QDPP(2)*TH( ARG )/(QDPP(8)*Z) )
      ENDIF
      PROFILING_STOP( 'f_beta_lrs' )
   END FUNCTION F_BETA_LRS

!***********************************************************************
! DESCRIPTION:
!>calculates following function:
!> \f$ 
!>     \frac{1}{8}(1-\tanh^2\frac{z}{2}+\frac{2}{z}\tanh\frac{z}{2}})
!>     -\frac{1}{z^2}\tanh^2\frac{z}{2}
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION F_BETA_LRS_REG( Z ) 
      QDPREAL               :: F_BETA_LRS_REG
      QDPREAL , INTENT(IN)  :: Z 
      QDPREAL               :: ARG 
      PROFILING_START( 'f_beta_lrs_reg' )

      ARG = Z/QDPP(2)
      ! use a taylor expansion for small arguments 
      IF ( ABS( Z ) < QDPP(1E-8) ) THEN
         F_BETA_LRS_REG=Z**4/QDPP(2880)-Z**6/QDPP(13440)
      ELSE
         ARG = TH( ARG ) 
         F_BETA_LRS_REG =( QDPP(0.125)-QDPP(0.125)*ARG**2 )&
                        +( QDPP(2)*ARG/(QDPP(8)*Z) ) &
                        -( ARG/Z )**2
      ENDIF
      PROFILING_STOP( 'f_beta_lrs_reg' )
   END FUNCTION F_BETA_LRS_REG

!***********************************************************************
! DESCRIPTION:
!>calculates first two derivatives of following function:
!> \f$ 
!>     \frac{1}{8}(1-\tanh^2\frac{z}{2}+\frac{2}{z}\tanh\frac{z}{2}})
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION F_BETA_LRS_DER( Z ) 
      QDPREAL               :: F_BETA_LRS_DER(2)
      QDPREAL , INTENT(IN)  :: Z 
      QDPREAL               :: ZP(2),T
      PROFILING_START('f_beta_lrs_der')

      T = TH( Z/QDPD(2) )
      ! For small arguments use taylor expansion 
      IF ( ABS( Z ) < ACC3 ) THEN
         ! for |Z|<1.E-6 the error of expansion and exact value is < 10^-39

         ! first derivative  w.r.t. Z
         ZP(1)= -Z/QDPP(12) + Z**3/QDPP(40) - (QDPP(17)*Z**5)/QDPP(3360)
               
         ! second derivative  w.r.t. Z
         ZP(2)= -QDPP(1)/QDPP(12) + QDPP(3)*Z**2/QDPP(40) - (QDPP(17)*Z**4)/QDPP(672)
      ELSE
         ! write functions in terms of Tanh[ Z/2 ] 
         ! first derivative  w.r.t. Z
         ZP(1)= ( QDPP( 1 )/( QDPP( 8 )*Z  ) ) &
              - ( T/QDPP( 8 ) )          &
              - T*(QDPP(1)/( QDPP( 2 )*Z) )**2  &
              - ( T**2/( QDPP( 8 )*Z ) ) &
              + ( T**3/( QDPP( 8 ) ) )

         ! second derivative  w.r.t. Z
         ZP(2)= -QDPP(0.0625)                     &
                -( QDPP(1)/(QDPP(2)*Z))**2        &
                +( T/QDPP(2) )*( QDPD(1)/Z )**3   &
                -( T/(QDPP(8)*Z) )                &
                +( T/QDPP(2) )**2                 &
                +( T/(QDPP(2)*Z) )**2             &
                +( T**3/(QDPP(8)*Z) )             &
                -( (QDPP(3)*T**4)/QDPP(16) )     
      ENDIF 

      F_BETA_LRS_DER(1) = ZP(1)
      F_BETA_LRS_DER(2) = ZP(2)

      PROFILING_STOP('f_beta_lrs_der')
   END FUNCTION F_BETA_LRS_DER

!***********************************************************************
! DESCRIPTION:
!>calculates first derivative of following function:
!> \f$ 
!>     \frac{1}{8}(1-\tanh^2\frac{z}{2}+\frac{2}{z}\tanh\frac{z}{2}})
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D1_F_BETA_LRS( Z ) 
      QDPREAL               :: D1_F_BETA_LRS
      QDPREAL , INTENT(IN)  :: Z 
      QDPREAL               :: T,DT
      PROFILING_START( 'd1_f_beta_lrs' )

      T = TH( Z/QDPD(2) )
      ! For small arguments use taylor expansion 
      IF ( ABS( Z ) < ACC3 ) THEN
         ! for |Z|<1.E-6 the error of expansion and exact value is < 10^-39

         ! first derivative  w.r.t. Z
         D1_F_BETA_LRS= -Z/QDPP(12) + Z**3/QDPP(40) - (QDPP(17)*Z**5)/QDPP(3360)
      ELSE
         ! write functions in terms of Tanh[ Z/2 ] 
         ! first derivative  w.r.t. Z
         D1_F_BETA_LRS= ( QDPP( 1 )/( QDPP( 8 )*Z  ) ) &
              - ( T/QDPP( 8 ) )          &
              - T*(QDPP(1)/( QDPP( 2 )*Z) )**2  &
              - ( T**2/( QDPP( 8 )*Z ) ) &
              + ( T**3/( QDPP( 8 ) ) )
      ENDIF 

      PROFILING_STOP( 'd1_f_beta_lrs' )
   END FUNCTION D1_F_BETA_LRS

!***********************************************************************
! DESCRIPTION:
!>calculates second derivative of following function:
!> \f$ 
!>     \frac{1}{8}(1-\tanh^2\frac{z}{2}+\frac{2}{z}\tanh\frac{z}{2}})
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D2_F_BETA_LRS( Z ) 
      QDPREAL               :: D2_F_BETA_LRS
      QDPREAL , INTENT(IN)  :: Z 
      QDPREAL               :: T
      PROFILING_START( 'd2_f_beta_lrs' )

      T = TH( Z/QDPD(2) )
      ! For small arguments use taylor expansion 
      IF ( ABS( Z ) < ACC3 ) THEN
         ! for |Z|<1.E-6 the error of expansion and exact value is < 10^-39
         ! second derivative  w.r.t. Z
         D2_F_BETA_LRS= -QDPP(1)/QDPP(12) + QDPP(3)*Z**2/QDPP(40) - (QDPP(17)*Z**4)/QDPP(672)
      ELSE
         ! second derivative  w.r.t. Z
         D2_F_BETA_LRS= -QDPP(0.0625)             &
                -( QDPP(1)/(QDPP(2)*Z))**2        &
                +( T/QDPP(2) )*( QDPD(1)/Z )**3   &
                -( T/(QDPP(8)*Z) )                &
                +( T/QDPP(2) )**2                 &
                +( T/(QDPP(2)*Z) )**2             &
                +( T**3/(QDPP(8)*Z) )             &
                -( (QDPP(3)*T**4)/QDPP(16) )     
      ENDIF 

      PROFILING_STOP( 'd2_f_beta_lrs' )
   END FUNCTION D2_F_BETA_LRS

!***********************************************************************
! DESCRIPTION:
!>calculates frst two dervatives of following function:
!> \f$ 
!>     \frac{1}{8}(1-\tanh^2\frac{z}{2}+\frac{2}{z}\tanh\frac{z}{2}})
!>     -\frac{1}{z^2}\tanh^2\frac{z}{2}
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION F_BETA_LRS_REG_DER( Z ) 
      QDPREAL               :: F_BETA_LRS_REG_DER(2)
      QDPREAL , INTENT(IN)  :: Z 
      QDPREAL               :: ZP(2),T,TZ
      PROFILING_START( 'beta2_reg_der' )

      ! For small arguments use taylor expansion 
      IF ( ABS( Z ) < ACC3 ) THEN
         ! first derivative  w.r.t. Z
         ZP(1)= Z**3/QDPP(720)-Z**5/QDPP(2240)
         ! second derivative  w.r.t. Z
         ZP(2)=Z**2/QDPP(240)-Z**4/QDPP(448)+QDPP(7)*Z**6/QDPP(10800)
      ELSE
         T = TH( Z/QDPD(2) )
         TZ = T/Z
         ! write functions in terms of Tanh[ Z/2 ] 
         ! first derivative  w.r.t. Z
         ZP(1)=                                   &
              - ( T/QDPP( 8 ) )                   &
              + ( T**3/QDPP( 8 ) )                &
              + ( QDPP(2)*TZ**2/Z  )              &
              - ( QDPP(5)/QDPP(4)*TZ/Z )          &
              + ( T*TZ**2 )                       &
              + ( QDPP( 1 )/( QDPP( 8 )*Z  ) )    &
              - ( T*TZ/QDPP( 8 ) )    

         ! second derivative  w.r.t. Z
         ZP(2)=                                   &
              -  QDPP(0.0625)                     &
              + ( T**2/QDPP( 4 ) )                &
              - (QDPP(3)/QDPP(16))*T**4           &
              - (QDPP(6)*(TZ**2/Z**2 ))           &
              + (QDPP(9)*TZ/(QDPP(2)*Z**2 ))      &
              - (QDPP(4)*TZ**3)                   &
              - (QDPP(3)/(QDPP(4)*Z**2 ))         &
              + (QDPP(9)/QDPP(4)*TZ**2)           &
              - (QDPP(3)/QDPP(2)*T**2*TZ**2)      &
              - (TZ/QDPP(8) )                     &
              + (TZ*T**2/QDPP(8) )      
      ENDIF 

      F_BETA_LRS_REG_DER(1) = ZP(1)
      F_BETA_LRS_REG_DER(2) = ZP(2)
      PROFILING_STOP( 'beta2_reg_der' )
   END FUNCTION F_BETA_LRS_REG_DER

!***********************************************************************
! DESCRIPTION:
!>calculates first dervatives of following function:
!> \f$ 
!>     \frac{1}{8}(1-\tanh^2\frac{z}{2}+\frac{2}{z}\tanh\frac{z}{2}})
!>     -\frac{1}{z^2}\tanh^2\frac{z}{2}
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D1_F_BETA_LRS_REG( Z ) 
      QDPREAL               :: D1_F_BETA_LRS_REG
      QDPREAL , INTENT(IN)  :: Z 
      QDPREAL               :: ZP(2),T,TZ
      PROFILING_START( 'd1_f_beta_lrs' )

      ! For small arguments use taylor expansion 
      IF ( ABS( Z ) < ACC3 ) THEN
         ! first derivative  w.r.t. Z
         D1_F_BETA_LRS_REG= Z**3/QDPP(720)-Z**5/QDPP(2240)
      ELSE
         T = TH( Z/QDPD(2) )
         TZ = T/Z
         ! write functions in terms of Tanh[ Z/2 ] 
         ! first derivative  w.r.t. Z
         D1_F_BETA_LRS_REG=                       &
              - ( T/QDPP( 8 ) )                   &
              + ( T**3/QDPP( 8 ) )                &
              + ( QDPP(2)*TZ**2/Z  )              &
              - ( QDPP(5)/QDPP(4)*TZ/Z )          &
              + ( T*TZ**2 )                       &
              + ( QDPP( 1 )/( QDPP( 8 )*Z  ) )    &
              - ( T*TZ/QDPP( 8 ) )    
      ENDIF 
      PROFILING_STOP( 'd1_f_beta_lrs' )
   END FUNCTION D1_F_BETA_LRS_REG

!***********************************************************************
! DESCRIPTION:
!>calculates second dervatives of following function:
!> \f$ 
!>     \frac{1}{8}(1-\tanh^2\frac{z}{2}+\frac{2}{z}\tanh\frac{z}{2}})
!>     -\frac{1}{z^2}\tanh^2\frac{z}{2}
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D2_F_BETA_LRS_REG( Z ) 
      QDPREAL               :: D2_F_BETA_LRS_REG
      QDPREAL , INTENT(IN)  :: Z 
      QDPREAL               :: T,TZ
      PROFILING_START( 'd2_f_beta_lrs' )

      ! For small arguments use taylor expansion 
      IF ( ABS( Z ) < ACC3 ) THEN
         ! first derivative  w.r.t. Z
         D2_F_BETA_LRS_REG = Z**2/QDPP(240)-Z**4/QDPP(448)+QDPP(7)*Z**6/QDPP(10800)
      ELSE
         T = TH( Z/QDPD(2) )
         TZ = T/Z
         ! write functions in terms of Tanh[ Z/2 ] 
         ! first derivative  w.r.t. Z
         D2_F_BETA_LRS_REG=                       &
              -  QDPP(0.0625)                     &
              + ( T**2/QDPP( 4 ) )                &
              - (QDPP(3)/QDPP(16))*T**4           &
              - (QDPP(6)*(TZ**2/Z**2 ))           &
              + (QDPP(9)*TZ/(QDPP(2)*Z**2 ))      &
              - (QDPP(4)*TZ**3)                   &
              - (QDPP(3)/(QDPP(4)*Z**2 ))         &
              + (QDPP(9)/QDPP(4)*TZ**2)           &
              - (QDPP(3)/QDPP(2)*T**2*TZ**2)      &
              - (TZ/QDPP(8) )                     &
              + (TZ*T**2/QDPP(8) )      
      ENDIF 
      PROFILING_STOP( 'd2_f_beta_lrs' )
   END FUNCTION D2_F_BETA_LRS_REG

!***********************************************************************
! DESCRIPTION:
!>calculates following function:
!> \f$ 
!>     \frac{1}{4}( \frac{\tanh\frac{z}{2}}}{z}
!>                 -\frac12(1-\tanh^2\frac{z}{2}) )
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION F_BETA_RRS( Z ) 
      QDPREAL               :: F_BETA_RRS
      QDPREAL , INTENT(IN)  :: Z 
      PROFILING_START( 'f_beta_rrs' )

      ! use a taylor expansion for small arguments 
      IF ( ABS( Z ) < QDPP(1E-8) ) THEN
         F_BETA_RRS = Z**2/QDPP(48) - Z**4/QDPP(240) + QDPP(17)*Z**6/QDPP(26880)
      ELSE
         F_BETA_RRS = ( QDPP(2)*F_THZ_OVER_Z( Z ) &
                    - ( QDPP(1)-TH( Z/QDPP(2) )**2 ) )/QDPP(8)
      ENDIF
      PROFILING_STOP( 'f_beta_rrs' )
   END FUNCTION F_BETA_RRS

!***********************************************************************
! DESCRIPTION:
!>calculates first derivative of following function:
!> \f$ 
!>     \frac{1}{4}( \frac{\tanh\frac{z}{2}}}{z}
!>                 -\frac12(1-\tanh^2\frac{z}{2}) )
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D1_F_BETA_RRS( Z ) 
      QDPREAL              :: D1_F_BETA_RRS
      QDPREAL, INTENT(IN)  :: Z 
      ! local  
      QDPREAL              :: TZ, T

      PROFILING_START( 'd1_f_beta_rrs' )

      ! use a taylor expansion for small arguments 
      IF ( ABS( Z ) < QDPP(1E-8) ) THEN
         D1_F_BETA_RRS=Z/QDPP(24)-Z**3/QDPP(60)+QDPP(17)*Z**5/QDPP(4480)
      ELSE
         TZ = F_THZ_OVER_Z( Z )
         T = TH( Z/QDPP(2) ) 
         !T/8 - T^3/8 - T/(4 z^2) + 1/(8 z) - T^2/(8 z)
         D1_F_BETA_RRS = (T-T**3)/QDPP(8)-(QDPP(2)*TZ/Z+T*TZ)/QDPP(8) &
                       + QDPP(1)/(QDPP(8)*Z)
      ENDIF
      PROFILING_STOP( 'd1_f_beta_rrs' )
   END FUNCTION D1_F_BETA_RRS

!***********************************************************************
! DESCRIPTION:
!>calculates second derivative of following function:
!> \f$ 
!>     \frac{1}{4}( \frac{\tanh\frac{z}{2}}}{z}
!>                 -\frac12(1-\tanh^2\frac{z}{2}) )
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D2_F_BETA_RRS( Z ) 
      QDPREAL              :: D2_F_BETA_RRS
      QDPREAL, INTENT(IN)  :: Z 
      ! local  
      QDPREAL              :: TZ, T

      PROFILING_START( 'd2_f_beta_rrs' )

      ! use a taylor expansion for small arguments 
      IF ( ABS( Z ) < QDPP(1E-8) ) THEN
         D2_F_BETA_RRS = QDPP(1)/QDPP(24)-Z**2/QDPP(20) &
                       + QDPP(17)*Z**4/QDPP(896) &
                       - QDPP(31)*Z**6/QDPP(6480)
      ELSE
         TZ = F_THZ_OVER_Z( Z )
         T = TH( Z/QDPP(2) ) 
         !1/16 - T^2/4 + (3 T^4)/16 + T/(2 z^3) - 1/(4 z^2) + T^2/(4 z^2) - T/( 8 z) + T^3/(8 z)
         D2_F_BETA_RRS = QDPP(1)/QDPP(16)       &
                       - T**2/QDPP( 4 )         &
                       + QDPP(3)*T**4/QDPP(16)  &
                       + TZ/(QDPP(2)*Z**2)      &
                       - QDPP(1)/(QDPP(4)*Z**2) &
                       + TZ**2/QDPP(4)          &
                       - TZ/QDPP(8)             &
                       + T**2*TZ/QDPP(8)      
      ENDIF
      PROFILING_STOP( 'd2_f_beta_rrs' )
   END FUNCTION D2_F_BETA_RRS

   FUNCTION F_BETA_RRS_DER( Z ) 
      QDPREAL               :: F_BETA_RRS_DER(2)
      QDPREAL , INTENT(IN)  :: Z 
      QDPREAL               :: T, TZ 
      PROFILING_START( 'f_beta_rrs_der' )

      IF ( ABS( Z ) < QDPP(1E-8) ) THEN
         F_BETA_RRS_DER(1)=Z/QDPP(24)-Z**3/QDPP(60)+QDPP(17)*Z**5/QDPP(4480)
         F_BETA_RRS_DER(2) = QDPP(1)/QDPP(24)-Z**2/QDPP(20) &
                       + QDPP(17)*Z**4/QDPP(896) &
                       - QDPP(31)*Z**6/QDPP(6480)
      ELSE
         TZ = F_THZ_OVER_Z( Z )
         T = TH( Z/QDPP(2) ) 
         !T/8 - T^3/8 - T/(4 z^2) + 1/(8 z) - T^2/(8 z)
         F_BETA_RRS_DER(1) = (T-T**3)/QDPP(8)-(QDPP(2)*TZ/Z+T*TZ)/QDPP(8) &
                       + QDPP(1)/(QDPP(8)*Z)
         F_BETA_RRS_DER(2) = QDPP(1)/QDPP(16)       &
                       - T**2/QDPP( 4 )         &
                       + QDPP(3)*T**4/QDPP(16)  &
                       + TZ/(QDPP(2)*Z**2)      &
                       - QDPP(1)/(QDPP(4)*Z**2) &
                       + TZ**2/QDPP(4)          &
                       - TZ/QDPP(8)             &
                       + T**2*TZ/QDPP(8)      
      ENDIF

      PROFILING_STOP( 'f_beta_rrs_der' )
   END FUNCTION F_BETA_RRS_DER

!***********************************************************************
! DESCRIPTION:
!>calculates this function:
!> \f$ 
!>     \frac{1}{2} ( \coth\frac{z}{2} -\frac{2}{z} )
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION COTH_REG( Z ) 
      QDPREAL               :: COTH_REG
      QDPREAL , INTENT(IN)  :: Z 
      QDPREAL               :: ARG 

      PROFILING_START( 'coth_reg' )

      ARG = Z/QDPP(2)
      ! use a taylor expansion for small arguments 
      IF ( ABS( Z ) < QDPP(0.3) ) THEN
         COTH_REG=ARG/QDPP(6)-ARG**3/QDPP(90)+ARG**5/QDPP(945)
      ELSE
         COTH_REG = (QDPP(1)/TH(ARG)-QDPP(1)/ARG)/QDPP(2)
      ENDIF
 
      PROFILING_STOP( 'coth_reg' )
   END FUNCTION COTH_REG

!***********************************************************************
! DESCRIPTION:
!>calculates derivative of this function:
!> \f$ 
!>     \frac{1}{2} ( \coth\frac{z}{2} -\frac{2}{z} )
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D1_COTH_REG( Z ) 
      QDPREAL               :: D1_COTH_REG
      QDPREAL , INTENT(IN)  :: Z 
      QDPREAL               :: ARG, CT 

      PROFILING_START( 'd1_coth_reg' )

      ARG = Z/QDPP(2)
      ! use a taylor expansion for small arguments 
      IF ( ABS( Z ) < QDPP(0.3) ) THEN
         D1_COTH_REG=QDPP(1)/QDPP(12)   &
                             -ARG**2/QDPP(60)    &
                             +ARG**4/QDPP(378)   &
                             -ARG**6/QDPP(2700) 
      ELSE
         CT=QDPP(1)/TH(ARG)
         D1_COTH_REG=( QDPP(1)/(ARG*ARG) -CT**2 + QDPP(1) )/QDPP(4)
      ENDIF
 
      PROFILING_STOP( 'd1_coth_reg' )
   END FUNCTION D1_COTH_REG

!***********************************************************************
! 
! DESCRIPTION:
!>calculates second derivative of 
!> \f$ 
!>    f''(z)= \frac14(\coth^3\frac{z}{2}-\coth\frac{z}{2}-
!>            \frac{8}{z^3})
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION D2_COTH_REG( Z ) 
      QDPREAL               :: D2_COTH_REG
      QDPREAL , INTENT(IN)  :: Z 
      QDPREAL               :: ARG , CT

      PROFILING_START( 'd2_coth_reg' )

      ARG = Z/QDPP(2)
      ! use a taylor expansion for small arguments 
      IF ( ABS( Z ) < QDPP(0.3) ) THEN
         D2_COTH_REG=-ARG/QDPP(60)      &
                     +ARG**3/QDPP(189)   &
                     -ARG**5/QDPP(900)   &
                     +ARG**7/QDPP(20790) 
      ELSE
         CT=QDPP(1)/TH(ARG)
         D2_COTH_REG=(-QDPP(1)/ARG**3 - CT + CT**3 )/QDPP(4) 
      ENDIF
 
      PROFILING_STOP( 'd2_coth_reg' )
   END FUNCTION D2_COTH_REG

!***********************************************************************
! 
! DESCRIPTION:
!>calculates taylor expansion of cos for low arguments Z
!> \f$ 
!>    f(z)  = 1 - Z^2/2
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION COS_TAYLOR( Z ) 
      QDPREAL               :: COS_TAYLOR
      QDPREAL , INTENT(IN)  :: Z 

      PROFILING_START( 'cos_taylor' )

      COS_TAYLOR = QDPP( 1 ) - Z**2/QDPP( 2 ) 

      PROFILING_STOP( 'cos_taylor' )
   END FUNCTION COS_TAYLOR

!***********************************************************************
! 
! DESCRIPTION:
!>calculates taylor expansion of sin for low arguments Z
!> \f$ 
!>    f(z)  = Z - Z^3/6
!> \f$ 
!> @param[in]  Z  (argument)
!***********************************************************************
   FUNCTION SIN_TAYLOR( Z ) 
      QDPREAL               :: SIN_TAYLOR
      QDPREAL , INTENT(IN)  :: Z 

      PROFILING_START( 'sin_taylor' )

      SIN_TAYLOR = Z - Z**3/QDPP( 6 ) 

      PROFILING_STOP( 'sin_taylor' )
   END FUNCTION SIN_TAYLOR

!****************************************************************************
!
! DESCRIPTION:
!>returns sign of a quad-precision real-valued variable
!  
!> @param[in]  X  (quad precision real)
! 
!****************************************************************************

   PURE FUNCTION SIGN_QD( X ) 
   INTEGER               :: SIGN_QD 
   QDPREAL , INTENT(IN)  :: X 
   ! quad library has no build in functionality for sign function
!#ifdef qd_emulate
      IF ( X > 0 ) THEN
         SIGN_QD = 1
      ELSE
         SIGN_QD =-1
      ENDIF
!#else
      !SIGN_QD = SIGN( 1._q, REAL( X, q) ) 
!#endif
   END FUNCTION SIGN_QD

!****************************************************************************
! DESCRIPTION:
!>calculates the sum of all components of a vector A
!  
!> @param[in]  A  (argument)
!****************************************************************************

   FUNCTION SUM_MPR( A  ) 
      USE prec
#ifdef qd_emulate
      USE qdmodule
#endif
      QDPREAL :: A(:)
      QDPREAL :: SUM_MPR, S
      INTEGER         :: I, N
   
      PROFILING_START('sum_mpr')
   
      S=0
      N = SIZE( A ) 
      DO I = 1, N
         S = S + A(I)
      ENDDO
      SUM_MPR = S 
   
      PROFILING_STOP('sum_mpr')
   END FUNCTION SUM_MPR 

END MODULE minimax_functions1D
