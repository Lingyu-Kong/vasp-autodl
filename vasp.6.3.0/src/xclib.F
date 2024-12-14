#include "symbol.inc"

MODULE xclib
!*******************************************************************
!
! the module xclib implements a number of common
! LDA exchange correlation functionals
! it has been moved into this file to allow inlining without 
! relying on inter file inlining (which is usually done only
! at rather high optimisation levels)
!
!*******************************************************************

  CONTAINS

!*******************************************************************
!
! In the following a couple of "standard" exchange
! and correlation functionals are implemented
! these usually use a.u. and Hartree-units 
! they calculate the energy density per electrons and the
! corresponding potential
!
! E_xc = \int eps_xc(rho(r)) rho(r) d^3 r
! v_xc = d (eps_xc(rho) rho)   / d rho 
!
!*******************************************************************


    FUNCTION ECCA(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Ceperley-Alder correlation energy as parametrised by Perdew/Zunger
! (see Phys.Rev. B23,5048 [1981], Appendix).
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      REAL(q),SAVE ::  A(2)=(/0.0622_q,0.0311_q/)
      REAL(q),SAVE ::  B(2)=(/-0.0960_q,-0.0538_q/)
!     REAL(q),SAVE ::  C(2)=(/0.0040_q,0.0014_q/)
      REAL(q),SAVE ::  C(2)=(/0.004038664055501747_q,0.001395274602717559_q/)
!     REAL(q),SAVE ::  D(2)=(/-0.0232_q,-0.0096_q/)
      REAL(q),SAVE ::  D(2)=(/-0.023264632546756681_q,-0.009602765503781227_q/)
      REAL(q),SAVE ::  G(2)=(/-0.2846_q,-0.1686_q/)
      REAL(q),SAVE ::  B1(2)=(/1.0529_q,1.3981_q/)
!     REAL(q),SAVE ::  B2(2)=(/0.3334_q,0.2611_q/)
      REAL(q),SAVE ::  B2(2)=(/0.333390000000000000_q, 0.261090000000000000_q/)
      IF (RS<=1.0_q) THEN
         RSL=LOG(RS)
         ECCA=A(IFLG)*RSL+B(IFLG)+C(IFLG)*RS*RSL+D(IFLG)*RS
      ELSE
         RSQ=SQRT(RS)
         ECCA=G(IFLG)/(1.0_q+B1(IFLG)*RSQ+B2(IFLG)*RS)
      END IF
      RETURN
    END FUNCTION ECCA

    FUNCTION VCCA(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Ceperley-Alder correlation potential as parametrised by Perdew/Zunger
! (see Phys.Rev. B23,5048 [1981], Appendix).
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION A(2),B(2),C(2),D(2),BT1(2),BT2(2)
      PARAMETER(X76=7.0_q/6.0_q, X43=4.0_q/3.0_q, &
     &  AP=0.03110_q*2.0_q, BP=-0.0480_q*2.0_q, CP=0.0020_q*2.0_q, DP=-0.0116_q*2.0_q, &
     &  AF=0.01555_q*2.0_q, BF=-0.0269_q*2.0_q, CF=0.0007_q*2.0_q, DF=-0.0048_q*2.0_q, &
     &  BP1=BP-AP/3.0_q, CP1=2.0_q*CP/3.0_q, DP1=(2.0_q*DP-CP)/3.0_q, &
     &  BF1=BF-AF/3.0_q, CF1=2.0_q*CF/3.0_q, DF1=(2.0_q*DF-CF)/3.0_q)
      SAVE A,B,C,D,BT1,BT2
      DATA A/AP,AF/,B/BP1,BF1/,C/CP1,CF1/,D/DP1,DF1/, &
     &     BT1/1.0529_q,1.3981_q/,BT2/0.3334_q,0.2611_q/
!KRESSE/FURTH---get a continous energy functional
      c(1)  = 0.004038664055501747_q * 2._q/3._q
      d(1)  =-0.023264632546756681_q * 2._q/3._q  -  0.004038664055501747_q / 3._q
      bt2(1)= 0.333390000000000000_q
      c(2)  = 0.001395274602717559_q * 2._q/3._q
      d(2)  =-0.009602765503781227_q * 2._q/3._q  -  0.001395274602717559_q / 3._q
      bt2(2)= 0.261090000000000000_q
!KRESSE/FURTH
      IF (RS<=1.0_q) THEN
         RSL=LOG(RS)
         VCCA=A(IFLG)*RSL+B(IFLG)+C(IFLG)*RS*RSL+D(IFLG)*RS
      ELSE
         RSQ=SQRT(RS)
         VCCA=ECCA(RS,IFLG)*(1.0_q+X76*BT1(IFLG)*RSQ+X43*BT2(IFLG)*RS) &
     &                     /(1.0_q+    BT1(IFLG)*RSQ+    BT2(IFLG)*RS)
      END IF
      RETURN
    END FUNCTION VCCA


    FUNCTION ECVO(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The Ceperley-Alder correlation energy as given by the Pade approximation
! technique of Vosko et al. (Can.J.Phys. 58,1200 [1980], eq.{4.4} with
! the parameters given in table 5, page 1207).
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION A(2),X02(2),B2(2),C2(2),Q2(2),XX02(2)
      PARAMETER(AP=0.0621814_q,X0P=-0.104980_q,BP=3.72744_q,CP=12.9352_q)
      PARAMETER(AF=0.0310907_q,X0F=-0.325000_q,BF=7.06042_q,CF=18.0578_q)
      PARAMETER(QP=6.1519908_q,QF=4.7309269_q)
      PARAMETER(XX0P=X0P*X0P+BP*X0P+CP,XX0F=X0F*X0F+BF*X0F+CF)
      SAVE A,X02,B2,C2,Q2,XX02
      DATA A/AP,AF/,X02/X0P,X0F/,B2/BP,BF/,C2/CP,CF/,Q2/QP,QF/, &
     &     XX02/XX0P,XX0F/
      X=SQRT(RS)
      XX=RS+B2(IFLG)*X+C2(IFLG)
      X0=X02(IFLG)
      B=B2(IFLG)
      C=C2(IFLG)
      QQ=Q2(IFLG)
      XX0=XX02(IFLG)
      ECVO=LOG((X-X0)*(X-X0)/XX)+2._q*(B+2._q*X0)/QQ*ATAN(QQ/(2._q*X+B))
      ECVO=-1._q*ECVO*B*X0/XX0+LOG(RS/XX)+2._q*B/QQ*ATAN(QQ/(2._q*X+B))
      ECVO=ECVO*A(IFLG)
      RETURN
    END FUNCTION ECVO

    FUNCTION VCVO(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The function ECVO(RS,IFLG)-RS/3.*d(ECVO(RS,IFLG))/d(RS) with function
! ECVO(RS,IFLG) given above (Ceperley-Alder potential derived from the
! approximation for ECVO of Vosko et al. discussed above).
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION A(2),X02(2),B2(2),C2(2),Q2(2),XX02(2)
      PARAMETER(AP=0.0621814_q,X0P=-0.104980_q,BP=3.72744_q,CP=12.9352_q)
      PARAMETER(AF=0.0310907_q,X0F=-0.325000_q,BF=7.06042_q,CF=18.0578_q)
      PARAMETER(QP=6.1519908_q,QF=4.7309269_q)
      PARAMETER(XX0P=X0P*X0P+BP*X0P+CP,XX0F=X0F*X0F+BF*X0F+CF)
      SAVE A,X02,B2,C2,Q2,XX02
      DATA A/AP,AF/,X02/X0P,X0F/,B2/BP,BF/,C2/CP,CF/,Q2/QP,QF/, &
     &     XX02/XX0P,XX0F/
      X=SQRT(RS)
      XX=RS+B2(IFLG)*X+C2(IFLG)
      X0=X02(IFLG)
      B=B2(IFLG)
      C=C2(IFLG)
      QQ=Q2(IFLG)
      XX0=XX02(IFLG)
      VCVO=-4._q*B*(1._q-X0*(B+2._q*X0)/XX0)/(QQ*QQ+(2._q*X+B)*(2._q*X+B))
      VCVO=VCVO-(2._q*X+B)/XX*(1._q-B*X0/XX0)-2._q*B*X0/XX0/(X-X0)+2._q/X
      VCVO=ECVO(RS,IFLG)-VCVO*A(IFLG)*X/6._q
      RETURN
    END FUNCTION VCVO
!jP---------------------------------------------------------------------------------
    FUNCTION ECVOIII(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The RPA correlation energy as given by the Pade approximation
! technique of Vosko et al. (Can.J.Phys. 58,1200 [1980], eq.{4.4} with
! the parameters given in table 5, page 1207).
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION A(2),X02(2),B2(2),C2(2),Q2(2),XX02(2)
      PARAMETER(AP=0.0621814_q,X0P=-0.409286_q,BP=13.0720_q,CP=42.7198_q)
      PARAMETER(AF=0.0310907_q,X0F=-0.743294_q,BF=20.1231_q,CF=101.578_q)
      PARAMETER(QP=0.0448999_q,QF=1.1716853_q)

      PARAMETER(XX0P=X0P*X0P+BP*X0P+CP,XX0F=X0F*X0F+BF*X0F+CF)
      SAVE A,X02,B2,C2,Q2,XX02
      DATA A/AP,AF/,X02/X0P,X0F/,B2/BP,BF/,C2/CP,CF/,Q2/QP,QF/, &
     &     XX02/XX0P,XX0F/
      X=SQRT(RS)
      XX=RS+B2(IFLG)*X+C2(IFLG)
      X0=X02(IFLG)
      B=B2(IFLG)
      C=C2(IFLG)
      QQ=Q2(IFLG)
      XX0=XX02(IFLG)
      ECVOIII=LOG((X-X0)*(X-X0)/XX)+2._q*(B+2._q*X0)/QQ*ATAN(QQ/(2._q*X+B))
      ECVOIII=-1._q*ECVOIII*B*X0/XX0+LOG(RS/XX)+2._q*B/QQ*ATAN(QQ/(2._q*X+B))
      ECVOIII=ECVOIII*A(IFLG)
      RETURN
    END FUNCTION ECVOIII

    FUNCTION VCVOIII(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The function ECVOIII(RS,IFLG)-RS/3.*d(ECVOIII(RS,IFLG))/d(RS) with function
! ECVOIII(RS,IFLG) given above (RPA potential derived from the
! approximation for ECVOIII of Vosko et al. discussed above).
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION A(2),X02(2),B2(2),C2(2),Q2(2),XX02(2)
      PARAMETER(AP=0.0621814_q,X0P=-0.409286_q,BP=13.0720_q,CP=42.7198_q)
      PARAMETER(AF=0.0310907_q,X0F=-0.743294_q,BF=20.1231_q,CF=101.578_q)
      PARAMETER(QP=0.0448999_q,QF=1.1716853_q)

      PARAMETER(XX0P=X0P*X0P+BP*X0P+CP,XX0F=X0F*X0F+BF*X0F+CF)
      SAVE A,X02,B2,C2,Q2,XX02
      DATA A/AP,AF/,X02/X0P,X0F/,B2/BP,BF/,C2/CP,CF/,Q2/QP,QF/, &
     &     XX02/XX0P,XX0F/
      X=SQRT(RS)
      XX=RS+B2(IFLG)*X+C2(IFLG)
      X0=X02(IFLG)
      B=B2(IFLG)
      C=C2(IFLG)
      QQ=Q2(IFLG)
      XX0=XX02(IFLG)
      VCVOIII=-4._q*B*(1._q-X0*(B+2._q*X0)/XX0)/(QQ*QQ+(2._q*X+B)*(2._q*X+B))
      VCVOIII=VCVOIII-(2._q*X+B)/XX*(1._q-B*X0/XX0)-2._q*B*X0/XX0/(X-X0)+2._q/X
      VCVOIII=ECVOIII(RS,IFLG)-VCVOIII*A(IFLG)*X/6._q
      RETURN
    END FUNCTION VCVOIII
!jP---------------------------------------------------------------------------------

    FUNCTION ECGL(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Gunnarson-Lundqvist correlation energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION C(2),R(2)
      PARAMETER(THIRD=1.0_q/3.0_q)
      SAVE C,R
      DATA C/0.0666_q,0.0406_q/,R/11.4_q,15.9_q/
      X=RS/R(IFLG)
      ECGL=-C(IFLG)*((1._q+X**3)*LOG(1._q+1._q/X)-THIRD+X*(0.5_q-X))
      RETURN
    END FUNCTION ECGL

    FUNCTION VCGL(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Gunnarson-Lundqvist correlation potential:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION C(2),R(2)
      SAVE C,R
      DATA C/0.0666_q,0.0406_q/,R/11.4_q,15.9_q/
      VCGL=-C(IFLG)*LOG(1._q+R(IFLG)/RS)
      RETURN
    END FUNCTION VCGL

    FUNCTION ECHL(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Hedin-Lundqvist correlation energy (J.Phys. C4,2064 [1971]):
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION C(2),R(2)
      PARAMETER(THIRD=1.0_q/3.0_q)
      SAVE C,R
      DATA C/0.045_q,0.0225_q/,R/21.0_q,52.917_q/
      X=RS/R(IFLG)
      ECHL=-C(IFLG)*((1._q+X**3)*LOG(1._q+1._q/X)-THIRD+X*(0.5_q-X))
      RETURN
    END FUNCTION ECHL

    FUNCTION VCHL(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Hedin-Lundqvist correlation potential (J.Phys. C4,2064 [1971]):
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION C(2),R(2)
      SAVE C,R
      DATA C/0.045_q,0.0225_q/,R/21.0_q,52.917_q/
      VCHL=-C(IFLG)*LOG(1._q+R(IFLG)/RS)
      RETURN
    END FUNCTION VCHL

    FUNCTION ECBH(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Barth-Hedin correlation energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION C(2),R(2)
      PARAMETER(THIRD=1.0_q/3.0_q)
      SAVE C,R
      DATA C/0.0504_q,0.0254_q/,R/30._q,75._q/
      X=RS/R(IFLG)
      ECBH=-C(IFLG)*((1._q+X**3)*LOG(1._q+1._q/X)-THIRD+X*(0.5_q-X))
      RETURN
    END FUNCTION ECBH

    FUNCTION VCBH(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Barth-Hedin correlation potential:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION C(2),R(2)
      SAVE C,R
      DATA C/0.0504_q,0.0254_q/,R/30._q,75._q/
      VCBH=-C(IFLG)*LOG(1._q+R(IFLG)/RS)
      RETURN
    END FUNCTION VCBH

      FUNCTION ECWI(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Wigner correlation energy (hopefully correct?):
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results --> warning: equals paramagnetic result
      DIMENSION CX(2),C(2),R(2)
      SAVE CX,C,R
      DATA CX/0.9163305865663_q,1.1545041946774_q/
      DATA C/7.8_q,7.8_q/,R/0.88_q,0.88_q/
      X=CX(IFLG)/R(IFLG)*(1._q+C(IFLG)/RS)
      ECWI=-CX(IFLG)/X/RS
      RETURN
    END FUNCTION ECWI

    FUNCTION VCWI(RS,IFLG)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Wigner correlation potential (hopefully correct?):
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results --> warning: equals paramagnetic result
      DIMENSION CX(2),C(2),R(2)
      SAVE CX,C,R
      DATA CX/0.9163305865663_q,1.1545041946774_q/
      DATA C/7.8_q,7.8_q/,R/0.88_q,0.88_q/
      X1=C(IFLG)/RS
      X2=1._q+X1
      X3=CX(IFLG)/R(IFLG)*X2
      B=1._q+1._q/X3
      F=1._q-X1/X2/(1._q+X3)
      F=1._q+F/3._q
      E=-CX(IFLG)*B/RS
      VCWI=(4._q/3._q-F*B)*CX(IFLG)/RS
      RETURN
    END FUNCTION VCWI

    FUNCTION EX(RS,IFLG,TREL)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Exchange energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
!     TREL  :  Relativistic correction or not (for details see
!              J.Phys. C12,2977(1979) )
      LOGICAL TREL
      DIMENSION CX(2)
      SAVE CX,CBETA
      DATA CX/0.9163305865663_q,1.1545041946774_q/,CBETA/0.0140_q/
      EX=-CX(IFLG)/RS
      IF (TREL) THEN
         B=CBETA/RS
         F=LOG(B+(SQRT(1+B*B)))/(B*B)
         F=(SQRT(1+B*B)/B)-F
!jF: the expression given above becomes numerically extremely instable for
!    very small values of B (small difference of two large numbers divided
!    by small number = noise) therefore use following for reasons of safety:
         IF (B.LT.1.E-5_q) F=B*(2._q/3._q-0.7_q*B*B)
         EX=(1._q-1.5_q*F*F)*EX
      END IF
      RETURN
    END FUNCTION EX

    FUNCTION VX(RS,IFLG,TREL)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Exchange potential:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
!     TREL  :  Relativistic correction or not (for details see
!              J.Phys. C12,2977(1979) )
      LOGICAL TREL
      DIMENSION CX(2)
      SAVE CX,CBETA
      DATA CX/1.2217741154217_q,1.5393389262365_q/,CBETA/0.0140_q/

      VX=-CX(IFLG)/RS
      IF (TREL) THEN
! Warning error in the paper of Bachelet et al. !!
         B=CBETA/RS
         F=LOG(B+(SQRT(1+B*B)))/B/SQRT(1+(B*B))
!        F=LOG(B+(SQRT(1+B*B)))/B/(1+(B*B))
!jF: potentially the expression given above becomes numerically instable for
!    very small values of B, therefore use following for reasons of safety:
         IF (B.LT.1.E-5_q) F=1._q-B*B*(2._q/3._q-B*B*31._q/30._q)
         VX=(-0.5_q+1.5_q*F)*VX
      END IF
      RETURN
    END FUNCTION VX
      
    FUNCTION FZ0(ZETA)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! Interpolation function between paramagnetic and ferromagnetic results
! for exchange energy and exchange potential. The parameter ZETA is:
! ZETA=(RHO[upspin] - RHO[downspin]) / (RHO[upspin] + RHO[downspin]).
      PARAMETER(C43=4._q/3._q,FAC=1.92366105093153632_q)
      Z=ABS(ZETA)
      IF (Z==0._q) THEN
         FZ0=0._q
      ELSE IF (Z>=1._q) THEN
         FZ0=1._q
      ELSE
         FZ0=(((1._q+Z)**C43)+((1._q-Z)**C43)-2._q)*FAC
      END IF
      RETURN
    END FUNCTION FZ0

    FUNCTION FZ1(ZETA)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The derivative dFZ0(ZETA)/d(ZETA), FZ0(ZETA) given above.
      PARAMETER(C13=1._q/3._q,FAC=2.56488140124204843_q)
      Z=ABS(ZETA)
      IF (Z==0._q) THEN
         FZ1=0._q
      ELSE IF (Z>=1._q) THEN
         FZ1=(2._q**C13)*FAC
      ELSE
         FZ1=(((1._q+Z)**C13)-((1._q-Z))**C13)*FAC
      END IF
      RETURN
    END FUNCTION FZ1

    FUNCTION ALPHA0(RS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The spin stiffness [d(EXC(RS,ZETA))**2/d**2(ZETA) |ZETA=0] as given by
! Vosko et al. (Can.J.Phys. 58,1200 [1980], eq.{4.4} with the parameters
! given on page 1209 [fitting to low-density values using eq.{4.7.}]).
! Warning: the values are multiplied by 1./(d(FZ0(ZETA))**2/d**2(ZETA))
!          at ZETA=0, FZ0(ZETA) given above.
      PARAMETER(X0=-0.0047584_q,B=1.13107_q,C=13.0045_q,QQ=7.12311_q)
      PARAMETER(XX0=X0*X0+B*X0+C,A=-0.019751631321681_q)
      X=SQRT(RS)
      XX=RS+B*X+C
      ALPHA0=LOG((X-X0)*(X-X0)/XX)+2._q*(B+2._q*X0)/QQ*ATAN(QQ/(2._q*X+B))
      ALPHA0=-1._q*ALPHA0*B*X0/XX0+LOG(RS/XX)+2._q*B/QQ*ATAN(QQ/(2._q*X+B))
      ALPHA0=ALPHA0*A
      RETURN
    END FUNCTION ALPHA0

    FUNCTION ALPHA1(RS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The function ALPHA0(RS)-RS/3.*dALPHA0(RS)/dRS, ALPHA0(RS) given above.
! Warning: the values are multiplied by 1./(d(FZ0(ZETA))**2/d**2(ZETA))
!          at ZETA=0, FZ0(ZETA) given above.
      PARAMETER(X0=-.0047584_q,B=1.13107_q,C=13.0045_q,QQ=7.12311_q)
      PARAMETER(XX0=X0*X0+B*X0+C,A=-0.0032919385536135_q)
      X=SQRT(RS)
      XX=RS+B*X+C
      ALPHA1=-4._q*B/QQ*(1._q-X0*(B+2._q*X0)/XX0)/(QQ*QQ+(2._q*X+B)*(2._q*X+B))
      ALPHA1=ALPHA1-(2._q*X+B)/XX*(1._q-B*X0/XX0)-2._q*B*X0/XX0/(X-X0)+2._q/X
      ALPHA1=ALPHA0(RS)-ALPHA1*A*X
      RETURN
    END FUNCTION ALPHA1
!----------------------------------------------------------------------
    FUNCTION ALPHA0_III(RS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The spin stiffness [d(EXC(RS,ZETA))**2/d**2(ZETA) |ZETA=0] as given by
! Vosko et al. (Can.J.Phys. 58,1200 [1980], eq.{4.4} with the parameters
! given on page 1209 [fitting to low-density values using eq.{4.7.}]).
! Warning: the values are multiplied by 1./(d(FZ0(ZETA))**2/d**2(ZETA))
!          at ZETA=0, FZ0(ZETA) given above.
      PARAMETER(X0=-0.228344_q,B=1.06835_q,C=11.4813_q,QQ=6.69207_q)
      PARAMETER(XX0=X0*X0+B*X0+C,A=-0.019751631321681_q)
      X=SQRT(RS)
      XX=RS+B*X+C
      ALPHA0_III=LOG((X-X0)*(X-X0)/XX)+2._q*(B+2._q*X0)/QQ*ATAN(QQ/(2._q*X+B))
      ALPHA0_III=-1._q*ALPHA0_III*B*X0/XX0+LOG(RS/XX)+2._q*B/QQ*ATAN(QQ/(2._q*X+B))
      ALPHA0_III=ALPHA0_III*A
      RETURN
    END FUNCTION ALPHA0_III

    FUNCTION ALPHA1_III(RS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! The function ALPHA0(RS)-RS/3.*dALPHA0(RS)/dRS, ALPHA0(RS) given above.
! Warning: the values are multiplied by 1./(d(FZ0(ZETA))**2/d**2(ZETA))
!          at ZETA=0, FZ0(ZETA) given above.
      PARAMETER(X0=-0.228344_q,B=1.06835_q,C=11.4813_q,QQ=6.69207_q)
      PARAMETER(XX0=X0*X0+B*X0+C,A=-0.0032919385536135_q)
      X=SQRT(RS)
      XX=RS+B*X+C
      ALPHA1_III=-4._q*B/QQ*(1._q-X0*(B+2._q*X0)/XX0)/(QQ*QQ+(2._q*X+B)*(2._q*X+B))
      ALPHA1_III=ALPHA1_III-(2._q*X+B)/XX*(1._q-B*X0/XX0)-2._q*B*X0/XX0/(X-X0)+2._q/X
      ALPHA1_III=ALPHA0_III(RS)-ALPHA1_III*A*X
      RETURN
    END FUNCTION ALPHA1_III
!----------------------------------------------------------------------

!----------------------------------------------------------------------
    SUBROUTINE CORPBE_LDA(RS,ZET,EC,VCUP,VCDN)
!----------------------------------------------------------------------
!  LDA part of the official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!  output: ec=lsd correlation energy from [a]
!        : vcup=lsd up correlation potential
!        : vcdn=lsd dn correlation potential
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, 
!     {\sl Generalized gradient approximation made simple}, sub.
!     to Phys. Rev.Lett. May 1996.
! [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
!     construction of a generalized gradient approximation:  The PW91
!     density functional}, submitted to Phys. Rev. B, Feb. 1996.
! [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      logical lgga
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet=0.06672455060314922_q,delt=bet/gamma)
      parameter(eta=1.e-12_q)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor_xc(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS)
      CALL gcor_xc(0.01554535_q,0.20548_q,14.1189_q,6.1977_q,3.3662_q, &
     &    0.62517_q,rtRS,EP,EPRS)
      CALL gcor_xc(0.0168869_q,0.11125_q,10.357_q,3.6231_q,0.88026_q, &
     &    0.49671_q,rtRS,ALFM,ALFRSM)
      ALFC = -ALFM
      Z4 = ZET**4
      F=((1._q+ZET)**THRD4+(1._q-ZET)**THRD4-2._q)/GAM
      EC = EU*(1._q-F*Z4)+EP*F*Z4-ALFM*F*(1._q-Z4)/FZZ
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LSD potential from [c](A1)
! ECRS = dEc/drs [c](A2)
! ECZET=dEc/dzeta [c](A3)
! FZ = dF/dzeta [c](A4)
      ECRS = EURS*(1._q-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1._q-Z4)/FZZ
      FZ = THRD4*((1._q+ZET)**THRD-(1._q-ZET)**THRD)/GAM
      ECZET = 4._q*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
     &        -(1._q-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3._q-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units
      EC=EC*2
      VCUP=VCUP*2
      VCDN=VCDN*2
      RETURN
    END SUBROUTINE CORPBE_LDA
    
    
!----------------------------------------------------------------------
! Judith Harl
! functionals for range-separated ACFDT (LDA - short range RPA):
! a bit akward at the moment since the range separation parameter
! is hard coded for now. Hopefully this will change ...
!----------------------------------------------------------------------

    SUBROUTINE CORPBE_LDA_RPA(RS,ZET,EC,VCUP,VCDN)
! RPA correlation energy for jellium 
! Perdew, Wang, prB 45 13244 (1992)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      logical lgga
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet=0.06672455060314922_q,delt=bet/gamma)
      parameter(eta=1.e-12_q)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor_xc_rpa(0.031091_q,0.082477_q,5.1486_q,1.6483_q,0.23647_q, &
     &    0.20614_q,rtrs,EU,EURS)
      CALL gcor_xc_rpa(0.015545_q,0.035374_q,6.4869_q,1.3083_q,0.15180_q, &
     &    0.082349_q,rtRS,EP,EPRS)
      CALL gcor_xc(0.016887_q,0.028829_q,10.357_q,3.6231_q,0.47990_q, &
     &    0.12279_q,rtRS,ALFM,ALFRSM)
      ALFC = -ALFM
      Z4 = ZET**4
      F=((1._q+ZET)**THRD4+(1._q-ZET)**THRD4-2._q)/GAM
      EC = EU*(1._q-F*Z4)+EP*F*Z4-ALFM*F*(1._q-Z4)/FZZ
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LSD potential from [c](A1)
! ECRS = dEc/drs [c](A2)
! ECZET=dEc/dzeta [c](A3)
! FZ = dF/dzeta [c](A4)
      ECRS = EURS*(1._q-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1._q-Z4)/FZZ
      FZ = THRD4*((1._q+ZET)**THRD-(1._q-ZET)**THRD)/GAM
      ECZET = 4._q*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
     &        -(1._q-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3._q-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units
      EC=EC*2
      VCUP=VCUP*2
      VCDN=VCDN*2
      RETURN
    END SUBROUTINE CORPBE_LDA_RPA


!  jH-new This routine is used to evaluate the RPA+ correction term, 
!  jH-new which is E^{LDA,QMC}_{c}-E^{LDA,RPA}_{c}

    SUBROUTINE CORPBE_LDA_RPA_PLUS(RS,ZET,EC,VCUP,VCDN)
! RPA correlation energy for jellium 
! Perdew, Wang, prB 45 13244 (1992)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      logical lgga
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet=0.06672455060314922_q,delt=bet/gamma)
      parameter(eta=1.e-12_q)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor_xc_rpa(0.031091_q,0.082477_q,5.1486_q,1.6483_q,0.23647_q, &
     &    0.20614_q,rtrs,EURPA,EURSRPA)
      CALL gcor_xc_rpa(0.015545_q,0.035374_q,6.4869_q,1.3083_q,0.15180_q, &
     &    0.082349_q,rtRS,EPRPA,EPRSRPA)
      CALL gcor_xc(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EUQM,EURSQM)
      CALL gcor_xc(0.01554535_q,0.20548_q,14.1189_q,6.1977_q,3.3662_q, &
     &    0.62517_q,rtRS,EPQM,EPRSQM)
      EU=EUQM-EURPA
      EP=EPQM-EPRPA
      EURS=EURSQM-EURSRPA
      EPRS=EPRSQM-EPRSRPA
! jH rpa spin stiffness
      CALL gcor_xc(0.016887_q,0.028829_q,10.357_q,3.6231_q,0.47990_q, &
     &    0.12279_q,rtRS,ALFMRPA,ALFRSMRPA)
! jH QMC spin stiffness
      CALL gcor_xc(0.0168869_q,0.11125_q,10.357_q,3.6231_q,0.88026_q, &
     &    0.49671_q,rtRS,ALFMQM,ALFRSMQM)     
      ALFCRPA = -ALFMRPA
      ALFCQM  = -ALFMQM
      Z4 = ZET**4
      F=((1._q+ZET)**THRD4+(1._q-ZET)**THRD4-2._q)/GAM
      ECRPA = EURPA*(1._q-F*Z4)+EPRPA*F*Z4-ALFMRPA*F*(1._q-Z4)/FZZ
      ECQM  =  EUQM*(1._q-F*Z4)+EPQM*F*Z4-ALFMQM*F*(1._q-Z4)/FZZ
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LSD potential from [c](A1)
! ECRS = dEc/drs [c](A2)
! ECZET=dEc/dzeta [c](A3)
! FZ = dF/dzeta [c](A4)
! jH Aufpassen: Potential 
      ECRSRPA = EURSRPA*(1._q-F*Z4)+EPRSRPA*F*Z4-ALFRSMRPA*F*(1._q-Z4)/FZZ
      ECRSQM  =  EURSQM*(1._q-F*Z4)+EPRSQM*F*Z4-ALFRSMQM*F*(1._q-Z4)/FZZ
      FZ = THRD4*((1._q+ZET)**THRD-(1._q-ZET)**THRD)/GAM
      ECZETRPA = 4._q*(ZET**3)*F*(EPRPA-EURPA+ALFMRPA/FZZ)+FZ*(Z4*EPRPA-Z4*EURPA &
     &        -(1._q-Z4)*ALFMRPA/FZZ)
      ECZETQM = 4._q*(ZET**3)*F*(EPQM-EUQM+ALFMQM/FZZ)+FZ*(Z4*EPQM-Z4*EUQM &
     &        -(1._q-Z4)*ALFMQM/FZZ) 
      COMMRPA = ECRPA -RS*ECRSRPA/3._q-ZET*ECZETRPA
      COMMQM  = ECQM  -RS*ECRSQM/3._q-ZET*ECZETQM
      VCUP = COMMQM + ECZETQM - (COMMRPA + ECZETRPA)
      VCDN = COMMQM - ECZETQM - (COMMRPA - ECZETRPA)

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units
!      EC=EC*2
       EC=(ECQM-ECRPA)*2
       VCUP=VCUP*2
       VCDN=VCDN*2
      RETURN
    END SUBROUTINE CORPBE_LDA_RPA_PLUS

    SUBROUTINE CORPBE_LDA_SR_0_15au(RS,ZET,EC,VCUP,VCDN)
!This routines gives a correlation energy functional defined by
!    Ec,SR[n0]:=   Ec,MC[n0]    -   Ec,RPA,LR[n0] 
!and   Ec,RPA,LR[n0] = Ec,RPA[n0]-Ec,RPA,sr[n0]
!where Ec,MC is the Perdew-Wang parametrization of the
!Monte Carlo correlation energies [1], Ec,RPA the P-W 
!parametrization of the RPA energy [1], and Ec,RPA,sr the 
!v = (1/r)Erfc(-mu r) short-range ACFDT-RPA energies for
!jellium calculated by Judith 
! [1] Perdew,Wang 1992 prB 45, 13244
! 
!   Just for the paramagnetic case
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      logical lgga
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet=0.06672455060314922_q,delt=bet/gamma)
      parameter(eta=1.e-12_q)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor_xc_rpa(0.031091_q,0.082477_q,5.1486_q,1.6483_q,0.23647_q, &
     &    0.20614_q,rtrs,EU,EURS)
      ERPATOT = EU
      ERPATOTRS = EURS
       CALL gcor_xc(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS)
      EMC     = EU
      EMCRS   = EURS 
       CALL gcor_xc_2(0.031091_q,-0.0089687_q,5.1486_q,1.3952503_q, & 
     & 0.2284487_q,0.0965737_q,1.6078169_q,rtrs,EU,EURS)
      ERPASR = EU
      ERPASRRS = EURS
! SR correlation energy so definiert das HEG limit
      EU = EMC - (ERPATOT - ERPASR)
      EURS = EMCRS -  (ERPATOTRS - ERPASRRS)

! jH thats enough

      EC=EU
      ECRS=EURS
      COMM=EC-S*ECRS/3._q 
      VCUP = COMM    
      VCDN = COMM

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units
      EC=EC*2
      VCUP=VCUP*2
      VCDN=VCDN*2
      RETURN
    END SUBROUTINE CORPBE_LDA_SR_0_15au

    SUBROUTINE CORPBE_LDA_SR_0_5A(RS,ZET,EC,VCUP,VCDN)
!This routines gives a correlation energy functional defined by
!    Ec,SR[n0]:=   Ec,MC[n0]    -   Ec,RPA,LR[n0] 
!and   Ec,RPA,LR[n0] = Ec,RPA[n0]-Ec,RPA,sr[n0]
!where Ec,MC is the Perdew-Wang parametrization of the
!Monte Carlo correlation energies [1], Ec,RPA the P-W 
!parametrization of the RPA energy [1], and Ec,RPA,sr the 
!v = (1/r)Erfc(-mu r) short-range ACFDT-RPA energies for
!jellium calculated by Judith 
! [1] Perdew,Wang 1992 prB 45, 13244
! 
!   Just for the paramagnetic case
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      logical lgga
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet=0.06672455060314922_q,delt=bet/gamma)
      parameter(eta=1.e-12_q)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor_xc_rpa(0.031091_q,0.082477_q,5.1486_q,1.6483_q,0.23647_q, &
     &    0.20614_q,rtrs,EU,EURS)
      ERPATOT = EU
      ERPATOTRS = EURS
       CALL gcor_xc(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS)
      EMC     = EU
      EMCRS     = EURS
       CALL gcor_xc_2(0.031091_q,-0.0106268_q,5.1486_q,1.4094821_q, & 
     & 1.2291950_q,0.2579335_q,1.8066972_q,rtrs,EU,EURS)
      ERPASR = EU
      ERPASRRS = EURS
      EU = EMC - (ERPATOT - ERPASR)
      EURS = EMCRS -  (ERPATOTRS - ERPASRRS)

! jH thats enough

      EC=EU
      ECRS=EURS
      COMM=EC-S*ECRS/3._q
      VCUP = COMM
      VCDN = COMM

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units
      EC=EC*2
      VCUP=VCUP*2
      VCDN=VCDN*2
      RETURN
    END SUBROUTINE CORPBE_LDA_SR_0_5A

    SUBROUTINE CORPBE_LDA_SR_1_0A(RS,ZET,EC,VCUP,VCDN)
!This routines gives a correlation energy functional defined by
!    Ec,SR[n0]:=   Ec,MC[n0]    -   Ec,RPA,LR[n0] 
!and   Ec,RPA,LR[n0] = Ec,RPA[n0]-Ec,RPA,sr[n0]
!where Ec,MC is the Perdew-Wang parametrization of the
!Monte Carlo correlation energies [1], Ec,RPA the P-W 
!parametrization of the RPA energy [1], and Ec,RPA,sr the 
!v = (1/r)Erfc(-mu r) short-range ACFDT-RPA energies for
!jellium calculated by Judith 
! [1] Perdew,Wang 1992 prB 45, 13244
! 
!   Just for the paramagnetic case
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      logical lgga
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet=0.06672455060314922_q,delt=bet/gamma)
      parameter(eta=1.e-12_q)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor_xc_rpa(0.031091_q,0.082477_q,5.1486_q,1.6483_q,0.23647_q, &
     &    0.20614_q,rtrs,EU,EURS)
      ERPATOT = EU
      ERPATOTRS = EURS
       CALL gcor_xc(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS)
      EMC     = EU
      EMCRS     = EURS
       CALL gcor_xc_2(0.031091_q,-0.0119450_q,5.1486_q,1.5999830_q, & 
     & 4.7029103_q,1.3622380_q,1.9748249_q,rtrs,EU,EURS)
      ERPASR = EU
      ERPASRRS = EURS
      EU = EMC -(ERPATOT - ERPASR)
      EURS = EMCRS -(ERPATOTRS - ERPASRRS)

! jH thats enough

      EC=EU
      ECRS=EURS
      COMM=EC-S*ECRS/3._q
      VCUP = COMM
      VCDN = COMM

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units
      EC=EC*2
      VCUP=VCUP*2
      VCDN=VCDN*2
      RETURN
    END SUBROUTINE CORPBE_LDA_SR_1_0A

    SUBROUTINE CORPBE_LDA_SR_2_0A(RS,ZET,EC,VCUP,VCDN)
!This routines gives a correlation energy functional defined by
!    Ec,SR[n0]:=   Ec,MC[n0]    -   Ec,RPA,LR[n0] 
!and   Ec,RPA,LR[n0] = Ec,RPA[n0]-Ec,RPA,sr[n0]
!where Ec,MC is the Perdew-Wang parametrization of the
!Monte Carlo correlation energies [1], Ec,RPA the P-W 
!parametrization of the RPA energy [1], and Ec,RPA,sr the 
!v = (1/r)Erfc(-mu r) short-range ACFDT-RPA energies for
!jellium calculated by Judith 
! [1] Perdew,Wang 1992 prB 45, 13244
! 
!   Just for the paramagnetic case
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      logical lgga
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet=0.06672455060314922_q,delt=bet/gamma)
      parameter(eta=1.e-12_q)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor_xc_rpa(0.031091_q,0.082477_q,5.1486_q,1.6483_q,0.23647_q, &
     &    0.20614_q,rtrs,EU,EURS)
      ERPATOT = EU
      ERPATOTRS = EURS
       CALL gcor_xc(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS)
      EMC     = EU
      EMCRS     = EURS
       CALL gcor_xc_2(0.031091_q,0.0127716_q,5.1486_q,2.5109626_q, & 
     & 15.3666137_q,10.2843432_q,2.1207238_q,rtrs,EU,EURS)
      ERPASR = EU
      ERPASRRS = EURS
      EU = EMC - (ERPATOT - ERPASR)
      EURS = EMCRS - (ERPATOTRS - ERPASRRS) 

! jH thats enough

      EC=EU
      ECRS=EURS
      COMM=EC-S*ECRS/3._q
      VCUP = COMM
      VCDN = COMM

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units
      EC=EC*2
      VCUP=VCUP*2
      VCDN=VCDN*2
      RETURN
    END SUBROUTINE CORPBE_LDA_SR_2_0A


!----------------------------------------------------------------------
    FUNCTION PBE_ALPHA(RS)
!----------------------------------------------------------------------
!  Alpha (spin stiffness)
!  from the official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!  spin stiffness 
!----------------------------------------------------------------------
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))

      rtrs=dsqrt(rs)
      CALL gcor_xc(0.0168869_q,0.11125_q,10.357_q,3.6231_q,0.88026_q, &
     &    0.49671_q,rtRS,ALFM,ALFRSM)

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units

      PBE_ALPHA = -ALFM/FZZ*2
      RETURN
    END FUNCTION PBE_ALPHA

! jH RPA_ALPHA new routine for RPA spin stiffness 
!----------------------------------------------------------------------
    FUNCTION RPA_ALPHA(RS)
!----------------------------------------------------------------------
!  Alpha (spin stiffness)
!  from the official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!  spin stiffness
!  jH this parametrisation has been already impl in CORPBE_LDA_RPA
!  jH parametrisation from Perdew, Wang, prB 45, 13244 (1992)
!----------------------------------------------------------------------
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))

      rtrs=dsqrt(rs)
      CALL gcor_xc(0.016887_q,0.028829_q,10.357_q,3.6231_q,0.47990_q, &
     &    0.12279_q,rtRS,ALFM,ALFRSM)

      ! the convention of the subroutines in this package
      ! is to return Rydberg energy units

      RPA_ALPHA = -ALFM/FZZ*2
      RETURN
    END FUNCTION RPA_ALPHA


!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
    SUBROUTINE GCOR_XC(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
! slimmed down version of GCOR used in PW91 routines, to interpolate
! LSD correlation energy, as given by (10) of
! J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
! K. Burke, May 11, 1996.
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      Q0 = -2._q*A*(1._q+A1*rtrs*rtrs)
      Q1 = 2._q*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1._q+1._q/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2._q*B2+rtrs*(3._q*B3+4._q*B4*rtrs))
      GGRS = -2._q*A*A1*Q2-Q0*Q3/(Q1*(1._q+Q1))
      RETURN
    END SUBROUTINE GCOR_XC

!----------------------------------------------------------------------
!#########################JUDITH###################################
!----------------------------------------------------------------------

    SUBROUTINE GCOR_XC_RPA(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
! This is the rs interpolation formula of Perdew/Wang Phys.Rev.B 45 (1992) 13244
! for the RPA case (p = 0.75). See FORMULA (10) in article 
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      Q0 = -2._q*A*(1._q+A1*rtrs*rtrs)
      Q1 = 2._q*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*sqrt(rtrs))))
      Q2 = DLOG(1._q+1._q/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2._q*B2+rtrs*(3._q*B3+7._q/2._q*B4*sqrt(rtrs)))
      GGRS = -2._q*A*A1*Q2-Q0*Q3/(Q1*(1._q+Q1))
      RETURN
      END SUBROUTINE GCOR_XC_RPA

!----------------------------------------------------------------------
!###########################JUDITH################################
!----------------------------------------------------------------------

    SUBROUTINE GCOR_XC_2(A,A1,B1,B2,B3,B4,p,rtrs,GG,GGRS)
! This is the rs interpolation formula of Perdew/Wang Phys.Rev.B 45 (1992) 13244
! for p = 2. See FORMULA (10) in article. This interpolation is used to
! parametrsize the short-range for v = (1/r)*Erfc(-mu r) RPA jellium energie --> Judith
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      Q0 = -2._q*A*(1._q+A1*rtrs*rtrs)
      Q1 = 2._q*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs**(2._q*p-1._q))))
      Q2 = DLOG(1._q+1._q/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2._q*B2+rtrs*(3._q*B3+(2._q*p+2._q)*B4*rtrs**(2._q*p-1._q)))
      GGRS = -2._q*A*A1*Q2-Q0*Q3/(Q1*(1._q+Q1))
      RETURN
      END SUBROUTINE GCOR_XC_2

!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
! Iann Gerber 18/11/04

!=====================================================================
! short range screened exchange functional
! EX_SR computes exchange-part with modified kernel
! interpolation formula of A. Savin and J. Toulouse 
!=====================================================================

    
    FUNCTION EX_SR(RS,RMU,IFLG)
!$ACC ROUTINE SEQ
#ifdef _OPENACC
      use openacc
#endif
      USE prec
      USE constant
      USE c2f_interface, ONLY : ERRF
      IMPLICIT REAL(q) (A-H,O-Z)
      REAL(q) :: RS        ! Wigner Seitz radius in a.u.
      REAL(q) :: RMU       ! inverse real space cutoff in a.u.
      INTEGER :: IFLG      ! paramagnetic (IFLG=1) or ferromagnetic (IFLG=2)
! Exchange energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION CX(2)
      SAVE CX
      DATA CX/0.9163305865663_q,1.1545041946774_q/
      EX_=-CX(IFLG)/RS

! another IG bug: IFLAG replaced by IFLG
      IF (IFLG==2) THEN
         QFAC=(6._q*PI*PI)**(1._q/3._q)
      ELSE
         QFAC=(3._q*PI*PI)**(1._q/3._q)
      ENDIF

      RHO = 3/(4._q*PI)/RS**3
      QF=QFAC*(RHO)**(1._q/3._q)
      A=RMU/2._q/QF
!      FRAC=(1 - (8._q/3._q)*A*( SQRT(PI)*ERRF(1/(2*A)) + (2*A-4*A*A*A)*EXP(-1/(4*A*A)) &
!     &        -3*A + 4*A*A*A ))
!      EX_SR=EX_*FRAC
! IG Simple formula replaced by 

! Test on the value of A
      IF (A < 1E-9_q) THEN
! Limit for small A
         EX_SR=EX_
      ELSEIF (A <= 100._q) THEN
! Intermediate Values of A
#ifdef _OPENACC
!(sm): only when compiling GPU code this construct will use standard ERF which has a GPU
!      implementation instead of the custom CPU implementation ERRF
      IF(ACC_ON_DEVICE(ACC_DEVICE_NOT_HOST)) THEN
         FRAC=(1.0_q - (8._q/3._q)*A*( SQRT(PI)*ERF(1._q/(2._q*A)) + &
     &        (2._q*A-4._q*A*A*A)*EXP(-1._q/(4._q*A*A))-3._q*A + 4._q*A*A*A ))
      ELSE
#endif
         FRAC=(1.0_q - (8._q/3._q)*A*( SQRT(PI)*ERRF(1._q/(2._q*A)) + &
     &        (2._q*A-4._q*A*A*A)*EXP(-1._q/(4._q*A*A))-3._q*A + 4._q*A*A*A ))
#ifdef _OPENACC
      ENDIF
#endif
         EX_SR=EX_*FRAC
      ELSEIF (A <= 1.E+9_q) THEN
! Development for large A
         FRAC=1._q/36._q/A/A
         EX_SR=EX_*FRAC
      ELSE
         EX_SR=0._q
      ENDIF
!      B=EXP(-1._q/(4._q*A*A))-1._q
!      C=2._q*A*A*B+1._q/2._q
!      FRAC=1.0_q - 8._q/3._q*A*(SQRT(PI)*ERRF(1._q/(2._q*A))+2._q*A*(B-C))
!      EX_SR=EX_*FRAC
      RETURN
    END FUNCTION EX_SR
      

    FUNCTION VX_SR(RS,RMU,IFLG)
!$ACC ROUTINE SEQ
      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
      REAL(q) :: RS        ! Wigner Seitz radius in a.u.
      REAL(q) :: RMU       ! inverse real space cutoff in a.u.
      INTEGER :: IFLG      ! paramagnetic (IFLG=1) or ferromagnetic (IFLG=2)
! Exchange energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION CX(2)
      SAVE CX
      DATA CX/0.9163305865663_q,1.1545041946774_q/
      VX_=-CX(IFLG)/RS

! another IG bug: IFLAG replaced by IFLG
      IF (IFLG==2) THEN
         QFAC=(6._q*PI*PI)**(1._q/3._q)
      ELSE
         QFAC=(3._q*PI*PI)**(1._q/3._q)
      ENDIF

      RHO = 3/(4._q*PI)/RS**3
      QF=QFAC*(RHO)**(1._q/3._q)
      A=RMU/2._q/QF

!      FRAC= 1._q - 32*A*A*A*A*EXP(-1/(4*A*A)) -8._q*A*A + (64._q/3._q)*A*A*A*A
!      VX_SR=EX_SR(RS,RMU,IFLG)+(1._q/3._q)*VX_*FRAC
! IG Simple formula replaced by 

! Test on the value of A
      IF (A < 1E-9_q) THEN
! Limit for small A
         VX_SR=(4._q/3._q)*VX_
      ELSEIF (A <= 100._q) THEN
! Intermediate Values of A
         FRAC= (1._q - 32._q*A*A*A*A*(EXP(-1._q/(4._q*A*A))-1._q) -8._q*A*A)
         VX_SR=EX_SR(RS,RMU,IFLG)+(1._q/3._q)*VX_*FRAC
      ELSEIF (A <= 1.E+9_q) THEN
! Development for large A
         VX_SR=-RMU/48._q/PI/A/A/A
      ELSE
         VX_SR=0._q
      ENDIF
     RETURN
   END FUNCTION VX_SR

   FUNCTION EC_SR(RS,RMU,EC) ! No spin polarization
      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
      REAL(q) :: RS        ! Wigner Seitz radius in a.u.
      REAL(q) :: RMU       ! inverse real space cutoff in a.u.

      SAVE U1,U2,V1,A,BET,GAM
      DATA U1/1.0270741452992294_q/
      DATA U2/-0.230160617208092_q/
      DATA V1/0.6196884832404359_q/
      DATA A/3.2581_q/
      DATA BET/163.44_q/
      DATA GAM/4.7125_q/

      D=32._q/3._q/PI
! changed to a.u. instead of A, iG please check
      RHO = 3/(4._q*PI)/RS**3

!     G(0) from Burke, Perdew & Ernzerhof
      GRS32=(GAM+RS)**3._q/2._q
      GRS12=SQRT(GAM+RS)
      G0=D*(GRS32+BET)*EXP(-A*GRS12)

      C1N=U1*RS + U2*RS*RS
      C1D=1._q+V1*RS
      C1=C1N/C1D
      C2D=0.5_q*PI*RHO*RHO*(G0-0.5_q)
      C2=EC/C2D

      SCALE=1._q + C1*RMU + C2*RMU*RMU

      EC_SR=EC/SCALE
      RETURN
    END FUNCTION EC_SR

    FUNCTION VC_SR(RS,RMU,VC,EC) ! No spin polarization
      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
      REAL(q) :: RS        ! Wigner Seitz radius in a.u.
      REAL(q) :: RMU       ! inverse real space cutoff in a.u.

      SAVE U1,U2,V1,A,BET,GAM
      DATA U1/1.0270741452992294_q/
      DATA U2/-0.230160617208092_q/
      DATA V1/0.6196884832404359_q/
      DATA A/3.2581_q/
      DATA BET/163.44_q/
      DATA GAM/4.7125_q/

      D=32._q/3._q/PI
! changed to a.u. instead of A, iG please check
      RHO = 3/(4._q*PI)/RS**3

!     G(0) from Burke, Perdew & Ernzerhof
      GRS32=(GAM+RS)**3._q/2._q
      GRS12=SQRT(GAM+RS)
      G0=D*(GRS32+BET)*EXP(-A*GRS12)

      C1N=U1*RS + U2*RS*RS
      C1D=1._q+V1*RS
      C1=C1N/C1D
      C2D=0.5_q*PI*RHO*RHO*(G0-0.5_q)
      C2=EC/C2D

      SCALE1=1._q + C1*RMU + C2*RMU*RMU

      DRS=-RS/3._q/RHO
      DG0=-0.5_q*D*EXP(-A*GRS12)*(A*BET-3._q*GAM-3._q*RS+A*GRS32)*DRS/GRS12

      DC1=(U1+2._q*U2*RS + U2*V1*RS*RS)*DRS/C1D/C1D
      DC2D=0.5_q*PI*RHO*(2._q*G0+RHO*DG0-1._q)

      DC2=(VC*C2D-EC*DC2D)/C2D/C2D

      SCALE2=DC1*RMU + DC2*RMU*RMU

      VC_SR=(VC*SCALE1-EC*SCALE2)/(SCALE1*SCALE1)
      RETURN
    END FUNCTION VC_SR

!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
            
!=====================================================================
! 
! Thomas-Fermi screened exchange functional
!  
!=====================================================================

    FUNCTION EX_SX(RS,RMU,IFLG)
      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
! Exchange energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION CX(2)
      SAVE CX
      DATA CX/0.9163305865663_q,1.1545041946774_q/
      EX_=-CX(IFLG)/RS

!     RMU is the Thomas Fermi wave vector supplied in a.u., RS in a.u. as well
      RHO = 3/(4._q*PI)/(RS)**3

! here is the "global" implementation Walter Wolf adopted in FLAPW
!      QF_GLOBAL=RMU**2*PI/4
!      A=RMU/QF_GLOBAL

      QF_LOCAL=(3._q*PI*PI*RHO)**(1._q/3._q)
      A=RMU/QF_LOCAL

      FRAC = 1 - (4._q/3._q) * A * ATAN(2._q/A) &
           & - (A*A/6._q) * ( 1._q - ( A*A/4._q  + 3._q) * LOG(1 + (4._q/(A*A) )))
      EX_SX=EX_*FRAC
      RETURN
    END FUNCTION EX_SX

    FUNCTION VX_SX(RS,RMU,IFLG)
      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
! Exchange energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      DIMENSION CX(2)
      SAVE CX
      DATA CX/0.9163305865663_q,1.1545041946774_q/

!     RMU is the Thomas Fermi wave vector supplied in a.u., RS in a.u. as well
      RHO = 3/(4._q*PI)/(RS)**3
      VX_=-CX(IFLG)/RS

! here is the "global" implementation Walter Wolf adopted in FLAPW
!      QF_GLOBAL=RMU**2*PI/4
!      A=RMU/QF_GLOBAL
!      FRAC = 1 - (4._q/3._q) * A * ATAN(2._q/A) &
!           & - (A*A/6._q) * ( 1._q - ( A*A/4._q  + 3._q) * LOG(1 + (4._q/(A*A) )))
!      VX_SX=VX_*FRAC

      QF_LOCAL=(3._q*PI*PI*RHO)**(1._q/3._q)
      A=RMU/QF_LOCAL

      FRAC = 1 + (1._q/2._q) * A*A - (1._q/8._q) * LOG( 1 + (4._q/(A*A) )) * A*A*A*A &
           &  - (1._q/2._q) * LOG( 1 + (4._q/(A*A))) *A*A
      VX_SX=EX_SX(RS,RMU,IFLG)+(1._q/3._q)*VX_*FRAC

      RETURN
    END FUNCTION VX_SX

! jP
!--------------------------------------------------------------------
! The routines for the evaluation of the spin-polarized 
! range-separated exchange-correlation energy in the LDA,
! LSDSR, ECORRLR, VCORRLR, EXCHANGELR, VEXCHANGELR, ECPW, and GPW,
! as well as the functions
! G0, G0D, DPOL, DPOLD, QRPA, QRPAD
! are downloaded from Paola Gori-Giorgi's webpage
! http://www.lct.jussieu.fr/pagesperso/gori/elegas.html
! [Phys. Rev. B 73, 155111 (2006)].
! They have been adapted to f90 by jP (04/25/08).
!--------------------------------------------------------------------
      SUBROUTINE lsdsr(rs,z,mu,excsr,vxcsrup,vxcsrdown)
!   Hartree atomic units used
!   for given density parameter 'rs', relative spin polarization 'z'= (nu -nd)/n
!   and cutoff parameter 'mu' 
!   gives the complementary  short-range exchange-correlation
!   energy  (i.e., xc energy of jellium minus xc energy of long-range
!   interacting electron gas) => 'excsr'
!   and the corresponding exchange-correlation potentials for
!   spin-up and spin-down electrons => 'vxcsrup','vxcsrdown'
!   from Paziani, Moroni, Gori-Giorgi, and Bachelet, cond-mat/0601343
        USE prec
        USE constant
        IMPLICIT NONE
        REAL(q) :: rs,z,mu,excsr,vxcsrup,vxcsrdown
        REAL(q) :: eclr,exlr,ec,ecd,ecz,ex
        REAL(q) :: vclrup,vclrdown,vxlrup,vxlrdown
        REAL(q) :: vxup,vxdown,vcup,vcdown
        REAL(q) :: alpha,cf

        alpha=(4.0_q/9.0_q/pi)**(1.0_q/3.0_q)
        cf=1.0_q/alpha
        
        ex=-3.0_q*cf/rs/8.0_q/pi*((1.0_q+z)**(4.0_q/3.0_q)+&
        &     (1.0_q-z)**(4.0_q/3.0_q))
        
        vxup=-(1.0_q+z)**(1.0_q/3.0_q)*(3.0_q/2.0_q/PI)**(2.0_q/3.0_q)/rs
        vxdown=-(1.0_q-z)**(1.0_q/3.0_q)*(3.0_q/2.0_q/PI)**(2.0_q/3.0_q)/rs
        
        call ecPW(rs,z,ec,ecd,ecz)
        vcup=ec-rs/3.0_q*ecd-(z-1.0_q)*ecz
        vcdown=ec-rs/3.0_q*ecd-(z+1.0_q)*ecz
        
        call exchangelr(rs,z,mu,exlr)
        call vexchangelr(rs,z,mu,vxlrup,vxlrdown)
        
        call ecorrlr(rs,z,mu,eclr)
        call vcorrlr(rs,z,mu,vclrup,vclrdown)
        
        excsr=ex+ec-(exlr+eclr)
        vxcsrup=vxup+vcup-(vxlrup+vclrup)
        vxcsrdown=vxdown+vcdown-(vxlrdown+vclrdown)

        RETURN
        END SUBROUTINE lsdsr

!----------------------------------------------------------------------
! Routine to calculate the spin-polarized
! short-range correlation energy 
! in the local density approximation -> see routine
! LSDSR; 
! Adaptation done by jP (04/25/08)
!----------------------------------------------------------------------
      SUBROUTINE ecorrsr(rs,z,mu,ecsr)
!   Hartree atomic units used
!   for given density parameter 'rs', relative spin polarization 'z'= (nu -nd)/n
!   and cutoff parameter 'mu' 
!   gives the complementary  short-range exchange-correlation
!   energy  (i.e., xc energy of jellium minus xc energy of long-range
!   interacting electron gas) => 'excsr'
!   and the corresponding exchange-correlation potentials for
!   spin-up and spin-down electrons => 'vxcsrup','vxcsrdown'
!   from Paziani, Moroni, Gori-Giorgi, and Bachelet, cond-mat/0601343
        USE prec
        USE constant
        IMPLICIT NONE
        REAL(q) :: rs,z,mu,ecsr
        REAL(q) :: eclr,ec,ecd,ecz

        
        call ecPW(rs,z,ec,ecd,ecz)
        call ecorrlr(rs,z,mu,eclr)
        
        ecsr=ec-eclr
        RETURN
      END SUBROUTINE ecorrsr

!----------------------------------------------------------------------
! Routine to calculate the spin-polarized
! short-range correlation potentials
! in the local density approximation -> see routine
! LSDSR; 
! Adaptations w.r.t. VASP done by jP (04/25/08)
!----------------------------------------------------------------------
      SUBROUTINE vcorrsr(rs,z,mu,vcsrup,vcsrdown)
!   Hartree atomic units used
!   for given density parameter 'rs', relative spin polarization 'z'= (nu -nd)/n
!   and cutoff parameter 'mu' 
!   gives the complementary  short-range exchange-correlation
!   energy  (i.e., xc energy of jellium minus xc energy of long-range
!   interacting electron gas) => 'excsr'
!   and the corresponding exchange-correlation potentials for
!   spin-up and spin-down electrons => 'vxcsrup','vxcsrdown'
!   from Paziani, Moroni, Gori-Giorgi, and Bachelet, cond-mat/0601343
        USE prec
        USE constant
        IMPLICIT NONE
        REAL(q) :: rs,z,mu,vcsrup,vcsrdown
        REAL(q) :: ec,ecd,ecz
        REAL(q) :: vclrup,vclrdown
        REAL(q) :: vcup,vcdown
        
        call ecPW(rs,z,ec,ecd,ecz)
        vcup=ec-rs/3.0_q*ecd-(z-1.0_q)*ecz
        vcdown=ec-rs/3.0_q*ecd-(z+1.0_q)*ecz
        
        call vcorrlr(rs,z,mu,vclrup,vclrdown)
        
        vcsrup=vcup-vclrup
        vcsrdown=vcdown-vclrdown
        RETURN
      END SUBROUTINE vcorrsr 
     
      SUBROUTINE ecorrlr(rs,z,mu,eclr)
!   Hartree atomic units used
!   for given density parameter rs, relative spin polarization z=(nu-nd)/n
!   and cutoff parameter mu 
!   gives the correlation energy of the LR gas
!    => eclr
        USE prec
        USE constant
        IMPLICIT NONE
        REAL(q) :: rs,z,mu,eclr,ec,ecd,ecz
        REAL(q) :: alpha,cf,phi
        REAL(q) :: d2anti,d3anti
        REAL(q) :: coe2,coe3,coe4,coe5
        REAL(q) :: a1,a2,a3,a4,b0
        REAL(q) :: q1a,q2a,q3a,t1a,t2a,t3a,adib

        alpha=(4.0_q/9.0_q/pi)**(1.0_q/3.0_q)
        cf=1.0_q/alpha
        
        phi=((1.0_q+z)**(2.0_q/3.0_q)+(1.0_q-z)**(2.0_q/3.0_q))/2.0_q
!  parameters from the fit
        adib   = 0.784949_q   
        q1a    = -0.388_q   
        q2a    = 0.676_q   
        q3a    = 0.547_q   
        t1a    = -4.95_q   
        t2a    = 1._q    
        t3a    = 0.31_q   
        
        b0=adib*rs
        
        d2anti=(q1a*rs+q2a*rs**2)*exp(-abs(q3a)*rs)/rs**2
        d3anti=(t1a*rs+t2a*rs**2)*exp(-abs(t3a)*rs)/rs**3
        
        coe2=-3._q/8._q/rs**3*(1._q-z**2)*(gg0(rs)-0.5_q)
        
        coe3=-(1._q-z**2)*gg0(rs)/(sqrt(2._q*pi)*rs**3)
        
        if(abs(z).eq.1._q) then
           
           coe4=-9._q/64._q/rs**3*(dpol(rs) &
                &        -cf**2*2**(5._q/3._q)/5._q/rs**2) 
           coe5=-9._q/40._q/(sqrt(2._q*pi)*rs**3)*dpol(rs)
           
        else
           
           coe4=-9._q/64._q/rs**3*(((1._q+z)/2._q)**2* &
                &        dpol(rs*(2/(1._q+z))**(1._q/3._q))+((1._q-z)/2._q)**2 &
                &        *dpol(rs*(2._q/(1._q-z))**(1._q/3._q))+ &
                &        (1.-z**2)*d2anti-cf**2/10._q*((1._q+z)**(8._q/3._q) &
                &        +(1.-z)**(8._q/3._q))/rs**2)
           
           coe5=-9._q/40._q/(sqrt(2._q*pi)*rs**3)*(((1._q+z)/2._q)**2 &
                &        *dpol(rs*(2._q/(1._q+z))**(1._q/3._q))+((1._q-z)/2._q)**2 &
                &        *dpol(rs*(2._q/(1._q-z))**(1._q/3._q))+(1._q-z**2)* &
                &        d3anti)
        endif
        
        call ecPW(rs,z,ec,ecd,ecz)
        
        a1=4._q*b0**6*coe3+b0**8*coe5
        a2=4._q*b0**6*coe2+b0**8*coe4+6._q*b0**4*ec
        a3=b0**8*coe3
        a4=b0**6*(b0**2*coe2+4._q*ec)
        
        eclr=(phi**3*Qrpa(mu*sqrt(rs)/phi)+a1*mu**3+a2*mu**4+a3*mu**5+ &
             &     a4*mu**6+b0**8*mu**8*ec)/((1._q+b0**2*mu**2)**4)
        
        RETURN
      END SUBROUTINE ecorrlr

      SUBROUTINE vcorrlr(rs,z,mu,vclrup,vclrdown)
!   Hartree atomic units used
!   for given density parameter rs, relative spin polarization z=(nu-nd)/n
!   and cutoff mu it gives the correlation LSD potential for LR interaction
!    => vclrup (spin-up electrons), vclrdown (spin-down electrons)
        USE prec
        USE constant
        IMPLICIT NONE
        REAL(q) :: rs,z,mu,eclr,eclrrs,eclrz,vclrup,vclrdown
        REAL(q) :: ec,ecd,ecz
        REAL(q) :: alpha,cf,phi
        REAL(q) :: d2anti,d3anti
        REAL(q) :: d2antid,d3antid,x
        REAL(q) :: coe2,coe3,coe4,coe5
        REAL(q) :: coe2rs,coe3rs,coe4rs,coe5rs
        REAL(q) :: coe2z,coe3z,coe4z,coe5z
        REAL(q) :: a1,a2,a3,a4,a5,b0,a1rs,a2rs,a3rs,a4rs,a5rs
        REAL(q) :: b0rs,a1z,a2z,a3z,a4z,a5z,b0z
        REAL(q) :: q1a,q2a,q3a,t1a,t2a,t3a,adib
        
        alpha=(4._q/9._q/pi)**(1._q/3._q)
        cf=1._q/alpha
        
        phi=((1._q+z)**(2._q/3._q)+(1._q-z)**(2._q/3._q))/2._q
!  parameters from the fit
        adib   = 0.784949_q   
        q1a    = -0.388_q   
        q2a    = 0.676_q   
        q3a    = 0.547_q   
        t1a    = -4.95_q   
        t2a    = 1._q    
        t3a    = 0.31_q   
        
        b0=adib*rs
        
        d2anti=(q1a+q2a*rs)*exp(-q3a*rs)/rs
        d3anti=(t1a+t2a*rs)*exp(-t3a*rs)/rs**2
        
        d2antid=-((q1a + q1a*q3a*rs + q2a*q3a*rs**2)/ &
             &    rs**2)*exp(-q3a*rs)
        d3antid=-((rs*t2a*(1 + rs*t3a) + t1a*(2 + rs*t3a))/ &
             &    rs**3)*exp(-rs*t3a)
        
        coe2=-3._q/8._q/rs**3*(1._q-z**2)*(gg0(rs)-0.5_q)
        coe2rs=-3._q/8._q/rs**3*(1._q-z**2)*gg0d(rs)+ &
             &     9._q/8._q/rs**4*(1._q-z**2)*(gg0(rs)-0.5_q)
        coe2z=-3._q/8._q/rs**3*(-2._q*z)*(gg0(rs)-0.5_q)

        coe3=-(1._q-z**2)*gg0(rs)/(sqrt(2._q*pi)*rs**3)
        coe3rs=-(1._q-z**2)*gg0d(rs)/(sqrt(2._q*pi)*rs**3)+ &
             &    3._q*(1._q-z**2)*gg0(rs)/(sqrt(2._q*pi)*rs**4) 
        coe3z=2._q*z*gg0(rs)/(sqrt(2._q*pi)*rs**3)
        
        if(abs(z).eq.1._q) then
           
           coe4=-9._q/64._q/rs**3*(dpol(rs) &
                &        -cf**2*2**(5._q/3._q)/5._q/rs**2)
           coe4rs=-3._q/rs*coe4-9._q/64._q/rs**3*(dpold(rs) &
                &        +2._q*cf**2*2**(5._q/3._q)/5._q/rs**3)
           coe4z=-9._q/64._q/rs**3*(dpol(rs)-rs/6._q*dpold(rs)-2._q*d2anti &
                &       -4._q/15._q/rs**2*cf**2*2._q**(5._q/3._q))*z
           coe5=-9._q/40._q/(sqrt(2._q*pi)*rs**3)*dpol(rs)
           coe5rs=-3._q/rs*coe5-9._q/40._q/(sqrt(2._q*pi)*rs**3)*dpold(rs)
           coe5z=-9._q/40._q/(sqrt(2._q*pi)*rs**3)*(dpol(rs)-rs/6._q* &
                &       dpold(rs)-2._q*d3anti)*z
           
        else
           
           coe4=-9._q/64._q/rs**3*(((1._q+z)/2._q)**2* &
                &        dpol(rs*(2/(1._q+z))**(1._q/3._q))+((1._q-z)/2._q)**2 &
                &        *dpol(rs*(2._q/(1._q-z))**(1._q/3._q))+ &
                &        (1.-z**2)*d2anti-cf**2/10._q*((1._q+z)**(8._q/3._q) &
                &        +(1.-z)**(8._q/3._q))/rs**2)
           coe4rs=-3._q/rs*coe4-9._q/64._q/rs**3*( &
                &        ((1._q+z)/2._q)**(5._q/3._q)*dpold(rs*(2/(1._q+z))** &
                &        (1._q/3._q))+((1._q-z)/2._q)**(5._q/3._q)* &
                &        dpold(rs*(2/(1._q-z))**(1._q/3._q))+(1._q-z**2)* &
                &        d2antid+cf**2/5._q*((1._q+z)**(8._q/3._q) &
                &        +(1._q-z)**(8._q/3._q))/rs**3)
           coe4z=-9._q/64._q/rs**3*(1._q/2._q*(1._q+z)* &
                &        dpol(rs*(2/(1._q+z))**(1._q/3._q))-1._q/2._q*(1._q-z)* &
                &        dpol(rs*(2/(1._q-z))**(1._q/3._q))-rs/6._q* &
                &        ((1._q+z)/2._q)**(2._q/3._q)*dpold(rs*(2/(1._q+z)) &
                &        **(1._q/3._q))+rs/6._q*((1._q-z)/2._q)**(2._q/3._q) &
                &        *dpold(rs*(2/(1._q-z))**(1._q/3._q))-2._q*z*d2anti- &
                &        4._q/15._q/rs**2*cf**2*((1._q+z)**(5._q/3._q)- &
                &        (1._q-z)**(5._q/3._q)))
           
           coe5=-9._q/40._q/(sqrt(2._q*pi)*rs**3)*(((1._q+z)/2._q)**2 &
                &        *dpol(rs*(2._q/(1._q+z))**(1._q/3._q))+((1._q-z)/2._q)**2 &
                &        *dpol(rs*(2._q/(1._q-z))**(1._q/3._q))+(1._q-z**2)* &
                &        d3anti)
           coe5rs=-3._q/rs*coe5-9._q/(40._q*sqrt(2._q*pi)*rs**3)*( &
                &        ((1._q+z)/2._q)**(5._q/3._q)*dpold(rs*(2/(1._q+z))** &
                &        (1._q/3._q))+((1._q-z)/2._q)**(5._q/3._q)* &
                &        dpold(rs*(2/(1._q-z))**(1._q/3._q))+(1._q-z**2)* &
                &        d3antid)
           coe5z=-9._q/40._q/(sqrt(2._q*pi)*rs**3)*(1._q/2._q*(1._q+z)* &
                &        dpol(rs*(2/(1._q+z))**(1._q/3._q))-1._q/2._q*(1._q-z)* &
                &        dpol(rs*(2/(1._q-z))**(1._q/3._q))-rs/6._q* &
                &        ((1._q+z)/2._q)**(2._q/3._q)*dpold(rs*(2/(1._q+z)) &
                &        **(1._q/3._q))+rs/6._q*((1._q-z)/2._q)**(2._q/3._q) &
                &        *dpold(rs*(2/(1._q-z))**(1._q/3._q))-2._q*z*d3anti)
           
        endif
        
        call ecPW(rs,z,ec,ecd,ecz)
        
        a1=4._q*b0**6*coe3+b0**8*coe5
        a1rs=24._q*adib*b0**5*coe3+4._q*b0**6*coe3rs+8._q*adib*b0**7* &
             &     coe5+b0**8*coe5rs
        a1z=4._q*b0**6*coe3z+b0**8*coe5z
        
        a2=4._q*b0**6*coe2+b0**8*coe4+6._q*b0**4*ec
        a2rs=24._q*adib*b0**5*coe2+4._q*b0**6*coe2rs+8._q*adib*b0**7* &
             &     coe4+b0**8*coe4rs+24._q*adib*b0**3*ec+6._q*b0**4*ecd
        a2z=4._q*b0**6*coe2z+b0**8*coe4z+6._q*b0**4*ecz
        
        a3=b0**8*coe3
        a3rs=8._q*adib*b0**7*coe3+b0**8*coe3rs
        a3z=b0**8*coe3z
        
        a4=b0**6*(b0**2*coe2+4._q*ec)
        a4rs=8._q*adib*b0**7*coe2+b0**8*coe2rs+24._q*adib*b0**5*ec+ &
             &     4._q*b0**6*ecd
        a4z=b0**6*(b0**2*coe2z+4._q*ecz)
        
        a5=b0**8*ec
        a5rs=8._q*adib*b0**7*ec+b0**8*ecd
        a5z=b0**8*ecz
        
        x=mu*sqrt(rs)/phi
        
        eclr=(phi**3*Qrpa(x)+a1*mu**3+a2*mu**4+a3*mu**5+ &
             &     a4*mu**6+a5*mu**8)/((1._q+b0**2*mu**2)**4)
        
        eclrrs=-4._q/(1._q+b0**2*mu**2)*2._q*adib*b0*mu**2*eclr+ &
             &     1._q/((1._q+b0**2*mu**2)**4)*(phi**2*mu/(2._q*sqrt(rs)) &
             &     *Qrpad(x)+ &
             &     a1rs*mu**3+a2rs*mu**4+a3rs*mu**5+a4rs*mu**6+a5rs*mu**8)
        
        if(z.eq.1._q) then
           vclrup=eclr-rs/3._q*eclrrs
           vclrdown=0._q
        elseif(z.eq.-1._q) then
           vclrup=0._q
           vclrdown=eclr-rs/3._q*eclrrs
        else
           
           eclrz=(phi**2*((1._q+z)**(-1._q/3._q)-(1._q-z)**(-1._q/3._q)) &
                &        *Qrpa(x)-phi*Qrpad(x)*mu*sqrt(rs)*((1._q+z)**(-1._q/3._q) &
                &        -(1._q-z)**(-1._q/3._q))/3._q+ &
                &        a1z*mu**3+a2z*mu**4+a3z*mu**5+ &
                &        a4z*mu**6+a5z*mu**8)/((1._q+b0**2*mu**2)**4)
           
           vclrup=eclr-rs/3._q*eclrrs-(z-1._q)*eclrz
           vclrdown=eclr-rs/3._q*eclrrs-(z+1._q)*eclrz
        endif
        RETURN
      END SUBROUTINE vcorrlr

      FUNCTION gg0(x)
        USE prec
        USE constant
!   on-top pair-distribution function
!   Gori-Giorgi and Perdew, PRB 64, 155102 (2001)
!   x -> rs
        IMPLICIT NONE
        REAL(q) :: gg0
        REAL(q) ::  C0f,D0f,E0f,F0f,x

        gg0              = 0.0_q
        C0f             = 0.0819306_q  
        D0f             = 0.752411_q    
        E0f             = -0.0127713_q   
        F0f             = 0.00185898_q   
        gg0=(1._q-(0.7317_q-D0f)*x+C0f*x**2+E0f*x**3+ &
             &     F0f*x**4)*exp(-abs(D0f)*x)/2._q
      RETURN
    END FUNCTION gg0

    FUNCTION gg0d(rs)
        USE prec
!   derivative of on-top pair-distribution function
!   Gori-Giorgi and Perdew, PRB 64, 155102 (2001)
        IMPLICIT NONE
        REAL(q) :: gg0d
        REAL(q) :: Bgg0,Cgg0,Dgg0,Egg0,Fgg0,rs

        gg0d             = 0.0_q
        Cgg0             = 0.0819306_q    
        Fgg0             = 0.752411_q     
        Dgg0             = -0.0127713_q   
        Egg0             = 0.00185898_q
        Bgg0             =0.7317_q-Fgg0
        gg0d=(-Bgg0+2*Cgg0*rs+3*Dgg0*rs**2+4*Egg0*rs**3)/2._q*exp(-Fgg0*rs) &
             &   - (Fgg0*(1 - Bgg0*rs + Cgg0*rs**2 + Dgg0*rs**3 + Egg0*rs**4))/ &
             &   2._q*exp(-Fgg0*rs)
        RETURN
      END FUNCTION gg0d

      FUNCTION dpol(rs)
        USE prec
        USE constant
        IMPLICIT NONE
        REAL(q) dpol
        REAL(q) :: cf,rs,p2p,p3p

        cf=(9._q*pi/4._q)**(1._q/3._q)
        p2p    = 0.04_q 
        p3p    = 0.4319_q   
        dpol=2._q**(5._q/3._q)/5._q*cf**2/rs**2*(1._q+(p3p-0.454555_q)*rs) &
             &     /(1._q+p3p*rs+p2p*rs**2)
      RETURN
    END FUNCTION dpol
    
    FUNCTION dpold(rs)
        USE prec
        USE constant
        IMPLICIT NONE
        REAL(q) dpold
        REAL(q) ::  cf,rs,p2p,p3p

        cf=(9._q*pi/4._q)**(1._q/3._q)
        p2p    = 0.04_q 
        p3p    = 0.4319_q   
        dpold=2._q**(5._q/3._q)/5._q*cf**2* &
             & (-2. + (0.454555 - 4.*p3p)*rs + &
             &    (-4.*p2p + &
             &       (0.90911 - 2.*p3p)*p3p)*rs**2 &
             &      + p2p*(1.363665 - 3.*p3p)* &
             &     rs**3)/ &
             &  (rs**3*(1. + p3p*rs + p2p*rs**2)**2)
      RETURN
    END FUNCTION dpold

    FUNCTION Qrpa(x)
      USE prec
      USE constant
      IMPLICIT NONE
        REAL(q) Qrpa
      REAL(q) :: a2,b2,c2,d2,x,Acoul

      Acoul=2._q*(log(2._q)-1._q)/pi**2
      a2              = 5.84605_q 
      c2              = 3.91744_q 
      d2              = 3.44851_q
      b2=d2-3._q/(2._q*pi*Acoul)*(4._q/(9._q*pi))**(1._q/3._q)
      Qrpa=Acoul*log((1._q+a2*x+b2*x**2+c2*x**3)/(1._q+a2*x+d2*x**2))
      RETURN
    END FUNCTION Qrpa

    FUNCTION Qrpad(x)
      USE prec
      USE constant
      IMPLICIT NONE
        REAL(q) Qrpad
      REAL(q) :: a2,b2,c2,d2,x,Acoul
      
      Acoul=2._q*(log(2._q)-1._q)/pi**2
      a2              = 5.84605_q 
      c2              = 3.91744_q 
      d2              = 3.44851_q
      b2=d2-3._q/(2._q*pi*Acoul)*(4._q/(9._q*pi))**(1._q/3._q)
      Qrpad=Acoul*((x*(b2*(2._q + a2*x) + &
           &      c2*x*(3._q + 2._q*a2*x) + &
           &      d2*(-2._q - a2*x + c2*x**3)))/ &
           &  ((1._q + a2*x + d2*x**2)* &
           &    (1._q + a2*x + b2*x**2 + c2*x**3)))
      RETURN
    END FUNCTION Qrpad
    
    SUBROUTINE exchangelr(rs,z,mu,exlr)
      USE prec
      USE constant
      USE c2f_interface, ONLY : ERRF
!   Hartree atomic units used
!   for given density parameter rs, relative spin polarization z=(nu-nd)/n
!   and cutoff mu it gives the exchange energy of the LR gas
!    => exlr
      IMPLICIT NONE
      REAL(q) :: rs,z,mu,exlr
      REAL(q) :: alpha,fx,y

      alpha=(4._q/9._q/pi)**(1._q/3._q)
      if(abs(z).eq.1._q) then
         y=mu*alpha*rs/2._q/2._q**(1._q/3._q)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + &
              &        sqrt(pi)*ERRF(1/(2.*y)))/pi)
         exlr=mu*fx
      else
         y=mu*alpha*rs/2._q/(1.+z)**(1._q/3._q)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + &
              &        sqrt(pi)*ERRF(1/(2.*y)))/pi)
         exlr=(1._q+z)*mu*fx/2._q
         y=mu*alpha*rs/2._q/(1.-z)**(1._q/3._q)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + &
              &        sqrt(pi)*ERRF(1/(2.*y)))/pi)
         exlr=exlr+(1._q-z)*mu*fx/2._q
      endif
      RETURN
    END SUBROUTINE exchangelr
    
    SUBROUTINE vexchangelr(rs,z,mu,vxlrup,vxlrdown)
      USE prec
      USE constant
      USE c2f_interface, ONLY : ERRF
!   Hartree atomic units used
!   for given density parameter rs, relative spin polarization z=(nu-nd)/n
!   and cutoff mu it gives the exchange LSD potential for LR interaction
!    => vxlrup (spin-up electrons), vxlrdown (spin-down electrons)
      IMPLICIT NONE
      REAL(q) :: rs,z,mu,vxlrup,vxlrdown
      REAL(q) :: alpha,fx,fx1,y,exlr,derrs,derz

      alpha=(4._q/9._q/pi)**(1._q/3._q)
      if(z.eq.1._q) then
         vxlrup=(rs*alpha*mu**2)/ &
              &   (2**(1._q/3._q)*pi) - (rs*alpha*mu**2)/(2**(1._q/3._q)*pi)* &
              &     exp(-2**(2._q/3._q)/(rs**2*alpha**2*mu**2)) - &
              &  (mu*ERRF(2**(1._q/3._q)/(rs*alpha*mu)))/sqrt(Pi)
         vxlrdown=0._q
      elseif(z.eq.-1._q) then
         vxlrdown=(rs*alpha*mu**2)/ &
              &   (2**(1._q/3._q)*pi) - (rs*alpha*mu**2)/(2**(1._q/3._q)*pi)* &
              &     exp(-2**(2._q/3._q)/(rs**2*alpha**2*mu**2)) - &
              &  (mu*ERRF(2**(1._q/3._q)/(rs*alpha*mu)))/sqrt(Pi)
         vxlrup=0._q
      else       
         y=mu*alpha*rs/2._q/(1.+z)**(1._q/3._q)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + &
              &        sqrt(pi)*ERRF(1/(2.*y)))/pi)
         fx1=(3._q*(1 + (-4._q + 4._q*exp(-1._q/(4._q*y**2)))*y**2))/pi
         derrs=1._q/4._q*(1._q+z)**(2._q/3._q)*mu**2*alpha*fx1
         derz=1._q/2._q*mu*fx-1._q/6._q*fx1*mu*y
         vxlrup=rs/3._q*derrs+(z-1._q)*derz
         vxlrdown=rs/3._q*derrs+(z+1._q)*derz
         
         y=mu*alpha*rs/2._q/(1.-z)**(1._q/3._q)
         fx=-((-3*y + 4*y**3 +(2*y - 4*y**3)*exp(-1./(4.*y**2)) + &
              &        sqrt(pi)*ERRF(1/(2.*y)))/pi)
         fx1=(3._q*(1 + (-4._q + 4._q*exp(-1._q/(4._q*y**2)))*y**2))/pi
         derrs=1._q/4._q*(1._q-z)**(2._q/3._q)*mu**2*alpha*fx1
         derz=-1._q/2._q*mu*fx+1._q/6._q*fx1*mu*y
         vxlrup=vxlrup+rs/3._q*derrs+(z-1._q)*derz
         vxlrdown=vxlrdown+rs/3._q*derrs+(z+1._q)*derz
         
         call exchangelr(rs,z,mu,exlr)
         vxlrup=exlr-vxlrup
         vxlrdown=exlr-vxlrdown
      endif
      RETURN
    END SUBROUTINE vexchangelr

!-----------------------------------------------------------------------
! correlation energy and its derivative w.r.t. rs and z at mu=infinity
! Perdew & Wang PRB 45, 13244 (1992)
!-----------------------------------------------------------------------
    SUBROUTINE ecPW(x,y,ec,ecd,ecz)
      USE prec
      USE constant
! in Hartree; ec=ec(rs,zeta)
! x -> rs; y -> zeta
!   ecd is d/drs ec
!   ecz is d/dz ec
      IMPLICIT NONE
      REAL(q) :: f02,ff,x,y,ec,ecd,ec0,ec0d,ec1,ec1d
      REAL(q) :: aaa,G,Gd,alfac,alfacd,ecz
      
      f02=4._q/(9._q*(2._q**(1._q/3._q)-1._q))
      
      ff=((1._q+y)**(4._q/3._q)+(1._q-y)**(4._q/3._q)- &
           &     2._q)/(2._q**(4._q/3._q)-2._q)
      
      aaa=(1._q-log(2._q))/pi**2
      call  GPW(x,aaa,0.21370_q,7.5957_q,3.5876_q,&
           &     1.6382_q,0.49294_q,G,Gd)
      ec0=G
      ec0d=Gd
      
      aaa=aaa/2._q
      call GPW(x,aaa,0.20548_q,14.1189_q,6.1977_q,&
           &     3.3662_q,0.62517_q,G,Gd)
      ec1=G
      ec1d=Gd
      call GPW(x,0.016887_q,0.11125_q,10.357_q,3.6231_q,&
           &     0.88026_q,0.49671_q,G,Gd)
      alfac=-G
      alfacd=-Gd
      
      ec=ec0+alfac*ff/f02*(1._q-y**4)+(ec1-ec0)*ff*y**4
      ecd=ec0d+alfacd*ff/f02*(1._q-y**4)+(ec1d-ec0d)* &
           &     ff*y**4
      ecz=alfac*(-4._q*y**3)*ff/f02+alfac*(1._q-y**4)/f02* &
           &     4._q/3._q*((1._q+y)**(1._q/3._q)-(1._q-y)**(1._q/3._q))/ &
           &     (2._q**(4._q/3._q)-2._q)+(ec1-ec0)*(4._q*y**3*ff+ &
           &     4._q/3._q*((1._q+y)**(1._q/3._q)-(1._q-y)**(1._q/3._q))/ &
           &     (2._q**(4._q/3._q)-2._q)*y**4)
      
      RETURN
    END SUBROUTINE ecPW
    
    SUBROUTINE GPW(x,Ac,alfa1,beta1,beta2,beta3,beta4,G,Gd)
      USE prec
!   Gd is d/drs G
      IMPLICIT NONE
      REAL(q) ::  G,Gd,Ac,alfa1,beta1,beta2,beta3,beta4,x
      G=-2._q*Ac*(1._q+alfa1*x)*dlog(1._q+1._q/(2._q* &
           &     Ac*(beta1*x**0.5_q+ &
           &     beta2*x+beta3*x**1.5_q+beta4*x**2)))
      Gd=(1._q+alfa1*x)*(beta2+beta1/(2._q*sqrt(x))+3._q*beta3* &
           &     sqrt(x)/2._q+2._q*beta4*x)/((beta1*sqrt(x)+beta2*x+ &
           &     beta3*x**(3._q/2._q)+beta4*x**2)**2*(1._q+1._q/ &
           &     (2._q*Ac*(beta1*sqrt(x)+beta2*x+beta3*x**(3._q/2._q)+ &
           &     beta4*x**2))))-2._q*Ac*alfa1*dlog(1._q+1._q/(2._q*Ac* &
           &     (beta1*sqrt(x)+beta2*x+beta3*x**(3._q/2._q)+ &
           &     beta4*x**2)))
      RETURN
    END SUBROUTINE GPW
!----------------------------------------------------------------------
!----------------------------------------------------------------------

END MODULE xclib
