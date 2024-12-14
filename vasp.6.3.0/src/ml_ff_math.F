#if defined(MPI) || defined(MPI_CHAIN)
#include "symbol.inc"

!****************************************************************************************************
! Module for atomic data
!****************************************************************************************************

      MODULE ATOM_DATA
        USE ML_FF_CONSTANT
        USE ML_FF_PREC
        IMPLICIT NONE

        INTEGER, PARAMETER :: MELEMENT=118

        CONTAINS

! Van-der-walls radius

        SUBROUTINE VDW_RADII (TYPE,R)
          IMPLICIT NONE
! Input variables
          CHARACTER(LEN=2), INTENT(IN) :: TYPE
! Output variables
          REAL(q), INTENT(OUT)         :: R
! Local variables
          INTEGER                      :: INTYP
          REAL(q)                      :: R0(1:MELEMENT)
          CHARACTER(LEN=2)             :: TYPE0(1:MELEMENT)   
! Data on Atomic radius
          DATA R0 /1.20                                                                                                                                                      ,1.40 &
                  ,2.20,1.90                                                                                                                        ,1.80,1.70,1.60,1.55,1.50,1.54 &
                  ,2.40,2.20                                                                                                                        ,2.10,2.10,1.95,1.80,1.80,1.88 &
                  ,2.80,2.20                                                                      ,2.30,2.15,2.05,2.05,2.05,2.05,2.00,2.00,2.00,2.10,2.10,2.10,2.05,1.90,1.90,2.02 &
                  ,2.90,2.55                                                                      ,2.40,2.30,2.15,2.10,2.05,2.05,2.00,2.05,2.10,2.20,2.20,2.25,2.20,2.10,2.10,2.15 &
                  ,3.00,2.70,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.50,2.25,2.20,2.10,2.05,2.00,2.00,2.05,2.10,2.05,2.20,2.30,2.30,2.30,2.30,2.30 &
                  ,2.40,2.40,2.40,2.40,2.40,2.30,2.40,2.40,2.40,2.40,2.40,2.40,2.40,2.40,2.40,2.40,2.40,2.40,2.40,2.40,2.40,2.40,2.40,2.40,2.40,2.40,2.40,2.40,2.40,2.40,2.40,2.40/
! Data on atomic symbols
          DATA TYPE0 /'H' ,'He'                                                                                                                                                       &
                     ,'Li','Be'                                                                                                                        ,'B' ,'C' ,'N' ,'O' ,'F' ,'Ne' &
                     ,'Na','Mg'                                                                                                                        ,'Al','Si','P' ,'S' ,'Cl','Ar' &
                     ,'K' ,'Ca'                                                                      ,'Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr' &
                     ,'Rb','Sr'                                                                      ,'Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I' ,'Xe' &
                     ,'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W' ,'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn' &
                     ,'Fr','Ra','Ac','Th','Pa','U' ,'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'/
          DO INTYP=1, MELEMENT
             IF(TYPE .EQ. TYPE0(INTYP)) THEN
                R=R0(INTYP)
                EXIT
             ENDIF
          ENDDO
        END SUBROUTINE VDW_RADII

      END MODULE ATOM_DATA

      MODULE ML_SPLINE

      CONTAINS
!**************** SUBROUTINE SPLCOF, SPLCOF_N0 *************************
! RCS:  $Id: ini.F,v 1.3 2002/08/14 13:59:39 kresse Exp $
!
!  Subroutine for calculating spline-coefficients
!  using the routines of the book 'numerical  recipes'
!  on input P(1,N) must contain x-values
!           P(2,N) must contain function-values
!  YP is the first derivatives at the first point
!  if >= 10^30 natural boundary-contitions (y''=0) are used
!
!  for point N always natural boundary-conditions are used in
!  SPLCOF, whereas SPLCOF_N0 assume 0 derivative at N
!  SPLCOF_NDER allows to specify a boundary condition
!  at both end points
!
!***********************************************************************
      SUBROUTINE SPLCOF_ML(P,N,NDIM,Y1P)
        USE ML_FF_PREC
        IMPLICIT REAL(q) (A-H,O-Z)
        REAL(q) :: P(:,:) !(NDIM,5)
        INTEGER :: N
        INTEGER :: NDIM
        REAL(q) :: Y1P
!
!     determination of spline coefficients
!     ------------------------------------
!     f = ((d*dx+c)*dx+b)*dx+a
!         between adjacent x - values
!
!     result
!     P-ARRAY
!     P(I,1) = X(I)
!     P(I,2) = A(I) = F(I)
!     P(I,3) = B(I)
!     P(I,4) = C(I)
!     P(I,5) = D(I)
!
        IF (Y1P> .99E30_q) THEN
           P(1,4)=0.0_q
           P(1,3)=0.0_q
        ELSE
           P(1,4)=-.5_q
           P(1,3)=(3._q/(P(2,1)-P(1,1)))*((P(2,2)-P(1,2))/ &
                  (P(2,1)-P(1,1))-Y1P)
        ENDIF
        DO I=2,N-1
           S=(P(I,1)-P(I-1,1))/(P(I+1,1)-P(I-1,1))
           R=S*P(I-1,4)+2._q
           P(I,4)=(S-1._q)/R
           P(I,3)=(6*((P(I+1,2)-P(I,2))/(P(I+1,1)-P(I,1))- &
                  (P(I,2)-P(I-1,2))/(P(I,1)-P(I-1,1)))/ &
                  (P(I+1,1)-P(I-1,1))-S*P(I-1,3))/R
        ENDDO
        P(N,4)=0.0_q
        P(N,3)=0.0_q
        DO I=N-1,1,-1
           P(I,4)=P(I,4)*P(I+1,4)+P(I,3)
        ENDDO
        DO I=1,N-1
           S= P(I+1,1)-P(I,1)
           R=(P(I+1,4)-P(I,4))/6
           P(I,5)=R/S
           P(I,4)=P(I,4)/2.0_q
           P(I,3)=(P(I+1,2)-P(I,2))/S-(P(I,4)+R)*S
        ENDDO
        RETURN
      END SUBROUTINE SPLCOF_ML
!
!  helper routine, which evaluates the spline fit at a specific
!  position
!
      SUBROUTINE SPLVAL_ML(X,F,FDER,P,NAC,NDIM)
        USE ML_FF_PREC
        IMPLICIT REAL(q) (A-H,O-Z)
        REAL(q) :: P(:,:) !(NDIM,5)
!  interval bisectioning
        I=1
        J=NAC
        IF (X   <P(I,1)) GO TO 60
        IF (X   <P(J,1)) GO TO 70
        K=J-1
        GOTO 90
60      K=1
        GOTO 90
70      K=(I+J)/2
        IF(I==K) GOTO 90
        IF (X   <P(K,1)) GO TO 80
        I=K
        GOTO 70
80      J=K
        GOTO 70
!
90      DX=X   -P(K,1)
        F   =((P(K,5)*DX+P(K,4))*DX+P(K,3))*DX+P(K,2)
        FDER=(3.0_q*P(K,5)*DX+2.0_q*P(K,4))*DX+P(K,3)
      END SUBROUTINE SPLVAL_ML
  
      END MODULE ML_SPLINE

!****************************************************************************************************
! Basic mathematical functions necessary for calculating SOAP kernel.
!****************************************************************************************************

      MODULE BASIS_SOAP
        USE ML_FF_CONSTANT
        USE ML_FF_PREC
        IMPLICIT NONE

        CONTAINS

! This subroutine construct radial spherical Bessel function basis sets.

        SUBROUTINE RJ_RSBESSEL (DR,FWIN,IWINDOW,LMAX,LVARTRAN,LWINDOW,MRB,NR,NRB,NVARTRAN,QGRID,QMAX,RB,RCUT,RG,RMEM,RNORM,WR)
          IMPLICIT NONE
! Input variables
          REAL(q), INTENT(IN) :: DR
          INTEGER, INTENT(IN) :: IWINDOW
          INTEGER, INTENT(IN) :: LMAX
          LOGICAL, INTENT(IN) :: LVARTRAN
          LOGICAL, INTENT(IN) :: LWINDOW
          INTEGER, INTENT(IN) :: MRB
          INTEGER, INTENT(IN) :: NR
          INTEGER, INTENT(IN) :: NVARTRAN
          REAL(q), INTENT(IN) :: RCUT
! Output variables
          REAL(q), INTENT(OUT) :: FWIN(0:,:) !(0:LMAX,1:MRB)
          INTEGER, INTENT(OUT) :: NRB(0:) !(0:LMAX)
          REAL(q), INTENT(OUT) :: QGRID(:,0:) !(1:MRB,0:LMAX)
          REAL(q), INTENT(OUT) :: QMAX
          REAL(q), INTENT(OUT) :: RB(:,:,0:) !(1:NR,1:MRB,0:LMAX)
          REAL(q), INTENT(OUT) :: RG(:) !(1:NR)
          REAL(q), INTENT(OUT) :: RMEM
          REAL(q), INTENT(OUT) :: RNORM(0:,:) !(0:LMAX,1:MRB)
          REAL(q), INTENT(OUT) :: WR
! Local variables
          INTEGER :: IR
          INTEGER :: IRB
          INTEGER :: L
          REAL(q) :: X
! Function
!          REAL(q) :: WINDOW_FUNCTION
! Set the temporary memory usage as zero
          RMEM=0.0_q
! Determine the Q-grids
          QMAX=1.0_q/(1E-12_q)
          DO L=0, LMAX
             CALL OPT_BEZERO(QGRID(1:,L),L,MRB)
             QMAX=MIN(QMAX,QGRID(MRB,L))
          ENDDO
          DO L=0, LMAX
             NRB(L)=0
             DO IRB=1, MRB
                IF(QGRID(IRB,L) .LE. QMAX) THEN
                   NRB(L)=NRB(L)+1
                ELSE
                   QGRID(IRB,L)=0.0_q
                ENDIF
             ENDDO
          ENDDO
! If LVARTRAN=.TRUE., the bessel function is spanned in the transformed radial coordinate.
          IF(LVARTRAN) THEN
             DO L=0, LMAX
                DO IRB=1, NRB(L)
                   QGRID(IRB,L)=QGRID(IRB,L)
                ENDDO
             ENDDO
          ELSE
! Otherwise, the function is spanned in the standard radial coordinate.
             QMAX=QMAX/RCUT
             DO L=0, LMAX
                DO IRB=1, NRB(L)
                   QGRID(IRB,L)=QGRID(IRB,L)/RCUT
                ENDDO
             ENDDO
          ENDIF
! Calculate window function
          IF(LWINDOW) THEN
             DO L=0, LMAX
                DO IRB=1, NRB(L)
                   FWIN(L,IRB)=WINDOW_FUNCTION (IWINDOW,QGRID(IRB,L),QMAX)
                ENDDO
             ENDDO
          ENDIF
! Calculate the radial resolution
          WR=2.0_q*PI/QMAX
! Determine the real space radial grids.
! If LVARTRAN=.TRUE., the radial coordinate is transformed by using the polynomial function.
          DO IR=1, NR
             RG(IR)=REAL(IR,q)*DR
          ENDDO
! Calculate the radial basis sets.
! If LVARTRAN=.TRUE., the bessel function is calculated on transformed radial grid.
          DO L=0, LMAX
             DO IRB=1, NRB(L)
                DO IR=1, NR
                   X=QGRID(IRB,L)*RG(IR)
                   CALL SBESSEL(X,RB(IR,IRB,L),L)
                ENDDO
             ENDDO
          ENDDO
! Normalization.
! If LVARTRAN=.TRUE., the normalization is done with respect to the integration over tne transformed radial coordinate.
          DO L=0, LMAX
             DO IRB=1, NRB(L)
                RNORM(L,IRB)=0.0_q
                DO IR=1, NR
                   RNORM(L,IRB)=RNORM(L,IRB)+(RB(IR,IRB,L)*RG(IR))**2*DR
                ENDDO
                RNORM(L,IRB)=SQRT(RNORM(L,IRB))
                RB(1:NR,IRB,L)=RB(1:NR,IRB,L)/RNORM(L,IRB)
             ENDDO
          ENDDO
        END SUBROUTINE RJ_RSBESSEL

!*******************************************************************
!  SUBROUTINE BEZERO
!  searches for NQ zeros j(qr)
!  i/o:
!         XNULL(NQ) result
!         L           quantum number l
!  great full spaghetti code (written by gK)
!********************************************************************

        SUBROUTINE OPT_BEZERO(XNULL,L,NQ)
          IMPLICIT NONE
          REAL(q) :: BJ1,BJ2
          REAL(q), PARAMETER :: BREAK= 1E-10_q
          REAL(q) :: DUMMY
          REAL(q) :: ETA
          INTEGER, INTENT(IN) :: L
          INTEGER :: N
          INTEGER, INTENT(IN) :: NQ
          REAL(q), PARAMETER :: STEP=.1_q
          REAL(q) :: SSTEP
          REAL(q) :: X,XX
          REAL(q), INTENT(OUT) :: XNULL(:) !(NQ)
! initialization
          X=STEP
          N=0
! entry point for next q_n
30        CALL SBESSE2(X, BJ1, DUMMY,  L)
! coarse search
10        X=X+STEP
          CALL SBESSE2(X, BJ2, DUMMY,  L)
! found one point
          IF(BJ1*BJ2 < 0) THEN
             ETA=0.0_q
! interval bisectioning
             SSTEP=STEP
             XX   =X
20           SSTEP=SSTEP/2
             IF (BJ1*BJ2 < 0) THEN
                XX=XX-SSTEP
             ELSE
                XX=XX+SSTEP
             ENDIF
             CALL SBESSE2(XX, BJ2, DUMMY,  L)
             IF (SSTEP > BREAK) GOTO 20
             N=N+1
             XNULL(N)=XX
             IF (N == NQ) RETURN
             GOTO 30
          ENDIF
          GOTO 10
        END SUBROUTINE

!*******************************************************************
!>
!>     caculates the spherical besselfunction and its derivate
!>     The computation is done depending on the value of the argument x
!>
!>      x>DEL: recursion j_0(x) = sin(x)
!>                       j_1(x) = cos(x)/x - sin(x)/x/x
!>                       [...]
!>                       j_n(x) = -j_{n-2}(x) + (2n+1)/x j_{n-1}(x)
!>
!>      x<DEL: taylor expansion
!>
!********************************************************************

        SUBROUTINE SBESSEL(X, BJ, LMAX)
          IMPLICIT NONE
          REAL(q),INTENT(IN) :: X !< argument x of the Bessel function j_n(x)
          REAL(q), INTENT(OUT) :: BJ !< value of the Bessel function j_n(x)
          INTEGER, INTENT(IN) :: LMAX !< index n of the Bessel function j_n(x)
          ! local variables
          REAL(q), PARAMETER :: DEL=1._q
          REAL(q) :: BJP
          INTEGER :: L
          REAL(q) :: TEMP
          IF ( X< DEL) THEN
             BJ=SBESSITER(X,LMAX)
          ELSE
             BJ = DSIN(X)/X
             BJP=(-DCOS(X)+BJ)/X
             DO L=2,LMAX+1
                TEMP= -BJ + (2*L-1)*BJP/ X
                BJ  = BJP
                BJP = TEMP
             ENDDO
          ENDIF
          RETURN
        END

!*******************************************************************
!>
!> Compute the bessel function and its derivative
!> The first derivative of a bessel function can be written as
!> j'_n(x) = - j_{n+1}(x) + n/x j_n(x)
!>
!*******************************************************************
       SUBROUTINE SBESSE2(X, BJ, BJP , LMAX)
         IMPLICIT NONE
         REAL(q), INTENT(IN) :: X !< argument x of the Bessel function j_n(x)
         REAL(q), INTENT(OUT) :: BJ !< value of the Bessel function j_n(x)
         REAL(q), INTENT(OUT) :: BJP !< value of the first derivative of the Bessel function j'_n(x)
         INTEGER, INTENT(IN) :: LMAX !< index n of the Bessel function j_n(x)
         ! local variables
         REAL(q), PARAMETER :: DEL=1._q
         INTEGER :: L
         REAL(q) :: TEMP
         IF ( X< DEL) THEN
            BJ  = SBESSITER(X,LMAX)
            BJP = SBESSITER(X,LMAX+1)
            BJP = -BJP+ BJ*LMAX/X
         ELSE
            BJ = DSIN(X)/X
            BJP=(-DCOS(X)+BJ)/X
            DO L=2,LMAX+1
               TEMP= -BJ + (2*L-1)*BJP/ X
               BJ  = BJP
               BJP = TEMP
            ENDDO
            BJP = -BJP+ BJ*LMAX/X
         ENDIF
         RETURN
       END

       FUNCTION SBESSITER(X,LMAX)
         USE PRINT2TERMINAL
         IMPLICIT NONE
         REAL(q)        :: X
         INTEGER        :: LMAX
         REAL(q)        :: SBESSITER
         INTEGER        :: L
         INTEGER        :: I
         REAL(q)        :: X2
         REAL(q)        :: FAK
         REAL(q)        :: SUM
         REAL(q)        :: D
         REAL(q),PARAMETER :: ABBRUCH = 1E-16_q
         INTEGER,PARAMETER :: ITERMAX = 40 

         X2= X*X/2
         FAK = 1._q
         DO L=1,LMAX
           FAK = FAK*X /(2*L+1)
         ENDDO

         SUM = 1._q
         D   = 1._q
         DO I=1,ITERMAX
           D = -D* X2/(I*(2*(LMAX+I)+1))
           IF ( ABS(D)< ABBRUCH) GOTO 230
           SUM = SUM + D
         ENDDO

         CALL ERROR_PRINT2TERMINAL(" &
              ERROR: SBESSELITER is not converging") 

 230     SBESSITER= FAK * SUM
         RETURN
       END

! Window function

       FUNCTION WINDOW_FUNCTION (IWINDOW,QGRID,QMAX)
         IMPLICIT NONE
         REAL(q), PARAMETER :: A11=5.0E-1_q                                                    ! Parameter for Hann window.
         REAL(q), PARAMETER :: A12=5.4E-1_q, A22=4.6E-1_q                                      ! Parameters for Hamming window.
         REAL(q), PARAMETER :: A13=4.2E-1_q, A23=5.0E-1_q, A33=0.8E-1_q                        ! Parameters for Blackman window.
         REAL(q), PARAMETER :: A14=0.355768_q ,A24=0.487396_q ,A34=0.144232_q ,A44=0.012604_q  ! Parameters for Nuttall window.
         REAL(q), PARAMETER :: A15=0.3635819_q,A25=0.4891775_q,A35=0.1365995_q,A45=0.0106411_q ! Parameters for Blackman-Nuttwall window.
         REAL(q), PARAMETER :: A16=0.35875_q  ,A26=0.48829_q  ,A36=0.14128_q  ,A46=0.01168_q   ! Parameters for Blackman-Harris window.
         INTEGER :: IWINDOW
         REAL(q), PARAMETER :: PI=3.14159265358979_q
         REAL(q) :: QGRID
         REAL(q) :: QMAX,QMAX_LOC
         REAL(q) :: WINDOW_FUNCTION
         QMAX_LOC=QMAX
! If IWINDOW=1, calculate Hann function.
         IF(IWINDOW .EQ. 1) THEN
            IF(QGRID .LE. QMAX_LOC) THEN
               WINDOW_FUNCTION=A11*(DCOS(PI*QGRID/QMAX_LOC)+1.0_q)
            ELSE
               WINDOW_FUNCTION=0.0_q
            ENDIF
! else if IWINDOW=2, calculate Hamming function.
         ELSE IF(IWINDOW .EQ. 2) THEN
            IF(QGRID .LE. QMAX_LOC) THEN
               WINDOW_FUNCTION=A12+A22*DCOS(PI*QGRID/QMAX_LOC)
            ELSE
               WINDOW_FUNCTION=A12-A22
            ENDIF
! else if IWINDOW=3, calculate Blackman function.
         ELSE IF(IWINDOW .EQ. 3) THEN
            IF(QGRID .LE. QMAX_LOC) THEN
               WINDOW_FUNCTION=A13+A23*DCOS(PI*QGRID/QMAX_LOC)+A33*DCOS(2.0_q*PI*QGRID/QMAX_LOC)
            ELSE
               WINDOW_FUNCTION=A13-A23+A33
            ENDIF
! else if IWINDOW=4, calculate Nuttall function.
         ELSE IF(IWINDOW .EQ. 4) THEN
            IF(QGRID .LE. QMAX_LOC) THEN
               WINDOW_FUNCTION=A14+A24*DCOS(PI*QGRID/QMAX_LOC)+A34*DCOS(2.0_q*PI*QGRID/QMAX_LOC)+A44*DCOS(3.0_q*PI*QGRID/QMAX_LOC)
            ELSE
               WINDOW_FUNCTION=A14-A24+A34-A44
            ENDIF
! else if IWINDOW=5, calculate Blackman-Nuttall function.
         ELSE IF(IWINDOW .EQ. 5) THEN
            IF(QGRID .LE. QMAX_LOC) THEN
               WINDOW_FUNCTION=A15+A25*DCOS(PI*QGRID/QMAX_LOC)+A35*DCOS(2.0_q*PI*QGRID/QMAX_LOC)+A45*DCOS(3.0_q*PI*QGRID/QMAX_LOC)
            ELSE
               WINDOW_FUNCTION=A15-A25+A35-A45
            ENDIF
! else if IWINDOW=6, calculate Blackman-Harris function.
         ELSE IF(IWINDOW .EQ. 6) THEN
            IF(QGRID .LE. QMAX_LOC) THEN
               WINDOW_FUNCTION=A16+A26*DCOS(PI*QGRID/QMAX_LOC)+A36*DCOS(2.0_q*PI*QGRID/QMAX_LOC)+A46*DCOS(3.0_q*PI*QGRID/QMAX_LOC)
            ELSE
               WINDOW_FUNCTION=A16-A26+A36-A46
            ENDIF
         ENDIF
         RETURN
       END

! Subroutine calculate the errors caused by the expansion of 
! atom distribution by the radial basis sets

        SUBROUTINE RAD_ERROR (DR,ERR,ICUT,LMAX,LMETRIC,LVARTRAN,MRB,NMETRIC,NR,NRB,NVARTRAN,RB,RCUT,RMETRIC,WION)
          USE ML_FF_CONSTANT
          USE ML_FF_PREC
          IMPLICIT NONE
! Input variables
          REAL(q), INTENT(IN)    :: DR
          INTEGER, INTENT(IN)    :: ICUT
          INTEGER, INTENT(IN)    :: LMAX
          LOGICAL, INTENT(IN)    :: LMETRIC
          LOGICAL, INTENT(IN)    :: LVARTRAN
          INTEGER, INTENT(IN)    :: MRB
          INTEGER, INTENT(IN)    :: NMETRIC
          INTEGER, INTENT(IN)    :: NR
          INTEGER, INTENT(IN)    :: NRB(0:) !(0:LMAX)
          INTEGER, INTENT(IN)    :: NVARTRAN
          REAL(q), INTENT(IN)    :: RB(:,:,0:) !(1:NR,1:MRB,0:LMAX)
          REAL(q), INTENT(IN)    :: RCUT
          REAL(q), INTENT(IN)    :: RMETRIC
          REAL(q), INTENT(IN)    :: WION
! In and output variables
          REAL(q), INTENT(INOUT) :: ERR
! Local variables
          REAL(q)                :: CL(1:MRB,0:LMAX)
          REAL(q)                :: F0(1:NR,0:LMAX)
          REAL(q)                :: F1(1:NR,0:LMAX)
          REAL(q)                :: FCUT(1:NR)
          REAL(q)                :: FMET(1:NR)
          INTEGER                :: IR
          INTEGER                :: IRB
          INTEGER                :: JR
          INTEGER                :: L
          REAL(q)                :: RG(1:NR) ! Grids in the original radial space. Please, note that its definition is different from that of FFM%RG.
          REAL(q)                :: RION(1:NR) ! Ionic position.
          REAL(q)                :: RIL(1:NR,0:LMAX)
          REAL(q)                :: RLIMIT
          REAL(q)                :: RR(1:NR)
          REAL(q)                :: RRI
          REAL(q)                :: RRJ(1:NR)
          REAL(q)                :: TG(1:NR) ! Grids in the transformed radial space.
! Constant.
          RLIMIT=0.5_q/AUTOA
! Determine the (transformed) radial mesh.
          IF(LVARTRAN) THEN
             DO IR=1, NR
                TG(IR)=REAL(IR,q)*DR
                RG(IR)=RCUT*TG(IR)**NVARTRAN
             ENDDO
          ELSE
             DO IR=1, NR
                TG(IR)=REAL(IR,q)*DR
                RG(IR)=TG(IR)
             ENDDO
          ENDIF
! Calculate cutoff function.
          IF(ICUT.EQ.1) THEN
             CALL BP0 (FCUT,NR,RCUT,RG)
          ELSE IF(ICUT.EQ.2) THEN
             CALL MO0 (FCUT,NR,RCUT,RG)
          ELSE IF(ICUT.EQ.3) THEN
             CALL RJ0 (FCUT,NR,RCUT,RG)
          ELSE IF(ICUT.EQ.4) THEN
             CALL WMC0 (FCUT,NR,RCUT,RG)
          ENDIF
! Calculate radial metric function.
          IF(LMETRIC) THEN
             CALL PMET0 (FMET,NMETRIC,NR,RG,RCUT,RMETRIC)
          ELSE
             FMET(1:NR)=1.0_q
          ENDIF
! Scan the radial direction.
! IR corresponds to the ionic position that is moved in the original radial space.
          DO IR=1, NR
!          IR=12
             IF(RG(IR) .GT. RLIMIT) THEN
                RRI=WION*RG(IR)**2
! Calculate modified spherical Bessel function within the original radial space.
! Because the modified spherical Bessel function easily diverse at large RG,
! the exponential functions that will be multiplied later are introduced 
! in the subroutine.
! Hence, the calculated modified spherical Bessel functions are not the
! original ones but the ones multiplied with exponential functions.
                DO JR=1, NR
                   RR(JR)=2.0_q*WION*RG(JR)*RG(IR)
                   RRJ(JR)=WION*RG(JR)**2
                ENDDO
                CALL IL0 (LMAX,NR,NR,RR,RRI,RRJ,RIL)
! Calculate the original functions.
                DO L=0, LMAX
                   DO JR=1, NR
! Calculate original radial functions
                      F0(JR,L)=4.0_q*PI*(WION/PI)**1.5*FCUT(IR)*FMET(JR)*RIL(JR,L)
                   ENDDO
                ENDDO
! Calculate the expanded functions.
                DO L=0, LMAX
! Calculate the expansion coefficients.
! If LMETRIC=.TRUE., this integration is done over the transformed radial coordinate.
                   DO IRB=1, NRB(L)
                      CL(IRB,L)=0.0_q
                      DO JR=1, NR
                         CL(IRB,L)=CL(IRB,L)+4.0_q*PI*(WION/PI)**1.5*FCUT(IR)*FMET(JR)*RB(JR,IRB,L)*RIL(JR,L)*TG(JR)**2*DR
                      ENDDO
                   ENDDO
! Reconstruct the radial function
                   F1(1:NR,L)=0.0_q
                   DO IRB=1, NRB(L)
                      DO JR=1, NR
                         F1(JR,L)=F1(JR,L)+CL(IRB,L)*RB(JR,IRB,L)
                      ENDDO
                   ENDDO
                ENDDO
! Evaluate the maximum error
                DO L=0, LMAX
                   DO JR=1, NR
                      ERR=MAX(ERR,ABS(F1(JR,L)-F0(JR,L)))
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
        END SUBROUTINE RAD_ERROR

! Subroutine calculating radial function for angular descriptor.
! The data on the radial meshes will be calculated.

        SUBROUTINE RAD_FUNC (DR,DRSPL,FNL,IBROAD,ICUT,LMAX,LMETRIC,LVARTRAN,MRB,MSPL, &
                   NMETRIC,NR,NRB,NSPL,NVARTRAN,QGRID,RB,RC,RCUT,RG,RMEM,RMETRIC,RNORM,USQ,WION,WR)
          USE ML_FF_STRUCT, ONLY : ML_MPI_PAR
          USE ML_SPLINE
          IMPLICIT NONE
! Input variables
          REAL(q), INTENT(IN) :: DR
          INTEGER, INTENT(IN) :: IBROAD
          INTEGER, INTENT(IN) :: ICUT
          INTEGER, INTENT(IN) :: LMAX
          LOGICAL, INTENT(IN) :: LMETRIC
          LOGICAL, INTENT(IN) :: LVARTRAN
          INTEGER, INTENT(IN) :: MRB
          INTEGER, INTENT(IN) :: MSPL
          INTEGER, INTENT(IN) :: NMETRIC
          INTEGER, INTENT(IN) :: NR
          INTEGER, INTENT(IN) :: NRB(0:) !(0:LMAX)
          INTEGER, INTENT(IN) :: NSPL
          INTEGER, INTENT(IN) :: NVARTRAN
          REAL(q), INTENT(IN) :: QGRID(:,0:) !(1:MRB,0:LMAX)
          REAL(q), INTENT(IN) :: RB(:,:,0:) !(1:NR,1:MRB,0:LMAX)
          REAL(q), INTENT(IN) :: RC(:) !(1:MRB)
          REAL(q), INTENT(IN) :: RCUT
          REAL(q), INTENT(IN) :: RG(:) !(1:NR)
          REAL(q), INTENT(IN) :: RMETRIC
          REAL(q), INTENT(IN) :: RNORM(0:,:) !(0:LMAX,1:MRB)
          REAL(q), INTENT(IN) :: USQ(:,:) !(1:MRB,1:MRB)
          REAL(q), INTENT(IN) :: WION
          REAL(q), INTENT(IN) :: WR
! Output variables
          REAL(q), INTENT(OUT) :: DRSPL
          REAL(q), INTENT(OUT) :: FNL(:,:,0:,:) !(1:MSPL,1:5,0:LMAX,1:MRB)
          REAL(q), INTENT(OUT) :: RMEM
! Local variables
          REAL(q) CR(1:NR,0:LMAX)
          REAL(q) DCR(1:NR,0:LMAX)
          REAL(q) DDOT
          REAL(q) DFCUT(1)
          REAL(q) DFMET(1)
          REAL(q) DFNL_LOC(1:MRB,0:LMAX)
          REAL(q) DEXP22(1)
          REAL(q) DRIL(1:NR,0:LMAX)
          REAL(q) DUMMY
          REAL(q) EXP21(1:NR)
          REAL(q) EXP22(1:MSPL)
          REAL(q) FCUT(1:MSPL)
          REAL(q) FMET(1:MAX(MSPL,NR))
          INTEGER I
          INTEGER IR
          INTEGER ISPL
          INTEGER L
          INTEGER N,NN
          REAL(q) PI4
          REAL(q) PREFAC
          REAL(q) RG_LOC(1:NR)
          REAL(q) RIL(1:NR,0:LMAX)
          REAL(q) RR(1:MAX(NR,MSPL))
          REAL(q) RRI
          REAL(q) RRJ(1:MAX(NR,MSPL))
          REAL(q) TG_LOC(1:NR)
          REAL(q) X
! Calculate the temporary memory usage
          RMEM=(REAL(MSPL,q)*8*2+ & ! FCUT and DFCUT
                REAL(MSPL,q)*8*2+ & ! FMET and DFMET
                REAL(NR,q)*REAL((LMAX+1),q)*8+ & ! CR
                REAL(NR,q)*REAL((LMAX+1),q)*8+ & ! DCR
                REAL(MRB,q)*REAL((LMAX+1),q)*8+ & ! DFNL_LOC
                REAL(NR,q)*REAL((LMAX+1),q)*8+ & ! DRIL
                REAL(NR,q)*8+ & ! EXP21
                REAL(MSPL,q)*8+ & !EXP22
                REAL(NR,q)*REAL((LMAX+1),q)*8+ & ! RIL
                REAL(MAX(NR,MSPL),q)*8 & ! RR
               )/1E6_q
! Constant
          PI4=4.0_q*PI
          PREFAC=PI4*(WION/PI)**1.5
! Set radial mesh for spline interpolation.
          DRSPL=RCUT/REAL(NSPL,q)
          DO L=0, LMAX
             DO N=1, NRB(L)
                DO ISPL=1, NSPL
                   FNL(ISPL,1,L,N)=ISPL*DRSPL
                ENDDO
             ENDDO
          ENDDO
! Calculate cutoff function on the spline mesh.
          IF(ICUT.EQ.1) THEN
             CALL BP0 (FCUT,NSPL,RCUT,FNL(1:,1,0,1))
             CALL BP1 (FCUT(1:),DFCUT(1:),1,RCUT,FNL(:1,1,0,1))
          ELSE IF(ICUT.EQ.2) THEN
             CALL MO0 (FCUT,NSPL,RCUT,FNL(1:,1,0,1))
             CALL MO1 (FCUT(1:),DFCUT(1:),1,RCUT,FNL(1:,1,0,1))
          ELSE IF(ICUT.EQ.3) THEN
             CALL RJ0 (FCUT,NSPL,RCUT,FNL(1:,1,0,1))
             CALL RJ1 (FCUT(1:),DFCUT(1:),1,RCUT,FNL(1:,1,0,1))
          ELSE IF(ICUT.EQ.4) THEN
             CALL WMC0 (FCUT,NSPL,RCUT,FNL(1:,1,0,1))
             CALL WMC1 (FCUT(1:),DFCUT(1:),1,RCUT,FNL(1:,1,0,1))
          ENDIF
! Calculate original radial grid and transformed one for numerical integration.
          IF(LVARTRAN) THEN
             DO IR=1, NR
                TG_LOC(IR)=REAL(IR,q)*DR
                RG_LOC(IR)=RCUT*TG_LOC(IR)**NVARTRAN
             ENDDO
          ELSE
             DO IR=1, NR
                TG_LOC(IR)=REAL(IR,q)*DR
                RG_LOC(IR)=TG_LOC(IR)
             ENDDO
          ENDIF
! If IBROAD /= 1, execute following procedure.
          IF(IBROAD .NE. 1) THEN
! Calculate metric function on the numerical integration grid.
             IF(LMETRIC) THEN
                CALL PMET0 (FMET,NMETRIC,NR,RG_LOC,RCUT,RMETRIC)
             ELSE
                FMET(1:NR)=1.0_q
             ENDIF
! Calculate the target ion distribution function on the spline grids.
! Because the modified spherical Bessel function easily diverse at large RG,
! the exponential functions that will be multiplied later are introduced 
! in the subroutine.
! Hence, the calculated modified spherical Bessel functions are not the
! original ones but the ones multiplied with exponential functions.
! Loop over spline grid.
             DO ISPL=1, NSPL
! Set numericl grids to compute modified spherical Bessel functions.
                RR(1:NR)=2.0_q*WION*FNL(ISPL,1,0,1)*RG_LOC(1:NR)
                RRI=WION*FNL(ISPL,1,0,1)**2
                RRJ(1:NR)=WION*RG_LOC(1:NR)**2
! Calculate modified spherical Bessel functions.
                CALL IL0 (LMAX,NR,NR,RR,RRI,RRJ,RIL)
! Execute numerical integrations.
                DO L=0, LMAX
                   DO IR=1, NR
                      CR(IR,L)=PREFAC*FCUT(ISPL)*FMET(IR)*RIL(IR,L)*DR*TG_LOC(IR)**2
                   ENDDO
                   DO N=1, NRB(L)
                      FNL(ISPL,2,L,N)=DDOT(NR,RB(1,N,L),1,CR(1,L),1)
                   ENDDO
                ENDDO
             ENDDO
! Calculate first derivative at the first grid.
! In this case, derivatives of original modified spherical Bessel functions are calculated,
! because their values at small RG only are necessary.
! Calculate Gauss function, exp(-RG_LOC^2/2/WION^2), on the numerical integration grid.
             CALL GAUSS0 (EXP21,NR,RG_LOC,WION)
! Calculate Gauss function, exp(-FNL^2/2/WION^2), on the spline grids.
             CALL GAUSS0 (EXP22,NSPL,FNL(1:,1,0,1),WION)
! Calculate derivative of Gauss function, exp(-RG_LOC^2/2/WION^2), on the numerical integration grid.
             CALL GAUSS1 (DEXP22(1:),EXP22(1:),1,FNL(1:,1,0,1),WION)
! Set radial component, FNL*RG_LOC.
             RR(1:NR)=2.0_q*WION*FNL(1,1,0,1)*RG_LOC(1:NR)
! Calculate derivative of modified spherical Bessel functions.
             CALL IL1 (DRIL,LMAX,NR,NR,RR,RIL)
! Then, execute the numerical integration over the (transformed) radial grid.
             DO L=0, LMAX
                DO IR=1, NR
                   DCR(IR,L)=PREFAC*(DFCUT(1)* FMET(IR)* EXP22(1)*EXP21(IR)* RIL(IR,L)+ &
                                      FCUT(1)* FMET(IR)*DEXP22(1)*EXP21(IR)* RIL(IR,L)+ &
                                      FCUT(1)* FMET(IR)* EXP22(1)*EXP21(IR)*DRIL(IR,L)*2.0_q*WION*RG_LOC(IR))* &
                                      DR*TG_LOC(IR)**2
                ENDDO
                DO N=1, NRB(L)
                   DFNL_LOC(N,L)=DDOT(NR,RB(1,N,L),1,DCR(1,L),1)
                ENDDO
             ENDDO
! Otherwise, simply interpolate the radial basis functions.
! This means that the atomic distribution is described by a sum non-broadened delta-functions.
! In this case, variable transformation is not executed.
          ELSE
! Calculate metric function
             IF(LMETRIC) THEN
                CALL PMET0 ( FMET   ,NMETRIC,NSPL,FNL(1:,1,0,1),RCUT,RMETRIC)
                CALL PMET1 (DFMET(1:),NMETRIC,   1,FNL(1:,1,0,1),RCUT,RMETRIC)
             ELSE
                FMET(1:NSPL) =1.0_q
                DFMET(1)=0.0_q
             ENDIF
! Calculate spherical Bessel functions on the spline grids
             DO L=0, LMAX
                DO N=1, NRB(L)
                   DO ISPL=1, NSPL
                      X=QGRID(N,L)*FNL(ISPL,1,L,N)
                      CALL SBESSEL(X,FNL(ISPL,2,L,N),L)
                      FNL(ISPL,2,L,N)=FCUT(ISPL)*FMET(ISPL)*FNL(ISPL,2,L,N)/RNORM(L,N)
                   ENDDO
                ENDDO
             ENDDO
             DO L=0, LMAX
                DO N=1, NRB(L)
                   X=QGRID(N,L)*FNL(1,1,L,N)
                   CALL SBESSE2(X,DUMMY,DFNL_LOC(N,L),L)
                   DFNL_LOC(N,L)=((DFCUT(1)*FMET(1)+FCUT(1)*DFMET(1))*DUMMY+QGRID(N,L)*FCUT(1)*FMET(1)*DFNL_LOC(N,L))/RNORM(L,N)
                ENDDO
             ENDDO
          ENDIF
! Determine the coefficients of spline interpolation
          DO L=0, LMAX
             DO N=1, NRB(L)
                CALL SPLCOF_ML (FNL(1:,1:,L,N),NSPL,NSPL,DFNL_LOC(N,L))
             ENDDO
          ENDDO
        END SUBROUTINE RAD_FUNC

! Subroutine calculating Gaussian function

       SUBROUTINE GAUSS0 (GAUSS,N,R,W)
         IMPLICIT NONE
! Input variables
         INTEGER, INTENT(IN) :: N
         REAL(q), INTENT(IN) :: R(:) !(1:N)
         REAL(q), INTENT(IN) :: W
! Output variable
         REAL(q), INTENT(OUT) :: GAUSS(:) !(1:N)
! Local variable
         INTEGER I
         DO I=1, N
            GAUSS(I)=EXP(-W*R(I)**2)
         ENDDO
       END SUBROUTINE GAUSS0

! Subroutine calculating 1st derivative of Gaussian function

       SUBROUTINE GAUSS1 (DGAUSS,GAUSS,N,R,W)
         IMPLICIT NONE
! Input variables
         INTEGER, INTENT(IN) :: N
         REAL(q), INTENT(IN) :: R(:) !(1:N)
         REAL(q), INTENT(IN) :: W
! Output variables
         REAL(q) DGAUSS(:) !(1:N)
         REAL(q) GAUSS(:) !(1:N)
! Local variable
         INTEGER I
         DO I=1, N
            GAUSS(I)=EXP(-W*R(I)**2)
            DGAUSS(I)=-2.0_q*W*R(I)*GAUSS(I)
         ENDDO
       END SUBROUTINE GAUSS1

! Subroutine for calculating spherical harmonics

       SUBROUTINE YLM0 (LMAX,MNSF,R0,XYZ0,YLM)
         IMPLICIT NONE
! Input variables
         INTEGER, INTENT(IN) :: LMAX
         INTEGER, INTENT(IN) :: MNSF
         REAL(q), INTENT(IN) :: XYZ0(:) !(1:3) ! Input geometry data
! Output variables
         REAL(q), INTENT(OUT) :: R0
         REAL(q), INTENT(OUT) :: YLM(:) !(1:MNSF) ! YLM
! Temporary data effective only in this subroutine
         REAL(q) EPS
         REAL(q) FAC(1:2*LMAX-1)
         REAL(q) FACTOR
         INTEGER I
         INTEGER MM,MMM
         INTEGER L
         INTEGER LL
         INTEGER LM,LM1,LM2,LM3,LMP,LMM
         REAL(q) PLM0
         REAL(q) RFPI
         REAL(q) RXY1(0:LMAX-1)
         REAL(q) SGN
         REAL(q) WW1,ZZ1
         REAL(q) WW2,ZZ2
         REAL(q) WW3,ZZ3
         REAL(q) XYZ1(1:3)
         REAL(q) ZR1,ZI1
! Constant
         EPS=1E-10_q
         RFPI=1.0_q/SQRT(4.0_q*PI)
! Spherical harmonics with l=0 and m=0
         IF(LMAX.LE.0) THEN
            RETURN
         ENDIF
! Factorial
         CALL FACTORIAL (FAC,2*LMAX-1)
         IF(LMAX.EQ.1) THEN
            R0=SQRT(XYZ0(1)*XYZ0(1)+XYZ0(2)*XYZ0(2)+XYZ0(3)*XYZ0(3))
            YLM(1)=RFPI
            RETURN
         ENDIF
! Angular conmponents
         R0=SQRT(XYZ0(1)*XYZ0(1)+XYZ0(2)*XYZ0(2)+XYZ0(3)*XYZ0(3))
         RXY1(0)=1.0_q
         IF(R0.LT.EPS) THEN
            XYZ1(1)=0.0_q
            XYZ1(2)=0.0_q
            XYZ1(3)=1.0_q
            RXY1(1)=0.0_q
         ELSE
            RXY1(1)=SQRT(XYZ0(1)*XYZ0(1)+XYZ0(2)*XYZ0(2))
            IF(RXY1(1)/R0.LT.EPS) THEN
               XYZ1(1)=1.0_q
               XYZ1(2)=0.0_q
               XYZ1(3)=SIGN(1.0_q,XYZ0(3))
               RXY1(1)=0.0_q
            ELSE
               XYZ1(1)=XYZ0(1)/RXY1(1)
               XYZ1(2)=XYZ0(2)/RXY1(1)
               XYZ1(3)=XYZ0(3)/R0
               RXY1(1)=RXY1(1)/R0
            ENDIF
         ENDIF
         IF(LMAX.GE.3) THEN
            DO L=2, LMAX-1
               RXY1(L)=RXY1(L-1)*RXY1(1)
            ENDDO
         ENDIF
! Derivatives of Legendre polynomials
! Calculation of Pl,m, sin(theta) * dPl,m/dx, sin2(theta) * d2Pl,m/dx2
! (l,m)=(0,0),(1,1)
         YLM(1)=1.0_q
         YLM(3)=-RXY1(1)
         IF(LMAX.GE.3) THEN
            YLM(7)=3.0_q*RXY1(2)
         ENDIF
! Calculation of sin(theta) * dPl,m/dx, sin2(theta) * d2Pl,m/dx2
! (l,m)=(l,l), l is larger than 2
         IF(LMAX.GE.4) THEN
            DO L=3, LMAX-1
               LM1=L**2+L+1
               LM2=(L-1)**2+L
               SGN=(-1.0_q)**MOD(L+1,2)
               YLM(LM1)=-REAL((2*L-1),q)*RXY1(1)*YLM(LM2)
            ENDDO
         ENDIF
! Calculation of sin(theta) * dPl,m/dx, sin2(theta) * d2Pl,m/dx2
! (l,m)=(l,l-1), l is larger than 0
         DO L=1, LMAX-1
            LM1=L**2+L
            LM2=(L-1)**2+L
            YLM(LM1)=XYZ1(3)*REAL((2*L-1),q)*YLM(LM2)
         ENDDO
! Calculation of sin(theta) * dPl,m/dx, sin2(theta) * d2Pl,m/dx2
! Remaining (l,m)
         DO MM=0, LMAX-3
            DO L=MM+2, LMAX-1
               LM1=L**2+MM+1
               LM2=(L-1)**2+MM+1
               LM3=(L-2)**2+MM+1
               YLM(LM1)=(XYZ1(3)*REAL((2*L-1),q)*YLM(LM2) &
                       -REAL((L+MM-1),q)*YLM(LM3))/REAL((L-MM),q)
            ENDDO
         ENDDO
         ZZ1=1.0_q
         WW1=0.0_q
! m=0 part
         DO LL=1, LMAX
            L=LL-1
            FACTOR=SQRT(REAL((2*L+1),q))*RFPI
            LM=L**2+1
            PLM0=YLM(LM)
            YLM(LM)=FACTOR*PLM0
         ENDDO
         DO MM=2, LMAX
            ZR1=ZZ1
            ZI1=WW1
            ZZ1=ZR1*XYZ1(1)-ZI1*XYZ1(2)
            WW1=ZR1*XYZ1(2)+ZI1*XYZ1(1)
            DO LL=MM, LMAX
               L=LL-1
               MMM=MM-1
               FACTOR=SQRT(REAL((2*L+1),q)*FAC(L-MMM+1)/(2.0_q*PI*FAC(L+MMM+1)))
               LMP=L**2+MM
               LMM=LMP+L
               SGN=(-1.0_q)**(MMM+1)
               PLM0=YLM(LMP)
               YLM(LMP)=-SGN*FACTOR*PLM0*ZZ1
               YLM(LMM)=-SGN*FACTOR*PLM0*WW1
            ENDDO
         ENDDO
       END SUBROUTINE YLM0

! 1st derivatives of real spherical harmonics

       SUBROUTINE YLM1 (EPS,LMAX,MNSF,XYZ0,COSD,COSP,R0,SIND,SIND_CRITERIA,SINP,YLM,YLMD,YLMP)
         IMPLICIT NONE
! Input variables
         REAL(q), INTENT(IN) :: EPS
         INTEGER, INTENT(IN) :: LMAX
         INTEGER, INTENT(IN) :: MNSF
         REAL(q), INTENT(IN) :: XYZ0(:) !(1:3) ! Input geometry data
! Output variables
         REAL(q), INTENT(OUT) :: COSD,COSP ! Output
         REAL(q), INTENT(OUT) :: R0
         REAL(q), INTENT(OUT) :: SIND,SIND_CRITERIA,SINP ! Output
         REAL(q), INTENT(OUT) :: YLM(:) !(1:MNSF) ! Output Ylm
         REAL(q), INTENT(OUT) :: YLMD(:) !(1:MNSF) ! Output dYlm/d(theta) or dYlm/dx|theta=0
         REAL(q), INTENT(OUT) :: YLMP(:) !(1:MNSF) ! Output dYlm/d(phi) or dYlm/dy|theta=0
! Temporary data effective only in this subroutine
         REAL(q) FAC(1:2*LMAX-1),FAC2(1:2*LMAX-1)
         REAL(q) FACTOR
         INTEGER I
         INTEGER L,LL
         INTEGER LM,LM1,LM2,LM3,LMM,LMP
         INTEGER MM,MMM
         REAL(q) PLM0,PLM1
         REAL(q) RFPI
         REAL(q) RXY1(0:LMAX-1)
         REAL(q) SGN
         REAL(q) WW1,ZZ1
         REAL(q) WW2,ZZ2
         REAL(q) WW3,ZZ3
         REAL(q) XYZ1(1:3),XYZ2(1:3),XYZ3(1:3)
         REAL(q) ZI1,ZI2,ZI3
         REAL(q) ZR1,ZR2,ZR3
! Constant
         RFPI=1.0_q/SQRT(4.0_q*PI)
! Spherical harmonics with l=0 and m=0
         IF(LMAX .LE. 0) THEN
            RETURN
         ENDIF
! Factorial
         CALL FACTORIAL (FAC,2*LMAX-1)
         CALL FACTORIAL2 (FAC2,2*LMAX-1)
! Angular conmponents
         IF(LMAX .EQ. 1) THEN
            YLM(1)=RFPI
            YLMD(1)=0.0_q
            YLMP(1)=0.0_q
            R0=SQRT(XYZ0(1)*XYZ0(1)+XYZ0(2)*XYZ0(2)+XYZ0(3)*XYZ0(3))
            IF(R0 .LT. EPS) THEN
               COSD=0.0_q
               SIND=0.0_q
               COSP=0.0_q
               SINP=0.0_q
            ELSE
               RXY1(0)=SQRT(XYZ0(1)*XYZ0(1)+XYZ0(2)*XYZ0(2))
               SIND_CRITERIA=RXY1(0)/R0
               IF(RXY1(0)/R0 .LT. EPS) THEN
                  XYZ1(3)=SIGN(1.0_q,XYZ0(3))
                  COSD=XYZ1(3)
                  SIND=0.0_q
                  COSP=0.0_q
                  SINP=0.0_q
               ELSE
                  XYZ1(1)=XYZ0(1)/RXY1(0)
                  XYZ1(2)=XYZ0(2)/RXY1(0)
                  XYZ1(3)=XYZ0(3)/R0
                  RXY1(0)=RXY1(0)/R0
                  COSD=XYZ1(3)
                  SIND=RXY1(0)
                  COSP=XYZ1(1)
                  SINP=XYZ1(2)
               ENDIF
            ENDIF
            RETURN
         ENDIF
         R0=SQRT(XYZ0(1)*XYZ0(1)+XYZ0(2)*XYZ0(2)+XYZ0(3)*XYZ0(3))
         RXY1(0)=1.0_q
         IF(R0 .LT. EPS) THEN
            XYZ1(1)=0.0_q
            XYZ1(2)=0.0_q
            XYZ1(3)=1.0_q
            XYZ2(1)=XYZ1(1)
            XYZ2(2)=XYZ1(2)
            XYZ2(3)=XYZ1(3)
            XYZ3(1)=XYZ1(1)
            XYZ3(2)=XYZ1(2)
            XYZ3(3)=XYZ1(3)
            RXY1(1)=0.0_q
            COSD=0.0_q
            SIND=0.0_q
            COSP=0.0_q
            SINP=0.0_q
         ELSE
            RXY1(1)=SQRT(XYZ0(1)*XYZ0(1)+XYZ0(2)*XYZ0(2))
            SIND_CRITERIA=RXY1(1)/R0
            IF(RXY1(1)/R0 .LT. EPS) THEN
               XYZ1(1)=1.0_q
               XYZ1(2)=0.0_q
               XYZ1(3)=SIGN(1.0_q,XYZ0(3))
               XYZ2(1)=XYZ1(2)
               XYZ2(2)=XYZ1(1)
               XYZ2(3)=XYZ1(3)
               XYZ3(1)=1.0_q/SQRT(2.0_q)
               XYZ3(2)=1.0_q/SQRT(2.0_q)
               XYZ3(3)=XYZ1(3)
               RXY1(1)=0.0_q
               COSD=XYZ1(3)
               SIND=0.0_q
               COSP=0.0_q
               SINP=0.0_q
            ELSE
               XYZ1(1)=XYZ0(1)/RXY1(1)
               XYZ1(2)=XYZ0(2)/RXY1(1)
               XYZ1(3)=XYZ0(3)/R0
               XYZ2(1)=XYZ1(1)
               XYZ2(2)=XYZ1(2)
               XYZ2(3)=XYZ1(3)
               XYZ3(1)=XYZ1(1)
               XYZ3(2)=XYZ1(2)
               XYZ3(3)=XYZ1(3)
               RXY1(1)=RXY1(1)/R0
               COSD=XYZ1(3)
               SIND=RXY1(1)
               COSP=XYZ1(1)
               SINP=XYZ1(2)
            ENDIF
         ENDIF
         IF(LMAX .GE. 3) THEN
            DO L=2, LMAX-1
               RXY1(L)=RXY1(L-1)*RXY1(1)
            ENDDO
         ENDIF
! Derivatives of Legendre polynomials
! Calculation of Pl,m, sin(theta) * dPl,m/dx, sin2(theta) * d2Pl,m/dx2
! (l,m)=(0,0),(1,1)
         YLM(1)=1.0_q
         YLMD(1)=0.0_q
         YLM(3)=-RXY1(1)
         YLMD(3)=XYZ1(3)
         IF(LMAX .GE. 3) THEN
            YLM(7)=3.0_q*RXY1(2)
            YLMD(7)=-6.0_q*RXY1(1)*XYZ1(3)
         ENDIF
! Calculation of sin(theta) * dPl,m/dx, sin2(theta) * d2Pl,m/dx2
! (l,m)=(l,l), l is larger than 2
         IF(LMAX .GE. 4) THEN
            DO L=3, LMAX-1
               LM1=L**2+L+1
               LM2=(L-1)**2+L
               SGN=(-1.0_q)**MOD(L+1,2)
               YLM(LM1)=-REAL((2*L-1),q)*RXY1(1)*YLM(LM2)
               YLMD(LM1)=SGN*REAL(L,q)*FAC2(L)*RXY1(L-1)*XYZ1(3)
            ENDDO
         ENDIF
! Calculation of sin(theta) * dPl,m/dx, sin2(theta) * d2Pl,m/dx2
! (l,m)=(l,l-1), l is larger than 0
         DO L=1, LMAX-1
            LM1=L**2+L
            LM2=(L-1)**2+L
            YLM(LM1)=XYZ1(3)*REAL((2*L-1),q)*YLM(LM2)
            YLMD(LM1)=REAL((2*L-1),q)*YLM(LM2)*RXY1(1) &
                     +XYZ1(3)*REAL((2*L-1),q)*YLMD(LM2)
         ENDDO
! Calculation of sin(theta) * dPl,m/dx, sin2(theta) * d2Pl,m/dx2
! Remaining (l,m)
         DO MM=0, LMAX-3
            DO L=MM+2, LMAX-1
               LM1=L**2+MM+1
               LM2=(L-1)**2+MM+1
               LM3=(L-2)**2+MM+1
               YLM(LM1)=(XYZ1(3)*REAL((2*L-1),q)*YLM(LM2) &
                      -REAL((L+MM-1),q)*YLM(LM3))/REAL((L-MM),q)
               YLMD(LM1)=(REAL((2*L-1),q)*YLM(LM2)*RXY1(1) &
                          +XYZ1(3)*REAL((2*L-1),q)*YLMD(LM2) &
                          -REAL((L+MM-1),q)*YLMD(LM3))/REAL((L-MM),q)
            ENDDO
         ENDDO
! Gradients of spherical harmonics
! m=0 part
         ZZ1=1.0_q
         WW1=0.0_q
         ZZ2=1.0_q
         WW2=0.0_q
         ZZ3=1.0_q
         WW3=0.0_q
         DO LL=1, LMAX
            L=LL-1
            FACTOR=SQRT(REAL((2*L+1),q))*RFPI
            LM=L**2+1 
            PLM0=YLM(LM)
            PLM1=YLMD(LM)
            YLM(LM)=FACTOR*PLM0
            YLMD(LM)=-FACTOR*PLM1
            YLMP(LM)=0.0_q
            IF(ABS(RXY1(1)) .LE. EPS) THEN
               YLMD(LM)=YLMD(LM)/R0
               YLMP(LM)=YLMD(LM)
            ENDIF
         ENDDO
         DO MM=2, LMAX
            ZR1=ZZ1
            ZI1=WW1
            ZR2=ZZ2
            ZI2=WW2
            ZR3=ZZ3
            ZI3=WW3
            ZZ1=ZR1*XYZ1(1)-ZI1*XYZ1(2)
            WW1=ZR1*XYZ1(2)+ZI1*XYZ1(1)
            ZZ2=ZR2*XYZ2(1)-ZI2*XYZ2(2)
            WW2=ZR2*XYZ2(2)+ZI2*XYZ2(1)
            ZZ3=ZR3*XYZ3(1)-ZI3*XYZ3(2)
            WW3=ZR3*XYZ3(2)+ZI3*XYZ3(1)
            DO LL=MM, LMAX
               L=LL-1
               MMM=MM-1
               FACTOR=SQRT(REAL((2*L+1),q)*FAC(L-MMM+1)/(2.0_q*PI*FAC(L+MMM+1)))
               LMP=L**2+MM
               LMM=LMP+L
               SGN=(-1.0_q)**(MMM+1)
               PLM0=YLM(LMP)
               PLM1=YLMD(LMP)
               YLM(LMP)=-SGN*FACTOR*PLM0*ZZ1
               YLM(LMM)=-SGN*FACTOR*PLM0*WW1
               YLMD(LMP)= SGN*FACTOR*PLM1*ZZ1
               YLMD(LMM)= SGN*FACTOR*PLM1*WW1
               IF(ABS(RXY1(1)) .LE. EPS) THEN
                  YLMD(LMP)=YLMD(LMP)/R0
                  YLMD(LMM)=YLMD(LMM)/R0
                  YLMP(LMP)=SGN*FACTOR*PLM1*ZZ2/R0
                  YLMP(LMM)=SGN*FACTOR*PLM1*WW2/R0
               ELSE
                  YLMP(LMP)=-REAL(MMM,q)*YLM(LMM)
                  YLMP(LMM)= REAL(MMM,q)*YLM(LMP)
               ENDIF
            ENDDO
         ENDDO
! Stop profiling
         RETURN
       END SUBROUTINE YLM1

! Factorial

       SUBROUTINE FACTORIAL (F,N)
         IMPLICIT NONE
! Input variable
         INTEGER, INTENT(IN) :: N
! Output variable
         REAL(q), INTENT(OUT) :: F(:) !(N)
! Local variable
         INTEGER I
         F(1)=1.0_q
         IF(N .GT. 1) THEN
            DO I=2, N
               F(I)=REAL((I-1),q)*F(I-1)
            ENDDO
         ENDIF
         RETURN
       END SUBROUTINE FACTORIAL

       SUBROUTINE FACTORIAL2(F,N)
         IMPLICIT NONE
! Input variable
         INTEGER, INTENT(IN) :: N
! Output variable
         REAL(q), INTENT(OUT) :: F(:) !(N)
         INTEGER I
         F(1)=1.0_q
         IF(N .GT. 1) THEN
            DO I=2, N
               F(I)=REAL((2*I-1),q)*F(I-1)
            ENDDO
         ENDIF
         RETURN
       END SUBROUTINE FACTORIAL2

! Behler-Parinello cutoff function

       SUBROUTINE BP0 (FCUT,N,RCUT,R)
         IMPLICIT NONE
! Input variables
         INTEGER, INTENT(IN) :: N
         REAL(q), INTENT(IN) :: R(:) !(1:N)
         REAL(q), INTENT(IN) :: RCUT
! Output variable
         REAL(q) FCUT(:) !(1:N)
! Local variable
         INTEGER I
         DO I=1, N
            IF(R(I) .LE. RCUT) THEN
               FCUT(I)=0.5_q*(DCOS(PI*R(I)/RCUT)+1.0_q)
            ELSE
               FCUT(I)=0.0_q
            ENDIF
         ENDDO
       END SUBROUTINE BP0

! 1st derivative of Behler-Parinello cutoff function

       SUBROUTINE BP1 (FCUT,DFCUT,N,RCUT,R)
         IMPLICIT NONE
! Input variables
         INTEGER, INTENT(IN) :: N
         REAL(q), INTENT(IN) :: R(:) !(1:N)
         REAL(q), INTENT(IN) :: RCUT
! Output variables
         REAL(q), INTENT(OUT) :: FCUT(:) !(1:N)
         REAL(q), INTENT(OUT) :: DFCUT(:) !(1:N)
! Local variable
         INTEGER I
         DO I=1, N
            IF(R(I) .LE. RCUT) THEN
               FCUT(I)=0.5_q*(DCOS(PI*R(I)/RCUT)+1.0_q)
               DFCUT(I)=0.5_q*(-PI/RCUT*DSIN(PI*R(I)/RCUT))
            ELSE
               FCUT(I)=0.0_q
               DFCUT(I)=0.0_q
            ENDIF
         ENDDO
       END SUBROUTINE BP1

! Cutoff function proposed by Miwa.

       SUBROUTINE MO0 (FCUT,N,RCUT,R)
         IMPLICIT NONE
! Input variables
         INTEGER, INTENT(IN) :: N
         REAL(q), INTENT(IN) :: RCUT
         REAL(q), INTENT(IN) :: R(:) !(1:N)
! Output variables
         REAL(q), INTENT(OUT) :: FCUT(:) !(1:N)
! Local variables
         INTEGER :: I
         REAL(q) :: X
         DO I=1, N
            X=4.0_q*R(I)/RCUT-3.0_q
            IF(X .LT. -1.0_q) THEN
               FCUT(I)=1.0_q
            ELSE IF(X .LT. 1.0_q) THEN
               FCUT(I)=0.25_q*(X**3-3.0_q*X+2.0_q)
            ELSE
               FCUT(I)=0.0_q
            ENDIF
         ENDDO
       END SUBROUTINE MO0

! 1st dericative of the cutoff function proposed by Miwa.

       SUBROUTINE MO1 (FCUT,DFCUT,N,RCUT,R)
         IMPLICIT NONE
! Input variables
         INTEGER, INTENT(IN) :: N
         REAL(q), INTENT(IN) :: RCUT
         REAL(q), INTENT(IN) :: R(:) !(1:N)
! Output variables
         REAL(q), INTENT(OUT) :: DFCUT(:) !(1:N)
         REAL(q), INTENT(OUT) :: FCUT(:) !(1:N)
! Local variables
         REAL(q) :: DXDR
         INTEGER :: I
         REAL(q) :: X
         DO I=1, N
            X=4.0_q*R(I)/RCUT-3.0_q
            DXDR=4.0_q/RCUT
            IF(X .LT. -1.0_q) THEN
               FCUT(I)=1.0_q
               DFCUT(I)=0.0_q
            ELSE IF(X .LT. 1.0_q) THEN
               FCUT(I)=0.25_q*(X**3-3.0_q*X+2.0_q)
               DFCUT(I)=0.75_q*(X**2-1.0_q)*DXDR
            ELSE
               FCUT(I)=0.0_q
               DFCUT(I)=0.0_q
            ENDIF
         ENDDO
       END SUBROUTINE MO1

! Cutoff function which gives more flat function within the cutoff radius.
! This might be better to be used with the metric method to make the meaning of metric clear.

       SUBROUTINE RJ0 (FCUT,N,RCUT,R)
         IMPLICIT NONE
! Input variables
         INTEGER, INTENT(IN) :: N
         REAL(q), INTENT(IN) :: RCUT
         REAL(q), INTENT(IN) :: R(:) !(1:N)
! Output variables
         REAL(q), INTENT(OUT) :: FCUT(:) !(1:N)
! Local variables
         REAL(q) :: DELTA
         INTEGER :: I
         REAL(q) :: RCUT1
         REAL(q) :: RCUT2
         REAL(q) :: X
         DELTA=0.5_q/AUTOA
         IF(RCUT-DELTA.GT.0.0_q) THEN
            RCUT2=RCUT
            RCUT1=RCUT2-DELTA
         ELSE
            RCUT2=RCUT
            RCUT1=0.5_q*RCUT2
         ENDIF
         DO I=1, N
            IF(R(I).LE.RCUT1) THEN
               FCUT(I)=1.0_q
            ELSE IF(R(I).LE.RCUT2) THEN
               FCUT(I)=(RCUT2**2-R(I)**2)**2*(RCUT2**2+2.0_q*R(I)**2-3.0_q*RCUT1**2)/(RCUT2**2-RCUT1**2)**3
            ELSE
               FCUT(I)=0.0_q
            ENDIF
         ENDDO
       END SUBROUTINE RJ0

! 1st dericative of the cutoff function proposed by R. Jinnouchi.

       SUBROUTINE RJ1 (FCUT,DFCUT,N,RCUT,R)
         IMPLICIT NONE
! Input variables
         INTEGER, INTENT(IN) :: N
         REAL(q), INTENT(IN) :: RCUT
         REAL(q), INTENT(IN) :: R(:) !(1:N)
! Output variables
         REAL(q), INTENT(OUT) :: DFCUT(:) !(1:N)
         REAL(q), INTENT(OUT) :: FCUT(:) !(1:N)
! Local variables
         REAL(q) :: DELTA
         REAL(q) :: DXDR
         INTEGER :: I
         REAL(q) :: RCUT1
         REAL(q) :: RCUT2
         REAL(q) :: X
         DELTA=0.5_q/AUTOA
         IF(RCUT-DELTA.GT.0.0_q) THEN
            RCUT2=RCUT
            RCUT1=RCUT2-DELTA
         ELSE
            RCUT2=RCUT
            RCUT1=0.5_q*RCUT2
         ENDIF
         DO I=1, N
            IF(R(I).LE.RCUT1) THEN
               FCUT(I)=1.0_q
               DFCUT(I)=0.0_q
            ELSE IF(R(I).LE.RCUT2) THEN
               FCUT(I)=(RCUT2**2-R(I)**2)**2*(RCUT2**2+2.0_q*R(I)**2-3.0_q*RCUT1**2)/(RCUT2**2-RCUT1**2)**3
               DFCUT(I)=12.0_q*R(I)*(RCUT2**2-R(I)**2)*(RCUT1**2-R(I)**2)/(RCUT2**2-RCUT1**2)**3
            ELSE
               FCUT(I)=0.0_q
               DFCUT(I)=0.0_q
            ENDIF
         ENDDO
       END SUBROUTINE RJ1

! Cutoff function proposed in M. J. Willatt, F. Musil and M. Ceriotti. PCCP, 20, 29661 (2018).
! This corresponds to the radial function used in the optimized radially scaling kernel.

       SUBROUTINE WMC0 (FCUT,N,RCUT,R)
         IMPLICIT NONE
! Input variables
         INTEGER, INTENT(IN)  :: N
         REAL(q), INTENT(IN)  :: RCUT
         REAL(q), INTENT(IN)  :: R(:) !(1:N)
! Output variables
         REAL(q), INTENT(OUT) :: FCUT(:) !(1:N)
! Local variables
         REAL(q)              :: C
         REAL(q)              :: DELTA
         REAL(q)              :: FCUT1
         REAL(q)              :: FCUT2
         INTEGER              :: I
         INTEGER              :: M
         REAL(q)              :: R0
         REAL(q)              :: RCUT1
         REAL(q)              :: RCUT2
         REAL(q)              :: X
! Setting the parameters.
         C=1.0_q
         DELTA=0.5_q/AUTOA
         M=7
         R0=2.0_q/AUTOA
         IF(RCUT-DELTA.GT.0.0_q) THEN
            RCUT2=RCUT
            RCUT1=RCUT2-DELTA
         ELSE
            RCUT2=RCUT
            RCUT1=0.5_q*RCUT2
         ENDIF
! Calculate the cufoff function.
         DO I=1, N
            IF(R(I).LE.RCUT1) THEN
               FCUT1=1.0_q
            ELSE IF(R(I).LE.RCUT2) THEN
               FCUT1=(RCUT2**2-R(I)**2)**2*(RCUT2**2+2.0_q*R(I)**2-3.0_q*RCUT1**2)/(RCUT2**2-RCUT1**2)**3
            ELSE
               FCUT1=0.0_q
            ENDIF
            FCUT2=C/(C+(R(I)/R0)**M)
            FCUT(I)=FCUT1*FCUT2
         ENDDO
       END SUBROUTINE WMC0

! Derivative of the cutoff function proposed in M. J. Willatt, F. Musil and M. Ceriotti. PCCP, 20, 29661 (2018).

       SUBROUTINE WMC1 (FCUT,DFCUT,N,RCUT,R)
         IMPLICIT NONE
! Input variables
         INTEGER, INTENT(IN) :: N
         REAL(q), INTENT(IN) :: RCUT
         REAL(q), INTENT(IN) :: R(:) !(1:N)
! Output variables
         REAL(q), INTENT(OUT) :: DFCUT(:) !(1:N)
         REAL(q), INTENT(OUT) :: FCUT(:) !(1:N)
! Local variables
         REAL(q) :: C
         REAL(q) :: DFCUT1
         REAL(q) :: DFCUT2
         REAL(q) :: DELTA
         REAL(q) :: DXDR
         REAL(q) :: FCUT1
         REAL(q) :: FCUT2
         INTEGER :: I
         INTEGER :: M
         REAL(q) :: R0
         REAL(q) :: RCUT1
         REAL(q) :: RCUT2
         REAL(q) :: X
! Set constants.
         C=1.0_q
         DELTA=0.5_q/AUTOA
         M=7
         R0=2.0_q/AUTOA
         IF(RCUT-DELTA.GT.0.0_q) THEN
            RCUT2=RCUT
            RCUT1=RCUT2-DELTA
         ELSE
            RCUT2=RCUT
            RCUT1=0.5_q*RCUT2
         ENDIF
! Calculate cutoff function.
         DO I=1, N
            IF(R(I).LE.RCUT1) THEN
               FCUT1=1.0_q
               DFCUT1=0.0_q
            ELSE IF(R(I).LE.RCUT2) THEN
               FCUT1=(RCUT2**2-R(I)**2)**2*(RCUT2**2+2.0_q*R(I)**2-3.0_q*RCUT1**2)/(RCUT2**2-RCUT1**2)**3
               DFCUT1=12.0_q*R(I)*(RCUT2**2-R(I)**2)*(RCUT1**2-R(I)**2)/(RCUT2**2-RCUT1**2)**3
            ELSE
               FCUT1=0.0_q
               DFCUT1=0.0_q
            ENDIF
            FCUT2=C/(C+(R(I)/R0)**M)
            DFCUT2=-C/(C+(R(I)/R0)**M)**2*M*(R(I)/R0)**(M-1)/R0
            FCUT(I)=FCUT1*FCUT2
            DFCUT(I)=DFCUT1*FCUT2+DFCUT2*FCUT1
         ENDDO
       END SUBROUTINE WMC1

! Polynomical radial metric function.
! Cutoff function is introduced also here because the resulted distribution function does not approach zero at R=RCUT.

       SUBROUTINE PMET0 (FMET,NMETRIC,N,R,RCUT,RMETRIC)
         IMPLICIT NONE
! Input variables
         INTEGER, INTENT(IN)  :: NMETRIC
         INTEGER, INTENT(IN)  :: N
         REAL(q), INTENT(IN)  :: R(:) !(1:N)
         REAL(q), INTENT(IN)  :: RCUT
         REAL(q), INTENT(IN)  :: RMETRIC
! Output variables
         REAL(q), INTENT(OUT) :: FMET(:) !(1:N)
! Local variables
         REAL(q)              :: FMET_HELP
         REAL(q)              :: FCUT_HELP
         INTEGER              :: I
         REAL(q)              :: X
         DO I=1, N
            X=4.0_q*R(I)/RCUT-3.0_q
            IF(X .LT. -1.0_q) THEN
               FCUT_HELP=1.0_q
            ELSE IF(X .LT. 1.0_q) THEN
               FCUT_HELP=0.25_q*(X**3-3.0_q*X+2.0_q)
            ELSE
               FCUT_HELP=0.0_q
            ENDIF
            FMET_HELP=(RMETRIC/(RMETRIC+R(I)))**NMETRIC
            FMET(I)=FMET_HELP*FCUT_HELP
!            FMET(I)=FMET_HELP
         ENDDO
       END SUBROUTINE PMET0

! Derivative of polynomical radial metric function

       SUBROUTINE PMET1 (DFMET,NMETRIC,N,R,RCUT,RMETRIC)
         IMPLICIT NONE
! Input variables
         INTEGER, INTENT(IN)  :: NMETRIC
         INTEGER, INTENT(IN)  :: N
         REAL(q), INTENT(IN)  :: R(:) !(1:N)
         REAL(q), INTENT(IN)  :: RCUT
         REAL(q), INTENT(IN)  :: RMETRIC
! Output variables
         REAL(q), INTENT(OUT) :: DFMET(:) !(1:N)
! Local variables
         REAL(q)              :: DFCUT_HELP
         REAL(q)              :: DFMET_HELP
         REAL(q)              :: DXDR
         REAL(q)              :: FCUT_HELP
         REAL(q)              :: FMET_HELP
         INTEGER              :: I
         REAL(q)              :: X
         DO I=1, N
            X=4.0_q*R(I)/RCUT-3.0_q
            DXDR=4.0_q/RCUT
            IF(X .LT. -1.0_q) THEN
               FCUT_HELP=1.0_q
               DFCUT_HELP=0.0_q
            ELSE IF(X .LT. 1.0_q) THEN
               FCUT_HELP=0.25_q*(X**3-3.0_q*X+2.0_q)
               DFCUT_HELP=0.75_q*(X**2-1.0_q)*DXDR
            ELSE
               FCUT_HELP=0.0_q
               DFCUT_HELP=0.0_q
            ENDIF
            FMET_HELP=(RMETRIC/(RMETRIC+R(I)))**NMETRIC
            DFMET_HELP=-REAL(NMETRIC,q)/RMETRIC*(RMETRIC/(RMETRIC+R(I)))**(NMETRIC+1)
            DFMET(I)=FCUT_HELP*DFMET_HELP+DFCUT_HELP*FMET_HELP
!            DFMET(I)=-REAL(NMETRIC,q)/RMETRIC*(RMETRIC/(RMETRIC+R(I)))**(NMETRIC+1)
         ENDDO
       END SUBROUTINE PMET1

! Modified spherical Bessel function

       SUBROUTINE IL0 (LMAX,N,M,R,RI,RJ,RIL)
         IMPLICIT NONE
! Input variables
         INTEGER, INTENT(IN)  :: LMAX
         INTEGER, INTENT(IN)  :: M
         INTEGER, INTENT(IN)  :: N
         REAL(q), INTENT(IN)  :: R(:) !(1:M)
         REAL(q), INTENT(IN)  :: RI
         REAL(q), INTENT(IN)  :: RJ(:) !(1:M)
! Output variable
         REAL(q), INTENT(OUT) :: RIL(:,0:) !(1:M,0:LMAX)
! Local variables
         REAL(q) :: EPS
         REAL(q) :: FAC2(0:LMAX)
         INTEGER :: I,L
! Constant
         EPS=1E-10_q
         CALL FACTORIAL2 (FAC2(0:),LMAX+1)
! LMAX=0
         IF(LMAX .EQ. 0) THEN
            DO I=1, N
               IF(R(I)**2 .GT. EPS) THEN
                  RIL(I,0)=0.5_q*(EXP(R(I)-RI-RJ(I))-EXP(-R(I)-RI-RJ(I)))/R(I)
               ELSE
                  RIL(I,0)=EXP(-(RI+RJ(I)))
               ENDIF
            ENDDO
! LMAX=1
         ELSE IF(LMAX .EQ. 1) THEN
            DO I=1, N
               IF(R(I)**2 .GT. EPS) THEN
                  RIL(I,0)=0.5_q*(EXP(R(I)-RI-RJ(I))-EXP(-R(I)-RI-RJ(I)))/R(I)
               ELSE
                  RIL(I,0)=EXP(-(RI+RJ(I)))
               ENDIF
               IF(R(I)**3 .GT. EPS) THEN
                  RIL(I,1)=(0.5_q*R(I)*(EXP(R(I)-RI-RJ(I))+EXP(-R(I)-RI-RJ(I)))-5D-1*(EXP(R(I)-RI-RJ(I))-EXP(-R(I)-RI-RJ(I))))/R(I)**2
               ELSE
                  RIL(I,1)=R(I)/FAC2(1)*EXP(-(RI+RJ(I)))
               ENDIF
            ENDDO
! LMAX > 1
         ELSE
            DO I=1, N
               IF(R(I)**2 .GT. EPS) THEN
                  RIL(I,0)=0.5_q*(EXP(R(I)-RI-RJ(I))-EXP(-R(I)-RI-RJ(I)))/R(I)
               ELSE
                  RIL(I,0)=EXP(-(RI+RJ(I)))
               ENDIF
               IF(R(I)**3 .GT. EPS) THEN
                  RIL(I,1)=(0.5_q*R(I)*(EXP(R(I)-RI-RJ(I))+EXP(-R(I)-RI-RJ(I)))-5D-1*(EXP(R(I)-RI-RJ(I))-EXP(-R(I)-RI-RJ(I))))/R(I)**2
               ELSE
                  RIL(I,1)=R(I)/FAC2(1)*EXP(-(RI+RJ(I)))
               ENDIF
            ENDDO
            DO L=2, LMAX, 1
               DO I=1, N
                  IF(R(I)**(L+2) .GT. EPS) THEN
                     RIL(I,L)=RIL(I,L-2)-(2*(L-1)+1)/R(I)*RIL(I,L-1)
                  ELSE
                     RIL(I,L)=R(I)**L/FAC2(L)*EXP(-(RI+RJ(I)))
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
       END SUBROUTINE IL0

! 1st derivative of modified spherical Bessel function

       SUBROUTINE IL1 (DRIL,LMAX,N,M,R,RIL)
         IMPLICIT NONE
! Input variables
         INTEGER, INTENT(IN) :: LMAX
         INTEGER, INTENT(IN) :: M
         INTEGER, INTENT(IN) :: N
         REAL(q), INTENT(IN) :: R(:) !(1:M)
! Output variables
         REAL(q) DRIL(:,0:) !(1:M,0:LMAX)
         REAL(q) RIL(:,0:) !(1:M,0:LMAX)
! Local variables
         REAL(q) :: EPS
         REAL(q) :: FAC2(0:LMAX)
         INTEGER :: I,L
! Constant
         EPS=1E-10_q
         CALL FACTORIAL2 (FAC2(0:),LMAX+1)
! LMAX=0
         IF(LMAX .EQ. 0) THEN
            DO I=1, N
               IF(R(I)**2 .GT. EPS) THEN
                  RIL(I,0)=DSINH(R(I))/R(I)
                  DRIL(I,0)=(R(I)*DCOSH(R(I))-DSINH(R(I)))/R(I)**2
               ELSE
                  RIL(I,0)=1.0_q
                  DRIL(I,0)=0.0_q
               ENDIF
            ENDDO
! LMAX=1
         ELSE IF (LMAX .EQ. 1) THEN
            DO I=1, N
               IF(R(I)**2 .GT. EPS) THEN
                  RIL(I,0)=DSINH(R(I))/R(I)
                  DRIL(I,0)=(R(I)*DCOSH(R(I))-DSINH(R(I)))/R(I)**2
               ELSE
                  RIL(I,0)=1.0_q
                  DRIL(I,0)=0.0_q
               ENDIF
               IF(R(I)**3 .GT. EPS) THEN
                  RIL(I,1)=(R(I)*DCOSH(R(I))-DSINH(R(I)))/R(I)**2
                  DRIL(I,1)=RIL(I,0)-3.0_q*RIL(I,1)/R(I)
                  DRIL(I,1)=(RIL(I,0)+2.0_q*DRIL(I,1))/3.0_q
               ELSE
                  RIL(I,1)=R(I)/FAC2(1)
                  DRIL(I,1)=1.0_q/FAC2(1)
               ENDIF
            ENDDO
! LMAX > 1
         ELSE
            DO I=1, N
               IF(R(I)**2 .GT. EPS) THEN
                  RIL(I,0)=DSINH(R(I))/R(I)
                  DRIL(I,0)=(R(I)*DCOSH(R(I))-DSINH(R(I)))/R(I)**2
               ELSE
                  RIL(I,0)=1.0_q
                  DRIL(I,0)=0.0_q
               ENDIF
               IF(R(I)**3 .GT. EPS) THEN
                  RIL(I,1)=(R(I)*DCOSH(R(I))-DSINH(R(I)))/R(I)**2
                  DRIL(I,1)=RIL(I,0)-3.0_q*RIL(I,1)/R(I)
                  DRIL(I,1)=(RIL(I,0)+2.0_q*DRIL(I,1))/3.0_q
               ELSE
                  RIL(I,1)=R(I)/FAC2(1)
                  DRIL(I,1)=1.0_q/FAC2(1)
               ENDIF
            ENDDO
            DO L=2, LMAX
               DO I=1, N
                  IF(R(I)**(L+2) .GT. EPS) THEN
                    RIL(I,L)=RIL(I,L-2)-(2*(L-1)+1)/R(I)*RIL(I,L-1)
                  ELSE
                    RIL(I,L)=R(I)**L/FAC2(L)
                  ENDIF
               ENDDO
            ENDDO
            DO L=2, LMAX
               DO I=1, N
                  IF(R(I)**(L+2) .GT. EPS) THEN
                     DRIL(I,L)=RIL(I,L-1)-(2*L+1)/R(I)*RIL(I,L)
                     DRIL(I,L)=(REAL(L,q)*RIL(I,L-1)+REAL((L+1),q)*DRIL(I,L))/REAL((2*L+1),q)
                  ELSE
                     DRIL(I,L)=L*R(I)**(L-1)/FAC2(L)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
       END SUBROUTINE IL1

      END MODULE BASIS_SOAP

!****************************************************************************************************
! Module to handle lattice parameters
!****************************************************************************************************

      MODULE LATTICE_ML
        USE ML_FF_PREC
        USE ML_FF_CONSTANT
        IMPLICIT NONE

        CONTAINS

! Calculate reciprocal latice vectors and cell volume from lattice
! vectors

        SUBROUTINE RECIPROCAL_LATTICE (A,B,OMEGA)
          IMPLICIT NONE
! Input variables
          REAL(q), INTENT(IN) :: A(:,:) !(1:3,1:3)
! Output variables
          REAL(q), INTENT(OUT) :: B(:,:) !(1:3,1:3)
          REAL(q), INTENT(OUT) :: OMEGA
! Local variables
          REAL(q) DPI
! Calculate reciprocal lattice vectors and cell volume
          DPI=2.0_q*PI
          OMEGA=(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
                +A(2,1)*(A(1,3)*A(3,2)-A(1,2)*A(3,3)) &
                +A(3,1)*(A(1,2)*A(2,3)-A(1,3)*A(2,2)))
          B(1,1)=DPI/OMEGA*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
          B(2,1)=DPI/OMEGA*(A(3,2)*A(1,3)-A(3,3)*A(1,2))
          B(3,1)=DPI/OMEGA*(A(1,2)*A(2,3)-A(1,3)*A(2,2))
          B(1,2)=DPI/OMEGA*(A(2,3)*A(3,1)-A(2,1)*A(3,3))
          B(2,2)=DPI/OMEGA*(A(3,3)*A(1,1)-A(3,1)*A(1,3))
          B(3,2)=DPI/OMEGA*(A(1,3)*A(2,1)-A(1,1)*A(2,3))
          B(1,3)=DPI/OMEGA*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
          B(2,3)=DPI/OMEGA*(A(3,1)*A(1,2)-A(3,2)*A(1,1))
          B(3,3)=DPI/OMEGA*(A(1,1)*A(2,2)-A(1,2)*A(2,1))
        END SUBROUTINE RECIPROCAL_LATTICE

      END MODULE LATTICE_ML

!****************************************************************************************************
! Module to make list of neighbours
!****************************************************************************************************

      MODULE NEIGHBOUR
        USE ML_FF_CONSTANT
        USE MPI_DATA
        USE ML_FF_PREC
        IMPLICIT NONE

        CONTAINS

! Set maximum of number of atoms for each element type within cutoff radius.
! MPI-parallel.

       SUBROUTINE SET_SIZE_EST (COMM_WORLD_X,A,B,LADD,LNIONS,MIONS,MITYP,MNEIB_EST,MTYP,NIONS,NITYP,NTYP,RCUT,XYZ)
         USE ML_FF_STRUCT, ONLY : ML_MPI_PAR
         IMPLICIT NONE
         TYPE(ML_MPI_PAR) COMM_WORLD_X
! Input variables
         REAL(q), INTENT(IN) :: A(:,:) !(1:3,1:3)
         REAL(q), INTENT(IN) :: B(:,:) !(1:3,1:3)
         INTEGER, INTENT(IN) :: LADD(:) !(1:MTYP)
         INTEGER, INTENT(IN) :: LNIONS(:,:) !(1:MITYP,1:MTYP)
         INTEGER, INTENT(IN) :: MIONS
         INTEGER, INTENT(IN) :: MITYP
         INTEGER, INTENT(IN) :: MTYP
         INTEGER, INTENT(IN) :: NIONS
         INTEGER, INTENT(IN) :: NITYP(:) !(1:MTYP)
         INTEGER, INTENT(IN) :: NTYP
         REAL(q), INTENT(IN) :: RCUT
         REAL(q), INTENT(IN) :: XYZ(:,:) !(1:3,1:MIONS)
! Output variables
         INTEGER, INTENT(INOUT) :: MNEIB_EST
! Local variables
         INTEGER IERR
         INTEGER INIONS
         INTEGER INIONS_HELP
         INTEGER MNEIB_TMP
         MNEIB_EST=0
         DO INIONS_HELP=1, NIONS, COMM_WORLD_X%NCPU
           INIONS=INIONS_HELP+COMM_WORLD_X%NODE_ME
           IF(INIONS.LE.NIONS) THEN
              CALL SET_SIZE_EST_POINT (A,B,INIONS,LADD,LNIONS,MIONS,MITYP,MNEIB_TMP,MTYP,NIONS,NITYP,NTYP,RCUT,XYZ)
              MNEIB_EST=MAX(MNEIB_EST,MNEIB_TMP)
            ENDIF
         ENDDO
! All-reduce MNEIB_EST
         CALL MPI_ALLREDUCE(MPI_IN_PLACE,MNEIB_EST,1,MPI_INTEGER,MPI_MAX,COMM_WORLD_X%MPI_COMM,IERR)
       END SUBROUTINE SET_SIZE_EST

! Set maximum of number of atoms within cutoff radius from certain atom

       SUBROUTINE SET_SIZE_EST_POINT (A,B,INIONS,LADD,LNIONS,MIONS,MITYP,MNEIB_EST,MTYP,NIONS,NITYP,NTYP,RCUT,XYZ)
         IMPLICIT NONE
! Input variables
         REAL(q), INTENT(IN) :: A(:,:) !(1:3,1:3)
         REAL(q), INTENT(IN) :: B(:,:) !(1:3,1:3)
         INTEGER, INTENT(IN) :: INIONS
         INTEGER, INTENT(IN) :: LADD(:) !(1:MTYP)
         INTEGER, INTENT(IN) :: LNIONS(:,:) !(1:MITYP,1:MTYP)
         INTEGER, INTENT(IN) :: MIONS
         INTEGER, INTENT(IN) :: MITYP
         INTEGER, INTENT(IN) :: MTYP
         INTEGER, INTENT(IN) :: NIONS
         INTEGER, INTENT(IN) :: NITYP(:) !(1:MTYP)
         INTEGER, INTENT(IN) :: NTYP
         REAL(q), INTENT(IN) :: RCUT
         REAL(q), INTENT(IN) :: XYZ(:,:) !(1:3,1:MIONS)
! Output variables
         INTEGER, INTENT(OUT) :: MNEIB_EST
! Local variables
         REAL(q) :: EPS
         INTEGER :: IRCX,IRCY,IRCZ
         INTEGER :: IXYZ
         INTEGER :: IMAX(1:3)
         INTEGER :: IMIN(1:3)
         INTEGER :: JNIONS
         INTEGER :: JNITYP
         INTEGER :: JNTYP
         INTEGER :: JNTYP0
         INTEGER :: JXYZ
         INTEGER :: MNEIB_TMP1(1:MTYP)
         INTEGER :: MNEIB_TMP2
         REAL(q) :: RR(1:3)
         REAL(q) :: RIJR
         REAL(q) :: XYZI(1:3)
         REAL(q) :: XYZJ(1:3)
         REAL(q) :: XYZIJR(1:3)
         REAL(q) :: XYZ_MAX(1:3,1:3)
         REAL(q) :: XYZ_MIN(1:3,1:3)
! Set constant.
         EPS=1E-4_q
! Set variables.
         RR(1:3)=SQRT(B(1,1:3)**2+B(2,1:3)**2+B(3,1:3)**2)
         XYZI(1:3)=XYZ(1:3,INIONS)
! Calculate the edges of spheres.
         DO IXYZ=1, 3
            DO JXYZ=1, 3
               XYZ_MAX(JXYZ,IXYZ)=XYZI(JXYZ)+RCUT*B(JXYZ,IXYZ)/RR(IXYZ)
               XYZ_MIN(JXYZ,IXYZ)=XYZI(JXYZ)-RCUT*B(JXYZ,IXYZ)/RR(IXYZ)
            ENDDO
         ENDDO
! Determin the lattice including the positions of edges.
         DO IXYZ=1, 3
            IMAX(IXYZ)=CEILING(0.5_q*(B(1,IXYZ)*XYZ_MAX(1,IXYZ)+B(2,IXYZ)*XYZ_MAX(2,IXYZ)+B(3,IXYZ)*XYZ_MAX(3,IXYZ))/PI)-1
            IMIN(IXYZ)=CEILING(0.5_q*(B(1,IXYZ)*XYZ_MIN(1,IXYZ)+B(2,IXYZ)*XYZ_MIN(2,IXYZ)+B(3,IXYZ)*XYZ_MIN(3,IXYZ))/PI)-1
         ENDDO
! Then, count the number of neighbours.
         MNEIB_EST=0
         DO JNTYP=1, NTYP
            JNTYP0=LADD(JNTYP)
            MNEIB_TMP2=0
            DO JNITYP=1, NITYP(JNTYP0)
               JNIONS=LNIONS(JNITYP,JNTYP0)
               XYZJ(1:3)=XYZ(1:3,JNIONS)
               DO IRCX=IMIN(1), IMAX(1), 1
                  DO IRCY=IMIN(2), IMAX(2), 1
                     DO IRCZ=IMIN(3), IMAX(3), 1
                        XYZIJR(1:3)=XYZJ(1:3)-XYZI(1:3)+IRCX*A(1:3,1)+IRCY*A(1:3,2)+IRCZ*A(1:3,3)
                        RIJR=SQRT(XYZIJR(1)**2+XYZIJR(2)**2+XYZIJR(3)**2)
                        IF((RIJR .LE. RCUT) .AND. (RIJR .GT. EPS)) then
                           MNEIB_TMP2=MNEIB_TMP2+1
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            MNEIB_EST=MAX(MNEIB_EST,MNEIB_TMP2)
         ENDDO
       END SUBROUTINE SET_SIZE_EST_POINT

! Set structural data on atoms surrounding certain atom

         SUBROUTINE SET_EST_DATA_POINT (INIT,A,B,INIONS,LNTYP,LNIONS_EST,MIONS,MROW_EST,MTYP, &
                                       NIONS,NNEIB_EST,NTYP,RCUT_TARGET,RCUT1,RCUT2,XYZ,XYZ_EST)
           IMPLICIT NONE
! Input variables
           REAL(q), INTENT(IN)  :: A(:,:) !(1:3,1:3)
           REAL(q), INTENT(IN)  :: B(:,:) !(1:3,1:3)
           INTEGER, INTENT(IN)  :: INIONS
           INTEGER, INTENT(IN)  :: INIT
           INTEGER, INTENT(IN)  :: LNTYP(:) !(1:MIONS)
           INTEGER, INTENT(IN)  :: MIONS
           INTEGER, INTENT(IN)  :: MROW_EST
           INTEGER, INTENT(IN)  :: MTYP
           INTEGER, INTENT(IN)  :: NIONS
           INTEGER, INTENT(IN)  :: NTYP
           REAL(q), INTENT(IN)  :: RCUT_TARGET
           REAL(q), INTENT(IN)  :: RCUT1
           REAL(q), INTENT(IN)  :: RCUT2
           REAL(q), INTENT(IN)  :: XYZ(:,:) !(1:3,1:MIONS)
! Output variables
           INTEGER, INTENT(OUT) :: LNIONS_EST(:,:) !(1:MROW_EST,1:MTYP)
           INTEGER, INTENT(OUT) :: NNEIB_EST(:) !(1:MTYP)
           REAL(q), INTENT(OUT) :: XYZ_EST(:,:,:) !(1:3,1:MROW_EST,1:MTYP)
! Local variables
           REAL(q)              :: EPS
           INTEGER              :: IMAX(1:3)
           INTEGER              :: IMIN(1:3)
           INTEGER              :: INEIB
           INTEGER              :: INTYP0
           INTEGER              :: IRCX
           INTEGER              :: IRCY
           INTEGER              :: IRCZ
           INTEGER              :: IXYZ
           INTEGER              :: JNIONS
           INTEGER              :: JNTYP0
           INTEGER              :: JXYZ
           INTEGER              :: LNIONS_EST_HELP(1:MROW_EST-1,1:MTYP)
           INTEGER              :: NNEIB_EST_HELP(1:MTYP)
           REAL(q)              :: RCUT_MAX
           REAL(q)              :: RCUT_MIN
           REAL(q)              :: RR(1:3)
           REAL(q)              :: RIJR
           REAL(q)              :: XYZI(1:3)
           REAL(q)              :: XYZJ(1:3)
           REAL(q)              :: XYZIJR(1:3)
           REAL(q)              :: XYZ_EST_HELP(1:3,1:MROW_EST-1,1:MTYP)
           REAL(q)              :: XYZ_MAX(1:3,1:3)
           REAL(q)              :: XYZ_MIN(1:3,1:3)
! Set maxmum and minimum cutoff radii.
           RCUT_MAX=MAX(RCUT1,RCUT2)
           RCUT_MIN=MIN(RCUT1,RCUT2)
! Set constant
           EPS=1E-4_q
! Initializatiion of data if necessary.
           IF(INIT.EQ.1) THEN
              LNIONS_EST=0
              LNIONS_EST_HELP=0
              NNEIB_EST=0
              NNEIB_EST_HELP=0
              XYZ_EST=0.0_q
              XYZ_EST_HELP=0.0_q
           ENDIF
! Set variables.
           INTYP0=LNTYP(INIONS)
           XYZI(1:3)=XYZ(1:3,INIONS)
           RR(1:3)=SQRT(B(1,1:3)**2+B(2,1:3)**2+B(3,1:3)**2)
           XYZI(1:3)=XYZ(1:3,INIONS)
! Calculate the edges of spheres.
           DO IXYZ=1, 3
              DO JXYZ=1, 3
                 XYZ_MAX(JXYZ,IXYZ)=XYZI(JXYZ)+RCUT_TARGET*B(JXYZ,IXYZ)/RR(IXYZ)
                 XYZ_MIN(JXYZ,IXYZ)=XYZI(JXYZ)-RCUT_TARGET*B(JXYZ,IXYZ)/RR(IXYZ)
              ENDDO
           ENDDO
! Determin the lattice including the positions of edges.
            DO IXYZ=1, 3
               IMAX(IXYZ)=CEILING(0.5_q*(B(1,IXYZ)*XYZ_MAX(1,IXYZ)+B(2,IXYZ)*XYZ_MAX(2,IXYZ)+B(3,IXYZ)*XYZ_MAX(3,IXYZ))/PI)-1
               IMIN(IXYZ)=CEILING(0.5_q*(B(1,IXYZ)*XYZ_MIN(1,IXYZ)+B(2,IXYZ)*XYZ_MIN(2,IXYZ)+B(3,IXYZ)*XYZ_MIN(3,IXYZ))/PI)-1
            ENDDO
! Then, make the configurational data.
           LNIONS_EST(1,INTYP0)=INIONS
           XYZ_EST(1:3,1,INTYP0)=0.0_q
           DO JNIONS=1, NIONS
              JNTYP0=LNTYP(JNIONS)
              XYZJ(1:3)=XYZ(1:3,JNIONS)
              DO IRCX=IMIN(1), IMAX(1), 1
                 DO IRCY=IMIN(2), IMAX(2), 1
                    DO IRCZ=IMIN(3), IMAX(3), 1
                       XYZIJR(1:3)=XYZJ(1:3)-XYZI(1:3)+IRCX*A(1:3,1)+IRCY*A(1:3,2)+IRCZ*A(1:3,3)
                       RIJR=SQRT(XYZIJR(1)**2+XYZIJR(2)**2+XYZIJR(3)**2)
                       IF((RIJR .LE. RCUT_MIN) .AND.  (RIJR .GT. EPS)) THEN
                          NNEIB_EST(JNTYP0)=NNEIB_EST(JNTYP0)+1
                          LNIONS_EST(NNEIB_EST(JNTYP0)+1,JNTYP0)=JNIONS
                          XYZ_EST(1:3,NNEIB_EST(JNTYP0)+1,JNTYP0)=XYZIJR(1:3)
                       ELSE IF((RIJR .LE. RCUT_TARGET) .AND.  (RIJR .GT. RCUT_MIN)) THEN
                          NNEIB_EST_HELP(JNTYP0)=NNEIB_EST_HELP(JNTYP0)+1
                          LNIONS_EST_HELP(NNEIB_EST_HELP(JNTYP0),JNTYP0)=JNIONS
                          XYZ_EST_HELP(1:3,NNEIB_EST_HELP(JNTYP0),JNTYP0)=XYZIJR(1:3)
                       ENDIF
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           IF(RCUT_TARGET.GT.RCUT_MIN) THEN
              DO JNTYP0=1, NTYP
                 DO INEIB=1, NNEIB_EST_HELP(JNTYP0)
                    NNEIB_EST(JNTYP0)=NNEIB_EST(JNTYP0)+1
                    LNIONS_EST(NNEIB_EST(JNTYP0)+1,JNTYP0)=LNIONS_EST_HELP(INEIB,JNTYP0)
                    XYZ_EST(1:3,NNEIB_EST(JNTYP0)+1,JNTYP0)=XYZ_EST_HELP(1:3,INEIB,JNTYP0)
                 ENDDO
              ENDDO
           ENDIF
         END SUBROUTINE SET_EST_DATA_POINT

      END MODULE NEIGHBOUR

!****************************************************************************************************
! Module for calculating angular descriptor
!****************************************************************************************************

      MODULE SOAP_KERNEL
        USE BASIS_SOAP
        USE ML_FF_CONSTANT
        USE MPI_DATA
        USE ML_FF_PREC
        IMPLICIT NONE

        CONTAINS

! Subroutine to calculate expansion coefficients

        SUBROUTINE D0CLM_SOAP (CLM,CLM_COUPLE,DRSPL,FNL,FWIN, &
                   INTYP0,LADD,LCOUPLE,LFLAG_VAR,LFLAG_VAR_SIC,LIONS_COUPLE,LMAX,LNIONS_EST,LSIC,LVAR,LVAR_SIC,LWINDOW, &
                   MIONS,MROW_EST,MRB,MSPL,MTYP,MROW_CLM,MROW_CLM_COUPLE,MROW_PS_SIC,MROW_PS_SIC_COUPLE, &
                   NNEIB_EST,NRB,NSPL,NTYP,PS_SIC,PS_SIC_COUPLE,RCOUPLE_IONS,RMEM,XYZ_EST)
          USE ML_SPLINE
          IMPLICIT NONE
! Input variables
          REAL(q), INTENT(IN)  :: DRSPL
          REAL(q), INTENT(IN)  :: FNL(:,:,0:,:,:) !(1:MSPL,1:5,0:LMAX,1:MRB,1:MTYP)
          REAL(q), INTENT(IN)  :: FWIN(0:,:) !(0:LMAX,1:MRB)
          INTEGER, INTENT(IN)  :: INTYP0
          INTEGER, INTENT(IN)  :: LADD(:) !(1:MTYP)
          LOGICAL, INTENT(IN)  :: LCOUPLE
          LOGICAL, INTENT(IN)  :: LFLAG_VAR(:,:,0:,:,:,:) !(1:MRB,1:MRB,0:LMAX,1:MTYP,1:MTYP,1:MTYP)
          LOGICAL, INTENT(IN)  :: LFLAG_VAR_SIC(:,:,0:,:,:) !(1:MRB,1:MRB,0:LMAX,1:MTYP,1:MTYP)
          LOGICAL, INTENT(IN)  :: LIONS_COUPLE(:) !(1:MIONS)
          INTEGER, INTENT(IN)  :: LMAX
          INTEGER, INTENT(IN)  :: LNIONS_EST(:,:) !(1:MROW_EST,1:MTYP)
          LOGICAL, INTENT(IN)  :: LSIC
          INTEGER, INTENT(IN)  :: LVAR(:,:,0:,:,:,:) !(1:MRB,1:MRB,0:LMAX,1:MTYP,1:MTYP,1:MTYP)
          INTEGER, INTENT(IN)  :: LVAR_SIC(:,:,0:,:,:) !(1:MRB,1:MRB,0:LMAX,1:MTYP,1:MTYP)
          LOGICAL, INTENT(IN)  :: LWINDOW
          INTEGER, INTENT(IN)  :: MIONS
          INTEGER, INTENT(IN)  :: MRB
          INTEGER, INTENT(IN)  :: MROW_CLM
          INTEGER, INTENT(IN)  :: MROW_CLM_COUPLE
          INTEGER, INTENT(IN)  :: MROW_EST
          INTEGER, INTENT(IN)  :: MROW_PS_SIC
          INTEGER, INTENT(IN)  :: MROW_PS_SIC_COUPLE
          INTEGER, INTENT(IN)  :: MSPL
          INTEGER, INTENT(IN)  :: MTYP
          INTEGER, INTENT(IN)  :: NNEIB_EST(:) !(1:MTYP)
          INTEGER, INTENT(IN)  :: NRB(0:) !(0:LMAX)
          INTEGER, INTENT(IN)  :: NSPL
          INTEGER, INTENT(IN)  :: NTYP
          REAL(q), INTENT(IN)  :: RCOUPLE_IONS(:) !(1:MIONS)
          REAL(q), INTENT(IN)  :: XYZ_EST(:,:,:) !(1:3,1:MROW_EST,1:MTYP)
! Output variables
          REAL(q), INTENT(OUT) :: CLM(:) !(1:MROW_CLM)
          REAL(q), INTENT(OUT) :: CLM_COUPLE(:) !(1:MROW_CLM_COUPLE)
          REAL(q), INTENT(OUT) :: PS_SIC(:) !(1:MROW_PS_SIC)
          REAL(q), INTENT(OUT) :: PS_SIC_COUPLE(:) !(1:MROW_PS_SIC_COUPLE)
          REAL(q), INTENT(OUT) :: RMEM
! Local variables
          REAL(q)              :: CLM_LOC(1:MRB)
          REAL(q)              :: F(1:MRB,0:LMAX)
          REAL(q)              :: FDER
          REAL(q)              :: FWIN_LOC
          INTEGER              :: IB
          INTEGER              :: ILM
          INTEGER              :: INEIB
          INTEGER              :: INEIB_TMP
          INTEGER              :: IVAR
          INTEGER              :: JB
          INTEGER              :: JNIONS
          INTEGER              :: JNTYP
          INTEGER              :: JNTYP0
          INTEGER              :: L
          INTEGER              :: M
          INTEGER              :: NAC
          REAL(q)              :: R
          REAL(q)              :: YLM(1:(LMAX+1)**2) ! Spherical function
! Record the memory usage for local variables
          RMEM=REAL((LMAX+1)*(MRB+LMAX+1),q)*8/1E6_q ! YLM and F
! Calculations of expansion coefficients for each element distribution
! Calculate invidual coefficient for 
          DO JNTYP=1, NTYP
             JNTYP0=LADD(JNTYP)
             DO INEIB=2, NNEIB_EST(JNTYP0)+1
                JNIONS=LNIONS_EST(INEIB,JNTYP0)
! Compute necessary spherical harmonics.
                CALL YLM0 (LMAX+1,(LMAX+1)**2,R,XYZ_EST(1:,INEIB,JNTYP0),YLM)
! Compute radial functions.
                NAC=INT(R/DRSPL)+1
                DO L=0, LMAX
!NEC$ novector
                   DO IB=1, NRB(L)
                      IF(LWINDOW) THEN
                         FWIN_LOC=FWIN(L,IB)
                      ELSE
                         FWIN_LOC=1.0_q
                      ENDIF
                      CALL SPLVAL_ML (R,F(IB,L),FDER,FNL(1:,1:,L,IB,JNTYP0),NAC,NSPL)
                      F(IB,L)=FWIN_LOC*F(IB,L)
                   ENDDO
                ENDDO
                DO L=0, LMAX
                   DO M=1, 2*L+1
                      ILM=L**2+M
! Compute the expansion coefficient.
                      DO IB=1, NRB(L)
                         CLM_LOC(IB)=F(IB,L)*YLM(ILM)
                      ENDDO
                      DO IB=1, NRB(L)
! Sum up the expansion coefficients.
                         IVAR=(LMAX+1)**2*MTYP*(IB-1)+(LMAX+1)**2*(JNTYP0-1)+ILM
                         CLM(IVAR)=CLM(IVAR)+RCOUPLE_IONS(JNIONS)*CLM_LOC(IB)
                         IF(LCOUPLE) THEN
                            IF(LIONS_COUPLE(JNIONS)) THEN
                               CLM_COUPLE(IVAR)=CLM_COUPLE(IVAR)+CLM_LOC(IB)
                            ENDIF
                         ENDIF
! Compute the self-interaction terms if necessary.
                         IF(LSIC) THEN
                            DO JB=IB, NRB(L)
                               IF(LFLAG_VAR_SIC(JB,IB,L,JNTYP0,INTYP0)) THEN
                                  IVAR=LVAR_SIC(JB,IB,L,JNTYP0,INTYP0)
                                  PS_SIC(IVAR)=PS_SIC(IVAR)+RCOUPLE_IONS(JNIONS)**2*CLM_LOC(IB)*CLM_LOC(JB)
! Compute the self-interaction coupling term if necessary.
                                  IF(LCOUPLE) THEN
                                     IF(LIONS_COUPLE(JNIONS)) THEN
                                        PS_SIC_COUPLE(IVAR)=PS_SIC_COUPLE(IVAR)+2.0_q*RCOUPLE_IONS(JNIONS)*CLM_LOC(IB)*CLM_LOC(JB)
                                     ENDIF
                                  ENDIF
                               ENDIF
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
        END SUBROUTINE D0CLM_SOAP

! Store l=m=0 component if necessary.
! This subroutine is used for calculations without using FFM%PS_ALL.
! Parallelized by scaLAPACK

        SUBROUTINE D0C00_SOAP (C00,C00_COUPLE,CLM,CLM_COUPLE,LADD,LCOUPLE,LMAX,LNORM1, &
                   MROW_C00,MROW_C00_COUPLE,MRB,MTYP,MROW_CLM,MROW_CLM_COUPLE,NRB,NTYP)
          IMPLICIT NONE
! Input variables
          REAL(q), INTENT(IN)  :: CLM(:) !(MROW_CLM)
          REAL(q), INTENT(IN)  :: CLM_COUPLE(:) !(MROW_CLM_COUPLE)
          INTEGER, INTENT(IN)  :: LADD(:) !(1:MTYP)
          LOGICAL, INTENT(IN)  :: LCOUPLE
          INTEGER, INTENT(IN)  :: LMAX
          LOGICAL, INTENT(IN)  :: LNORM1
          INTEGER, INTENT(IN)  :: MROW_C00
          INTEGER, INTENT(IN)  :: MROW_C00_COUPLE
          INTEGER, INTENT(IN)  :: MROW_CLM
          INTEGER, INTENT(IN)  :: MROW_CLM_COUPLE
          INTEGER, INTENT(IN)  :: MRB
          INTEGER, INTENT(IN)  :: MTYP
          INTEGER, INTENT(IN)  :: NRB(0:) !(0:LMAX)
          INTEGER, INTENT(IN)  :: NTYP
! Output variables
          REAL(q), INTENT(OUT) :: C00(:) !(1:MROW_C00)
          REAL(q), INTENT(OUT) :: C00_COUPLE(:) !(1:MROW_C00_COUPLE)
! Local variables
          REAL(q)              :: C00_COUPLE_SUM
          INTEGER              :: IB
          INTEGER              :: IVAR1
          INTEGER              :: IVAR2
          INTEGER              :: JNTYP
          INTEGER              :: JNTYP0
          INTEGER              :: JJNTYP0
          REAL(q)              :: RNORM
! LAPACK routine.
          REAL(q) :: DNRM2
! Initialization.
          C00=0.0_q
          IF(LCOUPLE) THEN
             C00_COUPLE=0.0_q
          ENDIF
! Calculate C00 and its norm.
          DO JNTYP=1, NTYP
             JNTYP0=LADD(JNTYP)
             DO IB=1, NRB(0)
                IVAR1=MRB*(JNTYP0-1)+IB
                IVAR2=(LMAX+1)**2*MTYP*(IB-1)+(LMAX+1)**2*(JNTYP0-1)+1
                C00(IVAR1)=CLM(IVAR2)
                IF(LCOUPLE) THEN
                   C00_COUPLE(IVAR1)=CLM_COUPLE(IVAR2)
                ENDIF
             ENDDO
          ENDDO
        END SUBROUTINE D0C00_SOAP

! Subroutine calculating angular descriptor part of SOAP

        SUBROUTINE D0PS_SOAP (AFILT,CLM,CLM_COUPLE,IAFILT,INTYP0, &
                   LAFILT,LADD,LCOUPLE,LFLAG_VAR,LFLAG_VAR_SIC,L_LNRB,LMAX,LNORM2,LSIC,LVAR,LVAR_SIC, &
                   MLNRB,MROW_PS,MROW_PS_COUPLE,MROW_PS_SIC,MROW_PS_SIC_COUPLE,MRB,MTYP,MROW_CLM,MROW_CLM_COUPLE, &
                   NLNRB,NRB,NRB_LNRB,NTYP,PS,PS_COUPLE,PS_SIC,PS_SIC_COUPLE,WVAR)
          IMPLICIT NONE
! Input variables
          REAL(q), INTENT(IN)  :: AFILT
          REAL(q), INTENT(IN)  :: CLM(:) !(1:MROW_CLM)
          REAL(q), INTENT(IN)  :: CLM_COUPLE(:) !(1:MROW_CLM_COUPLE)
          INTEGER, INTENT(IN)  :: IAFILT
          INTEGER, INTENT(IN)  :: INTYP0
          LOGICAL, INTENT(IN)  :: LAFILT
          INTEGER, INTENT(IN)  :: LADD(:) !(1:MTYP)
          LOGICAL, INTENT(IN)  :: LCOUPLE
          LOGICAL, INTENT(IN)  :: LFLAG_VAR(:,:,0:,:,:,:) !(1:MRB,1:MRB,0:LMAX,1:MTYP,1:MTYP,1:MTYP)
          LOGICAL, INTENT(IN)  :: LFLAG_VAR_SIC(:,:,0:,:,:) !(1:MRB,1:MRB,0:LMAX,1:MTYP,1:MTYP)
          INTEGER, INTENT(IN)  :: L_LNRB(:) !(1:MLNRB)
          INTEGER, INTENT(IN)  :: LMAX
          LOGICAL, INTENT(IN)  :: LNORM2
          LOGICAL, INTENT(IN)  :: LSIC
          INTEGER, INTENT(IN)  :: LVAR(:,:,0:,:,:,:) !(1:MRB,1:MRB,0:LMAX,1:MTYP,1:MTYP,1:MTYP)
          INTEGER, INTENT(IN)  :: LVAR_SIC(:,:,0:,:,:) !(1:MRB,1:MRB,0:LMAX,1:MTYP,1:MTYP)
          INTEGER, INTENT(IN)  :: MLNRB
          INTEGER, INTENT(IN)  :: MROW_CLM
          INTEGER, INTENT(IN)  :: MROW_CLM_COUPLE
          INTEGER, INTENT(IN)  :: MROW_PS
          INTEGER, INTENT(IN)  :: MROW_PS_COUPLE
          INTEGER, INTENT(IN)  :: MROW_PS_SIC
          INTEGER, INTENT(IN)  :: MROW_PS_SIC_COUPLE
          INTEGER, INTENT(IN)  :: MRB
          INTEGER, INTENT(IN)  :: MTYP
          INTEGER, INTENT(IN)  :: NLNRB
          INTEGER, INTENT(IN)  :: NRB(0:) !(0:LMAX)
          INTEGER, INTENT(IN)  :: NRB_LNRB(:) !(1:MLNRB)
          INTEGER, INTENT(IN)  :: NTYP
          REAL(q), INTENT(IN)  :: PS_SIC(:) !(1:MROW_PS_SIC)
          REAL(q), INTENT(IN)  :: PS_SIC_COUPLE(:) !(1:MROW_PS_SIC_COUPLE)
          REAL(q), INTENT(IN)  :: WVAR(:,:,0:,:,:,:) !(1:MRB,1:MRB,0:LMAX,1:MTYP,1:MTYP,1:MTYP)
! Output variables
          REAL(q), INTENT(OUT) :: PS(:) !(1:MROW_PS)
          REAL(q), INTENT(OUT) :: PS_COUPLE(:) !(1:MROW_PS_COUPLE)
! Local variables
          INTEGER              :: IB
          INTEGER              :: IERR
          INTEGER              :: ILM
          INTEGER              :: ILNR
          INTEGER              :: ILNRB
          INTEGER              :: INEIB
          INTEGER              :: IVAR1
          INTEGER              :: IVAR2
          INTEGER              :: IVAR3
          INTEGER              :: IVAR4
          INTEGER              :: JB
          INTEGER              :: JNTYP
          INTEGER              :: JNTYP0
          INTEGER              :: JJNTYP
          INTEGER              :: JJNTYP0
          INTEGER              :: KNTYP0
          INTEGER              :: L
          INTEGER              :: LOCVAR_ROW
          INTEGER              :: M
          REAL(q)              :: PI8
          REAL(q)              :: PREFAC(0:LMAX)
          REAL(q)              :: PS_LOC
          REAL(q)              :: PS_COUPLE_LOC
          REAL(q)              :: PS_COUPLE_SUM
          REAL(q)              :: RNORM
! Lapack routine.
          REAL(q) :: DNRM2
! Set constant variables
          PI8=8.0_q*PI**2
          DO L=0, LMAX
             IF(LAFILT) THEN
                IF(IAFILT.EQ.1) THEN
                   PREFAC(L)=SQRT(PI8/REAL((2*L+1),q))/SQRT(REAL((2*L+1),q))
                ELSE IF(IAFILT.EQ.2) THEN
                   PREFAC(L)=SQRT(PI8/REAL((2*L+1),q))/(1.0_q+AFILT*REAL((L*(L+1)),q)**2)**2
                ENDIF
             ELSE
                PREFAC(L)=SQRT(PI8/REAL((2*L+1),q))
             ENDIF
          ENDDO
! Calculate angular descriptor.
          PS=0.0_q
          IF(LCOUPLE) THEN
             PS_COUPLE=0.0_q
          ENDIF
          DO JNTYP=1, NTYP
             JNTYP0=LADD(JNTYP)
             DO JJNTYP=1, NTYP
                JJNTYP0=LADD(JJNTYP)
                DO ILNRB=1, NLNRB
                   L=L_LNRB(ILNRB)
                   IB=NRB_LNRB(ILNRB)
                   DO JB=IB, NRB(L)
                      IF(LFLAG_VAR(JB,IB,L,JJNTYP0,JNTYP0,INTYP0)) THEN
                         IVAR3=LVAR(JB,IB,L,JJNTYP0,JNTYP0,INTYP0)
                         PS_LOC=0.0_q
                         IF(LCOUPLE) THEN
                            PS_COUPLE_LOC=0.0_q
                         ENDIF
! Calculate the angular descriptor that contains self-interaction terms.
                         DO M=1, 2*L+1
                            ILM=L**2+M
                            IVAR1=(LMAX+1)**2*MTYP*(IB-1)+(LMAX+1)**2*( JNTYP0-1)+ILM
                            IVAR2=(LMAX+1)**2*MTYP*(JB-1)+(LMAX+1)**2*(JJNTYP0-1)+ILM
                            PS_LOC=PS_LOC+CLM(IVAR1)*CLM(IVAR2)
                            IF(LCOUPLE) THEN
                               PS_COUPLE_LOC=PS_COUPLE_LOC+CLM_COUPLE(IVAR1)*CLM(IVAR2)+CLM_COUPLE(IVAR2)*CLM(IVAR1)
                            ENDIF
                         ENDDO
! Self-interaction correction.
                         IF(LFLAG_VAR_SIC(JB,IB,L,JJNTYP0,INTYP0).AND.(JJNTYP0.EQ.JNTYP0).AND.LSIC) THEN
                            IVAR4=LVAR_SIC(JB,IB,L,JJNTYP0,INTYP0)
                            PS_LOC=PS_LOC-PS_SIC(IVAR4)
                            IF(LCOUPLE) THEN
                               PS_COUPLE_LOC=PS_COUPLE_LOC-PS_SIC_COUPLE(IVAR4)
                            ENDIF
                         ENDIF
! Calculate the angular descriptor.
                         PS(IVAR3)=WVAR(JB,IB,L,JJNTYP0,JNTYP0,INTYP0)*PREFAC(L)*PS_LOC
                         IF(LCOUPLE) THEN
                            PS_COUPLE(IVAR3)=WVAR(JB,IB,L,JJNTYP0,JNTYP0,INTYP0)*PREFAC(L)*PS_COUPLE_LOC
                         ENDIF
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
        END SUBROUTINE D0PS_SOAP

! Subroutine to normalize the descriptor

        SUBROUTINE NORMALIZATION_D0 (C00,C00_COUPLE,LSUPERVEC,LCOUPLE,LNORM1,LNORM2, &
                   MROW_C00,MROW_C00_COUPLE,MROW_PS,MROW_PS_COUPLE,NNVAR1,NNVAR2,PS, &
                   PS_COUPLE,RNORM_ALL,RNORM_C00,RNORM_PS,W1,W2)
          IMPLICIT NONE
! Input variables.
          LOGICAL, INTENT(IN)    :: LSUPERVEC
          LOGICAL, INTENT(IN)    :: LCOUPLE
          LOGICAL, INTENT(IN)    :: LNORM1
          LOGICAL, INTENT(IN)    :: LNORM2
          INTEGER, INTENT(IN)    :: MROW_C00
          INTEGER, INTENT(IN)    :: MROW_C00_COUPLE
          INTEGER, INTENT(IN)    :: MROW_PS
          INTEGER, INTENT(IN)    :: MROW_PS_COUPLE
          INTEGER, INTENT(IN)    :: NNVAR1
          INTEGER, INTENT(IN)    :: NNVAR2
          REAL(q), INTENT(IN)    :: W1
          REAL(q), INTENT(IN)    :: W2
! In and output variables.
          REAL(q), INTENT(INOUT) :: C00(:) !(1:MROW_C00)
          REAL(q), INTENT(INOUT) :: C00_COUPLE(:) !(1:MROW_C00_COUPLE)
          REAL(q), INTENT(INOUT) :: PS(:) !(1:MROW_PS)
          REAL(q), INTENT(INOUT) :: PS_COUPLE(:) !(1:MROW_PS_COUPLE)
! Output variables.
          REAL(q), INTENT(OUT)   :: RNORM_ALL
          REAL(q), INTENT(OUT)   :: RNORM_C00
          REAL(q), INTENT(OUT)   :: RNORM_PS
! Local variables.
          REAL(q)                :: ALL_COUPLE_SUM
          REAL(q)                :: C00_COUPLE_SUM
          INTEGER                :: IVAR
          REAL(q)                :: PS_COUPLE_SUM
          REAL(q)                :: RNORM
! Lapack routines.
          REAL(q) :: DNRM2
          IF(.NOT.LSUPERVEC) THEN
             IF(W1.GT.0.0_q) THEN
                RNORM_C00=DNRM2(NNVAR1,C00(1),1)
                IF(LNORM1) THEN
                   RNORM=1.0_q/RNORM_C00
                   CALL DSCAL (NNVAR1,RNORM,C00(1),1)
                   IF(LCOUPLE) THEN
                      CALL DSCAL (NNVAR1,RNORM,C00_COUPLE(1),1)
                      C00_COUPLE_SUM=0.0_q
                      DO IVAR=1, NNVAR1
                         C00_COUPLE_SUM=C00_COUPLE_SUM+C00_COUPLE(IVAR)*C00(IVAR)
                      ENDDO
                      DO IVAR=1, NNVAR1
                         C00_COUPLE(IVAR)=C00_COUPLE(IVAR)-C00(IVAR)*C00_COUPLE_SUM
                      ENDDO
                   ENDIF
                ENDIF
                RNORM_C00=RNORM_C00**2
             ENDIF
             IF(W2.GT.0.0_q) THEN
                RNORM_PS=DNRM2(NNVAR2,PS(1),1)
                IF(LNORM2) THEN
                   RNORM=1.0_q/RNORM_PS
                   CALL DSCAL (NNVAR2,RNORM,PS(1),1)
                   IF(LCOUPLE) THEN
                      CALL DSCAL (NNVAR2,RNORM,PS_COUPLE(1),1)
                      PS_COUPLE_SUM=0.0_q
                      DO IVAR=1, NNVAR2
                         PS_COUPLE_SUM=PS_COUPLE_SUM+PS_COUPLE(IVAR)*PS(IVAR)
                      ENDDO
                      DO IVAR=1, NNVAR2
                         PS_COUPLE(IVAR)=PS_COUPLE(IVAR)-PS(IVAR)*PS_COUPLE_SUM
                      ENDDO
                   ENDIF
                ENDIF
                RNORM_PS=RNORM_PS**2
             ENDIF
          ELSE
             RNORM_ALL=0.0_q
             IF(W1.GT.0.0_q) THEN
                RNORM=DNRM2(NNVAR1,C00(1),1)
                RNORM_ALL=RNORM_ALL+W1*RNORM**2
             ENDIF
             IF(W2.GT.0.0_q) THEN
                RNORM=DNRM2(NNVAR2,PS(1),1)
                RNORM_ALL=RNORM_ALL+W2*RNORM**2
             ENDIF
             RNORM_ALL=SQRT(RNORM_ALL)
             ALL_COUPLE_SUM=0.0_q
             IF(W1.GT.0.0_q) THEN
                IF(LNORM1) THEN
                   RNORM=SQRT(W1)/RNORM_ALL
                   CALL DSCAL (NNVAR1,RNORM,C00(1),1)
                   IF(LCOUPLE) THEN
                      CALL DSCAL (NNVAR1,RNORM,C00_COUPLE(1),1)
                      C00_COUPLE_SUM=0.0_q
                      DO IVAR=1, NNVAR1
                         C00_COUPLE_SUM=C00_COUPLE_SUM+C00_COUPLE(IVAR)*C00(IVAR)
                      ENDDO
                      ALL_COUPLE_SUM=ALL_COUPLE_SUM+C00_COUPLE_SUM
                   ENDIF
                ENDIF
             ENDIF
             IF(W2.GT.0.0_q) THEN
                IF(LNORM2) THEN
                   RNORM=SQRT(W2)/RNORM_ALL
                   CALL DSCAL (NNVAR2,RNORM,PS(1),1)
                   IF(LCOUPLE) THEN
                      CALL DSCAL (NNVAR2,RNORM,PS_COUPLE(1),1)
                      PS_COUPLE_SUM=0.0_q
                      DO IVAR=1, NNVAR2
                         PS_COUPLE_SUM=PS_COUPLE_SUM+PS_COUPLE(IVAR)*PS(IVAR)
                      ENDDO
                      ALL_COUPLE_SUM=ALL_COUPLE_SUM+PS_COUPLE_SUM
                   ENDIF
                ENDIF
             ENDIF
             IF(LCOUPLE) THEN
                IF(W1.GT.0.0_q) THEN
                   IF(LNORM1) THEN
                      DO IVAR=1, NNVAR1
                         C00_COUPLE(IVAR)=C00_COUPLE(IVAR)-C00(IVAR)*ALL_COUPLE_SUM
                      ENDDO
                   ENDIF
                ENDIF
                IF(W2.GT.0.0_q) THEN
                   IF(LNORM2) THEN
                      DO IVAR=1, NNVAR2
                         PS_COUPLE(IVAR)=PS_COUPLE(IVAR)-PS(IVAR)*ALL_COUPLE_SUM
                      ENDDO
                   ENDIF
                ENDIF
             ENDIF
             RNORM_ALL=RNORM_ALL**2
          ENDIF
        END SUBROUTINE NORMALIZATION_D0

! Calculate derivatives of expansion coefficients
! This subroutine is used in calculations using FFM%C00_ALL and FFM%PS_ALL
! Inside of this subroutine is not parallelized.

        SUBROUTINE D1CLM_SOAP (DCLM,DRSPL,FNL,FWIN,KNIONS,KNTYP0, &
                   LADD,LFLAG_VAR,LFLAG_VAR_SIC,LMAX,LNIONS_EST,LSIC,LVAR,LVAR_SIC,LWINDOW, &
                   MCOL_DCLM,MCOL_DPS_SIC,MIONS,MRB, &
                   MROW_DCLM,MROW_DPS_SIC,MROW_EST,MSPL,MTYP, &
                   NNEIB_EST,NRB,NSPL,NTYP, &
                   DPS_SIC,RCOUPLE_IONS,RMEM,XYZ_EST)
          USE ML_SPLINE
          IMPLICIT NONE
! Input variables
          REAL(q), INTENT(IN)  :: DRSPL
          REAL(q), INTENT(IN)  :: FNL(:,:,0:,:,:) !(1:MSPL,1:5,0:LMAX,1:MRB,1:MTYP)
          REAL(q), INTENT(IN)  :: FWIN(0:,:) !(0:LMAX,1:MRB)
          INTEGER, INTENT(IN)  :: KNIONS
          INTEGER, INTENT(IN)  :: KNTYP0
          INTEGER, INTENT(IN)  :: LADD(:) !(1:MTYP)
          LOGICAL, INTENT(IN)  :: LFLAG_VAR(:,:,0:,:,:,:) !(1:MRB,1:MRB,0:LMAX,1:MTYP,1:MTYP,1:MTYP)
          LOGICAL, INTENT(IN)  :: LFLAG_VAR_SIC(:,:,0:,:,:) !(1:MRB,1:MRB,0:LMAX,1:MTYP,1:MTYP)
          INTEGER, INTENT(IN)  :: LMAX
          INTEGER, INTENT(IN)  :: LNIONS_EST(:,:) !(1:MROW_EST,1:MTYP)
          LOGICAL, INTENT(IN)  :: LSIC
          INTEGER, INTENT(IN)  :: LVAR(:,:,0:,:,:,:) !(1:MRB,1:MRB,0:LMAX,1:MTYP,1:MTYP,1:MTYP)
          INTEGER, INTENT(IN)  :: LVAR_SIC(:,:,0:,:,:) !(1:MRB,1:MRB,0:LMAX,1:MTYP,1:MTYP)
          LOGICAL, INTENT(IN)  :: LWINDOW
          INTEGER, INTENT(IN)  :: MCOL_DCLM
          INTEGER, INTENT(IN)  :: MCOL_DPS_SIC
          INTEGER, INTENT(IN)  :: MIONS
          INTEGER, INTENT(IN)  :: MROW_DCLM
          INTEGER, INTENT(IN)  :: MROW_DPS_SIC
          INTEGER, INTENT(IN)  :: MROW_EST
          INTEGER, INTENT(IN)  :: MRB
          INTEGER, INTENT(IN)  :: MSPL
          INTEGER, INTENT(IN)  :: MTYP
          INTEGER, INTENT(IN)  :: NNEIB_EST(:) !(1:MTYP)
          INTEGER, INTENT(IN)  :: NRB(0:) !(0:LMAX)
          INTEGER, INTENT(IN)  :: NSPL
          INTEGER, INTENT(IN)  :: NTYP
          REAL(q), INTENT(IN)  :: RCOUPLE_IONS(:) !(1:MIONS)
          REAL(q), INTENT(IN)  :: XYZ_EST(:,:,:) !(1:3,1:MROW_EST,1:MTYP)
! Output variables
          REAL(q), INTENT(OUT) :: DCLM(:,:,:,:) !(1:3,1:MROW_DCLM,1:MTYP,1:MCOL_DCLM)
          REAL(q), INTENT(OUT) :: DPS_SIC(:,:,:,:) !(1:MROW_DPS_SIC,1:MCOL_DPS_SIC,1:3,1:MTYP)
          REAL(q), INTENT(OUT) :: RMEM
! Local variables
          REAL(q)              :: CLM_LOC1(1:MRB)
          REAL(q)              :: CLM_LOC2(1:MRB)
          REAL(q)              :: COSD
          REAL(q)              :: COSP
          REAL(q)              :: DCLM_LOC1(1:MRB,1:3)
          REAL(q)              :: DCLM_LOC2(1:MRB,1:3)
          REAL(q)              :: DD(1:3) ! Temporary variables used to calculate derivative of spherical harmonics.
          REAL(q)              :: DF(1:3) ! Temporary variable.
          REAL(q)              :: DFDD,DFDP,DFDR ! Temporary variables used to calculate derivative of spherical harmonics
          REAL(q)              :: DP(1:3) ! Temporary variable used to calculate derivative of spherical harmonics.
          REAL(q)              :: DRR(1:3) ! Temporary variable used to calculate derivative of sherical harmonics.
          REAL(q)              :: EPS ! Threshold prameter.
          REAL(q)              :: F(1:MRB,0:LMAX)
          REAL(q)              :: FDER(1:MRB,0:LMAX)
          REAL(q)              :: FWIN_LOC
          INTEGER              :: IB
          INTEGER              :: ILM
          INTEGER              :: INEIB
          INTEGER              :: INIONS
          INTEGER              :: INTYP
          INTEGER              :: INTYP0
          INTEGER              :: IVAR1
          INTEGER              :: IVAR2
          INTEGER              :: IXYZ
          INTEGER              :: JB,JNTYP,JNTYP0,JJNTYP,JJNTYP0 ! Same as above.
          INTEGER              :: L ! Temporary variable specifying angular momentum number
          INTEGER              :: M ! Temporary variable specifying magnetic momentum number
          INTEGER              :: NAC ! Location of the atomic position in the spline grid
          REAL(q)              :: R
          REAL(q)              :: SGN ! Temporary variable.
          REAL(q)              :: SIND
          REAL(q)              :: SIND_CRITERIA
          REAL(q)              :: SINP
          REAL(q)              :: YLM(1:(LMAX+1)**2) ! Spherical harmonics.
          REAL(q)              :: YLMD(1:(LMAX+1)**2) ! Derivative of spherical harmonics.
          REAL(q)              :: YLMP(1:(LMAX+1)**2) ! Derivative of spherical harmonics.
! Record dynamical memory usage for local variables
          RMEM=(REAL((LMAX+1),q)**2*8 &                              ! YLM
               +REAL((LMAX+1),q)**2*8 &                              ! YLMD
               +REAL((LMAX+1),q)**2*8 &                              ! YLMP
               )/1E6_q
! Set constant variables
          EPS=1E-10_q
! Calculate derivative of expansion coefficients.
          DCLM=0.0_q
          IF(LSIC) THEN
             DPS_SIC=0.0_q
          ENDIF
! Calculation of coefficients
          DO INTYP=1, NTYP
             INTYP0=LADD(INTYP)
! set up subcommunicator for NIONS
             IF(NNEIB_EST(INTYP0).GT.0) THEN
! MPI - parallel
                DO INEIB=2, NNEIB_EST(INTYP0)+1
                   INIONS=LNIONS_EST(INEIB,INTYP0)
! Compute derivatives of spherical harmonics.
                   CALL YLM1 (EPS,LMAX+1,(LMAX+1)**2,XYZ_EST(1:,INEIB,INTYP0),COSD,COSP,R,SIND,SIND_CRITERIA,SINP,YLM,YLMD,YLMP)
                   IF(ABS(SIND_CRITERIA).GT.EPS) THEN
                      DRR(1)= SIND*COSP
                      DD(1) = COSD*COSP/R
                      DP(1) =-SINP/SIND/R
                      DRR(2)=SIND*SINP
                      DD(2) =COSD*SINP/R
                      DP(2) =COSP/SIND/R
                      DRR(3)=COSD
                      DD(3) =-SIND/R
                      DP(3) =0.0_q
                   ENDIF
! Compute derivative of radial function.
                   NAC=INT(R/DRSPL)+1
                   DO L=0, LMAX
!NEC$ novector
                      DO IB=1, NRB(L)
                         IF(LWINDOW) THEN
                            FWIN_LOC=FWIN(L,IB)
                         ELSE
                            FWIN_LOC=1.0_q
                         ENDIF
                         CALL SPLVAL_ML (R,F(IB,L),FDER(IB,L),FNL(1:,1:,L,IB,KNTYP0),NAC,NSPL)
                         F(IB,L)   =FWIN_LOC*   F(IB,L)
                         FDER(IB,L)=FWIN_LOC*FDER(IB,L)
                      ENDDO
                   ENDDO
! Calculate derivatives of expansion coefficients.
                   DO L=0, LMAX
                      SGN=(-1.0_q)**L
                      DO M=1, 2*L+1
                         ILM=L**2+M
                         DO IB=1, NRB(L)
                            CLM_LOC1(IB)=SGN*F(IB,L)*YLM(ILM) ! C_nlm^ik (center is i)
                            CLM_LOC2(IB)=    F(IB,L)*YLM(ILM) ! C_nlm^ki (center is k)
                            IF(ABS(SIND_CRITERIA).GT.EPS) then
                               DFDD=   F(IB,L)*YLMD(ILM)
                               DFDP=   F(IB,L)*YLMP(ILM)
                               DFDR=FDER(IB,L)*YLM(ILM)
                               DO IXYZ=1, 3
                                  DF(IXYZ)=DRR(IXYZ)*DFDR+DD(IXYZ)*DFDD+DP(IXYZ)*DFDP
                               ENDDO
                            ELSE
                               IF(COSD .GE. 0.0_q) THEN
                                  DF(1)=    F(IB,L)*YLMD(ILM)
                                  DF(2)=    F(IB,L)*YLMP(ILM)
                                  DF(3)= FDER(IB,L)* YLM(ILM)
                               ELSE
                                  DF(1)=-   F(IB,L)*YLMD(ILM)
                                  DF(2)=-   F(IB,L)*YLMP(ILM)
                                  DF(3)=-FDER(IB,L)* YLM(ILM)
                               ENDIF
                            ENDIF
                            DO IXYZ=1, 3
                               DCLM_LOC1(IB,IXYZ)=-SGN*DF(IXYZ)  ! Derivative of C_nlm^ik with respect to Xk
                               DCLM_LOC2(IB,IXYZ)=-    DF(IXYZ)  ! Derivative of C_nlm^ki with respect to Xi
                            ENDDO
                         ENDDO
! Sum up the derivatives of expansion coefficients.
                         DO IXYZ=1, 3
                            DO IB=1, NRB(L)
                               IVAR1=(LMAX+1)**2*MTYP*(IB-1)+(LMAX+1)**2*(KNTYP0-1)+ILM
                               IVAR2=(LMAX+1)**2*MTYP*(IB-1)+(LMAX+1)**2*(INTYP0-1)+ILM
                               DCLM(IXYZ,IVAR1,INTYP0,INEIB)=DCLM(IXYZ,IVAR1,INTYP0,INEIB)+ & ! Copy dC_nlm^ki/dXk.
                               RCOUPLE_IONS(KNIONS)*DCLM_LOC1(IB,IXYZ)
                               DCLM(IXYZ,IVAR2,KNTYP0,    1)=DCLM(IXYZ,IVAR2,KNTYP0,    1)+ & ! Sum up dC_nlm^ik/dXi over i.
                               RCOUPLE_IONS(INIONS)*DCLM_LOC2(IB,IXYZ)
                            ENDDO
                         ENDDO
! Calculate derivatives of self-interaction corrections.
                         IF(LSIC) THEN
                            DO IXYZ=1, 3
                               DO IB=1, NRB(L)
                                  DO JB=IB, NRB(L)
                                     IF(LFLAG_VAR_SIC(JB,IB,L,KNTYP0,INTYP0)) THEN
                                        IVAR1=LVAR_SIC(JB,IB,L,KNTYP0,INTYP0)
                                        DPS_SIC(IVAR1,INEIB,IXYZ,INTYP0)= & ! Copy d(C_nlm^ki C_nlm^ki)/dXk.
                                        DPS_SIC(IVAR1,INEIB,IXYZ,INTYP0)+ &
                                        RCOUPLE_IONS(KNIONS)**2*(DCLM_LOC1(IB,IXYZ)*CLM_LOC1(JB)+DCLM_LOC1(JB,IXYZ)*CLM_LOC1(IB))
                                     ENDIF
                                     IF(LFLAG_VAR_SIC(JB,IB,L,INTYP0,KNTYP0)) THEN
                                        IVAR2=LVAR_SIC(JB,IB,L,INTYP0,KNTYP0)
                                        DPS_SIC(IVAR2,    1,IXYZ,KNTYP0)= & ! Sum up d(C_nlm^ik C_nlm^ik)/dXi over i.
                                        DPS_SIC(IVAR2,    1,IXYZ,KNTYP0)+ &
                                        RCOUPLE_IONS(INIONS)**2*(DCLM_LOC2(IB,IXYZ)*CLM_LOC2(JB)+DCLM_LOC2(JB,IXYZ)*CLM_LOC2(IB))
                                     ENDIF
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
        END SUBROUTINE D1CLM_SOAP

! Store l=m=0 component if necessary.
! This subroutine is used in calculations using FFM%C00_ALL.
! Except for Ceriotti's kernel, norm and inner products of DC00 are calculated.
! This subroutine is not parallelized.

        SUBROUTINE D1C00_SOAP (CLM,DC00,DCLM,LSUPERVEC,KNTYP0,LADD,LMAX,LNIONS_EST,LNIONS_EST_MAX, &
                   MCOL_CLM,MCOL_DC00,MCOL_DCLM,MRB,MROW_EST,MROW_EST_MAX,MROW_CLM,MROW_DC00,MROW_DCLM, &
                   MTYP,MTYP_ALLOCATE,NNEIB_EST,NNEIB_EST_MAX,NRB,NTYP)
          IMPLICIT NONE
! Input variables
          REAL(q), INTENT(IN)  :: CLM(:,:) !(1:MROW_CLM,1:MCOL_CLM)
          REAL(q), INTENT(IN)  :: DCLM(:,:,:,:) !(1:3,1:MROW_DCLM,1:MTYP,1:MCOL_DCLM)
          LOGICAL, INTENT(IN)  :: LSUPERVEC
          INTEGER, INTENT(IN)  :: KNTYP0
          INTEGER, INTENT(IN)  :: LADD(:) !(1:MTYP)
          INTEGER, INTENT(IN)  :: LMAX
          INTEGER, INTENT(IN)  :: LNIONS_EST(:,:) !(1:MROW_EST,1:MTYP)
          INTEGER, INTENT(IN)  :: LNIONS_EST_MAX(:,:) !(1:MROW_EST_MAX,1:MTYP)
          INTEGER, INTENT(IN)  :: MCOL_CLM
          INTEGER, INTENT(IN)  :: MCOL_DC00
          INTEGER, INTENT(IN)  :: MCOL_DCLM
          INTEGER, INTENT(IN)  :: MRB
          INTEGER, INTENT(IN)  :: MROW_CLM
          INTEGER, INTENT(IN)  :: MROW_DCLM
          INTEGER, INTENT(IN)  :: MROW_DC00
          INTEGER, INTENT(IN)  :: MROW_EST
          INTEGER, INTENT(IN)  :: MROW_EST_MAX
          INTEGER, INTENT(IN)  :: MTYP
          INTEGER, INTENT(IN)  :: MTYP_ALLOCATE
          INTEGER, INTENT(IN)  :: NNEIB_EST(:) !(1:MTYP)
          INTEGER, INTENT(IN)  :: NNEIB_EST_MAX(:) !(1:MTYP)
          INTEGER, INTENT(IN)  :: NRB(0:) !(0:LMAX)
          INTEGER, INTENT(IN)  :: NTYP
! Output variables
          REAL(q), INTENT(OUT) :: DC00(:,:,0:,:) !(1:MROW_DC00,1:MCOL_DC00,0:3,1:MTYP)
! Local variables 
          INTEGER              :: IB
          INTEGER              :: INEIB
          INTEGER              :: INIONS
          INTEGER              :: INTYP
          INTEGER              :: INTYP0
          INTEGER              :: IVAR1
          INTEGER              :: IVAR2
          INTEGER              :: IXYZ
          INTEGER              :: JNTYP
          INTEGER              :: JNTYP0
          INTEGER              :: JJNTYP0
          INTEGER, ALLOCATABLE :: LNIONS_EST_HELP(:,:)
          INTEGER, ALLOCATABLE :: NNEIB_EST_HELP(:)
!          INTEGER              :: LNIONS_EST_HELP(1:MAX(MROW_EST,MROW_EST_MAX),1:MTYP)
!          INTEGER              :: NNEIB_EST_HELP(1:MTYP)
          ! Allocate helping arrays 
          ALLOCATE(LNIONS_EST_HELP(1:MAX(MROW_EST,MROW_EST_MAX),1:MTYP_ALLOCATE))
          ALLOCATE(NNEIB_EST_HELP(1:MTYP_ALLOCATE))
! Initialization of DC00
          DC00=0.0_q
! Set lists.
! If LSUPERVEC=.FALSE., lists for *_MAX is not set. Therefore, we cannot use them.
          IF(.NOT.LSUPERVEC) THEN
             LNIONS_EST_HELP=LNIONS_EST
             NNEIB_EST_HELP=NNEIB_EST
! If LSUPERVEC=.TRUE., lists are set. We must use them.
          ELSE
             LNIONS_EST_HELP=LNIONS_EST_MAX
             NNEIB_EST_HELP=NNEIB_EST_MAX
          ENDIF
! Calculations of DC00.
          DO INTYP=1, NTYP
             INTYP0=LADD(INTYP)
             DO INEIB=1, NNEIB_EST_HELP(INTYP0)+1, 1
                INIONS=LNIONS_EST_HELP(INEIB,INTYP0)
                IF((INEIB.NE.1).OR.(INTYP0.EQ.KNTYP0)) THEN
                   DO JNTYP=1, NTYP
                      JNTYP0=LADD(JNTYP)
                      DO IB=1, NRB(0)
                         IVAR1=MRB*(JNTYP0-1)+IB
                         IVAR2=(LMAX+1)**2*MTYP*(IB-1)+(LMAX+1)**2*(JNTYP0-1)+1
                         DC00(IVAR1,INEIB,0,INTYP0)=CLM(IVAR2,INIONS)
                      ENDDO
                   ENDDO
                   IF(INEIB.LE.NNEIB_EST(INTYP0)+1) THEN
                      DO JNTYP=1, NTYP
                         JNTYP0=LADD(JNTYP)
                         DO IB=1, NRB(0)
                            IVAR1=MRB*(JNTYP0-1)+IB
                            IVAR2=(LMAX+1)**2*MTYP*(IB-1)+(LMAX+1)**2*(JNTYP0-1)+1
                            DO IXYZ=1, 3
                               DC00(IVAR1,INEIB,IXYZ,INTYP0)=DCLM(IXYZ,IVAR2,INTYP0,INEIB)
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
          ! Deallocate helping arrays
          DEALLOCATE(LNIONS_EST_HELP)
          DEALLOCATE(NNEIB_EST_HELP)
        END SUBROUTINE D1C00_SOAP

! Calculate angular descriptor within SOAP upto 1st derivatives.
! This subroutine is used in calculations using FFM%PS_ALL.
! This subroutine is not parallelized.

        SUBROUTINE D1PS_SOAP (AFILT,CLM,DCLM,DPS,DPS_SIC,IAFILT,LSUPERVEC,KNTYP0, &
                   LAFILT,LADD,LFLAG_VAR,LFLAG_VAR_SIC,L_LNRB,LMAX,LNIONS_EST,LNIONS_EST_MAX,LSIC,LVAR,LVAR_SIC, &
                   MCOL_CLM,MCOL_DCLM,MCOL_DPS,MCOL_DPS_SIC,MCOL_PS_SIC,MLNRB,MRB, &
                   MROW_EST,MROW_EST_MAX,MROW_CLM,MROW_DCLM,MROW_DPS,MROW_DPS_SIC,MROW_PS_SIC,MTYP,MTYP_ALLOCATE,&
                   NIONS,NLNRB,NNEIB_EST,NNEIB_EST_MAX,NRB,NRB_LNRB,NTYP,PS_SIC,WVAR)
          IMPLICIT NONE
! Input variables
          REAL(q), INTENT(IN)  :: AFILT
          REAL(q), INTENT(IN)  :: CLM(:,:) !(1:MROW_CLM,1:MCOL_CLM)
          REAL(q), INTENT(IN)  :: DCLM(:,:,:,:) !(1:3,1:MROW_DCLM,1:MTYP,1:MCOL_DCLM)
          REAL(q), INTENT(IN)  :: DPS_SIC(:,:,:,:) !(1:MROW_DPS_SIC,1:MCOL_DPS_SIC,1:3,1:MTYP)
          INTEGER, INTENT(IN)  :: KNTYP0
          INTEGER, INTENT(IN)  :: IAFILT
          LOGICAL, INTENT(IN)  :: LSUPERVEC
          LOGICAL, INTENT(IN)  :: LAFILT
          INTEGER, INTENT(IN)  :: LADD(:) !(1:MTYP)
          LOGICAL, INTENT(IN)  :: LFLAG_VAR(:,:,0:,:,:,:) !(1:MRB,1:MRB,0:LMAX,1:MTYP,1:MTYP,1:MTYP)
          LOGICAL, INTENT(IN)  :: LFLAG_VAR_SIC(:,:,0:,:,:) !(1:MRB,1:MRB,0:LMAX,1:MTYP,1:MTYP)
          INTEGER, INTENT(IN)  :: L_LNRB(:) !(1:MLNRB)
          INTEGER, INTENT(IN)  :: LMAX
          INTEGER, INTENT(IN)  :: LNIONS_EST(:,:) !(1:MROW_EST,1:MTYP)
          INTEGER, INTENT(IN)  :: LNIONS_EST_MAX(:,:) !(1:MROW_EST_MAX,1:MTYP)
          LOGICAL, INTENT(IN)  :: LSIC
          INTEGER, INTENT(IN)  :: LVAR(:,:,0:,:,:,:) !(1:MRB,1:MRB,0:LMAX,1:MTYP,1:MTYP,1:MTYP)
          INTEGER, INTENT(IN)  :: LVAR_SIC(:,:,0:,:,:) !(1:MRB,1:MRB,0:LMAX,1:MTYP,1:MTYP)
          INTEGER, INTENT(IN)  :: MCOL_CLM
          INTEGER, INTENT(IN)  :: MCOL_DCLM
          INTEGER, INTENT(IN)  :: MCOL_DPS
          INTEGER, INTENT(IN)  :: MCOL_DPS_SIC
          INTEGER, INTENT(IN)  :: MCOL_PS_SIC
          INTEGER, INTENT(IN)  :: MLNRB
          INTEGER, INTENT(IN)  :: MRB
          INTEGER, INTENT(IN)  :: MROW_CLM
          INTEGER, INTENT(IN)  :: MROW_DCLM
          INTEGER, INTENT(IN)  :: MROW_DPS
          INTEGER, INTENT(IN)  :: MROW_DPS_SIC
          INTEGER, INTENT(IN)  :: MROW_EST
          INTEGER, INTENT(IN)  :: MROW_EST_MAX
          INTEGER, INTENT(IN)  :: MROW_PS_SIC
          INTEGER, INTENT(IN)  :: MTYP
          INTEGER, INTENT(IN)  :: MTYP_ALLOCATE
          INTEGER, INTENT(IN)  :: NIONS
          INTEGER, INTENT(IN)  :: NLNRB
          INTEGER, INTENT(IN)  :: NNEIB_EST(:) !(1:MTYP)
          INTEGER, INTENT(IN)  :: NNEIB_EST_MAX(:) !(1:MTYP)
          INTEGER, INTENT(IN)  :: NRB(0:) !(0:LMAX)
          INTEGER, INTENT(IN)  :: NRB_LNRB(:) !(1:MLNRB)
          INTEGER, INTENT(IN)  :: NTYP
          REAL(q), INTENT(IN)  :: PS_SIC(:,:) !(1:MROW_PS_SIC,1:MCOL_PS_SIC)
          REAL(q), INTENT(IN)  :: WVAR(:,:,0:,:,:,:) !(1:MRB,1:MRB,0:LMAX,1:MTYP,1:MTYP,1:MTYP)
! Output variables
          REAL(q), INTENT(OUT) :: DPS(:,:,0:,:) !(1:MROW_DPS,1:MCOL_DPS,0:3,1:MTYP)
! Local variables
          REAL(q)              :: DPS_LOC(0:3)
          INTEGER              :: IB
          INTEGER              :: ILM
          INTEGER              :: ILNRB
          INTEGER              :: INEIB
          INTEGER              :: INIONS
          INTEGER              :: INTYP
          INTEGER              :: INTYP0
          INTEGER              :: IVAR1
          INTEGER              :: IVAR2
          INTEGER              :: IVAR3
          INTEGER              :: IVAR4
          INTEGER              :: IXYZ
          INTEGER              :: JB
          INTEGER              :: JNTYP
          INTEGER              :: JNTYP0
          INTEGER              :: JJNTYP
          INTEGER              :: JJNTYP0 ! Same as above.
!          INTEGER              :: LNIONS_EST_HELP(1:MAX(MROW_EST,MROW_EST_MAX),1:MTYP)
          INTEGER, ALLOCATABLE :: LNIONS_EST_HELP(:,:)
          INTEGER              :: LNTYP0
          INTEGER              :: LLNTYP0
          INTEGER              :: L ! Temporary variable specifying angular momentum number
          INTEGER              :: M ! Temporary variable specifying magnetic momentum number
!          INTEGER              :: NNEIB_EST_HELP(1:MTYP)
          INTEGER, ALLOCATABLE :: NNEIB_EST_HELP(:)
          REAL(q)              :: PI8 ! Constant 8*PI**2
          REAL(q)              :: PREFAC(0:LMAX)
          ! Allocate helping arrays
          ALLOCATE(LNIONS_EST_HELP(1:MAX(MROW_EST,MROW_EST_MAX),1:MTYP_ALLOCATE))
          ALLOCATE(NNEIB_EST_HELP(1:MTYP_ALLOCATE))
! Set lists.
! If LSUPERVEC=.FALSE., lists for *_MAX is not set. Therefore, we cannot use them.
          IF(.NOT.LSUPERVEC) THEN
             LNIONS_EST_HELP=LNIONS_EST
             NNEIB_EST_HELP=NNEIB_EST
! If LSUPERVEC=.TRUE., lists are set. We must use them.
          ELSE
             LNIONS_EST_HELP=LNIONS_EST_MAX
             NNEIB_EST_HELP=NNEIB_EST_MAX
          ENDIF
! Set constant variables.
          PI8=8.0_q*PI**2
          DO L=0, LMAX
             IF(LAFILT) THEN
                IF(IAFILT.EQ.1) THEN
                   PREFAC(L)=SQRT(PI8/REAL((2*L+1),q))/SQRT(REAL((2*L+1),q))
                ELSE IF(IAFILT.EQ.2) THEN
                   PREFAC(L)=SQRT(PI8/REAL((2*L+1),q))/(1.0_q+AFILT*REAL((L*(L+1)),q)**2)**2
                ENDIF
             ELSE
                PREFAC(L)=SQRT(PI8/REAL((2*L+1),q))
             ENDIF
          ENDDO
! Initialization
          DPS=0.0_q
! calculation of DPS
          DO INTYP=1, NTYP
             INTYP0=LADD(INTYP)
             DO INEIB=1, NNEIB_EST_HELP(INTYP0)+1, 1
                INIONS=LNIONS_EST_HELP(INEIB,INTYP0)
                IF((INEIB.NE.1).OR.(INTYP0.EQ.KNTYP0)) THEN
                   DO JNTYP=1, NTYP
                      JNTYP0=LADD(JNTYP)
                      DO JJNTYP=1, NTYP
                         JJNTYP0=LADD(JJNTYP)
                         DO ILNRB=1, NLNRB
                            L=L_LNRB(ILNRB)
                            IB=NRB_LNRB(ILNRB)
                            DO JB=IB, NRB(L)
                               IF(LFLAG_VAR(JB,IB,L,JJNTYP0,JNTYP0,INTYP0)) THEN
                                  IVAR1=LVAR(JB,IB,L,JJNTYP0,JNTYP0,INTYP0)
                                  DPS_LOC(0)=0.0_q
                                  DO M=1, 2*L+1
                                     ILM=L**2+M
                                     IVAR2=(LMAX+1)**2*MTYP*(IB-1)+(LMAX+1)**2*( JNTYP0-1)+ILM
                                     IVAR3=(LMAX+1)**2*MTYP*(JB-1)+(LMAX+1)**2*(JJNTYP0-1)+ILM
                                     DPS_LOC(0)=DPS_LOC(0)+CLM(IVAR2,INIONS)*CLM(IVAR3,INIONS)
                                  ENDDO
                                  IF(LFLAG_VAR_SIC(JB,IB,L,JJNTYP0,INTYP0).AND.LSIC.AND.(JJNTYP0.EQ.JNTYP0)) THEN
                                     IVAR4=LVAR_SIC(JB,IB,L,JJNTYP0,INTYP0)
                                     DPS_LOC(0)=DPS_LOC(0)-PS_SIC(IVAR4,INIONS)
                                  ENDIF
                                  DPS(IVAR1,INEIB,0,INTYP0)=DPS(IVAR1,INEIB,0,INTYP0)+WVAR(JB,IB,L,JJNTYP0,JNTYP0,INTYP0)*PREFAC(L)*DPS_LOC(0)
                               ENDIF
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                   IF(INEIB.LE.NNEIB_EST(INTYP0)+1) THEN
                      DO JNTYP=1, NTYP
                         JNTYP0=LADD(JNTYP)
                         DO JJNTYP=1, NTYP
                            JJNTYP0=LADD(JJNTYP)
                            DO ILNRB=1, NLNRB
                               L=L_LNRB(ILNRB)
                               IB=NRB_LNRB(ILNRB)
                               DO JB=IB, NRB(L)
                                  IF(LFLAG_VAR(JB,IB,L,JJNTYP0,JNTYP0,INTYP0)) THEN
                                     IVAR1=LVAR(JB,IB,L,JJNTYP0,JNTYP0,INTYP0)
                                     DO IXYZ=1, 3
                                        DPS_LOC(IXYZ)=0.0_q
                                     ENDDO
                                     DO M=1, 2*L+1
                                        ILM=L**2+M
                                        IVAR2=(LMAX+1)**2*MTYP*(IB-1)+(LMAX+1)**2*( JNTYP0-1)+ILM
                                        IVAR3=(LMAX+1)**2*MTYP*(JB-1)+(LMAX+1)**2*(JJNTYP0-1)+ILM
                                        DO IXYZ=1, 3
                                           DPS_LOC(IXYZ)=DPS_LOC(IXYZ)+(DCLM(IXYZ,IVAR3,INTYP0,INEIB)*CLM(IVAR2,INIONS)+ &
                                                                        DCLM(IXYZ,IVAR2,INTYP0,INEIB)*CLM(IVAR3,INIONS))
                                        ENDDO
                                     ENDDO
                                     IF(LFLAG_VAR_SIC(JB,IB,L,JJNTYP0,INTYP0).AND.LSIC.AND.(JJNTYP0.EQ.JNTYP0)) THEN
                                        DO IXYZ=1, 3
                                           IVAR4=LVAR_SIC(JB,IB,L,JJNTYP0,INTYP0)
                                           DPS_LOC(IXYZ)=DPS_LOC(IXYZ)-DPS_SIC(IVAR4,INEIB,IXYZ,INTYP0)
                                        ENDDO
                                     ENDIF
                                     DO IXYZ=1, 3
                                        DPS(IVAR1,INEIB,IXYZ,INTYP0)=DPS(IVAR1,INEIB,IXYZ,INTYP0)+WVAR(JB,IB,L,JJNTYP0,JNTYP0,INTYP0)*PREFAC(L)*DPS_LOC(IXYZ)
                                     ENDDO
                                  ENDIF
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
          ! Deallocate helping arrays
          DEALLOCATE(LNIONS_EST_HELP)
          DEALLOCATE(NNEIB_EST_HELP)
        END SUBROUTINE D1PS_SOAP

! Normalization of derivatives.

        SUBROUTINE NORMALIZATION_D1 (DC00,DPS,LSUPERVEC,KNTYP0,LADD,LNIONS_EST1,LNIONS_EST2,LNIONS_EST_MAX,LNITYP,LNORM1,LNORM2, &
                   MCOL_DC00,MCOL_DPS,MCOL_RNORM_ALL,MCOL_RNORM_C00,MCOL_RNORM_PS,MIONS,MITYP, &
                   MROW_EST1,MROW_EST2,MROW_EST_MAX,MROW_DC00,MROW_DPS,MROW_RNORM_C00,MROW_RNORM_PS,MTYP, &
                   NNEIB_EST1,NNEIB_EST2,NNEIB_EST_MAX,NNVAR1,NNVAR2,NTYP,RCUT1,RCUT2,RNORM_ALL,RNORM_C00,RNORM_PS,W1,W2)
          IMPLICIT NONE
! Input variables
          LOGICAL, INTENT(IN) :: LSUPERVEC
          INTEGER, INTENT(IN) :: KNTYP0
          INTEGER, INTENT(IN) :: LADD(:) !(1:MTYP)
          INTEGER, INTENT(IN) :: LNIONS_EST1(:,:) !(1:MROW_EST1,1:MTYP)
          INTEGER, INTENT(IN) :: LNIONS_EST2(:,:) !(1:MROW_EST2,1:MTYP)
          INTEGER, INTENT(IN) :: LNIONS_EST_MAX(:,:) !(1:MROW_EST_MAX,1:MTYP)
          INTEGER, INTENT(IN) :: LNITYP(:) !(1:MIONS)
          LOGICAL, INTENT(IN) :: LNORM1
          LOGICAL, INTENT(IN) :: LNORM2
          INTEGER, INTENT(IN) :: MCOL_DC00
          INTEGER, INTENT(IN) :: MCOL_DPS
          INTEGER, INTENT(IN) :: MCOL_RNORM_ALL
          INTEGER, INTENT(IN) :: MCOL_RNORM_C00
          INTEGER, INTENT(IN) :: MCOL_RNORM_PS
          INTEGER, INTENT(IN) :: MIONS
          INTEGER, INTENT(IN) :: MITYP
          INTEGER, INTENT(IN) :: MROW_EST1
          INTEGER, INTENT(IN) :: MROW_EST2
          INTEGER, INTENT(IN) :: MROW_EST_MAX
          INTEGER, INTENT(IN) :: MROW_DC00
          INTEGER, INTENT(IN) :: MROW_DPS
          INTEGER, INTENT(IN) :: MROW_RNORM_C00
          INTEGER, INTENT(IN) :: MROW_RNORM_PS
          INTEGER, INTENT(IN) :: MTYP
          INTEGER, INTENT(IN) :: NNEIB_EST1(:) !(1:MTYP)
          INTEGER, INTENT(IN) :: NNEIB_EST2(:) !(1:MTYP)
          INTEGER, INTENT(IN) :: NNEIB_EST_MAX(:) !(1:MTYP)
          INTEGER, INTENT(IN) :: NNVAR1
          INTEGER, INTENT(IN) :: NNVAR2(:) !(1:MTYP)
          INTEGER, INTENT(IN) :: NTYP
          REAL(q), INTENT(IN) :: RCUT1
          REAL(q), INTENT(IN) :: RCUT2
          REAL(q), INTENT(IN) :: RNORM_ALL(:,:) !(1:MROW_RNORM_ALL,1:MCOL_RNORM_ALL)
          REAL(q), INTENT(IN) :: RNORM_C00(:,:) !(1:MROW_RNORM_C00,1:MCOL_RNORM_C00)
          REAL(q), INTENT(IN) :: RNORM_PS(:,:) !(1:MROW_RNORM_PS,1:MCOL_RNORM_PS)
          REAL(q), INTENT(IN) :: W1
          REAL(q), INTENT(IN) :: W2
! In and output variables.
          REAL(q), INTENT(INOUT) :: DC00(:,:,0:,:) !(1:MROW_DC00,1:MCOL_DC00,0:3,1:MTYP)
          REAL(q), INTENT(INOUT) :: DPS(:,:,0:,:) !(1:MROW_DPS,1:MCOL_DPS,0:3,1:MTYP)
! Local variables.
          INTEGER             :: INEIB
          INTEGER             :: INIONS
          INTEGER             :: INITYP
          INTEGER             :: INTYP
          INTEGER             :: INTYP0
          INTEGER             :: IVAR
          INTEGER             :: IXYZ
          REAL(q)             :: PRODUCT_HELP1
          REAL(q)             :: PRODUCT_HELP2(1:3,1:MAX(MROW_EST1,MROW_EST2),1:MTYP)
          REAL(q)             :: RNORM_ALL_HELP
          REAL(q)             :: RNORM_C00_HELP
          REAL(q)             :: RNORM_PS_HELP
! LAPACK routine
          REAL(q)             :: DNRM2,DDOT
          ! Allocate helpint arrays
          ! Do calculations
          IF(.NOT.LSUPERVEC) THEN
             IF(W1.GT.0.0_q) THEN
                IF(LNORM1) THEN
                   DO INTYP=1, NTYP
                      INTYP0=LADD(INTYP)
                      DO INEIB=1, NNEIB_EST1(INTYP0)+1, 1
                         IF((INEIB.NE.1).OR.(INTYP0.EQ.KNTYP0)) THEN
                            INIONS=LNIONS_EST1(INEIB,INTYP0)
                            INITYP=LNITYP(INIONS)
                            RNORM_C00_HELP=1.0_q/SQRT(RNORM_C00(INITYP,INTYP0))
                            DO IXYZ=0, 3
                               CALL DSCAL(NNVAR1,RNORM_C00_HELP,DC00(1,INEIB,IXYZ,INTYP0),1)
                               IF (IXYZ.GT.0) THEN
                                  PRODUCT_HELP1=DDOT(NNVAR1,DC00(1,INEIB,0,INTYP0),1,DC00(1,INEIB,IXYZ,INTYP0),1)
                                  DO IVAR=1,NNVAR1
                                     DC00(IVAR,INEIB,IXYZ,INTYP0) = &
                                     DC00(IVAR,INEIB,IXYZ,INTYP0) - &
                                     DC00(IVAR,INEIB,0,INTYP0)*PRODUCT_HELP1
                                  ENDDO
                               ENDIF
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF
             ENDIF
             IF(W2.GT.0.0_q) THEN
                IF(LNORM2) THEN
                   DO INTYP=1, NTYP
                      INTYP0=LADD(INTYP)
                      DO INEIB=1, NNEIB_EST2(INTYP0)+1, 1
                         IF((INEIB.NE. 1).OR.(INTYP0.EQ.KNTYP0)) THEN
                            INIONS=LNIONS_EST2(INEIB,INTYP0)
                            INITYP=LNITYP(INIONS)
                            RNORM_PS_HELP=1.0_q/SQRT(RNORM_PS(INITYP,INTYP0))
                            DO IXYZ=0, 3
                               CALL DSCAL(NNVAR2(INTYP0),RNORM_PS_HELP,DPS(1,INEIB,IXYZ,INTYP0),1)
                               IF (IXYZ.GT.0) THEN
                                  PRODUCT_HELP1=DDOT(NNVAR2(INTYP0),DPS(1,INEIB,0,INTYP0),1,DPS(1,INEIB,IXYZ,INTYP0),1)
                                  DO IVAR=1,NNVAR2(INTYP0)
                                     DPS(IVAR,INEIB,IXYZ,INTYP0) = &
                                     DPS(IVAR,INEIB,IXYZ,INTYP0) - &
                                     DPS(IVAR,INEIB,0,INTYP0)*PRODUCT_HELP1
                                  ENDDO
                               ENDIF
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF
             ENDIF
          ELSE
! Normalize the super vectors and calculate products.
! The INEIB is listed so that the list for C00 is the same as that for PS inside MIN(RCUT1,RCUT2).
! By this reason, PRODUCT_HELP2 can be calculated simply by the following method.
             PRODUCT_HELP2=0.0_q
             IF(W1.GT.0.0_q) THEN
                IF(LNORM1) THEN
                   DO INTYP=1, NTYP
                      INTYP0=LADD(INTYP)
                      DO INEIB=1, NNEIB_EST_MAX(INTYP0)+1, 1
                         IF((INEIB.NE.1).OR.(INTYP0.EQ.KNTYP0)) THEN
                            INIONS=LNIONS_EST_MAX(INEIB,INTYP0)
                            INITYP=LNITYP(INIONS)
                            RNORM_ALL_HELP=SQRT(W1)/SQRT(RNORM_ALL(INITYP,INTYP0))
                            CALL DSCAL(NNVAR1,RNORM_ALL_HELP,DC00(1,INEIB,0,INTYP0),1)
                            IF(INEIB.LE.NNEIB_EST1(INTYP0)+1) THEN
                               DO IXYZ=1, 3
                                  CALL DSCAL(NNVAR1,RNORM_ALL_HELP,DC00(1,INEIB,IXYZ,INTYP0),1)
                                  IF (IXYZ.GT.0) THEN
                                     PRODUCT_HELP2(IXYZ,INEIB,INTYP0)=PRODUCT_HELP2(IXYZ,INEIB,INTYP0)+DDOT(NNVAR1,DC00(1,INEIB,0,INTYP0),1,DC00(1,INEIB,IXYZ,INTYP0),1)
                                  ENDIF
                               ENDDO
                            ENDIF
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF
             ENDIF
             IF(W2.GT.0.0_q) THEN
                IF(LNORM2) THEN
                   DO INTYP=1, NTYP
                      INTYP0=LADD(INTYP)
                      DO INEIB=1, NNEIB_EST_MAX(INTYP0)+1, 1
                         IF((INEIB.NE.1).OR.(INTYP0.EQ.KNTYP0)) THEN
                            INIONS=LNIONS_EST_MAX(INEIB,INTYP0)
                            INITYP=LNITYP(INIONS)
                            RNORM_ALL_HELP=SQRT(W2)/SQRT(RNORM_ALL(INITYP,INTYP0))
                            CALL DSCAL(NNVAR2(INTYP0),RNORM_ALL_HELP,DPS(1,INEIB,0,INTYP0),1)
                            IF(INEIB.LE.NNEIB_EST2(INTYP0)+1) THEN
                               DO IXYZ=1, 3
                                  CALL DSCAL(NNVAR2(INTYP0),RNORM_ALL_HELP,DPS(1,INEIB,IXYZ,INTYP0),1)
                                  IF (IXYZ.GT.0) THEN
                                     PRODUCT_HELP2(IXYZ,INEIB,INTYP0)=PRODUCT_HELP2(IXYZ,INEIB,INTYP0)+DDOT(NNVAR2(INTYP0),DPS(1,INEIB,0,INTYP0),1,DPS(1,INEIB,IXYZ,INTYP0),1)
                                  ENDIF
                               ENDDO
                            ENDIF
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF
             ENDIF
! Finalize the normalization.
             IF(W1.GT.0.0_q) THEN
                IF(LNORM1) THEN
                   DO INTYP=1, NTYP
                      INTYP0=LADD(INTYP)
                      DO INEIB=1, NNEIB_EST_MAX(INTYP0)+1, 1
                         IF((INEIB.NE.1).OR.(INTYP0.EQ.KNTYP0)) THEN
                            DO IXYZ=1, 3
                               DO IVAR=1,NNVAR1
                                  DC00(IVAR,INEIB,IXYZ,INTYP0) = &
                                  DC00(IVAR,INEIB,IXYZ,INTYP0) - &
                                  DC00(IVAR,INEIB,0,INTYP0)*PRODUCT_HELP2(IXYZ,INEIB,INTYP0)
                               ENDDO
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF
             ENDIF
             IF(W2.GT.0.0_q) THEN
                IF(LNORM2) THEN
                   DO INTYP=1, NTYP
                      INTYP0=LADD(INTYP)
                      DO INEIB=1, NNEIB_EST_MAX(INTYP0)+1, 1
                         IF((INEIB.NE. 1).OR.(INTYP0.EQ.KNTYP0)) THEN
                            DO IXYZ=1, 3
                               DO IVAR=1,NNVAR2(INTYP0)
                                  DPS(IVAR,INEIB,IXYZ,INTYP0) = &
                                  DPS(IVAR,INEIB,IXYZ,INTYP0) - &
                                  DPS(IVAR,INEIB,0,INTYP0)*PRODUCT_HELP2(IXYZ,INEIB,INTYP0)
                               ENDDO
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF
             ENDIF
          ENDIF
        END SUBROUTINE NORMALIZATION_D1

#ifndef ML_LEARNING_NOT_ACTIVE
! Calculate polynomical of inner product.
! This subroutine is used in calculations using distributed angular descriptor.
! ScaLAPACK-parallel.

        SUBROUTINE POLY_SOAP_INIT1 (CONTEXT_WORLD,DESCRIPTOR,INIT,LSUPERVEC,MX_LOC,MY_LOC,N,NX,NY,NHYP,SOAP1,SOAP2,W)
          USE ML_FF_STRUCT, ONLY : SCALAPACK_GRID, PARALLEL_DESCRIPTOR
          IMPLICIT NONE
! Descriptor
          TYPE (SCALAPACK_GRID)      :: CONTEXT_WORLD
          TYPE (PARALLEL_DESCRIPTOR) :: DESCRIPTOR
! Input variables
          INTEGER, INTENT(IN)  :: INIT
          LOGICAL, INTENT(IN)  :: LSUPERVEC
          INTEGER, INTENT(IN)  :: MX_LOC,MY_LOC
          INTEGER, INTENT(IN)  :: N
          INTEGER, INTENT(IN)  :: NX,NY
          INTEGER, INTENT(IN)  :: NHYP
          REAL(q), INTENT(IN)  :: SOAP1(:,:) !(1:MX_LOC,1:MY_LOC)
          REAL(q), INTENT(IN)  :: W
! Output variables
          REAL(q), INTENT(OUT) :: SOAP2(:,:) !(1:MX_LOC,1:MY_LOC)
! Local variables
          INTEGER              :: IX
          INTEGER              :: IY
          INTEGER              :: LOCVAR_COL
          INTEGER              :: LOCVAR_ROW
          REAL(q)              :: PREFAC
          INTEGER              :: TEST_COL
          INTEGER              :: TEST_ROW
! Constant.
          PREFAC=W/DBLE(N)
! Initialization if necessary
          DO IY=1, NY
#ifdef scaLAPACK
             TEST_COL=MOD((IY-1)/DESCRIPTOR%NB,CONTEXT_WORLD%NP_COL)
             IF(CONTEXT_WORLD%ME_COL.NE.TEST_COL) THEN
               CYCLE
             ENDIF
             LOCVAR_COL=GLOBAL2LOCAL(IY,CONTEXT_WORLD%NP_COL,DESCRIPTOR%NB)
#else
             LOCVAR_COL=IY
#endif
             DO IX=1, NX
#ifdef scaLAPACK
                TEST_ROW=MOD((IX-1)/DESCRIPTOR%MB,CONTEXT_WORLD%NP_ROW)
                IF(CONTEXT_WORLD%ME_ROW.NE.TEST_ROW) THEN
                  CYCLE
                ENDIF
                LOCVAR_ROW=GLOBAL2LOCAL(IX,CONTEXT_WORLD%NP_ROW,DESCRIPTOR%MB)
#else
                LOCVAR_ROW=IX
#endif
                IF(INIT.EQ.0) THEN
                   SOAP2(LOCVAR_ROW,LOCVAR_COL)=0.0_q
                ENDIF
                IF(.NOT.LSUPERVEC) THEN
                   SOAP2(LOCVAR_ROW,LOCVAR_COL)=SOAP2(LOCVAR_ROW,LOCVAR_COL)+PREFAC*SOAP1(LOCVAR_ROW,LOCVAR_COL)**NHYP
                ELSE
                   SOAP2(LOCVAR_ROW,LOCVAR_COL)=SOAP2(LOCVAR_ROW,LOCVAR_COL)+SOAP1(LOCVAR_ROW,LOCVAR_COL)
                ENDIF
             ENDDO
          ENDDO
        END SUBROUTINE POLY_SOAP_INIT1
#endif

#ifndef ML_LEARNING_NOT_ACTIVE
! Calculate polynomical of inner product.
! This subroutine is used in calculations using distributed angular descriptor.
! ScaLAPACK-parallel.

        SUBROUTINE POLY_SOAP_INIT2 (CONTEXT_WORLD,DESCRIPTOR,MX_LOC,MY_LOC,N,NX,NY,NHYP,SOAP)
          USE ML_FF_STRUCT, ONLY : SCALAPACK_GRID, PARALLEL_DESCRIPTOR
          IMPLICIT NONE
! Descriptor
          TYPE (SCALAPACK_GRID)      :: CONTEXT_WORLD
          TYPE (PARALLEL_DESCRIPTOR) :: DESCRIPTOR
! Input variables
          INTEGER, INTENT(IN)  :: MX_LOC,MY_LOC
          INTEGER, INTENT(IN)  :: N
          INTEGER, INTENT(IN)  :: NX,NY
          INTEGER, INTENT(IN)  :: NHYP
! Output variables
          REAL(q), INTENT(INOUT) :: SOAP(:,:) !(1:MX_LOC,1:MY_LOC)
! Local variables
          INTEGER              :: IX
          INTEGER              :: IY
          INTEGER              :: LOCVAR_COL
          INTEGER              :: LOCVAR_ROW
          INTEGER              :: TEST_COL
          INTEGER              :: TEST_ROW
! Initialization if necessary
          DO IY=1, NY
#ifdef scaLAPACK
             TEST_COL=MOD((IY-1)/DESCRIPTOR%NB,CONTEXT_WORLD%NP_COL)
             IF(CONTEXT_WORLD%ME_COL.NE.TEST_COL) THEN
               CYCLE
             ENDIF
             LOCVAR_COL=GLOBAL2LOCAL(IY,CONTEXT_WORLD%NP_COL,DESCRIPTOR%NB)
#else
             LOCVAR_COL=IY
#endif
             DO IX=1, NX
#ifdef scaLAPACK
                TEST_ROW=MOD((IX-1)/DESCRIPTOR%MB,CONTEXT_WORLD%NP_ROW)
                IF(CONTEXT_WORLD%ME_ROW.NE.TEST_ROW) THEN
                  CYCLE
                ENDIF
                LOCVAR_ROW=GLOBAL2LOCAL(IX,CONTEXT_WORLD%NP_ROW,DESCRIPTOR%MB)
#else
                LOCVAR_ROW=IX
#endif
                SOAP(LOCVAR_ROW,LOCVAR_COL)=SOAP(LOCVAR_ROW,LOCVAR_COL)**NHYP
             ENDDO
          ENDDO
        END SUBROUTINE POLY_SOAP_INIT2
#endif

! This subroutine is used in calculations using distributed angular descriptor.
! ScaLAPACK-parallel.

        SUBROUTINE POLY_SOAP1 (INIT,LSUPERVEC,MX_LOC,MY_LOC,N,NX,NY,NHYP,SOAP1,SOAP2,W)
          IMPLICIT NONE
! Input variables
          INTEGER, INTENT(IN)  :: INIT
          LOGICAL, INTENT(IN)  :: LSUPERVEC
          INTEGER, INTENT(IN)  :: MX_LOC,MY_LOC
          INTEGER, INTENT(IN)  :: N
          INTEGER, INTENT(IN)  :: NX,NY
          INTEGER, INTENT(IN)  :: NHYP
          REAL(q), INTENT(IN)  :: SOAP1(:,:) !(1:MX_LOC,1:MY_LOC)
          REAL(q), INTENT(IN)  :: W
! Output variables
          REAL(q), INTENT(INOUT) :: SOAP2(:,:) !(1:MX_LOC,1:MY_LOC)
! Local variables
          INTEGER              :: IX
          INTEGER              :: IY
          REAL(q)              :: PREFAC
! Constant.
          PREFAC=W/REAL(N,q)
! Initialization if necessary
          DO IY=1, NY
             DO IX=1, NX
                IF(INIT.EQ.0) THEN
                   SOAP2(IX,IY)=0.0_q
                ENDIF
                IF(.NOT.LSUPERVEC) THEN
                   SOAP2(IX,IY)=SOAP2(IX,IY)+PREFAC*SOAP1(IX,IY)**NHYP
                ELSE
                   SOAP2(IX,IY)=SOAP2(IX,IY)+SOAP1(IX,IY)
                ENDIF
             ENDDO
          ENDDO
        END SUBROUTINE POLY_SOAP1

! This subroutine is used in calculations using distributed angular descriptor.
! ScaLAPACK-parallel.

        SUBROUTINE POLY_SOAP2 (MX,MY,N,NX,NY,NHYP,SOAP)
          IMPLICIT NONE
! Input variables
          INTEGER, INTENT(IN)  :: MX,MY
          INTEGER, INTENT(IN)  :: N
          INTEGER, INTENT(IN)  :: NX,NY
          INTEGER, INTENT(IN)  :: NHYP
! Output variables
          REAL(q), INTENT(INOUT) :: SOAP(:,:) !(1:MX,1:MY)
! Local variables
          INTEGER              :: IX
          INTEGER              :: IY
! Initialization if necessary
          DO IY=1, NY
             DO IX=1, NX
                SOAP(IX,IY)=(SOAP(IX,IY)**NHYP)/REAL(N,q)
             ENDDO
          ENDDO
        END SUBROUTINE POLY_SOAP2


! This subroutine is used in calculations of the derivative of kernel with respect to the coupling constant.
! ScaLAPACK-parallel.

        SUBROUTINE POLY_SOAP_COUPLE1 (INIT,LSUPERVEC,MX,MY,MX_COUPLE,MY_COUPLE,NX,NY,NHYP,SOAP1,SOAP1_COUPLE,SOAP2_COUPLE,W)
          IMPLICIT NONE
! Input variables
          INTEGER, INTENT(IN)  :: INIT
          LOGICAL, INTENT(IN)  :: LSUPERVEC
          INTEGER, INTENT(IN)  :: MX       ,MY
          INTEGER, INTENT(IN)  :: MX_COUPLE,MY_COUPLE
          INTEGER, INTENT(IN)  :: NX,NY
          INTEGER, INTENT(IN)  :: NHYP
          REAL(q), INTENT(IN)  :: SOAP1(:,:) !(1:MX,1:MY)
          REAL(q), INTENT(IN)  :: SOAP1_COUPLE(:,:) !(1:MX_COUPLE,1:MY_COUPLE)
          REAL(q), INTENT(IN)  :: W
! Output variables
          REAL(q), INTENT(OUT) :: SOAP2_COUPLE(:,:) !(1:MX_COUPLE,1:MY_COUPLE)
! Local variables
          INTEGER              :: IX
          INTEGER              :: IY
          REAL(q)              :: PREFAC
! Constant.
          PREFAC=REAL(NHYP,q)*W
! Initialization if necessary
          DO IY=1, NY
             DO IX=1, NX
                IF(INIT.EQ.0) THEN
                   SOAP2_COUPLE(IX,IY)=0.0_q
                ENDIF
                IF(.NOT.LSUPERVEC) THEN
                   SOAP2_COUPLE(IX,IY)=SOAP2_COUPLE(IX,IY)+PREFAC*SOAP1_COUPLE(IX,IY)*SOAP1(IX,IY)**(NHYP-1)
                ELSE
                   SOAP2_COUPLE(IX,IY)=SOAP2_COUPLE(IX,IY)+SOAP1_COUPLE(IX,IY)
                ENDIF
             ENDDO
          ENDDO
        END SUBROUTINE POLY_SOAP_COUPLE1

! This subroutine is used in calculations of the derivative of kernel with respect to the coupling constant.
! ScaLAPACK-parallel.

        SUBROUTINE POLY_SOAP_COUPLE2 (MX,MY,MX_COUPLE,MY_COUPLE,NX,NY,NHYP,SOAP,SOAP_COUPLE)
          IMPLICIT NONE
! Input variables
          INTEGER, INTENT(IN)  :: MX       ,MY
          INTEGER, INTENT(IN)  :: MX_COUPLE,MY_COUPLE
          INTEGER, INTENT(IN)  :: NX,NY
          INTEGER, INTENT(IN)  :: NHYP
          REAL(q), INTENT(IN)  :: SOAP(:,:) !(1:MX,1:MY)
! Output variables
          REAL(q), INTENT(INOUT) :: SOAP_COUPLE(:,:) !(1:MX_COUPLE,1:MY_COUPLE)
! Local variables
          INTEGER              :: IX
          INTEGER              :: IY
          REAL(q)              :: PREFAC
! Constant.
          PREFAC=REAL(NHYP,q)
! Initialization if necessary
          DO IY=1, NY
             DO IX=1, NX
                SOAP_COUPLE(IX,IY)=PREFAC*SOAP_COUPLE(IX,IY)*(SOAP(IX,IY)**(NHYP-1))
             ENDDO
          ENDDO
        END SUBROUTINE POLY_SOAP_COUPLE2

        SUBROUTINE POLY_DSOAP (INTYP0,KNIONS,KNTYP0,LNIONS_EST,MROW_EST,MROW_DSOAP,MCOL_DSOAP,MROW_SOAP,MCOL_SOAP,NX,NY,DSOAP1,DSOAP2,SOAP,W)
          IMPLICIT NONE
! Input variables
          INTEGER, INTENT(IN) :: INTYP0
          INTEGER, INTENT(IN) :: KNIONS
          INTEGER, INTENT(IN) :: KNTYP0
          INTEGER, INTENT(IN) :: LNIONS_EST(:) !(1:MROW_EST)
          INTEGER, INTENT(IN) :: MROW_EST
          INTEGER, INTENT(IN) :: MROW_DSOAP
          INTEGER, INTENT(IN) :: MCOL_DSOAP
          INTEGER, INTENT(IN) :: MROW_SOAP 
          INTEGER, INTENT(IN) :: MCOL_SOAP
          INTEGER, INTENT(IN) :: NX
          INTEGER, INTENT(IN) :: NY
          REAL(q), INTENT(IN) :: DSOAP1(:,:) !(1:MROW_DSOAP,1:MCOL_DSOAP)
          REAL(q), INTENT(IN) :: SOAP(:,:) !(1:MROW_SOAP,1:MCOL_SOAP)
          REAL(q), INTENT(IN) :: W
! In and output variables
          REAL(q), INTENT(OUT) :: DSOAP2(:,:) !(1:MROW_DSOAP,1:MCOL_DSOAP)
! Local variables
          INTEGER :: INIONS
          INTEGER :: IX,IY
          DO IY=1, NY
             DO IX=1, NX
                IF((IX.NE.1).OR.(INTYP0.EQ.KNTYP0)) THEN
                   INIONS=LNIONS_EST(IX)
                   DSOAP2(IX,IY)=-SOAP(INIONS,IY)*DSOAP1(IX,IY)
                ELSE
                   DSOAP2(IX,IY)=0.0_q
                ENDIF
             ENDDO
          ENDDO
        END SUBROUTINE POLY_DSOAP

      END MODULE SOAP_KERNEL

#ifndef ML_LEARNING_NOT_ACTIVE
!****************************************************************************************************
! Module for regression
!****************************************************************************************************

      MODULE REGRESSION
        USE ML_FF_CONSTANT
        USE LAPACK_DATA
        USE LOGFILE
        USE MPI_DATA
        USE ML_FF_PREC
        IMPLICIT NONE

        CONTAINS

!****************************************************************************************************
! Bayesean linear regression using evidence approximation
! This subroutine is used for many-body fitting
!****************************************************************************************************

      SUBROUTINE BLEA_MB (LAPACKWORK,PAR_SUP_HANDLE,LOGF,CMAT,MROW_CMAT,NROW_CMAT,NCOL_CMAT,LOC_ROW_CMAT,LOC_COL_CMAT, &
                 DMAT_FMAT, &
                 DMAT_YMAT, &
                 FMAT,MROW_FMAT,NROW_FMAT,NCOL_FMAT, &
                 WMAT,NROW_WMAT, &
                 YMAT,MROW_YMAT,NROW_YMAT,LOC_ROW_YMAT,LOC_COL_YMAT, &
                 IFF,IREG,LBHEAD,NTYP,RMEM,SIGV,SIGV0,SIGW,SIGW0,IALGO_LINREG,NSTEP,ML_EPS_REG)
        USE ML_FF_CONSTANT
        USE MPI_DATA
        USE ML_FF_PREC
        USE ML_FF_STRUCT, ONLY: LAPACK_WORK, PARALLEL_SUPER, LOG_PAR, ML_IO_WRITE
        IMPLICIT NONE
! Input variables
        TYPE (LAPACK_WORK)     :: LAPACKWORK
        TYPE (PARALLEL_SUPER)  :: PAR_SUP_HANDLE
        TYPE (LOG_PAR)         :: LOGF
        REAL(q), INTENT(IN)    :: DMAT_FMAT(:) !(1:LOC_ROW_DMAT_FMAT) ! Scaling factor for FMAT
        REAL(q), INTENT(IN)    :: DMAT_YMAT(:) !(1:LOC_ROW_DMAT_YMAT) ! Scaling factor for YMAT
        INTEGER, INTENT(IN)    :: IALGO_LINREG ! Regression method
        INTEGER, INTENT(IN)    :: IREG ! Type of optimization method
        INTEGER, INTENT(IN)    :: LBHEAD(:) !(1:NTYP)
        INTEGER, INTENT(IN)    :: LOC_COL_CMAT
        INTEGER, INTENT(IN)    :: LOC_COL_YMAT
        INTEGER, INTENT(IN)    :: LOC_ROW_CMAT
        INTEGER, INTENT(IN)    :: LOC_ROW_YMAT
        INTEGER, INTENT(IN)    :: MROW_CMAT ! Row size of CMAT
        INTEGER, INTENT(IN)    :: MROW_FMAT ! Row size of FMAT
        INTEGER, INTENT(IN)    :: MROW_YMAT ! Row size of YMAT
        INTEGER, INTENT(IN)    :: NCOL_CMAT ! Column number of CMAT
        INTEGER, INTENT(IN)    :: NROW_CMAT ! Row number of CMAT
        INTEGER, INTENT(IN)    :: NCOL_FMAT ! Column number of FMAT
        INTEGER, INTENT(IN)    :: NROW_FMAT(:) !(1:NTYP) ! Row number of FMAT
        INTEGER, INTENT(IN)    :: NROW_WMAT(:) !(1:NTYP) ! Row number of WMAT
        INTEGER, INTENT(IN)    :: NROW_YMAT ! Row number of YMAT
        INTEGER, INTENT(IN)    :: NTYP ! The number of element
        INTEGER, INTENT(IN)    :: NSTEP ! Current time step
        REAL(q), INTENT(IN)    :: SIGV0,SIGW0 ! (Initial) regularization parameters
        REAL(q), INTENT(IN)    :: ML_EPS_REG ! Convergence parameter for optimization regularization parameters and weights
! Input/output variables
        REAL(q), INTENT(INOUT) :: FMAT(:,:,:) !(1:LOC_ROW_FMAT,1:LOC_COL_FMAT,1:NTYP) ! Design matrix
        REAL(q), INTENT(INOUT) :: SIGV,SIGW ! Precision parameters determined by evidence approximation
        REAL(q), INTENT(INOUT) :: YMAT(:,:) !(1:LOC_ROW_YMAT,1:LOC_COL_YMAT) ! Training data
! Output variables
        REAL(q), INTENT(OUT)   :: RMEM ! Local memory usage
        REAL(q), INTENT(OUT)   :: WMAT(:,:,:) !(1:LOC_ROW_WMAT,1:LOC_COL_WMAT,1:NTYP) ! Regression parameter
! In and output variables
        REAL(q), INTENT(INOUT) :: CMAT(:,:) !(1:LOC_ROW_CMAT,1:LOC_COL_CMAT) ! Covariance matrix
        LOGICAL, INTENT(INOUT) :: IFF ! Logical data determining the execution of FF generation
! Local variables
        REAL(q)                :: ANORM
        REAL(q)                :: CMAT0(1:LOC_ROW_CMAT,1:LOC_COL_CMAT) ! Covariance matrix
        REAL(q), ALLOCATABLE   :: CMAT1(:,:) ! Covariance submatrix
        REAL(q)                :: CRITERIA ! Convergence criteria
#ifndef scaLAPACK
        REAL(q), EXTERNAL      :: DDOT,DNRM2 ! LAPACK routines
#endif
        REAL(q)                :: EF ! Evidence function
        REAL(q)                :: EF_HELP ! Helping variable to compute the evidence function
        REAL(q)                :: EF_OLD ! Old evidence function used to check the convergence
        REAL(q)                :: SIGW_OLD ! Old SIGW to check for convergence
        REAL(q)                :: EIGENVALUE(1:MROW_CMAT) ! Eigenvalue of covariance matrix
        REAL(q)                :: EIGENVECTOR(1:LOC_ROW_CMAT,1:LOC_COL_CMAT) ! Eigen vector of covariance matrix
        REAL(q)                :: EPS1,EPS3,EPS4 ! Several threshold constants
        REAL(q)                :: EPS_REG
        REAL(q), SAVE          :: EPS_REG_LAST=0.0_q
        REAL(q)                :: GAM
        REAL(q), ALLOCATABLE   :: GMAT(:,:) ! G-vector
        INTEGER                :: ICOL
        INTEGER                :: IERR
        INTEGER                :: INFO
        INTEGER                :: INTYP0
        INTEGER                :: IROW
        INTEGER                :: ITE
        INTEGER                :: ITE_RESET
        INTEGER                :: JNTYP0
        LOGICAL                :: LEXIT
        INTEGER                :: LOCVAR_COL
        INTEGER                :: LOCVAR_ROW
        INTEGER                :: NITE
        INTEGER                :: NPROCS
        REAL(q)                :: QUOTIENT_OLD
        REAL(q)                :: RM(1:LOC_ROW_YMAT,1:LOC_COL_YMAT) ! Mean value
        INTEGER                :: TEST_COL
        INTEGER                :: TEST_ROW
        REAL(q)                :: VEC1(1:LOC_ROW_YMAT,1:LOC_COL_YMAT)
        REAL(q), ALLOCATABLE   :: VECONE(:,:)          ! vector filled with 1
        REAL(q), ALLOCATABLE   :: VEC3(:,:)
        REAL(q), ALLOCATABLE   :: WMATALL(:,:)
        REAL(q)                :: WMAX
        REAL(q)                :: WMIN
        REAL(q)                :: ZTR
        REAL(q)                :: ZTR_OLD
        REAL(q)                :: ZTR_MAX
! Variables needed for QR factorization and SVD
        INTEGER                :: RANK
        REAL(q), ALLOCATABLE   :: FMAT_HELP(:,:) ! Design matrix
        REAL(q), ALLOCATABLE   :: FMAT_TRANS_HELP(:,:) ! Design matrix
        REAL(q), ALLOCATABLE   :: YMAT_HELP(:,:) ! Training data
! Variables needed for QR
        INTEGER, ALLOCATABLE   :: JPVT(:)
! Variables needed for SVD
        REAL(q), ALLOCATABLE   :: SINGULAR_VALUES(:)
        REAL(q), ALLOCATABLE   :: SING_U(:,:)
        REAL(q), ALLOCATABLE   :: SING_VT(:,:)
        REAL(q), ALLOCATABLE   :: SIGMA(:,:)
        REAL(q), ALLOCATABLE   :: VSIG(:,:)
        REAL(q), ALLOCATABLE   :: B1_SVD(:,:)
! SET MANUALLY FOR NOW
        LOGICAL                :: LDGELSY ! (only for Lapack) if T use DGELSY for rank-deficient matrices, if F use DGELS
        LOGICAL                :: LDGELSD ! (only for Lapack) if T use DGELSD (driver using SVD), 
                                          ! if F use 'manual' implementation with DGESVD
        LDGELSY = .TRUE.
        LDGELSD = .FALSE.
        PROFILING_START('blea_mb')
! Record memory usage for non-allocated arrays
        RMEM=0.0_q
! Define scalapack descriptors
#ifdef scaLAPACK
        PAR_SUP_HANDLE%DESC_CMAT%M=MROW_FMAT
        PAR_SUP_HANDLE%DESC_CMAT%N=MROW_FMAT
        PAR_SUP_HANDLE%DESC_CMAT%MB=NBLOCK_SCALAPACK
        PAR_SUP_HANDLE%DESC_CMAT%NB=NBLOCK_SCALAPACK
        PAR_SUP_HANDLE%DESC_CMAT%IRSRC=0
        PAR_SUP_HANDLE%DESC_CMAT%ICSRC=0
        CALL SCALAPACK_SET_DESCRIPTOR(PAR_SUP_HANDLE%DESC_CMAT%DESC,PAR_SUP_HANDLE%DESC_CMAT%M,PAR_SUP_HANDLE%DESC_CMAT%N, &
             PAR_SUP_HANDLE%DESC_CMAT%MB,PAR_SUP_HANDLE%DESC_CMAT%NB,PAR_SUP_HANDLE%DESC_CMAT%IRSRC,PAR_SUP_HANDLE%DESC_CMAT%ICSRC, &
             PAR_SUP_HANDLE%CONTEXT_WORLD%ICTXT,PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL,  &
             PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_CMAT%LLD, &
             PAR_SUP_HANDLE%DESC_CMAT%LOC_M,PAR_SUP_HANDLE%DESC_CMAT%LOC_N)
        PAR_SUP_HANDLE%DESC_WMAT%M=MROW_CMAT
        PAR_SUP_HANDLE%DESC_WMAT%N=1
        PAR_SUP_HANDLE%DESC_WMAT%MB=NBLOCK_SCALAPACK
        PAR_SUP_HANDLE%DESC_WMAT%NB=1
        PAR_SUP_HANDLE%DESC_WMAT%IRSRC=0
        PAR_SUP_HANDLE%DESC_WMAT%ICSRC=0
        CALL SCALAPACK_SET_DESCRIPTOR(PAR_SUP_HANDLE%DESC_WMAT%DESC,PAR_SUP_HANDLE%DESC_WMAT%M,PAR_SUP_HANDLE%DESC_WMAT%N, &
             PAR_SUP_HANDLE%DESC_WMAT%MB,PAR_SUP_HANDLE%DESC_WMAT%NB,PAR_SUP_HANDLE%DESC_WMAT%IRSRC,PAR_SUP_HANDLE%DESC_WMAT%ICSRC, &
             PAR_SUP_HANDLE%CONTEXT_WORLD%ICTXT,PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL,  &
             PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_WMAT%LLD, &
             PAR_SUP_HANDLE%DESC_WMAT%LOC_M,PAR_SUP_HANDLE%DESC_WMAT%LOC_N)
        IF (IALGO_LINREG.EQ.2.OR.IALGO_LINREG.EQ.3.OR.IALGO_LINREG.EQ.4) THEN
           PAR_SUP_HANDLE%DESC_FMAT_HELP%M=MROW_CMAT
           PAR_SUP_HANDLE%DESC_FMAT_HELP%N=NCOL_FMAT
           PAR_SUP_HANDLE%DESC_FMAT_HELP%MB=NBLOCK_SCALAPACK
           PAR_SUP_HANDLE%DESC_FMAT_HELP%NB=NBLOCK_SCALAPACK
           PAR_SUP_HANDLE%DESC_FMAT_HELP%IRSRC=0
           PAR_SUP_HANDLE%DESC_FMAT_HELP%ICSRC=0
           CALL SCALAPACK_SET_DESCRIPTOR(PAR_SUP_HANDLE%DESC_FMAT_HELP%DESC,PAR_SUP_HANDLE%DESC_FMAT_HELP%M,PAR_SUP_HANDLE%DESC_FMAT_HELP%N, &
                PAR_SUP_HANDLE%DESC_FMAT_HELP%MB,PAR_SUP_HANDLE%DESC_FMAT_HELP%NB,PAR_SUP_HANDLE%DESC_FMAT_HELP%IRSRC,PAR_SUP_HANDLE%DESC_FMAT_HELP%ICSRC, &
                PAR_SUP_HANDLE%CONTEXT_WORLD%ICTXT,PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL, &
                PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_FMAT_HELP%LLD, &
                PAR_SUP_HANDLE%DESC_FMAT_HELP%LOC_M,PAR_SUP_HANDLE%DESC_FMAT_HELP%LOC_N)
           PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%M=NCOL_FMAT
           PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%N=NROW_CMAT
           PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%MB=NBLOCK_SCALAPACK
           PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%NB=NBLOCK_SCALAPACK
           PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%IRSRC=0
           PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%ICSRC=0
           CALL SCALAPACK_SET_DESCRIPTOR(PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%DESC,PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%M,PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%N, &
                PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%MB,PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%NB,PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%IRSRC,PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%ICSRC, &
                PAR_SUP_HANDLE%CONTEXT_WORLD%ICTXT,PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL, &
                PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%LLD, &
                PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%LOC_M,PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%LOC_N)
        ENDIF 
#endif
! Allocate arrays
#ifdef scaLAPACK
        IF(ALLOCATED(CMAT1)) THEN
           DEALLOCATE(CMAT1)
        ENDIF
        ALLOCATE(CMAT1(1:PAR_SUP_HANDLE%DESC_CMAT%LOC_M,1:PAR_SUP_HANDLE%DESC_CMAT%LOC_N))
        IF(ALLOCATED(GMAT)) THEN
           DEALLOCATE(GMAT)
        ENDIF
        ALLOCATE(GMAT(1:PAR_SUP_HANDLE%DESC_WMAT%LOC_M,1:PAR_SUP_HANDLE%DESC_WMAT%LOC_N))
        IF(ALLOCATED(VECONE)) THEN
           DEALLOCATE(VECONE)
        ENDIF
        ALLOCATE(VECONE(1:PAR_SUP_HANDLE%DESC_WMAT%LOC_M,1:PAR_SUP_HANDLE%DESC_WMAT%LOC_N))
        IF(ALLOCATED(VEC3)) THEN
           DEALLOCATE(VEC3)
        ENDIF
        ALLOCATE(VEC3(1:PAR_SUP_HANDLE%DESC_WMAT%LOC_M,1:PAR_SUP_HANDLE%DESC_WMAT%LOC_N))
        IF(ALLOCATED(WMATALL)) THEN
           DEALLOCATE(WMATALL)
        ENDIF
        ALLOCATE(WMATALL(1:PAR_SUP_HANDLE%DESC_WMAT%LOC_M,1:PAR_SUP_HANDLE%DESC_WMAT%LOC_N))
        IF (IALGO_LINREG.EQ.2.OR.IALGO_LINREG.EQ.3.OR.IALGO_LINREG.EQ.4) THEN
           IF(ALLOCATED(FMAT_HELP)) THEN
              DEALLOCATE(FMAT_HELP)
           ENDIF
           ALLOCATE(FMAT_HELP(1:PAR_SUP_HANDLE%DESC_FMAT_HELP%LOC_M,1:PAR_SUP_HANDLE%DESC_FMAT_HELP%LOC_N))
           IF(ALLOCATED(FMAT_TRANS_HELP)) THEN
              DEALLOCATE(FMAT_TRANS_HELP)
           ENDIF
           ALLOCATE(FMAT_TRANS_HELP(1:PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%LOC_M,1:PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%LOC_N))
           IF(ALLOCATED(YMAT_HELP)) THEN
              DEALLOCATE(YMAT_HELP)
           ENDIF
           ALLOCATE(YMAT_HELP(1:LOC_ROW_YMAT,1:LOC_COL_YMAT))
        ENDIF
#else
        IF(ALLOCATED(CMAT1)) THEN
           DEALLOCATE(CMAT1)
        ENDIF
        ALLOCATE(CMAT1(1:MROW_FMAT,1:MROW_FMAT))
        IF(ALLOCATED(GMAT)) THEN
           DEALLOCATE(GMAT)
        ENDIF
        ALLOCATE(GMAT(1:MROW_CMAT,1:1))
        IF(ALLOCATED(VECONE)) THEN
           DEALLOCATE(VECONE)
        ENDIF
        ALLOCATE(VECONE(1:MROW_CMAT,1:1))
        IF(ALLOCATED(VEC3)) THEN
           DEALLOCATE(VEC3)
        ENDIF
        ALLOCATE(VEC3(1:MROW_CMAT,1:1))
        IF(ALLOCATED(WMATALL)) THEN
           DEALLOCATE(WMATALL)
        ENDIF
        ALLOCATE(WMATALL(1:MROW_CMAT,1:1))
        IF (IALGO_LINREG.EQ.2.OR.IALGO_LINREG.EQ.3.OR.IALGO_LINREG.EQ.4) THEN
           IF(ALLOCATED(FMAT_HELP)) THEN
              DEALLOCATE(FMAT_HELP)
           ENDIF
           ALLOCATE(FMAT_HELP(1:MROW_CMAT,1:NCOL_FMAT))
           IF(ALLOCATED(FMAT_TRANS_HELP)) THEN
              DEALLOCATE(FMAT_TRANS_HELP)
           ENDIF
           ALLOCATE(FMAT_TRANS_HELP(1:NCOL_FMAT,1:NROW_CMAT))
           IF(ALLOCATED(YMAT_HELP)) THEN
              DEALLOCATE(YMAT_HELP)
           ENDIF
           ALLOCATE(YMAT_HELP(1:LOC_ROW_YMAT,1:LOC_COL_YMAT))
        ENDIF
#endif
! Set all helping arrays to 0 here
        WMAT=0.0_q
        WMATALL=0.0_q
        CMAT1=0.0_q
        GMAT=0.0_q
        VECONE=0.0_q
        VEC3=0.0_q
        CMAT=0.0_q
! Constant
        EPS1=1E-14_q
        EPS3=1E0_q
        EPS4=1E-9_q
        CRITERIA =3.0E-2_q
! set number of iterations
        IF(IREG.EQ.1) THEN
           NITE=1
        ELSE IF(IREG.EQ.2) THEN
           NITE=50
        ELSE IF(IREG.EQ.3) THEN
           NITE=50
        ENDIF
        VECONE=1.0_q
! Broadcast SIGW and SIGV to ensure they are the same on all cores
         CALL M_BCAST (PAR_SUP_HANDLE%COMM_WORLD,SIGV,0)
         CALL M_BCAST (PAR_SUP_HANDLE%COMM_WORLD,SIGW,0)
         CALL M_BCAST (PAR_SUP_HANDLE%COMM_WORLD,SIGV0,0)
         CALL M_BCAST (PAR_SUP_HANDLE%COMM_WORLD,SIGW0,0)
! Set the initial values for SIGV and SIGW
        IF(ABS(SIGW) .LT. EPS1) THEN
           SIGV=SIGV0
           SIGW=SIGW0
        ELSE IF(ABS(SIGV/SIGW) .LT. EPS3) THEN
           SIGV=SIGV0
           SIGW=SIGW0
        ENDIF
! Scale design matrix
        PROFILING_START('scale_blea_mb')
        DO INTYP0=1, NTYP
           DO ICOL=1, NCOL_FMAT
#ifdef scaLAPACK
              TEST_COL=MOD((ICOL-1)/PAR_SUP_HANDLE%DESC_FMAT_FFM%NB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL)
              IF(PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL.NE.TEST_COL) THEN
                 CYCLE
              ENDIF
              LOCVAR_COL=GLOBAL2LOCAL(ICOL,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_FMAT_FFM%NB)
#else
              LOCVAR_COL=ICOL
#endif
              DO IROW=1, NROW_FMAT(INTYP0)
#ifdef scaLAPACK
                 TEST_ROW=MOD((IROW-1)/PAR_SUP_HANDLE%DESC_FMAT_FFM%MB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW)
                 IF(PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW.NE.TEST_ROW) THEN
                    CYCLE
                 ENDIF
                 LOCVAR_ROW=GLOBAL2LOCAL(IROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%DESC_FMAT_FFM%MB)
#else
                 LOCVAR_ROW=IROW
#endif
                 FMAT(LOCVAR_ROW,LOCVAR_COL,INTYP0)=FMAT(LOCVAR_ROW,LOCVAR_COL,INTYP0)*DMAT_FMAT(LOCVAR_COL)
              ENDDO
           ENDDO
        ENDDO
        DO IROW=1, NROW_YMAT
#ifdef scaLAPACK
           TEST_ROW=MOD((IROW-1)/PAR_SUP_HANDLE%DESC_YMAT_FFM%MB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW)
           IF(PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW.NE.TEST_ROW) THEN
              CYCLE
           ENDIF
           LOCVAR_ROW=GLOBAL2LOCAL(IROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%DESC_FMAT_FFM%MB)
           LOCVAR_COL=GLOBAL2LOCAL(   1,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_FMAT_FFM%NB)
#else
           LOCVAR_ROW=IROW
           LOCVAR_COL=1
#endif
           YMAT(LOCVAR_ROW,LOCVAR_COL)=YMAT(LOCVAR_ROW,LOCVAR_COL)*DMAT_YMAT(LOCVAR_ROW)
        ENDDO
        PROFILING_STOP('scale_blea_mb')
!----------------------------------------------------------------------------------------------------
! Standard regression 
!----------------------------------------------------------------------------------------------------
        IF (IALGO_LINREG.EQ.2) THEN
           PROFILING_START('qr_blea_mb')
! Calculate regression parameters using QR factorization to solve the
! least squares problem ||Y-Phi w||=min;
! (P)DGELS* overwrites input matrices, so copy FMAT and YMAT into aux mat;
! To use DGELSY, we need to take the transpose of FMAT; (P)DGELS has the
! additional argument TRANS, but to make the Lapack and Scalapack implementations 
! similar, we don't use it and we take the transpose of FMAT in both cases
! First copy YMAT 
           YMAT_HELP=YMAT
! Then copy FMAT
           DO INTYP0=1, NTYP
              DO ICOL=1, NCOL_FMAT
#ifdef scaLAPACK
                 CALL PDCOPY (NROW_FMAT(INTYP0),FMAT(1,1,INTYP0),1,ICOL,PAR_SUP_HANDLE%DESC_FMAT_FFM%DESC,1, &
                      FMAT_HELP(1,1),1+LBHEAD(INTYP0),ICOL,PAR_SUP_HANDLE%DESC_FMAT_HELP%DESC,1)
#else
                 CALL DCOPY (NROW_FMAT(INTYP0),FMAT(1,ICOL,INTYP0),1,FMAT_HELP(1+LBHEAD(INTYP0),ICOL),1)
#endif
              ENDDO
           ENDDO
#ifdef scaLAPACK
! Transpose FMAT_HELP 
           CALL PDGEADD('T',NCOL_FMAT,NROW_CMAT,1.0_q,FMAT_HELP(1,1),1,1,PAR_SUP_HANDLE%DESC_FMAT_HELP%DESC, &
                        0.0_q,FMAT_TRANS_HELP(1,1),1,1,PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%DESC,IERR)
! First calculate the size of the work array, then calculate the regression coefficients
           CALL SCALAPACK_WORK_ARRAY_SIZE_PDGELS_INIT(PAR_SUP_HANDLE%SCAWORK)
           CALL PDGELS('N',NCOL_FMAT,NROW_CMAT,1,FMAT_TRANS_HELP(1,1),1,1,PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%DESC, &
                YMAT_HELP(1,1),1,1,PAR_SUP_HANDLE%DESC_YMAT_FFM%DESC, &
                PAR_SUP_HANDLE%SCAWORK%WORK(1),PAR_SUP_HANDLE%SCAWORK%LWORK,PAR_SUP_HANDLE%SCAWORK%INFO)
           CALL SCALAPACK_WORK_ARRAY_SIZE_PDGELS_FIN(PAR_SUP_HANDLE%SCAWORK)
           CALL SCALAPACK_WORK_ARRAY_PDGELS_INIT(PAR_SUP_HANDLE%SCAWORK)
           CALL PDGELS('N',NCOL_FMAT,NROW_CMAT,1,FMAT_TRANS_HELP(1,1),1,1,PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%DESC, &
                YMAT_HELP(1,1),1,1,PAR_SUP_HANDLE%DESC_YMAT_FFM%DESC, &
                PAR_SUP_HANDLE%SCAWORK%WORK(1),PAR_SUP_HANDLE%SCAWORK%LWORK,PAR_SUP_HANDLE%SCAWORK%INFO)
           CALL SCALAPACK_WORK_ARRAY_PDGELS_FIN(PAR_SUP_HANDLE%SCAWORK)
#else
           DO IROW=1, NROW_CMAT
              DO ICOL=1, NCOL_FMAT
                 FMAT_TRANS_HELP(ICOL,IROW) = FMAT_HELP(IROW,ICOL)
              ENDDO
           ENDDO
           IF (.NOT.LDGELSY) THEN
             CALL LAPACK_WORK_ARRAY_SIZE_DGELSY_INIT(LAPACKWORK)
             CALL DGELS('N',NCOL_FMAT,NROW_CMAT,1,FMAT_TRANS_HELP(1,1),NCOL_FMAT,YMAT_HELP(1,1),MROW_YMAT, &
                  LAPACKWORK%WORK,LAPACKWORK%LWORK,LAPACKWORK%INFO)
             CALL LAPACK_WORK_ARRAY_SIZE_DGELSY_FIN(LAPACKWORK)
             CALL LAPACK_WORK_ARRAY_DGELSY_INIT(LAPACKWORK)
             CALL DGELS('N',NCOL_FMAT,NROW_CMAT,1,FMAT_TRANS_HELP(1,1),NCOL_FMAT,YMAT_HELP(1,1),MROW_YMAT, &
                  LAPACKWORK%WORK,LAPACKWORK%LWORK,LAPACKWORK%INFO)
             CALL LAPACK_WORK_ARRAY_DGELSY_FIN(LAPACKWORK)
           ELSE
             IF (ALLOCATED(JPVT)) THEN
                DEALLOCATE(JPVT)
             ENDIF
             ALLOCATE(JPVT(NROW_CMAT))
             JPVT=0
             CALL LAPACK_WORK_ARRAY_SIZE_DGELSY_INIT(LAPACKWORK)
             CALL DGELSY (NCOL_FMAT,NROW_CMAT,1,FMAT_TRANS_HELP(1,1),NCOL_FMAT,YMAT_HELP(1,1),MROW_YMAT, &
                  JPVT,ML_EPS_REG,RANK,LAPACKWORK%WORK,LAPACKWORK%LWORK,LAPACKWORK%INFO)
             CALL LAPACK_WORK_ARRAY_SIZE_DGELSY_FIN(LAPACKWORK)
             CALL LAPACK_WORK_ARRAY_DGELSY_INIT(LAPACKWORK)
             CALL DGELSY (NCOL_FMAT,NROW_CMAT,1,FMAT_TRANS_HELP(1,1),NCOL_FMAT,YMAT_HELP(1,1),MROW_YMAT, &
                  JPVT,ML_EPS_REG,RANK,LAPACKWORK%WORK,LAPACKWORK%LWORK,LAPACKWORK%INFO)
             CALL LAPACK_WORK_ARRAY_DGELSY_FIN(LAPACKWORK)
           ENDIF
#endif
! Put results into WMAT
           DO INTYP0=1, NTYP
#ifdef scaLAPACK
              CALL PDCOPY(NROW_WMAT(INTYP0),YMAT_HELP(1,1),1+LBHEAD(INTYP0),1,PAR_SUP_HANDLE%DESC_WMAT%DESC, &
                   1,WMAT(1,1,INTYP0),1,1,PAR_SUP_HANDLE%DESC_WMAT_FFM%DESC,1)
#else
              CALL DCOPY (NROW_WMAT(INTYP0),YMAT_HELP(1+LBHEAD(INTYP0),1),1,WMAT(1,1,INTYP0),1)
#endif
           ENDDO
           ! Write regression results here
           LEXIT=.TRUE.
           CALL LOGLINE_REGR(LOGF,                   &
                             NSTEP,                  &
                             IALGO_LINREG,           &
                             1,                      &
                             0.0_q,                  &
                             0.0_q,                  &
                             0.0_q,                  &
                             0.0_q,                  &
                             0.0_q,                  &
                             0.0_q,                  &
                             LEXIT)
           
           IF(ALLOCATED(FMAT_HELP)) DEALLOCATE(FMAT_HELP)
           IF(ALLOCATED(FMAT_TRANS_HELP)) DEALLOCATE(FMAT_TRANS_HELP)
           IF(ALLOCATED(YMAT_HELP)) DEALLOCATE(YMAT_HELP)
           PROFILING_STOP('qr_blea_mb')
        ELSE IF (IALGO_LINREG.EQ.3.OR.IALGO_LINREG.EQ.4) THEN
           PROFILING_START('svd_blea_mb')
! Calculate regression parameters using SVD to solve the
! least squares problem ||Y-Phi w||=min;
! (P)DGESVD, DGELSD overwrite or destroy input matrices, so copy FMAT and YMAT into aux mat;
! To use (P)DGESVD, DGELSD we need to take the transpose of FMAT
! IALGO_LINREG=3 uses truncated SVD
! IALGO_LINREG=4 uses SVD plus Tikhonov regularization
! First copy YMAT
           YMAT_HELP=YMAT
! Then copy FMAT to FMAT_HELP
           DO INTYP0=1, NTYP
              DO ICOL=1, NCOL_FMAT
#ifdef scaLAPACK
                 CALL PDCOPY (NROW_FMAT(INTYP0),FMAT(1,1,INTYP0),1,ICOL,PAR_SUP_HANDLE%DESC_FMAT_FFM%DESC,1, &
                      FMAT_HELP(1,1),1+LBHEAD(INTYP0),ICOL,PAR_SUP_HANDLE%DESC_FMAT_HELP%DESC,1)
#else
                 CALL DCOPY (NROW_FMAT(INTYP0),FMAT(1,ICOL,INTYP0),1,FMAT_HELP(1+LBHEAD(INTYP0),ICOL),1)
#endif
              ENDDO
           ENDDO
! Allocate some helping arrays
           ALLOCATE(SINGULAR_VALUES(1:NROW_CMAT))
#ifdef scaLAPACK
! Transpose FMAT_HELP 
           CALL PDGEADD('T',NCOL_FMAT,NROW_CMAT,1.0_q,FMAT_HELP(1,1),1,1,PAR_SUP_HANDLE%DESC_FMAT_HELP%DESC, &
                        0.0_q,FMAT_TRANS_HELP(1,1),1,1,PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%DESC,IERR)
! Set descriptors for SING_U, SING_VT for SVD
           PAR_SUP_HANDLE%DESC_SING_U%M=NCOL_FMAT
           PAR_SUP_HANDLE%DESC_SING_U%N=NROW_CMAT !reduced SVD representation: U is rectangular matrix
           PAR_SUP_HANDLE%DESC_SING_U%MB=NBLOCK_SCALAPACK
           PAR_SUP_HANDLE%DESC_SING_U%NB=NBLOCK_SCALAPACK
           PAR_SUP_HANDLE%DESC_SING_U%IRSRC=0
           PAR_SUP_HANDLE%DESC_SING_U%ICSRC=0
           CALL SCALAPACK_SET_DESCRIPTOR(PAR_SUP_HANDLE%DESC_SING_U%DESC,PAR_SUP_HANDLE%DESC_SING_U%M,PAR_SUP_HANDLE%DESC_SING_U%N, &
              PAR_SUP_HANDLE%DESC_SING_U%MB,PAR_SUP_HANDLE%DESC_SING_U%NB,PAR_SUP_HANDLE%DESC_SING_U%IRSRC,PAR_SUP_HANDLE%DESC_SING_U%ICSRC, &
              PAR_SUP_HANDLE%CONTEXT_WORLD%ICTXT,PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL,  &
              PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_SING_U%LLD, &
              PAR_SUP_HANDLE%DESC_SING_U%LOC_M,PAR_SUP_HANDLE%DESC_SING_U%LOC_N)
           PAR_SUP_HANDLE%DESC_SING_VT%M=NROW_CMAT
           PAR_SUP_HANDLE%DESC_SING_VT%N=NROW_CMAT
           PAR_SUP_HANDLE%DESC_SING_VT%MB=NBLOCK_SCALAPACK
           PAR_SUP_HANDLE%DESC_SING_VT%NB=NBLOCK_SCALAPACK
           PAR_SUP_HANDLE%DESC_SING_VT%IRSRC=0
           PAR_SUP_HANDLE%DESC_SING_VT%ICSRC=0
           CALL SCALAPACK_SET_DESCRIPTOR(PAR_SUP_HANDLE%DESC_SING_VT%DESC,PAR_SUP_HANDLE%DESC_SING_VT%M,PAR_SUP_HANDLE%DESC_SING_VT%N, &
              PAR_SUP_HANDLE%DESC_SING_VT%MB,PAR_SUP_HANDLE%DESC_SING_VT%NB,PAR_SUP_HANDLE%DESC_SING_VT%IRSRC,PAR_SUP_HANDLE%DESC_SING_VT%ICSRC, &
              PAR_SUP_HANDLE%CONTEXT_WORLD%ICTXT,PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL,  &
              PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_SING_VT%LLD, &
              PAR_SUP_HANDLE%DESC_SING_VT%LOC_M,PAR_SUP_HANDLE%DESC_SING_VT%LOC_N)
! Deallocate FMAT_HELP before allocating U and VT
           DEALLOCATE(FMAT_HELP)
! Allocate U and VT for SVD
           ALLOCATE(SING_U(1:PAR_SUP_HANDLE%DESC_SING_U%LOC_M,1:PAR_SUP_HANDLE%DESC_SING_U%LOC_N))
           ALLOCATE(SING_VT(1:PAR_SUP_HANDLE%DESC_SING_VT%LOC_M,1:PAR_SUP_HANDLE%DESC_SING_VT%LOC_N))
! Then do the singular value decomposition
           CALL SCALAPACK_WORK_ARRAY_SIZE_PDGESVD_INIT(PAR_SUP_HANDLE%SCAWORK)
           CALL PDGESVD('V','V',NCOL_FMAT,NROW_CMAT,FMAT_TRANS_HELP(1,1),1,1,PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%DESC, &
                      SINGULAR_VALUES(1),SING_U,1,1,PAR_SUP_HANDLE%DESC_SING_U%DESC,SING_VT,1,1 , &
                      PAR_SUP_HANDLE%DESC_SING_VT%DESC,PAR_SUP_HANDLE%SCAWORK%WORK(1), &
                      PAR_SUP_HANDLE%SCAWORK%LWORK,PAR_SUP_HANDLE%SCAWORK%INFO)
           CALL SCALAPACK_WORK_ARRAY_SIZE_PDGESVD_FIN(PAR_SUP_HANDLE%SCAWORK)
           CALL SCALAPACK_WORK_ARRAY_PDGESVD_INIT(PAR_SUP_HANDLE%SCAWORK)
           CALL PDGESVD('V','V',NCOL_FMAT,NROW_CMAT,FMAT_TRANS_HELP(1,1),1,1,PAR_SUP_HANDLE%DESC_FMAT_TRANS_HELP%DESC, &
                      SINGULAR_VALUES(1),SING_U,1,1,PAR_SUP_HANDLE%DESC_SING_U%DESC,SING_VT,1,1 , &
                      PAR_SUP_HANDLE%DESC_SING_VT%DESC,PAR_SUP_HANDLE%SCAWORK%WORK(1), &
                      PAR_SUP_HANDLE%SCAWORK%LWORK,PAR_SUP_HANDLE%SCAWORK%INFO)
           CALL SCALAPACK_WORK_ARRAY_PDGESVD_FIN(PAR_SUP_HANDLE%SCAWORK)
! Condition number of the design matrix
           RANK=NROW_CMAT
! Only do the following step for truncated SVD
           IF (IALGO_LINREG.EQ.3) THEN 
! Determine the effective rank (SINGULAR_VALUES is a global array)
              DO IROW=1,NROW_CMAT
                 IF (SINGULAR_VALUES(IROW)/SINGULAR_VALUES(1).LT.ML_EPS_REG) THEN
                    RANK=IROW-1
                    EXIT
                 ENDIF
              ENDDO
           ENDIF
! Make descriptors for SIGMA and allocate SIGMA 
           PAR_SUP_HANDLE%DESC_SING_SIGMA%M=RANK
           PAR_SUP_HANDLE%DESC_SING_SIGMA%N=RANK
           PAR_SUP_HANDLE%DESC_SING_SIGMA%MB=NBLOCK_SCALAPACK
           PAR_SUP_HANDLE%DESC_SING_SIGMA%NB=NBLOCK_SCALAPACK
           PAR_SUP_HANDLE%DESC_SING_SIGMA%IRSRC=0
           PAR_SUP_HANDLE%DESC_SING_SIGMA%ICSRC=0
           CALL SCALAPACK_SET_DESCRIPTOR(PAR_SUP_HANDLE%DESC_SING_SIGMA%DESC,PAR_SUP_HANDLE%DESC_SING_SIGMA%M,PAR_SUP_HANDLE%DESC_SING_SIGMA%N, &
              PAR_SUP_HANDLE%DESC_SING_SIGMA%MB,PAR_SUP_HANDLE%DESC_SING_SIGMA%NB,PAR_SUP_HANDLE%DESC_SING_SIGMA%IRSRC,PAR_SUP_HANDLE%DESC_SING_SIGMA%ICSRC, &
              PAR_SUP_HANDLE%CONTEXT_WORLD%ICTXT,PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL, &
              PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_SING_SIGMA%LLD, &
              PAR_SUP_HANDLE%DESC_SING_SIGMA%LOC_M,PAR_SUP_HANDLE%DESC_SING_SIGMA%LOC_N)
           ALLOCATE(SIGMA(PAR_SUP_HANDLE%DESC_SING_SIGMA%LOC_M,PAR_SUP_HANDLE%DESC_SING_SIGMA%LOC_N))
           SIGMA=0.0_q
! Calculate SIGMA matrix (diagonal)
           DO IROW=1, RANK
              TEST_ROW=MOD((IROW-1)/PAR_SUP_HANDLE%DESC_SING_SIGMA%MB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW)
              TEST_COL=MOD((IROW-1)/PAR_SUP_HANDLE%DESC_SING_SIGMA%NB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL)
              IF(PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW.EQ.TEST_ROW.AND.PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL.EQ.TEST_COL) THEN
                 LOCVAR_ROW=GLOBAL2LOCAL(IROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%DESC_SING_SIGMA%MB)
                 LOCVAR_COL=GLOBAL2LOCAL(IROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_SING_SIGMA%NB)
! For truncated SVD we only take the inverse eigenvalues, for SVD plus Tikhonov we ad regularization parameters
                 IF (IALGO_LINREG.EQ.3) THEN
                    SIGMA(LOCVAR_ROW,LOCVAR_COL)=1.0_q/SINGULAR_VALUES(IROW)
                 ELSE IF (IALGO_LINREG.EQ.4) THEN
                    SIGMA(LOCVAR_ROW,LOCVAR_COL)=SINGULAR_VALUES(IROW)/(SINGULAR_VALUES(IROW)**2+SIGW0/SIGV0)
                 ENDIF
              ENDIF
           ENDDO
! Make descriptors for VSIG and allocate VSIG
           PAR_SUP_HANDLE%DESC_SING_VSIG%M=RANK
           PAR_SUP_HANDLE%DESC_SING_VSIG%N=RANK
           PAR_SUP_HANDLE%DESC_SING_VSIG%MB=NBLOCK_SCALAPACK
           PAR_SUP_HANDLE%DESC_SING_VSIG%NB=NBLOCK_SCALAPACK
           PAR_SUP_HANDLE%DESC_SING_VSIG%IRSRC=0
           PAR_SUP_HANDLE%DESC_SING_VSIG%ICSRC=0
           CALL SCALAPACK_SET_DESCRIPTOR(PAR_SUP_HANDLE%DESC_SING_VSIG%DESC,PAR_SUP_HANDLE%DESC_SING_VSIG%M,PAR_SUP_HANDLE%DESC_SING_VSIG%N, &
              PAR_SUP_HANDLE%DESC_SING_VSIG%MB,PAR_SUP_HANDLE%DESC_SING_VSIG%NB,PAR_SUP_HANDLE%DESC_SING_VSIG%IRSRC,PAR_SUP_HANDLE%DESC_SING_VSIG%ICSRC, &
              PAR_SUP_HANDLE%CONTEXT_WORLD%ICTXT,PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL, &
              PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_SING_VSIG%LLD, &
              PAR_SUP_HANDLE%DESC_SING_VSIG%LOC_M,PAR_SUP_HANDLE%DESC_SING_VSIG%LOC_N)
           ALLOCATE(VSIG(PAR_SUP_HANDLE%DESC_SING_VSIG%LOC_M,PAR_SUP_HANDLE%DESC_SING_VSIG%LOC_N))
! Calculate V times SiGMA (NB: PDGESVD returns VT)
           CALL PDGEMM('T','N',RANK,RANK,RANK,1.0_q, &
                      SING_VT,1,1,PAR_SUP_HANDLE%DESC_SING_VT%DESC, &
                      SIGMA,1,1,PAR_SUP_HANDLE%DESC_SING_SIGMA%DESC, &
                      0.0_q,VSIG,1,1,PAR_SUP_HANDLE%DESC_SING_VSIG%DESC)
! Make descriptors for b1=UT*b, where b is YMAT and allocate b1
           PAR_SUP_HANDLE%DESC_SING_B1_SVD%M=RANK
           PAR_SUP_HANDLE%DESC_SING_B1_SVD%N=1
           PAR_SUP_HANDLE%DESC_SING_B1_SVD%MB=NBLOCK_SCALAPACK
           PAR_SUP_HANDLE%DESC_SING_B1_SVD%NB=NBLOCK_SCALAPACK
           PAR_SUP_HANDLE%DESC_SING_B1_SVD%IRSRC=0
           PAR_SUP_HANDLE%DESC_SING_B1_SVD%ICSRC=0
           CALL SCALAPACK_SET_DESCRIPTOR(PAR_SUP_HANDLE%DESC_SING_B1_SVD%DESC,PAR_SUP_HANDLE%DESC_SING_B1_SVD%M,PAR_SUP_HANDLE%DESC_SING_B1_SVD%N, &
              PAR_SUP_HANDLE%DESC_SING_B1_SVD%MB,PAR_SUP_HANDLE%DESC_SING_B1_SVD%NB,PAR_SUP_HANDLE%DESC_SING_B1_SVD%IRSRC,PAR_SUP_HANDLE%DESC_SING_B1_SVD%ICSRC, &
              PAR_SUP_HANDLE%CONTEXT_WORLD%ICTXT,PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL, &
              PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_SING_B1_SVD%LLD, &
              PAR_SUP_HANDLE%DESC_SING_B1_SVD%LOC_M,PAR_SUP_HANDLE%DESC_SING_B1_SVD%LOC_N)
           ALLOCATE(B1_SVD(PAR_SUP_HANDLE%DESC_SING_B1_SVD%LOC_M,PAR_SUP_HANDLE%DESC_SING_B1_SVD%LOC_N))
! Calculate B1_SVD
           CALL PDGEMV('T',NCOL_FMAT,RANK,1.0_q,SING_U,1,1,PAR_SUP_HANDLE%DESC_SING_U%DESC, &
                       YMAT_HELP,1,1,PAR_SUP_HANDLE%DESC_YMAT_FFM%DESC,1,0.0_q, &
                       B1_SVD,1,1,PAR_SUP_HANDLE%DESC_SING_B1_SVD%DESC,1) 
! Calculate regression coefficients (here we reuse YMAT_HELP)
           YMAT_HELP=0.0_q
           CALL PDGEMV('N',RANK,RANK,1.0_q,VSIG,1,1,PAR_SUP_HANDLE%DESC_SING_VSIG%DESC, &
                       B1_SVD,1,1,PAR_SUP_HANDLE%DESC_SING_B1_SVD%DESC,1,0.0_q, &
                       YMAT_HELP,1,1,PAR_SUP_HANDLE%DESC_YMAT_FFM%DESC,1) 
#else
           DO IROW=1, NROW_CMAT
              DO ICOL=1, NCOL_FMAT
                 FMAT_TRANS_HELP(ICOL,IROW) = FMAT_HELP(IROW,ICOL)
              ENDDO
           ENDDO
           IF (.NOT.LDGELSD) THEN
              ALLOCATE(SING_U(1:NCOL_FMAT,1:NROW_CMAT)) !reduced SVD representation: U is rectangular matrix
              ALLOCATE(SING_VT(1:NROW_CMAT,1:NROW_CMAT))
! First do the singular value decomposition
              CALL LAPACK_WORK_ARRAY_SIZE_DGESVD_INIT(LAPACKWORK)
              CALL DGESVD('S','A',NCOL_FMAT,NROW_CMAT,FMAT_TRANS_HELP(1,1),NCOL_FMAT, &
                          SINGULAR_VALUES,SING_U,NCOL_FMAT,SING_VT,NROW_CMAT, &
                          LAPACKWORK%WORK(1),LAPACKWORK%LWORK,LAPACKWORK%INFO)
              CALL LAPACK_WORK_ARRAY_SIZE_DGESVD_FIN(LAPACKWORK)
              CALL LAPACK_WORK_ARRAY_DGESVD_INIT(LAPACKWORK)
              CALL DGESVD('S','A',NCOL_FMAT,NROW_CMAT,FMAT_TRANS_HELP(1,1),NCOL_FMAT, &
                          SINGULAR_VALUES,SING_U,NCOL_FMAT,SING_VT,NROW_CMAT, &
                          LAPACKWORK%WORK(1),LAPACKWORK%LWORK,LAPACKWORK%INFO)
              CALL LAPACK_WORK_ARRAY_DGESVD_FIN(LAPACKWORK)
              RANK=NROW_CMAT
! Do this only for truncated SVD
              IF (IALGO_LINREG.EQ.3)  THEN ! Truncated SVD
! Determine the effective rank
                DO IROW=1,NROW_CMAT
                   IF (SINGULAR_VALUES(IROW)/SINGULAR_VALUES(1).LT.ML_EPS_REG) THEN
                      RANK=IROW-1
                      EXIT
                   ENDIF
                ENDDO
              ENDIF
! Construct inverse of Sigma
              ALLOCATE(SIGMA(1:RANK,1:RANK))
              SIGMA=0.0d0
              DO IROW=1, RANK
! For truncated SVD we only take the inverse eigenvalues, for SVD plus Tikhonov we ad regularization parameters
                 IF (IALGO_LINREG.EQ.3) THEN
                    SIGMA(IROW,IROW) = 1.0_q/SINGULAR_VALUES(IROW)
                 ELSE IF (IALGO_LINREG.EQ.4) THEN
                    SIGMA(IROW,IROW) = SINGULAR_VALUES(IROW)/(SINGULAR_VALUES(IROW)**2+SIGW0/SIGV0)
                 ENDIF
              ENDDO
! Calculate product of V and Sigma (SING_VT is the transpose of V)
              ALLOCATE(VSIG(1:RANK,1:RANK))
              CALL DGEMM('T','N',RANK,RANK,RANK,1.0_q,SING_VT,NROW_CMAT,SIGMA, &
                         RANK,0.0_q,VSIG,RANK)
! Calculate b1=U^T*b, where b1 contains the first k entries and b is YMAT_HELP
              ALLOCATE(B1_SVD(1:RANK,1:1))
              CALL DGEMV('T',NCOL_FMAT,RANK,1.0_q,SING_U,NCOL_FMAT,YMAT_HELP,1,0.0_q,B1_SVD,1)
! Calculate solutions x, here we reuse YMAT_HELP
              YMAT_HELP=0.0_q
              CALL DGEMV('N',RANK,RANK,1.0_q,VSIG,RANK,B1_SVD,1,0.0_q,YMAT_HELP,1)
           ELSE 
              CALL LAPACK_WORK_ARRAY_SIZE_DGELSD_INIT(LAPACKWORK)
              CALL DGELSD (NCOL_FMAT,NROW_CMAT,1,FMAT_TRANS_HELP(1,1),NCOL_FMAT,YMAT_HELP(1,1),NROW_YMAT, &
                   SINGULAR_VALUES,ML_EPS_REG,RANK,LAPACKWORK%WORK,LAPACKWORK%LWORK,LAPACKWORK%IPIV,LAPACKWORK%INFO)
              CALL LAPACK_WORK_ARRAY_SIZE_DGELSD_FIN(LAPACKWORK)
              CALL LAPACK_WORK_ARRAY_DGELSD_INIT(LAPACKWORK)
              CALL DGELSD (NCOL_FMAT,NROW_CMAT,1,FMAT_TRANS_HELP(1,1),NCOL_FMAT,YMAT_HELP(1,1),NROW_YMAT, &
                   SINGULAR_VALUES,ML_EPS_REG,RANK,LAPACKWORK%WORK,LAPACKWORK%LWORK,LAPACKWORK%IPIV,LAPACKWORK%INFO)
              CALL LAPACK_WORK_ARRAY_DGELSD_FIN(LAPACKWORK)
           ENDIF
#endif
! Put results into WMAT
           DO INTYP0=1, NTYP
#ifdef scaLAPACK
              CALL PDCOPY(NROW_WMAT(INTYP0),YMAT_HELP(1,1),1+LBHEAD(INTYP0),1,PAR_SUP_HANDLE%DESC_WMAT%DESC, &
                   1,WMAT(1,1,INTYP0),1,1,PAR_SUP_HANDLE%DESC_WMAT_FFM%DESC,1)
#else
              CALL DCOPY (NROW_WMAT(INTYP0),YMAT_HELP(1+LBHEAD(INTYP0),1),1,WMAT(1,1,INTYP0),1)
#endif
           ENDDO
           ! Determine smallest and largest singular value
           WMIN=SINGULAR_VALUES(1)
           WMAX=SINGULAR_VALUES(1)
           DO IROW=1, MIN(NCOL_FMAT,NROW_CMAT)
              WMIN=MIN(WMIN,SINGULAR_VALUES(IROW))
              WMAX=MAX(WMAX,SINGULAR_VALUES(IROW))
           ENDDO
           ! Write regression results here
           LEXIT=.TRUE.
           CALL LOGLINE_REGR(LOGF,                   &
                             NSTEP,                  &
                             IALGO_LINREG,           &
                             1,                      &
                             0.0_q,                  &
                             0.0_q,                  &
                             0.0_q,                  &
                             WMAX / WMIN,            &
                             WMAX,                   &
                             0.0_q,                  &
                             LEXIT)
! Deallocate the stuff used here
           IF(ALLOCATED(FMAT_HELP)) DEALLOCATE(FMAT_HELP)
           IF(ALLOCATED(FMAT_TRANS_HELP)) DEALLOCATE(FMAT_TRANS_HELP)
           IF(ALLOCATED(SINGULAR_VALUES)) DEALLOCATE(SINGULAR_VALUES)
           IF(ALLOCATED(YMAT_HELP)) DEALLOCATE(YMAT_HELP)
           IF(ALLOCATED(B1_SVD)) DEALLOCATE(B1_SVD)
           IF(ALLOCATED(SIGMA)) DEALLOCATE(SIGMA)
           IF(ALLOCATED(SING_U)) DEALLOCATE(SING_U)
           IF(ALLOCATED(SING_VT)) DEALLOCATE(SING_VT)
           IF(ALLOCATED(VSIG)) DEALLOCATE(VSIG)
           PROFILING_STOP('svd_blea_mb')
        ELSE
! Bayesian linear regression
           PROFILING_START('init_cmat_blea_mb')
           DO INTYP0=1, NTYP
              DO JNTYP0=1, NTYP
#ifdef scaLAPACK
                 CALL PDGEMM('N','T',NROW_FMAT(INTYP0),NROW_FMAT(JNTYP0),NCOL_FMAT,1.0_q, &
                      FMAT(1,1,INTYP0),1,1,PAR_SUP_HANDLE%DESC_FMAT_FFM%DESC, &
                      FMAT(1,1,JNTYP0),1,1,PAR_SUP_HANDLE%DESC_FMAT_FFM%DESC, &
                      0.0_q,CMAT1(1,1),1,1,PAR_SUP_HANDLE%DESC_CMAT%DESC)
                 DO ICOL=1, NROW_FMAT(JNTYP0)
                    CALL PDCOPY(NROW_FMAT(INTYP0),CMAT1(1,1),1,ICOL,PAR_SUP_HANDLE%DESC_CMAT%DESC,1,CMAT0(1,1),1+LBHEAD(INTYP0),ICOL+LBHEAD(JNTYP0),PAR_SUP_HANDLE%DESC_CMAT_FFM%DESC,1)
                 ENDDO
#else
                 CALL DGEMM ('N','T',NROW_FMAT(INTYP0),NROW_FMAT(JNTYP0),NCOL_FMAT,1.0_q,FMAT(1,1,INTYP0),MROW_FMAT,FMAT(1,1,JNTYP0),MROW_FMAT,0.0_q,CMAT1(1,1),MROW_FMAT)
                 DO ICOL=1, NROW_FMAT(JNTYP0)
                    CALL DCOPY (NROW_FMAT(INTYP0),CMAT1(1,ICOL),1,CMAT0(1+LBHEAD(INTYP0),ICOL+LBHEAD(JNTYP0)),1)
                 ENDDO
#endif
              ENDDO
           ENDDO
           PROFILING_STOP('init_cmat_blea_mb')
           PROFILING_START('init_gmat_blea_mb')
           GMAT=0.0_q
#ifdef scaLAPACK
           DO INTYP0=1, NTYP
              CALL PDGEMV ('N',NROW_FMAT(INTYP0),NCOL_FMAT,1.0_q,FMAT(1,1,INTYP0),1,1,PAR_SUP_HANDLE%DESC_FMAT_FFM%DESC, &
                   YMAT(1,1),1,1,PAR_SUP_HANDLE%DESC_YMAT_FFM%DESC,1, &
                   0.0_q,GMAT(1,1),1+LBHEAD(INTYP0),1,PAR_SUP_HANDLE%DESC_WMAT%DESC,1)
           ENDDO
#else
           DO INTYP0=1, NTYP
              CALL DGEMV ('N',NROW_FMAT(INTYP0),NCOL_FMAT,1.0_q,FMAT(1,1,INTYP0),MROW_FMAT,YMAT(1,1),1,0_q,GMAT(1+LBHEAD(INTYP0),1),1)
           ENDDO
#endif
           PROFILING_STOP('init_gmat_blea_mb')
! Calculate eigenvalues of covariance matrix
           DO IROW=1,LOC_ROW_CMAT
              DO ICOL=1,LOC_COL_CMAT
                 CMAT(IROW,ICOL)=CMAT0(IROW,ICOL)
              ENDDO
           ENDDO
#ifdef scaLAPACK
           PROFILING_START('pdsyev_eigen_blea_mb')
           CALL SCALAPACK_WORK_ARRAY_SIZE_PDSYEV_INIT(PAR_SUP_HANDLE%SCAWORK)
           CALL PDSYEV('V','U',NROW_CMAT,CMAT(1,1),1,1,PAR_SUP_HANDLE%DESC_CMAT_FFM%DESC, &
                EIGENVALUE(1),EIGENVECTOR(1,1),1,1,PAR_SUP_HANDLE%DESC_CMAT_FFM%DESC, &
                PAR_SUP_HANDLE%SCAWORK%WORK(1),PAR_SUP_HANDLE%SCAWORK%LWORK,PAR_SUP_HANDLE%SCAWORK%INFO)
           CALL SCALAPACK_WORK_ARRAY_SIZE_PDSYEV_FIN(PAR_SUP_HANDLE%SCAWORK)
! Then, calculate the eigenvalues and eigenvectors of covariance matrix
           CALL SCALAPACK_WORK_ARRAY_PDSYEV_INIT(PAR_SUP_HANDLE%SCAWORK)
           CALL PDSYEV('V','U',NROW_CMAT,CMAT(1,1),1,1,PAR_SUP_HANDLE%DESC_CMAT_FFM%DESC, &
                EIGENVALUE(1),EIGENVECTOR(1,1),1,1,PAR_SUP_HANDLE%DESC_CMAT_FFM%DESC, &
                PAR_SUP_HANDLE%SCAWORK%WORK(1),PAR_SUP_HANDLE%SCAWORK%LWORK,PAR_SUP_HANDLE%SCAWORK%INFO)
           CALL SCALAPACK_WORK_ARRAY_PDSYEV_FIN(PAR_SUP_HANDLE%SCAWORK)
           PROFILING_STOP('pdsyev_eigen_blea_mb')
#else
           PROFILING_START('dsyev_eigen_blea_mb')
! First, calculate the size of work array
           CALL LAPACK_WORK_ARRAY_SIZE_DSYEV_INIT(LAPACKWORK)
           CALL DSYEV ('V','U',NROW_CMAT,CMAT(1,1),MROW_CMAT,EIGENVALUE(1),LAPACKWORK%WORK(1),LAPACKWORK%LWORK,LAPACKWORK%INFO)
           CALL LAPACK_WORK_ARRAY_SIZE_DSYEV_FIN(LAPACKWORK)
! Then, calculate the eigenvalues and eigenvectors of covariance matrix
           CALL LAPACK_WORK_ARRAY_DSYEV_INIT(LAPACKWORK)
           CALL DSYEV ('V','U',NROW_CMAT,CMAT(1,1),MROW_CMAT,EIGENVALUE(1),LAPACKWORK%WORK(1),LAPACKWORK%LWORK,LAPACKWORK%INFO)
           CALL LAPACK_WORK_ARRAY_DSYEV_FIN(LAPACKWORK)
           PROFILING_STOP ('dsyev_eigen_blea_mb')
#endif
! Set WMAX and WMIN
           WMAX=EIGENVALUE(1)
           WMIN=EIGENVALUE(1)
           DO IROW=2, NROW_CMAT, 1
              WMAX=MAX(WMAX,ABS(EIGENVALUE(IROW)))
!gK              WMIN=MIN(WMIN,ABS(EIGENVALUE(IROW)))
              WMIN=MIN(WMIN,EIGENVALUE(IROW))
           ENDDO
!----------------------------------------------------------------------------------------------------
! Optimization of SIGV and SIGW
!----------------------------------------------------------------------------------------------------
           ! Set convergence parameter EPS_REG to user defined input
           EPS_REG=MAX(ML_EPS_REG,EPS_REG_LAST/2)
           QUOTIENT_OLD=0.0_q
           ITE_RESET=10

           ! Initialize pivot array for linear equations here (PDGETRF, PDGETRI, PDGESV)
#ifdef scaLAPACK
           CALL SCALAPACK_WORK_ARRAY_PDGETRF_INIT (PAR_SUP_HANDLE%SCAWORK,NBLOCK_SCALAPACK,MROW_CMAT, &
                               PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW)
#else
           CALL LAPACK_WORK_ARRAY_DGETRF_INIT (LAPACKWORK,NROW_CMAT)
#endif
           DO ITE=1, NITE
              CALL EVIDENCE_APPROXIMATION()

              IF (ITE==1) THEN 
                 IF (EPS_REG_LAST.EQ.0.0_q) THEN
                    ! If this routine is called the first time one cannot expect a reliable error estimate ZTR
                    ! so we set ZTR_MAX to a large value
                    ZTR_MAX=1E+08_q
                 ELSE
                    ! We expect a reliable error estimate ZTR using SIGW/SIGV from previous calls
                    ZTR_MAX=ZTR*1.2_q
                 ENDIF
              ELSE
                 ZTR_MAX=MAX(ZTR_OLD*1.2_q,ZTR_MAX)
              ENDIF
! Check convergence of SIGW/SIGV and convergence of error estimate ZTR
              IF((ITE .NE. 1) .AND.(ABS(SIGW/SIGV-QUOTIENT_OLD).LT.CRITERIA*MAX(ABS(SIGW/SIGV),ABS(QUOTIENT_OLD))) & 
                 .AND. (ABS(ZTR-ZTR_OLD) .LT. CRITERIA*MAX(ABS(ZTR),ABS(ZTR_OLD)))) THEN
!              IF((ITE .NE. 1) .AND. (ABS(EF-EF_OLD) .LT. CRITERIA*MAX(ABS(EF),ABS(EF_OLD)))) THEN
!                   .AND. (ABS(SIGW-SIGW_OLD) .LT. 0.2*MAX(ABS(SIGW),ABS(SIGW_OLD))) ) THEN
                 LEXIT = .TRUE.
              ELSE
                 LEXIT = .FALSE.
              ENDIF
              ITE_RESET=ITE_RESET-1
              ! If the error estimate ZTR exceeds a threshold increase the regularization
              ! and avoid exit of the optimization loop
              IF ((ITE.NE.1).AND.(ABS(ZTR)>ABS(ZTR_MAX))) THEN
                 EPS_REG=EPS_REG*2.0_q
                 ITE_RESET=10
                 LEXIT = .FALSE.
              ELSE IF (ITE_RESET.LE.0) THEN
                 EPS_REG=EPS_REG*4.0_q
                 ITE_RESET=10
              ENDIF
              ! This checks if regularization becomes gigantic or if no
              ! convergence is reached in the iteration
              IF((ITE .EQ. NITE) .OR. (ABS(SIGV/SIGW) .LT. EPS4)) THEN
                 IF(IREG.EQ.2.OR.IREG.EQ.3) THEN
                    IFF=.FALSE.
                 ENDIF
                 LEXIT = .TRUE.
              ENDIF
              CALL M_BCAST (PAR_SUP_HANDLE%COMM_WORLD,LEXIT,0)
! We also need to broadcast IFF here, if the LEXIT from core 0 is enforced
              CALL M_BCAST (PAR_SUP_HANDLE%COMM_WORLD,IFF,0)
              CALL LOGLINE_REGR(LOGF,                   &
                                NSTEP,                  &
                                IALGO_LINREG,           &
                                ITE,                    &
                                1 / SIGW,               &
                                1 / SIGV,               &
                                ( SIGW / SIGV ) / WMAX, &
                                WMAX / WMIN,            &
                                WMAX,                   &
                                ZTR,                    &
                                LEXIT)
              QUOTIENT_OLD=SIGW/SIGV
! Only if LEXIT is .TRUE. we exit, it was broadcasted so no hangup can occur here
              IF (LEXIT) EXIT
           ENDDO
           ! Before exiting we need to calculate the inverse of the covariance matrix
           !  which is later needed to calculate the Bayesian error estimates
           CALL CALCULATE_INVERSE_OF_COVARIANCE_MATRIX()
#ifdef scaLAPACK
           CALL SCALAPACK_WORK_ARRAY_PDGETRF_FIN(PAR_SUP_HANDLE%SCAWORK)
#else
           CALL LAPACK_WORK_ARRAY_DGETRF_FIN(LAPACKWORK)
#endif
           EPS_REG_LAST=EPS_REG
!----------------------------------------------------------------------------------------------------
        ENDIF
! Rescaling
        PROFILING_START('rescale_blea_mb')
        DO INTYP0=1, NTYP
           DO ICOL=1, NCOL_FMAT
#ifdef scaLAPACK
              TEST_COL=MOD((ICOL-1)/PAR_SUP_HANDLE%DESC_FMAT_FFM%NB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL)
              IF(PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL.NE.TEST_COL) THEN
                 CYCLE
              ENDIF
              LOCVAR_COL=GLOBAL2LOCAL(ICOL,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_FMAT_FFM%NB)
#else
              LOCVAR_COL=ICOL
#endif
              DO IROW=1, NROW_FMAT(INTYP0)
#ifdef scaLAPACK
                 TEST_ROW=MOD((IROW-1)/PAR_SUP_HANDLE%DESC_FMAT_FFM%MB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW)
                 IF(PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW.NE.TEST_ROW) THEN
                    CYCLE
                 ENDIF
                 LOCVAR_ROW=GLOBAL2LOCAL(IROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%DESC_FMAT_FFM%MB)
#else
                 LOCVAR_ROW=IROW
#endif
                 FMAT(LOCVAR_ROW,LOCVAR_COL,INTYP0)=FMAT(LOCVAR_ROW,LOCVAR_COL,INTYP0)/DMAT_FMAT(LOCVAR_COL)
              ENDDO
           ENDDO
        ENDDO
        DO IROW=1, NROW_YMAT
#ifdef scaLAPACK
           TEST_ROW=MOD((IROW-1)/PAR_SUP_HANDLE%DESC_YMAT_FFM%MB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW)
           IF(PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW.NE.TEST_ROW) THEN
              CYCLE
           ENDIF
           LOCVAR_ROW=GLOBAL2LOCAL(IROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%DESC_FMAT_FFM%MB)
           LOCVAR_COL=GLOBAL2LOCAL(   1,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_FMAT_FFM%NB)
#else
           LOCVAR_ROW=IROW
           LOCVAR_COL=1
#endif
           YMAT(LOCVAR_ROW,LOCVAR_COL)=YMAT(LOCVAR_ROW,LOCVAR_COL)/DMAT_YMAT(LOCVAR_ROW)
        ENDDO
        PROFILING_STOP('rescale_blea_mb')
        PROFILING_STOP ('blea_mb')
1       FORMAT(100('/')/5X,'Start regression of many-body part'/100('=')/ &
                       5X,'ITE  Nbasis    1/SIGW^2       1/SIGV^2          EF',15X,'WMAX',10X,'WMAX/WMIN'/100('-'))
2       FORMAT(I7,I7,E16.6e2,E16.6e2,E16.6e2,E16.6e2,E16.6e2)
3       FORMAT(5X,'!!! Initialize noise- and precision parameters !!!')
4       FORMAT(5X,'Condition number of design matrix from SVD: ',E16.6e2)

        CONTAINS 

!****************************************************************************************************
!
!  subroutine to perform one step in the evidence approximation
!  input:  SIGV, SIGW, EPS_REG ...
!  output: update SIGV, SIGW
!
!****************************************************************************************************

           SUBROUTINE EVIDENCE_APPROXIMATION
              IMPLICIT NONE

              PROFILING_START('renew_cmat_blea_mb')

! Calculate covariance matrix
              CMAT=0.0_q
              DO ICOL=1, NCOL_CMAT
#ifdef scaLAPACK
                 CALL PDAXPY(NROW_CMAT,SIGV,CMAT0(1,1),1,ICOL,PAR_SUP_HANDLE%DESC_CMAT_FFM%DESC,1,CMAT(1,1),1,ICOL,PAR_SUP_HANDLE%DESC_CMAT_FFM%DESC,1)
#else
                 CALL DAXPY (NROW_CMAT,SIGV,CMAT0(1,ICOL),1,CMAT(1,ICOL),1)
#endif
              ENDDO
              DO ICOL=1, NCOL_CMAT
#ifdef scaLAPACK
                 TEST_ROW =MOD((ICOL-1)/PAR_SUP_HANDLE%DESC_CMAT_FFM%MB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW)
                 TEST_COL =MOD((ICOL-1)/PAR_SUP_HANDLE%DESC_CMAT_FFM%NB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL)
                 IF (PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL.NE.TEST_COL.OR.PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW.NE.TEST_ROW) THEN
                    CYCLE
                 ENDIF
                 LOCVAR_ROW=GLOBAL2LOCAL(ICOL,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%DESC_CMAT_FFM%MB)
                 LOCVAR_COL=GLOBAL2LOCAL(ICOL,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_CMAT_FFM%NB)
#else
                 LOCVAR_ROW=ICOL
                 LOCVAR_COL=ICOL
#endif
                 CMAT(LOCVAR_ROW,LOCVAR_COL)=CMAT(LOCVAR_ROW,LOCVAR_COL)+SIGW
              ENDDO
              PROFILING_STOP('renew_cmat_blea_mb')
! Solve linear equations
              PROFILING_START('linear_equations_blea_mb')
#ifdef scaLAPACK
              ! First do the LU factorization
              CALL PDGETRF(NROW_CMAT,NROW_CMAT,CMAT(1,1),1,1,PAR_SUP_HANDLE%DESC_CMAT_FFM%DESC,PAR_SUP_HANDLE%SCAWORK%IPIV,PAR_SUP_HANDLE%SCAWORK%INFO)
              ! We also need to copy GMAT to WMATALL, since the input vector is also the output vector
              ! WMATALL=SIGV*Phi^T*Y=SIGV*GMAT by SIGV
              CALL PDCOPY(NROW_CMAT,GMAT(1,1),1,1,PAR_SUP_HANDLE%DESC_WMAT%DESC,1,WMATALL(1,1),1,1,PAR_SUP_HANDLE%DESC_WMAT%DESC,1)
              ! Rescale matrix by SIGV
              CALL PDSCAL(NROW_CMAT,SIGV,WMATALL(1,1),1,1,PAR_SUP_HANDLE%DESC_WMAT%DESC,1)
              ! Now the actual solving of the equations is done here
              CALL PDGETRS('N',NROW_CMAT,1,CMAT(1,1),1,1,PAR_SUP_HANDLE%DESC_CMAT_FFM%DESC,PAR_SUP_HANDLE%SCAWORK%IPIV(1), &
                          WMATALL(1,1),1,1,PAR_SUP_HANDLE%DESC_WMAT%DESC,PAR_SUP_HANDLE%SCAWORK%INFO)
#else
              CALL DGETRF(NROW_CMAT,NROW_CMAT,CMAT(1,1),MROW_CMAT,LAPACKWORK%IPIV(1),LAPACKWORK%INFO)
              CALL DCOPY(NROW_CMAT,GMAT(1,1),1,WMATALL(1,1),1)
              CALL DSCAL(NROW_CMAT,SIGV,WMATALL(1,1),1)
              CALL DGETRS('N',NROW_CMAT,1,CMAT(1,1),MROW_CMAT,LAPACKWORK%IPIV(1), &
                          WMATALL(1,1),MROW_CMAT,LAPACKWORK%INFO)
#endif
              DO INTYP0=1, NTYP
#ifdef scaLAPACK
                 CALL PDCOPY(NROW_WMAT(INTYP0),WMATALL(1,1),1+LBHEAD(INTYP0),1,PAR_SUP_HANDLE%DESC_WMAT%DESC,1,WMAT(1,1,INTYP0),1,1,PAR_SUP_HANDLE%DESC_WMAT_FFM%DESC,1)
#else
                 CALL DCOPY (NROW_WMAT(INTYP0),WMATALL(1+LBHEAD(INTYP0),1),1,WMAT(1,1,INTYP0),1)
#endif
              ENDDO
! Estimate mean values
              RM=0.0_q
              DO INTYP0=1, NTYP
#ifdef scaLAPACK
                 CALL PDGEMV('T',NROW_FMAT(INTYP0),NCOL_FMAT,1.0_q,FMAT(1,1,INTYP0),1,1,PAR_SUP_HANDLE%DESC_FMAT_FFM%DESC,WMAT(1,1,INTYP0),1,1,PAR_SUP_HANDLE%DESC_WMAT_FFM%DESC,1,1.0_q,RM(1,1),1,1,PAR_SUP_HANDLE%DESC_YMAT_FFM%DESC,1)
#else
                 CALL DGEMV ('T',NROW_FMAT(INTYP0),NCOL_FMAT,1.0_q,FMAT(1,1,INTYP0),MROW_FMAT,WMAT(1,1,INTYP0),1,1.0_q,RM(1,1),1)
#endif
              ENDDO
              PROFILING_STOP('linear_equations_blea_mb')
! Calculate loss function
              PROFILING_START('renew_ztr_blea_mb')
              IF(ITE .NE. 1) THEN
                 EF_OLD=EF
                 SIGW_OLD=SIGW
                 ZTR_OLD=ZTR
              ELSE
                 SIGW_OLD=0.0_q
                 ZTR_OLD=0.0_q
              ENDIF

              DO IROW=1, NROW_YMAT
#ifdef scaLAPACK
                 TEST_ROW =MOD((IROW-1)/PAR_SUP_HANDLE%DESC_YMAT_FFM%MB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW)
                 TEST_COL =MOD((   1-1)/PAR_SUP_HANDLE%DESC_YMAT_FFM%NB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL)
                 IF (PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL.NE.TEST_COL.OR.PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW.NE.TEST_ROW) THEN
                    CYCLE
                 ENDIF
                 LOCVAR_ROW=GLOBAL2LOCAL(IROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%DESC_YMAT_FFM%MB)
                 LOCVAR_COL=GLOBAL2LOCAL(   1,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_YMAT_FFM%NB)
#else
                 LOCVAR_ROW=IROW
                 LOCVAR_COL=1
#endif
                 VEC1(LOCVAR_ROW,LOCVAR_COL)=RM(LOCVAR_ROW,LOCVAR_COL)-YMAT(LOCVAR_ROW,LOCVAR_COL)
              ENDDO
#ifdef scaLAPACK
              CALL PDNRM2(NROW_YMAT,ZTR,VEC1(1,1),1,1,PAR_SUP_HANDLE%DESC_YMAT_FFM%DESC,1)
              IF(PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL.EQ.0) THEN
                 NPROCS=1
              ELSE
                 NPROCS=0
              ENDIF
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,NPROCS,1,MPI_INTEGER,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,ZTR   ,1,MPI_REAL8  ,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
              ZTR=ZTR/REAL(NPROCS,q)
              ZTR=5.0E-1_q*ZTR**2
#else
              ZTR=5.0E-1_q*(DNRM2(NROW_YMAT,VEC1(1,1),1))**2
#endif
              PROFILING_STOP('renew_ztr_blea_mb')
!----------------------------------------------------------------------------------------------------
! Calculate evidence function and store marginal likelhood in EF
!----------------------------------------------------------------------------------------------------
              PROFILING_START('renew_ef_blea_mb')
              DO IROW=1, NROW_CMAT
#ifdef scaLAPACK
                 TEST_ROW =MOD((IROW-1)/PAR_SUP_HANDLE%DESC_WMAT%MB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW)
                 TEST_COL =MOD((   1-1)/PAR_SUP_HANDLE%DESC_WMAT%NB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL)
                 IF (PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL.NE.TEST_COL.OR.PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW.NE.TEST_ROW) THEN
                    CYCLE
                 ENDIF
                 LOCVAR_ROW=GLOBAL2LOCAL(IROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%DESC_WMAT%MB)
                 LOCVAR_COL=GLOBAL2LOCAL(   1,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_WMAT%NB)
#else
                 LOCVAR_ROW=IROW
                 LOCVAR_COL=1
#endif
!                 VEC3(LOCVAR_ROW,LOCVAR_COL)=-5.0E-1_q*LOG(SIGV*EIGENVALUE(IROW)+SIGW)
                 IF(EIGENVALUE(IROW)/WMAX .GT. EPS_REG) THEN
                    VEC3(LOCVAR_ROW,LOCVAR_COL)=-5.0E-1_q*LOG(SIGV*EIGENVALUE(IROW)+SIGW)
                 ELSE
!gK                    VEC3(LOCVAR_ROW,LOCVAR_COL)=0.0_q
                    VEC3(LOCVAR_ROW,LOCVAR_COL)=-5.0E-1_q*LOG(SIGW)
                 ENDIF
              ENDDO
#ifdef scaLAPACK
              CALL PDDOT(NROW_CMAT,EF,VECONE(1,1),1,1,PAR_SUP_HANDLE%DESC_WMAT%DESC,1,VEC3(1,1),1,1,PAR_SUP_HANDLE%DESC_WMAT%DESC,1)
              IF(PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL.EQ.0) THEN
                 NPROCS=1
              ELSE
                 NPROCS=0
              ENDIF
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,NPROCS,1,MPI_INTEGER,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,EF    ,1,MPI_REAL8  ,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
              EF=EF/REAL(NPROCS,q)
#else
              EF=DDOT(NROW_CMAT,VECONE(1,1),1,VEC3(1,1),1)
#endif
              EF=EF+0.5_q*REAL(NROW_CMAT,q)*LOG(SIGW)+0.5_q*REAL(NCOL_FMAT,q)*LOG(SIGV)-0.5_q*REAL(NCOL_FMAT,q)*LOG(2.0_q*PI)-SIGV*ZTR
#ifdef scaLAPACK
              CALL PDNRM2(NROW_CMAT,EF_HELP,WMATALL(1,1),1,1,PAR_SUP_HANDLE%DESC_WMAT%DESC,1)
              IF(PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL.EQ.0) THEN
                 NPROCS=1
              ELSE
                 NPROCS=0
              ENDIF
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,NPROCS ,1,MPI_INTEGER,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,EF_HELP,1,MPI_REAL8  ,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
              EF_HELP=EF_HELP/REAL(NPROCS,q)
              EF=EF-0.5_q*SIGW*EF_HELP**2
#else
              EF=EF-0.5_q*SIGW*(DNRM2(NROW_CMAT,WMATALL(1,1),1))**2
#endif
!----------------------------------------------------------------------------------------------------
! gamma in equation (C3) of PRB 100, 014105 (2019) sum_k lambda_k/(lambda_k + 1 / sigma_w^2)
! where lambda_k are the eigenvalues of Psî^T Psi / sigma_v^2
!----------------------------------------------------------------------------------------------------
              PROFILING_STOP('renew_ef_blea_mb')
! Calculate of gradients of evidence function gradient with respect to hyper-parameters
              PROFILING_START('renew_sig_blea_mb')
              DO IROW=1, NROW_CMAT
#ifdef scaLAPACK
                 TEST_ROW =MOD((IROW-1)/PAR_SUP_HANDLE%DESC_WMAT%MB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW)
                 TEST_COL =MOD((   1-1)/PAR_SUP_HANDLE%DESC_WMAT%NB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL)
                 if (PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL.NE.TEST_COL.OR.PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW.NE.TEST_ROW) THEN
                    CYCLE
                 ENDIF
                 LOCVAR_ROW=GLOBAL2LOCAL(IROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%DESC_WMAT%MB)
                 LOCVAR_COL=GLOBAL2LOCAL(   1,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_WMAT%NB)
#else
                 LOCVAR_ROW=IROW
                 LOCVAR_COL=1
#endif
!                 VEC3(LOCVAR_ROW,LOCVAR_COL)=EIGENVALUE(IROW)/(EIGENVALUE(IROW)+SIGW/SIGV)
                 IF(EIGENVALUE(IROW)/WMAX .GT. EPS_REG) THEN
                    VEC3(LOCVAR_ROW,LOCVAR_COL)=SIGV*EIGENVALUE(IROW)/(SIGV*EIGENVALUE(IROW)+SIGW)
                 ELSE
!gK                    VEC3(LOCVAR_ROW,LOCVAR_COL)=0.0_q
                    VEC3(LOCVAR_ROW,LOCVAR_COL)=SIGV*EIGENVALUE(IROW)/(SIGW)
                 ENDIF
              ENDDO
#ifdef scaLAPACK
              CALL PDDOT(NROW_CMAT,GAM,VECONE(1,1),1,1,PAR_SUP_HANDLE%DESC_WMAT%DESC,1,VEC3(1,1),1,1,PAR_SUP_HANDLE%DESC_WMAT%DESC,1)
              IF(PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL.EQ.0) THEN
                 NPROCS=1
              ELSE
                 NPROCS=0
              ENDIF
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,NPROCS ,1,MPI_INTEGER,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,GAM    ,1,MPI_REAL8  ,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
              GAM=GAM/REAL(NPROCS,q)
#else
              GAM=DDOT(NROW_CMAT,VECONE(1,1),1,VEC3(1,1),1)
#endif
!----------------------------------------------------------------------------------------------------
!  equation (C1) of PRB 100, 014105 (2019):  ANORM = sum_k |w_k|^2
!----------------------------------------------------------------------------------------------------
              DO IROW=1, NROW_CMAT
#ifdef scaLAPACK
                 TEST_ROW =MOD((IROW-1)/PAR_SUP_HANDLE%DESC_WMAT%MB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW)
                 TEST_COL =MOD((   1-1)/PAR_SUP_HANDLE%DESC_WMAT%NB,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL)
                 if (PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL.NE.TEST_COL.OR.PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW.NE.TEST_ROW) THEN
                    CYCLE
                 ENDIF
                 LOCVAR_ROW=GLOBAL2LOCAL(IROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW,PAR_SUP_HANDLE%DESC_WMAT%MB)
                 LOCVAR_COL=GLOBAL2LOCAL(   1,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_COL,PAR_SUP_HANDLE%DESC_WMAT%NB)
#else
                 LOCVAR_ROW=IROW
                 LOCVAR_COL=1
#endif
                 VEC3(LOCVAR_ROW,LOCVAR_COL)=WMATALL(LOCVAR_ROW,LOCVAR_COL)
              ENDDO
#ifdef scaLAPACK
              CALL PDNRM2(NROW_CMAT,ANORM,VEC3(1,1),1,1,PAR_SUP_HANDLE%DESC_WMAT%DESC,1)
              IF(PAR_SUP_HANDLE%CONTEXT_WORLD%ME_COL.EQ.0) THEN
                 NPROCS=1
              ELSE
                 NPROCS=0
              ENDIF
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,NPROCS ,1,MPI_INTEGER,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
              CALL MPI_ALLREDUCE(MPI_IN_PLACE,ANORM  ,1,MPI_REAL8  ,MPI_SUM,PAR_SUP_HANDLE%COMM_WORLD%MPI_COMM,IERR)
              ANORM=ANORM/REAL(NPROCS,q)
              ANORM=5.0E-1_q*ANORM**2
#else
              ANORM=5.0E-1_q*DNRM2(NROW_CMAT,VEC3(1,1),1)**2
#endif
             IF(IREG .EQ. 2) THEN
                 SIGV=0.5_q*(REAL(NCOL_FMAT,q)-GAM)/ZTR
                 SIGW=MAX(ABS(EPS_REG*WMAX)*SIGV,0.5_q*GAM/ANORM)
              ELSE IF(IREG.EQ.3) THEN
                 SIGW=MAX(ABS(EPS_REG*WMAX)*SIGV,0.5_q*GAM/ANORM)
              ENDIF
! SIGV and SIGW need to be broadcasted here to ensure they are the same on all cores
              CALL M_BCAST(PAR_SUP_HANDLE%COMM_WORLD,SIGV,0)
              CALL M_BCAST(PAR_SUP_HANDLE%COMM_WORLD,SIGW,0)

              PROFILING_STOP('renew_sig_blea_mb')

           END SUBROUTINE EVIDENCE_APPROXIMATION
              
!****************************************************************************************************
!  subroutine to calculate inverse of covariance matrix CMAT
!  this subroutine requires the pivoting information from LU decomposition 
!  (PAR_SUP_HANDLE%SCAWORK%IPIV) to be calculated before
!  input:  CMAT, PAR_SUP_HANDLE%SCAWORK%IPIV
!  output: CMAT
!****************************************************************************************************
           SUBROUTINE CALCULATE_INVERSE_OF_COVARIANCE_MATRIX
              IMPLICIT NONE
              PROFILING_START('inverse_cov_mat_blea_mb')
#ifdef scaLAPACK
              ! Here we calculate the optimal sizes of the work arrays
              CALL SCALAPACK_WORK_ARRAY_SIZE_PDGETRI_INIT (PAR_SUP_HANDLE%SCAWORK,NBLOCK_SCALAPACK,MROW_CMAT, &
                   PAR_SUP_HANDLE%CONTEXT_WORLD%ME_ROW,PAR_SUP_HANDLE%CONTEXT_WORLD%NP_ROW)
              CALL PDGETRI(NROW_CMAT,CMAT(1,1),1,1,PAR_SUP_HANDLE%DESC_CMAT_FFM%DESC,PAR_SUP_HANDLE%SCAWORK%IPIV(1),PAR_SUP_HANDLE%SCAWORK%WORK(1),PAR_SUP_HANDLE%SCAWORK%LWORK,PAR_SUP_HANDLE%SCAWORK%IWORK(1),PAR_SUP_HANDLE%SCAWORK%LIWORK,PAR_SUP_HANDLE%SCAWORK%INFO)
              CALL SCALAPACK_WORK_ARRAY_SIZE_PDGETRI_FIN(PAR_SUP_HANDLE%SCAWORK)
              ! Here is the actual inversion done
              CALL SCALAPACK_WORK_ARRAY_PDGETRI_INIT(PAR_SUP_HANDLE%SCAWORK)
              CALL PDGETRI(NROW_CMAT,CMAT(1,1),1,1,PAR_SUP_HANDLE%DESC_CMAT_FFM%DESC,PAR_SUP_HANDLE%SCAWORK%IPIV(1),PAR_SUP_HANDLE%SCAWORK%WORK(1),PAR_SUP_HANDLE%SCAWORK%LWORK,PAR_SUP_HANDLE%SCAWORK%IWORK(1),PAR_SUP_HANDLE%SCAWORK%LIWORK,PAR_SUP_HANDLE%SCAWORK%INFO)
              CALL SCALAPACK_WORK_ARRAY_PDGETRI_FIN(PAR_SUP_HANDLE%SCAWORK)
#else
              CALL LAPACK_WORK_ARRAY_SIZE_DGETRI_INIT(LAPACKWORK)
              CALL DGETRI(NROW_CMAT,CMAT(1,1),MROW_CMAT,LAPACKWORK%IPIV(1),LAPACKWORK%WORK(1),LAPACKWORK%LWORK,LAPACKWORK%INFO)
              CALL LAPACK_WORK_ARRAY_SIZE_DGETRI_FIN(LAPACKWORK)
              CALL LAPACK_WORK_ARRAY_DGETRI_INIT(LAPACKWORK)
              CALL DGETRI(NROW_CMAT,CMAT(1,1),MROW_CMAT,LAPACKWORK%IPIV(1),LAPACKWORK%WORK(1),LAPACKWORK%LWORK,LAPACKWORK%INFO)
              CALL LAPACK_WORK_ARRAY_DGETRI_FIN(LAPACKWORK)
#endif
              PROFILING_STOP('inverse_cov_mat_blea_mb')
           END SUBROUTINE CALCULATE_INVERSE_OF_COVARIANCE_MATRIX

         ! CONTAINS ends here
 
      END SUBROUTINE BLEA_MB

      END MODULE REGRESSION
#endif

#endif
