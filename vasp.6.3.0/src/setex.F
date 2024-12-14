#include "symbol.inc"
!***********************************************************************
!
! VASP calles the routines in the xclib module  only via the
! subroutines in the module SETEXM
! (SETEXM is the actual interface layer between exchange correlation
!  functionals and VASP)
! it stores for instances which exchange correlation type
! VASP uses, or which interpolation is used (LFCI)
!
!***********************************************************************

  MODULE setexm
    USE prec
    USE setexm_struct_def
    !
    ! LEXCH_TABLE specifies the exchange correlation type stored
    ! in the table EXCTAB
    ! the size of the table is specified by NSMA
    !
    INTEGER, PARAMETER :: NSMA=2000
    INTEGER, SAVE, PRIVATE :: LEXCH_TABLE=-1
    !
    ! the exchange stack can be used to save 
    ! the present exchange parameters temporarily
    !
    INTEGER, SAVE, PRIVATE :: ISTACK=0
    REAL(q), SAVE, PRIVATE :: EX_STACK(7,5)
    LOGICAL, SAVE, PRIVATE :: EX_LSTACK(4,5)

    CONTAINS

!******************* SUBROUTINE EXTYP *********************************
!
!  this subroutine interprets the string CEXCH
!  which determines the type of exchange correlation
!  and sets the integer LEXCH accordingly
!
!**********************************************************************

    SUBROUTINE EXTYP(CEXCH,LEXCH)
      USE prec
      USE main_mpi
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      CHARACTER (2) CEXCH
      INTEGER LEXCH
#ifdef libbeef
      LOGICAL, EXTERNAL :: BEEF_SET_TYPE
#endif

      LEXCH=-1
      IF (CEXCH=='  ') THEN
        LEXCH=0
      ELSE IF (CEXCH=='HL') THEN
        LEXCH=1
      ELSE IF (CEXCH=='PZ') THEN
        LEXCH=2
      ELSE IF (CEXCH=='CA') THEN
        LEXCH=2
      ELSE IF (CEXCH=='WI') THEN
        LEXCH=3
      ELSE IF (CEXCH=='PB') THEN
        LEXCH=4
      ELSE IF (CEXCH=='PW') THEN
        LEXCH=5
      ELSE IF (CEXCH=='LM') THEN
        LEXCH=6
      ELSE IF (CEXCH=='91') THEN
        LEXCH=7
      ELSE IF (CEXCH=='PE') THEN
        LEXCH=8
      ELSE IF (CEXCH=='RP') THEN
        LEXCH=9
      ELSE IF (CEXCH=='VW') THEN
        LEXCH=10
      ELSE IF (CEXCH=='B3') THEN
        LEXCH=11
      ELSE IF (CEXCH=='B5') THEN
        LEXCH=12
!aem the AM05 functional added
      ELSE IF (CEXCH=='AM') THEN
        LEXCH=13
!aem the AM05 functional added
! jP: adding PBEsol
      ELSE IF (CEXCH=='PS') THEN
        LEXCH=14
! jP: adding PBEsol
      ELSE IF (CEXCH=='CO') THEN
        LEXCH=100
! BEEF XC
      ELSE IF (CEXCH=='BF') THEN
#ifdef libbeef
        LEXCH=17
#if defined(MPI) || defined(MPI_CHAIN)
        IF (.NOT.BEEF_SET_TYPE(0, (COMM%NODE_ME==COMM%IONODE))) THEN
#else
        IF (.NOT.BEEF_SET_TYPE(0, .TRUE.)) THEN
#endif
           CALL vtutor%error('error setting beef type 0')
        ENDIF
#else
        CALL vtutor%error("VASP needs to be linked against libbeef for Bayesian error estimation &
           &functional support.\nlibbeef sources and binaries can be downloaded from &
           &suncat.stanford.edu")
#endif
! functionals for range-separated ACFDT (LDA - short range RPA):
! a bit akward at the moment since the range separation parameter
! is hard coded for now. Hopefully this will change ...
      ELSE IF (CEXCH=='RA') THEN
! jH-new RPA Perdew-Wang
        LEXCH=20
      ELSE IF (CEXCH=='03') THEN
        ! \mu = 0.3 A^-1
        LEXCH=21
      ELSE IF (CEXCH=='05') THEN
        ! \mu = 0.5 A^-1
        LEXCH=22
      ELSE IF (CEXCH=='10') THEN
        ! \mu = 1.0 A^-1
        LEXCH=23
      ELSE IF (CEXCH=='20') THEN
        ! \mu = 2.0 A^-1
        LEXCH=24
      ELSE IF (CEXCH=='PL') THEN
! jH-new RPAplus Perdew-Wang         
        LEXCH=30
!vdw jk
      ELSE IF (CEXCH=='RE') THEN
        LEXCH=40
      ELSE IF (CEXCH=='OR') THEN
        LEXCH=41
      ELSE IF (CEXCH=='BO') THEN
        LEXCH=42
      ELSE IF (CEXCH=='MK') THEN
        LEXCH=43
      ELSE IF (CEXCH=='ML') THEN
        LEXCH=44
!vdw jk
!vdw TbB
      ELSE IF (CEXCH=='CX') THEN
        LEXCH=45
!vdw TbB
#ifdef USELIBXC
      ELSE IF (CEXCH=='LI') THEN
        LEXCH=99
#endif
      ELSE IF (CEXCH/='--') THEN
        CALL vtutor%error("Error: This functional is not implemented.")
      ENDIF
      RETURN
    END SUBROUTINE EXTYP


#ifdef USELIBXC
!******************* SUBROUTINE EXTYP_LIBXC ***************************
!
!  this subroutine initializes the functionals LIBXC1 and LIBXC2
!  of the libxc library (if LIBXC1/LIBXC2 > 0)
!
!**********************************************************************
    SUBROUTINE EXTYP_LIBXC(LIBXC1,LIBXC2)
      USE base, ONLY : in_struct
      USE reader_tags
      USE fock_glb, ONLY : AEXX,LHFCALC
      USE tutor, ONLY : vtutor

      IMPLICIT NONE
      TYPE (in_struct) IO
      INTEGER :: IERR,N,LIBXC1,LIBXC2,LIBXC_TMP
      REAL(q) :: XC_DENS_TMP(1),XC_SIGMA_TMP(1),XC_LAPL_TMP(1),XC_TAU_TMP(1),XC_TAUW_TMP(1),XC_EXC_TMP(1)
      LOGICAL LOPEN
      CHARACTER (40) SZNAM1,SZNAM2
      CHARACTER (50) FUNCNAME
      TYPE (xc_f03_func_t) :: LIBXC_FUNC_TMP

      LIBXC_MASK=0
      XC_FAMILY=LIBXC_FAMILY_NONE
      XC_KIND=-1
      LIBXC_HYB=.FALSE.
      LIBXC_NOEXC=.FALSE.
      LTBOUNDLIBXC=.FALSE.

      IF (LIBXC1<=0.AND.LIBXC2<=0) RETURN

      IF (LIBXC1>0) THEN
         LIBXC_MASK(1)=LIBXC1
         CALL XC_F03_FUNC_INIT(LIBXC_FUNC_NSP(1),LIBXC1,1)
         CALL XC_F03_FUNC_INIT(LIBXC_FUNC_SPN(1),LIBXC1,2)
         LIBXC_INFO_NSP(1)=XC_F03_FUNC_GET_INFO(LIBXC_FUNC_NSP(1))
         LIBXC_INFO_SPN(1)=XC_F03_FUNC_GET_INFO(LIBXC_FUNC_SPN(1))
         XC_FAMILY(1)=XC_F03_FUNC_INFO_GET_FAMILY(LIBXC_INFO_SPN(1))
         FUNCNAME=XC_F03_FUNCTIONAL_GET_NAME(LIBXC1)
         CALL UPPER(FUNCNAME)
         IF (FUNCNAME(1:3) .EQ. 'HYB') LIBXC_HYB(1)=.TRUE.
         XC_KIND(1)=XC_F03_FUNC_INFO_GET_KIND(LIBXC_INFO_SPN(1))
         IF (LIBXC1==208) THEN
            CALL vtutor%error("Error: mBJ/TB09 can not be used via Libxc.")
         ENDIF
         IF (IAND(XC_F03_FUNC_INFO_GET_FLAGS(LIBXC_INFO_NSP(1)),XC_FLAGS_NEEDS_LAPLACIAN)==XC_FLAGS_NEEDS_LAPLACIAN) THEN
            CALL vtutor%error("Error: LIBXC1 is a Laplacian-dependent meta-GGA, but such functionals are not yet implemented in VASP.")
         ENDIF
         IF (IAND(XC_F03_FUNC_INFO_GET_FLAGS(LIBXC_INFO_NSP(1)),XC_FLAGS_HAVE_EXC)==0) THEN
            LIBXC_NOEXC(1)=.TRUE.
         ENDIF
      ENDIF
      IF (LIBXC2>0) THEN
         LIBXC_MASK(2)=LIBXC2
         CALL XC_F03_FUNC_INIT(LIBXC_FUNC_NSP(2),LIBXC2,1)
         CALL XC_F03_FUNC_INIT(LIBXC_FUNC_SPN(2),LIBXC2,2)
         LIBXC_INFO_NSP(2)=XC_F03_FUNC_GET_INFO(LIBXC_FUNC_NSP(2))
         LIBXC_INFO_SPN(2)=XC_F03_FUNC_GET_INFO(LIBXC_FUNC_SPN(2))
         XC_FAMILY(2)=XC_F03_FUNC_INFO_GET_FAMILY(LIBXC_INFO_SPN(2))
         FUNCNAME=XC_F03_FUNCTIONAL_GET_NAME(LIBXC2)
         CALL UPPER(FUNCNAME)
         IF (FUNCNAME(1:3) .EQ. 'HYB') LIBXC_HYB(2)=.TRUE.
         XC_KIND(2)=XC_F03_FUNC_INFO_GET_KIND(LIBXC_INFO_SPN(2))
         IF (LIBXC2==208) THEN
            CALL vtutor%error("Error: mBJ/TB09 can not be used via Libxc.")
         ENDIF
         IF (IAND(XC_F03_FUNC_INFO_GET_FLAGS(LIBXC_INFO_NSP(2)),XC_FLAGS_NEEDS_LAPLACIAN)==XC_FLAGS_NEEDS_LAPLACIAN) THEN
            CALL vtutor%error("Error: LIBXC2 is a Laplacian-dependent meta-GGA, but such functionals are not yet implemented in VASP.")
         ENDIF
         IF (IAND(XC_F03_FUNC_INFO_GET_FLAGS(LIBXC_INFO_NSP(2)),XC_FLAGS_HAVE_EXC)==0) THEN
            LIBXC_NOEXC(2)=.TRUE.
         ENDIF
      ENDIF

      CALL OPEN_INCAR_IF_FOUND(IO%IU5,LOPEN)
      SZNAM1=''
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'GGA', SZNAM1, 40, IERR, WRITEXMLINCAR)
      CALL STRIP(SZNAM1,N,'L'); CALL UPPER(SZNAM1)
      SZNAM2=''
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'METAGGA', SZNAM2, 40, IERR, WRITEXMLINCAR)
      CALL STRIP(SZNAM2,N,'L'); CALL UPPER(SZNAM2)
      CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, 'LTBOUNDLIBXC',LTBOUNDLIBXC , IERR, WRITEXMLINCAR)
      CALL CLOSE_INCAR_IF_FOUND(IO%IU5)

      IF (((XC_FAMILY(1).EQ.LIBXC_FAMILY_LDA) .OR. (XC_FAMILY(1).EQ.LIBXC_FAMILY_HYB_LDA) .OR. (XC_FAMILY(1).EQ.LIBXC_FAMILY_GGA) .OR. (XC_FAMILY(1).EQ.LIBXC_FAMILY_HYB_GGA) .OR. (XC_FAMILY(1).EQ.LIBXC_FAMILY_NONE)) .AND. &
     &    ((XC_FAMILY(2).EQ.LIBXC_FAMILY_LDA) .OR. (XC_FAMILY(2).EQ.LIBXC_FAMILY_HYB_LDA) .OR. (XC_FAMILY(2).EQ.LIBXC_FAMILY_GGA) .OR. (XC_FAMILY(2).EQ.LIBXC_FAMILY_HYB_GGA) .OR. (XC_FAMILY(2).EQ.LIBXC_FAMILY_NONE))) THEN
         IF (SZNAM1(1:5)/='LIBXC') THEN
            CALL vtutor%error("Error: LIBXC1 and LIBXC2 are LDA or GGA, but ""GGA=LIBXC"" in INCAR is missing.")
         ENDIF
         IF (SZNAM2(1:5)=='LIBXC') THEN
            CALL vtutor%error("Error: ""METAGGA=LIBXC"" in INCAR is specified, but neither LIBXC1 nor LIBXC2 are a meta-GGA.")
         ENDIF
      ENDIF
      IF ((XC_FAMILY(1).EQ.LIBXC_FAMILY_MGGA) .OR. (XC_FAMILY(1).EQ.LIBXC_FAMILY_HYB_MGGA) .OR. (XC_FAMILY(2).EQ.LIBXC_FAMILY_MGGA) .OR. (XC_FAMILY(2).EQ.LIBXC_FAMILY_HYB_MGGA)) THEN
         IF (SZNAM2(1:5)/='LIBXC') THEN
            CALL vtutor%error("Error: LIBXC1 and/or LIBXC2 are a meta-GGA, but ""METAGGA=LIBXC"" in INCAR is missing.")
         ENDIF
         IF (SZNAM1(1:5)=='LIBXC') THEN
            CALL vtutor%error("Error: ""GGA=LIBXC"" in INCAR is specified, but LIBXC1 and/or LIBXC2 are a meta-GGA.")
         ENDIF
      ENDIF

      IF ((XC_KIND(1).EQ.XC_EXCHANGE).AND.(XC_KIND(2).EQ.XC_EXCHANGE)) THEN
         CALL vtutor%error("Error: LIBXC1 and LIBXC2 are both exchange (_X_) functionals.")
      ELSEIF ((XC_KIND(1).EQ.XC_CORRELATION).AND.(XC_KIND(2).EQ.XC_CORRELATION)) THEN
         CALL vtutor%error("Error: LIBXC1 and LIBXC2 are both correlation (_C_) functionals.")
      ELSEIF ((XC_KIND(1).EQ.XC_EXCHANGE_CORRELATION).AND.(XC_KIND(2).EQ.XC_EXCHANGE_CORRELATION)) THEN
         CALL vtutor%error("Error: LIBXC1 and LIBXC2 are both exchange-correlation (_XC_) functionals.")
      ELSEIF ((XC_KIND(1).EQ.XC_EXCHANGE_CORRELATION).AND.(XC_KIND(2).NE.-1)) THEN
         CALL vtutor%error("Error: If LIBXC1 is an exchange-correlation (_XC_) functional, then LIBXC2 should not be specified.")
      ELSEIF ((XC_KIND(2).EQ.XC_EXCHANGE_CORRELATION).AND.(XC_KIND(1).NE.-1)) THEN
         CALL vtutor%error("Error: If LIBXC2 is an exchange-correlation (_XC_) functional, then LIBXC1 should not be specified.")
      ELSEIF (XC_KIND(1).EQ.XC_KINETIC) THEN
         CALL vtutor%error("Error: LIBXC1 is a kinetic (_K_) functional.")
      ELSEIF (XC_KIND(2).EQ.XC_KINETIC) THEN
         CALL vtutor%error("Error: LIBXC2 is a kinetic (_K_) functional.")
      ENDIF

      IF (LIBXC_HYB(1).OR.LIBXC_HYB(2)) THEN
         LHFCALC=.TRUE.
         AEXX=0.25_q
      ENDIF

      !A single point calculation with r2SCAN exchange to check that Libxc was compiled with the patch
      IF ((XC_FAMILY(1).EQ.LIBXC_FAMILY_MGGA) .OR. (XC_FAMILY(1).EQ.LIBXC_FAMILY_HYB_MGGA) .OR. (XC_FAMILY(2).EQ.LIBXC_FAMILY_MGGA) .OR. (XC_FAMILY(2).EQ.LIBXC_FAMILY_HYB_MGGA)) THEN
         LIBXC_TMP=497
         CALL XC_F03_FUNC_INIT(LIBXC_FUNC_TMP,LIBXC_TMP,1)
         XC_DENS_TMP=1._q
         XC_SIGMA_TMP=1._q
         XC_LAPL_TMP=0._q
         XC_TAUW_TMP=0.125_q*XC_SIGMA_TMP/XC_DENS_TMP
         XC_TAU_TMP=XC_TAUW_TMP-1000._q
         CALL XC_F03_MGGA_EXC(LIBXC_FUNC_TMP,1_8,XC_DENS_TMP,XC_SIGMA_TMP,XC_LAPL_TMP,XC_TAU_TMP,XC_EXC_TMP)
         IF (ABS(XC_EXC_TMP(1)) < 0.8689_q) THEN
            CALL vtutor%error("Error: Libxc was not compiled with the option --disable-fhc that is important for meta-GGA functionals (see the VASP installation guide for more details).")
         ENDIF
      ENDIF

      RETURN
    END SUBROUTINE EXTYP_LIBXC

!******************* SUBROUTINE SETPARAMS_LIBXC ***********************
!
!  this subroutine sets the parameters of functionals LIBXC1 and LIBXC2
!  of the libxc library (if LIBXC1/LIBXC2 > 0)
!
!**********************************************************************
    SUBROUTINE SETPARAMS_LIBXC
      USE base, ONLY : in_struct
      USE reader_tags
      USE incar_reader, ONLY: COUNT_ELEMENTS
      USE constant, ONLY : AUTOA
      USE fock_glb, ONLY : AEXX,HFSCREEN
      USE setexm_struct_def, ONLY : LIBXC_NPARAM,LIBXC_PARAM,LIBXC_PARAM_DESCR
      USE tutor, ONLY : vtutor

      IMPLICIT NONE
      TYPE (in_struct) IO
      INTEGER :: I,J,N,N_AGGAX,N_ALDAC,N_AGGAC,IERR
      REAL(q) PARAMTMP
      LOGICAL LOPEN
      CHARACTER (LEN=20) :: ICHAR,CHAR1

      LIBXC_NPARAM=0
      LIBXC_PARAM=0._q

      IF (LIBXC_MASK(1)<=0.AND.LIBXC_MASK(2)<=0) RETURN

      CALL OPEN_INCAR_IF_FOUND(IO%IU5,LOPEN)

      N_AGGAX=COUNT_ELEMENTS(INCAR_F,"AGGAX")
      N_ALDAC=COUNT_ELEMENTS(INCAR_F,"ALDAC")
      N_AGGAC=COUNT_ELEMENTS(INCAR_F,"AGGAC")
      IF ((XC_KIND(1).EQ.XC_EXCHANGE_CORRELATION).OR.(XC_KIND(2).EQ.XC_EXCHANGE_CORRELATION)) THEN
         IF ((N_AGGAX >= 1).OR.(N_ALDAC >= 1).OR.(N_AGGAC >= 1)) THEN
            CALL vtutor%error("Error: AGGAX, ALDAC, and AGGAC can not be specified when an exchange-correlation (_XC_) functional is chosen.")
         ENDIF
      ELSEIF ((LIBXC_HYB(1).EQV..TRUE.).OR.(LIBXC_HYB(2).EQV..TRUE.)) THEN
         IF (N_AGGAX >= 1) THEN
            CALL vtutor%error("Error: AGGAX can not be specified when a hybrid (HYB) functional is chosen.")
         ENDIF
      ENDIF

      DO I=1,2
         IF (LIBXC_MASK(I)>0) THEN

            LIBXC_NPARAM(I)=XC_F03_FUNC_INFO_GET_N_EXT_PARAMS(LIBXC_INFO_SPN(I))
            DO J=1,LIBXC_NPARAM(I)
               LIBXC_PARAM_DESCR(I,J)=XC_F03_FUNC_INFO_GET_EXT_PARAMS_DESCRIPTION(LIBXC_INFO_SPN(I),J-1)
               LIBXC_PARAM(I,J)=XC_F03_FUNC_INFO_GET_EXT_PARAMS_DEFAULT_VALUE(LIBXC_INFO_SPN(I),J-1)
            ENDDO

            IF ((LIBXC_MASK(I)==406) .OR. &
           &    (LIBXC_MASK(I)==412) .OR. &
           &    (LIBXC_MASK(I)==427) .OR. &
           &    (LIBXC_MASK(I)==428) .OR. &
           &    (LIBXC_MASK(I)==479) .OR. &
           &    (LIBXC_MASK(I)==480)) THEN
               LIBXC_PARAM(I,1)=AEXX
               CALL XC_F03_FUNC_SET_EXT_PARAMS(LIBXC_FUNC_NSP(I),LIBXC_PARAM(I,1:LIBXC_NPARAM(I)))
               CALL XC_F03_FUNC_SET_EXT_PARAMS(LIBXC_FUNC_SPN(I),LIBXC_PARAM(I,1:LIBXC_NPARAM(I)))
            ENDIF

            IF ((LIBXC_MASK(I)==427) .OR. &
           &    (LIBXC_MASK(I)==428) .OR. &
           &    (LIBXC_MASK(I)==479) .OR. &
           &    (LIBXC_MASK(I)==480)) THEN
               LIBXC_PARAM(I,2)=HFSCREEN*AUTOA
               LIBXC_PARAM(I,3)=HFSCREEN*AUTOA
               CALL XC_F03_FUNC_SET_EXT_PARAMS(LIBXC_FUNC_NSP(I),LIBXC_PARAM(I,1:LIBXC_NPARAM(I)))
               CALL XC_F03_FUNC_SET_EXT_PARAMS(LIBXC_FUNC_SPN(I),LIBXC_PARAM(I,1:LIBXC_NPARAM(I)))
            ENDIF

            DO J=1,LIBXC_NPARAM(I)
               WRITE(ICHAR,*) J
               IF (I==1) THEN
                  CHAR1 = 'LIBXC1_P'//TRIM(ADJUSTL(ICHAR))
               ELSEIF (I==2) THEN
                  CHAR1 = 'LIBXC2_P'//TRIM(ADJUSTL(ICHAR))
               ENDIF
               CALL PROCESS_INCAR(LOPEN, IO%IU0, IO%IU5, TRIM(CHAR1), PARAMTMP, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)
               IF (N>=1) LIBXC_PARAM(I,J)=PARAMTMP
            ENDDO
            IF (LIBXC_NPARAM(I)>=1) THEN
               CALL XC_F03_FUNC_SET_EXT_PARAMS(LIBXC_FUNC_NSP(I),LIBXC_PARAM(I,1:LIBXC_NPARAM(I)))
               CALL XC_F03_FUNC_SET_EXT_PARAMS(LIBXC_FUNC_SPN(I),LIBXC_PARAM(I,1:LIBXC_NPARAM(I)))
            ENDIF

         ENDIF
      ENDDO

      CALL CLOSE_INCAR_IF_FOUND(IO%IU5)

      RETURN
    END SUBROUTINE SETPARAMS_LIBXC
#endif

!***********************************************************************
!
! EX_MOD calculate the exchange energy density per particle eps_x
! possibly reduced by short range exchange hole or reduce by the amount 
! that is accounted for by the exact exchange
! that is 
! E_x = \int eps_x(rho(r)) rho(r) d^3 r
!
! VX_MOD calculates the exchange potential, which is the derivative
! of eps_x(rho(r)) rho(r) with respect to rho(r)
! v_x = d (eps_x(rho) rho) /d rho
! 
! in leading order these quantities are proportional to the density^(1/3)
!
!***********************************************************************


    FUNCTION EX_MOD(RS,IFLG,TREL)
      USE xclib
      USE constant
      IMPLICIT NONE
      REAL(q) EX_MOD, RS
      INTEGER IFLG
      LOGICAL TREL

      EX_MOD=EX(RS,IFLG,TREL)
      
      IF (LDASCREEN/=0 .OR. LDAX/=1 ) THEN
         IF (LDASCREEN==0) THEN
            EX_MOD=EX_MOD*LDAX
         ELSE IF (LUSE_THOMAS_FERMI) THEN
            EX_MOD=EX_MOD-EX_SX(RS,LDASCREEN*AUTOA,IFLG)*(1-LDAX)
         ELSE IF (LUSE_LONGRANGE_HF) THEN
         ! Iann Gerber: use EX_SR (for LDAX=0, i.e. AEXX=1)
!            EX_MOD=(1-LDAX)*EX(RS,IFLG,.FALSE.)-(EX_MOD-EX_SR(RS,LDASCREEN*AUTOA,IFLG))*(1-LDAX)
            EX_MOD=EX(RS,IFLG,.FALSE.)-(EX_MOD-EX_SR(RS,LDASCREEN*AUTOA,IFLG))*(1-LDAX)
         ELSE IF (LUSE_MODEL_HF) THEN
         ! model sX 
            EX_MOD=(EX(RS,IFLG,.FALSE.)-EX_SR(RS,LDASCREEN*AUTOA,IFLG))*LDAX
         ELSE 
         ! WPBE: LDAX*EX_SR + EX_LR = EX - (1-LDAX)*EX_SR
            EX_MOD=EX(RS,IFLG,.FALSE.)-EX_SR(RS,LDASCREEN*AUTOA,IFLG)*(1-LDAX)
         ENDIF
      ENDIF
    END FUNCTION EX_MOD

    FUNCTION VX_MOD(RS,IFLG,TREL)
      USE xclib
      USE constant
      IMPLICIT NONE
      REAL(q) VX_MOD, RS
      INTEGER IFLG
      LOGICAL TREL

      VX_MOD=VX(RS,IFLG,TREL)
      
      IF (LDASCREEN/=0 .OR. LDAX/=1 ) THEN
         IF (LDASCREEN==0) THEN
            VX_MOD=VX_MOD*LDAX
         ELSE IF (LUSE_THOMAS_FERMI) THEN
            VX_MOD=VX_MOD-VX_SX(RS,LDASCREEN*AUTOA,IFLG)*(1-LDAX)
         ELSE IF (LUSE_LONGRANGE_HF) THEN
         ! Iann Gerber: use VX_SR (for LDAX=0, i.e. AEXX=1)
            VX_MOD=VX(RS,IFLG,.FALSE.)-(VX_MOD-VX_SR(RS,LDASCREEN*AUTOA,IFLG))*(1-LDAX)
         ELSE IF (LUSE_MODEL_HF) THEN
         ! model sX 
            VX_MOD=(VX(RS,IFLG,.FALSE.)-VX_SR(RS,LDASCREEN*AUTOA,IFLG))*LDAX
         ELSE
         ! WPBE: LDAX*VX_SR + VX_LR = VX - (1-LDAX)*VX_SR
            VX_MOD=VX(RS,IFLG,.FALSE.)-VX_SR(RS,LDASCREEN*AUTOA,IFLG)*(1-LDAX)
         ENDIF
      ENDIF
    END FUNCTION VX_MOD

!***********************************************************************
!
! EC_MOD calculate the correlation energy density per particle eps_c
! that is 
! E_x = \int eps_c(rho(r)) rho(r) d^3 r
!
! VC_MOD calculates the exchange potential, which is the derivative
! of eps_c(rho(r)) rho(r) with respect to rho(r)
! v_c = d (eps_c(rho) rho) /d rho
! 
!
!***********************************************************************
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
! Iann Gerber 18/11/04
! Comments Joachim Paier jP

! Splitting the correlation part: meaning only take into account of the 
! short-range part if necessary
      FUNCTION EC_MOD(RS,ZETA,IFLG,TREL)
      USE xclib
      USE constant
      IMPLICIT NONE
      REAL(q) EC_MOD, RS
      REAL(q) ZETA              ! relative spin-polarization zeta := (nu-nd)/(nu+nd)
      REAL(q) EC_MOD_TMP
      INTEGER IFLG
      LOGICAL TREL

      EC_MOD=ECVO(RS,IFLG)
      IF (LRANGE_SEPARATED_CORR) THEN
         IF (LDASCREENC==0) THEN
            EC_MOD=EC_MOD
! jP here, we ask for the short-range LSDA correlation energy
         ELSE
            CALL ECORRSR(RS,ZETA,LDASCREENC*AUTOA,EC_MOD)
! jP: Routines of Paola Gori-Giorig written in a.u.
! jP: conversion to Rydberg units necessary for VASP
            EC_MOD=2.0_q*EC_MOD
            IF (LUSE_LONGRANGE_HF) THEN
! keep short range part of LDA correlation energy
!               EC_MOD=ECVO(RS,1)
!               EC_MOD=EC_SR(RS,LDASCREENC*AUTOA,EC_MOD)
               CALL ECORRSR(RS,ZETA,LDASCREENC*AUTOA,EC_MOD)
! jP: Routines of Paola Gori-Giorig written in a.u.
! jP: conversion to Rydberg units necessary for VASP
               EC_MOD=2.0_q*EC_MOD
! keep long range part of the LDA correlation energy
!                  EC_MOD=ECVO(RS,1)
!                  EC_MOD=EC_MOD-EC_SR(RS,LDASCREENC*AUTOA,EC_MOD)  
            ELSE
               CALL ECORRLR(RS,ZETA,LDASCREENC*AUTOA,EC_MOD)
! jP: Routines of Paola Gori-Giorig written in a.u.
! jP: conversion to Rydberg units necessary for VASP
               EC_MOD=2.0_q*EC_MOD            
            ENDIF
         ENDIF
      ENDIF
      END FUNCTION EC_MOD

      FUNCTION VC_MOD(RS,ZETA,IFLG,TREL)
      USE xclib
      USE constant
      IMPLICIT NONE
      REAL(q) VC_MOD, RS
      REAL(q)  ZETA, VCSRUP, VCSRDOWN, VCLRUP, VCLRDOWN
      INTEGER IFLG
      LOGICAL TREL
      
      VC_MOD=VCVO(RS,IFLG)
      IF (LRANGE_SEPARATED_CORR) THEN
         IF (LDASCREENC==0) THEN
            VC_MOD=VC_MOD
         ELSE
            CALL VCORRSR(RS,ZETA,LDASCREENC*AUTOA,VCSRUP,VCSRDOWN)
! jP: Routines of Paola Gori-Giorig written in a.u.
! jP: conversion to Rydberg units necessary for VASP
            VCSRUP=2.0_q*VCSRUP
            VCSRDOWN=2.0_q*VCSRDOWN
            IF (LUSE_LONGRANGE_HF) THEN
! keep short range part of LDA correlation potential
!               VC_MOD=VC_SR(RS,LDASCREENC*AUTOA,VCVO(RS,IFLG),ECVO(RS,IFLG))
               CALL VCORRSR(RS,ZETA,LDASCREENC*AUTOA,VCSRUP,VCSRDOWN)
! jP: Routines of Paola Gori-Giorig written in a.u.
! jP: conversion to Rydberg units necessary for VASP
               VCSRUP=2.0_q*VCSRUP
               VCSRDOWN=2.0_q*VCSRDOWN
            ELSE
! keep long range part of the LDA correlation potential
!              VC_MOD=VC_MOD-VC_SR(RS,LDASCREENC*AUTOA,VCVO(RS,IFLG),ECVO(RS,IFLG))
               CALL VCORRLR(RS,ZETA,LDASCREENC*AUTOA,VCLRUP,VCLRDOWN)
! jP: Routines of Paola Gori-Giorig written in a.u.
! jP: conversion to Rydberg units necessary for VASP
               VCSRUP=2.0_q*VCSRUP
               VCSRDOWN=2.0_q*VCSRDOWN
            ENDIF
         ENDIF
      ENDIF
      END FUNCTION VC_MOD
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------

!***********************************************************************
!
!  VASP interpolates the XC-energy density from a table (at least
!    the plane wave part)
!  the required table is generated here
!
!***********************************************************************

    SUBROUTINE SETUP_LDA_XC(ISPIN,IU6,IU0,IDIOT)
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE tutor, ONLY: vtutor
      USE ini
      IMPLICIT NONE

      INTEGER  ISPIN            ! spin

! arrays for tutor call
      INTEGER IU6,IU0,IDIOT
! temporary
      INTEGER N,NDUMMY
      REAL(q) AMARG
      CHARACTER (1) CSEL
      CHARACTER (2) CEXCH
      IF (LEXCH<0) THEN
         CALL vtutor%bug("internal ERROR in SETUP_LDA_XC: LEXCH has not been set up", __FILE__, __LINE__)
      ENDIF
!
! set the exchange correlation type for the internal table
! before the table is used XCTABLE_CHECK should be called
! to ascertain that LEXCH_TABLE is correct and equivalent to LEXCH
!

      CALL SET_EX_TABLE(LEXCH,NEXCH,EXCTAB%EXCTAB,EXCTAB%NEXCHF,EXCTAB%RHOEXC,IU0)

      AMARG=1E30_q  ! natural boundary conditions required
      CALL SPLCOF(EXCTAB%EXCTAB(1,1,1),EXCTAB%NEXCHF(2),NEXCH,AMARG)

      IF (ISPIN==2) THEN
         CALL SPLCOF(EXCTAB%EXCTAB(1,1,2),EXCTAB%NEXCHF(2),NEXCH,AMARG)
         CALL SPLCOF(EXCTAB%EXCTAB(1,1,3),EXCTAB%NEXCHF(2),NEXCH,AMARG)
         CALL SPLCOF(EXCTAB%EXCTAB(1,1,4),EXCTAB%NEXCHF(2),NEXCH,AMARG)
         CALL SPLCOF(EXCTAB%EXCTAB(1,1,5),EXCTAB%NEXCHF(2),NEXCH,AMARG)
         CALL SPLCOF(EXCTAB%EXCTAB(1,1,6),EXCTAB%NEXCHF(2),NEXCH,AMARG)
      ENDIF
      LEXCH_TABLE=LEXCH

      IF (IU6>=0) WRITE(IU6,7004)LEXCH,EXCTAB%RHOEXC(1),EXCTAB%NEXCHF(1), &
     &                     EXCTAB%RHOEXC(2),EXCTAB%NEXCHF(2)
 7004 FORMAT(' exchange correlation table for  LEXCH = ',I8/ &
     &       '   RHO(1)= ',F8.3,5X,'  N(1)  = ',I8/ &
     &       '   RHO(2)= ',F8.3,5X,'  N(2)  = ',I8)

!$ACC UPDATE DEVICE(EXCTAB) ASYNC(ACC_ASYNC_Q)
!$ACC UPDATE DEVICE(EXCTAB%RHOEXC,EXCTAB%NEXCHF,EXCTAB%EXCTAB) ASYNC(ACC_ASYNC_Q)
      RETURN
    END SUBROUTINE 


!************************ PROGRAM SET_EX_TABLE *************************
!
! set up the default xc-table
! i.e.
!  Ceperly Alder
!  with standard interpolation to spin
!  with relativistic correction
!
!***********************************************************************

    SUBROUTINE SET_EX_TABLE(LEXCH,N,EXCTAB,NEXCHF,RHOEXC,IU0)

      USE xclib
      USE string, ONLY: str
      USE tutor, ONLY: vtutor

      INTEGER LEXCH, N, IU0
      REAL(q) :: EXCTAB(N,5,6),RHOEXC(2)
      INTEGER :: NEXCHF(2)
    ! local
      LOGICAL TREL
      CHARACTER (13) CEXCH
      INTEGER LEXCH_LDA
      INTEGER J,I
      REAL(q) :: SLATER, RHOSMA, RHOMAX, RH, EXCP, DEXF, DECF, ALPHA, DRHO, ZETA, FZA, FZB

      TREL=.TRUE.
! convert LEXCH to local format
      IF (LEXCH==0) THEN
        LEXCH_LDA=0
      ELSE IF (LEXCH==1) THEN
    ! Hedin Lundquist
        LEXCH_LDA=4
      ELSE IF (LEXCH==2) THEN
    ! Ceperly-Alder
        LEXCH_LDA=1
      ELSE IF (LEXCH==3) THEN
    ! Wigner
        LEXCH_LDA=6
!vdw jk
      ELSE IF (LEXCH==8 .OR. LEXCH==9 .OR. (LEXCH.ge.40 .and. LEXCH.le.50)) THEN
!vdw jk
    ! Pade approximation of Perdew
        LEXCH_LDA=7
      ELSE IF (LEXCH==10) THEN
    ! Iann Gerber: LDA exchange + VWN correlation
        LEXCH_LDA=2
      ELSE IF (LEXCH==11) THEN
    ! Joachim Paier: B3LYP LDA part is VWNIII-correlation
        LEXCH_LDA=11
      ELSE IF (LEXCH==12) THEN
    ! Joachim Paier: B3LYP LDA part is VWN5-correlation
        LEXCH_LDA=12
!aem AM05 uses the PW LDA correlation (same as PBE)
      ELSE IF (LEXCH==13) THEN
        LEXCH_LDA=7
!aem end addition
! jP: PBEsol - LDA part identical to PBE
      ELSE IF (LEXCH==14) THEN
        LEXCH_LDA=7
! jP: PBEsol - LDA part identical to PBE
!BEEF xc
#ifdef libbeef
      ELSE IF (LEXCH==17) THEN
        LEXCH_LDA=7
#endif
      ELSE IF (LEXCH==100) THEN
        LEXCH_LDA=-1
! jH-new RPA Perdew Wang 
      ELSE IF (LEXCH==20) THEN
       LEXCH_LDA=20
! functionals for range-separated ACFDT (LDA - short range RPA)
      ELSE IF (LEXCH==21) THEN
        ! \mu = 0.3 A^-1
        LEXCH_LDA=21
      ELSE IF (LEXCH==22) THEN
        ! \mu = 0.3 A^-1
        LEXCH_LDA=22
      ELSE IF (LEXCH==23) THEN
        ! \mu = 0.3 A^-1
        LEXCH_LDA=23
      ELSE IF (LEXCH==24) THEN
        ! \mu = 0.3 A^-1
        LEXCH_LDA=24
! jH-new RPAplus Perdew Wang
      ELSE IF (LEXCH==30) THEN
        LEXCH_LDA=30
#ifdef USELIBXC
      ELSE IF (LEXCH==99) THEN
        LEXCH_LDA=-1
#endif
      ELSE
        LEXCH_LDA=1
      ENDIF

      IF (LEXCH_LDA==0) THEN
         CEXCH='  '
      ELSE IF (LEXCH_LDA==1) THEN
         CEXCH='Ceperly-Alder'
      ELSE IF (LEXCH_LDA==2) THEN
         CEXCH='VW'
      ELSE IF (LEXCH_LDA==3) THEN
         CEXCH='GL'
      ELSE IF (LEXCH_LDA==4) THEN
         CEXCH='HL'
      ELSE IF (LEXCH_LDA==5) THEN
         CEXCH='BH'
      ELSE IF (LEXCH_LDA==6) THEN
         CEXCH='WI'
      ELSE IF (LEXCH_LDA==7) THEN
         CEXCH='PB'
      ELSE IF (LEXCH_LDA==-1) THEN
         CEXCH='no ex-corr'
      ELSE IF (LEXCH_LDA==11) THEN
! Joachim Paier
         CEXCH='VWN3'
      ELSE IF (LEXCH_LDA==12) THEN
! Joachim Paier
         CEXCH='VWN5'
! jH-new RPA Perdew Wang 
      ELSE IF (LEXCH_LDA==20) THEN 
       CEXCH='RPA Perdew-Wang'
! functionals for range-separated ACFDT (LDA - short range RPA)
      ELSE IF (LEXCH_LDA==21) THEN
        ! \mu = 0.3 A^-1
        CEXCH='0.3'
      ELSE IF (LEXCH_LDA==22) THEN
        ! \mu = 0.3 A^-1
        CEXCH='0.5'
      ELSE IF (LEXCH_LDA==23) THEN
        ! \mu = 0.3 A^-1
        CEXCH='1.0'
      ELSE IF (LEXCH_LDA==24) THEN
        ! \mu = 0.3 A^-1
        CEXCH='2.0'
! jH-new RPAplus Perdew Wang 
      ELSE IF (LEXCH_LDA==30) THEN
        CEXCH='RPA+'
      ELSE
         CALL vtutor%bug("internal error in SET_EX_TABLE: Wrong exchange correlation type! " // &
            str(LEXCH_LDA), __FILE__, __LINE__)
      ENDIF
      

      IF (IU0>=0) THEN
         IF (LEXCH_LDA==7) THEN
            WRITE(IU0,*)'LDA part: xc-table for Pade appr. of Perdew'
         ELSE IF (LEXCH_LDA>20 .AND. LEXCH_LDA<25) THEN
            WRITE(IU0,'(A,A3,A)')'LDA part: xc-table for (LDA - short range RPA): hard-coded for mu= ',CEXCH,' A^-1'
         ELSE
            IF (LFCI==1) THEN 
               WRITE(IU0,*)'LDA part: xc-table for ',CEXCH, ', Vosko type interpolation para-ferro'
            ELSE
               WRITE(IU0,*)'LDA part: xc-table for ',CEXCH, ', standard interpolation'
            ENDIF
         ENDIF
      ENDIF

! Slater parameter
      SLATER=1._q
! standard interpolation   from para- to ferromagnetic corr
      RHOSMA=INT(0.5_q*1000)/1000._q
!#define hugeXCtable
#ifdef hugeXCtable
      RHOMAX=INT(2000.5_q*1000)/1000._q
      WRITE(0,*) '******* RHOMAX set to 2000 instead of 100 *****'
#else
      RHOMAX=INT(100.5_q*1000)/1000._q
#endif


      RHOEXC(1)=RHOSMA
      RHOEXC(2)=RHOMAX
      NEXCHF(1)=NSMA
      NEXCHF(2)=N

      IF (NSMA/=0) THEN
         RH=RHOSMA/NSMA/2
      ELSE
         RH=RHOMAX/N/100
      ENDIF
      J=1

      CALL EXCHG(LEXCH_LDA,RH,EXCP, &
     &           DEXF,DECF,ALPHA,SLATER,TREL)
        EXCTAB(J,1,1)=RH
        EXCTAB(J,1,2)=RH
        EXCTAB(J,1,3)=RH
        EXCTAB(J,1,4)=RH

        EXCTAB(J,2,1)=EXCP
        EXCTAB(J,2,2)=DEXF
        EXCTAB(J,2,3)=DECF
        EXCTAB(J,2,4)=ALPHA


      IF (NSMA/=0) THEN
         DRHO=RHOSMA/NSMA
         DO 100 I=1,NSMA-1
            J=I+1
            RH=DRHO*I
            CALL EXCHG(LEXCH_LDA,RH,EXCP, &
     &                 DEXF,DECF,ALPHA,SLATER,TREL)
          EXCTAB(J,1,1)=RH
          EXCTAB(J,1,2)=RH
          EXCTAB(J,1,3)=RH
          EXCTAB(J,1,4)=RH

          EXCTAB(J,2,1)=EXCP
          EXCTAB(J,2,2)=DEXF
          EXCTAB(J,2,3)=DECF
          EXCTAB(J,2,4)=ALPHA

  100    CONTINUE

         J=NSMA+1
         RH=RHOSMA
         CALL EXCHG(LEXCH_LDA,RH,EXCP, &
     &              DEXF,DECF,ALPHA,SLATER,TREL)
         EXCTAB(J,1,1)=RH
         EXCTAB(J,1,2)=RH
         EXCTAB(J,1,3)=RH
         EXCTAB(J,1,4)=RH

         EXCTAB(J,2,1)=EXCP
         EXCTAB(J,2,2)=DEXF
         EXCTAB(J,2,3)=DECF
         EXCTAB(J,2,4)=ALPHA
      ENDIF

      DRHO=(RHOMAX-RHOSMA)/(N-NSMA)
      DO 200 I=1,N-NSMA-1
         J=I+NSMA+1
         RH=DRHO*I+RHOSMA
         CALL EXCHG(LEXCH_LDA,RH,EXCP, &
     &              DEXF,DECF,ALPHA,SLATER,TREL)
         EXCTAB(J,1,1)=RH
         EXCTAB(J,1,2)=RH
         EXCTAB(J,1,3)=RH
         EXCTAB(J,1,4)=RH

         EXCTAB(J,2,1)=EXCP
         EXCTAB(J,2,2)=DEXF
         EXCTAB(J,2,3)=DECF
         EXCTAB(J,2,4)=ALPHA
  200 CONTINUE

      J=1
      ZETA=0
      FZA =.854960467080682810_q
      FZB =.854960467080682810_q
      EXCTAB(J,1,5)=ZETA
      EXCTAB(J,1,6)=ZETA
      EXCTAB(J,2,5)=FZA
      EXCTAB(J,2,6)=FZB

      DO 800 I=1,N-1
         ZETA=FLOAT(I)/FLOAT(N-1)
         FZA=FZ0(ZETA)/ZETA/ZETA
         FZB=FZ0(ZETA)/ZETA/ZETA
         J=I+1
         EXCTAB(J,1,5)=ZETA
         EXCTAB(J,1,6)=ZETA
         EXCTAB(J,2,5)=FZA
         EXCTAB(J,2,6)=FZB
  800 CONTINUE

      IF (.FALSE.) THEN
      DO 300 I=1,N
         WRITE(97,20) EXCTAB(I,1,1),EXCTAB(I,2,1)
  300 CONTINUE
      DO 400 I=1,N
         WRITE(97,20)  EXCTAB(I,1,2),EXCTAB(I,2,2)
  400 CONTINUE
      DO 500 I=1,N
         WRITE(97,20)  EXCTAB(I,1,3),EXCTAB(I,2,3)
  500 CONTINUE
      DO 600 I=1,N
         WRITE(97,20)  EXCTAB(I,1,4),EXCTAB(I,2,4)
  600 CONTINUE
      DO 700 I=1,N
         WRITE(97,20)  EXCTAB(I,1,5),EXCTAB(I,2,5)
  700 CONTINUE
      DO 710 I=1,N
         WRITE(97,20)  EXCTAB(I,1,6),EXCTAB(I,2,6)
  710 CONTINUE
   20 FORMAT((3(E24.16,2X)))
      ENDIF

      RETURN
   END SUBROUTINE SET_EX_TABLE


!**************** SUBROUTINE EXCHG *************************************
!
! EXCHG calculate the LDA part of xc energy density per particle eps_xc(rho)
! 
! E_xc = \int eps_xc(rho(r)) rho(r) d^3 r

! The subroutines calls various subroutines defined in xclib
! The function should be called with care, since the flag LEXCH_LDA
! is not compatible to LEXCH
!
! Before returning, the subroutine divides the results by the density^(1/3)
! which is the behaviour of the exchange density per particle e_x 
!
! this routine expects input in Angst, and the output is in eV/Angst
!
! Most of the called routines expect atomic units, hence
! EXCHG calculates the Wigner Seitz radius in a.u. and then calls
! the various exchange correlation routines (which are expected
! to return the exchange correlation energy as well as potential 
! in Rydberg)
!
!***********************************************************************

   SUBROUTINE EXCHG(LEXCH_LDA,RHO,EXCP,DEXF,DECF,ALPH,SLATER,TREL)
      USE constant
      USE xclib

      INTEGER LEXCH_LDA

      LOGICAL TREL
      REAL(q) :: RHO,EXCP,DEXF,DECF,ALPH,SLATER
      REAL(q) :: RHOP, RHOM, EXP, EXM, DELTA
   ! local
      REAL(q) :: RS, RH, ZETA, FZA, FZB, &
           RHOTHD, SE, ECLDA, ECD1LDA, ECD2LDA, ECLDA_MAG, A0

      IF (RHO==0) THEN
         EXCP=0._q
         DEXF=0._q
         DECF=0._q
         ALPH=0._q
         RETURN
      ENDIF
#ifdef test_exchange_potentials
! this small DO loop allows to test the consistency between 
! the exchange energy density per particle EX_MOD
! and the corresponding potential VX_MOD
      RHO=1.0
      DO
         DELTA=0.0001_q
         RHOP=RHO+DELTA
         RHOM=RHO-DELTA

         RHOTHD = RHOP**(1/3._q)
         RS = (3._q/(4._q*PI)/RHOP)**(1/3._q) /AUTOA
         EXP=EX_MOD(RS, 1, .FALSE.)

         RHOTHD = RHOM**(1/3._q)
         RS = (3._q/(4._q*PI)/RHOM)**(1/3._q) /AUTOA
         EXM=EX_MOD(RS, 1, .FALSE.)

         RHOTHD = RHO**(1/3._q)
         RS = (3._q/(4._q*PI)/RHO)**(1/3._q) /AUTOA

         WRITE(*,'(3F14.7)') RHO, (EXP*RHOP/AUTOA**3-EXM*RHOM/AUTOA**3)/2.0_q/(DELTA/AUTOA**3), VX_MOD(RS,1,.FALSE.)

         RHO=RHO+0.001
      ENDDO
#endif

      RHOTHD = RHO**(1/3._q)
      RS = (3._q/(4._q*PI)/RHO)**(1/3._q) /AUTOA
      IF (LEXCH_LDA==0) THEN
         EXCP=SLATER*EX_MOD(RS,1,TREL)
         DEXF=SLATER*EX_MOD(RS,2,TREL)-EXCP
         DECF=0._q
         ALPH=0._q
      ELSE IF (LEXCH_LDA==1) THEN
         EXCP=EX_MOD(RS,1,TREL)+ECCA(RS,1)*ALDAC
         DEXF=EX_MOD(RS,2,TREL)-EX_MOD(RS,1,TREL)
         DECF=ECCA(RS,2)*ALDAC-ECCA(RS,1)*ALDAC
         ALPH=0._q
      ELSE IF (LEXCH_LDA==2) THEN
! Iann Gerber 18/11/04
!        EXCP=EX_MOD(RS,1,TREL)+ECVO(RS,1)
         ZETA=0.0_q
         EXCP=EX_MOD(RS,1,TREL)+EC_MOD(RS,ZETA,1,TREL)*ALDAC
         DEXF=EX_MOD(RS,2,TREL)-EX_MOD(RS,1,TREL)
!        DECF=ECVO(RS,2)*ALDAC-ECVO(RS,1)*ALDAC
! hack: ZETA set to 1.0_q or 0.0_q, respectively
         DECF=EC_MOD(RS,1._q,2,TREL)*ALDAC-EC_MOD(RS,0._q,1,TREL)*ALDAC
         ALPH=0._q
! Iann Gerber 18/11/04; updated to scrLSDA by Joachim Paier (28/04/08)
      ELSE IF (LEXCH_LDA==3) THEN
         EXCP=EX_MOD(RS,1,TREL)+ECGL(RS,1)*ALDAC
         DEXF=EX_MOD(RS,2,TREL)-EX_MOD(RS,1,TREL)
         DECF=ECGL(RS,2)*ALDAC-ECGL(RS,1)*ALDAC
         ALPH=0._q
      ELSE IF (LEXCH_LDA==4) THEN
         EXCP=EX_MOD(RS,1,TREL)+ECHL(RS,1)*ALDAC
         DEXF=EX_MOD(RS,2,TREL)-EX_MOD(RS,1,TREL)
         DECF=ECHL(RS,2)*ALDAC-ECHL(RS,1)*ALDAC
         ALPH=0._q
      ELSE IF (LEXCH_LDA==5) THEN
         EXCP=EX_MOD(RS,1,TREL)+ECBH(RS,1)*ALDAC
         DEXF=EX_MOD(RS,2,TREL)-EX_MOD(RS,1,TREL)
         DECF=ECBH(RS,2)*ALDAC-ECBH(RS,1)*ALDAC
         ALPH=0._q
      ELSE IF (LEXCH_LDA==6) THEN
         EXCP=EX_MOD(RS,1,TREL)+ECWI(RS,1)*ALDAC
         DEXF=EX_MOD(RS,2,TREL)-EX_MOD(RS,1,TREL)
         DECF=ECWI(RS,2)*ALDAC-ECWI(RS,1)*ALDAC
         ALPH=0._q
      ELSE IF (LEXCH_LDA==7) THEN
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA)
         ZETA=1  ! ferromagnetic result
         CALL CORPBE_LDA(RS,ZETA,ECLDA_MAG,ECD1LDA,ECD2LDA)

         EXCP=EX_MOD(RS,1,.FALSE.)+ECLDA*ALDAC
         DEXF=EX_MOD(RS,2,.FALSE.)-EX_MOD(RS,1,.FALSE.)
         DECF=ECLDA_MAG*ALDAC-ECLDA*ALDAC
!        WRITE(77,'(5F14.7)') RS, ECCA(RS,1),ECLDA,ECCA(RS,2),ECLDA_MAG
         ALPH=0._q
      ELSE IF (LEXCH_LDA==11) THEN
! Joachim Paier
         EXCP=EX_MOD(RS,1,TREL)+ECVOIII(RS,1)*ALDAC
         DEXF=EX_MOD(RS,2,TREL)-EX_MOD(RS,1,TREL)
         DECF=ECVOIII(RS,2)*ALDAC-ECVOIII(RS,1)*ALDAC
         ALPH=0._q
! check is VO similar to CA
!        WRITE(77,'(5F14.7)') RS, ECCA(RS,2),ECVO(RS,2)
!        WRITE(78,'(5F14.7)') RS, ECCA(RS,1),ECVO(RS,1)
      ELSE IF (LEXCH_LDA==12) THEN
! Joachim Paier
         EXCP=EX_MOD(RS,1,TREL)+ECVO(RS,1)*ALDAC
         DEXF=EX_MOD(RS,2,TREL)-EX_MOD(RS,1,TREL)
         DECF=ECVO(RS,2)*ALDAC-ECVO(RS,1)*ALDAC
         ALPH=0._q
! check is VO similar to CA
!        WRITE(77,'(5F14.7)') RS, ECCA(RS,2),ECVO(RS,2)
!        WRITE(78,'(5F14.7)') RS, ECCA(RS,1),ECVO(RS,1)
! jH-new RPA with Perdew Wang para
      ELSE IF (LEXCH_LDA==20) THEN
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_RPA(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA)
         ZETA=1  ! ferromagnetic result
         CALL CORPBE_LDA_RPA(RS,ZETA,ECLDA_MAG,ECD1LDA,ECD2LDA)
         EXCP=EX_MOD(RS,1,.FALSE.)+ECLDA*ALDAC
         DEXF=EX_MOD(RS,2,.FALSE.)-EX_MOD(RS,1,.FALSE.)
         DECF=ECLDA_MAG*ALDAC-ECLDA*ALDAC
         ALPH=0._q
! Judith Harl
! functionals for range-separated ACFDT (LDA - short range RPA):
! a bit akward at the moment since the range separation parameter
! is hard coded for now. Hopefully this will change ...
      ELSE IF (LEXCH_LDA==21) THEN
! jh --- fuer den spin-polarisierten Fall noch nicht moeglich
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_SR_0_15au(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA)
         ZETA=0  ! auch paramagnetisch (nicht fuer spin-polarisiert)
         CALL CORPBE_LDA_SR_0_15au(RS,ZETA,ECLDA_MAG,ECD1LDA,ECD2LDA)
         EXCP=EX_MOD(RS,1,.FALSE.)+ECLDA*ALDAC
         DEXF=EX_MOD(RS,2,.FALSE.)-EX_MOD(RS,1,.FALSE.)
         DECF=ECLDA_MAG*ALDAC-ECLDA*ALDAC
         ALPH=0._q
      ELSE IF (LEXCH_LDA==22) THEN
! jh --- fuer den spin-polarisierten Fall noch nicht moeglich
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_SR_0_5A(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA)
         ZETA=0  ! auch paramagnetisch (nicht fuer spin-polarisiert)
         CALL CORPBE_LDA_SR_0_5A(RS,ZETA,ECLDA_MAG,ECD1LDA,ECD2LDA)
         EXCP=EX_MOD(RS,1,.FALSE.)+ECLDA*ALDAC
         DEXF=EX_MOD(RS,2,.FALSE.)-EX_MOD(RS,1,.FALSE.)
         DECF=ECLDA_MAG*ALDAC-ECLDA*ALDAC
         ALPH=0._q
      ELSE IF (LEXCH_LDA==23) THEN
! jh --- fuer den spin-polarisierten Fall noch nicht moeglich
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_SR_1_0A(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA)
         ZETA=0  ! auch paramagnetisch (nicht fuer spin-polarisiert)
         CALL CORPBE_LDA_SR_1_0A(RS,ZETA,ECLDA_MAG,ECD1LDA,ECD2LDA)
         EXCP=EX_MOD(RS,1,.FALSE.)+ECLDA*ALDAC
         DEXF=EX_MOD(RS,2,.FALSE.)-EX_MOD(RS,1,.FALSE.)
         DECF=ECLDA_MAG*ALDAC-ECLDA*ALDAC
         ALPH=0._q
      ELSE IF (LEXCH_LDA==24) THEN
! jh --- fuer den spin-polarisierten Fall noch nicht moeglich
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_SR_2_0A(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA)
         ZETA=0  ! auch paramagnetisch (nicht fuer spin-polarisiert)
         CALL CORPBE_LDA_SR_2_0A(RS,ZETA,ECLDA_MAG,ECD1LDA,ECD2LDA)
         EXCP=EX_MOD(RS,1,.FALSE.)+ECLDA*ALDAC
         DEXF=EX_MOD(RS,2,.FALSE.)-EX_MOD(RS,1,.FALSE.)
         DECF=ECLDA_MAG*ALDAC-ECLDA*ALDAC
         ALPH=0._q         
! jH-new RPA-plus = ELDAQMC-ELDARPA, exchange term is zero 
      ELSE IF (LEXCH_LDA==30) THEN
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_RPA_PLUS(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA)
         ZETA=1  ! ferromagnetic result
         CALL CORPBE_LDA_RPA_PLUS(RS,ZETA,ECLDA_MAG,ECD1LDA,ECD2LDA)
         EXCP=EX_MOD(RS,1,.FALSE.)+ECLDA*ALDAC
         DEXF=EX_MOD(RS,2,.FALSE.)-EX_MOD(RS,1,.FALSE.)
         DECF=ECLDA_MAG*ALDAC-ECLDA*ALDAC
         ALPH=0._q
      ELSE
         EXCP=0._q
         DEXF=0._q
         DECF=0._q
         ALPH=0._q
      ENDIF
      !
      ! for Perdew, Burke Ernzerhof we always use the 
      ! recommended  interpolation from nm to magnetic
      ! 
      IF (LEXCH_LDA==11) THEN
         A0=ALPHA0_III(RS)*ALDAC
         ALPH=DECF-A0
         DECF=A0
      ELSE IF (LEXCH_LDA==7) THEN
         A0=PBE_ALPHA(RS)*ALDAC
         ALPH=DECF-A0
         DECF=A0
! jH-new new switching for RPA 
      ELSE IF (LEXCH_LDA==20) THEN
         A0=RPA_ALPHA(RS)*ALDAC
         ALPH=DECF-A0
         DECF=A0
! jH-new new switching for RPAplus
      ELSE IF (LEXCH_LDA==30) THEN          
         A0=(PBE_ALPHA(RS)-RPA_ALPHA(RS))*ALDAC
         ALPH=DECF-A0
         DECF=A0
      ELSE IF (LFCI==1 .OR. LEXCH_LDA==12 ) THEN
         A0=ALPHA0(RS)*ALDAC
         ALPH=DECF-A0
         DECF=A0
      ENDIF

      EXCP=EXCP*RYTOEV/RHOTHD
      DEXF=DEXF*RYTOEV/RHOTHD
      DECF=DECF*RYTOEV/RHOTHD
      ALPH=ALPH*RYTOEV/RHOTHD

      RETURN
    END SUBROUTINE EXCHG

!************************ SUBROUTINE SET_LEXCH *************************
!
! in principle this is the only function that should be used
! to set the exchange correlation type
! however, presently the fock.F and pseudo.F modules access LEXCH
! directly
!
!***********************************************************************

    SUBROUTINE SET_LEXCH(LEXCH_SET)
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      IMPLICIT NONE
      INTEGER LEXCH_SET

      LEXCH=LEXCH_SET
!$ACC UPDATE DEVICE(LEXCH) ASYNC(ACC_ASYNC_Q)
    END SUBROUTINE SET_LEXCH


    SUBROUTINE PUSH_LEXCH(LEXCH_TEMP)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      IMPLICIT NONE
      INTEGER LEXCH_TEMP

      ISTACK=ISTACK+1
      IF (ISTACK>SIZE(EX_STACK,2)) THEN
         CALL vtutor%bug("internal ERROR in PUSH_LEXCH: push already exhausted " // str(ISTACK), __FILE__, __LINE__)
      ENDIF
      EX_STACK(1,ISTACK)=LEXCH
      LEXCH=LEXCH_TEMP
!$ACC UPDATE DEVICE(LEXCH) ASYNC(ACC_ASYNC_Q)
    END SUBROUTINE PUSH_LEXCH

    SUBROUTINE POP_LEXCH
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE tutor, ONLY: vtutor
      IF (ISTACK==0) THEN
         CALL vtutor%bug("internal ERROR in POP_LEXCH: push was not used", __FILE__, __LINE__)
      ENDIF
      LEXCH=EX_STACK(1,ISTACK)
      ISTACK=ISTACK-1
!$ACC UPDATE DEVICE(LEXCH) ASYNC(ACC_ASYNC_Q)
    END SUBROUTINE POP_LEXCH

!************************ SUBROUTINE SET_LEXCH *************************
!
! subroutine to temporarily use an alternative local 
! exchange correlation functional
!
!***********************************************************************

    SUBROUTINE PUSH_XC_TYPE(LEXCH_TEMP, LDAX_TEMP, LDAC_TEMP, GGAX_TEMP, GGAC_TEMP, LDASCREEN_TEMP)
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      USE fock_glb, ONLY: HFSCREEN, HFSCREENC, AEXX, L_MODEL_HF
      USE wpbe
      IMPLICIT NONE
      INTEGER LEXCH_TEMP
      REAL(q) LDAX_TEMP, LDAC_TEMP, LDASCREEN_TEMP, GGAX_TEMP, GGAC_TEMP

      ISTACK=ISTACK+1
      IF (ISTACK>SIZE(EX_STACK,2)) THEN
         CALL vtutor%bug("internal ERROR in PUSH_XC_TYPE: push already exhausted " // str(ISTACK), __FILE__, __LINE__)
      ENDIF
      ! gK: a comment is in place here
      ! if LDASCREEN is set to 0, xclib.F falls back to the non-range seperated functional, regardless of LUSE_MODEL_HF
      ! if this were not the case, one would need to set LUSE_MODEL_HF to .FALSE. 
      ! to recover the original PBE functional (see PUSH in paw.F)
      EX_STACK (1,ISTACK)=LEXCH
      EX_STACK (2,ISTACK)=LDAX
      EX_STACK (3,ISTACK)=LDASCREEN
      EX_STACK (4,ISTACK)=LFCI
      EX_STACK (5,ISTACK)=AGGAX
      EX_STACK (6,ISTACK)=AGGAC
      EX_STACK (7,ISTACK)=ALDAC
      EX_LSTACK(1,ISTACK)=LRANGE_SEPARATED_CORR
      EX_LSTACK(2,ISTACK)=LUSE_LONGRANGE_HF
      EX_LSTACK(3,ISTACK)=LUSE_THOMAS_FERMI
      EX_LSTACK(4,ISTACK)=LUSE_MODEL_HF
!      WRITE(*,'(A,4F8.3,2L5)') 'push',LDASCREEN, HFSCREEN, 1-LDAX, AEXX, LUSE_MODEL_HF, L_MODEL_HF

      LEXCH    =LEXCH_TEMP
      LDAX     =LDAX_TEMP
      ALDAC    =LDAC_TEMP
      LDASCREEN=LDASCREEN_TEMP
      AGGAX    =GGAX_TEMP
      AGGAC    =GGAC_TEMP
!gK 18.05.2021 also set some variables in fock_glb; this is way cleaner then the tinkering done before
! this was previously manually set in some places, but in a fairly inconsistent manner
! is was working most of the time, since the calling routine avoided to evaluate the Fock exchange
      AEXX = 1-LDAX
      HFSCREEN=LDASCREEN

      CALL INIT_WPBE
!$ACC UPDATE DEVICE(LEXCH, &
!$ACC&              LDAX, &
!$ACC&              LDASCREEN, &
!$ACC&              LFCI, &
!$ACC&              AGGAX, &
!$ACC&              AGGAC, &
!$ACC&              ALDAC, &
!$ACC&              LRANGE_SEPARATED_CORR, &
!$ACC&              LUSE_LONGRANGE_HF, &
!$ACC&              LUSE_THOMAS_FERMI, &
!$ACC&              LUSE_MODEL_HF) ASYNC(ACC_ASYNC_Q)
    END SUBROUTINE PUSH_XC_TYPE


    SUBROUTINE POP_XC_TYPE
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE wpbe
      USE tutor, ONLY: vtutor
      USE fock_glb, ONLY: HFSCREEN, HFSCREENC, AEXX, L_MODEL_HF
      IMPLICIT NONE

      IF (ISTACK==0) THEN
         CALL vtutor%bug("internal ERROR in POP_XC_TYPE: push was not used", __FILE__, __LINE__)
      ENDIF
!      WRITE(*,'(A,4F8.3,2L5)') 'before pop',LDASCREEN, HFSCREEN, 1-LDAX, AEXX, LUSE_MODEL_HF, L_MODEL_HF
      LEXCH                =EX_STACK (1,ISTACK)
      LDAX                 =EX_STACK (2,ISTACK)
      LDASCREEN            =EX_STACK (3,ISTACK)
      LFCI                 =EX_STACK (4,ISTACK)
      AGGAX                =EX_STACK (5,ISTACK)
      AGGAC                =EX_STACK (6,ISTACK)
      ALDAC                =EX_STACK (7,ISTACK)
      LRANGE_SEPARATED_CORR=EX_LSTACK(1,ISTACK)
      LUSE_LONGRANGE_HF    =EX_LSTACK(2,ISTACK)
      LUSE_THOMAS_FERMI    =EX_LSTACK(3,ISTACK)
      LUSE_MODEL_HF        =EX_LSTACK(4,ISTACK)
!gK 18.05.2021 also set some variables in fock_glb; this is way cleaner then the tinkering done before
      AEXX = 1-LDAX
      HFSCREEN=LDASCREEN
!      WRITE(*,'(A,4F8.3,2L5)') 'after pop',LDASCREEN, HFSCREEN, 1-LDAX, AEXX, LUSE_MODEL_HF, L_MODEL_HF

      ISTACK=ISTACK-1

      CALL INIT_WPBE
!$ACC UPDATE DEVICE(LEXCH, &
!$ACC&              LDAX, &
!$ACC&              LDASCREEN, &
!$ACC&              LFCI, &
!$ACC&              AGGAX, &
!$ACC&              AGGAC, &
!$ACC&              ALDAC, &
!$ACC&              LRANGE_SEPARATED_CORR, &
!$ACC&              LUSE_LONGRANGE_HF, &
!$ACC&              LUSE_THOMAS_FERMI, &
!$ACC&              LUSE_MODEL_HF) ASYNC(ACC_ASYNC_Q)
    END SUBROUTINE POP_XC_TYPE


!*********************** FUNCTION ISGGA  *******************************
!
! function that returns .TRUE. if a GGA functional has been selected
!
!***********************************************************************

    FUNCTION ISGGA()
      IMPLICIT NONE
      LOGICAL ISGGA
! Iann Gerber special case for splitting correlation functional LDA-only
      IF (LEXCH==10) THEN
         ISGGA=.FALSE.
      ELSE IF (LEXCH>3.AND.LEXCH<100) THEN
         ISGGA=.TRUE.
      ELSE
         ISGGA=.FALSE.
      ENDIF
    END FUNCTION ISGGA
      
!************************ FUNCTION ISLDAXC  ****************************
!
! function that returns .TRUE. if a local or semilocal xc functional is
! used
!
!***********************************************************************

    FUNCTION ISLDAXC()
      IMPLICIT NONE
      LOGICAL ISLDAXC
      IF (LEXCH>=0) THEN
         ISLDAXC=.TRUE.
      ELSE
         ISLDAXC=.FALSE.
      ENDIF
    END FUNCTION ISLDAXC

!************************ SUBROUTINE XCTABLE_CHECK *********************
!
! procedure that checks whether the exchange correlation table
! corresponds to the currently used exchange correlation type
! LEXCH
!
!***********************************************************************

    SUBROUTINE XCTABLE_CHECK
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      
      IF (LEXCH /= LEXCH_TABLE) THEN
         CALL vtutor%bug("internal ERROR in XCTABLE_CHECK: the exchange table is not properly set " &
            // "up " // str(LEXCH) // " " // str(LEXCH_TABLE), __FILE__, __LINE__)
      ENDIF
      
    END SUBROUTINE XCTABLE_CHECK


!*******************************************************************
!
!  calculate the exchange correlation potential
!  on a radial grid and the first 3 derivatives
!
!*******************************************************************

    SUBROUTINE EXCOR_DER_PARA(RHO, NDER, EXCA, TREL)
      USE constant
      IMPLICIT NONE

      LOGICAL TREL
      INTEGER NDER               ! number of derivatives
      REAL(q) RHO
      REAL(q) EXCA(4),EXCA_(4)
      REAL(q) RHOT,EXC0,EXC2,EXCD0,EXCD2,EPS
      ! this is the best compromise for densities between 1000 and 0.1
      REAL(q),PARAMETER :: DELTA=1E-3_q

      CALL EXCOR_PARA(RHO,EXCA(1),EXCA(2),TREL)

      IF (NDER>0) THEN
         EPS=DELTA*RHO

         RHOT=RHO-EPS
         CALL EXCOR_PARA(RHOT,EXC0,EXCD0,TREL)

         RHOT=RHO+EPS
         CALL EXCOR_PARA(RHOT,EXC2,EXCD2,TREL)
         ! 1st and 2nd derivative of energy
         EXCA_(2)=(EXC2-EXC0)/ (2*EPS)
         EXCA_(3)=(EXC2+EXC0-2*EXCA(1))/ (EPS*EPS)

         ! 2nd and 3nd derivative of potential=
         ! 1st and 2nd derivative of energy
         EXCA(3)=(EXCD2-EXCD0)/ (2*EPS)
         EXCA(4)=(EXCD2+EXCD0-2*EXCA(2))/ (EPS*EPS)
         ! WRITE(*,'(5E14.7)') EXCA(2),EXCA_(2),EXCA(3),EXCA_(3),EXCA(4)
      ENDIF

    END SUBROUTINE

    SUBROUTINE EXCOR_PARA(RHO,EXC,DEXC,TREL)
      USE constant
      USE xclib
      IMPLICIT NONE
      REAL(q) RHO,RS,EXC,DEXC
      REAL(q) ECLDA,ECDLDA,ECD2LDA,SK,T,EC,ECD,ECDD,ZETA
      LOGICAL TREL

      RS = ( (3._q/(4*PI)) / RHO)**(1/3._q) /AUTOA
#ifndef USELIBXC
      IF (LEXCH==-1) THEN
#else
      IF (LEXCH==-1.OR.LEXCH==99) THEN
#endif
         EXC = 0
         DEXC= 0
      !
      ! LDA part: Pade approximation to Ceperly Alder results
      !
!aem AM05 functional (LEXCH==13) added, it uses the PW correlation in the LDA part.
! jP: adding PBEsol
!vdw jk
      ELSE IF (LEXCH==8.OR.LEXCH==9.OR.LEXCH==13 .OR. LEXCH==14 .or. LEXCH==17 .or. (LEXCH.ge.40 .and. LEXCH.le.50)) THEN
!vdw jk
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)

         EXC = (EX_MOD(RS,1,.FALSE.)+ECLDA*ALDAC)*RHO*RYTOEV
         DEXC= (VX_MOD(RS,1,.FALSE.)+ECDLDA*ALDAC)*RYTOEV
! jH-new LEXCH==20 CORPBE_LDA_RPA
      ELSE IF (LEXCH==20)  THEN
         ZETA=0
         CALL CORPBE_LDA_RPA(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)

         EXC = (EX_MOD(RS,1,.FALSE.)+ECLDA*ALDAC)*RHO*RYTOEV
         DEXC= (VX_MOD(RS,1,.FALSE.)+ECDLDA*ALDAC)*RYTOEV
! jH-new LEXCH==21 CORPBE_LDA_SR_0_15au
      ELSE IF (LEXCH==21)  THEN
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_SR_0_15au(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)

         EXC = (EX_MOD(RS,1,.FALSE.)+ECLDA*ALDAC)*RHO*RYTOEV 
         DEXC= (VX_MOD(RS,1,.FALSE.)+ECDLDA*ALDAC)*RYTOEV
! jH-new LEXCH==22 CORPBE_LDA_SR_0_5A
      ELSE IF (LEXCH==22)  THEN
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_SR_0_5A(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)

         EXC = (EX_MOD(RS,1,.FALSE.)+ECLDA*ALDAC)*RHO*RYTOEV
         DEXC= (VX_MOD(RS,1,.FALSE.)+ECDLDA*ALDAC)*RYTOEV
! jH-new LEXCH==23 CORPBE_LDA_SR_1_0A
      ELSE IF (LEXCH==23)  THEN
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_SR_1_0A(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)

         EXC = (EX_MOD(RS,1,.FALSE.)+ECLDA*ALDAC)*RHO*RYTOEV
         DEXC= (VX_MOD(RS,1,.FALSE.)+ECDLDA*ALDAC)*RYTOEV
! jH-new LEXCH==24 CORPBE_LDA_SR_2_0A
      ELSE IF (LEXCH==24)  THEN
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_SR_2_0A(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)

         EXC = (EX_MOD(RS,1,.FALSE.)+ECLDA*ALDAC)*RHO*RYTOEV
         DEXC= (VX_MOD(RS,1,.FALSE.)+ECDLDA*ALDAC)*RYTOEV
! jH-new LEXCH==30 CORPBE_LDA_RPA_PLUS
      ELSE IF (LEXCH==30)  THEN 
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA_RPA_PLUS(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)            
         EXC = (EX_MOD(RS,1,.FALSE.)+ECLDA*ALDAC)*RHO*RYTOEV
         DEXC= (VX_MOD(RS,1,.FALSE.)+ECDLDA*ALDAC)*RYTOEV
      !  
      ! LDA: range separated exchange + range separated correlation
      !   
      ELSEIF (LEXCH==10) THEN
!         EXC = (EX_MOD(RS,1,TREL)+EC_MOD(RS,1,TREL)*ALDAC)*RHO*RYTOEV
!         DEXC= (VX_MOD(RS,1,TREL)+VC_MOD(RS,1,TREL)*ALDAC)*RYTOEV      
! Relativistic correction have been removed  
! jP be careful ... 
! jP: what does EC_MOD do in that case - LEXCH==10  ...
         EXC = (EX_MOD(RS,1,.FALSE.)+EC_MOD(RS,0._q,1,.FALSE.)*ALDAC)*RHO*RYTOEV
         DEXC= (VX_MOD(RS,1,.FALSE.)+VC_MOD(RS,0._q,1,.FALSE.)*ALDAC)*RYTOEV
      ELSEIF (LEXCH==11) THEN
! Joachim Paier
         EXC = (EX_MOD(RS,1,TREL)+ECVOIII(RS,1)*ALDAC)*RHO*RYTOEV
         DEXC= (VX_MOD(RS,1,TREL)+VCVOIII(RS,1)*ALDAC)*RYTOEV
         ! ECLDA=ECVO(RS,1)
         ! ECDLDA=VCVO(RS,1)
         ! jP: compare with CALL CORPBE_LDA(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)
      ELSEIF (LEXCH==12) THEN
! Joachim Paier
         EXC = (EX_MOD(RS,1,TREL)+ECVO(RS,1)*ALDAC)*RHO*RYTOEV
         DEXC= (VX_MOD(RS,1,TREL)+VCVO(RS,1)*ALDAC)*RYTOEV
         ! ECLDA=ECVO(RS,1)
         ! ECDLDA=VCVO(RS,1)
         ! jP: compare with CALL CORPBE_LDA(RS,ZETA,ECLDA,ECDLDA,ECD2LDA)
      ELSE IF (LEXCH==100) THEN
         EXC=0
         CALL COHSM1(RHO, DEXC)
      !
      ! LDA part: Perdew Zungers interpolation
      !
      ELSE
         EXC = (EX_MOD(RS,1,TREL)+ECCA(RS,1)*ALDAC)*RHO*RYTOEV
         DEXC= (VX_MOD(RS,1,TREL)+VCCA(RS,1)*ALDAC)*RYTOEV
      ENDIF

    END SUBROUTINE

!*******************************************************************
!
!  calculate the exchange correlation energy density
!  on a radial grid and the first 3 derivatives
!  for the spinpolarized case
!  as input the density (up and down) is required
!  derivatives with respect to up and down components are calculated
!
!*******************************************************************
      
    SUBROUTINE EXCOR_DER(RHOUP,RHODOWN,NDER,EXC,EXCD,EXCDD,EXCDDD,TREL)
      USE constant
      IMPLICIT NONE

      LOGICAL TREL
      INTEGER NDER               ! number of derivative
      REAL(q) RHOUP,RHODOWN      ! up and down component of density
      REAL(q) RHO(2),RHOIN(2)
      REAL(q) EXC,EXCD(2),EXCD_(2),EXCDD(2,2),EXCDD_(2,2),EXCDDD(2,2,2)
      REAL(q) TMP(-1:1,-1:1),T2(2,-1:1,-1:1)
      REAL(q) EPS(2)
      ! this is the best compromise for densities between 1000 and 0.1
      REAL(q),PARAMETER :: DELTA=1E-3_q
      INTEGER I,J

      EXCDD=0
      EXCDDD=0
      RHO(1)=RHOUP
      RHO(2)=RHODOWN
    ! function + derivative
      CALL EXCOR(RHO,EXC,EXCD,TREL)
      IF (NDER>0) THEN
        TMP(0,0) =EXC
        T2(:,0,0)=EXCD
        RHOIN=RHO
        ! calculate steps
        ! exc is approx rho^(4/3)
        ! the derivative thus rho^(1/3) and thus the derivative is
        ! of the order 1/rho
        EPS(1)=DELTA*MAX(RHO(1),RHO(2))
        EPS(2)=EPS(1)
        ! calculate all values on the 3x3 rectangle
        DO I=-1,1
          DO J=-1,1
             RHO(1)=RHOIN(1)+EPS(1)*I
             RHO(2)=RHOIN(2)+EPS(2)*J
             CALL EXCOR(RHO,TMP(I,J),T2(1,I,J),TREL)
          ENDDO
        ENDDO

    ! 1st and 2nd derivative of exchange correlation energy
    ! 1st derivative (EXCD_) is of course equal to EXCD
        EXCD_(1)   =(TMP(1,0)-TMP(-1,0)) / (2*EPS(1))
        EXCD_(2)   =(TMP(0,1)-TMP(0,-1)) / (2*EPS(2))
        EXCDD_(1,1)=(TMP(1,0)+TMP(-1,0)-2*TMP(0,0)) / (EPS(1)*EPS(1))
        EXCDD_(2,2)=(TMP(0,1)+TMP(0,-1)-2*TMP(0,0)) / (EPS(2)*EPS(2))
        EXCDD_(1,2)=((TMP(1,1)+TMP(-1,-1)-2*TMP(0,0))- &
                     (TMP(1,0)+TMP(-1,0)-2*TMP(0,0)) - &
                     (TMP(0,1)+TMP(0,-1)-2*TMP(0,0)))/(2*EPS(1)*EPS(2))
        EXCDD_(2,1)=EXCDD_(1,2)
    ! 1st derivative of potential = 2nd derivative of energy
        EXCDD (1,1)=(T2(1,1,0)-T2(1,-1,0)) / (2*EPS(1))
        EXCDD (2,2)=(T2(2,0,1)-T2(2,0,-1)) / (2*EPS(1))
        EXCDD (1,2)=(T2(1,0,1)-T2(1,0,-1)) / (2*EPS(1))
        EXCDD (2,1)=(T2(2,1,0)-T2(2,-1,0)) / (2*EPS(1))
    ! 2nd derivative of potential = 3rd derivative of energy
        EXCDDD(1,1,1)=(T2(1,1,0)+T2(1,-1,0)-2*T2(1,0,0)) / (EPS(1)*EPS(1))
        EXCDDD(1,2,2)=(T2(1,0,1)+T2(1,0,-1)-2*T2(1,0,0)) / (EPS(2)*EPS(2))
        EXCDDD(1,1,2)=((T2(1,1,1)+T2(1,-1,-1)-2*T2(1,0,0))- &
                       (T2(1,1,0)+T2(1,-1,0)-2*T2(1,0,0))- &
                       (T2(1,0,1)+T2(1,0,-1)-2*T2(1,0,0)))/(2*EPS(1)*EPS(2))
        EXCDDD(1,2,1)=EXCDDD(1,1,2)
        EXCDDD(2,1,1)=(T2(2,1,0)+T2(2,-1,0)-2*T2(2,0,0)) / (EPS(1)*EPS(1))
        EXCDDD(2,2,2)=(T2(2,0,1)+T2(2,0,-1)-2*T2(2,0,0)) / (EPS(2)*EPS(2))
        EXCDDD(2,1,2)=((T2(2,1,1)+T2(2,-1,-1)-2*T2(2,0,0)) - &
                       (T2(2,1,0)+T2(2,-1,0)-2*T2(2,0,0)) - &
                       (T2(2,0,1)+T2(2,0,-1)-2*T2(2,0,0)))/(2*EPS(1)*EPS(2))
        EXCDDD(2,2,1)=EXCDDD(2,1,2)
    !    WRITE(*,*) EXCDDD(2,1,2),EXCDDD(1,2,2),EXCDDD(2,2,2)
    !    WRITE(*,*) EXCDDD(1,1,2),EXCDDD(2,1,1),EXCDDD(1,1,1)
    ! plenty of cross checks are possible here
    !    WRITE(*,'(4E14.7)') EXCD,EXCD_
    !    WRITE(*,'(4E14.7)') EXCDD,EXCDD_
    !    WRITE(*,'(4E14.7)') EXCDDD
      ENDIF

    END SUBROUTINE EXCOR_DER


    SUBROUTINE EXCOR(RHO,EXC,EXCD,TREL)
      USE prec
      USE constant
      USE xclib

      IMPLICIT NONE
      LOGICAL TREL
      REAL(q) RHO(2),EXC,EXCD(2)
    ! local variables
      REAL(q) RS,ZETA,FZ,DFZ,EXP,EC,ECP,VXP,VCP,EXT,VX0,DVX,VX1,VX2,RH, &
              ECF,VXF,VCF,EXF,ALP,VALP,ZETA3,ZETA4

      RH=RHO(1)+RHO(2)
      RS = ( (3._q/(4*PI)) / RH)**(1/3._q) /AUTOA
      ZETA=(RHO(1)-RHO(2))/(RHO(1)+RHO(2))
      ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)
      ZETA3=(ZETA*ZETA)*ZETA
      ZETA4=(ZETA*ZETA)*(ZETA*ZETA)
#ifndef USELIBXC
      IF (LEXCH==-1) THEN
#else
      IF (LEXCH==-1.OR.LEXCH==99) THEN
#endif
         EXC=0
         EXCD=0
      !
      ! LDA part: Pade approximation to Ceperley Alder results
      !
!aem AM05 functional (LEXCH==13) added, it uses the PW correlation in the LDA part.
! jP: adding PBEsol
!vdw jk
      ELSE IF (LEXCH==8 .OR. LEXCH==9 .OR. LEXCH==13 .OR. LEXCH==14 .or. LEXCH==17 .or. (LEXCH.ge.40 .and. LEXCH.le.50)) THEN
!vdw jk
         FZ =FZ0(ZETA)          ! interpolation function for exchange from pm to fm
         DFZ=FZ1(ZETA)*SIGN(1._q,ZETA)
         EXP=EX_MOD(RS,1,.FALSE.)*RH*RYTOEV ; EXF=EX_MOD(RS,2,.FALSE.)*RH*RYTOEV
         VXP=VX_MOD(RS,1,.FALSE.)*RYTOEV    ; VXF=VX_MOD(RS,2,.FALSE.)*RYTOEV

         CALL CORPBE_LDA(RS,ZETA,EC,VCP,VCF)
         EC=EC*ALDAC*RYTOEV*RH ; VCP=VCP*ALDAC*RYTOEV ; VCF=VCF*ALDAC*RYTOEV
         
         VX0=VXP +(VXF-VXP)*FZ
         DVX=(EXF-EXP)*DFZ/RH
         EXC= EXP+(EXF-EXP)*FZ+EC

         EXCD(1)= VX0-DVX*(ZETA-1)+VCP
         EXCD(2)= VX0-DVX*(ZETA+1)+VCF
      !
      ! LDA: short range exchange + short range correlation
      !   
! jH-new noch fuer RPA und RPA-plus 
      ELSEIF (LEXCH==20) THEN
         FZ =FZ0(ZETA)          ! interpolation function for exchange from pm to fm
         DFZ=FZ1(ZETA)*SIGN(1._q,ZETA)
         EXP=EX_MOD(RS,1,.FALSE.)*RH*RYTOEV ; EXF=EX_MOD(RS,2,.FALSE.)*RH*RYTOEV
         VXP=VX_MOD(RS,1,.FALSE.)*RYTOEV    ; VXF=VX_MOD(RS,2,.FALSE.)*RYTOEV

         CALL CORPBE_LDA_RPA(RS,ZETA,EC,VCP,VCF)
         EC=EC*ALDAC*RYTOEV*RH ; VCP=VCP*ALDAC*RYTOEV ; VCF=VCF*ALDAC*RYTOEV
         
         VX0=VXP +(VXF-VXP)*FZ
         DVX=(EXF-EXP)*DFZ/RH
         EXC= EXP+(EXF-EXP)*FZ+EC

         EXCD(1)= VX0-DVX*(ZETA-1)+VCP
         EXCD(2)= VX0-DVX*(ZETA+1)+VCF
      ELSEIF (LEXCH==30) THEN
         FZ =FZ0(ZETA)          ! interpolation function for exchange from pm to fm
         DFZ=FZ1(ZETA)*SIGN(1._q,ZETA)
         EXP=EX_MOD(RS,1,.FALSE.)*RH*RYTOEV ; EXF=EX_MOD(RS,2,.FALSE.)*RH*RYTOEV
         VXP=VX_MOD(RS,1,.FALSE.)*RYTOEV    ; VXF=VX_MOD(RS,2,.FALSE.)*RYTOEV

         CALL CORPBE_LDA_RPA_PLUS(RS,ZETA,EC,VCP,VCF)
         EC=EC*ALDAC*RYTOEV*RH ; VCP=VCP*ALDAC*RYTOEV ; VCF=VCF*ALDAC*RYTOEV

         VX0=VXP +(VXF-VXP)*FZ
         DVX=(EXF-EXP)*DFZ/RH
         EXC= EXP+(EXF-EXP)*FZ+EC

         EXCD(1)= VX0-DVX*(ZETA-1)+VCP
         EXCD(2)= VX0-DVX*(ZETA+1)+VCF


      ELSEIF (LEXCH==10) THEN
! jP used interpolation: PBE (LEXCH==8)
         FZ =FZ0(ZETA)          ! interpolation function for exchange from pm to fm
         DFZ=FZ1(ZETA)*SIGN(1._q,ZETA)
         EXP=EX_MOD(RS,1,.FALSE.)*RH*RYTOEV ; EXF=EX_MOD(RS,2,.FALSE.)*RH*RYTOEV
         VXP=VX_MOD(RS,1,.FALSE.)*RYTOEV    ; VXF=VX_MOD(RS,2,.FALSE.)*RYTOEV

         EC= EC_MOD(RS,ZETA,2,TREL)*ALDAC*RYTOEV*RH 
         VCP=VC_MOD(RS,0._q,1,TREL)*ALDAC*RYTOEV 
         VCF=VC_MOD(RS,1._q,2,TREL)*ALDAC*RYTOEV
         
         VX0=VXP +(VXF-VXP)*FZ
         DVX=(EXF-EXP)*DFZ/RH
         EXC= EXP+(EXF-EXP)*FZ+EC

         EXCD(1)= VX0-DVX*(ZETA-1)+VCP
         EXCD(2)= VX0-DVX*(ZETA+1)+VCF

      ELSE IF (LEXCH==100) THEN
         EXC=0
         CALL COHSM1(RH, EXCD(1))
         EXCD(2)=EXCD(1)
      !
      ! LDA part: Perdew Zungers interpolation
      !
      ELSEIF (LEXCH==11) THEN
! Joachim Paier
         FZ =FZ0(ZETA)          ! interpolation function for exchange from pm to fm
         DFZ=FZ1(ZETA)*SIGN(1._q,ZETA)
         EXP=EX_MOD(RS,1,TREL)*RH*RYTOEV ; EXF=EX_MOD(RS,2,TREL)*RH*RYTOEV
         ECP=ECVOIII(RS,1)*RH*RYTOEV     ; ECF=ECVOIII(RS,2)*RH*RYTOEV
!         ECP=ECVO(RS,1)*RH*RYTOEV     ; ECF=ECVO(RS,2)*RH*RYTOEV
         VXP=VX_MOD(RS,1,TREL)*RYTOEV    ; VXF=VX_MOD(RS,2,TREL)*RYTOEV
         VCP=VCVOIII(RS,1)*RYTOEV        ; VCF=VCVOIII(RS,2)*RYTOEV
!         VCP=VCVO(RS,1)*RYTOEV        ; VCF=VCVO(RS,2)*RYTOEV

         ALP =ALPHA0_III(RS)*RH*RYTOEV
!         ALP =ALPHA0(RS)*RH*RYTOEV
         VALP=ALPHA1_III(RS)*RYTOEV
!         VALP=ALPHA1(RS)*RYTOEV

         VX0=VXP+(VXF-VXP)*FZ+  (VCP+(VALP+(VCF-VCP-VALP)*ZETA4)*FZ)*ALDAC
         DVX=((EXF-EXP)*DFZ  +(ALP+(ECF-ECP-ALP)*ZETA4)*DFZ*ALDAC+ &
              4*(ECF-ECP-ALP)*ZETA3 *FZ *ALDAC )/RH
        ! the more usual expression for this is
        ! ECP*( 1 - FZ Z4) +EP*FZ*Z4-ALP*F*(1._q-Z4)/(8/(9 gamma))
         EXC= EXP+(EXF-EXP)*FZ +ECP*ALDAC+(ALP+(ECF-ECP-ALP)*ZETA4)*FZ*ALDAC

         EXCD(1)= VX0-DVX*(ZETA-1)
         EXCD(2)= VX0-DVX*(ZETA+1)
      ELSEIF (LEXCH==12) THEN
! Joachim Paier
         FZ =FZ0(ZETA)          ! interpolation function for exchange from pm to fm
         DFZ=FZ1(ZETA)*SIGN(1._q,ZETA)
         EXP=EX_MOD(RS,1,TREL)*RH*RYTOEV ; EXF=EX_MOD(RS,2,TREL)*RH*RYTOEV
         ECP=ECVO(RS,1)*RH*RYTOEV        ; ECF=ECVO(RS,2)*RH*RYTOEV
         VXP=VX_MOD(RS,1,TREL)*RYTOEV    ; VXF=VX_MOD(RS,2,TREL)*RYTOEV
         VCP=VCVO(RS,1)*RYTOEV           ; VCF=VCVO(RS,2)*RYTOEV

         ALP =ALPHA0(RS)*RH*RYTOEV
         VALP=ALPHA1(RS)*RYTOEV

         VX0=VXP+(VXF-VXP)*FZ+  (VCP+(VALP+(VCF-VCP-VALP)*ZETA4)*FZ)*ALDAC
         DVX=((EXF-EXP)*DFZ  +(ALP+(ECF-ECP-ALP)*ZETA4)*DFZ*ALDAC+ &
              4*(ECF-ECP-ALP)*ZETA3 *FZ *ALDAC )/RH
        ! the more usual expression for this is
        ! ECP*( 1 - FZ Z4) +EP*FZ*Z4-ALP*F*(1._q-Z4)/(8/(9 gamma))
         EXC= EXP+(EXF-EXP)*FZ +ECP*ALDAC+(ALP+(ECF-ECP-ALP)*ZETA4)*FZ*ALDAC

         EXCD(1)= VX0-DVX*(ZETA-1)
         EXCD(2)= VX0-DVX*(ZETA+1)
      ELSE

         FZ =FZ0(ZETA)          ! interpolation function for exchange from pm to fm
         DFZ=FZ1(ZETA)*SIGN(1._q,ZETA)
         EXP=EX_MOD(RS,1,TREL)*RH*RYTOEV ; EXF=EX_MOD(RS,2,TREL)*RH*RYTOEV
         ECP=ECCA(RS,1)*RH*RYTOEV        ; ECF=ECCA(RS,2)*RH*RYTOEV
         VXP=VX_MOD(RS,1,TREL)*RYTOEV    ; VXF=VX_MOD(RS,2,TREL)*RYTOEV
         VCP=VCCA(RS,1)*RYTOEV           ; VCF=VCCA(RS,2)*RYTOEV

         IF (LFCI==1) THEN
            ALP =ALPHA0(RS)*RH*RYTOEV
            VALP=ALPHA1(RS)*RYTOEV
            
            VX0=VXP+(VXF-VXP)*FZ+  (VCP+(VALP+(VCF-VCP-VALP)*ZETA4)*FZ)*ALDAC
            DVX=((EXF-EXP)*DFZ  +(ALP+(ECF-ECP-ALP)*ZETA4)*DFZ*ALDAC+ &
                              4*(ECF-ECP-ALP)*ZETA3 *FZ *ALDAC )/RH
        ! the more usual expression for this is
        ! ECP*( 1 - FZ Z4) +EP*FZ*Z4-ALP*F*(1._q-Z4)/(8/(9 gamma))
            EXC= EXP+(EXF-EXP)*FZ +ECP*ALDAC+(ALP+(ECF-ECP-ALP)*ZETA4)*FZ*ALDAC
         ELSE
            VX0= VXP+(VXF-VXP)*FZ +VCP*ALDAC+(VCF-VCP)*FZ*ALDAC
            DVX=((EXF-EXP)*DFZ    +(ECF-ECP)*DFZ)/RH*ALDAC
            EXC= EXP+(EXF-EXP)*FZ +ECP*ALDAC+(ECF-ECP)*FZ*ALDAC
         ENDIF

         EXCD(1)= VX0-DVX*(ZETA-1)
         EXCD(2)= VX0-DVX*(ZETA+1)
      ENDIF

    END SUBROUTINE EXCOR

  END MODULE setexm
