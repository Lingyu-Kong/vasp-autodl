#include  "symbol.inc"

!************************ SUBROUTINE GGAALL *****************************
!
!  switch between different GGAs
!  presently only PW91, PBE and RPBE are implemented
!  (i.e. d exc / d rho and  d exc / d | grad rho | are calculated
!  directly
!  for other GGA functional finite differences are used to calculate
!  the required derivatives
!  LLDA allows to include the LDA contribution directly in this
!  routine.
!  This only works for PBE and RPBE and fails in all other cases
!
!************************************************************************

      SUBROUTINE GGAALL(D,DD, EXC,EXCD,EXCDD,LLDA)
!$ACC ROUTINE SEQ
      USE prec
      USE constant
      USE setexm_struct_def
#ifdef USELIBXC
      USE xc_f03_lib_m
#endif

      IMPLICIT REAL(q) (A-H,O-Z)
#ifdef USELIBXC
      REAL(q) DTMP(1),SIGMATMP(1),EXC1TMP(1),EXCD1TMP(1),EXCDD1TMP(1)
      REAL(q) SIGMAZERO(1),EXC1ZERO(1),EXCD1ZERO(1),EXCDD1ZERO(1)
#endif
      PARAMETER (DDELTA=1E-4_q)
      PARAMETER (THRD=1._q/3._q)
      LOGICAL LLDA
      LOGICAL:: FORCE_PBE=.FALSE.
!$ACC ROUTINE (GGA91_WB,EXCHPBE,CALC_EXCHWPBE_SP,CORunspPBE,B3LYPXCS,AM05, &
!$ACC& EXCHPBESOL,CALC_EXCHWPBEsol_SP,CORunspPBESOL,GGAEALL) SEQ
#ifdef USELIBXC
!$ACC ROUTINE (XC_F03_LDA_VXC,XC_F03_LDA_EXC_VXC,XC_F03_GGA_VXC,XC_F03_GGA_EXC_VXC) SEQ
#endif

      IF (LEXCH==7) THEN

        ! PW91 using the routines of Bird and White

        CALL GGA91_WB(D,DD,EXC,EXCD,EXCDD)
        EXC = 2*EXC / D
        EXCD =2*EXCD
        EXCDD=2*EXCDD
!        WRITE(*,'(10F14.7)') D, EXC,EXCD,EXCDD
!vdw jk
      ELSE IF (LEXCH==8 .OR. LEXCH==9 .OR. (LEXCH.GE.40 .AND. LEXCH.LE.50)) THEN
!vdw jk

        ! Perdew Burke Ernzerhof and revised functional

        IF (LEXCH==8) THEN
           ukfactor=1.0_q
        ELSE
           ukfactor=0.0_q
        ENDIF

        IF (D<=0) THEN
           EXC   = 0._q
           EXCD  = 0._q
           EXCDD = 0._q
           RETURN
        ENDIF

        DTHRD=exp(log(D)*THRD)
        RS=(0.75_q/PI)**THRD/DTHRD
        FK=(3._q*PI*PI)**THRD*DTHRD
        SK = SQRT(4.0_q*FK/PI)
        IF(D>1.E-10_q)THEN
           S=DD/(D*FK*2._q)
           T=DD/(D*SK*2._q)
        ELSE
           S=0.0_q
           T=0.0_q
        ENDIF
! Iann Gerber: range separating in GGA exchange
        IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
           CALL EXCHPBE(D,DTHRD,S,EXLDA,EXC,EXDLDA,EXCD,EXCDD, &
          &     ukfactor)
        ELSE
           CALL CALC_EXCHWPBE_SP(D,S, &
          &  EXLDA,EXDLDA,EXLDA_SR,EXLDA_LR,EXDLDA_SR,EXDLDA_LR, &
          &  EXWPBE_SR,EXWPBE_LR,EXWPBED_SR,EXWPBED_LR,EXWPBEDD_SR,EXWPBEDD_LR)
        ENDIF
! Iann Gerber: end modification

        CALL CORunspPBE(RS,ECLDA,ECDLDA,SK, &
             T,EC,ECD,ECDD,.TRUE.)

!        WRITE(*,'(10F14.7)') D,S,(EXC-EXLDA)/D,EXCD-EXDLDA,EXCDD
!        WRITE(*,'(10F14.7)') D,RS,SK,T,EC,ECD,ECDD

        IF (LLDA) THEN
        ! Add LDA contributions as well
           IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI  .OR. FORCE_PBE) THEN
              EXC  =EC  *AGGAC +(EXC-EXLDA)/D*AGGAX+EXLDA/D*LDAX+ECLDA  ! eventuell noch mit ALDAC multi
              EXCD =ECD *AGGAC +(EXCD-EXDLDA)*AGGAX+EXDLDA*LDAX +ECDLDA ! eventuell noch mit ALDAC multi
              EXCDD=ECDD*AGGAC + EXCDD*AGGAX
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
              ! Iann Gerber's functionals
                 ! gk replaced ECLD by ECLDA, 18.09.2018
                 EXC  =EC  *AGGAC+ EXWPBE_LR/D*AGGAX +EXWPBE_SR/D +ECLDA! ist EXWPBE nur diff
                 EXCD =ECD *AGGAC+ EXWPBED_LR*AGGAX+EXWPBED_SR+ECDLDA   !
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_LR*AGGAX +EXWPBEDD_SR !
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC  *AGGAC +EXWPBE_LR/D*AGGAX +ECLDA
                 EXCD =ECD *AGGAC +EXWPBED_LR*AGGAX+ECDLDA
                 EXCDD=ECDD*AGGAC +EXWPBEDD_LR*AGGAX
              ELSE
              ! wPBE
                 EXC  =EC  *AGGAC+ EXWPBE_SR/D*AGGAX +EXWPBE_LR/D +ECLDA
                 EXCD =ECD *AGGAC+ EXWPBED_SR*AGGAX+EXWPBED_LR+ECDLDA
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_SR*AGGAX +EXWPBEDD_LR
              ENDIF
           ENDIF
        ELSE
        ! Do not add LDA contributions
           IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
              EXC  =EC  *AGGAC +(EXC-EXLDA)/D*AGGAX
              EXCD =ECD *AGGAC +(EXCD-EXDLDA)*AGGAX
              EXCDD=ECDD*AGGAC + EXCDD*AGGAX
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
              ! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+(EXWPBE_LR-EXLDA_LR)/D*AGGAX  +(EXWPBE_SR-EXLDA_SR)/D
                 EXCD =ECD *AGGAC+(EXWPBED_LR-EXDLDA_LR)*AGGAX+(EXWPBED_SR-EXDLDA_SR)
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_LR*AGGAX +EXWPBEDD_SR
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC  *AGGAC+(EXWPBE_LR-EXLDA_LR)/D*AGGAX
                 EXCD =ECD *AGGAC+(EXWPBED_LR-EXDLDA_LR)*AGGAX
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_LR*AGGAX
              ELSE
              ! wPBE
                 EXC  =EC  *AGGAC+(EXWPBE_SR-EXLDA_SR)/D*AGGAX  +(EXWPBE_LR-EXLDA_LR)/D
                 EXCD =ECD *AGGAC+(EXWPBED_SR-EXDLDA_SR)*AGGAX+(EXWPBED_LR-EXDLDA_LR)
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_SR*AGGAX +EXWPBEDD_LR
              ENDIF
           ENDIF
        ENDIF

        ! Hartree -> Rydberg conversion
        EXC = 2*EXC
        EXCD =2*EXCD
        EXCDD=2*EXCDD

!        WRITE(*,'(10F14.7)') D, EXC,EXCD,EXCDD
! jP
      ELSE IF (LEXCH==11 .OR. LEXCH==12) THEN

        RHO1=D/2.0_q
        RHO2=D/2.0_q
        DRHO1=DD/2.0_q
        DRHO2=DD/2.0_q

        CALL  B3LYPXCS(RHO1,RHO2,DRHO1,DRHO2,DD,EXC,EXCD1,EXCDD1,EXCD2,EXCDD2,EXCDDA)


        EXC     = 2._q* EXC/D
        EXCD1   = 2._q* EXCD1
        EXCD2   = 2._q* EXCD2
        EXCDD1  = 2._q* EXCDD1
        EXCDD2  = 2._q* EXCDD2
        EXCDDA  = 2._q* EXCDDA

        EXCD    = EXCD1
        EXCDD   = EXCDD1+EXCDDA

!       WRITE(91,'(7F14.7)') RHO1, RHO2, DRHO1, DRHO2, EXC, EXCD, EXCDD

! jP
!aem AM05 added

      ELSE IF (LEXCH==13) THEN

        !aem DD sometimes comes in negative.
        !aem Since AM05 assumes that DD actually IS |grad rho|,
        !aem and thus DD always should be positive,
        !aem we need to use ABS(DD) as input to AM05.
        !aem However, my routine needs to give the same symmetries out as
        !aem PBE and PW91, therefor the sign-compensation of EXCDD.

        D1=D/2._q
        D2=D/2._q
        DD1=DD/2._q
        DD2=DD/2._q
        ider = 1
        CALL AM05(D1,D2,ABS(DD1),ABS(DD2),ider, &
             FXC,FXCLDA,VXCLDA1,VXCLDA2,FXCD1,FXCD2, &
             FXCDD1,FXCDD2)

        IF (LLDA) THEN
          EXC = FXC/D
          EXCD = FXCD1
        ELSE
          EXC = (FXC - FXCLDA)/D
          EXCD = FXCD1 - VXCLDA1
        ENDIF
        EXCDD = SIGN(1.0_q,DD)*FXCDD1

        ! Hartree -> Rydberg conversion
        EXC = 2*EXC
        EXCD =2*EXCD
        EXCDD=2*EXCDD



      ELSE IF (LEXCH==14) THEN

        ! Perdew Burke Ernzerhof and revised functional

           ukfactor=1.0_q

        IF (D<=0) THEN
           EXC   = 0._q
           EXCD  = 0._q
           EXCDD = 0._q
           RETURN
        ENDIF

        DTHRD=exp(log(D)*THRD)
        RS=(0.75_q/PI)**THRD/DTHRD
        FK=(3._q*PI*PI)**THRD*DTHRD
        SK = SQRT(4.0_q*FK/PI)
        IF(D>1.E-10_q)THEN
           S=DD/(D*FK*2._q)
           T=DD/(D*SK*2._q)
        ELSE
           S=0.0_q
           T=0.0_q
        ENDIF
! Iann Gerber: range separating in GGA exchange
        IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
           CALL EXCHPBESOL(D,DTHRD,S,EXLDA,EXC,EXDLDA,EXCD,EXCDD, &
          &     ukfactor)
        ELSE
           CALL CALC_EXCHWPBEsol_SP(D,S, &
          &  EXLDA,EXDLDA,EXLDA_SR,EXLDA_LR,EXDLDA_SR,EXDLDA_LR, &
          &  EXWPBE_SR,EXWPBE_LR,EXWPBED_SR,EXWPBED_LR,EXWPBEDD_SR,EXWPBEDD_LR)
        ENDIF
! Iann Gerber: end modification

        CALL CORunspPBESOL(RS,ECLDA,ECDLDA,SK, &
             T,EC,ECD,ECDD,.TRUE.)

!        WRITE(*,'(10F14.7)') D,S,(EXC-EXLDA)/D,EXCD-EXDLDA,EXCDD
!        WRITE(*,'(10F14.7)') D,RS,SK,T,EC,ECD,ECDD

        IF (LLDA) THEN
        ! Add LDA contributions as well
           IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI  .OR. FORCE_PBE) THEN
              EXC  =EC  *AGGAC +(EXC-EXLDA)/D*AGGAX+EXLDA/D*LDAX+ECLDA
              EXCD =ECD *AGGAC +(EXCD-EXDLDA)*AGGAX+EXDLDA*LDAX +ECDLDA
              EXCDD=ECDD*AGGAC + EXCDD*AGGAX
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
              ! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+ EXWPBE_LR/D*AGGAX +EXWPBE_SR/D +ECLDA
                 EXCD =ECD *AGGAC+ EXWPBED_LR*AGGAX+EXWPBED_SR+ECDLDA
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_LR*AGGAX +EXWPBEDD_SR
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC  *AGGAC +EXWPBE_LR/D*AGGAX +ECLDA
                 EXCD =ECD *AGGAC +EXWPBED_LR*AGGAX+ECDLDA
                 EXCDD=ECDD*AGGAC +EXWPBEDD_LR*AGGAX
              ELSE
              ! wPBE
                 EXC  =EC  *AGGAC+ EXWPBE_SR/D*AGGAX +EXWPBE_LR/D +ECLDA
                 EXCD =ECD *AGGAC+ EXWPBED_SR*AGGAX+EXWPBED_LR+ECDLDA
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_SR*AGGAX +EXWPBEDD_LR
              ENDIF
           ENDIF
        ELSE
        ! Do not add LDA contributions
           IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
              EXC  =EC  *AGGAC +(EXC-EXLDA)/D*AGGAX
              EXCD =ECD *AGGAC +(EXCD-EXDLDA)*AGGAX
              EXCDD=ECDD*AGGAC + EXCDD*AGGAX
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
              ! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+(EXWPBE_LR-EXLDA_LR)/D*AGGAX  +(EXWPBE_SR-EXLDA_SR)/D
                 EXCD =ECD *AGGAC+(EXWPBED_LR-EXDLDA_LR)*AGGAX+(EXWPBED_SR-EXDLDA_SR)
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_LR*AGGAX +EXWPBEDD_SR
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC  *AGGAC+(EXWPBE_LR-EXLDA_LR)/D*AGGAX
                 EXCD =ECD *AGGAC+(EXWPBED_LR-EXDLDA_LR)*AGGAX
                 EXCDD=ECDD*AGGAC+EXWPBEDD_LR*AGGAX
              ELSE
              ! wPBE
                 EXC  =EC  *AGGAC+(EXWPBE_SR-EXLDA_SR)/D*AGGAX  +(EXWPBE_LR-EXLDA_LR)/D
                 EXCD =ECD *AGGAC+(EXWPBED_SR-EXDLDA_SR)*AGGAX+(EXWPBED_LR-EXDLDA_LR)
                 EXCDD=ECDD*AGGAC+ EXWPBEDD_SR*AGGAX +EXWPBEDD_LR
              ENDIF
           ENDIF
        ENDIF

        ! Hartree -> Rydberg conversion
        EXC = 2*EXC
        EXCD =2*EXCD
        EXCDD=2*EXCDD

!        WRITE(*,'(10F14.7)') D, EXC,EXCD,EXCDD
! jP: end of adding PBEsol

! BEEF
#ifdef libbeef
      ELSE IF (LEXCH==17) THEN
        IF (D>1E-10_q) THEN
           DDSQ = DD*DD
           CALL BEEFX(D, DDSQ, EXBEEF, EXBEEFD, EXBEEFDD, LLDA)
           CALL BEEFLOCALCORR(D, DDSQ, ECBEEF, ECBEEFD, ECBEEFDD, LLDA)
           EXC = 2._q*(EXBEEF + ECBEEF)/D
           EXCD = 2._q*(EXBEEFD + ECBEEFD)
           EXCDD = 2._q*(EXBEEFDD + ECBEEFDD)*DD
        ELSE
           EXC   = 0._q
           EXCD  = 0._q
           EXCDD = 0._q
        ENDIF
#endif


      ELSE IF (LEXCH>20.AND.LEXCH<25) THEN
! functionals for range-separated ACFDT (LDA - short range RPA)
! dummy stub since these functionals are not gradient corrected
         EXC=0._q ; EXCD=0._q ; EXCDD=0._q


      ELSE IF (LEXCH==20.OR.LEXCH==30) THEN
         EXC=0._q ; EXCD=0._q ; EXCDD=0._q

#ifdef USELIBXC
      ELSE IF (LEXCH==99) THEN
         DTMP(1)=D
         SIGMATMP(1)=DD*DD
         SIGMAZERO(1)=0._q
         EXC   = 0._q
         EXCD  = 0._q
         EXCDD = 0._q
         CX = -(3._q/4._q)*(3._q/PI)**THRD
         EXLDA = CX*D**THRD
         EXDLDA = EXLDA*4._q/3._q

         DO IFUNC=1,2
            IF (LIBXC_MASK(IFUNC)>0) THEN
               IF ((XC_FAMILY(IFUNC) .EQ. LIBXC_FAMILY_LDA) .OR. (XC_FAMILY(IFUNC) .EQ. LIBXC_FAMILY_HYB_LDA)) THEN
                  IF (LIBXC_NOEXC(IFUNC)) THEN
                     CALL XC_F03_LDA_VXC(LIBXC_FUNC_NSP(IFUNC),1_8,DTMP,EXCD1TMP)
                     EXC1TMP = 0._q
                  ELSE
                     CALL XC_F03_LDA_EXC_VXC(LIBXC_FUNC_NSP(IFUNC),1_8,DTMP,EXC1TMP,EXCD1TMP)
                  ENDIF
                  IF ((XC_KIND(IFUNC) .EQ. XC_EXCHANGE) .AND. (LIBXC_HYB(IFUNC) .EQV. .FALSE.)) THEN
                     EXC   = EXC   + LDAX*EXC1TMP(1)
                     EXCD  = EXCD  + LDAX*EXCD1TMP(1)
                  ELSEIF (XC_KIND(IFUNC) .EQ. XC_CORRELATION) THEN
                     EXC   = EXC   + ALDAC*EXC1TMP(1)
                     EXCD  = EXCD  + ALDAC*EXCD1TMP(1)
                  ELSEIF ((XC_KIND(IFUNC) .EQ. XC_EXCHANGE_CORRELATION) .OR. (LIBXC_HYB(IFUNC) .EQV. .TRUE.)) THEN
                     EXC   = EXC   + EXC1TMP(1)
                     EXCD  = EXCD  + EXCD1TMP(1)
                  ENDIF
               ELSEIF ((XC_FAMILY(IFUNC) .EQ. LIBXC_FAMILY_GGA) .OR. (XC_FAMILY(IFUNC) .EQ. LIBXC_FAMILY_HYB_GGA)) THEN
                  IF (LIBXC_NOEXC(IFUNC)) THEN
                     CALL XC_F03_GGA_VXC(LIBXC_FUNC_NSP(IFUNC),1_8,DTMP,SIGMATMP,EXCD1TMP,EXCDD1TMP)
                     EXC1TMP = 0._q
                  ELSE
                     CALL XC_F03_GGA_EXC_VXC(LIBXC_FUNC_NSP(IFUNC),1_8,DTMP,SIGMATMP,EXC1TMP,EXCD1TMP,EXCDD1TMP)
                  ENDIF
                  IF ((XC_KIND(IFUNC) .EQ. XC_EXCHANGE) .AND. (LIBXC_HYB(IFUNC) .EQV. .FALSE.)) THEN
                     EXC   = EXC   + AGGAX*(EXC1TMP(1)-EXLDA) + LDAX*EXLDA
                     EXCD  = EXCD  + AGGAX*(EXCD1TMP(1)-EXDLDA) + LDAX*EXDLDA
                     EXCDD = EXCDD + AGGAX*EXCDD1TMP(1)*2._q*DD
                  ELSEIF (XC_KIND(IFUNC) .EQ. XC_CORRELATION) THEN
                     IF (LIBXC_NOEXC(IFUNC)) THEN
                        CALL XC_F03_GGA_VXC(LIBXC_FUNC_NSP(IFUNC),1_8,DTMP,SIGMAZERO,EXCD1ZERO,EXCDD1ZERO)
                        EXC1ZERO = 0._q
                     ELSE
                        CALL XC_F03_GGA_EXC_VXC(LIBXC_FUNC_NSP(IFUNC),1_8,DTMP,SIGMAZERO,EXC1ZERO,EXCD1ZERO,EXCDD1ZERO)
                     ENDIF
                     EXC   = EXC   + AGGAC*(EXC1TMP(1)-EXC1ZERO(1)) + ALDAC*EXC1ZERO(1)
                     EXCD  = EXCD  + AGGAC*(EXCD1TMP(1)-EXCD1ZERO(1)) + ALDAC*EXCD1ZERO(1)
                     EXCDD = EXCDD + AGGAC*EXCDD1TMP(1)*2._q*DD
                  ELSEIF ((XC_KIND(IFUNC) .EQ. XC_EXCHANGE_CORRELATION) .OR. (LIBXC_HYB(IFUNC) .EQV. .TRUE.)) THEN
                     EXC   = EXC   + EXC1TMP(1)
                     EXCD  = EXCD  + EXCD1TMP(1)
                     EXCDD = EXCDD + EXCDD1TMP(1)*2._q*DD
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         !Conversion to Rydbergs
         EXC   = 2._q*EXC
         EXCD  = 2._q*EXCD
         EXCDD = 2._q*EXCDD
#endif
      ELSE

        ! for all other functionals
        ! we use finite differences to calculate the required
        ! quantities
        ! presently no other functional are supported
        ! but in case one needs these routines

        DELTA=MIN(DDELTA,ABS(D)/100)
        D1=D-DELTA
        D2=D+DELTA
        CALL GGAEALL(D1,DD,1._q,1._q,VXC,EXC1)
        CALL GGAEALL(D2,DD,1._q,1._q,VXC,EXC2)
        EXCD=(EXC2*D2-EXC1*D1)/MAX((D2-D1),1E-10_q)

        DELTA=MIN(DDELTA,ABS(DD)/100)
        DD1=DD-DELTA
        DD2=DD+DELTA
        CALL GGAEALL(D,DD1,1._q,1._q,VXC,EXC1)
        CALL GGAEALL(D,DD2,1._q,1._q,VXC,EXC2)
        EXCDD=(EXC2*D-EXC1*D)/MAX((DD2-DD1),1E-10_q)
        CALL GGAEALL(D,DD,1._q,1._q,VXC,EXC)

!aem added in case somebody (like me) might want to use this routine
        ! Hartree -> Rydberg conversion
        EXC = 2*EXC
        EXCD =2*EXCD
        EXCDD=2*EXCDD
!aem end addition

      ENDIF
      RETURN
      END

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE GGAALL_SIMD(D_IN, DD, EXC, EXCD, EXCDD, LLDA, MASK0)
      USE prec
      USE constant
      USE setexm_struct_def
      USE simd
      USE tutor, ONLY: vtutor
      IMPLICIT REAL(q) (A-H, O-Z)
      TYPE(SIMD_REAL64), INTENT(IN)      :: D_IN, DD
      TYPE(SIMD_REAL64), INTENT(INOUT)   :: EXC, EXCD, EXCDD
      LOGICAL, INTENT(IN)                :: LLDA
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: D, RS, SK, S, T
      TYPE(SIMD_REAL64)                  :: EXLDA, EXDLDA, EXLDA_SR, EXLDA_LR, EXDLDA_SR, EXDLDA_LR
      TYPE(SIMD_REAL64)                  :: EXWPBE_SR, EXWPBE_LR, EXWPBED_SR, EXWPBED_LR, EXWPBEDD_SR, EXWPBEDD_LR
      TYPE(SIMD_REAL64)                  :: ECLDA, ECDLDA, EC, ECD, ECDD
      TYPE(SIMD_REAL64)                  :: DTHRD, UKFACTOR
      TYPE(SIMD_REAL64)                  :: RHO1, RHO2, DRHO1, DRHO2
      TYPE(SIMD_REAL64)                  :: EXCD1, EXCDD1, EXCD2, EXCDD2, EXCDDA
      TYPE(SIMD_REAL64)                  :: D1, D2, DD1, DD2, FXC, FXCLDA, VXCLDA1, VXCLDA2, FXCD1, FXCD2, FXCDD1, FXCDD2
      TYPE(SIMD_INT64)                   :: IDER
      TYPE(SIMD_MASK_REAL64)             :: MASK1
      INTEGER                            :: II
      LOGICAL                            :: TRUE_FOR_ANY
      LOGICAL, SAVE                      :: INIT_PRED = .TRUE., PRED_1, PRED_2, PRED_3, PRED_4, PRED_5, PRED_6, PRED_7
      REAL(q)                            :: TEMP_1, TEMP_2, TEMP_3, TEMP_4, TEMP_5, TEMP_6, TEMP_7, TEMP_8, TEMP_9
      REAL(q), PARAMETER                 :: DDELTA = 1.0E-4_q
      REAL(q), PARAMETER                 :: THRD = 0.333333333333333333333333_q ! 1._q/3._q
      REAL(q), PARAMETER                 :: P1 = 0.620350490899400086597382_q   ! (3._q/(4._q*pi))**(1._q/3._q)
      REAL(q), PARAMETER                 :: P2 = 3.093667726280135532590521_q   ! (3._q*pi*pi)**(1._q/3._q)
      REAL(q), PARAMETER                 :: P3 = 0.161620459673995492400067_q   ! 1._q/((3._q*pi*pi)**(1._q/3._q))
      REAL(q), PARAMETER                 :: P4 = 1.273239544735162764865777_q   ! (4._q/pi)
      REAL(q), PARAMETER                 :: P5 = 0.251928970342244873403814_q   ! 1._q/sqrt(p2 * p4)
      LOGICAL, PARAMETER                 :: FORCE_PBE = .FALSE.

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK0%X(II)) THEN
            D%X(II) = D_IN%X(II)
         ELSE
            D%X(II) = 1.0E-20_q
         ENDIF
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! LEXCH = 7
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (LEXCH .EQ. 7) THEN
         CALL GGA91_WB_SIMD(D, DD, EXC, EXCD, EXCDD, MASK0)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
         DO II = 0, SIMD_WIDTH_64BIT - 1
            IF (MASK0%X(II)) THEN
               EXC%X(II) = 2.0_q * EXC%X(II) / D%X(II)
               EXCD%X(II) = 2.0_q * EXCD%X(II)
               EXCDD%X(II) = 2.0_q * EXCDD%X(II)
            ENDIF
         ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! LEXCH = 8 | 9 | 14 | 40..50
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ELSE IF ((LEXCH .EQ. 8) .OR. (LEXCH .EQ. 9) .OR. (LEXCH .EQ. 14) .OR. ((LEXCH .GE. 40) .AND. (LEXCH .LE. 50))) THEN
! Perdew Burke Ernzerhof and revised functional
         IF (INIT_PRED) THEN
            PRED_1 = (LDASCREEN .EQ. 0.0_q) .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE
            PRED_2 = LLDA .AND. ((LDASCREEN .EQ. 0.0_q) .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE)
            PRED_3 = LLDA .AND. (.NOT. ((LDASCREEN .EQ. 0.0_q) .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE)) .AND. LUSE_LONGRANGE_HF
            PRED_4 = LLDA .AND. (.NOT. ((LDASCREEN .EQ. 0.0_q) .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE)) .AND. (.NOT. LUSE_LONGRANGE_HF)
            PRED_5 = (.NOT. LLDA) .AND. ((LDASCREEN .EQ. 0.0_q) .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE)
            PRED_6 = (.NOT. LLDA) .AND. (.NOT.((LDASCREEN .EQ. 0.0_q) .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE)) .AND. LUSE_LONGRANGE_HF
            PRED_7 = (.NOT. LLDA) .AND. (.NOT.((LDASCREEN .EQ. 0.0_q) .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE)) .AND. (.NOT. LUSE_LONGRANGE_HF)
            INIT_PRED = .FALSE.
         ENDIF

         TRUE_FOR_ANY = .FALSE.
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY)
         DO II = 0, SIMD_WIDTH_64BIT - 1
            MASK1%X(II) = MASK0%X(II)
            IF (MASK1%X(II)) THEN
               IF (D%X(II) .GT. 0.0_q) THEN
                  TRUE_FOR_ANY = .TRUE.
               ELSE
                  MASK1%X(II) = .FALSE.
                  D%X(II) = 1.0E-20_q
                  EXC%X(II) = 0.0_q
                  EXCD%X(II) = 0.0_q
                  EXCDD%X(II) = 0.0_q
               ENDIF
            ENDIF
         ENDDO

! return of no SIMD lane continues execution
         IF (.NOT. TRUE_FOR_ANY) RETURN

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2)
         DO II = 0, SIMD_WIDTH_64BIT - 1
            ! Perdew Burke Ernzerhof and revised functional
            IF ((LEXCH .EQ. 8) .OR. (LEXCH .EQ. 14)) THEN
               UKFACTOR%X(II) = 1.0_q
            ELSE
               UKFACTOR%X(II) = 0.0_q
            ENDIF
            DTHRD%X(II) = SIMD_POW(D%X(II), THRD)
            TEMP_1 = 1.0_q / DTHRD%X(II)
            RS%X(II) = P1 * TEMP_1
            SK%X(II) = SQRT(P2 * P4 * DTHRD%X(II))
            TEMP_2 = 1.0_q / D%X(II)
            IF(D%X(II) .GT. 1.0E-10_q)THEN
               S%X(II) = DD%X(II) * P3 * TEMP_1 * TEMP_2
               T%X(II) = DD%X(II) * P5 * SQRT(TEMP_1) * TEMP_2
            ELSE
               S%X(II) = 0.0_q
               T%X(II) = 0.0_q
            ENDIF
         ENDDO

! Iann Gerber: range separating in GGA exchange
         IF (PRED_1) THEN
            IF (LEXCH .EQ.14 ) THEN
               CALL EXCHPBESOL_SIMD(D, DTHRD, S, EXLDA, EXC, EXDLDA, EXCD, EXCDD, UKFACTOR, MASK1)
            ELSE
               CALL EXCHPBE_SIMD(D, DTHRD, S, EXLDA, EXC, EXDLDA, EXCD, EXCDD, UKFACTOR, MASK1)
            ENDIF
         ELSE
            CALL CALC_EXCHWPBE_SP_SIMD(D, S, &
            EXLDA, EXDLDA, EXLDA_SR, EXLDA_LR, EXDLDA_SR, EXDLDA_LR, &
            EXWPBE_SR, EXWPBE_LR, EXWPBED_SR, EXWPBED_LR, EXWPBEDD_SR, EXWPBEDD_LR, MASK1)
         ENDIF
! Iann Gerber: end modification
         CALL CORUNSPPBE_SIMD(RS, ECLDA, ECDLDA, SK, T, EC, ECD, ECDD, .TRUE., MASK1)

         IF (PRED_2) THEN
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, TEMP_3, TEMP_4)
            DO II = 0, SIMD_WIDTH_64BIT - 1
               TEMP_1 = 1.0_q / D%X(II)
               TEMP_2 = EC%X(II) * AGGAC + (EXC%X(II) - EXLDA%X(II)) * TEMP_1 * AGGAX + EXLDA%X(II) * TEMP_1 * LDAX + ECLDA%X(II) ! eventuell noch mit ALDAC multi
               TEMP_3 = ECD%X(II) * AGGAC + (EXCD%X(II) - EXDLDA%X(II)) * AGGAX + EXDLDA%X(II) * LDAX + ECDLDA%X(II) ! eventuell noch mit ALDAC multi
               TEMP_4 = ECDD%X(II) * AGGAC + EXCDD%X(II) * AGGAX
               IF (MASK1%X(II)) THEN
                  EXC%X(II) = 2.0_q * TEMP_2
                  EXCD%X(II) = 2.0_q * TEMP_3
                  EXCDD%X(II) = 2.0_q * TEMP_4
               ENDIF
            ENDDO
         ENDIF

         IF (PRED_3) THEN
         ! Iann Gerber's functionals
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, TEMP_3, TEMP_4)
            DO II = 0, SIMD_WIDTH_64BIT - 1
               TEMP_1 = 1.0_q / D%X(II)
               TEMP_2 = EC%X(II) * AGGAC + EXWPBE_LR%X(II) *TEMP_1 * AGGAX + EXWPBE_SR%X(II) * TEMP_1 + ECLD ! ist EXWPBE nur diff
               TEMP_3 = ECD%X(II) * AGGAC + EXWPBED_LR%X(II) * AGGAX + EXWPBED_SR%X(II) + ECDLDA%X(II)
               TEMP_4 = ECDD%X(II) * AGGAC + EXWPBEDD_LR%X(II) * AGGAX + EXWPBEDD_SR%X(II)
               IF (MASK1%X(II)) THEN
                  EXC%X(II) = 2.0_q * TEMP_2
                  EXCD%X(II) = 2.0_q * TEMP_3
                  EXCDD%X(II) = 2.0_q * TEMP_4
               ENDIF
            ENDDO
         ENDIF

         IF (PRED_4) THEN
         ! wPBE
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, TEMP_3, TEMP_4)
            DO II = 0, SIMD_WIDTH_64BIT - 1
               TEMP_1 = 1.0_q / D%X(II)
               TEMP_2 = EC%X(II) * AGGAC + EXWPBE_SR%X(II) * TEMP_1 * AGGAX + EXWPBE_LR%X(II) * TEMP_1 + ECLDA%X(II)
               TEMP_3 = ECD%X(II) * AGGAC + EXWPBED_SR%X(II) * AGGAX + EXWPBED_LR%X(II) + ECDLDA%X(II)
               TEMP_4 = ECDD%X(II) * AGGAC + EXWPBEDD_SR%X(II) * AGGAX + EXWPBEDD_LR%X(II)
               IF (MASK1%X(II)) THEN
                  EXC%X(II) = 2.0_q * TEMP_2
                  EXCD%X(II) = 2.0_q * TEMP_3
                  EXCDD%X(II) = 2.0_q * TEMP_4
               ENDIF
            ENDDO
         ENDIF

         IF (PRED_5) THEN
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, TEMP_3, TEMP_4)
            DO II = 0, SIMD_WIDTH_64BIT - 1
               TEMP_2 = EC%X(II) * AGGAC + (EXC%X(II) - EXLDA%X(II)) / D%X(II) * AGGAX
               TEMP_3 = ECD%X(II) * AGGAC + (EXCD%X(II) - EXDLDA%X(II)) * AGGAX
               TEMP_4 = ECDD%X(II) * AGGAC + EXCDD%X(II) * AGGAX
               IF (MASK1%X(II)) THEN
                  EXC%X(II) = 2.0_q * TEMP_2
                  EXCD%X(II) = 2.0_q * TEMP_3
                  EXCDD%X(II) = 2.0_q * TEMP_4
               ENDIF
            ENDDO
         ENDIF

         IF (PRED_6) THEN
         ! Iann Gerber's functionals
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, TEMP_3, TEMP_4)
            DO II = 0, SIMD_WIDTH_64BIT - 1
               TEMP_1 = 1.0_q / D%X(II)
               TEMP_2 = EC%X(II) * AGGAC + (EXWPBE_LR%X(II) - EXLDA_LR%X(II)) * TEMP_1 * AGGAX + (EXWPBE_SR%X(II) - EXLDA_SR%X(II)) * TEMP_1
               TEMP_3 = ECD%X(II) * AGGAC + (EXWPBED_LR%X(II) - EXDLDA_LR%X(II)) * AGGAX + (EXWPBED_SR%X(II) - EXDLDA_SR%X(II))
               TEMP_4 = ECDD%X(II) * AGGAC + EXWPBEDD_LR%X(II) * AGGAX + EXWPBEDD_SR%X(II)
               IF (MASK1%X(II)) THEN
                  EXC%X(II) = 2.0_q * TEMP_2
                  EXCD%X(II) = 2.0_q * TEMP_3
                  EXCDD%X(II) = 2.0_q * TEMP_4
               ENDIF
            ENDDO
         ENDIF

         IF (PRED_7) THEN
         ! wPBE
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, TEMP_3, TEMP_4)
            DO II = 0, SIMD_WIDTH_64BIT - 1
               TEMP_1 = 1.0_q / D%X(II)
               TEMP_2 = EC%X(II) * AGGAC + (EXWPBE_SR%X(II) - EXLDA_SR%X(II)) * TEMP_1 * AGGAX + (EXWPBE_LR%X(II) - EXLDA_LR%X(II)) * TEMP_1
               TEMP_3 = ECD%X(II) * AGGAC + (EXWPBED_SR%X(II) - EXDLDA_SR%X(II)) * AGGAX + (EXWPBED_LR%X(II) - EXDLDA_LR%X(II))
               TEMP_4 = ECDD%X(II) * AGGAC + EXWPBEDD_SR%X(II) * AGGAX + EXWPBEDD_LR%X(II)
               IF (MASK1%X(II)) THEN
                  EXC%X(II) = 2.0_q * TEMP_2
                  EXCD%X(II) = 2.0_q * TEMP_3
                  EXCDD%X(II) = 2.0_q * TEMP_4
               ENDIF
            ENDDO
         ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! LEXCH = 11 | 12
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ELSE IF ((LEXCH .EQ. 11) .OR. (LEXCH .EQ. 12)) THEN
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
         DO II = 0, SIMD_WIDTH_64BIT - 1
            RHO1%X(II) = 0.5_q * D%X(II)
            RHO2%X(II) = 0.5_q * D%X(II)
            DRHO1%X(II) = 0.5_q * DD%X(II)
            DRHO2%X(II) = 0.5_q * DD%X(II)
         ENDDO

         CALL B3LYPXCS_SIMD(RHO1, RHO2, DRHO1, DRHO2, DD, EXC, EXCD1, EXCDD1, EXCD2, EXCDD2, EXCDDA, MASK0)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
         DO II = 0, SIMD_WIDTH_64BIT - 1        
            EXCD1%X(II) = 2.0_q * EXCD1%X(II)
            EXCD2%X(II) = 2.0_q * EXCD2%X(II)
            EXCDD1%X(II) = 2.0_q * EXCDD1%X(II)
            EXCDD2%X(II) = 2.0_q * EXCDD2%X(II)
            EXCDDA%X(II) = 2.0_q * EXCDDA%X(II)
            IF (MASK0%X(II)) THEN
               EXC%X(II) = 2.0_q * EXC%X(II) / D%X(II)
               EXCD%X(II) = EXCD1%X(II)
               EXCDD%X(II) = EXCDD1%X(II) + EXCDDA%X(II)
            ENDIF
         ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! LEXCH = 13
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ELSE IF (LEXCH .EQ. 13) THEN
         !aem DD sometimes comes in negative.
         !aem Since AM05 assumes that DD actually IS |grad rho|,
         !aem and thus DD always should be positive,
         !aem we need to use ABS(DD) as input to AM05.
         !aem However, my routine needs to give the same symmetries out as
         !aem PBE and PW91, therefor the sign-compensation of EXCDD.
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
         DO II = 0, SIMD_WIDTH_64BIT - 1
            D1%X(II) = 0.5_q * D%X(II)
            D2%X(II) = D1%X(II)
            DD1%X(II) = ABS(0.5_q * DD%X(II))
            DD2%X(II) = DD1%X(II)
            IDER%X(II) = 1
         ENDDO
         
         CALL AM05_SIMD(D1, D2, DD1, DD2, IDER, FXC, FXCLDA, VXCLDA1, VXCLDA2, FXCD1, FXCD2, FXCDD1, FXCDD2, MASK0)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, TEMP_3, TEMP_4, TEMP_5, TEMP_6, TEMP_7, TEMP_8, TEMP_9, D_INV)
         DO II = 0, SIMD_WIDTH_64BIT - 1
            TEMP_1 = EXC%X(II)
            TEMP_2 = EXCD%X(II)
            TEMP_3 = EXCDD%X(II)
            TEMP_4 = FXC%X(II)
            TEMP_5 = FXCD1%X(II)
            TEMP_6 = FXCLDA%X(II)
            TEMP_7 = VXCLDA1%X(II)
            TEMP_8 = DD%X(II)
            TEMP_9 = FXCDD1%X(II)
            D_INV = 1.0_q / D%X(II)
            IF (MASK0%X(II)) THEN
               TEMP_1 = TEMP_4 * D_INV
               TEMP_2 = TEMP_5
               IF (.NOT. LLDA) THEN
                  TEMP_1 = TEMP_1 - TEMP_6 * D_INV
                  TEMP_2 = TEMP_2 - TEMP_7
               ENDIF
               TEMP_3 = TEMP_9
               IF (TEMP_8 .LT. 0.0_q) THEN
                  TEMP_3 = -TEMP_3
               ENDIF
               ! Hartree -> Rydberg conversion
               TEMP_1 = 2.0_q * TEMP_1
               TEMP_2 = 2.0_q * TEMP_2
               TEMP_3 = 2.0_q * TEMP_3
            ENDIF
            EXC%X(II) = TEMP_1
            EXCD%X(II) = TEMP_2
            EXCDD%X(II) = TEMP_3
         ENDDO
      ELSE IF ((LEXCH .GT. 20) .AND. (LEXCH .LT. 25)) THEN
         IF (MASK0%X(II)) THEN
            EXC%X(II) = 0.0_q
            EXCD%X(II) = 0.0_q
            EXCDD%X(II) = 0.0_q
         ENDIF
      ELSE IF ((LEXCH .EQ. 20) .OR. (LEXCH .EQ. 30)) THEN
         IF (MASK0%X(II)) THEN
            EXC%X(II) = 0.0_q
            EXCD%X(II) = 0.0_q
            EXCDD%X(II) = 0.0_q
         ENDIF
      ELSE
         ! not implemented: see GGAEALL
         CALL vtutor%bug("internal ERROR GGAEALL: Wrong LEXCH, scheme not implemented!", __FILE__, __LINE__)
      ENDIF
      RETURN
      END SUBROUTINE GGAALL_SIMD
#endif

!
!  Common interface to all gradient correction routines:
!  returns energy + potential (if defined)
!
      SUBROUTINE GGAEALL(RHO,DRHO,DABGRH,DDRHO,VXC,EXC)
!$ACC ROUTINE SEQ
      USE prec
      USE tutor, ONLY: vtutor
      IMPLICIT REAL(q) (A-H,O-Z)

NOACC CALL vtutor%bug('internal ERROR GGAEALL: Wrong LEXCH, scheme not implemented!', __FILE__, __LINE__)

      RETURN
      END


      SUBROUTINE GGA91_WB(D,DD,EXC,EXCD,EXCDD)
!$ACC ROUTINE SEQ
      USE prec
      USE setexm_struct_def
      IMPLICIT REAL(q) (A-H,O-Z)
      REAL(q),PARAMETER :: PI=3.1415926536_q
      REAL(q),PARAMETER :: THRD=0.333333333333333_q
!$ACC ROUTINE(EXCH2,GCOR,CORGGA1) SEQ
      IF (D<0) THEN
       EXC   = 0._q
       EXCD  = 0._q
       EXCDD = 0._q
       RETURN
      ENDIF

      G = 1.0_q
      RS=(0.75_q/(PI*D))**THRD
      FK=(3._q*PI*PI*D)**THRD
      SK = SQRT(4.0_q*FK/PI)
      IF(D>1.E-10_q)THEN
        S=DD/(D*FK*2._q)
        T=DD/(D*SK*2._q)
      ELSE
        S=0.0_q
        T=0.0_q
      ENDIF
      CALL EXCH2(D,S,EX,EXD,EXDD)
      CALL GCOR(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,1.00_q,RS,EU,EURS)
      EC = 0._q
      ECD = 0._q
      CALL CORGGA1(D,RS,T,EC1,EC1D,EC1DD, &
     &             FK,SK,G,EU,EURS,ECZET)
!      WRITE(*,'(10F14.6)') D,S,EX/D,EXD,EXDD
!      WRITE(*,'(10F14.6)') D,RS,SK,T,(EC+EC1)/D,ECD+EC1D,EC1DD
      EXC   = EX*AGGAX  +(EC+EC1)*AGGAC
      EXCD  = EXD*AGGAX +(ECD+EC1D)*AGGAC
      EXCDD = EXDD*AGGAX+EC1DD*AGGAC
      RETURN
      END

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE GGA91_WB_SIMD(D_IN, DD, EXC, EXCD, EXCDD, MASK0)
      USE prec
      USE setexm_struct_def
      USE simd
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE(SIMD_REAL64), INTENT(IN)      :: D_IN, DD
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EXC, EXCD, EXCDD
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: D, G, RS, FK, SK, S, T
      TYPE(SIMD_REAL64)                  :: EX, EXD, EXDD
      TYPE(SIMD_REAL64)                  :: EU, EURS
      TYPE(SIMD_REAL64)                  :: EC1, EC1D, EC1DD
      TYPE(SIMD_REAL64)                  :: ECZET
      TYPE(SIMD_MASK_REAL64)             :: MASK1
      INTEGER                            :: II
      LOGICAL                            :: TRUE_FOR_ANY
      REAL(q), PARAMETER                 :: PI = 3.1415926536_q
      REAL(q), PARAMETER                 :: THRD = 0.333333333333333_q
      REAL(q), PARAMETER                 :: P1 = 0.620350490899400086597382_q   ! (3._q/(4._q*pi))**(1._q/3._q)
      REAL(q), PARAMETER                 :: P2 = 3.093667726280135532590521_q   ! (3._q*pi*pi)**(1._q/3._q)
      REAL(q), PARAMETER                 :: P4 = 1.273239544735162764865777_q ! (4._q/pi)

      TRUE_FOR_ANY = .FALSE.
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         D%X(II) = 1.0E-20_q
         MASK1%X(II) = MASK0%X(II)
         IF (MASK1%X(II)) THEN
            D%X(II) = D_IN%X(II)
            IF (D_IN%X(II) .GE. 0.0_q) THEN
               TRUE_FOR_ANY = .TRUE.
            ELSE
               MASK1%X(II) = .FALSE.
            ENDIF
            EXC%X(II) = 0.0_q
            EXCD%X(II) = 0.0_q
            EXCDD%X(II) = 0.0_q
         ENDIF
      ENDDO

! return of no SIMD lane continues execution
      IF (.NOT. TRUE_FOR_ANY) RETURN

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(DTHRD)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         G%X(II) = 1.0_q
         DTHRD = SIMD_POW(D%X(II), THRD)
         RS%X(II) = P1 / DTHRD
         FK%X(II) = P2 * DTHRD
         SK%X(II) = SQRT(P4 * FK%X(II))
         IF(D%X(II) .GT. 1.0E-10_q) THEN
            S%X(II) = DD%X(II) / (D%X(II) * FK%X(II) * 2.0_q)
            T%X(II) = DD%X(II) / (D%X(II) * SK%X(II) * 2.0_q)
         ELSE
            S%X(II) = 0.0_q
            T%X(II) = 0.0_q
         ENDIF
      ENDDO

      CALL EXCH2_SIMD(D, S, EX, EXD, EXDD, MASK1)

      CALL GCOR_SIMD(0.0310907_q, 0.21370_q, 7.5957_q, 3.5876_q, 1.6382_q, 0.49294_q, 1.0_q, RS, EU, EURS, MASK1)

      CALL CORGGA1_SIMD(D, RS, T, EC1, EC1D, EC1DD, FK, SK, G, EU, EURS, ECZET, MASK1)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK1%X(II)) THEN
            EXC%X(II) = EX%X(II) * AGGAX + EC1%X(II) * AGGAC
            EXCD%X(II) = EXD%X(II) * AGGAX + EC1D%X(II) * AGGAC
            EXCDD%X(II) = EXDD%X(II) * AGGAX + EC1DD%X(II) * AGGAC
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE GGA91_WB_SIMD
#endif

!************************ SUBROUTINE GGAALL_GRID ***********************
!
! calculate the gradient corrected exchange correlation potential
! on a grid and store the derivates with respect to the density
! and with respect to the gradient
!
!***********************************************************************
#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES) && !defined(_OPENACC)
      SUBROUTINE GGAALL_GRID(DCHARG, DWORKG, DWORK, EXC, NP, OMEGA)
      USE prec
      USE constant
      USE setexm_struct_def
      USE simd
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: NP          ! number of grid points
      REAL(q), INTENT(IN)    :: OMEGA       ! volume
      REAL(q), INTENT(OUT)   :: EXC         ! exchange correlation energy
      REAL(q), INTENT(IN)    :: DCHARG(NP)  ! charge density
      REAL(q), INTENT(INOUT) :: DWORKG(NP)  ! entry: gradient of charge density
                                            ! exit:  derivative of energy w.r.t. gradient
      RGRID, INTENT(OUT)     :: DWORK(NP)   ! exit:  derivative of energy w.r.t. density
!     local
      INTEGER                :: I, II
      REAL(q)                :: RHO, EXCL, DEXC, DVXC
      TYPE(SIMD_REAL64)      :: VEXC, VEXCL, VDEXC, VDVXC, VDCHARG, VDWORKG, VTEMP_1, VTEMP_2
      TYPE(SIMD_MASK_REAL64) :: MASK0

      PROFILING_START('ggaall_grid')

      EXC = 0.0_q

      IF ((LEXCH .EQ. 7) .OR. (LEXCH .EQ. 8) .OR. (LEXCH .EQ. 9) .OR. (LEXCH .EQ. 11) .OR. (LEXCH .EQ. 13) .OR. (LEXCH .EQ. 14)) THEN

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(I, RHO, EXCL, DEXC, DVXC) &
!$OMP PRIVATE(II, VDCHARG, VDWORKG, VEXC, VEXCL, VDEXC, VDVXC, MASK0, VTEMP_1, VTEMP_2) &
!$OMP SHARED(NP, DCHARG, DWORKG, OMEGA, DWORK) &
!$OMP REDUCTION(+ : EXC)
         DO II = 0, SIMD_WIDTH_64BIT - 1
            VEXC%X(II) = 0.0_q
         ENDDO

!$OMP DO SCHEDULE(GUIDED)
         DO I = 1, NP, SIMD_WIDTH_64BIT
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
            DO II = 0, SIMD_WIDTH_64BIT - 1
               IF ((I + II) .LE. NP) THEN
                  MASK0%X(II) = .TRUE.
                  VDCHARG%X(II) = DCHARG(I + II)
                  VDWORKG%X(II) = DWORKG(I + II)
                  VTEMP_1%X(II) = VDCHARG%X(II) * AUTOA3
                  VTEMP_2%X(II) = VDWORKG%X(II) * AUTOA4
               ELSE
                  MASK0%X(II) = .FALSE.
               ENDIF
            ENDDO

            CALL GGAALL_SIMD(VTEMP_1, VTEMP_2, VEXCL, VDEXC, VDVXC, .FALSE., MASK0)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
            DO II = 0, SIMD_WIDTH_64BIT - 1
               IF (MASK0%X(II)) THEN
                  VEXC%X(II) = VEXC%X(II) + VEXCL%X(II) * VDCHARG%X(II) * RYTOEV * OMEGA
!                 DVXC=DVXC*RYTOEV*AUTOA
!                 store d f/ d (|d rho| ) / |d rho|  in DWORK
                  DWORK(I + II) = VDVXC%X(II) * RYTOEV * AUTOA / MAX(VDWORKG%X(II), 1.0E-10_q)
!                 store d f/ d rho  in DWORKG
                  DWORKG(I + II) = VDEXC%X(II) * RYTOEV
               ENDIF
            ENDDO
         ENDDO
!$OMP END DO

         DO II = 0, SIMD_WIDTH_64BIT - 1
            EXC = EXC + VEXC%X(II)
         ENDDO
!$OMP END PARALLEL
      ELSE
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP PRIVATE(I, RHO, EXCL, DEXC, DVXC) &
!$OMP SHARED(NP, DCHARG, DWORKG, OMEGA, DWORK) &
!$OMP REDUCTION(+ : EXC)
         DO I = 1, NP
            RHO = DCHARG(I)
            CALL GGAALL(RHO * AUTOA3, DWORKG(I) * AUTOA4, EXCL, DEXC, DVXC, .FALSE.)
            EXC = EXC + EXCL * RHO * RYTOEV * OMEGA
!           DVXC=DVXC*RYTOEV*AUTOA
!           store d f/ d (|d rho| ) / |d rho|  in DWORK
            DWORK(I) = DVXC * RYTOEV * AUTOA / MAX(DWORKG(I), 1.0E-10_q)
!           store d f/ d rho  in DWORKG
            DWORKG(I) = DEXC * RYTOEV
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      PROFILING_STOP('ggaall_grid')

      RETURN
      END SUBROUTINE
#else
      SUBROUTINE GGAALL_GRID(DCHARG, DWORKG, DWORK, EXC, NP, OMEGA)
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE prec
      USE constant
      IMPLICIT NONE
      INTEGER NP          ! number of grid points
      REAL(q) OMEGA       ! volume
      REAL(q) EXC         ! exchange correlation energy
      REAL(q) DCHARG(NP)  ! charge density
      REAL(q) DWORKG(NP)  ! entry: gradient of charge density
                          ! exit:  derivative of energy w.r.t. gradient
      RGRID   DWORK(NP)   ! exit:  derivative of energy w.r.t. density
    ! local
      INTEGER I
      REAL(q) RHO, EXCL, DEXC, DVXC
!$ACC ROUTINE (GGAALL) SEQ

      PROFILING_START('ggaall_grid')

      EXC=0
#ifdef _OPENACC
      IF (ACC_EXEC_ON) THEN

!$ACC PARALLEL LOOP REDUCTION(+:EXC) GANG VECTOR PRESENT(DCHARG,DWORKG,DWORK) &
!$ACC PRIVATE(RHO,EXCL,DEXC,DVXC) __IF_ASYNC__
      DO I=1,NP
         RHO= DCHARG(I)

         CALL GGAALL(RHO*AUTOA3,DWORKG(I)*AUTOA4,EXCL,DEXC,DVXC,.FALSE.)

         EXC=EXC+EXCL*RHO*RYTOEV*OMEGA
!        DVXC=DVXC*RYTOEV*AUTOA
        !  store d f/ d (|d rho| ) / |d rho|  in DWORK
         DWORK(I)  = DVXC *RYTOEV*AUTOA/ MAX(DWORKG(I),1.E-10_q)
        !  store d f/ d rho  in DWORKG
         DWORKG(I) = DEXC*RYTOEV
      ENDDO

      ELSE
#endif
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(NP, DCHARG, DWORKG, OMEGA, DWORK) &
!$OMP PRIVATE(I, RHO, EXCL, DEXC, DVXC) REDUCTION(+ : EXC)
      DO I=1,NP
         RHO= DCHARG(I)

         CALL GGAALL(RHO*AUTOA3,DWORKG(I)*AUTOA4,EXCL,DEXC,DVXC,.FALSE.)

         EXC=EXC+EXCL*RHO*RYTOEV*OMEGA
!        DVXC=DVXC*RYTOEV*AUTOA
        !  store d f/ d (|d rho| ) / |d rho|  in DWORK
         DWORK(I)  = DVXC *RYTOEV*AUTOA/ MAX(DWORKG(I),1.E-10_q)
        !  store d f/ d rho  in DWORKG
         DWORKG(I) = DEXC*RYTOEV
      ENDDO
!$OMP END PARALLEL DO
#ifdef _OPENACC
      ENDIF
#endif
      PROFILING_STOP('ggaall_grid')

      END SUBROUTINE
#endif

!=======================================================================
!
!  exchange according to Perdew and Wang 1991
!
!=======================================================================

      SUBROUTINE EXCH2(D,S,EX,EXD,EXDD)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
!  gga91 exchange for a spin-unpolarized electronic system
!  input d : density
!  input s:  abs(grad d)/(2*kf*d)
!  output:  exchange energy per electron (ex) and ites derivatives
!           w.r.t. d (exd) and dd (exdd)
      REAL(q),PARAMETER :: A1=0.19645_q,A2=0.27430_q,A3=0.15084_q,A4=100._q
      REAL(q),PARAMETER :: AX=-0.7385588_q,A=7.7956_q,B1=0.004_q
      REAL(q),PARAMETER :: THPITH=3.0936677262801_q
      THRD = 1._q/3._q
      THRD4 = 4._q/3._q
      FAC = AX*D**THRD
      S2 = S*S
      S3 = S2*S
      S4 = S2*S2
      P0 = 1.0_q/SQRT(1.0_q+A*A*S2)
      P1 = LOG(A*S+1.0_q/P0)
      P2 = EXP(-A4*S2)
      P3 = 1.0_q/(1.0_q+A1*S*P1+B1*S4)
      P4 = 1.0_q+A1*S*P1+(A2-A3*P2)*S2
      F = P3*P4
      EX = FAC*(F-1)*D
!     IF (S==0)
!    &         PRINT *,P1,P2,P3,P4,FAC,F
      P5 = 2.0_q*(S*(A2-A3*P2)+A3*A4*S3*P2-2.0_q*B1*S3)
      P6 = (A1*(P1+A*S*P0)+4.0_q*B1*S3)*((A2-A3*P2)*S2-B1*S4)
      FS = (P5*P3-P6*P3*P3)
      EXD = THRD4*FAC*(F-S*FS-1)
      EXDD = AX*FS*0.5_q/THPITH
      RETURN
      END

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE EXCH2_SIMD(D, S, EX, EXD, EXDD, MASK0)
      USE prec
      USE simd
      IMPLICIT REAL(q) (A-H,O-Z)
!  gga91 exchange for a spin-unpolarized electronic system
!  input d : density
!  input s:  abs(grad d)/(2*kf*d)
!  output:  exchange energy per electron (ex) and ites derivatives
!           w.r.t. d (exd) and dd (exdd)
      TYPE(SIMD_REAL64), INTENT(IN)      :: D, S
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EX, EXD, EXDD
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: A1 = 0.19645_q, A2 = 0.27430_q, A3 = 0.15084_q, A4 = 100.0_q
      REAL(q), PARAMETER                 :: A3_A4 = A3 * A4
      REAL(q), PARAMETER                 :: AX = -0.7385588_q, A = 7.7956_q, A_A = A * A
      REAL(q), PARAMETER                 :: B1 = 0.004_q, TWO_B1 = 2.0_q * B1, FOUR_B1 = 4.0_q * B1
      REAL(q), PARAMETER                 :: THPITH = AX * 0.5_q / 3.0936677262801_q
      REAL(q), PARAMETER                 :: THRD = 1.0_q / 3.0_q, THRD4 = 4.0_q * THRD

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(FAC, S2, S3, S4, P0, P1, P2, P7, P3, P8, P4, F, P9, P5, P6, FS)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         FAC = AX * SIMD_POW(D%X(II), THRD)
         S2 = S%X(II) * S%X(II)
         S3 = S2 * S%X(II)
         S4 = S2 * S2
         P0 = 1.0_q / SQRT(1.0_q + A_A * S2)
         P1 = SIMD_LOG(A * S%X(II) + 1.0_q / P0)
         P2 = SIMD_EXP(-A4 * S2)
         P7 = 1.0_q + A1 * S%X(II) * P1
         P3 = 1.0_q / (P7 + B1 * S4)
         P8 = A2 - A3 * P2
         P4 = P7 + P8 * S2
         F = P3 * P4
         P9 = FAC * (F - 1.0_q)
         P5 = 2.0_q * (S%X(II) * P8 + A3_A4 * S3 * P2 - TWO_B1 * S3)
         P6 = (A1 * (P1 + A * S%X(II) * P0) + FOUR_B1 * S3) * (P8 * S2 - B1 * S4)
         FS = (P5 - P6 * P3) * P3
         IF (MASK0%X(II)) THEN
            EX%X(II) = P9 * D%X(II)
            EXD%X(II) = THRD4 * (P9 - FAC * S%X(II) * FS)
            EXDD%X(II) = FS * THPITH
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE EXCH2_SIMD
#endif

!=======================================================================
!
!  GGA contribtuion to correlation Perdew and Wang 1991
!
!=======================================================================

      SUBROUTINE CORGGA1(D,RS,T,EC1,EC1D,EC1DD, &
     &                  FK,SK,G,EC,ECRS,ECZET)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
!  gga91 correlation
!  input rs: seitz radius
!  input t: abs(grad d)/(d*2.*ks)
!  output hn: nonlocal part of correlation energy
!  hnd,hndd : derivatives of hn w.r.t. d and dd
      REAL(q),PARAMETER :: XNU=15.75592_q,CC0=0.004235_q,CX=-0.001667212_q,ALF=0.09_q
      REAL(q),PARAMETER :: C1=0.002568_q,C2=0.023266_q,C3=7.389E-6_q,C4=8.723_q
      REAL(q),PARAMETER :: C5=0.472_q,C6=7.389E-2_q,A4=100.0_q
      REAL(q),PARAMETER :: THRD=0.33333333333333_q,SIXTH7=1.16666666666667_q
      REAL(q),PARAMETER :: PI=3.1415926536_q
      BET = XNU*CC0
      DELT = 2.0_q*ALF/BET
      G3 = G**3
      G4 = G3*G
      PON = -DELT*EC/(G3*BET)
      B = DELT/(EXP(PON)-1.0_q)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      T6 = T4*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.0_q+B*T2
      Q5 = 1.0_q+B*T2+B2*T4
      Q6 = C1+C2*RS+C3*RS2
      Q7 = 1.0_q+C4*RS+C5*RS2+C6*RS3
      CC = -CX + Q6/Q7
      R0 = (SK/FK)**2
      R1 = A4*R0*G4
      COEFF = CC-CC0-3.0_q*CX/7.0_q
      R2 = XNU*COEFF*G3
      R3 = EXP(-R1*T2)
      H0 = G3*(BET/DELT)*LOG(1.0_q+DELT*Q4*T2/Q5)
      H1 = R3*R2*T2
      H = H0+H1
!============================================================
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      H0T = 2.0_q*BET*T*(1.0_q+2.0_q*B*T2)/Q8
      H0B = -BET*T6*(2.0_q*B+B2*T2)/Q8
      H0RS = H0B*B*ECRS*(B+DELT)/BET
      H1T = 2.0_q*R3*R2*T*(1.0_q-R1*T2)
      CCRS = (C2+2._q*C3*RS)/Q7 - Q6*(C4+2._q*C5*RS+3._q*C6*RS2)/Q7**2
      R1RS = 100.0_q*R0/RS
      H1RS = XNU*T2*R3*(CCRS - COEFF*T2*R1RS)
! = = = = = = = = = = = =
      HT = H0T + H1T
      HRS = H0RS + H1RS
      EC1 = D*H
      EC1D = H-THRD*RS*HRS-SIXTH7*T*HT
      EC1DD = 0.5_q*HT/SK
      RETURN
      END

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE CORGGA1_SIMD(D, RS, T, EC1, EC1D, EC1DD, FK, SK, G, EC, ECRS, ECZET, MASK0)
      USE prec
      USE simd
      IMPLICIT REAL(q) (A-H,O-Z)
!  gga91 correlation
!  input rs: seitz radius
!  input t: abs(grad d)/(d*2.*ks)
!  output hn: nonlocal part of correlation energy
!  hnd,hndd : derivatives of hn w.r.t. d and dd
      TYPE(SIMD_REAL64), INTENT(IN)      :: D, RS, T, FK, SK, G, EC, ECRS, ECZET
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EC1, EC1D, EC1DD
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: XNU = 15.75592_q, CC0 = 0.004235_q, CX = -0.001667212_q, ALF = 0.09_q
      REAL(q), PARAMETER                 :: C1 = 0.002568_q, C2 = 0.023266_q, C3 = 7.389E-6_q, C4 = 8.723_q
      REAL(q), PARAMETER                 :: C5 = 0.472_q, C6 = 7.389E-2_q, A4 = 100.0_q
      REAL(q), PARAMETER                 :: THRD = 0.33333333333333_q, SIXTH7 = 1.16666666666667_q
      REAL(q), PARAMETER                 :: PI = 3.1415926536_q
      REAL(q), PARAMETER                 :: BET = XNU * CC0, BET_INV = 1.0_q / BET
      REAL(q), PARAMETER                 :: DELT = 2.0_q * ALF * BET_INV
      REAL(q), PARAMETER                 :: BET_DELT = BET / DELT
      REAL(q), PARAMETER                 :: CX_3_7 = CX * 3.0_q / 7.0_q

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(G3, G4, PON, B, B2, T2, T4, T6, RS2, RS3, Q4, Q5, Q6, Q7, CC, R0, R1, COEFF, R2, R3, H0, H1, H, Q8, H0T, H0B, H0RS, H1T, CCRS, R1RS, H1RS, HT, HRS)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         G3 = SIMD_POW(G%X(II), 3.0_q)
         G4 = G3 * G%X(II)
         PON = -DELT * EC%X(II) * BET_INV / G3
         B = DELT / (SIMD_EXP(PON) - 1.0_q)
         B2 = B * B
         T2 = T%X(II) * T%X(II)
         T4 = T2 * T2
         T6 = T4 * T2
         RS2 = RS%X(II) * RS%X(II)
         RS3 = RS2 * RS%X(II)
         Q4 = 1.0_q + B * T2
         Q5 = Q4 + B2 * T4
         Q6 = C1 + C2 * RS%X(II) + C3 * RS2
         Q7 = 1.0_q / (1.0_q + C4 * RS%X(II) + C5 * RS2 + C6 * RS3)
         CC = -CX + Q6 * Q7
         R0 = SK%X(II) / FK%X(II)
         R0 = R0 * R0
         R1 = A4 * R0 * G4
         COEFF = CC - CC0 - CX_3_7
         R2 = XNU * COEFF * G3
         R3 = SIMD_EXP(-R1 * T2)
         H0 = G3 * (BET_DELT) * LOG(1.0_q + DELT * Q4 * T2 / Q5)
         H1 = R3 * R2 * T2
         H = H0 + H1
!============================================================
         Q8 = 1.0_q / (Q5 * Q5 + DELT * Q4 * Q5 * T2)
         H0T = 2.0_q * BET * T%X(II) * (1.0_q + 2.0_q * B * T2) * Q8
         H0B = -BET * T6 * (2.0_q * B + B2 * T2) * Q8
         H0RS = H0B * B * ECRS%X(II) * (B + DELT) * BET_INV
         H1T = 2.0_q * R3 * R2 *T%X(II) * (1.0_q - R1 * T2)
         CCRS = (C2 + 2.0_q * C3 * RS%X(II)) * Q7 - Q6 * (C4 + 2.0_q * C5 * RS%X(II) + 3.0_q * C6 * RS2) * Q7 * Q7
         R1RS = 100.0_q * R0 / RS%X(II)
         H1RS = XNU * T2 * R3 * (CCRS - COEFF * T2 * R1RS)
!     = = = = = = = = = = = =
         HT = H0T + H1T
         HRS = H0RS + H1RS
         IF (MASK0%X(II)) THEN
            EC1%X(II) = D%X(II) * H
            EC1D%X(II) = H - THRD * RS%X(II) * HRS - SIXTH7 * T%X(II) * HT
            EC1DD%X(II) = 0.5_q * HT / SK%X(II)
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE CORGGA1_SIMD
#endif

!---------------------------------------------------------------------
! The P B E
! this routine was kindly supplied by Bjork Hammer
! it is essentially identical to Burkes routine but includes
! the revised P B E routines
!---------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE EXCHPBE(rho,rhothrd,s,exlda,expbe,exdlda,exd,exdd, &
     &                   ukfactor)
!$ACC ROUTINE SEQ
!----------------------------------------------------------------------
!  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
!  K Burkes modification of PW91 codes, May 14, 1996
!  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  INPUT rho : DENSITY
!  INPUT rhothrd : DENSITY^(1/3)
!  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
! [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
!     {\bf 40},  3399  (1989) (E).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Formulas:
!       e_x[unif]=ax*rho^(4/3)  [LDA]
! ax = -0.75*(3/pi)^(1/3)
!       e_x[PBE]=e_x[unif]*FxPBE(s)
!       FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
! uk, ul defined after [a](13)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      USE prec
      USE setexm_struct_def
      IMPLICIT REAL(q) (A-H,O-Z)
      parameter(thrd=1._q/3._q,thrd4=4._q/3._q)
      parameter(pi=3.14159265358979323846264338327950_q)
      parameter(ax=-0.738558766382022405884230032680836_q)
      parameter(um=0.2195149727645171_q,uk1=0.8040_q,ul1=um/uk1)
!$ACC ROUTINE (CX13) SEQ
!----------------------------------------------------------------------
! construct LDA exchange energy density
      exunif=AX*rhothrd
      exlda=exunif*rho
      exdlda=exunif*thrd4
!----------------------------------------------------------------------
! construct PBE enhancement factor
      S2 = S*S
!----------------------------------------------------------------------
!vdw jk
      IF (LEXCH .eq. 8) THEN
! These are the PBE96 and revPBE98 functionals
! scale uk with a factor
         uk = uk1*ukfactor
         ul = ul1/ukfactor
         P0=1._q+ul*S2
         FxPBE = 1._q+uk-uk/P0
         expbe = exlda*FxPBE
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
         Fs=2._q*um/(P0*P0)
      ELSE IF (LEXCH .eq. 9) THEN
! This is the RPBE functional [Hammer et al, PRB 59, 7413 (1999)]
         P0=exp(-ul1*S2)
         FxPBE = 1._q+uk1*(1.0_q-P0)
         expbe = exlda*FxPBE
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
         Fs=2._q*um*P0
      ELSE IF (LEXCH .eq. 40) THEN
!revPBE
         uk =1.2450_q
         ul = um/uk
         P0=1._q+ul*S2
         FxPBE = 1._q+uk-uk/P0
         expbe = exlda*FxPBE
         Fs=2._q*um/(P0*P0)
      ELSE IF (LEXCH .eq. 41) THEN
!opt PBE
         fix=-0.054732_q
         umm=0.175519_q
         uk=1.04804_q
         ul=umm/uk
         P0=1._q+ul*S2
         P1=exp(-ul*S2)
         FxPBE1 = 1._q+uk-(1._q+fix)*uk/P0
         FxRPBE = fix*uk*P1
         FxPBE = (FxPBE1+FxRPBE)
         expbe = exlda*FxPBE
         Fs=(1._q+fix)*2._q*umm/(P0*P0)-2._q*umm*P1*fix
      ELSE IF (LEXCH .eq. 42) THEN
!this is B88 exchange functional
!use 0.18333333333 0.220000000 to get optB88  exchange
       !  BETA = 0.0042_q
       !  x=s*2._q**(4._q/3._q)*(3._q*3.1415926535)**(1._q/3._q)
       !  arsinh(x) = LOG(  x + SQRT( x * x + 1.0_q))
       !  expbe = -BETA * RHOS**(4.0_q/3.0_q) * X*X/(1.0_q + 6.0_q * BETA * X * LOG(  X + SQRT( X * X + 1.0_q)) )
          sabs=abs(s)
          x_s=7.7955541793_q*sabs
!          P0=0.196447965_q*LOG(x_s+sqrt(x_s*x_s+1.0_q))
          P0=PARAM1*LOG(x_s+sqrt(x_s*x_s+1.0_q))
          P1=1._q+sabs*P0
!          FxPBE=(1._q+0.27429447_q*sabs*sabs/P1)
          FxPBE=(1._q+PARAM2*sabs*sabs/P1)
          expbe=exlda*FxPBE
          Fs=PARAM2*2._q/P1 - sabs*PARAM2/(P1**2)*(P0+sabs*PARAM1*7.7955541793_q/sqrt(1._q+x_s*x_s))
      ELSE IF (LEXCH .eq. 43) THEN
!Becke86MGC
!this is B86MGC with parameters
!PARAM1 is \mu
!PARAM2 is \kappa
!use 0.1234 and 1.00 to get optB86b  exchange
         P0=1._q+PARAM1/PARAM2*S2
         FxPBE = 1._q+PARAM1*S2/(P0**0.8_q)
         expbe = exlda*FxPBE
         Fs=2*PARAM1/(P0**0.8_q)-1.6_q*PARAM1*PARAM1/PARAM2*S2/(P0**(1.8_q))
      ELSE IF (LEXCH .eq. 44) THEN
!PerdewWang86 reparametrized by Murray
         S4 = S2*S2
         S6 = S2*S4
         fftnt=1._q/15._q
         P0=1._q+1.851_q*S2+17.33_q*S4+0.163_q*S6
         FxPBE = P0**fftnt
         expbe = exlda*FxPBE
         Fs=fftnt*FxPBE/P0*(3.7020_q+69.32_q*S2+0.978_q*S4)
!vdw jk
      ELSE IF (LEXCH .eq. 45) THEN
!vdw TbB
!CX13 by Berland and Hyldgaard
         CALL CX13(rho,rhothrd,s,exlda,expbe,exdlda,exd,exdd)
         RETURN
!vdw TbB
      ENDIF
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! calculate the partial derivatives of ex wrt n and |grad(n)|
!  0.3232409194=(3*pi^2)^(-1/3)
      exd =exunif*THRD4*(FxPBE-S2*Fs)
      exdd=0.5_q*ax*0.3232409194_q*S*Fs
      RETURN
      END
!----------------------------------------------------------------------
!----------------------------------------------------------------------

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE EXCHPBE_SIMD(RHO, RHOTHRD, S, EXLDA, EXPBE, EXDLDA, EXD, EXDD, UKFACTOR, MASK0)
!----------------------------------------------------------------------
!  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
!  K Burkes modification of PW91 codes, May 14, 1996
!  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  INPUT rho : DENSITY
!  INPUT rhothrd : DENSITY^(1/3)
!  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
! [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
!     {\bf 40},  3399  (1989) (E).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Formulas:
!       e_x[unif]=ax*rho^(4/3)  [LDA]
! ax = -0.75*(3/pi)^(1/3)
!       e_x[PBE]=e_x[unif]*FxPBE(s)
!       FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
! uk, ul defined after [a](13)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      USE prec
      USE setexm_struct_def
      USE simd
      IMPLICIT REAL(Q) (A-H,O-Z)
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHO, RHOTHRD, S, UKFACTOR
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EXLDA, EXPBE, EXDLDA, EXD, EXDD
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: S2, FXPBE, EXUNIF, FS
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: THRD = 1.0_q / 3.0_q, THRD4 = 4.0_q / 3.0_q
      REAL(q), PARAMETER                 :: PI = 3.14159265358979323846264338327950_q
      REAL(q), PARAMETER                 :: AX = -0.738558766382022405884230032680836_q
      REAL(q), PARAMETER                 :: UM = 0.2195149727645171_q, UK1 = 0.8040_q, UL1 = UM / UK1

!----------------------------------------------------------------------
! construct LDA exchange energy density
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         EXUNIF%X(II) = AX * RHOTHRD%X(II)
         IF (MASK0%X(II)) THEN
            EXLDA%X(II) = EXUNIF%X(II) * RHO%X(II)
            EXDLDA%X(II) = EXUNIF%X(II) * THRD4
         ENDIF
!----------------------------------------------------------------------
! construct PBE enhancement factor
         S2%X(II) = S%X(II) * S%X(II)
!----------------------------------------------------------------------
      ENDDO
!vdw jk
      IF (LEXCH .EQ. 8) THEN
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(UK, UL, P0)
         DO II = 0, SIMD_WIDTH_64BIT - 1
! These are the PBE96 and revPBE98 functionals
! scale uk with a factor
            UK = UK1 * UKFACTOR%X(II)
            UL = UL1 / UKFACTOR%X(II)
            P0 = 1.0_q / (1.0_q + UL * S2%X(II))
            FXPBE%X(II) = 1.0_q + UK - UK * P0
            IF (MASK0%X(II)) THEN
               EXPBE%X(II) = EXLDA%X(II) * FXPBE%X(II)
            ENDIF
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
            FS%X(II) = 2.0_q * UM * P0 * P0
         ENDDO
      ELSE IF (LEXCH .EQ. 9) THEN
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(P0)
         DO II = 0, SIMD_WIDTH_64BIT - 1
! This is the RPBE functional [Hammer et al, PRB 59, 7413 (1999)]
            P0 = SIMD_EXP(-UL1 * S2%X(II))
               FXPBE%X(II) = 1.0_q + UK1 * (1.0_q - P0)
            IF (MASK0%X(II)) THEN
               EXPBE%X(II) = EXLDA%X(II) * FXPBE%X(II)
            ENDIF
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
            FS%X(II) = 2.0_q * UM * P0
         ENDDO
      ELSE IF (LEXCH .EQ. 40) THEN
      ELSE IF (LEXCH .EQ. 41) THEN
      ELSE IF (LEXCH .EQ. 42) THEN
      ELSE IF (LEXCH .EQ. 43) THEN
      ELSE IF (LEXCH .EQ. 44) THEN
      ELSE IF (LEXCH .EQ. 45) THEN
      ENDIF
!vdw jk
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! calculate the partial derivatives of ex wrt n and |grad(n)|
!  0.3232409194=(3*pi^2)^(-1/3)
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK0%X(II)) THEN
            EXD%X(II) = EXUNIF%X(II) * THRD4 * (FXPBE%X(II) - S2%X(II) * FS%X(II))
            EXDD%X(II) = 0.5_q * AX * 0.3232409194_q * S%X(II) * FS%X(II)
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE EXCHPBE_SIMD
!----------------------------------------------------------------------
!----------------------------------------------------------------------
#endif
!######################################################################
!----------------------------------------------------------------------
![TbB 2014-02-24] Added cx13 by Berland and Hyldgaard
!                 (Phys. Rev. B 89, 035412 (2014))
!----------------------------------------------------------------------
      SUBROUTINE CX13(rho,rhothrd,s,exlda,expbe,exdlda,exd,exdd)
!$ACC ROUTINE SEQ
      use prec
      implicit none
      real(q) rho,rhothrd,exlda,expbe,exdlda,exd,exdd
      real(q) thrd,thrd4,pi,ax,mu,kappa,alpha,beta,mulv,a,b,c,m
      real(q) s,S2,S3,S4,S5,S6,G1,G2,G3,Fx,Fs,exunif,as6,base
      real(q) fs_rpw86,df_ds,df_rPW86_ds
      parameter(thrd=1._q/3._q,thrd4=4._q/3._q)
      parameter(pi=3.14159265358979323846264338327950_q)
      parameter(ax=-0.738558766382022405884230032680836_q)
      parameter(m=1.0_q/15.0_q,a=0.1234_q,b=17.33_q,c=0.163_q)
      parameter(mulv=.09434_q,alpha=0.021789_q,beta=1.15_q)
!----------------------------------------------------------------------
! construct LDA exchange energy density
      exunif=AX*rhothrd
      exlda=exunif*rho
      exdlda=exunif*thrd4
!----------------------------------------------------------------------
! construct enhancement factor
      s2 = s*s
      s3 = s2 * s
      s4 = s2 * s2
      s5 = s2 * s3
      s6 = s2 * s2 *s2
      as6 = alpha*s6
! Calculation of energy
      G1 = 1/(1+as6)
      G2 = as6/(beta+as6)
      G3 = alpha*s5/(beta+as6)
      base = (1.0_q+a/m*S2+b*S4+c*S6)
      fs_rpw86=base**m
      Fx = (1+mulv*s2)*G1 + G2*fs_rpw86
      expbe = exlda*Fx
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first derivatives of Fx w.r.t s.
      df_rpw86_ds=s*m*fs_rpw86/base*(2.0_q*a/m+4.0_q*b*S2+6.0_q*c*S4)
      Fs = 2*mulv*s*G1 - 6*alpha*s5*(1+mulv*s2)*G1**2 + &
      G2*df_rpw86_ds + 6*G3*(1-G2)*fs_rpw86
!----------------------------------------------------------------------
! calculate the partial derivatives of ex wrt n and |grad(n)|
!  0.3232409194=(3*pi^2)^(-1/3)
      exd =exunif*THRD4*(Fx-S*Fs)
      exdd=0.5_q*ax*0.3232409194_q*Fs
      RETURN
      END SUBROUTINE CX13


!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE CORunspPBE(RS,EC,VC,sk, &
     &                  T,H,DVC,ecdd,lgga)
!$ACC ROUTINE SEQ
!----------------------------------------------------------------------
!  Official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
!       : lgga=flag to do gga (0=>LSD only)
!       : lpot=flag to do potential (0=>energy only)
!  output: ec=lsd correlation energy from [a]
!        : vcup=lsd up correlation potential
!        : vcdn=lsd dn correlation potential
!        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
!        : dvcup=nonlocal correction to vcup
!        : dvcdn=nonlocal correction to vcdn
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
!$ACC ROUTINE(gcor2) SEQ
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
      CALL gcor2(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS)
      EC = EU
! check for zero energy, immediate return if true
      IF (EC==0.) THEN
        H=0; DVC=0; ecdd=0
        RETURN
      ENDIF
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LSD potential from [c](A1)
! ECRS = dEc/drs [c](A2)
! ECZET=dEc/dzeta [c](A3)
! FZ = dF/dzeta [c](A4)
      ECRS = EURS
      VC = EC -RS*ECRS/3._q
      if (.not.lgga) return
!----------------------------------------------------------------------
! PBE correlation energy
! G=phi(zeta), given after [a](3)
! DELT=bet/gamma
! B=A of [a](8)
      PON=-EC/(gamma)
      B = DELT/(DEXP(PON)-1._q)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1._q+B*T2
      Q5 = 1._q+B*T2+B2*T4
      H = (BET/DELT)*DLOG(1._q+DELT*Q4*T2/Q5)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
      T6 = T4*T2
      RSTHRD = RS/3._q
      FAC = DELT/B+1._q
      BEC = B2*FAC/(BET)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1._q+2._q*B*T2
      hB = -BET*B*T6*(2._q+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      hT = 2._q*BET*Q9/Q8
      DVC = H+HRS-7.0_q*T2*HT/6._q
      ecdd=0.5_q/sk*t*ht
      RETURN
      END
!----------------------------------------------------------------------

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE CORUNSPPBE_SIMD(RS, EC, VC, SK, T, H, DVC, ECDD, LGGA, MASK0)
!----------------------------------------------------------------------
!  Official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
!       : lgga=flag to do gga (0=>LSD only)
!       : lpot=flag to do potential (0=>energy only)
!  output: ec=lsd correlation energy from [a]
!        : vcup=lsd up correlation potential
!        : vcdn=lsd dn correlation potential
!        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
!        : dvcup=nonlocal correction to vcup
!        : dvcdn=nonlocal correction to vcdn
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
      USE setexm_struct_def
      USE simd
      IMPLICIT NONE
      TYPE(SIMD_REAL64), INTENT(IN)      :: RS, SK, T
      TYPE(SIMD_REAL64), INTENT(INOUT)   :: VC
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EC, H, DVC, ECDD
      LOGICAL, INTENT(IN)                :: LGGA
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: ECRS, RTRS
      TYPE(SIMD_MASK_REAL64)             :: MASK1
      REAL(q)                            :: PON, B, B2, T2, T4, T6, RS2, RS3
      REAL(q)                            :: RSTHRD, FAC, BEC, HB, HRS, HT, TEMP
      REAL(q)                            :: Q_0, Q_1, Q_2, Q_3, Q_4, Q_5, Q_8, Q_9
      INTEGER                            :: II
      REAL(q)                            :: TEMP_1, TEMP_2, TEMP_3, TEMP_4
      LOGICAL                            :: TRUE_FOR_ANY
! numbers for use in LSD energy spin-interpolation formula,  [c](9).
      REAL(q), PARAMETER                 :: GAMMA = 0.03109069086965489503494086371273_q ! (1._q-log(2._q))/pi**2
      REAL(q), PARAMETER                 :: BET = 0.06672455060314922_q ! coefficient in gradient expansion for correlation, [a](4).
      REAL(q), PARAMETER                 :: BETSOL = 0.046_q ! coefficient in gradient expansion for correlation, [a](4).
      REAL(q), PARAMETER                 :: DELT = BET / GAMMA
      REAL(q), PARAMETER                 :: DELTSOL = BETSOL / GAMMA
      REAL(q), PARAMETER                 :: THRD = 0.333333333333333333333333_q ! 1._q/3._q
      REAL(q), PARAMETER                 :: P_1 = 1.0_q / GAMMA                 ! 1._q/gamma
      REAL(q), PARAMETER                 :: P_2 = 1.0_q / BET                   ! 1._q/bet
      REAL(q), PARAMETER                 :: P_2SOL = 1.0_q / BETSOL             ! 1._q/betsol
      REAL(q), PARAMETER                 :: P_3 = 1.166666666666666740681535_q  ! 7._q/6._q

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EC=unpolarized LSD correlation energy
! ECRS=dEC/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         RTRS%X(II) = SQRT(RS%X(II))
      ENDDO

      CALL GCOR2_SIMD(0.0310907_q, 0.21370_q, 7.5957_q, 3.5876_q, 1.6382_q, 0.49294_q, RTRS, EC, ECRS, MASK0)

      TRUE_FOR_ANY = .FALSE.
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY)
      DO II = 0, SIMD_WIDTH_64BIT - 1
! check for zero energy, immediate return if true
         MASK1%X(II) = MASK0%X(II)
         IF (MASK1%X(II)) THEN
            IF (EC%X(II) .EQ. 0.0_q) THEN
               MASK1%X(II) = .FALSE.
               H%X(II) = 0.0_q
               DVC%X(II) = 0.0_q
               ECDD%X(II) = 0.0_q
            ENDIF
         ENDIF
         IF (MASK1%X(II)) THEN
! LSD potential from [c](A1)
! ECRS = dEc/drs [c](A2)
! ECZET=dEc/dzeta [c](A3)
! FZ = dF/dzeta [c](A4)
            VC%X(II) = EC%X(II) - RS%X(II) * ECRS%X(II) * THRD
            IF (LGGA) THEN
               TRUE_FOR_ANY = .TRUE.
            ELSE
               MASK1%X(II) = .FALSE.
            ENDIF
         ENDIF
      ENDDO

      IF (.NOT. TRUE_FOR_ANY) RETURN

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(PON, TEMP_1, TEMP_2, TEMP_3, TEMP_4, B, B2, T2, T4, RS2, RS3, Q_4, Q_5, T6, RSTHRD, FAC, BEC, Q_8, Q_9, HB, HRS, HT)
      DO II = 0,SIMD_WIDTH_64BIT - 1
!----------------------------------------------------------------------
! PBE correlation energy
! G=phi(zeta), given after [a](3)
! DELT=bet/gamma
! B=A of [a](8)
         PON = -EC%X(II) * P_1
         TEMP_1  =  SIMD_EXP(PON) - 1.0_q
         IF (LEXCH .EQ. 14) THEN
            TEMP_2 = DELTSOL
            TEMP_3 = P_2SOL
            TEMP_4 = BETSOL
         ELSE
            TEMP_2 = DELT
            TEMP_3 = P_2
            TEMP_4 = BET
         ENDIF
         B = TEMP_2 / TEMP_1
         B2 = B * B
         T2 = T%X(II) * T%X(II)
         T4 = T2 * T2
         RS2 = RS%X(II) * RS%X(II)
         RS3 = RS2 * RS%X(II)
         Q_4 = 1.0_q + B * T2
         Q_5 = 1.0_q + B * T2 + B2 * T4
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
         T6 = T4 * T2
         RSTHRD = RS%X(II) * THRD
         FAC = TEMP_1 + 1.0_q
         BEC = B2 * FAC * TEMP_3
         Q_8 = Q_5 * Q_5 + TEMP_2 * Q_4 * Q_5 * T2
         Q_9 = 1.0_q + 2.0_q * B * T2
         TEMP_1  =  1.0_q / Q_8
         HB = -TEMP_4 * B * T6 * (2.0_q + B * T2)  *  TEMP_1
         HRS = -RSTHRD * HB * BEC * ECRS%X(II)
         HT = 2.0_q * TEMP_4 * Q_9  *  TEMP_1
         IF (MASK1%X(II)) THEN
            H%X(II) = GAMMA * SIMD_LOG(1.0_q + TEMP_2 * Q_4 * T2 / Q_5)
            DVC%X(II) = H%X(II) + HRS - T2 * HT * P_3
            ECDD%X(II) = 0.5_q / SK%X(II) * T%X(II) * HT
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE CORUNSPPBE_SIMD
#endif


!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE CORPBE(RS,ZET,EC,VCUP,VCDN,g,sk, &
     &                  T,H,DVCUP,DVCDN,ecdd,lgga)
!$ACC ROUTINE SEQ
!----------------------------------------------------------------------
!  Official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
!       : lgga=flag to do gga (0=>LSD only)
!       : lpot=flag to do potential (0=>energy only)
!  output: ec=lsd correlation energy from [a]
!        : vcup=lsd up correlation potential
!        : vcdn=lsd dn correlation potential
!        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
!        : dvcup=nonlocal correction to vcup
!        : dvcdn=nonlocal correction to vcdn
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
!$ACC ROUTINE(GCOR2) SEQ
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
      CALL gcor2(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS)
      CALL gcor2(0.01554535_q,0.20548_q,14.1189_q,6.1977_q,3.3662_q, &
     &    0.62517_q,rtRS,EP,EPRS)
      CALL gcor2(0.0168869_q,0.11125_q,10.357_q,3.6231_q,0.88026_q, &
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
      if (.not.lgga) return
!----------------------------------------------------------------------
! PBE correlation energy
! G=phi(zeta), given after [a](3)
! DELT=bet/gamma
! B=A of [a](8)
      G3 = G**3
      PON=-EC/(G3*gamma)
      B = DELT/(DEXP(PON)-1._q)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1._q+B*T2
      Q5 = 1._q+B*T2+B2*T4
      H = G3*(BET/DELT)*DLOG(1._q+DELT*Q4*T2/Q5)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
      G4 = G3*G
      T6 = T4*T2
      RSTHRD = RS/3._q
      GZ=(((1._q+zet)**2+eta)**sixthm- &
     &((1._q-zet)**2+eta)**sixthm)/3._q
      FAC = DELT/B+1._q
      BG = -3._q*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1._q+2._q*B*T2
      hB = -BET*G3*B*T6*(2._q+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      hZ = 3._q*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
      hT = 2._q*BET*G3*Q9/Q8
      COMM = H+HRS-7.0_q*T2*HT/6._q
      PREF = HZ-GZ*T2*HT/G
      COMM = COMM-PREF*ZET
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
      ecdd=0.5_q/(sk*g)*t*ht
      RETURN
      END
!----------------------------------------------------------------------
! slimmed down version of GCOR used in PW91 routines, to interpolate
! LSD correlation energy, as given by (10) of
! J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
! K. Burke, May 11, 1996.
!----------------------------------------------------------------------
      SUBROUTINE GCOR2(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
DOESI IF (rtrs > 1.0E10_q) THEN; GG = 0.0_q; GGRS = 0.0_q; RETURN; END IF
      Q0 = -2._q*A*(1._q+A1*rtrs*rtrs)
      Q1 = 2._q*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1._q+1._q/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2._q*B2+rtrs*(3._q*B3+4._q*B4*rtrs))
      GGRS = -2._q*A*A1*Q2-Q0*Q3/(Q1*(1._q+Q1))
      RETURN
      END

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE GCOR2_SIMD(A, A1, B1, B2, B3, B4, RTRS, GG, GGRS, MASK0)
      USE prec
      USE simd
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE(SIMD_REAL64), INTENT(IN)       :: RTRS
      TYPE(SIMD_REAL64), INTENT(OUT)      :: GG, GGRS
      TYPE(SIMD_MASK_REAL64), INTENT(IN)  :: MASK0
      INTEGER                             :: II

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(Q0, Q1, Q2, Q3)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         Q0 = -2.0_q * A * (1.0_q + A1 * RTRS%X(II) * RTRS%X(II))
         Q1 = 2.0_q * A * RTRS%X(II) * (B1 + RTRS%X(II) * (B2 + RTRS%X(II) * (B3 + B4 * RTRS%X(II))))
         Q2 = SIMD_LOG(1.0_q + 1.0_q / Q1)
         Q3 = A * (B1 / RTRS%X(II) + 2.0_q * B2 + RTRS%X(II) * (3.0_q * B3 + 4.0_q * B4 * RTRS%X(II)))
         IF (MASK0%X(II)) THEN
            GG%X(II) = Q0 * Q2
            GGRS%X(II) = -2.0_q * A * A1 * Q2 - Q0 * Q3 / (Q1 * (1.0_q + Q1))
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE GCOR2_SIMD
#endif

!aem AM05 included
!----------------------------------------------------------------------
! Below follow routines needed for the am05 functional (called 'AM'
! in INCAR). Ann E. Mattsson April 2007, spin added May 2008.
!----------------------------------------------------------------------
! $Header:$
!**********************************************************************
! Armiento Mattsson am05 functional for exchange and correlation
! Version: 3_spin_xscss
!
! input
!   D1   electron spin-up density [bohr**(-3)]
!   D2   electron spin-down density [bohr**(-3)]
!   DD1  abs of gradient of spin-up density (D1) [bohr**(-4)]
!   DD2  abs of gradient of spin-down density (D2) [bohr**(-4)]
!   der  integer: 1 = calculate derivatives (FXCD{1,2} and FXDD{1,2})
!                 0 = don't calculate derivatives
!
! output (FXCLDA, VXCLDA{1,2} needed due to special VASP potential handling)
!   FXC       total exchange-correlation energy density [hartree]
!             = (2*D1*EX[2*D1] + 2*D2*EX[2*D2])/2 + (D1+D2)*EC[D1,D2]
!   FXCLDA    PW LDA exchange-correlation energy density [hartree]
!             = (2*D1*EXLDA[2*D1] + 2*D2*EXLDA[2*D2])/2 + (D1+D2)*ECLDA[D1,D2]
!   VXCLDA1   LDA exchange-correlation potential for spin-up [hartree]
!             = d(FXCLDA[D1,D2])/d(D1)
!   VXCLDA2   LDA exchange-correlation potential for spin-down [hartree]
!             = d(FXCLDA[D1,D2])/d(D2)
!   FXCD1     d(FXC[D1,D2])/d(D1) [hartree]
!   FXCD2     d(FXC[D1,D2])/d(D2) [hartree]
! In this version both the exchange and the correlation parts are
! dependent on spin-up/spin-down gradients (DD1/DD2) separately
!   FXCDD1     d(FXC[D1,D2])/d(DD1) [hartee/(bohr**(-1))]
!   FXCDD2     d(FXC[D1,D2])/d(DD2) [hartee/(bohr**(-1))]
!
! Citation request: when using this functional, please cite:
! "R. Armiento and A. E. Mattsson, PRB 72, 085108 (2005);
!  A. E. Mattsson and R. Armiento, PRB 79, 155101 (2009)."
!
! (The first paper for the AM05 functional, the second for the
! spin-polarized version)
!
! Copyright (c) 2005-2008, Rickard Armiento
!
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject
! to the following conditions:
!
! The above copyright notice and this permission notice shall be
! included in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
! OTHER DEALINGS IN THE SOFTWARE.
!
!**********************************************************************

      SUBROUTINE AM05(D1,D2,DD1,DD2,der, &
                      FXC,FXCLDA,VXCLDA1,VXCLDA2,FXCD1,FXCD2, &
                      FXCDD1,FXCDD2)
!$ACC ROUTINE SEQ
      USE prec
      implicit none

!     ** Input parameters
      REAL(q) D1, D2, DD1, DD2
      INTEGER der

!     ** Output parameters
      REAL(q) FXC, FXCLDA
      REAL(q) VXCLDA1, VXCLDA2
      REAL(q) FXCD1, FXCD2
      REAL(q) FXCDD1, FXCDD2

!     ** Constants and shorthand notation
      REAL(q) pi, g, a, c
      parameter (pi = 3.141592653589793238462643383279502884197_q)
      parameter (g = 0.8098_q, a = 2.804_q)
      parameter (c = 0.7168_q)

!     ** Local variables
      REAL(q) NTOT, KF1, KF2, su, sd, su2, sd2, RS, ZET
      REAL(q) exlda1, exlda2, vxlda1, vxlda2
      REAL(q) eclda, vclda1, vclda2
      REAL(q) fx, fc
      REAL(q) fxd1, fxd2, fxdd1, fxdd2
      REAL(q) fcd1, fcd2, fcdd1, fcdd2
      REAL(q) X1, X2, Xsos1, Xsos2
      REAL(q) F1, F2, Fsos1, Fsos2
      REAL(q) Hx1, Hx2, Hxsos1, Hxsos2
      REAL(q) Hc1, Hc2, Hcsos1, Hcsos2
      REAL(q) szsoz1, szsoz2, denom1, denom2
      REAL(q) zfac1, zfac2, zosn1, zosn2, w1, w2

!aem  ** Local dummies
      REAL(q) ECLDARS, ECZET, ALFC
!$ACC ROUTINE (CORLSD,AM05_LABERTW) SEQ

!     ** Cutoff
      if((D1 .le. 1.e-16_q) .and. (D2 .le. 1.e-16_q)) then
         NTOT=1.e-16_q
      else
         NTOT=D1+D2
      endif

!     *******************
!     LDA PW correlation
!     *******************
      RS = (3._q/(4._q*pi*NTOT))**(1._q/3._q)

!aem  Decided to use the already implemented CORLSD
!aem  which is the Perdew Wang correlation
!aem  (PRB 45, 13244 (1992)) that AM05 needs.
!aem  CORPBE_LDA in xclib.F is a short version of CORLSD
!aem  without the unneeded variables ECLDARS, ECZET, ALFC
      ZET = (D1-D2)/NTOT

      CALL CORLSD(RS,ZET,eclda,vclda1,vclda2,ECLDARS,ECZET,ALFC)

!      WRITE(*,'(2F14.7)') vclda1,vclda2

!     *******************
!     LDA exchange
!     *******************

      vxlda1 = -(3._q*(2._q*D1)/pi)**(1._q/3._q)
      vxlda2 = -(3._q*(2._q*D2)/pi)**(1._q/3._q)
      exlda1 = 3._q/4._q*vxlda1
      exlda2 = 3._q/4._q*vxlda2

!     **********************************
!     LDA exchange-correlation energy
!     and potential
!     **********************************

      FXCLDA = D1*exlda1 + D2*exlda2 + NTOT*eclda
      VXCLDA1 = vxlda1 + vclda1
      VXCLDA2 = vxlda2 + vclda2

!     ************************
!     Gradient contributions
!     ************************

!     ********************
!     Spin up
!     ********************

      KF1 = (3._q*pi*pi*(2._q*D1))**(1._q/3._q)
      if ((2._q*KF1*(2._q*D1)) .ge. 1.e-30_q) then
        su = (2._q*DD1)/(2._q*KF1*(2._q*D1))
      else if (DD1 .le. 1.e-30_q) then
        su = 0._q
      else if (DD1 .ge. 1.e-30_q) then
        su = 1.e30_q
      endif
      su2 = su**2

!     ** Interpolation index and

!     ********************
!     Correlation
!     ********************

!     ** Correlation refinement function Eqn. (12)

      if (su .ge. 1.e12_q) then
              X1 = 0._q
              Hc1 = g
      else if (su .le. 1.e-30_q) then
              X1 = 1._q
              Hc1 = 1._q
      else
          X1 = 1._q/(1._q + a*su2)
          Hc1 = X1 + g*(1._q - X1)
      endif

!     ********************
!     Exchange
!     ********************

!     ** Airy LAA refinement function
      CALL AM05_LABERTW(su**(3._q/2._q)/sqrt(24._q),w1)

!     ** am05_lambertw give back argument if it is < 1.0e-20
!     ** (1.0e-14)^{3/2} = 1.0e-21 => give  low s limit for z/s
!     ** zosn = normalized z/s
      if (su < 1.e-14_q) then
              zosn1 = 1._q
      else
              zosn1 = 24._q**(1._q/3._q)*w1**(2._q/3._q)/su
      end if
      zfac1 = su2*(zosn1*27._q/32._q/pi**2)**2

!     ** denom = denominator of Airy LAA refinement function
      denom1 = 1._q + c*su2*zosn1*(1._q + zfac1)**(1._q/4._q)
!     ** Airy exchange enhancement factor Eqn. (8)
      F1 = (c*su2 + 1._q)/denom1

!     ** Exchange refinement function Eqn. (12)
      Hx1 = X1 + (1._q - X1)*F1

!     ********************
!     Spin down
!     ********************

      KF2 = (3._q*pi*pi*(2._q*D2))**(1._q/3._q)
      if ((2._q*KF2*(2._q*D2)) .ge. 1.e-30_q) then
          sd = (2._q*DD2)/(2._q*KF2*(2._q*D2))
      else if (DD2 .le. 1.e-30_q) then
          sd = 0._q
      else if (DD2 .ge. 1.e-30_q) then
          sd = 1.e30_q
      endif
      sd2 = sd**2

!     ** Interpolation index and

!     ********************
!     Correlation
!     ********************

!     ** Correlation refinement function Eqn. (12)

      if (sd .ge. 1.e12_q) then
          X2 = 0._q
          Hc2 = g
      else if (sd .le.  1.e-30_q) then
          X2 = 1._q
          Hc2 = 1._q
      else
          X2 = 1._q/(1._q + a*sd2)
          Hc2 = X2 + g*(1._q - X2)
      endif

!     ********************
!     Exchange
!     ********************

!     ** Airy LAA refinement function
      CALL AM05_LABERTW(sd**(3._q/2._q)/sqrt(24._q),w2)

!     ** am05_lambertw give back argument if it is < 1.0e-20
!     ** (1.0e-14)^{3/2} = 1.0e-21 => give  low s limit for z/s
!     ** zosn = normalized z/s
      if (sd < 1.e-14_q) then
              zosn2 = 1._q
      else
              zosn2 = 24._q**(1._q/3._q)*w2**(2._q/3._q)/sd
      end if
      zfac2 = sd2*(zosn2*27._q/32._q/pi**2)**2

!     ** denom = denominator of Airy LAA refinement function
      denom2 = 1._q + c*sd2*zosn2*(1._q + zfac2)**(1._q/4._q)
!     ** Airy exchange enhancement factor Eqn. (8)
      F2 = (c*sd2 + 1._q)/denom2

!     ** Exchange refinement function Eqn. (12)
      Hx2 = X2 + (1._q - X2)*F2

!     ********************************
!     Exchange-correlation energy
!     ********************************

!     ** Exchange energy density, Ex = Integrate[fx]
      fx = D1*exlda1*Hx1 + D2*exlda2*Hx2

!     ** Correlation energy density, Ec = Integrate[fc]
      fc = eclda*(D1*Hc1+D2*Hc2)

      FXC = fx + fc

      if (der .eq. 0) return

!     *******************************
!     Potentials
!     *******************************

!     **************************
!     Spin up
!     **************************

      if ((D1 .le. 1.e-16_q) .AND. (su .le. 1.e-30_q)) then
          fxd1 = 0._q
          fxdd1 = 0._q
          fcd1 = vclda1*Hc2 + eclda*(Hc1 - Hc2)
          fcdd1 = 0._q
      else
          if (D1 .le. 1.e-16_q) then
              D1 = 1.e-16_q
          endif

!     ***************************
!         Exchange derivatives
!     ***************************

!     ** Interpolation index derivatives, 1/s dX/ds
          Xsos1 = -2._q*a*X1**2

!     ** Airy LAA refinement function derivatives, 1/s dF/ds
!     ** szsoz = s*(dz/ds)/z
          szsoz1 = 1._q/(1._q + w1)

          Fsos1 = c/denom1**2*(2._q - zosn1* &
     &            ((1._q - c*su2)*(1._q + zfac1)**(1._q/4._q) + &
     &             (1._q + c*su2)*(1._q + 3._q/2._q*zfac1)/ &
     &                      (1._q + zfac1)**(3._q/4._q)*szsoz1))

!     ** Refinement function derivatives, 1/s dHx/ds
!     ** We use that (1 - X) = a*X*s2
          Hxsos1 = (1._q - X1)*Fsos1 - (F1 - 1._q)*Xsos1

          fxd1 = vxlda1*Hx1 - 4._q/3._q*exlda1*su2*Hxsos1
          fxdd1 = exlda1*su*Hxsos1/2._q/KF1

!     *****************************
!        Correlation derivatives
!     *****************************

!     ** Correlation refinement function derivatives, 1/s dF/ds
          Hcsos1 = Xsos1*(1._q - g)

          fcd1 = vclda1/NTOT*(D1*Hc1 + D2*Hc2) &
     &           + eclda*D2/NTOT*(Hc1 - Hc2) &
     &           - 4._q/3._q*eclda*su2*Hcsos1
          fcdd1 = eclda*su*Hcsos1/2._q/KF1
      endif

!     **************************
!     Spin down
!     **************************

      if ((D2 .le. 1.e-16_q) .AND. (sd .le. 1.e-30_q)) then
          fxd2 = 0._q
          fxdd2 = 0._q
          fcd2 = vclda2*Hc1 + eclda*(Hc2 - Hc1)
          fcdd2 = 0._q
      else
          if (D2 .le.  1.e-16_q) then
             D2 = 1.e-16_q
          endif

!     ***************************
!         Exchange derivatives
!     ***************************

!     ** Interpolation index derivatives, 1/s dX/ds
          Xsos2 = -2._q*a*X2**2

!     ** Airy LAA refinement function derivatives, 1/s dF/ds
!     ** szsoz = s*(dz/ds)/z
          szsoz2 = 1._q/(1._q + w2)

          Fsos2 = c/denom2**2*(2._q - zosn2* &
     &            ((1._q - c*sd2)*(1._q + zfac2)**(1._q/4._q) + &
     &             (1._q + c*sd2)*(1._q + 3._q/2._q*zfac2)/ &
     &                      (1._q + zfac2)**(3._q/4._q)*szsoz2))

!     ** Refinement function derivatives, 1/s dHx/ds
!     ** We use that (1 - X) = a*X*s2
          Hxsos2 = (1._q - X2)*Fsos2 - (F2 - 1._q)*Xsos2

          fxd2 = vxlda2*Hx2 - 4._q/3._q*exlda2*sd2*Hxsos2
          fxdd2 = exlda2*sd*Hxsos2/2._q/KF2

!     *****************************
!        Correlation derivatives
!     *****************************

!     ** Correlation refinement function derivatives, 1/s dF/ds
          Hcsos2 = Xsos2*(1._q - g)

          fcd2 = vclda2/NTOT*(D1*Hc1+D2*Hc2) &
     &           + eclda*D1/NTOT*(Hc2 - Hc1) &
     &           - 4._q/3._q*eclda*sd2*Hcsos2
          fcdd2 = eclda*sd*Hcsos2/2._q/KF2
      endif

!     *************************************
!        Exchange-correlation derivatives
!     *************************************

      FXCD1 = fxd1 + fcd1
      FXCD2 = fxd2 + fcd2
      FXCDD1 = fxdd1 + fcdd1
      FXCDD2 = fxdd2 + fcdd2

      return
      end

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE AM05_SIMD(D1_IN, D2_IN, DD1, DD2, DER, FXC, FXCLDA, VXCLDA1, VXCLDA2, FXCD1, FXCD2, FXCDD1, FXCDD2, MASK0)
      USE prec
      USE simd
      IMPLICIT NONE
      TYPE(SIMD_REAL64), INTENT(IN)      :: D1_IN, D2_IN, DD1, DD2
      TYPE(SIMD_INT64), INTENT(IN)       :: DER
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64), INTENT(OUT)     :: FXC, FXCLDA, VXCLDA1, VXCLDA2, FXCD1, FXCD2, FXCDD1, FXCDD2
      TYPE(SIMD_REAL64)                  :: D1, D2
      TYPE(SIMD_REAL64)                  :: NTOT, RS, ZET, ECLDA, VCLDA1, VCLDA2, ECLDARS, ECZET, ALFC
      TYPE(SIMD_REAL64)                  :: KF1, KF2, SU, SD, W1, W2
      TYPE(SIMD_REAL64)                  :: ARG_1, ARG_2
      TYPE(SIMD_REAL64)                  :: ZOSN1, ZOSN2, ZFAC1, ZFAC2, DENOM1, DENOM2, F1, F2, X1, X2
      TYPE(SIMD_REAL64)                  :: HC1, HC2, HX1, HX2, EXLDA1, EXLDA2, VXLDA1, VXLDA2
      TYPE(SIMD_MASK_REAL64)             :: MASK1
      REAL(q)                            :: SU2, SD2, NTOT_INV, KF1_INV, KF2_INV
      REAL(q)                            :: FX, FC
      REAL(q)                            :: FXD1, FXD2, FXDD1, FXDD2
      REAL(q)                            :: FCD1, FCD2, FCDD1, FCDD2
      REAL(q)                            :: XSOS1, XSOS2
      REAL(q)                            :: FSOS1, FSOS2
      REAL(q)                            :: HXSOS1, HXSOS2
      REAL(q)                            :: HCSOS1, HCSOS2
      REAL(q)                            :: SZSOZ1, SZSOZ2
      REAL(q)                            :: TEMP_1, TEMP_2, TEMP_3, TEMP_4
      INTEGER                            :: II
      LOGICAL                            :: TRUE_FOR_ANY
      REAL(q), PARAMETER                 :: PI = 3.141592653589793238462643383279502884197_q
      REAL(q), PARAMETER                 :: G = 0.8098_q, A = 2.804_q, C = 0.7168_q
      REAL(q), PARAMETER                 :: THRD = 1.0_q / 3.0_q
      REAL(q), PARAMETER                 :: THRD2 = 2.0_q / 3.0_q
      REAL(q), PARAMETER                 :: P1 = (3.0_q / (4.0_q * PI)) ** THRD
      REAL(q), PARAMETER                 :: P2 = (6.0_q * PI * PI) ** THRD
      REAL(q), PARAMETER                 :: P3 = 1.0_q / SQRT(24.0_q)
      REAL(q), PARAMETER                 :: P4 = 24.0_q ** THRD
      REAL(q), PARAMETER                 :: P5 = 27.0_q / (32.0_q * PI * PI)
      REAL(q), PARAMETER                 :: P6 = (6.0_q / PI) ** THRD
      REAL(q), PARAMETER                 :: P7 = 4.0_q / 3.0_q
      
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, NTOT_INV)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         D1%X(II) = D1_IN%X(II)
         D2%X(II) = D2_IN%X(II)

!     ** CUTOFF
         IF((D1%X(II) .LE. 1.0E-16_q) .AND. (D2%X(II) .LE. 1.0E-16_q)) THEN
            TEMP_1 = 1.0E-16_q
         ELSE
            TEMP_1 = D1%X(II) + D2%X(II)
         ENDIF

         NTOT%X(II) = TEMP_1
         NTOT_INV = 1.0_q / TEMP_1

!     *******************
!     LDA PW CORRELATION (PART 1)
!     *******************
         RS%X(II) = P1 * SIMD_POW(NTOT_INV, THRD)

!AEM  DECIDED TO USE THE ALREADY IMPLEMENTED CORLSD
!AEM  WHICH IS THE PERDEW WANG CORRELATION
!AEM  (PRB 45, 13244 (1992)) THAT AM05 NEEDS.
!AEM  CORPBE_LDA IN XCLIB.F IS A SHORT VERSION OF CORLSD
!AEM  WITHOUT THE UNNEEDED VARIABLES ECLDARS, ECZET, ALFC
         ZET%X(II) = (D1%X(II) - D2%X(II)) * NTOT_INV

!     ************************
!     GRADIENT CONTRIBUTIONS (PART 1)
!     ************************

!     ********************
!     SPIN UP
!     ********************

         KF1%X(II) = P2 * SIMD_POW(D1%X(II), THRD)
         TEMP_1 = 4.0_q * KF1%X(II) * D1%X(II)
         IF (TEMP_1 .GE. 1.0E-30_q) THEN
            SU%X(II) = (2.0_q * DD1%X(II)) / TEMP_1
         ELSE IF (DD1%X(II) .LE. 1.0E-30_q) THEN
            SU%X(II) = 0.0_q
         ELSE IF (DD1%X(II) .GE. 1.0E-30_q) THEN
            SU%X(II) = 1.0E30_q
         ENDIF

!     ********************
!     SPIN DOWN
!     ********************

         KF2%X(II) = P2 * SIMD_POW(D2%X(II), THRD)
         TEMP_1 = 4.0_q * KF2%X(II) * D2%X(II)
         IF (TEMP_1 .GE. 1.0E-30_q) THEN
            SD%X(II) = (2.0_q * DD2%X(II)) / TEMP_1
         ELSE IF (DD2%X(II) .LE. 1.0E-30_q) THEN
            SD%X(II) = 0.0_q
         ELSE IF (DD2%X(II) .GE. 1.0E-30_q) THEN
            SD%X(II) = 1.0E30_q
         ENDIF

         ARG_1%X(II) = SIMD_POW(SU%X(II), 1.5_q) * P3
         ARG_2%X(II) = SIMD_POW(SD%X(II), 1.5_q) * P3
      ENDDO

!     ************************
!     GRADIENT CONTRIBUTIONS (PART 2)
!     ************************

!     *******************
!     LDA PW CORRELATION
!     *******************
      CALL CORLSD_SIMD(RS, ZET, ECLDA, VCLDA1, VCLDA2, ECLDARS, ECZET, ALFC, MASK0)

!     ********************
!     EXCHANGE
!     ********************

!     ** AIRY LAA REFINEMENT FUNCTION
      CALL AM05_LABERTW_SIMD(ARG_1, W1, MASK0)

!     ** AIRY LAA REFINEMENT FUNCTION
      CALL AM05_LABERTW_SIMD(ARG_2, W2, MASK0)

      TRUE_FOR_ANY = .FALSE.
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY) PRIVATE(SU2, SD2, TEMP_1, TEMP_2, FX, FC)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         SU2 = SU%X(II) * SU%X(II)
         SD2 = SD%X(II) * SD%X(II)

!     ** AM05_LAMBERTW GIVE BACK ARGUMENT IF IT IS < 1.0E-20
!     ** (1.0E-14)^{3/2} = 1.0E-21 => GIVE  LOW S LIMIT FOR Z/S
!     ** ZOSN = NORMALIZED Z/S
         IF (SU%X(II) < 1.0E-14_q) THEN
            TEMP_1 = 1.0_q
         ELSE
            TEMP_1 = P4 * SIMD_POW(W1%X(II), THRD2) / SU%X(II)
         END IF
         ZOSN1%X(II) = TEMP_1
         TEMP_2 = P5 * TEMP_1
         ZFAC1%X(II) = SU2 * TEMP_2 * TEMP_2
         
!     ** DENOM = DENOMINATOR OF AIRY LAA REFINEMENT FUNCTION
         DENOM1%X(II) = 1.0_q / (1.0_q + C * SU2 * TEMP_1 * SIMD_POW(1.0_q + ZFAC1%X(II), 0.25_q))
!     ** AIRY EXCHANGE ENHANCEMENT FACTOR EQN. (8)
         F1%X(II) = (C * SU2 + 1.0_q) * DENOM1%X(II)

!     ********************
!     CORRELATION
!     ********************

!     ** CORRELATION REFINEMENT FUNCTION EQN. (12)

         IF (SU%X(II) .GE. 1.0E12_q) THEN
            TEMP_1 = 0.0_q
            TEMP_2 = G
         ELSE IF (SU%X(II) .LE. 1.0E-30_q) THEN
            TEMP_1 = 1.0_q
            TEMP_2 = 1.0_q
         ELSE
            TEMP_1 = 1.0_q / (1.0_q + A * SU2)
            TEMP_2 = TEMP_1 + G * (1.0_q - TEMP_1)
         ENDIF
         X1%X(II) = TEMP_1
         HC1%X(II) = TEMP_2

!     ** EXCHANGE REFINEMENT FUNCTION EQN. (12)
         HX1%X(II) = TEMP_1 + (1.0_q - TEMP_1) * F1%X(II)
         
!     ** INTERPOLATION INDEX AND

!     ********************
!     CORRELATION
!     ********************

!     ** CORRELATION REFINEMENT FUNCTION EQN. (12)

         IF (SD%X(II) .GE. 1.0E12_q) THEN
            TEMP_1 = 0.0_q
            TEMP_2 = G
         ELSE IF (SD%X(II) .LE.  1.0E-30_q) THEN
            TEMP_1 = 1.0_q
            TEMP_2 = 1.0_q
         ELSE
            TEMP_1 = 1.0_q / (1.0_q + A * SD2)
            TEMP_2 = TEMP_1 + G * (1.0_q - TEMP_1)
         ENDIF
         X2%X(II) = TEMP_1
         HC2%X(II) = TEMP_2

!     ********************
!     EXCHANGE
!     ********************

!     ** AM05_LAMBERTW GIVE BACK ARGUMENT IF IT IS < 1.0E-20
!     ** (1.0E-14)^{3/2} = 1.0E-21 => GIVE  LOW S LIMIT FOR Z/S
!     ** ZOSN = NORMALIZED Z/S
         IF (SD%X(II) < 1.0E-14_q) THEN
            TEMP_1 = 1.0_q
         ELSE
            TEMP_1 = P4 * SIMD_POW(W2%X(II), THRD2) / SD%X(II)
         END IF
         ZOSN2%X(II) = TEMP_1
         TEMP_2 = P5 * TEMP_1
         ZFAC2%X(II) = SD2 * TEMP_2 * TEMP_2
         
!     ** DENOM = DENOMINATOR OF AIRY LAA REFINEMENT FUNCTION
         DENOM2%X(II) = 1.0_q / (1.0_q + C * SD2 * TEMP_1 * SIMD_POW(1.0_q + ZFAC2%X(II), 0.25_q))
!     ** AIRY EXCHANGE ENHANCEMENT FACTOR EQN. (8)
         F2%X(II) = (C * SD2 + 1.0_q) * DENOM2%X(II)

!     ** EXCHANGE REFINEMENT FUNCTION EQN. (12)
         HX2%X(II) = X2%X(II) + (1.0_q - X2%X(II)) * F2%X(II)

!     *******************
!     LDA EXCHANGE
!     *******************

         VXLDA1%X(II) = -P6 * SIMD_POW(D1%X(II), THRD)
         VXLDA2%X(II) = -P6 * SIMD_POW(D2%X(II), THRD)
         EXLDA1%X(II) = 0.75_q * VXLDA1%X(II)
         EXLDA2%X(II) = 0.75_q * VXLDA2%X(II)

!     ********************************
!     EXCHANGE-CORRELATION ENERGY
!     ********************************

!     ** EXCHANGE ENERGY DENSITY, EX = INTEGRATE[FX]
         FX = D1%X(II) * EXLDA1%X(II) * HX1%X(II) + D2%X(II) * EXLDA2%X(II) * HX2%X(II)

!     ** CORRELATION ENERGY DENSITY, EC = INTEGRATE[FC]
         FC = ECLDA%X(II) * (D1%X(II) * HC1%X(II) + D2%X(II) * HC2%X(II))

         IF (MASK0%X(II)) THEN
            FXC%X(II) = FX + FC

!     **********************************
!     LDA EXCHANGE-CORRELATION ENERGY
!     AND POTENTIAL
!     **********************************

            FXCLDA%X(II) = D1%X(II) * EXLDA1%X(II) + D2%X(II) * EXLDA2%X(II) + NTOT%X(II) * ECLDA%X(II)
            VXCLDA1%X(II) = VXLDA1%X(II) + VCLDA1%X(II)
            VXCLDA2%X(II) = VXLDA2%X(II) + VCLDA2%X(II)
         ENDIF

         MASK1%X(II) = MASK0%X(II)
         IF (MASK1%X(II)) THEN
            IF (DER%X(II) .NE. 0) THEN
               TRUE_FOR_ANY = .TRUE.
            ELSE
               MASK1%X(II) = .FALSE.
            ENDIF
         ENDIF
      ENDDO

      IF (.NOT. TRUE_FOR_ANY) RETURN

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(SU2, SD2, NTOT_INV, FXD1, FXDD1, FCD1, FCDD1, XSOS1, SZSOZ1, TEMP_1, TEMP_2, FSOS1, HXSOS1, KF1_INV, HCSOS1, FXD2, FXDD2, FCD2, FCDD2, XSOS2, SZSOZ2, FSOS2, HXSOS2, KF2_INV, HCSOS2)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         SU2 = SU%X(II) * SU%X(II)
         SD2 = SD%X(II) * SD%X(II)
         NTOT_INV = 1.0_q / NTOT%X(II)

!     *******************************
!     POTENTIALS
!     *******************************

!     **************************
!     SPIN UP
!     **************************

         IF ((D1%X(II) .LE. 1.0E-16_q) .AND. (SU%X(II) .LE. 1.0E-30_q)) THEN
            FXD1 = 0.0_q
            FXDD1 = 0.0_q
            FCD1 = VCLDA1%X(II) * HC2%X(II) + ECLDA%X(II) * (HC1%X(II) - HC2%X(II))
            FCDD1 = 0.0_q
         ELSE
            IF (D1%X(II) .LE. 1.0E-16_q) THEN
               D1%X(II) = 1.0E-16_q
            ENDIF

!     ***************************
!         EXCHANGE DERIVATIVES
!     ***************************

!     ** INTERPOLATION INDEX DERIVATIVES, 1/S DX/DS
            XSOS1 = -2.0_q * A * X1%X(II) * X1%X(II)

!     ** AIRY LAA REFINEMENT FUNCTION DERIVATIVES, 1/S DF/DS
!     ** SZSOZ = S*(DZ/DS)/Z
            SZSOZ1 = 1.0_q / (1.0_q + W1%X(II))
            TEMP_1 = SIMD_POW(1.0_q + ZFAC1%X(II), 0.25_q)
            TEMP_2 = TEMP_1 * TEMP_1 * TEMP_1
            FSOS1 = C * DENOM1%X(II) * DENOM1%X(II) * (2.0_q - ZOSN1%X(II) * &
     &           ((1.0_q - C * SU2) * TEMP_1 + &
     &            (1.0_q + C * SU2) * (1.0_q + 1.5_q * ZFAC1%X(II)) / TEMP_2 * SZSOZ1))
            
!     ** REFINEMENT FUNCTION DERIVATIVES, 1/S DHX/DS
!     ** WE USE THAT (1 - X) = A*X*S2
            HXSOS1 = (1.0_q - X1%X(II)) * FSOS1 - (F1%X(II) - 1.0_q) * XSOS1
            KF1_INV = 1.0_q / KF1%X(II)
            FXD1 = VXLDA1%X(II) * HX1%X(II) - P7 * EXLDA1%X(II) * SU2 * HXSOS1
            FXDD1 = 0.5_q * EXLDA1%X(II) * SU%X(II) * HXSOS1 * KF1_INV

!     *****************************
!        CORRELATION DERIVATIVES
!     *****************************

!     ** CORRELATION REFINEMENT FUNCTION DERIVATIVES, 1/S DF/DS
            HCSOS1 = XSOS1 * (1.0_q - G)
            FCD1 = VCLDA1%X(II) * NTOT_INV * (D1%X(II) * HC1%X(II) + D2%X(II) * HC2%X(II)) &
     &           + ECLDA%X(II) * D2%X(II) * NTOT_INV * (HC1%X(II) - HC2%X(II)) &
     &           - P7 * ECLDA%X(II) * SU2 * HCSOS1
            FCDD1 = 0.5_q * ECLDA%X(II) * SU%X(II) * HCSOS1 * KF1_INV
         ENDIF

!     **************************
!     SPIN DOWN
!     **************************

         IF ((D2%X(II) .LE. 1.0E-16_q) .AND. (SD%X(II) .LE. 1.0E-30_q)) THEN
            FXD2 = 0.0_q
            FXDD2 = 0.0_q
            FCD2 = VCLDA2%X(II) * HC1%X(II) + ECLDA%X(II) * (HC2%X(II) - HC1%X(II))
            FCDD2 = 0.0_q
         ELSE
            IF (D2%X(II) .LE.  1.0E-16_q) THEN
               D2%X(II) = 1.0E-16_q
            ENDIF
            
!     ***************************
!         EXCHANGE DERIVATIVES
!     ***************************

!     ** INTERPOLATION INDEX DERIVATIVES, 1/S DX/DS
            XSOS2 = -2.0_q * A * X2%X(II) * X2%X(II)

!     ** AIRY LAA REFINEMENT FUNCTION DERIVATIVES, 1/S DF/DS
!     ** SZSOZ = S*(DZ/DS)/Z
            SZSOZ2 = 1.0_q / (1.0_q + W2%X(II))
            TEMP_1 = SIMD_POW(1.0_q + ZFAC2%X(II), 0.25_q)
            TEMP_2 = TEMP_1 * TEMP_1 * TEMP_1
            FSOS2 = C * DENOM2%X(II) * DENOM2%X(II) * (2.0_q - ZOSN2%X(II) * &
     &           ((1.0_q - C * SD2) * TEMP_1 + &
     &            (1.0_q + C * SD2) * (1.0_q + 1.5_q * ZFAC2%X(II)) / TEMP_2 * SZSOZ2))

!     ** REFINEMENT FUNCTION DERIVATIVES, 1/S DHX/DS
!     ** WE USE THAT (1 - X) = A*X*S2
            HXSOS2 = (1.0_q - X2%X(II)) * FSOS2 - (F2%X(II) - 1.0_q) * XSOS2
            KF2_INV = 1.0_q / KF2%X(II)
            FXD2 = VXLDA2%X(II) * HX2%X(II) - P7 * EXLDA2%X(II) * SD2 * HXSOS2
            FXDD2 = 0.5_q * EXLDA2%X(II) * SD%X(II) * HXSOS2 * KF2_INV

!     *****************************
!        CORRELATION DERIVATIVES
!     *****************************

!     ** CORRELATION REFINEMENT FUNCTION DERIVATIVES, 1/S DF/DS
            HCSOS2 = XSOS2 * (1.0_q - G)
            FCD2 = VCLDA2%X(II) * NTOT_INV * (D1%X(II) * HC1%X(II) + D2%X(II) * HC2%X(II)) &
     &           + ECLDA%X(II) * D1%X(II) * NTOT_INV *(HC2%X(II) - HC1%X(II)) &
     &           - P7 * ECLDA%X(II) * SD2 * HCSOS2
            FCDD2 = 0.5_q * ECLDA%X(II) * SD%X(II) * HCSOS2 * KF2_INV
         ENDIF

!     *************************************
!        EXCHANGE-CORRELATION DERIVATIVES
!     *************************************

         IF (MASK1%X(II)) THEN
            FXCD1%X(II) = FXD1 + FCD1
            FXCD2%X(II) = FXD2 + FCD2
            FXCDD1%X(II) = FXDD1 + FCDD1
            FXCDD2%X(II) = FXDD2 + FCDD2
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE AM05_SIMD
#endif

!     ***********************************************
!       LambertW function.
!
!       Corless, Gonnet, Hare, Jeffrey, and Knuth (1996),
!         Adv. in Comp. Math. 5(4):329-359.
!       Implementation approach loosely inspired by the
!       GNU Octave version by N. N. Schraudolph, but this
!       implementation is only for real values and
!       principal branch.
!
!       Copyright (c) 2005, Rickard Armiento
!       All rights reserved.
!     ***********************************************

      SUBROUTINE AM05_LABERTW(z,result)
!$ACC ROUTINE SEQ
      USE prec
      USE tutor, ONLY : vtutor
      implicit none

!     input
      REAL(q) z
!     output
      REAL(q) result
!     local variables
      REAL(q) e,t,p
      INTEGER i

!     ** If z too low, go with the first term of the power expansion, z
      if( z .lt. 1.E-20_q) then
         result = z
         return
      endif

      e = exp(1._q)

!     ** Inital guess
      if( abs(z + 1._q/e) .gt. 1.45_q ) then
!        ** Asymptotic expansion at 0 and Inf
         result = log(z)
         result = result - log(result)
      else
!        ** Series expansion about -1/e to first order
         result = 1._q*sqrt(2._q*e*z + 2._q) - 1._q
      endif

!     ** Find result through iteration
      do i=1,10
         p = exp(result)
         t = result*p - z
         if( result .ne. -1._q ) then
            t = t/(p*(result + 1._q) - &
     &           0.5_q*(result + 2._q)*t/(result + 1._q))
         else
            t = 0._q
         endif
         result = result - t
         if(abs(t) < (2.48_q*1.E-14_q)*(1._q + abs(result))) then
            return
         endif
      enddo
!     ** This should never happen!;
NOACC call vtutor%error('am05_labertw: iteration limit reached. Should never happen: execution aborted.')
      return
      end

!aem AM05 included

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE AM05_LABERTW_SIMD(Z, RESULT, MASK0)
      USE prec
      USE simd
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      TYPE(SIMD_REAL64), INTENT(IN)      :: Z
      TYPE(SIMD_REAL64), INTENT(OUT)     :: RESULT
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_MASK_REAL64)             :: MASK1
      REAL(q)                            :: TEMP_1, T, P
      INTEGER                            :: I, II
      LOGICAL                            :: TRUE_FOR_ANY
      REAL(q), PARAMETER                 :: E = EXP(1.0_q), E_INV = 1.0_q / E
      
!     ** IF Z TOO LOW, GO WITH THE FIRST TERM OF THE POWER EXPANSION, Z
      TRUE_FOR_ANY = .FALSE.
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY) PRIVATE(TEMP_1)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         TEMP_1 = Z%X(II)
         IF (MASK0%X(II)) THEN
            IF (TEMP_1 .GE. 1.0E-20_q) THEN
               MASK1%X(II) = .TRUE.
               TRUE_FOR_ANY = .TRUE.
            ELSE
               MASK1%X(II) = .FALSE.
               RESULT%X(II) = TEMP_1
            ENDIF
         ENDIF
      ENDDO
      
      IF (.NOT. TRUE_FOR_ANY) RETURN

!     $OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1)
      DO II = 0, SIMD_WIDTH_64BIT - 1
!     ** INITAL GUESS
         IF (ABS(Z%X(II) + E_INV) .GT. 1.45_q) THEN
!        ** ASYMPTOTIC EXPANSION AT 0 AND INF
            TEMP_1 = SIMD_LOG(Z%X(II))
            TEMP_1 = TEMP_1 - SIMD_LOG(TEMP_1)
         ELSE
!        ** SERIES EXPANSION ABOUT -1/E TO FIRST ORDER
            TEMP_1 = 1.0_q * SQRT(2.0_q * E * Z%X(II) + 2.0_q) - 1.0_q
         ENDIF
         IF (MASK1%X(II)) THEN
            RESULT%X(II) = TEMP_1
         ENDIF
      ENDDO

!     ** FIND RESULT THROUGH ITERATION
      DO I = 1, 10
         TRUE_FOR_ANY = .FALSE.
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY) PRIVATE(P, T)
         DO II = 0, SIMD_WIDTH_64BIT - 1
            P = SIMD_EXP(RESULT%X(II))
            T = RESULT%X(II) * P - Z%X(II)
            IF (RESULT%X(II) .NE. -1.0_q ) THEN
               T =  T / (P * (RESULT%X(II) + 1.0_q) - 0.5_q * (RESULT%X(II) + 2.0_q) * T / (RESULT%X(II) + 1.0_q))
            ELSE
               T = 0.0_q
            ENDIF
            IF (MASK1%X(II)) THEN
               RESULT%X(II) = RESULT%X(II) - T
               IF (ABS(T) .GE. (2.48_q * 1.0E-14_q) * (1.0_q + ABS(RESULT%X(II)))) THEN
                  TRUE_FOR_ANY = .TRUE.
               ELSE
                  MASK1%X(II) = .FALSE.
               ENDIF
            ENDIF
         ENDDO
         IF (.NOT. TRUE_FOR_ANY) RETURN
      ENDDO
!     ** THIS SHOULD NEVER HAPPEN!;
      CALL vtutor%error("AM05_LABERTW_SIMD: ITERATION LIMIT REACHED. \n SHOULD NEVER HAPPEN: EXECUTION &
         &ABORTED.")
      RETURN
      END SUBROUTINE AM05_LABERTW_SIMD
#endif


!************************ SUBROUTINE GGASPIN *****************************
!
!  switch between different GGAs
!  presently only PBE and RPBE are implemented efficiently
!  (i.e. d exc / d rho and  d exc / d | grad rho | are calculated
!  directly
!  for other GGA functional finite differences are used to calculate
!  the required derivatives
!  LLDA allows to include the LDA contribution directly in this
!  routine.
!  This only works for PBE and RPBE and fails in all other cases
!
!***********************************************************************

      SUBROUTINE GGASPIN(D1,D2,DD1,DD2,DDA,EXC, &
         excd1,excd2,excq1,excq2,ecq,LLDA)
!$ACC ROUTINE SEQ

!     D1   density up
!     D2   density down
!     DD1  |gradient of density up|
!     DD2  |gradient of density down|
!     DDA  |gradient of the total density|
!     LLDA add lda contributions

      USE prec
      USE constant
      USE setexm_struct_def
      USE tutor, ONLY : vtutor
      IMPLICIT REAL(q) (A-H,O-Z)
#ifdef USELIBXC
      REAL(q) XC_DENS(2),XC_SIGMA(3),XC_EXC(1),XC_VRHO(2),XC_VSIGMA(3)
      REAL(q) XC_SIGMAZERO(3),XC_EXCZERO(1),XC_VRHOZERO(2),XC_VSIGMAZERO(3)
#endif
      LOGICAL LLDA
      PARAMETER (THRD=1._q/3._q,THRD4=4._q/3._q)
      LOGICAL:: FORCE_PBE=.FALSE.
!$ACC ROUTINE ( GGAXCS, EXCHPBE, CALC_EXCHWPBE_SP, CORPBE, B3LYPXCS, AM05, EXCHPBESOL, &
!$ACC& EXCHPBESOL, CALC_EXCHWPBEsol_SP, corpbesol ) SEQ
#ifdef USELIBXC
!$ACC ROUTINE ( XC_F03_LDA_VXC, XC_F03_LDA_EXC_VXC, XC_F03_GGA_VXC, XC_F03_GGA_EXC_VXC) SEQ
#endif
      ECQ=0

      IF (LEXCH==7) THEN
         CALL GGAXCS(D1,D2,DD1,DD2,EXC,EXCD1,EXCD2,EXCQ1,EXCQ2)
!
         D = D1 + D2
         EXC    = 2._q* EXC / D
         EXCD1  = 2._q* EXCD1
         EXCD2  = 2._q* EXCD2
         EXCQ1  = 2._q* EXCQ1
         EXCQ2  = 2._q* EXCQ2
      ELSE IF (LEXCH==8 .OR. LEXCH==9 .OR. (LEXCH.ge.18 .AND. LEXCH.le.50)) THEN
        IF (LEXCH==8) THEN
           ukfactor=1.0_q
        ELSE
           ukfactor=0.0_q
        ENDIF

        D=2*D1

#undef exchange_RHOU_RHOD

#ifdef  exchange_RHOU_RHOD
        D=D1+D2
#endif
        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN

           S=DD1/(D*FK)
#ifdef  exchange_RHOU_RHOD
           S=(DD1+DD2)/2/(D*FK)
#endif

        ELSE
           S=0.0_q
        ENDIF
        IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
           CALL EXCHPBE(D,DTHRD,S,EXLDA1,EXC1,EXDLDA1,EXCD1,EXCQ1, &
          &     ukfactor)
        ELSE
           CALL CALC_EXCHWPBE_SP(D,S, &
          &  EXLDA1,EXDLDA1,EXLDA_SR1,EXLDA_LR1,EXDLDA_SR1,EXDLDA_LR1, &
          &  EXWPBE_SR1,EXWPBE_LR1,EXWPBED_SR1,EXWPBED_LR1,EXWPBEDD_SR1,EXWPBEDD_LR1)
        ENDIF
!        WRITE(*,'(10F14.7)') D,S,(EXC1-EXLDA1)/D,EXCD1-EXDLDA1,EXCQ1
!        EXLDA1=0; EXC1=0; EXDLDA1=0; EXCD1=0; EXCQ1=0

        D=2*D2
#ifdef  exchange_RHOU_RHOD
        D=D1+D2
#endif

        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN
           S=DD2/(D*FK)
#ifdef  exchange_RHOU_RHOD
           S=(DD1+DD2)/2/(D*FK)
#endif
        ELSE
           S=0.0_q
        ENDIF

        IF (LDASCREEN==0._q.OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
           CALL EXCHPBE(D,DTHRD,S,EXLDA2,EXC2,EXDLDA2,EXCD2,EXCQ2, &
          &     ukfactor)
        ELSE
           CALL CALC_EXCHWPBE_SP(D,S, &
          &  EXLDA2,EXDLDA2,EXLDA_SR2,EXLDA_LR2,EXDLDA_SR2,EXDLDA_LR2, &
          &  EXWPBE_SR2,EXWPBE_LR2,EXWPBED_SR2,EXWPBED_LR2,EXWPBEDD_SR2,EXWPBEDD_LR2)
        ENDIF
!        EXLDA2=0; EXC2=0; EXDLDA2=0; EXCD2=0; EXCQ2=0
!        WRITE(*,'(10F14.7)') D,S,(EXC2-EXLDA2)/D,EXCD2-EXDLDA2,EXCQ2

        D=D1+D2
        DTHRD=exp(log(D)*THRD)
        RS=(0.75_q/PI)**THRD/DTHRD
        ZETA=(D1-D2)/D
        ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)

        FK=(3._q*PI*PI)**THRD*DTHRD
        SK = SQRT(4.0_q*FK/PI)

        G = (exp((2*THRD)*log(1._q+ZETA)) &
                +exp((2*THRD)*log(1._q-ZETA)))/2._q
        T = DDA/(D*2._q*SK*G)

        CALL corpbe(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA,G,SK, &
                     T,EC,ECD1,ECD2,ECQ,.TRUE.)
!        ECLDA=0 ; ECD1LDA=0 ; ECD2LDA=0; EC=0 ; ECD1=0 ; ECD2=0 ; ECQ=0
!        WRITE(*,'(10F14.7)') D,RS,SK,T,EC,ECD1,ECD2,ECQ

        IF (LLDA) THEN
        ! Add LDA contributions as well
           IF (LDASCREEN==0._q.OR. LUSE_THOMAS_FERMI  .OR. FORCE_PBE) THEN
              EXC  =EC  *AGGAC +(EXC1-EXLDA1+EXC2-EXLDA2)/(2*D)*AGGAX+(EXLDA1+EXLDA2)/(2*D)*LDAX+ECLDA
              EXCD1=ECD1*AGGAC +(EXCD1-EXDLDA1)*AGGAX +EXDLDA1*LDAX+ECD1LDA
              EXCD2=ECD2*AGGAC +(EXCD2-EXDLDA2)*AGGAX +EXDLDA2*LDAX+ECD2LDA
              EXCQ1=EXCQ1 *AGGAX
              EXCQ2=EXCQ2 *AGGAX
              ECQ  =ECQ *AGGAC
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
              ! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+(EXWPBE_LR1+EXWPBE_LR2)/(2*D)*AGGAX +(EXWPBE_SR1+EXWPBE_SR2)/(2*D)+ECLDA
                 EXCD1=ECD1*AGGAC+ EXWPBED_LR1*AGGAX +EXWPBED_SR1+ECD1LDA
                 EXCD2=ECD2*AGGAC+ EXWPBED_LR2*AGGAX +EXWPBED_SR2+ECD2LDA
                 EXCQ1=EXWPBEDD_LR1*AGGAX+EXWPBEDD_SR1
                 EXCQ2=EXWPBEDD_LR2*AGGAX+EXWPBEDD_SR2
                 ECQ  =ECQ*AGGAC
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC  *AGGAC+(EXWPBE_LR1+EXWPBE_LR2)/(2*D)*AGGAX+ECLDA
                 EXCD1=ECD1*AGGAC+ EXWPBED_LR1*AGGAX+ECD1LDA
                 EXCD2=ECD2*AGGAC+ EXWPBED_LR2*AGGAX+ECD2LDA
                 EXCQ1=EXWPBEDD_LR1*AGGAX
                 EXCQ2=EXWPBEDD_LR2*AGGAX
                 ECQ =ECQ*AGGAC
              ELSE
              ! wPBE
                 EXC  =EC  *AGGAC+(EXWPBE_SR1+EXWPBE_SR2)/(2*D)*AGGAX +(EXWPBE_LR1+EXWPBE_LR2)/(2*D)+ECLDA
                 EXCD1=ECD1*AGGAC+ EXWPBED_SR1*AGGAX +EXWPBED_LR1+ECD1LDA
                 EXCD2=ECD2*AGGAC+ EXWPBED_SR2*AGGAX +EXWPBED_LR2+ECD2LDA
                 EXCQ1=EXWPBEDD_SR1*AGGAX+EXWPBEDD_LR1
                 EXCQ2=EXWPBEDD_SR2*AGGAX+EXWPBEDD_LR2
                 ECQ  =ECQ*AGGAC
              ENDIF
           ENDIF
        ELSE
        ! Do not add LDA contribution
           IF (LDASCREEN==0._q.OR. LUSE_THOMAS_FERMI  .OR. FORCE_PBE) THEN
              EXC  =EC*  AGGAC +(EXC1-EXLDA1+EXC2-EXLDA2)/(2*D)*AGGAX
              EXCD1=ECD1*AGGAC +(EXCD1-EXDLDA1)*AGGAX
              EXCD2=ECD2*AGGAC +(EXCD2-EXDLDA2)*AGGAX
              EXCQ1=EXCQ1 *AGGAX
              EXCQ2=EXCQ2 *AGGAX
              ECQ  =ECQ *AGGAC
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
              ! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+(EXWPBE_LR1-EXLDA_LR1+EXWPBE_LR2-EXLDA_LR2)/(2*D)*AGGAX+ &
                &                 (EXWPBE_SR1-EXLDA_SR1+EXWPBE_SR2-EXLDA_SR2)/(2*D)
                 EXCD1=ECD1*AGGAC+(EXWPBED_LR1-EXDLDA_LR1)*AGGAX+(EXWPBED_SR1-EXDLDA_SR1)
                 EXCD2=ECD2*AGGAC+(EXWPBED_LR2-EXDLDA_LR2)*AGGAX+(EXWPBED_SR2-EXDLDA_SR2)
                 EXCQ1=EXWPBEDD_LR1*AGGAX+EXWPBEDD_SR1
                 EXCQ2=EXWPBEDD_LR2*AGGAX+EXWPBEDD_SR2
                 ECQ  =ECQ *AGGAC
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC  *AGGAC+(EXWPBE_LR1-EXLDA_LR1+EXWPBE_LR2-EXLDA_LR2)/(2*D)*AGGAX
                 EXCD1=ECD1*AGGAC+(EXWPBED_LR1-EXDLDA_LR1)*AGGAX
                 EXCD2=ECD2*AGGAC+(EXWPBED_LR2-EXDLDA_LR2)*AGGAX
                 EXCQ1=EXWPBEDD_LR1*AGGAX
                 EXCQ2=EXWPBEDD_LR2*AGGAX
                 ECQ =ECQ *AGGAC
              ELSE
              ! wPBE
                 EXC  =EC  *AGGAC+(EXWPBE_SR1-EXLDA_SR1+EXWPBE_SR2-EXLDA_SR2)/(2*D)*AGGAX+ &
                &                 (EXWPBE_LR1-EXLDA_LR1+EXWPBE_LR2-EXLDA_LR2)/(2*D)
                 EXCD1=ECD1*AGGAC+(EXWPBED_SR1-EXDLDA_SR1)*AGGAX+(EXWPBED_LR1-EXDLDA_LR1)
                 EXCD2=ECD2*AGGAC+(EXWPBED_SR2-EXDLDA_SR2)*AGGAX+(EXWPBED_LR2-EXDLDA_LR2)
                 EXCQ1=EXWPBEDD_SR1*AGGAX+EXWPBEDD_LR1
                 EXCQ2=EXWPBEDD_SR2*AGGAX+EXWPBEDD_LR2
                 ECQ  =ECQ *AGGAC
              ENDIF
           ENDIF
        ENDIF
        ! Hartree -> Rydberg conversion
        EXC    = 2*EXC
        EXCD1  = 2*EXCD1
        EXCD2  = 2*EXCD2
        EXCQ1  = 2*EXCQ1
        EXCQ2  = 2*EXCQ2
        ECQ    = 2*ECQ
! jP
      ELSE IF (LEXCH==11 .OR. LEXCH==12) THEN

        CALL  B3LYPXCS(D1,D2,DD1,DD2,DDA,EXC,EXCD1,EXCDD1,EXCD2,EXCDD2,EXCDDA)

        D = D1 + D2
        EXC    = 2.0_q*EXC/D
        EXCD1  = 2.0_q*EXCD1
        EXCD2  = 2.0_q*EXCD2
        EXCQ1  = 2.0_q*EXCDD1
        EXCQ2  = 2.0_q*EXCDD2
        ECQ    = 2.0_q*EXCDDA

!       WRITE(91,'(7F12.7)') D1, D2, EXC, EXCD1, EXCD2, EXCDD1, EXCDD2
! jP
!aem AM05 spin formulation added
      ELSE IF (LEXCH==13) THEN

        !aem DD sometimes comes in negative.
        !aem Since AM05 assumes that DD actually IS |grad rho|,
        !aem and thus DD always should be positive,
        !aem we need to use ABS(DD) as input to AM05.
        !aem However, my routine needs to give the same symmetries out as
        !aem PBE and PW91, therefor the sign-compensation of EXCDD.

        ider = 1
        CALL AM05(D1,D2,ABS(DD1),ABS(DD2),ider, &
                  FXC,FXCLDA,VXCLDA1,VXCLDA2,FXCD1,FXCD2, &
                  FXCDD1,FXCDD2)

        D = D1 + D2
        IF (LLDA) THEN
          EXC = FXC/D
          EXCD1 = FXCD1
          EXCD2 = FXCD2
        ELSE
          EXC = (FXC - FXCLDA)/D
          EXCD1 = FXCD1 - VXCLDA1
          EXCD2 = FXCD2 - VXCLDA2
        ENDIF
        EXCQ1 = SIGN(1.0_q,DD1)*FXCDD1
        EXCQ2 = SIGN(1.0_q,DD2)*FXCDD2
        ECQ = 0._q

        ! Hartree -> Rydberg conversion
        EXC    = 2._q* EXC
        EXCD1  = 2._q* EXCD1
        EXCD2  = 2._q* EXCD2
        EXCQ1  = 2._q* EXCQ1
        EXCQ2  = 2._q* EXCQ2
        ECQ    = 2._q* ECQ
!aem AM05 spin formulation added
! jP: adding PBEsol
!     Note that this is a plain copy-paste - but it offers
!     the possibility to add a screened PBEsol version later,
!     for a possible PBEsol hybrid-functional  - closely related to HSE
!     Further note: necessary condition: fitting the PBEsol-exchange hole!!!
      ELSE IF (LEXCH==14) THEN
           ukfactor=1.0_q


        D=2*D1

#undef exchange_RHOU_RHOD

#ifdef  exchange_RHOU_RHOD
        D=D1+D2
#endif
        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN

           S=DD1/(D*FK)
#ifdef  exchange_RHOU_RHOD
           S=(DD1+DD2)/2/(D*FK)
#endif

        ELSE
           S=0.0_q
        ENDIF
        IF (LDASCREEN==0._q .OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
           CALL EXCHPBESOL(D,DTHRD,S,EXLDA1,EXC1,EXDLDA1,EXCD1,EXCQ1, &
          &     ukfactor)
        ELSE
           CALL CALC_EXCHWPBEsol_SP(D,S, &
          &  EXLDA1,EXDLDA1,EXLDA_SR1,EXLDA_LR1,EXDLDA_SR1,EXDLDA_LR1, &
          &  EXWPBE_SR1,EXWPBE_LR1,EXWPBED_SR1,EXWPBED_LR1,EXWPBEDD_SR1,EXWPBEDD_LR1)
        ENDIF
!        WRITE(*,'(10F14.7)') D,S,(EXC1-EXLDA1)/D,EXCD1-EXDLDA1,EXCQ1
!        EXLDA1=0; EXC1=0; EXDLDA1=0; EXCD1=0; EXCQ1=0

        D=2*D2
#ifdef  exchange_RHOU_RHOD
        D=D1+D2
#endif

        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN
           S=DD2/(D*FK)
#ifdef  exchange_RHOU_RHOD
           S=(DD1+DD2)/2/(D*FK)
#endif
        ELSE
           S=0.0_q
        ENDIF

        IF (LDASCREEN==0._q.OR. LUSE_THOMAS_FERMI .OR. FORCE_PBE) THEN
           CALL EXCHPBESOL(D,DTHRD,S,EXLDA2,EXC2,EXDLDA2,EXCD2,EXCQ2, &
          &     ukfactor)
        ELSE
           CALL CALC_EXCHWPBEsol_SP(D,S, &
          &  EXLDA2,EXDLDA2,EXLDA_SR2,EXLDA_LR2,EXDLDA_SR2,EXDLDA_LR2, &
          &  EXWPBE_SR2,EXWPBE_LR2,EXWPBED_SR2,EXWPBED_LR2,EXWPBEDD_SR2,EXWPBEDD_LR2)
        ENDIF
!        EXLDA2=0; EXC2=0; EXDLDA2=0; EXCD2=0; EXCQ2=0
!        WRITE(*,'(10F14.7)') D,S,(EXC2-EXLDA2)/D,EXCD2-EXDLDA2,EXCQ2

        D=D1+D2
        DTHRD=exp(log(D)*THRD)
        RS=(0.75_q/PI)**THRD/DTHRD
        ZETA=(D1-D2)/D
        ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)

        FK=(3._q*PI*PI)**THRD*DTHRD
        SK = SQRT(4.0_q*FK/PI)

        G = (exp((2*THRD)*log(1._q+ZETA)) &
                +exp((2*THRD)*log(1._q-ZETA)))/2._q
        T = DDA/(D*2._q*SK*G)

        CALL corpbesol(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA,G,SK, &
                     T,EC,ECD1,ECD2,ECQ,.TRUE.)
!        ECLDA=0 ; ECD1LDA=0 ; ECD2LDA=0; EC=0 ; ECD1=0 ; ECD2=0 ; ECQ=0
!        WRITE(*,'(10F14.7)') D,RS,SK,T,EC,ECD1,ECD2,ECQ

        IF (LLDA) THEN
! Add LDA contributions as well
           IF (LDASCREEN==0._q.OR. LUSE_THOMAS_FERMI  .OR. FORCE_PBE) THEN
              EXC  =EC  *AGGAC +(EXC1-EXLDA1+EXC2-EXLDA2)/(2*D)*AGGAX+(EXLDA1+EXLDA2)/(2*D)*LDAX+ECLDA
              EXCD1=ECD1*AGGAC +(EXCD1-EXDLDA1)*AGGAX +EXDLDA1*LDAX+ECD1LDA
              EXCD2=ECD2*AGGAC +(EXCD2-EXDLDA2)*AGGAX +EXDLDA2*LDAX+ECD2LDA
              EXCQ1=EXCQ1 *AGGAX
              EXCQ2=EXCQ2 *AGGAX
              ECQ  =ECQ *AGGAC
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+(EXWPBE_LR1+EXWPBE_LR2)/(2*D)*AGGAX +(EXWPBE_SR1+EXWPBE_SR2)/(2*D)+ECLDA
                 EXCD1=ECD1*AGGAC+ EXWPBED_LR1*AGGAX +EXWPBED_SR1+ECD1LDA
                 EXCD2=ECD2*AGGAC+ EXWPBED_LR2*AGGAX +EXWPBED_SR2+ECD2LDA
                 EXCQ1=EXWPBEDD_LR1*AGGAX+EXWPBEDD_SR1
                 EXCQ2=EXWPBEDD_LR2*AGGAX+EXWPBEDD_SR2
                 ECQ  =ECQ*AGGAC
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC  =EC*AGGAC  +(EXWPBE_LR1+EXWPBE_LR2)/(2*D)*AGGAX+ECLDA
                 EXCD1=ECD1*AGGAC+ EXWPBED_LR1*AGGAX+ECD1LDA
                 EXCD2=ECD2*AGGAC+ EXWPBED_LR2*AGGAX+ECD2LDA
                 EXCQ1=EXWPBEDD_LR1*AGGAX
                 EXCQ2=EXWPBEDD_LR2*AGGAX
                 ECQ =ECQ*AGGAC
              ELSE
! wPBE
                 EXC  =EC  *AGGAC+(EXWPBE_SR1+EXWPBE_SR2)/(2*D)*AGGAX +(EXWPBE_LR1+EXWPBE_LR2)/(2*D)+ECLDA
                 EXCD1=ECD1*AGGAC+ EXWPBED_SR1*AGGAX +EXWPBED_LR1+ECD1LDA
                 EXCD2=ECD2*AGGAC+ EXWPBED_SR2*AGGAX +EXWPBED_LR2+ECD2LDA
                 EXCQ1=EXWPBEDD_SR1*AGGAX+EXWPBEDD_LR1
                 EXCQ2=EXWPBEDD_SR2*AGGAX+EXWPBEDD_LR2
                 ECQ  =ECQ*AGGAC
              ENDIF
           ENDIF
        ELSE
! Do not add LDA contribution
           IF (LDASCREEN==0._q.OR. LUSE_THOMAS_FERMI  .OR. FORCE_PBE) THEN
              EXC  =EC*  AGGAC +(EXC1-EXLDA1+EXC2-EXLDA2)/(2*D)*AGGAX
              EXCD1=ECD1*AGGAC +(EXCD1-EXDLDA1)*AGGAX
              EXCD2=ECD2*AGGAC +(EXCD2-EXDLDA2)*AGGAX
              EXCQ1=EXCQ1 *AGGAX
              EXCQ2=EXCQ2 *AGGAX
              ECQ  =ECQ *AGGAC
           ELSE
              IF (LUSE_LONGRANGE_HF) THEN
! Iann Gerber's functionals
                 EXC  =EC  *AGGAC+(EXWPBE_LR1-EXLDA_LR1+EXWPBE_LR2-EXLDA_LR2)/(2*D)*AGGAX+ &
                &                 (EXWPBE_SR1-EXLDA_SR1+EXWPBE_SR2-EXLDA_SR2)/(2*D)
                 EXCD1=ECD1*AGGAC+(EXWPBED_LR1-EXDLDA_LR1)*AGGAX+(EXWPBED_SR1-EXDLDA_SR1)
                 EXCD2=ECD2*AGGAC+(EXWPBED_LR2-EXDLDA_LR2)*AGGAX+(EXWPBED_SR2-EXDLDA_SR2)
                 EXCQ1=EXWPBEDD_LR1*AGGAX+EXWPBEDD_SR1
                 EXCQ2=EXWPBEDD_LR2*AGGAX+EXWPBEDD_SR2
                 ECQ  =ECQ *AGGAC
              ELSE IF (LUSE_MODEL_HF) THEN
                 EXC=EC*AGGAC+ (EXWPBE_LR1-EXLDA_LR1+EXWPBE_LR2-EXLDA_LR2)/(2*D)*AGGAX
                 EXCD1=ECD1*AGGAC+(EXWPBED_LR1-EXDLDA_LR1)*AGGAX
                 EXCD2=ECD2*AGGAC+(EXWPBED_LR2-EXDLDA_LR2)*AGGAX
                 EXCQ1=EXWPBEDD_LR1*AGGAX
                 EXCQ2=EXWPBEDD_LR2*AGGAX
                 ECQ =ECQ *AGGAC
              ELSE
! wPBE
                 EXC  =EC  *AGGAC+(EXWPBE_SR1-EXLDA_SR1+EXWPBE_SR2-EXLDA_SR2)/(2*D)*AGGAX+ &
                &                 (EXWPBE_LR1-EXLDA_LR1+EXWPBE_LR2-EXLDA_LR2)/(2*D)
                 EXCD1=ECD1*AGGAC+(EXWPBED_SR1-EXDLDA_SR1)*AGGAX+(EXWPBED_LR1-EXDLDA_LR1)
                 EXCD2=ECD2*AGGAC+(EXWPBED_SR2-EXDLDA_SR2)*AGGAX+(EXWPBED_LR2-EXDLDA_LR2)
                 EXCQ1=EXWPBEDD_SR1*AGGAX+EXWPBEDD_LR1
                 EXCQ2=EXWPBEDD_SR2*AGGAX+EXWPBEDD_LR2
                 ECQ  =ECQ *AGGAC
              ENDIF
           ENDIF
        ENDIF
! Hartree -> Rydberg conversion
        EXC    = 2*EXC
        EXCD1  = 2*EXCD1
        EXCD2  = 2*EXCD2
        EXCQ1  = 2*EXCQ1
        EXCQ2  = 2*EXCQ2
        ECQ    = 2*ECQ
! jH-new otherwise spin-polarized not possible
      ELSE IF (LEXCH==20 .OR. LEXCH==30) THEN
        EXC    = 0._q
        EXCD1  = 0._q
        EXCD2  = 0._q
        EXCQ1  = 0._q
        EXCQ2  = 0._q
        ECQ    = 0._q
! jP: adding PBEsol

!BEEF
#ifdef libbeef
      ELSE IF (LEXCH==17) THEN
        D = D1+D2
        IF (D>1E-10_q) THEN
           ZETA=(D1-D2)/D
           ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)
           IF (D1>1E-10_q) THEN
            DDSQ = 4._q*DD1*DD1
            CALL BEEFX(2._q*D1, DDSQ, EXBEEF1, EXBEEFD1, EXCQ1, LLDA)
            IF (DDSQ<1E-10_q) THEN
             EXBEEFD1 = 0._q
             EXCQ1 = 0._q
            END IF
           ELSE
            EXBEEF1 = 0._q
            EXBEEFD1 = 0._q
            EXCQ1 = 0._q
           END IF
           IF (D2>1E-10_q) THEN
            DDSQ = 4._q*DD2*DD2
            CALL BEEFX(2._q*D2, DDSQ, EXBEEF2, EXBEEFD2, EXCQ2, LLDA)
            IF (DDSQ<1E-10_q) THEN
             EXBEEFD2 = 0._q
             EXCQ2 = 0._q
            END IF
           ELSE
            EXBEEF2 = 0._q
            EXBEEFD2 = 0._q
            EXCQ2 = 0._q
           END IF
           CALL BEEFLOCALCORRSPIN(D, ZETA, DDA*DDA, ECBEEF, ECBEEFD1, ECBEEFD2, ECQ, LLDA)

           EXC = (EXBEEF1 + EXBEEF2 + 2._q*ECBEEF)/D
           EXCD1 = 2._q*(EXBEEFD1 + ECBEEFD1)
           EXCD2 = 2._q*(EXBEEFD2 + ECBEEFD2)
           EXCQ1 = 4._q*EXCQ1*DD1
           EXCQ2 = 4._q*EXCQ2*DD2
           ECQ = 2._q*ECQ*DDA
        ELSE
           EXC = 0._q
           EXCD1 = 0._q
           EXCD2 = 0._q
           EXCQ1 = 0._q
           EXCQ2 = 0._q
           ECQ = 0._q
        END IF
#endif

      ELSE IF (LEXCH>20.AND.LEXCH<25) THEN
! functionals for range-separated ACFDT (LDA - short range RPA)
! N.B. even if these functionals would support spin polarization,
! this would be a dummy stub since these functionals are not
! gradient corrected.
NOACC    CALL vtutor%error('GGA = 03 | 05 | 10 | 20  does not support spin polarization!')
#ifdef USELIBXC
      ELSE IF (LEXCH==99) THEN
         XC_DENS(1) = D1
         XC_DENS(2) = D2
         XC_SIGMA(1) = DD1*DD1
         XC_SIGMA(2) = 0.5_q*(DDA*DDA - DD1*DD1 - DD2*DD2)
         XC_SIGMA(3) = DD2*DD2
         XC_SIGMAZERO = 0._q
         EXC    = 0._q
         EXCD1  = 0._q
         EXCD2  = 0._q
         EXCQ1  = 0._q
         EXCQ2  = 0._q
         ECQ    = 0._q
         CX = -(3._q/4._q)*(3._q/PI)**THRD
         EXLDA = 0.5_q*CX*((2._q*D1)**THRD4 + (2._q*D2)**THRD4)/(D1+D2)
         EXDLDA1 = THRD4*CX*(2._q*D1)**THRD
         EXDLDA2 = THRD4*CX*(2._q*D2)**THRD

         DO IFUNC=1,2
            IF (LIBXC_MASK(IFUNC)>0) THEN
               IF ((XC_FAMILY(IFUNC) .EQ. LIBXC_FAMILY_LDA) .OR. (XC_FAMILY(IFUNC) .EQ. LIBXC_FAMILY_HYB_LDA)) THEN
                  IF (LIBXC_NOEXC(IFUNC)) THEN
                     CALL XC_F03_LDA_VXC(LIBXC_FUNC_SPN(IFUNC),1_8,XC_DENS,XC_VRHO)
                     XC_EXC = 0._q
                  ELSE
                     CALL XC_F03_LDA_EXC_VXC(LIBXC_FUNC_SPN(IFUNC),1_8,XC_DENS,XC_EXC,XC_VRHO)
                  ENDIF
                  IF ((XC_KIND(IFUNC) .EQ. XC_EXCHANGE) .AND. (LIBXC_HYB(IFUNC) .EQV. .FALSE.)) THEN
                     EXC    = EXC + LDAX*XC_EXC(1)
                     EXCD1  = EXCD1 + LDAX*XC_VRHO(1)
                     EXCD2  = EXCD2 + LDAX*XC_VRHO(2)
                  ELSEIF (XC_KIND(IFUNC) .EQ. XC_CORRELATION) THEN
                     EXC    = EXC + ALDAC*XC_EXC(1)
                     EXCD1  = EXCD1 + ALDAC*XC_VRHO(1)
                     EXCD2  = EXCD2 + ALDAC*XC_VRHO(2)
                  ELSEIF ((XC_KIND(IFUNC) .EQ. XC_EXCHANGE_CORRELATION) .OR. (LIBXC_HYB(IFUNC) .EQV. .TRUE.)) THEN
                     EXC    = EXC + XC_EXC(1)
                     EXCD1  = EXCD1 + XC_VRHO(1)
                     EXCD2  = EXCD2 + XC_VRHO(2)
                  ENDIF
               ELSEIF ((XC_FAMILY(IFUNC) .EQ. LIBXC_FAMILY_GGA) .OR. (XC_FAMILY(IFUNC) .EQ. LIBXC_FAMILY_HYB_GGA)) THEN
                  IF (LIBXC_NOEXC(IFUNC)) THEN
                     CALL XC_F03_GGA_VXC(LIBXC_FUNC_SPN(IFUNC),1_8,XC_DENS,XC_SIGMA,XC_VRHO,XC_VSIGMA)
                     XC_EXC = 0._q
                  ELSE
                     CALL XC_F03_GGA_EXC_VXC(LIBXC_FUNC_SPN(IFUNC),1_8,XC_DENS,XC_SIGMA,XC_EXC,XC_VRHO,XC_VSIGMA)
                  ENDIF
                  IF ((XC_KIND(IFUNC) .EQ. XC_EXCHANGE) .AND. (LIBXC_HYB(IFUNC) .EQV. .FALSE.)) THEN
                     EXC    = EXC + AGGAX*(XC_EXC(1)-EXLDA) + LDAX*EXLDA
                     EXCD1  = EXCD1 + AGGAX*(XC_VRHO(1)-EXDLDA1) + LDAX*EXDLDA1
                     EXCD2  = EXCD2 + AGGAX*(XC_VRHO(2)-EXDLDA2) + LDAX*EXDLDA2
                     EXCQ1  = EXCQ1 + AGGAX*(XC_VSIGMA(1)*2._q*DD1 - XC_VSIGMA(2)*DD1)
                     EXCQ2  = EXCQ2 + AGGAX*(XC_VSIGMA(3)*2._q*DD2 - XC_VSIGMA(2)*DD2)
                     ECQ    = ECQ + AGGAX*XC_VSIGMA(2)*DDA
                  ELSEIF (XC_KIND(IFUNC) .EQ. XC_CORRELATION) THEN
                     IF (LIBXC_NOEXC(IFUNC)) THEN
                        CALL XC_F03_GGA_VXC(LIBXC_FUNC_SPN(IFUNC),1_8,XC_DENS,XC_SIGMAZERO,XC_VRHOZERO,XC_VSIGMAZERO)
                        XC_EXCZERO = 0._q
                     ELSE
                        CALL XC_F03_GGA_EXC_VXC(LIBXC_FUNC_SPN(IFUNC),1_8,XC_DENS,XC_SIGMAZERO,XC_EXCZERO,XC_VRHOZERO,XC_VSIGMAZERO)
                     ENDIF
                     EXC    = EXC + AGGAC*(XC_EXC(1)-XC_EXCZERO(1)) + ALDAC*XC_EXCZERO(1)
                     EXCD1  = EXCD1 + AGGAC*(XC_VRHO(1)-XC_VRHOZERO(1)) + ALDAC*XC_VRHOZERO(1)
                     EXCD2  = EXCD2 + AGGAC*(XC_VRHO(2)-XC_VRHOZERO(2)) + ALDAC*XC_VRHOZERO(2)
                     EXCQ1  = EXCQ1 + AGGAC*(XC_VSIGMA(1)*2._q*DD1 - XC_VSIGMA(2)*DD1)
                     EXCQ2  = EXCQ2 + AGGAC*(XC_VSIGMA(3)*2._q*DD2 - XC_VSIGMA(2)*DD2)
                     ECQ    = ECQ + AGGAC*XC_VSIGMA(2)*DDA
                  ELSEIF ((XC_KIND(IFUNC) .EQ. XC_EXCHANGE_CORRELATION) .OR. (LIBXC_HYB(IFUNC) .EQV. .TRUE.)) THEN
                     EXC    = EXC + XC_EXC(1)
                     EXCD1  = EXCD1 + XC_VRHO(1)
                     EXCD2  = EXCD2 + XC_VRHO(2)
                     EXCQ1  = EXCQ1 + XC_VSIGMA(1)*2._q*DD1 - XC_VSIGMA(2)*DD1
                     EXCQ2  = EXCQ2 + XC_VSIGMA(3)*2._q*DD2 - XC_VSIGMA(2)*DD2
                     ECQ    = ECQ + XC_VSIGMA(2)*DDA
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         !Conversion to Rydbergs
         EXC   = 2._q*EXC
         EXCD1 = 2._q*EXCD1
         EXCD2 = 2._q*EXCD2
         EXCQ1 = 2._q*EXCQ1
         EXCQ2 = 2._q*EXCQ2
         ECQ   = 2._q*ECQ
#endif
      ELSE
NOACC    CALL vtutor%bug('internal ERROR GGASPIN: Wrong LEXCH, scheme not implemented!', __FILE__, __LINE__)
      ENDIF

      RETURN
      END

!************************ SUBROUTINE GGAXCS *****************************
!
!  calculates the PW 91 xc-functional for spin-dependent potentials
!  using the algorithm of J.A. White and D.M. Bird
!    (Phys.Rev.B 50,7 (1994) 4954)
!  presently finite differences are used to calculate
!  the required derivatives
!
!***********************************************************************

      SUBROUTINE GGAXCS(ro1,ro2,q1,q2,exc, &
     &   excd1,excd2,excq1,excq2)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      REAL(q),PARAMETER :: PI=3.1415926536_q
      REAL(q),PARAMETER :: THRD=0.333333333333333_q
!$ACC ROUTINE (FUN,DFRO1,DFRO2,DFQ1,DFQ2) SEQ

      ro   = ro1 + ro2
      exc  = fun(ro1,ro2,q1,q2)
!----------------------------------------------------------------------
!  The derivatives of f_xc necessary to
!  construct spin up and down potentials are
!  computed numerically
!----------------------------------------------------------------------
      excd1 = dfro1(ro1,ro2,q1,q2)
      excd2 = dfro2(ro1,ro2,q1,q2)
      excq1 = dfq1 (ro1,ro2,q1,q2)
      excq2 = dfq2 (ro1,ro2,q1,q2)

      RETURN
      END


!***********************************************************************
!     f_xc PW91 construction
!***********************************************************************
      function fun(ro1,ro2,q1,q2)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
!$ACC ROUTINE (EXPW,CPW) SEQ

      ro   = ro1 + ro2
      qq   = q1  + q2
      xi   = (ro1 - ro2)/ro

      call expw(ro1,q1,ex1)
      call expw(ro2,q2,ex2)
      call cpw(xi,ro,qq,ec)
      despin=0

      fx1  = ro1*ex1
      fx2  = ro2*ex2
      fcg  = ro*ec
      fsg  = ro*despin
      fun  = fx1 + fx2 + fcg + fsg

      return
      end

      function dfq1(ro1,ro2,q1,q2)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      parameter (delta=1E-2_q)
!$ACC ROUTINE (FUN) SEQ

      eps = min(delta,abs(0.001_q*q1))+1E-15_q
      x1  = q1 - eps
      x3  = q1 + eps
      f1  = fun(ro1,ro2,x1,q2)
      f3  = fun(ro1,ro2,x3,q2)
      dfq1 = (f3-f1) / (2*eps)
      return
      end

      function dfq2(ro1,ro2,q1,q2)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      parameter (delta=1E-2_q)
!$ACC ROUTINE (FUN) SEQ

      eps = min(delta,abs(0.001_q*q2))+1E-15_q
      x1  = q2 - eps
      x3  = q2 + eps
      f1  = fun(ro1,ro2,q1,x1)
      f3  = fun(ro1,ro2,q1,x3)
      dfq2 = (f3-f1) / (2*eps)
      return
      end

      function dfro1(ro1,ro2,q1,q2)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      parameter (delta=1E-4_q)
!$ACC ROUTINE (FUN) SEQ

      eps = min(delta,abs(0.001_q*ro1))+1E-15_q
      x1  = ro1 - eps
      x3  = ro1 + eps
      f1  = fun(x1,ro2,q1,q2)
      f3  = fun(x3,ro2,q1,q2)
      dfro1 = (f3-f1) / (2*eps)
      return
      end

      function dfro2(ro1,ro2,q1,q2)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      parameter (delta=1E-4_q)
!$ACC ROUTINE (FUN) SEQ

      eps = min(delta,abs(0.001_q*ro2))+1E-15_q
      x1  = ro2 - eps
      x3  = ro2 + eps
      f1  = fun(ro1,x1,q1,q2)
      f3  = fun(ro1,x3,q1,q2)
      dfro2 = (f3-f1) / (2*eps)
      return
      end

!***********************************************************************
!     GGA91 EXCHANGE
!  gga91 exchange for a spin-unpolarized electronic system
!  input d : density
!  input s:  abs(grad d)/(2*kf*d)
!  output:  exchange energy per electron (ex) and its derivatives
!           w.r.t. d (exd) and dd (exdd)
!***********************************************************************

      subroutine expw(rhos,rhops,ex)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)

      REAL(q),PARAMETER :: A1=0.19645_q,A2=0.27430_q,A3=0.15084_q,A4=100.0_q
      REAL(q),PARAMETER :: AX=-0.7385588_q,A=7.7956_q,B1=0.004_q
      REAL(q),PARAMETER :: THRD=0.33333333333333_q,THRD4=1.333333333333333_q
      REAL(q),PARAMETER :: THPITH=3.0936677262801_q

      rho=2._q*rhos
      rhop=2._q*rhops
      FAC = AX*rho**THRD

      w=0.16162046_q
      s=w*abs(rhop)/rho**(4._q/3._q)
      S2 = S*S
      S3 = S2*S
      S4 = S2*S2

      S2 = S*S
      S3 = S2*S
      S4 = S2*S2
      P0 = 1.0_q/SQRT(1.0_q+A*A*S2)
      P1 = LOG(A*S+1.0_q/P0)
      P2 = EXP(-A4*S2)
      P3 = 1.0_q/(1.0_q+A1*S*P1+B1*S4)
      P4 = 1.0_q+A1*S*P1+(A2-A3*P2)*S2
      F = P3*P4
      EX = FAC*(F-1)
      RETURN
      END

!***********************************************************************
!     GGA91 CORRELATION
!***********************************************************************
      subroutine cpw(xi,ro,qq,ecc)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
!$ACC ROUTINE(CORLSD) SEQ

      REAL(q),PARAMETER :: C1=0.002568_q,C2=0.023266_q,C3=7.389E-6_q,C4=8.723_q
      REAL(q),PARAMETER :: C5=0.472_q,C6=7.389E-2_q,A4=100.0_q
      REAL(q),PARAMETER :: XNU=15.75592_q,CC0=0.004235_q,CX=-0.001667212_q,ALF=0.09_q
      REAL(q),PARAMETER :: THRDM=-0.333333333333_q,THRD2=0.666666666667_q
      REAL(q),PARAMETER :: THRD=0.33333333333333_q,SIXTH7=1.16666666666667_q
      REAL(q),PARAMETER :: PI=3.1415926536_q
!
!     LOCAL CORRELATION is needed to compute EC:
!
      rs   = (4._q*pi*ro/3._q)**(-1._q/3._q)
      zet   = xi
      CALL CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
!
!     NONLOCAL CORRELATION:
!
      FK = 1.91915829_q/RS
      SK = SQRT(4.0_q*FK/PI)
      G = ((1.0_q+ZET)**THRD2+(1.0_q-ZET)**THRD2)/2.0_q
      t =  abs(qq)/(2.0_q*G*sk*ro)
!
      BET = XNU*CC0
      DELT = 2.0_q*ALF/BET
      G3 = G**3
      G4 = G3*G
      PON = -DELT*EC/(G3*BET)
      B = DELT/(EXP(PON)-1.0_q)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      T6 = T4*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.0_q+B*T2
      Q5 = 1.0_q+B*T2+B2*T4
      Q6 = C1+C2*RS+C3*RS2
      Q7 = 1.0_q+C4*RS+C5*RS2+C6*RS3
      CC = -CX + Q6/Q7
      R0 = (SK/FK)**2
      R1 = A4*R0*G4
      COEFF = CC-CC0-3.0_q*CX/7.0_q
      R2 = XNU*COEFF*G3
      R3 = EXP(-R1*T2)
      H0 = G3*(BET/DELT)*LOG(1.0_q+DELT*Q4*T2/Q5)
      H1 = R3*R2*T2
      H = H0+H1
!
      ECC=H
!
      RETURN
      END


      SUBROUTINE CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
!  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
!  INPUT: SEITZ RADIUS (RS), RELATIVE SPIN POLARIZATION (ZET)
!  OUTPUT: CORRELATION ENERGY PER ELECTRON (EC), UP- AND DOWN-SPIN
!     POTENTIALS (VCUP,VCDN), DERIVATIVES OF EC WRT RS (ECRS) & ZET (ECZET)
!  OUTPUT: CORRELATION CONTRIBUTION (ALFC) TO THE SPIN STIFFNESS
!     IMPLICIT REAL(q) (A-H,O-Z)
      REAL(q),PARAMETER :: GAM=0.5198421E0_q,FZZ=1.709921E0_q
      REAL(q),PARAMETER :: THRD=0.333333333333E0_q,THRD4=1.333333333333E0_q
!$ACC ROUTINE (GCOR) SEQ
      F = ((1.E0_q+ZET)**THRD4+(1.E0_q-ZET)**THRD4-2.E0_q)/GAM
      CALL GCOR(0.0310907E0_q,0.21370E0_q,7.5957E0_q,3.5876E0_q,1.6382E0_q, &
     &    0.49294E0_q,1.00E0_q,RS,EU,EURS)
      CALL GCOR(0.01554535E0_q,0.20548E0_q,14.1189E0_q,6.1977E0_q,3.3662E0_q, &
     &    0.62517E0_q,1.00E0_q,RS,EP,EPRS)
      CALL GCOR(0.0168869E0_q,0.11125E0_q,10.357E0_q,3.6231E0_q,0.88026E0_q, &
     &    0.49671E0_q,1.00E0_q,RS,ALFM,ALFRSM)
!  ALFM IS MINUS THE SPIN STIFFNESS ALFC
      ALFC = -ALFM
      Z4 = ZET**4
      EC = EU*(1.E0_q-F*Z4)+EP*F*Z4-ALFM*F*(1.E0_q-Z4)/FZZ
!  ENERGY DONE. NOW THE POTENTIAL:
      ECRS = EURS*(1.E0_q-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.E0_q-Z4)/FZZ
      FZ = THRD4*((1.E0_q+ZET)**THRD-(1.E0_q-ZET)**THRD)/GAM
      ECZET = 4.E0_q*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
     &        -(1.E0_q-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.E0_q-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      RETURN
      END

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE CORLSD_SIMD(RS, ZET, EC, VCUP, VCDN, ECRS, ECZET, ALFC, MASK0)
      USE prec
      USE simd
      IMPLICIT REAL(q) (A-H,O-Z)
!  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
!  INPUT: SEITZ RADIUS (RS), RELATIVE SPIN POLARIZATION (ZET)
!  OUTPUT: CORRELATION ENERGY PER ELECTRON (EC), UP- AND DOWN-SPIN
!     POTENTIALS (VCUP,VCDN), DERIVATIVES OF EC WRT RS (ECRS) & ZET (ECZET)
!  OUTPUT: CORRELATION CONTRIBUTION (ALFC) TO THE SPIN STIFFNESS
!     IMPLICIT REAL(q) (A-H,O-Z)
      TYPE(SIMD_REAL64), INTENT(IN)      :: RS, ZET
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EC, VCUP, VCDN, ECRS, ECZET, ALFC
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: EU, EURS, EP, EPRS, ALFM, ALFRSM
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: GAM = 0.5198421E0_q, GAM_INV = 1.0_q / GAM, FZZ = 1.709921E0_q, FZZ_INV = 1.0_q / FZZ
      REAL(q), PARAMETER                 :: THRD = 0.333333333333E0_q, THRD4 = 1.333333333333E0_q

      CALL GCOR_SIMD(0.0310907E0_q, 0.21370E0_q, 7.5957E0_q, 3.5876E0_q, 1.6382E0_q, 0.49294E0_q, 1.00E0_q, RS, EU, EURS, MASK0)

      CALL GCOR_SIMD(0.01554535E0_q, 0.20548E0_q, 14.1189E0_q, 6.1977E0_q, 3.3662E0_q, 0.62517E0_q, 1.00E0_q, RS, EP, EPRS, MASK0)

      CALL GCOR_SIMD(0.0168869E0_q, 0.11125E0_q, 10.357E0_q, 3.6231E0_q, 0.88026E0_q, 0.49671E0_q, 1.00E0_q, RS, ALFM, ALFRSM, MASK0)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, F, Z4, TEMP_3, TEMP_4, FZ, COMM)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         TEMP_1 = 1.0_q + ZET%X(II)
         TEMP_1 = TEMP_1 * TEMP_1
         TEMP_1 = TEMP_1 * TEMP_1
         TEMP_2 = 1.0_q - ZET%X(II)
         TEMP_2 = TEMP_2 * TEMP_2
         TEMP_2 = TEMP_2 * TEMP_2
         F = (TEMP_1 + TEMP_2 - 2.0_q) * GAM_INV
!  ALFM IS MINUS THE SPIN STIFFNESS ALFC
         TEMP_1 = ZET%X(II) * ZET%X(II)
         Z4 = TEMP_1 * TEMP_1
         TEMP_1 = F * Z4
         TEMP_2 = 1.0_q - TEMP_1
         TEMP_3 = (1.0_q - Z4) * FZZ_INV
         IF (MASK0%X(II)) THEN
            TEMP_4 = F * TEMP_3
            ALFC%X(II) = -ALFM%X(II)
            EC%X(II) = EU%X(II) * TEMP_2 + EP%X(II) * TEMP_1 - ALFM%X(II) * TEMP_4
!  ENERGY DONE. NOW THE POTENTIAL:
            ECRS%X(II) = EURS%X(II) * TEMP_2 + EPRS%X(II) * TEMP_1 - ALFRSM%X(II) * TEMP_4
         ENDIF
         FZ = THRD4 * (SIMD_POW(1.0_q + ZET%X(II), THRD) - SIMD_POW(1.0_q - ZET%X(II), THRD)) * GAM_INV
         IF (MASK0%X(II)) THEN
            ECZET%X(II) = 4.0_q * ZET%X(II) * ZET%X(II) * ZET%X(II) * F * (EP%X(II) - EU%X(II) + ALFM%X(II) * FZZ_INV) + &
     &           FZ * (Z4 * EP%X(II) - Z4 * EU%X(II) - ALFM%X(II) * TEMP_3)
            COMM = EC%X(II) - RS%X(II) * ECRS%X(II) * THRD - ZET%X(II) * ECZET%X(II)
            VCUP%X(II) = COMM + ECZET%X(II)
            VCDN%X(II) = COMM - ECZET%X(II)
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE CORLSD_SIMD
#endif


!=======================================================================
!
! original GCOR used in PW91 routines, to interpolate
! LSD correlation energy, as given by (10) of
! J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
! K. Burke, May 11, 1996.
!
!=======================================================================

      SUBROUTINE GCOR(A,A1,B1,B2,B3,B4,P,RS,GG,GGRS)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)

      P1 = P + 1.E0_q
      Q0 = -2.E0_q*A*(1.E0_q+A1*RS)
      RS12 = SQRT(RS)
      RS32 = RS12**3
      RSP = RS**P
      Q1 = 2.E0_q*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
      Q2 = LOG(1.E0_q+1.E0_q/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RS12+2.E0_q*B2+3.E0_q*B3*RS12+2.E0_q*B4*P1*RSP)
      GGRS = -2.E0_q*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
      RETURN
      END

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE GCOR_SIMD(A, A1, B1, B2, B3, B4, P, RS, GG, GGRS, MASK0)
      USE prec
      USE simd
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE(SIMD_REAL64), INTENT(IN)      :: RS
      TYPE(SIMD_REAL64), INTENT(OUT)     :: GG, GGRS
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      INTEGER                            :: II
      
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(P1, Q0, RS12, RS32, RSP, Q1, Q2, Q3)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         P1 = P + 1.0_q
         Q0 = -2.0_q * A * (1.0_q + A1 * RS%X(II))
         RS12 = SQRT(RS%X(II))
         RS32 = SIMD_POW(RS12, 3.0_q)
         RSP = SIMD_POW(RS%X(II), P)
         Q1 = 2.0_q * A * (B1 * RS12 + B3 * RS32 + (B2 + B4 * RSP) * RS%X(II))
         Q2 = SIMD_LOG(1.0_q + 1.0_q / Q1)
         Q3 = A * (B1 / RS12 + 2.0_q * B2 + 3.0_q * B3 * RS12 + 2.0_q * B4 * P1 * RSP)
         IF (MASK0%X(II)) THEN
            GG%X(II) = Q0 * Q2
            GGRS%X(II) = -2.0_q * A * A1 * Q2 - Q0 * Q3 / (Q1 * Q1 + Q1)
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE GCOR_SIMD
#endif


!***********************************************************************
!
! B3LYP exchange and correlation energy
! Note: for the Becke part only the gradient corrections are calculated
! whereas the correlation part includes the local (Wigner) term
!
! this subroutine calculates the derivative w.r.t. each argument
! of the B3LYP functional using finite differences
!
!***********************************************************************
     SUBROUTINE B3LYPXCS(RHO1,RHO2,DRHO1,DRHO2,DRHOA,EXC,EXCD1,EXCDD1, &
                         EXCD2,EXCDD2,EXCDDA)
!$ACC ROUTINE SEQ
     USE prec
     USE constant
     IMPLICIT NONE
     REAL(q) RHO1,RHO2,DRHO1,DRHO2,DRHOA
     REAL(q) EXC,EXCD1,EXCDD1,EXCD2,EXCDD2,EXCDDA
     REAL(q) EXCB3LYP

     REAL(q) DB3LYP_RHO1,DB3LYP_RHO2,DB3LYP_DRHO1,DB3LYP_DRHO2,DB3LYP_DRHOA
!$ACC ROUTINE (EXCB3LYP,DB3LYP_RHO1,DB3LYP_RHO2,DB3LYP_DRHO1,DB3LYP_DRHO2,DB3LYP_DRHOA) SEQ
!    B3LYP - Recipe:
!    EXC = EXC(LSDA) + a0 * (EX(HF) - EX(LSDA)) + aX * DELTA EX(B88) + aC * DELTA EC(LYP)
!
!    a0 = 0.20 ( == AEXX !)  : 0.2  * EX(LSDA)
!    aX = 0.72               : 0.72 * DELTA EX(B88)
!    aC = 0.81               : 0.81 * DELTA EC(LYP)
!

     IF (RHO1<0._q .OR. RHO2<0._q ) THEN
       EXC   = 0._q
       EXCD1  = 0._q
       EXCDD1 = 0._q
       EXCD2  = 0._q
       EXCDD2 = 0._q
       EXCDDA = 0._q
       RETURN
     ENDIF

! Exchange-correlation energy

     EXC=EXCB3LYP(RHO1,RHO2,DRHO1,DRHO2,DRHOA)

! Exchange correlation potential

     EXCD1=DB3LYP_RHO1(RHO1,RHO2,DRHO1,DRHO2,DRHOA)
     EXCD2=DB3LYP_RHO2(RHO1,RHO2,DRHO1,DRHO2,DRHOA)
     EXCDD1=DB3LYP_DRHO1(RHO1,RHO2,DRHO1,DRHO2,DRHOA)
     EXCDD2=DB3LYP_DRHO2(RHO1,RHO2,DRHO1,DRHO2,DRHOA)
     EXCDDA=DB3LYP_DRHOA(RHO1,RHO2,DRHO1,DRHO2,DRHOA)

!     WRITE(91,'(6F12.7)') RHO1, RHO2, EXCD1, EXCD2, EXCDD1, EXCDD2

     END SUBROUTINE B3LYPXCS

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE B3LYPXCS_SIMD(RHO1_IN, RHO2_IN, DRHO1, DRHO2, DRHOA, EXC, EXCD1, EXCDD1, EXCD2, EXCDD2, EXCDDA, MASK0)
      USE prec
      USE constant
      USE simd
      IMPLICIT NONE
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHO1_IN, RHO2_IN, DRHO1, DRHO2, DRHOA
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EXC, EXCD1, EXCDD1, EXCD2, EXCDD2, EXCDDA
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: RHO1, RHO2
      TYPE(SIMD_MASK_REAL64)             :: MASK1
      INTEGER                            :: II
      LOGICAL                            :: TRUE_FOR_ANY
!    B3LYP - Recipe:
!    EXC = EXC(LSDA) + a0 * (EX(HF) - EX(LSDA)) + aX * DELTA EX(B88) + aC * DELTA EC(LYP)
!
!    a0 = 0.20 ( == AEXX !)  : 0.2  * EX(LSDA)
!    aX = 0.72               : 0.72 * DELTA EX(B88)
!    aC = 0.81               : 0.81 * DELTA EC(LYP)
!
      TRUE_FOR_ANY = .FALSE.
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         MASK1%X(II) = MASK0%X(II)
         IF (MASK1%X(II)) THEN
            IF ((RHO1_IN%X(II) .GE. 0.0_q) .AND. (RHO2_IN%X(II) .GE. 0.0_q)) THEN
               TRUE_FOR_ANY = .TRUE.
               RHO1%X(II) = RHO1_IN%X(II)
               RHO2%X(II) = RHO2_IN%X(II)
            ELSE
               MASK1%X(II) = .FALSE.
               RHO1%X(II) = 0.0_q
               RHO2%X(II) = 0.0_q
               EXC%X(II) = 0.0_q
               EXCD1%X(II) = 0.0_q
               EXCDD1%X(II) = 0.0_q
               EXCD2%X(II) = 0.0_q
               EXCDD2%X(II) = 0.0_q
               EXCDDA%X(II) = 0.0_q
            ENDIF
         ENDIF
      ENDDO

      IF (.NOT. TRUE_FOR_ANY) RETURN

      CALL EXCB3LYP_SIMD(EXC, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK1)

      CALL DB3LYP_RHO1_SIMD(EXCD1, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK1)
      
      CALL DB3LYP_RHO2_SIMD(EXCD2, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK1)

      CALL DB3LYP_DRHO1_SIMD(EXCDD1, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK1)

      CALL DB3LYP_DRHO2_SIMD(EXCDD2, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK1)
 
      CALL DB3LYP_DRHOA_SIMD(EXCDDA, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK1)

      RETURN
      END SUBROUTINE B3LYPXCS_SIMD
#endif

!***********************************************************************
!
! B3LYP exchange and correlation energy
! Note: for the Becke part only the gradient corrections are calculated
! whereas the correlation part includes the local (Wigner) term
!
!***********************************************************************

      FUNCTION EXCB3LYP(RHO1,RHO2,DRHO1,DRHO2,DRHOA)
!$ACC ROUTINE SEQ
      USE prec
      USE setexm_struct_def

      IMPLICIT NONE
      REAL(q) EXCB3LYP

      REAL(q) :: RHO, RHO1,RHO2,DRHO1,DRHO2, DRHOA
      REAL(q) :: EXC, TEMP

      REAL(q) :: AX, AC
!$ACC ROUTINE (EXB88,ECLYP) SEQ

      AX=AGGAX
      AC=AGGAC

      EXCB3LYP = 0.0_q
      EXC=0._q
      RHO  = RHO1+RHO2

! Becke88 gradient corrections up electrons
      CALL EXB88(RHO1,DRHO1,TEMP)
      EXC = EXC + AX *TEMP
! Becke88 gradient corrections down electrons
      CALL EXB88(RHO2,DRHO2,TEMP)
      EXC = EXC + AX *TEMP

! LYP correlation
      CALL ECLYP(RHO1,RHO2,DRHO1,DRHO2,DRHOA,TEMP)
      EXC = EXC + AC*TEMP
      EXCB3LYP=EXC

      END FUNCTION EXCB3LYP

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE EXCB3LYP_SIMD(RES, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK0)
      USE prec
      USE setexm_struct_def
      USE simd
      IMPLICIT NONE
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHO1, RHO2, DRHO1, DRHO2, DRHOA
      TYPE(SIMD_REAL64), INTENT(OUT)     :: RES
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: TEMP, TEMP_RES
      INTEGER                            :: II

      CALL EXB88_SIMD(RHO1, DRHO1, TEMP, MASK0)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) 
      DO II = 0, SIMD_WIDTH_64BIT - 1
         TEMP_RES%X(II) = AGGAX * TEMP%X(II)
      ENDDO

      CALL EXB88_SIMD(RHO2, DRHO2, TEMP, MASK0)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) 
      DO II = 0, SIMD_WIDTH_64BIT - 1
         TEMP_RES%X(II) = TEMP_RES%X(II) +  AGGAX * TEMP%X(II)
      ENDDO

      CALL ECLYP_SIMD(RHO1, RHO2, DRHO1, DRHO2, DRHOA, TEMP, MASK0)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) 
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK0%X(II)) THEN
            RES%X(II) = TEMP_RES%X(II) + AGGAC * TEMP%X(II)
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE EXCB3LYP_SIMD
#endif

!**********************************************************************
!
! Becke 88 gradient correction to LDA exchange for a spin-unpolarized
! electronic system
! A.D. Becke, Phys. Rev. A, 38, 3098 (1988).
!
!**********************************************************************

      SUBROUTINE EXB88(RHOS,DRHOS,EX)
!$ACC ROUTINE SEQ
      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)

      REAL(q), PARAMETER :: BETA = 0.0042_q

      RHO43 = RHOS**(4.0_q/3.0_q)
      X = ABS(DRHOS)/RHO43
      X2 = X*X
      B6 = 6.0_q * BETA

      ! arsinh(x) = LOG(  x + SQRT( x * x + 1.0_q))
      DENOM = 1.0_q + B6*X*LOG(X + SQRT(X2 + 1.0_q))
      IF (RHOS .GE. 0.00001_q ) THEN
      EX = -BETA * RHO43 * X2/DENOM
      ELSE
      EX = 0.0_q
      END IF
      RETURN
      END SUBROUTINE EXB88

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE EXB88_SIMD(RHOS, DRHOS, EX, MASK0)
      USE prec
      USE constant
      USE simd
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHOS, DRHOS
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EX
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: BETA = -0.0042_q
      REAL(q), PARAMETER                 :: THRD4 = -4.0_q / 3.0_q 

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, TEMP_3, RHO43, X, X2, B6, DENOM)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         TEMP_1 = RHOS%X(II)
         TEMP_2 = DRHOS%X(II)
         TEMP_3 = EX%X(II)
         TEMP_3 = 0.0_q
         IF (TEMP_1 .GE. 0.00001_q) THEN
            RHO43 = SIMD_POW(TEMP_1, THRD4)
            X = ABS(TEMP_2) * RHO43
            X2 = X * X
            B6 = 6.0_q * BETA
            ! arsinh(x) = LOG(x + SQRT(x * x + 1.0_q))
            TEMP_2 = X + SQRT(1.0_q + X2)
            DENOM = 1.0_q - B6 * X * SIMD_LOG(TEMP_2)
            TEMP_3 = BETA * X2 / (RHO43 * DENOM)
         ENDIF
         IF (MASK0%X(II)) THEN
            EX%X(II) = TEMP_3
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE EXB88_SIMD
#endif

!**********************************************************************
!
! The correlation energy density of Lee, Yang and Parr
! B. Miehlich, A. Savin, H. Stoll and H. Preuss, Chem. Phys. Lett. 157,
! 200 (1989).
! C. Lee, W. Yang, and R. Parr, Phys. Rev. B 37, 785 (1988)
!
!**********************************************************************

      SUBROUTINE ECLYP(RHOU,RHOD,DRHOU,DRHOD,DRHOA,ELYP)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)

      REAL(q), PARAMETER :: CF=2.87123400018819152526_q
      REAL(q), PARAMETER :: ALYP=0.04918_q
      REAL(q), PARAMETER :: BLYP=0.132_q
      REAL(q), PARAMETER :: CLYP=0.2533_q
      REAL(q), PARAMETER :: DLYP=0.349_q

      RHO=RHOU+RHOD
      DRHO=DRHOA
      DRHOSQ=ABS(DRHO)*ABS(DRHO)
      DRHOUSQ=ABS(DRHOU)*ABS(DRHOU)
      DRHODSQ=ABS(DRHOD)*ABS(DRHOD)

      X =  CLYP* ( RHO**(-1.0_q/3.0_q) )

      Y =  DLYP* ( RHO**(-1.0_q/3.0_q) )

      Z =  1.0_q + DLYP * ( RHO**(-1.0_q/3.0_q) )

      OMEGA = EXP(-X)/Z * ( RHO**(-11.0_q/3.0_q) )

      DELTA = X + Y/Z

      SUMMAND_A=0
      SUMMAND_B=0
      SUMMAND_C=0
      SUMMAND_D=0
      SUMMAND_E=0

      IF (RHOU .GE. 0.00001_q .AND. RHOD .GE. 0.00001_q ) THEN
      ! Wigner term
      SUMMAND_A = -ALYP* 4._q/( 1._q + DLYP*RHO**(-1._q/3._q) ) * RHOU*RHOD/RHO
      ! second local term
      SUMMAND_B = 2._q**(11._q/3._q)*CF*( RHOU**(8._q/3._q) + RHOD**(8._q/3._q))
      ! non local terms
      SUMMAND_C = ( 47._q/18._q - 7._q/18._q * DELTA) * DRHOSQ - (5._q/2._q - 1._q/18._q * DELTA) * ( DRHOUSQ + DRHODSQ)

      SUMMAND_D = -(DELTA - 11._q)/9._q * ( (RHOU/RHO) * DRHOUSQ + (RHOD/RHO) * DRHODSQ)

      SUMMAND_E = (-2._q/3._q)*RHO*RHO*DRHOSQ + ( (2._q/3._q) *RHO*RHO - RHOU*RHOU) * DRHODSQ + &
                  ((2._q/3._q) *RHO*RHO - RHOD*RHOD) * DRHOUSQ

      ! "local" contribution (i.e Wigner term + term "B")
      ELYP =  SUMMAND_A - ALYP*BLYP*OMEGA* ( RHOU*RHOD * ( SUMMAND_B ) )

      ! "non-local" contribution
      ELYP = ELYP - ALYP*BLYP*OMEGA* ( RHOU*RHOD * ( SUMMAND_C + SUMMAND_D ) +SUMMAND_E )
      ELSE
      ELYP = 0.0_q
      END IF
      RETURN
      END SUBROUTINE ECLYP

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE ECLYP_SIMD(RHOU, RHOD, DRHOU, DRHOD, DRHOA, ELYP, MASK0)
      USE prec
      USE simd
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHOU, RHOD, DRHOU, DRHOD, DRHOA
      TYPE(SIMD_REAL64), INTENT(OUT)     :: ELYP
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: OMEGA, DELTA
      TYPE(SIMD_MASK_REAL64)             :: MASK1
      INTEGER                            :: II
      LOGICAL                            :: TRUE_FOR_ANY
      REAL(q), PARAMETER                 :: CF = 2.87123400018819152526_q
      REAL(q), PARAMETER                 :: ALYP = 0.04918_q
      REAL(q), PARAMETER                 :: BLYP = 0.132_q
      REAL(q), PARAMETER                 :: CLYP = 0.2533_q
      REAL(q), PARAMETER                 :: DLYP = 0.349_q
      REAL(q), PARAMETER                 :: THRD = 1.0_q / 3.0_q
      REAL(q), PARAMETER                 :: MTHRD = -1.0_q / 3.0_q
      REAL(q), PARAMETER                 :: P1 = -11.0_q / 3.0_q
      REAL(q), PARAMETER                 :: P2 = 2.0_q ** (11.0_q/3.0_q)
      REAL(q), PARAMETER                 :: P3 = 8.0_q / 3.0_q
      REAL(q), PARAMETER                 :: P4 = 47.0_q / 18.0_q
      REAL(q), PARAMETER                 :: P5 = 7.0_q / 18.0_q
      REAL(q), PARAMETER                 :: P6 = 1.0_q / 18.0_q
      REAL(q), PARAMETER                 :: P7 = 2.0_q / 3.0_q
      REAL(q), PARAMETER                 :: P8 = 1.0_q / 9.0_q

      TRUE_FOR_ANY = .FALSE.
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY) PRIVATE(TEMP_1, TEMP_2, RHO, RHO_MTHRD, X, Y, Z, RHO_INV, SUMMAND_A, SUMMAND_B, TEMP_3)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         TEMP_1 = RHOU%X(II)
         TEMP_2 = RHOD%X(II)
         RHO = TEMP_1 + TEMP_2
         RHO_MTHRD = SIMD_POW(RHO, MTHRD)
         X = CLYP * RHO_MTHRD
         Y = DLYP * RHO_MTHRD
         Z = 1.0_q / (1.0_q + DLYP * RHO_MTHRD)
         OMEGA%X(II) = SIMD_EXP(-X) * Z * SIMD_POW(RHO, P1)
         DELTA%X(II) = X + Y * Z
         MASK1%X(II) = MASK0%X(II)
         IF (MASK1%X(II)) THEN
            IF ((TEMP_1 .GE. 0.00001_q) .AND. (TEMP_2 .GE. 0.00001_q)) THEN
               TRUE_FOR_ANY = .TRUE.
               RHO_INV = 1.0_q / RHO
               ! Wigner term
               SUMMAND_A = -ALYP * 4.0_q / (1.0_q + DLYP * RHO_MTHRD) * TEMP_1 * TEMP_2 * RHO_INV
               ! second local term
               SUMMAND_B = P2 * CF * (SIMD_POW(TEMP_1, P3) + SIMD_POW(TEMP_2, P3))
               ! "local" contribution (i.e Wigner term + term "B")
               TEMP_3 = SUMMAND_A - ALYP * BLYP * OMEGA%X(II) * (TEMP_1 * TEMP_2 * SUMMAND_B)
            ELSE
               MASK1%X(II) = .FALSE.
               TEMP_3 = 0.0_q
            ENDIF
            ELYP%X(II) = TEMP_3
         ENDIF
      ENDDO

      IF (.NOT. TRUE_FOR_ANY) RETURN

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, RHO, DRHO, DRHOSQ, DRHOUSQ, DRHODSQ, RHO_INV, SUMMAND_C, SUMMAND_D, SUMMAND_E)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         TEMP_1 = RHOU%X(II)
         TEMP_2 = RHOD%X(II)
         RHO = TEMP_1 + TEMP_2
         DRHO = DRHOA%X(II)
         DRHOSQ = DRHO * DRHO
         DRHOUSQ = DRHOU%X(II) * DRHOU%X(II)
         DRHODSQ = DRHOD%X(II) * DRHOD%X(II)
         IF (MASK1%X(II)) THEN
            RHO_INV = 1.0_q / RHO
            ! non local terms
            SUMMAND_C = (P4 - P5 * DELTA%X(II)) * DRHOSQ - (2.5_q - P6 * DELTA%X(II)) * (DRHOUSQ + DRHODSQ)
            SUMMAND_D = -(DELTA%X(II) - 11.0_q) * P8 * (TEMP_1 * RHO_INV * DRHOUSQ + TEMP_2 * RHO_INV * DRHODSQ)              
            DELTA%X(II) = P7 * RHO * RHO
            SUMMAND_E = -DELTA%X(II) * DRHOSQ + (DELTA%X(II) - TEMP_1 * TEMP_1) * DRHODSQ + &
                        (DELTA%X(II) - TEMP_2 * TEMP_2) * DRHOUSQ
            ! "non-local" contribution
            ELYP%X(II) = ELYP%X(II) - ALYP * BLYP * OMEGA%X(II) * (TEMP_1 * TEMP_2 * (SUMMAND_C + SUMMAND_D) + SUMMAND_E)
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE ECLYP_SIMD
#endif

!**********************************************************************
!
! the following functions calculate the first derivative of
! the semi-local part of the B3LYP functional with respect to
! each argument
!
!**********************************************************************
      FUNCTION DB3LYP_DRHO1(RHO1,RHO2,DRHO1,DRHO2,DRHOA)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT NONE

      REAL(q), PARAMETER :: DELTA=1E-2_q
      REAL(q)            :: EPS, RHO1, RHO2, DRHO1, DRHO2,DRHOA, X1, X3, F1, F3
      REAL(q)            :: DB3LYP_DRHO1,EXCB3LYP
!$ACC ROUTINE(EXCB3LYP) SEQ

      EPS = MIN(DELTA,ABS(0.02_q*DRHO1))+1E-15_q
      X1  = DRHO1 - EPS
      X3  = DRHO1 + EPS
      F1  = EXCB3LYP(RHO1,RHO2,X1,DRHO2,DRHOA)
      F3  = EXCB3LYP(RHO1,RHO2,X3,DRHO2,DRHOA)
      DB3LYP_DRHO1 = (F3-F1) / (2*EPS)
      RETURN
      END FUNCTION DB3LYP_DRHO1

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE DB3LYP_DRHO1_SIMD(RES, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK0)
      USE prec
      USE simd
      IMPLICIT NONE
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHO1, RHO2, DRHO1, DRHO2, DRHOA
      TYPE(SIMD_REAL64), INTENT(OUT)     :: RES
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: EPS, X1, X3, F1, F3
      REAL(q)                            :: TEMP
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: DELTA = 1.0E-2_q

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         EPS%X(II) = MIN(DELTA, ABS(0.02_q * DRHO1%X(II))) + 1.0E-15_q
         X1%X(II) = DRHO1%X(II) - EPS%X(II)
         X3%X(II) = DRHO1%X(II) + EPS%X(II)
      ENDDO

      CALL EXCB3LYP_SIMD(F1, RHO1, RHO2, X1, DRHO2, DRHOA, MASK0)

      CALL EXCB3LYP_SIMD(F3, RHO1, RHO2, X3, DRHO2, DRHOA, MASK0)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK0%X(II)) THEN
            RES%X(II) = (F3%X(II) - F1%X(II)) / (2.0_q * EPS%X(II))
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DB3LYP_DRHO1_SIMD
#endif
!----------------------------------------------------------------------

      FUNCTION DB3LYP_DRHO2(RHO1,RHO2,DRHO1,DRHO2,DRHOA)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT NONE

      REAL(q), PARAMETER :: DELTA=1E-2_q
      REAL(q)            :: EPS, RHO1, RHO2, DRHO1, DRHO2,DRHOA, X1, X3, F1, F3
      REAL(q)            :: DB3LYP_DRHO2, EXCB3LYP
!$ACC ROUTINE(EXCB3LYP) SEQ

      EPS = MIN(DELTA,ABS(0.02_q*DRHO2))+1E-15_q
      X1  = DRHO2 - EPS
      X3  = DRHO2 + EPS
      F1  = EXCB3LYP(RHO1,RHO2,DRHO1,X1,DRHOA)
      F3  = EXCB3LYP(RHO1,RHO2,DRHO1,X3,DRHOA)
      DB3LYP_DRHO2 = (F3-F1) / (2*EPS)
      RETURN
      END FUNCTION DB3LYP_DRHO2

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE DB3LYP_DRHO2_SIMD(RES, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK0)
      USE prec
      USE simd
      IMPLICIT NONE
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHO1, RHO2, DRHO1, DRHO2, DRHOA
      TYPE(SIMD_REAL64), INTENT(OUT)     :: RES
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: EPS, X1, X3, F1, F3
      REAL(q)                            :: TEMP
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: DELTA = 1.0E-2_q

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         EPS%X(II) = MIN(DELTA, ABS(0.02_q * DRHO2%X(II))) + 1.0E-15_q
         X1%X(II) = DRHO2%X(II) - EPS%X(II)
         X3%X(II) = DRHO2%X(II) + EPS%X(II)
      ENDDO

      CALL EXCB3LYP_SIMD(F1, RHO1, RHO2, DRHO1, X1, DRHOA, MASK0)

      CALL EXCB3LYP_SIMD(F3, RHO1, RHO2, DRHO1, X3, DRHOA, MASK0)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK0%X(II)) THEN
            RES%X(II) = (F3%X(II) - F1%X(II)) / (2.0_q * EPS%X(II))
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DB3LYP_DRHO2_SIMD
#endif
!----------------------------------------------------------------------

      FUNCTION DB3LYP_DRHOA(RHO1,RHO2,DRHO1,DRHO2,DRHOA)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT NONE

      REAL(q), PARAMETER :: DELTA=1E-2_q
      REAL(q)            :: EPS, RHO1, RHO2, DRHO1, DRHO2,DRHOA, X1, X3, F1, F3
      REAL(q)            :: DB3LYP_DRHOA, EXCB3LYP
!$ACC ROUTINE(EXCB3LYP) SEQ

      EPS = MIN(DELTA,ABS(0.02_q*DRHOA))+1E-15_q
      X1  = DRHOA - EPS
      X3  = DRHOA + EPS
      F1  = EXCB3LYP(RHO1,RHO2,DRHO1,DRHO2,X1)
      F3  = EXCB3LYP(RHO1,RHO2,DRHO1,DRHO2,X3)
      DB3LYP_DRHOA = (F3-F1) / (2*EPS)
      RETURN
      END FUNCTION DB3LYP_DRHOA

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE DB3LYP_DRHOA_SIMD(RES, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK0)
      USE prec
      USE simd
      IMPLICIT NONE
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHO1, RHO2, DRHO1, DRHO2, DRHOA
      TYPE(SIMD_REAL64), INTENT(OUT)     :: RES
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: EPS, X1, X3, F1, F3
      REAL(q)                            :: TEMP
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: DELTA = 1.0E-2_q

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         EPS%X(II) = MIN(DELTA, ABS(0.02_q * DRHOA%X(II))) + 1.0E-15_q
         X1%X(II) = DRHOA%X(II) - EPS%X(II)
         X3%X(II) = DRHOA%X(II) + EPS%X(II)
      ENDDO

      CALL EXCB3LYP_SIMD(F1, RHO1, RHO2, DRHO1, DRHO2, X1, MASK0)

      CALL EXCB3LYP_SIMD(F3, RHO1, RHO2, DRHO1, DRHO2, X3, MASK0)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK0%X(II)) THEN
            RES%X(II) = (F3%X(II) - F1%X(II)) / (2.0_q * EPS%X(II))
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DB3LYP_DRHOA_SIMD
#endif
!----------------------------------------------------------------------

      FUNCTION DB3LYP_RHO1(RHO1,RHO2,DRHO1,DRHO2,DRHOA)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT NONE

      REAL(q), PARAMETER :: DELTA=1E-2_q
      REAL(q)            :: EPS, RHO1, RHO2, DRHO1, DRHO2,DRHOA, X1, X3, F1, F3
      REAL(q)            :: DB3LYP_RHO1, EXCB3LYP
!$ACC ROUTINE(EXCB3LYP) SEQ

      EPS = MIN(DELTA,ABS(0.02_q*RHO1))+1E-15_q
      X1  = RHO1 - EPS
      X3  = RHO1 + EPS
      F1  = EXCB3LYP(X1,RHO2,DRHO1,DRHO2,DRHOA)
      F3  = EXCB3LYP(X3,RHO2,DRHO1,DRHO2,DRHOA)
      DB3LYP_RHO1 = (F3-F1) / (2*EPS)
!      WRITE(82,*) RHO1, RHO2, DB3LYP_RHO1
      RETURN
      END FUNCTION DB3LYP_RHO1

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE DB3LYP_RHO1_SIMD(RES, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK0)
      USE prec
      USE simd
      IMPLICIT NONE
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHO1, RHO2, DRHO1, DRHO2, DRHOA
      TYPE(SIMD_REAL64), INTENT(OUT)     :: RES
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: EPS, X1, X3, F1, F3
      REAL(q)                            :: TEMP
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: DELTA = 1.0E-2_q

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         EPS%X(II) = MIN(DELTA, ABS(0.02_q * RHO1%X(II))) + 1.0E-15_q
         X1%X(II) = RHO1%X(II) - EPS%X(II)
         X3%X(II) = RHO1%X(II) + EPS%X(II)
      ENDDO

      CALL EXCB3LYP_SIMD(F1, X1, RHO2, DRHO1, DRHO2, DRHOA, MASK0)

      CALL EXCB3LYP_SIMD(F3, X3, RHO2, DRHO1, DRHO2, DRHOA, MASK0)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK0%X(II)) THEN
            RES%X(II) = (F3%X(II) - F1%X(II)) / (2.0_q * EPS%X(II))
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DB3LYP_RHO1_SIMD
#endif
!----------------------------------------------------------------------

      FUNCTION DB3LYP_RHO2(RHO1,RHO2,DRHO1,DRHO2,DRHOA)
!$ACC ROUTINE SEQ
      USE prec
      IMPLICIT NONE

      REAL(q), PARAMETER :: DELTA=1E-2_q
      REAL(q)            :: EPS, RHO1, RHO2, DRHO1, DRHO2,DRHOA, X1, X3, F1, F3
      REAL(q)            :: DB3LYP_RHO2, EXCB3LYP
!$ACC ROUTINE(EXCB3LYP) SEQ

      EPS = MIN(DELTA,ABS(0.02_q*RHO2))+1E-15_q
      X1  = RHO2 - EPS
      X3  = RHO2 + EPS
      F1  = EXCB3LYP(RHO1,X1,DRHO1,DRHO2,DRHOA)
      F3  = EXCB3LYP(RHO1,X3,DRHO1,DRHO2,DRHOA)
      DB3LYP_RHO2 = (F3-F1) / (2*EPS)
      RETURN
      END FUNCTION DB3LYP_RHO2

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE DB3LYP_RHO2_SIMD(RES, RHO1, RHO2, DRHO1, DRHO2, DRHOA, MASK0)
      USE prec
      USE simd
      IMPLICIT NONE
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHO1, RHO2, DRHO1, DRHO2, DRHOA
      TYPE(SIMD_REAL64), INTENT(OUT)     :: RES
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: EPS, X1, X3, F1, F3
      REAL(q)                            :: TEMP
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: DELTA = 1.0E-2_q

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         EPS%X(II) = MIN(DELTA, ABS(0.02_q * RHO2%X(II))) + 1.0E-15_q
         X1%X(II) = RHO2%X(II) - EPS%X(II)
         X3%X(II) = RHO2%X(II) + EPS%X(II)
      ENDDO
      
      CALL EXCB3LYP_SIMD(F1, RHO1, X1, DRHO1, DRHO2, DRHOA, MASK0)

      CALL EXCB3LYP_SIMD(F3, RHO1, X3, DRHO1, DRHO2, DRHOA, MASK0)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK0%X(II)) THEN
           RES%X(II) = (F3%X(II) - F1%X(II)) / (2.0_q * EPS%X(II))
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE DB3LYP_RHO2_SIMD
#endif

!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
!=====================================================================
!Functional wPBE short-range for calculation in separated interactions LR-SR
!with modified kernel Short-Range erf(\mu*r)/r in non-spin polarized case
! ref1: "Hybrid Functionals based on a screened Coulomb potential" PRB 2003, vol 118 p 8207
! ref2: "Generalized gradient approximation to the angle- and system-averaged exchange hole"
! JCP 1998, vol 109 p 3313
!Author: Iann Gerber
!Date : 01-04
!=====================================================================
      MODULE wpbe
      USE prec
      USE constant
      USE xclib
      USE setexm_struct_def
      IMPLICIT NONE
      PRIVATE :: LOGRHO,S
! Table for the range separated enhancement factor
      REAL(q), ALLOCATABLE, PRIVATE, SAVE :: FS_RS(:,:)
      REAL(q), ALLOCATABLE, PRIVATE, SAVE :: SSPLINES_RHO_S(:,:)
      REAL(q), ALLOCATABLE, PRIVATE, SAVE :: SSPLINES_S_RHO(:,:)
! Dimensions of the table
      INTEGER, PARAMETER, PRIVATE :: NRHO=2000
      REAL(q), PARAMETER, PRIVATE :: INV_NRHOMIN1=1._q/(NRHO-1)
      INTEGER, PARAMETER, PRIVATE :: NS=2000
      REAL(q), PARAMETER, PRIVATE :: INV_NSMIN1=1._q/(NS-1)
      INTEGER, PRIVATE, SAVE :: LOGRHO0

! Constants for enforcement of local Lieb-Oxford bound
! Value of the maximum of s without the rescaling
! define the value of MAXS in the table construction

      REAL(q), PARAMETER,PRIVATE :: STRANS=8.3_q
      REAL(q), PARAMETER,PRIVATE :: SMAX=8.5728844_q
      REAL(q), PARAMETER,PRIVATE :: SCONST=18.79622316_q

      REAL(q), PARAMETER, PRIVATE :: MAXLOGRHO= 22.5_q
      REAL(q), PARAMETER, PRIVATE :: MINLOGRHO=-33.0_q
!      REAL(q), PARAMETER, PRIVATE :: MAXS=(STRANS+3.0_q)
      REAL(q), PARAMETER, PRIVATE :: MAXS=(9.3_q)
      REAL(q), PARAMETER, PRIVATE :: MINS= 0.0_q
      REAL(q), PARAMETER, PRIVATE :: LOGRHOSTEP=(MAXLOGRHO-MINLOGRHO)/(1._q*(NRHO-1))
      REAL(q), PARAMETER, PRIVATE :: INV_LOGRHOSTEP=1._q/LOGRHOSTEP
      REAL(q), PARAMETER, PRIVATE :: SSTEP=MAXS/(1._q*(NS-1))
      REAL(q), PARAMETER, PRIVATE :: INV_SSTEP=1._q/SSTEP
!$ACC DECLARE CREATE(LOGRHO0,SSPLINES_RHO_S,SSPLINES_S_RHO,FS_RS)

      INTEGER, PRIVATE, SAVE :: LEXCH_WPBE_TABLE=-1
      REAL(q), PRIVATE, SAVE :: LDASCREEN_WPBE_TABLE=-1
      LOGICAL, PRIVATE, SAVE :: INITIALIZE_WPBE_TABLE=.TRUE.

      CONTAINS

!************************ SUBROUTINE INIT_WPBE *************************
!
! (Re)initialize the wpbe interpolation tables
!
!***********************************************************************

      SUBROUTINE INIT_WPBE
! Early exit if possible
      IF ((LEXCH/=8.AND.LEXCH/=9.AND.LEXCH/=14.AND.(LEXCH<18.OR.LEXCH>50)).OR. &
           LDASCREEN<=0.OR.LUSE_THOMAS_FERMI) RETURN
! Possibly reinitialize the table
      IF (.NOT.INITIALIZE_WPBE_TABLE) THEN
         IF (LDASCREEN_WPBE_TABLE/=LDASCREEN.OR.LEXCH_WPBE_TABLE/=LEXCH) INITIALIZE_WPBE_TABLE=.TRUE.
      ENDIF
! Deallocate the interpolation tables if allocated
      IF (INITIALIZE_WPBE_TABLE.AND.ALLOCATED(FS_RS)) DEALLOCATE(FS_RS,SSPLINES_RHO_S,SSPLINES_S_RHO)
! (Re)initialize the wpbe tables
      CALL INIT_WPBE_TABLE
      RETURN
      END SUBROUTINE INIT_WPBE


!************************ SUBROUTINE INIT_WPBE_TABLE *******************
!
!***********************************************************************

      SUBROUTINE INIT_WPBE_TABLE
      IMPLICIT NONE
      INTEGER IGRID,JGRID
      REAL(q) RHO
      REAL(q) FXWPBE_SR
! For testing purposes
      REAL(q) DFSR_DR,DFSR_DS,LGRHO,SV

! early exit if possible
      IF (.NOT.INITIALIZE_WPBE_TABLE) RETURN

!$OMP CRITICAL (omp_lck_wpbe_table)
!$    IF (INITIALIZE_WPBE_TABLE) THEN
! Allocate tables
         ALLOCATE(FS_RS(NRHO,NS),SSPLINES_RHO_S(NRHO,NS),SSPLINES_S_RHO(NRHO,NS))
! Fill the table containing the short range part of the
! range separated enhancement factor
         DO IGRID=1,NRHO
            IF ((LOGRHO(IGRID)>0._q).AND.((LOGRHO(IGRID)-LOGRHOSTEP)<0._q)) THEN
               LOGRHO0=IGRID
            ENDIF
            RHO=EXP(LOGRHO(IGRID))
            DO JGRID=1,NS
!gK: converted LDASCREEN to  atomic units
!LS chose between PBE (with EP exchange hole) and PBEsol (HJS exchange hole)
               IF (LEXCH==14) THEN
                  CALL EXCHWPBEsol_R(LDASCREEN*AUTOA,RHO,S(JGRID),FXWPBE_SR)
                  FS_RS(IGRID,JGRID)=FXWPBE_SR
               ELSE
                  CALL EXCHWPBE_R(LDASCREEN*AUTOA,RHO,S(JGRID),FXWPBE_SR)
                  FS_RS(IGRID,JGRID)=FXWPBE_SR
               ENDIF
            ENDDO
         ENDDO
! Calculate the 2nd derivative arrays
         CALL WPBE_SPLIE2
! Store the parameters for which the table was set up
         LEXCH_WPBE_TABLE=LEXCH
         LDASCREEN_WPBE_TABLE=LDASCREEN
! and flag the tables as initialized
         INITIALIZE_WPBE_TABLE = .FALSE.
!$    ENDIF
!$OMP END CRITICAL (omp_lck_wpbe_table)
!$ACC UPDATE DEVICE(LOGRHO0,FS_RS,SSPLINES_RHO_S,SSPLINES_S_RHO)

      RETURN
      END SUBROUTINE INIT_WPBE_TABLE


      FUNCTION LOGRHO(I)
!$ACC ROUTINE SEQ
      IMPLICIT NONE
      INTEGER I
      REAL(q) LOGRHO
      LOGRHO=MINLOGRHO+(I-1)*(MAXLOGRHO-MINLOGRHO)/(NRHO-1)
      END FUNCTION LOGRHO

      FUNCTION S(I)
!$ACC ROUTINE SEQ
      IMPLICIT NONE
      INTEGER I
      REAL(q) S
      S=MINS+(I-1)*(MAXS-MINS)/(NS-1)
! Rescaling the s values to ensure the Lieb-Oxford bound for s>8.3
!      IF (S>STRANS) S=SMAX-(SCONST/S**2)
      END FUNCTION S


      SUBROUTINE WPBE_SPLIE2
      IMPLICIT NONE
      INTEGER IGRID,JGRID
      REAL(q), ALLOCATABLE :: XTMP(:),YTMP(:),Y2TMP(:)

! Calculate SSPLINES_RHO_S: second derivative of FS_RS wrt. S at constant RHO
      ALLOCATE(XTMP(NS),YTMP(NS),Y2TMP(NS))
      DO IGRID=1,NRHO
         DO JGRID=1,NS
            YTMP(JGRID)=FS_RS(IGRID,JGRID)
            XTMP(JGRID)=S(JGRID)
         ENDDO
         CALL WPBE_SPLIN2(XTMP,YTMP,NS,1E30_q,1E30_q,Y2TMP)
         DO JGRID=1,NS
            SSPLINES_RHO_S(IGRID,JGRID)=Y2TMP(JGRID)
         ENDDO
      ENDDO
      DEALLOCATE(XTMP,YTMP,Y2TMP)
! Calculate SSPLINES_S_RHO: second derivative of FS_RS wrt. RHO at constant S
      ALLOCATE(XTMP(NRHO),YTMP(NRHO),Y2TMP(NRHO))
      DO IGRID=1,NS
         DO JGRID=1,NRHO
            YTMP(JGRID)=FS_RS(JGRID,IGRID)
            XTMP(JGRID)=LOGRHO(JGRID)
         ENDDO
         CALL WPBE_SPLIN2(XTMP,YTMP,NRHO,1E30_q,1E30_q,Y2TMP)
         DO JGRID=1,NRHO
            SSPLINES_S_RHO(IGRID,JGRID)=Y2TMP(JGRID)
         ENDDO
      ENDDO
      DEALLOCATE(XTMP,YTMP,Y2TMP)
      RETURN
      END SUBROUTINE WPBE_SPLIE2


      SUBROUTINE WPBE_SPLIN2(X,Y,N,YP1,YPN,Y2)
!$ACC ROUTINE SEQ
      IMPLICIT NONE
      INTEGER N
      REAL(q) X(N),Y(N),Y2(N)
      REAL(q) YP1,YPN

      INTEGER I,K
      REAL(q)  P,QN,SIG,UN,U(N)

! First point
      IF (YP1.GT.0.99E30_q) THEN
         Y2(1)=0; U(1)=0
      ELSE
!        Y2(1) =-0.5_q
         Y2(1) = 0.5_q
         U(1)  = (3._q/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
! Decomposition loop for the tridiagonal alg
      DO I=2,N-1
         SIG   = (X(I)-X(I-1))/(X(I+1)-X(I-1))
         P     = SIG*Y2(I-1)+2._q
         Y2(I) = (SIG-1._q)/P
         U(I)  = (6._q*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) &
        &        /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      ENDDO
! Last point
      IF (YPN.GT.0.99E30_q) THEN
         QN=0; UN=0
      ELSE
         QN = 0.5_q
         UN = (3._q/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N) = (UN-QN*U(N-1))/(QN*Y2(N-1)+1._q)
! Backsubstitution loop of the tridiagonal alg
      DO K=N-1,1,-1
         Y2(K) = Y2(K)*Y2(K+1)+U(K)
      ENDDO

      RETURN
      END SUBROUTINE WPBE_SPLIN2


      SUBROUTINE WPBE_SPLINE(X1,X2,Y)
!$ACC ROUTINE SEQ
      USE string, ONLY : str
      USE tutor, ONLY : vtutor
      IMPLICIT NONE

      REAL(q) X1,X2,Y
      REAL(q) D1X1,D1X2

      INTEGER J,K,JLO,JHI,KLO,KHI,JTMP,KTMP,KMIN

      REAL(q) YYTMPX2(6),Y2X2(6),LOGRHOTMP(6)
      REAL(q) YYTMPX1(6),Y2X1(6),STMP(6)
      REAL(q) H,A,B
      REAL(q) D1HI1,D1LO2,DIFFX1,DIFFX2
      REAL(q) D1A,D1B,H13,X_XJ_3_H,XJ1_X_3_H
      REAL(q) YX2,FX
!$ACC ROUTINE (WPBE_SPLIN2,LOGRHO) SEQ
! Rescaled the S parameter to lie in the bounds:
      IF (X2>STRANS) X2=SMAX-(SCONST/X2**2)
! Check whether X1 and X2 are within bounds
      IF (X1<LOGRHO(1).OR.X1>LOGRHO(NRHO)) THEN
         WRITE(*,*) 'WPBE_SPLINE: ERROR(1), LOG(RHO) out of bounds: ', &
        & X1,LOGRHO(1),LOGRHO(NRHO)
      ENDIF
      IF (X2<S(1).OR.X2>S(NS)) THEN
         WRITE(*,*) 'WPBE_SPLINE: ERROR(1), S out of bounds: ', &
        & X2,S(1),S(NS),EXP(X1)
      ENDIF
! Find bracketing indices at the rho axis
      IF (X1.GE.0._q) THEN
         JLO=INT(X1/LOGRHOSTEP)+LOGRHO0
      ELSE
         JLO=INT(ABS((X1-MINLOGRHO)/LOGRHOSTEP))+1
      ENDIF
      JHI=JLO+1

      IF ((LOGRHO(JLO).GT.X1).OR.(LOGRHO(JHI).LT.X1)) THEN
         IF ((LOGRHO(JLO).GT.X1)) THEN
            JLO=JLO-1
         ELSE
            JLO=JLO+1
         ENDIF
         JHI=JLO+1
      ENDIF

      IF ((LOGRHO(JLO).GT.X1).OR.(LOGRHO(JHI).LT.X1)) THEN
NOACC    CALL vtutor%error('WPBE_SPLINE: ERROR(2), LOG(RHO) out of bounds: ' // &
NOACC       str(X1) // ' ' // str(LOGRHO(JLO)) // ' ' // str(LOGRHO(JHI)) // &
NOACC       ' ' // str(JLO) // ' ' // str(JHI) // ' ' // str(LOGRHO0))
      ENDIF

! Find bracketing indices at the S axis
      KLO=INT((X2/SSTEP))+1
      KHI=KLO+1

      IF (S(KLO).GT.X2) THEN
         KHI=KLO
         KLO=KLO-1
      ENDIF

      IF ((S(KLO).GT.X2).OR.(S(KHI)-X2.LT.-1.0E-8_q)) THEN
NOACC    CALL vtutor%error('WPBE_SPLINE: ERROR(2), S out of bounds: ' // str(X2) // ' ' // &
NOACC       str(S(KLO)) // ' ' // str(S(KHI)) // ' ' // str(KLO) // ' ' // str(KHI))
         IF (S(KLO).GT.X2) THEN
            KLO=NS-4
            KHI=NS-3
         ENDIF
      ENDIF

#if defined(RETAIN_DEAD_CODE)
! Interpolate S on a section of the rows
      DO J=1,6
         ! Calculate function value
         H=S(KHI)-S(KLO)
         A=(S(KHI)-X2)/H
         B=(X2-S(KLO))/H

         JTMP=JLO-3+J
         YYTMPX1(J)=A*FS_RS(JTMP,KLO)+B*FS_RS(JTMP,KHI)+ &
        &  ((A*A*A-A)*SSPLINES_RHO_S(JTMP,KLO)+ &
        &   (B*B*B-B)*SSPLINES_RHO_S(JTMP,KHI))*(H*H)/6._q
      ENDDO
! Calculate a spline from these 6 function values
      DIFFX1=LOGRHOSTEP
      D1LO2 =((YYTMPX1(2)-YYTMPX1(1))/DIFFX1)
      D1HI1 =((YYTMPX1(6)-YYTMPX1(5))/DIFFX1)

      DO J=1,6
         LOGRHOTMP(J)=LOGRHO(JLO-3+J)
      ENDDO
      CALL WPBE_SPLIN2(LOGRHOTMP,YYTMPX1,6,D1LO2,D1HI1,Y2X1)
#endif

! Now do the same thing along the other dimension
      ! Take care of X2 close to zero
      IF (KLO.LT.3) THEN
         KMIN=2
         YYTMPX2(1)=1._q
         IF (KLO.LT.2) THEN
            KMIN=3
         ENDIF
      ELSE
         KMIN=1
      ENDIF

      DO K=KMIN,6
         ! Calculate function value
         H=LOGRHO(JHI)-LOGRHO(JLO)
         A=(LOGRHO(JHI)-X1)/H
         B=(X1-LOGRHO(JLO))/H

         KTMP=KLO-3+K

         YYTMPX2(K)=A*FS_RS(JLO,KTMP)+B*FS_RS(JHI,KTMP)+ &
        &  ((A*A*A-A)*SSPLINES_S_RHO(KTMP,JLO)+ &
        &   (B*B*B-B)*SSPLINES_S_RHO(KTMP,JHI))*(H*H)/6._q
      ENDDO

      IF (KMIN.GT.1) THEN
         IF (KMIN.EQ.2) THEN
            YYTMPX2(1)=YYTMPX2(3)
         ENDIF
         IF (KMIN.EQ.3) THEN
            YYTMPX2(2)=YYTMPX2(4)
            YYTMPX2(1)=YYTMPX2(5)
         ENDIF
      ENDIF

! Calculate a spline from these 6 function values
      DIFFX2=SSTEP
      D1LO2 =((YYTMPX2(2)-YYTMPX2(1))/DIFFX2)
      D1HI1 =((YYTMPX2(6)-YYTMPX2(5))/DIFFX2)

      IF (KMIN.EQ.1) THEN
         DO K=1,6
            STMP(K)=S(KLO-3+K)
         ENDDO
         CALL WPBE_SPLIN2(STMP,YYTMPX2,6,D1LO2,D1HI1,Y2X2)
      ELSE
         IF (KLO.EQ.2) THEN
            STMP(1)=S(KLO-1)-SSTEP
         ENDIF
         IF (KLO.EQ.1) THEN
            STMP(1)=S(KLO)-2._q*SSTEP
            STMP(2)=S(KLO)-SSTEP
         ENDIF
         DO K=KMIN,6
            STMP(K)=S(KLO-3+K)
         ENDDO
         CALL WPBE_SPLIN2(STMP,YYTMPX2,6,D1LO2,D1HI1,Y2X2)
      ENDIF
#if defined(RETAIN_DEAD_CODE)
! Interpolate the first spline for X1 (LOGRHO)
      H=LOGRHOSTEP
! Derivative in X1 direction
      D1A = -1._q/H
      D1B = -D1A
      H13 = (H*H*H)**(-1._q/3._q)

      XJ1_X_3_H=(LOGRHO(JHI)-X1)*(LOGRHO(JHI)-X1)*3._q*H13
      X_XJ_3_H =(X1-LOGRHO(JLO))*(X1-LOGRHO(JLO))*3._q*H13
! dF/dRHO
!      D1X1=(D1A*YYTMPX1(3)+D1B*YYTMPX1(4)+ &
!     &      (H*H/6._q)*((-XJ1_X_3_H-D1A)*Y2X1(3)+ &
!     &      (X_XJ_3_H-D1B)*Y2X1(4)))

! Interpolate the second spline for X2 (S)
#endif
      H=SSTEP
      A=(S(KHI)-X2)/H
      B=(X2-S(KLO))/H
! F(rho,s)
      YX2=A*YYTMPX2(3)+B*YYTMPX2(4)+ &
     &   ((A*A*A-A)*Y2X2(3)+(B*B*B-B)*Y2X2(4))*(H*H)/6._q
#if defined(RETAIN_DEAD_CODE)
! Derivative in X2 direction
      D1A = -1._q/H
      D1B = -D1A
      H13 = (H*H*H)**(-1._q/3._q)

      XJ1_X_3_H=(S(KHI)-X2)*(S(KHI)-X2)*3._q*H13
      X_XJ_3_H =(X2-S(KLO))*(X2-S(KLO))*3._q*H13
! dF/dS
!      D1X2=(D1A*YYTMPX2(3)+D1B*YYTMPX2(4)+ &
!     &      (H*H/6._q)*((-XJ1_X_3_H-D1A)*Y2X2(3)+ &
!     &      (X_XJ_3_H-D1B)*Y2X2(4)))
#endif
      Y=YX2

!     IF (Y.LT.0._q) THEN
!        WRITE(*,*) 'WPBE_SPLINE: Y<1',X1,X2,Y
!     ENDIF
      RETURN
      END SUBROUTINE WPBE_SPLINE

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE WPBE_SPLIN2_N6_SIMD(X, Y, YP1, YPN, Y2, MASK0)
      USE simd
      IMPLICIT NONE
      TYPE(SIMD_REAL64), INTENT(IN)      :: X(*), Y(*)
      TYPE(SIMD_REAL64), INTENT(IN)      :: YP1, YPN
      TYPE(SIMD_REAL64), INTENT(OUT)     :: Y2(*)
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: U(6), XN_SUB_XNM1(6), YN_SUB_YNM1(6)
      REAL(q)                            :: P, INV_P, QN, SIG, UN
      REAL(q)                            :: TEMP_1, TEMP_2, TEMP_3, TEMP_4, TEMP_5
      INTEGER                            :: I, II, K

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, TEMP_2, TEMP_3, TEMP_4, SIG, P, INV_P, QN, UN)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         DO K = 2, 6
            XN_SUB_XNM1(K)%X(II) = X(K)%X(II) - X(K - 1)%X(II)
            YN_SUB_YNM1(K)%X(II) = Y(K)%X(II) - Y(K - 1)%X(II)
         ENDDO
! First point
         IF (MASK0%X(II)) THEN
            Y2(1)%X(II) = 0.0_q
         ENDIF
         U(1)%X(II) = 0.0_q
         IF (YP1%X(II) .LE. 0.99E30_q) THEN
            IF (MASK0%X(II)) THEN
               Y2(1)%X(II) = 0.5_q
            ENDIF
            TEMP_1 = 1.0_q / XN_SUB_XNM1(2)%X(II)
            U(1)%X(II) = (3.0_q * TEMP_1) * (YN_SUB_YNM1(2)%X(II) * TEMP_1 - YP1%X(II))
         ENDIF
! Decomposition loop for the tridiagonal alg
         TEMP_2 = YN_SUB_YNM1(2)%X(II) / XN_SUB_XNM1(2)%X(II)
         DO I = 2, 5
            TEMP_1 = TEMP_2
            TEMP_3 = 1.0_q / XN_SUB_XNM1(I + 1)%X(II)
            TEMP_2 = YN_SUB_YNM1(I + 1)%X(II) * TEMP_3
            TEMP_4 = 1.0_q / (X(I + 1)%X(II) - X(I - 1)%X(II))
            SIG = XN_SUB_XNM1(I)%X(II) * TEMP_4
            P = SIG * Y2(I - 1)%X(II) + 2.0_q
            INV_P = 1.0_q / P
            IF (MASK0%X(II)) THEN
               Y2(I)%X(II) = (SIG - 1.0_q) * INV_P
            ENDIF
            U(I)%X(II) = (6.0_q * (TEMP_2 - TEMP_1) * TEMP_4 - SIG * U(I - 1)%X(II)) * INV_P
         ENDDO
! Last point
         QN = 0.0_q
         UN = 0.0_q
         IF (YPN%X(II).LE.0.99E30_q) THEN
            QN = 0.5_q
            UN = (3.0_q * TEMP_3) * (YPN%X(II) - TEMP_2)
         ENDIF
         Y2(6)%X(II) = (UN - QN * U(5)%X(II)) / (QN * Y2(5)%X(II) + 1.0_q)
! Backsubstitution loop of the tridiagonal alg
         DO K = 5, 1, -1
            IF (MASK0%X(II)) THEN
               Y2(K)%X(II) = Y2(K)%X(II) * Y2(K + 1)%X(II) + U(K)%X(II)
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE WPBE_SPLIN2_N6_SIMD


      SUBROUTINE WPBE_SPLINE_SIMD(X1_IN, X2_IN, Y, MASK0)
      USE simd
      IMPLICIT NONE
      TYPE(SIMD_REAL64), INTENT(IN)      :: X1_IN, X2_IN
      TYPE(SIMD_REAL64), INTENT(OUT)     :: Y
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: X1, X2
      TYPE(SIMD_REAL64)                  :: YYTMPX2(6), Y2X2(6), STMP(6)
      TYPE(SIMD_REAL64)                  :: D1HI1, D1LO2, S_KLO, S_KHI, LOGRHO_JLO, LOGRHO_JHI
      TYPE(SIMD_REAL64)                  :: TEMP_FS_RS_JLO(6), TEMP_FS_RS_JHI(6)
      TYPE(SIMD_REAL64)                  :: TEMP_SSPLINES_S_RHO_JLO(6), TEMP_SSPLINES_S_RHO_JHI(6)
      TYPE(SIMD_INT64)                   :: JLO, JHI, KLO, KHI
      TYPE(SIMD_MASK_REAL64)             :: MASK1
      LOGICAL                            :: TRUE_FOR_ANY
      REAL(q)                            :: D1X1, D1X2
      REAL(q)                            :: H, INV_H, A, B
      INTEGER                            :: II, J, K, KTMP
      LOGICAL                            :: PRED1, PRED2
      INTEGER, PARAMETER                 :: INT1 = 1, INT3 = 3
      REAL(q), PARAMETER                 :: P_1 = 0.166666666666666657414808_q ! 1._q/6._q

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(PRED1, PRED2)
      DO II = 0, SIMD_WIDTH_64BIT - 1
! Rescaled the S parameter to lie in the bounds:
         X1%X(II) = X1_IN%X(II)
         X2%X(II) = X2_IN%X(II)
         IF (X2%X(II) .GT. STRANS) X2%X(II) = SMAX - SCONST / (X2%X(II) * X2%X(II))
! Check whether X1 and X2 are within bounds
!        IF (X1.LT.LOGRHO(1).OR.X1.GT.LOGRHO(NRHO)) THEN
!           WRITE(*,*) 'WPBE_SPLINE: ERROR(1), LOG(RHO) out of bounds: ', &
!           &  X1,LOGRHO(1),LOGRHO(NRHO)
!        ENDIF
!        IF (X2.LT.S(1).OR.X2.GT.S(NS)) THEN
!           WRITE(*,*) 'WPBE_SPLINE: ERROR(1), S out of bounds: ', &
!           &  X2,S(1),S(NS),EXP(X1)
!        ENDIF
! Find bracketing indices at the rho axis
         IF (X1%X(II) .GE. 0.0_q) THEN
            JLO%X(II) = INT(X1%X(II) * INV_LOGRHOSTEP) + LOGRHO0
         ELSE
            JLO%X(II) = INT(ABS((X1%X(II) - MINLOGRHO) * INV_LOGRHOSTEP)) + 1
         ENDIF
         JHI%X(II) = JLO%X(II) + 1
         LOGRHO_JLO%X(II) = MINLOGRHO + (JLO%X(II) - 1) * (MAXLOGRHO - MINLOGRHO) * INV_NRHOMIN1
         LOGRHO_JHI%X(II) = MINLOGRHO + (JHI%X(II) - 1) * (MAXLOGRHO - MINLOGRHO) * INV_NRHOMIN1
         PRED1 = LOGRHO_JLO%X(II) .GT. X1%X(II)
         PRED2 = LOGRHO_JHI%X(II) .LT. X1%X(II)
         IF (PRED1) JLO%X(II) = JLO%X(II) - 1
         IF (PRED2) JLO%X(II) = JLO%X(II) + 1
         IF (PRED1 .OR. PRED2) THEN
            JHI%X(II) = JLO%X(II) + 1
            LOGRHO_JLO%X(II) = MINLOGRHO + (JLO%X(II) - 1) * (MAXLOGRHO - MINLOGRHO) * INV_NRHOMIN1
            LOGRHO_JHI%X(II) = MINLOGRHO + (JHI%X(II) - 1) * (MAXLOGRHO - MINLOGRHO) * INV_NRHOMIN1
         ENDIF

!        IF ((LOGRHO(JLO).GT.X1).OR.(LOGRHO(JHI).LT.X1)) THEN
!           WRITE(*,*) 'WPBE_SPLINE: ERROR(2), LOG(RHO) out of bounds: ', &
!           &  X1,LOGRHO(JLO),LOGRHO(JHI),JLO,JHI,LOGRHO0
!           STOP
!        ENDIF
         MASK1%X(II) = MASK0%X(II)
         IF (MASK1%X(II)) THEN
            IF ((LOGRHO_JLO%X(II) .GT. X1%X(II)) .OR. (LOGRHO_JHI%X(II) .LT. X1%X(II))) THEN
               MASK1%X(II) = .FALSE.
            ENDIF
         ENDIF

! Find bracketing indices at the S axis
         KLO%X(II) = INT(X2%X(II) * INV_SSTEP) + 1
         KHI%X(II) = KLO%X(II) + 1
         S_KLO%X(II) = MINS + (KLO%X(II) - 1) * (MAXS - MINS) * INV_NSMIN1
         IF (S_KLO%X(II) .GT. X2%X(II)) THEN
            KHI%X(II) = KLO%X(II)
            KLO%X(II) = KLO%X(II) - 1
            S_KLO%X(II) = MINS + (KLO%X(II) - 1) * (MAXS - MINS) * INV_NSMIN1
         ENDIF
         S_KHI%X(II) = MINS + (KHI%X(II) - 1) * (MAXS - MINS) * INV_NSMIN1
!        IF ((S(KLO).GT.X2).OR.(S(KHI)-X2.LT.-1.0E-8_q)) THEN
!           WRITE(*,*) 'WPBE_SPLINE: ERROR(2), S out of bounds: ', &
!           &  X2,S(KLO),S(KHI),KLO,KHI
!           STOP
!           IF (S(KLO).GT.X2) THEN
!              KLO=NS-4
!              KHI=NS-3
!           ENDIF
!        ENDIF

         IF ((S_KLO%X(II) .GT. X2%X(II)) .OR. ((S_KHI%X(II) - X2%X(II)) .LT. -1.0E-8_q)) THEN
            MASK1%X(II) = .FALSE.
         ENDIF

         DO K = 1, 6
            STMP(K)%X(II) = MINS + (KLO%X(II) - 3 + K - 1) * (MAXS - MINS) * INV_NSMIN1
         ENDDO
         IF (KLO%X(II) .EQ. 2) THEN
            STMP(1)%X(II) = MINS - SSTEP
         ENDIF
         IF (KLO%X(II) .EQ. 1) THEN
            STMP(1)%X(II) = MINS - 2.0_q * SSTEP
            STMP(2)%X(II) = MINS - SSTEP
         ENDIF
      ENDDO

      TRUE_FOR_ANY = .FALSE.
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (MASK1%X(II)) THEN
            TRUE_FOR_ANY = .TRUE.
         ENDIF
      ENDDO
      
      IF (.NOT. TRUE_FOR_ANY) RETURN
      

      DO II = 0, SIMD_WIDTH_64BIT - 1
! Now do the same thing along the other dimension
! Take care of X2 close to zero
         DO K = 1, 6
            KTMP = MAX(INT1, KLO%X(II) - INT3 + K)
            IF (MASK1%X(II)) THEN
               TEMP_FS_RS_JLO(K)%X(II) = FS_RS(JLO%X(II), KTMP)
               TEMP_FS_RS_JHI(K)%X(II) = FS_RS(JHI%X(II), KTMP)
               TEMP_SSPLINES_S_RHO_JLO(K)%X(II) = SSPLINES_S_RHO(KTMP, JLO%X(II))
               TEMP_SSPLINES_S_RHO_JHI(K)%X(II) = SSPLINES_S_RHO(KTMP, JHI%X(II))
            ENDIF
         ENDDO
      ENDDO

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(H, INV_H, A, B)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         DO K = 1, 6
! Calculate function value
            H = LOGRHO_JHI%X(II) - LOGRHO_JLO%X(II)
            INV_H = 1.0_q / H
            A = (LOGRHO_JHI%X(II) - X1%X(II)) * INV_H
            B = (X1%X(II) - LOGRHO_JLO%X(II)) * INV_H
            YYTMPX2(K)%X(II) = A * TEMP_FS_RS_JLO(K)%X(II) + B * TEMP_FS_RS_JHI(K)%X(II) + &
            &  ((A * A * A - A) * TEMP_SSPLINES_S_RHO_JLO(K)%X(II) + &
            &   (B * B * B - B) * TEMP_SSPLINES_S_RHO_JHI(K)%X(II)) * (H * H) * P_1
         ENDDO
         IF (KLO%X(II) .EQ. 2) THEN
            YYTMPX2(1)%X(II) = YYTMPX2(3)%X(II)
         ENDIF
         IF (KLO%X(II) .EQ. 1) THEN
            YYTMPX2(2)%X(II) = YYTMPX2(4)%X(II)
            YYTMPX2(1)%X(II) = YYTMPX2(5)%X(II)
         ENDIF
! Calculate a spline from these 6 function values
         D1LO2%X(II) = ((YYTMPX2(2)%X(II) - YYTMPX2(1)%X(II)) * INV_SSTEP)
         D1HI1%X(II) = ((YYTMPX2(6)%X(II) - YYTMPX2(5)%X(II)) * INV_SSTEP)
      ENDDO

      CALL WPBE_SPLIN2_N6_SIMD(STMP, YYTMPX2, D1LO2, D1HI1, Y2X2, MASK1)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(A, B)
      DO II = 0, SIMD_WIDTH_64BIT - 1
! Interpolate the second spline for X2 (S)
         A = (S_KHI%X(II) - X2%X(II)) * INV_SSTEP
         B = (X2%X(II) - S_KLO%X(II)) * INV_SSTEP
! F(rho,s)
         IF (MASK1%X(II)) THEN
            Y%X(II) = A * YYTMPX2(3)%X(II) + B * YYTMPX2(4)%X(II) + &
            &  ((A * A * A - A) * Y2X2(3)%X(II) + (B * B * B - B) * Y2X2(4)%X(II)) * (SSTEP * SSTEP) * P_1
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE WPBE_SPLINE_SIMD
#endif


!=====================================================================
!Functional wPBE long-range for calculation in separated interactions LR-SR
!with modified kernel Short-Range erf(\mu*r)/r
! ref1: "Hybrid Functionals based on a screened Coulomb potential" PRB 2003, vol 118 p 8207
! ref2: "Generalized gradient approximation to the angle- and system-averaged exchange hole"
! JCP 1998, vol 109 p 3313
! ref 3 jcp 2004 vol 120 7274
!Author: Iann Gerber
!Date : 05-04
!=====================================================================
      SUBROUTINE EXCHWPBE_R(OMEGA,RHO,SW,FXWPBE_SR)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  INPUT rho : DENSITY
!  INPUT OMEGA : SPLITTING RANGE
!  INPUT rhothrd : DENSITY^(1/3)
!  INPUT S:  (GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Formulas:
!       e_x[unif]=ax*rho^(4/3)  [LDA]
! ax = -0.75*(3/pi)^(1/3)
!       e_x[PBE]=e_x[unif]*FxwPBE(s,w,rho)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      USE prec
      USE constant
      USE c2f_interface, ONLY : ERRF,ERRFC
!      USE xclib
      IMPLICIT NONE
      REAL(q) OMEGA,RHO
      REAL(q) FXWPBE_SR
      REAL(q) F13,F43,F12,F14,F32,F34,F94,F98,F1516
      REAL(q) AX,UM,UK,UL,P0,FXPBE
      REAL(q) A,B,C,D,E
      REAL(q) HA1,HA2,HA3,HA4,HA5
      REAL(q) FC1,FC2,EA1,EA2,EA3,EA4,EA5,EA6,EA7,EA8,EB1,WCUT
      REAL(q) EGSCUT,EGA1,EGA2,EGA3
      REAL(q) EXPCUT,EXEI1,EXEI2,EXEI3,EXEI4
      REAL(q) PI2,SRPI,F89M,SREAL,STRA1S
      REAL(q) A2,A3,A4,A12,A32,A52,W,W2,W3,W4,W5,W6,W7,W8,XKF
      REAL(q) SW,S2,S3,S4,S5,S6
      REAL(q) HNUM,HDEN,H,HNU1S,HDE1S,H1S,F
      REAL(q) HSBW,HSBW2,HSBW3,HSBW4,HSBW6,HSBW12,HSBW32,HSBW52,HSBW72
      REAL(q) DHSB,DHSB2,DHSB3,DHSB4,DHSB5,DHSB12,DHSB32,DHSB52,DHSB72
      REAL(q) DHSB92,HA94,HA942,HA943,HA945,HA9412
      REAL(q) DHS,DHS2,DHS3,DHS4,DHS72,DHS92,DHSW,DHSW2,DHSW52,DHSW72
      REAL(q) GA,GB,EG,TM1,TM2,TM3,TM4,TM5,T10
      REAL(q) EXER,EXHA94,EIHA94,EXEI,T1,PN1,PN2
      REAL(q) F2,F3,F4,F5,F6,F7,F8,F9,T2T9
! Numerical factors
      PARAMETER(F13=1._q/3._q,F43=4._q/3._q,F12=0.5_q,F14=0.25_q)
      PARAMETER(F32=1.5_q,F34=0.75_q,F94=2.25_q,F98=1.125_q,F1516=0.9375_q)
! values for PBE enhancement factor calculation
      PARAMETER(UM=0.2195149727645171_q,UK=0.8040_q,UL=UM/UK)
! Constants  from the PBE hole
      PARAMETER(A=1.0161144_q,B=-0.37170836_q)
      PARAMETER(C=-0.077215461_q,D=0.57786348_q)
      PARAMETER(E=-0.051955731_q)
! Values for H(s)
      PARAMETER(HA1=0.00979681_q,HA2=0.0410834_q,HA3=0.187440_q)
      PARAMETER(HA4=0.00120824_q,HA5=0.0347188_q)
! Values for F(s)
      PARAMETER(FC1=6.4753871_q,FC2=0.47965830_q)
! Coefficients of the erfc(x) expansion (eb1 set later depending on wcut)
      PARAMETER(EA1=-1.128223946706117_q,EA2=1.452736265762971_q)
      PARAMETER(EA3=-1.243162299390327_q,EA4=0.971824836115601_q)
      PARAMETER(EA5=-0.568861079687373_q,EA6=0.246880514820192_q)
      PARAMETER(EA7=-0.065032363850763_q,EA8=0.008401793031216_q)
      PARAMETER(WCUT=14.0_q)
! Constants for polynomial expansion of EG for small s
      PARAMETER(EGSCUT=0.08_q,EGA1=-0.02628417880_q,EGA2=-0.07117647788_q)
      PARAMETER(EGA3=0.08534541323_q)
! Constants for large x in exp(x)*Ei(x)
      PARAMETER(EXPCUT=700_q,EXEI1=4.03640_q,EXEI2=1.15198_q)
      PARAMETER(EXEI3=5.03627_q,EXEI4=4.19160_q)

! General constants
      PI2=PI*PI
      SRPI=SQRT(PI)
      F89M=-8._q/9._q

!----------------------------------------------------------------------
! construct modified-PBE enhancement factor
!----------------------------------------------------------------------
!     INTERMEDIATE VARIABLES

! Calculate prelim variables
      XKF=(3._q*PI2*RHO)**F13
      A2=A*A ; A3=A2*A ; A4=A3*A ; A12=SQRT(A) ; A32=A12*A ; A52=A32*A
      W=OMEGA/XKF ; W2= W*W ; W3=W2*W ; W4=W2*W2 ; W5=W2*W3 ; W6=W3*W3
      W7=W6*W ; W8=W7*W

      S2=SW*SW
      S3=S2*SW
      S4=S2*S2
      S5=S4*SW
      S6=S4*S2

! Calculate H(s) and F(s) for the PBE hole
      HNUM=HA1*S2+HA2*S4
      HDEN=1.0_q+HA3*S4+HA4*S5+HA5*S6
      H=(HNUM)/(HDEN)
      HNU1S=2._q*HA1*SW+4._q*HA2*S3
      HDE1S=4.0_q*HA3*S3+5.0_q*HA4*S4+6.0_q*HA5*S5
      H1S=(HDEN*HNU1S-HNUM*HDE1S)/(HDEN*HDEN)
      F=FC1*H+FC2

! Set exponent of the Gaussian in the approximation of the erfc function
      IF (W<WCUT) THEN
         EB1= 1.455915450052607_q
      ELSE
         EB1=2.0_q
      ENDIF

! Calculate intermediate variables
      HSBW=S2*H+EB1*W2 ; HSBW2=HSBW*HSBW ; HSBW3=HSBW2*HSBW ;
      HSBW4=HSBW2*HSBW2 ; HSBW6=HSBW3*HSBW3
      HSBW12=SQRT(HSBW) ; HSBW32=HSBW12*HSBW ; HSBW52=HSBW32*HSBW
      HSBW72=HSBW52*HSBW

      DHSB=D+S2*H+EB1*W2 ; DHSB2=DHSB*DHSB ; DHSB3=DHSB2*DHSB
      DHSB4=DHSB2*DHSB2 ; DHSB5=DHSB4*DHSB
      DHSB12=SQRT(DHSB) ; DHSB32=DHSB12*DHSB ; DHSB52=DHSB32*DHSB
      DHSB72=DHSB52*DHSB ; DHSB92=DHSB72*DHSB

      HA94=F94*HSBW/A ; HA942=HA94*HA94 ; HA943=HA942*HA94
      HA945=HA943*HA942 ; HA9412=SQRT(HA94)

      DHS=D+S2*H ; DHS2=DHS*DHS ; DHS3=DHS2*DHS ; DHS4=DHS2*DHS2
      DHS72=DHS3*SQRT(DHS) ; DHS92=DHS72*DHS

      DHSW=DHS +W2 ; DHSW2=DHSW*DHSW ; DHSW52=SQRT(DHSW)*DHSW2
      DHSW72=DHSW52*DHSW

! Calculate G(s) using expansion for small s if necessary
      IF (SW>EGSCUT) THEN
         GA=SRPI*(15._q*E+6.0_q*C*(1.0_q+F*S2)*DHS + 4.0_q*B*(DHS2)&
        &   +8.0_q*A*(DHS3)) * (1.0_q/(16._q*DHS72))         &
        &   -F34*PI*SQRT(A)*EXP(F94*H*S2/A)* (1._q-ERRF(F32*SW*SQRT(H/A)))
         GB=F1516*SRPI*S2/DHS72
         EG=-(F34*PI+GA)/GB
      ELSE
         EG=EGA1+EGA2*S2+EGA3*S4
      ENDIF
! calculate the terms needed in any case

      TM2=(DHS2*B + DHS*C +2._q*E +DHS*S2*C*F +2._q*S2*EG )/2._q/DHS3
      TM3=-W*(4._q*DHSW2*B +6._q*DHSW*C + 15._q*E + 6.0_q*DHSW*S2*C*F + &
     &        15._q*S2*EG)/8._q/DHS/DHSW52
      TM4=-W3*(DHSW*C + 5._q*E + DHSW*S2*C*F + 5.0_q*S2*EG)/2._q/DHS2/DHSW52
      TM5=-W5*(E+S2*EG)/DHS3/DHSW52

! Calculate t10 unless that would generate a division by zero

      IF ((SW>0.0_q).OR.(W>0.0_q)) THEN
         T10=F12*A*LOG(HSBW/DHSB)
      ENDIF
! Calculate exp(x)*f(x) depending on the size of x

      IF (HA94<EXPCUT) THEN
         EXER=PI*EXP(HA94)*(ERRFC(HA9412))
         EXHA94=EXP(HA94)
         EIHA94=-EXPINT(1,HA94)
         EXEI=EXHA94*EIHA94
      ELSE
         EXER=PI*(1._q/(SRPI*HA9412)-1._q/(2._q*SQRT(PI*HA943))+ &
        &         3._q/(4._q*SQRT(PI*HA945)))
         EXEI=-(1._q/HA94)*(HA942+EXEI1*HA94+EXEI2)/            &
        &      (HA942+EXEI3*HA94+EXEI4)
      ENDIF
      IF (W==0.0_q) THEN
!Fall back to the PBE hole expression
         T1=-F12*A*EXEI
         IF (SW>0.0_q) THEN
            TM1=T1+T10
            FXWPBE_SR=F89M*(TM1+TM2)
          ELSE
            FXWPBE_SR=1._q
         ENDIF
      ELSE IF(W>WCUT) THEN
! Use simple gaussian approximation for large w
         TM1=-F12*A*(EXEI+LOG(DHSB)-LOG(HSBW))
         FXWPBE_SR=F89M*(TM1+TM2+TM3+TM4+TM5)
      ELSE
! For everything else use the full blown expression
!
! First calculate the polynomials for the first term
         PN1=-F32*EA1*A12*W + 27._q*EA3*W3/(8._q*A12)-243._q*EA5*W5/ &
        &     (32._q*A32) + 2187._q*EA7*W7/(128._q*A52)
         PN2=-A + F94*EA2*W2 - 81._q*EA4*W4/(16.0_q*A) + &
        &    729._q*EA6*W6/(64._q*A2) - 6561._q*EA8*W8/(256._q*A3)

! The first term is
         T1=F12*(PN1*EXER+PN2*EXEI)
! The factors for the main polynomials in w
         F2= F12*EA1*SRPI*A/DHSB12
         F3= F12*EA2*A/DHSB
         F4= EA3*SRPI*(-F98/HSBW12+F14*A/DHSB32)
         F5= EA4*(1._q/128._q)*(-144._q*(1._q/HSBW)+64._q*(1._q/DHSB2)*A)
         F6= EA5*(3._q*SRPI*(3._q*DHSB52*(9.0_q*HSBW-2._q*A)+4.0_q  &
        &            *HSBW32*A2))/(32._q*DHSB52*HSBW32*A)
         F7= EA6*(((32._q*A)/DHSB3 + (-36._q+(81._q*S2*H)/A)/HSBW2))&
        &      /32._q
         F8= EA7*(-3._q*SRPI*(-40._q*HSBW52*A3+9.0_q*DHSB72*(27._q  &
        &      *HSBW2-6.0_q*HSBW*A+4._q*A2)))/(128._q*DHSB72*HSBW52*A2)
         F9= (324._q*EA6*EB1*DHSB4*HSBW*A + EA8*(384._q*HSBW3*A3    &
        &        +DHSB4*(-729._q*HSBW2+324._q*HSBW*A-288._q*A2)))      &
        &        /(128._q*DHSB4*HSBW3*A2)

         T2T9= F2*W+F3*W2+F4*W3+F5*W4+F6*W5+F7*W6+F8*W7+f9*W8

! The final value of the first term for 0<omega<wcut is
         TM1= T1+ T2T9 +T10
         FXWPBE_SR=F89M*(TM1+TM2+TM3+TM4+TM5)
      ENDIF
!      write (*,*) RHO,S,FXWPBE_SR
      RETURN
      END SUBROUTINE

!LS
!
!=====================================================================
!Range-separated enhancementfactor for hybrid functionals
! ref1: "GGA model exchange holes for range-separated hybrids",
!        J.Chem.Phys. 128, 194105(2008)
! not yet implemented: stable version for \nu --> infinity
! J. Chem Theo. Comput. 5, 754 (2009)
!
!Date : 05-2009 ; Laurids Schimka
!=====================================================================
      SUBROUTINE EXCHWPBEsol_R(OMEGA,RHO,SW,FXWPBE_SR)
      USE prec
      USE constant
!      USE xclib
      IMPLICIT NONE
      REAL(q) OMEGA,RHO
      REAL(q) FXWPBE_SR
      REAL(q) F13,F43,F12,F14,F32,F34,F94,F98,F1516,F54
      REAL(q) ABAR,B,C,D,E
      REAL(q) HA2,HA3,HA4,HA5,HA6,HA7
      REAL(q) HB1,HB2,HB3,HB4,HB5,HB6,HB7,HB8,HB9
      REAL(q) PI2,SRPI,F89M,F25M,F415M,F65M,F45M,F125M
      REAL(q) F49M,F32M,F158M,F38M
      REAL(q) XKF
      REAL(q) SW,S2,S3,S4,S5,S6,S7,S8,S9
      REAL(q) nu,eta,lambda,zeta,chi
      REAL(q) lambda2,lambda3,lambda72
      REAL(q) chi2,chi3,chi5
      REAL(q) zeta12,eta12,nu2
      REAL(q) HNUM,HDEN,H,FBAR,EGBAR
      REAL(q) T1,T2,T3,T4,T5,T6



     ! Numerical factors
      PARAMETER(F13=1._q/3._q,F43=4._q/3._q,F12=0.5_q,F14=0.25_q)
      PARAMETER(F32=1.5_q,F34=0.75_q,F94=2.25_q,F98=1.125_q,F1516=0.9375_q)
      PARAMETER(F54=1.25_q)
     ! Constants  from the PBE(sol) hole
      PARAMETER(ABAR=0.757211_q,B=-0.106364_q,C=-0.118649_q,D=0.609650_q,E=-0.0477963_q)
! this parameters are for the PBE exchange hole
     ! Values for H(s)-PBE; values for PBEsol,B88 and B97x are also available
!      PARAMETER(HA2=0.0159941_q,HA3=0.0852995_q,HA4=-0.160368_q,HA5=0.152645_q)
!      PARAMETER(HA6=-0.0971263_q,HA7=0.0422061_q)
!      PARAMETER(HB1=5.33319_q,HB2=-12.4780_q,HB3=11.0988_q,HB4=-5.11013_q,HB5=1.71468_q)
!      PARAMETER(HB6=-0.610380_q,HB7=0.307555_q,HB8=-0.0770547_q,HB9=0.0334840_q)

     ! Parameters for PBEsol
      PARAMETER(HA2=0.0047333_q,HA3=0.0403304_q,HA4=-0.0574615_q,HA5=0.0435395_q)
      PARAMETER(HA6=-0.0216251_q,HA7=0.0063721_q)
      PARAMETER(HB1=8.52056_q,HB2=-13.9885_q,HB3=9.28583_q,HB4=-3.27287_q,HB5=0.843499_q)
      PARAMETER(HB6=-0.235543_q,HB7=0.0847074_q,HB8=-0.0171561_q,HB9=0.0050552_q)

     ! General constants
      PI2=PI*PI
      SRPI=SQRT(PI)
      F89M=-8._q/9._q;F25M=-2._q/5._q;F415M=-4._q/15._q;F65M=-6._q/5._q
      F45M=-4._q/5._q;F125M=-12._q/5._q;F49M=-4._q/9._q;F32M=-3._q/2._q
      F158M=-15._q/8._q;F38M=-3._q/8._q

!----------------------------------------------------------------------
! construct modified-PBE enhancement factor
!----------------------------------------------------------------------
!     INTERMEDIATE VARIABLES

! Calculate prelim variables
      XKF=(3._q*PI2*RHO)**F13  !(k_f)
      S2=SW*SW;S3=S2*SW;S4=S2*S2;S5=S4*SW;S6=S4*S2;S7=S6*SW;S8=S4*S4;S9=S8*SW
! Calculate H(s)
      HNUM=HA2*S2+HA3*S3+HA4*S4+HA5*S5+HA6*S6+HA7*S7
      HDEN=1._q+HB1*SW+HB2*S2+HB3*S3+HB4*S4+HB5*S5+HB6*S6+HB7*S7+HB8*S8+HB9*S9
      H=(HNUM)/(HDEN)
! Calculate Fbar(s)
      FBAR=1._q-1._q/(27._q*C)*S2/(1._q+S2/4._q)-1._q/(2._q*C)*S2*H
! Calculate chi,lambda,nu,zeta,eta
      nu=OMEGA/XKF;nu2=nu*nu
      zeta=S2*H;zeta12=SQRT(zeta)
      eta=ABAR+S2*H;eta12=SQRT(eta)
      lambda=D+S2*H;lambda2=lambda*lambda;lambda3=lambda2*lambda
      lambda72=lambda3*SQRT(lambda)
      chi=nu/(SQRT(lambda+nu2));chi2=chi*chi;chi3=chi2*chi;chi5=chi2*chi3
! Calculate EGbar(s)
      EGBAR=F25M*C*FBAR*lambda+F415M*B*lambda2+F65M*ABAR*lambda3+F45M*SRPI*lambda72+ &
     &F125M*lambda72*(zeta12-eta12)
! Calculate terms of FXWPBE_SR
      T1=ABAR+F49M*B/lambda*(1._q-chi)
      T2=F49M*C*FBAR/lambda2*(1._q+F32M*chi+F12*chi3)
      T3=F89M*EGBAR/lambda3*(1._q+F158M*chi+F54*chi3+F38M*chi5)
      T4=2._q*nu*(SQRT(zeta+nu2)-SQRT(eta+nu2))
      T5=2._q*zeta*LOG((nu+SQRT(zeta+nu2))/(nu+SQRT(lambda+nu2)))
      T6=-2._q*eta*LOG((nu+SQRT(eta+nu2))/(nu+SQRT(lambda+nu2)))
! Calculate FXWPBE_SR
      FXWPBE_SR=T1+T2+T3+T4+T5+T6
!      WRITE(77,'(5G16.7)') OMEGA,RHO,nu,SW,FXWPBE_SR
      RETURN
      END SUBROUTINE
!End LS

      FUNCTION EXPINT(N,X)
!
!       ============================================
!       Purpose: Compute exponential integral EN(x)
!       ============================================
!
      USE prec
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      REAL(q) EXPINT
      REAL(q) X,EULER,EPS,FPMIN
      REAL(q) B,C,D,H,A,DEL,FACT,PSI
      INTEGER N,MAXIT,I,II,NM1

      PARAMETER (EPS=1E-9_q,FPMIN=1E-30_q)
      PARAMETER (MAXIT=100)
      PARAMETER (EULER=0.577215664901532860606512_q)

      NM1=N-1
      IF ((N<0.OR.X<0.0_q).OR.(X==0._q.AND.(N==0.0_q.OR.N==1))) THEN
         CALL vtutor%stopCode()
      ELSE IF (N==0) THEN
         EXPINT=EXP(-X)/X
      ELSE IF (X==0.0_q) THEN
         EXPINT=1._q/NM1
      ELSE IF (X>1.0_q) THEN
         B=X+N
         C=1._q/FPMIN
         D=1.0_q/B
         H=D
         DO I=1,MAXIT
            A=-I*(NM1+I)
            B=B+2._q
            D=1.0_q/(A*D+B)
            C=B+A/C
            DEL=C*D
            H=H*DEL
            IF (ABS(DEL-1._q)< EPS) THEN
               EXPINT=H*EXP(-X)
               RETURN
            ENDIF
         ENDDO
         CALL vtutor%stopCode()
      ELSE
         IF (NM1/=0) THEN
            EXPINT=1.0_q/NM1
         ELSE
            EXPINT=-LOG(X)-EULER
         ENDIF
         FACT=1._q
         DO I=1,MAXIT
            FACT=-FACT*X/I
            IF (I/=NM1) THEN
               DEL=-FACT/(I-NM1)
            ELSE
               PSI=-EULER
               DO II=1,NM1
                  PSI=PSI+1._q/II
               ENDDO
               DEL=FACT*(LOG(X)+PSI)
            ENDIF
            EXPINT=EXPINT+DEL
            IF (ABS(DEL)<ABS(EXPINT)*EPS) RETURN
         ENDDO
         CALL vtutor%stopCode()
      ENDIF
      RETURN
      END FUNCTION EXPINT

      END MODULE wpbe

      SUBROUTINE CALC_EXCHWPBEsol_SP(D,S, &
     &  EXLDA,EXDLDA,EXLDA_SR,EXLDA_LR,EXDLDA_SR,EXDLDA_LR, &
     &  EXWPBE_SR,EXWPBE_LR,EXWPBED_SR,EXWPBED_LR,EXWPBEDD_SR,EXWPBEDD_LR)
!$ACC ROUTINE SEQ
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  INPUT D : DENSITY
!  INPUT S:  (GRAD rho)/(2*KF*rho), where kf=(6 pi^2 rho)^(1/3)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EXWPBE)
!  OUTPUT:  DERIVATIVE W.R.T. DENSITY (EXWPBED)
!  OUTPUT:  DERIVATIVE W.R.T. GRADIENT DENSITY (EXWPBEDD)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      USE prec
      USE constant
      USE wpbe
      IMPLICIT NONE

      REAL(q) F,FS,FSR,FLR
      REAL(q) D,LOGRHO,S,RS,SW

! Numerical derivatives test
      REAL(q) XKF,XKFMIN,XKFPLS,DD
      REAL(q) DDMIN,DDPLS
      REAL(q) SMIN,SPLS,DFDR,DFDDD
      REAL(q) DMIN,DPLS,FSR_MIN,FSR_PLS,XFKMIN,XFKPLS

      REAL(q) EXWPBE,EXWPBE_SR,EXWPBE_LR
      REAL(q) EXWPBED_SR,EXWPBED_LR,EXWPBEDD_SR,EXWPBEDD_LR
      REAL(q) EXLDA,EXDLDA,EXLDA_SR,EXLDA_LR,EXDLDA_SR,EXDLDA_LR
      REAL(q) EXPBE,EXPBED,EXPBEDD

! Parameters to calculate the PBE_SOL enhancement factor
      REAL(q), PARAMETER :: AX=-0.738558766382022405884230032680836_q
      REAL(q), PARAMETER :: UM=0.12345679012345679012_q  !LS
      REAL(q), PARAMETER :: UK=0.8040_q
      REAL(q), PARAMETER :: UL=UM/UK
!$ACC ROUTINE(WPBE_SPLINE) SEQ

! set up the table for F_SR(s)
NOACC CALL INIT_WPBE_TABLE

      SW=ABS(S)
      XKF=(3._q*PI*PI*D)**(1._q/3._q)

! Due to the fact that most of the times only the GGA corrections to the
! LDA exchange are wanted we will have to recalculate the LDA contributions
! in order to be able to subtract them from our PBE result.
! (Georg thinks this makes a lot of sense, but I guess there are not many
! people who would agree with this!)
      RS=(3._q/(4._q*PI)/D)**(1/3._q)
      EXLDA=EX(RS,1,.FALSE.)*D/2._q
      EXLDA_SR=EX_SR(RS,LDASCREEN*AUTOA,1)*D/2._q
      EXLDA_LR=EXLDA-EXLDA_SR
      EXDLDA=VX(RS,1,.FALSE.)/2._q
      EXDLDA_SR=VX_SR(RS,LDASCREEN*AUTOA,1)/2._q
      EXDLDA_LR=VX(RS,1,.FALSE.)/2._q-EXDLDA_SR

! PBE quantities
      F=1._q+UK-UK/(1._q+UL*SW*SW)
      FS=2._q*UM/(1._q+UL*SW*SW)/(1._q+UL*SW*SW)
      EXPBE=F*EXLDA
      EXPBED=4._q/3._q*(F-SW*SW*FS)*EXLDA/D
      EXPBEDD=0.5_q*AX*0.3232409194_q*SW*FS
! set some stuff to zero
      EXWPBE_SR=0
      EXWPBE_LR=0
      EXWPBED_SR=0
      EXWPBED_LR=0
      EXWPBEDD_SR=0
      EXWPBEDD_LR=0
! Interpolate the short range part of the enhancement
! coefficient from the tables
      IF (SW==0._q.OR.D<0._q) RETURN

      LOGRHO=LOG(D)
!    Cutoff criterion to enforce local Lieb-Oxford
!    this ensures that the enhancement factor does not exceed the
!    original one

      DD=SW*(2._q*XKF*D)
!LS
      CALL WPBE_SPLINE(LOGRHO,SW,FSR)
! get values from function directly (slower)
!      CALL EXCHWPBE_R(LDASCREEN*AUTOA,D,SW,FSR)
! Numerical derivatives with interpolated values:
! Firstly: derivatives w.r.t rho
      DMIN=D-0.0001_q*D
      LOGRHO=LOG(DMIN)
      XKFMIN=(3._q*PI*PI*DMIN)**(1._q/3._q)
      SMIN=DD/(DMIN*XKFMIN*2._q)
      CALL WPBE_SPLINE(LOGRHO,SMIN,FSR_MIN)

      DPLS=D+0.0001_q*D
      LOGRHO=LOG(DPLS)
      XKFPLS=(3._q*PI*PI*DPLS)**(1._q/3._q)
      SPLS=DD/(DPLS*XKFPLS*2._q)

      CALL WPBE_SPLINE(LOGRHO,SPLS,FSR_PLS)
      DFDR=(FSR_PLS-FSR_MIN)/2._q/D/0.0001_q
! Secondly: derivatives w.r.t gradient(rho)
      LOGRHO=LOG(D)
      DDMIN=DD-0.0001_q*DD
      DDPLS=DD+0.0001_q*DD
      SMIN=DDMIN/(D*XKF*2._q)
      SPLS=DDPLS/(D*XKF*2._q)

      CALL WPBE_SPLINE(LOGRHO,SMIN,FSR_MIN)
      CALL WPBE_SPLINE(LOGRHO,SPLS,FSR_PLS)
      DFDDD=(FSR_PLS-FSR_MIN)/2._q/DD/0.0001_q

! Find the complementary long range part
      FLR=F-FSR
! E_sr = fsr*Exlda
      EXWPBE_SR=FSR*EXLDA
      EXWPBE_LR=FLR*EXLDA
! dE_sr/drho = dfsr/drho*Exlda + fsr*dExlda/drho
      EXWPBED_SR=DFDR*EXLDA+4._q/3._q*FSR*EXLDA/D
      EXWPBED_LR=EXPBED-EXWPBED_SR
! dE_sr/d(grad rho) = dfsr/dgrad(rho)*Exlda
!gK correct the sign
      EXWPBEDD_SR=DFDDD*EXLDA *SIGN(1.0_q,S)
      EXWPBEDD_LR=(EXPBEDD-DFDDD*EXLDA) *SIGN(1.0_q,S)
!gK end
      RETURN
      END SUBROUTINE CALC_EXCHWPBEsol_SP

      SUBROUTINE CALC_EXCHWPBE_SP(D,S, &
     &  EXLDA,EXDLDA,EXLDA_SR,EXLDA_LR,EXDLDA_SR,EXDLDA_LR, &
     &  EXWPBE_SR,EXWPBE_LR,EXWPBED_SR,EXWPBED_LR,EXWPBEDD_SR,EXWPBEDD_LR)
!$ACC ROUTINE SEQ
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  INPUT D : DENSITY
!  INPUT S:  (GRAD rho)/(2*KF*rho), where kf=(6 pi^2 rho)^(1/3)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EXWPBE)
!  OUTPUT:  DERIVATIVE W.R.T. DENSITY (EXWPBED)
!  OUTPUT:  DERIVATIVE W.R.T. GRADIENT DENSITY (EXWPBEDD)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      USE prec
      USE constant
      USE wpbe
      IMPLICIT NONE

      REAL(q) F,FS,FSR,FLR
      REAL(q) D,LOGRHO,S,RS,SW

! Numerical derivatives test
      REAL(q) XKF,XKFMIN,XKFPLS,DD
      REAL(q) DDMIN,DDPLS
      REAL(q) SMIN,SPLS,DFDR,DFDDD
      REAL(q) DMIN,DPLS,FSR_MIN,FSR_PLS,XFKMIN,XFKPLS

      REAL(q) EXWPBE,EXWPBE_SR,EXWPBE_LR
      REAL(q) EXWPBED_SR,EXWPBED_LR,EXWPBEDD_SR,EXWPBEDD_LR
      REAL(q) EXLDA,EXDLDA,EXLDA_SR,EXLDA_LR,EXDLDA_SR,EXDLDA_LR
      REAL(q) EXPBE,EXPBED,EXPBEDD

! Parameters to calculate the PBE enhancement factor
      REAL(q), PARAMETER :: AX=-0.738558766382022405884230032680836_q
      REAL(q), PARAMETER :: UM=0.2195149727645171_q
      REAL(q), PARAMETER :: UK=0.8040_q
      REAL(q), PARAMETER :: UL=UM/UK
!$ACC ROUTINE(WPBE_SPLINE) SEQ

! set up the table for F_SR(s)
NOACC CALL INIT_WPBE_TABLE

      SW=ABS(S)
      XKF=(3._q*PI*PI*D)**(1._q/3._q)

! Due to the fact that most of the times only the GGA corrections to the
! LDA exchange are wanted we will have to recalculate the LDA contributions
! in order to be able to subtract them from our PBE result.
! (Georg thinks this makes a lot of sense, but I guess there are not many
! people who would agree with this!)
      RS=(3._q/(4._q*PI)/D)**(1/3._q)
      EXLDA=EX(RS,1,.FALSE.)*D/2._q
      EXLDA_SR=EX_SR(RS,LDASCREEN*AUTOA,1)*D/2._q
      EXLDA_LR=EXLDA-EXLDA_SR
      EXDLDA=VX(RS,1,.FALSE.)/2._q
      EXDLDA_SR=VX_SR(RS,LDASCREEN*AUTOA,1)/2._q
      EXDLDA_LR=VX(RS,1,.FALSE.)/2._q-EXDLDA_SR

! PBE quantities
      F=1._q+UK-UK/(1._q+UL*SW*SW)
      FS=2._q*UM/(1._q+UL*SW*SW)/(1._q+UL*SW*SW)
      EXPBE=F*EXLDA
      EXPBED=4._q/3._q*(F-SW*SW*FS)*EXLDA/D
      EXPBEDD=0.5_q*AX*0.3232409194_q*SW*FS
! set some stuff to zero
      EXWPBE_SR=0
      EXWPBE_LR=0
      EXWPBED_SR=0
      EXWPBED_LR=0
      EXWPBEDD_SR=0
      EXWPBEDD_LR=0
! Interpolate the short range part of the enhancement
! coefficient from the tables
      IF (SW==0._q.OR.D<0._q) RETURN

      LOGRHO=LOG(D)
!    Cutoff criterion to enforce local Lieb-Oxford
!    this ensures that the enhancement factor does not exceed the
!    original one

      DD=SW*(2._q*XKF*D)
      CALL WPBE_SPLINE(LOGRHO,SW,FSR)
! get values from function directly (slower)
!      CALL EXCHWPBE_R(LDASCREEN*AUTOA,D,SW,FSR)
! Numerical derivatives with interpolated values:
! Firstly: derivatives w.r.t rho
      DMIN=D-0.0001_q*D
      LOGRHO=LOG(DMIN)
      XKFMIN=(3._q*PI*PI*DMIN)**(1._q/3._q)
      SMIN=DD/(DMIN*XKFMIN*2._q)
      CALL WPBE_SPLINE(LOGRHO,SMIN,FSR_MIN)

      DPLS=D+0.0001_q*D
      LOGRHO=LOG(DPLS)
      XKFPLS=(3._q*PI*PI*DPLS)**(1._q/3._q)
      SPLS=DD/(DPLS*XKFPLS*2._q)

      CALL WPBE_SPLINE(LOGRHO,SPLS,FSR_PLS)
      DFDR=(FSR_PLS-FSR_MIN)/2._q/D/0.0001_q
! Secondly: derivatives w.r.t gradient(rho)
      LOGRHO=LOG(D)
      DDMIN=DD-0.0001_q*DD
      DDPLS=DD+0.0001_q*DD
      SMIN=DDMIN/(D*XKF*2._q)
      SPLS=DDPLS/(D*XKF*2._q)

      CALL WPBE_SPLINE(LOGRHO,SMIN,FSR_MIN)
      CALL WPBE_SPLINE(LOGRHO,SPLS,FSR_PLS)
      DFDDD=(FSR_PLS-FSR_MIN)/2._q/DD/0.0001_q

! Find the complementary long range part
      FLR=F-FSR
! E_sr = fsr*Exlda
      EXWPBE_SR=FSR*EXLDA
      EXWPBE_LR=FLR*EXLDA
! dE_sr/drho = dfsr/drho*Exlda + fsr*dExlda/drho
      EXWPBED_SR=DFDR*EXLDA+4._q/3._q*FSR*EXLDA/D
      EXWPBED_LR=EXPBED-EXWPBED_SR
! dE_sr/d(grad rho) = dfsr/dgrad(rho)*Exlda
!gK correct the sign
      EXWPBEDD_SR=DFDDD*EXLDA *SIGN(1.0_q,S)
      EXWPBEDD_LR=(EXPBEDD-DFDDD*EXLDA) *SIGN(1.0_q,S)
!gK end
      RETURN
      END SUBROUTINE CALC_EXCHWPBE_SP

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE EX_SR_SIMD(RS, RMU, IFLG, EX_SR, MASK0)
      USE prec
      USE constant
      USE simd
      IMPLICIT NONE
! Exchange energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      TYPE(SIMD_REAL64), INTENT(IN)      :: RS
      REAL(q), INTENT(IN)                :: RMU
      INTEGER, INTENT(IN)                :: IFLG
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EX_SR
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      REAL(q)                            :: CX, EX_, RHO, QFAC, QF, A, FRAC, Z, T, ERRF_RESULT
      REAL(q)                            :: TEMP_1, TEMP_2, TEMP_3
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: CX1 = 0.9163305865663_q
      REAL(q), PARAMETER                 :: CX2 = 1.1545041946774_q
      REAL(q), PARAMETER                 :: THRD = 0.333333333333333333333333_q    ! 1._q/3._q
      REAL(q), PARAMETER                 :: QFAC1 = 3.093667726280135532590521_q   ! (3._q*PI*PI)**(1._q/3._q)
      REAL(q), PARAMETER                 :: QFAC2 = 3.897777089720753540547094_q   ! (6._q*PI*PI)**(1._q/3._q)
      REAL(q), PARAMETER                 :: SQRT_PI = 1.772453850905515881919428_q ! sqrt(pi)
      REAL(q), PARAMETER                 :: P1 = 0.620350490899400086597382_q      ! (3._q/(4._q*pi))**(1._q/3._q)
      REAL(q), PARAMETER                 :: P_2 = 2.666666666666666518636930_q     ! 8._q/3._q
      REAL(q), PARAMETER                 :: P_3 = 0.027777777777777776235801_q     ! 1._q/36._q

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(QFAC, CX, TEMP_1, EX_, QF, A, TEMP_3, Z, T, ERRF_RESULT, TEMP_2, FRAC)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (IFLG .EQ. 1) THEN
            QFAC = QFAC1
            CX = CX1
         ELSE
            QFAC = QFAC2
            CX = CX2
         ENDIF
         TEMP_1 = 1.0_q / RS%X(II)
         EX_ = -CX * TEMP_1
! This version is faster and potentially more precise
         QF = QFAC * P1 * TEMP_1
         A = RMU / (2.0_q * QF)
         TEMP_3 = 0.0_q
! Test on the value of A
         IF (A .LT. 1.E-9_q) THEN
         ! Limit for small A
            TEMP_3 = EX_
         ENDIF
         TEMP_1 = 1.0_q / A
         IF ((A .GE. 1.0E-9_q) .AND. (A .LE. 100.0_q)) THEN
         ! Inline ERRF : potentially faster
            Z = 0.5_q * TEMP_1
            T = 1.0_q / (1.0_q + 0.5_q * Z)
            ERRF_RESULT = -Z * Z - 1.26551223_q + T * (1.00002368_q + T * (.37409196_q +  &
            &  T * (0.09678418_q + T * (-0.18628806_q + T * (0.27886807_q + T * (-1.13520398_q +  &
            &  T * (1.48851587_q + T * (-0.82215223_q + T * 0.17087277_q))))))))
            ERRF_RESULT = 1.0_q - T * SIMD_EXP(ERRF_RESULT)
         ! End of Inline ERRF
            TEMP_2 = 4.0_q * A * A *A
            FRAC = (1.0_q - P_2 * A * (SQRT_PI * ERRF_RESULT +  &
            &  (2.0_q * A - TEMP_2) * SIMD_EXP(-0.25_q * TEMP_1 * TEMP_1) - 3.0_q * A + TEMP_2))
            TEMP_3 = EX_ * FRAC
         ENDIF
         IF ((A .GT. 100.0_q) .AND. (A .LE. 1.0E+9_q)) THEN
         ! Development for large A
            FRAC = P_3 * TEMP_1
            TEMP_3 = EX_ * FRAC
         ENDIF
         IF (MASK0%X(II)) THEN
            EX_SR%X(II) = TEMP_3
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE EX_SR_SIMD


      SUBROUTINE VX_SR_SIMD(RS, RMU, IFLG, VX_SR, MASK0)
      USE prec
      USE constant
      USE simd
      IMPLICIT NONE
! Exchange energy:
!     IFLG=1:  Paramagnetic results
!     IFLG=2:  Ferromagnetic results
      TYPE(SIMD_REAL64), INTENT(IN)      :: RS
      REAL(q), INTENT(IN)                :: RMU
      INTEGER, INTENT(IN)                :: IFLG
      TYPE(SIMD_REAL64), INTENT(OUT)     :: VX_SR
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      REAL(q)                            :: CX, VX_, RHO, QFAC, QF, A, FRAC, Z, T, ERRF_RESULT, EX_SR_RESULT
      REAL(q)                            :: TEMP_1, TEMP_2, TEMP_3, TEMP_4, TEMP_5
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: CX1 = 0.9163305865663_q
      REAL(q), PARAMETER                 :: CX2 = 1.1545041946774_q
      REAL(q), PARAMETER                 :: THRD = 0.333333333333333333333333_q    ! 1._q/3._q
      REAL(q), PARAMETER                 :: QFAC1 = 3.093667726280135532590521_q   ! (3._q*PI*PI)**(1._q/3._q)
      REAL(q), PARAMETER                 :: QFAC2 = 3.897777089720753540547094_q   ! (6._q*PI*PI)**(1._q/3._q)
      REAL(q), PARAMETER                 :: SQRT_PI = 1.772453850905515881919428_q ! sqrt(pi)
      REAL(q), PARAMETER                 :: P_1 = 0.620350490899400086597382_q     ! (1._q/(4._q*pi))**(1._q/3._q)
      REAL(q), PARAMETER                 :: P_2 = 1.333333333333333259318465_q     ! 4._q/3._q
      REAL(q), PARAMETER                 :: P_3 = 2.666666666666666518636930_q     ! 8._q/3._q
      REAL(q), PARAMETER                 :: P_4 = 0.006631455962162306067009_q     ! 1._q/(48._q*pi)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(QFAC, CX, TEMP_1, VX_, QF, A, TEMP_5, Z, T, TEMP_2, TEMP_3, TEMP_4, FRAC, EX_SR_RESULT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         IF (IFLG .EQ. 1) THEN
            QFAC = QFAC1
            CX = CX1
         ELSE
            QFAC = QFAC2
            CX = CX2
         ENDIF
         TEMP_1 = 1.0_q / RS%X(II)
         VX_ = -CX * TEMP_1
! This version is faster and potentially more precise
         QF = QFAC * P_1 * TEMP_1
         A = RMU / (2.0_q * QF)
         TEMP_5 = 0.0_q
! Test on the value of A
         IF (A .LT. 1.0E-9_q) THEN
         ! Limit for small A
            TEMP_5 = P_2 * VX_
         ENDIF
         TEMP_1 = 1.0_q / A
         IF ((A .GE. 1.0E-9_q) .AND. (A .LE. 100.0_q)) THEN
         ! Intermediate Values of A
         ! Inline EX_SR : computations already carried out can be reused
         ! Inline ERRF
            Z = 0.5_q * TEMP_1
            T = 1.0_q / (1.0_q + 0.5_q * Z)
            ERRF_RESULT = -Z * Z - 1.26551223_q + T * (1.00002368_q + T * (.37409196_q +  &
            &  T * (0.09678418_q + T * (-0.18628806_q + T * (0.27886807_q + T * (-1.13520398_q +  &
            &  T * (1.48851587_q + T * (-0.82215223_q + T * 0.17087277_q))))))))
            ERRF_RESULT = 1.0_q - T * SIMD_EXP(ERRF_RESULT)
         ! End of Inline ERRF
            TEMP_2 = 4.0_q * A * A * A
            TEMP_3 = SIMD_EXP(-0.25_q * TEMP_1 * TEMP_1)
            FRAC = (1.0_q - P_3 * A * (SQRT_PI * ERRF_RESULT +  &
            &  (2.0_q * A - TEMP_2) * TEMP_3 - 3.0_q * A + TEMP_2))
            EX_SR_RESULT = VX_ * FRAC
         ! End of Inline EX_SR
            TEMP_2 = 4.0_q * A * A
            TEMP_4 = 2.0_q * TEMP_2
            FRAC = (1.0_q - TEMP_2 * TEMP_4 * (TEMP_3 - 1.0_q) - TEMP_4)
            TEMP_5 = EX_SR_RESULT + VX_ * FRAC * THRD
         ENDIF
         IF ((A .GT. 100.0_q) .AND. (A .LE. 1.0E+9_q)) THEN
         ! Development for large A
            TEMP_5 = -RMU * P_4 * TEMP_1 * TEMP_1 * TEMP_1
         ENDIF
         IF (MASK0%X(II)) THEN
            VX_SR%X(II) = TEMP_5
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE VX_SR_SIMD


      SUBROUTINE CALC_EXCHWPBE_SP_SIMD(D, S,  &
      &  EXLDA, EXDLDA, EXLDA_SR, EXLDA_LR, EXDLDA_SR, EXDLDA_LR,  &
      &  EXWPBE_SR, EXWPBE_LR, EXWPBED_SR, EXWPBED_LR, EXWPBEDD_SR, EXWPBEDD_LR, MASK0)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  INPUT D : DENSITY
!  INPUT S:  (GRAD rho)/(2*KF*rho), where kf=(6 pi^2 rho)^(1/3)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EXWPBE)
!  OUTPUT:  DERIVATIVE W.R.T. DENSITY (EXWPBED)
!  OUTPUT:  DERIVATIVE W.R.T. GRADIENT DENSITY (EXWPBEDD)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      USE prec
      USE constant
      USE wpbe
      USE simd
      IMPLICIT NONE
      TYPE(SIMD_REAL64), INTENT(IN)      :: D, S
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EXLDA, EXDLDA, EXLDA_SR, EXLDA_LR, EXDLDA_SR, EXDLDA_LR
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EXWPBE_SR, EXWPBE_LR, EXWPBED_SR, EXWPBED_LR, EXWPBEDD_SR, EXWPBEDD_LR
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
      TYPE(SIMD_REAL64)                  :: RS, SW, LOGRHO1, LOGRHO2, LOGRHO3, SMIN, SPLS, XKF
      TYPE(SIMD_REAL64)                  :: F, FSR, FSR_MIN, FSR_PLS, EXPBED, EXPBEDD, DD, DFDR, DFDDD, INV_D
      TYPE(SIMD_MASK_REAL64)             :: MASK1
      REAL(q)                            :: FS, FLR
! Numerical derivatives test
      REAL(q)                            :: XKFMIN, XKFPLS
      REAL(q)                            :: DDMIN, DDPLS
      REAL(q)                            :: DMIN, DPLS, XFKMIN, XFKPLS
      REAL(q)                            :: EXWPBE, EXPBE
      REAL(q)                            :: TEMP_1
      INTEGER                            :: II
      LOGICAL                            :: TRUE_FOR_ANY
! Parameters to calculate the PBE enhancement factor
      REAL(q), PARAMETER                 :: THRD = 0.333333333333333333333333_q
      REAL(q), PARAMETER                 :: AX = -0.738558766382022405884230032680836_q
      REAL(q), PARAMETER                 :: UM = 0.2195149727645171_q
      REAL(q), PARAMETER                 :: UMSOL = 0.12345679012345679012_q
      REAL(q), PARAMETER                 :: UK = 0.8040_q
      REAL(q), PARAMETER                 :: UL = UM / UK
      REAL(q), PARAMETER                 :: ULSOL = UMSOL / UK
      REAL(q), PARAMETER                 :: P_1 = 0.238732414637843004534545_q ! 3._q/(4._q*pi)
      REAL(q), PARAMETER                 :: P_2 = 1.333333333333333259318465_q ! 4._q/3._q

      CALL INIT_WPBE_TABLE

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT)
      DO II = 0, SIMD_WIDTH_64BIT - 1
! Due to the fact that most of the times only the GGA corrections to the
! LDA exchange are wanted we will have to recalculate the LDA contributions
! in order to be able to subtract them from our PBE result.
! (Georg thinks this makes a lot of sense, but I guess there are not many
! people who would agree with this!)
         RS%X(II) = SIMD_POW(P_1 / D%X(II), THRD)
! continued below
      ENDDO

      CALL EX_SR_SIMD(RS, LDASCREEN * AUTOA, 1, EXLDA_SR, MASK0)

      CALL VX_SR_SIMD(RS, LDASCREEN * AUTOA, 1, EXDLDA_SR, MASK0)

      TRUE_FOR_ANY = .FALSE.
!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) REDUCTION(.OR. : TRUE_FOR_ANY) PRIVATE(TEMP_1)
      DO II = 0, SIMD_WIDTH_64BIT - 1
         SW%X(II) = ABS(S%X(II))
         MASK1%X(II) = MASK0%X(II)
         IF (MASK1%X(II)) THEN
            TEMP_1 = 1.0_q / RS%X(II)
            EXLDA%X(II) = (-0.9163305865663_q * TEMP_1) * D%X(II) * 0.5_q
            EXDLDA%X(II) = (-1.2217741154217_q * TEMP_1) * 0.5_q
            EXLDA_SR%X(II) = EXLDA_SR%X(II) * D%X(II) * 0.5_q
            EXLDA_LR%X(II) = EXLDA%X(II) - EXLDA_SR%X(II)
            EXDLDA_SR%X(II) = EXDLDA_SR%X(II) * 0.5_q
            EXDLDA_LR%X(II) = EXDLDA%X(II) - EXDLDA_SR%X(II)
            EXWPBE_SR%X(II) = 0.0_q
            EXWPBE_LR%X(II) = 0.0_q
            EXWPBED_SR%X(II) = 0.0_q
            EXWPBED_LR%X(II) = 0.0_q
            EXWPBEDD_SR%X(II) = 0.0_q
            EXWPBEDD_LR%X(II) = 0.0_q
! Interpolate the short range part of the enhancement
! coefficient from the tables            
            IF ((SW%X(II) .NE. 0.0_q) .AND. (D%X(II) .GE. 0.0_q)) THEN
               TRUE_FOR_ANY = .TRUE.
            ELSE
               MASK1%X(II) = .FALSE.
            ENDIF
         ENDIF
      ENDDO

      IF (.NOT. TRUE_FOR_ANY) RETURN

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1, FS, EXPBE, DMIN, XKFMIN, DPLS, XKFPLS)
      DO II = 0, SIMD_WIDTH_64BIT - 1
! PBE quantities
         IF (LEXCH .EQ. 14) THEN
            TEMP_1 = 1.0_q / (1.0_q + ULSOL * SW%X(II) * SW%X(II))
            FS = 2.0_q * UMSOL * TEMP_1 * TEMP_1
         ELSE
            TEMP_1 = 1.0_q / (1.0_q + UL * SW%X(II) * SW%X(II))
            FS = 2.0_q * UM * TEMP_1 * TEMP_1
         ENDIF
         F%X(II) = 1.0_q + UK - UK * TEMP_1
         EXPBE = F%X(II) * EXLDA%X(II)
         INV_D%X(II) = 1.0_q / D%X(II)
         EXPBED%X(II) = P_2 * (F%X(II) - SW%X(II) * SW%X(II) * FS) * EXLDA%X(II) * INV_D%X(II)
         EXPBEDD%X(II) = 0.5_q * AX * 0.3232409194_q * SW%X(II) * FS
! Cutoff criterion to enforce local Lieb-Oxford
! this ensures that the enhancement factor does not exceed the
! original one
! get values from function directly (slower)
!     CALL EXCHWPBE_R(LDASCREEN*AUTOA,D,SW,FSR)
! Numerical derivatives with interpolated values:
! Firstly: derivatives w.r.t rho
         XKF%X(II) = SIMD_POW(3.0_q * PI * PI * D%X(II), THRD)
         DD%X(II) = SW%X(II) * (2.0_q * XKF%X(II) * D%X(II))
         DMIN = 0.9999_q * D%X(II)
         LOGRHO2%X(II) = SIMD_LOG(DMIN)
         XKFMIN = SIMD_POW(3.0_q * PI * PI * DMIN, THRD)
         SMIN%X(II) = DD%X(II) / (DMIN * XKFMIN * 2.0_q)
         DPLS = 1.0001_q * D%X(II)
         LOGRHO3%X(II) = SIMD_LOG(DPLS)
         XKFPLS = SIMD_POW(3.0_q * PI * PI * DPLS, THRD)
         SPLS%X(II) = DD%X(II) / (DPLS * XKFPLS * 2.0_q)
         LOGRHO1%X(II) = SIMD_LOG(D%X(II))
      ENDDO

      CALL WPBE_SPLINE_SIMD(LOGRHO1, SW, FSR, MASK1)

      CALL WPBE_SPLINE_SIMD(LOGRHO2, SMIN, FSR_MIN, MASK1)

      CALL WPBE_SPLINE_SIMD(LOGRHO3, SPLS, FSR_PLS, MASK1)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(FLR, DDMIN, DDPLS, TEMP_1)
      DO II=0,SIMD_WIDTH_64BIT-1
         DFDR%X(II) = (FSR_PLS%X(II) - FSR_MIN%X(II)) * 5000.0_q * INV_D%X(II)
! Find the complementary long range part
         FLR = F%X(II) - FSR%X(II)
! E_sr = fsr*Exlda
         IF (MASK1%X(II)) THEN
            EXWPBE_SR%X(II) = FSR%X(II) * EXLDA%X(II)
            EXWPBE_LR%X(II) = FLR * EXLDA%X(II)
         ENDIF
! Secondly: derivatives w.r.t gradient(rho)
         DDMIN = 0.9999_q * DD%X(II)
         DDPLS = 1.0001_q * DD%X(II)
         TEMP_1 = 0.5_q / (D%X(II) * XKF%X(II))
         SMIN%X(II) = DDMIN * TEMP_1
         SPLS%X(II) = DDPLS * TEMP_1
      ENDDO

      CALL WPBE_SPLINE_SIMD(LOGRHO1, SMIN, FSR_MIN, MASK1)

      CALL WPBE_SPLINE_SIMD(LOGRHO1, SPLS, FSR_PLS, MASK1)

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(TEMP_1)
      DO II=0,SIMD_WIDTH_64BIT-1
         DFDDD%X(II) = (FSR_PLS%X(II) - FSR_MIN%X(II)) * 5000.0_q / DD%X(II)
! dE_sr/drho = dfsr/drho*Exlda + fsr*dExlda/drho
         IF (MASK1%X(II)) THEN
            EXWPBED_SR%X(II) = DFDR%X(II) * EXLDA%X(II) + P_2 * FSR%X(II)*EXLDA%X(II) * INV_D%X(II)
            EXWPBED_LR%X(II) = EXPBED%X(II) - EXWPBED_SR%X(II)
! dE_sr/d(grad rho) = dfsr/dgrad(rho)*Exlda
         ENDIF
!gK correct the sign
         IF (S%X(II) .LT. 0.0_q) THEN
            TEMP_1 = -1.0_q
         ELSE
            TEMP_1 = 1.0_q
         ENDIF
         IF (MASK1%X(II)) THEN
            EXWPBEDD_SR%X(II) = DFDDD%X(II) * EXLDA%X(II) * TEMP_1
            EXWPBEDD_LR%X(II) = (EXPBEDD%X(II) - DFDDD%X(II) * EXLDA%X(II)) * TEMP_1
         ENDIF
      ENDDO
!gK end
      RETURN
      END SUBROUTINE CALC_EXCHWPBE_SP_SIMD
#endif

!----------------------------------------------------------------------
! jP: adding PBEsol
!     Note that this is a plain copy-paste - just the necessary constants are changed!
!     The quickest way to program
!     PBEsol is to take your existing PBE subroutine, and alter the two values
!     mu and beta:  In PBE, mu=0.219 and beta=0.067 these become mu=10/81
!     and beta=0.046 in PBEsol.   The routines below are simple modifications
!     of the PBE subroutines, to produce PBEsol.
!---------------------------------------------------------------------
! The P B E
! this routine was kindly supplied by Bjork Hammer
! it is essentially identical to Burkes routine but includes
! the revised P B E routines
!---------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
SUBROUTINE EXCHPBESOL(rho,rhothrd,s,exlda,expbe,exdlda,exd,exdd, &
     &                   ukfactor)
!$ACC ROUTINE SEQ
  !----------------------------------------------------------------------
  !  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
  !  K Burkes modification of PW91 codes, May 14, 1996
  !  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  !  INPUT rho : DENSITY
  !  INPUT rhothrd : DENSITY^(1/3)
  !  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
  !  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX)
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! References:
  ! [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
  ! [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
  !     {\bf 40},  3399  (1989) (E).
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! Formulas:
  !       e_x[unif]=ax*rho^(4/3)  [LDA]
  ! ax = -0.75*(3/pi)^(1/3)
  !       e_x[PBE]=e_x[unif]*FxPBE(s)
  !       FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
  ! uk, ul defined after [a](13)
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  USE prec
  IMPLICIT REAL(q) (A-H,O-Z)
  parameter(thrd=1._q/3._q,thrd4=4._q/3._q)
  parameter(pi=3.14159265358979323846264338327950_q)
  parameter(ax=-0.738558766382022405884230032680836_q)
  parameter(um=0.12345679012345679012_q,uk1=0.8040_q,ul1=um/uk1)

  !----------------------------------------------------------------------
  ! construct LDA exchange energy density
  exunif=AX*rhothrd
  exlda=exunif*rho
  exdlda=exunif*thrd4
  !----------------------------------------------------------------------
  ! construct PBE enhancement factor
  S2 = S*S
  !----------------------------------------------------------------------
  if (ukfactor.ne.0.0_q) then
     ! These are the PBE96 and revPBE98 functionals
     ! scale uk with a factor
     uk = uk1*ukfactor
     ul = ul1/ukfactor
     P0=1._q+ul*S2
     FxPBE = 1._q+uk-uk/P0
     expbe = exlda*FxPBE
     !----------------------------------------------------------------------
     !  ENERGY DONE. NOW THE POTENTIAL:
     !  find first derivatives of Fx w.r.t s.
     !  Fs=(1/s)*d FxPBE/ ds
     Fs=2._q*um/(P0*P0)
  else
     ! This is the RPBE functional [Hammer et al, PRB 59, 7413 (1999)]
     P0=exp(-ul1*S2)
     FxPBE = 1._q+uk1*(1.0_q-P0)
     expbe = exlda*FxPBE
     !----------------------------------------------------------------------
     !  ENERGY DONE. NOW THE POTENTIAL:
     !  find first derivatives of Fx w.r.t s.
     !  Fs=(1/s)*d FxPBE/ ds
     Fs=2._q*um*P0
  endif
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! calculate the partial derivatives of ex wrt n and |grad(n)|
  !  0.3232409194=(3*pi^2)^(-1/3)
  exd =exunif*THRD4*(FxPBE-S2*Fs)
  exdd=0.5_q*ax*0.3232409194_q*S*Fs
  RETURN
END SUBROUTINE EXCHPBESOL

#if defined(SIMD_OPTIMIZED) && defined(SIMD_TYPES)
      SUBROUTINE EXCHPBESOL_SIMD(RHO, RHOTHRD, S, EXLDA, EXPBE, EXDLDA, EXD, EXDD, UKFACTOR, MASK0)
  !----------------------------------------------------------------------
  !  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
  !  K Burkes modification of PW91 codes, May 14, 1996
  !  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  !  INPUT rho : DENSITY
  !  INPUT rhothrd : DENSITY^(1/3)
  !  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
  !  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX)
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! References:
  ! [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
  ! [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
  !     {\bf 40},  3399  (1989) (E).
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! Formulas:
  !       e_x[unif]=ax*rho^(4/3)  [LDA]
  ! ax = -0.75*(3/pi)^(1/3)
  !       e_x[PBE]=e_x[unif]*FxPBE(s)
  !       FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
  ! uk, ul defined after [a](13)
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
      USE prec
      USE simd
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE(SIMD_REAL64), INTENT(IN)      :: RHO, RHOTHRD, S, UKFACTOR
      TYPE(SIMD_REAL64), INTENT(OUT)     :: EXLDA, EXPBE, EXDLDA, EXD, EXDD
      TYPE(SIMD_MASK_REAL64), INTENT(IN) :: MASK0
!      TYPE(SIMD_REAL64)                  :: S2, FXPBE, EXUNIF, FS
      INTEGER                            :: II
      REAL(q), PARAMETER                 :: THRD = 1.0_q / 3.0_q, THRD4 = 4.0_q / 3.0_q
      REAL(q), PARAMETER                 :: PI = 3.14159265358979323846264338327950_q
      REAL(q), PARAMETER                 :: AX = -0.738558766382022405884230032680836_q
      REAL(q), PARAMETER                 :: UM = 0.12345679012345679012_q, UK1 = 0.8040_q, UL1 = UM / UK1
      REAL(q), PARAMETER                 :: P_1 = 0.5_q * AX * 0.3232409194_q

!$OMP SIMD OMP_SIMDLEN(SIMD_WIDTH_64BIT) PRIVATE(EXUNIF, S2, UK, UL, P0, FXPBE, FS)
      DO II = 0, SIMD_WIDTH_64BIT - 1
!----------------------------------------------------------------------
! construct LDA exchange energy density
         EXUNIF = AX * RHOTHRD%X(II)
         IF (MASK0%X(II)) THEN
            EXLDA%X(II) = EXUNIF * RHO%X(II)
            EXDLDA%X(II) = EXUNIF * THRD4
         ENDIF
!----------------------------------------------------------------------
! construct PBE enhancement factor
         S2 = S%X(II) * S%X(II)
!----------------------------------------------------------------------
         IF (UKFACTOR%X(II) .NE. 0.0_q) THEN
! These are the PBE96 and revPBE98 functionals
! scale uk with a factor
            UK = UK1 * UKFACTOR%X(II)
            UL = UL1 / UKFACTOR%X(II)
            P0 = 1.0_q / (1.0_q + UL * S2)
            FXPBE = 1.0_q + UK - UK * P0
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
            FS = 2.0_q * UM * P0 * P0
         ELSE
     ! This is the RPBE functional [Hammer et al, PRB 59, 7413 (1999)]
            P0 = SIMD_EXP(-UL1 * S2)
            FXPBE = 1.0_q + UK1 * (1.0_q - P0)
     !----------------------------------------------------------------------
     !  ENERGY DONE. NOW THE POTENTIAL:
     !  FIND FIRST DERIVATIVES OF FX W.R.T S.
     !  FS=(1/S)*D FXPBE/ DS
            FS = 2.0_q * UM * P0
         ENDIF
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! calculate the partial derivatives of ex wrt n and |grad(n)|
  !  0.3232409194=(3*pi^2)^(-1/3)
         IF (MASK0%X(II)) THEN
            EXPBE%X(II) = EXLDA%X(II) * FXPBE
            EXD%X(II) = EXUNIF * THRD4 * (FXPBE - S2 * FS)
            EXDD%X(II) = P_1 * S%X(II) * FS
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE EXCHPBESOL_SIMD
#endif
!----------------------------------------------------------------------
!----------------------------------------------------------------------

! jP ######################################################################
!----------------------------------------------------------------------
SUBROUTINE CORunspPBESOL(RS,EC,VC,sk, &
     &                  T,H,DVC,ecdd,lgga)
!$ACC ROUTINE SEQ
  !----------------------------------------------------------------------
  !  Official PBE correlation code. K. Burke, May 14, 1996.
  !  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
  !       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
  !       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
  !       : lgga=flag to do gga (0=>LSD only)
  !       : lpot=flag to do potential (0=>energy only)
  !  output: ec=lsd correlation energy from [a]
  !        : vcup=lsd up correlation potential
  !        : vcdn=lsd dn correlation potential
  !        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
  !        : dvcup=nonlocal correction to vcup
  !        : dvcdn=nonlocal correction to vcdn
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
  parameter(bet=0.046_q,delt=bet/gamma)
  parameter(eta=1.e-12_q)
!$ACC ROUTINE (GCOR2) SEQ
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
  CALL gcor2(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
       &    0.49294_q,rtrs,EU,EURS)
  EC = EU
  ! check for zero energy, immediate return if true
  IF (EC==0.) THEN
     H=0; DVC=0; ecdd=0
     RETURN
  ENDIF
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! LSD potential from [c](A1)
  ! ECRS = dEc/drs [c](A2)
  ! ECZET=dEc/dzeta [c](A3)
  ! FZ = dF/dzeta [c](A4)
  ECRS = EURS
  VC = EC -RS*ECRS/3._q
  if (.not.lgga) return
  !----------------------------------------------------------------------
  ! PBE correlation energy
  ! G=phi(zeta), given after [a](3)
  ! DELT=bet/gamma
  ! B=A of [a](8)
  PON=-EC/(gamma)
  B = DELT/(DEXP(PON)-1._q)
  B2 = B*B
  T2 = T*T
  T4 = T2*T2
  RS2 = RS*RS
  RS3 = RS2*RS
  Q4 = 1._q+B*T2
  Q5 = 1._q+B*T2+B2*T4
  H = (BET/DELT)*DLOG(1._q+DELT*Q4*T2/Q5)
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
  T6 = T4*T2
  RSTHRD = RS/3._q
  FAC = DELT/B+1._q
  BEC = B2*FAC/(BET)
  Q8 = Q5*Q5+DELT*Q4*Q5*T2
  Q9 = 1._q+2._q*B*T2
  hB = -BET*B*T6*(2._q+B*T2)/Q8
  hRS = -RSTHRD*hB*BEC*ECRS
  hT = 2._q*BET*Q9/Q8
  DVC = H+HRS-7.0_q*T2*HT/6._q
  ecdd=0.5_q/sk*t*ht
  RETURN
END SUBROUTINE CORunspPBESOL
!----------------------------------------------------------------------


!######################################################################
! jP ----------------------------------------------------------------------
SUBROUTINE CORPBESOL(RS,ZET,EC,VCUP,VCDN,g,sk, &
     &                  T,H,DVCUP,DVCDN,ecdd,lgga)
!$ACC ROUTINE SEQ
  !----------------------------------------------------------------------
  !  Official PBE correlation code. K. Burke, May 14, 1996.
  !  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
  !       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
  !       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
  !       : lgga=flag to do gga (0=>LSD only)
  !       : lpot=flag to do potential (0=>energy only)
  !  output: ec=lsd correlation energy from [a]
  !        : vcup=lsd up correlation potential
  !        : vcdn=lsd dn correlation potential
  !        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
  !        : dvcup=nonlocal correction to vcup
  !        : dvcdn=nonlocal correction to vcdn
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
  parameter(bet=0.046_q,delt=bet/gamma)
  parameter(eta=1.e-12_q)
!$ACC ROUTINE(GCOR2) SEQ
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
  CALL gcor2(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
       &    0.49294_q,rtrs,EU,EURS)
  CALL gcor2(0.01554535_q,0.20548_q,14.1189_q,6.1977_q,3.3662_q, &
       &    0.62517_q,rtRS,EP,EPRS)
  CALL gcor2(0.0168869_q,0.11125_q,10.357_q,3.6231_q,0.88026_q, &
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
  if (.not.lgga) return
  !----------------------------------------------------------------------
  ! PBE correlation energy
  ! G=phi(zeta), given after [a](3)
  ! DELT=bet/gamma
  ! B=A of [a](8)
  G3 = G**3
  PON=-EC/(G3*gamma)
  B = DELT/(DEXP(PON)-1._q)
  B2 = B*B
  T2 = T*T
  T4 = T2*T2
  RS2 = RS*RS
  RS3 = RS2*RS
  Q4 = 1._q+B*T2
  Q5 = 1._q+B*T2+B2*T4
  H = G3*(BET/DELT)*DLOG(1._q+DELT*Q4*T2/Q5)
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  ! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
  G4 = G3*G
  T6 = T4*T2
  RSTHRD = RS/3._q
  GZ=(((1._q+zet)**2+eta)**sixthm- &
       &((1._q-zet)**2+eta)**sixthm)/3._q
  FAC = DELT/B+1._q
  BG = -3._q*B2*EC*FAC/(BET*G4)
  BEC = B2*FAC/(BET*G3)
  Q8 = Q5*Q5+DELT*Q4*Q5*T2
  Q9 = 1._q+2._q*B*T2
  hB = -BET*G3*B*T6*(2._q+B*T2)/Q8
  hRS = -RSTHRD*hB*BEC*ECRS
  hZ = 3._q*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
  hT = 2._q*BET*G3*Q9/Q8
  COMM = H+HRS-7.0_q*T2*HT/6._q
  PREF = HZ-GZ*T2*HT/G
  COMM = COMM-PREF*ZET
  DVCUP = COMM + PREF
  DVCDN = COMM - PREF
  ecdd=0.5_q/(sk*g)*t*ht
  RETURN
END SUBROUTINE CORPBESOL
! jP: adding PBEsol
