#include "symbol.inc"
      MODULE setexm_struct_def
      USE prec
#ifdef USELIBXC
      USE xc_f03_lib_m
#endif
      IMPLICIT NONE
      !
      ! the table has two section one for low densities
      !   (1...NEXCHF(1),0...RHOEXC(1))
      ! and one for high densities
      !   (...NEXCHF(2),...RHOEXC(2))
      !
      INTEGER, PARAMETER :: NEXCH=4000

      TYPE exctable
        REAL(q) :: EXCTAB(NEXCH,5,6) ! table including spline coeff.
        REAL(q) :: RHOEXC(2)         ! maximal densities
        INTEGER :: NEXCHF(2)         ! number of points
        INTEGER :: LEXCH             ! type of exchange correlation
      END TYPE
      !
      ! string read from INCAR for GGA=XX entry
      !
      CHARACTER(LEN=2), SAVE :: SZGGA
      !
      ! LEXCH specifies the exchange correlation type VASP applies
      ! throughout the calculations
      ! it either corresponds to the GGA entry
      ! or is defaulted from the POTCAR files
      !
      INTEGER, SAVE :: LEXCH=-1
      !
      ! stores the exchange correlation interpolation table
      !
      TYPE (exctable), SAVE :: EXCTAB
      !
      ! interpolation of correlation from paramagnetic to 
      ! ferromagnetic case according to 
      ! Vosko, Wilk and Nusair, CAN. J. PHYS. 58, 1200 (1980)
      !
      INTEGER, SAVE :: LFCI=1
      !
      ! amount of LDA exchange
      ! this can be reduced e.g. if exact exchange is used
      ! (e.g. hybrid functionals)
      !
      REAL(q), SAVE :: LDAX=1
      REAL(q), SAVE :: ALDAC=1
      REAL(q), SAVE :: AGGAX=1
      REAL(q), SAVE :: AGGAC=1
      !
      ! screened LDA exchange parameter
      !
      REAL(q), SAVE :: LDASCREEN=0
      REAL(q), SAVE :: LDASCREENC=0
      !
      ! LRANGE_SEPARATED_CORR=.TRUE.  range separated LDA correlation
      ! LRANGE_SEPARATED_CORR=.FALSE. complete LDA correlation (default)
      !
      LOGICAL, SAVE :: LRANGE_SEPARATED_CORR=.FALSE.
      !
      ! LUSE_LONGRANGE_HF  short range LDA exchange interaction only
      ! long range contribution is done in HF
      ! the default is that HF treats short range, and LDA long range
      ! but using this flag the behavior can inverted
      ! the flag should be identical to LRHFCALC in fock.F
      !
      LOGICAL, SAVE :: LUSE_LONGRANGE_HF=.FALSE.
      !
      ! LUSE_THOMAS_FERMI Thomas Fermi screening in local exchange
      ! should be identical to L_THOMAS_FERMI in fock.F
      !
      LOGICAL, SAVE :: LUSE_THOMAS_FERMI=.FALSE.
      !
      ! LUSE_MODEL_HF no local exchange in the short range
      ! a fraction LDAX of the local exchange in the long range limit
      ! should be identical to L_MODEL_HF
      !
      LOGICAL, SAVE :: LUSE_MODEL_HF=.FALSE.
      !
      ! for vdW-DF functionals
      !
      REAL(q), SAVE :: PARAM1, PARAM2, PARAM3
      REAL(q), SAVE :: Zab_VDW
      LOGICAL, SAVE :: LUSE_VDW
      !
      ! rVV10
      !
      REAL(q), SAVE :: BPARAM
      REAL(q), SAVE :: CPARAM
      !
      ! libxc related variables
      !
      CHARACTER(LEN=35) :: LIBXC_VERSION='     F    Libxc'
#ifdef USELIBXC
      INTEGER, PARAMETER :: &
     &   LIBXC_FAMILY_UNKNOWN = -1, &
     &   LIBXC_FAMILY_NONE = 0, &
     &   LIBXC_FAMILY_LDA = 1, &
     &   LIBXC_FAMILY_GGA = 2, &
     &   LIBXC_FAMILY_MGGA = 4, &
     &   LIBXC_FAMILY_LCA = 8, &
     &   LIBXC_FAMILY_OEP = 16, &
     &   LIBXC_FAMILY_HYB_GGA = 32, &
     &   LIBXC_FAMILY_HYB_MGGA = 64, &
     &   LIBXC_FAMILY_HYB_LDA = 128
      INTEGER, SAVE :: LIBXC_MASK(2),XC_FAMILY(2),XC_KIND(2),LIBXC_NPARAM(2)
      DOUBLE PRECISION :: LIBXC_PARAM(2,200)
      CHARACTER(LEN=128) :: LIBXC_PARAM_DESCR(2,200)
      LOGICAL, SAVE :: LTBOUNDLIBXC, LIBXC_HYB(2), LIBXC_NOEXC(2)
      TYPE (xc_f03_func_t), SAVE :: LIBXC_FUNC_NSP(2),LIBXC_FUNC_SPN(2)
      TYPE (xc_f03_func_info_t), SAVE :: LIBXC_INFO_NSP(2),LIBXC_INFO_SPN(2)
#endif

!$ACC DECLARE CREATE(LEXCH, &
!$ACC&               EXCTAB, &
!$ACC&               LFCI, &
!$ACC&               LDAX, &
!$ACC&               ALDAC, &
!$ACC&               AGGAX, &
!$ACC&               AGGAC, &
!$ACC&               LDASCREEN, &
!$ACC&               LDASCREENC, &
!$ACC&               LRANGE_SEPARATED_CORR, &
!$ACC&               LUSE_LONGRANGE_HF, &
!$ACC&               LUSE_THOMAS_FERMI, &
!$ACC&               LUSE_MODEL_HF, &
!$ACC&               BPARAM, &
!$ACC&               CPARAM, &
!$ACC&               PARAM1, &
!$ACC&               PARAM2, &
!$ACC&               EXCTAB%RHOEXC, &
!$ACC&               EXCTAB%NEXCHF, &
!$ACC&               EXCTAB%EXCTAB)

#ifdef USELIBXC
!$ACC DECLARE CREATE(LIBXC_MASK, &
!$ACC&               XC_FAMILY, &
!$ACC&               LIBXC_HYB, &
!$ACC&               LIBXC_NOEXC, &
!$ACC&               XC_KIND, &
!$ACC&               LIBXC_FUNC_NSP, &
!$ACC&               LIBXC_FUNC_SPN, &
!$ACC&               LTBOUNDLIBXC)
#endif

      END MODULE setexm_struct_def
