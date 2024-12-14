!#define dotiming
#include "symbol.inc"

!***********************************************************************
!
!> this module implements the exact exchange operator (Hartree Fock)
!
! - until 01.10.03 written by Robin Hirschl
! - 08.10.03 cleaned up gK
! - 11.10.03 forces added by gK
! - 12.10.03 screened exchange added by gK
! - 21.10.03 one center onsite PAW corrections included gK
! - 30.10.03 core valence exchange added by gK
! - 32.02.04 screened exchange updated gK
! - 32.09.04 downsampling on HF grid added gK
! - 32.09.04 k-points downsampling added gK
! - 32.10.04 stress added gK
! - 18.04.05 modification of singularity correction to treat
!            k-points not included in the exchange kernel
!***********************************************************************


MODULE fock
  USE prec
  USE poscar
  USE nonl_high
  USE augfast
  USE wave_high
  USE fock_multipole
#ifdef shmem_bcast_buffer
  USE mpi_shmem
#endif
  USE fock_glb

  IMPLICIT NONE

CONTAINS
!**********************************************************************
!
!> read all variables related to exchange correlation treatment
!> from the INCAR file
!> this set both the variables in the setex module and
!> in the local fock module
!
!> @details @ref openmp :
!> In case openmp::omp_nthreads/=1 the default for NBLOCK_FOCK
!> (fock::nblock_fock) is changed to 2* openmp::omp_nthreads.
!
!**********************************************************************

  SUBROUTINE XC_FOCK_READER(IU5, IU0, IU6, SZPREC, ISIF, ISYM, IALGO, & 
      OMEGA, NTYP, NIONS, IMIX, AMIX, AMIX_MAG, BMIX, BMIX_MAG, LVTOT, &
      AFQMC_SET)
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE afqmc_struct, ONLY: afqmc_settings
      USE vaspxml
      USE base
      USE setexm
      USE full_kpoints
      USE vdw_df_helper, ONLY : SET_vdW_IO_UNITS
      USE openmp,        ONLY : omp_nthreads
      USE reader_tags
      USE tutor, ONLY: vtutor, isError, isAlert, isWarning, isAdvice, LHFOneCenter, ISYM2, &
          NMAXFOCKAE, DeprENCUTFOCK, PRECFOCK, vdWKlimes, HFParam, OEPdouble
#ifdef USELIBXC
      USE xc_f03_lib_m
#endif
      IMPLICIT NONE

      TYPE(afqmc_settings), INTENT(IN) :: AFQMC_SET
      INTEGER IU5, IU0, IU6
      CHARACTER(LEN=12) :: SZPREC
      INTEGER ISIF            !< calculate stress tensor or not (defaults to 0 for HF)
      INTEGER ISYM            !< symmetry
      INTEGER IALGO           !< applied algorithm
      REAL(q) :: OMEGA
      INTEGER :: NTYP, NIONS
      INTEGER :: IMIX
      REAL(q) :: AMIX, AMIX_MAG, BMIX, BMIX_MAG
      LOGICAL :: LVTOT
  ! local
      INTEGER IDUM, N, IERR
      REAL(q) RDUM
      COMPLEX(q) CDUM
      LOGICAL LOPEN,LDUM,LHFCALC_READ,USE_LIBXC
      LOGICAL, EXTERNAL :: CALCULATE_RESPONSE_FUNCTIONS, USE_OEP_IN_GW, LDMATRIX, FAST_FOCK, CALCULATE_RPA_FORCES
      CHARACTER (LEN=1)  :: CHARAC
      CHARACTER (LEN=40) :: SZNAM
#ifdef USELIBXC
      CHARACTER (LEN=256) :: LIBXC_NAME
      INTEGER :: LIBXC1,LIBXC2,NMAX,VMAJOR,VMINOR,VMICRO
#endif

      AEXX=0.0_q
      AGGAX=1.0_q
      AGGAC=1.0_q
      LDAX =1.0_q
      ALDAC=1.0_q
      HFSCREEN=0.0_q
      HFSCREENC=0.0_q
      HFRCUT=0.0_q
      PARAM1=0.1234000000_q
      PARAM2=1.00000000_q
      PARAM3=0._q
      LUSE_VDW=.FALSE.
      Zab_VDW=-0.8491_q
! rVV10
      BPARAM = 6.3_q 
      CPARAM = 0.0093_q
      
      LRSCOR=.FALSE.
      LHFCALC=CALCULATE_RESPONSE_FUNCTIONS().OR.LDMATRIX()
      LRHFCALC=.FALSE.
      L_MODEL_HF=.FALSE.
      L_THOMAS_FERMI=.FALSE.

      CALL OPEN_INCAR_IF_FOUND(IU5, LOPEN)
! HF calculations
      LHFCALC_READ=.FALSE.
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LHFCALC', LHFCALC, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)
      IF (LHFCALC .AND. IERR==0 .AND. N>=1) THEN
         LHFCALC_READ=.TRUE.
         ! if LHFCALC was READ (and only then) default AEXX to 0.25
         AEXX=0.25_q
         N=COUNT_ELEMENTS(INCAR_F,"METAGGA")
         IF (N>=1) CALL vtutor%error("XC_FOCK_READER: HF can not be combined with meta-GGAs.")
      ENDIF
! EXXOEP default is from chi.F
      EXXOEP=0
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'EXXOEP', EXXOEP, IERR, WRITEXMLINCAR)
      EXXOEP=MAX(0,MIN(EXXOEP,4))

! LSYMGRAD apply symmetry to restore symmetry of the gradient
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LSYMGRAD', LSYMGRAD, IERR, WRITEXMLINCAR)

! HF calculations one center terms
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LHFONE', LHFONE, IERR, WRITEXMLINCAR)

! force LHFCALC in some cases
      IF (EXXOEP>0 .OR. CALCULATE_RESPONSE_FUNCTIONS()) THEN
         LHFCALC=.TRUE.
      ENDIF
! long range correlation
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LRSCOR', LRSCOR, IERR, WRITEXMLINCAR)

! long range HF exclusively
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LRHFCALC', LRHFCALC, IERR, WRITEXMLINCAR)

! model HF
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LMODELHF', L_MODEL_HF, IERR, WRITEXMLINCAR)
      IF (L_MODEL_HF) THEN
         LRHFCALC=.FALSE.
      END IF

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LTHOMAS', L_THOMAS_FERMI, IERR, WRITEXMLINCAR)

         FOURORBIT=0
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'FOURORBIT', FOURORBIT, IERR, WRITEXMLINCAR)
      FOURORBIT=MAX(0,MIN(FOURORBIT,2))

      ENCUT4O=-1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ENCUT4O', ENCUT4O, IERR, WRITEXMLINCAR)

! If range separation, THOMAS fermi screening or four orbital terms are demanded 
! the calculation of the Fock exchange should be switched on as well
      IF (LRSCOR .OR. LRHFCALC .OR. L_MODEL_HF .OR. L_THOMAS_FERMI .OR. FOURORBIT>0) THEN
         LHFCALC=.TRUE.
      ENDIF

      IF (LHFCALC .AND. LHFONE) THEN
         CALL vtutor%write(isError, LHFOneCenter)
      ENDIF

!
! read in the GGA tag (enforces GGA, even if GGA was not used for PP creation)
!
      SZNAM='--'
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'GGA', SZNAM, 40, IERR, WRITEXMLINCAR)

      CALL STRIP(SZNAM,N,'L')
      CALL UPPER(SZNAM)
      SZGGA=SZNAM(1:2)

      CALL EXTYP(SZGGA,LEXCH)
! If we want a range separation in the LDA correlation then
! we have to switch to the Vosko-Wilk-Nusair treatment of
! LDA correlation
      IF (LRSCOR) CALL EXTYP('VW',LEXCH)
!
! read in LIBXC related tags
!
#ifdef USELIBXC
      NMAX=XC_F03_MAXIMUM_NAME_LENGTH()

      LIBXC1=-1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LIBXC1', LIBXC_NAME, NMAX, IERR, WRITEXMLINCAR, LCONTINUE=.TRUE., FOUNDNUMBER=N)
      IF (N>0) THEN
         CALL STRIP(LIBXC_NAME,N,'L')
         LIBXC1=XC_F03_FUNCTIONAL_GET_NUMBER(LIBXC_NAME(1:N))
      ENDIF
      IF (LIBXC1==-1) CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LIBXC1', LIBXC1, IERR, WRITEXMLINCAR, LCONTINUE=.TRUE.)

      LIBXC2=-1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LIBXC2', LIBXC_NAME, NMAX, IERR, WRITEXMLINCAR, LCONTINUE=.TRUE., FOUNDNUMBER=N)
      IF (N>0) THEN
         CALL STRIP(LIBXC_NAME,N,'L')
         LIBXC2=XC_F03_FUNCTIONAL_GET_NUMBER(LIBXC_NAME(1:N))
      ENDIF
      IF (LIBXC2==-1) CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LIBXC2', LIBXC2, IERR, WRITEXMLINCAR, LCONTINUE=.TRUE.)

      USE_LIBXC=.FALSE.
      SZNAM='--'
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'GGA', SZNAM, 40, IERR, WRITEXMLINCAR)
      CALL STRIP(SZNAM,N,'L')
      CALL UPPER(SZNAM)
      IF (SZNAM(1:5)=='LIBXC') USE_LIBXC=.TRUE.
      SZNAM='--'
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'METAGGA', SZNAM, 40, IERR, WRITEXMLINCAR)
      CALL STRIP(SZNAM,N,'L')
      CALL UPPER(SZNAM)
      IF (SZNAM(1:5)=='LIBXC') USE_LIBXC=.TRUE.

      IF ((LIBXC1==-1).AND.(LIBXC2==-1).AND.(USE_LIBXC.EQV..TRUE.)) THEN
         CALL vtutor%error("Error: Libxc is attempted to be used, but neither LIBXC1 nor LIBXC2 is specified.")
      ENDIF

      IF (((LIBXC1/=-1).OR.(LIBXC2/=-1)).AND.(USE_LIBXC.EQV..TRUE.)) THEN
         CALL XC_F03_VERSION(VMAJOR,VMINOR,VMICRO)
         WRITE(LIBXC_VERSION,'("     T    Libxc version ",I1,".",I1,".",I1)') VMAJOR, VMINOR, VMICRO
      ENDIF

      CALL EXTYP_LIBXC(LIBXC1,LIBXC2)
      IF (LHFCALC) THEN
         LHFCALC_READ=.TRUE.
         N=COUNT_ELEMENTS(INCAR_F,"METAGGA")
         IF (N>=1) CALL vtutor%error("XC_FOCK_READER: HF can not be combined with meta-GGAs.")
      ENDIF
#endif
!
! spin interpolation according to Vosko Wilk and Nusair
!
      LFCI=0
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'VOSKOWN', LFCI, IERR, WRITEXMLINCAR)
!
! parameters for x-functionals
!
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'PARAM1', PARAM1, IERR, WRITEXMLINCAR)

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'PARAM2', PARAM2, IERR, WRITEXMLINCAR)

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'PARAM3', PARAM3, IERR, WRITEXMLINCAR)

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LUSE_VDW', LUSE_VDW, IERR, WRITEXMLINCAR)

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'Zab_VDW', ZAB_VDW, IERR, WRITEXMLINCAR)
!
! read rVV10 parameters
!
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'BPARAM', BPARAM, IERR, WRITEXMLINCAR)

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'CPARAM', CPARAM, IERR, WRITEXMLINCAR)
!
! read in MODEL_GW flag
!
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MODEL_GW', MODEL_GW, IERR, WRITEXMLINCAR)
      IF (MODEL_GW>0) THEN
         LHFCALC=.TRUE.
         ! default the screening length to typical values for semiconductors
         HFSCREEN=2.0
         ! switch local exchange and correlation to Coulomb hole
         CALL SET_LEXCH(100)
      ENDIF
!
! read model screening constants
!
      MODEL_GW_EPS0=0.6_q*OMEGA/NIONS
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MODEL_EPS0', MODEL_GW_EPS0, IERR, WRITEXMLINCAR)

      MODEL_GW_ALPHA=1._q
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MODEL_ALPHA', MODEL_GW_ALPHA, IERR, WRITEXMLINCAR)
      !
      ! NKRED reduce k-points in Fock exchange
      !
      NKREDX=1
      NKREDY=1
      NKREDZ=1
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NKRED', NKREDX, IERR, WRITEXMLINCAR, FOUNDNUMBER=N)

      IF (N>=1) THEN
         NKREDY=NKREDX
         NKREDZ=NKREDX
      ENDIF

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NKREDX', NKREDX, IERR, WRITEXMLINCAR)

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NKREDY', NKREDY, IERR, WRITEXMLINCAR)

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NKREDZ', NKREDZ, IERR, WRITEXMLINCAR)

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'SHIFTRED', SHIFTRED, IERR, WRITEXMLINCAR)

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'HFKIDENT', HFKIDENT, IERR, WRITEXMLINCAR)

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ODDONLY', ODDONLY, IERR, WRITEXMLINCAR)

      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'EVENONLY', EVENONLY, IERR, WRITEXMLINCAR)

      IF (LHFCALC) THEN
!
! read lowest band to be included in HF term
!
         NBANDSGWLOW_FOCK=1
     
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBANDSGWLOW', NBANDSGWLOW_FOCK, IERR, WRITEXMLINCAR)
         NBANDSGWLOW_FOCK=MAX(NBANDSGWLOW_FOCK,1)
!
! reset ISIF to 0, and reread it from file
!
         ISIF=0
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ISIF', ISIF, IERR, WRITEXMLINCAR)
         IF (EXXOEP==0 .AND. .NOT. USE_OEP_IN_GW()) ISYM=3

         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ISYM', ISYM, IERR, WRITEXMLINCAR)
         IF (EXXOEP==0 .AND. .NOT. USE_OEP_IN_GW() .AND. ISYM == 2) THEN
            CALL vtutor%write(isAlert, ISYM2)
         ENDIF

! re-read BMIX tag (default is BMIX=0.001 for OEP methods)

         IF (EXXOEP/=0 .OR. USE_OEP_IN_GW()) THEN
            AMIX=0.4
            CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'AMIX', AMIX, IERR, WRITEXMLINCAR)

            AMIX_MAG=AMIX
            CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'AMIX_MAG', AMIX_MAG, IERR, WRITEXMLINCAR)

            BMIX=0.001_q
            CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'BMIX', BMIX, IERR, WRITEXMLINCAR)

            BMIX_MAG=BMIX
            CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'BMIX_MAG', BMIX_MAG, IERR, WRITEXMLINCAR)

            LVTOT=.TRUE.
            CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LVTOT', LVTOT, IERR, WRITEXMLINCAR)
         ELSE
! no Broyden mixer for HF, tends to blow up, simple Kerker mixer as default
! hence reread IMIX with different default
            IF ((IALGO>=60 .OR. IALGO<50) .AND. LHFCALC_READ)  THEN
               IMIX=1
               CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'IMIX', IMIX, IERR, WRITEXMLINCAR)
            ENDIF
         ENDIF
!
! HFLMAX maximum L quantum number for HF (LMAXFOCK is read alternatively)
!
         LMAX_FOCK=4
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'HFLMAX', LMAX_FOCK, IERR, WRITEXMLINCAR)

         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LMAXFOCK', LMAX_FOCK, IERR, WRITEXMLINCAR)

         ! default for LMAX_FOCKAE is set in chi.F (RESPONSE_READER)
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LMAXFOCKAE', LMAX_FOCKAE, IERR, WRITEXMLINCAR)

         ! values below -1 are not sensible (-1 no AE augmentation)
         LMAX_FOCKAE=MAX(-1,LMAX_FOCKAE)

         ! default for QMAX_FOCKAE is set in chi.F (RESPONSE_READER)
         ALLOCATE(QMAX_FOCKAE(NTYP))
         QMAX_FOCKAE=0
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'QMAXFOCKAE', QMAX_FOCKAE, NTYP, IERR, WRITEXMLINCAR)

         NMAX_FOCKAE=1
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NMAXFOCKAE', NMAX_FOCKAE, IERR, WRITEXMLINCAR)
         NMAX_FOCKAE=MIN(2,MAX(1,NMAX_FOCKAE))

         IF (NMAX_FOCKAE>=2) THEN
            CALL vtutor%write(isAlert, NMAXFOCKAE)
         ENDIF

         IF ((EXXOEP/=0 .OR. USE_OEP_IN_GW()).AND. LMAX_FOCKAE>=0) THEN
            LFOCKAEDFT=.TRUE.
         ELSE
            LFOCKAEDFT=.FALSE.
         ENDIF

         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LFOCKAEDFT', LFOCKAEDFT, IERR, WRITEXMLINCAR)
!
! read in HFALPHA from INCAR file, default value -1 (automatic determination)
!
         HFALPHA=-1

         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'HFALPHA', HFALPHA, IERR, WRITEXMLINCAR)
!
! read in MCALPHA from INCAR file, default value -1 (automatic determination)
!
         MCALPHA=0
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'MCALPHA', MCALPHA, IERR, WRITEXMLINCAR)
!
!
! read in HFSCREEN, screening length in Thomas-Fermi and Scuserias 
! screened exchange, as well as in the model GW
! Thomas-Fermi and model GW: q_TF
! Scuserias default: 0.6 A-1
!
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'HFSCREEN', HFSCREEN, IERR, WRITEXMLINCAR)
!
! screening length for correlation
!
         HFSCREENC=HFSCREEN
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'HFSCREENC', HFSCREENC, IERR, WRITEXMLINCAR)
!
! read in HFRCUT
! truncate the potential kernel at a finite radius R
!
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'HFRCUT', HFRCUT, IERR, WRITEXMLINCAR)
!
! read in AEXX amount of exact exchange and AGGAX (default 1-AEXX)
!
         ! PBE0 parameters
         ! C. Adamo, V. Barone, J. Chem. Phys. 110, 6158 (1999)
         ! M. Ernzerhof, G. E. Scuseria, J. Chem. Phys. 110, 5029 (1999)
         ! B3LYP parameters  AEXX=0.2 ; AGGAX=0.72 ; AGGAC=0.81

         IF (LRHFCALC .OR. L_THOMAS_FERMI .OR. MODEL_GW>0 .OR. EXXOEP/=0 .OR. USE_OEP_IN_GW()) THEN
            AEXX=1.0
         ENDIF
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'AEXX', AEXX, IERR, WRITEXMLINCAR)

         AGGAX=1-AEXX
         LDAX =1-AEXX

         ! for HF type calculations and full exchange, no correlation is usually included
         IF (AEXX==1.0) THEN
            ALDAC=0.0_q
            AGGAC=0.0_q
         ENDIF

         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'AGGAX', AGGAX, IERR, WRITEXMLINCAR)

         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ALDAX', LDAX, IERR, WRITEXMLINCAR)
!
! read in ENCUTFOCK cutoff for Fock exchange
! terminate if user selects it
!
         ENCUTFOCK=-1
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ENCUTFOCK', ENCUTFOCK, IERR, WRITEXMLINCAR)
!         IF (N>0) THEN
!           CALL vtutor%write(isError, DeprENCUTFOCK)
!         ENDIF

! HF calculations one center terms
         SZPRECFOCK='      '
         IF (FAST_FOCK()) SZPRECFOCK='fast'
         IF (SZPREC(1:1)=='s') SZPRECFOCK=SZPREC
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'PRECFOCK', SZPRECFOCK, 40, IERR, WRITEXMLINCAR)
         CALL STRIP(SZPRECFOCK,N,'L')
         CALL LOWER(SZPRECFOCK)
         IF (SZPRECFOCK(1:1)=='l') THEN
            SZPRECFOCK='low   '
         ELSE IF (SZPRECFOCK(1:1)=='m') THEN
            SZPRECFOCK='medium'
         ELSE IF (SZPRECFOCK(1:1)=='f') THEN
            SZPRECFOCK='fast'
         ELSE IF (SZPRECFOCK(1:1)=='n') THEN
            SZPRECFOCK='normal'
         ELSE IF (SZPRECFOCK(1:1)=='a') THEN
            SZPRECFOCK='accur.'
         ELSE IF (SZPRECFOCK(1:4)=='soft') THEN
            SZPRECFOCK='soft  '
         ELSE IF (SZPRECFOCK(1:1)=='s') THEN
            SZPRECFOCK='single'
         ELSE
            SZPRECFOCK='normal'
!            CALL vtutor%write(isAlert, PRECFOCK)
         ENDIF
!
! Use Adaptively Compressed Exchange (ACE) representation of
! the fock operator?
!
         LFOCKACE=.TRUE.
#ifndef fock_dblbuf
         LFOCKACE=.FALSE.
#endif
         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LFOCKACE', LFOCKACE, IERR, WRITEXMLINCAR)
#ifndef fock_dblbuf
         IF( LFOCKACE ) THEN
            CALL vtutor%error("The adaptively compressed exchange method LFOCKACE=.TRUE. requires &
               &to compile the code using -Dfock_dblbuf.")
         ENDIF
#endif
      ELSE IF (LHFONE) THEN
!
! read in AEXX amount of exact exchange and AGGAX (default 1-AEXX)
! one center PBE0 parameter
!
         AEXX=0.25

         CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'AEXX', AEXX, IERR, WRITEXMLINCAR)
      ENDIF
!
! amount of GGA correlation
!
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'AGGAC', AGGAC, IERR, WRITEXMLINCAR)
!
! amount of LDA correlation
!
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ALDAC', ALDAC, IERR, WRITEXMLINCAR)
!
! blocking factor in FOCK_ACC and FOCK_FORCE
!
      NBLOCK_FOCK=64
!$    IF (omp_nthreads>1) NBLOCK_FOCK=2*omp_nthreads
      CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBLOCK_FOCK', NBLOCK_FOCK, IERR, WRITEXMLINCAR)

      CALL CLOSE_INCAR_IF_FOUND(IU5)

      IF (LUSE_VDW) THEN
         CALL vtutor%write(isAdvice, vdWKlimes)
         CALL SET_vdW_IO_UNITS(IU0, IU6, ISIF)
      ENDIF
!
! set the corresponding variables in the LDA module (xclib.F)
! 
      LDASCREEN =HFSCREEN
      LDASCREENC=HFSCREENC
      LRANGE_SEPARATED_CORR=LRSCOR
      LUSE_LONGRANGE_HF=LRHFCALC
      LUSE_THOMAS_FERMI=L_THOMAS_FERMI
      LUSE_MODEL_HF=L_MODEL_HF

#ifdef USELIBXC
      CALL SETPARAMS_LIBXC
#endif

      IF (NKREDX==1 .AND. NKREDY==1 .AND. NKREDZ==1 .AND. .NOT. EVENONLY .AND. .NOT. ODDONLY ) & 
        HFKIDENT=.FALSE.

      IF (LHFCALC .AND. HFRCUT/=0 .AND. HFKIDENT) THEN
          CALL vtutor%write(isError, HFParam)
      ENDIF

      IF (LHFCALC) THEN
         CALL USE_FULL_KPOINTS
      ENDIF

      ! OEP methods do not work for low, fast or soft (double grid techniques)
      IF (EXXOEP>0) THEN
         IF ( SZPRECFOCK(1:1) /= SZPREC(1:1) .AND. SZPRECFOCK(1:1) /= 'm') THEN
            CALL vtutor%write(isWarning, OEPdouble)
            SZPRECFOCK='normal'
         ENDIF
      ENDIF

!$ACC UPDATE DEVICE(LEXCH, &
!$ACC&              LFCI, &
!$ACC&              LDAX, &
!$ACC&              ALDAC, &
!$ACC&              AGGAX, &
!$ACC&              AGGAC, &
!$ACC&              LDASCREEN, &
!$ACC&              LDASCREENC, &
!$ACC&              LRANGE_SEPARATED_CORR, &
!$ACC&              LUSE_LONGRANGE_HF, &
!$ACC&              LUSE_THOMAS_FERMI, &
!$ACC&              LUSE_MODEL_HF, &
!$ACC&              BPARAM, &
!$ACC&              CPARAM, &
!$ACC&              PARAM1, &
!$ACC&              PARAM2) ASYNC(ACC_ASYNC_Q)
    END SUBROUTINE XC_FOCK_READER


!**********************************************************************
!
!> write the Fock and exchange correlation related
!> parameters to the OUTCAR file
!
!**********************************************************************

  SUBROUTINE WRITE_FOCK(IU6)
      USE setexm
      IMPLICIT NONE

      INTEGER I,J,IU6
      REAL(q), EXTERNAL :: ENCUTGW_IN_CHI

      REAL(q) :: ENCUTGW
      CHARACTER (20) ICHAR

      ENCUTGW=ENCUTGW_IN_CHI()

      IF (LHFCALC) THEN
         IF (LMAX_FOCKAE==0) LFOCKAEDFT=.FALSE.
      ELSE
         LFOCKAEDFT=.FALSE.
      ENDIF

      IF (IU6>=0 .AND. LHFCALC) THEN
         WRITE(IU6,10 ) SZGGA,LEXCH,LIBXC_VERSION,LFCI,EXXOEP,LHFCALC,LSYMGRAD,SZPRECFOCK,LRHFCALC,LRSCOR,L_THOMAS_FERMI,L_MODEL_HF,LFOCKACE,ENCUT4O, &
           LMAX_FOCK,LMAX_FOCKAE,NMAX_FOCKAE,LFOCKAEDFT,NKREDX,NKREDY,NKREDZ,SHIFTRED,HFKIDENT,ODDONLY,EVENONLY,HFALPHA,MCALPHA,AEXX,HFSCREEN,HFSCREENC,HFRCUT, &
           LDAX,AGGAX,ALDAC,AGGAC,ENCUTFOCK,NBANDSGWLOW_FOCK,NBLOCK_FOCK
#ifdef USELIBXC
         IF ((LIBXC_MASK(1)>0).OR.(LIBXC_MASK(2)>0)) THEN
            WRITE(IU6,'(A)') ' Parameters of Libxc functionals:'
            DO I=1,2
               IF (LIBXC_MASK(I)>0) THEN
                  IF (I==1) THEN
                     WRITE(IU6,'(3A,I3,A)') '   LIBXC1(',TRIM(XC_F03_FUNCTIONAL_GET_NAME(LIBXC_MASK(I))),'=',LIBXC_MASK(I),'):'
                  ELSEIF (I==2) THEN
                     WRITE(IU6,'(3A,I3,A)') '   LIBXC2(',TRIM(XC_F03_FUNCTIONAL_GET_NAME(LIBXC_MASK(I))),'=',LIBXC_MASK(I),'):'
                  ENDIF
                  IF (LIBXC_NPARAM(I)==0) THEN
                     WRITE(IU6,'(A)') '     No parameter'
                  ELSE
                     DO J=1, LIBXC_NPARAM(I)
                        WRITE(ICHAR,*) J
                        WRITE(IU6,'(5A,E14.6)') '     P',TRIM(ADJUSTL(ICHAR)),': ',TRIM(LIBXC_PARAM_DESCR(I,J)),' =',LIBXC_PARAM(I,J)
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
            WRITE(IU6,*)
         ENDIF
#endif
      ENDIF
      IF (IU6>=0 .AND. .NOT.LHFCALC) THEN
         WRITE(IU6,100) SZGGA,LEXCH,LIBXC_VERSION,LFCI,LHFCALC,LHFONE,AEXX
#ifdef USELIBXC
         IF ((LIBXC_MASK(1)>0).OR.(LIBXC_MASK(2)>0)) THEN
            WRITE(IU6,'(A)') ' Parameters of Libxc functionals:'
            DO I=1,2
               IF (LIBXC_MASK(I)>0) THEN
                  IF (I==1) THEN
                     WRITE(IU6,'(3A,I3,A)') '   LIBXC1(',TRIM(XC_F03_FUNCTIONAL_GET_NAME(LIBXC_MASK(I))),'=',LIBXC_MASK(I),'):'
                  ELSEIF (I==2) THEN
                     WRITE(IU6,'(3A,I3,A)') '   LIBXC2(',TRIM(XC_F03_FUNCTIONAL_GET_NAME(LIBXC_MASK(I))),'=',LIBXC_MASK(I),'):'
                  ENDIF
                  IF (LIBXC_NPARAM(I)==0) THEN
                     WRITE(IU6,'(A)') '     No parameter'
                  ELSE
                     DO J=1, LIBXC_NPARAM(I)
                        WRITE(ICHAR,*) J
                        WRITE(IU6,'(5A,E14.6)') '     P',TRIM(ADJUSTL(ICHAR)),': ',TRIM(LIBXC_PARAM_DESCR(I,J)),' =',LIBXC_PARAM(I,J)
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
            WRITE(IU6,*)
         ENDIF
#endif
      ENDIF
      IF (IU6>=0 .AND. MODEL_GW>0) THEN
         WRITE(IU6,20) MODEL_GW,MODEL_GW_EPS0, MODEL_GW_ALPHA
20    FORMAT(' Model GW treatment:'  / &
           '   MODEL_GW   =',I6,  '    type of model GW         '/ & 
           '   MODEL_EPSO =',F6.2,  '    static dielectric constant'/&
           '   MODEL_ALPHA=',F6.2,  '    alpha'/)
      ENDIF

      IF (IU6>=0 .AND. LUSE_VDW) THEN
         WRITE(IU6,40) LUSE_VDW,Zab_VDW,PARAM1,PARAM2,PARAM3,BPARAM,CPARAM
40    FORMAT(' vdW DFT:'  / &
           '   LUSE_VDW   =',L6,  '    switch on vdW DFT'/ & 
           '   Zab_VDW    =',F6.4,  '    correlation parameter'/&
           '   PARAM1     =',F6.4,  /&
           '   PARAM2     =',F6.4,  /&
           '   PARAM3     =',F6.4,  /&
           '   -- rVV10 --',           /&
           '   BPARAM     =',F7.4,  /&
           '   CPARAM     =',F7.4,  /)
      ENDIF


      IF (ENCUTGW>=0 .AND. IU6>=0 .AND. EXXOEP /= 0) THEN
         WRITE(IU6,30) ENCUTGW
30    FORMAT(' Cutoff for OEP:'  / &
           '   ENCUTGW    =',F6.1,  '    cutoff applied to OEP potential'/)
      ENDIF
      
  10  FORMAT(' Exchange correlation treatment:'  / &
             '   GGA     =',A6,  '    GGA type'/ & 
             '   LEXCH   =',I6,  '    internal setting for exchange type'/ & 
             '   LIBXC   =',A35, / &
             '   VOSKOWN =',I6,  '    Vosko Wilk Nusair interpolation'/&
             '   EXXOEP  =',I6,  '    0=HF, 1=EXX-LHF (local Hartree Fock) 2=EXX OEP'/ &
             '   LHFCALC =',L6,  '    Hartree Fock is set to'/ &
             '   LSYMGRAD=',L6,  '    symmetrize gradient (conserves proper symmetry)'/ &
             '   PRECFOCK=',A6,  '    Normal, Fast or Accurate (Low or Medium for compatibility)'/ &
             '   LRHFCALC=',L6,  '    long range Hartree Fock'/ &
             '   LRSCOR  =',L6,  '    long range correlation only (use DFT for short range part)'/ &
             '   LTHOMAS =',L6,  '    Thomas Fermi screening in HF'/ &
             '   LMODELHF=',L6,  '    short range full HF, long range fraction AEXX'/ &
             '   LFOCKACE=',L6  ,'    use Adeptively-Compressed-Exchange operator'/ &
             '   ENCUT4O =',F6.1,'   cutoff for four orbital integrals eV'/ &
             '   LMAXFOCK=',I6  ,'    L truncation for augmentation on plane wave grid'  / &
             '   LMAXFOCKAE=',I4,'    L truncation for all-electron charge restoration on plane wave grid'  / &
             '   NMAXFOCKAE=',I4,'    number of basis functions for all-electron charge restoration' / &
             '   LFOCKAEDFT=',L6,  '  apply the AE augmentation even for DFT' /&
             '   NKREDX  =',I6  ,'    reduce k-point grid by'    / &
             '   NKREDY  =',I6  ,'    reduce k-point grid by'    / &
             '   NKREDZ  =',I6  ,'    reduce k-point grid by'    / &
             '   SHIFTRED=',L6  ,'    shift reduced grid of Gamma'    / &
             '   HFKIDENT=',L6  ,'    idential grid for each k-point'    / &
             '   ODDONLY =',L6  ,'    use only odd q-grid points'    / &
             '   EVENONLY=',L6  ,'    use only even q-grid points'    / &
             '   HFALPHA =',F10.4,' decay constant for conv. correction'/ &
             '   MCALPHA =',F10.4,' extent of test-charge in conv. correction in multipole expansion'/ &
             '   AEXX    =',F10.4,' exact exchange contribution' / &
             '   HFSCREEN=',F10.4,' screening length (either q_TF or 0.3 A-1)'/ &
             '   HFSCREENC=',F9.4,' screening length for correlation (either q_TF or 0.3 A-1)'/ &
             '   HFRCUT  =',F10.4,' spherical cutoff for potential kernel'/ &
             '   ALDAX   =',F10.4,' LDA exchange part'           / &
             '   AGGAX   =',F10.4,' GGA exchange part'           / &
             '   ALDAC   =',F10.4,' LDA correlation'             / &
             '   AGGAC   =',F10.4,' GGA correlation'             / &
             '   ENCUTFOCK=',F6.1,' apply spherical cutoff to Coloumb kernel'/ &
             '   NBANDSGWLOW=',I6  ,'    first orbital included in HF term'/ &
             '   NBLOCK_FOCK=',I6  ,'    blocking factor in FOCK_ACC'/ &
             )


  100 FORMAT(' Exchange correlation treatment:'  / &
             '   GGA     =',A6,  '    GGA type'/ & 
             '   LEXCH   =',I6,  '    internal setting for exchange type'/ & 
             '   LIBXC   =',A35, / &
             '   VOSKOWN =',I6,  '    Vosko Wilk Nusair interpolation'/&
             '   LHFCALC =',L6,  '    Hartree Fock is set to'/&
             '   LHFONE  =',L6,  '    Hartree Fock one center treatment'/ &
             '   AEXX    =',F10.4,' exact exchange contribution' / &
             )
       
    END SUBROUTINE WRITE_FOCK


!**********************************************************************
!
!>      XML_WRITE_XC_FOCK
!>   this subroutine writes the exchange parameters and the Fock
!>   parameters to a file
!
!**********************************************************************

    SUBROUTINE XML_WRITE_XC_FOCK
      USE setexm
      USE vaspxml
      IMPLICIT NONE
! local
      INTEGER IDUM, N
      REAL(q) RDUM
      COMPLEX(q) CDUM
      LOGICAL LDUM
      CHARACTER (1) :: CHARAC
 
      CALL XML_TAG("separator","electronic exchange-correlation")

      N=1
      CALL XML_INCAR('GGA','S',IDUM,RDUM,CDUM,LDUM,SZGGA,N)
      CALL XML_INCAR('VOSKOWN','I',LFCI,RDUM,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('LHFCALC','L',IDUM,RDUM,CDUM,LHFCALC,CHARAC,N)
      CALL XML_INCAR('PRECFOCK','S',IDUM,RDUM,CDUM,LDUM,SZPRECFOCK,N)
      CALL XML_INCAR('LSYMGRAD','L',IDUM,RDUM,CDUM,LSYMGRAD,CHARAC,N)
      CALL XML_INCAR('LHFONE','L',IDUM,RDUM,CDUM,LHFONE,CHARAC,N)
      CALL XML_INCAR('LRHFCALC','L',IDUM,RDUM,CDUM,LRHFCALC,CHARAC,N)
      CALL XML_INCAR('LTHOMAS','L',IDUM,RDUM,CDUM,L_THOMAS_FERMI,CHARAC,N)
      CALL XML_INCAR('LMODELHF','L',IDUM,RDUM,CDUM,L_MODEL_HF,CHARAC,N)
      CALL XML_INCAR('LFOCKACE','L',IDUM,RDUM,CDUM,LFOCKACE,CHARAC,N)
      CALL XML_INCAR('ENCUT4O','F',IDUM,ENCUT4O,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('EXXOEP','I',EXXOEP,RDUM,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('FOURORBIT','I',FOURORBIT,RDUM,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('AEXX','F',IDUM,AEXX,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('HFALPHA','F',IDUM,HFALPHA,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('MCALPHA','F',IDUM,MCALPHA,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('ALDAX','F',IDUM,LDAX,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('AGGAX','F',IDUM,AGGAX,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('ALDAC','F',IDUM,ALDAC,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('AGGAC','F',IDUM,AGGAC,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('NKREDX','I',NKREDX,RDUM,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('NKREDY','I',NKREDY,RDUM,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('NKREDZ','I',NKREDZ,RDUM,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('SHIFTRED','L',IDUM,RDUM,CDUM,SHIFTRED,CHARAC,N)
      CALL XML_INCAR('ODDONLY','L',IDUM,RDUM,CDUM,ODDONLY,CHARAC,N)
      CALL XML_INCAR('EVENONLY','L',IDUM,RDUM,CDUM,EVENONLY,CHARAC,N)
      CALL XML_INCAR('LMAXFOCK','I',LMAX_FOCK,RDUM,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('NMAXFOCKAE','I',NMAX_FOCKAE,RDUM,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('LFOCKAEDFT','L',IDUM,RDUM,CDUM,LFOCKAEDFT,CHARAC,N)
      CALL XML_INCAR('HFSCREEN','F',IDUM,HFSCREEN,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('HFSCREENC','F',IDUM,HFSCREENC,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('NBANDSGWLOW','I',NBANDSGWLOW_FOCK,RDUM,CDUM,LDUM,CHARAC,N)

      CALL XML_CLOSE_TAG
      CALL XML_TAG("separator","vdW DFT")

      CALL XML_INCAR('LUSE_VDW','L',IDUM,RDUM,CDUM,LUSE_VDW,CHARAC,N)
      CALL XML_INCAR('Zab_VDW','F',IDUM,Zab_VDW,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('PARAM1','F',IDUM,PARAM1,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('PARAM2','F',IDUM,PARAM2,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('PARAM3','F',IDUM,PARAM3,CDUM,LDUM,CHARAC,N)

      CALL XML_CLOSE_TAG

      CALL XML_TAG("separator","model GW")
      CALL XML_INCAR('MODEL_GW','I',MODEL_GW,RDUM,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('MODEL_EPS0','F',IDUM,MODEL_GW_EPS0,CDUM,LDUM,CHARAC,N)
      CALL XML_INCAR('MODEL_ALPHA','F',IDUM,MODEL_GW_ALPHA,CDUM,LDUM,CHARAC,N)
      CALL XML_CLOSE_TAG


    END SUBROUTINE XML_WRITE_XC_FOCK


!***********************************************************************
!
!> initialise the augmentation descriptors for the
!> fast augmentation
!> here is what is done:
!> ) setup grid for HF operator and some stuff for nonlocal terms ((c) Georg)
!>   ) this is alligned in real space with GRID (and GRID_US)
!>   ) in reciprocal space it uses only half of the G-vectors
!>     if gammareal is selected
!>   ) points are distributed among the NPAR processors in COMM_INB
!
!***********************************************************************

    SUBROUTINE SETUP_FOCK(T_INFO, P, WDES, GRID, LATT_CUR, LMDIM, SZPREC, IU6, IU0 )
#ifdef _OPENACC
      USE mcufft
#endif
      USE pseudo
      USE constant
      USE asa
      USE wave
      USE mgrid
      USE lattice
      USE full_kpoints
      USE wpbe
      USE tutor, ONLY: vtutor, isError, HF_NPAR, HFUSPP
      
      IMPLICIT NONE

      TYPE (type_info) T_INFO
      TYPE (potcar),TARGET::  P(T_INFO%NTYP)
      TYPE (wavedes), TARGET :: WDES !< wave function descriptor
      TYPE (grid_3d), TARGET :: GRID !< 3d- grid descriptor for wavefunctions
      TYPE (latt)    :: LATT_CUR     !< lattice parameters
      INTEGER        :: LMDIM        !< leading dimension of arrays like CDIJ
      INTEGER :: IU6, IU0            !< io units
      CHARACTER(LEN=12) :: SZPREC
    ! local
      INTEGER :: NT, L, NK, IND
      INTEGER :: N1, N2, N3
      REAL(q) :: G1, G2, G3, GIX, GIY, GIZ, ENERGI
      INTEGER :: NGX_FOCK, NGY_FOCK, NGZ_FOCK
      REAL(q) :: WFACT
    ! work array to perform 
      REAL(q) :: ENMAX
      REAL(q) :: SCALE
      LOGICAL, EXTERNAL :: CALCULATE_RESPONSE_FUNCTIONS,LDMATRIX

      IF (.NOT. LHFCALC .AND. .NOT. CALCULATE_RESPONSE_FUNCTIONS() .AND. .NOT. LDMATRIX()) RETURN
#ifdef MPI
!      IF (WDES%GRID%COMM%NCPU /=1 .OR. LUSE_PARALLEL_FFT) THEN
!         CALL vtutor%write(isError, HF_NPAR)
!      ENDIF
#endif
      IF (IU6>=0) CALL WFORCE(IU6)

      ! check if PPs are OK
      typ:  DO NT=1,T_INFO%NTYP        
         IF ((WDES%LOVERL) .AND. ( .NOT. ASSOCIATED( P(NT)%QPAW )) .AND. (P(NT)%PSDMAX/=0)) THEN
            CALL vtutor%write(isError, HFUSPP)
         ENDIF
      ENDDO typ
!==========================================================================
! generate the required GRID structures and the modified descriptors
!  GRID_FOCK and WDES_FOCK
!==========================================================================
      ! simple case, link the descriptors
      GRID_FOCK=> GRID
      WDES_FOCK=> WDES
      LWDES_FOCK=.FALSE.
      ! difficult case regenerate the grids
      IF ( SZPRECFOCK(1:1) /= SZPREC(1:1) .AND. SZPRECFOCK(1:1) /= 'm' .AND. SZPRECFOCK(1:6) /= 'single') THEN
          ! first allocate the required structures and initialise to default
          ALLOCATE(GRID_FOCK)
          ALLOCATE(WDES_FOCK)
          GRID_FOCK=GRID
          WDES_FOCK=WDES
          LWDES_FOCK=.TRUE.

          IF (SZPRECFOCK(1:1)=='a') THEN
             ! high precision do not allow for wrap around
             WFACT=4
          ELSE 
             ! medium-low precision allow for small wrap around
             WFACT=3
          ENDIF
          
          ! determine minimum required values for NGX, Y and NGZ
          IF (SZPRECFOCK(1:1)=='a' .OR. SZPRECFOCK(1:1)=='n') THEN
             NGX_FOCK =SQRT(WDES%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(1)/AUTOA))*WFACT+1.0_q
             NGY_FOCK =SQRT(WDES%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(2)/AUTOA))*WFACT+1.0_q
             NGZ_FOCK= SQRT(WDES%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(3)/AUTOA))*WFACT+1.0_q
          ELSE
             NGX_FOCK=0
             NGY_FOCK=0
             NGZ_FOCK=0
          ENDIF

          ! new version loops over all k-points in the entire BZ
          ! to determine the required NGX
          DO N3=1,GRID%NGZ_rd
             DO N2=1,GRID%NGY
                DO NK=1,WDES%NKPTS_FOR_GEN_LAYOUT
                   G3=(GRID%LPCTZ(N3)+WDES%VKPT(3,NK))
                   G2=(GRID%LPCTY(N2)+WDES%VKPT(2,NK))
                   DO N1=1,GRID%NGX_rd
                      G1=(GRID%LPCTX(N1)+WDES%VKPT(1,NK))
                      GIX= (G1*LATT_CUR%B(1,1)+G2*LATT_CUR%B(1,2)+G3*LATT_CUR%B(1,3)) *TPI
                      GIY= (G1*LATT_CUR%B(2,1)+G2*LATT_CUR%B(2,2)+G3*LATT_CUR%B(2,3)) *TPI
                      GIZ= (G1*LATT_CUR%B(3,1)+G2*LATT_CUR%B(3,2)+G3*LATT_CUR%B(3,3)) *TPI
                      ENERGI=HSQDTM*((GIX**2)+(GIY**2)+(GIZ**2))

                      ! exclude some components for gamma-only version (C(G)=C*(-G))
                      IF (GRID%NGZ/=GRID%NGZ_rd) THEN
                         IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)<0) CYCLE
                         IF (GRID%LPCTZ(N3)==0 .AND.GRID%LPCTY(N2)==0 .AND.GRID%LPCTX(N1)<0) CYCLE
                      ENDIF

                      IF(ENERGI<WDES%ENMAX) THEN
                         NGX_FOCK=MAX(NGX_FOCK,ABS(GRID%LPCTX(N1))*2+2)
                         NGY_FOCK=MAX(NGY_FOCK,ABS(GRID%LPCTY(N2))*2+2)
                         NGZ_FOCK=MAX(NGZ_FOCK,ABS(GRID%LPCTZ(N3))*2+2)
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

          ! in principle that should do
          ! however the previous algorithm sometimes leads to meshes that violate
          ! symmetry
          ! search for the next larger values that conserves symmetry
          WFACT=MAX(NGX_FOCK/LATT_CUR%ANORM(1),NGY_FOCK/LATT_CUR%ANORM(2),NGZ_FOCK/LATT_CUR%ANORM(3))+1E-6
          NGX_FOCK=WFACT*LATT_CUR%ANORM(1)
          NGY_FOCK=WFACT*LATT_CUR%ANORM(2)
          NGZ_FOCK=WFACT*LATT_CUR%ANORM(3)

          GRID_FOCK%NGPTAR(1)=NGX_FOCK
          GRID_FOCK%NGPTAR(2)=NGY_FOCK
          GRID_FOCK%NGPTAR(3)=NGZ_FOCK

          CALL FFTCHK(GRID_FOCK%NGPTAR)
          ! reinitialize the loop counters
          CALL INILGRD(GRID_FOCK%NGPTAR(1),GRID_FOCK%NGPTAR(2),GRID_FOCK%NGPTAR(3),GRID_FOCK)

          IF (IU6>=0) THEN
             WRITE(IU6,20) GRID_FOCK%NGPTAR
20           FORMAT(/' FFT grid for exact exchange (Hartree Fock) ',/ &
                  '  NGX =',I3,'; NGY =',I3,'; NGZ =',I3,/)
          ENDIF

          ! set WDES_FOCK%GRID to GRID_FOCK and generate the data layout for FFT's
          CALL GEN_LAYOUT(GRID_FOCK,WDES_FOCK, LATT_CUR%B, LATT_CUR%B, IU6, .TRUE.)
          CALL GEN_INDEX (GRID_FOCK,WDES_FOCK, LATT_CUR%B, LATT_CUR%B,-1, -1, .TRUE.)

          CALLMPI( MAPSET(GRID_FOCK))  ! generate the communication maps (patterns) for FFT
          ! now check whether it is equivalent to original layout in reciprocal space (G)
          ! this is required since the FFTEXT and FFTWAV would not work from the orbital arrays
          CALL CHECK_GEN_LAYOUT( WDES,  WDES_FOCK, WDES%NKPTS)
          !  init FFT (required if real to complex FFT is used)
          CALL FFTINI(WDES_FOCK%NINDPW(1,1),WDES_FOCK%NGVECTOR(1),WDES_FOCK%NKPTS,WDES_FOCK%NGDIM,GRID_FOCK)
          ! set plane wave index array INDPW between NKPTS and NKDIM
          CALL SET_INDPW_FOCK(GRID_FOCK,WDES_FOCK)
      ENDIF
!==========================================================================
! other initializations
!==========================================================================
      ! set the indexing arrays for k-points outside IRZ
      CALL SET_INDPW_FOCK(GRID,WDES)

      GRIDHF%COMM=>WDES_FOCK%COMM_INB
      CALL INILGRD(GRID_FOCK%NGX,GRID_FOCK%NGY,GRID_FOCK%NGZ,GRIDHF)
      CALL GEN_GRID_HF(GRIDHF, GRID_FOCK)

      !  generate the communication maps (patterns) for FFT
      CALLMPI( MAPSET(GRIDHF))

      ! create (and destroy) the FFT plans to generate "wisdom"
      CALL FFTGRIDPLAN(GRIDHF   )
      CALL FFTGRIDPLAN(GRID_FOCK)

#ifdef _OPENACC
      CALL ACC_CUFFT_MAKEPLANS(GRIDHF,   1)
      CALL ACC_CUFFT_MAKEPLANS(GRID_FOCK,1)
#endif

#ifdef gammareal
      IF (.NOT.GRIDHF%LREAL) THEN
         CALL vtutor%bug("internal error in SETUP_FOCK: GRIDHF is not real \n the gamma only version &
            &in combination with Hartree-Fock requires \n to compile the code using -DNGXhalf or " &
            // "-DNGZhalf (mpiversion)", __FILE__, __LINE__)
      ENDIF
#endif

      ! set up the structure (FAST_AUG_FOCK) to perform the fast augmentation 
      ! of the charge density

      IF (WDES%LOVERL .AND. LMAX_FOCK>=0) THEN
         CALL SETUP_FASTAUG( T_INFO, P, WDES, GRIDHF, LATT_CUR, LMDIM, SZPRECFOCK, & 
              FAST_AUG_FOCK, TRANS_MATRIX_FOCK, & 
              LMAX_FOCK, LMAX_FOCKAE, NMAX_FOCKAE, QMAX_FOCKAE, EXXOEP, IU6, IU0 )

         CALL RSPHER(GRIDHF, FAST_AUG_FOCK, LATT_CUR)
#ifdef shmem_rproj
         IF (FAST_AUG_FOCK%COMM_shmem%NODE_ME==1) FAST_AUG_FOCK%RPROJ= FAST_AUG_FOCK%RPROJ*SQRT(LATT_CUR%OMEGA)
         CALLMPI( MPI_barrier( FAST_AUG_FOCK%COMM_shmem%MPI_COMM, IND ))
#else
         FAST_AUG_FOCK%RPROJ= FAST_AUG_FOCK%RPROJ*SQRT(LATT_CUR%OMEGA)
#endif
         IF (IU6>=0) WRITE(IU6,*) 'Maximum index for augmentation-charges in exchange ',MAXVAL(FAST_AUG_FOCK%NLIMAX)

         ! set up AUG_DES describing data layout of one center terms
         CALL SETUP_AUG_DES(GRIDHF, WDES_FOCK, AUG_DES, FAST_AUG_FOCK )
      ELSE
         ALLOCATE(AUG_DES%LMMAX(SIZE(WDES%LMMAX)))
         ! to make live easier set LMMAX (number of channels) to zero for each type
         ! also set NPRO, NPROD etc. to zero
         AUG_DES%LMMAX=0
         AUG_DES%NPRO =0
         AUG_DES%NPROD=0
         AUG_DES%NPRO_TOT=0
      ENDIF

!      CALL SET_FOCK_ALPHA(GRIDHF%NGPTAR(1),GRIDHF%NGPTAR(2),GRIDHF%NGPTAR(3),LATT_CUR)

      ! set the decay constant for the exponential used to
      ! calculate convergence corrections
      ! this new setup yields stable values even for PRECFOCK = F
      ! 
      IF (HFALPHA==-1) THEN
         HFALPHA=6/(WDES%ENMAX /RYTOEV)/(2*PI/AUTOA)
      ENDIF

      IF (MCALPHA<0) THEN
         !     ideally nothing should depend on the value here: we started with ***/200 and
         !     then converged on this - seems to work:
         !        MCALPHA=(WDES%ENMAX /RYTOEV)/12
         MCALPHA=(WDES%ENMAX / RYTOEV)/18.9
         !        used in MP2 routines so they know that mcalpha was set here      
         MCALPHA_DEFAULT=.TRUE.

         IF (IU6>=0) THEN
            WRITE(IU6,"( / &
     &       ' MCALPHA set automatically to default value (depends on ENMAX)',/ & 
     &       '   MCALPHA =',F10.4,' extent of test-charge in conv. correction in multipole expansion'/)") MCALPHA
         ENDIF

      ENDIF

      ! precalculate and store convergence corrections
      IF (HFRCUT<0) THEN
         HFRCUT=0
         SCALE=KPOINTS_FULL%WTKPT(1)*NKREDX*NKREDY*NKREDZ
         IF (ODDONLY .OR. EVENONLY) SCALE=SCALE*2
         HFRCUT=(LATT_CUR%OMEGA/SCALE/4/PI*3)**(1.0_q/3.0_q)
         IF (IU0>=0) WRITE(IU0,*) 'HFRCUT set to (new)',HFRCUT
         IF (IU6>=0) WRITE(IU6,*) 'HFRCUT set to (new)',HFRCUT
      ENDIF

      ALLOCATE (FSG_STORE(WDES%NKPTS))
      CALL SET_FSG_STORE(GRIDHF, LATT_CUR, WDES)
!      DO NK=1,WDES%NKPTS
!         FSG_STORE(NK)=SET_FSG(GRIDHF, LATT_CUR, NK)
!     ENDDO

      IF (MCALPHA/=0) THEN 
         CALL FOCK_MULTIPOLE_CORR_SETUP(LATT_CUR, GRIDHF)
         CALL FOCK_MULTIPOLE_PROJ_SETUP(LATT_CUR, GRIDHF)
      ENDIF

      ! setup wpbe interpolation tables
      CALL INIT_WPBE

      IF (IU6>=0) THEN
         WRITE(IU6,*) ' SETUP_FOCK is finished'
         CALL WFORCE(IU6)
      ENDIF

    END SUBROUTINE SETUP_FOCK

!************************ SUBROUTINE RESETUP_FOCK_WDES *****************
!
!> this routine resets the WDES_FOCK and must be called when
!> WDES has been updated (for instance if the number of k-points)
!> has been changed
!
!***********************************************************************

    SUBROUTINE RESETUP_FOCK_WDES(WDES, LATT_CUR, LATT_INI, IU6)
      USE wave
      USE lattice
      USE full_kpoints
      IMPLICIT NONE

      TYPE (wavedes),TARGET :: WDES       !< wave function descriptor
      TYPE (latt)    :: LATT_CUR          !< lattice parameters
      TYPE (latt)    :: LATT_INI          !< initial lattice parameters (determining the basis set)
      INTEGER IU6
    ! local
      LOGICAL, EXTERNAL :: CALCULATE_RESPONSE_FUNCTIONS,LDMATRIX
      INTEGER :: NK
    

      IF (.NOT. LHFCALC .AND. .NOT. CALCULATE_RESPONSE_FUNCTIONS() .AND. .NOT. LDMATRIX()) RETURN

      IF ( LWDES_FOCK ) THEN
         CALL DEALLOCWDES(WDES_FOCK,LEXTEND=.TRUE.)

         WDES_FOCK=WDES
         CALL GEN_LAYOUT(GRID_FOCK,WDES_FOCK, LATT_CUR%B, LATT_INI%B, IU6, .TRUE.)
         CALL GEN_INDEX (GRID_FOCK,WDES_FOCK, LATT_CUR%B, LATT_INI%B, -1, -1, .TRUE.)

         ! now check whether it is equivalent to original layout
         ! this is required since the FFTEXT and FFTWAV would not work otherwise
         CALL CHECK_GEN_LAYOUT( WDES,  WDES_FOCK, WDES%NKPTS)
         
         ! set plane wave index array INDPW between NKPTS and NKDIM
         CALL SET_INDPW_FOCK(GRID_FOCK,WDES_FOCK)
      ELSE
         GRID_FOCK=> WDES%GRID
         WDES_FOCK=> WDES
      ENDIF
      IF (WDES%LOVERL .AND. LMAX_FOCK>=0) THEN
         CALL DEALLOC_AUG_DES(AUG_DES)
         CALL SETUP_AUG_DES(GRIDHF, WDES_FOCK, AUG_DES, FAST_AUG_FOCK )
      ENDIF
      CALL SET_INDPW_FOCK(WDES%GRID, WDES)

      ! re-determined convergence corrections

      DEALLOCATE(FSG_STORE)
      ALLOCATE (FSG_STORE(WDES%NKPTS))
      CALL SET_FSG_STORE(GRIDHF, LATT_CUR, WDES)
!      DO NK=1,WDES%NKPTS
!         FSG_STORE(NK)=SET_FSG(GRIDHF, LATT_CUR, NK)
!      ENDDO

    END SUBROUTINE RESETUP_FOCK_WDES

!************************ SUBROUTINE RESETUP_FOCK     ******************
!
!> this routine resets the FAST_AUG_FOCK structure
!> must be called when the lattice vectors have been change
!
!***********************************************************************

    SUBROUTINE RESETUP_FOCK(WDES, LATT_CUR)
      USE wave
      USE lattice
      IMPLICIT NONE

      TYPE (wavedes) :: WDES       !< wave function descriptor
      TYPE (latt)    :: LATT_CUR   !< lattice parameters
      LOGICAL, EXTERNAL :: CALCULATE_RESPONSE_FUNCTIONS,LDMATRIX
      INTEGER :: NK
      LOGICAL :: LREALLOCATE

      IF ((LHFCALC .OR. CALCULATE_RESPONSE_FUNCTIONS() .OR. LDMATRIX()) .AND. WDES%LOVERL .AND. LMAX_FOCK>=0 ) THEN
         CALL REAL_OPTLAY(GRIDHF, LATT_CUR, FAST_AUG_FOCK, .TRUE., LREALLOCATE, -1,-1)
         IF (LREALLOCATE) THEN
            ! reallocate real space projectors
            CALL NONLR_DEALLOC(FAST_AUG_FOCK)
            CALL NONLR_ALLOC(FAST_AUG_FOCK)
         END IF

         CALL RSPHER(GRIDHF, FAST_AUG_FOCK, LATT_CUR)
#ifdef shmem_rproj
         IF (FAST_AUG_FOCK%COMM_shmem%NODE_ME==1) FAST_AUG_FOCK%RPROJ= FAST_AUG_FOCK%RPROJ*SQRT(LATT_CUR%OMEGA)
         CALLMPI( M_barrier( FAST_AUG_FOCK%COMM_shmem ))
#else
         FAST_AUG_FOCK%RPROJ= FAST_AUG_FOCK%RPROJ*SQRT(LATT_CUR%OMEGA)
#endif

         ! recalculate convergence correction
         CALL SET_FSG_STORE(GRIDHF, LATT_CUR, WDES) 
!         DO NK=1,WDES%NKPTS
!            FSG_STORE(NK)=SET_FSG(GRIDHF, LATT_CUR, NK)
!         ENDDO

      ENDIF

    END SUBROUTINE RESETUP_FOCK
 
!************************ SUBROUTINE SET_FOCK_KPOINTS ******************
!
!> this routine resets the NKDIM entry in the WDES array to the
!> total number of k-points in the entire BZ
!> in this way the WDES and NONL structure will be allocated to allow
!> a handling of all k-points in the entire BZ
!> this is in principle no longer required since the FOCK routine
!> now calculates the wave function character using symmetry
!> whereas previously it had to do it using the non local projector functions
!
!***********************************************************************

  SUBROUTINE SET_FOCK_KPOINTS(NKDIM)
    USE prec
    USE full_kpoints
    IMPLICIT NONE
    INTEGER NKDIM
    REAL(q), POINTER :: VKPT(:,:)

    IF (.NOT. LHFCALC) RETURN

    CALL CHECK_FULL_KPOINTS
    NKDIM=KPOINTS_FULL%NKPTS

  END SUBROUTINE SET_FOCK_KPOINTS

!********************** SUBROUTINE SET_INDPW_FOCK  *********************
!
!> SET_INDPW_FOCK sets the INDPW array for the wavefunctions belonging to
!> k-points that have symmetry-equivalents k-points in the IBZ
!>
!>  the following arrays in the WDES array are initialised as well:
!>   NGVECTOR
!>   IGX, IGY, IGZ
!>   NINDPW
!
!***********************************************************************

    SUBROUTINE SET_INDPW_FOCK(GRID, WDES )
      USE prec
      USE mgrid
      USE wave
      USE full_kpoints
      IMPLICIT NONE

      TYPE (grid_3d) :: GRID
      TYPE (wavedes) :: WDES

      IF (.NOT. LHFCALC) RETURN
      
      CALL CHECK_FULL_KPOINTS
      CALL SET_INDPW_FULL(GRID, WDES, KPOINTS_FULL)

    END SUBROUTINE SET_INDPW_FOCK



!********************** SUBROUTINE SET_GFAC_VKPT **********************
!
!> setup the (possibly screened)  Coloumb kernel
!>
!> 4 pi e^2 / (G+k-q)**2 * screening
!>
!> or if LCONJG is set
!>
!> 4 pi e^2 / (G+k+q)**2 * screening
!>
!> additionally set the zero grid point value to the convergence
!> correction FSG if NK==NQ (see below SET_FSG)
!>
!> for maximal performance
!> the routine includes a weighting factor
!> 1/NRPLWV
!> and furthermore includes a q-point
!> weighting factor
!>    KPOINTS_FULL%WTKPT(NQ)*NKREDX*NKREDY*NKREDZ
!>
!> whereas the routine  SET_GFAC_WITHOUT_WEIGHT does not include
!> a k-point weight
!>
!> see below for more details
! 
!**********************************************************************

    SUBROUTINE SET_GFAC_VKPT(GRID, LATT_CUR, VKPT, NQ, FSG, POTFAK, ENCUT)
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE constant
      USE mgrid
      USE nonl_high
      USE lattice
      USE full_kpoints
      IMPLICIT NONE

      TYPE (grid_3d) GRID
      TYPE (latt) LATT_CUR
      REAL(q) :: VKPT(3)           !< coordinate of k-point in reciprocal lattice units
      INTEGER :: NQ                !< k-point at which contribution to exchange kernel is evaluated
      REAL(q) :: FSG
      REAL(q), OPTIONAL :: ENCUT
      REAL(q) :: POTFAK(GRID%MPLWV)
   ! local
      INTEGER    NI,NC,N1,N2,N3
      REAL(q) :: DKX,DKY,DKZ,GX,GY,GZ,GSQU,GSQUP,GQUAD,SCALE, RHOEFF, FOMEGAP2

      PROFILING_START('fock_set_gfac_vkpt')

      CALL CHECK_FULL_KPOINTS

      ! effective electron density in a.u.^-3
      RHOEFF=(HFSCREEN*HFSCREEN*AUTOA*AUTOA/(4._q*EXP(LOG(3._q/PI)/3._q)))**3
      ! plasma frequency
      FOMEGAP2=16._q*PI*RHOEFF/(AUTOA*AUTOA*AUTOA*AUTOA)

      ! all k-points in the full grid have equal weight equal to first one
      SCALE=KPOINTS_FULL%WTKPT(NQ)*NKREDX*NKREDY*NKREDZ *EDEPS/LATT_CUR%OMEGA/TPI**2
      IF (ODDONLY .OR. EVENONLY ) SCALE=SCALE*2

      ! set correct phase in FAST_AUG structure
      CALL PHASER_HF(GRID, LATT_CUR, FAST_AUG_FOCK,VKPT)

      DKX=VKPT(1)*LATT_CUR%B(1,1)+VKPT(2)*LATT_CUR%B(1,2)+VKPT(3)*LATT_CUR%B(1,3)
      DKY=VKPT(1)*LATT_CUR%B(2,1)+VKPT(2)*LATT_CUR%B(2,2)+VKPT(3)*LATT_CUR%B(2,3)
      DKZ=VKPT(1)*LATT_CUR%B(3,1)+VKPT(2)*LATT_CUR%B(3,2)+VKPT(3)*LATT_CUR%B(3,3)

!$ACC PARALLEL LOOP COLLAPSE(force:2) PRESENT(GRID,LATT_CUR,KPOINTS_FULL,POTFAK) &
!$ACC PRIVATE(N2,N3,NI,GX,GY,GZ,GSQU,GSQUP,GQUAD) __IF_ASYNC__
      col: DO NC=1,GRID%RC%NCOL
         N2=GRID%RC%I2(NC)
         N3=GRID%RC%I3(NC)
         row: DO N1=1,GRID%RC%NROW
            NI=N1+(NC-1)*GRID%RC%NROW
            GX=(GRID%LPCTX(N1)*LATT_CUR%B(1,1)+GRID%LPCTY(N2)* &
                 LATT_CUR%B(1,2)+GRID%LPCTZ(N3)*LATT_CUR%B(1,3))
            GY=(GRID%LPCTX(N1)*LATT_CUR%B(2,1)+GRID%LPCTY(N2)* &
                 LATT_CUR%B(2,2)+GRID%LPCTZ(N3)*LATT_CUR%B(2,3))
            GZ=(GRID%LPCTX(N1)*LATT_CUR%B(3,1)+GRID%LPCTY(N2)* &
                 LATT_CUR%B(3,2)+GRID%LPCTZ(N3)*LATT_CUR%B(3,3))
            GSQU=(DKX+GX)**2+(DKY+GY)**2+(DKZ+GZ)**2
            IF (GSQU<KPOINTS_FULL%VKPTMINDIST2/40) THEN
               POTFAK(NI)=FSG
            ELSE
               IF (MODEL_GW>0) THEN
                  GSQUP=GSQU*(TPI*TPI)
                  GQUAD=GSQUP*GSQUP
                  ! corrected by gK 04.04.2014 (seems MODEL_GW_ALPHA is usually set to 1.0_q)
                  ! although the Bechstedt et al. Solid State Comm. 84, 765 suggests to use 1.5_q
                  POTFAK(NI)=SCALE/(GSQU)/ &
                       (1+1/(1/(MODEL_GW_EPS0-1)+MODEL_GW_ALPHA*GSQUP/(HFSCREEN*HFSCREEN)+GQUAD/FOMEGAP2))
               ELSE IF (HFRCUT/=0) THEN
                  ! spherical cutoff on Coloumb kernel
                  ! see for instance C.A. Rozzi, PRB 73, 205119 (2006)
                  POTFAK(NI)=SCALE/(GSQU)*(1-COS(SQRT(GSQU)*TPI*HFRCUT))
               ELSE IF (HFSCREEN==0) THEN
                  ! Coloumb kernel
                  POTFAK(NI)=SCALE/(GSQU)
               ELSE IF (L_THOMAS_FERMI) THEN
                  ! screened Thomas Fermi kernel
                  POTFAK(NI)=SCALE/(GSQU+HFSCREEN*HFSCREEN/(TPI*TPI))
               ELSE IF (LRHFCALC) THEN
                  ! exponentially screened kernel corresponding to erfc( HFSCREEN r)/r 
                  POTFAK(NI)=SCALE/GSQU*EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN)))
               ELSE IF (L_MODEL_HF) THEN
                  ! AEXX * case(HFSCREEN==0) + (1-AEXX)* default
                  POTFAK(NI)=SCALE/GSQU*(1-(1-AEXX)*EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN))))
               ELSE
                  POTFAK(NI)=SCALE/GSQU*(1-EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN))))
               ENDIF
            ENDIF
            IF (PRESENT(ENCUT)) THEN
               IF (HSQDTM*TPI*TPI*GSQU>ENCUT) THEN
                  POTFAK(NI)=0
               ENDIF
            ENDIF
            IF (ENCUTFOCK>0) THEN
               IF (HSQDTM*TPI*TPI*GSQU>ENCUTFOCK) THEN
                  POTFAK(NI)=0
               ENDIF
            ENDIF
         ENDDO row
      ENDDO col

      CALL SETUNB_REAL(POTFAK,GRID)

      IF (L_MODEL_HF) THEN
!$ACC PARALLEL LOOP PRESENT(GRID,POTFAK) __IF_ASYNC__
         DO NC=1,GRID%RC%NP
            POTFAK(NC)=POTFAK(NC)*(1.0_q/GRID%NPLWV)
         ENDDO
      ELSE

!$ACC PARALLEL LOOP PRESENT(GRID,POTFAK) __IF_ASYNC__
         DO NC=1,GRID%RC%NP
            POTFAK(NC)=POTFAK(NC)*(AEXX/GRID%NPLWV)
         ENDDO
      ENDIF

      PROFILING_STOP('fock_set_gfac_vkpt')

    END SUBROUTINE SET_GFAC_VKPT

!
!> standard calling interface to the SET_GFAC routine
!> requires one to pass down two k-points for which the Coulomb kernel is required
!> it is assumed that the exchange operator is acting on an orbitals at psi_k
!> \sum_q w_q phi_q(r)   phi^*_q(r') psi_k(r') / | r-r'|
! 
    SUBROUTINE SET_GFAC(GRID, LATT_CUR, NK, NQ, FSG, POTFAK, LCONJG, ENCUT)
      USE lattice
      USE full_kpoints
      IMPLICIT NONE

      TYPE (grid_3d) GRID
      TYPE (latt) LATT_CUR
      INTEGER NK, NQ
      REAL(q) :: FSG
      LOGICAL, OPTIONAL :: LCONJG
      REAL(q), OPTIONAL :: ENCUT
      REAL(q) :: POTFAK(GRID%MPLWV)

      IF (PRESENT(ENCUT).AND. PRESENT(LCONJG)) THEN
         CALL SET_GFAC_VKPT(GRID, LATT_CUR, KPOINTS_FULL%VKPT(:,NK)+KPOINTS_FULL%VKPT(:,NQ), NQ, FSG, POTFAK, ENCUT)
      ELSE IF (PRESENT(ENCUT)) THEN
         CALL SET_GFAC_VKPT(GRID, LATT_CUR, KPOINTS_FULL%VKPT(:,NK)-KPOINTS_FULL%VKPT(:,NQ), NQ, FSG, POTFAK, ENCUT)
      ELSE IF (PRESENT(LCONJG)) THEN
         CALL SET_GFAC_VKPT(GRID, LATT_CUR, KPOINTS_FULL%VKPT(:,NK)+KPOINTS_FULL%VKPT(:,NQ), NQ, FSG, POTFAK)
      ELSE
         CALL SET_GFAC_VKPT(GRID, LATT_CUR, KPOINTS_FULL%VKPT(:,NK)-KPOINTS_FULL%VKPT(:,NQ), NQ, FSG, POTFAK)
      ENDIF

    END SUBROUTINE SET_GFAC

    SUBROUTINE SET_GFAC_WITHOUT_WEIGHT(GRID, LATT_CUR, NK, NQ, FSG, POTFAK, LCONJG)
      USE constant
      USE mgrid
      USE lattice
      USE full_kpoints
      IMPLICIT NONE

      TYPE (grid_3d) GRID
      TYPE (latt) LATT_CUR
      INTEGER NK, NQ
      REAL(q) :: FSG
      REAL(q) :: POTFAK(GRID%MPLWV)
      LOGICAL, OPTIONAL :: LCONJG
      ! local
      INTEGER NC
      REAL(q) SCALE

      IF (PRESENT(LCONJG)) THEN
         CALL SET_GFAC(GRID, LATT_CUR, NK, NQ, FSG, POTFAK, LCONJG)
      ELSE
         CALL SET_GFAC(GRID, LATT_CUR, NK, NQ, FSG, POTFAK)
      ENDIF

      SCALE=1/(KPOINTS_FULL%WTKPT(NQ)*NKREDX*NKREDY*NKREDZ)
      IF (ODDONLY .OR. EVENONLY ) SCALE=SCALE/2

      DO NC=1,GRID%RC%NP
         POTFAK(NC)=POTFAK(NC)*SCALE
      ENDDO

    END SUBROUTINE SET_GFAC_WITHOUT_WEIGHT

!********************** SUBROUTINE SET_GFAC_WAVEFUN *******************
!
!> setup the (possibly screened)  Coloumb kernel
!>
!> this version calculates the kernel within the spherical cutoff
!> sphere defined by the plane wave cutoff
!> it is required in the context of GW routines and BSE routines
! 
!**********************************************************************

    SUBROUTINE SET_GFAC_WAVEFUN(WDES1, LATT_CUR, FSG, POTFAK, ENCUT, ENCUTSOFT)
      USE constant
      USE wave
      USE lattice
      USE full_kpoints
      IMPLICIT NONE

      TYPE (wavedes1) WDES1   !< descriptor for a specific k-point
      TYPE (latt) :: LATT_CUR 
      REAL(q)     :: FSG
      REAL(q) :: POTFAK(WDES1%NGVECTOR)
      REAL(q), OPTIONAL :: ENCUT, ENCUTSOFT
   ! local
      INTEGER    NI,NP
      REAL(q) :: DKX,DKY,DKZ,GX,GY,GZ,GSQU,GSQUP,GQUAD,SCALE, SCALEFSG, RHOEFF, FOMEGAP2
      REAL(q) :: E

      PROFILING_START('fock_set_gfac_wavefun')

      ! effective electron density in a.u.^-3
      RHOEFF=(HFSCREEN*HFSCREEN*AUTOA*AUTOA/(4._q*EXP(LOG(3._q/PI)/3._q)))**3
      ! plasma frequency
      FOMEGAP2=16._q*PI*RHOEFF/(AUTOA*AUTOA*AUTOA*AUTOA)

      SCALE=EDEPS/LATT_CUR%OMEGA/TPI**2

      IF (KPOINTS_FULL%WTKPT(1)==0) THEN
         CALL vtutor%bug("internal error in SET_GFAC_WAVEFUN: division by zero in SCALEFSG, and FSG &
            &anyway most likely wrong \n needs to be fixed in CALCULATE_LOCAL_FIELD_PREPARE as well", &
            __FILE__, __LINE__)
      ENDIF
      
      SCALEFSG=1/(KPOINTS_FULL%WTKPT(1)*NKREDX*NKREDY*NKREDZ)
      IF (ODDONLY .OR. EVENONLY ) SCALEFSG=SCALEFSG/2

      DKX=WDES1%VKPT(1)*LATT_CUR%B(1,1)+ &
          WDES1%VKPT(2)*LATT_CUR%B(1,2)+ &
          WDES1%VKPT(3)*LATT_CUR%B(1,3)
      DKY=WDES1%VKPT(1)*LATT_CUR%B(2,1)+ &
          WDES1%VKPT(2)*LATT_CUR%B(2,2)+ &
          WDES1%VKPT(3)*LATT_CUR%B(2,3)
      DKZ=WDES1%VKPT(1)*LATT_CUR%B(3,1)+ &
          WDES1%VKPT(2)*LATT_CUR%B(3,2)+ &
          WDES1%VKPT(3)*LATT_CUR%B(3,3)

      NP=WDES1%NGVECTOR

      DO NI=1,NP
       
         GX=(WDES1%IGX(NI)*LATT_CUR%B(1,1)+WDES1%IGY(NI)* &
              LATT_CUR%B(1,2)+WDES1%IGZ(NI)*LATT_CUR%B(1,3))
         GY=(WDES1%IGX(NI)*LATT_CUR%B(2,1)+WDES1%IGY(NI)* &
              LATT_CUR%B(2,2)+WDES1%IGZ(NI)*LATT_CUR%B(2,3))
         GZ=(WDES1%IGX(NI)*LATT_CUR%B(3,1)+WDES1%IGY(NI)* &
              LATT_CUR%B(3,2)+WDES1%IGZ(NI)*LATT_CUR%B(3,3))
         
         GSQU=(DKX+GX)**2+(DKY+GY)**2+(DKZ+GZ)**2
         IF (GSQU<KPOINTS_FULL%VKPTMINDIST2/40) THEN
            POTFAK(NI)=FSG*SCALEFSG
         ELSE
            IF (MODEL_GW>0) THEN
               GSQUP=GSQU*(TPI*TPI)
               GQUAD=GSQUP*GSQUP
               POTFAK(NI)=SCALE/(GSQU)/ &
                    (1+1/(1/(MODEL_GW_EPS0-1)+MODEL_GW_ALPHA*GSQUP/(HFSCREEN*HFSCREEN)+GQUAD/FOMEGAP2))
            ELSE IF (HFRCUT/=0) THEN
               POTFAK(NI)=SCALE/(GSQU)*(1-COS(SQRT(GSQU)*TPI*HFRCUT))
            ELSE IF (HFSCREEN==0) THEN
               POTFAK(NI)=SCALE/(GSQU)
            ELSE IF (L_THOMAS_FERMI) THEN
               POTFAK(NI)=SCALE/(GSQU+HFSCREEN*HFSCREEN/(TPI*TPI))
            ELSE IF (LRHFCALC) THEN
               POTFAK(NI)=SCALE/GSQU*EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN)))
            ELSE IF (L_MODEL_HF) THEN
               ! AEXX * case(HFSCREEN==0) + (1-AEXX)* default
               POTFAK(NI)=SCALE/GSQU*(1-(1-AEXX)*EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN))))
            ELSE
               POTFAK(NI)=SCALE/GSQU*(1-EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN))))
            ENDIF
         ENDIF

         ! smooth cutoff function between  ENCUTSOFT and ENCUT
         E=HSQDTM*(GSQU*TPI**2)
         IF (PRESENT(ENCUT).AND. PRESENT(ENCUTSOFT)) THEN 
            IF (E>ENCUTSOFT) THEN
               POTFAK(NI)=POTFAK(NI)*(1+COS((E-ENCUTSOFT)/(ENCUT-ENCUTSOFT)*PI))/2
            ENDIF
         ENDIF
      ENDDO

      IF (L_MODEL_HF) THEN
         DO NI=1,NP
            POTFAK(NI)=POTFAK(NI)*(1.0_q/WDES1%GRID%NPLWV)
         ENDDO
      ELSE
         DO NI=1,NP
            POTFAK(NI)=POTFAK(NI)*(AEXX/WDES1%GRID%NPLWV)
         ENDDO
      ENDIF

      PROFILING_STOP('fock_set_gfac_wavefun')

    END SUBROUTINE SET_GFAC_WAVEFUN


!********************** SUBROUTINE SET_FSG_QDER ************************
!
!> calculate derivative of kernel with respect to q
!> can be done easily anyltically, but really no time to dvelve into
!> this for all functionals
!>
!> open issues: the derivative of FSG (the singularity correction)
!> w.r.t k is not included presently
!> I think it should be zero anyhow
!
!
!**********************************************************************

    SUBROUTINE SET_GFAC_QDER(GRID, LATT_CUR, NK, NQ, FSG, POTFAK, IDIR)
      USE constant
      USE mgrid
      USE nonl_high
      USE lattice
      USE full_kpoints
      IMPLICIT NONE

      TYPE (grid_3d) GRID
      TYPE (latt) LATT_CUR
      INTEGER NK, NQ
      REAL(q) :: FSG
      REAL(q) :: POTFAK(GRID%MPLWV)
      INTEGER :: IDIR                  !< direction
    ! local
      REAL(q) :: POTFAK_(GRID%MPLWV)
      REAL(q) :: DISPL(3)
      REAL(q), PARAMETER :: DIS=fd_displacement
      
      ! shift k-point slightly
      DISPL=0
      DISPL(IDIR)=DIS/TPI
      CALL KARDIR(1, DISPL(1), LATT_CUR%A)

      CALL SET_GFAC_VKPT(GRID, LATT_CUR, KPOINTS_FULL%VKPT(:,NK)-KPOINTS_FULL%VKPT(:,NQ)-DISPL, NQ, FSG, POTFAK_)
      CALL SET_GFAC_VKPT(GRID, LATT_CUR, KPOINTS_FULL%VKPT(:,NK)-KPOINTS_FULL%VKPT(:,NQ)+DISPL, NQ, FSG, POTFAK)

      POTFAK(1:GRID%RC%NP)=(POTFAK(1:GRID%RC%NP)-POTFAK_(1:GRID%RC%NP))/DIS/2

    END SUBROUTINE SET_GFAC_QDER


!********************** SUBROUTINE SET_FOCK_ALPHA ********************
!
!> the Posternak, Baldereschi convergence correction is defined
!> as the energy difference between a single exponential charge
!>    exp(-HFALPHA*(G+k)^2*2 pi 2)
!> and periodically repeated charges which are located
!> on a lattice reciprocal to the generating  k-point set
!> ideally we would like to have rather contracted charges
!> but there are limits imposed by the FFT grid
!
!**********************************************************************

   SUBROUTINE SET_FOCK_ALPHA(NGX,NGY,NGZ,LATT_CUR)
      USE prec
      USE lattice
      IMPLICIT NONE
      
      INTEGER NGX,NGY,NGZ
      TYPE (latt) LATT_CUR
      REAL(q) GSQU

      IF (HFALPHA==-1) THEN
         ! the minimum value for ALPHA is determine by the FFT grid
         ! seek shortest wave vector
         ! NGn * b(n) * a(n)/||a(n)||
         ! updated 24.10.2009 NGX/2-> NGX and summation in calculation of FSG
         ! from -NGX to +NGX
         GSQU=MIN(NGX/2/LATT_CUR%ANORM(1), &
                  NGY/2/LATT_CUR%ANORM(2), &
                  NGZ/2/LATT_CUR%ANORM(3))**2

         ! the larger GSQU is the smaller alpha can be
         ! corresponding to a more contracted charge
         ! the default below leads to errors of 0.5 meV
         ! compared to including an infinite number of G vector
         ! HFALPHA=1/GSQU/5 gK updated on 09.10.2008
         ! larger alpha results in more screened potential at large
         ! distances, this speeds up convergence
         HFALPHA=1/GSQU/2
      ENDIF

    END SUBROUTINE SET_FOCK_ALPHA

!********************** SUBROUTINE SET_FSG  ***************************
!
!> integratable singularity of exchange term:
!> calculate the average electrostatic potential prefactor divided
!> by the number of k-points
!> essentially
!>  4 pi e^2 / G^2 / NKPTS
!>
!> for the orbitals with k=k' and n=n' (i.e. those that are affected
!>  by divergence of the HF term in periodic boundary conditions)
!> this is the convergence accelerator suggested by F. Giggy and Baldereschi
!>
!>   sum_G,k epsilon/omega/pi^2 exp(-HFALPHA*(G+k)^2*2 pi 2)/ (G+k)^2
!>
!> if the generation scheme of the k-point mesh is known use a simple
!> monopole correction
!
!**********************************************************************

    FUNCTION SET_FSG(GRID, LATT_CUR, NK_ORIG)
      USE prec
      USE constant
      USE mgrid
      USE nonl_high
      USE lattice
      USE full_kpoints
      IMPLICIT NONE
      
      TYPE (grid_3d) GRID
      TYPE (latt) LATT_CUR
      INTEGER, OPTIONAL :: NK_ORIG
    ! local
      INTEGER NK,NQ
      REAL(q) :: SET_FSG
      REAL(q) :: HFALPHAP
      INTEGER    NC,N1,N2,N3, NN1, NN2, NN3, LPX, LPY, LPZ
      REAL(q) :: DKX,DKY,DKZ,GX,GY,GZ,GSQU,SCALE,FSG,FACTM,GX0,GY0,GZ0
      TYPE (latt) LATT_EWALD

      PROFILING_START('fock_set_fsg')

      CALL CHECK_FULL_KPOINTS
!
! spherical cutoff on Coloumb kernel
! see for instance C.A. Rozzi, PRB 73, 205119 (2006)
!
      IF (HFRCUT/=0) THEN
         SCALE=KPOINTS_FULL%WTKPT(1)*NKREDX*NKREDY*NKREDZ
         IF (ODDONLY .OR. EVENONLY) SCALE=SCALE*2
         FSG=HFRCUT*HFRCUT/2*EDEPS/LATT_CUR%OMEGA*SCALE
!
! simple version, that can be applied if
! a Monkhorst Pack or Gamma centered grid is used
      ELSE IF (KPOINTS_FULL%NKPTS_NON_ZERO==KPOINTS_FULL%NKPX*KPOINTS_FULL%NKPY*KPOINTS_FULL%NKPZ &
           .AND. (HFSCREEN==0.OR. MODEL_GW>0 ).AND. .NOT. ODDONLY .AND. .NOT. EVENONLY .AND. .NOT. HFKIDENT) THEN
         LATT_EWALD=LATT_CUR
         LATT_EWALD%A(:,1)=LATT_EWALD%A(:,1)*KPOINTS_FULL%NKPX/NKREDX
         LATT_EWALD%A(:,2)=LATT_EWALD%A(:,2)*KPOINTS_FULL%NKPY/NKREDY
         LATT_EWALD%A(:,3)=LATT_EWALD%A(:,3)*KPOINTS_FULL%NKPZ/NKREDZ
         CALL LATTIC(LATT_EWALD)
         CALL EWALD_MONO(FSG,1.0_q,LATT_EWALD)
         FSG=FSG*2
         IF (MODEL_GW>0) FSG=FSG/MODEL_GW_EPS0
      ELSE
!
! Giggy, Posternak, Baldereschi version
!
         IF (PRESENT(NK_ORIG)) THEN
            NK=NK_ORIG
         ELSE
            NK=1
         ENDIF

         FSG=0
         DO NQ=1,KPOINTS_FULL%NKPTS
            IF (KPOINTS_FULL%WTKPT(NQ)==0) CYCLE
            SCALE=KPOINTS_FULL%WTKPT(NQ)*NKREDX*NKREDY*NKREDZ *EDEPS/LATT_CUR%OMEGA/TPI**2
            IF (ODDONLY .OR. EVENONLY) SCALE=SCALE*2
            IF (SKIP_THIS_KPOINT_IN_FOCK(KPOINTS_FULL%VKPT(:,NQ)-KPOINTS_FULL%VKPT(:,NK))) CYCLE

#ifndef gammareal
            IF (.NOT. HFKIDENT) THEN
            DKX=(KPOINTS_FULL%VKPT(1,NK)-KPOINTS_FULL%VKPT(1,NQ))*LATT_CUR%B(1,1)+ &
                 (KPOINTS_FULL%VKPT(2,NK)-KPOINTS_FULL%VKPT(2,NQ))*LATT_CUR%B(1,2)+ &
                 (KPOINTS_FULL%VKPT(3,NK)-KPOINTS_FULL%VKPT(3,NQ))*LATT_CUR%B(1,3)
            DKY=(KPOINTS_FULL%VKPT(1,NK)-KPOINTS_FULL%VKPT(1,NQ))*LATT_CUR%B(2,1)+ &
                 (KPOINTS_FULL%VKPT(2,NK)-KPOINTS_FULL%VKPT(2,NQ))*LATT_CUR%B(2,2)+ &
                 (KPOINTS_FULL%VKPT(3,NK)-KPOINTS_FULL%VKPT(3,NQ))*LATT_CUR%B(2,3)
            DKZ=(KPOINTS_FULL%VKPT(1,NK)-KPOINTS_FULL%VKPT(1,NQ))*LATT_CUR%B(3,1)+ &
                 (KPOINTS_FULL%VKPT(2,NK)-KPOINTS_FULL%VKPT(2,NQ))*LATT_CUR%B(3,2)+ &
                 (KPOINTS_FULL%VKPT(3,NK)-KPOINTS_FULL%VKPT(3,NQ))*LATT_CUR%B(3,3)
            ELSE
            DKX=( -KPOINTS_FULL%VKPT(1,NQ))*LATT_CUR%B(1,1)+ &
                 (-KPOINTS_FULL%VKPT(2,NQ))*LATT_CUR%B(1,2)+ &
                 (-KPOINTS_FULL%VKPT(3,NQ))*LATT_CUR%B(1,3)
            DKY=( -KPOINTS_FULL%VKPT(1,NQ))*LATT_CUR%B(2,1)+ &
                 (-KPOINTS_FULL%VKPT(2,NQ))*LATT_CUR%B(2,2)+ &
                 (-KPOINTS_FULL%VKPT(3,NQ))*LATT_CUR%B(2,3)
            DKZ=( -KPOINTS_FULL%VKPT(1,NQ))*LATT_CUR%B(3,1)+ &
                 (-KPOINTS_FULL%VKPT(2,NQ))*LATT_CUR%B(3,2)+ &
                 (-KPOINTS_FULL%VKPT(3,NQ))*LATT_CUR%B(3,3)
            ENDIF
#else
            DKX=0._q; DKY=0._q; DKZ=0._q
#endif
            DO NN1=0,GRID%NGX,GRID%NGX
            DO NN2=0,GRID%NGY,GRID%NGY
            DO NN3=0,GRID%NGZ,GRID%NGZ
            col: DO NC=1,GRID%RC%NCOL
               N2=GRID%RC%I2(NC)
               N3=GRID%RC%I3(NC)
               LPY=MOD(GRID%LPCTY(N2)+NN2+GRID%NGY-1,2*GRID%NGY)-GRID%NGY+1
               LPZ=MOD(GRID%LPCTZ(N3)+NN3+GRID%NGZ-1,2*GRID%NGZ)-GRID%NGZ+1
               GX0=(LPY*LATT_CUR%B(1,2)+LPZ*LATT_CUR%B(1,3))+DKX
               GY0=(LPY*LATT_CUR%B(2,2)+LPZ*LATT_CUR%B(2,3))+DKY
               GZ0=(LPY*LATT_CUR%B(3,2)+LPZ*LATT_CUR%B(3,3))+DKZ

               row: DO N1=1,GRID%RC%NROW
                  LPX=MOD(GRID%LPCTX(N1)+NN1+GRID%NGX-1,2*GRID%NGX)-GRID%NGX+1

                  GX=(LPX*LATT_CUR%B(1,1)+GX0)
                  GY=(LPX*LATT_CUR%B(2,1)+GY0)
                  GZ=(LPX*LATT_CUR%B(3,1)+GZ0)
                  GSQU=(GX)**2+(GY)**2+(GZ)**2

                  FACTM=1._q
                  ! IF (GRID%LREAL .AND. (N1/=1 .AND. (N1/=GRID%RC%NROW))) FACTM=2._q
                  IF ( GRID%NGX /= GRID%NGX_rd .AND. (N1/=1 .AND. (N1/= GRID%NGX_rd))) FACTM=2._q
                  IF ( GRID%NGZ /= GRID%NGZ_rd .AND. (N3/=1 .AND. (N3/= GRID%NGZ_rd))) FACTM=2._q

                  IF (GSQU < KPOINTS_FULL%VKPTMINDIST2/40) THEN
                  ELSE
                     IF (MODEL_GW>0) THEN
                        FSG=FSG+FACTM*SCALE*EXP(-GSQU*(HFALPHA*TPI*TPI))/GSQU/MODEL_GW_EPS0
                     ELSE IF (HFSCREEN==0) THEN
                        ! Coulomb kernel
                        FSG=FSG+FACTM*SCALE*EXP(-GSQU*(HFALPHA*TPI*TPI))/GSQU
                     ELSE IF (L_THOMAS_FERMI) THEN
                        FSG=FSG+FACTM*SCALE*EXP(-GSQU*(HFALPHA*TPI*TPI))/(GSQU+HFSCREEN*HFSCREEN*(1.0_q/(TPI*TPI)))
                     ELSE IF (LRHFCALC) THEN
                        ! long range Coulomb kernel
                        FSG=FSG+FACTM*SCALE*EXP(GSQU*((-HFALPHA-1/(4*HFSCREEN*HFSCREEN))*TPI*TPI))/GSQU
                     ELSE IF (L_MODEL_HF) THEN
                        FSG=FSG+FACTM*SCALE/GSQU*EXP(-GSQU*(HFALPHA*TPI*TPI)) &
                               *(1-(1-AEXX)*EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN))))
                     ELSE
                        ! screened Coulomb kernel (short range)
                        FSG=FSG+FACTM*SCALE/GSQU*EXP(-GSQU*(HFALPHA*TPI*TPI)) &
                           *(1-EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN))))
                     ENDIF
                  ENDIF
               ENDDO row
            ENDDO col
            ENDDO
            ENDDO
            ENDDO
         ENDDO

         CALLMPI (M_sum_d(GRID%COMM,FSG,1))

         IF (MODEL_GW>0) THEN
            ! add selfenergy of an isolated erfc(r)  charge distribution / epsilon
            FSG=-FSG+EDEPS/(TPI**3)*2._q*HFALPHA*(PI/HFALPHA)**(3._q/2._q)/MODEL_GW_EPS0
         ELSE IF (HFSCREEN==0) THEN
            ! add selfenergy of an isolated erfc(r)  charge distribution
            FSG=-FSG+EDEPS/(TPI**3)*2._q*HFALPHA*(PI/HFALPHA)**(3._q/2._q)
         ELSE IF (L_THOMAS_FERMI) THEN
            FSG=0
         ELSE IF (LRHFCALC) THEN
            HFALPHAP=HFALPHA+1/HFSCREEN/HFSCREEN/4
            !
            ! we need the selfenergy of an erfc(r) charge distribution 
            ! for the long range Coulomb kernel
            !
            FSG=-FSG+EDEPS/(TPI**3)*2._q*HFALPHAP*(PI/HFALPHAP)**(3._q/2._q)
         ELSE IF (L_MODEL_HF) THEN
            !
            ! AEXX * case(HFSCREEN==0) + (1-AEXX)* default
            !
            HFALPHAP=HFALPHA+1/HFSCREEN/HFSCREEN/4
            FSG=-FSG+EDEPS/(TPI**3)*2._q*HFALPHA *(PI/HFALPHA )**(3._q/2._q)- &
                  (1-AEXX)*EDEPS/(TPI**3)*2._q*HFALPHAP*(PI/HFALPHAP)**(3._q/2._q)
         ELSE 
            HFALPHAP=HFALPHA+1/HFSCREEN/HFSCREEN/4
            !
            ! we need the selfenergy of an erfc(r) charge distribution in
            ! a screened short range Coloumb potential 
            !
            FSG=-FSG+EDEPS/(TPI**3)*2._q*HFALPHA *(PI/HFALPHA )**(3._q/2._q)- &
                 EDEPS/(TPI**3)*2._q*HFALPHAP*(PI/HFALPHAP)**(3._q/2._q)
         ENDIF
      ENDIF
      SET_FSG=FSG

      PROFILING_STOP('fock_set_fsg')

    END FUNCTION SET_FSG

!********************** SUBROUTINE SET_FSG_STORE **********************
!
! set FSG_STORE for all k-points calling SET_FSG
! this version parallelizes over k-points 
! depending on the selected k-point grid, it might also chose
! to calculate FSG only at the first k-point, and set the convergence
! correction for other k-points to the same value
!
!**********************************************************************
   
    SUBROUTINE SET_FSG_STORE(GRID, LATT_CUR, WDES)
      USE wave
      USE mgrid
      USE lattice
      USE full_kpoints
      IMPLICIT NONE
      
      TYPE (grid_3d) GRID
      TYPE (latt) LATT_CUR
      TYPE (wavedes),TARGET :: WDES       ! wave function descriptor
    ! local
      INTEGER NK

      FSG_STORE=0

      ! uniform grid, and no use of HFKIDENT, then FSG is identical for all k-points
      IF (KPOINTS_FULL%Mode>= MonkhorstPack .AND. KPOINTS_FULL%Mode <= GeneratingLattice .AND. &
          .NOT. HFKIDENT )  THEN
         FSG_STORE=SET_FSG(GRIDHF, LATT_CUR, 1)
      ELSE
      ! otherwise, one needs to calculate FSG for each k-point
         DO NK=1,WDES%NKPTS
         ! distribute round robin along KINTER and INTER
         ! no distribution within bands since that is utilized anyway in the
         ! GRID communicator
#ifdef MPI
            IF (MOD(NK-1,WDES%COMM_KINTER%NCPU*WDES%COMM_INTER%NCPU)== &
                (WDES%COMM_KINTER%NODE_ME-1)*WDES%COMM_INTER%NCPU+WDES%COMM_INTER%NODE_ME-1 ) THEN
              FSG_STORE(NK)=SET_FSG(GRIDHF, LATT_CUR, NK)
            ENDIF
#else
              FSG_STORE(NK)=SET_FSG(GRIDHF, LATT_CUR, NK)
#endif
         ENDDO

         CALLMPI( M_sum_d( WDES%COMM_KINTER, FSG_STORE , SIZE(FSG_STORE)))
         CALLMPI( M_sum_d( WDES%COMM_INTER, FSG_STORE , SIZE(FSG_STORE)))
      ENDIF

   END SUBROUTINE


!********************** SUBROUTINE SET_GFAC_DER ***********************
!
!> setup the derivative of the (possibly screened)
!> Coloumb kernel  e^2 /4 pi (G+k-q)**2 * (1-exp(  (G+k-q)**2/4/alpha^2))
!> corresponding the erfc(alpha*r)/r
!> additionally set the zero grid point value to FSG
!> if NK==NQ (see below SET_FSG)
! 
!**********************************************************************

    SUBROUTINE SET_GFAC_DER(GRID, LATT_CUR, NK, NQ, FSG, POTFAK)
      USE prec
      USE constant
      USE mgrid
      USE nonl_high
      USE lattice
      USE full_kpoints
      IMPLICIT NONE
      
      TYPE (grid_3d) GRID
      TYPE (latt) LATT_CUR
      INTEGER NK,NQ
      REAL(q) :: FSG(0:6)
      REAL(q) :: POTFAK(GRID%MPLWV,0:6)
   
      INTEGER    NI,NC,N1,N2,N3,I
      REAL(q) :: DKX,DKY,DKZ,GX,GY,GZ,GSQU,GSQU2,SCALE,EXP_GSQU_DER,FACTM

      PROFILING_START('fock_set_gfac_der')

      CALL CHECK_FULL_KPOINTS

      SCALE=KPOINTS_FULL%WTKPT(NQ)*NKREDX*NKREDY*NKREDZ *EDEPS/LATT_CUR%OMEGA/TPI**2
      IF (ODDONLY .OR. EVENONLY) SCALE=SCALE*2
#ifndef gammareal
      DKX=(KPOINTS_FULL%VKPT(1,NK)-KPOINTS_FULL%VKPT(1,NQ))*LATT_CUR%B(1,1)+ &
           (KPOINTS_FULL%VKPT(2,NK)-KPOINTS_FULL%VKPT(2,NQ))*LATT_CUR%B(1,2)+ &
           (KPOINTS_FULL%VKPT(3,NK)-KPOINTS_FULL%VKPT(3,NQ))*LATT_CUR%B(1,3)
      DKY=(KPOINTS_FULL%VKPT(1,NK)-KPOINTS_FULL%VKPT(1,NQ))*LATT_CUR%B(2,1)+ &
           (KPOINTS_FULL%VKPT(2,NK)-KPOINTS_FULL%VKPT(2,NQ))*LATT_CUR%B(2,2)+ &
           (KPOINTS_FULL%VKPT(3,NK)-KPOINTS_FULL%VKPT(3,NQ))*LATT_CUR%B(2,3)
      DKZ=(KPOINTS_FULL%VKPT(1,NK)-KPOINTS_FULL%VKPT(1,NQ))*LATT_CUR%B(3,1)+ &
           (KPOINTS_FULL%VKPT(2,NK)-KPOINTS_FULL%VKPT(2,NQ))*LATT_CUR%B(3,2)+ &
           (KPOINTS_FULL%VKPT(3,NK)-KPOINTS_FULL%VKPT(3,NQ))*LATT_CUR%B(3,3)
! set correct phase in FAST_AUG structure
      CALL PHASER_HF(GRID, LATT_CUR, FAST_AUG_FOCK,KPOINTS_FULL%VKPT(:,NK)-KPOINTS_FULL%VKPT(:,NQ))
#else
      DKX=0._q; DKY=0._q; DKZ=0._q
#endif

!$ACC PARALLEL LOOP COLLAPSE(force:2) PRESENT(GRID,LATT_CUR,POTFAK) &
!$ACC PRIVATE(N2,N3,NI,GX,GY,GZ,GSQU,GSQU2,FACTM,EXP_GSQU_DER) __IF_ASYNC__
      col: DO NC=1,GRID%RC%NCOL
         N2=GRID%RC%I2(NC)
         N3=GRID%RC%I3(NC)
         row: DO N1=1,GRID%RC%NROW
            NI=N1+(NC-1)*GRID%RC%NROW
            GX=(GRID%LPCTX(N1)*LATT_CUR%B(1,1)+GRID%LPCTY(N2)* &
                 LATT_CUR%B(1,2)+GRID%LPCTZ(N3)*LATT_CUR%B(1,3))+DKX
            GY=(GRID%LPCTX(N1)*LATT_CUR%B(2,1)+GRID%LPCTY(N2)* &
                 LATT_CUR%B(2,2)+GRID%LPCTZ(N3)*LATT_CUR%B(2,3))+DKY
            GZ=(GRID%LPCTX(N1)*LATT_CUR%B(3,1)+GRID%LPCTY(N2)* &
                 LATT_CUR%B(3,2)+GRID%LPCTZ(N3)*LATT_CUR%B(3,3))+DKZ
            GSQU=GX**2+GY**2+GZ**2
            GSQU2=GSQU*GSQU

            FACTM=1._q
            IF ( GRID%NGX /= GRID%NGX_rd .AND. (N1/=1 .AND. (N1/= GRID%NGX_rd))) FACTM=2._q
            IF ( GRID%NGZ /= GRID%NGZ_rd .AND. (N3/=1 .AND. (N3/= GRID%NGZ_rd))) FACTM=2._q

            IF ((GRID%LPCTX(N1)==0).AND.(GRID%LPCTY(N2)==0).AND.&
                 (GRID%LPCTZ(N3)==0) .AND. (NK==NQ)) THEN
               POTFAK(NI,0)=FSG(0)
               POTFAK(NI,1)=FSG(1)*FACTM
               POTFAK(NI,2)=FSG(2)*FACTM
               POTFAK(NI,3)=FSG(3)*FACTM
               POTFAK(NI,4)=FSG(4)*FACTM
               POTFAK(NI,5)=FSG(5)*FACTM
               POTFAK(NI,6)=FSG(6)*FACTM
            ELSE
               IF (HFSCREEN==0) THEN
                  ! derivative  -1/2 d (1/OMEGA/G^2) / d t_ij = G_i G_j/G^4 /OMEGA  - 1/OMEGA/G^2/2 delta_ij
                  POTFAK(NI,0)=SCALE/(GSQU)
                  EXP_GSQU_DER=SCALE/GSQU2
               ELSE IF (L_THOMAS_FERMI) THEN
                  POTFAK(NI,0)=SCALE/(GSQU+HFSCREEN*HFSCREEN/(TPI*TPI))
               ELSE IF (LRHFCALC) THEN
                  POTFAK(NI,0)=SCALE/GSQU*EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN)))
                  EXP_GSQU_DER=SCALE*EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN))) &
                              *(1/GSQU2+TPI*TPI/(4*HFSCREEN*HFSCREEN)/GSQU)
               ELSE IF (L_MODEL_HF) THEN
                  ! AEXX * case(HFSCREEN==0) + (1-AEXX)* default
                  POTFAK(NI,0)=SCALE/GSQU*(1-(1-AEXX)*EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN))))
               ELSE
                  POTFAK(NI,0)=(SCALE/GSQU )*(1-EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN))))
                  ! negative derivative of previous equation
                  EXP_GSQU_DER=(SCALE/GSQU2)*(1-EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN))))- &
                               (SCALE/GSQU )*(EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN))))* & 
                               (TPI*TPI/(4*HFSCREEN*HFSCREEN))
               ENDIF
               POTFAK(NI,1)=(GX*GX*EXP_GSQU_DER-POTFAK(NI,0)/2)*FACTM
               POTFAK(NI,2)=(GX*GY*EXP_GSQU_DER)*FACTM
               POTFAK(NI,3)=(GY*GY*EXP_GSQU_DER-POTFAK(NI,0)/2)*FACTM
               POTFAK(NI,4)=(GX*GZ*EXP_GSQU_DER)*FACTM
               POTFAK(NI,5)=(GY*GZ*EXP_GSQU_DER)*FACTM
               POTFAK(NI,6)=(GZ*GZ*EXP_GSQU_DER-POTFAK(NI,0)/2)*FACTM
            ENDIF
         ENDDO row
      ENDDO col

      IF (L_MODEL_HF) THEN
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(GRID,POTFAK) __IF_ASYNC__
         DO I=0,6
            DO NC=1,GRID%RC%NP
               POTFAK(NC,I)=POTFAK(NC,I)*(1.0_q/GRID%NPLWV)
            ENDDO
         ENDDO
      ELSE
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(GRID,POTFAK) __IF_ASYNC__
         DO I=0,6
            DO NC=1,GRID%RC%NP
               POTFAK(NC,I)=POTFAK(NC,I)*(AEXX/GRID%NPLWV)
            ENDDO
         ENDDO
      ENDIF
      
      DO I=0,6
         CALL SETUNB_REAL(POTFAK(1,I),GRID)
      ENDDO

      PROFILING_STOP('fock_set_gfac_der')

    END SUBROUTINE SET_GFAC_DER


!********************** SUBROUTINE SET_FSG_DER  ************************
!
!> integratable singularity of exchange term:
!> calculate the derivative of the average electrostatic potential prefactor
!> for the orbitals with k=k' and n=n' (i.e. those that are affected
!>  by divergence of the HF term in periodic boundary conditions)
!>
!>  Thesis gK Appendix B might be helpfull
!>  for the Coloum kernel this is:
!>  d / d t_ij  C/ Omega sum_G e(-alpha G^2) / G^2
!>       =  delta_ij C/ Omega sum_G e(-alpha G^2) / G^2 +
!>        + C/ Omega sum_G e(-alpha G^2) ( alpha / G^2 + 1/ G^4) G_i G_j
!
!**********************************************************************

    SUBROUTINE SET_FSG_DER(GRID, LATT_CUR, FSG)
      USE prec
      USE constant
      USE mgrid
      USE nonl_high
      USE lattice
      USE full_kpoints
      IMPLICIT NONE
      
      TYPE (grid_3d) GRID
      TYPE (latt) LATT_CUR
      INTEGER NK,NQ
      REAL(q) :: FSG(0:6)
      REAL(q) :: HFALPHAP
    ! local
      INTEGER    NC,N1,N2,N3, NN1, NN2, NN3, LPX, LPY, LPZ
      REAL(q) :: DKX,DKY,DKZ,GX,GY,GZ,GSQU,GSQU2,SCALE,EXP_GSQU_DER,FACTM
      INTEGER :: IDIR, JDIR, I
      REAL(q) :: FSG1, FSG2
      TYPE (latt) LATT_EWALD, LATT_FIN
      REAL(q) :: DIS=fd_displacement

      PROFILING_START('fock_set_fsg_der')

      CALL CHECK_FULL_KPOINTS
!
! this version is for molecules
! stress calculations (which require the derivative of FSG with respect to the volume)
! are not applicable
!
      IF (HFRCUT/=0) THEN
         FSG=0
         FSG(0)=HFRCUT*HFRCUT*EDEPS/LATT_CUR%OMEGA/2
!
! simple version, that can be applied for
! a simple Monkhorst Pack or Gamma centered grid
      ELSE IF (KPOINTS_FULL%NKPTS_NON_ZERO==KPOINTS_FULL%NKPX*KPOINTS_FULL%NKPY*KPOINTS_FULL%NKPZ .AND. HFSCREEN==0 &
         .AND. .NOT. ODDONLY .AND. .NOT. EVENONLY) THEN
         LATT_EWALD=LATT_CUR

         LATT_EWALD%A(:,1)=LATT_EWALD%A(:,1)*KPOINTS_FULL%NKPX/NKREDX
         LATT_EWALD%A(:,2)=LATT_EWALD%A(:,2)*KPOINTS_FULL%NKPY/NKREDY
         LATT_EWALD%A(:,3)=LATT_EWALD%A(:,3)*KPOINTS_FULL%NKPZ/NKREDZ

         CALL LATTIC(LATT_EWALD)
         CALL EWALD_MONO(FSG(0),1.0_q,LATT_EWALD)
         FSG(0)=FSG(0)*2

         I=0
         DO IDIR=1,3
            DO JDIR=1,IDIR
               I=I+1
               LATT_FIN=LATT_EWALD
               LATT_FIN%A(IDIR,:)=LATT_EWALD%A(IDIR,:)+DIS*LATT_EWALD%A(JDIR,:)
               CALL LATTIC(LATT_FIN)
               CALL EWALD_MONO(FSG1,1.0_q,LATT_FIN)
               LATT_FIN%A(IDIR,:)=LATT_EWALD%A(IDIR,:)-DIS*LATT_EWALD%A(JDIR,:)
               CALL LATTIC(LATT_FIN)
               CALL EWALD_MONO(FSG2,1.0_q,LATT_FIN)
               FSG(I)=(FSG1-FSG2)/2/DIS
            ENDDO
         ENDDO
      ELSE

         FSG=0
         NK=1
         DO NQ=1,KPOINTS_FULL%NKPTS
            SCALE=KPOINTS_FULL%WTKPT(NQ)*NKREDX*NKREDY*NKREDZ *EDEPS/LATT_CUR%OMEGA/TPI**2
            IF (ODDONLY .OR. EVENONLY) SCALE=SCALE*2
            IF (SKIP_THIS_KPOINT_IN_FOCK(KPOINTS_FULL%VKPT(:,NQ)-KPOINTS_FULL%VKPT(:,NK))) CYCLE

#ifndef gammareal
            DKX=(KPOINTS_FULL%VKPT(1,NK)-KPOINTS_FULL%VKPT(1,NQ))*LATT_CUR%B(1,1)+ &
                 (KPOINTS_FULL%VKPT(2,NK)-KPOINTS_FULL%VKPT(2,NQ))*LATT_CUR%B(1,2)+ &
                 (KPOINTS_FULL%VKPT(3,NK)-KPOINTS_FULL%VKPT(3,NQ))*LATT_CUR%B(1,3)
            DKY=(KPOINTS_FULL%VKPT(1,NK)-KPOINTS_FULL%VKPT(1,NQ))*LATT_CUR%B(2,1)+ &
                 (KPOINTS_FULL%VKPT(2,NK)-KPOINTS_FULL%VKPT(2,NQ))*LATT_CUR%B(2,2)+ &
                 (KPOINTS_FULL%VKPT(3,NK)-KPOINTS_FULL%VKPT(3,NQ))*LATT_CUR%B(2,3)
            DKZ=(KPOINTS_FULL%VKPT(1,NK)-KPOINTS_FULL%VKPT(1,NQ))*LATT_CUR%B(3,1)+ &
                 (KPOINTS_FULL%VKPT(2,NK)-KPOINTS_FULL%VKPT(2,NQ))*LATT_CUR%B(3,2)+ &
                 (KPOINTS_FULL%VKPT(3,NK)-KPOINTS_FULL%VKPT(3,NQ))*LATT_CUR%B(3,3)
#else
            DKX=0._q; DKY=0._q; DKZ=0._q
#endif
            DO NN1=0,GRID%NGX,GRID%NGX
            DO NN2=0,GRID%NGY,GRID%NGY
            DO NN3=0,GRID%NGZ,GRID%NGZ

            col: DO NC=1,GRID%RC%NCOL
               N2=GRID%RC%I2(NC)
               N3=GRID%RC%I3(NC)
               
               row: DO N1=1,GRID%RC%NROW
                  LPX=MOD(GRID%LPCTX(N1)+NN1+GRID%NGX-1,2*GRID%NGX)-GRID%NGX+1
                  LPY=MOD(GRID%LPCTY(N2)+NN2+GRID%NGY-1,2*GRID%NGY)-GRID%NGY+1
                  LPZ=MOD(GRID%LPCTZ(N3)+NN3+GRID%NGZ-1,2*GRID%NGZ)-GRID%NGZ+1

                  GX=(LPX*LATT_CUR%B(1,1)+LPY*LATT_CUR%B(1,2)+LPZ*LATT_CUR%B(1,3))+DKX
                  GY=(LPX*LATT_CUR%B(2,1)+LPY*LATT_CUR%B(2,2)+LPZ*LATT_CUR%B(2,3))+DKY
                  GZ=(LPX*LATT_CUR%B(3,1)+LPY*LATT_CUR%B(3,2)+LPZ*LATT_CUR%B(3,3))+DKZ
                  GSQU=GX**2+GY**2+GZ**2
                  GSQU2=GSQU*GSQU
                  
                  FACTM=1._q
                  ! IF (GRID%LREAL .AND. (N1/=1 .AND. (N1/=GRID%RC%NROW))) FACTM=2._q
                  IF ( GRID%NGX /= GRID%NGX_rd .AND. (N1/=1 .AND. (N1/= GRID%NGX_rd))) FACTM=2._q
                  IF ( GRID%NGZ /= GRID%NGZ_rd .AND. (N3/=1 .AND. (N3/= GRID%NGZ_rd))) FACTM=2._q
                  
                  IF ((GRID%LPCTX(N1)==0).AND.(GRID%LPCTY(N2)==0).AND.&
                       (GRID%LPCTZ(N3)==0) .AND. (NK==NQ)) THEN
                  ELSE
                     IF (HFSCREEN==0) THEN
                        ! Coloumb kernel
                        FSG(0)=FSG(0)+FACTM*SCALE*EXP(-GSQU*(HFALPHA*TPI*TPI))/GSQU
                        ! negative derivative of Coloumb kernel with  respect to G^2 (GSQU)
                        EXP_GSQU_DER =FACTM*SCALE*EXP(-GSQU*(HFALPHA*TPI*TPI))*(1/GSQU2+HFALPHA*TPI*TPI/GSQU)
                     ELSE IF (L_THOMAS_FERMI) THEN
                        FSG(0)=FSG(0)+FACTM*SCALE*EXP(-GSQU*(HFALPHA*TPI*TPI))/(GSQU+HFSCREEN*HFSCREEN*(1.0_q/(TPI*TPI)))
                        CALL vtutor%error("Stress not implemented")
                     ELSE IF (LRHFCALC) THEN
                        ! long range Coulomb kernel
                        FSG(0)=FSG(0)+FACTM*SCALE*EXP(GSQU*((-HFALPHA-1/(4*HFSCREEN*HFSCREEN))*TPI*TPI))/GSQU
!                       WRITE(*,*) 'Stress not implemented'
!                       STOP
                        ! negative derivative of kernel with  respect to G^2 (GSQU)
                        EXP_GSQU_DER =FACTM*SCALE*EXP(GSQU*((-HFALPHA-1/(4*HFSCREEN*HFSCREEN))*TPI*TPI)) &
                                    *(1/GSQU2+(HFALPHA+1/(4*HFSCREEN*HFSCREEN))*TPI*TPI/GSQU)

                     ELSE IF (L_MODEL_HF) THEN
                        FSG=FSG+FACTM*SCALE/GSQU*EXP(-GSQU*(HFALPHA*TPI*TPI)) &
                             *(1-(1-AEXX)*EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN))))
                        CALL vtutor%error("Stress not implemented")
                     ELSE
                        ! screened Coulomb kernel (short range)
                        FSG(0)=FSG(0)+FACTM*SCALE*EXP(-GSQU*(HFALPHA*TPI*TPI))/GSQU &
                             *(1-EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN))))
                        !
                        ! negative derivative of kernel with  respect to G^2 (GSQU)
                        EXP_GSQU_DER=FACTM*SCALE*EXP(-GSQU*(HFALPHA*TPI*TPI))*(1/GSQU2+(HFALPHA*TPI*TPI)/GSQU) &
                             *(1-EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN))))- &
                             FACTM*SCALE*EXP(-GSQU*(HFALPHA*TPI*TPI))/GSQU* &
                             EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN)))*(TPI*TPI/(4*HFSCREEN*HFSCREEN))
                     ENDIF
                     ! derivative of d G^2 / d t_ij = - 2 G_i G_j
                     FSG(1)=FSG(1)+GX*GX*EXP_GSQU_DER
                     FSG(2)=FSG(2)+GX*GY*EXP_GSQU_DER
                     FSG(3)=FSG(3)+GY*GY*EXP_GSQU_DER
                     FSG(4)=FSG(4)+GX*GZ*EXP_GSQU_DER
                     FSG(5)=FSG(5)+GY*GZ*EXP_GSQU_DER
                     FSG(6)=FSG(6)+GZ*GZ*EXP_GSQU_DER
                  ENDIF
               ENDDO row
            ENDDO col
            ENDDO
            ENDDO
            ENDDO
         ENDDO
         
         CALLMPI (M_sum_d(GRID%COMM,FSG,7))

         ! add contribution from one over volume dependency of FSG
         FSG(1)=FSG(1)-FSG(0)/2
         FSG(3)=FSG(3)-FSG(0)/2
         FSG(6)=FSG(6)-FSG(0)/2

      
         IF (HFSCREEN==0) THEN
            ! selfenergy of an isolated erfc(r) charge distribution
            FSG(0)=-FSG(0)+EDEPS/(TPI**3._q)*2._q*HFALPHA*(PI/HFALPHA)**(3._q/2._q)
            
         ELSE IF (L_THOMAS_FERMI) THEN
            ! Thomas Fermi, no analytical formula but small hence zero it
            FSG=0
         ELSE IF (LRHFCALC) THEN
            HFALPHAP=HFALPHA+1/HFSCREEN/HFSCREEN/4
            !
            ! selfenergy of an erfc(r) charge distribution 
            ! for the long range Coulomb kernel
            !
            FSG(0)=-FSG(0)+EDEPS/(TPI**3._q)*2._q*HFALPHAP*(PI/HFALPHAP)**(3._q/2._q)
         ELSE IF (L_MODEL_HF) THEN
            !
            ! AEXX * case(HFSCREEN==0) + (1-AEXX)* default
            !
            HFALPHAP=HFALPHA+1/HFSCREEN/HFSCREEN/4
            FSG(0)=-FSG(0)+EDEPS/(TPI**3)*2._q*HFALPHA *(PI/HFALPHA )**(3._q/2._q)- &
                 (1-AEXX)*EDEPS/(TPI**3)*2._q*HFALPHAP*(PI/HFALPHAP)**(3._q/2._q)
         ELSE 
            HFALPHAP=HFALPHA+1/HFSCREEN/HFSCREEN/4
            !
            ! selfenergy of an erfc(r) charge distribution in
            ! a screened short range Coloumb potential 
            !
            FSG(0)=-FSG(0)+EDEPS/(TPI**3._q)*2._q*HFALPHA *(PI/HFALPHA )**(3._q/2._q)- &
                 EDEPS/(TPI**3._q)*2._q*HFALPHAP*(PI/HFALPHAP)**(3._q/2._q)
         ENDIF
         
         FSG(1:6)=-FSG(1:6)
      ENDIF

      PROFILING_STOP('fock_set_fsg_der')

    END SUBROUTINE 

!***********************************************************************
!
!> Determine whether a specific q point (usually k1-k2) should
!> be used in the calculation of Bloch integrals such
!> as exact exchange
! 
!***********************************************************************

    FUNCTION SKIP_THIS_KPOINT_IN_FOCK(VKPT)
      USE full_kpoints
      REAL(q) VKPT(3)
      LOGICAL SKIP_THIS_KPOINT_IN_FOCK

      SKIP_THIS_KPOINT_IN_FOCK=.FALSE.
      IF (ODDONLY .AND. ABS(MODULO(NINT(VKPT(1)*KPOINTS_FULL%NKPX+ & 
                                        VKPT(2)*KPOINTS_FULL%NKPY+ & 
                                        VKPT(3)*KPOINTS_FULL%NKPZ),2))<=1E-6)  SKIP_THIS_KPOINT_IN_FOCK=.TRUE.
      IF (EVENONLY.AND. ABS(MODULO(NINT(VKPT(1)*KPOINTS_FULL%NKPX+ & 
                                        VKPT(2)*KPOINTS_FULL%NKPY+ & 
                                        VKPT(3)*KPOINTS_FULL%NKPZ+1),2))<=1E-6)  SKIP_THIS_KPOINT_IN_FOCK=.TRUE.
      IF (SHIFTRED) THEN
         IF (NKREDX==2 .AND. MOD(FLOOR(VKPT(1)*KPOINTS_FULL%NKPX+100.5_q),NKREDX)==0) SKIP_THIS_KPOINT_IN_FOCK=.TRUE.
         IF (NKREDY==2 .AND. MOD(FLOOR(VKPT(2)*KPOINTS_FULL%NKPY+100.5_q),NKREDY)==0) SKIP_THIS_KPOINT_IN_FOCK=.TRUE.
         IF (NKREDY==2 .AND. MOD(FLOOR(VKPT(3)*KPOINTS_FULL%NKPZ+100.5_q),NKREDZ)==0) SKIP_THIS_KPOINT_IN_FOCK=.TRUE.
      ELSE IF (NKREDX>=2 .OR. NKREDY>=2 .OR. NKREDZ>=2) THEN
         IF (MOD(FLOOR((VKPT(1)+32)*KPOINTS_FULL%NKPX+.5_q),NKREDX)/=0) SKIP_THIS_KPOINT_IN_FOCK=.TRUE.
         IF (MOD(FLOOR((VKPT(2)+32)*KPOINTS_FULL%NKPY+.5_q),NKREDY)/=0) SKIP_THIS_KPOINT_IN_FOCK=.TRUE.
         IF (MOD(FLOOR((VKPT(3)+32)*KPOINTS_FULL%NKPZ+.5_q),NKREDZ)/=0) SKIP_THIS_KPOINT_IN_FOCK=.TRUE.
      ENDIF
    END FUNCTION SKIP_THIS_KPOINT_IN_FOCK

!***********************************************************************
!
!> allocate the handle for HF type calculations
!>
!> @details @ref openmp : Under OpenMP the arrays
!>    fock_handle::omp_lck1 and fock_handle::omp_lck2
!> are allocated and initialized as OpenMP locks.
!> These are used to prevent data races in the access of 
!>    fock_handle::wq (\%CR) and fock_handle::ckappa
!> in wave_cacher::store_gw_acc.
!
!***********************************************************************

    SUBROUTINE ALLOCATE_FOCK_HANDLE( FH , LMDIM , NBANDSGW)
      USE ini
      TYPE (fock_handle) :: FH
      INTEGER LMDIM, NBANDSGW
!$    INTEGER i

!     ALLOCATE(FH)

      FH%LMDIM=LMDIM
      ALLOCATE( &
           FH%CXI(GRID_FOCK%MPLWV*WDES_FOCK%NRSPINORS, NBANDSGW, 2), &
           FH%CKAPPA(WDES_FOCK%NPROD, NBANDSGW, 2),  &
           FH%CDIJ(LMDIM,LMDIM,WDES_FOCK%NIONS,WDES_FOCK%NRSPINORS), &
           FH%CDLM(AUG_DES%NPROD*WDES_FOCK%NRSPINORS))

      CALL REGISTER_ALLOCATE(16._q*SIZE( FH%CXI), "fockhandle")
#ifdef gammareal
      CALL REGISTER_ALLOCATE(8._q*SIZE( FH%CKAPPA), "fockhandle")
      CALL REGISTER_ALLOCATE(8._q*SIZE( FH%CDIJ), "fockhandle")
      CALL REGISTER_ALLOCATE(8._q*SIZE( FH%CDLM), "fockhandle")
#else
      CALL REGISTER_ALLOCATE(16._q*SIZE( FH%CKAPPA), "fockhandle")
      CALL REGISTER_ALLOCATE(16._q*SIZE( FH%CDIJ), "fockhandle")
      CALL REGISTER_ALLOCATE(16._q*SIZE( FH%CDLM), "fockhandle")
#endif
      
      CALL SETWDES(WDES_FOCK, FH%WDESK, 0)

      CALL NEWWAV(FH%WQ, FH%WDESK, .TRUE.)
      FH%W1%CPROJ=>FH%CDLM(:)

      CALL CLEAR_FOCK_HANDLE( FH)

!$    ALLOCATE(FH%OMP_LCK1(NBANDSGW),FH%OMP_LCK2(NBANDSGW))
!$    DO i=1,NBANDSGW
!$       CALL OMP_INIT_LOCK(FH%OMP_LCK1(i))
!$       CALL OMP_INIT_LOCK(FH%OMP_LCK2(i))
!$    ENDDO

    END SUBROUTINE ALLOCATE_FOCK_HANDLE


    SUBROUTINE ALLOCATE_FOCK_HANDLE_NOACC( FH , LMDIM)
      USE ini
      TYPE (fock_handle) :: FH
      INTEGER LMDIM, NBANDSGW

      FH%LMDIM=LMDIM
      ALLOCATE( &
           FH%CDIJ(LMDIM,LMDIM,WDES_FOCK%NIONS,WDES_FOCK%NRSPINORS), &
           FH%CDLM(AUG_DES%NPROD*WDES_FOCK%NRSPINORS))

#ifdef gammareal
      CALL REGISTER_ALLOCATE(8._q*SIZE( FH%CDIJ), "fockhandle")
      CALL REGISTER_ALLOCATE(8._q*SIZE( FH%CDLM), "fockhandle")
#else
      CALL REGISTER_ALLOCATE(16._q*SIZE( FH%CDIJ), "fockhandle")
      CALL REGISTER_ALLOCATE(16._q*SIZE( FH%CDLM), "fockhandle")
#endif
      
      CALL SETWDES(WDES_FOCK, FH%WDESK, 0)

      CALL NEWWAV(FH%WQ, FH%WDESK, .TRUE.)
      FH%W1%CPROJ=>FH%CDLM(:)

      FH%CDIJ=0; FH%CDLM=0

    END SUBROUTINE ALLOCATE_FOCK_HANDLE_NOACC


!***********************************************************************
! 
!> @details @ref openmp : the OpenMP locks
!>    fock_handle::omp_lck1 and fock_handle::omp_lck2
!> are destroyed and the corresponding arrays deallocated.
!
!***********************************************************************

    SUBROUTINE DEALLOCATE_FOCK_HANDLE( FH)
      USE ini
      TYPE (fock_handle) :: FH
      TYPE (wavespin) W
      INTEGER LMDIM
!$    INTEGER i

      CALL DEREGISTER_ALLOCATE(16._q*SIZE( FH%CXI), "fockhandle")
#ifdef gammareal
      CALL DEREGISTER_ALLOCATE(8._q*SIZE( FH%CKAPPA), "fockhandle")
      CALL DEREGISTER_ALLOCATE(8._q*SIZE( FH%CDIJ), "fockhandle")
      CALL DEREGISTER_ALLOCATE(8._q*SIZE( FH%CDLM), "fockhandle")
#else
      CALL DEREGISTER_ALLOCATE(16._q*SIZE( FH%CKAPPA), "fockhandle")
      CALL DEREGISTER_ALLOCATE(16._q*SIZE( FH%CDIJ), "fockhandle")
      CALL DEREGISTER_ALLOCATE(16._q*SIZE( FH%CDLM), "fockhandle")
#endif

      DEALLOCATE(FH%CXI, FH%CKAPPA, FH%CDIJ, FH%CDLM)

      NULLIFY(FH%W1%CPROJ)
      CALL DELWAV(FH%WQ, .TRUE.)

!$    DO i=1,SIZE(FH%OMP_LCK1)
!$       CALL OMP_DESTROY_LOCK(FH%OMP_LCK1(i))
!$       CALL OMP_DESTROY_LOCK(FH%OMP_LCK2(i))
!$    ENDDO
!$    DEALLOCATE(FH%OMP_LCK1,FH%OMP_LCK2)

!     DEALLOCATE(FH)
      
    END SUBROUTINE DEALLOCATE_FOCK_HANDLE


    SUBROUTINE DEALLOCATE_FOCK_HANDLE_NOACC( FH)
      USE ini
      TYPE (fock_handle) :: FH
      TYPE (wavespin) W
      INTEGER LMDIM

#ifdef gammareal
      CALL DEREGISTER_ALLOCATE(8._q*SIZE( FH%CDIJ), "fockhandle")
      CALL DEREGISTER_ALLOCATE(8._q*SIZE( FH%CDLM), "fockhandle")
#else
      CALL DEREGISTER_ALLOCATE(16._q*SIZE( FH%CDIJ), "fockhandle")
      CALL DEREGISTER_ALLOCATE(16._q*SIZE( FH%CDLM), "fockhandle")
#endif

      DEALLOCATE(FH%CDIJ, FH%CDLM)

      NULLIFY(FH%W1%CPROJ)
      CALL DELWAV(FH%WQ, .TRUE.)

    END SUBROUTINE DEALLOCATE_FOCK_HANDLE_NOACC


    SUBROUTINE CLEAR_FOCK_HANDLE( FH )
      TYPE (fock_handle) :: FH

      FH%CXI   =0
      FH%CKAPPA=0
      FH%CDIJ  =0
      FH%CDLM  =0  
    END SUBROUTINE CLEAR_FOCK_HANDLE


!***********************************************************************
!
!> Determine the action of the Hamiltonian on all bands and
!> store the result in a set of wavefunctions WXI
! 
!***********************************************************************


  SUBROUTINE FOCK_ACC_ALL(GRID, LATT_CUR, NONLR_S, NONL_S, W, WXI, &
     &    LMDIM, LKINETIC, P, CQIJ, EXHF)
      USE prec
      USE wave_mpi
      USE wave
      USE lattice
      USE mpimy
      USE mgrid
      
      USE nonl_high
      USE pseudo

      IMPLICIT NONE

      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (potcar)      P(:)
      TYPE (wavespin), TARGET :: W
      TYPE (wavespin)    WXI
      INTEGER LMDIM
      LOGICAL LKINETIC
      INTEGER IU0
      OVERLAP            CQIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
      REAL(q) EXHF
!    local
      TYPE (wavedes), POINTER :: WDES
      TYPE (wavedes1)    WDES1          ! descriptor for one k-point
      TYPE (wavefun1)    W1             ! current wavefunction
      COMPLEX(q) CDCHF
      INTEGER ISP, NK, N, M, MM, NSTRIP, NSTRIP_ACT, NPOS, ISPINOR

      WDES=>W%WDES
      NSTRIP=((W%WDES%NSIM*2+W%WDES%NB_PAR-1)/W%WDES%NB_PAR)

      WXI%CPROJ=0
      WXI%CW=0

      CDCHF=0

      spin:  DO ISP=1,WDES%ISPIN
      kpoint: DO NK=1,WDES%NKPTS
#ifdef MPI
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
         WXI%CW(:,:,NK,ISP)=0   ! Zero accelerations
      ELSE
#endif
         IF (NONLR_S%LREAL) THEN
            CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES)
         ELSE
            CALL PHASE(WDES,NONL_S,NK)
         ENDIF

         CALL SETWDES(WDES,WDES1,NK); CALL SETWGRID_OLD(WDES1,GRID)
         DO NPOS=1,WDES%NBANDS,NSTRIP
            NSTRIP_ACT=MIN(WDES%NBANDS+1-NPOS,NSTRIP)
            IF (AEXX/=0) THEN
               CALL FOCK_ACC(GRID, LMDIM, LATT_CUR, W,  &
                 NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIP_ACT, &
                 WXI%CW(:,NPOS:,NK,ISP), P, CQIJ, CDCHF )
            ELSE
               WXI%CW(:,NPOS:NPOS+NSTRIP_ACT-1,NK,ISP)=0
            ENDIF
         ENDDO
         IF (LKINETIC) THEN
         DO NPOS=1,WDES%NBANDS
            DO ISPINOR=0,WDES1%NRSPINORS-1
               DO M=1,WDES1%NGVECTOR
                  MM=M+ISPINOR*WDES1%NGVECTOR
                  WXI%CW(MM,NPOS,NK,ISP)=WXI%CW(MM,NPOS,NK,ISP)+W%CW(MM,NPOS,NK,ISP)*WDES1%DATAKE(M,ISPINOR+1)
               ENDDO
            ENDDO
         ENDDO
         ENDIF
#ifdef MPI
      ENDIF
#endif
      END DO kpoint
      END DO spin

      CALLMPI( M_sum_z(WDES1%COMM_KIN,CDCHF,1))
      CALLMPI( M_sum_z(WDES1%COMM_KINTER,CDCHF,1))
      EXHF=CDCHF
      RETURN
    END SUBROUTINE FOCK_ACC_ALL


    SUBROUTINE FOCK_ALL( LATT_CUR, NONLR_S, NONL_S, W, &
     &    LMDIM, P, CQIJ, EXHF)
      USE prec
      USE wave_mpi
      USE wave
      USE lattice
      USE mpimy
      USE mgrid
      
      USE nonl_high
      USE pseudo

      IMPLICIT NONE

      TYPE (latt)        LATT_CUR
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (potcar)      P(:)
      TYPE (wavespin), TARGET :: W
      INTEGER LMDIM
      LOGICAL LKINETIC
      INTEGER IU0
      OVERLAP            CQIJ (LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
      REAL(q) EXHF
!    local
      TYPE (wavedes), POINTER :: WDES
      TYPE (wavedes1) WDES1          ! descriptor for one k-point
      COMPLEX(q) CDCHF
      INTEGER ISP, NK, NSTRIP, NSTRIP_ACT, NPOS, ISPINOR
      TYPE (wavefuna) WXI 
 
      WDES=>W%WDES
      NSTRIP=((W%WDES%NSIM*2+W%WDES%NB_PAR-1)/W%WDES%NB_PAR)
  

      CALL SETWDES(WDES,WDES1,0)
      CALL NEWWAVA(WXI, WDES1, NSTRIP)
 
      CDCHF=0
      
      spin:  DO ISP=1,WDES%ISPIN
      kpoint: DO NK=1,WDES%NKPTS
#ifdef MPI
         IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).EQ.WDES%COMM_KINTER%NODE_ME-1) THEN
#endif

         WXI%CW(:,:)=0      ! Zero accelerations
         IF (NONLR_S%LREAL) THEN
            CALL PHASER(W%WDES%GRID,LATT_CUR,NONLR_S,NK,WDES)
         ELSE
            CALL PHASE(WDES,NONL_S,NK)
         ENDIF

         CALL SETWDES(WDES,WDES1,NK); CALL SETWGRID_OLD(WDES1,W%WDES%GRID)
         DO NPOS=1,WDES%NBANDS,NSTRIP
            NSTRIP_ACT=MIN(WDES%NBANDS+1-NPOS,NSTRIP)
            IF (AEXX/=0) THEN
               CALL FOCK_ACC(W%WDES%GRID, LMDIM, LATT_CUR, W,  &
                 NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIP_ACT, &
                 WXI%CW, P, CQIJ, CDCHF )
            ENDIF
         ENDDO
#ifdef MPI
      ENDIF
#endif
      END DO kpoint
      END DO spin

      CALLMPI( M_sum_z(WDES1%COMM_KIN,CDCHF,1))
      CALLMPI( M_sum_z(WDES1%COMM_KINTER,CDCHF,1))
      EXHF=CDCHF
      RETURN
    END SUBROUTINE FOCK_ALL


!************************ SUBROUTINE FOCK_CHARGE ***********************
! 
!> calculate the total charge density (i.e. plane wave +
!> fast augmentation) from two orbitals
!>   w1 * w2^*
!
!***********************************************************************

  SUBROUTINE FOCK_CHARGE( W1, W2, GCHG, CRHOLM)
    TYPE (wavefun1) :: W1, W2  !< two wavefunctions
    GDEF ::  GCHG(:)           !< accumulated charged
    GDEF ::  CRHOLM(:)         !< temporary

    ! get the plane wave contribution to the charge
    CALL PW_CHARGE_TRACE(W1%WDES1, GCHG(1), W1%CR(1), W2%CR(1))

    ! add augmentation part to charge (if required)
    IF (W1%WDES1%LOVERL) THEN
       CALL DEPSUM_TWO_BANDS_RHOLM_TRACE(W1%CPROJ(:),W2%CPROJ(:), W1%WDES1, AUG_DES, &
            TRANS_MATRIX_FOCK, CRHOLM, 1._q, W1%WDES1%LOVERL)
      
       AUG_DES%RINPL=1._q ! multiplicator used by RACC0
       CALL RACC0_HF(FAST_AUG_FOCK, AUG_DES, CRHOLM(1), GCHG(1))
    ENDIF
  END SUBROUTINE FOCK_CHARGE


!************************ SUBROUTINE FOCK_CHARGE_MU ********************
!
!> blocked version working on a set of input orbitals W1
!
!***********************************************************************

  SUBROUTINE FOCK_CHARGE_MU( W1, W2, GCHG, CRHOLM)
#ifdef _OPENACC
    USE mopenacc
#endif
    USE openmp, ONLY : omp_nthreads
    TYPE (wavefun1) :: W1(:), W2  !< two wavefunctions
    GDEF ::  GCHG(:,:)            !< accumulated charged
    GDEF ::  CRHOLM(:,:)          !< temporary
    INTEGER N

    PROFILING_START('fock_charge_mu')

#ifndef _OPENACC
!$OMP PARALLEL DO PRIVATE(N)
    DO N=1,SIZE(W1)
    ! get the plane wave contribution to the charge
    CALL PW_CHARGE_TRACE(W1(N)%WDES1, GCHG(1,N), W1(N)%CR(1), W2%CR(1))

    ! add augmentation part to charge (if required)
    IF (W1(N)%WDES1%LOVERL) THEN
       CALL DEPSUM_TWO_BANDS_RHOLM_TRACE(W1(N)%CPROJ(:),W2%CPROJ(:), W1(N)%WDES1, AUG_DES, &
            TRANS_MATRIX_FOCK, CRHOLM(:,N), 1._q, W1(N)%WDES1%LOVERL)
    ENDIF
    ENDDO
!$OMP END PARALLEL DO
#else
     N=SIZE(W1)
     IF (ACC_EXEC_ON) THEN
        CALL PW_CHARGE_TRACE_MU_ACC(W1(1)%WDES1, GCHG(1,1), SIZE(GCHG,1), W1, W2%CR(1), N)
     ELSE
        CALL PW_CHARGE_TRACE_MU    (W1(1)%WDES1, GCHG(1,1), SIZE(GCHG,1), W1, W2%CR(1), N)
     ENDIF

    ! add augmentation part to charge (if required)
     CALL DEPSUM_TWO_BANDS_RHOLM_TRACE_MU(W1, W2%CPROJ(:), W1(1)%WDES1, AUG_DES, &
          TRANS_MATRIX_FOCK, CRHOLM(:,:), 1._q, W1(N)%WDES1%LOVERL, N)
#endif

    IF (W1(1)%WDES1%LOVERL) THEN
       AUG_DES%RINPL=1._q       ! multiplicator used by RACC0
!NOTE(sm): we do not replicate this change on the GPU because only RACC0MU_HF_ACC
!          would be using it. This is a specialized ACC routine, so we just remove
!          the multiplication there.
       IF (SIZE(CRHOLM,1) /= AUG_DES%NPROD*W1(1)%WDES1%NRSPINORS) THEN
          CALL vtutor%bug("internal error in VASP: FOCK_CHARGE_MU size mismatch " // str(SIZE(CRHOLM,&
             1)) // " " // str(AUG_DES%NPROD), __FILE__, __LINE__)
       ENDIF

       N=SIZE(W1)
#ifdef _OPENACC
       IF (ACC_EXEC_ON) THEN
          CALL RACC0MU_HF_ACC(FAST_AUG_FOCK, AUG_DES, CRHOLM(1,1), SIZE(CRHOLM,1), GCHG(1,1), SIZE(GCHG,1), N)
          PROFILING_STOP('fock_charge_mu')
          RETURN
       ENDIF
#endif
       CALL RACC0MU_HF(FAST_AUG_FOCK, AUG_DES, CRHOLM(1,1), SIZE(CRHOLM,1), GCHG(1,1), SIZE(GCHG,1), N)
    ENDIF

    PROFILING_STOP('fock_charge_mu')

  END SUBROUTINE FOCK_CHARGE_MU

  !> only augmentation part
  SUBROUTINE FOCK_CHARGE_AUG( W1, W2, GCHG, CRHOLM)
    TYPE (wavefun1) :: W1, W2  !< two wavefunctions
    GDEF ::  GCHG(:)           !< accumulated charged
    GDEF ::  CRHOLM(:)         !< temporary

    ! add augmentation part to charge (if required)
    IF (W1%WDES1%LOVERL) THEN
       CALL DEPSUM_TWO_BANDS_RHOLM_TRACE(W1%CPROJ(:),W2%CPROJ(:), W1%WDES1, AUG_DES, &
            TRANS_MATRIX_FOCK, CRHOLM, 1._q, W1%WDES1%LOVERL)
       
       AUG_DES%RINPL=1._q ! multiplicator used by RACC0
       CALL RACC0_HF(FAST_AUG_FOCK, AUG_DES, CRHOLM(1), GCHG(1))
    ENDIF
  END SUBROUTINE FOCK_CHARGE_AUG

!> conjugated version
!> calculate the total charge density (i.e. plane wave +
!> fast augmentation) from two wavefunctions without conjugation of the second wavefunction
!>   w1 * w2

  SUBROUTINE FOCK_CHARGE_NO_CONJG( W1, W2, GCHG, CRHOLM)
    TYPE (wavefun1) :: W1, W2  !< two wavefunctions
    GDEF ::  GCHG(:)           !< accumulated charged
    GDEF ::  CRHOLM(:)         !< temporary

    ! get the plane wave contribution to the charge
    CALL PW_CHARGE_TRACE_NO_CONJG(W1%WDES1, GCHG(1), W1%CR(1), W2%CR(1))

    ! add augmentation part to charge (if required)
    IF (W1%WDES1%LOVERL) THEN
       CALL DEPSUM_TWO_BANDS_RHOLM_NO_CONJG(W1%CPROJ(:),W2%CPROJ(:), W1%WDES1, AUG_DES, &
            TRANS_MATRIX_FOCK, CRHOLM, 1._q, W1%WDES1%LOVERL)
       
       AUG_DES%RINPL=1._q ! multiplicator used by RACC0
       CALL RACC0_HF(FAST_AUG_FOCK, AUG_DES, CRHOLM(1), GCHG(1))
    ENDIF
  END SUBROUTINE FOCK_CHARGE_NO_CONJG


!************************ SUBROUTINE FOCK_ACC  *************************
!
!> This is the main subroutine for calculating the Fock exchange
!> operator accelerations  F | psi_n> for a stripe of bands from
!> NPOS to NPOS+NSTRIP at the k-points NK and for spin ISP
!> the wavefunction onto which the action should be evaluated is
!> either taken from the original array W or from WIN_ORIG
!> is supplied
!> the convergence corrections are only correct if WIN_ORIG
! 
!**********************************************************************

!
!> currently two version are available: the first one blocks the global
!> gather opereration and collects a maximum of NBLOCK_FOCK states at
!> once. The second version does not block the global gather and collects
!> NSTRIPN*NB_PAR states at once. The latter can use a lot of memory,
!> especially when NB_PAR is large.
!

#ifndef non_blocked
!
!> blocked version
!
!> @details @ref openmp :
!> the loop over the bands on which the action of the Fock potential
!> is calculated (label: nband) is distributed over the available
!> threads in case NGLB is a multiple of openmp::omp_nthreads.
!
  SUBROUTINE FOCK_ACC(GRID_, LMDIM, LATT_CUR, W, NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIPN, &
       CH, P, CQIJ, CDCHF, WIN_ORIG, LSYMGRAD, EXHF_ACFDT, LCOMMENSURATE)
    USE ini
    USE sym_prec
    USE nonl_high
    USE wave_high
    USE lattice
    USE constant
    USE full_kpoints
    USE pseudo
    USE openmp, ONLY : omp_nthreads
    IMPLICIT NONE

! passed variables
    INTEGER LMDIM
    TYPE (grid_3d) GRID_                          ! not used (scheduled for removal)
    TYPE (latt) LATT_CUR
    TYPE (wavespin) W
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct)  NONL_S
    TYPE (potcar)       P(:)
    INTEGER NK,ISP,NPOS,NSTRIPN
    COMPLEX(q)       :: CH(:,:)                   !< accelerations in rec. space
    OVERLAP             CQIJ(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    COMPLEX(q) CDCHF                              !< double counting correction
    TYPE (wavefun1), OPTIONAL :: WIN_ORIG(:)      !< original wavefunctions distributed over bands
    LOGICAL,OPTIONAL :: LSYMGRAD
    REAL(q),OPTIONAL :: EXHF_ACFDT                !< difference between HF exchange, and ACFDT exchange
                                                  !< compare equation (12) in Harl et al. PRB 81 115126 (2010)
    LOGICAL, OPTIONAL:: LCOMMENSURATE

! local variables
    TYPE (wavedes1), TARGET :: WDESK, WDESQ, WDESQ_IRZ
    TYPE (wavefun1), TARGET :: WQ
    TYPE (wavefun1),ALLOCATABLE :: W1(:)
    TYPE (wavefun1),ALLOCATABLE :: WIN(:),WXI(:)
    REAL(q) :: WEIGHT
    REAL(q) :: FSG                                ! singularity correction setting V(G=0) to a finite value
    REAL(q) :: FSG_AEXX                           ! alternative treatment, by adding original orbital
    INTEGER ISPINOR
    INTEGER N, NS, NT, MQ, NP, NGLB, MM, ISP_IRZ
    INTEGER NQ, NQ_USED
    LOGICAL LSHIFT
    GDEF,      ALLOCATABLE :: GWORK(:,:)               ! fock pot in real sp
    GDEF,      ALLOCATABLE :: CRHOLM(:,:)              ! augmentation occupancy matrix
    GDEF,      ALLOCATABLE :: CDIJ(:,:,:,:)            ! D_lml'm'
    COMPLEX(q),ALLOCATABLE, TARGET :: CXI(:,:)         ! acc. in real space
    GDEF,      ALLOCATABLE, TARGET :: CDLM(:,:)        ! D_LM
    GDEF,      ALLOCATABLE, TARGET :: CKAPPA(:,:)      ! stores NL accelerations
    REAL(q),   ALLOCATABLE :: POTFAK(: __omp_arg(:))   ! 1/(G+dk)**2 (G)
    REAL(q)                :: EXX,FD
    TYPE( rotation_handle), POINTER :: ROT_HANDLE
    TYPE (wavespin) WHF
    COMPLEX(q) :: CWORK(W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)
    INTEGER, SAVE :: NWARN
    REAL(q) :: WEIGHT_Q
    REAL(q) :: EXHF
    INTEGER IBLK,NDO,NBLK,NBSTART,NBSTOP,NBSTORE
!$  INTEGER NSTRIP,NSTRIP_ACT,i

    GDEF, ALLOCATABLE :: CTMP1(:),CTMP2(:)
!$  INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM

    LOGICAL :: LCOMMENSURATE_STRIP

    PROFILING_START('fock_acc')

    IF (AEXX==0) THEN
       CH=0
       PROFILING_STOP('fock_acc')
       RETURN
    ENDIF
#ifdef MPI
    IF (MOD(NK-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1 .AND. NWARN<10) THEN
       NWARN=NWARN+1
       WRITE(*,*) 'fock contribution at unconventional k-point'
    END IF
#endif
! use temporarily another WDES
    WHF=W
    WHF%WDES => WDES_FOCK

    CALL CHECK_FULL_KPOINTS
    NULLIFY(ROT_HANDLE)

    IF (PRESENT(LCOMMENSURATE)) THEN
       LCOMMENSURATE_STRIP=LCOMMENSURATE
    ELSE
       LCOMMENSURATE_STRIP=.TRUE.
    ENDIF

    IF (LCOMMENSURATE_STRIP) THEN
       NDO=NSTRIPN*W%WDES%NB_PAR
    ELSE
       NDO=NSTRIPN
    ENDIF

    NBLK=NDO
    IF (NBLOCK_FOCK>0) NBLK=MIN(NBLK,NBLOCK_FOCK)

! allocate workspace to determine the action of the Hamiltonian
! on a strip of NDO bands in blocks of size NBLK
    CALL FOCK_ACC_ALLOCATE

! average electrostatic potential for k=k' and n=n'
    FSG=FSG_STORE(NK)
 
    IF (MCALPHA==0) THEN
       IF (L_MODEL_HF) THEN
          FSG_AEXX=FSG
       ELSE
          FSG_AEXX=FSG*AEXX
       ENDIF
! there are two in principle equivalent ways to perform
! the convergence correction
! either FSG_AEXX or FSG must be set
       IF (PRESENT(WIN_ORIG)) THEN
! numerically slightly less accurate version
! applies to any orbital
          FSG_AEXX=0
       ELSE
! this one is more accurate but works only for wavefunctions in the 
! set of wavefunctions spanned by W
          FSG=0
       ENDIF
    ELSE
! k=0 term in potfak must be free of other kinds of finite size errors
       FSG=0
       FSG_AEXX=0
    ENDIF
!==========================================================================
! initialise variables
!==========================================================================
    EXHF=0
    block: DO IBLK=1,NDO,NBLK

    NGLB=MIN(NDO-IBLK+1,NBLK)
    ! global index of the first and last band in this block, and the
    ! base index of the first band in the storeback array
    IF (LCOMMENSURATE_STRIP) THEN
       NBSTART=(NPOS-1)*W%WDES%NB_PAR+IBLK
       NBSTORE=IBLK
    ELSE
       NBSTART=NPOS+IBLK-1
       NBSTORE=NBSTART
    ENDIF
    NBSTOP=NBSTART+NGLB-1

    CDIJ=0; CDLM=0; CXI=0; CKAPPA=0
!==========================================================================
! fourier transform the bands for which the HF exchange need to be calculated
! to real space, then gather to WIN
!==========================================================================
    IF (PRESENT(WIN_ORIG)) THEN
!      CALL W1_GATHER_STRIP( WHF, IBLK, IBLK+NGLB-1, WIN_ORIG, WIN)
       CALL W1_IGATHER_STRIP( WHF, IBLK, IBLK+NGLB-1, WIN_ORIG, WIN)
    ELSE
!      CALL W1_GATHER_GLB_( WHF, NBSTART, NBSTOP, ISP, WIN)
       CALL W1_IGATHER_GLB( WHF, NBSTART, NBSTOP, ISP, WIN)
    ENDIF
    NQ_USED=0
!==========================================================================
!  loop over all q-points (index NQ)
!  sum_nq phi_nq mq (r') \int phi_nq mq(r) phi_nk mk(r) / (r-r') d3r 
!==========================================================================
    qpoints: DO NQ=1,KPOINTS_FULL%NKPTS
      IF( KPOINTS_FULL%WTKPT(NQ)==0 .OR. &
         (HFKIDENT.AND.SKIP_THIS_KPOINT_IN_FOCK(WHF%WDES%VKPT(:,NQ))) .OR. &
         (.NOT.HFKIDENT.AND.SKIP_THIS_KPOINT_IN_FOCK(KPOINTS_FULL%VKPT(:,NQ)-WHF%WDES%VKPT(:,NK)))) CYCLE
       NQ_USED=NQ_USED+1
       WEIGHT_Q=1
       IF (ALLOCATED(WEIGHT_K_POINT_PAIR_SMALL_GROUP) .AND. PRESENT(LSYMGRAD) ) THEN
          IF (LSYMGRAD) THEN
             IF (WEIGHT_K_POINT_PAIR_SMALL_GROUP(NK,NQ)==0) CYCLE
             WEIGHT_Q=WEIGHT_K_POINT_PAIR_SMALL_GROUP(NK,NQ)
          ENDIF
       ENDIF
          
       CALL SETWDES(WHF%WDES,WDESQ,NQ)
       CALL SETWDES(WHF%WDES,WDESQ_IRZ,KPOINTS_FULL%NEQUIV(NQ))

       ISP_IRZ=ISP
       IF (KPOINTS_FULL%SPINFLIP(NQ)==1) THEN
          ISP_IRZ=3-ISP
       ENDIF

   ! set POTFAK for this q and k point
       CALL SET_GFAC(GRIDHF,LATT_CUR,NK,NQ,FSG,POTFAK)

   ! loop over bands mq (occupied bands for present q-point on the local CPU)
       mband: DO MQ=1,WHF%WDES%NBANDS
          IF (ABS(WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ))<=1E-10) CYCLE mband
          IF ((MQ-1)*W%WDES%NB_PAR+W%WDES%NB_LOW<NBANDSGWLOW_FOCK) CYCLE mband


          IF (NQ<=WHF%WDES%NKPTS) THEN
             CALL W1_COPY(ELEMENT(WHF, WDESQ, MQ, ISP), WQ)
             CALL FFTWAV_W1(WQ)
          ELSE
#ifndef gammareal
           ! 
           ! symmetry must be considered if the wavefunctions for this
           ! k-point NQ (containing all k-points in the entire BZ) 
           ! are not stored in W
           !
             LSHIFT=.FALSE.
             IF ((ABS(KPOINTS_FULL%TRANS(1,NQ))>TINY) .OR. &
                 (ABS(KPOINTS_FULL%TRANS(2,NQ))>TINY) .OR. &
                 (ABS(KPOINTS_FULL%TRANS(3,NQ))>TINY)) LSHIFT=.TRUE.
             CALL W1_ROTATE_AND_FFT(WQ, ELEMENT(WHF, WDESQ_IRZ, MQ, ISP_IRZ), &
                  ROT_HANDLE, P, LATT_CUR, LSHIFT)
#endif
          ENDIF
!-----------------------------------------------------------------------------
! calculate fock potential and add to accelerations
!-----------------------------------------------------------------------------
      ! calculate charge phi_q nq(r) phi_k nk(r)
          CALL FOCK_CHARGE_MU( WIN(1:NGLB), WQ, GWORK, CRHOLM)

          WEIGHT=WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ)/GRIDHF%NPLWV*WEIGHT_Q
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP IF (MOD(NGLB,omp_nthreads)==0) &
!$OMP PRIVATE(N,EXX,NS,FD) &
!$OMP REDUCTION(+:EXHF)
          nband: DO N=1,NGLB
!$           i=1
             NS=NBSTART+N-1        ! global storage index of present band
      ! fft to reciprocal space
             CALL FFT3D(GWORK(1,N),GRIDHF,-1)
      ! model-GW set GFAC state dependent
             IF (MODEL_GW==2) THEN
!$              i=OMP_GET_THREAD_NUM()+1
                CALL MODEL_GW_SET_GFAC(GRIDHF, LATT_CUR, NK, NQ, KPOINTS_FULL%NEQUIV(NQ), &
                   NS, (MQ-1)*W%WDES%NB_PAR+W%WDES%NB_LOW, ISP, ISP_IRZ, FSG, POTFAK(1 __omp_arg(i)))
             ENDIF

             IF (MCALPHA/=0) THEN
                ! with finite size corrections:
                CALL APPLY_GFAC_MULTIPOLE(GRIDHF, GWORK(1,N), POTFAK(1 __omp_arg(i)))
             ELSE
                CALL APPLY_GFAC_EXCHANGE(GRIDHF, GWORK(1,N), POTFAK(1 __omp_arg(i)), EXX)
                IF (PRESENT(EXHF_ACFDT)) THEN
                   EXX=EXX*(0.5_q/GRIDHF%NPLWV)  ! divide by grid points
                   ! correct for self-interaction 
                   IF (NS==(MQ-1)*W%WDES%NB_PAR+W%WDES%NB_LOW .AND. NQ==NK) THEN
#ifdef MPI
                      IF (W%WDES%COMM_INB%NODE_ME==1) THEN
                         EXX=EXX+FSG_AEXX*0.5_q  ! one node adds corrections
                      ENDIF
#else
                      EXX=EXX+FSG_AEXX*0.5_q
#endif
                   ENDIF
                   ! use smaller occupancy
                   ! FD=MIN(WHF%FERTOT(NS,NK,ISP),WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ))
                   ! use one-electron occupancy of state at greater energy
                   IF (REAL(WHF%CELTOT(NS,NK,ISP),q) > REAL(WHF%CELEN(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ),q)) THEN
                      FD=WHF%FERTOT(NS,NK,ISP)
                   ELSE
                      FD=WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ)
                   ENDIF

                   EXHF=EXHF-EXX &
                         *WHF%WDES%RSPIN*WHF%WDES%WTKPT(NK)*WEIGHT_Q* &
                        (FD- & 
                         WHF%FERTOT(NS,NK,ISP)*WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ))
                ENDIF
             ENDIF

      ! back to real space to get  \int phi_q(r) phi_k(r) / (r-r') d3r
             CALL FFT3D(GWORK(1,N),GRIDHF,1)

      ! add to acceleration xi in real space
             CALL VHAMIL_TRACE(WDESK, GRID_FOCK, GWORK(1,N), WQ%CR(1), CXI(1,N), WEIGHT)
          ENDDO nband
!$OMP END PARALLEL DO
          IF (WHF%WDES%LOVERL) THEN
             ! calculate D_LM
             ! build the descriptor for RPRO1
             DO N=1,NGLB
                W1(N)%CPROJ => CDLM(:,N)
             ENDDO

             AUG_DES%RINPL=WEIGHT ! multiplicator for RPRO1
             CALL RPROMU_HF(FAST_AUG_FOCK, AUG_DES, W1, NGLB, GWORK(1,1), SIZE(GWORK,1))

!$OMP PARALLEL DO DEFAULT(SHARED) IF (MOD(NGLB,omp_nthreads)==0) PRIVATE(N,CDIJ)
             DO N=1,NGLB
                IF (WHF%WDES%NRSPINORS==2) CDLM(AUG_DES%NPRO+1:AUG_DES%NPRO*2, N)=CDLM(1:AUG_DES%NPRO, N)
      ! transform D_LM -> D_lml'm'
                CALL CALC_DLLMM_TRANS(WHF%WDES, AUG_DES, TRANS_MATRIX_FOCK, CDIJ, CDLM(:,N))
      ! add D_lml'm' to kappa_lm_N (sum over l'm')
                CALL OVERL_FOCK(WHF%WDES, LMDIM, CDIJ(1,1,1,1), WQ%CPROJ(1), CKAPPA(1,N),.TRUE.)
             ENDDO
!$OMP END PARALLEL DO
          ENDIF
       ENDDO mband
    ENDDO qpoints
    IF (((ODDONLY.OR. EVENONLY) .AND. NQ_USED*2 /=KPOINTS_FULL%NKPTS_NON_ZERO) .OR. &
        (.NOT. (ODDONLY.OR. EVENONLY) .AND. NQ_USED*NKREDX*NKREDY*NKREDZ /=KPOINTS_FULL%NKPTS_NON_ZERO)) THEN
       CALL vtutor%bug("internal error in FOCK_ACC: number of k-points incorrect " // str(NQ_USED) &
          // " " // str(KPOINTS_FULL%NKPTS_NON_ZERO) // " " // str(NKREDX) // " " // str(NKREDY) // &
          " " // str(NKREDZ), __FILE__, __LINE__)
    ENDIF
!-----------------------------------------------------------------------------
! end of main fock loop                                         
!-----------------------------------------------------------------------------
! collect CXI and CKAPPA
    IF (WHF%WDES%DO_REDIS) THEN
!      PROFILING_START('m_sum')
!      CALLMPI( M_sum_z(WDESK%COMM_INTER,CXI(1,1),GRID_FOCK%MPLWV*WHF%WDES%NRSPINORS*NGLB))
!      CALLMPI( M_sum_g(WDESK%COMM_INTER,CKAPPA(1,1),WHF%WDES%NPROD*NGLB))
!      PROFILING_STOP('m_sum')
       CALL W1_IREDUCE_GLB(WDESK,WXI,NBSTART,NBSTOP)
    END IF

! generate the descriptor for the original WDES
! only difference to WDESK is the FFT mesh (GRID and not GRIDHF)
    CALL SETWDES(W%WDES,WDESQ,NK)

! fourier transform local accelerations xi (only own bands)
!!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(N,NS,NT,CTMP1,CTMP2,CWORK,WEIGHT,ISPINOR,NP,MM) REDUCTION(+:CDCHF)
    fft_back:DO N=1,NGLB
       IF (MOD(NBSTART+N-1-1,W%WDES%NB_PAR)+1/=W%WDES%NB_LOW) CYCLE fft_back

       ! local storage indices
       NS=(NBSTART+N-1-1)/W%WDES%NB_PAR+1 ! source
       NT=(NBSTORE+N-1-1)/W%WDES%NB_PAR+1 ! target

!!!$     IF (MODULO(NT,OMP_GET_NUM_THREADS())/=OMP_GET_THREAD_NUM()) CYCLE fft_back

       CH(:,NT)=0

! add CKAPPA to CXI (full acceleration on band N now in CXI)
       IF (WHF%WDES%LOVERL) THEN
          IF (.NOT.PRESENT(WIN_ORIG)) THEN
             ! convergence correction non local part
             CTMP1(:)=FSG_AEXX*WHF%CPROJ(:,NS,NK,ISP)*WHF%FERWE(NS,NK,ISP)
             CALL OVERL1(WDESK, LMDIM, CQIJ(1,1,1,1), CQIJ(1,1,1,1), 0.0_q, CTMP1(1), CTMP2(1))
             CKAPPA(:,N)=CKAPPA(:,N)+CTMP2(:)
          ENDIF

          IF (NONLR_S%LREAL) THEN
             CWORK=0
             CALL RACC0(NONLR_S, WDESQ, CKAPPA(1,N), CWORK(1))
             DO ISPINOR=0,WDESQ%NRSPINORS-1
                CALL FFTEXT(WDESQ%NGVECTOR,WDESQ%NINDPW(1), &
                     CWORK(1+ISPINOR*WDESQ%GRID%MPLWV), &
                     CH(1+ISPINOR*WDESQ%NGVECTOR,NT),WDESQ%GRID,.TRUE.)
             ENDDO
          ELSE
             CALL VNLAC0(NONL_S, WDESK, CKAPPA(1,N), CH(1,NT))
          ENDIF
       ENDIF
       IF (.NOT.PRESENT(WIN_ORIG)) THEN
       ! convergence correction non local part
          CH(:,NT)=CH(:,NT)+FSG_AEXX*WHF%CW(:,NS,NK,ISP)*WHF%FERWE(NS,NK,ISP)
       ENDIF

       ! double counting hence subtract half the self energy
       ! and change sign since we have use e^2 to calculate the potential
       WEIGHT=WHF%FERWE(NS,NK,ISP)*WHF%WDES%WTKPT(NK)*0.5_q*WHF%WDES%RSPIN

       DO ISPINOR=0,WDESK%NRSPINORS-1
          CALL FFTEXT(WDESK%NGVECTOR,WDESK%NINDPW(1), &
               CXI(1+ISPINOR*GRID_FOCK%MPLWV,N), &
               CH(1+ISPINOR*WDESK%NGVECTOR,NT),GRID_FOCK,.TRUE.)

          DO NP=1,WDESK%NGVECTOR
             MM=NP+ISPINOR*WDESK%NGVECTOR
             CH(MM,NT)=-CH(MM,NT)
             CDCHF=CDCHF-CONJG(CH(MM,NT))*WHF%CW(MM,NS,NK,ISP)*WEIGHT
            ENDDO
       ENDDO
    ENDDO fft_back
!!!$OMP END PARALLEL

    ENDDO block

    IF (PRESENT(EXHF_ACFDT)) EXHF_ACFDT=EXHF_ACFDT+EXHF

! deallocate workspace
    CALL FOCK_ACC_DEALLOCATE

    PROFILING_STOP('fock_acc')

    CONTAINS
!> Allocate memory, we determine the action of the Hamiltonian
!> on NDO bands in blocks of size NLBK
    SUBROUTINE FOCK_ACC_ALLOCATE
    INTEGER ISTT,ISTATUS
    LOGICAL, SAVE :: LFIRST=.TRUE.
#ifdef shmem_bcast_buffer
    COMPLEX(q), POINTER :: TMP1(:,:),TMP2(:,:)
    GDEF , POINTER      :: TMP3(:,:)
#endif
! Register (and write out) the memory demands of FOCK_ACC before the actual allocation.
! We include only the part that scales with NBLK. This is done only the first time FOCK_ACC is called.
    IF (LFIRST) THEN
       CALL REGISTER_ALLOCATE(( &
           8._q*ndata*m_ GRIDHF%MPLWV+ &                        ! GWORK
           8._q*ndata*2*AUG_DES%NPROD+ &                        ! CRHOLM+CDLM
          16._q*      GRID_FOCK%MPLWV+ &                        ! CXI
           8._q*ndata*WHF%WDES%NPROD+ &                         ! CKAPPA
         (16._q*      WHF%WDES%NRPLWV+ &                        ! WIN%CW
          16._q*      WHF%WDES%GRID%MPLWV*WHF%WDES%NRSPINORS+ & ! WIN%CR
           8._q*ndata*WHF%WDES%NPROD) &                         ! WIN%CPROJ
#ifdef shmem_bcast_buffer
           /WHF%WDES%COMM_intra_node%NCPU &
#endif
           )*NBLK,'fock_wrk')

!      IF (WHF%WDES%COMM%IONODE==WHF%WDES%COMM%NODE_ME) &
!         CALL DUMP_ALLOCATE_TAG(8,'FOCK_ACC')

       LFIRST=.FALSE.
    ENDIF

    ALLOCATE( &
         GWORK(m_ GRIDHF%MPLWV,NBLK), &
         CRHOLM(AUG_DES%NPROD*WHF%WDES%NRSPINORS,NBLK), &
         W1(NBLK), &
         CDIJ(LMDIM,LMDIM,WHF%WDES%NIONS,WHF%WDES%NRSPINORS), &
         CDLM(AUG_DES%NPROD*WHF%WDES%NRSPINORS,NBLK), &
         STAT=ISTT)
    ISTATUS=ISTT

!$  IF (MODEL_GW/=2) THEN
!$     ALLOCATE(POTFAK(GRIDHF%MPLWV,omp_nthreads),STAT=ISTT)
!$  ELSE
       ALLOCATE(POTFAK(GRIDHF%MPLWV __omp_arg(1)),STAT=ISTT)
!$  ENDIF
    ISTATUS=ISTATUS+ISTT

    ALLOCATE(CTMP1(WHF%WDES%NPROD),CTMP2(WHF%WDES%NPROD),STAT=ISTT)
    ISTATUS=ISTATUS+ISTT

    CALL SETWDES(WHF%WDES,WDESQ,0)
    CALL NEWWAV(WQ, WDESQ, .TRUE., ISTT)
    ISTATUS=ISTATUS+ISTT

    CALL SETWDES(WHF%WDES,WDESK,NK)
    ALLOCATE(WIN(NBLK))
#ifdef shmem_bcast_buffer
    ! when the shared memory segments are already
    ! allocated we check whether they are large enough
    IF (LFOCK_SHM_ALLOC) THEN
       IF (FOCK_SHM_BCST1%desc%dim2<NBLK.OR. &
           FOCK_SHM_BCST2%desc%dim2<NBLK.OR. &
           FOCK_SHM_BCST3%desc%dim2<NBLK) THEN
          ! NBLK is too large, deallocate previous buffers
          CALL M_shmem_dealloc(FOCK_SHM_BCST1)
          CALL M_shmem_dealloc(FOCK_SHM_BCST2)
          CALL M_shmem_dealloc(FOCK_SHM_BCST3)
          ! mark SHMEM bcast buffers as unset
          LFOCK_SHM_ALLOC=.FALSE.
       ENDIF
    ENDIF
    IF (.NOT.LFOCK_SHM_ALLOC) THEN
       CALL M_shmem_alloc( &
      &   WHF%WDES%COMM_intra_node,WDESK%NRPLWV,NBLK,FOCK_SHM_BCST1,TMP1,ISTAT=ISTT)
       ISTATUS=ISTATUS+ISTT

       CALL M_shmem_alloc( &
      &   WHF%WDES%COMM_intra_node,WDESK%GRID%MPLWV*WDESK%NRSPINORS,NBLK,FOCK_SHM_BCST2,TMP2,ISTAT=ISTT)
       ISTATUS=ISTATUS+ISTT

!      CALL REGISTER_ALLOCATE(16._q*(WDESK%NRPLWV+WDESK%GRID%MPLWV*WDESK%NRSPINORS)*NBLK/WHF%WDES%COMM_intra_node%NCPU,'shmem_buf')
 
       CALL M_shmem_alloc( &
      &   WHF%WDES%COMM_intra_node,WDESK%NPROD,NBLK,FOCK_SHM_BCST3,TMP3,ISTAT=ISTT)
       ISTATUS=ISTATUS+ISTT

!      CALL REGISTER_ALLOCATE( 8._q*ndata*WDESK%NPROD*NBLK/WHF%WDES%COMM_intra_node%NCPU,'shmem_buf')

       LFOCK_SHM_ALLOC=.TRUE.
    ELSE
       TMP1(1:WDESK%NRPLWV,1:NBLK) => FOCK_SHM_BCST1%fptr_c(1:WDESK%NRPLWV*NBLK)
       TMP2(1:WDESK%GRID%MPLWV*WDESK%NRSPINORS,1:NBLK) => FOCK_SHM_BCST2%fptr_c(1:WDESK%GRID%MPLWV*WDESK%NRSPINORS*NBLK)
#ifdef gammareal
       TMP3(1:WDESK%NPROD,1:NBLK) => FOCK_SHM_BCST3%fptr_r(1:WDESK%NPROD*NBLK)
#else
       TMP3(1:WDESK%NPROD,1:NBLK) => FOCK_SHM_BCST3%fptr_c(1:WDESK%NPROD*NBLK)
#endif
    ENDIF

    DO N=1,NBLK
       WIN(N)%CW   =>TMP1(:,N)
       WIN(N)%CR   =>TMP2(:,N)
       WIN(N)%CPROJ=>TMP3(:,N)

       WIN(N)%WDES1=>WDESK
    ENDDO

    NULLIFY(TMP1,TMP2,TMP3)
#else
    DO N=1,NBLK
       CALL NEWWAV(WIN(N) , WDESK, .TRUE., ISTT)
       IF (ISTT/=0) EXIT
    ENDDO
    ISTATUS=ISTATUS+ISTT
#endif
    ALLOCATE( &
         CXI(GRID_FOCK%MPLWV*WHF%WDES%NRSPINORS,NBLK), &
         CKAPPA(WHF%WDES%NPROD,NBLK), &
         STAT=ISTT)
    ISTATUS=ISTATUS+ISTT

    ALLOCATE(WXI(NBLK))
    DO N=1,NBLK
       WXI(N)%CW   =>   CXI(:,N)
       WXI(N)%CPROJ=>CKAPPA(:,N)

       WXI(N)%WDES1=>WDESK
    ENDDO

    IF (ISTATUS/=0) THEN
!      IF (WHF%WDES%COMM%IONODE==WHF%WDES%COMM%NODE_ME) &
!         CALL DUMP_ALLOCATE_TAG(8,'FOCK_ACC')
       CALL vtutor%error("FOCK_ACC_ALLOCATE: ERROR: could not allocate enough workspace. Try reducing &
          &NBLOCK_FOCK.")
    ENDIF

    RETURN
    END SUBROUTINE FOCK_ACC_ALLOCATE


!> Deallocate memory, except for the shared memory buffers.
!> The latter are deallocated at the end of the run.
    SUBROUTINE FOCK_ACC_DEALLOCATE
    LOGICAL, SAVE :: LFIRST=.TRUE.
! Register deallocation of FOCK_ACC workspace (only on first call).
! When a shared memory segment is used for WIN this will not be deallocated.
    IF (LFIRST) THEN
       CALL DEREGISTER_ALLOCATE(( &
           8._q*ndata*m_ GRIDHF%MPLWV+ &                        ! GWORK
           8._q*ndata*2*AUG_DES%NPROD+ &                        ! CRHOLM+CDLM
          16._q*      GRID_FOCK%MPLWV+ &                        ! CXI
           8._q*ndata*WHF%WDES%NPROD  &                         ! CKAPPA
#ifndef shmem_bcast_buffer
        +(16._q*      WHF%WDES%NRPLWV+ &                        ! WIN%CW
          16._q*      WHF%WDES%GRID%MPLWV*WHF%WDES%NRSPINORS+ & ! WIN%CR
           8._q*ndata*WHF%WDES%NPROD) &                         ! WIN%CPROJ
#endif
           )*NBLK,'fock_wrk')

       LFIRST=.FALSE.
    ENDIF

    DEALLOCATE(CTMP1,CTMP2)

    DEALLOCATE(GWORK,CXI,CKAPPA,CRHOLM,CDIJ,CDLM,W1,POTFAK)
    CALL DEALLOCATE_ROT_HANDLE(ROT_HANDLE)
    CALL DELWAV(WQ,.TRUE.)
#ifdef shmem_bcast_buffer
    DO N=1,NBLK
       NULLIFY(WIN(N)%CW,WIN(N)%CR,WIN(N)%CPROJ,WIN(N)%WDES1)
    ENDDO
!    IF (LFOCK_SHM_ALLOC) THEN
!       CALL M_shmem_dealloc(FOCK_SHM_BCST1)
!       CALL M_shmem_dealloc(FOCK_SHM_BCST2)
!       CALL M_shmem_dealloc(FOCK_SHM_BCST3)
!       LFOCK_SHM_ALLOC=.FALSE.
!    ENDIF
#else
    DO N=1,NBLK
       CALL DELWAV(WIN(N) ,.TRUE.)
    ENDDO
#endif
    DEALLOCATE(WIN)

    DO N=1,NBLK
       NULLIFY(WXI(N)%CW,WXI(N)%CPROJ,WXI(N)%WDES1)
    ENDDO
    DEALLOCATE(WXI)

    RETURN
    END SUBROUTINE FOCK_ACC_DEALLOCATE

  END SUBROUTINE FOCK_ACC

#else
!
!> non blocked version
!
  SUBROUTINE FOCK_ACC(GRID_, LMDIM, LATT_CUR, W, &
       NONLR_S, NONL_S, NK, ISP, NPOS, NSTRIPN, &
       CH, P, CQIJ, CDCHF, WIN_ORIG, LSYMGRAD, EXHF_ACFDT)
    USE sym_prec
    USE nonl_high
    USE wave_high
    USE lattice
    USE constant
    USE full_kpoints
    USE pseudo
    IMPLICIT NONE

! passed variables
    INTEGER LMDIM
    TYPE (grid_3d) GRID_                          ! not used (scheduled for removal)
    TYPE (latt) LATT_CUR
    TYPE (wavespin) W
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct)  NONL_S
    TYPE (potcar)       P(:)
    INTEGER NK,ISP,NPOS,NSTRIPN
    COMPLEX(q)       :: CH(:,:)                   !< accelerations in rec. space
    OVERLAP             CQIJ(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
    COMPLEX(q) CDCHF                              !< double counting correction
    TYPE (wavefun1), OPTIONAL :: WIN_ORIG(:)      !< original wavefunctions distributed over bands
    LOGICAL,OPTIONAL :: LSYMGRAD
    REAL(q),OPTIONAL :: EXHF_ACFDT                !< difference between HF exchange, and ACFDT exchange
                                                  !< compare equation (12) in Harl et al. PRB 81 115126 (2010)
! local variables
    TYPE (wavedes1), TARGET :: WDESK, WDESQ, WDESQ_IRZ
    TYPE (wavefun1) :: W1, WQ
    TYPE (wavefun1),ALLOCATABLE :: WIN(:)
    REAL(q) :: WEIGHT
    REAL(q) :: FSG                                !< singularity correction setting V(G=0) to a finite value
    REAL(q) :: FSG_AEXX                           !< alternative treatment, by adding original orbital
    INTEGER ISPINOR
    INTEGER N, N_, MQ, NP, NGLB, MM, ISP_IRZ
    INTEGER NQ, NQ_USED
    LOGICAL LSHIFT
    GDEF, ALLOCATABLE      :: GWORK(:)         !< fock pot in real sp
    COMPLEX(q),ALLOCATABLE :: CXI(:,:)         !< acc. in real space
    GDEF,      ALLOCATABLE :: CKAPPA(:,:)      !< stores NL accelerations
    GDEF,      ALLOCATABLE :: CRHOLM(:)        !< augmentation occupancy matrix
    GDEF,      ALLOCATABLE :: CDIJ(:,:,:,:)    !< D_lml'm'
    GDEF,ALLOCATABLE,TARGET:: CDLM(:)          !< D_LM
    REAL(q),   ALLOCATABLE :: POTFAK(:)        !< 1/(G+dk)**2 (G)
    REAL(q)                :: EXX,FD
    TYPE( rotation_handle), POINTER :: ROT_HANDLE
    TYPE (wavespin) WHF
    COMPLEX(q) :: CWORK(W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)
    INTEGER, SAVE :: NWARN
    REAL(q)    :: WEIGHT_Q

    IF (AEXX==0) THEN
       CH=0
       RETURN
    ENDIF
#ifdef MPI
    IF (MOD(NK-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1 .AND. NWARN<10) THEN
       NWARN=NWARN+1
       WRITE(*,*) 'fock contribution at unconventional k-point'
    END IF
#endif
! use temporarily another WDES
    WHF=W
    WHF%WDES => WDES_FOCK

    CALL CHECK_FULL_KPOINTS
    NULLIFY(ROT_HANDLE)
! allocate memory, we have to do the action of the Hamiltonian acceleration on nstripn bands
! using strips of size m for the second band
    NGLB=NSTRIPN*W%WDES%NB_PAR

    ALLOCATE( &
         GWORK(m_ GRIDHF%MPLWV), &
         CXI(GRID_FOCK%MPLWV*WHF%WDES%NRSPINORS,NGLB), &
         CKAPPA(WHF%WDES%NPROD,NGLB),                   &
         CRHOLM(AUG_DES%NPROD*WHF%WDES%NRSPINORS), &
         CDIJ(LMDIM,LMDIM,WHF%WDES%NIONS,WHF%WDES%NRSPINORS), &
         CDLM(AUG_DES%NPROD*WHF%WDES%NRSPINORS), &
         POTFAK(GRIDHF%MPLWV))

    CALL SETWDES(WHF%WDES,WDESQ,0)
    CALL NEWWAV(WQ , WDESQ,.TRUE.)

    CALL SETWDES(WHF%WDES,WDESK,NK)
    ALLOCATE(WIN(NGLB))
    DO N=1,NGLB
       CALL NEWWAV(WIN(N) , WDESK, .TRUE.)
    ENDDO

! average electrostatic potential for k=k' and n=n'
    FSG=FSG_STORE(NK)
 
    IF (MCALPHA==0) THEN
       IF (L_MODEL_HF) THEN
          FSG_AEXX=FSG
       ELSE
          FSG_AEXX=FSG*AEXX
       ENDIF
! there are two in principle equivalent ways to perform
! the convergence correction
! either FSG_AEXX or FSG must be set
       IF (PRESENT(WIN_ORIG)) THEN
! numerically slightly less accurate version
! applies to any orbital
          FSG_AEXX=0
       ELSE
! this one is more accurate but works only for wavefunctions in the 
! set of wavefunctions spanned by W
          FSG=0
       ENDIF
       !     k=0 term in potfak must be free of other kinds of finite size errors
    ELSE
       FSG=0
       FSG_AEXX=0
    ENDIF
!==========================================================================
! initialise variables
!==========================================================================
    CDIJ=0; CDLM=0; CXI=0; CKAPPA=0; CH=0
!==========================================================================
! fourier transform the bands for which the HF exchange need to be calculated
! to real space, then gather to WIN
!==========================================================================
    IF (PRESENT(WIN_ORIG)) THEN
       CALL W1_GATHER_W1( WHF, NSTRIPN, WIN_ORIG, WIN)
    ELSE
       CALL W1_GATHER( WHF, NPOS, NPOS+NSTRIPN-1, ISP, WIN)
    ENDIF
    NQ_USED=0
!==========================================================================
!  loop over all q-points (index NQ)
!  sum_nq phi_nq mq (r') \int phi_nq mq(r) phi_nk mk(r) / (r-r') d3r 
!==========================================================================
    qpoints: DO NQ=1,KPOINTS_FULL%NKPTS
      IF( KPOINTS_FULL%WTKPT(NQ)==0 .OR. &
         (HFKIDENT.AND.SKIP_THIS_KPOINT_IN_FOCK(WHF%WDES%VKPT(:,NQ))) .OR. &
         (.NOT.HFKIDENT.AND.SKIP_THIS_KPOINT_IN_FOCK(KPOINTS_FULL%VKPT(:,NQ)-WHF%WDES%VKPT(:,NK)))) CYCLE
       NQ_USED=NQ_USED+1
       WEIGHT_Q=1
       IF (ALLOCATED(WEIGHT_K_POINT_PAIR_SMALL_GROUP) .AND. PRESENT(LSYMGRAD) ) THEN
          IF (LSYMGRAD) THEN
             IF (WEIGHT_K_POINT_PAIR_SMALL_GROUP(NK,NQ)==0) CYCLE
             WEIGHT_Q=WEIGHT_K_POINT_PAIR_SMALL_GROUP(NK,NQ)
          ENDIF
       ENDIF
          
       CALL SETWDES(WHF%WDES,WDESQ,NQ)
       CALL SETWDES(WHF%WDES,WDESQ_IRZ,KPOINTS_FULL%NEQUIV(NQ))

       ISP_IRZ=ISP
       IF (KPOINTS_FULL%SPINFLIP(NQ)==1) THEN
          ISP_IRZ=3-ISP
       ENDIF

   ! set POTFAK for this q and k point
       CALL SET_GFAC(GRIDHF,LATT_CUR,NK,NQ,FSG,POTFAK)

   ! loop over bands mq (occupied bands for present q-point on the local CPU)
       mband: DO MQ=1,WHF%WDES%NBANDS
          IF (ABS(WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ))<=1E-10) CYCLE mband
          IF ((MQ-1)*W%WDES%NB_PAR+W%WDES%NB_LOW<NBANDSGWLOW_FOCK) CYCLE mband


          IF (NQ<=WHF%WDES%NKPTS) THEN
             CALL W1_COPY(ELEMENT(WHF, WDESQ, MQ, ISP), WQ)
             CALL FFTWAV_W1(WQ)
          ELSE
#ifndef gammareal
           ! 
           ! symmetry must be considered if the wavefunctions for this
           ! k-point NQ (containing all k-points in the entire BZ) 
           ! are not stored in W
           !
             LSHIFT=.FALSE.
             IF ((ABS(KPOINTS_FULL%TRANS(1,NQ))>TINY) .OR. &
                 (ABS(KPOINTS_FULL%TRANS(2,NQ))>TINY) .OR. &
                 (ABS(KPOINTS_FULL%TRANS(3,NQ))>TINY)) LSHIFT=.TRUE.
             CALL W1_ROTATE_AND_FFT(WQ, ELEMENT(WHF, WDESQ_IRZ, MQ, ISP_IRZ), &
                  ROT_HANDLE, P, LATT_CUR, LSHIFT)
#endif
          ENDIF
!-----------------------------------------------------------------------------
! calculate fock potential and add to accelerations
!-----------------------------------------------------------------------------
      ! calculate charge phi_q nq(r) phi_k nk(r)
          nband: DO N=1,NGLB
             CALL FOCK_CHARGE( WIN(N), WQ, GWORK, CRHOLM)
      ! fft to reciprocal space
             CALL FFT3D(GWORK(1),GRIDHF,-1)
      ! model-GW set GFAC state dependent
             IF (MODEL_GW==2) &
                  CALL MODEL_GW_SET_GFAC(GRIDHF, LATT_CUR, NK, NQ, KPOINTS_FULL%NEQUIV(NQ), &
                  (NPOS-1)*W%WDES%NB_PAR+N, (MQ-1)*W%WDES%NB_PAR+W%WDES%NB_LOW, ISP, ISP_IRZ,FSG, POTFAK)
             IF (MCALPHA/=0) THEN
                ! with finite size corrections:
                CALL APPLY_GFAC_MULTIPOLE(GRIDHF, GWORK(1), POTFAK(1))
             ELSE
                CALL APPLY_GFAC_EXCHANGE(GRIDHF, GWORK(1), POTFAK(1), EXX)
                IF (PRESENT(EXHF_ACFDT)) THEN
                   EXX=EXX*(0.5_q/GRIDHF%NPLWV)  ! divide by grid points
                   N_=(NPOS-1)*WHF%WDES%NB_PAR+N ! global storage index of present band
                   ! correct for self-interaction 
                   IF (N_==(MQ-1)*W%WDES%NB_PAR+W%WDES%NB_LOW .AND. NQ==NK) THEN
#ifdef MPI
                      IF (W%WDES%COMM_INB%NODE_ME==1) THEN
                         EXX=EXX+FSG_AEXX*0.5_q  ! one node adds corrections
                      ENDIF
#else
                      EXX=EXX+FSG_AEXX*0.5_q
#endif
                   ENDIF
                   ! use smaller occupancy
                   ! FD=MIN(WHF%FERTOT(N_,NK,ISP),WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ))
                   ! use one-electron occupancy of state at greater energy
                   IF (REAL(WHF%CELTOT(N_,NK,ISP),q) > REAL(WHF%CELEN(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ),q)) THEN
                      FD=WHF%FERTOT(N_,NK,ISP)
                   ELSE
                      FD=WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ)
                   ENDIF

                   EXHF_ACFDT=EXHF_ACFDT-EXX &
                         *WHF%WDES%RSPIN*WHF%WDES%WTKPT(NK)*WEIGHT_Q* &
                        (FD- & 
                         WHF%FERTOT(N_,NK,ISP)*WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ))
                ENDIF
             ENDIF

      ! back to real space to get  \int phi_q(r) phi_k(r) / (r-r') d3r
             CALL FFT3D(GWORK(1),GRIDHF,1)

      ! add to acceleration xi in real space
             WEIGHT=WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ)/GRIDHF%NPLWV*WEIGHT_Q
             CALL VHAMIL_TRACE(WDESK, GRID_FOCK, GWORK(1), WQ%CR(1), CXI(1,N), WEIGHT)
                
             IF (WHF%WDES%LOVERL) THEN
! add to acceleration kappa
      ! calculate D_LM
                ! build the descriptor for RPRO1
                W1%CPROJ => CDLM(:)
                AUG_DES%RINPL=WEIGHT ! multiplicator for RPRO1
                CALL RPRO1_HF(FAST_AUG_FOCK,AUG_DES, W1, GWORK)
                IF (WHF%WDES%NRSPINORS==2) CDLM(AUG_DES%NPRO+1:AUG_DES%NPRO*2)=CDLM(1:AUG_DES%NPRO)
      ! transform D_LM -> D_lml'm'
                CALL CALC_DLLMM_TRANS(WHF%WDES, AUG_DES, TRANS_MATRIX_FOCK, CDIJ,CDLM)
      ! add D_lml'm' to kappa_lm_N (sum over l'm')
                CALL OVERL_FOCK(WHF%WDES, LMDIM, CDIJ(1,1,1,1), WQ%CPROJ(1), CKAPPA(1,N),.TRUE.)
             ENDIF
          ENDDO nband
       ENDDO mband
    ENDDO qpoints
    IF (((ODDONLY.OR. EVENONLY) .AND. NQ_USED*2 /=KPOINTS_FULL%NKPTS_NON_ZERO) .OR. &
        (.NOT. (ODDONLY.OR. EVENONLY) .AND. NQ_USED*NKREDX*NKREDY*NKREDZ /=KPOINTS_FULL%NKPTS_NON_ZERO)) THEN
       CALL vtutor%bug("internal error in FOCK_ACC: number of k-points incorrect " // str(NQ_USED) &
          // " " // str(KPOINTS_FULL%NKPTS_NON_ZERO) // " " // str(NKREDX) // " " // str(NKREDY) // &
          " " // str(NKREDZ), __FILE__, __LINE__)
    ENDIF
!-----------------------------------------------------------------------------
! end of main fock loop                                         
!-----------------------------------------------------------------------------
! collect CXI and CKAPPA
    IF (WHF%WDES%DO_REDIS) THEN
       CALLMPI( M_sum_z(WDESK%COMM_INTER,CXI(1,1),GRID_FOCK%MPLWV*WHF%WDES%NRSPINORS*NGLB))
       CALLMPI( M_sum_g(WDESK%COMM_INTER,CKAPPA(1,1),WHF%WDES%NPROD*NGLB))
    END IF

! generate the descriptor for the original WDES
! only difference to WDESK is the FFT mesh (GRID and not GRIDHF)
    CALL SETWDES(W%WDES,WDESQ,NK)

! fourier transform local accelerations xi (only own bands)
    fft_back:DO N=W%WDES%NB_LOW,NGLB,W%WDES%NB_PAR
       N_=(N-1)/W%WDES%NB_PAR+1

! add CKAPPA to CXI (full acceleration on band N now in CXI)
       IF (WHF%WDES%LOVERL) THEN
          IF (.NOT.PRESENT(WIN_ORIG)) THEN
             ! convergence correction non local part
             WQ%CPROJ(:)=FSG_AEXX*WHF%CPROJ(:,N_+NPOS-1,NK,ISP)*WHF%FERWE(N_+NPOS-1,NK,ISP)
             CALL OVERL1(WDESK, LMDIM, CQIJ(1,1,1,1), CQIJ(1,1,1,1), 0.0_q, WQ%CPROJ(1), WIN(1)%CPROJ(1))
             CKAPPA(:,N)=CKAPPA(:,N)+WIN(1)%CPROJ
          ENDIF

          IF (NONLR_S%LREAL) THEN
             CWORK=0
             CALL RACC0(NONLR_S, WDESQ, CKAPPA(1,N), CWORK(1))
             DO ISPINOR=0,WDESQ%NRSPINORS-1
                CALL FFTEXT(WDESQ%NGVECTOR,WDESQ%NINDPW(1), &
                     CWORK(1+ISPINOR*WDESQ%GRID%MPLWV), &
                     CH(1+ISPINOR*WDESQ%NGVECTOR,N_),WDESQ%GRID,.TRUE.)
             ENDDO
          ELSE
             CALL VNLAC0(NONL_S, WDESK, CKAPPA(1,N), CH(1,N_))
          ENDIF
       ENDIF
       IF (.NOT.PRESENT(WIN_ORIG)) THEN
       ! convergence correction non local part
          CH(:,N_)=CH(:,N_)+FSG_AEXX*WHF%CW(:,N_+NPOS-1,NK,ISP)*WHF%FERWE(N_+NPOS-1,NK,ISP)
       ENDIF

       ! double counting hence subtract half the self energy
       ! and change sign since we have use e^2 to calculate the potential
       WEIGHT=WHF%FERWE(N_+NPOS-1,NK,ISP)*WHF%WDES%WTKPT(NK)*0.5_q*WHF%WDES%RSPIN

       DO ISPINOR=0,WDESK%NRSPINORS-1
          CALL FFTEXT(WDESK%NGVECTOR,WDESK%NINDPW(1), &
               CXI(1+ISPINOR*GRID_FOCK%MPLWV,N), &
               CH(1+ISPINOR*WDESK%NGVECTOR,N_),GRID_FOCK,.TRUE.)

          DO NP=1,WDESK%NGVECTOR
             MM=NP+ISPINOR*WDESK%NGVECTOR
             CH(MM,N_)=-CH(MM,N_)
             CDCHF=CDCHF-CONJG(CH(MM,N_))*WHF%CW(MM,N_+NPOS-1,NK,ISP)*WEIGHT
            ENDDO
       ENDDO
    ENDDO fft_back
    DEALLOCATE(GWORK,CXI,CKAPPA,CRHOLM,CDIJ,CDLM,POTFAK)
    CALL DEALLOCATE_ROT_HANDLE(ROT_HANDLE)
    CALL DELWAV(WQ,.TRUE.)
    DO N=1,NGLB
       CALL DELWAV(WIN(N) ,.TRUE.)
    ENDDO
    DEALLOCATE(WIN)

  END SUBROUTINE FOCK_ACC
#endif

END MODULE fock

!************************ SUBROUTINE FOCK_CHARGE_NOINT *****************
! 
!> calculate the total charge density (i.e. plane wave +
!> fast augmentation) from two wavefunctions  w1 x w2^*
!> function without explicit interface
! 
!**********************************************************************

  SUBROUTINE FOCK_CHARGE_NOINT( W1, W2, GCHG, CRHOLM, NCRHOLM)
    USE fock
    TYPE (wavefun1) :: W1, W2  !< two wavefunctions
    GDEF ::  GCHG(*)           !< accumulated charged
    INTEGER :: NCRHOLM         !< size of CRHOLM
    GDEF ::  CRHOLM(NCRHOLM)   !< work array

    PROFILING_START('fock_charge_noint')

    ! get the plane wave contribution to the charge
    CALL PW_CHARGE_TRACE(W1%WDES1, GCHG(1), W1%CR(1), W2%CR(1))

    ! add augmentation part to charge (if required)
    IF (W1%WDES1%LOVERL) THEN
       CALL DEPSUM_TWO_BANDS_RHOLM_TRACE(W1%CPROJ(:),W2%CPROJ(:), W1%WDES1, AUG_DES, &
            TRANS_MATRIX_FOCK, CRHOLM, 1._q, W1%WDES1%LOVERL)
       
       AUG_DES%RINPL=1._q ! multiplicator used by RACC0
       CALL RACC0_HF(FAST_AUG_FOCK, AUG_DES, CRHOLM(1), GCHG(1))
    ENDIF

    PROFILING_STOP('fock_charge_noint')

  END SUBROUTINE FOCK_CHARGE_NOINT

  SUBROUTINE FOCK_CHARGE_NOINT_NOAE( W1, W2, GCHG, CRHOLM, NCRHOLM)
    USE fock

    TYPE (wavefun1) :: W1, W2  !< two wavefunctions
    GDEF ::  GCHG(*)           !< accumulated charged
    INTEGER :: NCRHOLM         !< size of CRHOLM
    GDEF ::  CRHOLM(NCRHOLM)   !< work array

    ! get the plane wave contribution to the charge
    CALL PW_CHARGE_TRACE(W1%WDES1, GCHG(1), W1%CR(1), W2%CR(1))

    ! add augmentation part to charge (if required)
    IF (W1%WDES1%LOVERL) THEN
       CALL DEPSUM_TWO_BANDS_RHOLM_TRACE_NOAE(W1%CPROJ(:),W2%CPROJ(:), W1%WDES1, AUG_DES, &
            TRANS_MATRIX_FOCK, CRHOLM, 1._q, W1%WDES1%LOVERL, FAST_AUG_FOCK%LMAX )
       
       AUG_DES%RINPL=1._q ! multiplicator used by RACC0
       CALL RACC0_HF(FAST_AUG_FOCK, AUG_DES, CRHOLM(1), GCHG(1))
    ENDIF
  END SUBROUTINE FOCK_CHARGE_NOINT_NOAE

!> calculate the total charge density (i.e. plane wave +
!> fast augmentation) from two wavefunctions without conjugation of the second wavefunction
!>   w1 * w2

  SUBROUTINE FOCK_CHARGE_NO_CONJG_NOINT( W1, W2, GCHG, CRHOLM, NCRHOLM)
    USE fock

    TYPE (wavefun1) :: W1, W2  !< two wavefunctions
    GDEF ::  GCHG(*)           !< accumulated charged
    INTEGER :: NCRHOLM         !< size of CRHOLM
    GDEF ::  CRHOLM(NCRHOLM)   !< temporary

    ! get the plane wave contribution to the charge
    CALL PW_CHARGE_TRACE_NO_CONJG(W1%WDES1, GCHG(1), W1%CR(1), W2%CR(1))

    ! add augmentation part to charge (if required)
    IF (W1%WDES1%LOVERL) THEN
       CALL DEPSUM_TWO_BANDS_RHOLM_NO_CONJG(W1%CPROJ(:),W2%CPROJ(:), W1%WDES1, AUG_DES, &
            TRANS_MATRIX_FOCK, CRHOLM, 1._q, W1%WDES1%LOVERL)
       
       AUG_DES%RINPL=1._q ! multiplicator used by RACC0
       CALL RACC0_HF(FAST_AUG_FOCK, AUG_DES, CRHOLM(1), GCHG(1))
    ENDIF
  END SUBROUTINE FOCK_CHARGE_NO_CONJG_NOINT

!************************ FOCK_CHARGE_ONE_CENTER_NOINT ****************
! 
!> calculate the total charge density (i.e. plane wave +
!> fast augmentation) from two wavefunctions  w1 x w2^*
!> also set the one center charge density
!> function without explicit interface
! 
!**********************************************************************

  SUBROUTINE FOCK_CHARGE_ONE_CENTER_NOINT( W1, W2, GCHG, H, CRHO_ONE_CENTER, CRHOLM, NCRHOLM)
    USE fock
    IMPLICIT NONE

    TYPE (wavefun1) :: W1, W2  !< two wavefunctions
    GDEF ::  GCHG(*)           !< accumulated charged
    TYPE (one_center_handle) :: H
    GDEF ::  CRHO_ONE_CENTER(H%TOTAL_ENTRIES) !< accumulated one center charge
    INTEGER :: NCRHOLM         !< size of CRHOLM
    GDEF ::  CRHOLM(NCRHOLM)   !< work array

    PROFILING_START('fock_charge_one_center_noint')

    ! get the plane wave contribution to the charge
    CALL PW_CHARGE_TRACE(W1%WDES1, GCHG(1), W1%CR(1), W2%CR(1))

    ! add augmentation part to charge (if required)
    IF (W1%WDES1%LOVERL) THEN
       CALL DEPSUM_TWO_BANDS_ONE_CTR_TR(W1%CPROJ(:),W2%CPROJ(:), W1%WDES1, AUG_DES, &
            TRANS_MATRIX_FOCK, CRHOLM, H, CRHO_ONE_CENTER, 1._q, W1%WDES1%LOVERL)

       AUG_DES%RINPL=1._q ! multiplicator used by RACC0
       CALL RACC0_HF(FAST_AUG_FOCK, AUG_DES, CRHOLM(1), GCHG(1))
    ENDIF

    PROFILING_STOP('fock_charge_one_center_noint')

  END SUBROUTINE FOCK_CHARGE_ONE_CENTER_NOINT



!************************ SUBROUTINE FOCK_QDER  ***********************
!
!> this subroutine calculates the action of the derivative of
!> the Fock part of the Hamiltonian acting onto a wavefunction
!>
!>  d H_Fock / d q psi
!>
!> it is required for the calculation of optical properties
!> presently the derivative with respect to the augmentation charges
!> is not implemented, but this contribution is expected to be rather
!> small
! 
!**********************************************************************

  SUBROUTINE FOCK_QDER(GRID_, LMDIM, LATT_CUR, W, &
       NONLR_S, NONLR_D, NONL_S, NONL_D, IDIR, NK, ISP, NPOS, NSTRIPN, &
       CH, P, CQIJ)
    USE sym_prec
    USE nonl_high
    USE wave_high
    USE lattice
    USE constant
    USE full_kpoints
    USE pseudo
    USE fock
    IMPLICIT NONE

! passed variables
    INTEGER LMDIM
    TYPE (grid_3d) GRID_                       ! not used (scheduled for removal)
    TYPE (latt) LATT_CUR
    TYPE (wavespin) W
    TYPE (nonlr_struct) NONLR_S, NONLR_D       !< second one for derivatives
    TYPE (nonl_struct)  NONL_S, NONL_D         !< second one for derivatives
    TYPE (potcar)      P(:)
    INTEGER :: IDIR                            !< derivative with respect to cart. component IDIR
    INTEGER NK,ISP,NPOS,NSTRIPN
    COMPLEX(q)       :: CH(W%WDES%NRPLWV,NSTRIPN) !< accelerations in rec. space
    OVERLAP             CQIJ(LMDIM,LMDIM,W%WDES%NIONS,W%WDES%NCDIJ)
! local variables
    TYPE (wavedes1), TARGET :: WDESK, WDESQ, WDESQ_IRZ
    TYPE (wavefun1) :: W1, WQ
    TYPE (wavefun1),ALLOCATABLE :: WIN(:)
    TYPE (wavefun1),ALLOCATABLE :: WTMP(:)
    LOGICAL, ALLOCATABLE :: LDO(:)
    REAL(q) :: WEIGHT
    REAL(q) :: FSG                             ! singularity correction
    INTEGER ISPINOR
    INTEGER N, N_, MQ, NP, NGLB, MM, ISP_IRZ
    INTEGER NQ, NQ_USED
    LOGICAL LSHIFT
    GDEF    :: GWORK (m_ GRIDHF%MPLWV,W%WDES%NRSPINORS) ! fock pot in real sp
    GDEF    :: GWORK2(m_ GRIDHF%MPLWV,W%WDES%NRSPINORS) ! fock pot in real sp
    COMPLEX(q),ALLOCATABLE :: CXI(:,:)         ! acc. in real space
    GDEF,      ALLOCATABLE :: CKAPPA1(:,:)     ! stores NL accelerations
    GDEF,      ALLOCATABLE :: CKAPPA2(:,:)     ! stores NL accelerations
    GDEF,TARGET,ALLOCATABLE:: CPROJD(:,:)      ! derivatives of proj operators with respect to q
    GDEF,      ALLOCATABLE :: CRHOLM(:)        ! augmentation occupancy matrix
    GDEF,      ALLOCATABLE :: CDIJ(:,:,:,:)    ! D_lml'm'
    GDEF,ALLOCATABLE,TARGET:: CDLM(:)          ! D_LM
    REAL(q),   ALLOCATABLE :: POTFAK(:)        ! Kernel usually 1/(G+dk)**2
    REAL(q),   ALLOCATABLE :: POTFAKQ(:)       ! derivative of Kernel w.r.t. k
    TYPE( rotation_handle), POINTER :: ROT_HANDLE
    TYPE (wavespin) WHF
    COMPLEX(q) :: CWORK(W%WDES%GRID%MPLWV*W%WDES%NRSPINORS)

    IF (AEXX==0) THEN
       CH=0
       RETURN
    ENDIF

#ifdef MPI
    IF (MOD(NK-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) THEN
       CH=0
       RETURN
    END IF
#endif
! use temporarily another WDES

    WHF=W
    WHF%WDES => WDES_FOCK

    CALL CHECK_FULL_KPOINTS
    NULLIFY(ROT_HANDLE)

! allocate memory, we have to do the acceleration on nstripn bands
! using strips of size m for the second band
    NGLB=NSTRIPN*W%WDES%NB_PAR

    ALLOCATE( &
         CXI(GRID_FOCK%MPLWV*WHF%WDES%NRSPINORS,NGLB), &
         CKAPPA1(WHF%WDES%NPROD,NGLB), &
         CKAPPA2(WHF%WDES%NPROD,NGLB), &
         CPROJD (WHF%WDES%NPROD,NGLB), &
         CRHOLM(AUG_DES%NPROD*WHF%WDES%NRSPINORS), &
         CDIJ(LMDIM,LMDIM,WHF%WDES%NIONS,WHF%WDES%NRSPINORS), &
         CDLM(AUG_DES%NPROD*WHF%WDES%NRSPINORS), &
         POTFAK(GRIDHF%MPLWV),POTFAKQ(GRIDHF%MPLWV), &
         LDO(NGLB))

    CALL SETWDES(WHF%WDES,WDESQ,0)
    CALL NEWWAV(WQ , WDESQ,.TRUE.)

    CALL SETWDES(WHF%WDES,WDESK,NK)
    ALLOCATE(WIN(NGLB))
    ALLOCATE(WTMP(NGLB))
    DO N=1,NGLB
       CALL NEWWAV(WIN(N) , WDESK,.TRUE.)
       ! WTMP is identical to WIN, except for CPROJ entry
       ! which will contain the derivative of the projectors w.r.t. k
       WTMP(N)=WIN(N)
       WTMP(N)%CPROJ => CPROJD(:,N)
    ENDDO

! average electrostatic potential for k=k' and n=n'
    FSG=FSG_STORE(NK)
!==========================================================================
! initialise variables
!==========================================================================
    CDIJ=0; CDLM=0; CXI=0; CKAPPA1=0; CKAPPA2=0
!==========================================================================
! fourier transform the bands for which the HF exchange need to be calculated
! to real space, then gather to WIN
!==========================================================================
    CALL W1_GATHER( WHF, NPOS, NPOS+NSTRIPN-1, ISP, WIN)

    ! calculate the derivative of the wave function character
    ! (could be done locally and merged, but this is simpler)
    ! generate the descriptor for the original WDES
    ! only difference to WDESK is the FFT mesh (GRID and not GRIDHF)
    CALL SETWDES(W%WDES,WDESQ,NK)
    
    IF (NONLR_S%LREAL) THEN
       LDO=.TRUE.
       CALL RPROMU(NONLR_D,WDESQ,WTMP, NGLB, LDO)
    ELSE
       DO N=1,NGLB
          CALL PROJ1(NONL_D,WDESQ, WTMP(N))
       ENDDO
    ENDIF

    NQ_USED=0
!==========================================================================
!  loop over all q-points (index NQ)
!  sum_nq phi_nq mq (r') \int phi_nq mq(r) phi_nk mk(r) / (r-r') d3r 
!==========================================================================
    qpoints: DO NQ=1,KPOINTS_FULL%NKPTS
      IF( KPOINTS_FULL%WTKPT(NQ)==0 .OR. &
         (HFKIDENT.AND.SKIP_THIS_KPOINT_IN_FOCK(WHF%WDES%VKPT(:,NQ))) .OR. &
         (.NOT.HFKIDENT.AND.SKIP_THIS_KPOINT_IN_FOCK(KPOINTS_FULL%VKPT(:,NQ)-WHF%WDES%VKPT(:,NK)))) CYCLE

       NQ_USED=NQ_USED+1

       CALL SETWDES(WHF%WDES,WDESQ,NQ)
       CALL SETWDES(WHF%WDES,WDESQ_IRZ,KPOINTS_FULL%NEQUIV(NQ))
       ISP_IRZ=ISP
       IF (KPOINTS_FULL%SPINFLIP(NQ)==1) THEN
          ISP_IRZ=3-ISP
       ENDIF

   ! set POTFAK for this q and k point
       CALL SET_GFAC(GRIDHF,LATT_CUR,NK,NQ,FSG,POTFAK)
       CALL SET_GFAC_QDER(GRIDHF,LATT_CUR,NK,NQ,0.0_q,POTFAKQ, IDIR)

   ! loop over bands mq (occupied bands for present q-point on the local CPU)
       mband: DO MQ=1,WHF%WDES%NBANDS
          IF (ABS(WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ))<=1E-10) CYCLE mband
          IF ((MQ-1)*W%WDES%NB_PAR+W%WDES%NB_LOW<NBANDSGWLOW_FOCK) CYCLE mband

          IF (NQ<=WHF%WDES%NKPTS) THEN
             CALL W1_COPY(ELEMENT(WHF, WDESQ, MQ, ISP), WQ)
             CALL FFTWAV_W1(WQ)
          ELSE
#ifndef gammareal
           ! 
           ! symmetry must be considered if the wavefunctions for this
           ! k-point NQ (containing all k-points in the entire BZ) 
           ! are not stored in W
           !
             LSHIFT=.FALSE.
             IF ((ABS(KPOINTS_FULL%TRANS(1,NQ))>TINY) .OR. &
                 (ABS(KPOINTS_FULL%TRANS(2,NQ))>TINY) .OR. &
                 (ABS(KPOINTS_FULL%TRANS(3,NQ))>TINY)) LSHIFT=.TRUE.
             CALL W1_ROTATE_AND_FFT(WQ, ELEMENT(WHF, WDESQ_IRZ, MQ, ISP_IRZ), &
                  ROT_HANDLE, P, LATT_CUR, LSHIFT)
#endif
          ENDIF
!-----------------------------------------------------------------------------
! calculate fock potential and add to accelerations
!-----------------------------------------------------------------------------
          nband: DO N=1,NGLB
! contribution sum_i int Q_ij(r) V(r) d3r <p_i | phi> -> CKAPPA1
! applies only for US PP and PAW
             IF (WHF%WDES%LOVERL) THEN
                CALL FOCK_CHARGE( WIN(N), WQ, GWORK(:,1), CRHOLM)
      ! fft to reciprocal space
                CALL FFT3D(GWORK(1,1),GRIDHF,-1)
                GWORK2=GWORK      ! save for later use
      ! multiply by 4 pi e^2/G^2 and divide by # of gridpoints to obtain potential
                CALL APPLY_GFAC(GRIDHF, GWORK(1,1), POTFAK(1))
      ! back to real space to get  \int phi_q(r) phi_k(r) / (r-r') d3r
                CALL FFT3D(GWORK(1,1),GRIDHF,1)
                WEIGHT=WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ)/GRIDHF%NPLWV
      ! calculate D_LM
                ! build the descriptor for RPRO1
                W1%CPROJ => CDLM(:)
                AUG_DES%RINPL=WEIGHT ! multiplicator for RPRO1
                CALL RPRO1_HF(FAST_AUG_FOCK,AUG_DES, W1, GWORK(:,1))
                IF (WHF%WDES%NRSPINORS==2) CDLM(AUG_DES%NPRO+1:AUG_DES%NPRO*2)=CDLM(1:AUG_DES%NPRO)
      ! transform D_LM -> D_lml'm'
                CALL CALC_DLLMM_TRANS(WHF%WDES, AUG_DES, TRANS_MATRIX_FOCK, CDIJ,CDLM)
      ! add D_lml'm' to kappa_lm_N (sum over l'm')
                CALL OVERL_FOCK(WHF%WDES, LMDIM, CDIJ(1,1,1,1), WQ%CPROJ(1), CKAPPA1(1,N),.TRUE.)
             ELSE
                GWORK2=0
             ENDIF

! contribution int d Q_ij(r)/ d q V(r) d3r <p_i | phi>
! presently only derivatives with respect to the change of the projection
! operators are implemented
! what is lacking is the derivative of Q_ij(r) w.r.t q
             GWORK=0
             CALL FOCK_CHARGE_AUG( WTMP(N), WQ, GWORK(:,1), CRHOLM)
      ! fft to reciprocal space
             CALL FFT3D(GWORK(1,1),GRIDHF,-1)
      ! multiply by 4 pi e^2/G^2 and divide by # of gridpoints to obtain potential
             CALL APPLY_GFAC(GRIDHF, GWORK(1,1), POTFAK(1))
      ! add derivative of Kernel with respect to q
             CALL APPLY_GFAC_QDER(GRIDHF, GWORK2(1,1),  GWORK(1,1), POTFAKQ(1))
      ! back to real space to get  \int phi_q(r) phi_k(r) / (r-r') d3r
             CALL FFT3D(GWORK(1,1),GRIDHF,1)
      ! add to acceleration xi in real space
             WEIGHT=WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ)/GRIDHF%NPLWV
             CALL VHAMIL_TRACE(WDESK, GRID_FOCK, GWORK(1,1), WQ%CR(1), CXI(1,N), WEIGHT)
                
             IF (WHF%WDES%LOVERL) THEN
! add to acceleration kappa
      ! calculate D_LM
                ! build the descriptor for RPRO1
                W1%CPROJ => CDLM(:)
                AUG_DES%RINPL=WEIGHT ! multiplicator for RPRO1
                CALL RPRO1_HF(FAST_AUG_FOCK,AUG_DES, W1, GWORK(:,1))
                IF (WHF%WDES%NRSPINORS==2) CDLM(AUG_DES%NPRO+1:AUG_DES%NPRO*2)=CDLM(1:AUG_DES%NPRO)
      ! transform D_LM -> D_lml'm'
                CALL CALC_DLLMM_TRANS(WHF%WDES, AUG_DES, TRANS_MATRIX_FOCK, CDIJ,CDLM)
      ! add D_lml'm' to kappa_lm_N (sum over l'm')
                CALL OVERL_FOCK(WHF%WDES, LMDIM, CDIJ(1,1,1,1), WQ%CPROJ(1), CKAPPA2(1,N),.TRUE.)
             ENDIF
          ENDDO nband
       ENDDO mband
    ENDDO qpoints
!-----------------------------------------------------------------------------
! end of main fock loop                                         
!-----------------------------------------------------------------------------
! collect CXI and CKAPPA1 and CKAPPA2
    IF (WHF%WDES%DO_REDIS) THEN
       CALLMPI( M_sum_z(WDESK%COMM_INTER,CXI(1,1),GRID_FOCK%MPLWV*WHF%WDES%NRSPINORS*NGLB))
       CALLMPI( M_sum_g(WDESK%COMM_INTER,CKAPPA1(1,1),WHF%WDES%NPROD*NGLB))
       CALLMPI( M_sum_g(WDESK%COMM_INTER,CKAPPA2(1,1),WHF%WDES%NPROD*NGLB))
    END IF

! generate the descriptor for the original WDES
! only difference to WDESK is the FFT mesh (GRID and not GRIDHF)
    CALL SETWDES(W%WDES,WDESQ,NK)

    CH=0
! fourier transform local accelerations xi (only own bands)
    fft_back:DO N=W%WDES%NB_LOW,NGLB,W%WDES%NB_PAR
       N_=(N-1)/W%WDES%NB_PAR+1

! add CKAPPA1 to CXI (full acceleration on band N now in CXI)
       IF (WHF%WDES%LOVERL) THEN
          IF (NONLR_S%LREAL) THEN
            ! contribution - r_idir | p_i > C1_i
             CWORK=0
             CALL RACC0(NONLR_D, WDESQ, CKAPPA1(1,N), CWORK(1))
             CWORK=-CWORK
             CALL RACC0(NONLR_S, WDESQ, CKAPPA2(1,N), CWORK(1))

             DO ISPINOR=0,WDESQ%NRSPINORS-1
                CALL FFTEXT(WDESQ%NGVECTOR,WDESQ%NINDPW(1), &
                     CWORK(1+ISPINOR*WDESQ%GRID%MPLWV), &
                     CH(1+ISPINOR*WDESQ%NGVECTOR,N_),WDESQ%GRID,.TRUE.)
             ENDDO
          ELSE
             ! contribution - | -i d p_i/ d k> C_i
             CALL VNLAC0(NONL_D, WDESK, CKAPPA1(1,N), CH(1,N_))
             CH(:,N_)=-CH(:,N_)
             CALL VNLAC0(NONL_S, WDESK, CKAPPA2(1,N), CH(1,N_))
          ENDIF
       ENDIF

       DO ISPINOR=0,WDESK%NRSPINORS-1
          CALL FFTEXT(WDESK%NGVECTOR,WDESK%NINDPW(1), &
               CXI(1+ISPINOR*GRID_FOCK%MPLWV,N), &
               CH(1+ISPINOR*WDESK%NGVECTOR,N_),GRID_FOCK,.TRUE.)

          DO NP=1,WDESK%NGVECTOR
             MM=NP+ISPINOR*WDESK%NGVECTOR
             CH(MM,N_)=-CH(MM,N_)
            ENDDO
       ENDDO
    ENDDO fft_back

    DEALLOCATE(CXI,CKAPPA1,CKAPPA2,CPROJD,CRHOLM,CDIJ,CDLM,POTFAK,POTFAKQ,LDO)
    CALL DEALLOCATE_ROT_HANDLE(ROT_HANDLE)
    CALL DELWAV(WQ,.TRUE.)
    DO N=1,NGLB
       CALL DELWAV(WIN(N) ,.TRUE.)
    ENDDO
    DEALLOCATE(WIN,WTMP)
  END SUBROUTINE FOCK_QDER


!***********************************************************************
!
!> interface function that returns true if Hartree-Fock is used
!> for EXX-OEP certain subroutine must not include the Hartree-Fock term
!> since the local potential contains effective contribution
!> - subrot
!> - linear_optics
!> - one center terms
!
!***********************************************************************

    FUNCTION USEFOCK_CONTRIBUTION()
      USE fock
      LOGICAL USEFOCK_CONTRIBUTION
      LOGICAL, EXTERNAL :: USE_OEP_IN_GW

      USEFOCK_CONTRIBUTION=LHFCALC .AND. ((EXXOEP==0 .AND. .NOT. USE_OEP_IN_GW()) &
           .OR. LHFCALC_FORCE)

    END FUNCTION USEFOCK_CONTRIBUTION

!***********************************************************************
!
!> interface function that returns true if Hartree-Fock is used
!> for the core-valence interaction
!
!***********************************************************************

    FUNCTION USEFOCK_ONECENTER()
      USE fock
      LOGICAL USEFOCK_ONECENTER

      USEFOCK_ONECENTER = LHFCALC .OR. LHFONE

    END FUNCTION USEFOCK_ONECENTER

!***********************************************************************
!
!> interface function that returns true if Hartree-Fock is used
!> for the one center AE term (one center HF treatment only)
!
!***********************************************************************

    FUNCTION USEFOCK_AE_ONECENTER()
      USE fock
      LOGICAL USEFOCK_AE_ONECENTER

      USEFOCK_AE_ONECENTER=LHFONE

    END FUNCTION USEFOCK_AE_ONECENTER

!***********************************************************************
!
!> function to apply one-center HF type Hamiltonians
!
!***********************************************************************

    SUBROUTINE APPLY_ONE_CENTER_AEXX
      USE setexm
      USE fock
      IF (LHFONE) CALL PUSH_XC_TYPE(LEXCH, 1-AEXX, ALDAC, 1-AEXX, AGGAC, 0.0_q)

    END SUBROUTINE APPLY_ONE_CENTER_AEXX

    SUBROUTINE RESTORE_ONE_CENTER_AEXX
      USE setexm
      USE fock
      IF (LHFONE) CALL POP_XC_TYPE

    END SUBROUTINE RESTORE_ONE_CENTER_AEXX

!***********************************************************************
!
!> interface function that returns the maximum L quantum number
!> for the one centre augmentation charges
!
!***********************************************************************

    FUNCTION FOCK_LMAXONECENTER()
      USE fock
      INTEGER  FOCK_LMAXONECENTER
      FOCK_LMAXONECENTER=LMAX_FOCK

    END FUNCTION FOCK_LMAXONECENTER

!***********************************************************************
!
!> interface function that returns true if the one center
!> contributions are calculated using
!
!***********************************************************************

    FUNCTION ONE_CENTER_NMAX_FOCKAE()
      USE fock
      INTEGER ONE_CENTER_NMAX_FOCKAE

      IF (LFOCKAEDFT) THEN
         ONE_CENTER_NMAX_FOCKAE=NMAX_FOCKAE
      ELSE
         ONE_CENTER_NMAX_FOCKAE=0
      ENDIF

    END FUNCTION ONE_CENTER_NMAX_FOCKAE

!***********************************************************************
!
!> switch off HF temporarily
!
!***********************************************************************

    SUBROUTINE PUSH_FOCK()
      USE fock
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      IF (LSTACK_FOCK) THEN
         CALL vtutor%bug("internal ERROR in PUSH_FOCK: push already used", __FILE__, __LINE__)
      ENDIF
      LHFCALC_STACK=LHFCALC
      LHFCALC=.FALSE.
      LSTACK_FOCK=.TRUE.

    END SUBROUTINE PUSH_FOCK

    SUBROUTINE POP_FOCK()
      USE fock
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      IF (.NOT. LSTACK_FOCK) THEN
         CALL vtutor%bug("internal ERROR in POP_FOCK: push was not used", __FILE__, __LINE__)
      ENDIF
      LHFCALC = LHFCALC_STACK
      LSTACK_FOCK=.FALSE.

    END SUBROUTINE POP_FOCK

!********************** SUBROUTINE APPLY_GFAC   ************************
!
!>  this subroutine multiplies the charge density by the potential kernel
!>  (essentially e^2 epsilon_0 / 4 pi (G +k-q)^2)
!
!***********************************************************************

   SUBROUTINE APPLY_GFAC(GRID, CWORK, POTFAK)
#ifdef _OPENACC
     USE mopenacc_struct_def
#endif
     USE prec
     USE mgrid
     IMPLICIT NONE

     TYPE (grid_3d) GRID
     INTEGER NP
     REAL(q)     :: POTFAK(GRID%MPLWV)
     COMPLEX(q)  :: CWORK(GRID%MPLWV)

     PROFILING_START('apply_gfac')

     ! multiply by 4 pi e^2/G^2 and divide by # of gridpoints to obtain potential
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(NP) SHARED(GRID, POTFAK, CWORK)
!$ACC PARALLEL LOOP PRESENT(GRID,CWORK,POTFAK) __IF_ASYNC__
     DO NP=1,GRID%RC%NP
        CWORK(NP)=POTFAK(NP)*CWORK(NP)
     ENDDO
NOACC !$OMP END PARALLEL DO

     PROFILING_STOP('apply_gfac')

   END SUBROUTINE APPLY_GFAC


   SUBROUTINE APPLY_GFAC_WEIGHT(GRID, CWORK, POTFAK, WEIGHT)
#ifdef _OPENACC
     USE mopenacc_struct_def
#endif
     USE prec
     USE mgrid
     IMPLICIT NONE

     TYPE (grid_3d) GRID
     INTEGER NP
     REAL(q)     :: POTFAK(GRID%MPLWV), WEIGHT
     COMPLEX(q)  :: CWORK(GRID%MPLWV)

     PROFILING_START('apply_gfac_weight')

     ! multiply by 4 pi e^2/G^2 and divide by # of gridpoints to obtain potential
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(NP) SHARED(GRID, POTFAK, CWORK, WEIGHT)
!$ACC PARALLEL LOOP PRESENT(GRID,CWORK,POTFAK) __IF_ASYNC__
     DO NP=1,GRID%RC%NP
        CWORK(NP)=POTFAK(NP)*CWORK(NP)*WEIGHT
     ENDDO
NOACC !$OMP END PARALLEL DO

     PROFILING_STOP('apply_gfac_weight')

   END SUBROUTINE APPLY_GFAC_WEIGHT


!
!>  this is abatched version of APPLY_GFAC that applies
!>  POTFAK to NSIM overlap densities CWORK(:,1:NSIM)
!

   SUBROUTINE APPLY_GFAC_WEIGHT_MU(GRID, CWORK, POTFAK, WEIGHT, NSIM)
#ifdef _OPENACC
     USE mopenacc_struct_def
#endif
     USE prec
     USE mgrid
     IMPLICIT NONE

     TYPE (grid_3d) GRID
     REAL(q)     :: POTFAK(GRID%MPLWV)
     COMPLEX(q)  :: CWORK(GRID%MPLWV,NSIM)
     REAL(q)     :: WEIGHT(NSIM)
     INTEGER     :: NSIM
     ! local
     INTEGER N,NP

     PROFILING_START('apply_gfac_weight_mu')

     ! multiply by 4 pi e^2/G^2 and divide by # of gridpoints to obtain potential
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(GRID,CWORK,POTFAK,WEIGHT) __IF_ASYNC__
     DO N=1,NSIM
        DO NP=1,GRID%RC%NP
           CWORK(NP,N)=POTFAK(NP)*CWORK(NP,N)*WEIGHT(N)
        ENDDO
     ENDDO

     PROFILING_STOP('apply_gfac_weight_mu')

   END SUBROUTINE APPLY_GFAC_WEIGHT_MU

!
!> this version also additionally returns the total energy contribution
!> sum_G rho_G POTFAK(G) rho*_G
!

   SUBROUTINE APPLY_GFAC_EXCHANGE(GRID, CWORK, POTFAK, EXCHANGE)
#ifdef _OPENACC
     USE mopenacc_struct_def
#endif
     USE prec
     USE mgrid
     USE tutor, ONLY: vtutor
     IMPLICIT NONE

     TYPE (grid_3d) GRID
     INTEGER NP
     REAL(q)     :: POTFAK(GRID%MPLWV)
     COMPLEX(q)  :: CWORK(GRID%MPLWV)
     REAL(q)     :: EXCHANGE

     PROFILING_START('apply_gfac_exchange')

     ! multiply by 4 pi e^2/G^2 and divide by # of gridpoints to obtain potential
     EXCHANGE=0
     IF (GRID%REAL2CPLX) THEN
        IF (.NOT. ASSOCIATED(GRID%FFTWEIGHT)) THEN
           CALL vtutor%bug("internal error in APPLY_GFAC_EXCHANGE: FFTWEIGHT is not associated", __FILE__, __LINE__)
        ENDIF
NOACC !$OMP PARALLEL DO DEFAULT(SHARED) &
NOACC !$OMP PRIVATE(NP) &
NOACC !$OMP REDUCTION(+:EXCHANGE)
!!DIR$ SIMD REDUCTION(+:EXCHANGE)
!$ACC PARALLEL LOOP PRESENT(POTFAK,CWORK,GRID,EXCHANGE) REDUCTION(+:EXCHANGE) __IF_ASYNC__
        DO NP=1,GRID%RC%NP
           EXCHANGE=EXCHANGE+POTFAK(NP)*CWORK(NP)*CONJG(CWORK(NP))*GRID%FFTWEIGHT(NP)
           CWORK(NP)=POTFAK(NP)*CWORK(NP)
        ENDDO
NOACC !$OMP END PARALLEL DO
     ELSE
NOACC !$OMP PARALLEL DO DEFAULT(SHARED) &
NOACC !$OMP PRIVATE(NP) &
NOACC !$OMP REDUCTION(+:EXCHANGE)
!!DIR$ SIMD REDUCTION(+:EXCHANGE)
!$ACC PARALLEL LOOP PRESENT(POTFAK,CWORK,EXCHANGE) REDUCTION(+:EXCHANGE) __IF_ASYNC__
        DO NP=1,GRID%RC%NP
           EXCHANGE=EXCHANGE+POTFAK(NP)*CWORK(NP)*CONJG(CWORK(NP))
           CWORK(NP)=POTFAK(NP)*CWORK(NP)
        ENDDO
NOACC !$OMP END PARALLEL DO
     ENDIF

     PROFILING_STOP('apply_gfac_exchange')

   END SUBROUTINE APPLY_GFAC_EXCHANGE

!
!> this is a batched version of APPLY_GFAC_EXCHANGE that applies
!> POTFAK to NSIM overlap densities CWORK(:,1:NSIM)
!

   SUBROUTINE APPLY_GFAC_EXCHANGE_MU(GRID, CWORK, POTFAK, EXCHANGE, NSIM)
#ifdef _OPENACC
     USE mopenacc_struct_def
#endif
     USE prec
     USE mgrid
     USE tutor, ONLY: vtutor
     IMPLICIT NONE

     TYPE (grid_3d) GRID
     REAL(q)     :: POTFAK(GRID%MPLWV)
     COMPLEX(q)  :: CWORK(GRID%MPLWV,NSIM)
     REAL(q)     :: EXCHANGE(NSIM)
     INTEGER     :: NSIM
     ! local
     INTEGER N,NP
     REAL(q) EXTMP

     PROFILING_START('apply_gfac_exchange_mu')

     ! multiply by 4 pi e^2/G^2 and divide by # of gridpoints to obtain potential
!$ACC KERNELS PRESENT(EXCHANGE) __IF_ASYNC__
     EXCHANGE=0
!$ACC END KERNELS
     IF (GRID%REAL2CPLX) THEN
        IF (.NOT. ASSOCIATED(GRID%FFTWEIGHT)) THEN
           CALL vtutor%bug("internal error in APPLY_GFAC_EXCHANGE_MU: FFTWEIGHT is not associated", __FILE__, __LINE__)
        ENDIF
!$ACC PARALLEL LOOP GANG PRESENT(POTFAK,CWORK,GRID,EXCHANGE) PRIVATE(EXTMP) __IF_ASYNC__
        DO N=1,NSIM
           EXTMP=0
!$ACC LOOP VECTOR REDUCTION(+:EXTMP)
           DO NP=1,GRID%RC%NP
              EXTMP=EXTMP+POTFAK(NP)*CWORK(NP,N)*CONJG(CWORK(NP,N))*GRID%FFTWEIGHT(NP)
              CWORK(NP,N)=POTFAK(NP)*CWORK(NP,N)
           ENDDO
           EXCHANGE(N)=EXTMP
        ENDDO
     ELSE
!$ACC PARALLEL LOOP GANG PRESENT(POTFAK,CWORK,GRID,EXCHANGE) PRIVATE(EXTMP) __IF_ASYNC__
        DO N=1,NSIM
           EXTMP=0
!$ACC LOOP VECTOR REDUCTION(+:EXTMP)
           DO NP=1,GRID%RC%NP
              EXTMP=EXTMP+POTFAK(NP)*CWORK(NP,N)*CONJG(CWORK(NP,N))
              CWORK(NP,N)=POTFAK(NP)*CWORK(NP,N)
           ENDDO
           EXCHANGE(N)=EXTMP
        ENDDO
     ENDIF

     PROFILING_STOP('apply_gfac_exchange_mu')

   END SUBROUTINE APPLY_GFAC_EXCHANGE_MU

#ifdef _OPENACC
   SUBROUTINE APPLY_GFAC_EXCHANGE_MU_ACC(GRID, CWORK, POTFAK, EXCHANGE, NSIM)
     USE mopenacc_struct_def
     USE prec
     USE mgrid
     USE tutor, ONLY: vtutor
     IMPLICIT NONE

     TYPE (grid_3d) :: GRID
     REAL(q)        :: POTFAK(GRID%MPLWV)
     COMPLEX(q)     :: CWORK(GRID%MPLWV,NSIM)
     REAL(q)        :: EXCHANGE(NSIM)
     INTEGER        :: NSIM
     ! local
     INTEGER, PARAMETER :: BLOCKSIZE=512
     REAL(q) :: TMP(BLOCKSIZE)
     REAL(q) :: EXTMP
     INTEGER :: N,NP,IBLOCK

     PROFILING_START('apply_gfac_exchange_mu_acc')

     ! multiply by 4 pi e^2/G^2 and divide by # of gridpoints to obtain potential
!$ACC KERNELS PRESENT(EXCHANGE) __IF_ASYNC__
     EXCHANGE=0
!$ACC END KERNELS
     IF (GRID%REAL2CPLX) THEN
        IF (.NOT. ASSOCIATED(GRID%FFTWEIGHT)) THEN
           CALL vtutor%bug("internal error in APPLY_GFAC_EXCHANGE_MU_ACC: FFTWEIGHT is not associated", __FILE__, __LINE__)
        ENDIF
!$ACC PARALLEL LOOP GANG PRESENT(POTFAK,CWORK,GRID,EXCHANGE) PRIVATE(TMP) __IF_ASYNC__
        DO IBLOCK=1,GRID%RC%NP,BLOCKSIZE
!$ACC CACHE(TMP(1:BLOCKSIZE))
!$ACC LOOP VECTOR
           DO NP=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,GRID%RC%NP)
              TMP(NP-IBLOCK+1)=POTFAK(NP)
           ENDDO
!$ACC LOOP SEQ
           DO N=1,NSIM
              EXTMP=0
!$ACC LOOP VECTOR REDUCTION(+:EXTMP)
              DO NP=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,GRID%RC%NP)
                 EXTMP=EXTMP+TMP(NP-IBLOCK+1)*REAL(CWORK(NP,N)*CONJG(CWORK(NP,N)))*GRID%FFTWEIGHT(NP)
                 CWORK(NP,N)=TMP(NP-IBLOCK+1)*CWORK(NP,N)
              ENDDO
!$ACC ATOMIC UPDATE
              EXCHANGE(N)=EXCHANGE(N)+EXTMP
           ENDDO
        ENDDO
     ELSE
!$ACC PARALLEL LOOP GANG PRESENT(POTFAK,CWORK,GRID,EXCHANGE) PRIVATE(TMP) __IF_ASYNC__
        DO IBLOCK=1,GRID%RC%NP,BLOCKSIZE
!$ACC CACHE(TMP(1:BLOCKSIZE))
!$ACC LOOP VECTOR
           DO NP=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,GRID%RC%NP)
              TMP(NP-IBLOCK+1)=POTFAK(NP)
           ENDDO
!$ACC LOOP SEQ
           DO N=1,NSIM
              EXTMP=0
!$ACC LOOP VECTOR REDUCTION(+:EXTMP)
              DO NP=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,GRID%RC%NP)
                 EXTMP=EXTMP+TMP(NP-IBLOCK+1)*REAL(CWORK(NP,N)*CONJG(CWORK(NP,N)))
                 CWORK(NP,N)=TMP(NP-IBLOCK+1)*CWORK(NP,N)
              ENDDO
!$ACC ATOMIC UPDATE
              EXCHANGE(N)=EXCHANGE(N)+EXTMP
           ENDDO
        ENDDO
     ENDIF

     PROFILING_STOP('apply_gfac_exchange_mu_acc')

   END SUBROUTINE APPLY_GFAC_EXCHANGE_MU_ACC
#endif
!
!> this version contracts over a second supplied density i.e.
!> sum_G rho^1_G POTFAK(G) rho*^2_G
!

   SUBROUTINE APPLY_GFAC_EXCHANGE_2(GRID, CWORK, CWORKG, POTFAK, EXCHANGE)

     USE prec
     USE mgrid
     USE tutor, ONLY: vtutor
     IMPLICIT NONE

     TYPE (grid_3d) GRID
     INTEGER NP
     REAL(q)     :: POTFAK(GRID%MPLWV)
     COMPLEX(q)  :: CWORK(GRID%MPLWV)
     COMPLEX(q)  :: CWORKG(GRID%MPLWV)
     COMPLEX(q)  :: EXCHANGE
     ! multiply by 4 pi e^2/G^2 and divide by # of gridpoints to obtain potential
     EXCHANGE=0
     IF (GRID%REAL2CPLX) THEN
        IF (.NOT. ASSOCIATED(GRID%FFTWEIGHT)) THEN
           CALL vtutor%bug("internal error in APPLY_GFAC_EXCHANGE: FFTWEIGHT is not associated", __FILE__, __LINE__)
        ENDIF
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(NP) &
!$OMP REDUCTION(+:EXCHANGE)
        DO NP=1,GRID%RC%NP
           EXCHANGE=EXCHANGE+POTFAK(NP)*CWORK(NP)*CONJG(CWORKG(NP))*GRID%FFTWEIGHT(NP)
           CWORK(NP)=POTFAK(NP)*CWORK(NP)
        ENDDO
!$OMP END PARALLEL DO
     ELSE
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(NP) &
!$OMP REDUCTION(+:EXCHANGE)
        DO NP=1,GRID%RC%NP
           EXCHANGE=EXCHANGE+POTFAK(NP)*CWORK(NP)*CONJG(CWORKG(NP))
           CWORK(NP)=POTFAK(NP)*CWORK(NP)
        ENDDO
!$OMP END PARALLEL DO
     ENDIF
   END SUBROUTINE APPLY_GFAC_EXCHANGE_2

!********************** SUBROUTINE EXCHANGE_GFAC   *********************
!
!>  this subroutine multiplies calculates
!>  sum_G e^2 / 4 pi (G +k-q)^2)  rho_G  rho^*_G
!
!***********************************************************************

   SUBROUTINE EXCHANGE_GFAC(GRID, CWORK, POTFAK, EXCHANGE)

     USE prec
     USE mgrid
     USE tutor, ONLY: vtutor
     IMPLICIT NONE

     TYPE (grid_3d) GRID
     INTEGER NP
     REAL(q)     :: POTFAK(GRID%MPLWV)
     COMPLEX(q)  :: CWORK(GRID%MPLWV)
     REAL(q)     :: EXCHANGE

     PROFILING_START('exchange_gfac')

     ! multiply by 4 pi e^2/G^2 and divide by # of gridpoints to obtain potential
     EXCHANGE=0
     IF (GRID%REAL2CPLX) THEN
        IF (.NOT. ASSOCIATED(GRID%FFTWEIGHT)) THEN
           CALL vtutor%bug("internal error in APPLY_GFAC_EXCHANGE: FFTWEIGHT is not associated", __FILE__, __LINE__)
        ENDIF
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(NP) &
!$OMP REDUCTION(+:EXCHANGE)
        DO NP=1,GRID%RC%NP
           EXCHANGE=EXCHANGE+POTFAK(NP)*CWORK(NP)*CONJG(CWORK(NP))*GRID%FFTWEIGHT(NP)
        ENDDO
     ELSE
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(NP) &
!$OMP REDUCTION(+:EXCHANGE)
        DO NP=1,GRID%RC%NP
           EXCHANGE=EXCHANGE+POTFAK(NP)*CWORK(NP)*CONJG(CWORK(NP))
        ENDDO
!$OMP END PARALLEL DO
     ENDIF

     PROFILING_STOP('exchange_gfac')

   END SUBROUTINE EXCHANGE_GFAC

!********************** SUBROUTINE APPLY_GFAC_WAVEFUN ******************
!
!>  this subroutine multiplies the charge density by the potential kernel
!>  within a cutoff sphere defined by a plane wave cutoff
!
!***********************************************************************

   SUBROUTINE APPLY_GFAC_WAVEFUN(WDESQ, CWORK, POTFAK)

     USE prec
     USE wave
     IMPLICIT NONE

     TYPE(wavedes1) :: WDESQ
     INTEGER NQ
     REAL(q)     :: POTFAK(WDESQ%NGVECTOR)
     COMPLEX(q)  :: CWORK(WDESQ%NGVECTOR)
   ! local
     INTEGER :: NP, I

     PROFILING_START('apply_gfac_wavefun')

     DO NP=1,WDESQ%NGVECTOR
        CWORK(NP)=POTFAK(NP)*CWORK(NP)
     ENDDO

     PROFILING_STOP('apply_gfac_wavefun')

   END SUBROUTINE APPLY_GFAC_WAVEFUN

!********************** SUBROUTINE APPLY_GFAC_QDER *********************
!
!>  this subroutine multiplies the charge density by the q derivative of
!>  the potential kernel and multiplies by -i
!
!***********************************************************************

   SUBROUTINE APPLY_GFAC_QDER(GRID, CWORK2, CWORK, POTFAK)

     USE prec
     USE mgrid
     IMPLICIT NONE
      
     TYPE (grid_3d) GRID
     INTEGER NP
     REAL(q)     :: POTFAK(GRID%MPLWV)
     COMPLEX(q)  :: CWORK2(GRID%MPLWV)    !< charge density
     COMPLEX(q)  :: CWORK(GRID%MPLWV)     !< result
     REAL(q)     :: CSUM

     PROFILING_START('apply_gfac_qder')

     DO NP=1,GRID%RC%NP
        CWORK(NP)=(POTFAK(NP)*CWORK2(NP))*(0.0_q,-1.0_q)+CWORK(NP)
     ENDDO

     PROFILING_STOP('apply_gfac_qder')

   END SUBROUTINE APPLY_GFAC_QDER


!********************** SUBROUTINE APPLY_GFAC_DER **********************
!
!>  this subroutine multiplies the charge density by the potential kernel
!>  (essentially e^2 epsilon_0 / 4 pi (G +k-q)^2)
!
!***********************************************************************

   SUBROUTINE APPLY_GFAC_DER(GRID, CWORK, POTFAK, DER, WEIGHT)
#ifdef _OPENACC
     USE mopenacc_struct_def
#endif
     USE prec
     USE mgrid
     IMPLICIT NONE
      
     TYPE (grid_3d) GRID
     INTEGER NP,I
     REAL(q)     :: POTFAK(GRID%MPLWV,0:6)
     COMPLEX(q)  :: CWORK(GRID%MPLWV)
     REAL(q)     :: DER(0:6), WEIGHT
     REAL(q)     :: A

     PROFILING_START('apply_gfac_der')

!$ACC PARALLEL LOOP GANG PRESENT(GRID,POTFAK,CWORK,DER) PRIVATE(A) __IF_ASYNC__
     DO I=0,6
        A=0
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(NP) SHARED(GRID,POTFAK,CWORK,I) REDUCTION(+:A)
!$ACC LOOP VECTOR REDUCTION(+:A)
        DO NP=1,GRID%RC%NP
           A=A+POTFAK(NP,I)*CWORK(NP)*CONJG(CWORK(NP))
        ENDDO
NOACC !$OMP END PARALLEL DO
        DER(I)=DER(I)+A/GRID%NPLWV*WEIGHT
     ENDDO

     PROFILING_STOP('apply_gfac_der')

   END SUBROUTINE APPLY_GFAC_DER

!
!>  this is a batched version of APPLY_GFAC_DER that applies
!>  POTFAK to NSIM overlap densities CWORK(:,1:NSIM)
!

   SUBROUTINE APPLY_GFAC_DER_MU(GRID, CWORK, POTFAK, DER, WEIGHT, NSIM)
#ifdef _OPENACC
     USE mopenacc_struct_def
#endif
     USE prec
     USE mgrid
     IMPLICIT NONE

     TYPE (grid_3d) GRID
     INTEGER N,NP,I
     REAL(q)     :: POTFAK(GRID%MPLWV,0:6)
     COMPLEX(q)  :: CWORK(GRID%MPLWV,NSIM)
     REAL(q)     :: DER(0:6), WEIGHT(NSIM)
     INTEGER     :: NSIM
     REAL(q)     :: A

     PROFILING_START('apply_gfac_der_mu')

!$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(GRID,POTFAK,CWORK,WEIGHT,DER) PRIVATE(A) __IF_ASYNC__
     DO N=1,NSIM
        DO I=0,6
           A=0
NOACC !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
NOACC !$OMP PRIVATE(NP) SHARED(GRID,POTFAK,CWORK,I,N) REDUCTION(+:A)
!$ACC LOOP VECTOR REDUCTION(+:A)
           DO NP=1,GRID%RC%NP
              A=A+POTFAK(NP,I)*CWORK(NP,N)*CONJG(CWORK(NP,N))
           ENDDO
NOACC !$OMP END PARALLEL DO
!$ACC ATOMIC UPDATE
           DER(I)=DER(I)+A/GRID%NPLWV*WEIGHT(N)
        ENDDO
     ENDDO

     PROFILING_STOP('apply_gfac_der_mu')

   END SUBROUTINE APPLY_GFAC_DER_MU

!********************** SUBROUTINE HFPARAMETERS ************************
!
!>  this subroutine returns the parameters for hybride functionals
!>  it is called in xcspin, xcgrad, setex and radial
!
!***********************************************************************

    SUBROUTINE HFPARAMETERS(BEXX)
      USE prec
      USE fock
      IMPLICIT NONE
      REAL(q) BEXX

      BEXX=AEXX

      IF (L_MODEL_HF) THEN
         BEXX=1.0
      ENDIF

    END SUBROUTINE HFPARAMETERS

!********************** SUBROUTINE RSPARAMETERS ************************
!
!>  this subroutine returns the parameters for range-separated
!>  functionals it is called by RAD_POT_EX_HAR in pawfock.F
!
!***********************************************************************

    SUBROUTINE RSPARAMETERS(RMU, LRHF, BRSEXX, BEXX)
      USE prec
      USE fock
      IMPLICIT NONE
      REAL(q) RMU
      LOGICAL LRHF
      REAL(q) :: BRSEXX  !< amount of range seperated exchange
      REAL(q) :: BEXX    !< amount of exact exchange

      RMU=HFSCREEN
      LRHF=LRHFCALC

      IF (L_MODEL_HF) THEN
         BRSEXX=(1-AEXX)
         BEXX=AEXX
      ELSEIF ((HFSCREEN/=0.AND.(.NOT.L_THOMAS_FERMI))) THEN
         BRSEXX=AEXX
         BEXX=0.0
      ELSE
         BRSEXX=0.0
         BEXX=AEXX
      ENDIF

    END SUBROUTINE RSPARAMETERS

!********************** SUBROUTINE OVERL_FOCK *************************
!
!> calcs D_lml'm'<p_l'm'|psi_m> and adds it to CRESUL (being kappa_n)
!> (m is not explicitly in the routine but take care when calling it)
!>
!> modified from OVERL in dfast.F
!>
!> @details @ref openmp :
!> the nested loop over \"types\" + \"ions-of-type\" is replaced by
!> a loop over \"all ions\" that is distributed over all available
!> threads.
!
!**********************************************************************

      SUBROUTINE OVERL_FOCK(WDES, LMDIM, CDIJ, CPROJ , CRESUL, LADD)
      USE prec
      USE wave
      
      IMPLICIT NONE
! passed variables
      TYPE (wavedes) WDES
      GDEF CRESUL(WDES%NPROD),CPROJ(WDES%NPROD)
      LOGICAL LOVERL
      INTEGER LMDIM
      GDEF CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NRSPINORS)
      LOGICAL LADD
! local variables
      INTEGER   :: ISPINOR, NPRO
      INTEGER   :: NIS,NT,LMMAXC,NI,L,LP 
!$    INTEGER   :: NPRO2
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS

      PROFILING_START('overl_fock')

      IF (.NOT. LADD) THEN
        CRESUL=0
      ENDIF

      spinor: DO ISPINOR=0,WDES%NRSPINORS-1            
         NPRO =ISPINOR *WDES%NPRO/2
#ifndef _OPENMP
         NIS =1
         DO NT=1,WDES%NTYP
            LMMAXC=WDES%LMMAX(NT)
            IF (LMMAXC==0) GOTO 100
            DO NI=NIS,WDES%NITYP(NT)+NIS-1
!NEC$ select_vector
               DO L =1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
!NEC$ outerloop_unroll(8)
                  DO LP=1,LMMAXC
                     CRESUL(L+NPRO)=CRESUL(L+NPRO)+ &
                          CDIJ(LP,L,NI,1+ISPINOR)* &
                          CPROJ(LP+NPRO)
                  ENDDO
               ENDDO
               NPRO = LMMAXC+NPRO
            ENDDO
100         NIS = NIS+WDES%NITYP(NT)
         ENDDO
#else
!$OMP PARALLEL DO IF (WDES%NIONS>=OMP_GET_NUM_THREADS()) DEFAULT(SHARED) &
!$OMP PRIVATE(NI,NT,LMMAXC,NPRO2,L,LP)
         DO NI=1,WDES%NIONS
            NT=WDES%ITYP(NI)
            LMMAXC=WDES%LMMAX(NT)
            IF (LMMAXC==0) CYCLE
            NPRO2=WDES%LMBASE(NI)+NPRO
!!$OMP PARALLEL DO IF (LMMAXC>=OMP_GET_NUM_THREADS()) & 
!!$OMP SCHEDULE(STATIC) DEFAULT(NONE) &
!!$OMP PRIVATE(L,LP) &
!!$OMP SHARED(LMMAXC,CRESUL,NPRO,CDIJ,NI,ISPINOR,CPROJ)
            DO L =1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
               DO LP=1,LMMAXC
                  CRESUL(L+NPRO2)=CRESUL(L+NPRO2)+ &
                       CDIJ(LP,L,NI,1+ISPINOR)* &
                       CPROJ(LP+NPRO2)
               ENDDO
            ENDDO
!!$OMP END PARALLEL DO
         ENDDO
!$OMP END PARALLEL DO
#endif
      ENDDO spinor
 
      PROFILING_STOP('overl_fock')

      RETURN
    END SUBROUTINE OVERL_FOCK

!**********************************************************************
!
!> F77 kernel to calculate the contribution to the ion decomposed
!> energy in the Fock routine
!
!**********************************************************************

    SUBROUTINE ECCP_NL_FOCK(LMDIM,LMMAXC,CDIJ0,CDIJ1, &
         CPROJ,CPROJ0,CPROJ1,ENL,WEIGHT)
      USE prec
      IMPLICIT NONE
      INTEGER LMMAXC, LMDIM
      REAL(q) ENL
      REAL(q) WEIGHT
      GDEF CDIJ0(LMDIM,LMDIM),CDIJ1(LMDIM,LMDIM)
      GDEF CPROJ(LMMAXC),CPROJ0(LMMAXC),CPROJ1(LMMAXC)
    ! local
      INTEGER L, LP
      REAL(q) ENL_ADD
!$    INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS

      PROFILING_START('eccp_nl_fock')

      ENL_ADD=0
!$OMP PARALLEL DO IF (LMMAXC>=OMP_GET_NUM_THREADS()) & 
!$OMP SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP PRIVATE(L,LP) &
!$OMP SHARED(LMMAXC,CPROJ,CDIJ0,CDIJ1,CPROJ0,CPROJ1) &
!$OMP REDUCTION(+:ENL_ADD)
      DO L=1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
         DO LP=1,LMMAXC
            ENL_ADD=ENL_ADD+REAL(CPROJ(LP)*(CDIJ0(LP,L)*GCONJG(CPROJ0(L))+ &
                                            CDIJ1(LP,L)*GCONJG(CPROJ1(L))))
         ENDDO
      ENDDO
!$OMP END PARALLEL DO
      ENL=ENL+ENL_ADD*WEIGHT

      PROFILING_STOP('eccp_nl_fock')

    END SUBROUTINE ECCP_NL_FOCK

!==========================================================================
! 
!> this subroutine adds a quadratic field to the potential
!> to correct for the errors cause by the periodic FFT
!> this is a "second" order correction and applicable only in the case
!> of isolated molecules
!>
!> WORK IN PROGRESS, BUT PROPABLY NOT VERY IMPORTANT
!
!==========================================================================


    SUBROUTINE QUADRATIC_FIELD_CORRECTION(GRID, GWORK, LATT_CUR) 
      USE prec
      USE lattice
      USE mgrid
      USE constant
      IMPLICIT NONE

      TYPE (grid_3d) GRID
      GDEF           GWORK(m_ GRID%MPLWV)    !< fock pot in real sp
      TYPE (latt)    LATT_CUR
    ! local
      INTEGER IDIP, NOUT, NOUTH, IDIR, INDMIN(3), NI, NC, N1, N2, N3, II
      REAL(q) QUADFAC, QUAD_FIELD, CUTOFF, POTCOR, XX
    !
    ! a smooth cutoff function is used for potentials and charge densities
    ! around the dividing plane
    ! WIDTH must be at least 1
    ! 4 grid points around the step function are usually a good choice
      REAL(q), PARAMETER :: WIDTH=4
    
  
      QUAD_FIELD=TPI/LATT_CUR%OMEGA*FELECT/3
   ! correct the potential
      dir: DO IDIP=1,3
         NOUT=GRID%NGPTAR(IDIP)
         NOUTH=NOUT/2
         QUADFAC=QUAD_FIELD*(LATT_CUR%ANORM(IDIP)/NOUT)**2
         IDIR=IDIP
         IF (GRID%RL%NFAST==3) THEN
            IDIR=MOD(IDIP,3)+1 ! mpi version: x-> N2, y-> N3, z-> N1
         ENDIF
         INDMIN(IDIR)=GRID%NGPTAR(IDIP)/2
         NI=0
         DO NC=1,GRID%RL%NCOL
            N2=GRID%RL%I2(NC)
            N3=GRID%RL%I3(NC)
            DO N1=1,GRID%RL%NROW
               NI=NI+1
               IF (IDIR==1) II=MOD(N1-INDMIN(1)+NOUT,NOUT)-NOUTH
               IF (IDIR==2) II=MOD(N2-INDMIN(2)+NOUT,NOUT)-NOUTH
               IF (IDIR==3) II=MOD(N3-INDMIN(3)+NOUT,NOUT)-NOUTH
               XX=ABS(ABS(II)-NOUTH)
               IF (XX > WIDTH) THEN
                  CUTOFF=1.0
               ELSE
                  CUTOFF=ABS(SIN(PI*XX/WIDTH/2))
               ENDIF
               POTCOR=-QUADFAC*II*II*CUTOFF*CUTOFF
               GWORK(NI)=GWORK(NI)+POTCOR
            ENDDO
         ENDDO
      ENDDO dir
    END SUBROUTINE QUADRATIC_FIELD_CORRECTION


!==========================================================================
!
!> standard calling interface to the SET_GFAC routine
!> requires one to pass down two k-points for which the Coulomb kernel is required
!> this one is without explicitly compiled interface
!
!==========================================================================

    SUBROUTINE SET_GFAC_NOINTER(GRID, LATT_CUR, NK, NQ, FSG, POTFAK)
      USE fock
      USE lattice
      USE full_kpoints
      IMPLICIT NONE

      TYPE (grid_3d) GRID
      TYPE (latt) LATT_CUR
      INTEGER NK, NQ
      REAL(q) :: FSG
      REAL(q) :: POTFAK(GRID%MPLWV)

      CALL SET_GFAC_VKPT(GRID, LATT_CUR, KPOINTS_FULL%VKPT(:,NK)-KPOINTS_FULL%VKPT(:,NQ), NQ, FSG, POTFAK)
    END SUBROUTINE SET_GFAC_NOINTER
