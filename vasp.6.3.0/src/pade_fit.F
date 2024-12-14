#include "symbol.inc"

!*********************************************************************
!
! Module pade_fit includes some routines for an analytic continuation 
! of diagonal elements of imaginary frequency-dependent Sigma(iw) to 
! the real axis, and calculate the QP energies and spectral functions
!
!*********************************************************************

MODULE pade_fit
   USE prec
   USE wave
   USE vaspxml
#ifdef qd_emulate
   USE qdmodule
#endif
   IMPLICIT NONE
   REAL(q),PARAMETER           :: ZMIN=0.45_q
   REAL(q),PARAMETER           :: ZMAX=1.35_q
   INTEGER,SAVE                :: NBANDS_OUT,KPTS_OUT,ISP_OUT
   INTEGER,ALLOCATABLE,SAVE    :: NBANDS_INDEX(:),KPTS_INDEX(:),ISP_INDEX(:)
   INTEGER,ALLOCATABLE,SAVE    :: TAG(:,:,:)              ! 1-SVD, 2-fractional
   COMPLEX(q),ALLOCATABLE,SAVE :: SIGMA(:,:,:)            ! self-energy at DFT eigenvalue
   COMPLEX(q),ALLOCATABLE,SAVE :: FRAC_COEFF_UP(:,:,:,:)  ! pade coefficients for fractional version
   COMPLEX(q),ALLOCATABLE,SAVE :: FRAC_COEFF_DN(:,:,:,:)  
   COMPLEX(q),ALLOCATABLE,SAVE :: SVD_COEFF_UP(:,:,:,:)   ! pade coefficients for SVD version
   COMPLEX(q),ALLOCATABLE,SAVE :: SVD_COEFF_DN(:,:,:,:)   
   COMPLEX(q),ALLOCATABLE,SAVE :: QP_ZEROS(:,:,:,:)       ! QP_ZEROS(2,NBANDSGW,NKPTS_IRZ,ISPIN)
   REAL(q),ALLOCATABLE,SAVE    :: QP_LINEAR(:,:,:,:)      ! QP_LINEAR(2,NBANDSGW,NKPTS_IRZ,ISPIN)
   REAL(q),ALLOCATABLE,SAVE    :: QP_DIFF(:,:,:)          ! QP_DIFF= QP(Fractional)-QP(SVD)
   LOGICAL,SAVE                :: LPADE_FIT               ! whether you want to show the fitted G and Sigma 

 
CONTAINS


!*********************************************************************
!
! reading related parameters from INCAR for pade fit
!
!*********************************************************************
SUBROUTINE QP_PADE_READER(W,IU5,IU0,IU6,NBANDSGW,NKPTS_IRZ)
  USE reader_tags
  IMPLICIT NONE
  TYPE (wavespin)   :: W
  INTEGER           :: IU5,IU0,IU6
  INTEGER           :: NBANDSGW                   !number of bands to be written to file
  INTEGER           :: NKPTS_IRZ                  !number of KPOINTS in the IBZ
  ! local
  INTEGER           :: I,IDUM,N,IERR
  REAL(q)           :: RDUM
  COMPLEX(q)        :: CDUM
  LOGICAL           :: LOPEN, LDUM
  CHARACTER(1)      :: CHARAC
  CHARACTER(40)     :: STRING
  CHARACTER(LEN=10) :: INCAR='INCAR'

  CALL OPEN_INCAR_IF_FOUND(IU5, LOPEN)

  NBANDS_OUT=NBANDSGW
  CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBANDS_OUT', NBANDS_OUT, IERR, WRITEXMLINCAR)
  IF (IU0>0 .AND. NBANDS_OUT>NBANDSGW) THEN
      WRITE(IU0,*)'Error: NBANDS_OUT should be <= NBANDSGW, We set NBANDS_OUT=NBANDSGW'
      NBANDS_OUT=NBANDSGW
  ENDIF
  IF (ALLOCATED(NBANDS_INDEX)) DEALLOCATE(NBANDS_INDEX)
  ALLOCATE(NBANDS_INDEX(NBANDS_OUT))

  DO I=1,NBANDS_OUT 
     NBANDS_INDEX(I)=I
  ENDDO
  CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'NBANDS_INDEX', NBANDS_INDEX, NBANDS_OUT, IERR, WRITEXMLINCAR)

  KPTS_OUT=NKPTS_IRZ  !by default
  CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'KPTS_OUT', KPTS_OUT, IERR, WRITEXMLINCAR)
  IF (IU0>0 .AND. KPTS_OUT>NKPTS_IRZ) THEN
      WRITE(IU0,*)'Error:KPTS_OUT should be <= NKPTS_IRZ, We set KPTS_OUT=NKPTS_IRZ'
      KPTS_OUT=NKPTS_IRZ
  ENDIF
  IF (ALLOCATED(KPTS_INDEX)) DEALLOCATE(KPTS_INDEX)
  ALLOCATE(KPTS_INDEX(KPTS_OUT))

  DO I=1,KPTS_OUT
     KPTS_INDEX(I)=I
  ENDDO
  CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'KPTS_INDEX', KPTS_INDEX, KPTS_OUT, IERR, WRITEXMLINCAR)

  ISP_OUT=W%WDES%ISPIN  !by default
  CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ISP_OUT', ISP_OUT, IERR, WRITEXMLINCAR)
  IF (IU0>0 .AND. ISP_OUT>W%WDES%ISPIN) THEN
      WRITE(IU0,*)'Error:ISP_OUT should be <= W%WDES%ISPIN, We set ISP_OUT=W%WDES%ISPIN'
      ISP_OUT=W%WDES%ISPIN
  ENDIF
  IF (ALLOCATED(ISP_INDEX)) DEALLOCATE(ISP_INDEX)
  ALLOCATE(ISP_INDEX(ISP_OUT))

  DO I=1,ISP_OUT
     ISP_INDEX(I)=I
  ENDDO
  CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'ISP_INDEX', ISP_INDEX, ISP_OUT, IERR, WRITEXMLINCAR)

  LPADE_FIT=.FALSE.
  CALL PROCESS_INCAR(LOPEN, IU0, IU5, 'LPADE_FIT', LPADE_FIT, IERR, WRITEXMLINCAR)

  CALL CLOSE_INCAR_IF_FOUND(IU5)
END SUBROUTINE QP_PADE_READER


!*********************************************************************
!
! allocate all global variables in this module
!
!*********************************************************************
SUBROUTINE ALLOCATE_ALL_PADE(NOMEGA,NOMEGA_SVD,NBANDSGW,NKPTS_IRZ,ISPIN)
    IMPLICIT NONE
    INTEGER         :: NOMEGA
    INTEGER         :: NOMEGA_SVD
    INTEGER         :: NBANDSGW   
    INTEGER         :: NKPTS_IRZ
    INTEGER         :: ISPIN

    IF ( ALLOCATED(QP_LINEAR) )  DEALLOCATE(QP_LINEAR)
    ALLOCATE( QP_LINEAR(2,NBANDSGW,NKPTS_IRZ,ISPIN) )
    QP_LINEAR=0.0_q

    IF ( ALLOCATED(QP_DIFF) )  DEALLOCATE(QP_DIFF)
    ALLOCATE( QP_DIFF(NBANDSGW,NKPTS_IRZ,ISPIN) )
    QP_DIFF=0.0_q

    IF ( ALLOCATED(QP_ZEROS) )  DEALLOCATE(QP_ZEROS)
    ALLOCATE( QP_ZEROS(2,NBANDSGW,NKPTS_IRZ,ISPIN) )
    QP_ZEROS=(0.0_q,0.0_q)

    IF ( ALLOCATED(SVD_COEFF_UP) )  DEALLOCATE(SVD_COEFF_UP)
    ALLOCATE( SVD_COEFF_UP(NOMEGA_SVD,NBANDSGW,NKPTS_IRZ,ISPIN) )
    SVD_COEFF_UP=(0.0_q,0.0_q)

    IF ( ALLOCATED(FRAC_COEFF_UP) )  DEALLOCATE(FRAC_COEFF_UP)
    ALLOCATE( FRAC_COEFF_UP(NOMEGA,NBANDSGW,NKPTS_IRZ,ISPIN) )
    FRAC_COEFF_UP=(0.0_q,0.0_q)

    IF ( ALLOCATED(TAG) )  DEALLOCATE(TAG)
    ALLOCATE( TAG(NBANDSGW,NKPTS_IRZ,ISPIN) )
    TAG=0

    IF ( ALLOCATED(SIGMA) )  DEALLOCATE(SIGMA)
    ALLOCATE( SIGMA(NBANDSGW,NKPTS_IRZ,ISPIN) )
    SIGMA=(0.0_q,0.0_q)

    IF (LPADE_FIT) THEN
       IF ( ALLOCATED(FRAC_COEFF_DN) )  DEALLOCATE(FRAC_COEFF_DN)
       ALLOCATE( FRAC_COEFF_DN(NOMEGA,NBANDSGW,NKPTS_IRZ,ISPIN) )
       FRAC_COEFF_DN=(0.0_q,0.0_q)
       IF ( ALLOCATED(SVD_COEFF_DN) )  DEALLOCATE(SVD_COEFF_DN)
       ALLOCATE( SVD_COEFF_DN(NOMEGA_SVD,NBANDSGW,NKPTS_IRZ,ISPIN) )
       SVD_COEFF_DN=(0.0_q,0.0_q)
    ENDIF
END SUBROUTINE ALLOCATE_ALL_PADE


!*********************************************************************
!
! deallocate all global variables in this module
!
!*********************************************************************
SUBROUTINE DEALLOCATE_ALL_PADE()
    USE chi_glb, ONLY: LOEP 
    IF ( ALLOCATED(QP_ZEROS) .AND. .NOT. LOEP    )  DEALLOCATE(QP_ZEROS)
    IF ( ALLOCATED(QP_DIFF)      )  DEALLOCATE(QP_DIFF)
    IF ( ALLOCATED(TAG)          )  DEALLOCATE(TAG)
   !IF ( ALLOCATED(SIGMA)        )  DEALLOCATE(SIGMA)    !Keep SIGMA for later use
    IF ( ALLOCATED(FRAC_COEFF_UP))  DEALLOCATE(FRAC_COEFF_UP)
    IF ( ALLOCATED(FRAC_COEFF_DN))  DEALLOCATE(FRAC_COEFF_DN)
    IF ( ALLOCATED(SVD_COEFF_UP) )  DEALLOCATE(SVD_COEFF_UP)
    IF ( ALLOCATED(SVD_COEFF_DN) )  DEALLOCATE(SVD_COEFF_DN)
END SUBROUTINE DEALLOCATE_ALL_PADE


!*********************************************************************
!
! Main routine to do the Pade fit for G0W0R, GW0R, and GWR
! using the imaginary time self-energy
!
!*********************************************************************
SUBROUTINE QP_PADE_FIT_FROM_TAU(W,WMEAN,IMAG_GRIDS,&
    SIGMAO,SIGMAU,NKPTS_IRZ,NBANDSGW,IO,LDUMP)
    USE mpimy
    USE minimax , ONLY :LOCAL_INDEX_CYCLIC 
    USE minimax_struct , ONLY :imag_grid_handle
    USE base , ONLY : in_struct
    IMPLICIT NONE
    TYPE(wavespin)         :: W                    ! wave function 
    TYPE(wavespin)         :: WMEAN                ! identical to W, except for CELTOT and FERTOT
    TYPE(imag_grid_handle) :: IMAG_GRIDS           ! grid handle
    GDEF                   :: SIGMAO(:,:,:,:)      ! diagonal of self-energy in tau domain
    GDEF                   :: SIGMAU(:,:,:,:)      ! diagonal of self-energy in tau domain
    INTEGER                :: NKPTS_IRZ            ! number of KPOINTS in the IBZ
    INTEGER                :: NBANDSGW             ! number of bands to be written to file
    TYPE(in_struct)        :: IO
    LOGICAL, OPTIONAL      :: LDUMP
    ! local variables 
    INTEGER,PARAMETER      :: NP=5             ! NP-point finite difference stencil
    INTEGER                :: IU6                  !for output OUTCAR
    INTEGER                :: IU0
    INTEGER                :: NOMEGA           ! number of (z_i,f_i) points,M
    INTEGER                :: NOMEGA_SVD       ! number of coefficients of SVD <=NOMEGA
    INTEGER                :: IND,NOMEGA_TMP
    INTEGER                :: I,J,K,L
    INTEGER                :: NB,NK,ISP
    REAL(q)                :: RETMP,IMTMP
    REAL(q),ALLOCATABLE    :: QP_SVD(:,:,:,:)     ! SVD version
    REAL(q),ALLOCATABLE    :: QP_FRAC(:,:,:,:)    ! fractional version
    COMPLEX(q),ALLOCATABLE :: IW_UP(:),SELF_IW_UP(:)
    COMPLEX(q),ALLOCATABLE :: IW_DN(:),SELF_IW_DN(:)
    COMPLEX(q)             :: EQP_TMP
    COMPLEX(q),ALLOCATABLE :: SIGMA_SVD(:,:,:)     
    COMPLEX(q),ALLOCATABLE :: SIGMA_FRAC(:,:,:)     
    COMPLEX(q),ALLOCATABLE :: ZEROS_SVD(:,:,:,:)  ! SVD version
    COMPLEX(q),ALLOCATABLE :: ZEROS_FRAC(:,:,:,:) ! fractional version
    INTEGER,ALLOCATABLE    :: I_FRAC(:,:,:,:) 
    ! for fourier transform from time to imaginary frequency
    INTEGER                :: NTAU
    REAL(q),ALLOCATABLE    :: FREQCOS(:)          !cosine imaginary frequecny
    REAL(q),ALLOCATABLE    :: FTCOS(:,:)          !cosine imaginary frequecny
    REAL(q),ALLOCATABLE    :: FTSIN(:,:)          !cosine imaginary frequecny
    REAL(q)                :: EFERMI(W%WDES%ISPIN) ! fermi energy
    INTEGER                :: NODE_ME=1
    INTEGER                :: IONODE=1
    INTEGER                :: NCPU=1
    PROFILING_START('QP_PADE_FIT_FROM_TAU')

    IU6=IO%IU6
    IU0=IO%IU0

    EFERMI = W%EFERMI
    ! include w=0 point
    NOMEGA=IMAG_GRIDS%NOMEGA
    NTAU=IMAG_GRIDS%T%NPOINTS
    ! frequency points for analytical continuation
    ! store frequency points 
    ALLOCATE( FREQCOS( NOMEGA ) ) 
    ALLOCATE( FTCOS( NOMEGA,NTAU ) ) 
    ALLOCATE( FTSIN( NOMEGA,NTAU ) ) 

    FREQCOS(1:NOMEGA) = IMAG_GRIDS%FER_RE(1:NOMEGA)
    ! restore points and matrices 
    DO I = 1, IMAG_GRIDS%NOMEGA
       FREQCOS( I ) = IMAG_GRIDS%FER_RE( I ) 
       FTCOS(I,1:NTAU) = IMAG_GRIDS%TO_FER_RE(I,1:NTAU)
       FTSIN(I,1:NTAU) = IMAG_GRIDS%TO_FER_RE_CONJG(I,1:NTAU)
    ENDDO

    ! SVD version only works for NOMEGA=even
    IF (MOD(NOMEGA,2) /= 0) THEN
      NOMEGA_SVD=NOMEGA-1
    ELSE
      NOMEGA_SVD=NOMEGA
    ENDIF

    ! allocate global variables
    CALL ALLOCATE_ALL_PADE(NOMEGA,NOMEGA_SVD,NBANDSGW,NKPTS_IRZ,W%WDES%ISPIN)

    ! allocate other temperary variables
    ALLOCATE( QP_SVD(2,NBANDSGW,NKPTS_IRZ,W%WDES%ISPIN) )
    QP_SVD=0.0_q
    ALLOCATE( QP_FRAC(2,NBANDSGW,NKPTS_IRZ,W%WDES%ISPIN) )
    QP_FRAC=0.0_q
    ALLOCATE( ZEROS_SVD(2,NBANDSGW,NKPTS_IRZ,W%WDES%ISPIN) )
    ZEROS_SVD=(0.0_q,0.0_q)
    ALLOCATE( ZEROS_FRAC(2,NBANDSGW,NKPTS_IRZ,W%WDES%ISPIN) )
    ZEROS_FRAC=(0.0_q,0.0_q)
    ALLOCATE( SIGMA_SVD(NBANDSGW,NKPTS_IRZ,W%WDES%ISPIN) )
    SIGMA_SVD=(0.0_q,0.0_q)
    ALLOCATE( SIGMA_FRAC(NBANDSGW,NKPTS_IRZ,W%WDES%ISPIN) )
    SIGMA_FRAC=(0.0_q,0.0_q)
    ALLOCATE( I_FRAC(2,NBANDSGW,NKPTS_IRZ,W%WDES%ISPIN) )
    I_FRAC=0

#ifdef MPI
    NODE_ME= W%WDES%COMM%NODE_ME
    IONODE = W%WDES%COMM%IONODE
    NCPU   = W%WDES%COMM%NCPU
#endif
  
    ALLOCATE( IW_UP(NOMEGA) )
    ALLOCATE( SELF_IW_UP(NOMEGA) )
    IW_UP=(0.0_q,0.0_q)
    SELF_IW_UP=(0.0_q,0.0_q)
    IF (LPADE_FIT) THEN
       ALLOCATE( IW_DN(NOMEGA) )
       ALLOCATE( SELF_IW_DN(NOMEGA) )
       IW_DN=(0.0_q,0.0_q)
       SELF_IW_DN=(0.0_q,0.0_q)
    ENDIF

    !=====================================================================!
    ! first, we pre-caculate SVD 
    !=====================================================================!
    DO ISP=1, W%WDES%ISPIN
      DO NK=1,NKPTS_IRZ
         DO NB=1,NBANDSGW
#ifdef MPI
          IND = LOCAL_INDEX_CYCLIC( NCPU, NBANDSGW, NODE_ME, NB ) 
          IF ( IND < 1 ) CYCLE 
#endif
             ! store self-energy to auxillary arrays, this includes zero point            
             DO I=1,NOMEGA
                ! set up self-energy including V_xc part 
                RETMP = 0.0_q
                IMTMP = 0.0_q
                DO K = 1, NTAU
                   RETMP = RETMP + (SIGMAU(K,NB,NK,ISP) - SIGMAO(K,NB,NK,ISP))*&
                           FTCOS(I,K)
                   IMTMP = IMTMP - (SIGMAO(K,NB,NK,ISP) + SIGMAU(K,NB,NK,ISP))*&
                           FTSIN(I,K)
                ENDDO
                IW_UP(I)     =CMPLX(EFERMI(ISP),FREQCOS(I),KIND=q)
                SELF_IW_UP(I)=CMPLX(RETMP,IMTMP,KIND=q)
                IF (LPADE_FIT) THEN
                   IW_DN(I)     =CMPLX(EFERMI(ISP),-FREQCOS(I),KIND=q)
                   SELF_IW_DN(I)=CMPLX(RETMP,-IMTMP,KIND=q)
                ENDIF
             ENDDO

             ! get QP shift by linearization around E_DFT for SVD case
             SVD_COEFF_UP(:,NB,NK,ISP)=(0.0_q,0.0_q)
             IF (MOD(NOMEGA,2)==0) THEN
                CALL PADE_LINEAR_COEFF(NOMEGA,IW_UP(1:NOMEGA),SELF_IW_UP(1:NOMEGA),&
                                  & SVD_COEFF_UP(1:NOMEGA,NB,NK,ISP),IU0)
             ELSE
                CALL PADE_SVD_COEFF(NOMEGA,IW_UP(1:NOMEGA),SELF_IW_UP(1:NOMEGA),&
                                  & NOMEGA_SVD,SVD_COEFF_UP(1:NOMEGA_SVD,NB,NK,ISP),IU0)
             ENDIF

             QP_SVD(:,NB,NK,ISP)=0.0_q
             SIGMA_SVD(NB,NK,ISP)=(0.0_q,0.0_q)
             CALL GET_QP_LINEAR_SVD(NOMEGA_SVD,W%CELTOT(NB,NK,ISP),&
                   & SVD_COEFF_UP(1:NOMEGA_SVD,NB,NK,ISP),SIGMA_SVD(NB,NK,ISP),&
                   & QP_SVD(:,NB,NK,ISP),IU0,NP,WMEAN%CELTOT(NB,NK,ISP))

             ! if fails, get QP shift by QR method
             IF ( MOD(NOMEGA,2)==0 .AND. (QP_SVD(2,NB,NK,ISP)>ZMAX .OR. QP_SVD(2,NB,NK,ISP)<ZMIN)) THEN
                SVD_COEFF_UP(:,NB,NK,ISP)=(0.0_q,0.0_q)
                CALL PADE_SVD_COEFF(NOMEGA,IW_UP(1:NOMEGA),SELF_IW_UP(1:NOMEGA),&
                                  & NOMEGA_SVD,SVD_COEFF_UP(1:NOMEGA_SVD,NB,NK,ISP),IU0)
                QP_SVD(:,NB,NK,ISP)=0.0_q
                SIGMA_SVD(NB,NK,ISP)=(0.0_q,0.0_q)
                CALL GET_QP_LINEAR_SVD(NOMEGA,W%CELTOT(NB,NK,ISP),&
                      & SVD_COEFF_UP(1:NOMEGA,NB,NK,ISP),SIGMA_SVD(NB,NK,ISP),&
                      & QP_SVD(:,NB,NK,ISP),IU0,NP,WMEAN%CELTOT(NB,NK,ISP))
             ENDIF

             ! get QP shift by finding the roots of QP equation using Brent-bisection method
             EQP_TMP=CMPLX(QP_SVD(1,NB,NK,ISP),0.0_q,KIND=q)
             ZEROS_SVD(:,NB,NK,ISP)=(0.0_q,0.0_q)
             CALL GET_ZEROS_BRENT_SVD(NOMEGA_SVD,ISP,NK,NB,EQP_TMP,IW_UP(1:NOMEGA_SVD),&
                      & SVD_COEFF_UP(1:NOMEGA_SVD,NB,NK,ISP),ZEROS_SVD(:,NB,NK,ISP),IU0,WMEAN%CELTOT(NB,NK,ISP))

             ! note that for the plot, we compare QR version and FRAC version
             ! LINEAR version almost yields the same as FRAC version
             IF (LPADE_FIT) THEN
                ! negative imaginary axis to negative real axis to avoid the branch cut
                SVD_COEFF_DN(:,NB,NK,ISP)=(0.0_q,0.0_q)
                CALL PADE_SVD_COEFF(NOMEGA,IW_DN(1:NOMEGA),SELF_IW_DN(1:NOMEGA),&
                                  & NOMEGA_SVD,SVD_COEFF_DN(1:NOMEGA_SVD,NB,NK,ISP),IU0)
                SVD_COEFF_UP(:,NB,NK,ISP)=(0.0_q,0.0_q)
                CALL PADE_SVD_COEFF(NOMEGA,IW_UP(1:NOMEGA),SELF_IW_UP(1:NOMEGA),&
                                  & NOMEGA_SVD,SVD_COEFF_UP(1:NOMEGA_SVD,NB,NK,ISP),IU0)
             ENDIF

         ENDDO !NB
      ENDDO !NK
    ENDDO !ISP

    IF (LPADE_FIT) THEN
       CALLMPI( M_sum_z(W%WDES%COMM,SVD_COEFF_DN,SIZE(SVD_COEFF_DN)) )
    ENDIF 
    CALLMPI( M_sum_z(W%WDES%COMM,SVD_COEFF_UP,SIZE(SVD_COEFF_UP)) )
    CALLMPI( M_sum_z(W%WDES%COMM,SIGMA_SVD,SIZE(SIGMA_SVD)) )
    CALLMPI( M_sum_z(W%WDES%COMM,ZEROS_SVD,SIZE(ZEROS_SVD)) )
    CALLMPI( M_sum_d(W%WDES%COMM,QP_SVD,SIZE(QP_SVD)) )


    !=====================================================================!
    ! second, we pre-caculate NOMEGA fractional version
    !=====================================================================!
    DO ISP=1, W%WDES%ISPIN
      DO NK=1,NKPTS_IRZ
         DO NB=1,NBANDSGW
#ifdef MPI
          IND = LOCAL_INDEX_CYCLIC( NCPU, NBANDSGW, NODE_ME, NB )
          IF ( IND < 1 ) CYCLE
#endif
             ! store self-energy to auxillary arrays, this includes zero point            
             DO I=1,NOMEGA
                ! set up self-energy including V_xc part 
                RETMP = 0.0_q
                IMTMP = 0.0_q
                DO K = 1, NTAU
                   RETMP = RETMP + (SIGMAU(K,NB,NK,ISP) - SIGMAO(K,NB,NK,ISP))*&
                           FTCOS(I,K)
                   IMTMP = IMTMP - (SIGMAO(K,NB,NK,ISP) + SIGMAU(K,NB,NK,ISP))*&
                           FTSIN(I,K)
                ENDDO
                IW_UP(I)     =CMPLX(EFERMI(ISP),FREQCOS(I),KIND=q)
                SELF_IW_UP(I)=CMPLX(RETMP,IMTMP,KIND=q)
                IF (LPADE_FIT) THEN
                   IW_DN(I)     =CMPLX(EFERMI(ISP),-FREQCOS(I),KIND=q)
                   SELF_IW_DN(I)=CMPLX(RETMP,-IMTMP,KIND=q)
                ENDIF
             ENDDO

             ! get QP shift by linearization around E_DFT for NOMEGA case
             NOMEGA_TMP=NOMEGA
             CALL CALCULATE_QP(NOMEGA_TMP,QP_FRAC,SIGMA_FRAC)

             ! if fails, get QP shift by linearization around E_DFT for NOMEGA-1 case
             IF ( QP_FRAC(2,NB,NK,ISP)>ZMAX .OR. QP_FRAC(2,NB,NK,ISP)<ZMIN) THEN
                I_FRAC(1,NB,NK,ISP)=1
                NOMEGA_TMP=NOMEGA-1
                CALL CALCULATE_QP(NOMEGA_TMP,QP_FRAC,SIGMA_FRAC)
             ENDIF

             ! get QP shift by finding the roots of QP equation using Brent-bisection method
             EQP_TMP=CMPLX(QP_FRAC(1,NB,NK,ISP),0.0_q,KIND=q)

             NOMEGA_TMP=NOMEGA
             DO WHILE(.TRUE.)
                FRAC_COEFF_UP(:,NB,NK,ISP)=(0.0_q,0.0_q)
                CALL PADE_FRAC_COEFF(NOMEGA_TMP,IW_UP(NOMEGA-NOMEGA_TMP+1:NOMEGA),&
                     & SELF_IW_UP(NOMEGA-NOMEGA_TMP+1:NOMEGA),FRAC_COEFF_UP(1:NOMEGA_TMP,NB,NK,ISP),IU0)
                ZEROS_FRAC(:,NB,NK,ISP)=(0.0_q,0.0_q)
                CALL GET_ZEROS_BRENT_FRAC(NOMEGA_TMP,ISP,NK,NB,EQP_TMP,&
                         & IW_UP(NOMEGA-NOMEGA_TMP+1:NOMEGA),FRAC_COEFF_UP(1:NOMEGA_TMP,NB,NK,ISP),&
                         & ZEROS_FRAC(:,NB,NK,ISP),IU0,WMEAN%CELTOT(NB,NK,ISP))
                IF (REAL(ZEROS_FRAC(2,NB,NK,ISP),q)>ZMAX .OR. REAL(ZEROS_FRAC(2,NB,NK,ISP),q)<ZMIN ) THEN
                   NOMEGA_TMP=NOMEGA_TMP-1
                   I_FRAC(2,NB,NK,ISP)=1
                ELSE
                   EXIT
                ENDIF
                IF (NOMEGA_TMP<MAX(NOMEGA-2,8) ) THEN
                   IF (IU0>=0) WRITE(IU0,'(A,3I6)') &
                       & "Warning: NOMEGA_TMP<MAX(NOMEGA-2,8) in GET_ZEROS_BRENT_FRAC",NB,NK,ISP
                   EXIT
                ENDIF
             ENDDO

             ! if still problematic, we report error and use original NOMEGA one 
             IF (NOMEGA_TMP<MAX(NOMEGA-2,8) .OR. & 
                & REAL(ZEROS_FRAC(2,NB,NK,ISP),q)>ZMAX .OR. REAL(ZEROS_FRAC(2,NB,NK,ISP),q)<ZMIN ) THEN
                IF (IU0>=0) WRITE(IU0,'(A,3I6)') "Warning: GET_ZEROS_BRENT_FRAC in Z ", NB,NK,ISP
                FRAC_COEFF_UP(:,NB,NK,ISP)=(0.0_q,0.0_q)
                CALL PADE_FRAC_COEFF(NOMEGA,IW_UP(1:NOMEGA),SELF_IW_UP(1:NOMEGA),&
                               & FRAC_COEFF_UP(1:NOMEGA,NB,NK,ISP),IU0)
                ZEROS_FRAC(:,NB,NK,ISP)=(0.0_q,0.0_q)
                CALL GET_ZEROS_BRENT_FRAC(NOMEGA,ISP,NK,NB,EQP_TMP,IW_UP(1:NOMEGA),&
                            & FRAC_COEFF_UP(1:NOMEGA,NB,NK,ISP),ZEROS_FRAC(:,NB,NK,ISP),IU0,WMEAN%CELTOT(NB,NK,ISP))
             ENDIF

             ! note that for the plot, we compare QR version and FRAC version
             ! LINEAR version almost yields the same as FRAC version
             IF (LPADE_FIT) THEN
                ! negative imaginary axis to negative real axis to avoid the branch cut
                FRAC_COEFF_DN(:,NB,NK,ISP)=(0.0_q,0.0_q)
                CALL PADE_FRAC_COEFF(NOMEGA,IW_DN(1:NOMEGA),SELF_IW_DN(1:NOMEGA),&
                         & FRAC_COEFF_DN(1:NOMEGA,NB,NK,ISP),IU0)
                FRAC_COEFF_UP(:,NB,NK,ISP)=(0.0_q,0.0_q)
                CALL PADE_FRAC_COEFF(NOMEGA,IW_UP(1:NOMEGA),SELF_IW_UP(1:NOMEGA),&
                         & FRAC_COEFF_UP(1:NOMEGA,NB,NK,ISP),IU0)
             ENDIF


         ENDDO !NB
      ENDDO !NK
    ENDDO !ISP

    IF (LPADE_FIT) THEN
       CALLMPI( M_sum_z(W%WDES%COMM,FRAC_COEFF_DN,SIZE(FRAC_COEFF_DN)) )
    ENDIF
    CALLMPI( M_sum_z(W%WDES%COMM,FRAC_COEFF_UP,SIZE(FRAC_COEFF_UP)) )
    CALLMPI( M_sum_z(W%WDES%COMM,SIGMA_FRAC,SIZE(SIGMA_FRAC)) )
    CALLMPI( M_sum_z(W%WDES%COMM,ZEROS_FRAC,SIZE(ZEROS_FRAC)) )
    CALLMPI( M_sum_d(W%WDES%COMM,QP_FRAC,SIZE(QP_FRAC)) )
    CALLMPI( M_sum_i(W%WDES%COMM,I_FRAC,SIZE(I_FRAC)) )


    !=====================================================================!
    ! third, we compare SVD and fractional version to select the 'best' one
    !=====================================================================!
    DO ISP=1, W%WDES%ISPIN  
      DO NK=1,NKPTS_IRZ
         DO NB=1,NBANDSGW
#ifdef MPI
          IND = LOCAL_INDEX_CYCLIC( NCPU, NBANDSGW, NODE_ME, NB )
          IF ( IND < 1 ) CYCLE
#endif

             ! get QP shift by linearization around E_DFT
             QP_LINEAR(:,NB,NK,ISP)=0.0_q

             QP_DIFF(NB,NK,ISP)=QP_FRAC(1,NB,NK,ISP)-QP_SVD(1,NB,NK,ISP)

             ! the criteria is that we try to keep Fractional version
             ! otherwise, we use the SVD version
             ! for linerization 
             IF (QP_FRAC(2,NB,NK,ISP)>=0.60_q .AND. &
               & QP_FRAC(2,NB,NK,ISP)<=1.0_q .AND. I_FRAC(1,NB,NK,ISP)==0) THEN
                QP_LINEAR(:,NB,NK,ISP)=QP_FRAC(:,NB,NK,ISP)
                SIGMA(NB,NK,ISP)=SIGMA_FRAC(NB,NK,ISP)
                TAG(NB,NK,ISP)=2
             ELSE IF (QP_SVD(2,NB,NK,ISP)>=0.60_q .AND. QP_SVD(2,NB,NK,ISP)<=1.0_q) THEN
                QP_LINEAR(:,NB,NK,ISP)=QP_SVD(:,NB,NK,ISP)
                SIGMA(NB,NK,ISP)=SIGMA_SVD(NB,NK,ISP)
                TAG(NB,NK,ISP)=1
             ELSE
                QP_LINEAR(:,NB,NK,ISP)=QP_FRAC(:,NB,NK,ISP)
                SIGMA(NB,NK,ISP)=SIGMA_FRAC(NB,NK,ISP)
                TAG(NB,NK,ISP)=2
             ENDIF

             ! for Zeros
             IF ( REAL(ZEROS_FRAC(2,NB,NK,ISP),q)>=0.60_q .AND. &
                & REAL(ZEROS_FRAC(2,NB,NK,ISP),q)<=1.0_q .AND. I_FRAC(2,NB,NK,ISP)==0) THEN
                QP_ZEROS(:,NB,NK,ISP)=ZEROS_FRAC(:,NB,NK,ISP)
             ELSE IF (REAL(ZEROS_SVD(2,NB,NK,ISP),q)>=0.60_q .AND. REAL(ZEROS_SVD(2,NB,NK,ISP),q)<=1.0_q) THEN
                QP_ZEROS(:,NB,NK,ISP)=ZEROS_SVD(:,NB,NK,ISP)
             ELSE
                QP_ZEROS(:,NB,NK,ISP)=ZEROS_FRAC(:,NB,NK,ISP)
             ENDIF

         ENDDO !NB
      ENDDO !NK
    ENDDO !ISP

    CALLMPI( M_sum_z(W%WDES%COMM,SIGMA,SIZE(SIGMA)) )
    CALLMPI( M_sum_d(W%WDES%COMM,QP_LINEAR,SIZE(QP_LINEAR)) )
    CALLMPI( M_sum_d(W%WDES%COMM,QP_DIFF,SIZE(QP_DIFF)) )
    CALLMPI( M_sum_z(W%WDES%COMM,QP_ZEROS,SIZE(QP_ZEROS)) )
    CALLMPI( M_sum_i(W%WDES%COMM,TAG,SIZE(TAG)) )

    IF (NODE_ME==1) THEN
      ! write out the QP_shift to OUTCAR
      CALL WRITE_QP_SHIFT(W,NKPTS_IRZ,NBANDSGW,TAG,IU6)
      ! write out the Pade fitted sigma(w) for specified bands and kpoints in vasprun.xml file
      IF (LPADE_FIT) THEN
         DO ISP=1,ISP_OUT
           CALL WRITE_SIGMA_XML(WMEAN%CELTOT,NOMEGA,NOMEGA_SVD,IW_UP,IW_DN,&
                  & EFERMI(ISP_INDEX(ISP)),2001,ISP_INDEX(ISP),IU0)
         ENDDO
      ENDIF
    ENDIF  

    CALL DEALLOCATE_ALL_PADE()
    DEALLOCATE(QP_SVD, ZEROS_SVD, SIGMA_SVD)
    DEALLOCATE(QP_FRAC, ZEROS_FRAC, SIGMA_FRAC, I_FRAC)
    DEALLOCATE(IW_UP,SELF_IW_UP)
    IF (LPADE_FIT) DEALLOCATE(IW_DN,SELF_IW_DN)
    DEALLOCATE(FREQCOS, FTCOS, FTSIN)

    PROFILING_STOP('QP_PADE_FIT_FROM_TAU')

CONTAINS


!*********************************************************************
!
! get QP shift by linearization around E_DFT for NOMEGA_TMP case
!
!*********************************************************************
SUBROUTINE CALCULATE_QP(NOMEGA_TMP,QP_BB,SIGMA_BB)
    IMPLICIT NONE
    INTEGER     :: NOMEGA_TMP
    REAL(q)     :: QP_BB(:,:,:,:)
    COMPLEX(q)  :: SIGMA_BB(:,:,:)

    FRAC_COEFF_UP(:,NB,NK,ISP)=(0.0_q,0.0_q)
    CALL PADE_FRAC_COEFF(NOMEGA_TMP,IW_UP(NOMEGA-NOMEGA_TMP+1:NOMEGA),&
         & SELF_IW_UP(NOMEGA-NOMEGA_TMP+1:NOMEGA),FRAC_COEFF_UP(1:NOMEGA_TMP,NB,NK,ISP),IU0)
    QP_BB(:,NB,NK,ISP)=0.0_q
    CALL GET_QP_LINEAR_FRAC(NOMEGA_TMP,W%CELTOT(NB,NK,ISP),&
          & IW_UP(NOMEGA-NOMEGA_TMP+1:NOMEGA),FRAC_COEFF_UP(1:NOMEGA_TMP,NB,NK,ISP),&
          & SIGMA_BB(NB,NK,ISP),QP_BB(:,NB,NK,ISP),IU0,NP,WMEAN%CELTOT(NB,NK,ISP))

END SUBROUTINE CALCULATE_QP

END SUBROUTINE QP_PADE_FIT_FROM_TAU



!*********************************************************************
!
! Same above but using the imaginary frequency self-energy 
!
!*********************************************************************
SUBROUTINE QP_PADE_FIT(W,WMEAN,FREQCOS,&
          & SIGMA_COS,SIGMA_SIN,NKPTS_IRZ,NBANDSGW,IU6,IU0)
    USE mpimy
    USE minimax, ONLY: LOCAL_INDEX_CYCLIC
    IMPLICIT NONE
    TYPE(wavespin)  :: W                       ! wave function 
    TYPE(wavespin)  :: WMEAN                   ! identical to W, except for CELTOT and FERTOT
    REAL(q)         :: FREQCOS(:)              ! cosine imaginary frequecny
    COMPLEX(q)      :: SIGMA_COS(:,:,:,:)      ! diagonal of cos transformed Self-energy
    COMPLEX(q)      :: SIGMA_SIN(:,:,:,:)      ! diagonal of sin transformed Self-energy
    INTEGER         :: NKPTS_IRZ               ! number of KPOINTS in the IBZ
    INTEGER         :: NBANDSGW                ! number of bands to be written to file
    INTEGER         :: IU6                     ! for output OUTCAR
    INTEGER         :: IU0
    ! local variables 
    INTEGER,PARAMETER      :: NP=5             ! NP-point finite difference stencil
    REAL(q)         :: EFERMI(W%WDES%ISPIN)    ! fermi energy
    INTEGER                :: NOMEGA           ! number of (z_i,f_i) points,M
    INTEGER                :: NOMEGA_SVD       ! number of coefficients of SVD <=NOMEGA
    INTEGER                :: IND,NOMEGA_TMP
    INTEGER                :: I,J,K,L
    INTEGER                :: NB,NK,ISP
    REAL(q)                :: RETMP,IMTMP
    REAL(q),ALLOCATABLE    :: QP_SVD(:,:,:,:)     ! SVD version
    REAL(q),ALLOCATABLE    :: QP_FRAC(:,:,:,:)    ! fractional version
    COMPLEX(q),ALLOCATABLE :: IW_UP(:),SELF_IW_UP(:)
    COMPLEX(q),ALLOCATABLE :: IW_DN(:),SELF_IW_DN(:)
    COMPLEX(q)             :: EQP_TMP
    COMPLEX(q),ALLOCATABLE :: SIGMA_SVD(:,:,:)     
    COMPLEX(q),ALLOCATABLE :: SIGMA_FRAC(:,:,:)     
    COMPLEX(q),ALLOCATABLE :: ZEROS_SVD(:,:,:,:)  ! SVD version
    COMPLEX(q),ALLOCATABLE :: ZEROS_FRAC(:,:,:,:) ! fractional version
    INTEGER,ALLOCATABLE    :: I_FRAC(:,:,:,:) 
    INTEGER                :: NODE_ME =1
    INTEGER                :: IONODE=1
    INTEGER                :: NCPU=1

    NOMEGA=SIZE(SIGMA_COS,1)
    EFERMI = W%EFERMI

    ! SVD version only works for NOMEGA=even
    IF (MOD(NOMEGA,2) /= 0) THEN
      NOMEGA_SVD=NOMEGA-1
    ELSE
      NOMEGA_SVD=NOMEGA
    ENDIF

    ! allocate global variables
    CALL ALLOCATE_ALL_PADE(NOMEGA,NOMEGA_SVD,NBANDSGW,NKPTS_IRZ,W%WDES%ISPIN)

    ! allocate other temperary variables
    ALLOCATE( QP_SVD(2,NBANDSGW,NKPTS_IRZ,W%WDES%ISPIN) )
    QP_SVD=0.0_q
    ALLOCATE( QP_FRAC(2,NBANDSGW,NKPTS_IRZ,W%WDES%ISPIN) )
    QP_FRAC=0.0_q
    ALLOCATE( ZEROS_SVD(2,NBANDSGW,NKPTS_IRZ,W%WDES%ISPIN) )
    ZEROS_SVD=(0.0_q,0.0_q)
    ALLOCATE( ZEROS_FRAC(2,NBANDSGW,NKPTS_IRZ,W%WDES%ISPIN) )
    ZEROS_FRAC=(0.0_q,0.0_q)
    ALLOCATE( SIGMA_SVD(NBANDSGW,NKPTS_IRZ,W%WDES%ISPIN) )
    SIGMA_SVD=(0.0_q,0.0_q)
    ALLOCATE( SIGMA_FRAC(NBANDSGW,NKPTS_IRZ,W%WDES%ISPIN) )
    SIGMA_FRAC=(0.0_q,0.0_q)
    ALLOCATE( I_FRAC(2,NBANDSGW,NKPTS_IRZ,W%WDES%ISPIN) )
    I_FRAC=0

#ifdef MPI
    NODE_ME= W%WDES%COMM%NODE_ME
    IONODE = W%WDES%COMM%IONODE
    NCPU   = W%WDES%COMM%NCPU
#else
    NODE_ME= 1
    IONODE = 1
    NCPU   = 1
#endif
  
    ALLOCATE( IW_UP(NOMEGA) )
    ALLOCATE( SELF_IW_UP(NOMEGA) )
    IW_UP=(0.0_q,0.0_q)
    SELF_IW_UP=(0.0_q,0.0_q)
    IF (LPADE_FIT) THEN
       ALLOCATE( IW_DN(NOMEGA) )
       ALLOCATE( SELF_IW_DN(NOMEGA) )
       IW_DN=(0.0_q,0.0_q)
       SELF_IW_DN=(0.0_q,0.0_q)
    ENDIF

    !=====================================================================!
    ! first, we pre-caculate SVD 
    !=====================================================================!
    DO ISP=1, W%WDES%ISPIN
      DO NK=1,NKPTS_IRZ
         DO NB=1,NBANDSGW
#ifdef MPI
          IND = LOCAL_INDEX_CYCLIC( NCPU, NBANDSGW, NODE_ME, NB )
          IF ( IND < 1 ) CYCLE
#endif
             ! store self-energy to auxillary arrays            
             DO I=1,NOMEGA
                RETMP= SIGMA_COS(I,NB,NK,ISP)
                ! "-1" sign comes from G_unocc(paper)=-G_unocc(VASP)
                IMTMP=-SIGMA_SIN(I,NB,NK,ISP)
                IW_UP(I)     =CMPLX(EFERMI(ISP),FREQCOS(I),KIND=q)
                SELF_IW_UP(I)=CMPLX(RETMP,IMTMP,KIND=q)
                IF (LPADE_FIT) THEN
                   IW_DN(I)     =CMPLX(EFERMI(ISP),-FREQCOS(I),KIND=q)
                   SELF_IW_DN(I)=CMPLX(RETMP,-IMTMP,KIND=q)
                ENDIF
             ENDDO

             ! get QP shift by linearization around E_DFT for SVD case
             SVD_COEFF_UP(:,NB,NK,ISP)=(0.0_q,0.0_q)
             IF (MOD(NOMEGA,2)==0) THEN
                CALL PADE_LINEAR_COEFF(NOMEGA,IW_UP(1:NOMEGA),SELF_IW_UP(1:NOMEGA),&
                                  & SVD_COEFF_UP(1:NOMEGA,NB,NK,ISP),IU0)
             ELSE
                CALL PADE_SVD_COEFF(NOMEGA,IW_UP(1:NOMEGA),SELF_IW_UP(1:NOMEGA),&
                                  & NOMEGA_SVD,SVD_COEFF_UP(1:NOMEGA_SVD,NB,NK,ISP),IU0)
             ENDIF

             QP_SVD(:,NB,NK,ISP)=0.0_q
             SIGMA_SVD(NB,NK,ISP)=(0.0_q,0.0_q)
             CALL GET_QP_LINEAR_SVD(NOMEGA_SVD,W%CELTOT(NB,NK,ISP),&
                   & SVD_COEFF_UP(1:NOMEGA_SVD,NB,NK,ISP),SIGMA_SVD(NB,NK,ISP),&
                   & QP_SVD(:,NB,NK,ISP),IU0,NP,WMEAN%CELTOT(NB,NK,ISP))

             ! if fails, get QP shift by QR method
             IF ( MOD(NOMEGA,2)==0 .AND. (QP_SVD(2,NB,NK,ISP)>ZMAX .OR. QP_SVD(2,NB,NK,ISP)<ZMIN)) THEN
                SVD_COEFF_UP(:,NB,NK,ISP)=(0.0_q,0.0_q)
                CALL PADE_SVD_COEFF(NOMEGA,IW_UP(1:NOMEGA),SELF_IW_UP(1:NOMEGA),&
                                  & NOMEGA_SVD,SVD_COEFF_UP(1:NOMEGA_SVD,NB,NK,ISP),IU0)
                QP_SVD(:,NB,NK,ISP)=0.0_q
                SIGMA_SVD(NB,NK,ISP)=(0.0_q,0.0_q)
                CALL GET_QP_LINEAR_SVD(NOMEGA,W%CELTOT(NB,NK,ISP),&
                      & SVD_COEFF_UP(1:NOMEGA,NB,NK,ISP),SIGMA_SVD(NB,NK,ISP),&
                      & QP_SVD(:,NB,NK,ISP),IU0,NP,WMEAN%CELTOT(NB,NK,ISP))
             ENDIF

             ! get QP shift by finding the roots of QP equation using Brent-bisection method
             EQP_TMP=CMPLX(QP_SVD(1,NB,NK,ISP),0.0_q,KIND=q)
             ZEROS_SVD(:,NB,NK,ISP)=(0.0_q,0.0_q)
             CALL GET_ZEROS_BRENT_SVD(NOMEGA_SVD,ISP,NK,NB,EQP_TMP,IW_UP(1:NOMEGA_SVD),&
                      & SVD_COEFF_UP(1:NOMEGA_SVD,NB,NK,ISP),ZEROS_SVD(:,NB,NK,ISP),IU0,WMEAN%CELTOT(NB,NK,ISP))

             ! note that for the plot, we compare QR version and FRAC version
             ! LINEAR version almost yields the same as FRAC version
             IF (LPADE_FIT) THEN
                ! negative imaginary axis to negative real axis to avoid the branch cut
                SVD_COEFF_DN(:,NB,NK,ISP)=(0.0_q,0.0_q)
                CALL PADE_SVD_COEFF(NOMEGA,IW_DN(1:NOMEGA),SELF_IW_DN(1:NOMEGA),&
                                  & NOMEGA_SVD,SVD_COEFF_DN(1:NOMEGA_SVD,NB,NK,ISP),IU0)
                SVD_COEFF_UP(:,NB,NK,ISP)=(0.0_q,0.0_q)
                CALL PADE_SVD_COEFF(NOMEGA,IW_UP(1:NOMEGA),SELF_IW_UP(1:NOMEGA),&
                                  & NOMEGA_SVD,SVD_COEFF_UP(1:NOMEGA_SVD,NB,NK,ISP),IU0)
             ENDIF

         ENDDO !NB
      ENDDO !NK
    ENDDO !ISP

    IF (LPADE_FIT) THEN
       CALLMPI( M_sum_z(W%WDES%COMM,SVD_COEFF_DN,SIZE(SVD_COEFF_DN)) )
    ENDIF
    CALLMPI( M_sum_z(W%WDES%COMM,SVD_COEFF_UP,SIZE(SVD_COEFF_UP)) )
    CALLMPI( M_sum_z(W%WDES%COMM,SIGMA_SVD,SIZE(SIGMA_SVD)) )
    CALLMPI( M_sum_z(W%WDES%COMM,ZEROS_SVD,SIZE(ZEROS_SVD)) )
    CALLMPI( M_sum_d(W%WDES%COMM,QP_SVD,SIZE(QP_SVD)) )


    !=====================================================================!
    !second, we pre-caculate NOMEGA fractional version
    !=====================================================================!
    DO ISP=1, W%WDES%ISPIN
      DO NK=1,NKPTS_IRZ
         DO NB=1,NBANDSGW
#ifdef MPI
          IND = LOCAL_INDEX_CYCLIC( NCPU, NBANDSGW, NODE_ME, NB )
          IF ( IND < 1 ) CYCLE
#endif
             ! store self-energy to auxillary arrays            
             DO I=1,NOMEGA
                RETMP= SIGMA_COS(I,NB,NK,ISP)
                ! "-1" sign comes from G_unocc(paper)=-G_unocc(VASP)
                IMTMP=-SIGMA_SIN(I,NB,NK,ISP)
                IW_UP(I)     =CMPLX(EFERMI(ISP),FREQCOS(I),KIND=q)
                SELF_IW_UP(I)=CMPLX(RETMP,IMTMP,KIND=q)
                IF (LPADE_FIT) THEN
                   IW_DN(I)     =CMPLX(EFERMI(ISP),-FREQCOS(I),KIND=q)
                   SELF_IW_DN(I)=CMPLX(RETMP,-IMTMP,KIND=q)
                ENDIF
             ENDDO

             ! get QP shift by linearization around E_DFT for NOMEGA case
             NOMEGA_TMP=NOMEGA
             CALL CALCULATE_QP(NOMEGA_TMP,QP_FRAC,SIGMA_FRAC)

             ! if fails, get QP shift by linearization around E_DFT for NOMEGA-1 case
             IF ( QP_FRAC(2,NB,NK,ISP)>ZMAX .OR. QP_FRAC(2,NB,NK,ISP)<ZMIN) THEN
                I_FRAC(1,NB,NK,ISP)=1
                NOMEGA_TMP=NOMEGA-1
                CALL CALCULATE_QP(NOMEGA_TMP,QP_FRAC,SIGMA_FRAC)
             ENDIF

             ! get QP shift by finding the roots of QP equation using Brent-bisection method
             EQP_TMP=CMPLX(QP_FRAC(1,NB,NK,ISP),0.0_q,KIND=q)

             NOMEGA_TMP=NOMEGA
             DO WHILE(.TRUE.)
                FRAC_COEFF_UP(:,NB,NK,ISP)=(0.0_q,0.0_q)
                CALL PADE_FRAC_COEFF(NOMEGA_TMP,IW_UP(NOMEGA-NOMEGA_TMP+1:NOMEGA),&
                     & SELF_IW_UP(NOMEGA-NOMEGA_TMP+1:NOMEGA),FRAC_COEFF_UP(1:NOMEGA_TMP,NB,NK,ISP),IU0)
                ZEROS_FRAC(:,NB,NK,ISP)=(0.0_q,0.0_q)
                CALL GET_ZEROS_BRENT_FRAC(NOMEGA_TMP,ISP,NK,NB,EQP_TMP,&
                         & IW_UP(NOMEGA-NOMEGA_TMP+1:NOMEGA),FRAC_COEFF_UP(1:NOMEGA_TMP,NB,NK,ISP),&
                         & ZEROS_FRAC(:,NB,NK,ISP),IU0,WMEAN%CELTOT(NB,NK,ISP))
                IF (REAL(ZEROS_FRAC(2,NB,NK,ISP),q)>ZMAX .OR. REAL(ZEROS_FRAC(2,NB,NK,ISP),q)<ZMIN ) THEN
                   NOMEGA_TMP=NOMEGA_TMP-1
                   I_FRAC(2,NB,NK,ISP)=1
                ELSE
                   EXIT
                ENDIF
                IF (NOMEGA_TMP<MAX(NOMEGA-2,8) ) THEN
                   IF (IU0>=0) WRITE(IU0,'(A,3I6)') &
                       & "Warning: NOMEGA_TMP<MAX(NOMEGA-2,8) in GET_ZEROS_BRENT_FRAC",NB,NK,ISP
                   EXIT
                ENDIF
             ENDDO

             ! if still problematic, we report error and use original NOMEGA one 
             IF (NOMEGA_TMP<MAX(NOMEGA-2,8) .OR. & 
                & REAL(ZEROS_FRAC(2,NB,NK,ISP),q)>ZMAX .OR. REAL(ZEROS_FRAC(2,NB,NK,ISP),q)<ZMIN ) THEN
                IF (IU0>=0) WRITE(IU0,'(A,3I6)') "Warning: GET_ZEROS_BRENT_FRAC in Z ", NB,NK,ISP
                FRAC_COEFF_UP(:,NB,NK,ISP)=(0.0_q,0.0_q)
                CALL PADE_FRAC_COEFF(NOMEGA,IW_UP(1:NOMEGA),SELF_IW_UP(1:NOMEGA),&
                               & FRAC_COEFF_UP(1:NOMEGA,NB,NK,ISP),IU0)
                ZEROS_FRAC(:,NB,NK,ISP)=(0.0_q,0.0_q)
                CALL GET_ZEROS_BRENT_FRAC(NOMEGA,ISP,NK,NB,EQP_TMP,IW_UP(1:NOMEGA),&
                            & FRAC_COEFF_UP(1:NOMEGA,NB,NK,ISP),ZEROS_FRAC(:,NB,NK,ISP),IU0,WMEAN%CELTOT(NB,NK,ISP))
             ENDIF

             ! note that for the plot, we compare QR version and FRAC version
             ! LINEAR version almost yields the same as FRAC version
             IF (LPADE_FIT) THEN
                !negative imaginary axis to negative real axis to avoid the branch cut
                FRAC_COEFF_DN(:,NB,NK,ISP)=(0.0_q,0.0_q)
                CALL PADE_FRAC_COEFF(NOMEGA,IW_DN(1:NOMEGA),SELF_IW_DN(1:NOMEGA),&
                         & FRAC_COEFF_DN(1:NOMEGA,NB,NK,ISP),IU0)
                FRAC_COEFF_UP(:,NB,NK,ISP)=(0.0_q,0.0_q)
                CALL PADE_FRAC_COEFF(NOMEGA,IW_UP(1:NOMEGA),SELF_IW_UP(1:NOMEGA),&
                         & FRAC_COEFF_UP(1:NOMEGA,NB,NK,ISP),IU0)
             ENDIF


         ENDDO !NB
      ENDDO !NK
    ENDDO !ISP

    IF (LPADE_FIT) THEN
       CALLMPI( M_sum_z(W%WDES%COMM,FRAC_COEFF_DN,SIZE(FRAC_COEFF_DN)))
    ENDIF
    CALLMPI( M_sum_z(W%WDES%COMM,FRAC_COEFF_UP,SIZE(FRAC_COEFF_UP)) )
    CALLMPI( M_sum_z(W%WDES%COMM,SIGMA_FRAC,SIZE(SIGMA_FRAC)) )
    CALLMPI( M_sum_z(W%WDES%COMM,ZEROS_FRAC,SIZE(ZEROS_FRAC)) )
    CALLMPI( M_sum_d(W%WDES%COMM,QP_FRAC,SIZE(QP_FRAC)) )
    CALLMPI( M_sum_i(W%WDES%COMM,I_FRAC,SIZE(I_FRAC)) )


    !=====================================================================!
    ! third, we compare SVD and fractional version to select the 'best' one
    !=====================================================================!
    DO ISP=1, W%WDES%ISPIN  
      DO NK=1,NKPTS_IRZ
         DO NB=1,NBANDSGW
#ifdef MPI
          IND = LOCAL_INDEX_CYCLIC( NCPU, NBANDSGW, NODE_ME, NB )
          IF ( IND < 1 ) CYCLE
#endif

             !get QP shift by linearization around E_DFT
             QP_LINEAR(:,NB,NK,ISP)=0.0_q

             QP_DIFF(NB,NK,ISP)=QP_FRAC(1,NB,NK,ISP)-QP_SVD(1,NB,NK,ISP)

             ! the criteria is that we try to keep Fractional version
             ! otherwise, we use the SVD version
             ! for linerization 
             IF (QP_FRAC(2,NB,NK,ISP)>=0.60_q .AND. &
               & QP_FRAC(2,NB,NK,ISP)<=1.0_q .AND. I_FRAC(1,NB,NK,ISP)==0) THEN
                QP_LINEAR(:,NB,NK,ISP)=QP_FRAC(:,NB,NK,ISP)
                SIGMA(NB,NK,ISP)=SIGMA_FRAC(NB,NK,ISP)
                TAG(NB,NK,ISP)=2
             ELSE IF (QP_SVD(2,NB,NK,ISP)>=0.60_q .AND. QP_SVD(2,NB,NK,ISP)<=1.0_q) THEN
                QP_LINEAR(:,NB,NK,ISP)=QP_SVD(:,NB,NK,ISP)
                SIGMA(NB,NK,ISP)=SIGMA_SVD(NB,NK,ISP)
                TAG(NB,NK,ISP)=1
             ELSE
                QP_LINEAR(:,NB,NK,ISP)=QP_FRAC(:,NB,NK,ISP)
                SIGMA(NB,NK,ISP)=SIGMA_FRAC(NB,NK,ISP)
                TAG(NB,NK,ISP)=2
             ENDIF

             ! for Zeros
             IF ( REAL(ZEROS_FRAC(2,NB,NK,ISP),q)>=0.60_q .AND. &
                & REAL(ZEROS_FRAC(2,NB,NK,ISP),q)<=1.0_q .AND. I_FRAC(2,NB,NK,ISP)==0) THEN
                QP_ZEROS(:,NB,NK,ISP)=ZEROS_FRAC(:,NB,NK,ISP)
             ELSE IF (REAL(ZEROS_SVD(2,NB,NK,ISP),q)>=0.60_q .AND. REAL(ZEROS_SVD(2,NB,NK,ISP),q)<=1.0_q) THEN
                QP_ZEROS(:,NB,NK,ISP)=ZEROS_SVD(:,NB,NK,ISP)
             ELSE
                QP_ZEROS(:,NB,NK,ISP)=ZEROS_FRAC(:,NB,NK,ISP)
             ENDIF

         ENDDO !NB
      ENDDO !NK
    ENDDO !ISP

    CALLMPI( M_sum_z(W%WDES%COMM,SIGMA,SIZE(SIGMA)) )
    CALLMPI( M_sum_d(W%WDES%COMM,QP_LINEAR,SIZE(QP_LINEAR)) )
    CALLMPI( M_sum_d(W%WDES%COMM,QP_DIFF,SIZE(QP_DIFF)) )
    CALLMPI( M_sum_z(W%WDES%COMM,QP_ZEROS,SIZE(QP_ZEROS)))
    CALLMPI( M_sum_i(W%WDES%COMM,TAG,SIZE(TAG)))

    IF (NODE_ME==1) THEN
      ! write out the QP_shift to OUTCAR
      CALL WRITE_QP_SHIFT(W,NKPTS_IRZ,NBANDSGW,TAG,IU6)
      ! write out the Pade fitted sigma(w) for specified bands and kpoints in vasprun.xml file
      IF (LPADE_FIT) THEN
         DO ISP=1,ISP_OUT
           CALL WRITE_SIGMA_XML(WMEAN%CELTOT,NOMEGA,NOMEGA_SVD,IW_UP,IW_DN,&
                  & EFERMI(ISP_INDEX(ISP)),2001,ISP_INDEX(ISP),IU0)
         ENDDO
      ENDIF
    ENDIF  

    CALL DEALLOCATE_ALL_PADE()
    DEALLOCATE(QP_SVD, ZEROS_SVD, SIGMA_SVD)
    DEALLOCATE(QP_FRAC, ZEROS_FRAC, SIGMA_FRAC, I_FRAC)
    DEALLOCATE(IW_UP,SELF_IW_UP)
    IF (LPADE_FIT) DEALLOCATE(IW_DN,SELF_IW_DN)


CONTAINS


!*********************************************************************
!
! get QP shift by linearization around E_DFT for NOMEGA_TMP case
!
!*********************************************************************
SUBROUTINE CALCULATE_QP(NOMEGA_TMP,QP_BB,SIGMA_BB)
    IMPLICIT NONE
    INTEGER     :: NOMEGA_TMP
    REAL(q)     :: QP_BB(:,:,:,:)
    COMPLEX(q)  :: SIGMA_BB(:,:,:)

    FRAC_COEFF_UP(:,NB,NK,ISP)=(0.0_q,0.0_q)
    CALL PADE_FRAC_COEFF(NOMEGA_TMP,IW_UP(NOMEGA-NOMEGA_TMP+1:NOMEGA),&
         & SELF_IW_UP(NOMEGA-NOMEGA_TMP+1:NOMEGA),FRAC_COEFF_UP(1:NOMEGA_TMP,NB,NK,ISP),IU0)
    QP_BB(:,NB,NK,ISP)=0.0_q
    CALL GET_QP_LINEAR_FRAC(NOMEGA_TMP,W%CELTOT(NB,NK,ISP),&
          & IW_UP(NOMEGA-NOMEGA_TMP+1:NOMEGA),FRAC_COEFF_UP(1:NOMEGA_TMP,NB,NK,ISP),&
          & SIGMA_BB(NB,NK,ISP),QP_BB(:,NB,NK,ISP),IU0,NP,WMEAN%CELTOT(NB,NK,ISP))

END SUBROUTINE CALCULATE_QP

END SUBROUTINE QP_PADE_FIT


!*********************************************************************
!
! write self-energy to vasprun.xml for G0W0(iw) in screened_2e.F
! 2/3 of the odd part of the SIGMA_IMAG is replaced by SIGMA_SIN 
!*********************************************************************
SUBROUTINE WRITE_SIGMA_IW_XML(ISP,NK,NB,OMEGA,IMAG_GRIDS,SIGMA_IMAG,SIGMA_SIN,WMEAN_CELTOT,NELM)
    USE minimax_struct 
    USE tutor, ONLY: vtutor, isAlert, GWkweight
    IMPLICIT NONE
    INTEGER    :: ISP,NK,NB
    REAL(q)    :: OMEGA(:)
    TYPE(imag_grid_handle)  :: IMAG_GRIDS   !time and frequency grid handle
    COMPLEX(q) :: SIGMA_IMAG(:)             !complex self-energy at cos grids
    REAL(q)    :: SIGMA_SIN(:)              !odd part of self-energy at sin grids
    COMPLEX(q) :: WMEAN_CELTOT
    INTEGER    :: NELM
    !local
    INTEGER              :: I,NO,IDUM,NOMEGA
    REAL(q),ALLOCATABLE  :: SELF_IW(:,:),P(:,:),F(:)
    REAL(q)              :: FDER, RE
    CHARACTER(len=100)   :: TXT
    COMPLEX(q)           :: CDUM
    LOGICAL              :: LDUM
    CHARACTER(1)         :: CHARAC

    IF (SIZE(OMEGA) /= IMAG_GRIDS%NOMEGA) THEN
       CALL vtutor%error("ERROR in WRITE_SIGMA_IW_XML "//str(SIZE(OMEGA))//" "//str( IMAG_GRIDS%NOMEGA) )
    ENDIF
    NOMEGA=IMAG_GRIDS%NOMEGA

    ALLOCATE(P(2*NOMEGA,5)) ; P=0.0_q
    ALLOCATE(F(NOMEGA)) ; F=0.0_q

    ALLOCATE(SELF_IW(3,NOMEGA))

    !spline interpolate SIGMA_SIN with open boundary conditions
    P=0.0_q
    DO I=1,NOMEGA
       P(I,1)=-IMAG_GRIDS%FER_IM(NOMEGA-I+1)
       P(I,2)=-SIGMA_SIN(NOMEGA-I+1)
       P(I+NOMEGA,1)=IMAG_GRIDS%FER_IM(I)
       P(I+NOMEGA,2)=SIGMA_SIN(I)
    ENDDO

    !first, get the coefficients of spline fit
    CALL SPLCOF(P,2*IMAG_GRIDS%NOMEGA,2*IMAG_GRIDS%NOMEGA,10E30_q)
    !second, spline-interpolate G_SIN_SIN to cos grids
    F=0.0_q
    DO I=1,IMAG_GRIDS%NOMEGA
       CALL SPLVAL(IMAG_GRIDS%FER_RE(I),F(I),FDER,P,2*IMAG_GRIDS%NOMEGA,2*IMAG_GRIDS%NOMEGA)
    ENDDO

    !now we overwrite the first 2/3 of SIGMA_SIN data points with
    !interpolated data
    DO I=1,MAX(1,2*IMAG_GRIDS%NOMEGA/3)
       RE=REAL(SIGMA_IMAG(I),q)
       SIGMA_IMAG(I)=CMPLX(RE,F(I),KIND=q)
    ENDDO

    IF (ISMEMBER(ISP,ISP_INDEX) .AND. ISMEMBER(NK,KPTS_INDEX) .AND. ISMEMBER(NB,NBANDS_INDEX)) THEN
        DO NO=1,NOMEGA
           SELF_IW(1,NO)=OMEGA(NO)
           SELF_IW(2,NO)=REAL(SIGMA_IMAG(NO)+WMEAN_CELTOT,q)
           SELF_IW(3,NO)=AIMAG(SIGMA_IMAG(NO))
        ENDDO
        WRITE(TXT,"(4(A,I0))") " ISP= ", ISP, " NK= ", NK, " NB= ", NB, " iter= ", NELM
        CALL XML_INCAR(TRIM(TXT)//" KS basis HF_MEAN= ",'F',IDUM, REAL(WMEAN_CELTOT,q),CDUM,LDUM,CHARAC,1)
        CALL XML_VECARRAY("selfenergy along imaginary axis"//TRIM(TXT))
        CALL XML_ARRAY_REAL(SELF_IW,"(F17.10,' ')")
        CALL XML_CLOSE_TAG
    ENDIF

    DEALLOCATE(SELF_IW,P,F)

END SUBROUTINE WRITE_SIGMA_IW_XML
!*********************************************************************


!*********************************************************************
FUNCTION ISMEMBER(A,B)
    IMPLICIT NONE
    INTEGER :: A,I
    INTEGER :: B(:)
    LOGICAL :: ISMEMBER

    ISMEMBER=.FALSE.
    DO I=1,SIZE(B)
       IF (B(I)==A) THEN
           ISMEMBER=.TRUE.
       ENDIF
    ENDDO
END FUNCTION ISMEMBER
!*********************************************************************


!*********************************************************************
!
! write self-energy to vasprun.xml for ISP channel
!
!*********************************************************************
SUBROUTINE WRITE_SIGMA_XML(WMEAN_CELTOT,NOMEGA,NOMEGA_SVD,IW_UP,IW_DN,EFERMI,NGRIDS,ISP,IU0)
    IMPLICIT NONE
    COMPLEX(q) :: WMEAN_CELTOT(:,:,:)
    INTEGER    :: NOMEGA,NOMEGA_SVD
    COMPLEX(q) :: IW_UP(NOMEGA)
    COMPLEX(q) :: IW_DN(NOMEGA)
    REAL(q)    :: EFERMI
    INTEGER    :: NGRIDS 
    INTEGER    :: ISP,IU0
    !local
    INTEGER    :: I,J,NK,NB

    DO J=1,KPTS_OUT
       NK=KPTS_INDEX(J)
       DO I=1,NBANDS_OUT
          NB=NBANDS_INDEX(I)
          CALL WRITE_G_SIGMA_SVD(NOMEGA_SVD,SVD_COEFF_UP(1:NOMEGA_SVD,NB,NK,ISP),&
                  & SVD_COEFF_DN(1:NOMEGA_SVD,NB,NK,ISP),EFERMI,NGRIDS,NK,NB,ISP,&
                  & IU0,WMEAN_CELTOT(NB,NK,ISP))
          CALL WRITE_G_SIGMA_FRAC(NOMEGA,IW_UP(1:NOMEGA),FRAC_COEFF_UP(1:NOMEGA,NB,NK,ISP),&
                  & IW_DN(1:NOMEGA),FRAC_COEFF_DN(1:NOMEGA,NB,NK,ISP),&
                  & EFERMI,NGRIDS,NK,NB,ISP,IU0,WMEAN_CELTOT(NB,NK,ISP))
       ENDDO
    ENDDO

END SUBROUTINE WRITE_SIGMA_XML
!*********************************************************************


!*********************************************************************
!
! this helper routine is used to do an analytical continuation by Pade fit 
! in the old G0W0 code using Merzuk imaginary grids, fractional version
!
!*********************************************************************
SUBROUTINE QP_PADE_FIT_FRAC(NB,NK,ISP,NOMEGA,OMEGA,SIGMA_IMAG,EDFT,EFERMI,COEFF,&
                     & SIGMA_QP,QP_LIN,QP_ZERO,IU0,L_FRAC,WMEAN_CELTOT)
    IMPLICIT NONE
    INTEGER         :: NB,NK,ISP
    INTEGER         :: NOMEGA
    REAL(q)         :: OMEGA(NOMEGA)           ! imaginary frequency
    COMPLEX(q)      :: SIGMA_IMAG(NOMEGA)      ! self-energy
    REAL(q)         :: EDFT                    ! CELNEW
    REAL(q)         :: EFERMI                  ! fermi energy
    COMPLEX(q)      :: COEFF(NOMEGA)           ! pade coefficients
    COMPLEX(q)      :: SIGMA_QP                ! self-energy at E_DFT
    REAL(q)         :: QP_LIN(2)               ! linearized QP energy and Z
    COMPLEX(q)      :: QP_ZERO(2)              ! QP energy and Z of Zeros solution
    INTEGER         :: IU0
    LOGICAL         :: L_FRAC(2)
    COMPLEX(q)      :: WMEAN_CELTOT
    ! local variables 
    INTEGER,PARAMETER      :: NP=5             ! NP-point finite difference stencil
    INTEGER                :: NOMEGA_TMP
    REAL(q)                :: Z1,Z2,KK
    REAL(q)                :: QP_TMP1(2), QP_TMP2(2)
    COMPLEX(q)             :: C_DFT,EQP_TMP,E1,E2
    COMPLEX(q)             :: SIGMA_TMP1,SIGMA_TMP2
    COMPLEX(q),ALLOCATABLE :: IW(:)

    ALLOCATE(IW(NOMEGA))

    C_DFT=CMPLX(EDFT,0.0_q,KIND=q)
    IW=CMPLX(EFERMI,OMEGA,KIND=q)

    ! get QP shift by linearization around E_DFT
    ! intial try for both NOMEGA and NOMEGA-1 cases, and select the 'best' solution
    QP_LIN(:)=0.0_q

    ! for NOMEGA-1 cases
    NOMEGA_TMP=NOMEGA-1
    COEFF(:)=(0.0_q,0.0_q)
    CALL PADE_FRAC_COEFF(NOMEGA_TMP,IW(NOMEGA-NOMEGA_TMP+1:NOMEGA),&
         & SIGMA_IMAG(NOMEGA-NOMEGA_TMP+1:NOMEGA),COEFF(1:NOMEGA_TMP),IU0)
    QP_TMP2(:)=0.0_q
    CALL GET_QP_LINEAR_FRAC(NOMEGA_TMP,C_DFT,IW(NOMEGA-NOMEGA_TMP+1:NOMEGA),&
          & COEFF(1:NOMEGA_TMP),SIGMA_TMP2,QP_TMP2(:),IU0,NP,WMEAN_CELTOT)

    ! for NOMEGA cases
    NOMEGA_TMP=NOMEGA
    COEFF(:)=(0.0_q,0.0_q)
    CALL PADE_FRAC_COEFF(NOMEGA_TMP,IW(NOMEGA-NOMEGA_TMP+1:NOMEGA),&
         & SIGMA_IMAG(NOMEGA-NOMEGA_TMP+1:NOMEGA),COEFF(1:NOMEGA_TMP),IU0)
    QP_TMP1(:)=0.0_q
    CALL GET_QP_LINEAR_FRAC(NOMEGA_TMP,C_DFT,IW(NOMEGA-NOMEGA_TMP+1:NOMEGA),&
          & COEFF(1:NOMEGA_TMP),SIGMA_TMP1,QP_TMP1(:),IU0,NP,WMEAN_CELTOT)

    ! the criteria is that if Z fator for NOMEGA cases is "reasonable", say [0.6,1.0)
    ! otherwise choose the one that is closer to 0.8_q
    IF (0.6_q<=QP_TMP1(2) .AND. QP_TMP1(2)<1.0_q ) THEN
      QP_LIN(:)=QP_TMP1(:)
      SIGMA_QP=SIGMA_TMP1
    ELSE
      IF (ABS(QP_TMP1(2)-0.8_q) <= ABS(QP_TMP2(2)-0.8_q)) THEN
        QP_LIN(:)=QP_TMP1(:)
        SIGMA_QP=SIGMA_TMP1
      ELSE
        QP_LIN(:)=QP_TMP2(:)
        SIGMA_QP=SIGMA_TMP2
      ENDIF
    ENDIF

    ! get QP shift by finding the roots of QP equation using Brent-bisection method
    EQP_TMP=CMPLX(QP_LIN(1),0.0_q,KIND=q)

    NOMEGA_TMP=NOMEGA
    DO WHILE(.TRUE.)
       COEFF(:)=(0.0_q,0.0_q)
       CALL PADE_FRAC_COEFF(NOMEGA_TMP,IW(NOMEGA-NOMEGA_TMP+1:NOMEGA),&
            & SIGMA_IMAG(NOMEGA-NOMEGA_TMP+1:NOMEGA),COEFF(1:NOMEGA_TMP),IU0)
       QP_ZERO=(0.0_q,0.0_q)
       CALL GET_ZEROS_BRENT_FRAC(NOMEGA_TMP,ISP,NK,NB,EQP_TMP,&
            & IW(NOMEGA-NOMEGA_TMP+1:NOMEGA),COEFF(1:NOMEGA_TMP),QP_ZERO(:),IU0,WMEAN_CELTOT)
       IF (REAL(QP_ZERO(2),q)>ZMAX .OR. REAL(QP_ZERO(2),q)<ZMIN ) THEN
          NOMEGA_TMP=NOMEGA_TMP-1
          L_FRAC(2)=.FALSE.
       ELSE
          EXIT
       ENDIF
       IF (NOMEGA_TMP<MAX(NOMEGA-2,8) ) THEN
          IF (IU0>=0) WRITE(IU0,'(A,3I6)') "Warning: NOMEGA_TMP<MAX(NOMEGA-2,8) in GET_ZEROS_BRENT_FRAC",NB,NK,ISP
          EXIT
       ENDIF
    ENDDO

    ! if still problematic, we report error and use original one
    IF (NOMEGA_TMP<MAX(NOMEGA-2,8) .OR. REAL(QP_ZERO(2),q)>ZMAX .OR. REAL(QP_ZERO(2),q)<ZMIN ) THEN
       IF (IU0>=0) WRITE(IU0,'(A,3I6)') "Warning: GET_ZEROS_BRENT_FRAC in Z", NB,NK,ISP
       COEFF=(0.0_q,0.0_q)
       CALL PADE_FRAC_COEFF(NOMEGA,IW(1:NOMEGA),SIGMA_IMAG(1:NOMEGA),COEFF(1:NOMEGA),IU0)
       QP_ZERO=(0.0_q,0.0_q)
       CALL GET_ZEROS_BRENT_FRAC(NOMEGA,ISP,NK,NB,EQP_TMP,&
            & IW(1:NOMEGA),COEFF(1:NOMEGA),QP_ZERO(:),IU0, WMEAN_CELTOT)
    ENDIF

    DEALLOCATE(IW)
END SUBROUTINE QP_PADE_FIT_FRAC
!*********************************************************************


!*********************************************************************
!
! this helper routine is used to do an analytical continuation by Pade fit 
! in the old G0W0 code using Merzuk imaginary grids, SVD version
!
!*********************************************************************
SUBROUTINE QP_PADE_FIT_SVD(NB,NK,ISP,NOMEGA,OMEGA,SIGMA_IMAG,EDFT,EFERMI,&
                    & NOMEGA_SVD,COEFF,SIGMA_QP,QP_LIN,QP_ZERO,IU0,L_SVD,WMEAN_CELTOT)
    IMPLICIT NONE
    INTEGER         :: NB,NK,ISP
    INTEGER         :: NOMEGA
    REAL(q)         :: OMEGA(NOMEGA)           ! imaginary frequency
    COMPLEX(q)      :: SIGMA_IMAG(NOMEGA)      ! self-energy
    REAL(q)         :: EDFT                    ! CELNEW
    REAL(q)         :: EFERMI                  ! fermi energy
    INTEGER         :: NOMEGA_SVD
    COMPLEX(q)      :: COEFF(NOMEGA_SVD)       ! pade coefficients
    COMPLEX(q)      :: SIGMA_QP                ! self-energy at E_DFT
    REAL(q)         :: QP_LIN(2)               ! linearized QP energy and Z
    COMPLEX(q)      :: QP_ZERO(2)              ! QP energy and Z of Zeros solution
    INTEGER         :: IU0 
    LOGICAL         :: L_SVD
    COMPLEX(q)      :: WMEAN_CELTOT
    ! local variables 
    INTEGER,PARAMETER      :: NP=5             ! NP-point finite difference stencil
    INTEGER                :: NOMEGA_TMP
    REAL(q)                :: Z1,Z2,KK
    COMPLEX(q)             :: C_DFT,EQP_TMP,E1,E2
    COMPLEX(q),ALLOCATABLE :: IW(:)

    ALLOCATE(IW(NOMEGA))

    C_DFT=CMPLX(EDFT,0.0_q,KIND=q)
    IW=CMPLX(EFERMI,OMEGA,KIND=q)

    ! get QP shift by linearization around E_DFT for SVD case
    COEFF(:)=(0.0_q,0.0_q)
    IF (MOD(NOMEGA,2)==0) THEN
       CALL PADE_LINEAR_COEFF(NOMEGA,IW(1:NOMEGA),SIGMA_IMAG(1:NOMEGA),&
                         & COEFF(1:NOMEGA),IU0)
    ELSE
       CALL PADE_SVD_COEFF(NOMEGA,IW(1:NOMEGA),SIGMA_IMAG(1:NOMEGA),&
                     & NOMEGA_SVD,COEFF(1:NOMEGA_SVD),IU0)
    ENDIF

    QP_LIN(:)=0.0_q
    SIGMA_QP=(0.0_q,0.0_q)
    CALL GET_QP_LINEAR_SVD(NOMEGA_SVD,C_DFT,&
          & COEFF(1:NOMEGA_SVD),SIGMA_QP,QP_LIN(:),IU0,NP,WMEAN_CELTOT)

    ! if fails, get QP shift by QR method
    IF ( MOD(NOMEGA,2)==0 .AND. (QP_LIN(2)>ZMAX .OR. QP_LIN(2)<ZMIN)) THEN
       COEFF(:)=(0.0_q,0.0_q)
       CALL PADE_SVD_COEFF(NOMEGA,IW(1:NOMEGA),SIGMA_IMAG(1:NOMEGA),&
                     & NOMEGA_SVD,COEFF(1:NOMEGA_SVD),IU0)
       QP_LIN(:)=0.0_q
       SIGMA_QP=(0.0_q,0.0_q)
       CALL GET_QP_LINEAR_SVD(NOMEGA,C_DFT,COEFF(1:NOMEGA),SIGMA_QP,QP_LIN(:),IU0,NP,WMEAN_CELTOT)
    ENDIF

    ! get QP shift by finding the roots of QP equation using Brent-bisection method
    EQP_TMP=CMPLX(QP_LIN(1),0.0_q,KIND=q)
    QP_ZERO=(0.0_q,0.0_q)

    COEFF(:)=(0.0_q,0.0_q)
    IF (MOD(NOMEGA,2)==0) THEN
       CALL PADE_LINEAR_COEFF(NOMEGA,IW(1:NOMEGA),SIGMA_IMAG(1:NOMEGA),&
                         & COEFF(1:NOMEGA),IU0)
    ELSE
       CALL PADE_SVD_COEFF(NOMEGA,IW(1:NOMEGA),SIGMA_IMAG(1:NOMEGA),&
                     & NOMEGA_SVD,COEFF(1:NOMEGA_SVD),IU0)
    ENDIF
    QP_ZERO=(0.0_q,0.0_q)
    CALL GET_ZEROS_BRENT_SVD(NOMEGA_SVD,ISP,NK,NB,EQP_TMP,&
             & IW(1:NOMEGA_SVD),COEFF(1:NOMEGA_SVD),QP_ZERO(:),IU0,WMEAN_CELTOT)

    DEALLOCATE(IW)
END SUBROUTINE QP_PADE_FIT_SVD
!*********************************************************************


!*********************************************************************
!
! write QP energies to OUTCAR 
!
!*********************************************************************
SUBROUTINE WRITE_QP_SHIFT(W,NKPTS_IRZ,NBANDSGW,TAG,IU6)
    IMPLICIT NONE
    TYPE(wavespin)    :: W    
    INTEGER           :: NKPTS_IRZ,NBANDSGW
    INTEGER           :: TAG(:,:,:)
    INTEGER           :: ISP,NK,NB,IU6

    DO ISP=1, W%WDES%ISPIN
       IF (W%WDES%ISPIN==2 .AND. IU6>=0) WRITE(IU6,'(/A,I1)') ' spin component ',ISP
       DO NK=1,NKPTS_IRZ
          IF (W%WDES%WTKPT(NK)==0) CYCLE ! zero k-point weighted points are useless right now in the GW code
          IF (IU6>=0) THEN
             WRITE(IU6,777) NK,W%WDES%VKPT(1,NK),W%WDES%VKPT(2,NK),W%WDES%VKPT(3,NK)
          ENDIF
          DO NB=1,NBANDSGW
             !write the output to OUTCAR file         
             IF (IU6>=0)  WRITE(IU6,888) NB,REAL(W%CELTOT(NB,NK,ISP),q),REAL(SIGMA(NB,NK,ISP),q),QP_LINEAR(1,NB,NK,ISP),&
                           & QP_LINEAR(2,NB,NK,ISP),REAL(QP_ZEROS(1,NB,NK,ISP),q), REAL(QP_ZEROS(2,NB,NK,ISP),q),&
                           & W%FERTOT(NB,NK,ISP)*W%WDES%RSPIN,AIMAG(QP_ZEROS(1,NB,NK,ISP)), QP_DIFF(NB,NK,ISP), TAG(NB,NK,ISP)
          ENDDO
       ENDDO
    ENDDO

777  FORMAT(/' k-point ',I3,' :',3X,3F10.4/  &
  & " band No.  KS-energies   sigma(KS)    QP-e(linear)    Z         QP-e(zeros)     Z        occupation    Imag(E_QP)    QP_DIFF TAG"/)
888  FORMAT(3X,I4,3X,9(F10.4,3X),I1)

END SUBROUTINE WRITE_QP_SHIFT
!*********************************************************************


!*********************************************************************
!
! help routine to write self-energy to vasprun.xml for SVD version 
!
!*********************************************************************
SUBROUTINE WRITE_G_SIGMA_SVD(NOMEGA,COEFF_UP,COEFF_DN,EFERMI,NOUT,NK,NB,ISP,IU0,WMEAN_CELTOT)
    IMPLICIT NONE
    INTEGER                 :: NOMEGA
    COMPLEX(q)              :: COEFF_UP(NOMEGA)
    COMPLEX(q)              :: COEFF_DN(NOMEGA)
    REAL(q)                 :: EFERMI
    INTEGER                 :: NOUT,NK,NB,ISP,IU0
    COMPLEX(q)              :: WMEAN_CELTOT
    ! local
    INTEGER                 :: NO
    REAL(q)                 :: RFREQ,RMIN,RMAX
    COMPLEX(q),ALLOCATABLE  :: RFREQ_UP(:)
    COMPLEX(q),ALLOCATABLE  :: RFREQ_DN(:)
    COMPLEX(q),ALLOCATABLE  :: SELF_W_UP(:)
    COMPLEX(q),ALLOCATABLE  :: SELF_W_DN(:)
    COMPLEX(q),ALLOCATABLE  :: CFREQ_UP(:)
    COMPLEX(q),ALLOCATABLE  :: SELF_UP(:)
    REAL(q),ALLOCATABLE     :: SIGMA_W_FIT(:,:)
    REAL(q),ALLOCATABLE     :: G_W_FIT(:,:)
    REAL(q),ALLOCATABLE     :: SIGMA_IW_FIT(:,:)
    CHARACTER(len=100)      :: TXT

    WRITE(TXT,"(3(A,I0),A)") " ISP=",ISP," NK=", NK, " NB=", NB, " SVD "

    ALLOCATE( RFREQ_UP(NOUT),  SELF_W_UP(NOUT) )
    ALLOCATE( RFREQ_DN(NOUT),  SELF_W_DN(NOUT) )
    ALLOCATE( CFREQ_UP(NOUT),  SELF_UP(NOUT)   )
    ALLOCATE( SIGMA_W_FIT(3,2*NOUT-1) ) 
    ALLOCATE( G_W_FIT(3,2*NOUT-1)     ) 
    ALLOCATE( SIGMA_IW_FIT(3,NOUT)    ) 

    RMIN=-100.0_q
    RMAX= 100.0_q
    DO NO=1,NOUT
       RFREQ=0.0_q+(RMAX-0.0_q)/(NOUT-1)*(NO-1)
       RFREQ_UP(NO)=CMPLX(RFREQ,0.0_q,KIND=q)
       RFREQ=0.0_q+300.0_q/(NOUT-1)*(NO-1)
       CFREQ_UP(NO)=CMPLX(0.0_q,RFREQ,KIND=q)
    ENDDO

    DO NO=1,NOUT
       RFREQ=RMIN+(0.0_q-RMIN)/(NOUT-1)*(NO-1)
       RFREQ_DN(NO)=CMPLX(RFREQ,0.0_q,KIND=q)
    ENDDO

    SIGMA_W_FIT(1,1:NOUT)=  REAL(RFREQ_DN,q)
    SIGMA_W_FIT(1,NOUT+1:2*NOUT-1)=  REAL(RFREQ_UP(2:NOUT),q)
    SIGMA_IW_FIT(1,:)=AIMAG(CFREQ_UP)

    CALL PADE_SVD_EVAL(SIZE(COEFF_UP),COEFF_UP,SIZE(RFREQ_UP),RFREQ_UP,SELF_W_UP,IU0,WMEAN_CELTOT)
    CALL PADE_SVD_EVAL(SIZE(COEFF_DN),COEFF_DN,SIZE(RFREQ_DN),RFREQ_DN,SELF_W_DN,IU0,WMEAN_CELTOT)

    SIGMA_W_FIT(2,1:NOUT)= REAL(SELF_W_DN,q)
    SIGMA_W_FIT(2,NOUT+1:2*NOUT-1)=REAL(SELF_W_UP(2:NOUT),q)
    SIGMA_W_FIT(3,1:NOUT)=AIMAG(SELF_W_DN)
    SIGMA_W_FIT(3,NOUT+1:2*NOUT-1)=AIMAG(SELF_W_UP(2:NOUT))

    CALL SIGMA_W_TO_G_W(SIGMA_W_FIT,G_W_FIT)

    CALL PADE_SVD_EVAL(SIZE(COEFF_UP),COEFF_UP,SIZE(CFREQ_UP),CFREQ_UP+EFERMI,SELF_UP,IU0,WMEAN_CELTOT)
    SIGMA_IW_FIT(2,:)= REAL(SELF_UP,q)
    SIGMA_IW_FIT(3,:)=AIMAG(SELF_UP)

    CALL XML_VECARRAY("FIT_SIGMA_W along real axis"//TRIM(TXT))
    CALL XML_ARRAY_REAL( SIGMA_W_FIT )
    CALL XML_CLOSE_TAG

    CALL XML_VECARRAY("FIT_G_W along real axis"//TRIM(TXT))
    CALL XML_ARRAY_REAL( G_W_FIT )
    CALL XML_CLOSE_TAG

    CALL XML_VECARRAY("FIT_SIGMA_IW along imaginary axis"//TRIM(TXT))
    CALL XML_ARRAY_REAL( SIGMA_IW_FIT )
    CALL XML_CLOSE_TAG

    DEALLOCATE( SIGMA_W_FIT, G_W_FIT, SIGMA_IW_FIT ) 
    DEALLOCATE( RFREQ_UP, SELF_W_UP ) 
    DEALLOCATE( RFREQ_DN, SELF_W_DN ) 
    DEALLOCATE( CFREQ_UP, SELF_UP ) 
END SUBROUTINE WRITE_G_SIGMA_SVD


!*********************************************************************
!
! help routine to write self-energy to vasprun.xml for fractional version 
!
!*********************************************************************
SUBROUTINE WRITE_G_SIGMA_FRAC(NOMEGA,IW_UP,COEFF_UP,IW_DN,COEFF_DN,&
                 & EFERMI,NOUT,NK,NB,ISP,IU0,WMEAN_CELTOT)
    IMPLICIT NONE
    INTEGER                 :: NOMEGA
    COMPLEX(q)              :: IW_UP(NOMEGA),COEFF_UP(NOMEGA)
    COMPLEX(q)              :: IW_DN(NOMEGA),COEFF_DN(NOMEGA)
    REAL(q)                 :: EFERMI
    INTEGER                 :: NOUT,NK,NB,ISP,IU0
    COMPLEX(q)              :: WMEAN_CELTOT
    ! local
    INTEGER                 :: NO
    REAL(q)                 :: RFREQ,RMIN,RMAX
    COMPLEX(q),ALLOCATABLE  :: RFREQ_UP(:)
    COMPLEX(q),ALLOCATABLE  :: RFREQ_DN(:)
    COMPLEX(q),ALLOCATABLE  :: SELF_W_UP(:)
    COMPLEX(q),ALLOCATABLE  :: SELF_W_DN(:)
    COMPLEX(q),ALLOCATABLE  :: CFREQ_UP(:)
    COMPLEX(q),ALLOCATABLE  :: SELF_UP(:)
    REAL(q),ALLOCATABLE     :: SIGMA_W_FIT(:,:)
    REAL(q),ALLOCATABLE     :: G_W_FIT(:,:)
    REAL(q),ALLOCATABLE     :: SIGMA_IW_FIT(:,:)
    CHARACTER(len=100)      :: TXT

    WRITE(TXT,"(3(A,I0),A)") " ISP=",ISP," NK=", NK, " NB=", NB, " FRAC "

    ALLOCATE( RFREQ_UP(NOUT),  SELF_W_UP(NOUT) )
    ALLOCATE( RFREQ_DN(NOUT),  SELF_W_DN(NOUT) )
    ALLOCATE( CFREQ_UP(NOUT),  SELF_UP(NOUT)   )
    ALLOCATE( SIGMA_W_FIT(3,2*NOUT-1) ) 
    ALLOCATE( G_W_FIT(3,2*NOUT-1)     ) 
    ALLOCATE( SIGMA_IW_FIT(3,NOUT)    ) 

    RMIN=-100.0_q
    RMAX= 100.0_q
    DO NO=1,NOUT
       RFREQ=0.0_q+(RMAX-0.0_q)/(NOUT-1)*(NO-1)
       RFREQ_UP(NO)=CMPLX(RFREQ,0.0_q,KIND=q)
       RFREQ=0.0_q+300.0_q/(NOUT-1)*(NO-1)
       CFREQ_UP(NO)=CMPLX(0.0_q,RFREQ,KIND=q)
    ENDDO

    DO NO=1,NOUT
       RFREQ=RMIN+(0.0_q-RMIN)/(NOUT-1)*(NO-1)
       RFREQ_DN(NO)=CMPLX(RFREQ,0.0_q,KIND=q)
    ENDDO
 
    !fitst collect the frequencies
    SIGMA_W_FIT(1,1:NOUT)=  REAL(RFREQ_DN,q)
    SIGMA_W_FIT(1,NOUT+1:2*NOUT-1)=  REAL(RFREQ_UP(2:NOUT),q)
    SIGMA_IW_FIT(1,:)=AIMAG(CFREQ_UP)

    CALL PADE_FRAC_EVAL(SIZE(IW_UP),IW_UP,COEFF_UP,NOUT,RFREQ_UP,SELF_W_UP,IU0,WMEAN_CELTOT)
    CALL PADE_FRAC_EVAL(SIZE(IW_DN),IW_DN,COEFF_DN,NOUT,RFREQ_DN,SELF_W_DN,IU0,WMEAN_CELTOT)

    SIGMA_W_FIT(2,1:NOUT)= REAL(SELF_W_DN,q)
    SIGMA_W_FIT(2,NOUT+1:2*NOUT-1)=REAL(SELF_W_UP(2:NOUT),q)
    SIGMA_W_FIT(3,1:NOUT)=AIMAG(SELF_W_DN)
    SIGMA_W_FIT(3,NOUT+1:2*NOUT-1)=AIMAG(SELF_W_UP(2:NOUT))

    CALL SIGMA_W_TO_G_W(SIGMA_W_FIT,G_W_FIT)

    CALL PADE_FRAC_EVAL(SIZE(IW_UP),IW_UP,COEFF_UP,NOUT,CFREQ_UP+EFERMI,SELF_UP,IU0,WMEAN_CELTOT)
    SIGMA_IW_FIT(2,:)= REAL(SELF_UP,q)
    SIGMA_IW_FIT(3,:)=AIMAG(SELF_UP)

    CALL XML_VECARRAY("FIT_SIGMA_W along real axis"//TRIM(TXT))
    CALL XML_ARRAY_REAL( SIGMA_W_FIT )
    CALL XML_CLOSE_TAG

    CALL XML_VECARRAY("FIT_G_W along real axis"//TRIM(TXT))
    CALL XML_ARRAY_REAL( G_W_FIT )
    CALL XML_CLOSE_TAG

    CALL XML_VECARRAY("FIT_SIGMA_IW along imaginary axis"//TRIM(TXT))
    CALL XML_ARRAY_REAL( SIGMA_IW_FIT )
    CALL XML_CLOSE_TAG

    DEALLOCATE( SIGMA_W_FIT, G_W_FIT, SIGMA_IW_FIT ) 
    DEALLOCATE( RFREQ_UP, SELF_W_UP ) 
    DEALLOCATE( RFREQ_DN, SELF_W_DN ) 
    DEALLOCATE( CFREQ_UP, SELF_UP ) 
END SUBROUTINE WRITE_G_SIGMA_FRAC


!*********************************************************************
!
! calculation of G(w) from Sigma(w)
! G(w)=1/(w-Sigma(w))
!
!*********************************************************************
SUBROUTINE SIGMA_W_TO_G_W(SIGMA_W,GREEN_W)
    IMPLICIT NONE
    REAL(q)    :: SIGMA_W(:,:), GREEN_W(:,:)
    ! local
    COMPLEX(q) :: OMEGA, SELF, G, ONE
    INTEGER    :: I

    ONE=CMPLX(1.0_q,0.0_q,q)
    DO I=1,SIZE(SIGMA_W,2)
       OMEGA=CMPLX(SIGMA_W(1,I),0.0_q,KIND=q)
       SELF=CMPLX(SIGMA_W(2,I),SIGMA_W(3,I),KIND=q)
       G=ONE/(OMEGA-SELF)
       GREEN_W(1,I)=SIGMA_W(1,I)
       GREEN_W(2,I)=REAL(G)
       GREEN_W(3,I)=AIMAG(G)
    ENDDO
END SUBROUTINE SIGMA_W_TO_G_W


!*********************************************************************
!
! calculation of Sigma(iw) from G(iw)
! Sigma(iw)=iw+EF-1/G(iw)
!
!*********************************************************************
SUBROUTINE G_IW_TO_SIGMA_IW(FREQ,G,SIGMA,EFERMI)
    IMPLICIT NONE
    REAL(q)     :: FREQ        !freqency
    COMPLEX(q)  :: SIGMA       !self-energy SIGMA
    COMPLEX(q)  :: G           !green function G
    REAL(q)     :: EFERMI
    SIGMA=cmplx(0.0,FREQ,q)+EFERMI-1.0_q/G
END SUBROUTINE G_IW_TO_SIGMA_IW
!*********************************************************************


!*********************************************************************
!
! calculate the QP energies by linearization around DFT eigenvalues
! Fractional continuation version
!
!*********************************************************************
SUBROUTINE GET_QP_LINEAR_FRAC(NOMEGA,EDFT,IW,COEFF,SIGMA,QP,IU0,NP,WMEAN_CELTOT)
    USE tutor, ONLY: vtutor, isAlert, GWkweight
    IMPLICIT NONE
    INTEGER               :: NOMEGA
    COMPLEX(q)            :: EDFT
    COMPLEX(q),INTENT(IN) :: IW(NOMEGA)   
    COMPLEX(q),INTENT(IN) :: COEFF(NOMEGA)  
    COMPLEX(q)            :: SIGMA        ! self-energy at DFT energy (including T+V_H)
    REAL(q)               :: QP(2)        ! approximate QP energy
    INTEGER               :: IU0
    INTEGER               :: NP           ! NP-point ifferentiation stencils
    COMPLEX(q)            :: WMEAN_CELTOT
    ! local 
    REAL(q)               :: KK
    COMPLEX(q)            :: E_DFT(NP),H_DFT(NP)
    COMPLEX(q),PARAMETER  :: DELTA=(0.1_q,0.0_q) !0.1 seems best

    !quick return
    IF (NP/=3 .AND. NP/=5 .AND. NP/=7 ) THEN
      CALL vtutor%error("ERROR: stencils NP in GET_QP_LINEAR " // str(NP))
    ENDIF    
    
    SELECT CASE (NP)
      CASE(3)
        E_DFT(1)=EDFT
        E_DFT(2)=E_DFT(1)+DELTA  
        E_DFT(3)=E_DFT(1)-DELTA  

        CALL PADE_FRAC_EVAL(NOMEGA,IW,COEFF,NP,E_DFT,H_DFT,IU0,WMEAN_CELTOT)
        SIGMA=H_DFT(1)

        KK=REAL(H_DFT(2)-H_DFT(3),q)/(2.0_q*REAL(DELTA,q))
        QP(2)=1.0_q/(1.0_q-KK)
        QP(1)=REAL(SIGMA-EDFT,q)*QP(2)+REAL(EDFT,q)

      CASE(5)
        E_DFT(1)=EDFT
        E_DFT(2)=E_DFT(1)+DELTA  
        E_DFT(3)=E_DFT(1)+DELTA*2.0_q
        E_DFT(4)=E_DFT(1)-DELTA  
        E_DFT(5)=E_DFT(1)-DELTA*2.0_q  

        CALL PADE_FRAC_EVAL(NOMEGA,IW,COEFF,NP,E_DFT,H_DFT,IU0,WMEAN_CELTOT)
        SIGMA=H_DFT(1)

        KK=REAL(-1.0_q*H_DFT(3)+8.0_q*H_DFT(2)-8.0_q*H_DFT(4)+1.0_q*H_DFT(5),q)/(12.0_q*REAL(DELTA,q))
        QP(2)=1.0_q/(1.0_q-KK)
        QP(1)=REAL(SIGMA-EDFT,q)*QP(2)+REAL(EDFT,q)

      CASE(7)
        E_DFT(1)=EDFT
        E_DFT(2)=E_DFT(1)+DELTA  
        E_DFT(3)=E_DFT(1)+DELTA*2.0_q
        E_DFT(4)=E_DFT(1)+DELTA*3.0_q
        E_DFT(5)=E_DFT(1)-DELTA  
        E_DFT(6)=E_DFT(1)-DELTA*2.0_q  
        E_DFT(7)=E_DFT(1)-DELTA*3.0_q  

        CALL PADE_FRAC_EVAL(NOMEGA,IW,COEFF,NP,E_DFT,H_DFT,IU0,WMEAN_CELTOT)
        SIGMA=H_DFT(1)

        KK=REAL(0.016666667_q*H_DFT(4)-0.15_q*H_DFT(3)+0.75_q*H_DFT(2)-0.75_q*H_DFT(5)+&
               & 0.15_q*H_DFT(6)-0.016666667_q*H_DFT(7),q)/REAL(DELTA,q)
        QP(2)=1.0_q/(1.0_q-KK)
        QP(1)=REAL(SIGMA-EDFT,q)*QP(2)+REAL(EDFT,q)
     END SELECT

END SUBROUTINE GET_QP_LINEAR_FRAC


!*********************************************************************
!
! calculate the QP energies by linearization around DFT eigenvalues
! Pade SVD verstion 
!
!*********************************************************************
SUBROUTINE GET_QP_LINEAR_SVD(NOMEGA,EDFT,COEFF,SIGMA,QP,IU0,NP,WMEAN_CELTOT)
    USE tutor, ONLY: vtutor, isAlert, GWkweight
    IMPLICIT NONE
    INTEGER               :: NOMEGA
    COMPLEX(q)            :: EDFT
    COMPLEX(q),INTENT(IN) :: COEFF(NOMEGA)
    COMPLEX(q)            :: SIGMA        ! self-energy at DFT energy (including T+V_H)
    REAL(q)               :: QP(2)        ! approximate QP energy
    INTEGER               :: IU0
    INTEGER               :: NP           ! NP-point ifferentiation stencils
    COMPLEX(q)            :: WMEAN_CELTOT
    ! local 
    REAL(q)               :: KK
    COMPLEX(q)            :: E_DFT(NP),H_DFT(NP)
    COMPLEX(q),PARAMETER  :: DELTA=(0.1_q,0.0_q) ! 0.1 seems best

    !quick return
    IF (NP/=3 .AND. NP/=5 .AND. NP/=7 ) THEN
      CALL vtutor%error("ERROR: stencils NP in GET_QP_LINEAR_FRAC "//str( NP) )
    ENDIF

    SELECT CASE (NP)
      CASE(3)
        E_DFT(1)=EDFT
        E_DFT(2)=E_DFT(1)+DELTA
        E_DFT(3)=E_DFT(1)-DELTA

        CALL PADE_SVD_EVAL(NOMEGA,COEFF,NP,E_DFT,H_DFT,IU0,WMEAN_CELTOT)
        SIGMA=H_DFT(1)

        KK=REAL(H_DFT(2)-H_DFT(3),q)/(2.0_q*REAL(DELTA,q))
        QP(2)=1.0_q/(1.0_q-KK)
        QP(1)=REAL(SIGMA-EDFT,q)*QP(2)+REAL(EDFT,q)

      CASE(5)
        E_DFT(1)=EDFT
        E_DFT(2)=E_DFT(1)+DELTA
        E_DFT(3)=E_DFT(1)+DELTA*2.0_q
        E_DFT(4)=E_DFT(1)-DELTA
        E_DFT(5)=E_DFT(1)-DELTA*2.0_q

        CALL PADE_SVD_EVAL(NOMEGA,COEFF,NP,E_DFT,H_DFT,IU0,WMEAN_CELTOT)
        SIGMA=H_DFT(1)

        KK=REAL(-1.0_q*H_DFT(3)+8.0_q*H_DFT(2)-8.0_q*H_DFT(4)+1.0_q*H_DFT(5),q)/(12.0_q*REAL(DELTA,q))
        QP(2)=1.0_q/(1.0_q-KK)
        QP(1)=REAL(SIGMA-EDFT,q)*QP(2)+REAL(EDFT,q)

      CASE(7)
        E_DFT(1)=EDFT
        E_DFT(2)=E_DFT(1)+DELTA
        E_DFT(3)=E_DFT(1)+DELTA*2.0_q
        E_DFT(4)=E_DFT(1)+DELTA*3.0_q
        E_DFT(5)=E_DFT(1)-DELTA
        E_DFT(6)=E_DFT(1)-DELTA*2.0_q
        E_DFT(7)=E_DFT(1)-DELTA*3.0_q

        CALL PADE_SVD_EVAL(NOMEGA,COEFF,NP,E_DFT,H_DFT,IU0,WMEAN_CELTOT)
        SIGMA=H_DFT(1)

        KK=REAL(0.016666667_q*H_DFT(4)-0.15_q*H_DFT(3)+0.75_q*H_DFT(2)-0.75_q*H_DFT(5)+&
              & 0.15_q*H_DFT(6)-0.016666667_q*H_DFT(7),q)/REAL(DELTA,q)
        QP(2)=1.0_q/(1.0_q-KK)
        QP(1)=REAL(SIGMA-EDFT,q)*QP(2)+REAL(EDFT,q)
     END SELECT

END SUBROUTINE GET_QP_LINEAR_SVD


!*********************************************************************
!
! calculate the QP energies by searching the root of E_QP=\Sigma(E_QP)
! via Brent scheme for FRAC version
!
!*********************************************************************
SUBROUTINE GET_ZEROS_BRENT_FRAC(NOMEGA,ISP,NK,NB,EDFT,IW,COEFF,QP,IU0,WMEAN_CELTOT)
    IMPLICIT NONE
    INTEGER               :: NOMEGA
    INTEGER               :: ISP,NK,NB
    COMPLEX(q)            :: EDFT          ! here is the sulution of linearized version
    COMPLEX(q),INTENT(IN) :: IW(NOMEGA)   
    COMPLEX(q),INTENT(IN) :: COEFF(NOMEGA)  
    COMPLEX(q)            :: QP(2)
    INTEGER               :: IU0
    COMPLEX(q)            :: WMEAN_CELTOT
    ! local 
    COMPLEX(q),PARAMETER  :: DELTA=(0.1_q,0.0_q)
    COMPLEX(q)            :: E_TMP(3),H_TMP(3)
    INTEGER,PARAMETER     :: NUM=100
    INTEGER               :: NTOT,I,J
    REAL(q)               :: X1,X2,DE,EPS
    REAL(q), DIMENSION(:),POINTER :: XB1,XB2
    REAL(q), ALLOCATABLE  :: ROOT(:),DIST(:)
    CHARACTER(len=100)    :: TXT
   
    NULLIFY(XB1,XB2)

    DE=2.0_q ; NTOT=0
    DO
      X1=REAL(EDFT,q)-DE
      X2=REAL(EDFT,q)+DE
      CALL FIND_BRACKET(IW,COEFF,X1,X2,NUM,XB1,XB2,NTOT,FUNC_FRAC,WMEAN_CELTOT)
      IF (NTOT==0) THEN
        DE=DE+0.5
      ELSE
        EXIT
      ENDIF
      IF (DE>6.0_q) EXIT
    ENDDO

    IF (NTOT==0) THEN
       IF (IU0>=0) WRITE(IU0,'(A,3I5)') &
                    & "internal error in VASP: NTOT==0,GET_ZEROS_BRENT_FRAC failed", NB,NK,ISP
       QP(1)=(-1000.0_q,-1000.0_q)
       QP(2)=(-1000.0_q,-1000.0_q)
       RETURN
    ENDIF

    ALLOCATE(ROOT(NTOT),DIST(NTOT))

    DO I=1,NTOT
      ROOT(I)=ROOT_FIND_BRENT(IW,COEFF,XB1(I),XB2(I),0.001_q,FUNC_FRAC,WMEAN_CELTOT)
      IF (NTOT==1 .AND. ABS(ROOT(1)-43210.0_q)<0.0001_q) THEN
          QP(1)=(-1000.0_q,-1000.0_q)
          QP(2)=(-1000.0_q,-1000.0_q)
          RETURN
      ENDIF
      DIST(I)=ABS(ROOT(I)-REAL(EDFT,q))
    ENDDO

    J=1
    EPS=10000.0_q
    DO I=1,NTOT
       IF (DIST(I)<EPS) THEN 
          EPS=DIST(I)  
          J=I
       ENDIF
    ENDDO

#ifdef debug
444 format("(A,2X,",I0,"F9.3,4X,A,F9.3)")
555 format(3(A,I5,2X),4X,A,F9.3)
    WRITE(TXT,444) NTOT
    WRITE(*,555) "NK=",NK, "NB=",NB, "NROOT=",NTOT, "E_QP(linear)=", REAL(EDFT,q)
    WRITE(*,trim(TXT)) "TOTAL ROOTS=", ROOT, "FINAL_QP=", ROOT(J)
    WRITE(*,*)
#endif

    ! determine the Z factor
    E_TMP=0.0_q ; H_TMP=0.0_q
    E_TMP(1)=ROOT(J)
    E_TMP(2)=E_TMP(1)+DELTA
    E_TMP(3)=E_TMP(1)-DELTA

    CALL PADE_FRAC_EVAL(NOMEGA,IW,COEFF,3,E_TMP,H_TMP,IU0,WMEAN_CELTOT)
    QP(1)=H_TMP(1)
    QP(2)=1.0_q/(1-REAL(H_TMP(2)-H_TMP(3),q)/(2.0_q*REAL(DELTA,q)) )
 
    DEALLOCATE(ROOT,DIST)
    IF (ASSOCIATED(XB1)) DEALLOCATE(XB1)
    IF (ASSOCIATED(XB2)) DEALLOCATE(XB2)
END SUBROUTINE GET_ZEROS_BRENT_FRAC


!*********************************************************************
!
! calculate the QP energies by searching the root of E_QP=\Sigma(E_QP)
! via Brent scheme for SVD version 
!
!*********************************************************************
SUBROUTINE GET_ZEROS_BRENT_SVD(NOMEGA,ISP,NK,NB,EDFT,IW,COEFF,QP,IU0,WMEAN_CELTOT)
    IMPLICIT NONE
    INTEGER               :: NOMEGA
    INTEGER               :: ISP,NK,NB
    COMPLEX(q)            :: EDFT            ! here is the sulution of linearized version
    COMPLEX(q),INTENT(IN) :: IW(NOMEGA)   
    COMPLEX(q),INTENT(IN) :: COEFF(NOMEGA)  
    COMPLEX(q)            :: QP(2)
    INTEGER               :: IU0
    COMPLEX(q)            :: WMEAN_CELTOT
    ! local 
    COMPLEX(q),PARAMETER  :: DELTA=(0.1_q,0.0_q)
    COMPLEX(q)            :: E_TMP(3),H_TMP(3)
    INTEGER,PARAMETER     :: NUM=100
    INTEGER               :: NTOT,I,J
    REAL(q)               :: X1,X2,DE,EPS
    REAL(q), DIMENSION(:),POINTER :: XB1,XB2
    REAL(q), ALLOCATABLE  :: ROOT(:),DIST(:)
    CHARACTER(len=100)    :: TXT
   
    NULLIFY(XB1,XB2)

    DE=2.0_q ; NTOT=0
    DO  
      X1=REAL(EDFT,q)-DE
      X2=REAL(EDFT,q)+DE
      CALL FIND_BRACKET(IW,COEFF,X1,X2,NUM,XB1,XB2,NTOT,FUNC_SVD,WMEAN_CELTOT)
      IF (NTOT==0) THEN
        DE=DE+0.5
      ELSE
        EXIT
      ENDIF
      IF (DE>6.0_q) EXIT
    ENDDO

    IF (NTOT==0) THEN
       IF (IU0>=0) WRITE(IU0,'(A,3I5)') &
                    & "internal error in VASP: NTOT==0,GET_ZEROS_BRENT_FRAC failed", NB,NK,ISP
       QP(1)=(-1000.0_q,-1000.0_q)
       QP(2)=(-1000.0_q,-1000.0_q)
       RETURN
    ENDIF

    ALLOCATE(ROOT(NTOT),DIST(NTOT))

    DO I=1,NTOT
      ROOT(I)=ROOT_FIND_BRENT(IW,COEFF,XB1(I),XB2(I),0.001_q,FUNC_SVD,WMEAN_CELTOT)
      IF (NTOT==1 .AND. ABS(ROOT(1)-43210.0_q)<0.0001_q) THEN
          QP(1)=(-1000.0_q,-1000.0_q)
          QP(2)=(-1000.0_q,-1000.0_q)
          RETURN
      ENDIF
      DIST(I)=ABS(ROOT(I)-REAL(EDFT,q))
    ENDDO

    J=1
    EPS=10000.0_q
    DO I=1,NTOT
       IF (DIST(I)<EPS) THEN 
          EPS=DIST(I)  
          J=I
       ENDIF
    ENDDO

#ifdef debug
444 format("(A,2X,",I0,"F9.3,4X,A,F9.3)")
555 format(3(A,I5,2X),4X,A,F9.3)
    WRITE(TXT,444) NTOT
    WRITE(*,555) "NK=",NK, "NB=",NB, "NROOT=",NTOT, "E_QP(linear)=", REAL(EDFT,q)
    WRITE(*,trim(TXT)) "TOTAL ROOTS=", ROOT, "FINAL_QP=", ROOT(J)
    WRITE(*,*)
#endif

    ! determine the Z factor
    E_TMP=0.0_q ; H_TMP=0.0_q
    E_TMP(1)=ROOT(J)
    E_TMP(2)=E_TMP(1)+DELTA
    E_TMP(3)=E_TMP(1)-DELTA

    CALL PADE_SVD_EVAL(SIZE(COEFF),COEFF,3,E_TMP,H_TMP,IU0,WMEAN_CELTOT)
    QP(1)=H_TMP(1)
    QP(2)=1.0_q/(1-REAL(H_TMP(2)-H_TMP(3),q)/(2.0_q*REAL(DELTA,q)) )
 
    DEALLOCATE(ROOT,DIST)
    IF (ASSOCIATED(XB1)) DEALLOCATE(XB1)
    IF (ASSOCIATED(XB2)) DEALLOCATE(XB2)
END SUBROUTINE GET_ZEROS_BRENT_SVD
!*********************************************************************


!*********************************************************************
!
! find the bracket to search the root
! Given a function fx defined on the interval from x1-x2 subdivide the 
! interval into n equally spaced segments, and search for zero crossings 
! of the function. nb is input as the maximum number of roots sought, and 
! is reset to the number of bracketing pairs xb1(1:nb),xb2(1:nb) that 
! are found.  
!
!*********************************************************************
SUBROUTINE FIND_BRACKET(IW,COEFF,X1,X2,N,XB1,XB2,NB,FUNC,WMEAN_CELTOT)
    IMPLICIT NONE
    COMPLEX(q),INTENT(IN) :: IW(:)
    COMPLEX(q),INTENT(IN) :: COEFF(:)
    REAL(q), INTENT(IN)   :: X1,X2
    INTEGER, INTENT(IN)   :: N
    INTEGER, INTENT(OUT)  :: NB
    REAL(q), DIMENSION(:), POINTER :: XB1,XB2
    COMPLEX(q)            :: WMEAN_CELTOT
    REAL(q)    :: FUNC    ! function as input
    ! local
    INTEGER :: I
    REAL(q) :: DX
    REAL(q), DIMENSION(0:N) :: F,X
    LOGICAL, DIMENSION(1:N) :: MASK
    LOGICAL, SAVE :: INIT=.TRUE.

    IF (INIT) THEN
      INIT=.FALSE.
      NULLIFY(XB1,XB2)
    END IF

    IF (ASSOCIATED(XB1)) DEALLOCATE(XB1)
    IF (ASSOCIATED(XB2)) DEALLOCATE(XB2)
   
    DX=(X2-X1)/N
    X=X1+DX*ARTH_I(0,1,N+1)
    DO I=0,N
       F(I)=FUNC(SIZE(COEFF),IW(1:SIZE(COEFF)),COEFF,X(I),WMEAN_CELTOT)
    END DO

    MASK=(F(1:N)*F(0:N-1) <= 0.0)
    
    NB=COUNT(MASK)
    ALLOCATE(XB1(NB),XB2(NB))
    XB1(1:NB)=PACK(X(0:N-1),MASK)
    XB2(1:NB)=PACK(X(1:N),MASK)

END SUBROUTINE FIND_BRACKET


!*********************************************************************
!
! help routine for FIND_BRACKET
!
!*********************************************************************
FUNCTION ARTH_I(FIRST,INCREMENT,N)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: FIRST,INCREMENT,N
    INTEGER, DIMENSION(N) :: ARTH_I
    !LOCAL
    INTEGER, PARAMETER    :: NPAR_ARTH=16, NPAR2_ARTH=8
    INTEGER :: K,K2,TEMP
    
    IF (N > 0) ARTH_I(1)=FIRST
    IF (N <= NPAR_ARTH) THEN
       DO K=2,N
          ARTH_I(K)=ARTH_I(K-1)+INCREMENT
       END DO
    ELSE
       DO K=2,NPAR2_ARTH
          ARTH_I(K)=ARTH_I(K-1)+INCREMENT
       END DO
       TEMP=INCREMENT*NPAR2_ARTH
       K=NPAR2_ARTH
       DO
         IF (K >= N) EXIT
         K2=K+K
         ARTH_I(K+1:MIN(K2,N))=TEMP+ARTH_I(1:MIN(K,N-K))
         TEMP=TEMP+TEMP
         K=K2
       END DO
    END IF
END FUNCTION ARTH_I


!*********************************************************************
!
! find the root via Brent scheme
!
!*********************************************************************
FUNCTION ROOT_FIND_BRENT(IW,COEFF,X1,X2,TOL,FUNC,WMEAN_CELTOT)
    IMPLICIT NONE
    COMPLEX(q),INTENT(IN) :: IW(:)
    COMPLEX(q),INTENT(IN) :: COEFF(:)
    REAL(q), INTENT(IN)   :: X1,X2,TOL
    COMPLEX(q)            :: WMEAN_CELTOT
    REAL(q) :: ROOT_FIND_BRENT
    REAL(q) :: FUNC   ! function as input
    ! local
    INTEGER, PARAMETER :: ITMAX=100
    REAL(q), PARAMETER :: EPS=EPSILON(X1)
    INTEGER(q) :: ITER
    REAL(q)    :: A,B,C,D,E,FA,FB,FC
    REAL(q)    :: P,Q1,R,S,TOL1,XM

    A=X1
    B=X2
    FA=FUNC(SIZE(COEFF),IW(1:SIZE(COEFF)),COEFF,A,WMEAN_CELTOT)
    FB=FUNC(SIZE(COEFF),IW(1:SIZE(COEFF)),COEFF,B,WMEAN_CELTOT)
    IF ((FA > 0.0 .AND. FB > 0.0) .OR. (FA < 0.0 .AND. FB < 0.0)) THEN
        WRITE(*,*) 'ERROR: root must be bracketed for ROOT_FIND_BRENT'
        ROOT_FIND_BRENT=43210.0_q
        RETURN
    ENDIF    
    C=B
    FC=FB
    DO ITER=1,ITMAX
       IF ((FB > 0.0 .AND. FC > 0.0) .OR. (FB < 0.0 .AND. FC < 0.0)) THEN
                C=A
                FC=FA
                D=B-A
                E=D
       END IF
       IF (ABS(FC) < ABS(FB)) THEN
                A=B
                B=C
                C=A
                FA=FB
                FB=FC
                FC=FA
       END IF
       TOL1=2.0_q*EPS*ABS(B)+0.5_q*TOL
       XM=0.5_q*(C-B)
       IF (ABS(XM) <= TOL1 .OR. FB == 0.0) THEN
                ROOT_FIND_BRENT=B
                RETURN
       END IF
       IF (ABS(E) >= TOL1 .AND. ABS(FA) > ABS(FB)) THEN
                S=FB/FA
                IF (A == C) THEN
                        P=2.0_q*XM*S
                        Q1=1.0_q-S
                ELSE
                        Q1=FA/FC
                        R=FB/FC
                        P=S*(2.0_q*XM*Q1*(Q1-R)-(B-A)*(R-1.0_q))
                        Q1=(Q1-1.0_q)*(R-1.0_q)*(S-1.0_q)
                END IF
                IF (P > 0.0) Q1=-Q1
                P=ABS(P)
                IF (2.0_q*P  <  MIN(3.0_q*XM*Q1-ABS(TOL1*Q1),ABS(E*Q1))) THEN
                        E=D
                        D=P/Q1
                ELSE
                        D=XM
                        E=D
                END IF
       ELSE
          D=XM
          E=D
       END IF
       A=B
       FA=FB
       B=B+MERGE(D,SIGN(TOL1,XM), ABS(D) > TOL1 )
       FB=FUNC(SIZE(COEFF),IW(1:SIZE(COEFF)),COEFF,B,WMEAN_CELTOT)
    END DO
    WRITE(*,*) 'ERROR: ROOT_FIND_BRENT: exceeded maximum iterations: 100'
    ROOT_FIND_BRENT=43210.0_q
    RETURN
END FUNCTION ROOT_FIND_BRENT


!*********************************************************************
!
! FUNC_FRAC=SIGMA(w)-w for FRAC version
!
!*********************************************************************
FUNCTION FUNC_FRAC(N,IZ,A,X,WMEAN_CELTOT)
    IMPLICIT NONE
    INTEGER,   INTENT(IN)   :: N        ! number of imaginary data sets (iz(i),ih(i),i=1,...N))
    COMPLEX(q),INTENT(IN)   :: IZ(N)    ! imaginary frequencies iz(i)
    COMPLEX(q),INTENT(IN)   :: A(N)     ! coefficients of pade
    REAL(q),   INTENT(IN)   :: X
    COMPLEX(q)              :: WMEAN_CELTOT
    REAL(q)                 :: FUNC_FRAC
    ! local
    INTEGER    :: IU0
    COMPLEX(q) :: RZ(1),RH(1)

    RZ(1)=CMPLX(X,0.0_q,q)
    CALL PADE_FRAC_EVAL(N,IZ,A,1,RZ(1),RH(1),IU0,WMEAN_CELTOT)
    FUNC_FRAC=REAL(RH(1),q)-X
END FUNCTION FUNC_FRAC


!*********************************************************************
!
! FUNC_SVD=SIGMA(w)-w for SVD version 
!
!*********************************************************************
FUNCTION FUNC_SVD(N,IZ,A,X,WMEAN_CELTOT)
    IMPLICIT NONE
    INTEGER,   INTENT(IN)   :: N        ! number of imaginary data sets (iz(i),ih(i),i=1,...N))
    COMPLEX(q),INTENT(IN)   :: IZ(N)    ! imaginary frequencies iz(i)  !not used
    COMPLEX(q),INTENT(IN)   :: A(N)     ! coefficients of pade
    REAL(q),   INTENT(IN)   :: X
    COMPLEX(q)              :: WMEAN_CELTOT
    REAL(q)                 :: FUNC_SVD
    ! local
    INTEGER    :: IU0
    COMPLEX(q) :: RZ(1),RH(1)

    RZ(1)=CMPLX(X,0.0_q,q)
    CALL PADE_SVD_EVAL(N,A,1,RZ,RH,IU0,WMEAN_CELTOT)
    FUNC_SVD=REAL(RH(1),q)-X
END FUNCTION FUNC_SVD


!**********************************************************************
!
! this subroutine is implemented by following the idea from             
! G. A. Baker Jr., Essentials of Pade Approximants                      
! (Academic Press, New York, 1975),Chapter 8                            
! See also related works in                                             
! K. Lee and K. Chang, Phys. Rev. B 54, R8285 (1996)                    
! H. J. Vidberg and J. W. Serene, J, Low, Temp. Phys. 29, 179 (1977)    
!
!**********************************************************************
SUBROUTINE ACONT_PADE_FRAC(N,IZ,IH, M,RZ,RH,IU0,WMEAN_CELTOT)
    USE tutor, ONLY: vtutor, isAlert, GWkweight
    IMPLICIT NONE
    INTEGER,   INTENT(IN)   :: N        ! number of imaginary data sets (iz(i),ih(i),i=1,...N))
    COMPLEX(q),INTENT(IN)   :: IZ(N)    ! imaginary frequencies iz(i)
    COMPLEX(q),INTENT(IN)   :: IH(N)    ! function evaluated at imaginary frequencies iz(i)
    INTEGER,   INTENT(IN)   :: M        ! number of real data sets (rz(i),rh(i),i=1,...M))
    COMPLEX(q),INTENT(IN)   :: RZ(M)    ! real frequencies for which analytic continuation is to be performed
    COMPLEX(q),INTENT(OUT)  :: RH(M)    ! function evaluated at real frequencies rz(i)
    INTEGER,   INTENT(IN)   :: IU0
    COMPLEX(q),OPTIONAL     :: WMEAN_CELTOT  ! mean-field HF eigenvalue
    ! local
    COMPLEX(q), ALLOCATABLE :: ATMP(:), BTMP(:), A(:), G(:, :)
    COMPLEX(q)              :: CTMP
    INTEGER                 :: I,J

    IF (N<2) THEN
       CALL vtutor%bug("internal ERROR in VASP: N should be >=2 in ACONT_PADE_FIT(N,IZ,IH,M,RZ,RH,"&
          // "IU0) Stop!", __FILE__, __LINE__)
    END IF

    ALLOCATE(ATMP(N), BTMP(0:N), A(N), G(N,N))

    A(1)=IH(1)     ! A(1)=G(1,1)=F(1)  
    G(1,:)=IH(:)   ! (G_1(Z(i))=F(i), i=1,N) =>  G(1,1:N)=F(1:N) 
    DO I=2,N
       DO J=2,N
          CTMP=(IZ(J)-IZ(I-1))*G(I-1,J)
          IF ( ABS(CTMP) < 1.0E-6  ) THEN
             G(I,J)=0.0_q
          ELSE
             G(I,J) = ( A(I-1)-G(I-1,J) ) / CTMP
          ENDIF 
       ENDDO
       A(I)=G(I,I)   
    ENDDO

    DO I=1,M
       CTMP=RZ(I)
       BTMP(:)=CMPLX(1.0_q, 0.0_q,q)
       BTMP(0)=CMPLX(0.0_q, 0.0_q,q)
       ATMP(1)=A(1)
       DO J=2,N
          ATMP(J)=A(J)*(CTMP-IZ(J-1)) 
       ENDDO
       CALL PADE_EVALUATE(N, ATMP, BTMP, RH(I))
    ENDDO

    ! add mean-field HF eigenvalue
    IF (PRESENT(WMEAN_CELTOT)) THEN
      RH=RH+WMEAN_CELTOT
    ENDIF

    DEALLOCATE(ATMP, BTMP, A, G)
END SUBROUTINE ACONT_PADE_FRAC


!*********************************************************************
!
! obtain the coefficients A(:) of pade fit for FRAC version
!
!*********************************************************************
SUBROUTINE PADE_FRAC_COEFF(N,IZ,IH,A,IU0)
    USE tutor, ONLY: vtutor, isAlert, GWkweight
    IMPLICIT NONE
    INTEGER,   INTENT(IN)   :: N        ! number of imaginary data sets (iz(i),ih(i),i=1,...N))
    COMPLEX(q),INTENT(IN)   :: IZ(N)    ! imaginary frequencies iz(i)
    COMPLEX(q),INTENT(IN)   :: IH(N)    ! function evaluated at imaginary frequencies iz(i)
    COMPLEX(q),INTENT(OUT)  :: A(N)     ! coefficients of pade fit
    INTEGER,   INTENT(IN)   :: IU0
    ! local
    COMPLEX(q), ALLOCATABLE :: G(:,:)
    COMPLEX(q)              :: CTMP
    INTEGER                 :: I,J

    IF (N<2) THEN
       CALL vtutor%bug("internal ERROR in VASP: N should be >=2 in PADE_FIT_COEFF(N,IZ,IH,A,IU0) Stop!", &
          __FILE__, __LINE__)
    END IF
    ALLOCATE(G(N,N))

    A(1)=IH(1)     ! A(1)=G(1,1)=F(1)  
    G(1,:)=IH(:)   ! (G_1(Z(i))=F(i), i=1,N) =>  G(1,1:N)=F(1:N) 
    DO I=2,N
       DO J=2,N
          CTMP=(IZ(J)-IZ(I-1))*G(I-1,J)
          IF ( ABS(CTMP) < 1.0E-8  ) THEN
             G(I,J)=0.0_q
          ELSE
             G(I,J) = ( A(I-1)-G(I-1,J) ) / CTMP
          ENDIF 
       ENDDO
       A(I)=G(I,I)   
    ENDDO

    DEALLOCATE(G)
END SUBROUTINE PADE_FRAC_COEFF


!*********************************************************************
!
! Evaluation of the pade_fit at specified real frequency
! PADE_FRAC_COEFF should be called before you call PADE_FRAC_EVAL
! note: for interpolation: IZ=(EF+iw) to RZ=(EF+iw')
!       for analytic continuation, IZ=(EF+iw) to RZ=w'
!
!*********************************************************************
SUBROUTINE PADE_FRAC_EVAL(N,IZ,A,M,RZ,RH,IU0,WMEAN_CELTOT)
    IMPLICIT NONE
    INTEGER,   INTENT(IN)   :: N        ! number of imaginary data sets (iz(i),ih(i),i=1,...N))
    COMPLEX(q),INTENT(IN)   :: IZ(N)    ! imaginary frequencies iz(i)
    COMPLEX(q),INTENT(IN)   :: A(N)     ! coefficients of pade
    INTEGER,   INTENT(IN)   :: M        ! number of real data sets (rz(i),rh(i),i=1,...M))
    COMPLEX(q),INTENT(IN)   :: RZ(M)    ! real frequencies for which analytic continuation is to be performed
    COMPLEX(q),INTENT(OUT)  :: RH(M)    ! function evaluated at real frequencies rz(i)
    INTEGER,   INTENT(IN)   :: IU0
    COMPLEX(q),OPTIONAL     :: WMEAN_CELTOT  ! mean-field HF eigenvalue
    ! local
    COMPLEX(q), ALLOCATABLE :: ATMP(:), BTMP(:)
    COMPLEX(q)              :: CTMP
    INTEGER                 :: I,J

    ALLOCATE(ATMP(N), BTMP(0:N))

    DO I=1,M
       CTMP=RZ(I)
       BTMP(:)=CMPLX(1.0_q, 0.0_q,q)
       BTMP(0)=CMPLX(0.0_q, 0.0_q,q)
       ATMP(1)=A(1)
       DO J=2,N
          ATMP(J)=A(J)*(CTMP-IZ(J-1)) 
       ENDDO
       CALL PADE_EVALUATE(N, ATMP, BTMP, RH(I))
    ENDDO

    ! add mean-field HF eigenvalue
    IF (PRESENT(WMEAN_CELTOT)) THEN
      RH=RH+WMEAN_CELTOT
    ENDIF

    DEALLOCATE(ATMP, BTMP)
END SUBROUTINE PADE_FRAC_EVAL


!*********************************************************************
SUBROUTINE PADE_EVALUATE(N, A, B, FUN)
    IMPLICIT NONE
    INTEGER,   INTENT(IN)  :: N
    COMPLEX(q),INTENT(IN)  :: A(N)
    COMPLEX(q),INTENT(IN)  :: B(0:N)
    COMPLEX(q),INTENT(OUT) :: FUN
    ! local
    INTEGER                :: I
    FUN=B(N)
    DO I=N,1,-1
       FUN=B(I-1)+A(I)/FUN
    ENDDO
END SUBROUTINE PADE_EVALUATE


!*********************************************************************
!
! obtain the coefficients of pade fit  AX=B
! in this routine, it requires NOMEGA >= NOMEGA_SVD
! A(NOMEGA,NOMEGA_SVD) * X(NOMEGA_SVD) = B(NOMEGA)
! NOTE: double precision by calling LAPACK
!
!*********************************************************************
SUBROUTINE PADE_SVD_COEFF_q(NOMEGA,IZ,IH,NOMEGA_SVD,X,IU0)
    USE tutor, ONLY: vtutor, isAlert, GWkweight
    IMPLICIT NONE
    INTEGER,   INTENT(IN)   :: NOMEGA              ! number of imaginary data sets 
    COMPLEX(q),INTENT(IN)   :: IZ(NOMEGA)          ! imaginary frequencies iz(i)
    COMPLEX(q),INTENT(IN)   :: IH(NOMEGA)          ! function evaluated at iz(i)
    INTEGER,   INTENT(IN)   :: NOMEGA_SVD          ! number of SVD coefficients
    COMPLEX(q),INTENT(OUT)  :: X(NOMEGA_SVD)       ! solution of AX=B
    INTEGER,   INTENT(IN)   :: IU0
    ! local
    COMPLEX(q), ALLOCATABLE :: A(:,:),B(:)
    REAL(q)                 :: WEIGHT
    INTEGER                 :: I,J,N

    ! quick check 
    IF (MOD(NOMEGA_SVD,2) /=0 ) THEN
      CALL vtutor%bug("internal error in PADE_SVD_COEFF"//str( NOMEGA_SVD), __FILE__, __LINE__)
    ENDIF
    
    ALLOCATE(A(NOMEGA,NOMEGA_SVD))     ! A matrix
    ALLOCATE(B(NOMEGA))                ! B verctor
    A=(0.0_q,0.0_q)
    B=(0.0_q,0.0_q)

    N=NOMEGA_SVD/2
    DO I=1,NOMEGA
       DO J=1,NOMEGA_SVD,2
         A(I,J)=IZ(NOMEGA-I+1)**(J/2)
         A(I,J+1)=-A(I,J)*IH(NOMEGA-I+1)
       ENDDO
       WEIGHT=0.0_q
       DO J=1,NOMEGA_SVD
          WEIGHT=WEIGHT+ABS(A(I,J))
       ENDDO
       WEIGHT=WEIGHT/REAL(NOMEGA_SVD,q)
       ! determine A matrix    
       A(I,:)=A(I,:)/WEIGHT
       ! determine B vector
       B(I)=IH(NOMEGA-I+1)*(IZ(NOMEGA-I+1)**N)/WEIGHT
    ENDDO

    ! solve AX=B by QR in double precision by calling LAPACK
    CALL SOLVE_LSQ_QR(NOMEGA,NOMEGA_SVD,A,B,X,IU0)

    DEALLOCATE(A,B)
END SUBROUTINE PADE_SVD_COEFF_q


!*********************************************************************
!
! obtain the coefficients of pade fit  AX=B
! in this routine, it requires NOMEGA >= NOMEGA_SVD
! A(NOMEGA,NOMEGA_SVD) * X(NOMEGA_SVD) = B(NOMEGA)
! NOTE: quadruple precision
!
!*********************************************************************
SUBROUTINE PADE_SVD_COEFF(NOMEGA,IZ,IH,NOMEGA_SVD,X,IU0)
#ifdef qd_emulate
   USE qdmodule
#endif
    USE base, ONLY: TOREAL
    USE tutor, ONLY: vtutor, isAlert, GWkweight
    IMPLICIT NONE
    INTEGER,   INTENT(IN)   :: NOMEGA              ! number of imaginary data sets 
    COMPLEX(q),INTENT(IN)   :: IZ(NOMEGA)          ! imaginary frequencies iz(i)
    COMPLEX(q),INTENT(IN)   :: IH(NOMEGA)          ! function evaluated at iz(i)
    INTEGER,   INTENT(IN)   :: NOMEGA_SVD          ! number of SVD coefficients
    COMPLEX(q),INTENT(OUT)  :: X(NOMEGA_SVD)       ! solution of AX=B
    INTEGER,   INTENT(IN)   :: IU0
    ! local
    QDPCOMPLEX, ALLOCATABLE :: A(:,:),B(:),XX(:)
    COMPLEX(q)              :: CD
    QDPCOMPLEX              :: CQD
    QDPREAL                 :: R1,R2,R3
    QDPREAL                 :: WEIGHT
    INTEGER                 :: I,J,N

    ! quick check 
    IF (MOD(NOMEGA_SVD,2) /=0 ) THEN
      CALL vtutor%bug("internal error in PADE_SVD_COEFF"//str( NOMEGA_SVD), __FILE__, __LINE__)
    ENDIF
    
    ALLOCATE(A(NOMEGA,NOMEGA_SVD))     ! A matrix
    ALLOCATE(B(NOMEGA))                ! B verctor
    ALLOCATE(XX(NOMEGA_SVD))           ! X verctor
    A=0
    B=0
    XX=0

    N=NOMEGA_SVD/2
    DO I=1,NOMEGA
       DO J=1,NOMEGA_SVD,2
         A(I,J)=IZ(NOMEGA-I+1)**(J/2)
         R1 = REAL(IH(NOMEGA-I+1),q)
         R2 = AIMAG(IH(NOMEGA-I+1) )
         !A(I,J+1)=-A(I,J)*QDPC( REAL(IH(NOMEGA-I+1),q), AIMAG( IH(NOMEGA-I+1) ))
         A(I,J+1)=-A(I,J)*QDPC( R1, R2 )
       ENDDO
       WEIGHT=0
       DO J=1,NOMEGA_SVD
          WEIGHT=WEIGHT+ABS(A(I,J))
       ENDDO
       WEIGHT=WEIGHT/REAL(NOMEGA_SVD,q)
       ! determine A matrix    
       A(I,:)=A(I,:)/WEIGHT
       ! determine B vector
       !B(I)=IH(NOMEGA-I+1)*(IZ(NOMEGA-I+1)**N)/WEIGHT
       CD = IH(NOMEGA-I+1)*(IZ(NOMEGA-I+1)**N)
       R1 = REAL( CD , q)
       R2 = AIMAG( CD )
       B(I)= QDPC( R1, R2 )/WEIGHT
    ENDDO

    ! solve AX=B by QR in quadruple precision
    CALL QQRSOLVE(A,NOMEGA,NOMEGA_SVD,B,XX)
    DO J=1,NOMEGA_SVD
       CD = XX(J)
       R1 = REAL( CD , q)
       R2 = AIMAG( CD )
       X(J)=CMPLX(TOREAL(R1),TOREAL(R2),q)
    ENDDO

    DEALLOCATE(A,B,XX)
END SUBROUTINE PADE_SVD_COEFF
!*********************************************************************


!*********************************************************************
!
! evaluate the Pade polynminal at RZ
! PADE_SVD_COEFF should be called before you call PADE_SVD_EVAL
! note: for interpolation: IZ=(EF+iw) to RZ=(EF+iw')
!       for analytic continuation, IZ=(EF+iw) to RZ=w'
!
!*********************************************************************
SUBROUTINE PADE_SVD_EVAL(NOMEGA_SVD, X, M,RZ, RH, IU0, WMEAN_CELTOT)
    USE tutor, ONLY: vtutor, isAlert, GWkweight
    IMPLICIT NONE
    INTEGER,   INTENT(IN)  :: NOMEGA_SVD
    COMPLEX(q),INTENT(IN)  :: X(NOMEGA_SVD) ! SVD coefficients
    INTEGER,   INTENT(IN)  :: M             ! number of real data sets (rz(i),rh(i),i=1,...M))
    COMPLEX(q),INTENT(IN)  :: RZ(M)         ! real frequencies for which analytic continuation is performed
    COMPLEX(q),INTENT(OUT) :: RH(M)         ! function evaluated at real frequencies rz(i)
    INTEGER,   INTENT(IN)  :: IU0
    COMPLEX(q),OPTIONAL    :: WMEAN_CELTOT  ! mean-field HF eigenvalue
    ! local
    INTEGER                :: I,J,N
    COMPLEX(q)             :: UP,DN
    COMPLEX(q),ALLOCATABLE :: P(:),QQ(:)

    IF (MOD(NOMEGA_SVD,2) /= 0) THEN
       CALL vtutor%error( "NOMEGA_SVD in PADE_SVD_EVAL should be an even number "//str(NOMEGA_SVD) )
    ENDIF

    IF (SIZE(RZ) /= SIZE(RH) ) THEN
       CALL vtutor%error(" SIZE(RZ) /= SIZE(RH) "//str(SIZE(RZ))//" "//str( SIZE(RH) ))
    ENDIF    

    N=NOMEGA_SVD/2

    ALLOCATE(P(N),QQ(N))
    P=(0.0,0.0)
    QQ=(0.0,0.0)
    P=X(1:NOMEGA_SVD:2)
    QQ=X(2:NOMEGA_SVD:2)

    DO I=1,SIZE(RZ)
      UP=(0.0,0.0)
      DN=(0.0,0.0)
      DN=DN+RZ(I)**N
      DO J=1,N
        UP=UP+P(J)*(RZ(I)**(J-1))
        DN=DN+QQ(J)*(RZ(I)**(J-1))
      ENDDO
      RH(I)=UP/DN
    ENDDO

    ! add mean-field HF eigenvalue
    IF (PRESENT(WMEAN_CELTOT)) THEN
      RH=RH+WMEAN_CELTOT
    ENDIF

    DEALLOCATE(P,QQ)
END SUBROUTINE PADE_SVD_EVAL


!*********************************************************************
!
! Solve AX=B, with A(M,N), X(N) and B(M) using QR by LAPACK
! NOTE that M>=N and this routine does not modify A and B
!
!*********************************************************************
SUBROUTINE SOLVE_LSQ_QR(M,N,A,B,X,IU0)
    USE tutor, ONLY: vtutor, isAlert, GWkweight
    IMPLICIT NONE
    INTEGER,    INTENT(IN)  :: M,N
    COMPLEX(q), INTENT(IN)  :: A(M,N)
    COMPLEX(q), INTENT(IN)  :: B(M)
    COMPLEX(q), INTENT(OUT) :: X(N)
    INTEGER, INTENT(IN)     :: IU0
    !LOCAL
    INTEGER                 :: LDA,LDB,LWORK,INFO
    COMPLEX(q),ALLOCATABLE  :: WORK(:),A2(:,:),B2(:)

    ! quick return
    IF (M<N) THEN
      CALL vtutor%bug('internal error in VASP: SOLVE_LSQ_QR: M <N in SOLVE_LSQ_QR '//str(M)//" "//str(N), __FILE__, __LINE__)
    ENDIF

    ! back up A of B, this routine do not touch A and B
    ALLOCATE(A2(M,N))
    A2=A
    ALLOCATE(B2(M))
    B2=B

    ! set up input parameters 
    LDA=M
    LDB=MAX(M,N)

    ! obtain optimal working array sizes
    ALLOCATE(WORK(1))
    LWORK=-1
    CALL ZGELS('N',M,N,1,A2,LDA,B2,LDB,WORK,LWORK,INFO)
    LWORK=INT(WORK(1))
    IF ( LWORK <= 0 )  THEN
       CALL vtutor%bug('internal error in VASP: SOLVE_LSQ_QR: illegal size from ZGELS '//str(LWORK), __FILE__, __LINE__)
    ENDIF
    ! allocate working array with optimal size 
    DEALLOCATE(WORK)
    ALLOCATE(WORK(LWORK))

    ! actural LSQ QR sovle AX=B comes now 
    CALL ZGELS('N',M,N,1,A2,LDA,B2,LDB,WORK,LWORK,INFO)
    IF ( INFO/=0) THEN
       CALL vtutor%bug('internal error in VASP: SOLVE_LSQ_QR: ZGELS reports'//str(INFO), __FILE__, __LINE__)
    ENDIF

    X(1:N)=B2(1:N)
 
    DEALLOCATE(A2,B2)
    DEALLOCATE(WORK)

END SUBROUTINE SOLVE_LSQ_QR
!*********************************************************************


!*********************************************************************
!
! obtain the coefficients of pade fit  AX=B
! A(NOMEGA,NOMEGA) * X(NOMEGA) = B(NOMEGA)
!
!*********************************************************************
SUBROUTINE PADE_LINEAR_COEFF(NOMEGA,IZ,IH,X,IU0)
    USE tutor, ONLY: vtutor, isAlert, GWkweight
    IMPLICIT NONE
    INTEGER,   INTENT(IN)   :: NOMEGA         ! number of imaginary data sets 
    COMPLEX(q),INTENT(IN)   :: IZ(NOMEGA)     ! imaginary frequencies iz(i)
    COMPLEX(q),INTENT(IN)   :: IH(NOMEGA)     ! function evaluated at iz(i)
    COMPLEX(q),INTENT(OUT)  :: X(NOMEGA)      ! solution of AX=B
    INTEGER,   INTENT(IN)   :: IU0
    ! local
    COMPLEX(q), ALLOCATABLE :: A(:,:),B(:)
    REAL(q)                 :: WEIGHT
    INTEGER                 :: I,J,N

    ! quick check 
    IF (MOD(NOMEGA,2) /=0 ) THEN
      CALL vtutor%bug( "internal error in PADE_LINEAR_COEFF "//str(NOMEGA), __FILE__, __LINE__)
    ENDIF
    
    ALLOCATE(A(NOMEGA,NOMEGA))     ! A matrix
    ALLOCATE(B(NOMEGA))            ! B verctor
    A=(0.0_q,0.0_q)
    B=(0.0_q,0.0_q)

    ! determine A matrix    
    N=NOMEGA/2
    DO I=1,NOMEGA
       DO J=1,NOMEGA,2
         A(I,J)=IZ(NOMEGA-I+1)**(J/2)
         A(I,J+1)=-A(I,J)*IH(NOMEGA-I+1)
       ENDDO
       WEIGHT=0.0_q
       DO J=1,NOMEGA
          WEIGHT=WEIGHT+ABS(A(I,J))
       ENDDO
       WEIGHT=WEIGHT/REAL(NOMEGA,q)
       ! determine A matrix    
       A(I,:)=A(I,:)/WEIGHT
       ! determine B vector
       B(I)=IH(NOMEGA-I+1)*(IZ(NOMEGA-I+1)**N)/WEIGHT
    ENDDO

    ! second, solve AX=B by LU decopisition
    CALL LINEAR_SOLVE('N',NOMEGA,A,B,X,IU0)

    DEALLOCATE(A,B)
END SUBROUTINE PADE_LINEAR_COEFF
!*********************************************************************


!*********************************************************************
!
! solve AX=B by Lapack iteratively
!   A(N,N) = (in) general matrix A         
!   B(N) = (in) the right hand side matrix B  
!   X(N) = (out) the solution matrix X(N)    
!   (NOTE: this routine does not modify A and B)
! 
!*********************************************************************
SUBROUTINE LINEAR_SOLVE(TRANS,N,A,B,X,IU0)
    USE tutor, ONLY: vtutor, isAlert, GWkweight
    IMPLICIT NONE
    CHARACTER(LEN=*)           :: TRANS          ! 'N', 'T', 'C'
    INTEGER, INTENT(in)        :: N
    COMPLEX(q), INTENT(in)     :: A(N,N)
    COMPLEX(q), INTENT(in)     :: B(N)
    COMPLEX(q), INTENT(out)    :: X(N)
    INTEGER, INTENT(in)        :: IU0
    ! local
    INTEGER, PARAMETER         :: MAXITER=4
    INTEGER                    :: LDA          ! The leading dimension of the array A
    INTEGER                    :: LDB          ! The leading dimension of the array B
    INTEGER                    :: INFO         !  = 0:  successful exit
    INTEGER                    :: I
    INTEGER    ,ALLOCATABLE    :: IPIV(:)      ! The pivot indices from DGETRF
    COMPLEX(q) ,ALLOCATABLE    :: ATMP(:,:),DX(:)

    ALLOCATE( ATMP(N,N),DX(N),IPIV(N) )
    LDA=N
    LDB=N

    ! first LU decomposition of A=LU
    ATMP=A
    CALL ZGETRF(N,N,ATMP,LDA,IPIV,INFO)
    
    ! solve AX=B
    X=B
    IF (INFO ==0 ) THEN
       CALL ZGETRS(TRANS, N, 1, ATMP, LDA, IPIV, X, LDB, INFO)
    ENDIF

    ! check for convergence.
    IF(INFO /= 0 .AND. IU0>=0 ) THEN
       CALL vtutor%bug( "internal error in VASP: LINEAR_SOLVE did not work properly when calling ZGETRS " &
          //str(INFO), __FILE__, __LINE__)
    ENDIF

    ! check how accurate solution is and iterate 
    DO I=1,MAXITER
      !calculate DX=AX-B
      DX=B
      CALL ZGEMV('N',N,N,(1.0_q,0.0_q),A,N,X,1,(-1.0_q,0.0_q),DX,1)
      ! IF(IU0>=0 ) WRITE(IU0,*) "ERROR= ",REAL(DOT_PRODUCT(DX,DX),q)
      ! solve again with LU decoposited ATMP and same IPIV
      CALL ZGETRS(TRANS, N, 1, ATMP, LDA, IPIV, DX, LDB, INFO)
      X=X-DX
    ENDDO

    DEALLOCATE(ATMP,DX,IPIV)
END SUBROUTINE LINEAR_SOLVE
!*********************************************************************


!*********************************************************************
!
! Pade SVD verstion, combine like ACONT_PADE_FRAC
! 
!*********************************************************************
SUBROUTINE ACONT_PADE_SVD(NOMEGA,IZ,IH,M,RZ,RH,IU0,NOMEGA_SVD,WMEAN_CELTOT)
    IMPLICIT NONE
    INTEGER,   INTENT(IN)   :: NOMEGA        ! number of imaginary data sets 
    COMPLEX(q),INTENT(IN)   :: IZ(NOMEGA)    ! imaginary frequencies iz(i)
    COMPLEX(q),INTENT(IN)   :: IH(NOMEGA)    ! function evaluated at imaginary frequencies iz(i)
    INTEGER,   INTENT(IN)   :: M             ! number of real data sets (rz(i),rh(i),i=1,...M))
    COMPLEX(q),INTENT(IN)   :: RZ(M)         ! real frequencies for which analytic continuation is done
    COMPLEX(q),INTENT(OUT)  :: RH(M)         ! function evaluated at real frequencies rz(i)
    INTEGER,   INTENT(IN)   :: IU0
    INTEGER,   INTENT(IN)   :: NOMEGA_SVD
    COMPLEX(q),OPTIONAL     :: WMEAN_CELTOT  ! mean-field HF eigenvalue
    ! local
    COMPLEX(q), ALLOCATABLE :: X(:)

    ALLOCATE(X(NOMEGA_SVD)) ; X=(0.0_q,0.0_q)

    CALL PADE_SVD_COEFF(NOMEGA,IZ,IH,NOMEGA_SVD,X,IU0)

    IF (PRESENT(WMEAN_CELTOT)) THEN 
       CALL PADE_SVD_EVAL(NOMEGA_SVD,X,M,RZ,RH,IU0,WMEAN_CELTOT) 
    ELSE
       CALL PADE_SVD_EVAL(NOMEGA_SVD,X,M,RZ,RH,IU0) 
    ENDIF

    DEALLOCATE(X) 
END SUBROUTINE ACONT_PADE_SVD
!*********************************************************************


!*********************************************************************
!
! linear Pade verstion, combine like ACONT_PADE_FRAC
! 
!*********************************************************************
SUBROUTINE ACONT_PADE_LINEAR(NOMEGA,IZ,IH,M,RZ,RH,IU0,WMEAN_CELTOT)
    USE tutor, ONLY: vtutor, isAlert, GWkweight
    IMPLICIT NONE
    INTEGER,   INTENT(IN)   :: NOMEGA        ! number of imaginary data sets (iz(i),ih(i),i=1,...N))
    COMPLEX(q),INTENT(IN)   :: IZ(NOMEGA)    ! imaginary frequencies iz(i)
    COMPLEX(q),INTENT(IN)   :: IH(NOMEGA)    ! function evaluated at imaginary frequencies iz(i)
    INTEGER,   INTENT(IN)   :: M             ! number of real data sets 
    COMPLEX(q),INTENT(IN)   :: RZ(M)         ! real frequencies for which analytic continuation is done
    COMPLEX(q),INTENT(OUT)  :: RH(M)         ! function evaluated at real frequencies rz(i)
    INTEGER,   INTENT(IN)   :: IU0
    COMPLEX(q),OPTIONAL     :: WMEAN_CELTOT  ! mean-field HF eigenvalue
    ! local
    COMPLEX(q), ALLOCATABLE :: A(:,:), B(:), X(:)
    REAL(q)                 :: WEIGHT
    INTEGER                 :: I,J,N

    IF (MOD(NOMEGA,2) /= 0) THEN
       CALL vtutor%error("NOMEGA in ACONT_PADE_LINEAR should be an even number: "//str(NOMEGA))
    ENDIF

    ALLOCATE(A(NOMEGA,NOMEGA))
    ALLOCATE(B(NOMEGA))
    ALLOCATE(X(NOMEGA))
    A=(0.0_q,0.0_q)
    B=(0.0_q,0.0_q)
    X=(0.0_q,0.0_q)
    
    ! determine A matrix    
    N=NOMEGA/2
    DO I=1,NOMEGA
       DO J=1,NOMEGA,2
         A(I,J)=IZ(NOMEGA-I+1)**(J/2)
         A(I,J+1)=-A(I,J)*IH(NOMEGA-I+1)
       ENDDO
       WEIGHT=0.0_q
       DO J=1,NOMEGA
          WEIGHT=WEIGHT+ABS(A(I,J))
       ENDDO
       WEIGHT=WEIGHT/REAL(NOMEGA,q)
       ! determine A matrix    
       A(I,:)=A(I,:)/WEIGHT
       ! determine B vector
       B(I)=IH(NOMEGA-I+1)*(IZ(NOMEGA-I+1)**N)/WEIGHT
    ENDDO

    CALL LINEAR_SOLVE('N',NOMEGA,A,B,X,IU0)

    IF (PRESENT(WMEAN_CELTOT)) THEN 
       CALL PADE_SVD_EVAL(NOMEGA,X,M,RZ,RH,IU0,WMEAN_CELTOT) 
    ELSE
       CALL PADE_SVD_EVAL(NOMEGA,X,M,RZ,RH,IU0) 
    ENDIF

    DEALLOCATE(A,B,X) 
END SUBROUTINE ACONT_PADE_LINEAR
!*********************************************************************


!*********************************************************************
!
! the integrated analytic continuation by Pade fit
! 
!*********************************************************************
SUBROUTINE ACONT_PADE_FIT(NOMEGA,IZ,IH,M,RZ,RH,IU0,NOMEGA_SVD,OPTION,WMEAN_CELTOT)
    IMPLICIT NONE
    INTEGER,   INTENT(IN)   :: NOMEGA        ! number of imaginary data sets 
    COMPLEX(q),INTENT(IN)   :: IZ(NOMEGA)    ! imaginary frequencies iz(i)
    COMPLEX(q),INTENT(IN)   :: IH(NOMEGA)    ! function evaluated at imaginary frequencies iz(i)
    INTEGER,   INTENT(IN)   :: M             ! number of real data sets (rz(i),rh(i),i=1,...M))
    COMPLEX(q),INTENT(IN)   :: RZ(M)         ! real frequencies for which analytic continuation is done
    COMPLEX(q),INTENT(OUT)  :: RH(M)         ! function evaluated at real frequencies rz(i)
    INTEGER,   INTENT(IN)   :: IU0
    INTEGER,   INTENT(IN)   :: NOMEGA_SVD
    INTEGER,   INTENT(IN)   :: OPTION
    COMPLEX(q),OPTIONAL     :: WMEAN_CELTOT  ! mean-field HF eigenvalue

    IF (PRESENT(WMEAN_CELTOT)) THEN
       SELECT CASE(OPTION)
         CASE(1)
            CALL ACONT_PADE_FRAC(NOMEGA,IZ,IH,M,RZ,RH,IU0,WMEAN_CELTOT)
         CASE(2)
            CALL ACONT_PADE_SVD(NOMEGA,IZ,IH,M,RZ,RH,IU0,NOMEGA_SVD,WMEAN_CELTOT)
         CASE(3)
            CALL ACONT_PADE_LINEAR(NOMEGA,IZ,IH,M,RZ,RH,IU0,WMEAN_CELTOT)
       END SELECT
    ELSE
       SELECT CASE(OPTION)
         CASE(1)
            CALL ACONT_PADE_FRAC(NOMEGA,IZ,IH,M,RZ,RH,IU0)
         CASE(2)
            CALL ACONT_PADE_SVD(NOMEGA,IZ,IH,M,RZ,RH,IU0,NOMEGA_SVD)
         CASE(3)
            CALL ACONT_PADE_LINEAR(NOMEGA,IZ,IH,M,RZ,RH,IU0)
       END SELECT
    ENDIF

END SUBROUTINE ACONT_PADE_FIT
!*********************************************************************


!*********************************************************************
!
! calculate the QR decomposition of A in quadruple precision
!
!*********************************************************************
SUBROUTINE QQRDC(A, M, N, QQ, R, P, PIVOTING)
    USE minimax_functions1D, ONLY: SIGN_QD
    IMPLICIT NONE
    QDPCOMPLEX    :: A(:,:)
    INTEGER              :: M, N
    QDPCOMPLEX    :: QQ(:,:), R(:,:)
    INTEGER, ALLOCATABLE :: P(:)
    LOGICAL              :: PIVOTING
    ! local
    QDPCOMPLEX, ALLOCATABLE :: H(:,:)
    QDPCOMPLEX, ALLOCATABLE :: V(:,:)
    QDPCOMPLEX, ALLOCATABLE :: Y(:)
    QDPREAL, ALLOCATABLE    :: NORMS(:)
    QDPCOMPLEX              :: ALPHA
    INTEGER           :: I, J, K
    INTEGER           :: PIVCOL
    QDPREAL           :: TMP, R1, R2
    QDPCOMPLEX        :: TMPC

    IF (PIVOTING) THEN
      P(1:N) = (/ (I, I=1,N) /)
    END IF

    ALLOCATE(H(M,M))
    ALLOCATE(V(M,M))
    ALLOCATE(Y(M))
    ALLOCATE(NORMS(N))

    R = A
    QQ = 0
    DO I=1,M
       QQ(I,I) = 1
    ENDDO

    DO I=1,N
      IF (PIVOTING) THEN
        ! determine pivot column
        NORMS = 0
        PIVCOL = I
        DO J=I,N
          !NORMS(J) = SUM(R(:,J) * CONJG(R(:,J)))
          NORMS(J) = 0 
          DO K = 1, SIZE( R,1) 
             TMP = QDPC2QDPREAL( R( K, J )*CONJG( R(K,J) ) )
             NORMS(J) = NORMS(J) + TMP
          ENDDO

          IF (NORMS(J) > NORMS(PIVCOL)) THEN
            PIVCOL = J
          END IF
        END DO

        ! swap columns
        IF (PIVCOL /= I) THEN
          R(:,I:PIVCOL:PIVCOL-I) = R(:,PIVCOL:I:I-PIVCOL)
          P(I:PIVCOL:PIVCOL-I) = P(PIVCOL:I:I-PIVCOL)
        END IF
      END IF

      Y = 0
      Y(I:M) = R(I:M,I)

      ! the sign of alpha is in principle arbitrary; if Y(I) is real, then
      ! usually -SIGN(Y(I)) is used, in order to avoid loss of precision in
      ! the substraction later on; in our case however, Y(I) is complex, so
      ! it has no defined sign; I therefore use -SIGN(REAL(Y(I))) here
      !ALPHA = SQRT(SUM(Y*Y)) * (-QDPP(1)) * SIGN_QD( QDPC2QDPREAL( Y(I) ) )
      TMP = 0 
      DO K = 1, SIZE( Y ) 
         TMP = TMP + ABS( QDPC2QDPREAL( Y(K)*Y(K) ) )
      ENDDO
      R1 = -SIGN_QD( QDPC2QDPREAL( Y(I) ) )
      R2 = 0
      ALPHA = QDPC( SQRT( TMP ),QDPP(0) ) * QDPC( R1, R2 )

      ! compute the reflection vector, store it in the leftmost column
      ! of the matrix V; this in necessary because we will later need a
      ! matrix of the form MATMUL_QDC(V, TRANSPOSE_QDC(V))
      V = 0
      V(:,I) = Y
      V(I,I) = V(I,I) - ALPHA

      ! compute the Householder matrix
      H = 0
      DO J=1,M
         H(J,J) = 1
      ENDDO
      TMPC = 0 
      DO K = 1, SIZE( V,1 ) 
         TMPC=TMPC+V(K,I)*V(K,I)
      ENDDO
      !H = H - QDPP(2)*MATMUL_QDC(V, TRANSPOSE_QDC(V)) / SUM(V(:,I) * V(:,I))
      H = H - QDPP(2)*MATMUL_QDC(V, TRANSPOSE_QDC(V)) / TMPC

      R = MATMUL_QDC(H, R)
      QQ = MATMUL_QDC(QQ, H)
    END DO

    ! invert the permutation
    P(P) = (/ (I, I=1,N) /)

    DEALLOCATE(H,V,Y,NORMS)
END SUBROUTINE QQRDC
!*********************************************************************


!*********************************************************************
!
! solve an upper triangular system of equations, using backwards
! substitution, in quadruple precision
! writen by Manuel Grumet
!
!*********************************************************************
SUBROUTINE QBSSOLVE(R, N, B, X)
    IMPLICIT NONE
    QDPCOMPLEX :: R(:,:)
    QDPCOMPLEX :: B(:), X(:)
    INTEGER           :: N,I
    INTEGER           :: K
    QDPCOMPLEX :: TMPC

    DO I=N,1,-1
       TMPC=0
       DO K = I+1, N
         TMPC=TMPC+R(I,K)*X(K)
       ENDDO
       !X(I) = (B(I) - SUM(R(I,I+1:N) * X(I+1:N))) / R(I,I)
       X(I) = (B(I) - TMPC) / R(I,I)
    END DO
END SUBROUTINE QBSSOLVE


!*********************************************************************
!
! solve a system of equations using QQRDC in quadruple precision
!
!*********************************************************************
SUBROUTINE QQRSOLVE(A, M, N, B, X)
    IMPLICIT NONE
    QDPCOMPLEX :: A(:,:), B(:), X(:)
    INTEGER           :: M, N
    ! local
    QDPCOMPLEX, ALLOCATABLE :: QQ(:,:), R(:,:)
    QDPCOMPLEX, ALLOCATABLE :: QB(:)
    INTEGER, ALLOCATABLE :: P(:)
    INTEGER :: I

    ALLOCATE(QQ(M,M))
    ALLOCATE(R(M,N))
    ALLOCATE(P(N))
    ALLOCATE(QB(M))

    ! A = Q R P^T
    CALL QQRDC(A, M, N, QQ, R, P, .TRUE.)

    ! A x = b
    ! Q R P^T x = b
    ! R P^T x = Q^T b
    QB = MATMULX_QDC(TRANSPOSE_QDC(QQ), B)

    ! R P^T x = Q^T b
    ! P^T x = R^-1 Q^T b
    CALL QBSSOLVE(R, N, QB(1:N), X)

    ! x = P P^T x
    X = X(P)

    ! write euclidean matrix norm to stdout
    ! WRITE (*,"(E14.6)",ADVANCE="NO") &
    !   SQRT(SUM(ABS(R)**2))

    ! write norm of residual vector to stdout
    ! WRITE (*,"(E14.6)",ADVANCE="NO") &
    !   SQRT(REAL(SUM(QB(N+1:M)*CONJG(QB(N+1:M)))))

    DEALLOCATE(QQ,R,P,QB)
END SUBROUTINE QQRSOLVE

   !*****************************
   ! matrix multiplication for qdc
   !*****************************
   FUNCTION MATMUL_QDC( A, B )
       QDPCOMPLEX :: A(:,:), B(:,:)
       QDPCOMPLEX :: MATMUL_QDC(SIZE(A,1),SIZE(B,2))
       ! logical 
       INTEGER    :: I, J, K
       QDPCOMPLEX :: R(SIZE(A,1),SIZE(B,2))

#ifndef qd_emulate
       MATMUL_QDC = MATMUL( A, B ) 
       RETURN
#endif
       
       IF ( SIZE( A,2 ) /= SIZE( B, 1 ) ) THEN
          WRITE(*,*)'ERROR, MATMUL_QDC has wrong shape',SIZE(A,2),SIZE(B,1)
       ENDIF
       
       R = 0 
       DO I = 1, SIZE( A, 1 )      
          DO J = 1, SIZE( B, 2 )      
             DO K = 1, SIZE( B, 1 )      
                R(I,J) = R(I,J) + A(I,K)*B(K,J)
             ENDDO
          ENDDO
       ENDDO
       MATMUL_QDC = R
       
   END FUNCTION MATMUL_QDC 

   !*****************************
   ! matrix multiplication for qdc
   !*****************************
   FUNCTION MATMULX_QDC( A, B )
       QDPCOMPLEX :: A(:,:), B(:)
       QDPCOMPLEX :: MATMULX_QDC(SIZE(A,1))
       ! logical 
       INTEGER    :: I, J, K
       QDPCOMPLEX :: R(SIZE(A,1))

#ifndef qd_emulate
       MATMULX_QDC = MATMUL( A, B ) 
       RETURN
#endif
       
       IF ( SIZE( A,2 ) /= SIZE( B, 1 ) ) THEN
          WRITE(*,*)'ERROR, MATMULX_QDC has wrong shape',SIZE(A,2),SIZE(B,1)
       ENDIF
       
       R = 0 
       DO I = 1, SIZE( A, 1 )      
          DO K = 1, SIZE( B, 1 )      
             R(I) = R(I) + A(I,K)*B(K)
          ENDDO
       ENDDO
       MATMULX_QDC = R
       
   END FUNCTION MATMULX_QDC 

   !*****************************
   ! matrix transposition
   !*****************************
   FUNCTION TRANSPOSE_QDC( A )
       QDPCOMPLEX :: A(:,:)
       QDPCOMPLEX :: TRANSPOSE_QDC(SIZE(A,2),SIZE(A,1))
       ! logical 
       INTEGER    :: I, J

#ifndef qd_emulate
       TRANSPOSE_QDC = TRANSPOSE( A )
       RETURN
#endif
       
       DO I = 1, SIZE( A, 1 )      
          DO J = 1, SIZE( A, 2 )      
              TRANSPOSE_QDC(J,I) = A( I, J )
          ENDDO
       ENDDO
       
   END FUNCTION TRANSPOSE_QDC 

END MODULE pade_fit
