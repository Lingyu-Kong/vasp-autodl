!#define dotiming
!#define debug
#include "symbol.inc"

MODULE choleski
  USE prec
  USE dfast
CONTAINS
!************************ SUBROUTINE ORTHCH ****************************
! RCS:  $Id: choleski2.F,v 1.5 2002/04/16 07:28:38 kresse Exp $
!
! this subroutine orthonormalizes a set of complex (wave-)functions
! using a Choleski-decomposition of the overlap matrix (O = L L^H)
! in conjunction with inversion of the result of this decomposition
! (U --> U^-1). If we have arbitrary functions {|cfin_i>} on input,
! we have first to build up the overlap matrix OVL_ij=<cfin_i|cfin_j>
! then we have to decompose it (OVL_ij=L_ik U_kj), have to invert
! U_ij and then to form the output set |cfout_i>=U^-1_ji|cfin_j>. As
! one can show easily it holds: <cfout_i|cfout_j>=delta_ij !! Voila!
!
!> @details @ref openmp :
!> NSTRIP is multiplied with openmp::omp_nthreads.
!> This improves the performance of the ::orth1 and ::orth1_distri
!> calls under OpenMP.
!
!***********************************************************************

      SUBROUTINE ORTHCH(WDES,W, LOVERL,LMDIM,CQIJ,NKSTART,NKSTOP)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE scala
      USE wave_high
      USE wave
      USE wave_mpi
      USE tutor, ONLY: ORTHONORMALIZATION_FAILED, argument, isError
      USE openmp, ONLY : omp_nthreads
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavespin) W
      TYPE (wavedes)  WDES
      TYPE (wavedes1) WDES1
      INTEGER, OPTIONAL :: NKSTART,NKSTOP

      INTEGER :: MY_NKSTART,MY_NKSTOP
      LOGICAL   LOVERL
      OVERLAP   CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      GDEF,ALLOCATABLE,TARGET :: CPROW(:,:),COVL(:,:)
! redistributed plane wave coefficients
      COMPLEX(q), POINTER :: CW_RED(:,:)
      GDEF   , POINTER :: CPROW_RED(:,:),CPROJ_RED(:,:)
      LOGICAL DO_REDIS
      TYPE (REDIS_PW_CTR),POINTER :: H_PW

      LOGICAL IS_ALREADY_REDIS   ! are the plane wave coefficients already redistributed
      LOGICAL LASYNC_            ! ansyncronous exchange of pw coefficients

      PROFILING_START('orthch')

#ifdef _OPENACC
!$ACC WAIT IF(ACC_EXEC_ON)
      PUSH_ACC_EXEC_ON(.TRUE.)
      ACC_ASYNC_Q = ACC_ASYNC_ASYNC

      CALL ACC_COPYIN_TYPED_VAR(WDES)
      CALL ACC_COPYIN_TYPED_VAR(W)

!$ACC ENTER DATA COPYIN(CQIJ) __IF_ASYNC__
#endif

#ifdef MPI
      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE
      NCPU=WDES%COMM_INTER%NCPU ! number of procs involved in band dis.
#else
      NCPU=1
      NODE_ME=0
      IONODE=0
#endif
!-----------------------------------------------------------------------
! determine whether redistribution is required
!-----------------------------------------------------------------------
      IF (NCPU /= 1) THEN

        DO_REDIS=.TRUE.
        NRPLWV_RED=WDES%NRPLWV/NCPU
        NPROD_RED =WDES%NPROD /NCPU
        LASYNC_=LASYNC
        ! it is possible that the bands are distributed of plane wave coefficients upon entry
        IS_ALREADY_REDIS=W%OVER_BAND
        ! of course no need to do it in overlap with calculations
        IF (IS_ALREADY_REDIS) LASYNC_=.FALSE.
      ELSE

        DO_REDIS=.FALSE.
        NRPLWV_RED=WDES%NRPLWV
        NPROD_RED =WDES%NPROD

      ENDIF
      NB_TOT=WDES%NB_TOT
      NBANDS=WDES%NBANDS
      NSTRIP=NSTRIP_STANDARD
!$    NSTRIP=NSTRIP*omp_nthreads

! allocate work space
      ALLOCATE(CPROW(WDES%NPROD,NBANDS))

      IF (.NOT.__IF_ACC_OFF__(LscaAWARE)) THEN
         ALLOCATE(COVL(NB_TOT,NB_TOT))
      ELSE
         CALL INIT_scala(WDES%COMM_KIN, NB_TOT)
         ALLOCATE(COVL(SCALA_NP(), SCALA_NQ()))
      ENDIF
!$ACC ENTER DATA CREATE(WDES1,CPROW,COVL) __IF_ASYNC__

      MY_NKSTART=1
      MY_NKSTOP=WDES%NKPTS
      IF(PRESENT(NKSTART))MY_NKSTART=NKSTART
      IF(PRESENT(NKSTOP))MY_NKSTOP=NKSTOP
!=======================================================================
      spin:    DO ISP=1,WDES%ISPIN
      kpoints: DO NK=MY_NKSTART,MY_NKSTOP
#ifdef MPI
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
!=======================================================================
      IF (__IF_ACC_OFF__(LscaAWARE)) CALL INIT_scala(WDES%COMM_KIN,WDES%NB_TOTK(NK,ISP))
      CALL SETWDES(WDES,WDES1,NK)
!   get pointers for redistributed wavefunctions
!   I can not guarantee that this works with all f90 compilers
!   please see comments in wave_mpi.F
      IF (NCPU /= 1) THEN
        CALL SET_WPOINTER(CW_RED,   NRPLWV_RED, NB_TOT, W%CW(1,1,NK,ISP))
        CALL SET_GPOINTER(CPROJ_RED, NPROD_RED, NB_TOT, W%CPROJ(1,1,NK,ISP))
        CALL SET_GPOINTER(CPROW_RED, NPROD_RED, NB_TOT, CPROW(1,1))
      ELSE
        CW_RED    => W%CW(:,:,NK,ISP)
        CPROJ_RED => W%CPROJ(:,:,NK,ISP)
        CPROW_RED => CPROW(:,:)
      ENDIF

!   set number of wavefunctions after redistribution
      NPL = WDES1%NPL
      NPRO= WDES1%NPRO

      CALLMPI( SET_NPL_NPRO(WDES1, NPL, NPRO))

      NPRO_O=NPRO
      IF (.NOT. LOVERL) NPRO_O=0
!=======================================================================
!  calculate overlap matrix (only upper triangle is needed):
!=======================================================================
      IF (DO_REDIS .AND. LASYNC_) THEN
         CALL REDIS_PW_ALLOC(WDES, NSTRIP, H_PW)
         DO NPOS=1,NSTRIP
           CALL REDIS_PW_START(WDES, W%CW(1,NPOS,NK,ISP), NPOS, H_PW)
        ENDDO
      ENDIF

      CALL OVERL(WDES1,LOVERL,LMDIM,CQIJ(1,1,1,ISP),W%CPROJ(1,1,NK,ISP),CPROW(1,1))
    ! redistribute everything

      statmg
      IF (DO_REDIS) THEN
        CALL REDIS_PROJ(WDES1, NBANDS, CPROW(1,1))
        stotmgt("RD_PROJ")
        CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
        stotmgt("RD_PROJ")
        IF (.NOT. LASYNC_ .AND. .NOT. IS_ALREADY_REDIS)  CALL REDIS_PW  (WDES1, NBANDS, W%CW   (1,1,NK,ISP))
        stotmgt("RD_PW  ")
      ENDIF

      IF (.NOT.__IF_ACC_OFF__(LscaAWARE)) THEN
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(COVL) __IF_ASYNC__
         DO N=1,NB_TOT
            DO I=1,NB_TOT
               COVL(I,N)=zero
            ENDDO
         ENDDO
       ENDIF
!
! there is a strange bug in the PII optimized blas DGEMM, which seems
! to access in certain instances data beyond the last array element
! if a matrix is multiplied with a vector (second call to ORTH1)
! to work around this I avoid calling ORTH1 with NB_TOT-NPOS+1=1
      DO NPOS=1,NBANDS-NSTRIP,NSTRIP
        IF (DO_REDIS .AND. LASYNC_) THEN
        DO NPOS_=NPOS,NPOS+NSTRIP-1
          CALL REDIS_PW_STOP (WDES, W%CW(1,NPOS_,NK,ISP), NPOS_, H_PW)
          IF (NPOS_+NSTRIP<=NBANDS) &
          CALL REDIS_PW_START(WDES, W%CW(1,NPOS_+NSTRIP,NK,ISP), NPOS_+NSTRIP, H_PW)
        ENDDO
        ENDIF

        NPOS_RED  =(NPOS-1)*NCPU+1
        NSTRIP_RED=NSTRIP*NCPU

        IF (.NOT.__IF_ACC_OFF__(LscaAWARE)) THEN
           CALL ORTH1('U',CW_RED(1,1),CW_RED(1,NPOS_RED),CPROJ_RED(1,1), &
                CPROW_RED(1,NPOS_RED),NB_TOT, &
                NPOS_RED,NSTRIP_RED,NPL,NPRO_O,NRPLWV_RED,NPROD_RED,COVL(1,1))
        ELSE
           CALL ORTH1_DISTRI('U',CW_RED(1,1),CW_RED(1,NPOS_RED),CPROJ_RED(1,1), &
                CPROW_RED(1,NPOS_RED),NB_TOT, &
                NPOS_RED,NSTRIP_RED,NPL,NPRO_O,NRPLWV_RED,NPROD_RED,COVL(1,1), & 
                WDES%COMM_KIN, WDES%NB_TOTK(NK,ISP))
        ENDIF
      ENDDO

      IF (DO_REDIS .AND. LASYNC_) THEN
      DO NPOS_=NPOS,NBANDS
          CALL REDIS_PW_STOP (WDES, W%CW(1,NPOS_,NK,ISP), NPOS_, H_PW)
      ENDDO
      ENDIF
      
      NPOS_RED  =(NPOS-1)*NCPU+1
      NSTRIP_RED=(NBANDS-NPOS+1)*NCPU

      IF (.NOT.__IF_ACC_OFF__(LscaAWARE)) THEN
         CALL ORTH1('U',CW_RED(1,1),CW_RED(1,NPOS_RED),CPROJ_RED(1,1), &
              CPROW_RED(1,NPOS_RED),NB_TOT, &
              NPOS_RED,NSTRIP_RED,NPL,NPRO_O,NRPLWV_RED,NPROD_RED,COVL(1,1))
      ELSE
         CALL ORTH1_DISTRI('U',CW_RED(1,1),CW_RED(1,NPOS_RED),CPROJ_RED(1,1), &
              CPROW_RED(1,NPOS_RED),NB_TOT, &
              NPOS_RED,NSTRIP_RED,NPL,NPRO_O,NRPLWV_RED,NPROD_RED,COVL(1,1), &  
              WDES%COMM_KIN, WDES%NB_TOTK(NK,ISP))
      ENDIF

      IF (DO_REDIS .AND. LASYNC_) CALL REDIS_PW_DEALLOC(H_PW)

      stotmgt("IN_PROD")
      IF (.NOT.__IF_ACC_OFF__(LscaAWARE)) THEN
         CALLMPI( M_sum_g(WDES%COMM_KIN,COVL(1,1),NB_TOT*NB_TOT))
      ENDIF
      stotmgt("GLB_SUM")

#ifdef debug
      io_begin
      CALL DUMP_HAM( "overlap",WDES, COVL)
      io_end
#endif
!=======================================================================
! Choleski-decomposition of the overlap matrix + inversion of the result
! calling LAPACK-routines ZPOTRF (decomposition) and ZTRTRI (inversion):
!=======================================================================
      PROFILING_START('choleski')

      IF (.NOT.__IF_ACC_OFF__(LscaAWARE)) THEN
         IF (__IF_ACC_OFF__(LscaLAPACK.AND.LscaLU)) THEN
            CALL pPOTRF_TRTRI(WDES%COMM_KIN, COVL(1,1),NB_TOT, WDES%NB_TOTK(NK,ISP))
            CALLMPI( M_sum_g(WDES%COMM_KIN,COVL(1,1), NB_TOT*NB_TOT))
         ELSE
            INFO=0
#ifdef gammareal
            CALL __DPOTRF__ &
#else
            CALL __ZPOTRF__ &
#endif
                 & ('U',WDES%NB_TOTK(NK,ISP),COVL(1,1),NB_TOT,INFO)
            IF (INFO/=0) THEN
               CALL vtutor%write(isError, ORTHONORMALIZATION_FAILED, argument(ival=[1, NK, ISP]))
            ENDIF
#ifdef gammareal
            CALL __DTRTRI__ &
#else
            CALL __ZTRTRI__ &
#endif
                 & ('U','N',WDES%NB_TOTK(NK,ISP),COVL(1,1),NB_TOT,INFO)
            IF (INFO/=0) THEN
               CALL vtutor%write(isError, ORTHONORMALIZATION_FAILED, argument(ival=[2, NK, ISP]))
            ENDIF
         ENDIF
      ELSE
         CALL BG_pPOTRF_TRTRI(COVL(1,1), WDES%NB_TOTK(NK,ISP), INFO)
         IF (INFO/=0) THEN
            CALL vtutor%write(isError, ORTHONORMALIZATION_FAILED, argument(ival=[3, NK, ISP]))
         ENDIF
      ENDIF

      PROFILING_STOP('choleski')

      stotmgt("ZPOTRF")
!=======================================================================
!  construct the orthogonal set:
!=======================================================================
      IF (.NOT.__IF_ACC_OFF__(LscaAWARE)) THEN
         CALL LINCOM('U',CW_RED,CPROJ_RED,COVL(1,1), &
              WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP),NPL,NPRO,NRPLWV_RED,NPROD_RED,NB_TOT, &
              CW_RED,CPROJ_RED)
      ELSE
         CALL LINCOM_DISTRI('U',CW_RED(1,1),CPROJ_RED(1,1),COVL(1,1), &
              WDES%NB_TOTK(NK,ISP),NPL,NPRO,NRPLWV_RED,NPROD_RED,NB_TOT, & 
              WDES%COMM_KIN, NBLK)
      ENDIF

      stotmgt("LINCOM")

     !  back redistribution
      IF (DO_REDIS) THEN
        CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
        IF (LASYNC_ .OR. IS_ALREADY_REDIS) THEN
           W%OVER_BAND=.TRUE.
        ELSE
           ! if the routine was entered with already redistributed wave-functions
           ! then no need to redistribute them back
           CALL REDIS_PW  (WDES1, NBANDS, W%CW   (1,1,NK,ISP))
        ENDIF
      ENDIF
      stotmgt("RD_PW  ")
!=======================================================================
      ENDDO kpoints
      ENDDO spin
!=======================================================================
!$ACC EXIT DATA DELETE(CPROW,COVL) __IF_ASYNC__
      DEALLOCATE(CPROW,COVL)

#ifdef _OPENACC
      CALL ACC_DELETE_TYPED_VAR(WDES1)

      CALL ACC_DELETE_TYPED_VAR(WDES)
      CALL ACC_COPYOUT_TYPED_VAR(W)
!$ACC EXIT DATA DELETE(CQIJ) __IF_ASYNC__

!$ACC WAIT IF(ACC_EXEC_ON)
      ACC_ASYNC_Q = ACC_ASYNC_SYNC

      POP_ACC_EXEC_ON
#endif

      PROFILING_STOP('orthch')

      RETURN
    END SUBROUTINE ORTHCH


!************************ SUBROUTINE ORTHCH_DUAL ***********************
! RCS:  $Id: choleski2.F,v 1.5 2002/04/16 07:28:38 kresse Exp $
!
! this subroutine determines a dual set of orbitals with the
! property 
!  < phi_i | phi_j> = delta_ij
! the dual set W_DUAL is passed down as a pointer to a wavespin variable
! if NC potential are used W_DUAL is simply linked to W
! otherwise if W_DUAL is not associated, the required arrays are
! allocated and determined (make sure to NULLIFY W_DUAL before 
!   calling the routined)
! if it is already associated, the routine recalculates W_DUAL from W but
! assumes that the proper allocation was done before
!
!***********************************************************************

    SUBROUTINE ORTHCH_DUAL(WDES,W,W_DUAL, LOVERL,LMDIM,CQIJ)
      USE prec
      USE scala
      USE wave_high
      USE wave
      USE wave_mpi
      IMPLICIT NONE

      TYPE (wavespin), TARGET  ::  W
      TYPE (wavespin), POINTER ::  W_DUAL
      TYPE (wavedes)  WDES
      TYPE (wavedes1) WDES1

      LOGICAL   LOVERL
      INTEGER   LMDIM
      OVERLAP   CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      GDEF,ALLOCATABLE,TARGET :: COVL(:,:)
! redistributed plane wave coefficients
      TYPE (wavefuna)    WA             ! array pointers
      TYPE (wavefuna)    WA_DUAL        ! array pointers
      INTEGER   NODE_ME, IONODE, NCPU
      INTEGER   NB_TOT, NBANDS, NSTRIP, NSTRIP_ACT, ISP, NK, N, I, INFO, NPOS, NPOS_RED, NSTRIP_RED
! scaAWARE is often a little bit tricky to implement
! routine was tested and currenty LscaAWARE_LOCAL is set to LscaAWARE
! see below
      LOGICAL :: LscaAWARE_LOCAL=.FALSE.
#ifdef MPI
      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE
      NCPU=WDES%COMM_INTER%NCPU ! number of procs involved in band dis.
#else
      NCPU=1
      NODE_ME=0
      IONODE=0
#endif
      IF (.NOT. LOVERL) THEN
         W_DUAL => W
         RETURN
      ENDIF

! if there are any subspected problems with duality comment out
! the following line (or set LscaAWARE_LOCAL =.FALSE.)
      LscaAWARE_LOCAL=LscaAWARE

      IF (.NOT. ASSOCIATED(W_DUAL)) THEN
         ALLOCATE(W_DUAL)
         CALL ALLOCW( WDES, W_DUAL)
!#define test_w_dual
#ifdef test_w_dual
         W_DUAL%CPROJ=W%CPROJ
         W_DUAL%CW   =W%CW
         W_DUAL%CELTOT=W%CELTOT
         W_DUAL%FERTOT=W%FERTOT
#endif
      ENDIF

#ifndef test_w_dual
      W_DUAL%CPROJ=W%CPROJ
      W_DUAL%CW   =W%CW
      W_DUAL%CELTOT=W%CELTOT
      W_DUAL%FERTOT=W%FERTOT
#endif
!-----------------------------------------------------------------------
! determine whether redistribution is required
!-----------------------------------------------------------------------
      NB_TOT=WDES%NB_TOT
      NBANDS=WDES%NBANDS
      NSTRIP=NSTRIP_STANDARD

      IF (.NOT. LscaAWARE_LOCAL) THEN
         ALLOCATE(COVL(NB_TOT,NB_TOT))
      ELSE
         CALL INIT_scala(WDES%COMM_KIN, NB_TOT)
         ALLOCATE(COVL(SCALA_NP(), SCALA_NQ()))
      ENDIF
!=======================================================================
      spin:    DO ISP=1,WDES%ISPIN
      kpoints: DO NK=1,WDES%NKPTS
#ifdef MPI
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
!=======================================================================
      IF (LscaAWARE_LOCAL) CALL INIT_scala(WDES%COMM_KIN,WDES%NB_TOTK(NK,ISP))

      CALL SETWDES(WDES,WDES1,NK)
      WA     =ELEMENTS(W, WDES1, ISP)
      WA_DUAL=ELEMENTS(W_DUAL, WDES1, ISP)
!=======================================================================
!  calculate overlap matrix (only upper triangle is needed):
!=======================================================================

      statmg
      IF (WDES%DO_REDIS) THEN
         CALL REDISTRIBUTE_PROJ(WA)
         CALL REDISTRIBUTE_PW(WA)
         CALL REDISTRIBUTE_PROJ(WA_DUAL)
         CALL REDISTRIBUTE_PW(WA_DUAL)
      ENDIF

      IF (.NOT. LscaAWARE_LOCAL) THEN
         DO N=1,NB_TOT
            DO I=1,NB_TOT
               COVL(I,N)=(0._q,0._q)
            ENDDO
         ENDDO
      ENDIF
!
      DO NPOS=1,NBANDS,NSTRIP
        NSTRIP_ACT=MIN(WDES%NBANDS+1-NPOS,NSTRIP)
        NPOS_RED  =(NPOS-1)*NCPU+1
        NSTRIP_RED=NSTRIP_ACT*NCPU

        IF (.NOT. LscaAWARE_LOCAL) THEN
           ! W^+ x W_DUAL (W_DUAL=W usually)
           CALL ORTH1('U',WA%CW_RED(1,1),WA_DUAL%CW_RED(1,NPOS_RED),WA%CPROJ_RED(1,1), &
                WA_DUAL%CPROJ_RED(1,1),NB_TOT, &
                NPOS_RED,NSTRIP_RED, WDES1%NPL_RED,0,WDES1%NRPLWV_RED,WDES1%NPROD_RED,COVL(1,1))
        ELSE
           CALL ORTH1_DISTRI('U',WA%CW_RED(1,1),WA_DUAL%CW_RED(1,NPOS_RED),WA%CPROJ_RED(1,1), &
                WA_DUAL%CPROJ_RED(1,1),NB_TOT, &
                NPOS_RED,NSTRIP_RED, WDES1%NPL_RED,0,WDES1%NRPLWV_RED,WDES1%NPROD_RED,COVL(1,1), & 
                WDES%COMM_KIN, WDES%NB_TOTK(NK,ISP))
        ENDIF
      ENDDO
      
      stotmgt("IN_PROD")
      IF (.NOT. LscaAWARE_LOCAL) THEN
         CALLMPI( M_sum_g(WDES%COMM_KIN,COVL(1,1),NB_TOT*NB_TOT))
      ENDIF
      stotmgt("GLB_SUM")

#if defined(debug) || defined test_w_dual
      io_begin
      CALL DUMP_HAM( "overlap",WDES, COVL)
      io_end
#endif
!=======================================================================
! Choleski-decomposition of the overlap matrix + inversion of the result
! calling LAPACK-routines ZPOTRF (decomposition) and ZPOTRI (inversion):
!=======================================================================
      IF (.NOT. LscaAWARE_LOCAL) THEN
         IF (LscaLAPACK .AND. LscaLU ) THEN
            CALL pPOTRF_POTRI(WDES%COMM_KIN, COVL(1,1),NB_TOT, WDES%NB_TOTK(NK,ISP))
            CALLMPI( M_sum_g(WDES%COMM_KIN,COVL(1,1), NB_TOT*NB_TOT))
         ELSE
            INFO=0
#ifdef gammareal
            CALL DPOTRF &
#else
            CALL ZPOTRF &
#endif
                 & ('U',WDES%NB_TOTK(NK,ISP),COVL(1,1),NB_TOT,INFO)
            IF (INFO/=0) THEN
               CALL vtutor%error("LAPACK: Routine ZPOTRF failed! " // str(INFO) // " " // str(NK) // &
                  " " // str(ISP))
            ENDIF
#ifdef gammareal
            CALL DPOTRI &
#else
            CALL ZPOTRI &
#endif
                 & ('U',WDES%NB_TOTK(NK,ISP),COVL(1,1),NB_TOT,INFO)
            IF (INFO/=0) THEN
               CALL vtutor%error("LAPACK: Routine ZTRTRI failed! " // str(INFO) // " " // str(NK) // &
                  " " // str(ISP))
            ENDIF
         ENDIF
         ! add lower triangle not calculated by POTRI
         DO N=1,NB_TOT
            DO I=N+1,NB_TOT
               COVL(I,N)=GCONJG(COVL(N,I))
            ENDDO
         ENDDO
      ELSE
         CALL BG_pPOTRF_POTRI(COVL(1,1), WDES%NB_TOTK(NK,ISP), INFO)
         IF (INFO/=0) THEN
            CALL vtutor%error("LAPACK: Routine ZTRTRI failed! " // str(INFO) // " " // str(NK) // " &
               &" // str(ISP))
         ENDIF
      ENDIF
#ifdef debug
      io_begin
      CALL DUMP_HAM( "transformation",WDES, COVL)
      io_end
#endif

      stotmgt("ZPOTRF")
!=======================================================================
!  construct the dual set by inplace transformation of W_DUAL
!=======================================================================
      IF (.NOT. LscaAWARE_LOCAL) THEN
         CALL LINCOM('F',WA_DUAL%CW_RED,WA_DUAL%CPROJ_RED,COVL(1,1), &
              WDES%NB_TOTK(NK,ISP),WDES%NB_TOTK(NK,ISP), & 
              WDES1%NPL_RED, WDES1%NPRO_RED, WDES1%NRPLWV_RED, WDES1%NPROD_RED, NB_TOT, &
              WA%CW_RED, WA%CPROJ_RED)
      ELSE
         ! POTRF set only the upper triangle of the transformation matrix
         ! add the lower triangle using DISTRI_HERM
         CALL LINCOM_DISTRI_HERM('F',WA_DUAL%CW_RED(1,1),WA_DUAL%CPROJ_RED(1,1),COVL(1,1), &
              WDES%NB_TOTK(NK,ISP), & 
              WDES1%NPL_RED, WDES1%NPRO_RED, WDES1%NRPLWV_RED, WDES1%NPROD_RED, NB_TOT, & 
              WDES%COMM_KIN, NBLK)
      ENDIF

      stotmgt("LINCOM")

     !  back redistribution
      IF (WDES%DO_REDIS) THEN
        CALL REDISTRIBUTE_PROJ(WA)
        CALL REDISTRIBUTE_PW  (WA)
        CALL REDISTRIBUTE_PROJ(WA_DUAL)
        CALL REDISTRIBUTE_PW  (WA_DUAL)
      ENDIF
      stotmgt("RD_PW  ")
!=======================================================================
      ENDDO kpoints
      ENDDO spin
!=======================================================================
      DEALLOCATE(COVL)

      RETURN
      END SUBROUTINE ORTHCH_DUAL


!************************ SUBROUTINE DUMP_S ****************************
!
! dump the overlap matrix between the occupied states
! this is essentially identical to the previous routine
!
!***********************************************************************


      SUBROUTINE DUMP_S(WDES,W, LOVERL,LMDIM,CQIJ)
      USE prec
      USE scala
      USE wave_high
      USE wave
      USE wave_mpi
      USE dfast
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavespin) W
      TYPE (wavedes)  WDES
      TYPE (wavedes1) WDES1

      LOGICAL   LOVERL
      OVERLAP   CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
      GDEF,ALLOCATABLE,TARGET :: CPROW(:,:),COVL(:,:)

      ! redistributed plane wave coefficients
      COMPLEX(q), POINTER :: CW_RED(:,:)
      GDEF   , POINTER :: CPROW_RED(:,:),CPROJ_RED(:,:)
      LOGICAL DO_REDIS
      REAL(q) :: NELECT
#ifdef MPI
      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE
      NCPU=WDES%COMM_INTER%NCPU ! number of procs involved in band dis.
#else
      NCPU=1
      NODE_ME=0
      IONODE=0
#endif
      IF (NCPU /= 1) THEN

         DO_REDIS=.TRUE.
         NRPLWV_RED=WDES%NRPLWV/NCPU
         NPROD_RED =WDES%NPROD /NCPU
      ELSE
         DO_REDIS=.FALSE.
         NRPLWV_RED=WDES%NRPLWV
         NPROD_RED =WDES%NPROD
      ENDIF

      NB_TOT=WDES%NB_TOT
      NBANDS=WDES%NBANDS

      NSTRIP=NSTRIP_STANDARD

! allocate work space
      ALLOCATE(CPROW(WDES%NPROD,NBANDS),COVL(NB_TOT,NB_TOT))

      NELECT=0

      spin:    DO ISP=1,WDES%ISPIN
      kpoints: DO NK=1,WDES%NKPTS
#ifdef MPI
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
      CALL SETWDES(WDES,WDES1,NK)
! get pointers for redistributed wavefunctions
      IF (NCPU /= 1) THEN
         CALL SET_WPOINTER(CW_RED,   NRPLWV_RED, NB_TOT, W%CW(1,1,NK,ISP))
         CALL SET_GPOINTER(CPROJ_RED, NPROD_RED, NB_TOT, W%CPROJ(1,1,NK,ISP))
         CALL SET_GPOINTER(CPROW_RED, NPROD_RED, NB_TOT, CPROW(1,1))
      ELSE
         CW_RED    => W%CW(:,:,NK,ISP)
         CPROJ_RED => W%CPROJ(:,:,NK,ISP)
         CPROW_RED => CPROW(:,:)
      ENDIF

! set number of wavefunctions after redistribution
      NPL = WDES1%NPL
      NPRO= WDES1%NPRO

      CALLMPI( SET_NPL_NPRO(WDES1, NPL, NPRO))

      NPRO_O=NPRO
      IF (.NOT. LOVERL) NPRO_O=0

      CALL OVERL(WDES1, LOVERL,LMDIM,CQIJ(1,1,1,ISP), W%CPROJ(1,1,NK,ISP),CPROW(1,1))

! redistribute everything
      IF (DO_REDIS) THEN
         CALL REDIS_PROJ(WDES1, NBANDS, CPROW(1,1))
         CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
         CALL REDIS_PW  (WDES1, NBANDS, W%CW   (1,1,NK,ISP))
      ENDIF

! calculate overlap
      DO N=1,NB_TOT
         DO I=1,NB_TOT
            COVL(I,N)=(0._q,0._q)
         ENDDO
      ENDDO

      DO NPOS=1,NBANDS-NSTRIP,NSTRIP
         NPOS_RED  =(NPOS-1)*NCPU+1
         NSTRIP_RED=NSTRIP*NCPU

         CALL ORTH1('U',CW_RED(1,1),CW_RED(1,NPOS_RED),CPROJ_RED(1,1), &
              CPROW_RED(1,NPOS_RED),NB_TOT, &
              NPOS_RED,NSTRIP_RED,NPL,NPRO_O,NRPLWV_RED,NPROD_RED,COVL(1,1))
      ENDDO

      
      NPOS_RED  =(NPOS-1)*NCPU+1
      NSTRIP_RED=(NBANDS-NPOS+1)*NCPU

      CALL ORTH1('U',CW_RED(1,1),CW_RED(1,NPOS_RED),CPROJ_RED(1,1), &
           CPROW_RED(1,NPOS_RED),NB_TOT, &
           NPOS_RED,NSTRIP_RED,NPL,NPRO_O,NRPLWV_RED,NPROD_RED,COVL(1,1))

      DO N=1,NB_TOT
         NELECT=NELECT+ COVL(N, N)*W%FERWE(N,NK,ISP)*WDES%RSPIN*WDES%WTKPT(NK)
      ENDDO

      CALLMPI( M_sum_g(WDES%COMM_KIN,COVL(1,1),NB_TOT*NB_TOT))

! back redistribution
      IF (DO_REDIS) THEN
        CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
        CALL REDIS_PW  (WDES1, NBANDS, W%CW   (1,1,NK,ISP))
      ENDIF

      WRITE(*,*) 'NK=',NK,'ISP=',ISP
      NPL2=MIN(10,NB_TOT)
      DO  N1=1,NPL2
         WRITE(*,1)N1,(REAL( COVL(N1,N2) ,KIND=q) ,N2=1,NPL2)
      ENDDO
      WRITE(*,*)
#ifndef gammareal
      DO N1=1,NPL2
         WRITE(*,2)N1,(AIMAG(COVL(N1,N2)),N2=1,NPL2)
      ENDDO
#endif
      WRITE(*,*)
      WRITE(*,'(20E10.3)') 0,(REAL( COVL(N1,N1)-1 ,KIND=q) ,N1=1,NPL2)
      WRITE(*,*)
1     FORMAT(1I2,3X,20F9.5)
2     FORMAT(1I2,3X,20E9.1)

      ENDDO kpoints
      ENDDO spin

      WRITE(*,'(A,F14.8)') 'total number of electrons NELECT=',NELECT
      DEALLOCATE(CPROW,COVL)

      RETURN
      END SUBROUTINE

!                              __________________
!_____________________________/ CALC_PAW_OVERLAP \______________________________
!
!> @brief Calculates the S-overlap between two wavefunctions WA1 and WA2.
!>
!> The wavefunction descriptors associated with WA1 and WA2 have to be
!> compatible, preferably identical. WA1 and WA2 may refer to the same storage
!> location. This is detected and MPI redistribution is handled correctly for
!> this case.
!>
!> @param[in] WA1             Bra state
!> @param[in] WA2             Ket state
!> @param[in] CQIJ            Augmentation charges
!> @param[in,out] OVERLAP_MAT Calculated overlap matrix
!> @param[in] ADD             If true, result is added to OVERLAP_MAT instead
!>                            of overwriting it.
!> @param[in] ONLY_AUG        If true, only augmentation part is calculated.
!_______________________________________________________________________________
      SUBROUTINE CALC_PAW_OVERLAP(WA1, WA2, CQIJ, OVERLAP_MAT, ADD, ONLY_AUG)
         USE wave_struct_def, ONLY: wavefuna, &
                                    wavedes1
         USE wave_high,       ONLY: OVERL, &
                                    REDISTRIBUTE_PW, &
                                    REDISTRIBUTE_PROJ, &
                                    NEWWAVA_PROJ, &
                                    DELWAVA_PROJ
         USE wave_mpi,        ONLY: SET_WPOINTER
         USE tutor,           ONLY: vtutor
         USE string,          ONLY: str
         IMPLICIT NONE
   
         TYPE(wavefuna),      INTENT(IN)     :: WA1, WA2
         OVERLAP,             INTENT(IN)     :: CQIJ(:, :, :, :)
         GDEF,                INTENT(INOUT)  :: OVERLAP_MAT(:, :)
         LOGICAL, OPTIONAL,   INTENT(IN)     :: ADD
         LOGICAL, OPTIONAL,   INTENT(IN)     :: ONLY_AUG
   
         INTEGER        :: NB_TOT
         INTEGER        :: STRIP
         INTEGER        :: NPL
         LOGICAL        :: WA_SAME
         TYPE(wavefuna) :: WOVL
         TYPE(wavedes1) :: WDES1
         LOGICAL        :: MY_ADD, MY_ONLY_AUG
   
         MY_ONLY_AUG = .FALSE.; IF (PRESENT(ONLY_AUG)) MY_ONLY_AUG = ONLY_AUG
         MY_ADD = .FALSE.; IF (PRESENT(ADD)) MY_ADD = ADD
   
         WA_SAME = ASSOCIATED(WA1%CW, WA2%CW)
   
         WDES1 = WA1%WDES1
         NB_TOT = WDES1%NB_TOT
   
         IF (ANY(SHAPE(OVERLAP_MAT) /= NB_TOT)) &
            CALL vtutor%bug("CALC_PAW_OVERLAP: Shape mismatch: " // &
               str(NB_TOT) // ", " // str(SIZE(OVERLAP_MAT, 1)) // ", " // str(SIZE(OVERLAP_MAT, 2)), &
               __FILE__, __LINE__)
   
         ! Allocate only projections
         CALL NEWWAVA_PROJ(WOVL, WDES1)
   
         IF (MY_ADD) THEN
#if defined(MPI) || defined(MPI_CHAIN)
            OVERLAP_MAT = OVERLAP_MAT / WDES1%COMM%NCPU
#endif
         ELSE
            OVERLAP_MAT = 0
         ENDIF
   
         ! Reuse plane-wave components of WA2
         WOVL%CW=>WA2%CW
         IF (WDES1%DO_REDIS) THEN
            CALL SET_WPOINTER(WOVL%CW_RED, WDES1%NRPLWV_RED, NB_TOT, WOVL%CW(1, 1))
         ELSE
            WOVL%CW_RED=>WOVL%CW
         ENDIF
         
         ! Compute \sum_j Q_ij < p_j | W_P >
         CALL OVERL(WDES1, WDES1%LOVERL, SIZE(CQIJ, 1), CQIJ(1, 1, 1, 1), WA2%CPROJ(1, 1), WOVL%CPROJ(1, 1))
   
         ! Redistribute plane-wave components
         IF (.NOT. MY_ONLY_AUG) THEN
            CALL REDISTRIBUTE_PW(WA1)
            IF (.NOT. WA_SAME) CALL REDISTRIBUTE_PW(WOVL)
         ENDIF
         ! Redistribute projections
         CALL REDISTRIBUTE_PROJ(WA1)
         CALL REDISTRIBUTE_PROJ(WOVL)
   
         NPL = WDES1%NPL_RED
         ! NPL == 0 skips calculation of plane-wave part
         IF (MY_ONLY_AUG) NPL = 0
   
         ! Compute < WA1_m | S | WA2_n >
         DO STRIP = 1, NB_TOT - NSTRIP_STANDARD_GLOBAL, NSTRIP_STANDARD_GLOBAL
            CALL ORTH2( &
               WA1%CW_RED(1, 1), WOVL%CW_RED(1, STRIP), WA1%CPROJ_RED(1, 1), WOVL%CPROJ_RED(1, STRIP), &
               NB_TOT, STRIP, NSTRIP_STANDARD_GLOBAL, NPL, &
               WDES1%NPRO_O_RED, WDES1%NRPLWV_RED, WDES1%NPROD_RED, OVERLAP_MAT(1, 1))
         ENDDO
   
         CALL ORTH2( &
            WA1%CW_RED(1, 1), WOVL%CW_RED(1, STRIP), WA1%CPROJ_RED(1, 1), WOVL%CPROJ_RED(1, STRIP), &
            NB_TOT, STRIP, NB_TOT-STRIP+1, NPL, &
            WDES1%NPRO_O_RED, WDES1%NRPLWV_RED, WDES1%NPROD_RED, OVERLAP_MAT(1, 1))
   
         CALLMPI( M_sum_g(WDES1%COMM_KIN, OVERLAP_MAT(1, 1), NB_TOT * NB_TOT))
   
         ! Return to original plane-wave distribution
         IF (.NOT. MY_ONLY_AUG) THEN
            CALL REDISTRIBUTE_PW(WA1)
            IF (.NOT. WA_SAME) CALL REDISTRIBUTE_PW(WA2)
         ENDIF
         ! Return to original projection distribution
         CALL REDISTRIBUTE_PROJ(WA1)
   
         NULLIFY(WOVL%CW)
         CALL DELWAVA_PROJ(WOVL)
      END SUBROUTINE

END MODULE
