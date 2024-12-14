#include "symbol.inc"
!***********************************************************************
!
!***********************************************************************

      MODULE fock_frc

      USE prec
      USE fock
#ifdef _OPENACC
      USE mopenacc
#endif
      IMPLICIT NONE

      PUBLIC :: FOCK_FORCE

      PRIVATE

#ifdef _OPENACC
      INTEGER :: ACC_ASYNC_FFT_AND_GATHER_Q=1
#endif
      CONTAINS

!************************ SUBROUTINE FOCK_FORCE ************************
!
! Calculate the contribution to the Hellmann-Feynman forces and
! the stress tensors (if LSIF=.TRUE.) due to Fock exchange
!
!> @details @ref openmp :
!> the loop over the occupied orbitals owned by the particular
!> MPI-rank (label: mband) is distributed over all available OpenMP
!> threads.
!
!***********************************************************************

      SUBROUTINE FOCK_FORCE( W,LATT_CUR,NONLR_S,NONL_S,P,LMDIM,FORHF,SIFHF,LSIF,IU0,IU6 )
!$ACC ROUTINE(NI_GLOBAL) SEQ
#ifdef _OPENACC
      USE cufft
      USE openacc
#endif
      USE wave_high
      USE nonl_high
      USE lattice
      USE pseudo
      USE full_kpoints
      USE sym_prec

      TYPE (wavespin)     :: W
      TYPE (latt)         :: LATT_CUR
      TYPE (nonlr_struct) :: NONLR_S
      TYPE (nonl_struct)  :: NONL_S
      TYPE (potcar)       :: P(NONLR_S%NTYP)

      INTEGER :: LMDIM

      REAL(q) :: FORHF(3,NONLR_S%NIONS)
      REAL(q) :: SIFHF(3,3)

      LOGICAL :: LSIF

      INTEGER :: IU0,IU6

      ! local variables
      TYPE (wavespin) :: WHF
      TYPE (wavefun1) :: WQ
      TYPE (wavefun1), ALLOCATABLE :: W1(:)
      TYPE (wavefun1), ALLOCATABLE :: WIN(:)
      TYPE (wavedes1), TARGET :: WDESK,WDESQ,WDESQ_IRZ

      TYPE (nonlr_struct), ALLOCATABLE :: FAST_AUG(:)

      GDEF, ALLOCATABLE :: CPROJKXYZ(:,:,:),CPROJXYZ(:,:,:)
      GDEF, ALLOCATABLE :: CRHOLM(:,:)
#if defined(_OPENACC) && defined(gammareal)
      GDEF, POINTER, CONTIGUOUS :: GWORK(:,:)
      GDEF, ALLOCATABLE, TARGET :: GWORK_FWD(:,:)
      GDEF, ALLOCATABLE, TARGET :: GWORK_BCK(:,:)
#else
      GDEF, ALLOCATABLE :: GWORK(:,:)
#endif
#ifndef _OPENACC
      GDEF, ALLOCATABLE :: CDIJ(:,:,:,:,:)
#endif
      GDEF, ALLOCATABLE, TARGET:: CDLM0(:,:),CDLM(:,:)

      REAL(q), ALLOCATABLE :: POTFAK(:,:)
      REAL(q), ALLOCATABLE :: FSG(:)

      REAL(q) :: WEIGHT_Q,WEIGHT
      REAL(q) :: RTMP

      REAL(q), ALLOCATABLE :: ENL(:,:)
      REAL(q) :: SIF(0:6)

      INTEGER :: ISP,NK,NB,NQ,MQ,ISP_IRZ,N,ISPINOR,NI,NT,LMMAXC,NPRO,NIP
      INTEGER :: NB_TOT,NB_TOTK,NBLK,NBLOCK,NBLOCK_ACT,NPOS,NGLB
      INTEGER :: NDIR,IDIR,I,J
      LOGICAL :: LSKIP,LSHIFT
      REAL(q) :: FD

      TYPE (rotation_handle), POINTER :: ROT_HANDLE

#ifdef _OPENACC
      LOGICAL :: DO_BATCHED_FFT=.FALSE.
      INTEGER :: ierr,NS
      REAL(q), ALLOCATABLE :: GFAC_MU_WEIGHTS(:)
#ifndef gammareal
      INTEGER :: fftdims(3)
      INTEGER :: batchedPlan
#else
      INTEGER :: batchedPlanD2Z(2),batchedPlanZ2D(2)
#endif
#endif
      PROFILING_START('fock_force')

      ! early exit if possible
      IF ((.NOT.W%WDES%LOVERL.AND..NOT.LSIF).OR.MODEL_GW>0.OR.AEXX==0) THEN
         FORHF=0; SIFHF=0
         PROFILING_STOP('fock_force')
         RETURN
      ENDIF

#ifdef _OPENACC_TST
      PUSH_ACC_EXEC_ON(.TRUE.)
!$ACC WAIT
      ACC_ASYNC_Q=ACC_ASYNC_SYNC

      CALL ACC_COPYIN_TYPED_VAR(W)
      CALL ACC_COPYIN_TYPED_VAR(NONLR_S)
      CALL ACC_COPYIN_TYPED_VAR(NONL_S)
      CALL ACC_COPYIN_TYPED_VAR(P)
      CALL ACC_COPYIN_TYPED_VAR(WDES_FOCK)
      CALL ACC_COPYIN_TYPED_VAR(AUG_DES)
      CALL ACC_COPYIN_TYPED_VAR(FAST_AUG_FOCK)
!$ACC ENTER DATA COPYIN(LATT_CUR,TRANS_MATRIX_FOCK) __IF_ASYNC__
#endif

      CALL CHECK_FULL_KPOINTS
      NULLIFY(ROT_HANDLE)

      WHF=W
      WHF%WDES=>WDES_FOCK
#ifdef _OPENACC
      CALL ACC_COPYIN_TYPED_VAR(WHF)
#endif

      NDIR=3; IF (LSIF) NDIR=9

      ! determine the number of (partially) occupied states
      NB_TOT=0
      DO ISP=1,WHF%WDES%ISPIN
         DO NK=1,WHF%WDES%NKPTS
            DO NB=1,WHF%WDES%NB_TOT
               IF (ABS(W%FERTOT(NB,NK,ISP))>1E-8_q) NB_TOT=MAX(NB_TOT,NB)
            ENDDO
         ENDDO
      ENDDO

      NBLK=NBLOCK_FOCK
      NBLK=MIN(NBLK,NB_TOT)

      CALL WRK_ALLOCATE
!$ACC WAIT IF(ACC_EXEC_ON)

!$ACC KERNELS PRESENT(FORHF,SIF) __IF_ASYNC__
      FORHF=0; SIF=0
!$ACC END KERNELS

      spn: DO ISP=1,WHF%WDES%ISPIN
      kpt: DO NK=1,WHF%WDES%NKPTS

         ! set all NK dependent stuff
         CALL PREAMBLE_K ; IF (LSKIP) CYCLE kpt

         ! run over all bands at NK in chunks of size NBLOCK_ACT
         band: DO NPOS=1,NB_TOTK,NBLOCK
            NBLOCK_ACT=MIN(NB_TOTK-NPOS+1,NBLOCK)

#ifdef _OPENACC
            DO_BATCHED_FFT=(NBLOCK_ACT==NBLOCK).AND.ACC_EXEC_ON
#endif
            CALL FFT_AND_GATHER
            qpt: DO NQ=1,KPOINTS_FULL%NKPTS

               ! set all NQ dependent stuff
               CALL PREAMBLE_Q ; IF (LSKIP) CYCLE qpt

#ifndef _OPENACC
!$            CALL SET_ROT_HANDLE(P,LATT_CUR,WDESQ,ROT_HANDLE)

!$OMP PARALLEL PRIVATE(WQ) REDUCTION(+:SIF) REDUCTION(-:FORHF)
!$            CALL NEWWAV(WQ,WDESQ,.TRUE.)
!$OMP DO SCHEDULE(STATIC) FIRSTPRIVATE(GWORK,CRHOLM,AUG_DES,POTFAK,W1,CDLM0,CDLM,CDIJ) &
!$OMP PRIVATE(MQ,LSHIFT,N,NGLB,WEIGHT,IDIR,ENL,ISPINOR,NI,NT,LMMAXC,NPRO,NIP,RTMP)
#endif
               mband: DO MQ=1,WHF%WDES%NBANDS
                  IF (ABS(WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ))<=1E-10_q .OR. &
                      (MQ-1)*W%WDES%NB_PAR+W%WDES%NB_LOW<NBANDSGWLOW_FOCK) CYCLE mband

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

                     CALL W1_ROTATE_AND_FFT(WQ, ELEMENT(WHF, WDESQ_IRZ, MQ, ISP_IRZ), ROT_HANDLE, P, LATT_CUR, LSHIFT)
#endif
                  ENDIF

                  ! calculate charge phi_q nq(r) phi_k nk(r)
#if defined(_OPENACC) && defined(gammareal)
                  GWORK=>GWORK_FWD
#endif
                  CALL FOCK_CHARGE_MU(WIN(1:NBLOCK_ACT),WQ,GWORK,CRHOLM)

#ifdef _OPENACC
                  IF (DO_BATCHED_FFT) THEN
                     PROFILING_START('batched_fft')
                     ! fft to reciprocal space
#ifdef gammareal
!$ACC WAIT(ACC_ASYNC_Q) ASYNC(ACC_ASYNC_Q+1)
                     DO N=1,NBLOCK
                        NS=MOD(N,2)+1
                        ierr = cufftSetStream(batchedPlanD2Z(NS),acc_get_cuda_stream(ACC_ASYNC_Q+NS-1))
!$ACC HOST_DATA USE_DEVICE(GWORK_FWD,GWORK_BCK)
                        ierr = ierr + cufftExecD2Z(batchedPlanD2Z(NS),GWORK_FWD(1,N),GWORK_BCK(1,N))
!$ACC END HOST_DATA
                     ENDDO
!$ACC WAIT(ACC_ASYNC_Q+1) ASYNC(ACC_ASYNC_Q)
                     GWORK=>GWORK_BCK
#else
                     ierr = cufftSetStream(batchedPlan,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(GWORK)
                     ierr = ierr + cufftExecZ2Z(batchedPlan,GWORK,GWORK,CUFFT_FORWARD)
!$ACC END HOST_DATA
#endif
                     PROFILING_STOP('batched_fft')
                  ELSE
                     DO N=1,NBLOCK_ACT
                        ! fft to reciprocal space
                        CALL FFT3D(GWORK(1,N),GRIDHF,-1)
                     ENDDO
                  ENDIF

                  CALL CALC_GFAC_MU_WEIGHTS

                  IF (LSIF) CALL APPLY_GFAC_DER_MU(GRIDHF,GWORK(1,1),POTFAK(1,0),SIF(0),GFAC_MU_WEIGHTS,NBLOCK_ACT)

                  IF (WHF%WDES%LOVERL) THEN
                     ! multiply by 4 pi e^2/G^2 and divide by # of gridpoints to obtain potential
                     CALL APPLY_GFAC_WEIGHT_MU(GRIDHF,GWORK(1,1),POTFAK(1,0),GFAC_MU_WEIGHTS,NBLOCK_ACT)

                     IF (DO_BATCHED_FFT) THEN
                        ! back to real space to get  \int phi_q(r) phi_k(r) / (r-r') d3r
                        PROFILING_START('batched_fft')
#ifdef gammareal
!$ACC WAIT(ACC_ASYNC_Q) ASYNC(ACC_ASYNC_Q+1)
                        DO N=1,NBLOCK_ACT
                           NS=MOD(N,2)+1
                           ierr = cufftSetStream(batchedPlanZ2D(NS),acc_get_cuda_stream(ACC_ASYNC_Q+NS-1))
!$ACC HOST_DATA USE_DEVICE(GWORK_FWD,GWORK_BCK)
                           ierr = ierr + cufftExecZ2D(batchedPlanZ2D(NS),GWORK_BCK(1,N),GWORK_FWD(1,N))
!$ACC END HOST_DATA
                        ENDDO
!$ACC WAIT(ACC_ASYNC_Q+1) ASYNC(ACC_ASYNC_Q)
                        GWORK=>GWORK_FWD
#else
!$ACC HOST_DATA USE_DEVICE(GWORK)
                        ierr = ierr + cufftExecZ2Z(batchedPlan,GWORK,GWORK,CUFFT_INVERSE)
!$ACC END HOST_DATA
#endif
                        PROFILING_STOP('batched_fft')
                     ELSE
                        DO N=1,NBLOCK_ACT
                           ! back to real space to get  \int phi_q(r) phi_k(r) / (r-r') d3r
                           CALL FFT3D(GWORK(1,N),GRIDHF,1)
                        ENDDO
                     ENDIF
                  ENDIF
#else
                  nband: DO N=1,NBLOCK_ACT

                     NGLB=NPOS+N-1
                     ! fft to reciprocal space
                     CALL FFT3D(GWORK(1,N),GRIDHF,-1)

                     WEIGHT=WHF%WDES%RSPIN*WHF%WDES%WTKPT(NK)*WHF%FERTOT(NGLB,NK,ISP)*  &
                          WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ)*WEIGHT_Q

                     IF (LSIF) CALL APPLY_GFAC_DER(GRIDHF,GWORK(1,N),POTFAK(1,0),SIF(0),WEIGHT)

                     IF (WHF%WDES%LOVERL) THEN
                        ! multiply by 4 pi e^2/G^2 and divide by # of gridpoints to obtain potential
                        ! and multiply in the k-point and Fermi-weight (gK 19.06.2020)
                        ! this is cleaner and more transparent
                        CALL APPLY_GFAC_WEIGHT(GRIDHF,GWORK(1,N),POTFAK(1,0),WEIGHT)
                        ! back to real space to get  \int phi_q(r) phi_k(r) / (r-r') d3r
                        CALL FFT3D(GWORK(1,N),GRIDHF,1)
                     ENDIF
                  ENDDO nband
#endif
                  IF (WHF%WDES%LOVERL) THEN
                     ! multiplicative factor used in RPROMU_HF_*
                     AUG_DES%RINPL=1.0_q/GRIDHF%NPLWV
!$ACC UPDATE DEVICE(AUG_DES%RINPL) __IF_ASYNC__
#ifdef _OPENACC
                     IF (ACC_EXEC_ON) THEN
                        CALL RPROMU_HF_ACC(FAST_AUG_FOCK,AUG_DES,CDLM0,SIZE(CDLM0,1),NBLOCK_ACT,GWORK(1,1),SIZE(GWORK,1))
                     ELSE
#endif
                        ! workspace used by RPROMU_HF
                        DO N=1,NBLOCK_ACT
                           W1(N)%CPROJ => CDLM0(:,N)
                        ENDDO
                        CALL RPROMU_HF(FAST_AUG_FOCK,AUG_DES,W1,NBLOCK_ACT,GWORK(1,1),SIZE(GWORK,1))
#ifdef _OPENACC
                     ENDIF
#endif
                     IF (WHF%WDES%NRSPINORS==2) THEN
                        DO N=1,NBLOCK_ACT
                           CALL __GCOPY__(AUG_DES%NPRO,CDLM0(1,N),1,CDLM0(AUG_DES%NPRO+1,N),1)
                        ENDDO
                     ENDIF

                     DO IDIR=1,NDIR
#ifdef _OPENACC
                        IF (ACC_EXEC_ON) THEN
                           CALL RPROMU_HF_ACC(FAST_AUG(IDIR),AUG_DES,CDLM,SIZE(CDLM,1),NBLOCK_ACT,GWORK(1,1),SIZE(GWORK,1))
                        ELSE
#endif
                           ! workspace used by RPROMU_HF
                           DO N=1,NBLOCK_ACT
                              W1(N)%CPROJ => CDLM(:,N)
                           ENDDO
                           CALL RPROMU_HF(FAST_AUG(IDIR),AUG_DES,W1,NBLOCK_ACT,GWORK(1,1),SIZE(GWORK,1))
#ifdef _OPENACC
                        ENDIF
#endif
                        IF (WHF%WDES%NRSPINORS==2) THEN
                           DO N=1,NBLOCK_ACT
                              CALL __GCOPY__(AUG_DES%NPRO,CDLM(1,N),1,CDLM(AUG_DES%NPRO+1,N),1)
                           ENDDO
                        ENDIF
#ifndef _OPENACC
                        ENL=0
                        DO N=1,NBLOCK_ACT

                           CALL CALC_DLLMM_TRANS(WHF%WDES,AUG_DES,TRANS_MATRIX_FOCK,CDIJ(:,:,:,:,1),CDLM0(:,N))

                           CALL CALC_DLLMM_TRANS(WHF%WDES,AUG_DES,TRANS_MATRIX_FOCK,CDIJ(:,:,:,:,2),CDLM(:,N))

                           NGLB=N+NPOS-1

                           DO ISPINOR=0,WHF%WDES%NRSPINORS-1
                              DO NI=1,WHF%WDES%NIONS
                                 NT=WHF%WDES%ITYP(NI)

                                 LMMAXC=WHF%WDES%LMMAX(NT)
                                 IF (LMMAXC==0) CYCLE

                                 NPRO=WHF%WDES%LMBASE(NI)+ISPINOR*WHF%WDES%NPRO/2

                                 CALL ECCP_NL_FOCK(LMDIM,LMMAXC,CDIJ(1,1,NI,1+ISPINOR,1),CDIJ(1,1,NI,1+ISPINOR,2), &
                                      WQ%CPROJ(NPRO+1),CPROJKXYZ(NPRO+1,N,IDIR),WIN(N)%CPROJ(NPRO+1),ENL(N,NI),1.0_q)
                              ENDDO
                           ENDDO

                        ENDDO
#else
#if PGI_BETWEEN_20_9_AND_21_3
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
#endif
                        CALL DLLMM_TRANS_ECCP_NL_FOCK_ACC
#endif
                        PROFILING_START('contraction')

                        IF (IDIR<=3) THEN
!$ACC PARALLEL LOOP GANG PRIVATE(RTMP,NIP) PRESENT(WHF,WHF%WDES,FORHF,ENL) __IF_ASYNC__
                           DO NI=1,WHF%WDES%NIONS
                              NIP=NI_GLOBAL(NI,WHF%WDES%COMM_INB)
                              RTMP=0
!$ACC LOOP VECTOR REDUCTION(+:RTMP)
                              DO N=1,NBLOCK_ACT
                                 RTMP=RTMP-ENL(N,NI)
                              ENDDO
!$ACC ATOMIC UPDATE
                              FORHF(IDIR,NIP)=FORHF(IDIR,NIP)+RTMP
                           ENDDO
                        ELSE
!$ACC PARALLEL LOOP GANG PRIVATE(RTMP) PRESENT(WHF,WHF%WDES,SIF,ENL) __IF_ASYNC__
                           DO NI=1,WHF%WDES%NIONS
                              RTMP=0
!$ACC LOOP VECTOR REDUCTION(+:RTMP)
                              DO N=1,NBLOCK_ACT
                                 RTMP=RTMP+ENL(N,NI)
                              ENDDO
!$ACC ATOMIC UPDATE
                              SIF(IDIR-3)=SIF(IDIR-3)+RTMP
                           ENDDO
                        ENDIF

                        PROFILING_STOP('contraction')

                     ENDDO
                  ENDIF

               ENDDO mband
#ifndef _OPENACC
!$OMP END DO
!$             CALL DELWAV(WQ,.TRUE.)
!$OMP END PARALLEL
#endif
            ENDDO qpt
         ENDDO band

      ENDDO kpt
      ENDDO spn

      CALLMPI( M_sum_d(WHF%WDES%COMM_KINTER,FORHF(1,1),NONLR_S%NIONS*3))
      CALLMPI( M_sum_d(WHF%WDES%COMM_KINTER,SIF,7))

      CALLMPI( M_sum_d(WDESK%COMM_KIN,FORHF(1,1),NONLR_S%NIONS*3))
      CALLMPI( M_sum_d(WDESK%COMM_KIN,SIF,7))

!$ACC UPDATE SELF(FORHF,SIF) __IF_ASYNC__
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)

      FORHF=-FORHF

      IDIR=0
      DO I=1,3
         DO J=1,I
            IDIR=IDIR+1
            SIFHF(I,J)=SIF(IDIR)
            SIFHF(J,I)=SIF(IDIR)
         ENDDO
      ENDDO

      CALL WRK_DEALLOCATE
#ifdef _OPENACC
      CALL ACC_DELETE_TYPED_VAR(WHF)
#endif
      CALL WNULLIFY(WHF)

#ifdef _OPENACC_TST
      CALL ACC_DELETE_TYPED_VAR(W)
      CALL ACC_DELETE_TYPED_VAR(NONLR_S)
      CALL ACC_DELETE_TYPED_VAR(NONL_S)
      CALL ACC_DELETE_TYPED_VAR(P)
      CALL ACC_DELETE_TYPED_VAR(AUG_DES)
      CALL ACC_DELETE_TYPED_VAR(WDES_FOCK)
      CALL ACC_DELETE_TYPED_VAR(FAST_AUG_FOCK)
      CALL ACC_DELETE_TYPED_VAR(KPOINTS_FULL)
!$ACC EXIT DATA DELETE(LATT_CUR,TRANS_MATRIX_FOCK) __IF_ASYNC__
!$ACC WAIT
      ACC_ASYNC_Q=ACC_ASYNC_SYNC
      POP_ACC_EXEC_ON
#endif

      PROFILING_STOP('fock_force')

!***********************************************************************
!***********************************************************************
!
! Internal subroutines: begin
!
!***********************************************************************
!***********************************************************************
      CONTAINS

!************************ SUBROUTINE WRK_ALLOCATE **********************
!
! Allocate workspace and setup the derivatives of the augmentation
! charges w.r.t. the ionic positions and the lattice vectors
! (if LSIF=.TRUE.).
!
!***********************************************************************

      SUBROUTINE WRK_ALLOCATE
      USE tutor, ONLY: vtutor, isAlert, FockForce
      ! local variables
      TYPE (latt) :: LATT_FIN1,LATT_FIN2

      REAL(q) :: DIS,DISPL1(3,NONLR_S%NIONS),DISPL2(3,NONLR_S%NIONS)
      INTEGER :: N,IDIR,I,J

      INTEGER :: ISTATUS,ISTT

      PROFILING_START('wrk_allocate')

#if defined(_OPENACC) && defined(gammareal)
      ALLOCATE(GWORK_FWD(m_ GRIDHF%MPLWV,NBLK),GWORK_BCK(m_ GRIDHF%MPLWV,NBLK),STAT=ISTATUS)
!$ACC ENTER DATA CREATE(GWORK_FWD,GWORK_BCK) __IF_ASYNC__
#else
      ALLOCATE(GWORK(m_ GRIDHF%MPLWV,NBLK),STAT=ISTATUS)
!$ACC ENTER DATA CREATE(GWORK) __IF_ASYNC__
#endif

      ALLOCATE(POTFAK(GRIDHF%MPLWV,0:NDIR-3),FSG(0:NDIR-3),STAT=ISTT)
      ISTATUS=ISTATUS+ISTT

!$ACC ENTER DATA CREATE(POTFAK,FORHF,SIF) __IF_ASYNC__

! average electrostatic potential prefactor for k=k' and n=n'
      IF (LSIF) THEN
         CALL SET_FSG_DER(GRIDHF,LATT_CUR,FSG)
      ELSE
         FSG(0)=SET_FSG(GRIDHF,LATT_CUR)
      ENDIF

!$ACC ENTER DATA CREATE(WDESQ) __IF_ASYNC__
      CALL SETWDES(WHF%WDES,WDESQ,0)
#if defined(_OPENACC) || !defined(_OPENMP)
!$ACC ENTER DATA CREATE(WQ)  __IF_ASYNC__
      CALL NEWWAV(WQ,WDESQ,.TRUE.,ISTT)
      ISTATUS=ISTATUS+ISTT
#endif

!$ACC ENTER DATA CREATE(WDESK) __IF_ASYNC__
      CALL SETWDES(WHF%WDES,WDESK,0)

      ALLOCATE(WIN(NBLK),W1(NBLK))
!$ACC ENTER DATA CREATE(WIN(:)) __IF_ASYNC__
      DO N=1,NBLK
         CALL NEWWAV(WIN(N),WDESK,.TRUE.,ISTT)
         IF (ISTT/=0) EXIT
      ENDDO
      ISTATUS=ISTATUS+ISTT

      IF (WHF%WDES%LOVERL) THEN
         ALLOCATE(CRHOLM(AUG_DES%NPROD*WHF%WDES%NRSPINORS,NBLK), &
                  CDLM0(AUG_DES%NPROD*WHF%WDES%NRSPINORS,NBLK), &
                  CDLM(AUG_DES%NPROD*WHF%WDES%NRSPINORS,NBLK), &
                  CPROJXYZ(WHF%WDES%NPROD,WHF%WDES%NBANDS,NDIR), &
                  CPROJKXYZ(WHF%WDES%NPROD,NBLK,NDIR), &
                  ENL(NBLK,NONLR_S%NIONS),STAT=ISTT)
         ISTATUS=ISTATUS+ISTT
#ifndef _OPENACC
         ALLOCATE(CDIJ(LMDIM,LMDIM,WHF%WDES%NIONS,WHF%WDES%NRSPINORS,2),STAT=ISTT)
         ISTATUS=ISTATUS+ISTT
#endif
         ALLOCATE(FAST_AUG(NDIR),STAT=ISTT)
!!$ACC ENTER DATA CREATE(FAST_AUG(:)) __IF_ASYNC__
         ISTATUS=ISTATUS+ISTT

! setup first derivative of augmentation charges with respect to ionic positions
         PUSH_ACC_EXEC_ON(.FALSE.)

         DIS=fd_displacement

         DO IDIR=1,3
            CALL COPY_FASTAUG(FAST_AUG_FOCK,FAST_AUG(IDIR))

            DISPL1=0
            DISPL1(IDIR,:)=-DIS

            DISPL2=0
            DISPL2(IDIR,:)= DIS

            CALL RSPHER_ALL(GRIDHF,FAST_AUG(IDIR),LATT_CUR,LATT_CUR,LATT_CUR,DISPL1,DISPL2,1)

            ! phase factor identical to FAST_AUG_FOCK
            IF (ASSOCIATED(FAST_AUG(IDIR)%CRREXP)) DEALLOCATE(FAST_AUG(IDIR)%CRREXP)
            FAST_AUG(IDIR)%CRREXP=>FAST_AUG_FOCK%CRREXP

#ifdef shmem_rproj
            IF (FAST_AUG(IDIR)%COMM_shmem%NODE_ME==1) FAST_AUG(IDIR)%RPROJ=FAST_AUG(IDIR)%RPROJ*SQRT(LATT_CUR%OMEGA)*(1._q/(2._q*DIS))
            CALLMPI( M_barrier(FAST_AUG(IDIR)%COMM_shmem))
#else
            FAST_AUG(IDIR)%RPROJ=FAST_AUG(IDIR)%RPROJ*SQRT(LATT_CUR%OMEGA)*(1._q/(2._q*DIS))
#endif
         ENDDO

! setup first derivative of augmentation charges with respect to lattice vectors
         IF (LSIF) THEN
            IDIR=3
            DO I=1,3
               DO J=1,I
                  IDIR=IDIR+1
                  LATT_FIN1%A=LATT_CUR%A
                  LATT_FIN2%A=LATT_CUR%A
                  LATT_FIN1%A(I,:)=LATT_CUR%A(I,:)+DIS*LATT_CUR%A(J,:)
                  LATT_FIN2%A(I,:)=LATT_CUR%A(I,:)-DIS*LATT_CUR%A(J,:)

                  CALL LATTIC(LATT_FIN1)
                  CALL LATTIC(LATT_FIN2)

                  CALL COPY_FASTAUG(FAST_AUG_FOCK,FAST_AUG(IDIR))

                  DISPL1=0
                  CALL RSPHER_ALL(GRIDHF,FAST_AUG(IDIR),LATT_FIN2,LATT_FIN1,LATT_CUR,DISPL1,DISPL1,1,LOMEGA=.TRUE.)

                  ! phase factor identical to FAST_AUG_FOCK
                  IF (ASSOCIATED(FAST_AUG(IDIR)%CRREXP)) DEALLOCATE(FAST_AUG(IDIR)%CRREXP)
                  FAST_AUG(IDIR)%CRREXP=>FAST_AUG_FOCK%CRREXP

#ifdef shmem_rproj
                 IF (FAST_AUG(IDIR)%COMM_shmem%NODE_ME==1) FAST_AUG(IDIR)%RPROJ=FAST_AUG(IDIR)%RPROJ*(1._q/(2._q*DIS))
                 CALLMPI( M_barrier(FAST_AUG(IDIR)%COMM_shmem))
#else
                 FAST_AUG(IDIR)%RPROJ=FAST_AUG(IDIR)%RPROJ*(1._q/(2._q*DIS))
#endif
               ENDDO
            ENDDO
         ENDIF

#ifdef _OPENACC
         POP_ACC_EXEC_ON
         DO IDIR=1,NDIR
            CALL ACC_COPYIN_TYPED_VAR(FAST_AUG(IDIR))
         ENDDO
!$ACC ENTER DATA CREATE(CRHOLM,CDLM0,CDLM,CPROJXYZ,CPROJKXYZ,ENL) __IF_ASYNC__
#endif
      ENDIF

      CALLMPI( M_sum_i(WHF%WDES%COMM,ISTATUS,1))
      IF (ISTATUS/=0) THEN
! return if there was not enough memory: the forces will be bogus
! but we might still end up with a usable WAVECAR file.
         CALL vtutor%write(isAlert, FockForce)
         PROFILING_STOP('fock_force')
         RETURN
      ENDIF

#ifdef _OPENACC
#ifdef gammareal
      ierr = cufftPlan3d(batchedPlanD2Z(1),GRIDHF%NGPTAR(3),GRIDHF%NGPTAR(2),GRIDHF%NGPTAR(1),CUFFT_D2Z)
      ierr = cufftPlan3d(batchedPlanD2Z(2),GRIDHF%NGPTAR(3),GRIDHF%NGPTAR(2),GRIDHF%NGPTAR(1),CUFFT_D2Z)
      ierr = cufftPlan3d(batchedPlanZ2D(1),GRIDHF%NGPTAR(3),GRIDHF%NGPTAR(2),GRIDHF%NGPTAR(1),CUFFT_Z2D)
      ierr = cufftPlan3d(batchedPlanZ2D(2),GRIDHF%NGPTAR(3),GRIDHF%NGPTAR(2),GRIDHF%NGPTAR(1),CUFFT_Z2D)
#else
      fftdims(1) = GRIDHF%NGPTAR(3)
      fftdims(2) = GRIDHF%NGPTAR(2)
      fftdims(3) = GRIDHF%NGPTAR(1)
      ierr = cufftPlanMany(batchedPlan,3,fftdims, &
                                         fftdims,1,m_ GRIDHF%MPLWV,&
                                         fftdims,1,m_ GRIDHF%MPLWV,CUFFT_Z2Z,NBLK)
#endif
      ALLOCATE(GFAC_MU_WEIGHTS(NBLK))
!$ACC ENTER DATA CREATE(GFAC_MU_WEIGHTS) __IF_ASYNC__
#endif
      PROFILING_STOP('wrk_allocate')

      RETURN
      END SUBROUTINE WRK_ALLOCATE


!************************ SUBROUTINE WRK_DEALLOCATE ********************
!
! Deallocate workspace
!
!***********************************************************************

      SUBROUTINE WRK_DEALLOCATE
      ! local variables
      INTEGER :: N,IDIR

      PROFILING_START('wrk_deallocate')

#ifdef _OPENACC
!$ACC EXIT DATA DELETE(GFAC_MU_WEIGHTS) __IF_ASYNC__
      DEALLOCATE(GFAC_MU_WEIGHTS)
#ifdef gammareal
      ierr = cufftDestroy(batchedPlanD2Z(1))
      ierr = cufftDestroy(batchedPlanD2Z(2))
      ierr = cufftDestroy(batchedPlanZ2D(1))
      ierr = cufftDestroy(batchedPlanZ2D(2))
#else
      ierr = ierr + cufftDestroy(batchedPlan)
#endif
#endif
      IF (WHF%WDES%LOVERL) THEN
         DO N=1,NBLK
            NULLIFY(W1(N)%CPROJ)
         ENDDO

!$ACC EXIT DATA DELETE(CRHOLM,CDLM0,CDLM,CPROJXYZ,CPROJKXYZ,ENL) __IF_ASYNC__
         DEALLOCATE(CRHOLM,CDLM0,CDLM,CPROJXYZ,CPROJKXYZ,ENL)
#ifndef _OPENACC
         DEALLOCATE(CDIJ)
#endif
#ifdef _OPENACC
         DO IDIR=1,NDIR
            CALL ACC_DELETE_TYPED_VAR(FAST_AUG(IDIR))
         ENDDO
#endif
         PUSH_ACC_EXEC_ON(.FALSE.)

         DO IDIR=1,NDIR
            NULLIFY(FAST_AUG(IDIR)%CRREXP)
            CALL NONLR_DEALLOC(FAST_AUG(IDIR))
         ENDDO

         POP_ACC_EXEC_ON

!!!$ACC EXIT DATA DELETE(FAST_AUG) __IF_ASYNC__
         DEALLOCATE(FAST_AUG)
      ENDIF
#ifdef _OPENACC
      CALL ACC_DELETE_TYPED_VAR(WDESQ)
      CALL ACC_DELETE_TYPED_VAR(WDESK)
#endif
#if defined(_OPENACC) || !defined(_OPENMP)
      CALL DELWAV(WQ,.TRUE.)
!$ACC EXIT DATA DELETE(WQ) __IF_ASYNC__
#endif
      DO N=1,NBLK
         CALL DELWAV(WIN(N),.TRUE.)
      ENDDO
!$ACC EXIT DATA DELETE(WIN(:)) __IF_ASYNC__
      DEALLOCATE(WIN,W1)

#if defined(_OPENACC) && defined(gammareal)
!$ACC EXIT DATA DELETE(GWORK_FWD,GWORK_BCK) __IF_ASYNC__
      DEALLOCATE(GWORK_FWD,GWORK_BCK)
      NULLIFY(GWORK)
#else
!$ACC EXIT DATA DELETE(GWORK) __IF_ASYNC__
      DEALLOCATE(GWORK)
#endif

!$ACC EXIT DATA DELETE(POTFAK,FORHF,SIF) __IF_ASYNC__
      DEALLOCATE(POTFAK)

      CALL DEALLOCATE_ROT_HANDLE(ROT_HANDLE)

      PROFILING_STOP('wrk_deallocate')

      RETURN
      END SUBROUTINE WRK_DEALLOCATE


!************************ SUBROUTINE PREAMBLE_K ************************
!
! Initialize all k-point dependent quantities
!
!***********************************************************************

      SUBROUTINE PREAMBLE_K

#ifdef MPI
      IF (MOD(NK-1,WHF%WDES%COMM_KINTER%NCPU).NE.WHF%WDES%COMM_KINTER%NODE_ME-1) THEN
         LSKIP=.TRUE. ; RETURN
      ENDIF
#endif
      LSKIP=.FALSE.

      NB_TOTK=MIN(NB_TOT, W%WDES%NB_TOTK(NK,ISP))
      NBLOCK =MIN(NB_TOTK,NBLK)

      CALL SETWDES(WHF%WDES,WDESK,NK)

! first derivative of wavefunction character with respect to all ionic positions
      IF (WHF%WDES%LOVERL) THEN
         IF (NONLR_S%LREAL) THEN
            CALL PHASER(W%WDES%GRID,LATT_CUR,NONLR_S,NK,W%WDES)
            CALL RPROXYZ(W%WDES%GRID,NONLR_S,P,LATT_CUR,W,W%WDES,ISP,NK,CPROJXYZ)
            IF (LSIF) CALL RPROLAT_DER(W%WDES%GRID,NONLR_S,P,LATT_CUR,W,W%WDES,ISP,NK,CPROJXYZ(1,1,4))
         ELSE
            CALL PHASE(W%WDES,NONL_S,NK)
            CALL PROJXYZ(NONL_S,W%WDES,W,LATT_CUR,ISP,NK,CPROJXYZ(:,:,1:3))
            IF (LSIF) CALL PROJLAT_DER(P,NONL_S,W%WDES,W,LATT_CUR,ISP,NK,CPROJXYZ(1,1,4))
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE PREAMBLE_K
#define nonblocking_bcast
#ifndef nonblocking_bcast
!************************ SUBROUTINE FFT_AND_GATHER ********************
!
!***********************************************************************
      SUBROUTINE FFT_AND_GATHER
      ! local variables
      INTEGER :: NI,N,NB_LOCAL,IDIR,NP

      PROFILING_START('fft_and_gather')

#ifdef _OPENACC
      CALL ACC_SET_ASYNC_Q(ACC_ASYNC_FFT_AND_GATHER_Q)
#endif

      DO N=NPOS,NPOS+NBLOCK_ACT-1
         IF (MOD(N-1,WHF%WDES%NB_PAR)+1==WHF%WDES%NB_LOW) THEN
            NI=N-NPOS+1 ; NB_LOCAL=1+(N-1)/WHF%WDES%NB_PAR
            CALL W1_COPY(ELEMENT(WHF,WIN(NI)%WDES1,NB_LOCAL,ISP),WIN(NI))
            CALL FFTWAV_W1(WIN(NI))
         ENDIF
      ENDDO

! copy the derivatives of the wave function characters into CPROJK array
      IF (WHF%WDES%LOVERL) THEN
! copy the derivatives of the wave function characters into CPROJKXYZ
!$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(WHF,CPROJKXYZ,CPROJXYZ) PRIVATE(NI,NB_LOCAL) __IF_ASYNC__
         DO IDIR=1,NDIR
            DO N=NPOS,NPOS+NBLOCK_ACT-1
               IF (MOD(N-1,WHF%WDES%NB_PAR)+1==WHF%WDES%NB_LOW) THEN
                  NI=N-NPOS+1 ; NB_LOCAL=1+(N-1)/WHF%WDES%NB_PAR
!$ACC LOOP VECTOR
                  DO NP=1,WHF%WDES%NPROD
                     CPROJKXYZ(NP,NI,IDIR)=CPROJXYZ(NP,NB_LOCAL,IDIR)
                  ENDDO
               ENDIF
           ENDDO
        ENDDO
      ENDIF
#ifdef MPI
! distribute WIN and CPROJKXYZ to all nodes
      redis: IF (WHF%WDES%COMM_INTER%NCPU>1) THEN
         DO N=NPOS,NPOS+NBLOCK_ACT-1
            NI=N-NPOS+1

            CALL M_bcast_z_from(WHF%WDES%COMM_INTER,WIN(NI)%CR(1), &
           &     SIZE(WIN(NI)%CR),MOD(N-1,WHF%WDES%NB_PAR)+1)

            IF (WHF%WDES%LOVERL) THEN
#ifndef gammareal
               CALL M_bcast_z_from(WHF%WDES%COMM_INTER,WIN(NI)%CPROJ(1), &
              &     SIZE(WIN(NI)%CPROJ),MOD(N-1,WHF%WDES%NB_PAR)+1)

               DO IDIR=1,NDIR
                  CALL M_bcast_z_from(WHF%WDES%COMM_INTER,CPROJKXYZ(:,NI,IDIR), &
                 &     SIZE(CPROJKXYZ,1),MOD(N-1,WHF%WDES%NB_PAR)+1)
               ENDDO
#else
               CALL M_bcast_d_from(WHF%WDES%COMM_INTER,WIN(NI)%CPROJ(1), &
              &     SIZE(WIN(NI)%CPROJ),MOD(N-1,WHF%WDES%NB_PAR)+1)

               DO IDIR=1,NDIR
                  CALL M_bcast_d_from(WHF%WDES%COMM_INTER,CPROJKXYZ(:,NI,IDIR), &
                 &     SIZE(CPROJKXYZ,1),MOD(N-1,WHF%WDES%NB_PAR)+1)
               ENDDO
#endif
            ENDIF
         ENDDO
      ENDIF redis
#endif
#ifdef _OPENACC
!$ACC WAIT(ACC_ASYNC_FFT_AND_GATHER_Q) IF(ACC_EXEC_ON)
      ACC_ASYNC_Q=ACC_ASYNC_SYNC
#endif
      PROFILING_STOP('fft_and_gather')

      RETURN
      END SUBROUTINE FFT_AND_GATHER
#else
!!#define interleave_communication
#ifdef interleave_communication
!************************ SUBROUTINE FFT_AND_GATHER ********************
!
!***********************************************************************
      SUBROUTINE FFT_AND_GATHER
#ifdef USENCCL
      USE nccl2for
#endif
      ! local variables
      INTEGER :: NI,N,NB_LOCAL,IDIR,NP
#ifdef MPI
      INTEGER, ALLOCATABLE :: requests(:)
      INTEGER :: nrequests
#endif
#ifdef USENCCL
      TYPE (ncclResult) :: ncclRes
#endif

      PROFILING_START('fft_and_gather')

#ifdef _OPENACC
      CALL ACC_SET_ASYNC_Q(ACC_ASYNC_FFT_AND_GATHER_Q)
#endif
#ifdef MPI
      IF (WHF%WDES%COMM_INTER%NCPU>1) THEN
         ALLOCATE(requests((NDIR+2)*NBLOCK_ACT))
         nrequests=0
#ifdef USENCCL
!$ACC WAIT(ACC_ASYNC_FFT_AND_GATHER_Q) IF(ACC_EXEC_ON)
         ncclRes = ncclGroupStart()
#endif
      ENDIF
#endif

      DO N=NPOS,NPOS+NBLOCK_ACT-1
         NI=N-NPOS+1 ; NB_LOCAL=1+(N-1)/WHF%WDES%NB_PAR

         IF (MOD(N-1,WHF%WDES%NB_PAR)+1==WHF%WDES%NB_LOW) THEN
            CALL W1_COPY(ELEMENT(WHF,WIN(NI)%WDES1,NB_LOCAL,ISP),WIN(NI))
            CALL FFTWAV_W1(WIN(NI))
         ENDIF
#ifdef MPI
         IF (WHF%WDES%COMM_INTER%NCPU>1) THEN
            nrequests=nrequests+1
            CALL M_ibcast_z_from(WHF%WDES%COMM_INTER,WIN(NI)%CR(1), &
           &     SIZE(WIN(NI)%CR),MOD(N-1,WHF%WDES%NB_PAR)+1,requests(nrequests))
         ENDIF
#endif

         IF (WHF%WDES%LOVERL) THEN
#ifdef MPI
            IF (WHF%WDES%COMM_INTER%NCPU>1) THEN
               nrequests=nrequests+1
#ifndef gammareal
               CALL M_ibcast_z_from(WHF%WDES%COMM_INTER,WIN(NI)%CPROJ(1), &
              &     SIZE(WIN(NI)%CPROJ),MOD(N-1,WHF%WDES%NB_PAR)+1,requests(nrequests))
#else
               CALL M_ibcast_d_from(WHF%WDES%COMM_INTER,WIN(NI)%CPROJ(1), &
              &     SIZE(WIN(NI)%CPROJ),MOD(N-1,WHF%WDES%NB_PAR)+1,requests(nrequests))
#endif
            ENDIF
#endif
            DO IDIR=1,NDIR
               IF (MOD(N-1,WHF%WDES%NB_PAR)+1==WHF%WDES%NB_LOW) THEN
!$ACC PARALLEL LOOP PRESENT(WHF,CPROJKXYZ,CPROJXYZ) __IF_ASYNC__
                  DO NP=1,WHF%WDES%NPROD
                     CPROJKXYZ(NP,NI,IDIR)=CPROJXYZ(NP,NB_LOCAL,IDIR)
                  ENDDO
               ENDIF
#ifdef MPI
               IF (WHF%WDES%COMM_INTER%NCPU>1) THEN
                  nrequests=nrequests+1
#ifndef gammareal
                  CALL M_ibcast_z_from(WHF%WDES%COMM_INTER,CPROJKXYZ(:,NI,IDIR), &
                 &     SIZE(CPROJKXYZ,1),MOD(N-1,WHF%WDES%NB_PAR)+1,requests(nrequests))
#else         
                  CALL M_ibcast_d_from(WHF%WDES%COMM_INTER,CPROJKXYZ(:,NI,IDIR), &
                 &     SIZE(CPROJKXYZ,1),MOD(N-1,WHF%WDES%NB_PAR)+1,requests(nrequests))
#endif
               ENDIF
#endif
            ENDDO
         ENDIF
      ENDDO
#ifdef MPI
      IF (WHF%WDES%COMM_INTER%NCPU>1) THEN
#if defined(USENCCL) && defined (_OPENACC)
         ncclRes = ncclGroupEnd()
         IF ( WHF%WDES%COMM_INTER%LUSENCCL ) THEN
!$ACC WAIT(ACC_ASYNC_FFT_AND_GATHER_Q) IF(ACC_EXEC_ON)
         ELSE
#endif
            CALL M_waitall(nrequests,requests(1))
#if defined(USENCCL) && defined (_OPENACC)
         ENDIF
#endif
         DEALLOCATE(requests)
      ENDIF
#endif
#ifdef _OPENACC
      ACC_ASYNC_Q=ACC_ASYNC_SYNC
#endif
      PROFILING_STOP('fft_and_gather')

      RETURN
      END SUBROUTINE FFT_AND_GATHER
#else
      SUBROUTINE FFT_AND_GATHER
      ! local variables
      INTEGER :: NI,N,NB_LOCAL,IDIR,NP
#ifdef MPI
      INTEGER, ALLOCATABLE :: requests(:)
      INTEGER :: nrequests
#endif
#ifdef USENCCL
      TYPE (ncclResult) :: ncclRes
#endif

      PROFILING_START('fft_and_gather')

#ifdef _OPENACC
      CALL ACC_SET_ASYNC_Q(ACC_ASYNC_FFT_AND_GATHER_Q)
#endif

      DO N=NPOS,NPOS+NBLOCK_ACT-1
         IF (MOD(N-1,WHF%WDES%NB_PAR)+1==WHF%WDES%NB_LOW) THEN
            NI=N-NPOS+1 ; NB_LOCAL=1+(N-1)/WHF%WDES%NB_PAR
            CALL W1_COPY(ELEMENT(WHF,WIN(NI)%WDES1,NB_LOCAL,ISP),WIN(NI))
            CALL FFTWAV_W1(WIN(NI))
         ENDIF
      ENDDO

! copy the derivatives of the wave function characters into CPROJK array
      IF (WHF%WDES%LOVERL) THEN
! copy the derivatives of the wave function characters into CPROJKXYZ
!$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(WHF,CPROJKXYZ,CPROJXYZ) PRIVATE(NI,NB_LOCAL) __IF_ASYNC__
         DO IDIR=1,NDIR
            DO N=NPOS,NPOS+NBLOCK_ACT-1
               IF (MOD(N-1,WHF%WDES%NB_PAR)+1==WHF%WDES%NB_LOW) THEN
                  NI=N-NPOS+1 ; NB_LOCAL=1+(N-1)/WHF%WDES%NB_PAR
!$ACC LOOP VECTOR
                  DO NP=1,WHF%WDES%NPROD
                     CPROJKXYZ(NP,NI,IDIR)=CPROJXYZ(NP,NB_LOCAL,IDIR)
                  ENDDO
               ENDIF
           ENDDO
        ENDDO
      ENDIF
#ifdef MPI
! distribute WIN and CPROJKXYZ to all nodes
      redis: IF (WHF%WDES%COMM_INTER%NCPU>1) THEN
         ALLOCATE(requests((NDIR+2)*NBLOCK_ACT))
         nrequests=0
#ifdef USENCCL
!!!$ACC WAIT(ACC_ASYNC_FFT_AND_GATHER_Q) IF(ACC_EXEC_ON)
         ncclRes = ncclGroupStart()
#endif
         DO N=NPOS,NPOS+NBLOCK_ACT-1
            NI=N-NPOS+1

            nrequests=nrequests+1
            CALL M_ibcast_z_from(WHF%WDES%COMM_INTER,WIN(NI)%CR(1), &
           &     SIZE(WIN(NI)%CR),MOD(N-1,WHF%WDES%NB_PAR)+1,requests(nrequests))

            IF (WHF%WDES%LOVERL) THEN
               nrequests=nrequests+1
#ifndef gammareal
               CALL M_ibcast_z_from(WHF%WDES%COMM_INTER,WIN(NI)%CPROJ(1), &
              &     SIZE(WIN(NI)%CPROJ),MOD(N-1,WHF%WDES%NB_PAR)+1,requests(nrequests))

               DO IDIR=1,NDIR
                  nrequests=nrequests+1
                  CALL M_ibcast_z_from(WHF%WDES%COMM_INTER,CPROJKXYZ(:,NI,IDIR), &
                 &     SIZE(CPROJKXYZ,1),MOD(N-1,WHF%WDES%NB_PAR)+1,requests(nrequests))
               ENDDO
#else
               CALL M_ibcast_d_from(WHF%WDES%COMM_INTER,WIN(NI)%CPROJ(1), &
              &     SIZE(WIN(NI)%CPROJ),MOD(N-1,WHF%WDES%NB_PAR)+1,requests(nrequests))

               DO IDIR=1,NDIR
                  nrequests=nrequests+1
                  CALL M_ibcast_d_from(WHF%WDES%COMM_INTER,CPROJKXYZ(:,NI,IDIR), &
                 &     SIZE(CPROJKXYZ,1),MOD(N-1,WHF%WDES%NB_PAR)+1,requests(nrequests))
               ENDDO
#endif
            ENDIF
         ENDDO
#if defined(USENCCL) && defined (_OPENACC)
         ncclRes = ncclGroupEnd()
         IF ( WHF%WDES%COMM_INTER%LUSENCCL ) THEN
!$ACC WAIT(ACC_ASYNC_FFT_AND_GATHER_Q) IF(ACC_EXEC_ON)
         ELSE
#endif
            CALL M_waitall(nrequests,requests(1))
#if defined(USENCCL) && defined (_OPENACC)
         ENDIF
#endif
         DEALLOCATE(requests)
      ENDIF redis
#endif
      PROFILING_STOP('fft_and_gather')

      RETURN
      END SUBROUTINE FFT_AND_GATHER
#endif
#endif
!************************ SUBROUTINE PREAMBLE_Q ************************
!
! Initialize all q-point dependent quantities
!
!***********************************************************************

      SUBROUTINE PREAMBLE_Q

      IF (KPOINTS_FULL%WTKPT(NQ)==0.OR.(HFKIDENT.AND.SKIP_THIS_KPOINT_IN_FOCK(WHF%WDES%VKPT(:,NQ))) .OR. &
          (.NOT.HFKIDENT.AND.SKIP_THIS_KPOINT_IN_FOCK(KPOINTS_FULL%VKPT(:,NQ)-WHF%WDES%VKPT(:,NK)))) THEN
          LSKIP=.TRUE. ; RETURN
      ENDIF

      IF (ALLOCATED(WEIGHT_K_POINT_PAIR_SMALL_GROUP).AND.LSYMGRAD) THEN
         IF (WEIGHT_K_POINT_PAIR_SMALL_GROUP(NK,NQ)==0) THEN
            LSKIP=.TRUE. ; RETURN
         ENDIF
         WEIGHT_Q=WEIGHT_K_POINT_PAIR_SMALL_GROUP(NK,NQ)
      ELSE
         WEIGHT_Q=1
      ENDIF

      LSKIP=.FALSE.

      CALL SETWDES(WHF%WDES,WDESQ,NQ)

      CALL SETWDES(WHF%WDES,WDESQ_IRZ,KPOINTS_FULL%NEQUIV(NQ))

      ISP_IRZ=ISP; IF (KPOINTS_FULL%SPINFLIP(NQ)==1) ISP_IRZ=3-ISP

      ! set POTFAK for this q and k point
      IF (LSIF) THEN
         CALL SET_GFAC_DER(GRIDHF,LATT_CUR,NK,NQ,FSG,POTFAK)
      ELSE
         CALL SET_GFAC(GRIDHF,LATT_CUR,NK,NQ,FSG(0),POTFAK(1,0))
      ENDIF

      RETURN
      END SUBROUTINE PREAMBLE_Q


!************************ SUBROUTINE CALC_GFAC_DER_MU ******************
!
!***********************************************************************
#ifdef _OPENACC
      SUBROUTINE CALC_GFAC_MU_WEIGHTS
      ! local variables
      INTEGER :: N,NGLB
      REAL(q) :: FD

      DO N=1,NBLOCK_ACT
         NGLB=N+NPOS-1

         FD=WHF%FERTOT(NGLB,NK,ISP)*WHF%FERWE(MQ,KPOINTS_FULL%NEQUIV(NQ),ISP_IRZ)
         GFAC_MU_WEIGHTS(N)=WHF%WDES%RSPIN*WHF%WDES%WTKPT(NK)*FD*WEIGHT_Q

      ENDDO
!$ACC UPDATE DEVICE(GFAC_MU_WEIGHTS) __IF_ASYNC__

      RETURN
      END SUBROUTINE CALC_GFAC_MU_WEIGHTS
#endif

!************************ SUBROUTINE DLLMM_TRANS_ECCP_NL_FOCK_ACC ******
!
!***********************************************************************
#ifdef _OPENACC
      SUBROUTINE DLLMM_TRANS_ECCP_NL_FOCK_ACC
      ! local variables
      GDEF :: CDIJ0,CDIJ,GTMP
      INTEGER :: ISPINOR,N,NI,NT,LMMAXC,NPRO,NAUG,L,LP,LM,NGLB

      PROFILING_START('dllmm_trans_eccp_nl_fock_acc')

!$ACC KERNELS PRESENT(ENL) __IF_ASYNC__
      ENL=0
!$ACC END KERNELS

      DO ISPINOR=0,WHF%WDES%NRSPINORS-1
!$ACC PARALLEL LOOP COLLAPSE(2) GANG __IF_ASYNC__ &
!$ACC PRESENT(WHF,AUG_DES,TRANS_MATRIX_FOCK,CDLM0,CDLM,WQ,CPROJKXYZ,WIN,ENL) &
!$ACC PRIVATE(NT,LMMAXC,NPRO,NAUG,GTMP,L,LP,NGLB)
         DO N=1,NBLOCK_ACT
         DO NI=1,WHF%WDES%NIONS
            NT=WHF%WDES%ITYP(NI)

            LMMAXC=WHF%WDES%LMMAX(NT)
            IF (LMMAXC==0) CYCLE

            NPRO=WHF%WDES%LMBASE(NI)+ISPINOR*WHF%WDES%NPRO/2
            NAUG=AUG_DES%LMBASE(NI)+ISPINOR*AUG_DES%NPRO

            GTMP=0
!$ACC LOOP COLLAPSE(2) VECTOR PRIVATE(CDIJ0,CDIJ,LM) REDUCTION(+:GTMP)
            DO L =1,LMMAXC
            DO LP=1,LMMAXC
               CDIJ0=0; CDIJ=0
#ifdef bug1
!$ACC LOOP VECTOR REDUCTION(+:CDIJ0,CDIJ)
#endif
               DO LM=1,AUG_DES%LMMAX(NT)
                  CDIJ0=CDIJ0+TRANS_MATRIX_FOCK(LP,L,LM,NT)*CDLM0(LM+NAUG,N)
                  CDIJ =CDIJ +TRANS_MATRIX_FOCK(LP,L,LM,NT)*CDLM (LM+NAUG,N)
               ENDDO

               GTMP=GTMP+WQ%CPROJ(LP+NPRO)*(CDIJ0*GCONJG(CPROJKXYZ(L+NPRO,N,IDIR))+CDIJ*GCONJG(WIN(N)%CPROJ(L+NPRO)))

            ENDDO
            ENDDO

            NGLB=NPOS+N-1
            ENL(N,NI)=ENL(N,NI)+REAL(GTMP,KIND=q)

         ENDDO
         ENDDO
      ENDDO

      PROFILING_STOP('dllmm_trans_eccp_nl_fock_acc')

      RETURN
      END SUBROUTINE DLLMM_TRANS_ECCP_NL_FOCK_ACC
#endif
!***********************************************************************
!***********************************************************************
!
! Internal subroutines: end
!
!***********************************************************************
!***********************************************************************

      END SUBROUTINE FOCK_FORCE

      END MODULE fock_frc
