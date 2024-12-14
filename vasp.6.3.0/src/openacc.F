#include "symbol.inc"
!***********************************************************************
!
!***********************************************************************

  MODULE mopenacc
#ifdef _OPENACC
  USE prec

  USE mopenacc_struct_def

  USE mgrid_struct_def
  USE wave_struct_def
  USE nonlr_struct_def
  USE nonl_struct_def
  USE pseudo_struct_def
  USE mkpoints_struct_def

  USE string, ONLY: str
  USE tutor, ONLY: vtutor

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ACC_EXEC_ALLOW,ACC_EXEC_ON,ACC_PUSH_EXEC_ON,ACC_POP_EXEC_ON
  PUBLIC :: ACC_ASYNC_Q,ACC_ASYNC_SYNC,ACC_ASYNC_ASYNC,ACC_SET_ASYNC_Q,ACC_SYNC_ASYNC_Q
  PUBLIC :: ACC_COPYIN_TYPED_VAR,ACC_CREATE_TYPED_VAR,ACC_DELETE_TYPED_VAR,ACC_COPYOUT_TYPED_VAR
  PUBLIC :: ACC_IS_PRESENT
#ifdef USENCCL
  PUBLIC :: LUSENCCL
#endif

  PUBLIC :: INIT_ACC

  INTEGER :: ACC_MAX_STREAMS=1024

  INTEGER, PARAMETER :: MAXLEVEL=20
  INTEGER :: ACC_EXEC_ON_LEVEL=0
  LOGICAL :: ACC_EXEC_ON_STACK(MAXLEVEL)=.FALSE.

#ifdef ACC_DEBUG
  LOGICAL, PUBLIC :: ACC_DISABLE_SET_ASYNC_Q=.FALSE.
#endif

  INTERFACE ACC_COPYIN_TYPED_VAR
    MODULE PROCEDURE COPYIN_NONL_STRUCT
    MODULE PROCEDURE COPYIN_NONLR_STRUCT
    MODULE PROCEDURE COPYIN_NONLR_PROJ
    MODULE PROCEDURE COPYIN_LAYOUT
    MODULE PROCEDURE COPYIN_GRID_MAP
    MODULE PROCEDURE COPYIN_GRID_3D
    MODULE PROCEDURE COPYIN_TRANSIT
    MODULE PROCEDURE COPYIN_WAVEDES
    MODULE PROCEDURE COPYIN_WAVEDES1
    MODULE PROCEDURE COPYIN_WAVESPIN
    MODULE PROCEDURE COPYIN_WAVEFUNA
    MODULE PROCEDURE COPYIN_WAVEFUN1
    MODULE PROCEDURE COPYIN_POTCARS
    MODULE PROCEDURE COPYIN_KPOINTS_FULL
  END INTERFACE ACC_COPYIN_TYPED_VAR

  INTERFACE ACC_CREATE_TYPED_VAR
    MODULE PROCEDURE CREATE_WAVESPIN
  END INTERFACE ACC_CREATE_TYPED_VAR

  INTERFACE ACC_DELETE_TYPED_VAR
    MODULE PROCEDURE DELETE_NONL_STRUCT
    MODULE PROCEDURE DELETE_NONLR_STRUCT
    MODULE PROCEDURE DELETE_NONLR_PROJ
    MODULE PROCEDURE DELETE_LAYOUT
    MODULE PROCEDURE DELETE_GRID_MAP
    MODULE PROCEDURE DELETE_GRID_3D
    MODULE PROCEDURE DELETE_TRANSIT
    MODULE PROCEDURE DELETE_WAVEDES
    MODULE PROCEDURE DELETE_WAVEDES1
    MODULE PROCEDURE DELETE_WAVESPIN
    MODULE PROCEDURE DELETE_WAVEFUNA
    MODULE PROCEDURE DELETE_WAVEFUN1
    MODULE PROCEDURE DELETE_POTCARS
    MODULE PROCEDURE DELETE_KPOINTS_FULL
  END INTERFACE ACC_DELETE_TYPED_VAR

  INTERFACE ACC_COPYOUT_TYPED_VAR
    MODULE PROCEDURE COPYOUT_WAVESPIN
  END INTERFACE ACC_COPYOUT_TYPED_VAR

  INTERFACE ACC_SYNC_ASYNC_Q
    MODULE PROCEDURE SYNC_ASYNC_Q_ALL
    MODULE PROCEDURE SYNC_ASYNC_Q_RANGE
    MODULE PROCEDURE SYNC_ASYNC_Q_MASKED
  END INTERFACE ACC_SYNC_ASYNC_Q

  CONTAINS

!****************** SUBROUTINE INIT_ACC ********************************
!
!***********************************************************************

  SUBROUTINE INIT_ACC(IO)
    USE base
    USE main_mpi
!$  USE openmp
    USE openacc, ONLY : ACC_GET_DEVICE_TYPE,ACC_DEVICE_NVIDIA,ACC_DEVICE_HOST, &
                        ACC_GET_NUM_DEVICES,ACC_SET_DEVICE_NUM, &
                        ACC_INIT

    TYPE(in_struct) IO
    ! local variables
    INTEGER :: DEVICE_NUM,NUM_DEVICES,DEVICE_TYPE
    INTEGER :: NBLOCK_RANKS_PER_DEVICE=1

    PROFILING_START('init_acc')

    DEVICE_NUM=0

    DEVICE_TYPE=ACC_GET_DEVICE_TYPE()

    IF (DEVICE_TYPE==ACC_DEVICE_HOST) THEN
       ACC_EXEC_ALLOW=.FALSE.
       PROFILING_STOP('init_acc')
       RETURN
    ENDIF

#ifdef MPI
    IF (.NOT. CUDA_AWARE_SUPPORT()) CALL vtutor%error("INIT_ACC: your MPI is not CUDA-aware, STOPPING")
#endif

    ACC_EXEC_ALLOW=.TRUE.

    IF (DEVICE_TYPE==ACC_DEVICE_NVIDIA) THEN
       NUM_DEVICES=ACC_GET_NUM_DEVICES(ACC_DEVICE_NVIDIA)
#ifdef MPI
       DEVICE_NUM=MOD((COMM%NODE_ME-1)/NBLOCK_RANKS_PER_DEVICE,NUM_DEVICES)
#endif
!$OMP PARALLEL
       CALL ACC_SET_DEVICE_NUM(DEVICE_NUM,ACC_DEVICE_NVIDIA)
!$OMP END PARALLEL
    ENDIF

!$OMP PARALLEL
    CALL ACC_INIT(DEVICE_TYPE)
!$OMP END PARALLEL

    IF (IO%IU0>=0) WRITE(IO%IU0,'(" OpenACC runtime initialized ... ",I4," GPUs detected")') NUM_DEVICES

#if defined(MPI) && defined(USENCCL)
    CALL M_init_nccl(COMM_INTER)
    CALL M_init_nccl(COMM_INB)
    CALL M_init_nccl(COMM_KIN)
    CALL M_init_nccl(COMM_KINTER)
#endif

#ifdef debug
    CALL ACC_SET_WAIT_ON_EXCEPTION(.TRUE.)
#endif

    PROFILING_STOP('init_acc')
#ifdef MPI
    CONTAINS
       LOGICAL FUNCTION CUDA_AWARE_SUPPORT()
       USE iso_c_binding

       INTERFACE
         INTEGER(c_int) FUNCTION MPIX_Query_cuda_support() BIND(C, name="MPIX_Query_cuda_support")
         END FUNCTION
       END INTERFACE

       CHARACTER(LEN=1) :: ENVVAR_VALUE
       INTEGER :: ENVVAR_STAT

       ! This should tell us if MPI is CUDA-aware
       CUDA_AWARE_SUPPORT = MPIX_Query_cuda_support() == 1
       ! However, for OpenMPI some env variables can still deactivate it even though the previous
       ! check was positive
       CALL GET_ENVIRONMENT_VARIABLE("OMPI_MCA_mpi_cuda_support", ENVVAR_VALUE, STATUS=ENVVAR_STAT)
       IF (ENVVAR_STAT==0 .AND. ENVVAR_VALUE=='0') CUDA_AWARE_SUPPORT = .FALSE.
       CALL GET_ENVIRONMENT_VARIABLE("OMPI_MCA_opal_cuda_support", ENVVAR_VALUE, STATUS=ENVVAR_STAT)
       IF (ENVVAR_STAT==0 .AND. ENVVAR_VALUE=='0') CUDA_AWARE_SUPPORT = .FALSE.
       ! Just in case we might be non-OpenMPI, and their MPIX_Query_cuda_support behaves similarly
       CALL GET_ENVIRONMENT_VARIABLE("MV2_USE_CUDA", ENVVAR_VALUE, STATUS=ENVVAR_STAT)
       IF (ENVVAR_STAT==0 .AND. ENVVAR_VALUE=='0') CUDA_AWARE_SUPPORT = .FALSE.
       CALL GET_ENVIRONMENT_VARIABLE("MPICH_RDMA_ENABLED_CUDA", ENVVAR_VALUE, STATUS=ENVVAR_STAT)
       IF (ENVVAR_STAT==0 .AND. ENVVAR_VALUE=='0') CUDA_AWARE_SUPPORT = .FALSE.
       CALL GET_ENVIRONMENT_VARIABLE("PMPI_GPU_AWARE", ENVVAR_VALUE, STATUS=ENVVAR_STAT)
       IF (ENVVAR_STAT==0) CUDA_AWARE_SUPPORT =(ENVVAR_VALUE == '1')

       RETURN
       END FUNCTION CUDA_AWARE_SUPPORT
#endif
  END SUBROUTINE INIT_ACC


#ifdef debug
  SUBROUTINE SLEEP_ON_EXCEPTION()
     INTERFACE
     FUNCTION GETPID() BIND(C, NAME='getpid')
        USE ISO_C_BINDING
        INTEGER(C_INT) :: GETPID
     END FUNCTION
     END INTERFACE

        PRINT *, "EXCEPTION ENCOUNTERED! waiting for the debugger to be attached to PID ", getpid()
        DO WHILE (.TRUE.)
           CALL SLEEP(1)
        END DO
  END SUBROUTINE SLEEP_ON_EXCEPTION

  SUBROUTINE ACC_SET_WAIT_ON_EXCEPTION(ENABLE)
    USE ISO_C_BINDING
    IMPLICIT NONE
    LOGICAL :: ENABLE
    INTEGER :: IERR
    INTERFACE
    FUNCTION ATEXIT(FPTR) BIND(C, NAME='atexit')
       USE ISO_C_BINDING
       TYPE(C_FUNPTR), VALUE, INTENT(IN) :: FPTR
       INTEGER(c_int) :: ATEXIT
    END FUNCTION
    END INTERFACE

    TYPE(C_FUNPTR) :: fptr
    fptr = c_funloc(SLEEP_ON_EXCEPTION)
!$OMP PARALLEL
    IF ( ENABLE ) IERR = ATEXIT(fptr)
!$OMP END PARALLEL

  END SUBROUTINE ACC_SET_WAIT_ON_EXCEPTION
#endif
!****************** SUBROUTINE COPYIN_NONL_STRUCT **********************
!
!***********************************************************************

  SUBROUTINE COPYIN_NONL_STRUCT(VAR)
    TYPE(nonl_struct) :: VAR
    ! early exit if possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%NITYP, &
!$ACC&                  VAR%ITYP, &
!$ACC&                  VAR%LMMAX, &
!$ACC&                  VAR%LMBASE, &
!$ACC&                  VAR%POSION, &
!$ACC&                  VAR%VKPT_SHIFT, &
!$ACC&                  VAR%CQFAK, &
!$ACC&                  VAR%QPROJ, &
!$ACC&                  VAR%CREXP) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE COPYIN_NONL_STRUCT


!****************** SUBROUTINE COPYIN_NONLR_STRUCT *********************
!
!***********************************************************************

  SUBROUTINE COPYIN_NONLR_STRUCT(VAR)
    TYPE(nonlr_struct) :: VAR
    ! local variables
    INTEGER :: I
    ! early exit if possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%NITYP, &
!$ACC&                  VAR%ITYP, &
!$ACC&                  VAR%LMAX, &
!$ACC&                  VAR%LMMAX, &
!$ACC&                  VAR%CHANNELS, &
!$ACC&                  VAR%PSRMAX, &
!$ACC&                  VAR%RSMOOTH, &
!$ACC&                  VAR%POSION, &
!$ACC&                  VAR%VKPT_SHIFT, &
!$ACC&                  VAR%NLIMAX, &
!$ACC&                  VAR%NLI, &
!$ACC&                  VAR%LMBASE, &
!$ACC&                  VAR%NLIBASE, &
!$ACC&                  VAR%RPROJ, &
!$ACC&                  VAR%CRREXP) ASYNC(ACC_ASYNC_Q)
    ! derived type pointers
    DO I=1,SIZE(VAR%BETA)
       CALL ACC_COPYIN_TYPED_VAR(VAR%BETA(I))
    ENDDO
  END SUBROUTINE COPYIN_NONLR_STRUCT


!****************** SUBROUTINE COPYIN_NONLR_PROJ ***********************
!
!***********************************************************************

  SUBROUTINE COPYIN_NONLR_PROJ(VAR)
    TYPE(nonlr_proj) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%PSPRNL, &
!$ACC&                  VAR%LPS) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE COPYIN_NONLR_PROJ


!****************** SUBROUTINE COPYIN_LAYOUT ***************************
!
!***********************************************************************

  SUBROUTINE COPYIN_LAYOUT(VAR)
    TYPE(layout) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%I2, &
!$ACC&                  VAR%I3, &
!$ACC&                  VAR%INDEX) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE COPYIN_LAYOUT


!****************** SUBROUTINE COPYIN_GRID_MAP *************************
!
!***********************************************************************

  SUBROUTINE COPYIN_GRID_MAP(VAR)
    TYPE(grid_map) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%N, &
!$ACC&                  VAR%PTR, &
!$ACC&                  VAR%RMT, &
!$ACC&                  VAR%TBL, &
!$ACC&                  VAR%NI, &
!$ACC&                  VAR%PTRI, &
!$ACC&                  VAR%TBLI, &
!$ACC&                  VAR%RMTI) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE COPYIN_GRID_MAP


!****************** SUBROUTINE COPYIN_GRID_3D **************************
!
!***********************************************************************

  SUBROUTINE COPYIN_GRID_3D(VAR)
    TYPE(grid_3d) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%LPCTX, &
!$ACC&                  VAR%LPCTY, &
!$ACC&                  VAR%LPCTZ, &
!$ACC&                  VAR%LPCTX_, &
!$ACC&                  VAR%LPCTY_, &
!$ACC&                  VAR%LPCTZ_, &
!$ACC&                  VAR%FFTSCA, &
!$ACC&                  VAR%FFTWEIGHT, &
!$ACC&                  VAR%NINDPWCONJG, &
!$ACC&                  VAR%IND_IN_SPHERE, &
!$ACC&                  VAR%COMM, &
!$ACC&                  VAR%COMM_KIN, &
!$ACC&                  VAR%COMM_KINTER) ASYNC(ACC_ASYNC_Q)

    ! native type pointer if associated
    IF (ASSOCIATED(VAR%RL_FFT)) THEN
!!$ACC ENTER DATA COPYIN(VAR%RL_FFT) ASYNC(ACC_ASYNC_Q)
       CALL ACC_COPYIN_TYPED_VAR(VAR%RL_FFT)
       CALL ACC_ATTACH_ASYNC(VAR%RL_FFT,ACC_ASYNC_Q)
    ENDIF
    ! derived type pointers

    CALL ACC_COPYIN_TYPED_VAR(VAR%RC)
!!    CALL ACC_ATTACH_ASYNC(VAR%RC,ACC_ASYNC_Q)

    CALL ACC_COPYIN_TYPED_VAR(VAR%IN)
!!    CALL ACC_ATTACH_ASYNC(VAR%IN,ACC_ASYNC_Q)

    CALL ACC_COPYIN_TYPED_VAR(VAR%RL)
!!    CALL ACC_ATTACH_ASYNC(VAR%RL,ACC_ASYNC_Q)

    CALL ACC_COPYIN_TYPED_VAR(VAR%RC_IN)
!!    CALL ACC_ATTACH_ASYNC(VAR%RC_IN,ACC_ASYNC_Q)

    CALL ACC_COPYIN_TYPED_VAR(VAR%IN_RL)
!!    CALL ACC_ATTACH_ASYNC(VAR%IN_RL,ACC_ASYNC_Q)
  END SUBROUTINE COPYIN_GRID_3D


!****************** SUBROUTINE COPYIN_TRANSIT **************************
!
!***********************************************************************

  SUBROUTINE COPYIN_TRANSIT(VAR)
    TYPE(transit) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%IND1, &
!$ACC&                  VAR%INDCOL) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE COPYIN_TRANSIT


!****************** SUBROUTINE COPYIN_WAVEDES **************************
!
!***********************************************************************

  SUBROUTINE COPYIN_WAVEDES(VAR)
    TYPE(wavedes) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%NPLWKP, &
!$ACC&                  VAR%NGVECTOR, &
!$ACC&                  VAR%NGVECTOR_POS, &
!$ACC&                  VAR%NPLWKP_TOT, &
!$ACC&                  VAR%NB_TOTK, &
!$ACC&                  VAR%PL_INDEX, &
!$ACC&                  VAR%PL_COL, &
!$ACC&                  VAR%NPRO_POS, &
!$ACC&                  VAR%LMMAX, &
!$ACC&                  VAR%LMBASE, &
!$ACC&                  VAR%NITYP, &
!$ACC&                  VAR%ITYP, &
!$ACC&                  VAR%NT_GLOBAL, &
!$ACC&                  VAR%VKPT, &
!$ACC&                  VAR%WTKPT, &
!$ACC&                  VAR%NINDPW, &
!$ACC&                  VAR%LUSEINV, &
!$ACC&                  VAR%NINDPW_INV, &
!$ACC&                  VAR%FFTSCA, &
!$ACC&                  VAR%MAP_TO_FULL, &
!$ACC&                  VAR%IGX, &
!$ACC&                  VAR%IGY, &
!$ACC&                  VAR%IGZ, &
!$ACC&                  VAR%DATAKE, &
!$ACC&                  VAR%AT_GAMMA) ASYNC(ACC_ASYNC_Q)
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%COMM)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_INTER)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_INTER) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_INB)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_INB) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_KINTER)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_KINTER) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_KIN)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_KIN) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_SHMEM)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_SHMEM) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_intra_node)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_intra_node) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_inter_node)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_inter_node) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%GRID)) THEN
       CALL ACC_COPYIN_TYPED_VAR(VAR%GRID)
       CALL ACC_ATTACH_ASYNC(VAR%GRID,ACC_ASYNC_Q)
    ENDIF
  END SUBROUTINE COPYIN_WAVEDES


!****************** SUBROUTINE COPYIN_WAVEDES1 *************************
!
!***********************************************************************

  SUBROUTINE COPYIN_WAVEDES1(VAR)
    TYPE(wavedes1) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%PL_INDEX, &
!$ACC&                  VAR%PL_COL, &
!$ACC&                  VAR%NPRO_POS, &
!$ACC&                  VAR%LMMAX, &
!$ACC&                  VAR%LMBASE, &
!$ACC&                  VAR%NITYP, &
!$ACC&                  VAR%ITYP, &
!$ACC&                  VAR%NT_GLOBAL, &
!$ACC&                  VAR%NINDPW, &
!$ACC&                  VAR%NINDPW_INV, &
!$ACC&                  VAR%FFTSCA, &
!$ACC&                  VAR%MAP_TO_FULL, &
!$ACC&                  VAR%IGX, &
!$ACC&                  VAR%IGY, &
!$ACC&                  VAR%IGZ, &
!$ACC&                  VAR%VKPT, &
!$ACC&                  VAR%DATAKE) ASYNC(ACC_ASYNC_Q)
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%COMM)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_INTER)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_INTER) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_INB)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_INB) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_KINTER)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_KINTER) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_KIN)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_KIN) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_SHMEM)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_SHMEM) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_intra_node)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_intra_node) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_inter_node)) THEN
!$ACC ENTER DATA COPYIN(VAR%COMM_inter_node) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%GRID)) THEN
       CALL ACC_COPYIN_TYPED_VAR(VAR%GRID)
       CALL ACC_ATTACH_ASYNC(VAR%GRID,ACC_ASYNC_Q)
    ENDIF
  END SUBROUTINE COPYIN_WAVEDES1


!****************** SUBROUTINE COPYIN_WAVESPIN *************************
!
!***********************************************************************

  SUBROUTINE COPYIN_WAVESPIN(VAR)
    TYPE(wavespin) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!note the order of CELEN and CELTOT, since CELEN is a subarray of CELTOT,
!CELTOT needs to be created first, the same holds for FERWE and FERTOT
!$ACC ENTER DATA COPYIN(VAR%CPTWFP, &
!$ACC&                  VAR%CPROJ, &
!$ACC&                  VAR%FERTOT, &
!$ACC&                  VAR%FERWE, &
!$ACC&                  VAR%AUXTOT, &
!$ACC&                  VAR%CELTOT, &
!$ACC&                  VAR%CELEN) ASYNC(ACC_ASYNC_Q)
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%WDES)) THEN
       CALL ACC_COPYIN_TYPED_VAR(VAR%WDES)
       CALL ACC_ATTACH_ASYNC(VAR%WDES,ACC_ASYNC_Q)
    ENDIF
  END SUBROUTINE COPYIN_WAVESPIN


!****************** SUBROUTINE COPYIN_WAVEFUNA *************************
!
!***********************************************************************

  SUBROUTINE COPYIN_WAVEFUNA(VAR)
    TYPE(wavefuna) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!$ACC ENTER DATA COPYIN(VAR%CPTWFP, &
!$ACC&                  VAR%CPROJ, &
!$ACC&                  VAR%CW_RED, &
!$ACC&                  VAR%CPROJ_RED) ASYNC(ACC_ASYNC_Q)
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%WDES1)) THEN
       CALL ACC_COPYIN_TYPED_VAR(VAR%WDES1)
       CALL ACC_ATTACH_ASYNC(VAR%WDES1,ACC_ASYNC_Q)
    ENDIF
  END SUBROUTINE COPYIN_WAVEFUNA


!****************** SUBROUTINE COPYIN_WAVEFUN1 *************************
!
!***********************************************************************

  SUBROUTINE COPYIN_WAVEFUN1(VAR)
    TYPE(wavefun1) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers if associated
    IF (ASSOCIATED(VAR%CPTWFP)) THEN
!$ACC ENTER DATA COPYIN(VAR%CPTWFP) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%CPROJ)) THEN
!$ACC ENTER DATA COPYIN(VAR%CPROJ) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%CR)) THEN
!$ACC ENTER DATA COPYIN(VAR%CR) ASYNC(ACC_ASYNC_Q)
    ENDIF
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%WDES1)) THEN
       CALL ACC_COPYIN_TYPED_VAR(VAR%WDES1)
       CALL ACC_ATTACH_ASYNC(VAR%WDES1,ACC_ASYNC_Q)
    ENDIF
  END SUBROUTINE COPYIN_WAVEFUN1


!****************** SUBROUTINE COPYIN_POTCARS **************************
!
!***********************************************************************

  SUBROUTINE COPYIN_POTCARS(VAR)
    TYPE(potcar) :: VAR(:)
    INTEGER :: NT
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers if associated
    ! ATTENTION THIS IS NOT COMPLETE, BUT ALL THAT IS NEEDED WITHIN POTLOK
    DO NT = 1, SIZE(VAR)
!$ACC ENTER DATA COPYIN(VAR(NT)%PSP, &
!$ACC&                  VAR(NT)%PSPCOR, &
!$ACC&                  VAR(NT)%PSPRHO, &
!$ACC&                  VAR(NT)%LPS) ASYNC(ACC_ASYNC_Q)
    ENDDO
  END SUBROUTINE COPYIN_POTCARS


!****************** SUBROUTINE COPYIN_KPOINTS_FULL *********************
!
!***********************************************************************

  SUBROUTINE COPYIN_KPOINTS_FULL(VAR)
    TYPE(skpoints_full) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type variables
!$ACC ENTER DATA COPYIN(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers if associated
    IF (ASSOCIATED(VAR%NG_INDEX)) THEN
!$ACC ENTER DATA COPYIN(VAR%NG_INDEX) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%PHASE)) THEN
!$ACC ENTER DATA COPYIN(VAR%PHASE) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%LINV)) THEN
!$ACC ENTER DATA COPYIN(VAR%LINV) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%ROTMAP)) THEN
!$ACC ENTER DATA COPYIN(VAR%ROTMAP) ASYNC(ACC_ASYNC_Q)
    ENDIF
#ifndef oldsym
    IF (ASSOCIATED(VAR%RSSYMOP)) THEN
!$ACC ENTER DATA COPYIN(VAR%RSSYMOP) ASYNC(ACC_ASYNC_Q)
    ENDIF
#endif
  END SUBROUTINE COPYIN_KPOINTS_FULL


!****************** SUBROUTINE CREATE_WAVESPIN *************************
!
!***********************************************************************

  SUBROUTINE CREATE_WAVESPIN(VAR)
    TYPE(wavespin) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type variables
!$ACC ENTER DATA CREATE(VAR) ASYNC(ACC_ASYNC_Q)
    ! native type pointers
!note the order of CELEN and CELTOT, since CELEN is a subarray of CELTOT,
!CELTOT needs to be created first, the same holds for FERWE and FERTOT
!$ACC ENTER DATA CREATE(VAR%CPTWFP, &
!$ACC&                  VAR%CPROJ, &
!$ACC&                  VAR%FERTOT, &
!$ACC&                  VAR%FERWE, &
!$ACC&                  VAR%AUXTOT, &
!$ACC&                  VAR%CELTOT, &
!$ACC&                  VAR%CELEN) ASYNC(ACC_ASYNC_Q)
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%WDES)) THEN
       CALL ACC_COPYIN_TYPED_VAR(VAR%WDES)
       CALL ACC_ATTACH_ASYNC(VAR%WDES,ACC_ASYNC_Q)
    ENDIF
  END SUBROUTINE CREATE_WAVESPIN


!****************** SUBROUTINE DELETE_NONL_STRUCT **********************
!
!***********************************************************************

  SUBROUTINE DELETE_NONL_STRUCT(VAR)
    TYPE(nonl_struct) :: VAR
    ! early exit if possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%NITYP, &
!$ACC&                 VAR%ITYP, &
!$ACC&                 VAR%LMMAX, &
!$ACC&                 VAR%LMBASE, &
!$ACC&                 VAR%POSION, &
!$ACC&                 VAR%VKPT_SHIFT, &
!$ACC&                 VAR%CQFAK, &
!$ACC&                 VAR%QPROJ, &
!$ACC&                 VAR%CREXP) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_NONL_STRUCT


!****************** SUBROUTINE DELETE_NONLR_STRUCT *********************
!
!***********************************************************************

  SUBROUTINE DELETE_NONLR_STRUCT(VAR)
    TYPE(nonlr_struct) :: VAR
    ! logical variables
    INTEGER :: I
    ! early exit if possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! derived type pointers
    DO I=1,SIZE(VAR%BETA)
       CALL ACC_DELETE_TYPED_VAR(VAR%BETA(I))
    ENDDO
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%NITYP, &
!$ACC&                 VAR%ITYP, &
!$ACC&                 VAR%LMAX, &
!$ACC&                 VAR%LMMAX, &
!$ACC&                 VAR%CHANNELS, &
!$ACC&                 VAR%PSRMAX, &
!$ACC&                 VAR%RSMOOTH, &
!$ACC&                 VAR%POSION, &
!$ACC&                 VAR%VKPT_SHIFT, &
!$ACC&                 VAR%NLIMAX, &
!$ACC&                 VAR%NLI, &
!$ACC&                 VAR%LMBASE, &
!$ACC&                 VAR%NLIBASE, &
!$ACC&                 VAR%RPROJ, &
!$ACC&                 VAR%CRREXP) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_NONLR_STRUCT


!****************** SUBROUTINE DELETE_NONLR_PROJ ***********************
!
!***********************************************************************

  SUBROUTINE DELETE_NONLR_PROJ(VAR)
    TYPE(nonlr_proj) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%PSPRNL, &
!$ACC&                 VAR%LPS) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_NONLR_PROJ


!****************** SUBROUTINE DELETE_LAYOUT ***************************
!
!***********************************************************************

  SUBROUTINE DELETE_LAYOUT(VAR)
    TYPE(layout) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%I2, &
!$ACC&                 VAR%I3, &
!$ACC&                 VAR%INDEX) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_LAYOUT


!****************** SUBROUTINE DELETE_GRID_MAP *************************
!
!***********************************************************************

  SUBROUTINE DELETE_GRID_MAP(VAR)
    TYPE(grid_map) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%N, &
!$ACC&                 VAR%PTR, &
!$ACC&                 VAR%RMT, &
!$ACC&                 VAR%TBL, &
!$ACC&                 VAR%NI, &
!$ACC&                 VAR%PTRI, &
!$ACC&                 VAR%TBLI, &
!$ACC&                 VAR%RMTI) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_GRID_MAP


!****************** SUBROUTINE DELETE_GRID_3D **************************
!
!***********************************************************************

  SUBROUTINE DELETE_GRID_3D(VAR)
    TYPE(grid_3d) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! derived type pointers
    CALL ACC_DELETE_TYPED_VAR(VAR%RC)
!!    IF (.NOT.ACC_IS_PRESENT(VAR%RC,1)) CALL ACC_DETACH_ASYNC(VAR%RC,ACC_ASYNC_Q)

    CALL ACC_DELETE_TYPED_VAR(VAR%IN)
!!    IF (.NOT.ACC_IS_PRESENT(VAR%IN,1)) CALL ACC_DETACH_ASYNC(VAR%IN,ACC_ASYNC_Q)

    CALL ACC_DELETE_TYPED_VAR(VAR%RL)
!!    IF (.NOT.ACC_IS_PRESENT(VAR%RL,1)) CALL ACC_DETACH_ASYNC(VAR%RL,ACC_ASYNC_Q)

    CALL ACC_DELETE_TYPED_VAR(VAR%RC_IN)
!!    IF (.NOT.ACC_IS_PRESENT(VAR%RC_IN,1)) CALL ACC_DETACH_ASYNC(VAR%RC_IN,ACC_ASYNC_Q)

    CALL ACC_DELETE_TYPED_VAR(VAR%IN_RL)
!!    IF (.NOT.ACC_IS_PRESENT(VAR%IN_RL,1)) CALL ACC_DETACH_ASYNC(VAR%IN_RL,ACC_ASYNC_Q)

    ! native type pointer if associated
    IF (ASSOCIATED(VAR%RL_FFT)) THEN
!!$ACC EXIT DATA DELETE(VAR%RL_FFT) ASYNC(ACC_ASYNC_Q)
       CALL ACC_DELETE_TYPED_VAR(VAR%RL_FFT)
       IF (.NOT.ACC_IS_PRESENT(VAR%RL_FFT,1)) CALL ACC_DETACH_ASYNC(VAR%RL_FFT,ACC_ASYNC_Q)
    ENDIF
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%LPCTX, &
!$ACC&                 VAR%LPCTY, &
!$ACC&                 VAR%LPCTZ, &
!$ACC&                 VAR%LPCTX_, &
!$ACC&                 VAR%LPCTY_, &
!$ACC&                 VAR%LPCTZ_, &
!$ACC&                 VAR%FFTSCA, &
!$ACC&                 VAR%FFTWEIGHT, &
!$ACC&                 VAR%NINDPWCONJG, &
!$ACC&                 VAR%IND_IN_SPHERE, &
!$ACC&                 VAR%COMM, &
!$ACC&                 VAR%COMM_KIN, &
!$ACC&                 VAR%COMM_KINTER) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_GRID_3D


!****************** SUBROUTINE COPYIN_TRANSIT **************************
!
!***********************************************************************

  SUBROUTINE DELETE_TRANSIT(VAR)
    TYPE(transit) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%IND1, &
!$ACC&                  VAR%INDCOL) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_TRANSIT


!****************** SUBROUTINE DELETE_WAVEDES **************************
!
!***********************************************************************

  SUBROUTINE DELETE_WAVEDES(VAR)
    TYPE(wavedes) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%COMM)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_INTER)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_INTER) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_INB)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_INB) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_KINTER)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_KINTER) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_KIN)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_KIN) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_SHMEM)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_SHMEM) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_intra_node)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_intra_node) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_inter_node)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_inter_node) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%GRID)) THEN
       CALL ACC_DELETE_TYPED_VAR(VAR%GRID)
       CALL ACC_DETACH_ASYNC(VAR%GRID,ACC_ASYNC_Q)
    ENDIF
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%NPLWKP, &
!$ACC&                 VAR%NGVECTOR, &
!$ACC&                 VAR%NGVECTOR_POS, &
!$ACC&                 VAR%NPLWKP_TOT, &
!$ACC&                 VAR%NB_TOTK, &
!$ACC&                 VAR%PL_INDEX, &
!$ACC&                 VAR%PL_COL, &
!$ACC&                 VAR%NPRO_POS, &
!$ACC&                 VAR%LMMAX, &
!$ACC&                 VAR%LMBASE, &
!$ACC&                 VAR%NITYP, &
!$ACC&                 VAR%ITYP, &
!$ACC&                 VAR%NT_GLOBAL, &
!$ACC&                 VAR%VKPT, &
!$ACC&                 VAR%WTKPT, &
!$ACC&                 VAR%NINDPW, &
!$ACC&                 VAR%LUSEINV, &
!$ACC&                 VAR%NINDPW_INV, &
!$ACC&                 VAR%FFTSCA, &
!$ACC&                 VAR%MAP_TO_FULL, &
!$ACC&                 VAR%IGX, &
!$ACC&                 VAR%IGY, &
!$ACC&                 VAR%IGZ, &
!$ACC&                 VAR%DATAKE, &
!$ACC&                 VAR%AT_GAMMA) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_WAVEDES


!****************** SUBROUTINE DELETE_WAVEDES1 *************************
!
!***********************************************************************

  SUBROUTINE DELETE_WAVEDES1(VAR)
    TYPE(wavedes1) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%COMM)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_INTER)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_INTER) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_INB)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_INB) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_KINTER)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_KINTER) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_KIN)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_KIN) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_SHMEM)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_SHMEM) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_intra_node)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_intra_node) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%COMM_inter_node)) THEN
!$ACC EXIT DATA DELETE(VAR%COMM_inter_node) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%GRID)) THEN
       CALL ACC_DELETE_TYPED_VAR(VAR%GRID)
       CALL ACC_DETACH_ASYNC(VAR%GRID,ACC_ASYNC_Q)
    ENDIF
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%PL_INDEX, &
!$ACC&                 VAR%PL_COL, &
!$ACC&                 VAR%NPRO_POS, &
!$ACC&                 VAR%LMMAX, &
!$ACC&                 VAR%LMBASE, &
!$ACC&                 VAR%NITYP, &
!$ACC&                 VAR%ITYP, &
!$ACC&                 VAR%NT_GLOBAL, &
!$ACC&                 VAR%NINDPW, &
!$ACC&                 VAR%NINDPW_INV, &
!$ACC&                 VAR%FFTSCA, &
!$ACC&                 VAR%MAP_TO_FULL, &
!$ACC&                 VAR%IGX, &
!$ACC&                 VAR%IGY, &
!$ACC&                 VAR%IGZ, &
!$ACC&                 VAR%VKPT, &
!$ACC&                 VAR%DATAKE) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_WAVEDES1


!****************** SUBROUTINE DELETE_WAVESPIN *************************
!
!***********************************************************************

  SUBROUTINE DELETE_WAVESPIN(VAR)
    TYPE(wavespin) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%WDES)) THEN
       CALL ACC_DELETE_TYPED_VAR(VAR%WDES)
       CALL ACC_DETACH_ASYNC(VAR%WDES,ACC_ASYNC_Q)
    ENDIF
    ! native type pointers
!note the order of CELEN and CELTOT, since CELEN is a subarray of CELTOT,
!CELEN needs to be freed first, the same holds for FERWE and FERTOT.
!$ACC EXIT DATA DELETE(VAR%CPTWFP, &
!$ACC&                 VAR%CPROJ, &
!$ACC&                 VAR%FERWE, &
!$ACC&                 VAR%FERTOT, &
!$ACC&                 VAR%AUXTOT, &
!$ACC&                 VAR%CELEN, &
!$ACC&                 VAR%CELTOT) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_WAVESPIN


!****************** SUBROUTINE DELETE_WAVEFUNA *************************
!
!***********************************************************************

  SUBROUTINE DELETE_WAVEFUNA(VAR)
    TYPE(wavefuna) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%WDES1)) THEN
       CALL ACC_DELETE_TYPED_VAR(VAR%WDES1)
       CALL ACC_DETACH_ASYNC(VAR%WDES1,ACC_ASYNC_Q)
    ENDIF
    ! native type pointers
!$ACC EXIT DATA DELETE(VAR%CPTWFP, &
!$ACC&                  VAR%CPROJ, &
!$ACC&                  VAR%CW_RED, &
!$ACC&                  VAR%CPROJ_RED) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_WAVEFUNA


!****************** SUBROUTINE DELETE_WAVEFUN1 *************************
!
!***********************************************************************

  SUBROUTINE DELETE_WAVEFUN1(VAR)
    TYPE(wavefun1) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type pointers if associated
    IF (ASSOCIATED(VAR%CPTWFP)) THEN
!$ACC EXIT DATA DELETE(VAR%CPTWFP) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%CPROJ)) THEN
!$ACC EXIT DATA DELETE(VAR%CPROJ) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%CR)) THEN
!$ACC EXIT DATA DELETE(VAR%CR) ASYNC(ACC_ASYNC_Q)
    ENDIF
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%WDES1)) THEN
       CALL ACC_DELETE_TYPED_VAR(VAR%WDES1)
       CALL ACC_DETACH_ASYNC(VAR%WDES1,ACC_ASYNC_Q)
    ENDIF
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_WAVEFUN1


!****************** SUBROUTINE DELETE_POTCARS **************************
!
!***********************************************************************

  SUBROUTINE DELETE_POTCARS(VAR)
    TYPE(potcar) :: VAR(:)
    INTEGER :: NT
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type pointers
    ! ATTENTION THIS IS NOT COMPLETE, BUT ALL THAT IS NEEDED WITHIN POTLOK
    DO NT = 1, SIZE(VAR)
!$ACC EXIT DATA DELETE(VAR(NT)%PSP, &
!$ACC&                 VAR(NT)%PSPCOR, &
!$ACC&                 VAR(NT)%PSPRHO, &
!$ACC&                 VAR(NT)%LPS) ASYNC(ACC_ASYNC_Q)
    ENDDO
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_POTCARS


!****************** SUBROUTINE DELETE_KPOINTS_FULL *********************
!
!***********************************************************************

  SUBROUTINE DELETE_KPOINTS_FULL(VAR)
    TYPE(skpoints_full) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! native type pointers if associated
    IF (ASSOCIATED(VAR%NG_INDEX).AND.ACC_IS_PRESENT(VAR%NG_INDEX,1)) THEN
!$ACC EXIT DATA DELETE(VAR%NG_INDEX) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%PHASE).AND.ACC_IS_PRESENT(VAR%PHASE,1)) THEN
!$ACC EXIT DATA DELETE(VAR%PHASE) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%LINV)) THEN
!$ACC EXIT DATA DELETE(VAR%LINV) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IF (ASSOCIATED(VAR%ROTMAP)) THEN
!$ACC EXIT DATA DELETE(VAR%ROTMAP) ASYNC(ACC_ASYNC_Q)
    ENDIF
#ifndef oldsym
    IF (ASSOCIATED(VAR%RSSYMOP)) THEN
!$ACC EXIT DATA DELETE(VAR%RSSYMOP) ASYNC(ACC_ASYNC_Q)
    ENDIF
#endif
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE DELETE_KPOINTS_FULL


!****************** SUBROUTINE COPYOUT_WAVESPIN ************************
!
!***********************************************************************

  SUBROUTINE COPYOUT_WAVESPIN(VAR)
    TYPE(wavespin) :: VAR
    ! early exit of possible
    IF (.NOT.ACC_EXEC_ON) RETURN
    ! derived type pointers if associated
    IF (ASSOCIATED(VAR%WDES)) THEN
       CALL ACC_DELETE_TYPED_VAR(VAR%WDES)
       CALL ACC_DETACH_ASYNC(VAR%WDES,ACC_ASYNC_Q)
    ENDIF
    ! native type pointers
!note the order of CELEN and CELTOT, since CELEN is a subarray of CELTOT,
!CELEN needs to be freed first, the same holds for FERWE and FERTOT.
!$ACC EXIT DATA DELETE( VAR%FERWE, &
!$ACC&                  VAR%CELEN) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA COPYOUT(VAR%CPTWFP, &
!$ACC&                  VAR%CPROJ, &
!$ACC&                  VAR%FERTOT, &
!$ACC&                  VAR%AUXTOT, &
!$ACC&                  VAR%CELTOT) ASYNC(ACC_ASYNC_Q)
    ! native type variables
!$ACC EXIT DATA DELETE(VAR) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE COPYOUT_WAVESPIN


!****************** SUBROUTINE ACC_SET_ASYNC_Q *************************
!
!***********************************************************************

  SUBROUTINE ACC_SET_ASYNC_Q(VAR)
    INTEGER :: VAR
#ifdef ACC_DEBUG
    IF (ACC_DISABLE_SET_ASYNC_Q) RETURN
#endif
    IF (ACC_EXEC_ON) ACC_ASYNC_Q=MOD(VAR-1,ACC_MAX_STREAMS)+1
  END SUBROUTINE ACC_SET_ASYNC_Q


!****************** SUBROUTINE ACC_SYNC_ASYNC_Q ************************
!
!***********************************************************************

  SUBROUTINE SYNC_ASYNC_Q_ALL
!!#define use_cudaDeviceSynchronize
#ifdef use_cudaDeviceSynchronize
    USE cudafor, ONLY : cudaDeviceSynchronize,cudaGetErrorString,cudaGetLastError,cudaSuccess
    INTEGER :: ierr
    IF (ACC_EXEC_ON) THEN
       ierr = cudaDeviceSynchronize()
       IF (ierr/=cudaSuccess) WRITE(*,*) 'SYNC_ASYNC_Q_ALL: internal error:',cudaGetErrorString(cudaGetLastError())
    ENDIF
    RETURN
#endif
!$ACC WAIT IF(ACC_EXEC_ON)
  END SUBROUTINE SYNC_ASYNC_Q_ALL

  SUBROUTINE SYNC_ASYNC_Q_RANGE(NQ)
!!#define use_cudaStreamSynchronize
#ifdef use_cudaStreamSynchronize
    USE openacc, ONLY : acc_get_cuda_stream
    USE cudafor, ONLY : cudaStreamSynchronize,cudaGetErrorString,cudaGetLastError,cudaSuccess
    INTEGER :: ierr
#endif
    INTEGER :: NQ
    ! local variables
    INTEGER :: I
    IF (.NOT.ACC_EXEC_ON) RETURN
    DO I=1,NQ
#ifdef use_cudaStreamSynchronize
       ierr = cudaStreamSynchronize(acc_get_cuda_stream(I))
       IF (ierr/=cudaSuccess) WRITE(*,*) 'SYNC_ASYNC_Q_RANGE: internal error:',cudaGetErrorString(cudaGetLastError())
#else
!$ACC WAIT(I)
#endif
    ENDDO
  END SUBROUTINE SYNC_ASYNC_Q_RANGE

  SUBROUTINE SYNC_ASYNC_Q_MASKED(NQ,MASK)
!!#define use_cudaStreamSynchronize
#ifdef use_cudaStreamSynchronize
    USE openacc, ONLY : acc_get_cuda_stream
    USE cudafor, ONLY : cudaStreamSynchronize,cudaGetErrorString,cudaGetLastError,cudaSuccess
    INTEGER :: ierr
#endif
    INTEGER :: NQ
    LOGICAL :: MASK(:)
    ! local variables
    INTEGER :: I
    IF (.NOT.ACC_EXEC_ON) RETURN
    DO I=1,NQ
       IF (MASK(I)) THEN
#ifdef use_cudaStreamSynchronize
          ierr = cudaStreamSynchronize(acc_get_cuda_stream(I))
          IF (ierr/=cudaSuccess) WRITE(*,*) 'SYNC_ASYNC_Q_MASKED: internal error:',cudaGetErrorString(cudaGetLastError())
#else
!$ACC WAIT(I)
#endif
       ENDIF
    ENDDO
  END SUBROUTINE SYNC_ASYNC_Q_MASKED


!****************** SUBROUTINE ACC_PUSH_EXEC_ON ************************
!
!***********************************************************************

  SUBROUTINE ACC_PUSH_EXEC_ON(VAR)
    LOGICAL :: VAR
!$ACC WAIT
    IF (ACC_EXEC_ON_LEVEL==MAXLEVEL) THEN
       CALL vtutor%error("ACC_PUSH_EXEC_ON: ERROR: stack is full")
    ENDIF
    ACC_EXEC_ON_LEVEL=ACC_EXEC_ON_LEVEL+1
    ACC_EXEC_ON_STACK(ACC_EXEC_ON_LEVEL)=ACC_EXEC_ON
    ACC_EXEC_ON=VAR.AND.ACC_EXEC_ALLOW
  END SUBROUTINE ACC_PUSH_EXEC_ON


!****************** SUBROUTINE ACC_POP_EXEC_ON *************************
!
!***********************************************************************

  SUBROUTINE ACC_POP_EXEC_ON
!$ACC WAIT
    IF (ACC_EXEC_ON_LEVEL==0) THEN
       CALL vtutor%error("ACC_POP_EXEC_ON: ERROR: stack is empty")
    ENDIF
    ACC_EXEC_ON=ACC_EXEC_ON_STACK(ACC_EXEC_ON_LEVEL)
    ACC_EXEC_ON_LEVEL=ACC_EXEC_ON_LEVEL-1
  END SUBROUTINE ACC_POP_EXEC_ON
#else
  CONTAINS
  SUBROUTINE acc_dummy
    WRITE(*,*)'Im a DEC compiler so I need this line'
  END SUBROUTINE acc_dummy
#endif
  END MODULE mopenacc

#ifdef _OPENACC
!***********************************************************************
!
! This module contains the infrastructure to create and destroy
! cuFFT plans for grid_3d structures.
!
!***********************************************************************

  MODULE mcufft
  USE prec
  USE mgrid_struct_def, ONLY : grid_3d

  USE cufft
  USE openacc, ONLY : ACC_ASYNC_SYNC,ACC_ASYNC_NOVAL,acc_get_cuda_stream
  USE mopenacc, ONLY : ACC_EXEC_ALLOW,ACC_ASYNC_ASYNC

  USE string, ONLY: str
  USE tutor, ONLY: vtutor

  USE iso_c_binding

  IMPLICIT NONE

!!  PRIVATE

  PUBLIC :: ACC_CUFFT_MAKEPLANS,ACC_CUFFT_MAKEPLAN,ACC_CUFFT_GETPLAN,ACC_CUFFT_DESTROYPLAN
  PUBLIC :: acc_cufft_plan

  TYPE acc_cufft_plan
    TYPE(grid_3d), POINTER :: GRID
    INTEGER :: stream
    INTEGER :: nbatch
    INTEGER :: CUFFT_3D_Z2Z,CUFFT_3D_D2Z,CUFFT_3D_Z2D
    INTEGER :: CUFFT_RC_Z2Z,CUFFT_IN_Z2Z
    INTEGER :: CUFFT_RL_Z2Z,CUFFT_RL_D2Z,CUFFT_RL_Z2D
    CHARACTER(LEN=:), ALLOCATABLE :: tag
#if PGI_BETWEEN_20_9_AND_XX_X
    LOGICAL :: LFORCE_ACC_WAIT
#endif
  END TYPE acc_cufft_plan

  TYPE(acc_cufft_plan), ALLOCATABLE :: plans(:)

  INTEGER :: createdPlans = 0

#if PGI_BEFORE_XX_X
  INTEGER :: CUFFT_VER = -1
#endif
#if PGI_BETWEEN_20_9_AND_XX_X
  INTEGER :: DEVICE_CACHE_SIZE = -1
#endif

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE planAssign
  END INTERFACE

  INTERFACE OPERATOR(==)
    MODULE PROCEDURE planCompare
  END INTERFACE

  INTEGER(c_int), BIND(c,name="CUFFT_WORKAREA_MINIMAL") :: CUFFT_WORKAREA_MINIMAL

  INTERFACE
    INTEGER(c_int) FUNCTION cufftXtSetWorkAreaPolicy(plan, policy, workSize) BIND(c, name="cufftXtSetWorkAreaPolicy")
        IMPORT
        INTEGER(c_int), VALUE :: plan
        INTEGER(c_int), VALUE :: policy
        INTEGER(c_size_t)     :: workSize
    END FUNCTION cufftXtSetWorkAreaPolicy
  END INTERFACE

  CONTAINS

!****************** SUBROUTINE ACC_CUFFT_MAKEPLANS *********************
!
!***********************************************************************

  SUBROUTINE ACC_CUFFT_MAKEPLANS(GRID,NSTREAMS)
    TYPE(grid_3d), INTENT(IN), TARGET :: GRID
    INTEGER, INTENT(IN), OPTIONAL :: NSTREAMS
    ! local variables
    TYPE(acc_cufft_plan) :: new
    INTEGER :: n,i,ierr

    PROFILING_START('acc_cufft_makeplans')

    IF (.NOT.ACC_EXEC_ALLOW) THEN
       PROFILING_STOP('acc_cufft_makeplans')
       RETURN
    ENDIF

    ! Always create a plan for ACC_ASYNC_SYNC
    CALL ACC_CUFFT_MAKEPLAN(GRID,INT(ACC_ASYNC_SYNC))

    ! And for ACC_ASYNC_ASYNC
    CALL ACC_CUFFT_MAKEPLAN(GRID,ACC_ASYNC_ASYNC)

    ! Then for additional streams if requested
    N=1; IF (PRESENT(NSTREAMS)) N=NSTREAMS
    DO I=1,N
      CALL ACC_CUFFT_MAKEPLAN(GRID,I)
    ENDDO

    ! And last but not least for ACC_ASYNC_NOVAL
    CALL ACC_CUFFT_MAKEPLAN(GRID,INT(ACC_ASYNC_NOVAL))

    PROFILING_STOP('acc_cufft_makeplans')

  END SUBROUTINE ACC_CUFFT_MAKEPLANS


!****************** SUBROUTINE ACC_CUFFT_MAKEPLAN **********************
!
!***********************************************************************

  SUBROUTINE ACC_CUFFT_MAKEPLAN(GRID,ISTREAM,NBATCH,NSTRIDE,TAG)
    TYPE(grid_3d), INTENT(IN), TARGET      :: GRID
    INTEGER, INTENT(IN)                    :: ISTREAM
    INTEGER, INTENT(IN), OPTIONAL          :: NBATCH
    INTEGER, INTENT(IN), OPTIONAL          :: NSTRIDE
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: TAG
    ! local variables
    TYPE(acc_cufft_plan) :: NEW
    INTEGER :: NG(3),BATCHES(3),LD,IERR
#define use_cufft_workarea_minimal
#ifdef use_cufft_workarea_minimal
    INTEGER(KIND=c_size_t) :: ISIZE
    INTEGER :: IDUMMY(1)
#endif

    PROFILING_START('acc_cufft_makeplan')

    IF (.NOT.ACC_EXEC_ALLOW) THEN
       PROFILING_STOP('acc_cufft_makeplan')
       RETURN
    ENDIF

    IERR=0

    NEW%GRID   =>GRID
    NEW%STREAM = ISTREAM

    NEW%NBATCH = 1 ; IF (PRESENT(NBATCH)) NEW%NBATCH = NBATCH

    ! when we already have a plan for this GRID+ISTREAM+NBATCH we may exit here
    IF (planFind(NEW)/=0) THEN
       PROFILING_STOP('acc_cufft_makeplan')
       RETURN
    ENDIF

#if PGI_BEFORE_XX_X
    IF (CUFFT_VER == -1 ) IERR = cufftGetVersion(CUFFT_VER)
#endif

    NEW%TAG = "" ; IF (PRESENT(TAG)) NEW%TAG = TAG

    NG=NEW%GRID%NGPTAR

    IF (NEW%GRID%RL%NFAST==1 .AND. NEW%GRID%RL_FFT%NFAST==1) THEN
       ! serial fft
       IF (NEW%NBATCH==1) THEN
          ! create single fft plans
          IF (NEW%GRID%REAL2CPLX) THEN
             ! 3d real-to-complex
             IERR = IERR + cufftPlan3d(NEW%CUFFT_3D_D2Z,NG(3),NG(2),NG(1),CUFFT_D2Z)
             IERR = IERR + cufftSetStream(NEW%CUFFT_3D_D2Z,acc_get_cuda_stream(NEW%STREAM))
             ! 3d complex-to-real
             IERR = IERR + cufftPlan3d(NEW%CUFFT_3D_Z2D,NG(3),NG(2),NG(1),CUFFT_Z2D)
             IERR = IERR + cufftSetStream(NEW%CUFFT_3D_Z2D,acc_get_cuda_stream(NEW%STREAM))
          ELSE
             ! 3d complex-to-complex
             IERR = IERR + cufftPlan3d(NEW%CUFFT_3D_Z2Z,NG(3),NG(2),NG(1),CUFFT_Z2Z)
             IERR = IERR + cufftSetStream(NEW%CUFFT_3D_Z2Z,acc_get_cuda_stream(NEW%STREAM))
          ENDIF
       ELSE
          ! size required for in-place FFTs (in complex words)
          LD=NEW%GRID%MPLWV
          ! allow for custom strided array access
          IF (PRESENT(NSTRIDE)) LD=NSTRIDE

          ! create batched fft plans
          IF (NEW%GRID%REAL2CPLX) THEN
             ! batched 3d real-to-complex
             IERR = IERR + cufftPlanMany(NEW%CUFFT_3D_D2Z,3,[NG(3),NG(2),NG(1)], &
                                         [NG(3),NG(2), NG(1)+2   ],1,2*LD, &
                                         [NG(3),NG(2),(NG(1)+2)/2],1,  LD, &
                                         CUFFT_D2Z,NEW%NBATCH)
             IERR = IERR + cufftSetStream(NEW%CUFFT_3D_D2Z,acc_get_cuda_stream(NEW%STREAM))
             ! batched 3d real-to-complex
             IERR = IERR + cufftPlanMany(NEW%CUFFT_3D_Z2D,3,[NG(3),NG(2),NG(1)], &
                                         [NG(3),NG(2),(NG(1)+2)/2],1,  LD, &
                                         [NG(3),NG(2), NG(1)+2   ],1,2*LD, &
                                         CUFFT_Z2D,NEW%NBATCH)
             IERR = IERR + cufftSetStream(NEW%CUFFT_3D_Z2D,acc_get_cuda_stream(NEW%STREAM))
          ELSE
             ! batched 3d complex-to-complex
#ifdef use_cufft_workarea_minimal
             IERR = IERR + cufftCreate(NEW%CUFFT_3D_Z2Z)
             ! write(*,*) 'ierror cufftCreate :',ierr

             IERR = IERR + cufftSetAutoAllocation(NEW%CUFFT_3D_Z2Z,0)
             ! write(*,*) 'ierror SetAutoAlloc:',ierr

             IERR = IERR + cufftMakePlanMany(NEW%CUFFT_3D_Z2Z,3,[NG(3),NG(2),NG(1)], &
                                             [NG(3),NG(2),NG(1)],1,LD, &
                                             [NG(3),NG(2),NG(1)],1,LD, &
                                             CUFFT_Z2Z,NEW%NBATCH,ISIZE)
             ! write(*,*) 'ierror MakePlanMany:',ierr,isize

             IERR = IERR + cufftXtSetWorkAreaPolicy(NEW%CUFFT_3D_Z2Z,CUFFT_WORKAREA_MINIMAL,ISIZE)
             ! write(*,*) 'ierror SetWorkAreaPolicy:',ierr,cufft_workarea_minimal
#if PGI_BEFORE_XX_X
             IF (CUFFT_VER < 10500) THEN
!$ACC ENTER DATA CREATE(IDUMMY)
!$ACC HOST_DATA USE_DEVICE(IDUMMY)
             IERR = IERR + cufftSetWorkArea(NEW%CUFFT_3D_Z2Z,IDUMMY)
!$ACC END HOST_DATA
!$ACC EXIT DATA DELETE(IDUMMY)
             ! write(*,*) 'ierror SetWorkArea:',ierr
             ENDIF
#endif
             IERR = IERR + cufftSetStream(NEW%CUFFT_3D_Z2Z,acc_get_cuda_stream(NEW%STREAM))
             ! write(*,*) 'ierror SetStream:',ierr
#else
             IERR = IERR + cufftPlanMany(NEW%CUFFT_3D_Z2Z,3,[NG(3),NG(2),NG(1)], &
                                         [NG(3),NG(2),NG(1)],1,LD,[NG(3),NG(2),NG(1)],1,LD,CUFFT_Z2Z,NEW%NBATCH)
             IERR = IERR + cufftSetStream(NEW%CUFFT_3D_Z2Z,acc_get_cuda_stream(NEW%STREAM))
#endif
          ENDIF
       ENDIF
    ELSE
       IF (NEW%NBATCH/=1) THEN
          CALL vtutor%bug("ACC_CUFFT_MAKEPLAN: batched parallel ffts not supported yet",__FILE__,__LINE__)
       ENDIF

       ! parallel fft
       BATCHES(1)=NEW%GRID%RC%NCOL
       BATCHES(2)=NEW%GRID%IN%NCOL
       BATCHES(3)=NEW%GRID%RL_FFT%NCOL

       IF (batches(1)>0) THEN
          ! fft on data with reciprocal space layout
          IERR = IERR + cufftPlanMany(NEW%CUFFT_RC_Z2Z,1,NG(1), &
                                      NG(1),1,NG(1),NG(1),1,NG(1),CUFFT_Z2Z,BATCHES(1))
          IERR = IERR + cufftSetStream(NEW%CUFFT_RC_Z2Z,acc_get_cuda_stream(NEW%STREAM))
       ENDIF
       IF (BATCHES(2)>0) THEN
          ! fft on data with intermediate layout
          IERR = IERR + cufftPlanMany(NEW%CUFFT_IN_Z2Z,1,NG(2), &
                                      NG(2),BATCHES(2),1,NG(2),BATCHES(2),1,CUFFT_Z2Z,BATCHES(2))
          IERR = IERR + cufftSetStream(NEW%CUFFT_IN_Z2Z,acc_get_cuda_stream(NEW%STREAM))
       ENDIF
       IF (BATCHES(3)>0) THEN
          ! fft on data with real space layout
          IF (NEW%GRID%REAL2CPLX) THEN
             ! real-to-complex
             IERR = IERR + cufftPlanMany(NEW%CUFFT_RL_D2Z,1,NG(3), &
                                         NG(3),1,NG(3)+2,NG(3),1,(NG(3)+2)/2,CUFFT_D2Z,BATCHES(3))
             IERR = IERR + cufftSetStream(NEW%CUFFT_RL_D2Z,acc_get_cuda_stream(NEW%STREAM))
             ! complex-to-real
             IERR = IERR + cufftPlanMany(NEW%CUFFT_RL_Z2D,1,NG(3), &
                                         NG(3),1,(NG(3)+2)/2,NG(3),1,NG(3)+2,CUFFT_Z2D,BATCHES(3))
             IERR = IERR + cufftSetStream(NEW%CUFFT_RL_Z2D,acc_get_cuda_stream(NEW%STREAM))
          ELSE
             ! complex-to-complex
             IERR = IERR + cufftPlanMany(NEW%CUFFT_RL_Z2Z,1,NG(3), &
                                         NG(3),1,NG(3),NG(3),1,NG(3),CUFFT_Z2Z,BATCHES(3))
             IERR = IERR + cufftSetStream(NEW%CUFFT_RL_Z2Z,acc_get_cuda_stream(NEW%STREAM))
          ENDIF
       ENDIF
    ENDIF

#if PGI_BETWEEN_20_9_AND_XX_X
    CALL L2CACHE_OPT_WAR(NEW)
#endif

    ! In case something went wrong along the way ...
    IF (IERR/=0) THEN
       CALL vtutor%error("ACC_CUFFT_MAKEPLAN: could not create plan")
    ENDIF

    CALL planAdd(NEW)

    PROFILING_STOP('acc_cufft_makeplan')

  END SUBROUTINE ACC_CUFFT_MAKEPLAN


!****************** SUBROUTINE ACC_CUFFT_DESTROYPLAN *******************
!
!***********************************************************************

  SUBROUTINE ACC_CUFFT_DESTROYPLAN(GRID,ISTREAM,NBATCH,TAG)
    TYPE(grid_3d), INTENT(IN), POINTER     :: GRID
    INTEGER, INTENT(IN), OPTIONAL          :: ISTREAM
    INTEGER, INTENT(IN), OPTIONAL          :: NBATCH
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: TAG
    ! local variables
    TYPE(acc_cufft_plan), ALLOCATABLE :: tempPlans(:)
    INTEGER :: i,j,ierr
    LOGICAL :: keepPlan

    IF (.NOT.ACC_EXEC_ALLOW) RETURN

    j=0; ierr=0

    ALLOCATE(tempPlans(createdPlans))

    DO i=1,createdPlans
       keepPlan = .FALSE.

       IF (.NOT.ASSOCIATED(GRID,plans(i)%GRID)) keepPLan = .TRUE.

       IF (PRESENT(ISTREAM)) THEN
          IF (plans(i)%stream/=ISTREAM) keepPlan = .TRUE.
       ENDIF

       IF (PRESENT(NBATCH)) THEN
          IF (plans(i)%NBATCH/=NBATCH) keepPlan = .TRUE.
       ENDIF

       IF (PRESENT(TAG)) keepPlan = (plans(i)%TAG/=TAG)

       IF (keepPlan) THEN
          j=j+1; tempPlans(j)=plans(i)
       ELSE
          CALL DESTROYPLAN
       ENDIF
    ENDDO

    IF (ierr/=0) THEN
       CALL vtutor%error("ACC_CUFFT_DESTROYPLAN: could not destroy plans")
    ENDIF

    DEALLOCATE(plans)

    IF (j/=0) THEN
       ALLOCATE(plans(j))
       DO i=1,j
          plans(i)=tempPlans(i)
       ENDDO
    ENDIF

    DEALLOCATE(tempPlans)

    createdPlans=j

    CONTAINS

    SUBROUTINE DESTROYPLAN
       IF (GRID%RL%NFAST==1) THEN
          ! serial fft
          IF (GRID%REAL2CPLX) THEN
             ierr = ierr + cufftDestroy(plans(i)%CUFFT_3D_D2Z)
             ierr = ierr + cufftDestroy(plans(i)%CUFFT_3D_Z2D)
          ELSE
             ierr = ierr + cufftDestroy(plans(i)%CUFFT_3D_Z2Z)
          ENDIF
       ELSE
          IF (GRID%RC%NCOL>0) &
             ierr = ierr + cufftDestroy(plans(i)%CUFFT_RC_Z2Z)
          IF (GRID%IN%NCOL>0) &
             ierr = ierr + cufftDestroy(plans(i)%CUFFT_IN_Z2Z)
          IF (GRID%RL_FFT%NCOL>0) THEN
             IF (GRID%REAL2CPLX) THEN
                ierr = ierr + cufftDestroy(plans(i)%CUFFT_RL_D2Z)
                ierr = ierr + cufftDestroy(plans(i)%CUFFT_RL_Z2D)
             ELSE
                ierr = ierr + cufftDestroy(plans(i)%CUFFT_RL_Z2Z)
             ENDIF
          ENDIF
       ENDIF
    END SUBROUTINE DESTROYPLAN

  END SUBROUTINE ACC_CUFFT_DESTROYPLAN


!****************** FUNCTION ACC_CUFFT_GETPLAN *************************
!
!***********************************************************************

  FUNCTION ACC_CUFFT_GETPLAN(GRID,ISTREAM,NBATCH,RES) RESULT(ISTAT)
    TYPE(grid_3d), INTENT(IN), TARGET :: GRID
    INTEGER, INTENT(IN) :: ISTREAM
    INTEGER, INTENT(IN) :: NBATCH
    TYPE(acc_cufft_plan), INTENT(OUT) :: RES
    INTEGER, INTENT(OUT) :: ISTAT
    ! local variables
    TYPE(acc_cufft_plan) :: q
    INTEGER :: i

    q%GRID   =>GRID
    q%stream = ISTREAM
    q%nbatch = NBATCH

    ISTAT=0

    ! search matching plan
    i=planFind(q)

    IF (i/=0) THEN
       ! return matching plan
       RES=plans(i)
    ELSE
      ! no match found: return in error
      ISTAT=-1
    ENDIF
  END FUNCTION ACC_CUFFT_GETPLAN


!****************** SUBROUTINE PLANADD *********************************
!
!***********************************************************************

  SUBROUTINE planAdd(a)
    TYPE(acc_cufft_plan), INTENT(IN) :: a
    ! local variables
    INTEGER :: i

    ! Early exit if plan already exists
    IF (planFind(a)/=0) RETURN

    ! Find storage space for a new plan
    i=planAlloc()

    ! Add plan
    plans(i)=a
  END SUBROUTINE planAdd


!****************** FUNCTION PLANALLOC *********************************
!
!***********************************************************************

  FUNCTION planAlloc
    INTEGER :: planAlloc
    ! local variables
    TYPE(acc_cufft_plan), ALLOCATABLE :: tempPlans(:)
    INTEGER :: i,oldSize

    ! on first call the plans array is not allocated
    IF (.NOT.ALLOCATED(plans)) THEN
       ALLOCATE(plans(1))
       NULLIFY(plans(1)%grid)
       createdPlans=0

    ! when the plans array is full
    ELSEIF(createdPlans==SIZE(plans)) THEN
       oldSize=SIZE(plans)
       ALLOCATE(tempPlans(oldsize))
       DO i=1,oldSize
          tempPlans(i)=plans(i)
       ENDDO
       DEALLOCATE(plans)
       ALLOCATE(plans(oldsize+1))
       DO i=1,oldSize
          plans(i)=tempPlans(i)
       ENDDO
       DEALLOCATE(tempPlans)
       NULLIFY(plans(oldsize+1)%grid)

    ! error .. this may not happen
    ELSEIF(createdPlans>SIZE(plans)) THEN
       STOP
    ENDIF

    ! advance the number of created plans by one
    createdPlans=createdPlans+1

    ! and pass the index back
    planAlloc=createdPlans
  END FUNCTION planAlloc


!****************** FUNCTION PLANFIND **********************************
!
!***********************************************************************

  FUNCTION planFind(a)
    TYPE(acc_cufft_plan), INTENT(IN) :: a
    INTEGER :: planFind
    ! local variables
    INTEGER :: i

    planFind=0
    DO i=1,createdPlans
       IF (plans(i)==a) THEN
          planFind=i; RETURN
       ENDIF
    ENDDO
  END FUNCTION planFind


!****************** SUBROUTINE PLANASSIGN ******************************
!
!***********************************************************************

  SUBROUTINE planAssign(a,b)
    TYPE(acc_cufft_plan), INTENT(OUT) :: a
    TYPE(acc_cufft_plan), INTENT(IN)  :: b

    a%GRID            =>b%GRID
    a%stream          = b%stream
    a%nbatch          = b%nbatch
    a%tag             = b%tag

    a%CUFFT_3D_Z2Z    = b%CUFFT_3D_Z2Z
    a%CUFFT_3D_D2Z    = b%CUFFT_3D_D2Z
    a%CUFFT_3D_Z2D    = b%CUFFT_3D_Z2D

    a%CUFFT_RC_Z2Z    = b%CUFFT_RC_Z2Z
    a%CUFFT_IN_Z2Z    = b%CUFFT_IN_Z2Z
    a%CUFFT_RL_Z2Z    = b%CUFFT_RL_Z2Z
    a%CUFFT_RL_D2Z    = b%CUFFT_RL_D2Z
    a%CUFFT_RL_Z2D    = b%CUFFT_RL_Z2D

#if PGI_BETWEEN_20_9_AND_XX_X
    a%LFORCE_ACC_WAIT = b%LFORCE_ACC_WAIT
#endif
  END SUBROUTINE planAssign


!****************** FUNCTION PLANCOMPARE *******************************
!
!***********************************************************************

  FUNCTION planCompare(a,b) RESULT(res)
    TYPE(acc_cufft_plan), INTENT(IN) :: a,b
    LOGICAL, INTENT(OUT) :: res

    res=.FALSE.

    IF (.NOT.ASSOCIATED(a%GRID,b%GRID)) RETURN
    IF ( a%stream/=b%stream) RETURN
    IF ( a%nbatch/=b%nbatch) RETURN

    res=.TRUE.
  END FUNCTION planCompare

#if PGI_BETWEEN_20_9_AND_XX_X
  SUBROUTINE L2CACHE_OPT_WAR(PLAN)
    USE mopenacc_struct_def
    TYPE(acc_cufft_plan) :: PLAN
    INTEGER :: IERR
    INTEGER :: dataSize

    PLAN%LFORCE_ACC_WAIT = .FALSE.

    IF ( PLAN%nbatch == 1 ) RETURN

    IF ( CUFFT_VER == -1 ) IERR = cufftGetVersion(CUFFT_VER)
    IF ( CUFFT_VER < 10300 .OR. CUFFT_VER > 10501 ) RETURN

    IF ( .NOT. IS_DIM_AFFECTED(PLAN%GRID%NGPTAR(1)) ) RETURN
    IF ( .NOT. IS_DIM_AFFECTED(PLAN%GRID%NGPTAR(2)) ) RETURN
    IF ( .NOT. IS_DIM_AFFECTED(PLAN%GRID%NGPTAR(3)) ) RETURN

    dataSize = PRODUCT(PLAN%GRID%NGPTAR) * PLAN%nbatch * 16
    IF ( DEVICE_CACHE_SIZE == -1 ) DEVICE_CACHE_SIZE=GET_DEVICE_CACHE_SIZE()
    IF ( DEVICE_CACHE_SIZE > dataSize ) RETURN

    PLAN%LFORCE_ACC_WAIT = .TRUE.

    CONTAINS
      LOGICAL FUNCTION IS_DIM_AFFECTED(X)
        INTEGER :: X
        INTEGER, PARAMETER, DIMENSION(30) :: LIST=(/2, 4, 8, 16, 32, 64, 128, &
                                                    3, 9, 81, 243, 729, 2187, 6561, &
                                                    5, 25, 125, 625, 3125, &
                                                    6, 36, 216, 1296, 7776, &
                                                    7, 49, 343, 2401, &
                                                   11, 121/)
        IS_DIM_AFFECTED = ANY(LIST==X)
        RETURN
      END FUNCTION

      INTEGER FUNCTION GET_DEVICE_CACHE_SIZE()
        USE cudafor, only : cudaGetDeviceProperties, cudaDeviceProp
        USE openacc, only : ACC_DEVICE_NVIDIA, acc_get_device_num
        INTEGER :: IERR
        type (cudaDeviceProp) :: prop

        IERR =cudaGetDeviceProperties(prop, acc_get_device_num(ACC_DEVICE_NVIDIA))
        GET_DEVICE_CACHE_SIZE = prop%l2CacheSize
        IF ( prop%major >= 8 ) GET_DEVICE_CACHE_SIZE=GET_DEVICE_CACHE_SIZE/2
        RETURN
      END FUNCTION GET_DEVICE_CACHE_SIZE

  END SUBROUTINE L2CACHE_OPT_WAR
#endif

  END MODULE mcufft


!***********************************************************************
!
! This module only holds a few variables related to the use of cublas
! and a subroutine to initialize the cublasHandle ACC_CUBLASHANDLE.
!
! The actual cublas calls are wrapped in ACC_* subroutines that live
! outside of mcublas (see immediately below the aforementioned module).
!
!***********************************************************************

  MODULE mcublas
  USE prec

  USE cublas
  USE openacc,  ONLY : acc_get_cuda_stream

  USE mopenacc, ONLY : ACC_EXEC_ALLOW,ACC_EXEC_ON,ACC_ASYNC_Q

  USE tutor, ONLY: vtutor

  IMPLICIT NONE

  TYPE(cublasHandle) :: ACC_CUBLASHANDLE
!$OMP THREADPRIVATE(ACC_CUBLASHANDLE)

  INTEGER(4) :: CUBLAS_VERSION=0

  CONTAINS

  INTEGER FUNCTION CHAR_TO_OP(MODE)
    CHARACTER :: MODE
    SELECT CASE (MODE)
      CASE ('N','n')
        CHAR_TO_OP = CUBLAS_OP_N
      CASE ('T','t')
        CHAR_TO_OP = CUBLAS_OP_T
      CASE ('C','c')
        CHAR_TO_OP = CUBLAS_OP_C
      CASE DEFAULT
        CALL vtutor%bug("CHAR_TO_OP: invalid argument:" // MODE,__FILE__,__LINE__)
    END SELECT
  END FUNCTION CHAR_TO_OP


  INTEGER FUNCTION CHAR_TO_FILL(MODE)
    CHARACTER :: MODE
    SELECT CASE (MODE)
      CASE ('U','u')
        CHAR_TO_FILL = CUBLAS_FILL_MODE_UPPER
      CASE ('L','l')
        CHAR_TO_FILL = CUBLAS_FILL_MODE_LOWER
      !CASE ('F','f')
      !  CHAR_TO_FILL = CUBLAS_FILL_MODE_FULL
      CASE DEFAULT
        CALL vtutor%bug("CHAR_TO_FILL: invalid argument:" // MODE,__FILE__,__LINE__)
    END SELECT
  END FUNCTION CHAR_TO_FILL


  INTEGER FUNCTION CHAR_TO_SIDE(MODE)
    CHARACTER :: MODE
    SELECT CASE (MODE)
      CASE ('L','l')
        CHAR_TO_SIDE = CUBLAS_SIDE_LEFT
      CASE ('R','r')
        CHAR_TO_SIDE = CUBLAS_SIDE_RIGHT
      CASE DEFAULT
        CALL vtutor%bug("CHAR_TO_SIDE: invalid argument:" // MODE,__FILE__,__LINE__)
    END SELECT
  END FUNCTION CHAR_TO_SIDE


  INTEGER FUNCTION CHAR_TO_DIAG(MODE)
    CHARACTER :: MODE
    SELECT CASE (MODE)
      CASE ('N','n')
        CHAR_TO_DIAG = CUBLAS_DIAG_NON_UNIT
      CASE ('U','u')
        CHAR_TO_DIAG = CUBLAS_DIAG_UNIT
      CASE DEFAULT
        CALL vtutor%bug("CHAR_TO_DIAG: invalid argument:" // MODE,__FILE__,__LINE__)
    END SELECT
  END FUNCTION CHAR_TO_DIAG


  SUBROUTINE ACC_INIT_CUBLAS()
    INTEGER(4) :: IERR

    PROFILING_START('acc_init_cublas')

    IF (.NOT.ACC_EXEC_ALLOW) THEN
       PROFILING_STOP('acc_init_cublas')
       RETURN
    ENDIF

!$OMP PARALLEL
    IERR = cublasCreate(ACC_CUBLASHANDLE)
    IF (IERR/=0) THEN
       CALL vtutor%error("ERROR setting up cublasHandle")
    ENDIF
!$OMP END PARALLEL

    CALL GET_CUBLAS_VERSION()

    PROFILING_STOP('acc_init_cublas')

  END SUBROUTINE ACC_INIT_CUBLAS


  SUBROUTINE ACC_FINISH_CUBLAS()
    INTEGER(4) :: IERR

    IF (.NOT.ACC_EXEC_ALLOW) RETURN

!$OMP PARALLEL
    IERR = cublasDestroy(ACC_CUBLASHANDLE)
    IF (IERR/=0) THEN
       CALL vtutor%error("ERROR releasing cublasHandle")
    ENDIF
!$OMP END PARALLEL
  END SUBROUTINE ACC_FINISH_CUBLAS

  SUBROUTINE GET_CUBLAS_VERSION()
    USE tutor, ONLY: vtutor
    INTEGER(4) :: IERR

!$OMP PARALLEL
!$OMP MASTER
    IERR = cublasGetVersion(ACC_CUBLASHANDLE,CUBLAS_VERSION)
    IF (IERR/=0) THEN
       CALL vtutor%bug("ERROR in cublasGetVersion", __FILE__, __LINE__)
    ENDIF
!$OMP END MASTER
!$OMP END PARALLEL
  END SUBROUTINE GET_CUBLAS_VERSION

  END MODULE mcublas


  SUBROUTINE ACC_DCOPY(N,DX,INCX,DY,INCY)
    USE mcublas
    REAL(q) :: DX(*),DY(*)
    INTEGER :: N,INCX,INCY
    ! local variables
    INTEGER(4) :: IERR

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL DCOPY(N,DX,INCX,DY,INCY)
       RETURN
    ENDIF
#ifdef use_kernels_construct
    IF (INCX==1.AND.INCY==1) THEN
!$ACC KERNELS PRESENT(DX(1:N),DY(1:N)) ASYNC(ACC_ASYNC_Q)
       DY(1:N)=DX(1:N)
!$ACC END KERNELS
       RETURN
    ENDIF
#endif
    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(DX,DY)
    IERR = cublasDcopy_v2(ACC_CUBLASHANDLE,N,DX,INCX,DY,INCY)
!$ACC END HOST_DATA
  END SUBROUTINE ACC_DCOPY


  SUBROUTINE ACC_ZCOPY(N,ZX,INCX,ZY,INCY)
    USE mcublas
    COMPLEX(q) :: ZX(*),ZY(*)
    INTEGER :: N,INCX,INCY
    ! local variables
    INTEGER(4) :: IERR

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL ZCOPY(N,ZX,INCX,ZY,INCY)
       RETURN
    ENDIF
#ifdef use_kernels_construct
    IF (INCX==1.AND.INCY==1) THEN
!$ACC KERNELS PRESENT(ZX(1:N),ZY(1:N)) ASYNC(ACC_ASYNC_Q)
       ZY(1:N)=ZX(1:N)
!$ACC END KERNELS
       RETURN
    ENDIF
#endif
    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(ZX,ZY)
    IERR = cublasZcopy_v2(ACC_CUBLASHANDLE,N,ZX,INCX,ZY,INCY)
!$ACC END HOST_DATA
  END SUBROUTINE ACC_ZCOPY


  SUBROUTINE ACC_DAXPY(N,DA,DX,INCX,DY,INCY)
    USE mcublas
    REAL(q) :: DX(*),DY(*)
    REAL(q) :: DA
    INTEGER :: N,INCX,INCY
    ! local variables
    INTEGER(4) :: IERR

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL DAXPY(N,DA,DX,INCX,DY,INCY)
       RETURN
    ENDIF

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(DX,DY)
    IERR = cublasDaxpy_v2(ACC_CUBLASHANDLE,N,DA,DX,INCX,DY,INCY)
!$ACC END HOST_DATA
  END SUBROUTINE ACC_DAXPY


  SUBROUTINE ACC_ZAXPY(N,ZA,ZX,INCX,ZY,INCY)
    USE mcublas
    COMPLEX(q) :: ZX(*),ZY(*)
    COMPLEX(q) :: ZA
    INTEGER :: N,INCX,INCY
    ! local variables
    INTEGER(4) :: IERR

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL ZAXPY(N,ZA,ZX,INCX,ZY,INCY)
       RETURN
    ENDIF

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(ZX,ZY)
    IERR = cublasZaxpy_v2(ACC_CUBLASHANDLE,N,ZA,ZX,INCX,ZY,INCY)
!$ACC END HOST_DATA
  END SUBROUTINE ACC_ZAXPY


  SUBROUTINE ACC_DSCAL(N,DA,DX,INCX)
    USE mcublas
    REAL(q) :: DX(*)
    REAL(q) :: DA
    INTEGER :: N,INCX
    ! local variables
    INTEGER(4) :: IERR

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL DSCAL(N,DA,DX,INCX)
       RETURN
    ENDIF

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(DX)
    IERR = cublasDscal_v2(ACC_CUBLASHANDLE,N,DA,DX,INCX)
!$ACC END HOST_DATA
  END SUBROUTINE ACC_DSCAL


  SUBROUTINE ACC_ZSCAL(N,ZA,ZX,INCX)
    USE mcublas
    COMPLEX(q) :: ZX(*)
    COMPLEX(q) :: ZA
    INTEGER :: N,INCX
    ! local variables
    INTEGER(4) :: IERR

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL ZSCAL(N,ZA,ZX,INCX)
       RETURN
    ENDIF

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(ZX)
    IERR = cublasZscal_v2(ACC_CUBLASHANDLE,N,ZA,ZX,INCX)
!$ACC END HOST_DATA
  END SUBROUTINE ACC_ZSCAL


  SUBROUTINE ACC_ZDSCAL(N,DA,ZX,INCX)
    USE mcublas
    COMPLEX(q) :: ZX(*)
    REAL(q) :: DA
    INTEGER :: N,INCX
    ! local variables
    INTEGER(4) :: IERR

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL ZDSCAL(N,DA,ZX,INCX)
       RETURN
    ENDIF

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(ZX)
    IERR = cublasZdscal_v2(ACC_CUBLASHANDLE,N,DA,ZX,INCX)
!$ACC END HOST_DATA
  END SUBROUTINE ACC_ZDSCAL


  REAL(q) FUNCTION ACC_DDOT(N,DX,INCX,DY,INCY)
    USE mcublas
    REAL(q) :: DX(*),DY(*)
    INTEGER :: N,INCX,INCY
    ! local variables
    INTEGER(4) :: IERR

    IF (.NOT.ACC_EXEC_ON) THEN
       ACC_DDOT = DDOT(N,DX,INCX,DY,INCY)
       RETURN
    ENDIF

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(DX,DY)
    IERR = cublasDdot_v2(ACC_CUBLASHANDLE,N,DX,INCX,DY,INCY,ACC_DDOT)
!$ACC END HOST_DATA
  END FUNCTION ACC_DDOT


  COMPLEX(q) FUNCTION ACC_ZDOTC(N,ZX,INCX,ZY,INCY)
    USE mcublas
    COMPLEX(q) :: ZX(*),ZY(*)
    INTEGER :: N,INCX,INCY
    ! local variables
    INTEGER(4) :: IERR

    IF (.NOT.ACC_EXEC_ON) THEN
       ACC_ZDOTC = ZDOTC(N,ZX,INCX,ZY,INCY)
       RETURN
    ENDIF

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(ZX,ZY)
    IERR = cublasZdotc_v2(ACC_CUBLASHANDLE,N,ZX,INCX,ZY,INCY,ACC_ZDOTC)
!$ACC END HOST_DATA
  END FUNCTION ACC_ZDOTC


  SUBROUTINE ACC_DGEMV(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    USE mcublas
    REAL(q) :: A(LDA,*),X(*),Y(*)
    REAL(q) :: ALPHA,BETA
    INTEGER :: M,N,LDA,INCX,INCY
    CHARACTER(1) :: MODE
    ! local variables
    REAL(q), ALLOCATABLE :: YTMP(:)
    INTEGER :: IDX,IDY,I
    LOGICAL :: LGEMM

    INTEGER(4) :: IERR

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL DGEMV(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
       RETURN
    ENDIF

    IF (MODE=='N'.OR.MODE=='n') THEN
       IDX=N*INCX
       IDY=M*INCY
    ELSE
       IDX=M*INCX
       IDY=N*INCY
    ENDIF

    IF ( CUBLAS_VERSION < 10100 ) THEN
! Use DGEMM or not? This should be made conditional upon a performance test.
    LGEMM=.FALSE.

    IF (LGEMM.AND.INCX==1.AND.INCY==1) THEN
       IF (ABS(BETA)>0._q) THEN
!$ACC ENTER DATA COPYIN(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
       ELSE
!$ACC ENTER DATA CREATE(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
       ENDIF
       IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
       IF (MODE=='N'.OR.MODE=='n') THEN
!$ACC HOST_DATA USE_DEVICE(Y,A,X)
          IERR = cublasDgemm_v2(ACC_CUBLASHANDLE,CUBLAS_OP_N,CUBLAS_OP_N,M,1,N,ALPHA,A,LDA,X,N,BETA,Y,IDY)
!$ACC END HOST_DATA
       ELSE
!$ACC HOST_DATA USE_DEVICE(Y,A,X)
          IERR = cublasDgemm_v2(ACC_CUBLASHANDLE,CHAR_TO_OP(MODE),CUBLAS_OP_N,N,1,M,ALPHA,A,LDA,X,M,BETA,Y,IDY)
!$ACC END HOST_DATA
       ENDIF
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
       RETURN
    ENDIF

    IF (MODE=='T'.OR.MODE=='t') THEN
       IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
       IF (ABS(BETA)>0._q) THEN
          ALLOCATE(YTMP(IDY))
!$ACC ENTER DATA COPYIN(Y(1:IDY)) CREATE(YTMP) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,X,YTMP)
          DO I=1,IDY,INCY
             IERR = cublasDdot_v2(ACC_CUBLASHANDLE,M,A(:,I),1,X,INCX,YTMP(I))
          ENDDO
!$ACC END HOST_DATA
!$ACC PARALLEL LOOP PRESENT(YTMP,Y) ASYNC(ACC_ASYNC_Q)
          DO I=1,IDY,INCY
             Y(I) = ALPHA* YTMP(I) + BETA*Y(I)
          ENDDO
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) DELETE(YTMP) ASYNC(ACC_ASYNC_Q)
          DEALLOCATE(YTMP)
       ELSE
          ALLOCATE(YTMP(IDY))
!$ACC ENTER DATA COPYIN(Y(1:IDY)) CREATE(YTMP) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,X,YTMP)
          DO I=1,IDY,INCY
             IERR = cublasDdot_v2(ACC_CUBLASHANDLE,M,A(:,I),1,X,INCX,YTMP(I))
          ENDDO
!$ACC END HOST_DATA
!$ACC PARALLEL LOOP PRESENT(YTMP,Y) ASYNC(ACC_ASYNC_Q)
          DO I=1,IDY,INCY
             Y(I) = ALPHA* YTMP(I)
          ENDDO
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) DELETE(YTMP) ASYNC(ACC_ASYNC_Q)
          DEALLOCATE(YTMP)
       ENDIF
       RETURN
    ENDIF
    ENDIF

    IF (ABS(BETA)>0._q) THEN
!$ACC ENTER DATA COPYIN(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
    ELSE
!$ACC ENTER DATA CREATE(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(Y,A,X)
    IERR = cublasDgemv_v2(ACC_CUBLASHANDLE,CHAR_TO_OP(MODE),M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!$ACC END HOST_DATA
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE ACC_DGEMV


  SUBROUTINE ACC_ZGEMV(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    USE mcublas
    COMPLEX(q) :: A(LDA,*),X(*),Y(*)
    COMPLEX(q) :: ALPHA,BETA
    INTEGER :: M,N,LDA,INCX,INCY
    CHARACTER(1) :: MODE
    ! local variables
    COMPLEX(q), ALLOCATABLE :: YTMP(:)
    INTEGER :: IDX,IDY,I
    LOGICAL :: LGEMM

    INTEGER(4) :: IERR

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL ZGEMV(MODE,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
       RETURN
    ENDIF

    IF (MODE=='N'.OR.MODE=='n') THEN
       IDX=N*INCX
       IDY=M*INCY
    ELSE
       IDX=M*INCX
       IDY=N*INCY
    ENDIF

    IF ( CUBLAS_VERSION < 10100 ) THEN
! Use ZGEMM or not? This should be made conditional upon a performance test.
    LGEMM=.FALSE.

    IF (LGEMM.AND.INCX==1.AND.INCY==1) THEN
       IF (ABS(BETA)>0._q) THEN
!$ACC ENTER DATA COPYIN(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
       ELSE
!$ACC ENTER DATA CREATE(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
       ENDIF
       IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
       IF (MODE=='N'.OR.MODE=='n') THEN
!$ACC HOST_DATA USE_DEVICE(Y,A,X)
          IERR = cublasZgemm_v2(ACC_CUBLASHANDLE,CUBLAS_OP_N,CUBLAS_OP_N,M,1,N,ALPHA,A,LDA,X,N,BETA,Y,IDY)
!$ACC END HOST_DATA
       ELSE
!$ACC HOST_DATA USE_DEVICE(Y,A,X)
          IERR = cublasZgemm_v2(ACC_CUBLASHANDLE,CHAR_TO_OP(MODE),CUBLAS_OP_N,N,1,M,ALPHA,A,LDA,X,M,BETA,Y,IDY)
!$ACC END HOST_DATA
       ENDIF
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
       RETURN
    ENDIF

    IF (MODE=='C'.OR.MODE=='c') THEN
       IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
       IF (ABS(BETA)>0._q) THEN
          ALLOCATE(YTMP(IDY))
!$ACC ENTER DATA COPYIN(Y(1:IDY)) CREATE(YTMP) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,X,YTMP)
          DO I=1,IDY,INCY
             IERR = cublasZdotc_v2(ACC_CUBLASHANDLE,M,A(:,I),1,X,INCX,YTMP(I))
          ENDDO
!$ACC END HOST_DATA
!$ACC PARALLEL LOOP PRESENT(YTMP,Y) ASYNC(ACC_ASYNC_Q)
          DO I=1,IDY,INCY
             Y(I) = ALPHA* YTMP(I) + BETA*Y(I)
          ENDDO
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) DELETE(YTMP) ASYNC(ACC_ASYNC_Q)
          DEALLOCATE(YTMP)
       ELSE
          ALLOCATE(YTMP(IDY))
!$ACC ENTER DATA COPYIN(Y(1:IDY)) CREATE(YTMP) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,X,YTMP)
          DO I=1,IDY,INCY
             IERR = cublasZdotc_v2(ACC_CUBLASHANDLE,M,A(:,I),1,X,INCX,YTMP(I))
          ENDDO
!$ACC END HOST_DATA
!$ACC PARALLEL LOOP PRESENT(YTMP,Y) ASYNC(ACC_ASYNC_Q)
          DO I=1,IDY,INCY
             Y(I) = ALPHA* YTMP(I)
          ENDDO
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) DELETE(YTMP) ASYNC(ACC_ASYNC_Q)
          DEALLOCATE(YTMP)
       ENDIF
       RETURN
    ENDIF
    IF (MODE=='T'.OR.MODE=='t') THEN
       IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
       IF (ABS(BETA)>0._q) THEN
          ALLOCATE(YTMP(IDY))
!$ACC ENTER DATA COPYIN(Y(1:IDY)) CREATE(YTMP) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,X,YTMP)
          DO I=1,IDY,INCY
             IERR = cublasZdotu_v2(ACC_CUBLASHANDLE,M,A(:,I),1,X,INCX,YTMP(I))
          ENDDO
!$ACC END HOST_DATA
!$ACC PARALLEL LOOP PRESENT(YTMP,Y) ASYNC(ACC_ASYNC_Q)
          DO I=1,IDY,INCY
             Y(I) = ALPHA* YTMP(I) + BETA*Y(I)
          ENDDO
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) DELETE(YTMP) ASYNC(ACC_ASYNC_Q)
          DEALLOCATE(YTMP)
       ELSE
          ALLOCATE(YTMP(IDY))
!$ACC ENTER DATA COPYIN(Y(1:IDY)) CREATE(YTMP) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,X,YTMP)
          DO I=1,IDY,INCY
             IERR = cublasZdotu_v2(ACC_CUBLASHANDLE,M,A(:,I),1,X,INCX,YTMP(I))
          ENDDO
!$ACC END HOST_DATA
!$ACC PARALLEL LOOP PRESENT(YTMP,Y) ASYNC(ACC_ASYNC_Q)
          DO I=1,IDY,INCY
             Y(I) = ALPHA* YTMP(I)
          ENDDO
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) DELETE(YTMP) ASYNC(ACC_ASYNC_Q)
          DEALLOCATE(YTMP)
       ENDIF
       RETURN
    ENDIF
    ENDIF

    IF (ABS(BETA)>0._q) THEN
!$ACC ENTER DATA COPYIN(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
    ELSE
!$ACC ENTER DATA CREATE(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
    ENDIF
    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(Y,A,X)
    IERR = cublasZgemv_v2(ACC_CUBLASHANDLE,CHAR_TO_OP(MODE),M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!$ACC END HOST_DATA
!$ACC EXIT DATA COPYOUT(Y(1:IDY)) ASYNC(ACC_ASYNC_Q)
  END SUBROUTINE ACC_ZGEMV


  SUBROUTINE ACC_DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    USE mcublas
    REAL(q) :: A(LDA,*),B(LDB,*),C(LDC,*)
    REAL(q) :: ALPHA,BETA
    INTEGER :: M,N,K,LDA,LDB,LDC
    CHARACTER(1) :: TRANSA,TRANSB
    ! local variables
    INTEGER(4) :: IERR

    IF (N<1) RETURN

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
       RETURN
    ENDIF

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(A,B,C)
    IERR = cublasDgemm_v2(ACC_CUBLASHANDLE,CHAR_TO_OP(TRANSA),CHAR_TO_OP(TRANSB),M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!$ACC END HOST_DATA
  END SUBROUTINE ACC_DGEMM


  SUBROUTINE ACC_ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    USE mcublas
    COMPLEX(q) :: A(LDA,*),B(LDB,*),C(LDC,*)
    COMPLEX(q) :: ALPHA,BETA
    INTEGER :: M,N,K,LDA,LDB,LDC
    CHARACTER(1) :: TRANSA,TRANSB
    ! local variables
    INTEGER(4) :: IERR

    IF (N<1) RETURN

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
       RETURN
    ENDIF

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(A,B,C)
    IERR = cublasZgemm_v2(ACC_CUBLASHANDLE,CHAR_TO_OP(TRANSA),CHAR_TO_OP(TRANSB),M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!$ACC END HOST_DATA
  END SUBROUTINE ACC_ZGEMM


  SUBROUTINE ACC_DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
    USE mcublas
    REAL(q) :: A(LDA,*),B(LDB,*)
    REAL(q) :: ALPHA
    INTEGER :: LDA,LDB,M,N
    CHARACTER(1) :: DIAG,SIDE,TRANSA,UPLO
    ! local variables
    INTEGER(4) :: IERR

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
       RETURN
    ENDIF

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(A,B)
    IERR = cublasDtrmm_v2(ACC_CUBLASHANDLE,CHAR_TO_SIDE(SIDE),CHAR_TO_FILL(UPLO),CHAR_TO_OP(TRANSA),CHAR_TO_DIAG(DIAG),M,N,ALPHA,A,LDA,B,LDB,B,LDB)
!$ACC END HOST_DATA
  END SUBROUTINE ACC_DTRMM


  SUBROUTINE ACC_ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
    USE mcublas
    COMPLEX(q) :: A(LDA,*),B(LDB,*)
    COMPLEX(q) :: ALPHA
    INTEGER :: LDA,LDB,M,N
    CHARACTER(1) :: DIAG,SIDE,TRANSA,UPLO
    ! local variables
    INTEGER(4) :: IERR

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
       RETURN
    ENDIF

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))
!$ACC HOST_DATA USE_DEVICE(A,B)
    IERR = cublasZtrmm_v2(ACC_CUBLASHANDLE,CHAR_TO_SIDE(SIDE),CHAR_TO_FILL(UPLO),CHAR_TO_OP(TRANSA),CHAR_TO_DIAG(DIAG),M,N,ALPHA,A,LDA,B,LDB,B,LDB)
!$ACC END HOST_DATA
  END SUBROUTINE ACC_ZTRMM


!***********************************************************************
!
! This module only holds a few variables related to the use of
! cusolverDn and a subroutine to initialize and destroy the
! cusolverDnHandle ACC_CUSOLVERLHANDLE.
!
! The actual cusolverDn calls are wrapped in ACC_* subroutines that
! live outside of mcusolver (see immediately below the aforementioned
! module).
!
!***********************************************************************

  MODULE mcusolver
  USE prec

  USE cusolverDn, ONLY : cusolverDnHandle, &
                         cusolverDnCreate, &
                         cusolverDnDestroy, &
                         cusolverDnSetStream, &
                         cusolverDnDpotrf_bufferSize, &
                         cusolverDnZpotrf_bufferSize, &
                         cusolverDnDgetrf_bufferSize, &
                         cusolverDnZgetrf_bufferSize, &
                         cusolverDnDsygvd_bufferSize, &
                         cusolverDnZhegvd_bufferSize, &
                         cusolverDnDsyevd_bufferSize, &
                         cusolverDnZheevd_bufferSize, &
                         cusolverDnDpotrf, &
                         cusolverDnZpotrf, &
                         cusolverDnDgetrf, &
                         cusolverDnDgetrs, &
                         cusolverDnZgetrf, &
                         cusolverDnZgetrs, &
                         cusolverDnDsyevd, &
                         cusolverDnZheevd, &
                         CUSOLVER_EIG_MODE_NOVECTOR, &
                         CUSOLVER_EIG_MODE_VECTOR, &
                         CUSOLVER_EIG_TYPE_1, &
                         CUSOLVER_EIG_TYPE_2, &
                         CUSOLVER_EIG_TYPE_3

  USE cublas, ONLY :     CUBLAS_FILL_MODE_UPPER, &
                         CUBLAS_FILL_MODE_LOWER, &
                         CUBLAS_OP_N, &
                         CUBLAS_OP_T, &
                         CUBLAS_OP_C

  USE openacc,  ONLY : acc_get_cuda_stream,c_devloc

  USE mopenacc, ONLY : ACC_EXEC_ALLOW,ACC_EXEC_ON,ACC_ASYNC_Q

  USE string, ONLY: str
  USE tutor, ONLY: vtutor

  IMPLICIT NONE

  TYPE(cusolverDnHandle) :: ACC_CUSOLVERHANDLE
!$OMP THREADPRIVATE(ACC_CUSOLVERHANDLE)

  INTERFACE
  INTEGER(c_int) FUNCTION cusolverDnZhegvd(handle,itype,jobz,uplo,n,A,lda,B,ldb,W,work,lwork,devInfo) BIND(C, name="cusolverDnZhegvd")
    USE ISO_C_BINDING
    USE cusolverDn, ONLY : cusolverDnHandle
    IMPLICIT NONE
    TYPE(cusolverDnHandle),VALUE::handle
    INTEGER(c_int),VALUE::itype
    INTEGER(c_int),VALUE::jobz
    INTEGER(c_int),VALUE::uplo
    INTEGER(c_int),VALUE::n
    !dir$ ignore_tkr(trd) A
    DOUBLE COMPLEX, DIMENSION(*) :: A
    INTEGER(c_int),VALUE::lda
    !dir$ ignore_tkr(trd) B
    DOUBLE COMPLEX, DIMENSION(*) :: B
    INTEGER(c_int),VALUE::ldb
    !dir$ ignore_tkr(trd) W
    DOUBLE PRECISION, DIMENSION(*)::W
    !dir$ ignore_tkr(trd) work
    DOUBLE COMPLEX, DIMENSION(*) :: work
    INTEGER(c_int),VALUE::lwork
    !dir$ ignore_tkr(trd) devInfo
    INTEGER, DIMENSION(1)::devInfo
  END FUNCTION
  END INTERFACE

  INTERFACE
  INTEGER(c_int) FUNCTION cusolverDnDsygvd(handle,itype,jobz,uplo,n,A,lda,B,ldb,W,work,lwork,devInfo) BIND(C, name="cusolverDnDsygvd")
    USE ISO_C_BINDING
    USE cusolverDn, ONLY : cusolverDnHandle
    IMPLICIT NONE
    TYPE(cusolverDnHandle),VALUE::handle
    INTEGER(c_int),VALUE::itype
    INTEGER(c_int),VALUE::jobz
    INTEGER(c_int),VALUE::uplo
    INTEGER(c_int),VALUE::n
    !dir$ ignore_tkr(trd) A
    DOUBLE PRECISION, DIMENSION(*) :: A
    INTEGER(c_int),VALUE::lda
    !dir$ ignore_tkr(trd) B
    DOUBLE PRECISION, DIMENSION(*) :: B
    INTEGER(c_int),VALUE::ldb
    !dir$ ignore_tkr(trd) W
    DOUBLE PRECISION, DIMENSION(*)::W
    !dir$ ignore_tkr(trd) work
    DOUBLE PRECISION, DIMENSION(*) :: work
    INTEGER(c_int),VALUE::lwork
    !dir$ ignore_tkr(trd) devInfo
    INTEGER, DIMENSION(1)::devInfo
  END FUNCTION
  END INTERFACE

  CONTAINS

  SUBROUTINE ACC_INIT_CUSOLVER()
    INTEGER(4) :: IERR

    PROFILING_START('acc_init_cusolver')

    IF (.NOT.ACC_EXEC_ALLOW) THEN
       PROFILING_STOP('acc_init_cusolver')
       RETURN
    ENDIF

!$OMP PARALLEL
    IERR = cusolverDnCreate(ACC_CUSOLVERHANDLE)
    IF (IERR/=0) THEN
       CALL vtutor%error("ERROR setting up cusolverHandle")
    ENDIF
!$OMP END PARALLEL

    PROFILING_STOP('acc_init_cusolver')

  END SUBROUTINE ACC_INIT_CUSOLVER


  SUBROUTINE ACC_FINISH_CUSOLVER()
     INTEGER(4) :: IERR

     IF (.NOT.ACC_EXEC_ALLOW) RETURN

!$OMP PARALLEL
     IERR = cusolverDnDestroy(ACC_CUSOLVERHANDLE)
     IF (IERR/=0) THEN
        CALL vtutor%error("ERROR releasing cusolverHandle")
     ENDIF
!$OMP END PARALLEL
  END SUBROUTINE ACC_FINISH_CUSOLVER

  END MODULE mcusolver


  SUBROUTINE ACC_DPOTRF(UPLO,N,A,LDA,INFO)
    USE mcusolver
    IMPLICIT NONE
    REAL(q) :: A(LDA,*)
    INTEGER :: N,LDA,INFO
    CHARACTER(1) :: UPLO
    ! local variables
    REAL(q), ALLOCATABLE :: WORK(:)
    INTEGER :: LWORK
    INTEGER(4) :: IERR,UPLO_,INFO_

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL DPOTRF(UPLO,N,A,LDA,INFO)
       RETURN
    ENDIF
#ifdef magma_dpotrf
    CALL magmaf_dpotrf_gpu(UPLO,N,c_devloc(A),LDA,INFO)
    RETURN
#endif
    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

!$ACC HOST_DATA USE_DEVICE(A)
    IERR = IERR + cusolverDnDpotrf_bufferSize(ACC_CUSOLVERHANDLE,UPLO_,N,A,LDA,LWORK)
!$ACC END HOST_DATA

    ALLOCATE(WORK(LWORK))
!$ACC ENTER DATA CREATE(WORK,INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,WORK,INFO_)
    IERR = IERR + cusolverDnDpotrf(ACC_CUSOLVERHANDLE,UPLO_,N,A,LDA,WORK,LWORK,INFO_)
!$ACC END HOST_DATA

!$ACC EXIT DATA DELETE(WORK) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA COPYOUT(INFO_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK)

    INFO=INFO_*100+IERR
  END SUBROUTINE ACC_DPOTRF


  SUBROUTINE ACC_ZPOTRF(UPLO,N,A,LDA,INFO)
    USE mcusolver
    IMPLICIT NONE
    COMPLEX(q) :: A(LDA,*)
    INTEGER :: N,LDA,INFO
    CHARACTER(1) :: UPLO
    ! local variables
    COMPLEX(q), ALLOCATABLE :: WORK(:)
    INTEGER :: LWORK
    INTEGER(4) :: IERR,UPLO_,INFO_

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL ZPOTRF(UPLO,N,A,LDA,INFO)
       RETURN
    ENDIF
#ifdef magma_zpotrf
    CALL magmaf_zpotrf_gpu(UPLO,N,c_devloc(A),LDA,INFO)
    RETURN
#endif
    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

!$ACC HOST_DATA USE_DEVICE(A)
    IERR = IERR + cusolverDnZpotrf_bufferSize(ACC_CUSOLVERHANDLE,UPLO_,N,A,LDA,LWORK)
!$ACC END HOST_DATA

    ALLOCATE(WORK(LWORK))
!$ACC ENTER DATA CREATE(WORK,INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,WORK,INFO_)
    IERR = IERR + cusolverDnZpotrf(ACC_CUSOLVERHANDLE,UPLO_,N,A,LDA,WORK,LWORK,INFO_)
!$ACC END HOST_DATA

!$ACC EXIT DATA DELETE(WORK) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA COPYOUT(INFO_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK)

    INFO=INFO_*100+IERR
  END SUBROUTINE ACC_ZPOTRF


  SUBROUTINE ACC_DGETRF(M,N,A,LDA,IPIV,INFO)
    USE mcusolver
    IMPLICIT NONE
    REAL(q) :: A(LDA,*)
    INTEGER :: M,N,LDA,INFO,IPIV(*)
    ! local variables
    REAL(q), ALLOCATABLE :: WORK(:)
    INTEGER :: LWORK
    INTEGER(4) :: IERR,INFO_

    IF (.NOT.ACC_EXEC_ON) THEN
      CALL DGETRF(M,N,A,LDA,IPIV,INFO)
      RETURN
    ENDIF

    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

!$ACC HOST_DATA USE_DEVICE(A)
    IERR = IERR + cusolverDnDgetrf_bufferSize(ACC_CUSOLVERHANDLE,M,N,A,LDA,LWORK)
!$ACC END HOST_DATA

    ALLOCATE(WORK(LWORK))
!$ACC ENTER DATA CREATE(WORK,INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,WORK,IPIV,INFO_)
    IERR = IERR + cusolverDnDgetrf(ACC_CUSOLVERHANDLE,M,N,A,LDA,WORK,IPIV,INFO_)
!$ACC END HOST_DATA

!$ACC EXIT DATA DELETE(WORK) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA COPYOUT(INFO_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK)

    INFO=INFO_*100+IERR
  END SUBROUTINE ACC_DGETRF


  SUBROUTINE ACC_ZGETRF(M,N,A,LDA,IPIV,INFO)
    USE mcusolver
    IMPLICIT NONE
    COMPLEX(q) :: A(LDA,*)
    INTEGER :: M,N,LDA,INFO,IPIV(*)
    ! local variables
    COMPLEX(q), ALLOCATABLE :: WORK(:)
    INTEGER :: LWORK
    INTEGER(4) :: IERR,INFO_

    IF (.NOT.ACC_EXEC_ON) THEN
      CALL ZGETRF(M,N,A,LDA,IPIV,INFO)
      RETURN
    ENDIF

    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

!$ACC HOST_DATA USE_DEVICE(A)
    IERR = IERR + cusolverDnZgetrf_bufferSize(ACC_CUSOLVERHANDLE,M,N,A,LDA,LWORK)
!$ACC END HOST_DATA

    ALLOCATE(WORK(LWORK))
!$ACC ENTER DATA CREATE(WORK,INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,WORK,IPIV,INFO_)
    IERR = IERR + cusolverDnZgetrf(ACC_CUSOLVERHANDLE,M,N,A,LDA,WORK,IPIV,INFO_)
!$ACC END HOST_DATA

!$ACC EXIT DATA DELETE(WORK) ASYNC(ACC_ASYNC_Q)
!$ACC EXIT DATA COPYOUT(INFO_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK)

    INFO=INFO_*100+IERR
  END SUBROUTINE ACC_ZGETRF


  SUBROUTINE ACC_DGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
    USE mcusolver
    IMPLICIT NONE
    REAL(q) :: A(LDA,*),B(LDB,*)
    INTEGER :: N,NRHS,LDA,LDB,INFO,IPIV(*)
    CHARACTER(1) :: TRANS
    ! local variables
    INTEGER(4) :: IERR,TRANS_,INFO_

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL DGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
       RETURN
    ENDIF

    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (TRANS=='N'.OR.TRANS=='n') THEN
       TRANS_=CUBLAS_OP_N
    ELSE
       TRANS_=CUBLAS_OP_T
    ENDIF

!$ACC HOST_DATA USE_DEVICE(A,IPIV,B,INFO_)
    IERR = IERR + cusolverDnDgetrs(ACC_CUSOLVERHANDLE,TRANS_,N,NRHS,A,LDA,IPIV,B,LDB,INFO_)
!$ACC END HOST_DATA

!$ACC EXIT DATA COPYOUT(INFO_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)

    INFO=INFO_*100+IERR
  END SUBROUTINE ACC_DGETRS


  SUBROUTINE ACC_ZGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
    USE mcusolver
    IMPLICIT NONE
    COMPLEX(q) :: A(LDA,*),B(LDB,*)
    INTEGER :: N,NRHS,LDA,LDB,INFO,IPIV(*)
    CHARACTER(1) :: TRANS
    ! local variables
    INTEGER(4) :: IERR,TRANS_,INFO_

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL ZGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
       RETURN
    ENDIF

    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (TRANS=='N'.OR.TRANS=='n') THEN
       TRANS_=CUBLAS_OP_N
    ELSEIF (TRANS=='T'.OR.TRANS=='t') THEN
       TRANS_=CUBLAS_OP_T
    ELSE
       TRANS_=CUBLAS_OP_C
    ENDIF

!$ACC HOST_DATA USE_DEVICE(A,IPIV,B,INFO_)
    IERR = IERR + cusolverDnZgetrs(ACC_CUSOLVERHANDLE,TRANS_,N,NRHS,A,LDA,IPIV,B,LDB,INFO_)
!$ACC END HOST_DATA

!$ACC EXIT DATA COPYOUT(INFO_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)

    INFO=INFO_*100+IERR
  END SUBROUTINE ACC_ZGETRS


  SUBROUTINE ACC_DTRTRI(UPLO,DIAG,N,A,LDA,INFO)
    USE mcublas
    IMPLICIT NONE
    REAL(q) :: A(LDA,*)
    INTEGER :: N,LDA,INFO
    CHARACTER(1) :: UPLO,DIAG
    ! local variables
    REAL(q), ALLOCATABLE :: WORK(:)
    INTEGER :: I,J
    INTEGER(4) :: IERR

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL DTRTRI(UPLO,DIAG,N,A,LDA,INFO)
       RETURN
    ENDIF
#ifdef magma_dtrtri
    CALL magmaf_dtrtri_gpu(UPLO,DIAG,N,c_devloc(A),LDA,INFO)
    RETURN
#endif
    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       ALLOCATE(WORK(N*N))
!$ACC ENTER DATA CREATE(WORK) ASYNC(ACC_ASYNC_Q)
!$ACC KERNELS PRESENT(WORK) ASYNC(ACC_ASYNC_Q)
       WORK=0
       DO I=1,N
          WORK(I+(I-1)*N)=1
       ENDDO
!$ACC END KERNELS

!$ACC HOST_DATA USE_DEVICE(A,WORK)
       IERR = IERR + cublasDtrsm_v2(ACC_CUBLASHANDLE,CUBLAS_SIDE_LEFT,CUBLAS_FILL_MODE_UPPER, &
                                    CUBLAS_OP_N,CUBLAS_DIAG_NON_UNIT,N,N,done,A,LDA,WORK,N)
!$ACC END HOST_DATA

!$ACC KERNELS PRESENT(A,WORK) ASYNC(ACC_ASYNC_Q)
       DO I=1,N
          DO J=1,N
             A(J,I)=WORK(J+(I-1)*N)
          ENDDO
       ENDDO
!$ACC END KERNELS

!$ACC EXIT DATA DELETE(WORK) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
       DEALLOCATE(WORK)
    ENDIF

    INFO=IERR
  END SUBROUTINE ACC_DTRTRI


  SUBROUTINE ACC_ZTRTRI(UPLO,DIAG,N,A,LDA,INFO)
    USE mcublas
    IMPLICIT NONE
    COMPLEX(q) :: A(LDA,*)
    INTEGER :: N,LDA,INFO
    CHARACTER(1) :: UPLO,DIAG
    ! local variables
    COMPLEX(q), ALLOCATABLE :: WORK(:)
    INTEGER :: I,J
    INTEGER(4) :: IERR

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL ZTRTRI(UPLO,DIAG,N,A,LDA,INFO)
       RETURN
    ENDIF
#ifdef magma_ztrtri
    CALL magmaf_ztrtri_gpu(UPLO,DIAG,N,c_devloc(A),LDA,INFO)
    RETURN
#endif

    IERR = cublasSetStream(ACC_CUBLASHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       ALLOCATE(WORK(N*N))
!$ACC ENTER DATA CREATE(WORK) ASYNC(ACC_ASYNC_Q)
!$ACC KERNELS PRESENT(WORK) ASYNC(ACC_ASYNC_Q)
       WORK=0
       DO I=1,N
          WORK(I+(I-1)*N)=1
       ENDDO
!$ACC END KERNELS

!$ACC HOST_DATA USE_DEVICE(A,WORK)
       IERR = IERR + cublasZtrsm_v2(ACC_CUBLASHANDLE,CUBLAS_SIDE_LEFT,CUBLAS_FILL_MODE_UPPER, &
                                    CUBLAS_OP_N,CUBLAS_DIAG_NON_UNIT,N,N,zone,A,LDA,WORK,N)
!$ACC END HOST_DATA

!$ACC KERNELS PRESENT(A,WORK) ASYNC(ACC_ASYNC_Q)
       DO I=1,N
          DO J=1,N
             A(J,I)=WORK(J+(I-1)*N)
          ENDDO
       ENDDO
!$ACC END KERNELS

!$ACC EXIT DATA DELETE(WORK) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
       DEALLOCATE(WORK)
    ENDIF

    INFO=IERR
  END SUBROUTINE ACC_ZTRTRI


  SUBROUTINE ACC_DSYGV(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,INFO)
    USE mcusolver
    IMPLICIT NONE
    REAL(q) :: A(LDA,*),B(LDB,*),W(*),WORK(*)
    INTEGER :: ITYPE,N,LDA,LDB,LWORK,INFO
    CHARACTER(1) :: JOBZ,UPLO
    ! local variables
    REAL(q), ALLOCATABLE :: WORK_(:)
    INTEGER :: ITYPE_,JOBZ_,UPLO_,LWORK_
    INTEGER(4) :: IERR,INFO_

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL DSYGV(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,INFO)
       RETURN
    ENDIF
#ifdef magma_dsygv
!   CALL magmaf_dsygv_gpu
    RETURN
#endif
    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (ITYPE==1) THEN
       ITYPE_=CUSOLVER_EIG_TYPE_1
    ELSEIF (ITYPE==2) THEN
       ITYPE_=CUSOLVER_EIG_TYPE_2
    ELSE
       ITYPE_=CUSOLVER_EIG_TYPE_3
    ENDIF

    IF (JOBZ=='N'.OR.JOBZ=='n') THEN
       JOBZ_=CUSOLVER_EIG_MODE_NOVECTOR
    ELSE
       JOBZ_=CUSOLVER_EIG_MODE_VECTOR
    ENDIF

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

!$ACC ENTER DATA CREATE(W(1:N)) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,B,W)
    IERR = IERR + cusolverDnDsygvd_bufferSize(ACC_CUSOLVERHANDLE,ITYPE_,JOBZ_,UPLO_,N,A,LDA,B,LDB,W,LWORK_)
!$ACC END HOST_DATA
    ALLOCATE(WORK_(LWORK_))
!$ACC ENTER DATA CREATE(WORK_,INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,B,W,WORK_,INFO_)
    IERR = IERR + cusolverDnDsygvd(ACC_CUSOLVERHANDLE,ITYPE_,JOBZ_,UPLO_,N,A,LDA,B,LDB,W,WORK_,LWORK_,INFO_)
!$ACC END HOST_DATA
!$ACC EXIT DATA COPYOUT(W(1:N),INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC EXIT DATA DELETE(WORK_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK_)

    INFO=INFO_*100+IERR
  END SUBROUTINE ACC_DSYGV


  SUBROUTINE ACC_ZHEGV(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,RWORK,INFO)
    USE mcusolver
    IMPLICIT NONE
    COMPLEX(q) :: A(LDA,*),B(LDB,*),WORK(*)
    REAL(q) :: W(*),RWORK(*)
    INTEGER :: ITYPE,N,LDA,LDB,LWORK,INFO
    CHARACTER(1) :: JOBZ,UPLO
    ! local variables
    COMPLEX(q), ALLOCATABLE :: WORK_(:)
    INTEGER :: ITYPE_,JOBZ_,UPLO_,LWORK_
    INTEGER(4) :: IERR,INFO_

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL ZHEGV(ITYPE,JOBZ,UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,RWORK,INFO)
       RETURN
    ENDIF
#ifdef magma_zhegv
!   CALL magmaf_zhegv_gpu
    RETURN
#endif
    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (ITYPE==1) THEN
       ITYPE_=CUSOLVER_EIG_TYPE_1
    ELSEIF (ITYPE==2) THEN
       ITYPE_=CUSOLVER_EIG_TYPE_2
    ELSE
       ITYPE_=CUSOLVER_EIG_TYPE_3
    ENDIF

    IF (JOBZ=='N'.OR.JOBZ=='n') THEN
       JOBZ_=CUSOLVER_EIG_MODE_NOVECTOR
    ELSE
       JOBZ_=CUSOLVER_EIG_MODE_VECTOR
    ENDIF

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

!$ACC ENTER DATA CREATE(W(1:N)) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,B,W)
       IERR = IERR + cusolverDnZhegvd_bufferSize(ACC_CUSOLVERHANDLE,ITYPE_,JOBZ_,UPLO_,N,A,LDA,B,LDB,W,LWORK_)
!$ACC END HOST_DATA

       ALLOCATE(WORK_(LWORK_))
!$ACC ENTER DATA CREATE(WORK_,INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,B,W,WORK_,INFO_)
    IERR = IERR + cusolverDnZhegvd(ACC_CUSOLVERHANDLE,ITYPE_,JOBZ_,UPLO_,N,A,LDA,B,LDB,W,WORK_,LWORK_,INFO_)
!$ACC END HOST_DATA
!$ACC EXIT DATA COPYOUT(W(1:N),INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC EXIT DATA DELETE(WORK_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK_)

    INFO=INFO_*100+IERR
  END SUBROUTINE ACC_ZHEGV


  SUBROUTINE ACC_DSYEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
    USE mcusolver
    IMPLICIT NONE
    REAL(q) :: A(LDA,*),W(*),WORK(*)
    INTEGER :: N,LDA,LWORK,INFO
    CHARACTER(1) :: JOBZ,UPLO
    ! local variables
    REAL(q), ALLOCATABLE :: WORK_(:)
    INTEGER :: JOBZ_,UPLO_,LWORK_
    INTEGER(4) :: IERR,INFO_

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL DSYEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
       RETURN
    ENDIF
#ifdef magma_dsyev
!   CALL magmaf_dsyev_gpu
    RETURN
#endif
    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (JOBZ=='N'.OR.JOBZ=='n') THEN
       JOBZ_=CUSOLVER_EIG_MODE_NOVECTOR
    ELSE
       JOBZ_=CUSOLVER_EIG_MODE_VECTOR
    ENDIF

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

!$ACC ENTER DATA CREATE(W(1:N)) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,W)
    IERR = IERR + cusolverDnDsyevd_bufferSize(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,LWORK_)
!$ACC END HOST_DATA
    ALLOCATE(WORK_(LWORK_))
!$ACC ENTER DATA CREATE(WORK_,INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,W,WORK_,INFO_)
    IERR = IERR + cusolverDnDsyevd(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,WORK_,LWORK_,INFO_)
!$ACC END HOST_DATA
!$ACC EXIT DATA COPYOUT(W(1:N),INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC EXIT DATA DELETE(WORK_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK_)

    INFO=INFO_*100+IERR
  END SUBROUTINE ACC_DSYEV


  SUBROUTINE ACC_ZHEEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)
    USE mcusolver
    IMPLICIT NONE
    COMPLEX(q) :: A(LDA,*),WORK(*)
    REAL(q) :: W(*),RWORK(*)
    INTEGER :: N,LDA,LWORK,INFO
    CHARACTER(1) :: JOBZ,UPLO
    ! local variables
    COMPLEX(q), ALLOCATABLE :: WORK_(:)
    INTEGER :: JOBZ_,UPLO_,LWORK_
    INTEGER(4) :: IERR,INFO_

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL ZHEEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,INFO)
       RETURN
    ENDIF
#ifdef magma_zheev
!   CALL magmaf_zheev_gpu
    RETURN
#endif
    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (JOBZ=='N'.OR.JOBZ=='n') THEN
       JOBZ_=CUSOLVER_EIG_MODE_NOVECTOR
    ELSE
       JOBZ_=CUSOLVER_EIG_MODE_VECTOR
    ENDIF

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

!$ACC ENTER DATA CREATE(W(1:N)) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,W)
    IERR = IERR + cusolverDnZheevd_bufferSize(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,LWORK_)
!$ACC END HOST_DATA
    ALLOCATE(WORK_(LWORK_))
!$ACC ENTER DATA CREATE(WORK_,INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,W,WORK_,INFO_)
    IERR = IERR + cusolverDnZheevd(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,WORK_,LWORK_,INFO_)
!$ACC END HOST_DATA
!$ACC EXIT DATA COPYOUT(W(1:N),INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC EXIT DATA DELETE(WORK_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK_)

    INFO=INFO_*100+IERR
  END SUBROUTINE ACC_ZHEEV


  SUBROUTINE ACC_DSYEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)
    USE mcusolver
    IMPLICIT NONE
    REAL(q) :: A(LDA,*),W(*),WORK(*),Z(LDZ,*),VL,VU,ABSTOL
    INTEGER :: IWORK(*),IFAIL(*),N,LDA,IL,IU,M,LDZ,LWORK,INFO
    CHARACTER(1) :: JOBZ,RANGE,UPLO
    ! local variables
    REAL(q), ALLOCATABLE :: WORK_(:)
!$ACC DECLARE DEVICE_RESIDENT(WORK_)
    INTEGER :: JOBZ_,UPLO_,LWORK_,I
    INTEGER(4) :: IERR,INFO_

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL DSYEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO)
       RETURN
    ENDIF
#ifdef magma_dsyevx
!   CALL magmaf_dsyevx_gpu
    RETURN
#endif
    IF (RANGE/='A'.AND.RANGE/='a') THEN
       CALL vtutor%error("ACC_DSYEVX: ERROR: RANGE= " // RANGE // " not implemented")
    ENDIF

    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (JOBZ=='N'.OR.JOBZ=='n') THEN
       JOBZ_=CUSOLVER_EIG_MODE_NOVECTOR
    ELSE
       JOBZ_=CUSOLVER_EIG_MODE_VECTOR
    ENDIF

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

!$ACC ENTER DATA CREATE(W(1:N)) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,W)
    IERR = IERR + cusolverDnDsyevd_bufferSize(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,LWORK_)
!$ACC END HOST_DATA
    ALLOCATE(WORK_(LWORK_))
!$ACC ENTER DATA CREATE(INFO_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,W,WORK_,INFO_)
    IERR = IERR + cusolverDnDsyevd(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,WORK_,LWORK_,INFO_)
!$ACC END HOST_DATA
!$ACC EXIT DATA COPYOUT(W(1:N),INFO_) ASYNC(ACC_ASYNC_Q)

    IF (JOBZ=='V'.OR.JOBZ=='v') THEN
!$ACC KERNELS PRESENT(A,Z) ASYNC(ACC_ASYNC_Q)
       DO I=1,N
          Z(1:N,I)=A(1:N,I)
       ENDDO
!$ACC END KERNELS
    ENDIF

!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK_)

    INFO=INFO_*100+IERR
  END SUBROUTINE ACC_DSYEVX


  SUBROUTINE ACC_ZHEEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)
    USE mcusolver
    IMPLICIT NONE
    COMPLEX(q) :: A(LDA,*),Z(LDZ,*),WORK(*)
    REAL(q) :: W(*),RWORK(*),VL,VU,ABSTOL
    INTEGER :: IWORK(*),IFAIL(*),N,LDA,IL,IU,M,LDZ,LWORK,INFO
    CHARACTER(1) :: JOBZ,RANGE,UPLO
    ! local variables
    COMPLEX(q), POINTER :: WORK_(:)
    INTEGER :: JOBZ_,UPLO_,LWORK_,I
    INTEGER(4) :: IERR,INFO_
#ifdef magma_zheevx
    COMPLEX(q) :: WA(N,N),WZ(N,M)
#endif

    IF (.NOT.ACC_EXEC_ON) THEN
       CALL ZHEEVX(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,RWORK,IWORK,IFAIL,INFO)
       RETURN
    ENDIF
#ifdef magma_zheevx
!$ACC HOST_DATA USE_DEVICE(A,Z)
    CALL magmaf_zheevx_gpu( &
   &   JOBZ,RANGE,UPLO,N,c_devloc(A),LDA,VL,VU,IL,IU,ABSTOL,M,W,c_devloc(Z),LDZ,WA,N,WZ,N,WORK,     -1,RWORK,IWORK,IFAIL,INFO)
!$ACC END HOST_DATA
    LWORK_=WORK(1); ALLOCATE(WORK_(LWORK_))
!$ACC HOST_DATA USE_DEVICE(A,Z)
    CALL magmaf_zheevx_gpu( &
   &   JOBZ,RANGE,UPLO,N,c_devloc(A),LDA,VL,VU,IL,IU,ABSTOL,M,W,c_devloc(Z),LDZ,WA,N,WZ,N,WORK_,LWORK_,RWORK,IWORK,IFAIL,INFO)
!$ACC END HOST_DATA
    DEALLOCATE(WORK_)
    RETURN
#endif
    IF (RANGE/='A'.AND.RANGE/='a') THEN
       CALL vtutor%error("ACC_ZHEEVX: ERROR: RANGE= " // RANGE // " not implemented")
    ENDIF

    IERR = cusolverDnSetStream(ACC_CUSOLVERHANDLE,acc_get_cuda_stream(ACC_ASYNC_Q))

    IF (JOBZ=='N'.OR.JOBZ=='n') THEN
       JOBZ_=CUSOLVER_EIG_MODE_NOVECTOR
    ELSE
       JOBZ_=CUSOLVER_EIG_MODE_VECTOR
    ENDIF

    IF (UPLO=='U'.OR.UPLO=='u') THEN
       UPLO_=CUBLAS_FILL_MODE_UPPER
    ELSE
       UPLO_=CUBLAS_FILL_MODE_LOWER
    ENDIF

!$ACC ENTER DATA CREATE(W(1:N)) ASYNC(ACC_ASYNC_Q)
!$ACC HOST_DATA USE_DEVICE(A,W)
    IERR = IERR + cusolverDnZheevd_bufferSize(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,LWORK_)
!$ACC END HOST_DATA

    ALLOCATE(WORK_(LWORK_))
!$ACC ENTER DATA CREATE(INFO_,WORK_) ASYNC(ACC_ASYNC_Q)

!$ACC HOST_DATA USE_DEVICE(A,W,WORK_,INFO_)
    IERR = IERR + cusolverDnZheevd(ACC_CUSOLVERHANDLE,JOBZ_,UPLO_,N,A,LDA,W,WORK_,LWORK_,INFO_)
!$ACC END HOST_DATA
!$ACC EXIT DATA COPYOUT(W(1:N),INFO_) ASYNC(ACC_ASYNC_Q)

    IF (JOBZ=='V'.OR.JOBZ=='v') THEN
!$ACC KERNELS PRESENT(A,Z) ASYNC(ACC_ASYNC_Q)
       DO I=1,N
          Z(1:N,I)=A(1:N,I)
       ENDDO
!$ACC END KERNELS
    ENDIF

!$ACC EXIT DATA DELETE(WORK_) ASYNC(ACC_ASYNC_Q)
!$ACC WAIT(ACC_ASYNC_Q)
    DEALLOCATE(WORK_)

    INFO=INFO_*100+IERR
  END SUBROUTINE ACC_ZHEEVX


!****************** SUBROUTINE SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX *******
!
! This routine takes one complex argument and accepts this as two real
! arguments for the real and imaginary parts, respectively. It then
! performs two separate OpenACC atomic adds, one for each part.
!
!***********************************************************************

  SUBROUTINE SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX(SPLIT_CMPLX,TO_ADD)
!$ACC ROUTINE SEQ
    USE prec
    REAL(q),DIMENSION(2) :: SPLIT_CMPLX
    COMPLEX(q) :: TO_ADD
!$ACC ATOMIC UPDATE
    SPLIT_CMPLX(1)=SPLIT_CMPLX(1)+REAL(TO_ADD) ! real part
!$ACC END ATOMIC
!$ACC ATOMIC UPDATE
    SPLIT_CMPLX(2)=SPLIT_CMPLX(2)+AIMAG(TO_ADD) ! imaginary part
!$ACC END ATOMIC
  END SUBROUTINE SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX


!***********************************************************************
!
! The infrastructure in this module allows subroutines that are executed
! on the GPU to raise critical exceptions to be acted  upon returning
! after returning CPU side.
! It only supports isBug and isError types of exceptions and does not
! support passing information (arguments) from inside the GPU kernel.
! Furthermore only one exception can be logged: when ACC_VTUTOR_LOG
! is called multiple times it will store the instance with the largest
! IdAndSeverity = Id + IDMAX * SEVERITY (see ACC_VTUTOR_LOG).
!
!***********************************************************************

  MODULE tutoracc
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: ACC_VTUTOR_LOG, ACC_VTUTOR_CHECK

  INTEGER :: IdAndSeverity = -1
  INTEGER :: IDMAX = 1000
!$ACC DECLARE COPYIN(IdAndSeverity,IDMAX)

  CONTAINS

  SUBROUTINE ACC_VTUTOR_LOG(SEVERITY,ID)
!$ACC ROUTINE SEQ
    INTEGER :: SEVERITY,ID
!$ACC ATOMIC UPDATE
    IdAndSeverity = MAX(IdAndSeverity, SEVERITY*IDMAX + ID)
  END SUBROUTINE ACC_VTUTOR_LOG


  SUBROUTINE ACC_VTUTOR_CHECK
    USE mopenacc_struct_def, ONLY : ACC_EXEC_ON
    USE tutor, ONLY : vtutor, VMBJ_BrentInitialPositions, VMBJ_TooManyIterations
    USE string, ONLY : str
    ! local variables
    CHARACTER(LEN=:), ALLOCATABLE :: FILENAME
    INTEGER :: ID,SEVERITY

!$ACC WAIT IF(ACC_EXEC_ON)
!$ACC UPDATE SELF(IdAndSeverity) IF(ACC_EXEC_ON)
    IF (IdAndSeverity==-1) RETURN

    SEVERITY = IdAndSeverity/IDMAX
    ID = IdAndSeverity - SEVERITY*IDMAX

    SELECT CASE (ID)
    CASE(VMBJ_BrentInitialPositions)
       FILENAME = 'metagga.F'
    CASE(VMBJ_TooManyIterations)
       FILENAME = 'metagga.F'
    CASE DEFAULT
       FILENAME = 'UnkownFile'
    END SELECT

    CALL vtutor%write(SEVERITY, ID, filename=FILENAME)

    ! currently this point will never be reached because we have
    ! only implemented fatal exceptions ... this might change in
    ! future though
    IdAndSeverity = -1
!$ACC UPDATE DEVICE(IdAndSeverity) IF(ACC_EXEC_ON)
    RETURN
  END SUBROUTINE ACC_VTUTOR_CHECK

  END MODULE tutoracc
#endif

#ifdef CUDA_PINNED_MEMORY
  MODULE mcudafor
  USE prec
  USE iso_c_binding, ONLY : C_PTR,C_SIZE_T,C_LOC

  USE cudafor, ONLY : cudaHostRegister, cudaHostRegisterPortable, &
                      cudaHostUnregister, cudaPointerGetAttributes, cudaPointerAttributes

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: PIN_MEMORY,UNPIN_MEMORY,C_LOC,C_SIZE_T

  CONTAINS

  SUBROUTINE PIN_MEMORY(CPTR,ISIZE,ISTATUS)
    TYPE(C_PTR), VALUE :: CPTR
    INTEGER(KIND=C_SIZE_T) :: ISIZE
    INTEGER, OPTIONAL :: ISTATUS
    ! local variables
    INTEGER :: IERR

    IERR = cudaHostRegister(CPTR,ISIZE,cudaHostRegisterPortable)

    IF (PRESENT(ISTATUS)) ISTATUS=IERR
    RETURN
  END SUBROUTINE PIN_MEMORY


  SUBROUTINE UNPIN_MEMORY(CPTR,ISTATUS)
    TYPE(C_PTR), VALUE :: CPTR
    INTEGER, OPTIONAL :: ISTATUS
    ! local variables
    TYPE(cudaPointerAttributes) :: ATTR
    INTEGER :: IERR

    IERR = cudaPointerGetAttributes(ATTR,CPTR)

    IF (IERR==0) IERR = cudaHostUnregister(CPTR)

    IF (PRESENT(ISTATUS)) ISTATUS=IERR
    RETURN
  END SUBROUTINE UNPIN_MEMORY

  END MODULE mcudafor
#endif
!***********************************************************************
!
! OpenACC related Doxygen documentation page.
!
!***********************************************************************

!> @page openacc OpenACC
!! @tableofcontents
!! @section general All changes
!! All datatypes, subroutines, and functions that have been changed
!! to port VASP to GPU under OpenACC.
!!
!! \li      subrot::add_gamma_from_file
!! \li       hamil::eccp
!! \li       hamil::eccp_tau
!! \li       hamil::eccp_vec
!! \li       david::eddav
!! \li david_inner::eddav_inner
!! \li      subrot::eddiag
!! \li      subrot::eddiag_exact
!! \li    rmm_diis::eddrmm
!! \li       steep::edstep
!! \li            ::elmin
!! \li       hamil::hamilt_local
!! \li       hamil::hamilt_local_tau
!! \li       hamil::hamiltmu
!! \li       hamil::hamiltmu_tau
!! \li       hamil::hamiltmu_vec
!! \li            ::kinhamil
!! \li            ::kinhamil_tau
!! \li            ::kinhamil_vec
!! \li     locproj::lprj_proall
!! \li            ::potlok
!! \li        meta::potlok_metagga
!! \li morbitalmag::set_dd_magatom
!! \li        pawm::set_dd_paw
!! \li        meta::set_kineden
!! \li            ::setdij_
!! \li            ::setdij_avec_
!! \li     msphpro::sphpro_fast
!! \li morbitalmag::vectorpot
!! \li        wave::wvreal
!!
!! @section active Ported
!! \li       hamil::eccp
!! \li       hamil::eccp_tau
!! \li       hamil::eccp_vec
!! \li      subrot::eddiag
!! \li    rmm_diis::eddrmm
!! \li            ::elmin
!! \li       hamil::hamilt_local
!! \li       hamil::hamilt_local_tau
!! \li       hamil::hamiltmu
!! \li       hamil::hamiltmu_tau
!! \li       hamil::hamiltmu_vec
!! \li            ::kinhamil
!! \li            ::kinhamil_tau
!! \li            ::kinhamil_vec
!! \li        wave::wvreal
!!
!! @section inactive OpenACC deactivated
!! \li      subrot::add_gamma_from_file
!! \li       david::eddav
!! \li david_inner::eddav_inner
!! \li      subrot::eddiag_exact
!! \li       steep::edstep
!! \li     locproj::lprj_proall
!! \li         pot::potlok
!! \li        meta::potlok_metagga
!! \li morbitalmag::set_dd_magatom
!! \li        pawm::set_dd_paw
!! \li        meta::set_kineden
!! \li            ::setdij_
!! \li            ::setdij_avec_
!! \li     msphpro::sphpro_fast
!! \li morbitalmag::vectorpot
