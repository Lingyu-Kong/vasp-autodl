#include "symbol.inc"
MODULE nonl_high
  USE prec
  USE nonlr
  USE nonl
  IMPLICIT NONE
!***********************************************************************
!
!  this module contains all high level routines to calculate
!  the non local projection operators
!  either 
!  ) on a plane wave grid
!  ) on a equally spaced grid
!
!***********************************************************************
  CONTAINS
!************************** SUBROUTINE PROALL***************************
!
! this subroutine
! calculates the scalar product of the current wavefunctions
! stored in W with the projectors
!
!  C_lme ion,n = < b_lme ion | psi_n >
!
! and stores the result in W%CPROJ (wave function character)
!
!***********************************************************************

  SUBROUTINE PROALL(GRID,LATT_CUR,NONLR_S,NONL_S,W)
#ifdef _OPENACC
    USE mopenacc
#endif
    USE poscar
    USE lattice
    
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W
  ! local
    INTEGER NK

    PROFILING_START('proall')

#ifdef _OPENACC
!$ACC WAIT IF(ACC_EXEC_ON)
    PUSH_ACC_EXEC_ON(.TRUE.)

    ACC_ASYNC_Q = ACC_ASYNC_ASYNC

    CALL ACC_COPYIN_TYPED_VAR(W)
    CALL ACC_COPYIN_TYPED_VAR(GRID)
    CALL ACC_COPYIN_TYPED_VAR(NONL_S)
    CALL ACC_COPYIN_TYPED_VAR(NONLR_S)
!$ACC ENTER DATA COPYIN(LATT_CUR) __IF_ASYNC__
#endif

    DO NK=1,W%WDES%NKPTS
#ifdef MPI
       IF (MOD(NK-1,W%WDES%COMM_KINTER%NCPU).NE.W%WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
       IF (NONLR_S%LREAL) THEN
          CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,W%WDES)
          CALL RPRO(NONLR_S,W%WDES,W,GRID,NK)
       ELSE
          CALL PHASE(W%WDES,NONL_S,NK)
          CALL PROJ(NONL_S,W%WDES,W,NK)
       ENDIF
    ENDDO

#ifdef _OPENACC
!$ACC WAIT IF(ACC_EXEC_ON)
    CALL ACC_COPYOUT_TYPED_VAR(W)
    CALL ACC_DELETE_TYPED_VAR(GRID)
    CALL ACC_DELETE_TYPED_VAR(NONL_S)
    CALL ACC_DELETE_TYPED_VAR(NONLR_S)
!$ACC EXIT DATA DELETE(LATT_CUR) __IF_ASYNC__

    ACC_ASYNC_Q = ACC_ASYNC_SYNC

    POP_ACC_EXEC_ON
#endif

    PROFILING_STOP('proall')

    RETURN
  END SUBROUTINE PROALL


!************************** SUBROUTINE PROALL***************************
!
! this subroutine calculates the scalar product of the current wavefunctions
! stored in W1(:) with the projectors
!
!  C_lme ion,n = < b_lme ion | psi_n >
!
! and stores the result in W1(:)%CPROJ (wave function character)
!
!***********************************************************************
  

  SUBROUTINE W1_PROJALL(WDES1, W1, NONLR_S, NONL_S, NMAX)
#ifdef _OPENACC
    USE mopenacc_struct_def
#endif
    IMPLICIT NONE
    TYPE (wavedes1) :: WDES1
    TYPE (wavefun1) :: W1(:)
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    INTEGER, OPTIONAL :: NMAX
  ! local
    INTEGER NMAX_, NP
!$  INTEGER NSTRIP,NSTRIP_ACT
!$  INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS

    PROFILING_START('w1_projall')

    IF (PRESENT(NMAX)) THEN
       NMAX_=NMAX
    ELSE
       NMAX_=SIZE(W1)
    ENDIF
 
    IF ( NONLR_S%LREAL ) THEN
#ifdef _OPENACC
       IF (ACC_EXEC_ON) THEN
          CALL RPROMU_ACC(NONLR_S, WDES1, W1, NMAX_, W1%LDO)
          PROFILING_STOP('w1_projall')
          RETURN
       ENDIF
#endif
       IF (NMAX_ >1 ) THEN
          CALL RPROMU(NONLR_S, WDES1, W1, NMAX_, W1%LDO)
       ELSE
          DO NP=1,NMAX_
             IF (.NOT. W1(NP)%LDO) CYCLE
             CALL RPRO1(NONLR_S, WDES1, W1(NP))
          ENDDO
       ENDIF
    ELSE
#ifdef _OPENACC
       IF (ACC_EXEC_ON) THEN
          CALL PROJMU_ACC(NONL_S, WDES1, W1, NMAX_, W1%LDO)
          PROFILING_STOP('w1_projall')
          RETURN
       ENDIF
#endif
       DO NP=1,NMAX_
          IF (.NOT. W1(NP)%LDO) CYCLE
          CALL PROJ1(NONL_S,WDES1,W1(NP))
       ENDDO
    ENDIF

    PROFILING_STOP('w1_projall')

  END SUBROUTINE W1_PROJALL


  SUBROUTINE W1_PROJ(W1, NONLR_S, NONL_S)
#ifdef _OPENACC
    USE mopenacc_struct_def
#endif
    IMPLICIT NONE
    TYPE (wavefun1) ::  W1
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct)  NONL_S
  ! local
    IF (NONLR_S%LREAL) THEN
#ifdef _OPENACC
       IF (ACC_EXEC_ON) THEN
          CALL RPRO1_ACC(NONLR_S,W1%WDES1,W1)
       ELSE
#endif
          CALL RPRO1(NONLR_S,W1%WDES1,W1)
#ifdef _OPENACC
       ENDIF
#endif
    ELSE
       CALL PROJ1(NONL_S,W1%WDES1,W1)
    ENDIF
  END SUBROUTINE W1_PROJ

END MODULE nonl_high
