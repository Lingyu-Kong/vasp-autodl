!#define nonlr_single
#define avoid_reduction
!#define IntelKNL
#include "symbol.inc"
!***********************************************************************
!
!  this module contains all the routines required to support 
!  real space presentation of the non local projection operators
!  on an equally spaced grid
!
!***********************************************************************
MODULE nonlr
  USE prec
  USE wave
  USE mgrid
  USE nonlr_struct_def
#ifdef _OPENACC
  USE mopenacc
#endif

  INTERFACE
     SUBROUTINE RACC0(NONLR_S, WDES1, CPROJ_LOC, CRACC)
       USE nonlr_struct_def
       USE wave
       TYPE (nonlr_struct) NONLR_S
       TYPE (wavedes1)     WDES1
       COMPLEX(q) CRACC
       GDEF       CPROJ_LOC
     END SUBROUTINE RACC0
  END INTERFACE

  INTERFACE
     SUBROUTINE RACC0_HF(NONLR_S, WDES1, CPROJ_LOC, CRACC)
       USE nonlr_struct_def
       USE wave
       TYPE (nonlr_struct) NONLR_S
       TYPE (wavedes1)     WDES1
       GDEF   CRACC
       GDEF   CPROJ_LOC
     END SUBROUTINE RACC0_HF
  END INTERFACE


  INTERFACE
     SUBROUTINE RACC0_REAL(NONLR_S, WDES1, CPROJ_LOC, CRACC)
       USE nonlr_struct_def
       USE wave
       TYPE (nonlr_struct) NONLR_S
       TYPE (wavedes1)     WDES1
       REAL(q) CRACC
       GDEF   CPROJ_LOC
     END SUBROUTINE RACC0_REAL
  END INTERFACE


  INTERFACE
     SUBROUTINE RACC0MU(NONLR_S, WDES1, CPROJ_LOC, CRACC, LD, NSIM, LDO)
       USE nonlr_struct_def
       USE wave
       TYPE (nonlr_struct) NONLR_S
       TYPE (wavedes1)     WDES1
       COMPLEX(q) CRACC
       INTEGER    LD
       GDEF       CPROJ_LOC
       INTEGER    NSIM
       LOGICAL LDO(*)
     END SUBROUTINE RACC0MU
  END INTERFACE

  INTERFACE
     SUBROUTINE RACC0MU_HF(NONLR_S, WDES1, CPROJ_LOC, LD1, CRACC, LD2, NSIM)
       USE nonlr_struct_def
       USE wave
       TYPE (nonlr_struct) NONLR_S
       TYPE (wavedes1)     WDES1
       GDEF       CRACC
       INTEGER    LD1,LD2
       GDEF       CPROJ_LOC
       INTEGER    NSIM
     END SUBROUTINE RACC0MU_HF
  END INTERFACE

CONTAINS


!****************** SUBROUTINE NONLR_SETUP ****************************
!
! NONLR_SETUP is the base initialisation routine
! o it sets the number of types and the number of ions
! o it links the positions descriptors to the T_INFO structure
! o it links the tables to the pseudpotential structure
!
! before the data structure can be used in actual calculations
! the following calls must be made:
! ) REAL_OPTLAY  determine the number of grid point in the
!                real space cutoff spheres
! ) NONLR_ALLOC  allocate the required tables 
! ) SPHERE       set the tables for the fast evaluation of the
!                non local projetors
! 
!***********************************************************************

  SUBROUTINE  NONLR_SETUP(NONLR_S,T_INFO,P,LREAL,LSPIRAL,COMM)
    USE pseudo
    USE poscar
    USE mpimy
    IMPLICIT NONE


    TYPE (nonlr_struct) NONLR_S
    TYPE (type_info)   T_INFO
    TYPE (potcar)      P(T_INFO%NTYP)
    LOGICAL LREAL
    LOGICAL LSPIRAL
    TYPE (communic), TARGET, OPTIONAL :: COMM
    INTEGER, EXTERNAL :: MAXL1
    ! local var
    INTEGER NT


    NONLR_S%LREAL  =LREAL
    NONLR_S%NK     =0
    NONLR_S%NTYP   =T_INFO%NTYP
    NONLR_S%NIONS  =T_INFO%NIONS
    NONLR_S%IRALLOC=-1
    NONLR_S%IRMAX  =-1
    NONLR_S%NITYP  =>T_INFO%NITYP
    NONLR_S%ITYP   =>T_INFO%ITYP
    NONLR_S%POSION =>T_INFO%POSION
    NONLR_S%LSPIRAL=LSPIRAL

    ALLOCATE(NONLR_S%LMAX  (NONLR_S%NTYP), &
         NONLR_S%LMMAX (NONLR_S%NTYP), &
         NONLR_S%CHANNELS(NONLR_S%NTYP), &
         NONLR_S%PSRMAX(NONLR_S%NTYP), &
         NONLR_S%RSMOOTH(NONLR_S%NTYP), &
         NONLR_S%BETA  (NONLR_S%NTYP))

    NULLIFY(NONLR_S%NLIMAX, NONLR_S%NLI, NONLR_S%RPROJ, NONLR_S%CRREXP, NONLR_S%VKPT_SHIFT)

    DO NT=1,T_INFO%NTYP
       NONLR_S%LMAX(NT)    = MAXL1(P(NT))
       NONLR_S%LMMAX(NT)   =P(NT)%LMMAX
       NONLR_S%CHANNELS(NT)=P(NT)%LMAX
       NONLR_S%PSRMAX(NT)  =P(NT)%PSRMAX
       NONLR_S%RSMOOTH(NT) =0
       NONLR_S%BETA(NT)%PSPRNL=>P(NT)%PSPRNL
       NONLR_S%BETA(NT)%LPS   =>P(NT)%LPS
    ENDDO

    NULLIFY(NONLR_S%LMBASE,NONLR_S%NLIBASE)

    NONLR_S%SELECTED_ION=-1

#ifdef shmem_rproj
    NONLR_S%COMM_shmem => COMM
#endif
    RETURN
  END SUBROUTINE NONLR_SETUP


!****************** SUBROUTINE NONLR_ALLOC *****************************
!
! allocate required work arrays
! this function can be called only after NONLR_SETUP
! since it requires the number of data points in the real space cutoff
! spheres (NONLR_S%IRMAX and NONLR_S%IRALLOC)
!
!> @details @ref openmp :
!> Allocate
!>      nonlr_struct::nlimax,
!>      nonlr_struct::nli,
!>      nonlr_struct::rproj,
!>      nonlr_struct::crrexp,
!>      nonlr_struct::nlibase,
!> with the additional dimension openmp::omp_nthreads_nonlr_rspace
!> (slowest index).
!
!***********************************************************************

  SUBROUTINE NONLR_ALLOC(NONLR_S)
    USE pseudo
    USE ini
    USE mpi_shmem
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT NONE

    TYPE (nonlr_struct), TARGET :: NONLR_S
    ! local variables
    INTEGER NIONS,IRMAX
!$  INTEGER, EXTERNAL :: OMP_GET_MAX_THREADS

    NIONS = NONLR_S%NIONS
    IRMAX = NONLR_S%IRMAX

    IF (NONLR_S%LREAL) THEN
       IF (NONLR_S%IRMAX<0.OR.NONLR_S%IRALLOC<0) THEN
          CALL vtutor%bug("internal ERROR in NONLR_ALLOC: IRMAX or IRALLOC is not set " // &
             str(NONLR_S%IRMAX) // " " // str(NONLR_S%IRALLOC) // "\n call REAL_OPTLAY before NONLR_ALLOC", &
             __FILE__, __LINE__)
       ENDIF
       ALLOCATE(NONLR_S%NLIMAX(NIONS  __NOACC_omp_arg(omp_nthreads_nonlr_rspace)), &
          NONLR_S%NLI   (IRMAX,NIONS  __NOACC_omp_arg(omp_nthreads_nonlr_rspace)))
#ifdef shmem_rproj
       CALL M_shmem_alloc(NONLR_S%COMM_shmem,NONLR_S%IRALLOC  __NOACC_omp_arg(omp_nthreads_nonlr_rspace),NONLR_S%RPROJ_shmem,NONLR_S%RPROJ)
#else
       ALLOCATE(NONLR_S%RPROJ (NONLR_S%IRALLOC  __NOACC_omp_arg(omp_nthreads_nonlr_rspace)))
#endif
#ifdef shmem_rproj
       CALL REGISTER_ALLOCATE(8._q*SIZE(NONLR_S%RPROJ)/NONLR_S%COMM_shmem%NCPU, "nonlr-proj")
#else
       CALL REGISTER_ALLOCATE(8._q*SIZE(NONLR_S%RPROJ), "nonlr-proj")
#endif

#ifndef gammareal
       IF (.NOT.NONLR_S%LSPIRAL ) THEN
          ALLOCATE(NONLR_S%CRREXP(IRMAX,NIONS,1 __NOACC_omp_arg(omp_nthreads_nonlr_rspace)))
          CALL REGISTER_ALLOCATE(16._q*SIZE(NONLR_S%CRREXP), "nonlr-proj")
       ELSE
          ALLOCATE(NONLR_S%CRREXP(IRMAX,NIONS,2 __NOACC_omp_arg(omp_nthreads_nonlr_rspace)))
          CALL REGISTER_ALLOCATE(16._q*SIZE(NONLR_S%CRREXP), "nonlr-proj")
       ENDIF
#endif
       ALLOCATE(NONLR_S%LMBASE(NIONS+1),NONLR_S%NLIBASE(NIONS+1 __NOACC_omp_arg(omp_nthreads_nonlr_rspace)))
    ENDIF
    RETURN
  END SUBROUTINE NONLR_ALLOC

!****************** SUBROUTINE NONLR_DEALLOC ***************************
!
! deallocate all work arrays
! but leave the links to the pseudopotentials and ions open
! such that the projectors can be reinitialized by a call to 
! NONLR_SETUP, NONLR_ALLOC and SPHERE
!
!***********************************************************************

  SUBROUTINE  NONLR_DEALLOC(NONLR_S)
    USE ini
    USE mpi_shmem
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT NONE

    TYPE (nonlr_struct), TARGET :: NONLR_S

    IF (NONLR_S%LREAL) THEN
#ifdef shmem_rproj
       CALL DEREGISTER_ALLOCATE(8._q*SIZE(NONLR_S%RPROJ)/NONLR_S%COMM_shmem%NCPU, "nonlr-proj")
#else
       CALL DEREGISTER_ALLOCATE(8._q*SIZE(NONLR_S%RPROJ), "nonlr-proj")
#endif
       DEALLOCATE(NONLR_S%NLIMAX,NONLR_S%NLI)
#ifdef shmem_rproj
       CALL M_shmem_dealloc(NONLR_S%RPROJ_shmem)
#else
       DEALLOCATE(NONLR_S%RPROJ)
#endif
       IF (ASSOCIATED(NONLR_S%CRREXP)) THEN
          CALL DEREGISTER_ALLOCATE(16._q*SIZE(NONLR_S%CRREXP), "nonlr-proj")
          DEALLOCATE(NONLR_S%CRREXP)
       ENDIF
    ENDIF

    NULLIFY(NONLR_S%NLIMAX, NONLR_S%NLI, NONLR_S%RPROJ,NONLR_S%CRREXP)
    RETURN
  END SUBROUTINE NONLR_DEALLOC

!****************** SUBROUTINE NONLR_DESTROY ***************************
!
! destroy the  links to the pseudopotentials and release all resources
! allocated by the NONLR_S descriptor
! the subroutine performs all the operations performed by 
! NONLR_DEALLOC
! and destroys all other links as well
!
!***********************************************************************

  SUBROUTINE NONLR_DESTROY(NONLR_S)
    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    ! local variables

    IF (NONLR_S%LREAL) THEN
       IF (ASSOCIATED(NONLR_S%NLI))    DEALLOCATE(NONLR_S%NLI)
       IF (ASSOCIATED(NONLR_S%CRREXP)) DEALLOCATE(NONLR_S%CRREXP)
#ifdef shmem_rproj
       IF (ASSOCIATED(NONLR_S%RPROJ))  CALL M_shmem_dealloc(NONLR_S%RPROJ_shmem)
#else
       IF (ASSOCIATED(NONLR_S%RPROJ))  DEALLOCATE(NONLR_S%RPROJ)
#endif
       IF (ASSOCIATED(NONLR_S%NLIMAX)) DEALLOCATE(NONLR_S%NLIMAX)
    ENDIF
    NULLIFY(NONLR_S%NLIMAX, NONLR_S%NLI, NONLR_S%RPROJ,NONLR_S%CRREXP)

    DEALLOCATE(NONLR_S%LMAX, &
         NONLR_S%LMMAX, &
         NONLR_S%CHANNELS, &
         NONLR_S%PSRMAX, &
         NONLR_S%RSMOOTH, &
         NONLR_S%BETA)

    DEALLOCATE(NONLR_S%LMBASE, &
         NONLR_S%NLIBASE)

    RETURN
  END SUBROUTINE NONLR_DESTROY

!****************** SUBROUTINE NONLR_ASSIGN ****************************
!
! This routine accomplishes the same thing as the assignment:
!
!   NONLR_LHS = NONLR_RHS
!
! would, but allows for the use of "fast"-memory (in the sense that
! the fm_* arrays are excluded from the assignment).
!
!***********************************************************************

  SUBROUTINE NONLR_ASSIGN(NONLR_LHS,NONLR_RHS)
    IMPLICIT NONE
    TYPE (nonlr_struct) NONLR_LHS,NONLR_RHS

    NONLR_LHS%LREAL       = NONLR_RHS%LREAL
    NONLR_LHS%NTYP        = NONLR_RHS%NTYP
    NONLR_LHS%NIONS       = NONLR_RHS%NIONS
    NONLR_LHS%SELECTED_ION= NONLR_RHS%SELECTED_ION
    NONLR_LHS%IRMAX       = NONLR_RHS%IRMAX
    NONLR_LHS%IRALLOC     = NONLR_RHS%IRALLOC
    NONLR_LHS%NK          = NONLR_RHS%NK
    NONLR_LHS%LSPIRAL     = NONLR_RHS%LSPIRAL

    NONLR_LHS%NITYP       =>NONLR_RHS%NITYP
    NONLR_LHS%ITYP        =>NONLR_RHS%ITYP
    NONLR_LHS%LMAX        =>NONLR_RHS%LMAX
    NONLR_LHS%LMMAX       =>NONLR_RHS%LMMAX
    NONLR_LHS%CHANNELS    =>NONLR_RHS%CHANNELS
    NONLR_LHS%PSRMAX      =>NONLR_RHS%PSRMAX
    NONLR_LHS%RSMOOTH     =>NONLR_RHS%RSMOOTH
    NONLR_LHS%POSION      =>NONLR_RHS%POSION
    NONLR_LHS%VKPT_SHIFT  =>NONLR_RHS%VKPT_SHIFT
    NONLR_LHS%BETA        =>NONLR_RHS%BETA

    NONLR_LHS%NLIMAX      =>NONLR_RHS%NLIMAX

    NONLR_LHS%LMBASE      =>NONLR_RHS%LMBASE
    NONLR_LHS%NLIBASE     =>NONLR_RHS%NLIBASE

    NONLR_LHS%NLI         =>NONLR_RHS%NLI
    NONLR_LHS%RPROJ       =>NONLR_RHS%RPROJ
    NONLR_LHS%CRREXP      =>NONLR_RHS%CRREXP

#ifdef shmem_rproj
    NONLR_LHS%COMM_shmem  =>NONLR_RHS%COMM_shmem
#endif
    RETURN
  END SUBROUTINE NONLR_ASSIGN

!****************** SUBROUTINE REAL_OPTLAY *****************************
!
! determine 
! NONLR_S%IRMAX and NONLR_S%IRALLOC in the non local PP structure
!
! for the parallel version the subroutine also optimizes
! the layout (i.e. data distribution) of the data points (columns)
! in real space such that the total number of grid points
! is the same on all nodes (this requires an update of the GRID structure)
! if LNOREDIS is set the data layout is, however,
! not allowed to change in the GRID structure
!
! NOTE: if the data layout is updated in the GRID structure the parallel
! fast Fourier transformation must be reinitialised
! furthermore  the routine signals to the calling routine whether 
! a reallocation of the work arrays in NONLR_S is required
! via NONLR_ALLOC. This is done be checking the current setting
! of NONLR_S%IRMAX and NONLR_S%IRALLOC
!
!> @details @ref openmp : Partition each real space projection operator
!> into openmp::omp_nthreads_nonlr_rspace parts. The work on these
!> parts will be distributed over the available OpenMP threads lateron.
!> @n
!> Currently OMP_NUM_THREADS/=1 requires NCORE=1! This means
!> the fast index of the real space grid of the non local projection
!> operators is the x-index (GRID\%RL\%NFAST/=3).
!> @n
!> For openmp::omp_nonlr_planewise = .TRUE. (default) the projectors are
!> distributed over (x,y)-planes in a round robin fashion.
!> In case openmp::omp_nonlr_planewise = .FALSE. the projectors are
!> distributed over (z,y)-columns.
!
!***********************************************************************

  SUBROUTINE REAL_OPTLAY(GRID,LATT_CUR,NONLR_S,LNOREDIS, &
       LREALLOCATE,IU6,IU0)

    USE lattice
    USE constant
    USE pseudo
    USE openmp, ONLY : omp_nthreads_nonlr_rspace,omp_nonlr_planewise

    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (grid_3d)      GRID
    TYPE (latt)         LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    ! local work arrays
    INTEGER, ALLOCATABLE :: USED_ROWS(:,:) ! counts how many elements
    ! must be allocated for one row
    INTEGER, ALLOCATABLE :: REDISTRIBUTION_INDEX(:)
    LOGICAL  LNOREDIS   ! no redistribution allowed
    LOGICAL  LREALLOCATE

!$  INTEGER idim

    LREALLOCATE=.FALSE.

    IF (.NOT. NONLR_S%LREAL) RETURN
!=======================================================================
! loop over all ions
!=======================================================================
    NLIIND=0
    IRMAX =0
    NIS=1

    IF (GRID%RL%NFAST==3) THEN
       ALLOCATE(USED_ROWS(GRID%NGX,GRID%NGY), &
            REDISTRIBUTION_INDEX(GRID%NGX*GRID%NGY) )
       USED_ROWS=0
       IRALLOC=0  ! number of real space proj on local node
    ENDIF

!$OMP PARALLEL DO DEFAULT(PRIVATE) FIRSTPRIVATE(NIS) &
!$OMP SHARED(omp_nthreads_nonlr_rspace,NONLR_S,GRID,LATT_CUR,USED_ROWS,omp_nonlr_planewise) &
!$OMP REDUCTION(max:IRMAX,NLIIND,IRALLOC)
!$  omp: DO idim=1,omp_nthreads_nonlr_rspace

!$  IRMAX=MAX(IRMAX,0); NLIIND=MAX(NLIIND,0); IRALLOC=MAX(IRALLOC,0)

    type: DO NT=1,NONLR_S%NTYP
       IF (NONLR_S%LMMAX(NT)==0) GOTO 600
       LMMAXC=NONLR_S%LMMAX(NT)
       ions: DO NI=NIS,NONLR_S%NITYP(NT)+NIS-1

          IF (NONLR_S%SELECTED_ION>0 .AND. NONLR_S%SELECTED_ION/=NI) CYCLE
          ARGSC=NPSRNL/NONLR_S%PSRMAX(NT)

!=======================================================================
! find lattice points contained within the cutoff-sphere
! this loop might be done in scalar unit
!=======================================================================
          F1=1._q/GRID%NGX
          F2=1._q/GRID%NGY
          F3=1._q/GRID%NGZ

!-----------------------------------------------------------------------
! restrict loop to points contained within a cubus around the ion
!-----------------------------------------------------------------------
!sh
!          D1= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(1)*GRID%NGX
!          D2= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(2)*GRID%NGY
!          D3= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(3)*GRID%NGZ
          D1= (NONLR_S%PSRMAX(NT)+NONLR_S%RSMOOTH(NT))*LATT_CUR%BNORM(1)*GRID%NGX
          D2= (NONLR_S%PSRMAX(NT)+NONLR_S%RSMOOTH(NT))*LATT_CUR%BNORM(2)*GRID%NGY
          D3= (NONLR_S%PSRMAX(NT)+NONLR_S%RSMOOTH(NT))*LATT_CUR%BNORM(3)*GRID%NGZ

          N3LOW= INT(NONLR_S%POSION(3,NI)*GRID%NGZ-D3+10*GRID%NGZ+.99_q)-10*GRID%NGZ
          N2LOW= INT(NONLR_S%POSION(2,NI)*GRID%NGY-D2+10*GRID%NGY+.99_q)-10*GRID%NGY
          N1LOW= INT(NONLR_S%POSION(1,NI)*GRID%NGX-D1+10*GRID%NGX+.99_q)-10*GRID%NGX

          N3HI = INT(NONLR_S%POSION(3,NI)*GRID%NGZ+D3+10*GRID%NGZ)-10*GRID%NGZ
          N2HI = INT(NONLR_S%POSION(2,NI)*GRID%NGY+D2+10*GRID%NGY)-10*GRID%NGY
          N1HI = INT(NONLR_S%POSION(1,NI)*GRID%NGX+D1+10*GRID%NGX)-10*GRID%NGX
!-----------------------------------------------------------------------
! loop over cubus
! MPI version z ist the fast index
!-----------------------------------------------------------------------
          IF (GRID%RL%NFAST==3) THEN
          IND=1

          DO N2=N2LOW,N2HI
             X2=(N2*F2-NONLR_S%POSION(2,NI))
             N2P=MOD(N2+10*GRID%NGY,GRID%NGY)

!!$           IF (MOD(N2P,omp_nthreads_nonlr_rspace)/=(idim-1)) CYCLE

             DO N1=N1LOW,N1HI
                X1=(N1*F1-NONLR_S%POSION(1,NI))
                N1P=MOD(N1+10*GRID%NGX,GRID%NGX)

                NCOL=GRID%RL%INDEX(N1P,N2P)

                DO N3=N3LOW,N3HI
                   X3=(N3*F3-NONLR_S%POSION(3,NI))

                   X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
                   Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
                   Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)

                   D=SQRT(X*X+Y*Y+Z*Z)
                   ARG=(D*ARGSC)+1
                   NADDR=INT(ARG)

!sh                IF (NADDR<NPSRNL) THEN
                   IF (D<(NPSRNL-1)/ARGSC+NONLR_S%RSMOOTH(NT)) THEN

                      IND=IND+1
                      USED_ROWS(N1P+1,N2P+1)=USED_ROWS(N1P+1,N2P+1)+LMMAXC

!$                    IF (MOD(NCOL,omp_nthreads_nonlr_rspace)/=(idim-1)) CYCLE

                      ! if on local processor add to IRALLOC
                      IF (NCOL/=0) IRALLOC=IRALLOC+LMMAXC
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
          ELSE
!-----------------------------------------------------------------------
! loop over cubus around one ion
! conventional version x is fast index
!-----------------------------------------------------------------------
          IND=1
          DO N3=N3LOW,N3HI
             X3=(N3*F3-NONLR_S%POSION(3,NI))
             N3P=MOD(N3+10*GRID%NGZ,GRID%NGZ)

!$           IF (MOD(N3P,omp_nthreads_nonlr_rspace)/=(idim-1).AND.omp_nonlr_planewise) CYCLE

             DO N2=N2LOW,N2HI
                X2=(N2*F2-NONLR_S%POSION(2,NI))
                N2P=MOD(N2+10*GRID%NGY,GRID%NGY)

                NCOL=GRID%RL%INDEX(N2P,N3P)

!$              IF (MOD(NCOL,omp_nthreads_nonlr_rspace)/=(idim-1).AND.(.NOT.omp_nonlr_planewise)) CYCLE

                DO N1=N1LOW,N1HI
                   X1=(N1*F1-NONLR_S%POSION(1,NI))

                   X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
                   Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
                   Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)

                   D=SQRT(X*X+Y*Y+Z*Z)
                   ARG=(D*ARGSC)+1
                   NADDR=INT(ARG)

!sh                IF (NADDR<NPSRNL) THEN
                   IF (D<(NPSRNL-1)/ARGSC+NONLR_S%RSMOOTH(NT)) THEN
                      N1P=MOD(N1+10*GRID%NGX,GRID%NGX)
                      NCHECK=N1P+(NCOL-1)*GRID%NGX+1
                      IF (NCHECK /= 1+N1P+GRID%NGX*(N2P+GRID%NGY* N3P)) THEN
                         CALL vtutor%bug("REAL_OPT: internal ERROR: " // str(N1P) // " " // str(N2P) &
                            // " " // str(N3P) // " " // str(NCOL), __FILE__, __LINE__)
                      ENDIF
                      IND=IND+1
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
          ENDIF

          INDMAX=IND-1
          IRMAX =MAX(IRMAX,INDMAX)
          NLIIND=NONLR_S%LMMAX(NT)*INDMAX+NLIIND
       ENDDO ions
600    NIS = NIS+NONLR_S%NITYP(NT)
    ENDDO type

!$  ENDDO omp
!$OMP END PARALLEL DO

!=======================================================================
! now redistribute rows in MPI version
!=======================================================================
#ifdef MPI
    IF (GRID%RL%NFAST==3) THEN
       ! first check whether everything is ok
       NLISUM=SUM(USED_ROWS)
       IF (NLISUM /= NLIIND) THEN
          CALL vtutor%bug("REAL_OPTLAY: internal error (1) " // str(NLISUM) // " " // str(NLIIND), __FILE__, __LINE__)
       ENDIF
       ! setup redistribution index
       NCOL_TOT=GRID%NGY*GRID%NGX
       DO I=1,NCOL_TOT
          REDISTRIBUTION_INDEX(I)=I
       ENDDO
       IF (.NOT. LNOREDIS .AND. GRID%COMM%NCPU>1 ) THEN
          WRITE(*,*) 'resort distribution'
          CALL SORT_REDIS(NCOL_TOT,USED_ROWS(1,1),REDISTRIBUTION_INDEX(1))
          CALL REAL_OPTLAY_GRID(GRID,     REDISTRIBUTION_INDEX,USED_ROWS(1,1),IRALLOC)
       ENDIF
       IF (IRMAX >  NONLR_S%IRMAX) THEN
          ! IRMAX  is the maximum global number, could be improved !!!!
          NONLR_S%IRMAX   =IRMAX  *1.1
#ifdef memalign
          NONLR_S%IRMAX=__align_DIM(NONLR_S%IRMAX)
#endif
          LREALLOCATE=.TRUE.
       ENDIF
       IF( IRALLOC > NONLR_S%IRALLOC) THEN
          ! more safety on parallel machines increase by 20 %
          NONLR_S%IRALLOC =IRALLOC*1.2
          LREALLOCATE=.TRUE.
       ENDIF

       IALLOC_MAX=  IRALLOC
       IALLOC_MIN= -IRALLOC
       CALLMPI (M_max_i( GRID%COMM, IALLOC_MAX, 1))
       CALLMPI (M_max_i( GRID%COMM, IALLOC_MIN, 1))
       CALLMPI (M_sum_i( GRID%COMM, IRALLOC, 1))
       IALLOC_MIN=-IALLOC_MIN

       AKBYTES=1024/8  ! conversion from 8 byte words to kbytes
!$     AKBYTES=AKBYTES/omp_nthreads_nonlr_rspace

       IF (.NOT. LNOREDIS .AND. GRID%COMM%NCPU>1  .AND. IU6>=0) &
            WRITE(IU6,*)'redistribution in real space done'
       IF (.NOT. LNOREDIS .AND. GRID%COMM%NCPU>1 .AND. IU0>=0) &
            WRITE(IU6,*)'redistribution in real space done'
       IF (IU6>=0) &
            WRITE(IU6,1) IRALLOC/AKBYTES,IALLOC_MAX/AKBYTES,IALLOC_MIN/AKBYTES
1      FORMAT(/' real space projection operators:'/ &
            '  total allocation   :',F14.2,' KBytes'/ &
            '  max/ min on nodes  :',2F14.2/)

       IF (NLISUM /= NLIIND) THEN
          CALL vtutor%bug("REAL_OPTLAY: internal error (2) " // str(IRALLOC) // " " // str(NLIIND), __FILE__, __LINE__)
       ENDIF
       DEALLOCATE(USED_ROWS,REDISTRIBUTION_INDEX)
    ELSE
#endif
!=======================================================================
! serial version
!=======================================================================
! to avoid too often reallocation increase values by 10 %
       IF (IRMAX >  NONLR_S%IRMAX) THEN
          NONLR_S%IRMAX   =IRMAX*1.1
#ifdef memalign
          NONLR_S%IRMAX=__align_DIM(NONLR_S%IRMAX)
#endif
          LREALLOCATE=.TRUE.
       ENDIF

       IF( NLIIND > NONLR_S%IRALLOC) THEN
          NONLR_S%IRALLOC =NLIIND*1.1
          LREALLOCATE=.TRUE.
       ENDIF
#ifdef MPI
    ENDIF
#endif
  END SUBROUTINE REAL_OPTLAY


#ifdef MPI
!
! step through all columns and distribute them onto proc.
! in the manner 1 ... NCPU - NCPU ... 1 - 1 ... NCPU - etc.
!
  SUBROUTINE REAL_OPTLAY_GRID(GRID,REDISTRIBUTION_INDEX,USED_ROWS,IRALLOC)
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT REAL(q) (A-H,O-Z)

    TYPE (grid_3d)     GRID
    INTEGER REDISTRIBUTION_INDEX(GRID%NGX*GRID%NGY), &
         USED_ROWS(GRID%NGX*GRID%NGY)
    LOGICAL LUP

!$  INTEGER idim

    GRID%RL%INDEX= 0
    GRID%RL%NCOL = 0

    NODE_TARGET=0  ! NODE onto which column has to go
    LUP=.TRUE.     ! determines whether NODE_TARGET is increased or decreased
    IRALLOC=0

    NCOL_TOT=GRID%NGY*GRID%NGX

!$OMP PARALLEL DO DEFAULT(PRIVATE) &
!$OMP SHARED(omp_nthreads_nonlr_rspace,NCOL_TOT,REDISTRIBUTION_INDEX,GRID) &
!$OMP REDUCTION(max:IRALLOC)
!$  omp: DO idim=1,omp_nthreads_nonlr_rspace

!$  IRALLOC=MAX(IRALLOC,0)

    DO NCOL=1,NCOL_TOT
       IND_REDIS=REDISTRIBUTION_INDEX(NCOL)
       N2=MOD(IND_REDIS-1,GRID%NGX)+1         ! x index (is fast)
       N3=   (IND_REDIS-1)/GRID%NGX+1         ! y index
       IF (LUP) THEN
          IF (NODE_TARGET == GRID%COMM%NCPU) THEN
             LUP=.FALSE.
          ELSE
             NODE_TARGET=NODE_TARGET+1
          ENDIF
       ELSE
          IF (NODE_TARGET == 1) THEN
             LUP=.TRUE.
          ELSE
             NODE_TARGET=NODE_TARGET-1
          ENDIF
       ENDIF

       IND_ON_CPU=(NCOL-1)/GRID%COMM%NCPU+1

       ! element on local node
       ! set up required elements
       IF (NODE_TARGET == GRID%COMM%NODE_ME) THEN
          GRID%RL%NCOL=GRID%RL%NCOL+1
          IF (IND_ON_CPU /= GRID%RL%NCOL) THEN
             CALL vtutor%bug("REAL_OPTLAY: internal error(3) " // str(GRID%COMM%NODE_ME) // " " // &
                str(IND_ON_CPU) // " " // str(GRID%RL%NCOL), __FILE__, __LINE__)
          ENDIF

!$        IF (MOD(IND_ON_CPU,omp_nthreads_nonlr_rspace)/=(idim-1)) CYCLE

          GRID%RL%INDEX(N2-1,N3-1)=IND_ON_CPU
          GRID%RL%I2(IND_ON_CPU)=N2 ! I2 contains x index
          GRID%RL%I3(IND_ON_CPU)=N3 ! I3      the y index
          IRALLOC=IRALLOC+USED_ROWS(NCOL)
       ENDIF
    ENDDO

!$  ENDDO omp
!$OMP END PARALLEL DO

  END SUBROUTINE REAL_OPTLAY_GRID
#endif


!****************** SUBROUTINE RSPHER  *********************************
!
!  subroutine RSPHER calculates the spherical harmonics multiplied
!  by the radial projector function in real space
!  the result is the real space projection operator NONLR_S
!    RPROJ = 1/Omega ^(1/2) Xi(r-R(N)) Y_lm(r-R(N) exp(i k r-R(N))
!
!  all ions can be displaced by a constant shift to allow
!  the evaluation of the first derivative of the projector functions
!
!  RSPHER is the simple interface
!  whereas RSPHER_ALL allows for finite difference calculations
!
!  IDISPL:
!     0 set projector function
!     1 use finite differences to calculate the derivative of
!       the projector function with respect to the specified displacement
!  LOMEGA: use 1/volume scaling (required for FAST_AUG)
!
!***********************************************************************

  SUBROUTINE RSPHER(GRID,NONLR_S, LATT_CUR )
    USE mpimy
    USE lattice
    USE constant

    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    TYPE (wavedes)     WDES
    INTEGER NK
    REAL(q)   DISPL(3,NONLR_S%NIONS)

    DISPL=0
    CALL RSPHER_ALL(GRID,NONLR_S, LATT_CUR, LATT_CUR, LATT_CUR, &
         DISPL,DISPL, 0)
    RETURN
  END SUBROUTINE RSPHER


!****************** SUBROUTINE RSPHER_ALL ******************************
!
!  subroutine RSPHER_ALL calculates the spherical harmonics multiplied
!  by the radial projector function in real space
!  the result is the real space projection operator NONLR_S
!    RPROJ = 1/Omega ^(1/2) Xi(r-R(N)) Y_lm(r-R(N) exp(i k r-R(N))
!
!  all ions can be displaced by a constant shift to allow
!  the evaluation of the first derivative of the projector functions
!
!  IDISPL:
!     0 set projector function
!     1 use finite differences to calculate the derivative of
!       the projector function with respect to the specified displacement
!  LOMEGA: use 1/volume scaling (required for FAST_AUG)
!
!> @details @ref openmp : Partition each real space projection operator
!> into openmp::omp_nthreads_nonlr_rspace parts. The work on these
!> parts will be distributed over the available OpenMP threads lateron.
!> @n
!> Currently OMP_NUM_THREADS/=1 requires NCORE=1. Consequently
!> the fast index of the real space grid of the non local projection
!> operators is the x-index (GRID\%RL\%NFAST/=3).
!> @n
!> For openmp::omp_nonlr_planewise = .TRUE. (default) the projectors are
!> distributed over (x,y)-planes in a round robin fashion.
!> In case openmp::omp_nonlr_planewise = .FALSE. the projectors are
!> distributed over (z,y)-columns.
!> @n@n
!> In addition to NONLR_S\%RPROJ (nonlr_struct::rproj), this routine
!> sets NONLR_S\%LMBASE (nonlr_struct::lmbase) and NONLR_S\%NLIBASE
!> (nonlr_struct::nlibase).
!> The latter two arrays allows us to contract nested loops over
!> \"types\" + \"ions-of-typ\" to single loops over \"all ions\".
!> In the OpenMP version this is used in many routines
!> (see for instance nonlr::rpromu_hf).
!
!***********************************************************************

  SUBROUTINE RSPHER_ALL(GRID,NONLR_S,LATT_FIN1, LATT_FIN2, LATT_CUR, &
       DISPL1, DISPL2, IDISPL, LOMEGA)
    USE pseudo
    USE mpimy
    USE lattice
    USE constant
    USE asa
    USE openmp, ONLY : omp_nthreads_nonlr_rspace,omp_nonlr_planewise
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (nonlr_struct) NONLR_S
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR,LATT_FIN1,LATT_FIN2,LATT_FIN
    INTEGER NK
    INTEGER IDISPL      ! 0 no finite differences, 1 finite differences
    LOGICAL, OPTIONAL :: LOMEGA
    ! work arrays
    REAL(q),ALLOCATABLE :: DIST(:),XS(:),YS(:),ZS(:),VPS(:),YLM(:,:),VYLM(:)
    REAL(q) :: DISPL1(3,NONLR_S%NIONS),DISPL2(3,NONLR_S%NIONS)
    REAL(q) :: DISPL(3)
    TYPE (smoothing_handle) :: SH
    INTEGER :: ISH

!$  INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM,OMP_GET_MAX_THREADS
!$  INTEGER idim,ndim

    PROFILING_START('rspher_all')

#ifdef shmem_rproj
    ! before we zero NONLR_S%RPROJ we have to
    ! make sure no rank in COMM_shmem is still using it
    CALLMPI( M_barrier( NONLR_S%COMM_shmem ))
    ! only the first rank sets to zero
    IF (NONLR_S%COMM_shmem%NODE_ME==1) NONLR_S%RPROJ=0
    ! all have to wait until the first rank is
    ! finished zeroing NONLR_S%RPROJ
    CALLMPI( M_barrier( NONLR_S%COMM_shmem ))
#else
    NONLR_S%RPROJ=0
#endif

    LYDIM=MAXVAL(NONLR_S%LMAX)
    LMYDIM=(LYDIM+1)**2          ! number of lm pairs

    LMMAX =MAXVAL(NONLR_S%LMMAX) ! number of nlm indices in the non local potential
    IRMAX=NONLR_S%IRMAX

    CALL RSPHER_SMOOTH( SH, NONLR_S , GRID, LATT_CUR )

    ALLOCATE(DIST(IRMAX),XS(IRMAX),YS(IRMAX),ZS(IRMAX),VPS(IRMAX),YLM(IRMAX,LMYDIM),VYLM(IRMAX*LMMAX))

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

!$OMP PARALLEL DO NUM_THREADS(omp_nthreads_nonlr_rspace) &
!!$OMP PARALLEL DO &
!$OMP PRIVATE(idim, &
!$OMP         LYDIM,DIST,XS,YS,ZS,VPS,YLM,VYLM,ISH,NLIIND,NIS,NT,NI,ARGSC,F1,F2,F3,D1,D2,D3,N3LOW,N2LOW,N1LOW,N3HI,N2HI,N1HI, &
!$OMP         IDIS,LATT_FIN,DISPL,IND,N2,X2,N2P,N1,X1,N1P,NCOL,N3,X3,XC,YC,ZC,D,ARG,NADDR,X,Y,Z,N3P,ZZ,YY,XX,INDMAX, &
!$OMP         LMIND,L,FAKT,I,REM,LL,MMAX,LMBASE,LM,IBAS) &
!$OMP SHARED(ndim,LMYDIM,LMMAX,IRMAX,SH,NONLR_S,GRID,LATT_CUR,IDISPL,LATT_FIN1,LATT_FIN2,DISPL1,DISPL2,LOMEGA,omp_nonlr_planewise)
!$  omp: DO idim=1,ndim
    smooth: DO ISH=1,SH%N
!=======================================================================
! loop over all ions
!=======================================================================
      NLIIND=0
      NIS=1
      type: DO NT=1,NONLR_S%NTYP
         IF (NONLR_S%LMMAX(NT)==0) GOTO 600
         ions: DO NI=NIS,NONLR_S%NITYP(NT)+NIS-1
 
            IF (NONLR_S%SELECTED_ION>0 .AND. NONLR_S%SELECTED_ION/=NI) THEN
               NONLR_S%NLIMAX(NI __NOACC_omp_arg(:))=0
               CYCLE
            ENDIF

            ARGSC=NPSRNL/NONLR_S%PSRMAX(NT)
!=======================================================================
! find lattice points contained within the cutoff-sphere
! this loop might be done in scalar unit
!=======================================================================
            F1=1._q/GRID%NGX
            F2=1._q/GRID%NGY
            F3=1._q/GRID%NGZ

!-----------------------------------------------------------------------
! restrict loop to points contained within a cubus around the ion
!-----------------------------------------------------------------------
!sh
!            D1= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(1)*GRID%NGX
!            D2= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(2)*GRID%NGY
!            D3= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(3)*GRID%NGZ
            D1= (NONLR_S%PSRMAX(NT)+NONLR_S%RSMOOTH(NT))*LATT_CUR%BNORM(1)*GRID%NGX
            D2= (NONLR_S%PSRMAX(NT)+NONLR_S%RSMOOTH(NT))*LATT_CUR%BNORM(2)*GRID%NGY
            D3= (NONLR_S%PSRMAX(NT)+NONLR_S%RSMOOTH(NT))*LATT_CUR%BNORM(3)*GRID%NGZ
 
            N3LOW= INT(NONLR_S%POSION(3,NI)*GRID%NGZ-D3+10*GRID%NGZ+.99_q)-10*GRID%NGZ
            N2LOW= INT(NONLR_S%POSION(2,NI)*GRID%NGY-D2+10*GRID%NGY+.99_q)-10*GRID%NGY
            N1LOW= INT(NONLR_S%POSION(1,NI)*GRID%NGX-D1+10*GRID%NGX+.99_q)-10*GRID%NGX
 
            N3HI = INT(NONLR_S%POSION(3,NI)*GRID%NGZ+D3+10*GRID%NGZ)-10*GRID%NGZ
            N2HI = INT(NONLR_S%POSION(2,NI)*GRID%NGY+D2+10*GRID%NGY)-10*GRID%NGY
            N1HI = INT(NONLR_S%POSION(1,NI)*GRID%NGX+D1+10*GRID%NGX)-10*GRID%NGX
 
            VYLM= 0
 
            dis: DO IDIS=-ABS(IDISPL),ABS(IDISPL),2
 
               IF (IDIS==-1) THEN
                  LATT_FIN=LATT_FIN1
                  DISPL=DISPL1(:,NI)
               ELSE IF (IDIS==1) THEN
                  LATT_FIN=LATT_FIN2
                  DISPL=DISPL2(:,NI)
               ELSE
                  LATT_FIN=LATT_CUR
                  DISPL=0
               ENDIF
!-----------------------------------------------------------------------
! loop over cubus
! MPI version z ist the fast index
!-----------------------------------------------------------------------
               IF (GRID%RL%NFAST==3) THEN
               IND=1
 
               DO N2=N2LOW,N2HI
                  X2=(N2*F2-NONLR_S%POSION(2,NI))
                  N2P=MOD(N2+10*GRID%NGY,GRID%NGY)

!!$                IF (MOD(N2P,ndim)/=(idim-1)) CYCLE

                  DO N1=N1LOW,N1HI
                     X1=(N1*F1-NONLR_S%POSION(1,NI))
                     N1P=MOD(N1+10*GRID%NGX,GRID%NGX)
 
                     NCOL=GRID%RL%INDEX(N1P,N2P)
                     IF (NCOL==0) CYCLE ! not on local node go on
                     IF (GRID%RL%I2(NCOL) /= N1P+1 .OR. GRID%RL%I3(NCOL) /= N2P+1) THEN
                        CALL vtutor%bug("RSPHER: internal ERROR: " // str(GRID%RL%I2(NCOL)) // " "&
                            // str(N1P+1) // " " // str(GRID%RL%I3(NCOL)) // " " // str(N2P+1), &
                            __FILE__, __LINE__)
                     ENDIF

!$                   IF (MOD(NCOL,ndim)/=(idim-1)) CYCLE

!OCL SCALAR
                     DO N3=N3LOW,N3HI
                        X3=(N3*F3-NONLR_S%POSION(3,NI))
 
                        XC= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
                        YC= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
                        ZC= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)
 
                        D=SQRT(XC*XC+YC*YC+ZC*ZC)
                        ARG=(D*ARGSC)+1
                        NADDR=INT(ARG)
!sh                     IF (NADDR<NPSRNL) THEN
                        IF (D<(NPSRNL-1)/ARGSC+NONLR_S%RSMOOTH(NT)) THEN
                           X= (X1+SH%X1(ISH))*LATT_FIN%A(1,1)+(X2+SH%X2(ISH))*LATT_FIN%A(1,2)+(X3+SH%X3(ISH))*LATT_FIN%A(1,3)
                           Y= (X1+SH%X1(ISH))*LATT_FIN%A(2,1)+(X2+SH%X2(ISH))*LATT_FIN%A(2,2)+(X3+SH%X3(ISH))*LATT_FIN%A(2,3)
                           Z= (X1+SH%X1(ISH))*LATT_FIN%A(3,1)+(X2+SH%X2(ISH))*LATT_FIN%A(3,2)+(X3+SH%X3(ISH))*LATT_FIN%A(3,3)
!sh end
                           N3P=MOD(N3+10*GRID%NGZ,GRID%NGZ)
                           NONLR_S%NLI (IND,NI __NOACC_omp_arg(idim)) =1+N3P+ GRID%NGZ*(NCOL-1)
 
                           ZZ=Z-DISPL(3)
                           YY=Y-DISPL(2)
                           XX=X-DISPL(1)
                           ! the calculation of the | R(ion)-R(mesh)+d | for displaced ions
                           ! was done using the well known formula  | R+d | = | R | + d . R/|R|
                           ! this improves the stability of finite differences considerable
                           !IF (D<1E-4_q) THEN
                           !  DIST(IND)=1E-4_q
                           !ELSE
                           !  DIST(IND)=MAX(D-IDIS*(DISX*X+DISY*Y+DISZ*Z)/D,1E-10_q)
                           !ENDIF
                           DIST(IND)=MAX(SQRT(XX*XX+YY*YY+ZZ*ZZ),1E-10_q)
 
                           XS(IND)  =XX/DIST(IND)
                           YS(IND)  =YY/DIST(IND)
                           ZS(IND)  =ZZ/DIST(IND)
                           IND=IND+1
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
               ELSE
!-----------------------------------------------------------------------
! loop over cubus around one ion
! conventional version x is fast index
!-----------------------------------------------------------------------
               IND=1
               DO N3=N3LOW,N3HI
                  X3=(N3*F3-NONLR_S%POSION(3,NI))
                  N3P=MOD(N3+10*GRID%NGZ,GRID%NGZ)

!$                IF (MOD(N3P,ndim)/=(idim-1).AND.omp_nonlr_planewise) CYCLE

                  DO N2=N2LOW,N2HI
                     X2=(N2*F2-NONLR_S%POSION(2,NI))
                     N2P=MOD(N2+10*GRID%NGY,GRID%NGY)
 
                     NCOL=GRID%RL%INDEX(N2P,N3P)

!$                   IF (MOD(NCOL,ndim)/=(idim-1).AND.(.NOT.omp_nonlr_planewise)) CYCLE

                     DO N1=N1LOW,N1HI
                        X1=(N1*F1-NONLR_S%POSION(1,NI))
 
                        XC= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
                        YC= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
                        ZC= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)
 
                        D=SQRT(XC*XC+YC*YC+ZC*ZC)
!sh
!                        ARG=(D*ARGSC)+1
!                        NADDR=INT(ARG)
!                        IF (NADDR<NPSRNL) THEN
                        IF (D<(NPSRNL-1)/ARGSC+NONLR_S%RSMOOTH(NT)) THEN
                           X= (X1+SH%X1(ISH))*LATT_FIN%A(1,1)+(X2+SH%X2(ISH))*LATT_FIN%A(1,2)+(X3+SH%X3(ISH))*LATT_FIN%A(1,3)
                           Y= (X1+SH%X1(ISH))*LATT_FIN%A(2,1)+(X2+SH%X2(ISH))*LATT_FIN%A(2,2)+(X3+SH%X3(ISH))*LATT_FIN%A(2,3)
                           Z= (X1+SH%X1(ISH))*LATT_FIN%A(3,1)+(X2+SH%X2(ISH))*LATT_FIN%A(3,2)+(X3+SH%X3(ISH))*LATT_FIN%A(3,3)
!shend

                           N1P=MOD(N1+10*GRID%NGX,GRID%NGX)
                           NONLR_S%NLI (IND,NI __NOACC_omp_arg(idim)) =N1P+(NCOL-1)*GRID%NGX+1
                           IF (NONLR_S%NLI (IND,NI __NOACC_omp_arg(idim)) /= 1+N1P+GRID%NGX*(N2P+GRID%NGY* N3P)) THEN
                              CALL vtutor%bug("RSHPER internal ERROR: " // str(N1P) // " " // &
                                 str(N2P) // " " // str(N3P) // " " // str(NCOL), __FILE__, __LINE__)
                           ENDIF

                           ZZ=Z-DISPL(3)
                           YY=Y-DISPL(2)
                           XX=X-DISPL(1)
                           ! the calculation of the | R(ion)-R(mesh)+d | for displaced ions
                           ! was done using the well known formula  | R+d | = | R | + d . R/|R|
                           ! this improves the stability of finite differences considerable
                           !IF (D<1E-4_q) THEN
                           !  DIST(IND)=1E-4_q
                           !ELSE
                           !  DIST(IND)=MAX(D-IDIS*(DISX*X+DISY*Y+DISZ*Z)/D,1E-10_q)
                           !ENDIF
                           DIST(IND)=MAX(SQRT(XX*XX+YY*YY+ZZ*ZZ),1E-10_q)
 
                           XS(IND)  =XX/DIST(IND)
                           YS(IND)  =YY/DIST(IND)
                           ZS(IND)  =ZZ/DIST(IND)
                           IND=IND+1
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
               ENDIF
!-----------------------------------------------------------------------
!  compare maximum index with INDMAX
!-----------------------------------------------------------------------
               INDMAX=IND-1
               IF (INDMAX>NONLR_S%IRMAX) THEN
                  CALL vtutor%bug("internal ERROR: RSPHER: NONLR_S%IRMAX must be increased to " // &
                     str(INT(INDMAX*1.1_q)), __FILE__, __LINE__)
               ENDIF
               NONLR_S%NLIMAX(NI __NOACC_omp_arg(idim))=INDMAX

               IF (INDMAX==0) CYCLE ions
!=======================================================================
! now calculate the tables containing the spherical harmonics
! multiplied by the pseudopotential
!=======================================================================
               LYDIM=NONLR_S%LMAX(NT)
               CALL SETYLM(LYDIM,INDMAX,YLM,XS,YS,ZS)
 
               LMIND=1
               l_loop: DO L=1,NONLR_S%CHANNELS(NT)
!-----------------------------------------------------------------------
! interpolate the non-local pseudopotentials
! and multiply by (LATT_CUR%OMEGA)^(1/2)
! interpolation is done here using spline-fits this inproves the
! numerical stability of the forces the MIN operation takes care
! that the index is between  1 and NPSRNL
!-----------------------------------------------------------------------
                  FAKT= SQRT(LATT_FIN%OMEGA)
                  IF (PRESENT (LOMEGA)) THEN
                     IF (LOMEGA) FAKT=LATT_FIN%OMEGA
                  ENDIF                   
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(I,REM)
                  DO IND=1,INDMAX
                     I  =MIN(INT(DIST(IND)*ARGSC)+1,NPSRNL-1)
                     REM=DIST(IND)-NONLR_S%BETA(NT)%PSPRNL(I,1,L)
                     VPS(IND)=(NONLR_S%BETA(NT)%PSPRNL(I,2,L)+REM*(NONLR_S%BETA(NT)%PSPRNL(I,3,L)+ &
                          &         REM*(NONLR_S%BETA(NT)%PSPRNL(I,4,L)+REM*NONLR_S%BETA(NT)%PSPRNL(I,5,L))))*FAKT
                  ENDDO
 
                  LL=NONLR_S%BETA(NT)%LPS(L)
                  MMAX=2*LL
 
                  ! invert sign for first displacement
                  IF (IDIS==-1) THEN
                     DO IND=1,INDMAX
                        VPS(IND)=-VPS(IND)
                     ENDDO
                  ENDIF
 
                  LMBASE=LL**2+1
 
                  DO LM=0,MMAX
                     DO IND=1,INDMAX
                        IBAS = (LMIND-1+LM)*INDMAX
                        VYLM(IBAS+IND)=VYLM(IBAS+IND)+VPS(IND)*YLM(IND,LM+LMBASE)
                     ENDDO
                  ENDDO
 
                  LMIND=LMIND+MMAX+1
               ENDDO l_loop
 
               IF (LMIND-1/=NONLR_S%LMMAX(NT)) THEN
                  CALL vtutor%bug("internal ERROR: RSPHER: NONLR_S%LMMAX is wrong " // str(LMIND-1) &
                     // " " // str(NONLR_S%LMMAX(NT)), __FILE__, __LINE__)
               ENDIF
 
               IF ( NONLR_S%LMMAX(NT)*INDMAX+NLIIND > NONLR_S%IRALLOC) THEN
                       CALL vtutor%bug("internal ERROR RSPHER: running out of buffer " // &
                          str(NLIIND) // " " // str(INDMAX) // " " // str(NONLR_S%LMMAX(NT)) // " " &
                          // str(NT) // " " // str(NONLR_S%IRALLOC) // "\nnonlr.F:Out of buffer RSPHER", &
                          __FILE__, __LINE__)
               ENDIF

            ENDDO dis
!-----------------------------------------------------------------------
! finally store the coefficients
!-----------------------------------------------------------------------
            DO LMIND=1,NONLR_S%LMMAX(NT)
#ifdef shmem_rproj
               IF (MOD(LMIND-1,NONLR_S%COMM_shmem%NCPU)/=NONLR_S%COMM_shmem%NODE_ME-1) CYCLE
#endif
               DO IND=1,INDMAX
                  IBAS = (LMIND-1)*INDMAX
!sh
                  NONLR_S%RPROJ(IND+IBAS+NLIIND __NOACC_omp_arg(idim))=NONLR_S%RPROJ(IND+IBAS+NLIIND __NOACC_omp_arg(idim))+VYLM(IND+IBAS)*SH%WEIGHT(ISH)
               ENDDO
            ENDDO

            NLIIND= NONLR_S%LMMAX(NT)*INDMAX+NLIIND
!=======================================================================
! end of loop over ions
!=======================================================================
         ENDDO ions
  600    NIS = NIS+NONLR_S%NITYP(NT)
      ENDDO type
    ENDDO smooth

!=======================================================================
! Here we set LMBASE and NLIBASE. These arrays of offsets allow
! us to contract nested loops over types+ions_of_type to single loops
! over all ions (used lateron in for instance RPROMU_HF).
!=======================================================================
    NONLR_S%NLIBASE(: __NOACC_omp_arg(idim))=0
!$  IF (idim==1) THEN
       NONLR_S%LMBASE(1)=0
!$  ENDIF
    DO NI=1,NONLR_S%NIONS
       NT=NONLR_S%ITYP(NI)
!$     IF (idim==1) THEN
           NONLR_S%LMBASE(NI+1)=NONLR_S%LMBASE(NI)+NONLR_S%LMMAX(NT)
!$     ENDIF
       IF (NONLR_S%SELECTED_ION>0 .AND. NONLR_S%SELECTED_ION/=NI) THEN
          NONLR_S%NLIBASE(NI+1 __NOACC_omp_arg(idim))=NONLR_S%NLIBASE(NI __NOACC_omp_arg(idim))
       ELSE
          NONLR_S%NLIBASE(NI+1 __NOACC_omp_arg(idim))=NONLR_S%NLIBASE(NI __NOACC_omp_arg(idim))+NONLR_S%LMMAX(NT)*NONLR_S%NLIMAX(NI __NOACC_omp_arg(idim))
       ENDIF
    ENDDO
!$  ENDDO omp
!$OMP END PARALLEL DO

    DEALLOCATE(DIST,XS,YS,ZS,VPS,YLM,VYLM)

    CALL RSPHER_SMOOTH_DEALLOCATE( SH )
#ifdef shmem_rproj
    ! and let's make sure everybody is done
    CALLMPI( M_barrier( NONLR_S%COMM_shmem ))
#endif
!$ACC UPDATE DEVICE(NONLR_S%NLIBASE,NONLR_S%LMBASE,NONLR_S%RPROJ,NONLR_S%NLIMAX,NONLR_S%NLI) __IF_ASYNC__ IF_PRESENT
    PROFILING_STOP('rspher_all')

  END SUBROUTINE RSPHER_ALL


!****************** SUBROUTINE RSPHER_SMOOTH ***************************
!
! subroutine to smooth the real space projector functions
! using a simple real space method
!
!***********************************************************************

  SUBROUTINE RSPHER_SMOOTH( SH, NONLR_S, GRID, LATT_CUR)
    USE lattice
    USE constant
    IMPLICIT NONE
    TYPE (smoothing_handle) SH
    TYPE (nonlr_struct) NONLR_S
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    
    INTEGER :: N1, N2, N3, I1, I2, I3, NT
    INTEGER, PARAMETER :: IREFINE=4
    REAL(q) :: F1, F2, F3, X, Y, Z, D
    INTEGER, PARAMETER :: NMAX=2000
    REAL(q) :: X1(NMAX), X2(NMAX), X3(NMAX), WEIGHT(NMAX)
    REAL(q) :: RSMOOTH

! currently only one fixed value for RSMOOTH is allowed
! check that they are all the same

    RSMOOTH=NONLR_S%RSMOOTH(1)
 
    DO NT=1,SIZE(NONLR_S%RSMOOTH)
       IF (RSMOOTH/=NONLR_S%RSMOOTH(NT)) THEN
          CALL vtutor%bug("RSPHER_SMOOTH: internal error only one value for RSMOOTH allowed " // &
             str(RSMOOTH) // " " // str(NONLR_S%RSMOOTH(NT)) // " " // str(NT), __FILE__, __LINE__)
       ENDIF
    ENDDO
 
    IF (RSMOOTH==0) THEN
       SH%N=1
 
       ALLOCATE(SH%WEIGHT(SH%N))
       ALLOCATE(SH%X1(SH%N), SH%X2(SH%N), SH%X3(SH%N))
 
       SH%N=1
       SH%X1=0
       SH%X2=0
       SH%X3=0
       SH%WEIGHT=1
    ELSE
       N1= RSMOOTH*LATT_CUR%BNORM(1)*IREFINE*GRID%NGX
       N2= RSMOOTH*LATT_CUR%BNORM(2)*IREFINE*GRID%NGY
       N3= RSMOOTH*LATT_CUR%BNORM(3)*IREFINE*GRID%NGZ
 
       SH%N=0
 
       DO I1=-N1,N1
          F1=REAL(I1,q)/GRID%NGX/IREFINE
          DO I2=-N2,N2
             F2=REAL(I2,q)/GRID%NGY/IREFINE
             DO I3=-N3,N3
                F3=REAL(I3,q)/GRID%NGZ/IREFINE
                X= F1*LATT_CUR%A(1,1)+F2*LATT_CUR%A(1,2)+F3*LATT_CUR%A(1,3)
                Y= F1*LATT_CUR%A(2,1)+F2*LATT_CUR%A(2,2)+F3*LATT_CUR%A(2,3)
                Z= F1*LATT_CUR%A(3,1)+F2*LATT_CUR%A(3,2)+F3*LATT_CUR%A(3,3)
                
                D=SQRT(X*X+Y*Y+Z*Z)
                IF (D<RSMOOTH) THEN
                   SH%N=SH%N+1
                   IF( (SH%N) > NMAX) THEN
                      CALL vtutor%bug("internal error in RSPHER_SMOOTH: increase NMAX", __FILE__, __LINE__)
                   ENDIF
                   X1(SH%N)=F1
                   X2(SH%N)=F2
                   X3(SH%N)=F3
!                   WEIGHT(SH%N)= (COS(D/RSMOOTH*PI)+1)/2
                   WEIGHT(SH%N)= EXP(-4*(D/RSMOOTH)**2)
!                   WRITE(*,'(5F14.7)') D,F1,F2,F3

                ENDIF
             ENDDO
          ENDDO
       ENDDO
       ALLOCATE(SH%WEIGHT(SH%N))
       ALLOCATE(SH%X1(SH%N), SH%X2(SH%N), SH%X3(SH%N))
 
       SH%X1=X1(1:SH%N)
       SH%X2=X2(1:SH%N)
       SH%X3=X3(1:SH%N)
       SH%WEIGHT=WEIGHT(1:SH%N)/SUM(WEIGHT(1:SH%N))
       DO I1=1,SH%N
          X= SH%X1(I1)*LATT_CUR%A(1,1)+SH%X2(I1)*LATT_CUR%A(1,2)+SH%X3(I1)*LATT_CUR%A(1,3)
          Y= SH%X1(I1)*LATT_CUR%A(2,1)+SH%X2(I1)*LATT_CUR%A(2,2)+SH%X3(I1)*LATT_CUR%A(2,3)
          Z= SH%X1(I1)*LATT_CUR%A(3,1)+SH%X2(I1)*LATT_CUR%A(3,2)+SH%X3(I1)*LATT_CUR%A(3,3)
          D=SQRT(X*X+Y*Y+Z*Z)
!          WRITE(*,'(5F14.7)') D,WEIGHT(I1),SH%X1(I1),SH%X2(I1),SH%X3(I1)
       ENDDO
       WRITE(*,*) 'grid is refined by ',IREFINE,SH%N
    END IF

  END SUBROUTINE RSPHER_SMOOTH


  SUBROUTINE RSPHER_SMOOTH_DEALLOCATE( SH )

    TYPE (smoothing_handle) :: SH

    DEALLOCATE(SH%WEIGHT)
    DEALLOCATE(SH%X1, SH%X2, SH%X3)
  
  END SUBROUTINE RSPHER_SMOOTH_DEALLOCATE


!****************** SUBROUTINE PHASER **********************************
!
! subroutine PHASER
! recalculates the phase factor for the real-space projectors
!
!> @details @ref openmp : Partition each real space projection operator
!> into openmp::omp_nthreads_nonlr_rspace parts. The work on these
!> parts will be distributed over the available OpenMP threads lateron.
!> @n
!> Currently OMP_NUM_THREADS/=1 requires NCORE=1. Consequently
!> the fast index of the real space grid of the non local projection
!> operators is the x-index (GRID\%RL\%NFAST/=3).
!> @n
!> For openmp::omp_nonlr_planewise = .TRUE. (default) the projectors are
!> distributed over (x,y)-planes in a round robin fashion.
!> In case openmp::omp_nonlr_planewise = .FALSE. the projectors are
!> distributed over (z,y)-columns.
!
!***********************************************************************

  SUBROUTINE PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES)
    USE lattice
    USE mpimy
    USE constant
    USE pseudo
    USE openmp, ONLY : omp_nthreads_nonlr_rspace,omp_nonlr_planewise
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes)     WDES

!$  INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
!$  INTEGER idim,ndim

    PROFILING_START('phaser')

    NONLR_S%NK=NK
!$ACC UPDATE DEVICE(NONLR_S%NK) __IF_ASYNC__

#ifdef gammareal
    PROFILING_STOP('phaser')
    RETURN
#else
NOACC !$  ndim=SIZE(NONLR_S%CRREXP,4)
DOACC !$  ndim=1

!$OMP PARALLEL DO NUM_THREADS(omp_nthreads_nonlr_rspace) &
!!$OMP PARALLEL DO &
!$OMP PRIVATE(idim, &
!$OMP         VKX,VKY,VKZ,QX,QY,QZ,NSPINORS,ISPINOR,NIS,NT,NI,ARGSC,F1,F2,F3,D1,D2,D3,N3LOW,N2LOW,N1LOW,N3HI,N2HI,N1HI, &
!$OMP         IND,N2,X2,N2P,N1,X1,N1P,NCOL,N3,X3,X,Y,Z,D,ARG,NADDR,N3P) &
!$OMP SHARED(ndim,GRID,LATT_CUR,NONLR_S,NK,WDES,omp_nonlr_planewise)
!$ omp: DO idim=1,ndim
!-----------------------------------------------------------------------
! k-point in kartesian coordiantes
!-----------------------------------------------------------------------
    VKX= WDES%VKPT(1,NK)*LATT_CUR%B(1,1)+WDES%VKPT(2,NK)*LATT_CUR%B(1,2)+WDES%VKPT(3,NK)*LATT_CUR%B(1,3)
    VKY= WDES%VKPT(1,NK)*LATT_CUR%B(2,1)+WDES%VKPT(2,NK)*LATT_CUR%B(2,2)+WDES%VKPT(3,NK)*LATT_CUR%B(2,3)
    VKZ= WDES%VKPT(1,NK)*LATT_CUR%B(3,1)+WDES%VKPT(2,NK)*LATT_CUR%B(3,2)+WDES%VKPT(3,NK)*LATT_CUR%B(3,3)

!-----------------------------------------------------------------------
! spin spiral propagation vector in cartesian coordinates
! is simply zero when LSPIRAL=.FALSE.
!-----------------------------------------------------------------------
    QX= (WDES%QSPIRAL(1)*LATT_CUR%B(1,1)+WDES%QSPIRAL(2)*LATT_CUR%B(1,2)+WDES%QSPIRAL(3)*LATT_CUR%B(1,3))/2
    QY= (WDES%QSPIRAL(1)*LATT_CUR%B(2,1)+WDES%QSPIRAL(2)*LATT_CUR%B(2,2)+WDES%QSPIRAL(3)*LATT_CUR%B(2,3))/2
    QZ= (WDES%QSPIRAL(1)*LATT_CUR%B(3,1)+WDES%QSPIRAL(2)*LATT_CUR%B(3,2)+WDES%QSPIRAL(3)*LATT_CUR%B(3,3))/2

!=======================================================================
! Loop over NSPINORS: here only in case of spin spirals NRSPINOR=2
!=======================================================================
    IF (NONLR_S%LSPIRAL) THEN 
       NSPINORS=2
    ELSE
       NSPINORS=1
    ENDIF

    spinor: DO ISPINOR=1,NSPINORS
!=======================================================================
! loop over all ions
!=======================================================================
       NIS=1

!OCL SCALAR
       type: DO NT=1,NONLR_S%NTYP
          IF (NONLR_S%LMMAX(NT)==0) GOTO 600
          ions: DO NI=NIS,NONLR_S%NITYP(NT)+NIS-1

             IF (NONLR_S%SELECTED_ION>0 .AND. NONLR_S%SELECTED_ION/=NI) CYCLE

             IF (ASSOCIATED(NONLR_S%VKPT_SHIFT)) THEN
                VKX= (WDES%VKPT(1,NK)+NONLR_S%VKPT_SHIFT(1,NI))*LATT_CUR%B(1,1)+ & 
                     (WDES%VKPT(2,NK)+NONLR_S%VKPT_SHIFT(2,NI))*LATT_CUR%B(1,2)+ & 
                     (WDES%VKPT(3,NK)+NONLR_S%VKPT_SHIFT(3,NI))*LATT_CUR%B(1,3)
                VKY= (WDES%VKPT(1,NK)+NONLR_S%VKPT_SHIFT(1,NI))*LATT_CUR%B(2,1)+ & 
                     (WDES%VKPT(2,NK)+NONLR_S%VKPT_SHIFT(2,NI))*LATT_CUR%B(2,2)+ & 
                     (WDES%VKPT(3,NK)+NONLR_S%VKPT_SHIFT(3,NI))*LATT_CUR%B(2,3)
                VKZ= (WDES%VKPT(1,NK)+NONLR_S%VKPT_SHIFT(1,NI))*LATT_CUR%B(3,1)+ & 
                     (WDES%VKPT(2,NK)+NONLR_S%VKPT_SHIFT(2,NI))*LATT_CUR%B(3,2)+ &
                     (WDES%VKPT(3,NK)+NONLR_S%VKPT_SHIFT(3,NI))*LATT_CUR%B(3,3)
             ENDIF
!-----------------------------------------------------------------------
! check some quantities
!-----------------------------------------------------------------------
             ARGSC=NPSRNL/NONLR_S%PSRMAX(NT)

             F1=1._q/GRID%NGX
             F2=1._q/GRID%NGY
             F3=1._q/GRID%NGZ

!-----------------------------------------------------------------------
! restrict loop to points contained within a cubus around the ion
!-----------------------------------------------------------------------
!sh
!            D1= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(1)*GRID%NGX
!            D2= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(2)*GRID%NGY
!            D3= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(3)*GRID%NGZ
             D1= (NONLR_S%PSRMAX(NT)+NONLR_S%RSMOOTH(NT))*LATT_CUR%BNORM(1)*GRID%NGX
             D2= (NONLR_S%PSRMAX(NT)+NONLR_S%RSMOOTH(NT))*LATT_CUR%BNORM(2)*GRID%NGY
             D3= (NONLR_S%PSRMAX(NT)+NONLR_S%RSMOOTH(NT))*LATT_CUR%BNORM(3)*GRID%NGZ

             N3LOW= INT(NONLR_S%POSION(3,NI)*GRID%NGZ-D3+GRID%NGZ+.99_q)-GRID%NGZ
             N2LOW= INT(NONLR_S%POSION(2,NI)*GRID%NGY-D2+GRID%NGY+.99_q)-GRID%NGY
             N1LOW= INT(NONLR_S%POSION(1,NI)*GRID%NGX-D1+GRID%NGX+.99_q)-GRID%NGX

             N3HI = INT(NONLR_S%POSION(3,NI)*GRID%NGZ+D3)
             N2HI = INT(NONLR_S%POSION(2,NI)*GRID%NGY+D2)
             N1HI = INT(NONLR_S%POSION(1,NI)*GRID%NGX+D1)
!-----------------------------------------------------------------------
! loop over cubus
! MPI version z ist the fast index
!-----------------------------------------------------------------------
             IF (GRID%RL%NFAST==3) THEN
             IND=1

             DO N2=N2LOW,N2HI
                X2=(N2*F2-NONLR_S%POSION(2,NI))
                N2P=MOD(N2+10*GRID%NGY,GRID%NGY)

!!$              IF (MOD(N2P,ndim)/=(idim-1)) CYCLE

                DO N1=N1LOW,N1HI
                   X1=(N1*F1-NONLR_S%POSION(1,NI))
                   N1P=MOD(N1+10*GRID%NGX,GRID%NGX)

                   NCOL=GRID%RL%INDEX(N1P,N2P)
                   IF (NCOL==0) CYCLE
                   IF (GRID%RL%I2(NCOL) /= N1P+1 .OR. GRID%RL%I3(NCOL) /= N2P+1) THEN
                      CALL vtutor%bug("internal ERROR PHASER: " // str(GRID%RL%I2(NCOL)) // " " // &
                         str(N1P+1) // " " // str(GRID%RL%I3(NCOL)) // " " // str(N2P+1), __FILE__, __LINE__)
                   ENDIF

!$                 IF (MOD(NCOL,ndim)/=(idim-1)) CYCLE

                   DO N3=N3LOW,N3HI
                      X3=(N3*F3-NONLR_S%POSION(3,NI))

                      X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
                      Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
                      Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)

                      D=SQRT(X*X+Y*Y+Z*Z)
                      ARG=(D*ARGSC)+1
                      NADDR=INT(ARG)

!sh                   IF (NADDR<NPSRNL) THEN
                      IF (D<(NPSRNL-1)/ARGSC+NONLR_S%RSMOOTH(NT)) THEN
                         NONLR_S%CRREXP(IND,NI,ISPINOR __NOACC_omp_arg(idim))=EXP(CITPI*(X*(VKX-QX)+Y*(VKY-QY)+Z*(VKZ-QZ)))
                         IND=IND+1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             ELSE
!-----------------------------------------------------------------------
! loop over cubus around one ion
! conventional version x is fast index
!-----------------------------------------------------------------------
             IND=1
             DO N3=N3LOW,N3HI
                X3=(N3*F3-NONLR_S%POSION(3,NI))
                N3P=MOD(N3+10*GRID%NGZ,GRID%NGZ)

!$              IF (MOD(N3P,ndim)/=(idim-1).AND.omp_nonlr_planewise) CYCLE

                DO N2=N2LOW,N2HI
                   X2=(N2*F2-NONLR_S%POSION(2,NI))
                   N2P=MOD(N2+10*GRID%NGY,GRID%NGY)

                   NCOL=GRID%RL%INDEX(N2P,N3P)

!$                 IF (MOD(NCOL,ndim)/=(idim-1).AND.(.NOT.omp_nonlr_planewise)) CYCLE

                   DO N1=N1LOW,N1HI
                      X1=(N1*F1-NONLR_S%POSION(1,NI))

                      X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
                      Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
                      Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)

                      D=SQRT(X*X+Y*Y+Z*Z)
                      ARG=(D*ARGSC)+1
                      NADDR=INT(ARG)

!sh                   IF (NADDR<NPSRNL) THEN
                      IF (D<(NPSRNL-1)/ARGSC+NONLR_S%RSMOOTH(NT)) THEN
                         NONLR_S%CRREXP(IND,NI,ISPINOR __NOACC_omp_arg(idim))=EXP(CITPI*(X*(VKX-QX)+Y*(VKY-QY)+Z*(VKZ-QZ)))
                         IND=IND+1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             ENDIF
!=======================================================================
! end of loop over ions
!=======================================================================
          ENDDO ions
600       NIS = NIS+NONLR_S%NITYP(NT)
       ENDDO type

       ! conjugate phase alteration for spin down: -q/2 -> q/2
       QX=-QX
       QY=-QY
       QZ=-QZ
    ENDDO spinor
!$ ENDDO omp
!$OMP END PARALLEL DO
#endif

!$ACC UPDATE DEVICE(NONLR_S%CRREXP) __IF_ASYNC__

    PROFILING_STOP('phaser')

    RETURN
  END SUBROUTINE PHASER

!****************** SUBROUTINE PHASERR  ********************************
!
! recalculates the phase factor e^(ik r-R_i) times x,y,z
! the cartesian direction, selecting x, y or z,  is supplied by an 
! index IDIR
!
!> @details @ref openmp : Partition each real space projection operator
!> into openmp::omp_nthreads_nonlr_rspace parts. The work on these
!> parts will be distributed over the available OpenMP threads lateron.
!> @n
!> Currently OMP_NUM_THREADS/=1 requires NCORE=1. Consequently
!> the fast index of the real space grid of the non local projection
!> operators is the x-index (GRID\%RL\%NFAST/=3).
!> @n
!> For openmp::omp_nonlr_planewise = .TRUE. (default) the projectors are
!> distributed over (x,y)-planes in a round robin fashion.
!> In case openmp::omp_nonlr_planewise = .FALSE. the projectors are
!> distributed over (z,y)-columns.
!
!***********************************************************************

  SUBROUTINE PHASERR(GRID,LATT_CUR,NONLR_S,NK,WDES,IDIR)
    USE lattice
    USE mpimy
    USE constant
    USE pseudo
    USE openmp, ONLY : omp_nthreads_nonlr_rspace,omp_nonlr_planewise
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes)     WDES
    INTEGER IDIR

    REAL (q) XX(3)

!$  INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
!$  INTEGER idim,ndim

    NONLR_S%NK=NK

    IF (.NOT. ASSOCIATED(NONLR_S%CRREXP)) THEN
       CALL vtutor%bug("internal error in PHASERR: CRREXP is not allocated", __FILE__, __LINE__)
    ENDIF

NOACC !$  ndim=SIZE(NONLR_S%CRREXP,4)
DOACC !$  ndim=1

!$OMP PARALLEL DO NUM_THREADS(omp_nthreads_nonlr_rspace) &
!!$OMP PARALLEL DO &
!$OMP PRIVATE(idim, &
!$OMP         VKX,VKY,VKZ,QX,QY,QZ,NSPINORS,ISPINOR,NIS,NT,NI,ARGSC,F1,F2,F3,D1,D2,D3,N3LOW,N2LOW,N1LOW,N3HI,N2HI,N1HI, &
!$OMP         IND,N2,X2,N2P,N1,X1,N1P,NCOL,N3,X3,X,Y,Z,XX,D,ARG,NADDR,N3P) &
!$OMP SHARED(ndim,GRID,LATT_CUR,NONLR_S,NK,WDES,IDIR,omp_nonlr_planewise)
!$ omp: DO idim=1,ndim
!-----------------------------------------------------------------------
! k-point in kartesian coordiantes
!-----------------------------------------------------------------------
    VKX= WDES%VKPT(1,NK)*LATT_CUR%B(1,1)+WDES%VKPT(2,NK)*LATT_CUR%B(1,2)+WDES%VKPT(3,NK)*LATT_CUR%B(1,3)
    VKY= WDES%VKPT(1,NK)*LATT_CUR%B(2,1)+WDES%VKPT(2,NK)*LATT_CUR%B(2,2)+WDES%VKPT(3,NK)*LATT_CUR%B(2,3)
    VKZ= WDES%VKPT(1,NK)*LATT_CUR%B(3,1)+WDES%VKPT(2,NK)*LATT_CUR%B(3,2)+WDES%VKPT(3,NK)*LATT_CUR%B(3,3)

!-----------------------------------------------------------------------
! spin spiral propagation vector in cartesian coordinates
! is simply zero when LSPIRAL=.FALSE.
!-----------------------------------------------------------------------
    QX= (WDES%QSPIRAL(1)*LATT_CUR%B(1,1)+WDES%QSPIRAL(2)*LATT_CUR%B(1,2)+WDES%QSPIRAL(3)*LATT_CUR%B(1,3))/2
    QY= (WDES%QSPIRAL(1)*LATT_CUR%B(2,1)+WDES%QSPIRAL(2)*LATT_CUR%B(2,2)+WDES%QSPIRAL(3)*LATT_CUR%B(2,3))/2
    QZ= (WDES%QSPIRAL(1)*LATT_CUR%B(3,1)+WDES%QSPIRAL(2)*LATT_CUR%B(3,2)+WDES%QSPIRAL(3)*LATT_CUR%B(3,3))/2

!=======================================================================
! Loop over NSPINORS: here only in case of spin spirals NRSPINOR=2
!=======================================================================
    IF (NONLR_S%LSPIRAL) THEN 
       NSPINORS=2
    ELSE
       NSPINORS=1
    ENDIF

    spinor: DO ISPINOR=1,NSPINORS
       NIS=1

!OCL SCALAR
       type: DO NT=1,NONLR_S%NTYP
          IF (NONLR_S%LMMAX(NT)==0) GOTO 600
          ions: DO NI=NIS,NONLR_S%NITYP(NT)+NIS-1

             IF (NONLR_S%SELECTED_ION>0 .AND. NONLR_S%SELECTED_ION/=NI) CYCLE
             ARGSC=NPSRNL/NONLR_S%PSRMAX(NT)

             IF (ASSOCIATED(NONLR_S%VKPT_SHIFT)) THEN
                VKX= (WDES%VKPT(1,NK)+NONLR_S%VKPT_SHIFT(1,NI))*LATT_CUR%B(1,1)+ & 
                     (WDES%VKPT(2,NK)+NONLR_S%VKPT_SHIFT(2,NI))*LATT_CUR%B(1,2)+ & 
                     (WDES%VKPT(3,NK)+NONLR_S%VKPT_SHIFT(3,NI))*LATT_CUR%B(1,3)
                VKY= (WDES%VKPT(1,NK)+NONLR_S%VKPT_SHIFT(1,NI))*LATT_CUR%B(2,1)+ & 
                     (WDES%VKPT(2,NK)+NONLR_S%VKPT_SHIFT(2,NI))*LATT_CUR%B(2,2)+ & 
                     (WDES%VKPT(3,NK)+NONLR_S%VKPT_SHIFT(3,NI))*LATT_CUR%B(2,3)
                VKZ= (WDES%VKPT(1,NK)+NONLR_S%VKPT_SHIFT(1,NI))*LATT_CUR%B(3,1)+ & 
                     (WDES%VKPT(2,NK)+NONLR_S%VKPT_SHIFT(2,NI))*LATT_CUR%B(3,2)+ &
                     (WDES%VKPT(3,NK)+NONLR_S%VKPT_SHIFT(3,NI))*LATT_CUR%B(3,3)
             ENDIF
!=======================================================================
! find lattice points contained within the cutoff-sphere
! this loop might be done in scalar unit
!=======================================================================
             F1=1._q/GRID%NGX
             F2=1._q/GRID%NGY
             F3=1._q/GRID%NGZ

!-----------------------------------------------------------------------
! restrict loop to points contained within a cubus around the ion
!-----------------------------------------------------------------------
!sh
!            D1= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(1)*GRID%NGX
!            D2= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(2)*GRID%NGY
!            D3= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(3)*GRID%NGZ
             D1= (NONLR_S%PSRMAX(NT)+NONLR_S%RSMOOTH(NT))*LATT_CUR%BNORM(1)*GRID%NGX
             D2= (NONLR_S%PSRMAX(NT)+NONLR_S%RSMOOTH(NT))*LATT_CUR%BNORM(2)*GRID%NGY
             D3= (NONLR_S%PSRMAX(NT)+NONLR_S%RSMOOTH(NT))*LATT_CUR%BNORM(3)*GRID%NGZ

             N3LOW= INT(NONLR_S%POSION(3,NI)*GRID%NGZ-D3+GRID%NGZ+.99_q)-GRID%NGZ
             N2LOW= INT(NONLR_S%POSION(2,NI)*GRID%NGY-D2+GRID%NGY+.99_q)-GRID%NGY
             N1LOW= INT(NONLR_S%POSION(1,NI)*GRID%NGX-D1+GRID%NGX+.99_q)-GRID%NGX

             N3HI = INT(NONLR_S%POSION(3,NI)*GRID%NGZ+D3)
             N2HI = INT(NONLR_S%POSION(2,NI)*GRID%NGY+D2)
             N1HI = INT(NONLR_S%POSION(1,NI)*GRID%NGX+D1)

!-----------------------------------------------------------------------
! loop over cubus
! MPI version z ist the fast index
!-----------------------------------------------------------------------
             IF (GRID%RL%NFAST==3) THEN
             IND=1

             DO N2=N2LOW,N2HI
                X2=(N2*F2-NONLR_S%POSION(2,NI))
                N2P=MOD(N2+10*GRID%NGY,GRID%NGY)

!!$              IF (MOD(N2P,ndim)/=(idim-1)) CYCLE

                DO N1=N1LOW,N1HI
                   X1=(N1*F1-NONLR_S%POSION(1,NI))
                   N1P=MOD(N1+10*GRID%NGX,GRID%NGX)

                   NCOL=GRID%RL%INDEX(N1P,N2P)
                   IF (NCOL==0) CYCLE
                   IF (GRID%RL%I2(NCOL) /= N1P+1 .OR. GRID%RL%I3(NCOL) /= N2P+1) THEN
                      CALL vtutor%bug("internal ERROR PHASER: " // str(GRID%RL%I2(NCOL)) // " " // &
                         str(N1P+1) // " " // str(GRID%RL%I3(NCOL)) // " " // str(N2P+1), __FILE__, __LINE__)
                   ENDIF

!$                 IF (MOD(NCOL,ndim)/=(idim-1)) CYCLE

                   DO N3=N3LOW,N3HI
                      X3=(N3*F3-NONLR_S%POSION(3,NI))

                      X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
                      Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
                      Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)
                      XX(1)=X
                      XX(2)=Y
                      XX(3)=Z

                      D=SQRT(X*X+Y*Y+Z*Z)
                      ARG=(D*ARGSC)+1
                      NADDR=INT(ARG)

!sh                   IF (NADDR<NPSRNL) THEN
                      IF (D<(NPSRNL-1)/ARGSC+NONLR_S%RSMOOTH(NT)) THEN
                         NONLR_S%CRREXP(IND,NI,ISPINOR __NOACC_omp_arg(idim))=EXP(CITPI*(X*(VKX-QX)+Y*(VKY-QY)+Z*(VKZ-QZ)))*XX(IDIR)
                         IND=IND+1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             ELSE
!-----------------------------------------------------------------------
! loop over cubus around one ion
! conventional version x is fast index
!-----------------------------------------------------------------------
             IND=1
             DO N3=N3LOW,N3HI
                X3=(N3*F3-NONLR_S%POSION(3,NI))
                N3P=MOD(N3+10*GRID%NGZ,GRID%NGZ)

!$              IF (MOD(N3P,ndim)/=(idim-1).AND.omp_nonlr_planewise) CYCLE

                DO N2=N2LOW,N2HI
                   X2=(N2*F2-NONLR_S%POSION(2,NI))
                   N2P=MOD(N2+10*GRID%NGY,GRID%NGY)

                   NCOL=GRID%RL%INDEX(N2P,N3P)

!$                 IF (MOD(NCOL,ndim)/=(idim-1).AND.(.NOT.omp_nonlr_planewise)) CYCLE

                   DO N1=N1LOW,N1HI
                      X1=(N1*F1-NONLR_S%POSION(1,NI))

                      X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
                      Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
                      Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)
                      XX(1)=X
                      XX(2)=Y
                      XX(3)=Z

                      D=SQRT(X*X+Y*Y+Z*Z)
                      ARG=(D*ARGSC)+1
                      NADDR=INT(ARG)

!sh                   IF (NADDR<NPSRNL) THEN
                      IF (D<(NPSRNL-1)/ARGSC+NONLR_S%RSMOOTH(NT)) THEN
                         NONLR_S%CRREXP(IND,NI,ISPINOR __NOACC_omp_arg(idim))=EXP(CITPI*(X*(VKX-QX)+Y*(VKY-QY)+Z*(VKZ-QZ)))*XX(IDIR)
                         IND=IND+1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             ENDIF
!=======================================================================
! end of loop over ions
!=======================================================================
          ENDDO ions
600       NIS = NIS+NONLR_S%NITYP(NT)
       ENDDO type

       ! conjugate phase alteration for spin down: -q/2 -> q/2
       QX=-QX
       QY=-QY
       QZ=-QZ
    ENDDO spinor
!$ ENDDO omp
!$OMP END PARALLEL DO

    RETURN
  END SUBROUTINE PHASERR

!****************** SUBROUTINE PHASER_HF  ******************************
!
! subroutine PHASER_HF
! recalculates the phase factor for the real-space projectors
! for this version the k-point coordinate is explicitly supplied
!
!> @details @ref openmp : Partition each real space projection operator
!> into openmp::omp_nthreads_nonlr_rspace parts. The work on these
!> parts will be distributed over the available OpenMP threads lateron.
!> @n
!> Currently OMP_NUM_THREADS/=1 requires NCORE=1. Consequently
!> the fast index of the real space grid of the non local projection
!> operators is the x-index (GRID\%RL\%NFAST/=3).
!> @n
!> For openmp::omp_nonlr_planewise = .TRUE. (default) the projectors are
!> distributed over (x,y)-planes in a round robin fashion.
!> In case openmp::omp_nonlr_planewise = .FALSE. the projectors are
!> distributed over (z,y)-columns.
!
!***********************************************************************

  SUBROUTINE PHASER_HF(GRID,LATT_CUR,NONLR_S,VK)
    USE lattice
    USE mpimy
    USE constant
    USE pseudo
    USE openmp, ONLY : omp_nthreads_nonlr_rspace,omp_nonlr_planewise
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    REAL(q)            VK(3)
!$  INTEGER idim,ndim

    PROFILING_START('phaser_hf')

#ifdef gammareal
    PROFILING_STOP('phaser_hf')
    RETURN
#else
NOACC !$  ndim=SIZE(NONLR_S%CRREXP,4)
DOACC !$  ndim=1

!$OMP PARALLEL DO NUM_THREADS(omp_nthreads_nonlr_rspace) &
!$OMP PRIVATE(idim,NSPINORS,ISPINORS,NIS,NT,NI,ARGSC,F1,F2,F3,D1,D2,D3,N3LOW,N2LOW,N1LOW, &
!$OMP N3HI,N2HI,N1HI,IND,N2,X2,N2P,N1,X1,N1P,NCOL,N3,X3,X,Y,Z,D,ARG,NADDR,N3P) &
!$OMP SHARED(ndim,GRID,LATT_CUR,NONLR_S,VK,omp_nonlr_planewise)
!$ omp: DO idim=1,ndim
!=======================================================================
! Loop over NSPINORS: here only in case of spin spirals NRSPINOR=2
!=======================================================================
    IF (NONLR_S%LSPIRAL) THEN 
       NSPINORS=2
    ELSE
       NSPINORS=1
    ENDIF

    spinor: DO ISPINOR=1,NSPINORS

!=======================================================================
! loop over all ions
!=======================================================================
       NIS=1

!OCL SCALAR
       type: DO NT=1,NONLR_S%NTYP
          IF (NONLR_S%LMMAX(NT)==0) GOTO 600
          ions: DO NI=NIS,NONLR_S%NITYP(NT)+NIS-1

             IF (NONLR_S%SELECTED_ION>0 .AND. NONLR_S%SELECTED_ION/=NI) CYCLE
             ARGSC=NPSRNL/NONLR_S%PSRMAX(NT)

!=======================================================================
! find lattice points contained within the cutoff-sphere
! this loop might be done in scalar unit
!=======================================================================
             F1=1._q/GRID%NGX
             F2=1._q/GRID%NGY
             F3=1._q/GRID%NGZ
!-----------------------------------------------------------------------
! restrict loop to points contained within a cubus around the ion
!-----------------------------------------------------------------------
!sh
!            D1= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(1)*GRID%NGX
!            D2= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(2)*GRID%NGY
!            D3= NONLR_S%PSRMAX(NT)*LATT_CUR%BNORM(3)*GRID%NGZ
             D1= (NONLR_S%PSRMAX(NT)+NONLR_S%RSMOOTH(NT))*LATT_CUR%BNORM(1)*GRID%NGX
             D2= (NONLR_S%PSRMAX(NT)+NONLR_S%RSMOOTH(NT))*LATT_CUR%BNORM(2)*GRID%NGY
             D3= (NONLR_S%PSRMAX(NT)+NONLR_S%RSMOOTH(NT))*LATT_CUR%BNORM(3)*GRID%NGZ

             N3LOW= INT(NONLR_S%POSION(3,NI)*GRID%NGZ-D3+GRID%NGZ+.99_q)-GRID%NGZ
             N2LOW= INT(NONLR_S%POSION(2,NI)*GRID%NGY-D2+GRID%NGY+.99_q)-GRID%NGY
             N1LOW= INT(NONLR_S%POSION(1,NI)*GRID%NGX-D1+GRID%NGX+.99_q)-GRID%NGX

             N3HI = INT(NONLR_S%POSION(3,NI)*GRID%NGZ+D3)
             N2HI = INT(NONLR_S%POSION(2,NI)*GRID%NGY+D2)
             N1HI = INT(NONLR_S%POSION(1,NI)*GRID%NGX+D1)

!-----------------------------------------------------------------------
! loop over cubus
! MPI version z ist the fast index
!-----------------------------------------------------------------------
             IF (GRID%RL%NFAST==3) THEN
             IND=1

             DO N2=N2LOW,N2HI
                X2=(N2*F2-NONLR_S%POSION(2,NI))
                N2P=MOD(N2+10*GRID%NGY,GRID%NGY)

!!$              IF (MOD(N2P,ndim)/=(idim-1)) CYCLE

                DO N1=N1LOW,N1HI
                   X1=(N1*F1-NONLR_S%POSION(1,NI))
                   N1P=MOD(N1+10*GRID%NGX,GRID%NGX)

                   NCOL=GRID%RL%INDEX(N1P,N2P)
                   IF (NCOL==0) CYCLE
                   IF (GRID%RL%I2(NCOL) /= N1P+1 .OR. GRID%RL%I3(NCOL) /= N2P+1) THEN
                      CALL vtutor%bug("internal ERROR PHASER: " // str(GRID%RL%I2(NCOL)) // " " // &
                         str(N1P+1) // " " // str(GRID%RL%I3(NCOL)) // " " // str(N2P+1), __FILE__, __LINE__)
                   ENDIF

!$                 IF (MOD(NCOL,ndim)/=(idim-1)) CYCLE

                   DO N3=N3LOW,N3HI
                      X3=(N3*F3-NONLR_S%POSION(3,NI))

                      X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
                      Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
                      Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)

                      D=SQRT(X*X+Y*Y+Z*Z)
                      ARG=(D*ARGSC)+1
                      NADDR=INT(ARG)

!sh                   IF (NADDR<NPSRNL) THEN
                      IF (D<(NPSRNL-1)/ARGSC+NONLR_S%RSMOOTH(NT)) THEN
                         NONLR_S%CRREXP(IND,NI,ISPINOR __NOACC_omp_arg(idim))=EXP(CITPI*(X1*VK(1)+X2*VK(2)+X3*VK(3)))
                         IND=IND+1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             ELSE
!-----------------------------------------------------------------------
! loop over cubus around one ion
! conventional version x is fast index
!-----------------------------------------------------------------------
             IND=1
             DO N3=N3LOW,N3HI
                X3=(N3*F3-NONLR_S%POSION(3,NI))
                N3P=MOD(N3+10*GRID%NGZ,GRID%NGZ)

!$              IF (MOD(N3P,ndim)/=(idim-1).AND.omp_nonlr_planewise) CYCLE

                DO N2=N2LOW,N2HI
                   X2=(N2*F2-NONLR_S%POSION(2,NI))
                   N2P=MOD(N2+10*GRID%NGY,GRID%NGY)

                   NCOL=GRID%RL%INDEX(N2P,N3P)

!$                 IF (MOD(NCOL,ndim)/=(idim-1).AND.(.NOT.omp_nonlr_planewise)) CYCLE

                   DO N1=N1LOW,N1HI
                      X1=(N1*F1-NONLR_S%POSION(1,NI))

                      X= X1*LATT_CUR%A(1,1)+X2*LATT_CUR%A(1,2)+X3*LATT_CUR%A(1,3)
                      Y= X1*LATT_CUR%A(2,1)+X2*LATT_CUR%A(2,2)+X3*LATT_CUR%A(2,3)
                      Z= X1*LATT_CUR%A(3,1)+X2*LATT_CUR%A(3,2)+X3*LATT_CUR%A(3,3)

                      D=SQRT(X*X+Y*Y+Z*Z)
                      ARG=(D*ARGSC)+1
                      NADDR=INT(ARG)

!sh                   IF (NADDR<NPSRNL) THEN
                      IF (D<(NPSRNL-1)/ARGSC+NONLR_S%RSMOOTH(NT)) THEN
                         NONLR_S%CRREXP(IND,NI,ISPINOR __NOACC_omp_arg(idim))=EXP(CITPI*(X1*VK(1)+X2*VK(2)+X3*VK(3)))
                         IND=IND+1
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             ENDIF
!=======================================================================
! end of loop over ions
!=======================================================================
          ENDDO ions
600       NIS = NIS+NONLR_S%NITYP(NT)
       ENDDO type

    ENDDO spinor
!$ ENDDO omp
!$OMP END PARALLEL DO
#endif

!$ACC UPDATE DEVICE(NONLR_S%CRREXP) __IF_ASYNC__

    PROFILING_STOP('phaser_hf')

    RETURN
  END SUBROUTINE PHASER_HF

!****************** SUBROUTINE NONLR_SET_SINGLE_ION ********************
!
! generate a non local projector  for a single ion
! reuse the data structure of the full projector whereever that
! is possible
! NOTE: deallocation must go through  NONLR_DEALLOC_SINGLE_ION
! since otherwise the original data structure is destroyed
! the first projector contains the conventional projector, whereas the second
! one is the derivative of the projector
!
!> @details @ref openmp : 
!> Allocate
!>      nonlr_struct::nlimax,
!>      nonlr_struct::rproj,
!> with the additional dimension openmp::omp_nthreads_nonlr_rspace
!> (slowest index).
!
!***********************************************************************

  SUBROUTINE  NONLR_SET_SINGLE_ION(GRID,LATT_CUR, NONLR_S, NONLR_ION, NONLR_IOND, ION, IDIR)
#ifdef _OPENACC
    USE mopenacc
#endif
    USE lattice
    USE mpimy
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT NONE

    TYPE (grid_3d)      GRID
    TYPE (latt)         LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonlr_struct) NONLR_ION  ! descriptor for single ion projector
    TYPE (nonlr_struct) NONLR_IOND ! derivative of projector
    LOGICAL LREALLOCATE
    INTEGER ION, IDIR
    REAL(q)  DISPL(3,NONLR_S%NIONS)
    REAL(q), PARAMETER :: DIS=fd_displacement
#ifdef shmem_rproj
    TYPE(communic), TARGET, SAVE :: DUMMY
    DUMMY%NCPU=1; DUMMY%NODE_ME=1
#endif

    PUSH_ACC_EXEC_ON(.FALSE.)

    NONLR_ION%LREAL=.FALSE.
    NONLR_IOND%LREAL=.FALSE.

    IF (.NOT. NONLR_S%LREAL) THEN
       POP_ACC_EXEC_ON
       RETURN
    ENDIF

    ! copy entire descriptor
!   NONLR_ION=NONLR_S
    CALL NONLR_ASSIGN(NONLR_ION,NONLR_S)
    ! select one specific ion
    NONLR_ION%SELECTED_ION=ION
    ! determine IRALLOC
    NONLR_ION%IRALLOC=-1
    CALL REAL_OPTLAY(GRID,LATT_CUR,NONLR_ION,.TRUE.,LREALLOCATE,-1,-1)

    ! allocate required arrays
    ALLOCATE( &
         NONLR_ION%LMBASE (NONLR_ION%NIONS+1), &
         NONLR_ION%NLIMAX (NONLR_ION%NIONS   __NOACC_omp_arg(omp_nthreads_nonlr_rspace)), &
         NONLR_ION%NLIBASE(NONLR_ION%NIONS+1 __NOACC_omp_arg(omp_nthreads_nonlr_rspace)), &
         NONLR_ION%RPROJ  (NONLR_ION%IRALLOC __NOACC_omp_arg(omp_nthreads_nonlr_rspace)))

    ! copy entire descriptor
!   NONLR_IOND=NONLR_S
    CALL NONLR_ASSIGN(NONLR_IOND,NONLR_S)
    ! select one specific ion
    NONLR_IOND%SELECTED_ION=ION
    ! determine IRALLOC
    NONLR_IOND%IRALLOC=-1
    CALL REAL_OPTLAY(GRID,LATT_CUR,NONLR_IOND,.TRUE.,LREALLOCATE,-1,-1)

    ! allocate required arrays
    ALLOCATE( &
         NONLR_IOND%LMBASE (NONLR_IOND%NIONS+1), &
         NONLR_IOND%NLIMAX (NONLR_IOND%NIONS   __NOACC_omp_arg(omp_nthreads_nonlr_rspace)), &
         NONLR_IOND%NLIBASE(NONLR_IOND%NIONS+1 __NOACC_omp_arg(omp_nthreads_nonlr_rspace)), &
         NONLR_IOND%RPROJ  (NONLR_IOND%IRALLOC __NOACC_omp_arg(omp_nthreads_nonlr_rspace)))

    ! setup single ion projector
    DISPL=0
#ifdef shmem_rproj
    ! the following will effectively switch off the use
    ! of shared memory for the NONLR_ION%RPROJ
    NONLR_ION%COMM_shmem=>DUMMY
#endif
    CALL RSPHER_ALL(GRID,NONLR_ION,LATT_CUR,LATT_CUR,LATT_CUR, DISPL,DISPL,0)

    ! calculate first derivative of projector
    DISPL(IDIR,:)=DIS
#ifdef shmem_rproj
    ! the following will effectively switch off the use
    ! of shared memory for the NONLR_IOND%RPROJ
    NONLR_IOND%COMM_shmem=>DUMMY
#endif
    CALL RSPHER_ALL(GRID,NONLR_IOND,LATT_CUR,LATT_CUR,LATT_CUR, -DISPL,DISPL,1)
    NONLR_IOND%RPROJ= NONLR_IOND%RPROJ*(0.5_q/DIS)

#ifdef _OPENACC
    POP_ACC_EXEC_ON

    CALL ACC_COPYIN_TYPED_VAR(NONLR_ION)
    CALL ACC_COPYIN_TYPED_VAR(NONLR_IOND)
#endif
    RETURN
  END SUBROUTINE NONLR_SET_SINGLE_ION


  SUBROUTINE  NONLR_DEALLOC_SINGLE_ION(NONLR_S)
#ifdef _OPENACC
    USE mopenacc
#endif
    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S

    IF (NONLR_S%LREAL) THEN
#ifdef _OPENACC
       CALL ACC_DELETE_TYPED_VAR(NONLR_S)
#endif
       DEALLOCATE(NONLR_S%LMBASE,NONLR_S%NLIMAX,NONLR_S%NLIBASE,NONLR_S%RPROJ)
    ENDIF

    NULLIFY(NONLR_S%LMBASE,NONLR_S%NLIMAX,NONLR_S%NLIBASE,NONLR_S%NLI,NONLR_S%RPROJ,NONLR_S%CRREXP)
    RETURN
  END SUBROUTINE NONLR_DEALLOC_SINGLE_ION


!****************** SUBROUTINE NONLR_ALLOC_CRREXP **********************
!
! allocate the CRREXP array
!
!> @details @ref openmp : 
!> Allocate
!>      nonlr_struct::crrexp
!> with the additional dimension openmp::omp_nthreads_nonlr_rspace
!> (slowest index).
!
!***********************************************************************

  SUBROUTINE  NONLR_ALLOC_CRREXP(NONLR_S)
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT NONE

    TYPE (nonlr_struct), TARGET :: NONLR_S
    ! local variables
    INTEGER NIONS,IRMAX

    NIONS = NONLR_S%NIONS
    IRMAX = NONLR_S%IRMAX

    NULLIFY(NONLR_S%CRREXP)
    IF (.NOT.NONLR_S%LSPIRAL ) THEN
       ALLOCATE(NONLR_S%CRREXP(IRMAX,NIONS,1 __NOACC_omp_arg(omp_nthreads_nonlr_rspace)))
    ELSE
       ALLOCATE(NONLR_S%CRREXP(IRMAX,NIONS,2 __NOACC_omp_arg(omp_nthreads_nonlr_rspace)))
    ENDIF
    RETURN
  END SUBROUTINE NONLR_ALLOC_CRREXP


!****************** SUBROUTINE NONLR_DEALLOC_CRREXP **********************
!
! deallocate the CRREXP array
!
!***********************************************************************

  SUBROUTINE  NONLR_DEALLOC_CRREXP(NONLR_S)
    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    IF (ASSOCIATED(NONLR_S%CRREXP)) THEN
       DEALLOCATE(NONLR_S%CRREXP)
       NULLIFY(NONLR_S%CRREXP)
    ENDIF
    RETURN
  END SUBROUTINE NONLR_DEALLOC_CRREXP


!****************** SUBROUTINE RPRO1 *********************************
!
! this subroutine calculates the scalar product of one wavefunction with
! all projector functions in real space
! thesis gK Equ. (10.36)
!
!> @details @ref openmp : each real space projection operator is
!> partitioned into openmp::omp_nthreads_nonlr_rspace parts.
!> The work on these parts is done by openmp::omp_nthreads_nonlr_rspace
!> threads. The work is distributed by means of an explicit loop
!> (OMP PARALLEL DO, label omp) and a reduction of CPROJ.
!> @n
!> Alternatively, instead of the loop over the real space parts
!> of the projectors, the loop over ions may be distributed over
!> the available threads.
!> To this end the precompiler flag <TT>avoid_reduction</TT> must
!> be defined.
!> This avoids the reduction of CPROJ.
!> @n@n
!> Under OpenMP the nested loop over \"types\" + \"ions-of-typ\"
!> is replaced by a single loop over \"all ions\".
!
!*********************************************************************
#if !defined(IntelKNL) || !defined(_OPENMP)
  SUBROUTINE RPRO1(NONLR_S, WDES1, W1)
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1

  ! local
    INTEGER :: IP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NIS, NT, LMMAXC, &
         NI, INDMAX, IND, LM
    REAL(q) :: RP, SUMR, SUMI
    COMPLEX(q) CTMP
    REAL(qn),PARAMETER :: ONE=1,ZERO=0
    INTEGER, PARAMETER :: NLM=__align_DIM(101)

!$  INTEGER i,ndim
!$  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM

#ifdef avoidalloc
    REAL(qn) :: WORK(NONLR_S%IRMAX*2),TMP(NLM,2)
    GDEF    :: CPROJ(WDES1%NPRO_TOT)
#else
    REAL(qn),ALLOCATABLE:: WORK(:),TMP(:,:)
    GDEF,ALLOCATABLE   :: CPROJ(:)
#endif
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: WORK,TMP
#endif

    PROFILING_START('rpro1')

#ifndef avoidalloc
    ALLOCATE(CPROJ(WDES1%NPRO_TOT),WORK(NONLR_S%IRMAX*2),TMP(NLM,2))
#endif

    IF (WDES1%NK /= NONLR_S%NK) &
       CALL vtutor%bug("RPRO1: PHASE not properly set up " // str(WDES1%NK) &
          // " " // str(NONLR_S%NK),__FILE__,__LINE__)

    CPROJ=0

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

#ifndef avoid_reduction
!$OMP PARALLEL DO &
!$OMP NUM_THREADS(omp_nthreads_nonlr_rspace) &
!$OMP PRIVATE(i,LMBASE,ISPIRAL,ISPINOR,NLIIND,NIS,NT,LMMAXC,NI,IND,INDMAX,IP,WORK,CTMP,LM,SUMR,SUMI,RP,TMP) &
!$OMP SHARED(NONLR_S,WDES1,W1) &
!$OMP REDUCTION(+:CPROJ)
#else
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(NI,NT,LMMAXC,INDMAX,LMBASE,NLIIND,WORK,IND,IP,CTMP,LM,SUMR,SUMI,RP,TMP)
#endif
!$  omp: DO i=1,ndim
!=======================================================================
! loop over ions
!=======================================================================
    LMBASE= 0

    ISPIRAL = 1
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
#ifndef _OPENMP
       NLIIND= 0
       NIS=1

       typ: DO NT=1,NONLR_S%NTYP
          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) GOTO 600

          ion: DO NI=NIS,NONLR_S%NITYP(NT)+NIS-1
!=======================================================================
!  extract the relevant points for this ion
!=======================================================================
             INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
             IF (INDMAX == 0) GOTO 100
#else
#ifdef avoid_reduction
!$OMP DO
#endif
       ion: DO NI=1,NONLR_S%NIONS
          NT=NONLR_S%ITYP(NI)

          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) CYCLE ion

          INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
          IF (INDMAX==0) CYCLE ion

          LMBASE=NONLR_S%LMBASE(NI)+ISPINOR*NONLR_S%LMBASE(NONLR_S%NIONS+1)
          NLIIND=NONLR_S%NLIBASE(NI __NOACC_omp_arg(i))
#endif
             IF (ASSOCIATED(NONLR_S%CRREXP)) THEN
#ifdef CRREXP_inline
                DO IND=1,INDMAX
                   IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
#ifdef gammareal
                   WORK(IND)      = REAL( W1%CR(IP) ,KIND=qn)
#else
                   CTMP=  W1%CR(IP)*NONLR_S%CRREXP(IND,NI,ISPIRAL __NOACC_omp_arg(i))
                   WORK(IND)      = REAL( CTMP ,KIND=qn)
                   WORK(IND+NONLR_S%IRMAX)=AIMAG(CTMP)
#endif
                ENDDO
#else
                CALL CRREXP_MUL_WAVE( INDMAX, NONLR_S%CRREXP(1,NI,ISPIRAL __NOACC_omp_arg(i)),NONLR_S%NLI(1,NI __NOACC_omp_arg(i)), &
                           W1%CR(ISPINOR*WDES1%GRID%MPLWV+1),WORK(1),WORK(NONLR_S%IRMAX+1))
#endif
             ELSE
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP)
                DO IND=1,INDMAX
                   IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                   WORK(IND)      = REAL( W1%CR(IP) ,KIND=qn)
                ENDDO
             ENDIF
!=======================================================================
! loop over composite indexes L,M
!=======================================================================
#if defined(do_loops) || defined(pro_loop)
             DO LM=1,LMMAXC
                SUMR=0
                SUMI=0
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(RP) REDUCTION(+:SUMR,SUMI)
                DO IND=1,INDMAX
                   RP   =  NONLR_S%RPROJ(IND+ (LM-1)*INDMAX+NLIIND __NOACC_omp_arg(i))
                   SUMR = SUMR + WORK(IND)       * RP
#ifndef gammareal
                   SUMI = SUMI + WORK(IND+NONLR_S%IRMAX) * RP
#endif
                ENDDO
                CPROJ(LM+LMBASE)=GREAL(CMPLX( SUMR , SUMI ,KIND=q) *WDES1%RINPL)
             ENDDO
#else
             CALL DGEMVn( 'T' , INDMAX, LMMAXC, ONE , NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), &
                  INDMAX, WORK(1) , 1 , ZERO ,  TMP(1,1), 1)
#ifndef gammareal
             CALL DGEMVn( 'T' , INDMAX, LMMAXC, ONE , NONLR_S%RPROJ(1+NLIIND  __NOACC_omp_arg(i)), &
                  INDMAX, WORK(1+NONLR_S%IRMAX) , 1 , ZERO ,  TMP(1,2), 1)
#endif

             l_loop: DO LM=1,LMMAXC
#ifdef gammareal
#ifndef _OPENMP
                CPROJ(LM+LMBASE)=TMP(LM,1)*WDES1%RINPL
#else
                CPROJ(LM+LMBASE)=CPROJ(LM+LMBASE)+TMP(LM,1)*WDES1%RINPL
#endif
#else
                SUMR=TMP(LM,1)
                SUMI=TMP(LM,2)
#ifndef _OPENMP
                CPROJ(LM+LMBASE)=GREAL(CMPLX( SUMR , SUMI ,KIND=q) *WDES1%RINPL)
#else
                CPROJ(LM+LMBASE)=CPROJ(LM+LMBASE)+GREAL(CMPLX( SUMR , SUMI ,KIND=q) *WDES1%RINPL)
#endif
#endif
             ENDDO l_loop
#endif
#ifndef _OPENMP
100          LMBASE= LMMAXC+LMBASE
             NLIIND= LMMAXC*INDMAX+NLIIND
          ENDDO ion

600       NIS = NIS+NONLR_S%NITYP(NT)
       ENDDO typ
#else
       ENDDO ion
#ifdef avoid_reduction
!$OMP END DO
#endif
#endif
       IF (NONLR_S%LSPIRAL) ISPIRAL=2
    ENDDO spinor
!$  ENDDO omp
#ifdef avoid_reduction
!$OMP END PARALLEL
#else
!$OMP END PARALLEL DO
#endif

    ! distribute the projected wavefunctions to nodes
    CALL DIS_PROJ(WDES1,CPROJ(1),W1%CPROJ(1))

#ifndef avoidalloc
    DEALLOCATE(CPROJ,WORK,TMP)
#endif

    PROFILING_STOP('rpro1')

    RETURN
  END SUBROUTINE RPRO1
#else
  SUBROUTINE RPRO1(NONLR_S, WDES1, W1)
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1

  ! local
    INTEGER :: IP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NIS, NT, LMMAXC, &
         NI, INDMAX, IND, LM
    REAL(q) :: RP, SUMR, SUMI
    COMPLEX(q) CTMP
!$  INTEGER i,ndim
!$  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM

#ifdef avoidalloc
    REAL(qn) :: WORK(NONLR_S%IRMAX,2)
    GDEF    :: CPROJ(WDES1%NPRO_TOT)
#else
    REAL(qn),ALLOCATABLE:: WORK(:,:),TMP(:,:)
    GDEF,ALLOCATABLE   :: CPROJ(:)
#endif
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: WORK,CPROJ
#endif

    PROFILING_START('rpro1')

#ifndef avoidalloc
    ALLOCATE(CPROJ(WDES1%NPRO_TOT),WORK(NONLR_S%IRMAX,2))
#endif

    IF (WDES1%NK /= NONLR_S%NK) &
       CALL vtutor%bug("RPRO1: PHASE not properly set up " // str(WDES1%NK) &
          // " " // str(NONLR_S%NK),__FILE__,__LINE__)

    CPROJ=0

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

#ifndef avoid_reduction
!$OMP PARALLEL DO &
!$OMP NUM_THREADS(omp_nthreads_nonlr_rspace) &
!$OMP PRIVATE(i,LMBASE,ISPIRAL,ISPINOR,NLIIND,NT,LMMAXC,NI,IND,INDMAX,IP,WORK,CTMP) &
!$OMP SHARED(NONLR_S,WDES1,W1) &
!$OMP REDUCTION(+:CPROJ)
#else
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(NI,NT,LMMAXC,INDMAX,LMBASE,NLIIND,WORK,IND,IP,CTMP)
#endif
!$  omp: DO i=1,ndim
!=======================================================================
! loop over ions
!=======================================================================
    LMBASE= 0

    ISPIRAL = 1
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
#ifdef avoid_reduction
!$OMP DO
#endif
       ion: DO NI=1,NONLR_S%NIONS
          NT=NONLR_S%ITYP(NI)

          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) CYCLE ion

          INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
          IF (INDMAX==0) CYCLE ion

          LMBASE=NONLR_S%LMBASE(NI)+ISPINOR*NONLR_S%LMBASE(NONLR_S%NIONS+1)
          NLIIND=NONLR_S%NLIBASE(NI __NOACC_omp_arg(i))

             IF (ASSOCIATED(NONLR_S%CRREXP)) THEN
#ifdef CRREXP_inline
                DO IND=1,INDMAX
                   IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
#ifdef gammareal
                   WORK(IND,1)= REAL( W1%CR(IP) ,KIND=qn)
#else
                   CTMP=  W1%CR(IP)*NONLR_S%CRREXP(IND,NI,ISPIRAL __NOACC_omp_arg(i))
                   WORK(IND,1)= REAL( CTMP ,KIND=qn)
                   WORK(IND,2)=AIMAG(CTMP)
#endif
                ENDDO
#else
                CALL CRREXP_MUL_WAVE( INDMAX, NONLR_S%CRREXP(1,NI,ISPIRAL __NOACC_omp_arg(i)),NONLR_S%NLI(1,NI __NOACC_omp_arg(i)), &
                           W1%CR(ISPINOR*WDES1%GRID%MPLWV+1),WORK(1,1),WORK(1,2))
#endif
             ELSE
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP)
                DO IND=1,INDMAX
                   IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                   WORK(IND,1)= REAL( W1%CR(IP) ,KIND=qn)
                ENDDO
             ENDIF

             CALL GEMMV_AND_SCATTER(INDMAX,LMMAXC,NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), &
                                    WORK(1,1),WORK(1,2),CPROJ(LMBASE+1),WDES1%RINPL)
       ENDDO ion
#ifdef avoid_reduction
!$OMP END DO
#endif
       IF (NONLR_S%LSPIRAL) ISPIRAL=2
    ENDDO spinor
!$  ENDDO omp
#ifdef avoid_reduction
!$OMP END PARALLEL
#else
!$OMP END PARALLEL DO
#endif

    ! distribute the projected wavefunctions to nodes
    CALL DIS_PROJ(WDES1,CPROJ(1),W1%CPROJ(1))

#ifndef avoidalloc
    DEALLOCATE(CPROJ,WORK)
#endif

    PROFILING_STOP('rpro1')

    RETURN
  END SUBROUTINE RPRO1
#endif
#ifdef _OPENACC
  SUBROUTINE RPRO1_ACC(NONLR_S, WDES1, W1)
    USE openmp, ONLY : omp_nthreads_acc
    IMPLICIT NONE

    TYPE (nonlr_struct) :: NONLR_S
    TYPE (wavedes1)     :: WDES1
    TYPE (wavefun1)     :: W1

  ! local
    INTEGER :: NP, IP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NT, LMMAXC, NI, INDMAX, IND, LM
    GDEF :: CTMP

!$  INTEGER i,ndim

#ifdef avoidalloc
    GDEF :: GWORK(NONLR_S%IRMAX),CPROJ(WDES1%NPRO_TOT)
#else
    GDEF, ALLOCATABLE :: GWORK(:),CPROJ(:)
#endif

    PROFILING_START('rpro1_acc')

#ifndef avoidalloc
    ALLOCATE(GWORK(NONLR_S%IRMAX),CPROJ(WDES1%NPRO_TOT))
#endif
!$ACC ENTER DATA CREATE(CPROJ) __IF_ASYNC__

    IF (WDES1%NK /= NONLR_S%NK) &
       CALL vtutor%bug("RPRO1_ACC: PHASE not properly set up " // str(WDES1%NK) &
          // " " // str(NONLR_S%NK),__FILE__,__LINE__)

!$ACC KERNELS PRESENT(CPROJ) __IF_ASYNC__
    CPROJ=0
!$ACC END KERNELS

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

!$  omp: DO i=1,ndim
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1

       ISPIRAL=1; IF (NONLR_S%LSPIRAL) ISPIRAL=ISPINOR+1

!=======================================================================
! loop over ions
!=======================================================================
!$ACC PARALLEL LOOP GANG PRESENT(CPROJ,WDES1,WDES1%GRID,NONLR_S,W1) &
!$ACC PRIVATE(NT,LMMAXC,INDMAX,LMBASE,NLIIND,GWORK) __IF_ASYNC__
       ion: DO NI=1,NONLR_S%NIONS
          NT=NONLR_S%ITYP(NI)

          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) CYCLE ion

          INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
          IF (INDMAX==0) CYCLE ion

          LMBASE=NONLR_S%LMBASE(NI)+ISPINOR*NONLR_S%LMBASE(NONLR_S%NIONS+1)
          NLIIND=NONLR_S%NLIBASE(NI __NOACC_omp_arg(i))

!$ACC LOOP VECTOR PRIVATE(IP)
          DO IND=1,INDMAX
             IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
#ifdef gammareal
             GWORK(IND)=REAL(W1%CR(IP),q)
#else
             GWORK(IND)=W1%CR(IP)*NONLR_S%CRREXP(IND,NI,ISPIRAL __NOACC_omp_arg(i))
#endif
          ENDDO

!$ACC LOOP PRIVATE(CTMP)
          DO LM=1,LMMAXC
             CTMP=0
!$ACC LOOP VECTOR REDUCTION(+:CTMP)
             DO IND=1,INDMAX
                   CTMP=CTMP+GWORK(IND)*NONLR_S%RPROJ(IND+(LM-1)*INDMAX+NLIIND __NOACC_omp_arg(i))
             ENDDO
             CPROJ(LM+LMBASE)=CPROJ(LM+LMBASE)+CTMP*WDES1%RINPL
          ENDDO

       ENDDO ion

    ENDDO spinor
!$  ENDDO omp

!   PROFILING_START('rpro1_dis')

    ! distribute the projected wavefunctions to nodes
    CALL DIS_PROJ(WDES1,CPROJ(1),W1%CPROJ(1))

!   PROFILING_STOP('rpro1_dis')

!$ACC EXIT DATA DELETE(CPROJ) __IF_ASYNC__
#ifndef avoidalloc
    DEALLOCATE(CPROJ,GWORK)
#endif

    PROFILING_STOP('rpro1_acc')

    RETURN

  END SUBROUTINE RPRO1_ACC
#endif

!****************** SUBROUTINE RPRO1_HF ******************************
!
! this subroutine calculates the scalar product of one wavefunction with
! all projector functions in real space
! thesis gK Equ. (10.36)
! cannot use RPRO1 because of possible type real of GCR
!
!> @details @ref openmp : each real space projection operator is
!> partitioned into openmp::omp_nthreads_nonlr_rspace parts.
!> The work on these parts is done by openmp::omp_nthreads_nonlr_rspace
!> threads. The work is distributed by means of an explicit loop
!> (OMP PARALLEL DO, label omp) and a reduction of CPROJ.
!> @n
!> Alternatively, instead of the loop over the real space parts
!> of the projectors, the loop over ions may be distributed over
!> the available threads. 
!> To this end the precompiler flag <TT>avoid_reduction</TT> must
!> be defined.
!> This avoids the reduction of CPROJ.
!> @n@n
!> Under OpenMP the nested loop over \"types\" + \"ions-of-typ\"
!> is replaced by a single loop over \"all ions\".
!
!*********************************************************************
#if !defined(IntelKNL) || !defined(_OPENMP)
  SUBROUTINE RPRO1_HF(NONLR_S,WDES1,W1,GCR)
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1
    GDEF :: GCR(*)

    INTEGER :: IP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NIS, NT, LMMAXC, &
         NI, INDMAX, IND, LM
    REAL(q) :: RP, SUMR, SUMI
    COMPLEX(q) CTMP 
    REAL(qn),PARAMETER :: ONE=1,ZERO=0
    INTEGER, PARAMETER :: NLM=__align_DIM(101)

!$  INTEGER i,ndim
!$  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM

#ifdef avoidalloc
    REAL(qn) :: WORK(NONLR_S%IRMAX*2),TMP(NLM,2)
    GDEF    :: CPROJ(WDES1%NPRO_TOT)
#else
    REAL(qn),ALLOCATABLE:: WORK(:),TMP(:,:)
    GDEF,ALLOCATABLE   :: CPROJ(:)
#endif
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: WORK,TMP
#endif

    PROFILING_START('rpro1_hf')

#ifndef avoidalloc
    ALLOCATE(CPROJ(WDES1%NPRO_TOT),WORK(NONLR_S%IRMAX*2),TMP(NLM,2))
#endif

    CPROJ=0

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

#ifndef avoid_reduction
!$OMP PARALLEL DO &
!$OMP NUM_THREADS(omp_nthreads_nonlr_rspace) DEFAULT(NONE) &
!$OMP PRIVATE(i,WORK,TMP,LMBASE,ISPIRAL,ISPINOR,NLIIND,NIS,NT,LMMAXC,NI,INDMAX,IND,IP,LM,SUMR,SUMI,RP) &
!$OMP SHARED(NONLR_S,WDES1,W1,GCR) &
!$OMP REDUCTION(+:CPROJ)
#else
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(NI,NT,LMMAXC,INDMAX,LMBASE,NLIIND,WORK,IND,IP,LM,SUMR,SUMI,RP,TMP)
#endif
!$  omp: DO i=1,ndim
!=======================================================================
! loop over ions
!=======================================================================
    LMBASE= 0

    ISPIRAL = 1
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
#ifndef _OPENMP
       NLIIND= 0
       NIS=1

       typ: DO NT=1,NONLR_S%NTYP
          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) GOTO 600

          ion: DO NI=NIS,NONLR_S%NITYP(NT)+NIS-1
!=======================================================================
!  extract the relevant points for this ion
!=======================================================================
             INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
             IF (INDMAX == 0) GOTO 100
#else
#ifdef avoid_reduction
!$OMP DO
#endif
       ion: DO NI=1,NONLR_S%NIONS
          NT=NONLR_S%ITYP(NI)

          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) CYCLE ion

          INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
          IF (INDMAX==0) CYCLE ion

          LMBASE=NONLR_S%LMBASE(NI)+ISPINOR*NONLR_S%LMBASE(NONLR_S%NIONS+1)
          NLIIND=NONLR_S%NLIBASE(NI __NOACC_omp_arg(i))
#endif
             IF (ASSOCIATED(NONLR_S%CRREXP)) THEN
                CALL CRREXP_MUL_GWAVE( INDMAX, NONLR_S%CRREXP(1,NI,ISPIRAL __NOACC_omp_arg(i)),NONLR_S%NLI(1,NI __NOACC_omp_arg(i)), &
                           GCR(ISPINOR*WDES1%GRID%MPLWV+1),WORK(1),WORK(NONLR_S%IRMAX+1))
             ELSE
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP)
                DO IND=1,INDMAX
                   IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                   WORK(IND)      = REAL( GCR(IP) ,KIND=qn)
                ENDDO
             ENDIF
!=======================================================================
! loop over composite indexes L,M
!=======================================================================
#if defined(do_loops) || defined(pro_loop)
             DO LM=1,LMMAXC
                SUMR=0
                SUMI=0
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(RP) REDUCTION(+:SUMR,SUMI)
                DO IND=1,INDMAX
                   RP   =  NONLR_S%RPROJ(IND+ (LM-1)*INDMAX+NLIIND __NOACC_omp_arg(i))
                   SUMR = SUMR + WORK(IND)       * RP
#ifndef gammareal
                   SUMI = SUMI + WORK(IND+NONLR_S%IRMAX) * RP
#endif
                ENDDO
                CPROJ(LM+LMBASE)=GREAL(CMPLX( SUMR , SUMI ,KIND=q) *WDES1%RINPL)
             ENDDO
#else
             CALL DGEMVn( 'T' , INDMAX, LMMAXC, ONE , NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), &
                  INDMAX, WORK(1) , 1 , ZERO ,  TMP(1,1), 1)
#ifndef gammareal
             CALL DGEMVn( 'T' , INDMAX, LMMAXC, ONE , NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), &
                  INDMAX, WORK(1+NONLR_S%IRMAX) , 1 , ZERO ,  TMP(1,2), 1)
#endif

             l_loop: DO LM=1,LMMAXC
#ifdef gammareal
#ifndef _OPENMP
                CPROJ(LM+LMBASE)=TMP(LM,1)*WDES1%RINPL
#else
                CPROJ(LM+LMBASE)=CPROJ(LM+LMBASE)+TMP(LM,1)*WDES1%RINPL
#endif
#else
                SUMR=TMP(LM,1)
                SUMI=TMP(LM,2)
#ifndef _OPENMP
                CPROJ(LM+LMBASE)=GREAL(CMPLX( SUMR , SUMI ,KIND=q) *WDES1%RINPL)
#else
                CPROJ(LM+LMBASE)=CPROJ(LM+LMBASE)+GREAL(CMPLX( SUMR , SUMI ,KIND=q) *WDES1%RINPL)
#endif
#endif
             ENDDO l_loop
#endif
#ifndef _OPENMP
100          LMBASE= LMMAXC+LMBASE
             NLIIND= LMMAXC*INDMAX+NLIIND
          ENDDO ion

600       NIS = NIS+NONLR_S%NITYP(NT)
       ENDDO typ
#else
       ENDDO ion
#ifdef avoid_reduction
!$OMP END DO
#endif
#endif
       IF (NONLR_S%LSPIRAL) ISPIRAL=2
    ENDDO spinor
!$  ENDDO omp
#ifdef avoid_reduction
!$OMP END PARALLEL
#else
!$OMP END PARALLEL DO
#endif

    ! distribute the projected wavefunctions to nodes
    CALL DIS_PROJ(WDES1,CPROJ(1),W1%CPROJ(1))

#ifndef avoidalloc
    DEALLOCATE(CPROJ,WORK,TMP)
#endif

    PROFILING_STOP('rpro1_hf')

    RETURN
  END SUBROUTINE RPRO1_HF
#else
  SUBROUTINE RPRO1_HF(NONLR_S,WDES1,W1,GCR)
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1
    GDEF :: GCR(*)

    INTEGER :: IP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NT, LMMAXC, &
         NI, INDMAX, IND

!$  INTEGER i,ndim
!$  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM

#ifdef avoidalloc
    REAL(qn) :: WORK(NONLR_S%IRMAX,2)
    GDEF     :: CPROJ(WDES1%NPRO_TOT)
#else
    REAL(qn),ALLOCATABLE:: WORK(:,:)
    GDEF,ALLOCATABLE    :: CPROJ(:)
#endif
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: WORK,CPROJ
#endif

    PROFILING_START('rpro1_hf')

#ifndef avoidalloc
    ALLOCATE(CPROJ(WDES1%NPRO_TOT),WORK(NONLR_S%IRMAX,2))
#endif

    CPROJ=0

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

#ifndef avoid_reduction
!$OMP PARALLEL DO &
!$OMP NUM_THREADS(omp_nthreads_nonlr_rspace) DEFAULT(NONE) &
!$OMP PRIVATE(i,WORK,LMBASE,ISPIRAL,ISPINOR,NLIIND,NT,LMMAXC,NI,INDMAX,IND,IP) &
!$OMP SHARED(NONLR_S,WDES1,W1,GCR) &
!$OMP REDUCTION(+:CPROJ)
#else
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(NI,NT,LMMAXC,INDMAX,LMBASE,NLIIND,WORK,IND,IP)
#endif
!$  omp: DO i=1,ndim
!=======================================================================
! loop over ions
!=======================================================================
    LMBASE= 0

    ISPIRAL = 1
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
#ifdef avoid_reduction
!$OMP DO
#endif
       ion: DO NI=1,NONLR_S%NIONS
          NT=NONLR_S%ITYP(NI)

          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) CYCLE ion

          INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
          IF (INDMAX==0) CYCLE ion

          LMBASE=NONLR_S%LMBASE(NI)+ISPINOR*NONLR_S%LMBASE(NONLR_S%NIONS+1)
          NLIIND=NONLR_S%NLIBASE(NI __NOACC_omp_arg(i))

             IF (ASSOCIATED(NONLR_S%CRREXP)) THEN
                CALL CRREXP_MUL_GWAVE( INDMAX, NONLR_S%CRREXP(1,NI,ISPIRAL __NOACC_omp_arg(i)),NONLR_S%NLI(1,NI __NOACC_omp_arg(i)), &
                           GCR(ISPINOR*WDES1%GRID%MPLWV+1),WORK(1,1),WORK(1,2))
             ELSE
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP)
                DO IND=1,INDMAX
                   IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                   WORK(IND,1) = REAL( GCR(IP) ,KIND=qn)
                ENDDO
             ENDIF

             CALL GEMMV_AND_SCATTER(INDMAX,LMMAXC,NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), &
                                    WORK(1,1),WORK(1,2),CPROJ(LMBASE+1),WDES1%RINPL)
       ENDDO ion
#ifdef avoid_reduction
!$OMP END DO
#endif
       IF (NONLR_S%LSPIRAL) ISPIRAL=2
    ENDDO spinor
!$  ENDDO omp
#ifdef avoid_reduction
!$OMP END PARALLEL
#else
!$OMP END PARALLEL DO
#endif

    ! distribute the projected wavefunctions to nodes
    CALL DIS_PROJ(WDES1,CPROJ(1),W1%CPROJ(1))

#ifndef avoidalloc
    DEALLOCATE(CPROJ,WORK)
#endif

    PROFILING_STOP('rpro1_hf')

    RETURN
  END SUBROUTINE RPRO1_HF
#endif


!****************** SUBROUTINE RPRO1_REAL ****************************
!
! this subroutine calculates the scalar product of one wavefunction with
! all projector functions in real space
! thesis gK Equ. (10.36)
! this routine is essentially identical to RPRO1 above, but GCR
! is defined as REAL in this subroutine
!
! for details on OMP see previous routine
!
!*********************************************************************

#if !defined(IntelKNL) || !defined(_OPENMP)
  SUBROUTINE RPRO1_REAL(NONLR_S,WDES1,W1,GCR)
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1
    REAL(q) :: GCR(*)

    INTEGER :: IP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NIS, NT, LMMAXC, &
         NI, INDMAX, IND, LM
    REAL(q) :: RP, SUMR, SUMI
    COMPLEX(q) CTMP 
    REAL(qn),PARAMETER :: ONE=1,ZERO=0
    INTEGER, PARAMETER :: NLM=__align_DIM(101)

!$  INTEGER i,ndim
!$  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM

#ifdef avoidalloc
    REAL(qn) :: WORK(NONLR_S%IRMAX*2),TMP(NLM,2)
    GDEF    :: CPROJ(WDES1%NPRO_TOT)
#else
    REAL(qn),ALLOCATABLE:: WORK(:),TMP(:,:)
    GDEF,ALLOCATABLE   :: CPROJ(:)
#endif
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: WORK,TMP
#endif

    PROFILING_START('rpro1_hf')

#ifndef avoidalloc
    ALLOCATE(CPROJ(WDES1%NPRO_TOT),WORK(NONLR_S%IRMAX*2),TMP(NLM,2))
#endif

    CPROJ=0

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

#ifndef avoid_reduction
!$OMP PARALLEL DO &
!$OMP NUM_THREADS(omp_nthreads_nonlr_rspace) DEFAULT(NONE) &
!$OMP PRIVATE(i,WORK,TMP,LMBASE,ISPIRAL,ISPINOR,NLIIND,NIS,NT,LMMAXC,NI,INDMAX,IND,IP,LM,SUMR,SUMI,RP) &
!$OMP SHARED(NONLR_S,WDES1,W1,GCR) &
!$OMP REDUCTION(+:CPROJ)
#else
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(NI,NT,LMMAXC,INDMAX,LMBASE,NLIIND,IP,WORK,IND,LM,SUMR,SUMI,RP,TMP)
#endif
!$  omp: DO i=1,ndim
!=======================================================================
! loop over ions
!=======================================================================
    LMBASE= 0

    ISPIRAL = 1
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
#ifndef _OPENMP
       NLIIND= 0
       NIS=1

       typ: DO NT=1,NONLR_S%NTYP
          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) GOTO 600

          ion: DO NI=NIS,NONLR_S%NITYP(NT)+NIS-1
!=======================================================================
!  extract the relevant points for this ion
!=======================================================================
             INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
             IF (INDMAX == 0) GOTO 100
#else
#ifdef avoid_reduction
!$OMP DO
#endif
       ion: DO NI=1,NONLR_S%NIONS
          NT=NONLR_S%ITYP(NI)

          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) CYCLE ion

          INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
          IF (INDMAX==0) CYCLE ion

          LMBASE=NONLR_S%LMBASE(NI)+ISPINOR*NONLR_S%LMBASE(NONLR_S%NIONS+1)
          NLIIND=NONLR_S%NLIBASE(NI __NOACC_omp_arg(i))
#endif
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP)
          DO IND=1,INDMAX
             IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
             WORK(IND)      = REAL( GCR(IP) ,KIND=qn)
          ENDDO
!=======================================================================
! loop over composite indexes L,M
!=======================================================================
#if defined(do_loops) || defined(pro_loop)
             DO LM=1,LMMAXC
                SUMR=0
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(RP) REDUCTION(+:SUMR,SUMI)
                DO IND=1,INDMAX
                   RP   =  NONLR_S%RPROJ(IND+ (LM-1)*INDMAX+NLIIND __NOACC_omp_arg(i))
                   SUMR = SUMR + WORK(IND)       * RP
                ENDDO
                CPROJ(LM+LMBASE)=SUMR*WDES1%RINPL
             ENDDO
#else
             CALL DGEMVn( 'T' , INDMAX, LMMAXC, ONE , NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), &
                  INDMAX, WORK(1) , 1 , ZERO ,  TMP(1,1), 1)

             l_loop: DO LM=1,LMMAXC
#ifndef _OPENMP
                CPROJ(LM+LMBASE)=TMP(LM,1)*WDES1%RINPL
#else
                CPROJ(LM+LMBASE)=CPROJ(LM+LMBASE)+TMP(LM,1)*WDES1%RINPL
#endif
             ENDDO l_loop
#endif
#ifndef _OPENMP
100          LMBASE= LMMAXC+LMBASE
             NLIIND= LMMAXC*INDMAX+NLIIND
          ENDDO ion

600       NIS = NIS+NONLR_S%NITYP(NT)
       ENDDO typ
#else
       ENDDO ion
#ifdef avoid_reduction
!$OMP END DO
#endif
#endif
       IF (NONLR_S%LSPIRAL) ISPIRAL=2
    ENDDO spinor
!$  ENDDO omp
#ifdef avoid_reduction
!$OMP END PARALLEL
#else
!$OMP END PARALLEL DO
#endif

    ! distribute the projected wavefunctions to nodes
    CALL DIS_PROJ(WDES1,CPROJ(1),W1%CPROJ(1))

#ifndef avoidalloc
    DEALLOCATE(CPROJ,WORK,TMP)
#endif

    PROFILING_STOP('rpro1_hf')

    RETURN
  END SUBROUTINE RPRO1_REAL
#else
  SUBROUTINE RPRO1_REAL(NONLR_S,WDES1,W1,GCR)
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1
    REAL(q) :: GCR(*)

    INTEGER :: IP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NT, LMMAXC, &
         NI, INDMAX, IND

!$  INTEGER i,ndim
!$  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM

#ifdef avoidalloc
    REAL(qn) :: WORK(NONLR_S%IRMAX,2)
    GDEF     :: CPROJ(WDES1%NPRO_TOT)
#else
    REAL(qn),ALLOCATABLE:: WORK(:,:)
    GDEF,ALLOCATABLE    :: CPROJ(:)
#endif
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: WORK,CPROJ
#endif

    PROFILING_START('rpro1_hf')

#ifndef avoidalloc
    ALLOCATE(CPROJ(WDES1%NPRO_TOT),WORK(NONLR_S%IRMAX,2))
#endif

    CPROJ=0

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

#ifndef avoid_reduction
!$OMP PARALLEL DO &
!$OMP NUM_THREADS(omp_nthreads_nonlr_rspace) DEFAULT(NONE) &
!$OMP PRIVATE(i,WORK,LMBASE,ISPIRAL,ISPINOR,NLIIND,NT,LMMAXC,NI,INDMAX,IND,IP) &
!$OMP SHARED(NONLR_S,WDES1,W1,GCR) &
!$OMP REDUCTION(+:CPROJ)
#else
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(NI,NT,LMMAXC,INDMAX,LMBASE,NLIIND,WORK,IND,IP)
#endif
!$  omp: DO i=1,ndim
!=======================================================================
! loop over ions
!=======================================================================
    LMBASE= 0

    ISPIRAL = 1
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
#ifdef avoid_reduction
!$OMP DO
#endif
       ion: DO NI=1,NONLR_S%NIONS
          NT=NONLR_S%ITYP(NI)

          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) CYCLE ion

          INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
          IF (INDMAX==0) CYCLE ion

          LMBASE=NONLR_S%LMBASE(NI)+ISPINOR*NONLR_S%LMBASE(NONLR_S%NIONS+1)
          NLIIND=NONLR_S%NLIBASE(NI __NOACC_omp_arg(i))

!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP)
          DO IND=1,INDMAX
             IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
             WORK(IND,1) = REAL( GCR(IP) ,KIND=qn)
          ENDDO
          CALL DGEMMV_AND_SCATTER(INDMAX,LMMAXC,NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), &
               WORK(1,1),CPROJ(LMBASE+1),WDES1%RINPL)
       ENDDO ion
#ifdef avoid_reduction
!$OMP END DO
#endif
       IF (NONLR_S%LSPIRAL) ISPIRAL=2
    ENDDO spinor
!$  ENDDO omp
#ifdef avoid_reduction
!$OMP END PARALLEL
#else
!$OMP END PARALLEL DO
#endif
    ! distribute the projected wavefunctions to nodes
    CALL DIS_PROJ(WDES1,CPROJ(1),W1%CPROJ(1))

#ifndef avoidalloc
    DEALLOCATE(CPROJ,WORK)
#endif

    PROFILING_STOP('rpro1_hf')

    RETURN
  END SUBROUTINE RPRO1_REAL
#endif

!****************** SUBROUTINE RPROMU ********************************
!
!  this subroutine  calculates the projection of a set of
!  bands onto the
!  real space projection operators
!
!> @details @ref openmp : each real space projection operator is
!> partitioned into openmp::omp_nthreads_nonlr_rspace parts.
!> The work on these parts is done by openmp::omp_nthreads_nonlr_rspace
!> threads. The work is distributed by means of an explicit loop
!> (OMP PARALLEL DO, label omp) and a reduction of CPROJ.
!> @n
!> Alternatively, instead of the loop over the real space parts
!> of the projectors, the loop over ions may be distributed over
!> the available threads.
!> To this end the precompiler flag <TT>avoid_reduction</TT> must
!> be defined.
!> This avoids the reduction of CPROJ.
!> @n@n
!> Under OpenMP the nested loop over \"types\" + \"ions-of-typ\"
!> is replaced by a single loop over \"all ions\".
!
!*********************************************************************
#if !defined(IntelKNL) || !defined(_OPENMP)
  SUBROUTINE RPROMU(NONLR_S, WDES1, W1, NSIM, LDO)
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT NONE

    INTEGER NSIM
    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1(NSIM)
    LOGICAL            LDO(NSIM)
  ! local
    INTEGER :: NP, IP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NIS, NT, LMMAXC, &
         NI, INDMAX, IND, IND0, NPFILL, LM
    REAL(q) :: SUMR, SUMI
    COMPLEX(q) CTMP 
    REAL(qn),PARAMETER :: ONE=1,ZERO=0
    INTEGER, PARAMETER :: NLM=__align_DIM(101)

!$  INTEGER i,ndim
!$  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM

#ifdef avoidalloc
    REAL(qn) :: WORK(ndata*NONLR_S%IRMAX*NSIM),TMP(NLM, ndata*NSIM)
    GDEF    :: CPROJ(WDES1%NPRO_TOT,NSIM)
#else
    REAL(qn),ALLOCATABLE:: WORK(:),TMP(:,:)
    GDEF,ALLOCATABLE   :: CPROJ(:,:)
#endif
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: WORK,TMP
#endif

    PROFILING_START('rpromu')

#ifndef avoidalloc
    ALLOCATE(CPROJ(WDES1%NPRO_TOT,NSIM),WORK(ndata*NONLR_S%IRMAX*NSIM),TMP(NLM, ndata*NSIM))
#endif
    IF (WDES1%NK /= NONLR_S%NK) &
       CALL vtutor%bug("RPROMU: PHASE not properly set up " // str(WDES1%NK) &
          // " " // str(NONLR_S%NK),__FILE__,__LINE__)

    CPROJ=0

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

#ifndef avoid_reduction
!$OMP PARALLEL DO &
!$OMP NUM_THREADS(omp_nthreads_nonlr_rspace) &
!$OMP PRIVATE(i,WORK,TMP,LMBASE,ISPIRAL,ISPINOR,NLIIND,NIS,NT,LMMAXC,NI,INDMAX,IND0,NPFILL,NP,IND,IP,CTMP,LM,SUMR,SUMI) &
!$OMP SHARED(NONLR_S,WDES1,W1,NSIM,LDO) &
!$OMP REDUCTION(+:CPROJ)
#else
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(NI,NT,LMMAXC,INDMAX,LMBASE,NLIIND,IND0,NPFILL,NP,WORK,IND,IP,CTMP,LM,SUMR,SUMI,TMP)
#endif
!$  omp: DO i=1,ndim
!=======================================================================
! loop over ions
!=======================================================================
    LMBASE= 0

    ISPIRAL = 1
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
#ifndef _OPENMP
       NLIIND= 0
       NIS=1

       typ: DO NT=1,NONLR_S%NTYP
          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) GOTO 600

          ion: DO NI=NIS,NONLR_S%NITYP(NT)+NIS-1
!=======================================================================
!  extract the relevant points for this ion
!=======================================================================
             INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
             IF (INDMAX == 0) GOTO 100
#else
#ifdef avoid_reduction
!$OMP DO
#endif
       ion: DO NI=1,NONLR_S%NIONS
          NT=NONLR_S%ITYP(NI)

          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) CYCLE ion

          INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
          IF (INDMAX==0) CYCLE ion

          LMBASE=NONLR_S%LMBASE(NI)+ISPINOR*NONLR_S%LMBASE(NONLR_S%NIONS+1)
          NLIIND=NONLR_S%NLIBASE(NI __NOACC_omp_arg(i))
#endif
             IND0=0
             NPFILL=0

             DO NP=1,NSIM

                IF (LDO(NP)) THEN
                   IF (ASSOCIATED(NONLR_S%CRREXP)) THEN
#ifdef CRREXP_inline
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP,CTMP)
                      DO IND=1,INDMAX
                         IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                         CTMP=  W1(NP)%CR(IP)*NONLR_S%CRREXP(IND,NI,ISPIRAL __NOACC_omp_arg(i))
                         WORK(IND+IND0) = REAL( CTMP ,KIND=q)
                         WORK(IND+(NONLR_S%IRMAX+IND0))=AIMAG(CTMP)
                      ENDDO
#else
                      CALL CRREXP_MUL_WAVE( INDMAX, NONLR_S%CRREXP(1,NI,ISPIRAL __NOACC_omp_arg(i)), NONLR_S%NLI(1,NI __NOACC_omp_arg(i)), &
                           W1(NP)%CR(ISPINOR*WDES1%GRID%MPLWV+1),WORK(IND0+1),WORK(IND0+(NONLR_S%IRMAX+1)))
#endif
                   ELSE
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP)
                      DO IND=1,INDMAX
                         IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                         WORK(IND+IND0) = REAL( W1(NP)%CR(IP) ,KIND=qn)
                      ENDDO
                   ENDIF

                   IND0=IND0+ndata * NONLR_S%IRMAX
                   NPFILL=NPFILL+1
                ENDIF

             ENDDO

!=======================================================================
! loop over composite indexes L,M
!=======================================================================
             PROFILING_START('rpromu_gemm')
#ifdef RPROMU_DGEMV
             DO NP=0,ndata*NPFILL-1
                CALL DGEMVn( 'T' , INDMAX, LMMAXC, ONE , NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), &
                     INDMAX, WORK(1+NONLR_S%IRMAX*NP) , 1 , ZERO ,  TMP(1,NP+1), 1)
             ENDDO
#else
             CALL DGEMMn( 'T', 'N' , LMMAXC,  ndata*NPFILL, INDMAX, ONE, &
                  NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), INDMAX, WORK(1), NONLR_S%IRMAX, &
                  ZERO,  TMP(1,1), NLM )
#endif
             PROFILING_STOP('rpromu_gemm')

             IND0=0
             DO NP=1,NSIM
                IF (LDO(NP)) THEN
                   l_loop: DO LM=1,LMMAXC
#ifdef gammareal
#ifndef _OPENMP
                      CPROJ(LM+LMBASE,NP)=TMP(LM,1+IND0)*WDES1%RINPL
#else
                      CPROJ(LM+LMBASE,NP)=CPROJ(LM+LMBASE,NP)+TMP(LM,1+IND0)*WDES1%RINPL
#endif
#else
                      SUMR=TMP(LM,1+IND0)
                      SUMI=TMP(LM,2+IND0)
#ifndef _OPENMP
                      CPROJ(LM+LMBASE,NP)=GREAL(CMPLX( SUMR , SUMI ,KIND=q) *WDES1%RINPL)
#else
                      CPROJ(LM+LMBASE,NP)=CPROJ(LM+LMBASE,NP)+GREAL(CMPLX( SUMR , SUMI ,KIND=q) *WDES1%RINPL)
#endif
#endif
                   ENDDO l_loop
                   IND0=IND0+ndata
                ENDIF
             ENDDO
#ifndef _OPENMP
100          LMBASE= LMMAXC+LMBASE
             NLIIND= LMMAXC*INDMAX+NLIIND
          ENDDO ion


600       NIS = NIS+NONLR_S%NITYP(NT)
       ENDDO typ
#else
       ENDDO ion
#ifdef avoid_reduction
!$OMP END DO
#endif
#endif
       IF (NONLR_S%LSPIRAL) ISPIRAL=2
    ENDDO spinor
!$  ENDDO omp
#ifdef avoid_reduction
!$OMP END PARALLEL
#else
!$OMP END PARALLEL DO
#endif

    ! distribute the projected wavefunctions to nodes
    DO NP=1,NSIM
       IF (LDO(NP)) THEN
          CALL DIS_PROJ(WDES1,CPROJ(1,NP),W1(NP)%CPROJ(1))
       ENDIF
    ENDDO

#ifndef avoidalloc
    DEALLOCATE(CPROJ,WORK,TMP)
#endif

    PROFILING_STOP('rpromu')

    RETURN
  END SUBROUTINE RPROMU
#else
  SUBROUTINE RPROMU(NONLR_S, WDES1, W1, NSIM, LDO)
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT NONE

    INTEGER NSIM
    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)     WDES1
    TYPE (wavefun1)     W1(NSIM)
    LOGICAL             LDO(NSIM)
  ! local
    INTEGER :: NP, IP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NT, LMMAXC, &
         NI, INDMAX, IND
    COMPLEX(q) CTMP

!$  INTEGER i,ndim
!$  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM

#ifdef avoidalloc
    REAL(qn) :: WORK(NONLR_S%IRMAX,2)
    GDEF     :: CPROJ(WDES1%NPRO_TOT,NSIM)
#else
    REAL(qn),ALLOCATABLE:: WORK(:,:)
    GDEF,ALLOCATABLE    :: CPROJ(:,:)
#endif
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: WORK
#endif

    PROFILING_START('rpromu')

#ifndef avoidalloc
    ALLOCATE(CPROJ(WDES1%NPRO_TOT,NSIM),WORK(NONLR_S%IRMAX,2))
#endif
    IF (WDES1%NK /= NONLR_S%NK) &
       CALL vtutor%bug("RPROMU: PHASE not properly set up " // str(WDES1%NK) &
          // " " // str(NONLR_S%NK),__FILE__,__LINE__)

    CPROJ=0

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

#ifndef avoid_reduction
!$OMP PARALLEL DO &
!$OMP NUM_THREADS(omp_nthreads_nonlr_rspace) &
!$OMP PRIVATE(i,WORK,LMBASE,ISPIRAL,ISPINOR,NLIIND,NT,LMMAXC,NI,INDMAX,NP,IND,IP,CTMP) &
!$OMP SHARED(NONLR_S,WDES1,W1,NSIM,LDO) &
!$OMP REDUCTION(+:CPROJ)
#else
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(NI,NT,LMMAXC,INDMAX,LMBASE,NLIIND,NP,WORK,IND,IP,CTMP)
#endif
!$  omp: DO i=1,ndim
!=======================================================================
! loop over ions
!=======================================================================
    LMBASE= 0

    ISPIRAL = 1
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
#ifdef avoid_reduction
!$OMP DO
#endif
       ion: DO NI=1,NONLR_S%NIONS
          NT=NONLR_S%ITYP(NI)

          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) CYCLE ion

          INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
          IF (INDMAX==0) CYCLE ion

          LMBASE=NONLR_S%LMBASE(NI)+ISPINOR*NONLR_S%LMBASE(NONLR_S%NIONS+1)
          NLIIND=NONLR_S%NLIBASE(NI __NOACC_omp_arg(i))

             DO NP=1,NSIM
                IF (.NOT.LDO(NP)) CYCLE

!               PROFILING_START('rpromu_gather')

                IF (ASSOCIATED(NONLR_S%CRREXP)) THEN
#ifdef CRREXP_inline
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP,CTMP)
                   DO IND=1,INDMAX
                      IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                      CTMP=  W1(NP)%CR(IP)*NONLR_S%CRREXP(IND,NI,ISPIRAL __NOACC_omp_arg(i))
                      WORK(IND,1)= REAL(CTMP,KIND=q)
                      WORK(IND,2)=AIMAG(CTMP)
                   ENDDO
#else
                   CALL CRREXP_MUL_WAVE( INDMAX, NONLR_S%CRREXP(1,NI,ISPIRAL __NOACC_omp_arg(i)),NONLR_S%NLI(1,NI __NOACC_omp_arg(i)), &
                        W1(NP)%CR(ISPINOR*WDES1%GRID%MPLWV+1),WORK(1,1),WORK(1,2))
#endif
                ELSE
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP)
                   DO IND=1,INDMAX
                      IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                      WORK(IND,1) = REAL( W1(NP)%CR(IP) ,KIND=qn)
                   ENDDO
                ENDIF

!               PROFILING_STOP('rpromu_gather')
!               PROFILING_START('rpromu_gemm')

                CALL GEMMV_AND_SCATTER(INDMAX,LMMAXC,NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), &
                                       WORK(1,1),WORK(1,2),CPROJ(LMBASE+1,NP),WDES1%RINPL)

!               PROFILING_STOP('rpromu_gemm')
             ENDDO
       ENDDO ion
#ifdef avoid_reduction
!$OMP END DO
#endif
       IF (NONLR_S%LSPIRAL) ISPIRAL=2
    ENDDO spinor
!$  ENDDO omp
#ifdef avoid_reduction
!$OMP END PARALLEL
#else
!$OMP END PARALLEL DO
#endif

!   PROFILING_START('rpromu_dis')

    ! distribute the projected wavefunctions to nodes
    DO NP=1,NSIM
       IF (LDO(NP)) THEN
          CALL DIS_PROJ(WDES1,CPROJ(1,NP),W1(NP)%CPROJ(1))
       ENDIF
    ENDDO

!   PROFILING_STOP('rpromu_dis')

#ifndef avoidalloc
    DEALLOCATE(CPROJ,WORK)
#endif

    PROFILING_STOP('rpromu')

    RETURN
  END SUBROUTINE RPROMU
#endif
#ifdef _OPENACC
  SUBROUTINE RPROMU_ACC(NONLR_S, WDES1, W1, NSIM, LDO)
    USE openmp, ONLY : omp_nthreads_acc
    IMPLICIT NONE

    TYPE (nonlr_struct) :: NONLR_S
    TYPE (wavedes1)     :: WDES1
    TYPE (wavefun1)     :: W1(NSIM)
    INTEGER             :: NSIM
    LOGICAL             :: LDO(NSIM)

  ! local
    INTEGER :: NP, IP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NT, LMMAXC, NI, INDMAX, IND, LM, LMP
    GDEF :: CTMP

    INTEGER, PARAMETER :: BLOCKSIZE=512
    INTEGER :: IBLOCK

!$ACC ROUTINE(SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX) SEQ

!$  INTEGER i,ndim

#ifdef avoidalloc
    GDEF :: GWORK(BLOCKSIZE),CPROJ(WDES1%NPRO_TOT,NSIM)
#else
    GDEF, ALLOCATABLE :: GWORK(:),CPROJ(:,:)
#endif

    PROFILING_START('rpromu_acc')

#ifndef avoidalloc
    ALLOCATE(GWORK(BLOCKSIZE),CPROJ(WDES1%NPRO_TOT,NSIM))
#endif
!$ACC ENTER DATA CREATE(CPROJ) __IF_ASYNC__

    IF (WDES1%NK /= NONLR_S%NK) &
       CALL vtutor%bug("RPROMU: PHASE not properly set up " // str(WDES1%NK) &
          // " " // str(NONLR_S%NK),__FILE__,__LINE__)

!$ACC KERNELS PRESENT(CPROJ) __IF_ASYNC__
    CPROJ=0
!$ACC END KERNELS

!$ACC ENTER DATA COPYIN(LDO) __IF_ASYNC__

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

!$  omp: DO i=1,ndim
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1

       ISPIRAL=1; IF (NONLR_S%LSPIRAL) ISPIRAL=ISPINOR+1

!=======================================================================
! loop over ions
!=======================================================================
!$ACC PARALLEL LOOP COLLAPSE(3) GANG PRESENT(CPROJ,WDES1,WDES1%GRID,NONLR_S,W1,LDO) &
!$ACC PRIVATE(NT,LMMAXC,INDMAX,LMBASE,NLIIND,IBLOCK,GWORK) __IF_ASYNC__
       block: DO IBLOCK=0,NONLR_S%IRMAX-1,BLOCKSIZE
       ion: DO NI=1,NONLR_S%NIONS
          nsim: DO NP=1,NSIM

             NT=NONLR_S%ITYP(NI); LMMAXC=NONLR_S%LMMAX(NT)
             INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))

             IF (LMMAXC==0 .OR. INDMAX==0 .OR. .NOT.LDO(NP)) CYCLE

             LMBASE=NONLR_S%LMBASE(NI)+ISPINOR*NONLR_S%LMBASE(NONLR_S%NIONS+1)
             NLIIND=NONLR_S%NLIBASE(NI __NOACC_omp_arg(i))

!$ACC CACHE(GWORK(1:BLOCKSIZE))
!$ACC LOOP VECTOR PRIVATE(IP)
             DO IND=IBLOCK+1,MIN(IBLOCK+BLOCKSIZE,INDMAX)
                IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
#ifdef gammareal
                GWORK(IND-IBLOCK)=REAL(W1(NP)%CR(IP),q)
#else
                GWORK(IND-IBLOCK)=W1(NP)%CR(IP)*NONLR_S%CRREXP(IND,NI,ISPIRAL __NOACC_omp_arg(i))
#endif
             ENDDO

!$ACC LOOP SEQ PRIVATE(CTMP,LMP)
             DO LM=1,LMMAXC
                CTMP=0
                LMP=(LM-1)*INDMAX+NLIIND
!$ACC LOOP VECTOR REDUCTION(+:CTMP)
                DO IND=IBLOCK+1,MIN(IBLOCK+BLOCKSIZE,INDMAX)
                   CTMP=CTMP+GWORK(IND-IBLOCK)*NONLR_S%RPROJ(IND+LMP __NOACC_omp_arg(i))
                ENDDO
#ifdef gammareal
!$ACC ATOMIC UPDATE
                CPROJ(LM+LMBASE,NP)=CPROJ(LM+LMBASE,NP)+CTMP*WDES1%RINPL
#else
                CALL SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX(CPROJ(LM+LMBASE,NP),CTMP*WDES1%RINPL)
#endif
             ENDDO


          ENDDO nsim
       ENDDO ion
       ENDDO block

    ENDDO spinor
!$  ENDDO omp

!   PROFILING_START('rpromu_dis')

    ! distribute the projected wavefunctions to nodes
!!!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(NP) &
!!!$OMP NUM_THREADS(omp_nthreads_acc)
    DO NP=1,NSIM
       IF (LDO(NP)) CALL DIS_PROJ(WDES1,CPROJ(1,NP),W1(NP)%CPROJ(1))
    ENDDO
!!!$OMP END PARALLEL DO

!   PROFILING_STOP('rpromu_dis')

!$ACC EXIT DATA DELETE(LDO,CPROJ) __IF_ASYNC__
#ifndef avoidalloc
    DEALLOCATE(CPROJ,GWORK)
#endif

    PROFILING_STOP('rpromu_acc')

    RETURN
  END SUBROUTINE RPROMU_ACC
#endif

!****************** SUBROUTINE RPROMU_HF *****************************
!
! essentially a copy of the previous routine 
! -  with W1(NP)%CR(i) -> GCR(i,NP)
! -  phase factor check removed
! -  LDO array removed
!
!> @details @ref openmp : each real space projection operator is
!> partitioned into openmp::omp_nthreads_nonlr_rspace parts.
!> The work on these parts is done by openmp::omp_nthreads_nonlr_rspace
!> threads. The work is distributed by means of an explicit loop
!> (OMP PARALLEL DO, label omp) and a reduction of CPROJ.
!> @n
!> Alternatively, instead of the loop over the real space parts
!> of the projectors, the loop over ions may be distributed over
!> the available threads.
!> To this end the precompiler flag <TT>avoid_reduction</TT> must
!> be defined.
!> This avoids the reduction of CPROJ.
!> @n@n
!> Under OpenMP the nested loop over \"types\" + \"ions-of-typ\"
!> is replaced by a single loop over \"all ions\".
!
!*********************************************************************
#if !defined(IntelKNL) || !defined(_OPENMP)
  SUBROUTINE RPROMU_HF(NONLR_S, WDES1, W1, NSIM, GCR, LD)
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT NONE

    INTEGER NSIM
    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1(NSIM)
    INTEGER :: LD  ! leading dimension of GCR
    GDEF :: GCR(LD, NSIM)
  ! local
    INTEGER :: NP, IP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NIS, NT, LMMAXC, &
         NI, INDMAX, IND, IND0, NPFILL, LM
    REAL(q) :: SUMR, SUMI
    REAL(qn),PARAMETER :: ONE=1,ZERO=0
    INTEGER, PARAMETER :: NLM=__align_DIM(101)

!$  INTEGER i,ndim
!$  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM

#ifdef avoidalloc
    REAL(qn) :: WORK(ndata*NONLR_S%IRMAX*NSIM),TMP(NLM, ndata*NSIM)
    GDEF     :: CPROJ(WDES1%NPRO_TOT,NSIM)
#else
    REAL(qn),ALLOCATABLE:: WORK(:),TMP(:,:)
    GDEF,ALLOCATABLE    :: CPROJ(:,:)
#endif
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: WORK,TMP
#endif

    PROFILING_START('rpromu_hf')

#ifndef avoidalloc
    ALLOCATE(CPROJ(WDES1%NPRO_TOT,NSIM),WORK(ndata*NONLR_S%IRMAX*NSIM),TMP(NLM, ndata*NSIM))
#endif

    CPROJ=0

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

#ifndef avoid_reduction
!$OMP PARALLEL DO &
!$OMP NUM_THREADS(omp_nthreads_nonlr_rspace) &
!$OMP PRIVATE(i,WORK,TMP,LMBASE,ISPIRAL,ISPINOR,NLIIND,NIS,NT,LMMAXC,NI,INDMAX,IND0,NPFILL,NP,IND,IP,LM,SUMR,SUMI) &
!$OMP SHARED(NONLR_S,WDES1,W1,GCR,NSIM) &
!$OMP REDUCTION(+:CPROJ)
#else
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(NI,NT,LMMAXC,INDMAX,LMBASE,NLIIND,IND0,NPFILL,NP,WORK,IND,IP,TMP,LM,SUMR,SUMI)
#endif
!$  omp: DO i=1,ndim
!=======================================================================
! loop over ions
!=======================================================================
    LMBASE= 0

    ISPIRAL = 1
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
#ifndef _OPENMP
       NLIIND= 0
       NIS=1

       typ: DO NT=1,NONLR_S%NTYP
          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) GOTO 600

          ion: DO NI=NIS,NONLR_S%NITYP(NT)+NIS-1
!=======================================================================
!  extract the relevant points for this ion
!=======================================================================
             INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
             IF (INDMAX == 0) GOTO 100
#else
#ifdef avoid_reduction
!$OMP DO
#endif
       ion: DO NI=1,NONLR_S%NIONS
          NT=NONLR_S%ITYP(NI)

          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) CYCLE ion

          INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
          IF (INDMAX==0) CYCLE ion

          LMBASE=NONLR_S%LMBASE(NI)+ISPINOR*NONLR_S%LMBASE(NONLR_S%NIONS+1)
          NLIIND=NONLR_S%NLIBASE(NI __NOACC_omp_arg(i))
#endif
!            PROFILING_START('rpromu_hf_gather')

             IND0=0
             NPFILL=0

             DO NP=1,NSIM

                   IF (ASSOCIATED(NONLR_S%CRREXP)) THEN
                      CALL CRREXP_MUL_GWAVE( INDMAX, NONLR_S%CRREXP(1,NI,ISPIRAL __NOACC_omp_arg(i)),NONLR_S%NLI(1,NI __NOACC_omp_arg(i)), &
                           GCR(ISPINOR*WDES1%GRID%MPLWV+1,NP),WORK(IND0+1),WORK(IND0+(NONLR_S%IRMAX+1)))
                   ELSE
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP)
                      DO IND=1,INDMAX
                         IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                         WORK(IND+IND0) = REAL( GCR(IP, NP) ,KIND=qn)
                      ENDDO
                   ENDIF

                   IND0=IND0+ndata * NONLR_S%IRMAX
                   NPFILL=NPFILL+1

             ENDDO

!            PROFILING_STOP('rpromu_hf_gather')

!=======================================================================
! loop over composite indexes L,M
!=======================================================================
!            PROFILING_START('rpromu_hf_gemm')
#ifdef RPROMU_DGEMV
             DO NP=0,ndata*NPFILL-1
                CALL DGEMVn( 'T' , INDMAX, LMMAXC, ONE , NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), &
                     INDMAX, WORK(1+NONLR_S%IRMAX*NP) , 1 , ZERO ,  TMP(1,NP+1), 1)
             ENDDO
#else
             CALL DGEMMn( 'T', 'N' , LMMAXC,  ndata*NPFILL, INDMAX, ONE, &
                  NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), INDMAX, WORK(1), NONLR_S%IRMAX, &
                  ZERO,  TMP(1,1), NLM )
#endif
!            PROFILING_STOP('rpromu_hf_gemm')
!            PROFILING_START('rpromu_hf_scatter')

             IND0=0
             DO NP=1,NSIM
                   l_loop: DO LM=1,LMMAXC
#ifdef gammareal
#ifndef _OPENMP
                      CPROJ(LM+LMBASE,NP)=TMP(LM,1+IND0)*WDES1%RINPL
#else
                      CPROJ(LM+LMBASE,NP)=CPROJ(LM+LMBASE,NP)+TMP(LM,1+IND0)*WDES1%RINPL
#endif
#else
                      SUMR=TMP(LM,1+IND0)
                      SUMI=TMP(LM,2+IND0)
#ifndef _OPENMP
                      CPROJ(LM+LMBASE,NP)=GREAL(CMPLX( SUMR , SUMI ,KIND=q) *WDES1%RINPL)
#else
                      CPROJ(LM+LMBASE,NP)=CPROJ(LM+LMBASE,NP)+GREAL(CMPLX( SUMR , SUMI ,KIND=q) *WDES1%RINPL)
#endif
#endif
                   ENDDO l_loop
                   IND0=IND0+ndata
             ENDDO

!            PROFILING_STOP('rpromu_hf_scatter')
#ifndef _OPENMP
100          LMBASE= LMMAXC+LMBASE
             NLIIND= LMMAXC*INDMAX+NLIIND
          ENDDO ion

600       NIS = NIS+NONLR_S%NITYP(NT)
       ENDDO typ
#else
       ENDDO ion
#ifdef avoid_reduction
!$OMP END DO
#endif
#endif
       IF (NONLR_S%LSPIRAL) ISPIRAL=2
    ENDDO spinor
!$  ENDDO omp
#ifdef avoid_reduction
!$OMP END PARALLEL
#else
!$OMP END PARALLEL DO
#endif

!   PROFILING_START('rpromu_hf_dis')

    ! distribute the projected wavefunctions to nodes
    DO NP=1,NSIM
       CALL DIS_PROJ(WDES1,CPROJ(1,NP),W1(NP)%CPROJ(1))
    ENDDO

!   PROFILING_STOP('rpromu_hf_dis')

#ifndef avoidalloc
    DEALLOCATE(CPROJ,WORK,TMP)
#endif

    PROFILING_STOP('rpromu_hf')

    RETURN
  END SUBROUTINE RPROMU_HF
#else
  SUBROUTINE RPROMU_HF(NONLR_S, WDES1, W1, NSIM, GCR, LD)
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT NONE

    INTEGER NSIM
    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)    WDES1
    TYPE (wavefun1)    W1(NSIM)
    INTEGER :: LD  ! leading dimension of GCR
    GDEF :: GCR(LD, NSIM)
  ! local
    INTEGER :: NP, IP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NT, LMMAXC, &
         NI, INDMAX, IND

!$  INTEGER i,ndim
!$  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM

#ifdef avoidalloc
    REAL(qn) :: WORK(NONLR_S%IRMAX,2)
    GDEF     :: CPROJ(WDES1%NPRO_TOT,NSIM)
#else
    REAL(qn),ALLOCATABLE:: WORK(:,:)
    GDEF,ALLOCATABLE    :: CPROJ(:,:)
#endif
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: WORK
#endif

    PROFILING_START('rpromu_hf')

#ifndef avoidalloc
    ALLOCATE(CPROJ(WDES1%NPRO_TOT,NSIM),WORK(NONLR_S%IRMAX,2))
#endif

    CPROJ=0

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

#ifndef avoid_reduction
!$OMP PARALLEL DO &
!$OMP NUM_THREADS(omp_nthreads_nonlr_rspace) &
!$OMP PRIVATE(i,WORK,LMBASE,ISPIRAL,ISPINOR,NLIIND,NT,LMMAXC,NI,INDMAX,NP,IND,IP) &
!$OMP SHARED(NONLR_S,WDES1,W1,GCR,NSIM) &
!$OMP REDUCTION(+:CPROJ)
#else
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(NI,NT,LMMAXC,INDMAX,LMBASE,NLIIND,NP,WORK,IND,IP)
#endif
!$  omp: DO i=1,ndim
!=======================================================================
! loop over ions
!=======================================================================
    LMBASE= 0

    ISPIRAL = 1
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
#ifdef avoid_reduction
!$OMP DO
#endif
       ion: DO NI=1,NONLR_S%NIONS
          NT=NONLR_S%ITYP(NI)

          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) CYCLE ion

             INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
             IF (INDMAX==0) CYCLE ion

             LMBASE=NONLR_S%LMBASE(NI)+ISPINOR*NONLR_S%LMBASE(NONLR_S%NIONS+1)
             NLIIND=NONLR_S%NLIBASE(NI __NOACC_omp_arg(i))

             DO NP=1,NSIM
!               PROFILING_START('rpromu_hf_gather')

                IF (ASSOCIATED(NONLR_S%CRREXP)) THEN
                   CALL CRREXP_MUL_GWAVE( INDMAX, NONLR_S%CRREXP(1,NI,ISPIRAL __NOACC_omp_arg(i)),NONLR_S%NLI(1,NI __NOACC_omp_arg(i)), &
                        GCR(ISPINOR*WDES1%GRID%MPLWV+1,NP),WORK(1,1),WORK(1,2))
                ELSE
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP)
                   DO IND=1,INDMAX
                      IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                      WORK(IND,1) = REAL( GCR(IP, NP) ,KIND=qn)
                   ENDDO
                ENDIF

!               PROFILING_STOP('rpromu_hf_gather')
!               PROFILING_START('rpromu_hf_gemm')

                CALL GEMMV_AND_SCATTER(INDMAX,LMMAXC,NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), &
                                       WORK(1,1),WORK(1,2),CPROJ(LMBASE+1,NP),WDES1%RINPL)

!               PROFILING_STOP('rpromu_hf_gemm')
             ENDDO
       ENDDO ion
#ifdef avoid_reduction
!$OMP END DO
#endif
       IF (NONLR_S%LSPIRAL) ISPIRAL=2
    ENDDO spinor
!$  ENDDO omp
#ifdef avoid_reduction
!$OMP END PARALLEL
#else
!$OMP END PARALLEL DO
#endif

!   PROFILING_START('rpromu_hf_dis')

    ! distribute the projected wavefunctions to nodes
    DO NP=1,NSIM
       CALL DIS_PROJ(WDES1,CPROJ(1,NP),W1(NP)%CPROJ(1))
    ENDDO

!   PROFILING_STOP('rpromu_hf_dis')

#ifndef avoidalloc
    DEALLOCATE(CPROJ,WORK)
#endif

    PROFILING_STOP('rpromu_hf')

    RETURN
  END SUBROUTINE RPROMU_HF
#endif

#ifdef _OPENACC
  SUBROUTINE RPROMU_HF_ACC(NONLR_S, WDES1, CDLM, LDC, NSIM, GCR, LD)
    USE openmp, ONLY : omp_nthreads_acc
    IMPLICIT NONE

    TYPE (nonlr_struct) :: NONLR_S
    TYPE (wavedes1)     :: WDES1
    GDEF                :: CDLM(LDC,NSIM)
    GDEF                :: GCR(LD,NSIM)
    INTEGER             :: LDC             ! leading dimension of CDLM
    INTEGER             :: LD              ! leading dimension of GCR
    INTEGER             :: NSIM

  ! local
    INTEGER :: NP, IP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NT, LMMAXC, NI, INDMAX, IND, LM
    GDEF :: CTMP

    INTEGER, PARAMETER :: BLOCKSIZE=512
    INTEGER :: IBLOCK

!$  INTEGER i,ndim

#ifdef avoidalloc
    GDEF :: GWORK(BLOCKSIZE),CPROJ(WDES1%NPRO_TOT,NSIM)
#else
    GDEF, ALLOCATABLE :: GWORK(:),CPROJ(:,:)
#endif

    PROFILING_START('rpromu_hf_acc')

#ifndef avoidalloc
    ALLOCATE(GWORK(BLOCKSIZE),CPROJ(WDES1%NPRO_TOT,NSIM))
#endif

!$ACC ENTER DATA CREATE(CPROJ) ASYNC(ACC_ASYNC_Q)
!$ACC KERNELS PRESENT(CPROJ) ASYNC(ACC_ASYNC_Q)
    CPROJ=0
!$ACC END KERNELS

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

!$  omp: DO i=1,ndim
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1

       ISPIRAL=1; IF (NONLR_S%LSPIRAL) ISPIRAL=ISPINOR+1

!=======================================================================
! loop over ions
!=======================================================================
!$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(CPROJ,WDES1,WDES1%GRID,NONLR_S,GCR) &
!$ACC PRIVATE(NT,LMMAXC,INDMAX,LMBASE,NLIIND,IBLOCK) ASYNC(ACC_ASYNC_Q)
       ion: DO NI=1,NONLR_S%NIONS
          DO NP=1,NSIM

             NT=NONLR_S%ITYP(NI)

             LMMAXC=NONLR_S%LMMAX(NT)
             IF (LMMAXC==0) CYCLE ! ion

             INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
             IF (INDMAX==0) CYCLE ! ion

             LMBASE=NONLR_S%LMBASE(NI)+ISPINOR*NONLR_S%LMBASE(NONLR_S%NIONS+1)
             NLIIND=NONLR_S%NLIBASE(NI __NOACC_omp_arg(i))

!$ACC LOOP SEQ PRIVATE(GWORK)
             block: DO IBLOCK=1,INDMAX,BLOCKSIZE
!$ACC CACHE(GWORK(1:BLOCKSIZE))
!$ACC LOOP WORKER VECTOR PRIVATE(IP)
                DO IND=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,INDMAX)
                   IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
#ifdef gammareal
                   GWORK(IND-IBLOCK+1)=GCR(IP,NP)
#else
                   GWORK(IND-IBLOCK+1)=GCR(IP,NP)*NONLR_S%CRREXP(IND,NI,ISPIRAL __NOACC_omp_arg(i))
#endif
                ENDDO

!$ACC LOOP WORKER PRIVATE(CTMP)
                DO LM=1,LMMAXC
                   CTMP=0
!$ACC LOOP VECTOR REDUCTION(+:CTMP)
                   DO IND=IBLOCK,MIN(IBLOCK+BLOCKSIZE-1,INDMAX)
                      CTMP=CTMP+GWORK(IND-IBLOCK+1)*NONLR_S%RPROJ(IND+(LM-1)*INDMAX + NLIIND __NOACC_omp_arg(i))
                   ENDDO
                   CPROJ(LM+LMBASE,NP)=CPROJ(LM+LMBASE,NP)+CTMP*WDES1%RINPL
                ENDDO

             ENDDO block

          ENDDO
       ENDDO ion

    ENDDO spinor
!$  ENDDO omp

#ifdef MPI
    IF (WDES1%COMM_INB%NCPU>1) THEN
    ! distribute the projected wavefunctions to nodes
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
       DO NP=1,NSIM
          CALL DIS_PROJ(WDES1,CPROJ(1,NP),CDLM(1,NP))
       ENDDO
!$OMP END PARALLEL DO
    ELSE
#endif
    ! in this case a simple copy suffices
!$ACC KERNELS PRESENT(CDLM,WDES1,CPROJ) ASYNC(ACC_ASYNC_Q)
       CDLM(1:WDES1%NPRO,1:NSIM)=CPROJ(1:WDES1%NPRO,1:NSIM)
!$ACC END KERNELS
#ifdef MPI
    ENDIF
#endif

!$ACC EXIT DATA DELETE(CPROJ) ASYNC(ACC_ASYNC_Q)
#ifndef avoidalloc
    DEALLOCATE(CPROJ,GWORK)
#endif

    PROFILING_STOP('rpromu_hf_acc')

    RETURN

  END SUBROUTINE RPROMU_HF_ACC
#endif

!****************** SUBROUTINE RPRO **********************************
!
!  this subroutine  calculates the projection of all bands onto the
!  real space projection operators doing a set of
!  bands at the same time
!  
!*********************************************************************

  SUBROUTINE RPRO(NONLR_S,WDES,W,GRID,NK)
    IMPLICIT NONE
    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes)     WDES
    TYPE (wavedes1)    WDES1
    TYPE (wavespin)    W
    TYPE (grid_3d)     GRID
    INTEGER NK

    CALL RPRO_ISP(NONLR_S,WDES,W,GRID,0,NK)

  END SUBROUTINE RPRO


  SUBROUTINE RPRO_ISP(NONLR_S,WDES,W,GRID,ISP_SWITCH,NK)
    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W
    TYPE (grid_3d)     GRID
    INTEGER ISP_SWITCH       ! 0 all spin components, 1 or 2 only selected
    INTEGER NK
  ! local variables
    TYPE (wavedes1)  WDES1
    TYPE (wavefun1), TARGET :: W1(WDES%NSIM)
    LOGICAL :: LDO(WDES%NSIM)
    INTEGER :: NSIM, NT, N, NPL, NGVECTOR, ISP_START, ISP_END, ISP, NUP, NN, NNP, ISPINOR

    PROFILING_START('rpro_isp')

    LDO=.TRUE.
    NSIM = WDES%NSIM
    DO NT=1,NONLR_S%NTYP
       IF (NONLR_S%LMMAX(NT)/=0) GOTO 300
    ENDDO
    ! shortcut for NC potentials
    PROFILING_STOP('rpro_isp')
    RETURN

300 CONTINUE
    ! setup descriptor
!$ACC ENTER DATA CREATE(WDES1) __IF_ASYNC__
    CALL SETWDES(WDES,WDES1,NK)
    CALL SETWGRID_OLD(WDES1,GRID)

!$ACC ENTER DATA CREATE(W1(:)) __IF_ASYNC__
    DO N=1,NSIM
       CALL NEWWAV_R(W1(N),WDES1)
    ENDDO

    NPL=WDES%NPLWKP(NK)
    NGVECTOR=WDES%NGVECTOR(NK)

    IF (ISP_SWITCH==1 .OR. ISP_SWITCH==2) THEN
       ISP_START=ISP_SWITCH
       ISP_END  =ISP_SWITCH
    ELSE
       ISP_START=1
       ISP_END  =WDES%ISPIN
    ENDIF

    DO ISP=ISP_START,ISP_END
       DO N=1,WDES%NBANDS,NSIM
          NUP=MIN(N+NSIM-1,WDES%NBANDS)
          DO NN=N,NUP
             NNP=NN-N+1
             CALL SETWAV(W,W1(NNP),WDES1,NN,ISP)
             DO ISPINOR=0,WDES%NRSPINORS-1
                CALL FFTWAV(NGVECTOR,WDES%NINDPW(1,NK),W1(NNP)%CR(1+ISPINOR*WDES1%GRID%MPLWV),W1(NNP)%CW(1+ISPINOR*NGVECTOR),GRID)
             ENDDO
          ENDDO
#ifdef _OPENACC
          IF (ACC_EXEC_ON) THEN
             CALL RPROMU_ACC(NONLR_S, WDES1, W1, NUP-N+1, LDO)
             CYCLE
          ENDIF
#endif
          IF (NSIM/=1) THEN
             CALL RPROMU(NONLR_S,WDES1,W1,NUP-N+1,LDO)
          ELSE
             CALL RPRO1(NONLR_S,WDES1,W1(1))
          ENDIF
       ENDDO
    ENDDO

    DO N=1,NSIM
       CALL DELWAV_R(W1(N))
    ENDDO
!$ACC EXIT DATA DELETE(W1(:)) __IF_ASYNC__
#ifdef _OPENACC
    CALL ACC_DELETE_TYPED_VAR(WDES1)
#endif
    PROFILING_STOP('rpro_isp')

    RETURN
  END SUBROUTINE RPRO_ISP


!****************** SUBROUTINE RACCT *********************************
!
!  this subroutine  calculates the non local part of the gradient for
!  all bands.
!  it is only for performance testing
!
!*********************************************************************

  SUBROUTINE RACCT(NONLR_S,WDES,W,GRID,CDIJ,CQIJ,ISP,LMDIM, NK)
    USE mpimy
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes)     WDES
    TYPE (wavedes1)    WDES1
    TYPE (wavespin)    W
    TYPE (grid_3d)     GRID
    OVERLAP  CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ),CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    ! local variables
    TYPE (wavefun1)  W1(WDES%NSIM)
    LOGICAL ::       LDO(WDES%NSIM)
    COMPLEX(q),ALLOCATABLE :: CWORK(:,:),CWORK2(:)
    REAL(q) :: EVALUE(WDES%NSIM)

    LDO=.TRUE.
    NSIM = WDES%NSIM
    DO NT=1,NONLR_S%NTYP
       IF (NONLR_S%LMMAX(NT)/=0) GOTO 300
    ENDDO
    RETURN

300 CONTINUE
    ALLOCATE(CWORK(GRID%MPLWV,NSIM),CWORK2(GRID%MPLWV))


    ! setup descriptor
    CALL SETWDES(WDES,WDES1,NK); CALL SETWGRID_OLD(WDES1,GRID)

    NPL=WDES%NPLWKP(NK)
    NGVECTOR=WDES%NGVECTOR(NK)

    DO ISP=1,WDES%ISPIN
       DO N=1,WDES%NBANDS,NSIM
          NUP=MIN(N+NSIM-1,WDES%NBANDS)
          CWORK=0
          DO NN=N,NUP
             NNP=NN-N+1
             CALL SETWAV(W,W1(NNP),WDES1,NN,ISP)
             EVALUE(NNP)=W%CELEN(N,1,ISP)
          ENDDO
          CALL RACCMU(NONLR_S,WDES1,W1, LMDIM,CDIJ(1,1,1,ISP),CQIJ(1,1,1,ISP),EVALUE,CWORK(1,1), &
               WDES1%GRID%MPLWV*WDES1%NRSPINORS, NSIM, LDO)
          DO NN=N,NUP
             NNP=NN-N+1
             DO  ISPINOR=0,WDES%NRSPINORS-1
                CALL FFTEXT(NGVECTOR,WDES%NINDPW(1,NK),CWORK(1+ISPINOR*GRID%MPLWV,NNP),CWORK2(1+ISPINOR*NGVECTOR),GRID,.FALSE.)
             ENDDO
          ENDDO


       ENDDO
    ENDDO
    DEALLOCATE(CWORK,CWORK2)

    RETURN
  END SUBROUTINE RACCT


!****************** SUBROUTINE RLACC *********************************
!
!  subroutine for calculating the non local contribution of
!  the Hamiltonian, using real space projection scheme
!  the result of the wavefunction projected on the projection operatores
!  must be given in CPROJ
!  the result is added to  CRACC
!                -----
!*********************************************************************

  SUBROUTINE RACC(NONLR_S, W1, CDIJ, CQIJ, ISP, EVALUE,  CRACC)
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (nonlr_struct) NONLR_S
    TYPE (wavefun1)    W1
    COMPLEX(q)  CRACC(:)
    OVERLAP  CDIJ(:,:,:,:),CQIJ(:,:,:,:)
    INTEGER  ISP
  ! work arrays
    GDEF :: CRESUL(W1%WDES1%NPROD)

    CALL OVERL1(W1%WDES1, SIZE(CDIJ,1), CDIJ(1,1,1,ISP), CQIJ(1,1,1,ISP), EVALUE, W1%CPROJ(1),CRESUL(1))
    CALL RACC0(NONLR_S, W1%WDES1, CRESUL(1), CRACC(1))

    RETURN
  END SUBROUTINE RACC


!****************** SUBROUTINE RACCMU *******************************
!
!  subroutine for calculating the non local contribution of
!  the Hamiltonian, using real space projection scheme
!  for a set of bands simultaneously
!  the result is added to  CRACC
!                -----
!*********************************************************************

!
! scheduled for removal
!
  SUBROUTINE RACCMU(NONLR_S,WDES1,W1, &
       &     LMDIM,CDIJ,CQIJ,EVALUE, CRACC,LD, NSIM, LDO)
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)     WDES1
    TYPE (wavedes)      WDES
    TYPE (wavefun1)     W1(NSIM)

    COMPLEX(q) CRACC(LD, NSIM)
    OVERLAP    CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
         CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    REAL(q)    EVALUE(NSIM)
    LOGICAL    LDO(NSIM)
    ! work arrays
    GDEF :: CRESUL(WDES1%NPROD,NSIM)

    PROFILING_START('raccmu')

DOACC CALL ACC_SYNC_ASYNC_Q
!$ACC ENTER DATA CREATE(CRESUL) IF(ACC_EXEC_ON)
    DO NP=1,NSIM
       IF (LDO(NP)) THEN

          SET_ACC_ASYNC_Q(NP)

          CALL OVERL1(WDES1, LMDIM,CDIJ,CQIJ, EVALUE(NP), W1(NP)%CPROJ(1),CRESUL(1,NP))
       ENDIF
    ENDDO
#ifdef _OPENACC
    IF (ACC_EXEC_ON) THEN
       CALL RACC0MU_ACC(NONLR_S,WDES1,CRESUL(1,1),CRACC(1,1),LD, NSIM,LDO)
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
       ACC_ASYNC_Q=ACC_ASYNC_SYNC
!$ACC EXIT DATA DELETE(CRESUL) WAIT(ACC_ASYNC_Q)
       PROFILING_STOP('raccmu')
       RETURN
    ENDIF
#endif
    IF (NSIM/=1) THEN
       CALL RACC0MU(NONLR_S,WDES1,CRESUL(1,1),CRACC(1,1),LD, NSIM,LDO)
    ELSE
       CALL RACC0(NONLR_S,WDES1,CRESUL(1,1),CRACC(1,1))
    ENDIF

    PROFILING_STOP('raccmu')

    RETURN
  END SUBROUTINE RACCMU


  SUBROUTINE RACCMU_(NONLR_S, WDES1, W1, CDIJ, CQIJ, ISP, EVALUE, CRACC)
    USE openmp, ONLY : omp_nthreads_acc
    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)     WDES1
    TYPE (wavefun1)     W1(:)
    OVERLAP CDIJ(:,:,:,:), CQIJ(:,:,:,:)
    INTEGER :: ISP
    REAL(q)    EVALUE(:)
    COMPLEX(q) CRACC(:,:)
    ! work arrays
    GDEF :: CRESUL(WDES1%NPROD,SIZE(W1))
    INTEGER :: NP

    PROFILING_START('raccmu_')

#ifdef _OPENACC
    IF (ACC_EXEC_ON) THEN
       CALL ACC_SYNC_ASYNC_Q

!$ACC ENTER DATA CREATE(CRESUL) COPYIN(EVALUE) IF(ACC_EXEC_ON)
       CALL OVERL_MU(WDES1,SIZE(CDIJ,1),CDIJ(1,1,1,ISP),CQIJ(1,1,1,ISP),EVALUE(1),W1(1),CRESUL(1,1),SIZE(W1))
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)

!!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc) SCHEDULE(STATIC)
!!    DO NP=1,SIZE(W1)
!!       IF (W1(NP)%LDO) THEN
!!
!!          SET_ACC_ASYNC_Q(NP)
!!
!!          CALL OVERL1(WDES1,SIZE(CDIJ,1),CDIJ(1,1,1,ISP),CQIJ(1,1,1,ISP),EVALUE(NP),W1(NP)%CPROJ(1),CRESUL(1,NP))
!!       ENDIF
!!    ENDDO
!!!$OMP END PARALLEL DO

       CALL RACC0MU_ACC(NONLR_S,WDES1,CRESUL(1,1),CRACC(1,1),SIZE(CRACC,1),SIZE(W1),W1%LDO)

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
       ACC_ASYNC_Q=ACC_ASYNC_SYNC
!$ACC EXIT DATA DELETE(CRESUL,EVALUE)
       PROFILING_STOP('raccmu_')
       RETURN
    ENDIF
#endif

    DO NP=1,SIZE(W1)
       IF (W1(NP)%LDO) THEN
          CALL OVERL1(WDES1,SIZE(CDIJ,1),CDIJ(1,1,1,ISP),CQIJ(1,1,1,ISP),EVALUE(NP),W1(NP)%CPROJ(1),CRESUL(1,NP))
       ENDIF
    ENDDO

    IF (SIZE(W1)/=1) THEN
       CALL RACC0MU(NONLR_S,WDES1,CRESUL(1,1),CRACC(1,1),SIZE(CRACC,1),SIZE(W1),W1%LDO)
    ELSE
       CALL RACC0(NONLR_S,WDES1,CRESUL(1,1),CRACC(1,1))
    ENDIF

    PROFILING_STOP('raccmu_')

    RETURN
  END SUBROUTINE RACCMU_


!
! same as before but for complex EVALUE
! scheduled for removal
!
  SUBROUTINE RACCMU_C(NONLR_S,WDES1,W1, &
       &     LMDIM,CDIJ,CQIJ,CEVALUE, CRACC,LD, NSIM, LDO)
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)     WDES1
    TYPE (wavedes)      WDES
    TYPE (wavefun1)     W1(NSIM)

    COMPLEX(q) CRACC(LD, NSIM)
    OVERLAP    CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
         CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    COMPLEX(q) CEVALUE(NSIM)
    LOGICAL    LDO(NSIM)
    ! work arrays
    GDEF :: CRESUL(WDES1%NPROD,NSIM)

    PROFILING_START('raccmu_c')

!$ACC ENTER DATA CREATE(CRESUL) IF(ACC_EXEC_ON)
    DO NP=1,NSIM
       IF (LDO(NP)) THEN

          SET_ACC_ASYNC_Q(NP)

          CALL OVERL1_C(WDES1, LMDIM,CDIJ,CQIJ, CEVALUE(NP), W1(NP)%CPROJ(1),CRESUL(1,NP))
       ENDIF
    ENDDO
#ifdef _OPENACC
    IF (ACC_EXEC_ON) THEN
       CALL RACC0MU_ACC(NONLR_S,WDES1,CRESUL(1,1),CRACC(1,1),LD, NSIM,LDO)
!$ACC EXIT DATA DELETE(CRESUL) __IF_ASYNC__
       PROFILING_STOP('raccmu_c')
       RETURN
    ENDIF
#endif
    IF (NSIM/=1) THEN
       CALL RACC0MU(NONLR_S,WDES1,CRESUL(1,1),CRACC(1,1),LD, NSIM,LDO)
    ELSE
       CALL RACC0(NONLR_S,WDES1,CRESUL(1,1),CRACC(1,1))
    ENDIF

    PROFILING_STOP('raccmu_c')

    RETURN
  END SUBROUTINE RACCMU_C


  SUBROUTINE RACCMU_C_(NONLR_S, WDES1, W1, &
       &     CDIJ, CQIJ, ISP, CEVALUE, CRACC)
    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)     WDES1
    TYPE (wavefun1)     W1(:)
    OVERLAP CDIJ(:,:,:,:), CQIJ(:,:,:,:)
    INTEGER :: ISP
    COMPLEX(q) CEVALUE(:)
    COMPLEX(q) CRACC(:,:)
    ! work arrays
    GDEF :: CRESUL(WDES1%NPROD,SIZE(W1))
    INTEGER :: NP

DOACC CALL ACC_SYNC_ASYNC_Q
!$ACC ENTER DATA CREATE(CRESUL) IF(ACC_EXEC_ON)
    DO NP=1,SIZE(W1)
       IF (W1(NP)%LDO) THEN

          SET_ACC_ASYNC_Q(NP)

          CALL OVERL1_C(WDES1, SIZE(CDIJ,1),CDIJ(1,1,1,ISP),CQIJ(1,1,1,ISP), CEVALUE(NP), W1(NP)%CPROJ(1),CRESUL(1,NP))
       ENDIF
    ENDDO
#ifdef _OPENACC
    IF (ACC_EXEC_ON) THEN
       CALL RACC0MU_ACC(NONLR_S,WDES1,CRESUL(1,1),CRACC(1,1),SIZE(CRACC,1), SIZE(W1),W1%LDO)
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
       ACC_ASYNC_Q=ACC_ASYNC_SYNC
!$ACC EXIT DATA DELETE(CRESUL) WAIT(ACC_ASYNC_Q)
       RETURN
    ENDIF
#endif
    IF (SIZE(W1)/=1) THEN
       CALL RACC0MU(NONLR_S,WDES1,CRESUL(1,1),CRACC(1,1),SIZE(CRACC,1), SIZE(W1),W1%LDO)
    ELSE
       CALL RACC0(NONLR_S,WDES1,CRESUL(1,1),CRACC(1,1))
    ENDIF

    RETURN
  END SUBROUTINE RACCMU_C_

! here CDIJ and QCIJ are always complex
  SUBROUTINE RACCMU_CCDIJ(NONLR_S,WDES1,W1, &
       &     LMDIM,CDIJ,CQIJ,EVALUE, CRACC,LD, NSIM, LDO)
    IMPLICIT COMPLEX(q) (C)
    IMPLICIT REAL(q) (A-B,D-H,O-Z)

    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)    WDES1
    TYPE (wavedes)    WDES
    TYPE (wavefun1)    W1(NSIM)

    COMPLEX(q) CRACC(LD, NSIM)
    COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
         CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
    REAL(q)    EVALUE(NSIM)
    LOGICAL    LDO(NSIM)
    ! work arrays
    GDEF :: CRESUL(WDES1%NPROD,NSIM)

    DO NP=1,NSIM
       IF (LDO(NP)) THEN
          CALL OVERL1_CCDIJ(WDES1, LMDIM,CDIJ,CQIJ, EVALUE(NP), W1(NP)%CPROJ(1),CRESUL(1,NP))
       ENDIF
    ENDDO
    IF (NSIM/=1) THEN
       CALL RACC0MU(NONLR_S,WDES1,CRESUL(1,1),CRACC(1,1),LD, NSIM,LDO)
    ELSE
       CALL RACC0(NONLR_S,WDES1,CRESUL(1,1),CRACC(1,1))
    ENDIF

    RETURN
  END SUBROUTINE RACCMU_CCDIJ


!****************** SUBROUTINE RNLPR ***********************************
! subroutine for calculating the non-local energy per ion
!
! E(ION,k) = SUM(BAND,L,M) Z(ION,BAND,L,M,k) CONJG( Z(ION,BAND,L,M,k))
!
!> @details @ref openmp :
!> under OpenMP the nested loop over \"types\" + \"ions-of-typ\"
!> is replaced by a single loop over \"all ions\", and is
!> distributed over the available threads (OMP PARALLEL DO).
!
!***********************************************************************

  SUBROUTINE RNLPR(GRID,NONLR_S, P, LATT_FIN1, LATT_FIN2, LATT_CUR, W, WDES, &
       &    LMDIM, CDIJ, CQIJ, ENL)
    USE pseudo
    USE mpimy
    USE lattice
    USE constant
    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    TYPE (potcar)      P(NONLR_S%NTYP)
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W,WTMP
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR,LATT_FIN1,LATT_FIN2
    INTEGER LMDIM
    OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    OVERLAP CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    REAL(q) ENL(NONLR_S%NIONS)
  ! local
    INTEGER NK, ISP, N, ISPINOR, ISPINOR_, LBASE, LBASE_, NT, LMMAXC, &
         NI, L, LP
    REAL(q) EVALUE, WEIGHT, RTMP
    GDEF,ALLOCATABLE,TARGET :: CPROW(:,:,:,:)
    REAL(q) DISPL(3,NONLR_S%NIONS)

    ALLOCATE(CPROW(WDES%NPROD,WDES%NBANDS,WDES%NKPTS,WDES%ISPIN))
!$ACC ENTER DATA CREATE(CPROW) __IF_ASYNC__
!=======================================================================
!  calculate the projection operator
!=======================================================================
    NK=1
    DISPL=0
    CALL RSPHER_ALL(GRID,NONLR_S,LATT_FIN1,LATT_FIN2,LATT_CUR,DISPL,DISPL, 1)

!$ACC ENTER DATA CREATE(ENL) __IF_ASYNC__
!$ACC KERNELS PRESENT(ENL) __IF_ASYNC__
    ENL=0
!$ACC END KERNELS

!$ACC ENTER DATA CREATE(WTMP) __IF_ASYNC__
    WTMP=W
    WTMP%CPROJ => CPROW  ! relink the CPROJ array to temporary workspace

    kpoint: DO NK=1,WDES%NKPTS
#ifdef MPI
       IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
       CALL PHASER(GRID,LATT_CUR,NONLR_S, NK,WDES)
       CALL RPRO(NONLR_S,WDES,WTMP,GRID,NK)

!$ACC PARALLEL LOOP COLLAPSE(force:5) GANG PRESENT(WDES,CDIJ,W,W%CPROJ,CPROW,CQIJ,ENL,W%CELEN,W%FERWE) &
!$ACC PRIVATE(RTMP,EVALUE,WEIGHT,LBASE,LBASE_,NT,LMMAXC) __IF_ASYNC__
       spin: DO ISP=1,WDES%ISPIN
!=======================================================================
!  sum up to give the non-local energy per ion
!=======================================================================

          band: DO N=1,WDES%NBANDS
             EVALUE=W%CELEN(N,NK,ISP)
             WEIGHT=WDES%WTKPT(NK)*W%FERWE(N,NK,ISP)*WDES%RSPIN

             spinor: DO ISPINOR=0,WDES%NRSPINORS-1
                DO ISPINOR_=0,WDES%NRSPINORS-1
NOACC !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(NI,NT,LMMAXC,LBASE,LBASE_,L,LP,RTMP) &
NOACC !$OMP SHARED(WDES,ENL,WEIGHT,W,CPROW,N,NK,ISP,ISPINOR,ISPINOR_,CDIJ,CQIJ,EVALUE)
                   DO NI=1,WDES%NIONS
                      NT=WDES%ITYP(NI)
                      LMMAXC=WDES%LMMAX(NT)
                      IF (LMMAXC==0) CYCLE
                      LBASE =WDES%LMBASE(NI)+ISPINOR *WDES%NPRO/2
                      LBASE_=WDES%LMBASE(NI)+ISPINOR_*WDES%NPRO/2

                      RTMP = 0.0_q
!$ACC LOOP COLLAPSE(2) REDUCTION(+:RTMP)
                      DO L=1 ,LMMAXC
                         DO LP=1,LMMAXC
                            RTMP=RTMP+WEIGHT*W%CPROJ(LBASE_+LP,N,NK,ISP)*GCONJG(CPROW(LBASE+L,N,NK,ISP))* &
                                 (CDIJ(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)-EVALUE*CQIJ(LP,L,NI,ISP+ISPINOR_+2*ISPINOR))
                         ENDDO
                      ENDDO
!$ACC ATOMIC UPDATE
                      ENL(NI)=ENL(NI)+RTMP
                   ENDDO
NOACC !$OMP END PARALLEL DO
                ENDDO
             ENDDO spinor

          ENDDO band
       ENDDO spin

    ENDDO kpoint
!$ACC EXIT DATA COPYOUT(ENL) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)
!$ACC EXIT DATA DELETE(WTMP) __IF_ASYNC__

    CALLMPI( M_sum_d( WDES%COMM_KINTER, ENL(1),NONLR_S%NIONS) )
!$ACC EXIT DATA DELETE(CPROW) __IF_ASYNC__
    DEALLOCATE(CPROW)
    RETURN
  END SUBROUTINE RNLPR


!****************** SUBROUTINE STRNLR **********************************
!
!  subroutine for calculating the non-local contributions to stress,
!  use central differences
!  all components to the stress tensor are calculated
!  except if ISIF = 1
!
!  uncomment CTEST-lines if you want to test finit-differences
!
!***********************************************************************

  SUBROUTINE STRNLR(GRID,NONLR_S,P,LATT_CUR,W, &
      &    CDIJ,CQIJ, ISIF,FNLSIF)
    USE pseudo
    USE mpimy
    USE lattice
    USE constant
    IMPLICIT NONE

    TYPE (grid_3d)     GRID
    TYPE (nonlr_struct) NONLR_S
    TYPE (potcar)      P(NONLR_S%NTYP)
    TYPE (wavespin)    W
    TYPE (latt)        LATT_CUR
    OVERLAP CDIJ(:,:,:,:)
    OVERLAP CQIJ(:,:,:,:)
    INTEGER ISIF                         ! which componets
    REAL(q) FNLSIF(3,3)                ! result stress tensor
  ! local
    TYPE (latt)        LATT_FIN1,LATT_FIN2
    INTEGER :: IDIR, JDIR, I, J, NI
  ! magnitude used to distort lattice
    REAL(q) :: DIS=fd_displacement
    REAL(q) ::  ENL(NONLR_S%NIONS)

    PROFILING_START('strnlr')

    !TEST
    !      DIS=1E-3
    ! 1000 DIS=DIS/2
    !TEST
!=======================================================================
! initialise non-local forces to zero
!=======================================================================
    FNLSIF=0
!=======================================================================
! calculate the contribution to the energy from the nonlocal
! pseudopotential for elongation of each basis-vector
!=======================================================================
    DO IDIR=1,3
       DO JDIR=1,3

          LATT_FIN1=LATT_CUR
          LATT_FIN2=LATT_CUR
          IF (ISIF==1) THEN
             !  only isotrop pressure
             DO I=1,3
                DO J=1,3
                   LATT_FIN1%A(I,J)=LATT_CUR%A(I,J)*(1+DIS/3)
                   LATT_FIN2%A(I,J)=LATT_CUR%A(I,J)*(1-DIS/3)
                ENDDO
             ENDDO
          ELSE
             !  all directions
             DO I=1,3
                LATT_FIN1%A(IDIR,I)=LATT_CUR%A(IDIR,I)+DIS*LATT_CUR%A(JDIR,I)
                LATT_FIN2%A(IDIR,I)=LATT_CUR%A(IDIR,I)-DIS*LATT_CUR%A(JDIR,I)
             ENDDO
          ENDIF
          CALL LATTIC(LATT_FIN1)
          CALL LATTIC(LATT_FIN2)

          CALL RNLPR(GRID,NONLR_S,P,LATT_FIN1,LATT_FIN2,LATT_CUR,W,W%WDES, &
               &    SIZE(CDIJ,1),CDIJ,CQIJ,ENL)

          DO NI=1,NONLR_S%NIONS
             FNLSIF(IDIR,JDIR)=FNLSIF(IDIR,JDIR)+ENL(NI)
          ENDDO
          !
          !  only isotrop pressure terminate loop
          !
          IF (ISIF==1) THEN
             FNLSIF(2,2)= FNLSIF(1,1)
             FNLSIF(3,3)= FNLSIF(1,1)
             GOTO 400 ! terminate (not very clean but who cares)
          ENDIF

       ENDDO
    ENDDO
!=======================================================================
! calculation finished  scale pressure
!=======================================================================
400 CONTINUE

    CALLMPI( M_sum_d(W%WDES%COMM_KIN, FNLSIF, 9))
    FNLSIF=FNLSIF/DIS
    !TEST
    !      WRITE(*,'(E10.3,3E14.7)')DIS,((FNLSIF(I,J),I=1,3),J=1,3)
    !      IF (DIS>1E-10) GOTO 1000
    !TEST

    PROFILING_STOP('strnlr')

    RETURN
  END SUBROUTINE STRNLR


!****************** SUBROUTINE FORNLR **********************************
!
!  subroutine for calculating the non local contribution
!  to the forces acting onto the ions (using simple central finite
!  differences)
!  uncomment CTEST-lines if you want to test finit-differences
!
!> @details @ref openmp :
!> under OpenMP the nested loop over \"types\" + \"ions-of-typ\"
!> is replaced by a single loop over \"all ions\", and is
!> distributed over the available threads (OMP PARALLEL DO).
!
!***********************************************************************

  SUBROUTINE FORNLR(GRID, NONLR_S, P, LATT_CUR, W, &
       &    CDIJ, CQIJ, DISPL0, FORNL)
    USE pseudo
    USE mpimy
    USE lattice
    USE constant

    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    TYPE (potcar)      P(NONLR_S%NTYP)
    TYPE (wavespin), TARGET ::    W
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR

    REAL(q) FORNL(3,NONLR_S%NIONS)
    REAL(q) DISPL0(3,NONLR_S%NIONS)
    OVERLAP CDIJ(:,:,:,:)
    OVERLAP CQIJ(:,:,:,:)
  ! local
    TYPE (wavespin)    WTMP
    TYPE (wavedes), POINTER :: WDES
    REAL(q) ENL(NONLR_S%NIONS), EVALUE, WEIGHT, RTMP
    INTEGER IDIR, NK, ISP, N, ISPINOR, ISPINOR_, LBASE, LBASE_, NT, &
         LMMAXC, NI, L, LP, NIP
    REAL(q) DISPL1(3,NONLR_S%NIONS),DISPL2(3,NONLR_S%NIONS)
    ! allocate required work space
    GDEF,ALLOCATABLE,TARGET :: CPROW(:,:,:,:)
    GDEF,POINTER :: CPROT(:,:,:,:)
  ! magnitude used for finite differences
    REAL(q) :: DIS=fd_displacement

    PROFILING_START('fornlr')

    !TEST
    !      DIS=1E-3
    ! 1000 DIS=DIS/2
    !TEST

    WDES=>W%WDES

    ALLOCATE(CPROW(WDES%NPROD,WDES%NBANDS,WDES%NKPTS,WDES%ISPIN))
!$ACC ENTER DATA CREATE(CPROW) __IF_ASYNC__
!=======================================================================
! initialise non-local forces to zero
!=======================================================================
    FORNL=0
!=======================================================================
! calculate the contribution to the force from the nonlocal
! projection functions for displacement X using central (semianalytical)
! finite differences (about 9 digits precision)
!=======================================================================
    dir: DO IDIR=1,3
!$ACC ENTER DATA CREATE(ENL) __IF_ASYNC__
!$ACC KERNELS PRESENT(ENL) __IF_ASYNC__
       ENL=0
!$ACC END KERNELS

       DISPL1=DISPL0
       DISPL1(IDIR,:)= DISPL0(IDIR,:)-DIS

       DISPL2=DISPL0
       DISPL2(IDIR,:)= DISPL0(IDIR,:)+DIS
       CALL RSPHER_ALL(GRID,NONLR_S,LATT_CUR,LATT_CUR,LATT_CUR, DISPL1, DISPL2, 1)

!$ACC ENTER DATA CREATE(WTMP) __IF_ASYNC__
       WTMP=W
       WTMP%CPROJ => CPROW       ! relink the CPROJ array to temporary workspace

       kpoint: DO NK=1,WDES%NKPTS
#ifdef MPI
          IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
          CALL PHASER(GRID,LATT_CUR,NONLR_S, NK,WDES)
          CALL RPRO(NONLR_S,WDES,WTMP,GRID,NK)

!$ACC PARALLEL LOOP COLLAPSE(force:5) PRESENT(WDES,CDIJ,W,W%CPROJ,CPROW,CQIJ,ENL,W%CELEN,W%FERWE) &
!$ACC PRIVATE(RTMP,EVALUE,WEIGHT,LBASE,LBASE_,NT,LMMAXC) __IF_ASYNC__
          spin: DO ISP=1,WDES%ISPIN

             band: DO N=1,WDES%NBANDS
                EVALUE=W%CELEN(N,NK,ISP)
                WEIGHT=WDES%WTKPT(NK)*W%FERWE(N,NK,ISP)*WDES%RSPIN

                spinor: DO ISPINOR=0,WDES%NRSPINORS-1
                   DO ISPINOR_=0,WDES%NRSPINORS-1
NOACC !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(NI,NT,LMMAXC,LBASE,LBASE_,L,LP,RTMP) &
NOACC !$OMP SHARED(WDES,CDIJ,CQIJ,NK,N,ISP,ISPINOR,ISPINOR_,EVALUE,W,CPROW,ENL,WEIGHT)
                      DO NI=1,WDES%NIONS
                         NT=WDES%ITYP(NI)
                         LMMAXC=WDES%LMMAX(NT)
                         IF (LMMAXC==0) CYCLE
                         LBASE =WDES%LMBASE(NI)+ISPINOR *WDES%NPRO/2
                         LBASE_=WDES%LMBASE(NI)+ISPINOR_*WDES%NPRO/2
                         RTMP=0
!$ACC LOOP COLLAPSE(2) REDUCTION(+:RTMP)
                         DO L=1 ,LMMAXC
                            DO LP=1,LMMAXC
                               RTMP=RTMP+WEIGHT*W%CPROJ(LBASE_+LP,N,NK,ISP)*GCONJG(CPROW(LBASE+L,N,NK,ISP))* &
                                    (CDIJ(LP,L,NI,ISP+ISPINOR_+2*ISPINOR)-EVALUE*CQIJ(LP,L,NI,ISP+ISPINOR_+2*ISPINOR))
                            ENDDO
                         ENDDO
!$ACC ATOMIC UPDATE
                         ENL(NI)=ENL(NI)+RTMP
                      ENDDO 
NOACC !$OMP END PARALLEL DO
                   ENDDO
                ENDDO spinor
             ENDDO band
          ENDDO spin
       ENDDO kpoint
!$ACC EXIT DATA COPYOUT(ENL) IF(ACC_EXEC_ON) WAIT(ACC_ASYNC_Q)
!$ACC EXIT DATA DELETE(WTMP) __IF_ASYNC__

       DO NI=1,WDES%NIONS
          NIP=NI_GLOBAL(NI, WDES%COMM_INB)
          FORNL(IDIR,NIP)=FORNL(IDIR,NIP)-ENL(NI)/DIS
       ENDDO
    ENDDO dir

    CALLMPI( M_sum_d(WDES%COMM, FORNL(1,1),NONLR_S%NIONS*3))
!$ACC EXIT DATA DELETE(CPROW) __IF_ASYNC__
    DEALLOCATE(CPROW)

    !TEST
    !      WRITE(*,'(4E20.12)') DIS
    !      WRITE(*,'("n",3F15.9 )') FORNL
    !      GOTO 1000
    !TEST

    PROFILING_STOP('fornlr')

    RETURN
  END SUBROUTINE FORNLR


!****************** SUBROUTINE RPROXYZ *********************************
!
! this subroutine calculates the first order change of the 
! wave function character upon moving the ions for one selected k-point
! and spin component
! the results are stored in CPROJXYZ
! mind that either the bra or the kat can vary, therefore
! a factor two has to be included
!
!***********************************************************************

  SUBROUTINE RPROXYZ(GRID, NONLR_S, P, LATT_CUR, W, WDES, ISP, NK, CPROJXYZ)
    USE pseudo
    USE mpimy
    USE lattice
    USE constant

    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    TYPE (potcar)      P(NONLR_S%NTYP)
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W,WTMP
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    INTEGER            ISP, NK
    GDEF :: CPROJXYZ(WDES%NPROD, WDES%NBANDS, 3)

    !-----some temporary arrays
    REAL(q) :: DISPL(3,NONLR_S%NIONS), DIS
    GDEF,ALLOCATABLE,TARGET :: CPROW(:,:,:,:)
    INTEGER :: IDIR

    PROFILING_START('rproxyz')

    ALLOCATE(CPROW(WDES%NPROD,WDES%NBANDS,WDES%NKPTS,WDES%ISPIN))
!$ACC ENTER DATA CREATE(CPROW) __IF_ASYNC__

!$ACC KERNELS PRESENT(CPROW) __IF_ASYNC__
    CPROW(:,:,NK,ISP)=0
!$ACC END KERNELS

!$ACC ENTER DATA CREATE(WTMP) __IF_ASYNC__
    WTMP=W
    WTMP%CPROJ => CPROW  ! relink the CPROJ array to temporary workspace

    DIS=fd_displacement

    DO IDIR=1,3
       ! operator= beta(r-R-dis)
       DISPL=0; DISPL(IDIR,:)= DIS
       ! this includes the factor two since a displacement + and - is performed by RSPHER_ALL
       !   2 d projector / d dis   = projector(R+dis) - projector(R-dis) / dis
       CALL RSPHER_ALL(GRID,NONLR_S,LATT_CUR,LATT_CUR,LATT_CUR,-DISPL,DISPL,1)

       CALL RPRO_ISP(NONLR_S,WDES,WTMP,GRID,ISP,NK)
!$ACC KERNELS PRESENT(CPROJXYZ,CPROW) __IF_ASYNC__
       CPROJXYZ(:,:,IDIR)=CPROW(:,:,NK,ISP)/DIS
!$ACC END KERNELS
    ENDDO

    DISPL=0
    CALL RSPHER_ALL(GRID,NONLR_S,LATT_CUR,LATT_CUR,LATT_CUR,DISPL,DISPL,0)

!$ACC EXIT DATA DELETE(WTMP,CPROW) __IF_ASYNC__
    DEALLOCATE(CPROW)

    PROFILING_STOP('rproxyz')

    RETURN
  END SUBROUTINE RPROXYZ


!****************** SUBROUTINE RPROLAT_DER *****************************
!
! this subroutine calculates the first order change of the 
! wave function character upon changing the lattice
!
!***********************************************************************

  SUBROUTINE RPROLAT_DER(GRID, NONLR_S, P, LATT_CUR, W, WDES, ISP, NK, CPROJXYZ)
    USE pseudo
    USE mpimy
    USE lattice
    USE constant

    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    TYPE (potcar)      P(NONLR_S%NTYP)
    TYPE (wavedes)     WDES
    TYPE (wavespin)    W,WTMP
    TYPE (grid_3d)     GRID
    TYPE (latt)        LATT_CUR
    INTEGER            ISP, NK
    GDEF :: CPROJXYZ(WDES%NPROD, WDES%NBANDS, 6)

    !-----some temporary arrays
    REAL(q) :: DISPL(3,NONLR_S%NIONS), DIS
    GDEF,ALLOCATABLE,TARGET :: CPROW(:,:,:,:)
    INTEGER :: IDIR, JDIR, IJDIR
    TYPE (latt)        LATT_FIN1,LATT_FIN2

    PROFILING_START('rprolat_der')

    ALLOCATE(CPROW(WDES%NPROD,WDES%NBANDS,WDES%NKPTS,WDES%ISPIN))
!$ACC ENTER DATA CREATE(CPROW) __IF_ASYNC__

!$ACC KERNELS PRESENT(CPROW) __IF_ASYNC__
    CPROW(:,:,NK,ISP)=0
!$ACC END KERNELS

!$ACC ENTER DATA CREATE(WTMP) __IF_ASYNC__
    WTMP=W
    WTMP%CPROJ => CPROW  ! relink the CPROJ array to temporary workspace

    DIS=fd_displacement

    IJDIR=0
    DO IDIR=1,3
       DO JDIR=1,IDIR
          IJDIR=IJDIR+1
          LATT_FIN1=LATT_CUR
          LATT_FIN1%A(IDIR,:)=LATT_CUR%A(IDIR,:)+DIS*LATT_CUR%A(JDIR,:)

          LATT_FIN2=LATT_CUR
          LATT_FIN2%A(IDIR,:)=LATT_CUR%A(IDIR,:)-DIS*LATT_CUR%A(JDIR,:)

          CALL LATTIC(LATT_FIN1)
          CALL LATTIC(LATT_FIN2)

          DISPL=0
          ! this includes the factor two since a displacement + and - is performed by RSPHER_ALL
          !   2 d projector / d dis   = projector(R+dis) - projector(R-dis) / dis
          CALL RSPHER_ALL(GRID,NONLR_S,LATT_FIN2,LATT_FIN1,LATT_CUR,DISPL,DISPL,1)

          CALL RPRO_ISP(NONLR_S,WDES,WTMP,GRID,ISP,NK)
!$ACC KERNELS PRESENT(CPROJXYZ,CPROW) __IF_ASYNC__
          CPROJXYZ(:,:,IJDIR)=CPROW(:,:,NK,ISP)/DIS
!$ACC END KERNELS
       ENDDO
    ENDDO

    DISPL=0
    CALL RSPHER_ALL(GRID,NONLR_S,LATT_CUR,LATT_CUR,LATT_CUR,DISPL,DISPL,0)

!$ACC EXIT DATA DELETE(WTMP,CPROW) __IF_ASYNC__
    DEALLOCATE(CPROW)

    PROFILING_STOP('rprolat_der')

    RETURN
  END SUBROUTINE RPROLAT_DER

END MODULE nonlr


!****************** SUBROUTINE RACC0 ********************************
!
! this subroutine calculates a linear combination of
! projection operatores in real space
! the result is added to CRACC
!
!> @details @ref openmp : each real space projection operator is
!> partitioned into openmp::omp_nthreads_nonlr_rspace parts.
!> The work on these parts is done by openmp::omp_nthreads_nonlr_rspace
!> threads. The work is distributed by means of an explicit loop
!> (OMP PARALLEL DO, label omp).
!
!*********************************************************************

  SUBROUTINE RACC0(NONLR_S,WDES1,CPROJ_LOC,CRACC)
    USE nonlr_struct_def
    USE wave
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)     WDES1
    COMPLEX(q) CRACC(WDES1%NRSPINORS*WDES1%GRID%MPLWV)
    GDEF       CPROJ_LOC(WDES1%NPROD)

  ! local
    REAL(q) RP
    INTEGER IP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NIS, NT, LMMAXC, NI, INDMAX, L, IND, LM
    COMPLEX(q)          :: CTMP
    REAL(qn),PARAMETER  :: ONE=1,ZERO=0
#ifdef avoidalloc
    REAL(qn) :: WORK(NONLR_S%IRMAX*2),TMP(101,2)
    GDEF    :: CPROJ(WDES1%NPRO_TOT)
#else 
    REAL(qn),ALLOCATABLE :: WORK(:),TMP(:,:)
    GDEF,ALLOCATABLE    :: CPROJ(:)
#endif
!$  INTEGER i,ndim
!$  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM
#ifdef memalign
!!DIR$ ASSUME_ALIGNED CRACC: __align_BTS
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: WORK,TMP
#endif

#ifndef avoidalloc
    ALLOCATE(WORK(NONLR_S%IRMAX*2),TMP(101,2),CPROJ(WDES1%NPRO_TOT))
#endif

    PROFILING_START('racc0')

    IF (WDES1%NK /= NONLR_S%NK) &
       CALL vtutor%bug("RACC0: PHASE not properly set up " // str(WDES1%NK) &
          // " " // str(NONLR_S%NK),__FILE__,__LINE__)

    ! merge projected wavefunctions from all nodes (if distributed over 
    !   plane wave coefficients)
    CALL MRG_PROJ(WDES1,CPROJ(1),CPROJ_LOC(1))

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

!$OMP PARALLEL DO NUM_THREADS(omp_nthreads_nonlr_rspace) &
!$OMP PRIVATE(i,WORK,TMP,LMBASE,ISPIRAL,ISPINOR,NLIIND,NIS,NT,LMMAXC,NI,INDMAX,L,CTMP,IND,LM,RP,IP) &
!$OMP SHARED(NONLR_S,WDES1,CPROJ,CRACC)
!$  omp: DO i=1,ndim
!=======================================================================
! loop over ions
!=======================================================================
    LMBASE= 0

    ISPIRAL = 1
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1

       NLIIND= 0
       NIS=1

       typ: DO NT=1,NONLR_S%NTYP
          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) GOTO 600

          ion: DO  NI=NIS,NONLR_S%NITYP(NT)+NIS-1
             INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
             IF (INDMAX == 0) GOTO 100
!=======================================================================
! set TMP
!=======================================================================
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(CTMP)
             DO L=1,LMMAXC
                CTMP= CPROJ(LMBASE+L)*WDES1%RINPL
                TMP(L,1)= REAL( CTMP ,KIND=qn)
#ifndef gammareal
                TMP(L,2)=AIMAG(CTMP)
#endif
             ENDDO
!=======================================================================
! calculate SUM(LM=1,NONLR_S%LMMAX) NONLR_S%RPROJ(K,LM) * TMP(LM)
!=======================================================================
#ifdef do_loops

             ! initialise accelerations to 0 (real &  imaginary part seperately)

!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
             DO IND=1,INDMAX
                WORK(IND)=0
#ifndef gammareal
                WORK(IND+NONLR_S%IRMAX)=0
#endif
             ENDDO

             DO LM=1,LMMAXC
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(RP)
                DO IND=1,INDMAX
                   RP  =  NONLR_S%RPROJ(IND+(LM-1)*INDMAX+NLIIND __NOACC_omp_arg(i))
                   WORK(IND)      =WORK(IND)      +TMP(LM,1)*RP
#ifndef gammareal
                   WORK(IND+NONLR_S%IRMAX)=WORK(IND+NONLR_S%IRMAX)+TMP(LM,2)*RP
#endif
                ENDDO
             ENDDO
#else
             CALL DGEMVn( 'N' , INDMAX, LMMAXC, ONE , NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), &
                  INDMAX, TMP(1,1) , 1 , ZERO , WORK(1), 1)
#ifndef gammareal
             CALL DGEMVn( 'N' , INDMAX, LMMAXC, ONE , NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), &
                  INDMAX, TMP(1,2) , 1 , ZERO , WORK(1+NONLR_S%IRMAX), 1)
#endif

#endif
!=======================================================================
!  add the non local contribution to the accelerations in real space
!=======================================================================
             IF (ASSOCIATED(NONLR_S%CRREXP)) THEN
#ifdef CRREXP_inline
!DIR$ IVDEP
!NEC$ list_vector
!OCL NOVREC
!$OMP SIMD PRIVATE(IP)
                DO IND=1,INDMAX
                   IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                   CRACC(IP)= CRACC(IP)+ &
                        GCMPLX( WORK(IND) , WORK(IND+NONLR_S%IRMAX) ,q) *CONJG(NONLR_S%CRREXP(IND,NI,ISPIRAL __NOACC_omp_arg(i)))
                ENDDO
#else
                CALL CRREXP_MUL_WORK_ADD(INDMAX,NONLR_S%CRREXP(1,NI,ISPIRAL __NOACC_omp_arg(i)),NONLR_S%NLI(1,NI __NOACC_omp_arg(i)), &
                          WORK(1),WORK(1+NONLR_S%IRMAX),CRACC(1+ISPINOR*WDES1%GRID%MPLWV))
#endif
             ELSE
!DIR$ IVDEP
!NEC$ list_vector
!OCL NOVREC
!$OMP SIMD PRIVATE(IP)
                DO IND=1,INDMAX
                   IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                   CRACC(IP)= CRACC(IP)+WORK(IND)
                ENDDO
             ENDIF

100          LMBASE= LMMAXC+LMBASE
             NLIIND= LMMAXC*INDMAX+NLIIND
          ENDDO ion
600       NIS = NIS+NONLR_S%NITYP(NT)
       ENDDO typ
       IF (NONLR_S%LSPIRAL) ISPIRAL=2
    ENDDO spinor
!$  ENDDO omp
!$OMP END PARALLEL DO

#ifndef avoidalloc
    DEALLOCATE(WORK,TMP,CPROJ)
#endif

    PROFILING_STOP('racc0')

    RETURN
  END SUBROUTINE RACC0


!****************** SUBROUTINE RACC0_HF ********************************
!
! this subroutine calculates a linear combination of
! projection operatores in real space
! the only difference is the CRACC is defined as GDEF whereas it is
! defined COMPLEX in the previous version
!
!> @details @ref openmp : each real space projection operator is
!> partitioned into openmp::omp_nthreads_nonlr_rspace parts.
!> The work on these parts is done by openmp::omp_nthreads_nonlr_rspace
!> threads. The work is distributed by means of an explicit loop
!> (OMP PARALLEL DO, label omp).
!
!*********************************************************************

  SUBROUTINE RACC0_HF(NONLR_S, WDES1, CPROJ_LOC, CRACC)
    USE nonlr_struct_def
    USE wave
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT NONE
    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)     WDES1
    GDEF   CRACC(WDES1%NRSPINORS*WDES1%GRID%MPLWV)
    GDEF   CPROJ_LOC(WDES1%NPROD)

    ! work array
    REAL(q) RP
    INTEGER IP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NIS, NT, LMMAXC, NI, INDMAX, L, IND, LM
    COMPLEX(q)          :: CTMP
    REAL(qn),PARAMETER :: ONE=1,ZERO=0
#ifdef avoidalloc
    REAL(qn) :: WORK(NONLR_S%IRMAX*2),TMP(101,2)
    GDEF     :: CPROJ(WDES1%NPRO_TOT)
#else
    REAL(qn),ALLOCATABLE :: WORK(:),TMP(:,:)
    GDEF,ALLOCATABLE     :: CPROJ(:)
#endif
!$  INTEGER i,ndim
!$  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM

#ifndef avoidalloc
    ALLOCATE(WORK(NONLR_S%IRMAX*2),TMP(101,2),CPROJ(WDES1%NPRO_TOT))
#endif

    PROFILING_START('racc0_hf')

    ! merge projected wavefunctions from all nodes (if distributed over 
    !   plane wave coefficients)
    CALL MRG_PROJ(WDES1,CPROJ(1),CPROJ_LOC(1))

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

!$OMP PARALLEL DO NUM_THREADS(omp_nthreads_nonlr_rspace) DEFAULT(NONE) &
!$OMP PRIVATE(i,WORK,TMP,LMBASE,ISPIRAL,ISPINOR,NLIIND,NIS,NT,LMMAXC,NI,INDMAX,L,CTMP,IND,LM,RP,IP) &
!$OMP SHARED(ndim,NONLR_S,WDES1,CPROJ,CRACC)
!$  omp: DO i=1,ndim
!=======================================================================
! loop over ions
!=======================================================================
    LMBASE= 0

    ISPIRAL = 1
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1

       NLIIND= 0
       NIS=1

       typ: DO NT=1,NONLR_S%NTYP
          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) GOTO 600

          ion: DO  NI=NIS,NONLR_S%NITYP(NT)+NIS-1
             INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
             IF (INDMAX == 0) GOTO 100
!=======================================================================
! set TMP
!=======================================================================
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(CTMP)
             DO L=1,LMMAXC
                CTMP= CPROJ(LMBASE+L)*WDES1%RINPL
                TMP(L,1)= REAL( CTMP ,KIND=qn)
#ifndef gammareal
                TMP(L,2)=AIMAG(CTMP)
#endif
             ENDDO
!=======================================================================
! calculate SUM(LM=1,NONLR_S%LMMAX) NONLR_S%RPROJ(K,LM) * TMP(LM)
!=======================================================================
#ifdef do_loops

             ! initialise accelerations to 0 (real &  imaginary part seperately)

!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
             DO IND=1,INDMAX
                WORK(IND)=0
#ifndef gammareal
                WORK(IND+NONLR_S%IRMAX)=0
#endif
             ENDDO

             DO LM=1,LMMAXC
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(RP)
                DO IND=1,INDMAX
                   RP  =  NONLR_S%RPROJ(IND+(LM-1)*INDMAX+NLIIND __NOACC_omp_arg(i))
                   WORK(IND)      =WORK(IND)      +TMP(LM,1)*RP
#ifndef gammareal
                   WORK(IND+NONLR_S%IRMAX)=WORK(IND+NONLR_S%IRMAX)+TMP(LM,2)*RP
#endif
                ENDDO
             ENDDO
#else
             CALL DGEMVn( 'N' , INDMAX, LMMAXC, ONE , NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), &
                  INDMAX, TMP(1,1) , 1 , ZERO , WORK(1), 1)
#ifndef gammareal
             CALL DGEMVn( 'N' , INDMAX, LMMAXC, ONE , NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), &
                  INDMAX, TMP(1,2) , 1 , ZERO , WORK(1+NONLR_S%IRMAX), 1)
#endif

#endif
!=======================================================================
!  add the non local contribution to the accelerations in real space
!=======================================================================
             IF (ASSOCIATED(NONLR_S%CRREXP)) THEN
#ifdef CRREXP_inline
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP)
                DO IND=1,INDMAX
                   IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                   CRACC(IP)= CRACC(IP)+ &
                        GCMPLX( WORK(IND) , WORK(IND+NONLR_S%IRMAX) ,q) *CONJG(NONLR_S%CRREXP(IND,NI,ISPIRAL __NOACC_omp_arg(i)))
                ENDDO
#else
                CALL CRREXP_MUL_WORK_GADD(INDMAX,NONLR_S%CRREXP(1,NI,ISPIRAL __NOACC_omp_arg(i)),NONLR_S%NLI(1,NI __NOACC_omp_arg(i)), &
                          WORK(1),WORK(1+NONLR_S%IRMAX),CRACC(1+ISPINOR*WDES1%GRID%MPLWV))
#endif
             ELSE
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP)
                DO IND=1,INDMAX
                   IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                   CRACC(IP)= CRACC(IP)+WORK(IND)
                ENDDO
             ENDIF


100          LMBASE= LMMAXC+LMBASE
             NLIIND= LMMAXC*INDMAX+NLIIND
          ENDDO ion
600       NIS = NIS+NONLR_S%NITYP(NT)
       ENDDO typ

       IF (NONLR_S%LSPIRAL) ISPIRAL=2
    ENDDO spinor
!$  ENDDO omp
!$OMP END PARALLEL DO

#ifndef avoidalloc
    DEALLOCATE(WORK,TMP,CPROJ)
#endif

    PROFILING_STOP('racc0_hf')

    RETURN
  END SUBROUTINE RACC0_HF


!****************** SUBROUTINE RACC0_REAL ******************************
!
! this subroutine calculates a linear combination of
! projection operatores in real space
! the only difference is the CRACC is defined as REAL whereas it is
! defined COMPLEX in RACC0 and GDEF in RACC0_HF
! this is required from the fast RPAR, GWR routines
!
! for details on OMP see previous routine
!
!*********************************************************************

  SUBROUTINE RACC0_REAL(NONLR_S, WDES1, CPROJ_LOC, CRACC)
    USE nonlr_struct_def
    USE wave
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT NONE
    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)     WDES1
    REAL(q) CRACC(WDES1%NRSPINORS*WDES1%GRID%MPLWV)
    GDEF    CPROJ_LOC(WDES1%NPROD)

    ! work array
    REAL(q) RP
    INTEGER IP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NIS, NT, LMMAXC, NI, INDMAX, L, IND, LM
    COMPLEX(q)          :: CTMP
    REAL(qn),PARAMETER :: ONE=1,ZERO=0
#ifdef avoidalloc
    REAL(qn) :: WORK(NONLR_S%IRMAX*2),TMP(101,2)
    GDEF     :: CPROJ(WDES1%NPRO_TOT)
#else
    REAL(qn),ALLOCATABLE :: WORK(:),TMP(:,:)
    GDEF,ALLOCATABLE     :: CPROJ(:)
#endif
!$  INTEGER i,ndim
!$  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM

#ifndef avoidalloc
    ALLOCATE(WORK(NONLR_S%IRMAX*2),TMP(101,2),CPROJ(WDES1%NPRO_TOT))
#endif

    PROFILING_START('racc0_real')

    ! merge projected wavefunctions from all nodes (if distributed over 
    !   plane wave coefficients)
    CALL MRG_PROJ(WDES1,CPROJ(1),CPROJ_LOC(1))

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

!$OMP PARALLEL DO NUM_THREADS(omp_nthreads_nonlr_rspace) DEFAULT(NONE) &
!$OMP PRIVATE(i,WORK,TMP,LMBASE,ISPIRAL,ISPINOR,NLIIND,NIS,NT,LMMAXC,NI,INDMAX,L,CTMP,IND,LM,RP,IP) &
!$OMP SHARED(ndim,NONLR_S,WDES1,CPROJ,CRACC)
!$  omp: DO i=1,ndim
!=======================================================================
! loop over ions
!=======================================================================
    LMBASE= 0

    ISPIRAL = 1
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1

       NLIIND= 0
       NIS=1

       typ: DO NT=1,NONLR_S%NTYP
          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) GOTO 600

          ion: DO  NI=NIS,NONLR_S%NITYP(NT)+NIS-1
             INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
             IF (INDMAX == 0) GOTO 100
!=======================================================================
! set TMP
!=======================================================================
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(CTMP)
             DO L=1,LMMAXC
                CTMP= CPROJ(LMBASE+L)*WDES1%RINPL
                TMP(L,1)= REAL( CTMP ,KIND=qn)
             ENDDO
!=======================================================================
! calculate SUM(LM=1,NONLR_S%LMMAX) NONLR_S%RPROJ(K,LM) * TMP(LM)
!=======================================================================
#ifdef do_loops
             ! initialise accelerations to 0 (real &  imaginary part seperately)
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
             DO IND=1,INDMAX
                WORK(IND)=0
             ENDDO

             DO LM=1,LMMAXC
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(RP)
                DO IND=1,INDMAX
                   RP  =  NONLR_S%RPROJ(IND+(LM-1)*INDMAX+NLIIND __NOACC_omp_arg(i))
                   WORK(IND)      =WORK(IND)      +TMP(LM,1)*RP
                ENDDO
             ENDDO
#else
             CALL DGEMVn( 'N' , INDMAX, LMMAXC, ONE , NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), &
                  INDMAX, TMP(1,1) , 1 , ZERO , WORK(1), 1)

#endif
!=======================================================================
!  add the non local contribution to the accelerations in real space
!=======================================================================
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP)
             DO IND=1,INDMAX
                IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                CRACC(IP)= CRACC(IP)+WORK(IND)
             ENDDO

100          LMBASE= LMMAXC+LMBASE
             NLIIND= LMMAXC*INDMAX+NLIIND
          ENDDO ion
600       NIS = NIS+NONLR_S%NITYP(NT)
       ENDDO typ

       IF (NONLR_S%LSPIRAL) ISPIRAL=2
    ENDDO spinor
!$  ENDDO omp
!$OMP END PARALLEL DO

#ifndef avoidalloc
    DEALLOCATE(WORK,TMP,CPROJ)
#endif

    PROFILING_STOP('racc0_real')

    RETURN
  END SUBROUTINE RACC0_REAL


!****************** SUBROUTINE RACC0MU *******************************
!
! this subroutine calculates a set of linear combination of
! projection operatores in real space
! the result is added to CRACC
!               -----
!
!> @details @ref openmp : each real space projection operator is
!> partitioned into openmp::omp_nthreads_nonlr_rspace parts.
!> The work on these parts is done by openmp::omp_nthreads_nonlr_rspace
!> threads. The work is distributed by means of an explicit loop
!> (OMP PARALLEL DO, label omp).
!
!*********************************************************************

  SUBROUTINE RACC0MU(NONLR_S, WDES1, CPROJ_LOC, CRACC, LD, NSIM, LDO)
    USE nonlr_struct_def
    USE wave
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)     WDES1
    INTEGER LD                        ! leading dimension of CRACC
    INTEGER NSIM                      ! do NSIM bands at a time
    COMPLEX(q) CRACC(LD,NSIM)         ! result in real space
    GDEF   CPROJ_LOC(WDES1%NPROD,NSIM)! wave function character
    LOGICAL LDO(NSIM)                 ! which bands are included

  ! local
    INTEGER, PARAMETER  :: NLM=__align_DIM(101)
    REAL(qn),PARAMETER  :: ONE=1,ZERO=0
    INTEGER NP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NT, LMMAXC, NI, &
         INDMAX, IND0, NPFILL, L, IND, IP
    COMPLEX(q)          :: CTMP

!$  INTEGER i,ndim
!$  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS

#ifdef avoidalloc
    GDEF     :: CPROJ(WDES1%NPRO_TOT,NSIM)
    REAL(qn) :: WORK(ndata*NSIM*NONLR_S%IRMAX),TMP(NLM,ndata*2*NSIM)
#else
    REAL(qn),ALLOCATABLE:: WORK(:),TMP(:,:)
    GDEF,ALLOCATABLE    :: CPROJ(:,:)
#endif
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: WORK,TMP
#endif

    PROFILING_START('racc0mu')

#ifndef avoidalloc
    ALLOCATE(CPROJ(WDES1%NPRO_TOT,NSIM),WORK(ndata*NSIM*NONLR_S%IRMAX),TMP(NLM,ndata*2*NSIM))
#endif

    IF (WDES1%NK /= NONLR_S%NK) &
       CALL vtutor%bug("RACC0MU: PHASE not properly set up " // str(WDES1%NK) &
          // " " // str(NONLR_S%NK),__FILE__,__LINE__)

    ! merge projected wavefunctions from all nodes
    DO NP=1,NSIM
       IF (LDO(NP)) THEN
          CALL MRG_PROJ(WDES1,CPROJ(1,NP),CPROJ_LOC(1,NP))
       ENDIF
    ENDDO

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

!$OMP PARALLEL DO NUM_THREADS(omp_nthreads_nonlr_rspace) DEFAULT(NONE) &
!$OMP PRIVATE(i,WORK,TMP,LMBASE,ISPIRAL,ISPINOR,NLIIND,NT,LMMAXC,NI,INDMAX,IND0,NPFILL,NP,L,CTMP,IND,IP) &
!$OMP SHARED(ndim,NONLR_S,WDES1,CPROJ,CRACC,NSIM,LDO)
!$  omp: DO i=1,ndim
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1

       ISPIRAL=1; IF (NONLR_S%LSPIRAL) ISPIRAL=ISPINOR+1

!=======================================================================
! loop over ions
!=======================================================================

       ion: DO NI=1,NONLR_S%NIONS
          NT=NONLR_S%ITYP(NI)

          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) CYCLE ion

          INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
          IF (INDMAX==0) CYCLE ion

          LMBASE=NONLR_S%LMBASE(NI)+ISPINOR*NONLR_S%LMBASE(NONLR_S%NIONS+1)
          NLIIND=NONLR_S%NLIBASE(NI __NOACC_omp_arg(i))

!=======================================================================
! set TMP
!=======================================================================
          IND0=0
          NPFILL=0
          DO NP=1,NSIM
             IF (LDO(NP)) THEN
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(CTMP)
                DO L=1,LMMAXC
                   CTMP= CPROJ(LMBASE+L,NP)*WDES1%RINPL
                   TMP(L,1+IND0)= REAL( CTMP ,KIND=qn)
#ifndef gammareal
                   TMP(L,2+IND0)=AIMAG(CTMP)
#endif
                ENDDO
                IND0=IND0+ndata
                NPFILL=NPFILL+1
             ENDIF

          ENDDO
!=======================================================================
! calculate SUM(LM=1,NONLR_S%LMMAX) NONLR_S%RPROJ(K,LM) * TMP(LM)
!=======================================================================
          PROFILING_START('gemm')
#ifdef RACCMU_DGEMV
          DO NP=0,ndata*NPFILL-1
             CALL DGEMVn( 'N' , INDMAX, LMMAXC,  ONE , NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), &
                  INDMAX,  TMP(1,1+NP) , 1 , ZERO ,  WORK(1+NONLR_S%IRMAX*NP), 1)
          ENDDO
#else
          CALL DGEMMn( 'N' , 'N', INDMAX, ndata*NPFILL, LMMAXC, ONE, &
               NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), INDMAX, TMP(1,1) , NLM , &
               ZERO , WORK(1), NONLR_S%IRMAX)
#endif
          PROFILING_STOP('gemm')
!=======================================================================
!  add the non local contribution to the accelerations
!=======================================================================
          IND0=0
          DO NP=1,NSIM
             IF (LDO(NP)) THEN
                IF (ASSOCIATED(NONLR_S%CRREXP)) THEN
#ifdef CRREXP_inline
!DIR$ IVDEP
!NEC$ list_vector
!OCL NOVREC
!$OMP SIMD PRIVATE(IP)
                   DO IND=1,INDMAX
                      IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                      CRACC(IP,NP)= CRACC(IP,NP)+ &
                           GCMPLX( WORK(IND+IND0), WORK(IND+(NONLR_S%IRMAX+IND0)) ,q)* &
                           CONJG(NONLR_S%CRREXP(IND,NI,ISPIRAL __NOACC_omp_arg(i) ))
                   ENDDO
#else
                   CALL CRREXP_MUL_WORK_ADD(INDMAX,NONLR_S%CRREXP(1,NI,ISPIRAL __NOACC_omp_arg(i)),NONLR_S%NLI(1,NI __NOACC_omp_arg(i)), &
                       WORK(1+IND0),WORK(1+IND0+NONLR_S%IRMAX),CRACC(1+ISPINOR*WDES1%GRID%MPLWV,NP))
#endif
                ELSE
!DIR$ IVDEP
!NEC$ list_vector
!OCL NOVREC
!$OMP SIMD PRIVATE(IP)
                   DO IND=1,INDMAX
                      IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                      CRACC(IP,NP)= CRACC(IP,NP)+WORK(IND+IND0)
                   ENDDO
                ENDIF

                IND0=IND0+ndata * NONLR_S%IRMAX
             ENDIF
          ENDDO
       ENDDO ion

    ENDDO spinor
!$  ENDDO omp
!$OMP END PARALLEL DO

#ifndef avoidalloc
    DEALLOCATE(CPROJ,WORK,TMP)
#endif

    PROFILING_STOP('racc0mu')

    RETURN
  END SUBROUTINE RACC0MU

#ifdef _OPENACC
  SUBROUTINE RACC0_ACC(NONLR_S, WDES1, CPROJ_LOC, CRACC)
    USE mopenacc_struct_def
    USE nonlr_struct_def
    USE wave
    USE openmp, ONLY : omp_nthreads_acc
    IMPLICIT NONE

!$ACC ROUTINE(SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX) SEQ

    TYPE (nonlr_struct) :: NONLR_S
    TYPE (wavedes1)     :: WDES1
    COMPLEX(q)          :: CRACC(WDES1%NRSPINORS*WDES1%GRID%MPLWV)
    GDEF                :: CPROJ_LOC(WDES1%NPROD) ! wave function character

  ! local
    INTEGER :: NP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NT, LMMAXC, NI, INDMAX, L, IND, IP
    COMPLEX(q) :: CTMP

    INTEGER :: BLOCKSIZE=128*1024
    INTEGER :: IBLOCK

!$  INTEGER i,ndim

#ifdef avoidalloc
    GDEF :: CPROJ(WDES1%NPRO_TOT)
#else
    GDEF, ALLOCATABLE :: CPROJ(:)
#endif

    PROFILING_START('racc0_acc')

#ifndef avoidalloc
    ALLOCATE(CPROJ(WDES1%NPRO_TOT))
#endif

!$ACC ENTER DATA CREATE(CPROJ) ASYNC(ACC_ASYNC_Q)

    IF (WDES1%NK /= NONLR_S%NK) &
       CALL vtutor%bug("RACC0_ACC: PHASE not properly set up " // str(WDES1%NK) &
          // " " // str(NONLR_S%NK),__FILE__,__LINE__)

    ! merge projected wavefunctions from all nodes
    CALL MRG_PROJ(WDES1,CPROJ(1),CPROJ_LOC(1))

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

!$  omp: DO i=1,ndim

!$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(CPROJ,WDES1,WDES1%GRID,NONLR_S,CRACC) &
!$ACC PRIVATE(ISPIRAL,NT,LMMAXC,INDMAX,LMBASE,NLIIND,IBLOCK) ASYNC(ACC_ASYNC_Q)
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
!=======================================================================
! loop over ions
!=======================================================================

       ion: DO NI=1,NONLR_S%NIONS
          NT=NONLR_S%ITYP(NI)

          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) CYCLE ion

          INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
          IF (INDMAX==0) CYCLE ion

          LMBASE=NONLR_S%LMBASE(NI)+ISPINOR*NONLR_S%LMBASE(NONLR_S%NIONS+1)
          NLIIND=NONLR_S%NLIBASE(NI __NOACC_omp_arg(i))
          ISPIRAL=1; IF (NONLR_S%LSPIRAL) ISPIRAL=ISPINOR+1

          DO IBLOCK=0,INDMAX/BLOCKSIZE
!$ACC LOOP VECTOR PRIVATE(IP,CTMP)
             DO IND=IBLOCK*BLOCKSIZE+1,MIN((IBLOCK+1)*BLOCKSIZE,INDMAX)

                CTMP=0
!$ACC LOOP REDUCTION(+:CTMP)
                DO L=1,LMMAXC
!$ACC CACHE(CPROJ(LMBASE:LMBASE+LMMAXC))
                   CTMP=CTMP+CPROJ(L+LMBASE)*NONLR_S%RPROJ(IND+(L-1)*INDMAX+NLIIND __NOACC_omp_arg(i))
                ENDDO
#ifndef gammareal
                CTMP=CTMP*CONJG(NONLR_S%CRREXP(IND,NI,ISPIRAL __NOACC_omp_arg(i)))
#endif
                IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                CALL SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX(CRACC(IP),CTMP*WDES1%RINPL)
             ENDDO
          ENDDO
       ENDDO ion

    ENDDO spinor
!$  ENDDO omp
!$ACC EXIT DATA DELETE(CPROJ) ASYNC(ACC_ASYNC_Q)

#ifndef avoidalloc
    DEALLOCATE(CPROJ)
#endif

    PROFILING_STOP('racc0_acc')

    RETURN
  END SUBROUTINE RACC0_ACC


  SUBROUTINE RACC0MU_ACC(NONLR_S, WDES1, CPROJ_LOC, CRACC, LD, NSIM, LDO)
    USE mopenacc
    USE nonlr_struct_def
    USE wave
    USE openmp, ONLY : omp_nthreads_acc
    IMPLICIT NONE

!$ACC ROUTINE(SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX) SEQ

    TYPE (nonlr_struct) :: NONLR_S
    TYPE (wavedes1)     :: WDES1
    INTEGER             :: LD                          ! leading dimension of CRACC
    INTEGER             :: NSIM                        ! do NSIM bands at a time
    COMPLEX(q)          :: CRACC(LD,NSIM)              ! result in real space
    GDEF                :: CPROJ_LOC(WDES1%NPROD,NSIM) ! wave function character
    LOGICAL             :: LDO(NSIM)                   ! which bands are included

  ! local
    INTEGER :: NP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NT, LMMAXC, NI, INDMAX, L, IND, IP
    COMPLEX(q) :: CTMP

    INTEGER :: BLOCKSIZE=128*1024
    INTEGER :: IBLOCK

!$  INTEGER i,ndim

#ifdef avoidalloc
    GDEF :: CPROJ(WDES1%NPRO_TOT,NSIM)
#else
    GDEF, ALLOCATABLE :: CPROJ(:,:)
#endif

    PROFILING_START('racc0mu_acc')

#ifndef avoidalloc
    ALLOCATE(CPROJ(WDES1%NPRO_TOT,NSIM))
#endif

!$ACC ENTER DATA CREATE(CPROJ)

    IF (WDES1%NK /= NONLR_S%NK) &
       CALL vtutor%bug("RACC0MU: PHASE not properly set up " // str(WDES1%NK) &
          // " " // str(NONLR_S%NK),__FILE__,__LINE__)

    ! merge projected wavefunctions from all nodes
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(NP) &
!$OMP NUM_THREADS(omp_nthreads_acc)
    DO NP=1,NSIM
       IF (LDO(NP)) THEN
          SET_ACC_ASYNC_Q(NP)
          CALL MRG_PROJ(WDES1,CPROJ(1,NP),CPROJ_LOC(1,NP))
       ENDIF
    ENDDO
!$OMP END PARALLEL DO

    CALL ACC_SYNC_ASYNC_Q(NSIM,LDO)
    ACC_ASYNC_Q=ACC_ASYNC_ASYNC

!$ACC ENTER DATA COPYIN(LDO) __IF_ASYNC__

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

!$  omp: DO i=1,ndim

!$ACC PARALLEL LOOP COLLAPSE(3) GANG PRESENT(CPROJ,WDES1,WDES1%GRID,NONLR_S,CRACC) &
!$ACC PRIVATE(ISPIRAL,NT,LMMAXC,INDMAX,LMBASE,NLIIND,IBLOCK) __IF_ASYNC__
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
!=======================================================================
! loop over ions
!=======================================================================

       ion: DO NI=1,NONLR_S%NIONS
          DO NP=1,NSIM
             IF (.NOT.LDO(NP)) CYCLE
             NT=NONLR_S%ITYP(NI)

             LMMAXC=NONLR_S%LMMAX(NT)
             IF (LMMAXC==0) CYCLE !ion

             INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
             IF (INDMAX==0) CYCLE !ion

             LMBASE=NONLR_S%LMBASE(NI)+ISPINOR*NONLR_S%LMBASE(NONLR_S%NIONS+1)
             NLIIND=NONLR_S%NLIBASE(NI __NOACC_omp_arg(i))
             ISPIRAL=1; IF (NONLR_S%LSPIRAL) ISPIRAL=ISPINOR+1

             DO IBLOCK=0,INDMAX/BLOCKSIZE
!$ACC LOOP VECTOR PRIVATE(IP,CTMP)
             DO IND=IBLOCK*BLOCKSIZE+1,MIN((IBLOCK+1)*BLOCKSIZE,INDMAX)

                CTMP=0
!$ACC LOOP REDUCTION(+:CTMP)
                DO L=1,LMMAXC
!$ACC CACHE(CPROJ(1+LMBASE:LMBASE+LMMAXC,NP))
                   CTMP=CTMP+CPROJ(L+LMBASE,NP)*NONLR_S%RPROJ(IND+(L-1)*INDMAX+NLIIND __NOACC_omp_arg(i))
                ENDDO
#ifndef gammareal
                CTMP=CTMP*CONJG(NONLR_S%CRREXP(IND,NI,ISPIRAL __NOACC_omp_arg(i)))
#endif
                IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                CALL SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX(CRACC(IP,NP),CTMP*WDES1%RINPL)
             ENDDO
             ENDDO
          ENDDO
       ENDDO ion

    ENDDO spinor
!$  ENDDO omp
!$ACC EXIT DATA DELETE(LDO,CPROJ) __IF_ASYNC__

#ifndef avoidalloc
    DEALLOCATE(CPROJ)
#endif

    PROFILING_STOP('racc0mu_acc')

    RETURN
  END SUBROUTINE RACC0MU_ACC


  SUBROUTINE RACC0MU_HF_ACC(NONLR_S, WDES1, CPROJ_LOC, LD1, CRACC, LD2, NSIM)
    USE mopenacc_struct_def
    USE nonlr_struct_def
    USE wave
    USE openmp, ONLY : omp_nthreads_acc
    IMPLICIT NONE

!$ACC ROUTINE(SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX) SEQ

    TYPE (nonlr_struct) :: NONLR_S
    TYPE (wavedes1)     :: WDES1
    GDEF                :: CPROJ_LOC(LD1,NSIM)         ! wave function character
    GDEF                :: CRACC(LD2,NSIM)             ! result in real space
    INTEGER             :: LD1,LD2                     ! leading dimension of CRACC
    INTEGER             :: NSIM                        ! do NSIM bands at a time
 
  ! local
    INTEGER :: NP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NT, LMMAXC, NI, INDMAX, L, IND, IP
    GDEF :: CTMP

    COMPLEX(q) :: CRREXP

!$  INTEGER i,ndim

#ifdef avoidalloc
    GDEF :: CPROJ(WDES1%NPRO_TOT,NSIM)
#else
    GDEF, ALLOCATABLE :: CPROJ(:,:)
#endif

    PROFILING_START('racc0mu_hf_acc')

#ifndef avoidalloc
    ALLOCATE(CPROJ(WDES1%NPRO_TOT,NSIM))
#endif

!$ACC ENTER DATA CREATE(CPROJ) ASYNC(ACC_ASYNC_Q)

#ifdef MPI
    IF (WDES1%COMM_INB%NCPU>1) THEN
    ! merge projected wavefunctions from all nodes
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(NP) NUM_THREADS(omp_nthreads_acc)
       DO NP=1,NSIM
          CALL MRG_PROJ(WDES1,CPROJ(1,NP),CPROJ_LOC(1,NP))
       ENDDO
!$OMP END PARALLEL DO
    ELSE
#endif
    ! in this case a simple copy suffices
!$ACC KERNELS PRESENT(WDES1,CPROJ,CPROJ_LOC) __IF_ASYNC__
         CPROJ(1:WDES1%NPRO,1:NSIM)=CPROJ_LOC(1:WDES1%NPRO,1:NSIM)
!$ACC END KERNELS
#ifdef MPI
    ENDIF
#endif

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

!$  omp: DO i=1,ndim

    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
       ISPIRAL=1; IF (NONLR_S%LSPIRAL) ISPIRAL=ISPINOR+1
!$ACC PARALLEL LOOP COLLAPSE(2) GANG VECTOR PRESENT(CPROJ,WDES1,WDES1%GRID,NONLR_S,CRACC) &
!$ACC PRIVATE(NT,LMMAXC,INDMAX,LMBASE,NLIIND,IP,CRREXP) ASYNC(ACC_ASYNC_Q)
!=======================================================================
! loop over ions
!=======================================================================
       ion: DO NI=1,NONLR_S%NIONS
          DO IND=1,NONLR_S%IRMAX

             NT=NONLR_S%ITYP(NI); LMMAXC=NONLR_S%LMMAX(NT)
             INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))

             IF ( IND > INDMAX .OR. LMMAXC == 0 .OR. INDMAX == 0 ) CYCLE !ion

             LMBASE=NONLR_S%LMBASE(NI)+ISPINOR*NONLR_S%LMBASE(NONLR_S%NIONS+1)
             NLIIND=NONLR_S%NLIBASE(NI __NOACC_omp_arg(i))+IND

             IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
#ifndef gammareal
             CRREXP=CONJG(NONLR_S%CRREXP(IND,NI,ISPIRAL __NOACC_omp_arg(i)))
#endif

!$ACC LOOP SEQ PRIVATE(CTMP)
             DO NP=1,NSIM
                CTMP=0
!$ACC LOOP SEQ
                DO L=1,LMMAXC
                   CTMP=CTMP+CPROJ(L+LMBASE,NP)*NONLR_S%RPROJ((L-1)*INDMAX+NLIIND __NOACC_omp_arg(i))
                ENDDO
!NOTE(sm): Original (CPU) version multiplies CTMP with WDES1%RINPL here. The only call
!          call site for RACC0MU_HF is FOCK_CHARGE_MU, so we just remove
!          the multiplication there and spare the ACC UPDATE there.
#ifndef gammareal
                CALL SPLIT_CMPLX_ATOMIC_ADD_FROM_CMPLX(CRACC(IP,NP),CTMP*CRREXP) !*WDES1%RINPL)
#else
!$ACC ATOMIC UPDATE
                CRACC(IP,NP)=CRACC(IP,NP) + CTMP !*WDES1%RINPL
#endif
             ENDDO
          ENDDO
       ENDDO ion

    ENDDO spinor
!$  ENDDO omp
!$ACC EXIT DATA DELETE(CPROJ) ASYNC(ACC_ASYNC_Q)

#ifndef avoidalloc
    DEALLOCATE(CPROJ)
#endif

    PROFILING_STOP('racc0mu_hf_acc')

    RETURN
  END SUBROUTINE RACC0MU_HF_ACC
#endif

!****************** SUBROUTINE RACC0MU_HF ****************************
!
! exact copy of the previous routine with CRACC defined as GDEF
! phase factor test is also removed
!
!> @details @ref openmp : each real space projection operator is
!> partitioned into openmp::omp_nthreads_nonlr_rspace parts.
!> The work on these parts is done by openmp::omp_nthreads_nonlr_rspace
!> threads. The work is distributed by means of an explicit loop
!> (OMP PARALLEL DO, label omp).
!
!*********************************************************************

  SUBROUTINE RACC0MU_HF(NONLR_S, WDES1, CPROJ_LOC, LD1, CRACC, LD2, NSIM)
    USE nonlr_struct_def
    USE wave
    USE openmp, ONLY : omp_nthreads_nonlr_rspace
    IMPLICIT NONE

    TYPE (nonlr_struct) NONLR_S
    TYPE (wavedes1)     WDES1
    INTEGER LD1,LD2                   ! leading dimension of CPROJ_LOC and CRACC
    INTEGER NSIM                      ! do NSIM bands at a time
    GDEF   CPROJ_LOC(LD1,NSIM)        ! wave function character
    GDEF   CRACC(LD2,NSIM)            ! result in real space

  ! local
    INTEGER, PARAMETER  :: NLM=__align_DIM(101)
    REAL(qn),PARAMETER  :: ONE=1,ZERO=0
    INTEGER NP, LMBASE, ISPIRAL, ISPINOR, NLIIND, NIS, NT, LMMAXC, NI, &
         INDMAX, IND0, NPFILL, L, IND, IP
    COMPLEX(q)          :: CTMP

!$  INTEGER i,ndim
!$  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS

#ifdef avoidalloc
    GDEF    :: CPROJ(WDES1%NPRO_TOT,NSIM)
    REAL(qn) :: WORK(ndata*NSIM*NONLR_S%IRMAX),TMP(NLM,ndata*2*NSIM)
#else
    REAL(qn),ALLOCATABLE:: WORK(:),TMP(:,:)
    GDEF,ALLOCATABLE    :: CPROJ(:,:)
#endif
#ifdef memalign
!DIR$ ATTRIBUTES ALIGN: __align_BTS :: WORK,TMP
#endif

    PROFILING_START('racc0mu_hf')

#ifndef avoidalloc
    ALLOCATE(CPROJ(WDES1%NPRO_TOT,NSIM),WORK(ndata*NSIM*NONLR_S%IRMAX),TMP(NLM,ndata*2*NSIM))
#endif

!   PROFILING_START('racc0mu_hf_dis')

    ! merge projected wavefunctions from all nodes
    DO NP=1,NSIM
       CALL MRG_PROJ(WDES1,CPROJ(1,NP),CPROJ_LOC(1,NP))
    ENDDO

!   PROFILING_STOP('racc0mu_hf_dis')

NOACC !$  ndim=SIZE(NONLR_S%NLIMAX,2)
DOACC !$  ndim=1

!$OMP PARALLEL DO NUM_THREADS(omp_nthreads_nonlr_rspace) &
!$OMP PRIVATE(i,WORK,TMP,LMBASE,ISPIRAL,ISPINOR,NLIIND,NIS,NT,LMMAXC,NI,INDMAX,IND0,NPFILL,NP,L,CTMP,IND,IP) &
!$OMP SHARED(NONLR_S,WDES1,CPROJ,CRACC,NSIM)
!$  omp: DO i=1,ndim
!=======================================================================
! loop over ions
!=======================================================================
    LMBASE= 0

    ISPIRAL = 1
    spinor: DO ISPINOR=0,WDES1%NRSPINORS-1
#ifndef _OPENMP
       NLIIND= 0
       NIS=1

       typ: DO NT=1,NONLR_S%NTYP
          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) GOTO 600

          ion: DO NI=NIS,NONLR_S%NITYP(NT)+NIS-1
             INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
             IF (INDMAX == 0) GOTO 100
#else
!#ifndef avoid_reduction
!!$OMP PARALLEL DO &
!!$OMP DEFAULT(SHARED) &
!!$OMP PRIVATE(NI,NT,LMMAXC,INDMAX,LMBASE,NLIIND,IND0,NPFILL,L,CTMP,TMP,NP,WORK,IND,IP) &
!!$OMP REDUCTION(+:CRACC)
!#endif
       ion: DO NI=1,NONLR_S%NIONS
          NT=NONLR_S%ITYP(NI)

          LMMAXC=NONLR_S%LMMAX(NT)
          IF (LMMAXC==0) CYCLE ion

          INDMAX=NONLR_S%NLIMAX(NI __NOACC_omp_arg(i))
          IF (INDMAX==0) CYCLE ion

          LMBASE=NONLR_S%LMBASE(NI)+ISPINOR*NONLR_S%LMBASE(NONLR_S%NIONS+1)
          NLIIND=NONLR_S%NLIBASE(NI __NOACC_omp_arg(i))
#endif
!=======================================================================
! set TMP
!=======================================================================
!            PROFILING_START('racc0mu_hf_gather')

             IND0=0
             NPFILL=0
             DO NP=1,NSIM
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(CTMP)
                   DO L=1,LMMAXC
                      CTMP= CPROJ(LMBASE+L,NP)*WDES1%RINPL
                      TMP(L,1+IND0)= REAL( CTMP ,KIND=qn)
#ifndef gammareal
                      TMP(L,2+IND0)=AIMAG(CTMP)
#endif
                   ENDDO
                   IND0=IND0+ndata
                   NPFILL=NPFILL+1
             ENDDO

!            PROFILING_STOP('racc0mu_hf_gather')
!=======================================================================
! calculate SUM(LM=1,NONLR_S%LMMAX) NONLR_S%RPROJ(K,LM) * TMP(LM)
!=======================================================================
!            PROFILING_START('racc0mu_hf_gemm')
#ifdef RACCMU_DGEMV
             DO NP=0,ndata*NPFILL-1
                CALL DGEMVn( 'N' , INDMAX, LMMAXC,  ONE , NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), &
                     INDMAX,  TMP(1,1+NP) , 1 , ZERO ,  WORK(1+NONLR_S%IRMAX*NP), 1)
             ENDDO
#else
             CALL DGEMMn( 'N' , 'N', INDMAX, ndata*NPFILL, LMMAXC, ONE, &
                  NONLR_S%RPROJ(1+NLIIND __NOACC_omp_arg(i)), INDMAX, TMP(1,1) , NLM , &
                  ZERO , WORK(1), NONLR_S%IRMAX)
#endif
!            PROFILING_STOP('racc0mu_hf_gemm')
!=======================================================================
!  add the non local contribution to the accelerations
!=======================================================================
!            PROFILING_START('racc0mu_hf_scatter')

             IND0=0
             DO NP=1,NSIM
                   IF (ASSOCIATED(NONLR_S%CRREXP)) THEN
#ifdef CRREXP_inline
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP)
                      DO IND=1,INDMAX
                         IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                         CRACC(IP,NP)= CRACC(IP,NP)+ &
                              GCMPLX( WORK(IND+IND0), WORK(IND+(NONLR_S%IRMAX+IND0)) ,q)* &
                              CONJG(NONLR_S%CRREXP(IND,NI,ISPIRAL __NOACC_omp_arg(i)))
                      ENDDO
#else
                      CALL CRREXP_MUL_WORK_GADD(INDMAX,NONLR_S%CRREXP(1,NI,ISPIRAL __NOACC_omp_arg(i)),NONLR_S%NLI(1,NI __NOACC_omp_arg(i)), &
                          WORK(1+IND0),WORK(1+IND0+NONLR_S%IRMAX),CRACC(1+ISPINOR*WDES1%GRID%MPLWV,NP))
#endif
                   ELSE
!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP)
                      DO IND=1,INDMAX
                         IP=NONLR_S%NLI(IND,NI __NOACC_omp_arg(i))+ISPINOR*WDES1%GRID%MPLWV
                         CRACC(IP,NP)= CRACC(IP,NP)+WORK(IND+IND0)
                      ENDDO
                   ENDIF

                   IND0=IND0+ndata * NONLR_S%IRMAX
             ENDDO

!            PROFILING_STOP('racc0mu_hf_scatter')
#ifndef _OPENMP
100          LMBASE= LMMAXC+LMBASE
             NLIIND= LMMAXC*INDMAX+NLIIND
          ENDDO ion

600       NIS = NIS+NONLR_S%NITYP(NT)
       ENDDO typ
#else
       ENDDO ion
!#ifndef avoid_reduction
!!$OMP END PARALLEL DO
!#endif
#endif
       IF (NONLR_S%LSPIRAL) ISPIRAL=2
    ENDDO spinor
!$  ENDDO omp
!$OMP END PARALLEL DO
#ifndef avoidalloc
    DEALLOCATE(CPROJ,WORK,TMP)
#endif

    PROFILING_STOP('racc0mu_hf')

    RETURN
  END SUBROUTINE RACC0MU_HF


!***********************************************************************
!
! small f77 helper  routines to 
! multiply with phasefactor and divide into real and imaginary part
!
!***********************************************************************
  
  SUBROUTINE CRREXP_MUL_WAVE( INDMAX, CRREXP, NLI, CR, WORK1, WORK2)
    USE prec
    IMPLICIT NONE
    INTEGER INDMAX
    COMPLEX(q) :: CRREXP(INDMAX)
    INTEGER    :: NLI(INDMAX)
    COMPLEX(q) :: CR(*)
    REAL(qn)   :: WORK1(INDMAX), WORK2(INDMAX)
    COMPLEX(qn):: CTMP
  ! local
    INTEGER IND, IP

    PROFILING_START('crrexp_mul_wave')

!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP,CTMP)
    DO IND=1,INDMAX
       IP=NLI(IND)
       CTMP=    CR(IP)*CRREXP(IND)
       WORK1(IND) = REAL( CTMP ,KIND=q)
       WORK2(IND)=  AIMAG(CTMP)
    ENDDO

    PROFILING_STOP('crrexp_mul_wave')

  END SUBROUTINE CRREXP_MUL_WAVE


  SUBROUTINE CRREXP_MUL_GWAVE( INDMAX, CRREXP, NLI, CR, WORK1, WORK2)
    USE prec
    IMPLICIT NONE
    INTEGER INDMAX
    COMPLEX(q) :: CRREXP(INDMAX)
    INTEGER    :: NLI(INDMAX)
    GDEF       :: CR(*)
    REAL(qn)   :: WORK1(INDMAX), WORK2(INDMAX)
    COMPLEX(qn):: CTMP
  ! local
    INTEGER IND, IP

!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP,CTMP)
    DO IND=1,INDMAX
       IP=NLI(IND)
       CTMP=    CR(IP)*CRREXP(IND)
       WORK1(IND) = REAL( CTMP ,KIND=q)
       WORK2(IND)=  AIMAG(CTMP)
    ENDDO
  END SUBROUTINE CRREXP_MUL_GWAVE


  SUBROUTINE CRREXP_MUL_WORK_ADD( INDMAX, CRREXP, NLI, WORK1, WORK2, CR)
    USE prec
    IMPLICIT NONE
    INTEGER INDMAX
    COMPLEX(q) :: CRREXP(INDMAX)
    INTEGER    :: NLI(INDMAX)
    REAL(qn)   :: WORK1(INDMAX), WORK2(INDMAX)
    COMPLEX(q) :: CR(*)
! local
    COMPLEX(q) :: CTMP
    GDEFN :: CTMPN
    INTEGER IND,IP

    PROFILING_START('crrexp_mul_work_add')

!DIR$ IVDEP
!NEC$ ivdep
!OCL NOVREC
!$OMP SIMD PRIVATE(IP,CTMP,CTMPN)
    DO IND=1,INDMAX
       IP=NLI(IND)
       CTMPN=GCMPLX( WORK1(IND), WORK2(IND), q)
       CTMP =CONJG(CRREXP(IND))
       CR(IP)=CR(IP)+CTMPN*CTMP
    ENDDO

    PROFILING_STOP('crrexp_mul_work_add')

  END SUBROUTINE CRREXP_MUL_WORK_ADD


  SUBROUTINE CRREXP_MUL_WORK_GADD( INDMAX, CRREXP, NLI, WORK1, WORK2, CR)
    USE prec
    IMPLICIT NONE
    INTEGER INDMAX
    COMPLEX(q), INTENT(IN) :: CRREXP(INDMAX)
    INTEGER,    INTENT(IN) :: NLI(INDMAX)
    REAL(qn),   INTENT(IN) :: WORK1(INDMAX), WORK2(INDMAX)
    GDEF       :: CR(*)
! local
    COMPLEX(q) :: CTMP
    GDEFN :: CTMPN
    INTEGER IND,IP
!!DIR$ IVDEP
!NEC$ ivdep
!!OCL NOVREC
!$OMP SIMD PRIVATE(IP,CTMP,CTMPN)
    DO IND=1,INDMAX
       IP=NLI(IND)
       CTMPN=GCMPLX( WORK1(IND), WORK2(IND), q)
       CTMP =CONJG(CRREXP(IND))
       CR(IP)=CR(IP)+CTMPN*CTMP
    ENDDO
  END SUBROUTINE CRREXP_MUL_WORK_GADD


  SUBROUTINE GEMMV_AND_SCATTER(INDMAX,LMMAXC,A,X1,X2,CP,RINPL)
    USE prec
    IMPLICIT NONE
    INTEGER,    INTENT(IN) :: INDMAX,LMMAXC
    REAL(qn),   INTENT(IN) :: A(INDMAX,LMMAXC)
    REAL(qn),   INTENT(IN) :: X1(INDMAX),X2(INDMAX)
    REAL(qn),   INTENT(IN) :: RINPL
    GDEF                   :: CP(*)
    !local variable
    INTEGER I,J
    REAL(qn) R1,R2

#ifdef gammareal
    DO I=1,LMMAXC
      R1=0
      DO J=1,INDMAX
        R1 = R1 + A(J,I)*X1(J)
      ENDDO
      CP(I)=CP(I)+R1*RINPL
    ENDDO
#else
    DO I=1,LMMAXC
      R1=0
      R2=0
      DO J=1,INDMAX
        R1 = R1 + A(J,I)*X1(J)
        R2 = R2 + A(J,I)*X2(J)
      ENDDO
      CP(I)=CP(I)+CMPLX(R1*RINPL,R2*RINPL,KIND=q)
    ENDDO
#endif
  END SUBROUTINE GEMMV_AND_SCATTER

!
! this version assumes that only X1 is setup, X2 is not passed to the routines
!

  SUBROUTINE DGEMMV_AND_SCATTER(INDMAX,LMMAXC,A,X1,CP,RINPL)
    USE prec
    IMPLICIT NONE
    INTEGER,    INTENT(IN) :: INDMAX,LMMAXC
    REAL(qn),   INTENT(IN) :: A(INDMAX,LMMAXC)
    REAL(qn),   INTENT(IN) :: X1(INDMAX)
    REAL(qn),   INTENT(IN) :: RINPL
    GDEF                   :: CP(*)
    !local variable
    INTEGER I,J
    REAL(qn) R1

    DO I=1,LMMAXC
      R1=0
      DO J=1,INDMAX
        R1 = R1 + A(J,I)*X1(J)
      ENDDO
      CP(I)=CP(I)+R1*RINPL
    ENDDO
  END SUBROUTINE DGEMMV_AND_SCATTER
