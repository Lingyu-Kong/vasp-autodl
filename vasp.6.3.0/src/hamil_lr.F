#define debug
#include "symbol.inc"
MODULE hamil_lr
  USE prec
CONTAINS
!************************ SUBROUTINE LR_HAMIL ***************************
!
! this subroutine evaluates the "first order perturbation"
!    |xi>  = H(1) - e(0) S(1) - e(1) S(0) | phi>
! when an ion moves
! H(1), epsilon(1) and S(1) are the first order change of the
! Hamiltonian, eigenenergies and overlap, respectively
! more specifically:
!    |xi>  = V(1) + sum_ij | p_i> D(1)_ij <p_j | phi > 
!                 + sum_ij | d/dR p_i> D(0)_ij - e(0) Q(0)_ij  <p_j | phi >
!                 + sum_ij | p_i> D(0)_ij - e(0) Q(0)_ij < d/dR p_j | phi >
!                 - sum_ij | p_i> e(1) Q(0)_ij <p_j | phi >
!                 - e(1) |phi>
! for details see also HAMILMU_LR
! the matrix D(0) is supplied in the array CDIJ0
! the direction for which the derivative is calculated is supplied 
! by an ion index ION and the direction IDIR (in cartesian coordinates)
! e(1) is returned in WXI%CELEN
! < d/dR p_j | phi > is returned in WXI%CPROJ,
! whereas WXI%CW is set to <G | xi>
!
! furthermore the routine calculates the first order change of
! the onsite occupancy matrix for the considered ion (CRHODE), 
!  sum_nk,n <phi| (| d/dR p_i>   <p_j |
!                 +| p_i>  < d/dR p_j |)| phi>
! and the term
!         <phi|(sum_ij | d/dR p_i> Q(0)_ij  <p_j | 
!              +sum_ij | p_i> Q(0)_ij < d/dR p_j |) |phi>
! (change of norm) for each wavefunction
! this term is returned in WXI%FERWE (somewhat unclean)
!
!***********************************************************************

  SUBROUTINE LR_HAMIL(GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W0, WXI, WDES, &
        LMDIM,CDIJ, CDIJ0, CQIJ, CRHODE, ION, IDIR,  SV, CSHIFT, RMS, ICOUEV, LSUBTRACT)
#ifdef _OPENACC
    USE mopenacc
#endif
    USE prec
    USE wave
    USE wave_high
    USE base
    USE lattice
    USE mpimy
    USE mgrid
    USE nonl_high
    USE hamil
    USE constant
    USE wave_mpi

    IMPLICIT NONE

    TYPE (grid_3d)     GRID
    TYPE (info_struct) INFO
    TYPE (latt)        LATT_CUR
    TYPE (nonlr_struct) NONLR_S
    TYPE (nonl_struct) NONL_S
    TYPE (wavespin)    W0, WXI
    TYPE (wavedes)     WDES

    RGRID   SV(DIMREAL(GRID%MPLWV),WDES%NCDIJ) ! local potential
    INTEGER LMDIM  
    OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    OVERLAP CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    OVERLAP CDIJ0(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)
    OVERLAP CRHODE(LMDIM,LMDIM,WDES%NCDIJ)
    REAL(q) RMS        ! magnitude of the residual vector
    INTEGER ICOUEV     ! number of H | phi> evaluations
    INTEGER ION, IDIR
    LOGICAL LSUBTRACT  ! subtract - sum_ij | p_i> e(1) Q(0)_ij <p_j | phi > - e(1) |phi>
    REAL(q) CSHIFT
! local work arrays and structures
    TYPE (wavedes1)  WDES1           ! descriptor for one k-point
    TYPE (wavefun1)  W0_1(WDES%NSIM) ! current wavefunction
    TYPE (wavefun1)  WTMP(WDES%NSIM) ! see below
    TYPE (nonlr_struct) NONLR_ION,NONLR_IOND

    REAL(q) :: ORTH (WDES%NSIM), ORTH_
    REAL(q) :: FNORM(WDES%NSIM),FNORM_

    COMPLEX(q), ALLOCATABLE :: CF(:,:)
    GDEF, TARGET, ALLOCATABLE :: CPROJ(:,:)

    INTEGER :: NSIM                  ! number of bands treated simultaneously
    INTEGER :: NODE_ME, IONODE
    INTEGER :: NP, ISP, NK, NPL, NGVECTOR, NB_DONE, N, IDUMP, ISPINOR, LD, M, MM
    INTEGER :: NB(WDES%NSIM)         ! contains a list of bands currently optimized
    REAL(q) :: EVALUE(WDES%NSIM)     ! eigenvalue during optimization
    REAL(q) :: DS(WDES%NSIM)         ! change of norm related to update of |p_i>
    LOGICAL :: LDO(WDES%NSIM)        ! band finished
    LOGICAL :: LSTOP

      PROFILING_START('lr_hamil')

#ifdef _OPENACC_TST
      PUSH_ACC_EXEC_ON(.TRUE.)

      ACC_ASYNC_Q=ACC_ASYNC_ASYNC

      CALL ACC_COPYIN_TYPED_VAR(GRID)
      CALL ACC_COPYIN_TYPED_VAR(W0)
      CALL ACC_COPYIN_TYPED_VAR(WXI)
      CALL ACC_COPYIN_TYPED_VAR(WDES)
      CALL ACC_COPYIN_TYPED_VAR(NONL_S)
      CALL ACC_COPYIN_TYPED_VAR(NONLR_S)

!$ACC ENTER DATA COPYIN(CDIJ0,CQIJ,LATT_CUR,CDIJ,SV) __IF_ASYNC__
!$ACC ENTER DATA CREATE(CRHODE) __IF_ASYNC__
#endif

      NODE_ME=0
      IONODE =0
#ifdef MPI
      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE
#endif
!=======================================================================
!  INITIALISATION:
! maximum  number of bands treated simultaneously
!=======================================================================
      RMS   =0
      ICOUEV=0
      NSIM=WDES%NSIM

      ALLOCATE(CF(WDES%NRPLWV,NSIM),CPROJ(WDES%NPROD,NSIM))
!$ACC ENTER DATA CREATE(WDES1,W0_1(:),WTMP(:),CF,CPROJ,DS,ORTH,FNORM) __IF_ASYNC__
      LD=WDES%NRPLWV

      CALL SETWDES(WDES, WDES1, 0)

      DO NP=1,NSIM
         CALL NEWWAV_R(W0_1(NP),WDES1)
      ENDDO

      CALL NONLR_SET_SINGLE_ION(GRID, LATT_CUR, NONLR_S, NONLR_ION, NONLR_IOND, ION, IDIR)

!$ACC KERNELS PRESENT(CRHODE,WXI) __IF_ASYNC__
      CRHODE=0
      WXI%CW=0
      WXI%CPROJ=0
!$ACC END KERNELS
!=======================================================================
      spin:    DO ISP=1,WDES%ISPIN
      kpoints: DO NK=1,WDES%NKPTS
#ifdef MPI
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
#endif
      CALL SETWDES(WDES,WDES1,NK); CALL SETWGRID_OLD(WDES1,GRID) 

      NPL=WDES1%NPL
      NGVECTOR=WDES1%NGVECTOR

      IF (INFO%LREAL) THEN
        CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES)
        ! not very clean but these descriptors use same phase factors as NONLR_S
        NONLR_ION%NK=NK
        NONLR_IOND%NK=NK
!$ACC UPDATE DEVICE(NONLR_ION%NK,NONLR_IOND%NK) __IF_ASYNC__
      ELSE
        CALL PHASE(WDES,NONL_S,NK)
      ENDIF

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)

!=======================================================================
      NB_DONE=0                 ! index of the bands allready optimised
      bands: DO

        NB=0                    ! empty the list of bands, which are optimized currently
!
!  check the NB list, whether there is any empty slot
!  fill in a not yet optimized wavefunction into the slot
!
        IDUMP=0

        PROFILING_START('newband')

        newband: DO NP=1,NSIM
        IF (NB_DONE < WDES%NBANDS ) THEN

           SET_ACC_ASYNC_Q(NP)

           NB_DONE=NB_DONE+1
           N     =NB_DONE
           NB(NP)=NB_DONE
           ICOUEV=ICOUEV+1

           CALL SETWAV(W0,W0_1(NP),WDES1,N,ISP)  ! fill band N into W0_1(NP)
#ifdef MPI
           IF (NODE_ME /= IONODE) IDUMP=0
#endif
           ! fft to real space
           DO ISPINOR=0,WDES%NRSPINORS-1
              CALL FFTWAV(NGVECTOR, WDES%NINDPW(1,NK),W0_1(NP)%CR(1+ISPINOR*GRID%MPLWV),W0_1(NP)%CW(1+ISPINOR*NGVECTOR),GRID)
           ENDDO
           ! WTMP is identical to W0_1, except for the CPROJ component
           ! which will contain the derivative of the projectors with
           ! respect to the displaced ion, after the call HAMILMU_LR
           WTMP(NP)=W0_1(NP)
           WTMP(NP)%CPROJ => CPROJ(:,NP)

#ifdef _OPENACC
!!!$ACC ENTER DATA COPYIN(WTMP(NP)%CW,WTMP(NP)%CPROJ,WTMP(NP)%WDES1) __IF_ASYNC__
           CALL ACC_COPYIN_TYPED_VAR(WTMP(NP))
#endif
           EVALUE(NP)=W0_1(NP)%CELEN
        ENDIF
        ENDDO newband

        PROFILING_STOP('newband')

!=======================================================================
! if the NB list is now empty end the bands DO loop
!=======================================================================
        LSTOP=.TRUE.
        LDO  =.FALSE.
        DO NP=1,NSIM
           IF ( NB(NP) /= 0 ) THEN
              LSTOP  =.FALSE.
              LDO(NP)=.TRUE.     ! band not finished yet
           ENDIF
        ENDDO
        IF (LSTOP) EXIT bands

!=======================================================================
! determine gradient and store it
!=======================================================================
        !  store H | psi > temporarily
        !  to have uniform stride for result array
        CALL HAMILTMU_LR(WDES1,W0_1,WTMP,NONLR_S,NONLR_ION,NONLR_IOND,NONL_S, &
                & GRID,  INFO%LREAL, EVALUE, LATT_CUR, &
                & LMDIM,CDIJ(1,1,1,ISP),CDIJ0(1,1,1,ISP),CQIJ(1,1,1,ISP), CRHODE(1,1,ISP),DS, &
                & SV(1,ISP), CF(1,1),LD, NSIM, LDO, WDES%WTKPT(NK), ION, IDIR, LSUBTRACT, CSHIFT)

        PROFILING_START('i2')

        i2: DO NP=1,NSIM
           N=NB(NP); IF (.NOT. LDO(NP)) CYCLE i2
#ifdef _OPENACC
           SET_ACC_ASYNC_Q(NP)
           CALL KERNEL1(NP,FNORM(NP),ORTH(NP))
#else
           ORTH_ =0
           FNORM_=0

           DO ISPINOR=0,WDES%NRSPINORS-1
              DO M=1,NGVECTOR
                 MM=M+ISPINOR*NGVECTOR
                 IF (LSUBTRACT) THEN
                    CF(MM,NP)=CF(MM,NP)-W0_1(NP)%CELEN*W0_1(NP)%CW(MM)
                 ENDIF
                 IF (WDES%LSPIRAL.AND.WDES%DATAKE(M,ISPINOR+1,NK)>INFO%ENINI) CF(MM,NP)=0
                 FNORM_ =FNORM_+CF(MM,NP)*CONJG(CF(MM,NP))
                 ORTH_  =ORTH_ +CF(MM,NP)*CONJG(W0_1(NP)%CW(MM))
                 WXI%CW(MM,N,NK,ISP)=CF(MM,NP)
              ENDDO
           ENDDO

           CALLMPI( M_sum_2(WDES%COMM_INB, FNORM_, ORTH_))

           FNORM(NP)=FNORM_ ; ORTH(NP)=ORTH_
#endif
           CALL __GCOPY__(WDES%NPROD,WTMP(NP)%CPROJ(1),1,WXI%CPROJ(1,N,NK,ISP),1)
!!!$ACC KERNELS PRESENT(WXI,W0_1(NP),DS(NP)) __IF_ASYNC__
!!           WXI%CPROJ(:,N,NK,ISP)=WTMP(NP)%CPROJ
!$ACC KERNELS PRESENT(WXI,W0_1,DS) __IF_ASYNC__
           WXI%CELEN(N,NK,ISP)=W0_1(NP)%CELEN
           WXI%FERWE(N,NK,ISP)=DS(NP)
!$ACC END KERNELS

!!!$ACC EXIT DATA DELETE(WTMP(NP)%CW,WTMP(NP)%CPROJ,WTMP(NP)%WDES1) __IF_ASYNC__
#ifdef _OPENACC
           CALL ACC_DELETE_TYPED_VAR(WTMP(NP))
#endif
!=======================================================================
! move onto the next block of bands
!=======================================================================
        ENDDO i2

        DO NP=1,NSIM
           N=NB(NP); IF (.NOT. LDO(NP)) CYCLE

           SET_ACC_ASYNC_Q(NP)
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)

           IF (ABS(ORTH(NP))>1E-4.AND.LSUBTRACT) THEN
              WRITE(0,*)'HAMIL_LR internal error: the vector H(1)-e(1) S(1) |phi(0)> is not orthogonal to |phi(0)>',NK,N,ORTH(NP)
!              STOP
           ENDIF

           IF (IDUMP==2) THEN
#if PGI_BEFORE_XX_X
!$ACC UPDATE SELF(W0_1(NP)%CELEN,WXI%CELEN(N:N,NK:NK,ISP:ISP),WXI%FERWE(N:N,NK:NK,ISP:ISP)) IF(ACC_EXEC_ON)
#else
!$ACC UPDATE SELF(W0_1(NP)%CELEN,WXI%CELEN(N,NK,ISP),WXI%FERWE(N,NK,ISP)) IF(ACC_EXEC_ON)
#endif
               WRITE(*,'(I3,E11.4,"R ",E11.4,"E ",E11.4,"O ",2E14.7)') N,SQRT(ABS(FNORM(NP))),REAL(W0_1(NP)%CELEN,q),ORTH(NP),REAL(WXI%CELEN(N,NK,ISP)),WXI%FERWE(N,NK,ISP)
           ENDIF

           RMS=RMS+WDES%RSPIN*WDES%WTKPT(NK)*W0%FERWE(N,NK,ISP)*SQRT(ABS(FNORM(NP)))/WDES%NB_TOT
        ENDDO

        PROFILING_STOP('i2')

!=======================================================================
      ENDDO bands
      ENDDO kpoints
      ENDDO spin
!=======================================================================

      CALLMPI( M_sum_d(WDES%COMM_INTER, RMS, 1))
      CALLMPI( M_sum_d(WDES%COMM_KINTER, RMS, 1))
#ifdef realmode
      CALLMPI( M_sum_d(WDES%COMM_INTER, CRHODE, LMDIM*LMDIM*WDES%NCDIJ))
      CALLMPI( M_sum_d(WDES%COMM_KINTER, CRHODE, LMDIM*LMDIM*WDES%NCDIJ))
#else
      CALLMPI( M_sum_d(WDES%COMM_INTER, CRHODE, LMDIM*LMDIM*WDES%NCDIJ*2))
      CALLMPI( M_sum_d(WDES%COMM_KINTER, CRHODE, LMDIM*LMDIM*WDES%NCDIJ*2))
#endif

!$ACC WAIT IF(ACC_EXEC_ON)
DOACC ACC_ASYNC_Q=ACC_ASYNC_ASYNC

      DO NP=1,NSIM
         CALL DELWAV_R(W0_1(NP))
      ENDDO
!$ACC EXIT DATA DELETE(W0_1(:),WTMP(:),CF,CPROJ,DS,ORTH,FNORM) __IF_ASYNC__
      DEALLOCATE(CF, CPROJ)

      CALL  NONLR_DEALLOC_SINGLE_ION(NONLR_ION)
      CALL  NONLR_DEALLOC_SINGLE_ION(NONLR_IOND)

#ifdef _OPENACC
      CALL ACC_DELETE_TYPED_VAR(WDES1)
#endif
#ifdef _OPENACC_TST
!$ACC EXIT DATA DELETE(CDIJ0,CQIJ,LATT_CUR,CDIJ,SV) __IF_ASYNC__
!$ACC EXIT DATA COPYOUT(CRHODE) __IF_ASYNC__

      CALL ACC_DELETE_TYPED_VAR(GRID)
      CALL ACC_DELETE_TYPED_VAR(WDES)
      CALL ACC_COPYOUT_TYPED_VAR(WXI)
      CALL ACC_DELETE_TYPED_VAR(W0)
      CALL ACC_DELETE_TYPED_VAR(NONL_S)
      CALL ACC_DELETE_TYPED_VAR(NONLR_S)

!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
      ACC_ASYNC_Q=ACC_ASYNC_SYNC

      POP_ACC_EXEC_ON
#endif

      PROFILING_STOP('lr_hamil')
#ifdef _OPENACC
      CONTAINS

         SUBROUTINE KERNEL1(NP,FNORM,ORTH)
         REAL(q) :: FNORM,ORTH
         INTEGER :: NP
         ! local
         INTEGER :: ISPINOR,M,MM
!$ACC KERNELS PRESENT(ORTH,FNORM) __IF_ASYNC__
         ORTH =0
         FNORM=0
!$ACC END KERNELS

         IF (LSUBTRACT) THEN
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(WDES,CF,W0_1,WDES,WXI,ORTH,FNORM) PRIVATE(MM) REDUCTION(+:FNORM,ORTH) __IF_ASYNC__
            DO ISPINOR=0,WDES%NRSPINORS-1
               DO M=1,NGVECTOR
                  MM=M+ISPINOR*NGVECTOR
                  CF(MM,NP)=CF(MM,NP)-W0_1(NP)%CELEN*W0_1(NP)%CW(MM)
                  IF (WDES%LSPIRAL.AND.WDES%DATAKE(M,ISPINOR+1,NK)>INFO%ENINI) CF(MM,NP)=0
                  FNORM =FNORM+CF(MM,NP)*CONJG(CF(MM,NP))
                  ORTH  =ORTH +CF(MM,NP)*CONJG(W0_1(NP)%CW(MM))
                  WXI%CW(MM,N,NK,ISP)=CF(MM,NP)
               ENDDO
            ENDDO
         ELSE
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(WDES,CF,W0_1,WDES,WXI,ORTH,FNORM) PRIVATE(MM) REDUCTION(+:FNORM,ORTH) __IF_ASYNC__
            DO ISPINOR=0,WDES%NRSPINORS-1
               DO M=1,NGVECTOR
                  MM=M+ISPINOR*NGVECTOR
                  IF (WDES%LSPIRAL.AND.WDES%DATAKE(M,ISPINOR+1,NK)>INFO%ENINI) CF(MM,NP)=0
                  FNORM =FNORM+CF(MM,NP)*CONJG(CF(MM,NP))
                  ORTH  =ORTH +CF(MM,NP)*CONJG(W0_1(NP)%CW(MM))
                  WXI%CW(MM,N,NK,ISP)=CF(MM,NP)
               ENDDO
            ENDDO
         ENDIF

         CALLMPI( M_sum_d(WDES%COMM_INB, FNORM, 1))
         CALLMPI( M_sum_d(WDES%COMM_INB, ORTH,  1))
!$ACC UPDATE SELF(ORTH,FNORM) __IF_ASYNC__
         RETURN
         END SUBROUTINE KERNEL1
#endif
    END SUBROUTINE LR_HAMIL


!************************* SUBROUTINE HAMILTMU_LR *********************
!
! this subroutine calculates the first order change of H
! acting onto a set of wavefuntions
! the  wavefunction must be given in reciprocal space C and real
! space CR
! CH contains the result
!      |xi>  = V(1) | phi > + | p_i> d/dR D_ij <p_j | phi > 
!             + sum_ij |d p_i/d R > D_ij - e(0) Q_ij <p_j | phi >
!             + sum_ij | p_i> D_ij - e(0) Q_ij <d p_j/d R | phi >
!                 - sum_ij | p_i> e(1) Q_ij <p_j | phi >
!
! V(1) is the first order change of the local potential  (SV)
! D(1) is the first order change of the PAW strenght     (CDIJ)
! D(0) is the original strength                          (CDIJ0)
! e(0) is the zero order eigen energy
! e(1) is the first order change of the eigen energy
!
! e(1) is evaluated during the calculation of |xi>:
! e(1) = <phi| V(1) | phi > + <phi | p_i> d/dR D_ij <p_j | phi >
!        +  sum_ij <phi |d p_i/d R > D_ij - e Q_ij <p_j | phi > + c.c.
!
! the direction for which the derivative is calculated is supplied 
! by an ion index ION and the direction IDIR (in cartesian coordinates)
!
! furthermore the routine calculates the first order change of
! the onsite occupancy matrix for the considered ion
!  sum_nk,n <phi| (| d/dR p_i>   <p_j |
!                 +| p_i>  < d/dR p_j |)| phi>
! this is not particularly elegant but the only place
! where this quantity can be calculated without further
! complications
!
! NOTE: the calling routine has to subtract e(1) |phi> to get the
! correct vector xi
!
!***********************************************************************

    SUBROUTINE HAMILTMU_LR( &
         WDES1,W0_1,WTMP,NONLR_S,NONLR_ION,NONLR_IOND,NONL_S, &
         GRID, LREAL, EVALUE0, LATT_CUR, &
         LMDIM,CDIJ,CDIJ0,CQIJ,CRHODE,DS, &
         SV,CH,LD, NSIM, LDO, WTKPT, ION, IDIR, LSUBTRACT, CSHIFT)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE mpimy
      USE mgrid
      USE wave
      USE nonl_high
      USE lattice
      USE hamil

      IMPLICIT NONE

      INTEGER NSIM,NP,LD
      INTEGER LMDIM, NGVECTOR, ISPINOR, ISPINOR_, MM, MM_ 
      TYPE (grid_3d)     GRID
      TYPE (nonlr_struct) NONLR_S,NONLR_ION, NONLR_IOND
      TYPE (nonl_struct) NONL_S
      TYPE (wavefun1)    W0_1(NSIM)
      TYPE (wavefun1)    WTMP(NSIM)
      TYPE (wavedes1)    WDES1
      TYPE (latt)        LATT_CUR

      RGRID      SV(DIMREAL(GRID%MPLWV),WDES1%NRSPINORS*WDES1%NRSPINORS) ! local potential
      OVERLAP    CDIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
                 CQIJ(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS), &
                 CDIJ0(LMDIM,LMDIM,WDES1%NIONS,WDES1%NRSPINORS*WDES1%NRSPINORS)
      OVERLAP    CRHODE(LMDIM,LMDIM,WDES1%NRSPINORS*WDES1%NRSPINORS)
      COMPLEX(q) CH(LD,NSIM)
      REAL(q)    EVALUE0(NSIM)
      LOGICAL LREAL
      LOGICAL LDO(NSIM)
      INTEGER ION, IDIR
      REAL(q) WTKPT
      REAL(q) CSHIFT
      LOGICAL LSUBTRACT
! local variables
      COMPLEX(q), ALLOCATABLE :: CWORK1(:,:)
      COMPLEX(q) :: EVALUE0_(NSIM), CTMP(NSIM), CE
      REAL(q) :: DS(NSIM), EVALUE1(NSIM), RINPLW, WEIGHT
      INTEGER :: LMBASE, LMBASE_, NIS, LMMAXC, NI, L, LP, NT, M

      PROFILING_START('hamiltmu_lr')

      ALLOCATE(CWORK1(GRID%MPLWV*WDES1%NRSPINORS,NSIM)) 
!$ACC ENTER DATA CREATE(CWORK1,CTMP,EVALUE1) IF(ACC_EXEC_ON)

      RINPLW=1._q/GRID%NPLWV
      NGVECTOR=WDES1%NGVECTOR

      EVALUE0_=EVALUE0+CMPLX(0.0_q,CSHIFT,q)
!=======================================================================
! calculate the local contribution (result in CWORK1)
!=======================================================================
      DO NP=1,NSIM
         IF ( LDO(NP) ) THEN
#ifdef _OPENACC
            SET_ACC_ASYNC_Q(NP)
            CALL KERNEL1(NP,W0_1(NP)%CELEN)
#else
            CE=0
            CWORK1(:,NP)=0

            DO ISPINOR =0,WDES1%NRSPINORS-1
            DO ISPINOR_=0,WDES1%NRSPINORS-1
               DO M=1,GRID%RL%NP
                  MM =M+ISPINOR *GRID%MPLWV
                  MM_=M+ISPINOR_*GRID%MPLWV
                  CWORK1(MM,NP)=  CWORK1(MM,NP)+(SV(M,1+ISPINOR_+2*ISPINOR) *W0_1(NP)%CR(MM_)*RINPLW)
                  CE=CE + CONJG(W0_1(NP)%CR(MM)) *(SV(M,1+ISPINOR_+2*ISPINOR) *W0_1(NP)%CR(MM_)*RINPLW)
               ENDDO
            ENDDO
            ENDDO

            CALLMPI( M_sum_z(WDES1%COMM_INB, CE, 1))

            W0_1(NP)%CELEN=CE
#endif
         ENDIF
      ENDDO
!=======================================================================
! non-local contribution in real-space
!=======================================================================
      IF (LREAL) THEN
         ! contribution | d p_i/ d R > D_ij - e(0) Q_ij < p_j |
         IF (CSHIFT==0) THEN
            CALL RACCMU(NONLR_IOND,WDES1,W0_1,LMDIM,CDIJ0,CQIJ,EVALUE0,CWORK1,GRID%MPLWV*WDES1%NRSPINORS,NSIM,LDO)
         ELSE
            CALL RACCMU_C(NONLR_IOND,WDES1,W0_1,LMDIM,CDIJ0,CQIJ,EVALUE0_,CWORK1,GRID%MPLWV*WDES1%NRSPINORS,NSIM,LDO)
         ENDIF
         ! contribution | p_i > D_ij - e(0) Q_ij < d p_j/ d R |
#ifdef _OPENACC
         IF(ACC_EXEC_ON) THEN
            CALL RPROMU_ACC(NONLR_IOND,WDES1,WTMP, NSIM, LDO)
         ELSE
#endif
            CALL RPROMU(NONLR_IOND,WDES1,WTMP, NSIM, LDO)
#ifdef _OPENACC
         ENDIF
#endif
         IF (CSHIFT==0) THEN
            CALL RACCMU(NONLR_ION,WDES1,WTMP,LMDIM,CDIJ0,CQIJ,EVALUE0,CWORK1,GRID%MPLWV*WDES1%NRSPINORS,NSIM,LDO)
         ELSE
! note(mM): RACCMU_C will start putting work into NSIM asynchronous execution queues. The work started
!   in RPROMU_ACC above was put into the ACC_ASYNC_ASYNC queue, so we have to wait for this to finish.
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_EXEC_ON)
            CALL RACCMU_C(NONLR_ION,WDES1,WTMP,LMDIM,CDIJ0,CQIJ,EVALUE0_,CWORK1,GRID%MPLWV*WDES1%NRSPINORS,NSIM,LDO)
         ENDIF

         ! non local contributions to e(1)
         DO NP=1,NSIM
            IF ( LDO(NP) ) THEN

               SET_ACC_ASYNC_Q(NP)

               CALL ECCP_NL_ALL(WDES1,W0_1(NP),WTMP(NP),CDIJ0,CQIJ,EVALUE0(NP),CTMP(NP))
!$ACC KERNELS PRESENT(W0_1,CTMP) __IF_ASYNC__
               W0_1(NP)%CELEN=W0_1(NP)%CELEN+2*REAL(CTMP(NP),q)
!$ACC END KERNELS
               CALL ECCP_NL_ALL(WDES1,W0_1(NP),W0_1(NP),CDIJ,CQIJ,0.0_q,CTMP(NP))
!$ACC KERNELS PRESENT(W0_1,CTMP) __IF_ASYNC__
               W0_1(NP)%CELEN=W0_1(NP)%CELEN+CTMP(NP)
               W0_1(NP)%CELEN=REAL(W0_1(NP)%CELEN,q)
!$ACC END KERNELS

               IF (LSUBTRACT) THEN
!$ACC KERNELS PRESENT(EVALUE1,W0_1) __IF_ASYNC__
                  EVALUE1(NP)=W0_1(NP)%CELEN
!$ACC END KERNELS
               ELSE
!$ACC KERNELS PRESENT(EVALUE1,W0_1) __IF_ASYNC__
                  EVALUE1(NP)=0
!$ACC END KERNELS
               ENDIF
            ENDIF
         ENDDO

         ! contribution | p_i > d D_ij / d R - e(1) Q_ij < p_j |
         CALL RACCMU(NONLR_S,WDES1,W0_1,LMDIM,CDIJ,CQIJ,EVALUE1,CWORK1,GRID%MPLWV*WDES1%NRSPINORS,NSIM,LDO)

         DO NP=1,NSIM
            IF ( LDO(NP) ) THEN

               SET_ACC_ASYNC_Q(NP)

               DO ISPINOR=0,WDES1%NRSPINORS-1
                  CALL FFTEXT(NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV,NP),CH(1+ISPINOR*NGVECTOR,NP),GRID,.FALSE.)
               ENDDO
            ENDIF
         ENDDO
!=======================================================================
! calculate the non local contribution in reciprocal space
!=======================================================================
      ELSE
         DO NP=1,NSIM
            IF ( LDO(NP) ) THEN

               SET_ACC_ASYNC_Q(NP)

!$ACC KERNELS PRESENT(CH) __IF_ASYNC__
               CH(:,NP)=0
!$ACC END KERNELS
!!! test_
!!!$ACC ENTER DATA COPYIN(EVALUE0_(NP)) __IF_ASYNC__
!!! test_
               ! contribution | d p_i/ d R > D_ij - epsilon Q_ij < p_j |
               IF (CSHIFT==0) THEN
! test_
#if PGI_BEFORE_XX_X
!$ACC ENTER DATA COPYIN(EVALUE0(NP:NP) ) __IF_ASYNC__
#else
!$ACC ENTER DATA COPYIN(EVALUE0(NP) ) __IF_ASYNC__
#endif
! test_
                  CALL VNLACC_DER(NONL_S,W0_1(NP),CDIJ0,CQIJ,EVALUE0(NP),CH(:,NP), &
                       LATT_CUR, ION, IDIR )
               ELSE
! test_
#if PGI_BEFORE_XX_X
!$ACC ENTER DATA COPYIN(EVALUE0_(NP:NP)) __IF_ASYNC__
#else
!$ACC ENTER DATA COPYIN(EVALUE0_(NP)) __IF_ASYNC__
#endif
! test_
                  CALL VNLACC_DER_C(NONL_S,W0_1(NP),CDIJ0,CQIJ,EVALUE0_(NP),CH(:,NP), &
                       LATT_CUR, ION, IDIR )
               ENDIF
               ! contribution | p_i > D_ij - epsilon Q_ij < d p_j/ d R |
               CALL PROJ_DER(NONL_S, WDES1, WTMP(NP), LATT_CUR, ION, IDIR)
               IF (CSHIFT==0) THEN
                  CALL VNLACC_DER(NONL_S,WTMP(NP),CDIJ0,CQIJ,EVALUE0(NP),CH(:,NP), &
                       LATT_CUR, ION, 0 )
! test_
#if PGI_BEFORE_XX_X
!$ACC EXIT DATA DELETE(EVALUE0(NP:NP) ) __IF_ASYNC__
#else
!$ACC EXIT DATA DELETE(EVALUE0(NP) ) __IF_ASYNC__
#endif
! test_
               ELSE
                  CALL VNLACC_DER_C(NONL_S,WTMP(NP),CDIJ0,CQIJ,EVALUE0_(NP),CH(:,NP), &
                       LATT_CUR, ION, 0 )
! test_
#if PGI_BEFORE_XX_X
!$ACC EXIT DATA DELETE(EVALUE0_(NP:NP)) __IF_ASYNC__
#else
!$ACC EXIT DATA DELETE(EVALUE0_(NP)) __IF_ASYNC__
#endif
! test_
               ENDIF
!!! test_
!!!$ACC EXIT DATA DELETE(EVALUE0_(NP)) __IF_ASYNC__
!!! test_
               ! non local contributions to e(1)
               CALL ECCP_NL_ALL(WDES1,W0_1(NP),WTMP(NP),CDIJ0,CQIJ,EVALUE0(NP),CTMP(NP))
!$ACC KERNELS PRESENT(W0_1,CTMP) __IF_ASYNC__
               W0_1(NP)%CELEN=W0_1(NP)%CELEN+2*REAL(CTMP(NP),q)
!$ACC END KERNELS
               CALL ECCP_NL_ALL(WDES1,W0_1(NP),W0_1(NP),CDIJ,CQIJ,0.0_q,CTMP(NP))
!$ACC KERNELS PRESENT(W0_1,CTMP) __IF_ASYNC__
               W0_1(NP)%CELEN=W0_1(NP)%CELEN+CTMP(NP)
               W0_1(NP)%CELEN=REAL(W0_1(NP)%CELEN,q)
!$ACC END KERNELS
               IF (LSUBTRACT) THEN
!$ACC KERNELS PRESENT(EVALUE1,W0_1) __IF_ASYNC__
                  EVALUE1(NP)=W0_1(NP)%CELEN
!$ACC END KERNELS
               ELSE
!$ACC KERNELS PRESENT(EVALUE1,W0_1) __IF_ASYNC__
                  EVALUE1(NP)=0
!$ACC END KERNELS
               ENDIF
               ! contribution | p_i > d D_ij / d R - e(1) Q_ij < p_j |
               CALL VNLACC_ADD(NONL_S,W0_1(NP),CDIJ,CQIJ,1,EVALUE1(NP),CH(:,NP))

               DO ISPINOR=0,WDES1%NRSPINORS-1
                  CALL FFTEXT(NGVECTOR,WDES1%NINDPW(1),CWORK1(1+ISPINOR*WDES1%GRID%MPLWV,NP),CH(1+ISPINOR*NGVECTOR,NP),GRID,.TRUE.)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!=======================================================================
! calculate the first order change of the onsite occupancy matrix
!=======================================================================
      DO NP=1,NSIM
         IF ( LDO(NP) ) THEN

            SET_ACC_ASYNC_Q(NP)

            CALL ECCP_NL_ALL(WDES1,W0_1(NP),WTMP(NP),CQIJ,CQIJ,0.0_q,CTMP(NP))
!$ACC KERNELS PRESENT(DS,CTMP) __IF_ASYNC__
            DS(NP)=REAL(CTMP(NP),q)*2
!$ACC END KERNELS
         ENDIF
      ENDDO

#ifdef _OPENACC
      CALL ACC_SYNC_ASYNC_Q(NSIM,LDO)
      ACC_ASYNC_Q=ACC_ASYNC_ASYNC
#endif

      DO NP=1,NSIM
         IF ( LDO(NP) ) THEN
            spinor: DO ISPINOR =0,WDES1%NRSPINORS-1
            DO ISPINOR_=0,WDES1%NRSPINORS-1

               LMBASE =ISPINOR *(WDES1%NPRO/2)
               LMBASE_=ISPINOR_*(WDES1%NPRO/2)

               NIS   =1
               typ:  DO NT=1,WDES1%NTYP
                  LMMAXC=WDES1%LMMAX(NT)
                  IF (LMMAXC==0) GOTO 210

                  DO NI=NIS,WDES1%NITYP(NT)+NIS-1
                  ! is this ion the one we seek
                     IF (NI_GLOBAL(NI, WDES1%COMM_INB)==ION) THEN

                        WEIGHT=WDES1%RSPIN*W0_1(NP)%FERWE*WTKPT
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CRHODE,W0_1,WTMP) __IF_ASYNC__
                        DO L=1,LMMAXC
!DIR$ IVDEP
!OCL NOVREC
                        DO LP=1,LMMAXC
                           CRHODE(LP,L,ISPINOR_+2*ISPINOR+1)=CRHODE(LP,L,ISPINOR_+2*ISPINOR+1)+ &
                                WEIGHT*W0_1(NP)%CPROJ(L+LMBASE)*GCONJG(WTMP(NP)%CPROJ(LP+LMBASE_))+ &
                                WEIGHT*WTMP(NP)%CPROJ(L+LMBASE)*GCONJG(W0_1(NP)%CPROJ(LP+LMBASE_))
                        ENDDO
                        ENDDO
                     ENDIF

                     LMBASE = LMMAXC+LMBASE
                     LMBASE_= LMMAXC+LMBASE_
                  ENDDO

  210             NIS = NIS+WDES1%NITYP(NT)
               ENDDO typ
            ENDDO
            ENDDO spinor

         ENDIF
      ENDDO

!$ACC EXIT DATA DELETE(EVALUE1,CTMP,CWORK1) __IF_ASYNC__
      DEALLOCATE(CWORK1)

      PROFILING_STOP('hamiltmu_lr')
#ifdef _OPENACC
      CONTAINS

         SUBROUTINE KERNEL1(NP,CE)
         COMPLEX(q) :: CE
         INTEGER    :: NP
         ! local
         INTEGER :: ISPINOR,ISPINOR_,M,MM,MM_
!$ACC KERNELS PRESENT(CE,CWORK1) __IF_ASYNC__
         CE=0
         CWORK1(:,NP)=0
!$ACC END KERNELS
         DO ISPINOR =0,WDES1%NRSPINORS-1
!$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(CE,CWORK1,SV,WDES1,W0_1,GRID) PRIVATE(MM,MM_) REDUCTION(+:CE) __IF_ASYNC__
           DO ISPINOR_=0,WDES1%NRSPINORS-1
              DO M=1,GRID%RL%NP
                 MM =M+ISPINOR *GRID%MPLWV
                 MM_=M+ISPINOR_*GRID%MPLWV
                 CWORK1(MM,NP)=  CWORK1(MM,NP)+(SV(M,1+ISPINOR_+2*ISPINOR) *W0_1(NP)%CR(MM_)*RINPLW)
                 CE=CE + CONJG(W0_1(NP)%CR(MM)) *(SV(M,1+ISPINOR_+2*ISPINOR) *W0_1(NP)%CR(MM_)*RINPLW)
              ENDDO
           ENDDO
        ENDDO
        CALLMPI( M_sum_z(WDES1%COMM_INB, CE, 1))
        RETURN
        END SUBROUTINE KERNEL1
#endif
    END SUBROUTINE HAMILTMU_LR

  END MODULE hamil_lr
