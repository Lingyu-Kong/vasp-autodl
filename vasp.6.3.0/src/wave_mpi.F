#include "symbol.inc"
!************************************************************************
! RCS:  $Id: wave_mpi.F,v 1.3 2001/04/03 10:43:18 kresse Exp $
!
!  this module contains the routines required to communicate
!  wavefunctions and projected wavefunctions
!
!  there are also two quite tricky routines SET_WPOINTER
!  which allow to generate pointer to F77-sequenced arrays
!  There is no guarantee that this will work on all computers, but
!  currently it seems to be ok
!
!***********************************************************************

  MODULE wave_mpi
#ifdef _OPENACC
    USE mopenacc
#endif
    USE prec
    USE mpimy
!
! I need an interface block here because I pass the
! first element of a pointer to CW
! most F90 comiler do not allow such a construct
!
      INTERFACE
      SUBROUTINE SET_WPOINTER(CW_P, N1, N2, CW)
      USE prec
      INTEGER N1,N2
      COMPLEX(q), POINTER :: CW_P(:,:)
      COMPLEX(q), TARGET  :: CW
      END SUBROUTINE
      END INTERFACE

      INTERFACE
      SUBROUTINE SET_GPOINTER(CW_P, N1, N2, CW)
      USE prec
      INTEGER N1,N2
      GDEF, POINTER :: CW_P(:,:)
      GDEF, TARGET  :: CW
      END SUBROUTINE
      END INTERFACE


      ! each communcation package needs a unique identifier (ICOMM)
      ! this is handled by the global variable ICOMM_BASE_HANDLE
      ! which is incremented / decremented by ICOMM_INCREMENT
      ! whenever a redistribution handle is allocated

      INTEGER, PARAMETER :: ICOMM_BASE=10000, ICOMM_INCREMENT=1000
      INTEGER, SAVE :: ICOMM_BASE_HANDLE=ICOMM_BASE

      ! is assyncronous communication allowed or not

      LOGICAL :: LASYNC=.FALSE.

      TYPE REDIS_PW_CTR
         INTEGER :: NB                   ! number of bands to be redistributed
         INTEGER :: NBANDS               ! maximum band index
         INTEGER :: NEXT                 ! next vacant slot
         COMPLEX(q), POINTER :: CW(:,:)  ! storage for redistribution
         INTEGER,POINTER :: BAND(:)      ! bands that are currently redistributed
         INTEGER,POINTER :: SREQUEST(:,:)! handles for outstanding send requests
         INTEGER,POINTER :: RREQUEST(:,:)! handles for outstanding receive requests
         INTEGER :: ICOMM                ! identifier for communication
         TYPE(communic),POINTER :: COMM  ! communication handle
      END TYPE REDIS_PW_CTR

    CONTAINS
!************************ SUBROUTINE REDIS_PW_ALLOC ********************
! 
! This subroutine allocates the storage required for asynchronous
! redistribution of wave function coeffients
! the redistribution can be initiated for up to NSIM bands 
! at the same time
!
!***********************************************************************

    SUBROUTINE REDIS_PW_ALLOC(WDES, NSIM, H)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)   :: WDES
      INTEGER :: NSIM                  ! number of bands done simultaneausly
      TYPE (redis_pw_ctr),POINTER :: H ! handle
#ifdef MPI
      INTEGER NRPLWV

      NRPLWV=WDES%NRPLWV

      IF (MOD(NRPLWV,WDES%NB_PAR) /= 0) THEN
        CALL vtutor%bug("REDIS_PW_ALLOC: internal error(1) " // str(NRPLWV) // " " // str(WDES%NB_PAR), __FILE__, __LINE__)
      ENDIF

      ALLOCATE(H)
!$ACC ENTER DATA CREATE(H) __IF_ASYNC__
      H%NB=NSIM
      H%NBANDS=WDES%NBANDS
      H%COMM=>WDES%COMM_INTER
      H%ICOMM=ICOMM_BASE_HANDLE

      ALLOCATE(H%SREQUEST(H%COMM%NCPU,NSIM), H%RREQUEST(H%COMM%NCPU,NSIM))
      ALLOCATE(H%CW(NRPLWV,NSIM),H%BAND(NSIM))
!$ACC ENTER DATA CREATE(H%CW) __IF_ASYNC__
      H%BAND=0

      H%NEXT=1

      ICOMM_BASE_HANDLE=ICOMM_BASE_HANDLE+ICOMM_INCREMENT
#endif
    END SUBROUTINE REDIS_PW_ALLOC


    SUBROUTINE REDIS_PW_DEALLOC(H)
      USE prec
      USE wave
      IMPLICIT NONE

      INTEGER :: NSIM                  ! number of bands done simultaneausly
      TYPE (redis_pw_ctr),POINTER :: H ! handle
      INTEGER :: IND
#ifdef MPI
      IF (.NOT.ASSOCIATED(H)) THEN
        CALL vtutor%bug("REDIS_PW_DEALLOC: internal error(1)", __FILE__, __LINE__)
      ENDIF

      DO IND=1, H%NB
         IF ( H%BAND(IND) /= 0) THEN
            CALL vtutor%bug("internal error: REDIS_PW_DEALLOC not all lines STOPED", __FILE__, __LINE__)
         ENDIF
      ENDDO
!$ACC EXIT DATA DELETE(H%CW) __IF_ASYNC__
      DEALLOCATE(H%CW,H%BAND)
      DEALLOCATE(H%SREQUEST, H%RREQUEST)
!$ACC EXIT DATA DELETE(H) __IF_ASYNC__
      DEALLOCATE(H)

      ICOMM_BASE_HANDLE=ICOMM_BASE_HANDLE-ICOMM_INCREMENT
      IF (ICOMM_BASE_HANDLE < ICOMM_BASE) THEN
        CALL vtutor%bug("REDIS_PW_DEALLOC: internal error(2) " // str(ICOMM_BASE_HANDLE) // " " // &
           str(ICOMM_BASE), __FILE__, __LINE__)
      ENDIF
#endif
    END SUBROUTINE REDIS_PW_DEALLOC

  END MODULE wave_mpi

!=======================================================================
!
! this routine returns a pointer to an SEQUENCED F77 like array
! with a given storage convention
! 1. dimension is N1 second one N2
!
!=======================================================================

      SUBROUTINE SET_WPOINTER(CW_P, N1, N2, CW)
      USE prec
      IMPLICIT NONE
      INTEGER N1,N2
      COMPLEX(q), POINTER :: CW_P(:,:)
      COMPLEX(q), TARGET  :: CW(N1,N2)
      CW_P => CW(:,:)

      END SUBROUTINE

      SUBROUTINE SET_GPOINTER(CW_P, N1, N2, CW)
      USE prec
      IMPLICIT NONE
      INTEGER N1,N2
      GDEF, POINTER :: CW_P(:,:)
      GDEF, TARGET  :: CW(N1,N2)
      CW_P => CW(:,:)

      END SUBROUTINE

!=======================================================================
!  distribute projector part of wavefunctions over all nodes
!=======================================================================

      SUBROUTINE DIS_PROJ(WDES1,CPROJ,CPROJ_LOCAL)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)    WDES1
      GDEF CPROJ(WDES1%NPRO_TOT)
      GDEF CPROJ_LOCAL(WDES1%NPRO)
#ifdef MPI
      INTEGER NPRO,NT,LMMAXC,NPRO_POS,L,NI,NPRO_TOT
!
! quick copy if possible
      IF (WDES1%COMM_INB%NCPU==1) THEN
!$ACC KERNELS PRESENT(CPROJ_LOCAL,WDES1,CPROJ) __IF_ASYNC__
         CPROJ_LOCAL(1:WDES1%NPRO)=CPROJ(1:WDES1%NPRO)
!$ACC END KERNELS
         RETURN
      ENDIF
      CALLMPI( M_sum_g(WDES1%COMM_INB, CPROJ, WDES1%NPRO_TOT))

!$ACC PARALLEL LOOP PRESENT(WDES1,CPROJ,CPROJ_LOCAL) PRIVATE(NT,LMMAXC,NPRO_POS,NPRO,L) __IF_ASYNC__
      DO NI=1,WDES1%NIONS
         NT=WDES1%ITYP(NI)
         LMMAXC=WDES1%LMMAX(NT)
         IF (LMMAXC/=0) THEN
            NPRO_POS=WDES1%NPRO_POS(NI)
            NPRO=WDES1%LMBASE(NI)
!$ACC LOOP
            DO L=1,LMMAXC
               CPROJ_LOCAL(L+NPRO)=CPROJ(L+NPRO_POS)
            ENDDO
          ENDIF
      ENDDO

      IF (WDES1%LNONCOLLINEAR) THEN
         
         NPRO_TOT= WDES1%NPRO_TOT/2
!$ACC PARALLEL LOOP PRESENT(WDES1,CPROJ,CPROJ_LOCAL) PRIVATE(NT,LMMAXC,NPRO_POS,NPRO,L) __IF_ASYNC__
         DO NI=1,WDES1%NIONS
            NT=WDES1%ITYP(NI)
            LMMAXC=WDES1%LMMAX(NT)
            IF (LMMAXC/=0) THEN
               NPRO_POS=WDES1%NPRO_POS(NI)
               NPRO=WDES1%LMBASE(NI)+WDES1%LMBASE(WDES1%NIONS+1)
!$ACC LOOP
               DO L=1,LMMAXC
                  CPROJ_LOCAL(L+NPRO)=CPROJ(L+NPRO_POS+NPRO_TOT)
               ENDDO
             ENDIF
         ENDDO

      ENDIF
#else
!$ACC KERNELS PRESENT(CPROJ_LOCAL,WDES1,CPROJ) __IF_ASYNC__
      CPROJ_LOCAL(1:WDES1%NPRO)=CPROJ(1:WDES1%NPRO)
!$ACC END KERNELS
#endif
      END SUBROUTINE

!=======================================================================
!  merge projector part wavefunctions from all nodes
!  definitely not the most efficient implementation
!  but it works
!=======================================================================

      SUBROUTINE MRG_PROJ(WDES1,CPROJ,CPROJ_LOCAL)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)    WDES1
      GDEF CPROJ(WDES1%NPRO_TOT)
      GDEF CPROJ_LOCAL(WDES1%NPRO)
#ifdef MPI
      INTEGER NPRO,NT,LMMAXC,NPRO_POS,L,NI,NPRO_TOT
!
! quick copy if possible
      IF (WDES1%COMM_INB%NCPU==1) THEN
!$ACC KERNELS PRESENT(WDES1,CPROJ,CPROJ_LOCAL) __IF_ASYNC__
         CPROJ(1:WDES1%NPRO)=CPROJ_LOCAL(1:WDES1%NPRO)
!$ACC END KERNELS
         RETURN
      ENDIF
!$ACC KERNELS PRESENT(WDES1,CPROJ) __IF_ASYNC__
      CPROJ(1:WDES1%NPRO_TOT)=0
!$ACC END KERNELS

!$ACC PARALLEL LOOP PRESENT(WDES1,WDES1%LMBASE,WDES1%NPRO_POS,WDES1%LMMAX,WDES1%ITYP,CPROJ,CPROJ_LOCAL) &
!$ACC& PRIVATE(NT,LMMAXC,NPRO,NPRO_POS,L) __IF_ASYNC__
      DO NI=1,WDES1%NIONS
         NT = WDES1%ITYP(NI)
         LMMAXC=WDES1%LMMAX(NT)
         IF (LMMAXC/=0) THEN
            NPRO = WDES1%LMBASE(NI)
            NPRO_POS=WDES1%NPRO_POS(NI)
!$ACC LOOP
            DO L=1,LMMAXC
               CPROJ(L+NPRO_POS)=CPROJ_LOCAL(L+NPRO)
            ENDDO
         ENDIF
      ENDDO

      IF (WDES1%LNONCOLLINEAR) THEN

         NPRO_TOT= WDES1%NPRO_TOT/2
!$ACC PARALLEL LOOP PRESENT(WDES1,CPROJ,CPROJ_LOCAL) PRIVATE(NT,LMMAXC,NPRO,NPRO_POS,L) &
!$ACC __IF_ASYNC__
         DO NI=1,WDES1%NIONS
            NT=WDES1%ITYP(NI)
            LMMAXC=WDES1%LMMAX(NT)
            IF (LMMAXC/=0) THEN
               NPRO = WDES1%LMBASE(NI)+WDES1%LMBASE(WDES1%NIONS+1)
               NPRO_POS=WDES1%NPRO_POS(NI)
!$ACC LOOP
               DO L=1,LMMAXC
                  CPROJ(L+NPRO_POS+NPRO_TOT)=CPROJ_LOCAL(L+NPRO)
               ENDDO
            ENDIF
         ENDDO

      ENDIF
      CALLMPI( M_sum_g(WDES1%COMM_INB,CPROJ, WDES1%NPRO_TOT))
#else
!$ACC KERNELS PRESENT(CPROJ,WDES1,CPROJ_LOCAL) __IF_ASYNC__
      CPROJ(1:WDES1%NPRO)=CPROJ_LOCAL(1:WDES1%NPRO)
!$ACC END KERNELS
#endif

      END SUBROUTINE

!=======================================================================
!  distribute plane wave part of one wavefunction over all COMM_INB nodes
!  the supplied wavefunction must have the standard serial layout and
!  must be defined on the node COMM_INB%IONODE
!=======================================================================

      SUBROUTINE DIS_PW(WDES1,CW,CW_LOCAL)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)    WDES1
      COMPLEX(q) CW(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NPL)
#ifdef MPI
      INTEGER NC,IND,I
      INTEGER NPL_TOT,NPL

      CALLMPI( M_bcast_z(WDES1%COMM_INB, CW, WDES1%NPL_TOT))

#ifdef  MPI_barrier_after_bcast
      CALLMPI( MPI_barrier( WDES1%COMM_INB%MPI_COMM, I))
#endif

      IND=0
      DO NC=1,WDES1%NCOL
        DO I=1,WDES1%PL_COL(NC)
          CW_LOCAL(IND+I)=CW(WDES1%PL_INDEX(NC)+I)
        ENDDO
        IND=IND+WDES1%PL_COL(NC)
      ENDDO

      IF ( WDES1%LNONCOLLINEAR ) THEN
         IND=0
         NPL=WDES1%NGVECTOR
         NPL_TOT=WDES1%NPL_TOT / 2
         DO NC=1,WDES1%NCOL
            DO I=1,WDES1%PL_COL(NC)
               CW_LOCAL(IND+I+NPL)=CW(WDES1%PL_INDEX(NC)+I+NPL_TOT)
            ENDDO
            IND=IND+WDES1%PL_COL(NC)
         ENDDO
      ENDIF

      IF (IND*WDES1%NRSPINORS/=WDES1%NPL) THEN
        CALL vtutor%bug("internal error in DIS_PW: " // str(IND) // " " // str(WDES1%NPL), __FILE__, __LINE__)
      ENDIF
#else
      CW_LOCAL(1:WDES1%NPL)=CW(1:WDES1%NPL)
#endif
      END SUBROUTINE


!=======================================================================
!  merge plane wave part of wavefunctions from all COMM_INB nodes
!  using an index array 
!  only used when reading orbitals in Gamma only
!=======================================================================

      SUBROUTINE DIS_PW_GAMMA(WDES1,CW,CW_LOCAL, INDEX, LCONJG)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)    WDES1
      COMPLEX(q) CW(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NPL)
      LOGICAL :: LCONJG(WDES1%NGVECTOR)
      INTEGER :: INDEX(WDES1%NGVECTOR)
#ifdef MPI
      INTEGER NC,IND,I
      INTEGER NPL_TOT,NPL

      CALLMPI( M_bcast_z(WDES1%COMM_INB, CW, WDES1%NPL_TOT))

#ifdef  MPI_barrier_after_bcast
      CALLMPI( MPI_barrier( WDES1%COMM_INB%MPI_COMM, I))
#endif

      DO IND=1,WDES1%NGVECTOR
         IF (LCONJG(IND)) THEN
            CW_LOCAL(IND)=CONJG(CW(INDEX(IND)))
         ELSE
            CW_LOCAL(IND)=CW(INDEX(IND))
         ENDIF
      ENDDO

      IF ( WDES1%LNONCOLLINEAR ) THEN
         NPL=WDES1%NGVECTOR
         NPL_TOT=WDES1%NPL_TOT / 2
         DO IND=1,WDES1%NGVECTOR
            IF (LCONJG(IND)) THEN
               CW_LOCAL(IND+NPL)=CONJG(CW(INDEX(IND)+NPL_TOT))
            ELSE
               CW_LOCAL(IND+NPL)=CW(INDEX(IND)+NPL_TOT)
            ENDIF
         ENDDO
      ENDIF

#else
      CW_LOCAL(1:WDES1%NPL)=CW(1:WDES1%NPL)
#endif
      END SUBROUTINE


!=======================================================================
! distribute one specific band to nodes
! only required after i.e. reading wavefunctions
! it is sufficient if wavefunctions are defined on master node
! but mind, that all nodes must call this communication routine
!=======================================================================

      SUBROUTINE DIS_PW_BAND(WDES1, NB, CW, CW_LOCAL)
      USE prec
      USE wave
      IMPLICIT NONE
      INTEGER NB,NB_L
      TYPE (wavedes1)    WDES1
      COMPLEX(q) CW(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NRPLWV,WDES1%NBANDS)
      INTEGER ierror
#ifdef MPI

! first copy wavefunction to all nodes (it would be sufficient to
!     copy it to master nodes in each inband-communicator, but
!     I am not sure who the master node is)
      CALLMPI( M_bcast_z(WDES1%COMM_KIN, CW, WDES1%NPL_TOT))
#ifdef  MPI_barrier_after_bcast
      CALLMPI( MPI_barrier( WDES1%COMM_KIN%MPI_COMM, ierror ))
#endif

      NB_L=NB_LOCAL(NB,WDES1)
      IF ( MOD(NB-1,WDES1%NB_PAR)+1 == WDES1%NB_LOW) THEN
        CALL DIS_PW(WDES1,CW,CW_LOCAL(1,NB_L))
      ENDIF
#else
      CALL DIS_PW(WDES1,CW,CW_LOCAL(1,NB))
#endif
      END SUBROUTINE

!
! new version uses an index array and works even if parallel or seriell FFT is used
!

      SUBROUTINE DIS_PW_BAND_GAMMA(WDES1, NB, CW, CW_LOCAL, INDEX, LCONJG)
      USE prec
      USE wave
      IMPLICIT NONE
      INTEGER NB,NB_L
      TYPE (wavedes1)    WDES1
      COMPLEX(q) CW(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NRPLWV,WDES1%NBANDS)
      LOGICAL :: LCONJG(WDES1%NGVECTOR)
      INTEGER :: INDEX(WDES1%NGVECTOR)
      INTEGER ierror
#ifdef MPI

! first copy wavefunction to all nodes (it would be sufficient to
!     copy it to master nodes in each inband-communicator, but
!     I am not sure who the master node is)
      CALLMPI( M_bcast_z(WDES1%COMM_KIN, CW, WDES1%NPL_TOT))
#ifdef  MPI_barrier_after_bcast
      CALLMPI( MPI_barrier( WDES1%COMM_KIN%MPI_COMM, ierror ))
#endif

      NB_L=NB_LOCAL(NB,WDES1)
      IF ( MOD(NB-1,WDES1%NB_PAR)+1 == WDES1%NB_LOW) THEN
         CALL DIS_PW_GAMMA(WDES1,CW,CW_LOCAL(1,NB_L), INDEX, LCONJG)
      ENDIF
#else
      ! serial version is always ok
      CALL DIS_PW(WDES1,CW,CW_LOCAL(1,NB))
#endif
      END SUBROUTINE


!=======================================================================
!  merge plane wave part of wavefunctions from all COMM_INB nodes
!=======================================================================

      SUBROUTINE MRG_PW(WDES1,CW,CW_LOCAL)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)    WDES1
      COMPLEX(q) CW(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NPL)
#ifdef MPI
      INTEGER NC,IND,I
      INTEGER NPL_TOT,NPL

      CW=0
      IND=0
      DO NC=1,WDES1%NCOL
        DO I=1,WDES1%PL_COL(NC)
          CW(WDES1%PL_INDEX(NC)+I)=CW_LOCAL(IND+I)
        ENDDO
        IND=IND+WDES1%PL_COL(NC)
      ENDDO

      IF ( WDES1%LNONCOLLINEAR ) THEN
         IND=0
         NPL=WDES1%NGVECTOR
         NPL_TOT=WDES1%NPL_TOT / 2
         DO NC=1,WDES1%NCOL
            DO I=1,WDES1%PL_COL(NC)
               CW(WDES1%PL_INDEX(NC)+I+NPL_TOT)=CW_LOCAL(IND+I+NPL)
            ENDDO
            IND=IND+WDES1%PL_COL(NC)
         ENDDO
      ENDIF

      IF (IND*WDES1%NRSPINORS/=WDES1%NPL) THEN
        CALL vtutor%bug("internal error MRG_PW: " // str(IND) // " " // str(WDES1%NPL), __FILE__, __LINE__)
      ENDIF

      CALLMPI( M_sum_z(WDES1%COMM_INB, CW(1), WDES1%NPL_TOT))
#else
      CW(1:WDES1%NPL)=CW_LOCAL(1:WDES1%NPL)
#endif
      END SUBROUTINE


!=======================================================================
!  merge plane wave part of wavefunctions from all COMM_INB nodes
!  using an index array
!  only used when writing orbitals in Gamma only
!=======================================================================

      SUBROUTINE MRG_PW_GAMMA(WDES1, CW, CW_LOCAL, INDEX, LCONJG)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)    WDES1
      COMPLEX(q) CW(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NPL)
      LOGICAL :: LCONJG(WDES1%NGVECTOR)
      INTEGER :: INDEX(WDES1%NGVECTOR)
#ifdef MPI
      INTEGER IND
      INTEGER NPL_TOT,NPL

      CW=0
      DO IND=1,WDES1%NGVECTOR
         IF (LCONJG(IND)) THEN
            CW(INDEX(IND))=CONJG(CW_LOCAL(IND))
         ELSE
            CW(INDEX(IND))=CW_LOCAL(IND)
         ENDIF
      ENDDO

      IF ( WDES1%LNONCOLLINEAR ) THEN
         NPL=WDES1%NGVECTOR
         NPL_TOT=WDES1%NPL_TOT / 2
         DO IND=1,WDES1%NGVECTOR
            IF (LCONJG(IND)) THEN
               CW(INDEX(IND)+NPL_TOT)=CONJG(CW_LOCAL(IND+NPL))
            ELSE
               CW(INDEX(IND)+NPL_TOT)=CW_LOCAL(IND+NPL)
            ENDIF
         ENDDO
      ENDIF

      CALLMPI( M_sum_z(WDES1%COMM_INB, CW(1), WDES1%NPL_TOT))
#else
      CW(1:WDES1%NPL)=CW_LOCAL(1:WDES1%NPL)
#endif
      END SUBROUTINE

!=======================================================================
! merge  band NB from all nodes into a local copy CW
! only required for i.e. writing wavefunctions
!=======================================================================

      SUBROUTINE MRG_PW_BAND(WDES1, NB, CW, CW_LOCAL)
      USE prec
      USE wave
      IMPLICIT NONE
      INTEGER NB,NB_L
      TYPE (wavedes1)    WDES1
      COMPLEX(q) CW(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NRPLWV,WDES1%NBANDS)
#ifdef MPI
      CW=0
      NB_L=NB_LOCAL(NB,WDES1)
! if wavefunction is located on this node merge it to CW
      IF ( MOD(NB-1,WDES1%NB_PAR)+1 == WDES1%NB_LOW) THEN
        CALL MRG_PW(WDES1,CW,CW_LOCAL(1,NB_L))
      ENDIF
! then merge over all nodes using inter-band communication
      CALLMPI( M_sum_z(WDES1%COMM_INTER, CW(1), WDES1%NPL_TOT))
#else
      CALL MRG_PW(WDES1,CW,CW_LOCAL(1,NB))
#endif
      END SUBROUTINE

!
! new version uses an index array and works even if parallel or seriell FFT is used
!

      SUBROUTINE MRG_PW_BAND_GAMMA(WDES1, NB, CW, CW_LOCAL, INDEX, LCONJG)
      USE prec
      USE wave
      IMPLICIT NONE
      INTEGER NB,NB_L

      TYPE (wavedes1)    WDES1
      COMPLEX(q) CW(WDES1%NPL_TOT)
      COMPLEX(q) CW_LOCAL(WDES1%NRPLWV,WDES1%NBANDS)
      LOGICAL :: LCONJG(WDES1%NGVECTOR)
      INTEGER :: INDEX(WDES1%NGVECTOR)
#ifdef MPI
      CW=0
      NB_L=NB_LOCAL(NB,WDES1)
! if wavefunction is located on this node merge it to CW
      IF ( MOD(NB-1,WDES1%NB_PAR)+1 == WDES1%NB_LOW) THEN
         CALL MRG_PW_GAMMA(WDES1,CW,CW_LOCAL(1,NB_L), INDEX, LCONJG)
      ENDIF
! then merge over all nodes using inter-band communication
      CALLMPI( M_sum_z(WDES1%COMM_INTER, CW(1), WDES1%NPL_TOT))
#else
      ! serial version is always ok
      CALL MRG_PW(WDES1,CW,CW_LOCAL(1,NB))
#endif
      END SUBROUTINE

!=======================================================================
! merge eigenvalues from all nodes
!=======================================================================

      SUBROUTINE MRG_CEL(WDES,W)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)    WDES
      TYPE (wavedes1)   WDES1
      TYPE (wavespin)   W
      INTEGER I,NK,NB_GLOBAL,NB,NCEL

      PUSH_ACC_EXEC_ON(.FALSE.)
#ifdef MPI
!
! first zero out all components from other nodes
!
      spin:   DO I=1,WDES%ISPIN
      kpoint: DO NK=1,WDES%NKPTS
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
         DO NB_GLOBAL=1,WDES%NB_TOT
            W%CELTOT(NB_GLOBAL,NK,I)=0
         END DO
      ELSE
         CALL SETWDES(WDES,WDES1,0)
         band: DO NB_GLOBAL=1,WDES%NB_TOT
            NB=NB_LOCAL(NB_GLOBAL,WDES1)
            IF (NB==0) THEN
               W%CELTOT(NB_GLOBAL,NK,I)=0
            ENDIF
         ENDDO band
      ENDIF
      ENDDO kpoint
      ENDDO spin
!
! then merge CELTOT from all nodes
!
      NCEL=WDES%NB_TOT*WDES%NKPTS*WDES%ISPIN
      CALLMPI( M_sum_z( WDES%COMM_KINTER, W%CELTOT(1,1,1), NCEL))
      CALLMPI( M_sum_z( WDES%COMM_INTER,  W%CELTOT(1,1,1), NCEL))
#endif
      ! finally set eigenvalues of states beyond WDES%NB_TOTK to large values
      CALL WFSET_HIGH_CELEN(W)

      POP_ACC_EXEC_ON

      END SUBROUTINE


      SUBROUTINE MRG_CELTOT(WDES,W)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)    WDES
      TYPE (wavedes1)   WDES1
      TYPE (wavespin)   W
      INTEGER I,NK,NB_GLOBAL,NB,NCEL
!$ACC ROUTINE(NB_LOCAL) SEQ
#ifdef MPI
!
! first zero out all components from other nodes
!
!$ACC ENTER DATA CREATE(WDES1) __IF_ASYNC__
      CALL SETWDES(WDES,WDES1,0)

!$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(WDES,WDES1,W) __IF_ASYNC__
      spin:   DO I=1,WDES%ISPIN
      kpoint: DO NK=1,WDES%NKPTS
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
!$ACC LOOP VECTOR
         DO NB_GLOBAL=1,WDES%NB_TOT
            W%CELTOT(NB_GLOBAL,NK,I)=0
         END DO
      ELSE
!$ACC LOOP VECTOR
         band: DO NB_GLOBAL=1,WDES%NB_TOT
            NB=NB_LOCAL(NB_GLOBAL,WDES1)
            IF (NB==0) THEN
               W%CELTOT(NB_GLOBAL,NK,I)=0
            ENDIF
         ENDDO band
      ENDIF
      ENDDO kpoint
      ENDDO spin

#ifdef _OPENACC
    CALL ACC_DELETE_TYPED_VAR(WDES1)
#endif
!
! then merge CELTOT from all nodes
!
      NCEL=WDES%NB_TOT*WDES%NKPTS*WDES%ISPIN
      CALLMPI( M_sum_z( WDES%COMM_KINTER, W%CELTOT(1,1,1), NCEL))
      CALLMPI( M_sum_z( WDES%COMM_INTER,  W%CELTOT(1,1,1), NCEL))
#endif

      ! finally set eigenvalues of states beyond WDES%NB_TOTK to large values
      CALL WFSET_HIGH_CELEN(W)

      END SUBROUTINE


!=======================================================================
! merge partial occupancies from all nodes
!=======================================================================

      SUBROUTINE MRG_FERWE(WDES,W)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)    WDES
      TYPE (wavedes1)   WDES1
      TYPE (wavespin)   W
      INTEGER I,NK,NB_GLOBAL,NB,NCEL
!$ACC ROUTINE(NB_LOCAL) SEQ
#ifdef MPI
!
! first zero out all components from other nodes
!
!$ACC ENTER DATA CREATE(WDES1) __IF_ASYNC__
      CALL SETWDES(WDES,WDES1,0)

!$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(WDES,WDES1,W) __IF_ASYNC__
      spin:   DO I=1,WDES%ISPIN
      kpoint: DO NK=1,WDES%NKPTS
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
!$ACC LOOP VECTOR
         DO NB_GLOBAL=1,WDES%NB_TOT
            W%FERTOT(NB_GLOBAL,NK,I)=0
         END DO
      ELSE
!$ACC LOOP VECTOR
         band: DO NB_GLOBAL=1,WDES%NB_TOT
           NB=NB_LOCAL(NB_GLOBAL,WDES1)
           IF (NB==0) THEN
             W%FERTOT(NB_GLOBAL,NK,I)=0
           ENDIF
         ENDDO band
      ENDIF
      ENDDO kpoint
      ENDDO spin

#ifdef _OPENACC
    CALL ACC_DELETE_TYPED_VAR(WDES1)
#endif
!
! then merge FERTOT from all nodes
!
      NCEL=WDES%NB_TOT*WDES%NKPTS*WDES%ISPIN
      CALLMPI( M_sum_d( WDES%COMM_KINTER, W%FERTOT(1,1,1), NCEL))
      CALLMPI( M_sum_d( WDES%COMM_INTER,  W%FERTOT(1,1,1), NCEL))
#endif
      END SUBROUTINE


!=======================================================================
! merge auxilary array from all nodes
!=======================================================================

      SUBROUTINE MRG_AUX(WDES,W)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)    WDES
      TYPE (wavedes1)   WDES1
      TYPE (wavespin)   W
      INTEGER I,NK,NB_GLOBAL,NB,NCEL
!$ACC ROUTINE(NB_LOCAL) SEQ
#ifdef MPI
!
! first zero out all components from other nodes
!
!$ACC ENTER DATA CREATE(WDES1) __IF_ASYNC__
      CALL SETWDES(WDES,WDES1,0)

!$ACC PARALLEL LOOP COLLAPSE(2) GANG PRESENT(WDES,WDES1,W) __IF_ASYNC__
      spin:   DO I=1,WDES%ISPIN
      kpoint: DO NK=1,WDES%NKPTS
      IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
!$ACC LOOP VECTOR
         DO NB_GLOBAL=1,WDES%NB_TOT
            W%AUXTOT(NB_GLOBAL,NK,I)=0
         END DO
      ELSE
!$ACC LOOP VECTOR
         band: DO NB_GLOBAL=1,WDES%NB_TOT
           NB=NB_LOCAL(NB_GLOBAL,WDES1)
           IF (NB==0) THEN
             W%AUXTOT(NB_GLOBAL,NK,I)=0
           ENDIF
         ENDDO band
      ENDIF
      ENDDO kpoint
      ENDDO spin

#ifdef _OPENACC
    CALL ACC_DELETE_TYPED_VAR(WDES1)
#endif
!
! then merge AUXTOT from all nodes
!
      NCEL=WDES%NB_TOT*WDES%NKPTS*WDES%ISPIN
      CALLMPI( M_sum_d( WDES%COMM_KINTER, W%AUXTOT(1,1,1), NCEL))
      CALLMPI( M_sum_d( WDES%COMM_INTER,  W%AUXTOT(1,1,1), NCEL))
#endif
      END SUBROUTINE

!************************ SUBROUTINE REDIS_PW **************************
!
! redistribute plane wave coefficients from over band to over
! plane wave coefficient or vice versa
! this operation is done in place to reduce storage demands
!
! mind that if the routine is called twice the original distribution
! is obtained
!
!***********************************************************************

      SUBROUTINE REDIS_PW(WDES1, NBANDS, CW)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)  WDES1
      INTEGER NBANDS
      COMPLEX(q) :: CW(WDES1%NRPLWV*NBANDS)
#if defined(T3D_SMA)
      INTEGER MALLOC_DONE
      INTEGER, EXTERNAL :: ISHM_CHECK
      COMMON /SHM/ MALLOC_DONE, PBUF
      POINTER ( PBUF, CWORK ) ;  COMPLEX(q):: CWORK(*)
#else
      COMPLEX(q), ALLOCATABLE :: CWORK(:)
#endif
! local variables
      INTEGER :: NRPLWV,N,NB,INFO
#ifdef _OPENACC
      LOGICAL :: ACC_ACTIVE
      ACC_ACTIVE=ACC_EXEC_ON.AND.ACC_IS_PRESENT(CW)
#endif
#ifdef MPI
      PROFILING_START('redis_pw')

      ! quick return if possible
      IF (WDES1%COMM_INTER%NCPU ==1) THEN
         PROFILING_STOP('redis_pw')
         RETURN
      ENDIF

      NRPLWV=WDES1%NRPLWV

      IF (MOD(NRPLWV,WDES1%NB_PAR) /= 0) THEN
        CALL vtutor%bug("REDIS_PW: internal error(1) " // str(NRPLWV) // " " // str(WDES1%NB_PAR), __FILE__, __LINE__)
      ENDIF
#ifdef T3D_SMA
      IF (ISHM_CHECK(NRPLWV*2) == 0) THEN
         CALL vtutor%error("REDIS_PW: shmem workspace unsufficient")
      ENDIF
#else
      ALLOCATE( CWORK(NRPLWV))
DOESI CWORK=(0.0_q, 0.0_q)
!$ACC DATA CREATE(CWORK) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
#endif
      DO NB=0,NBANDS-1
         CALL M_alltoall_z(WDES1%COMM_INTER, NRPLWV, CW(NB*NRPLWV+1), CWORK(1))
         CALL __ZCOPY__( NRPLWV, CWORK(1), 1, CW(NB*NRPLWV+1), 1)
      ENDDO
#ifndef T3D_SMA
!$ACC WAIT(ACC_ASYNC_Q) IF(ACC_ACTIVE)
!$ACC END DATA
      DEALLOCATE(CWORK)
#endif
      PROFILING_STOP('redis_pw')
#endif
      END SUBROUTINE

!************************ SUBROUTINE REDIS_PW_OVER_BANDS ***************
!
! redistribute all coefficients if and only if LOVER_BAND
! is .TRUE.
! after this call the bands are distributed over nodes
! and LOVER_BAND is .FALSE.
!
!***********************************************************************

      SUBROUTINE REDIS_PW_OVER_BANDS(WDES, W)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)     WDES
      TYPE (wavespin)    W

#ifdef MPI
      PROFILING_START('redis_pw_over_bands')

      IF ( W%OVER_BAND ) THEN
         CALL REDIS_PW_ALL(WDES, W)
      ENDIF

      PROFILING_STOP('redis_pw_over_bands')
#endif
      END SUBROUTINE REDIS_PW_OVER_BANDS

!************************ SUBROUTINE REDIS_PW_ALL **********************
!
! redistribute all plane wave coeffcients, change the order 
! W%OVER_BAND
! calling the routine twice restores the original order
!
!***********************************************************************

      SUBROUTINE REDIS_PW_ALL(WDES, W)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)     WDES
      TYPE (wavespin)    W

#ifdef MPI
      TYPE (wavedes1)    WDES1          ! descriptor for one k-point
      INTEGER NCPU, ISP, NK

      PROFILING_START('redis_pw_all')
     
      NCPU   =WDES%COMM_INTER%NCPU ! number of procs involved in band dis.
      IF (NCPU /= 1) THEN
         DO ISP=1,WDES%ISPIN
         DO NK=1,WDES%NKPTS
            IF (MOD(NK-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) CYCLE
            CALL SETWDES(WDES,WDES1,NK)
            CALL REDIS_PW  (WDES1, WDES%NBANDS, W%CW   (1,1,NK,ISP))
         ENDDO
         ENDDO

         W%OVER_BAND=.NOT. W%OVER_BAND
      ENDIF

      PROFILING_STOP('redis_pw_all')
#endif
      END SUBROUTINE REDIS_PW_ALL

!************************ SUBROUTINE REDIS_PW_START ********************
!
! redistribute plane wave coefficients from over band to over
! plane wave coefficient or vice versa asynchronously 
! one band is initiated
!
!***********************************************************************

      SUBROUTINE REDIS_PW_START(WDES, CW, BANDINDEX, H)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE wave_mpi
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)  WDES
      INTEGER BANDINDEX
      COMPLEX(q) :: CW(WDES%NRPLWV)
      TYPE(redis_pw_ctr) :: H
! local variables
      INTEGER :: NRPLWV
      INTEGER :: IND,N

#ifdef MPI
      PROFILING_START('redis_pw_start')

      IF (BANDINDEX>H%NBANDS) THEN
         WRITE(0,*) 'REDIS_PW_START: internal error, non existing band',BANDINDEX
      ENDIF
      NRPLWV=WDES%NRPLWV

      IND=H%NEXT
      H%NEXT=H%NEXT+1 ; IF (H%NEXT > H%NB) H%NEXT=1

      IF (H%BAND(IND)/=0) THEN
         CALL vtutor%bug("REDIS_PW_START: internal error, slot NEXT not empty " // str(H%NEXT), __FILE__, __LINE__)
      ENDIF
      ! fill in band index
      H%BAND(IND)=BANDINDEX
#ifdef no_async
!$ACC PARALLEL LOOP PRESENT(CW,H%CW) __IF_ASYNC__
      DO N=1,NRPLWV
         H%CW(N,IND)=CW(N)
      ENDDO
#else
      CALL M_alltoall_d_async(H%COMM, 2*NRPLWV, CW(1), H%CW(1,IND), &
           H%ICOMM+IND, H%SREQUEST(1,IND), H%RREQUEST(1,IND) )
#endif
      PROFILING_STOP('redis_pw_start')
#endif
      END SUBROUTINE

!************************ SUBROUTINE REDIS_PW_STOP  ********************
!
! finish redistribution of  plane wave coefficients from over band to
! over plane wave coefficient or vice versa.
!
!***********************************************************************

      SUBROUTINE REDIS_PW_STOP(WDES, CW, BANDINDEX, H)
      USE prec
      USE wave_mpi
      USE wave
      IMPLICIT NONE

      TYPE (wavedes)  WDES
      INTEGER BANDINDEX
      COMPLEX(q) :: CW(WDES%NRPLWV)
      TYPE(redis_pw_ctr) :: H
! local variables
      INTEGER :: NRPLWV
      INTEGER :: IND,N
#ifdef MPI
      PROFILING_START('redis_pw_stop')

      IF (BANDINDEX>H%NBANDS) THEN
         WRITE(0,*) 'REDIS_PW_STOP: internal error, non existing band',BANDINDEX
      ENDIF

      NRPLWV=WDES%NRPLWV

      ! usually the required band is in the slot NEXT, however
      ! possibly we have to search for it

      IND=H%NEXT
      IF (H%BAND(IND) /= BANDINDEX) THEN
         DO IND=1, H%NB
            IF ( H%BAND(IND) == BANDINDEX) EXIT
         ENDDO
         IND=MIN(IND,H%NB)
         IF  (H%BAND(IND) /= BANDINDEX) THEN
            CALL vtutor%bug("REDIS_PW_STOP: internal error, can not find " // str(BANDINDEX) // " " &
               // str(H%BAND), __FILE__, __LINE__)
         ENDIF
      ENDIF
      ! fill in zero into the band index array
      H%BAND(IND)=0
#ifdef no_async
      CALL M_alltoall_d_async(H%COMM, 2*NRPLWV, H%CW(1,IND), CW(1) , &
           H%ICOMM+IND, H%SREQUEST(1,IND), H%RREQUEST(1,IND) )
      CALL M_alltoall_wait(H%COMM, H%SREQUEST(1,IND), H%RREQUEST(1,IND) )
#else
      ! wait for send and receive to finish
      CALL M_alltoall_wait(H%COMM, H%SREQUEST(1,IND), H%RREQUEST(1,IND) )
!$ACC PARALLEL LOOP PRESENT(CW,H%CW) __IF_ASYNC__
      DO N=1,NRPLWV
         CW(N)=H%CW(N,IND)
      ENDDO
#endif
      PROFILING_STOP('redis_pw_stop')
#endif
      END SUBROUTINE


!************************ SUBROUTINE REDIS_PROJ ************************
!
! redistribute projector part of wave function from over band to over
! plane wave coefficient or vice versa
! this operation is done in place to reduce storage demands
!
! mind that if the routine is called twice the original distribution
! is obtained
!
!***********************************************************************

      SUBROUTINE REDIS_PROJ(WDES1, NBANDS, CPROJ)
#ifdef _OPENACC
      USE mopenacc
#endif
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)  WDES1
      INTEGER NBANDS
      GDEF :: CPROJ(WDES1%NPROD*NBANDS)
#if defined(T3D_SMA)
      INTEGER MALLOC_DONE
      INTEGER, EXTERNAL :: ISHM_CHECK
      COMMON /SHM/ MALLOC_DONE, PBUF
      POINTER ( PBUF, CWORK ) ;  GDEF:: CWORK(*)
#else
      GDEF, ALLOCATABLE :: CWORK(:)
#endif

#ifdef gammareal
      INTEGER, PARAMETER :: MCOMP=1
#else
      INTEGER, PARAMETER :: MCOMP=2
#endif
! local variables
#ifdef MPI
      INTEGER :: NPROD,N,NB,INFO
#ifdef _OPENACC
      LOGICAL :: ACC_ACTIVE
      ACC_ACTIVE=ACC_EXEC_ON.AND.ACC_IS_PRESENT(CPROJ)
#endif
      PROFILING_START('redis_proj')

      IF (WDES1%COMM_INTER%NCPU==1 .OR. WDES1%NPROD==0) THEN
         PROFILING_STOP('redis_proj')
         RETURN
      ENDIF

      NPROD=WDES1%NPROD

      IF (MOD(NPROD,WDES1%NB_PAR) /= 0) THEN
        CALL vtutor%bug("REDIS_PW: NPROD not dividable by NB_PAR:" // str(NPROD) // " " // str(WDES1%NB_PAR), __FILE__, __LINE__)
      ENDIF
#ifdef T3D_SMA
      IF (ISHM_CHECK(NPROD*MCOMP) == 0) THEN
         CALL vtutor%bug("REDIS_PROJ: insufficient shmem workspace", __FILE__, __LINE__)
      ENDIF
#else
      ALLOCATE( CWORK(NPROD))
!$ACC ENTER DATA CREATE(CWORK) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
DOESI CWORK=0
#endif
      DO NB=0,NBANDS-1
#ifdef gammareal
         CALL M_alltoall_d(WDES1%COMM_INTER, NPROD, CPROJ(NB*NPROD+1), CWORK(1))
#else
         CALL M_alltoall_z(WDES1%COMM_INTER, NPROD, CPROJ(NB*NPROD+1), CWORK(1))
#endif
         CALL __GCOPY__( NPROD, CWORK(1), 1, CPROJ(NB*NPROD+1), 1)
      ENDDO
#ifndef T3D_SMA
!$ACC EXIT DATA DELETE(CWORK) IF(ACC_ACTIVE) ASYNC(ACC_ASYNC_Q)
      DEALLOCATE(CWORK)
#endif
      PROFILING_STOP('redis_proj')
#endif
      END SUBROUTINE


!************************ SUBROUTINE SET_NPL_NPRO **********************
!
! set the local number of plane waves on a node after redistribution
! of plane wave coefficients and projected wavefunctions
! on entry the global number of plane waves must be given
!
!***********************************************************************

      SUBROUTINE SET_NPL_NPRO(WDES1, NPL, NPRO)
      USE prec
      USE wave
      IMPLICIT NONE

      TYPE (wavedes1)  WDES1
      INTEGER NPL,NPRO,NPL_REM,NPRO_REM,I
! local variables
      NPL_REM =NPL
      NPRO_REM=NPRO
#ifdef MPI
      NPL =WDES1%NRPLWV/WDES1%COMM_INTER%NCPU
      NPRO=WDES1%NPROD /WDES1%COMM_INTER%NCPU
      DO I=1,WDES1%COMM_INTER%NODE_ME
        NPL =MIN(NPL,NPL_REM)
        NPRO=MIN(NPRO,NPRO_REM)
        NPL_REM  =NPL_REM -NPL
        NPRO_REM =NPRO_REM-NPRO
      ENDDO
#else
      NPL=WDES1%NRPLWV
      NPRO=WDES1%NPROD
#endif

      END SUBROUTINE

!***********************************************************************
!
! Copy wavefunction coefficients from each k-point communicator to
! replicated copies in other communicators
! Should use a gather operation
!
!***********************************************************************

    SUBROUTINE KPAR_SYNC_WAVEFUNCTIONS(WDES,W)
      USE prec
      USE wave
      USE wave_mpi
      USE mpimy
      IMPLICIT NONE
      TYPE (wavedes)     WDES
      TYPE (wavespin)    W
      INTEGER :: ISP,K,inode

#ifdef MPI
      IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
         DO ISP=1,WDES%ISPIN
            DO K=1,WDES%NKPTS
               inode=MOD(K-1,WDES%COMM_KINTER%NCPU)+1
! the product of the number of bands and number of orbitals is possible exceeding INTEGER (KIND=4)
               CALLMPI( M_bcast_z8_from(WDES%COMM_KINTER, W%CW(1,1,K,ISP),SIZE(W%CW,1,KIND=selected_int_kind(15))*SIZE(W%CW,2,KIND=selected_int_kind(15)),inode) )
#ifdef gammareal
               CALLMPI( M_bcast_d8_from(WDES%COMM_KINTER, W%CPROJ(1,1,K,ISP),SIZE(W%CPROJ,1,KIND=selected_int_kind(15))*SIZE(W%CPROJ,2,KIND=selected_int_kind(15)),inode) )
#else
               CALLMPI( M_bcast_z8_from(WDES%COMM_KINTER, W%CPROJ(1,1,K,ISP),SIZE(W%CPROJ,1,KIND=selected_int_kind(15))*SIZE(W%CPROJ,2,KIND=selected_int_kind(15)),inode) )
#endif
            END DO
         END DO
      ENDIF
#endif
    END SUBROUTINE KPAR_SYNC_WAVEFUNCTIONS

!***********************************************************************
!
! Copy eigenvalues from each k-point communicator to
! replicated copies in other communicators
!
!***********************************************************************

    SUBROUTINE KPAR_SYNC_CELTOT(WDES,W)
! Not the same as a MRG_CEL since we do not sum over within bands
      USE prec
      USE wave
      USE wave_mpi
      USE mpimy
      IMPLICIT NONE
      TYPE (wavedes) :: WDES
      TYPE (wavespin) :: W
      INTEGER :: ISP,K

#ifdef MPI
      IF (WDES%COMM_KINTER%NCPU.GT.1) THEN

         PUSH_ACC_EXEC_ON(.FALSE.)

         DO ISP=1,WDES%ISPIN
            DO K=1,WDES%NKPTS
               IF (MOD(K-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
                  W%CELTOT(:,K,ISP)=0.
               END IF
            END DO
         END DO
! Reduce and broadcast for simplicity. Really an allgather.
         CALLMPI( M_sum_master_z( WDES%COMM_KINTER, W%CELTOT, SIZE(W%CELTOT) ) )
         CALLMPI( M_bcast_z( WDES%COMM_KINTER, W%CELTOT, SIZE(W%CELTOT) ) )

         POP_ACC_EXEC_ON

      ENDIF
#endif
    END SUBROUTINE KPAR_SYNC_CELTOT

!***********************************************************************
!
! Copy occupancies from each k-point communicator to
! replicated copies in other communicators
!
!***********************************************************************

    SUBROUTINE KPAR_SYNC_FERTOT(WDES,W)
      USE prec
      USE wave
      USE wave_mpi
      USE mpimy
      IMPLICIT NONE
      TYPE (wavedes) :: WDES
      TYPE (wavespin) :: W
      INTEGER :: ISP,K

#ifdef MPI
      IF (WDES%COMM_KINTER%NCPU.GT.1) THEN

         PUSH_ACC_EXEC_ON(.FALSE.)

         DO ISP=1,WDES%ISPIN
            DO K=1,WDES%NKPTS
               IF (MOD(K-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
                  W%FERTOT(:,K,ISP)=0.
               END IF
            END DO
         END DO
! Reduce and broadcast for simplicity. Really an allgather.
         CALLMPI( M_sum_master_d( WDES%COMM_KINTER, W%FERTOT, SIZE(W%FERTOT) ) )
         CALLMPI( M_bcast_d( WDES%COMM_KINTER, W%FERTOT, SIZE(W%FERTOT) ) )

         POP_ACC_EXEC_ON

      ENDIF
#endif
    END SUBROUTINE KPAR_SYNC_FERTOT

!***********************************************************************
!
! Copy auxilary array from each k-point communicator to
! replicated copies in other communicators
!
!***********************************************************************

    SUBROUTINE KPAR_SYNC_AUXTOT(WDES,W)
      USE prec
      USE wave
      USE wave_mpi
      USE mpimy
      IMPLICIT NONE
      TYPE (wavedes) :: WDES
      TYPE (wavespin) :: W
      INTEGER :: ISP,K

#ifdef MPI
      IF (WDES%COMM_KINTER%NCPU.GT.1) THEN

         PUSH_ACC_EXEC_ON(.FALSE.)

         DO ISP=1,WDES%ISPIN
            DO K=1,WDES%NKPTS
               IF (MOD(K-1,WDES%COMM_KINTER%NCPU).NE.WDES%COMM_KINTER%NODE_ME-1) THEN
                  W%AUXTOT(:,K,ISP)=0.
               END IF
            END DO
         END DO
! Reduce and broadcast for simplicity. Really an allgather.
         CALLMPI( M_sum_master_d( WDES%COMM_KINTER, W%AUXTOT, SIZE(W%AUXTOT) ) )
         CALLMPI( M_bcast_d( WDES%COMM_KINTER, W%AUXTOT, SIZE(W%AUXTOT) ) )

         POP_ACC_EXEC_ON

      ENDIF
#endif
    END SUBROUTINE KPAR_SYNC_AUXTOT

!***********************************************************************
!
! Sync orbitals, eigenvalues and occupancies
!
!***********************************************************************

    SUBROUTINE KPAR_SYNC_ALL(WDES,W)
      USE wave
      IMPLICIT NONE
      TYPE (wavedes) :: WDES
      TYPE (wavespin) :: W
#ifdef MPI
      IF (WDES%COMM_KINTER%NCPU.GT.1) THEN
         CALL KPAR_SYNC_FERTOT(WDES,W)
         CALL KPAR_SYNC_CELTOT(WDES,W)
         CALL KPAR_SYNC_AUXTOT(WDES,W)
         CALL KPAR_SYNC_WAVEFUNCTIONS(WDES,W)
      ENDIF
#endif
    END SUBROUTINE KPAR_SYNC_ALL
