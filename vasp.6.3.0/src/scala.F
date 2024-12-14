#include "symbol.inc"
!=======================================================================
! RCS:  $Id: scala.F,v 1.5 2003/06/27 13:22:22 kresse Exp kresse $
!
! Module containing wrapper for scaLAPACK
! written by Gilles de Wijs (gD) and Georg Kresse (gK)
! modified to run on any number of nodes by Dario Alfe
! BG_ routines were added by IBM fellow (bG)
!
! 10.09.2014   gK:
! - removed the T3D support (broken anyway most likely)
! - BG_ routines are now always compiled in
! - added some comments
! - added double check that BG_INIT_scala was called before calling
!    compute routines (simply not sure this is done properly in VASP)
! - removed a lot of unused variables
! - some clean up, but a lot of off and creapy stuff is still around
!   specifically the redundance between BG_INIT_SCALA and INIT_SCALA
!   is odd
!
!=======================================================================

#ifndef scaLAPACK


 MODULE scala
   USE prec
   USE mpimy

   IMPLICIT NONE

   LOGICAL, PUBLIC :: LscaLAPACK = .FALSE.
   LOGICAL, PUBLIC :: LscaLU     = .FALSE.  ! use parallel LU  decomposition
   LOGICAL, PUBLIC :: LscaAWARE  = .FALSE.  ! do not store sub space matrix
   LOGICAL, PUBLIC :: LscaAWARE_read  = .FALSE.  ! LscaAWARE was read from file

  CONTAINS
    SUBROUTINE  INIT_SCALAAWARE( NB_TOT, NRPLWV, COMM )
      INTEGER :: NB_TOT, NRPLWV
      TYPE (communic) COMM
    END SUBROUTINE INIT_SCALAAWARE

    SUBROUTINE pPOTRF_TRTRI (COMM, AMATIN, NDIM, N)
      TYPE (communic) COMM
      INTEGER NDIM,N            ! NxN matrix to be distributed
      GDEF    AMATIN(NDIM,N)    ! input/output matrix
    END SUBROUTINE pPOTRF_TRTRI

    SUBROUTINE pPOTRF_POTRI (COMM, AMATIN, NDIM, N)
      TYPE (communic) COMM
      INTEGER NDIM,N            ! NxN matrix to be distributed
      GDEF    AMATIN(NDIM,N)    ! input/output matrix
    END SUBROUTINE pPOTRF_POTRI

    SUBROUTINE pDSSYEX_ZHEEVX(COMM,AMATIN,W, NDIM, N, ABSTOL, AMATIN_SCALA)
      TYPE (communic) COMM
      INTEGER NDIM,N                        ! NxN matrix to be distributed
      GDEF    AMATIN(NDIM,N)                ! input/output matrix
      REAL(q) W(N)                          ! eigenvalues
      REAL(q), OPTIONAL ::  ABSTOL          ! tolerance
      GDEF, OPTIONAL, TARGET ::  AMATIN_SCALA(N*N) ! alternative input/output matrix (scaLAPACK format)
    END SUBROUTINE pDSSYEX_ZHEEVX

    SUBROUTINE BG_pDSSYEX_ZHEEVX(COMM,A, W, N)
      TYPE (communic) COMM
      INTEGER N                             ! NxN matrix to be distributed
      GDEF    A(1,1)                        ! input/output matrix
      REAL(q) W(N)                          ! eigenvalues
    END SUBROUTINE BG_pDSSYEX_ZHEEVX

    SUBROUTINE BG_pPOTRF_TRTRI(A, N, INFO)
      TYPE (communic) COMM
      INTEGER N          ! NxN matrix to be distributed
      GDEF    A(1,1)     ! input/output matrix
      INTEGER INFO
    END SUBROUTINE BG_pPOTRF_TRTRI

    SUBROUTINE BG_pPOTRF_POTRI(A, N, INFO)
      TYPE (communic) COMM
      INTEGER N          ! NxN matrix to be distributed
      GDEF    A(1,1)     ! input/output matrix
      INTEGER INFO
    END SUBROUTINE BG_pPOTRF_POTRI

    SUBROUTINE BG_CHANGE_DIAGONALE(N,A,IU0)
      INTEGER N,IU0
      GDEF    A(*)
    END SUBROUTINE BG_CHANGE_DIAGONALE

    SUBROUTINE INIT_scala(COMM, N)
      IMPLICIT NONE
      TYPE (communic) :: COMM
      INTEGER N
    END SUBROUTINE INIT_scala

    FUNCTION SCALA_NP()
      INTEGER SCALA_NP
      SCALA_NP=1
    END FUNCTION SCALA_NP

    FUNCTION SCALA_NQ()
      INTEGER SCALA_NQ
      SCALA_NQ=1
    END FUNCTION SCALA_NQ

    SUBROUTINE INIT_scala_t3d
    END SUBROUTINE INIT_scala_t3d

 END MODULE

#else

 MODULE scala
    USE prec
    USE mpimy
#ifdef ELPA
    USE elpa
    LOGICAL, PUBLIC :: LELPA      = .TRUE.    ! use ELPA
    INTEGER,SAVE,PRIVATE :: ELPA_MPI_COMM_ROWS,ELPA_MPI_COMM_COLS,ELPA_MPI_COMM_ALL
#else
    LOGICAL, PUBLIC :: LELPA      = .FALSE.   ! can't use ELPA
#endif

    LOGICAL, PUBLIC :: LscaLAPACK = .TRUE.
    LOGICAL, PUBLIC :: LscaLU     = .FALSE.   ! use parallel LU  decomposition
    LOGICAL, PUBLIC :: LscaAWARE  = .FALSE.   ! use parallel LU  decomposition
    LOGICAL, PUBLIC :: LscaAWARE_read  = .FALSE. ! LscaAWARE was read from file

!
! customize if required
!
      INTEGER,PRIVATE,PARAMETER :: NB=16      ! blocking factor for distribution of matrices
      ! P4 optimal, larger values were even slightly better (160) but still
      ! slower on a Gigabit cluster than ZHEEVX

      ! customization of matrix diagonalization

      INTEGER,PRIVATE :: NCLUST=48            ! maximum cluster of eigenvector
      REAL(q),PRIVATE :: ABSTOL_DEF=1e-10_q   ! specifies eigenvector orthogonality tolerance
      REAL(q),PRIVATE :: ORFAC=-1.e0_q        ! controls reorthogonalisation of eigenvectors
      REAL(qs),PRIVATE:: ORFAC_single=-1.e0_qs! controls reorthogonalisation of eigenvectors

      !  see manpage for pSSYEVX pZHEEVX:
      !    NCLUST is the number of eigenvalues in the largest cluster,
      !    where a cluster is defined as a set of close eigenvalues:
      !    {W(K),...,W(K+NCLUST-1)|W(J+1)<= W(J)+orfac*norm(A)}
!
! end customization
!
      INTEGER,PARAMETER,PRIVATE :: DTYPE_=1, &
                           CTXT_=2, M_=3, N_=4, MB_=5, NB_=6,  &
                           RSRC_=7, CSRC_ =8, LLD_=9, DLEN_=9

! DESCA( DTYPE_ )The descriptor type.  In most cases,
!                DTYPE_ = BLOCK_CYCLIC_2D = 1
! DESCA( CTXT_ ) The BLACS context handle, indicating
!                the BLACS process grid A is distribu-
!                ted over. The context itself is glo-
!                bal, but the handle (the integer
!                value) may vary.
! DESCA( M_ )    The number of rows in the global
!                array A.
! DESCA( N_ )    The number of columns in the global
!                array A.
! DESCA( MB_ )   The blocking factor used to distribute
!                the rows of the array.
! DESCA( NB_ )   The blocking factor used to distribute
!                the columns of the array.
! DESCA( RSRC_ ) The process row over which the first
!                row of the array A is distributed.
! DESCA( CSRC_ ) The process column over which the
!                first column of the array A is
!                distributed.
! DESCA( LLD_ )  The leading dimension of the local
!                array.  LLD_A >= MAX(1,LOCr(M_A)).
! DLEN_          dimension of DESCA

!
! INIT_scala initializes the following variables
!
! usually VASP uses a "light-weighted" implementation of scalapack, in the
! sense that the matrices are stored in the conventional seriel manner
! and only when the scaLAPACK routines are called, the matrices are
! stored distributed and then treated by scaLAPACK
! this has some disatvantages, for instance the storage requirements
! are huge
! if the variable  LscaAWARE is set in the INCAR, however
! VASP uses distributed matrices from the outside
! however, even in this case, scaLAPACK is essentially completely
! "unvisible" to VASP, and all scaLAPACK related routines
! should go in here
!
! in several place (bse.F, chi_GG.F) that strategy did not work out,
! and the implementations is directly in the main VASP code
!
! most routines in this module use a global array descriptor DESCSTD
! which is set up by INIT_SCALA. INIT_SCALA can be called
! repeatedly if the matrix dimensions or the processor grid change
! however, all subroutines that have an argument DESCA
! are save and they can be used in different contexts
!
      TYPE scalapack_des
         INTEGER :: MAGIC=0        ! tells whether MODULE is initialized
         INTEGER :: N              ! dimension of matrix currently treated
                                   ! in all routines checked against the actual dimension
         INTEGER :: MPI_COMM       ! mpi communicator used for scaLAPACK
         INTEGER :: NPROCS         ! number of PEs
         INTEGER :: ICTXT          ! context handle of grid (usually MPI context)
         INTEGER :: NP             ! number of rows on the processor
         INTEGER :: NQ             ! number of cols on the processor
         INTEGER :: NPROW,NPCOL    ! processor grid dimensions
         INTEGER :: MYROW,MYCOL    ! processor coordinates in ps. grid
         INTEGER :: LWWORK,LIWORK,LRWORK,MALLOCPQ ! sizes of scalapack workarrays
      END TYPE scalapack_des

      ! global sca-lapack descriptor
      ! used as "default" in VASP

      TYPE (scalapack_des), SAVE, PRIVATE :: GSD
      INTEGER,SAVE   :: DESCSTD( DLEN_ ) ! distributed matrix descriptor array

! few notes: many of the variables are redundant with DESCSTS
!
! GSD%NP= NUMROC(GSD%N,DESCSTD(MB_),GSD%MYROW,0,GSD%NPROW)
! GSD%NQ= NUMROC(GSD%N,DESCSTD(NB_),GSD%MYCOL,0,GSD%NPCOL)
! CALL BLACS_GRIDINFO(DESCSTD(CTXT_),GSD%NPROW,GSD%NPCOL,GSD%MYROW,GSD%MYCOL) to determine GSD%NPROW,...GSD%MYCOL
!
! not sure how lightweighted the evaluation is
!

! size of arrays defined with GDEF
#ifdef gammareal
      INTEGER, PARAMETER,PRIVATE :: MCOMP=1
#else
      INTEGER, PARAMETER,PRIVATE :: MCOMP=2
#endif

 CONTAINS
!=======================================================================
!
! this subroutine sets LscaAWARE to a default value if not read from
! INCAR
!
!=======================================================================
    SUBROUTINE  INIT_SCALAAWARE( NB_TOT, NRPLWV, COMM )
      IMPLICIT NONE

      INTEGER :: NB_TOT, NRPLWV
      TYPE (communic) COMM
    ! local
      REAL(q) :: STORAGE_HAM, STORAGE_WAVE

#ifdef gammareal
      STORAGE_HAM =8._q*NB_TOT*NB_TOT
      STORAGE_WAVE=16._q*NB_TOT*NRPLWV/ COMM%NCPU
#else
      STORAGE_HAM =16._q*NB_TOT*NB_TOT
      STORAGE_WAVE=16._q*NB_TOT*NRPLWV/ COMM%NCPU
#endif
      IF (.NOT. LscaAWARE_read ) THEN
         IF (STORAGE_HAM>STORAGE_WAVE) THEN
            LscaAWARE=LscaLAPACK
         ELSE
            LscaAWARE=.FALSE.
         ENDIF
      ENDIF

    END SUBROUTINE

!=======================================================================
!
! function to determine the blocking factor used in scaLAPACK
!
!=======================================================================

    SUBROUTINE  QUERRY_SCALA_NB( NB_QUERRY)
      NB_QUERRY=NB
    END SUBROUTINE QUERRY_SCALA_NB

!=======================================================================
! helper function, factorizes NPROCS into NPROW*NPCOL
! with restriction that NPROW ~ NPCOL and NPROW >= NPCOL
!=======================================================================
    SUBROUTINE  FERMAT_RAZOR( NPROCS, NPROW, NPCOL )
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE
       INTEGER, INTENT(IN)  :: NPROCS
       INTEGER, INTENT(OUT) :: NPROW, NPCOL
       ! local
       INTEGER              :: NPRIMES
       INTEGER, POINTER     :: PRIMES(:) => NULL()
       INTEGER, POINTER     :: ALL_PRIMES(:) => NULL()
       INTEGER, POINTER     :: NEXP(:) => NULL()
       INTEGER              :: I,IP
       INTEGER              :: NODD
       INTEGER              :: A, B
       PROFILING_START( 'FERMAT_RAZOR' )

       IF ( NPROCS == 1 ) THEN
          NPROW = 1
          NPCOL = 1
          PROFILING_STOP( 'FERMAT_RAZOR' )
       ENDIF

       ! obtain all primes up to SQRT( NPROCS  )
       CALL ERATOSTHENES_SIEVE( NPROCS, PRIMES, ALL_PRIMES )
       NPRIMES = SIZE( PRIMES )
       ! obtain prime factorization of NPROCS
       CALL PRIME_FACTORIZATION( NPROCS, NEXP, PRIMES )

       ! first check if NPROCS is a perfect square
       IF ( IS_PERFECT_SQUARE( NPROCS ) ) THEN
          NPROW = PERFECT_SQUARE( NPROCS )
          NPCOL = PERFECT_SQUARE( NPROCS )
       ELSE
       ! the idea is to factorize NPROCS = A*B
       ! with restriction |A - B| -> 0
       ! fermat sieve works well if N is odd
       NODD = NPROCS/(PRIMES(1)**NEXP(1))

       ! in case NPROCS is even
       IF ( NODD == 1 ) THEN

          NPROW = PRIMES(1)**( NEXP(1)/2 )
          NPCOL = PRIMES(1)**( NEXP(1)/2+MOD(NEXP(1),2) )
       ! in case NPROCS is not even
       ! use fermats sieve on NODD
       ELSE

          ! perfect square nothing to do
          IF ( IS_PERFECT_SQUARE( NODD ) ) THEN
             NPROW = PERFECT_SQUARE( NODD )
             NPCOL = PERFECT_SQUARE( NODD )
          ELSE
             IF ( IS_PRIME( NODD, ALL_PRIMES ) ) THEN
                NPROW = 1
                NPCOL = NODD
             ELSE
                ! Fermat's way to factorize an odd integer:
                ! NODD = A**2- B**2 = ( A - B ) * ( A + B )
                ! and write A**2 - N = B**2
                ! and increase A by 1 until A**2 - N is a perfect square
                ! this gives B and finally the factorization
                A = CEILING( SQRT( REAL(NODD) ) )
                B = A*A - NODD
                DO
                   IF ( IS_PERFECT_SQUARE( B ) ) EXIT
                   A = A+1
                   B = A*A-NODD
                ENDDO
                NPROW = A+PERFECT_SQUARE( B )
                NPCOL = A-PERFECT_SQUARE( B )
             ENDIF
          ENDIF
          ! distribute remaining factors of 2
          ! such that NPROW is roughly NPCOL
          DO I = 1, NEXP(1)
             IF ( NPROW < NPCOL ) THEN
                NPROW = NPROW*2
             ELSE
                NPCOL = NPCOL*2
             ENDIF
          ENDDO
       ENDIF

       ENDIF ! perfect square

       IF ( NPROW*NPCOL /= NPROCS ) THEN
          CALL vtutor%error("ERROR, FERMAT_RAZOR failed " // str(NPROCS) // " " // str(NPROW) // " "&
              // str(NPCOL))
       ENDIF

       DEALLOCATE( PRIMES )
       DEALLOCATE( ALL_PRIMES )
       DEALLOCATE( NEXP )
       NULLIFY( PRIMES )
       NULLIFY( ALL_PRIMES )
       NULLIFY( NEXP )

       PROFILING_STOP( 'FERMAT_RAZOR' )

       CONTAINS

       ! Eratosthenes_Sieve to obtain primes up to N
       SUBROUTINE ERATOSTHENES_SIEVE(N, PRIM_SQRT, PRIMES)
          IMPLICIT NONE
          INTEGER, INTENT(IN)  :: N
          INTEGER, POINTER     :: PRIM_SQRT(:)
          INTEGER              :: I, J, MAXPRIME, STAT
          LOGICAL              :: A(N)
          INTEGER, POINTER     :: PRIMES(:)

          MAXPRIME = FLOOR(SQRT(REAL(N)))
          A = .TRUE.

          DO I=2,MAXPRIME
             J = I*I
             DO
                A(J) = .FALSE.
                J = J + I ; IF ( J .GT. N ) EXIT
             ENDDO
          ENDDO !I

          ALLOCATE( PRIMES( COUNT(A)-1 ) )
          J = 1
          DO I=2,N ! SKIP 1
             IF ( .NOT. A(I) ) CYCLE
             PRIMES( J ) = I
             J = J + 1 ; IF ( J > SIZE(PRIMES) ) EXIT
          ENDDO


          J = 0
          DO I = 1, SIZE( PRIMES )
             IF ( PRIMES( I )*PRIMES( I ) > N ) THEN
                J = J+1
                EXIT
             ENDIF
             J = J+1
          ENDDO

          ALLOCATE( PRIM_SQRT( J ) )
          DO I = 1, J
             PRIM_SQRT( I ) = PRIMES( I )
          ENDDO
       END SUBROUTINE ERATOSTHENES_SIEVE

       ! obtains prime factoriztion of N
       SUBROUTINE PRIME_FACTORIZATION( N, NEXP, PRIMES )
          IMPLICIT NONE
          INTEGER          :: N
          INTEGER, POINTER :: NEXP(:)
          INTEGER          :: PRIMES(:)
          ! local
          INTEGER          :: NPRIMES
          INTEGER          :: IP
          INTEGER          :: M

          NPRIMES=SIZE( PRIMES )
          ALLOCATE( NEXP( NPRIMES ) )
          NEXP = 0
          ! 1 is not considered as prime here
          M = N
          !loop over all primes in PRIMES
          prime: DO IP = 1 , NPRIMES
             ! find exponent
             DO
                ! if this prime is a divisor
                IF( MOD( M, PRIMES( IP ) ) == 0 ) THEN
                   M = M/PRIMES( IP )
                   ! exponent of factor
                   NEXP( IP ) = NEXP( IP ) + 1
                ELSE
                   EXIT
                ENDIF
             ENDDO
          ENDDO prime
       END SUBROUTINE PRIME_FACTORIZATION

       ! returns true if B is prime
       FUNCTION IS_PRIME( B, PRIMES )
          LOGICAL :: IS_PRIME
          INTEGER :: PRIMES(:)
          INTEGER :: B
          ! local
          INTEGER :: I
          INTEGER :: ILIMIT

          ILIMIT = SIZE( PRIMES )

          DO I = 1, ILIMIT
             IF ( B == PRIMES( I ) ) EXIT
          ENDDO

          IF ( I > ILIMIT ) THEN
             IS_PRIME = .FALSE.
          ELSE
             IS_PRIME = .TRUE.
          ENDIF
       END FUNCTION IS_PRIME

       ! returns true if B is perfect square
       FUNCTION IS_PERFECT_SQUARE( B )
          LOGICAL :: IS_PERFECT_SQUARE
          INTEGER :: B
          ! local
          INTEGER :: I
          INTEGER :: ILIMIT

          ILIMIT =  CEILING( SQRT( REAL(B) ) )
          DO I = 1, ILIMIT
             IF ( I*I == B ) EXIT
          ENDDO
          IF ( I > ILIMIT ) THEN
             IS_PERFECT_SQUARE = .FALSE.
          ELSE
             IS_PERFECT_SQUARE = .TRUE.
          ENDIF
       END FUNCTION IS_PERFECT_SQUARE

       ! returns perfect square of B
       FUNCTION PERFECT_SQUARE( B )
          INTEGER :: PERFECT_SQUARE
          INTEGER :: B
          ! local
          INTEGER :: I
          INTEGER :: ILIMIT

          ILIMIT =  CEILING( SQRT( REAL(B) ) )
          DO I = 1, ILIMIT
             IF ( I*I == B ) EXIT
          ENDDO
          PERFECT_SQUARE = I
       END FUNCTION PERFECT_SQUARE

    END SUBROUTINE FERMAT_RAZOR

!=======================================================================
!
! messy scalapack initialization, calculate all required workspace
!
! the first routine
!    SUBROUTINE INIT_scala(myCOMM, N )
! sets up a global structure for scaLAPACK communication
! this is sufficiently flexible for most parts of VASP
!
! However in the more complicated cases, such as GW, it might
! be necessary to use different contextes for instance
! subdividing the nodes into different communication universes
! in that case
!   SUBROUTINE INIT_scala_DESC(myCOMM, N , DESCA, GS)
! should be used with local DESCA and GS (scalapack_des) descriptors
!
!=======================================================================

    SUBROUTINE INIT_scala(myCOMM, N )
      TYPE(communic) :: myCOMM
      INTEGER N

      CALL INIT_scala_DESC(myCOMM, N , DESCSTD, GSD)

    END SUBROUTINE INIT_scala

!
! here is the context safe version
!

    SUBROUTINE INIT_scala_DESC(myCOMM, N , DESCA, GS)

      USE main_mpi            ! to get world communicator
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      INTEGER              :: DESCA( DLEN_ ) ! distributed matrix descriptor array
      TYPE (scalapack_des) :: GS             ! descriptor that should be initialized by the routine
      TYPE(communic) :: myCOMM               ! VASP MPI communicator
      INTEGER N                              ! dimension of matrices to be handled
! local
      INTEGER INFO,NEIG
      INTEGER NP0,NQ0,NN
      INTEGER ierror

      INTRINSIC :: MAX,MIN
      INTEGER,EXTERNAL :: NUMROC,ICEIL
#ifdef oldversion
      INTEGER :: FACTORS(4), PWR(4)
      DATA FACTORS /2, 3, 5, 7/   ! a few prime factors likely to be be used
      INTEGER :: NNPCOL, NNPROW, MR, NDUM, NMIN, I
#endif
      EXTERNAL BLACS_PINFO
      EXTERNAL DESCINIT
#ifdef ELPA
      INTEGER ISUCCESS
#endif

      ! if gs%magic is set and the communicator changes or
      ! the size N of the matrix changes, then reset the GS structure
      IF (GS%MAGIC == 1789231) THEN
         IF (myCOMM%MPI_COMM /= GS%MPI_COMM .OR. GS%N /=N) THEN
            ! exit the BLACS context first
            CALL BLACS_GRIDEXIT(GS%ICTXT)
            ! require re-initialization
            GS%MAGIC=0
         ENDIF
      ENDIF

! if GS%MAGIC is not set to the right value reinitialize the GS structure
      IF (GS%MAGIC /= 1789231) THEN

      GS%MAGIC=1789231
      GS%N=N
      GS%MPI_COMM=myCOMM%MPI_COMM

! number of nodes (for one image)
      GS%NPROCS = myCOMM%NCPU

# ifdef oldversion
! determine processor grid, according to naive guess
      PWR=0
      MR = GS%NPROCS
      GS%NPROW=1; GS%NPCOL=1
      DO I=1,4
         DO
           IF( MOD(MR,FACTORS(I)) == 0 ) THEN
              PWR(I) = PWR(I) + 1
              MR = MR/FACTORS(I)
           ELSE
              EXIT
           ENDIF
         ENDDO
         NNPROW = PWR(I)/2
         NNPCOL = PWR(I) - NNPROW
         IF(MOD(I,2)==0) THEN         ! to make the grid more uniformly distributed (hopefully)
           NDUM = NNPCOL
           NNPCOL = NNPROW
           NNPROW = NDUM
         ENDIF
         GS%NPROW = GS%NPROW*FACTORS(I)**NNPROW
         GS%NPCOL = GS%NPCOL*FACTORS(I)**NNPCOL
      ENDDO
      IF( GS%NPROW*GS%NPCOL < GS%NPROCS )THEN  ! final adjustment to the total number of procs
         NMIN = MIN(GS%NPROW,GS%NPCOL)
         NDUM = GS%NPROCS/GS%NPROW/GS%NPCOL
         IF(GS%NPROW==NMIN)THEN          ! again, hoping that in this way the grid is more uniform
           GS%NPROW = GS%NPROW * NDUM   ! you might think to a better algorithm ( I dont think it
         ELSE                             ! really matters )
           GS%NPCOL = GS%NPCOL * NDUM
         ENDIF
      ENDIF
#else
! determine precessor grid using Fermat Sieve
      CALL FERMAT_RAZOR( GS%NPROCS, GS%NPROW, GS%NPCOL )
#endif

! make processor grids
      CALL MPI_barrier( myCOMM%MPI_COMM, ierror )

      CALL PROCMAP( myCOMM, GS%ICTXT, 2, GS%NPROW, GS%NPCOL)

! calculate local size of matrices

      CALL BLACS_GRIDINFO( GS%ICTXT, GS%NPROW, GS%NPCOL, GS%MYROW, GS%MYCOL )

      GS%NP = NUMROC(N,NB,GS%MYROW,0,GS%NPROW)   ! get number of rows on proc
      GS%NQ = NUMROC(N,NB,GS%MYCOL,0,GS%NPCOL)   ! get number of cols on proc

      CALL DESCINIT(DESCA,N,N,NB,NB,0,0,GS%ICTXT,MAX( 1, GS%NP ),INFO ) ! setup descriptor
      IF (INFO.NE.0) THEN
        CALL vtutor%bug("INIT_SCALA: DESCA, DESCINIT, INFO: " // str(INFO), __FILE__, __LINE__)
      ENDIF

#ifdef ELPA
! if one of the ranks has no data we can not use ELPA
      ISUCCESS=0; IF (GS%NP==0.OR.GS%NQ==0) ISUCCESS=1
      CALLMPI( M_sum_i(myCOMM,ISUCCESS,1) )
      IF (ISUCCESS>0) THEN
         LELPA=.FALSE.
         IF (myCOMM%NODE_ME==myCOMM%IONODE) &
            WRITE(*,*) 'ELPA: falling back to scaLAPACK'
      ELSE
         ELPA_MPI_COMM_ALL=myCOMM%MPI_COMM
      ENDIF
#endif

! calculate scalapack workspace and allocate, pDSSYEX_ZHEEVX only

   ! iwork
      GS%LIWORK=6*MAX(N,GS%NPROW*GS%NPCOL+1,4)
   ! wwork
      NN = MAX( N, NB, 2 )
      NEIG = N                ! number of eigenvectors requested
      IF ((DESCA(MB_).NE.NB).OR.(DESCA(NB_).NE.NB) .OR. &
          (DESCA(RSRC_).NE.0) .OR.(DESCA( CSRC_).NE.0)) THEN
        CALL vtutor%bug("INIT_SCALA: pSSYEZS_ETC, ERROR", __FILE__, __LINE__)
      ENDIF
      NP0 = NUMROC( NN, NB, 0, 0, GS%NPROW )
      NQ0 = MAX( NUMROC( NEIG, NB, 0, 0, GS%NPCOL ), NB )
#ifdef gammareal
      GS%LRWORK=1
      GS%LWWORK=5*N+MAX(5*NN,NP0*NQ0)+ICEIL(NEIG,GS%NPROW*GS%NPCOL)*NN+2*NB*NB+ &
         NCLUST*N
#else
      GS%LRWORK=4*N+MAX( 5*NN, NP0 * NQ0 )+ICEIL( NEIG, GS%NPROW*GS%NPCOL)*NN+ &
         NCLUST*N
      GS%LWWORK=N + MAX((NP0 +NQ0 +NB ) * NB, 3)
#endif

      GS%MALLOCPQ=MAX(GS%NP*GS%NQ,1)

! gK: dump global variables to see what has been set
!      WRITE(*,*) GS%NPROCS         ,' number of PEs'
!      WRITE(*,*) GS%ICTXT          ,' context handle of grid'
!      WRITE(*,*) GS%NPROW,GS%NPCOL ,' processor grid dimensions'
!      WRITE(*,*) GS%MYROW,GS%MYCOL ,' processor coordinates in ps. grid'
!      WRITE(*,*) GS%NP             ,' number of rows on the processor'
!      WRITE(*,*) GS%NQ             ,' number of cols on the processor'
!      WRITE(*,*) GS%LWWORK,GS%LIWORK,GS%LRWORK,GS%MALLOCPQ ,' sizes of scalapack workarrays'

      ENDIF

    END SUBROUTINE INIT_scala_DESC


!=======================================================================
!
! any routine should either call INIT_scala or CHECK_scala
! CHECK_scala, gracefully aborts if the GSD%N is incompatible
! to the supplied value N
!
!=======================================================================

    SUBROUTINE CHECK_scala (N, STRING)
      USE tutor, ONLY: vtutor
      INTEGER N
      CHARACTER (LEN=*) :: STRING
      IF (N /= GSD%N) THEN
         CALL vtutor%bug("CHECK_scala: BG_INIT_SCALA_BG was not called beforehand for " &
            // "this k-point \n CHECK_scala was called by " // STRING, __FILE__, __LINE__)
      ENDIF

    END SUBROUTINE CHECK_scala

!=======================================================================
!
! routine to querry the GSD%NP and GSD%NQ variables
!
!=======================================================================

    FUNCTION SCALA_NP()
      INTEGER SCALA_NP
      SCALA_NP=MAX(GSD%NP,1)
    END FUNCTION SCALA_NP


    FUNCTION SCALA_NQ()
      INTEGER SCALA_NQ
      SCALA_NQ=MAX(GSD%NQ,1)
    END FUNCTION SCALA_NQ

!=======================================================================
!
!     this subroutine maps processors onto a grid
!     and allows to perform scaLAPACK operations on a subset of nodes
!     taken from BLACS example code (written by Clint Whaley 7/26/94)
!     modified by gD and gK
!
!=======================================================================

    SUBROUTINE PROCMAP(VCOMM,CONTEXT, MAPPING, NPROW, NPCOL)
      USE main_mpi     ! to get MPI communication handle
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE
      TYPE(communic) :: VCOMM

!     .. Scalar Arguments ..
      INTEGER :: CONTEXT, MAPPING, NPROW, NPCOL
      INTEGER :: IMAP(NPROW, NPCOL)
      INTEGER ierror
!     ..
!
!  Purpose
!  =======
!  PROCMAP maps NPROW*NPCOL onto an VASP-Communicator
!  It is a nice system independent implementation
!
!  Arguments
!  =========
!
!  CONTEXT      (output) INTEGER
!               This integer is used by the BLACS to indicate a context.
!               A context is a universe where messages exist and do not
!               interact with other contexts messages.  The context includes
!               the definition of a grid, and each processs coordinates in it.
!
!  MAPPING      (input) INTEGER
!               Way to map processes to grid.  Choices are:
!               1 : row-major natural ordering
!               2 : column-major natural ordering
!
!  NPROW        (input) INTEGER
!               The number of process rows the created grid
!               should have.
!
!  NPCOL        (input) INTEGER
!               The number of process columns the created grid
!               should have.
!
!  =====================================================================
!
      INTEGER,EXTERNAL :: BLACS_PNUM
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDMAP
      INTEGER TMPCONTXT, NPROCS, I, J, K
      INTEGER BLACS_IDS(VCOMM%NCPU), MY_BLACS_ID
      INTEGER NPROW_TMP, NPCOL_TMP, MYROW_TMP, MYCOL_TMP
!
!     See how many processes there are in the system
!     NPROCS here refers to all processors being used

      NPROCS=VCOMM%NCPU
!
!     Temporarily map all processes into 1 x NPROCS grid
!
      TMPCONTXT=VCOMM%MPI_COMM
      CALL MPI_barrier( VCOMM%MPI_COMM, ierror )
      CALL BLACS_GRIDINIT( TMPCONTXT, 'Row', 1, NPROCS )
      CALL MPI_barrier( VCOMM%MPI_COMM, ierror )

!
!     If we want a row-major natural ordering
!
      CALL BLACS_GRIDINFO( TMPCONTXT, NPROW_TMP, NPCOL_TMP, MYROW_TMP, MYCOL_TMP)

      MY_BLACS_ID= BLACS_PNUM(TMPCONTXT, MYROW_TMP, MYCOL_TMP)
      CALL MPI_ALLGATHER(MY_BLACS_ID,1,MPI_INTEGER, &
     &     BLACS_IDS(1),1,MPI_INTEGER,VCOMM%MPI_COMM,IERROR)
      IF (IERROR.NE.0) CALL vtutor%error('ERROR allgather in procmap' // str(ierror))

      K=1
      IF (MAPPING .EQ. 1) THEN
         DO I = 1, NPROW
            DO J = 1, NPCOL
               IMAP(I, J) = BLACS_IDS( K )
               K = K + 1
            END DO
         END DO
!
!     If we want a column-major natural ordering
!
      ELSE IF (MAPPING .EQ. 2) THEN
         DO J = 1, NPCOL
            DO I = 1, NPROW
               IMAP(I, J) = BLACS_IDS( K )
               K = K + 1
            END DO
         END DO
      ELSE
         CALL vtutor%error("INIT_SCALA, PROCMAP: unknown mapping, STOP")
      END IF
!
!     Free temporary context
!
      CALL MPI_barrier( VCOMM%MPI_COMM, ierror )
      CALL BLACS_GRIDEXIT(TMPCONTXT)
!      WRITE (*,*) 'MAPPING',NPCOL,NPROW
!
!     Apply the new mapping to form desired context
!
      CONTEXT=VCOMM%MPI_COMM
      CALL BLACS_GRIDMAP( CONTEXT, IMAP, NPROW, NPROW, NPCOL )

!      WRITE(*,'("local map",16I3)') ((IMAP(I,J),I=1,NPROW),J=1,NPCOL)

      RETURN
    END SUBROUTINE PROCMAP

!=======================================================================
!
! the following routines rely on the global GSD communicator
! they however re-initialize GSD as well, so should be safe
! against changes of the matrix dimension
!
!=======================================================================
!
! this subroutine determines U^-1: A= U+ U and invert U
!   1. calls INIT_scala
!   2. calls the setup of distributed the distributed matrix
!   3. calls PDPOTRF and PDTRTRI (gamma) or PZPOTRF and PZTRTRI
!   4. calls reconstruction of distributed data into patched matrix
! note: the sum over processors of the patched matrix is not realised
!       in this module
!
!=======================================================================


    SUBROUTINE pPOTRF_TRTRI (COMM, AMATIN, NDIM, N)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER NDIM           ! leading dimension of AMATIN
      INTEGER N              ! NxN matrix to be distributed
      GDEF    AMATIN(NDIM,N) ! input/output matrix
! local variables
      INTEGER INFO

      GDEF, ALLOCATABLE ::  A(:)

      INTEGER,EXTERNAL :: NUMROC
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER,EXTERNAL :: BLACS_PNUM

      PROFILING_START('ppotrf_trtri')

      CALL INIT_scala(COMM, N)

!-----------------------------------------------------------------------
! allocate workarray (on T3D allocated in init_T3D)
!-----------------------------------------------------------------------
      ALLOCATE(A(GSD%MALLOCPQ))
!-----------------------------------------------------------------------
! do actual calculation
!-----------------------------------------------------------------------
   calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN
!       get parts of the matrix into the local arrays
        CALL DISTRI(AMATIN,NDIM,N,A,DESCSTD)

        INFO=0
#ifdef gammareal
        CALL PDPOTRF( 'U', N, A, 1, 1, DESCSTD, INFO )
#else
        CALL PZPOTRF( 'U', N, A, 1, 1, DESCSTD, INFO )
#endif

        IF (INFO.NE.0) THEN
          CALL vtutor%error("pPOTRF_TRTRI, POTRF, INFO: " // str(INFO) // "\n STOP")
        ENDIF
#ifdef gammareal
         CALL PDTRTRI &
#else
         CALL PZTRTRI &
#endif
     &    ('U', 'N', N, A, 1, 1, DESCSTD, INFO)
        IF (INFO.NE.0) THEN
          CALL vtutor%error("pPOTRF_TRTRI, TRTRI, INFO: " // str(INFO))
        ENDIF

      CALL RECON(AMATIN,NDIM,N,A,DESCSTD)
   ENDIF calc

!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      DEALLOCATE(A)

      PROFILING_STOP('ppotrf_trtri')

      RETURN
    END SUBROUTINE pPOTRF_TRTRI


!=======================================================================
!
! this subroutine determines A^-1: A= U+ U and invert (U+ U)
!   1. calls INIT_scala
!   2. calls the setup of distributed the distributed matrix
!   3. calls PDPOTRF and PDTRTRI (gamma) or PZPOTRF and PZTRTRI
!   4. calls reconstruction of distributed data into patched matrix
! note: the sum over processors of the patched matrix is not realised
!       in this module
!
!=======================================================================


    SUBROUTINE pPOTRF_POTRI (COMM, AMATIN, NDIM, N)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER NDIM           ! leading dimension of AMATIN
      INTEGER N              ! NxN matrix to be distributed
      GDEF    AMATIN(NDIM,N) ! input/output matrix
! local variables
      INTEGER INFO

      GDEF, ALLOCATABLE ::  A(:)

      INTEGER,EXTERNAL :: NUMROC
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER,EXTERNAL :: BLACS_PNUM

      PROFILING_START('ppotrf_potri')

      CALL INIT_scala(COMM, N)
!-----------------------------------------------------------------------
! allocate workarray
!-----------------------------------------------------------------------
      ALLOCATE(A(GSD%MALLOCPQ))
!-----------------------------------------------------------------------
! do actual calculation
!-----------------------------------------------------------------------
   calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN
!       get parts of the matrix into the local arrays
        CALL DISTRI(AMATIN,NDIM,N,A,DESCSTD)

        INFO=0
#ifdef gammareal
        CALL PDPOTRF( 'U', N, A, 1, 1, DESCSTD, INFO )
#else
        CALL PZPOTRF( 'U', N, A, 1, 1, DESCSTD, INFO )
#endif

        IF (INFO.NE.0) THEN
          CALL vtutor%error("pPOTRF_TRTRI, POTRF, INFO: " // str(INFO) // "\n STOP")
        ENDIF
#ifdef gammareal
         CALL PDPOTRI &
#else
         CALL PZPOTRI &
#endif
     &    ('U', N, A, 1, 1, DESCSTD, INFO)
        IF (INFO.NE.0) THEN
          CALL vtutor%error("pPOTRF_POTRI, POTRI, INFO: " // str(INFO))
        ENDIF

      CALL RECON(AMATIN,NDIM,N,A,DESCSTD)
   ENDIF calc
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      DEALLOCATE(A)

      PROFILING_STOP('ppotrf_potri')

      RETURN
    END SUBROUTINE pPOTRF_POTRI


!=======================================================================
!
! call to pZHEEVX respectively pDSYEVX
! i.e. diagonalization of a symmetric or hermitian matrix
!
!=======================================================================

    SUBROUTINE pDSSYEX_ZHEEVX(COMM, AMATIN, W, NDIM, N, ABSTOL_, AMATIN_SCALA)
#ifdef ELPA
      USE iso_c_binding
#endif
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N                          ! NxN matrix to be distributed
      INTEGER NDIM                       ! leading dimension of matrix
      GDEF    AMATIN(NDIM,N)             ! input/output matrix
      REAL(q) W(N)                       ! eigenvalues
      REAL(q), OPTIONAL      ::  ABSTOL_ ! tolerance
      GDEF, OPTIONAL, TARGET ::  AMATIN_SCALA(:) ! alternative input/output matrix (scaLAPACK format)
! local variables
      INTEGER ICLUSTR(2*COMM%NCPU)
      INTEGER IFAIL(N)
      REAL(q) :: GAP(COMM%NCPU)
      REAL(q) :: ZERO=0.0_q
      INTEGER INFO
      INTEGER M,NZ
      INTEGER I2
      REAL(q) :: ABSTOL

      GDEF, ALLOCATABLE, TARGET :: A(:) ! A matrix to be diagonalized (allocatable)
      GDEF, ALLOCATABLE, TARGET :: Z(:) ! Z eigenvector matrix
      GDEF, POINTER        ::     AP(:) ! A matrix to be diagonalized (pointer) either to A or AMATIN_SCALA
      REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
      GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
      INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array

      INTEGER,EXTERNAL :: NUMROC,ICEIL
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER NODE_ME,IONODE

#ifdef ELPA
      CLASS(elpa_t), POINTER :: ET
      INTEGER :: ISUCCESS

      GDEF, POINTER :: EP(:,:), EZ(:,:)
#endif

      PROFILING_START('pdssyex_zheevx')

      CALL INIT_scala(COMM, N)

      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE

      IF ((GSD%NPROW*GSD%NPCOL) > COMM%NCPU) THEN
         CALL vtutor%error("pDSSYEX_ZHEEVX: too many processors " // str(GSD%NPROW*GSD%NPCOL) // " "&
             // str(COMM%NCPU))
      ENDIF
!-----------------------------------------------------------------------
! allocation
!-----------------------------------------------------------------------
      IF (PRESENT(AMATIN_SCALA)) THEN
         IF (SIZE(AMATIN_SCALA)<GSD%MALLOCPQ) CALL vtutor%bug("pDSSYEX_ZHEEVX: AMATIN_SCALA too small: " &
               // str(SIZE(AMATIN_SCALA)) // " " // str(GSD%MALLOCPQ),__FILE__,__LINE__)
         AP=>AMATIN_SCALA
      ELSE
         ALLOCATE(A(GSD%MALLOCPQ))
         AP=>A
      ENDIF
      ALLOCATE(Z(GSD%MALLOCPQ))
      ALLOCATE(IWORK(GSD%LIWORK))
      ALLOCATE(WWORK(GSD%LWWORK))
      ALLOCATE(RWORK(GSD%LRWORK))

      IF (PRESENT(ABSTOL_)) THEN
         ABSTOL=ABSTOL_
      ELSE
         ABSTOL=ABSTOL_DEF
      ENDIF
!-----------------------------------------------------------------------
! do calculation
!-----------------------------------------------------------------------
      calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN

         !     get parts of the matrix into the local arrays
         IF (.NOT. PRESENT(AMATIN_SCALA)) THEN
            CALL DISTRI(AMATIN, NDIM, N, A, DESCSTD)
         ENDIF
         INFO=0

#ifdef ELPA
        ! call ELPA 2-stage solver
        IF (LELPA) THEN
           IF (ELPA_INIT(20171201) /= ELPA_OK) &
              CALL vtutor%error("ELPA_INIT: ELPA API version (20171201) not suppported by your ELPA library")

           CALL c_f_pointer(c_loc(AP),EP,[GSD%NP,GSD%NQ])
           CALL c_f_pointer(c_loc( Z),EZ,[GSD%NP,GSD%NQ])

           ET => ELPA_ALLOCATE()

           CALL ET%SET("na",  N, ISUCCESS)
           CALL ET%SET("nev", N, ISUCCESS)

           CALL ET%SET("local_nrows", GSD%NP, ISUCCESS)
           CALL ET%SET("local_ncols", GSD%NQ, ISUCCESS)

           CALL ET%SET("nblk", NB, ISUCCESS)

           CALL ET%SET("mpi_comm_parent", ELPA_MPI_COMM_ALL, ISUCCESS)
           CALL ET%SET("process_row", GSD%MYROW, ISUCCESS)
           CALL ET%SET("process_col", GSD%MYCOL, ISUCCESS)

           ISUCCESS = ET%SETUP()

           CALL ET%SET("solver", ELPA_SOLVER_2STAGE, ISUCCESS)
           CALL ET%EIGENVECTORS(EP, W, EZ, ISUCCESS)

           CALL ELPA_DEALLOCATE(ET)
           CALL ELPA_UNINIT()

           NULLIFY(EP,EZ)
        ENDIF
#endif
        ! call scaLAPACK routine
        IF (.NOT. LELPA) THEN
#ifdef gammareal
         CALL PDSYEVX( 'V', 'A', 'U', N, AP(1), 1, 1, DESCSTD, ZERO, ZERO, 13,  &
              -13, ABSTOL, M, NZ, W, ORFAC, Z, 1, 1, DESCSTD, WWORK, &
              GSD%LWWORK, IWORK, GSD%LIWORK, IFAIL, ICLUSTR, GAP, INFO )
#else
         CALL PZHEEVX( 'V', 'A', 'U', N, AP(1), 1, 1, DESCSTD, ZERO, ZERO, 13,  &
              -13, ABSTOL, M, NZ, W, ORFAC, Z, 1, 1, DESCSTD, &
              WWORK, GSD%LWWORK,  RWORK, GSD%LRWORK, IWORK, GSD%LIWORK, &
              IFAIL, ICLUSTR, GAP, INFO )
#endif
         IF (M.NE.N) THEN
            CALL vtutor%error("GSD%LWWORK " // str(GSD%LWWORK) // " " // str(GSD%LRWORK) // " " // &
               str(GSD%LIWORK) // " " // str(DESCSTD(M_)) // "\nERROR in subspace rotation PDSYEVX/ &
               &PZHEEVX: not enough eigenvalues found " // str(M) // " " // str(N))
         ENDIF
         IF (NZ.NE.N) THEN
            CALL vtutor%error("ERROR in subspace rotation PDSYEVX/ PZHEEVX: not enough eigenvectors &
               &computed " // str(M) // " " // str(N))
         ENDIF
         DO I2=1,N
            IF (IFAIL(I2).NE.0) THEN
               CALL vtutor%error("ERROR in subspace rotation PDSYEVX/ PZHEEVX: I2,IFAIL= " // str(I2) &
                  // " " // str(IFAIL(I2)))
            ENDIF
         ENDDO
         IF (INFO<0) THEN
            CALL vtutor%error("ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
         ELSE IF (INFO>0) THEN
            IF (MOD(INFO,2).NE.0) THEN
               CALL vtutor%error("ERROR eigenvector not converged PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
            ELSEIF (MOD(INFO/2,2).NE.0) THEN
               ! this error condition happens very often, but does not imply that
               ! the result is not usefull
               ! IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,'(1X,A,I5)') "WARNING in PDSYEVX/ PZHEEVX: eigenvectors not reorthogonalized IFAIL= ",INFO
            ELSE
               CALL vtutor%error("ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
            ENDIF
         ENDIF
        ENDIF

        IF (PRESENT(AMATIN_SCALA)) THEN
           AMATIN_SCALA(1:GSD%MALLOCPQ)=Z(1:GSD%MALLOCPQ)
        ELSE
           CALL RECON(AMATIN, NDIM, N, Z, DESCSTD)
        ENDIF

     ENDIF calc
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      IF (.NOT.PRESENT(AMATIN_SCALA)) THEN
         DEALLOCATE(A)
      ENDIF
      DEALLOCATE(Z)
      DEALLOCATE(IWORK)
      DEALLOCATE(WWORK)
      DEALLOCATE(RWORK)

      PROFILING_STOP('pdssyex_zheevx')

      RETURN
    END SUBROUTINE pDSSYEX_ZHEEVX

!=======================================================================
!
! call to pZHEEVD respectively pDSYEVD
! i.e. diagonalization of a symmetric or hermitian matrix
!
!=======================================================================

    SUBROUTINE pDSYEV_ZHEEVD(COMM, AMATIN, W, NDIM, N, AMATIN_SCALA)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N                          ! NxN matrix to be distributed
      INTEGER NDIM                       ! leading dimension of matrix
      GDEF    AMATIN(NDIM,N)             ! input/output matrix
      REAL(q) W(N)                       ! eigenvalues
      GDEF, OPTIONAL, TARGET ::  AMATIN_SCALA(:) ! alternative input/output matrix (scaLAPACK format)
! local variables
      GDEF, ALLOCATABLE,TARGET ::  A(:) ! A matrix to be diagonalized (allocatable)
      GDEF, POINTER  ::  AP(:)          ! A matrix to be diagonalized (pointer) either to A or AMATIN_SCALA
      INTEGER NODE_ME,IONODE

      CALL INIT_scala(COMM, N)
      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE

      IF ((GSD%NPROW*GSD%NPCOL) > COMM%NCPU) THEN
         CALL vtutor%error("pDSYEV_ZHEEVD: too many processors " // str(GSD%NPROW*GSD%NPCOL) // " "&
             // str(COMM%NCPU))
      ENDIF
!-----------------------------------------------------------------------
! allocation
!-----------------------------------------------------------------------
      IF (PRESENT(AMATIN_SCALA)) THEN
         AP=>AMATIN_SCALA
      ELSE
         ALLOCATE(A(GSD%MALLOCPQ))
         AP=>A
      ENDIF
!-----------------------------------------------------------------------
! do calculation
!-----------------------------------------------------------------------
      calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN

         !     get parts of the matrix into the local arrays (AP=>A)
         IF (.NOT. PRESENT(AMATIN_SCALA)) THEN
            CALL DISTRI(AMATIN, NDIM, N, A(1), DESCSTD)
         ENDIF
        
         CALL PDSYEV_ZHEEVD_DESC(AP(:), W, N, DESCSTD, COMM)
! if desired this can be also replaced by the standard matrix diagonalization
!         CALL PDSYEV_ZHEEV_DESC(AP(:), W, N, DESCSTD, COMM)

         IF (PRESENT(AMATIN_SCALA)) THEN
            AMATIN_SCALA=AP
         ELSE
            CALL RECON(AMATIN, NDIM, N, A(1), DESCSTD)
         ENDIF

     ENDIF calc
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      IF (.NOT.PRESENT(AMATIN_SCALA)) THEN
         DEALLOCATE(A)
      ENDIF

      RETURN
    END SUBROUTINE
    
!=======================================================================
!
! clean the matrix elements along the diagonal (force them to be real)
! as it is done in the seriel version
!
!=======================================================================


    SUBROUTINE BG_CHANGE_DIAGONALE(N, A, IU0)
      IMPLICIT NONE
      GDEF    A(*) ! distributed matrix
      INTEGER N    ! dimension of the matrix
      INTEGER IU0  ! IO unit for error report
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO
#ifndef gammareal
      CALL CHECK_scala (N, 'BG_CHANGE_DIAGONALE' )

      PROFILING_START('bg_change_diagonale')

      CALL BLACS_GRIDINFO(DESCSTD(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCSTD(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCSTD(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCSTD(NB_)
        J1RES=MIN(DESCSTD(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCSTD(MB_)
            I1RES=MIN(DESCSTD(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( (DESCSTD(MB_)*MYROW+NPROW*(I1-1)+I2) == (DESCSTD(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                 IF (ABS(AIMAG(A(IROW_JCOL)))>1E-2_q .AND. IU0>=0) THEN
                    WRITE(IU0,*)'WARNING: Sub-Space-Matrix is not hermitian subr', &
                    &                 AIMAG(A(IROW_JCOL)),DESCSTD(NB_)*MYCOL+NPCOL*(J1-1)+J2
                 ENDIF
                A(IROW_JCOL) = REAL( A(IROW_JCOL) ,KIND=q)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCSTD(LLD_)
        ENDDO
      ENDDO

      PROFILING_STOP('bg_change_diagonale')
#endif
      RETURN
    END SUBROUTINE BG_CHANGE_DIAGONALE

!=======================================================================
!
! matrix diagonalization scaLAPACK aware version
! this version is virtually identical to pDSSYEX_ZHEEVX
! when called with the additional AMATIN_SCALA argument
! i.e.
!  BG_pDSSYEX_ZHEEVX(COMM, A, W, N) could be replaced by
!  pDSSYEX_ZHEEVX(COMM, ADUMMY, W, NDUMMY, N, AMATIN_SCALA=A)
! A is a two dimension array here, however
!
!=======================================================================

    SUBROUTINE BG_pDSSYEX_ZHEEVX(COMM, A, W, N)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N               ! NxN matrix to be distributed
      GDEF    A(MAX(GSD%NP,1),MAX(GSD%NQ,1))  ! input/output matrix
      REAL(q) W(N)            ! eigenvalues
! local variables
      INTEGER ICLUSTR(2*COMM%NCPU)
      INTEGER IFAIL(N)
      REAL(q) :: GAP(COMM%NCPU)
      REAL(q) :: ZERO=0.0_q
      INTEGER INFO
      INTEGER M,NZ
      INTEGER I2

      GDEF    Z(MAX(GSD%NP,1),MAX(GSD%NQ,1))
      REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
      GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
      INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
      INTEGER,EXTERNAL :: NUMROC,ICEIL
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDINFO,  &
               BLACS_GRIDEXIT, BLACS_EXIT
      INTEGER NODE_ME,IONODE
#ifdef ELPA
      CLASS(elpa_t), POINTER :: ET
      INTEGER :: ISUCCESS
#endif

      PROFILING_START('bg_pdssyex_zheevx')

      CALL CHECK_scala (N, 'BG_pDSSYEX_ZHEEVX' )

      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE

      IF ((GSD%NPROW*GSD%NPCOL) > COMM%NCPU) THEN
        CALL vtutor%error("pDSSYEX_ZHEEVX: too many processors " // str(GSD%NPROW*GSD%NPCOL) // " " &
           // str(COMM%NCPU))
      ENDIF
!-----------------------------------------------------------------------
! allocation (on T3D allocated in init_T3D)
!-----------------------------------------------------------------------
      ALLOCATE(IWORK(GSD%LIWORK))
      ALLOCATE(WWORK(GSD%LWWORK))
      ALLOCATE(RWORK(GSD%LRWORK))
!-----------------------------------------------------------------------
! do calculation
!-----------------------------------------------------------------------
   calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN

      INFO=0

#ifdef ELPA
      ! call ELPA 2-stage solver
      IF (ELPA) THEN
           IF (ELPA_INIT(20171201) /= ELPA_OK) &
              CALL vtutor%error("ELPA_INIT: ELPA API version (20171201) not suppported by your ELPA library")

           ET => ELPA_ALLOCATE()

           CALL ET%SET("na",  N, ISUCCESS)
           CALL ET%SET("nev", N, ISUCCESS)

           CALL ET%SET("local_nrows", GSD%NP, ISUCCESS)
           CALL ET%SET("local_ncols", GSD%NQ, ISUCCESS)

           CALL ET%SET("nblk", NB, ISUCCESS)

           CALL ET%SET("mpi_comm_parent", ELPA_MPI_COMM_ALL, ISUCCESS)
           CALL ET%SET("process_row", GSD%MYROW, ISUCCESS)
           CALL ET%SET("process_col", GSD%MYCOL, ISUCCESS)

           ISUCCESS = ET%SETUP()

           CALL ET%SET("solver", ELPA_SOLVER_2STAGE, ISUCCESS)
           CALL ET%EIGENVECTORS(A, W, Z, ISUCCESS)

           CALL ELPA_DEALLOCATE(ET)
           CALL ELPA_UNINIT()
      ENDIF
#endif
      ! call scaLAPACK routine
      IF (.NOT. LELPA) THEN
#ifdef gammareal
        CALL PDSYEVX( 'V', 'A', 'U', N, A, 1, 1, DESCSTD, ZERO, ZERO, 13,  &
                    -13, ABSTOL_DEF, M, NZ, W, ORFAC, Z, 1, 1, DESCSTD, WWORK, &
                    GSD%LWWORK, IWORK, GSD%LIWORK, IFAIL, ICLUSTR, GAP, INFO )
#else
        CALL PZHEEVX( 'V', 'A', 'U', N, A, 1, 1, DESCSTD, ZERO, ZERO, 13,  &
                    -13, ABSTOL_DEF, M, NZ, W, ORFAC, Z, 1, 1, DESCSTD, &
                    WWORK, GSD%LWWORK,  RWORK, GSD%LRWORK, IWORK, GSD%LIWORK, &
                    IFAIL, ICLUSTR, GAP, INFO )
#endif
        IF (M.NE.N) THEN
          CALL vtutor%error("GSD%LWWORK " // str(GSD%LWWORK) // " " // str(GSD%LRWORK) // " " // &
             str(GSD%LIWORK) // " " // str(DESCSTD(M_)) // "\nERROR in subspace rotation PDSYEVX: not &
             &enough eigenvalues found " // str(M) // " " // str(N))
        ENDIF
        IF (NZ.NE.N) THEN
          CALL vtutor%error("ERROR in subspace rotation PDSYEVX: not enough eigenvectors computed " &
             // str(M) // " " // str(N))
        ENDIF
        DO I2=1,N
          IF (IFAIL(I2).NE.0) THEN
            CALL vtutor%error("ERROR in subspace rotation PDSYEVX: I2,IFAIL= " // str(I2) // " " // &
               str(IFAIL(I2)))
          ENDIF
        ENDDO
        IF (INFO<0) THEN
           CALL vtutor%error("ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
        ELSE IF (INFO>0) THEN
           IF (MOD(INFO,2).NE.0) THEN
              CALL vtutor%error("ERROR eigenvector not converged PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
           ELSEIF (MOD(INFO/2,2).NE.0) THEN
              ! this error condition happens very often, but does not imply that
              ! the result is not usefull
              ! IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,'(1X,A,I5)') "WARNING in PDSYEVX/ PZHEEVX: eigenvectors not reorthogonalized IFAIL= ",INFO
           ELSE
              CALL vtutor%error("ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
           ENDIF
        ENDIF
      ENDIF

      A=Z

   ENDIF calc
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      DEALLOCATE(IWORK)
      DEALLOCATE(WWORK)
      DEALLOCATE(RWORK)

      PROFILING_STOP('bg_pdssyex_zheevx')

      RETURN
    END SUBROUTINE BG_pDSSYEX_ZHEEVX

!=======================================================================
!
! matrix diagonalization scaLAPACK aware version
! this version is virtually identical to pDSSYED_ZHEEVD
! when called with the additional AMATIN_SCALA argument
! i.e.
!  BG_pDSSYED_ZHEEVD(COMM, A, W, N) could be replaced by
!  pDSSYED_ZHEEVD(COMM, ADUMMY, W, NDUMMY, N, AMATIN_SCALA=A)
! A is a two dimension array here, however
!
!=======================================================================

    SUBROUTINE BG_pDSYEV_ZHEEVD(COMM, A, W, N)
      USE string, ONLY: str
      USE tutor, ONLY: vtutor
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N               ! NxN matrix to be distributed
      GDEF    A(MAX(GSD%NP,1),MAX(GSD%NQ,1))  ! input/output matrix
      REAL(q) W(N)            ! eigenvalues
! local variables
      INTEGER NODE_ME,IONODE

      CALL CHECK_scala (N, 'BG_pDSYEV_ZHEEVD' )

      NODE_ME=COMM%NODE_ME
      IONODE =COMM%IONODE

      IF ((GSD%NPROW*GSD%NPCOL) > COMM%NCPU) THEN
        CALL vtutor%error("pDSYEV_ZHEEVD: too many processors " // str(GSD%NPROW*GSD%NPCOL) // " " &
           // str(COMM%NCPU))
      ENDIF
!-----------------------------------------------------------------------
! do calculation
!-----------------------------------------------------------------------
      calc: IF( GSD%MYROW < GSD%NPROW .AND. GSD%MYCOL < GSD%NPCOL ) THEN

        CALL PDSYEV_ZHEEVD_DESC2(A, W, N, DESCSTD, COMM)
      ENDIF calc
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
      RETURN
    END SUBROUTINE
    
!=======================================================================
!
! matrix LU decompositions scaLAPACK aware version
!  Determine A = U+ U (POTRF) and invert U (TRTRI)
!
!=======================================================================

    SUBROUTINE BG_pPOTRF_TRTRI(A, N, INFO)
      IMPLICIT NONE

      INTEGER N                  ! NxN matrix to be distributed
      GDEF    A(MAX(GSD%NP,1),MAX(GSD%NQ,1))     ! input/output matrix
      INTEGER INFO

      PROFILING_START('bg_ppotrf_trtri')

      INFO=0

      CALL CHECK_scala (N, 'BG_pPOTRF_TRTRI' )

#ifdef gammareal
      CALL PDPOTRF &
#else
      CALL PZPOTRF &
#endif
        & ('U',N, A(1,1),1,1,DESCSTD,INFO)
      IF (INFO/=0) RETURN
#ifdef gammareal
      CALL PDTRTRI &
#else
      CALL PZTRTRI &
#endif
          & ('U','N',N,A(1,1),1,1,DESCSTD,INFO)

      PROFILING_STOP('bg_ppotrf_trtri')

      END SUBROUTINE BG_pPOTRF_TRTRI

!=======================================================================
!
! matrix LU decompositions scaLAPACK aware version
!  Determine A = U+ U (POTRF) and invert A = U+U (POTRI)
!
!=======================================================================

    SUBROUTINE BG_pPOTRF_POTRI(A, N, INFO)
      IMPLICIT NONE

      INTEGER N               ! NxN matrix to be distributed
      GDEF    A(MAX(GSD%NP,1),MAX(GSD%NQ,1))     ! input/output matrix
      INTEGER INFO

      PROFILING_START('bg_ppotrf_potri')

      CALL CHECK_scala (N, 'BG_pPOTRF_POTRI' )
#ifdef gammareal
         CALL PDPOTRF &
#else
         CALL PZPOTRF &
#endif
          & ('U',N, A(1,1),1,1,DESCSTD,INFO)
#ifdef gammareal
         CALL PDPOTRI &
#else
         CALL PZPOTRI &
#endif
          & ('U',N, A(1,1), 1, 1, DESCSTD, INFO)

      PROFILING_STOP('bg_ppotrf_potri')

     END SUBROUTINE BG_pPOTRF_POTRI



!=======================================================================
!
! the following routines are called with descriptors DESCA
! the routines can be called from any context or subgroup
! of MPI processed as long as DESCA is properly set up
! they also to not refer to GSD
! in that sense they are much safer to use, but matter
! of fact the caller must "know" scalapack
!
! for the time being the check however that the passed matrix
! dimensions are consistent with those set up in GSD
! by calling CHECK_scala
! this should be removed after more extensive testing
!
!=======================================================================

    SUBROUTINE DISTRI(AMATIN, NDIM, N, A, DESCA)
      USE tutor, ONLY: vtutor
      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      INTEGER NDIM,N
      GDEF    AMATIN(NDIM,N),A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      IF (N>NDIM) THEN
         CALL vtutor%bug("DISTRI: leading dimension of matrix too small", __FILE__, __LINE__)
      END IF

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     ITEST=0

!     setup distributed matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              A(IROW+(JCOL-1)*DESCA(LLD_))= &
               AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                      DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DWRITE "DISTRI",ITEST,MYCOL,MYROW

      RETURN
    END SUBROUTINE DISTRI

    ! same as above using a single precision input array

    SUBROUTINE DISTRI_SINGLE(AMATIN, NDIM, N, A,DESCA)
      USE tutor, ONLY: vtutor
      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      INTEGER NDIM,N
      GDEFS   AMATIN(NDIM,N)
      GDEF    A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      IF (N>NDIM) THEN
         CALL vtutor%bug("internal error in scala.F: leading dimension of matrix too small", __FILE__, __LINE__)
      ENDIF

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     ITEST=0

!     setup distributed matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              A(IROW+(JCOL-1)*DESCA(LLD_))= &
               AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                      DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DWRITE "DISTRI",ITEST,MYCOL,MYROW

      RETURN
    END SUBROUTINE DISTRI_SINGLE

    ! same as above using a single precision input array

    SUBROUTINE DISTRI_SINGLE_SINGLE(AMATIN, NDIM, N, A,DESCA)
      USE tutor, ONLY: vtutor
      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      INTEGER NDIM,N
      GDEFS   AMATIN(NDIM,N)
      GDEFS   A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      IF (N>NDIM) THEN
         CALL vtutor%bug("internal error in scala.F: leading dimension of matrix too small", __FILE__, __LINE__)
      ENDIF

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     ITEST=0

!     setup distributed matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              A(IROW+(JCOL-1)*DESCA(LLD_))= &
               AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                      DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DWRITE "DISTRI",ITEST,MYCOL,MYROW

      RETURN
    END SUBROUTINE DISTRI_SINGLE_SINGLE


!=======================================================================
!
! merge distributed matrix into one large matrix
!
!=======================================================================

    SUBROUTINE RECON(AMATIN,NDIM, N,A,DESCA)

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      INTEGER NDIM,N
      GDEF    AMATIN(NDIM, N),A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

!     clear the matrix
      AMATIN=0._q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                     DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=   &
              A(IROW+(JCOL-1)*DESCA(LLD_))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON

    ! same as above using a single precision output array

    SUBROUTINE RECON_SINGLE(AMATIN,NDIM, N,A,DESCA)

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      INTEGER NDIM,N
      GDEFS   AMATIN(NDIM, N)
      GDEF    A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

!     clear the matrix
      AMATIN=0.0_q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                     DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=   &
              A(IROW+(JCOL-1)*DESCA(LLD_))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON_SINGLE


    ! same as above using a single precision output array

    SUBROUTINE RECON_SINGLE_SINGLE(AMATIN,NDIM, N,A,DESCA)

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      INTEGER NDIM,N
      GDEFS   AMATIN(NDIM, N)
      GDEFS   A(*)
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

!     clear the matrix
      AMATIN=0.0_q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          IROW=0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW=IROW+1
              AMATIN(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2,    &
                     DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=   &
              A(IROW+(JCOL-1)*DESCA(LLD_))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON_SINGLE_SINGLE

!=======================================================================
!
! matrix diagonalization scaLAPACK aware version
! this version expects a proper scaLAPACK descriptor to
! note that this routine does stop when eigenvectors canot be
! determined
!
!=======================================================================

     SUBROUTINE PDSYEVX_ZHEEVX_DESC(A, W, N, DESCA_, COMM, JOBIN, ABSTOL_)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               ! dimension of matrix
       GDEF    A(:)            ! input/output matrix
       REAL(q) W(N)            ! eigenvalues
       INTEGER DESCA_( DLEN_ ) ! distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
       REAL(q), OPTIONAL :: ABSTOL_
! local variables
       REAL(q) :: ABSTOL
       INTEGER ICLUSTR(2*COMM%NCPU)
       INTEGER IFAIL(N)
       REAL(q) :: GAP(COMM%NCPU)
       REAL(q) :: ZERO=0.0_q
       INTEGER INFO
       INTEGER M,NZ
       INTEGER I2

       GDEF, ALLOCATABLE ::     Z(:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors

       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1), IWORK(1))
DOESI  RWORK = 0.0_q; WWORK = 0; IWORK = 0
       ALLOCATE(Z(SIZE(A,KIND=qi8)))
DOESI  Z = 0

      IF (PRESENT(ABSTOL_)) THEN
         ABSTOL=ABSTOL_
      ELSE
         ABSTOL=ABSTOL_DEF
      ENDIF

#ifdef gammareal
       CALL PDSYEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
            -13, ABSTOL, M, NZ, W, ORFAC, Z(1), 1, 1, DESCA, WWORK, &
            LWWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )
       RWORK(1)=1
       WWORK(1)=WWORK(1)+NCLUST*N
#else
       CALL PZHEEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
            -13, ABSTOL, M, NZ, W, ORFAC, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, &
            IFAIL, ICLUSTR, GAP, INFO )
       RWORK(1)=RWORK(1)+NCLUST*N
#endif
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q)))
       LIWORK= MAX(1,IWORK(1))
       LRWORK= MAX(1,NINT(RWORK(1)))

       DEALLOCATE(RWORK, WWORK, IWORK)
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))
DOESI  RWORK = 0.0_q; WWORK = 0; IWORK = 0

#ifdef gammareal
       CALL PDSYEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
            -13, ABSTOL, M, NZ, W, ORFAC, Z(1), 1, 1, DESCA, WWORK, &
            LWWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )
#else
       CALL PZHEEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
            -13, ABSTOL, M, NZ, W, ORFAC, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, &
            IFAIL, ICLUSTR, GAP, INFO )
#endif
       IF (M.NE.N) THEN
          CALL vtutor%error("LWWORK " // str(LWWORK) // " " // str(LRWORK) // " " // str(LIWORK) // &
             " " // str(DESCA(M_)) // "\nERROR in subspace rotation PDSYEVX/ PZHEEVX: not enough &
             &eigenvalues found " // str(M) // " " // str(N))
       ENDIF
       IF (NZ.NE.N) THEN
          CALL vtutor%error("ERROR in subspace rotation PDSYEVX/ PZHEEVX: not enough eigenvectors &
             &computed " // str(M) // " " // str(N))
       ENDIF
       DO I2=1,N
          IF (IFAIL(I2).NE.0) THEN
             CALL vtutor%error("ERROR in subspace rotation PDSYEVX/ PZHEEVX: I2,IFAIL= " // str(I2) &
                // " " // str(IFAIL(I2)))
          ENDIF
       ENDDO
       IF (INFO<0) THEN
          CALL vtutor%error("ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
       ELSE IF (INFO>0) THEN
          IF (MOD(INFO,2).NE.0) THEN
             CALL vtutor%error("ERROR eigenvector not converged PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
          ELSEIF (MOD(INFO/2,2).NE.0) THEN
          ! this error condition happens very often, but does not imply that
          ! the result is not usefull
          !  IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,'(1X,A,I5)') "WARNING in PDSYEVX/ PZHEEVX: eigenvectors not reorthogonalized IFAIL= ",INFO
          ELSE
             CALL vtutor%error("ERROR in diagonalization PDSYEVX/ PZHEEVX: IFAIL= " // str(INFO))
          ENDIF
       ENDIF

       A=Z
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

    END SUBROUTINE PDSYEVX_ZHEEVX_DESC

!=======================================================================
!
! matrix diagonalization scaLAPACK aware version
! this version expects a proper scaLAPACK descriptor to
! note that this routine does stop when eigenvectors
!
!=======================================================================

     SUBROUTINE PSSYEVX_CHEEVX_DESC(A, W, N, DESCA_, COMM, JOBIN, ABSTOL_)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               ! dimension of matrix
       GDEFS    A(:)           ! input/output matrix
       REAL(qs) W(N)           ! eigenvalues
       INTEGER DESCA_( DLEN_ ) ! distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
       REAL(q), OPTIONAL :: ABSTOL_
! local variables
       REAL(qs) :: ABSTOL
       INTEGER ICLUSTR(2*COMM%NCPU)
       INTEGER IFAIL(N)
       REAL(qs) :: GAP(COMM%NCPU)
       REAL(qs) :: ZERO=0.0_q
       INTEGER INFO
       INTEGER M,NZ
       INTEGER I2

       GDEFS, ALLOCATABLE ::    Z(:)
       REAL(qs), ALLOCATABLE::  RWORK(:) ! work array
       GDEFS, ALLOCATABLE   ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors

       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1), IWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))

      IF (PRESENT(ABSTOL_)) THEN
         ABSTOL=ABSTOL_
      ELSE
         ABSTOL=ABSTOL_DEF
      ENDIF

#ifdef gammareal
       CALL PSSYEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
            -13, ABSTOL, M, NZ, W, ORFAC_single, Z(1), 1, 1, DESCA, WWORK, &
            LWWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )
       RWORK(1)=1
       WWORK(1)=WWORK(1)+NCLUST*N
#else
       CALL PCHEEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
            -13, ABSTOL, M, NZ, W, ORFAC_single, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, &
            IFAIL, ICLUSTR, GAP, INFO )
       RWORK(1)=RWORK(1)+NCLUST*N
#endif
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=qs)))
       LIWORK= MAX(1,IWORK(1))
       LRWORK= MAX(1,NINT(RWORK(1)))

       DEALLOCATE(RWORK, WWORK, IWORK)
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))

#ifdef gammareal
       CALL PSSYEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
            -13, ABSTOL, M, NZ, W, ORFAC_single, Z(1), 1, 1, DESCA, WWORK, &
            LWWORK, IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )
#else
       CALL PCHEEVX( JOB, 'A', 'U', N, A(1), 1, 1, DESCA, ZERO, ZERO, 13,  &
            -13, ABSTOL, M, NZ, W, ORFAC_single, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, &
            IFAIL, ICLUSTR, GAP, INFO )
#endif
       IF (M.NE.N) THEN
          CALL vtutor%error("LWWORK " // str(LWWORK) // " " // str(LRWORK) // " " // str(LIWORK) // &
             " " // str(DESCA(M_)) // "\nERROR in subspace rotation PSSYEVX/ PCHEEVX: not enough &
             &eigenvalues found " // str(M) // " " // str(N))
       ENDIF
       IF (NZ.NE.N) THEN
          CALL vtutor%error("ERROR in subspace rotation PSSYEVX/ PCHEEVX: not enough eigenvectors &
             &computed " // str(M) // " " // str(N))
       ENDIF
       DO I2=1,N
          IF (IFAIL(I2).NE.0) THEN
             CALL vtutor%error("ERROR in subspace rotation PSSYEVX/ PCHEEVX: I2,IFAIL= " // str(I2) &
                // " " // str(IFAIL(I2)))
          ENDIF
       ENDDO
       IF (INFO<0) THEN
          CALL vtutor%error("ERROR in diagonalization PSSYEVX/ PCHEEVX: IFAIL= " // str(INFO))
       ELSE IF (INFO>0) THEN
          IF (MOD(INFO,2).NE.0) THEN
             CALL vtutor%error("ERROR eigenvector not converged PSSYEVX/ PCHEEVX: IFAIL= " // str(INFO))
          ELSEIF (MOD(INFO/2,2).NE.0) THEN
          ! this error condition happens very often, but does not imply that
          ! the result is not usefull
!             IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,'(1X,A,I5)') "WARNING in PSSYEVX/ PCHEEVX: eigenvectors not reorthogonalized IFAIL= ",INFO
          ELSE
             CALL vtutor%error("ERROR in diagonalization PSSYEVX/ PCHEEVX: IFAIL= " // str(INFO))
          ENDIF
       ENDIF

       A=Z
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

     END SUBROUTINE PSSYEVX_CHEEVX_DESC


!=======================================================================
!
! matrix diagonalization using divide and conquer scaLAPACK
! this version expects a proper scaLAPACK descriptor to
!
!=======================================================================

     SUBROUTINE PDSYEV_ZHEEVD_DESC(A, W, N, DESCA_, COMM, JOBIN)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               ! dimension of matrix
       GDEF    A(:)            ! input/output matrix
       REAL(q) W(N)            ! eigenvalues
       INTEGER DESCA_( DLEN_ )  ! distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO

       GDEF, ALLOCATABLE ::     Z(:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       PROFILING_START( 'PDSYEV_ZHEEVD_DESC' )
       IF ( SIZE( A , KIND=qi8 ) == 0_qi8 ) THEN
          CALL vtutor%error("PDSYEV_ZHEEVD_DESC reports: Rank " // str(COMM%NODE_ME) // " has no data")
       ENDIF

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1), IWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))

#ifdef gammareal
       CALL PDSYEVD( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, IWORK, LIWORK, INFO )
       RWORK(1)=1
       WWORK(1)=WWORK(1)*2 ! often crashes without this mod
#else
       CALL PZHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
#endif
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q)))
       LIWORK= MAX(1,IWORK(1))
       LRWORK= MAX(1,NINT(RWORK(1)))

       DEALLOCATE(RWORK, WWORK, IWORK)
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))

#ifdef gammareal
       CALL PDSYEVD( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, IWORK, LIWORK, INFO )
#else
       CALL PZHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
#endif

       A=Z
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

       IF (INFO<0) THEN
          CALL vtutor%error("ERROR in diagonalization PDSYEVD/ PZHEEVD: IFAIL= " // str(INFO))
       ELSE IF (INFO>0) THEN
          IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,*) "WARNING: PDSYEVD/ PZHEEVD failed using fallback PDSYEV/ PZHEEV"
          CALL PDSYEV_ZHEEV_DESC(A, W, N, DESCA, COMM)
       ENDIF
! just in case, PDSYEVD does seem to be buggy and sometimes not all nodes
! have the right eigenvalues in particular if matrices are small and stored
! on a single core
       CALLMPI ( M_bcast_d_from(COMM , W(1), SIZE(W), 1))

       PROFILING_STOP( 'PDSYEV_ZHEEVD_DESC' )
    END SUBROUTINE PDSYEV_ZHEEVD_DESC

!=======================================================================
!
! matrix diagonalization using divide and conquer scaLAPACK
! same as above, but complex version regardless of gammareal
!
!=======================================================================

     SUBROUTINE PZHEEVD_DESC(A, W, N, DESCA_, COMM, JOBIN)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               ! dimension of matrix
       COMPLEX(q) A(:)            ! input/output matrix
       REAL(q) W(N)            ! eigenvalues
       INTEGER DESCA_( DLEN_ )  ! distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO

       COMPLEX(q), ALLOCATABLE ::     Z(:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       COMPLEX(q), ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       PROFILING_START( 'PZHEEVD_DESC' )
       IF ( SIZE( A , KIND=qi8 ) == 0_qi8 ) THEN
          CALL vtutor%error("PZHEEVD_DESC reports: Rank " // str(COMM%NODE_ME) // " has no data")
       ENDIF

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1), IWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))

       CALL PZHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q)))
       LIWORK= MAX(1,IWORK(1))
       LRWORK= MAX(1,NINT(RWORK(1)))

       DEALLOCATE(RWORK, WWORK, IWORK)
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))

       CALL PZHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )

       A=Z
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

       IF (INFO<0) THEN
          CALL vtutor%error("ERROR in diagonalization PZHEEVD: IFAIL= " // str(INFO))
       ELSE IF (INFO>0) THEN
          IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,*) "WARNING: PZHEEVD failed using fallback PZHEEV"
          CALL PZHEEV_DESC(A, W, N, DESCA, COMM)
       ENDIF
! just in case, PDSYEVD does seem to be buggy and sometimes not all nodes
! have the right eigenvalues in particular if matrices are small and stored
! on a single core
       CALLMPI ( M_bcast_d_from(COMM , W(1), SIZE(W), 1))

       PROFILING_STOP( 'PZHEEVD_DESC' )
    END SUBROUTINE PZHEEVD_DESC

!=======================================================================
!
! matrix diagonalization using divide and conquer scaLAPACK
! this version uses a two-dimensional input array
!
!=======================================================================

     SUBROUTINE PDSYEV_ZHEEVD_DESC2(A, W, N, DESCA_, COMM, JOBIN)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               ! dimension of matrix
       GDEF    A(:,:)          ! input/output matrix
       REAL(q) W(N)            ! eigenvalues
       INTEGER DESCA_( DLEN_ ) ! distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO

       GDEF, ALLOCATABLE ::     Z(:,:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1), IWORK(1))
       ALLOCATE(Z(SIZE(A,1), SIZE(A,2)))

#ifdef gammareal
       CALL PDSYEVD( JOB, 'U', N, A(1,1), 1, 1, DESCA, &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK, IWORK, LIWORK, INFO )
       RWORK(1)=1
       WWORK(1)=WWORK(1)*2 ! often crashes without this mod
#else
       CALL PZHEEVD( JOB, 'U', N, A(1,1), 1, 1, DESCA,  &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
#endif
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q)))
       LIWORK= MAX(1,IWORK(1))
       LRWORK= MAX(1,NINT(RWORK(1)))

       DEALLOCATE(RWORK, WWORK, IWORK)
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))

#ifdef gammareal
       CALL PDSYEVD( JOB, 'U', N, A(1,1), 1, 1, DESCA, &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK, IWORK, LIWORK, INFO )
#else
       CALL PZHEEVD( JOB, 'U', N, A(1,1), 1, 1, DESCA,  &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
#endif

       A=Z
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

       IF (INFO<0) THEN
          CALL vtutor%error("ERROR in diagonalization PDSYEVD/ PZHEEVD: IFAIL= " // str(INFO))
       ELSE IF (INFO>0) THEN
          IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,*) "WARNING: PDSYEVD/ PZHEEVD failed using fallback PDSYEV/ PZHEEV"
          CALL PDSYEV_ZHEEV_DESC2(A, W, N, DESCA, COMM)
       ENDIF
! just in case, PDSYEVD does seem to be buggy and sometimes not all nodes
! have the right eigenvalues in particular if matrices are small and stored in
! single core
       CALLMPI ( M_bcast_d_from(COMM , W(1), SIZE(W), 1))


    END SUBROUTINE PDSYEV_ZHEEVD_DESC2

!=======================================================================
!
! matrix diagonalization using divide and conquer scaLAPACK
! single precision version
!
!=======================================================================

     SUBROUTINE PSSYEV_CHEEVD_DESC(A, W, N, DESCA_, COMM, JOBIN)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N                ! dimension of matrix
       GDEFS    A(:)            ! input/output matrix
       REAL(qs) W(N)            ! eigenvalues
       INTEGER DESCA_( DLEN_ )  ! distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO

       GDEFS, ALLOCATABLE ::     Z(:)
       REAL(qs), ALLOCATABLE ::  RWORK(:) ! work array
       GDEFS, ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER, ALLOCATABLE ::  IWORK(:) ! integer work array
       INTEGER :: LWWORK, LIWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       PROFILING_START( 'PSSYEV_CHEEVD_DESC' )
       IF ( SIZE( A , KIND=qi8 ) == 0_qi8 ) THEN
          CALL vtutor%error("PSSYEV_CHEEVD_DESC reports: Rank " // str(COMM%NODE_ME) // " has no data")
       ENDIF

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LIWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1), IWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))

#ifdef gammareal
       CALL PSSYEVD( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, IWORK, LIWORK, INFO )
       RWORK(1)=1
       WWORK(1)=WWORK(1)*2 ! often crashes without this mod
#else
       CALL PCHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
#endif
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=qs)))
       LIWORK= MAX(1,IWORK(1))
       LRWORK= MAX(1,NINT(RWORK(1)))

       ! allocate optimal work space
       DEALLOCATE(RWORK, WWORK, IWORK)
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK), IWORK(LIWORK))

#ifdef gammareal
       CALL PSSYEVD( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, IWORK, LIWORK, INFO )
#else
       CALL PCHEEVD( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, IWORK, LIWORK, INFO )
#endif

       A=Z
       DEALLOCATE(Z)
       DEALLOCATE(IWORK, WWORK, RWORK)

       IF (INFO<0) THEN
          CALL vtutor%error("ERROR in diagonalization PSSYEVD/ PCHEEVD: IFAIL= " // str(INFO))
       ELSE IF (INFO>0) THEN
          IF (COMM%NODE_ME==COMM%IONODE) WRITE(*,*) "WARNING: PSSYEVD/ PCHEEVD failed using fallback PSSYEV/ PCHEEV"
          CALL PSSYEV_CHEEV_DESC(A, W, N, DESCA, COMM)
       ENDIF
! just in case, PSSYEVD does seem to be buggy and sometimes not all nodes
! have the right eigenvalues in particular if matrices are small and stored
! on a single core
       CALLMPI ( M_bcast_s_from(COMM , W(1), SIZE(W), 1))

       PROFILING_STOP( 'PSSYEV_CHEEVD_DESC' )
    END SUBROUTINE PSSYEV_CHEEVD_DESC

!=======================================================================
!
! matrix diagonalization using standard scaLAPACK  diagonalization
! usually not very fast but robust
!
!=======================================================================

     SUBROUTINE PDSYEV_ZHEEV_DESC(A, W, N, DESCA_, COMM, JOBIN)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               ! dimension of matrix
       GDEF    A(:)            ! input/output matrix
       REAL(q) W(N)            ! eigenvalues
       INTEGER DESCA_( DLEN_ ) ! distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO
       GDEF, ALLOCATABLE ::     Z(:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER :: LWWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))

#ifdef gammareal
       CALL PDSYEV( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, INFO )
       RWORK(1)=1
#else
       CALL PZHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
#endif
!jF:  added N to LWWORK ("heuristic"!) since sometimes still crashes happened
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q))) + N
       LRWORK= MAX(1,NINT(RWORK(1)))*2  ! crashes for small matrices otherwise

       DEALLOCATE(RWORK, WWORK )
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK))

#ifdef gammareal
       CALL PDSYEV( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, INFO )
#else
       CALL PZHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
#endif
       IF (INFO.NE.0) THEN
          CALL vtutor%error("ERROR in diagonalization PDSYEV/ PZHEEV: IFAIL= " // str(INFO))
       ENDIF

       A=Z
       DEALLOCATE(Z)
       DEALLOCATE(WWORK, RWORK)

    END SUBROUTINE PDSYEV_ZHEEV_DESC

!=======================================================================
!
! matrix diagonalization using standard scaLAPACK  diagonalization
! complex version regardless of gammareal
!
!=======================================================================

    SUBROUTINE PZHEEV_DESC(A, W, N, DESCA_, COMM, JOBIN)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               ! dimension of matrix
       COMPLEX(q)    A(:)            ! input/output matrix
       REAL(q) W(N)            ! eigenvalues
       INTEGER DESCA_( DLEN_ ) ! distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO
       COMPLEX(q), ALLOCATABLE ::     Z(:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       COMPLEX(q), ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER :: LWWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))

       CALL PZHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
!jF:  added N to LWWORK ("heuristic"!) since sometimes still crashes happened
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q))) + N
       LRWORK= MAX(1,NINT(RWORK(1)))*2  ! crashes for small matrices otherwise

       DEALLOCATE(RWORK, WWORK )
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK))

       CALL PZHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
       IF (INFO.NE.0) THEN
          CALL vtutor%error("ERROR in diagonalization PZHEEV: IFAIL= " // str(INFO))
       ENDIF

       A=Z
       DEALLOCATE(Z)
       DEALLOCATE(WWORK, RWORK)

    END SUBROUTINE PZHEEV_DESC

!=======================================================================
!
! matrix diagonalization using standard scaLAPACK  diagonalization
! this version uses a two-dimensional input array
!
!=======================================================================

     SUBROUTINE PDSYEV_ZHEEV_DESC2(A, W, N, DESCA_, COMM, JOBIN)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               ! dimension of matrix
       GDEF    A(:,:)          ! input/output matrix
       REAL(q) W(N)            ! eigenvalues
       INTEGER DESCA_( DLEN_ ) ! distributed matrix descriptor array
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO
       GDEF, ALLOCATABLE ::     Z(:,:)
       REAL(q), ALLOCATABLE ::  RWORK(:) ! work array
       GDEF, ALLOCATABLE    ::  WWORK(:) ! work array
       INTEGER :: LWWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1))
       ALLOCATE(Z(SIZE(A,1), SIZE(A,2)))

#ifdef gammareal
       CALL PDSYEV( JOB, 'U', N, A(1,1), 1, 1, DESCA, &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK, INFO )
       RWORK(1)=1
#else
       CALL PZHEEV( JOB, 'U', N, A(1,1), 1, 1, DESCA,  &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
#endif
!jF:  added N to LWWORK ("heuristic"!) since sometimes still crashes happened
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=q))) + N
       LRWORK= MAX(1,NINT(RWORK(1)))*2  ! crashes for small matrices otherwise

       DEALLOCATE(RWORK, WWORK )
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK))

#ifdef gammareal
       CALL PDSYEV( JOB, 'U', N, A(1,1), 1, 1, DESCA, &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK, INFO )
#else
       CALL PZHEEV( JOB, 'U', N, A(1,1), 1, 1, DESCA,  &
            W, Z(1,1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
#endif
       IF (INFO.NE.0) THEN
          CALL vtutor%error("ERROR in diagonalization PDSYEV/ PZHEEV: IFAIL= " // str(INFO))
       ENDIF

       A=Z
       DEALLOCATE(Z)
       DEALLOCATE(WWORK, RWORK)

    END SUBROUTINE PDSYEV_ZHEEV_DESC2


!=======================================================================
!
! matrix diagonalization using standard scaLAPACK  diagonalization
! single precision version
!
!=======================================================================

     SUBROUTINE PSSYEV_CHEEV_DESC(A, W, N, DESCA_, COMM, JOBIN)
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               ! dimension of matrix
       GDEFS   A(:)            ! input/output matrix
       REAL(qs) W(N)           ! eigenvalues
       INTEGER :: DESCA_( DLEN_ )
       TYPE(communic) :: COMM
       CHARACTER (LEN=1), OPTIONAL :: JOBIN
! local variables
       INTEGER INFO
       GDEFS, ALLOCATABLE ::     Z(:)
       REAL(qs), ALLOCATABLE :: RWORK(:) ! work array
       GDEFS, ALLOCATABLE   ::  WWORK(:) ! work array
       INTEGER :: LWWORK, LRWORK
       INTEGER :: DESCA( DLEN_ )
       CHARACTER (LEN=1) :: JOB

       JOB='V'  ! default is eigenvalues and eigenvectors
       IF (PRESENT(JOBIN)) JOB=JOBIN

       DESCA=DESCA_
       DESCA(M_)=N
       DESCA(N_)=N

       LWWORK=-1
       LRWORK=-1
       ALLOCATE(RWORK(1), WWORK(1))
       ALLOCATE(Z(SIZE(A,KIND=qi8)))

#ifdef gammareal
       CALL PSSYEV( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, INFO )
       RWORK(1)=1
#else
       CALL PCHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
#endif
!jF:  added N to LWWORK ("heuristic"!) since sometimes still crashes happened
       LWWORK= MAX(1,NINT(REAL(WWORK(1),KIND=qs))) + N
       LRWORK= MAX(1,NINT(RWORK(1)))*2  ! crashes for small matrices otherwise

       DEALLOCATE(RWORK, WWORK )
       ALLOCATE(RWORK(LRWORK), WWORK(LWWORK))

#ifdef gammareal
       CALL PSSYEV( JOB, 'U', N, A(1), 1, 1, DESCA, &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK, INFO )
#else
       CALL PCHEEV( JOB, 'U', N, A(1), 1, 1, DESCA,  &
            W, Z(1), 1, 1, DESCA, &
            WWORK, LWWORK,  RWORK, LRWORK, INFO )
#endif
       IF (INFO.NE.0) THEN
          CALL vtutor%error("ERROR in diagonalization PSSYEV/ PCHEEV: IFAIL= " // str(INFO))
       ENDIF

       A=Z
!-----------------------------------------------------------------------
! deallocation
!-----------------------------------------------------------------------
       DEALLOCATE(Z)
       DEALLOCATE(WWORK, RWORK)

    END SUBROUTINE PSSYEV_CHEEV_DESC


!=======================================================================
!
! the routine makes an (approximately unitary) transformation matrix
! exactly unitary (or orthogonal)
!
! if the orbitals (observing <n|m>= delta_mn )
! are transformed by a matrix A, they will no longer be
! exactly orthonormal
!    |n'> = |n> A_nn'  (see LINCOM)
! but rather observe
!    <n'| m'> =  A^+_n'm A_mm' = (A+ A)_n'm'
!
! The routine calculate S=  A+ A then determines
!    S= U+ U  (U upper triangular matrix)
! inverts U and updates the transformation matrix to
!    A <-  A U^-1
! one can easily show that
!   (A U^-1)^+ (A U^-1) = 1
!
!=======================================================================

     SUBROUTINE MAKE_UNITARY_DESC(A, N, DESCA )
       USE string, ONLY: str
       USE tutor, ONLY: vtutor
       IMPLICIT NONE

       INTEGER N               ! dimension of matrix
       GDEF    A(:)            ! input/output matrix
       INTEGER DESCA( DLEN_ )  ! distributed matrix descriptor array
       INTEGER :: INFO
! local variables
       GDEF, ALLOCATABLE ::     Z(:)

       ALLOCATE(Z(SIZE(A,KIND=qi8)))

       ! calculate S= A+ A
       CALL PGGEMM( trans, 'N', N, N, N, one,  &
          A(1), 1, 1, DESCA, &
          A(1), 1, 1, DESCA, zero,  &
          Z(1), 1, 1, DESCA )

       ! determine S= U+ U  (U upper triangular matrix)
#ifdef gammareal
       CALL PDPOTRF &
#else
       CALL PZPOTRF &
#endif
          & ('U', N, Z(1), 1, 1, DESCA, INFO)

       IF (INFO/=0) THEN
          CALL vtutor%bug("internal error in VASP: MAKE_UNITARY_DESC PZPOTRF failed " // str(INFO), __FILE__, __LINE__)
       ENDIF
       ! invert U <- U^-1
#ifdef gammareal
       CALL PDTRTRI &
#else
       CALL PZTRTRI &
#endif
          & ('U','N', N, Z(1), 1, 1, DESCA, INFO)
       IF (INFO/=0) THEN
          CALL vtutor%bug("internal error in VASP: MAKE_UNITARY_DESC PZTRTRI failed " // str(INFO), __FILE__, __LINE__)
       ENDIF
! triangular update of matrix A <- A U^-1
! note that the lower part of Z is set so DGEMM does not work here
#ifdef gammareal
       CALL PDTRMM &
#else
       CALL PZTRMM &
#endif
     &                ('R', 'U', 'N', 'N' , N, N, one, &
     &                 Z, 1, 1, DESCA, A, 1, 1, DESCA)

       DEALLOCATE(Z)

     END SUBROUTINE MAKE_UNITARY_DESC

!=======================================================================
!
! add diagonal elements AD (usually eigenvalues) to the supplied
! distributed matrix
!
!=======================================================================


    SUBROUTINE ADD_TO_DIAGONALE(N, A, AD, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      COMPLEX(q) A(:)        ! distributed matrix
      COMPLEX(q) AD(:)       ! diagonal elements to be added to A
      INTEGER N              ! dimension of the matrix
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                A(IROW_JCOL) = A(IROW_JCOL)+AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE ADD_TO_DIAGONALE

!
! same as above for a real valued function AD
!

    SUBROUTINE ADD_TO_DIAGONALE_REAL(N, A, AD, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      GDEF    A(:)       ! distributed matrix
      REAL(q) AD(:)      ! diagonal elements to be added to A
      INTEGER N          ! dimension of the matrix
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                A(IROW_JCOL) = A(IROW_JCOL)+AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE ADD_TO_DIAGONALE_REAL


!
! same as above, but set diagonal part
!

    SUBROUTINE SET_DIAGONALE_REAL(N, A, AD, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      GDEF    A(:)       ! distributed matrix
      REAL(q) AD(:)      ! diagonal elements to be subtracted from A
      INTEGER N          ! dimension of the matrix
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                IF (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2 <= SIZE(AD) ) THEN
                   A(IROW_JCOL) = AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
                ENDIF
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SET_DIAGONALE_REAL


!=======================================================================
!
! determine diagonal elements AD of the supplied
! distributed matrix  A
! important note: global sum is required after the call !!
!
!=======================================================================


    SUBROUTINE DETERMINE_DIAGONALE(N, A, AD, DESCA)
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      COMPLEX(q) A(:)    ! distributed matrix
      COMPLEX(q) AD(:)   ! on exit: diagonal elements
      INTEGER N          ! dimension of the matrix
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AD=0

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=A(IROW_JCOL)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE DETERMINE_DIAGONALE

!
! this version is for a GDEF input array (REAL(q) for gammaonly)
!

    SUBROUTINE DETERMINE_DIAGONALE_GDEF(N, A, AD, DESCA)
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      GDEF    A(:)       ! distributed matrix
      GDEF   AD(:)       ! on exit: diagonal elements
      INTEGER N          ! dimension of the matrix
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO
      PROFILING_START( 'determine_diagonale_gdef' )

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AD=0

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF(  (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) == &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                AD(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)=A(IROW_JCOL)
                GOTO 100 ! next column, since we have found the right row
              ENDIF
            ENDDO
          ENDDO
100       JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      PROFILING_STOP( 'determine_diagonale_gdef' )
      RETURN
    END SUBROUTINE DETERMINE_DIAGONALE_GDEF


!=========================================================================
! computes the taylor expansion of ( f_r - f_s ) / ( e_r - e_s )
! by expanding f_s around e_r, that is
! f_r ~ f_s + f'_s( e_r - e_s ) + f''_s( e_r - e_s)^2 / 2
!=========================================================================
    FUNCTION TAYLOR_FR_FS_OVER_ER_ES( F_R, F_S, E_R, E_S, BETA )
        REAL(q),INTENT(IN) :: F_R
        REAL(q),INTENT(IN) :: F_S
        REAL(q),INTENT(IN) :: E_R
        REAL(q),INTENT(IN) :: E_S
        REAL(q),INTENT(IN) :: BETA
        REAL(q) :: TAYLOR_FR_FS_OVER_ER_ES
        ! local

        TAYLOR_FR_FS_OVER_ER_ES = BETA*F_S*(1-F_S)*&
          (BETA*(E_R-E_S)/2*(1-2*F_S)-1)

    END FUNCTION TAYLOR_FR_FS_OVER_ER_ES

!=========================================================================
! computes the taylor expansion of ( f_r*e_r - f_s*e_s ) / ( e_r - e_s )
! by expanding f_s around e_r, that is
! f_r ~ g_s - beta*f_s( 1 - f_s )*( e_r = e_m * beta /2 ( e_m - e_n ) * (1-2*f_s)
!=========================================================================
    FUNCTION TAYLOR_FR_ER_FS_ES_OVER_ER_ES( F_R, F_S, E_R, E_S, BETA )
        REAL(q),INTENT(IN) :: F_R
        REAL(q),INTENT(IN) :: F_S
        REAL(q),INTENT(IN) :: E_R
        REAL(q),INTENT(IN) :: E_S
        REAL(q),INTENT(IN) :: BETA
        REAL(q) :: TAYLOR_FR_ER_FS_ES_OVER_ER_ES
        ! local

        TAYLOR_FR_ER_FS_ES_OVER_ER_ES = F_S - BETA * F_S * (1-F_S)*&
           ( E_R + BETA*E_S * ( E_S - E_R ) * ( 1 - 2*F_S)/2 )
    END FUNCTION TAYLOR_FR_ER_FS_ES_OVER_ER_ES

!=======================================================================
!
! divide distribute matrix by eigenvalue difference
!
!   (       |       )         (                |                )
!   (       |       )         ( A_ij*(f_i-f_j) | A_ia*(f_i-f_a) )
!   (  A_ij | A_ia  )         ( -------------- | -------------- )
!   (       |       )         (    e_i -e_j    |    e_i -e_a    )
!   (_______|_______)  ----\  (________________|________________)
!   (       |       )  ----/  (                |                )
!   (       |       )         ( A_ai*(f_a-f_i) | A_ab*(f_a-f_b) )
!   ( A_ai  | A_ab  )         ( -------------- | -------------- )
!   (       |       )         (    e_a -e_i    |    e_a -e_b    )
!   (       |       )         (                |                )
!
! this is usually required in standard perturbation theory
! the operation is only done in the occupied-unoccupied part
! of the matrix
! the first version is used to construct the density matrix
! from a given perturbation
! the second version is used to construct an approximately
! unitary rotation matrix
!
!=======================================================================

    SUBROUTINE DIVIDE_BY_EIGENVALUEDIFF(N, A, CELTOT, FERTOT, DESCA, &
      LFINITE_T, AUXTOT )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )! distributed matrix descriptor array
      GDEF    A(:)          ! distributed matrix
      COMPLEX(q) CELTOT(:)  ! eigenvalues
      REAL(q)    FERTOT(:)  ! occupancies
      INTEGER N             ! dimension of the matrix
      LOGICAL :: LFINITE_T  ! finite temperature mode
      REAL(q)    AUXTOT(:)  ! derivatives of occupancies
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
               IROW_JCOL = IROW_JCOL+1

               IF ( LFINITE_T ) THEN
                  CALL SET_A_RS_AT_T( A(IROW_JCOL),        &
                     FERTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2),        &
                     FERTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2),        &
                     REAL(CELTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2),q),&
                     REAL(CELTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2),q),&
                     AUXTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2),        &
                     AUXTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2))
               ELSE
                  CALL SET_A_RS_AT_T0( A(IROW_JCOL),       &
                     FERTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2),        &
                     FERTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2),        &
                     REAL(CELTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2),q),&
                     REAL(CELTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2),q) )

               ENDIF

            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO
      RETURN

      CONTAINS

      ! ------------------------------------------------------------------------
      ! gives matrix element at T=0, F_R and F_S has integer values 0 and 1 only
      ! ------------------------------------------------------------------------
      SUBROUTINE SET_A_RS_AT_T0( A, FI, FJ, EI, EJ )
         GDEF, INTENT(INOUT) :: A
         REAL(q),INTENT(IN)  :: FI
         REAL(q),INTENT(IN)  :: FJ
         REAL(q),INTENT(IN)  :: EI
         REAL(q),INTENT(IN)  :: EJ
         ! local
         REAL(q), PARAMETER  :: DIFMAX=0.001_q
         REAL(q) :: DIFF

         ! identical occupancies, clear rotation matrix
         DIFF=( FI- FJ )

         IF (ABS(DIFF)<DIFMAX) THEN
            A = 0
         ELSE
            ! (e_i - e_j)/(f_i-f_j)
            DIFF = ( EI - EJ ) / DIFF

            ! huge rotation element, should never occur except for tiny
            ! smearing parameters of the order of DIFMAX
            IF (ABS(DIFF)<DIFMAX) THEN
               A = 0
            ELSE
               A = A /DIFF
            ENDIF
         END IF

      END SUBROUTINE SET_A_RS_AT_T0

      ! ------------------------------------------------------------------------
      ! gives matrix element at T>0, FI and FJ between 0 and 1
      ! ------------------------------------------------------------------------
      SUBROUTINE SET_A_RS_AT_T( A, FI, FJ, EI, EJ, DFI, DFJ )
         GDEF, INTENT(INOUT) :: A
         REAL(q),INTENT(IN)  :: FI
         REAL(q),INTENT(IN)  :: FJ
         REAL(q),INTENT(IN)  :: EI
         REAL(q),INTENT(IN)  :: EJ
         REAL(q),INTENT(IN)  :: DFI
         REAL(q),INTENT(IN)  :: DFJ
         ! local
         REAL(q), PARAMETER  :: DIFMAX=0.001_q

         IF (ABS( EI - EJ )<DIFMAX) THEN
            ! this was tested for small transition energies against equation below
            A = A * (DFI + DFJ)/2
         ELSE
            A =  A* ( FI- FJ ) / ( EI - EJ )
         END IF

      END SUBROUTINE SET_A_RS_AT_T

    END SUBROUTINE DIVIDE_BY_EIGENVALUEDIFF


!=======================================================================
!
! divide matrix A by eigenvalue difference times eigenvalue in the
! uocc_occ block ( i,j \in occ ; a,b \in uocc), remove uocc_uocc block
!
!   (       |       )         (          |           )
!   (       |       )         (          |  e_i*A_ia )
!   (  A_ij | A_ia  )         (   A_ij   |  -------- )
!   (       |       )         (          |  e_i-e_a  )
!   (_______|_______)  ----\  (__________|___________)  ---    A'
!   (       |       )  ----/  (          |           )  ---
!   (       |       )         ( e_i*A_ai |     0     )
!   ( A_ai  | A_ab  )         ( -------- |    0 0    )
!   (       |       )         ( e_i-e_a  |     0     )
!   (       |       )         (          |           )
!
! in index notation
!
! in index notation
!       A' = A_ij * ( e_i f_i - e_j f_j) / (e_i - e_j)
!          = A_ij/2 [ f_i+f_j + (e_i+e_j) (f_i-f_j) /(e_i - e_j)]
!
! at finite temperature one has 0 < f_i,  f_j < 1
! so careful approximation of the fraction has to be done
! here f is assumed to be a Fermi function
!
! and add it to matrix B:
! B=B+A'
!=======================================================================

    SUBROUTINE DIVIDE_BY_EIGENVALUEDIFFxEIGENVAL(N, A, B, CELTOT, FERTOT, DESCA,&
      LFINITE_T, AUXTOT)
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      GDEF    A(:)          ! distributed Hamilton matrix
      GDEF    B(:)          ! distributed result matrix
      COMPLEX(q) CELTOT(:)  ! eigenvalues
      REAL(q)    FERTOT(:)  ! occupancies
      REAL(q)    AUXTOT(:)  ! occupancies
      INTEGER N             ! dimension of the matrix
      LOGICAL LFINITE_T
      REAL(q) :: EFERMI
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTRINSIC  MIN
      INTEGER    I, J
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
               IROW_JCOL = IROW_JCOL+1
               ! global column index is DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
               ! global row index DESCA(MB_)*MYROW+NPROW*(I1-1)+I2

               J = DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
               I = DESCA(MB_)*MYROW+NPROW*(I1-1)+I2

               ! compute matrix elements
               IF ( LFINITE_T ) THEN
                  CALL MATRIX_ELEMENT_AT_T( B(IROW_JCOL), A(IROW_JCOL),&
                      FERTOT(J), FERTOT(I), REAL(CELTOT(J),q), REAL(CELTOT(I),q), &
                      AUXTOT(I), AUXTOT(J))
               ELSE
                  CALL MATRIX_ELEMENT_AT_T0( B(IROW_JCOL), A(IROW_JCOL),&
                       FERTOT(J), FERTOT(I), REAL(CELTOT(J),q), REAL(CELTOT(I),q)  )
               ENDIF

            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      CONTAINS
      ! ------------------------------------------------------------------------
      ! gives matrix element at T=0, FI and FJ have integer values 0 and 1 only
      ! ------------------------------------------------------------------------
      SUBROUTINE MATRIX_ELEMENT_AT_T0( B, A, FJ, FI, EJ, EI )
         GDEF,INTENT(INOUT) :: B
         GDEF,INTENT(IN) :: A
         REAL(q), INTENT(IN) :: FJ, FI
         REAL(q), INTENT(IN) :: EJ, EI
         !local
         REAL(q), PARAMETER  :: DIFMAX=0.001_q

         IF (ABS(EJ-EI)<DIFMAX) THEN
            ! not the exact solution but a sensible simple approximation at T=0
            B = B +  A/2 * ( FI +FJ )
         ELSE
            B = B + A * ( FJ*EJ - FI*EI ) / (EJ-EI)
         END IF
      END SUBROUTINE MATRIX_ELEMENT_AT_T0

      ! ------------------------------------------------------------------------
      ! the matrix element at T>0 is obtained from
      ! sum_n i w_n g( iw_n) g( iw_n ) e^{-iw_n eta}
      ! = 1/2 ( sign( eta) -1) + (f_n e_n - f_m e_m ) / ( e_n - e_m )
      ! ------------------------------------------------------------------------
      SUBROUTINE MATRIX_ELEMENT_AT_T( B, A, FJ, FI, EJ, EI, DFJ, DFI )
         GDEF,INTENT(INOUT) :: B
         GDEF,INTENT(IN) :: A
         REAL(q), INTENT(IN) :: FJ, FI
         REAL(q), INTENT(IN) :: EJ, EI
         REAL(q), INTENT(IN) :: DFJ, DFI
         !local
         REAL(q), PARAMETER  :: DIFMAX=0.001_q
! the actual Matsubara sum is 
! 
! 1/beta \sum_n ( iw_n ) G(iw_n + mu ) F G( iw_n + mu) e^{-i w_n eta)
!    fj ej - fi ei
! = ---------------- 
!       ej - ei 
!
! This result follows from the evaluation of the series for ~ -w_n^2 
! and the series for ~ -I w_n^3 Exp[-iw_n eta] (eta>0)
! The term proportional to -w_n^2 gives the result -1/2, while the term 
! proportional to -I w_n^3 Exp[-iw_n 0-] yiels +1/2.
!
! For terms on the diagonal ei ~ ej 
! make a Taylor exansion of fj ej  around ei 
! 
! fj ej ~ [fi + fi' (ej-ei)  ][ ei + (ej-ei) ]
!       ~ fi ei + fi' (ej -ei ) ei + fi (ej-ei) 
!
! and insert this in the first term of the result giving 
! 
!    fi ei + fi' (ej -ei ) ei + fi (ej-ei) - fi ei
! = -----------------------------------------------  = fi' ei + fi 
!                        ej - ei 
!
! symmetrizing yields 
!
!    fi + fj    fi' ei + fj' ej
! =  ------- +  ---------------
!       2              2
! 
! to include off-diagonal terms the seond term is symmetrized as follows
! not this is valid since ei ~ ej holds true
!
!    fi + fj    fi' ei + fj' ej + fi' ej + fj' ei
! =  ------- +  ---------------------------------
!       2                       4
! 
         ! energies close to each other
         IF (ABS(EJ-EI)<DIFMAX) THEN
            !B = B +  A* ( FI+FJ + EI * DFI + EJ*DFJ ) / 2
            ! taylor expansion if occupancies are close (tested against below)
            B = B +  A* ( FI +FJ + (EI +EJ) * (DFI + DFJ)/2) / 2
         ELSE
            !B = B + A*( FJ*EJ - FI*EI )/(EJ-EI)
            B = B + A*(( FJ*EJ - FI*EI )/(EJ-EI))
         ENDIF
       END SUBROUTINE MATRIX_ELEMENT_AT_T

    END SUBROUTINE DIVIDE_BY_EIGENVALUEDIFFxEIGENVAL

!=======================================================================
!
! divide distribute matrix by eigenvalue difference
! this is usually required in standard perturbation theory
! the operation is only done in the occupied-unoccupied part
! of the matrix
!
!=======================================================================

    SUBROUTINE DIVIDE_BY_EIGENVALUEDIFF_SIGNED(N, A, CELTOT, DESCA, LFINITE_T )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      GDEF    A(:)          ! distributed matrix
      GDEF  CELTOT(:)  ! eigenvalues
      INTEGER N             ! dimension of the matrix
      LOGICAL :: LFINITE_T
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      REAL(q), PARAMETER  :: DIFMAX=0.001_q
      REAL(q) :: DIFCELL
      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1

              DIFCELL=REAL(CELTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)- &
                           CELTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2))
              IF (ABS(DIFCELL)<DIFMAX) THEN
                 A(IROW_JCOL) = 0
              ELSE
                 ! the equations above also involve
                 ! the difference between the occupancies, for the density matrix they would drop out
                 ! in other words: rotations in the (un)occupied-(un)occupied  block are irrelevant
                 ! for energies
                 ! tests indicate that inclusion of DIFCELL makes indeed little difference for energies
                 A(IROW_JCOL) = A(IROW_JCOL)/DIFCELL
                 IF (ABS(A(IROW_JCOL))>0.4_q) THEN
                    ! damp rotation matrix as in Loewdin case
                    A(IROW_JCOL)=A(IROW_JCOL)/ABS(ABS(A(IROW_JCOL)))*0.4_q
                 ENDIF
              ENDIF
              ! set diagonal to 1
              IF (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2==DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) THEN
                 A(IROW_JCOL)=1
              ENDIF


            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE DIVIDE_BY_EIGENVALUEDIFF_SIGNED

!=======================================================================
!
! calculate |H_ai|^2 / (epsilon_a-epsilon_i) *(f_a-f_i)
! both the matrix H_ai and the eigenvalues, and the occupancies are passed
! down
! the calculation is only performed for |epsilon_a - epsilon_i| >DIFFMAX
!
!=======================================================================

    SUBROUTINE SUM_SECOND_ORDER(N, A, CELTOT, FERTOT, CSUM, DESCA, LFINITE_T )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      GDEF    A(:)          ! distributed matrix
      REAL(q) CELTOT(:)     ! eigenvalues epsilon
      REAL(q) FERTOT(:)     ! occupancies
      COMPLEX(q) CSUM       ! final result |H_ai|^2 / (epsilon_a-epsilon_i)
      INTEGER N             ! dimension of the matrix
      LOGICAL :: LFINITE_T  ! finite temperature mode
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      REAL(q), PARAMETER  :: DIFMAX=0.001_q
      REAL(q) :: DIFCELL
      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CSUM=0
      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              DIFCELL=REAL(CELTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)- &
                          CELTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2))
              IF (ABS(DIFCELL)>=DIFMAX) THEN
                 CSUM=CSUM+A(IROW_JCOL)*GCONJG(A(IROW_JCOL))/DIFCELL* &
                 (FERTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)-FERTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2))
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SUM_SECOND_ORDER

!=======================================================================
!
! clear occupied-occupied and unoccupied-unoccupied block of
! a square matrix
!
!=======================================================================

    SUBROUTINE CLEAR_OCCOCC_UNOCCUNOCC(N, A, LAST_OCCUPIED, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      GDEF    A(:)          ! distributed matrix
      INTEGER N             ! dimension of the matrix
      INTEGER LAST_OCCUPIED ! last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) >LAST_OCCUPIED .AND.  &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)<= LAST_OCCUPIED) .OR. &
                  ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)<=LAST_OCCUPIED .AND. &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)>  LAST_OCCUPIED)) THEN
                     ! do nothing
               ELSE
                  ! remove the element from the perturbation
                  A(IROW_JCOL)=0
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE CLEAR_OCCOCC_UNOCCUNOCC

!=======================================================================
!
! clear occupied-unoccupied and unoccupied-occupied blocks of
! a square matrix
!
!=======================================================================

    SUBROUTINE CLEAR_OCCUNOCC_UNOCCOCC(N, A, LAST_OCCUPIED, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      GDEF    A(:)          ! distributed matrix
      INTEGER N             ! dimension of the matrix
      INTEGER LAST_OCCUPIED ! last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) >LAST_OCCUPIED .AND.  &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)<= LAST_OCCUPIED) .OR. &
                  ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)<=LAST_OCCUPIED .AND. &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)>  LAST_OCCUPIED)) THEN
                  ! remove the element from the perturbation
                  A(IROW_JCOL)=0
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE CLEAR_OCCUNOCC_UNOCCOCC


!=======================================================================
!
! clear occupied-unoccupied and unoccupied-occupied blocks of
! a square matrix using supplied fermi-weights
! specifically the supplied matrix is multiplied by
! ferwe(ncol)*ferwe(nrow) + (1-ferwe(ncol))*(1-ferwe(nrow))
!
!=======================================================================

    SUBROUTINE CLEAR_OCCUNOCC_UNOCCOCC_FER(N, A, FERTOT, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      COMPLEX(q) A(:)        ! distributed matrix
      INTEGER N              ! dimension of the matrix
      REAL(q)    FERTOT(:)   ! occupancies
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              A(IROW_JCOL)=A(IROW_JCOL)* &
                 (FERTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)* FERTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)+ &
                  (1-FERTOT(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2))*(1-FERTOT(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)))
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE CLEAR_OCCUNOCC_UNOCCOCC_FER

!=======================================================================
!
! set occupied-occupied block of a square matrix to identity matrix
!
!=======================================================================

    SUBROUTINE SET_OCCOCC_IDENTITY(N, A, LAST_OCCUPIED, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      GDEF    A(:)          ! distributed matrix
      INTEGER N             ! dimension of the matrix
      INTEGER LAST_OCCUPIED ! last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)<= LAST_OCCUPIED .AND.  &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)<= LAST_OCCUPIED) THEN
                  ! remove the element from the perturbation
                  A(IROW_JCOL)=0._q
                  IF ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) .EQ. (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)) THEN
                     A(IROW_JCOL)=1._q
                  ENDIF
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SET_OCCOCC_IDENTITY

!=======================================================================
!
! set occupied-occupied block to a square matrix
! n+1        0             0           0
! 0          n             0           0
! 0          0             n-1         0
! 0          0             0           1
! this will preserve the order of the states during diagonalization
! of the density matrix
!
!=======================================================================

    SUBROUTINE SET_OCCOCC_INCREASING(N, A, LAST_OCCUPIED, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      GDEF    A(:)          ! distributed matrix
      INTEGER N             ! dimension of the matrix
      INTEGER LAST_OCCUPIED ! last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)<= LAST_OCCUPIED .AND.  &
                   (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)<= LAST_OCCUPIED) THEN
                  ! remove the element from the perturbation
                  A(IROW_JCOL)=0._q
                  IF ((DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) .EQ. (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)) THEN
                     A(IROW_JCOL)=LAST_OCCUPIED+1-(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)
                  ENDIF
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SET_OCCOCC_INCREASING

!=======================================================================
!
! conjugate the upper triangle of a matrix
! this routine constructs from the change of the density matrix
! the unitary transformation for the orbitals
! one can also simply take the lower triangle and change the sign
! in the upper
!
!=======================================================================

    SUBROUTINE MINUS_UPPER(N, A, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array
      GDEF    A(:)          ! distributed matrix
      INTEGER N             ! dimension of the matrix
      INTEGER LAST_OCCUPIED ! last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL

          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
              IROW_JCOL = IROW_JCOL+1
              IF( (DESCA(MB_)*MYROW+NPROW*(I1-1)+I2) < &
                  (DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) ) THEN
                   A(IROW_JCOL)=-(A(IROW_JCOL))
              ENDIF
            ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE MINUS_UPPER


!=======================================================================
!
! multiply distribute matrix by eigenvalue difference
!   1/(iw - epsilon_i) and  1/(iw - epsilon_a)
! from the left and right
!
!=======================================================================


    SUBROUTINE MULTIPLY_BY_IVEIGENVALUE(N, A, CELTOTINV, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  ! distributed matrix descriptor array
      COMPLEX(q) A(:)         ! distributed matrix
      COMPLEX(q) CELTOTINV(:) ! multiplication factors e.g. 1/(iw-eigenvalue)
      INTEGER N               ! dimension of the matrix
      INTEGER LAST_OCCUPIED   ! last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*CELTOTINV(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) &
                                           *CELTOTINV(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE MULTIPLY_BY_IVEIGENVALUE

!
! use only real part of transformation
!
    SUBROUTINE MULTIPLY_BY_IVEIGENVALUE_REAL(N, A, CELTOTINV, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  ! distributed matrix descriptor array
      COMPLEX(q) A(:)         ! distributed matrix
      COMPLEX(q) CELTOTINV(:) ! multiplication factors e.g. 1/(iw-eigenvalue)
      INTEGER N               ! dimension of the matrix
      INTEGER LAST_OCCUPIED   ! last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*REAL(CELTOTINV(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2) &
                                           *CELTOTINV(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2),q)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE MULTIPLY_BY_IVEIGENVALUE_REAL


!=======================================================================
!
! right multiply matrix by a vector (most likely a BLACS routines
! exists to do this
!  A <- A R    or  A_ij = A_ij R_j
!
!=======================================================================
!
    SUBROUTINE RIGHT_MULTIPLY_BY_R(N, A, R, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  ! distributed matrix descriptor array
      GDEF A(:)               ! distributed matrix
      REAL(q) R(:)            ! multiplication factors R
      INTEGER N               ! dimension of the matrix
      INTEGER LAST_OCCUPIED   ! last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*R(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RIGHT_MULTIPLY_BY_R


    SUBROUTINE RIGHT_MULTIPLY_BY_C(N, A, R, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  ! distributed matrix descriptor array
      COMPLEX(q) A(:)         ! distributed matrix
      COMPLEX(q) R(:)         ! multiplication factors R
      INTEGER N               ! dimension of the matrix
      INTEGER LAST_OCCUPIED   ! last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*R(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RIGHT_MULTIPLY_BY_C

    SUBROUTINE LEFT_MULTIPLY_BY_R(N, A, R, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  ! distributed matrix descriptor array
      GDEF A(:)               ! distributed matrix
      REAL(q) R(:)            ! multiplication factors R
      INTEGER N               ! dimension of the matrix
      INTEGER LAST_OCCUPIED   ! last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*R(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE LEFT_MULTIPLY_BY_R

    SUBROUTINE LEFT_MULTIPLY_BY_C(N, A, R, DESCA )
      IMPLICIT NONE

      INTEGER DESCA( DLEN_ )  ! distributed matrix descriptor array
      COMPLEX(q) A(:)               ! distributed matrix
      COMPLEX(q) R(:)            ! multiplication factors R
      INTEGER N               ! dimension of the matrix
      INTEGER LAST_OCCUPIED   ! last band treated as occupied
 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW_JCOL
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      JCOL = 0
      ! loop over local columns
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          IROW_JCOL = JCOL
          ! loop over local rows
          DO I1=1,NP,DESCA(MB_)
             I1RES=MIN(DESCA(MB_),NP-I1+1)
             DO I2=1,I1RES
                IROW_JCOL = IROW_JCOL+1
                A(IROW_JCOL) = A(IROW_JCOL)*R(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)
             ENDDO
          ENDDO
          JCOL = JCOL + DESCA(LLD_)
       ENDDO
      ENDDO

      RETURN
    END SUBROUTINE LEFT_MULTIPLY_BY_C

!=======================================================================
!
! distribute slice of the global matrix AMATIN to a distributed matrix A
! as input the routine expects several columns of the matrix
! to be distributed
! the columns of the input matrix are COLUMN_LOW up to COLUMN_HIGH
! all rows are passed to the routine
!
!=======================================================================

    SUBROUTINE DISTRI_SLICE(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)
      INTEGER NDIM                                   ! first dimension of matrix
      INTEGER N                                      ! rank of matrix (not used except by CHECK_scala)
      INTEGER COLUMN_LOW, COLUMN_HIGH                ! starting and end row supplied in AMATIN
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  ! input matrix (globally known)
      GDEF    A(*)                                   ! on return: output distributed matrix
      INTEGER DESCA( DLEN_ )                         ! distributed matrix descriptor array
    ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      PROFILING_START('distri_slice')

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     WRITE(*,*) 'distri',NDIM, N

!     setup distributed matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          ! global column index: block_size*processor coordinate +grid dimension*block number + index in block
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN

             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
          ! global row index: block_size*processor coordinate +grid dimension*block number + index in block
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   A(IROW+(JCOL-1)*DESCA(LLD_))= &
                   AMATIN(ROW, COL-COLUMN_LOW+1)
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO
      DWRITE "DISTRI",ITEST,MYCOL,MYROW

      PROFILING_STOP('distri_slice')

      RETURN
    END SUBROUTINE DISTRI_SLICE

!
! distribute matrix and add hermitian conjugated elements
!
!
    SUBROUTINE DISTRI_SLICE_HERM(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)
      INTEGER NDIM                                   ! first dimension of matrix
      INTEGER N                                      ! rank of matrix (not used except by CHECK_scala)
      INTEGER COLUMN_LOW, COLUMN_HIGH                ! starting and end row supplied in AMATIN
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  ! input matrix (globally known)
      GDEF    A(*)                                   ! on return: output distributed matrix
      INTEGER DESCA( DLEN_ )                         ! distributed matrix descriptor array
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index
      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      PROFILING_START('distri_slice_herm')

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     set conjugated part first
      IROW=0
      DO I1=1,NP,DESCA(MB_)
         I1RES=MIN(DESCA(MB_),NP-I1+1)
         DO I2=1,I1RES
            IROW=IROW+1
          ! global row index: block_size*processor coordinate +grid dimension*block number + index in block
            ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
            IF (ROW>=COLUMN_LOW .AND. ROW<=COLUMN_HIGH) THEN
               JCOL=0
               DO J1=1,NQ,DESCA(NB_)
                  J1RES=MIN(DESCA(NB_),NQ-J1+1)
                  DO J2=1,J1RES
                     JCOL=JCOL+1
          ! global column index: block_size*processor coordinate +grid dimension*block number + index in block
                     COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
                       A(IROW+(JCOL-1)*DESCA(LLD_))=GCONJG(AMATIN(COL, ROW-COLUMN_LOW+1))
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

!     non conjugated part
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          ! global column index: block_size*processor coordinate +grid dimension*block number + index in block
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
          ! global row index: block_size*processor coordinate +grid dimension*block number + index in block
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                     A(IROW+(JCOL-1)*DESCA(LLD_))=AMATIN(ROW, COL-COLUMN_LOW+1)
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      PROFILING_STOP('distri_slice_herm')

      RETURN
    END SUBROUTINE DISTRI_SLICE_HERM

!
! distribute matrix and add hermitian conjugated elements
! in this version it is assumed that the lower triangle of the input matrix is set
! this requires some extra precaution that previously properly set elements
! are not overwritten
! only the lower part of the input matrix is used as source
!
    SUBROUTINE DISTRI_SLICE_HERM_L(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)
      INTEGER NDIM                                   ! first dimension of matrix
      INTEGER N                                      ! rank of matrix (not used except by CHECK_scala)
      INTEGER COLUMN_LOW, COLUMN_HIGH                ! starting and end row supplied in AMATIN
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  ! input matrix (globally known)
      GDEF    A(*)                                   ! on return: output distributed matrix
      INTEGER DESCA( DLEN_ )                         ! distributed matrix descriptor array
! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index
      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)
      DWRITE "DISTRI,NP,NQ",NP,NQ," MB,NB",DESCA(MB_),DESCA(NB_)," LLD",DESCA(LLD_)

!     set conjugated part
      IROW=0
      DO I1=1,NP,DESCA(MB_)
         I1RES=MIN(DESCA(MB_),NP-I1+1)
         DO I2=1,I1RES
            IROW=IROW+1
          ! global row index: block_size*processor coordinate +grid dimension*block number + index in block
            ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
            IF (ROW>=COLUMN_LOW .AND. ROW<=COLUMN_HIGH) THEN
               JCOL=0
               DO J1=1,NQ,DESCA(NB_)
                  J1RES=MIN(DESCA(NB_),NQ-J1+1)
                  DO J2=1,J1RES
                     JCOL=JCOL+1
          ! global column index: block_size*processor coordinate +grid dimension*block number + index in block
                     COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
                     IF (COL>=ROW) A(IROW+(JCOL-1)*DESCA(LLD_))=GCONJG(AMATIN(COL, ROW-COLUMN_LOW+1))
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

!     non conjugated part
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          ! global column index: block_size*processor coordinate +grid dimension*block number + index in block
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
          ! global row index: block_size*processor coordinate +grid dimension*block number + index in block
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   IF (ROW>=COL) A(IROW+(JCOL-1)*DESCA(LLD_))=AMATIN(ROW, COL-COLUMN_LOW+1)
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE DISTRI_SLICE_HERM_L

!=======================================================================
!
! reconstruct a slice of the originally distributed matrix AMATIN
!
!=======================================================================

    SUBROUTINE RECON_SLICE(AMATIN, NDIM, N, A, DESCA, COLUMN_LOW, COLUMN_HIGH)

      INTEGER DESCA( DLEN_ )             ! distributed matrix descriptor
      INTEGER COLUMN_LOW, COLUMN_HIGH    ! first and last column of the matrix to be merged
      INTEGER NDIM                       ! first dimension of matrix AMATIN
      INTEGER N                          ! rank of matrix (not used except by CHECK_scala)
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)
                                         ! on return: matrix stored in serial layout between column
                                         !   COLUMN_LOW and COLUMN_HIGH
      GDEF    A(*)                       ! distributed matrix
   ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      PROFILING_START('recon_slice')

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AMATIN=0.0_q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   AMATIN(ROW,COL-COLUMN_LOW+1)= &
                     A(IROW+(JCOL-1)*DESCA(LLD_))
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      PROFILING_STOP('recon_slice')

      RETURN
    END SUBROUTINE RECON_SLICE

!
! complex version of the previous routine
!

    SUBROUTINE RECON_SLICE_C(AMATIN, NDIM, N, A, DESCA, COLUMN_LOW, COLUMN_HIGH)

      INTEGER DESCA( DLEN_ )             ! distributed matrix descriptor
      INTEGER COLUMN_LOW, COLUMN_HIGH    ! first and last column of the matrix to be merged
      INTEGER NDIM                       ! first dimension of matrix AMATIN
      INTEGER N                          ! rank of matrix (not used except by CHECK_scala)
      COMPLEX(q) AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)
                                         ! on return: matrix stored in serial layout between column
                                         !   COLUMN_LOW and COLUMN_HIGH
      COMPLEX(q) A(*)                    ! distributed matrix
   ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AMATIN=0.0_q

!     rebuild patched matrix
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   AMATIN(ROW,COL-COLUMN_LOW+1)= &
                     A(IROW+(JCOL-1)*DESCA(LLD_))
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON_SLICE_C

!
! similar version as above but add complex conjugated elements
! in the lower triangle  in AMATIN
!
    SUBROUTINE RECON_SLICE_HERM(AMATIN,NDIM,N,A,DESCA, COLUMN_LOW, COLUMN_HIGH)

      INTEGER DESCA( DLEN_)              ! distributed matrix descriptor array
      INTEGER COLUMN_LOW, COLUMN_HIGH    ! first and last column of the matrix to be merged
      INTEGER NDIM                       ! first dimension of matrix AMATIN
      INTEGER N                          ! rank of matrix (not used except by CHECK_scala)
      GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)
                                         ! on return: matrix stored in serial layout between column
      GDEF    A(*)                       ! distributed matrix
   ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,IROW,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2
      INTEGER ROW,COL            ! global column and row index

      INTRINSIC  MIN
      INTEGER    NUMROC
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL CHECK_scala (N, 'RECON_SLICE_HERM' )

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)
      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      AMATIN=0.0_q

!     set conjugated part first restrict to lower triangle
      IROW=0
      DO I1=1,NP,DESCA(MB_)
         I1RES=MIN(DESCA(MB_),NP-I1+1)
         DO I2=1,I1RES
            IROW=IROW+1
          ! global row index: block_size*processor coordinate +grid dimension*block number + index in block
            ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
            IF (ROW>=COLUMN_LOW .AND. ROW<=COLUMN_HIGH) THEN
               JCOL=0
               DO J1=1,NQ,DESCA(NB_)
                  J1RES=MIN(DESCA(NB_),NQ-J1+1)
                  DO J2=1,J1RES
                     JCOL=JCOL+1
          ! global column index: block_size*processor coordinate +grid dimension*block number + index in block
                     COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
                     IF (ROW<COL) &
                     AMATIN(COL, ROW-COLUMN_LOW+1)= &
                        GCONJG(A(IROW+(JCOL-1)*DESCA(LLD_)))
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO

!     now the non conjugated part but restrict to upper triangle
      JCOL=0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL=JCOL+1
          COL=DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2
          IF (COL>=COLUMN_LOW .AND. COL<=COLUMN_HIGH) THEN
             IROW=0
             DO I1=1,NP,DESCA(MB_)
                I1RES=MIN(DESCA(MB_),NP-I1+1)
                DO I2=1,I1RES
                   IROW=IROW+1
                   ROW=DESCA(MB_)*MYROW+NPROW*(I1-1)+I2
                   IF (ROW<=COL) &
                      AMATIN(ROW,COL-COLUMN_LOW+1)= &
                            A(IROW+(JCOL-1)*DESCA(LLD_))
                ENDDO
             ENDDO
          ENDIF
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE RECON_SLICE_HERM

!=======================================================================
!
! shift clear off diagonal components below NLOW
!
!=======================================================================

    SUBROUTINE SCA_SHIFT_BANDS_LOW(N, NLOW, A, SHIFT, DESCA )
      IMPLICIT NONE

      INTEGER N              ! actual dimension of the matrix
      INTEGER NLOW           ! first band index not to be shifted
      GDEF    A(*)           ! distributed matrix
      REAL(q) SHIFT          ! magnitude of shift
      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array

 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2, JCOL_GLOBAL, IROW_GLOBAL

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      ! loop over local columns
      JCOL = 0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL = JCOL + 1
          JCOL_GLOBAL=(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)

          ! loop over local rows
          IROW = 0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
               IROW = IROW + 1
               IROW_GLOBAL=(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)

               IF (IROW_GLOBAL < NLOW .OR. JCOL_GLOBAL < NLOW ) THEN
                  IF (IROW_GLOBAL== JCOL_GLOBAL) THEN
                     A(IROW+ (JCOL-1)*DESCA(LLD_) )=A(IROW+ (JCOL-1)*DESCA(LLD_) )+SHIFT
                  ELSE
                     A(IROW+ (JCOL-1)*DESCA(LLD_) )=0.0_q
                  ENDIF
               ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SCA_SHIFT_BANDS_LOW


    SUBROUTINE SCA_SHIFT_BANDS_HIGH(N, NHIGH, A, SHIFT, DESCA )
      IMPLICIT NONE

      INTEGER N              ! dimension of the matrix
      INTEGER NHIGH          ! first band index not to be shifted
      GDEF    A(*)           ! distributed matrix
      REAL(q) SHIFT          ! magnitude of shift
      INTEGER DESCA( DLEN_ ) ! distributed matrix descriptor array

 ! local
      INTEGER NP,NQ
      INTEGER I1RES,J1RES,JCOL
      INTEGER MYROW, MYCOL, NPROW, NPCOL
      INTEGER I1,I2,J1,J2, JCOL_GLOBAL, IROW_GLOBAL

      INTRINSIC  MIN
      INTEGER    NUMROC,IROW
      EXTERNAL   NUMROC,BLACS_GRIDINFO

      CALL BLACS_GRIDINFO(DESCA(CTXT_),NPROW,NPCOL,MYROW,MYCOL)

      NP = NUMROC(N,DESCA(MB_),MYROW,0,NPROW)
      NQ = NUMROC(N,DESCA(NB_),MYCOL,0,NPCOL)

      ! loop over local columns
      JCOL = 0
      DO J1=1,NQ,DESCA(NB_)
        J1RES=MIN(DESCA(NB_),NQ-J1+1)
        DO J2=1,J1RES
          JCOL = JCOL + 1
          JCOL_GLOBAL=(DESCA(NB_)*MYCOL+NPCOL*(J1-1)+J2)

          ! loop over local rows
          IROW = 0
          DO I1=1,NP,DESCA(MB_)
            I1RES=MIN(DESCA(MB_),NP-I1+1)
            DO I2=1,I1RES
               IROW = IROW + 1
               IROW_GLOBAL=(DESCA(MB_)*MYROW+NPROW*(I1-1)+I2)

               IF (IROW_GLOBAL > NHIGH .OR. JCOL_GLOBAL > NHIGH ) THEN
                  IF (IROW_GLOBAL== JCOL_GLOBAL) THEN
                     A(IROW+ (JCOL-1)*DESCA(LLD_) )=A(IROW+ (JCOL-1)*DESCA(LLD_) )+SHIFT
                  ELSE
                     A(IROW+ (JCOL-1)*DESCA(LLD_) )=0.0_q
                  ENDIF
               ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
    END SUBROUTINE SCA_SHIFT_BANDS_HIGH

!=======================================================================
!
! performs a pseudo inverse of CHI1D
! CHI1D is passed as pointer to 1D array
!
! the following lines can be used to reshape a 2D array into 1D form
!      USE, INTRINSIC :: ISO_C_BINDING
!      GDEFCPTR, ALLOCATABLE   :: CHI1D(:)
!      CHI1D = RESHAPE( CHI_INV, (/MY_NROWS*MY_NCOLS/) )
!
!=======================================================================
   SUBROUTINE PSEUDO_INVERSE_DESC( CHI1D, NP, DESC, THRESHHOLD, IU6 )
      USE tutor, ONLY: vtutor
      USE string, ONLY: str
      GDEF, POINTER        :: CHI1D(:)
      INTEGER              :: NP               !dimension of submatrix
      INTEGER              :: DESC( DLEN_ )  ! distributed matrix descriptor array
      REAL(q)              :: THRESHHOLD
      INTEGER              :: IU6
      ! local
      INTEGER              :: J
      INTEGER              :: PINFO                !scaLAPACK routine info variable
      INTEGER              :: LWORK, LRWORK        !for scaLAPACK routine
      REAL(q)              :: E(NP)
      GDEF   ,ALLOCATABLE  :: EV(:)
      GDEF   ,ALLOCATABLE  :: V(:)
      GDEF, ALLOCATABLE    :: WORK(:)
#ifndef gammareal
      COMPLEX(q), ALLOCATABLE :: CWORK(:)             !needed for complex computation
#endif
      ALLOCATE(EV(SIZE(CHI1D)))
      ALLOCATE(V(SIZE(CHI1D)))
      ! initialize eigenvalues and eigenvectors
      E = zero
      EV = zero
      !
      !first get the optimal working arrays
      !
      ALLOCATE(WORK(1))
#ifdef gammareal
      !for real values only
      CALL PDSYEV( 'V', 'U', NP, CHI1D(1), 1, 1, DESC, E(1),&
           EV(1), 1, 1, DESC, WORK(1), -1, PINFO)
      IF ( PINFO /=0 ) THEN
         CALL vtutor%error("ERROR in PSEUDO_INVERSE_DESC for NP=" // str(NP) // " 1st call of &
            &PDSYEV returns " // str(PINFO))
      ENDIF
#else
      !for a hermitian matrix
      ALLOCATE(CWORK(1))
      CALL PZHEEV( 'V', 'U', NP, CHI1D(1), 1, 1, DESC, E(1), &
           EV(1), 1, 1, DESC, WORK(1), -1, CWORK(1), -1, PINFO )
      IF ( PINFO /=0 ) THEN
         CALL vtutor%error("ERROR in PSEUDO_INVERSE_DESC for NP=" // str(NP) // " 1st call of &
            &PZHEEV returns " // str(PINFO))
      ENDIF
#endif

      !LWORK is now the optimal size of WORK
      LWORK=WORK(1)
      IF ( WORK(1)  /= 0 ) THEN   !if it's not zero reallocate WORK array
         DEALLOCATE(WORK)
         ALLOCATE(WORK(LWORK))
      ENDIF

      !diagonalize
#ifdef gammareal
      !real version
      CALL PDSYEV( 'V', 'U', NP, CHI1D(1), 1, 1, DESC, E(1),&
           EV(1), 1, 1, DESC, WORK(1), LWORK , PINFO)
      IF ( PINFO < 0 ) THEN
         CALL vtutor%error("ERROR in PSEUDO_INVERSE_DESC for NP " // str(NP) // " 2nd call of &
            &PDSYEV returns " // str(PINFO))
      ENDIF
#else
      !in the complex version do the same with the CWORK array
      LRWORK=CWORK(1)
      IF ( LRWORK  /= 0 ) THEN   !if it's not zero reallocate WORK array
        DEALLOCATE(CWORK)
        ALLOCATE(CWORK(LRWORK))
      ENDIF
      !hermitian version
      CALL PZHEEV( 'V', 'U', NP, CHI1D(1), 1, 1, DESC, E(1),&
           EV(1), 1, 1, DESC, WORK(1), LWORK, CWORK(1), LRWORK, PINFO )
      IF ( PINFO < 0 ) THEN
         CALL vtutor%error("ERROR in PSEUDO_INVERSE_DESC for NP " // str(NP) // " 2nd call of &
            &PZHEEV returns " // str(PINFO))
      ENDIF
      DEALLOCATE( CWORK )
#endif

      IF ( PINFO > 0 ) THEN
         WRITE(*,'(" Eigenvector #",I4," may not be converged" )') PINFO
      ENDIF

      IF (IU6>=0) WRITE(IU6,'(A)') ' eigenvalues of the response function'
      IF (IU6>=0) WRITE(IU6,'(2X,10F12.6)') E(1:NP)

      ! invert eigenvalues based on passed threshhold
      DO J=1,NP
         IF (ABS(E(J))<THRESHHOLD) THEN
            E(J)=0
         ELSE
            ! uncomment the following line to represent original matrix
            E(J)=1/E(J)
         ENDIF
      ENDDO

      ! form diagonal matrix with eigenvalues
      V = zero
      CALL ADD_TO_DIAGONALE_REAL(NP, V, E,  DESC )
      ! re-shape back to 2D array
       ! EIG * V,  multiply with rotation matrix
      CALL PGGEMM( 'N', 'N', NP, NP, NP, one,  &
         EV(1), 1, 1, DESC,              &
         V(1), 1, 1, DESC,         &
         zero, CHI1D(1), 1, 1, DESC )
      V = CHI1D
       ! 1/EIG * V,  multiply with rotation matrix
      CALL PGGEMM( 'N', trans, NP, NP, NP, one,  &
         EV(1), 1, 1, DESC,              &
         V(1), 1, 1, DESC,         &
         zero, CHI1D(1), 1, 1, DESC )

      DEALLOCATE( WORK )
      DEALLOCATE( EV )
      DEALLOCATE( V )
   END SUBROUTINE PSEUDO_INVERSE_DESC


END MODULE
#endif


#ifdef scaLAPACK

!========================== ORTH1_DISTRI ===============================
!
! ORTH1_DISTRI is essentially identical to ORTH1
! it calculates a slice or stripe of the matrix
!   O(I,J) =  <CPTWFP(I) | 1+ Q | CFW(J) >
!
! and THEN DISTRIBUTES THE SLICE in a scaLAPACK compatible manner
! using DISTRI_SLICE
! ORTH1 requires the storage of a global NB_TOT time NB_TOT matrix
! whereas ORTH1_DISTRI only the storage of the distributed matrix
!
! arguments are identical to ORTH1, except for additional communicator
!
!=======================================================================

     SUBROUTINE ORTH1_DISTRI(CSEL,CPTWFP,CFW,CPROJ,CPROW, &
           NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,MY_COMM,NBANDSK)
        USE prec
        USE scala
        USE tutor, ONLY: vtutor
        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(*)
        CHARACTER*(*) CSEL

        PROFILING_START('orth1_distri')

        CALL CHECK_scala (NBANDSK, 'ORTH1_DISTRI' )

        COVL=0
        CALL ORTH1_NOSUBINDEX(CSEL,CPTWFP(1,1),CFW(1,1),CPROJ(1,1), &
           CPROW(1,1),NBANDS, &
           NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL(1,1))

        PROFILING_START('m_sum_g')
        CALLMPI( M_sum_g(MY_COMM,COVL(1,1),NBANDS*NSTRIP))
        PROFILING_STOP('m_sum_g')

#ifdef ELPA
        IF (LELPA) THEN
           IF (CSEL=='U') THEN
              CALL DISTRI_SLICE_HERM(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCSTD,NPOS,NPOS+NSTRIP-1)
           ELSE
              CALL vtutor%bug("internal error in VASP: ORTH1_DISTRI is only tested for CSEL=U try " &
                 // "DISTRI_SLICE_HERM_L", __FILE__, __LINE__)
           ENDIF
        ENDIF
#endif
        IF (.NOT. LELPA) THEN
           CALL DISTRI_SLICE(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCSTD,NPOS,NPOS+NSTRIP-1)
        ENDIF

        PROFILING_STOP('orth1_distri')

      END SUBROUTINE ORTH1_DISTRI

!
! same as above, this version always adds the Hermitian conjugated elements
! to the matrix
! costs little extra time, and is more secure, in the sense that
! we can use the constructed matrix in places where the lower triangle is required as well
!

      SUBROUTINE ORTH1_DISTRI_HERM(CSEL,CPTWFP,CFW,CPROJ,CPROW, &
           NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,MY_COMM,NBANDSK)
        USE prec
        USE scala
        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(*)
        CHARACTER*(*) CSEL

        CALL CHECK_scala (NBANDSK, 'ORTH1_DISTRI_HERM' )

        COVL=0
        CALL ORTH1_NOSUBINDEX(CSEL,CPTWFP(1,1),CFW(1,1),CPROJ(1,1), &
           CPROW(1,1),NBANDS, &
           NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL(1,1))

        CALLMPI( M_sum_g(MY_COMM,COVL(1,1),NBANDS*NSTRIP))

        IF (CSEL=='U') THEN
           CALL DISTRI_SLICE_HERM(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCSTD,NPOS,NPOS+NSTRIP-1)
        ELSE
           CALL DISTRI_SLICE_HERM_L(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCSTD,NPOS,NPOS+NSTRIP-1)
        ENDIF
      END SUBROUTINE ORTH1_DISTRI_HERM


!========================== ORTH1_DISTRI_HERM_DESC =====================
!
! same as above, this version requires an additional argument, the descriptor DESCA
! of the matrix
!
!=======================================================================

      SUBROUTINE ORTH1_DISTRI_HERM_DESC(CSEL,CPTWFP,CFW,CPROJ,CPROW, &
           NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,MY_COMM,NBANDSK, DESCA)
        USE prec
        USE scala
        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(*)
        CHARACTER*(*) CSEL
        INTEGER DESCA(*)

        COVL=0
        CALL ORTH1_NOSUBINDEX(CSEL,CPTWFP(1,1),CFW(1,1),CPROJ(1,1), &
           CPROW(1,1),NBANDS, &
           NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL(1,1))

        CALLMPI( M_sum_g(MY_COMM,COVL(1,1),NBANDS*NSTRIP))

        IF (CSEL=='U') THEN
           CALL DISTRI_SLICE_HERM(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCA,NPOS,NPOS+NSTRIP-1)
        ELSE
           CALL DISTRI_SLICE_HERM_L(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCA,NPOS,NPOS+NSTRIP-1)
        ENDIF
      END SUBROUTINE ORTH1_DISTRI_HERM_DESC


!========================== ORTH1_DISTRI_DAVIDSON ======================
!
! ORTH1_DISTRI_DAVIDSON is identical to ORTH1_DISTRI
! but also changes the intermediate COVL matrix as required in
! the davidson.F routine (requires vector R)
!
!=======================================================================

      SUBROUTINE ORTH1_DISTRI_DAVIDSON(CSEL,CPTWFP,CFW,CPROJ,CPROW, &
           NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,R,MY_COMM, NBANDSK )
        USE prec
        USE scala

        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK, I
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(*)
        REAL(q) :: R(NSTRIP)
        CHARACTER*(*) CSEL

        CALL CHECK_scala (NBANDSK, 'ORTH1_DISTRI_DAVIDSON' )

        COVL=0
        CALL ORTH1_NOSUBINDEX(CSEL,CPTWFP(1,1),CFW(1,1),CPROJ(1,1), &
           CPROW(1,1),NBANDS, &
           NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL(1,1))

        COVL(NPOS:NPOS+NSTRIP-1,1:NSTRIP)=0
        IF (MY_COMM%NODE_ME==MY_COMM%IONODE) THEN
           DO I=1,NSTRIP
              COVL(NPOS-1+I,I)=R(I)
           ENDDO
        ENDIF

        CALLMPI( M_sum_g(MY_COMM,COVL(1,1),NBANDS*NSTRIP))

#ifdef ELPA
        IF (LELPA) THEN
           CALL DISTRI_SLICE_HERM(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCSTD,NPOS,NPOS+NSTRIP-1)
        ENDIF
#endif
        IF (.NOT. LELPA) THEN
           CALL DISTRI_SLICE(COVL(1,1),NBANDS,NBANDSK,COVL2(1),DESCSTD,NPOS,NPOS+NSTRIP-1)
        ENDIF

      END SUBROUTINE ORTH1_DISTRI_DAVIDSON

!========================== LINCOM_DISTRI ==============================
!
! LINCOM_DISTRI does essentially the same as LINCOM but
! accepts a distributed matrix COVL2 instead of a global
! NB_TOT times NB_TOT matrix
! the required matrix COVL is reconstructed on the fly in memory
! slice by slice
! the disadvantage is that the transformed wavefunctions need
! to be stored in a temporary array
!
!=======================================================================

      SUBROUTINE LINCOM_DISTRI(MODE,CF,CPROF,COVL,NIN,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN, MY_COMM, NBLK)
        USE prec
        USE scala

        IMPLICIT NONE
        CHARACTER*1 MODE
        TYPE(communic) :: MY_COMM
        INTEGER     NIN, NPL, NPRO, NPLDIM, NPROD, LDTRAN, NBLK
        COMPLEX(q) ::   CF(NPLDIM,NIN),CF_RESULT(NPLDIM,NIN)
        GDEF        CPROF(NPROD,NIN),CPROF_RESULT(NPROD,NIN)
        GDEF        COVL(*)
      ! local
        GDEF        COVL_GLOBAL(LDTRAN,NBLK)
        INTEGER     NPOS, NOUT

DOESI   CF_RESULT=(0.0_q, 0.0_q); CPROF_RESULT=0; COVL_GLOBAL=0
        PROFILING_START('lincom_distri')

        CALL CHECK_scala (NIN, 'LINCOM_DISTRI' )

      ! reconstruct matrix block by block and transform wavefunctions
        DO NPOS=1,NIN-NBLK,NBLK
          NOUT=NBLK
          CALL RECON_SLICE(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCSTD,NPOS,NPOS+NOUT-1)

          PROFILING_START('m_sum_g')
          CALLMPI( M_sum_g(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT))
          PROFILING_STOP('m_sum_g')

          CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))
        ENDDO

        NOUT=NIN-NPOS+1
        CALL RECON_SLICE(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCSTD,NPOS,NPOS+NOUT-1)

        PROFILING_START('m_sum_g')
        CALLMPI( M_sum_g(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT))
        PROFILING_STOP('m_sum_g')

        CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))

        CF   =CF_RESULT
        CPROF=CPROF_RESULT

        PROFILING_STOP('lincom_distri')

      END SUBROUTINE LINCOM_DISTRI

!
! similar version as above but add the Hermitian conjugated elements
! in the lower triangle
! can be used for "full" mode (MODE="F") if the transformation
! is Hermitian and only upper triangle is set (ZPOTRI does that)
!
      SUBROUTINE LINCOM_DISTRI_HERM(MODE,CF,CPROF,COVL,NIN,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN, MY_COMM, NBLK)
        USE prec
        USE scala

        IMPLICIT NONE
        CHARACTER*1 MODE
        TYPE(communic) :: MY_COMM
        INTEGER     NIN, NPL, NPRO, NPLDIM, NPROD, LDTRAN, NBLK
        COMPLEX(q) ::   CF(NPLDIM,NIN),CF_RESULT(NPLDIM,NIN)
        GDEF        CPROF(NPROD,NIN),CPROF_RESULT(NPROD,NIN)
        GDEF        COVL(*)
      ! local
        GDEF        COVL_GLOBAL(LDTRAN,NBLK)
        INTEGER     NPOS, NOUT

        PROFILING_START('lincom_distri_herm')

        CALL CHECK_scala (NIN, 'LINCOM_DISTRI_HERM' )

      ! reconstruct matrix block by block and transform wavefunctions
        DO NPOS=1,NIN-NBLK,NBLK
          NOUT=NBLK
          CALL RECON_SLICE_HERM(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCSTD,NPOS,NPOS+NOUT-1)

          CALLMPI( M_sum_g(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT))

          CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))
        ENDDO

        NOUT=NIN-NPOS+1
        CALL RECON_SLICE_HERM(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCSTD,NPOS,NPOS+NOUT-1)

        CALLMPI( M_sum_g(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT))

        CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))

        CF   =CF_RESULT
        CPROF=CPROF_RESULT

        PROFILING_STOP('lincom_distri_herm')

      END SUBROUTINE LINCOM_DISTRI_HERM


!========================== LINCOM_DISTRI_DESC =========================
!
! LINCOM_DISTRI_DESC essentially the same as LINCOM_DISTRI
! but additionally requires the distributed matrix descriptor
!
!=======================================================================

      SUBROUTINE LINCOM_DISTRI_DESC(MODE,CF,CPROF,COVL,NIN,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN, MY_COMM, NBLK, DESCA)
        USE prec
        USE scala

        IMPLICIT NONE
        CHARACTER*1 MODE
        TYPE(communic) :: MY_COMM
        INTEGER     NIN, NPL, NPRO, NPLDIM, NPROD, LDTRAN, NBLK
        COMPLEX(q) ::   CF(NPLDIM,NIN),CF_RESULT(NPLDIM,NIN)
        GDEF        CPROF(NPROD,NIN),CPROF_RESULT(NPROD,NIN)
        GDEF        COVL(*)
        INTEGER DESCA(*)                               ! distributed matrix descriptor array
      ! local
        GDEF        COVL_GLOBAL(LDTRAN,NBLK)
        INTEGER     NPOS, NOUT

      ! reconstruct matrix block by block and transform wavefunctions
        DO NPOS=1,NIN-NBLK,NBLK
          NOUT=NBLK
          CALL RECON_SLICE(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCA  ,NPOS,NPOS+NOUT-1)

          CALLMPI( M_sum_g(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT))

          CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))
        ENDDO

        NOUT=NIN-NPOS+1
        CALL RECON_SLICE(COVL_GLOBAL(1,1),LDTRAN,NIN,COVL(1),DESCA  ,NPOS,NPOS+NOUT-1)

        CALLMPI( M_sum_g(MY_COMM,COVL_GLOBAL(1,1),LDTRAN*NOUT))

        CALL LINCOM_SLICE(MODE,CF(1,1),CPROF(1,1),COVL_GLOBAL(1,1), &
                NIN,NPOS,NOUT,NPL,NPRO,NPLDIM,NPROD,LDTRAN, &
                CF_RESULT(1,NPOS),CPROF_RESULT(1,NPOS))

        CF   =CF_RESULT
        CPROF=CPROF_RESULT

      END SUBROUTINE LINCOM_DISTRI_DESC

!=======================================================================
!
! routine to call DISTRI_SLICE without explicitly checking
! the interface
!
!=======================================================================
      SUBROUTINE DISTRI_SLICE_NOINT(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)
        USE prec
        USE scala
        IMPLICIT NONE

        INTEGER NDIM                                   ! first dimension of matrix
        INTEGER N                                      ! rank of matrix
        INTEGER COLUMN_LOW, COLUMN_HIGH                ! starting and end row supplied in AMATIN
        GDEF    AMATIN(NDIM,COLUMN_HIGH-COLUMN_LOW+1)  ! input matrix (globally known)
        GDEF    A(*)                                   ! output distributed matrix
        INTEGER DESCA(*)                               ! distributed matrix descriptor array

        CALL DISTRI_SLICE(AMATIN,NDIM,N,A,DESCA,COLUMN_LOW, COLUMN_HIGH)

      END SUBROUTINE DISTRI_SLICE_NOINT


      SUBROUTINE BG_pDSSYEX_ZHEEVX_NOINT(COMM,A,W,N)
      USE prec
      USE scala
      IMPLICIT NONE

      TYPE (communic) COMM
      INTEGER N               ! NxN matrix to be distributed
      GDEF    A(*)            ! input/output matrix
      REAL(q) W(N)            ! eigenvalues

      CALL BG_pDSSYEX_ZHEEVX(COMM,A(1),W,N)

      END SUBROUTINE BG_pDSSYEX_ZHEEVX_NOINT


      SUBROUTINE SCA_SHIFT_BANDS_LOW_NOINT(N, NLOW, A, SHIFT, DESCA )
      USE scala
      IMPLICIT NONE

      INTEGER N              ! actual dimension of the matrix
      INTEGER NLOW           ! first band index not to be shifted
      GDEF    A(*)           ! distributed matrix
      REAL(q) SHIFT          ! magnitude of shift
      INTEGER DESCA(*)       ! distributed matrix descriptor array

      CALL SCA_SHIFT_BANDS_LOW(N, NLOW, A, SHIFT, DESCA )
      END SUBROUTINE

      SUBROUTINE SCA_SHIFT_BANDS_HIGH_NOINT(N, NHIGH, A, SHIFT, DESCA )
      USE scala
      IMPLICIT NONE

      INTEGER N              ! actual dimension of the matrix
      INTEGER NHIGH          ! first band index not to be shifted
      GDEF    A(*)           ! distributed matrix
      REAL(q) SHIFT          ! magnitude of shift
      INTEGER DESCA(*)       ! distributed matrix descriptor array

      CALL SCA_SHIFT_BANDS_HIGH(N, NHIGH, A, SHIFT, DESCA )
      END SUBROUTINE


!=======================================================================
!
! dummy routines if scaLAPACK is not compiled in
!
!=======================================================================
#else
      SUBROUTINE ORTH1_DISTRI(CSEL,CPTWFP,CFW,CPROJ,CPROW,NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,MY_COMM, NBANDSK)
        USE prec
        USE scala

        IMPLICIT NONE

        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        TYPE(communic) :: MY_COMM
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(1,1)
        CHARACTER*(*) CSEL

        WRITE(0,*) "internal error in scala.F:  ORTH1_DISTRI is not available"

      END SUBROUTINE ORTH1_DISTRI

      SUBROUTINE ORTH1_DISTRI_DAVIDSON(CSEL,CPTWFP,CFW,CPROJ,CPROW,NBANDS,NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2,R,MY_COMM, NBANDSK)
        USE prec
        USE scala

        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(1,1)
        REAL(q) :: R(NSTRIP)
        CHARACTER*(*) CSEL

        WRITE(0,*) "internal error in scala.F:  ORTH1_DISTRI_DAVIDSON is not available"

      END SUBROUTINE

      SUBROUTINE ORTH1_DISTRI_DAVMAT(CSEL,CPTWFP,CFW,CPROJ,CPROW,NBANDS, &
           NPOS,NSTRIP,NPL,NPRO,NPLDIM,NPROD,COVL2, &
           CHAM_CORRECTION,MY_COMM, NBANDSK )
        USE prec
        USE scala

        IMPLICIT NONE

        TYPE(communic) :: MY_COMM
        INTEGER :: NPRO, NPL, NSTRIP, NPLDIM, NPOS, NPROD, NBANDS, NBANDSK
        COMPLEX(q) :: CPTWFP(NPLDIM,NBANDS), CFW(NPLDIM,NSTRIP)
        GDEF      CPROJ(NPROD,NBANDS)
        GDEF      CPROW(NPROD,NSTRIP)
        GDEF      COVL(NBANDS,NSTRIP), COVL2(1,1)
        GDEF      CHAM_CORRECTION(NSTRIP,NSTRIP)
        CHARACTER*(*) CSEL

        WRITE(0,*) "internal error in scala.F:  ORTH1_DISTRI_DAVIDSON is not available"

      END SUBROUTINE ORTH1_DISTRI_DAVMAT

      SUBROUTINE LINCOM_DISTRI(MODE,CF,CPROF,COVL,NIN,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN, MY_COMM, NBLK)
        USE prec
        USE scala

        IMPLICIT NONE
        CHARACTER*1 MODE
        TYPE(communic) :: MY_COMM
        INTEGER     NIN, NPL, NPRO, NPLDIM, NPROD, LDTRAN, NBLK
        COMPLEX(q) ::   CF(NPLDIM,NIN),CF_RESULT(NPLDIM,NIN)
        GDEF        CPROF(NPROD,NIN),CPROF_RESULT(NPROD,NIN)
        GDEF        COVL(1,1)

        WRITE(0,*) "internal error in scala.F:  LINCOM_DISTRI is not available"

      END SUBROUTINE LINCOM_DISTRI


      SUBROUTINE LINCOM_DISTRI_HERM(MODE,CF,CPROF,COVL,NIN,NPL, &
     &           NPRO,NPLDIM,NPROD,LDTRAN, MY_COMM, NBLK)
        USE prec
        USE scala

        IMPLICIT NONE
        CHARACTER*1 MODE
        TYPE(communic) :: MY_COMM
        INTEGER     NIN, NPL, NPRO, NPLDIM, NPROD, LDTRAN, NBLK
        COMPLEX(q) ::   CF(NPLDIM,NIN),CF_RESULT(NPLDIM,NIN)
        GDEF        CPROF(NPROD,NIN),CPROF_RESULT(NPROD,NIN)
        GDEF        COVL(1,1)

        WRITE(0,*) "internal error in scala.F:  LINCOM_DISTRI is not available"

      END SUBROUTINE LINCOM_DISTRI_HERM
#endif




