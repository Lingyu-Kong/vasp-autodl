#include "symbol.inc"
!************************ SUBROUTINE STUFAK ****************************
! RCS:  $Id: stufak.F,v 1.1 2000/11/15 08:13:54 kresse Exp $
!
! this subroutine calculates the structure factor on the grid of
! reciprocal lattice vectors
! cstrf(g) = sum over ions (-exp(ig.r)) where r is the position of the
! ion
!***********************************************************************

      SUBROUTINE STUFAK(GRIDC,T_INFO,CSTRF)
#ifdef _OPENACC
      USE mopenacc_struct_def
#endif
      USE prec

      USE mpimy
      USE mgrid
      USE poscar
      USE constant
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (type_info)   T_INFO

      COMPLEX(q) CSTRF(GRIDC%MPLWV,T_INFO%NTYP)

      PROFILING_START('stufak')

! loop over all types of atoms
      NIS=1
      typ: DO NT=1,T_INFO%NTYP
      CALL STUFAK_ONE(GRIDC,T_INFO%NITYP(NT),T_INFO%POSION(1,NIS),T_INFO%VCA(NT),CSTRF(1,NT))
      NIS=NIS+T_INFO%NITYP(NT)

      ENDDO typ

!$ACC UPDATE DEVICE(CSTRF) IF_PRESENT __IF_ASYNC__
      PROFILING_STOP('stufak')

      RETURN
      END


!************************ SUBROUTINE STUFAK_ONE ************************
!
! this subroutine calculates the structure factor on the grid of
! for one species (i.e. partial structure factor)
! cstrf(g) = sum over ions (-exp(ig.r)) where r is the position of the
! ion
!
!***********************************************************************

      SUBROUTINE STUFAK_ONE(GRIDC,NIONS,POSION,VCA,CSTRF)
      USE prec

      USE mpimy
      USE mgrid
      USE poscar
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      REAL(q) POSION(3,NIONS)
      REAL(q) VCA

      TYPE (grid_3d)     GRIDC

      COMPLEX(q) CSTRF(GRIDC%RC%NP)

      CSTRF=0

#ifndef _OPENMP
      ion: DO NI=1,NIONS
!=======================================================================
! loop over all grid points
!=======================================================================
#ifdef vector
!-----------------------------------------------------------------------
! simple version for vector machines
!-----------------------------------------------------------------------
         DO N=1,GRIDC%RC%NP
            N1= MOD((N-1),GRIDC%RC%NROW) +1
            NC= (N-1)/GRIDC%RC%NROW+1
            N2= GRIDC%RC%I2(NC)
            N3= GRIDC%RC%I3(NC)

            G3=POSION(3,NI)*GRIDC%LPCTZ(N3)
            G2=POSION(2,NI)*GRIDC%LPCTY(N2)
            G1=POSION(1,NI)*GRIDC%LPCTX(N1)

            CSTRF(N)=CSTRF(N)+EXP(-CITPI*(G3+G2+G1))*VCA
         ENDDO
#else
!-----------------------------------------------------------------------
! more envolved version which is faster on most (scalar) machines
! and includes support for parallel machines
!-----------------------------------------------------------------------
         CX =EXP(-CITPI*POSION(1,NI))
         G1 =POSION(1,NI)*xmin(GRIDC%NGX)

         col: DO NC=1,GRIDC%RC%NCOL
            N=(NC-1)*GRIDC%RC%NROW+1

            N2= GRIDC%RC%I2(NC)
            N3= GRIDC%RC%I3(NC)
            G2=POSION(2,NI)*GRIDC%LPCTY(N2)
            G3=POSION(3,NI)*GRIDC%LPCTZ(N3)
            CE=EXP(-CITPI*(G3+G2+G1))*VCA
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
            DO N1P=0,GRIDC%RC%NROW-1
               SETN1fromN1P
               CSTRF(N+N1)=CSTRF(N+N1)+CE
               CE=CE*CX
            ENDDO
         ENDDO col
#endif
!-----------------------------------------------------------------------
!  next ion
!-----------------------------------------------------------------------
      ENDDO ion
#else
!-----------------------------------------------------------------------
! OpenMP version: the loop with the largest count (col) has been OMP
! parallelized, and the loop over ions has been moved inside it to
! maximize the amount of work inside the parallel region.
! This comes at the cost of an increase in EXP function calls inside the
! col loop (related to CX), but should yield less cache misses in the
! N1P loop ...
!-----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(N,N2,N3,NI,G1,G2,G3,CX,CE,N1,N1P)
      col: DO NC=1,GRIDC%RC%NCOL
         N=(NC-1)*GRIDC%RC%NROW+1

         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)

         ion:  DO NI=1,NIONS

            G2=POSION(2,NI)*GRIDC%LPCTY(N2)
            G3=POSION(3,NI)*GRIDC%LPCTZ(N3)
            G1=POSION(1,NI)*xmin(GRIDC%NGX)

            CX=EXP(-CITPI*POSION(1,NI))
            CE=EXP(-CITPI*(G3+G2+G1))*VCA
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
            DO N1P=0,GRIDC%RC%NROW-1
              SETN1fromN1P
              CSTRF(N+N1)=CSTRF(N+N1)+CE
              CE=CE*CX
            ENDDO

          ENDDO ion
      ENDDO col
!$OMP END PARALLEL DO
#endif
      RETURN
      END

