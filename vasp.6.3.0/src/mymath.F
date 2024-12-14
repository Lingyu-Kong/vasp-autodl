       MODULE MYMATH
         USE prec
         USE constant
         USE poscar
         IMPLICIT NONE
         INCLUDE "gadget.inc"

       CONTAINS

       FUNCTION coordreallocate(p,m)
         !c little usefull tool taken from Numerical Recipies
         !c reallocates pointer preserving its size
         TYPE(coordinate),DIMENSION(:),POINTER :: p,coordreallocate
         INTEGER,INTENT(IN) :: m
         INTEGER :: mold

         ALLOCATE(coordreallocate(m))
         IF (.NOT. ASSOCIATED(p)) RETURN
         !ALLOCATE(coordreallocate(m))
         mold=SIZE(p)
         mold=min(mold,m)
         coordreallocate(1:mold)=p(1:mold)
         DEALLOCATE(p)
       END FUNCTION

       FUNCTION qreallocate(p,m,n)
         !c little usefull tool taken from Numerical Recipies
         !c reallocates pointer preserving its size
         REAL(q),DIMENSION(:,:),POINTER :: p,qreallocate
         INTEGER :: shp(2)
         INTEGER,INTENT(IN) :: m,n
         INTEGER :: mold,nold

         ALLOCATE(qreallocate(m,n))
         IF (.NOT. ASSOCIATED(p)) RETURN
         shp=SHAPE(p)
         mold=shp(1)
         nold=shp(2)
         mold=min(mold,m)
         nold=min(nold,n)
         qreallocate(1:mold,1:nold)=p(1:mold,1:nold)
         DEALLOCATE(p)
       END FUNCTION

       FUNCTION ireallocate(p,m,n)
         !c little usefull tool taken from Numerical Recipies
         !c reallocates pointer preserving its size
         INTEGER,DIMENSION(:,:),POINTER :: p,ireallocate
         INTEGER :: shp(2)
         INTEGER,INTENT(IN) :: m,n
         INTEGER :: mold,nold

         ALLOCATE(ireallocate(m,n))
         IF (.NOT. ASSOCIATED(p)) RETURN
         shp=SHAPE(p)
         mold=shp(1)
         nold=shp(2)
         mold=min(mold,m)
         nold=min(nold,n)
         ireallocate(1:mold,1:nold)=p(1:mold,1:nold)
         DEALLOCATE(p)
       END FUNCTION

       FUNCTION OUTERPROD(NDIM,avect,bvect)
       !c returns outerproduct of two vectors
         INTEGER :: NDIM
         INTEGER :: i,j
         REAL(q) :: avect(NDIM)
         REAL(q) :: bvect(NDIM)
         REAL(q) :: OUTERPROD(NDIM,NDIM)

         DO i=1,NDIM
           DO j=1,NDIM
             OUTERPROD(i,j)=avect(i)*bvect(j)
           ENDDO
         ENDDO
       END FUNCTION

       FUNCTION OUTERPROD2(NDIM,avect,bvect)
       !c returns outerproduct of two vectors
         INTEGER :: NDIM
         INTEGER :: avect(NDIM)
         INTEGER :: i,j
         REAL(q) :: bvect(NDIM)
         REAL(q) :: OUTERPROD2(NDIM,NDIM)

         DO i=1,NDIM
           DO j=1,NDIM
             OUTERPROD2(i,j)=avect(i)*bvect(j)
           ENDDO
         ENDDO
       END FUNCTION


       FUNCTION INNERPROD(NDIM,avect,bvect)
       !c returns innerproduct of two vectors
         INTEGER :: NDIM
         REAL(q) :: avect(NDIM)
         REAL(q) :: bvect(NDIM)
         REAL(q) :: INNERPROD
         
         INNERPROD=SUM(avect*bvect)
       END FUNCTION
       
       FUNCTION MATVECTMUL(NDIM,mat,vect)
         INTEGER :: i,NDIM
         REAL(q) :: mat(NDIM,NDIM)
         REAL(q) :: vect(NDIM)
         REAL(q) :: MATVECTMUL(NDIM)

         DO i=1,NDIM
           MATVECTMUL(i)=SUM(mat(i,:)*vect)
         ENDDO
       END FUNCTION

       FUNCTION VECTORSIZE(NDIM,VECT)
         INTEGER :: NDIM
         REAL(q) :: VECT(NDIM)
         REAL(q) :: SVECT(NDIM)
         REAL(q) :: VECTORSIZE

         SVECT=VECT*VECT
         VECTORSIZE=SUM(SVECT)
         VECTORSIZE=VECTORSIZE**(0.5)
       END FUNCTION

       FUNCTION CROSSPROD(NDIM,A,B)
         INTEGER :: NDIM
         REAL(q) :: A(NDIM),B(NDIM)
         REAL(q) :: X,Y,Z
         REAL(q) :: CROSSPROD(NDIM)

         X=A(2)*B(3)-B(2)*A(3)
         Y=A(3)*B(1)-B(3)*A(1)
         Z=A(1)*B(2)-B(1)*A(2)
         CROSSPROD=(/X,Y,Z/)
       END FUNCTION

       SUBROUTINE SVDVALVEC(imat,dim,d,u,vt)
         INTEGER :: dim
         REAL(q) :: imat(dim,dim)
         REAL(q) :: d(dim)         !c eigenvalues of the matrix
         REAL(q) :: e(dim-1)
         REAL(q) :: tauq(dim)
         REAL(q) :: taup(dim)
         !REAL(q),ALLOCATABLE :: work(:)
         REAL(q) :: work(4*dim)
         REAL(q) :: vt(dim,dim)
         REAL(q) :: u(dim,dim)     !c eigenvectors of the matrix
         REAL(q) :: c(1,1)         
         INTEGER :: lwork
         INTEGER :: info
         INTEGER :: n
         EXTERNAL :: DGEBRD
         EXTERNAL :: DORGBR
         EXTERNAL :: DBDSQR

         lwork=4*dim
         !ALLOCATE(work(lwork))
        
         CALL DGEBRD(dim,dim,imat,dim,d,e,tauq,taup,work,lwork,info)
         vt=imat
         u=imat
         CALL DORGBR('P',dim,dim,dim,vt,dim,taup,work,lwork,info)
         CALL DORGBR('Q',dim,dim,dim,u,dim,tauq,work,lwork,info)
         CALL DBDSQR('U',dim,dim,dim,0,d,e,vt,dim,u,dim,c,dim,work,info)
         !DEALLOCATE(work)
         
       END SUBROUTINE
        
       SUBROUTINE EIGVAL_GENERAL(imat,dim,d)
       !c computes eigenvalues of general square matrix
         INTEGER :: dim
         REAL(q),INTENT(in) :: imat(dim,dim)
         REAL(q),allocatable :: jmat(:,:)
         REAL(q) :: d(dim)         !c eigenvalues of the matrix
         REAL(q) :: di(dim)        !c eigenvalues of the matrix, imaginary part
         REAL(q) :: work(4*dim)
         REAL(q) :: vt(dim,dim)
         REAL(q) :: u(dim,dim)     !c eigenvectors of the matrix 
         INTEGER :: lwork 
         INTEGER :: info
         EXTERNAL :: DGEEV
        
         allocate(jmat(dim,dim))
         jmat=imat 
         lwork=4*dim 
         !c dgeev overwrites the inputed matrix!!!
         jmat=imat
         CALL DGEEV( 'N', 'N', dim, jmat, dim, d, di, u, dim, vt, &
            &                  dim, WORK, LWORK, INFO )
         deallocate(jmat) 
       END SUBROUTINE EIGVAL_GENERAL

       SUBROUTINE ZEIGVAL_GENERAL(imat,N,w)
       !c computes eigenvalues of general square matrix
         INTEGER :: N,lda
         INTEGER :: lwork,ldvl,ldvr
         INTEGER :: info
         COMPLEX(q),INTENT(in) :: imat(N,N)
         COMPLEX(q),allocatable :: A(:,:)
         COMPLEX(q) :: w(N)      !c eigenvalues of the matrix
         COMPLEX(q):: vl(N,N),vr(N,N)
         COMPLEX(q):: work(4*N)
         REAL(q):: RWORK(2*N)
         EXTERNAL :: ZGEEV

         lda=N
         lwork=4*N
         ldvl=N;ldvr=N
         allocate(A(N,N))
         !w=cmplx(0.)
         !vl=cmplx(0.);vr=cmplx(0.)
         !allocate(A(lda,N))
         !allocate(work(lwork))
         !allocate(vl(lvdl,N));allocate(vr(lvdr,N))
         !allocate(RWORK(2*N))
         !c dgeev overwrites the inputed matrix!!!
         A=imat
         CALL ZGEEV( 'N', 'N', N, A, lda, w, vl, ldvl, vr, &
            &                  ldvr, WORK, LWORK,RWORK, INFO )
         !CALL CGEEV( 'N', 'N', N, A, lda, w,  ldvl,  &
         !   &                  ldvr, WORK, LWORK,RWORK, INFO )
         !deallocate(vr,vl,work,A,rwork)
         deallocate(A)
         !write(*,*) 'INFO',INFO
       END SUBROUTINE ZEIGVAL_GENERAL

       SUBROUTINE ZEIGVAL_GENERAL2(imat,N,W)
         USE string, ONLY: str
         USE tutor, ONLY: vtutor
         INTEGER :: i
         INTEGER :: N
         INTEGER :: LWORK,IWORK,LRWORK
         INTEGER :: info
         COMPLEX(q),ALLOCATABLE :: A(:,:)
         COMPLEX(q),ALLOCATABLE:: work(:),RWORK(:)
         EXTERNAL :: ZGESDD
         REAL(q) :: qdummy
         COMPLEX(q) :: imat(N,N)
         REAL(q) :: S(N)
         COMPLEX(q) :: U(N,N),VT(N,N)
         COMPLEX(q) :: W(N)

         LWORK=2*N*N+2*N+N
         LRWORK=N*(5*N+7)
         IWORK=8*N
         ALLOCATE(WORK(LWORK),RWORK(LRWORK))
         allocate(A(N,N))
         A=imat
         WORK=0.;RWORK=0.

         CALL ZGESDD('N',N,N,A,N,S,U,N,VT,N,WORK,LWORK,RWORK,IWORK,INFO)

         IF (INFO==0) THEN
           DO i=1,N
             w(i)=CMPLX(S(i),0._q,q)
           ENDDO
         ELSE
           CALL vtutor%error("ZEIGVAL_GENERAL2(ZGESDD): " // str(INFO))
         ENDIF
         DEALLOCATE(A,RWORK,WORK)
       END SUBROUTINE   ZEIGVAL_GENERAL2

       SUBROUTINE INVERSE_Z(A,N)
         USE string, ONLY: str
         USE tutor, ONLY: vtutor
         INTEGER :: N,lda,M
         INTEGER :: lwork,ldvl,ldvr
         INTEGER :: info
         COMPLEX(q) :: A(N,N)
         INTEGER :: IPIV(N)
         COMPLEX(q):: work(4*N)
         EXTERNAL :: ZGETRF,ZGETRI

         LWORK=4*N !c this should be optimized using ILAENV
         M=N
         LDA=N

         CALL ZGETRF( M, N, A, LDA, IPIV, INFO )

         IF (INFO==0) THEN
           CALL ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
           IF (INFO .NE. 0) THEN
             CALL vtutor%error("INVERSE_Z(ZGETRI): " // str(INFO))
           ENDIF
         ELSE
           CALL vtutor%error("INVERSE_Z(ZGETRF): " // str(INFO))
         ENDIF 
       END SUBROUTINE INVERSE_Z

       SUBROUTINE SVD_INVERSE_Z(A,N)
         !USE c2f_interface, ONLY : ERRF
         USE string, ONLY: str
         USE tutor, ONLY: vtutor
         INTEGER :: i,j
         INTEGER :: N,LDA
         INTEGER :: lwork,ldvl,ldvr,IWORK,LRWORK
         INTEGER :: info
         COMPLEX(q) :: A(N,N)
         COMPLEX(q),ALLOCATABLE:: work(:),RWORK(:)
         EXTERNAL :: ZGEMM,ZGESDD
         REAL(q),PARAMETER :: tiny_my=1e-6
         REAL(q) :: qdummy
         COMPLEX(q) :: imat(N,N)
         COMPLEX(q) :: zdumm0=CMPLX(0._q,0._q,q),zdumm1=CMPLX(1._q,0._q,q)
         REAL(q) :: S(N)
         COMPLEX(q) :: U(N,N),VT(N,N)
         !REAL(q), EXTERNAL :: ERRF

         LWORK=2*N*N+2*N+N
         LRWORK=N*(5*N+7)
         IWORK=8*N
         ALLOCATE(WORK(LWORK),RWORK(LRWORK))
         WORK=0.;RWORK=0.
         LDA=N

         CALL ZGESDD('A',N,N,A,N,S,U,N,VT,N,WORK,LWORK,RWORK,IWORK,INFO)

         IF (INFO==0) THEN
           A=CMPLX(0._q,0._q,q)
           DO i=1,N
             qdummy=S(i)
             IF (ABS(qdummy)>tiny_my) THEN
               U(:,i)=U(:,i)/qdummy
             ELSE
               U(:,i)=0.
             ENDIF
           ENDDO

           CALL ZGEMM('N','N',N,N,N,zdumm1,U,N,VT,N, &
                &    zdumm0,A,N)

         ELSE
           CALL vtutor%error("INVERSE2_Z(ZGESDD): " // str(INFO))
         ENDIF
         DEALLOCATE(RWORK,WORK)
       END SUBROUTINE SVD_INVERSE_Z

       SUBROUTINE SVD_INVERSE_POSITIVE_Z(A,N,tiny_my)
         !USE c2f_interface, ONLY : ERRF
         USE string, ONLY: str
         USE tutor, ONLY: vtutor
         INTEGER :: i,j
         INTEGER :: N,LDA
         INTEGER :: lwork,ldvl,ldvr,IWORK,LRWORK
         INTEGER :: info
         COMPLEX(q) :: A(N,N)
         COMPLEX(q),ALLOCATABLE:: work(:),RWORK(:)
         EXTERNAL :: ZGEMM,ZGESDD
         REAL(q) :: tiny_my
         REAL(q) :: qdummy
         COMPLEX(q) :: imat(N,N)
         COMPLEX(q) :: zdumm0=CMPLX(0._q,0._q,q),zdumm1=CMPLX(1._q,0._q,q)
         REAL(q) :: S(N)
         COMPLEX(q) :: U(N,N),VT(N,N)
         !REAL(q), EXTERNAL :: ERRF

         LWORK=2*N*N+2*N+N
         LRWORK=N*(5*N+7)
         IWORK=8*N
         ALLOCATE(WORK(LWORK),RWORK(LRWORK))
         WORK=0.;RWORK=0.
         LDA=N

         CALL ZGESDD('A',N,N,A,N,S,U,N,VT,N,WORK,LWORK,RWORK,IWORK,INFO)

         IF (INFO==0) THEN
           A=CMPLX(0._q,0._q,q)
           DO i=1,N
             qdummy=S(i)
             IF (qdummy>tiny_my) THEN
               U(:,i)=U(:,i)/qdummy
             ELSE
               U(:,i)=0.
             ENDIF
           ENDDO

           !CALL ZGEMM('N','N',N,N,N,zdumm1,U,N,VT,N, &
           !     &    zdumm0,A,N)
           CALL ZGEMM('C','C',N,N,N,zdumm1,VT,N,U ,N, &
                &    zdumm0,A,N)

         ELSE
           CALL vtutor%error("INVERSE2_Z(ZGESDD): " // str(INFO))
         ENDIF

         A=transpose(A)
         A=CONJG(A)

         DEALLOCATE(RWORK,WORK)
       END SUBROUTINE SVD_INVERSE_POSITIVE_Z

       SUBROUTINE SVD_INVERSE_AND_SHIFT_Z(A,N,P)
         USE c2f_interface, ONLY : ERRF
         USE string, ONLY: str
         USE tutor, ONLY: vtutor
         INTEGER :: i,j
         INTEGER :: N,LDA
         INTEGER :: lwork,ldvl,ldvr,IWORK,LRWORK
         INTEGER :: info
         COMPLEX(q) :: A(N,N)
         COMPLEX(q),ALLOCATABLE:: work(:),RWORK(:)
         EXTERNAL :: ZGEMM,ZGESDD,ZGESVD
         REAL(q) :: P
         REAL(q) :: OneOverP
         REAL(q),PARAMETER :: tiny_my=1e-6
         REAL(q),PARAMETER :: SQRTPi=SQRT(Pi)
         REAL(q) :: qdummy
         COMPLEX(q) :: imat(N,N)
         COMPLEX(q) :: zdumm0=CMPLX(0._q,0._q,q),zdumm1=CMPLX(1._q,0._q,q)
         REAL(q) :: S(N)
         COMPLEX(q) :: U(N,N),VT(N,N)
         !REAL(q), EXTERNAL :: ERRF
         LOGICAL :: LZGESDD=.FALSE. !c SDD version is usually faster but crashes 
                                    !c for some dimensions (bug in LAPACK) 

         IF (ABS(P)>0._q) THEN
           OneOverP=1._q/P
         ELSE
            CALL vtutor%error("Error(SVD_INVERSE_AND_SHIFT_Z): P=0")
         ENDIF

         !c minimal dimensions as required by ZGESDD 
         !LWORK=2*N*N+2*N+N
         LWORK=N*(2*N+3)
         LRWORK=N*(5*N+7)
         IWORK=8*N

         !c minimal dimensions as required by ZGESVD
         !LWORK=3*N
         !LRWORK=5*N
         !c are always smaller than those for ZGESDD

         ALLOCATE(WORK(LWORK),RWORK(LRWORK))
         WORK=0._q;RWORK=0._q
         LDA=N

         IF (LZGESDD) THEN
           CALL ZGESDD('A',N,N,A,N,S,U,N,VT,N,WORK,LWORK,RWORK,IWORK,INFO)
         ELSE
           CALL ZGESVD('A','A',N,N,A,N,S,U,N,VT,N,WORK,LWORK,RWORK,INFO)
         END IF


         IF (INFO==0) THEN
           A=CMPLX(0._q,0._q,q)
           DO i=1,N
             !qdummy=S(i)

             qdummy=S(i)-1._q
             IF (qdummy<0.) THEN
               !qdummy=0.5*((-qdummy)**4)*pi**0.5
               !qdummy=-(ERRF(qdummy))**0.25
               qdummy=0.5*((-qdummy)**P)*SQRTPi
               qdummy=-(ERRF(qdummy))**OneOverP
             ENDIF
             qdummy=qdummy+1._q

             IF (ABS(qdummy)>tiny_my) THEN
               U(:,i)=U(:,i)/qdummy
             ELSE
               U(:,i)=CMPLX(0._q,0._q,q)
             ENDIF
           ENDDO

           !CALL ZGEMM('N','N',N,N,N,zdumm1,U,N,VT,N, &
           !     &    zdumm0,A,N)
           CALL ZGEMM('C','C',N,N,N,zdumm1,VT,N,U ,N, &
                &    zdumm0,A,N)


         ELSE
           CALL vtutor%error("INVERSE2_Z(ZGESDD): " // str(INFO))
         ENDIF
         DEALLOCATE(RWORK,WORK)
       END SUBROUTINE SVD_INVERSE_AND_SHIFT_Z

       SUBROUTINE SVD_TIM_reconstruct(AT,T,A,N)
         USE c2f_interface, ONLY : ERRF
         USE string, ONLY: str
         USE tutor, ONLY: vtutor
         INTEGER :: i,j
         INTEGER :: N,LDA
         INTEGER :: lwork,ldvl,ldvr,IWORK,LRWORK
         INTEGER :: info
         COMPLEX(q) :: AT(N,N),T(N,N),A(N,N)
         COMPLEX(q),ALLOCATABLE:: work(:),RWORK(:)
         EXTERNAL :: ZGEMM,ZGESDD
         REAL(q),PARAMETER :: tiny_my=1e-6
         REAL(q) :: qdummy
         COMPLEX(q) :: imat(N,N)
         COMPLEX(q) :: zdumm0=CMPLX(0._q,0._q,q),zdumm1=CMPLX(1._q,0._q,q),zdummM1=CMPLX(-1._q,0._q,q)
         REAL(q) :: S(N)
         COMPLEX(q) :: U(N,N),VT(N,N)
         !REAL(q), EXTERNAL :: ERRF

         LWORK=2*N*N+2*N+N
         LRWORK=N*(5*N+7)
         IWORK=8*N
         ALLOCATE(WORK(LWORK),RWORK(LRWORK))
         WORK=0.;RWORK=0.
         LDA=N

         CALL ZGESDD('A',N,N,AT,N,S,U,N,VT,N,WORK,LWORK,RWORK,IWORK,INFO)

         IF (INFO==0) THEN
           AT=CMPLX(0._q,0._q,q)
           A=CMPLX(0._q,0._q,q)
           DO i=1,N
             !qdummy=S(i)

             qdummy=S(i)-1._q
             IF (qdummy<0.) THEN
               qdummy=0.5*((-qdummy)**4)*pi**0.5
               qdummy=-(ERRF(qdummy))**0.25
             ENDIF
             qdummy=qdummy+1._q

             !A(:,i)=U(:,i)*(S(i)-1._q)

             IF (ABS(qdummy)>tiny_my) THEN
               AT(:,i)=U(:,i)/qdummy
               !T(:,i)=U(:,i)*(qdummy-1._q)/REAL(alpha(i,i))
               A(:,i)=U(:,i)*(qdummy-1._q)
             ELSE
               AT(:,i)=0.
               A(:,i)=0.
             ENDIF
           ENDDO

           U=AT
           AT=0.
           CALL ZGEMM('N','N',N,N,N,zdumm1,U,N,VT,N, &
                &    zdumm0,AT,N)

           U=A
           A=0.
           CALL ZGEMM('N','N',N,N,N,zdumm1,U,N,VT,N, &
                &    zdumm0,A,N)

           VT=T
           CALL INVERSE_Z(VT,N)
           !CALL SVD_INVERSE_Z(VT,N)

           U=A
           A=0.
           CALL ZGEMM('N','N',N,N,N,zdummM1,U,N,VT,N, &
                &    zdumm0,A,N)

         ELSE
           CALL vtutor%error("INVERSE2_Z(ZGESDD): " // str(INFO))
         ENDIF
         DEALLOCATE(RWORK,WORK)
       END SUBROUTINE SVD_TIM_reconstruct

       SUBROUTINE INVERSE_SYM_D(A,N)
         USE string, ONLY: str
         USE tutor, ONLY: vtutor
         INTEGER :: N,LDA,i,j
         INTEGER :: info
         REAL(q) :: A(N,N)
         CHARACTER :: UPLO
         EXTERNAL :: DPOTRF,DPOTRI 
       
         UPLO='U'
         LDA=N

         !c factorization
         CALL DPOTRF(UPLO,N,A,LDA,INFO)
         IF (INFO.EQ.0) THEN
           CALL DPOTRI(UPLO,N,A,LDA,INFO)
           IF (INFO .NE. 0) THEN
             CALL vtutor%error("INVERSE_SYM_D(DPOTRI): " // str(INFO))
           ENDIF
         ELSE
           CALL vtutor%error("INVERSE_SYM_D(DPOTRF): " // str(INFO))
         ENDIF

         !A=A+TRANSPOSE(A)
         DO i=1,N
           DO j=i+1,N
             A(j,i)=A(i,j)
           ENDDO
         ENDDO
       END SUBROUTINE INVERSE_SYM_D



       SUBROUTINE SVDINVERSE(imat,dim,ninfo)
       !c generalized inverse of a square matrix
         INTEGER :: dim
         REAL(q) :: imat(dim,dim)
         REAL(q) :: d(dim)         !c eigenvalues of the matrix
         REAL(q) :: u(dim,dim)     !c eigenvectors of the matrix
         REAL(q) :: vt(dim,dim)
         REAL(q) :: maxdd
         INTEGER :: ninfo
         INTEGER :: n


         CALL SVDVALVEC(imat,dim,d,u,vt) !c calculates evectors and evalues
         maxdd=MAXVAL(ABS(d))
         ninfo=0
         imat=0_q
         DO n=1,dim
           IF (ABS(d(n)/maxdd)>1e-04) THEN
             imat(n,n)=1/d(n)
           ELSE
             ninfo=ninfo+1
           ENDIF
         ENDDO
         imat=MATMUL(imat,TRANSPOSE(u))
         imat=MATMUL(TRANSPOSE(vt),imat)
       END SUBROUTINE


       FUNCTION GIVE_RANK(imat,dim) RESULT(ninfo)
       !c generalized inverse of a square matrix
         INTEGER :: dim
         REAL(q) :: imat(dim,dim)
         REAL(q) :: d(dim)         !c eigenvalues of the matrix
         REAL(q) :: u(dim,dim)     !c eigenvectors of the matrix
         REAL(q) :: vt(dim,dim)
         REAL(q) :: maxdd
         INTEGER :: ninfo
         INTEGER :: n


         CALL SVDVALVEC(imat,dim,d,u,vt) !c calculates evectors and evalues
         maxdd=MAXVAL(ABS(d))
         ninfo=0
         DO n=1,dim
           IF (ABS(d(n)/maxdd)<1e-04) THEN
             ninfo=ninfo+1
           ENDIF
         ENDDO
       END FUNCTION

       SUBROUTINE THREETOONE(N,TVECTOR,M,OVECTOR)
       !c transforms three-column format into
       !c one-column vector
         INTEGER :: i,j,N,M         !M>=N!!!
         REAL(q) :: TVECTOR(3,N)    ! vector in three-column format
         REAL(q) :: OVECTOR(M)   ! vector in one-column format

         OVECTOR=0
         DO i=1,N
           DO j=1,3
             OVECTOR(3*i+j-3)=TVECTOR(j,i)
           ENDDO
         ENDDO
       END SUBROUTINE
       
       SUBROUTINE ONETOTHREE(N,TVECTOR,OVECTOR)
       !c transforms one-column format into
       !c three-column vector
         INTEGER :: i,j,N
         REAL(q) :: TVECTOR(3,N)    ! vector in three-column format
         REAL(q) :: OVECTOR(3*N)   ! vector in one-column format
         
         DO i=1,N
           DO j=1,3
             TVECTOR(j,i)=OVECTOR(3*i+j-3)
           ENDDO
         ENDDO
       END SUBROUTINE

       SUBROUTINE Min_Image_Cyclic(prims1,prims2)
         REAL(q) :: prims1,prims2

         DO
           IF (prims2-prims1>PI) THEN
             prims2=prims2-2*PI
           ELSE
             EXIT
           ENDIF
         ENDDO
         DO
           IF (prims2-prims1<=-PI) THEN
             prims2=prims2+2*PI
           ELSE
             EXIT
           ENDIF
         ENDDO
       END SUBROUTINE Min_Image_Cyclic
       
       SUBROUTINE DECYCLE(prims1,prims2,ICOORDINATES)
         TYPE(coordstructure) :: ICOORDINATES
         REAL(q),DIMENSION(:) :: prims1,prims2
         INTEGER :: i

         DO i=1,ICOORDINATES%NUMINTERNALS
           IF (ICOORDINATES%COORDSTRUCT(i)%TAG=='A ' .OR. &
               ICOORDINATES%COORDSTRUCT(i)%TAG=='T ') THEN
             DO
               IF ((prims2(i)-prims1(i))>PI) THEN
                 prims2(i)=prims2(i)-2*PI
               ELSE
                 EXIT
               ENDIF
             ENDDO
             DO
               IF ((prims2(i)-prims1(i))<=-PI) THEN
                 prims2(i)=prims2(i)+2*PI
               ELSE
                 EXIT
               ENDIF
             ENDDO
           ENDIF
         ENDDO
       END SUBROUTINE DECYCLE

       SUBROUTINE REPAIR_TORSIONS(COORDINATES)
         TYPE(coordstructure) :: COORDINATES
         INTEGER :: i

         DO i=1,COORDINATES%NUMINTERNALS
           IF (COORDINATES%COORDSTRUCT(i)%TAG=='T ') THEN
             DO
               IF (COORDINATES%COORDSTRUCT(i)%VALUE>PI) THEN
                 COORDINATES%COORDSTRUCT(i)%VALUE=COORDINATES%COORDSTRUCT(i)%VALUE-2*PI
               ELSE
                 EXIT
               ENDIF
             ENDDO
             DO
               IF (COORDINATES%COORDSTRUCT(i)%VALUE<=-PI) THEN
                 COORDINATES%COORDSTRUCT(i)%VALUE=COORDINATES%COORDSTRUCT(i)%VALUE+2*PI
               ELSE
                 EXIT
               ENDIF
             ENDDO
           ENDIF
         ENDDO
       END SUBROUTINE REPAIR_TORSIONS

       SUBROUTINE NORMALIZE_MATROW(mat)
         REAL(q),DIMENSION(:,:) :: mat
         REAL(q) :: norm
         INTEGER :: i,dim1,dim2

         dim1=SIZE(mat(:,1))
         dim2=SIZE(mat(1,:))
         DO i=1,dim1
           norm=VECTORSIZE(dim2,mat(i,:))
           IF (norm > 1e-5) THEN
             mat(i,:)=mat(i,:)/norm
           ELSE
             mat(i,:)=0
           ENDIF
         END DO
       END SUBROUTINE NORMALIZE_MATROW

       SUBROUTINE ORTHO_NORMALIZE(mat)
         REAL(q),DIMENSION(:,:) :: mat
         REAL(q) :: norm
         INTEGER :: i,j,k,dim1,dim2

         dim1=SIZE(mat(:,1))
         dim2=SIZE(mat(1,:))
         CALL NORMALIZE_MATROW(mat)
         DO i=2,dim1
           DO j=1,i-1
             mat(i,:)=mat(i,:)-SUM(mat(i,:)*mat(j,:))*mat(j,:)
             norm=VECTORSIZE(dim2,mat(i,:))
             IF (norm>1e-5) THEN
               mat(i,:)=mat(i,:)/norm
             ELSE
               mat(i,:)=0
             END IF
           ENDDO
         ENDDO
       END SUBROUTINE ORTHO_NORMALIZE

       SUBROUTINE iconst_bmat(mat,ICOORDINATES,T_INFO)
         TYPE(coordstructure) :: ICOORDINATES
         TYPE(type_info) :: T_INFO
         REAL(q),DIMENSION(:,:) :: mat
         REAL(q),ALLOCATABLE :: bmat(:,:),cmat(:,:),dmat(:,:)
         REAL(q) :: norm
         INTEGER :: i,j,k,dim1,dim2

         dim1=SIZE(mat(:,1))
         dim2=SIZE(mat(1,:))
         DO i=1,T_INFO%NIONS
           DO j=1,3
             IF ( .NOT. T_INFO%LSFOR(j,i)) mat(:,3*(i-1)+j)=0
           ENDDO
         ENDDO

         ALLOCATE(bmat(dim1,dim2),cmat(dim1,dim2),dmat(dim1,dim1))
         cmat=0._q
         DO i=1,dim1
           !IF (ICOORDINATES%COORDSTRUCT(i)%STATUS/=1) THEN
           IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0) THEN
             cmat(i,:)=mat(i,:)
           ENDIF
         ENDDO
         CALL ORTHO_NORMALIZE(cmat)
         dmat=MATMUL(mat,TRANSPOSE(cmat))
         DO i=1,dim1
           dmat(i,i)=0._q
         ENDDO
         bmat=MATMUL(dmat,cmat)
         mat=mat-bmat
         DEALLOCATE(dmat,cmat,bmat)
       END SUBROUTINE iconst_bmat

       SUBROUTINE where_shortest(DPOS,lattmat,COORDSTRUCT)
       !c determines the minimal image translations for the
       !c internal coordinate
         TYPE(coordinate) :: COORDSTRUCT
         REAL(q) :: DPOS(:,:)            !c positions in direct
         REAL(q) :: lattmat(3,3)         !c direct and reciprocal lattice

         COORDSTRUCT%WHERE=0
         SELECT CASE(COORDSTRUCT%TAG)
           CASE ('R ')
             COORDSTRUCT%WHERE(1,:)=(/0,0,0/)
             COORDSTRUCT%WHERE(2,:)=shortest_dist(DPOS,lattmat,COORDSTRUCT%WHAT(1),COORDSTRUCT%WHAT(2))
           CASE ('W ')
             COORDSTRUCT%WHERE(1,:)=(/0,0,0/)
             COORDSTRUCT%WHERE(2,:)=shortest_dist(DPOS,lattmat,COORDSTRUCT%WHAT(1),COORDSTRUCT%WHAT(2))
           CASE ('WP')
             COORDSTRUCT%WHERE(1,:)=(/0,0,0/)
             COORDSTRUCT%WHERE(2,:)=shortest_dist(DPOS,lattmat,COORDSTRUCT%WHAT(1),COORDSTRUCT%WHAT(2))
           CASE ('Q ')
             COORDSTRUCT%WHERE(1,:)=(/0,0,0/)
             COORDSTRUCT%WHERE(2,:)=shortest_dist(DPOS,lattmat,COORDSTRUCT%WHAT(1),COORDSTRUCT%WHAT(2))
           CASE ('A ')
             COORDSTRUCT%WHERE(1,:)=shortest_dist(DPOS,lattmat,COORDSTRUCT%WHAT(2),COORDSTRUCT%WHAT(1))
             COORDSTRUCT%WHERE(2,:)=(/0,0,0/)
             COORDSTRUCT%WHERE(3,:)=shortest_dist(DPOS,lattmat,COORDSTRUCT%WHAT(2),COORDSTRUCT%WHAT(3))
           CASE ('T ')
             COORDSTRUCT%WHERE(1,:)=shortest_dist(DPOS,lattmat,COORDSTRUCT%WHAT(2),COORDSTRUCT%WHAT(1))
             COORDSTRUCT%WHERE(2,:)=(/0,0,0/)
             COORDSTRUCT%WHERE(3,:)=shortest_dist(DPOS,lattmat,COORDSTRUCT%WHAT(2),COORDSTRUCT%WHAT(3))
             COORDSTRUCT%WHERE(4,:)=shortest_dist(DPOS,lattmat,COORDSTRUCT%WHAT(3),COORDSTRUCT%WHAT(4))+&
             &COORDSTRUCT%WHERE(3,:)
           CASE ('M ')
             COORDSTRUCT%WHERE(1,:)=(/0,0,0/)
             COORDSTRUCT%WHERE(2,:)=shortest_dist(DPOS,lattmat,COORDSTRUCT%WHAT(1),COORDSTRUCT%WHAT(2))
             COORDSTRUCT%WHERE(3,:)=shortest_dist(DPOS,lattmat,COORDSTRUCT%WHAT(2),COORDSTRUCT%WHAT(3))
             COORDSTRUCT%WHERE(3,:)=COORDSTRUCT%WHERE(3,:)+COORDSTRUCT%WHERE(2,:)
           CASE ('B ')  !c should be done better - in fact the nearest midpoint-midpoint distance is desired
             COORDSTRUCT%WHERE(1,:)=(/0,0,0/)
             COORDSTRUCT%WHERE(2,:)=shortest_dist(DPOS,lattmat,COORDSTRUCT%WHAT(1),COORDSTRUCT%WHAT(2))
             COORDSTRUCT%WHERE(3,:)=shortest_dist(DPOS,lattmat,COORDSTRUCT%WHAT(1),COORDSTRUCT%WHAT(3))
             COORDSTRUCT%WHERE(4,:)=shortest_dist(DPOS,lattmat,COORDSTRUCT%WHAT(3),COORDSTRUCT%WHAT(4))
             COORDSTRUCT%WHERE(4,:)=COORDSTRUCT%WHERE(4,:)+COORDSTRUCT%WHERE(3,:)
           CASE ('P ')
             COORDSTRUCT%WHERE(1,:)=shortest_dist(DPOS,lattmat,COORDSTRUCT%WHAT(2),COORDSTRUCT%WHAT(1))
             COORDSTRUCT%WHERE(2,:)=(/0,0,0/)
             COORDSTRUCT%WHERE(3,:)=shortest_dist(DPOS,lattmat,COORDSTRUCT%WHAT(2),COORDSTRUCT%WHAT(3))
             COORDSTRUCT%WHERE(4,:)=shortest_dist(DPOS,lattmat,COORDSTRUCT%WHAT(3),COORDSTRUCT%WHAT(4))+&
             &COORDSTRUCT%WHERE(3,:)
           CASE('E ','F ','G ')
             COORDSTRUCT%WHERE(1,:)=(/0,0,0/);COORDSTRUCT%WHERE(2,:)=(/0,0,0/);
             COORDSTRUCT%WHERE(3,:)=(/0,0,0/);COORDSTRUCT%WHERE(4,:)=(/0,0,0/);
         END SELECT
       END SUBROUTINE where_shortest

       FUNCTION shortest_dist(DPOS,lattmat,what1,what2) RESULT(where)
       !c determines the minimal image translation for the
       !c separation between two atoms
         REAL(q) :: DPOS(:,:)  ! positions in direct
         REAL(q) :: CPOS_1(3),CPOS_2(3)
         REAL(q) :: lattmat(3,3)         !c direct and reciprocal lattice
         REAL(q) :: sdist,dist
         INTEGER :: where(3),trans(3)
         INTEGER :: what1,what2
         INTEGER :: i,j,k

         CPOS_1=DPOS(:,what1)
         CPOS_1=MATMUL(CPOS_1,lattmat)
         CPOS_2=DPOS(:,what2)
         CPOS_2=MATMUL(CPOS_2,lattmat)
         sdist=VECTORSIZE(3,CPOS_2-CPOS_1)
         where=(/0,0,0/)
         dist=0.0
         DO i=-1,1
           DO j=-1,1
             DO k=-1,1
               trans=(/0,0,0/)+(/i,j,k/)
               CPOS_2=DPOS(:,what2)+trans
               CPOS_2=MATMUL(CPOS_2,lattmat)
               dist=VECTORSIZE(3,CPOS_2-CPOS_1)
               IF (dist<sdist) THEN
                 sdist=dist
                 where=(/i,j,k/)
               ENDIF
             ENDDO
           ENDDO
         ENDDO
       END FUNCTION shortest_dist

       SUBROUTINE const_utrans(utrans,ICOORDINATES,CONSTRAINTS)
         TYPE(coordstructure) :: ICOORDINATES
         REAL(q),DIMENSION(:,:) :: utrans
         INTEGER:: CONSTRAINTS(:)
         INTEGER,ALLOCATABLE :: dummy(:),cummy(:)
         INTEGER :: i,dim1,dim2

         CONSTRAINTS=1
         dim1=SIZE(utrans(:,1))
         dim2=SIZE(utrans(1,:))
         ALLOCATE(dummy(dim2))
         ALLOCATE(cummy(dim2))
         dummy=0
         cummy=0
         DO i=1,ICOORDINATES%NUMINTERNALS
           IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==0) dummy(i)=1
           IF (ICOORDINATES%COORDSTRUCT(i)%STATUS==-1) cummy(i)=1 
         ENDDO

         DO i=1,dim1
           IF (SUM(dummy*(ABS(utrans(i,:))))>1e-4) THEN
             CONSTRAINTS(i)=0
           ENDIF
         ENDDO

         DO i=1,dim1
           IF (SUM(cummy*(ABS(utrans(i,:))))>1e-4) THEN
             CONSTRAINTS(i)=-1
           ENDIF
         ENDDO

         DEALLOCATE(cummy)
         DEALLOCATE(dummy)
       END SUBROUTINE const_utrans

       SUBROUTINE KINETIC_E(EKIN,NIONS,NTYP,ITYP,POMASS,POTIM,A,V)
          USE base
          USE lattice
          USE ini
          USE chain
          INTEGER :: NI,NT,NTYP,NIONS
          REAL(q) :: V(3,NIONS),VTMP(3)
          REAL(q) :: A(3,3),POMASS(NTYP)
          INTEGER ::  ITYP(NIONS)
          REAL(q) :: UL,UT,FACT,EKIN,POTIM

          IF (POTIM==0) THEN
            EKIN=0.0_q
            RETURN
          ENDIF

          !c set unit length,unit time ...
          UL=1E-10_q
          UT=POTIM*1E-15_q

          !c this factor converts  (atomic mass* UL**2/UT**2) to UE (eV)
          FACT=(AMTOKG/EVTOJ)*(UL/UT)**2

          !c convert to cartesian coordinetes and calculated EKIN
          EKIN=0
          DO NI=1,NIONS
            VTMP(1)=   V(1,NI)
            VTMP(2)=   V(2,NI)
            VTMP(3)=   V(3,NI)
            CALL  DIRKAR(1,VTMP,A)
            NT=ITYP(NI)
            EKIN=EKIN+ (VTMP(1)**2+VTMP(2)**2+VTMP(3)**2)*POMASS(NT)
          ENDDO
          EKIN= EKIN*FACT/2._q
        END SUBROUTINE KINETIC_E

       SUBROUTINE KINETIC_E_Lat(EKIN,AMASS,POTIM,V)
          USE base
          USE lattice
          USE ini
          USE chain
          REAL(q) :: V(3,3),V_tmp(3,3)
          REAL(q) :: AMASS
          REAL(q) :: UL,UT,FACT,EKIN,POTIM
          INTEGER :: i,j

          IF (POTIM==0) THEN
            EKIN=0.0_q
            RETURN
          ENDIF

          !c set unit length,unit time ...
          UL=1E-10_q
          UT=POTIM*1E-15_q

          !c this factor converts  (atomic mass* UL**2/UT**2) to UE (eV)
          FACT=(AMTOKG/EVTOJ)*(UL/UT)**2

          V_tmp=MATMUL(V,TRANSPOSE(V))
          EKIN=V_tmp(1,1)+V_tmp(2,2)+V_tmp(3,3)
          EKIN=EKIN*AMASS*FACT/2._q
         

!           !c convert to kartesian coordinetes and calculated EKIN
!           EKIN=0._q
!           DO i=1,3 
!             DO j=1,3
!               !!upper traingle only!!!
!               !!IF (j<=i) THEN
!                 EKIN=EKIN+V(i,j)**2*AMASS
!               !!ENDIF
!             ENDDO
!           ENDDO
!           EKIN= EKIN*FACT/2._q
        END SUBROUTINE KINETIC_E_Lat

       SUBROUTINE init_velocities(V,T_INFO,DYN,LATT_CUR,NDEGREES_OF_FREEDOM,LSCALE,LCMASS)
       !c initialization of velocities for the atomic positions
         USE lattice
         USE ini
         TYPE(type_info) :: T_INFO
         TYPE(dynamics) :: DYN
         TYPE(latt) :: LATT_CUR
         REAL(q) :: V(:,:)
         REAL(q) :: CMASS(3,1)
         REAL(q) :: RNULL,BMP,FACT,UL,UT,EKIN,SCALE,TEMPER
         INTEGER :: KP,NT,NDEGREES_OF_FREEDOM,i,j
         LOGICAL :: LSCALE !c rescale velocity so that the resulting T is Tsoll?
         LOGICAL :: LCMASS !c remove shift of center of mass?
         INTEGER :: Ltxyz(3,1)

         RNULL=0._q
         UL =1E-10_q
         UT =DYN%POTIM*1E-15_q
         FACT= (AMTOKG/EVTOJ)*(UL/UT)**2


         !c velocities in cartesian coords 
         V=0._q
         DO KP=1,T_INFO%NIONS
           NT=T_INFO%ITYP(KP)
           BMP=SQRT(DYN%TEMP*BOLKEV/(T_INFO%POMASS(NT)*FACT))
           !IF (T_INFO%LSFOR(1,KP)) V(1,KP)=boltzmann_distribution(RNULL,BMP) 
           !IF (T_INFO%LSFOR(2,KP)) V(2,KP)=boltzmann_distribution(RNULL,BMP)
           !IF (T_INFO%LSFOR(3,KP)) V(3,KP)=boltzmann_distribution(RNULL,BMP)
           V(1,KP)=boltzmann_distribution(RNULL,BMP)
           V(2,KP)=boltzmann_distribution(RNULL,BMP)
           V(3,KP)=boltzmann_distribution(RNULL,BMP)
         ENDDO

         CALL KARDIR(T_INFO%NIONS,V,LATT_CUR%B)

         Ltxyz=1
         !c selective dynamics applies to fractional and NOT cartesian coords!!!
         DO KP=1,T_INFO%NIONS
           IF (.NOT. T_INFO%LSFOR(1,KP)) THEN
             V(1,KP)=0._q
             Ltxyz(1,1)=0
           ENDIF
           IF (.NOT. T_INFO%LSFOR(2,KP)) THEN
             V(2,KP)=0._q
             Ltxyz(2,1)=0
           ENDIF
           IF (.NOT. T_INFO%LSFOR(3,KP)) THEN
             V(3,KP)=0._q
             Ltxyz(3,1)=0
           ENDIF
         ENDDO


         IF (LCMASS) THEN
           CMASS=0._q
           CALL GIVE_CMASS(T_INFO,V,CMASS)
           !c it don't remove CMASS if there is a fixed direction
           DO i=1,3
             IF (Ltxyz(i,1)==0)  CMASS(i,1)=0._q
           ENDDO
           DO i=1,T_INFO%NIONS
             DO j=1,3
                V(j,i)=V(j,i)-CMASS(j,1)
             ENDDO
           ENDDO
         ENDIF

         CALL KINETIC_E(EKIN,T_INFO%NIONS,T_INFO%NTYP,T_INFO%ITYP,T_INFO%POMASS, &
               & DYN%POTIM,LATT_CUR%A,V)

         TEMPER=2*EKIN/BOLKEV/NDEGREES_OF_FREEDOM

         IF (EKIN==0) THEN
           SCALE=0
         ELSE
           SCALE=SQRT(DYN%TEMP/TEMPER)
         ENDIF

         IF (LSCALE) V(:,:)=V(:,:)*SCALE
       END SUBROUTINE init_velocities

       FUNCTION boltzmann_velocity(DYN,LATT_CUR,MASS) RESULT(V)
        !c return boltzmann distributed velocities according to the current
        !c temperature in the array V.
        !c note: no scaling and no drift correction will be performed
          USE lattice
          USE ini
          REAL(q) :: V(3)
          TYPE(dynamics), INTENT(IN) :: DYN
          TYPE(latt), INTENT(IN) :: LATT_CUR
          REAL(q), INTENT(IN) :: MASS

          INTEGER :: J
          REAL(q) :: UL, UT
          REAL(q) :: MU, SIGMA

          UL = 1E-10_q
          UT = DYN%POTIM*1E-15_q
          MU = 0.0_q

          !c the standard deviation is sqrt(k_B*T / m)
          SIGMA = SQRT(DYN%TEMP*BOLK / (MASS*AMTOKG)) * UT / UL
          DO J = 1, 3
            !c the function boltzmann_distribution actually returns a
            !c normal distributed variable X with E(X)=MU and V(X)=SIGMA^2
            V(J) = boltzmann_distribution(MU,SIGMA)
          END DO

          CALL KARDIR(1,V,LATT_CUR%B)
       END FUNCTION


       SUBROUTINE boltzmann_velocities(T_INFO,DYN,LATT_CUR,V)
        !c return boltzmann distributed velocities according to the current
        !c temperature in the array V.
        !c note: no scaling and no drift correction will be performed
          USE lattice
          USE ini
          TYPE(type_info), INTENT(IN) :: T_INFO
          TYPE(dynamics), INTENT(IN) :: DYN
          TYPE(latt), INTENT(IN) :: LATT_CUR
          REAL(q), INTENT(OUT) :: V(3,T_INFO%NIONS)

          INTEGER :: I, IT, J
          REAL(q) :: UL, UT
          REAL(q) :: MU, SIGMA

          UL = 1E-10_q
          UT = DYN%POTIM*1E-15_q
          MU = 0.0_q

          V = 0.0_q
          DO I = 1, T_INFO%NIONS
            IT = T_INFO%ITYP(I)
            !c the standard deviation is sqrt(k_B*T / m)
            SIGMA = SQRT(DYN%TEMP*BOLK / (T_INFO%POMASS(IT)*AMTOKG)) * UT / UL
            DO J = 1, 3
              !c the function boltzmann_distribution actually returns a
              !c normal distributed variable X with E(X)=MU and V(X)=SIGMA^2
              !IF (T_INFO%LSFOR(J,I)) V(J,I) = boltzmann_distribution(MU,SIGMA)
              V(J,I) = boltzmann_distribution(MU,SIGMA)
            END DO
          END DO

          CALL KARDIR(T_INFO%NIONS,V,LATT_CUR%B)
          DO I = 1, T_INFO%NIONS
            DO J = 1, 3
              IF (.NOT. T_INFO%LSFOR(J,I)) V(J,I) = 0._q
            END DO
          END DO

       END SUBROUTINE

       SUBROUTINE friction_Forces_stoch(V,T_INFO,DYN,LATT_CUR,gamma)
       !c contribution of stochastic friction forces to acceleration (Langevin dynamics)
       !c we actually compute (1/2*a*dt^2)
         USE lattice
         USE ini
         TYPE(type_info) :: T_INFO
         TYPE(dynamics) :: DYN
         TYPE(latt) :: LATT_CUR
         REAL(q) :: V(:,:) !c stochastic force 
         REAL(q) :: CMASS(3,1)
         REAL(q) :: RNULL,BMP,FACT
         INTEGER :: KP,NI,NT,i,j
         REAL(q) :: GAMMA(T_INFO%NTYP) ! friction coefs in ps^(-1)
        
         RNULL=0._q
        
         FACT=1e12*EVTOJ*AMTOKG/1E-15
    
         V=0._q
         DO KP=1,T_INFO%NIONS
           NT=T_INFO%ITYP(KP)
           BMP=2*gamma(NT)*BOLKEV*DYN%TEMP*T_INFO%POMASS(NT)/DYN%POTIM
           !c force in J/m:
           BMP=BMP*FACT
           BMP=SQRT(BMP)
           V(1,KP)=boltzmann_distribution(RNULL,BMP)
           V(2,KP)=boltzmann_distribution(RNULL,BMP)
           V(3,KP)=boltzmann_distribution(RNULL,BMP)
         ENDDO

!          !c transform force to eV/A
!          V=V/EVTOJ*1E-10_q
! 
!          !c now compute scaled acceleration that is compatible with that
!          !c in DYN%D2C:
!          FACT=(DYN%POTIM**2)*EVTOJ/AMTOKG *1E-10_q

         FACT=(DYN%POTIM**2)/AMTOKG *1E-20_q

         DO KP=1,T_INFO%NIONS
           NT=T_INFO%ITYP(KP)
           V(1,KP)=V(1,KP)*FACT/2/T_INFO%POMASS(NT)
           V(2,KP)=V(2,KP)*FACT/2/T_INFO%POMASS(NT)
           V(3,KP)=V(3,KP)*FACT/2/T_INFO%POMASS(NT)
         ENDDO

         !c scalled acceleration for fractional coordinates
         CALL KARDIR(T_INFO%NIONS,V,LATT_CUR%B)

         !c take care of fixed coordinates
         DO I = 1, T_INFO%NIONS
            DO J = 1, 3
              IF (.NOT. T_INFO%LSFOR(J,I)) V(J,I) = 0._q
            END DO
         END DO
       END SUBROUTINE friction_Forces_stoch

      SUBROUTINE friction_Forces_Lstoch(V,DYN,GAMMA_L,AMASS)
       !c contribution of stochastic friction forces to acceleration (Langevin dynamics)
       !c we actually compute (1/2*a*dt^2)
         USE lattice
         USE ini
         TYPE(dynamics) :: DYN
         REAL(q) :: V(3,3) !c stochastic force 
         REAL(q) :: RNULL,BMP,FACT
         INTEGER :: i,j
         REAL(q) :: GAMMA_L ! friction coefs in ps^(-1)
         REAL(q) :: AMASS ! mass for lattice parameters
        
         RNULL=0._q        
         FACT=1e12*EVTOJ*AMTOKG/1E-15
         BMP=2*GAMMA_L*BOLKEV*DYN%TEMP*AMASS/DYN%POTIM
         BMP=BMP*FACT
         BMP=SQRT(BMP)
      
         !c friction force in J/m
         V=0._q
         DO i=1,3
           DO j=1,3
             V(i,j)=boltzmann_distribution(RNULL,BMP) 
           ENDDO
         ENDDO

         !c force in eV/A
         V=V/EVTOJ*1E-10_q

         !c compute one half of acceleration
         FACT=(DYN%POTIM**2)*EVTOJ/AMTOKG *1E-10_q
         V=V*FACT/2/Amass 
       END SUBROUTINE friction_Forces_Lstoch

       SUBROUTINE friction_Forces_momentum(V,T_INFO,DYN,gamma)
       !c contribution of momentum dependent friction forces to acceleration (Langevin dynamics)
       !c for the atomic positions
       !c we actually compute (1/2*a*dt^2)
         USE lattice
         USE ini
         TYPE(type_info) :: T_INFO
         TYPE(dynamics) :: DYN
         REAL(q) :: V(:,:) !c momentum dependent friction force 
         INTEGER :: KP,NT
         REAL(q) :: GAMMA(T_INFO%NTYP) ! friction coefs in ps^(-1)
  
         !acceleration in (fractional coord.)/fs/ps
         V=0._q
         DO KP=1,T_INFO%NIONS
           NT=T_INFO%ITYP(KP)
           IF (T_INFO%LSFOR(1,KP)) V(1,KP)=DYN%VEL(1,KP)*GAMMA(NT)*DYN%POTIM
           IF (T_INFO%LSFOR(2,KP)) V(2,KP)=DYN%VEL(2,KP)*GAMMA(NT)*DYN%POTIM
           IF (T_INFO%LSFOR(3,KP)) V(3,KP)=DYN%VEL(3,KP)*GAMMA(NT)*DYN%POTIM
         ENDDO

         !c now compute one-half of the acceleration in (fractional coord.)/fs^2
         V=-1*V*1e-3/2
       END SUBROUTINE friction_Forces_momentum

       SUBROUTINE friction_Forces_Lmomentum(V,AVEL,DYN,GAMMA_L)
       !c contribution of momentum dependent friction forces to acceleration (Langevin dynamics)
       !c for the lattice components
       !c we actually compute (1/2*a*dt^2)
         USE lattice
         USE ini
         TYPE(dynamics) :: DYN
         REAL(q) :: AVEL(3,3) !c lattice velocities
         REAL(q) :: V(:,:) !c momentum dependent friction force 
         INTEGER :: i,j
         REAL(q) :: GAMMA_L ! friction coefs in ps^(-1)
  
         !acceleration in (fractional coord.)/fs/ps
         V=0._q
         DO i=1,3   
           DO j=1,3   
             !!IF (i>=j) THEN
               V(i,j)=AVEL(i,j)*GAMMA_L*DYN%POTIM
             !!ENDIF
           ENDDO
         ENDDO

         !c now compute one-half of the acceleration in A/fs^2
         V=-1*V*1e-3/2
       END SUBROUTINE friction_Forces_Lmomentum

       SUBROUTINE langevin_temperature(EKIN1,EKIN2,TEIN1,TEIN2,T_INFO,DYN,A,GAMMA)
       !c contribution of momentum dependent friction forces to acceleration (Langevin dynamics)
       !c we actually compute (1/2*a*dt^2)
         USE lattice
         USE ini
         TYPE(type_info) :: T_INFO
         TYPE(dynamics) :: DYN
         REAL(q) :: VTMP(3) !c momentum dependent friction force 
         REAL(q) :: A(3,3)
         INTEGER :: NI,NT,I
         REAL(q) :: GAMMA(T_INFO%NTYP) ! friction coefs in ps^(-1)
         REAL(q) :: UL,UT,FACT
         REAL(q) :: ekin1,ekin2,tein1,tein2
         INTEGER :: count1,count2
  
         UL=1E-10_q
         UT=DYN%POTIM*1E-15_q

          !c this factor converts  (atomic mass* UL**2/UT**2) to UE (eV)
         FACT=(AMTOKG/EVTOJ)*(UL/UT)**2

         ekin1=0._q;ekin2=0._q
         count1=0;count2=0

         DO NI=1,T_INFO%NIONS
            VTMP(1)=   DYN%VEL(1,NI)
            VTMP(2)=   DYN%VEL(2,NI)
            VTMP(3)=   DYN%VEL(3,NI)
            CALL  DIRKAR(1,VTMP,A)
            NT=T_INFO%ITYP(NI)
            IF (GAMMA(NT)>0._q) THEN
              DO I=1,3
                IF (T_INFO%LSFOR(I,NI)) THEN 
                  count1=count1+1
                  EKIN1=EKIN1+(VTMP(I))**2*T_INFO%POMASS(NT)
                ENDIF
              ENDDO
            ELSE
              DO I=1,3
                IF (T_INFO%LSFOR(I,NI)) THEN 
                  count2=count2+1
                  EKIN2=EKIN2+(VTMP(I))**2*T_INFO%POMASS(NT)
                ENDIF
              ENDDO
            ENDIF
          ENDDO

          !write(*,*) 'test_count:',count1,count2
          tein1=0._q;tein2=0._q
          EKIN1= EKIN1*FACT/2._q
          EKIN2= EKIN2*FACT/2._q
          IF (count1 .GT. 0) TEIN1=2*EKIN1/BOLKEV/count1
          IF (count2 .GT. 0) TEIN2=2*EKIN2/BOLKEV/count2
       END SUBROUTINE langevin_temperature

       SUBROUTINE init_Avelocities(V,DYN,mass)
       !c initialization of velocities for the lattice parameters
         USE lattice
         USE ini
         TYPE(dynamics) :: DYN
         REAL(q) :: V(3,3)
         REAL(q) :: CMASS(3,1)
         REAL(q) :: RNULL,BMP,FACT,UL,UT,EKIN,SCALE,TEMPER
         INTEGER :: i,j
         REAL(q) :: mass
        
         RNULL=0._q
         UL =1E-10_q
         UT =DYN%POTIM*1E-15_q
         FACT= (AMTOKG/EVTOJ)*(UL/UT)**2
      
         V=0._q
         !write(*,*) 'V_0',V
         DO i=1,3
           DO j=1,3
             !!upper traingle only!!!
             !IF (i>=j) THEN
               BMP=SQRT(DYN%TEMP*BOLKEV/(mass*FACT))
               V(i,j)=boltzmann_distribution(RNULL,BMP)
             !END IF
           ENDDO 
         ENDDO
 
       END SUBROUTINE init_Avelocities

       FUNCTION boltzmann_distribution(RNULL,WIDTH)
         USE random_seeded
         REAL(q) :: num1,num2
         REAL(q) :: boltzmann_distribution,WIDTH,RNULL

         CALL RANE_ION(num1)
         num2=0._q
         DO
           CALL RANE_ION(num2)
           num2=ABS(num2)
           IF (num2 .GT. 1E-08) EXIT
         ENDDO
         IF (num2 .GT. 1._q) num2=1._q

         boltzmann_distribution= COS( TPI*num1 ) * SQRT( 2._q*ABS(LOG(num2)) )
         boltzmann_distribution = WIDTH * boltzmann_distribution  +  RNULL
       END FUNCTION boltzmann_distribution

       FUNCTION gasdev(idum)
         USE random_seeded
         integer :: idum
         real(q) :: fac,rsq,v1,v2,gasdev
         real(q) :: ran2
         integer, save :: iset=0
         real(q), SAVE :: gset

        if(idum.lt.0) iset=0
        if(iset.eq.0) then
          rsq=10._q
          DO
            CALL RANE_ION(v1)
            CALL RANE_ION(v2)
            v1 = 2.*v1 - 1._q
            v2 = 2.*v2 - 1._q
            rsq = v1**2 + v2**2
            IF (rsq .LE. 1._q .AND. rsq .GT. 0._q) EXIT
          ENDDO
          fac = sqrt(-2.*log(rsq)/rsq)
          gset = v1*fac
          gasdev = v2*fac
          iset = 1
        else
          gasdev = gset
          iset = 0
        end if
        return
      END FUNCTION gasdev



      SUBROUTINE GIVE_CMASS(T_INFO,X,CMASS)
          USE base
          USE lattice
          TYPE(type_info) :: T_INFO
          REAL(q) :: X(3,T_INFO%NIONS),X_(3,T_INFO%NIONS)
          REAL(q) :: CMASS(3,1)
          REAL(q) :: MASS
          REAL(q) :: TOTALMASS
          INTEGER :: i,j
          INTEGER :: Ltxyz(3)

          !c are there any fixed atoms?
          Ltxyz=1
          DO i=1,T_INFO%NIONS
            DO j=1,3
              IF (.NOT. T_INFO%LSFOR(j,i) .AND. Ltxyz(j)==1) Ltxyz(j)=0
            ENDDO
          ENDDO

          TOTALMASS=0._q
          MASS=0._q
          X_=0._q
          DO i=1,T_INFO%NIONS
            MASS=(T_INFO%POMASS(T_INFO%ITYP(i)))
            TOTALMASS=TOTALMASS+MASS
            DO j=1,3
              X_(j,i)=X(j,i)*MASS
            ENDDO
          ENDDO
          DO i=1,3
            CMASS(i,:)=SUM(X_(i,:))*Ltxyz(i)/TOTALMASS
          ENDDO
        END SUBROUTINE GIVE_CMASS

        SUBROUTINE put_in_box(NIONS,X)
          REAL(q) :: X(3,NIONS)
          INTEGER :: NIONS
          INTEGER :: i,j

          DO i=1,3
            DO j=1,NIONS
              DO
                IF (X(i,j) .GE. 1._q) THEN
                  X(i,j)=X(i,j)-1._q
                ELSE IF (X(i,j) .LT. 0._q) THEN
                  X(i,j)=X(i,j)+1._q
                ELSE
                  EXIT
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        END SUBROUTINE put_in_box

        SUBROUTINE minimize_difference(X,Y,NIONS)
          REAL(q) :: X(3,NIONS)
          REAL(q) :: Y(3,NIONS)
          INTEGER :: NIONS
          INTEGER :: i,j
  
          DO i=1,3
            DO j=1,NIONS
              DO
                IF (X(i,j)-Y(i,j) .GT. 0.5_q) THEN
                  X(i,j)=X(i,j)-1._q
                ELSE IF (X(i,j)-Y(i,j) .LE. -0.5_q) THEN
                  X(i,j)=X(i,j)+1._q
                ELSE
                  EXIT
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        END SUBROUTINE minimize_difference

!         FUNCTION ERRF(x)
!           INTEGER :: SIG
!           REAL(q) :: x,errf,t,y
!           REAL(q),PARAMETER :: a1 =  0.254829592
!           REAL(q),PARAMETER :: a2 = -0.284496736
!           REAL(q),PARAMETER :: a3 =  1.421413741
!           REAL(q),PARAMETER :: a4 = -1.453152027
!           REAL(q),PARAMETER :: a5 =  1.061405429
!           REAL(q),PARAMETER :: p  =  0.3275911
! 
!           SIG = 1
!           IF (x < 0.) SIG = -1
!           x = ABS(x)
!           t = 1._q/(1._q + p*x)
!           y = 1._q - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*EXP(-x*x)
!           errf=SIG*y
!         END FUNCTION ERRF

   SUBROUTINE min_imageV(n,dxn)
    !c minimum image convention for vector
     USE prec
     INTEGER ::n,i
     REAL(q) :: dxn(n)
  
      !write(*,*) 'dx1',dxn 
      DO i=1,n
        DO 
          IF (dxn(i) .GT. 0.5_q) THEN
            dxn(i)=dxn(i)-1._q
          ELSE IF (dxn(i) .LE. -0.5_q) THEN
            dxn(i)=dxn(i)+1._q
          ELSE 
            EXIT
          END IF
        END DO
      ENDDO
    END SUBROUTINE min_imageV

    SUBROUTINE min_imageV_shift(n,dxn,shift)
    !c minimum image convention for vector
     USE prec
     INTEGER ::n,i
     REAL(q) :: dxn(n)
     INTEGER :: shift(n)

     shift=0
      !write(*,*) 'dx1',dxn 
      DO i=1,n
        DO
          IF (dxn(i) .GT. 0.5_q) THEN
            dxn(i)=dxn(i)-1._q
            shift(i)=-1
          ELSE IF (dxn(i) .LE. -0.5_q) THEN
            dxn(i)=dxn(i)+1._q
            shift(i)=1
          ELSE
            EXIT
          END IF
        END DO
      ENDDO
    END SUBROUTINE min_imageV_shift



   SUBROUTINE index_table2d(n,nn,p,q)
     !nn=n*(n+1)/2  
     INTEGER, intent(in) :: n,nn
     INTEGER :: i, F,G, gg,ff
     INTEGER,intent(out) :: p(nn),q(nn)
    
     p=0;q=0
     F=n
     G=1
     gg=G
     ff=n-1
     DO i=1,nn
       IF (i .GT. F) THEN
         G=G+1
         F=F+ff
         ff=ff-1
         gg=G
       END IF
       p(i)=G
       q(i)=gg
       gg=gg+1
     END DO
  
   END SUBROUTINE index_table2d

   SUBROUTINE gauss_chebyshevTB(nb, omega, weight)
   !c not tested yet!!!
     implicit none
     integer, intent(in)    :: nb
     REAL(q), intent(inout) :: omega(nb+1), weight(nb+1)
     integer                :: i, j,nnb
     REAL(q)                :: t1, t2, t3

     omega=0._q
     weight=0._q

     t1=pi/(4._q*nb)
     t2=2._q*t1
     DO i=1,nb
       j=nb+1-i
       t3=t1*(2._q*i-1._q)
       omega(i+1) = 1._q/tan(t3)
       weight(i+1) = t2/(sin(t3))**2
     ENDDO
   
   END SUBROUTINE gauss_chebyshevTB

   FUNCTION VERBALIZE_DFTMETHOD() RESULT(DFTmethod)
   !c determine the DFT method used in calculation,
   !c identification of only those methods which are supported by 
   !c at least one vdW correction method is attempted
     USE setexm
     USE fock
     USE setxcmeta
     CHARACTER (LEN=10) :: DFTmethod

     DFTmethod='unknown'
     !c SCAN metaGGA
     IF (ID_METAGGA==50) THEN
       DFTmethod='scan'
     ELSE IF (ID_METAGGA==35) THEN
       DFTmethod='m06l'
     ELSE IF (ID_METAGGA==25) THEN
       DFTmethod='tpss'
     !ELSE IF ((LEXCH==4) .AND. (ID_METAGGA==-1)) THEN
     !  !c PB
     !  IF (LDAX==1.0) THEN
     !    DFTmethod='bp'
     !  ENDIF
     ELSE IF ((LEXCH==8) .AND. (ID_METAGGA==-1)) THEN
       !c standard PBE
       IF (LDAX==1.0_q) THEN
         DFTmethod='pbe'
       !c PBE0
       ELSE IF (AEXX==0.25_q .AND. AGGAX==0.75_q .AND. AGGAC==1.0_q .AND. ALDAC==1.0_q .AND. HFSCREEN==0.0_q) THEN
         DFTmethod='pbe0'
       !c HSE03
       ELSE IF (AEXX==0.25_q .AND. AGGAX==0.75_q .AND. AGGAC==1.0_q .AND. ALDAC==1.0_q .AND. HFSCREEN==0.3_q) THEN
         DFTmethod='hse03'
       !c HSE06
       ELSE IF (AEXX==0.25_q .AND. AGGAX==0.75_q .AND. AGGAC==1.0_q .AND. ALDAC==1.0_q .AND. HFSCREEN==0.2_q) THEN
         DFTmethod='hse06'
       ELSE IF (AEXX==1.00_q .AND. AGGAC==0.0_q .AND. ALDAC==0.0_q) THEN
         DFTmethod='hf'
       ENDIF
     ELSE IF ((LEXCH==9) .AND. (ID_METAGGA==-1)) THEN
       IF (LDAX==1.0_q) THEN
         DFTmethod='rpbe'
       ENDIF
     ELSE IF ((LEXCH==11) .AND. (ID_METAGGA==-1)) THEN
       !c B3LYP
       IF (AEXX==0.20_q .AND. AGGAX==0.72_q .AND. AGGAC==0.81_q .AND. ALDAC==0.19_q) THEN
         DFTmethod='b3lyp'
       ENDIF
     ELSE IF ((LEXCH==14) .AND. (ID_METAGGA==-1)) THEN
       IF (LDAX==1.0_q) THEN
         DFTmethod='pbesol'
       ENDIF
     ELSE IF ((LEXCH==40) .AND. (ID_METAGGA==-1)) THEN
       IF (LDAX==1.0_q) THEN
         DFTmethod='revpbe'
       ENDIF
     ENDIF

   END FUNCTION VERBALIZE_DFTMETHOD 

   FUNCTION GREATEST_COMMON_DIVISOR(NA, NB)
     INTEGER GREATEST_COMMON_DIVISOR
     INTEGER NA, NB
   ! local
     INTEGER IA, IB, ITEMP
     
     IA = NA
     IB = NB
1    IF (IB.NE.0) THEN
        ITEMP = IA
        IA = IB
        IB = MOD(ITEMP, IB)
        GOTO 1
     END IF
     GREATEST_COMMON_DIVISOR = IA
     RETURN
   END FUNCTION GREATEST_COMMON_DIVISOR

   END MODULE
