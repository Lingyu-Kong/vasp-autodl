
!> @brief This module provides a basic implementation of the
!> Direct inversion of the iterative subspace or pulay mixing
!> for the accelaration of SCF type calculations
!> by Henrique Miranda
!>
!> The idea is to solve a non-linear problem of the type:
!>   F(p) = p
!> where p is a m dimensional vector.
!> The procedure goes as follows, we start with an intial guess x_1,in and evaluate:
!>   F(p_1,in) = p_1,out
!> the p_2,in for the next iteration is obtained by linear mixing:
!>   p_2,in = alpha*p_1,out + (1-alpha)*p_1,in   (1)
!> Once we have the input and output vectors of nmax iterations we start using Pulay mixing:
!> In Pulay mixing we assume that a new trial vector can be generated as a linear combination of
!> previous nmax trial vectors:
!>   p_i+1 = sum_i^nmax c_i * p_i   (2)
!> The c_i weights are determined by solving the following linear system of equations:
!>   Bc=b
!> Where:
!>     B_ij = <p_i-p_i-1|p_j-p_j-1> and B_nmax+1,* = -1 and B_*,nmax+1 = -1
!>  and
!>     b_i = 0 with 1<i<nmax and b_nmax+1 = -1
!> the resulting vector c_i with 1<i<nmax gives the weights to use in (2)
!> and c_nmax+1 is the norm of the residual vector.
!> We only use pulay mixing once we have nmax > 2, the first step is done by a special linear step (see diis_pulay)
!> To ensure B does not become simular we add DIAG_EPS to the diagonal elements.
!> If the pulay jump leads to a change of p that is larger than the one in the previous step
!> then we take a simple linear step as in (1).
!> Following http://dx.doi.org/10.1016/j.cplett.2015.06.029 we restart the history
!> once the number of vectors p in the history is equal to nmax
module diis
    use prec
    implicit none
    real(q),parameter :: TOL = 1e-8_q !> Tolerance to check whether the weights sum to one
    real(q),parameter :: DIAG_EPS = 0.01 !> Number to add to the diagonal of B to prevent it to become singular

    !> @brief DIIS datatype
    type :: diis_type
        integer :: n !< current number of vectors
        integer :: nmax !< maximum number of steps in the history
        integer :: m !< size of the vector
        logical :: debug !< whether to use debug mode
        real(q) :: alpha !< mixing factor for linear mixer
        real(q) :: dw !< maximum residual
        real(q),allocatable :: p0(:) !< step for linear mixing
        real(q),allocatable :: p(:,:) !< iteration history
        real(q),allocatable :: b(:,:) !< B matrix to solve the linear problem
        real(q),allocatable :: c(:) !< coefficients obtained from the linear problem
    end type

    contains
    !> @brief Initialize DIIS object
    subroutine diis_init(self,m,nmax,alpha)
        type(diis_type) :: self
        real(q) :: alpha
        integer :: m, nmax
        self%n = -1
        self%nmax = nmax
        self%m = m
        self%alpha = alpha
        self%debug = .false.
        self%dw = 1e8
        allocate(self%p0(self%m))
        allocate(self%p(self%m,0:self%nmax))
        allocate(self%c(self%nmax+1))
        allocate(self%b(self%nmax+1,self%nmax+1))
    end subroutine diis_init

    !> @brief Save this step for linear mixing
    subroutine diis_save(self,p)
        type(diis_type) :: self
        real(q) :: p(self%m)
        self%p0 = p
    end subroutine diis_save

    !> @brief Add a vector to the history
    subroutine diis_append(self,p)
        type(diis_type) :: self
        real(q) :: p(self%m)
        ! local variables
        integer :: i

        ! This is based of from http://dx.doi.org/10.1016/j.cplett.2015.06.029
        if (self%n==self%nmax) call diis_reset(self)
        if (self%n<self%nmax) then
            ! if the number of iterations is smaller than iteration
            ! history we just need to add a new vector
            self%n = self%n + 1
            self%p(:,self%n)  = p
        else
            ! move all the vectors one level
            ! TODO: this should not involve memory copies!
            do i=0,self%nmax-1
                self%p(:,i) = self%p(:,i+1)
            enddo
            self%p(:,self%nmax)  = p
        endif
    end subroutine diis_append

    !> @brief Reset history
    subroutine diis_reset(self)
        type(diis_type) :: self
        ! local variables
        integer :: i, shift
        shift = min(self%n,self%nmax)
        self%n=1
        ! move all the vectors to the begining
        ! TODO: this should not involve memory copies!
        do i=0,self%n
            self%p(:,i) = self%p(:,shift-self%n+i)
        enddo
    end subroutine diis_reset

    !> @brief Edit the current vector in history
    subroutine diis_edit(self,p)
        type(diis_type) :: self
        real(q) :: p(self%m)
        if (self%n==-1) then
           self%p(:,0) = p
        else if (self%n<self%nmax) then
           self%p(:,self%n) = p
        else
           self%p(:,self%nmax) = p
        endif
    end subroutine diis_edit

    !> @brief Print the iteration history
    !> Usually this is not practical as the vectors are too large
    subroutine diis_print(self)
        type(diis_type) :: self
        ! local variables
        integer :: i
        do i=1,min(self%nmax,self%n)
            write(*,*) i,self%p(:,i)
        enddo
    end subroutine diis_print

    !> @brief Perform a step using simple linear mixing
    subroutine diis_simple(self,new_p)
        type(diis_type) :: self
        real(q),intent(out) :: new_p(self%m)
        ! local variables
        integer :: step
        step = min(self%n,self%nmax)
        new_p = self%alpha*self%p(:,step) + (1.0_q-self%alpha)*self%p0
        self%dw = maxval(abs(self%p(:,step)-self%p0))
    end subroutine diis_simple

    !> @brief Perform a step using pulay mixing
    subroutine diis_pulay(self,new_p,ierr)
        type(diis_type) :: self
        real(q),intent(out) :: new_p(self%m)
        integer :: ierr
        ! local variables
        integer :: i,j,info
        integer :: bsize,step
        integer :: ipiv(self%nmax+1)
        real(q) :: b(self%nmax+1,self%nmax+1)
        real(q) :: prod, dw
        ierr = 0
        bsize = min(self%n+1,self%nmax+1)
        if (self%n<2) then
            ! in case there is not enough history to start the pulay method we take a linear step
            step = min(self%n,self%nmax)
            new_p = self%alpha*self%p(:,step) + (1.0_q-self%alpha)*self%p(:,step-1)
            self%dw = maxval(abs(new_p-self%p0))
            return
        endif
        ! Compute B
        self%b = 0
        do i=1,bsize-1
            do j=1,bsize-1
                if (j>i) cycle
                prod = dot_product(self%p(:,j)-self%p(:,j-1),&
                                   self%p(:,i)-self%p(:,i-1))
                self%b(j,i) = prod
                self%b(i,j) = prod
            enddo
            self%b(bsize,i) = -1
            self%b(i,bsize) = -1
            self%b(bsize,bsize) = 0
            self%b(i,i) = self%b(i,i) + DIAG_EPS
        enddo
        ! compute B and b
        self%c = 0
        self%c(min(self%n+1,self%nmax+1)) = -1
        ! solve Bc = b
        b = self%b
        call dgesv(bsize,1,b,self%nmax+1,ipiv,self%c,self%nmax+1,info)
        if (info/=0) then
            if (self%debug) write(*,*) 'dgesv error', info
            ierr = 2
            return
        endif
        if (abs(1.0_q-sum(self%c(:bsize-1)))>TOL) then
            if (self%debug) write(*,*) 'the weights must sum to one', sum(self%c)
            ierr = 3
            return
        endif
        ! print after solving
        if (self%debug) then
            write(*,*) 'iter',self%n
            do i=1,self%nmax+1
                write(*,'(8F8.4)') self%b(:,i), self%c(i)
            enddo
        endif
        ! compute step
        new_p = 0
        do i=1,bsize-1
            new_p = new_p + self%p(:,i)*self%c(i)
        enddo
        dw = maxval(abs(new_p-self%p0))
        ! if the residual is larger then previous step then we use linear mixing
        if (dw>self%dw) then
            call diis_simple(self,new_p)
        else
            self%dw = dw
        endif
    end subroutine diis_pulay

    !> @brief External interface to obtain the next guess for the next step
    !> In the first steps (until nmax) a simple linear mixing is used,
    !> afterwards the pulay mixing is used
    subroutine diis_step(self,new_p)
        type(diis_type) :: self
        real(q),intent(out) :: new_p(self%m)
        ! local variables
        integer :: ierr
        logical :: pulay
        pulay = .true.
        if (self%n<1) then
            write(*,*) 'cannot make step first must call diis_append'
            call exit(0)
        endif
        call diis_pulay(self,new_p,ierr)
        if (pulay) then
            if (self%debug) write(*,*) 'pulay mixing'
        else
            if (self%debug) write(*,*) 'simple mixing'
        endif
    end subroutine diis_step

    !> @brief Free the DIIS structure
    subroutine diis_free(self)
        type(diis_type) :: self
        deallocate(self%p)
        deallocate(self%b)
        deallocate(self%c)
    end subroutine diis_free

end module diis
