!NEC$ options "-no-ftrace"

#include "symbol.inc"
module random_seeded
!===========================================================================
! RCS:  $Id: random.F,v 1.1 2000/11/15 08:13:54 kresse Exp $
!
! This random number generator originally appeared in Toward a Universal
! Random Number Generator by George Marsaglia and Arif Zaman.
! Florida State University Report: FSU-SCRI-87-50 (1987)
!
! It was later modified by F. James and published in A Review of Pseudo-
! random Number Generators
!
! Some final small modifications have been done by J. Furthmueller
! Technical University of Vienna, November 1993
!
! THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
!       (However, a newly discovered technique can yield
!         a period of 10^600. But that is still in the development stage.)
!
! It passes ALL of the tests for random number generators and has a period
!   of 2^144, is completely portable (gives bit identical results on all
!   machines with at least 24-bit mantissas in the floating point
!   representation).
!
! The algorithm is a combination of a Fibonacci sequence (with lags of 97
!   and 33, and operation "subtraction plus one, modulo one") and an
!   "arithmetic sequence" (using subtraction).
!
! On a Vax 11/780, this random number generator can produce a number in
!    13 microseconds.
! (Note by J. Furthmueller: in 2.5 microseconds on a IBM RS6000/Model 580)
!========================================================================
    use prec, only: q
    use constant

    implicit none

    real(q), parameter :: twopi = 6.283185307179586_q
    real(q), parameter :: cd = 7654321.0_q / 16777216.0_q
    real(q), parameter :: cm = 16777213.0_q / 16777216.0_q
    integer, parameter :: max_ij = 31328, max_kl = 30081
    integer, parameter :: countlarge = 1000000

    type, public :: random_number_generator
        logical, private :: initialized = .false.
        real(q), private :: u(97)
        real(q), private :: c = 362436.0_q / 16777216.0_q
        integer, private :: i97 = 97, j97 = 33

#ifdef __NEC_TUNE__
!NEC_RP: CHK_SZ is choosen by initial values of I97 and J97 ensuring that
!        I97_L and J97_L are always different
        integer, private :: rand_blksz = 33
        integer, private :: rand_cnt = 97
        real(q), private :: rand_vec(97)
        integer, private :: c_i = 362436
        integer, private :: cd_i = 7654321
        integer, private :: cm_i = 16777213
        real(q), private :: c_div = 1.0_q / 16777216.0_q
#endif // __NEC_TUNE__

    contains
        generic, public :: random => random_real_scalar, random_real_array
        generic, public :: uniform => uniform_real_scalar, uniform_real_array, &
            uniform_complex_scalar, uniform_complex_array
        procedure, private :: initialize
        procedure, private :: random_real_scalar, random_real_array
        procedure, private :: uniform_real_scalar, uniform_real_array, &
            uniform_complex_scalar, uniform_complex_array
    end type random_number_generator

#ifndef __NEC_TUNE__
    type(random_number_generator) gen_wave, gen_ion
#else // __NEC_TUNE__
    type(random_number_generator), public :: gen_wave, gen_ion
#endif // __NEC_TUNE__

    private
    public :: rane_ion, rang_ion, rane, rang

contains

    real(q) function rane(reset)
        logical, intent(in), optional :: reset
        if (present(reset)) gen_wave%initialized = .not.reset
        call gen_wave%random(rane)
    end function rane

    real(q) function rang(center, width)
        ! This subroutine produces a normal distribution (Gaussian distribution)
        real(q), intent(in) :: center, width
        call gen_wave%uniform(rang, center, width)
    end function rang

    subroutine rane_ion(rng_ion, put, get)
        real(q), intent(out) :: rng_ion
        integer, intent(in), optional :: put(3)
        integer, intent(out), optional :: get(3)
        integer, save :: rng_state(3) = [54217137, 0, 0]
        integer, parameter :: seed1_max = 900000000
            ! The first state stores the seeds ij * (max_kl + 1) + kl, the other two
            ! the number of random numbers generated so far:
            ! rng_state(3) * countlarge + rng_state(2)
        if (present(get)) then
            ! just return the RNG seed (and initialize if needed)
            get = rng_state
        else if (present(put)) then
            ! ensure that seed(1) doesn't exceed max. sensible value
            rng_state(1) = mod(put(1), seed1_max)
            rng_state(2) = mod(put(2), countlarge)
            rng_state(3) = put(3) + put(2) / countlarge
            ! reseed and reinitialise if SEED is provided
            call reset(gen_ion, rng_state)
        else
            ! return random number if no special action is taken
            rng_state(2) = mod(rng_state(2) + 1, countlarge)
            if (rng_state(2) == 0) rng_state(3) = rng_state(3) + 1
            call gen_ion%random(rng_ion)
        end if
    end subroutine rane_ion

    real(q) function rang_ion(center, width)
        ! This subroutine produces a normal distribution (Gaussian distribution)
        real(q), intent(in) :: center, width
        call gen_ion%uniform(rang_ion, center, width)
    end function rang_ion

    subroutine reset(generator, rng_state)
        class(random_number_generator), intent(out) :: generator
        integer, intent(in) :: rng_state(3)
        integer seed_ij, seed_kl, ii, jj
        real(q) skip
        !c regenerate pair of seeds used originally
        seed_ij = rng_state(1) / (max_kl + 1)
        seed_kl = rng_state(1) - (max_kl + 1) * seed_ij
        call generator%initialize(seed_ij, seed_kl)
        ! skip rng_state3 * countlong + rng_state2 numbers (presumably already used
        ! in a previous run) on initialization
        do ii = 1, rng_state(3)
            do jj = 1, countlarge
                call generator%random(skip)
            end do
        end do
        do ii = 1, rng_state(2)
            call generator%random(skip)
        end do
    end subroutine reset

    subroutine initialize(generator, seed_ij, seed_kl)
        ! This is the initialization routine for the random number generator RANMAR()
        ! NOTE: The seed variables can have values between:    0 <= IJ <= 31328
        !                                                      0 <= KL <= 30081
        ! The random number sequences created by these two seeds are of sufficient
        ! length to complete an entire calculation with. For example, if several
        ! different groups are working on different parts of the same calculation,
        ! each group could be assigned its own IJ seed. This would leave each group
        ! with 30000 choices for the second seed. That is to say, this random
        ! number generator can create 900 million different subsequences -- with
        ! each subsequence having a length of approximately 10^30.
        !
        ! Use IJ = 1802 & KL = 9373 to test the random number generator. The
        ! subroutine RANMAR should be used to generate 20000 random numbers.
        ! Then display the next six random numbers generated multiplied by 4096*4096
        ! If the random number generator is working properly, the random numbers
        ! should be:
        !           6533892.0  14220222.0  7275067.0
        !           6172232.0  8354498.0   10633180.0
        use tutor, only: vtutor
        class(random_number_generator), intent(out) :: generator
        integer, intent(in) :: seed_ij, seed_kl
        if ((seed_ij < 0).or.(seed_ij > max_ij).or.(seed_kl < 0).or.(seed_kl > max_kl)) then
            call vtutor%error("The first random number seed must have a value between 0 and 31328 \n The &
               &second seed must have a value between 0 and 30081")
        end if
        call initialize_u(seed_ij, seed_kl, generator%u)
        generator%initialized = .true.
    end subroutine initialize

    pure subroutine initialize_u(seed_ij, seed_kl, u)
        integer, intent(in) :: seed_ij, seed_kl
        real(q), intent(out) :: u(97)
        integer i, j, k, l, m, ii, jj
        real(q) s, t
        !
        i = mod(seed_ij / 177, 177) + 2
        j = mod(seed_ij, 177) + 2
        k = mod(seed_kl / 169, 178) + 1
        l = mod(seed_kl, 169)
        !
        do ii = 1, size(u)
            s = 0.0_q
            t = 0.5_q
            do jj = 1, 24
                m = mod(mod(i * j, 179) * k, 179)
                i = j
                j = k
                k = m
                l = mod(53 * l + 1, 169)
                if (mod(l * m, 64) >= 32) then
                    s = s + t
                end if
                t = 0.5_q * t
            end do
            u(ii) = s
        end do
    end subroutine initialize_u

    subroutine random_real_scalar(generator, scalar)
        ! This is the random number generator proposed by George Marsaglia in
        ! Florida State University Report: FSU-SCRI-87-50
        ! It was slightly modified by F. James to produce an array of pseudorandom
        ! numbers.
        class(random_number_generator), intent(inout) :: generator
        real(q), intent(out) :: scalar
        integer, parameter :: seed_ij = 1802, seed_kl = 9373
        real(q) uni
#ifndef __NEC_TUNE__
        if (.not.generator%initialized) call generator%initialize(seed_ij, seed_kl)
        uni = generator%u(generator%i97) - generator%u(generator%j97)
        if (uni < 0.0_q) uni = uni + 1.0_q
        generator%u(generator%i97) = uni
        generator%i97 = generator%i97 - 1
        if (generator%i97 == 0) generator%i97 = 97
        generator%j97 = generator%j97 - 1
        if (generator%j97 == 0) generator%j97 = 97
        generator%c = generator%c - cd
        if (generator%c < 0.0_q) generator%c = generator%c + cm
        uni = uni - generator%c
        if (uni < 0.0_q) uni = uni + 1.0_q
        scalar = uni
#else // __NEC_TUNE__
        integer :: i, c_il, i97_l, j97_l

        if (generator%rand_cnt < generator%rand_blksz) then
!NEC_RP: if precomputed random number exists use one of those
           generator%rand_cnt = generator%rand_cnt + 1
           scalar = generator%rand_vec(generator%rand_cnt)
        else
!NEC_RP: if precomputed random number exists use one of those
           if (.not.generator%initialized) call generator%initialize(seed_ij, seed_kl)
!NEC$ ivdep
           do i = 1, generator%rand_blksz
              c_il = generator%c_i - i * generator%cd_i
              c_il = modulo(c_il, generator%cm_i)
              generator%c = c_il * generator%c_div
              i97_l = generator%i97-(i-1)
              if (i97_l <= 0) i97_l = i97_l + 97
              j97_l = generator%j97-(i-1)
              if (j97_l <= 0) j97_l = j97_l + 97
              uni = generator%u(i97_l) - generator%u(j97_l)
              if (uni < 0.0_q) uni = uni + 1.0_q
              generator%u(i97_l) = uni
              uni = uni - generator%c
              if (uni < 0.0_q) uni = uni + 1.0_q
              generator%rand_vec(i) = uni
           end do
           generator%i97 = generator%i97 - generator%rand_blksz
           if (generator%i97 <= 0) generator%i97 = generator%i97 + 97
           generator%j97 = generator%j97 - generator%rand_blksz
           if (generator%j97 <= 0) generator%j97 = generator%j97 + 97
           generator%c_i = generator%c_i - generator%rand_blksz * generator%cd_i
           generator%c_i = modulo(generator%c_i, generator%cm_i)
!NEC_RP: return first new random number
           generator%rand_cnt = 1
           scalar = generator%rand_vec(generator%rand_cnt)
        end if
#endif // __NEC_TUNE__
    end subroutine random_real_scalar

    subroutine random_real_array(generator, array)
        class(random_number_generator), intent(inout) :: generator
        real(q), intent(out) :: array(:)
        integer ii
        do ii = 1, size(array)
            call generator%random(array(ii))
        end do
    end subroutine random_real_array

    subroutine uniform_real_scalar(generator, scalar, center, width)
        ! This subroutine produces a normal distribution (Gaussian distribution)
        ! using the built in fortran random number generator
        ! which is initialized with a random seed for tbdyn
        class(random_number_generator), intent(inout) :: generator
        real(q), intent(out) :: scalar
        real(q), intent(in), optional :: center, width
        real(q), parameter :: threshold = 1.0e-8_q
        real(q) first, second
        !
        call generator%random(first)
        do
            ! the random number generator includes zero, which leads to a problem
            ! in the logarithm, hence we replace it with a nonzero value
            call generator%random(second)
            if (second > threshold) exit
        end do
        scalar = cos(twopi * first) * sqrt(-2.0_q * log(second))
        if (present(width)) scalar = scalar * width
        if (present(center)) scalar = scalar + center
    end subroutine uniform_real_scalar

    subroutine uniform_real_array(generator, array, center, width)
        class(random_number_generator), intent(inout) :: generator
        real(q), intent(out) :: array(:)
        real(q), intent(in), optional :: center, width
        integer ii
        do ii = 1, size(array)
            call generator%uniform(array(ii), center, width)
        end do
    end subroutine uniform_real_array

    subroutine uniform_complex_scalar(generator, scalar, center, width)
        class(random_number_generator), intent(inout) :: generator
        complex(q), intent(out) :: scalar
        real(q), intent(in), optional :: center, width
        real(q) re_part, im_part
        call generator%uniform(re_part, center, width)
        call generator%uniform(im_part, width=width) ! only real part shifted
        scalar = cmplx(re_part, im_part, kind=q)
    end subroutine uniform_complex_scalar

    subroutine uniform_complex_array(generator, array, center, width)
        class(random_number_generator), intent(inout) :: generator
        complex(q), intent(out) :: array(:)
        real(q), intent(in), optional :: center, width
        integer ii
        do ii = 1, size(array)
            call generator%uniform(array(ii), center, width)
        end do
    end subroutine uniform_complex_array

end module random_seeded
