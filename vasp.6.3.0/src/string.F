#include "symbol.inc"

module string

    use prec, only: q,qd
    USE base, ONLY:  TOREAL

    implicit none

    integer, parameter :: maxLength = 255
    real(q), parameter :: precision = 1e-3_q
    real(q), parameter :: minFloat = 10_q * precision, maxFloat = 1.0_q / precision

    !! If you pass this random string to the search routines they will ignore whitespace
    !! before newlines
    character(len=*), parameter :: any_newline = "esY8zP9FxW"

    !!
    !! Convert a numerical or logical variable into a string
    !!
    !! The routines provided in this module convert a given variable into a string.
    !! The focus of the routines is on stability and reasonable defaults with the
    !! target of quickly inspecting whether a variable exhibits reasonable values
    !! for example in the context of debugging or error messages.
    !!
    !! For the real and complex quantities, there is also the option to switch
    !! between scientific and floating point output and for the logicals you can
    !! print in the long form '.TRUE.' and '.FALSE.'. For all quantities, you
    !! can also give a format string to display the variable, but the default
    !! is used instead of the given format whenever the format would result in
    !! an error or an overflow.
    !!
    interface str
        module procedure logicalDefault, logicalCompact, logicalFormat, &
            integerFormat, integerArray, realDefault, realFormat, &
            realArray, real2dArray, complexDefault, complexFormat, &
            complexArray, complex2dArray
#if !defined(noQuadPrecision) || defined(qd_emulate)
        module procedure qrealDefault
#endif
    end interface str

    !! Find or count all instances of one or more search patterns.
    !!
    !! Initialize this type by setting *escape* and *ignore_between* if you want to
    !! use them.
    !! *escape* is a character that will lead to the directly preceeding character
    !! to be ignored, e.g., you can set escape="\" to ignore all "\n" when
    !! searching for "n".
    !! *ignore_between* is a character that will lead to all matches to be ignored
    !! between two instances of the pattern. A typical usage would be to exclude
    !! matches inside a string by setting ignore_between='"'.
    !! One additional feature is to pass the *any_newline* constant defined in this
    !! module as search string, which will also ignore all matches where a newline
    !! is escaped, but some spaces are between the escape and the newline character.
    !!
    !! The main routines of the class are *find_in* to obtain the indices of all
    !! instances of the search pattern and *count_in* to count the amount of all
    !! instances. Both routines can be used with either a scalar or an array search
    !! pattern.
    !!
    type, public :: indexer
        character(len=:), allocatable :: escape, ignore_between
    contains
        procedure find_in_scalar
        procedure find_in_array
        procedure count_in_scalar
        procedure count_in_array
        procedure interpret_special_characters
        generic :: find_in => find_in_scalar, find_in_array
        generic :: count_in => count_in_scalar, count_in_array
    end type indexer

    type count_in_status
        logical :: is_escape = .false., active_escape = .false.
        logical :: active_escape_with_whitespace = .false.
        logical :: is_ignore = .false., active_ignore = .false.
        logical is_newline, is_space
    end type count_in_status

    interface useScientificFormat
        module procedure useScientificFormatReal, useScientificFormatComplex
    end interface useScientificFormat

    private
    public str, string_from_file, string_from_unit, parseNewline, lowercase, uppercase, any_newline

contains

    function string_from_file(filename, ierr) result (content)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ierr
        character(len=:), allocatable :: content
        integer file_unit
        open(newunit=file_unit,file=filename,status='old',form='unformatted',access='stream',iostat=ierr)
        if (ierr /= 0) then
            content = ""
            return
        end if
        content = string_from_unit(file_unit, ierr)
        close(file_unit)
    end function string_from_file

    function string_from_unit(file_unit, ierr) result (content)
        integer, intent(in) :: file_unit
        integer, intent(out) :: ierr
        character(len=:), allocatable :: content
        character(len=40) access_
        integer filesize
        inquire(unit=file_unit, size=filesize)
        allocate(character(len=filesize) :: content)
        read(file_unit, pos=1, iostat=ierr) content
    end function string_from_unit

    pure subroutine find_in_scalar(this, text, search, indices)
        class(indexer), intent(in) :: this
        character(len=*), intent(in) :: text, search
        integer, allocatable, intent(out) :: indices(:)
        call this%find_in(text, [search], indices)
    end subroutine find_in_scalar

    pure subroutine find_in_array(this, text, search, indices)
        class(indexer), intent(in) :: this
        character(len=*), intent(in) :: text, search(:)
        integer, allocatable, intent(out) :: indices(:)
        integer num_elements
        !
        call this%count_in(text, search, num_elements)
        allocate(indices(num_elements))
        call this%count_in(text, search, num_elements, indices)
    end subroutine find_in_array

    pure subroutine count_in_scalar(this, text, search, num_elements, indices)
        class(indexer), intent(in) :: this
        character(len=*), intent(in) :: text, search
        integer, intent(out) :: num_elements
        integer, intent(out), optional :: indices(:)
        call this%count_in(text, [search], num_elements, indices)
    end subroutine count_in_scalar

    pure subroutine count_in_array(this, text, search, num_elements, indices)
        ! return number of elements found and optionally the indices
        class(indexer), intent(in) :: this
        character(len=*), intent(in) :: text, search(:)
        integer, intent(out) :: num_elements
        integer, intent(out), optional :: indices(:)
        type(count_in_status) status_
        integer first
        num_elements = 0
        do first = 1, len(text)
            call this%interpret_special_characters(text(first:), status_)
            if (found_any_matches(search, text(first:), status_)) then
                num_elements = num_elements + 1
                if (present(indices)) indices(num_elements) = first
            end if
            call update_active_regions(status_)
        end do
    end subroutine count_in_array

    pure subroutine interpret_special_characters(this, text, status_)
        class(indexer), intent(in) :: this
        character(len=*), intent(in) :: text
        type(count_in_status), intent(inout) :: status_
        if (allocated(this%ignore_between)) &
            status_%is_ignore = compare(this%ignore_between, text, status_%active_escape)
        if (allocated(this%escape)) &
            status_%is_escape = compare(this%escape, text, status_%active_escape)
        status_%is_newline = compare(new_line('n'), text, status_%active_escape_with_whitespace)
        status_%is_space = compare(" ", text, .false.)
    end subroutine interpret_special_characters

    pure logical function found_any_matches(search, text, status_) result (res)
        character(len=*), intent(in) :: text, search(:)
        type(count_in_status), intent(in) :: status_
        logical match(size(search))
        where (search /= any_newline)
            match = compare(search, text, status_%active_escape)
        elsewhere
            match = status_%is_newline
        end where
        res = any(match).and..not.status_%active_ignore
    end function found_any_matches

    elemental logical function compare(pattern, text, active_escape)
        character(len=*), intent(in) :: pattern, text
        logical, intent(in) :: active_escape
        integer length
        length = max(1, len_trim(pattern))
        if (length > len(text)) then
            compare = .false.
        else
            compare = (text(1:length) == trim(pattern)).and..not.active_escape
        end if
    end function compare

    pure subroutine update_active_regions(status_)
        type(count_in_status), intent(inout) :: status_
        if (status_%is_ignore) status_%active_ignore = .not.status_%active_ignore
        status_%active_escape_with_whitespace = status_%is_escape .or. &
            (status_%active_escape_with_whitespace.and.status_%is_space)
        status_%active_escape = status_%is_escape
    end subroutine update_active_regions

    pure function lowercase(str) result (res)
        character(len=*), intent(in) :: str
        character(len=len(str)) :: res
        integer, parameter :: shift = ichar('A') - ichar('a')
        integer i, ich
        res = str
        do i = 1,len(str)
            ich = ichar(res(i:i))
            if ((ich >= ichar('A')).and.(ich <= ichar('Z'))) res(i:i) = achar(ich - shift)
        end do
    end function lowercase

    pure function uppercase(str) result (res)
        character(len=*), intent(in) :: str
        character(len=len(str)) :: res
        integer, parameter :: shift = ichar('A') - ichar('a')
        integer i, ich
        res = str
        do i = 1,len(str)
            ich = ichar(res(i:i))
            if ((ich >= ichar('a')).and.(ich <= ichar('z'))) res(i:i) = achar(ich + shift)
        end do
    end function uppercase

    !! Replace \n with new_line('n')
    pure function parseNewline(str) result (res)
        character(len=*), intent(in) :: str
        character(len=:), allocatable :: res, tmp
        character(*), parameter :: escapeCharacter = '\n'
        integer pos
        res = str
        ! some compiler settings replace escape characters automatically
        if (escapeCharacter == new_line('n')) return
        do
            pos = index(res, escapeCharacter)
            if (pos == 0) return
            ! temporary variable seems necessary for older gfortran
            tmp = res(1:pos-1) // new_line('n') // res(pos+2:len(res))
            res = tmp
        end do
    end function parseNewline

    pure function logicalDefault(lvar) result (res)
        logical, intent(in) :: lvar
        character(len=:), allocatable :: res
        if (lvar) then
            res = 'T'
        else
            res = 'F'
        end if
    end function logicalDefault

    pure function logicalCompact(lvar, compact) result (res)
        logical, intent(in) :: lvar
        logical, intent(in) :: compact
        character(len=:), allocatable :: res
        if (compact) then
            res = str(lvar)
        else
            if (lvar) then
                res = '.TRUE.'
            else
                res = '.FALSE.'
            end if
        end if
    end function logicalCompact

    pure function logicalFormat(lvar, formatString) result (res)
        logical, intent(in) :: lvar
        character(len=*), intent(in) :: formatString
        character(len=:), allocatable :: res
        character(len=maxLength) tmp
        integer ierr
        write(tmp, formatString, iostat=ierr) lvar
        if (ierr == 0) then
            res = trim(tmp)
        else
            res = str(lvar)
        end if
    end function logicalFormat

    pure function integerFormat(ivar, formatString) result (res)
        integer, intent(in) :: ivar
        character(len=*), intent(in), optional :: formatString
        character(len=*), parameter :: defaultFormat = '(i0)'
        character(len=:), allocatable :: res
        character(len=maxLength) tmp
        integer ierr
        if (present(formatString)) &
            write(tmp, formatString, iostat=ierr) ivar
        if (.not.present(formatString).or.ierr /= 0.or.tmp(1:1) == '*') &
            write(tmp, defaultFormat) ivar
        res = trim(tmp)
    end function integerFormat

    pure function integerArray(iarr, formatString) result (res)
        integer, intent(in) :: iarr(:)
        character(len=*), intent(in), optional :: formatString
        character(len=:), allocatable :: res
        integer ii
        res = ""
        do ii = 1, size(iarr)
            if (ii > 1) res = res // " "
            res = res // str(iarr(ii), formatString)
        end do
    end function integerArray

    pure function realDefault(rvar, scientific) result (res)
        real(q), intent(in) :: rvar
        logical, intent(in), optional :: scientific
        character(len=:), allocatable :: res
        logical scientific_
        call convertSpecialInput(rvar, res)
        if (allocated(res)) return
        scientific_ = useScientificFormat(rvar, scientific)
        if (scientific_) then
            res = realScientific(rvar)
        else
            res = realFloatingPoint(rvar)
        end if
    end function realDefault

    pure function realFormat(rvar, formatString) result (res)
        real(q), intent(in) :: rvar
        character(len=*), intent(in) :: formatString
        character(len=:), allocatable :: res
        character(len=maxLength) tmp
        integer ierr
        write(tmp, formatString, iostat=ierr) rvar
        if (ierr == 0.and.tmp(1:1) /= '*') then
            res = trim(tmp)
        else
            res = str(rvar)
        end if
    end function realFormat

    pure function realArray(rarr, scientific) result (res)
        real(q), intent(in) :: rarr(:)
        logical, intent(in), optional :: scientific
        character(len=:), allocatable :: res
        logical scientific_
        integer ii
        scientific_ = any(useScientificFormat(rarr, scientific))
        res = ""
        do ii = 1, size(rarr)
            if (ii > 1) res = res // " "
            res = res // str(rarr(ii), scientific_)
        end do
    end function realArray

    pure function real2dArray(rarr, scientific) result (res)
        real(q), intent(in) :: rarr(:,:)
        logical, intent(in), optional :: scientific
        character(len=:), allocatable :: res
        logical scientific_
        integer ii, jj
        scientific_ = any(useScientificFormat(rarr, scientific))
        res = ""
        do jj = 1, size(rarr, 2)
            do ii = 1, size(rarr, 1)
                if (ii /= 1.or.jj /= 1) res = res // " "
                res = res // str(rarr(ii, jj), scientific_)
            end do
        end do
    end function real2dArray

    pure subroutine convertSpecialInput(rvar, res)
        use ieee_arithmetic, only: ieee_is_normal
        real(q), intent(in) :: rvar
        character(len=:), allocatable, intent(out) :: res
        if (.not.ieee_is_normal(rvar)) then
            call convertAbnormalInput(rvar, res)
        else if (abs(rvar) < tiny(rvar)) then
            res = '0'
        end if
    end subroutine convertSpecialInput

    pure subroutine convertAbnormalInput(rvar, res)
        use ieee_arithmetic, only: ieee_is_nan, ieee_is_negative
        real(q), intent(in) :: rvar
        character(len=:), allocatable, intent(out) :: res
        if (ieee_is_nan(rvar)) then
            res = 'NaN'
        else if (ieee_is_negative(rvar)) then
            res = '-Infinity'
        else
            res = 'Infinity'
        end if
    end subroutine convertAbnormalInput

    elemental logical function useScientificFormatReal(rvar, scientific) result (scientific_)
        real(q), intent(in) :: rvar
        logical, intent(in), optional :: scientific
        real(q) absValue
        absValue = abs(rvar)
        if (present(scientific)) then
            scientific_ = scientific
        else
            scientific_ = absValue > dzero .and. &
                (minFloat > abs(rvar) .or. abs(rvar) > maxFloat)
        end if
    end function useScientificFormatReal

    pure function realScientific(rvar) result (res)
        real(q), intent(in) :: rvar
        character(len=:), allocatable :: res, expString
        character(len=maxLength) tmp
        integer exponent
        write(tmp, '(es11.3e3)') rvar
        tmp = trim(adjustl(tmp))
        exponent = index(tmp, 'E')
        if (exponent > 0) then
            expString = tmp(exponent+1:len_trim(tmp))
        else
            ! some compilers always use two digits for the exponential part
            exponent = scan(tmp, '-+', back = .true.) - 1
            expString = tmp(exponent+1:len_trim(tmp))
        end if
        expString = removeRedundantZeros(expString, backwards = .false.)
        if (expString /= '') expString = 'e' // expString
        res = removeRedundantZeros(tmp(1:exponent), backwards = .true.) // expString
    end function realScientific

    pure function realFloatingPoint(rvar) result (res)
        real(q), intent(in) :: rvar
        character(len=:), allocatable :: res, formatString
        character(len=maxLength) tmp
        integer lastNonZero
        ! factor two to account for rounding
        if (2.0_q * abs(rvar) < precision) then
            res = '0'
        else
            formatString = '(f' // str(maxLength) // '.3)'
            write(tmp, formatString) rvar
            res = removeRedundantZeros(trim(adjustl(tmp)), backwards = .true.)
        end if
    end function realFloatingPoint

    pure function removeRedundantZeros(stringOfNumber, backwards) result (res)
        character(len=*), intent(in) :: stringOfNumber
        logical, intent(in) :: backwards
        character(len=:), allocatable :: res
        character(*), parameter :: nonZero = "123456789."
        integer firstNonZero, signOffset
        signOffset = 0
        if (stringOfNumber(1:1) == '-') signOffset = 1
        firstNonZero = scan(stringOfNumber(1+signOffset:), nonZero, backwards)
        if (firstNonZero == 0) then
            res = ''
        else
            firstNonZero = firstNonZero + signOffset
            if (stringOfNumber(firstNonZero:firstNonZero) == '.') &
                firstNonZero = firstNonZero - 1
            if (backwards) then
                res = stringOfNumber(:firstNonZero)
            else
                res = stringOfNumber(:signOffset) // stringOfNumber(firstNonZero:)
            end if
        end if
    end function removeRedundantZeros

    pure function complexDefault(cvar, scientific) result (res)
        complex(q), intent(in) :: cvar
        logical, intent(in), optional :: scientific
        character(len=:), allocatable :: res
        logical scientific_
        scientific_ = useScientificFormat(cvar, scientific)
        if (significantImagPart(cvar, scientific_)) then
            res = convertRealAndImagPart(cvar, scientific_)
        else
            res = str(real(cvar), scientific_)
        end if
    end function complexDefault

    pure function convertRealAndImagPart(cvar, scientific) result (res)
        complex(q), intent(in) :: cvar
        logical, intent(in) :: scientific
        character(len=:), allocatable :: res
        res = '(' // str(real(cvar), scientific) // ', ' &
            // str(aimag(cvar), scientific) // ')'
    end function convertRealAndImagPart

    pure function complexFormat(cvar, formatString) result (res)
        complex(q), intent(in) :: cvar
        character(len=*), intent(in) :: formatString
        character(len=:), allocatable :: res
        character(len=maxLength) tmp
        integer ierr
        write(tmp, formatString, iostat=ierr) cvar
        if (ierr == 0.and.tmp(1:1) /= '*') then
            res = trim(tmp)
        else
            res = str(cvar)
        end if
    end function complexFormat

    pure function complexArray(carr, scientific) result (res)
        complex(q), intent(in) :: carr(:)
        logical, intent(in), optional :: scientific
        character(len=:), allocatable :: res
        logical sigImagPart, scientific_
        integer ii
        scientific_ = any(useScientificFormat(carr, scientific))
        sigImagPart = any(significantImagPart(carr, scientific_))
        if (sigImagPart) then
            res = ""
            do ii = 1, size(carr)
                if (ii > 1) res = res // " "
                res = res // convertRealAndImagPart(carr(ii), scientific_)
            end do
        else
            res = str(real(carr), scientific_)
        end if
    end function complexArray

    pure function complex2dArray(carr, scientific) result (res)
        complex(q), intent(in) :: carr(:,:)
        logical, intent(in), optional :: scientific
        character(len=:), allocatable :: res
        logical sigImagPart, scientific_
        integer ii, jj
        scientific_ = any(useScientificFormat(carr, scientific))
        sigImagPart = any(significantImagPart(carr, scientific_))
        if (sigImagPart) then
            res = ""
            do jj = 1, size(carr, 2)
                do ii = 1, size(carr, 1)
                    if (ii /= 1.or.jj /= 1) res = res // " "
                    res = res // convertRealAndImagPart(carr(ii, jj), scientific_)
                end do
            end do
        else
            res = str(real(carr), scientific_)
        end if
    end function complex2dArray

    elemental logical function significantImagPart(cvar, scientific)
        use ieee_arithmetic, only: ieee_is_normal
        complex(q), intent(in) :: cvar
        logical, intent(in) :: scientific
        significantImagPart = scientific .or. (abs(aimag(cvar)) >= precision) &
            .or..not.ieee_is_normal(aimag(cvar))
    end function significantImagPart

    elemental logical function useScientificFormatComplex(cvar, scientific) result (scientific_)
        use ieee_arithmetic, only: ieee_is_normal
        complex(q), intent(in) :: cvar
        logical, intent(in), optional :: scientific
        real(q) testValue
        if (present(scientific)) then
            scientific_ = scientific
            return
        end if
        associate(re => real(cvar), im => aimag(cvar))
            if (ieee_is_normal(re).and.ieee_is_normal(im)) then
                testValue = max(abs(re), abs(im))
            else if (ieee_is_normal(re)) then
                testValue = re
            else if (ieee_is_normal(im)) then
                testValue = im
            else
                testValue = done
            end if
            scientific_ = useScientificFormat(testValue)
        end associate
    end function useScientificFormatComplex

#if !defined(noQuadPrecision) || defined(qd_emulate)
    function qrealDefault(qrvar, scientific) result (res)
#ifdef qd_emulate
        USE qdmodule
#endif
        QDPREAL, intent(in) :: qrvar
        logical, intent(in), optional :: scientific
        character(len=:), allocatable :: res
        res = str(TOREAL(qrvar), scientific)
    end function qrealDefault
#endif

end module string
