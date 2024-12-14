#include "symbol.inc"
module tutor
    use prec, only: q
    use iso_fortran_env, only: output_unit, error_unit

    implicit none

    enum, bind(c)
        enumerator :: never, isBug, isError, isAlert, isWarning, isAdvice
    end enum
!$acc declare copyin(isBug, isError)

    enum, bind(c)
        ! NOTE: add new error codes to the end to avoid screwing up the tests
        enumerator :: UndefinedCode
        enumerator :: RealSpaceWithoutOptimization, SplineInterpolateProjectors, &
            VASP_4_4, NoRealSpaceAndYouCould, NoRealSpaceAndYouShould, IALGO8, &
            RealSpaceNoMoreRecommended, WrongOptimizationRealSpace, LongLattice, &
            DifferentXCGradTypes, NumberOfElectrons, NBANDSchanged, &
            DifferentLDAXCTypes, DifferentXCFiles, AM05spin, NoLDAU, &
            DoubleCounting, NoWannier, LREALA, DifferentSlaterXCTypes, &
            PartialDOS, NoOptics, Position, CenterOfMassDrift, EnforcedLDA, &
            DipolEpsilon, DipolCubic, ALGO_A_ISMEAR, ModelGW, HFUSPP, &
            LHFLinearResponse, LinearReponsePOSCAR, LinearResponseLREAL, &
            FockForce, HF_NPAR, NPAREfficiency, NPAR, LHFOneCenter, HFParam, &
            OEPdouble, NMAXFOCKAE, ISYM2, LSPECTRALinefficient, LSPECTRALGWset, &
            LPOTOK, XiSingular, AEXXzero, GWoptics, GW_HFoptics, ENCUTGW_LHF, &
            DeprENCUTFOCK, PRECFOCK, KPOINTS_HF, KPOINTS_PEAD, KPOINTS_INTER, &
            GWkweight, GWnoreal, ACFDTRnotSupported, PEADnoVirtuals, PEADlargeField, &
            PEAD_NCOREnot1, PEADPolNoMetal, PEADEpsNoMetal, &
            LINTET, GRDSHFT_NOTET, GammaK, KPOINTS_IO, NoInversion, &
            FFTGridIsNotSufficient, MetaGGAandForces, MetaGGAResponse, &
            SmallCSHIFT, TetNoMetal, POTIMlarge, GRDSHFTnotInIBZ, GRDSHFTK, &
            NoOcc, NoOddGrid, GammaC, ISYM_MD, BSEantires, BSEconjg, BSEncv, &
            Wmissing, AuxWAVECAR, SpinSpiral, HighestBandsOcccupied, LVTOT, &
            vdWKlimes, NoKPAR, FERWEOptics, OmegaPar, TauPar, &
            TauParCPU, ImagGridParForLTMP2, chiGG_NOmega, OmegaGridOMEGAMIN, WannNotFound, &
            TargetNotSet, NTARGET_STATES, NCRPALOW, CRPADisentangle, LCAO_INCAR, &
            SHMEM_CRPA, OmegaGridOmegaTL, NOmega32, OmegaGrid140, NUPDOWN, &
            NPAR_IALGO8, SHMEMError, KProjection, Gamma, Gamma2, DifferentRelXCTypes, &
            Alternant, LBONEandLZORA, RandomSeed, PEADnowIsMetal, &
            NBANDlessNWAN, ELPHtesting, GDESMatBlockRed, GDESMatSize, &
            NOMEGA_DUMPtooLarge, SigmaTMBPT, ALLOC_WAVEFUN, &
            FermiFiniteTemp, BlochBasisUsed, TargetOut, NCRPA_BANDSout, &
            LOCPROJnotSet, ALGO_GW, ALGO_GW0, ALGO_G0W0, ALGO_scGW0, &
            ALGO_scGW, POTCARmissing, SymprecDependence, &
            ALLOC_RESPONSE, ALLOC_RESPONSE_DISTR, MAXMEM_LOW, AUTOSET_MEM, &
            RPAFORCES_COMPAT, RPAENCUTGWSOFT, VMBJ_BrentInitialPositions, VMBJ_TooManyIterations, &
            ORTHONORMALIZATION_FAILED
        enumerator :: numTopic
    end enum
!$acc declare copyin(VMBJ_BrentInitialPositions, VMBJ_TooManyIterations)

    type logger
        integer :: unitOut = output_unit, unitErr = error_unit, unitBug = error_unit
        integer :: stopOn = isError
        logical :: discardOutputToNegativeUnit = .false.
    contains
        procedure :: bug => logBug
        procedure :: error => logError
        procedure :: alert => logAlert
        procedure :: warning => logWarning
        procedure :: advice => logAdvice
        procedure shouldStop
        procedure stopCode
        procedure write
        procedure printMessageToUnits
    end type logger

    ! Note that the global logger vtutor needs to be initialized before its first use
    ! otherwise it will use the default output and error unit
    type(logger) :: vtutor

    character(*), parameter :: headerBug = '&
        &                _     ____    _    _    _____     _ \n&
        &               | |   |  _ \  | |  | |  / ____|   | |\n&
        &               | |   | |_) | | |  | | | |  __    | |\n&
        &               |_|   |  _ <  | |  | | | | |_ |   |_|\n&
        &                _    | |_) | | |__| | | |__| |    _ \n&
        &               (_)   |____/   \____/   \_____|   (_)\n\n'
    character(*), parameter :: footerBug = '\n\n&
        &If you are not a developer, you should not encounter this problem.\n&
        &Please submit a bug report.\n\n'

    character(*), parameter :: headerError = '\n&
        &EEEEEEE  RRRRRR   RRRRRR   OOOOOOO  RRRRRR      ###     ###     ###\n&
        &E        R     R  R     R  O     O  R     R     ###     ###     ###\n&
        &E        R     R  R     R  O     O  R     R     ###     ###     ###\n&
        &EEEEE    RRRRRR   RRRRRR   O     O  RRRRRR       #       #       # \n&
        &E        R   R    R   R    O     O  R   R                          \n&
        &E        R    R   R    R   O     O  R    R      ###     ###     ###\n&
        &EEEEEEE  R     R  R     R  OOOOOOO  R     R     ###     ###     ###\n\n'
    character(*), parameter :: footerError = '\n\n&
        &  ---->  I REFUSE TO CONTINUE WITH THIS SICK JOB ... BYE!!! <----\n\n'

    character(*), parameter :: headerAlert = '\n&
        &      W    W    AA    RRRRR   N    N  II  N    N   GGGG   !!!\n&
        &      W    W   A  A   R    R  NN   N  II  NN   N  G    G  !!!\n&
        &      W    W  A    A  R    R  N N  N  II  N N  N  G       !!!\n&
        &      W WW W  AAAAAA  RRRRR   N  N N  II  N  N N  G  GGG   ! \n&
        &      WW  WW  A    A  R   R   N   NN  II  N   NN  G    G     \n&
        &      W    W  A    A  R    R  N    N  II  N    N   GGGG   !!!\n\n'
    character(*), parameter :: footerAlert = '\n\n'

    character(*), parameter :: headerWarning = ' WARNING: '
    character(*), parameter :: footerWarning = ''

    character(*), parameter :: headerAdvice = '\n&
        &          ----> ADVICE to this user running VASP <----\n\n'
    character(*), parameter :: footerAdvice = '\n\n'


    character(*), parameter :: horizontalLineSegment = '-', verticalLineSegment = '|'

    type messageFormatter
        integer :: tabWidth = 0, padding = 5, totalWidth = 77
        logical :: useBox = .true.
    contains
        procedure formatMessage
        procedure horizontalLineOfBox
        procedure wrapMessage
        procedure determineWidth
        procedure formatLine
    end type messageFormatter

    type argument
        integer, allocatable :: ival(:)
        logical, allocatable :: lval(:)
        real(q), allocatable :: rval(:)
    end type argument

    type messageSplitter
        integer :: first = 1, last = 0
        character(len=:), allocatable :: message
    contains
        procedure finished
        procedure nextPart
    end type messageSplitter

    private
    public logger, argument, vtutor, never, isBug, isError, isAlert, isWarning, isAdvice, &
        numTopic, RealSpaceWithoutOptimization, SplineInterpolateProjectors, VASP_4_4, &
        NoRealSpaceAndYouCould, NoRealSpaceAndYouShould, IALGO8, DifferentXCFiles, &
        RealSpaceNoMoreRecommended, WrongOptimizationRealSpace, LongLattice, &
        DifferentXCGradTypes, NumberOfElectrons, NBANDSchanged, DifferentLDAXCTypes, &
        AM05spin, NoLDAU, DoubleCounting, NoWannier, LREALA, DifferentSlaterXCTypes, &
        PartialDOS, NoOptics, Position, CenterOfMassDrift, EnforcedLDA, &
        DipolEpsilon, DipolCubic, ALGO_A_ISMEAR, ModelGW, HFUSPP, ISYM_MD, &
        LHFLinearResponse, LinearReponsePOSCAR, LinearResponseLREAL, FockForce, &
        HF_NPAR, NPAREfficiency, NPAR, LHFOneCenter, HFParam, OEPdouble, NMAXFOCKAE, &
        ISYM2, LSPECTRALinefficient, LSPECTRALGWset, LPOTOK, XiSingular, AEXXzero, &
        GWoptics, GW_HFoptics, ENCUTGW_LHF, DeprENCUTFOCK, PRECFOCK, KPOINTS_HF, &
        KPOINTS_PEAD, KPOINTS_INTER, GWkweight, GWnoreal, ACFDTRnotSupported, &
        PEADnoVirtuals, PEADlargeField, PEAD_NCOREnot1, PEADPolNoMetal, PEADEpsNoMetal, &
        PEADnowIsMetal, LINTET, GRDSHFT_NOTET, GammaK, KPOINTS_IO, NoInversion, &
        FFTGridIsNotSufficient, MetaGGAandForces, MetaGGAResponse, SmallCSHIFT, &
        TetNoMetal, POTIMlarge, GRDSHFTnotInIBZ, GRDSHFTK, NoOcc, NoOddGrid, &
        GammaC, BSEantires, BSEconjg, BSEncv, Wmissing, AuxWAVECAR, SpinSpiral, &
        HighestBandsOcccupied, LVTOT, vdWKlimes, NoKPAR, FERWEOptics, &
        OmegaPar, TauPar, TauParCPU, ImagGridParForLTMP2, chiGG_NOmega, &
        OmegaGridOMEGAMIN, WannNotFound, TargetNotSet, NTARGET_STATES, NCRPALOW, &
        CRPADisentangle, LCAO_INCAR, SHMEM_CRPA, OmegaGridOmegaTL, NOmega32, &
        OmegaGrid140, NUPDOWN, NPAR_IALGO8, SHMEMError, KProjection, Gamma, Gamma2, &
        DifferentRelXCTypes, Alternant, LBONEandLZORA, RandomSeed, &
        NBANDlessNWAN, ELPHtesting, GDESMatBlockRed, GDESMatSize, &
        NOMEGA_DUMPtooLarge, SigmaTMBPT, FermiFiniteTemp, &
        BlochBasisUsed, TargetOut, NCRPA_BANDSout, LOCPROJnotSet, ALGO_GW, &
        ALGO_GW0, ALGO_G0W0, ALGO_scGW0, ALGO_scGW,POTCARmissing, SymprecDependence, &
        vtutor_error, ALLOC_RESPONSE, ALLOC_RESPONSE_DISTR, MAXMEM_LOW, AUTOSET_MEM, &
        ALLOC_WAVEFUN, RPAFORCES_COMPAT, RPAENCUTGWSOFT, deactivateVtutorOutputForTest, vtutor_bug, &
        VMBJ_BrentInitialPositions, VMBJ_TooManyIterations, ORTHONORMALIZATION_FAILED

#ifdef MPI
    include 'mpif.h'
#endif

contains

    subroutine vtutor_bug(message, filename, linenumber)
        character(len=*), intent(in) :: message, filename
        integer, intent(in) :: linenumber
        call vtutor%bug(message, filename, linenumber)
    end subroutine vtutor_bug

    subroutine vtutor_error(message)
        character(len=*), intent(in) :: message
        call vtutor%error(message)
    end subroutine vtutor_error

    subroutine logBug(this, message, filename, linenumber)
        class(logger), intent(in) :: this
        character(len=*), intent(in) :: message, filename
        integer, intent(in) :: linenumber
        type(messageFormatter) formatter
        character(len=:), allocatable :: whereInCode, formattedMessage
        whereInCode = createWhereInCodeInfo(filename, linenumber)
        formattedMessage = formatter%formatMessage(headerBug // whereInCode // message // footerBug)
        call this%printMessageToUnits(formattedMessage, [this%unitBug])
        if (this%stopOn >= isBug) call this%stopCode(isBug)
    end subroutine logBug

    function createWhereInCodeInfo(filename, linenumber) result(res)
        use string, only: str
        character(len=*), intent(in) :: filename
        integer, intent(in) :: linenumber
        character(len=:), allocatable :: res
        if (linenumber > 0) then
           res = "internal error in: " // trim(filename) // "  at line: " // str(linenumber) // "\n\n"
        else
           res = "internal error in: " // trim(filename) // "  at line: " // "UnknownLine" // "\n\n"
        endif
    end function createWhereInCodeInfo

    subroutine logError(this, message, dont_stop)
        class(logger), intent(in) :: this
        character(len=*), intent(in) :: message
        logical, intent(in), optional :: dont_stop
        type(messageFormatter) formatter
        character(len=:), allocatable :: formattedMessage
        formattedMessage = formatter%formatMessage(headerError // message // footerError)
        call this%printMessageToUnits(formattedMessage, [this%unitOut, this%unitErr])
        if (this%shouldStop(isError, dont_stop)) call this%stopCode(isError)
    end subroutine logError

    subroutine logAlert(this, message, dont_stop)
        class(logger), intent(in) :: this
        character(len=*), intent(in) :: message
        logical, intent(in), optional :: dont_stop
        type(messageFormatter) formatter
        character(len=:), allocatable :: formattedMessage
        formattedMessage = formatter%formatMessage(headerAlert // message // footerAlert)
        call this%printMessageToUnits(formattedMessage, [this%unitOut, this%unitErr])
        if (this%shouldStop(isAlert, dont_stop)) call this%stopCode(isAlert)
    end subroutine logAlert

    subroutine logWarning(this, message, dont_stop)
        class(logger), intent(in) :: this
        character(len=*), intent(in) :: message
        logical, intent(in), optional :: dont_stop
        type(messageFormatter) formatter
        character(len=:), allocatable :: formattedMessage
        formatter = messageFormatter(tabWidth = 5, padding = 0, totalWidth = 78, &
            useBox = .false.)
        formattedMessage = formatter%formatMessage(headerWarning // message // footerWarning)
        call this%printMessageToUnits(formattedMessage, [this%unitOut])
        if (this%shouldStop(isWarning, dont_stop)) call this%stopCode(isWarning)
    end subroutine logWarning

    subroutine logAdvice(this, message, dont_stop)
        class(logger), intent(in) :: this
        character(len=*), intent(in) :: message
        logical, intent(in), optional :: dont_stop
        type(messageFormatter) formatter
        character(len=:), allocatable :: formattedMessage
        formattedMessage = formatter%formatMessage(headerAdvice // message // footerAdvice)
        call this%printMessageToUnits(formattedMessage, [this%unitOut, this%unitErr])
        if (this%shouldStop(isAdvice, dont_stop)) call this%stopCode(isAdvice)
    end subroutine logAdvice

    function formatMessage(formatter, message) result(res)
        use string, only: parseNewLine
        class(messageFormatter), intent(in) :: formatter
        character(len=*), intent(in) :: message
        character(len=:), allocatable :: res
        res = formatter%horizontalLineOfBox()
        res = res // formatter%wrapMessage(parseNewline(message))
        res = res // formatter%horizontalLineOfBox()
    end function formatMessage

    function horizontalLineOfBox(formatter) result (res)
        class(messageFormatter), intent(in) :: formatter
        character(len=:), allocatable :: res
        if (formatter%useBox) then
            res = ' ' // repeat(horizontalLineSegment, formatter%totalWidth) // new_line('n')
        else
            res = ''
        end if
    end function horizontalLineOfBox

    function wrapMessage(formatter, message) result(res)
        class(messageFormatter), intent(in) :: formatter
        character(len=*), intent(in) :: message
        character(len=:), allocatable :: res, part
        type(messageSplitter) splitter
        logical firstLine
        integer width
        res = ''
        firstLine = .true.
        allocate(character(len=len(message)) :: splitter%message)
        splitter%message = message
        do while (.not.splitter%finished())
            width = formatter%determineWidth(firstLine)
            part = splitter%nextPart(width)
            res = res // formatter%formatLine(firstLine, part)
            firstLine = .false.
        end do
    end function wrapMessage

    logical function finished(splitter)
        class(messageSplitter), intent(in) :: splitter
        finished = splitter%first > len(splitter%message)
    end function finished

    function nextPart(splitter, width) result (res)
        class(messageSplitter), intent(inout) :: splitter
        integer, intent(in) :: width
        character(len=:), allocatable :: res
        integer, parameter :: skipSpace = 2
        integer ii
        character char
        do ii = splitter%first, splitter%first + width
            if (ii > len(splitter%message)) then
                splitter%last = len(splitter%message)
                exit
            end if
            char = splitter%message(ii:ii)
            if (char == ' '.or.char == new_line('n')) splitter%last = ii - 1
            if (char == new_line('n')) exit
        end do
        res = splitter%message(splitter%first:splitter%last)
        splitter%first = splitter%last + skipSpace
    end function nextPart

    integer function determineWidth(formatter, firstLine) result(width)
        class(messageFormatter), intent(in) :: formatter
        logical, intent(in) :: firstLine
        width = formatter%totalWidth - 2 * formatter%padding
        if (.not.firstLine) width = width - formatter%tabWidth
    end function determineWidth

    function formatLine(formatter, firstLine, part) result(res)
        class(messageFormatter), intent(in) :: formatter
        logical, intent(in) :: firstLine
        character(len=*), intent(in) :: part
        character(len=:), allocatable :: res
        integer numSpace
        if (formatter%useBox) then
            res = repeat(' ', formatter%padding) // part
            numSpace = formatter%totalWidth - len(res)
            res = verticalLineSegment // res // repeat(' ', numSpace) // verticalLineSegment
        else
            res = part
            if (.not.firstLine) res = repeat(' ', formatter%tabWidth) // res
        end if
        res = res // new_line('n')
    end function formatLine

    subroutine printMessageToUnits(this, message, units)
        class(logger), intent(in) :: this
        character(len=*), intent(in) :: message
        integer, intent(in) :: units(:)
        integer iunit
        do iunit = 1, size(units)
            if (this%discardOutputToNegativeUnit .and. units(iunit) < 0) cycle
            if (unitNotOpen(units(iunit))) cycle
            !$omp critical(omp_wrt_stdout)
            write(units(iunit), '(a)') message
            flush(units(iunit))
            !$omp end critical(omp_wrt_stdout)
        end do
    end subroutine printMessageToUnits

    logical function unitNotOpen(iunit)
        integer, intent(in) :: iunit
        logical unitOpen
        inquire(iunit, opened = unitOpen)
        unitNotOpen = .not.unitOpen
    end function unitNotOpen
#ifdef MPI
    subroutine stopMpi(errCode)
        integer, intent(in), optional :: errCode
        integer errCode_, request, ierr
        logical complete
        errCode_ = 0
        if (present(errCode)) errCode_ = errCode
        if (errCode_ /= 0) then
            call MPI_ibarrier(MPI_comm_world, request, ierr)
            complete = tryWaitingForOtherCPU(request, ierr)
            if (ierr /= MPI_success.or..not.complete) &
                call MPI_abort(MPI_comm_world, 1, ierr)
        end if
        call MPI_finalize(ierr)
        if (ierr /= MPI_success) &
            call MPI_abort(MPI_comm_world, 1, ierr)
    end subroutine stopMpi

    logical function tryWaitingForOtherCPU(request, ierr) result(complete)
        integer, intent(in) :: request
        integer, intent(inout) :: ierr
        integer try
        if (ierr /= MPI_success) return
        do try = 1, 30
            call MPI_test(request, complete, MPI_STATUS_IGNORE, ierr)
            if (ierr /= MPI_success.or.complete) exit
            call sleep(1)
        end do
    end function tryWaitingForOtherCPU
#endif

    logical function shouldStop(this, severity, dont_stop)
        class(logger), intent(in) :: this
        integer, intent(in) :: severity
        logical, intent(in), optional :: dont_stop
        shouldStop = this%stopOn >= severity
        if (present(dont_stop)) then
            shouldStop = .not.dont_stop.and.shouldStop
        end if
    end function shouldStop

    subroutine stopCode(this, errCode)
        class(logger), intent(in) :: this
        integer, intent(in), optional :: errCode
        logical errStop
        errStop = .false.
        if (present(errCode)) errStop = errCode /= never
        CALLMPI( stopMpi(errCode) )
        if (errStop) stop 1
        stop
    end subroutine stopCode

    subroutine write(this, severity, topic, arg, filename, linenumber)
        class(logger), intent(in) :: this
        integer, intent(in) :: severity, topic
        type(argument), intent(in), optional :: arg
        character(len=*), intent(in), optional :: filename
        integer, intent(in), optional :: linenumber
        character(len=:), allocatable :: message, fnam
        integer line
        message = generateMessage(topic, arg)
        select case (severity)
        case (isBug)
            fnam = "UnkownFile" ; if (present(filename)) fnam = filename
            line = -1 ; if (present(linenumber)) line = linenumber
            call this%bug(message, fnam, line)
        case (isError)
            call this%error(message)
        case (isAlert)
            call this%alert(message)
        case (isWarning)
            call this%warning(message)
        case (isAdvice)
            call this%advice(message)
        end select
    end subroutine write

    function generateMessage(topic, arg) result(res)
        use string, only: str
        integer, intent(in) :: topic
        type(argument), intent(in), optional :: arg
        character(len=:), allocatable :: res, work
        select case (topic)
        case(RealSpaceWithoutOptimization)
            res = 'The real-space projection scheme for the treatment of the non-local &
                &pseudopotentials has been switched on but on file POTCAR I have &
                &not found any entries, signalling me that you have never performed a &
                &real-space optimization of all non-local projectors! BE WARNED that &
                &a calculation using the real-space projection scheme together with &
                &non-local projectors which were not real-space optimized (according &
                &to the plane-wave cutoff used in this calculation) might give more &
                &or less inaccurate results. I hope you know what you are doing.'
        case(SplineInterpolateProjectors)
            res = 'The projection operators are constructed using a spline interpolation. &
                &This is needed for the calculation of chemical shifts by means of &
                &linear response but yields slightly different total energies compared &
                &to the default setup of the projectors (NLSPLINE=.FALSE.). '
        case(VASP_4_4)
            res = 'You are running vasp.4.5 in the vasp.4.4 compatibility mode. &
                &vasp.4.5 has some numerical improvements, which are not applied in this &
                &mode (charges at unbalanced lattice vectors are no longer zeroed, &
                &PAW augmentation charges are integrated more accurately).'
        case(NoRealSpaceAndYouCould)
            res = 'You have a (more or less) ''large supercell'' and for larger cells &
                &it might be more efficient to use real-space projection operators. &
                &Therefore, try LREAL= Auto in the INCAR file.\n&
                &Mind: If you want to do very accurate calculations, keep the &
                &reciprocal projection scheme (i.e. LREAL=.FALSE.).'
        case(NoRealSpaceAndYouShould)
            res = 'You have a (more or less) ''large supercell'' and for larger cells &
                &it might be more efficient to use real-space projection operators. &
                &Therefore, try LREAL= Auto in the INCAR file.\n&
                &Mind: For very accurate calculation, you might also keep the &
                &reciprocal projection scheme (i.e. LREAL=.FALSE.).'
        case(IALGO8)
            res = 'Recently, Corning got a patent for the Teter Allan Payne algorithm. &
                &Therefore, VASP.4.5 does not support IALGO=8 any longer. &
                &A much faster algorithm, IALGO=38, is now implemented in VASP. &
                &This algorithm is a blocked Davidson-like method and as reliable as &
                &IALGO=8 used to be.\n&
                &For ultimate performance, IALGO=48 is still the method of choice &
                &-- SO MUCH ABOUT PATENTS :)'
        case(RealSpaceNoMoreRecommended)
            res = 'You have a (more or less) ''small supercell'' and for smaller cells &
                &it is recommended to use the reciprocal-space projection scheme! &
                &The real-space optimization is not efficient for small cells and it &
                &is also less accurate ...\n&
                &Therefore, set LREAL=.FALSE. in the INCAR file.'
        case(WrongOptimizationRealSpace)
            res = 'One real-space projector is optimized for E = ' // str(arg%rval(1)) // &
                ' eV but you are using a cutoff of ENMAX = ' // str(arg%rval(2)) // &
                ' eV (QCUT = ' // str(arg%rval(3)) // ' a.u.). This makes no sense. &
                &Reoptimize the projector with the a.u. value given above.'
        case(LongLattice)
            res = 'One of the lattice vectors is very long (>50 A), but AMIN is rather large. &
                &This can spoil convergence since charge sloshing might occur &
                &along the long lattice vector. If problems with convergence are &
                &observed, try to decrease AMIN to a smaller value (e.g. 0.01).\n&
                &Note: This warning only applies if the self-consistency cycle is used.'
        case(DifferentXCGradTypes, DifferentLDAXCTypes)
            res = 'You have built up your multi-ion-type POTCAR file out of POTCAR &
                &files with incompatible specifications for the XC-types used to &
                &generate the pseudopotential. This makes no sense at all!! What &
                &I found is LEXCH = ' // str(arg%ival(1)) // ' for atom types <= ' // &
                str(arg%ival(3)-1) // ' but LEXCH = ' // str(arg%ival(2)) // ' was &
                &found for atom type = ' // str(arg%ival(3)) // '. Use identical &
                &XC-functionals for the pseudopotential generation for all atom types, please ... !'
        case(NumberOfElectrons)
            res = 'The number of bands is not sufficient to hold all electrons. &
                &I am sorry, but this calculation does not make sense. &
                &If you want to calculate only selected states, you need to set &
                &either EFERMI oder EREF. EREF specifies around which energy states &
                &are supposed to be calculated.\nI found NBANDS = ' // &
                str(arg%ival(2)) // ', NELECT = ' // str(arg%ival(1)) // '.'
        case(NBANDSchanged)
            res = 'The number of bands has been changed from the values supplied &
                &in the INCAR file. This is a result of running the parallel version. &
                &The orbitals not found in the WAVECAR file will be initialized with &
                &random numbers, which is usually adequate. For correlated calculations, &
                &however, you should redo the groundstate calculation.\nI found NBANDS = ' // &
                &str(arg%ival(1)) // '. Now, NBANDS = ' // str(arg%ival(2)) // '.'
        case(DifferentXCFiles)
            res = 'The XC type on the POTCAR file is not compatible with the XC type on the EXHCAR file.'
        case(AM05spin)
            res = 'The AM05 functional for spin-polarized calculations is not supported. &
                &Maybe in a future release...'
        case(NoLDAU)
            res = 'VASP supports LDA+U only for PAW potentials and not for US or NC &
                &pseudopotentials. Please restart with appropriate potentials.'
        case(DoubleCounting)
            res = 'LDA+U in combination with GW calculations causes a double &
                &counting of unidentifiable many-body perturbation terms.'
        case(NoWannier)
            res = 'VASP supports the calculation of maximally-localized Wannier &
                &functions only in the Gamma-only version of the code.'
        case(LREALA)
            res = 'LREAL=A is not well tested yet, please use LREAL=O &
                &or test your results very carefully.'
        case(DifferentSlaterXCTypes)
            res = 'You have built up your multi-ion-type POTCAR file out of POTCAR &
                &files with incompatible specifications for the XC-types used to &
                &generate the pseudopotential. This makes no sense at all!! What &
                &I found is slater parameter = ' // str(arg%rval(1)) // ' for atom &
                &types <= ' // str(arg%ival(1)-1) // ' but slater parameter = ' // &
                str(arg%rval(2)) // ' was found for atom type = ' // str(arg%ival(1)) // &
                '. Use identical slater parameters in the exchange functionals for &
                &the pseudopotential generation for all atom types, please ... !'
        case(DifferentRelXCTypes)
            res = 'You have built up your multi-ion-type POTCAR file out of POTCAR &
                &files with incompatible specifications for the XC-types used to &
                &generate the pseudopotential. This makes no sense at all!! What &
                &I found is that the flag which switches on/off the relativistic &
                &corrections has been set to ' // str(arg%lval(1), compact = .false.) // &
                ' for atom types <= ' // str(arg%ival(1) - 1) // ' but it was &
                &set to ' // str(arg%lval(2), compact = .false.) // ' for atom type no. ' // &
                str(arg%ival(1)) // '. Use identical XC-functionals for ' // &
                'the pseudopotential generation for all atom types, please ...!'
        case(PartialDOS)
            res = 'The partial DOS and the PROCAR file are not evaluated for NPAR/=1. &
                &Please rerun with NPAR=1.'
        case(NoOptics)
            res = 'The file OPTICS can be written only for NPAR=1. &
                &Optical properties can be found in vasprun.xml and OUTCAR, however. &
                &If you need the OPTICS file, rerun with NPAR=1.'
        case(Position)
            res = 'The distance between some ions is very small. &
                &Please check the nearest-neighbor list in the OUTCAR file.\n&
                &I HOPE YOU KNOW WHAT YOU ARE DOING! '
        case(CenterOfMassDrift)
            res = 'The initial velocities result in a center-of-mass drift but &
                &there must be no drift. The drift will be removed! '
        case(EnforcedLDA)
            res = 'You enforced a specific xc type in the INCAR file but &
                &a different type was found in the POTCAR file.\n&
                &I HOPE YOU KNOW WHAT YOU ARE DOING!'
        case(DipolEpsilon)
            res = 'You specified a dielectric constant different from 1 (EPSILON). &
                &This indicates that you perform calculations for the bulk. &
                &At the same time, LDIPOL is set to .TRUE., a mode which does not &
                &work properly for bulk systems!\n&
                &I HOPE YOU KNOW WHAT YOU ARE DOING!'
        case(DipolCubic)
            res = 'LDIPOL = .TRUE. must be selected only for cubic supercells, since &
                &the quadrupole corrections are currently only implemented for this &
                &specific geometry.\n&
                &Please change your box, or set LDIPOL=.FALSE.'
        case(ALGO_A_ISMEAR)
            res = 'ALGO=A and IALGO=5X tend to fail with the tetrahedron method &
                &(e.g. Bloechl''s method, ISMEAR=-5, is not variational). &
                &Please switch to IMSEAR=0-n, except for DOS calculations. &
                &For DOS calculations use IALGO=53 after pre-converging with ISMEAR>=0.\n&
                &I HOPE YOU KNOW WHAT YOU ARE DOING!'
        case(ModelGW)
            res = 'Model GW is not a variational method. &
                &The conjugate-gradient algorithm therefore does not work. &
                &Please switch to IALGO=53-55.\n&
                &I HOPE YOU KNOW WHAT YOU ARE DOING!'
        case(HFUSPP)
            res = 'HF is only implemented for the PAW method and for norm-conserving &
                &pseudopotentials.\n&
                &One of the POTCAR data sets corresponds to an US pseudopotential.'
        case(LHFLinearResponse)
            res = 'HF is not implemented in the linear-response routines. &
                &Please use the finite-difference versions.'
        case(LinearReponsePOSCAR)
            res = 'The POSCAR file is problematic for the linear-response routines &
                &since the displacement is not along Cartesian directions. &
                &Please send this POSCAR to vasp support, we will try to fix the issue. &
                &You can use finite differences for the time being.'
        case(LinearResponseLREAL)
            res = 'Linear response in principle supports LREAL = Auto (or LREAL=.TRUE.). &
                &However, in some cases we have found that the real-space projection &
                &is not sufficiently accurate.\n&
                &I strongly recommend that you switch it off by setting LREAL=.FALSE.'
        case(FockForce)
            res = 'ALLOCATE for the force calculations using Hartree-Fock routines failed. &
                &Your forces are not correct!!!\n&
                &However, the WAVECAR file might be useful for continuation.'
        case(HF_NPAR)
            res = 'HF is only implemented for NPAR = number of nodes. ' // &
                'Otherwise, communication between nodes would be required during the &
                &calculation of the exchange potential, which would be utterly slow.'
        case(NPAREfficiency)
            res = 'For optimal performance we recommend to set \n&
                &  NCORE = 2 up to number-of-cores-per-socket\n' // &
                'NCORE specifies how many cores store one orbital (NPAR=cpu/NCORE). &
                &This setting can greatly improve the performance of VASP for DFT. &
                &The default, NCORE=1 might be grossly inefficient &
                &on modern multi-core architectures or massively parallel machines. &
                &Do your own testing! More info at https://www.vasp.at/wiki/index.php/NCORE\n&
                &Unfortunately you need to use the default for GW and RPA calculations &
                &(for HF NCORE is supported but not extensively tested yet).'
        case(NPAR)
            res = 'VASP internal routines have requested a change of the k-point set. ' // &
                'Unfortunately, this is only possible if NPAR=number of nodes. &
                &Please remove the tag NPAR from the INCAR file and restart the &
                &calculation.'
        case(LHFOneCenter)
            res = 'HF one-center treatment and a full HF treatment are not compatible. &
                &Please either set LHFCALC to .FALSE. or remove the LHFONE tag &
                &from the INCAR file.'
        case(HFParam)
            res = 'A combination of HF parameters was found in the INCAR file that is &
                &not supported, e.g. HFRCUT and HFKIDENT.'
        case(OEPdouble)
            res = 'The EXX-OEP and LHF methods do not support a different grid for HF &
                &and the local Hamiltonian.\nPRECFOCK was set to Normal.'
        case(NMAXFOCKAE)
            res = 'You use NMAXFOCKAE=2 in the INCAR file. Be careful:\n&
                &Adding two functions to mimic the all-electron charge density on the &
                &plane-wave grid requires very fine FFT grids. Please set NGX manually &
                &or reconsider your choice.'
        case(ISYM2)
            res = 'You have selected ISYM=2 for a HF type calculation. This will &
                &symmetrize the charge density but not the exchange potential. &
                &I suggest to use ISYM=3 instead. This uses symmetry to obtain ' // &
                'the orbitals at all k-points in the Brillouin zone, but does not &
                &apply symmetry to the Hartree potential directly. &
                &The resultant charge might have lower symmetry than the crystal &
                &but, at least, Hartree and exchange are fully compatible.'
        case(ISYM_MD)
            res = 'ISYM>0 for a molecular dynamics simmulation is not recommended.'
        case(LSPECTRALinefficient)
            res = 'LSPECTRAL is very inefficient if not sufficient memory is available. &
                &The available memory per core should be set by MAXMEM = XXXX. &
                &The default for MAXMEM is presently XXXX=2800 (MByte). ' // &
                'If this has already been done and if the warning prevails, &
                &you should either decrease NOMEGA or decrease ENCUTGW, or distribute ' // &
                'the calculations onto more nodes.\n&
                &Note: The present calculations might lead to excessive memory paging!!'
        case(LSPECTRALGWset)
            res = 'The default for LSPECTRALGW is now LSPECTRALGW=LSPECTRAL. &
                &This will change QP energies but increases the numerical precision &
                &of the QP slightly (for instance, Z factors are usually more accurate). &
                &You can use the old version by specifying LSPECTRALGW=.FALSE.'
        case(LPOTOK)
            res = 'Presently the RPA-OEP routine requires reading in a POT file &
                &containing an initial guess for the local potential. This could be &
                &generated by using the KLI or EXX-OEP routines (EXXOEP=1 or 2) and &
                &setting LVTOT=.TRUE. in order to write the POT file. &
                &Supply this file by setting ICHARG=4.'
        case(XiSingular)
            res = 'The response function can not be inverted. &
                &This indicates that the number of linearly independent transitions &
                &is smaller than the rank of the response function. &
                &Try to increase the number of k-points, or try to decrease ENCUTGW.'
        case(AEXXzero)
            res = 'You calculate the electron-hole interaction using AEXX=0. &
                &This is hardly what you want to do.\n&
                &Maybe you have forgotten to set AEXX in the INCAR file &
                &(For LHFCALC=.TRUE., the default is AEXX=0.25 but if LHFCALC is not &
                &set, the default is AEXX=0.0).'
        case(GWoptics)
            res = 'The derivative of the wavefunctions with respect to k (WAVEDER) &
                &can not be found. You should redo the groundstate calculation &
                &using LOPTICS=.TRUE. in order to write the WAVEDER file. &
                &However for metals, the present setting is ok.'
        case(GW_HFoptics)
            res = 'Presently LREAL is not supported for the calculation of optical &
                &properties and HF type Hamiltonians.\nPlease use LREAL=.FALSE.'
        case(ENCUTGW_LHF)
            res = 'ENCUTGW is set for a local HF calculation. &
                &This might cause very slow convergence with respect to &
                &ENCUTGW, since even the Hartree potential is cut off. &
                &Use with extreme care (or set ENCUTGW=0 in the INCAR file).'
        case(DeprENCUTFOCK)
            res = 'ENCUTFOCK is no longer supported!\n&
                &Use PRECFOCK to select the mode for HF type calculations:\n&
                &PRECFOCK= L low (coarse grid for HF, normal augmentation charge)\n&
                &PRECFOCK= M medium (normal grid for HF, normal augmentation charge)\n&
                &PRECFOCK= F fast (coarse grid for HF, soft augmentation charge)\n&
                &PRECFOCK= N normal (PREC=N grid for HF, soft augmentation charge)\n&
                &PRECFOCK= A accurate (PREC=A grid for HF, soft augmentation charge)\n\n&
                &L is equivalent to ENCUTFOCK=0 in vasp.5.2.2.\n&
                &M is equivalent to vasp.5.2.2 if ENCUTFOCK is not set.\n&
                &M&L cause significant noise in the forces and are no longer recommended.\n&
                &N is recommended for routine calculations, as there is little noise to be expected.\n&
                &A is now recommended for very accurate calculations.\n&
                &F is now recommended for quick calculations (even phonons often ok, &
                & speedup between 2 and 4).'
        case(PRECFOCK)
            res = 'Use PRECFOCK to select the mode for HF type calculations:\n&
                &PRECFOCK= L low (coarse grid for HF, normal augmentation charge)\n&
                &PRECFOCK= M medium (normal grid for HF, normal augmentation charge)\n&
                &PRECFOCK= F fast (coarse grid for HF, soft augmentation charge)\n&
                &PRECFOCK= N normal (PREC=N grid for HF, soft augmentation charge)\n&
                &PRECFOCK= A accurate (PREC=A grid for HF, soft augmentation charge)\n\n&
                &L is equivalent to ENCUTFOCK=0 in vasp.5.2.2.\n&
                &M is equivalent to vasp.5.2.2 if ENCUTFOCK is not set.\n&
                &M&L cause significant noise in the forces and are no longer recommended.\n&
                &N is recommended for routine calculations, as there is little noise to be expected.\n&
                &A is now recommended for very accurate calculations.\n&
                &F is now recommended for quick calculations (even phonons often ok, &
                & speedup between 2 and 4).'
        case(KPOINTS_HF)
            res = 'Your generating k-point grid is not commensurate with the symmetry &
                &of the lattice. This can cause slow convergence with respect &
                &to k-points for HF-type calculations.\n&
                &Suggested SOLUTIONS:\n&
                & ) If not already the case, use automatic k-point generation,\n&
                & ) and/or shift your grid to Gamma (G) (e.g. required for\n&
                &   hexagonal, fcc, and fcc-orthorhombic lattices).\n&
                & ) and/or set KBLOWUP = .FALSE.'
        case(KPOINTS_PEAD)
            res = 'Your generating k-point grid is not commensurate with the symmetry &
                &of the lattice. This does not sit well in combination with the &
                &PEAD routines, sorry ...\n&
                &Suggested SOLUTIONS:\n&
                & ) If not already the case, use automatic k-point generation,\n&
                & ) and/or shift your grid to Gamma (G) (e.g. required for\n&
                &   hexagonal, fcc, and fcc-orthorhombic lattices).\n&
                & ) and/or set KBLOWUP = .FALSE.'
        case(KPOINTS_INTER)
            res = 'Your generating k-point grid is not commensurate with the symmetry &
                &of the lattice. Presently, it is not possible to write the velocity ' // &
                'operator in the extended zone scheme for such k-point grids.\n&
                &Suggested SOLUTIONS:\n&
                & ) If not already the case, use automatic k-point generation,\n&
                & ) and/or shift your grid to Gamma (G) (e.g. required for\n&
                &   hexagonal, fcc, and fcc-orthorhombic lattices).\n&
                & ) and/or set KBLOWUP = .FALSE.'
        case(GWkweight)
            res = 'The weights do not sum up to one for the summation over q in the &
                &GW routine. To save your life, they have been renormalized to 1. &
                &This circumstance can occur under certain conditions. &
                &Please check your results carefully.'
        case(GWnoreal)
            res = 'Presently, the real-space GW code only supports Gamma-point-only &
                &calculations. Please use the standard version for arbitrary k-point &
                &grids.'
        case(ACFDTRnotSupported)
            res = 'Presently, the real-space GW code requires scaLAPACK. &
                &This code is designed for massively parallel computations. &
                &Please use the standard code for routine calculations.'
        case(PEADnoVirtuals)
            res = 'Since you included no unoccupied states, the PEAD routines can not &
                &estimate whether EFIELD_PEAD exceeds the critical field strength. &
                &Please check your results carefully.'
        case(PEADlargeField)
            res = 'One or more components of EFIELD_PEAD are too large for comfort. In all &
                &probability, you are too near to the onset of Zener tunneling.\n' // &
                '\n     e |E dot A_1| = ' // str(arg%rval(1), '(F8.5)') // ' > 1/10, E_g/N_1 = ' // str(arg%rval(4), '(F8.5)') // &
                '\n     e |E dot A_2| = ' // str(arg%rval(2), '(F8.5)') // ' > 1/10, E_g/N_2 = ' // str(arg%rval(5), '(F8.5)') // &
                '\n     e |E dot A_3| = ' // str(arg%rval(3), '(F8.5)') // ' > 1/10, E_g/N_3 = ' // str(arg%rval(6), '(F8.5)') // &
                '\n\nPossible SOLUTIONS:\n&
                & ) Choose a smaller electric field.\n&
                & ) Use a less dense grid of k-points.'
        case(PEAD_NCOREnot1)
            res = 'The PEAD routines do not work for NCORE /= 1 (or NPAR /= NCPU). &
                &Several features internally rely on the PEAD routines, so even if you &
                &have not specified LPEAD=.TRUE., try restarting your job with NCORE=1.'
        case(PEADPolNoMetal)
            res = 'The calculation of the macroscopic polarization by means of the &
                &Berry-phase expressions (LCALCPOL=.TRUE.) requires your system to be &
                &insulating. This does not seem to be the case.\n&
                &VASP will skip this part of your job, sorry...'
        case(PEADEpsNoMetal)
            res = 'The calculation of changes in the macroscopic polarization due to the &
                &application of a finite electric field, by means of the PEAD method &
                &(LCALCEPS=.TRUE. or EFIELD_PEAD/=0), requires your system to be &
                &insulating. This does not seem to be the case.\n&
                &VASP will skip this part of your job, sorry...'
        case(PEADnowIsMetal)
            res = 'The application of a finite electric field, by means of the PEAD &
                &method (LCALCEPS=.TRUE. or EFIELD_PEAD/=0), requires your system to be &
                &insulating. This does not seem to be the case (any longer). You might try &
                &to reduce the size of the electric field and/or the smearing width.'
        case(LINTET)
            res = 'Changing the smearing to ISMEAR=1 (the default) for those parts &
                &of this job that use k-points specified in line-mode.'
        case(GRDSHFT_NOTET)
            res = 'The linear tetrahedron method can not be used together with the &
                &grid-shift interpolation routine.'
        case(GammaK)
            res = 'You are using the Gamma-point-only version with more than one k-point &
                &or some other non-Gamma k-point.'
        case(KPOINTS_IO)
            res = 'Error reading KPOINTS file.\n&
                &The error occurred at line: ' // str(arg%ival(1)) // '.'
        case(NoInversion)
            res = 'Full k-point grid generated. Inversion symmetry is not applied.'
        case(FFTGridIsNotSufficient)
            res = 'Your FFT grids (NGX,NGY,NGZ) are not sufficient for an accurate &
                &calculation. Thus, the results might be wrong. Good settings for NGX, &
                &NGY and NGZ are ' // str(arg%ival(1)) // ', ' // str(arg%ival(2)) // &
                ' and ' // str(arg%ival(3)) //  ', respectively.\n&
                &Mind: This setting results in a small but reasonable wrap-around error. &
                &It is also necessary to adjust these values to the FFT routines you use.'
        case(MetaGGAandForces)
            res = 'You have switched METAGGA and/or ASPHERICAL radial calculations on &
                &but ions are moved based on forces. &
                &METAGGA and LASPH change only the exchange and correlation energy. &
                &You might want to reconsider your choices.\n&
                &Current values: LASPH=' // str(arg%lval(1), compact = .false.) // &
                ' and LMETAGGA=' // str(arg%lval(2), compact = .false.) // '.'
        case(MetaGGAResponse)
            res = 'You have switched on METAGGA and linear-response routines. &
                &This combination is currently not supported. &
                &If you have selected LOPTICS=.TRUE., VASP will continue but some &
                &matrix elements [H,r] are neglected (you might try to use &
                &LPEAD=.TRUE. which works irrespective of H). &
                &In all other cases, VASP will stop. SOOO sorry...'
        case(SmallCSHIFT)
            res = 'You are using a small complex shift. This can yield unreliable results &
               &when the Kramers-Kronig relation is used for systems with spin-orbit &
               &interactions.\nPlease increase CSHIFT in the INCAR file.'
        case(TetNoMetal)
            res = 'Integration of the imaginary part of the dielectric function by means &
               &of the tetrahedron method is not supported for metallic systems. &
               &VASP will skip the calculation of the optical properties.\n&
               &Please rerun with ISMEAR > -4.'
        case(POTIMlarge)
            res = 'Your timestep is larger than 0.1 Angst.\n&
                &For finite differences, this really does not make sense. I will &
                &reset POTIM to 0.015. I recommend to use 0.01 to 0.02 for finite &
                &differences.'
        case(GRDSHFTnotInIBZ)
            res = 'The grid-shift interpolation procedure is not able to find this point in the IBZ.'
        case(GRDSHFTK)
            res = 'The grid-shift interpolation procedure only works for odd KINTER.'
        case(NoOcc)
            res = 'The energy of the highest-lying unoccupied state is less than the &
                &lowest lying occupied states. Are you sure your system contains &
                &unoccupied states?'
        case(NoOddGrid)
            res = 'The grid-shift interpolation procedure does not work for even grid &
                &divisions.'
        case(GammaC)
            res = 'The k-points need to be Gamma centered.'
        case(BSEantires)
            res = 'BSE calculations beyond the Tamm-Dankoff approximation require you &
                &to set the flag LORBITALREAL=.TRUE. in the groundstate as well as &
                &in the BSE calculation.\n&
                &Please redo the groundstate calculation with this flag set.'
        case(BSEconjg)
            res = 'BSE calculations beyond the Tamm-Dankoff approximation require you &
                &to set the flag LORBITALREAL=.TRUE. in the groundstate as well as &
                &in the BSE calculation.\n&
                &Please redo the groundstate calculation with this flag set.'
        case(BSEncv)
            res = 'The number of conduction band and valence band pairs is zero. &
                &This indicates that NBANDSV or NBANDSO are not properly set up. &
                &The VASP default is to set them to the number of occupied orbitals. &
                &Try to set them manually.\n&
                &It is also possible that OMEGAMAX is too small; please increase it.'
        case(Wmissing)
            res = 'The file WXXXX.tmp for XXXX=' // str(arg%ival(1)) // ' can not be read. &
                &Rerun the GW calculations with LRPA= .TRUE. and do not change &
                &ENCUTGW or NBANDS between runs.'
        case(AuxWAVECAR)
            res = 'An auxiliary WAVECAR file was read.\n&
                &W will be calculated using these auxiliary wavefunctions. &
                &The WAVEDER.XXX file must match the WAVECAR file!'
        case(SpinSpiral)
            res = 'To represent the spin spiral you requested, with a kinetic &
                &energy cutoff of ENINI= ' // str(arg%rval(1)) // ' eV, &
                &choose ENMAX > ' // str(arg%rval(3)) // ' eV. &
                &Currently, ENMAX= ' // str(arg%rval(2)) // ' eV.'
        case(HighestBandsOcccupied)
            res = 'Your highest band is occupied at some k-points! Unless you are &
                &performing a calculation for an insulator or semiconductor, without &
                &unoccupied bands, you have included TOO FEW BANDS!! Please increase &
                &the parameter NBANDS in file INCAR to ensure that the highest band &
                &is unoccupied at all k-points. It is always recommended to &
                &include a few unoccupied bands to accelerate the convergence of &
                &molecular-dynamics runs (even for insulators or semiconductors), &
                &since the presence of unoccupied bands improves wavefunction &
                &prediction and helps to suppress ''band-crossings''.' // &
                fermiWeights(arg%ival, arg%rval)
        case(LVTOT)
            res = 'For LVTOT=.TRUE., VASP.5.x writes the TOTAL local potential to &
                &the file LOCPOT. If you want only the Hartree contributions, use &
                &LVHAR=.TRUE. instead. '
        case(vdWKlimes)
            res = 'You have switched on vdW-DFT.\n&
                &This routine was written and supplied by Jiri Klimes. &
                &We recommend that you carefully read and cite the following &
                &publications:\n&
                &  J. Klimes, D.R. Bowler, A. Michelides\n&
                &     J. Phys.: Cond Matt. 22 022201 (2010)\n&
                &  J. Klimes, D.R. Bowler, A. Michelides\n&
                &     Phys. Rev. B. 83, 195131 (2011)\n&
                &and references therein.'
        case(NoKPAR)
            res = 'You have enabled k-point parallelism (KPAR>1) but this does not &
                &work for the interpolation routines.'
        case(FERWEOptics)
            res = 'Integration of the imaginary part of the dielectric function by &
                &means of fixed occupancies (ISMEAR=-2) is not supported. VASP will &
                &skip the calculation of the optical properties.'
        case(OmegaPar)
            res = 'NOMEGAPAR must be smaller than NOMEGA and the used number of MPI &
                & ranks. I am setting NOMEGAPAR to '//str(arg%ival(1))//', the largest possible value.'
        case(TauPar)
            res = 'NTAUPAR must be smaller than NOMEGA and the used number of MPI &
                & ranks. I am setting NTAUPAR to '//str(arg%ival(1))//', the largest possible value.'
        case(TauParCPU)
            res = 'The number of MPI ranks ('//str(arg%ival(1))//') is not&
                  & divisible by NTAUPAR ('//str(arg%ival(2))//') and/or NOMEGAPAR &
                  &('//str(arg%ival(3))//').'
        case(ImagGridParForLTMP2)
            res = 'The present laplace transformed MP2 implentation does not&
                & support parallelization over imaginary time and frequency &
                & points. NTAUPAR and NOMEGAPAR are set to 1 in the following.'
        case(chiGG_NOmega)
            res = 'Presently, the real-time/real-space code requires a minimum of &
                &3 imaginary times.\nPlease set NOMEGA to a sensible number. &
                &Typically 8 imaginary times suffice.'
        case(OmegaGridOMEGAMIN)
            res = 'You have a very small band gap. This results in slow convergence &
                &of the energy with respect to NOMEGA. &
                &Usually, VASP defaults OMEGAMIN to the bandgap which guarantees an &
                &optimal grid for OMEGAGRID=140. &
                &The gap was found to be OMEGAMIN=' // str(arg%rval(1)) // '. &
                &Now, OMEGAMIN is set to 0.02 instead, which effectively removes some &
                &excitations from the spectrum.\n&
                &You can set OMEGAMIN manually in the INCAR, however.'
        case(WannNotFound)
            res = 'The file wannier90.win can not be found.'
        case(TargetNotSet)
            res = 'You have not defined target states for the projection onto the &
                &subspace. VASP will continue and use all states of the localized &
                &basis as targets.\n&
                &You can set target states in the INCAR file with NTARGET_STATES'
        case(NTARGET_STATES)
            res = 'NTARGET_STATES contains more non-zero values than target states &
                &available. VASP will continue and neglect screening effects of &
                &ALL Wannier states! If this is not desired NTARGET_STATES has to &
                &be set properly in the INCAR.'
        case(NCRPALOW)
            res = 'NCRPALOW is not set in INCAR.\n&
                &NCRPALOW should correspond to the lowest band represented by &
                &the localized basis.'
        case(CRPADisentangle)
            res = 'CRPA-disentanglement method needs a minimal basis with the same size &
                &as the target space. Your chosen targets span a ' // str(arg%ival(2)) // &
                '-dimensional sub space in a space with ' // str(arg%ival(1)) // &
                ' dimensions. Change either NTARGET_STATES or the localized basis set. &
                &Alternatively, use the standard CRPA or the weighted method.'
        case(LCAO_INCAR)
            res = 'You have not specified an ion number WANPROJ_I and/or an angular &
                &momentum WANPROJ_L in the INCAR file. This is mandatory if you &
                &want to use LCAOs!!'
        case(SHMEM_CRPA)
            res = 'NCSHMEM>1 is not supported for CRPA calculations.\n&
                &Please set NCSHMEM=1 for all CRPA calculations.'
        case(OmegaGridOmegaTL)
            res = 'The maximum frequency in the frequency grid is smaller than the &
                &largest excitation energy. This can cause errors in the &
                &correlation energy. Please increase NOMEGA to make sure that the &
                &error is sufficiently small.'
        case(NOmega32)
            res = 'OMEGAGRIRD = 140 does not allow for more than 32 grid points. &
                &For this method, the error decays exponentially with the number of &
                &grid points. Typically, NOMEGA=8-16 is sufficient even for small-gap &
                &systems. Please do your own testing.'
        case(OmegaGrid140)
            res = 'Presently, the real-time/real-space code requires either &
                &OMEGAGRID=140 (minimax grids) or OMEGAGRID=150 (least square grids).'
        case(NUPDOWN)
            res = 'NUPDOWN is presently not supported because it can result in a &
                &different EFERMI for up and down electrons. Switching it off now &
                &and continuing... Please check that occupancies are still ok by, e.g., ' // &
                'performing a one-step KS calculation.'
        case(NPAR_IALGO8)
            res = 'IALGO = 8 supports only parallelism over k-points (KPAR) and &
                &distribution of orbitals over cores.\n&
                &Specifically, NPAR must be set to 1 and NCORE should not be set in &
                &the INCAR file.'
        case(SHMEMError)
            res = 'Can not allocate shared memory. &
                &The simple solution is to run using NCSHMEM=1. &
                &However, this can increase the memory requirements significantly. &
                &If possible, increase the available shared memory. &
                &The available memory can be set in the kernel. &
                &/sbin/sysctl kernel{max,all,mni} can be used to determine the &
                &available memory. ''max'' often needs to be increased (half the physical &
                &memory is often a good choice).'
        case(KProjection)
            res = 'The projection LKPROJ is presently only implemented for NCORE=1. &
                &Please set NCORE=1 (or NPAR to the total number of cores) &
                &and restart the calculations.'
        case(Gamma)
            res = 'The file GAMMA is not compatible with the present calculation. &
                &The first line must specify the number of ions, bands and k-points.'
        case(Gamma2)
            res = 'The file GAMMA can not be read. Each set must start with the k-point &
                &number the start and final band for which the correction of the &
                &density matrix is supplied. NK=' // str(arg%ival(1)) // '.'
        case(Alternant)
            res = 'Could not find all extrema of error function. Try to increase the &
                &number of sampling points. Number of points used in this run: &
                &NMAXALT=' // str(arg%ival(1)) // '.'
        case(LBONEandLZORA)
            res = 'You have set both LBONE=.TRUE. as well as LZORA=.TRUE.. The former will &
                &take precedence, i.e., LZORA will now be set to .FALSE.'
        case(RandomSeed)
            res = 'The Random Number Generator may be seeded as follows:\n&
                &  0 <= SEED1 < 900000000\n&
                &  0 <= SEED2 < 1000000\n&
                &  0 <= SEED3\n&
                &Your current RANDOM_SEED:' // str(arg%ival(1)) // ' ' // &
                str(arg%ival(2)) // ' ' // str(arg%ival(3)) // '.'
        case(NBANDlessNWAN)
            res = 'The total number of bands (' // str(arg%ival(1)) // ') &
                &is smaller than the available number of local projection centers (' // &
                str(arg%ival(2)) // '). For generating Wannier functions, this &
                &does not make sense.\n Please increase the total number of &
                &bands by setting a higher value for NBANDS in the INCAR file. &
                &Alternatively, you might want to lower the number of local &
                &functions. Consult the documentation of the LOCPROJ flag for &
                &further information.'
        case(ELPHtesting)
            res = 'Electron-Phonon testing mode is active (LELPH_TESTING = .TRUE.). &
                &Unwanted behavior might occur. I highly recommend disabling it.'
        case(GDESMatBlockRed)
            res = 'One or more MPI groups includes at least one CPU that does not &
                &carry data. This would cause ScaLAPACK to crash and is due to &
                &the automatically selected processor grid and blocking factor NB. &
                &To overcome this problem NB has been divided by 2.\n&
                &However, this will slow down your calculation and it is suggested &
                &to change NTAUPAR or NBANDS or both in the INCAR in order to change &
                &the processor grid (NPROW, NPCOL) such that (NCOL/NTAUPAR-1)*MB &
                &< NBANDS as well as (NROW/NTAUPAR-1)*MB < NBANDS holds true.&
                &Values for this run are:&
                &\nNTAUPAR = ' // str(arg%ival(5)) // &
                '\nNBANDS  = ' // str(arg%ival(4)) // &
                '\nNPCOL   = ' // str(arg%ival(1)) // &
                '\nNPROW   = ' // str(arg%ival(2)) // &
                '\nMB      = ' // str(arg%ival(3))
        case(GDESMatSize)
            res = 'One or more MPI groups includes at least one CPU that does not &
                &carry data. This would cause ScaLAPACK to crash and is due to &
                &the automatically selected processor grid for ScaLAPACK.\n&
                &You can influence the processor grid (NPROW, NPCOL) by changing &
                &NTAUPAR or NBANDS or both such that (NPCOL/NTAUPAR-1)*MB < NBANDS &
                &as well as (NPROW/NTAUPAR-1)*MB < NBANDS. Values for this run are:&
                &\nNTAUPAR = ' // str(arg%ival(5)) // &
                '\nNBANDS  = ' // str(arg%ival(4)) // &
                '\nNPCOL   = ' // str(arg%ival(1)) // &
                '\nNPROW   = ' // str(arg%ival(2)) // &
                '\nMB      = ' // str(arg%ival(3))
        case(NOMEGA_DUMPtooLarge)
            res = 'I found a value for NOMEGA_DUMP > NOMEGA, will use NOMEGA_DUMP &
                &= NOMEGA instead.'
        case(SigmaTMBPT)
            res = 'The temperature for the Matsubara axis differs from the one found &
                &in the WAVECAR. This might cause severe integration errors for &
                &GW and RPA calculations. TMBPT should not be used in production &
                &runs.\n&
                &The inverse temperature (beta) is set via the k-point smearing &
                &beta=1/SIGMA.'
        case(FermiFiniteTemp)
            res = 'Finite temperature many-body pertubation algorithms work only &
                &with Fermi smearing weights (ISMEAR=-1).'
        case(BlochBasisUsed)
            res = 'You have not specified any basis functions for your CRPA run. &
                &VASP will use Bloch functions. You can define a basis set with the ' // &
                'LOCPROJ tag or using the VASP2WANNIER90 interface.'
        case(TargetOut)
            res = 'You selected more target states than present for this run.\n&
                &Check your NTARGET_STATES setting.'
        case(NCRPA_BANDSout)
            res = 'You selected more CRPA bands than present for this run.\n&
                &Check your NCRPA_BANDS setting.'
        case(LOCPROJnotSet)
            res = 'You have not specified any basis functions with the LOCRPOJ tag &
                &for your CRPA calculation!'
        case(ALGO_GW)
            res = 'ALGO = GW will be deprecated in future releases. Please use &
                &ALGO = EVGW from now on.'
        case(ALGO_GW0)
            res = 'ALGO = GW0 will be deprecated in future releases. Please use &
                &ALGO = EVGW0 from now on.'
        case(ALGO_G0W0)
            res = 'ALGO = G0W0 will be deprecated in future releases. Please use &
                &ALGO = EVGW0 ; NELM=1 from now on.'
        case(ALGO_scGW0)
            res = 'ALGO = scGW0 will be deprecated in future releases. Please use &
                &ALGO = QPGW0 from now on.'
        case(ALGO_scGW)
            res = 'ALGO = scGW will be deprecated in future releases. Please use &
                &ALGO = QPGW from now on.'
        case(POTCARmissing)
            res = 'Cannot find POTCAR file!'
        case(SymprecDependence)
            res = 'The determination of the symmetry of your systems shows a strong &
                &dependence on the tolerance parameter SYMPREC. This is often the result &
                &of "noise" in the structural parameters (lattice vectors and/or &
                &ionic positions). Set NWRITE=3 for more details.'
        case(ALLOC_RESPONSE)
            res = 'Could not allocate body of response function on mpi rank 0 of size: '&
                &//str(arg%ival(1))//' MB. Reducing NOMEGAPAR or using more computing ' // &
                'nodes might solve this problem.'
        case(ALLOC_RESPONSE_DISTR)
            res = 'Could not allocate body of response function on mpi rank 0 of size: '&
                &//str(arg%ival(1))//' MB. Reducing NTAUPAR or using more computing ' // &
                'nodes might solve this problem.'
        case(ALLOC_WAVEFUN)
            res = 'Could not allocate wavefunction on mpi rank 0 of size: '&
                &//str(arg%ival(1))//' MB. Reducing ENCUT and/or using more ' // &
                'computing nodes might solve this problem.'
        case(MAXMEM_LOW)
            res = 'This job will probably crash, due to insufficient memory &
                &available. Available memory per mpi rank: '//str(arg%ival(1))//' MB, &
                &required memory: '//str(arg%ival(2))//' MB. Reducing ' // &
                'NTAUPAR or using more computing nodes might solve this problem.'
        case(AUTOSET_MEM)
            res = 'Failed to automatically set available memory for this job. &
                &Please set MAXMEM by hand, because the default value of '//&
                &str(arg%ival(1))//' MB might be grossly inefficient.'
        case(RPAFORCES_COMPAT)
            res = 'Currently RPA forces (LRPAFORCES = .TRUE.) can be calculated &
                &only with the formalism of euclidean quantum field theory &
                &(imaginary time Greens functions). That is, only ALGO = &
                &RPAR[K] | G0W0R[K] or similar options are supported.'
        case(RPAENCUTGWSOFT)
            res = 'The default for ENCUTGWSOFT is different in this version of VASP. &
                & If you wish to obtain identical results as using pre vasp.6.3 &
                & please manually set ENCUTGWSOFT =' // str(arg%rval(1)) // ' in the INCAR file.'
        case(VMBJ_BrentInitialPositions)
            res = 'VMBJ_BRENT: Unable to find initial points.'
        case(VMBJ_TooManyIterations)
            res = 'VMBJ_BRENT: Too many iterations.'
        case(ORTHONORMALIZATION_FAILED)
            select case (arg%ival(1))
            case (1)
                work = 'LAPACK: Routine ZPOTRF failed!'
            case (2)
                work = 'LAPACK: Routine ZTRTRI failed!'
            case (3)
                work = 'scaLAPACK: Routine ZPOTRF ZTRTRI failed!'
            case default
                work = 'unknown argument'
            end select
            res = 'Orbital orthonormalization failed in the inversion of matrix\n' &
                // work // ' kpoint: ' // str(arg%ival(2)) // ' spin: ' // str(arg%ival(3)) // '\n\n&
                &Possible solutions for this issue are:\n&
                &- Please check whether the atoms are too close to each other.\n&
                &- Decreasing POTIM might help during relaxations or MDs.\n&
                &- Deleting CHGCAR or WAVECAR if incompatible with the structure.\n&
                &- Try calculating without symmetrizing (ISYM = 0 or ISYM = -1).\n&
                &- As last resort: using a different PAW for some of the atoms.'
        case default
            res = 'undefined tutor message ' // str(topic)
        end select
    end function generateMessage

    function fermiWeights(kptIndex, weight) result(res)
        use string, only: str
        integer, allocatable, intent(in) :: kptIndex(:)
        real(q), allocatable, intent(in) :: weight(:)
        character(len=:), allocatable :: res
        integer iband
        logical writeWeights
        writeWeights = .false.
        if (allocated(kptIndex) .and. allocated(weight)) then
            writeWeights = size(kptIndex) == size(weight)
        end if
        if (writeWeights) then
            res = '\n\nIn the following, all k-points are listed with the &
                &Fermi weights of the highest band given in paranthesis:\n'
            do iband = 1, size(kptIndex) - 1
                res = res // str(kptIndex(iband), '(3X,I5)') // ' (' // &
                    str(weight(iband), '(F8.5)') // ')\n'
            end do
            iband = size(kptIndex)
            res = res // 'The total occupancy of band no. ' // str(kptIndex(iband)) // &
                ' is ' // str(weight(iband)) // ' electrons...'
        else
            res = ''
        end if
    end function fermiWeights

    subroutine deactivateVtutorOutputForTest()
        vtutor%discardOutputToNegativeUnit = .true.
        vtutor%unitOut = -1
        vtutor%unitErr = -1
        vtutor%unitBug = -1
    end subroutine deactivateVtutorOutputForTest

end module tutor
