!#define debug
#include "symbol.inc"
!************************************************************************
!
!  this module contains part of the routines required to perform LCAO
!  written by Roman Wahl
!  04/2020 Re-writting of the root finding algorithm by Henrique Miranda
!***********************************************************************
  
  MODULE LCAO

    USE prec
    USE nonlr
    USE mgrid
    USE paw

    IMPLICIT NONE

    TYPE atoms
       TYPE (rgrid) :: FR                ! radial grid
       INTEGER :: OCCUPIED_VALENCE_STATES
       INTEGER :: ATOMLMMAX              ! number of LM-entries for occupied wavefunctions
       INTEGER :: LEXCH                  ! exchange type
       INTEGER :: LMAX                   ! number of entries in LPS for occupied channels
       INTEGER :: LMMAX                  ! number of LM-entries for LPS
       INTEGER :: MAXL                   ! maximum l quantum number appearing in the atom (depends in QPAW)
       INTEGER :: NPLIN                  ! number of points on lin grid (at this point equal to NPRSNL
       INTEGER, ALLOCATABLE :: NBF(:)        ! l-dep number of bessel basis functions
       INTEGER, ALLOCATABLE :: PSRMAX(:)
       INTEGER, ALLOCATABLE :: IRCHECK(:)    ! corresponding index on the r-grid for RCHECK
       INTEGER, ALLOCATABLE :: N(:)          ! main quantum number (since this a pseudopotential description it starts with 1)
                                         ! -1 no used
       INTEGER, ALLOCATABLE :: POTCAR_NODES(:)      ! # of allowed nodes (Kalium!!!)
       INTEGER, ALLOCATABLE :: LCAO_NODES(:)      ! # of allowed nodes (Kalium!!!)
       INTEGER, ALLOCATABLE :: LPS(:)    ! copy of LPS of PP
       INTEGER, ALLOCATABLE :: L(:,:)    ! L(occupied_valence_states,3): column 1 --> l quantum number
                                         !                                      2 --> start in array PP%LPS
                                         !                                      3 --> end in array PP%LPS
                                         !                                      4 --> position of wavefunction in WPS 
                                         !                                            important if double occupied
                                         ! necessary for solving wavefunctions
       REAL(q) :: RCOVI                  ! covalent radius, used for the confinement potential
       REAL(q) :: RCOVI_SHIFT            ! shift used for confinment
       REAL(q) :: VPSRMAX
       REAL(q) :: DEXCCORE
       REAL(q) :: ZVALF
       REAL(q), ALLOCATABLE :: RCHECK(:)     ! radius where inward and outward calculated wavefunctions are matched
       REAL(q), ALLOCATABLE :: WFRMAXL(:)    ! real space cutoff for the wave functions
       REAL(q) :: WFRMAX                 ! real space cutoff for the wave functions
       REAL(q) :: SUPGRIDMAX             ! cut off of support grid (equal to highest WFRMAX)
       REAL(q), ALLOCATABLE :: OCC(:)        ! occupancy for this channel
       REAL(q), ALLOCATABLE :: POTCAR_E(:)   ! initial energy for this channel (from POTCAR)
       REAL(q), ALLOCATABLE :: E(:)          ! energy for this channel
       REAL(q), ALLOCATABLE :: POT(:)        ! self consistency potential
       REAL(q), ALLOCATABLE :: POTC(:)       ! extended core potential
       REAL(q), ALLOCATABLE :: POTAE(:)      ! ae potential
       REAL(q), ALLOCATABLE :: POTAEC(:)     ! ae core-potential
       REAL(q), ALLOCATABLE :: POTPS(:)      ! local PP on r-grid for atomic reference configuration (valence only)
       REAL(q), ALLOCATABLE :: POTPSC(:)     ! local PP on r-grid (core only), V_H[\tilde{n}_Zc]
       REAL(q), ALLOCATABLE :: AUG(:,:)      ! L-dependent augmentation charge on r-grid
       REAL(q), ALLOCATABLE :: WFCT_L(:,:,:) ! calculated radial pseudo-wavefunctions on log grid
       REAL(q), ALLOCATABLE :: WFCT_G(:,:,:) ! Interpolating function of the Fourier transform of WFCT_L
       REAL(q), ALLOCATABLE :: WFCT(:,:,:)   ! calculated radial pseudo-wavefunctions
       REAL(q), ALLOCATABLE :: WKIN(:,:,:)   ! second derivative of wavefunctions for kinetic energy
       REAL(q), ALLOCATABLE :: WPS(:,:)      ! partial radial pseudo-wavefunctions (gridsize of POTCAR!)
       REAL(q), ALLOCATABLE :: WAE(:,:)      ! partial radial ae-wavefunctions (gridsize of POTCAR!)
       REAL(q), ALLOCATABLE :: RHO(:)        ! self consistent density
       REAL(q), ALLOCATABLE :: RHOPS(:)      ! frozen core pseudo partial density (1:FNMAX)
       REAL(q), ALLOCATABLE :: RHOAE(:)      ! frozen core all electron partial density (1:PP%R%NMAX)
       REAL(q), ALLOCATABLE :: PROJLOG(:,:)  ! projection operator on full logarithmic grid
       REAL(q), ALLOCATABLE :: DIJ(:,:)      ! strength
       REAL(q), ALLOCATABLE :: QIJ(:,:)      !
       REAL(q), ALLOCATABLE :: QION(:,:)     !
       REAL(q), ALLOCATABLE :: DION(:,:)     !
       REAL(q), ALLOCATABLE :: B(:,:)        ! coefficients for besselization
       REAL(q), ALLOCATABLE :: A(:,:,:)      ! coefficients for the orthonormal bessel basis
       REAL(q), ALLOCATABLE :: Q(:,:)        ! q's for bessel functions
       CHARACTER(LEN=2) ELEMENT              ! Name of element
    END TYPE atoms

    !> Container for some arguments to be passed to the WAVEDERIVENL function
    TYPE wavederivenl_args
       TYPE(atoms), POINTER :: CATOM                  !< atom datatype defined above
       REAL(q), POINTER :: POT(:)        !< self consistency potential
       INTEGER, POINTER :: IL            !< index of the quantum number
    END TYPE

    TYPE (atoms), ALLOCATABLE, TARGET, SAVE :: ATOM_LCAO(:)  ! atom descriptor

    REAL(q), ALLOCATABLE :: RCUT(:)      ! real space cutoffs of atoms

    INTEGER, PARAMETER :: NGMAX=1000     ! number of points in interpolating spline function WFCT_G

  CONTAINS

    SUBROUTINE LCAO_INIT(P,INFO,IU0,IU6)
      USE pseudo
      USE constant
      USE prec
      USE cl
      USE base 
      USE ini
      USE lattice
      USE poscar
      USE wave
      USE setexm
      USE nonlr_struct_def
#ifdef VASP_HDF5
      USE vhdf5
#endif
      TYPE (potcar),TARGET :: P(:)
      TYPE (info_struct) INFO

      INTEGER :: IU0, IU6
      !
      ! local variables
      !
      TYPE (potcar), POINTER    :: PP
      TYPE (atoms), POINTER     :: CATOM
      TYPE (rgrid)              :: FR

      INTEGER, PARAMETER :: MAX_ITER = 500
      INTEGER :: I, J, ION, ISPIN, NDIM, YLMMAX
      INTEGER :: ITER, MINITER

      REAL(q), ALLOCATABLE :: EOLDOUT(:), POT(:,:,:), POTPS_SAVE(:)
      REAL(q), ALLOCATABLE :: RHOOUT(:,:,:), RHOIN(:,:,:), RHOTMP(:,:,:)
      REAL(q), ALLOCATABLE :: CRHODEIN(:,:), CRHODEOUT(:,:), CRHODETMP(:,:)
      REAL(q), ALLOCATABLE :: RHOLM(:)
      REAL(q) :: RCOVI,  RCHECK, EDCT, EDCH, EDCT1, EDCTXC, EDCHXC, EDCT1XC
      REAL(q) :: ALPHA                ! mixing parameter for linear mixing
      REAL(q) :: R, ETOT, ETOTOLD, ENCORR, ENMAX
      CHARACTER(30) :: FORM

      REAL(q) R0,R1,oo,VADD

      ALLOCATE(ATOM_LCAO(1:SIZE(P)))

      ENMAX=INFO%ENMAX

      IF (IU0>=0) THEN
         WRITE(IU0,*)'==========================================='
         WRITE(IU0,*)'           LCAO calculation'
         WRITE(IU0,*)'==========================================='
      ENDIF
      IF (IU6>=0) THEN
         WRITE(IU6,*)
         WRITE(IU6,*)'Starting ATOMIC calculations'
         WRITE(IU6,*)
         WRITE(IU6,'(A,20(X,A2))')'Following elements found in POTCAR:',P(:)%ELEMENT
         WRITE(IU6,*)
      ENDIF
      !      
      !  dim. of spin = 0 -> no magnetic field
      !
      ISPIN = 1
      ITER = 0
      !============================================================================
      !
      !
      ! cycle through all elements found in POTCAR
      !
      !
      !============================================================================
      !
      ! reads in the atomic specific cut off radii, should find its place in the 
      ! INCAR file, same ordering as in POTCAR required
      ! 
      ALLOCATE(RCUT(SIZE(P))); RCUT=10._q
      
      ! OPEN(123,FILE='RCUTCAR')
      ! DO I=1,SIZE(P)
      !    READ(123,*)RCUT(I)
      ! ENDDO
      ! CLOSE(123)
      
!      print *,rcuts

      DO ION=1,SIZE(P)   
         PP=> P(ION)
         CATOM => ATOM_LCAO(ION)
         IF (IU6>=0) THEN
            WRITE(IU6,'(/A,A2)') 'Calculating: ',PP%ELEMENT
            WRITE(IU6,*)'Starting DFT run on atom'
         ENDIF
         !
         CALL PUSH_XC_TYPE(PP%LEXCH,1.0_q,1.0_q,1.0_q,1.0_q, 0.0_q)
         ALLOCATE(POTPS_SAVE(SIZE(PP%POTPS))); POTPS_SAVE=PP%POTPS
         !     
         ! Mixing parameter
         !
         !  with 0.5 Co_pv_GW LDA does not convrge !
         !
         ALPHA=0.1d0
         !
         ! Minimum iterations with fixed atomic occupancies (better for convergence)
         !
         MINITER=3
         !
         ! determine extended grid by imposing the constraint that the outermost
         ! grid point is larger 200 A
         !
         CALL DETERMINE_GRIDS(PP,CATOM,IU0,IU6)
         !
         ! Confinement potential parameter (should be the covalent radius)
         !
         ! Al = 118 pm = 1.18 A
         ! Si = 111 pm = 1.11 A
         ! Ti = 132 pm = 1.32 A
         !     
         ! Have a look at which radius wavefunction is nearly zero
         ! and whether we need confining potential
         !
         ! CATOM%RCOVI = (1._q/(1.85_q*1.18_q))
         !
         CATOM%RCOVI = 0.0_q
!        CATOM%RCOVI = (1._q/(0.8*1.11_q))
         CATOM%RCOVI_SHIFT = 0.0_q
!        CATOM%RCOVI_SHIFT =CATOM%FR%R(CATOM%FR%NMAX)**2*CATOM%RCOVI
         !
         ! seek number of occupied one electron states 
         ! using occupancies in PP descriptor and allocate
         ! arrays which depend on the l-quantum numbers
         !
         CALL SEEK_OCCUPIED_VALENCE_STATES(PP,CATOM,IU0,IU6)
         !
         ! Quantity which is important for the density arrays due to the LM's 
         !
         YLMMAX = (CATOM%MAXL + 1)**2
         !
         ! Allocate temporary arrays for calculations on each element
         !
         ALLOCATE(RHOOUT(CATOM%FR%NMAX,YLMMAX,1), RHOIN(CATOM%FR%NMAX,YLMMAX,1),RHOTMP(CATOM%FR%NMAX,YLMMAX,1))
         ALLOCATE(CRHODEOUT(CATOM%LMMAX,CATOM%LMMAX),CRHODEIN(CATOM%LMMAX,CATOM%LMMAX),CRHODETMP(CATOM%LMMAX,CATOM%LMMAX))
         ALLOCATE(POT(CATOM%FR%NMAX,YLMMAX,1),EOLDOUT(1:CATOM%OCCUPIED_VALENCE_STATES),RHOLM(CATOM%LMMAX*CATOM%LMMAX))
         !
         ! get initial potential for atom on extended grid
         !
         CALL EXTEND_ARRAYS(PP,CATOM,IU0,IU6)
         !
         ! now we have to transform the projection operators to a logarithmic grid, 
         ! due to the fact, that they are stored on a linear grid in the POTCAR file
         ! and the atom solver needs it on a logarithmic grid
         !
         CALL LCAO_PROJ_LINTOLOG(PP,CATOM,IU0,IU6)
         !
         ! calculate atomic occupation number in the form used by RAD_POT
         !
         CRHODEIN=0._q
         CALL LCAO_SET_CRHODE_ATOM(CRHODEIN(:,:), CATOM, PP)
         RHOLM=0._q
         CALL LCAO_TRANS_RHOLM(CRHODEIN(:,:),RHOLM(:), CATOM)
         !
         ! Set the matching point of the wavefunctions to the maximum r for the non local 
         ! contributions (HAS TO BEEN CHECKED WHICH IS BEST CHOICE: 1, 1.5)
         !
         CATOM%IRCHECK(:) =  CATOM%PSRMAX(:)
         DO I=1,CATOM%OCCUPIED_VALENCE_STATES
            CATOM%RCHECK(I)=CATOM%FR%R(CATOM%IRCHECK(I))
         ENDDO

         CALL LCAO_SET_DIJ(CATOM,PP,CRHODEIN,RHOLM,EDCH,EDCT1,EDCHXC,EDCT1XC,IU0, IU6)

         EOLDOUT = CATOM%E
         CALL CALCULATE_WAVEFUNCTIONS(CATOM,0,IU0,IU6)

         ! Save the initial wavefunctions
         CATOM%WFCT_L(:,3,:) = CATOM%WFCT_L(:,2,:)

         RHOIN=0._q
         CALL LCAO_CHARGE(CATOM, RHOIN(:,1,1),1,IU0,IU6)
        
         CRHODEOUT=CRHODEIN

         ETOT = 0
         !
         ! Self consistency loop
         !
         IF (IU0>=0) WRITE(IU0,*) 'solving self-consistent calculation'
         DO ITER=1,MAX_ITER

            RHOOUT = 0._q
            CALL LCAO_CHARGE(CATOM, RHOOUT(:,1,1),1,IU0,IU6)
            !
            !
            !
            IF (ITER>MINITER) THEN
               CRHODEOUT=0._q
               CALL LCAO_RECAL_CRHODE(CATOM,CRHODEOUT)
            ENDIF
            !
            ! density mixing
            !
            RHOTMP = RHOIN*(1._q-ALPHA) + RHOOUT*ALPHA
            !
            ! mixing of the occupancies of the augmentation channels
            ! 
            CRHODETMP = CRHODEIN*(1._q-ALPHA) + CRHODEOUT*ALPHA
            !
            ! Copy density due to the fact that we need it for mixing with 
            ! augmentation
            !
            RHOIN = RHOTMP
            CRHODEIN = CRHODETMP
            !
            RHOLM=0._q
            CALL LCAO_TRANS_RHOLM(CRHODETMP(:,:),RHOLM(:), CATOM)
            !     
            ! augment charge
            !
            CALL RAD_AUG_CHARGE(RHOTMP(:,:,1), CATOM%FR, RHOLM, CATOM%LMAX, CATOM%LPS,  &
                 CATOM%MAXL, CATOM%AUG, PP%QPAW)
            !
            ! Calculate potential
            !
            CATOM%POT = 0._q
            POT = 0._q
            CALL RAD_POT(CATOM%FR, ISPIN, 0, 0, .FALSE.,  &
                 RHOTMP, CATOM%RHOPS, CATOM%POT, POT, EDCT, EDCTXC)

            DO I=1,CATOM%FR%NMAX
               CATOM%POT(I) =-POT(I,1,1)/(2._q*SQRT(PI)) &
                    - CATOM%FR%R(I)**2._q*CATOM%RCOVI + CATOM%RCOVI_SHIFT

!              R1=RCUT(ION); R0=R1-1._q; oo=150._q
!              VADD=MIN(EXP(-(R1-R0)/(CATOM%FR%R(I)-R0))/(R1-CATOM%FR%R(I)),oo)
!              IF (CATOM%FR%R(I)< R0) VADD=0._q
!              IF (CATOM%FR%R(I)>=R1) VADD=oo
!              CATOM%POT(I) =-POT(I,1,1)/(2._q*SQRT(PI))-VADD
            ENDDO

            IF (ITER>MINITER) CALL LCAO_SET_DIJ(CATOM,PP,CRHODETMP,RHOLM,EDCH,EDCT1,EDCHXC,EDCT1XC,IU0, IU6)

            CALL CALCULATE_WAVEFUNCTIONS(CATOM,ITER,IU0,IU6)

            IF (CATOM%OCCUPIED_VALENCE_STATES>9+1) THEN
               WRITE(FORM,'(A6,I2,A6)') '(I3,A,',CATOM%OCCUPIED_VALENCE_STATES+2,'F15.8)'         
            ELSE
               WRITE(FORM,'(A6,I1,A6)') '(I3,A,',CATOM%OCCUPIED_VALENCE_STATES+2,'F15.8)'         
            ENDIF

            ETOTOLD = ETOT
            ETOT = 0
            DO I=1,CATOM%OCCUPIED_VALENCE_STATES
               ETOT = ETOT + CATOM%OCC(I)*(2*CATOM%L(I,1)+1)*(CATOM%E(I)+CATOM%RCOVI_SHIFT)
            ENDDO
            !
            ! DEXCCORE only shift see paw.F
            !
            ETOT = ETOT+EDCT+EDCH-EDCT1-CATOM%DEXCCORE
            IF (IU6>=0.AND.(MOD(ITER,10)==0.OR.ITER==1)) WRITE(IU6,'(I3,A2,2F15.8)')ITER,': ',ETOT,DABS(ETOT-ETOTOLD)
            !
            ! 1.E-5 is to small for some elements 
            !
            IF(DABS(ETOT-ETOTOLD)<1E-7.AND.ITER>MINITER) EXIT
            IF (MOD(ITER,1)==0) THEN
               IF(IU6>0) WRITE(IU6,'(I4,A2,2F18.10)') ITER,': ',ETOT,DABS(ETOT-ETOTOLD)
               IF(IU0>0) WRITE(IU0,'(I4,A2,2F18.10)') ITER,': ',ETOT,DABS(ETOT-ETOTOLD)
            ENDIF
         ENDDO
         IF (ITER>MAX_ITER) call vtutor%error('LCAO self-consistent calculation failed')
         !
         ! Calculate total energy of atom by summation of eigenvalues and 
         ! double counting correction
         !
         ETOT = 0
         DO I=1,CATOM%OCCUPIED_VALENCE_STATES
            ETOT = ETOT + CATOM%OCC(I)*(2*CATOM%L(i,1)+1)*(CATOM%E(I)+CATOM%RCOVI_SHIFT)
         ENDDO
         !
         !
         !
         IF (CATOM%OCCUPIED_VALENCE_STATES+4>9) THEN
            WRITE(FORM,'(A1,I2,A7)') '(',CATOM%OCCUPIED_VALENCE_STATES+1,'E20.9)'         
         ELSE
            WRITE(FORM,'(A1,I1,A7)') '(',CATOM%OCCUPIED_VALENCE_STATES+1,'E20.9)'         
         ENDIF
         IF (IU6>=0) THEN
            WRITE(IU6,'(/A,F15.6)')' sum of eigenenergies : ',ETOT
            WRITE(IU6, '(A,F15.6)')' dc corr. E~          : ',EDCT
            WRITE(IU6, '(A,F15.6)')' DEXCCORE             : ',CATOM%DEXCCORE
            WRITE(IU6, '(A,F15.6)')' dc corr. E^          : ',EDCH
            WRITE(IU6, '(A,F15.6)')' dc corr. E~1         : ',EDCT1
            WRITE(IU6, '(A,F15.6)')' TOTAL ENERGIE        : ',ETOT+EDCT+EDCH-EDCT1-CATOM%DEXCCORE
            WRITE(IU6,'(A,F15.6/)')' Diff to POTCAR       : ',DABS( ETOT+EDCT+EDCH-EDCT1-CATOM%DEXCCORE+PP%EATOM)
         ENDIF
         IF (IU0>=0) THEN
            WRITE(IU0,'(A,F15.9)')' Finished DFT run with diff. to POTCAR of: ',DABS( ETOT+EDCT+EDCH-EDCT1-CATOM%DEXCCORE+PP%EATOM)
            DO I=1,CATOM%OCCUPIED_VALENCE_STATES
               WRITE(IU0,*)'eigenenergy: ',CATOM%E(I)
            ENDDO
         ENDIF
         !
         ! set cutoff of support grid 
         !
         CATOM%SUPGRIDMAX=MAXVAL(RCUT)
         CALL DET_CUTOFF(CATOM,RCUT(ION),IU0,IU6)
         ! 
         ! approximate radial functions by Bessel functions
         !
         ALLOCATE(CATOM%WFCT(CATOM%NPLIN,5,CATOM%OCCUPIED_VALENCE_STATES), &
              CATOM%WKIN(CATOM%NPLIN,5,CATOM%OCCUPIED_VALENCE_STATES))
         !
         CALL LCAO_FIND_OPT_BESBASIS_LIN(CATOM,ENMAX,ENCORR,IU0,IU6)
         ALLOCATE(CATOM%WFCT_G(NGMAX,5,CATOM%OCCUPIED_VALENCE_STATES))
         !
         CALL LCAO_WFCT_L_TO_G(CATOM,ENMAX,IU0,IU6)
         !
         ! deallocate local arrays
         !
         DEALLOCATE(RHOOUT, RHOIN, RHOTMP)
         DEALLOCATE(CRHODEOUT,CRHODEIN,CRHODETMP)
         DEALLOCATE(POT,EOLDOUT,RHOLM)
         !
         PP%POTPS=POTPS_SAVE; DEALLOCATE(POTPS_SAVE)
         CALL POP_XC_TYPE
         !
      ENDDO

      IF (IU0>=0) THEN
         WRITE(IU0,*) 'Atomic calculation finished in '//str(ITER)//' iterations'
         WRITE(IU0,*)
      ENDIF

#ifdef VASP_HDF5
      ! write data in hdf5
      CALL VH5_WRITE_LCAO(IH5OUTFILEID, ATOM_LCAO)
#endif

    END SUBROUTINE LCAO_INIT

!**********************************************************************
!
!      LCAO_FIND_OPT_BESBASIS_LIN
!
! routine which finds the optimal number of basis functions according to
! the supplied energy cutoff
!
! Finally we get a almost "smooth" LCAO-basis, very similar to the original 
! wavefunctions, but without any kinks or other nasty things.
!
!**********************************************************************

    SUBROUTINE LCAO_FIND_OPT_BESBASIS_LIN(CATOM,ENMAX_ENTRY,ENCORR,IU0,IU6)

      USE constant

      TYPE (atoms), POINTER :: CATOM
      REAL (q) :: ENMAX_ENTRY
      REAL (q) :: ENCORR
      INTEGER :: IU0,IU6
      ! local variables
      REAL(q), ALLOCATABLE :: BBASIS(:,:,:),KBBASIS(:,:,:),KBBASIS_PEN(:,:,:)
      REAL(q), ALLOCATABLE :: E_BESBAS(:)

      INTEGER :: NBF(1:CATOM%OCCUPIED_VALENCE_STATES)
      INTEGER :: MAXNBF

      INTEGER :: L,I

      REAL(q) :: ENMAX,DERIV
      !
      ! from now on one could use a smaller grid with catom%irmax gridpoints, because
      ! bessel functions have beyond this point zero value and derivative.
      !
      ! this is NOT TESTED or implemented until now
      !
      ENMAX=ENMAX_ENTRY/2._q ! /1.5
      DO L=1,CATOM%OCCUPIED_VALENCE_STATES
         ! HM: I lower 2._q*ENMAX to ENMAX otherwise the calculations fail when the cutoff is too large
         CALL DET_QBFZERO(CATOM,L,ENMAX,NBF(L),IU0,IU6)
         CATOM%NBF(L) = NBF(L)
      ENDDO

      IF (IU0>=0) WRITE(IU0,*)'Performing besselization of the atomic orbitals'
      IF (IU6>=0) THEN
         WRITE(IU6,*)'BESSELIZATION'
         WRITE(IU6,*)'============='
         WRITE(IU6,*)
      ENDIF
      
      MAXNBF = MAXVAL(NBF)

      
      ALLOCATE(E_BESBAS(CATOM%OCCUPIED_VALENCE_STATES))
      ALLOCATE(BBASIS(CATOM%OCCUPIED_VALENCE_STATES,MAXNBF,CATOM%NPLIN))
      ALLOCATE(KBBASIS(CATOM%OCCUPIED_VALENCE_STATES,MAXNBF,CATOM%NPLIN))
      ALLOCATE(KBBASIS_PEN(CATOM%OCCUPIED_VALENCE_STATES,MAXNBF,CATOM%NPLIN))
      ALLOCATE(CATOM%Q(MAXNBF+1,CATOM%OCCUPIED_VALENCE_STATES),CATOM%A(MAXNBF+1,MAXNBF,CATOM%OCCUPIED_VALENCE_STATES))
      ALLOCATE(CATOM%B(MAXNBF,CATOM%OCCUPIED_VALENCE_STATES))

      ENCORR = 0._q
      E_BESBAS =0._q
      CATOM%A=0._q
      CATOM%B=0._q
      CATOM%Q=0._q
      BBASIS = 0._q
      KBBASIS=0._q
      KBBASIS_PEN=0._q

      DO L=1,CATOM%OCCUPIED_VALENCE_STATES

         CALL SETUP_LDEP_BESBASIS_LIN(CATOM,NBF(L),L,CATOM%Q(:,L),CATOM%A(:,:,L),BBASIS,KBBASIS,KBBASIS_PEN,ENMAX,IU0,IU6)

         IF (IU6>=0) THEN
            WRITE(IU6,*)'LOGARITHMIC'
            WRITE(IU6,*)'ENM: ',ENMAX
            WRITE(IU6,*)'NBF: ',NBF(L)
            WRITE(IU6,'(A,I2)')' L = ',CATOM%L(L,1)
            WRITE(IU6,*)'------'
            WRITE(IU6,'(A,I2,A)')' Using ',NBF(L),' LCBF basis functions'
            WRITE(IU6,'(A /, (10F10.4))') ' reciprocal lattice vectors ',CATOM%Q(1:NBF(L)+1,L)
            WRITE(IU6,'(A /, (10F10.4))') ' kinetic energy             ',(CATOM%Q(1:NBF(L)+1,L)*CATOM%Q(1:NBF(L)+1,L))*HSQDTM
            WRITE(IU6,*) 
         ENDIF

         CALL LCAO_WFCT_TO_BESBASIS_LIN(CATOM,NBF(L),L,BBASIS,KBBASIS,KBBASIS_PEN,CATOM%B(:,L),E_BESBAS,IU0,IU6)

         IF (IU6>=0) THEN
            WRITE(IU6,'(A,F15.9)')' Full-DFT energy: ',CATOM%E(L)
            WRITE(IU6,'(A,F15.9)')' LCBF energy:     ',E_BESBAS(L)
            WRITE(IU6,'(A,F15.9)')' Difference:      ',CATOM%E(L)-E_BESBAS(L)
            WRITE(IU6,*) 
         ENDIF
         IF (IU0>=0) WRITE(IU0,'(A,F15.9)')' LCBF energy:     ',E_BESBAS(L)

         !
         ! add the energy gained by confining the atom in a smaller sphere by smaller rcut to INFO%EALLAT
         ! 
         ! mind that this is only somehow a first order correction, due to the fact that we 
         ! calculcate the besselized wavefunctions by a single shot on the full relaxed atomic 
         ! potential. Therefore we won't get the atomic limit if the confinement is rather strong
         ! 
         ENCORR = ENCORR + CATOM%OCC(L)*(CATOM%E(L)-E_BESBAS(L))*(2*CATOM%L(L,1)+1)
         
         CALL LCAO_WFCT_TO_LIN_BESBASIS(CATOM,NBF(L),L,BBASIS,KBBASIS,CATOM%B(:,L),IU0,IU6)

         IF (CATOM%L(L,1)/=1) THEN
            DERIV=0._q
         ELSE
            DERIV=(-2._q*CATOM%WFCT(3,2,L)+16._q*CATOM%WFCT(2,2,L))/(12*(CATOM%WFCT(2,1,L)-CATOM%WFCT(1,1,L)))
         ENDIF
         CALL SPLCOF(CATOM%WFCT(1,1,L),CATOM%NPLIN,CATOM%NPLIN,DERIV)
      ENDDO
      
      DEALLOCATE(BBASIS,KBBASIS,KBBASIS_PEN)
      DEALLOCATE(E_BESBAS)


    END SUBROUTINE LCAO_FIND_OPT_BESBASIS_LIN

!*******************************************************************
!> @brief set up type FR (rgrid) of ATOM so that the last grid
!>   point is beyond 200 Angstroem. This should be large enough.
!>   Hard coded.
!*******************************************************************

  SUBROUTINE DETERMINE_GRIDS(PP, CATOM, IU0, IU6)

    USE prec
    
    TYPE (potcar), POINTER :: PP
    TYPE (atoms), POINTER  :: CATOM

    INTEGER :: IU0, IU6
    INTEGER :: FNMAX, I
    REAL(q), PARAMETER :: RMAX = 200._q
    REAL(q) :: FREND, R

    FREND = PP%R%REND 
    FNMAX = PP%R%NMAX

    IF (FNMAX==500.OR.FNMAX==2000) THEN
       ! full grids normally used in the pseudopotential generator
       IF (IU0>=0) WRITE(IU0,*)' Full grid in POTCAR'
    ELSE
       DO WHILE (FREND < RMAX)
          FREND = FREND*PP%R%D
          FNMAX = FNMAX + 1
       ENDDO
    ENDIF
    !
    ! Allocate arrays depending on the radial grid
    !
    ALLOCATE(CATOM%FR%R(FNMAX),CATOM%FR%SI(FNMAX),CATOM%POTPS(FNMAX), &
         CATOM%POTPSC(FNMAX), CATOM%RHOPS(FNMAX), CATOM%RHO(FNMAX), & 
         CATOM%POT(FNMAX), & 
         CATOM%POTC(FNMAX), CATOM%POTAE(FNMAX),CATOM%POTAEC(FNMAX), &
         CATOM%RHOAE(FNMAX))
    
    CATOM%FR%RSTART = PP%R%RSTART
    CATOM%FR%RMAX = PP%R%RMAX
    CATOM%FR%H = PP%R%H
    CATOM%FR%D = PP%R%D
    CATOM%FR%REND = FREND
    CATOM%FR%NMAX = FNMAX
    CATOM%FR%R(1:PP%R%NMAX) = PP%R%R(1:PP%R%NMAX)
    
    DO I=PP%R%NMAX+1,CATOM%FR%NMAX
       CATOM%FR%R(I) = CATOM%FR%R(I-1)*CATOM%FR%D
    ENDDO
    !
    ! align grid to avoid rounding errors and set up the weights for
    ! the simpson integration in CATOM%FR%SI
    !
    CALL RAD_ALIGN(CATOM%FR,CATOM%FR%RMAX)

    IF (IU6>=0) THEN
       WRITE(IU6,'(A)')
       WRITE(IU6,'(A)') 'Properties of new grid'
       WRITE(IU6,'(A)') '----------------------'
       WRITE(IU6,'(A,I4,A,F6.2,A,I4,A,F6.2,A)') 'Grid extended from ',PP%R%NMAX,' gridpoints (', &
            PP%R%REND,' A) to ',CATOM%FR%NMAX,' gridpoints (',CATOM%FR%REND,' A)'

       WRITE(IU6,'(A,F15.9)') 'RSTART: ', CATOM%FR%RSTART
       WRITE(IU6,'(A,F15.9)') 'RMAX:   ', CATOM%FR%RMAX
       WRITE(IU6,'(A,F15.9)') 'H:      ', CATOM%FR%H
       WRITE(IU6,'(A,F15.9)') 'D:      ', CATOM%FR%D
       WRITE(IU6,'(A,F15.9)') 'REND:   ', CATOM%FR%REND
       WRITE(IU6,'(A,I5)')    'NMAX:   ', CATOM%FR%NMAX
    ENDIF

  END SUBROUTINE DETERMINE_GRIDS

!*******************************************************************
!
!      SEEK_OCCUPIED_VALENCE_STATES
! determines a lot of different stuff like quantum numbers,
! occupation,... , and allocates all necessary arrays
! 
!*******************************************************************

  SUBROUTINE SEEK_OCCUPIED_VALENCE_STATES(PP, CATOM, IU0, IU6 )

    USE prec

    TYPE (potcar), POINTER :: PP
    TYPE (atoms), POINTER :: CATOM
    
    INTEGER :: IU0, IU6
    INTEGER :: I, J, K, L
    INTEGER :: CHANNEL, OCCUPIED_VALENCE_STATES
    !
    ! determine number of occupied valence states for array 
    ! allocation
    !
    CATOM%ELEMENT = PP%ELEMENT

    OCCUPIED_VALENCE_STATES = 0
    
    DO CHANNEL=1,PP%LMAX
       IF (PP%QATO(CHANNEL,CHANNEL) /=0) THEN
          OCCUPIED_VALENCE_STATES=OCCUPIED_VALENCE_STATES+1
       ENDIF
    ENDDO
    CATOM%OCCUPIED_VALENCE_STATES=OCCUPIED_VALENCE_STATES
    
    ALLOCATE(CATOM%N(CATOM%OCCUPIED_VALENCE_STATES), CATOM%L(CATOM%OCCUPIED_VALENCE_STATES,4), &
         CATOM%OCC(CATOM%OCCUPIED_VALENCE_STATES), CATOM%E(CATOM%OCCUPIED_VALENCE_STATES), &
         CATOM%POTCAR_E(CATOM%OCCUPIED_VALENCE_STATES), &
         CATOM%WFCT_L(CATOM%FR%NMAX,5,CATOM%OCCUPIED_VALENCE_STATES), & 
         CATOM%PSRMAX(CATOM%OCCUPIED_VALENCE_STATES),CATOM%NBF(CATOM%OCCUPIED_VALENCE_STATES))
    !
    ! determine pseudo-n, l quantum number and the occupation of each channel
    !
    CATOM%N = -1
    OCCUPIED_VALENCE_STATES = 0
    DO CHANNEL = 1,PP%LMAX
       IF (PP%QATO(CHANNEL,CHANNEL) /= 0) THEN
          OCCUPIED_VALENCE_STATES = OCCUPIED_VALENCE_STATES + 1
          CATOM%L(OCCUPIED_VALENCE_STATES,1) = PP%LPS(CHANNEL)
          CATOM%L(OCCUPIED_VALENCE_STATES,4) = CHANNEL
          IF (OCCUPIED_VALENCE_STATES-1 > 0) THEN
             IF (CATOM%L(OCCUPIED_VALENCE_STATES,1) == CATOM%L(OCCUPIED_VALENCE_STATES-1,1)) THEN
                CATOM%N(OCCUPIED_VALENCE_STATES) = CATOM%N(OCCUPIED_VALENCE_STATES-1) + 1
                CATOM%L(OCCUPIED_VALENCE_STATES,2) = CATOM%L(OCCUPIED_VALENCE_STATES-1,2)
             ELSE
                CATOM%N(OCCUPIED_VALENCE_STATES) = 1
                CATOM%L(OCCUPIED_VALENCE_STATES,2) = CHANNEL
             ENDIF
          ELSE
             CATOM%N(OCCUPIED_VALENCE_STATES) = 1
             CATOM%L(OCCUPIED_VALENCE_STATES,2) = CHANNEL
          ENDIF
          CATOM%OCC(OCCUPIED_VALENCE_STATES) = PP%QATO(CHANNEL,CHANNEL)
          CATOM%E(OCCUPIED_VALENCE_STATES) = PP%E(CHANNEL)! - CATOM%RCOVI_SHIFT
          CATOM%POTCAR_E(OCCUPIED_VALENCE_STATES) = PP%E(CHANNEL)
       ENDIF
    ENDDO
    !
    ! Only if all energy entries are all zero exit, because elements like Rb_pv have a energy of
    ! approx 0.00x in the l=2 channel so in the POTCAR the entry is zero. 
    ! in l
    !
    DO CHANNEL=1,OCCUPIED_VALENCE_STATES
       IF(PP%E(CHANNEL)/=0._q) EXIT
    ENDDO
    IF (CHANNEL==OCCUPIED_VALENCE_STATES+1) THEN
       CALL vtutor%error("No energies are in the POTCAR file --> get new one with them \n EXTING")
    ENDIF
    !
    ! Determine end numbers
    !
    DO CHANNEL=1,CATOM%OCCUPIED_VALENCE_STATES
       DO I=CATOM%L(CHANNEL,2),PP%LMAX
          IF(PP%LPS(I)/=CATOM%L(CHANNEL,1)) THEN
             CATOM%L(CHANNEL,3)=I-1
             EXIT
          ELSE
             CATOM%L(CHANNEL,3)=I
          ENDIF
       ENDDO
    ENDDO
    !
    ! Set LMAX to maxium number in CATOM%L(:,3), otherwise
    ! a projector of a high l quantum number can be used, 
    ! where the wavefunction is not calculated for the 
    ! minimal basis set
    !
    CATOM%LMAX=MAXVAL(CATOM%L(1:CATOM%OCCUPIED_VALENCE_STATES,3))
    CATOM%MAXL=MAXVAL(CATOM%L(:,1))

    IF ( ASSOCIATED(PP%QPAW) ) THEN
       CATOM%MAXL=2*CATOM%MAXL
    ENDIF

    CATOM%ATOMLMMAX=0
    DO I=1,OCCUPIED_VALENCE_STATES
       CATOM%ATOMLMMAX = CATOM%ATOMLMMAX + 2*CATOM%L(I,1)+1
    ENDDO
    !
    ! Allocate l-dependent quantities
    !
    ALLOCATE(CATOM%DIJ(CATOM%LMAX,CATOM%LMAX), CATOM%DION(CATOM%LMAX, CATOM%LMAX), &
         CATOM%QIJ(CATOM%LMAX,CATOM%LMAX), CATOM%QION(CATOM%LMAX, CATOM%LMAX), CATOM%LPS(CATOM%LMAX), &
         CATOM%POTCAR_NODES(CATOM%OCCUPIED_VALENCE_STATES),CATOM%LCAO_NODES(CATOM%OCCUPIED_VALENCE_STATES), &
         CATOM%WPS(CATOM%FR%NMAX,CATOM%LMAX), &
         CATOM%WAE(CATOM%FR%NMAX,CATOM%LMAX),  &
         CATOM%WFRMAXL(CATOM%OCCUPIED_VALENCE_STATES), &
         CATOM%IRCHECK(CATOM%OCCUPIED_VALENCE_STATES),CATOM%RCHECK(CATOM%OCCUPIED_VALENCE_STATES), &
         CATOM%AUG(CATOM%FR%NMAX,0:CATOM%MAXL), CATOM%PROJLOG(CATOM%FR%NMAX,CATOM%LMAX))

    CATOM%LPS(:) = PP%LPS(1:CATOM%LMAX)
    CATOM%LMMAX=0
    DO I=1,CATOM%LMAX
       CATOM%LMMAX = CATOM%LMMAX + 2*CATOM%LPS(I)+1
    ENDDO

    IF (IU6>=0) THEN
       WRITE(IU6,'(A)')
       WRITE(IU6,'(A)') 'l-properties of potential'
       WRITE(IU6,'(A)') '----------------------'
       WRITE(IU6,'(A,I5)') 'LMAX:      ', CATOM%LMAX
       WRITE(IU6,'(A,I5)') 'LMMAX:     ', CATOM%LMMAX
       WRITE(IU6,'(A,I5)') 'ATOMLMMAX: ', CATOM%ATOMLMMAX
       WRITE(IU6,'(A,I5)') 'MAXL:      ', CATOM%MAXL
       WRITE(IU6,'(A)')
    ENDIF
    !
    ! Determine number of nodes of pseudo wavefunction in the POTCAR file (still workaround)
    !
    DO I=1,CATOM%OCCUPIED_VALENCE_STATES
       DO J=1,CATOM%LMAX
          IF (CATOM%L(I,4)==J) THEN
             CATOM%POTCAR_NODES(I)=0
             DO K=2,PP%R%NMAX
                IF(PP%WPS(K,J)*PP%WPS(K-1,J)<0._q) CATOM%POTCAR_NODES(I)=CATOM%POTCAR_NODES(I)+1
             ENDDO
          ENDIF
       ENDDO
    ENDDO

    IF (IU6>=0) THEN
       WRITE(IU6,'(A)')'    l  |  lbeg  lend  wave nodes'
       WRITE(IU6,'(A)')' ------+-------------------------'
       DO CHANNEL=1,OCCUPIED_VALENCE_STATES
          WRITE(IU6,'(I5,A,4I5)') CATOM%L(CHANNEL,1),'  |',CATOM%L(CHANNEL,2:4), CATOM%POTCAR_NODES(CHANNEL)
       ENDDO
    ENDIF

    CATOM%QION(:,:) = PP%QION(1:CATOM%LMAX,1:CATOM%LMAX)
    CATOM%DION(:,:) = PP%DION(1:CATOM%LMAX,1:CATOM%LMAX)
    CATOM%ZVALF = PP%ZVALF

    IF (IU6>=0) THEN
       WRITE(IU6,*)
       WRITE(IU6,'(A,F5.2)')'Valence: ',CATOM%ZVALF
       WRITE(IU6,*)
       WRITE(IU6,'(A)')'        |  n    E     occ.'
       WRITE(IU6,'(A)')' -------+----------------------'
       DO CHANNEL = 1, OCCUPIED_VALENCE_STATES
          WRITE(IU6,'(A,I2,A,I3,F8.2,F9.4)')' l = ',CATOM%L(CHANNEL,1),' |',CATOM%N(CHANNEL), &
               CATOM%E(CHANNEL)+CATOM%RCOVI_SHIFT,CATOM%OCC(CHANNEL)
       ENDDO
       WRITE(IU6,*)
! test
       DO CHANNEL=1,CATOM%LMAX
          WRITE(IU6,'(20F14.7)') CATOM%DION(:,CHANNEL)
       ENDDO
! test
    ENDIF

  END SUBROUTINE SEEK_OCCUPIED_VALENCE_STATES
 
!*******************************************************************
!
! EXTEND_ARRAYS
!
!  copies potentials from PP to CATOMS and extends POTPSC with
!  Coulomb term, and POTPS with a trail 1/r.
! 
!*******************************************************************

  SUBROUTINE EXTEND_ARRAYS(PP,CATOM, IU0, IU6)
    
    USE prec
    USE constant
    USE radial
!   USE setexm
    
    TYPE (potcar), POINTER :: PP
    TYPE (atoms), POINTER :: CATOM
    
    INTEGER :: I, IU0, IU6
    CHARACTER(len=10) :: FORM
    REAL(q) :: A0,A1
#define include_this
#ifdef include_this
    INTEGER LYMAX,LMMAX,K
    REAL(q)  RHOLM(PP%LMMAX*PP%LMMAX)
    OVERLAP  CRHODE(PP%LMMAX,PP%LMMAX,1)
    REAL(q),ALLOCATABLE :: POT(:,:,:), RHO(:,:,:)
    REAL(q) DOUBLEPS_OLD,EXCG_OLD,DEXC

    INTEGER, EXTERNAL :: MAXL1

    LYMAX =MAXL1(PP)*2
    LMMAX=(LYMAX+1)**2

    PP%VPSRMAX=-PP%POTPS(PP%R%NMAX)+PP%POTPSC(PP%R%NMAX)
    ALLOCATE(POT(PP%R%NMAX,LMMAX,1),RHO(PP%R%NMAX,LMMAX,1))
 
    ! set RHOLM and CRHODE to atomic occupancies
    CALL SET_CRHODE_ATOM(CRHODE(:,:,1), PP)
    ! transform CRHODE to llp,LM
    RHOLM=0
    CALL TRANS_RHOLM( CRHODE(:,:,1), RHOLM(:), PP )
    ! get the spin up density
    CRHODE=CRHODE/2

    RHO=0; PP%POTPS=0
    CALL RAD_CHARGE( RHO(:,:,1), PP%R, RHOLM, PP%LMAX, PP%LPS, PP%WPS )
    CALL RAD_AUG_CHARGE(  RHO(:,:,1), PP%R, RHOLM, PP%LMAX, PP%LPS,  &
         LYMAX, PP%AUG, PP%QPAW )

!   CALL PUSH_XC_TYPE(PP%LEXCH,1.0_q,1.0_q,1.0_q,1.0_q, 0.0_q)
    CALL RAD_POT( PP%R, 1, 1, 1, .FALSE.,  &
            RHO, PP%RHOPS, PP%POTPS, POT, DOUBLEPS_OLD, EXCG_OLD)

    DO K=1,PP%R%NMAX
       PP%POTPS(K)=-POT(K,1,1)/(2*SQRT(PI))
    ENDDO

    CALL RAD_CORE_XC( PP%R, PP%RHOAE, DEXC)

    CATOM%DEXCCORE=PP%DEXCCORE+DEXC

!   CALL POP_XC_TYPE

    DEALLOCATE(POT,RHO)
#endif
 
    CATOM%WFCT_L = 0._q
    DO I=1,CATOM%OCCUPIED_VALENCE_STATES
       CATOM%WFCT_L(:,1,I) = CATOM%FR%R(:)
    ENDDO

    CATOM%LEXCH = PP%LEXCH
    CATOM%POTPSC = 0._q
    CATOM%POTPSC(1:PP%R%NMAX) = PP%POTPSC(1:PP%R%NMAX)
    CATOM%POTC = 0._q
    CATOM%POTC(1:PP%R%NMAX) = PP%POTPSC(1:PP%R%NMAX)
    !
    ! shift has to be done due to the fact that the potential of the valence electrons
    ! is recalculated before routine LCAO is called. Propably the fact that on the small grid
    ! the density of valence electrons is not fully described, a shift occurs. In VPSRMAX the original
    ! difference at the boundary between POTPS and POTPSC from the POTCAR file is saved. Due to the fact
    ! that POTPSC is not recalculated it is possible to shift the potential back.
    !
    CATOM%POTPS = 0._q
    CATOM%POTPS(1:PP%R%NMAX) = PP%POTPS(1:PP%R%NMAX)-PP%VPSRMAX-PP%POTPS(PP%R%NMAX)+PP%POTPSC(PP%R%NMAX)
    CATOM%POT = 0._q
    CATOM%POT(1:PP%R%NMAX) = CATOM%POTPS(1:PP%R%NMAX)
    CATOM%POTAE = 0._q
    CATOM%POTAE(1:PP%R%NMAX) = PP%POTAE(1:PP%R%NMAX)-PP%VPSRMAX-PP%POTPS(PP%R%NMAX)+PP%POTPSC(PP%R%NMAX)
    CATOM%RHO = 0._q
    CATOM%RHO(1:PP%R%NMAX) = PP%RHOPS(1:PP%R%NMAX)
    CATOM%RHOPS = 0._q
    CATOM%RHOPS(1:PP%R%NMAX) = PP%RHOPS(1:PP%R%NMAX)
    CATOM%RHOAE = 0._q
    CATOM%RHOAE(1:PP%R%NMAX) = PP%RHOAE(1:PP%R%NMAX)

    CATOM%AUG=0._q
    DO I=0,CATOM%MAXL
       CATOM%AUG(1:PP%R%NMAX,I) = PP%AUG(1:PP%R%NMAX,I)
    ENDDO
    !
    ! POTPS is extended by -1/r+a0/exp(x)**a1
    !    coefficients are choosen to match the two last
    !    points of the known potentials
    ! POTPSC by Coulomb term
    !
    A1=LOG((CATOM%POTPS(PP%R%NMAX)+PP%ZVALF*FELECT/CATOM%FR%R(PP%R%NMAX))/(CATOM%POTPS(PP%R%NMAX-1)+PP%ZVALF*FELECT/CATOM%FR%R(PP%R%NMAX-1)))/(CATOM%FR%R(PP%R%NMAX-1)-CATOM%FR%R(PP%R%NMAX))
    A0=(CATOM%POTPS(PP%R%NMAX-1)+PP%ZVALF*FELECT/CATOM%FR%R(PP%R%NMAX-1))*exp(CATOM%FR%R(PP%R%NMAX-1))**a1
    

    A0 = (CATOM%POTPS(PP%R%NMAX))*CATOM%FR%R(PP%R%NMAX)

    DO I = PP%R%NMAX+1,CATOM%FR%NMAX
       CATOM%POTC(I) = -PP%ZVALF*FELECT/CATOM%FR%R(I)
       CATOM%POT(I) = A0/CATOM%FR%R(I)
       CATOM%RHO(I)=0._q         
       CATOM%AUG(I,:)=0._q
    ENDDO
    
    DO I = 1,CATOM%FR%NMAX
       CATOM%POT(I) =  CATOM%POT(I) - CATOM%FR%R(I)**2._q*CATOM%RCOVI + CATOM%RCOVI_SHIFT!/CATOM%FR%R(I)
    ENDDO

    IF (IU0>=0) WRITE(IU0,*)'extending potentials, augmentation charges and core density to full grid'
    !
    ! copy pseudo and all electron partial waves to arrays, necessary for on site terms
    !
    CATOM%WPS=0._q
    CATOM%WAE=0._q
    CATOM%WPS(1:PP%R%NMAX,:) = PP%WPS(1:PP%R%NMAX,1:CATOM%LMAX)
    CATOM%WAE(1:PP%R%NMAX,:) = PP%WAE(1:PP%R%NMAX,1:CATOM%LMAX)


#ifdef debug
    WRITE(FORM,'(A1,I1,A7)') '(',CATOM%LMAX+1,'F20.10)'         

    OPEN(123,file='wps.dat')
    WRITE(123,FORM)(CATOM%FR%R(I),CATOM%WPS(I,1:CATOM%LMAX),I=1,CATOM%FR%NMAX)
    CLOSE(123)

    OPEN(1234,FILE='arrays_extended.dat')
    WRITE(1234,'(8F20.10)') (CATOM%FR%R(I),CATOM%POT(I),CATOM%POTAE(I),CATOM%POTC(I),-CATOM%POT(I)+CATOM%POTC(I), &
         CATOM%AUG(I,0),CATOM%RHOPS(I),CATOM%RHO(I),I=1,CATOM%FR%NMAX)
    CLOSE(1234)
#endif
    
  END SUBROUTINE EXTEND_ARRAYS

!***********************************************************************
!
!     LCAO_PROJ_LINTOLOG
! Transforms the projection operator stored on a linear
! grid to a logarithmic grid with continuation beyond
! the stored values of the POTCAR file
!
!***********************************************************************

  SUBROUTINE LCAO_PROJ_LINTOLOG(PP,CATOM,IU0,IU6)
    
    USE prec
    USE cl

    TYPE (potcar), POINTER :: PP
    TYPE (atoms), POINTER :: CATOM
    
    INTEGER I, J, K, PMAX
    INTEGER IU0, IU6
      
    REAL(q) :: TEMP(CATOM%FR%NMAX,2)

    REAL(q) :: DUMMY
    CHARACTER(30) :: FORM
    !      
    !     Calculation of the spline coefficients is already done in 
    !     pseudo.F and saved in the PSPNRL array. 
    !
#ifdef debug
    IF (CATOM%LMAX+1>9) THEN
       WRITE(FORM,'(A1,I2,A7)') '(',CATOM%LMAX+1,'E20.10)'         
    ELSE
       WRITE(FORM,'(A1,I1,A7)') '(',CATOM%LMAX+1,'E20.10)'         
    ENDIF

    OPEN(1234,FILE='projlin.dat')
    DO I=1,NPSRNL
       WRITE(1234,FORM) PP%PSPRNL_ORIG(I,1,1),PP%PSPRNL_ORIG(I,2,1:CATOM%LMAX)*PP%PSPRNL_ORIG(I,1,1)
    ENDDO
    CLOSE(1234)
#endif
    IF (IU6>=0) WRITE(IU6,*)'calculating proj. operator on log. grid for ',PP%ELEMENT
    !
    !      Search for the point in the logarithmis grid which corresponds to the 
    !      last point on the linear grid
    !
    CATOM%PROJLOG = 0._q
      
    DO J=1,CATOM%LMAX    
       DO PMAX=1,PP%R%NMAX
          IF(CATOM%FR%R(PMAX) > PP%PSPRNL_ORIG(NPSRNL,1,J)) EXIT
       ENDDO
       PMAX=PMAX-1

       DO I=1,PMAX
          CALL SPLVAL(CATOM%FR%R(I),CATOM%PROJLOG(I,J),DUMMY,PP%PSPRNL_ORIG(1,1,J),NPSRNL,NPSRNL)
          !
          ! Due to the storage layout of the projection operators
          ! we have to multiply with r (radial SE)
          !
          CATOM%PROJLOG(I,J) = CATOM%PROJLOG(I,J)*CATOM%FR%R(I)
       ENDDO

    ENDDO
    !
    ! Determine matching point, which is the r where the 
    ! influence of the projector vanishes
    !
    CATOM%PSRMAX=0
    DO J=1,CATOM%OCCUPIED_VALENCE_STATES
       DO K=CATOM%L(J,2),CATOM%L(J,3)
          DO I=CATOM%FR%NMAX,1,-1
             IF (CATOM%PROJLOG(I,K)>1.E-10) THEN
                CATOM%PSRMAX(J)=MAX(CATOM%PSRMAX(J),I)
                EXIT
             ENDIF
          ENDDO
       ENDDO
    ENDDO

#ifdef debug
    IF (CATOM%LMAX+1>9) THEN
       WRITE(FORM,'(A1,I2,A7)') '(',CATOM%LMAX+1,'E20.10)'         
    ELSE
       WRITE(FORM,'(A1,I1,A7)') '(',CATOM%LMAX+1,'E20.10)'         
    ENDIF

    OPEN(1234,FILE='projlog.dat')
    DO I=1,CATOM%FR%NMAX
       WRITE(1234,FORM) CATOM%FR%R(I),CATOM%PROJLOG(I,:)
    ENDDO
    CLOSE(1234)
#endif      
    
  END SUBROUTINE LCAO_PROJ_LINTOLOG
 
!*******************************************************************
!
! set the array CRHODE such that it corresponds
! to the atomic reference occupancies
! 
!*******************************************************************

    SUBROUTINE LCAO_SET_CRHODE_ATOM(CRHODE, CATOM, PP)
      USE prec
      USE pseudo
      REAL(q) CRHODE(:,:)

      TYPE (atoms),POINTER:: CATOM
      TYPE (potcar), POINTER :: PP
      INTEGER LM, LL, LHI, LOW, MMAX, L, LP, M

      CRHODE=0

      LOW=1
      LM =1
      block: DO
         LL=CATOM%LPS(LOW)
         ! search block with same L
         DO LHI=LOW,CATOM%LMAX
            IF (LL/=CATOM%LPS(LHI)) EXIT
         ENDDO
         LHI=LHI-1
         MMAX=2*LL+1

         DO L =LOW,LHI
         DO LP=LOW,LHI
            DO M =0,MMAX-1
               CRHODE(LM+(L-LOW)*MMAX+M,LM+(LP-LOW)*MMAX+M)=PP%QATO(L,LP)
            ENDDO
         ENDDO
         ENDDO
      
         ! set new LOW value and LM value and go on
         LM=LM+(LHI-LOW+1)*MMAX
         LOW=LHI+1
         IF (LOW > CATOM%LMAX) EXIT block
      ENDDO block
    END SUBROUTINE LCAO_SET_CRHODE_ATOM

!*******************************************************************
!
!  transform the real part of the occupancies RHO(lm,l'm') 
!  to RHO(ll',L,M) using Clebsch Gordan coefficients'
!
!' RHO(ll',LM) =  sum C(LM,ll',mm')  RHO(lm,l'm')
!  where C(LM,ll',mm') = \int Y_LM Y_lm Y_l'm' d Omega
!
!  the storage layout of RHO(llp,LM) is somewhat akward
!  for each l lp pair, (2l+1) (2lp+1) elements must be stored
!  they are stored in the order
!    Lmin,M=0 ... Lmin,M=2*Lmin+1
!     ...
!    Lmax,M=0 ... Lmax,M=2*Lmax+1
!  where Lmin and Lmax are given by the triangular rule
!  Lmin = | l-l' |  and Lmax = | l+l'|
!  certain elements in this array will be always zero, because
!  the sum rule allows only L=Lmin,Lmin+2,...,Lmax
!
!*******************************************************************
    SUBROUTINE LCAO_TRANS_RHOLM( RHOLLMM, RHOLM, CATOM)
      USE pseudo
      USE asa
      TYPE (atoms),POINTER :: CATOM

      REAL(q) RHOLLMM(:,:)   ! net augmentation charge
      REAL(q) RHOLM(:)       ! local charge for each L,M
    ! local varible
      INTEGER CH1,CH2,LL,LLP,LM,LMP,LMINDX,ISTART,IEND,IC,M,MP
      INTEGER IBASE,JBASE,LMIN,LMAX
      REAL(q) FAKT
   ! loop over all channels (l,epsilon)
      IBASE=0

      LM=1
      DO CH1=1,CATOM%LMAX
      LMP=LM
      DO CH2=CH1,CATOM%LMAX

         ! quantum numbers l and lp of these two channels
         LL =CATOM%LPS(CH1)
         LLP=CATOM%LPS(CH2)

         CALL YLM3LOOKUP(LL,LLP,LMINDX)

         ! Lmin and Lmax
         LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)
         
         ! transform coefficients
         FAKT=1
         IF (CH1 /= CH2) THEN
            FAKT=FAKT*2
         ENDIF
         ! JS(IC) is pointing to L*L+M+1, we must subtract LMIN*LMIN
         JBASE=IBASE-LMIN*LMIN
         
         RHOLM(IBASE+1:IBASE+(2*LL+1)*(2*LLP+1))=0
#ifdef debug
         WRITE(0,*) 'CALC_RHOLM: RHOLLMM is',CH1,CH2
         DO MP=1,2*LLP+1
            WRITE(0,'(10F10.7)') (RHOLLMM(LM+M-1,LMP+MP-1),M=1,2*LL+1)
         ENDDO
#endif
         DO M =1,2*LL+1
            DO MP=1,2*LLP+1
               LMINDX=LMINDX+1
               
               ISTART=INDCG(LMINDX)
               IEND  =INDCG(LMINDX+1)
               DO  IC=ISTART,IEND-1
                  RHOLM(JS(IC)+JBASE)= RHOLM(JS(IC)+JBASE)+ YLM3(IC)*FAKT* &
                       &         RHOLLMM(LM+M-1,LMP+MP-1)
               ENDDO
            ENDDO
         ENDDO

#ifdef debug
      WRITE(0,*) 'CALC_RHOLM: augmentation charges are',CH1,CH2
      DO MP=LMIN,LMAX
         WRITE(0,'(I3,10F10.7)') MP,(RHOLM(JBASE+MP**2+M),M=1,MP*2+1)
      ENDDO
#endif

      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO


    END SUBROUTINE LCAO_TRANS_RHOLM

!*******************************************************************
!
!    LCAO_SET_DIJ
! Set up DIJ including on site terms
!
!*******************************************************************
    SUBROUTINE LCAO_SET_DIJ(CATOM,PP,CRHODE,RHOLM,EDCH,EDCT1,EDCHXC,EDCT1XC,IU0, IU6)

      USE prec
      USE constant
      USE radial
      USE pseudo 

      TYPE (potcar), POINTER :: PP
      TYPE (atoms), POINTER :: CATOM     
      
      REAL(q) :: DIJ(CATOM%LMMAX,CATOM%LMMAX)
      REAL(q) :: DIJSMALL(CATOM%LMAX,CATOM%LMAX)
      REAL(q) :: POT(CATOM%FR%NMAX,(CATOM%MAXL + 1)**2,1)
      REAL(q) :: FPOT(CATOM%FR%NMAX)
      REAL(q) :: CRHODE(:,:), RHOLM(:)
      REAL(q) :: DLM(CATOM%LMMAX*CATOM%LMMAX)
      REAL(q) :: RHO(CATOM%FR%NMAX,(CATOM%MAXL + 1)**2,1)
      REAL(q) :: EDCH, EDCT1, EDCHXC, EDCT1XC
      INTEGER :: IU0, IU6
      ! local variables
      INTEGER, EXTERNAL :: MAXL_AUG
      INTEGER :: I,J, IOFF, JOFF,LYMAX

      DLM = 0._q
      CATOM%DIJ = CATOM%DION
      CATOM%QIJ = CATOM%QION

      FPOT(:) = (-CATOM%POT(:) + CATOM%POTC(:))

      CALL LCAO_RAD_POT_WEIGHT(CATOM%FR,FPOT)
      !
      ! \int v^\tilde_eff(r) Q^\hat^L_ij(r) d^3r
      !
      CALL LCAO_RAD_AUG_PROJ(CATOM,FPOT,CATOM%FR,CATOM%DIJ, CATOM%AUG, PP%QPAW(:,:,0),IU0,IU6 )
      !
      ! On site pseudo and all electron terms, wavefunctions contribute through the occupancy
      ! of each augmentation channel CRHODE.
      !
      ! First pseudo partial part
      !
      ! <phi^~_i|v^~_eff^1 - veff^~^1_atom|phi^~_j> - \int (v^~_eff^1 - veff^~^1_atom )Q^\hat^L_ij(r) d^3r
      !
      ! Recalculate n^~^1
      !
      RHO = 0._q
      CALL LCAO_CHARGE(CATOM, RHO(:,1,1), 3, IU0, IU6, CRHODE,PP%R%NMAX)
      !     
      ! augment charge rho=n^~^1 + n^~
      !
      CALL RAD_AUG_CHARGE(RHO(:,:,1), PP%R, RHOLM, CATOM%LMAX, CATOM%LPS,  &
           CATOM%MAXL, CATOM%AUG, PP%QPAW)
      !
      ! Calculate effective potential
      !
      FPOT = 0._q
      POT=0._q
      CALL RAD_POT(PP%R, 1, 0, 0, .FALSE., RHO, PP%RHOPS, FPOT, POT, EDCT1, EDCT1XC)
      !
      ! shift has to be applied, because not the full charge
      ! is in the sphere
      !
      POT(1:PP%R%NMAX,1,1)=-POT(1:PP%R%NMAX,1,1)/(2._q*SQRT(PI))-PP%VPSRMAX-PP%POTPS(PP%R%NMAX)+PP%POTPSC(PP%R%NMAX)
#ifdef debug
      OPEN(123,FILE='potps.dat')
      DO I=1,PP%R%NMAX
         WRITE(123,'(4F20.10)') CATOM%FR%R(i),POT(i,1,1),CATOM%POTPS(i),POT(i,1,1)-CATOM%POTPS(i)
      ENDDO
      CLOSE(123)
#endif


      DO I=1,PP%R%NMAX
         POT(I,1,1)=POT(I,1,1)-CATOM%POTPS(I) &
              - CATOM%FR%R(I)**2._q*CATOM%RCOVI + CATOM%RCOVI_SHIFT!/CATOM%FR%R(I)
      ENDDO
      !
      ! weight potential
      !
      CALL RAD_POT_WEIGHT(PP%R,1,0,POT)
      !
      ! -<phi^~_i|v^~_eff^1 - veff^~^1_atom|phi^~_j> 
      !
      CALL RAD_PROJ(POT(:,:,1),PP%R,-1._q,DLM,CATOM%LMAX,CATOM%LPS,CATOM%WPS)
      !
      ! - \int (v^~_eff^1 - veff^~^1_atom )Q^\hat^L_ij(r) d^3r
      !
      CALL RAD_AUG_PROJ(POT(:,:,1), PP%R, DLM, CATOM%LMAX, CATOM%LPS, &
           CATOM%MAXL, CATOM%AUG, PP%QPAW )
      !
      ! Second all electronpartial part
      !
      ! <phi_i|v_eff^1 - veff^1_atom|phi_j> - \int (v_eff^1 - veff^1_atom )Q^\hat^L_ij(r) d^3r
      !
      ! Recalculate n^1
      !
      CALL LCAO_CHARGE(CATOM, RHO(:,1,1), 2, IU0, IU6, CRHODE,PP%R%NMAX)
      !
      ! Calculate effective potential
      !
      FPOT = 0._q
      POT = 0._q
      CALL RAD_POT(PP%R, 1, 0, 0, .FALSE., RHO, PP%RHOAE, FPOT, POT, EDCH, EDCHXC)
      !
      ! shift has to be applied, because not the full charge
      ! is in the sphere
      !
      POT(1:PP%R%NMAX,1,1)=-POT(1:PP%R%NMAX,1,1)/(2._q*SQRT(PI))-PP%VPSRMAX-PP%POTPS(PP%R%NMAX)+PP%POTPSC(PP%R%NMAX)
#ifdef debug
      OPEN(123,FILE='potae.dat')
      DO I=1,PP%R%NMAX
         WRITE(123,'(4F20.10)') CATOM%FR%R(i),POT(i,1,1),CATOM%POTAE(i),POT(i,1,1)-CATOM%POTAE(i)
      ENDDO
      CLOSE(123)
#endif
      DO I=1,CATOM%FR%NMAX
         POT(I,1,1)=POT(I,1,1)-CATOM%POTAE(I) &
              - CATOM%FR%R(I)**2._q*CATOM%RCOVI + CATOM%RCOVI_SHIFT!/CATOM%FR%R(I)
      ENDDO
      !
      ! weight potential
      !
      CALL RAD_POT_WEIGHT(PP%R,1,0,POT)
      !
      ! +<phi_i|v_eff^1 - v_eff^1_atom|phi_j> 
      !
      CALL RAD_PROJ(POT(:,:,1),PP%R,1._q,DLM,CATOM%LMAX,CATOM%LPS,CATOM%WAE)
      !
      ! map linear orderer DLM to size of DIJSMALL
      !

      DIJ = 0
      CALL LCAO_TRANS_DLM(DIJ,DLM,CATOM)

      DIJSMALL=0._q
      IOFF = 0
      DO I=1,CATOM%LMAX
         JOFF = 0
         DO J=1,CATOM%LMAX
            IF(CATOM%LPS(I)/=CATOM%LPS(J)) CYCLE
            DIJSMALL(I,J) = DIJ(I+IOFF,J+JOFF)
            JOFF = JOFF + 2*CATOM%LPS(J) 
         ENDDO
         IOFF = IOFF + 2*CATOM%LPS(I) 
      ENDDO
      !
      ! add contributions to DIJ
      !
      CATOM%DIJ = CATOM%DIJ+DIJSMALL

    END SUBROUTINE LCAO_SET_DIJ

!*******************************************************************
!
!  LCAO_RAD_POT_WEIGHT
!     same as RAD_POT_WEIGHT but without m-component
!  because POT will be required only for integration we
!  multiply now with the weights
!
!*******************************************************************

  SUBROUTINE LCAO_RAD_POT_WEIGHT(R,POT)
    
    REAL(q) :: POT(:)     ! radial potential
    TYPE (rgrid) :: R
    INTEGER I
    
    DO I=1,R%NMAX
       POT(I)=POT(I)*R%SI(I)
    ENDDO
  END SUBROUTINE LCAO_RAD_POT_WEIGHT

!*******************************************************************
!
!  LCAO_RAD_AUG_PROJ
!     same as RAD_AUG_PROJ but without m-component
!  calculate the integral
!   D(ll'LM) =  \int V(r,L,M) Q(r,L) Q(PAW,ll' L) dr
!  on a radial grid
!  Q(r,L) are the L-dependent 1-normalized compensation charges
!  the potential is given by
!       V(r) =  \sum_lm pot_lm(r) * Y_lm(r)
!  and   pot_lm(r) is stored in POT(2l+1+m,..)
!
!*******************************************************************

  SUBROUTINE LCAO_RAD_AUG_PROJ(CATOM,POT,R,DIJ,AUG,QPAW,IU0,IU6 )
    
    USE constant
    USE prec
    
    TYPE (atoms), POINTER :: CATOM
    TYPE (rgrid) R
    
    INTEGER :: IU0, IU6
    
    REAL(q) :: POT(:)
    REAL(q) :: DIJ(:,:)
    REAL(q) :: AUG(:,0:)    ! 1-normalized L-dep compensation charge
    REAL(q) :: QPAW(:,:)
    !
    ! local variables
    !
    REAL(q) :: RHOLMT,SUM
    INTEGER I, J, K

    !-----------------------------------------------------------------------
    ! first calculate \int V(L,M) Q(L,M)
    !-----------------------------------------------------------------------
    RHOLMT=0._q
    SUM=0._q

    DO I=1,R%NMAX
       SUM=SUM+POT(I)*AUG(I,0)
    ENDDO
    RHOLMT = SUM

    !-----------------------------------------------------------------------
    ! than multiply with QPAW(llp, L) and add the DLM
    !-----------------------------------------------------------------------

    DO I=1,CATOM%OCCUPIED_VALENCE_STATES
       IF (CATOM%N(I)>1) CYCLE
       DO J=CATOM%L(I,2),CATOM%L(I,3)
          DO K=CATOM%L(I,2),CATOM%L(I,3)
             DIJ(J,K) = DIJ(J,K) + RHOLMT*QPAW(J,K)
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE LCAO_RAD_AUG_PROJ

!*******************************************************************
!
!  LCAO_CHARGE calculates all three different kind of charges
!  depending on ITYPE.
!  ITYPE:  1. pseudo-charge density
!          2. all electron partial density
!          3. pseudo partial density
!  For the latter two the occupancy of each augmentation channel
!  has to be supplied (CRHODE), otherwise --> stopping.
!
!*******************************************************************

  SUBROUTINE LCAO_CHARGE(CATOM, RHO, ITYPE, IU0, IU6, CRHODE, PPRMAX)

      USE constant
      USE prec

      TYPE (atoms), POINTER :: CATOM

      REAL(q) :: RHO(:)
      REAL(q), OPTIONAL :: CRHODE(:,:)
      INTEGER :: IU0, IU6
! local variables
      INTEGER :: I, J, N, ITYPE
      INTEGER :: IOFF, JOFF
      INTEGER, OPTIONAL :: PPRMAX
      REAL(q) :: SCALE

      SCALE=1._q/(2*SQRT(PI))
      RHO=0._q

      IF (ITYPE>1.AND.(.NOT.PRESENT(CRHODE).OR..NOT.PRESENT(PPRMAX))) THEN
         CALL vtutor%error("calling LCAO_CHARGE without CRHODE")
      ENDIF

      SELECT CASE (ITYPE)
         CASE(1)
            DO N=1,CATOM%OCCUPIED_VALENCE_STATES
               DO I=1,CATOM%FR%NMAX
                  RHO(I)=RHO(I)+CATOM%WFCT_L(I,2,N)*CATOM%WFCT_L(I,2,N)*CATOM%OCC(N)*(2*CATOM%L(N,1)+1)*SCALE
               ENDDO
            ENDDO
         CASE(2)
            IOFF = 0
            DO I=1,CATOM%LMAX
               JOFF = 0
               DO J=1,CATOM%LMAX
                  DO N=1,PPRMAX
                     RHO(N)=RHO(N)+CATOM%WAE(N,I)*CATOM%WAE(N,J)*CRHODE(I+IOFF,J+JOFF)*(2*CATOM%LPS(I)+1)*SCALE
                  ENDDO
                  JOFF = JOFF+2*CATOM%LPS(J)
               ENDDO
               IOFF = IOFF+2*CATOM%LPS(I)
            ENDDO
         CASE(3)
            IOFF = 0
            DO I=1,CATOM%LMAX
               JOFF = 0
               DO J=1,CATOM%LMAX
                  DO N=1,PPRMAX
                     RHO(N)=RHO(N)+CATOM%WPS(N,I)*CATOM%WPS(N,J)*CRHODE(I+IOFF,J+JOFF)*(2*CATOM%LPS(I)+1)*SCALE
                  ENDDO
                  JOFF = JOFF+2*CATOM%LPS(J)
               ENDDO
               IOFF = IOFF+2*CATOM%LPS(I)
            ENDDO
      END SELECT

    END SUBROUTINE LCAO_CHARGE

!*******************************************************************
!
!  transform D(llp,L,M) to the representation D(lm,l'm')
!  using Clebsch Gordan coefficients and add to another array
!
!  D(lm,l'm') =  sum C(LM,ll',mm') D(llp,LM)
!  where C(LM,ll',mm') = \int Y_LM Y_lm Y_l'm' d Omega
!
!  the storage layout of D(llp,LM) is somewhat akward see above
!
!*******************************************************************


    SUBROUTINE LCAO_TRANS_DLM( DLLMM, DLM, CATOM)
      USE pseudo
      USE asa
      USE constant
      TYPE (atoms), POINTER :: CATOM
      REAL(q) DLLMM(:,:)   ! net augmentation charge
      REAL(q) DLM(:)       ! local charge for each L,M
    ! local varible
      INTEGER CH1,CH2,LL,LLP,LM,LMP,LMINDX,ISTART,IEND,IC,M,MP
      INTEGER IBASE,JBASE,LMIN,LMAX,INMIN,INMAX

   ! loop over all channels (l,epsilon)
      IBASE=0

      LM=1
      DO CH1=1,CATOM%LMAX
      LMP=LM
      DO CH2=CH1,CATOM%LMAX

   ! quantum numbers l and lp of these two channels
      LL =CATOM%LPS(CH1)
      LLP=CATOM%LPS(CH2)

      CALL YLM3LOOKUP(LL,LLP,LMINDX)
   ! Lmin and Lmax
      LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)

   ! JS(IC) is pointing to L*L+M+1, we must subtract LMIN*LMIN
      JBASE=IBASE-LMIN*LMIN
#ifdef debug
      DO M=LMIN,LMAX
         WRITE(0,'(3I3,10F14.7)') CH1,CH2,M, &
               (DLM(JBASE+M*M+MP),MP=1,2*M+1)
      ENDDO
#endif

      INMIN=1000
      INMAX =-1000

      DO M =1,2*LL+1
      DO MP=1,2*LLP+1
         LMINDX=LMINDX+1

         ISTART=INDCG(LMINDX)
         IEND  =INDCG(LMINDX+1)
         DO  IC=ISTART,IEND-1
            INMIN=MIN(JS(IC)+JBASE,INMIN)
            INMAX=MAX(JS(IC)+JBASE,INMAX)
            DLLMM(LM+M-1,LMP+MP-1)=DLLMM(LM+M-1,LMP+MP-1)+ &
                            DLM(JS(IC)+JBASE)*YLM3(IC)
         ENDDO
      ENDDO
      ENDDO
  ! fill symmetric components (CH2 is >= CH1)
      IF  (CH1 /= CH2) THEN
         DO M =1,2*LL+1
         DO MP=1,2*LLP+1
            DLLMM(LMP+MP-1,LM+M-1)=DLLMM(LM+M-1,LMP+MP-1)
         ENDDO
         ENDDO
      ENDIF

      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO

    END SUBROUTINE LCAO_TRANS_DLM

!*******************************************************************
!
!> @brief calculations the atomic orbitals for the occupied valence states
!> by inward and outward integration, and matching the first
!> derivative at the maximum r of non local contributions,
!> furthermore the number of nodes is also checked to be right
!
!*******************************************************************
  SUBROUTINE CALCULATE_WAVEFUNCTIONS(CATOM,NITER,IU0,IU6)

    USE prec
    USE constant
    USE cl
    USE root_find, ONLY : ZBRAK, ZBRENT
    USE iso_c_binding

    TYPE (atoms), POINTER :: CATOM
    INTEGER :: NITER
    INTEGER :: IU0, IU6

    ! local variables
    REAL(q),PARAMETER :: TOL=1e-12
    INTEGER :: I, J, K, KJ, KI, IS
    INTEGER,TARGET :: IL
    INTEGER :: ITER, NODES, NERR
    REAL(q),TARGET :: POT(CATOM%FR%NMAX)
    REAL(q) :: WAVEWORK(CATOM%FR%NMAX)
    REAL(q) :: NORM, NORMW, W
    REAL(q) :: A,B,C,FA,FB,FC,DX
    REAL(q), ALLOCATABLE :: WAVEPROJ(:), TEMP(:,:)

    TYPE(C_PTR) :: ARGS
    TYPE(wavederivenl_args),TARGET :: ARGS_LOCAL
    REAL(q) :: X1,X2,X1_TMP,X2_TMP
    ! The following parameters have been carefully tuned by running
    ! all the POTCARS and timing the LCAO calculation.
    ! When changing them make sure that all the POTCARS still work and there is not a big performance penalty
    REAL(q),PARAMETER :: MARGIN=0.5 ! Size of the interval
    !> Number of intervals in which to look for zeros.
    !> energy range is [E-NTRY*MARGIN:E+NTRY*MARGIN]
    INTEGER,PARAMETER :: NTRY=16
    INTEGER,PARAMETER :: N0=10 !< Number of divisions of interval
    !> Number of divisions of the divisions of the interval.
    !> A zero can be found in a interval corresponding to MARGIN can be up to N0*NSAMPLING
    INTEGER,PARAMETER :: NSAMPLING=100
    INTEGER,PARAMETER :: NZEROS_MAX=5 !< Maximum number of zeros that can be found in interval
    INTEGER :: NB,N
    INTEGER :: ZEROS_FOUND, ICLOSEST_ZERO
    REAL(q) :: DELTA, EDIFF, CLOSEST_ZERO
    REAL(q) :: XB1(NZEROS_MAX),XB2(NZEROS_MAX)
    REAL(q) :: E_ZEROS(NZEROS_MAX)

#ifdef debug
    WRITE(IU6,'(A,F6.3)')'    RCOVI: ',CATOM%RCOVI
    WRITE(IU6,'(A,F6.3)')'   RCHECK: ',CATOM%RCHECK
    WRITE(IU6,'(A,I4)')'   NMAX: ',CATOM%FR%NMAX
    WRITE(IU6,'(A,I4)')'   IRCHECK: ',CATOM%IRCHECK
    WRITE(IU6,'(A,I4)')'   LDIM: ',CATOM%LMAX
#endif

    
    POT = 0._q
    DO I=1,CATOM%FR%NMAX
       POT(I)=(CATOM%POTC(I) - CATOM%POT(I))
    ENDDO

!#define debug
    DO IL=1,CATOM%OCCUPIED_VALENCE_STATES
       !
       ! Starting value is the energy given in POTCAR
       !
       C = CATOM%E(IL)
       !
       ! It is important to avoid calculation of an l=3 state
       ! until all other states are sufficiently converged,
       ! otherwise calculation diverges.
       ! HM: Had to increase from 5 to 8 due to LDA/Ac POTCAR
       !
       IF((CATOM%L(IL,1)==3).AND.NITER<8.AND.NITER>0) THEN
          CYCLE
       ENDIF

       !WRITE(*,*) NITER, str(IL)//'/'//str(CATOM%OCCUPIED_VALENCE_STATES), CATOM%L(IL,1), CATOM%POTCAR_NODES(IL), CATOM%E(IL)

       ARGS_LOCAL%CATOM=>CATOM
       ARGS_LOCAL%POT=>POT
       ARGS_LOCAL%IL=>IL
       ARGS = C_LOC(ARGS_LOCAL)
       !A = CATOM%E(IL)
       !FA = FUNC(A,3,ARGS)
       !WRITE(*,*) A, FA

       !A  = CATOM%E(IL)
       !CALL WAVEDERIVENL(CATOM,A,POT,IL,KJ,FA,NODES,NERR)
       !OPEN(UNIT=99,FILE='START_'//str(IL)//'_'//str(NITER))
       !DO I=1,CATOM%FR%NMAX
       !   WRITE(99,*) CATOM%WFCT_L(I,2,IL)
       !ENDDO
       !CLOSE(99)

       ! The problem is that the radial Schrodinger equation has to be solved
       ! using the shooting method (Shooting to a fixed point).
       ! We choose a radius and perform inwards and outward integration and match the value
       ! at this fitting point. The derivative DX1 and DX2 however might not match.
       ! We need to find the energies at which the derivative matches (DX1-DX2==0).
       !
       ! This is a matter of finding the zeros of the function in a robust way.
       !

       ZEROS_FOUND = 0
       DO I=1,NTRY  ! Loop over intervals
          DO IS=1,2 ! Loop over lower or higher interval
             ! Define intervals in which to look for zeros
             IF (I==1) THEN
                ! Since the first two intervals are adjacent we use only one
                X1 = CATOM%E(IL) - I*MARGIN
                X2 = CATOM%E(IL) + I*MARGIN
                N  = N0*2
             ELSE
                IF (IS==1) THEN
                   ! Look in lower interval
                   X1 = CATOM%E(IL) - I*MARGIN
                   X2 = CATOM%E(IL) - (I-1)*MARGIN
                ELSE
                   ! Look in higher interval
                   X1 = CATOM%E(IL) + (I-1)*MARGIN
                   X2 = CATOM%E(IL) + I*MARGIN
                ENDIF
                N  = N0
             ENDIF

             !WRITE(*,*) I,IS,X1,X2
             ! Bracket zero
             ! An energy interval is used to find zeros by the rule of changing sign.
             ! This energy interval DX=(X2-X1)/N0 might be too large and miss some zeros.
             ! to increse the sampling we shift the bounds by DELTA that is smaller than DX
             DX = (X2-X1)/N0
             DO J=0,NSAMPLING-1 ! Loop over sampling
                DELTA = DX/NSAMPLING*J ! here we compute the DELTA
                !WRITE(*,*) 'NSAMPLING',J,DX,DELTA
                NB = NZEROS_MAX
                ! Find zeros
                CALL ZBRAK(FUNC,3,ARGS,X1+DELTA,X2+DELTA,N,XB1,XB2,NB)
                ! If we did not find any zeros then we move to the next sampling point
                IF (NB==0) CYCLE
                ! We did find some zeros
                DO K=1,NB ! Loop over all the zeros that were found
                   !WRITE(*,*) '========'
                   !WRITE(*,*) I,XB1(K),XB2(K)
                   ! Call bisection method to find the zero
                   X1_TMP=XB1(K)
                   X2_TMP=XB2(K)
                   !CALL BISECT(FUNC,3,ARGS,X1_TMP,X2_TMP,TOL,C,FC) ! Timing is 1563s for all 622 POTCARS
                   CALL ZBRENT(FUNC,3,ARGS,X1_TMP,X2_TMP,TOL,C,FC)  ! Timing is 1022s for all 622 POTCARS
                   ! Compute solution
                   CALL WAVEDERIVENL(CATOM,C,POT,IL,KJ,FC,NODES,NERR)
                   ! We found a zero. Check if its really a zero,
                   ! can be a divergent point in which case FC is huge
                   IF (ABS(FC)<1e-3) THEN
                     ! The number of nodes is roughly correct so we exit
                     ! We allow to have +- 1 node with respect to the reference in the POTCAR file
                     IF (ABS(NODES-CATOM%POTCAR_NODES(IL))<2) GOTO 42
                     !IF (NODES==CATOM%POTCAR_NODES(IL)) GOTO 42
                     ! The number of nodes is not correct but we keep track of it
                     IF (ZEROS_FOUND>=NZEROS_MAX) CYCLE
                     ZEROS_FOUND = ZEROS_FOUND + 1
                     E_ZEROS(ZEROS_FOUND) = C
                   ENDIF
                ENDDO
             ENDDO ! Loop over sampling

             IF (I==1) CYCLE ! As mentioned above, the first interval is special
          ENDDO ! Loop over lower and higher interval
       ENDDO ! Loop over intervals
       ! We reached here so it means we failed
       NB=0
       42 CONTINUE

       ! There is nothing we can do!
       IF (NB==0) THEN
#ifdef debug
          ! Write a list of zeros we have found
          DO I=1,ZEROS_FOUND
             CALL WAVEDERIVENL(CATOM,E_ZEROS(I),POT,IL,KJ,FC,NODES,NERR)
             WRITE(*,*) I,E_ZEROS(I),FC,NODES,CATOM%IRCHECK(IL),KJ
             OPEN(UNIT=99,FILE='WF_'//str(I))
             DO J=1,CATOM%FR%NMAX
                WRITE(99,*) CATOM%WFCT_L(J,2,IL)
             ENDDO
             CLOSE(99)
          ENDDO
          ! Write behaviour of the derivative
          OPEN(UNIT=99,FILE='ZEROS')
          X1 = CATOM%E(IL) - NTRY*MARGIN
          X2 = CATOM%E(IL) + NTRY*MARGIN
          N = N0*NSAMPLING*(2*NTRY+1)/1000
          DX=(X2-X1)/N
          C=X1
          WRITE(*,*) X1,X2
          DO I=1,N
             CALL WAVEDERIVENL(CATOM,C,POT,IL,KJ,FC,NODES,NERR)
             WRITE(99,*) C,FC,NODES
             C=C+DX
          ENDDO
          CLOSE(99)
#endif
          CALL vtutor%bug('Found no zeros in interval '//str(X1)//':'//str(X2)//&
                          '.\nThis should not happen, please report indicating which POTCAR you are using.',__FILE__,__LINE__)
       ENDIF

       ! In this check we allow the difference between expected and computed nodes to be 1
       IF (ABS(NODES-CATOM%POTCAR_NODES(IL))>1) THEN
#ifdef debug
          WRITE(*,*) 'Wrong number of nodes NEXPECT: '//str(CATOM%POTCAR_NODES(IL))//' NFOUND: '//str(NODES)
          WRITE(*,*) 'final point', CATOM%IRCHECK(IL), KJ, C, FC
          OPEN(UNIT=99,FILE='FIN_'//str(IL)//'_'//str(NITER))
          DO I=1,CATOM%FR%NMAX
             WRITE(99,*) CATOM%WFCT_L(I,2,IL)
          ENDDO
          CLOSE(99)
#endif
          CALL vtutor%bug('Found no zeros in interval '//str(X1)//':'//str(X2)//&
                          '.\nThis should not happen, please report indicating which POTCAR you are using.',__FILE__,__LINE__)

       ENDIF
       CATOM%LCAO_NODES(IL) = NODES

       !
       ! append at the end the exponential tail
       !
       DO I=KJ+1,CATOM%FR%NMAX
          W=EXP(SQRT(-2._q*(C+CATOM%RCOVI_SHIFT))*(CATOM%FR%R(KJ)-CATOM%FR%R(I)))
          IF (W<(1.E-35_q)/ABS(CATOM%WFCT_L(KJ,2,IL))) THEN
             CATOM%WFCT_L(I:CATOM%FR%NMAX,2,IL)=0
             EXIT
          ENDIF
          CATOM%WFCT_L(I,2,IL)=W*CATOM%WFCT_L(KJ,2,IL)
       ENDDO

       CATOM%E(IL) = C
       !
       ! orthogonalize to the second projector
       ! 
       ALLOCATE(WAVEPROJ(CATOM%L(IL,3)-CATOM%L(IL,2)+1),TEMP(CATOM%FR%NMAX,CATOM%L(IL,3)-CATOM%L(IL,2)+1))       
       !
       ! Normalize wavefunction with the overlap operator and the augmentation charge
       ! 
       !
       WAVEPROJ=0._q

       DO I=CATOM%L(IL,2),CATOM%L(IL,3)
          J = I+1-CATOM%L(IL,2)
          TEMP=0._q
          TEMP(:,J) = CATOM%WFCT_L(:,2,IL)*CATOM%PROJLOG(:,I)
          CALL SIMPI(CATOM%FR,TEMP(:,J),WAVEPROJ(J))          
       ENDDO

       NORM = 0._q
       DO I=CATOM%L(IL,2),CATOM%L(IL,3)
          DO J=CATOM%L(IL,2),CATOM%L(IL,3)
             NORM = NORM + WAVEPROJ(I+1-CATOM%L(IL,2))*CATOM%QION(I,J)* WAVEPROJ(J+1-CATOM%L(IL,2))
          ENDDO
       ENDDO

       TEMP = 0._q
       NORMW = 0._q
       TEMP(:,1) = CATOM%WFCT_L(:,2,IL)*CATOM%WFCT_L(:,2,IL)
       CALL SIMPI(CATOM%FR,TEMP(:,1),NORMW)

       NORM = DSQRT(1._q/(NORM+NORMW))

       CATOM%WFCT_L(:,2,IL) = CATOM%WFCT_L(:,2,IL)*NORM

       DEALLOCATE(TEMP,WAVEPROJ)
    ENDDO ! loop over occupied states

  END SUBROUTINE CALCULATE_WAVEFUNCTIONS

!>@brief Wrapper to WAVEDERIVNL so that it can be used in the bisection routine
  FUNCTION FUNC(X,NARGS,ARGS) RESULT(FX)
    USE iso_c_binding
    REAL(q), INTENT(IN) :: X
    INTEGER, INTENT(IN) :: NARGS
    TYPE(C_PTR), VALUE, INTENT(IN) :: ARGS
    REAL(q) :: FX
    ! local variables
    TYPE (wavederivenl_args),POINTER :: ARGS_LOCAL
    INTEGER :: KJ,NODES,NERR

    CALL C_F_POINTER(ARGS,ARGS_LOCAL)
    CALL WAVEDERIVENL(ARGS_LOCAL%CATOM,X,ARGS_LOCAL%POT,ARGS_LOCAL%IL,KJ,FX,NODES,NERR)
  END FUNCTION FUNC

#ifdef VASP_HDF5
!>@brief Write the LCAO atomic information to hdf5
  subroutine vh5_write_lcao(fileid, atom_lcao, group, subgroup)
    use hdf5
    use vhdf5_base
    use vhdf5
    integer(HID_T), intent(in) :: fileid !< hdf5 file handle
    type (atoms), intent(in)  :: atom_lcao(:) !< Array containing the LCAO information
    character(len=*), intent(in), optional :: group !< group name at the lower level (default: GRP_RESULTS)
    character(len=*), intent(in), optional :: subgroup !< subgroup name (default: GRP_PROJECTORS)
    !-----------------------------------------------------
    integer(HID_T) :: groupid, subgroupid
    integer :: iatom
    character(len=2) :: element
    character(len=40) :: my_group, my_subgroup
    !-----------------------------------------------------
    my_group    = GRP_RESULTS; if (present(group))    my_group    = trim(group)
    my_subgroup = GRP_LCAO;    if (present(subgroup)) my_subgroup = trim(subgroup)
    VH5_CHECK(vh5_group_open_or_create(fileid,  trim(my_group), groupid))
    VH5_CHECK(vh5_group_open_or_create(groupid, trim(my_subgroup), subgroupid))

    do iatom=1,size(atom_lcao)
       element = atom_lcao(iatom)%element
       ! Write potential
       ! Write number of nodes in potcar
       VH5_CHECK(vh5_write(subgroupid,element//"/potcar_nodes",atom_lcao(iatom)%potcar_nodes(:)))
       ! Write number of nodes in LCAO
       VH5_CHECK(vh5_write(subgroupid,element//"/lcao_nodes",atom_lcao(iatom)%lcao_nodes(:)))
       ! Write atomic wavefunctions in the logarithmic grid
       VH5_CHECK(vh5_write(subgroupid,element//"/wfct_l",atom_lcao(iatom)%wfct_l(:,:3,:)))
       ! Write potcar pseudo and ae wavefunctions in the logarithmic grid
       VH5_CHECK(vh5_write(subgroupid,element//"/potcar_r",atom_lcao(iatom)%fr%r(:)))
       VH5_CHECK(vh5_write(subgroupid,element//"/potcar_wps",atom_lcao(iatom)%wps(:,:)))
       VH5_CHECK(vh5_write(subgroupid,element//"/potcar_wae",atom_lcao(iatom)%wae(:,:)))
       ! This is the mapping from LCAO to the POTCAR (because only occupied orbitals are used for LCAO)
       VH5_CHECK(vh5_write(subgroupid,element//"/lcao_l",atom_lcao(iatom)%l(:,:)))
       ! Write atomic energies
       VH5_CHECK(vh5_write(subgroupid,element//"/lcao_energy",atom_lcao(iatom)%e(:)))
       ! Write atomic energies from POTCAR
       VH5_CHECK(vh5_write(subgroupid,element//"/potcar_energy",atom_lcao(iatom)%potcar_e(:)))
    enddo
    !
    ! close/deallocate
    !
    VH5_CHECK(vh5_group_close_writing(subgroupid))
    VH5_CHECK(vh5_group_close_writing(groupid))
    VH5_CHECK(vh5_file_flush(fileid))
  end subroutine vh5_write_lcao
#endif

!****************************************************************************
!
!> @brief Perfroms the inward and outward integration for a given 
!> l-channel and calculate the difference in the first derivative at CATOM%IRCHECK
!
!****************************************************************************
  SUBROUTINE WAVEDERIVENL(CATOM,E,POT,IL,KJ,DX,NODES,NERR)
    
    USE prec
    USE cl
    
    TYPE (atoms), POINTER :: CATOM

    INTEGER :: I
    INTEGER :: IL !< l quantum number for which to compute
    INTEGER :: KI, KJ, NODES, IU0, IU6, NERR
    REAL(q) :: POT(CATOM%FR%NMAX)

    !local variables
    REAL(q) :: E,  DX, ACONV, DX1, DX2
    REAL(q) :: A1(CATOM%FR%NMAX), A2(CATOM%FR%NMAX), TMP(CATOM%FR%NMAX)

    REAL(q), ALLOCATABLE :: WAVEPROJ(:),TEMP(:,:)
    REAL(q) :: NORM, NORMW, INTEGRAL
    INTEGER :: j

    A1 = 0._q
    !
    ! performs the inward integration
    !
    CALL INWINT(E,CATOM%L(IL,1),CATOM%FR,POT,CATOM%IRCHECK(IL)-3,KJ,NODES,A1,A2)
    A2 = 0._q
    !
    ! solves the schroedinger equation with a nonlocal contribution
    !
    CALL NLSOLVE(CATOM,E,IL,POT,A2,NODES,IU0,IU6,NERR)
    !
    ! write both parts of wavefunction in array
    !
    KI=CATOM%IRCHECK(IL)

    ACONV = A2(KI)/A1(KI)
    DO I=1,KI
       CATOM%WFCT_L(I,2,IL) = A2(I)*SIGN(1._q,ACONV)
       !CATOM%WFCT_L(I,2,IL) = A2(I)
    ENDDO
    DO I=KI+1, CATOM%FR%NMAX
       CATOM%WFCT_L(I,2,IL) = A1(I)*ABS(ACONV)
       !CATOM%WFCT_L(I,2,IL) = A1(I)*ACONV
    ENDDO
    !
    ! orthoganlize to the second projector
    ! 
    ALLOCATE(WAVEPROJ(CATOM%L(IL,3)-CATOM%L(IL,2)+1),TEMP(CATOM%FR%NMAX,10))       
       
    WAVEPROJ=0._q
    
    DO I=CATOM%L(IL,2),CATOM%L(IL,3)
       J = I+1-CATOM%L(IL,2)
       TEMP=0._q
       TEMP(:,J) = CATOM%WFCT_L(:,2,IL)*CATOM%PROJLOG(:,I)
       CALL SIMPI(CATOM%FR,TEMP(:,J),WAVEPROJ(J))          
    ENDDO
    
    NORM = 0._q
    DO I=CATOM%L(IL,2),CATOM%L(IL,3)
       DO J=CATOM%L(IL,2),CATOM%L(IL,3)
          NORM = NORM + WAVEPROJ(I+1-CATOM%L(IL,2))*CATOM%QIJ(I,J)* WAVEPROJ(J+1-CATOM%L(IL,2))
       ENDDO
    ENDDO
    
    TEMP = 0._q
    NORMW = 0._q
    TEMP(:,1) = CATOM%WFCT_L(:,2,IL)*CATOM%WFCT_L(:,2,IL)
    CALL SIMPI(CATOM%FR,TEMP(:,1),NORMW)
    
    NORM = DSQRT(1._q/(NORM+NORMW))
    
    CATOM%WFCT_L(:,2,IL) = CATOM%WFCT_L(:,2,IL)*NORM
    
    DEALLOCATE(WAVEPROJ,Temp)
    !
    ! calculate the logarithmic derivatives DX1 and DX2 and the 
    ! difference
    !
    DX1 = (A1(KI-2)+8._q*A1(KI+1)-8._q*A1(KI-1)-A1(KI+2))/&
          CATOM%FR%R(KI)/CATOM%FR%H/12._q/A1(KI)
    DX2 = (A2(KI-2)+8._q*A2(KI+1)-8._q*A2(KI-1)-A2(KI+2))/&
          CATOM%FR%R(KI)/CATOM%FR%H/12._q/A2(KI)

    DX = DX1 - DX2

    NODES=0
    DO I=CATOM%FR%NMAX,2,-1
       IF(CATOM%WFCT_L(I,2,IL)*CATOM%WFCT_L(I-1,2,IL)<0._q) NODES=NODES+1
    ENDDO

  END SUBROUTINE WAVEDERIVENL

!***********************************************************************
!  SUBROUTINE NLSOLVE
!  subroutine solves  the Schroedinger equation for a non-local
!  generalized Vanderbild pseudopotential based on a generalisation
!  of the algorithm given by X. Gonze et.al Phys Rev B44, 8503 (1991)
!   PROJ    projection operatores x_i (in Ry)
!   E       Energy in Ry
!   VR      local Potential in Hartree
!
!***********************************************************************

  SUBROUTINE NLSOLVE(CATOM,E,IL,POT,AHOM,NODES,IU0,IU6,NERR)

    USE cl
    USE prec
    USE constant
    
    TYPE (atoms) :: CATOM

    INTEGER :: IU0, IU6
    INTEGER :: I, J, K, IL, KI, NODES, NERR
    INTEGER :: LMAX, IFAIL
        
    REAL(q) :: E
    REAL(q) :: AHOM(CATOM%FR%NMAX), POT(CATOM%FR%NMAX)
    REAL(q) :: B(CATOM%FR%NMAX), AIHOM(CATOM%FR%NMAX,CATOM%LMAX)
    REAL(q) :: TMP(CATOM%FR%NMAX), ATMP(CATOM%LMAX)

    REAL(q) :: A(CATOM%L(IL,3)-CATOM%L(IL,2)+1), AA(CATOM%L(IL,3)-CATOM%L(IL,2)+1,CATOM%L(IL,3)-CATOM%L(IL,2)+1)
    REAL(q) :: V(CATOM%L(IL,3)-CATOM%L(IL,2)+1)

    KI = 0
    !
    ! for outward integration use three more points, needed for derivatives
    !
    !-----------------------------------------------------------------------
    !  first solve the homogenous schroedinger-equation  (T+V-E) phi  = 0
    !  and the inhomogenous differential equations       (T+V-E) phi_i= x_i
    !-----------------------------------------------------------------------
    !
    ! homogenous
    !
    AHOM = 0._q
    B = 0._q

    NERR = 0
    CALL OUTINT(E,CATOM%L(IL,1),CATOM%FR,POT,KI,NODES,AHOM,B,NERR,RC=CATOM%RCHECK(IL)*CATOM%FR%D**3._q,TOL=1E-8_q)
    AHOM = AHOM/(2._q*RYTOEV)
    !
    ! inhomogenous
    !
    AIHOM = 0._q
    DO I=CATOM%L(IL,2),CATOM%L(IL,3)
       CALL OUTINT(E,CATOM%L(IL,1),CATOM%FR,POT,KI,NODES,AIHOM(:,I),B,NERR,&
                   RC=CATOM%RCHECK(IL)*CATOM%FR%D**3._q,IH=CATOM%PROJLOG(:,I),TOL=1E-8_q)
       AIHOM(:,I)=AIHOM(:,I)/(2._q*RYTOEV)
    ENDDO

    !-----------------------------------------------------------------------
    !  now calculate the vector a_i   = -  (D_ij- E Q_ij ) <x_i|phi>
    !-----------------------------------------------------------------------
    DO I=CATOM%L(IL,2),CATOM%L(IL,3)
       TMP(:) = 0._q
       TMP(:) = CATOM%PROJLOG(:,I)*AHOM(:)
       CALL SIMPI(CATOM%FR,TMP,ATMP(I))
    ENDDO

    A = 0._q
    DO I=CATOM%L(IL,2),CATOM%L(IL,3)
       DO J=CATOM%L(IL,2),CATOM%L(IL,3)
          A(I-CATOM%L(IL,2)+1) = A(I-CATOM%L(IL,2)+1) - (CATOM%DIJ(I,J) - E*CATOM%QIJ(I,J))*ATMP(J)
       ENDDO
    ENDDO
    !-----------------------------------------------------------------------
    !  calculate the matrix    aa_ik = (D_ij- E Q_ij ) <x_j|phi_k> +delta_ik
    !-----------------------------------------------------------------------

    ATMP = 0._q
    AA = 0._q


    DO K=CATOM%L(IL,2),CATOM%L(IL,3)
       DO J=CATOM%L(IL,2),CATOM%L(IL,3)
          TMP=0._q
          TMP(:) = CATOM%PROJLOG(:,J)*AIHOM(:,K)
          CALL SIMPI(CATOM%FR,TMP,ATMP(J))
       ENDDO
       DO I=CATOM%L(IL,2),CATOM%L(IL,3)
          IF (K==I) AA(I-CATOM%L(IL,2)+1,K-CATOM%L(IL,2)+1)= AA(I-CATOM%L(IL,2)+1,K-CATOM%L(IL,2)+1)+1._q
          DO J=CATOM%L(IL,2),CATOM%L(IL,3)
             AA(I-CATOM%L(IL,2)+1,K-CATOM%L(IL,2)+1) = AA(I-CATOM%L(IL,2)+1,K-CATOM%L(IL,2)+1) + (CATOM%DIJ(I,J) - E*CATOM%QIJ(I,J))*ATMP(J)
          ENDDO
       ENDDO
    ENDDO
    !-----------------------------------------------------------------------
    ! for using lapack we have to copy the calculated entries to the
    ! beginning of the matrix AA and vector A, and for obtaining the solution
    ! in the end also AIHOM
    !-----------------------------------------------------------------------

    LMAX = CATOM%L(IL,3)-CATOM%L(IL,2)+1
    !-----------------------------------------------------------------------
    !  solve the linear equation  aa_ik v_k = a_i
    !
    !  uses lapack routines dgetrf and dgetrs      
    !-----------------------------------------------------------------------
    V=0._q
    IFAIL=0

    CALL DGETRF(LMAX,LMAX,AA,LMAX,V,IFAIL)
    IF (IFAIL/=0) THEN
       CALL vtutor%error("NLSOLVE: lapack routine DGETRF FAILED ! IFAIL: " // str(IFAIL) //&
                         "\nCheck if you are using the latest pseudopotential.")
    ENDIF
    IFAIL=0
    CALL DGETRS('N',LMAX,1,AA,LMAX,V,A,LMAX,IFAIL)
    IF (IFAIL/=0) THEN
       CALL vtutor%error("NLSOLVE: lapack routine DGETRS FAILED ! IFAIL: " // str(IFAIL) //&
                         "\nCheck if you are using the latest pseudopotential.")
    ENDIF
   
    V(:)=A(:)
    !-----------------------------------------------------------------------
    !  calculate the solution
    !-----------------------------------------------------------------------

    DO I=1,LMAX
       DO J=1,CATOM%FR%NMAX
          AHOM(J) = AHOM(J) + V(I)*AIHOM(J,I+CATOM%L(IL,2)-1)
       ENDDO
    ENDDO
    RETURN
    
  END SUBROUTINE NLSOLVE

!*******************************************************************
!
!     LCAO_RECAL_CRHODE
! Recalculate the occupation numbers with the newly calculated
! wavefunctions
!
!*******************************************************************
  SUBROUTINE LCAO_RECAL_CRHODE(CATOM,CRHODE)

    USE cl

    TYPE (atoms), POINTER :: CATOM

    REAL(q) :: CRHODE(:,:)
    ! local variables

    INTEGER :: I, J, K, L
    INTEGER :: IINT, JINT
    REAL(q) :: TEMP(CATOM%FR%NMAX),WAVEPROJ(2)
    CRHODE=0._q
    !
    ! Calculation of the array is performed in blocks
    !
    IINT = 0
    DO I=1,CATOM%LMAX
       JINT = 0
       DO J=1,CATOM%LMAX
          DO K=1,CATOM%OCCUPIED_VALENCE_STATES
             IF(CATOM%L(K,1)/=CATOM%LPS(I).OR.CATOM%L(K,1)/=CATOM%LPS(J)) CYCLE
             !
             ! Calculates for a given pair of projectors the intgeral
             ! with the wavefuncitons 
             !      
             WAVEPROJ=0._q
             TEMP(:) = CATOM%WFCT_L(:,2,K)*CATOM%PROJLOG(:,I)
             CALL SIMPI(CATOM%FR,TEMP(:),WAVEPROJ(1))
             TEMP(:) = CATOM%WFCT_L(:,2,K)*CATOM%PROJLOG(:,J)
             CALL SIMPI(CATOM%FR,TEMP(:),WAVEPROJ(2))
             DO L=1,(2*CATOM%L(K,1)+1)
                CRHODE(IINT+L,JINT+L) = CRHODE(IINT+L,JINT+L) + WAVEPROJ(1)*WAVEPROJ(2)*CATOM%OCC(K)
             ENDDO
          ENDDO
          JINT = JINT + 2*CATOM%LPS(J)+1
       ENDDO
       IINT = IINT + 2*CATOM%LPS(I)+1
    ENDDO
    
  END SUBROUTINE LCAO_RECAL_CRHODE

!****************************************************************
!
!         DET_CUTOFF
!
! Determines the cutoff on the radial grid. The desired cut-off has
! to be supplied in the INCAR. With this user defined point,
! the next available gridpoint on the logarithmic grid is used.
!
! Mind that the largest cutoff of all occupied states is 
! taken as the cutoff for the atom under consideration.
!
! Additionally the number of gridpoints on the linear grid,
! where the besselization is performed, is determined. It 
! has turned out that it is a fairly very save choice to 
! use 100 gridpoints/Angstroem. 
! 
! CHANGES:
!   - changed the behaviour for the determination of the outermost
!     gridpoint. Now it is the supplied one. Choosing the
!     the nearest point on the logarithmic grid is not desireable
!     because in the range of 4 Angstroem and more only one or two
!     gridpoints are in the space of 1 Angstroem. So no real choice is 
!     possible.
!
!****************************************************************

    SUBROUTINE DET_CUTOFF(CATOM,RCUT,IU0,IU6)

      TYPE (atoms), POINTER :: CATOM
      INTEGER :: IU0, IU6
      REAL(q) :: RCUT
      !
      ! local variables
      ! 
      INTEGER :: I,J,NMAX
      INTEGER, PARAMETER :: NPA = 100 ! NPA...number of points per Angstroem, 100 seems to be a secure choice here
      !
      ! Find cutoff point
      !
      CATOM%WFRMAX=0._q
      CATOM%WFRMAXL=0._q

      DO I=1,CATOM%OCCUPIED_VALENCE_STATES
         CATOM%WFRMAX = RCUT
         CATOM%WFRMAXL(I)=CATOM%WFRMAX
      ENDDO


      IF (IU6>=0) THEN
         WRITE(IU6,*)
         WRITE(IU6,*)'radial cutoff and # of points:'
         WRITE(IU6,*)'------------------------------'
         WRITE(IU6,'(A,A2,A,F8.4)') ' Maximum cut off for ',CATOM%ELEMENT,' : ',CATOM%WFRMAX
      ENDIF
      IF (IU0>=0) WRITE(IU0,'(A,A2,A,F8.4)') ' Maximum cut off for ',CATOM%ELEMENT,' : ',CATOM%WFRMAX
      !
      ! determine number of necessary grid points on the linear lattice
      !
      NMAX=INT(CATOM%SUPGRIDMAX*NPA)
      ! 
      ! check whether we have an odd number of points, necessary for
      ! later simpson integration. Will converge badly otherwise, and 
      ! will be in principle wrong.
      !
      NMAX=NMAX+(1-MOD(NMAX,2))
      CATOM%NPLIN=NMAX
      !
      IF (IU6>=0) THEN
         WRITE(IU6,'(A,I3,A,I5)')' # of points (',NPA,'/Angstr.):',CATOM%NPLIN
         WRITE(IU6,*)
      ENDIF

    END SUBROUTINE DET_CUTOFF

!**********************************************************************
!
!       DET_QBFZERO
! helper routine which determines the number of basis functions with
! the cutoff energy ENMAX
!
!**********************************************************************

    SUBROUTINE DET_QBFZERO(CATOM,L,ENMAX,NBF,IU0,IU6)

      USE constant

      TYPE (atoms), POINTER :: CATOM
      
      INTEGER :: L
      REAL(q) :: ENMAX
      INTEGER :: NBF
      INTEGER :: IU0,IU6
! local
      INTEGER :: I
      REAL(q) :: QBFZERO(101)
      
      DO I=1,100
         CALL AUG_BEZERO(QBFZERO(:),CATOM%L(L,1),I+1)
         QBFZERO(:)=QBFZERO(:)/CATOM%WFRMAXL(L)
         IF (QBFZERO(I+1)*QBFZERO(I+1)*HSQDTM>ENMAX) EXIT
      ENDDO
      IF(I==100) THEN
         CALL vtutor%error('DET_QBFZERO, number of basis functions is to small for this energy increase')
      ENDIF
      NBF=(I-1)

    END SUBROUTINE DET_QBFZERO

!**********************************************************************
!
!       SETUP_LDEP_BESBASIS_LIN
! routine to set up the basis functions for a specific l-quantum 
! number on the linear grid used for approximating the atomic orbitals. 
! 
! First SETUP_LDEP_BESBASIS_LOG has to be called, because coefficients
! QBFZERO, A, are needed from there.
!
!**********************************************************************

    SUBROUTINE SETUP_LDEP_BESBASIS_LIN(CATOM,NBF,L,QBFZERO,A,BBASIS,KBBASIS,KBBASIS_PEN,ENMAX,IU0,IU6)

      USE radial
      USE constant

      TYPE (atoms), POINTER :: CATOM

      INTEGER, PARAMETER :: M=128        ! # of gridpoints for gauss-legendre integration 
                                         ! has to be bigger then 32, because number of besselfunctions
                                         ! can be large, therefore rapid oszillations occur in the 
                                         ! interval.

      INTEGER :: IU0,IU6                 ! I/O
      INTEGER :: NBF                     ! # of basis functions
      INTEGER :: L                       ! l-quantum number

      REAL(q) :: QBFZERO(:)              ! q's where j_l(q*RMAX) = 0
      REAL(q) :: A(:,:)                  ! coefficients matrix
      REAL(q) :: BBASIS(:,:,:)           ! basis functions on logarithmic grid
      REAL(q) :: KBBASIS(:,:,:)          ! kinetic term of basis functions on logarithmic grid
      REAL(q) :: KBBASIS_PEN(:,:,:)      ! kinetic term of basis functions on logarithmic grid
! local
      REAL(q) :: WE(M),RA(M)             ! weights and abscisses of gauss-legendre integration
      REAL(q) :: BJ(M,NBF+1)             ! values of bessel functions on grid of gauss-legendre
      REAL(q) :: PROJ(1:NBF)             ! projections of wavefunctions onto current one
      REAL(q) :: LINDIST
      REAL(q) :: BJ2,BJP,BJP2
      REAL(q) :: WFCT,WFCT2,WKIN,WKIN_PEN,ENMAX,NORM

      INTEGER :: I,J,N,K,ITYPE,IFAIL

      CHARACTER(30) :: FORM,FILEOUT

      EXTERNAL GAUSSI2                   ! gauss legendre routine for weights and abscisses

      PROJ=0._q
      !
      ! Calculate q's for j_l(q*RMAX) = 0
      !
      CALL AUG_BEZERO(QBFZERO(:),CATOM%L(L,1),NBF+1)
      QBFZERO(:)=QBFZERO(:)/CATOM%WFRMAXL(L)
      !
      ! set up gauss-legendre weights and abscisses
      !      
      ITYPE=0
      CALL GAUSSI(GAUSSI2,0.0_q,CATOM%WFRMAXL(L),ITYPE,M,WE,RA,IFAIL)
      !
      ! calculate for given QBFZERO y-value for besselfunctions
      ! on the abscisses of gauss-legendre
      !
      DO I=1,NBF+1
         DO J=1,M
            CALL SBESSEL(QBFZERO(I)*RA(J),BJ(J,I),CATOM%L(L,1))
         ENDDO
      ENDDO
      !
      ! Initialize basis function 
      !
      nbasis: DO I=1,NBF
         !
         ! Calculate the two last terms of the series. The first is set to one
         ! and the last one so that the derivative at RMAX is zero. Mind that 
         ! we have to multiply the first derivative by the corresponding QBFZERO 
         ! due to the chain rule and the fact that j_l(x) is replaced by
         ! j_l(q_i*r)
         !
         ! !!!!!!!
         !
         ! seems that A(I+1,I) is always one when A(I,I) is one, why ?
         ! not hard coded any more
         ! 
         A(I,I)=1._q
         CALL SBESSE2(QBFZERO(I)*CATOM%WFRMAXL(L),BJ2,BJP,CATOM%L(L,1))
         A(I+1,I)=(-BJP*QBFZERO(I))! *A(I,1) comment in if A(I,J) /= 1._q
         CALL SBESSE2(QBFZERO(I+1)*CATOM%WFRMAXL(L),BJ2,BJP,CATOM%L(L,1))
         A(I+1,I)=A(I+1,I)/(BJP*QBFZERO(I+1))
         !
         ! hardcoded here, otherwise uncomment 4 lines above and comment
         ! the following statement
         ! 
         !            A(I+1,I)=1._q
         IF (I>=2) THEN
            !
            ! Calculate projections of actual wavefunction first
            !
            PROJ=0._q
            DO J=1,NBF-1
               NORM=0._q
               DO K=1,M
                  WFCT=0._q
                  DO N=1,J+1
                     WFCT = WFCT + (A(N,J)*BJ(K,N))
                  ENDDO
                  NORM = NORM + WFCT*(A(I,I)*BJ(K,I) + A(I+1,I)*BJ(K,I+1))*RA(K)**2._q*WE(K)
               ENDDO
               PROJ(J)=NORM
            ENDDO
            !
            ! calculate coefficients
            !
            DO J=1,NBF
               DO K=1,NBF-1
                  A(J,I) = A(J,I) - PROJ(K)*A(J,K)
               ENDDO
            ENDDO
         ENDIF
         !
         ! Calculate Norm 
         !
         NORM=0._q
         DO J=1,M
            WFCT=0._q
            DO K=1,I+1
               WFCT = WFCT + (A(K,I)*BJ(J,K))
            ENDDO
            NORM = NORM + (WFCT*RA(J))**2._q*WE(J)
         ENDDO
         A(:,I) = A(:,I)/SQRT(NORM)
         !
         ! Bring basis functions and corresponding T|sum_j A(j,i)*j_l(q_j*r) to logarithmic grid
         !
         BBASIS(L,I,:)=0._q
         KBBASIS(L,I,:)=0._q
         KBBASIS_PEN(L,I,:)=0._q

         LINDIST=CATOM%SUPGRIDMAX/(CATOM%NPLIN-1)
         DO J=1,CATOM%NPLIN
            IF((J-1)*LINDIST>CATOM%WFRMAXL(L)) EXIT
            WFCT = 0._q
            WKIN = 0._q
            WKIN_PEN = 0._q
            DO N=1,I+1
               CALL SBESSEL(QBFZERO(N)*(J-1)*LINDIST,WFCT2,CATOM%L(L,1))
               WFCT = WFCT + A(N,I)*WFCT2
               WKIN = WKIN + A(N,I)*WFCT2*QBFZERO(N)**2._q
               IF (QBFZERO(N)**2._q*HSQDTM-ENMAX >ENMAX) THEN
                  CALL vtutor%bug("internal error in SETUP_LDEP_BESBASIS_LIN: exceeding bounds", __FILE__, __LINE__)
               ENDIF

               IF (QBFZERO(N)**2._q*HSQDTM>ENMAX) THEN
                  WKIN_PEN = WKIN_PEN + A(N,I)*WFCT2*QBFZERO(N)**2._q*(1._q + TAN((QBFZERO(N)**2._q*HSQDTM-ENMAX)/ENMAX*PI/2._q))
                  
               ELSE
                  WKIN_PEN = WKIN_PEN + A(N,I)*WFCT2*QBFZERO(N)**2._q
               ENDIF
            ENDDO
            WKIN = WKIN*HSQDTM
            WKIN_PEN = WKIN_PEN*HSQDTM
            BBASIS(L,I,J)=WFCT
            KBBASIS(L,I,J)=WKIN
            KBBASIS_PEN(L,I,J)=WKIN_PEN
         ENDDO
      ENDDO nbasis

#ifdef debug
      !
      ! dump into file(s)
      !
      IF (NBF+1>9) THEN
         WRITE(FORM,'(A1,I2,A7)') '(',NBF+1,'F20.9)'         
      ELSE
         WRITE(FORM,'(A1,I1,A7)') '(',NBF+1,'F20.9)'         
      ENDIF
      WRITE(FILEOUT,'(A9,I1,A4)')'besbasis_',CATOM%L(L,1),'.dat'

      LINDIST=CATOM%WFRMAX/(CATOM%NPLIN-1)
      OPEN(123,FILE=FILEOUT)
      DO J=1,CATOM%NPLIN
         WRITE(123,FORM)(J-1)*LINDIST ,BBASIS(L,1:NBF,J)
         IF ((J-1)*LINDIST>CATOM%WFRMAXL(L)) EXIT
      ENDDO
      CLOSE(123)

      WRITE(FILEOUT,'(A10,I1,A4)')'kbesbasis_',CATOM%L(L,1),'.dat'

      LINDIST=CATOM%WFRMAX/(CATOM%NPLIN-1)
      OPEN(123,FILE=FILEOUT)
      DO J=1,CATOM%NPLIN
         WRITE(123,FORM)(J-1)*LINDIST ,KBBASIS(L,1:NBF,J)
         IF ((J-1)*LINDIST>CATOM%WFRMAXL(L)) EXIT
      ENDDO

      CLOSE(123)

      WRITE(FILEOUT,'(A13,I1,A4)')'kbesbasispen_',CATOM%L(L,1),'.dat'

      LINDIST=CATOM%WFRMAX/(CATOM%NPLIN-1)
      OPEN(123,FILE=FILEOUT)
      DO J=1,CATOM%NPLIN
         WRITE(123,FORM)(J-1)*LINDIST ,KBBASIS_PEN(L,1:NBF,J)
         IF ((J-1)*LINDIST>CATOM%WFRMAXL(L)) EXIT
      ENDDO

      CLOSE(123)

#endif

    END SUBROUTINE SETUP_LDEP_BESBASIS_LIN
 
!**********************************************************************
!
!      LCAO_WFCT_TO_BESBASIS_LIN
!
! similar to LCAO_WFCT_TO_BESBASIS but performs everything on the linear
! grid.
! 
!**********************************************************************

    SUBROUTINE LCAO_WFCT_TO_BESBASIS_LIN(CATOM,NBF,L,BBASIS,KBBASIS,KBBASIS_PEN,COEFF,E,IU0,IU6)     

      TYPE (atoms), POINTER :: CATOM

      REAL(q) :: BBASIS(:,:,:)
      REAL(q) :: KBBASIS(:,:,:)
      REAL(q) :: KBBASIS_PEN(:,:,:)
      REAL(q) :: E(:)
      REAL(q) :: COEFF(:)
      INTEGER :: L
      INTEGER :: NBF
      INTEGER :: IU0,IU6
! lapack
      REAL(q) :: EV(NBF)
      REAL(q) :: WORK(1:3*NBF+1)
      INTEGER :: IWORK(1:5*NBF+3)
! local
      REAL(q) :: H(NBF,NBF)
      REAL(q) :: S(NBF,NBF)
      REAL(q) :: PROJ(CATOM%LMAX,NBF)
      REAL(q) :: TEMP(CATOM%nplin)
      REAL(q) :: WFCTEMP,WKINTEMP,NORM
      
      REAL(q) :: SPLARRAY(CATOM%FR%NMAX,5)
      REAL(q) :: LINPOT(CATOM%NPLIN)
      REAL(q) :: LINPROJ(CATOM%NPLIN,CATOM%LMAX)
      REAL(q) :: DERIV,LINDIST,TMP
      REAL(q) :: R(CATOM%NPLIN)
      REAL(q) :: SI(CATOM%NPLIN)
      CHARACTER(LEN=13) :: OUTFILE

      INTEGER :: IERR
      INTEGER :: I,J,M,N,LL

      SPLARRAY=0._q
      !
      ! spline potentials
      !
      SPLARRAY(:,1)=CATOM%FR%R(:)
      SPLARRAY(:,2)=CATOM%POTC(:)-CATOM%POT(:)

      
      DERIV=(SPLARRAY(2,2)-SPLARRAY(1,2))/(SPLARRAY(2,1)-SPLARRAY(1,1))

      CALL SPLCOF(SPLARRAY(1,1),CATOM%FR%NMAX,CATOM%FR%NMAX,DERIV)

      LINDIST=CATOM%SUPGRIDMAX/(CATOM%NPLIN-1)

      R(1)=0._q
      LINPOT(1)=CATOM%POTC(1)-CATOM%POT(1)

      DO J=2,CATOM%NPLIN
         R(J)=(J-1)*LINDIST
         CALL SPLVAL(R(J),LINPOT(J),TMP,SPLARRAY(:,:),CATOM%FR%NMAX,CATOM%FR%NMAX)
      ENDDO

      CALL LIN_SET_SIMP(CATOM%NPLIN,SI,LINDIST)
      !
      ! spline projectors
      !
      DO I=CATOM%L(L,2),CATOM%L(L,3)

         SPLARRAY(:,2)=CATOM%PROJLOG(:,I)
      
         DERIV=(SPLARRAY(2,2)-SPLARRAY(1,2))/(SPLARRAY(2,1)-SPLARRAY(1,1))

         CALL SPLCOF(SPLARRAY(1,1),CATOM%FR%NMAX,CATOM%FR%NMAX,DERIV)
      
         LINPROJ(1,I)=0._q

         DO J=2,CATOM%NPLIN
            LINPROJ(J,I)=(J-1)*LINDIST
            CALL SPLVAL(R(J),LINPROJ(J,I),TMP,SPLARRAY(:,:),CATOM%FR%NMAX,CATOM%FR%NMAX)
         ENDDO
      ENDDO

      H=0._q
      S=0._q
      PROJ=0._q
      TEMP=0._q

      DO I=1,NBF
         DO J=CATOM%L(L,2),CATOM%L(L,3)
            TEMP(:) = LINPROJ(:,J)*BBASIS(L,I,:)*R(:)
            CALL LIN_SIMPI(CATOM%NPLIN,R,SI,TEMP,PROJ(J,I))
         ENDDO
      ENDDO

      DO I=1,NBF
         DO J=1,NBF
            !
            ! explicit calculation of the overlap matrix is absolutely
            ! necessary, because due to rounding errors, which develope
            ! at high NBF, which we also get in other matrix part,
            ! the routine will not be stable if overlap is set to unity.
            !
            TEMP(:) =  BBASIS(L,I,:)*BBASIS(L,J,:)*r(:)**2._q
            CALL lin_SIMPI(catom%nplin,R,si,TEMP(:),S(I,J))
            
            TEMP(:) = (BBASIS(L,I,:)*KBBASIS_PEN(L,J,:) + BBASIS(L,I,:)*linpot(:)*BBASIS(L,J,:))*R(:)**2._q
            CALL lin_SIMPI(catom%nplin,R,si,TEMP(:),H(I,J))
            DO M=CATOM%L(L,2),CATOM%L(L,3)
               DO N=CATOM%L(L,2),CATOM%L(L,3)
                  H(I,J) = H(I,J) + CATOM%DIJ(M,N)*PROJ(M,I)*PROJ(N,J)
                  S(I,J) = S(I,J) + CATOM%QIJ(M,N)*PROJ(M,I)*PROJ(N,J)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DO I=1,NBF
         DO J=I+1,NBF
            H(I,J)=(H(I,J)+H(J,I))/2
            H(J,I)=H(I,J)
         ENDDO
      ENDDO
      !
      ! generalized symmetric eigenvalue problem solved by lapack routine
      !
      IERR=0
      CALL DSYGV(1,'V','U',NBF,H,NBF,S,NBF,EV, WORK, 3*NBF+1, IERR)
      
      IF (IERR/=0) THEN
         CALL vtutor%bug("internal error in LCAO_WFCT_TO_BESBASIS_LIN: DSYGV failed \n IERR: " &
            // str(IERR), __FILE__, __LINE__)
      ENDIF

      LL=0
      DO I=1,L
         IF (CATOM%L(I,1)==CATOM%L(L,1)) LL=LL+1
      ENDDO
!     E(L)=EV(1)
      E(L)=EV(LL)

      !
      ! recalculate the norm of the full wavefunction
      !
      NORM=0._q
      
      DO N=1,CATOM%nplin
         WFCTEMP=0._q
         DO I=1,NBF
!           WFCTEMP = WFCTEMP + BBASIS(L,I,N)*H(I,1)
            WFCTEMP = WFCTEMP + BBASIS(L,I,N)*H(I,LL)
         ENDDO
         TEMP(N)=WFCTEMP
      ENDDO
      
      DO I=CATOM%L(L,2),CATOM%L(L,3)
         CALL lin_SIMPI(catom%nplin,r,si,TEMP(:)*linPROJ(:,I)*R(:),PROJ(I,1))
      ENDDO

      CALL lin_SIMPI(catom%nplin,R,si,TEMP(:)*TEMP(:)*R(:)**2._q,NORM)
      
      DO I=CATOM%L(L,2),CATOM%L(L,3)
         DO J=CATOM%L(L,2),CATOM%L(L,3)
            NORM = NORM + PROJ(I,1)*CATOM%QIJ(I,J)*PROJ(J,1)
         ENDDO
      ENDDO
!     COEFF(:)=H(:,1)/SQRT(NORM)
      COEFF(:SIZE(H,1))=H(:,LL)/SQRT(NORM)

    END SUBROUTINE LCAO_WFCT_TO_BESBASIS_LIN

    SUBROUTINE LIN_SET_SIMP(NMAX,SI,LINDIST)

      REAL(q) :: SI(:)
      REAL(q) :: LINDIST
      INTEGER :: NMAX,K
      
      SI=0._q

      IF(MOD(NMAX,2)/=1) THEN
         CALL vtutor%error("even number if points needed in simpson integration")
      ENDIF

      DO K=NMAX,3,-2
         SI(K)=    LINDIST/3._q +SI(K)
         SI(K-1)=4*LINDIST/3._q
         SI(K-2)=  LINDIST/3._q
      ENDDO

    END SUBROUTINE LIN_SET_SIMP

    SUBROUTINE LIN_SIMPI(NMAX,R,SI,F,RES)
      
      REAL(q) :: R(:),SI(:),F(:)
      REAL(q) :: RES

      INTEGER :: NMAX,K
      
      RES=0._q

      DO K=1,NMAX
         RES=RES+F(K)*SI(K)
      ENDDO

    END SUBROUTINE LIN_SIMPI

!**********************************************************************
!
!      LCAO_WFCT_TO_LIN_BESBASIS
!
! routine which uses the b-basis (LC of bessel functions on linear grid)
! and coefficients calculated by LCAO_WFCT_TO_BESBASIS to calculate
! the basis functions on the linear grid.
!
!**********************************************************************

    SUBROUTINE LCAO_WFCT_TO_LIN_BESBASIS(CATOM,NBF,L,BBASIS,KBBASIS,COEFF,IU0,IU6)
      
      TYPE (atoms), POINTER :: CATOM

      REAL(q) :: BBASIS(:,:,:)
      REAL(q) :: KBBASIS(:,:,:)
      REAL(q) :: COEFF(:)

      INTEGER :: NBF, L, IU0,IU6
! local
      INTEGER :: N, I
      REAL(q) :: WFCT, WKIN, LINDIST, DERIV, PREFACT
      !
      ! Determine the right sign in comparison to the numerical
      ! wavefunction, use second point, because first is for l>0 
      ! equal zero.
      !
      WFCT=0._q
      WKIN=0._q
      DO I=1,NBF
         WFCT = WFCT + BBASIS(L,I,2)*COEFF(I)
         WKIN = WKIN + KBBASIS(L,I,2)*COEFF(I)
      ENDDO
      
      IF (WFCT*CATOM%WFCT_L(2,2,L)>0._q) THEN
         PREFACT=1._q
      ELSE
         PREFACT=-1._q
      ENDIF

      LINDIST=CATOM%SUPGRIDMAX/(CATOM%NPLIN-1)
      DO N=1,CATOM%NPLIN
         WFCT=0._q
         WKIN=0._q
         DO I=1,NBF
            WFCT = WFCT + BBASIS(L,I,N)*COEFF(I)
            WKIN = WKIN + KBBASIS(L,I,N)*COEFF(I)
         ENDDO
         CATOM%WFCT(N,1,L)=(N-1)*LINDIST
         CATOM%WKIN(N,1,L)=(N-1)*LINDIST
         CATOM%WFCT(N,2,L)=PREFACT*WFCT
         CATOM%WKIN(N,2,L)=PREFACT*WKIN
      ENDDO

    END SUBROUTINE LCAO_WFCT_TO_LIN_BESBASIS

!**********************************************************************
!
!      LCAO_WFCT_L_TO_G
!
!**********************************************************************

    SUBROUTINE LCAO_WFCT_L_TO_G(CATOM,ENMAX,IU0,IU6)

      USE constant
      USE radial, ONLY : BESSEL_TRANSFORM_RADIAL_FUNCTION

      TYPE (atoms), POINTER :: CATOM
      REAL(q) :: ENMAX
      INTEGER :: IU0,IU6
      ! local variables 
      REAL(q), ALLOCATABLE :: FTMP(:),FGTMP(:,:)
      INTEGER NMAX,I

      NMAX=SIZE(CATOM%WFCT_G,1)

      ALLOCATE(FTMP(CATOM%FR%NMAX),FGTMP(NMAX,5))

      DO I=1,CATOM%OCCUPIED_VALENCE_STATES
         FTMP(:)=CATOM%WFCT_L(:,2,I)/(CATOM%FR%R(:)) !**(ATOM_LCAO(NT)%L(I,1)+1))

         CALL BESSEL_TRANSFORM_RADIAL_FUNCTION( &
        &   CATOM%L(I,1),CATOM%FR,FTMP,REAL(SQRT(2._q*ENMAX/HSQDTM)/NMAX,KIND=q),FGTMP)

         CATOM%WFCT_G(:,:,I)=FGTMP(:,:)
      ENDDO
      DEALLOCATE(FTMP,FGTMP)

    END SUBROUTINE LCAO_WFCT_L_TO_G
 END MODULE LCAO

