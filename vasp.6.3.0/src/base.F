#include "symbol.inc"
!************************************************************************
! RCS:  $Id: base.F,v 1.2 2001/02/20 14:44:56 kresse Exp $
!
!> this module contains some control data structures for VASP
!
!***********************************************************************
      MODULE PREC
      INTEGER, PARAMETER :: q =SELECTED_REAL_KIND(10)
#ifndef noQuadPrecision
      INTEGER, PARAMETER :: qd=SELECTED_REAL_KIND(30)
#else
      INTEGER, PARAMETER :: qd=SELECTED_REAL_KIND(10)
#endif
      INTEGER, PARAMETER :: qs=SELECTED_REAL_KIND(5)
      INTEGER, PARAMETER :: qi8=SELECTED_INT_KIND(15)
!
!> this parameter controls the step width in for numerical differentiation
!> in some VASP routines
!> 1E-5 is very reliable yielding at least 7 digits in all contributions
!> to the forces
!> 1E-4, however, is better suited for second derivatives
!> for reasons of consistentcy with previous versions 1E-5
      REAL(q), PARAMETER :: fd_displacement=1E-5
      END MODULE

      MODULE BASE
      USE prec
!
!> header type information
!
        TYPE info_struct
!only  INFO
! mode information
        LOGICAL LREAL                !< real space projection/ reciprocal space proj.
        LOGICAL LOVERL               !< vanderbilt type PP read in ?
        LOGICAL LCORE                !< any partial core read in ?
        LOGICAL LCHCON               !< charge density constant during run
        LOGICAL LCHCOS               !< allways same as above in cur. impl.
        LOGICAL LONESW               !< use all bands simultaneous
        LOGICAL LONESW_AUTO          !< switch automatically between LONESW and DIIS
        LOGICAL LPRECONDH            !< precondition the subspace rotation matrix
        LOGICAL LDAVID               !< use block davidson
        LOGICAL LEXACT_DIAG          !< use exact diagonalization
        LOGICAL LRMM                 !< use RMM-DIIS algorithm
        LOGICAL LORTHO               !< orthogonalize
        LOGICAL LABORT,LSOFT         !< soft / hard stop
        LOGICAL LSTOP                !< stop with this iteration
        LOGICAL LPOTOK               !< local potential ok ?
        LOGICAL LMIX                 !< mixing done
        LOGICAL LMETAGGA             !< metagga
        LOGICAL LASPH                !< aspherical radial PAW
! exchange correlation spin
        LOGICAL LXLDA                !< calculate LDA exchange only in POTLOK
        INTEGER ISPIN                !< spin 1=no 2 =yes
        REAL(q)    RSPIN             !< 2 for spinpolar. 1 for non spinpol.
! electronic relaxation
        INTEGER ISTART               !< how to start up
        INTEGER ICHARG               !< initial charge density
        INTEGER INIWAV               !< how to initialize wavefunctions
        INTEGER INICHG               !< how to initialize charge
        INTEGER NELM                 !< maximal number of el-steps in SCC
        INTEGER NELMALL              !< maximal number of el-steps in direct optimization
        INTEGER NELMIN               !< minimal number of el-steps
        INTEGER NELMDL               !< number of delay el-steps
        INTEGER IALGO                !< algorithm for el-relax
        INTEGER IALGO_COMPAT         !< compatibility flag for IALGO
                                     !< MOD(IALGO_COMPAT,2) = 1 use old DIIS algorithm
        INTEGER NDAV                 !< number of steps in RMM-DIIS
        REAL(q)    TIME              !< timestep for el (IALGO>50)
        REAL(q)    WEIMIN,EBREAK,DEPER  !< control tags for elm
        REAL(q)    EDIFF                !< accuracy for electronic relaxation
        LOGICAL LDIAG                !< level reordering allowed or not
        LOGICAL LSUBROT              !< sub space rotation to optimize rotation matrix
        LOGICAL LPDIAG               !< sub space rotation before iterat. diag.
        LOGICAL LCDIAG               !< recalculate eigenvalues after iterat. diag.
! cutoff information
        REAL(q)    ENMAX             !< cutoff for calculations
        REAL(q)    ENINI             !< cutoff during delay
        REAL(q)    ENAUG             !< cutoff for augmentation charges
! some important things
        REAL(q)    EALLAT            !< total energy of all atoms
        REAL(q) NELECT               !< number of electrons
        REAL(q) NUP_DOWN             !< spin multiplicity
        INTEGER NBANDTOT             !< total number of bands
        INTEGER MCPU                 !< max number of proc (dimensioned)
        INTEGER NCPU                 !< actual number of proc
        LOGICAL LCORR                !< correction to forces
        INTEGER IQUASI               !< eigenvalue/occupation number corrections
        INTEGER TURBO                !< turbo mode
        INTEGER IFOLD                !< folded eigenproblem: 1=LFOLD,2=LFOLDHpsi
        INTEGER IRESTART             !< whether to restart: 2=restart with 2 optimized vectors
        INTEGER IHARMONIC            !< harmonic Ritz values
        INTEGER NREBOOT              !< number of reboots
        INTEGER NMIN                 !< reboot dimension
        REAL(q) EREF                 !< reference energy to select bands
        LOGICAL NLSPLINE             !< use spline interpolation to construct projection operator
! GW algorithms
        INTEGER ICHIREAL             !< which GW algorithm , > 0 selects space time algos
        LOGICAL LGW                  !< perform GW calculations
        LOGICAL LCHI                 !< calculate response functions
        LOGICAL LscQPGW              !< iterate wavefunctions (selfconsistent GW)
        LOGICAL LGW0                 !< leave initial wavefunctions unmodified for calculations of W
        LOGICAL LG0W0                !< leave initial wavefunctions unmodified for calculations of W
        LOGICAL LCRPA                !< swithes on CRPA in GW routines
        LOGICAL LBSE                 !< solve BSE or Cassida equations (TD-DFT)
        LOGICAL LACFDT               !< ACFDT calculation (adiabatic connection, fluctuation dissipation theorem)
        LOGICAL LOEP                 !< OEP calculation: calculate the RPA-OEP potential
        LOGICAL LEXX                 !< OEP calculation: calculate the EXX-OEP potential
        LOGICAL LHFCORRECT           !< HFCORRECT calculation (adiabatic connection, fluctuation dissipation theorem)
        INTEGER EXXOEP               !< Energy cutoff settings for OEP/LHF methods
        LOGICAL LGWNO                !< GW natural orbitals (optimized for description of Coloumb hole)
        LOGICAL LCORBSE              !< solve BSE or Cassida equations and calculate correlation energy
        LOGICAL LQPBSE               !< solve BSE or Cassida equations and calculate self-energy
        LOGICAL L2E4W                !< Calculate 2-electron 4-wannier-orbital integrals
        LOGICAL L2E4W_ALL            !< all terms not just one-center terms

! characters allways last
        CHARACTER*40 SZNAM1          !< header of INCAR
        CHARACTER*12 SZPREC          !< precision information
      END TYPE

      TYPE in_struct
!only  IO
        LOGICAL LOPEN                !< files open at startup
        INTEGER IU0                  !< unit for error
        INTEGER IU6                  !< unit for stdout
        INTEGER IU5                  !< unit for stdin
        INTEGER NWRITE               !< how much information is written out
        INTEGER IDIOT                !< how much information is written out
        INTEGER ICMPLX               !< size of a complex item upon IO
        INTEGER MRECL                !< maximal size of record length
        LOGICAL LREALD               !< no LREAL read in
        LOGICAL LMUSIC               !< jF (just a joke)
        LOGICAL LFOUND               !< WAVECAR exists ?
        LOGICAL LWAVE                !< write WAVECAR
        LOGICAL LCHARG               !< write CHGCAR
        LOGICAL LVTOT                !< write total local potential to LOCPOT
        LOGICAL LVHAR                !< write Hartree potential to LOCPOT
        LOGICAL LPDENS               !< write partial density (charge density for one band)
        INTEGER LORBIT               !< write orbit/dos
        LOGICAL LELF                 !< write elf
        LOGICAL LOPTICS              !< calculate/write optical matrix elements
        LOGICAL LPETIM               !< timing information
        INTEGER IUVTOT               !< unit for local potential
        LOGICAL INTERACTIVE          !< vasp runs interactive
        INTEGER IRECLW               !< record lenght for WAVECAR
        LOGICAL LDOWNSAMPLE          !< read WAVECAR of denser k-grid
        LOGICAL LH5                  !< general switch triggering h5 in/output
        LOGICAL LWAVEH5              !< write WAVECAR to h5
        LOGICAL LCHARGH5             !< write CHGCAR to h5
        LOGICAL LWAVEDERF            !< write file WAVEDERF
      END TYPE

      TYPE mixing
!only MIX
        INTEGER IUBROY               !< unit for broyden mixer
        REAL(q)     AMIX             !< mixing parameter A
        REAL(q)     BMIX             !< mixing parameter B
        REAL(q)     AMIX_MAG         !< mixing parameter A for magnetization
        REAL(q)     BMIX_MAG         !< mixing parameter B for magnetization
        REAL(q)     AMIN             !< minimal mixing parameter A
        REAL(q)     WC               !< weight factor for Johnsons method
        INTEGER  IMIX                !< type of mixing
        INTEGER  INIMIX              !< initial mixing matrix
        INTEGER  MIXPRE              !< form of metric for mixing
        LOGICAL  LRESET              !< reset mixer on next call (set when ions move)
        LOGICAL  HARD_RESET          !< force hard reset of mixer (force full reset regardless of MAXMIX)
        INTEGER  MAXMIX              !< maximum number of mixing steps (if positive LRESET does not apply)
        INTEGER  NEIG                !< number of eigenvalues
        INTEGER  MREMOVE             !< how many vectors are removed once
          ! iteration depth is reached
        REAL(q) :: EIGENVAL(512)     !< eigenvalues of dielectric matrix
        REAL(q) :: AMEAN             !< mean eigenvalue
        LOGICAL :: MIXFIRST          !< mix before diagonalization (or after)
      END TYPE

      TYPE symmetry
!only  SYMM
        INTEGER,POINTER:: ROTMAP(:,:,:) ! 
        REAL(q),POINTER:: TAU(:,:)      ! jF
        REAL(q),POINTER:: TAUROT(:,:)   ! jF
        REAL(q),POINTER:: WRKROT(:)     ! jF
        REAL(q),POINTER:: PTRANS(:,:)   ! jF
        REAL(q),POINTER:: MAGROT(:,:)   ! jF
        INTEGER,POINTER:: INDROT(:)     ! jF
        INTEGER :: ISYMOP(3,3,48)       !< Space group symmetry operation
        INTEGER :: IGRPOP(3,3,48)       !< Rotation part of the sapce group symmetry operation
        REAL(q) :: GTRANS(3,48)         !< Translation part of the space group symmetry operation
        INTEGER :: INVMAP(48)           !< Map to the inverse elements of each group element
        REAL(q) :: AP(3,3)              !< Lattice parameters of the primitive cell
        INTEGER :: NPCELL               !< Number of primitive cells
        INTEGER ISYM                    !< symmetry on/of
        INTEGER NROT                    !< number of rotations
        INTEGER NPTRANS                 !< number of primitive translations
      END TYPE

      TYPE prediction
!only PRED
        INTEGER IWAVPR               !< prediction of wavefunctions
        INTEGER INIPRE               !< initialized yes/no
        INTEGER IPRE                 !< what was done in wavefunction predic.
        INTEGER IUDIR                !< unit for prediction of wavefunction
        INTEGER ICMPLX               !< size of complex word
        REAL(q)  ALPHA,BETA
      END TYPE

      TYPE dipol
!only DIP
      INTEGER IDIPCO                 !< direction (0 no dipol corrections)
        LOGICAL LCOR_DIP             !< correct potential
        REAL(q) :: POSCEN(3)         !< position of center
        REAL(q) :: DIPOLC(3)         !< calculated dipol
        REAL(q) :: QUAD              !< trace of quadrupol
        INTEGER :: INDMIN(3)         !< position of minimum
        REAL(q) :: EDIPOL,EMONO,E_ION_EXTERN
      REAL(q),POINTER :: FORCE(:,:)
       REAL(q) :: VACUUM(2)          !< vacuum level
      END TYPE

      TYPE smear_struct
!only SMEAR_LOOP
        INTEGER ISMCNT               !
        REAL(q)    SMEARS(200)       !
      END TYPE


      TYPE paco_struct
!only PACO
        INTEGER NPACO                !< number of grid points for pair corr.
        REAL(q) APACO                !< cutoff
        REAL(q),ALLOCATABLE :: SIPACO(:,:) ! accumulated partial pair correlation function, second index is flattened
        INTEGER SMEANP               !< number of configurations stored
        REAL(q) OMEGA                !< average volume times number of accumulated timesteps
      END TYPE                       
                             
      TYPE energy
        REAL(q)    :: TOTENMGGA  = 0.0_q  !< total energy for METAGGA calculation
        REAL(q)    :: TOTENASPH  = 0.0_q  !< total energy for aspherical GGA
        REAL(q)    :: EBANDSTR   = 0.0_q  !< bandstructure energy
        REAL(q)    :: DENC       = 0.0_q  !< -1/2 hartree (d.c.)
        REAL(q)    :: XCENC      = 0.0_q  !< -V(xc)+E(xc) (d.c.)
        REAL(q)    :: EXCG       = 0.0_q  !< E(xc) (LDA+GGA)
        REAL(q)    :: EXCM       = 0.0_q  !< E(xc) (metaGGA)
        REAL(q)    :: EXLDA      = 0.0_q  !< LDA excchange energy
        REAL(q)    :: ECLDA      = 0.0_q  !< LDA correlation energy
        REAL(q)    :: EXGGA      = 0.0_q  !< GGA exchange energy
        REAL(q)    :: ECGGA      = 0.0_q  !< GGA correlation energy
        REAL(q)    :: EXHF       = 0.0_q  !< Hartree-Fock exchange energy
        REAL(q)    :: EXHF_ACFDT = 0.0_q  !< difference between HF energy, and exchange energy in ACFDT
        REAL(q)    :: EDOTP      = 0.0_q  !< Electric field \dot Polarization
        REAL(q)    :: TEWEN      = 0.0_q  !< Ewald energy
        REAL(q)    :: PSCENC     = 0.0_q  !< alpha Z (V(q->0) Z)
        REAL(q)    :: EENTROPY   = 0.0_q  !< Entropy term
        REAL(q)    :: PAWPS      = 0.0_q  !< paw double counting corrections
        REAL(q)    :: PAWAE      = 0.0_q  !< paw double counting corrections
        REAL(q)    :: PAWPSG     = 0.0_q  !< paw xc energies (LDA+GGA)
        REAL(q)    :: PAWAEG     = 0.0_q  !< paw xc energies (LDA+GGA)
        REAL(q)    :: PAWCORE    = 0.0_q  !< exchange correlation energy of core (LDA+GGA)
        REAL(q)    :: PAWPSM     = 0.0_q  !< paw xc energies (metaGGA)
        REAL(q)    :: PAWAEM     = 0.0_q  !< paw xc energies (metaGGA)
        REAL(q)    :: PAWCOREM   = 0.0_q  !< exchange correlation energy of core (metaGGA)
        REAL(q)    :: PAWPSAS    = 0.0_q  !< paw xc energies (aspherical)
        REAL(q)    :: PAWAEAS    = 0.0_q  !< paw xc energies (aspherical)
        COMPLEX(q) :: CVZERO     = (0.0_q, 0.0_q)  !< average local potential
        REAL(q)    :: ECGWGM     = 0.0_q  !< GW Galitskii-Migdal correlation energy
        REAL(q)    :: EKLGG0     = 0.0_q  !< Klein-contribution: Tr( G_0^-1 . G - 1 )
        REAL(q)    :: EKLLOG     = 0.0_q  !< Klein-contribution: Tr[Ln( G . G_0^-1 )]
        REAL(q)    :: ETRGHF     = 0.0_q  !< trace of HF hamiltonian with denstiy matrix
        REAL(q)    :: ELOG1G     = 0.0_q  !< Tr Log ( 1 - G S )
        REAL(q)    :: ELOGG0     = 0.0_q  !< Tr Log ( G_0 )
        REAL(q)    :: ERPA_CUT   = 0.0_q  !< RPA corr. energy with highest cutoff caculated
        REAL(q)    :: ERPA_INF   = 0.0_q  !< RPA corr. energy extrapolated to infinte basis set cutoff
        REAL(q)    :: ESCPC      = 0.0_q  !< self consistent potential correction double counting
      END TYPE

      CHARACTER (LEN=10) :: INCAR='INCAR'

      CONTAINS
!
!> small subroutine which tries to give good dimensions for 1 dimension
!
      SUBROUTINE MAKE_STRIDE (N)
      INTEGER N,NEW
      INTEGER, PARAMETER :: NGOOD=16

      NEW=(N+NGOOD)/NGOOD
      NEW=NEW*NGOOD+1
      N=NEW

      END SUBROUTINE

!****************************************************************************
! DESCRIPTION:
!>returns a double-precision real from a possible emulated quad precision
!>real  
!> @param[in]  A  (argument)
!****************************************************************************
   FUNCTION TOREAL( A ) 
      USE prec
#ifdef qd_emulate
      USE qdmodule
#endif 
      REAL(q) :: TOREAL
      QDPREAL :: A
      TOREAL = A 
   END FUNCTION TOREAL

      END MODULE
