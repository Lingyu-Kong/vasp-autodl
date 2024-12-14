#include "symbol.inc"
MODULE pseudo_struct_def
  USE prec
  USE radial_struct_def
  !
  ! pseudo potential discription include file
  ! only included if MODULES are not supported
  !only P
  INTEGER, PARAMETER:: NEKERR=100
  INTEGER, PARAMETER:: NPSPTS=1000
  INTEGER, PARAMETER:: NPSNL=100,NPSRNL=100
  !
  ! structure to support radial grid
  ! some of the entries are updated for instance by SET_PAW_ATOM_POT
  ! in this case the original date read from file are stored in entries _ORIG
  !
  TYPE potcar
     REAL(q) ZVALF,POMASS,RWIGS  ! valence, mass, wigner seitz radius
     REAL(q) ZVALF_ORIG          ! original valence
     ! ZVAL might be reset when core level shifts are calculated
     REAL(q) EATOM               ! atomic energy
     REAL(q) EATOM_CORRECTED     ! atomic energy with convergence correction
     REAL(q) EATOM_ORIG
     REAL(q) EATOM_CORRECTED_ORIG
     REAL(q) EGGA(4)             ! atomic energy for GGA
     REAL(q) ENMAXA,ENMINA,EAUG  ! energy cutoffs
     REAL(q) QOPT1,QOPT2         ! real space optimized for this cutoffs
     REAL(q) EKECUT(NEKERR),EKEERR(NEKERR) ! description of cutofferrors
     INTEGER LEXCH             ! exchange type
     REAL(q) PSGMAX              ! maximal G for local potential
     REAL(q) PSMAXN              ! maximal G for non local potential
     REAL(q) PSRMAX              ! maximal r for non local contrib.  (in fact rmax=PSRMAX/NPSNL*(NPSNL-1))
     REAL(q) PSDMAX              ! maximal r for augmentation charge (in fact rmax=PSDMAX/NPSNL*(NPSNL-1))
     REAL(q) RDEP                ! outermost grid point on radial grid used in LDA+U
                                 ! usually equivalent to largest matching radius in PSCTR
     REAL(q) RCUTRHO             ! cutoff radius for pseudo charge dens.
     REAL(q) RCUTATO             ! cutoff radius for atomic charge dens.
     REAL(q),POINTER :: RHOSPL(:,:) => NULL() ! pseudopot. charge dens. in real space spline coeff.
     REAL(q),POINTER :: ATOSPL(:,:) => NULL() ! atomic charge dens. in real space spline coeff.
     REAL(q),POINTER :: USESPL(:,:) => NULL() ! atomic charge dens. in real space spline coeff.
     REAL(q) USEZ
     REAL(q) USECUT
     REAL(q),POINTER :: PSP(:,:) => NULL() ! local pseudopotential in rec. space
     REAL(q),POINTER :: PSPCOR(:) => NULL()! partial core information in rec. space
     REAL(q),POINTER :: PSPTAU(:) => NULL()! partial kinetic energy density  information in rec. space
     REAL(q),POINTER :: PSPTAUVAL(:) => NULL()! kinetic energy density of valence electrons information in rec. space
     REAL(q),POINTER :: PSPRHO(:) => NULL()! atomic pseudo charge density in rec. space
     REAL(q),POINTER :: PSPNL(:,:) => NULL()  ! non local proj. rec. space
     REAL(q),POINTER :: PSPNL_SPLINE(:,:,:) => NULL() ! non local proj. rec. space spline fit
     REAL(q),POINTER :: PSPRNL(:,:,:) => NULL() ! non local proj. real space
     REAL(q),POINTER :: PSPRNL_ORIG(:,:,:) => NULL()
     REAL(q),POINTER :: DION(:,:) => NULL()! non local strength
     REAL(q),POINTER :: DION_ORIG(:,:) => NULL()
     REAL(q),POINTER :: QION(:,:) => NULL()! spherical augmentation charge
     REAL(q),POINTER :: QTOT(:,:) => NULL()! total charge in each channel
     REAL(q),POINTER :: QPAW(:,:,:) => NULL()! integrated augmentation charge 
             ! stores essentially WAE*WAE-WPS*WPS r^l
     COMPLEX(q),POINTER :: JPAW(:,:,:,:) => NULL()! 
             ! stores essentially <WAE| j(r) Y_lm(r) | WAE> -<WPS| j(r) Y_lm(r) | WPS>
     REAL(q),POINTER :: QPAW_FOCK(:,:,:,:) => NULL()
             ! similar to QPAW but stores coefficients that determine how much of 
             ! AUG_FOCK is added
     REAL(q),POINTER :: QATO(:,:) => NULL()! initial occupancies (in atom)
     REAL(q),POINTER :: QDEP(:,:,:) => NULL() ! L-dependent augmentation charges on regular grid
     REAL(q),POINTER :: QDEP_FOCK(:,:,:) => NULL() ! L-dependent charge with zero moment on regular grid
                                 ! for technical reasons this also includes all elements in QDEP
                                 ! see comments in fast_aug.F
     REAL(q),POINTER :: NABLA(:,:,:) => NULL()! atomic augmentation matrix elements of nabla operator
     REAL(q) PSCORE              ! is equal to V(q) + 4 Pi / q^2
     REAL(q) ESELF               ! self energy of pseudized ion (usually not used)
     INTEGER,POINTER :: LPS(:) => NULL()   ! L quantum number for each PP
     INTEGER,POINTER :: NLPRO(:) => NULL() ! unused
     INTEGER,POINTER :: NDEP(:,:) => NULL()! number of augmentation channels per ll'
     REAL(q),POINTER :: E(:) => NULL()     ! linearization energies
     REAL(q),POINTER :: E_ORIG(:) => NULL()! linearization energies as stored on POTCAR
     ! quantities defined on radial grid
     TYPE (rgrid)    :: R        ! radial grid
     REAL(q),POINTER :: POTAE(:) => NULL() ! frozen core AE potential on r-grid for atomic reference configuration (valence only) as read from file (must not be updated to get correct core states)
     REAL(q),POINTER :: POTAE_ORIG(:) => NULL()
     REAL(q),POINTER :: POTAE_XCUPDATED(:) => NULL()! as above, might be updated if xc potential changes
     REAL(q),POINTER :: POTPS(:) => NULL() ! local PP on r-grid for atomic reference configuration (valence only)
     REAL(q),POINTER :: POTPS_ORIG(:) => NULL()
     REAL(q),POINTER :: POTPSC(:) => NULL()! local PP on r-grid (core only), V_H[\tilde{n}_Zc]
     REAL(q),POINTER :: POTPSC_ORIG(:) => NULL()
     REAL(q),POINTER :: KLIC(:) => NULL()  ! averaged local exchange potential (KLI) on r-grid
     REAL(q),POINTER :: RHOAE(:) => NULL() ! frozen core charge rho(r)r^2 on r-grid
     REAL(q),POINTER :: TAUAE(:) => NULL() ! kinetic energy density of core
     REAL(q),POINTER :: RHOPS(:) => NULL() ! frozen pseudo partial core charge rho(r)r^2 on r-grid
     REAL(q),POINTER :: TAUPS(:) => NULL() ! kinetic energy density of partial core
     REAL(q),POINTER :: WAE(:,:) => NULL() ! ae valence wavefunction on r-grid
     REAL(q),POINTER :: WPS(:,:) => NULL() ! pseudo valence wavefunction on r-grid
     REAL(q),POINTER :: AUG(:,:) => NULL() ! L-dependent augmentation charge on r-grid
     REAL(q),POINTER :: AUG_SOFT(:,:) => NULL()  ! L-dependent augmentation charge on r-grid (soft version)
     REAL(q),POINTER :: AUG_FOCK(:,:,:) => NULL()! L-dependent charge with zero moment on r-grid
     REAL(q)         :: DEXCCORE ! exchange correlation energy of frozen core
     REAL(q)         :: DEXCCORE_ORIG ! exchange correlation energy of frozen core
     REAL(q)         :: DEXCCOREM! exchange correlation energy of frozen core (MetaGGA)
     ! relaxed core stuff
     REAL(q),POINTER :: C(:,:) => NULL()     ! partial wave expansion coefficients
     REAL(q),POINTER :: BETA(:,:) => NULL()  ! projectors on radial grid
     REAL(q),POINTER :: RHOAE00(:) => NULL() ! spherical component of valence AE charge density
     REAL(q),POINTER :: RHOPS00(:) => NULL() ! spherical component of valence PS charge density
     REAL(q),POINTER :: RHOPSPW(:) => NULL() ! spherical component of valence PS charge density from PW grid
     REAL(q),POINTER :: V00_AE(:) => NULL()  ! spherical component of AE potential
     REAL(q),POINTER :: V00_PS(:) => NULL()  ! spherical component of PS potential
     REAL(q),POINTER :: WKINAE(:,:) => NULL()! T | \psi_i >
     REAL(q),POINTER :: WKINPS(:,:) => NULL()! T | \tilde{\psi}_i >
     REAL(q),POINTER :: DIJ(:,:) => NULL()   ! < \psi_i | T+V |\psi_j > - < \tilde{\psi}_i | T+\tilde{V} |\tilde{\psi}_j >
     REAL(q),POINTER :: QIJ(:,:) => NULL()   ! < \psi_i | \psi_j > - < \tilde{\psi}_i | \tilde{\psi}_j >
     REAL(q),POINTER :: POTAEC(:) => NULL()  ! AE core potential 
     REAL(q)         :: AVERAGEPOT(3) ! average local potential (AE,PS,PW)
     REAL(q)         :: VPSRMAX    ! local PP at boundary
     REAL(q),POINTER :: CLEV(:) => NULL()    ! core state eigenenergies
     REAL(q)         :: ECORE(2)   ! core contributions to total energy
     REAL(q)         :: VCA        ! weight of this potential (can be overwritten by INCAR VCA tag)
     ! atomic stuff
     REAL(q),POINTER :: ATOMIC_J(:) => NULL()   ! relativistic quantum number
     REAL(q),POINTER :: ATOMIC_E(:) => NULL()   ! eigenenergy
     REAL(q),POINTER :: ATOMIC_OCC(:) => NULL() ! occupation number
     INTEGER,POINTER :: ATOMIC_N(:) => NULL()   ! main quantum number
     INTEGER,POINTER :: ATOMIC_L(:) => NULL()   ! angular momentum
     ! most integers are on the end (alignment ?)
     INTEGER      LDIM           ! (leading) dimension for l channels (>=LMAX)
     INTEGER      LMAX           ! total number of l-channels for non local PP
     INTEGER      LMDIM          ! (leading) dimension for lm channels (>=LMMAX)
     INTEGER      LMMAX          ! total number nlm-channels for non local PP
     INTEGER      LDIM2          ! dimension for augmentation arryas
     LOGICAL      LREAL          ! real space optimized ?
     LOGICAL      LUNSCR         ! partial core has been unscreened
     INTEGER      LMAX_CALC      ! maximum L for onsite terms in PAW
     ! might be overwritten by LMAXPAW line in the INCAR file
     REAL(q)      ZCORE          ! charge of core
     CHARACTER*40 SZNAMP         ! header
     CHARACTER*2  ELEMENT        ! Name of element
     REAL(q), POINTER :: OPTPROJ(:,:,:) => NULL()   ! optimized projectors for each quantum number L
  END TYPE potcar
END MODULE pseudo_struct_def
