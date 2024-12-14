#include "symbol.inc"
MODULE wave_struct_def
  USE prec
  USE mpimy
  USE mgrid_struct_def

  !
  !  structure required to support storage of wavefunction for
  !  several kpoints and spins
  !
  TYPE wavedes
     REAL(q) RSPIN                 ! spin multiplicity
     REAL(q) ENMAX                 ! energy cutoff
     INTEGER NRSPINORS             ! number of spinors (1 for collinear, 2 for non collinear)
     INTEGER NGDIM                 ! first dimension of any array related to the plane wave basis
     INTEGER NRPLWV                ! first dimension of wavefunction array
     ! collinear:  NRPLWV=NGDIM, noncollinear:  NRPLWV=2*NGDIM
     INTEGER NRPLWV_RED            ! local number of coefficients in wave function array after data redistribution
     INTEGER NPROD                 ! first dimension of projected wave array
     INTEGER NPRO                  ! local number of elements in projected wave array
     INTEGER NPRO_TOT              ! total number of elements (summed over all nodes)
     ! NPRO, NPROD, and NPRO_TOT are all doubled in the non collinear version
     INTEGER NPROD_RED             ! dimension of projected wave array after redistribution
     INTEGER NBANDS                ! local number of bands
     INTEGER NB_TOT                ! total number bands
     INTEGER NB_PAR                ! distribution over bands (number of bands done in parallel )= WDES%COMM_INTER%NCPU
     INTEGER NSIM                  ! band blocking (mainly for seriel version)
     INTEGER NB_LOW                ! lowest band index in global
     INTEGER NKDIM                 ! total number of k-points in the entire Brillouin zone (BZ)
     ! required for HF calculations (otherwise equal to NKPTS)
     INTEGER NKPTS_FOR_GEN_LAYOUT  ! number of k-points used for the generation of the data layout
     ! this must not change when the number of k-point changes
     INTEGER NKPTS                 ! number of k-points in the irreducable wedge of the BZ (IBZ)
     INTEGER ISPIN                 ! number of spins
     INTEGER NCDIJ                 ! dimension of arrays like CDIJ, CQIJ
     INTEGER NIONS                 ! number of ions stored locally 
     INTEGER NTYP                  ! number of types stored locally
     TYPE (grid_3d), POINTER ::GRID! pointer to a grid if FFT's are required
     INTEGER,POINTER :: NPLWKP(:)  ! number of coefficients for each k-point and band per node
     INTEGER,POINTER :: NGVECTOR(:)! number of G-vectors in the basis for each k-point per node
     ! collinear: NPLWKP= NGVECTOR, noncollinear NPLWKP = 2*NGVECTOR 
     ! NGVECTOR is the same for collinear and non collinear calculations
     ! (summed over nodes, doubled in the non collinear case)
     INTEGER,POINTER :: NGVECTOR_POS(:)! sum of NGVECTOR up to (but not including) the current node
     INTEGER,POINTER :: NPLWKP_TOT(:)  ! total number of coefficients in plane wave array at each k-points
     INTEGER,POINTER :: NB_TOTK(:,:)! number of bands to be calculated for each k-point and spin
     ! possibly smaller than NB_TOT
     INTEGER         :: NCOL       ! number of columns
     INTEGER,POINTER,CONTIGUOUS :: PL_INDEX(:,:) ! index a column would have in serial version
     INTEGER,POINTER,CONTIGUOUS :: PL_COL(:,:)! number of plane wave in this column
     INTEGER,POINTER ::NPRO_POS(:) ! for each atom, start index of entries in CPROJ in serial version
     INTEGER,POINTER :: LMMAX(:)   ! total number of NLM quantum numbers for each type
     INTEGER,POINTER :: LMBASE(:)  !
     INTEGER,POINTER :: NITYP(:)   ! number of ions stored locally for each type
     INTEGER,POINTER :: ITYP(:)    ! type for each ion
     INTEGER,POINTER ::NT_GLOBAL(:)! global type index for this type 
     REAL(q),POINTER :: VKPT(:,:)  ! coordinate of k-point
     REAL(q),POINTER :: WTKPT(:)   ! symmetry weight-factor for each k-point
     INTEGER,POINTER,CONTIGUOUS :: NINDPW(:,:)! index to the FFT box for each pw comp and k-point
     LOGICAL,POINTER,CONTIGUOUS :: LUSEINV(:) ! for each k-point decides whether reduced G grid can be used (compare AT_GAMMA)
     INTEGER,POINTER,CONTIGUOUS :: NINDPW_INV(:,:)! index to the FFT box for each pw comp and k-point to the G vector -G-k
     REAL(q),POINTER,CONTIGUOUS :: FFTSCA(:,:,:)  ! scaling index if plane wave coefficients are reduced (LUSEINV .TRUE.)
     INTEGER,POINTER,CONTIGUOUS :: MAP_TO_FULL(:,:,:) ! map from half-grid mode to full grid mode
     INTEGER,POINTER,CONTIGUOUS :: IGX(:,:)   ! x index of each pw comp and k-point
     INTEGER,POINTER,CONTIGUOUS :: IGY(:,:)   ! y index of each pw comp and k-point
     INTEGER,POINTER,CONTIGUOUS :: IGZ(:,:)   ! z index of each pw comp and k-point
     REAL(q),POINTER,CONTIGUOUS :: DATAKE(:,:,:) ! kinetic energy for each plane wave
     ! last index labels up and down components
     ! of the spinor in case of spin spirals
     REAL(q) QSPIRAL(3)            ! propagation vector of spin spiral
     TYPE(communic),POINTER  :: COMM,COMM_INTER,COMM_INB
     TYPE(communic),POINTER  :: COMM_KINTER,COMM_KIN
     TYPE(communic),POINTER  :: COMM_SHMEM,COMM_intra_node,COMM_inter_node
     REAL(q) SAXIS(3)              ! quantisation axis of the spin operator
! TODO maybe AT_GAMMA should be removed
     LOGICAL,POINTER :: AT_GAMMA(:)! indicates that a k-point corresponds to gamma
                                   ! selects special treatment
     LOGICAL LORBITALREAL          ! special treatment at gamma
     LOGICAL LOVERL                ! overlap required
     LOGICAL DO_REDIS              ! data redistribution required
     LOGICAL LNONCOLLINEAR         ! noncollinear calculations
     LOGICAL LSORBIT               ! spin orbit coupling
     LOGICAL LGAMMA                ! gamma point only, projected wavefunction character is REAL
                                   ! this is only .TRUE. if precompiler flag gammareal is define 
     LOGICAL LSPIRAL               ! calculate spin spirals?
     LOGICAL LZEROZ                ! set m_z to zero in SET_CHARGE?
     INTEGER NBANDSLOW             ! lowest band to be optimized (-1 no restrictions)
     INTEGER NBANDSHIGH            ! highest band to be optimized (-1 no restrictions)
  END TYPE wavedes

  !
  ! description for one k point
  ! contains also all information required for simple calculations in real space
  !
  TYPE wavedes1
     !only WDES1
     REAL(q) RSPIN                 ! spin multiplicity
     INTEGER NRSPINORS             ! number of spinors (1 for collinear, 2 for non collinear)
     INTEGER NGDIM                 ! first dimension of any array related to the basis
     INTEGER NRPLWV                ! first dimension of wavefunction array (stores coefficients in rec. space)
     ! collinear:  NRPLWV=NGDIM, noncollinear:  NRPLWV=2*NGDIM
     INTEGER NRPLWV_RED            ! local number of electrons in wavefunctionarray after data redistribution
     INTEGER NPROD                 ! first dimension of projected wave array
     INTEGER NPRO                  ! local number of elements in projected wave array
     INTEGER NPRO_O_RED            ! local number of elements for overlap (0 for normconserving pot) after resdistribution
     INTEGER NPRO_TOT              ! total number of elements
     INTEGER NPROD_RED             ! local number of elements in projected wave array after redistribution
     ! NPRO, NPROD, and NPRO_TOT are all doubled in the non collinear version
     INTEGER NPRO_RED              ! local number of elements in projected wave array after redistribution
     INTEGER NBANDS                ! bands
     INTEGER NB_TOT                ! total number bands
     INTEGER NB_PAR                ! distribution over bands (number of bands done in parallel )= WDES%COMM_INTER%NCPU
     INTEGER NSIM                  ! band blocking (mainly for serial version)
     INTEGER NB_LOW                ! lowest band index in global 
     INTEGER NPL                   ! number of plane waves coefficients (local)
     INTEGER NPL_RED               ! number of plane waves coefficients after data redistribution
     INTEGER NGVECTOR              ! number of G-vectors in the basis (local)
     ! collinear: NGVECTOR == NPL, noncollinear 2*NGVECTOR == NPL
     INTEGER NGVECTOR_POS          ! sum of NGVECTOR up to (but not including) the current node
     INTEGER NPL_TOT               ! total number of plane waves (global)
     INTEGER NB_TOTK(2)            ! number of bands to be calculated
     INTEGER NIONS                 ! number of ions stored locally
     INTEGER NTYP                  ! number of types stored locally
     REAL(q) RINPL                 ! inverse of total number of plane waves
     TYPE (grid_3d), POINTER ::GRID=>NULL() ! pointer to a grid if FFT's are required
     INTEGER NK                    ! k-point number (required for non-local proj.)
     ! few things which are only required in parallel version
     INTEGER         :: NCOL       ! number of columns
     INTEGER,POINTER,CONTIGUOUS :: PL_INDEX(:)=>NULL() ! index a column would have in serial version
     INTEGER,POINTER,CONTIGUOUS :: PL_COL(:)=>NULL()   ! number of plane waves in this column
     INTEGER,POINTER :: NPRO_POS(:)=>NULL() ! index CPROJ would have in serial version
     INTEGER,POINTER :: LMMAX(:)=>NULL()    ! total number of NLM quantum numbers for each type
     INTEGER,POINTER :: LMBASE(:)=>NULL()   !
     INTEGER,POINTER :: NITYP(:)=>NULL()    ! number of ions stored locally for each type
     INTEGER,POINTER :: ITYP(:)=>NULL()     ! type for each ion
     INTEGER,POINTER ::NT_GLOBAL(:)=>NULL() ! global type index for this type 
     INTEGER,POINTER,CONTIGUOUS :: NINDPW(:)=>NULL()   ! index to the FFT box for each pw comp
     LOGICAL                    :: LUSEINV    ! for each k-point decides whether reduced G grid can be used (compare AT_GAMMA)
     INTEGER,POINTER,CONTIGUOUS :: NINDPW_INV(:)=>NULL() ! index to the FFT box for each pw comp and k-point to the G vector -G-k
     REAL(q),POINTER,CONTIGUOUS :: FFTSCA(:,:)=>NULL()   ! scaling index if plane wave coefficients are reduced (LUSEINV .TRUE.)
     INTEGER,POINTER,CONTIGUOUS :: MAP_TO_FULL(:,:)=>NULL()  ! map from half-grid mode to full grid mode
     INTEGER,POINTER,CONTIGUOUS :: IGX(:)=>NULL()      ! x index of each pw comp and k-point
     INTEGER,POINTER,CONTIGUOUS :: IGY(:)=>NULL()      ! y index of each pw comp and k-point
     INTEGER,POINTER,CONTIGUOUS :: IGZ(:)=>NULL()      ! z index of each pw comp and k-point
     REAL(q),POINTER :: VKPT(:)=>NULL()     ! k-point
     REAL(q)         :: WTKPT      ! symmetry weight-factor for this k-point
     REAL(q),POINTER,CONTIGUOUS :: DATAKE(:,:)=>NULL() ! kinetic energy for each plane wave
     ! last index labels up and down components
     ! of the spinor in case of spin spirals
     REAL(q) QSPIRAL(3)            ! propagation vector of spin spiral
     TYPE(communic),POINTER  :: COMM=>NULL() 
     TYPE(communic),POINTER  :: COMM_INTER=>NULL() ,COMM_INB=>NULL() 
     TYPE(communic),POINTER  :: COMM_KINTER=>NULL() ,COMM_KIN=>NULL() 
     TYPE(communic),POINTER  :: COMM_SHMEM=>NULL() ,COMM_intra_node=>NULL() ,COMM_inter_node=>NULL() 
     REAL(q) SAXIS(3)              ! quantisation axis of the spin operator
     LOGICAL AT_GAMMA              ! indicates that a k-point corresponds to gamma
                                   ! selects special treatment
     LOGICAL LORBITALREAL          ! special treatment at gamma
     LOGICAL LOVERL                ! overlap required
     LOGICAL DO_REDIS              ! data redistribution required
     LOGICAL LNONCOLLINEAR         ! allows one to turn on noncollinear calculations
     LOGICAL LSORBIT               ! spin orbit coupling
     LOGICAL LGAMMA                ! gamma point only, projected wavefunction character is REAL
     LOGICAL LSPIRAL               ! calculate spin spirals?
     LOGICAL LZEROZ                ! set m_z to zero in SET_CHARGE?
  END TYPE wavedes1

  !
  ! structure required to store a set of orbitals including band index, k-point and spin
  ! in most places in VASP structures of this type are defined as TYPE(wavespin) :: W
  !
  TYPE wavespin
     TYPE(wavedes),POINTER:: WDES         ! descriptor for orbital
     COMPLEX(q),POINTER           :: CPTWFP(:,:,:,:)
     GDEF      ,POINTER           :: CPROJ(:,:,:,:)
     REAL(q),   POINTER:: FERWE(:,:,:)    ! fermi-weight for each band stored on local node
                                   ! local band index, k-point, spin
     REAL(q),   POINTER:: AUX  (:,:,:)    ! auxilary
                                   ! local band index, k-point, spin
     COMPLEX(q),POINTER:: CELEN(:,:,:)    ! eigenvalues
                                   ! local band index, k-point, spin
     REAL(q),   POINTER:: FERTOT(:,:,:)   ! global array for fermi-weights
                                   ! global band index, k-point, spin
     REAL(q),   POINTER:: AUXTOT(:,:,:)   ! global array for auxilary quantities
                                   ! global band index, k-point, spin
     COMPLEX(q),POINTER:: CELTOT(:,:,:)   ! global array for eigenvalues
                                   ! global band index, k-point, spin
     LOGICAL           :: OVER_BAND       ! distribution over bands or not
     ! store fermi energy
     REAL(q), POINTER :: EFERMI(:) ! fermi energy for each spin
  END TYPE wavespin

  !
  ! structure required to store a set of orbitals including and band index, k-point
  ! this type is hardly used in VASP, but was required during early porting to F90
  !
  TYPE wavefun
     TYPE(wavedes),POINTER:: WDES       ! descriptor for orbital
     COMPLEX(q),POINTER,CONTIGUOUS:: CPTWFP(:,:,:) ! plane wave coefficients of orbital (orbital) C_G,n,k
                                   ! indices: plane wave, local band index, k-point
     GDEF      ,POINTER,CONTIGUOUS:: CPROJ(:,:,:)  ! orbital character <p_alpha| phi_n,k>
                                   ! indices: projectors, local band index, k-point
     COMPLEX(q),POINTER,CONTIGUOUS:: CR(:,:,:)     ! orbital in real space
                                   ! hardly ever allocated
     REAL(q),   POINTER:: FERWE(:,:)    ! fermi-weight for each band
                                   ! local band index, k-point
     REAL(q),   POINTER:: AUX  (:,:)    ! auxilary
                                   ! local band index, k-point
     COMPLEX(q),POINTER:: CELEN(:,:)    ! eigenvalues
                                   ! local band index, k-point
     REAL(q),   POINTER:: FERTOT(:,:)   ! global array for fermi-weights
                                   ! global band index, k-point
     REAL(q),   POINTER:: AUXTOT(:,:)   ! global array for auxilary quantities
                                   ! global band index, k-point
     COMPLEX(q),POINTER:: CELTOT(:,:)   ! global array for eigenvalues
                                   ! global band index, k-point
     LOGICAL, POINTER  :: LOPT(:,:,:)   ! optimize this orbital
     LOGICAL           :: OVER_BAND     ! distribution over bands or not
  END TYPE wavefun

  !
  ! structure required to store one orbital (no spin no bands)
  !
  TYPE wavefun1
     TYPE(wavedes1), POINTER         :: WDES1     => NULL()  ! descriptor for orbital
     COMPLEX(q), POINTER, CONTIGUOUS :: CPTWFP(:) => NULL()  ! plane wave coefficients of a single orbital
     GDEF      , POINTER, CONTIGUOUS :: CPROJ(:)  => NULL()  ! orbital character <p_alpha| phi>
     COMPLEX(q), POINTER, CONTIGUOUS :: CR(:)     => NULL()  ! orbital in real space
     REAL(q)           :: FERWE       ! fermi-weight for each band
     REAL(q)           :: AUX         ! auxilary
     COMPLEX(q)        :: CELEN       ! eigenvalues
     INTEGER           :: NB          ! band index if it applies (not used)
     INTEGER           :: ISP         ! spin index if it applies (not used)
     LOGICAL           :: LDO         ! initialised and operational
  END TYPE wavefun1

  !
  ! structure required to store a set of orbitals (no spin and k dependency)
  !
  TYPE wavefuna
     TYPE(wavedes1), POINTER         :: WDES1          => NULL() ! descriptor for orbital
     COMPLEX(q), POINTER, CONTIGUOUS :: CPTWFP(:,:)    => NULL() ! plane wave coefficients of a set of orbitals
     GDEF      , POINTER, CONTIGUOUS :: CPROJ(:,:)     => NULL() ! orbital character <p_alpha| phi>
     COMPLEX(q), POINTER, CONTIGUOUS :: CW_RED(:,:)    => NULL() ! redistributed orbitals
     GDEF      , POINTER, CONTIGUOUS :: CPROJ_RED(:,:) => NULL() ! orbital character redistributed
     COMPLEX(q), POINTER, CONTIGUOUS :: CR(:,:)        => NULL() ! orbital in real space
     REAL(q),    POINTER :: FERWE(:)  ! fermi-weight for each band
     REAL(q),    POINTER :: AUX  (:)  ! auxilary
     COMPLEX(q), POINTER :: CELEN(:)  ! eigenvalues
     INTEGER             :: ISP       ! spin index if apply (not used)
     INTEGER             :: FIRST_DIM ! first dimension, if acessed as a multidimensional array
     LOGICAL             :: OVER_BAND ! distribution over bands or not (not used)
  END TYPE wavefuna

END MODULE wave_struct_def
