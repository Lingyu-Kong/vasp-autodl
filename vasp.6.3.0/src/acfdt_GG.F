#include "symbol.inc"
! use this flag for old (automatic) NOMEGA distribution

MODULE ACFDT_GG
  USE chi_base
  USE fock 
  USE acfdt
  IMPLICIT NONE

  INTEGER, PRIVATE :: NE=8
  REAL(q), PARAMETER, PRIVATE :: STEP_BETWEEN_FREQUENCIES=1.05
! alternative spanning the same frequency range but with only 4 points
!  INTEGER, PARAMETER, PRIVATE :: NE=4
!  REAL(q) :: STEP_BETWEEN_FREQUENCIES=1.120576983388_q

  REAL(q), PRIVATE :: DATAKEMAX                   ! maximum energy used in kernel
! Fermi-wave vector for HEG kernel
! set some default to avoid crashing 
!
  REAL(q), PRIVATE :: KFERMI=2.0_q

  !BLACS related variables and function
  !PROCS.. number of processors, ME... processor number, NPROW... number of rows in process grid
  !NPCOL... number of columns in process grid, myrow,mycol... my coordinates in process grid
  !MB... blocking factor of rows for block cyclic distribution
  !MB... blocking factor of columns for block cyclic distribution      
  INTEGER, SAVE :: PROCS,ME
  ! data distribution for block cyclic case
  INTEGER, SAVE :: NPROW,NPCOL,MYROW,MYCOL

  ! data distribution for column wise data distribution
  INTEGER, SAVE :: NPROW0,NPCOL0,MYROW0,MYCOL0

  ! blocking factors
  INTEGER, SAVE :: MB=64, NB=64
  INTEGER, SAVE :: DESCA(9),DESCB(9)
  INTEGER, EXTERNAL :: NUMROC, INDXG2P
  !numroc... function that returns the number of rows or columns of the LOCAL ARRAY if you enter
  !the corresponding properties of the GLOBAL ARRAY
#ifdef scaLAPACK

CONTAINS

!********************************* XI_ACFDT_SETUP   *************************************
!
! setup the energy cutoffs at which the correlation energy is 
! calculated 
! 
!****************************************************************************************

  SUBROUTINE XI_ACFDT_SETUP_GG( COR, ENCUTGW, ENCUTGWSOFT)
    USE constant
    TYPE (correlation), POINTER :: COR
    REAL(q) :: ENCUTGW, ENCUTGWSOFT
    INTEGER :: I

    ! use single shot method with SCK for RPA 
    IF ( LSCK ) THEN
       NE = 1 
    ELSE
       NE = 8 
    ENDIF 

    ALLOCATE(COR)
    COR%NE=NE
    ALLOCATE(COR%ENCUTGW(NE))
    ALLOCATE(COR%ENCUTGWSOFT(NE))
    ALLOCATE(COR%CORRELATION(NE))
    ALLOCATE(COR%CORRSOSEX(NE))
    ALLOCATE(COR%CORRELATION_K(NE))
    ALLOCATE(COR%CORRSOSEX_K(NE))
    ALLOCATE(COR%CORRMP2DIR(NE))
    ALLOCATE(COR%CORRMP2EX(NE))
    ALLOCATE(COR%CORRMP2DIR_K(NE))
    ALLOCATE(COR%CORRMP2EX_K(NE))
    ALLOCATE(COR%CORRELATION_LAMBDA(NE,0:NLAMBDA))

    COR%ENCUTGW(1)    =ENCUTGW
    COR%ENCUTGWSOFT(1)=ENCUTGWSOFT
    DO I=2,NE
       COR%ENCUTGW(I)    =COR%ENCUTGW(I-1)/STEP_BETWEEN_FREQUENCIES
       COR%ENCUTGWSOFT(I)=COR%ENCUTGWSOFT(I-1)/STEP_BETWEEN_FREQUENCIES
       IF (COR%ENCUTGWSOFT(I)<=0) COR%ENCUTGWSOFT(I)=-1
    ENDDO

    COR%CORRELATION_LAMBDA=0
    COR%CORRELATION=0
    COR%CORRSOSEX=0
    COR%CORRELATION_K=0
    COR%CORRSOSEX_K=0
    COR%CORRMP2DIR=0
    COR%CORRMP2EX=0
    COR%CORRMP2DIR_K=0
    COR%CORRMP2EX_K=0
    
  END SUBROUTINE XI_ACFDT_SETUP_GG


  SUBROUTINE XI_ACFDT_DEALLOCATE_GG( COR )
    TYPE (correlation), POINTER :: COR

    DEALLOCATE(COR%ENCUTGW)
    DEALLOCATE(COR%CORRELATION)
    DEALLOCATE(COR%CORRSOSEX)
    DEALLOCATE(COR%CORRELATION_K)
    DEALLOCATE(COR%CORRSOSEX_K)
    DEALLOCATE(COR%CORRMP2DIR)
    DEALLOCATE(COR%CORRMP2EX)
    DEALLOCATE(COR%CORRMP2DIR_K)
    DEALLOCATE(COR%CORRMP2EX_K)
    DEALLOCATE(COR%CORRELATION_LAMBDA)
    DEALLOCATE(COR)
    NULLIFY(COR)

  END SUBROUTINE XI_ACFDT_DEALLOCATE_GG

!********************************* CALCULATE_RPA_CORRELATION_ENERGY *********************
!
! determine the RPA correlation energy
! this subroutine integrates the entire calculation in one
! routine to conserve memory and speed things up
! 
!****************************************************************************************
SUBROUTINE CALCULATE_RPA_CORRELATION_ENERGY( WDES, CHI, NROWS, NCOLS, WGWQ, &
    IMAG_GRIDS, LATT_CUR, COR, NQ , IDIR_MAX, IO)
    USE base
    USE constant
    USE mpimy
    USE ini
    USE minimax_struct, ONLY : imag_grid_handle
    USE minimax, ONLY : LOCAL_INDEX_TO_GLOBAL 
    IMPLICIT NONE
    TYPE (wavedes)          :: WDES
    TYPE (responsefunction) :: CHI                  !response function
    INTEGER                 :: NROWS, NCOLS         !number of local rows and cols of CHI
    TYPE (wavedes1)         :: WGWQ
    TYPE (imag_grid_handle) :: IMAG_GRIDS           !imaginary grids 
    TYPE (latt)             :: LATT_CUR
    TYPE (correlation)      :: COR                  !correlation energy
    INTEGER                 :: NQ                   !current q-point
    INTEGER                 :: IDIR_MAX             !maximum number of directions
    TYPE (in_struct)        :: IO                 
  ! local
    GDEF, POINTER, DIMENSION(:,:) :: CHI_WORK

    TYPE (communic)         :: COMM_BETWEENOMEGA      !communicator between tau
    TYPE (communic)         :: COMM_INOMEGA           !communicator in tau
    INTEGER                 :: NOMEGA               !number of frequency points
    INTEGER                 :: NOMEGA_SIMULTANEOUS  !# of freq. points on one core

    INTEGER                 :: NROWS_NEW, NCOLS_NEW !new size of redistributed chi 
    INTEGER                 :: NBLOCK                 !block size
    INTEGER                 :: NODE_CRITICAL          !node containing cut off
    INTEGER ::  IDIR, NP_NEW, I,IP, ILAMBDA, NK, NCURR
    INTEGER :: I_, J
    COMPLEX(q) :: SUM, SUMMP2, SUMSOSEX, SUMMP2EX, SUM_sr, SUMMP2_sr
    REAL(q) :: LAMBDA
    INTEGER :: NOMEGA_GLOBAL
    INTEGER :: NOMEGA_IN_ROOT_GROUP                   ! # of frequencies in root group
    INTEGER :: NOMEGA_CURRENT
    LOGICAL :: LDMODE
    INTEGER, SAVE :: NK_LAST=-1

    ! use single shot method with SCK for RPA 
    IF ( LSCK ) THEN
       NE = 1 
    ELSE
       NE = 8 
    ENDIF 

    ! total number of points
    NOMEGA = IMAG_GRIDS%NOMEGA

    NOMEGA_SIMULTANEOUS = IMAG_GRIDS%B%NPOINTS_IN_GROUP
    COMM_INOMEGA = IMAG_GRIDS%B%COMM_IN_GROUP
    COMM_BETWEENOMEGA = IMAG_GRIDS%B%COMM_BETWEEN_GROUPS

    IF (CHI%LREAL .AND. .NOT. CHI%LREALSTORE) THEN
       CALL vtutor%bug("internal error in CALCULATE_RPA_CORRELATION_ENERGY: \n or the Gamma point " &
          // "version CHI%LREALSTORE must be set", __FILE__, __LINE__)
    ENDIF

    PROFILING_START('calculate_rpa_correlation_energy')

    !each group needs to loop over the same # of frequencies internally,
    !otherwise when initializing BLACS in LNTRACE_OF_CHI the routine 
    !BLACS_PINFO in GEN_PROC_GROUP_GRIDS_BC waits for the CPUs in other group
    !to call itself. 
    !since the root group contains always the maximum # of frequencies we 
    !use NOMEGA_SIMULTANEOUS from the root node for the internal # of frequencies for all groups
    NOMEGA_IN_ROOT_GROUP = 0
    IF ( WDES%COMM_KIN%NODE_ME == 1 ) NOMEGA_IN_ROOT_GROUP = NOMEGA_SIMULTANEOUS
    CALLMPI( M_bcast_i( WDES%COMM_KIN, NOMEGA_IN_ROOT_GROUP, 1 ) ) 

    !scaLAPACK block size and processor grid is set here
    IF ( CHI%NQ == 1 .AND. NK_LAST < 0 ) THEN
       CALL BLOCKSIZE_AND_PROC_GRID( NROWS, NCOLS, COMM_INOMEGA%NCPU, IO%IU0)
    ENDIF 
    NK_LAST = NK_LAST +1
 
    ! Merzuk's version changed the sign of RESPONSEFUN here
    ! I moved this over to the place where the RESPONSEFUN is first used

    !--------------------------------------------------------------------------------------------
    lamb: DO ILAMBDA=0,NLAMBDA                       !optionally loop over lambda parameter
    LAMBDA=1.0_q*(ILAMBDA+1)/MAX(1,NLAMBDA)          !set current lambda value
       !-----------------------------------------------------------------------------------------
       cutoff: DO I=1, NE                          !loop over different energy cutoffs
       !-----------------------------------------------------------------------------------------
          !initialize correlation energy for the current energy cutoff # I:
          COR%CORRELATION_K(I)=0
          COR%CORRMP2DIR_K(I)=0
          omega_sim: DO NCURR = 1, NOMEGA_IN_ROOT_GROUP  !loop over frequencies within group
          !--------------------------------------------------------------------------------------
             IF ( NCURR <= NOMEGA_SIMULTANEOUS ) THEN
                !this is the global omega point treated currently in each group 
                NOMEGA_GLOBAL = LOCAL_INDEX_TO_GLOBAL( NCURR, COMM_BETWEENOMEGA%NODE_ME, &
                    COMM_BETWEENOMEGA%NCPU, NOMEGA)
                NOMEGA_CURRENT = NCURR
                LDMODE = .FALSE.
             ELSE
                !group having less points than the root group use their largest local frequency
                NOMEGA_GLOBAL = LOCAL_INDEX_TO_GLOBAL( NOMEGA_SIMULTANEOUS, COMM_BETWEENOMEGA%NODE_ME, &
                    COMM_BETWEENOMEGA%NCPU, NOMEGA)

                NOMEGA_CURRENT = NOMEGA_SIMULTANEOUS
                LDMODE = .TRUE. 
             ENDIF

             !-----------------------------------------------------------------------------------
             ! directions only need to be considered for the head 
             IF ( CHI%LGAMMA ) THEN
             directions: DO IDIR = 1, IDIR_MAX          !loop over the tree directions  
             !-----------------------------------------------------------------------------------
                !for the long-wave limit we copy head and wings to RESPONSEFUN
                CALL BODY_FROM_WING_GG( CHI, IDIR, NOMEGA_CURRENT, NROWS, NCOLS, COMM_INOMEGA)

                !CHI is distributed in groups to all processors
                !each group contains CHI at one frequency point. 
                !multiplication with the Coulomb potential needs to be done carefully,
                !CHI will be used for GW calculations, therefore
                !we copy the matrix to CHI_WORK and calculate the RPA correlation energy
                !using CHI_WORK 

                !first we allocate CHI_WORK in the same way as CHI
                ALLOCATE(CHI_WORK(NROWS,NCOLS)) 
                !multiply with v^1/2 from left and right
                !furthermore set the correct energy cutoff
                CHI_WORK = 0
                CALL XI_LOCAL_FIELD_ACFDT_GG( CHI_WORK, CHI, NROWS, NCOLS,                 &
                   COMM_INOMEGA, WDES%COMM_KIN, LATT_CUR, COR%ENCUTGW(I), COR%ENCUTGWSOFT(I),  &
                   NOMEGA, NOMEGA_SIMULTANEOUS, NOMEGA_CURRENT, WGWQ , NP_NEW, NODE_CRITICAL,      &
                   .FALSE., I, IO%IU0 )
                !scaLAPACK block size and processor grid is set here
                CALL BLOCKSIZE_AND_PROC_GRID( NROWS, NCOLS, COMM_INOMEGA%NCPU,-1)

                !use scaLAPACK to compute eigenvalues of CHI_WORK
                CALL COL2BC_AND_RPA(CHI_WORK, NROWS, NCOLS, NP_NEW, COR,                   &
                   COMM_INOMEGA, COMM_BETWEENOMEGA, WDES%COMM_KIN, NOMEGA_CURRENT, NOMEGA_GLOBAL,      &
                   NOMEGA_SIMULTANEOUS, IMAG_GRIDS%BOS_RE_WEIGHT, I, IO,                                &
                   CHI%LREALSTORE, LDMODE , CHI%LGAMMA, IDIR_MAX ) 

                DEALLOCATE ( CHI_WORK)     
             !-----------------------------------------------------------------------------------
             ENDDO directions  
             ELSE 
             !-----------------------------------------------------------------------------------
                !CHI is distributed in groups to all processors
                !each group contains CHI at one frequency point. 
                !multiplication with the Coulomb potential needs to be done carefully,
                !CHI will be used for GW calculations, therefore
                !we copy the matrix to CHI_WORK and calculate the RPA correlation energy
                !using CHI_WORK 
                
                !first we allocate CHI_WORK in the same way as CHI
                ALLOCATE(CHI_WORK(NROWS,NCOLS)) 
                CHI_WORK = 0

                !multiply with v^1/2 from left and right
                !furthermore set the correct energy cutoff
                CALL XI_LOCAL_FIELD_ACFDT_GG( CHI_WORK, CHI, NROWS, NCOLS,                 &
                   COMM_INOMEGA, WDES%COMM_KIN, LATT_CUR, COR%ENCUTGW(I), COR%ENCUTGWSOFT(I),  &
                   NOMEGA, NOMEGA_SIMULTANEOUS, NOMEGA_CURRENT, WGWQ , NP_NEW, NODE_CRITICAL,      &
                   .FALSE., I, IO%IU0 )

                !scaLAPACK block size and processor grid is set here
                CALL BLOCKSIZE_AND_PROC_GRID( NROWS, NCOLS, COMM_INOMEGA%NCPU,-1)

                !use scaLAPACK to compute eigenvalues of CHI_WORK
                CALL COL2BC_AND_RPA(CHI_WORK, NROWS, NCOLS, NP_NEW, COR,                   &
                   COMM_INOMEGA, COMM_BETWEENOMEGA, WDES%COMM_KIN, NOMEGA_CURRENT, NOMEGA_GLOBAL,      &
                   NOMEGA_SIMULTANEOUS, IMAG_GRIDS%BOS_RE_WEIGHT, I, IO,                                &
                   CHI%LREALSTORE, LDMODE, CHI%LGAMMA ) 
                 
                DEALLOCATE ( CHI_WORK)     
             !-----------------------------------------------------------------------------------
             ENDIF
          !--------------------------------------------------------------------------------------
          ENDDO omega_sim
       !-----------------------------------------------------------------------------------------
       ENDDO cutoff
       !-----------------------------------------------------------------------------------------
   
       !finally communicate the result between nodes (note that COR%CORR... is complex)
       CALLMPI( M_sum_d(COMM_BETWEENOMEGA, COR%CORRELATION_K, 2*NE ))  !RPA
       CALLMPI( M_sum_d(COMM_BETWEENOMEGA, COR%CORRMP2DIR_K, 2*NE ))   !MP2

       !scale appropriately  
       COR%CORRELATION_K=COR%CORRELATION_K/2*KPOINTS_ORIG%WTKPT(NQ) 
       COR%CORRMP2DIR_K =COR%CORRMP2DIR_K/2*KPOINTS_ORIG%WTKPT(NQ) 
       !save the result ot CORELATION and MP2DIR      
       COR%CORRELATION=COR%CORRELATION  +COR%CORRELATION_K
       COR%CORRMP2DIR =COR%CORRMP2DIR   +COR%CORRMP2DIR_K
                 
    !--------------------------------------------------------------------------------------------
    ENDDO lamb
   !---------------------------------------------------------------------------------------------

    PROFILING_STOP('calculate_rpa_correlation_energy')

 RETURN
END SUBROUTINE CALCULATE_RPA_CORRELATION_ENERGY

!******************************* XI_LOCAL_FIELD_ACFDT_GG ********************************
!
! this subroutine truncates the response function by multiplying
! with a truncated  Coulomb kernel smoothly going from 1/G^2 to zero
! between ENCUTGWSOFT and ENCUTGW
!
!****************************************************************************************

  SUBROUTINE XI_LOCAL_FIELD_ACFDT_GG( CHI_WORK, CHI, NROWS, NCOLS, &
    COMM_INOMEGA, GLOBALCOMM, LATT_CUR, ENCUT, ENCUTSOFT, &
    NOMEGA, NOMEGA_SIMULTANEOUS, NCURR, WGWQ , MAXINDEX, &
    NODE_CRITICAL, LSR,LDD,INU)
    USE constant
    USE fock
    USE mpimy
    USE c2f_interface, ONLY : ERRF
    IMPLICIT NONE
    GDEF,POINTER  :: CHI_WORK(:,:)
    TYPE (responsefunction) :: CHI
    TYPE (latt) LATT_CUR
    INTEGER :: NCOLS, NROWS
    TYPE(communic) :: COMM_INOMEGA     !communicator in one group
    TYPE(communic) :: GLOBALCOMM     !global communicator
    REAL(q) :: ENCUT, ENCUTSOFT
    INTEGER :: NOMEGA
    INTEGER :: NOMEGA_SIMULTANEOUS   !# of frequencies in one group
    INTEGER :: NCURR                  !current freuquency point
    TYPE (wavedes1) :: WGWQ
    INTEGER :: MAXINDEX              ! number of grid points using truncated Coulomb kernel
    INTEGER :: NODE_CRITICAL         ! node id containing cutoff
    LOGICAL :: LSR
    INTEGER :: LDD,INU
  ! local
    INTEGER    NI, NP, NI_, N1, I, J, I_
    REAL(q) :: DKX, DKY, DKZ, GX, GY, GZ, GSQU, SCALE, POTFAK, E
    INTEGER :: NI_LOC, NOLD_LOC, NODE_RECV, IERROR
    REAL(q), ALLOCATABLE :: DATAKE(:)
    INTEGER, ALLOCATABLE :: OLD_INDEX(:)
    INTEGER :: ISEND(4) , IRECV(4)
    GDEF, ALLOCATABLE    :: CHI_RECV(:)
    GDEF, ALLOCATABLE    :: CHI_SEND(:)
    REAL(q) :: DFUN, SFUN
    ! neu
    REAL(q) :: OMEGBK,QC
    REAL(q) :: GAMMASCALE

    ! multipole correction
    INTEGER    NJ,NK
    COMPLEX(q) :: TSUM(10)
    COMPLEX(q) :: POTCORRECTION(10)
    !for communication
    INTEGER    :: NI_GLOB, NSEND
    INTEGER    :: NODE_SEND

    REAL(q) :: Q1, Q2

    IF (ENCUTSOFT>=0) THEN
       Q1=SQRT(ENCUTSOFT/HSQDTM)
       Q2=SQRT(ENCUT/HSQDTM)
    ENDIF
    ! e^2/ volume
    ! mK scales by EDEPS/LATT_CUR%OMEGA here
    SCALE=EDEPS/LATT_CUR%OMEGA
!    SCALE=1.0_q
!=======================================================================
! first set up the truncated Coulomb kernel
! smoothly going from 1/G^2 to zero between ENCUTGWSOFT and ENCUTGW
!=======================================================================
    DKX=(WGWQ%VKPT(1))*LATT_CUR%B(1,1)+ &
        (WGWQ%VKPT(2))*LATT_CUR%B(1,2)+ &
        (WGWQ%VKPT(3))*LATT_CUR%B(1,3)
    DKY=(WGWQ%VKPT(1))*LATT_CUR%B(2,1)+ &
        (WGWQ%VKPT(2))*LATT_CUR%B(2,2)+ &
        (WGWQ%VKPT(3))*LATT_CUR%B(2,3)
    DKZ=(WGWQ%VKPT(1))*LATT_CUR%B(3,1)+ &
        (WGWQ%VKPT(2))*LATT_CUR%B(3,2)+ &
        (WGWQ%VKPT(3))*LATT_CUR%B(3,3)

    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2

    ALLOCATE( DATAKE(NROWS), OLD_INDEX(NROWS))
    DO NI = 1, NROWS
     OLD_INDEX(NI) = NI
    ENDDO

    !treated frequency point
    N1=NCURR

    DATAKEMAX=0.0
    MAXINDEX=0
    DO NI=1,NP
       NI_=NI
       IF (WGWQ%LGAMMA) NI_=(NI-1)/2+1
       
       GX=(WGWQ%IGX(NI_)*LATT_CUR%B(1,1)+WGWQ%IGY(NI_)* &
            LATT_CUR%B(1,2)+WGWQ%IGZ(NI_)*LATT_CUR%B(1,3))
       GY=(WGWQ%IGX(NI_)*LATT_CUR%B(2,1)+WGWQ%IGY(NI_)* &
            LATT_CUR%B(2,2)+WGWQ%IGZ(NI_)*LATT_CUR%B(2,3))
       GZ=(WGWQ%IGX(NI_)*LATT_CUR%B(3,1)+WGWQ%IGY(NI_)* &
            LATT_CUR%B(3,2)+WGWQ%IGZ(NI_)*LATT_CUR%B(3,3))
          
       GSQU=(DKX+GX)**2+(DKY+GY)**2+(DKZ+GZ)**2
       

       IF (ABS(GSQU)<G2ZERO) THEN
          ! head and wing
          POTFAK=SCALE
          ! if HFRCUT is set, a fixed spherical cutoff is used
          ! resulting in a finite values at G=0
          ! since the head stores the q^2, the final contribution is zero
          ! from the head
          IF (HFRCUT/=0) THEN
             POTFAK=0
          ! if LRHFCALC.AND.LRSCOR the Coulomb kernel is replaced by v_{LR}=v*exp(-q^2/(4*mu^2))
          ! for q=0  v_{LR}=1, this introduces an error for small mu because of the finite q-grid
          ! In order to correct this error, v_{LR}(q=0) is set to average of v_{LR} in a sphere  
          ! with volume corresponding to BZ volume/ # q-points       
          ELSE IF (LRHFCALC.AND.LRSCOR) THEN
             CALL CELVOL(LATT_CUR%B(1,1),LATT_CUR%B(1,2),LATT_CUR%B(1,3),OMEGBK)
             QC=TPI*(0.75_q/PI*OMEGBK/KPOINTS_FULL%NKPTS)**(1._q/3._q)
             GAMMASCALE=8*PI*HFSCREEN*HFSCREEN* &
            &   (-EXP(-QC*QC/4/HFSCREEN/HFSCREEN)*QC+HFSCREEN*SQRT(PI)*ERRF(QC/2/HFSCREEN))/ &
            &   (OMEGBK*TPI*TPI*TPI/KPOINTS_FULL%NKPTS)
             POTFAK=POTFAK*GAMMASCALE
          ! if LSR the Coulomb kernel is replaced by v_{SR}=v*[1-exp(-q^2/(4*mu^2))]
          ! for q=0  v_{SR}=0, this introduces an error for small mu because of the finite q-grid
          ! In order to correct this error, v_{SR}(q=0) is set to average of v_{SR} in a sphere
          ! with volume corresponding to BZ volume/ # q-points
          ELSE IF (LSR) THEN
             CALL CELVOL(LATT_CUR%B(1,1),LATT_CUR%B(1,2),LATT_CUR%B(1,3),OMEGBK)
             QC=TPI*(0.75_q/PI*OMEGBK/KPOINTS_FULL%NKPTS)**(1._q/3._q)
             GAMMASCALE=8*PI*HFSCREEN*HFSCREEN* &
            &   (-EXP(-QC*QC/4/HFSCREEN/HFSCREEN)*QC+HFSCREEN*SQRT(PI)*ERRF(QC/2/HFSCREEN))/ &
            &   (OMEGBK*TPI*TPI*TPI/KPOINTS_FULL%NKPTS)
             POTFAK=POTFAK*(1-GAMMASCALE)
          ENDIF

          ! switch off standard convergence correction
          IF(MCALPHA/=0) THEN
             POTFAK=0
          ENDIF

       ELSE
          ! the factor 1/(2 pi)^2 is required to obtain proper reciprocal
          ! lattice vector lenght
          POTFAK=SCALE/(GSQU*TPI**2)
          IF (HFRCUT/=0) THEN
             ! spherical cutoff on Coloumb kernel
             ! see for instance C.A. Rozzi, PRB 73, 205119 (2006)
!             POTFAK=POTFAK*(1-COS(SQRT(GSQU)*TPI*HFRCUT)*EXP(-(SQRT(GSQU)*TPI*HFRCUT)**2*HFRCUT_SMOOTH))
!test
             CALL DELSTP(3,SQRT(GSQU)*TPI*HFRCUT/10,DFUN,SFUN)
             POTFAK=POTFAK*(SFUN-0.5)*2
!test
          ELSE IF (LSR) THEN
             ! use the short range part of the Coulomb kernel
             POTFAK=POTFAK*(1._q-EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN))))
          ELSE IF (LRHFCALC.AND.LRSCOR) THEN
             ! use the long range part of the Coulomb kernel
             POTFAK=POTFAK*EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN)))
          ENDIF
       ENDIF

       ! smooth cutoff function between  ENCUTSOFT and ENCUT
       E=HSQDTM*(GSQU*TPI**2)
       IF (ENCUT>=0 .AND. E>ENCUT) THEN
          POTFAK=0
       ELSE

          MAXINDEX=MAXINDEX+1
          OLD_INDEX(MAXINDEX)=NI
          
          IF (ENCUTSOFT>=0 .AND. E>ENCUTSOFT) THEN
           ! cosine window
             IF ( .NOT. LSCK) THEN
                POTFAK=POTFAK*(1+COS((E-ENCUTSOFT)/(ENCUT-ENCUTSOFT)*PI))/2
           ! squeezed Coulomb kernel
             ELSE 
                POTFAK=POTFAK*SQUEEZED_COULOMB_KERNEL(SQRT(GSQU)*TPI, Q1, Q2 )
             ENDIF
          ENDIF
       ENDIF

! POTFAK in fock has an additional factor. Since the corrections were originally devised
! for the fock routine this ensures consistency
       IF(MCALPHA/=0) THEN
          POTFAK=POTFAK*(1.0_Q/WGWQ%GRID%NPLWV)
       ENDIF

       DATAKE(NI)=SQRT(POTFAK)
       ! maximum kinetic energy
       IF(DATAKE(NI)>DATAKEMAX) THEN
          DATAKEMAX=DATAKE(NI)
       ENDIF

    ENDDO

    ! setup of multipole corrections
    IF (MCALPHA/=0) THEN
       CALL FOCK_MULTIPOLE_CORR_SETUP(LATT_CUR, GRIDHF)
       CALL GW_MULTIPOLE_PROJ_SETUP( WGWQ, LATT_CUR, DATAKE,DATAKEMAX )
    ENDIF

!=======================================================================
! now multiply the response function with the kernel
! from left and right hand side
!=======================================================================

    !save chi to chi_work
    ! change sign of RESPONSEFUN  here 
    !CHI_WORK=0._q
    IF (CHI%LREALSTORE) THEN
          CHI_WORK(1:NROWS,1:NCOLS) = -CHI%RESPONSER(   1 : NROWS , 1 : NCOLS , N1 )
    ELSE
          CHI_WORK(1:NROWS,1:NCOLS) = -CHI%RESPONSEFUN( 1 : NROWS , 1 : NCOLS , N1 )
    ENDIF

!WRITE(500+GLOBALCOMM%NODE_ME,'(" REAL PART CHI_WORK for NOMEGA_SIMULTANEOUS",I2 )')NCURR
!DO I_ = 1 , NROWS
!   WRITE(500+GLOBALCOMM%NODE_ME,'(32F10.4)')&
!      REAL(CHI_WORK(I_,1:NCOLS  ))
!ENDDO

    !kill elements larger than NP
    IF ( NP < NROWS ) THEN
     DO NI = NP + 1, NROWS
      DATAKE(NI) = 0
     ENDDO
    ENDIF

    !DATAKE(:) has the coulomb potenial, this needs to be multiplied with CHI_WORK
    !and symmetrized
    DO NI=1,NROWS  !from the left
       POTFAK=DATAKE( NI )
       CHI_WORK(NI,1:NCOLS)=CHI_WORK(NI,1:NCOLS)*POTFAK
    ENDDO

    !from right
    DO NI=1,NCOLS
       POTFAK=DATAKE( NI + (COMM_INOMEGA%NODE_ME-1)*NCOLS )
       !IF ( NI + (COMM_INOMEGA%NODE_ME-1)*NCOLS > NP ) POTFAK = 0 
       CHI_WORK(1:NROWS,NI)=CHI_WORK(1:NROWS,NI)*POTFAK
    ENDDO


    ! multipole corrections:
    IF(MCALPHA/=0) THEN 

      CALL vtutor%error("Sorry MCALPHA not implemented")

       DO NI=1,NP
          
          TSUM=0.0
          DO NK=1,10
             DO NJ=1,NP
                TSUM(NK)=TSUM(NK)+CHI_WORK(NI,NJ)*CONJG(P_PROJ2(NJ,NK))                
             ENDDO
          ENDDO

          POTCORRECTION=(MATMUL(AAH_MAT,TSUM))

          DO NK=1,10
             DO NJ=1,NP
                CHI_WORK(NI,NJ)=CHI_WORK(NI,NJ)+P_PROJ2(NJ,NK)*POTCORRECTION(NK)
             ENDDO
          ENDDO
          
       ENDDO
       
       
       DO NI=1,NP
          
          TSUM=0.0
          DO NK=1,10
             DO NJ=1,NP
                TSUM(NK)=TSUM(NK)+CHI_WORK(NJ,NI)*P_PROJ2(NJ,NK)
             ENDDO
          ENDDO

          POTCORRECTION=(MATMUL(AAH_MAT,TSUM))

          DO NK=1,10
             DO NJ=1,NP
                CHI_WORK(NJ,NI)=CHI_WORK(NJ,NI)+CONJG(P_PROJ2(NJ,NK))*POTCORRECTION(NK)
             ENDDO
          ENDDO
          
       ENDDO
! When using MCALPHA POTFAK had an additional normalisation factor to ensure consistency with fock.F
! acfdt does not need this extra factor so we have to take it out again:
       CHI_WORK=CHI_WORK*(WGWQ%GRID%NPLWV)
    ENDIF

   
    !smooth energy cutoff done by replacing the NI^th coulmn by the OLD_INDEX(NI)^th column 
    !and analogously for the row. 
#ifdef debug
    IF (INU>=0 ) WRITE(*,'(" Shrinking started for cut ",2I4)')LDD
#endif 

    !The replacement of the rows is straight forward 
    !and is be done directly by the following loop
    !IF (GLOBALCOMM%NODE_ME.eq.1) WRITE (*,*) 'test CHI_WORK ', SUM(CHI_WORK)
    DO NI=1,MAXINDEX
      CHI_WORK(NI,1:NCOLS)=CHI_WORK(OLD_INDEX(NI),1:NCOLS)
      IF (NI.ne.OLD_INDEX(NI)) CHI_WORK(OLD_INDEX(NI),:) =0._q
    ENDDO
    !added to test the dMP2 result 
    !CHI_WORK(MAXINDEX+1: ,:)=0._q
    !IF (GLOBALCOMM%NODE_ME.eq.1) WRITE (*,*) 'test CHI_WORK ', SUM(CHI_WORK)

    !for the columns this is not straight forward
    !since the whole matrix is distributed in stripes of NCOLS
    !among the nodes in each frequency group.

    !send columns between nodes
    ALLOCATE(CHI_SEND( NROWS ) )  
    ALLOCATE(CHI_RECV( NROWS ) ) 
    CHI_RECV = 0 
    CHI_SEND = 0 

    !determine which node in group has the critical stripe.
    DO NI = 1, COMM_INOMEGA%NCPU 
       IF (MAXINDEX > NI*NCOLS) CYCLE
       IF ( (NI-1)*NCOLS < MAXINDEX .AND. MAXINDEX <= NI*NCOLS) THEN 
          NODE_CRITICAL = NI 
          EXIT
       ENDIF
    ENDDO

#ifdef debug
 WRITE(100+WGWQ%COMM%NODE_ME,'("CUTOFF #",I4," NODE",I2," CRIT.NODE",I2," NROWS",I4," NP_NEW",I4)')&
      LDD, WGWQ%COMM%NODE_ME,NODE_CRITICAL, NROWS,MAXINDEX 
#endif
    CALL MPI_BARRIER(COMM_INOMEGA%MPI_COMM,IERROR)   !each node should start at the same time from here
    !==============================================================================================
    !we loop over the cpus in one group and let other node send the appropriate columns
    DO NODE_RECV = 1, COMM_INOMEGA%NCPU
    !==============================================================================================
        loc_col : DO NI = 1 , NCOLS                      !search column which should be replaced
        !------------------------------------------------------------------------------------------
 
           NI_GLOB = NI + (NODE_RECV-1)*NCOLS   !global index for current node NODE_RECV
  
 
           !find node NODE_SEND containing OLD_INDEX(NI_GLOB)
           DO NODE_SEND = 1, COMM_INOMEGA%NCPU
             IF ( 1 + ( NODE_SEND-1)*NCOLS <= OLD_INDEX(NI_GLOB) .AND. &
                  OLD_INDEX(NI_GLOB) <= NODE_SEND*NCOLS) EXIT
           ENDDO
 
           NSEND = OLD_INDEX(NI_GLOB) - (NODE_SEND-1)*NCOLS   !local column index on sending node
 
           IF ( NSEND < 1 ) THEN 
              CALL vtutor%error("Cannot send column with negative index " // str(NODE_CRITICAL) // "&
                 & " // str(NODE_RECV) // " " // str(NODE_SEND) // " " // str(NI_GLOB) // " " // &
                 str(OLD_INDEX(NI_GLOB)) // " " // str(NSEND) // " " // str(NCOLS))
           ENDIF
 
           
           !if old column is on the same node simply replace it 
           IF ( NODE_SEND == COMM_INOMEGA%NODE_ME .AND. NODE_RECV == COMM_INOMEGA%NODE_ME ) THEN 
#if debug  
              WRITE(100+WGWQ%COMM%NODE_ME,&
                  '( " COL ON SAME CPU #",I4," COL #",2I4, "REPLACED BY COL #", I4,10F10.5)')&
                  NODE_RECV, NI,NI_GLOB, NSEND,CHI_WORK(1:10,NSEND)
#endif
              CHI_WORK(1:NROWS, NI) = CHI_WORK(1:NROWS, NSEND)
              !jK added
              IF (NSEND.ne.NI) CHI_WORK(:,NSEND) = 0._q

           !otherwise tell NODE_SEND to send the appropriate column NSEND to NODE_RECV
           ELSEIF ( NODE_RECV < NODE_SEND ) THEN
              IF ( NODE_SEND == COMM_INOMEGA%NODE_ME ) THEN
#if debug  
                 WRITE(100+WGWQ%COMM%NODE_ME,'( " COL #",I4, "SENT TO COL", I4 , " ON CPU #", I4)')&
                 NSEND,NI,NODE_RECV
#endif
                 CHI_SEND(1:NROWS) = CHI_WORK(1:NROWS, NSEND)
                 CHI_WORK(:,NSEND) = 0._q

                 !send this column to node NODE_RECV              
                 IF (CHI%LREALSTORE) THEN
                    CALLMPI ( M_send_d( COMM_INOMEGA, NODE_RECV, CHI_SEND, NROWS) )
                 ELSE
                    CALLMPI ( M_send_z( COMM_INOMEGA, NODE_RECV, CHI_SEND, NROWS) )
                 ENDIF

              ELSEIF ( NODE_RECV == COMM_INOMEGA%NODE_ME ) THEN 
#if debug  
                 WRITE(100+WGWQ%COMM%NODE_ME,'( " RECV COL #", I4 ," FROM CPU #"'//&
                    ',I4, " FOR COL", I4)') NSEND, NODE_SEND, NI
#endif
                 IF (CHI%LREALSTORE) THEN
                    CALLMPI ( M_recv_d( COMM_INOMEGA, NODE_SEND, CHI_RECV, NROWS) )
                 ELSE
                    CALLMPI ( M_recv_z( COMM_INOMEGA, NODE_SEND, CHI_RECV, NROWS) )
                 ENDIF
                     
                 !and replace CHI_RECV 
                 CHI_WORK(1:NROWS, NI ) = CHI_RECV(1:NROWS)
              ENDIF
           ENDIF

        !mpi barrier here 
        CALL MPI_BARRIER(COMM_INOMEGA%MPI_COMM,IERROR)   
        !------------------------------------------------------------------------------------------
        ENDDO loc_col                                !local columns  
    !==============================================================================================
    ENDDO                                    !cpus in group
    !==============================================================================================
    DEALLOCATE( CHI_SEND )
    DEALLOCATE( CHI_RECV )

    DEALLOCATE(DATAKE, OLD_INDEX)
    RETURN
 END SUBROUTINE XI_LOCAL_FIELD_ACFDT_GG


!****************************************************************************************
!
! redistributes column wise distributed CHI_WORK  2D block-cyclic wise 
! also the RPA and MP2 correlation energy is determined.
!
!****************************************************************************************

SUBROUTINE COL2BC_AND_RPA(CHI_WORK, NROWS, NCOLS, NP_NEW, COR, &
    COMM_INOMEGA, COMM_BETWEENOMEGA, GLOBALCOMM, NOMEGA_CURRENT, &
    NOMEGA_GLOBAL, NOMEGA_SIMULTANEOUS, OMEGAWEIGHT, NCUT, IO ,&
    LREALC , LDMODE, LGAMMA, IDIR_MAX)
    USE base
    USE constant
    USE mpimy 
    USE ini
    USE scala
    USE chi_glb, ONLY: LLTDMP2 
    IMPLICIT NONE
    GDEF, POINTER, DIMENSION(:,:)    :: CHI_WORK    !matrix to be redistributed
    INTEGER                 :: NROWS, NCOLS         !number of local rows and cols of CHI
    INTEGER                 :: NP_NEW               !dimension of submatrix
    TYPE (correlation)      :: COR                  !correlation energy
    TYPE (communic)         :: COMM_INOMEGA         !communicator in tau
    TYPE (communic)         :: COMM_BETWEENOMEGA    !communicator between tau
    TYPE (communic)         :: GLOBALCOMM           !global communicator
    INTEGER                 :: NOMEGA_CURRENT       !current frequency point in group
    INTEGER                 :: NOMEGA_GLOBAL        !global frequency point treated by group  
    INTEGER                 :: NOMEGA_SIMULTANEOUS  !frequencies on one node
    REAL(q)                 :: OMEGAWEIGHT(:)       !integeration weights
    INTEGER                 :: NCUT                 !current cut off 
    TYPE (in_struct)        :: IO                   !IO 
    LOGICAL                 :: LREALC               !real or complex calculation
    LOGICAL                 :: LDMODE               !dummy mode?
    LOGICAL                 :: LGAMMA               !doing longwave limit?
    INTEGER, OPTIONAL       :: IDIR_MAX             !dimensions for longwave limit 
    !local                
    GDEF, POINTER, DIMENSION(:,:)  :: CHI_TMP(:,:)  !redistributed chi
    INTEGER                 :: I,J                  !some loop variables 
    INTEGER                 :: NROWS_NEW, NCOLS_NEW !new size of redistributed chi 
    INTEGER                 :: NROWS_USED, NCOLS_USED 
    INTEGER                 :: IERROR,INFO          !error variables
    INTEGER                 :: ICTXA                !contextes for BLACS (column major)
    INTEGER                 :: ICTXB                !contextes for BLACS (block cyclic)
    REAL(q), ALLOCATABLE    :: E(:)                 !eigenvectors and global eigenvalues
    INTEGER                 :: PINFO                !scaLAPACK routine info variable
    INTEGER                 :: LWORK, LRWORK        !for scaLAPACK routine 
    INTEGER                 :: LLD_A, LLD_B         !leading (row) dimensions for desciptors  
    GDEF, ALLOCATABLE       :: WORK(:), EV( :, : )  !working array and eigenvectors 
    COMPLEX(q), ALLOCATABLE :: RWORK(:)             !needed for complex computation 
    REAL(q)                 :: RTMP                 !auxilairy 
    REAL(q)                 :: SUM2                 !stores sum of squares of CHI_WORK elements 
                                                    !to obtain dMP2 contribution
    COMPLEX(q)              :: CHI_TRACE            !trace of CHI_WORK to get RPA energy
    LOGICAL                 :: LRPADIAG             !use diagonalisation or Cholesky
    REAL*8, external        :: PDLATRA  
    COMPLEX*16, external    :: PZLATRA  

    PROFILING_START( 'col2bc_and_rpa' )

    ! if linear term is calculated use diagonalization routine
    ! usually linear term is not calculated, so Cholesky decomposition is used
    ! by default
    LRPADIAG = LDMP1 

    !LRPADIAG=.FALSE.               ! get energy from Cholesky for RPA and from \sum |A_ij|^2 for dMP2 
    !LRPADIAG=.TRUE.               ! get energy from diagonalisation 


    IF ( LGAMMA .AND. .NOT.(PRESENT(IDIR_MAX)) ) THEN
          CALL vtutor%error("ERROR in COL2BC_AND_RPA: Longwave wanted but IDIR_MAX not present")
    ENDIF

    !---------------------------------------------------------------------------------------
    !initialize 1D column-major processor grid
!    ICTXA=GLOBALCOMM%MPI_COMM
!    CALL GEN_PROC_GROUP_GRIDS( COMM_BETWEENOMEGA%NODE_ME, COMM_INOMEGA%NCPU ,&
!       COMM_BETWEENOMEGA, ICTXA)

    ! row major ordering (but column-major ordering is identical if there is only 1 column)
    NPROW0=1
    NPCOL0=COMM_INOMEGA%NCPU
    CALL PROCMAP( COMM_INOMEGA, ICTXA, 2, NPROW0, NPCOL0)
    CALL BLACS_GRIDINFO( ICTXA, NPROW0, NPCOL0, MYROW0, MYCOL0) !get the ids for scaLAPACK

    !define corresponding array descriptor DESCA
    !since CHI_WORK is distributed blockwisely,  
    !the leading dimension of rows is simply NROWS, 
    !check this 
    LLD_A = MAX( 1, NUMROC(NROWS, NCOLS, MYROW0, 0, NPROW0) )
    IF ( LLD_A /= NROWS ) THEN 
       CALL vtutor%error("ERROR in COL2BC_AND_RPA: , LLD_A has wrong size for NODE " // &
          str(GLOBALCOMM%NODE_ME) // " " // str(LLD_A) // " " // str(NROWS))
    ENDIF 
    CALL DESCINIT(DESCA, NROWS , NROWS, NCOLS, NCOLS, 0, 0, ICTXA, &
       LLD_A, INFO )
    !Eventually take care for errors 
    IF ( INFO /= 0 ) THEN 
       CALL vtutor%error("Error in 1D col-major initialization " // str(INFO))
    ENDIF

!    ICTXB=GLOBALCOMM%MPI_COMM
!    !initialize 2D block-cyclic processor grid
!    CALL GEN_PROC_GROUP_GRIDS_BC( COMM_BETWEENOMEGA%NCPU, &
!                                  COMM_BETWEENOMEGA, COMM_INOMEGA, ICTXB)
    ! set NPROW,NPCOL globally
    CALL FERMAT_RAZOR( COMM_INOMEGA%NCPU, NPROW, NPCOL )
    ! generate processor grid
    CALL PROCMAP( COMM_INOMEGA, ICTXB, 1, NPROW, NPCOL)
    CALL BLACS_GRIDINFO( ICTXB, NPROW, NPCOL, MYROW, MYCOL) !get the ids for scaLAPACK

    !leading dimension of local array for block cyclic distribution
    LLD_B = MAX( 1 , NUMROC(NROWS, MB, MYROW, 0, NPROW) )

    !initialize the array descriptor for block cyclic distribution
    CALL DESCINIT(DESCB, NROWS , NROWS, MB, NB, 0, 0, ICTXB, &
       LLD_B , INFO )
    IF ( INFO /= 0 ) THEN 
       CALL vtutor%error("Error in 2D block-cyclic initialization " // str(INFO))
    ENDIF

    !check if CHI_WORK has proper size
    NROWS_NEW = NUMROC(DESCB(3),DESCB(5),MYROW,0,NPROW)
    NCOLS_NEW = NUMROC(DESCB(4),DESCB(6),MYCOL,0,NPCOL) 
 
    !the following is only neccessary if we are not in dummy mode
    IF( .NOT. LDMODE ) THEN

       !SUM2=0._q
       !DO I=1,SIZE(CHI_WORK,1)
       !  DO J=1,SIZE(CHI_WORK,2)
       !   SUM2=SUM2+CHI_WORK(I,J)*GCONJG(CHI_WORK(I,J))
      !   ENDDO
      ! ENDDO
      !IF (GLOBALCOMM%NODE_ME.eq.1) WRITE(*,*) 'MP2 from very old matrix ' , SUM2*0.5


       !allocate storage for redistributed matrix
       ALLOCATE(CHI_TMP(NROWS_NEW,NCOLS_NEW))
       CHI_TMP = 0 

       !now call the BLACS redistribution routine 
#ifdef gammareal
        CALL PDGEMR2D( NROWS, NROWS, CHI_WORK, 1, 1, DESCA, & 
                                     CHI_TMP, 1, 1, DESCB, &
                       ICTXB)
#else
        CALL PZGEMR2D( NROWS, NROWS, CHI_WORK, 1, 1, DESCA, & 
                                     CHI_TMP, 1, 1, DESCB, &
                       ICTXB)
#endif

#ifdef debug
       !test the new block-cyclic distribution 
       WRITE(100+GLOBALCOMM%NODE_ME + (NOMEGA_GLOBAL-1)*10,'( "REAL" ,3I4)')NROWS,NCOLS, NP_NEW
       DO I = 1 , NROWS
         WRITE(100+GLOBALCOMM%NODE_ME + (NOMEGA_GLOBAL-1)*10,'(32F10.5)') REAL(CHI_WORK(I,1:NCOLS))
       ENDDO
       WRITE(200+GLOBALCOMM%NODE_ME + (NOMEGA_GLOBAL-1)*10,'( "IMAG" ,3I4)')NROWS,NCOLS,NP_NEW
       DO I = 1 , NROWS
         WRITE(200+GLOBALCOMM%NODE_ME + (NOMEGA_GLOBAL-1)*10,'(32F10.5)') IMAG(CHI_WORK(I,1:NCOLS))
       ENDDO
#endif    

       !SUM2=0._q
       !DO I=1,SIZE(CHI_TMP ,1)
       !  DO J=1,SIZE(CHI_TMP ,2)
       !   SUM2=SUM2+CHI_TMP(I,J)*GCONJG(CHI_TMP(I,J))
       !  ENDDO
       !ENDDO
       !IF (GLOBALCOMM%NODE_ME.eq.1) WRITE(*,*) 'MP2 from old matrix    ', SUM2*0.5

       !first deallocate CHI_WORK
       DEALLOCATE(CHI_WORK)
       !and point CHI_WORK to CHI_TMP
       CHI_WORK => CHI_TMP

#ifdef debug
       WRITE(100+GLOBALCOMM%NODE_ME,'( "block cyclic ",2I4 )')NROWS_NEW,NCOLS_NEW
       DO I = 1 , NROWS_NEW
          WRITE(100+GLOBALCOMM%NODE_ME,'(32F10.5)') REAL(CHI_WORK(I,1:NCOLS_NEW))
       ENDDO
#endif 

       IF (.NOT. LRPADIAG) THEN 
          !!!!!!!
          ! obtain the dMP2 energy here
          !!!!!!!
          ! gK: actual matrix size of CHI_WORK is obviously NP_NEW
          ! is this different than NROWS, well I guess
          NROWS_USED = NUMROC(NP_NEW,DESCB(5),MYROW,0,NPROW)
          NCOLS_USED = NUMROC(NP_NEW,DESCB(6),MYCOL,0,NPCOL) 

          SUM2=0._q
          DO J=1,NCOLS_USED
             DO I=1,NROWS_USED
                SUM2=SUM2+CHI_WORK(I,J)*GCONJG(CHI_WORK(I,J))
             ENDDO
          ENDDO
          CALLMPI( M_sum_d(COMM_INOMEGA, SUM2 , 1 ))
       ENDIF

       IF (LRPADIAG) THEN
          !allocate eigenvalue and eigenvector arrays
          ALLOCATE( E( NROWS) ) ;  E = 0
          ALLOCATE( EV( NROWS_NEW , NCOLS_NEW ) ) ; EV = 0 

          !first get the optimal size of working array
          ALLOCATE(WORK(1))
          IF(LREALC) THEN
             !for real values only  
             CALL PDSYEV( 'N', 'L', NP_NEW, CHI_WORK(1,1), 1, 1, DESCB, E(1),&
                  EV(1,1), 1, 1, DESCB, WORK(1), -1, PINFO)
          ELSE
             !for a hermitian matrix 
             ALLOCATE(RWORK(1))
             CALL PZHEEV( 'N', 'L', NP_NEW, CHI_WORK(1,1), 1, 1, DESCB, E(1), &
                  EV(1,1), 1, 1, DESCB, WORK(1), -1, RWORK(1), -1, PINFO )
          ENDIF 

          IF ( PINFO /=0 ) THEN
             CALL vtutor%error("ERROR in COL2BC_AND_RPA for NCUT " // str(NCUT) // " 1st call of &
                &PDSYEV returns " // str(PINFO))
          ENDIF

          !LWORK is now the optimal size of WORK 
          LWORK=WORK(1)
          IF ( WORK(1)  /= 0 ) THEN   !if it's not zero reallocate WORK array
             DEALLOCATE(WORK)
             ALLOCATE(WORK(LWORK))
          ENDIF
          IF (.NOT. LREALC) THEN
             !in the complex version do the same with the RWORK array
             LRWORK=RWORK(1)
             IF ( LRWORK  /= 0 ) THEN   !if it's not zero reallocate WORK array
               DEALLOCATE(RWORK)
               ALLOCATE(RWORK(LRWORK))
             ENDIF
          ENDIF

       !timing is of interest
!       CALL START_TIMING("PD")


          !IF (GLOBALCOMM%NODE_ME.eq.1) WRITE (*,*) 'chi_work size ', NROWS_NEW, NCOLS_NEW, SIZE(CHI_WORK,1), SIZE(CHI_WORK,2), NROWS, NCOLS, NP_NEW !, NPROW, NPCOL, MYROW, MYCOL

          !diagonalize 
          IF (LREALC) THEN
             !real version
             CALL PDSYEV( 'N', 'L', NP_NEW, CHI_WORK(1,1), 1, 1, DESCB, E(1),&
                  EV(1,1), 1, 1, DESCB, WORK(1), LWORK , PINFO)
          ELSE
             !hermitian version 
             CALL PZHEEV( 'N', 'L', NP_NEW, CHI_WORK(1,1), 1, 1, DESCB, E(1),&
                  EV(1,1), 1, 1, DESCB, WORK(1), LWORK, RWORK(1), LRWORK, PINFO )
          ENDIF 
          IF ( PINFO < 0 ) THEN
             CALL vtutor%error("ERROR in COL2BC_AND_RPA for NCUT " // str(NCUT) // " 2nd call of &
                &PDSYEV returns " // str(PINFO))
          ENDIF
          IF ( PINFO > 0 ) THEN
             IF( IO%IU0>=0) THEN
                WRITE(*,'("ERROR in COL2BC_AND_RPA for NCUT",I2," 2nd call of PDSYEV returns ",I4)')&
                   NCUT,PINFO
                WRITE(*,'(" Eigenvector #",I4," may not be converged" )') PINFO
             ENDIF
          ENDIF
       ELSE
          ! use Cholesky decomposition
          !calculate trace
          !add one to the diagonal
          CHI_TRACE=DTRACE_AND_ADD( NP_NEW, CHI_WORK, 1, 1, DESCB) 
          IF (LREALC) THEN
             !real version 
             !call Cholesky
             CALL PDPOTRF( 'U', NP_NEW, CHI_WORK(1,1), 1, 1, DESCB, PINFO)
          ELSE
             !complex version
             !call Cholesky
             CALL PZPOTRF( 'U', NP_NEW, CHI_WORK(1,1), 1, 1, DESCB, PINFO)
          ENDIF
          !do log of diagonal
          RTMP=2._q*DLOG_DIAG( NP_NEW, CHI_WORK, 1, 1, DESCB)
          !energy obtained as 2\sum_i ln(A)_{ii}
          !IF (GLOBALCOMM%NODE_ME.eq.1)  WRITE (*,*) 'RPA from element add     ', RTMP-CHI_TRACE
          !IF (GLOBALCOMM%NODE_ME.eq.1)  WRITE (*,*) 'dMP2 from element add     ', 0.5*SUM2

          RTMP=(RTMP-CHI_TRACE)* OMEGAWEIGHT( NOMEGA_GLOBAL )
          IF ( LGAMMA ) RTMP = RTMP / IDIR_MAX
          COR%CORRELATION_K( NCUT ) = COR%CORRELATION_K( NCUT ) + DBLE(RTMP)

          SUM2=0.5_q*SUM2* OMEGAWEIGHT( NOMEGA_GLOBAL )
          IF ( LGAMMA ) SUM2 = SUM2 / IDIR_MAX
          COR%CORRMP2DIR_K( NCUT ) = COR%CORRMP2DIR_K( NCUT ) - SUM2
 
       ENDIF

!       CALL STOP_TIMING("PD",IO%IU6,"PDSYEVD")

       !exit blacs  
       CALL BLACS_GRIDEXIT( ICTXB )  
       CALL BLACS_GRIDEXIT( ICTXA )  

       IF (LRPADIAG) THEN
          !this is the trace (note that E has wrong sign)    
          DO I=1, NROWS
             IF ( 1 + E(I) > 1E-8_q  ) THEN
                !RPA (Pines & Noizeres respectively Gell-Mann & Brueckner)
                RTMP = ( LOG( 1 + E( I ) ) - E( I ) ) * OMEGAWEIGHT( NOMEGA_GLOBAL ) 
                IF ( LGAMMA ) RTMP = RTMP / IDIR_MAX
                COR%CORRELATION_K( NCUT ) = COR%CORRELATION_K( NCUT ) + RTMP
   
                ! linear term
                IF ( LDMP1 ) THEN
                   RTMP = ( E( I ) ) * OMEGAWEIGHT( NOMEGA_GLOBAL )
                !direct Moeller Plessett term of second order
                ELSE
                   RTMP = ( E( I )*E( I ) / 2 ) * OMEGAWEIGHT( NOMEGA_GLOBAL )
                ENDIF
                IF ( LGAMMA ) RTMP = RTMP / IDIR_MAX
                COR%CORRMP2DIR_K( NCUT )  = COR%CORRMP2DIR_K( NCUT ) - RTMP
             ELSE 
                IF (COMM_INOMEGA%NODE_ME==COMM_INOMEGA%IONODE) &
                   WRITE(*,'(" Error in COL2BC_AND_RPA, eigenvalue #",I4," is negative '//&
                   'for group",I3,":",F12.5)') I, COMM_BETWEENOMEGA%NODE_ME, E( I )
             ENDIF
          ENDDO
          !RTMP=0._q
          !DO I=1,NROWS
          !   RTMP=RTMP+E(I)*E(I)
          !ENDDO
          !IF (GLOBALCOMM%NODE_ME.eq.1) WRITE (*,*) 'dMP2 from diagonalisation ',  0.5*RTMP
          !IF (GLOBALCOMM%NODE_ME.eq.1)  WRITE (*,*) 'dMP2 from element add     ', 0.5*SUM2
          DEALLOCATE(EV,E)
       ENDIF

    ELSE !dummy mode
       !exit blacs  
       CALL BLACS_GRIDEXIT( ICTXB )  
       CALL BLACS_GRIDEXIT( ICTXA )  

    ENDIF

    ! in case laplace transformed direct MP2 was done, the RPA energy does not make sense
    IF( LLTDMP2 ) COR%CORRELATION_K( NCUT ) = 0 

    PROFILING_STOP( 'col2bc_and_rpa' )
    RETURN
END SUBROUTINE COL2BC_AND_RPA


!****************************************************************************************
!
!This routine initializes a 1D process grid which is used for BLACS, i.e
!1 ROW and COMM_INOMEGA%NCPU COLS. 
!This is done for NGROUP groups by assinging an individual CONTEXT(I)
!for each group individually, (I=1,...,NGROUP)
!
!****************************************************************************************

SUBROUTINE GEN_PROC_GROUP_GRIDS(NGROUP, NODE_CRITICAL,&
                                COMM_BETWEENOMEGA,&
                                CONTEXT)
   IMPLICIT NONE
   INTEGER                :: NGROUP,NODE_CRITICAL
   TYPE(communic)         :: COMM_BETWEENOMEGA
   TYPE(communic)         :: COMM_INOMEGA
   INTEGER                :: CONTEXT        ! upon entry the MPI id
                                            ! upon exit the generated BLACS id
   !local
   INTEGER, ALLOCATABLE   :: USERMAP(:,:)
   INTEGER                :: I,J
    !----------------------------------------------------------------------
    !BLACS, get info about processors
    ! CALL BLACS_PINFO(ME,PROCS)

    CALL vtutor%bug("internal error in VASP: GEN_PROC_GROUP_GRIDS broken for some BLACS implementations", __FILE__, __LINE__)

    CALL MPI_barrier(CONTEXT, I )  ! defensive barrier
    CALL MPI_COMM_SIZE(CONTEXT, PROCS, I)  ! get number of processors

    !we need to group processors into subgrids, where
    !each subgrid contains CHI at a ceartain frequency point 
    !this is best done with a grid map of the usual node id to the sub grids
    ALLOCATE(USERMAP( 1 , NODE_CRITICAL  ))
    USERMAP=0

    !----------------------------------------------------------------------
    !loop over all groups
    DO I = 1, NGROUP
       !a 1D processor grid is not optimal for large matrices
       !probably a 2D block cyclic distribution is better, however I don't care
       NPROW0 = 1                   
       NPCOL0 = NODE_CRITICAL

       IF ( COMM_BETWEENOMEGA%NODE_ME == I ) THEN            !define the mapping for
          DO J = 1, NODE_CRITICAL                            !the Ith group using this
             USERMAP(1,J)= J-1 + (I-1)*NODE_CRITICAL         !mapping
          ENDDO

         IF ( COMM_INOMEGA%NODE_ME <= NODE_CRITICAL ) THEN
!         CALL BLACS_GET(0, 0, CONTEXT)
          CALL BLACS_GRIDMAP( CONTEXT, USERMAP, 1, NPROW0, NPCOL0 )  !map to the sub grids
          CALL BLACS_GRIDINFO( CONTEXT, NPROW0, NPCOL0, MYROW0, MYCOL0) !get the ids for scaLAPACK
         ENDIF
#if debug
       WRITE(*,'(" CPU ",I2," MYROW0,MYCOL0",2I2, "NPROW0,NPCOL0",2I2)')ME,MYROW0,MYCOL0,&
                  NPROW0, NPCOL0
#endif
      ENDIF

    ENDDO
    !----------------------------------------------------------------------
   DEALLOCATE(USERMAP)
   RETURN
END SUBROUTINE GEN_PROC_GROUP_GRIDS

!****************************************************************************************
!
!This routine initializes a 2D process grid which is used for BLACS
!the grid is set in the same way as in BLOCKSIZE_AND_PROC_GRID.
!This is done for NGROUP groups by assinging an individual CONTEXT(I)
!for each group individually, (I=1,...,NGROUP)
!
!****************************************************************************************

SUBROUTINE GEN_PROC_GROUP_GRIDS_BC(NGROUP, &
                                COMM_BETWEENOMEGA, COMM_INOMEGA,&
                                CONTEXT)
   USE scala, ONLY: FERMAT_RAZOR
   IMPLICIT NONE
   INTEGER                :: NGROUP
   TYPE(communic)         :: COMM_BETWEENOMEGA
   TYPE(communic)         :: COMM_INOMEGA
   INTEGER                :: CONTEXT        ! upon entry the MPI id
                                            ! upon exit the generated BLACS id
   !local
   INTEGER, ALLOCATABLE   :: USERMAP(:,:)
   INTEGER                :: I,J,IG, NCPU_GROUP
    !----------------------------------------------------------------------
    !BLACS, get info about processors
    ! CALL BLACS_PINFO(ME,PROCS)

    CALL vtutor%bug("internal error in VASP: GEN_PROC_GROUP_GRIDS broken for some BLACS implementations", __FILE__, __LINE__)

    CALL MPI_barrier(CONTEXT, I )  ! defensive barrier
    CALL MPI_COMM_SIZE(CONTEXT, PROCS, I)  ! get number of processors

    NCPU_GROUP = COMM_INOMEGA%NCPU 
  
    CALL FERMAT_RAZOR( NCPU_GROUP, NPROW, NPCOL ) 

    !we need to group processors into subgrids, where
    !each subgrid contains CHI at a ceartain frequency point 
    !this is best done with a grid map of the usual node id to the sub grids
    ALLOCATE(USERMAP( NPROW , NPCOL )) 
 
    !----------------------------------------------------------------------
    !loop over all groups
    DO IG = 1, NGROUP
      IF ( COMM_BETWEENOMEGA%NODE_ME == IG ) THEN               !define the mapping for

       DO I = 1, NPROW
        DO J = 1, NPCOL                                       !the IGth group uses this
          USERMAP(I,J)= J-1 + (I-1)*NPCOL + (IG-1)*NCPU_GROUP !mapping
        ENDDO
       ENDDO
    
!      CALL BLACS_GET(0, 0, CONTEXT)     
       CALL BLACS_GRIDMAP( CONTEXT, USERMAP, NPROW, NPROW, NPCOL )  !map to the sub grids
       CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL) !get the ids for scaLAPACK
 
#ifdef debug   
       WRITE(*,'(" CPU ",I2," MYROW,MYCOL",2I2, " NPROW,NPCOL",2I2)')ME,MYROW,MYCOL,&
                  NPROW, NPCOL
#endif
      ENDIF

    ENDDO 
    !----------------------------------------------------------------------
   DEALLOCATE(USERMAP)

   RETURN
END SUBROUTINE GEN_PROC_GROUP_GRIDS_BC



!****************************************************************************************
!
! determines block size for scaLAPACK and the processor grid
!
!****************************************************************************************

SUBROUTINE BLOCKSIZE_AND_PROC_GRID( NROWS, NCOLS, NCPU_GROUP, INU )
    USE scala, ONLY: FERMAT_RAZOR 
    IMPLICIT NONE
    INTEGER                 :: NROWS,NCOLS
    INTEGER                 :: NCPU_GROUP
    INTEGER                 :: INU

    CALL FERMAT_RAZOR( NCPU_GROUP, NPROW, NPCOL )     

    !write which processor grid is used
    IF ( INU >=0 ) THEN
        WRITE(*,'(" Using a ",I3,"x",I3," processor grid for a ",I5," x ",I3,"*",I5," matrix")')&
           NPROW, NPCOL, NROWS, NCPU_GROUP,NCOLS
    ENDIF

    !block size is minimum of blocking factor of column major distribution and MB
    MB = MIN( NCOLS , MB ) 
    NB = MIN( NCOLS , NB ) 

    IF ( INU >= 0 ) THEN
       WRITE(*,'(" distribution is done in blocks of ",I3," x ", I3)')MB, NB
    ENDIF

   
    RETURN

END SUBROUTINE BLOCKSIZE_AND_PROC_GRID

!******************************* LIN_REG_GG**********************************************
!
! perform linear regression of correlation energy versus
! 1/energy^(3/2)
!
!***************************************************************************************

  SUBROUTINE LIN_REG_GG(COR, E, INOUT)
  USE base 
  TYPE (correlation), POINTER :: COR
  TYPE( energy ) :: E
  INTEGER :: INOUT
 
  REAL(q) :: SX, SY, SXY, SX2, SXM, SYM, SXYM, SX2M
  REAL(q) :: AREG, BREG, AREGM, BREGM
  INTEGER :: I, ILAMBDA
  REAL(q) :: LAMBDA

    SX = 0.0   
    SY = 0.0
    SXY = 0.0
    SX2 = 0.0
    SXM = 0.0   
    SYM = 0.0
    SXYM = 0.0
    SX2M = 0.0
  
    AREG = 0.0
    BREG = 0.0
    AREGM = 0.0
    BREGM = 0.0

    IF ( LDMP1 ) THEN
    IF (INOUT>=0) WRITE(INOUT,'("      cutoff energy     smooth cutoff   RPA   correlation   Hartree contr. to MP1")')
    ELSE
    IF (INOUT>=0) WRITE(INOUT,'("      cutoff energy     smooth cutoff   RPA   correlation   Hartree contr. to MP2")')
    ENDIF

    IF (INOUT>=0) WRITE(INOUT,'("---------------------------------------------------------------------------------")')       


    ! set correlation energy with highest cutoff 
    IF ( L2ORDER ) THEN
       E%ERPA_CUT=REAL(COR%CORRMP2DIR(1),q)
    ELSE
       E%ERPA_CUT=REAL(COR%CORRELATION(1),q)
    ENDIF

    DO I=1,COR%NE
      IF (INOUT>=0) WRITE(INOUT,'(" ",2F18.3,2F20.10)') COR%ENCUTGW(I),COR%ENCUTGWSOFT(I), &
                 REAL(COR%CORRELATION(I),q), REAL(COR%CORRMP2DIR(I),q) 
      IF ( COR%NE==1 ) CYCLE   
      SX = SX+ 1/COR%ENCUTGW(I)**1.5
      SY = SY+ REAL(COR%CORRELATION(I),q)
      SXY = SXY+ 1/COR%ENCUTGW(I)**1.5*REAL(COR%CORRELATION(I),q)
      SX2 = SX2 + (1/COR%ENCUTGW(I)**1.5)**2
      SXM = SXM+ 1/COR%ENCUTGW(I)**1.5
      SYM = SYM+ REAL(COR%CORRMP2DIR(I),q)
      SXYM = SXYM+ 1/COR%ENCUTGW(I)**1.5*REAL(COR%CORRMP2DIR(I),q)
      SX2M = SX2M + (1/COR%ENCUTGW(I)**1.5)**2      
    ENDDO

    IF ( .NOT. COR%NE==1 ) THEN
        BREG =  (SY*SX - COR%NE*SXY)/(SX*SX - COR%NE*SX2) 
        AREG =  1/(COR%NE+0.0)*(SY - BREG*SX)
        BREGM = (SYM*SXM - COR%NE*SXYM)/(SXM*SXM - COR%NE*SX2M) 
        AREGM =  1/(COR%NE+0.0)*(SYM - BREGM*SXM) 

        IF (INOUT >=0) WRITE(INOUT,'("  linear regression    ")')  
        IF (INOUT >=0) WRITE(INOUT,'("  converged value                    ", 2F20.10)') AREG, AREGM
        ! store value from regression 
        E%ERPA_INF = AREG
        !   IF (INOUT >=0) WRITE(INOUT,'("  slope                              ", 2F20.10)') BREG, BREGM
    ELSE 
        E%ERPA_INF = REAL(COR%CORRELATION(1),q)
    ENDIF

    IF (NLAMBDA>=1) THEN
       IF (INOUT >=0) WRITE(INOUT,'("  converged value lambda= ",F7.4,"     ", 2F20.10)') 0.0_q,0.0_q
       DO ILAMBDA=0,NLAMBDA
       LAMBDA=1.0_q*(ILAMBDA+1)/MAX(1,NLAMBDA)
       SX = 0.0   
       SY = 0.0
       SXY = 0.0
       SX2 = 0.0

       DO I=1,COR%NE
          SX = SX+ 1/COR%ENCUTGW(I)**1.5
          SY = SY+ REAL(COR%CORRELATION_LAMBDA(I,ILAMBDA),q)
          SXY = SXY+ 1/COR%ENCUTGW(I)**1.5*REAL(COR%CORRELATION_LAMBDA(I,ILAMBDA),q)
          SX2 = SX2 + (1/COR%ENCUTGW(I)**1.5)**2
       ENDDO

       BREG =  (SY*SX - COR%NE*SXY)/(SX*SX - COR%NE*SX2) 
       AREG =  1/(COR%NE+0.0)*(SY - BREG*SX)

       IF (INOUT >=0) WRITE(INOUT,'("  converged value lambda= ",F7.4,"     ", 2F20.10)') LAMBDA,AREG
       ENDDO
    ENDIF


  END SUBROUTINE LIN_REG_GG

!***************************************************************************************
! copies HEAD, WING and CWING to RESPONSEFUN of CHI.
! this routine is the analogue of BODY_FROM_WING in chi_base.F 
! for a blockwise distributed CHI 
!***************************************************************************************

   SUBROUTINE BODY_FROM_WING_GG( CHI, IDIR, NOMEGA_CURRENT, NROWS, NCOLS, COMM_IN ) 
    IMPLICIT NONE
    TYPE (responsefunction) :: CHI            ! full response function
    INTEGER                 :: IDIR           ! current direction 
    INTEGER                 :: NOMEGA_CURRENT ! current frequency point in group
    TYPE (communic)         :: COMM_IN        ! communicator in group 
    INTEGER                 :: NROWS          ! number of rows in CHI%RESPONSEFUN
    INTEGER                 :: NCOLS          ! number of columns in CHI%RESPONSEFUN
 
    IF (.NOT.CHI%LGAMMA) THEN
       CALL vtutor%bug("internal error in BODY_FROM_WING_GG: LGAMMA is not set", __FILE__, __LINE__)
    ENDIF
     
    IF (CHI%LREALSTORE) THEN
       IF ( COMM_IN%NODE_ME*NCOLS > SIZE( CHI%CWINGR,1) )THEN 
          WRITE(*,*)' internal error in BODY_FROM_WING_GG: size of CWINGR too small',&
             COMM_IN%NODE_ME*NCOLS,  SIZE( CHI%CWINGR,1)
          RETURN
       ENDIF

       !cwings are distributed blockwise with blocking factor NCOLS
       CHI%RESPONSER(1,1:NCOLS,NOMEGA_CURRENT)=& 
          CHI%CWINGR(1+(COMM_IN%NODE_ME-1)*NCOLS:COMM_IN%NODE_ME*NCOLS,IDIR,NOMEGA_CURRENT)
        
       !only master node of group gets the head and wing
       IF ( COMM_IN%NODE_ME == 1 ) THEN
          CHI%RESPONSER(:,1,NOMEGA_CURRENT)=CHI%WINGR(:,IDIR,NOMEGA_CURRENT)
          CHI%RESPONSER(1,1,NOMEGA_CURRENT)=CHI%HEAD(IDIR,IDIR,NOMEGA_CURRENT)
       ENDIF

    ELSE
       IF ( COMM_IN%NODE_ME*NCOLS > SIZE( CHI%CWING,1) )THEN 
          WRITE(*,*)' internal error in BODY_FROM_WING_GG: size of CWING too small',&
             COMM_IN%NODE_ME*NCOLS,  SIZE( CHI%CWING,1)
          RETURN
       ENDIF 

       !cwings are distributed blockwise with blocking factor NCOLS
       CHI%RESPONSEFUN(1,1:NCOLS,NOMEGA_CURRENT)=& 
          CHI%CWING(1+(COMM_IN%NODE_ME-1)*NCOLS:COMM_IN%NODE_ME*NCOLS,IDIR,NOMEGA_CURRENT)
        
       !only master node of group gets the head and wing
       IF ( COMM_IN%NODE_ME == 1 ) THEN
          CHI%RESPONSEFUN(:,1,NOMEGA_CURRENT)=CHI%WING(:,IDIR,NOMEGA_CURRENT)
          CHI%RESPONSEFUN(1,1,NOMEGA_CURRENT)=CHI%HEAD(IDIR,IDIR,NOMEGA_CURRENT)
       ENDIF
    ENDIF
 

   END SUBROUTINE BODY_FROM_WING_GG

!************************************ DETERMINE_VXV ************************************
!
! this subroutine truncates the response function by multiplying from left and right
! with a truncated Coulomb kernel smoothly going from 1/G^2 to zero
! between ENCUTGWSOFT and ENCUTGW
!
!****************************************************************************************

  SUBROUTINE DETERMINE_VXV( CHI_WORK, CHI, NROWS, NCOLS, &
    COMM_INOMEGA, GLOBALCOMM, LATT_CUR, ENCUT, ENCUTSOFT, &
    NOMEGA, NOMEGA_SIMULTANEOUS, NCURR, WGWQ , MAXINDEX, &
    NODE_CRITICAL, LSR,LDD,INU)
    USE constant
    USE fock
    USE mpimy
    USE c2f_interface, ONLY : ERRF
    IMPLICIT NONE
    GDEF,POINTER  :: CHI_WORK(:,:)
    TYPE (responsefunction) :: CHI
    TYPE (latt) LATT_CUR
    INTEGER :: NCOLS, NROWS
    TYPE(communic) :: COMM_INOMEGA     !communicator in one group
    TYPE(communic) :: GLOBALCOMM     !global communicator
    REAL(q) :: ENCUT, ENCUTSOFT
    INTEGER :: NOMEGA
    INTEGER :: NOMEGA_SIMULTANEOUS   !# of frequencies in one group
    INTEGER :: NCURR                  !current freuquency point
    TYPE (wavedes1) :: WGWQ
    INTEGER :: MAXINDEX              ! number of grid points using truncated Coulomb kernel
    INTEGER :: NODE_CRITICAL         ! node id containing cutoff
    LOGICAL :: LSR
    INTEGER :: LDD,INU
  ! local
    INTEGER    NI, NP, NI_, N1, I, J, I_
    REAL(q) :: DKX, DKY, DKZ, GX, GY, GZ, GSQU, SCALE, POTFAK, E
    INTEGER :: NI_LOC, NOLD_LOC, NODE_RECV, IERROR
    REAL(q), ALLOCATABLE :: DATAKE(:)
    INTEGER, ALLOCATABLE :: OLD_INDEX(:)
    INTEGER :: ISEND(4) , IRECV(4)
    GDEF, ALLOCATABLE    :: CHI_RECV(:)
    GDEF, ALLOCATABLE    :: CHI_SEND(:)
    REAL(q) :: DFUN, SFUN
    ! neu
    REAL(q) :: OMEGBK,QC
    REAL(q) :: GAMMASCALE

    ! multipole correction
    INTEGER    NJ,NK
    COMPLEX(q) :: TSUM(10)
    COMPLEX(q) :: POTCORRECTION(10)
    !for communication
    INTEGER    :: NI_GLOB, NSEND
    INTEGER    :: NODE_SEND

    ! e^2/ volume 
    SCALE=EDEPS/LATT_CUR%OMEGA
!=======================================================================
! first set up the truncated Coulomb kernel
! smoothly going from 1/G^2 to zero between ENCUTGWSOFT and ENCUTGW
!=======================================================================
    DKX=(WGWQ%VKPT(1))*LATT_CUR%B(1,1)+ &
        (WGWQ%VKPT(2))*LATT_CUR%B(1,2)+ &
        (WGWQ%VKPT(3))*LATT_CUR%B(1,3)
    DKY=(WGWQ%VKPT(1))*LATT_CUR%B(2,1)+ &
        (WGWQ%VKPT(2))*LATT_CUR%B(2,2)+ &
        (WGWQ%VKPT(3))*LATT_CUR%B(2,3)
    DKZ=(WGWQ%VKPT(1))*LATT_CUR%B(3,1)+ &
        (WGWQ%VKPT(2))*LATT_CUR%B(3,2)+ &
        (WGWQ%VKPT(3))*LATT_CUR%B(3,3)

    NP=WGWQ%NGVECTOR
    IF (WGWQ%LGAMMA) NP=NP*2

    ALLOCATE( DATAKE(NROWS), OLD_INDEX(NROWS))
    DO NI = 1, NROWS
     OLD_INDEX(NI) = NI
    ENDDO

    !treated frequency point
    N1=NCURR

    DATAKEMAX=0.0
    MAXINDEX=0
    DO NI=1,NP
       NI_=NI
       IF (WGWQ%LGAMMA) NI_=(NI-1)/2+1
       
       GX=(WGWQ%IGX(NI_)*LATT_CUR%B(1,1)+WGWQ%IGY(NI_)* &
            LATT_CUR%B(1,2)+WGWQ%IGZ(NI_)*LATT_CUR%B(1,3))
       GY=(WGWQ%IGX(NI_)*LATT_CUR%B(2,1)+WGWQ%IGY(NI_)* &
            LATT_CUR%B(2,2)+WGWQ%IGZ(NI_)*LATT_CUR%B(2,3))
       GZ=(WGWQ%IGX(NI_)*LATT_CUR%B(3,1)+WGWQ%IGY(NI_)* &
            LATT_CUR%B(3,2)+WGWQ%IGZ(NI_)*LATT_CUR%B(3,3))
          
       GSQU=(DKX+GX)**2+(DKY+GY)**2+(DKZ+GZ)**2
       

       IF (ABS(GSQU)<G2ZERO) THEN
          ! head and wing
          POTFAK=SCALE
          ! if HFRCUT is set, a fixed spherical cutoff is used
          ! resulting in a finite values at G=0
          ! since the head stores the q^2, the final contribution is zero
          ! from the head
          IF (HFRCUT/=0) THEN
             POTFAK=0
          ! if LRHFCALC.AND.LRSCOR the Coulomb kernel is replaced by v_{LR}=v*exp(-q^2/(4*mu^2))
          ! for q=0  v_{LR}=1, this introduces an error for small mu because of the finite q-grid
          ! In order to correct this error, v_{LR}(q=0) is set to average of v_{LR} in a sphere  
          ! with volume corresponding to BZ volume/ # q-points       
          ELSE IF (LRHFCALC.AND.LRSCOR) THEN
             CALL CELVOL(LATT_CUR%B(1,1),LATT_CUR%B(1,2),LATT_CUR%B(1,3),OMEGBK)
             QC=TPI*(0.75_q/PI*OMEGBK/KPOINTS_FULL%NKPTS)**(1._q/3._q)
             GAMMASCALE=8*PI*HFSCREEN*HFSCREEN* &
            &   (-EXP(-QC*QC/4/HFSCREEN/HFSCREEN)*QC+HFSCREEN*SQRT(PI)*ERRF(QC/2/HFSCREEN))/ &
            &   (OMEGBK*TPI*TPI*TPI/KPOINTS_FULL%NKPTS)
             POTFAK=POTFAK*GAMMASCALE
          ! if LSR the Coulomb kernel is replaced by v_{SR}=v*[1-exp(-q^2/(4*mu^2))]
          ! for q=0  v_{SR}=0, this introduces an error for small mu because of the finite q-grid
          ! In order to correct this error, v_{SR}(q=0) is set to average of v_{SR} in a sphere
          ! with volume corresponding to BZ volume/ # q-points
          ELSE IF (LSR) THEN
             CALL CELVOL(LATT_CUR%B(1,1),LATT_CUR%B(1,2),LATT_CUR%B(1,3),OMEGBK)
             QC=TPI*(0.75_q/PI*OMEGBK/KPOINTS_FULL%NKPTS)**(1._q/3._q)
             GAMMASCALE=8*PI*HFSCREEN*HFSCREEN* &
            &   (-EXP(-QC*QC/4/HFSCREEN/HFSCREEN)*QC+HFSCREEN*SQRT(PI)*ERRF(QC/2/HFSCREEN))/ &
            &   (OMEGBK*TPI*TPI*TPI/KPOINTS_FULL%NKPTS)
             POTFAK=POTFAK*(1-GAMMASCALE)
          ENDIF

          ! switch off standard convergence correction
          IF(MCALPHA/=0) THEN
             POTFAK=0
          ENDIF

       ELSE
          ! the factor 1/(2 pi)^2 is required to obtain proper reciprocal
          ! lattice vector lenght
          POTFAK=SCALE/(GSQU*TPI**2)
          IF (HFRCUT/=0) THEN
             ! spherical cutoff on Coloumb kernel
             ! see for instance C.A. Rozzi, PRB 73, 205119 (2006)
!             POTFAK=POTFAK*(1-COS(SQRT(GSQU)*TPI*HFRCUT)*EXP(-(SQRT(GSQU)*TPI*HFRCUT)**2*HFRCUT_SMOOTH))
!test
             CALL DELSTP(3,SQRT(GSQU)*TPI*HFRCUT/10,DFUN,SFUN)
             POTFAK=POTFAK*(SFUN-0.5)*2
!test
          ELSE IF (LSR) THEN
             ! use the short range part of the Coulomb kernel
             POTFAK=POTFAK*(1._q-EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN))))
          ELSE IF (LRHFCALC.AND.LRSCOR) THEN
             ! use the long range part of the Coulomb kernel
             POTFAK=POTFAK*EXP(-GSQU*(TPI*TPI/(4*HFSCREEN*HFSCREEN)))
          ENDIF
       ENDIF

       ! smooth cutoff function between  ENCUTSOFT and ENCUT
       E=HSQDTM*(GSQU*TPI**2)
       IF (ENCUT>=0 .AND. E>ENCUT) THEN
          POTFAK=0
       ELSE

          MAXINDEX=MAXINDEX+1
          OLD_INDEX(MAXINDEX)=NI
          
          IF (ENCUTSOFT>=0 .AND. E>ENCUTSOFT) THEN
             POTFAK=POTFAK*(1+COS((E-ENCUTSOFT)/(ENCUT-ENCUTSOFT)*PI))/2
          ENDIF
       ENDIF

! POTFAK in fock has an additional factor. Since the corrections were originally devised
! for the fock routine this ensures consistency
       IF(MCALPHA/=0) THEN
          POTFAK=POTFAK*(1.0_Q/WGWQ%GRID%NPLWV)
       ENDIF

       DATAKE(NI)=SQRT(POTFAK)
       ! maximum kinetic energy
       IF(DATAKE(NI)>DATAKEMAX) THEN
          DATAKEMAX=DATAKE(NI)
       ENDIF

    ENDDO

    ! setup of multipole corrections
    IF (MCALPHA/=0) THEN
       CALL FOCK_MULTIPOLE_CORR_SETUP(LATT_CUR, GRIDHF)
       CALL GW_MULTIPOLE_PROJ_SETUP( WGWQ, LATT_CUR, DATAKE,DATAKEMAX )
    ENDIF

!=======================================================================
! now multiply the response function with the kernel
! from left and right hand side
!=======================================================================

    !save chi to chi_work
    IF (CHI%LREALSTORE) THEN
          CHI_WORK(1:NROWS,1:NCOLS) = CHI%RESPONSER(   1 : NROWS , 1 : NCOLS , N1 )
    ELSE
          CHI_WORK(1:NROWS,1:NCOLS) = CHI%RESPONSEFUN( 1 : NROWS , 1 : NCOLS , N1 )
    ENDIF

!#ifdef debug
!dump CHI%RESPONSER here for first cutoff only!
IF ( LDD == 1 ) THEN
   WRITE(200+GLOBALCOMM%NODE_ME,*) ' X for omega point ',N1 
   DO NI=1,MIN(NROWS,16)
      WRITE(200+GLOBALCOMM%NODE_ME,'(16F12.6)')CHI_WORK(NI,1:MIN(NCOLS,16))
   ENDDO
IF ( GLOBALCOMM%NODE_ME==1) THEN
   WRITE(99,*) ' V for omega point ',N1
ENDIF
ENDIF
!#endif

    !kill elements larger than NP
    IF ( NP < NROWS ) THEN
       DO NI = NP + 1, NROWS
          DATAKE(NI) = 0
       ENDDO
    ENDIF

    !DATAKE(:) has the coulomb potenial, this needs to be multiplied with CHI_WORK
    !and symmetrized
    !from the left
    DO NI=1,NROWS  
       !DATAKE contains the square root of V!
       POTFAK=DATAKE( NI )*DATAKE( NI )
!#ifdef debug
!dump diagonal of V
IF ( LDD == 1 ) THEN
IF ( GLOBALCOMM%NODE_ME==1) THEN
   WRITE(99,'(2F20.10)')POTFAK
ENDIF 
ENDIF
!#endif 
       CHI_WORK(NI,1:NCOLS)=CHI_WORK(NI,1:NCOLS)*POTFAK
    ENDDO
    !from right
    DO NI=1,NCOLS
       !DATAKE contains the square root of V!
       POTFAK=DATAKE( NI + (COMM_INOMEGA%NODE_ME-1)*NCOLS )**2
       CHI_WORK(1:NROWS,NI)=CHI_WORK(1:NROWS,NI)*POTFAK
    ENDDO


    ! multipole corrections:
    IF(MCALPHA/=0) THEN 
      CALL vtutor%error("Sorry MCALPHA not implemented")
    ENDIF
   
    !smooth energy cutoff done by replacing the NI^th coulmn by the OLD_INDEX(NI)^th column 
    !and analogously for the row. 
#ifdef debug
    IF (INU>=0 ) WRITE(*,'(" Shrinking started for cut ",2I4)')LDD
#endif 

    !The replacement of the rows is straight forward 
    DO NI=1,MAXINDEX
       CHI_WORK(NI,1:NCOLS)=CHI_WORK(OLD_INDEX(NI),1:NCOLS)
    ENDDO

    !for the columns this is not straight forward
    !since the whole matrix is distributed in stripes of NCOLS
    !among the nodes in each frequency group.

    !send columns between nodes
    ALLOCATE(CHI_SEND( NROWS ) )  
    ALLOCATE(CHI_RECV( NROWS ) ) 
    CHI_RECV = 0 
    CHI_SEND = 0 

    !determine which node in group has the critical stripe.
    DO NI = 1, COMM_INOMEGA%NCPU 
       IF (MAXINDEX > NI*NCOLS) CYCLE
       IF (((NI-1)*NCOLS < MAXINDEX).AND.(MAXINDEX <= NI*NCOLS)) THEN
          NODE_CRITICAL = NI 
          EXIT
       ENDIF
    ENDDO

#ifdef debug
 WRITE(100+WGWQ%COMM%NODE_ME,'("CUTOFF #",I4," NODE",I2," CRIT.NODE",I2," NROWS",I4," NP_NEW",I4)')&
      LDD, WGWQ%COMM%NODE_ME,NODE_CRITICAL, NROWS,MAXINDEX 
#endif
    CALL MPI_BARRIER(COMM_INOMEGA%MPI_COMM,IERROR)   !each node should start at the same time from here
    !==============================================================================================
    !we loop over the cpus in one group and let other node send the appropriate columns
    DO NODE_RECV = 1, COMM_INOMEGA%NCPU
    !==============================================================================================
        loc_col : DO NI = 1 , NCOLS                      !search column which should be replaced
        !------------------------------------------------------------------------------------------
 
           NI_GLOB = NI + (NODE_RECV-1)*NCOLS   !global index for current node NODE_RECV
  
 
           !find node NODE_SEND containing OLD_INDEX(NI_GLOB)
           DO NODE_SEND = 1, COMM_INOMEGA%NCPU
             IF ( 1 + ( NODE_SEND-1)*NCOLS <= OLD_INDEX(NI_GLOB) .AND. &
                  OLD_INDEX(NI_GLOB) <= NODE_SEND*NCOLS) EXIT
           ENDDO
 
           NSEND = OLD_INDEX(NI_GLOB) - (NODE_SEND-1)*NCOLS   !local column index on sending node
 
           IF ( NSEND < 1 ) THEN 
              CALL vtutor%error("Cannot send column with negative index " // str(NODE_CRITICAL) // "&
                 & " // str(NODE_RECV) // " " // str(NODE_SEND) // " " // str(NI_GLOB) // " " // &
                 str(OLD_INDEX(NI_GLOB)) // " " // str(NSEND) // " " // str(NCOLS))
           ENDIF
 
           
           !if old column is on the same node simply replace it 
           IF ( NODE_SEND == COMM_INOMEGA%NODE_ME .AND. NODE_RECV == COMM_INOMEGA%NODE_ME ) THEN 
#if debug  
              WRITE(100+WGWQ%COMM%NODE_ME,&
                  '( " COL ON SAME CPU #",I4," COL #",2I4, "REPLACED BY COL #", I4,10F10.5)')&
                  NODE_RECV, NI,NI_GLOB, NSEND,CHI_WORK(1:10,NSEND)
#endif
              CHI_WORK(1:NROWS, NI) = CHI_WORK(1:NROWS, NSEND)

           !otherwise tell NODE_SEND to send the appropriate column NSEND to NODE_RECV
           ELSEIF ( NODE_RECV < NODE_SEND ) THEN
              IF ( NODE_SEND == COMM_INOMEGA%NODE_ME ) THEN
#if debug  
                 WRITE(100+WGWQ%COMM%NODE_ME,'( " COL #",I4, "SENT TO COL", I4 , " ON CPU #", I4)')&
                 NSEND,NI,NODE_RECV
#endif
                 CHI_SEND(1:NROWS) = CHI_WORK(1:NROWS, NSEND)

                 !send this column to node NODE_RECV              
                 IF (CHI%LREALSTORE) THEN
                    CALLMPI ( M_send_d( COMM_INOMEGA, NODE_RECV, CHI_SEND, NROWS) )
                 ELSE
                    CALLMPI ( M_send_z( COMM_INOMEGA, NODE_RECV, CHI_SEND, NROWS) )
                 ENDIF

              ELSEIF ( NODE_RECV == COMM_INOMEGA%NODE_ME ) THEN 
#if debug  
                 WRITE(100+WGWQ%COMM%NODE_ME,'( " RECV COL #", I4 ," FROM CPU #"'//&
                    ',I4, " FOR COL", I4)') NSEND, NODE_SEND, NI
#endif
                 IF (CHI%LREALSTORE) THEN
                    CALLMPI ( M_recv_d( COMM_INOMEGA, NODE_SEND, CHI_RECV, NROWS) )
                 ELSE
                    CALLMPI ( M_recv_z( COMM_INOMEGA, NODE_SEND, CHI_RECV, NROWS) )
                 ENDIF
                     
                 !and replace CHI_RECV 
                 CHI_WORK(1:NROWS, NI ) = CHI_RECV(1:NROWS)
              ENDIF
           ENDIF

        !mpi barrier here 
        CALL MPI_BARRIER(COMM_INOMEGA%MPI_COMM,IERROR)   
        !------------------------------------------------------------------------------------------
        ENDDO loc_col                                !local columns  
    !==============================================================================================
    ENDDO                                    !cpus in group
    !==============================================================================================
    DEALLOCATE( CHI_SEND )
    DEALLOCATE( CHI_RECV )

    DEALLOCATE(DATAKE, OLD_INDEX)

    !save the result V.X.V to CHI
    !save chi to chi_work
    IF (CHI%LREALSTORE) THEN
       CHI%RESPONSER(   1 : NROWS , 1 : NCOLS , N1 ) = CHI_WORK(1:NROWS,1:NCOLS)
    ELSE
       CHI%RESPONSEFUN( 1 : NROWS , 1 : NCOLS , N1 ) = CHI_WORK(1:NROWS,1:NCOLS) 
    ENDIF

#ifdef debug
!dump CHI%RESPONSER here for first cutoff only!
IF ( LDD == 1 ) THEN
   WRITE(300+GLOBALCOMM%NODE_ME,*) ' X for omega point ',N1 
   DO NI=1,MIN(NROWS,16)
      WRITE(300+GLOBALCOMM%NODE_ME,'(16F12.6)')CHI_WORK(NI,1:MIN(NCOLS,16))
   ENDDO
ENDIF
#endif 
    RETURN
 END SUBROUTINE DETERMINE_VXV

!************************************ STORE_VXV_TO_CHI **********************************
!
! stores the second order screened potential W^(2) = V.X.V to CHI 
! DEPRECATED
!
!****************************************************************************************
SUBROUTINE STORE_VXV_TO_CHI( WDES, CHI, NROWS, NCOLS, WGWQ, &
    COMM_BETWEENOMEGA, COMM_INOMEGA, LATT_CUR, NOMEGA, &
    NOMEGA_START,NOMEGA_SIMULTANEOUS, COR, NQ , IDIR_MAX, NCUTOFF, IO)
    USE base
    USE constant
    USE mpimy
    USE ini
    IMPLICIT NONE
    TYPE (wavedes)          :: WDES
    TYPE (responsefunction) :: CHI                  !response function
    INTEGER                 :: NROWS, NCOLS         !number of local rows and cols of CHI
    TYPE (wavedes1)         :: WGWQ
    TYPE (communic)         :: COMM_BETWEENOMEGA    !communicator between tau
    TYPE (communic)         :: COMM_INOMEGA         !communicator in tau
    TYPE (latt)             :: LATT_CUR
    INTEGER                 :: NOMEGA_START         !number of frequency points
    INTEGER                 :: NOMEGA               !number of frequency points
    INTEGER                 :: NOMEGA_SIMULTANEOUS  !# of freq. points on one core
    TYPE (correlation)      :: COR                  !correlation energy
    INTEGER                 :: NQ                   !current q-point
    INTEGER                 :: IDIR_MAX             !maximum number of directions
    INTEGER                 :: NCUTOFF              !cutoff index for stoarge
    TYPE (in_struct)        :: IO                 
  ! local
    GDEF, POINTER, DIMENSION(:,:) :: CHI_WORK
    INTEGER                 :: NROWS_NEW, NCOLS_NEW !new size of redistributed chi 
    INTEGER                 :: NBLOCK                 !block size
    INTEGER                 :: NODE_CRITICAL          !node containing cut off
    INTEGER ::  IDIR, NP_NEW, I,IP, NK, NCURR
    INTEGER :: NOMEGA_GLOBAL
    INTEGER :: NOMEGA_IN_ROOT_GROUP                   ! # of frequencies in root group
    INTEGER :: NOMEGA_CURRENT
    LOGICAL :: LDMODE
    
    IF (CHI%LREAL .AND. .NOT. CHI%LREALSTORE) THEN
       CALL vtutor%bug("internal error in STORE_XVX_TO_CHI: \n or the Gamma point version " &
          // "CHI%LREALSTORE must be set", __FILE__, __LINE__)
    ENDIF

    IF (NCUTOFF > NE .OR. NCUTOFF <= 0 ) THEN
       CALL vtutor%error("cutoff index negative or too large! " // str(NE) // " " // str(NCUTOFF))
    ENDIF

    !each group needs to loop over the same # of frequencies internally,
    !otherwise when initializing BLACS in LNTRACE_OF_CHI the routine 
    !BLACS_PINFO in GEN_PROC_GROUP_GRIDS_BC waits for the CPUs in other group
    !to call itself. 
    !since the root group contains always the maximum # of frequencies we 
    !use NOMEGA_SIMULTANEOUS from the root node for the internal # of frequencies for all groups
    NOMEGA_IN_ROOT_GROUP = 0
    IF ( WDES%COMM_KIN%NODE_ME == 1 ) NOMEGA_IN_ROOT_GROUP = NOMEGA_SIMULTANEOUS
    CALLMPI( M_bcast_i( WDES%COMM_KIN, NOMEGA_IN_ROOT_GROUP, 1 ) ) 

    !scaLAPACK block size and processor grid is set here
    IF ( CHI%NQ == 1 ) THEN
       CALL BLOCKSIZE_AND_PROC_GRID( NROWS, NCOLS, COMM_INOMEGA%NCPU, IO%IU0 )
    ENDIF 
 
    omega_sim: DO NCURR = 1, NOMEGA_IN_ROOT_GROUP  !loop over frequencies within group
    !--------------------------------------------------------------------------------------
       IF ( NCURR <= NOMEGA_SIMULTANEOUS ) THEN
          !this is the global omega point treated currently in each group 
          NOMEGA_GLOBAL = NOMEGA_START + NCURR - 1 
          NOMEGA_CURRENT = NCURR
          LDMODE = .FALSE.
       ELSE
          !group having less points than the root group use their largest local frequency
          NOMEGA_GLOBAL = NOMEGA_START + NOMEGA_SIMULTANEOUS - 1 
          NOMEGA_CURRENT = NOMEGA_SIMULTANEOUS
          LDMODE = .TRUE. 
       ENDIF

       !-----------------------------------------------------------------------------------
       ! directions only need to be considered for the head 
       IF ( CHI%LGAMMA ) THEN
       directions: DO IDIR = 1, IDIR_MAX          !loop over the tree directions  
       !-----------------------------------------------------------------------------------
          !for the long-wave limit we copy head and wings to RESPONSEFUN
          CALL BODY_FROM_WING_GG( CHI, IDIR, NOMEGA_CURRENT, NROWS, NCOLS, COMM_INOMEGA)

          !CHI is distributed in groups to all processors
          !each group contains CHI at one frequency point. 
          !multiplication with the Coulomb potential needs to be done carefully,
          !CHI will be used for GW calculations, therefore
          !we copy the matrix to CHI_WORK and calculate the RPA correlation energy
          !using CHI_WORK 

          !first we allocate CHI_WORK in the same way as CHI
          ALLOCATE(CHI_WORK(NROWS,NCOLS)) 

!test: compute V.X.V and store the result to CHI
          CHI_WORK = 0
          CALL DETERMINE_VXV( CHI_WORK, CHI, NROWS, NCOLS,                 &
             COMM_INOMEGA, WDES%COMM_KIN, LATT_CUR, COR%ENCUTGW(NCUTOFF), &
             COR%ENCUTGWSOFT(NCUTOFF),  &
             NOMEGA, NOMEGA_SIMULTANEOUS, NOMEGA_CURRENT, WGWQ , NP_NEW, NODE_CRITICAL,      &
             .FALSE., NCUTOFF, IO%IU0 )
!test 
          DEALLOCATE ( CHI_WORK)     
       !-----------------------------------------------------------------------------------
       ENDDO directions  
       ELSE 
       !-----------------------------------------------------------------------------------
          !CHI is distributed in groups to all processors
          !each group contains CHI at one frequency point. 
          !multiplication with the Coulomb potential needs to be done carefully,
          !CHI will be used for GW calculations, therefore
          !we copy the matrix to CHI_WORK and calculate the RPA correlation energy
          !using CHI_WORK 
          
          !first we allocate CHI_WORK in the same way as CHI
          ALLOCATE(CHI_WORK(NROWS,NCOLS)) 
!test: compute V.X.V
          CHI_WORK = 0
          CALL DETERMINE_VXV( CHI_WORK, CHI, NROWS, NCOLS,                 &
             COMM_INOMEGA, WDES%COMM_KIN, LATT_CUR, COR%ENCUTGW(NCUTOFF), &
             COR%ENCUTGWSOFT(NCUTOFF),  &
             NOMEGA, NOMEGA_SIMULTANEOUS, NOMEGA_CURRENT, WGWQ , NP_NEW, NODE_CRITICAL,      &
             .FALSE., NCUTOFF, IO%IU0 )
!test
          DEALLOCATE ( CHI_WORK)     
       !-----------------------------------------------------------------------------------
       ENDIF
    !--------------------------------------------------------------------------------------
    ENDDO omega_sim

 RETURN
END SUBROUTINE STORE_VXV_TO_CHI

!****************************************************************************************
!
!> adds 1 to distributed matrix A
!
!****************************************************************************************
   FUNCTION DTRACE_AND_ADD( N, A, IA, JA, DESCA )
      REAL(q)      DTRACE_AND_ADD
      INTEGER      IA, JA, N
      INTEGER      DESCA( * )
      GDEF         A( * )
      INTEGER      BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_, &
                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER  ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1, &
                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER      ICURCOL, ICURROW, II, IOFFA, J, JB, JJ, JN, &
                   LDA, LL, MYCOL, MYROW, NPCOL, NPROW
      GDEF         TRACE
      EXTERNAL     BLACS_GRIDINFO, INFOG2L, ZGSUM2D, DGSUM2D
      INTEGER      ICEIL
      EXTERNAL     ICEIL
      INTRINSIC    MIN, MOD
      
      PROFILING_START( 'dtrace_and_add' ) 
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      
      TRACE = 0.0_q
      IF( N.EQ.0 ) THEN
         DTRACE_AND_ADD = TRACE
         RETURN
      END IF
      
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, II, JJ, ICURROW, ICURCOL )
      
      JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
      JB = JN-JA+1
      LDA = DESCA( LLD_ )
      IOFFA = II + ( JJ - 1 ) * LDA
      
      !    Handle first diagonal block separately
      
      IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
         DO LL = IOFFA, IOFFA + (JB-1)*(LDA+1), LDA+1
            TRACE = TRACE + A( LL )
             A(LL)=A(LL)+1.0_q
          ENDDO
      END IF
      
      IF( MYROW.EQ.ICURROW ) IOFFA = IOFFA + JB
      IF( MYCOL.EQ.ICURCOL ) IOFFA = IOFFA + JB*LDA
      ICURROW = MOD( ICURROW+1, NPROW )
      ICURCOL = MOD( ICURCOL+1, NPCOL )
      
      !    Loop over the remaining block of columns
      
      DO J = JN+1, JA+N-1, DESCA( NB_ )
         JB = MIN( JA+N-J, DESCA( NB_ ) )
      
         IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
            DO LL = IOFFA, IOFFA + (JB-1)*(LDA+1), LDA+1
               TRACE = TRACE + A( LL )
               A(LL)=A(LL)+1.0_q
            ENDDO
         END IF
      
         IF( MYROW.EQ.ICURROW ) IOFFA = IOFFA + JB
         IF( MYCOL.EQ.ICURCOL ) IOFFA = IOFFA + JB*LDA
         ICURROW = MOD( ICURROW+1, NPROW )
         ICURCOL = MOD( ICURCOL+1, NPCOL )
      ENDDO
#ifdef gammareal
      CALL DGSUM2D( DESCA( CTXT_ ), 'All', ' ', 1, 1, TRACE, 1, -1, 1 )
#else
      CALL ZGSUM2D( DESCA( CTXT_ ), 'All', ' ', 1, 1, TRACE, 1, -1, 1 )
#endif
      DTRACE_AND_ADD = TRACE
      
      PROFILING_STOP( 'dtrace_and_add' ) 
   END FUNCTION DTRACE_AND_ADD

!****************************************************************************************
!
!> calculates ln of distributed, cholesky decomposed  matrix A
!
!****************************************************************************************
   FUNCTION DLOG_DIAG( N, A, IA, JA, DESCA )
      REAL(q)      DLOG_DIAG
      INTEGER      IA, JA, N
      INTEGER      DESCA( * )
      GDEF         A( * )
      INTEGER      BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_, &
                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER  ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1, &
                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER      ICURCOL, ICURROW, II, IOFFA, J, JB, JJ, JN, &
                   LDA, LL, MYCOL, MYROW, NPCOL, NPROW
      GDEF         TRACE
      EXTERNAL     BLACS_GRIDINFO, INFOG2L, ZGSUM2D, DGSUM2D
      INTEGER      ICEIL
      EXTERNAL     ICEIL
      INTRINSIC    MIN, MOD
      
      PROFILING_START( 'dlog_diag' )   

      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
      
      TRACE = 0.0_q
      IF( N.EQ.0 ) THEN
         DLOG_DIAG = TRACE
         RETURN
      END IF
      
      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL, II, JJ, ICURROW, ICURCOL )
      
      JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
      JB = JN-JA+1
      LDA = DESCA( LLD_ )
      IOFFA = II + ( JJ - 1 ) * LDA
      
      !    Handle first diagonal block separately
      
      IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
         DO LL = IOFFA, IOFFA + (JB-1)*(LDA+1), LDA+1
            IF ( REAL(A(LL),q) > 10*EPSILON(1._q) .AND. GAIMAG(A(LL)) < 1.E-8_q ) THEN
               TRACE = TRACE + LOG( REAL(A( LL ),q) )
               !  A(LL)=A(LL)+1.0_q
            ELSE 
               WRITE(*,'("Internal error in VASP, eigenvalue #",I4," is negative or complex ",2F20.10)')&
                 LL, A( LL )
            ENDIF
         ENDDO
      END IF
      
      IF( MYROW.EQ.ICURROW ) IOFFA = IOFFA + JB
      IF( MYCOL.EQ.ICURCOL ) IOFFA = IOFFA + JB*LDA
      ICURROW = MOD( ICURROW+1, NPROW )
      ICURCOL = MOD( ICURCOL+1, NPCOL )
      
      !    Loop over the remaining block of columns
      
      DO J = JN+1, JA+N-1, DESCA( NB_ )
         JB = MIN( JA+N-J, DESCA( NB_ ) )
      
         IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
            DO LL = IOFFA, IOFFA + (JB-1)*(LDA+1), LDA+1
                !A(LL)=A(LL)+1.0_q
               IF ( REAL(A(LL),q) > 10*EPSILON(1._q) .AND. GAIMAG(A(LL)) < 1.E-8_q ) THEN
                  TRACE = TRACE + LOG( REAL(A( LL ),q) )
                  !  A(LL)=A(LL)+1.0_q
               ELSE 
                  WRITE(*,'("Internal error in VASP 2, eigenvalue #",I4," is negative or complex ",2F20.10)')&
                    LL, A( LL )
               ENDIF
            ENDDO
         END IF
      
         IF( MYROW.EQ.ICURROW ) IOFFA = IOFFA + JB
         IF( MYCOL.EQ.ICURCOL ) IOFFA = IOFFA + JB*LDA
         ICURROW = MOD( ICURROW+1, NPROW )
         ICURCOL = MOD( ICURCOL+1, NPCOL )
      ENDDO
      
#ifdef gammareal
      CALL DGSUM2D( DESCA( CTXT_ ), 'All', ' ', 1, 1, TRACE, 1, -1, 1 )
#else
      CALL ZGSUM2D( DESCA( CTXT_ ), 'All', ' ', 1, 1, TRACE, 1, -1, 1 )
#endif
      DLOG_DIAG = TRACE
      PROFILING_STOP( 'dlog_diag' )   
   END FUNCTION DLOG_DIAG

#endif

END MODULE ACFDT_GG
