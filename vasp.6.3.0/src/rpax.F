#include "symbol.inc"
MODULE rpax
      USE prec
      USE fock
      USE chi_base      
      USE lattice
      USE wpot
      USE OMP_LIB
      IMPLICIT NONE
!**********************************************************************
! The central quantity is:
! 
!  < v_k1,n1 v'_k2,n2 | W | c'_k4,n4 c_k3,n3 > =
!  int_d3r d3r' c_k1-q,n3(r')  v'*_k2,n2(r') c'_k2+q,n4(r)  v*_k1,n1(r)  W(r',r)
!
! scaLAPACK
! FTOD=Fourier-transformed overlap density
!**********************************************************************
#if defined (gammareal) || !defined (MPI) || !defined (scaLAPACK)
#else
      PRIVATE :: CALC_FTOD,REDISTRIBUTE_FTOD_GRID,INIT_BLACS_COLS,INIT_BLACS_GRID
#endif
      COMPLEX(q) , ALLOCATABLE, PUBLIC, SAVE :: FTOD_PW(:,:,:), FTOD_PW_AI(:,:,:,:,:,:)
      COMPLEX(q) , ALLOCATABLE, PUBLIC :: FTOD_PW_IA(:,:,:,:,:,:), FTOD_PW_IJ(:,:,:,:,:,:)
      COMPLEX(q) :: E_MP2, E_MP2_C
      COMPLEX(q) :: E_CCSD, E_CCSD_X
      GDEF      , ALLOCATABLE, PUBLIC :: FTOD_OC(:,:,:), FTOD_OC_AB(:,:,:,:,:,:), FTOD_OC_AI(:,:,:,:,:,:)
      GDEF , ALLOCATABLE, PUBLIC :: FTOD_OC_IA(:,:,:,:,:,:)
      GDEF , ALLOCATABLE, PUBLIC :: FTOD_OC_IJ(:,:,:,:,:,:)
      GDEF , ALLOCATABLE, PUBLIC ::  OC_IJ_TMP(:,:,:,:), PW_IJ_TMP(:,:,:,:)
      GDEF , ALLOCATABLE, PUBLIC ::  OC_IA_TMP(:,:,:,:), PW_IA_TMP(:,:,:,:), OC_IA_TMP2(:,:,:,:), PW_IA_TMP2(:,:,:,:)
      GDEF , ALLOCATABLE, PUBLIC ::  OC_AI_TMP(:,:,:,:), PW_AI_TMP(:,:,:,:)
      GDEF , ALLOCATABLE, PUBLIC :: TE4O(:,:), VVOO(:,:,:,:), OVOV(:,:,:,:), JLIK(:,:,:,:), &
         KLIJ(:,:,:,:), OOVV(:,:,:,:), VV(:,:), VOVO(:,:,:,:), &
         OOOV(:,:,:,:), VV2(:,:), C_NV(:,:), CHI_GAI(:,:,:,:,:,:), AIG(:,:,:), GAI_TMP(:,:,:,:)
      GDEFS , ALLOCATABLE, PUBLIC :: T2(:,:,:,:,:,:), T2_N_TMP(:,:,:,:,:), &
         OOVV_S(:,:,:,:), VVOO_S(:,:,:,:), &
         T2_TMP(:,:,:,:), VOVO2_S(:,:,:,:), &
         VOVO_S(:,:,:,:), VOVO3_S(:,:,:,:), &
         OOVV2_S(:,:,:,:), VV_S(:,:), VV2_S(:,:), &
         VOOV_S(:,:,:,:),  OOOV_S(:,:,:,:)     
      INTEGER, PUBLIC, SAVE :: NGVECTOR, NHVECTOR
      integer, ALLOCATABLE, PUBLIC, SAVE :: KPTS_MKPTS(:)
      integer, ALLOCATABLE, PUBLIC, SAVE :: PROCS_KPTS(:)
      integer, ALLOCATABLE, PUBLIC, SAVE :: MKPTS_KPTS(:)
      integer, PUBLIC, SAVE :: MY_NKPTS
      COMPLEX(qs) :: mix
      INTEGER, PUBLIC, SAVE  :: CONTXT_COLS,CONTXT_GRID
      INTEGER, PUBLIC, SAVE :: NPROW_GRID, ncc
      INTEGER, PUBLIC, SAVE :: VBMAX
      INTEGER, PUBLIC, SAVE :: PROCS, ME
      INTEGER, PRIVATE, SAVE :: NPROW,NPCOL,MYROW,MYCOL,NBLOCKAB
      INTEGER :: TETS, TET2S, TETSORT, iteration, NCLENGTH, NCLENGTH_N, GKQ
      !Scalapack array descriptor for RESPF__ before redistribution
      integer, dimension(9)   :: desc_FTOD_PW_br, desc_FTOD_OC_br
      !Scalapack array descriptor for RESPF__ after redistribution
      integer, dimension(9)   :: desc_FTOD_PW, desc_FTOD_OC
      !BLACS related variables and function
      !PROCS.. number of processors, ME... processor number, NPROW... number of rows in process grid
      !NPCOL... number of columns in process grid, myrow,mycol... my coordinates in process grid
      !mb... blocking size of rows for block cyclic distribution
      !nb... blocking size of columns for block cyclic distribution       
      
      REAL(q) :: NFLOAT4O,tetfpo,tetcom,tetr2d,tetvv,tetcc,tetvc,tetcontv
      INTEGER, EXTERNAL :: NUMROC, BLACS_PNUM
      TYPE(one_center_handle), POINTER, PUBLIC, SAVE :: H
      LOGICAL :: SINGLES,RING,LCCD,CANONICAL
      LOGICAL, PRIVATE :: MP2, T2READ
      CHARACTER(LEN=10) ::  DIR_APP
      INTEGER :: contr_ftod_time,time_fakic,time_con_akic
      
      CONTAINS

#if defined(gammareal) || !defined(MPI) || !defined (scaLAPACK)

      SUBROUTINE CALCULATE_RPAX(P,WDES,W,LATT_CUR,T_INFO,IO,KPOINTS,WGW,ENCUTGW, ENCUTGWSOFT, LMAXMP2,NBANDSGW)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE base
         USE mkpoints
         USE mpimy
         USE ini
         USE choleski
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(wavedes) WGW
         TYPE(wavespin) W        
         TYPE(type_info) T_INFO
         TYPE(potcar) P(T_INFO%NTYP)
         TYPE(latt) LATT_CUR
         INTEGER LMAXMP2,NBANDSGW
         REAL(q) :: ENCUTGW, ENCUTGWSOFT,change,denom
         TYPE (in_struct) IO
         TYPE (kpoints_struct) KPOINTS
         INTEGER :: NI,NJ,NR,NS,KI,KJ,KR,KS,KQ,KQ_,direction,KR_,NA
         integer :: TWOE4ORBITAL_COLS,TWOE4ORBITAL_ROWS,niteration,nwfup
         integer :: time_array1(8),time_array2(8),ems
         COMPLEX(q) :: E_TOT_CCSD, E_TOT_CCSD_X

         
      END SUBROUTINE CALCULATE_RPAX

#else

!***********************************************************************
!***********************************************************************

      SUBROUTINE CALCULATE_RPAX(P,WDES,W,LATT_CUR,T_INFO,IO,KPOINTS,WGW,ENCUTGW, ENCUTGWSOFT, LMAXMP2,NBANDSGW)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE base
         USE mkpoints
         USE mpimy
         USE ini
         USE choleski
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(wavedes) WGW
         TYPE(wavespin) W        
         TYPE(type_info) T_INFO
         TYPE(potcar) P(T_INFO%NTYP)
         TYPE(latt) LATT_CUR
         INTEGER LMAXMP2,NBANDSGW
         REAL(q) :: ENCUTGW, ENCUTGWSOFT,change,denom
         TYPE (in_struct) IO
         TYPE (kpoints_struct) KPOINTS
         INTEGER :: NI,NJ,NR,NS,KI,KJ,KR,KS,KQ,KQ_,direction,KR_,NA
         integer :: TWOE4ORBITAL_COLS,TWOE4ORBITAL_ROWS,niteration,nwfup
         integer :: time_array1(8),time_array2(8),ems
         COMPLEX(q) :: E_TOT_CCSD, E_TOT_CCSD_X

!#ifdef gammareal
!   ncc=1
!#else
         ncc=2
!#endif
!         IF (IO%IU0>=0) write(*,*)'ccsd setup:'
!         IF (IO%IU0>=0) write(*,*)'q:',q,'qs:',qs,'NBANDSGW=',NBANDSGW
         !CALL MPI_TYPE_CREATE(MPI_COMPLEX, 4, MPI_COMPLEX4)
         NBLOCKAB=14
         MP2=.FALSE.
         CANONICAL=.TRUE.
         SINGLES=.FALSE.
         LCCD=.FALSE.  !double-checked with older version using diamond till 3x3x3 kpoints
         RING=.TRUE.  !double-checked with ACFDT version using diamond till 3x3x3 kpoints
         IF (RING) SINGLES=.FALSE.
         IF (LCCD) SINGLES=.FALSE.
         
         CALL CHECK_FULL_KPOINTS ! all set up properly ?
!         CALL CHECK_SHIFTED_KPOINTS(WDES,W,KPOINTS)
!         IF ((IO%IU0>0) .and. (SHIFTED_KPOINTS)) write(IO%IU0,*)'You are using a shifted k-mesh.'
         !Initialize the 1D process grid
         CALL INIT_BLACS_COLS()
         CALL INIT_MY_KPOINTS(WDES,WGW)
         !Calculate the Response functions. calc_ftod also calls INIT_BLACS_GRID
         !IF (IO%IU0>=0)
         CALL CALC_FTOD(WDES,WGW,W,P,T_INFO,LATT_CUR,LMAXMP2,ENCUTGW,ENCUTGWSOFT,IO,FSG_STORE(1),0,NBANDSGW)

         !IF (NIT==1) CALL LOAD_T2(WDES,NBANDSGW)

         CALL SETUP_T2AMPLITUDES(WDES,NBANDSGW)

         !CALL CALC_FOCKMATRIX(WDES,W,WGW)

         !one contxt for every process
         !CALL INIT_BLACS_CONTXT_0()

         !FTOD_OC_IJ(:,:,:,:,:,:)
         CALL CALC_EFOCK(WDES,WGW,IO)


         !Allocate the TWOE4ORBITAL and TWOE4ORBITAL_X matrix
         
         IF (IO%IU0>=0) write(*,*)''
         IF (IO%IU0>=0) WRITE(IO%IU0,*)'Calculating T2AMPLITUDES'
         IF (IO%IU0>=0) write(*,*)''

         !IF (.not. T2READ) 
         T2=(0.0_qs,0.0_qs)
        ! T2_N=(0.0_qs,0.0_qs)
         !T1=(0.0_qs,0.0_qs)
         !T1_N=(0.0_qs,0.0_qs)
         E_TOT_CCSD=zero
         E_TOT_CCSD_X=zero

         !WRITE (DIR_APP  , "('T2',I1,I1)") MOD(WDES%COMM%NODE_ME/10,10),MOD(WDES%COMM%NODE_ME,10)
         !open(unit = 7,file = DIR_APP)
         !READ(7,*)T2

         CALL APPEND_CDER_TO_FTOD(WDES,W,LATT_CUR,1)
         

         DO GKQ=1,KPOINTS_ORIG%NKPTS !WDES%NKPTS 
         IF (IO%IU0>=0) WRITE(*,*)'Calculating RPA+SOSEX corr. energy from kpoint ',GKQ
         IF (IO%IU0>=0) WRITE(IO%IU6,*)'Calculating RPA+SOSEX corr. energy from kpoint ',GKQ
         T2=(0.0_qs,0.0_qs)
         DO ITERATION=1,12
            CALL START_TIMING("G")
!            write(*,*)'tets1'
            CALL CALC_T1_T2_N(WDES,WGW,W,ITERATION,NBANDSGW)
!            write(*,*)'tets2'
!            CALL UPDATE_T2(WDES,WGW,W,IO,NBANDSGW,ITERATION)
	    IF (IO%IU0>=0) WRITE(*,*)'iteration:',iteration
            CALL CALC_CCSD_ENERGY(WDES,WGW,W,IO,NBANDSGW)

            CALL STOP_TIMING("G",IO%IU6,"LCCD")
         enddo
	    CALL CALC_CCSD_ENERGY(WDES,WGW,W,IO,NBANDSGW)
            E_TOT_CCSD=E_TOT_CCSD+E_CCSD*KPOINTS_ORIG%WTKPT(GKQ)/WDES%WTKPT(GKQ)
            E_TOT_CCSD_X=E_TOT_CCSD_X+E_CCSD_X*KPOINTS_ORIG%WTKPT(GKQ)/WDES%WTKPT(GKQ)

!         CALL CALC_FINAL_T2(WDES,WGW,W,IO,ITERATION,NBANDSGW)

         ENDDO

         IF (IO%IU0>=0) WRITE(*,*)'Final RPA+SOSEX corr. energy ',E_TOT_CCSD
         IF (IO%IU0>=0) WRITE(*,*)'Final SOSEX energy ',-E_TOT_CCSD_X
 
         IF (IO%IU0>=0) WRITE(IO%IU6,*)'Final RPA+SOSEX corr. energy ',E_TOT_CCSD
         IF (IO%IU0>=0) WRITE(IO%IU6,*)'Final SOSEX energy ',-E_TOT_CCSD_X


         !IF (IO%IU0>=0) WRITE(IO%IU6,*)'E MP2=',E_MP2

         !CALL MAKE_BRUECKNER_WF(WDES,W, LATT_CUR, LOVERL,CQIJ,LMDIM)
         !CALL ORTHCH(WDES,W, LOVERL,LMDIM,CQIJ)


!        WRITE (DIR_APP  , "('T2',I1,I1)") MOD(WDES%COMM%NODE_ME/10,10),MOD(WDES%COMM%NODE_ME,10)
!        open(unit = 7,file = DIR_APP)
!        write(7,*)T2
         
!        CLOSE(7)

!         open(unit = 7,file = "EV_SHIFT")

         IF (IO%IU0>=0) THEN
!            DO KI=1,WDES%NKPTS
               !write(7,*)KI,":",WDES%VKPT(:,KI)
!            DO NI=1,VBMAX
!               write(7,*)KI,NI
!               write(7,*)K_KI(NI,NI,KI)
!            ENDDO
!            DO NI=1,(PROCS*WDES%NBANDS-VBMAX)
!               write(7,*)KI,NI+VBMAX
!               write(7,*)K_AC(NI,NI,KI)
!            ENDDO
!            ENDDO
         ENDIF
!         CLOSE(7)

!         open(unit = 7,file = "T1")

!         IF (IO%IU0>=0) THEN
!            DO KI=1,WDES%NKPTS
!               !write(7,*)KI,":",WDES%VKPT(:,KI)
!            DO NI=1,VBMAX
!               write(7,*)KI,NI,":"
!               !write(7,*)T1(:,NI,KI)
!            DO NA=1,(PROCS*WDES%NBANDS-VBMAX)
!               write(7,*)NA+VBMAX
!               write(7,*)T1(NA,NI,KI)
!            ENDDO
!            ENDDO
!            ENDDO
!         ENDIF

         call BLACS_GRIDEXIT(contxt_cols)
         call BLACS_GRIDEXIT(contxt_grid)
         
      END SUBROUTINE CALCULATE_RPAX


      !***********************************************************************
      !This routine calculates the fourier-transformed overlap integrals(FTOD) <i|-G|a>
      !and <j|G|b>*4*pi*e^2/(G+q)^2 and also calls the redistribution routine.
      !This is done for the plane-wave part, as well as for the one-center terms.
      !Note, that in the gamma-only version it is enough to save <i|-G|a>*SQRT(4*pi*e^2/(G+q)^2)
      !for the plane-wave part, because (<i|-G|a>)*=<i|G|a>.
      !However, this approximation cannot be made for the one-center terms in the gamma-only
      !version, because of practical reasons.
      !The <i|G|a> quantities are most important for the construction of the two electron
      !four orbital integrals <ij|ab>.  <ij|ab>=4pi e^2 \sum_G <i|-G|a> <j|G|b> /(G+q)^2
      !***********************************************************************

      SUBROUTINE CALC_FTOD(WDES,WGW,W,P,T_INFO,LATT_CUR,LMAXMP2,ENCUTGW,ENCUTGWSOFT,IO,FSG,IT,NBANDSGW)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE mkpoints
         USE lattice
         USE constant
         USE base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(wavedes) WGW
         TYPE(wavespin) W
         TYPE(type_info) T_INFO
         TYPE(potcar) P(T_INFO%NTYP)
         TYPE(latt) LATT_CUR
         INTEGER LMAXMP2
         TYPE (in_struct) IO
         ! local variables
         TYPE(wavespin) WHF
         TYPE(wavedes1) WGWQ
         TYPE(wavedes1) WDESKI,WDESKA,WDESKB
         TYPE(wavefun1), ALLOCATABLE :: WI(:), WA(:), WB(:)
         INTEGER KQ,KI,KA,KB,KI_IN_FULL_ORIG,KQ_,MKQ
         INTEGER NSTRIP, NSTRIPA, ISP, NFFT, IT,NBANDSGW
         INTEGER NBI,NBA,NBAA, i, nbi_start, nbi_end,rnba
         INTEGER NP ! number of plane waves for the overlap density i*(r)a(r) GCHGIA(r)
         REAL(q) :: FSG ! singularity correction
         REAL(q) :: ENCUTGW,ENCUTGWSOFT 
         REAL(q) :: POTFAK(GRIDHF%MPLWV) 
         COMPLEX(q) CPHASE(GRIDHF%MPLWV)
         COMPLEX(q) CPHASE2(GRIDHF%MPLWV)
         LOGICAL LPHASE  
         LOGICAL LPHASE2
         COMPLEX(q) :: GWORK(MAX(GRIDHF%MPLWV,WGW%GRID%MPLWV)) !work array for calculating GCHGIA
         COMPLEX(q), ALLOCATABLE :: GCHGIA(:,:,:)  ! charge
         GDEF      , ALLOCATABLE :: CRHOIA(:,:)    ! one-center charge
         GDEF      , ALLOCATABLE :: CRHOIB(:,:)    ! one-center charge
         GDEF      , ALLOCATABLE :: CRHOLM(:)      ! augmentation occupancy matrix
         COMPLEX(q), ALLOCATABLE :: tmp_FTOD_PW(:,:,:,:)
         GDEF      , ALLOCATABLE :: tmp_FTOD_OC(:,:,:,:)
         Real(qs) :: mem_req
         
         CALL CHECK_FULL_KPOINTS
         
         IF (LMAXMP2>=0) THEN
           CALL SET_UP_ONE_CENTER_H(WDES,P,T_INFO,LMAXMP2,H)
         ENDIF
         !FSG=zero
         WHF=W
         WHF%WDES => WDES_FOCK
         NSTRIP=30
         CALL SETWDES(WHF%WDES,WDESKI,0)
         CALL SETWDES(WHF%WDES,WDESKA,0)
         CALL SETWDES(WHF%WDES,WDESKB,0)
         
         VBMAX=LAST_FILLED_XI_NOMOD(W,1,1)
         
         DO ISP=2,WDES%ISPIN
            VBMAX=MAX(LAST_FILLED_XI_NOMOD(W,1,ISP),VBMAX)
         ENDDO
         
         ALLOCATE(WI(NSTRIP),WA(NSTRIP),WB(NSTRIP))
         DO NBI=1,NSTRIP !VBMAX
            CALL NEWWAV(WI(NBI),WDESKI,.TRUE.)
         ENDDO
         DO NBA=1,NSTRIP
            CALL NEWWAV(WA(NBA),WDESKA,.TRUE.)
         ENDDO
         DO NBA=1,NSTRIP
            CALL NEWWAV(WB(NBA),WDESKB,.TRUE.)
         ENDDO
         
         NGVECTOR=MAXVAL(WGW%NGVECTOR(:))
         NGVECTOR=NGVECTOR+1
         NHVECTOR=0
         IF (ASSOCIATED(H)) THEN
            NHVECTOR=H%TOTAL_ENTRIES
         ENDIF
         ! Make sure that NGVECTOR is suited for the BLACS process grid
         IF (MOD(NGVECTOR,NPROW_GRID)/=0) THEN
            NGVECTOR=NGVECTOR+(NPROW_GRID-MOD(NGVECTOR,NPROW_GRID))
         ENDIF
         ! Make sure that NHVECTOR is suited for the BLACS process grid
         IF (MOD(NHVECTOR,NPROW_GRID)/=0) THEN
            NHVECTOR=NHVECTOR+(NPROW_GRID-MOD(NHVECTOR,NPROW_GRID))
         ENDIF
         
         
         !Initialize the 2D process grid
         CALL INIT_BLACS_GRID(WDES)
         !write(*,*)'before setupf ftod'
         
         !STOP
         IF (IT==0) THEN
            IF (ALLOCATED(FTOD_PW)) THEN
               DEALLOCATE(FTOD_PW)
            ENDIF
            CALL SETUP_FTOD(WDES,NBANDSGW)
         ENDIF
         
         
         !write out an estimate for the required memory
         mem_req=2.0_q*(PROCS*PROCS*WDES%NBANDS*WDES%NBANDS)*16.0_q/1024.0_q/1024.0_q/1024.0_q !TWOE4ORBITAL arrays
         mem_req=mem_req+NGVECTOR*PROCS*WDES%NBANDS*WDES%NBANDS*WDES%NKPTS*WDES%NKPTS*2.0_q*16.0_q/1024.0_q/1024.0_q/1024.0_q 
         mem_req=mem_req+(WDES%NBANDS-VBMAX)*VBMAX*(WDES%NBANDS-VBMAX)*VBMAX*WDES%NKPTS*WDES%NKPTS*WDES%NKPTS*16.0_q/1024.0_q/1024.0_q/1024.0_q 
         IF (LMAXMP2>=0) mem_req=mem_req+NHVECTOR*PROCS*WDES%NBANDS*VBMAX*WDES%NKPTS*WDES%NKPTS*2.0_q*16.0_q/1024.0_q/1024.0_q/1024.0_q
         IF (IO%IU0>0) THEN
            WRITE(*,"(A,E10.3,A)") 'For LCCD calculations approximately',mem_req,'GB RAM will be required.' 
         ENDIF
         !Setup the descriptors for the distributed TWOE4ORBITAL matrix   
         !and allocate the TWOE4ORBITAL and TWOE4ORBITAL_X matrices
         
         ALLOCATE(tmp_FTOD_PW(NGVECTOR,WDES%NBANDS,NSTRIP,ncc)) !VBMAX
         IF (ASSOCIATED(H)) THEN
            ALLOCATE(tmp_FTOD_OC(NHVECTOR,WDES%NBANDS,NSTRIP,2)) !VBMAX
         ENDIF
         
         IF (IO%IU0>0) THEN
            WRITE(IO%IU0,*)
            WRITE(IO%IU0,*)'Calculating fourier transformed overlap densities:'
         ENDIF

         call BLACS_GRIDINFO(CONTXT_COLS, NPROW, NPCOL, MYROW, MYCOL)
         
         IF (ASSOCIATED(H)) THEN
            ALLOCATE(CRHOIA(NHVECTOR, NSTRIP),CRHOIB(NHVECTOR, NSTRIP))
            FTOD_OC=zero
         ENDIF
         
         spin: DO ISP=1,WDES%ISPIN
         kqloop: DO KQ=1,WDES%NKPTS
            !KQ=KPOINT_IN_FULL_GRID(-WDES%VKPT(:,MKQ),KPOINTS_FULL)
!            IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KQ)/=0.00_q)) CYCLE
            IF (IO%IU0>=0) WRITE(IO%IU0,*)
            IF (IO%IU0>=0) THEN
               IF (WDES%ISPIN==1) THEN
                 WRITE(IO%IU0,'("NQ=",I4,3F10.4,", ")') KQ,WDES%VKPT(:,KQ)
               ELSE
                 WRITE(IO%IU0,'("NQ=",I4,3F10.4,", ",A1,A1,", ")') KQ,WDES%VKPT(:,KQ),ISP
               ENDIF
            ENDIF
            CALL SETWDES(WGW,WGWQ,KQ)
         
            NP=WGWQ%NGVECTOR
            IF (NP>NGVECTOR) THEN
               WRITE(*,*)'Internal error in "Calc_2orbital_response": NP larger than NGVECTOR'
               EXIT
            ENDIF
            ALLOCATE(GCHGIA(NP,NSTRIP,2),CRHOLM(AUG_DES%NPRO*WDES%NRSPINORS))
            !write(*,*)'procs=',procs
            kiloop: DO KI=1,WDES%NKPTS
!               IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KI)==0)) CYCLE
               tmp_FTOD_PW=(0._q,0._q)
               IF (ASSOCIATED(H)) tmp_FTOD_OC=zero
               
               CALL GWPROGRESS(IO%IU0, KI,WDES%NKPTS,KQ,WDES%NKPTS)
               KI_IN_FULL_ORIG=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,KI),KPOINTS_FULL_ORIG)
               ! collect NSTRIP nbi bands at k_i
               CALL SETWDES(WHF%WDES,WDESKI,KI)
               
               DO NBI_start=1,PROCS*WDES%NBANDS,NSTRIP
                  call BLACS_GRIDINFO(CONTXT_COLS, NPROW, NPCOL, MYROW, MYCOL)
                  NBI_end=NBI_start+MIN(PROCS*WDES%NBANDS-NBI_start+1,NSTRIP)-1
                  
                  CALL W1_GATHER_GLB(WHF,NBI_start,NBI_end,ISP,WI)
                  ! k_b = k_i - k_q - G
                  KB=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KQ)+WDES%VKPT(:,KI),KPOINTS_FULL)
!                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KB)==0)) Write(*,*)'error in calc_ftod with shifted k-mesh'
                  ! k_a = k_i + k_q - G
                  KA=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KI)-WDES%VKPT(:,KQ),KPOINTS_FULL)
!                  IF ((SHIFTED_KPOINTS) .and. (WDES%WTKPT(KA)==0)) Write(*,*)'error in calc_ftod with shifted k-mesh'
               
                  CALL SETWDES(WHF%WDES,WDESKA,KA)
               
                  ! CPHASE(r) = e^iGr, where G = k_i - k_q - k_b
                  CALL SETPHASE(WDES%VKPT(:,KI)-WDES%VKPT(:,KQ)-WDES%VKPT(:,KA),GRIDHF,CPHASE,LPHASE)
               
                  CALL SET_GFAC_WITHOUT_WEIGHT(GRIDHF,LATT_CUR,KI,KA,FSG,POTFAK)
                  ! 1/(G+q)**2
               
                  IF (ENCUTGW /= ENCUTGWSOFT .AND. ENCUTGW > 0 &
                   .AND. ENCUTGWSOFT > 0) THEN
                     CALL SET_GFAC_WAVEFUN(WGWQ,LATT_CUR,FSG,POTFAK,ENCUTGW,ENCUTGWSOFT)
                  ELSE
                     CALL SET_GFAC_WAVEFUN(WGWQ,LATT_CUR,FSG,POTFAK)
                  ENDIF
               
#ifdef gammareal
                  POTFAK=SQRT(POTFAK)
#endif

                  ! loop over all bands
                  DO NBA=1,WDES%NBANDS,NSTRIP
                     NSTRIPA=MIN(WDES%NBANDS+1-NBA,NSTRIP)
                  ! FFT{psi_a} to real space
                  DO NBAA=1,NSTRIPA !copy and fourier transform NSTRIP wave functions
                     CALL W1_COPY( ELEMENT(WHF,WDESKA,NBA+NBAA-1,ISP),WA(NBAA))
                     CALL FFTWAV_W1(WA(NBAA))
                  ENDDO
                  ! loop over valence bands only
                  DO NBI=1,MIN(PROCS*WDES%NBANDS-NBI_start+1,NSTRIP)
                     GCHGIA=0
                     
                     IF (ASSOCIATED(H)) THEN
                        CRHOIA=zero
                        CRHOIB=zero
                        CRHOLM=zero
                     ENDIF
                     !loop over all bands in NSTRIP 
                     DO NBAA=1,NSTRIPA
		            CALL LOC2GLOB((NBA+NBAA-1),ME,WDES%NB_TOT,PROCS,1,rnba)
                            IF (((NBI+NBI_START-1)>VBMAX) .and. ((RNBA)> VBMAX)) CYCLE
                     ! calculate rho(r)=psi_i(r)* psi_a(r) for,
                     ! one center terms and, on the plane wave grid.
                        IF (ASSOCIATED(H)) THEN
                           CALL FOCK_CHARGE_ONE_CENTER_NOINT( WI(NBI),WA(NBAA),&
                             GWORK(1),H,CRHOIA(1,NBAA), CRHOLM,  SIZE(CRHOLM))
                        ELSE
                           CALL FOCK_CHARGE_NOINT( WI(NBI),WA(NBAA), GWORK(1), &
                             CRHOLM, SIZE(CRHOLM))
                        ENDIF
                        ! Set phase e^iGr, where G = k_i - k_q - k_b
                        IF (LPHASE) THEN
                           CALL APPLY_PHASE( GRIDHF, CPHASE(1), GWORK(1), &
                             GWORK(1))
                        ENDIF
                   
                     ! FFT{rho} to reciprocal space
                        CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), &
                          GWORK(1),GCHGIA(1,NBAA,1),WGWQ%GRID,.FALSE.)
                        NFFT=NFFT+1
                     ! multiply with potential factor
                     
                        CALL APPLY_GFAC_WAVEFUN(WGWQ,GCHGIA(1,NBAA,1), &
                         POTFAK(1))
                         
                     ENDDO !NBAA (loop over all bands in NSTRIP)
                     
                     IF (ASSOCIATED(H)) THEN
                        CALL APPLY_PHASE_ONE_CENTER(WHF%WDES, H, CRHOIA(:,:), &
                          WHF%WDES%VKPT(:,KI)-WHF%WDES%VKPT(:,KA))
                     ENDIF
                     
                     IF (ASSOCIATED(H)) THEN
                        ! use CRHOIA as temporary work array
                        CALL APPLY_ONE_CENTER_H( WHF%WDES, H, CRHOIA(:,:), CRHOIB(:,:), NSTRIPA)
                     ENDIF
                     
                     DO NBAA=1,NSTRIPA
		            CALL LOC2GLOB((NBA+NBAA-1),ME,WDES%NB_TOT,PROCS,1,rnba)
 		            IF (((NBI+NBI_START-1)>VBMAX) .and. ((RNBA)> VBMAX)) CYCLE
                        !copy response functions to RESPF_PW and RESPF_OC
#ifdef gammareal
                           tmp_FTOD_PW(1:NP,NBAA+NBA-1,NBI,1)=(GCHGIA(1:NP,NBAA,1))*SQRT(1.0_q/GRIDHF%NPLWV)
#else
                           tmp_FTOD_PW(1:NP,NBAA+NBA-1,NBI,1)=(GCHGIA(1:NP,NBAA,1))
#endif
                        IF (ASSOCIATED(H)) THEN
                           tmp_FTOD_OC(1:H%TOTAL_ENTRIES,NBAA+NBA-1,NBI,1)=(CRHOIA(1:H%TOTAL_ENTRIES,NBAA))
                        ENDIF

                     ENDDO
                     
!!!!!!!!!!!!!!!!! GAMMA-only version !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef gammareal
                                     
                     ! calculate rho(r)=psi_i(r)* psi_a(r) for,
                     ! (1._q,0._q) center terms and, on the plane wave grid.
                     IF (ASSOCIATED(H)) THEN
                        CRHOIB=0
                        CRHOLM=0
                        !loop over all bands in NSTRIP 
                        DO NBAA=1,NSTRIPA
			    CALL LOC2GLOB((NBA+NBAA-1),ME,WDES%NB_TOT,PROCS,1,rnba)
			   IF (((NBI+NBI_START-1)>VBMAX) .and.( (RNBA)> VBMAX)) CYCLE
                           CALL FOCK_CHARGE_ONE_CENTER_NOINT( WA(NBAA),WI(NBI),&
                             GWORK(1),H,CRHOIA(1,NBAA), CRHOLM,  SIZE(CRHOLM))                
                        ! FFT{rho} to reciprocal space
                           CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), &
                             GWORK(1),GCHGIA(1,NBAA,1),WGWQ%GRID,.FALSE.)
                           NFFT=NFFT+1
                        ! multiply with potential factor
                        ENDDO !NBAA (loop over all bands in NSTRIP)
                        DO NBAA=1,NSTRIPA
                          !copy response functions to RESPF_PW and RESPF_OC
                           IF (ASSOCIATED(H)) THEN
                              tmp_FTOD_OC(1:H%TOTAL_ENTRIES,NBAA+NBA-1,NBI,2)=(CRHOIA(1:H%TOTAL_ENTRIES,NBAA))    !000000000000000000000
                           ENDIF
                        ENDDO
                     ENDIF
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!GAMMA only version !!!!!!!!!!!!!!!!!!!!!!!!

                     ENDDO !NBI (loop over nstrip nbi bands only)
                  ENDDO !NBA (loop over all bands)
               
#ifdef gammareal
               
#else
               
                  CALL SETWDES(WHF%WDES,WDESKB,KB)        
               
                  ! k_i - k_a = k_q + G
                  CALL PHASER_HF(GRIDHF,LATT_CUR,FAST_AUG_FOCK,WDES%VKPT(:,KB)-WDES%VKPT(:,KI))
               
                  ! CPHASE(r) = e^iGr, where G = k_i - k_q - k_a
                  CALL SETPHASE(WDES%VKPT(:,KB)-WDES%VKPT(:,KQ)-WDES%VKPT(:,KI),GRIDHF,CPHASE,LPHASE)
               
                  ! loop over all bands
                  DO NBA=1,WDES%NBANDS,NSTRIP
                     NSTRIPA=MIN(WDES%NBANDS+1-NBA,NSTRIP)
                     ! FFT{psi_a} to real space
                     DO NBAA=1,NSTRIPA !copy and fourier transform NSTRIP wave functions
                        CALL W1_COPY( ELEMENT(WHF,WDESKB,NBA+NBAA-1,ISP),WB(NBAA))
                        CALL FFTWAV_W1(WB(NBAA))
                     ENDDO
                     ! loop over valence bands only
                     DO NBI=1,MIN(PROCS*WDES%NBANDS-NBI_start+1,NSTRIP)
                        GCHGIA=0
                     
                        IF (ASSOCIATED(H)) THEN
                           CRHOIB=0
                           CRHOLM=0
                        ENDIF
                        !loop over all bands in NSTRIP
                        DO NBAA=1,NSTRIPA
			   CALL LOC2GLOB((NBA+NBAA-1),ME,WDES%NB_TOT,PROCS,1,rnba)
			   IF (((NBI+NBI_START-1)>VBMAX) .and. ((RNBA)> VBMAX)) CYCLE
                        ! GCHGIA number 2 for X-changed waves
                        ! calculate rho(r)=psi_i(r)* psi_a(r) for,
                        ! one center terms and, on the plane wave grid.
                           IF (ASSOCIATED(H)) THEN
                              CALL FOCK_CHARGE_ONE_CENTER_NOINT( WB(NBAA),WI(NBI),&
                                GWORK(1),H,CRHOIB(1,NBAA), CRHOLM,  SIZE(CRHOLM))
                           ELSE
                              CALL FOCK_CHARGE_NOINT( WB(NBAA),WI(NBI), GWORK(1), &
                               CRHOLM, SIZE(CRHOLM))
                           ENDIF
                           ! Set phase e^iGr, where G = k_i - k_q - k_a
                           IF (LPHASE) THEN
                           CALL APPLY_PHASE( GRIDHF, CPHASE(1), GWORK(1), &
                             GWORK(1) )
                              !IF (KQ==1) WRITE(*,*)'error: no apply_phase need for kq=1'
                           ENDIF  

                        ! FFT{rho} to reciprocal space
                           CALL FFTEXT(WGWQ%NGVECTOR, WGWQ%NINDPW(1), &
                             GWORK(1),GCHGIA(1,NBAA,2),WGWQ%GRID,.FALSE.)
                           NFFT=NFFT+1
                        ! multiply with potential factor                        
                        ENDDO !NBAA (loop over all bands in NSTRIP)
                     
                        IF (ASSOCIATED(H)) THEN
                           CALL APPLY_PHASE_ONE_CENTER(WHF%WDES, H, CRHOIB(:,:), &
                             WHF%WDES%VKPT(:,KB)-WHF%WDES%VKPT(:,KI))
                        ENDIF
                     
                        DO NBAA=1,NSTRIPA
                           !copy response functions to RESPF_PW and RESPF_OC
			   CALL LOC2GLOB((NBA+NBAA-1),ME,WDES%NB_TOT,PROCS,1,rnba)
                          IF (((NBI+NBI_START-1)>VBMAX) .and. ((RNBA)> VBMAX)) CYCLE
                           tmp_FTOD_PW(1:NP,NBAA+NBA-1,NBI,2)=(GCHGIA(1:NP,NBAA,2))*(1.0_q/GRIDHF%NPLWV)
                           IF (ASSOCIATED(H)) THEN
                              tmp_FTOD_OC(1:H%TOTAL_ENTRIES,NBAA+NBA-1,NBI,2)=CRHOIB(1:H%TOTAL_ENTRIES,NBAA)
                           ENDIF

                        ENDDO
                     ENDDO !NBI (loop over nstrip bands only)
                  ENDDO !NBA loop over all bands
               
#endif
                  CALL REDISTRIBUTE_FTOD_GRID(W,WDES,KI,KQ,ISP,tmp_FTOD_PW,tmp_FTOD_OC,nbi_start,nbi_end)
               ENDDO !NBI_start loop over all bands
            ENDDO kiloop
         
            DEALLOCATE(GCHGIA,CRHOLM)
         
         ENDDO kqloop
         ENDDO spin
         
         IF (ALLOCATED(CRHOIA)) DEALLOCATE(CRHOIA)
         IF (ALLOCATED(CRHOIB)) DEALLOCATE(CRHOIB)
         IF (ALLOCATED(tmp_FTOD_OC)) DEALLOCATE(tmp_FTOD_OC)
         IF (ALLOCATED(tmp_FTOD_PW)) DEALLOCATE(tmp_FTOD_PW)
         
         DO NBI=1,NSTRIP
            CALL DELWAV(WI(NBI),.TRUE.)
         ENDDO
         DO NBA=1,NSTRIP
            CALL DELWAV(WA(NBA),.TRUE.)
         ENDDO
         DO NBA=1,NSTRIP
            CALL DELWAV(WB(NBA),.TRUE.)
         ENDDO
         DEALLOCATE(WI,WA,WB)
         
         
         RETURN
      END SUBROUTINE CALC_FTOD

      !***********************************************************************
      !This subroutine redistributes the fourier-transformed overlap integrals <i|G|a>
      !from the column-only process grid (CONTXT) to the quadratic 
      !process grid (CONTXT_GRID)
      !***********************************************************************

      SUBROUTINE REDISTRIBUTE_FTOD_GRID(W,WDES,KI,KQ,ISP,tmp_FTOD_PW,tmp_FTOD_OC,nbi_start,nbi_end)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE (wavespin) W
         INTEGER :: KA,KI,KQ,ISP, nbi_start, nbi_end
         GDEF, TARGET :: tmp_FTOD_OC(:,:,:,:)
         COMPLEX(q), TARGET :: tmp_FTOD_PW(:,:,:,:)
         INTEGER :: FTOD_PW_rows, FTOD_PW_cols,cc,FTOD_OC_rows, FTOD_OC_cols
         INTEGER :: NBI,FTOD_PW_rows_br,NBA
         
         KA=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KI)-WDES%VKPT(:,KQ),KPOINTS_FULL)
         
         ! Prepare array descriptors for ScaLAPACK
         call BLACS_GRIDINFO(CONTXT_GRID, NPROW, NPCOL, MYROW, MYCOL)
!         MB=MIN(((NGVECTOR)/NPROW),50)   !Row blocking size
!         FTOD_PW_rows = numroc(NGVECTOR,mb,MYROW,0,NPROW)
!         desc_FTOD_PW(3) = NGVECTOR       ! global number of rows
!         desc_FTOD_PW(5) = mb       ! row block size
!         desc_FTOD_PW(9) = MAX(1,FTOD_PW_rows) ! leading dimension of local array    

         desc_FTOD_PW_br(1) = 1              ! descriptor type
         desc_FTOD_PW_br(2) = contxt_cols         ! blacs context
         desc_FTOD_PW_br(3) = NGVECTOR    ! global number of rows
         desc_FTOD_PW_br(4) = (PROCS*WDES%NBANDS) ! global number of cols
         desc_FTOD_PW_br(5) = NGVECTOR     ! row block size
         desc_FTOD_PW_br(6) = 1             ! col block size
         desc_FTOD_PW_br(7) = 0              ! initial process row
         desc_FTOD_PW_br(8) = 0              ! initial process col
         desc_FTOD_PW_br(9) = MAX(1,(NGVECTOR)) ! leading dimension of local array
         !Distribute the fourier-transformed overlap integrals
         
         desc_FTOD_PW=desc_FTOD_PW_br
         desc_FTOD_PW(2)=contxt_cols
         desc_FTOD_PW(6)=(PROCS*WDES%NBANDS)
         desc_FTOD_PW(8) = PROCS_KPTS(KI)
         
         DO NBI=nbi_start,nbi_end
            DO cc=1,ncc
               CALL PZGEMR2D(NGVECTOR,(PROCS*WDES%NBANDS),tmp_FTOD_PW(1,1,(NBI-nbi_start+1),cc),1,1,&
                desc_FTOD_PW_br,FTOD_PW(1,1,cc),1,1,desc_FTOD_PW,contxt_grid)
               IF (ME==PROCS_KPTS(KI)) THEN
                  DO NBA=1,(PROCS*WDES%NBANDS)
                     IF (FILLED_MP2_ORBITAL(W%FERTOT(NBI,KI,1))) THEN
                        IF (FILLED_MP2_ORBITAL(W%FERTOT(NBA,KA,1))) FTOD_PW_IJ(:,NBA,NBI,KQ,MKPTS_KPTS(KI),cc)= &
                           (FTOD_PW(:,NBA,cc))
                        IF (EMPTY_MP2_ORBITAL(W%FERTOT(NBA,KA,1))) FTOD_PW_IA(:,NBA-VBMAX,NBI,KQ,MKPTS_KPTS(KI),cc)= &
                           (FTOD_PW(:,NBA,cc))
                     ELSE
                        IF (FILLED_MP2_ORBITAL(W%FERTOT(NBA,KA,1))) FTOD_PW_AI(:,NBA,NBI-VBMAX,KQ,MKPTS_KPTS(KI),cc)= &
                           (FTOD_PW(:,NBA,cc))
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         
         
!#ifdef gammareal         
!         desc_FTOD_PW(9) = MAX(1,m_ desc_FTOD_PW(9)) ! leading dimension of local array
!         desc_FTOD_PW(5) = m_ desc_FTOD_PW(5)
!         desc_FTOD_PW(3) = m_ desc_FTOD_PW(3)
!#else

!#endif
         !ONE-center part RESPF_OC
           
         IF (ASSOCIATED(H)) THEN
         
            ! Prepare array descriptors for ScaLAPACK
         
            desc_FTOD_OC_br(1) = 1            ! descriptor type
            desc_FTOD_OC_br(2) = contxt_cols       ! blacs context
            desc_FTOD_OC_br(3) = NHVECTOR     ! global number of rows
            desc_FTOD_OC_br(4) = (PROCS*WDES%NBANDS) ! global number of cols
            desc_FTOD_OC_br(5) = NHVECTOR     ! row block size
            desc_FTOD_OC_br(6) = 1            ! col block size
            desc_FTOD_OC_br(7) = 0            ! initial process row
            desc_FTOD_OC_br(8) = 0            ! initial process col
            desc_FTOD_OC_br(9) = MAX(1,(NHVECTOR)) ! leading dimension of local array

         desc_FTOD_OC=desc_FTOD_OC_br
         desc_FTOD_OC(2)=contxt_cols
         desc_FTOD_OC(6)=(PROCS*WDES%NBANDS)
         desc_FTOD_OC(8) = PROCS_KPTS(KI)

            DO NBI=nbi_start,nbi_end
               DO cc=1,2
!#ifdef gammareal
!                  CALL PDGEMR2D(NHVECTOR,(PROCS*WDES%NBANDS),tmp_FTOD_OC(1,1,NBI-nbi_start+1,cc),1,1,&
!                   desc_FTOD_OC_br,FTOD_OC(1,1,cc),1,1,desc_FTOD_OC,contxt_grid)
!#else
                  CALL PZGEMR2D(NHVECTOR,(PROCS*WDES%NBANDS),tmp_FTOD_OC(1,1,NBI-nbi_start+1,cc),1,1,&
                   desc_FTOD_OC_br,FTOD_OC(1,1,cc),1,1,desc_FTOD_OC,contxt_grid)
!#endif

               IF (ME==PROCS_KPTS(KI)) THEN
                  DO NBA=1,(PROCS*WDES%NBANDS)
                     IF (FILLED_MP2_ORBITAL(W%FERTOT(NBI,KI,1))) THEN
                        IF (FILLED_MP2_ORBITAL(W%FERTOT(NBA,KA,1))) FTOD_OC_IJ(:,NBA,NBI,KQ,MKPTS_KPTS(KI),cc)= &
                           FTOD_OC(:,NBA,cc)
                        IF (EMPTY_MP2_ORBITAL(W%FERTOT(NBA,KA,1))) FTOD_OC_IA(:,NBA-VBMAX,NBI,KQ,MKPTS_KPTS(KI),cc)= &
                           FTOD_OC(:,NBA,cc)
                     ELSE
                        IF (FILLED_MP2_ORBITAL(W%FERTOT(NBA,KA,1))) FTOD_OC_AI(:,NBA,NBI-VBMAX,KQ,MKPTS_KPTS(KI),cc)= &
                           FTOD_OC(:,NBA,cc)
                        IF (EMPTY_MP2_ORBITAL(W%FERTOT(NBA,KA,1))) FTOD_OC_AB(:,NBA-VBMAX,NBI-VBMAX,KQ,MKPTS_KPTS(KI),cc)= &
                           FTOD_OC(:,NBA,cc)
                     ENDIF
                  ENDDO
               ENDIF
               
               ENDDO
            ENDDO

         ENDIF

      END SUBROUTINE REDISTRIBUTE_FTOD_GRID

      !***********************************************************************
      !This routine initializes a process grid which is used for the 
      !pre-calculation of the fourier-transformed overlap integrals <i|G|a>
      !This process grid has NPROCS(number of processors used) columns and 1 row.
      !The assigned context variable is called CONTXT_COLS.
      !*********************************************************************** 

      SUBROUTINE INIT_BLACS_COLS()
         implicit none           
         INTEGER :: a_PRCS, i
         REAL :: NPCOL_TMP
         
         !first we create a column-only process grid in order to
         !calculate the <i|-G|a> and <j|G|b> quantities
         call BLACS_PINFO(ME,PROCS)
         nprow=1
         npcol=PROCS
         call BLACS_PINFO(ME,PROCS)
         call BLACS_GET     (0, 0, CONTXT_COLS)
         
         call BLACS_GRIDINIT(CONTXT_COLS, 'R', NPROW, NPCOL)
         call BLACS_GRIDINFO(CONTXT_COLS, NPROW, NPCOL, MYROW, MYCOL)
        
         !since, a quadratic process grid is going to be needed for
         !the matrix-matrix multiplications, its optimal row and column
         !numbers are estimated for the given number of processors
        a_PRCS=CEILING(SQRT(Real(PROCS)))
        IF (a_PRCS==SQRT(Real(PROCS))) THEN
           NPROW_GRID=a_PRCS
           NPCOL_TMP=a_PRCS
        ENDIF
        IF (a_PRCS/=SQRT(Real(PROCS))) THEN
           DO i=1,CEILING(SQRT(Real(PROCS)))
              NPCOL_TMP=PROCS/Real(i)
              IF ((NPCOL_TMP-INT(NPCOL_TMP))==0) THEN
                 NPROW_GRID=i
              ENDIF
           ENDDO
           NPCOL_TMP=PROCS/NPROW_GRID
        ENDIF
!   NPROW_GRID=PROCS
!   NPCOL_TMP=1
        IF (ABS(NPROW_GRID-a_PRCS)>(a_PRCS/2.0)) THEN
           WRITE(*,*)'The allocated number of CPUs does not allow for the use of an efficient process grid.'
           WRITE(*,*)'Suggested number of CPUs: 4, 16, 32, ....'
        ENDIF
      END SUBROUTINE INIT_BLACS_COLS

      !***********************************************************************
      !This routine initializes the process grid, which is used for the
      !matrix-matrix multiplications and their
      !block-cyclic data distribution. Note that the context variable
      !CONTXT_GRID is used for this grid.
      !The routine tries to create the most quadratic grid possible for the given
      !number of processors.
      !*********************************************************************** 

      SUBROUTINE INIT_BLACS_GRID(WDES)
         implicit none           
         TYPE(wavedes) WDES
         
         call BLACS_PINFO(ME,PROCS)
         
         nprow=NPROW_GRID
         npcol=PROCS/NPROW
         IF (MOD(PROCS/NPROW,1)/=0) THEN
            WRITE(*,*)'internal error in INIT_BLACS_GRID: Bad process grid'
         ENDIF
         
         IF (NGVECTOR>(PROCS*WDES%NBANDS)) THEN
            IF (NPROW_GRID<(PROCS/NPROW)) THEN
               NPCOL=NPROW_GRID
               NPROW=PROCS/NPROW_GRID
               NPROW_GRID=NPROW
            ENDIF
         ENDIF         
         IF (NGVECTOR<(PROCS*WDES%NBANDS)) THEN
            IF (NPROW_GRID>(PROCS/NPROW)) THEN
               NPCOL=NPROW_GRID
               NPROW=PROCS/NPROW_GRID
               NPROW_GRID=NPROW
            ENDIF
         ENDIF
         IF ((MYROW==0) .AND. (MYCOL==0)) THEN
           !WRITE(*,*)'You are using',PROCS,' processors.'
           WRITE(*,'(A,I3,A,I3,A)')'The allocated processors form a',NPROW_GRID,'x',NPCOL,' grid.'
         ENDIF

         call BLACS_PINFO(ME,PROCS)
         call BLACS_GET     (0, 0, CONTXT_GRID)
         call BLACS_GRIDINIT(CONTXT_GRID, 'R', NPROW, NPCOL)
         call BLACS_GRIDINFO(CONTXT_GRID, NPROW, NPCOL, MYROW, MYCOL)

      END SUBROUTINE INIT_BLACS_GRID


      SUBROUTINE INIT_MY_KPOINTS(WDES,WGW)
         implicit none
         TYPE(wavedes) WDES
         TYPE(wavedes) WGW
         INTEGER :: MNK, ATPROCS

         call BLACS_PINFO(ME,PROCS)

         MY_NKPTS=0
         ATPROCS=-1
         DO MNK=1,WDES%NKPTS
            ATPROCS=MOD((ATPROCS+1),PROCS)
            IF (ATPROCS==ME) MY_NKPTS=MY_NKPTS+1        
         ENDDO
         ALLOCATE(KPTS_MKPTS(MY_NKPTS))
         ALLOCATE(PROCS_KPTS(WDES%NKPTS))
         ALLOCATE(MKPTS_KPTS(WDES%NKPTS))

         MKPTS_KPTS(:)=0
         MY_NKPTS=0
         ATPROCS=-1
         DO MNK=1,WDES%NKPTS
            ATPROCS=MOD((ATPROCS+1),PROCS)
            IF (ATPROCS==ME) MY_NKPTS=MY_NKPTS+1
            IF (ATPROCS==ME) KPTS_MKPTS(MY_NKPTS)=MNK
            IF (ATPROCS==ME) MKPTS_KPTS(MNK)=MY_NKPTS
            PROCS_KPTS(MNK)=ATPROCS
         ENDDO

      END SUBROUTINE INIT_MY_KPOINTS
      
      !***********************************************************************
      !Allocate the matrices T2AMPLITUDES(i,j,a,b)
      !in the block-cyclic data distribution
      !and setup its descriptor(desc_T2AMPLITUDES),
      !which is required by the scaLAPACK routines
      !***********************************************************************
      
      SUBROUTINE SETUP_T2AMPLITUDES(WDES,NBANDSGW)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER MB,NB,T2AMPLITUDES_ROWS,T2AMPLITUDES_COLS,NBANDSGW
         CALL BLACS_GRIDINFO(CONTXT_GRID, NPROW, NPCOL, MYROW, MYCOL)
         NPROW=NPROW_GRID
         NPCOL=PROCS/NPROW
         
         NCLENGTH=(NBANDSGW-VBMAX)*(PROCS*WDES%NBANDS-VBMAX)+(NBANDSGW-VBMAX)*(PROCS*WDES%NBANDS-NBANDSGW)
         NCLENGTH_N=(NBANDSGW-VBMAX)*(PROCS*WDES%NBANDS-VBMAX)

         !IF (ALLOCATED(T2)) THEN
         !   DEALLOCATE(T2)
         !ENDIF
         !IF (ALLOCATED(T2_N)) THEN
         !   DEALLOCATE(T2_N)
         !ENDIF

         ALLOCATE(T2((PROCS*WDES%NBANDS-VBMAX),(PROCS*WDES%NBANDS-VBMAX),&
           VBMAX,VBMAX,WDES%NKPTS,MY_NKPTS))
         ALLOCATE(CHI_GAI(NGVECTOR,(PROCS*WDES%NBANDS-VBMAX),VBMAX,WDES%NKPTS,MY_NKPTS,ncc))
         ALLOCATE(GAI_TMP(NGVECTOR,(PROCS*WDES%NBANDS-VBMAX),VBMAX,WDES%NKPTS))
         ALLOCATE(AIG((PROCS*WDES%NBANDS-VBMAX),VBMAX,NGVECTOR))
!         ALLOCATE(T2_N((NBANDSGW-VBMAX)*(PROCS*WDES%NBANDS-VBMAX),&
!           VBMAX,VBMAX,WDES%NKPTS,MY_NKPTS))
!         ALLOCATE(PT2_N(NCLENGTH,&
!           VBMAX,VBMAX,WDES%NKPTS,MY_NKPTS))
!         ALLOCATE(T2_N_TMP((NBANDSGW-VBMAX)*(PROCS*WDES%NBANDS-VBMAX),&
!           VBMAX,VBMAX,WDES%NKPTS,WDES%NKPTS))
         ALLOCATE(T2_TMP((PROCS*WDES%NBANDS-VBMAX),(PROCS*WDES%NBANDS-VBMAX),&
            VBMAX,VBMAX))

         !ALLOCATE(CHI_CKAI_KTMP2((PROCS*WDES%NBANDS-VBMAX),VBMAX,(PROCS*WDES%NBANDS-VBMAX),VBMAX,WDES%NKPTS))
!         ALLOCATE(K_KI(VBMAX,VBMAX,WDES%NKPTS))
!         ALLOCATE(K_MP2_KI(VBMAX,VBMAX,WDES%NKPTS))
!         ALLOCATE(K_AC((PROCS*WDES%NBANDS-VBMAX),(PROCS*WDES%NBANDS-VBMAX),WDES%NKPTS))
!         ALLOCATE(K_KC((VBMAX),(PROCS*WDES%NBANDS-VBMAX),WDES%NKPTS))
!         ALLOCATE(L_KI(VBMAX,VBMAX,WDES%NKPTS))
!         ALLOCATE(L_AC((PROCS*WDES%NBANDS-VBMAX),(PROCS*WDES%NBANDS-VBMAX),WDES%NKPTS))
         ALLOCATE(VV_S((PROCS*WDES%NBANDS-VBMAX),(PROCS*WDES%NBANDS-VBMAX)))
         ALLOCATE(VV2_S((PROCS*WDES%NBANDS-VBMAX),(PROCS*WDES%NBANDS-VBMAX)))
         ALLOCATE(VV((PROCS*WDES%NBANDS-VBMAX),(PROCS*WDES%NBANDS-VBMAX)))
         ALLOCATE(VV2((PROCS*WDES%NBANDS-VBMAX),(PROCS*WDES%NBANDS-VBMAX)))
         ALLOCATE(TE4O(VBMAX,VBMAX))
!         ALLOCATE(T1((PROCS*WDES%NBANDS-VBMAX),VBMAX,WDES%NKPTS))
!         ALLOCATE(T1_N((PROCS*WDES%NBANDS-VBMAX),VBMAX,WDES%NKPTS))
!         ALLOCATE(T1_T((PROCS*WDES%NBANDS-VBMAX),VBMAX,WDES%NKPTS))
!         ALLOCATE(VVOO_BA((PROCS*WDES%NBANDS-VBMAX),(PROCS*WDES%NBANDS-VBMAX),VBMAX,VBMAX,WDES%NKPTS,MY_NKPTS))
         
         ALLOCATE(VVOO((PROCS*WDES%NBANDS-VBMAX),(PROCS*WDES%NBANDS-VBMAX),VBMAX,VBMAX))
         ALLOCATE(OOVV(VBMAX,VBMAX,(PROCS*WDES%NBANDS-VBMAX),(PROCS*WDES%NBANDS-VBMAX)))
         ALLOCATE(OOVV_S(VBMAX,VBMAX,(PROCS*WDES%NBANDS-VBMAX),(PROCS*WDES%NBANDS-VBMAX)))
         ALLOCATE(OOVV2_S(VBMAX,VBMAX,(PROCS*WDES%NBANDS-VBMAX),(PROCS*WDES%NBANDS-VBMAX)))
         ALLOCATE(VVOO_S((PROCS*WDES%NBANDS-VBMAX),(PROCS*WDES%NBANDS-VBMAX),VBMAX,VBMAX))
         ALLOCATE(VOVO_S((PROCS*WDES%NBANDS-VBMAX),VBMAX,(PROCS*WDES%NBANDS-VBMAX),VBMAX))
         ALLOCATE(VOVO((PROCS*WDES%NBANDS-VBMAX),VBMAX,(PROCS*WDES%NBANDS-VBMAX),VBMAX))
         ALLOCATE(VOVO2_S((PROCS*WDES%NBANDS-VBMAX),VBMAX,(PROCS*WDES%NBANDS-VBMAX),VBMAX))
         ALLOCATE(VOVO3_S((PROCS*WDES%NBANDS-VBMAX),VBMAX,(PROCS*WDES%NBANDS-VBMAX),VBMAX))
         ALLOCATE(OVOV(VBMAX,(PROCS*WDES%NBANDS-VBMAX),VBMAX,(PROCS*WDES%NBANDS-VBMAX)))
         ALLOCATE(OOOV(VBMAX,VBMAX,VBMAX,(PROCS*WDES%NBANDS-VBMAX)))
         ALLOCATE(OOOV_S(VBMAX,VBMAX,VBMAX,(PROCS*WDES%NBANDS-VBMAX)))   
         ALLOCATE(VOOV_S((PROCS*WDES%NBANDS-VBMAX),VBMAX,VBMAX,(PROCS*WDES%NBANDS-VBMAX)))
         ALLOCATE(KLIJ(VBMAX,VBMAX,VBMAX,VBMAX))
         ALLOCATE(JLIK(VBMAX,VBMAX,VBMAX,VBMAX))
!         T2AMPLITUDES=zero
!         T2AMPLITUDES_new=zero

      END SUBROUTINE SETUP_T2AMPLITUDES

      !***********************************************************************
      !Allocate the arrays RESPF_PW and 
      !RESPF_OC in the block-cyclic data distribution
      !and setup their descriptors(desc_RESPF_PW and desc_RESPF_OC),
      !which are required by the scalaLAPACK routines
      !***********************************************************************
      
      SUBROUTINE SETUP_FTOD(WDES,NBANDSGW)
         IMPLICIT NONE
         TYPE(wavedes) WDES      
         INTEGER NBANDSGW

         call BLACS_GRIDINFO(CONTXT_GRID, NPROW, NPCOL, MYROW, MYCOL)
         !Blocking size for block-cyclic distribution of RESPF_PW          
         
         desc_FTOD_PW(1) = 1              ! descriptor type
         desc_FTOD_PW(2) = CONTXT_GRID      ! blacs context
         desc_FTOD_PW(3) = NGVECTOR       ! global number of rows
         desc_FTOD_PW(4) = (PROCS*WDES%NBANDS)   ! global number of cols
         desc_FTOD_PW(5) = NGVECTOR       ! row block size
         desc_FTOD_PW(6) = (PROCS*WDES%NBANDS)       ! column block size
         desc_FTOD_PW(7) = 0              ! initial process row
         desc_FTOD_PW(8) = 0              ! initial process column
         desc_FTOD_PW(9) = MAX(1,NGVECTOR) ! leading dimension of local array

         
         ALLOCATE(FTOD_PW(NGVECTOR,(PROCS*WDES%NBANDS),ncc))
         ALLOCATE(FTOD_PW_IA(NGVECTOR,(PROCS*WDES%NBANDS-VBMAX),VBMAX,WDES%NKPTS,MY_NKPTS,ncc))
         ALLOCATE(FTOD_PW_AI(NGVECTOR,VBMAX,(PROCS*WDES%NBANDS-VBMAX),WDES%NKPTS,MY_NKPTS,ncc))
         ALLOCATE(FTOD_PW_IJ(NGVECTOR,VBMAX,VBMAX,WDES%NKPTS,MY_NKPTS,ncc))
         !ALLOCATE(FTOD_PW_IJ_TMP(NGVECTOR,VBMAX,VBMAX,WDES%NKPTS))
         
         ALLOCATE(PW_AI_TMP(NGVECTOR,VBMAX,(PROCS*WDES%NBANDS-VBMAX),WDES%NKPTS))
         ALLOCATE(PW_IA_TMP(NGVECTOR,(PROCS*WDES%NBANDS-VBMAX),VBMAX,WDES%NKPTS))
         ALLOCATE(PW_IA_TMP2(NGVECTOR,(PROCS*WDES%NBANDS-VBMAX),VBMAX,WDES%NKPTS))
         ALLOCATE(PW_IJ_TMP(NGVECTOR,VBMAX,VBMAX,WDES%NKPTS))
	 
         IF (ASSOCIATED(H)) THEN         
            !Blocking size for block-cyclic distribution of RESPF_PW    
            
            desc_FTOD_OC(1) = 1              ! descriptor type
            desc_FTOD_OC(2) = CONTXT_GRID        ! blacs context
            desc_FTOD_OC(3) = NHVECTOR       ! global number of rows
            desc_FTOD_OC(4) = (PROCS*WDES%NBANDS)   ! global number of cols
            desc_FTOD_OC(5) = NHVECTOR             ! row block size
            desc_FTOD_OC(6) = (PROCS*WDES%NBANDS)             ! column block size
            desc_FTOD_OC(7) = 0              ! initial process row
            desc_FTOD_OC(8) = 0              ! initial process column
            desc_FTOD_OC(9) = MAX(1,NHVECTOR) ! leading dimension of local array
            
            
            ALLOCATE(FTOD_OC(NHVECTOR,(PROCS*WDES%NBANDS),2))
            ALLOCATE(FTOD_OC_IA(NHVECTOR,(PROCS*WDES%NBANDS-VBMAX),VBMAX,WDES%NKPTS,MY_NKPTS,2))
            ALLOCATE(FTOD_OC_AI(NHVECTOR,VBMAX,(PROCS*WDES%NBANDS-VBMAX),WDES%NKPTS,MY_NKPTS,2))
            ALLOCATE(FTOD_OC_IJ(NHVECTOR,VBMAX,VBMAX,WDES%NKPTS,MY_NKPTS,2))
            ALLOCATE(OC_IJ_TMP(NHVECTOR,VBMAX,VBMAX,WDES%NKPTS))
            ALLOCATE(OC_AI_TMP(NHVECTOR,VBMAX,(PROCS*WDES%NBANDS-VBMAX),WDES%NKPTS))
            ALLOCATE(OC_IA_TMP(NHVECTOR,(PROCS*WDES%NBANDS-VBMAX),VBMAX,WDES%NKPTS))
            ALLOCATE(OC_IA_TMP2(NHVECTOR,(PROCS*WDES%NBANDS-VBMAX),VBMAX,WDES%NKPTS))

         ENDIF

      END SUBROUTINE SETUP_FTOD
      
      !***********************************************************************
      ! TRANSFORM LOCAL TO GLOBAL ARRAY INDEX FOR BLOCK-CYCLIC ARRAY
      ! DISTRIBUTION
      !***********************************************************************
      SUBROUTINE LOC2GLOB(li,p,n,np,nb,gi)
         IMPLICIT NONE
         INTEGER :: li   ! local index
         INTEGER :: p    ! index in processor grid (either MYROW or MYCOL)
         INTEGER :: n    ! global array dimension
         INTEGER :: np   ! processor array dimension
         INTEGER :: nb   ! blocking size
         INTEGER :: gi   ! global index
         INTEGER :: litmp   
 
         litmp = li-1
         gi=(((litmp/nb)*np)+p)*nb+mod(litmp,nb)+1
         RETURN
      END SUBROUTINE LOC2GLOB

      subroutine APPEND_CDER_TO_FTOD(WDES,W,LATT_CUR,dir)
         USE constant
         USE full_kpoints
         USE mkpoints
         USE wave
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE (wavespin) W
         TYPE(latt) LATT_CUR
         integer :: dir !xyz direction of the dipolemoments
         integer :: mni,mna,mki,mkq,mng,mrng,mrna,mrni
         INTEGER :: FTOD_ROWS,TWOE4ORBITAL_COLS
         GDEF :: CDER_BETWEEN_STATE_IA(3),CDER_BETWEEN_STATE_AI(3)
         real(qs) :: ediff_th


         do mki=1,MY_NKPTS

         MKQ=KPOINT_IN_FULL_GRID(W%WDES%VKPT(:,KPTS_MKPTS(mki))-W%WDES%VKPT(:,KPTS_MKPTS(mki)),KPOINTS_FULL) !just to find the gamma-point
            do mni=1,VBMAX
               do mna=1,(PROCS*WDES%NBANDS-VBMAX)

                  CDER_BETWEEN_STATE_IA(:)=zero
                  CDER_BETWEEN_STATE_AI(:)=zero

                  CALL  CDER_BETWEEN_STATES_ROTATED( &
                   CDER_BETWEEN_STATE_IA,LATT_CUR,KPTS_MKPTS(mki), 1, mna+VBMAX, mni)
                        !CALL  CDER_BETWEEN_STATES_ROTATED( &
                      !CDER_BETWEEN_STATE_AI,LATT_CUR,mki, 1, mni, mrna)
                      
                  FTOD_PW_IA(ngvector,mna,mni,mkq,mki,1)=(CDER_BETWEEN_STATE_IA(dir)*EDEPS/LATT_CUR%OMEGA)
                  FTOD_PW_IA(ngvector,mna,mni,mkq,mki,2)=conjg(CDER_BETWEEN_STATE_IA(dir))
                  
                  FTOD_PW_AI(ngvector,mni,mna,mkq,mki,1)=conjg((CDER_BETWEEN_STATE_IA(dir)*EDEPS/LATT_CUR%OMEGA))
                  FTOD_PW_AI(ngvector,mni,mna,mkq,mki,2)=(CDER_BETWEEN_STATE_IA(dir))

               enddo
            enddo
         enddo

      end subroutine APPEND_CDER_TO_FTOD  

      Subroutine CALC_EFOCK(WDES,WGW,IO)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(wavedes) WGW
         TYPE (in_struct) IO
         INTEGER :: KQ,KI,KJ,KQ_,NI,NJ
!    write(*,*)'me,my_kpts',me,my_kpts,PROCS_KP
         TE4O=zero

         call BLACS_GRIDINFO(CONTXT_COLS, NPROW, NPCOL, MYROW, MYCOL)
         DO KI=1,WDES%NKPTS
!       IF (PROCS_KPTS(KI)==ME) THEN
!          CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, VBMAX*VBMAX*WDES%NKPTS, FTOD_PW_IJ(1,1,1,1,MKPTS_KPTS(KI),1), NGVECTOR)
!          FTOD_PW_IJ_TMP(:,:,:,:)=FTOD_PW_IJ(:,:,:,:,MKPTS_KPTS(KI),1)
!       ELSE
!          CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, VBMAX*VBMAX*WDES%NKPTS, FTOD_PW_IJ_TMP(1,1,1,1), NGVECTOR,0,PROCS_KPTS(KI))
!       ENDIF
         CALL BCAST2ALL_FTOD_IJ(WDES,KI,1,PW_IJ_TMP,OC_IJ_TMP)
         
         DO KJ=1,MY_NKPTS
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
            KQ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KI)-WDES%VKPT(:,KPTS_MKPTS(KJ)),KPOINTS_FULL)

!       CALL ZGEMM(trans,'n',(1),(1),&
!                 (NGVECTOR),(1._q,0._q),FTOD_PW_IJ(1,NJ,NI,KQ,KI,1),(NGVECTOR),&
!                 FTOD_PW_IJ(1,NI,NJ,KQ,MY_KPTS(KJ),2),(NGVECTOR),&
!                 (1._q,0._q), TE4O(1,1),(1))

!            CALL ZGEMM(trans,'n',(1),(1),&
!                 (NGVECTOR),(1._q,0._q),PW_IJ_TMP(1,NJ,NI,KQ),(NGVECTOR),&
!                 FTOD_PW_IJ(1,NI,NJ,KQ,KJ,2),(NGVECTOR),&
!                 (1._q,0._q), TE4O(1,1),(1))
                 
            IF (ASSOCIATED(H)) THEN
!            CALL ZGEMM(trans,'n',(1),(1),&
!                 (NHVECTOR),(1._q,0._q),OC_IJ_TMP(1,NJ,NI,KQ),(NHVECTOR),&
!                 FTOD_OC_IJ(1,NI,NJ,KQ,KJ,2),(NHVECTOR),&
!                 (1._q,0._q), TE4O(1,1),(1))
            ENDIF

               CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
                PW_IJ_TMP(1,NJ,NI,KQ),OC_IJ_TMP(1,NJ,NI,KQ), &
                (1)*1,&
                FTOD_PW_IJ(1,NI,NJ,KQ,KJ,2),FTOD_OC_IJ(1,NI,NJ,KQ,KJ,2),&
                (1)*1,TE4O(1,1),(1.0_q,0.0_q))

         ENDDO
         ENDDO
         ENDDO
         ENDDO
         
         TE4O(1,1)=GCONJG(TE4O(1,1)*KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(1))
         CALLMPI( M_sum_z(WGW%COMM_INTER, TE4O(1,1), 1))
         !IF (IO%IU0>=0) write(*,*)'EFOCK',TE4O(1,1)
         !IF (IO%IU0>=0) write(*,*)'VBMAX=',VBMAX

      END SUBROUTINE CALC_EFOCK

      SUBROUTINE CALC_CCSD_ENERGY(WDES,WGW,W,IO,NBANDSGW,MKB)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(wavedes) WGW
         TYPE(wavespin) W 
         TYPE (in_struct) IO
         INTEGER, OPTIONAL :: MKB
         INTEGER :: KQ,KI,KJ,KQ_,NI,NJ,KA,KB,NA,NB,NC,MKI
         INTEGER :: NBANDSGW
         GDEFS :: ETMP,KW3
         
         KW3=KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(1)
         
         E_CCSD=zero
         E_CCSD_X=zero
         ETMP=zero

         DO KB=1,WDES%NKPTS

!           IF (.not. CANONICAL) THEN
!              IF (.not. MP2) THEN
!              CALL CGEMM('C','n',1,1, &
!                (VBMAX)*(PROCS*WDES%NBANDS-VBMAX),(2._qs,0._qs),F_AI(1,1,KB),&
!                (VBMAX)*(PROCS*WDES%NBANDS-VBMAX), &
!                T1(1,1,KB),(VBMAX)*(PROCS*WDES%NBANDS-VBMAX), &
!                (1._qs,0._qs), ETMP,1)
!              ENDIF
!           ENDIF

            CALL BCAST2ALL_FTOD_AI(WDES,KB,2,PW_AI_TMP,OC_AI_TMP)
            DO KA=1,MY_NKPTS
!            DO KQ=1,WDES%NKPTS
               KQ=KPOINT_IN_FULL_GRID(-WDES%VKPT(:,GKQ),KPOINTS_FULL)
               OOVV=zero
               KI=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KPTS_MKPTS(KA))- &
                  WDES%VKPT(:,KQ),KPOINTS_FULL)
               KJ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KB)+ &
                  WDES%VKPT(:,KQ),KPOINTS_FULL)
               KQ_=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KPTS_MKPTS(KA))- &
                  WDES%VKPT(:,KJ),KPOINTS_FULL)
               DO NA=1,(PROCS*WDES%NBANDS-VBMAX)
               DO NB=1,(PROCS*WDES%NBANDS-VBMAX)
                  CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
                   PW_AI_TMP(1,1,NB,KQ_),OC_AI_TMP(1,1,NB,KQ_), &
                   VBMAX,&
                   FTOD_PW_AI(1,1,NA,KQ_,KA,1),FTOD_OC_AI(1,1,NA,KQ_,KA,1),&
                   VBMAX,OOVV(1,1,NA,NB),zero)
               ENDDO
               ENDDO

               DO NI=1,VBMAX
               DO NJ=1,VBMAX
                  NC=0
                  DO NB=1,(PROCS*WDES%NBANDS-VBMAX)
                  DO NA=1,(PROCS*WDES%NBANDS-VBMAX) !(NBANDSGW-VBMAX)
 !                    NC=NC+1
                     E_CCSD_X=E_CCSD_X+(OOVV(NI,NJ,NA,NB))* &
                         (T2(NA,NB,NI,NJ,KB,KA))*KPOINTS_FULL%WTKPT(1) &
                        *KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(1)
                  ENDDO
                  ENDDO
!                  IF (NBANDSGW<PROCS*WDES%NBANDS) THEN
!                     DO NB=1,(NBANDSGW-VBMAX)
!                     DO NA=(NBANDSGW-VBMAX)+1,(PROCS*WDES%NBANDS-VBMAX)
!                        NC=NC+1
!                        E_CCSD_X=E_CCSD_X+(OOVV(NI,NJ,NA,NB))* &
!                           (T2(NC,NI,NJ,KI,KJ,KA))*KPOINTS_FULL%WTKPT(1) &
!                           *KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(1)
!                     ENDDO
!                     ENDDO
!                  ENDIF
               ENDDO
               ENDDO


               OOVV=-CONJG((0.5_qs,0.0_qs)*OOVV)

               DO NA=1,(PROCS*WDES%NBANDS-VBMAX)
               DO NB=1,(PROCS*WDES%NBANDS-VBMAX)
                  CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
                   FTOD_PW_AI(1,1,NA,KQ,KA,1),FTOD_OC_AI(1,1,NA,KQ,KA,1),&
                   VBMAX,&
                   PW_AI_TMP(1,1,NB,KQ),OC_AI_TMP(1,1,NB,KQ), &
                   VBMAX,OOVV(1,1,NA,NB),one)
               ENDDO
               ENDDO

               OOVV=CONJG((2.0_qs,0.0_qs)*OOVV)

!               CALL SORT_O1V1O2V2_V2V1O2O1(WDES,OVOV, &
!                 VVOO)

!              VVOO=VVOO

               DO NI=1,VBMAX
               DO NJ=1,VBMAX
                  NC=0
                  DO NB=1,(PROCS*WDES%NBANDS-VBMAX)
                  DO NA=1,(PROCS*WDES%NBANDS-VBMAX) !(NBANDSGW-VBMAX)
!                     NC=NC+1
                     E_CCSD=E_CCSD+OOVV(NI,NJ,NA,NB)* &
                         (T2(NA,NB,NI,NJ,KB,KA))*KPOINTS_FULL%WTKPT(1) &
                        *KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(1)
                  ENDDO
                  ENDDO
!                  IF (NBANDSGW<PROCS*WDES%NBANDS) THEN
!                  DO NB=1,(NBANDSGW-VBMAX)
!                  DO NA=(NBANDSGW-VBMAX)+1,(PROCS*WDES%NBANDS-VBMAX)
!                     NC=NC+1
!                     E_CCSD=E_CCSD+OOVV(NI,NJ,NA,NB)* &
!                         (T2(NC,NI,NJ,KI,KJ,KA))*KPOINTS_FULL%WTKPT(1) &
!                        *KPOINTS_FULL%WTKPT(1)*KPOINTS_FULL%WTKPT(1)
!                  ENDDO
!                  ENDDO
!                  ENDIF
               ENDDO
               ENDDO
               
!            ENDDO
            ENDDO
         ENDDO
         CALL M_sum_z(WGW%COMM_INTER, E_CCSD, 1)
         CALL M_sum_z(WGW%COMM_INTER, E_CCSD_X, 1)
         E_CCSD=E_CCSD+ETMP*KPOINTS_FULL%WTKPT(1)
         E_CCSD_X=E_CCSD_X+ETMP*KPOINTS_FULL%WTKPT(1)
         IF (IO%IU0>=0) WRITE(IO%IU0,*)'RPA+SOSEX corr. energy ',E_CCSD
         IF (IO%IU0>=0) WRITE(IO%IU6,*)'RPA+SOSEX corr. energy ',E_CCSD
         IF (IO%IU0>=0) WRITE(IO%IU0,*)'SOSEX energy ',-E_CCSD_X
         IF (IO%IU0>=0) WRITE(IO%IU6,*)'SOSEX energy ',-E_CCSD_X
         IF (ITERATION==1) E_MP2=E_CCSD
         

      END SUBROUTINE CALC_CCSD_ENERGY

      Subroutine BCAST2ALL_FTOD_IJ(WDES,KI,MNCC,PW_IJ_TMP,OC_IJ_TMP)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: KI,MNCC
         GDEF :: PW_IJ_TMP(:,:,:,:), OC_IJ_TMP(:,:,:,:)

            call BLACS_PINFO(ME,PROCS)

            IF (PROCS_KPTS(KI)==ME) THEN
               !PW_IJ_TMP(:,:,:,:)=one
               !write(*,*)'test1'
               !PW_IJ_TMP(:,:,:,:)=FTOD_PW_IJ(:,:,:,:,MKPTS_KPTS(KI),1)
               !write(*,*)'test2',FTOD_PW_IJ(1,1,1,1,MKPTS_KPTS(KI),MIN(NCC,MNCC)),CONTXT_COLS,NGVECTOR,VBMAX*VBMAX*WDES%NKPTS
               CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, VBMAX*VBMAX*WDES%NKPTS, FTOD_PW_IJ(1,1,1,1,MKPTS_KPTS(KI),MIN(NCC,MNCC)), NGVECTOR)
               PW_IJ_TMP(:,:,:,:)=FTOD_PW_IJ(:,:,:,:,MKPTS_KPTS(KI),MIN(NCC,MNCC)) !FTOD_PW_IJ(1,1,1,1,MKPTS_KPTS(KI),MIN(NCC,MNCC))
            ELSE
               CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, VBMAX*VBMAX*WDES%NKPTS, PW_IJ_TMP(1,1,1,1), NGVECTOR,0,PROCS_KPTS(KI))
            ENDIF
            
            IF (ASSOCIATED(H)) THEN
            
               IF (PROCS_KPTS(KI)==ME) THEN
!#ifdef gammareal
!                  CALL DGEBS2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, VBMAX*VBMAX*WDES%NKPTS, FTOD_OC_IJ(1,1,1,1,MKPTS_KPTS(KI),MNCC), NHVECTOR)
!                  OC_IJ_TMP(:,:,:,:)=FTOD_OC_IJ(:,:,:,:,MKPTS_KPTS(KI),MNCC)
!#else
                  CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, VBMAX*VBMAX*WDES%NKPTS, FTOD_OC_IJ(1,1,1,1,MKPTS_KPTS(KI),MNCC), NHVECTOR)
                  OC_IJ_TMP(:,:,:,:)=FTOD_OC_IJ(:,:,:,:,MKPTS_KPTS(KI),MNCC)
!#endif
               ELSE
!#ifdef gammareal
!                  CALL DGEBR2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, VBMAX*VBMAX*WDES%NKPTS, OC_IJ_TMP(1,1,1,1), NHVECTOR,0,PROCS_KPTS(KI))
!#else
                  CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, VBMAX*VBMAX*WDES%NKPTS, OC_IJ_TMP(1,1,1,1), NHVECTOR,0,PROCS_KPTS(KI))
!#endif
               ENDIF
            
            ENDIF
      
      END Subroutine BCAST2ALL_FTOD_IJ


      SUBROUTINE ORDER_FTOD_IA(WDES,KI,MNCC,PW_IA_TMP,OC_IA_TMP)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: KI,MNCC, ATPRO
         GDEF :: PW_IA_TMP(:,:,:,:), OC_IA_TMP(:,:,:,:)
         INTEGER :: ORDERLIST(PROCS),tag,ierror,myrank
         INTEGER :: req(PROCS),recreq
         INTEGER :: status(MPI_STATUS_SIZE, PROCS)
         
         tag=1
         call BLACS_PINFO(ME,PROCS)
         CALL MPI_COMM_RANK(WDES%comm, myrank, ierror)
         me=myrank
         ORDERLIST=-1
         IF (ME==0) THEN
            ORDERLIST(1)=KI
            DO ATPRO=1,PROCS-1
               !write(*,*)'Im:',ME,'and wait for ATPRO to send',ATPRO
               CALL MPI_RECV( ORDERLIST(ATPRO+1), 1, MPI_INTEGER, ATPRO, tag, &
                  WDES%COMM, status(1,ATPRO), ierror) !req(ATPRO), ierror )
!               write(*,*)'I and wait for ATPRO to send',ATPRO
               !CALL MPI_WAIT(req(ATPRO), status(1,ATPRO), ierror)
            ENDDO
         ELSE
            !write(*,*)'Im:',ME,'and send to 0'
            call MPI_ISEND(KI, 1, MPI_INTEGER, 0, tag, &
                 WDES%COMM,req(ME+1), ierror)
            CALL MPI_WAIT(req(ME+1), status(1,ME+1), ierror)
         ENDIF

         !write(*,*)'test3,me',me,ORDERLIST
         
         IF (ME==0) THEN
            CALL IGEBS2D(CONTXT_COLS, 'All', 'i-ring', PROCS, 1, ORDERLIST(1), PROCS)
         ELSE
            CALL IGEBR2D(CONTXT_COLS, 'All', 'i-ring', PROCS, 1, ORDERLIST(1), PROCS, 0,0)
         ENDIF
         !write(*,*)'test',ORDERLIST

         DO ATPRO=1,PROCS
            !write(*,*)'atpro:',ATPRO
            IF ((PROCS_KPTS(ORDERLIST(ATPRO))==ME) .and. ((ATPRO-1)/=ME)) THEN
               !write(*,*)'test2',ORDERLIST
               !write(*,*)'Im:',ME,'and send to',ATPRO
               call MPI_ISEND(FTOD_PW_IA(1,1,1,1,MKPTS_KPTS(ORDERLIST(ATPRO)),MIN(NCC,MNCC)), &
                           NGVECTOR*(PROCS*WDES%NBANDS-VBMAX)*VBMAX*WDES%NKPTS, MPI_COMPLEX16, ATPRO-1, tag, &
                           WDES%COMM, req(ME+1), ierror)
               CALL MPI_WAIT(req(ME+1), status(1,ME+1), ierror)
            ENDIF
            IF ((PROCS_KPTS(KI)/=ME) .and. ((ATPRO-1)==ME)) THEN
               call MPI_RECV(PW_IA_TMP(1,1,1,1), &
                           NGVECTOR*(PROCS*WDES%NBANDS-VBMAX)*VBMAX*WDES%NKPTS, MPI_COMPLEX16, PROCS_KPTS(KI), tag, &
                           WDES%COMM, status(1,ME+1), ierror) !, recreq, ierror)
            ENDIF
         ENDDO
         IF ((PROCS_KPTS(KI)==ME)) THEN
            PW_IA_TMP(:,:,:,:)=FTOD_PW_IA(:,:,:,:,MKPTS_KPTS(KI),MIN(NCC,MNCC))
         ENDIF
         !write(*,*)'finish',orderlist
         !IF (PROCS_KPTS(KI)/=ME) CALL MPI_WAIT(recreq, status(1,ME), ierror)

!           PW_IA_TMP(:,:,:,:)=zero
!           OC_IA_TMP(:,:,:,:)=zero

!            IF (PROCS_KPTS(KI)==ME) THEN
!               CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, (PROCS*WDES%NBANDS-VBMAX)*VBMAX*WDES%NKPTS, FTOD_PW_IA(1,1,1,1,MKPTS_KPTS(KI),MIN(NCC,MNCC)), NGVECTOR)
!               PW_IA_TMP(:,:,:,:)=FTOD_PW_IA(:,:,:,:,MKPTS_KPTS(KI),MIN(NCC,MNCC))
!            ELSE
!               CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, (PROCS*WDES%NBANDS-VBMAX)*VBMAX*WDES%NKPTS, PW_IA_TMP(1,1,1,1), NGVECTOR,0,PROCS_KPTS(KI))
!            ENDIF
            
!            IF (ASSOCIATED(H)) THEN
!               IF (PROCS_KPTS(KI)==ME) THEN
!                  CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, (PROCS*WDES%NBANDS-VBMAX)*VBMAX*WDES%NKPTS, FTOD_OC_IA(1,1,1,1,MKPTS_KPTS(KI),MNCC), NHVECTOR)
!                  OC_IA_TMP(:,:,:,:)=FTOD_OC_IA(:,:,:,:,MKPTS_KPTS(KI),MNCC)
!               ELSE
!                  CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, (PROCS*WDES%NBANDS-VBMAX)*VBMAX*WDES%NKPTS, OC_IA_TMP(1,1,1,1), NHVECTOR,0,PROCS_KPTS(KI))
!               ENDIF
            
!            ENDIF
      
      END Subroutine ORDER_FTOD_IA
      
      SUBROUTINE ORDER_T2(WDES,KA,KB,XY)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: KA,KB,MNCC, ATPRO
         GDEFS :: XY(:,:,:,:)
         INTEGER :: ORDERLIST(PROCS,2),tag,ierror,myrank
         INTEGER :: req(PROCS),recreq
         INTEGER :: status(MPI_STATUS_SIZE, PROCS)
         
         tag=1
         call BLACS_PINFO(ME,PROCS)
         CALL MPI_COMM_RANK(WDES%comm, myrank, ierror)
         me=myrank
         ORDERLIST=-1
         IF (ME==0) THEN
            ORDERLIST(1,1)=KA
	    ORDERLIST(1,2)=KB
            DO ATPRO=1,PROCS-1
               !write(*,*)'Im:',ME,'and wait for ATPRO to send',ATPRO
               CALL MPI_RECV( ORDERLIST(ATPRO+1,1), 1, MPI_INTEGER, ATPRO, tag, &
                  WDES%COMM, status(1,ATPRO), ierror) !req(ATPRO), ierror )
!               write(*,*)'I and wait for ATPRO to send',ATPRO
               !CALL MPI_WAIT(req(ATPRO), status(1,ATPRO), ierror)
            ENDDO
         ELSE
            !write(*,*)'Im:',ME,'and send to 0'
            call MPI_ISEND(KA, 1, MPI_INTEGER, 0, tag, &
                 WDES%COMM,req(ME+1), ierror)
            CALL MPI_WAIT(req(ME+1), status(1,ME+1), ierror)
         ENDIF
	 IF (ME==0) THEN
            ORDERLIST(1,1)=KA
	    ORDERLIST(1,2)=KB
            DO ATPRO=1,PROCS-1
               !write(*,*)'Im:',ME,'and wait for ATPRO to send',ATPRO
               CALL MPI_RECV( ORDERLIST(ATPRO+1,2), 1, MPI_INTEGER, ATPRO, tag, &
                  WDES%COMM, status(1,ATPRO), ierror) !req(ATPRO), ierror )
!               write(*,*)'I and wait for ATPRO to send',ATPRO
               !CALL MPI_WAIT(req(ATPRO), status(1,ATPRO), ierror)
            ENDDO
         ELSE
            !write(*,*)'Im:',ME,'and send to 0'
            call MPI_ISEND(KB, 1, MPI_INTEGER, 0, tag, &
                 WDES%COMM,req(ME+1), ierror)
            CALL MPI_WAIT(req(ME+1), status(1,ME+1), ierror)
         ENDIF

         !write(*,*)'test3,me',me,ORDERLIST
         
         IF (ME==0) THEN
            CALL IGEBS2D(CONTXT_COLS, 'All', 'i-ring', PROCS*2, 1, ORDERLIST(1,1), PROCS*2)
         ELSE
            CALL IGEBR2D(CONTXT_COLS, 'All', 'i-ring', PROCS*2, 1, ORDERLIST(1,1), PROCS*2, 0,0)
         ENDIF
         !write(*,*)'test',ORDERLIST

         DO ATPRO=1,PROCS
            !write(*,*)'atpro:',ATPRO
            IF ((PROCS_KPTS(ORDERLIST(ATPRO,1))==ME) .and. ((ATPRO-1)/=ME)) THEN
               !write(*,*)'test2',ORDERLIST
               !write(*,*)'Im:',ME,'and send to',ATPRO
               call MPI_ISEND(T2(1,1,1,1,ORDERLIST(ATPRO,2),MKPTS_KPTS(ORDERLIST(ATPRO,1))), &
                           (PROCS*WDES%NBANDS-VBMAX)*(PROCS*WDES%NBANDS-VBMAX)*VBMAX*VBMAX, MPI_COMPLEX8, ATPRO-1, tag, &
                           WDES%COMM, req(ME+1), ierror)
               CALL MPI_WAIT(req(ME+1), status(1,ME+1), ierror)
            ENDIF
            IF ((PROCS_KPTS(KA)/=ME) .and. ((ATPRO-1)==ME)) THEN
               call MPI_RECV(XY(1,1,1,1), &
                           (PROCS*WDES%NBANDS-VBMAX)*(PROCS*WDES%NBANDS-VBMAX)*VBMAX*VBMAX, MPI_COMPLEX8, PROCS_KPTS(KA), tag, &
                           WDES%COMM, status(1,ME+1), ierror) !, recreq, ierror)
            ENDIF
         ENDDO
         IF ((PROCS_KPTS(KA)==ME)) THEN
            XY(:,:,:,:)=T2(:,:,:,:,KB,MKPTS_KPTS(KA)) !FTOD_PW_IA(:,:,:,:,MKPTS_KPTS(KI),MIN(NCC,MNCC))
         ENDIF
         !write(*,*)'finish',orderlist
         !IF (PROCS_KPTS(KI)/=ME) CALL MPI_WAIT(recreq, status(1,ME), ierror)

!           PW_IA_TMP(:,:,:,:)=zero
!           OC_IA_TMP(:,:,:,:)=zero

!            IF (PROCS_KPTS(KI)==ME) THEN
!               CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, (PROCS*WDES%NBANDS-VBMAX)*VBMAX*WDES%NKPTS, FTOD_PW_IA(1,1,1,1,MKPTS_KPTS(KI),MIN(NCC,MNCC)), NGVECTOR)
!               PW_IA_TMP(:,:,:,:)=FTOD_PW_IA(:,:,:,:,MKPTS_KPTS(KI),MIN(NCC,MNCC))
!            ELSE
!               CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, (PROCS*WDES%NBANDS-VBMAX)*VBMAX*WDES%NKPTS, PW_IA_TMP(1,1,1,1), NGVECTOR,0,PROCS_KPTS(KI))
!            ENDIF
            
!            IF (ASSOCIATED(H)) THEN
!               IF (PROCS_KPTS(KI)==ME) THEN
!                  CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, (PROCS*WDES%NBANDS-VBMAX)*VBMAX*WDES%NKPTS, FTOD_OC_IA(1,1,1,1,MKPTS_KPTS(KI),MNCC), NHVECTOR)
!                  OC_IA_TMP(:,:,:,:)=FTOD_OC_IA(:,:,:,:,MKPTS_KPTS(KI),MNCC)
!               ELSE
!                  CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, (PROCS*WDES%NBANDS-VBMAX)*VBMAX*WDES%NKPTS, OC_IA_TMP(1,1,1,1), NHVECTOR,0,PROCS_KPTS(KI))
!               ENDIF
            
!            ENDIF
      
      END Subroutine ORDER_T2

      Subroutine BCAST2ALL_FTOD_IA(WDES,KI,MNCC,PW_IA_TMP,OC_IA_TMP)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: KI,MNCC
         GDEF :: PW_IA_TMP(:,:,:,:), OC_IA_TMP(:,:,:,:)
            
            PW_IA_TMP(:,:,:,:)=zero
            OC_IA_TMP(:,:,:,:)=zero
            call BLACS_PINFO(ME,PROCS)

            IF (PROCS_KPTS(KI)==ME) THEN
               CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, (PROCS*WDES%NBANDS-VBMAX)*VBMAX*WDES%NKPTS, FTOD_PW_IA(1,1,1,1,MKPTS_KPTS(KI),MIN(NCC,MNCC)), NGVECTOR)
               PW_IA_TMP(:,:,:,:)=FTOD_PW_IA(:,:,:,:,MKPTS_KPTS(KI),MIN(NCC,MNCC))
            ELSE
               CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, (PROCS*WDES%NBANDS-VBMAX)*VBMAX*WDES%NKPTS, PW_IA_TMP(1,1,1,1), NGVECTOR,0,PROCS_KPTS(KI))
            ENDIF
            
            IF (ASSOCIATED(H)) THEN
            
               IF (PROCS_KPTS(KI)==ME) THEN
!#ifdef gammareal
!                  CALL DGEBS2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, (PROCS*WDES%NBANDS-VBMAX)*VBMAX*WDES%NKPTS, FTOD_OC_IA(1,1,1,1,MKPTS_KPTS(KI),MNCC), NHVECTOR)
!                  OC_IA_TMP(:,:,:,:)=FTOD_OC_IA(:,:,:,:,MKPTS_KPTS(KI),MNCC)
!#else
                  CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, (PROCS*WDES%NBANDS-VBMAX)*VBMAX*WDES%NKPTS, FTOD_OC_IA(1,1,1,1,MKPTS_KPTS(KI),MNCC), NHVECTOR)
                  OC_IA_TMP(:,:,:,:)=FTOD_OC_IA(:,:,:,:,MKPTS_KPTS(KI),MNCC)
!#endif
               ELSE
!#ifdef gammareal
!                  CALL DGEBR2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, (PROCS*WDES%NBANDS-VBMAX)*VBMAX*WDES%NKPTS, OC_IA_TMP(1,1,1,1), NHVECTOR,0,PROCS_KPTS(KI))
!#else
                  CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, (PROCS*WDES%NBANDS-VBMAX)*VBMAX*WDES%NKPTS, OC_IA_TMP(1,1,1,1), NHVECTOR,0,PROCS_KPTS(KI))
!#endif
               ENDIF
            
            ENDIF
      
      END Subroutine BCAST2ALL_FTOD_IA

      Subroutine BCAST2ALL_FTOD_AI(WDES,KA,MNCC,PW_AI_TMP,OC_AI_TMP)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: KA,MNCC
         GDEF :: PW_AI_TMP(:,:,:,:), OC_AI_TMP(:,:,:,:)
            
            PW_AI_TMP(:,:,:,:)=zero
            OC_AI_TMP(:,:,:,:)=zero
            call BLACS_PINFO(ME,PROCS)

            IF (PROCS_KPTS(KA)==ME) THEN
               CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, VBMAX*(PROCS*WDES%NBANDS-VBMAX)*WDES%NKPTS, FTOD_PW_AI(1,1,1,1,MKPTS_KPTS(KA),MIN(NCC,MNCC)), NGVECTOR)
               PW_AI_TMP(:,:,:,:)=FTOD_PW_AI(:,:,:,:,MKPTS_KPTS(KA),MIN(NCC,MNCC))
            ELSE
               CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, VBMAX*(PROCS*WDES%NBANDS-VBMAX)*WDES%NKPTS, PW_AI_TMP(1,1,1,1), NGVECTOR,0,PROCS_KPTS(KA))
            ENDIF
            
            IF (ASSOCIATED(H)) THEN
            
               IF (PROCS_KPTS(KA)==ME) THEN
!#ifdef gammareal
!                  CALL DGEBS2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, VBMAX*(PROCS*WDES%NBANDS-VBMAX)*WDES%NKPTS, FTOD_OC_AI(1,1,1,1,MKPTS_KPTS(KA),MNCC), NHVECTOR)
!                  OC_AI_TMP(:,:,:,:)=FTOD_OC_AI(:,:,:,:,MKPTS_KPTS(KA),MNCC)
!#else
                  CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, VBMAX*(PROCS*WDES%NBANDS-VBMAX)*WDES%NKPTS, FTOD_OC_AI(1,1,1,1,MKPTS_KPTS(KA),MNCC), NHVECTOR)
                  OC_AI_TMP(:,:,:,:)=FTOD_OC_AI(:,:,:,:,MKPTS_KPTS(KA),MNCC)
!#endif
               ELSE
!#ifdef gammareal
!                  CALL DGEBR2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, VBMAX*(PROCS*WDES%NBANDS-VBMAX)*WDES%NKPTS, OC_AI_TMP(1,1,1,1), NHVECTOR,0,PROCS_KPTS(KA))
!#else
                  CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, VBMAX*(PROCS*WDES%NBANDS-VBMAX)*WDES%NKPTS, OC_AI_TMP(1,1,1,1), NHVECTOR,0,PROCS_KPTS(KA))
!#endif
               ENDIF
            
            ENDIF
      
      END Subroutine BCAST2ALL_FTOD_AI
      
      SUBROUTINE BCAST2ALL_CHI_GAI(WDES,KA,MNCC,GAI_TMP,OC_GAI_TMP)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: KA,MNCC
         GDEF :: GAI_TMP(:,:,:,:)
         GDEF, OPTIONAL :: OC_GAI_TMP(:,:,:,:)
         integer :: time_array1(8),time_array2(8),ems
            
            call date_and_time (values=time_array1)
            
            call BLACS_PINFO(ME,PROCS)

            IF (PROCS_KPTS(KA)==ME) THEN
               CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, (PROCS*WDES%NBANDS-VBMAX)*(VBMAX)*WDES%NKPTS, CHI_GAI(1,1,1,1,MKPTS_KPTS(KA),MIN(NCC,MNCC)), NGVECTOR)
               GAI_TMP(:,:,:,:)=CHI_GAI(:,:,:,:,MKPTS_KPTS(KA),MIN(NCC,MNCC))
            ELSE
               CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NGVECTOR, (PROCS*WDES%NBANDS-VBMAX)*(VBMAX)*WDES%NKPTS,GAI_TMP(1,1,1,1), NGVECTOR,0,PROCS_KPTS(KA))
            ENDIF
            
            IF (ASSOCIATED(H)) THEN
            
               IF (PROCS_KPTS(KA)==ME) THEN
!                  CALL ZGEBS2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, (PROCS*WDES%NBANDS-VBMAX)*(PROCS*WDES%NBANDS-VBMAX)*WDES%NKPTS, FTOD_OC_AB(1,1,1,1,MKPTS_KPTS(KA),MNCC), NHVECTOR)
!                  OC_AB_TMP(:,:,:,:)=FTOD_OC_AB(:,:,:,:,MKPTS_KPTS(KA),MNCC)
               ELSE
!                  CALL ZGEBR2D(CONTXT_COLS, 'All', 'i-ring', NHVECTOR, (PROCS*WDES%NBANDS-VBMAX)*(PROCS*WDES%NBANDS-VBMAX)*WDES%NKPTS, OC_AB_TMP(1,1,1,1), NHVECTOR,0,PROCS_KPTS(KA))
               ENDIF
            
            ENDIF
            
            call date_and_time (values=time_array2)
            ems=((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
                time_array2(8)-time_array1(8))
            TETS=TETS+ems
            
      END SUBROUTINE BCAST2ALL_CHI_GAI

      Subroutine BCAST2ALL_T2(WDES,KA,T2_MTMP)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: KA, CHIL, CHIM
         GDEFS :: T2_MTMP(:,:,:,:,:)
         integer :: time_array1(8),time_array2(8),ems

            call BLACS_PINFO(ME,PROCS)
            call date_and_time (values=time_array1)

            IF (PROCS_KPTS(KA)==ME) THEN
               CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (WDES%NB_TOT-VBMAX)*(WDES%NB_TOT-VBMAX), (VBMAX*VBMAX*WDES%NKPTS), &
                  T2(1,1,1,1,1,MKPTS_KPTS(KA)),(WDES%NB_TOT-VBMAX)*(WDES%NB_TOT-VBMAX))
               T2_MTMP(:,:,:,:,:)=T2(:,:,:,:,:,MKPTS_KPTS(KA))
            ELSE
               CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (WDES%NB_TOT-VBMAX)*(WDES%NB_TOT-VBMAX), (VBMAX*VBMAX*WDES%NKPTS), &
                  T2_MTMP(1,1,1,1,1), (WDES%NB_TOT-VBMAX)*(WDES%NB_TOT-VBMAX),0,PROCS_KPTS(KA))
            ENDIF

            call date_and_time (values=time_array2)
            ems=((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
                time_array2(8)-time_array1(8))
            TET2S=TET2S+ems

      END Subroutine BCAST2ALL_T2

!      Subroutine BCAST2ALL_T2_N(WDES,NBANDSGW,KA,T2_MTMP)
!         use mkpoints
!         use base
!         IMPLICIT NONE
!         TYPE(wavedes) WDES
!         INTEGER :: KA, CHIL, CHIM, NBANDSGW
!         GDEFS :: T2_MTMP(:,:,:,:,:)
!         integer :: time_array1(8),time_array2(8),ems

!            call BLACS_PINFO(ME,PROCS)
!            call date_and_time (values=time_array1)

!            IF (PROCS_KPTS(KA)==ME) THEN
!               CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (NBANDSGW-VBMAX)*(PROCS*WDES%NBANDS-VBMAX), (VBMAX*VBMAX*WDES%NKPTS*WDES%NKPTS), &
!                  T2_N(1,1,1,1,1,MKPTS_KPTS(KA)), (NBANDSGW-VBMAX)*(PROCS*WDES%NBANDS-VBMAX))
!               T2_MTMP(:,:,:,:,:)=T2_N(:,:,:,:,:,MKPTS_KPTS(KA))
!            ELSE
!               CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (NBANDSGW-VBMAX)*(PROCS*WDES%NBANDS-VBMAX),(VBMAX*VBMAX*WDES%NKPTS*WDES%NKPTS), &
!                  T2_MTMP(1,1,1,1,1), (NBANDSGW-VBMAX)*(PROCS*WDES%NBANDS-VBMAX),0,PROCS_KPTS(KA))
!            ENDIF

!            call date_and_time (values=time_array2)
!            ems=((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
!                time_array2(8)-time_array1(8))
!            TET2S=TET2S+ems

!      END Subroutine BCAST2ALL_T2_N
      
      Subroutine BCAST2ALL_T2_KJ_KA(WDES,KJ,KA,T2_MTMP)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: KA, CHIL, CHIM, KJ
         GDEFS :: T2_MTMP(:,:,:,:,:)
         integer :: time_array1(8),time_array2(8),ems

            call BLACS_PINFO(ME,PROCS)
            call date_and_time (values=time_array1)

            IF (PROCS_KPTS(KA)==ME) THEN
!               CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (PROCS*WDES%NBANDS-VBMAX), (PROCS*WDES%NBANDS-VBMAX)*(VBMAX*VBMAX*WDES%NKPTS), &
!                 T2(1,1,1,1,1,KJ,MKPTS_KPTS(KA)), (PROCS*WDES%NBANDS-VBMAX))
!               T2_MTMP(:,:,:,:,:)=T2(:,:,:,:,:,KJ,MKPTS_KPTS(KA))
            ELSE
!               CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (PROCS*WDES%NBANDS-VBMAX), (PROCS*WDES%NBANDS-VBMAX)*(VBMAX*VBMAX*WDES%NKPTS), &
!                  T2_MTMP(1,1,1,1,1), (PROCS*WDES%NBANDS-VBMAX),0,PROCS_KPTS(KA))
            ENDIF

            call date_and_time (values=time_array2)
            ems=((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
                time_array2(8)-time_array1(8))
            TET2S=TET2S+ems

      END Subroutine BCAST2ALL_T2_KJ_KA

      Subroutine BCAST2ALL_T2_KI_KJ_KA(W,WDES,KI,KJ,KA,T2_MTMP,NBANDSGW)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(wavespin) W
         INTEGER :: KA, CHIL, CHIM, KJ, KI,KB,NI,NJ,NBANDSGW,KQ,NA,NB
         GDEFS :: T2_MTMP(:,:,:,:)
         integer :: time_array1(8),time_array2(8),ems

            call BLACS_PINFO(ME,PROCS)
            call date_and_time (values=time_array1)

            IF (NBANDSGW/=WDES%NB_TOT) THEN

            KQ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KA)- &
               WDES%VKPT(:,KI),KPOINTS_FULL)
            KB=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KJ)- &
               WDES%VKPT(:,KQ),KPOINTS_FULL)

            CALL BCAST2ALL_FTOD_AI(WDES,KB,2,PW_AI_TMP,OC_AI_TMP)

            ENDIF

            IF (PROCS_KPTS(KA)==ME) THEN
            
               VV_S=(0.0_qs,0.0_qs)
               VVOO_S=(0.0_qs,0.0_qs)
               DO NI=1,VBMAX
               DO NJ=1,VBMAX
                  CALL WRITE_T2_2_VV(WDES,NBANDSGW,VV_S,T2(1,NI,NJ,KI,KJ,MKPTS_KPTS(KA)))
                  VVOO_S(:,:,NI,NJ)=VVOO_S(:,:,NI,NJ)+VV_S(:,:)
               ENDDO
               ENDDO
               
               DO NA=(NBANDSGW-VBMAX)+1,(WDES%NB_TOT-VBMAX)
               DO NB=(NBANDSGW-VBMAX)+1,(WDES%NB_TOT-VBMAX)
                  CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
                    FTOD_PW_AI(1,1,NA,KQ,MKPTS_KPTS(KA),1),FTOD_OC_AI(1,1,NA,KQ,MKPTS_KPTS(KA),1),&
                    VBMAX,&
                    PW_AI_TMP(1,1,NB,KQ),OC_AI_TMP(1,1,NB,KQ), &
                    VBMAX,TE4O(1,1),zero)
                  DO NI=1,VBMAX
                  DO NJ=1,VBMAX
                     VVOO_S(NA,NB,NI,NJ)=VVOO_S(NA,NB,NI,NJ)+(TE4O(NI,NJ))/&
                        (W%CELTOT(NJ,KJ,1)+ &
                        W%CELTOT(NI,KI,1)- &
                        W%CELTOT(NB+VBMAX,KB,1)- &
                        W%CELTOT(NA+VBMAX,KA,1))  !bcjk
                  ENDDO
                  ENDDO
               ENDDO
               ENDDO


               CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (PROCS*WDES%NBANDS-VBMAX), (PROCS*WDES%NBANDS-VBMAX)*(VBMAX*VBMAX), &
                  VVOO_S(1,1,1,1), (PROCS*WDES%NBANDS-VBMAX))
               T2_MTMP=VVOO_S
            ELSE
               CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (PROCS*WDES%NBANDS-VBMAX), (PROCS*WDES%NBANDS-VBMAX)*(VBMAX*VBMAX), &
                  T2_MTMP(1,1,1,1), (PROCS*WDES%NBANDS-VBMAX),0,PROCS_KPTS(KA))
            ENDIF

            call date_and_time (values=time_array2)
            ems=((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
                time_array2(8)-time_array1(8))
            TET2S=TET2S+ems

      END Subroutine BCAST2ALL_T2_KI_KJ_KA

      Subroutine BCAST2ALL_AS_T2_KI_KJ_KA(W,WDES,KJ,KK,KB,T2_MTMP,NBANDSGW)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(wavespin) W
         INTEGER :: KB, CHIL, CHIM, KJ, KK,NK,NJ,NC,NB,KC,NBANDSGW,KQ,KQ_
         GDEFS :: T2_MTMP(:,:,:,:)
         integer :: time_array1(8),time_array2(8),ems

            call BLACS_PINFO(ME,PROCS)
            call date_and_time (values=time_array1)
            
            IF (NBANDSGW/=WDES%NB_TOT) THEN
               KQ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KB)- &
                  WDES%VKPT(:,KJ),KPOINTS_FULL)
               KC=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KK)- &
                  WDES%VKPT(:,KQ),KPOINTS_FULL)
               KQ_=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KB)- &
                  WDES%VKPT(:,KK),KPOINTS_FULL)
               CALL BCAST2ALL_FTOD_AI(WDES,KC,2,PW_AI_TMP,OC_AI_TMP)
               
            IF (PROCS_KPTS(KB)==ME) THEN
               T2_MTMP=(0.0_qs,0.0_qs)
               IF (.NOT. RING) THEN
               DO NB=(NBANDSGW-VBMAX)+1,(WDES%NB_TOT-VBMAX)
               DO NC=(NBANDSGW-VBMAX)+1,(WDES%NB_TOT-VBMAX)
                  CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
                    PW_AI_TMP(1,1,NC,KQ_),OC_AI_TMP(1,1,NC,KQ_), &
                    VBMAX,&
                    FTOD_PW_AI(1,1,NB,KQ_,MKPTS_KPTS(KB),1),FTOD_OC_AI(1,1,NB,KQ_,MKPTS_KPTS(KB),1),&
                    VBMAX,TE4O(1,1),zero)
                  T2_MTMP(NB,NC,:,:)=-CONJG(TE4O(:,:)) !bcjk
               ENDDO
               ENDDO
               VV_S=(0.0_qs,0.0_qs)
               DO NK=1,VBMAX
               DO NJ=1,VBMAX
                  CALL WRITE_T2_2_VV(WDES,NBANDSGW,VV_S,T2(1,NK,NJ,KK,KJ,MKPTS_KPTS(KB)))
                  T2_MTMP(:,:,NJ,NK)=T2_MTMP(:,:,NJ,NK)-VV_S(:,:)
               ENDDO
               ENDDO
               ENDIF

               VV_S=(0.0_qs,0.0_qs)
               DO NK=1,VBMAX
               DO NJ=1,VBMAX
                  CALL WRITE_T2_2_VV(WDES,NBANDSGW,VV_S,T2(1,NJ,NK,KJ,KK,MKPTS_KPTS(KB)))
                  T2_MTMP(:,:,NJ,NK)=T2_MTMP(:,:,NJ,NK)+(2.0_qs,0.0_qs)*VV_S(:,:)
               ENDDO
               ENDDO
               
               DO NC=(NBANDSGW-VBMAX)+1,(WDES%NB_TOT-VBMAX)
               DO NB=(NBANDSGW-VBMAX)+1,(WDES%NB_TOT-VBMAX)
                  CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
                    FTOD_PW_AI(1,1,NB,KQ,MKPTS_KPTS(KB),1),FTOD_OC_AI(1,1,NB,KQ,MKPTS_KPTS(KB),1),&
                    VBMAX,&
                    PW_AI_TMP(1,1,NC,KQ),OC_AI_TMP(1,1,NC,KQ), &
                    VBMAX,TE4O(1,1),zero)
                  DO NK=1,VBMAX
                  DO NJ=1,VBMAX
                     T2_MTMP(NB,NC,NJ,NK)=(T2_MTMP(NB,NC,NJ,NK)+((2.0_qs,0.0_qs)*TE4O(NJ,NK)))/&
                        (W%CELTOT(NJ,KJ,1)+ &
                        W%CELTOT(NK,KK,1)- &
                        W%CELTOT(NB+VBMAX,KB,1)- &
                        W%CELTOT(NC+VBMAX,KC,1))  !bcjk
                  ENDDO
                  ENDDO
               ENDDO
               ENDDO

            ENDIF
            ENDIF



            IF (PROCS_KPTS(KB)==ME) THEN
               CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (WDES%NB_TOT-VBMAX)*(WDES%NB_TOT-VBMAX),(VBMAX*VBMAX), &
                  T2_MTMP(1,1,1,1), (WDES%NB_TOT-VBMAX)*(WDES%NB_TOT-VBMAX))
                !T2_MTMP=T2_MTMP
            ELSE
               CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (WDES%NB_TOT-VBMAX)*(WDES%NB_TOT-VBMAX),(VBMAX*VBMAX), &
                  T2_MTMP(1,1,1,1), (WDES%NB_TOT-VBMAX)*(WDES%NB_TOT-VBMAX),0,PROCS_KPTS(KB))
            ENDIF

            call date_and_time (values=time_array2)
            ems=((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
                time_array2(8)-time_array1(8))
            TET2S=TET2S+ems

      END Subroutine BCAST2ALL_AS_T2_KI_KJ_KA

      Subroutine BCAST2ALL_AS_T2_KJ_KA(WDES,KJ,KA,T2_MTMP)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: KA, CHIL, CHIM, KJ, KK, KB, NJ, NK
         GDEFS :: T2_MTMP(:,:,:,:,:)
         integer :: time_array1(8),time_array2(8),ems

            call BLACS_PINFO(ME,PROCS)
            call date_and_time (values=time_array1)

            IF (PROCS_KPTS(KA)==ME) THEN
               DO KK=1,WDES%NKPTS
               DO NK=1,VBMAX
               DO NJ=1,VBMAX
!                 IF (.NOT. RING) T2_MTMP(:,:,NK,NJ,KK)=-T2(:,:,NK,NJ,KK,KJ,MKPTS_KPTS(KA))+(2.0_qs,0.0_qs)*T2(:,:,NJ,NK,KJ,KK,MKPTS_KPTS(KA))
!                 IF (RING) T2_MTMP(:,:,NK,NJ,KK)=(2.0_qs,0.0_qs)*T2(:,:,NJ,NK,KJ,KK,MKPTS_KPTS(KA))
               ENDDO
               ENDDO
               ENDDO
            
!               CALL CGEBS2D(CONTXT_COLS, 'All', 'i-ring', (PROCS*WDES%NBANDS-VBMAX), (PROCS*WDES%NBANDS-VBMAX)*(VBMAX*VBMAX*WDES%NKPTS), &
!                 T2_MTMP(1,1,1,1,1), (PROCS*WDES%NBANDS-VBMAX))
!               !T2_MTMP(:,:,:,:,:)=T2(:,:,:,:,:,KJ,MKPTS_KPTS(KA))
            ELSE
!               CALL CGEBR2D(CONTXT_COLS, 'All', 'i-ring', (PROCS*WDES%NBANDS-VBMAX), (PROCS*WDES%NBANDS-VBMAX)*(VBMAX*VBMAX*WDES%NKPTS), &
!                  T2_MTMP(1,1,1,1,1), (PROCS*WDES%NBANDS-VBMAX),0,PROCS_KPTS(KA))
            ENDIF

            call date_and_time (values=time_array2)
            ems=((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
                time_array2(8)-time_array1(8))
            TET2S=TET2S+ems

      END Subroutine BCAST2ALL_AS_T2_KJ_KA
      

      SUBROUTINE CONTR_FTOD(NGV,NHV,PW_XY,OC_XY,LNXY,PW_XY_,OC_XY_, &
                            LNXY_,XY2E4ORB,beta)
      IMPLICIT NONE
      INTEGER :: NGV,NHV,LNXY,LNXY_
      GDEF :: XY2E4ORB,PW_XY,OC_XY, &
              PW_XY_,OC_XY_, beta
      integer :: time_array1(8),time_array2(8),ems
!#ifdef gammareal
      call date_and_time (values=time_array1)
!#else
            CALL ZGEMM(trans,'n',(LNXY),(LNXY_),&
                 (NGV),(1._q,0._q),PW_XY,(NGV),&
                 PW_XY_,(NGV),&
                 beta, XY2E4ORB,(LNXY))
!#endif

            IF (ASSOCIATED(H)) THEN
!#ifdef gammareal

!#else
            CALL ZGEMM(trans,'n',(LNXY),(LNXY_),&
                 (NHV),(1._q,0._q),OC_XY,(NHV),&
                 OC_XY_,(NHV),&
                 (1._q,0._q), XY2E4ORB,(LNXY))
!#endif
            ENDIF
            
      call date_and_time (values=time_array2)
      ems=((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
                time_array2(8)-time_array1(8))
      contr_ftod_time=contr_ftod_time+ems

      END SUBROUTINE
      
      SUBROUTINE SORT_O1V1O2V2_V2V1O2O1(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEF :: XY(:,:,:,:),XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2
         
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NO1=1,VBMAX
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NO2=1,VBMAX
            XY2(NV1,NV2,NO1,NO2)=XY(NO2,NV2,NO1,NV1)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
          
      END SUBROUTINE SORT_O1V1O2V2_V2V1O2O1

      SUBROUTINE SORT_O1V1O2V2_V1O2V2O1(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEF :: XY(:,:,:,:),XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2
         
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NO1=1,VBMAX
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NO2=1,VBMAX
            XY2(NV1,NO2,NV2,NO1)=XY(NO1,NV1,NO2,NV2)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         
      END SUBROUTINE SORT_O1V1O2V2_V1O2V2O1

      SUBROUTINE SORT_V1V2O1O2_O1O2V1V2(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:,:),XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2
         
         DO NO2=1,VBMAX
         DO NO1=1,VBMAX
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)
            XY2(NO1,NO2,NV1,NV2)=XY(NV1,NV2,NO1,NO2)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         
      END SUBROUTINE SORT_V1V2O1O2_O1O2V1V2

      SUBROUTINE SORT_V1V2O1O2_O1O2V2V1(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:,:),XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2

         DO NO2=1,VBMAX
         DO NO1=1,VBMAX
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)
            XY2(NO1,NO2,NV2,NV1)=XY(NV1,NV2,NO1,NO2)
         ENDDO
         ENDDO
         ENDDO
         ENDDO

      END SUBROUTINE SORT_V1V2O1O2_O1O2V2V1

      SUBROUTINE SORT_V1V2O1O2_O2O1V2V1(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:,:),XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2

         DO NO2=1,VBMAX
         DO NO1=1,VBMAX
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)
            XY2(NO2,NO1,NV2,NV1)=XY(NV1,NV2,NO1,NO2)
         ENDDO
         ENDDO
         ENDDO
         ENDDO

      END SUBROUTINE SORT_V1V2O1O2_O2O1V2V1
      
      SUBROUTINE SORT_O1O2V1V2_V2O2O1V1(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:,:),XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2

         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NO2=1,VBMAX
         DO NO1=1,VBMAX
            XY2(NV2,NO2,NO1,NV1)=XY(NO1,NO2,NV1,NV2)
         ENDDO
         ENDDO
         ENDDO
         ENDDO

      END SUBROUTINE SORT_O1O2V1V2_V2O2O1V1

      SUBROUTINE SORT_O1O2V1V2_V1O1V2O2(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:,:),XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2

         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NO2=1,VBMAX
         DO NO1=1,VBMAX
            XY2(NV1,NO1,NV2,NO2)=XY(NO1,NO2,NV1,NV2)
         ENDDO
         ENDDO
         ENDDO
         ENDDO

      END SUBROUTINE SORT_O1O2V1V2_V1O1V2O2
      
      SUBROUTINE SORT_V1V2_V2V1(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:),XY2(:,:)
         INTEGER :: NV1,NV2

         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)
            XY2(NV2,NV1)=XY(NV1,NV2)
         ENDDO
         ENDDO
         
      END SUBROUTINE SORT_V1V2_V2V1

      SUBROUTINE SORT_V1O1_O1V1(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:),XY2(:,:)
         INTEGER :: NV1,NO1

         DO NO1=1,(VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)
            XY2(NO1,NV1)=XY(NV1,NO1)
         ENDDO
         ENDDO
         
      END SUBROUTINE SORT_V1O1_O1V1
      
      SUBROUTINE SORT_O1V1_V1O1(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:),XY2(:,:)
         INTEGER :: NV1,NO1

         DO NO1=1,(VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)
            XY2(NV1,NO1)=XY(NO1,NV1)
         ENDDO
         ENDDO
         
      END SUBROUTINE SORT_O1V1_V1O1

      SUBROUTINE SORT_O1V1V2O2_V1V2O1O2(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEF :: XY(:,:,:,:)
         GDEFS :: XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2
         
         DO NO2=1,VBMAX
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NO1=1,VBMAX
            XY2(NV1,NV2,NO1,NO2)=XY(NO1,NV1,NV2,NO2)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         
      END SUBROUTINE SORT_O1V1V2O2_V1V2O1O2
      
      SUBROUTINE SORT_V2V1O1O2_V1V2O1O2(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:,:)
         GDEFS :: XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2
         
         DO NO1=1,VBMAX
         DO NO2=1,VBMAX
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
            XY2(NV1,NV2,NO1,NO2)=XY(NV2,NV1,NO1,NO2)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         
      END SUBROUTINE SORT_V2V1O1O2_V1V2O1O2

      SUBROUTINE SORT_V2V1O1O2_V2O2V1O1(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEF :: XY(:,:,:,:)
         GDEFS :: XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2
         
         DO NO1=1,VBMAX
         DO NO2=1,VBMAX
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
            XY2(NV2,NO2,NV1,NO1)=XY(NV2,NV1,NO1,NO2)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         
      END SUBROUTINE SORT_V2V1O1O2_V2O2V1O1
      
      SUBROUTINE SORT_V1V2O1O2_V2O1V1O2(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:,:)
         GDEFS :: XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2
         
         DO NO1=1,VBMAX
         DO NO2=1,VBMAX
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)      
            XY2(NV2,NO1,NV1,NO2)=XY(NV1,NV2,NO1,NO2)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         
      END SUBROUTINE SORT_V1V2O1O2_V2O1V1O2

      SUBROUTINE SORT_V1V2O1O2_V2O1V1O2_KTMP(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:,:,:)
         GDEFS :: XY2(:,:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2,KK
         
         DO KK=1,WDES%NKPTS
         DO NO1=1,VBMAX
         DO NO2=1,VBMAX
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)      
            XY2(NV2,NO1,NV1,NO2,KK)=XY(NV1,NV2,NO1,NO2,KK)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         
      END SUBROUTINE SORT_V1V2O1O2_V2O1V1O2_KTMP
      
      SUBROUTINE SORT_V1V2O1O2_V2V1O2O1(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:,:)
         GDEFS :: XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2,KK
         
         DO NO1=1,VBMAX
         DO NO2=1,VBMAX
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)      
            XY2(NV2,NV1,NO2,NO1)=XY(NV1,NV2,NO1,NO2)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         
      END SUBROUTINE SORT_V1V2O1O2_V2V1O2O1
      
      SUBROUTINE SORT_T2_V1V2O1O2_V2O1V1O2(WDES,KK,KI,KA,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2,KK, KI, KA
         
         DO NO1=1,VBMAX
         DO NO2=1,VBMAX
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)      
!            XY2(NV2,NO1,NV1,NO2)=T2(NV1,NV2,NO1,NO2,KK,KI,KA)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         
      END SUBROUTINE SORT_T2_V1V2O1O2_V2O1V1O2
      
      SUBROUTINE SORT_AS_T2_V1V2O1O2_V2O1V1O2(WDES,KK,KI,KA,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2,KK, KI, KA
         
         DO NO1=1,VBMAX
         DO NO2=1,VBMAX
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)      
!            XY2(NV2,NO1,NV1,NO2)=2.0_qs*T2(NV1,NV2,NO2,NO1,KI,KK,KA)- &
!              T2(NV1,NV2,NO1,NO2,KK,KI,KA)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         
      END SUBROUTINE SORT_AS_T2_V1V2O1O2_V2O1V1O2

      SUBROUTINE SORT_V1O1V2O2_V1V2O1O2(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:,:)
         GDEFS :: XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2
         
             DO NO1=1,VBMAX
             DO NO2=1,VBMAX
             DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
             DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)  
            XY2(NV1,NV2,NO1,NO2)=XY(NV1,NO1,NV2,NO2)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         
      END SUBROUTINE SORT_V1O1V2O2_V1V2O1O2
      
      SUBROUTINE SORT_V1O1V2O2_V2O1V1O2(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:,:)
         GDEFS :: XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2
         
             DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)
             DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)  
            XY2(NV1,:,NV2,:)=XY(NV2,:,NV1,:)
         ENDDO
         ENDDO
         
      END SUBROUTINE SORT_V1O1V2O2_V2O1V1O2
      
      SUBROUTINE SORT_V1V2O1O2_V2O2V1O1(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:,:)
         GDEFS :: XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2
         integer :: time_array1(8),time_array2(8),ems
         
         call date_and_time (values=time_array1)
         DO NO1=1,VBMAX
         DO NO2=1,VBMAX
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)      
            XY2(NV2,NO2,NV1,NO1)=XY(NV1,NV2,NO1,NO2)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
         call date_and_time (values=time_array2)
            ems=((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
                time_array2(8)-time_array1(8))
         TETSORT=TETSORT+ems
            
      END SUBROUTINE SORT_V1V2O1O2_V2O2V1O1

      SUBROUTINE SORT_V1V2O1O2_V1O1V2O2(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:,:)
         GDEFS :: XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2
         integer :: time_array1(8),time_array2(8),ems
         
         call date_and_time (values=time_array1)
         DO NO1=1,VBMAX
         DO NO2=1,VBMAX
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
            XY2(:,NO1,NV2,NO2)=XY(:,NV2,NO1,NO2)
         ENDDO
         ENDDO
         ENDDO
         call date_and_time (values=time_array2)
            ems=((time_array2(6)*60+time_array2(7))*1000-(time_array1(6)*60+time_array1(7))*1000+&
                time_array2(8)-time_array1(8))
         TETSORT=TETSORT+ems
            
      END SUBROUTINE SORT_V1V2O1O2_V1O1V2O2

      SUBROUTINE SORT_V1V2O1O2_V1V2O2O1(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:,:)
         GDEFS :: XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2
         integer :: time_array1(8),time_array2(8),ems
         
         DO NO1=1,VBMAX
         DO NO2=1,VBMAX
            XY2(:,:,NO2,NO1)=XY(:,:,NO1,NO2)
         ENDDO
         ENDDO
            
      END SUBROUTINE SORT_V1V2O1O2_V1V2O2O1

      SUBROUTINE SORT_O1V1V2O2_V2O2V1O1(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEF :: XY(:,:,:,:)
         GDEFS :: XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2
         
         DO NO2=1,VBMAX
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)      
         DO NO1=1,VBMAX
            XY2(NV2,NO2,NV1,NO1)=XY(NO1,NV1,NV2,NO2)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      
      END SUBROUTINE SORT_O1V1V2O2_V2O2V1O1

      SUBROUTINE SORT_O2V1V2O1_V1O1V2O2(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:,:)
         GDEFS :: XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2
         
         DO NO2=1,VBMAX
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)      
         DO NO1=1,VBMAX
            XY2(NV1,NO1,NV2,NO2)=XY(NO2,NV1,NV2,NO1)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      
      END SUBROUTINE SORT_O2V1V2O1_V1O1V2O2
      
      SUBROUTINE SORT_V2O1V1O2_V1V2O1O2(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:,:)
         GDEFS :: XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2
         
         DO NO2=1,VBMAX
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)      
         DO NO1=1,VBMAX
            XY2(NV1,NV2,NO1,NO2)=XY(NV2,NO1,NV1,NO2)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      
      END SUBROUTINE SORT_V2O1V1O2_V1V2O1O2

      SUBROUTINE SORT_V1V2O1O2_O1V1V2O2(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:,:)
         GDEFS :: XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2

         DO NO2=1,VBMAX 
         DO NO1=1,VBMAX
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)
            XY2(NO1,NV1,NV2,NO2)=XY(NV1,NV2,NO1,NO2)
         ENDDO
         ENDDO
         ENDDO
         ENDDO

      END SUBROUTINE SORT_V1V2O1O2_O1V1V2O2

      SUBROUTINE SORT_V1V2O1O2_O2V1V2O1(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:,:)
         GDEFS :: XY2(:,:,:,:)
         INTEGER :: NV1,NO1,NV2,NO2

         DO NO2=1,VBMAX 
         DO NO1=1,VBMAX
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(PROCS*WDES%NBANDS-VBMAX)
            XY2(NO2,NV1,NV2,NO1)=XY(NV1,NV2,NO1,NO2)
         ENDDO
         ENDDO
         ENDDO
         ENDDO

      END SUBROUTINE SORT_V1V2O1O2_O2V1V2O1

      SUBROUTINE SORT_JLIK_KLIJ(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEF :: XY(:,:,:,:),XY2(:,:,:,:)
         INTEGER :: NK,NL,NI,NJ
      
         DO NK=1,VBMAX
         DO NL=1,VBMAX
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
            XY2(NK,NL,NI,NJ)=XY(NJ,NL,NI,NK)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      
      END SUBROUTINE SORT_JLIK_KLIJ

      SUBROUTINE WRITE_T2_2_VV(WDES,NBANDSGW,MVV,T)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: T(*)
         GDEFS :: MVV(:,:)
         INTEGER :: NV2,NV1,NC,NBANDSGW

         NC=0
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(NBANDSGW-VBMAX)
            NC=NC+1
            MVV(NV1,NV2)=T(NC)
         ENDDO
         ENDDO
         IF (NBANDSGW<PROCS*WDES%NBANDS) THEN
         DO NV2=1,(NBANDSGW-VBMAX)
         DO NV1=(NBANDSGW-VBMAX)+1,(PROCS*WDES%NBANDS-VBMAX)
            NC=NC+1
            MVV(NV1,NV2)=T(NC)
         ENDDO
         ENDDO
         ENDIF
         
       END SUBROUTINE WRITE_T2_2_VV

      SUBROUTINE WRITE_T2_N_2_VV(WDES,NBANDSGW,VV,T)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: T(*)
         GDEFS :: VV(:,:)
         INTEGER :: NV2,NV1,NC,NBANDSGW

         NC=0
         DO NV2=1,(PROCS*WDES%NBANDS-VBMAX)
         DO NV1=1,(NBANDSGW-VBMAX)
            NC=NC+1
            VV(NV1,NV2)=T(NC)
         ENDDO
         ENDDO
         
       END SUBROUTINE WRITE_T2_N_2_VV
       

      SUBROUTINE SORT_CO1O2_O1O2C(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:)
         GDEFS :: XY2(:,:,:)
         INTEGER :: NC,NO1,NO2

         DO NO2=1,VBMAX 
         DO NO1=1,VBMAX
         DO NC=1,NCLENGTH
            XY2(NO1,NO2,NC)=XY(NC,NO1,NO2)
         ENDDO
         ENDDO
         ENDDO

      END SUBROUTINE SORT_CO1O2_O1O2C

      SUBROUTINE SORT_CO1O2_CO2O1(WDES,XY,XY2)
         IMPLICIT NONE
         TYPE(wavedes) WDES
         GDEFS :: XY(:,:,:)
         GDEFS :: XY2(:,:,:)
         INTEGER :: NC,NO1,NO2

         DO NO2=1,VBMAX 
         DO NO1=1,VBMAX
            XY2(:,NO1,NO2)=XY(:,NO2,NO1)
         ENDDO
         ENDDO

      END SUBROUTINE SORT_CO1O2_CO2O1

      SUBROUTINE UPDATE_T2(WDES,WGW,W,IO,NBANDSGW,IT)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(wavedes) WGW
         TYPE(wavespin) W 
         TYPE (in_struct) IO
         !INTEGER :: MKB
         INTEGER :: KQ,KI,KJ,KQ_,NI,NJ,KA,KB,NA,NB,NBANDSGW,NC,NC_N,IT
         COMPLEX(q) :: E_CCSD
!	 COMPLEX(qs) :: mix
         
!         mix=(0.0_qs,0.0_qs)
!         IF ((RING) .and. (.not. (IT==1))) mix=(0.5_qs,0.0_qs)
	 

         KQ=GKQ

         DO KB=1,WDES%NKPTS
            DO KA=1,MY_NKPTS
               KI=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KQ)+ &
                  WDES%VKPT(:,KPTS_MKPTS(KA)),KPOINTS_FULL)
               KJ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KB)- &
                  WDES%VKPT(:,KQ),KPOINTS_FULL)
               DO NI=1,VBMAX
               DO NJ=1,VBMAX
                  DO NB=1,(PROCS*WDES%NBANDS-VBMAX)
                  DO NA=1,(PROCS*WDES%NBANDS-VBMAX)
!                     T2(NA,NB,NI,NJ,KB,KA)=T2(NA,NB,NI,NJ,KB,KA)*(mix)+(VVOO_BA(NA,NB,NI,NJ,KB,KA)*((1.0_qs,0.0_qs)-mix))&
!                         /(W%CELTOT(NI,KI,1)+ &
!                        W%CELTOT(NJ,KJ,1)- &
!                        W%CELTOT(NA+VBMAX,KPTS_MKPTS(KA),1)- &
!                        W%CELTOT(NB+VBMAX,KB,1))
                  ENDDO
                  ENDDO
               ENDDO
               ENDDO
            ENDDO
         ENDDO
         
      END SUBROUTINE UPDATE_T2

      SUBROUTINE APPLY_DENOM_2_IJAB(W,WDES,KI,KJ,KA,KB,XY)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(wavedes) WGW
         TYPE(wavespin) W 
         TYPE (in_struct) IO
         GDEFS :: XY(:,:,:,:)
         INTEGER :: KQ,KI,KJ,KQ_,NI,NJ,KA,KB,NA,NB,NBANDSGW,NC,NC_N
         
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
            DO NA=1,(WDES%NB_TOT-VBMAX)
            DO NB=1,(WDES%NB_TOT-VBMAX)
               XY(NI,NJ,NA,NB)=XY(NI,NJ,NA,NB)/(W%CELTOT(NI,KI,1)+ &
                        W%CELTOT(NJ,KJ,1)- &
                        W%CELTOT(NA+VBMAX,KA,1)- &
                        W%CELTOT(NB+VBMAX,KB,1))
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         
      END SUBROUTINE APPLY_DENOM_2_IJAB

      SUBROUTINE APPLY_DENOM_2_ABIJ(W,WDES,KA,KB,KI,KJ,XY)
         use mkpoints
         use base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(wavedes) WGW
         TYPE(wavespin) W 
         TYPE (in_struct) IO
         GDEFS :: XY(:,:,:,:)
         INTEGER :: KQ,KI,KJ,KQ_,NI,NJ,KA,KB,NA,NB,NBANDSGW,NC,NC_N
         
         DO NI=1,VBMAX
         DO NJ=1,VBMAX
            DO NA=1,(WDES%NB_TOT-VBMAX)
            DO NB=1,(WDES%NB_TOT-VBMAX)
               XY(NA,NB,NI,NJ)=XY(NA,NB,NI,NJ)/(W%CELTOT(NI,KI,1)+ &
                        W%CELTOT(NJ,KJ,1)- &
                        W%CELTOT(NA+VBMAX,KA,1)- &
                        W%CELTOT(NB+VBMAX,KB,1))
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         
      END SUBROUTINE APPLY_DENOM_2_ABIJ

      SUBROUTINE FORM_CHI_GAI(W,WDES,NBANDSGW)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE mkpoints
         USE lattice
         USE constant
         USE base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(wavespin) W
         INTEGER KI,KJ,KA,KK,NI,NJ,NK,NL,KQ,KQ_,KC,KD,KL,NA, &
           ND,NBANDSGW,NE,NKPTS_TMP
         COMPLEX(q) :: KW
         
!        write(*,*)'test1a'
         call BLACS_PINFO(ME,PROCS)
!        write(*,*)'test1b'
         CHI_GAI=zero

        KW=KPOINTS_FULL%WTKPT(1)*(1.0_q,0.0_q)
!        KW=one

        IF (.NOT. LCCD) THEN

           DO KD=1,WDES%NKPTS
              !CALL BCAST2ALL_T2(WDES,KD,T2_TMP)
	
!               DO KQ=1,WDES%NKPTS
                  KQ=GKQ
                  !KPOINT_IN_FULL_GRID(-WDES%VKPT(:,GKQ),KPOINTS_FULL)

!                  KD=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KL)+ &
!                     WDES%VKPT(:,KQ),KPOINTS_FULL)

                  CALL BCAST2ALL_FTOD_AI(WDES,KD,1,PW_AI_TMP,OC_AI_TMP)
                  DO NA=1,(WDES%NB_TOT-VBMAX)
                  DO NI=1,VBMAX
                     PW_IA_TMP(:,NA,NI,:)=PW_AI_TMP(:,NI,NA,:)
                  ENDDO
                  ENDDO
                  
!                 KD=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KL)- &
!                     WDES%VKPT(:,KQ),KPOINTS_FULL)
                  
                 CALL BCAST2ALL_FTOD_AI(WDES,KD,2,PW_AI_TMP,OC_AI_TMP)
		 
               NKPTS_TMP=MY_NKPTS
               IF ((PROCS*MY_NKPTS)<WDES%NKPTS) NKPTS_TMP=NKPTS_TMP+1
               DO KA=1,NKPTS_TMP !MY_NKPTS
                  CALL ORDER_T2(WDES,KD,KPTS_MKPTS(KA),T2_TMP)
		  IF (KA>MY_NKPTS) CYCLE
		  
!                KD=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KL)+ &
!                  WDES%VKPT(:,KQ),KPOINTS_FULL)
!                KQ_=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KD)- &
!                  WDES%VKPT(:,KL),KPOINTS_FULL)
                KI=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KPTS_MKPTS(KA))- &
                  WDES%VKPT(:,KQ),KPOINTS_FULL)

                KQ_=KPOINT_IN_FULL_GRID(-WDES%VKPT(:,KI)+ &
                   WDES%VKPT(:,KPTS_MKPTS(KA)),KPOINTS_FULL)
                VVOO_S(:,:,:,:)=(0.0_qs,0.0_qs)

               VVOO_S(:,:,:,:)=T2_TMP(:,:,:,:) !,KPTS_MKPTS(KA))
               CALL SORT_V1V2O1O2_V1O1V2O2(WDES,VVOO_S,VOVO2_S) !dlai
               VOVO=CONJG(VOVO2_S)

               DO NA=1,(WDES%NB_TOT-VBMAX)
               DO NI=1,VBMAX
                   AIG(NA,NI,:)=(PW_AI_TMP(:,NI,NA,KQ_)) !dlg
               ENDDO
               ENDDO
               
                  CALL ZGEMM('t','n',(NGVECTOR), &
                    (PROCS*WDES%NBANDS-VBMAX)*VBMAX,(PROCS*WDES%NBANDS-VBMAX)*(VBMAX), &
                    KW,AIG(1,1,1),(PROCS*WDES%NBANDS-VBMAX)*(VBMAX), &
                    VOVO(1,1,1,1),(PROCS*WDES%NBANDS-VBMAX)*(VBMAX), &
                    (1.0_q,0.0_q), CHI_GAI(1,1,1,KI,KA,2),NGVECTOR)

!               VOVO=(VOVO2_S)
               !KQ_=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KI)- &
               !    WDES%VKPT(:,KPTS_MKPTS(KA)),KPOINTS_FULL)
               KI=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KPTS_MKPTS(KA))+ &
                    WDES%VKPT(:,KQ),KPOINTS_FULL)

               VVOO_S(:,:,:,:)=T2(:,:,:,:,KD,KA)
               CALL SORT_V1V2O1O2_V2O2V1O1(WDES,VVOO_S,VOVO2_S) !dlai
               VOVO=(VOVO2_S)

               DO NI=1,VBMAX
               DO NA=1,(WDES%NB_TOT-VBMAX)
                   AIG(NA,NI,:)=(PW_IA_TMP(:,NA,NI,KQ)) !dlg
               ENDDO
               ENDDO
                  CALL ZGEMM('t','n',(NGVECTOR), &
                    (PROCS*WDES%NBANDS-VBMAX)*VBMAX,(PROCS*WDES%NBANDS-VBMAX)*(VBMAX), &
                    KW,AIG(1,1,1),(PROCS*WDES%NBANDS-VBMAX)*(VBMAX), &
                    VOVO(1,1,1,1),(PROCS*WDES%NBANDS-VBMAX)*(VBMAX), &
                    (1.0_q,0.0_q), CHI_GAI(1,1,1,KI,KA,1),NGVECTOR)
               ENDDO !KA

!               ENDDO !KQ
           ENDDO !KL
        
        ENDIF !LCCD
      END SUBROUTINE FORM_CHI_GAI

      SUBROUTINE CONTR_CHI_GAI(W,WDES,NBANDSGW,MKB)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE mkpoints
         USE lattice
         USE constant
         USE base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(wavespin) W
         INTEGER KI,KJ,KA,KK,NI,NJ,NK,NL,KQ,KQ_,KC,KD,KL,NA, &
           ND,NBANDSGW,NE,NC,KB,NB,NKPTS_TMP
         INTEGER :: MKB
         COMPLEX(qs) :: KW
         
!        write(*,*)'test1a'
         call BLACS_PINFO(ME,PROCS)
         !PT2_N=(0.0_qs,0.0_qs)
!        write(*,*)'test1b'

        KW=KPOINTS_FULL%WTKPT(1)
        GAI_TMP=zero 
!RETURN
        KB=MKB 
!        write(*,*)'bcasting'
               CALL BCAST2ALL_CHI_GAI(WDES,KB,2,GAI_TMP)
               !CALL BCAST2ALL_FTOD_AI(WDES,KB,2,PW_AI_TMP,OC_AI_TMP)
            
!               DO KI=1,WDES%NKPTS
               KQ=GKQ
               NKPTS_TMP=MY_NKPTS
               IF ((PROCS*MY_NKPTS)<WDES%NKPTS) NKPTS_TMP=NKPTS_TMP+1
               DO KA=1,NKPTS_TMP !MY_NKPTS

               IF (KA>MY_NKPTS) THEN
!                  write(*,*),KA,KQ
                  KI=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KPTS_MKPTS(MY_NKPTS))+ &
                                 WDES%VKPT(:,KQ),KPOINTS_FULL)
!                  write(*,*),KI,KA
               ELSE
                  KI=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KPTS_MKPTS(KA))+ &
                     WDES%VKPT(:,KQ),KPOINTS_FULL)
               ENDIF
!               CALL BCAST2ALL_FTOD_IA(WDES,KI,1,PW_IA_TMP,OC_IA_TMP)
!               write(*,*)ME,'ordering',KI
               CALL ORDER_FTOD_IA(WDES,KI,1,PW_IA_TMP,OC_IA_TMP)

!               write(*,*)ME,'received',KI

               IF (KA>MY_NKPTS) CYCLE
!               DO KA=1,MY_NKPTS

!               KQ=KPOINT_IN_FULL_GRID(-WDES%VKPT(:,KPTS_MKPTS(KA))+ &
!                  WDES%VKPT(:,KI),KPOINTS_FULL)
               KJ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KB)- &
                     WDES%VKPT(:,KQ),KPOINTS_FULL)
               VOVO=zero


!               KQ_=KPOINT_IN_FULL_GRID(-WDES%VKPT(:,KI)+ &
!                  WDES%VKPT(:,KPTS_MKPTS(KA)),KPOINTS_FULL)
               VVOO_S(:,:,:,:)=(0.0_qs,0.0_qs)
               
               !CHI_GAI(1,NA,NI,KI,KA,1)
               !GAI_TMP(1,NB,NJ,KJ)
               
                CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
                  CHI_GAI(1,1,1,KI,KA,1),FTOD_OC_AI(1,1,NA,KQ_,KA,1), &
                  VBMAX*(WDES%NB_TOT-VBMAX), &
                  GAI_TMP(1,1,1,KJ),OC_AI_TMP(1,1,ND,KQ_), &
                  VBMAX*(WDES%NB_TOT-VBMAX),VOVO(1,1,1,1),(0.0_q,0.0_q)) !aibj

                CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
                  PW_IA_TMP(1,1,1,KQ),OC_IA_TMP(1,1,1,KQ), &
                  VBMAX*(WDES%NB_TOT-VBMAX), &
                  GAI_TMP(1,1,1,KJ),OC_AI_TMP(1,1,ND,KQ), &
                  VBMAX*(WDES%NB_TOT-VBMAX),VOVO(1,1,1,1),(2.0_q,0.0_q)) !aibj 

                  !something seems to be wrong with GAI_TMP
                VOVO_S=CONJG(VOVO)
                CALL SORT_V1O1V2O2_V1V2O1O2(WDES,VOVO_S,VVOO_S)
                CALL APPLY_DENOM_2_ABIJ(W,WDES,KPTS_MKPTS(KA),KB,KI,KJ,VVOO_S)
                T2(:,:,:,:,KB,KA)=T2(:,:,:,:,KB,KA)+VVOO_S(:,:,:,:)*(2.0_qs,0.0_qs)*((1.0_qs,0.0_qs)-mix)
               DO NI=1,VBMAX
               DO NJ=1,VBMAX
                  DO NB=1,(WDES%NB_TOT-VBMAX)
!                     T2(:,NB,NI,NJ,KB,KA)=T2(:,NB,NI,NJ,KB,KA)+(VOVO_S(:,NI,NB,NJ))*(2.0_qs,0.0_qs)*((1.0_qs,0.0_qs)-mix)
                      !VVOO_BA(:,NB,NI,NJ,KB,KA)=VVOO_BA(:,NB,NI,NJ,KB,KA)+(VOVO_S(:,NI,NB,NJ))*(2.0_qs,0.0_qs)
                      
                  ENDDO
               ENDDO
               ENDDO

               ENDDO !KA
               !ENDDO !KI

               CALL BCAST2ALL_CHI_GAI(WDES,KB,1,GAI_TMP)
!               KQ=KPOINT_IN_FULL_GRID(-WDES%VKPT(:,GKQ),KPOINTS_FULL)
               KJ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KB)- &
                  WDES%VKPT(:,KQ),KPOINTS_FULL)
               CALL BCAST2ALL_FTOD_IA(WDES,KJ,2,PW_IA_TMP,OC_IA_TMP)

               DO KA=1,MY_NKPTS

!               KQ=KPOINT_IN_FULL_GRID(-WDES%VKPT(:,KB)+ &
!                  WDES%VKPT(:,KJ),KPOINTS_FULL)
               KI=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KPTS_MKPTS(KA))+ &
                     WDES%VKPT(:,KQ),KPOINTS_FULL)

               !VVO(:,:,:,:)=(0.0_qs,0.0_qs)
               
               !CHI_GAI(1,NA,NI,KI,KA,1)
               !GAI_TMP(1,NB,NJ,KJ)
               
!                CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
!                  GAI_TMP(1,1,1,KJ),OC_AI_TMP(1,1,ND,KQ_), &
!                  VBMAX*(WDES%NB_TOT-VBMAX), &
!                  CHI_GAI(1,1,1,KI,KA,2),FTOD_OC_AI(1,1,NA,KQ_,KA,2), &
!                  VBMAX*(WDES%NB_TOT-VBMAX),VOVO(1,1,1,1),zero) !bjai

!		CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
!                  PW_IA_TMP(1,1,1,KQ),OC_IA_TMP(1,1,1,KQ), &
!                  VBMAX*(WDES%NB_TOT-VBMAX), &
!                  CHI_GAI(1,1,1,KI,KA,2),FTOD_OC_AI(1,1,NA,KQ_,KA,2), &
!                  VBMAX*(WDES%NB_TOT-VBMAX),VOVO(1,1,1,1),one) !bjai
		  
!                CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
!                  PW_IA_TMP(1,1,1,KQ),OC_IA_TMP(1,1,1,KQ), &
!                  VBMAX*(WDES%NB_TOT-VBMAX), &
!                  CHI_GAI(1,1,1,KI,KA,2),FTOD_OC_AI(1,1,NA,KQ_,KA,2), &
!                  VBMAX*(WDES%NB_TOT-VBMAX),VOVO(1,1,1,1),zero) !bjai
                 !VOVO=zero
		 
                CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
                  CHI_GAI(1,1,1,KI,KA,1),FTOD_OC_AI(1,1,NA,KQ_,KA,1), &
                  VBMAX*(WDES%NB_TOT-VBMAX), &
                  PW_IA_TMP(1,1,1,KQ),OC_IA_TMP(1,1,1,KQ), &
                  VBMAX*(WDES%NB_TOT-VBMAX),VOVO(1,1,1,1),zero) !bjai
                  
                  VOVO_S=CONJG(VOVO)
                  CALL SORT_V1O1V2O2_V1V2O1O2(WDES,VOVO_S,VVOO_S)
                  CALL APPLY_DENOM_2_ABIJ(W,WDES,KPTS_MKPTS(KA),KB,KI,KJ,VVOO_S)
                  T2(:,:,:,:,KB,KA)=T2(:,:,:,:,KB,KA)+VVOO_S(:,:,:,:)*(2.0_qs,0.0_qs)*((1.0_qs,0.0_qs)-mix)
               DO NI=1,VBMAX
               DO NJ=1,VBMAX
                   DO NB=1,(WDES%NB_TOT-VBMAX)
!                      VVOO_BA(:,NB,NI,NJ,KB,KA)=VVOO_BA(:,NB,NI,NJ,KB,KA)+(VOVO_S(:,NI,NB,NJ))*(2.0_qs,0.0_qs)
                   ENDDO
               ENDDO
               ENDDO

               ENDDO !KA
               !ENDDO !KI

        
      END SUBROUTINE CONTR_CHI_GAI


      Subroutine CALC_T1_T2_N(WDES,WGW,W,ITERATION,NBANDSGW)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE mkpoints
         USE lattice
         USE constant
         USE base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(wavedes) WGW
         TYPE(wavespin) W
         INTEGER ITERATION,NBANDSGW
         INTEGER KI,KJ,KA,KB,NI,NJ,NA,NB,KQ,KQ_,NC,KK
         integer :: time_array1(8),time_array2(8),ems
!         COMPLEX(qs) :: mix

         mix=(0.5_qs,0.0_qs)
         IF (iteration==1) mix=(0.0_qs,0.0_qs)
         !T2_N=(0.0_qs,0.0_qs)
         !PT2_N=(0.0_qs,0.0_qs)
!         VVOO_BA=(0.0_qs,0.0_qs)

         TETSORT=0
         contr_ftod_time=0
         time_fakic=0
         time_con_akic=0

         !write(*,*)'forming chi_gai'
         CALL  FORM_CHI_GAI(W,WDES,NBANDSGW)
         !write(*,*)'contracting'
         DO KB=1,WDES%NKPTS
        ! CALL CONTR_CHI_GAI(W,WDES,NBANDSGW,KB)
         
        !if (ME==0) WRITE(*,*)']'
         VV_S=(0.0_qs,0.0_qs)

         !if (ME==0) WRITE(*,*)'calculating P{T2}'
         !CALL P_T2_N(WDES,NBANDSGW) !ok
         !if (ME==0) WRITE(*,*)'TETSORT=',TETSORT
         !if (ME==0) WRITE(*,*)'TETS(ab)=',TETS
         !if (ME==0) WRITE(*,*)'TET2S(T2)=',TET2S

!         DO KB=1,WDES%NKPTS
            CALL BCAST2ALL_FTOD_AI(WDES,KB,1,PW_AI_TMP,OC_AI_TMP)

            DO KA=1,MY_NKPTS
               KQ=GKQ
!            DO KQ=1,WDES%NKPTS
               KI=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KPTS_MKPTS(KA))+ &
                  WDES%VKPT(:,KQ),KPOINTS_FULL)
               KJ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KB)- &
                  WDES%VKPT(:,KQ),KPOINTS_FULL)
                  
               CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
                PW_AI_TMP(1,1,1,KQ),OC_AI_TMP(1,1,1,KQ), &
                (PROCS*WDES%NBANDS-VBMAX)*VBMAX,&
                FTOD_PW_AI(1,1,1,KQ,KA,2),FTOD_OC_AI(1,1,1,KQ,KA,2),&
                (PROCS*WDES%NBANDS-VBMAX)*VBMAX,OVOV(1,1,1,1),zero)
               CALL SORT_O1V1O2V2_V2V1O2O1(WDES,OVOV, &
                 VVOO)
               !VVOO_BA(:,:,:,:,KB,KA)=VVOO_BA(:,:,:,:,KB,KA)+VVOO(:,:,:,:)
               VVOO_S=VVOO
               CALL APPLY_DENOM_2_ABIJ(W,WDES,KPTS_MKPTS(KA),KB,KI,KJ,VVOO_S)
               T2(:,:,:,:,KB,KA)=T2(:,:,:,:,KB,KA)*mix+((1.0_qs,0.0_qs)-mix)*VVOO_S(:,:,:,:)
               !VVOO_BA(:,:,:,:,KB,KA)=VVOO(:,:,:,:)
!                DO NI=1,VBMAX
!                DO NJ=1,VBMAX
!                   NC=0
!                DO NB=1,(PROCS*WDES%NBANDS-VBMAX)
!                DO NA=1,(PROCS*WDES%NBANDS-VBMAX) !(NBANDSGW-VBMAX)
!                   NC=NC+1
!                   T2_N(NC,NI,NJ,KB,KA)=T2_N(NC,NI,NJ,KB,KA)+VVOO(NA,NB,NI,NJ)
!                ENDDO
!                ENDDO
!                ENDDO
!                ENDDO
!            ENDDO
            ENDDO

         ENDDO
         DO KB=1,WDES%NKPTS
            CALL CONTR_CHI_GAI(W,WDES,NBANDSGW,KB)
         ENDDO
         !if (ME==0) WRITE(*,*)'TETS(ab)=',TETS
         !if (ME==0) WRITE(*,*)'TET2S(T2)=',TET2S
         !if (iteration==9) T2(23:(PROCS*WDES%NBANDS-VBMAX),10:(PROCS*WDES%NBANDS-VBMAX),:,:,:,:,:)=zero
         !if (iteration==9) T2_N(23:(PROCS*WDES%NBANDS-VBMAX),10:(PROCS*WDES%NBANDS-VBMAX),:,:,:,:,:)=zero
 
      END Subroutine CALC_T1_T2_N

      Subroutine CALC_FINAL_T2(WDES,WGW,W,IO,ITERATION,NBANDSGW)
         USE prec
         USE poscar
         USE pseudo
         USE wave_high
         USE full_kpoints
         USE mkpoints
         USE lattice
         USE constant
         USE base
         IMPLICIT NONE
         TYPE(wavedes) WDES
         TYPE(wavedes) WGW
         TYPE(wavespin) W
         TYPE (in_struct) IO
         INTEGER ITERATION,NBANDSGW
         INTEGER KI,KJ,KA,KB,NI,NJ,NA,NB,KQ,KQ_,NC,KK
         integer :: time_array1(8),time_array2(8),ems

         E_CCSD=zero
         E_CCSD_X=zero
         CALL  FORM_CHI_GAI(W,WDES,NBANDSGW)
!         DO KB=1,WDES%NKPTS
         
!         VVOO_BA=(0.0_qs,0.0_qs)
         TETSORT=0
         
         CALL CONTR_CHI_GAI(W,WDES,NBANDSGW,KB)

         CALL BCAST2ALL_FTOD_AI(WDES,KB,2,PW_AI_TMP,OC_AI_TMP)
         DO KB=1,WDES%NKPTS
         DO KA=1,MY_NKPTS
            KQ=KPOINT_IN_FULL_GRID(-WDES%VKPT(:,GKQ),KPOINTS_FULL)
            KI=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KPTS_MKPTS(KA))- &
               WDES%VKPT(:,KQ),KPOINTS_FULL)
            KJ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KB)+ &
               WDES%VKPT(:,KQ),KPOINTS_FULL)
            DO NA=1,WDES%NB_TOT-VBMAX
            DO NB=1,WDES%NB_TOT-VBMAX
               CALL CONTR_FTOD(NGVECTOR,NHVECTOR, &
                FTOD_PW_AI(1,1,NA,KQ,KA,1),FTOD_OC_AI(1,1,NA,KQ,KA,1), &
                (VBMAX),&
                PW_AI_TMP(1,1,NB,KQ),OC_AI_TMP(1,1,NB,KQ), &
                (VBMAX),TE4O(1,1),zero)
               VVOO(NA,NB,:,:)=TE4O(:,:)
            ENDDO
            ENDDO
!            VVOO_BA(:,:,:,:,KB,KA)=VVOO_BA(:,:,:,:,KB,KA)+(VVOO(:,:,:,:))
         ENDDO
         ENDDO

         TET2S=0
         
           CALL UPDATE_T2(WDES,WGW,W,IO,NBANDSGW,ITERATION)
         DO KB=1,WDES%NKPTS
         DO KA=1,MY_NKPTS
          KQ=KPOINT_IN_FULL_GRID(-WDES%VKPT(:,GKQ),KPOINTS_FULL)
          KI=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KPTS_MKPTS(KA))- &
            WDES%VKPT(:,KQ),KPOINTS_FULL)
          KJ=KPOINT_IN_FULL_GRID(WDES%VKPT(:,KB)+ &
             WDES%VKPT(:,KQ),KPOINTS_FULL)
          DO NI=1,VBMAX
          DO NJ=1,VBMAX
!            CALL WRITE_T2_2_VV(WDES,NBANDSGW,VV_S,T2(1,NI,NJ,KB,KA))
            DO NA=1,WDES%NB_TOT-VBMAX
            DO NB=1,WDES%NB_TOT-VBMAX
               IF ((NA>(NBANDSGW-VBMAX)) .and. (NB>(NBANDSGW-VBMAX))) CYCLE
!               VVOO_BA(NA,NB,NI,NJ,KB,KA)=VV_S(NA,NB)
            ENDDO
            ENDDO
            ENDDO
            ENDDO
         ENDDO
         ENDDO

         CALL CALC_CCSD_ENERGY(WDES,WGW,W,IO,NBANDSGW)
!         ENDDO
         
         TETSORT=0


           IF (IO%IU0>=0) WRITE(IO%IU0,*)'Final E_CCSD=',E_CCSD

         
         CALL M_sum_z(WGW%COMM_INTER, E_CCSD, 1)
         CALL M_sum_z(WGW%COMM_INTER, E_CCSD_X, 1)
         
         IF (IO%IU0>=0) WRITE(IO%IU0,*)'Final E_CCSD=',E_CCSD
         IF (IO%IU0>=0) WRITE(IO%IU6,*)'Final E CCSD=',E_CCSD
         IF (IO%IU0>=0) WRITE(IO%IU0,*)'Final E_CCSD_X=',E_CCSD_X
         IF (IO%IU0>=0) WRITE(IO%IU6,*)'Final E CCSD_X=',E_CCSD_X

      END Subroutine CALC_FINAL_T2
      
      SUBROUTINE LOAD_T2(WDES,NBANDSGW)
         USE prec
         USE poscar
         IMPLICIT NONE
         TYPE(wavedes) WDES
         INTEGER :: NBANDSGW

!         ALLOCATE(T2_N((PROCS*WDES%NBANDS-VBMAX),(PROCS*WDES%NBANDS-VBMAX),VBMAX,VBMAX,WDES%NKPTS,WDES%NKPTS,MY_NKPTS))
!        ALLOCATE(T2((NBANDSGW-VBMAX),(NBANDSGW-VBMAX),VBMAX,VBMAX,WDES%NKPTS,WDES%NKPTS,MY_NKPTS))
!        WRITE (DIR_APP  , "('T2',I1,I1)") MOD(WDES%COMM%NODE_ME/10,10),MOD(WDES%COMM%NODE_ME,10)
!        INQUIRE(FILE=DIR_APP, EXIST=T2READ) 
!        IF (T2READ) THEN
!           open(unit = 7,file = DIR_APP, ERR=1983)
!           READ(7,*)T2
         
!           T2_N(1:(NBANDSGW-VBMAX),1:(NBANDSGW-VBMAX),1:VBMAX,1:VBMAX,1:WDES%NKPTS,1:WDES%NKPTS,1:MY_NKPTS)= &
!              T2(1:(NBANDSGW-VBMAX),1:(NBANDSGW-VBMAX),1:VBMAX,1:VBMAX,1:WDES%NKPTS,1:WDES%NKPTS,1:MY_NKPTS)
!           CLOSE(7)
            
!           write(*,*)'T2 read'
!        ENDIF
!1983    CONTINUE
!         DEALLOCATE(T2)
!        ALLOCATE(T2((PROCS*WDES%NBANDS-VBMAX),(PROCS*WDES%NBANDS-VBMAX),VBMAX,VBMAX,WDES%NKPTS,WDES%NKPTS,MY_NKPTS))
!        T2=T2_N
      
      END SUBROUTINE LOAD_T2

#endif
END MODULE rpax
