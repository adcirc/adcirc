C******************************************************************************
C PADCIRC VERSION 45.12 03/17/2006                                            *
C  last changes in this file VERSION 45.08                                    *
C                                                                             *
C This module declares all global variables that are not exclusive to the 3D  *
C routines.  The 3D exclusive variables are declared in global_3DVS           *
C                                                                             *
C******************************************************************************

C 
      MODULE GLOBAL
      USE SIZES
C...
C...SET GLOBAL PARAMETER CONSTANTS
C...

C.....nominal density of water RHOWAT0
      REAL(SZ), PARAMETER ::  RhoWat0=1000.D0

C.....Sigma T value of reference density
      REAL(SZ), PARAMETER ::  SigT0=RHOWAT0-1000.D0 

C.....PI and degrees to radians conversions
      REAL(8), PARAMETER  ::  PI=3.141592653589793D0
      REAL(8), PARAMETER  ::  DEG2RAD = PI/180.D0
      REAL(8), PARAMETER  ::  RAD2DEG = 180.D0/PI

C.....Days to seconds conversion
      REAL(8), PARAMETER :: day2sec=24.d0*3600.d0
      REAL(8), PARAMETER :: sec2day=1.d0/day2sec

C.....parameters used in barrier overflow 
      REAL(SZ), PARAMETER ::  BARMIN=0.04D0

C     jgf46.21 Added support for IBTYPE=52.
      LOGICAL :: IBTYPE52 = .FALSE.
      REAL(SZ),ALLOCATABLE, TARGET ::   ElevDisc(:)
      REAL(SZ) FluxSettlingTime
      INTEGER  FluxSettlingIT

C      real(sz) rampriver,drampriver,rampriver1,rampriver2

C     jgf46.08 Fine grained ramp functions (jgf46.21 split flux b.c.s)
      REAL(SZ) RampExtFlux,DRampExtFlux ! Ramp for external flux b.c.s
      REAL(SZ) RampIntFlux,DRampIntFlux ! Ramp for internal flux b.c.s
      REAL(SZ) RampElev,DRampElev    ! Ramp for elevation boundary conditions.
      REAL(SZ) RampTip,DRampTip      ! Ramp for tidal potential
      REAL(SZ) RampMete,DRampMete    ! Ramp for wind and atmospheric pressure 
      REAL(SZ) RampWRad,DRampWRad    ! Ramp for wave radiation stress
C
C     jgf46.10 user-controlled warning, dump and error levels on elevations
      REAL(SZ) WarnElev       ! elevation at which a warning is issued
      INTEGER iWarnElevDump   ! 1 to dump warning elevs to fort.69
      LOGICAL WarnElevDump    ! .True. if iWarnElevDump.eq.1
      INTEGER WarnElevDumpLimit ! max dumps; ADCIRC terminates if exceeded 
      INTEGER WarnElevDumpCounter ! counter for number of elev dumps to fort.69
      REAL(SZ) ErrorElev      ! ADCIRC terminates if this elev is exceeded
C     
      INTEGER ScreenUnit      ! jgf46.19 I/O unit where screen output is sent

      type HOTSTART_t
        character(10)     :: filename
        integer           :: lun
      end type HOTSTART_t

      type(HOTSTART_t):: hss

      INTEGER,SAVE        :: NP_G, NE_G  ! global number of nodes, elements
C     global number of elevation, velocity, and meteorological stations
      INTEGER, SAVE       :: NSTAE_G, NSTAV_G, NSTAM_G
      INTEGER, SAVE       :: NSTAC_G     ! number of concentration stations
      integer, save       :: COMM        ! MPI communicator.
      integer, save       :: FileFmtVersion      ! File format version.
      integer, save       :: FileFmtMajor = 1
      integer, save       :: FileFmtMinor = 1
      integer, save       :: FileFmtRev   = 0

C...
C...DECLARE ALL ARRAYS
C...
#ifdef CMPI
      INTEGER, ALLOCATABLE ::   IDUMY(:)
      REAL(SZ),ALLOCATABLE ::   DUMY1(:),DUMY2(:)
#endif
      REAL(SZ),ALLOCATABLE, TARGET ::  ETA1(:), ETA2(:), UU2(:), VV2(:)
#ifndef NELMAX
      REAL(SZ),ALLOCATABLE, TARGET ::  ETAMAX(:) ! v46.50 sb 11/11/2006
#endif
#ifndef NUVMAX
      REAL(SZ),ALLOCATABLE, TARGET ::  UMAX(:) ! v46.50 sb 11/11/2006
#endif
      REAL(SZ),ALLOCATABLE, TARGET ::  ET00(:), UU00(:), VV00(:)
      REAL(SZ),ALLOCATABLE, TARGET ::  RMP00(:), RMU00(:), RMV00(:)
      REAL(SZ),ALLOCATABLE, TARGET ::  Pr2(:), WVNXOut(:), WVNYOut(:)
#ifndef NPRMAX
      REAL(SZ),ALLOCATABLE, TARGET ::  PrMin(:) ! v46.50 sb 11/11/2006
#endif
#ifndef NWVMAX
      REAL(SZ),ALLOCATABLE, TARGET ::  WVNOutMax(:) ! v46.50 sb 11/11/2006
#endif
      REAL(SZ),ALLOCATABLE, TARGET ::  CC00(:)
      REAL(SZ),ALLOCATABLE, TARGET ::  EtaDisc(:) ! vjp 9/2006 added for hotstart file
      REAL(SZ),ALLOCATABLE, TARGET ::   CH1(:)
      REAL(SZ),ALLOCATABLE ::   ETAS(:), ETA0(:),UU0(:),VV0(:),ETAS0(:)
      REAL(SZ),ALLOCATABLE ::   UU1(:),VV1(:)
      REAL(SZ),ALLOCATABLE ::   QX1(:),QX2(:),QY1(:),QY2(:)
      REAL(SZ),ALLOCATABLE ::   QX0(:),QY0(:)
      REAL(SZ),ALLOCATABLE ::   TK0(:),TK2(:)
      REAL(SZ),ALLOCATABLE ::   DP(:),SFAC(:)
      REAL(SZ),ALLOCATABLE ::   MOM_LV_X(:),MOM_LV_Y(:),GWCE_LV(:)
      REAL(SZ),ALLOCATABLE ::   CORIF(:)
      REAL(SZ),ALLOCATABLE ::   TPK(:),FFT(:)
      REAL(SZ),ALLOCATABLE ::   FACET(:),ETRF(:)
      REAL(SZ),ALLOCATABLE ::   ESBIN1(:),ESBIN2(:)
      REAL(SZ),ALLOCATABLE ::   TEMP_LV_A(:,:),TEMP_LV_B(:,:)
      REAL(SZ),ALLOCATABLE ::   QN0(:),QN1(:),QN2(:)
      REAL(SZ),ALLOCATABLE ::   EN0(:),EN1(:),EN2(:)
      REAL(SZ),ALLOCATABLE ::   BNDLEN2O3(:)
      REAL(SZ),ALLOCATABLE ::   CSII(:),SIII(:)
      REAL(SZ),ALLOCATABLE ::   ZNGIF1(:),ZNGIF2(:),ZNGIF3(:)
      REAL(SZ),ALLOCATABLE ::   QNAM(:,:),QNPH(:,:)
      REAL(SZ),ALLOCATABLE ::   ENAM(:,:),ENPH(:,:)
      REAL(SZ),ALLOCATABLE ::   QNIN1(:),QNIN2(:)
      REAL(SZ),ALLOCATABLE ::   ENIN1(:),ENIN2(:)
      REAL(SZ),ALLOCATABLE ::   CSI(:),SII(:)
      REAL(SZ),ALLOCATABLE ::   STAIE1(:),STAIE2(:),STAIE3(:)
      REAL(8), ALLOCATABLE ::   XEV(:),YEV(:),SLEV(:),SFEV(:)
      REAL(SZ),ALLOCATABLE ::   STAIV1(:),STAIV2(:),STAIV3(:)
      REAL(8), ALLOCATABLE ::   XEC(:),YEC(:),SLEC(:),SFEC(:)
      REAL(SZ),ALLOCATABLE ::   STAIC1(:),STAIC2(:),STAIC3(:)
      REAL(8), ALLOCATABLE ::   XEM(:),YEM(:),SLEM(:),SFEM(:)
      REAL(SZ),ALLOCATABLE ::   STAIM1(:),STAIM2(:),STAIM3(:)
      REAL(SZ),ALLOCATABLE ::   TRANS_LV_B(:),TRANS_LV_A(:)
      REAL(SZ),ALLOCATABLE ::   SOURSIN(:)
      REAL(SZ),ALLOCATABLE ::   WSX1(:),WSY1(:),PR1(:)
      REAL(SZ),ALLOCATABLE ::   WSX2(:),WSY2(:)
      REAL(SZ),ALLOCATABLE ::   WVNX1(:),WVNY1(:),PRN1(:)
      REAL(SZ),ALLOCATABLE ::   WVNX2(:),WVNY2(:),PRN2(:)
      REAL(SZ),ALLOCATABLE ::   RSNX1(:),RSNY1(:),RSNX2(:),RSNY2(:)
#ifndef NRSMAX
      REAL(SZ),ALLOCATABLE, TARGET :: RSNMax(:) ! v46.50 sb 11/11/2006
#endif
      REAL(SZ),ALLOCATABLE ::   TK(:)
      REAL(8), ALLOCATABLE ::   EMO(:,:),EFA(:,:)
      REAL(8), ALLOCATABLE ::   XEL(:),YEL(:),SLEL(:),SFEL(:)
      REAL(8) ,ALLOCATABLE ::   AREAS(:), TotalArea(:)
      REAL(SZ),ALLOCATABLE ::   XVELAV(:),YVELAV(:),XVELVA(:),YVELVA(:)
      REAL(SZ),ALLOCATABLE ::   ELAV(:),ELVA(:)
      REAL(SZ),ALLOCATABLE ::   AUV11(:),AUV12(:),AUV13(:),AUV14(:)
      REAL(SZ),ALLOCATABLE ::   DUU1(:),DUV1(:),DVV1(:),BSX1(:),BSY1(:)
      REAL(SZ),ALLOCATABLE ::   TIP1(:),TIP2(:)
      REAL(SZ),ALLOCATABLE ::   SALTAMP(:,:),SALTPHA(:,:)
      REAL(SZ),ALLOCATABLE ::   OBCCOEF(:,:),COEF(:,:)
      REAL(SZ),ALLOCATABLE ::   WKSP(:),RPARM(:)
      REAL(SZ),ALLOCATABLE ::   ABD(:,:),ZX(:)
      REAL(SZ),ALLOCATABLE ::   LSXX(:),LSXY(:),LSYX(:),LSYY(:)

      INTEGER,ALLOCATABLE, TARGET ::    NOFF(:), NODECODE(:)
      INTEGER,ALLOCATABLE ::    ME2GW(:)
      INTEGER,ALLOCATABLE ::    NBV(:),LBCODEI(:)
      INTEGER,ALLOCATABLE ::    NEleZNG(:)
      INTEGER,ALLOCATABLE ::    NNODECODE(:)
      INTEGER,ALLOCATABLE ::    NOFFOLD(:) 
      INTEGER,ALLOCATABLE ::    NIBCNT(:)
      INTEGER,ALLOCATABLE ::    NM(:,:)
      INTEGER,ALLOCATABLE ::    NNEIGH(:),MJU(:),NODELE(:)
      INTEGER,ALLOCATABLE ::    NEITAB(:,:),NEITABELE(:,:)
      INTEGER,ALLOCATABLE ::    NIBNODECODE(:)
      INTEGER,ALLOCATABLE ::    LBArray_Pointer(:)
      INTEGER,ALLOCATABLE ::    NNC(:)
      INTEGER,ALLOCATABLE ::    NNE(:)
      INTEGER,ALLOCATABLE ::    NNV(:)
      INTEGER,ALLOCATABLE ::    NNM(:)
      INTEGER,ALLOCATABLE ::    IWKSP(:),IPARM(:),IPV(:)
      INTEGER,ALLOCATABLE ::    NVDLL(:),NBD(:)
      INTEGER,ALLOCATABLE ::    NBDV(:,:)
      INTEGER,ALLOCATABLE ::    NVELL(:)
      INTEGER,ALLOCATABLE ::    NBVV(:,:)
C
CCC START ADDITION BY MCF FOR NETCDFIO
C
      INTEGER,ALLOCATABLE ::    IBTYPENC(:)
C      
CCC FINISH ADDITION BY MCF FOR NETCDFIO
C
C.....for buoyancy forcing
      REAL(SZ),ALLOCATABLE ::   DASigT(:), DASalt(:), DATemp(:)
      REAL(SZ),ALLOCATABLE ::   VIDBCPDXOH(:), VIDBCPDYOH(:)

C.....for scalar transport
      REAL(SZ),ALLOCATABLE ::   DAConc(:)

C.....for internal barrier boundaries with flowthrough pipes
      REAL(SZ),ALLOCATABLE ::   BARLANHTR(:),BARLANCFSPR(:)
      REAL(SZ),ALLOCATABLE ::   BARINHTR(:),BARINCFSBR(:),BARINCFSPR(:)
      REAL(SZ),ALLOCATABLE ::   PIPEHTR(:),PIPECOEFR(:),PIPEDIAMR(:)
      REAL(SZ),ALLOCATABLE ::   BARLANHT(:),BARLANCFSP(:)
      REAL(SZ),ALLOCATABLE ::   FFF(:),FFACE(:)
      REAL(SZ),ALLOCATABLE ::   BTRAN3(:),BTRAN4(:),BTRAN5(:)
      REAL(SZ),ALLOCATABLE ::   BTRAN6(:),BTRAN7(:),BTRAN8(:)
      REAL(SZ),ALLOCATABLE ::   BARINHT(:),BARINCFSB(:),BARINCFSP(:)
      REAL(SZ),ALLOCATABLE ::   PIPEHT(:),PIPECOEF(:),PIPEDIAM(:)
      REAL(SZ),ALLOCATABLE ::   RBARWL1AVG(:),RBARWL2AVG(:)
      REAL(SZ),ALLOCATABLE ::   RPIPEWL1AVG(:),RPIPEWL2AVG(:)
      INTEGER, ALLOCATABLE ::   IBCONN(:),IBCONNR(:),NTRAN1(:),NTRAN2(:)

C.....for wet/dry
      REAL(SZ)                  HABSMIN,HOFF   

cjjw...added 2 lines 
C.....for matrix conditioning,global to allow for reading in from hotstart file
      REAL(SZ)                  EP !jgf45.08 EP global for predictor-corrector
C...
C...DECLARE COMMON BLOCKS
C...
      INTEGER NTSTEPS,ITMV
      REAL(SZ) DT,FMV
      REAL(8) TIMEBEG
      COMMON /MEANSQ/ TIMEBEG,DT,FMV,NTSTEPS,ITMV
C
      INTEGER NHARFR
      COMMON /LSQFREQS/ NHARFR
C     
      INTEGER NP,NOLICA,NOLIFA,NSCREEN,IHOT,ICS,IDen

C...
C...DECLARE REAL(8) AND CHAR VARIABLES, EQUIVALENCES
C...
      REAL(8) STATIM,REFTIM,DTDP
      REAL(8) AVGXY,DIF1R,DIF2R,DIF3R
      REAL(8) AEMIN,AE,AA,A1,A2,A3,X1,X2,X3,X4,Y1,Y2,Y3,Y4
      REAL(8) XL0,XL1,XL2,YL0,YL1,YL2,SLAM0,SFEA0
      REAL(8) WREFTIM,WTIMED,WTIME2,WTIME1,WTIMINC,QTIME1,QTIME2
      REAL(8) FTIMINC,ETIMINC,RSTIME1,RSTIME2,RSTIMINC
      REAL(8) DELX,DELY,DIST,DELDIST,DELETA
      REAL(8),ALLOCATABLE :: AMIG(:),AMIGT(:),FAMIG(:)
      REAL(8),ALLOCATABLE :: PER(:),PERT(:),FPER(:)
      REAL(8),ALLOCATABLE :: FREQ(:),FF(:),FACE(:)
      REAL(8),ALLOCATABLE :: SLAM(:),SFEA(:),X(:),Y(:)

      CHARACTER(80)             :: RUNDES, RUNID, AGRID, AGRID2, AFRIC
      CHARACTER(4)              :: RDES4(20),RID4(20),AID4(20)
      CHARACTER(8)              :: RDES8(10),RID8(10),AID8(10)
      CHARACTER(10)             :: ALPHA
      CHARACTER(5),ALLOCATABLE  :: TIPOTAG(:),BOUNTAG(:),FBOUNTAG(:)

c.ral v45.02 to fix compiler problem, see also change in read_input.F
c      EQUIVALENCE (RDES4(1),RDES8(1),RUNDES), (RID4(1),RID8(1),RUNID),
c     &            (AID4(1),AID8(1),AGRID)

C...
C...DECLARE AND INITIALIZE LOGICAL VARIABLES 
C...  jgf45.06 Default choices are set to .TRUE. in read_input.F
C     jgf45.07 The definitions of many of these are in timestep.F
      LOGICAL :: C2DDI               = .FALSE.
      LOGICAL :: C3D                 = .FALSE.
      LOGICAL :: C3DDSS              = .FALSE.
      LOGICAL :: C3DVS               = .FALSE.
      LOGICAL :: C2D_BTrans          = .FALSE.
      LOGICAL :: C2D_PTrans          = .FALSE.
      LOGICAL :: CBaroclinic         = .FALSE.
      LOGICAL :: CSmag_Eh            = .FALSE.
      LOGICAL :: CGWCE_Orig          = .FALSE.
      LOGICAL :: CGWCE_New           = .FALSE.
      LOGICAL :: CGWCE_Lump          = .FALSE.
      LOGICAL :: CGWCE_LS_KGQ        = .FALSE.
      LOGICAL :: CGWCE_LS_2PartQ     = .FALSE.
      LOGICAL :: CGWCE_LS_2PartV     = .FALSE.
      LOGICAL :: CGWCE_LS_2PartSQ    = .FALSE.
      LOGICAL :: CGWCE_LS_2PartSV    = .FALSE.
      LOGICAL :: CGWCE_Advec_NC      = .FALSE.
      LOGICAL :: CGWCE_Advec_C1      = .FALSE.
      LOGICAL :: CGWCE_Advec_C2      = .FALSE.
      LOGICAL :: CME_Orig            = .FALSE.
      LOGICAL :: CME_New_NC          = .FALSE.
      LOGICAL :: CME_New_C1          = .FALSE.
      LOGICAL :: CME_New_C2          = .FALSE.
      LOGICAL :: CPRECOR             = .FALSE.
      LOGICAL :: CME_LS_IBPQ         = .FALSE.
      LOGICAL :: CME_LS_IBPV         = .FALSE.
      LOGICAL :: CME_LS_IBPSQ        = .FALSE.
      LOGICAL :: CME_LS_IBPSV        = .FALSE.
      LOGICAL :: CME_LS_2PartQ       = .FALSE.
      LOGICAL :: CME_LS_2PartV       = .FALSE.
      LOGICAL :: CME_AreaInt_Orig    = .FALSE.
      LOGICAL :: CME_AreaInt_Corr    = .FALSE.
      LOGICAL :: CTIP                = .FALSE.
      LOGICAL :: CHARMV              = .FALSE.
C...
C...EXPLICITLY DECLARE ADDITIONAL VARIABLES
C...
      INTEGER  NODEDRYMIN,NODEWETMIN,IBTYPE
      INTEGER  IM,JGW,JKI,JME
      INTEGER  JNMM,KMIN,N1,N2,N3,NABOUT
      INTEGER  NBFR,NBOU,NBVI,NBVJ,NCOR,NE,NE2,NP2
      INTEGER  NEIMIN,NEIMAX,NETA,NFFR
      INTEGER  NFLUXB,NFLUXF,NFLUXIB,NFLUXRBC
      INTEGER  NFLUXGBC
      INTEGER  NFLUXIBP,NPIPE
      INTEGER  NFOVER,NHG,NHY,NOLICAT,NOPE,NOUTC 
      INTEGER  NOUTE,NSPOOLE,NOUTV,NSPOOLV,NPRBI
      INTEGER  NRAMP,NRS,NSTAE,NSTARTDRY,NSTAV,NT,NTCYFE
      INTEGER  NTCYFV,NTCYSE,NTCYSV,NTIF,NTIP,NTRSPE,NTRSPV
      INTEGER  NVEL,NVELEXT,NVELME,NWLAT,NWLON,NWS
      INTEGER  IBSTART, ICHA, ICSTP, IDSETFLG, IER
      INTEGER  IESTP
      INTEGER  IFNLCAT, IFNLCT, IFNLFA
      INTEGER  IFWIND, IGCP, IGEP, IGPP, IGVP, IGWP, IHABEG
      INTEGER  IHARIND, IHOTSTP, ILUMP, IMHS, IPSTP
      INTEGER  IREFYR,IREFMO,IREFDAY,IREFHR,IREFMIN, ISLDIA
      INTEGER  ITIME, ITEMPSTP, ITEST, ITHAF, ITHAS
      INTEGER  ITHS, ITITER, ITMAX, IVSTP, IWSTP, IWTIME, IWTIMEP
      INTEGER  IWYR, J12, J13, J21, J23, J31, J32
      INTEGER  JJ, KEMAX, KVMAX, LRC, LUMPT
      INTEGER  MBW, MDF, NBDJ, NBNCTOT, NBW
      INTEGER  NCBND, NCTOT, NDRY, NDSETSC
      INTEGER  NCChange !jgf45.06 flag indicating wetting / drying has occurred
      INTEGER  NDSETSE, NDSETSV, NDSETSW, NEle, NHAGE, NHAGV
      INTEGER  NHAINC, NHASE, NHASV, NHSINC, NHSTAR

      INTEGER  NNBB, NNBB1, NNBB2
      INTEGER  NNBB1R, NNBB1L, NNBB2R, NNBB2L ! sb46.28.sb05.05 11/01/2006
      INTEGER  NNBB1WN, NNBB2WN               ! sb46.28.sb05.05 WN stands for wet neighbors 11/01/2006
      INTEGER  NOUTGC, NOUTGE, NOUTGV, NOUTGW, NOUTM
      INTEGER  NSCOUC, NSCOUE, NSCOUGC, NSCOUGE, NSCOUGV
      INTEGER  NSCOUGW, NSCOUM, NSCOUV, NSPOOLC, NSPOOLGC, NSPOOLGE
      INTEGER  NSPOOLGV, NSPOOLGW, NSPOOLM
      INTEGER  NSTAC, NSTAM, NTCYFC, NTCYFGC, NTCYFGE
      INTEGER  NTCYFGV, NTCYFGW
      INTEGER  NTCYFM, NTCYSC, NTCYSGC, NTCYSGE
      INTEGER  NTCYSGV, NTCYSGW, NTCYSM
      INTEGER  NTRSPC, NTRSPM, NUMITR, NW, NWET, NWSEGWI, NWSGGWI
      INTEGER  NP_GLOBAL

      REAL(SZ) ADVECX, ADVECY, AGIRD, ARG
      REAL(SZ) ARG1, ARG2, ARGJ, ARGJ1, ARGJ2
      REAL(SZ) AUV21, AUV22
      REAL(SZ) BARAVGWT, BEDSTR, BNDLEN2O3NC
      REAL(SZ) CELERITY, COND, CONVCR
      REAL(SZ) DT2, DTO2, DTOHPP
      REAL(SZ) ETIME1, ETIME2, ETRATIO
      REAL(SZ) G, GO3, GO2, GA00, GB00A00, GC00, GDTO2, GFAO2

      REAL(SZ) P11, P22, P33
      REAL(SZ) QFORCEI, QFORCEJ

      REAL(SZ) QUNORM, QUTemp, QVNORM, RAMP, RAMP1, RAMP2, RBARWL
      REAL(SZ) RBARWL1,RBARWL1F, RBARWL2, RBARWL2F, RFF, RFF1, RFF2
      REAL(SZ) SADVDTO3
      REAL(SZ) TADVODT
      REAL(SZ) THAF, THAS, THENALLDSSSTUP, TIMEIT, TOUTFC
      REAL(SZ) TKWET
      REAL(SZ) VEL, VELABS
      REAL(SZ) VELMAX
      REAL(SZ) A00,B00,C00,ANGINN
      REAL(SZ) CORI,COSTHETA,COSTHETA1,COSTSET,CROSS,CROSS1 
      REAL(SZ) DAY,DOTVEC,DRAMP,DUM1,DUM2
      REAL(SZ) H0 
      REAL(SZ) RNDAY,THETA,THETA1
      REAL(SZ) XL
      REAL(SZ) VECNORM,VL1X,VL1Y,VL2X,VL2Y,WLATMAX
      REAL(SZ) WLONMIN,WLATINC,WLONINC, VELMIN

      REAL(8)  RNP_GLOBAL
      REAL(8)  REFSEC   ! required to run in either 32-bit or 64-bit


! CF Start addition for NETCDF I/O 6/21/06
!     -----------------
!     Attribute general
!     -----------------
      CHARACTER(80) :: title, institution, source, history, base_date
      CHARACTER(80) :: references, comments, host, convention, contact
      CHARACTER(7),  ALLOCATABLE  :: STATNUMB(:,:)
      CHARACTER(50), ALLOCATABLE  :: STATNAME(:,:)
      CHARACTER(7),  ALLOCATABLE  :: STATNUMBV(:,:)
      CHARACTER(50), ALLOCATABLE  :: STATNAMEV(:,:)
      INTEGER, PARAMETER :: SNUMLEN = 7
      INTEGER, PARAMETER :: SNAMLEN = 50
      CHARACTER(7),  ALLOCATABLE  :: STATNUMBM(:,:)
      CHARACTER(50), ALLOCATABLE  :: STATNAMEM(:,:)
! CF Finish addition for NETCDF I/O

C     jgf48.03 Arrays for mapping subdomain nodes and elements to full
C     domain ... used in globalio module but allocated and populated
C     in messenger module
      integer,allocatable,save,target :: nodes_lg(:) 
      integer,allocatable,save,target :: imap_el_lg(:)
      integer,allocatable,save,target :: imap_stae_lg(:)
      integer,allocatable,save,target :: imap_stav_lg(:)
      integer,allocatable,save,target :: imap_stam_lg(:)
      integer,allocatable,save,target :: imap_stac_lg(:)
      integer, save            :: float_type
C     jgf48.03 Variable, arrays, and datastructure that are used by both 
C     globalio and writer modules
      integer, parameter :: BUFSIZE_MAX = 131072
      real(sz), save :: buf(BUFSIZE_MAX)
      real(sz), save :: resultBuf(BUFSIZE_MAX)
C   
C     jgf48.03 the following data structure is used to describe
C     the size and shape of the output data. 
      type OutputDataDescript_t
        integer           :: specifier  ! output format from fort.15 (NOUTE,etc)
        REAL(SZ)          :: initial_value 
        ! num_items_per_record is number of "tuples", i.e.
        ! it is equal to 1 if scalar, 2 if 2D vector, 3 if 3D vector, etc
        integer           :: num_items_per_record
        ! array length for the full domain
        integer           :: num_fd_records
        ! array length for this domain (in serial, refers to full domain,
        ! in parallel, refers to the subdomain processed by this PE)
        integer           :: num_records_this
        integer,  pointer :: imap(:)   ! local-to-global mapping
        integer,  pointer :: iarray(:)
        REAL(SZ), pointer :: array(:)  ! scalar data or 1st component of vector
        REAL(SZ), pointer :: array2(:) ! 2nd component of vector
        REAL(SZ), pointer :: array3(:) ! 3rd component of vector
        REAL(SZ), pointer :: array_g(:)  ! full domain data scalar etc
        REAL(SZ), pointer :: array2_g(:) ! full domain 2nd vec component
        REAL(SZ), pointer :: array3_g(:) ! full domain 3rd vec component
      end type OutputDataDescript_t
      INTEGER ITAU0P ! position in tau0 output file for variable tau0
C
C     jgf48.03 Variables shared between messenger and writer modules
      INTEGER,SAVE ::  REALTYPE, DBLETYPE
      ! COMMUNICATORS FOR GLOBAL FILE WRITING
      INTEGER,SAVE, ALLOCATABLE ::  COMM_WRITER(:)
      ! = 0 if this is a compute proc; > 0 if this is a writer node
      INTEGER,SAVE :: WRITER_ID
      ! Following signals are used in messenger.F and writer.F
      INTEGER,SAVE :: SIG_VAL   = 1
      INTEGER,SAVE :: SIG_WRITE = 100
      INTEGER,SAVE :: SIG_TERM  = 999
C     End of variables shared between messenger and writer modules
C

C-------------------end of data declarations----------------------------------C

      CONTAINS
C...
C...Allocate space for Arrays dimensioned by MNE and MNP
C...
      SUBROUTINE ALLOC_MAIN1a()
      IMPLICIT NONE
C
      ALLOCATE ( SLAM(MNP),SFEA(MNP),X(MNP),Y(MNP))
      ALLOCATE ( DP(MNP),SFAC(MNP)) 
      ALLOCATE ( NM(MNE,3)) 
      ALLOCATE ( ETAS(MNP))
      ALLOCATE ( GWCE_LV(MNP))
      ALLOCATE ( UU1(MNP),VV1(MNP))
      ALLOCATE ( QX1(MNP),QY1(MNP))
      ALLOCATE ( NNODECODE(MNP))
      ALLOCATE ( NOFF(MNE),NOFFOLD(MNE)) 
      ALLOCATE ( NNEIGH(MNP),MJU(MNP),NODELE(MNP))
      ALLOCATE ( NIBNODECODE(MNP))
      ALLOCATE ( CH1(MNP),TRANS_LV_B(MNP),TRANS_LV_A(MNP))
      ALLOCATE ( SOURSIN(MNP))
      ALLOCATE ( TK(MNP))
      ALLOCATE ( AREAS(MNE), TotalArea(MNP))
      ALLOCATE ( UU2(MNP),VV2(MNP))
      ALLOCATE ( QX2(MNP),QY2(MNP))
      ALLOCATE ( ETA1(MNP),ETA2(MNP), EtaDisc(MNP))
#ifndef NELMAX
      ALLOCATE ( ETAMAX(MNP)) ! v46.50 sb 11/11/2006
#endif
#ifndef NUVMAX
      ALLOCATE ( UMAX(MNP)) ! v46.50 sb 11/11/2006
#endif
      ALLOCATE ( CORIF(MNP))
      ALLOCATE ( MOM_LV_X(MNP),MOM_LV_Y(MNP))
      ALLOCATE ( LBArray_Pointer(MNP))
      ALLOCATE ( CSI(MNP),SII(MNP))
      ALLOCATE ( NODECODE(MNP))
      ALLOCATE ( NIBCNT(MNP) )
      ALLOCATE ( DASigT(MNP),DATemp(MNP),DASalt(MNP))
      ALLOCATE ( VIDBCPDXOH(MNP), VIDBCPDYOH(MNP))
      ALLOCATE ( DAConc(MNP))
      ALLOCATE ( LSXX(MNP),LSXY(MNP),LSYX(MNP),LSYY(MNP))
ckmd	 Added in parameters for the pc algorithm 
      ALLOCATE ( ETA0(MNP), UU0(MNP), VV0(MNP))
      ALLOCATE ( QX0(MNP), QY0(MNP), ETAS0(MNP))
      ALLOCATE ( TK0(MNP), TK2(MNP))

#ifdef CMPI
      ALLOCATE ( IDUMY(1), DUMY1(1), DUMY2(1) )
#endif

#ifdef CVEC
      ALLOCATE ( TEMP_LV_A(MNE,3),TEMP_LV_B(MNE,3) )
#endif

c     IF (C3DDSS) THEN
c        ALLOCATE(AUV13(MNP),AUV14(MNP))
c     ENDIF

C     jgf These were inside the C3D if-block below, but since they are
C     used in the 2DDI section of timestep.F, they must be allocated in
C     the non-C3D case.

      ALLOCATE(AUV11(MNP),AUV12(MNP))

      IF (C3D) THEN
         ALLOCATE(DUU1(MNP),DUV1(MNP),DVV1(MNP),
     &        BSX1(MNP),BSY1(MNP))
      ENDIF

      EtaDisc(:) = 0.0d0
#ifndef NELMAX
      ETAMAX(:) = -99999.d0 ! v46.50 sb 11/11/2006
#endif
#ifndef NUVMAX
      UMAX(:) = 0.d0 ! v46.50 sb 11/11/2006
#endif

      RETURN
      END SUBROUTINE


C...
C...Allocate space for Arrays needed to determine neighbor tables   
C...
      SUBROUTINE ALLOC_MAIN1b()
      IMPLICIT NONE
C  Neighbor Table

      ALLOCATE ( NEITAB(MNP,MNEI))
      ALLOCATE ( NEITABELE(MNP,MNEI))
C
      RETURN
      END SUBROUTINE


C...
C...Allocate space for Arrays dimensioned by MNOPE and MNETA
C...     
      SUBROUTINE ALLOC_MAIN2()
      IMPLICIT NONE
C    
      ALLOCATE ( ESBIN1(MNETA),ESBIN2(MNETA))
      ALLOCATE ( NBDV(MNOPE,MNETA))
      ALLOCATE ( NVDLL(MNOPE),NBD(MNETA))
      ALLOCATE ( EMO(MNBFR,MNETA),EFA(MNBFR,MNETA))
C      
      RETURN
      END SUBROUTINE


C...
C...Allocate space for nonperiodic zero and nonzero normal flow
C...boundary arrays including barriers
C...
      SUBROUTINE ALLOC_MAIN3()
      IMPLICIT NONE
C
      ALLOCATE ( QN0(MNVEL),QN1(MNVEL),QN2(MNVEL))
      ALLOCATE ( EN0(MNVEL),EN1(MNVEL),EN2(MNVEL))
      ALLOCATE ( ENIN1(MNVEL),ENIN2(MNVEL))
      ALLOCATE ( NBV(MNVEL),LBCODEI(MNVEL))
      ALLOCATE ( BNDLEN2O3(MNVEL))
      ALLOCATE ( ME2GW(MNVEL))
      ALLOCATE ( CSII(MNVEL),SIII(MNVEL))
      ALLOCATE ( ZNGIF1(MNVEL),ZNGIF2(MNVEL),ZNGIF3(MNVEL))
      ALLOCATE ( NEleZNG(MNVEL))
      ALLOCATE ( BARLANHT(MNVEL),BARLANCFSP(MNVEL))
      ALLOCATE ( BARLANHTR(MNVEL),BARLANCFSPR(MNVEL))
      ALLOCATE ( BARINHT(MNVEL),BARINCFSB(MNVEL),BARINCFSP(MNVEL))
      ALLOCATE ( PIPEHT(MNVEL),PIPECOEF(MNVEL),PIPEDIAM(MNVEL))
      ALLOCATE ( IBCONN(MNVEL))
      ALLOCATE ( BARINHTR(MNVEL),BARINCFSBR(MNVEL),BARINCFSPR(MNVEL))
      ALLOCATE ( PIPEHTR(MNVEL),PIPECOEFR(MNVEL),PIPEDIAMR(MNVEL))
      ALLOCATE ( IBCONNR(MNVEL),NTRAN1(MNVEL),NTRAN2(MNVEL))
      ALLOCATE ( BTRAN3(MNVEL),BTRAN4(MNVEL),BTRAN5(MNVEL))
      ALLOCATE ( BTRAN6(MNVEL),BTRAN7(MNVEL),BTRAN8(MNVEL))
      ALLOCATE ( RBARWL1AVG(MNVEL),RBARWL2AVG(MNVEL))
      ALLOCATE ( RPIPEWL1AVG(MNVEL),RPIPEWL2AVG(MNVEL))
      ALLOCATE ( QNIN1(MNVEL),QNIN2(MNVEL))
      ALLOCATE ( NBVV(MNBOU,0:MNVEL))
      ALLOCATE ( NVELL(MNBOU))

CCC START ADDITION BY MCF FOR NETCDFIO

      ALLOCATE ( IBTYPENC(MNBOU))

CCC FINISH ADDITION BY MCF FOR NETCDFIO

C      
C     jgf46.21 Added support for IBTYPE=52.
      ALLOCATE ( ElevDisc(MNVEL))

      RETURN
      END SUBROUTINE


C...
C...Allocate space for tidal potential terms 
C...
      SUBROUTINE ALLOC_MAIN4a()
      IMPLICIT NONE
C
      ALLOCATE ( TPK(MNTIF),AMIGT(MNTIF),FFT(MNTIF) )
      ALLOCATE ( FACET(MNTIF),PERT(MNTIF),ETRF(MNTIF) )
      ALLOCATE ( TIPOTAG(MNTIF) )
      
      IF ( CTIP ) THEN
         ALLOCATE( TIP1(MNP),TIP2(MNP))
      ENDIF  
C      
      RETURN
      END SUBROUTINE


C...     
C...Allocate space for Earth load/self-attraction tide 
C...     
      SUBROUTINE ALLOC_MAIN4b()
      IMPLICIT NONE
C
      ALLOCATE ( SALTAMP(MNTIF,MNP),SALTPHA(MNTIF,MNP) )
C
      RETURN
      END SUBROUTINE


C...     
C...Allocate space for Arrays dimensioned by MNBFR      
C...     
      SUBROUTINE ALLOC_MAIN5()
      IMPLICIT NONE
C
      ALLOCATE ( AMIG(MNBFR),PER(MNBFR))
      ALLOCATE ( FF(MNBFR),FACE(MNBFR))
      ALLOCATE ( BOUNTAG(MNBFR) )
C      
      RETURN
      END SUBROUTINE


C...
C...Allocate space for periodic normal flow boundary conditions
C...
      SUBROUTINE ALLOC_MAIN6()
      IMPLICIT NONE
C
      ALLOCATE ( QNAM(MNFFR,MNVEL),QNPH(MNFFR,MNVEL))
      ALLOCATE ( ENAM(MNFFR,MNVEL),ENPH(MNFFR,MNVEL))
      ALLOCATE ( FBOUNTAG(MNFFR))
      ALLOCATE ( FAMIG(MNFFR), FFF(MNFFR), FFACE(MNFFR), FPER(MNFFR) )
C      
      RETURN
      END SUBROUTINE


C...
C...Allocate space for arrays used for station elevation output
C...
      SUBROUTINE ALLOC_MAIN7()
      IMPLICIT NONE
C
      ALLOCATE ( NNE(MNSTAE),ET00(MNSTAE))
      ALLOCATE ( STAIE1(MNSTAE),STAIE2(MNSTAE),STAIE3(MNSTAE))
      ALLOCATE ( XEL(MNSTAE),YEL(MNSTAE),SLEL(MNSTAE),SFEL(MNSTAE))
C
      RETURN
      END SUBROUTINE


C...     
C...Allocate space for arrays used for station velocity output
C...
      SUBROUTINE ALLOC_MAIN8()
      IMPLICIT NONE
C
      ALLOCATE ( XEV(MNSTAV),YEV(MNSTAV),SLEV(MNSTAV),SFEV(MNSTAV))
      ALLOCATE ( NNV(MNSTAV))
      ALLOCATE ( UU00(MNSTAV),VV00(MNSTAV))
      ALLOCATE ( STAIV1(MNSTAV),STAIV2(MNSTAV),STAIV3(MNSTAV))
C      
      RETURN
      END SUBROUTINE


C...
C...Allocate space for arrays used for station concentration output
C...
      SUBROUTINE ALLOC_MAIN9()
      IMPLICIT NONE
C
      ALLOCATE ( XEC(MNSTAC),YEC(MNSTAC),SLEC(MNSTAC),SFEC(MNSTAC))
      ALLOCATE ( NNC(MNSTAC))
      ALLOCATE ( CC00(MNSTAC))
      ALLOCATE ( STAIC1(MNSTAC),STAIC2(MNSTAC),STAIC3(MNSTAC))
C      
      RETURN
      END SUBROUTINE


C...
C...Allocate space for arrays used for station meteorological output
C...
      SUBROUTINE ALLOC_MAIN10()
      IMPLICIT NONE
C
      ALLOCATE ( XEM(MNSTAM),YEM(MNSTAM),SLEM(MNSTAM),SFEM(MNSTAM))
      ALLOCATE ( NNM(MNSTAM))
      ALLOCATE ( RMU00(MNSTAM),RMV00(MNSTAM),RMP00(MNSTAM))
      ALLOCATE ( STAIM1(MNSTAM),STAIM2(MNSTAM),STAIM3(MNSTAM))
C
      RETURN
      END SUBROUTINE


C...
C...Allocate space for arrays needed by GWCE matrix and iterative solver   
C...
      SUBROUTINE ALLOC_MAIN11()
      IMPLICIT NONE

C  Arrays used by JCG iterative solver
      
      ALLOCATE( OBCCOEF(MNETA,MNEI-1),COEF(MNP,MNEI))
      ALLOCATE( IWKSP(3*MNP),WKSP(4*MNP+400) )
      ALLOCATE( IPARM(12),RPARM(12) )
C
      RETURN
      END SUBROUTINE


C...
C...Allocate space for wind forcing   
C...
      SUBROUTINE ALLOC_MAIN12()
      IMPLICIT NONE
C
      ALLOCATE ( WSX1(MNP),WSY1(MNP),PR1(MNP) )
      ALLOCATE ( WSX2(MNP),WSY2(MNP),PR2(MNP) )
      ALLOCATE ( WVNX1(MNP),WVNY1(MNP),PRN1(MNP) )
      ALLOCATE ( WVNX2(MNP),WVNY2(MNP),PRN2(MNP) )
      ALLOCATE ( RSNX1(MNP),RSNY1(MNP),RSNX2(MNP),RSNY2(MNP) )
      ALLOCATE ( WVNXOUT(MNP),WVNYOUT(MNP) )

#ifndef NPRMAX
      ALLOCATE ( PrMin(MNP) ) ! v46.50 sb 11/11/2006
      PrMin(:) = 99999.d0
#endif
#ifndef NWVMAX
      ALLOCATE ( WVNOUTMax(MNP) ) ! v46.50 sb 11/11/2006
      WVNOUTMAX(:) = 0.d0
#endif
#ifndef NRSMAX
      ALLOCATE ( RSNMax(MNP) )               ! v46.50 sb 11/11/2006
      RSNMax(:) = 0.d0
#endif

C
      RETURN
      END SUBROUTINE

C...
C...Allocate space for harmonic analysis means and variance
C...calculations, this is in global data because the variables are
C...used in main source, outside of HA analysis subroutines.  This
C...should probably be changed.
C...
      SUBROUTINE ALLOC_MAIN14()
      IMPLICIT NONE
C
      ALLOCATE ( XVELAV(MNP),YVELAV(MNP),XVELVA(MNP),YVELVA(MNP) )
      ALLOCATE ( ELAV(MNP),ELVA(MNP) )
C
      RETURN
      END SUBROUTINE

      function VERSION_NUMBER(Major, Minor, Rev) result(vn)

      implicit none
      integer   :: vn, major, minor, rev

      vn = ior(ior(ishft(major,20),ishft(minor,10)),rev)

      return
      end function

      function CMP_VERSION_NUMBERS(a,b) result(match)

      implicit none
      integer  a,b
      logical  :: match

      match = ishft(a,-10) == ishft(b,-10)

      return
      end function
      


      END MODULE GLOBAL
