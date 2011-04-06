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
      USE VERSION
      USE KDTREE2_MODULE   !v49.48.01 tcm fast search algorithm
      IMPLICIT NONE
      SAVE
C...
C...SET GLOBAL PARAMETER CONSTANTS
C...

C.....nominal density of water RHOWAT0
      REAL(SZ), PARAMETER ::  RhoWat0=1000.D0

C.....Sigma T value of reference density
      REAL(SZ), PARAMETER ::  SigT0=RHOWAT0-1000.D0 

C.....Background Atmospheric Pressure in mb
      REAL(SZ), PARAMETER :: PRBCKGRND = 1013.0d0   !tcm 20100617 v49.16 Added

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

C... v49.48 tcm -- added for coupling with STWAVE (NRS=4)
      REAL(8)  :: STARTWAVE   ! Start of Wave Service
      REAL(8)  :: ENDWAVE     ! End of Wave Service

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
      INTEGER :: screenUnit = 6 ! jgf46.19 I/O unit where screen output is sent
C
C     Variables related to hotstarting.
      type HOTSTART_t
        character(10)     :: filename
        integer           :: lun
      end type HOTSTART_t

      type(HOTSTART_t):: hss

      INTEGER             :: NP_G, NE_G  ! global number of nodes, elements
C     global number of elevation, velocity, and meteorological stations
      INTEGER             :: NSTAE_G, NSTAV_G, NSTAM_G
      INTEGER             :: NSTAC_G     ! number of concentration stations

      INTEGER             :: COMM        ! MPI communicator.
      INTEGER             :: FileFmtVersion      ! File format version.
!      INTEGER             :: FileFmtMajor = 1
!      INTEGER             :: FileFmtMinor = 2
!      INTEGER             :: FileFmtRev   = 0
C
C     jgf49.44: These variables represent fulldomain arrays and are 
C     used in hstart.F for reading fulldomain hotstart files as well as in
C     write_output.F/writeHotstart and 
C     writer.F/writeHotstart_through_HSwriter
C     for writing fulldomain hotstart files. 
      REAL(SZ), ALLOCATABLE, TARGET :: ETA1_g(:)
      REAL(SZ), ALLOCATABLE, TARGET :: ETA2_g(:)
      REAL(SZ), ALLOCATABLE, TARGET :: EtaDisc_g(:)
      REAL(SZ), ALLOCATABLE, TARGET :: UU2_g(:)
      REAL(SZ), ALLOCATABLE, TARGET :: VV2_g(:)
      REAL(SZ), ALLOCATABLE, TARGET :: CH1_g(:)
      INTEGER, ALLOCATABLE, TARGET :: NodeCode_g(:)
      INTEGER, ALLOCATABLE :: NNodeCode_g(:)
      INTEGER, ALLOCATABLE, TARGET :: NOFF_g(:)
C...
C...DECLARE ALL ARRAYS
C...

c..RJW merged 09/02/2008 Casey 071219: Added the following variable declaration for 3D mass balance.
C             RESELEM is the number of processors on which each element appears.
      INTEGER, ALLOCATABLE :: RESELEM(:)
c.
#ifdef CMPI
      INTEGER, ALLOCATABLE ::   IDUMY(:)
      REAL(SZ),ALLOCATABLE ::   DUMY1(:),DUMY2(:)
#endif
      REAL(SZ),ALLOCATABLE, TARGET ::  ETA1(:), ETA2(:), UU2(:), VV2(:)
#ifndef NELMAX
      REAL(SZ),ALLOCATABLE, TARGET ::  ETAMAX(:) ! v46.50 sb 11/11/2006
         REAL(SZ),ALLOCATABLE, TARGET ::  HOT_ETAMAX(:) 
         REAL(SZ),ALLOCATABLE, TARGET ::  HOT_ETAMAX_G(:) 
#endif
#ifndef NUVMAX
      REAL(SZ),ALLOCATABLE, TARGET ::  UMAX(:) ! v46.50 sb 11/11/2006
         REAL(SZ),ALLOCATABLE, TARGET ::  HOT_UMAX(:) 
         REAL(SZ),ALLOCATABLE, TARGET ::  HOT_UMAX_G(:) 
#endif
      REAL(SZ),ALLOCATABLE, TARGET ::  ET00(:), UU00(:), VV00(:)
      REAL(SZ),ALLOCATABLE, TARGET ::  RMP00(:), RMU00(:), RMV00(:)
      REAL(SZ),ALLOCATABLE, TARGET ::  RMICE00(:), CICEOUT(:) !v49.64.01 tcm -added for ice concentration
      REAL(SZ),ALLOCATABLE, TARGET ::  Pr2(:), WVNXOut(:), WVNYOut(:)
#ifndef NPRMAX
      REAL(SZ),ALLOCATABLE, TARGET ::  PrMin(:) ! v46.50 sb 11/11/2006
         REAL(SZ),ALLOCATABLE, TARGET ::  HOT_PrMin(:) 
         REAL(SZ),ALLOCATABLE, TARGET ::  HOT_PrMin_G(:) 
#endif
#ifndef NWVMAX
      REAL(SZ),ALLOCATABLE, TARGET ::  WVNOutMax(:) ! v46.50 sb 11/11/2006
         REAL(SZ),ALLOCATABLE, TARGET ::  HOT_WVNOutMax(:) 
         REAL(SZ),ALLOCATABLE, TARGET ::  HOT_WVNOutMax_G(:) 
#endif
#ifdef CSWAN
Casey 090302: Added the following arrays for output of radiation stress gradients.
      REAL(SZ),ALLOCATABLE, TARGET :: RSNXOUT(:), RSNYOUT(:)
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
      REAL(8), ALLOCATABLE, TARGET :: XEV(:) ! x coords of current vel stations
      REAL(8), ALLOCATABLE, TARGET :: YEV(:) ! y coords of current vel stations
      REAL(8), ALLOCATABLE, TARGET :: SLEV(:)! longitudes of current vel sta
      REAL(8), ALLOCATABLE, TARGET :: SFEV(:)! latitudes of current vel sta

      REAL(SZ),ALLOCATABLE ::   STAIV1(:),STAIV2(:),STAIV3(:)
      REAL(8), ALLOCATABLE ::   XEC(:),YEC(:),SLEC(:),SFEC(:)
      REAL(SZ),ALLOCATABLE ::   STAIC1(:),STAIC2(:),STAIC3(:)
      REAL(8), ALLOCATABLE, TARGET :: XEM(:) ! x coords of met stations
      REAL(8), ALLOCATABLE, TARGET :: YEM(:) ! y coords of met stations
      REAL(8), ALLOCATABLE, TARGET :: SLEM(:) ! longitudes of met stations
      REAL(8), ALLOCATABLE, TARGET :: SFEM(:) ! latitudes of met stations

      REAL(SZ),ALLOCATABLE ::   STAIM1(:),STAIM2(:),STAIM3(:)
      REAL(SZ),ALLOCATABLE ::   TRANS_LV_B(:),TRANS_LV_A(:)
      REAL(SZ),ALLOCATABLE ::   SOURSIN(:)
      REAL(SZ),ALLOCATABLE ::   WSX1(:),WSY1(:),PR1(:)
      REAL(SZ),ALLOCATABLE ::   WSX2(:),WSY2(:)
      REAL(SZ),ALLOCATABLE ::   WVNX1(:),WVNY1(:),PRN1(:)
      REAL(SZ),ALLOCATABLE ::   WVNX2(:),WVNY2(:),PRN2(:)
      REAL(SZ),ALLOCATABLE ::   RSNX1(:),RSNY1(:),RSNX2(:),RSNY2(:)
      REAL(SZ),ALLOCATABLE ::   CICE1(:),CICE2(:)   !v49.64.01 tcm -- added 
                                                    !for ice concentration
#ifndef NRSMAX
      REAL(SZ),ALLOCATABLE, TARGET :: RSNMax(:) ! v46.50 sb 11/11/2006
      REAL(SZ),ALLOCATABLE, TARGET :: HOT_RSNMax(:) 
      REAL(SZ),ALLOCATABLE, TARGET :: HOT_RSNMax_g(:) 
#endif
      REAL(SZ),ALLOCATABLE :: TK(:)
      REAL(8),ALLOCATABLE :: EMO(:,:),EFA(:,:)
      REAL(8),ALLOCATABLE,TARGET :: XEL(:)  ! x-coords of elevation stations
      REAL(8),ALLOCATABLE,TARGET :: YEL(:)  ! y-coords of elevation stations
      REAL(8),ALLOCATABLE,TARGET :: SLEL(:) ! latitudes of elevation stations
      REAL(8),ALLOCATABLE,TARGET :: SFEL(:) ! longitudes of elevation stations
      REAL(8),ALLOCATABLE ::   AREAS(:), TotalArea(:)
      REAL(SZ),ALLOCATABLE ::   AUV11(:),AUV12(:),AUV13(:),AUV14(:)
      REAL(SZ),ALLOCATABLE, TARGET ::   DUU1(:),DUV1(:),DVV1(:)
      REAL(SZ),ALLOCATABLE, TARGET ::   BSX1(:),BSY1(:)
      REAL(SZ),ALLOCATABLE ::   TIP1(:),TIP2(:)
      REAL(SZ),ALLOCATABLE ::   SALTAMP(:,:),SALTPHA(:,:)
      REAL(SZ),ALLOCATABLE ::   OBCCOEF(:,:),COEF(:,:)
      REAL(SZ),ALLOCATABLE ::   COEFD(:) !jgf48.4619 from Seizo (Lumped GWCE)
      REAL(SZ),ALLOCATABLE ::   WKSP(:),RPARM(:)
      REAL(SZ),ALLOCATABLE ::   ABD(:,:),ZX(:)
      REAL(SZ),ALLOCATABLE ::   LSXX(:),LSXY(:),LSYX(:),LSYY(:)

C    kmd48.33bc - added variables for the levels of no motion and top temperature boundary condition
      REAL(SZ),ALLOCATABLE ::   LNM_BC1(:), LNM_BC2(:), LNM_BC(:) ! level of no motion
      REAL(SZ),ALLOCATABLE ::   q_heat1(:), q_heat2(:), HFLUX(:)  ! top temperature
      REAL(SZ),ALLOCATABLE, TARGET :: Sponge(:)  ! sponge information for boundary conditions

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
      INTEGER,ALLOCATABLE ::    IBTYPE(:)  !jgf49.32 now ibtype is an array
      INTEGER,ALLOCATABLE ::    IBTYPEE(:) !jgf49.32 added ibtypee
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
      REAL(SZ) DT
C
      INTEGER NP,NOLICA,NOLIFA,IHOT,ICS,IDen
      INTEGER :: NSCREEN = 1

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
C    kmd48.33 - variable added
      REAL(8) DTDPHS  ! time step for hot start

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
C    kmd48.33bc  - logical statement for hot start using different time step
      LOGICAL :: CHOTHS              = .FALSE.
  
C     jgf48.4627 When only meteorological output is requested, skip past
C     the GWCE and Momentum equations, so that ADCIRC runs faster and
C     only calculates the requested quantities.
      LOGICAL :: METONLY             = .FALSE.

C...
C...EXPLICITLY DECLARE ADDITIONAL VARIABLES
C...
      INTEGER  NODEDRYMIN,NODEWETMIN
      INTEGER  IM,JGW,JKI,JME
      INTEGER  JNMM,KMIN,N1,N2,N3
      INTEGER :: NABOUT = 0
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
      INTEGER  IBSTART, ICSTP, IDSETFLG, IER
      INTEGER  IESTP
      INTEGER  IFNLCAT, IFNLCT, IFNLFA
      INTEGER  IFWIND, IGCP, IGEP, IGPP, IGVP, IGWP
      INTEGER  IHOTSTP, ILUMP, IMHS, IPSTP
      INTEGER  IREFYR,IREFMO,IREFDAY,IREFHR,IREFMIN, ISLDIA
      INTEGER  ITIME, ITEMPSTP, ITEST 
      INTEGER  ITHS, ITITER, ITMAX, IVSTP, IWSTP, IWTIME, IWTIMEP
      INTEGER  IWYR, J12, J13, J21, J23, J31, J32
      INTEGER  JJ, KEMAX, KVMAX, LRC, LUMPT
      INTEGER  MBW, MDF, NBDJ, NBNCTOT, NBW
      INTEGER  NCBND, NCTOT, NDRY, NDSETSC
      INTEGER  NCChange !jgf45.06 flag indicating wetting / drying has occurred
      INTEGER  NDSETSE, NDSETSV, NDSETSW, NEle
      INTEGER  NHSINC, NHSTAR
#ifdef CSWAN
Casey 090302: Added a variable for output of radiation stress gradients.
      INTEGER  IGRadS
#endif

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

C    kmd48.33bc - added variables for 3D boundary conditions
      INTEGER  RES_BC_FLAG  ! boundary condition flag for salinity and temperature
      INTEGER  BCFLAG_LNM   ! boundary condition flag for levels of no motion
      INTEGER  BCFLAG_TEMP  ! boundary condition flag for heat flux (top temperature)
      REAL(SZ) RBCTIMEINC, RBCTIME1, RBCTIME2  ! time information for level of no motion
      REAL(SZ) SBCTIMEINC, SBCTIME1, SBCTIME2  ! time information for salinity 
      REAL(SZ) TBCTIMEINC, TBCTIME1, TBCTIME2  ! time information for temperature
      REAL(SZ) TTBCTIMEINC, TTBCTIME1, TTBCTIME2  ! time information for heat flux (top temperature)
      REAL(SZ) BCSTATIM, RBCRATIO
      REAL(SZ) SBCSTATIM, SBCRATIO
      REAL(SZ) TBCSTATIM, TBCRATIO
      REAL(SZ) TTBCSTATIM, TTBCRATIO
      REAL(SZ) SPONGEDIST  ! distance information for the sponge layer
      LOGICAL :: OUTPUTSPONGE  =  .FALSE.

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
      REAL(SZ) THENALLDSSSTUP, TIMEIT, TOUTFC
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
      CHARACTER(80) :: title       ! description of the dataset
      CHARACTER(80) :: institution ! where the original data was produced
      CHARACTER(80) :: source      ! method of production of the original data
      CHARACTER(80) :: history     ! provides an audit trail for modifications to the original data 
      CHARACTER(80) :: base_date 
      CHARACTER(80) :: comments    ! miscellaneous information about the data
      CHARACTER(80) :: host
      CHARACTER(80) :: convention
      CHARACTER(80) :: contact     ! user's institution/ e-mail address
      CHARACTER(80) :: references  ! published/web-based references that describe the data
      CHARACTER(7),  ALLOCATABLE :: STATNUMB(:,:)
      CHARACTER(50), ALLOCATABLE, TARGET  :: STATNAME(:)
      CHARACTER(7),  ALLOCATABLE :: STATNUMBV(:,:)
      CHARACTER(50), ALLOCATABLE, TARGET  :: STATNAMEV(:)
      INTEGER, PARAMETER :: SNUMLEN = 7
      INTEGER, PARAMETER :: SNAMLEN = 50
      CHARACTER(7),  ALLOCATABLE  :: STATNUMBM(:,:)
      CHARACTER(50), ALLOCATABLE, TARGET  :: STATNAMEM(:)
      CHARACTER(50), ALLOCATABLE, TARGET  :: STATNAMEC(:)
      
      INTEGER NCIESTP, NCIPSTP,NCIVSTP,NCIWSTP
      INTEGER NCIGEP,NCIGPP, NCIGVP,NCIGWP
      INTEGER NCCOUE,NCCOUM,NCCOUV,NCCOUW
      INTEGER NCCOUGE,NCCOUGW,NCCOUGV
      
! CF Finish addition for NETCDF I/O

! TCM v49.64.01 Addition for Ice Concentration Fields
      INTEGER :: NCICE,IGIP,IICESTP
      REAL(8) :: CICE_TIMINC,CICE_TIME1,CICE_TIME2

C     jgf48.03 Arrays for mapping subdomain nodes and elements to full
C     domain ... used in globalio module but allocated and populated
C     in messenger module
      integer,allocatable,target :: nodes_lg(:) 
      integer,allocatable,target :: imap_el_lg(:)
      integer,allocatable,target :: imap_stae_lg(:)
      integer,allocatable,target :: imap_stav_lg(:)
      integer,allocatable,target :: imap_stam_lg(:)
      integer,allocatable,target :: imap_stac_lg(:)
      integer                    :: float_type
C     jgf48.03 Variable, arrays, and datastructure that are used by both 
C     globalio and writer modules
      integer, parameter :: BUFSIZE_MAX = 131072
      real(sz) :: buf(BUFSIZE_MAX)
      real(sz) :: resultBuf(BUFSIZE_MAX)
      integer :: integerBuffer(BUFSIZE_MAX)
      integer :: integerResultBuffer(BUFSIZE_MAX)
C   
C     jgf48.03 the following data structure is used to describe
C     the size and shape of the output data. 
C     TCM48.4618 Added ConsiderWetDry, alternate_value and field_name elements to the structure
      type OutputDataDescript_t
        CHARACTER(12)     :: file_name ! name of output file
        character(20)     :: field_name !Name of the field associated with the data structure
        integer           :: specifier  ! output format from fort.15 (NOUTE,etc)
        REAL(SZ)          :: initial_value 
        integer           :: int_initial_value
        logical           :: ConsiderWetDry ! Flag for considering nodes wet/dry status
        real(sz)          :: alternate_value  ! Value that can be used as an alternate value
        REAL(SZ), pointer :: x_coord(:) ! x-coord or longitude (station or mesh)
        REAL(SZ), pointer :: y_coord(:) ! y-coord or lat (station or mesh)
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
        REAL(SZ), pointer :: array2D(:,:) ! subdomain matrix
        INTEGER, pointer  :: iarray_g(:) ! full domain integer
        REAL(SZ), pointer :: array_g(:)  ! full domain data scalar etc
        REAL(SZ), pointer :: array2_g(:) ! full domain 2nd vec component
        REAL(SZ), pointer :: array3_g(:) ! full domain 3rd vec component
        REAL(SZ), pointer :: array2D_g(:,:) ! fulldomain matrix
        REAL(SZ), pointer :: hotstart(:) ! from hotstart (used by min max)
        REAL(SZ), pointer :: hotstart_g(:) ! from hotstart (used by min max)
      end type OutputDataDescript_t
      INTEGER ITAU0P ! position in tau0 output file for variable tau0
C
C     jgf48.03 Variables shared between messenger and writer modules
      INTEGER ::  REALTYPE, DBLETYPE
      ! COMMUNICATORS FOR GLOBAL FILE WRITING
      INTEGER, ALLOCATABLE ::  COMM_WRITER(:)
      ! COMMUNICATORS FOR HOTSTART FILE WRITING !st3
      INTEGER, ALLOCATABLE ::  COMM_WRITEH(:), COMM_HSLEEP(:)

      ! = 0 if this is a compute proc; > 0 if this is a writer node
      INTEGER :: WRITER_ID
      ! Following signals are used in messenger.F and writer.F
      INTEGER :: SIG_VAL   = 1
      INTEGER :: SIG_WRITE = 100
      INTEGER :: SIG_TERM  = 999  ! signal to terminate adc_writer
!  tcm v49.48.01
      integer :: SIG_PAUSE = 99 ! signal to pause adc_writer

C     End of variables shared between messenger and writer modules
C
C     jgf49.44: Log levels, in order from largest amount of log messages
C     written (DEBUG) to fewest log messages written (ERROR). Compared
C     to the value of NABOUT when determining which messages to write
C     to the screen or to log files.
      INTEGER, PARAMETER :: DEBUG = -1  ! write all messages and echo input
      INTEGER, PARAMETER :: ECHO = 0    ! echo input, plus write all non-debug
      INTEGER, PARAMETER :: INFO = 1    ! don't echo input; write all non-debug
      INTEGER, PARAMETER :: WARNING = 2 ! don't echo input; write only warn/err
      INTEGER, PARAMETER :: ERROR = 3   ! don't echo input; only fatal msgs
C
      CHARACTER(len=10), dimension(5) :: logLevelNames  
      CHARACTER(len=50), dimension(50) :: messageSources ! subroutine names
      CHARACTER(len=1024) :: scratchMessage ! used for formatted messages
      INTEGER :: sourceNumber ! index into messageSources for current sub

C     v49.48.01 tcm 20110126 Variables related to kdtree searches
      REAL(SZ), ALLOCATABLE :: BCXY(:,:), RMAX(:)
      INTEGER :: SRCHDP=100
      TYPE(KDTREE2), POINTER :: TREE
      TYPE(KDTREE2_RESULT), ALLOCATABLE :: KDRESULTS(:)

C-------------------end of data declarations----------------------------------C

      CONTAINS

C--------------------------------------------------------------------
C     S U B R O U T I N E    I N I T   L O G G I N G    
C--------------------------------------------------------------------
C     jgf49.44: Initialize the names for the logging levels and the counter
C     for the current subroutine.
C--------------------------------------------------------------------
      SUBROUTINE initLogging()
      IMPLICIT NONE 
C 
      logLevelNames(1) = "DEBUG"
      logLevelNames(2) = "ECHO"
      logLevelNames(3) = "INFO"
      logLevelNames(4) = "WARNING"
      logLevelNames(5) = "ERROR"
      sourceNumber = 0 
C--------------------------------------------------------------------
      END SUBROUTINE initLogging
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     S U B R O U T I N E    S C R E E N   
C--------------------------------------------------------------------
C     jgf49.44: General purpose subroutine to write a message to
C     the screen with a certain "logging level", and subject to the
C     user's selection of where to write screen output. The logging 
C     level is controlled by NABOUT from the fort.15 file. The actual
C     destination of messages written to the screen is controlled by
C     NSCREEN from the fort.15 file. 
C
C     In parallel, only the processor with rank 0 actually writes
C     the message.
C
C     This subroutine assumes that the global variable "caller" has 
C     been set to the name of the subroutine calling it. Therefore,
C     the setMessageSource subroutine must be called at the beginning 
C     of the subroutine that calls this one, and unsetMessageSource
C     must be called at the end. 
C--------------------------------------------------------------------
      SUBROUTINE screen(level, message)
      IMPLICIT NONE 
      INTEGER, intent(in) :: level
      CHARACTER(*), intent(in) :: message
      INTEGER j  ! loop counter for stack
C 
      IF (myProc.eq.0) THEN
         IF (NSCREEN.ne.0) THEN
            IF (level.ge.NABOUT) THEN
#ifdef FULL_STACK
               write(screenUnit,331,advance="no") 
     7             trim(logLevelNames(level+2)),
     &             (trim(messageSources(j)),j=1,sourceNumber) 
               write(screenUnit,332) trim(message)
#else
               write(screenUnit,333) trim(logLevelNames(level+2)),
     &             trim(messageSources(sourceNumber)), trim(message)
#endif
#ifdef FLUSH_MESSAGES
               flush(screenUnit)
#endif
            ENDIF
         ENDIF
      ENDIF
331   FORMAT(A,": ",A,50(:,"->",A))   
332   FORMAT(": ",A)   
333   FORMAT(A,": ",A,": ",A)   
C--------------------------------------------------------------------
      END SUBROUTINE screen
C--------------------------------------------------------------------


C--------------------------------------------------------------------
C     S U B R O U T I N E    L O G   M E S S A G E   
C--------------------------------------------------------------------
C     jgf49.44: General purpose subroutine to write a message to
C     the fort.16 file. In parallel, processors of all ranks will
C     write the message to their own subdomain fort.16 files.
C
C     This subroutine assumes that the global variable "caller" has 
C     been set to the name of the subroutine calling it. Therefore,
C     the setMessageSource subroutine must be called at the beginning 
C     of the subroutine that calls this one, and unsetMessageSource
C     must be called at the end. 
C--------------------------------------------------------------------
      SUBROUTINE logMessage(level, message)
      IMPLICIT NONE 
      INTEGER, intent(in) :: level
      CHARACTER(*), intent(in) :: message
      INTEGER j  ! loop counter for stack
C
      IF (level.ge.NABOUT) THEN
#ifdef FULL_STACK
         write(16,331,advance="no") trim(logLevelNames(level+2)),
     &       (trim(messageSources(j)),j=1,sourceNumber) 
         write(16,332) trim(message)
#else
         write(16,333) trim(logLevelNames(level+2)),
     &       trim(messageSources(sourceNumber)), trim(message)
#endif
#ifdef FLUSH_MESSAGES
               flush(16)
#endif
      ENDIF
331   FORMAT(A,": ",A,50(:,"->",A))   
332   FORMAT(": ",A)   
333   FORMAT(A,": ",A,": ",A)   
C--------------------------------------------------------------------
      END SUBROUTINE logMessage
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     S U B R O U T I N E   A L L    M E S S A G E   
C--------------------------------------------------------------------
C     jgf49.44: General purpose subroutine to write a message to
C     both the screen and to the fort.16 log file.
C--------------------------------------------------------------------
      SUBROUTINE allMessage(level, message)
      IMPLICIT NONE 
      INTEGER, intent(in) :: level
      CHARACTER(*), intent(in) :: message
C 
      call logMessage(level, message)
      call screen(level, message)
C--------------------------------------------------------------------
      END SUBROUTINE allMessage
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     S U B R O U T I N E   S E T   M E S S A G E   S O U R C E   
C--------------------------------------------------------------------
C     jgf49.44: Sets the name of the subroutine that is writing
C     log and/or screen messages. Must use at the start of any subroutine
C     that calls screen, logMessage, or allMessage.
C--------------------------------------------------------------------
      SUBROUTINE setMessageSource(source)
      IMPLICIT NONE 
      CHARACTER(*), intent(in) :: source 
C
      sourceNumber = sourceNumber + 1
      messageSources(sourceNumber) = source
C--------------------------------------------------------------------
      END SUBROUTINE setMessageSource
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     S U B R O U T I N E   U N S E T   M E S S A G E   S O U R C E   
C--------------------------------------------------------------------
C     jgf49.44: Removes the name of the subroutine that is no longer 
C     writing log and/or screen messages. Must use at the end of 
C     any subroutine that calls screen, logMessage, or allMessage.
C--------------------------------------------------------------------
      SUBROUTINE unsetMessageSource()
      IMPLICIT NONE 
C
      sourceNumber = sourceNumber - 1
C--------------------------------------------------------------------
      END SUBROUTINE unsetMessageSource
C--------------------------------------------------------------------


C--------------------------------------------------------------------
C                    S U B R O U T I N E 
C     A L L O C A T E   F U L L   D O M A I N   I O   A R R A Y S 
C--------------------------------------------------------------------
C     jgf49.44: Allocates memory for the full domain arrays for i/o
C     purposes; these arrays are both read from and written to hotstart 
C     files. Arrays that are only written are allocated in the 
C     subroutine that does the writing. 
C--------------------------------------------------------------------
      SUBROUTINE allocateFullDomainIOArrays()
      IMPLICIT NONE 
C 
      call setMessageSource("allocateFullDomainIOArrays")
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ALLOCATE(ETA1_g(NP_G))
      ALLOCATE(ETA2_g(NP_G))
      ALLOCATE(UU2_g(NP_G))
      ALLOCATE(VV2_g(NP_G))
      IF (IM.eq.10) THEN
         ALLOCATE(CH1_g(NP_G))
      ENDIF
      ALLOCATE(EtaDisc_g(NP_G))
      ALLOCATE(NodeCode_g(NP_G))
      ALLOCATE(NOFF_g(NE_G))
      ALLOCATE(NNodeCode_g(NP_G))
C
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C--------------------------------------------------------------------
      END SUBROUTINE allocateFullDomainIOArrays
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     S U B R O U T I N E      B I N A R Y   R E A D   2 D 
C--------------------------------------------------------------------
C     jgf49.44: Reads a real(sz) array with the specified length
C     from a binary hotstart file; takes the file counter as an 
C     input argument and returns the new value of the file counter
C     at the end of the read. Its called "2D" even though it reads
C     a 1D fortran array because the array represents a 2D physical
C     quantity in an unstructured mesh. 
C--------------------------------------------------------------------
      SUBROUTINE binaryRead2D(array, length, lun, counter)
      USE SIZES, ONLY : SZ
      IMPLICIT NONE
      REAL(SZ), intent(out), dimension(:) :: array ! the data to read from the hotstart file
      INTEGER, intent(in) :: length  ! the number of values to read
      INTEGER, intent(in) :: lun     ! fortran logical unit number to read from 
      INTEGER, intent(inout) :: counter ! i/o position in the file 
C 
      INTEGER i      ! array index
C      
      call setMessageSource("binaryRead2D")
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
       DO i=1, length
          READ(lun,REC=counter) array(i)
          counter = counter + 1
       END DO
C
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C--------------------------------------------------------------------
      END SUBROUTINE binaryRead2D
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     S U B R O U T I N E      B I N A R Y   W R I T E   2 D 
C--------------------------------------------------------------------
C     jgf49.49.01: Writes a real(sz) array with the specified length
C     to a binary hotstart file; takes the file counter as an 
C     input argument and returns the new value of the file counter
C     at the end of the write.
C--------------------------------------------------------------------
      SUBROUTINE binaryWrite2D(array, length, lun, counter)
      USE SIZES, ONLY : SZ
      IMPLICIT NONE
      REAL(SZ), intent(in), dimension(:) :: array ! data to write
      INTEGER, intent(in) :: length  ! the number of values to write
      INTEGER, intent(in) :: lun     ! fortran logical unit number to write to 
      INTEGER, intent(inout) :: counter ! i/o position in the file 
C 
      INTEGER i      ! array index
C      
      call setMessageSource("binaryWrite2D")
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
       DO i=1, length
          WRITE(lun,REC=counter) array(i)
          counter = counter + 1
       END DO
C
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C--------------------------------------------------------------------
      END SUBROUTINE binaryWrite2D
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     S U B R O U T I N E      B I N A R Y   I N T   R E A D   2 D 
C--------------------------------------------------------------------
C     jgf49.44: Same as binaryRead2D except it reads an integer.
C--------------------------------------------------------------------
      SUBROUTINE binaryIntRead2D(array, length, lun, counter)
      USE SIZES, ONLY : SZ
      IMPLICIT NONE
      INTEGER, intent(out), dimension(:) :: array ! the data to read from the hotstart file
      INTEGER, intent(in) :: length  ! the number of values to read
      INTEGER, intent(in) :: lun     ! fortran logical unit number to read from 
      INTEGER, intent(inout) :: counter ! i/o position in the file 
C 
      INTEGER i      ! array index
C      
      call setMessageSource("binaryIntRead2D")
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
       DO i=1, length
          READ(lun,REC=counter) array(i)
          counter = counter + 1
       END DO
C
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C--------------------------------------------------------------------
      END SUBROUTINE binaryIntRead2D
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     S U B R O U T I N E      B I N A R Y   R E A D   3 D 
C--------------------------------------------------------------------
C     jgf49.44: Same as binaryRead2D except it reads data with 
C     both horizontal and vertical lengths. Also used for reading
C     harmonic analysis arrays from binary hotstart file.
C--------------------------------------------------------------------
      SUBROUTINE binaryRead3D(array, iLength, jLength, lun,
     &                        counter)
      USE SIZES, ONLY : SZ
      IMPLICIT NONE
      REAL(SZ), intent(out), dimension(:,:) :: array ! data to read from the hotstart file
      INTEGER, intent(in) :: iLength  ! the number of horiz values to read
      INTEGER, intent(in) :: jLength  ! the number of layers
      INTEGER, intent(in) :: lun     ! fortran logical unit number to read from 
      INTEGER, intent(inout) :: counter ! i/o position in the file 
C 
      INTEGER i,j      ! array indices
C      
      call setMessageSource("binaryRead3D")
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
       DO j=1, jLength
          DO i=1, iLength
            READ(lun,REC=counter) array(i,j)
            counter = counter + 1
          END DO
       END DO
C
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C--------------------------------------------------------------------
      END SUBROUTINE binaryRead3D
C--------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E 
C            R E A D  A N D  M A P  T O  S U B D O M A I N  2 D  
C-----------------------------------------------------------------------
C     jgf49.48.04 Reads a fulldomain array in binary format and maps
C     to a subdomain array.
C-----------------------------------------------------------------------
      subroutine readAndMapToSubdomain2D(sd_array, lun, counter)
      IMPLICIT NONE
C
      REAL(SZ), intent(out), dimension(:) :: sd_array ! subdomain array
      INTEGER, intent(in) :: lun     ! fortran logical unit number to read from 
      INTEGER, intent(inout) :: counter ! i/o position in the file 
      INTEGER :: sd_node_number
      INTEGER :: fd_node_number
C 
      REAL(SZ), allocatable, dimension(:) :: fd_array ! full domain array
C
      call setMessageSource("readAndMapToSubdomain2D")
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ! serial
      IF ((MNPROC.eq.1).or.(READ_LOCAL_HOT_START_FILES.eqv..true.)) THEN
         CALL binaryRead2D(sd_array, NP, lun, counter)
      ! parallel
      ELSE
         allocate(fd_array(NP_G))
         CALL binaryRead2D(fd_array, NP_G, lun, counter)
         ! loop over subdomain nodes
         DO sd_node_number=1,np
            ! get the corresponding fulldomain node number 
            fd_node_number = ABS(NODES_LG(sd_node_number))
            ! fill in the subdomain arrays with the corresponding
            ! fulldomain values
            sd_array(sd_node_number) = fd_array(fd_node_number)
         END DO
         deallocate(fd_array)
      ENDIF
C
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE readAndMapToSubdomain2D
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E 
C            R E A D  A N D  M A P  T O  S U B D O M A I N  I N T  2 D  
C-----------------------------------------------------------------------
C     jgf49.48.04 Reads a fulldomain array in binary format and maps
C     to a subdomain array.
C-----------------------------------------------------------------------
      subroutine readAndMapToSubdomainInt2D(sd_array, lun, counter)
      IMPLICIT NONE
C
      INTEGER, intent(out), dimension(:) :: sd_array ! subdomain array
      INTEGER, intent(in) :: lun     ! fortran logical unit number to read from 
      INTEGER, intent(inout) :: counter ! i/o position in the file 
      INTEGER :: sd_node_number
      INTEGER :: fd_node_number
C 
      INTEGER, allocatable, dimension(:) :: fd_array ! full domain array
C
      call setMessageSource("readAndMapToSubdomainInt2D")
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ! serial
      IF ((MNPROC.eq.1).or.(READ_LOCAL_HOT_START_FILES.eqv..true.)) THEN
         CALL binaryIntRead2D(sd_array, NP, lun, counter)
      ! parallel
      ELSE
         allocate(fd_array(NP_G))
         CALL binaryIntRead2D(fd_array, NP_G, lun, counter)
         ! loop over subdomain nodes
         DO sd_node_number=1,np
            ! get the corresponding fulldomain node number 
            fd_node_number = ABS(NODES_LG(sd_node_number))
            ! fill in the subdomain arrays with the corresponding
            ! fulldomain values
            sd_array(sd_node_number) = fd_array(fd_node_number)
         END DO
         deallocate(fd_array)
      ENDIF
C
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE readAndMapToSubdomainInt2D
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     S U B R O U T I N E 
C            R E A D  A N D  M A P  E L E   
C                 T O  S U B D O M A I N  I N T  2 D  
C-----------------------------------------------------------------------
C     jgf49.48.04 Reads a fulldomain array in binary format and maps
C     to a subdomain array (for elemental mappings as opposed to node
C     mapping).
C-----------------------------------------------------------------------
      subroutine readAndMapEleToSubdomainInt2D(sd_array, lun, counter)
      IMPLICIT NONE
C
      INTEGER, intent(out), dimension(:) :: sd_array ! subdomain array
      INTEGER, intent(in) :: lun     ! fortran logical unit number to read from 
      INTEGER, intent(inout) :: counter ! i/o position in the file 
      INTEGER :: sd_ele_number
      INTEGER :: fd_ele_number
C 
      INTEGER, allocatable, dimension(:) :: fd_array ! full domain array
C
      call setMessageSource("readAndMapEleToSubdomainInt2D")
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ! serial
      IF ((MNPROC.eq.1).or.(READ_LOCAL_HOT_START_FILES.eqv..true.)) THEN
         CALL binaryIntRead2D(sd_array, NE, lun, counter)
      ! parallel
      ELSE
         allocate(fd_array(NE_G))
         CALL binaryIntRead2D(fd_array, NE_G, lun, counter)
         ! loop over subdomain nodes
         DO sd_ele_number=1,ne
            ! get the corresponding fulldomain ele number 
            fd_ele_number = ABS(IMAP_EL_LG(sd_ele_number))
            ! fill in the subdomain arrays with the corresponding
            ! fulldomain values
            sd_array(sd_ele_number) = fd_array(fd_ele_number)
         END DO
         deallocate(fd_array)
      ENDIF
C
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE readAndMapEleToSubdomainInt2D
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C     S U B R O U T I N E 
C            R E A D  A N D  M A P  T O  S U B D O M A I N  3 D  
C-----------------------------------------------------------------------
C     jgf49.48.04 Reads a fulldomain array in binary format and maps
C     to a subdomain array.
C-----------------------------------------------------------------------
      subroutine readAndMapToSubdomain3D(sd_array, nfen, lun, counter)
      IMPLICIT NONE
C
      REAL(SZ), intent(out), dimension(:,:) :: sd_array ! subdomain array
      INTEGER, intent(in) :: nfen    ! number of vertical layers  
      INTEGER, intent(in) :: lun     ! fortran logical unit number to read from 
      INTEGER, intent(inout) :: counter ! i/o position in the file 
      INTEGER :: sd_node_number
      INTEGER :: fd_node_number
C 
      REAL(SZ), allocatable, dimension(:,:) :: fd_array ! full domain array
C
      call setMessageSource("readAndMapToSubdomain3D")
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      ! serial
      IF ((MNPROC.eq.1).or.(READ_LOCAL_HOT_START_FILES.eqv..true.)) THEN
         CALL binaryRead3D(sd_array, NP, NFEN, lun, counter)
      ! parallel
      ELSE
         allocate(fd_array(NP_G,NFEN))
         CALL binaryRead3D(fd_array, NP_G, NFEN, lun, counter)
         ! loop over subdomain nodes
         DO sd_node_number=1,np
            ! get the corresponding fulldomain node number 
            fd_node_number = ABS(NODES_LG(sd_node_number))
            ! fill in the subdomain arrays with the corresponding
            ! fulldomain values
            sd_array(sd_node_number,:) = fd_array(fd_node_number,:)
         END DO
         deallocate(fd_array)
      ENDIF
C
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      END SUBROUTINE readAndMapToSubdomain3D
C-----------------------------------------------------------------------

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

C    kmd48.33bc - arrays for baroclinic
      ALLOCATE (q_heat1(MNP), q_heat2(MNP), HFLUX(MNP) )
      ALLOCATE (sponge(MNP) )

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

C... v49.48.02 tcm -- Variables related to kdtree searches
C...    which are deallocated at the end of read_input.F
      ALLOCATE( BCXY(2,MNE), RMAX(MNE) )

      RETURN
      END SUBROUTINE ALLOC_MAIN1a


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
      END SUBROUTINE ALLOC_MAIN1b


C...
C...Allocate space for Arrays dimensioned by MNOPE and MNETA
C...     
      SUBROUTINE ALLOC_MAIN2()
      IMPLICIT NONE
C    
      ALLOCATE ( ESBIN1(MNETA),ESBIN2(MNETA))
      ALLOCATE ( NBDV(MNOPE,MNETA))
      ALLOCATE ( NVDLL(MNOPE),NBD(MNETA), IBTYPEE(MNOPE))
      ALLOCATE ( EMO(MNBFR,MNETA),EFA(MNBFR,MNETA))
C    kmd48.33bc - arrays for boundary conditions
      IF (CBAROCLINIC) THEN
         ALLOCATE ( LNM_BC1(MNETA), LNM_BC2(MNETA) )
         ALLOCATE ( LNM_BC(MNETA) )
      END IF 
C      
      RETURN
      END SUBROUTINE ALLOC_MAIN2


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
      ALLOCATE ( NVELL(MNBOU), IBTYPE(MNBOU))
C      
C     jgf46.21 Added support for IBTYPE=52.
      ALLOCATE ( ElevDisc(MNVEL))
C
C     jgf49.58 Explicitly initialize variables that would otherwise be
C     used without initialization.
      EN2(:) = 0.0D0

      RETURN
      END SUBROUTINE ALLOC_MAIN3


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
      END SUBROUTINE ALLOC_MAIN4a


C...     
C...Allocate space for Earth load/self-attraction tide 
C...     
      SUBROUTINE ALLOC_MAIN4b()
      IMPLICIT NONE
C
      ALLOCATE ( SALTAMP(MNTIF,MNP),SALTPHA(MNTIF,MNP) )
C
      RETURN
      END SUBROUTINE ALLOC_MAIN4b


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
      END SUBROUTINE ALLOC_MAIN5


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
      END SUBROUTINE ALLOC_MAIN6


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
      END SUBROUTINE ALLOC_MAIN7


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
      END SUBROUTINE ALLOC_MAIN8


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
      END SUBROUTINE ALLOC_MAIN9


C...
C...Allocate space for arrays used for station meteorological output
C...
      SUBROUTINE ALLOC_MAIN10()
      IMPLICIT NONE
C
      ALLOCATE ( XEM(MNSTAM),YEM(MNSTAM),SLEM(MNSTAM),SFEM(MNSTAM))
      ALLOCATE ( NNM(MNSTAM))
      ALLOCATE ( RMU00(MNSTAM),RMV00(MNSTAM),RMP00(MNSTAM))
      if (NCICE.NE.0) then
         ALLOCATE ( RMICE00(MNSTAM) )    ! v49.64.01 tcm added for ice concentration
      endif
      ALLOCATE ( STAIM1(MNSTAM),STAIM2(MNSTAM),STAIM3(MNSTAM))
C
      RETURN
      END SUBROUTINE ALLOC_MAIN10


C...
C...Allocate space for arrays needed by GWCE matrix and iterative solver   
C...

      SUBROUTINE ALLOC_MAIN11_LUMPED()
      IMPLICIT NONE
C     jgf48.4619: Array used in lumped GWCE (from Seizo)
      ALLOCATE( COEFD(MNP) ) 
      RETURN
      END SUBROUTINE ALLOC_MAIN11_LUMPED




      SUBROUTINE ALLOC_MAIN11()
      IMPLICIT NONE

C  Arrays used by JCG iterative solver
      
      ALLOCATE( OBCCOEF(MNETA,MNEI-1),COEF(MNP,MNEI))
      ALLOCATE( IWKSP(3*MNP),WKSP(4*MNP+400) )
      ALLOCATE( IPARM(12),RPARM(12) )
C
      RETURN
      END SUBROUTINE ALLOC_MAIN11


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
      if (NCICE.NE.0) then
         ALLOCATE ( CICE1(MNP),CICE2(MNP),CICEOUT(MNP) ) !tcm v49.64.01 added
                                                         ! for ice concentration
      ENDIF
      ALLOCATE ( WVNXOUT(MNP),WVNYOUT(MNP) )
      WSX1(:) = 0.0d0
      WSY1(:) = 0.0d0
      PR1(:) = 0.0d0
      WSX2(:) = 0.0d0
      WSY2(:) = 0.0d0
      PR2(:) = 0.0d0
      IF (NCICE.NE.0) then
         CICE1(:) = 0.0d0 !tcm v49.64.01 added for ice concentration
         CICE2(:) = 0.0d0 !tcm v49.64.01 added for ice concentration
         CICEOUT(:) = 0.0d0 !tcm v49.64.01 added for ice concentration
      ENDIF
#ifdef CSWAN
Casey 090302: Added the next line for output of radiation stress gradients.
      ALLOCATE ( RSNXOUT(MNP),RSNYOUT(MNP) )
#endif

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
      END SUBROUTINE ALLOC_MAIN12

C-------------------------------------------------------------------

      function VERSION_NUMBER(Major, Minor, Rev) result(vn)

      implicit none
      integer   :: vn, major, minor, rev

      vn = ior(ior(ishft(major,20),ishft(minor,10)),rev)

      return
      end function VERSION_NUMBER

C-------------------------------------------------------------------

      function CMP_VERSION_NUMBERS(a,b) result(match)

      implicit none
      integer  a,b
      logical  :: match

      match = ishft(a,-10) == ishft(b,-10)

      return
      end function CMP_VERSION_NUMBERS
C-------------------------------------------------------------------

      !=================================================================
      !=================================================================
      !=================================================================
      !      =====                                            =====
      !      =====            SUBROUTINE PARSE                =====
      !      =====                                            =====
      !=================================================================
      !=================================================================
      !=================================================================

      !=================================================================
      ! This subroutines parse station names so they can passed on to
      ! the NetCDF I/O module of ADCIRC and be written out in the netCDF
      ! output files.
      !
      ! Revision history:
      !
      ! Date      Programmer                       Description of change
      ! ----      ----------                       ---------------------
      ! 12/16/07  Cristina Forbes, PSGS @ UNC-IMS  Wrote parser
      !                                            based on character
      !                                            manipulation subroutines
      !                                            developed at the
      !                                            University of Oslo, Norway
      ! v49.29    Jason Fleming  Fixes to avoid placing leading exclamation
      !                          points in the station names, and to make use
      !                          of trailing exclamation points as station name
      !                          terminators or comment separators.
      !=================================================================
      subroutine parse (line, lvar)
      IMPLICIT NONE
      CHARACTER(132) LINE,L
      CHARACTER(50) LVAR(3)
      INTEGER LS, L1, LEN1, I
      INTEGER NAMSTART  ! position of beginning of station name   
      INTEGER NAMEND    ! position of end of station name
      INTEGER BANGEND   ! location of terminating "!", if any
C
      LS=1
      L=TRIM(ADJUSTL(LINE))
      LEN1=LEN(L)

!     EXTRACT STATION DATA FROM LINE

      DO I=1,2
!        DO I=1,3
         L1=LEN(TRIM(ADJUSTL(L(1:INDEX(L," ")-1))))
         LVAR(I)=ADJUSTL(L(LS:L1))
         L=ADJUSTL(L(L1+1:LEN(L)))
      END DO
C
C     Find start of station name, ignoring leading "!" if it is there
      NAMSTART=1
      IF ( L(1:1).EQ."!" ) THEN
         NAMSTART=2
      ENDIF
C
C     Find end of station name, ignoring the first occurrence of an
C     embedded "!" as well as any characters to the right of that, since
C     they are assumed to be a comment.
      NAMEND=LEN(L)
      BANGEND=INDEX(L(NAMSTART:NAMEND),"!")
      IF ( BANGEND.NE.0 ) THEN
         NAMEND = BANGEND + (NAMSTART-2)
      ENDIF
      LVAR(3)=ADJUSTL(L(NAMSTART:NAMEND)) !jgf49.29: trim leading "!" and spaces
      !=================================================================
      END SUBROUTINE parse
      !=================================================================

      !=================================================================
        SUBROUTINE i2a(in2t, res)
C     jgf49.31 Subroutine to convert an integer number to a string. Moved
C     here from read_input.F to make it more generally available in different
C     modules.
      !=================================================================
        IMPLICIT NONE
!// Arguments
        INTEGER, INTENT(IN)            :: in2t
        CHARACTER(LEN=*), INTENT(OUT)  :: res
!// Local variables
        INTEGER                        :: i,j,k
        CHARACTER(LEN=80)              :: sbuf
        sbuf = ' '
        res = ' '
        i = in2t
        k = 1
        DO
        j = MOD(i,10)
        sbuf(k:k) = ACHAR(j+48)
        i = i / 10
        IF (i <= 0) THEN
        EXIT
        END IF
        k = k + 1
        END DO
        k = LEN_TRIM(sbuf)
        IF (k>1) THEN
        j = k
        i = 1
        DO i = 1, k
         res(i:i) = sbuf(j:j)
         j = j - 1
        END DO
        ELSE
        res(1:1) = sbuf(1:1)
        END IF
      !=================================================================
        END SUBROUTINE i2a
      !=================================================================
C
      !=================================================================
      FUNCTION front_trim(buf) RESULT(res)
C     jgf49.31 Subroutine to return the index of the first non-space 
C     character of a string.
      !=================================================================
      IMPLICIT NONE
      !// Arguments
      CHARACTER(LEN=*), INTENT(IN)   :: buf
      INTEGER                        :: res
      !// Local variables
      INTEGER                        :: i, lng
      lng = LEN_TRIM(buf)
      DO i = 1, lng
         IF (buf(i:i) .NE. ' ') THEN
            res = i
            RETURN
         END IF
      END DO
      !=================================================================
      END FUNCTION front_trim
      !=================================================================
C
      !=================================================================
      FUNCTION a2d(buf) RESULT(res)
C     jgf49.31 Subroutine to convert a double precision real number to
C     a string. Moved here from read_input.F to make it more generally
C     available to different parts of ADCIRC.
      !=================================================================
      IMPLICIT NONE
      !// Arguments
        CHARACTER(LEN=*)               :: buf
        DOUBLE PRECISION               :: res
!// Local variables
        INTEGER                        :: i,j,k,itmp, foffs, lng
        INTEGER                        :: p, q
        DOUBLE PRECISION               :: dtmp
        LOGICAL                        :: neg
        LOGICAL                        :: exponential
        neg = .FALSE.
        exponential = .FALSE.
        lng = LEN_TRIM(buf)
        foffs = front_trim(buf)
        k = INDEX(buf,'.')
        p = INDEX(buf,'e')
        q = INDEX(buf,'E')
        IF(p /= -1 .AND. p /= 0) THEN
        exponential = .TRUE.
        END IF
        IF(q /= -1 .AND. q /= 0) THEN
        exponential = .TRUE.
        p = q
        END IF
        IF (k /= -1 .AND. k /= 0) THEN
!// We have a floating point number
        itmp = 0
!// Get the integer part of the number
        DO i = foffs, k - 1
        IF(buf(i:i) .EQ. '-') THEN
        neg = .TRUE.
        CONTINUE
        END IF
        IF(buf(i:i) .GE. '0' .AND. buf(i:i) .LE. '9' .AND.
     &  buf(i:i) .NE. ' ') THEN
        itmp = itmp * 10
        itmp = itmp + (IACHAR(buf(i:i)) - 48)
        END IF
        END DO
        res = DBLE(itmp)
        dtmp = 0.
        q = 0
        IF(.NOT. exponential) THEN
!// We do not have an exponential number
        DO i = LEN_TRIM(buf), k+1, -1
        itmp = (IACHAR(buf(i:i)) - 48)
        dtmp = dtmp + FLOAT(itmp)
        dtmp = dtmp / 10.
        END DO
        ELSE
!// We have an exponential number
        DO i = p-1, k+1, -1
        itmp = (IACHAR(buf(i:i)) - 48)
        dtmp = dtmp + FLOAT(itmp)
        dtmp = dtmp / 10.
        END DO
        q = a2i(buf(p+1:LEN_TRIM(buf)))
        END IF
        res = res + dtmp
        IF(exponential) THEN
        res = res * 10.**q
        END IF
        IF(neg) THEN
        res = res * (-1)
        END IF
        ELSE
!// We have an integer
        res = DBLE(a2i(buf))
        END IF
      !=================================================================
        END FUNCTION a2d
      !=================================================================
C
      !=================================================================
      FUNCTION a2f(buf) RESULT(res)
C     jgf49.31 Subroutine to convert a single precision real number to
C     a string. Moved here from read_input.F to make it more generally
C     available to different parts of ADCIRC.
      !=================================================================
      IMPLICIT NONE
      !// Arguments
      CHARACTER(LEN=*)               :: buf
      REAL                           :: res
      !// Local variables
      INTEGER                        :: i,j,k,itmp, foffs, lng
      INTEGER                        :: p, q
      REAL                           :: rtmp
      LOGICAL                        :: neg
      LOGICAL                        :: exponential
      neg = .FALSE.
      exponential = .FALSE.
      lng = LEN_TRIM(buf)
      foffs = front_trim(buf)
      k = INDEX(buf,'.')
      p = INDEX(buf,'e')
      q = INDEX(buf,'E')
      IF(p /= -1 .AND. p /= 0) THEN
         exponential = .TRUE.
      END IF
      IF(q /= -1 .AND. q /= 0) THEN
         exponential = .TRUE.
         p = q
      END IF
      IF (k /= -1 .AND. k /= 0) THEN
         !// We have a floating point number
         itmp = 0
         !// Get the integer part of the number
         DO i = foffs, k - 1
            IF(buf(i:i) .EQ. '-') THEN
               neg = .TRUE.
               CONTINUE
            END IF
            IF(buf(i:i) .GE. '0' .AND. buf(i:i) .LE. '9' .AND.
     &           buf(i:i) .NE. ' ') THEN
               itmp = itmp * 10
               itmp = itmp + (IACHAR(buf(i:i)) - 48)
            END IF
         END DO
         res = DBLE(itmp)
         rtmp = 0.
         q = 0
         IF(.NOT. exponential) THEN
!// We do not have an exponential number
            DO i = LEN_TRIM(buf), k+1, -1
               itmp = (IACHAR(buf(i:i)) - 48)
               rtmp = rtmp + FLOAT(itmp)
               rtmp = rtmp / 10.
            END DO
         ELSE
!// We have an exponential number
            DO i = p-1, k+1, -1
               itmp = (IACHAR(buf(i:i)) - 48)
               rtmp = rtmp + FLOAT(itmp)
               rtmp = rtmp / 10.
            END DO
            q = a2i(buf(p+1:LEN_TRIM(buf)))
         END IF
         res = res + rtmp
         IF(exponential) THEN
            res = res * 10.**q
         END IF
         IF(neg) THEN
            res = res * (-1)
         END IF
      ELSE
         !// We have an integer
         res = REAL(a2i(buf))
      END IF
      !=================================================================
      END FUNCTION a2f
      !=================================================================
C
      !=================================================================
      FUNCTION a2i(buf) RESULT(res)
C     jgf49.31 Subroutine to convert a string to an integer.
C     Moved here from read_input.F to make it more generally
C     available to different parts of ADCIRC.
      !=================================================================
      IMPLICIT NONE
      !// Arguments
      CHARACTER(LEN=*), INTENT(IN)   :: buf
      INTEGER                        :: res
      !// Local variables
      INTEGER                        :: i, foffs, lng
      LOGICAL                        :: neg
      neg = .FALSE.
      lng = LEN_TRIM(buf)
      foffs = front_trim(buf)
      res = 0
      DO i = foffs, lng
         IF(buf(i:i) .EQ. '-') THEN
            neg = .TRUE.
           CONTINUE
         END IF
         IF(buf(i:i) .GE. '0' .AND. buf(i:i) .LE. '9' .AND.
     &        buf(i:i) .NE. ' ') THEN
            res = res * 10
            res = res + (IACHAR(buf(i:i)) - 48)
         END IF
      END DO
      IF(neg) THEN
         res = res * (-1)
      END IF
      !=================================================================
      END FUNCTION a2i
      !=================================================================
C     ------------------------------------------------------------------
C     ------------------------------------------------------------------
      END MODULE GLOBAL
C     ------------------------------------------------------------------
C     ------------------------------------------------------------------
