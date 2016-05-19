C******************************************************************************
C PADCIRC VERSION 45.12 03/17/2006                                            *
C  last changes in this file VERSION 45.08                                    *
C... TCM V50.66.01 --ADDING TIME DEPENDENT BATHYMETRY                         *
C                 ADDED NDDT,BTIMINC,DP1,DP2,BTIME1,BTIME2                    *
C                                                                             *
C This module declares all global variables that are not exclusive to the 3D  *
C routines.  The 3D exclusive variables are declared in global_3DVS           *
C                                                                             *
C******************************************************************************

C
      MODULE GLOBAL
      USE VERSION
      USE SIZES, ONLY : SZ, MNP, MNE, MNBFR, MNETA, MNSTAE, MNSTAV, 
     &   MNEI, MNTIF, MNVEL, MYPROC, MNFFR, MNSTAC, MNSTAM
      USE KDTREE2_MODULE   !v49.48.01 tcm fast search algorithm
      IMPLICIT NONE
      SAVE
C...
C...SET GLOBAL PHYSICAL CONSTANTS
C...
      REAL(sz) :: g = 9.80665 !default gravitational acceleration
C.....nominal density of water RHOWAT0
      REAL(SZ), PARAMETER ::  RhoWat0=1000.D0
C.....Sigma T value of reference density
      REAL(SZ), PARAMETER ::  SigT0=RHOWAT0-1000.D0
C.....Background Atmospheric Pressure in mb
      REAL(SZ), PARAMETER :: PRBCKGRND = 1013.0d0   !tcm 20100617 v49.16 Added
      REAL(8) :: rhoAir = 1.15d0         ! kg/m^3
      REAL(8) :: windReduction = 0.9d0
      REAL(8) :: one2ten = 0.8928d0
      REAL(8) :: ten2one = 1.d0/0.8928d0
Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
      REAL(SZ) :: WaveWindMultiplier = 1.D0
C.....PI etc
      REAL(8), PARAMETER ::  PI=3.141592653589793D0
      REAL(8), PARAMETER ::  DEG2RAD = PI/180.D0  ! degrees to radians
      REAL(8), PARAMETER ::  RAD2DEG = 180.D0/PI  ! radians to degrees
      REAL(8), PARAMETER ::  e=2.718281828459045d0
C.....Days to seconds conversion
      REAL(8), PARAMETER :: day2sec=24.d0*3600.d0
      REAL(8), PARAMETER :: sec2day=1.d0/day2sec
      ! jgf50.38.03: Additional useful constants
      REAL(8), PARAMETER :: Rearth = 6378206.4d0 ! radius of earth (m)
      REAL(8), PARAMETER :: nm2m = 1852.d0 ! nautical miles to meters
      REAL(8), PARAMETER :: m2nm = 1.d0/nm2m ! meters to nautical miles
      REAL(8), PARAMETER :: kt2ms = nm2m / 3600.0 ! knots to m/s
      REAL(8), PARAMETER :: ms2kt = 1.d0 / kt2ms ! m/s to knots
      REAL(8), PARAMETER :: omega = 2.0d0*pi / 86164.2d0
      REAL(8), PARAMETER :: mb2pa = 100.d0
C     jgf46.21 Added support for IBTYPE=52.
      REAL(SZ),ALLOCATABLE, TARGET ::   ElevDisc(:)
      REAL(SZ) FluxSettlingTime
      INTEGER  FluxSettlingIT
C... v49.48 tcm -- added for coupling with STWAVE (NRS=4)
      REAL(8)  :: STARTWAVE   ! Start of Wave Service
      REAL(8)  :: ENDWAVE     ! End of Wave Service
      LOGICAL  :: CPL2STWAVE = .FALSE.  !Coupled to STWAVE (NRS=4)
      LOGICAL  :: Flag_ElevError = .FALSE.  
C      real(sz) rampriver,drampriver,rampriver1,rampriver2
C     jgf46.08 Fine grained ramp functions (jgf46.21 split flux b.c.s)
      REAL(SZ) RampExtFlux,DRampExtFlux ! Ramp for external flux b.c.s
      REAL(SZ) RampIntFlux,DRampIntFlux ! Ramp for internal flux b.c.s
      REAL(SZ) RampElev,DRampElev    ! Ramp for elevation boundary conditions.
      REAL(SZ) RampTip,DRampTip      ! Ramp for tidal potential
      REAL(SZ) RampMete,DRampMete    ! Ramp for wind and atmospheric pressure
      REAL(SZ) RampWRad,DRampWRad    ! Ramp for wave radiation stress
Corbitt 1203022: Added Zachs Fix for Assigning a Start Time to Mete Ramping
      REAL(SZ) DUnRampMete           ! Time to not have ramped Mete
C
C     jgf46.10 user-controlled warning, dump and error levels on elevations
      REAL(SZ) WarnElev       ! elevation at which a warning is issued
      INTEGER iWarnElevDump   ! 1 to dump warning elevs to fort.69
      LOGICAL WarnElevDump    ! .True. if iWarnElevDump.eq.1
      INTEGER WarnElevDumpLimit ! max dumps; ADCIRC terminates if exceeded
      INTEGER WarnElevDumpCounter ! counter for number of elev dumps to fort.69
      REAL(SZ) ErrorElev      ! ADCIRC terminates if this elev is exceeded
      INTEGER :: screenUnit = 6 ! jgf46.19 I/O unit where screen output is sent
      INTEGER :: Terminate_LocalProc = 0 ! zc - allow a single processor to initiate code shutdown
C
C     Variables related to hotstarting.
      type HOTSTART_t
        character(20)     :: filename
        integer           :: lun
      end type HOTSTART_t

      type(HOTSTART_t):: hss

      !-----------------------------------------------------------------
      !        F U L L   D O M A I N   V A R I A B L E S   
      !                  O N L Y   U S E D   
      !        I N   P A R A L L E L   E X E C U T I O N
      !-----------------------------------------------------------------
      ! jgf51.21.25: The variables in this section are only used in
      ! parallel execution but they are declared here (rather than in
      ! the messenger module) because they are used in many different 
      ! modules, including harm.F, write_output.F, globalio.F et al.
      !-----------------------------------------------------------------
      integer :: comm        ! MPI communicator.
      integer :: np_g, ne_g  ! global (full domain) number of nodes, elements

C     jgf48.03 Arrays for mapping subdomain nodes and elements to full
C     domain ... used in globalio module but allocated and populated
C     in messenger module
      integer,allocatable,target :: nodes_lg(:)
      integer,allocatable,target :: imap_el_lg(:) 
C
C     global number of elevation, velocity, and meteorological stations
      integer :: nstae_g, nstav_g, nstam_g
      integer :: nstac_g     ! number of concentration stations
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
      ! full domain min and max files in parallel
      !tcm v51.20.01 added variable for time of max (time is referenced off TimeLoc)
      REAL(SZ),ALLOCATABLE, TARGET :: HOT_ETAMAX_G(:)
      REAL(SZ),ALLOCATABLE, TARGET :: ETAMAX_Time(:) 
      REAL(SZ),ALLOCATABLE, TARGET :: ETAMAX_Time_G(:) 
      REAL(SZ),ALLOCATABLE, TARGET :: HOT_UMAX_G(:)
      REAL(SZ),ALLOCATABLE, TARGET :: UMAX_Time(:) 
      REAL(SZ),ALLOCATABLE, TARGET :: UMAX_Time_G(:) 
      REAL(SZ),ALLOCATABLE, TARGET :: HOT_PrMin_G(:)
      REAL(SZ),ALLOCATABLE, TARGET :: PrMin_Time(:) 
      REAL(SZ),ALLOCATABLE, TARGET :: PrMin_Time_G(:) 
      REAL(SZ),ALLOCATABLE, TARGET :: HOT_WVNOutMax_G(:)
      REAL(SZ),ALLOCATABLE, TARGET :: WVNOutMax_Time(:) 
      REAL(SZ),ALLOCATABLE, TARGET :: WVNOutMax_Time_G(:) 
      REAL(SZ),ALLOCATABLE, TARGET :: HOT_RSNMax_g(:)
      REAL(SZ),ALLOCATABLE, TARGET :: RSNMax_Time(:) 
      REAL(SZ),ALLOCATABLE, TARGET :: RSNMax_Time_G(:) 
     
#if defined CSWAN || defined ADCSWAN
      !...Global SWAN Arrays
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_HSOut_g(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_DIROut_g(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_TM01Out_g(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_TPSOut_g(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_WINDXOut_g(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_WINDYOut_g(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_TM02Out_g(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_TMM10Out_g(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_HSMaxOut_g(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_DIRMaxOut_g(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_TM01MaxOut_g(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_TPSMaxOut_g(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_WINDMaxOut_g(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_TM02MaxOut_g(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_TMM10MaxOut_g(:)
#endif
      
      ! mappings from subdomain to fulldomain for stations
      integer,allocatable,target :: imap_stae_lg(:)
      integer,allocatable,target :: imap_stav_lg(:)
      integer,allocatable,target :: imap_stam_lg(:)
      integer,allocatable,target :: imap_stac_lg(:)
C
C     jgf48.03 Variable, arrays, and datastructure that are used by both
C     globalio and writer modules
      integer, parameter :: BUFSIZE_MAX = 131072
      real(sz) :: buf(BUFSIZE_MAX)
      real(sz) :: resultBuf(BUFSIZE_MAX)
      integer :: integerBuffer(BUFSIZE_MAX)
      integer :: integerResultBuffer(BUFSIZE_MAX)
      integer                    :: float_type
C
C     jgf48.03 Variables shared between messenger and writer modules
      integer ::  realtype, dbletype
      !
      ! there is a separate mpi communicator for each dedicated
      ! writer processor; their communicator ids are stored in the 
      ! following array
      integer, allocatable ::  comm_writer(:)
      ! communicators for hotstart file writing !st3
      integer, allocatable ::  comm_writeh(:), comm_hsleep(:)
      !
      ! = 0 if this is a compute proc; > 0 if this is a writer node
      integer :: writer_id
      ! following signals are used in messenger.f and writer.f
      integer, parameter :: sig_val   = 1
      integer, parameter :: sig_write = 100
      integer, parameter :: sig_mesh = 14   ! signal to read mesh for xdmf output
      integer, parameter :: sig_term  = 999  ! signal to terminate adc_writer
!  tcm v49.48.01
      integer, parameter :: sig_pause = 99 ! signal to pause adc_writer
      !----------------------------------------------------------------
C     End of variables shared between messenger and other modules
      !----------------------------------------------------------------
C     
      INTEGER             :: FileFmtVersion      ! File format version.

c..RJW merged 09/02/2008 Casey 071219: Added the following variable declaration for 3D mass balance.
C             RESELEM is the number of processors on which each element appears.
      INTEGER, ALLOCATABLE :: RESELEM(:)
c.
#ifdef CMPI
      INTEGER, ALLOCATABLE ::   IDUMY(:)
      REAL(SZ),ALLOCATABLE ::   DUMY1(:),DUMY2(:)
#endif
      !jgf50.44: Added variable to record the time that the min or max
      ! has occurred
      REAL(SZ),ALLOCATABLE, TARGET ::  ETA1(:), ETA2(:), UU2(:), VV2(:)
      REAL(SZ),ALLOCATABLE, TARGET ::  ETAMAX(:) ! v46.50 sb 11/11/2006
      REAL(SZ),ALLOCATABLE, TARGET ::  UMAX(:) ! v46.50 sb 11/11/2006
      REAL(SZ),ALLOCATABLE, TARGET ::  ET00(:), UU00(:), VV00(:)
      REAL(SZ),ALLOCATABLE, TARGET ::  RMP00(:), RMU00(:), RMV00(:)
      REAL(SZ),ALLOCATABLE, TARGET ::  RMICE00(:), CICEOUT(:) !v49.64.01 tcm -added for ice concentration
      REAL(SZ),ALLOCATABLE, TARGET ::  Pr2(:), WVNXOut(:), WVNYOut(:)
      REAL(SZ),ALLOCATABLE, TARGET ::  PrMin(:) ! v46.50 sb 11/11/2006
      REAL(SZ),ALLOCATABLE, TARGET ::  WVNOutMax(:) ! v46.50 sb 11/11/2006
      
#if defined CSWAN || defined ADCSWAN     
      !...Local SWAN Arrays
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_HSOut(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_DIROut(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_TM01Out(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_TPSOut(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_WINDXOut(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_WINDYOut(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_TM02Out(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_TMM10Out(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_HSMaxOut(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_DIRMaxOut(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_TM01MaxOut(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_TPSMaxOut(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_WINDMaxOut(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_TM02MaxOut(:)
      REAL(SZ), ALLOCATABLE, TARGET      :: SWAN_TMM10MaxOut(:)
#endif      

!  tcm v50.75 removed ifdef cswan to allow for use whenever nrs=3 or nrs=4
!#ifdef CSWAN
Casey 090302: Added the following arrays for output of radiation stress gradients.
      REAL(SZ),ALLOCATABLE, TARGET :: RSNXOUT(:), RSNYOUT(:)
!#endif

#if defined CSWAN || defined ADCSWAN      
Cobell 20120510: Added logicals for turning on/off swan output files   
      LOGICAL                      :: SWAN_OutputHS = .TRUE.
      LOGICAL                      :: SWAN_OutputDIR = .TRUE.
      LOGICAL                      :: SWAN_OutputTM01 = .FALSE.
      LOGICAL                      :: SWAN_OutputTPS = .TRUE.
      LOGICAL                      :: SWAN_OutputWIND = .FALSE.
      LOGICAL                      :: SWAN_OutputTM02 = .FALSE.
      LOGICAL                      :: SWAN_OutputTMM10 = .TRUE.
      LOGICAL                      :: SWAN_OutputAgg(7)
#endif

Cobell - Added for Time Varying Weirs, Set to default values
      LOGICAL                      :: USE_TVW = .FALSE.     !..From namelist, used time varying weirs
      CHARACTER(80)                :: TVW_FILE = 'fort.tvw' !..From namelist, file for time varying weirs
      REAL(SZ),ALLOCATABLE,TARGET  :: TVW(:)                !..Current elevation of weirs, for output


Cobell -  Added the following for output of ESL at end of run
      REAL(SZ),ALLOCATABLE, TARGET :: ESLONOFF(:)

      REAL(SZ),ALLOCATABLE, TARGET ::  CC00(:)
      REAL(SZ),ALLOCATABLE, TARGET ::  EtaDisc(:) ! vjp 9/2006 added for hotstart file
      REAL(SZ),ALLOCATABLE, TARGET ::   CH1(:)
      REAL(SZ),ALLOCATABLE ::   ETAS(:), ETA0(:),UU0(:),VV0(:),ETAS0(:)
      REAL(SZ),ALLOCATABLE ::   UU1(:),VV1(:)
      REAL(SZ),ALLOCATABLE ::   QX1(:),QX2(:),QY1(:),QY2(:)
      REAL(SZ),ALLOCATABLE ::   QX0(:),QY0(:)
      REAL(SZ),ALLOCATABLE ::   TK0(:),TK2(:)
      REAL(SZ),ALLOCATABLE ::  DP1(:),DP2(:) !TCM V50.66.01 FOR TIME VARYING BATHYMETRY
      REAL(SZ),ALLOCATABLE ::   MOM_LV_X(:),MOM_LV_Y(:),GWCE_LV(:)
      REAL(SZ),ALLOCATABLE ::   CORIF(:)
      REAL(SZ),ALLOCATABLE ::   TPK(:),FFT(:)
      REAL(SZ),ALLOCATABLE ::   FACET(:),ETRF(:)
      REAL(SZ),ALLOCATABLE ::   FFF(:),FFACE(:)
      REAL(SZ),ALLOCATABLE ::   ESBIN1(:),ESBIN2(:)
      REAL(SZ),ALLOCATABLE ::   TEMP_LV_A(:,:),TEMP_LV_B(:,:)
      REAL(SZ),ALLOCATABLE ::   QN0(:),QN1(:),QN2(:)
      REAL(SZ),ALLOCATABLE ::   EN0(:),EN1(:),EN2(:)
      REAL(SZ),ALLOCATABLE ::   QNAM(:,:),QNPH(:,:)
      REAL(SZ),ALLOCATABLE ::   ENAM(:,:),ENPH(:,:)
      REAL(SZ),ALLOCATABLE ::   QNIN1(:),QNIN2(:)
      REAL(SZ),ALLOCATABLE ::   ENIN1(:),ENIN2(:)
      !
      ! spatial interpolating factors for elevation stations
      real(sz),allocatable, target :: staie1(:),staie2(:),staie3(:)
      real(8), allocatable, target :: xev(:) ! x coords of current vel stations
      real(8), allocatable, target :: yev(:) ! y coords of current vel stations
      real(8), allocatable, target :: slev(:)! longitudes of current vel sta
      real(8), allocatable, target :: sfev(:)! latitudes of current vel sta

      ! spatial interpolating factors for velocity stations
      real(sz),allocatable, target ::   staiv1(:),staiv2(:),staiv3(:)
      real(8), allocatable, target ::   xec(:),yec(:),slec(:),sfec(:)
      ! spatial interpolating factors for concentration stations
      real(sz),allocatable ::   staic1(:),staic2(:),staic3(:)

      ! spatial interpolating factors for meteorological stations
      real(sz),allocatable, target ::   staim1(:),staim2(:),staim3(:)
      real(8), allocatable, target :: xem(:) ! x coords of met stations
      real(8), allocatable, target :: yem(:) ! y coords of met stations
      real(8), allocatable, target :: slem(:) ! longitudes of met stations
      real(8), allocatable, target :: sfem(:) ! latitudes of met stations

      REAL(SZ),ALLOCATABLE ::   TRANS_LV_B(:),TRANS_LV_A(:)
      REAL(SZ),ALLOCATABLE ::   SOURSIN(:)
      REAL(SZ),ALLOCATABLE ::   WSX1(:),WSY1(:),PR1(:)
      REAL(SZ),ALLOCATABLE ::   WSX2(:),WSY2(:)
      REAL(SZ),ALLOCATABLE ::   WVNX1(:),WVNY1(:),PRN1(:)
      REAL(SZ),ALLOCATABLE ::   WVNX2(:),WVNY2(:),PRN2(:)
      REAL(SZ),ALLOCATABLE ::   RSNX1(:),RSNY1(:),RSNX2(:),RSNY2(:)
      REAL(SZ),ALLOCATABLE ::   CICE1(:),CICE2(:)   !v49.64.01 tcm -- added
                                                    !for ice concentration
      REAL(SZ),ALLOCATABLE, TARGET :: RSNMax(:) ! v46.50 sb 11/11/2006

      REAL(SZ),ALLOCATABLE :: TK(:)
      REAL(8),ALLOCATABLE :: EMO(:,:),EFA(:,:)
      REAL(8),ALLOCATABLE,TARGET :: XEL(:)  ! x-coords of elevation stations
      REAL(8),ALLOCATABLE,TARGET :: YEL(:)  ! y-coords of elevation stations
      REAL(8),ALLOCATABLE,TARGET :: SLEL(:) ! latitudes of elevation stations
      REAL(8),ALLOCATABLE,TARGET :: SFEL(:) ! longitudes of elevation stations

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
      INTEGER  NBFR,NCOR,NE2,NP2

C    kmd48.33bc - added variables for the levels of no motion and top temperature boundary condition
      REAL(SZ),ALLOCATABLE ::   LNM_BC1(:), LNM_BC2(:), LNM_BC(:) ! level of no motion
      REAL(SZ),ALLOCATABLE ::   q_heat1(:), q_heat2(:), HFLUX(:)  ! top temperature
      REAL(SZ),ALLOCATABLE, TARGET :: Sponge(:)  ! sponge information for boundary conditions
      INTEGER,ALLOCATABLE, TARGET ::    NOFF(:), NODECODE(:)
      ! jgf52.08.08: Enable analyst to remove NOFF from wetting and
      ! drying algorithm; .true. by default to conform with prior
      ! ADCIRC versions. 
      LOGICAL :: noffActive = .true.
      INTEGER,ALLOCATABLE ::    NNODECODE(:)
      INTEGER,ALLOCATABLE ::    NOFFOLD(:)
      INTEGER,ALLOCATABLE ::    NIBCNT(:)
      INTEGER,ALLOCATABLE ::    NIBNODECODE(:)
      !
      ! lists of elements containing each station 
      integer, allocatable, target :: nnc(:) ! concentration station elements
      integer, allocatable, target :: nne(:) ! elevation station elements
      integer, allocatable, target :: nnv(:) ! velocity station elements
      integer, allocatable, target :: nnm(:) ! meteorological station elements
      !
      INTEGER,ALLOCATABLE ::    IWKSP(:),IPARM(:),IPV(:)
C
C.....for buoyancy forcing
      REAL(SZ),ALLOCATABLE ::   DASigT(:), DASalt(:), DATemp(:)
      REAL(SZ),ALLOCATABLE ::   VIDBCPDXOH(:), VIDBCPDYOH(:)
C.....for scalar transport
      REAL(SZ),ALLOCATABLE ::   DAConc(:)
C.....for wet/dry
      REAL(SZ)                  HABSMIN,HOFF
cjjw...added 2 lines
C.....for matrix conditioning,global to allow for reading in from hotstart file
      REAL(SZ)                  EP !jgf45.08 EP global for predictor-corrector

      REAL(SZ) DT
C
      INTEGER NOLICA,NOLIFA,IHOT,IDen
      INTEGER :: NSCREEN = 1
C...
C...DECLARE REAL(8) AND CHAR VARIABLES, EQUIVALENCES
C...
      REAL(8) STATIM,REFTIM,DTDP
      REAL(8) WREFTIM,WTIMED,WTIME2,WTIME1,WTIMINC,QTIME1,QTIME2
      REAL(8) FTIMINC,ETIMINC,RSTIME1,RSTIME2,RSTIMINC

      REAL(8) DIST,DELDIST,DELETA
      REAL(8),ALLOCATABLE :: AMIG(:),AMIGT(:),FAMIG(:)
      REAL(8),ALLOCATABLE :: PER(:),PERT(:),FPER(:)
      REAL(8),ALLOCATABLE :: FREQ(:),FF(:),FACE(:)


      REAL(8) DTDPHS  !kmd48.33 - added time step for hot start
      CHARACTER(80)             :: RUNDES, RUNID, AFRIC
      CHARACTER(4)              :: RDES4(20),RID4(20)
      CHARACTER(8)              :: RDES8(10),RID8(10),AID8(10)
      CHARACTER(10)             :: ALPHA ! periodic constituent name
      CHARACTER(10),ALLOCATABLE :: ELEVALPHA(:) ! periodic constituent name
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
      ! This is .true. if the time step size of the current run is
      ! different from the time step size in the hotstart file; also
      ! true for IDEN=1, 3D baroclinic runs 
      LOGICAL :: CHOTHS              = .FALSE.

! kmd - Evan's updates for rivers above MSL
      LOGICAL :: River_above_MSL     = .FALSE.

C     jgf48.4627 When only meteorological output is requested, skip past
C     the GWCE and Momentum equations, so that ADCIRC runs faster and
C     only calculates the requested quantities.
      LOGICAL :: METONLY             = .FALSE.
C
C     jgf51.14: Moved the parameters that control output timing 
C     from read_input.F to here so that they can be used to recalculate
C     starting and ending time step for output when the time step
C     size has changed after a hotstart. 
      REAL(SZ) TOUTSC
      REAL(SZ) TOUTSE,TOUTFE
      REAL(SZ) TOUTSV,TOUTFV
      REAL(SZ) TOUTFGC, TOUTFGE, TOUTFGV, TOUTFGW, TOUTF_TVW
      REAL(SZ) TOUTFM, TOUTSGC, TOUTSGE, TOUTSGV, TOUTSGW, TOUTSM, TOUTS_TVW

      REAL(SZ),ALLOCATABLE ::   RPIPEWL1AVG(:),RPIPEWL2AVG(:)
      INTEGER  NODEDRYMIN,NODEWETMIN
      INTEGER  IM
      INTEGER :: NABOUT = 0
      INTEGER  NFFR
      INTEGER  NFOVER,NHG,NHY,NOLICAT,NOUTC
      INTEGER  NOUTE,NSPOOLE,NOUTV,NSPOOLV
      INTEGER  NRAMP,NRS,NSTAE,NSTARTDRY,NSTAV,NT,NTCYFE
      INTEGER  NTCYFV,NTCYSE,NTCYSV,NTIF,NTIP,NTRSPE,NTRSPV
      INTEGER  NTCYS_TVW,NTCYF_TVW
      INTEGER  NWLAT,NWLON,NWS
      INTEGER  ICSTP, IDSETFLG, IER

Corbitt 120322: Create Elemental Versions of IFNLCAT and IFNLCT
      INTEGER  IFNLCAT, IFNLCT, IFNLFA, IFNLCTE, IFNLCATE
      INTEGER  IFWIND
      ! 
      ! fulldomain (a.k.a. global) file position counters
      integer, target :: igcp ! concentration
      integer, target :: igep ! elevation
      integer, target :: igpp ! barometric pressure
      integer, target :: igvp ! water velocity
      integer, target :: igwp ! wind velocity
      integer, target :: itau0p ! position in tau0 output file for variable tau0
      integer, target :: igsp ! sponge layer output
      integer, target :: igdp ! debug elevation
      ! station file position counters
      integer, target :: ivstp ! file position in velocity station file
      integer, target :: iwstp ! file position in wind velocity station file
      integer, target :: ipstp ! position in barometric pressure station file
      integer, target :: iestp ! file position in station elevation file 
      integer, target :: ig_tvw ! file position in weir elevation file
      integer, target :: nodeCodeOutFilePos
      integer, target :: noffOutFilePos
      ! jgf52.08.03: Some files (e.g., initiallydry.63 and endrisinginun.63)
      ! are only written once, but are written using the same subroutines
      ! that write time varying output files; the output subroutine increments
      ! the file position counter; so we have to have a fake file position 
      ! counter even for output files that are only written once.
      integer, target :: fauxFilePosition 

      INTEGER  IHOTSTP, ILUMP, IMHS
      INTEGER  IREFYR,IREFMO,IREFDAY,IREFHR,IREFMIN, ISLDIA
      INTEGER  ITIME, ITEMPSTP, ITEST
      INTEGER  ITHS, ITITER, ITMAX    
      INTEGER IWTIME, IWTIMEP
      
      INTEGER  IWYR, J12, J13, J21, J23, J31, J32
      INTEGER  JJ, KEMAX, KVMAX, LRC, LUMPT
      INTEGER  MBW, MDF, NBDJ, NBNCTOT, NBW
      INTEGER  NCBND, NCTOT, NDRY, NDSETSC
      INTEGER  NCChange !jgf45.06 flag indicating wetting / drying has occurred
      INTEGER  NDSETSE, NDSETSV, NDSETSW, NDSETS_TVW,NEle
      INTEGER  NHSINC, NHSTAR
!  tcm v50.75 removed ifdef cswan to allow for use whenever nrs=3 or nrs=4
!#ifdef CSWAN 
Casey 090302: Added a variable for output of radiation stress gradients.
      INTEGER, TARGET :: IGRadS
!#endif
      INTEGER  NNBB, NNBB1, NNBB2
      INTEGER  NNBB1R, NNBB1L, NNBB2R, NNBB2L ! sb46.28.sb05.05 11/01/2006
      INTEGER  NNBB1WN, NNBB2WN               ! sb46.28.sb05.05 WN stands for wet neighbors 11/01/2006
      INTEGER  NOUTGC, NOUTGE, NOUTGV, NOUTGW, NOUTM, NOUT_TVW
      INTEGER, TARGET :: NSCOUC, NSCOUE, NSCOUGC, NSCOUGE, NSCOUGV ! output spool counters
      INTEGER, TARGET :: NSCOUGW, NSCOURS, NSCOUGP, NSCOUM, NSCOUV ! output spool counters
      INTEGER, TARGET :: NSCOUGI, NSCOUGT, NSCOUGS  ! output spool counters
      INTEGER, TARGET :: NSCOU_TVW ! output spool counters
      INTEGER, TARGET :: NSCOUP ! output spool counter for fort.71 (bar. press. sta.)
      INTEGER, TARGET :: NSCOUI ! output spool counter for fort.91 (ice sta.)     
      integer, target :: nodeCodeOutSpoolCounter
      integer, target :: noffOutSpoolCounter
      integer, target :: initiallyDrySpoolCounter
      
      INTEGER :: NSPOOLC, NSPOOLGC, NSPOOLGE ! output increments
      INTEGER :: NSPOOLGV, NSPOOLGW, NSPOOLM, NSPOOL_TVW ! output increments
      INTEGER  NSTAC, NSTAM, NTCYFC, NTCYFGC, NTCYFGE
      INTEGER  NTCYFGV, NTCYFGW
      INTEGER  NTCYFM, NTCYSC, NTCYSGC, NTCYSGE
      INTEGER  NTCYSGV, NTCYSGW, NTCYSM
      INTEGER  NTRSPC, NTRSPM, NUMITR, NW, NWET, NWSEGWI, NWSGGWI
      INTEGER  NP_GLOBAL
C     kmd48.33bc - added variables for 3D boundary conditions
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

      REAL(SZ) RIVBCTIMINC, RIVBCTIME1, RIVBCTIME2  ! time information for river - baroclinic
      REAL(SZ) BCRivRATIO, RIVBCSTATIM

      REAL(SZ) ADVECX, ADVECY, ARG
      REAL(SZ) ARG1, ARG2, ARGJ, ARGJ1, ARGJ2
      REAL(SZ) AUV21, AUV22
      REAL(SZ) BEDSTR, BNDLEN2O3NC
      REAL(SZ) CELERITY, COND, CONVCR
      REAL(SZ) DT2, DTO2, DTOHPP
      REAL(SZ) ETIME1, ETIME2, ETRATIO
      REAL(SZ) GO3, GO2, GA00, GB00A00, GC00, GDTO2, GFAO2

      REAL(SZ) P11, P22, P33
      REAL(SZ) QFORCEI, QFORCEJ
      REAL(SZ) QUNORM, QUTemp, QVNORM, RAMP, RAMP1, RAMP2, RBARWL
      REAL(SZ) RBARWL1,RBARWL1F, RBARWL2, RBARWL2F, RFF, RFF1, RFF2
      REAL(SZ) SADVDTO3
      REAL(SZ) TADVODT
      REAL(SZ) TIMEIT, TOUTFC
      REAL(SZ) TKWET
      REAL(SZ) VEL, VELABS
      REAL(SZ) VELMAX

      REAL(SZ) A00,B00,C00
      REAL(SZ) CORI
      REAL(SZ) DAY,DRAMP,DUM1,DUM2
      REAL(SZ) H0
      REAL(SZ) RNDAY
      REAL(SZ) WLATMAX
      REAL(SZ) WLONMIN,WLATINC,WLONINC, VELMIN
      REAL(8)  RNP_GLOBAL
      REAL(8)  REFSEC   ! required to run in either 32-bit or 64-bit

      LOGICAL :: useNetCDF = .false. ! .true. if netcdf is used at all
      LOGICAL :: useNetCDFOutput = .false. ! only for output, not hotstart
      LOGICAL :: useXDMF = .false.
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

#ifdef ADCNETCDF
      LOGICAL, PARAMETER :: netcdf_avail = .true. 
#else
      LOGICAL, PARAMETER :: netcdf_avail = .false. ! .true. if NetCDF support was compiled in
#endif
      
! CF Finish addition for NETCDF I/O

! TCM v49.64.01 Addition for Ice Concentration Fields
      INTEGER :: NCICE
      INTEGER, TARGET :: IGIP ! position in fulldomain Ice output file
      INTEGER, TARGET :: IICESTP ! position in station Ice output file
      REAL(8) :: CICE_TIMINC,CICE_TIME1,CICE_TIME2

! TCM V50.66.01 ADDITIONS FOR TIME VARYING BATHYMETRY
      INTEGER :: NDDT  !FLAG TO USE TIME VARYING BATHYMETRY
      REAL(8) :: BTIME1,BTIME2  !TIMES FOR TIME VARYING BATHYMETRY Records
      REAL(8) :: BTIME_END  !time that 
      !variables related to global bathymetry output file
      INTEGER, TARGET :: IGBP ! position in global bathy output file
      INTEGER :: NDSETSB,NTCYFGB,NTCYSGB,NTRSPB
      INTEGER, TARGET :: NSCOUGB
      !variables related to station bathymetry output file
      INTEGER :: NTCYSB,NTCYFB
      INTEGER, TARGET :: NSCOUB
      INTEGER, TARGET :: IBSTP  ! position in station bathy output file      
      REAL(8) :: BTIMINC ! TIME INCREMENT (SECONDS) FOR BATHYMETRY CHANGES USED FOR TIME VARYING BATHYMETRY
      REAL(8) :: BCHGTIMINC ! time increment (seconds) over which bathymetry changes during a btimeinc interval
      REAL(SZ),ALLOCATABLE, TARGET ::  DP00(:)

! TCM v52.30.01 ADDITIONS FOR USING SMAGORINSKY eddy MODEL OPTIONS
      REAL(SZ) :: SMAG_UPPER_LIM = 100.d0  ! m^2/s (default upper limit)
      REAL(SZ) :: SMAG_LOWER_LIM = 1.0d-8  ! m^2/s (default lower limit)
      LOGICAL :: SMAG_COMP_FLAG = .FALSE.   !BY DEFAULT DO NOT USE SMAG
      ! 
      ! all info needed for self describing dataset in XDMF
      type xdmfMetaData_t
         integer*8 :: xdmfFortranObj  ! object associated with the file
         integer :: depthID    ! index of the depth attribute
         integer :: geometryID ! index of the Geometry data in the output file
         real(8), allocatable :: tempCoord(:,:) ! holds Geometry data
         integer :: topologyID ! index of the Topology data in the output file
         integer, allocatable :: xdmf_nm(:,:)  ! 0-offset connectivity array
         logical :: createdIDs        ! .true. if infoIDs have been created 
         character(80) :: variable_name 
         integer :: variable_name_id
         character(80) :: long_name 
         integer :: long_name_id
         character(80) :: standard_name
         integer :: standard_name_id
         character(80) :: coordinates
         integer :: coordinates_id
         character(80) :: units
         integer :: units_id
         character(80) :: positive
         integer :: positive_id
         real(sz), allocatable :: data_array3(:,:)
      end type xdmfMetaData_t
      !
      ! info required for dedicated writers
      type writerData_t
         logical :: bufferInitialized  !.true. if memory has been allocated
         logical :: bufferLoaded       !.true. if buffered data must be written
         real(sz) :: myTime
         integer :: myTimeStep
      end type writerData_t
C
C     jgf48.03 the following data structure is used to describe
C     the size and shape of the output data.
C     TCM48.4618 Added ConsiderWetDry, alternate_value and field_name 
C     elements to the structure
      type OutputDataDescript_t
        integer           :: lun
        CHARACTER(1024)   :: file_name ! name of output file
        character(20)     :: field_name !Name of the field associated with the data structure
        integer           :: specifier  ! output format from fort.15 (NOUTE,etc)
        REAL(SZ)          :: initial_value
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
        logical :: isInteger ! .true. if this is an integer array
        INTEGER, pointer  :: iarray_g(:) ! full domain integer
        REAL(SZ), pointer :: array_g(:)  ! full domain data scalar etc
        REAL(SZ), pointer :: array2_g(:) ! full domain 2nd vec component
        REAL(SZ), pointer :: array3_g(:) ! full domain 3rd vec component
        REAL(SZ), pointer :: array2D_g(:,:) ! fulldomain matrix
        REAL(SZ), pointer :: hotstart(:) ! from hotstart (used by min max)
        REAL(SZ), pointer :: hotstart_g(:) ! from hotstart (used by min max)
        logical :: writeFlag ! .true. if the data should be written during this run
        logical :: isStation ! .true. if the data represent a recording station
        logical :: divideByDepth ! .true. if the solution must be divided by depth b/f output
        integer, allocatable :: writerFormats(:) ! formats that dedicated writers can write
        logical :: useWriter ! .true. if a dedicated writer should be used to write these data
        integer :: file_extension  ! end of the file name       
        character(len=32) :: file_basename  ! usually 'fort' but sometimes 'rads' etc
        real(sz), pointer :: interped_array(:) ! array that stations spatially interpolate
        real(sz), pointer :: interped_array2(:) ! array2 that stations spatially interpolate
        integer :: startTimeStep ! time step after which output should begin
        integer :: endTimeStep   ! last time step when output should occur    
        integer :: outputTimeStepIncrement ! num time steps between outputs          
        integer, pointer :: spoolCounter ! counts time steps between outputs
        integer, pointer :: filepos ! keeps track of position in output files
        integer, pointer :: elements(:) ! array of elements that contain the stations
        real(sz), pointer :: interp_fac1(:) ! array of interp. factors for stations
        real(sz), pointer :: interp_fac2(:) ! array of interp. factors for stations        
        real(sz), pointer :: interp_fac3(:) ! array of interp. factors for stations
        logical :: initialized ! .true. if initialization has occurred (netcdf, xdmf)
        logical :: minmax_timestamp ! .true. if a time stamp will be written to a max file
        logical :: readMaxMin ! .true. if max/min file should be read upon hotstart
        type(xdmfMetaData_t) :: xdmfMD ! used if this data will be written to XDMF
        type(writerData_t) :: wMD ! used if this data will be written by dedicated writer
      end type OutputDataDescript_t
      !
      ! inunThresh is used to set a minimum inundation depth for 
      ! accumulating inundationTime for detailed inundation output
      real(sz) :: inunThresh = 0.6d0
      !
      ! inundationOutput can be set to .true. in read_input()
      logical :: inundationOutput = .false. 

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
      CHARACTER(len=1024) :: scratchFormat  ! used for Fortran format strings
      INTEGER :: sourceNumber ! index into messageSources for current sub

C     v49.48.01 tcm 20110126 Variables related to kdtree searches
      INTEGER :: SRCHDP=100
      TYPE(KDTREE2), POINTER :: TREE
      TYPE(KDTREE2_RESULT), ALLOCATABLE :: KDRESULTS(:)

      CONTAINS

C--------------------------------------------------------------------
C     S U B R O U T I N E    I N I T   L O G G I N G
C--------------------------------------------------------------------
C     jgf49.44: Initialize the names for the logging levels and the counter
C     for the current subroutine.
C--------------------------------------------------------------------
      SUBROUTINE initLogging()
      IMPLICIT NONE
C...
C...  OPEN STATEMENT FOR UNIT 16 OUTPUT FILE (ADCIRC LOG FILE)
C...
      sourceNumber = 0
      logLevelNames(1) = "DEBUG"
      logLevelNames(2) = "ECHO"
      logLevelNames(3) = "INFO"
      logLevelNames(4) = "WARNING"
      logLevelNames(5) = "ERROR"

C--------------------------------------------------------------------
      END SUBROUTINE initLogging
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C     S U B R O U T I N E    O P E N   L O G   F I L E
C--------------------------------------------------------------------
C     jgf50.65: Open the log file; this must be called after make dirname
C     so that we know where to put this log file.
C--------------------------------------------------------------------
      SUBROUTINE openLogFile()
      use sizes, only : localDir
      IMPLICIT NONE
C...
C...  OPEN STATEMENT FOR UNIT 16 OUTPUT FILE (ADCIRC LOG FILE)
C...
C     cms51.06: moved fort.33 file from read_input.F to here
#ifdef CMPI
      OPEN(16,FILE=trim(localdir)//'/'//'fort.16', ACTION='WRITE',
     &    STATUS='REPLACE')
      OPEN(33,FILE=trim(localdir)//'/'//'fort.33', ACTION='WRITE',
     &    STATUS='REPLACE')
#else
      OPEN(16,FILE='fort.16', ACTION='WRITE', STATUS='REPLACE')
      OPEN(33,FILE='fort.33', ACTION='WRITE', STATUS='REPLACE')
#endif
C--------------------------------------------------------------------
      END SUBROUTINE openLogFile
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
      SUBROUTINE screenMessage(level, message)
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
     &             trim(logLevelNames(level+2)),
     &             (trim(messageSources(j)),j=1,sourceNumber)
               write(screenUnit,332) trim(message)
#else
               write(screenUnit,333) trim(logLevelNames(level+2)),
     &             trim(messageSources(sourceNumber)), trim(message)
#endif
!#ifdef FLUSH_MESSAGES
               flush(screenUnit)
!#endif
            ENDIF
         ENDIF
      ENDIF
331   FORMAT(A,": ",A,50(:,"->",A))
332   FORMAT(": ",A)
333   FORMAT(A,": ",A,": ",A)
C--------------------------------------------------------------------
      END SUBROUTINE screenMessage
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
!#ifdef FLUSH_MESSAGES
               flush(16)
!#endif
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
      call screenMessage(level, message)
      call logMessage(level, message)

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

C-----------------------------------------------------------------------
C     S U B R O U T I N E   O P E N  F I L E  F O R  R E A D
C-----------------------------------------------------------------------
C     jgf50.16 Added general subroutine for opening an existing input
C     file for reading. Includes error checking.
C-----------------------------------------------------------------------
      SUBROUTINE openFileForRead(lun, filename, errorIO)
      IMPLICIT NONE
      INTEGER, intent(in) :: lun   ! fortran logical unit number
      CHARACTER(*), intent(in) :: filename ! full pathname of file
      INTEGER, intent(out) :: errorIO  ! zero if the file opened successfully
      LOGICAL :: fileFound    ! .true. if the file is present
C
      call setMessageSource("openFileForRead")
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      errorIO = 0
C
C     Check to see if file exists
      write(scratchMessage,21) lun
  21  format("Searching for file to open on unit ",i0,"...")
      call logMessage(INFO,trim(scratchMessage))
      inquire(FILE=trim(filename),EXIST=fileFound)
      if (fileFound.eqv..false.) then
         write(scratchMessage,23) trim(filename)
  23     format("The file '",A,"' was not found.")
         call allMessage(INFO,scratchMessage)
         errorIO = 1
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         return  ! file not found
      else
         write(scratchMessage,24) trim(filename)
  24     format("The file '",A,"' was found. The file will be opened.")
         call logMessage(INFO,trim(scratchMessage))
      endif
C
C     Open existing file
      OPEN(lun,FILE=trim(filename),STATUS='OLD',
     & ACTION='READ',IOSTAT=errorIO)
      if (errorIO.ne.0) then
          write(scratchMessage,25) trim(filename)
  25      format("Could not open the file '",A,"'.")
          call allMessage(ERROR,trim(scratchMessage))
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
          call allMessage(DEBUG,"Return.")
#endif
          call unsetMessageSource()
          return  ! file found but could not be opened
      else
         write(scratchMessage,26) trim(filename)
  26     format("The file '",A,"' was opened successfully.")
         call logMessage(INFO,trim(scratchMessage))
      endif

#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
C-----------------------------------------------------------------------
      END SUBROUTINE openFileForRead
C-----------------------------------------------------------------------

C     ----------------------------------------------------------------
C      F U N C T I O N   S P H E R I C A L   D I S T A N C E
C     ----------------------------------------------------------------
C     jgf49.1001 Function to get the distance along the surface of
C     a sphere (the earth's surface in this case).
C     ----------------------------------------------------------------
      REAL(SZ) FUNCTION sphericalDistance(dx, dy, y1, y2)
      IMPLICIT NONE
      REAL(SZ), intent(in) :: DX    ! longitude distance in radians
      REAL(SZ), intent(in) :: DY    ! latitude distance in radians
      REAL(SZ), intent(in) :: y1    ! degrees latitude of starting point
      REAL(SZ), intent(in) :: y2    ! degrees latitude of ending point
C
C     compute the distances based on haversine formula for
C     distance along a sphere
      sphericalDistance = Rearth *
     &  ( 2.0d0*ASIN(
     &     sqrt(sin(DY/2.0d0)**2.0d0
     &      + cos(y1*DEG2RAD)*cos(y2*DEG2RAD)*sin(DX/2.0d0)**2.0d0
     &   )))

      RETURN
C     ----------------------------------------------------------------
      END FUNCTION sphericalDistance
C     ----------------------------------------------------------------

C...
C...Allocate space for Arrays dimensioned by MNE and MNP
C...
      SUBROUTINE ALLOC_MAIN1a()
      IMPLICIT NONE
C  
      !jgfdebug TODO: rationalize initialization of these variables
      ALLOCATE ( ETAS(MNP))
      ALLOCATE ( GWCE_LV(MNP))
      ALLOCATE ( UU1(MNP),VV1(MNP))
      uu1(:) = 0.d0
      vv1(:) = 0.d0
      ALLOCATE ( QX1(MNP),QY1(MNP))
      qx1(:) = 0.d0
      qy1(:) = 0.d0
      ALLOCATE ( NNODECODE(MNP))
      nnodecode(:) = 0
      ALLOCATE ( NOFF(MNE),NOFFOLD(MNE))
      noff(:) = 0
      noffold(:) = 0
      ALLOCATE ( NIBNODECODE(MNP))
      ALLOCATE ( CH1(MNP),TRANS_LV_B(MNP),TRANS_LV_A(MNP))
      ALLOCATE ( SOURSIN(MNP))
      ALLOCATE ( TK(MNP))

      ALLOCATE ( UU2(MNP),VV2(MNP))
      uu2(:) = 0.d0
      vv2(:) = 0.d0
      ALLOCATE ( QX2(MNP),QY2(MNP))
      qx2(:) = 0.d0
      qy2(:) = 0.d0
      ALLOCATE ( ETA1(MNP),ETA2(MNP))
      eta1(:) = 0.d0
      eta2(:) = 0.d0

      ALLOCATE ( CORIF(MNP))
      ALLOCATE ( MOM_LV_X(MNP),MOM_LV_Y(MNP))

      ALLOCATE ( NODECODE(MNP))
      ALLOCATE ( NIBCNT(MNP) )
      ALLOCATE ( DASigT(MNP),DATemp(MNP),DASalt(MNP))
      ALLOCATE ( VIDBCPDXOH(MNP), VIDBCPDYOH(MNP))
      ALLOCATE ( DAConc(MNP))
      ALLOCATE ( LSXX(MNP),LSXY(MNP),LSYX(MNP),LSYY(MNP))
ckmd   Added in parameters for the pc algorithm
      ALLOCATE ( ETA0(MNP), UU0(MNP), VV0(MNP))
      ALLOCATE ( QX0(MNP), QY0(MNP), ETAS0(MNP))
      ALLOCATE ( TK0(MNP), TK2(MNP))

C    kmd48.33bc - arrays for baroclinic
      ALLOCATE (q_heat1(MNP), q_heat2(MNP), HFLUX(MNP) )
      ALLOCATE (sponge(MNP) )

C    Allocate array for ESL i/o - only needed if
C    loadEleSlopeLim.eqv..true
      ALLOCATE (ESLONOFF(MNP))
      ESLONOFF = 0.d0

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

      ALLOCATE(EtaDisc(MNP))
      EtaDisc(:) = 0.0d0

      ALLOCATE ( ETAMAX(MNP)) ! v46.50 sb 11/11/2006
      ALLOCATE ( ETAMAX_TIME(MNP)) ! tcm v51.20.01 
      ALLOCATE ( UMAX(MNP)) ! v46.50 sb 11/11/2006
      ALLOCATE ( UMAX_TIME(MNP)) ! tcm v51.20.01 
      ETAMAX(:) = -99999.d0 ! v46.50 sb 11/11/2006
      EtaMax_Time(:) = 0.d0
      UMAX(:) = 0.d0 ! v46.50 sb 11/11/2006
      UMax_Time(:) = 0.d0





      RETURN
      END SUBROUTINE ALLOC_MAIN1a





C...
C...Allocate space for Arrays dimensioned by MNOPE and MNETA
C...
      SUBROUTINE ALLOC_MAIN2()
      IMPLICIT NONE
C
      ALLOCATE ( ESBIN1(MNETA),ESBIN2(MNETA))

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
      ALLOCATE ( QNIN1(MNVEL),QNIN2(MNVEL))
C
C     jgf46.21 Added support for IBTYPE=52.
      ALLOCATE ( ElevDisc(MNVEL))
C
C     jgf49.58 Explicitly initialize variables that would otherwise be
C     used without initialization.
      EN2(:) = 0.0D0

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
      ALLOCATE ( EMO(MNBFR,MNETA),EFA(MNBFR,MNETA))
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

c...  tcm v50.66.01 added for time varying bathymetry
      IF (NDDT.NE.0) ALLOCATE (DP00(MNSTAE) )
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
      WSX1(:) = 0.0d0
      WSY1(:) = 0.0d0
      PR1(:) = 0.0d0
      WSX2(:) = 0.0d0
      WSY2(:) = 0.0d0
      PR2(:) = 0.0d0

      ALLOCATE ( RSNX1(MNP),RSNY1(MNP),RSNX2(MNP),RSNY2(MNP) )
      RSNX1(:) = 0.0d0 !jgf51.46: Added initialization of RSNX1, RSNY1
      RSNY1(:) = 0.0d0
      RSNX2(:) = 0.d0
      RSNY2(:) = 0.d0

      if (NCICE.NE.0) then
         ALLOCATE ( CICE1(MNP),CICE2(MNP),CICEOUT(MNP) ) !tcm v49.64.01 added
                                                        ! for ice concentration
         CICE1(:) = 0.0d0 !tcm v49.64.01 added for ice concentration
         CICE2(:) = 0.0d0 !tcm v49.64.01 added for ice concentration
         CICEOUT(:) = 0.0d0 !tcm v49.64.01 added for ice concentration
      ENDIF

      ALLOCATE ( WVNXOUT(MNP),WVNYOUT(MNP) )
      WVNXOUT(:) = 0.d0
      WVNYOUT(:) = 0.d0

!  tcm v50.75 removed ifdef cswan to allow for use whenever nrs=3 or nrs=4
      IF ((ABS(NRS).EQ.3).or.(ABS(NRS).EQ.4)) then
!#ifdef CSWAN
Casey 090302: Added the next line for output of radiation stress gradients.
         ALLOCATE ( RSNXOUT(MNP),RSNYOUT(MNP) )
         RSNXOUT(:) = 0.d0  !tcm v50.75 added initialization to zero
         RSNYOUT(:) = 0.d0
      ENDIF
!#endif
      ALLOCATE ( PrMin(MNP)) ! v46.50 sb 11/11/2006
      PrMin(:) = 99999.d0
      ALLOCATE (PrMin_Time(MNP) ) !tcm v51.20.01
      PrMin_Time(:) = 0.d0
      ALLOCATE ( WVNOUTMax(MNP)) ! v46.50 sb 11/11/2006
      WVNOUTMAX(:) = 0.d0
      ALLOCATE (WVNOutMax_Time(MNP) ) !tcm v51.20.01
      WVNOutMax_Time(:) = 0.d0
      ALLOCATE ( RSNMax(MNP))               ! v46.50 sb 11/11/2006
      RSNMax(:) = 0.d0
      ALLOCATE (RSNMax_Time(MNP) ) !tcm v51.20.01
      RSNMax_Time(:) = 0.d0

C
      RETURN
      END SUBROUTINE ALLOC_MAIN12

! TCM V50.66.01 ADDITIONS FOR TIME VARYING BATHYMETRY (ALLOC_MAIN13)
C...
C...Allocate space for TIME VARYING BATHYMETRY
C...
      SUBROUTINE ALLOC_MAIN13()
      IMPLICIT NONE
      
      ALLOCATE ( DP1(MNP), DP2(MNP) )
      DP1(:) = 0.D0
      DP2(:) = 0.D0

      END SUBROUTINE ALLOC_MAIN13

      !-----------------------------------------------------------------
      !                      S U B R O U T I N E   
      !                G E T    F O R M A T    I N T E G E R
      !-----------------------------------------------------------------
      ! Sets file format integer that corresponds to the file 
      ! format string.
      !-----------------------------------------------------------------
      subroutine getFormatInteger(typeString, typeHolder)
      implicit none
      character(len=20), intent(in) :: typeString
      integer, intent(out) :: typeHolder
      !
      select case(trim(typeString))
      case('OFF')
         typeHolder = 0 
      case('ASCII')
         typeHolder = 1
      case('BINARY')
         typeHolder = 2 
      case('NETCDF3')
         typeHolder = 3 
      case('SPARSE_ASCII')
         typeHolder = 4 
      case('NETCDF4')
         typeHolder = 5 
      case('XDMF')
         typeHolder = 7
      case default
         write(6,'("WARNING: Unrecognized format ",i0,".")') trim(typeString)
      end select
      !----------------------------------------------------------------------
      end subroutine getFormatInteger
      !----------------------------------------------------------------------


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
      ! v50.36    Chris Massey   Added logic to parse to handle tabs and commas
      !                          and not just spaces.
      !=================================================================
      subroutine parse (line, lvar)
      IMPLICIT NONE
      CHARACTER(132) LINE,L
      CHARACTER(50) LVAR(3)
      INTEGER LS, L1, LEN1, I
      INTEGER L2,L3,L4,PB
      INTEGER NAMSTART  ! position of beginning of station name
      INTEGER NAMEND    ! position of end of station name
      INTEGER BANGEND   ! location of terminating "!", if any
      CHARACTER(len=1) :: tab_character
C
      call setMessageSource("parse")
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      tab_character = achar(9)
      LS=1
      L=TRIM(ADJUSTL(LINE))
      CALL COMPACT(L)
      LEN1=LEN(L)
!     EXTRACT STATION DATA FROM LINE
      L2 = 0
      L3 = 0
      L4 = 0
      DO I=1,2
!        DO I=1,3
         if (index(L," ").ne.0) then
            L2=LEN(TRIM(ADJUSTL(L(1:INDEX(L," ")-1))))          !space
         endif
         if (index(L,tab_character).ne.0) then
            L3=LEN(Trim(adjustl(L(1:index(L,tab_character)-1))))     !tab
         endif
         if (index(L,",").ne.0) then
            L4=LEN(trim(adjustl(L(1:index(L,",")-1))))          !comma
         endif
         !Determine starting position based on which delimiter was found first
         L1 = L2
         PB = 1
         IF(L1.GT.L3.AND.L3.NE.0) THEN
            L1 = L3
            PB = 2
         ENDIF
         IF(L1.GT.L4.AND.L4.NE.0) THEN
           L1=L4
           PB = 2
         ENDIF

         LVAR(I)=ADJUSTL(L(LS:L1))
         L=ADJUSTL(L(L1+1:LEN(L)))
         CALL COMPACT(L)

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
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      !=================================================================
      END SUBROUTINE parse
      !=================================================================

      !=================================================================
       SUBROUTINE COMPACT(STR)
C     tcm50.36 SUBROUTINE TO CONVERT MULTIPLE SPACES AND TABS TO SINGLE
C              SPACES; DELETES CONTROL CHARACTERS;REMOVES INITIAL SPACES.
C              Adapted from stringmod.f90 from the website:
C              http://www.gbenthien.net/strings/index.html
      !=================================================================
       IMPLICIT NONE

       INTEGER :: I,K,ISP,LENSTR,ICH
       CHARACTER(LEN=*):: STR
       CHARACTER(LEN=1):: CH
       CHARACTER(LEN=LEN_TRIM(STR)):: OUTSTR

      call setMessageSource("compact")
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
       STR=ADJUSTL(STR)
       LENSTR=LEN_TRIM(STR)
       OUTSTR=' '
       ISP=0
       K=0

       DO I=1,LENSTR
          CH=STR(I:I)
          ICH=IACHAR(CH)

          SELECT CASE(ICH)

          CASE(9,32)     ! SPACE OR TAB CHARACTER
             IF(ISP==0) THEN
                K=K+1
                OUTSTR(K:K)=' '
             END IF
             ISP=1

          CASE(33:)      ! NOT A SPACE, QUOTE, OR CONTROL CHARACTER
             K=K+1
             OUTSTR(K:K)=CH
             ISP=0

          END SELECT

       END DO

       STR=ADJUSTL(OUTSTR)
C
#if defined(GLOBAL_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      !=================================================================
       END SUBROUTINE COMPACT
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
      res = 1
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
