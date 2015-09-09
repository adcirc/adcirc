C***********************************************************************
C PADCIRC VERSION 45.12 03/17/2006                                     *
C last changes in this file VERSION 45.05                              *
C                                                                      *
C***********************************************************************
C
       MODULE WIND
       USE SIZES, ONLY : SZ, LOCALDIR
       USE GLOBAL, ONLY : ScreenUnit, setMessageSource, allMessage,
     &    screenMessage, logMessage, unsetMessageSource, 
     &    DEBUG, ECHO, INFO, WARNING, ERROR, scratchMessage,
     &    Rearth, pi, deg2rad, rad2deg, e, PRBCKGRND, g, rhoWat0,
     &    sphericalDistance, rhoAir, one2ten, windReduction,
     &    openFileForRead
#ifdef CMPI
      use messenger, only : msg_fini
#endif
      integer, parameter :: allowable_nws(17)
     &      = (/0,1,2,3,4,5,6,7,8,10,11,12,15,16,19,20,29/)
      logical nwsOK     ! true if NWS value in fort.15 is allowed

Casey 110518: Added for Mark Powell's sectorial wind drag.
      REAL(SZ) :: EyeLatR(3) = 0.D0
      REAL(SZ) :: EyeLonR(3) = 0.D0
      LOGICAL  :: FoundEye
C
      REAL(8), PARAMETER :: TWOPI=PI*2.D0
      REAL(8), PARAMETER :: HFPI=PI/2.D0
      REAL(8), PARAMETER :: mperdeg  = Rearth * pi / 180.0d0

      INTEGER StormNumber !jgf46.28 For Holland Wind model wind multiplier
      REAL(8) WindRefTime !jgf46.29 seconds since beginning of year, this
                          !corresponds to time=0 of the simulation
      REAL(8) BLAdj       !jgf46.32 boundary layer adjustment for Holland model

      ! jgf49.1001 Variables for use in NWS=29, for embedding an asymmetric
      ! vortex from NWS19 into an OWI basin scale met field from NAM
      ! using NWS12
      REAL(SZ) vortexRMW  ! Rmax of vortex data, nautical miles
      REAL(SZ) vortexLat  ! @center of vortex, degrees west latitude
      REAL(SZ) vortexLon  ! @center of vortex, degrees north longitude
      REAL(SZ), ALLOCATABLE :: vortexWVNX2(:) ! u wind vel from vortex model
      REAL(SZ), ALLOCATABLE :: vortexWVNY2(:) ! v wind vel from vortex model
      REAL(SZ), ALLOCATABLE :: vortexPRN2(:)  ! bar. press. from vortex model
      !
      ! radial extent of unblended vortex winds, normalized to current Rmax
      REAL(SZ) pureVortex
      !
      ! radial distance to point where vortex effects insignificant,
      ! normalized to current Rmax
      REAL(SZ) pureBackground

      REAL(SZ) :: WindDragLimit          ! maximum wind drag value
      CHARACTER(len=10) :: DragLawString ! controls wind drag formulation

      ! used for vis/diagnostics on nws19 and nws20
      LOGICAL :: writeFullCircleRmaxes = .false.
      LOGICAL :: writeRadialVandP = .false.
      LOGICAL :: writeSpatialVandP = .false.
      LOGICAL :: writeSpatialuvp = .false.
      LOGICAL :: writeNodaluvp = .false.
      LOGICAL :: moving_grid = .false.
      
      !switch to turn on or off geostrophic balance in GAHM
      ! on (default): Coriolis term included, phifactors will be calculated before being used 
      ! off:parameter is set to 'true', phifactors will be set to constant 1
      LOGICAL :: geostrophic_switch = .true.
      INTEGER :: geofactor = 1  !turn on or off gostrophic balance
      
      ! jgf50.38.05: used to keep track of HWind files and data
      type hWindData_t
         ! input data from fort.22
         character(1024)   :: file_name    ! full path to HWind file
         real(sz)          :: cycleTime    ! time (sec) relative to ADCIRC cold start or REFTIM
         real(sz)          :: Pc           ! storm central pressure (mb) or -1 when using p-w rel
         real(sz)          :: rampVal      ! ramping multiplier for this dataset
         ! input data from hwind files themselves
         real(sz)          :: dh           ! grid spacing (m) (assumed uniform)
         integer           :: nh           ! number of grid points (assumed same x and y dirs)
         real(sz)          :: cLon         ! longitude of center of storm (deg)
         real(sz)          :: cLat         ! latitude of center of storm (deg)
         real(sz), allocatable :: u(:,:)   ! east wind velocity (m/s)
         real(sz), allocatable :: v(:,:)   ! north wind velocity (m/s)
         ! derived data
         integer           :: indexLimit   ! data extent in index space (-indexLimit:indexLimit)
         real(sz)          :: vmax         ! maximum wind speed in the dataset (m/s)
         real(sz)          :: rmax         ! distance from storm center to vmax (nm)
         real(sz)          :: rmaxAngle    ! +CCW from east (deg) angle to rmax
         logical           :: loaded       ! .true. if we have loaded the data from file
      end type hWindData_t
      type(hWindData_t), allocatable, target :: hWindFiles(:) ! array of info about HWind files
      ! jgf50.38.05: Added the following for processing HWind data
      real(sz), allocatable :: wvnx_work(:) ! working array for wind u-velocity
      real(sz), allocatable :: wvny_work(:) ! working array for wind v-velocity
      real(8), allocatable :: x_mercator(:) ! mesh lons in mercator, origin at storm center
      real(8), allocatable :: y_mercator(:) ! mesh lats in mercator, origin at storm center
      integer, allocatable :: x_indices(:)  ! i index of each mesh node in hwind grid
      integer, allocatable :: y_indices(:)  ! j index of each mesh node in hwind grid
      real(sz), allocatable :: w1(:)        ! interp. weight for q11
      real(sz), allocatable :: w2(:)        ! interp. weight for q21
      real(sz), allocatable :: w3(:)        ! interp. weight for q22
      real(sz), allocatable :: w4(:)        ! interp. weight for q12
      real(sz), allocatable :: windSpeeds(:) ! u,v magnitude (m/s); used to find Vmax and Rmax
      logical, allocatable  :: inside(:)    ! true if a node is in the hwind bounding box
      integer               :: numFiles     ! number of hwind files listed in fort.22
      character(1024)       :: metComment   ! comment line at top of fort.22
      real(sz)              :: hWindMultiplier ! all hwind u,v multiplied by this
      character(32)         :: pressureWindRelationship ! dvorak, etc
      integer               :: currentCycle ! number of most recent
                                            ! applicable hwind file
      real(sz)              :: angleNow     ! +CCW from east (degrees) from storm center rmaxNow
                                            ! (deg) used by aswip to analyze radial vel

      ! tcm 51.06.02: used to keep track of GFDL Met files and data
      type gfdl_Data_t
         ! input data from fort.22
         character(1024)   :: file_name    ! full path to GFDL met files
         real(sz)          :: cycleTime    ! time (sec) relative to ADCIRC cold start or REFTIM
         real(sz)          :: rampVal          ! ramping multiplier for this dataset
         real(sz)          :: max_extrap_dist  ! maximum distance for extrapolating values (meters)
         integer           :: numFiles         !total number of GFDL met Files to be read
         ! input data from gfdl met files themselves
         integer           :: nh           ! number of grid points
         ! derived data
         logical           :: loaded       ! .true. if we have loaded the data from file
      end type gfdl_Data_t
      type(gfdl_Data_t), allocatable, target :: gfdl_Files(:) ! array of info about GFDL Met files
      real(sz) :: GFDL_max_extrap_dist     ! maximum distance for extrapolating values (meters)
      real(sz) :: GFDL_WindMultiplier      ! all GFDL wind u,v's multiplied by this
      !
      ! jgf52.51.21: Fixes for NWS10. Changed from NWLAT and NWLON to 
      ! latb and lonb and moved from read_input() to wind module and 
      ! gave them parameter status. 
      integer, parameter :: latb = 1536 !190 change to gfs t1534 grid
      integer, parameter :: lonb = 3072 !384 change to gfs t1534 grid
      ! numerical values at four corners of gaussian grid
      integer, allocatable ::  n00(:),n10(:),n11(:),n01(:)
      ! interpolation factors from mesh vertex to four corners of 
      ! gaussian grid
      real(sz), allocatable :: d00(:),d10(:),d11(:),d01(:)
      real(sz), allocatable :: colrab(:),dummy(:),gclat(:),gclon(:)
      ! met values on gaussian grid vertices
      real(sz), allocatable ::  ug(:),vg(:),pg(:)
      integer :: nws10DatasetCounter ! indicates which dataset to read

C------------------------end of data declarations-----------------------
      CONTAINS
C
C***********************************************************************
C                                                                      *
C   THE FOLLOWING SUBROUTINES READ IN AND IN SOME CASES INTERPOLATE    *
C   ONTO THE ADCIRC GRID WIND AND PRESSURE FIELDS IN VARIOUS INPUT     *
C   FORMATS.                                                           *
C                                                                      *
C   ALL WIND SPEEDS ARE CONVERTED TO M/S AND ALL PRESSURES TO M OF H20 *
C   BEFORE THEY ARE RETURNED.                                          *
C                                                                      *
C***********************************************************************

C-----------------------------------------------------------------------

C     ----------------------------------------------------------------
C       I N I T   W I N D   M O D U L E
C     ----------------------------------------------------------------
C
C     jgf50.32 Subroutine to initialize meteorological parameters to their
C     default values.
C
C     ----------------------------------------------------------------
      SUBROUTINE initWindModule()

      ! jgf50.32: The wind drag limit was originally set to 0.002d0
      ! in v46 ... the higher value of 0.0035d0 was used in ~v48 and
      ! is now the default.
      WindDragLimit = 0.0035d0

      ! jgf50.32: If DragLawString is not set in the fort.15 input file,
      ! it will remain 'default'; this will cause 'Garratt' to be used
      ! in FUNCTION WindDrag, and 'IceCube' to be used in
      ! FUNCTION WindIceDrag.
      
      DragLawString = 'default'
#ifdef POWELL
      ! jgf50.32 Provide backward compatibility with hardcoded POWELL
      ! support.
      DragLawString = 'Powell'
#endif

      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE initWindModule
C     ----------------------------------------------------------------


C     ----------------------------------------------------------------
C      F U N C T I O N   W I N D   D R A G
C     ----------------------------------------------------------------
C
C     jgf46.00 Function to calculate wind drag coefficient from the
C     windspeed.
C
C     ----------------------------------------------------------------
      REAL(SZ) FUNCTION WindDrag(WindSpeed, NodeNumber)
      !Casey 110518: Added variables for Mark Powell's sector-based wind drag.
      USE GLOBAL, ONLY: DEG2RAD
      USE MESH, ONLY : SFEA, SLAM
      REAL(SZ) WindSpeed
      ! jgf50.32: Even if present, NodeNumber is ignored by Garratt formulation
      INTEGER,OPTIONAL :: NodeNumber
C
      INTRINSIC :: ATAN2
      INTRINSIC :: MOD
      INTEGER   :: IS
      LOGICAL   :: FoundSector
      REAL(SZ)  :: Dir1
      REAL(SZ)  :: Dir2
      REAL(SZ)  :: Drag(3)
      REAL(SZ)  :: NodeDirection
      REAL(SZ)  :: NodeLat
      REAL(SZ)  :: NodeLon
      REAL(SZ)  :: StormDirection
      REAL(SZ)  :: Weight(3)
C
      SELECT CASE(TRIM(DragLawString))
C
      CASE("Garratt","GARRATT","garratt","default")
         WindDrag = 0.001d0*(0.75d0+0.067d0*WindSpeed)
         IF(WindDrag.gt.WindDragLimit) WindDrag=WindDragLimit
C
Casey 110518: Added Mark Powell's sector-based wind drag.
      CASE("Powell","POWELL","powell")
C.. Check whether we have previously found the eye.
         IF((EyeLatR(1).EQ.0.D0).OR.(EyeLonR(1).EQ.0.D0).OR.
     &      (EyeLatR(2).EQ.0.D0).OR.(EyeLonR(2).EQ.0.D0).OR.
     &      (EyeLatR(3).EQ.0.D0).OR.(EyeLonR(3).EQ.0.D0))THEN
            FoundEye = .FALSE.
         ENDIF
C.. If the eye has not been found then apply Garratt.
         WindDrag = 0.001D0*(0.75D0+0.067D0*WindSpeed)
         IF(.NOT.FoundEye)THEN
            FoundSector = .FALSE.
            Weight = 0.D0
            IF(WindDrag.GT.WindDragLimit) WindDrag = WindDragLimit
         ELSE
C.. We have found the eye!
C.. Compute the storm's direction.  We use the convention of northward
C.. being zero degrees, and then cycling clockwise.
            Dir1 = ATAN2(EyeLatR(2)-EyeLatR(1),EyeLonR(2)-EyeLonR(1))
            Dir2 = ATAN2(EyeLatR(3)-EyeLatR(2),EyeLonR(3)-EyeLonR(2))
            StormDirection = ATAN2(SIN(Dir1)+SIN(Dir2),COS(Dir1)+COS(Dir2))
            StormDirection = StormDirection/DEG2RAD
            StormDirection = 90.D0 - StormDirection
            IF(StormDirection.LT.0.D0)THEN
               StormDirection = StormDirection + 360.D0
            ENDIF
C.. Compute the direction of node relative to the eye.
            NodeLon = SLAM(NodeNumber)/DEG2RAD
            NodeLat = SFEA(NodeNumber)/DEG2RAD
            NodeDirection = ATAN2(NodeLat-EyeLatR(3),NodeLon-EyeLonR(3))
            NodeDirection = NodeDirection/DEG2RAD
            NodeDirection = 90.D0 - NodeDirection
            IF(NodeDirection.LT.0.D0)THEN
               NodeDirection = NodeDirection + 360.D0
            ENDIF
C.. Compute the weights for the sector in which the node is located.
C..     0- 40 : Transition from left (3) to right (1).
C..    40-130 : Right (1).
C..   130-170 : Transition from right (1) to rear (2).
C..   170-220 : Rear (2).
C..   220-260 : Transition from rear (2) to left (3).
C..   260-  0 : Left (3).
            FoundSector = .FALSE.
            DO IS=1,6
               IF(IS.EQ.1)THEN
                  Dir1 = MOD(StormDirection+  0.D0,360.D0)
                  Dir2 = MOD(StormDirection+ 40.D0,360.D0)
               ELSEIF(IS.EQ.2)THEN
                  Dir1 = MOD(StormDirection+ 40.D0,360.D0)
                  Dir2 = MOD(StormDirection+130.D0,360.D0)
               ELSEIF(IS.EQ.3)THEN
                  Dir1 = MOD(StormDirection+130.D0,360.D0)
                  Dir2 = MOD(StormDirection+170.D0,360.D0)
               ELSEIF(IS.EQ.4)THEN
                  Dir1 = MOD(StormDirection+170.D0,360.D0)
                  Dir2 = MOD(StormDirection+220.D0,360.D0)
               ELSEIF(IS.EQ.5)THEN
                  Dir1 = MOD(StormDirection+220.D0,360.D0)
                  Dir2 = MOD(StormDirection+260.D0,360.D0)
               ELSEIF(IS.EQ.6)THEN
                  Dir1 = MOD(StormDirection+260.D0,360.D0)
                  Dir2 = MOD(StormDirection+  0.D0,360.D0)
               ENDIF
               IF(Dir1.GT.Dir2)THEN
                  IF((Dir1.LE.NodeDirection).AND.(NodeDirection.LE.360.D0))THEN
                     Dir2 = Dir2 + 360.D0
                  ELSEIF((0.D0.LE.NodeDirection).AND.(NodeDirection.LE.Dir2))THEN
                     Dir1 = Dir1 - 360.D0
                  ENDIF
               ENDIF
               IF((Dir1.LE.NodeDirection).AND.(NodeDirection.LE.Dir2))THEN
                  IF(MOD(IS,2).EQ.0)THEN
                     FoundSector = .TRUE.
                     Weight = 0.D0
                     Weight(IS/2) = 1.D0
                  ELSE
                     FoundSector = .TRUE.
                     Weight = 0.D0
                     IF(IS.EQ.1)THEN
                        IW1 = 3
                        IW2 = 1
                     ELSEIF(IS.EQ.3)THEN
                        IW1 = 1
                        IW2 = 2
                     ELSEIF(IS.EQ.5)THEN
                        IW1 = 2
                        IW2 = 3
                     ENDIF
                     Weight(IW1) = 1.D0 + (0.D0-1.D0)/(Dir2-Dir1)*(NodeDirection-Dir1)
                     Weight(IW2) = 0.D0 + (1.D0-0.D0)/(Dir2-Dir1)*(NodeDirection-Dir1)
                  ENDIF
               ENDIF
            ENDDO
C.. Apply the wind drag limit in case of emergency.
            IF(.NOT.FoundSector)THEN
               Weight = 0.D0
               IF(WindDrag.GT.WindDragLimit) WindDrag = WindDragLimit
            ELSE
C.. Determine the wind drag for sector 1 (right).
               Drag(1) = WindDrag
               IF(WindSpeed.LE.35.0)THEN
                  IF(Drag(1).GT.0.0020)THEN
                     Drag(1) = 0.0020
                  ENDIF
               ELSEIF(WindSpeed.LE.45.0)THEN
                  Drag(1) = 0.0020 + (0.0030-0.0020)/(45.0-35.0)*(WindSpeed-35.0)
               ELSE
                  Drag(1) = 0.0030
               ENDIF
C.. Determine the wind drag for sector 2 (rear).
               Drag(2) = WindDrag
               IF(WindSpeed.LE.35.0)THEN
                  IF(Drag(2).GT.0.0020)THEN
                     Drag(2) = 0.0020
                  ENDIF
               ELSEIF(WindSpeed.LE.45.0)THEN
                  Drag(2) = 0.0020 + (0.0010-0.0020)/(45.0-35.0)*(WindSpeed-35.0)
               ELSE
                  Drag(2) = 0.0010
               ENDIF
C.. Determine the wind drag for sector 3 (left).
               Drag(3) = WindDrag
               IF(Drag(3).GT.0.0018)THEN
                  IF(WindSpeed.LE.25.0)THEN
                     Drag(3) = 0.0018
                  ELSEIF(WindSpeed.LE.30.0)THEN
                     Drag(3) = 0.0018 + (0.0045-0.0018)/(30.0-25.0)*(WindSpeed-25.0)
                  ELSEIF(WindSpeed.LE.45.0)THEN
                     Drag(3) = 0.0045 + (0.0010-0.0045)/(45.0-30.0)*(WindSpeed-30.0)
                  ELSE
                     Drag(3) = 0.0010
                  ENDIF
               ENDIF
C.. Apply a weighted average.
               WindDrag = Weight(1)*Drag(1) + Weight(2)*Drag(2) + Weight(3)*Drag(3)
            ENDIF
         ENDIF

      CASE DEFAULT
         WRITE(16,*) 'ERROR: Wind drag law not recognized:'
         WRITE(16,'(A10)') DragLawString
         WRITE(16,*) 'Execution will now be terminated.'
#ifdef CMPI
         call msg_fini()
#endif
         STOP
C
      END SELECT
C
      RETURN
C     ----------------------------------------------------------------
      END FUNCTION WindDrag
C     ----------------------------------------------------------------


C     ----------------------------------------------------------------
C      F U N C T I O N   W I N D  I C E  D R A G
C     ----------------------------------------------------------------
C
C     tcm v49.64.01 -- Added a function to modify wind drag coefficents
C          based on ice effects from a percentage of ice coverage
C          value (ice concentration limits 0.0 to 100.0)
C          If the percentage of ice is less than 1.0% then no
C          modifications are done to the incoming wind drag
C          coefficient values.
C
C     ----------------------------------------------------------------
      REAL(SZ) FUNCTION WindIceDrag(WindDrag,PercentIce)
      REAL(SZ) PercentIce
      REAL(SZ) pic,WindDrag,IceDrag
C
      SELECT CASE(TRIM(DragLawString))
C
      CASE("RaysIce")
C........If there is virtually no ice, then
C........just use the non-ice drag values
         if (PercentIce .lt. 1.d0) then
            WindIceDrag = WindDrag
         else
            pic = PercentIce*0.01d0
            IceDrag = (0.125d0 + 0.5d0*pic*(1.0d0-pic))*0.01d0
            WindIceDrag = max(IceDrag,WindDrag)
         endif

      CASE("IceCube","default")
C........This formula is a cubic function such that
C....... When pic = 0 the drag coefficient is equal to
C....... the minimum value for Garratt (.000075), when pic = 50
C....... it has a value of 0.0025 and a zero gradient, and
C....... when pic = 100 it has a value of 0.00125.
C....... This represents a smoother transition between Garratt
C....... and Ray's Ice for low ice values and low winds

         IF (PercentIce.lt.0.d0) then
            WindIceDrag = WindDrag
         ELSE
            pic = PercentIce*0.01d0
            IceDrag = (0.075d0 + 0.75d0*pic - 0.9d0*pic*pic
     &              + 0.2d0*pic*pic*pic)*0.01d0
            WindIceDrag = max(IceDrag,WindDrag)
         ENDIF

      CASE DEFAULT
         WRITE(16,*) 'ERROR: Ice drag law not recognized:'
         WRITE(16,'(A10)') DragLawString
         WRITE(16,*) 'Execution will now be terminated.'
#ifdef CMPI
         call msg_fini()
#endif
         STOP
C
      END SELECT
      IF(WindIceDrag.gt.WindDragLimit) WindIceDrag=WindDragLimit
C
      RETURN
C     ----------------------------------------------------------------
      END FUNCTION WindIceDrag
C     ----------------------------------------------------------------

C     ----------------------------------------------------------------
C     S U B R O U T I N E   G E T   B L E N D   F A C T O R
C     ----------------------------------------------------------------
C     jgf49.1001 Subroutine to get the blending factor for combining
C     NWS19 and NWS12 (from NAM) met data.
C     ----------------------------------------------------------------
      SUBROUTINE getBlendFactor(i, SLAM, SFEA, bf)
      IMPLICIT NONE
      INTEGER, intent(in) :: i  ! node number
      REAL(SZ), intent(in) :: SLAM(:)  ! nodal longitude, degrees
      REAL(SZ), intent(in) :: SFEA(:)  ! nodal latitude, degrees
      REAL(SZ), intent(out):: bf       ! blend factor, 0 to 1
      REAL(SZ) XCOOR ! longitude of node in degrees
      REAL(SZ) YCOOR ! latitude of node in degrees
C     distance from center of vortex to mesh node:
      REAL(SZ) DX    ! longitude distance in radians
      REAL(SZ) DY    ! latitude distance in radians
      REAL(SZ) rad   ! point-to-point distance over a sphere, in meters
      REAL(SZ) vortexOnly ! radial distance to edge of pure vortex, meters
      REAL(SZ) backgroundOnly ! radial dist to undisturbed background, meters
C
      XCOOR=SLAM(I)*RAD2DEG
      YCOOR=SFEA(I)*RAD2DEG
C
      DX=(XCOOR-vortexLon)*DEG2RAD
      DY=(YCOOR-vortexLat)*DEG2RAD
C
C     compute the distance from the center of the storm to the mesh
C     node along a sphere
      rad = sphericalDistance(DX, DY, vortexLat, YCOOR)
C
C     compute the radial distances to the various regions, in meters
      vortexOnly = pureVortex * vortexRMW
     &           * 1.852000003180799d0 * 1000.0d0
      backgroundOnly = pureBackground * vortexRMW
     &           * 1.852000003180799d0 * 1000.0d0
C     determine which interpolation region this mesh node falls in
      IF ( rad.le.vortexOnly ) bf = 1.0d0
      IF ( rad.ge.backgroundOnly ) bf = 0.0d0
C     interpolate linearly based on radial distance in the blended region
      IF ( (rad.gt.vortexOnly).and.(rad.lt.backgroundOnly) ) THEN
         bf = (rad-vortexOnly)/(backgroundOnly-vortexOnly)
Cjgfdebug         WRITE(16,*) "bf is ", bf
         IF ( (bf.gt.1.0d0).or.(bf.lt.0.0d0) ) THEN
            bf = 0.0d0
         ENDIF
      ENDIF
Cjgfdebug      WRITE(16,*) vortexLat, vortexLon, vortexRMW, rad, bf
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE getBlendFactor
C     ----------------------------------------------------------------

C***********************************************************************
C                                                                      *
C   Convert time from year,month,day,hour,min,sec into seconds since   *
C   the beginning of the year.                                         *
C                                                                      *
C***********************************************************************

      SUBROUTINE TIMECONV(IYR,IMO,IDAY,IHR,IMIN,SEC,TIMESEC,
     &     MyProc, NScreen, ScreenUnit)
      IMPLICIT NONE
      INTEGER IYR,IMO,IDAY,IHR,IMIN,ILEAP
      INTEGER MyProc, NScreen, ScreenUnit
      REAL*8 TIMESEC,SEC
C
      TIMESEC = (IDAY-1)*86400 + IHR*3600 + IMIN*60 + SEC
      IF(IMO.GE.2)  TIMESEC = TIMESEC + 31*86400
      ILEAP = (IYR/4)*4
      IF((ILEAP.EQ.IYR).AND.(IMO.GE.3)) TIMESEC = TIMESEC + 29*86400
      IF((ILEAP.NE.IYR).AND.(IMO.GE.3)) TIMESEC = TIMESEC + 28*86400
      IF(IMO.GE.4)  TIMESEC = TIMESEC + 31*86400
      IF(IMO.GE.5)  TIMESEC = TIMESEC + 30*86400
      IF(IMO.GE.6)  TIMESEC = TIMESEC + 31*86400
      IF(IMO.GE.7)  TIMESEC = TIMESEC + 30*86400
      IF(IMO.GE.8)  TIMESEC = TIMESEC + 31*86400
      IF(IMO.GE.9)  TIMESEC = TIMESEC + 31*86400
      IF(IMO.GE.10) TIMESEC = TIMESEC + 30*86400
      IF(IMO.GE.11) TIMESEC = TIMESEC + 31*86400
      IF(IMO.EQ.12) TIMESEC = TIMESEC + 30*86400
      IF(IMO.GT.12) THEN
         IF (NScreen.ne.0.and.MyProc.eq.0) THEN
            WRITE(ScreenUnit,*)
     &           'FATAL ERROR IN SUBROUTINE TIMECONV - MONTH > 12 '
         ENDIF
         WRITE(16,*) 'FATAL ERROR IN SUBROUTINE TIMECONV - MONTH > 12 '
#ifdef CMPI
         call msg_fini()
#endif
         STOP
      ENDIF
      RETURN
      END SUBROUTINE TIMECONV

C***********************************************************************
C                                                                      *
C   READ IN AND INTERPOLATE ONTO THE ADCIRC GRID WIND FIELDS FROM U.S. *
C   NAVY FLEET NUMERIC WIND FILES.                                     *
C                                                                      *
C   NOTE: The ADCIRC grid information consists only of the Lon and Lat *
C   of the nodes.  THE LONS AND LATS MUST BE IN RADIANS!               *
C                                                                      *
C                                                                      *
C   NWLAT = MAXIMUM NUMBER OF LATITUDES IN FLEET NUMERIC WIND FILE     *
C            SET = 1 IF FLEET NUMERIC WIND FILE NOT IN USE             *
C   NWLON = MAXIMUM NUMBER OF LONGITUDES IN FLEET NUMERIC WIND FILE    *
C            SET = 1 IF FLEET NUMERIC WIND FILE NOT IN USE             *
C                                                                      *
C                        R.L. 4/17/96                                  *
C                                                                      *
C   R.L. 4/2/01  changed MNWLAT,MNWLON in ALLOCATE statement to        *
C                NWLAT,NWLON                                           *
C***********************************************************************

      SUBROUTINE NWS3GET(X,Y,SLAM,SFEA,WVNX,WVNY,IWTIME,IWYR,WTIMED,NP,
     &     NWLON,NWLAT,WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,
     &     NScreen,ScreenUnit)
      USE SIZES
      IMPLICIT NONE
      INTEGER, SAVE :: FIRSTCALL = 0
      INTEGER IWTIME,IWYR,IWMO,IWDAY,IWHR,NP,NWLON,NWLAT,ICS,I,J
      INTEGER NScreen,ScreenUnit
      REAL*8 WTIMED
      REAL*8 X(*),Y(*),SLAM(*),SFEA(*),XCOOR,YCOOR
      INTEGER  LATIND1,LATIND2,LONIND1,LONIND2
      REAL(SZ) WLATMAX,WLONMIN,WLATINC,WLONINC,WSPEED,WDIR
      REAL(SZ) WLATM,WLONM,XWRATIO,YWRATIO
      REAL(SZ),ALLOCATABLE,SAVE :: WVXFN(:,:),WVYFN(:,:),PRN(:,:)
      REAL(SZ) WVNX(*),WVNY(*)
C
      IF (FIRSTCALL.EQ.0) THEN
         FIRSTCALL = 1
         ALLOCATE ( WVXFN(NWLAT,NWLON),WVYFN(NWLAT,NWLON),
     &             PRN(NWLAT,NWLON) )
      ENDIF
C
      READ(22,*) IWTIME
      IWYR = IWTIME/1000000
      IWMO = IWTIME/10000 - IWYR*100
      IWDAY = IWTIME/100 - IWYR*10000 - IWMO*100
      IWHR = IWTIME - IWYR*1000000 - IWMO*10000 - IWDAY*100
      CALL TIMECONV(IWYR,IWMO,IWDAY,IWHR,0,0.0D0,WTIMED,
     &     MyProc,NScreen,ScreenUnit)
C
      DO I=1,NWLAT
         READ(22,*) (WVXFN(I,J),J=1,NWLON)
      END DO
      DO I=1,NWLAT
         READ(22,*) (WVYFN(I,J),J=1,NWLON)
      END DO
C
      DO I=1,NWLAT              !CONVERT TO X AND Y COMPONENTS
         DO J=1,NWLON
            WSPEED=WVXFN(I,J)
            WDIR=WVYFN(I,J)*DEG2RAD
            WVXFN(I,J)=-WSPEED*SIN(WDIR)
            WVYFN(I,J)=-WSPEED*COS(WDIR)
         END DO
      END DO

      DO I=1,NP                 !INTERPOLATE TO ADCIRC GRID
         IF(ICS.EQ.2) THEN
            YCOOR=SFEA(I)*RAD2DEG
            XCOOR=SLAM(I)*RAD2DEG
         ENDIF
         IF(ICS.EQ.1) THEN
            YCOOR=Y(I)
            XCOOR=X(I)
         ENDIF
         LATIND2=(WLATMAX-YCOOR)/WLATINC + 1
         IF(LATIND2.EQ.NWLAT) LATIND2=LATIND2-1
         LATIND1=LATIND2 + 1
         LONIND1=(XCOOR-WLONMIN)/WLONINC + 1
         IF(LONIND1.EQ.NWLON) LONIND1=LONIND1-1
         LONIND2=LONIND1+1
         WLONM = WLONMIN + (LONIND1-1)*WLONINC
         WLATM = WLATMAX - (LATIND1-1)*WLATINC
         XWRATIO=(XCOOR-WLONM)/WLONINC
         YWRATIO=(YCOOR-WLATM)/WLATINC
C
         WVNX(I) = WVXFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVXFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVXFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVXFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         WVNY(I) = WVYFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVYFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVYFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVYFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
      END DO
C
      RETURN
      END SUBROUTINE NWS3GET



C***********************************************************************
C                                                                      *
C   Read onto the ADCIRC grid wind fields from the PBL-JAG model       *
C                                                                      *
C   Output from this subroutine is U,V (M/S) and P (M H20) on the      *
C   ADCIRC grid.                                                       *
C                                                                      *
C   The background pressure is assumed to be 1013 Mbars                *
C                                                                      *
C                           R.L.11/06/96                               *
C   R.L.09/04/00 added RHOWAT0 to call                                 *
C   R.L. 4/2/01  changed MNP dimensions to *                           *
C   R.L. 3/15/03 accounted for PRN=0                                   *   !RAL0315+ OK
C***********************************************************************

      SUBROUTINE NWS4GET(WVNX,WVNY,PRN,NP,RHOWAT0,G)
      USE SIZES
      IMPLICIT NONE
      INTEGER   NP,I,NHG
      REAL(SZ)  WVNX(*),WVNY(*),PRN(*)
      REAL(SZ)  RHOWAT0,RHOWATG,G
      CHARACTER*80 PBLJAGF
C
      RHOWATG=RHOWAT0*G
      DO I=1,NP
        WVNX(I)=0.d0
        WVNY(I)=0.d0
        PRN(I)=101300.d0/RHOWATG
      END DO
 170  READ(22,'(A80)') PBLJAGF
      IF(PBLJAGF(2:2).EQ.'#') GOTO 170
 171  READ(PBLJAGF,'(I8,5E13.5)') NHG,WVNX(NHG),WVNY(NHG),
     &                            PRN(NHG)
C
C     jgf46.02 From now on, wind files must contain data that are
C     appropriate for the time increment listed within the
C     file. Therefore, the following two lines were commented out.
C     WVNX(NHG)=WVNX(NHG)*1.04d0*0.5144d0 !CONVERT 30-MIN WINDS IN
C     WVNY(NHG)=WVNY(NHG)*1.04d0*0.5144d0 !KNOTS TO 10-MIN WIND IN M/S
C     jgf46.02 Added the following two lines.
      WVNX(NHG)=WVNX(NHG)*0.5144d0 !CONVERT KNOTS TO  M/S
      WVNY(NHG)=WVNY(NHG)*0.5144d0
      PRN(NHG)=100.d0*PRN(NHG)/RHOWATG !CONVERT MILLIBARS TO M OF WATER
      IF(PRN(NHG).EQ.0.) PRN(NHG)=101300.d0/RHOWATG                       !RAL0315+ OK
      READ(22,'(A80)') PBLJAGF
      IF(PBLJAGF(2:2).NE.'#') GOTO 171
      RETURN
      END SUBROUTINE


C***********************************************************************
C                                                                      *
C   Read in and interpolate onto the ADCIRC grid wind and pressure     *
C   fields from a meteorological file on a rectangular grid (either in *
C   Longitude, Latitude or Cartesian coordinates, consistent with the  *
C   ADCIRC grid coordinates).  If the ADCIRC grid is in Lon,Lat these  *
C   MUST BE IN RADIANS!                                                *
C                                                                      *
C   It is assumed that the meteorological grid is set up so that y     *
C   (e.g., latitude) varies from north (k=1) to south (k=NWLAT) and x  *
C   (e.g., longitude) varies from west (j=1) to east (j=NWLON).        *
C                                                                      *
C   The spatial extents of the meteorological grid must be consistent  *
C   with the ADCIRC model domain.  For example, if ADCIRC uses negative*
C   longitude values to indicate locations W of the Greenwich meridian,*
C   the meteorological file must be similarly organized.  Any grid that*
C   crosses the Greenwich Meridian should be organized so that the seam*
C   occurs @ 180 deg longitude. Therefore, the meteorological and      *
C   ADCIRC grids should use negative longitudes W of the Greenwich     *
C   Meridian and positive longitudes to the E.                         *
C                                                                      *
C                                                                      *
C   NOTE:  It is assumed that the met file data is oriented so that    *
C          the outer loop is on latitude and the inner loop is on      *
C          longitude.  For example:                                    *
C          line 1             lat 1,     lon 1                         *
C          line 2             lat 1,     lon 2                         *
C            .                                                         *
C          line nwlon         lat 1,     lon nwlon                     *
C          line nwlon+1       lat 2,     lon 1                         *
C          line nwlon+2       lat 2,     lon 2                         *
C            .                                                         *
C          line 2*nwlon       lat 2,     lon nwlon                     *
C          line 2*nwlon+1     lat 3,     lon 1                         *
C          line 2*nwlon+2     lat 3,     lon 2                         *
C            .                                                         *
C          line nwlon*nwlat   lat nwlat, lon nwlon                     *
C                                                                      *
C   NOTE:  It is assumed that he met file data is oriented so that     *
C          latitude varies from the northern most value (lat 1) to the *
C          southern most value (lat nwlat) and longitude varies in an  *
C          easterly direction (e.g. from 0 to 360 where positive       *
C          longitudes are angles measured easterly of the GM.          *
C                                                                      *
C   NOTE:  For the global AVN grid running from 0.5 - 359.5 deg lon    *
C          and 90 - -90 deg lat in 1 degree increments, NWLAT=181 and  *
C          NWLON=360 yielding a total number of entries in the file    *
C          of 65160.                                                   *
C                                                                      *
C   NOTE:  It is assumed that wind velocity is in EAST,NORTH components*
C          in M/2 and pressure is in N/M^2                             *
C                                                                      *
C   NOTE:  WLATMAX,WLONMIN,WLATINC,WLONINC should be in deg.           *
C                                                                      *
C   NOTE:  This should wrap if XCOORD > WLONMIN+NWLON*WLONINC  or      *
C          XCOORD < WLONMIN                                            *
C                                                                      *
C                                                                      *
C   MNWLAT = MAXIMUM NUMBER OF LATITUDES IN WIND FILE                  *
C            SET = 1 IF FLEET NUMERIC WIND FILE NOT IN USE             *
C   MNWLON = MAXIMUM NUMBER OF LONGITUDES IN WIND FILE                 *
C            SET = 1 IF FLEET NUMERIC WIND FILE NOT IN USE             *
C                                                                      *
C                           R.L. 4/13/99                               *
C                           R.L.09/04/00 added RHOWAT0 to call         *
C   R.L.09/04/00 added RHOWAT0 to call                                 *
C   R.L. 4/2/01  changed MNWLAT,MNWLON in ALLOCATE statement to        *
C                NWLAT,NWLON                                           *
C   R.L. 8/10/05 eliminated adding 360 to negative longitudes to match *
C                AVN model grid setup.  User is now required to provide*
C                met and ADCIRC grid that correspond in space without  *
C                adjusting the longitude.  Also the input variable     *
C                order has been changed to U,V,P to be consistent with *
C                other NWS input formats.                              *
C***********************************************************************

      SUBROUTINE NWS6GET(X,Y,SLAM,SFEA,WVNX,WVNY,PRESS,NP,NWLON,NWLAT,
     &     WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
      USE SIZES

      IMPLICIT NONE
      INTEGER, SAVE :: FIRSTCALL = 0
      INTEGER NP,NWLON,NWLAT,I,J,ICS
      REAL(SZ) RHOWAT0,RHOWATG,G
      INTEGER  LATIND1,LATIND2,LONIND1,LONIND2
      REAL(SZ) WLATMAX,WLONMIN,WLATINC,WLONINC,XWRATIO,YWRATIO
      REAL(SZ) WLATM,WLONM
      REAL*8 X(*),Y(*),SLAM(*),SFEA(*),XCOOR,YCOOR
      REAL(SZ) WVNX(*),WVNY(*),PRESS(*)
      REAL(SZ),SAVE,ALLOCATABLE :: WVXFN(:,:),WVYFN(:,:),PRN(:,:)
C
      IF (FIRSTCALL.EQ.0) THEN
         FIRSTCALL = 1
         ALLOCATE ( WVXFN(NWLAT,NWLON),WVYFN(NWLAT,NWLON),
     &        PRN(NWLAT,NWLON) )
      ENDIF
C
      RHOWATG=RHOWAT0*G
      DO I=1,NWLAT
         DO J=1,NWLON
            READ(22,*) WVXFN(I,J),WVYFN(I,J),PRN(I,J)
         END DO
      END DO

      DO I=1,NP                 !INTERPOLATE TO ADCIRC GRID
         IF(ICS.EQ.2) THEN
            YCOOR=SFEA(I)*RAD2DEG
            XCOOR=SLAM(I)*RAD2DEG
         ENDIF
         IF(ICS.EQ.1) THEN
            YCOOR=Y(I)
            XCOOR=X(I)
         ENDIF
         LATIND2=(WLATMAX-YCOOR)/WLATINC + 1
         IF(LATIND2.EQ.NWLAT) LATIND2=LATIND2-1
         LATIND1=LATIND2 + 1
         LONIND1=(XCOOR-WLONMIN)/WLONINC + 1
         LONIND2=LONIND1 + 1
C
         WLONM = WLONMIN + (LONIND1-1)*WLONINC
         WLATM = WLATMAX - (LATIND1-1)*WLATINC
         XWRATIO=(XCOOR-WLONM)/WLONINC
         YWRATIO=(YCOOR-WLATM)/WLATINC
C
         IF(LONIND1.EQ.0) LONIND1=NWLON
         IF(LONIND1.EQ.NWLON) LONIND2=1
C
         WVNX(I) = WVXFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVXFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVXFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVXFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         WVNY(I) = WVYFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVYFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVYFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVYFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         PRESS(I) = PRN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + PRN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + PRN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + PRN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         PRESS(I) = PRESS(I)/RHOWATG

      END DO
C
      RETURN
      END SUBROUTINE NWS6GET

C     ----------------------------------------------------------------
C      S U B R O U T I N E     N W S 7 G E T
C     ----------------------------------------------------------------
C
C     jgf46.01 Subroutine to get a set of surface wind stresses and
C     barometric pressure on a rectangular grid (either in Longitude,
C     Latitude or Cartesian coordinates, consistent with the ADCIRC grid
C     coordinates) from unit 22 and interpolate them in space onto the
C     ADCIRC grid. If the ADCIRC grid is in Lon, Lat these must be in
C     radians.
C
C     It is assumed that the meteorological grid is set up so that y
C     (e.g., latitude) varies from north (k=1) to south (k=NWLAT) and x
C     (e.g., longitude) varies from west (j=1) to east (j=NWLON).
C
C     The spatial extents of the meteorological grid must be consistent
C     with the ADCIRC model domain.  For example, if ADCIRC uses
C     negative longitude values to indicate locations W of the Greenwich
C     meridian, the meteorological file must be similarly organized.
C     Any grid that crosses the Greenwich Meridian should be organized
C     so that the seam occurs @ 180 deg longitude. Therefore, the
C     meteorological and ADCIRC grids should use negative longitudes W
C     of the Greenwich Meridian and positive longitudes to the E.
C
C   NOTE:  It is assumed that the met file data is oriented so that
C          the outer loop is on latitude and the inner loop is on
C          longitude.  For example:
C          line 1             lat 1,     lon 1
C          line 2             lat 1,     lon 2
C            .
C          line nwlon         lat 1,     lon nwlon
C          line nwlon+1       lat 2,     lon 1
C          line nwlon+2       lat 2,     lon 2
C            .
C          line 2*nwlon       lat 2,     lon nwlon
C          line 2*nwlon+1     lat 3,     lon 1
C          line 2*nwlon+2     lat 3,     lon 2
C            .
C          line nwlon*nwlat   lat nwlat, lon nwlon
C
C   NOTE:  It is assumed that he met file data is oriented so that
C          latitude varies from the northern most value (lat 1) to the
C          southern most value (lat nwlat) and longitude varies in an
C          easterly direction (e.g. from 0 to 360 where positive
C          longitudes are angles measured easterly of the GM.
C
C   NOTE:  It is assumed that wind stress is in EAST, NORTH components
C          in m/s and pressure is in N/m^2
C
C   NOTE:  WLATMAX,WLONMIN,WLATINC,WLONINC should be in deg.
C
C   NOTE:  This should wrap if XCOORD > WLONMIN+NWLON*WLONINC  or
C          XCOORD < WLONMIN
C
C     ----------------------------------------------------------------
      SUBROUTINE NWS7GET(X,Y,SLAM,SFEA,WVNX,WVNY,PRESS,NP,NWLON,NWLAT,
     &     WLATMAX,WLONMIN,WLATINC,WLONINC,ICS,RHOWAT0,G)
      IMPLICIT NONE
      INTEGER, intent(in) :: NP     ! number of nodes
      REAL(8), intent(in), dimension(NP) :: X
      REAL(8), intent(in), dimension(NP) :: Y
      REAL(8), intent(in), dimension(NP) :: SLAM
      REAL(8), intent(in), dimension(NP) :: SFEA
      REAL(SZ), intent(out), dimension(NP) :: WVNX
      REAL(SZ), intent(out), dimension(NP) :: WVNY
      REAL(SZ), intent(out), dimension(NP) :: PRESS
      INTEGER, intent(in) :: NWLON  ! max. # of longitudes in stress file
      INTEGER, intent(in) :: NWLAT  ! max. # of latitutes in stress file
      REAL(SZ), intent(in) :: WLATMAX
      REAL(SZ), intent(in) :: WLONMIN
      REAL(SZ), intent(in) :: WLATINC
      REAL(SZ), intent(in) :: WLONINC
      INTEGER, intent(in) :: ICS    ! coord.sys., 1=cartesian, 2=spherical
      REAL(SZ), intent(in):: RHOWAT0! ref. dens. of water
      REAL(SZ), intent(in):: G      ! gravitational constant

      INTEGER I                     ! node loop counter
      INTEGER K, J                  ! latitude, longitude loop counters
      REAL(SZ) RHOWATG             ! ref. dens. of water * grav. constant
      INTEGER  LATIND1,LATIND2,LONIND1,LONIND2
      REAL(SZ) WLATM,WLONM,XWRATIO,YWRATIO
      REAL(8) XCOOR,YCOOR
      REAL(SZ), SAVE, ALLOCATABLE :: WVXFN(:,:),WVYFN(:,:),PRN(:,:)
      LOGICAL, SAVE :: MemoryAllocated = .False.
C
      IF (.not.MemoryAllocated) THEN
         ALLOCATE ( WVXFN(NWLAT,NWLON),WVYFN(NWLAT,NWLON),
     &        PRN(NWLAT,NWLON) )
         MemoryAllocated = .True.
      ENDIF
C
      RHOWATG=RHOWAT0*G
      DO K=1,NWLAT
         DO J=1,NWLON
            READ(22,*) WVXFN(K,J),WVYFN(K,J),PRN(K,J)
         END DO
      END DO
C
      DO I=1,NP                 !INTERPOLATE TO ADCIRC GRID
         IF(ICS.EQ.2) THEN
            YCOOR=SFEA(I)*RAD2DEG
            XCOOR=SLAM(I)*RAD2DEG
         ENDIF
         IF(ICS.EQ.1) THEN
            YCOOR=Y(I)
            XCOOR=X(I)
         ENDIF
         LATIND2=(WLATMAX-YCOOR)/WLATINC + 1
         IF(LATIND2.EQ.NWLAT) LATIND2=LATIND2-1
         LATIND1=LATIND2 + 1
         LONIND1=(XCOOR-WLONMIN)/WLONINC + 1
         LONIND2=LONIND1 + 1
C
         WLONM = WLONMIN + (LONIND1-1)*WLONINC
         WLATM = WLATMAX - (LATIND1-1)*WLATINC
         XWRATIO=(XCOOR-WLONM)/WLONINC
         YWRATIO=(YCOOR-WLATM)/WLATINC
C
         IF(LONIND1.EQ.0) LONIND1=NWLON
         IF(LONIND1.EQ.NWLON) LONIND2=1
C
         WVNX(I) = WVXFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVXFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVXFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVXFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         WVNY(I) = WVYFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVYFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVYFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVYFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         PRESS(I) = PRN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + PRN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + PRN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + PRN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         PRESS(I) = PRESS(I)/RHOWATG
      END DO
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE NWS7GET
C     ----------------------------------------------------------------
   

      !----------------------------------------------------------------
      !         S U B R O U T I N E   N W S  1 0  I N I T
      !----------------------------------------------------------------
      ! jgf52.51.21: Decoupled the NWS10 initialization from the 
      ! rest of the subroutine to correct numerical issues associated
      ! with uninitialized variables.
      !----------------------------------------------------------------
      subroutine nws10Init(flon,flat)
      use global, only : ihot
      use mesh, only : np
      implicit none
      real(sz), intent(in) :: flon(:) ! adcirc mesh longitudes (rad)
      real(sz), intent(in) :: flat(:) ! adcirc mesh latitudes (rad)
      real(sz) :: gdlon
      integer :: j, jj ! loop counters
      !
      ! gaussian met values 
      allocate ( ug(latb*lonb),vg(latb*lonb),pg(latb*lonb) )
      ! four corner met values in grid cell associated with each adcirc mesh vertex 
      allocate ( n00(np),n10(np),n11(np),n01(np) )
      ! four interpolating factors associated with each adcirc mesh vertex
      allocate ( d00(np),d10(np),d11(np),d01(np) )
      ! 
      allocate ( colrab(latb),dummy(latb),gclat(latb),gclon(lonb) )
C
C...  Set up the Gaussian grid and determine the interpolating factors 
C     for the ADCIRC grid.
      call glats(latb/2,colrab,dummy,dummy,dummy)
      do j=1,latb/2 
         gclat(j)=colrab(j)
         jj=latb-j+1
         gclat(jj)=pi-colrab(j)
      enddo
      gdlon=twopi/lonb
      do j=1,lonb
         gclon(j)=gdlon*(j-1)
      end do
      call g2rini(flon,flat)
      !
      ! Set the counter value. According to the ADCIRC documentation,
      ! fort.200 is read upon hotstart, but not at coldstart. 
      if (ihot.eq.0) then
         nws10DatasetCounter = 1
      else
         nws10DatasetCounter = 0
      endif
         
      !----------------------------------------------------------------
      end subroutine nws10Init
      !----------------------------------------------------------------

C***********************************************************************
C  Subroutine to compute the factors to interpolate from a global      *
C  Gaussian Lat/Lon grid with T126 resolution (GRIB Grid type 126)     *
C  onto another grid.                                                  *
C                                                                      *
C  The new grid is a series of longitude and latitude points contained *
C  in the FLON and FLAT arrays with a total number of points NP        *
C                                                                      *
C       modified from the original G2RINI by R.L. 4/17/96              *
C    
C     jgf51.52.21: Removed duplicate variables; moved grid locations
C     and interpolating factors from subroutine arguments to module
C     variables; brought log messages up to date with current code.
C
C***********************************************************************
      subroutine g2rini(flon,flat)

      use mesh, only : np
      implicit none
      real(8), intent(in) :: flat(:),flon(:) ! adcirc mesh coordinates
      integer :: N,I,LON,LONP1,LAT,LATP1
      real(8) :: DLAT,DLON,FLONWORK,COLAT,DDLAT,XLAT,DFLAT,DFLAT1,
     &     DDLON,XLON,DFLON,DFLON1
C
C...Compute estimated DLAT, true DLON for Gaussian grid
      DLAT=PI/REAL(LATB-1,KIND=SZ)
      DLON=TWOPI/REAL(LONB,KIND=SZ)
      N=0

C...Loop through all the nodes in the grid to be interpolated onto and
C.....compute the interpolating factors.

      DO I=1,NP
C.....   Compute initial guess of which lon value FLON(I) is in the Gaussian file
C......  Check that this value is reasonable.
         flonwork=flon(i)
         if(flonwork.lt.0.) flonwork=flonwork+twopi
         lon=flonwork/dlon + 1
         lonp1=lon+1
         if(lon.eq.lonb) lonp1=1 !circle condition
         if((lon.lt.1).or.(lon.gt.lonb)) then
            call allMessage(ERROR,'Mesh node outside gaussian grid.')
            write(scratchMessage,
     &       '("i=",i0," lon=",f15.8," dlon=",f15.8," flon(i)=",f15.8)')
     &       i, lon, dlon, flon(i)
            call allMessage(ERROR, scratchMessage)
            call windTerminate()
         endif

C.....   Compute initial guess of which lat value FLAT(I) is in the Gaussian file
C....... Check that this value is reasonable.
         colat=hfpi-flat(i)
         lat=colat/dlat + 1
         if(lat.eq.latb) lat=lat-1
         latp1=lat+1
         if((lat.lt.1).or.(lat.gt.latb)) then
            call allMessage(ERROR,'Mesh node outside gaussian grid.')
            write(scratchMessage,
     &       '("i=",i0," lat=",f15.8," dlat=",f15.8," flat(i)=",f15.8)')
     &       i, lat, dlat, flat(i)
            call allMessage(ERROR, scratchMessage)
            call windTerminate()
         endif

5        CONTINUE
         IF((COLAT.GE.GCLAT(LAT)).AND.(COLAT.LE.GCLAT(LATP1))) GO TO 9
         IF(COLAT.LT.GCLAT(LAT)) THEN
            LATP1=LAT
            LAT=LAT-1
            IF(LAT.LE.0) THEN
               LAT=1
               LATP1=2
               GOTO 9
            ENDIF
            GOTO 5
         ENDIF
         IF(COLAT.GT.GCLAT(LATP1)) THEN
            LAT=LAT+1
            LATP1=LAT+1
            IF(LAT.GE.LATB ) THEN
               LAT=LATB-1
               LATP1=LATB
               GOTO 9
            ENDIF
            GOTO 5
         ENDIF

 9       CONTINUE
         DDLAT=GCLAT(LATP1)-GCLAT(LAT)
         XLAT=GCLAT(LAT)
         DFLAT1=(COLAT-XLAT)/DDLAT
         IF(LAT.EQ.1) DFLAT1=MAX(0.d0,DFLAT1) !MODIFY THIS FOR POLAR POINTS
         IF(LATP1.EQ.LATB) DFLAT1=MIN(1.d0,DFLAT1) !MODIFY THIS FOR POLAR POINTS
         DFLAT=1.d0-DFLAT1
         DDLON=DLON
         XLON=GCLON(LON)
         DFLON1=(FLONWORK-XLON)/DDLON
         DFLON=1.d0-DFLON1
         N=N+1
         D00(N)=DFLON*DFLAT
         D10(N)=DFLON1*DFLAT
         D11(N)=DFLON1*DFLAT1
         D01(N)=DFLON*DFLAT1
         N00(N)=LON+(LAT-1)*LONB
         N10(N)=LONP1+(LAT-1)*LONB
         N11(N)=LONP1+(LATP1-1)*LONB
         N01(N)=LON+(LATP1-1)*LONB

      END DO
      !----------------------------------------------------------------
      end subroutine g2rini
      !----------------------------------------------------------------


C***********************************************************************
C                                                                      *
C   Read in and interpolate onto the ADCIRC grid wind fields from U.S. *
C   National Weather Service AVN model SFLUX meteorological files.     *
C                                                                      *
C   The input files are in binary and have been created by the GRIB    *
C   unpacking program unpkgrb1.f to extract only the U 10M, V 10M, and *
C   surface P fields.    THE BINARY INPUT HAS BEEN ELIMINATED!!!!      *
C   The input files are in ASCII and contain surface P, U 10M and V 10M*
C   fields.                                                            *
C                                                                      *
C   The SFLUX files utilize a global Gaussian Lon/Lat grid which is    *
C   constructed in these subroutines.                                  *
C                                                                      *
C   NOTE: The ADCIRC grid information consists only of the Lon and Lat *
C   of the nodes.  THE LONS AND LATS MUST BE IN RADIANS!               *
C                                                                      *
C   Output from this subroutine is U,V (M/S) and P (M H20) on the      *
C   ADCIRC grid.                                                       *
C                                                                      *
C   MNWLAT = LATB = 190    FOR GAUSSIAN GRID                           *
C   MNWLON = LONB = 384    FOR GAUSSIAN GRID                           *
C                                                                      *
C                           R.L. 4/14/99                               *
C                           R.L.09/04/00 added RHOWAT0 to call         *
C   R.L. 4/2/01  changed MNWLAT,MNWLON in ALLOCATE statement to        *
C                LATB,LONB; eliminated MNWP as a dimension             *
C***********************************************************************
      subroutine nws10get(flon,flat,ull,vll,pll)
      use sizes, only : gblinputdir
      use global, only : wtiminc
      use mesh, only : np
      implicit none
      REAL(sz), intent(in) :: flat(:) ! mesh latitudes (radians)
      REAL(sz), intent(in) :: flon(:) ! mesh longitudes (radians)
      real(sz), intent(out) :: ull(:),vll(:),pll(:) ! mesh met values
      !
      real(sz) :: rhowatg ! density of water * gravitational acceleration
      real(sz) :: p1,p2,p3,p4 ! pressure at 4 corners of gaussian grid cell 
      real(sz) :: u1,u2,u3,u4 ! east velocity at 4 corners of gaussian grid cell
      real(sz) :: v1,v2,v3,v4 ! north velocity at 4 corners of gaussian grid cell
      character(len=8) :: fname1 ! name of GFS file to read
      integer :: iext        ! numerical extension to file name
      integer :: ios = 0     ! i/o status for opening the fort.xxx files
      integer :: i, n        ! loop counters
C
      rhowatg=rhowat0*g
      ! 
      ! Check to be sure that wtiminc is evenly divisible.
      if (modulo(int(wtiminc),3600).ne.0) then
         call allMessage(ERROR,
     &      'For NWS=10, WTIMINC must be evenly divisible by 3600.')
         call windTerminate()
      endif
      ! 
      ! Form the data file name.
      iext = 200 + nws10DatasetCounter * int(wtiminc/3600.d0) 
      write(fname1,'("fort.",i3)') iext    
      !
      ! jgf51.52.21: Fixed name of met file with full path for use in 
      ! parallel operation.  
      call openFileForRead(iext,trim(gblinputdir)//'/'//trim(fname1),ios)
      if (ios.ne.0) then
         call windTerminate()
      endif
C
C...  Read the ASCII data file
      do i=1,lonb*latb
         read(iext,*) pg(i),ug(i),vg(i)
      enddo
      close(iext)
C
C.....Go from the Gaussian grid to the ADCIRC grid
C.....Convert pressure from N/M^2 to M of H20
      do n=1,np
         ! grab the numerical values at the 4 corners of the grid 
         ! cell containing this mesh node
         p1=pg(n00(n))
         p2=pg(n10(n))
         p3=pg(n11(n))
         p4=pg(n01(n))
         u1=ug(n00(n))
         u2=ug(n10(n))
         u3=ug(n11(n))
         u4=ug(n01(n))
         v1=vg(n00(n))
         v2=vg(n10(n))
         v3=vg(n11(n))
         v4=vg(n01(n))
         ! spatially interpolate from the 4 values to the mesh vertex
         pll(n)=p1*d00(n)+p2*d10(n)+p3*d11(n)+p4*d01(n)
         ull(n)=u1*d00(n)+u2*d10(n)+u3*d11(n)+u4*d01(n)
         vll(n)=v1*d00(n)+v2*d10(n)+v3*d11(n)+v4*d01(n)
         pll(n)=pll(n)/rhowatg
      end do
C
      nws10DatasetCounter = nws10DatasetCounter + 1
      !----------------------------------------------------------------      
      end subroutine nws10get
      !----------------------------------------------------------------

C***********************************************************************
C  Subroutine to compute the latutudes in a Global Gaussian Lat/Lon    *
C  grid with T126 resolution (GRIB Grid type 126).                     *
C                                                                      *
C      modified from the original GLATS by R.L. 4/24/96                *
C      modified from T126 to T1534 by Y.F 9/1/2015       
C
C      yf51.52.25: Wiki link: https://en.wikipedia.org/wiki/Gaussian_grid
C***********************************************************************
      SUBROUTINE GLATS(LGGHAF,COLRAD,WGT,WGTCS,RCS2)
      USE SIZES
      IMPLICIT NONE
      REAL(SZ) COLRAD(*),WGT(*),WGTCS(*),RCS2(*)
      INTEGER LGGHAF,L2,K,K1,ITER
      REAL(SZ) SI,SCALEVAL,RL2,DRAD,RAD,P1,P2,EPS,PHI,X,W,SN,RC
      real(sz) :: dradz ! initial guess that represents a distance between gridpoints (radian)
      integer :: multiple ! yf51.52.25: reduce initial guess for DRAD for T1534 grid  
C
      EPS=1.d-6
C     EPS=1.d-12
C     PRINT 101
C     101  FORMAT ('0 I   COLAT   COLRAD     WGT', 12X, 'WGTCS',
CCCC  1 10X, 'ITER  RES')
C
      SI = 1.0d0
      L2=2*LGGHAF
      RL2=L2
      SCALEVAL = 2.0d0/(RL2*RL2)
      K1=L2-1
      !yf51.52.25 DRADZ = PI / 360.d0
      !yf51.52.25: In equation L2/(190+2), L2 (1536) is a number of 
      ! latitude from T1534 grid and 190 is ones from the original 
      ! T126 grid. +2 is polar points which didn't include in T126 grid.
      ! The reason why I used the equation and the integer "multiple" 
      ! are I want to use the original T126 latitudinal number (as the
      ! base number) for the future update of GFS.
      multiple = L2/(190+2)
      dradz = pi/(multiple*360.d0)
      RAD = 0.0
      DO 1000 K=1,LGGHAF
         ITER=0
         DRAD=DRADZ
 1       CALL POLY(L2,RAD,P2)
 2       P1 =P2
         ITER=ITER+1
         RAD=RAD+DRAD
         CALL POLY(L2,RAD,P2)
         IF(SIGN(SI,P1).EQ.SIGN(SI,P2)) GO TO 2
         IF(DRAD.LT.EPS)GO TO 3
         RAD=RAD-DRAD
         DRAD = DRAD * 0.25d0
         GO TO 1
 3       CONTINUE
         COLRAD(K)=RAD
         PHI = RAD * 180.d0 / PI
         CALL POLY(K1,RAD,P1)
         X = COS(RAD)
         W = SCALEVAL * (1.0d0 - X*X)/ (P1*P1)
         WGT(K) = W
         SN = SIN(RAD)
         W=W/(SN*SN)
         WGTCS(K) = W
         RC=1.d0/(SN*SN)
         RCS2(K) = RC
         CALL POLY(L2,RAD,P1)
C     PRINT 102,K,PHI,COLRAD(K),WGT(K),WGTCS(K),ITER,P1
C     102  FORMAT(1H ,I2,2X,F6.2,2X,F10.7,2X,E13.7,2X,E13.7,2X,I4,2X,D13.7)
 1000 CONTINUE
c     PRINT 100,LGGHAF
c     100  FORMAT(1H ,'SHALOM FROM 0.0 GLATS FOR ',I3)
      RETURN
      END SUBROUTINE GLATS


C***********************************************************************
C  Subroutine used by GLATS.                 
C     yf51.52.25: The subroutine POLY calculates the polynomial used in
C     Gaussian quadrature.
C***********************************************************************
      SUBROUTINE POLY(N,RAD,P)
      USE SIZES
      IMPLICIT NONE
      INTEGER N,I
      REAL(SZ) RAD,P,X,Y1,Y2,Y3,G
C
      X = COS(RAD)
      Y1 = 1.0d0
      Y2=X
      DO 1 I=2,N
         G=X*Y2
         Y3=G-Y1+G-(G-Y1)/REAL(I,KIND=SZ)
         Y1=Y2
         Y2=Y3
 1    CONTINUE
      P=Y3
      RETURN
      END SUBROUTINE POLY


C***********************************************************************
C                                                                      *
C   Read in and interpolate onto the ADCIRC grid wind fields from U.S. *
C   National Weather Service ETA-29 model that have been stripped down *
C   and given to us by NOAA.                                           *
C                                                                      *
C   The input files are in binary and have been created by NOAA and    *
C   contain only the U 10M, V 10M, (M/S) and surface P fields (mbars). *
C                                                                      *
C   The ETA-29 model uses an E grid and therefore the U and V          *
C   components are not oriented along lines of constant latitute and   *
C   longitude. These must be converted to be useful in ADCIRC.         *
C                                                                      *
C   NOTE: The ADCIRC grid information consists only of the Lon and Lat *
C   of the nodes.  THE LONS AND LATS MUST BE IN RADIANS!               *
C                                                                      *
C   Output from this subroutine is U,V (M/S) and P (M H20) on the      *
C   ADCIRC grid.                                                       *
C                                                                      *
C   MNWLAT = LATB = 271    FOR ETA-29 GRID                             *
C   MNWLON = LONB = 181    FOR ETA-29 GRID                             *
C                                                                      *
C                           R.L. 1/11/97                               *
C   R.L.09/04/00 added RHOWAT0 to call                                 *
C   R.L. 4/02/01  elminiated MNWP as a dimension                       *
C   R.L. 9/14/01  changed MNWLAT,MNWLON in ALLOCATE statement to       *
C                271,181                                               *
C***********************************************************************

      SUBROUTINE NWS11GET(NWSEGWI,IDSETFLG,FLON,FLAT,ULL,VLL,PLL,NP,
     &     RHOWAT0,G)
      USE SIZES
      IMPLICIT NONE
      INTEGER,SAVE  ::  ICALL = 0
      INTEGER NWSEGWI,IDSETFLG,NP,I,IEXT,IDIG1,IDIG2,IDIG3,KERR,N
      INTEGER IYEAR,IMONTH,IDAY,IHOUR
      REAL*8 RHOWATG100,FLONDEG,FLATDEG
      REAL(SZ) P1,P2,P3,U1,U2,U3,V1,V2,V3,UE29,VE29,CBETAU,SBETAU,G
      REAL(SZ) RHOWAT0
      REAL(SZ) ULL(*),VLL(*),PLL(*)
      REAL*8 FLAT(*),FLON(*)
C
      INTEGER,SAVE,ALLOCATABLE ::  N1(:),N2(:),N3(:)
      REAL(SZ),SAVE,ALLOCATABLE :: D1(:),D2(:),D3(:),BETAU(:)
      REAL(SZ),SAVE,ALLOCATABLE :: UE(:),VE(:),PE(:)
C
      CHARACTER*1 FNAME2(8)
      CHARACTER*8 FNAME1
      EQUIVALENCE (FNAME1,FNAME2)
      LOGICAL FOUND
C
      IF (ICALL.EQ.0) THEN
         ICALL = 1
         ALLOCATE ( N1(MNP),N2(MNP),N3(MNP) )
         ALLOCATE ( D1(MNP),D2(MNP),D3(MNP),BETAU(MNP) )
         ALLOCATE ( UE(181*271),VE(181*271),PE(181*271) )
      ENDIF
C
      RHOWATG100=RHOWAT0*G*100.d0

C...  The first time the subroutine is called, setup the interpolating factors
C...  between the Eta-29 grid and the ADCIRC grid.

      IF((NWSEGWI.EQ.0).AND.(IDSETFLG.EQ.0)) THEN
         if (myproc == 0) then
          WRITE(screenunit,*) 'Computing ETA29 met field interp factors'
         endif
         DO I=1,NP
            flondeg=rad2deg*flon(i)
            flatdeg=rad2deg*flat(i)
            CALL E29SEARCH(I,FLONDEG,FLATDEG,N1(I),N2(I),N3(I),
     &           D1(I),D2(I),D3(I),betau(i))
         END DO
         RETURN
      ENDIF

C...  Figure out the met data file name

      FNAME1='fort.   '
      IEXT=200 + NWSEGWI
      IDIG1=IEXT/100
      IDIG2=(IEXT-100*IDIG1)/10
      IDIG3=(IEXT-100*IDIG1-10*IDIG2)
      FNAME2(6)=CHAR(IDIG1+48)
      FNAME2(7)=CHAR(IDIG2+48)
      FNAME2(8)=CHAR(IDIG3+48)

C...  If appropriate, enter, locate and open the met data file

 1010 FORMAT(' File ',A8,' WAS NOT FOUND!  FATAL ERROR',/)
 1011 FORMAT(' File ',A8,' WAS FOUND!  Opening & Processing file',/)

      if (myproc == 0) WRITE(screenunit,*) '  '
      INQUIRE(FILE=FNAME1,EXIST=FOUND)
      IF(FOUND) GOTO 32
      if (myproc == 0) WRITE(screenunit,1010) FNAME1
      WRITE(16,1010) FNAME1
#ifdef CMPI
      call msg_fini()
#endif
      STOP
 32   if (myproc == 0) WRITE(screenunit,1011) FNAME1
      IF((NWSEGWI.EQ.0).OR.(IDSETFLG.EQ.1)) OPEN(IEXT,FILE=FNAME1,
     &status='old',access='sequential',form='unformatted',iostat=kerr)

C...  Read the met data file

      READ(IEXT,END=1100) IYEAR,IMONTH,IDAY,IHOUR
      READ(IEXT,END=1100) UE,VE,PE

      IF(NWSEGWI.EQ.0) THEN     !If the first file, read until the end
         DO I=2,IDSETFLG
            READ(IEXT,END=1100) IYEAR,IMONTH,IDAY,IHOUR
            READ(IEXT,END=1100) UE,VE,PE
         ENDDO
      ENDIF

 1100 IF(IDSETFLG.EQ.8) CLOSE(IEXT)

C.....Interpolate onto ADCIRC grid
C.....Convert velocity from the E grid reference to a lat/lon reference
C.....Convert pressure from millibars to N/M^2 to M of H20

      DO N=1,NP
         P1=PE(N1(N))
         P2=PE(N2(N))
         P3=PE(N3(N))
         U1=UE(N1(N))
         U2=UE(N2(N))
         U3=UE(N3(N))
         V1=VE(N1(N))
         V2=VE(N2(N))
         V3=VE(N3(N))
         UE29=U1*D1(N)+U2*D2(N)+U3*D3(N)
         VE29=V1*D1(N)+V2*D2(N)+V3*D3(N)
         CBETAU=COS(BETAU(N))
         SBETAU=SIN(BETAU(N))
         ULL(N)=UE29*CBETAU - VE29*SBETAU
         VLL(N)=UE29*SBETAU + VE29*CBETAU
         PLL(N)=P1*D1(N)+P2*D2(N)+P3*D3(N)
         PLL(N)=PLL(N)/RHOWATG100
      END DO

      RETURN
      END SUBROUTINE



C***********************************************************************
C  Subroutine to find where a given lon,lat falls in the Eta29 grid,   *
C     determine the interpolating factors to interpolate Eta29 fields  *
C     to that position, and finally to compute the angle to rotate the *
C     Eta29 velocity field to get to a lon, lat coordinated system.    *
C                                                                      *
C                    Written by R.L.       1/12/98                     *
C***********************************************************************

      subroutine e29search(node,FLON,FLAT,NN1,NN2,NN3,DD1,DD2,DD3,betau)
      implicit none
      integer nn1,nn2,nn3,node,icode,nwlon,nwlat,ifflag
      integer i,j,im2,jm2,n,ia,ja,na,ib,jb,nb,ic,jc,nc,id,jd,nd,
     &  ie,je,ne,ig,jg,ng,if
      real(sz) dd1,dd2,dd3,betau,ri,x1,x2,x3,x4,y1,y2,y3,y4
      real(sz) aemin,areas,a1,a2,a3,aa,ae,lambda
      real(8) lamda0,phi0,rphi0,cphi0,sphi0,tphi0,dlamda,dphi,rdlamda,
     &       rdphi,rflat,tflat,sflat,cflat,a,rlamar,cphiicrlamda,phiarg,
     &       rphii,rlamda,ri1,ri2,rj,dgtora,flon,flat
      real(sz) lamda,lamdaa,lamdab,lamdac,lamdad,lamdae,lamdag
      real(sz) phi,phia,phib,phic,phid,phie,phig
c
      icode=0
      nwlon=181
      nwlat=271
      dgtora=deg2rad
      lamda0=-97.0d0
      phi0=41.0d0
      rphi0=dgtora*phi0
      cphi0=cos(rphi0)
      sphi0=sin(rphi0)
      tphi0=tan(rphi0)
      dlamda=7.d0/36.d0
      dphi=5.d0/27.d0
      rdlamda=dgtora*dlamda
      rdphi=dgtora*dphi
c
      rflat=flat*dgtora
        tflat=tan(rflat)
      sflat=sin(rflat)
      cflat=cos(rflat)

c     compute the position of the closest node in the E29 grid

      a=flon-lamda0
      rlamar=cos(a*dgtora)
      cphiicrlamda=(rlamar+tflat*tphi0)*cflat*cphi0
      phiarg=sflat
      rphii=asin((phiarg-sphi0*cphiicrlamda)/cphi0)
      rlamda=acos(cphiicrlamda/cos(rphii))
      if(flon.lt.lamda0) rlamda=-rlamda
c
      ri2=(rlamda/rdlamda+nwlon+1)/2.
      ri1=(rlamda/rdlamda+nwlon)/2.
      rj=rphii/rdphi+(nwlat+1)/2
      j=(rj+0.5d0)
      ri=ri1
      if(mod(j,2).eq.0) ri=ri2
      i=(ri+0.5d0)

c     if (myproc == 0) then
c       write(screenunit,*) "lamda, phi = ",flon,flat
c       write(screenunit,*) "ri1, ri2, ri, rj = ",ri1,ri2,ri,rj
c       write(screenunit,*) "i, j = ",i,j
c     endif

      if ((rj.lt.1).or.(rj.gt.nwlat)) then
c        write(333,*) 'ADCIRC grid node ',node,
c     &             ' falls outside of the ETA 29 grid'
        icode=1
        NN1=1
        NN2=1
        NN3=1
        DD1=0
        DD2=0
        DD3=0
        return
      endif

      if (mod(j,2).eq.0) then
         if ((ri.lt.1).or.(ri.gt.(nwlon+0.5d0))) then
c          write(333,*) 'ADCIRC grid node ',node,
c     &                 ' falls outside of the ETA 29 grid'
            icode=1
            NN1=1
            NN2=1
            NN3=1
            DD1=0
            DD2=0
            DD3=0
            return
         endif
      endif

      if (mod(j,2).ne.0) then
         if ((ri.lt.0.5).or.(ri.gt.nwlon)) then
c           write(333,*) 'ADCIRC grid node ',node,
c     &                 ' falls outside of the ETA 29 grid'
            icode=1
            NN1=1
            NN2=1
            NN3=1
            DD1=0
            DD2=0
            DD3=0
            return
         endif
      endif

c     compute the coordinates of the closest Eta29 grid node

      jm2=(nwlat+1)/2
      im2=nwlon*2
      call e29calc(i,j,lamda,phi,n)

c     compute the coordinates of neighbor node "a" (located SW of closest node)

      if ((i.eq.1).and.(mod(j,2).eq.0)) then
         ia=i
         ja=j-2
      else
         ia=i
         if(mod(j,2).eq.0) ia=i-1
         ja=j-1
      endif
c                                 this neighbor lies outside of Eta29 grid
      if ((ia.lt.1).or.(ja.lt.1)) then
         na=0
      else
         call e29calc(ia,ja,lamdaa,phia,na)
      endif

c     compute the coordinates of neighbor node "b" (located W of closest node)

      ib=i-1
      jb=j
      if (ib.lt.1) then         !this neighbor lies outside of Eta29 grid
         nb=0
      else
         call e29calc(ib,jb,lamdab,phib,nb)
      endif

c     compute the coordinates of neighbor node "c" (located NW of closest node)

      if ((i.eq.1).and.(mod(j,2).eq.0)) then
         ic=i
         jc=j+2
      else
         ic=ia
         jc=j+1
      endif
c                                    this neighbor lies outside of Eta29 grid
      if ((ic.lt.1).or.(jc.gt.nwlat)) then
         nc=0
      else
         call e29calc(ic,jc,lamdac,phic,nc)
      endif

c     compute the coordinates of neighbor node "d" (located NE of closest node)

      if ((i.eq.181).and.(mod(j,2).ne.0)) then
         id=i
         jd=j+2
      else
         id=ic+1
         jd=j+1
      endif
c                                    this neighbor lies outside of Eta29 grid
      if ((id.gt.nwlon).or.(jd.gt.nwlat)) then
         nd=0
      else
         call e29calc(id,jd,lamdad,phid,nd)
      endif

c     compute the coordinates of neighbor node "e" (located E of closest node)

      ie=i+1
      je=j
      if (ie.gt.nwlon) then     !this neighbor lies outside of Eta29 grid
         ne=0
      else
         call e29calc(ie,je,lamdae,phie,ne)
      endif

c     compute the coordinates of neighbor node "g" (located SE of closest node)

      if ((i.eq.181).and.(mod(j,2).ne.0)) then
         ig=i
         jg=j-2
      else
         ig=id
         jg=j-1
      endif
c                                    this neighbor lies outside of Eta29 grid
      if ((ig.gt.nwlon).or.(jg.lt.1)) then
         ng=0
      else
         call e29calc(ig,jg,lamdag,phig,ng)
      endif

c     if (myproc == 0) then
c      write(screenunit,*) 'closest E29 node i,j = ',n,i,j,lamda,phi
c      if(na.eq.0) write(screenunit,*) 'point a falls outside of Eta29 grid'
c      if(na.ne.0) write(screenunit,*) 'point a   = ',na,ia,ja,lamdaa,phia
c      if(nb.eq.0) write(screenunit,*) 'point b falls outside of Eta29 grid'
c      if(nb.ne.0) write(screenunit,*) 'point b   = ',nb,ib,jb,lamdab,phib
c      if(nc.eq.0) write(screenunit,*) 'point c falls outside of Eta29 grid'
c      if(nc.ne.0) write(screenunit,*) "point c   = ",nc,ic,jc,lamdac,phic
c      if(nd.eq.0) write(screenunit,*) 'point d falls outside of Eta29 grid'
c      if(nd.ne.0) write(screenunit,*) "point d   = ",nd,id,jd,lamdad,phid
c      if(ne.eq.0) write(screenunit,*) 'point e falls outside of Eta29 grid'
c      if(ne.ne.0) write(screenunit,*) "point e   = ",ne,ie,je,lamdae,phie
c      if(ng.eq.0) write(screenunit,*) 'point g falls outside of Eta29 grid'
c      if(ng.ne.0) write(screenunit,*) "point g   = ",ng,ig,jg,lamdag,phig
c     endif

      NN1=1
      NN2=1
      NN3=1
      DD1=0
      DD2=0
      DD3=0
      X1=lamda
      X4=flon
      Y1=phi
      Y4=flat
      ifflag=0
      AEMIN=99999.d0

c     test if the point is in triangle ij - b - a

      if ((na.ne.0).and.(nb.ne.0)) then
         X2=lamdab
         X3=lamdaa
         Y2=phib
         Y3=phia
         AREAS=ABS((X1-X3)*(Y2-Y3)+(X3-X2)*(Y1-Y3))
         A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
         A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
         A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
         AA=ABS(A1)+ABS(A2)+ABS(A3)
         AE=ABS(AA-AREAS)/AREAS
c     write(333,*) "AE = ",AE
         IF((AE.LT.1.0d-5).AND.(AE.LT.AEMIN)) THEN
            AEMIN=AE
            NN1=n
            NN2=nb
            NN3=na
            DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
            DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
            DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
            call betaucalc(i,j,DD1,ib,jb,DD2,ia,ja,DD3,betau)
            ifflag=ifflag+1
c     write(333,*) 'position found in triangle ij - b - a'
         ENDIF
      endif

c     if along the west boundary, test if the point is in triangle ij - c - a

      if((i.eq.1).and.(mod(j,2).ne.0)) then
         if((na.ne.0).and.(nc.ne.0)) then
            X2=lamdac
            X3=lamdaa
            Y2=phic
            Y3=phia
            AREAS=ABS((X1-X3)*(Y2-Y3)+(X3-X2)*(Y1-Y3))
            A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
            A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
            A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
            AA=ABS(A1)+ABS(A2)+ABS(A3)
            AE=ABS(AA-AREAS)/AREAS
c     write(333,*) "AE = ",AE
            IF((AE.LT.1.0d-5).AND.(AE.LT.AEMIN)) THEN
               NN1=n
               NN2=nc
               NN3=na
               DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
               DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
               DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
               call betaucalc(i,j,DD1,ic,jc,DD2,ia,ja,DD3,betau)
               ifflag=ifflag+1
c     write(333,*) 'position found in triangle ij - c - a'
            ENDIF
         endif
      endif

c     test if the point is in triangle ij - c - b

      if((nb.ne.0).and.(nc.ne.0)) then
         X2=lamdac
         X3=lamdab
         Y2=phic
         Y3=phib
         AREAS=ABS((X1-X3)*(Y2-Y3)+(X3-X2)*(Y1-Y3))
         A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
         A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
         A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
         AA=ABS(A1)+ABS(A2)+ABS(A3)
         AE=ABS(AA-AREAS)/AREAS
c     write(333,*) "AE = ",AE
         IF((AE.LT.1.0d-5).AND.(AE.LT.AEMIN)) THEN
            NN1=n
            NN2=nc
            NN3=nb
            DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
            DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
            DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
            call betaucalc(i,j,DD1,ic,jc,DD2,ib,jb,DD3,betau)
            ifflag=ifflag+1
c     write(333,*) 'position found in triangle ij - c - b'
         ENDIF
      endif

c     test if the point is in triangle ij - d - c

      if((nc.ne.0).and.(nd.ne.0)) then
         X2=lamdad
         X3=lamdac
         Y2=phid
         Y3=phic
         AREAS=ABS((X1-X3)*(Y2-Y3)+(X3-X2)*(Y1-Y3))
         A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
         A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
         A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
         AA=ABS(A1)+ABS(A2)+ABS(A3)
         AE=ABS(AA-AREAS)/AREAS
c     write(333,*) "AE = ",AE
         IF((AE.LT.1.0d-5).AND.(AE.LT.AEMIN)) THEN
            NN1=n
            NN2=nd
            NN3=nc
            DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
            DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
            DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
            call betaucalc(i,j,DD1,id,jd,DD2,ic,jc,DD3,betau)
            ifflag=ifflag+1
c     write(333,*) 'position found in triangle ij - d - c'
         ENDIF
      endif

c     if along the east boundary, test if the point is in triangle ij - g - d

      if((i.eq.181).and.(mod(j,2).eq.0)) then
         if((nd.ne.0).and.(ng.ne.0)) then
            X2=lamdag
            X3=lamdad
            Y2=phig
            Y3=phid
            AREAS=ABS((X1-X3)*(Y2-Y3)+(X3-X2)*(Y1-Y3))
            A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
            A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
            A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
            AA=ABS(A1)+ABS(A2)+ABS(A3)
            AE=ABS(AA-AREAS)/AREAS
c     write(333,*) "AE = ",AE
            IF((AE.LT.1.0d-5).AND.(AE.LT.AEMIN)) THEN
               NN1=n
               NN2=ng
               NN3=nd
               DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
               DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
               DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
               call betaucalc(i,j,DD1,ig,jg,DD2,id,jd,DD3,betau)
               ifflag=ifflag+1
c     write(333,*) 'position found in triangle ij - g - d'
            ENDIF
         endif
      endif

c     test if the point is in triangle ij - e - d

      if((nd.ne.0).and.(ne.ne.0)) then
         X2=lamdae
         X3=lamdad
         Y2=phie
         Y3=phid
         AREAS=ABS((X1-X3)*(Y2-Y3)+(X3-X2)*(Y1-Y3))
         A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
         A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
         A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
         AA=ABS(A1)+ABS(A2)+ABS(A3)
         AE=ABS(AA-AREAS)/AREAS
c     write(333,*) "AE = ",AE
         IF((AE.LT.1.0d-5).AND.(AE.LT.AEMIN)) THEN
            NN1=n
            NN2=ne
            NN3=nd
            DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
            DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
            DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
            call betaucalc(i,j,DD1,ie,je,DD2,id,jd,DD3,betau)
            ifflag=ifflag+1
c     write(333,*) 'position found in triangle ij - e - d'
         ENDIF
      endif

c     test if the point is in triangle ij - g - e

      if((ne.ne.0).and.(ng.ne.0)) then
         X2=lamdag
         X3=lamdae
         Y2=phig
         Y3=phie
         AREAS=ABS((X1-X3)*(Y2-Y3)+(X3-X2)*(Y1-Y3))
         A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
         A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
         A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
         AA=ABS(A1)+ABS(A2)+ABS(A3)
         AE=ABS(AA-AREAS)/AREAS
c     write(333,*) "AE = ",AE
         IF((AE.LT.1.0d-5).AND.(AE.LT.AEMIN)) THEN
            NN1=n
            NN2=ng
            NN3=ne
            DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
            DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
            DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
            call betaucalc(i,j,DD1,ig,jg,DD2,ie,je,DD3,betau)
            ifflag=ifflag+1
c     write(333,*) 'position found in triangle ij - g - e'
         ENDIF
      endif

c     test if the point is in triangle ij - a - g

      if((na.ne.0).and.(ng.ne.0)) then
         X2=lamdaa
         X3=lamdag
         Y2=phia
         Y3=phig
         AREAS=ABS((X1-X3)*(Y2-Y3)+(X3-X2)*(Y1-Y3))
         A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
         A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
         A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
         AA=ABS(A1)+ABS(A2)+ABS(A3)
         AE=ABS(AA-AREAS)/AREAS
c     write(333,*) "AE = ",AE
         IF((AE.LT.1.0d-5).AND.(AE.LT.AEMIN)) THEN
            NN1=n
            NN2=na
            NN3=ng
            DD1=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS
            DD2=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS
            DD3=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS
            call betaucalc(i,j,DD1,ia,ja,DD2,ig,jg,DD3,betau)
            ifflag=ifflag+1
c     write(333,*) 'position found in triangle ij - a - g'
         ENDIF
      endif

c      if(ifflag.eq.0) then
c         write(333,*) 'position not found'
c        if (myproc == 0) then
c         write(screenunit,*) 'position not found in subroutine E29SEARCH'
c        endif
c        icode=3
c      else
c        if (myproc == 0) then
c         write(screenunit,*) 'i,j,NN1,NN2,NN3,DD1,DD2,DD3'
c        endif
c         write(333,999) i,j,NN1,NN2,NN3,DD1,DD2,DD3,betau/dgtora
c 999     format(5I8,1x,3E13.6)
c      endif

      return
      end subroutine



C***********************************************************************
C  Subroutine to compute the longititude and latitude of a given i,j   *
C       position in the Eta29 grid.                                    *
C                                                                      *
C                    Written by R.L.       1/11/98                     *
C***********************************************************************

      subroutine e29calc(i,j,lamda,phi,n)
      implicit none
      integer i,j,n,nwlon,nwlat,im2,jm2,i1,i2,i1p1,i1m1,i2p1,i2m1,
     &     i3p1,i3m1
      real(sz) lamda,phi,phii,dlon,dlat,dlnt,arg,betau1,betau2,betau3
      real(8) lamda0,phi0,rphi0,cphi0,sphi0,tphi0,dlamda,dphi,rdlamda,
     &     rdphi,a,rlamar,phiarg,rlamda,dgtora
c
      nwlon=181
      nwlat=271
      dgtora=deg2rad
      lamda0=-97.0d0
      phi0=41.0d0
      rphi0=dgtora*phi0
      cphi0=cos(rphi0)
      sphi0=sin(rphi0)
      tphi0=tan(rphi0)
      dlamda=7.d0/36.d0
      dphi=5.d0/27.d0
      rdlamda=dgtora*dlamda
      rdphi=dgtora*dphi
c
      jm2=(nwlat+1)/2
      im2=nwlon*2
c
      phii=rdphi*real(j-jm2,kind=8)
      i1=2*i-1
      i2=2*i
      if(mod(j,2).ne.0) then
         rlamda=rdlamda*real(i2-nwlon,kind=8)
      else
         rlamda=rdlamda*real(i1-nwlon,kind=8)
      endif
      phiarg= sin(phii)*cphi0+cos(phii)*sphi0*cos(rlamda)
      if(phiarg.gt.1.0d0) phiarg=1.0d0
      if(phiarg.lt.-1.0d0) phiarg=-1.0d0
      phi=asin(phiarg)
      rlamar= cos(phii)*cos(rlamda)/(cos(phi)*cphi0)-tan(phi)*tphi0
      if(rlamar.gt.1.0d0) rlamar=1.0d0
      if(rlamar.lt.-1.d0) rlamar=-1.d0
      a=acos(rlamar)/dgtora
      if(rlamda.le.0.) then
         lamda=lamda0-a
      else
         lamda=lamda0+a
      endif
      phi=phi/dgtora
      n=nwlon*(j-1)+i
C
      return
      end subroutine


C***********************************************************************
C  Subroutine to compute the conversion angle between the E29 velocity *
C       field and a lon,lat coordinate system.                         *
C                                                                      *
C                    Written by R.L.       1/12/98                     *
C***********************************************************************

      subroutine betaucalc(i1,j1,dd1,i2,j2,dd2,i3,j3,dd3,betau)
      implicit none
      integer i1,j1,i2,j2,i3,j3,n,i1p1,i1m1,i2p1,i2m1,i3p1,i3m1
      real(sz) dd1,dd2,dd3,betau
      real(sz) lamda,lamdap1,lamdam1,phi,phip1,phim1,dlon,dlat,
     &     dlnt,arg,betau1,betau2,betau3,dgtora
c
      dgtora=deg2rad
c
      if(i1.ne.181) then
         i1p1=i1+1
      else
         i1p1=i1
      endif
      if(i1.ne.1) then
         i1m1=i1-1
      else
         i1m1=i1
      endif
      call e29calc(i1,j1,lamda,phi,n)
      call e29calc(i1p1,j1,lamdap1,phip1,n)
      call e29calc(i1m1,j1,lamdam1,phim1,n)
      dlon=(lamdap1-lamdam1)*cos(phi*dgtora)
      dlat=phip1-phim1
      dlnt=sqrt(dlon*dlon+dlat*dlat)
      arg=dlat/dlnt
      if(arg.gt.1.d0) arg=1.d0
      if(arg.lt.-1.d0) arg=-1.d0
      betau1=asin(arg)
c
      if(i2.ne.181) then
         i2p1=i2+1
      else
         i2p1=i2
      endif
c
      if(i2.ne.1) then
         i2m1=i2-1
      else
         i2m1=i2
      endif
c
      call e29calc(i2,j2,lamda,phi,n)
      call e29calc(i2p1,j2,lamdap1,phip1,n)
      call e29calc(i2m1,j2,lamdam1,phim1,n)
      dlon=(lamdap1-lamdam1)*cos(phi*dgtora)
      dlat=phip1-phim1
      dlnt=sqrt(dlon*dlon+dlat*dlat)
      arg=dlat/dlnt
      if(arg.gt.1.d0) arg=1.d0
      if(arg.lt.-1.d0) arg=-1.d0
      betau2=asin(arg)
c
      if(i3.ne.181) then
         i3p1=i3+1
      else
         i3p1=i3
      endif
c
      if(i3.ne.1) then
         i3m1=i3-1
      else
         i3m1=i3
      endif
c
      call e29calc(i3,j3,lamda,phi,n)
      call e29calc(i3p1,j3,lamdap1,phip1,n)
      call e29calc(i3m1,j3,lamdam1,phim1,n)
      dlon=(lamdap1-lamdam1)*cos(phi*dgtora)
      dlat=phip1-phim1
      dlnt=sqrt(dlon*dlon+dlat*dlat)
      arg=dlat/dlnt
      if(arg.gt.1.d0) arg=1.d0
      if(arg.lt.-1.d0) arg=-1.d0
      betau3=asin(arg)
      betau=dd1*betau1+dd2*betau2+dd3*betau3
C
      return
      end subroutine

C     ----------------------------------------------------------------
C      S U B R O U T I N E   H O L L A N D  G E T
C     ----------------------------------------------------------------
C
C     jgf46.05 Subroutine to calculate wind velocity at nodes from
C     the Holland Wind model.
C
C     The format statement takes into account whether the track data is
C     hindcast/nowcast (BEST) or forecast (OFCL).
C
C     The first line in the file MUST be a hindcast, since the central
C     pressure and the RMW are carried forward from hindcasts into
C     forecasts. So there needs to be at least one hindcast to carry the data
C     forward.
C
C     Assumes spherical coordinates (ICS=2 in fort.15 file).
C
C     Based on bob's NWS67GET (See below).
C
C     ----------------------------------------------------------------
      SUBROUTINE HollandGet(X,Y,SLAM,SFEA,WVNX,WVNY,PRESS,NP,
     &     ICS,RHOWAT0,G,TIMELOC,NSCREEN,ScreenUnit)
      USE SIZES, ONLY : MyProc, LOCALDIR
      IMPLICIT NONE
      INTEGER, intent(in) :: NP,ICS
      REAL(SZ), intent(in) :: RHOWAT0,G,TIMELOC
      REAL(8), intent(in), dimension(NP) :: X, Y, SLAM, SFEA
      REAL(SZ), intent(out), dimension(NP) :: WVNX,WVNY
      REAL(SZ), intent(out), dimension(NP) :: PRESS
      INTEGER, intent(in) :: NSCREEN
      INTEGER, intent(in) :: ScreenUnit
      INTEGER I, J
      REAL(SZ),SAVE,ALLOCATABLE :: RAD(:),DX(:),DY(:)
      REAL(SZ),SAVE,ALLOCATABLE :: XCOOR(:),YCOOR(:) !lat and lon of nodes
      REAL(SZ),SAVE,ALLOCATABLE :: V_r(:)            !wind sp. as fn. of dist
      REAL(SZ),SAVE,ALLOCATABLE :: THETA(:)
      REAL(SZ) :: TVX,TVY,RRP,RMW,A,B,WTRATIO
      REAL(SZ) :: TransSpdX, TransSpdY
      REAL(SZ) :: TM, cpress, lon, lat, spd
      REAL(SZ) :: ts        ! storm translation speed, m/s
      REAL(SZ) :: WindMultiplier         ! for storm 2 in LPFS ensemble
      REAL(SZ) :: centralPressureDeficit ! difference btw ambient and cpress
C
      REAL(8) :: omega, coriolis
      REAL(8) :: mperdeg

      LOGICAL, SAVE :: FIRSTCALL = .True.
C
      mperdeg  = Rearth * pi / 180.0d0
      omega = 2.0d0*pi / 86164.2d0
C
      IF (FIRSTCALL) THEN
         FIRSTCALL = .False.
         ALLOCATE (RAD(NP),DX(NP),DY(NP),XCOOR(NP),YCOOR(NP))
         ALLOCATE (V_r(NP),THETA(NP))
C
C     The subroutine only works for ICS=2 (spherical coordinates)
         DO I=1,NP
            XCOOR(I)=SLAM(I)*RAD2DEG
            YCOOR(I)=SFEA(I)*RAD2DEG
         END DO
      ENDIF
C
C     Get data for this time step.
      CALL GetHollandStormData(lat,lon,cpress,spd,rrp,rmw,tvx,tvy,TIMELOC,
     &     nscreen,screenunit)
C
C     jgf50.32: If we are using sector-based wind drag, record the location
C     of the center of the storm.
      IF ((lat.ne.EyeLatR(3)).or.(lon.ne.EyeLonR(3))) THEN
         EyeLatR(1) = EyeLatR(2)
         EyeLonR(1) = EyeLonR(2)
         EyeLatR(2) = EyeLatR(3)
         EyeLonR(2) = EyeLonR(3)
         EyeLatR(3) = lat
         EyeLonR(3) = lon
         FoundEye = .true.
      ENDIF
C
C     Calculate and limit central pressure deficit; some track files
C     (e.g., Charley 2004) may have a central pressure greater than the
C     ambient pressure that this subroutine assumes

C     jgf51.14: Apply a factor of 100.0 to convert PRBCKGRND from
C     mb to pa for use in this subroutine.
      centralPressureDeficit = PRBCKGRND * 100.d0 - cpress
      IF ( centralPressureDeficit .lt. 100.d0 ) THEN
         centralPressureDeficit = 100.d0
      ENDIF
C
C     jgf46.29 Subtract the translational speed of the storm from the
C     observed max wind speed to avoid distortion in the Holland curve
C     fit. The translational speed will be added back later.
      ts=sqrt(tvx*tvx+tvy*tvy)
      spd=spd-ts
C
C     Convert wind speed from 10 meter altitude (which is what the
C     NHC forecast contains) to wind speed at the top of the atmospheric
C     boundary layer (which is what the Holland curve fit requires).
      spd=spd/BLAdj
C
C     Calculate Holland parameters and limit the result to its appropriate
C     range.
      B = rhoAir*e*(spd**2.d0)/(centralPressureDeficit)
      IF (B.lt.1.0d0) B=1.0d0
      IF (B.gt.2.5d0) B=2.5d0
C
C     Calculate Holland A parameter. (jgf46.32jgf4 commented out)
C      A = (RMW*1000.d0)**B
C
#ifdef DEBUG_HOLLAND
      WRITE(16,4321) B
 4321 FORMAT(/,2x,'Holland B parameter is ',e16.8)
#endif
C     jgf46.28 If we are running storm 2 in the Lake Pontchartrain
C     Forecast System ensemble, the final wind speeds should be
C     multiplied by 1.2.
      IF (StormNumber.eq.2) THEN
         WindMultiplier=1.2d0
      ELSE
         WindMultiplier=1.0d0
      ENDIF
C
C     Calculate wind velocity and pressure at each node.
      DO I=1,NP
         DX(I)=(XCOOR(I)-lon)*DEG2RAD
         DY(I)=(YCOOR(I)-lat)*DEG2RAD
         THETA(I)=ATAN2(DY(I),DX(I))
C   RJW v48.45
C     compute the distances based on haversine formula for distance along a sphere
        rad(i) = sphericalDistance(DX(I),DY(I),lat, ycoor(i))
C       Rad(i)=Rearth*(2.0d0*ASIN(sqrt(sin(DY(I)/2.0d0)**2.0d0+
C     & cos(lat*DEG2RAD)*cos(YCOOR(i)*DEG2RAD)*sin(DX(I)/2.0d0)**2.0d0)))
C compute coriolis
        coriolis = 2.0d0 * omega * sin(YCOOR(I)*DEG2RAD)

        PRESS(I)=(cpress+(centralPressureDeficit)*
     &        EXP(-(RMW*1000.d0/RAD(I))**B)) / (RHOWAT0*G)
         V_r(I) = sqrt(
     &        (RMW*1000.d0/RAD(I))**B *
     &        EXP(1.d0-(RMW*1000.d0/RAD(I))**B)*spd**2.d0
     &        + (RAD(I)**2.d0)*(CORIOLIS**2.d0)/4.d0
     &        )
     &        - RAD(I)*CORIOLIS/2.d0
C
C     jgf46.31 Determine translation speed that should be added to final
C     storm wind speed. This is tapered to zero as the storm wind tapers
C     to zero toward the eye of the storm and at long distances from the
C     storm.
         TransSpdX = (abs(V_r(I))/spd)*TVX
         TransSpdY = (abs(V_r(I))/spd)*TVY
C     Apply mutliplier for Storm2 in LPFS ensemble.
         V_r(I) = V_r(I) * WindMultiplier
C
C     Find the velocity components.
         WVNX(I)=-V_r(I)*SIN(THETA(I))
         WVNY(I)= V_r(I)*COS(THETA(I))
C
C     jgf46.19 Convert wind velocity from top of atmospheric boundary
C     layer (which is what the Holland curve fit produces) to wind
C     velocity at 10m above the earth's surface (factor of
C     0.7).
         WVNX(I)=WVNX(I)*BLAdj
         WVNY(I)=WVNY(I)*BLAdj
C
C     jgf46.21 Also convert from 1 minute averaged winds to 10
C     minute averaged winds (0.88 factor).
         WVNX(I)=WVNX(I)*one2ten
         WVNY(I)=WVNY(I)*one2ten
C
C     jgf46.31 Add the storm translation speed.
         WVNX(I)=WVNX(I)+TransSpdX
         WVNY(I)=WVNY(I)+TransSpdY
C
C     jgf46.31 Set the wind velocities to zero outside the last closed
C     isobar.
C         IF (RAD(I).gt.rrp) THEN
C            WVNX(I)=0.0d0
C            WVNY(I)=0.0d0
C         ENDIF
C
      END DO
C     ----------------------------------------------------------------
      END SUBROUTINE HollandGet
C     ----------------------------------------------------------------

C     ----------------------------------------------------------------
C      S U B R O U T I N E   G E T  H O L L A N D  S T O R M  D A T A
C     ----------------------------------------------------------------
C
C     jgf46.08 Subroutine to support HollandGet. Gets the next line from
C     the file, skipping lines that are time repeats. Interpolates in
C     time if we are between wind data points. Does conversions to the
C     proper units. Uses old values of central pressure and RMW if the
C     line is a forecast, since forecasts do not have that data in them.
C     Assumes longitude is WEST longitude, latitude is NORTH latitude.
C
C     ----------------------------------------------------------------
C
      SUBROUTINE GetHollandStormData(LatOut,LonOut,CPressOut,SpdOut,
     &     RRPOut,RMWOut,TVXOut,TVYOut,TIMELOC,NScreen,ScreenUnit)
      USE SIZES, ONLY : MyProc, GBLINPUTDIR
      USE VORTEX, ONLY : uvtrans
      IMPLICIT NONE
      REAL(SZ),intent(in) :: TIMELOC
      INTEGER,intent(in) :: nscreen, screenunit
      REAL(SZ),intent(out) :: LatOut, LonOut, CPressOut
      REAL(SZ),intent(out) :: SpdOut, RRPOut, RMWOut, TVXOut, TVYOut
C
      INTEGER, ALLOCATABLE, SAVE :: iYear(:),iMth(:),iDay(:),iHr(:)
      INTEGER, ALLOCATABLE, SAVE :: iLat(:),iLon(:)
      CHARACTER(1), ALLOCATABLE, SAVE :: ns(:), ew(:)
      INTEGER, ALLOCATABLE, SAVE :: iSpd(:),iCPress(:),iRRP(:),iRMW(:)
      REAL(SZ), ALLOCATABLE, SAVE :: Lat(:),Lon(:),Spd(:)
      REAL(SZ), ALLOCATABLE, SAVE :: CPress(:),RRP(:),RMW(:)
      REAL(SZ), ALLOCATABLE, SAVE :: TVX(:), TVY(:)      ! jgf46.32jgf9
      CHARACTER(len=4), ALLOCATABLE, SAVE :: CastType(:) !hindcast,forecast
      INTEGER, ALLOCATABLE, SAVE :: iFcstInc(:) ! hours
      REAL(SZ), ALLOCATABLE, SAVE :: FcstInc(:) ! seconds
      REAL(SZ),SAVE :: WTRATIO          !time ratio used for interpolation
      REAL(8), ALLOCATABLE, SAVE :: CastTime(:) ! seconds since start of year
      INTEGER,SAVE :: iNowcastCPress, iNowcastRRP, iNowcastRMW
      LOGICAL,SAVE :: FIRSTCALL = .True.
      INTEGER,SAVE :: i  ! Current array counter for fort.22 file
      INTEGER,SAVE :: nl ! Number of lines in the fort.22 file
      INTEGER,SAVE :: pl ! populated length of Holland Data array
      INTEGER :: j ! loop counter
      INTEGER :: ios ! return code for an i/o operation
C     ------------------------------------------------------
C     BEGIN Code executed upon first call to this subroutine
C     ------------------------------------------------------
      IF (FIRSTCALL) THEN
C
C     Determine the number of lines in the file.
         nl=0
         call openFileForRead(22,TRIM(GBLINPUTDIR)//'/'//'fort.22',ios)
         if ( ios.gt.0 ) then
            call allMessage(ERROR,
     &       "The symmetric vortex parameter file was not found. "
     &       //"ADCIRC terminating.")
            call windTerminate()
         endif

         DO
            READ(UNIT=22,FMT='(A170)',END=8888)
            nl=nl+1
         ENDDO
 8888    REWIND(22)
C
C     Dimension the arrays according to the number of lines in the file,
C     this will be greater than or equal to the size of the array we need
C     (probably greater because of the repeated lines that we throw away)
         ALLOCATE(iYear(nl),iMth(nl),iDay(nl),iHr(nl),iLat(nl),iLon(nl),
     &        iSpd(nl),iCpress(nl),iRRP(nl),iRMW(nl),iFcstInc(nl))
         ALLOCATE(ns(nl),ew(nl))
         ALLOCATE(Lat(nl),Lon(nl),Spd(nl),CPress(nl),RRP(nl),RMW(nl),
     &        FcstInc(nl),TVX(nl),TVY(nl))
         ALLOCATE(CastType(nl))
         ALLOCATE(CastTime(nl))
C
C     Now read the data into the arrays. The first
C     line must be a hindcast/nowcast.
         i=1
C
         DO
C     Get another line of data from the file and check to see if the
C     line represents a new point in time, or is a repeated time
C     point. Repeated time points occur in hindcasts for the purpose of
C     describing winds in the quadrants of the storm. We don't use the
C     quadrant-by-quadrant wind data. Repeated time data occur in the
C     forecast because the time data is just the time that the forecast
C     was made. The important parameter in the forecast file is the
C     forecast increment.
            READ(UNIT=22,FMT=228,END=9999)
     &           iYear(i),iMth(i),iDay(i),iHr(i),
     &           CastType(i),iFcstInc(i),iLat(i),ns(i),iLon(i),ew(i),
     &           iSpd(i),iCPress(i),iRRP(i),iRMW(i)

C                yr,mo,dy,hr, ,type, inc,  lat,NS,  lon,EW,  spd,   pc,
 228        format(8x,i4,i2,i2,i2,6x,a4,2x,i3,1x,i4,a1,2x,i4,a1,2x,i3,2x,i4,
     &     47x,i3,2x,i3)
C              RRP,  RMW

C
C
            SELECT CASE(trim(CastType(i)))
C           ------------
            CASE("BEST")     ! nowcast/hindcast
C           ------------
C     Check to see if this is a repeated line. If so, go directly to the
C     next line without any processing.
               IF (i.gt.1) THEN
               IF ( iYear(i).eq.iYear(i-1).and.iMth(i).eq.iMth(i-1).and.
     &              iDay(i).eq.iDay(i-1).and.iHr(i).eq.iHr(i-1)) THEN
                  CYCLE
               ENDIF
               ENDIF
C
C     Save the central pressure, radius of last closed isobar, and
C     radius to max wind for use in forecasts
               iNowcastCPress=iCPress(i)
               iNowcastRMW=iRMW(i)
               iNowcastRRP=iRRP(i)
C
C     Determine the time of this hindcast in seconds since the beginning
C     of the year.
               CALL TimeConv(iYear(i),iMth(i),iDay(i),iHr(i),0,0.d0,
     &              CastTime(i),MyProc,NScreen,ScreenUnit)
C
C     Determine the CastTime in seconds since the beginning of the simulation.
               CastTime(i)=CastTime(i)-WindRefTime
               FcstInc(i)=iFcstInc(i)
C
C           ------------
            CASE("OFCL")        ! forecast
C           ------------
C     Check to see if this is a repeated line (i.e., a forecast that
C     coincides with the nowcast, or a repeated forecast). If so, go
C     directly to the next line without any processing.
               IF (i.gt.1 ) THEN
               IF ((iFcstInc(i).eq.0.and.
     &              (iYear(i).eq.iYear(i-1).and.iMth(i).eq.iMth(i-1)
     &              .and.iDay(i).eq.iDay(i-1).and.iHr(i).eq.iHr(i-1)))
     &              .or.
     &             (iFcstInc(i).ne.0.and.iFcstInc(i).eq.iFcstInc(i-1)))
     &              THEN
                  CYCLE
               ENDIF
               ENDIF
               FcstInc(i) = iFcstInc(i)
C
C     Determine the time of this forecast in seconds since the beginning
C     of the year.
               IF ( iFcstInc(i).eq.0 ) THEN
                  CALL TimeConv(iYear(i),iMth(i),iDay(i),iHr(i),0,0.d0,
     &                 CastTime(i),MyProc,NScreen,ScreenUnit)
                  CastTime(i)=CastTime(i)-WindRefTime
               ELSE
                  FcstInc(i) = FcstInc(i) * 3600.d0 ! convert hours to seconds
                  CastTime(i) = CastTime(i-1) +
     &                 ( FcstInc(i) - FcstInc(i-1) )
               ENDIF
C
C     jgf48.4637 If the forecast values of central pressure or RMW are
C     zero (and they will be if unless the user has filled them in, because
C     the NHC does not forecast these parameters), exit with a fatal error.
               IF ( (iCPress(i).eq.0).or.(iRMW(i).eq.0) ) THEN
                  WRITE(16,1000)   ! unit 22 Holland Storm File
                  WRITE(16,1121)   ! cpress and/or Rmax were zero
                  WRITE(16,1131)   ! must both be non zero
                  IF (NScreen.ne.0.and.MyProc.eq.0) THEN
                     WRITE(ScreenUnit,1000) ! unit 22 Holland Storm File
                     WRITE(ScreenUnit,1021) CastType(i) ! cpress, Rmax were 0
                     WRITE(ScreenUnit,1131) ! describe valid input
                  ENDIF
#ifdef CMPI
                  CALL MSG_FINI()
#endif
                  STOP
               ENDIF
C
            CASE DEFAULT        ! unrecognized
               WRITE(16,1000)   ! unit 22 Holland Storm File
               WRITE(16,1021) CastType(i),MyProc ! contains invalid name
               WRITE(16,1031)   ! describe valid input
               WRITE(16,1041)   ! tell which column failed
               IF (NScreen.ne.0.and.MyProc.eq.0) THEN
                  WRITE(ScreenUnit,1000) ! unit 22 Holland Storm File
                  WRITE(ScreenUnit,1025) CastType(i) ! contains invalid name
                  WRITE(ScreenUnit,1031) ! describe valid input
                  WRITE(ScreenUnit,1041) ! tell which column failed
               ENDIF
               STOP
            END SELECT
C
C     Convert integers to reals.
            Lat(i) = iLat(i)
            Lon(i) = iLon(i)
            Spd(i) = iSpd(i)
            CPress(i) = iCPress(i)
            RRP(i) = iRRP(i)
            RMW(i) = iRMW(i)
C
C     Convert units.
            Lat(i) = Lat(i) / 10.d0 ! convert 10ths of degs to degs
            Lon(i) = Lon(i) / 10.d0 ! convert 10ths of degs to degs
            if ( ew(i).eq.'W' ) then
               Lon(i) = -1.d0 * Lon(i) ! negative b/c WEST longitude
            endif
            if ( ns(i).eq.'S' ) then
               Lat(i) = -1.d0 * Lat(i) ! negative b/c SOUTH latitude
            endif
            CPress(i) = CPress(i) * 100.d0 ! convert mbar to Pa
            RRP(i) = RRP(i) * 1.852000003180799d0 * 1000.0d0 ! convert nm to m
            RMW(i) = RMW(i) * 1.852000003180799d0 ! convert nm to km
            Spd(i) = Spd(i) * 0.51444444d0 ! convert kts to m/s
C
#ifdef DEBUG_HOLLAND
            WRITE(16,1244) CastTime(i),Lat(i),Lon(i),
     &           Spd(i),CPress(i),RRP(i),RMW(i),WindRefTime
            if ( i.gt.1 ) then
               write(16,2355) FcstInc(i-1)
            endif
            IF (NScreen.ne.0.and.MyProc.eq.0) THEN
               WRITE(ScreenUnit,1244) CastTime(i),Lat(i),Lon(i),
     &              Spd(i),CPress(i),RRP(i),RMW(i),WindRefTime
            ENDIF
 1244       FORMAT('CastTime ',e16.8,' Lat ',f6.2,' Lon ',f6.2,
     &           /,'Spd ',f8.2,' CPress ',f10.2,' RRP ',f12.2,
     &           ' RMW ',f8.2, ' WindRefTime ',e16.8)
 2355       format('FcstInc(i-1) ',e16.8)
#endif
C
C     Save the number of non-repeated lines from the fort.22 file, this
C     is the populated length of the array.
            pl=i
C
C     Increment array counter
            i=i+1
C
         ENDDO
 9999    CLOSE(22)
C
C     Calculate storm translation velocities based on change in position, then
C     convert degrees/time to m/s
         ! initialize storm translation speeds
         TVX = 0.0
         TVY = 0.0
         DO i=2, pl
C     Calculate storm translation velocities based on change in position,
            ! approximate u and v translation velocities
            call uvTrans(lat(i-1),lon(i-1),
     &                    lat(i),lon(i),
     &                    CastTime(i-1),CastTime(i),
     &                    tvx(i),tvy(i))
     
         ENDDO
         ! extrapolate later storm translation speeds back to the
         ! first storm position
         IF ( pl.ge.2 ) THEN
            TVX(1) = TVX(2)
            TVY(1) = TVY(2)
         ENDIF
C
C     Determine the correspondence between the current simulation time and
C     the fort.22 file.
         i=2
         DO
            IF (TIMELOC.ge.CastTime(i-1).and.TIMELOC.lt.CastTime(i)) THEN
               EXIT
            ELSE
               i=i+1
               IF (i.gt.pl) THEN
                  call allMessage(ERROR,
     &             "The Storm Hindcast/Forecast Input File (unit 22) "
     &             //"does not contain times/dates that correspond "
     &             //"to the ADCIRC current model time. "
     &             //" ADCIRC terminating.")
                  call windTerminate()
               ENDIF
            ENDIF
         ENDDO
         FIRSTCALL = .False.
      ENDIF
C     ----------------------------------------------------------
C     END Code executed only upon first call to this subroutine
C     ----------------------------------------------------------
C
C
C     ----------------------------------------------------------
C     BEGIN Code executed on every call to this subroutine
C     ----------------------------------------------------------
C
C     If time exceeds the next hindcast/nowcast/forecast time, increment the
C     array counter.
      IF (TIMELOC.gt.CastTime(i)) THEN
         i=i+1
C         jgf51.14: Check to see that we haven't gone off the end of 
C         meteorological data.
          IF (i.gt.pl) THEN
             call allMessage(ERROR,'The simulation time has extended '//
     &          'beyond the end of the meteorological dataset.')
             call windTerminate()
          ENDIF
      ENDIF
C
C     Interpolate w.r.t. time
      WTRATIO=(TIMELOC-CastTime(i-1))/(CastTime(i)-CastTime(i-1))
      CPressOut = CPress(i-1) + WTRATIO*(CPress(i)-CPress(i-1))
      LonOut = Lon(i-1) + WTRATIO*(Lon(i)-lon(i-1))
      LatOut = Lat(i-1) + WTRATIO*(Lat(i)-lat(i-1))
      SpdOut = Spd(i-1) + WTRATIO*(Spd(i)-Spd(i-1))
      RRPOut = RRP(i-1) + WTRATIO*(RRP(i)-RRP(i-1))
      RMWOut = RMW(i-1) + WTRATIO*(RMW(i)-RMW(i-1))
      TVXOut = TVX(i-1) + WTRATIO*(TVX(i)-TVX(i-1))
      TVYOut = TVY(i-1) + WTRATIO*(TVY(i)-TVY(i-1))
C
C     ----------------------------------------------------------
C     END Code executed on every call to this subroutine
C     ----------------------------------------------------------
 10   format(8x,i4,i2,i2,i2,6x,a4,7x,i3,4x,i3,3x,i3,2x,i4,52x,i3)
 12   format(8x,i4,i2,i2,i2,6x,a4,2x,i3,2x,i3,4x,i3,3x,i3,2x,i4,52x,i3)
 14   format(8x,i4,i2,i2,i2,6x,a4,2x,i3,1x,i4,3x,i4,3x,i3,2x,i4,47x,i3,
     &     2x,i3)
 1000 FORMAT('ERROR: The Storm Hindcast/Forecast Input File (unit 22)')
 1021 FORMAT('contains invalid TECH identifier: ',A4,' on PROC=',I4)
 1025 FORMAT('contains invalid TECH identifier: ',A4)
 1031 FORMAT('Valid TECH input is BEST (hindcast) or OFCL (forecast).')
 1041 FORMAT('Check the 5th column of the fort.22 file.')

 1121 FORMAT('contains invalid data for central pressure or Rmax.')
 1131 FORMAT('These values must both be nonzero.')
C
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE GetHollandStormData
C     ----------------------------------------------------------------


C     ----------------------------------------------------------------
C      S U B R O U T I N E   N W S  6 7  G E T
C     ----------------------------------------------------------------
C
C     jgf46.02 Subroutine written by Brian Blanton to calculate wind
C     velocity at nodes from the Holland Wind model. THIS SUBROUTINE IS
C     NOT CURRENTLY CALLED FROM ANYWHERE, AND IS INCLUDED HERE AS
C     REFERENCE MATERIAL. Code may be written in the future to call this
C     subroutine.
C
C     From Brian:
C
C     "This NWS67GET doesn't read a standard track file but rather it
C     reads a fort.22 file that contains along-track parameters.
C
C     time      lon      lat      dp          du       dv     RMW    B
C                                [Pa]        [m/s]    [m/s]    km
C      0.000  -89.600   18.792 11500.000    0.000    8.200    56   1.40
C      6.000  -89.600   20.563 11500.000    0.000    8.200    56   1.40
C     12.000  -89.600   22.334 11500.000    0.000    8.200    56   1.40
C     18.000  -89.600   24.105 11500.000    0.000    8.200    56   1.40
C     24.000  -89.600   25.876 11500.000    0.000    8.200    56   1.40
C     30.000  -89.600   27.648 11500.000    0.000    8.200    56   1.40
C     36.000  -89.600   29.419 11500.000    0.000    8.200    56   1.40
C     42.000  -89.600   31.190 11500.000    0.000    8.200    56   1.40
C
C     time isn't actually used.
C     dp is press diff from ambient
C     du,dv is storm translation vel
C     not sure why RMW is in km, but it made sense at the time."
C
C     ----------------------------------------------------------------
C23456
      SUBROUTINE NWS67GET  (X,Y,SLAM,SFEA,WVNX,WVNY,PRESS,NP,ICS,
     &     RHOWAT0,G,TIMELOC,WTIMINC,WTIME1,WTIME2)
      USE SIZES
      IMPLICIT NONE
      REAL(8), PARAMETER :: Ambient_Pressure=101300.d0      !   Pascals!!
      REAL(8), PARAMETER :: Coriolis=7.287588493541447D-05  ! 1/sec
      REAL(8), PARAMETER :: Rho_Air = 1.15d0                ! kg/ m/m/m
      INTEGER, SAVE :: FIRSTCALL = 0
      INTEGER NP,I,J,ICS
      REAL(SZ) RHOWAT0,G,WTIME1,WTIME2,TIMELOC,WTIMINC
      REAL(SZ) X(*),Y(*),SLAM(*),SFEA(*)
      REAL(SZ) WVNX(*),WVNY(*),PRESS(*)
! LOCAL VARIABLES
      REAL(SZ),SAVE,ALLOCATABLE :: RAD(:),DX(:),DY(:),XCOOR (:),YCOOR(:)
      REAL(SZ),SAVE,ALLOCATABLE :: EARG(:)
      REAL(SZ),SAVE,ALLOCATABLE :: VG(:),THETA(:)
      REAL(SZ),SAVE :: XC1,YC1,DP1,TVX1,TVY1,RMW1,A1,B1
      REAL(SZ),SAVE :: XC2,YC2,DP2,TVX2,TVY2,RMW2,A2,B2
      REAL(SZ) :: XC,YC,DP,TVX,TVY,RMW,A,B,WTRATIO
      REAL(SZ) :: Central_Pressure,T1,T2,T3,D,N,TT,TM
      REAL(SZ),SAVE :: RHOWATG
      IF (FIRSTCALL.EQ.0) THEN
         FIRSTCALL = 1
         RHOWATG=RHOWAT0*G
         ALLOCATE (RAD(NP),DX(NP),DY(NP),XCOOR(NP),YCOOR(NP),EARG (NP))
         ALLOCATE (VG(NP),THETA(NP))
!         IF(ICS.EQ.2) THEN
            DO I=1,NP
            XCOOR(I)=SLAM(I)*RAD2DEG
            YCOOR(I)=SFEA(I)*RAD2DEG
            END DO
!         ENDIF
!         IF(ICS.EQ.1) THEN
!            DO I=1,NP
!               XCOOR(i)=X(i)
!                  YCOOR(i)=Y(i)
!            END DO
!         ENDIF
         READ(22,*) TT,XC1,YC1,DP1,TVX1,TVY1,RMW1,B1
         READ(22,*) TT,XC2,YC2,DP2,TVX2,TVY2,RMW2,B2
         OPEN(667,file='fort.22.adc')
      ENDIF
      IF(TIMELOC.GT.WTIME2) THEN
         WTIME1=WTIME2
         WTIME2=WTIME2+WTIMINC
         XC1=XC2
         YC1=YC2
         DP1=DP2
         TVX1=TVX2
         TVY1=TVY2
         RMW1=RMW2
         B1=B2
         READ(22,*) TT,XC2,YC2,DP2,TVX2,TVY2,RMW2,B2
      END IF
      WTRATIO=(TIMELOC-WTIME1)/WTIMINC
      DP  = DP1  + WTRATIO*(DP2-DP1)
      XC  = XC1  + WTRATIO*(XC2-XC1)
      YC  = YC1  + WTRATIO*(YC2-YC1)
      TVX = TVX1 + WTRATIO*(TVX2-TVX1)
      TVY = TVY1 + WTRATIO*(TVY2-TVY1)
      RMW = RMW1 + WTRATIO*(RMW2-RMW1)
      B   = B1   + WTRATIO*(B2-B1)

      Central_Pressure=Ambient_Pressure-DP
      A=(RMW*1000)**B
      DX=XCOOR-XC
      DY=YCOOR-YC
      THETA=ATAN2(DY,DX)
      RAD=SQRT(DX*DX+DY*DY)*100.d0*1000.d0  ! into meters!!
      EARG=EXP(-A/RAD**B)

      ! WVNX,WVNY,PRESS are returned!!
      DO I=1,NP
         PRESS(I)=(Central_Pressure + (DP)*EARG(I))/RHOWATG
         D=Rho_Air*RAD(I)**B
         N=A*B*DP*EARG(I)
         T1=N/D
         T2=((CORIOLIS**2.d0)/4.d0)*RAD(I)**2.d0
         T3=RAD(I)*CORIOLIS/2.d0
         VG(I)=sqrt(T1+T2)-T3
         TM=(RAD(I)/(RMW*1000))-5.
         IF (TM.LT.0.)THEN
            TM=0.
         ELSEIF (TM.GT.1.)THEN
            TM=1.
         END IF
         TM=-TM/2
         TM=COS(TM*PI)
         WVNX(I)=-VG(I)*SIN(THETA(I))+TVX*TM
         WVNY(I)= VG(I)*COS(THETA(I))+TVY*TM
      END DO

C      write(667,'(12(f12.4,1x))')TIMELOC,XC,YC,TVX,TVY,DP,RMW,B
C     ----------------------------------------------------------------
      END SUBROUTINE NWS67GET
C     ----------------------------------------------------------------



      !=================================================================
      !=================================================================
      !=================================================================
      !      =====                                           =====
      !      =====            SUBROUTINE NWS19GET            =====
      !      =====                                           =====
      !=================================================================
      !=================================================================
      !=================================================================

      !=================================================================
      ! This subroutine takes a wind file in best track format
      ! that has been preprocessed by the ASWIP program and
      ! filled out forecast increments from the intended start of the simulation
      ! creates an asymmetric hurricane vortex, returning a wind field
      ! (wvnx,wvny) and atmospheric pressure-induced water surface
      ! elevation (press) on-the-fly for each model time step.
      !
      ! On input:
      !    slam     Longitude at model nodal points (radians)
      !    sfea     Latitude  at model nodal points (radians)
      !    np       Number of nodal points in model grid
      !    TIMELOC     Model time (seconds)
      !    ics      Coordinate system selection parameter
      !                ics = 1: cartesian coordinates
      !                ics = 2: spherical coordinates
      !                Note: Subroutine valid only for ics = 2
      !
      ! On output:
      !    wvnx     x component of wind velocity at nodal points (m/s)
      !    wvny     y component of wind velocity at nodal points (m/s)
      !    press    Atmospheric pressure-induced water surface
      !             elevation at nodal points (m)
      !
      ! Revision history:
      !    Date         Programmer                 Description of change
      !    ----         ----------                 ---------------------
      !    06/19/06     Cristina Forbes, UNC-CEP   Wrote original code
      !    07/09/06     Craig  Mattocks, UNC-CEP   Code clean up
      !    07/20/06     Jason G. Fleming UNC-IMS   Add into v46.16
      !    07/20/06     Cristina Forbes  (wrote)   v46.16: comment out debug
      !                 Jason G. Fleming (merged)  msg, check for all 4 radii
      !    08/15/06     Cristina Forbes, UNC-CEP   if missing values or radii
      !                                            = 0 -> set winds = 0  and
      !                                            press = P0, comment out
      !                                            fort.300 output
      !    08/25/06     Cristina Forbes, UNC-CEP   Implemented hotstart
      !                                            capability
      !    08/25/06     Craig  Mattocks, UNC-CEP   Applied wind reduction:
      !                                            top of SFC layer --> SFC
      !    09/12/06     Craig  Mattocks, UNC-CEP   Subtracted translational
      !                                            wind speed from Vmax.
      !    07/05/07     Jason G. Fleming           Merged into v47.04, added
      !                                            ScreenUnit
      !    12/12/08     Cristina Forbes, UNC-IMS   Comment out debug statements
      !    04/21/08     Cristina Forbes, UNC-IMS   Added capability to compute
      !                                            winds in all hemispheres
      !    05/12/09     Robert Weaver, UNC-IMS     Modified damping of translational
      !                                            velocity based on ratio of V/Vmax
      !    05/19/09     Cristina Forbes, UNC-IMS   Implemented old bug fix: units conversion
      !    05/22/09     Cristina Forbes, UNC-IMS   Reverted damping of translational
      !                                            velocity back to the original formulation
      !                                            developed by Craig Mattocks
      !    05/26/09    Cristina Forbes, UNC-IMS    Fixed implementation of Rick Luettich
      !                                            V/Vmax tapering formulation (as in NWS=8)
      !                                            for future experimentation - not activated
      !    7/8/09      Robert Weaver               rewrite to accept new input file w/ rmax
      !                                            and Holland B
      !    7/21/09     Robert Weaver               Modify to let user choose which isotach to use
      !                                            when submitting Rmax
      !                                            NWS=19
      !
      !
      !
      ! Reference:
      !
      !     This NWS=19 option is based on the work of
      !     the following publications:
      !
      !     Mattocks, C. and C. Forbes (2008): A real-time, event-triggered
      !     storm surge forecasting system for the State of North Carolina,
      !     Ocean Modelling, 25, pp. 95-119, 10.1016/j.ocemod.2008.06.008
      !
      !     Mattocks, C., C. Forbes and L. Ran (2006) Design and implementation
      !     of a real-time storm surge and flood forecasting capability for the
      !     State of North Carolina. UNC-CEP Technical Report, November 30, 2006, 103 pp.
      !
      !     http://www.unc.edu/~cmattock/ncfs/NCFSreport.pdf
      !=================================================================
      SUBROUTINE NWS19GET(slam,sfea, wvnx,wvny, press, np, TIMELOC, ics)
      USE VORTEX, ONLY : newVortex, setRmaxes, fitRmaxes, uvp,
     &                   getLatestRmax, getLatestAngle
      IMPLICIT NONE

      INTEGER  :: np, ics
      INTEGER  :: ient,itpc
      INTEGER  :: i, nwi
      INTEGER  :: num_entry
      REAL(sz) :: TIMELOC, time1, wtratio
      REAL( 8), SAVE, ALLOCATABLE :: xcoor(:),ycoor(:)
      REAL( 8) :: slam(np),sfea(np)
      REAL( 8) :: wvnx(np),wvny(np), press(np)
      REAL(sz) :: RhoWatG

      REAL,DIMENSION(:),ALLOCATABLE,SAVE::lat,lon
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE::ilat,ilon,atcfRMW
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: ilat1,ilon1,atcfRMW1
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: ispd,icpress
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: ispd1,icpress1
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: iyear,imth,iday,ihr
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE ::iyear1,imth1,iday1,ihr1
      REAL(sz), SAVE ::  lat1, lon1,  spd1,  cpress1,  rmw1
      REAL(sz), SAVE ::  dt, told, cLat1, cLon1
      CHARACTER(LEN=4),DIMENSION(:),ALLOCATABLE, SAVE:: type! hindcast/nowcast or forecast?
      CHARACTER(LEN=4),DIMENSION(:),ALLOCATABLE, SAVE:: type1 !hindcast/nowcast or forecast?
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE::iFcstInc1, iFcstInc ! hours between forecasts
      INTEGER , SAVE :: firstCall = 0
      LOGICAL , SAVE :: firstTime = .TRUE.
      REAL(sz),DIMENSION(:),ALLOCATABLE,SAVE :: HollB,HollB1
      INTEGER ,DIMENSION(:),ALLOCATABLE,SAVE :: advr,advr1
      REAL(sz), SAVE :: Vmax
      REAL(sz), SAVE :: Vr
      REAL(sz), SAVE, DIMENSION(4) :: r
      REAL(sz), SAVE :: Pn
      REAL(sz), SAVE :: Pc
      REAL(sz), SAVE :: cLat
      REAL(sz), SAVE :: cLon
      REAL(sz), SAVE :: cHollB
      REAL(sz) , DIMENSION(4), SAVE :: crmaxw
      REAL(sz), SAVE :: timeOffset
!-------------------------------
      INTEGER,SAVE :: icyc
!------------------------------
      INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: cycle_num
      INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: isotachs_per_cycle
      INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: isotachs_per_cycle1
      INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: ipn,ipn1,cycle_num1

      REAL(sz) , DIMENSION(:,:), ALLOCATABLE,SAVE :: rmaxw,rmaxw1
      REAL(sz) , DIMENSION(:), ALLOCATABLE,SAVE::CycleTime,CycleTime1
      REAL(sz) , DIMENSION(:), ALLOCATABLE,SAVE::uTrans1,uTrans
      REAL(sz) , DIMENSION(:), ALLOCATABLE,SAVE::vTrans1,vTrans
      REAL(sz), DIMENSION(:), ALLOCATABLE, SAVE :: h_speed
      INTEGER , DIMENSION(:,:), ALLOCATABLE,SAVE::quadflag,quadflag1
      INTEGER ,DIMENSION(:),ALLOCATABLE,SAVE :: ivr, dir,speed
      REAL(sz), SAVE :: dx,dy
      REAL(sz) :: VmaxBL
      REAL(sz) :: stormMotionU
      REAL(sz) :: stormMotionV
      REAL(sz) :: stormMotion
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: ivr1, dir1, speed1
      INTEGER , DIMENSION(:,:),ALLOCATABLE, SAVE :: ir,ir1
      CHARACTER(LEN = 10),DIMENSION(:),ALLOCATABLE,SAVE::name,name1
      REAL(sz),SAVE :: timeOld,timeNew
      REAL(sz) :: uTransNow, vTransNow ! time-interpolated overland speed, kts
      CHARACTER(1) ew,ns
      INTEGER, SAVE :: num_cycles ! num unique date/times in the fort.22
C
      call setMessageSource("nws19get")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      timeNew = TIMELOC
      rhoWatG = rhoWat0 * g
      !----------------------------------------
      ! Read (lon,lat) at nodal points and
      ! transform them from radians --> degrees
      !----------------------------------------
      IF (firstTime) THEN
         firstTime = .FALSE.
         ALLOCATE (xcoor(np), ycoor(np))

         DO i=1,np
            xcoor(i) = slam(i) * rad2deg
            ycoor(i) = sfea(i) * rad2deg
         END DO
         timeOffset=TIMELOC
      END IF   ! end of first Time data allocations and such
      !
      ! the mesh changes at each time step if we are writing out
      ! the full circle Rmaxes
      if ((writeFullCircleRmaxes.eqv..true.).or.
     &           (writeRadialVandP.eqv..true.).or.
     &           (writeSpatialVandP.eqv..true.)) then
         xcoor(:) = slam(:) * rad2deg
         ycoor(:) = sfea(:) * rad2deg
      endif
      !---------------
      ! Initialization
      ! it's the first time the routine is called so read
      ! in all variables and set the desired radii and
      ! such
      ! jgf49.10: Replaced opening and closing of fort.22 file with
      ! rewind, since opening is already done in the calling routine;
      ! this also makes it easier to vary the file name (needed
      ! for NWS29).
      !---------------
      IF (firstCall .EQ. 0) THEN  ! first time the routine is called
         i=0
         DO      ! determine the number of entries in the file
            READ(22,*,END=991)
            i=i+1
         ENDDO
991      continue
         num_entry=i
         rewind(22)
         i=num_entry

         allocate(advr1(i),iyear1(i),imth1(i),iday1(i),ihr1(i),type1(i))
         allocate(advr(i),iyear(i),imth(i),iday(i),ihr(i),type(i))
         allocate(iFcstInc1(i),ilat1(i),ilon1(i), ispd1(i))
         allocate(iFcstInc(i) ,ilat(i) ,ilon(i) , ispd(i))
         allocate(lat(i),lon(i))
         allocate(icpress1(i), ivr1(i),ir1(i,4), rmaxw1(i,4), ipn1(i))
         allocate(icpress(i) , ivr(i) ,ir(i,4) , rmaxw(i,4) , ipn(i))
         allocate(dir1(i),speed1(i),cycle_num1(i))
         allocate(dir(i) ,speed(i) ,cycle_num(i))
         allocate(isotachs_per_cycle1(i),CycleTime1(i))
         allocate(isotachs_per_cycle(i) ,CycleTime(i))
         allocate(quadflag1(i,4),name1(i))
         allocate(quadflag(i,4) ,name(i))
         allocate(uTrans1(i), vTrans1(i),HollB1(i),atcfRMW1(i) )
         allocate(uTrans(i) , vTrans(i) ,HollB(i) ,atcfRMW(i) )
         allocate(h_speed(i))

         ! Read parsed NHC advisory in modified best-track format

26       FORMAT(3x, i3, 2x, i4, 3i2, 6x, a4, 2(2x,i3),a1, 3x, i3,a1, 2x,
     &         i3,2x, i4, 6x, i3, 7x, 4(i4,2x), i4, 8x, i3,
     &          27x, 2(i3,2x), a10, 2x, i4, 2x, i4,
     &         2x, 4(i1,2x),2x, 4(f6.1,2x), 2x, f8.4) !,2(2x,f10.4))

         ient=1

         READ(22,26) advr1(ient), iyear1(ient),imth1(ient),iday1(ient),
     &    ihr1(ient), type1(ient),iFcstInc1(ient), ilat1(ient),ns,
     &    ilon1(ient),ew, ispd1(ient), icpress1(ient), ivr1(ient),
     &    (ir1(ient,i),i=1,4), ipn1(ient), atcfRMW1(ient),
     &    dir1(ient),speed1(ient), name1(ient),
     &    cycle_num1(ient),isotachs_per_cycle1(ient),
     &    (quadflag1(ient,i),i=1,4),(rmaxw1(ient,i),i=1,4),HollB1(ient)

         icyc=1

         isotachs_per_cycle(icyc) = isotachs_per_cycle1(ient)
         advr(icyc)  =  advr1(ient)
         iyear(icyc) =  iyear1(ient)
         imth(icyc)  =  imth1(ient)
         iday(icyc)  =  iday1(ient)
         ihr(icyc)   =  ihr1(ient)
         type(icyc)  =  type1(ient)
         iFcstInc(icyc)= iFcstInc1(ient)
         ilat(icyc)  =  ilat1(ient)
         ilon(icyc)  =  ilon1(ient)
         ispd(icyc)  =  ispd1(ient)
         icpress(icyc) =  icpress1(ient)
         ivr(icyc)   =  ivr1(ient)
         ipn(icyc)   =  ipn1(ient)
         atcfRMW(icyc) = atcfRMW1(ient)
         !    dir(icyc)   =  dir1(ient)
         h_speed(icyc) =  real(speed1(ient))
         name(icyc)  =  name1(ient)
         cycle_num(icyc) = cycle_num1(ient)
         HollB(icyc) = HollB1(ient)

         ! switch sin and cos b/c these are compass directions
         uTrans(icyc) = sin(dir1(ient)*deg2rad)*speed1(ient)
         vTrans(icyc) = cos(dir1(ient)*deg2rad)*speed1(ient)

         CycleTime(icyc) =  iFcstInc1(ient) * 3600.d0 + timeOffset
!         The logic for the code is set for W and N so if E or S
!         multiply lats and lons by -1
!
         IF(ew.EQ.'E')THEN
            ilon(icyc)=(-1)*ilon(icyc)
         ENDIF
         IF(ns.EQ.'S')THEN
            ilat(icyc)=(-1)*ilat(icyc)
         ENDIF
         lon(icyc) = (-1.0d0) * ilon(icyc) * 0.1d0
         lat(icyc)=ilat(icyc) * 0.1d0

         ! PICK OUT THE CYCLE DATA THAT THE USER WANTS TO USE
         do i=1,4  ! loop quads
           if (quadflag1(ient,i)==1) then
              rmaxw(icyc,i) = rmaxw1(ient,i)
              quadflag(icyc,i) = quadflag1(ient,i)
           endif
         enddo
! ---------------------------------------------------------------------------
         DO ient=2,num_entry  ! LOOP through entris (lines) in input file

            READ(22,26) advr1(ient), iyear1(ient),imth1(ient),
     &        iday1(ient),
     &        ihr1(ient), type1(ient),iFcstInc1(ient), ilat1(ient),ns,
     &        ilon1(ient),ew, ispd1(ient), icpress1(ient), ivr1(ient),
     &        (ir1(ient,i),i=1,4), ipn1(ient), atcfRMW1(ient),
     &        dir1(ient),speed1(ient), name1(ient),
     &        cycle_num1(ient),isotachs_per_cycle1(ient),
     &        (quadflag1(ient,i),i=1,4),(rmaxw1(ient,i),i=1,4),
     &        HollB1(ient)

            if ( cycle_num1(ient) == cycle_num1(ient-1)) then
               ! do nothing
            else
               IF ( iFcstInc1(ient) == 0 .AND. iFcstInc(icyc) == 0) then
                  icyc = icyc
               ELSE
                  icyc=icyc+1
               ENDIF
               isotachs_per_cycle(icyc) = isotachs_per_cycle1(ient)

               advr(icyc)  =  advr1(ient)
               iyear(icyc) =  iyear1(ient)
               imth(icyc)  =  imth1(ient)
               iday(icyc)  =  iday1(ient)
               ihr(icyc)   =  ihr1(ient)
               type(icyc)  =  type1(ient)
               iFcstInc(icyc)= iFcstInc1(ient)
               ilat(icyc)  =  ilat1(ient)
               ilon(icyc)  =  ilon1(ient)
               ispd(icyc)  =  ispd1(ient)
               icpress(icyc) =  icpress1(ient)
               ivr(icyc)   =  ivr1(ient)
               ipn(icyc)   =  ipn1(ient)
               atcfRMW(icyc) = atcfRMW1(ient)
               dir(icyc)   =  dir1(ient)
               h_speed(icyc) =  speed1(ient)
               name(icyc)  =  name1(ient)
               cycle_num(icyc) = cycle_num1(ient)
               HollB(icyc) = HollB1(ient)

              ! switch sin and cos b/c these are compass directions
              uTrans(icyc)=sin(dir1(ient)*deg2rad)*speed1(ient)
              vTrans(icyc)=cos(dir1(ient)*deg2rad)*speed1(ient)

               CycleTime(icyc) =  iFcstInc(icyc) * 3600.d0 + timeOffset

!              The logic for the code is set for W and N so if E or S
!              multiply lats and lons by -1
               IF(ew.EQ.'E')THEN
                  ilon(icyc)=(-1)*ilon(icyc)
               ENDIF
               IF(ns.EQ.'S')THEN
                  ilat(icyc)=(-1)*ilat(icyc)
               ENDIF
               lon(icyc) = (-1.0d0) * ilon(icyc) * 0.1d0
               lat(icyc) = ilat(icyc) * 0.1d0

            endif

            ! PICK OUT THE CYCLE DATA THAT THE USER WANTS TO USE
            do i=1,4 ! loop through quads
              if (quadflag1(ient,i)==1) then
                rmaxw(icyc,i) = rmaxw1(ient,i)
                quadflag(icyc,i) = quadflag1(ient,i)
              endif
            enddo ! loop through quads

         ENDDO ! loop through entries (lines ) in input file
         write(scratchMessage,12) icyc
 12      format('There are ',i2,' cycles in the fort.22 input file')
         CALL logMessage(INFO,scratchMessage)
         num_cycles = icyc
         icyc=2

      END IF   ! (firstCall .EQ. 0)

      !--------------------------------------
      !
      ! This part of code executes each time step
      !
      !--------------------------------------

C     If time exceeds the next hindcast/nowcast/forecast time, increment the
C     array counter.
C     jgf50.30 Changed .ge. to .gt. so that we wouldn't start looking at the
C     next cycle until we are in it ... this should prevent problems when
C     the ADCIRC run ends right at the end of the fort.22 ... otherwise,
C     we are running off the end of the data.
      IF (TIMELOC.gt.CycleTime(icyc)) THEN
         IF (icyc.gt.num_cycles) THEN
            CALL allMessage(WARNING,
     & "ADCIRC simulation time is later than the last data in fort.22.")
            CALL allMessage(WARNING,
     &      "The simulation has run out of meteorological data.")
         ELSE
            icyc=icyc+1
         ENDIF
      ENDIF

      ! Interpolate NHC forecast interval in
      ! time to obtain values at model time.
      wtratio=(TIMELOC-CycleTime(icyc-1))/
     &                  (CycleTime(icyc)-CycleTime(icyc-1))
      !
      ! Perform time interpolation, transform variables from integers
      ! to real numbers for hurricane vortex calcualtions.
      Vmax  =  1.d0*(ispd(icyc-1) +
     &                         wtratio * (ispd(icyc)-ispd(icyc-1)))
      Pn    =  1.d0*(ipn(icyc-1)  +
     &                         wtratio * (ipn(icyc)-ipn(icyc-1)))
      Pc    =  1.d0*(icpress(icyc-1) +
     &                         wtratio*(icpress(icyc)-icpress(icyc-1)))
      cLat   =  lat(icyc-1)  +
     &                         wtratio * (lat(icyc)-lat(icyc-1))
      cLon  =  lon(icyc-1)  +
     &                         wtratio * (lon(icyc)-lon(icyc-1))
      cHollB   =  HollB(icyc-1) +
     &                         wtratio * (HollB(icyc)-HollB(icyc-1))
      do i=1,4
         crmaxw(i)=rmaxw(icyc-1,i) +
     &                          wtratio*(rmaxw(icyc,i)-rmaxw(icyc-1,i))
      enddo
      uTransNow = uTrans(icyc-1) + wtratio
     &                  * (uTrans(icyc)-utrans(icyc-1))
      vTransNow = vTrans(icyc-1) + wtratio
     &                  * (vTrans(icyc)-vTrans(icyc-1))

      ! jgf49.0803 Make the current Rmax and center location
      ! available at the module level for use in NWS29
      vortexLat = cLat
      vortexLon = cLon
      vortexRMW = 0.0d0
      !jgf50.32: For sector-based wind drag, record the center
      !of the storm.
      IF ((cLat.ne.EyeLatR(3)).or.(cLon.ne.EyeLonR(3))) THEN
         EyeLatR(1) = EyeLatR(2)
         EyeLonR(1) = EyeLonR(2)
         EyeLatR(2) = EyeLatR(3)
         EyeLonR(2) = EyeLonR(3)
         EyeLatR(3) = cLat
         EyeLonR(3) = cLon
         FoundEye = .true.
      ENDIF

      ! arithmetic mean of the Rmax in the 4 storm quadrants
      do i=1,4
         vortexRMW = vortexRMW + 0.25d0*crmaxw(i)
      end do

      !-------------------------------
      ! Create a new asymmetric hurricane vortex.
      !
      ! Note: Subtract translational speed from Vmax, then
      ! scale (Vmax - Vt) and Vr up to the top of the surface,
      ! where the cylcostrophic wind balance is valid.
      !-------------------------------------------------------
      stormMotion = 1.5d0*(SQRT(uTransNow**2.d0+vTransNow**2.d0))**0.63d0
      VmaxBL = (Vmax-stormMotion)/windReduction
      call newVortex(Pn,Pc,cLat,cLon,VmaxBL)
      call setRmaxes(crmaxw)
      call fitRmaxes()

      !-------------------------------------------------------------
      ! Calculate wind and pressure fields at model nodal points.
      !
      ! Note: the asymmetric vortex wind speed is reduced from the
      ! top of the surface layer to the surface, then converted from
      ! a 1-minute (max sustained) to a 10-minute average prior to
      ! adding the translational velocity in subroutine uvp.
      !-------------------------------------------------------------
      stormMotionU=sin(dir(icyc)/rad2deg)*stormMotion
      stormMotionV=cos(dir(icyc)/rad2deg)*stormMotion
      DO i=1,np
         CALL uvp(ycoor(i),xcoor(i),stormMotionU,stormMotionV,
     &        wvnx(i),wvny(i), press(i)) !
         !-------------------------------------------
         ! Convert atmospheric pressure (Pascals) to
         ! atmospheric pressure-induced water surface
         ! elevation (meters).
         !-------------------------------------------
         press(i) = press(i) / RhoWatG

         if (writeFullCircleRmaxes.eqv..true.) then
            write(444,*) getLatestRmax(), getLatestAngle()
         endif

      END DO
      !-----------------------------------
      ! NHC advisory best-track i/o format
      !-----------------------------------
      firstCall = 1
      told=TIMELOC

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
999   RETURN
C***********************************************************************
      END SUBROUTINE NWS19GET
C***********************************************************************




      !=================================================================
      !=================================================================
      !=================================================================
      !      =====                                           =====
      !      =====            SUBROUTINE NWS20GET            =====
      !      =====                                           =====
      !=================================================================
      !=================================================================
      !=================================================================

      !=================================================================
      ! This subroutine takes a wind file in best track format
      ! that has been preprocessed by the ASWIP program and
      ! filled out forecast increments from the intended start of the simulation
      ! creates an asymmetric hurricane vortex, returning a wind field
      ! (wvnx,wvny) and atmospheric pressure-induced water surface
      ! elevation (press) on-the-fly for each model time step.
      !
      ! On input:
      !    slam     Longitude at model nodal points (radians)
      !    sfea     Latitude  at model nodal points (radians)
      !    np       Number of nodal points in model grid
      !    TIMELOC     Model time (seconds)
      !    ics      Coordinate system selection parameter
      !                ics = 1: cartesian coordinates
      !                ics = 2: spherical coordinates
      !                Note: Subroutine valid only for ics = 2
      !
      ! On output:
      !    wvnx     x component of wind velocity at nodal points (m/s)
      !    wvny     y component of wind velocity at nodal points (m/s)
      !    press    Atmospheric pressure-induced water surface
      !             elevation at nodal points (m)
      !
      ! Revision history:
      !    Date         Programmer                 Description of change
      !    ----         ----------                 ---------------------
      !    06/19/06     Cristina Forbes, UNC-CEP   Wrote original code
      !    07/09/06     Craig  Mattocks, UNC-CEP   Code clean up
      !    07/20/06     Jason G. Fleming UNC-IMS   Add into v46.16
      !    07/20/06     Cristina Forbes  (wrote)   v46.16: comment out debug
      !                 Jason G. Fleming (merged)  msg, check for all 4 radii
      !    08/15/06     Cristina Forbes, UNC-CEP   if missing values or radii
      !                                            = 0 -> set winds = 0  and
      !                                            press = P0, comment out
      !                                            fort.300 output
      !    08/25/06     Cristina Forbes, UNC-CEP   Implemented hotstart
      !                                            capability
      !    08/25/06     Craig  Mattocks, UNC-CEP   Applied wind reduction:
      !                                            top of SFC layer --> SFC
      !    09/12/06     Craig  Mattocks, UNC-CEP   Subtracted translational
      !                                            wind speed from Vmax.
      !    07/05/07     Jason G. Fleming           Merged into v47.04, added
      !                                            ScreenUnit
      !    12/12/08     Cristina Forbes, UNC-IMS   Comment out debug statements
      !    04/21/08     Cristina Forbes, UNC-IMS   Added capability to compute
      !                                            winds in all hemispheres
      !    05/12/09     Robert Weaver, UNC-IMS     Modified damping of translational
      !                                            velocity based on ratio of V/Vmax
      !    05/19/09     Cristina Forbes, UNC-IMS   Implemented old bug fix: units conversion
      !    05/22/09     Cristina Forbes, UNC-IMS   Reverted damping of translational
      !                                            velocity back to the original formulation
      !                                            developed by Craig Mattocks
      !    05/26/09    Cristina Forbes, UNC-IMS    Fixed implementation of Rick Luettich
      !                                            V/Vmax tapering formulation (as in NWS=8)
      !                                            for future experimentation - not activated
      !    7/8/09      Robert Weaver               rewrite to accept new input file w/ rmax
      !                                            and Holland B
      !    7/21/09     Robert Weaver               Modify to let user choose which isotach to use
      !                                            when submitting Rmax
      !                                            NWS=19
      !    5/20/13     Jie Gao, UNC-IMS            Added nws=20 option to use quadrant-varying Holland B, 
      !                                            and PhiFactor, and occasionally 
      !                                            Vpseudo=Vr if violations occur when taking out 
      !                                            translational speed; 
      !                                            use full gradient wind equation
      !                                            introduce a linear-weighted combination method to use
      !                                            multiple isotachs 
      !                                             
      !
      !
      ! Reference:
      !
      !     This NWS=19 option is based on the work of
      !     the following publications:
      !
      !     Mattocks, C. and C. Forbes (2008): A real-time, event-triggered
      !     storm surge forecasting system for the State of North Carolina,
      !     Ocean Modelling, 25, pp. 95-119, 10.1016/j.ocemod.2008.06.008
      !
      !     Mattocks, C., C. Forbes and L. Ran (2006) Design and implementation
      !     of a real-time storm surge and flood forecasting capability for the
      !     State of North Carolina. UNC-CEP Technical Report, November 30, 2006, 103 pp.
      !
      !     http://www.unc.edu/~cmattock/ncfs/NCFSreport.pdf
      !=================================================================
      SUBROUTINE NWS20GET(slam,sfea, wvnx,wvny, press, np, TIMELOC, ics)
      USE VORTEX, ONLY : setVortex, fitRmaxes4,  uvpr, spInterp,
     &                   Rmaxes4, QuadFlag4, QuadIr4, Bs4, VmBL4  ! jgfdebug
      USE GLOBAL, ONLY : kt2ms, nm2m,deg2rad,omega, Rearth, m2nm
      IMPLICIT NONE

      INTEGER  :: np, ics
      INTEGER  :: ient,itpc
      INTEGER  :: i, nwi,j
      INTEGER  :: num_entry
      REAL(sz) :: TIMELOC, time1, wtratio
      REAL( 8), SAVE, ALLOCATABLE :: xcoor(:),ycoor(:)
      REAL( 8) :: slam(np),sfea(np)
      REAL( 8) :: wvnx(np),wvny(np), press(np)
      REAL(sz) :: RhoWatG
      REAL(sz) :: corio                       ! Coriolis force (1/s)

      REAL,DIMENSION(:),ALLOCATABLE,SAVE::lat,lon
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE::ilat,ilon,atcfRMW
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: ilat1,ilon1,atcfRMW1
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: ispd,icpress
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: ispd1,icpress1
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: iyear,imth,iday,ihr
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE ::iyear1,imth1,iday1,ihr1
      REAL(sz), SAVE ::  lat1, lon1,  spd1,  cpress1,  rmw1
      REAL(sz), SAVE ::  dt, told, cLat1, cLon1
      CHARACTER(LEN=4),DIMENSION(:),ALLOCATABLE, SAVE:: type! hindcast/nowcast or forecast?
      CHARACTER(LEN=4),DIMENSION(:),ALLOCATABLE, SAVE:: type1 !hindcast/nowcast or forecast?
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE::iFcstInc1, iFcstInc ! hours between forecasts
      INTEGER , SAVE :: firstCall = 0
      LOGICAL , SAVE :: firstTime = .TRUE.
      REAL(sz),DIMENSION(:),ALLOCATABLE,SAVE :: HollB,HollB1
      REAL(sz),DIMENSION(:,:),ALLOCATABLE,SAVE :: VmaxesBL1,HollBs1
      REAL(sz),DIMENSION(:,:,:),ALLOCATABLE,SAVE :: HollBs,VmaxesBL
      INTEGER ,DIMENSION(:),ALLOCATABLE,SAVE :: advr,advr1
      REAL(sz), SAVE :: Vmax
      REAL(sz), SAVE :: Vr
      REAL(sz), SAVE, DIMENSION(4) :: r
      REAL(sz), SAVE :: Pn
      REAL(sz), SAVE :: Pc
      REAL(sz), SAVE :: cLat
      REAL(sz), SAVE :: cLon
      REAL(sz), SAVE :: cHollB
      ! save parameters for each node 
      ! change them to allocatable, save Aug. 2013, Jie
      REAL(sz), ALLOCATABLE, SAVE :: cHollBs1(:), cHollBs2(:)
      REAL(sz), ALLOCATABLE, SAVE :: cPhiFactor1(:), cPhiFactor2(:)
      REAL(sz), ALLOCATABLE, SAVE :: cVmwBL1(:), cVmwBL2(:)
      REAL(sz), ALLOCATABLE, SAVE :: crmaxw1(:), crmaxw2(:)
      REAL(sz), ALLOCATABLE, SAVE :: crmaxwTrue1(:), crmaxwTrue2(:)
      REAL(sz), ALLOCATABLE, SAVE :: cHollBs(:)
      REAL(sz), ALLOCATABLE, SAVE :: cPhiFactor(:)
      REAL(sz), ALLOCATABLE, SAVE :: cVmwBL(:)
      REAL(sz), ALLOCATABLE, SAVE :: crmaxw(:)
      REAL(sz), ALLOCATABLE, SAVE :: crmaxwTrue(:)
      REAL(sz), ALLOCATABLE, SAVE :: aziangle(:), dist(:)
      REAL(sz), SAVE :: timeOffset
!-------------------------------
      INTEGER,SAVE :: icyc
      INTEGER,SAVE :: isot
!------------------------------
      INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: cycle_num
      INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: isotachs_per_cycle
      INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: isotachs_per_cycle1
      INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: ipn,ipn1,cycle_num1

      REAL(sz) , DIMENSION(:,:), ALLOCATABLE,SAVE :: rmaxw1
      REAL(sz) , DIMENSION(:,:,:), ALLOCATABLE,SAVE :: rmaxw
      REAL(sz) , DIMENSION(:), ALLOCATABLE,SAVE::CycleTime,CycleTime1
      REAL(sz) , DIMENSION(:), ALLOCATABLE,SAVE::uTrans1,uTrans
      REAL(sz) , DIMENSION(:), ALLOCATABLE,SAVE::vTrans1,vTrans
      REAL(sz), DIMENSION(:), ALLOCATABLE, SAVE :: h_speed
      INTEGER , DIMENSION(:,:), ALLOCATABLE,SAVE::quadflag1
      INTEGER , DIMENSION(:,:,:), ALLOCATABLE,SAVE::quadflag
      INTEGER ,DIMENSION(:),ALLOCATABLE,SAVE :: ivr, dir,speed
      REAL(sz), SAVE :: dx,dy
      REAL(sz) :: VmaxBL
      REAL(sz) :: stormMotionU
      REAL(sz) :: stormMotionV
      REAL(sz) :: stormMotion
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: ivr1, dir1, speed1
      INTEGER , DIMENSION(:,:),ALLOCATABLE, SAVE :: ir1
      INTEGER , DIMENSION(:,:,:),ALLOCATABLE, SAVE :: ir
      CHARACTER(LEN = 10),DIMENSION(:),ALLOCATABLE,SAVE::name,name1
      REAL(sz),SAVE :: timeOld,timeNew
      REAL(sz) :: uTransNow, vTransNow ! time-interpolated overland speed, kts
      REAL(sz) :: dirNow ! Jie
      CHARACTER(1) ew,ns
      INTEGER, SAVE :: num_cycles ! num unique date/times in the fort.22
      
      call setMessageSource("nws20get")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      timeNew = TIMELOC
      rhoWatG = rhoWat0 * g
      !----------------------------------------
      ! Read (lon,lat) at nodal points and
      ! transform them from radians --> degrees
      !----------------------------------------
      IF (firstTime) THEN
         firstTime = .FALSE.
         ALLOCATE (xcoor(np), ycoor(np))
         ALLOCATE (cHollBs1(np), cHollBs2(np))
         ALLOCATE (cPhiFactor1(np), cPhiFactor2(np))
         ALLOCATE (cVmwBL1(np), cVmwBL2(np))
         ALLOCATE (crmaxw1(np), crmaxw2(np))
         ALLOCATE (crmaxwTrue1(np), crmaxwTrue2(np))
         ALLOCATE (cHollBs(np))
         ALLOCATE (cPhiFactor(np))
         ALLOCATE (cVmwBL(np))
         ALLOCATE (crmaxw(np))
         ALLOCATE (crmaxwTrue(np))
         ALLOCATE (aziangle(np), dist(np))

         DO i=1,np
            xcoor(i) = slam(i) * rad2deg
            ycoor(i) = sfea(i) * rad2deg
         END DO
         timeOffset=TIMELOC
      END IF   ! end of first Time data allocations and such
      !
      ! the mesh changes at each time step if we are writing out
      ! the full circle Rmaxes
      if ((writeFullCircleRmaxes.eqv..true.).or.
     &           (writeRadialVandP.eqv..true.).or.
     &           (writeSpatialVandP.eqv..true.).or.
     &           (writeSpatialuvp.eqv..true.)) then
         xcoor(:) = slam(:) * rad2deg
         ycoor(:) = sfea(:) * rad2deg
      endif
      !---------------
      ! Initialization
      ! it's the first time the routine is called so read
      ! in all variables and set the desired radii and
      ! such
      ! jgf49.10: Replaced opening and closing of fort.22 file with
      ! rewind, since opening is already done in the calling routine;
      ! this also makes it easier to vary the file name (needed
      ! for NWS29).
      !---------------
      IF (firstCall .EQ. 0) THEN  ! first time the routine is called
         i=0
         DO      ! determine the number of entries in the file
            READ(22,*,END=992)
            i=i+1
         ENDDO
992      continue
         num_entry=i
         rewind(22)
         i=num_entry

         allocate(advr1(i),iyear1(i),imth1(i),iday1(i),ihr1(i),type1(i))
         allocate(advr(i),iyear(i),imth(i),iday(i),ihr(i),type(i))
         allocate(iFcstInc1(i),ilat1(i),ilon1(i), ispd1(i))
         allocate(iFcstInc(i) ,ilat(i) ,ilon(i) , ispd(i))
         allocate(lat(i),lon(i))
         allocate(icpress1(i), ivr1(i),ir1(i,4), rmaxw1(i,4), ipn1(i))
         allocate(icpress(i), ivr(i) ,ir(i,4,4) , rmaxw(i,4,4) , ipn(i))
         allocate(dir1(i),speed1(i),cycle_num1(i))
         allocate(dir(i) ,speed(i) ,cycle_num(i))
         allocate(isotachs_per_cycle1(i),CycleTime1(i))
         allocate(isotachs_per_cycle(i) ,CycleTime(i))
         allocate(quadflag1(i,4),name1(i))
         allocate(quadflag(i,4,4),name(i))
         allocate(uTrans1(i), vTrans1(i),HollB1(i),atcfRMW1(i) )
         allocate(uTrans(i) , vTrans(i) ,HollB(i) ,atcfRMW(i) )
         allocate(HollBs1(i,4)) 
         allocate(HollBs(i,4,4)) 
         allocate(h_speed(i))
         allocate(VmaxesBL1(i,4), VmaxesBL(i,4,4))

         ! jgf: initialize arrays to zero
         ir = 0.0d0
         rmaxw = 0.0d0
         quadflag = 0 
         HollBs = 0.0d0
         VmaxesBL = 0.0d0
         
         ! Read parsed NHC advisory in modified best-track format

27       FORMAT(3x, i3, 2x, i4, 3i2, 6x, a4, 2(2x,i3),a1, 3x, i3,a1, 2x,
     &         i3,2x, i4, 6x, i3, 7x, 4(i4,2x), i4, 8x, i3,
     &          27x, 2(i3,2x), a10, 2x, i4, 2x, i4,
     &         2x, 4(i1,2x),2x, 4(f6.1,2x),1x, 8(f8.4,2x), f8.4)

         ient = 1

         READ(22,27) advr1(ient), iyear1(ient),imth1(ient),iday1(ient),
     &    ihr1(ient), type1(ient),iFcstInc1(ient), ilat1(ient),ns,
     &    ilon1(ient),ew, ispd1(ient), icpress1(ient), ivr1(ient),
     &    (ir1(ient,i),i=1,4), ipn1(ient), atcfRMW1(ient),
     &    dir1(ient),speed1(ient), name1(ient),
     &    cycle_num1(ient),isotachs_per_cycle1(ient),
     &    (quadflag1(ient,i),i=1,4),(rmaxw1(ient,i),i=1,4),HollB1(ient),
     &    (HollBs1(ient,i),i=1,4),(VmaxesBL1(ient,i),i=1,4)

         icyc=1
         isot=1
         isotachs_per_cycle(icyc) = isotachs_per_cycle1(ient)
         advr(icyc)  =  advr1(ient)
         iyear(icyc) =  iyear1(ient)
         imth(icyc)  =  imth1(ient)
         iday(icyc)  =  iday1(ient)
         ihr(icyc)   =  ihr1(ient)
         type(icyc)  =  type1(ient)
         iFcstInc(icyc)= iFcstInc1(ient)
         ilat(icyc)  =  ilat1(ient)
         ilon(icyc)  =  ilon1(ient)
         ispd(icyc)  =  ispd1(ient)
         icpress(icyc) =  icpress1(ient)
         ivr(icyc)   =  ivr1(ient)
         ipn(icyc)   =  ipn1(ient)
         atcfRMW(icyc) = atcfRMW1(ient)
         !    dir(icyc)   =  dir1(ient)
         h_speed(icyc) =  real(speed1(ient))
         name(icyc)  =  name1(ient)
         cycle_num(icyc) = cycle_num1(ient)
         HollB(icyc) = HollB1(ient)

         ! switch sin and cos b/c these are compass directions
         uTrans(icyc) = sin(dir1(ient)*deg2rad)*speed1(ient)
         vTrans(icyc) = cos(dir1(ient)*deg2rad)*speed1(ient)

         CycleTime(icyc) =  iFcstInc1(ient) * 3600.d0 + timeOffset
!         The logic for the code is set for W and N so if E or S
!         multiply lats and lons by -1
!
         IF(ew.EQ.'E')THEN
            ilon(icyc)=(-1)*ilon(icyc)
         ENDIF
         IF(ns.EQ.'S')THEN
            ilat(icyc)=(-1)*ilat(icyc)
         ENDIF
         lon(icyc) = (-1.0d0) * ilon(icyc) * 0.1d0
         lat(icyc)=ilat(icyc) * 0.1d0

         ! PICK OUT THE CYCLE DATA THAT THE USER WANTS TO USE
         do i=1,4  ! loop quads
           if (quadflag1(ient,i)==1) then
              ir(icyc,i, isot) = ir1(ient,i)
              rmaxw(icyc,i,isot) = rmaxw1(ient,i)
              quadflag(icyc,i,isot) = quadflag1(ient,i)
              HollBs(icyc,i,isot) = HollBs1(ient,i)
              VmaxesBL(icyc,i,isot) = VmaxesBL1(ient,i)
           endif
         enddo
         
! ---------------------------------------------------------------------------
         DO ient=2,num_entry  ! LOOP through entris (lines) in input file

         READ(22,27) advr1(ient), iyear1(ient),imth1(ient),iday1(ient),
     &    ihr1(ient), type1(ient),iFcstInc1(ient), ilat1(ient),ns,
     &    ilon1(ient),ew, ispd1(ient), icpress1(ient), ivr1(ient),
     &    (ir1(ient,i),i=1,4), ipn1(ient), atcfRMW1(ient),
     &    dir1(ient),speed1(ient), name1(ient),
     &    cycle_num1(ient),isotachs_per_cycle1(ient),
     &    (quadflag1(ient,i),i=1,4),(rmaxw1(ient,i),i=1,4),HollB1(ient),
     &    (HollBs1(ient,i),i=1,4),(VmaxesBL1(ient,i),i=1,4)
     
            if ( cycle_num1(ient) == cycle_num1(ient-1)) then
                isot=isot+1 ! same icyc, next isotach   
            else
                isot=1 ! initialize isotach #
               IF ( iFcstInc1(ient) == 0 .AND. iFcstInc(icyc) == 0) then
                  icyc = icyc
               ELSE
                  icyc = icyc+1
               ENDIF
               isotachs_per_cycle(icyc) = isotachs_per_cycle1(ient)

               advr(icyc)  =  advr1(ient)
               iyear(icyc) =  iyear1(ient)
               imth(icyc)  =  imth1(ient)
               iday(icyc)  =  iday1(ient)
               ihr(icyc)   =  ihr1(ient)
               type(icyc)  =  type1(ient)
               iFcstInc(icyc)= iFcstInc1(ient)
               ilat(icyc)  =  ilat1(ient)
               ilon(icyc)  =  ilon1(ient)
               ispd(icyc)  =  ispd1(ient)
               icpress(icyc) =  icpress1(ient)
               ivr(icyc)   =  ivr1(ient)
               ipn(icyc)   =  ipn1(ient)
               atcfRMW(icyc) = atcfRMW1(ient)
               dir(icyc)   =  dir1(ient)
               h_speed(icyc) =  speed1(ient)
               name(icyc)  =  name1(ient)
               cycle_num(icyc) = cycle_num1(ient)
               HollB(icyc) = HollB1(ient)           

              ! switch sin and cos b/c these are compass directions
              uTrans(icyc)=sin(dir1(ient)*deg2rad)*speed1(ient)
              vTrans(icyc)=cos(dir1(ient)*deg2rad)*speed1(ient)

               CycleTime(icyc) =  iFcstInc(icyc) * 3600.d0 + timeOffset

!              The logic for the code is set for W and N so if E or S
!              multiply lats and lons by -1
               IF(ew.EQ.'E')THEN
                  ilon(icyc)=(-1)*ilon(icyc)
               ENDIF
               IF(ns.EQ.'S')THEN
                  ilat(icyc)=(-1)*ilat(icyc)
               ENDIF
               lon(icyc) = (-1.0d0) * ilon(icyc) * 0.1d0
               lat(icyc) = ilat(icyc) * 0.1d0

            endif

            ! PICK OUT THE CYCLE DATA THAT THE USER WANTS TO USE
            do i=1,4 ! loop through quads
              if (quadflag1(ient,i)==1) then
                ir(icyc,i, isot) = ir1(ient,i)
                rmaxw(icyc,i,isot) = rmaxw1(ient,i)
                quadflag(icyc,i,isot) = quadflag1(ient,i)
                HollBs(icyc,i,isot) = HollBs1(ient,i)
                VmaxesBL(icyc,i,isot) = VmaxesBL1(ient,i)
              endif
            enddo ! loop through quads

         ENDDO ! loop through entries (lines ) in input file
         write(scratchMessage,13) icyc
 13      format('There are ',i2,' cycles in the fort.22 input file')
         CALL logMessage(INFO,scratchMessage)
         num_cycles = icyc
         icyc=2
      END IF   ! (firstCall .EQ. 0)

      !--------------------------------------
      !
      ! This part of code executes each time step
      !
      !--------------------------------------

C     If time exceeds the next hindcast/nowcast/forecast time, increment the
C     array counter.
C     jgf50.30 Changed .ge. to .gt. so that we wouldn't start looking at the
C     next cycle until we are in it ... this should prevent problems when
C     the ADCIRC run ends right at the end of the fort.22 ... otherwise,
C     we are running off the end of the data.
      
      IF (TIMELOC.gt.CycleTime(icyc)) THEN
         IF (icyc.gt.num_cycles) THEN
            CALL allMessage(WARNING,
     & "ADCIRC simulation time is later than the last data in fort.22.")
            CALL allMessage(WARNING,
     &      "The simulation has run out of meteorological data.")
         ELSE
         icyc=icyc+1
         ENDIF
      ENDIF

      ! Interpolate NHC forecast interval in
      ! time to obtain values at model time.
      ! obtain the latest hurricane center location
      wtratio=(TIMELOC-CycleTime(icyc-1))/
     &                  (CycleTime(icyc)-CycleTime(icyc-1))
      cLat   = lat(icyc-1)  +
     &                         wtratio * (lat(icyc)-lat(icyc-1))
      cLon  =  lon(icyc-1)  +
     &                         wtratio * (lon(icyc)-lon(icyc-1))
      ! locate each node based on the latest hurricane center location   
         do i=1,np
              dx = deg2rad * Rearth * (xcoor(i) - cLon)*COS(deg2rad*cLat)
              dy = deg2rad * Rearth * (ycoor(i) - cLat)
              dist(i) = SQRT(dx*dx + dy*dy)* m2nm
              aziangle(i) = 360.0d0 + rad2deg * ATAN2(dx,dy)
            IF (aziangle(i) > 360.d0) aziangle(i) = aziangle(i) - 360.d0
         enddo 
         
 ! set parameters at icyc-1
 !        call setQuadFlag4(reshape(quadflag(icyc-1,:,:),(/4,4/)))
 !        call setQuadIr4(reshape(ir(icyc-1,:,:),(/4,4/)))
 !        call setRmaxes4(reshape(rmaxw(icyc-1,:,:),(/4,4/)))
 !        call setShapeParameters4(reshape(HollBs(icyc-1,:,:),(/4,4/)))
 !        call setVmaxesBL4(reshape(VmaxesBL(icyc-1,:,:),(/4,4/)))
 
         QuadFlag4(2:5,1:4) = quadflag(icyc-1,1:4,1:4)
         QuadIr4(2:5,1:4) = real(ir(icyc-1,1:4,1:4))
         Rmaxes4(2:5,1:4) = rmaxw(icyc-1,1:4,1:4)
         Bs4(2:5,1:4) = HollBs(icyc-1,1:4,1:4)
         VmBL4(2:5,1:4) = VmaxesBL(icyc-1,1:4,1:4)
         
         call fitRmaxes4()
!699      FORMAT(5(f8.4,2x)) 
!700      FORMAT(a25, 2x, 3(f12.2,2x))
 !Jie debug GAHM         
 !        WRITE(16,*) 'Previous: Spatially-interpolated storm parameters',
 !    &                   ' at icyc-1= ', icyc-1, 
 !    &                   ' CycleTime=', CycleTime(icyc-1) 
         do i=1,np
             crmaxw1(i) = spInterp(aziangle(i),dist(i),1)
             crmaxwTrue1(i) = spInterp(aziangle(i),1.d0,1) 
             !an artificial number 1.0 is chosen to ensure only
             !rmax from the highest isotach is picked
             cHollBs1(i)= spInterp(aziangle(i),dist(i),2)
             cVmwBL1(i) = spInterp(aziangle(i),dist(i),3)
 !Jie debug GAHM       
 !            WRITE(16,699) crmaxw1(i),  crmaxwTrue1(i),
 !    &                   cHollBs1(i) ,cVmwBL1(i)
         end do
         
         
         ! set parameters at icyc
!         call setVortex(ipn(icyc),icpress(icyc),cLat,cLon
!     $                  ,isotachs_per_cycle(icyc))
!         call setQuadFlag4(reshape(quadflag(icyc,:,:),(/4,4/)))
!         call setQuadIr4(reshape(ir(icyc,:,:),(/4,4/)))
!         call setRmaxes4(reshape(rmaxw(icyc,:,:),(/4,4/)))
!         call setShapeParameters4(reshape(HollBs(icyc,:,:),(/4,4/)))
!         call setVmaxesBL4(reshape(VmaxesBL(icyc,:,:),(/4,4/)))  
      
         QuadFlag4(2:5,1:4) = quadflag(icyc,1:4,1:4)
         QuadIr4(2:5,1:4) = real(ir(icyc,1:4,1:4))
         Rmaxes4(2:5,1:4) = rmaxw(icyc,1:4,1:4)
         Bs4(2:5,1:4) = HollBs(icyc,1:4,1:4)
         VmBL4(2:5,1:4) = VmaxesBL(icyc,1:4,1:4)
         call fitRmaxes4()
 !Jie debug GAHM         
 !        WRITE(16,*) 'After: Spatially-interpolated storm parameters',
 !    &                   ' at icyc= ', icyc, 
 !    &                   ' CycleTime=', CycleTime(icyc)          
         do i=1,np
             crmaxw2(i) = spInterp(aziangle(i),dist(i),1)
             crmaxwTrue2(i) = spInterp(aziangle(i),1.d0,1) 
             !an artificial number 1.0 is chosen to ensure only
             !rmax from the highest isotach is picked             
             cHollBs2(i)= spInterp(aziangle(i),dist(i),2)
             cVmwBL2(i) = spInterp(aziangle(i),dist(i),3)
              !Jie debug GAHM       
  !           WRITE(16,699) crmaxw2(i),  crmaxwTrue2(i),
  !   &                   cHollBs2(i) ,cVmwBL2(i)
         end do
      ! Perform time interpolation, transform variables from integers
      ! to real numbers for hurricane vortex calcualtions.

      Pn    =  1.d0*(ipn(icyc-1)  +
     &                         wtratio*(ipn(icyc)-ipn(icyc-1)))
      Pc    =  1.d0*(icpress(icyc-1) +
     &                         wtratio*(icpress(icyc)-icpress(icyc-1)))
      crmaxw   =  crmaxw1(:) +
     &                         wtratio*(crmaxw2(:)-crmaxw1(:))
      crmaxwTrue   =  crmaxwTrue1(:) +
     &                         wtratio*(crmaxwTrue2(:)-crmaxwTrue1(:))
      cHollBs  =  cHollBs1(:) +
     &                         wtratio*(cHollBs2(:)-cHollBs1(:))
      cVmwBL   =  cVmwBL1(:) +
     &                         wtratio*(cVmwBL2(:)-cVmwBL1(:))
      corio = 2.0d0 * omega * SIN(deg2rad*cLat)
      
!Jie debug GAHM         
!         WRITE(16,700) 'Temporal interpolation:',
!     &           CycleTime(icyc-1), TIMELOC, CycleTime(icyc) 
     
      do i=1,np
      cPhiFactor(i) =  1 + 
     &                       cVmwBL(i)*kt2ms*crmaxw(i)*nm2m*corio/
     &                       (cHollBs(i)*((cVmwBL(i)*kt2ms)**2+
     &                       cVmwBL(i)*kt2ms*crmaxw(i)*nm2m*corio))
!Jie debug GAHM       
!             WRITE(16,699) crmaxw(i),  crmaxwTrue(i),
!     &                   cHollBs(i) ,cVmwBL(i), cPhiFactor(i)
      enddo
      

     
     
      uTransNow = uTrans(icyc-1) + wtratio
     &                  * (uTrans(icyc)-utrans(icyc-1))
      vTransNow = vTrans(icyc-1) + wtratio
     &                  * (vTrans(icyc)-vTrans(icyc-1))
     
      if (writeFullCircleRmaxes.eqv..true.) then
            do i=1,np
              write(444,*)  crmaxw(i), aziangle(i),dist(i)
            end do
      endif     
      ! jgf49.0803 Make the current Rmax and center location
      ! available at the module level for use in NWS29
      vortexLat = cLat
      vortexLon = cLon
      vortexRMW = 0.0d0
      !jgf50.32: For sector-based wind drag, record the center
      !of the storm.
      IF ((cLat.ne.EyeLatR(3)).or.(cLon.ne.EyeLonR(3))) THEN
         EyeLatR(1) = EyeLatR(2)
         EyeLonR(1) = EyeLonR(2)
         EyeLatR(2) = EyeLatR(3)
         EyeLonR(2) = EyeLonR(3)
         EyeLatR(3) = cLat
         EyeLonR(3) = cLon
         FoundEye = .true.
      ENDIF

      ! arithmetic mean of the Rmax in the 4 storm quadrants
      !do i=1,4
      !   vortexRMW = vortexRMW + 0.25d0*crmaxw(i)
      !end do

      !-------------------------------
      ! Create a new asymmetric hurricane vortex.
      !
      ! Note: Subtract translational speed from Vmax, then
      ! scale (Vmax - Vt) and Vr up to the top of the surface,
      ! where the cylcostrophic wind balance is valid.
      !-------------------------------------------------------
      !-------------------------------------------------------------
      ! Calculate wind and pressure fields at model nodal points.
      !
      ! Note: the asymmetric vortex wind speed is reduced from the
      ! top of the surface layer to the surface, then converted from
      ! a 1-minute (max sustained) to a 10-minute average prior to
      ! adding the translational velocity in subroutine uvp.
      !-------------------------------------------------------------
      stormMotion = 1.5d0*(SQRT(uTransNow**2.d0+vTransNow**2.d0))**0.63d0
      dirNow = rad2deg * ATAN2(uTransNow, vTransNow)
      if (dirNow .lt. 0.d0) dirNow = dirNow + 360.d0
      stormMotionU = sin(dirNow/rad2deg)*stormMotion
      stormMotionV = cos(dirNow/rad2deg)*stormMotion
      call setVortex(Pn,Pc,cLat,cLon)

      DO i=1,np
         CALL uvpr(dist(i),aziangle(i),crmaxw(i),crmaxwTrue(i),
     &        cHollBs(i),cVmwBL(i),cPhiFactor(i),stormMotionU,
     &        stormMotionV,geofactor,wvnx(i),wvny(i),press(i))
         !-------------------------------------------
         ! Convert atmospheric pressure (Pascals) to
         ! atmospheric pressure-induced water surface
         ! elevation (meters).
         !-------------------------------------------
         press(i) = press(i) / RhoWatG

      END DO

      !-----------------------------------
      ! NHC advisory best-track i/o format
      !-----------------------------------
      firstCall = 1
      told=TIMELOC

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
999   RETURN
C***********************************************************************
      END SUBROUTINE NWS20GET
C***********************************************************************



C     ----------------------------------------------------------------
C      S U B R O U T I N E     N W S  1 5  I N I T
C     ----------------------------------------------------------------
C     jgf50.38.05: This subroutine checks input data, allocates memory
C     verifies that input files are all present, etc.
C     ----------------------------------------------------------------
      subroutine NWS15INIT(timeloc)
      use sizes, only : globaldir, mnp
      use global, only : ihot, nws, openFileForRead, nabout
      implicit none
C
      real(8), intent(in) :: timeloc ! adcirc time in seconds since coldstart
C
      logical :: fileFound   ! true if the file could be found
      integer :: errorio     ! .gt. 1 if there was an i/o error
      integer :: currentLine ! line from fort.22 being processed, used in error msgs
      logical :: cycleTimeFound ! true if the adcirc time falls into one
                                ! of the hwind time intervals
      integer :: i           ! loop counter
C
      call setMessageSource("nws15init")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      errorio = 0
      call openFileForRead(22,trim(globaldir)//"/fort.22",errorio)
      if (errorio.gt.0) call windTerminate()
      !
      ! read fort.22 to determine how many HWind files there are
      currentLine = 1
      do i=1,3
         ! throw header away this time
         read(22,'(A)',end=246,err=248,iostat=errorio) scratchMessage
         currentLine = currentLine + 1
      end do
      numFiles = 0
      do   ! loop until we run out of data, just counting the files this time
         read(22,'(A)',end=123,err=248,iostat=errorio) scratchMessage
         currentLine = currentLine + 1
         numFiles = numFiles + 1
      end do
 123  continue ! jump to here when we've run out of lines in the fort.22
      allocate(hWindFiles(numFiles))
      !
      ! now go back to the beginning of the file and read in the data
      rewind(22)
      currentLine = 1
      read(22,'(A)',end=246,err=248,iostat=errorio) metComment
      call logMessage(ECHO,"metComment: "//trim(metComment))
      currentLine = currentLine + 1
      read(22,*,end=246,err=248,iostat=errorio) hWindMultiplier
      currentLine = currentLine + 1
      write(scratchMessage,'("hWindMultiplier=",E15.8)') hWindMultiplier
      call logMessage(ECHO,scratchMessage)
      read(22,'(A)',end=248,err=248,iostat=errorio) scratchMessage
      currentLine = currentLine + 1
      pressureWindRelationship
     &      = scratchMessage(1:index(adjustl(scratchMessage)," "))
      call logMessage(ECHO,"pressureWindRelationship: "
     &     //trim(pressureWindRelationship))
      do i=1,numFiles
         read(22,*,end=246,err=248,iostat=errorio)
     &     hWindFiles(i)%cycleTime,
     &     hWindFiles(i)%Pc,
     &     hWindFiles(i)%rampVal,
     &     hwindFiles(i)%file_name
         write(scratchMessage,
     &   '("time=",E15.8," Pc=",E15.8," ramp=",E15.8," file_name=",A)')
     &       hWindFiles(i)%cycleTime, hWindFiles(i)%Pc,
     &       hWindFiles(i)%rampVal, trim(hwindFiles(i)%file_name)
         call logMessage(ECHO,scratchMessage)
         currentLine = currentLine + 1
      end do
      close(22)
      !
      ! log the number of HWind files that were specified
      write(scratchMessage,
     &    '("There were ",I3," HWind files specified.")') numFiles
      call logMessage(INFO,scratchMessage)
      !
      ! need at least 2 files so we have something to interpolate between
      ! in time
      if (numFiles.lt.2) then
         call allMessage(ERROR,"Only one HWind file was specified.")
         call allMessage(ERROR,
     &    "ADCIRC requires at least two HWind files for interpolation.")
         call windTerminate()
      endif
      !
      ! check to make sure that the cycleTimes are
      ! monotonically increasing
      do i=1,numFiles-1
         if (hWindFiles(i)%cycleTime.ge.hWindFiles(i+1)%cycleTime) then
            write(scratchMessage,
     &    '("Cycle time ",I3," must be less cycle time ",I3,".")')i,i+1
            call allMessage(ERROR,scratchMessage)
            call windTerminate()
         endif
      end do
      !
      ! check to make sure the specified pressure-wind relationship
      ! is supported by ADCIRC
      select case(trim(pressureWindRelationship))
      case("dvorak","knaffzehr","specifiedPc","background")
         call logMessage(INFO,"ADCIRC will use the "//
     &   trim(pressureWindRelationship)//" pressure-wind relationship.")
      case default
         call allMessage(ERROR,"The "//trim(pressureWindRelationship)//
     &      " pressure-wind relationship is not supported by ADCIRC.")
         call windTerminate()
      end select
      !
      ! check to make sure that all the HWind files are actually there;
      ! we assume that the file names contain the full path to the file;
      ! if the file name contains forward slashes, it must be surrounded
      ! by quotation marks
      do i=1,numFiles-1
         fileFound = .false.
         inquire(file=hWindFiles(i)%file_name,exist=fileFound)
         if (fileFound.eqv..false.) then
            call allMessage(ERROR,"The HWind file "//
     &          trim(hWindFiles(i)%file_name)//" was not found.")
            call windTerminate()
         endif
      end do
      call logMessage(INFO,"All specified HWind files were found.")
      hWindFiles(:)%loaded = .false.
      !
      ! allocate memory for mesh coordinates in mercator projection
      allocate(x_mercator(mnp),y_mercator(mnp))
      ! allocate memory for lower left HWind grid indices assigned to each mesh node
      allocate(x_indices(mnp),y_indices(mnp))
      ! allocate memory for the interpolation weights
      allocate(w1(mnp),w2(mnp),w3(mnp),w4(mnp))
      ! allocate memory for wind u,v velocity work arrays
      allocate(wvnx_work(mnp),wvny_work(mnp))
      ! allocate memory for wind speed; this is used to find Vmax and Rmax
      allocate(windSpeeds(mnp))
      ! allocate memory to keep track of which mesh nodes are inside
      ! the hwind bounding box
      allocate(inside(mnp))
      !
      ! relate timing in the files to current ADCIRC time
      hWindFiles(:)%cycleTime = hWindFiles(:)%cycleTime * 3600.d0 ! convert hours to seconds
      if  ((nws.lt.0).and.(ihot.ne.0)) then
         ! if times are relative to hotstart time, add the hotstart time
         ! so that we can consider the times relative to the coldstart
         ! time from here on
         hWindFiles(:)%cycleTime = hWindFiles(:)%cycleTime + timeloc
      endif
      !
      ! locate the cycle that is relevant to current ADCIRC time
      cycleTimeFound = .false.
      do i=1,numFiles-1
         if ( (hWindFiles(i)%cycleTime.le.timeloc).and.
     &        (hWindFiles(i+1)%cycleTime.gt.timeloc) ) then
            cycleTimeFound = .true.
            exit
         endif
      end do
      ! if we couldn't find a relevant time interval, game over
      if ( cycleTimeFound.eqv..false. ) then
         if ( timeloc.lt.hWindFiles(1)%cycleTime ) then
            call allMessage(ERROR,
     & "The current ADCIRC time is earlier than the first HWind file.")
         endif
      if ( timeloc.gt.hWindFiles(numFiles)%cycleTime ) then
            call allMessage(ERROR,
     &    "The current ADCIRC time is later than the last HWind file.")
         endif
         call windTerminate()
      endif
      currentCycle = i
      !
      ! load up the data for this time interval
      call loadHWindFile(currentCycle)
      call loadHWindFile(currentCycle+1)
      !
      ! initialize eye locations for sector based wind drag (Powell)
      EyeLatR = 0.d0
      EyeLonR = 0.d0

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
C
C     We jump to this section if there was an error reading a file.
246   call allMessage(ERROR,"Unexpectedly reached end-of-file.") ! END jumps here
248   call allMessage(ERROR,"I/O error during file access.") !  ERR jumps here
      call allMessage(INFO,
     &   "Check the fort.16 file for more information.")
      if (nabout.gt.0) then
         call allMessage(INFO,
     &  "Reducing the value of NABOUT to 0"
     &  //" will maximize the information written to the fort.16 file,"
     &  //" which may aid in troubleshooting this issue.")
      endif
      write(scratchMessage,'("Could not read line ",I2," of fort.22.")')
     &  currentLine
      call allMessage(ERROR,scratchMessage)
      write(scratchMessage,'(A,I3,A)')
     &       'The value of the i/o error flag was ',errorio,'.'
      call allMessage(ERROR,scratchMessage)
      CALL windTerminate()
C     ----------------------------------------------------------------
      END SUBROUTINE NWS15INIT
C     ----------------------------------------------------------------


C***********************************************************************
C      S U B R O U T I N E     N W S  1 5  G E T
C***********************************************************************
C      jgf50.38.05: Developing this subroutine to read in HWind formatted
C      files. The fort.22 will consist of 2 sections: a header section
C      that has parameters for the dataset as a whole, and then the list
C      of times, in hours, and the corresponding file names that
C      correspond to those times.
C
C      HWind-specific fort.22 file format is as follows:
C      ! first line is a comment line, max length 1024 characters
C      1.0    ! 2nd line is a velocity magnitude multiplier
C      dvorak ! 3rd line: one word for the pressure-wind relationship
C      0.0  -1  0.0  al092011_0828_1330  ! time (hours), Pc (mb), ramp mult, filename
C      6.0  -1  0.5  al092011_0828_1930
C      12.0 -1  1.0  al092011_0829_0130
C      (etc)
C
C      HWind File Format notes:
C         a. the (u,v) data are on a regular grid
C         b. the regular grid is a mercator projection with origin at storm center
C         c. the mercator grid spacing is in meters and is uniform in x and y (dx=dy)
C         d. the dimensions (nx,ny) of the mercator grid are equal (nx=ny)
C         e. the grid dimensions change from snapshot to snapshot
C            for example, the first shapshot may be 161x161
C                         while the 2nd snapshot may be 121x121 (etc)
C         f. sequential hwind snapshots will not be evenly spaced in time
C         g. hwind data do not contain barometric pressure information
C***********************************************************************
      subroutine NWS15GET(wvnx, wvny, press, timeloc)
      use sizes, only : mnp
      use global, only : mb2pa, m2nm, nm2m
      use mesh, only : slam, sfea, np, dp
      implicit none
      real(sz), intent(out) :: wvnx(:)  ! wind u-velocity (m/s)
      real(sz), intent(out) :: wvny(:)  ! wind v-velocity (m/s)
      real(sz), intent(out) :: press(:) ! barometric pressure (mH2O)
      real(8), intent(in) :: timeloc    ! adcirc time, seconds since coldstart
C
      real(sz) :: wtratio  ! time weighting factor for time interpolation
      real(sz) :: cLonNow  ! time interpolated longitude at center of storm
      real(sz) :: cLatNow  ! time interpolated latitude at center of storm
      real(sz) :: rampValNow  ! time interpolated met ramp value
      real(sz) :: k_0      ! constant used in mercator projection
      real(sz) :: y_mercator_0 ! mercator y coord at center of storm
      real(sz) :: rmaxNow  ! (nm) used in estimating central pressure
      real(sz) :: pc       ! (mbar) used to find pressure field
      real(sz) :: B        ! Holland cyclone shape parameter
      real(sz) :: dist     ! (m) distance from storm center to mesh node
      real(sz) :: dxNow,dyNow ! (nm) x and y distance from storm center rmaxNow
      real(sz) :: vmaxNow  ! (m/s) time interpolated fulldomain vmax
      real(sz) :: oldDX, oldDY ! (nm) distances to rmax in last hwind cycle
      real(sz) :: newDX, newDY ! (nm) distances to rmax in next hwind cycle
      real(sz) :: angle    ! (nm) angle to location of mesh vmax (on subdomain in parallel)
      real(sz) :: maxSpeed ! (m/s) max wind speed on mesh (subdomain in parallel)
      integer  :: maxNode  ! mesh node on which maxSpeed occurs
      integer  :: i        ! loop counter
C
      call setMessageSource("nws15get")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      !
      ! check to see if the current ADCIRC time is still in the interval of
      ! the HWind cycles we are currently accessing ... if the ADCIRC
      ! time has gone beyond this interval, load up the next file
      if ( timeloc.gt.hWindFiles(currentCycle+1)%cycleTime) then
         ! "unload" the stale hwind data
         deallocate(hWindFiles(currentCycle)%u)
         deallocate(hWindFiles(currentCycle)%v)
         hWindFiles(currentCycle)%loaded = .false.
         ! update cycle number
         currentCycle = currentCycle + 1
         ! load up the new data, unless we've gone past the last file
         if ((currentCycle+1).le.numFiles) then
            call loadHWindFile(currentCycle+1)
         else
            call allMessage(ERROR,
     &  "The ADCIRC time has exceeded the time of the last HWind file.")
            write(scratchMessage,
     &    '("ADCIRC time=",E15.3," (sec), last HWind time=",E15.3,".")')
     &      timeloc, hWindFiles(numFiles)%cycleTime
            call windTerminate()
         endif
      endif
      !
      ! compute the time weighting factor for time interpolation
      wtratio = (timeloc - hWindFiles(currentCycle)%cycleTime) /
     &          ( hWindFiles(currentCycle+1)%cycleTime -
     &            hWindFiles(currentCycle)%cycleTime )
      !
      ! time-interpolate the location of the storm center
      cLonNow = hWindFiles(currentCycle)%cLon +
     &     wtratio * (hWindFiles(currentCycle+1)%cLon -
     &                 hWindFiles(currentCycle)%cLon )
      cLatNow = hWindFiles(currentCycle)%cLat +
     &     wtratio * (hWindFiles(currentCycle+1)%cLat -
     &                 hWindFiles(currentCycle)%cLat )
      !
      ! set eye location for sector based wind drag
      if ( (cLatNow.ne.EyeLatR(3)).or.(cLonNow.ne.EyeLonR(3))) then
         EyeLatR(1) = EyeLatR(2)
         EyeLonR(1) = EyeLonR(2)
         EyeLatR(2) = EyeLatR(3)
         EyeLonR(2) = EyeLonR(3)
         EyeLatR(3) = cLatNow
         EyeLonR(3) = cLonNow
         FoundEye = .TRUE.
      endif
      !write(scratchMessage,*) "EyeLatR: ",(EyeLatR(i),i=1,3) !jgfdebug
      !call logMessage(DEBUG,scratchMessage) !jgfdebug
      !write(scratchMessage,*) "EyeLonR: ",(EyeLatR(i),i=1,3) !jgfdebug
      !call logMessage(DEBUG,scratchMessage) !jgfdebug
      !
      ! time-interpolate the vmax, rmax, and angle to rmax
      ! vmaxNow represents the current vmax on the whole domain
      ! ... not just  on this particular subdomain (in parallel)
      vmaxNow = hWindFiles(currentCycle)%vmax +
     &     wtratio * (hWindFiles(currentCycle+1)%vmax -
     &                 hWindFiles(currentCycle)%vmax )
      rmaxNow = hWindFiles(currentCycle)%rmax +
     &     wtratio * (hWindFiles(currentCycle+1)%rmax -
     &                 hWindFiles(currentCycle)%rmax )
      oldDX = cos(hWindFiles(currentCycle)%rmaxAngle*deg2rad)
      oldDY = sin(hWindFiles(currentCycle)%rmaxAngle*deg2rad)
      newDX = cos(hWindFiles(currentCycle+1)%rmaxAngle*deg2rad)
      newDY = sin(hWindFiles(currentCycle+1)%rmaxAngle*deg2rad)
      dxNow = oldDX + wtratio * (newDX - oldDX)
      dyNow = oldDY + wtratio * (newDY - oldDY)
      angleNow = 360.0d0 + rad2deg * ATAN2(dxNow,dyNow)
      if ( angleNow > 360.d0) angleNow = angleNow - 360.d0
      !
      ! time-interpolate the ramp value
      rampValNow = hWindFiles(currentCycle)%rampval +
     &     wtratio * (hWindFiles(currentCycle+1)%rampval -
     &                 hWindFiles(currentCycle)%rampval )
      !
      ! transform the ADCIRC mesh into mercator with origin at the
      ! actual current storm center
      k_0 = cos(deg2rad*cLatNow)
      x_mercator = k_0 * Rearth * (slam - deg2rad*cLonNow)
      y_mercator_0 = k_0 *
     &    Rearth * log(tan(0.25d0*pi+0.5d0*deg2rad*cLatNow))
      y_mercator = k_0 * Rearth * log(tan(0.25d0*pi+0.5d0*sfea))
     &    - y_mercator_0
      !
      ! initialize arrays
      wvnx = 0.d0
      wvny = 0.d0
      !------------------------------------------------------------
      ! interpolate from most recent HWind dataset onto ADCIRC mesh
      !------------------------------------------------------------
      call hWindInterpolateVelocity(wvnx, wvny, currentCycle,
     &    rampValNow, (1.d0-wtratio) )
      !------------------------------------------------------------
      ! interpolate from next (future) HWind dataset onto ADCIRC mesh
      !------------------------------------------------------------
      call hWindInterpolateVelocity(wvnx, wvny, currentCycle+1,
     &    rampValNow, wtratio)
      !
      ! determine maximum interpolated wind speeds on adcirc mesh
      ! (or on this subdomain mesh if we are running in parallel)
      windSpeeds(:) = sqrt(wvnx(:)**2+wvny(:)**2)
      !
      ! if this run is for wind analysis
      ! write out Rmax (nautical miles) and azimuthal angle (degrees),
      if (writeFullCircleRmaxes.eqv..true.) then
         maxSpeed = maxval(windSpeeds)
         maxNode = maxloc(windSpeeds,1)
         dist = m2nm * sphericalDistance(slam(maxNode)-deg2rad*cLonNow,
     &      sfea(maxNode)-deg2rad*cLatNow,
     &      cLatNow, rad2deg*sfea(maxNode))
         angle = 360.0d0 + rad2deg
     &       * ATAN2((slam(maxNode)-deg2rad*cLonNow),
     &                (sfea(maxNode)-deg2rad*cLatNow))
         if ( angle > 360.d0) angle = angle - 360.d0
         write(444,*) dist, angle
      endif
      !
      ! use max wind speed to estimate minimum central pressure
      ! using the specified pressure-wind relationship
      select case(trim(pressureWindRelationship))
         case("dvorak")
            pc = 1015.d0 - (vmaxNow/3.92d0)**(1.0d0/0.644d0);
         case("knaffzehr")
            pc = 1010.d0 - (vmaxNow/2.3d0)**(1.0d0/0.76d0);
         case("specifiedPc")
            ! time-interpolate the specified central pressure value
            pc = hWindFiles(currentCycle)%Pc +
     &           wtratio * (hWindFiles(currentCycle+1)%Pc -
     &                      hWindFiles(currentCycle)%Pc )
         case("background")
            pc = 1013.d0
         case default
            ! should be unreachable; already quality checked in nws15init
            call allMessage(ERROR,"The "
     &          //trim(pressureWindRelationship)//
     &        " pressure-wind relationship is not supported by ADCIRC.")
            call windTerminate()
      end select
      !
      ! use central pressure and max wind speed to estimate the Holland B value
      B = vmaxNow**2*RhoAir*EXP(1.d0)/((1013.d0 - Pc)*mb2pa)
      B = MAX( MIN(B,2.5d0), 1.0d0) ! limit B to range 1.0->2.5
      !
      ! using Vmax, Rmax, Pc, and B ... estimate the pressure
      ! field using Holland's curve fit ... even compute pressure outside
      ! the bounding box of the hwind data, just so we don't get a
      ! destabilizing step change in the barometric pressure
      select case(trim(pressureWindRelationship))
         case("dvorak","knaffzehr","specifiedPc")
            do i=1,np
               dist = sphericalDistance(slam(i)-deg2rad*cLonNow,
     &               sfea(i)-deg2rad*cLatNow, cLatNow, rad2deg*sfea(i))
               ! the following assumes the rmax is uniform all the way
               ! around the storm
               press(i) = (Pc + (1013.d0 - Pc)
     &               * EXP(-(rmaxNow*nm2m/dist)**B))
            end do
         case("background")
            press = 1013.d0
         case default
            ! should be unreachable
      end select
      !
      ! apply ramp value from fort.22 to pressure field
      press = ( press - 1013.d0 ) * rampValNow + 1013.d0
      !
      ! convert pressure from mb to mH2O
      press = mb2pa * press / ( RhoWat0 * g )
      ! convert from 1 minute avg winds to 10 minute averaged winds
      ! TODO: confirm that HWind data are 1 minute averaged winds
      wvnx = wvnx * one2ten
      wvnx = wvnx * one2ten
C
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C     ----------------------------------------------------------------
      END SUBROUTINE NWS15GET
C     ----------------------------------------------------------------



C     ----------------------------------------------------------------
C      S U B R O U T I N E    L O A D  H W I N D  F I L E
C     ----------------------------------------------------------------
C
C     jgf50.38.05: Subroutine to load up the wind vectors data
C     from a single HWind file.
C
C     ----------------------------------------------------------------
      subroutine loadHWindFile(cn)
      use global, only : openFileForRead, nabout, m2nm
      implicit none
      integer, intent(in) :: cn ! cycle number of file to load
C
      real(sz), allocatable    :: windMag(:,:) ! hwind speeds, used to find hwind vmax
      real(sz)                 :: dx, dy ! distances (m) to location of vmax
      integer                  :: maxIndices(2) !i,j indices of hwind vmax
      real(sz)                 :: throwAwayCoordinate ! dummy variable
      integer                  :: numValues ! total number of u,v velocity pairs
      integer                  :: indexLimit ! abs(array bounds)
      integer                  :: errorio ! nonzero if there was an i/o error
      integer                  :: currentLine ! used in error messages
      integer                  :: dataset ! used in error messages
      integer                  :: counter ! used to count uvWind array values
      integer                  :: i, j    ! loop counters
      character(1)             :: paren
      character(1)             :: comma
C
      call setMessageSource("loadHWindFile")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      call logMessage(INFO,"Loading the HWind file '"
     &                     //trim(hWindFiles(cn)%file_name)//"'.")
      call openFileForRead(220,trim(hWindFiles(cn)%file_name), errorio)
      if (errorio.gt.0) then
         call windTerminate()
      endif
      currentLine = 1
      dataset = 0
      !
      ! SURFACE WIND COMPONENTS FOR HURRICANE  name_0825_22_30
      read(220,'(A)',end=246,err=248,iostat=errorio) scratchMessage
      call logMessage(ECHO,"HWind comment line: "//trim(scratchMessage))
      !
      ! DX=DY= 6.02280 KILOMETERS.
      currentLine = currentLine + 1
      read(220,'(7X,F8.5,A)',end=246,err=248,iostat=errorio)
     &   hWindFiles(cn)%dh
      write(scratchMessage,'("Mercator grid spacing is ",F8.5," km.")')
     &    hWindFiles(cn)%dh
      call logMessage(ECHO,scratchMessage)
      !
      ! STORM CENTER LOCALE IS -77.2730 EAST LONGITUDE and  27.2950 NORTH LATITUDE ... STORM CENTER IS AT (X,Y)=(0,0)
      currentLine = currentLine + 1
      read(220,'(23X,F9.4,18X,F9.4,A)',end=246,err=248,iostat=errorio)
     &   hWindFiles(cn)%cLon, hWindFiles(cn)%cLat ! degrees
      write(scratchMessage,'("center lon=",E15.8," center lat=",E15.8)')
     &   hWindFiles(cn)%cLon, hWindFiles(cn)%cLat
      call logMessage(ECHO,scratchMessage)
      ! dimensions
      do i=1,4
         currentLine = currentLine + 1
         read(220,'(A)',end=246,err=248,iostat=errorio) scratchMessage
         call logMessage(ECHO,"Coordinate: "//trim(scratchMessage))
         currentLine = currentLine + 1
         read(220,*,end=246,err=248,iostat=errorio) hWindFiles(cn)%nh
         write(scratchMessage,'("HWind dimension read: ",i3,".")')
     &      hWindFiles(cn)%nh
         call logMessage(ECHO,scratchMessage)
         ! we won't need the actual coordinates
         currentLine = currentLine + 1
         read(220,'(6F13.4)',end=246,err=248,iostat=errorio)
     &      (throwAwayCoordinate,j=1,hWindFiles(cn)%nh)
      end do
      !
      ! read the actual vector data
      !
      ! SURFACE WIND COMPONENTS ... M/S ... COMPLEX ARRAY W=(U,V)
      currentLine = currentLine + 1
      read(220,'(A)',end=246,err=248,iostat=errorio) scratchMessage
      currentLine = currentLine + 1
      call logMessage(ECHO,"HWind sfc wind line: "
     &      //trim(scratchMessage))
      currentLine = currentLine + 1
      read(220,'(A)',end=246,err=248,iostat=errorio) scratchMessage
      currentLine = currentLine + 1
      call logMessage(ECHO,"HWind sfc wind dimensions line: "
     &       //trim(scratchMessage))
      indexLimit = (hWindFiles(cn)%nh - 1) / 2
      allocate(hWindFiles(cn)
     &    %u(-indexLimit:indexLimit,-indexLimit:indexLimit))
      allocate(hWindFiles(cn)
     &    %v(-indexLimit:indexLimit,-indexLimit:indexLimit))
      allocate(windMag(-indexLimit:indexLimit,-indexLimit:indexLimit))
      ! read all the wind vector values
      currentLine = currentLine + 1
      dataset = 1
      do j=-indexLimit,indexLimit
            ! ref bottom of http://www.aoml.noaa.gov/hrd/Storm_pages/grid.html
         read(220,'(2(A1,F13.5,A1,F13.5,A1))',end=246,err=248,iostat=errorio)
     &       (paren, hWindFiles(cn)%u(i,j), comma,
     &        hWindFiles(cn)%v(i,j), paren,
     &        i=-indexLimit, indexLimit)
         dataset = dataset + 1
      end do
      close(220)
      call logMessage(INFO,"Finished loading the HWind file '"
     &                     //trim(hWindFiles(cn)%file_name)//"'.")
      !
      ! convert mercator mesh spacing dh from km to m
      hWindFiles(cn)%dh = hWindFiles(cn)%dh * 1000.d0
      !
      ! set index limit according to values found in the file
      hWindFiles(cn)%indexLimit = indexLimit
      !
      ! find the maximum wind speed in this dataset
      windMag(:,:) = sqrt(hWindFiles(cn)%u(:,:)**2
     &    + hWindFiles(cn)%v(:,:)**2)
      hWindFiles(cn)%vmax = maxval(windMag)
      maxIndices = maxloc(windMag) ! returns indices as if array is 1-indexed
      !
      ! now calculate rmax in nautical miles
      dx = (real(maxindices(1)-(indexLimit+1)))*hWindFiles(cn)%dh
      dy = (real(maxindices(2)-(indexLimit+1)))*hWindFiles(cn)%dh
      hWindFiles(cn)%rmax = m2nm * sqrt(dx**2 + dy**2)
      !
      ! calculate angle at which rmax occurs
      hWindFiles(cn)%rmaxAngle = 360.0d0 + rad2deg * ATAN2(dx,dy)
      if (hWindFiles(cn)%rmaxAngle > 360.d0) then
         hWindFiles(cn)%rmaxAngle = hwindFiles(cn)%rmaxAngle - 360.d0
      endif
      deallocate(windMag)
      !
      ! set the state variable
      hWindFiles(cn)%loaded = .true.
C
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return
C
C     We jump to this section if there was an error reading a file.
246   call allMessage(ERROR,"Unexpectedly reached end-of-file.") ! END jumps here
248   call allMessage(ERROR,"I/O error during file access.") !  ERR jumps here
      call allMessage(ERROR,
     &   "Check the fort.16 file for more information.")
           if (nabout.gt.0) then
         call allMessage(INFO,
     &  "Reducing the value of NABOUT to 0"
     &  //" will maximize the information written to the fort.16 file,"
     &  //" which may aid in troubleshooting this issue.")
      endif
      if ( dataset.eq.0 ) then
         write(scratchMessage,
     &    '("Could not read line ",I2," from the file ",A,".")')
     &    currentLine, trim(hWindFiles(cn)%file_name)
         call allMessage(ERROR,scratchMessage)
      else
         write(scratchMessage,
     &     '("Could not read dataset ",I2," from the file ",A,".")')
     &    dataset, trim(hWindFiles(cn)%file_name)
         call allMessage(ERROR,scratchMessage)
      endif
      write(scratchMessage,'(A,I5,A)')
     &       'The value of the i/o error flag was ',errorio,'.'
      call allMessage(ERROR,scratchMessage)
      CALL windTerminate()
C
C     ----------------------------------------------------------------
      end subroutine loadHWindFile
C     ----------------------------------------------------------------


C     ----------------------------------------------------------------
C      S U B R O U T I N E
C          H W I N D   I N T E R P O L A T E   V E L O C I T Y
C     ----------------------------------------------------------------
C
C     jgf50.38.05: Interpolates two hwind velocity datasets in
C     space and time onto an adcirc mesh.
C
C     ----------------------------------------------------------------
      subroutine hWindInterpolateVelocity(wvnx, wvny, cn, rv, tw)
      use mesh, only : np
      implicit none
      real(sz), intent(out) :: wvnx(:) ! wind u-velocity
      real(sz), intent(out) :: wvny(:) ! wind v-velocity
      integer, intent(in)   :: cn      ! cycle number, or hwind file number
      real(sz), intent(in)  :: rv      ! ramp value from hwind file
      real(sz), intent(in)  :: tw      ! time weighting for time interpolation
C
      real(sz)              :: dh      ! readability shorthand for hwind dh
      integer               :: il      ! readability shorthand for hwind integerLimit
      real(sz), pointer     :: hwu(:,:)! readability shorthand for hwind u
      real(sz), pointer     :: hwv(:,:)! readability shorthand for hwind v
      real(sz)              :: oneOverW2 ! 1/w**2 ; const in interp weights
      integer               :: i       ! loop counter
C
      call setMessageSource("hwindinterpolatevelocity")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      dh = hWindFiles(cn)%dh
      hwu => hWindFiles(cn)%u
      hwv => hWindFiles(cn)%v
      il = hWindFiles(cn)%indexLimit
      !
      ! each mesh point falls within a box of the hwind data;
      ! we form an array of indices that represent the HWind grid indices
      ! of the lower left corners of the boxes that each mesh point falls in
      x_indices = floor(x_mercator/dh)
      y_indices = floor(y_mercator/dh)
      !
      ! calculate bilinear interpolation weights
      oneOverW2 = 1.d0 / real(dh**2)
      where ( (x_indices.ge.-il).and.(y_indices.ge.-il).and.
     &        (x_indices.le.(il-1)).and.(y_indices.le.(il-1)) )
         ! the adcirc node is inside the bounding box of the HWind data;
         ! we will calculate interpolation weights for these nodes
         inside = .true.
         w1 = ((dh*(x_indices+1) - x_mercator)
     &       *  (dh*(y_indices+1) - y_mercator)) * oneOverW2
         w2 = ((x_mercator - dh*x_indices)
     &       * (dh*(y_indices+1) - y_mercator) ) * oneOverW2
         w3 = ((x_mercator - dh*x_indices)
     &       * (y_mercator - dh*y_indices) ) * oneOverW2
         w4 = ((dh*(x_indices+1) - x_mercator )
     &       * (y_mercator - dh*y_indices) ) * oneOverW2
      elsewhere
         inside = .false.
      end where
      !
      ! apply interpolation weights to perform bilinear interpolation
      ! of wind velocities
      do i=1,np
         if (inside(i).eqv..false.) then
            cycle
         endif
         wvnx_work(i) = w1(i)*hwu(x_indices(i),y_indices(i))
     &           + w2(i)*hwu(x_indices(i)+1,y_indices(i))
     &           + w3(i)*hwu(x_indices(i)+1,y_indices(i)+1)
     &           + w4(i)*hwu(x_indices(i),y_indices(i)+1)
         wvny_work(i) = w1(i)*hwv(x_indices(i),y_indices(i))
     &           + w2(i)*hwv(x_indices(i)+1,y_indices(i))
     &           + w3(i)*hwv(x_indices(i)+1,y_indices(i)+1)
     &           + w4(i)*hwv(x_indices(i),y_indices(i)+1)
      end do
      !
      ! apply ramp value, time weighting, and user specified
      ! wind multiplier; add the contribution of this dataset to the total
      ! at this mesh node
      where (inside.eqv..true.)
         wvnx = wvnx + ( wvnx_work * rv * tw * hWindMultiplier )
         wvny = wvny + ( wvny_work * rv * tw * hWindMultiplier )
      end where
C
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C
C     ----------------------------------------------------------------
      end subroutine hWindInterpolateVelocity
C     ----------------------------------------------------------------

C
C ---------------------------------------------------------------------
C  subroutine to initalize GFDL inputs    
C       Written By:  Chris Massey, USACE-ERDC-CHL, Vicksburg, MS
C       Date:  2/7/2013
C ---------------------------------------------------------------------
C    Opens the fort.22 file
C          Line 1 -- ASCII Text Header
C          Line 2 -- Windmultiplier Value (real number)
C          Line 3 -- Maximum Extrapolation distance
C          Lines 4-end -- cycleTime  rampValue  filename
C

       subroutine init_gfdl(timeloc)
       use sizes, only : globaldir, sz
       use global, only : ihot, nws, nabout, openFileForRead,
     &                setMessageSource,unsetMessageSource,allMessage,
     &                logMessage,scratchMessage,DEBUG,ECHO,INFO,
     &                ERROR,WARNING,PRBCKGRND,RHOWAT0,G,
     &                wvnx1,wvny1,prn1,wvnx2,wvny2,prn2
      use mesh, only : np, x, y, slam0, sfea0
      implicit none

      real(sz), intent(in) :: timeloc ! adcirc time in seconds since coldstart
C
      logical :: fileFound   ! true if the file could be found
      integer :: errorio     ! .gt. 1 if there was an i/o error
      integer :: currentLine ! line from fort.22 being processed, used in error msgs
      logical :: cycleTimeFound ! true if the adcirc time falls into one
                                ! of the gfdl time intervals
      integer :: i           ! loop counter

C
      call setMessageSource("gfdl_init")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      errorio = 0
      call openFileForRead(22,trim(globaldir)//"/fort.22",errorio)
      if (errorio.gt.0) call windTerminate()
      !
      ! read fort.22 to determine how many GFDL files there are
      currentLine = 1
      do i=1,3
         ! throw header away this time
         read(22,'(A)',end=246,err=248,iostat=errorio) scratchMessage
         currentLine = currentLine + 1
      end do
      numFiles = 0      
      do   ! loop until we run out of data, just counting the files this time
         read(22,'(A)',end=123,err=248,iostat=errorio) scratchMessage
         currentLine = currentLine + 1
         numFiles = numFiles + 1
      end do
 123  continue ! jump to here when we've run out of lines in the fort.22
      allocate( gfdl_files(numFiles) )

      !
      ! now go back to the beginning of the file and read in the data
      rewind(22)
      currentLine = 1
      read(22,'(A)',end=246,err=248,iostat=errorio) metComment
      call logMessage(ECHO,"metComment: "//trim(metComment))
      currentLine = currentLine + 1
      read(22,*,end=246,err=248,iostat=errorio) GFDL_WindMultiplier
      currentLine = currentLine + 1
      write(scratchMessage,'("WindMultiplier=",E15.8)') GFDL_WindMultiplier
      call logMessage(ECHO,scratchMessage)
      read(22,*,end=246,err=248,iostat=errorio) GFDL_max_extrap_dist
      currentLine = currentLine + 1
      write(scratchMessage,
     &  '("GFDL Maximum Extrapolation Distance=",E15.8)') GFDL_max_extrap_dist
      call logMessage(ECHO,scratchMessage)

      do i=1,numFiles
         read(22,*,end=246,err=248,iostat=errorio)
     &     GFDL_Files(i)%cycleTime,
     &     GFDL_Files(i)%rampVal,
     &     GFDL_Files(i)%file_name

         GFDL_Files(i)%max_extrap_dist = GFDL_max_extrap_dist
         GFDL_Files(i)%numFiles = numFiles
         write(scratchMessage,
     &   '("time=",E15.8," ramp=",E15.8," file_name=",A)')
     &       GFDL_Files(i)%cycleTime, 
     &       GFDL_Files(i)%rampVal, trim(GFDL_Files(i)%file_name)
         call logMessage(ECHO,scratchMessage)
         currentLine = currentLine + 1
      end do
      close(22)
 
      ! log the number of GFDL Met files that were specified
      write(scratchMessage,
     &    '("There were ",i0," GFDL Met files specified.")') numFiles
      call logMessage(INFO,scratchMessage)
      !
      ! need at least 2 files so we have something to interpolate between
      ! in time
      if (numFiles.lt.2) then
         call allMessage(ERROR,"Only one GFDL Met file was specified.")
         call allMessage(ERROR,
     &    "ADCIRC requires at least two GFDL Met files for interpolation.")
         call windTerminate()
      endif
      !
      ! check to make sure that the cycleTimes are
      ! monotonically increasing
      do i=1,numFiles-1
         if (GFDL_Files(i)%cycleTime.ge.GFDL_Files(i+1)%cycleTime) then
            write(scratchMessage,
     &  '("Cycle time ",i0," must be less cycle time ",i0,".")')i,i+1
            call allMessage(ERROR,scratchMessage)
            call windTerminate()
         endif
      end do

      !
      ! check to make sure that all the GFDL met files are actually there;
      ! we assume that the file names contain the full path to the file;
      ! if the file name contains forward slashes, it must be surrounded
      ! by quotation marks
      do i=1,numFiles-1
         fileFound = .false.
         inquire(file=GFDL_Files(i)%file_name,exist=fileFound)
         if (fileFound.eqv..false.) then
            call allMessage(ERROR,"The GFDL Met file "//
     &          trim(GFDL_Files(i)%file_name)//" was not found.")
            call windTerminate()
         endif
      GFDL_Files(i)%loaded = .false.
      end do

      write(scratchMessage,'("All specified GFDL Met Files were found.")')
      call logMessage(INFO,scratchMessage)
      GFDL_Files(:)%loaded = .false.

      !
      ! relate timing in the files to current ADCIRC time
      GFDL_Files(:)%cycleTime = GFDL_Files(:)%cycleTime * 3600.d0 ! convert hours to seconds
      if  ((nws.lt.0).and.(ihot.ne.0)) then
         ! if times are relative to hotstart time, add the hotstart time
         ! so that we can consider the times relative to the coldstart
         ! time from here on
         GFDL_Files(:)%cycleTime = GFDL_Files(:)%cycleTime + timeloc
      endif

      !
      ! locate the cycle that is relevant to current ADCIRC time
      cycleTimeFound = .false.
      do i=1,numFiles-1
         if (( GFDL_Files(i)%cycleTime.le.timeloc).and.
     &        (GFDL_Files(i+1)%cycleTime.gt.timeloc) ) then
            cycleTimeFound = .true.
            exit
         endif
      end do

      currentCycle = i

      ! if we couldn't find a relevant time interval, then insert a blank wind field
      if (cycleTimeFound.eqv. .false.) then
         if ( timeloc.lt.GFDL_Files(1)%cycleTime ) then

! ... tcm 20130801 v51.06.10 replaced error exiting with insertion of blank snaps
            write(scratchMessage,'("The current ADCIRC time ",E15.8,
     &        " (days) is earlier than the first GFDL Met file time: ",
     &        E15.8, " (days). Inserting a blank wind snap.")')
     &        timeloc/86400.d0,GFDL_Files(1)%cycleTime/86400.d0
            call logMessage(WARNING,scratchMessage)
            currentCycle = 1
         endif

         if ( timeloc.gt.GFDL_Files(numFiles)%cycleTime ) then

            write(scratchMessage,'("The current ADCIRC time ",
     &        E15.8," (days) is later than the last GFDL Met file time: ",
     &        E15.8, " (days). Inserting a blank wind snap.")')
     &         timeloc/86400.d0,GFDL_Files(numFiles)%cycleTime/86400.d0
            call logMessage(WARNING,scratchMessage)
            currentCycle = numFiles-1
         endif

      endif

      if (cycleTimeFound.eqv. .true.) then
      ! load up the data for this time interval
         call GET_GFDL(NP,WVNX1,WVNY1,PRN1,x,y,
     &                GFDL_Files(currentCycle)%file_name,
     &                SLAM0,SFEA0,
     &                GFDL_Files(currentCycle)%max_extrap_dist)
         call GET_GFDL(NP,WVNX2,WVNY2,PRN2,x,y,
     &                GFDL_Files(currentCycle+1)%file_name,
     &                SLAM0,SFEA0,
     &                GFDL_Files(currentCycle)%max_extrap_dist)
      else
         WVNX1(:) = 0.d0
         WVNY1(:) = 0.d0
         PRN1(:) = PRBCKGRND*100.D0/(RHOWAT0*G)  !SET TO BACKGROUND PRESSURE 1013.0D0 ADJUSTED for ADCIRC units
      endif

      !If no errors then return now

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      !cleans up the setMessageSource set at the beginning of the subroutine
      call unsetMessageSource()
       
      return

C
C     We jump to this section if there was an error reading a file.
246   call allMessage(ERROR,"Unexpectedly reached end-of-file.") ! END jumps here
248   call allMessage(ERROR,"I/O error during file access.") !  ERR jumps here
      call allMessage(INFO,
     &   "Check the fort.16 file for more information.")
      if (nabout.gt.0) then
         call allMessage(INFO,
     &  "Reducing the value of NABOUT to 0"
     &  //" will maximize the information written to the fort.16 file,"
     &  //" which may aid in troubleshooting this issue.")
      endif
      write(scratchMessage,'("Could not read line ",i0," of fort.22.")')
     &  currentLine
      call allMessage(ERROR,scratchMessage)
      write(scratchMessage,'(A,i0,A)')
     &       'The value of the i/o error flag was ',errorio,'.'
      call allMessage(ERROR,scratchMessage)
      CALL windTerminate()
      call unsetMessageSource()
C ---------------------------------------------------------------------      
      end subroutine init_gfdl
C ---------------------------------------------------------------------
C

C
C ---------------------------------------------------------------------
C    subroutine to read in an ascii GFDL met file and interpolate
C       u-,v- winds and pressure onto the adcirc mesh.
C       Written By:  Chris Massey, USACE-ERDC-CHL, Vicksburg, MS
C       Date:  2/7/2013
C ---------------------------------------------------------------------
C
C   Each ASCII GFDL met file contains one or more nested grid data
C   where the nested grids are allowed to change in time.  Coarse grid
C   data is not stored where finer nest data is given.
C   Format of the file:
C   Line 1:  Number of grid cells (f10.4) NCELLS
C   Lines 2-NCELLS+1:  Have ten columns of data formatted as 10f10.4
C         1. u (m/sec)
C         2. v (m/sec)
C         3. Temperature  (K)
C         4. mixing ratio(kg/kg)
C         5. storm accum precipitation (cm)
C         6. sea level pressure (hPa)
C         7. longitude (decimal deg)
C         8. latitude (decimal deg)
C         9. hurricane hour
C        10. nest number (this is not always present)
C
       SUBROUTINE GET_GFDL(NADCPTS,ADC_U_WIND,ADC_V_WIND,ADC_PRESS,
     &                     ADC_X,ADC_Y,FNAME,SLAM0,SFEA0,
     &                     MAX_INTERP_DIST)
       USE SIZES
       USE MESH, ONLY : CPP
       USE GLOBAL, ONLY : DEG2RAD,PRBCKGRND,G,RHOWAT0,
     &        setMessageSource,unsetMessageSource,allMessage,
     &        logMessage,scratchMessage,DEBUG,ECHO,INFO,ERROR
       USE KDTREE2_MODULE     ! USED FOR FAST SEARCHING
       
       IMPLICIT NONE
       REAL(SZ), ALLOCATABLE :: XY_GFDL(:,:),XLOC(:)
       REAL(SZ), ALLOCATABLE :: PRESS_GFDL(:)
       REAL(SZ), ALLOCATABLE :: U_WIND_GFDL(:),V_WIND_GFDL(:)
       REAL(SZ), ALLOCATABLE :: WTS(:)
       INTEGER, ALLOCATABLE :: LOCS(:),II2(:)

       REAL(SZ) :: LON_GFDL,LAT_GFDL,TOL,DIST
       REAL(SZ) :: RPTS_GFDL,UTMP,VTMP,TTMP,MTMP,SPTMP,PTMP
       REAL(SZ) :: LONGTMP,LATTMP,HTMP,NTMP,XSTA,YSTA,MAXDIST
       REAL(SZ), INTENT(IN) :: SLAM0,SFEA0
       REAL(SZ), INTENT(IN) :: MAX_INTERP_DIST
       REAL(SZ), INTENT(OUT) :: ADC_U_WIND(NADCPTS),ADC_V_WIND(NADCPTS)
       REAL(SZ), INTENT(OUT) :: ADC_PRESS(NADCPTS)
       REAL(SZ), INTENT(IN)  :: ADC_X(NADCPTS),ADC_Y(NADCPTS)
!       REAL(SZ), PARAMETER ::  PI=3.141592653589793D0
!       REAL(SZ), PARAMETER ::  DEG2RAD = PI/180.D0  ! DEGREES TO RADIANS

       INTEGER :: NPTS_GFDL,SRCHDP,LUN,ERRORIO,ITC,I,IEK,IP,INCELL
       INTEGER, INTENT(IN) :: NADCPTS
       CHARACTER(240),INTENT(IN) :: FNAME
       TYPE(KDTREE2), POINTER :: TREE
       TYPE(KDTREE2_RESULT), ALLOCATABLE :: KDRESULTS(:)

      call setMessageSource("get_gfdl")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

       SRCHDP = 6  !THE NUMBER OF NEAREST NEIGHBORS TO FIND
       TOL = 1.0D-6  !TOLERANCE FOR ESTIMATING ZERO DISTANCE BETWEEN POINTS
       IP = 3        !POWER VALUE USED IN THE INVERESE WEIGHTED DISTANCE FORMULA
       LUN = 287  !FILE UNIT NUMBER TO OPEN FOR THE GFDL FILE

C     INVERSE WEIGHTED DISTANCE FORMULA
C     P_I = SUM( {J=1 TO SRCHDP} (P_J/(D_I_J)^IP) )/SUM( {J=1 TO SRCHDP} (1/(D_I_J)^IP ) )
C Where
C   P_J IS THE VALUE AT LOCATION J
C   D_I_J - DISTANCE BETWEEN POINT I AND POINT J
C   IP = POWER TO RAISE DISTANCE (THE LARGER IP THE MORE INFLUENCE CLOSER POINTS HAVE)

        call logMessage(INFO,"Loading the GFDL Met file '"
     &                     //trim(fname)//"'.")
       call openFileForRead(lun,trim(fname), errorio)
       if (errorio.gt.0) then
          call windTerminate()
       endif

!       OPEN(LUN,FILE=TRIM(FNAME),STATUS='OLD',ACTION='READ',
!     &            IOSTAT=ERRORIO)
!       IF (ERRORIO.NE.0) THEN
!          WRITE(*,25) TRIM(FNAME)
!       ENDIF
!25      FORMAT("COULD NOT OPEN THE FILE '",A,"'.")


       ! READ THE FIRST LINE OF THE GFDL OUTPUT WHICH IS THE NUMBER OF POINTS IN THIS RECORD
!       READ(LUN,'(F10.4)') RPTS_GFDL  
!       NPTS_GFDL = INT(RPTS_GFDL)
        read(LUN,'(I10)') NPTS_GFDL
        rpts_gfdl = real(npts_gfdl,sz)

!       tcm 20130507 -- changed the read above due to changes in GFDL ASCII files

       ALLOCATE( XY_GFDL(2,NPTS_GFDL) )
       ALLOCATE( PRESS_GFDL(NPTS_GFDL) )
       ALLOCATE( U_WIND_GFDL(NPTS_GFDL), V_WIND_GFDL(NPTS_GFDL) )

! READ IN EACH GFDL POINT LOCATION INFORMATION
       DO I=1,NPTS_GFDL
!          READ(LUN,'(10F10.4)') UTMP,VTMP,TTMP,MTMP,SPTMP,PTMP,
!     &                      LONGTMP,LATTMP,HTMP,NTMP
          READ(LUN,'(9F10.4)') UTMP,VTMP,TTMP,MTMP,SPTMP,PTMP,
     &                      LONGTMP,LATTMP,HTMP
          LON_GFDL = LONGTMP*DEG2RAD !CONVERT TO RADIANS
          LAT_GFDL = LATTMP*DEG2RAD  !CONVERT TO RADIANS
          U_WIND_GFDL(I) = UTMP
          V_WIND_GFDL(I) = VTMP
          PRESS_GFDL(I) = PTMP
          CALL CPP(XY_GFDL(1,I),XY_GFDL(2,I),LON_GFDL,LAT_GFDL,
     &             SLAM0,SFEA0)

       ENDDO !I
         
       CLOSE(LUN)  !CLOSE THIS GFDL DATA FILE

       write(scratchMessage,'("Read ",i0," GFDL Met Points.")')
     &    NPTS_GFDL
      call logMessage(ECHO,scratchMessage)

!     IF A MAXIMUM DISTANCE ISN'T SPECIFIED, THEN DETERMINE AN ESTIMATE
!     TO GFDL GRID SPACING IN ORDER TO USE FOR MAXIMUM DISTANCE (maxdist = 6*cell spacing)
      IF (MAX_INTERP_DIST == -99999.D0 ) THEN
          ALLOCATE( XLOC(NPTS_GFDL),II2(NPTS_GFDL) )

          MAXDIST = MINVAL(XY_GFDL(1,:))  !MINIMUM X VALUE
          XLOC(:) = XY_GFDL(1,:)
          II2(:) = 0
          WHERE( ABS(XLOC-MAXDIST) < TOL ) II2 = 1
          INCELL = MAX(SUM(II2)-1,1)
          MAXDIST = (MAXVAL(XY_GFDL(2,:))-MINVAL(XY_GFDL(2,:)))/REAL(INCELL)  !ESTIMATED CELL SPACING
          MAXDIST = 6.0D0*MAXDIST  !MAXIMUM DISTANCE TO LOOK AWAY
                   
          DEALLOCATE(XLOC)
          DEALLOCATE(II2)
       ELSE
          MAXDIST = MAX_INTERP_DIST
       ENDIF
       
       write(scratchMessage,'("Maximum GFDL Interpolation Distance = ",e24.16)'),maxdist
       call logMessage(ECHO,scratchMessage)

       
C... ALLOCATE SPACE FOR KDTREE SEARCH
C...  BE SURE THE MAXIMUM SEARCH DEPTH IS NOT LARGER THAN
C...  THE NUMBER OF ELEMENTS BEING KEPT
        IF (NADCPTS.LT.SRCHDP) SRCHDP = NADCPTS

C...  CREATE THE SEARCH TREE
      TREE => KDTREE2_CREATE(XY_GFDL,REARRANGE=.TRUE.,SORT=.TRUE.)


C...  ALLOCATE SPACE FOR THE SEARCH RESULTS FROM THE TREE
C...  THIS SPACE WILL BE DEALLOCATED LATER IN THE SUBROUTINE
      ALLOCATE(KDRESULTS(SRCHDP))

      ALLOCATE (WTS(SRCHDP) )
      ALLOCATE (LOCS(SRCHDP) )

! INITIALIZE ADCIRC VARIABLES
      ADC_PRESS(:) = PRBCKGRND  !SET TO BACKGROUND PRESSURE 1013.0D0
      ADC_U_WIND(:) = 0.D0
      ADC_V_WIND(:) = 0.D0

!     LOOP OVER EACH ADCIRC NODE
      DO I=1,NADCPTS      
         XSTA = ADC_X(I)  !ADCIRC NODE X-COORDINATE
         YSTA = ADC_Y(I)  !ADCIRC NODE Y-COORDIANTE

!        FIND THE SRCHDP CLOSEST POINTS TO THIS NODE
         CALL KDTREE2_N_NEAREST(TP=TREE,QV=(/XSTA,YSTA/),
     &                  NN=SRCHDP,RESULTS=KDRESULTS)

!
!     FOR SIMPLE NEAREST NEIGHBOR INTERPOLATION
!         ITC = 1
!         IEK = KDRESULTS(ITC)%IDX  !LOCATION OF CLOSEST POINT
!         ADC_PRESS(I) = PRESS_GFDL(IEK)
!         ADC_U_WIND(I) = U_WIND_GFDL(IEK)
!         ADC_V_WIND(I) = V_WIND_GFDL(IEK)

!     FOR INVERSE WEIGHTED INTERPOLATION
         !DETERMINE THE WEIGHTS
         DO ITC = 1,SRCHDP
            IEK = KDRESULTS(ITC)%IDX  !LOCATION OF ITC CLOSEST POINT
            DIST = SQRT(KDRESULTS(ITC)%DIS) !DISTANCE FROM THIS POINT TO THE LOCATION OF THE NODE
            IF (DIST > TOL) THEN  !CHECK TO SEE HOW FAR AWAY
               IF (DIST <= MAXDIST) THEN
                  WTS(ITC) = 1.D0/(DIST**IP)
               ELSE  !POINT IS TOO FAR AWAY TO BE USED
                  WTS(ITC) = -99999.D0
               ENDIF   
            ELSE  !POINT IS CLOSE ENOUGH TO USE AS DIRECT VALUE
               WTS(ITC) = -1.0D0
            ENDIF
            LOCS(ITC) = IEK
         ENDDO  !ITC

         !COMPUTED THE INVERSE DISTANCE WEIGHTED AVERAGE         
         UTMP = 0.D0
         VTMP = 0.D0
         PTMP = 0.D0
         TTMP = 0.D0
         ITC = 1;
         DO WHILE (ITC <= SRCHDP)
            ! IF POINT IS LOCATED AT AN ORIGINAL GFDL GRID THEN USE THAT VALUE ONLY
            IF (WTS(ITC) == -1.D0) THEN
               UTMP = U_WIND_GFDL(LOCS(ITC))
               VTMP = V_WIND_GFDL(LOCS(ITC))
               PTMP = PRESS_GFDL(LOCS(ITC))
               TTMP = 1.0D0
               ITC = SRCHDP + 20
            ELSEIF (WTS(ITC) == -99999.D0 ) THEN  !THIS LOCATION IS TOO FAR AWAY TO USE
               ITC = ITC + 1  
            ELSE  !COMPUTE INVERSE WEIGHTED AVERAGING
               UTMP = UTMP + U_WIND_GFDL(LOCS(ITC))*WTS(ITC)
               VTMP = VTMP + V_WIND_GFDL(LOCS(ITC))*WTS(ITC)
               PTMP = PTMP + PRESS_GFDL(LOCS(ITC))*WTS(ITC)  !SUM OF THE VALUES/WEIGHTS
               TTMP = TTMP + WTS(ITC)     ! SUM OF THE WEIGHTS
               ITC = ITC + 1
            ENDIF
         ENDDO  !ITC

         !SET THE NODAL VALUES TO THE WEIGHTED AVERAGE
         IF (TTMP.NE.0.D0) THEN
            ADC_U_WIND(I) = UTMP/TTMP
            ADC_V_WIND(I) = VTMP/TTMP
            ADC_PRESS(I) = (PTMP/TTMP)*100.D0/(RHOWAT0*G)   !CONVERT MILLIBARS TO METERS OF WATER
         ELSE  !NO POINT WAS CLOSE ENOUGH TO INTERPOLATE
            ADC_U_WIND(I) = 0.D0
            ADC_V_WIND(I) = 0.D0
            ADC_PRESS(I) = PRBCKGRND*100.D0/(RHOWAT0*G)
         ENDIF
         
      ENDDO !I LOOP OVER ADCIRC NODES

! CLEAN UP MEMORY
      DEALLOCATE(KDRESULTS)
      DEALLOCATE(U_WIND_GFDL)
      DEALLOCATE(V_WIND_GFDL)
      DEALLOCATE(PRESS_GFDL)
      DEALLOCATE(XY_GFDL)
      DEALLOCATE(WTS)
      DEALLOCATE(LOCS)

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
       !cleans up the setMessageSource set at the beginning of the subroutine
       call unsetMessageSource()

       RETURN            
C ---------------------------------------------------------------------
       END SUBROUTINE GET_GFDL
C ---------------------------------------------------------------------

C********************************************************************

C ---------------------------------------------------------------------
C
      subroutine nws16get(timeloc,wvnx,wvny,press)
      use mesh, only : x,y,np,slam0,sfea0
      use global, only : setMessageSource,unsetMessageSource,allMessage,
     &        logMessage,scratchMessage,DEBUG,ECHO,INFO,ERROR,
     &        WARNING,PRBCKGRND,RHOWAT0,G

      implicit none
      integer :: i,numFiles,currentCycle
      logical :: cycleTimeFound
      real(8), intent(in) :: timeloc   ! adcirc time, seconds since coldstart
      real(sz), intent(out) :: wvnx(:)  ! wind u-velocity (m/s)
      real(sz), intent(out) :: wvny(:)  ! wind v-velocity (m/s)
      real(sz), intent(out) :: press(:)  !barometric pressure (mH20)

      call setMessageSource("nws16get")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      numFiles = GFDL_Files(1)%numFiles  !This is total number of files

      !
      ! locate the cycle that is relevant to current ADCIRC time
      cycleTimeFound = .false.
      do i=1,numFiles-1
         if (( GFDL_Files(i)%cycleTime.le.timeloc).and.
     &        (GFDL_Files(i+1)%cycleTime.gt.timeloc) ) then
            cycleTimeFound = .true.
            exit
         endif
      end do
      currentCycle = i

      if (cycleTimeFound.eqv. .false.) then
         if ( timeloc.lt.GFDL_Files(1)%cycleTime ) then

! ... tcm 20130801 v51.06.10 replaced error exiting with insertion of blank snaps
            write(scratchMessage,'("The current ADCIRC time ",
     &        E15.8," (days) is earlier than the first GFDL Met file time: ",
     &        E15.8, " (days). Inserting a blank wind snap.")')
     &        timeloc/86400.d0,GFDL_Files(1)%cycleTime/86400.d0
            call logMessage(WARNING,scratchMessage)
            currentCycle = 1
         endif

         if ( timeloc.gt.GFDL_Files(numFiles)%cycleTime ) then

            write(scratchMessage,'("The current ADCIRC time ",
     &        E15.8, " (days) is later than the last GFDL Met file time: ",
     &        E15.8, " (days). Inserting a blank wind snap.")')
     &         timeloc/86400.d0,GFDL_Files(numFiles)%cycleTime/86400.d0
            call logMessage(WARNING,scratchMessage)
            currentCycle = numFiles-1
         endif

      endif


      if (cycleTimeFound.eqv. .true.) then
         ! load up the data for this time interval

         write(scratchMessage,'("At ADCIRC time: ",
     &      E15.8," (days), ready to read next GFDL Met File with time: ",
     &      E15.8," (days).")') timeloc/86400.d0,
     &         GFDL_Files(currentCycle)%cycleTime/86400.d0

         call logMessage(INFO,scratchMessage)

         call GET_GFDL(NP,WVNX,WVNY,PRess,x,y,
     &                GFDL_Files(currentCycle)%file_name,
     &                SLAM0,SFEA0,
     &                GFDL_Files(currentCycle)%max_extrap_dist)
      else
         WVNX(:) = 0.d0
         WVNY(:) = 0.d0
         PRess(:) = PRBCKGRND*100.D0/(RHOWAT0*G)  !SET TO BACKGROUND PRESSURE 1013.0D0 ADJUSTED for ADCIRC units
      endif

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
        !cleans up the setMessageSource set at the beginning of the subroutine
        call unsetMessageSource()
      
        return
        end subroutine nws16get
C
C ---------------------------------------------------------------------


C***********************************************************************
C                                                                      *
C   End of subroutines to read wind and pressure fields                *
C                                                                      *
C***********************************************************************


C***********************************************************************
C                                                                      *
C   Read onto the ADCIRC grid radiation stress fields in the PBL-JAG   *
C   (hurricane) model format.                                          *
C                                                                      *
C                                                                      *
C                           R.L.05/12/99                               *
C***********************************************************************

      SUBROUTINE RSGET(RSNX,RSNY,NP)
      USE SIZES
      IMPLICIT NONE
      INTEGER NP,I,NHG
      REAL(SZ) RSNX(*),RSNY(*)
      CHARACTER*80 PBLJAGF
C
      DO I=1,NP
         RSNX(I)=0.d0
         RSNY(I)=0.d0
      END DO
 170  READ(23,'(A80)') PBLJAGF
      IF(PBLJAGF(2:2).EQ.'#') GOTO 170
 171  READ(PBLJAGF,'(I8,5E13.5)') NHG,RSNX(NHG),RSNY(NHG)
      READ(23,'(A80)') PBLJAGF
      IF(PBLJAGF(2:2).NE.'#') GOTO 171
C
      RETURN
      END SUBROUTINE


C----------------------------------------------------------------------
C...  jgf50.38.05: Subroutine to terminate the run cleanly.
C...
C----------------------------------------------------------------------
      SUBROUTINE windTerminate()
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      IMPLICIT NONE
C
      call setMessageSource("windTerminate")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      call allMessage(ERROR,"ADCIRC terminating.")
#ifdef CMPI
      call msg_fini()
#endif
      stop
C

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE windTerminate
C----------------------------------------------------------------------




      END MODULE

