C***********************************************************************
C PADCIRC VERSION 45.12 03/17/2006                                     *
C last changes in this file VERSION 45.05                              *
C                                                                      *
C***********************************************************************
C
       MODULE WIND
       USE SIZES, ONLY : LOCALDIR
       USE GLOBAL, ONLY : ScreenUnit, setMessageSource, allMessage,
     &    screenMessage, logMessage, unsetMessageSource,
     &    DEBUG, ECHO, INFO, WARNING, ERROR, scratchMessage,
     &    sphericalDistance, openFileForRead, RNDAY, cice_timinc, 
     &    cice_time1, cice_time2
#ifdef DATETIME
     &   ,basedatetime
#endif
       USE CONSTANTS, ONLY: pi, deg2rad, rad2deg, e, PRBCKGRND, g, 
     &                      rhoWat0, Rearth, rhoAir, one2ten,
     &                      windReduction
#ifdef CMPI
      use messenger, only : msg_fini
#endif
#ifdef DATETIME
      use datetime_module, only: datetime, timedelta, strptime
#endif
#if defined(DATETIME) && (defined(GRIB2API) || defined(ADCNETCDF))
      integer, parameter :: allowable_nws(21)
     &      = (/0,1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,19,20,29,30/)
#elif defined(ADCNETCDF)
      integer, parameter :: allowable_nws(20)
     &      = (/0,1,2,3,4,5,6,7,8,10,11,12,13,15,16,17,19,20,29,30/)
#else
C... v50.xx sm -- nws=17 for NUOPC coupling hwrf data
      integer, parameter :: allowable_nws(19)
     &      = (/0,1,2,3,4,5,6,7,8,10,11,12,15,16,17,19,20,29,30/)
#endif
      logical nwsOK     ! true if NWS value in fort.15 is allowed

      real(8), allocatable :: wvnx1(:) ! wind velocity in x-dir at each node at prev ts
      real(8), allocatable :: wvny1(:) ! wind velocity in y-dir at each node at prev ts
      real(8), allocatable :: wvnx2(:) ! wind velocity in x-dir at each node at new ts
      real(8), allocatable :: wvny2(:) ! wind velocity in y-dir at each node at new ts
      real(8), allocatable :: prn1(:) ! barometric presssure at each node at prev ts
      real(8), allocatable :: prn2(:) ! barometric pressure at each node at new ts
      real(8), allocatable :: w(:,:)  ! static weights for main nws13 grid

      real(8) :: wtime1 ! time associated with previous met dataset
      real(8) :: wtime2 ! time associated with new met dataset
      real(8) :: wtime1i ! intermediate NWS13 interpolation timestep
      real(8) :: wtime2i ! intermediate NWS13 interpolation timestep
      real(8) :: wtiminc ! time increment between met datasets (seconds)
      real(8) :: wtimed ! time in seconds associated with a dataset
      real(8) :: wreftim ! met reference time

Casey 110518: Added for Mark Powell's sectorial wind drag.
      REAL(8)  :: EyeLat(3) = 0.D0
      REAL(8)  :: EyeLon(3) = 0.D0
      LOGICAL  :: FoundEye
      REAL(8),allocatable :: PWeight(:,:)
      Logical,allocatable :: foundSector(:)
C
      ! related to NWS=11
      integer :: IDSETFLG
      integer :: IWTIME
      integer :: IWTIMEP
      integer :: IWYR
      integer :: NWSEGWI
      integer :: NWSGGWI

      integer  :: NWLON,NWLAT
      real(8) :: WLATMAX,WLONMIN,WLATINC,WLONINC
      integer  :: IREFYR,IREFMO,IREFDAY,IREFHR,IREFMIN

      REAL(8), PARAMETER :: TWOPI=PI*2.D0
      REAL(8), PARAMETER :: HFPI=PI/2.D0
      REAL(8), PARAMETER :: mperdeg  = Rearth * pi / 180.0d0
      real(8) :: prbckgrnd_mH2O, RhoWat0g
      real(8) :: airwaterdensityrat

      INTEGER StormNumber !jgf46.28 For Holland Wind model wind multiplier
      REAL(8) WindRefTime !jgf46.29 seconds since beginning of year, this
                          !corresponds to time=0 of the simulation
      REAL(8) BLAdj       !jgf46.32 boundary layer adjustment for Holland model

      ! jgf49.1001 Variables for use in NWS=29, for embedding an asymmetric
      ! vortex from NWS19 into an OWI basin scale met field from NAM
      ! using NWS12
      REAL(8) vortexRMW  ! Rmax of vortex data, nautical miles
      REAL(8) vortexLat  ! @center of vortex, degrees west latitude
      REAL(8) vortexLon  ! @center of vortex, degrees north longitude
      REAL(8), ALLOCATABLE :: vortexWVNX2(:) ! u wind vel from vortex model
      REAL(8), ALLOCATABLE :: vortexWVNY2(:) ! v wind vel from vortex model
      REAL(8), ALLOCATABLE :: vortexPRN2(:)  ! bar. press. from vortex model
      !
      ! radial extent of unblended vortex winds, normalized to current Rmax
      REAL(8) pureVortex
      !
      ! radial distance to point where vortex effects insignificant,
      ! normalized to current Rmax
      REAL(8) pureBackground
      !
      ! jgf50.32: The wind drag limit was originally set to 0.002d0
      ! in v46 ... the higher value of 0.0035d0 was used in ~v48 and
      ! is now the default.
      REAL(8) :: WindDragLimit = 0.0035d0   ! maximum wind drag value

      ! jgf50.32: If DragLawString is not set in the fort.15 input file,
      ! it will remain 'default'; this will cause 'Garratt' to be used
      ! in FUNCTION WindDrag, and 'IceCube' to be used in
      ! FUNCTION WindIceDrag.
      INTEGER,PARAMETER :: DRAGLAW_GARRATT = 99990
      INTEGER,PARAMETER :: DRAGLAW_POWELL  = 99991
      INTEGER,PARAMETER :: DRAGLAW_SWELL   = 99992
#ifndef POWELL
      CHARACTER(len=10) :: DragLawString = 'default' ! controls wind drag formulation
      INTEGER :: dragLawType = DRAGLAW_GARRATT
#else
      ! jgf50.32 Provide backward compatibility with hardcoded POWELL
      ! support.
#warning The flag -DPOWELL is deprecated and will be removed in a future release. Please use the &metcontrol name list functionality to avoid issues
      
      CHARACTER(len=10) :: DragLawString = 'Powell'
      INTEGER :: dragLawType = DRAGLAW_POWELL
#endif
      logical :: garatt, swell, powell
      !
      ! used for vis/diagnostics on nws19 and nws20
      LOGICAL :: writeFullCircleRmaxes = .false.
      LOGICAL :: writeRadialVandP = .false.
      LOGICAL :: writeSpatialVandP = .false.
      LOGICAL :: writeSpatialuvp = .false.
      LOGICAL :: writeNodaluvp = .false.

      !switch to turn on or off geostrophic balance in GAHM
      ! on (default): Coriolis term included, phifactors will be calculated before being used
      ! off:parameter is set to 'true', phifactors will be set to constant 1
      LOGICAL :: geostrophic_switch = .true.
      INTEGER :: geofactor = 1  !turn on or off gostrophic balance

      ! jgf50.38.05: used to keep track of HWind files and data
      type hWindData_t
         ! input data from fort.22
         character(1024)   :: file_name    ! full path to HWind file
         real(8)          :: cycleTime    ! time (sec) relative to ADCIRC cold start or REFTIM
         real(8)          :: Pc           ! storm central pressure (mb) or -1 when using p-w rel
         real(8)          :: rampVal      ! ramping multiplier for this dataset
         ! input data from hwind files themselves
         real(8)          :: dh           ! grid spacing (m) (assumed uniform)
         integer           :: nh           ! number of grid points (assumed same x and y dirs)
         real(8)          :: cLon         ! longitude of center of storm (deg)
         real(8)          :: cLat         ! latitude of center of storm (deg)
         real(8), allocatable :: u(:,:)   ! east wind velocity (m/s)
         real(8), allocatable :: v(:,:)   ! north wind velocity (m/s)
         ! derived data
         integer           :: indexLimit   ! data extent in index space (-indexLimit:indexLimit)
         real(8)          :: vmax         ! maximum wind speed in the dataset (m/s)
         real(8)          :: rmax         ! distance from storm center to vmax (nm)
         real(8)          :: rmaxAngle    ! +CCW from east (deg) angle to rmax
         logical           :: loaded       ! .true. if we have loaded the data from file
      end type hWindData_t
      type(hWindData_t), allocatable, target :: hWindFiles(:) ! array of info about HWind files
      ! jgf50.38.05: Added the following for processing HWind data
      real(8), allocatable :: wvnx_work(:) ! working array for wind u-velocity
      real(8), allocatable :: wvny_work(:) ! working array for wind v-velocity
      real(8), allocatable :: x_mercator(:) ! mesh lons in mercator, origin at storm center
      real(8), allocatable :: y_mercator(:) ! mesh lats in mercator, origin at storm center
      integer, allocatable :: x_indices(:)  ! i index of each mesh node in hwind grid
      integer, allocatable :: y_indices(:)  ! j index of each mesh node in hwind grid
      real(8), allocatable :: w1(:)        ! interp. weight for q11
      real(8), allocatable :: w2(:)        ! interp. weight for q21
      real(8), allocatable :: w3(:)        ! interp. weight for q22
      real(8), allocatable :: w4(:)        ! interp. weight for q12
      real(8), allocatable :: windSpeeds(:) ! u,v magnitude (m/s); used to find Vmax and Rmax
      logical, allocatable  :: inside(:)    ! true if a node is in the hwind bounding box
      integer               :: numFiles     ! number of hwind files listed in fort.22
      character(1024)       :: metComment   ! comment line at top of fort.22
      real(8)              :: hWindMultiplier ! all hwind u,v multiplied by this
      character(32)         :: pressureWindRelationship ! dvorak, etc
      integer               :: currentCycle ! number of most recent
                                            ! applicable hwind file
      real(8)              :: angleNow     ! +CCW from east (degrees) from storm center rmaxNow
                                            ! (deg) used by aswip to analyze radial vel

      ! tcm 51.06.02: used to keep track of GFDL Met files and data
      type gfdl_Data_t
         ! input data from fort.22
         character(1024)   :: file_name    ! full path to GFDL met files
         real(8)          :: cycleTime    ! time (sec) relative to ADCIRC cold start or REFTIM
         real(8)          :: rampVal          ! ramping multiplier for this dataset
         real(8)          :: max_extrap_dist  ! maximum distance for extrapolating values (meters)
         integer           :: numFiles         !total number of GFDL met Files to be read
         ! input data from gfdl met files themselves
         integer           :: nh           ! number of grid points
         ! derived data
         logical           :: loaded       ! .true. if we have loaded the data from file
      end type gfdl_Data_t
      type(gfdl_Data_t), allocatable, target :: gfdl_Files(:) ! array of info about GFDL Met files
      real(8) :: GFDL_max_extrap_dist     ! maximum distance for extrapolating values (meters)
      real(8) :: GFDL_WindMultiplier      ! all GFDL wind u,v's multiplied by this
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
      real(8), allocatable :: d00(:),d10(:),d11(:),d01(:)
      real(8), allocatable :: colrab(:),dummy(:),gclat(:),gclon(:)
      ! met values on gaussian grid vertices
      real(8), allocatable ::  ug(:), vg(:), pg(:), ig(:)
      integer :: nws10DatasetCounter ! indicates which dataset to read

      ! jgf52.30.09: Needed for vortex models with storm centers described
      ! in positive degrees east longitude when mesh nodes end up being
      ! defined in negative degrees east longitude
      logical :: entireMeshIsWestOfGMT ! true if all longitudes are negative

      ! jgf: variables needed for water level correction data
      ! from control file
      character(len=2048) :: dynamicWaterLevelCorrectionFileName='null' ! water level correction file
      real(8) :: dynamicWaterLevelCorrectionRampStart = 0.d0 ! seconds after run start
      real(8) :: dynamicWaterLevelCorrectionRampEnd = 0.d0   ! seconds after run start
      integer :: dynamicWaterLevelCorrectionSkipSnaps = 0 ! number of correction data sets to skip
      character(len=10) :: dynamicWaterLevelCorrectionRampReferenceTime='hotstart' ! coldstart or hotstart
      real(8) :: dynamicWaterLevelCorrectionMultiplier = 1.d0 ! multiply each correction value by this
      ! limited to this module
      integer, parameter :: dynamicWaterLevelCorrectionUn = 21 ! i/o unit number for water lvl correction
      integer :: dynamicWaterLevelCorrectionBlankSnaps = 0 ! counter for initial blank snaps
      integer :: dynamicWaterLevelCorrectionDataSetCounter=0 ! counts real datasets read from file
      logical :: dynamicWaterLevelCorrectionIsConstant=.false. ! true if only one dataset provided
      real(8), allocatable :: dynamicWaterLevelCorrectionsOld_m(:) ! from file in m
      real(8), allocatable :: dynamicWaterLevelCorrectionsNew_m(:) ! from file in m
      real(8), allocatable :: dynamicWaterLevelCorrectionFD(:) ! hold fulldomain data in parallel
      character(len=2048) :: dynamicWaterLevelCorrectionComment = 'null' ! comment from correction file
      real(8) :: dynamicWaterLevelCorrectionTimeIncrement = 0.d0 ! correction time increment (s)
      real(8) :: dynamicWaterLevelCorrectionTimeNew = 0.d0 ! time of latest dataset (s)
      real(8) :: dynamicWaterLevelCorrectionTimeOld = 0.d0 ! time of previous dataset (s)
      real(8) :: dynamicWaterLevelCorrectionDefaultValue = 0.d0 ! for otherwise unspecified nodes (m)
      integer :: dynamicWaterLevelCorrectionCurrentLine ! used in error messages
      logical :: dynamicWaterLevelCorrectionFileEnded = .false. ! true when all data have been read
#ifdef DATETIME
      integer :: ub, ubc ! size of weightsp and indp interpolant weights
      real(8),allocatable :: weightsp(:,:,:) 
      integer,allocatable  :: indp(:,:,:) 
      character(len=200) :: Pfile, Wfile, Pinv, Winv, Cinv, Cfile
      character(len=200) :: Wfile1, Pinv1, Winv1
      INTEGER :: PfileNCID, WfileNCID, CfileNCID, Wfile1NCID
      character(len=200) :: Pvar, Uvar, Vvar, Cvar, Tvar, Lonvar, Latvar
      character(len=200) :: Tdim, Londim, Latdim, Tformat
      type(datetime) :: CurDT, refdate, stepdate
      logical :: grb2flag
      integer :: NT, NTC
#endif
      real(8) :: wtime1_12 ! time associated with previous met dataset
      real(8) :: wtime2_12 ! time associated with new met dataset
      real(8) :: wtiminc_12 ! time increment between met datasets (seconds)
      logical  :: change_14 ! true when the NWS = 14 data has been updated
      logical,allocatable :: hasdata(:) ! if OWI data exists at this vertex       
      !...Function interface for wind drag
      !   Used so that the selection of wind drag
      !   is only performed once. Each wind drag
      !   function has the same function signiture
      ABSTRACT INTERFACE
        FUNCTION p_computeWindDrag(WindSpeed,NodeNumber) RESULT(Drag)
            REAL(8),INTENT(IN)  :: WindSpeed
            INTEGER,INTENT(IN)  :: NodeNumber
            REAL(8)             :: Drag
        END FUNCTION p_computeWindDrag
      END INTERFACE
      PROCEDURE(p_computeWindDrag),POINTER :: f_computeWindDrag => null()

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
C     jgf50.32 Subroutine to initialize meteorological parameters to their
C     default values.
C     jgf53.dev: Now allocates memory, whereas default values are now set
C     in the module variable declaration section above.
C     This subroutine must be called after read_input() so that it can
C     use the values of np and g.
C     ----------------------------------------------------------------
      subroutine initWindModule()
      use mesh, only : np
      use global, only: NWS
      implicit none
      
      call setMessageSource("initWindModule")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      allocate ( wvnx1(np),wvny1(np) )
      allocate ( wvnx2(np),wvny2(np) )
      allocate ( prn1(np) )
      allocate ( prn2(np) )

      if (NWS.eq.-14) allocate ( HasData(np) )
C     WJP: this is the value used in the surface stress computation  
C     (can contol through the rhoAir parameter in namelist)
      airwaterdensityrat = rhoAir / rhoWat0

C.... tcm v49.16 20100617 added
C.... convert background pressure from millibars to meters of water
      rhoWat0g = rhoWat0 * g
      PRBCKGRND_MH2O = 100.0D0 * PRBCKGRND / rhoWat0g
      prn1 = PRBCKGRND_MH2O
      prn2 = PRBCKGRND_MH2O

      ! dw
      wvnx1 = 0.0
      wvny1 = 0.0
      wvnx2 = 0.0
      wvny2 = 0.0

      CALL checkWindDragType()
      CALL mapWindDragFunctionPointer()

      if (NWS.eq.0) return

      ! Make logicals here if NWS is not 0 
      garatt = .false.; swell = .false.; powell = .false.
      SELECT CASE(TRIM(DragLawString))
         CASE("Garratt","GARRATT","garratt","default")
             garatt = .true.
             write(16,*) 'Garratt drag law used' 
#ifdef CSWAN
         CASE("Swell","SWELL","swell")
             swell   = .true.
             write(16,*) 'Swell drag law used' 
#endif
         CASE("Powell","POWELL","powell")
             if (NWS.eq.12.or.NWS.eq.29.or.NWS.eq.30) then
                powell  = .true.
                allocate(Pweight(3,np)); Pweight = 0d0
                allocate(foundSector(np)); foundSector = .false.
                write(16,*) 'Powell drag law used'
             else
                garatt = .true.
                write(16,*) 'Garratt drag law used because NWS ', 
     &                      'incompatible with Powell'
             endif
         CASE DEFAULT
            WRITE(16,*) 'ERROR: Wind drag law not recognized:'
            WRITE(16,'(A10)') DragLawString
            WRITE(16,*) 'Execution will now be terminated.'
#ifdef CMPI
            call msg_fini()
#endif
            CALL EXIT(1)
C
      END SELECT

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

!     ----------------------------------------------------------------
      END SUBROUTINE initWindModule
!     ----------------------------------------------------------------

C----------------------------------------------------------------------
C     Sets an integer for the drag law to be used later in checks
C     This is simpler than checking a string
C----------------------------------------------------------------------
      SUBROUTINE checkWindDragType()
C----------------------------------------------------------------------
          USE GLOBAL,ONLY: toLowercase
          IMPLICIT NONE
          CHARACTER(LEN(dragLawString)) :: dragLawLowercase
          draglawlowercase = toLowercase(dragLawString)
          SELECT CASE(TRIM(draglawlowercase))
          CASE("default","garratt")
            draglawtype = DRAGLAW_GARRATT
          CASE("powell")
            draglawtype = DRAGLAW_POWELL
          CASE("swell")
            draglawtype = DRAGLAW_SWELL
          CASE DEFAULT
            draglawtype = DRAGLAW_GARRATT
          END SELECT
C----------------------------------------------------------------------
      END SUBROUTINE checkWindDragType
C----------------------------------------------------------------------

!     ----------------------------------------------------------------
!>    @brief Maps the function pointer used in the calculation of
!>    wind drag. This is used so that the logical statements required
!>    to select the wind drag are only executed once at startup
!     ----------------------------------------------------------------
      SUBROUTINE mapWindDragFunctionPointer()
!     ----------------------------------------------------------------
          IMPLICIT NONE
#ifdef POWELL
          f_computeWindDrag => PowellWindDrag
#else
          SELECT CASE(DragLawType)
          CASE(DRAGLAW_GARRATT)
              f_computeWindDrag => CappedGarrattWindDrag
          CASE(DRAGLAW_POWELL)
              f_computeWindDrag => PowellWindDrag
          CASE(DRAGLAW_SWELL)
#ifndef CSWAN
              call allMessage(WARNING,"SWAN Model not enabled."//
     &              "Drag formulation will be Garratt.")
              f_computeWindDrag => CappedGarrattWindDrag
#else
              f_computeWindDrag => SwanWindDrag
#endif
          CASE DEFAULT
              call allMessage(ERROR,"Wind drag law not recognized"//TRIM(DragLawString))
#ifdef CMPI
              call msg_fini()
#endif
              CALL EXIT(1)
          END SELECT
#endif
!     ----------------------------------------------------------------
      END SUBROUTINE mapWindDragFunctionPointer
!     ----------------------------------------------------------------

!     ---------------------------------------------------------------
!>    @brief Computes wind drag using a function with the signiture:
!>    result = function(windSpeed,nodeNumber)
!>    @param[in] WindSpeed wind magnitude
!>    @param[in] NodeNumber adcirc mesh node index
!>    @result wind drag coefficient
!     ----------------------------------------------------------------
      REAL(8) FUNCTION WindDrag(WindSpeed, NodeNumber)
        USE GLOBAL,ONLY: windDragOut
        IMPLICIT NONE
        REAL(8) :: WindSpeed
        INTEGER  :: NodeNumber
        WindDrag = f_computeWindDrag(windSpeed,NodeNumber)
        windDragOut(nodeNumber) = windDrag
        RETURN
!     ----------------------------------------------------------------
      END FUNCTION WindDrag
!     ----------------------------------------------------------------

#ifdef CSWAN
!     ----------------------------------------------------------------
!>    @brief Use the wind drag coefficients calculated inside SWAN
!>    for wind drag inside of adcirc
!>    @param[in] WindSpeed wind magnitude
!>    @param[in] NodeNumber adcirc mesh node index
!>    @result wind drag coefficient
!     ----------------------------------------------------------------
       FUNCTION SwanWindDrag(WindSpeed,NodeNumber) RESULT(Drag)
!     ----------------------------------------------------------------
           USE GLOBAL,ONLY: Swan_WDragCo
           IMPLICIT NONE
           REAL(8),INTENT(IN)  :: WindSpeed
           INTEGER,INTENT(IN)  :: NodeNumber
           REAL(8)             :: Drag
           Drag = SWAN_WDragCo(NodeNumber)
!     ----------------------------------------------------------------
       END FUNCTION SwanWindDrag
!     ----------------------------------------------------------------
#endif

!     ----------------------------------------------------------------
!>    @brief Function applying the Garratt wind drag formulation.
!>
!>    Wind drag is calculated using Garratt, 1977 by:
!>    \f[
!>       d = \frac{\left( 0.75 + 0.067 \left| \vec{w} \right| \right)}{1000}
!>    \f]
!>
!>    @param[in] WindSpeed wind magnitude
!>    @param[in] NodeNumber adcirc mesh node index
!>    @result wind drag coefficient
!     ----------------------------------------------------------------
       FUNCTION GarrattWindDrag(WindSpeed,NodeNumber) RESULT(Drag)
!     ----------------------------------------------------------------
           IMPLICIT NONE
           REAL(8),INTENT(IN)  :: WindSpeed
           INTEGER,INTENT(IN)  :: NodeNumber
           REAL(8)             :: Drag
           Drag = 0.001d0*(0.75d0+0.067d0*WindSpeed)
       END FUNCTION GarrattWindDrag

!     ----------------------------------------------------------------
!>    @brief Applies the Garratt Wind Drag capped at WindDragLimit.
!>    By default the cap applied is 0.0035, however this can be adjusted
!>    using the metControl namelist parameter WindDragLimit
!>    @param[in] WindSpeed wind magnitude
!>    @param[in] NodeNumber adcirc mesh node index
!>    @result wind drag coefficient
!     ----------------------------------------------------------------
       FUNCTION CappedGarrattWindDrag(WindSpeed,NodeNumber) RESULT(WindDrag)
!     ----------------------------------------------------------------
           IMPLICIT NONE
           REAL(8),INTENT(IN)  :: WindSpeed
           INTEGER,INTENT(IN)  :: NodeNumber
           REAL(8)             :: WindDrag
           WindDrag = MIN(GarrattWindDrag(WindSpeed,NodeNumber),WindDragLimit)
!     ----------------------------------------------------------------
       END FUNCTION CappedGarrattWindDrag
!     ----------------------------------------------------------------

!     ----------------------------------------------------------------
!>    @brief Applies Mark Powell's sector-based wind drag. This
!>    applies different wind drag values based upon the sector of the
!>    storm a node lies within. The is done to account for the direction
!>    of wave travel and the coincident direction of the wind vector
!>    @param[in] WindSpeed wind magnitude
!>    @param[in] NodeNumber adcirc mesh node index
!>    @result wind drag coefficient
!     ----------------------------------------------------------------
      FUNCTION PowellWindDrag(WindSpeed,NodeNumber) RESULT(PowellDrag)
!     ----------------------------------------------------------------
          USE GLOBAL, ONLY: DEG2RAD
          USE MESH, ONLY : SFEA, SLAM
          IMPLICIT NONE
          INTRINSIC           :: ATAN2
          INTRINSIC           :: MOD
          REAL(8),INTENT(IN)  :: WindSpeed
          INTEGER,INTENT(IN)  :: NodeNumber
          INTEGER             :: IS
          INTEGER             :: IW1,IW2,IW3
          LOGICAL             :: FoundSector
          REAL(8)             :: Dir1
          REAL(8)             :: Dir2
          REAL(8)             :: Drag(3)
          REAL(8)             :: NodeDirection
          REAL(8)             :: NodeLat
          REAL(8)             :: NodeLon
          REAL(8)             :: StormDirection
          REAL(8)             :: Weight(3)
          REAL(8)             :: PowellDrag
          REAL(8)             :: rawGarrattWindDrag

!.. Check whether we have previously found the eye.
         IF((EyeLat(1).EQ.0.D0).OR.(EyeLon(1).EQ.0.D0).OR.
     &      (EyeLat(2).EQ.0.D0).OR.(EyeLon(2).EQ.0.D0).OR.
     &      (EyeLat(3).EQ.0.D0).OR.(EyeLon(3).EQ.0.D0))THEN
            FoundEye = .FALSE.
         ENDIF
!.. If the eye has not been found then apply Garratt.
         IF(.NOT.FoundEye)THEN
            PowellDrag = CappedGarrattWindDrag(WindSpeed,NodeNumber)
            FoundSector = .FALSE.
            Weight = 0.D0
         ELSE
!.. We have found the eye!
!.. Compute the storm's direction.  We use the convention of northward
!.. being zero degrees, and then cycling clockwise.
            Dir1 = ATAN2(EyeLat(2)-EyeLat(1),EyeLon(2)-EyeLon(1))
            Dir2 = ATAN2(EyeLat(3)-EyeLat(2),EyeLon(3)-EyeLon(2))
            StormDirection = ATAN2(SIN(Dir1)+SIN(Dir2),COS(Dir1)+COS(Dir2))
            StormDirection = StormDirection/DEG2RAD
            StormDirection = 90.D0 - StormDirection
            IF(StormDirection.LT.0.D0)THEN
               StormDirection = StormDirection + 360.D0
            ENDIF
!.. Compute the direction of node relative to the eye.
            NodeLon = SLAM(NodeNumber)/DEG2RAD
            NodeLat = SFEA(NodeNumber)/DEG2RAD
            NodeDirection = ATAN2(NodeLat-EyeLat(3),NodeLon-EyeLon(3))
            NodeDirection = NodeDirection/DEG2RAD
            NodeDirection = 90.D0 - NodeDirection
            IF(NodeDirection.LT.0.D0)THEN
               NodeDirection = NodeDirection + 360.D0
            ENDIF
!.. Compute the weights for the sector in which the node is located.
!..     0- 40 : Transition from left (3) to right (1).
!..    40-130 : Right (1).
!..   130-170 : Transition from right (1) to rear (2).
!..   170-220 : Rear (2).
!..   220-260 : Transition from rear (2) to left (3).
!..   260-  0 : Left (3).
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


!.. Apply garratt wind drag in case of emergency
            IF(.NOT.FoundSector)THEN
               PowellDrag = CappedGarrattWindDrag(WindSpeed,NodeNumber)
               Weight = 0.D0
               RETURN
            ELSE
               rawGarrattWindDrag = GarrattWindDrag(WindSpeed,NodeNumber)
!.. Determine the wind drag for sector 1 (right).
               Drag(1) = rawGarrattWindDrag
               IF(WindSpeed.LE.35.0D0)THEN
                  IF(Drag(1).GT.0.0020D0)THEN
                     Drag(1) = 0.0020D0
                  ENDIF
               ELSEIF(WindSpeed.LE.45.0D0)THEN
                  Drag(1) = 0.0020D0 + (0.0030D0-0.0020D0)/(45.0D0-35.0D0)*(WindSpeed-35.0D0)
               ELSE
                  Drag(1) = 0.0030D0
               ENDIF
!.. Determine the wind drag for sector 2 (rear).
               Drag(2) = rawGarrattWindDrag
               IF(WindSpeed.LE.35.0D0)THEN
                  IF(Drag(2).GT.0.0020D0)THEN
                     Drag(2) = 0.0020D0
                  ENDIF
               ELSEIF(WindSpeed.LE.45.0D0)THEN
                  Drag(2) = 0.0020D0 + (0.0010D0-0.0020D0)/(45.0D0-35.0D0)*(WindSpeed-35.0D0)
               ELSE
                  Drag(2) = 0.0010D0
               ENDIF
!.. Determine the wind drag for sector 3 (left).
               Drag(3) = rawGarrattWindDrag
               IF(Drag(3).GT.0.0018D0)THEN
                  IF(WindSpeed.LE.25.0D0)THEN
                     Drag(3) = 0.0018D0
                  ELSEIF(WindSpeed.LE.30.0D0)THEN
                     Drag(3) = 0.0018D0 + (0.0045D0-0.0018D0)/(30.0D0-25.0D0)*(WindSpeed-25.0D0)
                  ELSEIF(WindSpeed.LE.45.0D0)THEN
                     Drag(3) = 0.0045D0 + (0.0010D0-0.0045D0)/(45.0D0-30.0D0)*(WindSpeed-30.0D0)
                  ELSE
                     Drag(3) = 0.0010D0
                  ENDIF
               ENDIF
!.. Apply a weighted average.
               PowellDrag = Weight(1)*Drag(1) + Weight(2)*Drag(2) + Weight(3)*Drag(3)
            ENDIF
         ENDIF
!     ----------------------------------------------------------------
      END FUNCTION PowellWindDrag
!     ----------------------------------------------------------------

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
      REAL(8) FUNCTION WindIceDrag(WindDrag,PercentIce)
      implicit none
      REAL(8) :: PercentIce
      REAL(8) :: pic, WindDrag, IceDrag, Cform
      REAL(8) :: Cskin = 1.5d-3, Cform_max = 3.67d-3, beta = 0.6d0
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

      CASE("IceCube")
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
      
      CASE("Lupkes","default","Powell")
C........This formula is a function that explicitly
C....... combines air-sea drag (Garratt) with sea-ice skin drag and
C....... sea-ice form drag components with area fraction weights 
C....... see: Lüpkes, et al. (2012), JGR: doi:10.1029/2012JD017630. &
C....... Joyce et al. (2019), OM: doi:10.1016/j.ocemod.2019.101421
         IF (PercentIce.lt.0.d0) then
            WindIceDrag = WindDrag
         ELSE
            pic = min(PercentIce*0.01d0,1d0)
            Cform = Cform_max*(1d0-pic)**beta
            WindIceDrag = WindDrag*(1d0 - pic) + (Cskin + Cform)*pic
         ENDIF

      CASE DEFAULT
         WRITE(16,*) 'ERROR: Ice drag law not recognized:'
         WRITE(16,'(A10)') DragLawString
         WRITE(16,*) 'Execution will now be terminated.'
#ifdef CMPI
         call msg_fini()
#endif
         CALL EXIT(1)
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
      SUBROUTINE getBlendFactor(i, bf)
      USE MESH, ONLY : SLAM, SFEA
      IMPLICIT NONE
      INTEGER, intent(in) :: i  ! node number
      REAL(8), intent(out):: bf       ! blend factor, 0 to 1
      REAL(8) XCOOR ! longitude of node in degrees
      REAL(8) YCOOR ! latitude of node in degrees
C     distance from center of vortex to mesh node:
      REAL(8) DX    ! longitude distance in radians
      REAL(8) DY    ! latitude distance in radians
      REAL(8) rad   ! point-to-point distance over a sphere, in meters
      REAL(8) vortexOnly ! radial distance to edge of pure vortex, meters
      REAL(8) backgroundOnly ! radial dist to undisturbed background, meters
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
         bf = (backgroundOnly-rad)/(backgroundOnly-vortexOnly)
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
      SUBROUTINE TIMECONV(IYR,IMO,IDAY,IHR,IMIN,SEC,TIMESEC)
      IMPLICIT NONE
      INTEGER IYR,IMO,IDAY,IHR,IMIN,ILEAP
      REAL*8 TIMESEC,SEC
C
      call setMessageSource("timeConv")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
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
         call allMessage(ERROR,
     &     'Fatal error in subroutine timeconv: month > 12.')
         call windTerminate()

      ENDIF
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE TIMECONV
C-----------------------------------------------------------------------

C----------------------------------------------------------------------
C                    S U B R O U T I N E
C       G E T   M E T E O R O L O G I C A L   F O R C I N G
C----------------------------------------------------------------------
C     Called from timestep() to perform i/o, time interpolation,
C     wind reduction, canopy effects, drag modification due to ice,
C     computation of stress from velocity, units conversion, and
C     any other procedures necessary for producing met data at a
C     particular point in time.
C----------------------------------------------------------------------
      subroutine getMeteorologicalForcing(nws, wsx2, wsy2, pr2,
     &   wvnxout, wvnyout, timeloc)
      use owiwind,only : nws12get, findStormCenter    ! sb added 09/xx/2006     
      use global, only : ciceout, h0, eta2, ncice, rampmete,
     &   sponge, res_bc_flag, IDEN, SigTS, cice1, cice2, IFSPROTS
      use constants, only: g
      use nodalAttributes, only : applyDirectionalWindReduction,
     &   loadCanopyCoef, loadDirEffRLen, applyCanopyCoefficient
      use mesh, only : dp, ics, np, uvectmp, vvectmp, DRVMAP2DSPVEC
#ifdef ADCNETCDF
      USE OWIWIND_NETCDF, ONLY: NWS13GET
#endif
      implicit none
      integer, intent(in) :: nws
      real(8), intent(out) :: wsx2(:)
      real(8), intent(out) :: wsy2(:)
      real(8), intent(out) :: pr2(:)
      real(8), intent(out) :: wvnxout(:)
      real(8), intent(out) :: wvnyout(:)
      real(8), intent(in) :: timeloc
      real(8) :: windx ! wind vel (m/s) in x-dir
      real(8) :: windy ! wind vel (m/s) in y-dir

      real(8) :: prdiff
      real(8) :: vortexPRDIFF
      real(8) :: wDragCo
      real(8) :: windMag
      real(8) :: wtratio ! used to time interpolated between met datasets
      real(8) :: cicenow, icedragco, ctratio, icemag, icex, icey 
      integer :: nhg ! node number temp
      real(8) :: bf ! blend factor for multiple sets of wind data
      integer :: i ! node loop counter
      real(8) :: wtratio_12, wtratio_14 ! used to time interpolated between met datasets

      call setMessageSource("getMeteorologicalForcing")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif


C------------------------MET FORCING---------------------------------------
C
C...  UPDATE THE WIND STRESS AND SURFACE PRESSURE AND READ IN NEW VALUES
C...  FROM UNIT 22.  APPLY Ramp FUNCTION.
C
C  tcm v49.16 20100617 -- Changed pressure ramping so that we apply the ramp
C     to the difference between the background pressure and the forced pressure,
C     then add that ramped difference back to the background pressure
C
C     No wind, radiation stress or atmospheric pressure forcings are used.
      IF(NWS.EQ.1) THEN
         DO I=1,NP
            READ(22,*) NHG,WSX2(I),WSY2(I),PR2(I) !read in
            WSX2(I)=RampMete*WSX2(I) ! apply ramp function
            WSY2(I)=RampMete*WSY2(I)
C  tcm v49.16 20100617
!            PR2(I)=RampMete*PR2(I)
            PRDIFF = RampMete*(PR2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            wvnxout(i)=WSX2(i) ! for met recording station output
            wvnyout(i)=WSY2(i)
         END DO
      ENDIF

C     Wind stress and atmospheric pressure are read in at all grid nodes
C     at a time interval that does not equal the model time
C     step. Interpolation in time is used to synchronize the wind and
C     pressure information with the model time step.
      IF(ABS(NWS).EQ.2) THEN
         IF(TimeLoc.GT.WTIME2) THEN !determine if met file time incr. is exceeded
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
               WVNX1(I)=WVNX2(I) ! move current data to old
               WVNY1(I)=WVNY2(I)
               PRN1(I)=PRN2(I)
               READ(22,*) NHG,WVNX2(I),WVNY2(I),PRN2(I) ! read in
            END DO
         ENDIF
         WTRatio=(TimeLoc-WTIME1)/WTIMINC ! interpolate
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WSX2(I) = RampMete*WindX !apply ramp
            WSY2(I) = RampMete*WindY
C  tcm v49.16 20100617
!            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                          -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            wvnxout(i)=WSX2(i) !for met recording sta. output
            wvnyout(i)=WSY2(i)
         END DO
      ENDIF

C     Wind velocity in US Navy Fleet Numeric format interpolated in
C     space onto the ADCIRC grid and in time to synchronize the wind and
C     pressure information with the model time step. Garratt's formula
C     is used to compute wind stress from the wind velocity.
      IF(NWS.EQ.3) THEN
         IF(TimeLoc.GT.WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
            END DO
            CALL NWS3GET()
         ENDIF
         WTRatio=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete*airwaterdensityrat*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*airwaterdensityrat*WDragCo*WindY*WindMag
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
         END DO
      ENDIF

C     Wind velocity and atmospheric pressure are read in (PBL/JAG
C     format) at selected ADCIRC grid nodes. Interpolation in time is
C     used to synchronize the wind and pressure information with the
C     model time step. Garratt's formula is used to compute wind stress
C     from wind velocity.
      IF(ABS(NWS).EQ.4) THEN
         IF(TimeLoc.GT.WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
               PRN1(I)=PRN2(I)
            END DO
            CALL NWS4GET()
         ENDIF
         WTRatio=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete*airwaterdensityrat*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*airwaterdensityrat*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
!            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
         END DO
      ENDIF

C     Wind velocity and atmospheric pressure are read in at all grid
C     nodes. Interpolation in time is used to synchronize the wind and
C     pressure information with the model time step. Garratt's formula
C     is used to compute wind stress from wind velocity.
      IF(ABS(NWS).EQ.5) THEN
         IF(TimeLoc.GT.WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
               PRN1(I)=PRN2(I)
               READ(22,*) NHG,WVNX2(I),WVNY2(I),PRN2(I)
            END DO
         ENDIF
         WTRatio=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete*airwaterdensityrat*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*airwaterdensityrat*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
!            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
         END DO
      ENDIF

C     Wind velocity and atmospheric pressure are read in for a
C     rectangular grid (either in Longitude, Latitude or Cartesian
C     coordinates, consistent with the grid coordinates) and
C     interpolated in space onto the ADCIRC grid and in time to
C     synchronize the wind and pressure information with the model time
C     step. Garratt's formula is used to compute wind stress from the
C     wind velocity.
      IF(NWS.EQ.6) THEN
         IF(TimeLoc.GT.WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
               PRN1(I)=PRN2(I)
            END DO
            NWSGGWI=NWSGGWI+1
            CALL NWS6GET()
         ENDIF
         WTRatio=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete*airwaterdensityrat*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*airwaterdensityrat*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
!            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
         END DO
      ENDIF

C     jgf46.01 New option to read in surface wind stress and atmospheric
C     pressure for a rectangular grid (either in Longitude, Latitude or
C     Cartesian coordinates, consistent with the grid coordinates) and
C     interpolate in space onto the ADCIRC grid. Interpolation in time
C     is used to synchronize the wind and pressure information with the
C     model time step.
      IF(ABS(NWS).EQ.7) THEN
         IF(TimeLoc.GT.WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
               PRN1(I)=PRN2(I)
            END DO
            CALL NWS7GET()
         ENDIF
         WTRatio=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WSX2(I) = RampMete*WindX !apply ramp
            WSY2(I) = RampMete*WindY
C  tcm v49.16 20100617
!            PR2(I)=RampMete*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            wvnxout(i)=WSX2(i) !for met recording sta. output
            wvnyout(i)=WSY2(i)
         END DO
      ENDIF

C     jgf46.02 New option to read in hurricane locations and generate
C     generate hurricane winds from the Holland Wind Model.
      IF(ABS(NWS).EQ.8) THEN
         call HollandGet(wvnx2,wvny2,prn2,timeloc)
         DO I=1,NP
            WindX = WVNX2(I)
            WindY = WVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete*airwaterdensityrat*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*airwaterdensityrat*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
!            PR2(I)=RampMete*PRN2(I)
            PRDIFF = RampMete*(PRN2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
         ENDDO
      ENDIF
C
C     Wind velocity (10 m) and atmospheric pressure are read in from a
C     sequence of National Weather Service (NWS) Aviation (AVN) model
C     output files. Each AVN file is assumed to contain data on a
C     Gaussian longitude, latitude grid at a single time. Consecutive
C     files in the sequence are separated by N hours in time. Garratt's
C     formula is used to compute wind stress from the wind velocity.
      IF(NWS.EQ.10) THEN
         IF(TimeLoc.GT.WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
               PRN1(I)=PRN2(I)
            END DO
            IF (NCICE.eq.10) THEN 
               cice1 = cice2
               call nws10get(wvnx2,wvny2,prn2,cice2) 
            ELSE
               call nws10get(wvnx2,wvny2,prn2) 
            ENDIF
         ENDIF
         WTRatio=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF (NCICE.eq.10) then
               Cicenow = 100d0*( CICE1(I) + 
     &                           WTRatio*( CICE2(I)-CICE1(I) ) )
               WDragCo = WindIceDrag(WDragCo,Cicenow)
            ELSEIF (NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete*airwaterdensityrat*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*airwaterdensityrat*WDragCo*WindY*WindMag
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                          -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
         END DO
      ENDIF
C
C     Wind velocity (10 m) and atmospheric pressure are read in from a
C     sequence of stripped down National Weather Service (NWS) ETA 29km
C     model output files. Each ETA file is assumed to contain data on an
C     E grid for a single day (8 data sets, one every 3 hours, beginning
C     @ 03:00 and continuing through 24:00 of the given day). The wind
C     data is converted to an east-west, north-south coordinate system
C     inside ADCIRC. Garratt's formula is used to compute wind stress
C     from the wind velocity.
      IF(NWS.EQ.11) THEN
         IF(TimeLoc.GT.WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
               PRN1(I)=PRN2(I)
            END DO
            IDSETFLG=IDSETFLG+1
            IF(IDSETFLG.GT.8) THEN
               NWSEGWI=NWSEGWI+1
               IDSETFLG=1
            ENDIF
            CALL NWS11GET()
         ENDIF
         WTRatio=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete*airwaterdensityrat*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*airwaterdensityrat*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                          -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
         END DO
      ENDIF

C...  sb46.28sb01 NWS=12 reads in raw OWI files 09/xx/2006
C     WJP: Added for NWS = 14 (GRIB2 e.g. CFS)
      IF (ABS(NWS).EQ.12.or.ABS(NWS).eq.14) THEN
         IF (TimeLoc.GT.WTIME2) THEN
            CHANGE_14 = .true.
            WTIME1 = WTIME2
            WTIME2 = WTIME2 + WTIMINC
            DO I = 1,NP
              WVNX1(I) = WVNX2(I)
              WVNY1(I) = WVNY2(I)
              PRN1(I)  = PRN2(I)
            END DO
            IF (ABS(NWS).eq.12) then
               CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,PRBCKGRND)
#ifdef DATETIME
#if defined(GRIB2API) || defined(ADCNETCDF)
            ELSE
               IF ((NCICE.eq.14.or.NCICE.eq.13).and.
     &             TimeLoc.GT.CICE_TIME2) THEN
                  CICE_TIME1 = CICE_TIME2
                  CICE_TIME2 = CICE_TIME2 + CICE_TIMINC
                  CICE1 = CICE2
                  CALL NWS14GET(WVNX2,WVNY2,PRN2,CICE2)
               ELSE
                  CALL NWS14GET(WVNX2,WVNY2,PRN2) 
               ENDIF
#endif
#endif
            ENDIF
            if (draglawtype.EQ.draglaw_powell) then
               call findStormCenter(eyeLat, eyeLon, foundEye)
            endif
         ELSE
            CHANGE_14 = .false.
         ENDIF
         WTRATIO = (TimeLoc-WTIME1)/WTIMINC
         IF(NCICE.NE.0)THEN
             CTRATIO = (TimeLoc-CICE_TIME1)/CICE_TIMINC
         ENDIF
         IF (NWS.eq.-14) then
            IF (TimeLoc.GT.WTIME2_12) THEN
               ! Now get OWI file data (will overwrite where data
               ! exists)
               WTIME1_12 = WTIME2_12
               WTIME2_12 = WTIME2_12 + WTIMINC_12
               ! Don't do this is we have already for NWS = 14 data
               ! on same time snap
               IF (.not.CHANGE_14) THEN
                  DO I = 1,NP
                     IF (.not.HasData(I)) cycle
                     WVNX1(I) = WVNX2(I)
                     WVNY1(I) = WVNY2(I)
                     PRN1(I)  = PRN2(I)
                  END DO
               ENDIF
               CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,PRBCKGRND,HasData)
               if (draglawtype.EQ.draglaw_powell) then
                  call findStormCenter(eyeLat, eyeLon, foundEye)
               endif
            ENDIF
            WTRATIO_12 = (TimeLoc-WTIME1_12)/WTIMINC_12
            WTRATIO_14 = WTRATIO
         ENDIF
         DO I = 1,NP
            IF (NWS.eq.-14) then
               WTRATIO = WTRATIO_14
               IF (HasData(I)) WTRATIO = WTRATIO_12
            ENDIF
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            !Casey 110518: Enable Mark Powell's sector-based wind drag.
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
              CALL ApplyDirectionalWindReduction(I, WDragCo,
     &             WindMag, DP(I), ETA2(I), WindX, WindY)
              WindMag = SQRT(WindX*WindX+WindY*WindY)
              !Casey 110518: Enable Mark Powell's sector-based wind drag.
              WDragCo = WindDrag(WindMag,I)
            ENDIF
            IF (LoadCanopyCoef) then
               CALL ApplyCanopyCoefficient(I,WindX,WindY)
            ENDIF
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF (NCICE.eq.14) then
               Cicenow = 100d0*(CICE1(I) + 
     &                          CTRATIO*(CICE2(I)-CICE1(I)))
               WDragCo = WindIceDrag(WDragCo,Cicenow)
               ciceout(i) = Cicenow ! GML added 20210727
            elseif (NCICE.ne.0.and.NCICE.ne.13) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            endif
            ! WJP: added varying sea surface density
            if (abs(IDEN).ge.5) then
               if (SigTS(I) < -1d2) then
                  airwaterdensityrat = rhoAir/RhoWat0
               else
                  airwaterdensityrat = rhoAir/(SigTS(I) + RhoWat0)
               endif
            endif
            WSX2(I) = RampMete*airwaterdensityrat*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete*airwaterdensityrat*WDRAGCO*WINDY*WINDMAG
            PRDIFF  = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                          -PRBCKGRND_MH2O)
#ifdef NOIVB
            PR2(I)  = PRBCKGRND_MH2O
#else
            PR2(I)  = PRBCKGRND_MH2O + PRDIFF
#endif
            WVNXOUT(I) = RampMete*WINDX
            WVNYOUT(I) = RampMete*WINDY

          END DO
      ENDIF

#ifdef ADCNETCDF
Casey 180318: Added for NWS=13, ARC modified 190108
        IF(ABS(NWS).EQ.13) THEN
          IF(TimeLoc.GT.WTIME2I) THEN
            WTIME1I=WTIME2I
            WTIME2I=WTIME2I+WTIMINC
            DO I=1,NP
              WVNX1(I)=WVNX2(I)
              WVNY1(I)=WVNY2(I)
              PRN1(I)=PRN2(I)
            END DO
            CALL NWS13GET(WTIME2I,WVNX2,WVNY2,PRN2,
     &          WTIME2,EyeLon,EyeLat,FoundEye,W)
!            write(6,*) "timeloc", timeloc, "wtiminc", wtiminc, "inside",
!     &        wtime1i, wtime2i, wtime1, wtime2
          ENDIF
          WTRATIO=(TimeLoc-WTIME1I)/WTIMINC
          DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
              CALL ApplyDirectionalWindReduction(I, WDragCo,
     &             WindMag, DP(I), ETA2(I), WindX, WindY)
              WindMag = SQRT(WindX*WindX+WindY*WindY)
              WDragCo = WindDrag(WindMag,I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete*0.001293d0*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete*0.001293d0*WDRAGCO*WINDY*WINDMAG
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                          -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WINDX
            WVNYOUT(I)=RampMete*WINDY
          ENDDO
        ENDIF
#endif

C     RJW ffpl Merged:
!     rjw added nws = 19: asymmetric hurricane winds
      IF(NWS.EQ.19) THEN
         CALL NWS19GET(wvnx2, wvny2, prn2, timeloc)
         DO I=1,NP
            WindX = WVNX2(I)
            WindY = WVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
            WSX2(I) = RampMete*airwaterdensityrat*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*airwaterdensityrat*WDragCo*WindY*WindMag
C  tcm v49.16 20100617C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
!            PR2(I)=RampMete*PRN2(I)
            PRDIFF = RampMete*(PRN2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
         END DO
      ENDIF
C
C     jie added nws = 20: generalized asymmetric vortex model
      IF(NWS.EQ.20) THEN
         CALL NWS20GET(wvnx2, wvny2, prn2, timeloc)
         DO I=1,NP
            WindX = WVNX2(I)
            WindY = WVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
            WSX2(I) = RampMete*airwaterdensityrat*WDragCo*WindX*WindMag
            WSY2(I) = RampMete*airwaterdensityrat*WDragCo*WindY*WindMag
C  tcm v49.16 20100617C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
!            PR2(I)=RampMete*PRN2(I)
            PRDIFF = RampMete*(PRN2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
         END DO
      ENDIF
C
C     jgf49.1001 Added NWS29 for embedding an asymmetric vortex from
C     NWS19 into an OWI basin scale met field from NWS12 (derived from NAM).
      IF((ABS(NWS).EQ.29).or.(ABS(NWS).eq.30)) THEN
         ! bring in next set of OWI met data
         IF(TimeLoc.GT.WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
               WVNX1(I)=WVNX2(I)
               WVNY1(I)=WVNY2(I)
               PRN1(I)=PRN2(I)
            END DO
            CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,PRBCKGRND)
            if(draglawtype.EQ.draglaw_powell)then
               call findStormCenter(eyeLat, eyeLon, foundEye)
            endif
         ENDIF
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         ! bring in next set of asymmetric vortex met data to separate arrays
         select case(abs(nws))
         case(29)
            CALL NWS19GET(vortexWVNX2,vortexWVNY2,vortexPRN2,timeloc)
         case(30)
            CALL NWS20GET(vortexWVNX2,vortexWVNY2,vortexPRN2,timeloc)
         case default
            ! should be unreachable
            call logMessage(ERROR,'NWS must equal 29 or 30.')
         end select
         DO I=1,NP
            ! compute wind stress due to background met at this node
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
            ! jgf49.1001 NAM winds already contain wind reduction
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
            WSX2(I) = RampMete*airwaterdensityrat*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete*airwaterdensityrat*WDRAGCO*WINDY*WINDMAG
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                          -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WINDX
            WVNYOUT(I)=RampMete*WINDY
            ! compute wind stress due to vortex met at this node
            WindX = vortexWVNX2(I)
            WindY = vortexWVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
            CALL getBlendFactor(I, bf)
            ! blend the wind stresses and barometric pressures
            WSX2(I) = bf*RampMete*airwaterdensityrat*WDragCo*WindX*WindMag
     &               +(1.0d0-bf)*WSX2(I)
            WSY2(I) = bf*RampMete*airwaterdensityrat*WDragCo*WindY*WindMag
     &               +(1.0d0-bf)*WSY2(I)
            vortexPRDIFF = RampMete*(vortexPRN2(I)-PRBCKGRND_MH2O)
            vortexPRN2(I) = PRBCKGRND_MH2O + vortexPRDIFF
            PR2(I)= bf*vortexPRN2(I)+(1.0d0-bf)*PR2(I)
            WVNXOUT(I)=bf*RampMete*WindX+(1.0d0-bf)*WVNXOUT(I)
            WVNYOUT(I)=bf*RampMete*WindY+(1.0d0-bf)*WVNYOUT(I)
         END DO
      ENDIF

      !
      ! jgf50.38.05: Added NWS=15 for reading HWind data
       IF(ABS(NWS).EQ.15) THEN
          CALL NWS15GET(WVNX2,WVNY2,PRN2,timeloc)
          DO I=1,NP
             windx = wvnx2(i)
             windy = wvny2(i)
             windMag = sqrt(windx**2 + windy**2)
             WDragCo = WindDrag(windMag, i)
             IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(i, WDragCo,
     &              windMag, DP(I), ETA2(I),
     &              windx, windy)
                WindMag = SQRT(windx**2+windy**2)
                WDragCo = WindDrag(WindMag, i)
             ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(i,windx,windy)
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
             WSX2(I) = RampMete*airwaterdensityrat*WDRAGCO*WINDX*WINDMAG
             WSY2(I) = RampMete*airwaterdensityrat*WDRAGCO*WINDY*WINDMAG
             PR2(I) = PRBCKGRND_MH2O + RampMete*(PRN2(I)-PRBCKGRND_MH2O)
             WVNXOUT(I)=RampMete*WINDX
             WVNYOUT(I)=RampMete*WINDY
          END DO
       ENDIF
C.... tcm v51.06.02 added for GFDL Met Data
        IF(ABS(NWS).EQ.16) THEN
          IF(TimeLoc.GT.WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
              WVNX1(I)=WVNX2(I)
              WVNY1(I)=WVNY2(I)
              PRN1(I)=PRN2(I)
            END DO
            CALL NWS16GET(timeloc,WVNX2,WVNY2,PRN2)
          ENDIF
          WTRATIO=(TimeLoc-WTIME1)/WTIMINC
          DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            !Casey 110518: Enable Mark Powell's sector-based wind drag.
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
              CALL ApplyDirectionalWindReduction(I, WDragCo,
     &             WindMag, DP(I), ETA2(I), WindX, WindY)
              WindMag = SQRT(WindX*WindX+WindY*WindY)
              !Casey 110518: Enable Mark Powell's sector-based wind drag.
              WDragCo = WindDrag(WindMag,I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete*airwaterdensityrat*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete*airwaterdensityrat*WDRAGCO*WINDY*WINDMAG
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                          -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete*WINDX
            WVNYOUT(I)=RampMete*WINDY
          END DO
        ENDIF
C
C    kmd48.33bc - added in information for the sponge layer
       IF (RES_BC_FLAG.GT.0) THEN
         DO I=1,NP
            WSX2(I)=sponge(I)*WSX2(I)
            WSY2(I)=sponge(I)*WSY2(I)
         END DO
       END IF
         
!PV BEG:: NWS=17
!... v50.xx sm, v55.xx PV, -- NWS=17 for NUOPC coupling
      IF(ABS(NWS).EQ.17) THEN
         IF(TimeLoc.GT.WTIME2) THEN
            WTIME1=WTIME2
            WTIME2=WTIME2+WTIMINC
            DO I=1,NP
              WVNX1(I)=WVNX2(I)
              WVNY1(I)=WVNY2(I)
              PRN1(I)=PRN2(I)
            END DO
         ENDIF
C++ GML added 20210727 ncice=17
         if ((NCICE.ne. 0).and.
     &       TimeLoc.GT.CICE_TIME2) THEN
             CICE_TIME1 = CICE_TIME2
             CICE_TIME2 = CICE_TIME2 + CICE_TIMINC
             CICE1 = CICE2
         endif
C++
         WTRatio=(TimeLoc-WTIME1)/WTIMINC
C++ GML added 20210727
         IF(NCICE.NE.0)THEN
             CTRATIO = (TimeLoc-CICE_TIME1)/CICE_TIMINC
         ENDIF
C++
         if (WTRatio .GT. 1.0) then
             WTRatio = 1.0
             print *, '[Warning:] WTRatio is .GE. 1.0 and Set to 1.0'
         endif

         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
           IF(NCICE.NE.0) THEN
               Cicenow = 100d0*(CICE1(I) + CTRATIO*(CICE2(I)-CICE1(I)))   ! GML added 20210727
               WDragCo = WindIceDrag(WDragCo,Cicenow)                     ! GML added 20210727
               ciceout(i) = Cicenow                                       ! GML added 20210727
C               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete*airwaterdensityrat*WDragCo*WindX*WindMag   !P.V units: m^2/s^2 aka
            WSY2(I) = RampMete*airwaterdensityrat*WDragCo*WindY*WindMag   !    per water unit mass (1/RhoWat0)
C  tcm v49.16 20100617
            PRDIFF = RampMete*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            !PV Do we need to output the real wind velocities or the ramped ones?
            WVNXOUT(I)=RampMete*WindX
            WVNYOUT(I)=RampMete*WindY
         END DO
      ENDIF
!PV END:: NWS=17

C      DW/WJP: rotate wind stresses if coordinate transform
       IF ( IFSPROTS .eq. 1 ) THEN
          UVECTMP(1:NP) = WSX2(1:NP)
          VVECTMP(1:NP) = WSY2(1:NP)
          CALL DRVMAP2DSPVEC( WSX2, WSY2,
     &              UVECTMP, VVECTMP, NP, FWD = .TRUE. )
       ENDIF

C--------------------END MET FORCING---------------------------------------

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

C----------------------------------------------------------------------
      end subroutine getMeteorologicalForcing
C----------------------------------------------------------------------


C----------------------------------------------------------------------
C                   S U B R O U T I N E
C    C O L D  S T A R T   M E T E O R O L O G I C A L   F O R C I N G
C----------------------------------------------------------------------
C     Sets up meteorological i/o for a coldstart.
C----------------------------------------------------------------------
      subroutine coldstartMeteorologicalForcing(nws)
      use sizes, only : gblinputdir, inputdir
      use constants, only: g, mb2pa
      use global, only : ciceout, h0, eta2, ncice, rampmete,
     &   sponge, res_bc_flag, statim, wsx1, wsy1, wsx2, wsy2, pr2, pr1,
     &   cice1, cice2
      use owiwind,only : nws12init, nws12get, findStormCenter    ! sb added 09/xx/2006          
      use nodalAttributes, only : applyDirectionalWindReduction,
     &   loadCanopyCoef, loadDirEffRLen, applyCanopyCoefficient
      use mesh, only : dp, ics, np
#ifdef ADCNETCDF
Casey 180318: Added for NWS=13
      USE OWIWIND_NETCDF, ONLY: NWS13INIT,NWS13GET
#endif
#ifdef DATETIME
      use datetime_module, only: datetime
#endif
      implicit none
      integer, intent(in) :: nws
      real(8) :: windx ! wind vel (m/s) in x-dir
      real(8) :: windy ! wind vel (m/s) in y-dir
      real(8) :: wtratio ! used to time interpolated between met datasets
      integer :: nhg ! node number temp
      real(8) :: bf ! blend factor for multiple sets of wind data
      integer :: i ! node loop counter

      call setMessageSource("coldStartMeteorologicalForcing")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif


      WSX2(:) = 0.d0
      WSY2(:) = 0.D0
      PR2(:) = PRBCKGRND_MH2O

      WSX1(:) = 0.d0
      WSY1(:) = 0.d0
!      PR1(:) = PRBCKGRND_MH2O ! -- DW for 517
!++ GML
      if (NCICE.ne. 0) then
      cice1(:) = 0.d0
      cice2(:) = 0.D0
      endif
!++

      IF(NWS.EQ.1) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
      ENDIF

      IF(ABS(NWS).EQ.2) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
      ENDIF

      IF(NWS.EQ.3) THEN
C....    tcm_v49.04 20091124 -- Changed from local inputdir to global inputdir
C....                           to correspond with having only a global wind file.
!         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22')

         DO WHILE(.TRUE.)
            CALL NWS3GET()
            IF(IWYR.NE.IREFYR) THEN
               IWTIMEP=IWTIME
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
               END DO
               CYCLE
            ENDIF
            IF(WTIMED.LE.WREFTIM) THEN
               IWTIMEP=IWTIME
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
               END DO
               CYCLE
            ENDIF
            EXIT
         ENDDO
         write(scratchMessage,'(a,i15)') 'FOUND WIND DATA AT TIME= ',IWTIMEP
         call allMessage(INFO,scratchMessage)
         write(scratchMessage,'(a,i15)') 'FOUND WIND DATA AT TIME= ',IWTIME
         call allMessage(INFO,scratchMessage)
         WTIME2=WTIMED-WREFTIM  !CAST INTO MODEL TIME REFERENCE
         WTIME1=WTIME2-WTIMINC
      ENDIF

      IF(ABS(NWS).EQ.4) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         WTIME2=WTIME1+WTIMINC
         CALL NWS4GET()
         wvnx1 = wvnx2
         wvny1 = wvny2
         CALL NWS4GET()
      ENDIF

      IF(ABS(NWS).EQ.5) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
      ENDIF

      IF(NWS.EQ.6) THEN
C....    tcm_v49.04 20091124 -- Changed from local inputdir to global inputdir
C....                           to correspond with having only a global wind file.
!         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22')
         CALL NWS6GET()
         wvnx1 = wvnx2
         wvny1 = wvny2
         CALL NWS6GET()
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
      ENDIF

C     jgf46.00 Added option to directly apply surface stress without any
C     other correction factors.
C     xyc52.30 corrected NWS=7 (wind stress input provided on a
C     rectangular grid, same as NWS=6) 
      IF(ABS(NWS).EQ.7) THEN
C         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
C         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
C         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22')
         CALL NWS7GET()
         wvnx1 = wvnx2
         wvny1 = wvny2
         CALL NWS7GET()
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
      ENDIF

C     jgf46.02 New option to read in hurricane locations and generate
C     generate hurricane winds from the Holland Wind Model.
      IF(ABS(NWS).EQ.8) THEN
         CALL HollandGet(WVNX2,WVNY2,PRN2,statim*86400.d0)
      ENDIF

C     RJW Merged:
!    rjw added nws = 19: asymmetric hurricane winds v2.0
      IF(NWS.EQ.19) THEN
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22',STATUS='OLD')
        CALL NWS19GET(WVNX2,WVNY2,PRN2,STATIM*86400.D0)
      ENDIF

!    jie added nws = 20: generalized asymmetric vortex model
      IF(NWS.EQ.20) THEN
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22',STATUS='OLD')
        CALL NWS20GET(WVNX2,WVNY2,PRN2,STATIM*86400.D0)
      ENDIF

C.....sb46.28sb01 NWS=-12 and 12 was added to deal with raw OWI files.  09/xx/2006
      ! jgf49.1001 Added NWS=29, asymmetric vortex wind (NWS19) embedded in
      ! OWI (NWS12) basin scale wind field.
      IF ((ABS(NWS).EQ.12).OR.(ABS(NWS).EQ.29).OR.(ABS(NWS).EQ.30)) THEN
        CALL NWS12INIT(WVNX1,WVNY1,PRN1,NP,PRBCKGRND)
        CALL NWS12GET(WVNX1,WVNY1,PRN1,NP,PRBCKGRND)
        if(draglawtype.EQ.draglaw_powell)then
            call findStormCenter(eyeLat, eyeLon, foundEye)
        endif
        CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,PRBCKGRND)
        if(draglawtype.EQ.draglaw_powell)then
            call findStormCenter(eyeLat, eyeLon, foundEye)
        endif
        WTIME1 = STATIM*86400.D0
        WTIME2 = WTIME1 + WTIMINC
      ENDIF
      !
      ! jgf49.1001 Added NWS=29, asymmetric vortex wind (NWS19) embedded in
      ! OWI (NWS12) basin scale wind field.
      ! jgf: Added NWS=30 (GAHM+OWI)
      IF ((abs(NWS).EQ.29).or.(abs(NWS).eq.30)) THEN
         ! now blend in the vortex met field
         ALLOCATE(vortexWVNX2(NP),vortexWVNY2(NP),vortexPRN2(NP))
         if (abs(nws).eq.29) then
            OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'NWS_19_fort.22',
     &           STATUS='OLD')
            CALL NWS19GET(vortexWVNX2,vortexWVNY2,vortexPRN2,STATIM*86400.D0)
         endif
         if (abs(nws).eq.30) then
            OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'NWS_20_fort.22',
     &           STATUS='OLD')
            CALL NWS20GET(vortexWVNX2,vortexWVNY2,vortexPRN2,STATIM*86400.D0)
         endif
      ENDIF
      !
      ! jgf51.52.21: Fixed NWS10 and cleaned up interface.
      IF(NWS.EQ.10) THEN
         wtime1 = statim*86400.d0
         wtime2 = wtime1 + wtiminc
         call nws10init()
         wvnx1 = 0.d0
         wvny1 = 0.d0
         prn1  = prbckgrnd * mb2pa / ( rhowat0 * g)
         ! read from first available file (not fort.200) and spatially
         ! interpolate onto adcirc mesh
         IF (NCICE.eq.10) THEN
            cice1 = 0.d0
            call nws10get(wvnx2,wvny2,prn2,cice2) 
         ELSE 
            call nws10get(wvnx2,wvny2,prn2) 
         ENDIF
      ENDIF

      IF(NWS.EQ.11) THEN
         WTIME1=STATIM*86400.D0
         WTIME2=WTIME1+WTIMINC
         NWSEGWI=0
         IDSETFLG=0
         call allMessage(INFO,
     &     'The e29 met grid interpolating factors are being computed.')
         CALL NWS11GET()        !JUST COMPUTE INTERPOLATING FACTORS
         call allMessage(INFO,
     &      'Finished computing e29 interpolating factors.')
         NWSEGWI=1
         IDSETFLG=1
         wvnx1 = wvnx2
         wvny1 = wvny2
         CALL NWS11GET()        !NOW INTERPOLATE 1st WIND FIELD
      ENDIF
      !
      ! 
#ifdef DATETIME
#if defined(GRIB2API) || defined(ADCNETCDF)
      IF (ABS(NWS).eq.14) THEN
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
         ! Set datetime and make/check the grib header
         CALL NWS14INIT()
         ! Initialising winds (pressure is already initialized)
         ! Get the first wind snap
         if (NCICE.eq.14) then
            CICE_TIME1 = WTIME1
            CICE_TIME2 = CICE_TIME1 + CICE_TIMINC 
            CALL NWS14GET(WVNX1,WVNY1,PRN1,CICE1)
         else
            CALL NWS14GET(WVNX1,WVNY1,PRN1)
         endif
         ! Get owi winds (should only overwrite where they exist) 
         if (NWS.eq.-14) then
            CALL NWS12INIT(WVNX1,WVNY1,PRN1,NP,PRBCKGRND)
            CALL NWS12GET(WVNX1,WVNY1,PRN1,NP,PRBCKGRND,HasData)
            WTIME1_12 = WTIME1
            WTIME2_12 = WTIME1_12 + WTIMINC_12
         endif
         ! Get the second wind snap
         if (NCICE.eq.14) then
            CALL NWS14GET(WVNX2,WVNY2,PRN2,CICE2)
         else
            CALL NWS14GET(WVNX2,WVNY2,PRN2)
         endif
         ! Get owi winds (should only overwrite where they exist) 
         if (NWS.eq.-14) then
            CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,PRBCKGRND,HasData)
         endif
      ENDIF
#endif
#endif

#ifdef ADCNETCDF
Casey 180318: Added for NWS=13, ARC modified 190108
      IF(ABS(NWS).EQ.13)THEN
         WTIME1=STATIM*86400.D0
         WTIME2=WTIME1+WTIMINC
         WTIME2I=WTIME2
         WTIME1I=WTIME1
        CALL NWS13INIT(0.D0,W)
        CALL NWS13GET(WTIME1,WVNX1,WVNY1,PRN1,
     &      WTIME2,EyeLon,EyeLat,FoundEye,W)
        CALL NWS13GET(WTIME2,WVNX2,WVNY2,PRN2,
     &      WTIME2,EyeLon,EyeLat,FoundEye,W)
      ENDIF
#endif

      !
      ! jgf50.38.05: Added NWS=15,-15 for reading HWind data
      IF (ABS(NWS).EQ.15) THEN
         CALL NWS15INIT(STATIM*86400.0)
      ENDIF
C
C....tcm v51.06.02: Added NWS=16,-16 for reading ASCII GFDL Met data
      IF (ABS(NWS).EQ.16) THEN
         CALL INIT_GFDL(STATIM*86400.0)
      ENDIF
      !
      ! jgf52.14: Implementing Casey's fix for uninitialized
      ! wvnx1 and wvny1 arrays; they are not used/needed/initialized
      ! by the parametric vortex meteorological models but are used in
      ! the padcswan_init() subroutine to initialize the SWAN wind
      ! velocities. Failure to initialize them causes NaNs in the
      ! solution.
      if ( (abs(nws).eq.8).or.(abs(nws).eq.19).or.(abs(nws).eq.20)) then
         WVNX1 = WVNX2
         WVNY1 = WVNY2
      endif

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C----------------------------------------------------------------------
      end subroutine coldstartMeteorologicalForcing
C----------------------------------------------------------------------



C----------------------------------------------------------------------
C                   S U B R O U T I N E
C    H O T S T A R T   M E T E O R O L O G I C A L   F O R C I N G
C----------------------------------------------------------------------
C     Sets up meteorological i/o for a hotstart.
C----------------------------------------------------------------------
      subroutine hotstartMeteorologicalForcing(nws, wsx2, wsy2, pr2,
     &   wvnxout, wvnyout, timeloc, rampmete2)
      use sizes, only : gblinputdir, inputdir
      use constants, only: g
      use global, only : ciceout, h0, eta2, ncice, rampmete,
     &   sponge, res_bc_flag, statim, dtdp, dtdphs, iths,
     &   IDEN, SigTS, cice1, cice2, IFSPROTS
      use owiwind,only : nws12init, nws12get, findStormCenter    ! sb added 09/xx/2006
      use nodalAttributes, only : applyDirectionalWindReduction,
     &   loadCanopyCoef, loadDirEffRLen, applyCanopyCoefficient
      use mesh, only : dp, ics, np, uvectmp, vvectmp, DRVMAP2DSPVEC
#ifdef ADCNETCDF
      USE OWIWIND_NETCDF, ONLY: NWS13INIT,NWS13GET
#endif
      implicit none
      integer, intent(in) :: nws
      real(8), intent(out) :: wsx2(:)
      real(8), intent(out) :: wsy2(:)
      real(8), intent(out) :: pr2(:)
      real(8), intent(out) :: wvnxout(:)
      real(8), intent(out) :: wvnyout(:)
      real(8), intent(in) :: timeloc
      real(8), intent(in) :: rampMete2   ! Ramp for wind and atmospheric pressure @ITHS
      real(8) :: timeit
      real(8) :: windx ! wind vel (m/s) in x-dir
      real(8) :: windy ! wind vel (m/s) in y-dir
      real(8) :: prdiff
      real(8) :: wDragCo
      real(8) :: windMag
      real(8) :: wtratio ! used to time interpolated between met datasets
      real(8) :: cicenow, ctratio 
      integer :: nhg ! node number temp
      real(8) :: bf ! blend factor for multiple sets of wind data
      real(8) :: vortexPRDIFF
      integer :: it
      integer :: i ! node loop counter
      integer :: j
      real(8) :: wtratio_12, wtratio_14 ! used to time interpolated between met datasets
      integer :: TimeLoc_NM ! timeloc to the nearest minute

      call setMessageSource("hotStartMeteorologicalForcing")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

C
C......RESTART WIND AND PRESSURE INFORMATION

      !c DW, for nonwind and NEM
      SELECT CASE(NWS)
      CASE (0,17)
         WVNX1(:) = 0.d0
         WVNY1(:) = 0.d0
         WVNX2(:) = 0.d0
         WVNY2(:) = 0.d0
         PRN1(:)=PRBCKGRND_MH2O
         PRN2(:)=PRBCKGRND_MH2O
         WSX2(:) = 0.d0
         WSY2(:) = 0.D0
         PR2(:) = PRBCKGRND_MH2O
      END SELECT

C     No wind, radiation stress or atmospheric pressure forcings are used.
      IF(NWS.EQ.1) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         DO J=1,ITHS
            DO I=1,NP
               READ(22,*) NHG,WSX2(I),WSY2(I),PR2(I)
            END DO
         END DO
         DO I=1,NP
            WSX2(I)=RampMete2*WSX2(I)
            WSY2(I)=RampMete2*WSY2(I)
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PR2(I)
            PRDIFF = RampMete2*(PR2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

C     Wind stress and atmospheric pressure are read in at all grid nodes
C     at a time interval that does not equal the model time
C     step. Interpolation in time is used to synchronize the wind and
C     pressure information with the model time step.
      IF(NWS.EQ.2) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)  !TCM v49.02 Added read for rec. 1
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         DO IT=1,ITHS
            TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
            IF(TIMEIT.GT.WTIME2) THEN
               WTIME1=WTIME2
               WTIME2=WTIME2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
                  PRN1(I)=PRN2(I)
                  READ(22,*) NHG,WVNX2(I),WVNY2(I),PRN2(I)
               END DO
            ENDIF
         END DO
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WSX2(I) = RampMete2*WINDX
            WSY2(I) = RampMete2*WINDY
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                          -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

      IF(NWS.EQ.-2) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = TimeLoc
         WTIME2 = WTIME1 + WTIMINC
         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         DO I=1,NP
            WSX2(I) = RampMete2*WVNX1(I)
            WSY2(I) = RampMete2*WVNY1(I)
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PRN1(I)
            PRDIFF = RampMete2*(PRN1(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

C     Wind velocity in US Navy Fleet Numeric format interpolated in
C     space onto the ADCIRC grid and in time to synchronize the wind and
C     pressure information with the model time step. Garratt's formula
C     is used to compute wind stress from the wind velocity.
      IF(NWS.EQ.3) THEN
C....    tcm_v49.04 20091124 -- Changed from local inputdir to global inputdir
C....                           to correspond with having only a global wind file.
!         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22')
         DO WHILE(.TRUE.)
            CALL NWS3GET()
            IF(IWYR.NE.IREFYR) THEN
               IWTIMEP=IWTIME
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
               END DO
               CYCLE
            ENDIF
            IF(WTIMED.LE.WREFTIM) THEN
               IWTIMEP=IWTIME
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
               END DO
               CYCLE
            ENDIF
            EXIT
         ENDDO
         write(scratchMessage,'(a,f15.7)') 'FOUND WIND DATA AT TIME= ',IWTIMEP
         call allMessage(INFO,scratchMessage)
         write(scratchMessage,'(a,f15.7)') 'FOUND WIND DATA AT TIME= ',IWTIME
         call allMessage(INFO,scratchMessage)
         WTIME2=WTIMED-WREFTIM  !CAST INTO MODEL TIME REFERENCE
         WTIME1=WTIME2-WTIMINC
         DO IT=1,ITHS
            TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
            IF(TIMEIT.GT.WTIME2) THEN
               WTIME1=WTIME2
               WTIME2=WTIME2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
               END DO
               CALL NWS3GET()
               write(scratchMessage,'(a,f15.7)')
     &            'Wind file advanced to time = ',iwtime
               call allMessage(INFO,scratchMessage)
            ENDIF
         END DO
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP              !INTERPOLATE IN TIME
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG=SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I)=RampMete2*airwaterdensityrat*WDRAGCO*WINDX*WINDMAG
            WSY2(I)=RampMete2*airwaterdensityrat*WDRAGCO*WINDY*WINDMAG
         END DO
      ENDIF

C     Wind velocity and atmospheric pressure are read in (PBL/JAG
C     format) at selected ADCIRC grid nodes. Interpolation in time is
C     used to synchronize the wind and pressure information with the
C     model time step. Garratt's formula is used to compute wind stress
C     from wind velocity.
      IF(NWS.EQ.4) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
         CALL NWS4GET()  !TCM v49.02 Added read for rec. 1
         wvnx1 = wvnx2
         wvny1 = wvny2
         CALL NWS4GET()
         DO IT=1,ITHS
            TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
            IF(TIMEIT.GT.WTIME2) THEN
               WTIME1=WTIME2
               WTIME2=WTIME2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
                  PRN1(I)=PRN2(I)
               END DO
               CALL NWS4GET()
            ENDIF
         END DO
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

      IF(NWS.EQ.-4) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = TimeLoc
         WTIME2 = WTIME1 + WTIMINC
         CALL NWS4GET()
         wvnx1 = wvnx2
         wvny1 = wvny2
         CALL NWS4GET()
         DO I=1,NP
            WINDX = WVNX1(I)
            WINDY = WVNY1(I)
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PRN1(I)
            PRDIFF = RampMete2*(PRN1(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

C     Wind velocity and atmospheric pressure are read in at all grid
C     nodes. Interpolation in time is used to synchronize the wind and
C     pressure information with the model time step. Garratt's formula
C     is used to compute wind stress from wind velocity.
      IF(NWS.EQ.5) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         WTIME2 = WTIME1 + WTIMINC
         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)  !TCM v49.02 Added read for rec. 1
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         DO IT=1,ITHS
            TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
            IF(TIMEIT.GT.WTIME2) THEN
               WTIME1=WTIME2
               WTIME2=WTIME2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
                  PRN1(I)=PRN2(I)
                  READ(22,*) NHG,WVNX2(I),WVNY2(I),PRN2(I)
               END DO
            ENDIF
         END DO
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
C     jjw-42.06j wrote, jgf46.00 added following two lines:
            WVNXOUT(I)=RampMete2*WINDX
            WVNYOUT(I)=RampMete2*WINDY
         END DO
      ENDIF

      IF(NWS.EQ.-5) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = TimeLoc
         WTIME2 = WTIME1 + WTIMINC
         READ(22,*) (NHG,WVNX1(I),WVNY1(I),PRN1(I),I=1,NP)
         READ(22,*) (NHG,WVNX2(I),WVNY2(I),PRN2(I),I=1,NP)
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC !jgf46.00 added
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.01 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PRN1(I)
            PRDIFF = RampMete2*(PRN1(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete2*WINDX !jgf46.00 added  ! tcm 20100618 v49.16 changed RampMete to RampMete2
            WVNYOUT(I)=RampMete2*WINDY !jgf46.00 added  ! tcm 20100618 v49.16 changed RampMete to RampMete2
         END DO
      ENDIF

C     Wind velocity and atmospheric pressure are read in for a
C     rectangular grid (either in Longitude, Latitude or Cartesian
C     coordinates, consistent with the grid coordinates) and
C     interpolated in space onto the ADCIRC grid and in time to
C     synchronize the wind and pressure information with the model time
C     step. Garratt's formula is used to compute wind stress from the
C     wind velocity.
      IF(NWS.EQ.6) THEN
C....    tcm_v49.04 20091124 -- Changed from local inputdir to global inputdir
C....                           to correspond with having only a global wind file.
!         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22')
C     The following 3 lines are a hardwire to allow a non standard met
C     file to be read in at time zero in a hot start.  They should be
C     eliminated or commented out for normal operation
c     OPEN(199,FILE=TRIM(INPUTDIR)//'/'//'fort.199')
c     READ(199,*) (NHG,PRN1(I),WVNX1(I),WVNY1(I),I=1,NP)
c     CLOSE(199)
C     The following CALL statement should be uncommented for normal operation
         CALL NWS6GET()
         wvnx1 = wvnx2
         wvny1 = wvny2
         CALL NWS6GET()
         WTIME1=TimeLoc
         WTIME2=WTIME1+WTIMINC
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

C     jgf46.01 New option to read in surface wind stress and atmospheric
C     pressure for a rectangular grid (either in Longitude, Latitude or
C     Cartesian coordinates, consistent with the grid coordinates) and
C     interpolate in space onto the ADCIRC grid. Interpolation in time
C     is used to synchronize the wind and pressure information with the
C     model time step.
C
      IF(NWS.EQ.7) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         CALL NWS7GET()  !TCMv49.02 Changed (WVNX2,WVNY2,PRN2) to (WVNX1,WVNY1,PRN1)
         wvnx1 = wvnx2
         wvny1 = wvny2
         CALL NWS7GET()
         WTIME1=TimeLoc
         WTIME2=WTIME1+WTIMINC
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WSX2(I) = RampMete2*WindX !apply ramp
            WSY2(I) = RampMete2*WindY
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            wvnxout(i)=WSX2(i) !for met recording sta. output
            wvnyout(i)=WSY2(i)
         END DO
      ENDIF

      IF(NWS.EQ.-7) THEN
         OPEN(22,FILE=TRIM(INPUTDIR)//'/'//'fort.22')
         WTIME1 = STATIM*86400.D0
         CALL NWS7GET()
         wvnx1 = wvnx2
         wvny1 = wvny2
         CALL NWS7GET()
         WTIME1=TimeLoc
         WTIME2=WTIME1+WTIMINC
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
C     jgfdebug46.01 How do we convert these marine wind stresses into
C     directional land surface stresses?
         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WSX2(I) = RampMete2*WindX !apply ramp
            WSY2(I) = RampMete2*WindY
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            wvnxout(i)=WSX2(i) !for met recording sta. output
            wvnyout(i)=WSY2(i)
         END DO
      ENDIF

C     jgf46.02 New option to read in hurricane locations and generate
C     generate hurricane winds from the Holland Wind Model.
      IF(ABS(NWS).EQ.8) THEN
         CALL HollandGet(WVNX2,WVNY2,PRN2,timeloc)
         DO I=1,NP
            WindX = WVNX2(I)
            WindY = WVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete2*airwaterdensityrat*WDragCo*WindX*WindMag
            WSY2(I) = RampMete2*airwaterdensityrat*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PRN2(I)
            PRDIFF = RampMete2*(PRN2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete2*WindX
            WVNYOUT(I)=RampMete2*WindY
         ENDDO
      ENDIF

!     rjw added nws = 19: asymmetric hurricane winds v2.0

      IF(ABS(NWS).EQ.19) THEN
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22',STATUS='OLD')
         CALL NWS19GET(WVNX2,WVNY2,PRN2,TimeLoc)
         DO I=1,NP
            WindX = WVNX2(I)
            WindY = WVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete2*airwaterdensityrat*WDragCo*WindX*WindMag
            WSY2(I) = RampMete2*airwaterdensityrat*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PRN2(I)
            PRDIFF = RampMete2*(PRN2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete2*WindX
            WVNYOUT(I)=RampMete2*WindY
         ENDDO
      ENDIF

!     jie added nws = 20: generalized asymmetric vortex model

      IF(ABS(NWS).EQ.20) THEN
         OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'fort.22',STATUS='OLD')
         CALL NWS20GET(WVNX2,WVNY2,PRN2,TimeLoc)
         DO I=1,NP
            WindX = WVNX2(I)
            WindY = WVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete2*airwaterdensityrat*WDragCo*WindX*WindMag
            WSY2(I) = RampMete2*airwaterdensityrat*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*PRN2(I)
            PRDIFF = RampMete2*(PRN2(I)-PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete2*WindX
            WVNYOUT(I)=RampMete2*WindY
         ENDDO
      ENDIF
C

C     Wind velocity (10 m) and atmospheric pressure are read in from a
C     sequence of National Weather Service (NWS) Aviation (AVN) model
C     output files. Each AVN file is assumed to contain data on a
C     Gaussian longitude, latitude grid at a single time. Consecutive
C     files in the sequence are separated by N hours in time. Garratt's
C     formula is used to compute wind stress from the wind velocity.
      IF(NWS.EQ.10) THEN
         WTIME1 = TimeLoc
         WTIME2 = WTIME1 + WTIMINC
         !
         ! jgf51.52.21: Fixed NWS10 and cleaned up interface.
         call nws10init()
         ! read from fort.200 and spatially interpolate onto adcirc mesh
         IF (NCICE.eq.10) THEN
            ! read from fort.200 and spatially interpolate onto adcirc mesh
            call nws10get(wvnx1,wvny1,prn1,cice1)          
            ! read from first file after fort.200 and spatially interpolate onto adcirc mesh
            call nws10get(wvnx2,wvny2,prn2,cice2)
         else
            ! read from fort.200 and spatially interpolate onto adcirc mesh
            call nws10get(wvnx1,wvny1,prn1)          
            ! read from first file after fort.200 and spatially interpolate onto adcirc mesh
            call nws10get(wvnx2,wvny2,prn2)
         endif
         !
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF (NCICE.eq.10) then
               Cicenow = 100d0*( CICE1(I) + 
     &                           WTRatio*( CICE2(I)-CICE1(I) ) )
               WDragCo = WindIceDrag(WDragCo,Cicenow)
            ELSEIF (NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDY*WINDMAG
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF

C     Wind velocity (10 m) and atmospheric pressure are read in from a
C     sequence of stripped down National Weather Service (NWS) ETA 29km
C     model output files. Each ETA file is assumed to contain data on an
C     E grid for a single day (8 data sets, one every 3 hours, beginning
C     @ 03:00 and continuing through 24:00 of the given day). The wind
C     data is converted to an east-west, north-south coordinate system
C     inside ADCIRC. Garratt's formula is used to compute wind stress
C     from the wind velocity.
      IF(NWS.EQ.11) THEN
         WTIME1=TimeLoc
         WTIME2=WTIME1+WTIMINC
         NWSEGWI=0
         IDSETFLG=0
         call allMessage(INFO,
     &      'The e29 met grid interpolating factors are being computed.')
         CALL NWS11GET()        !JUST COMPUTE INTERPOLATING FACTORS
         call allMessage(INFO,
     &      'Finished computing e29 interpolating factors.')
         NWSEGWI=0
         IDSETFLG=8
         CALL NWS11GET()        !NOW COMPUTE HOTSTART WIND FILED
         NWSEGWI=1
         IDSETFLG=1
         wvnx1 = wvnx2
         wvny1 = wvny2
         CALL NWS11GET()        !NOW COMPUTE NEXT WIND FIELD
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDY*WINDMAG
C  tcm v49.16 20100617
!            PR2(I)=RampMete2*(PRN1(I)+WTRATIO*(PRN2(I)-PRN1(I)))
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
         END DO
      ENDIF
C
C.....sb46.28sb01 NWS=12,-12 were added to deal with raw OWI files.  09/xx/2006
      ! jgf: Combined NWS=12 and NWS=-12 and added abs(NWS)=29 and abs(NWS)=30
      IF((abs(NWS).eq.12).or.(abs(NWS).eq.29).or.(abs(NWS).eq.30)) THEN
         WTIME1 = TimeLoc      ! nws=-12, wind data relative to hot start time
         if (nws.gt.0) then
            WTIME1 = STATIM*86400.D0 ! wind data relative to cold start time
         endif
         WTIME2 = WTIME1 + WTIMINC
         CALL NWS12INIT(WVNX1,WVNY1,PRN1,NP,PRBCKGRND)  ! This just initializes some variables and might set WVNX1,etc... if skipping ahead.
         CALL NWS12GET(WVNX1,WVNY1,PRN1,NP,PRBCKGRND)  !TCM v49.02 Added read for rec. 1
         CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,PRBCKGRND)
          ! fast forward if wind data are relative to cold start time
         if (nws.gt.0) then
            DO IT=1,ITHS
               TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
               IF(TIMEIT.GT.WTIME2) THEN
                  WTIME1=WTIME2
                  WTIME2=WTIME2+WTIMINC
                  DO I=1,NP
                     WVNX1(I)=WVNX2(I)
                     WVNY1(I)=WVNY2(I)
                     PRN1(I)=PRN2(I)
                  END DO
                  CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,PRBCKGRND)
               ENDIF
            END DO
         endif
         WTRATIO=(TimeLoc-WTIME1)/WTIMINC
         DO I=1,NP
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            WSX2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDY*WINDMAG
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                          -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete2*WINDX
            WVNYOUT(I)=RampMete2*WINDY
         END DO
      ENDIF

#ifdef ADCNETCDF
Casey 180318: Added for NWS=13, ARC modified 190108
      IF(ABS(NWS).EQ.13)THEN
         WTIME1 = TimeLoc
         WTIME2 = WTIME1 + WTIMINC
         WTIME2I = WTIME2
         WTIME1I = WTIME1
        CALL NWS13INIT(TIMELOC,W)
        CALL NWS13GET(WTIME1,WVNX1,WVNY1,PRN1,
     &      WTIME2,EyeLon,EyeLat,FoundEye,W)
        CALL NWS13GET(WTIME2,WVNX2,WVNY2,PRN2,
     &      WTIME2,EyeLon,EyeLat,FoundEye,W)
        WTRATIO=(TimeLoc-WTIME1)/WTIMINC
        DO I=1,NP
          WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
          WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
          WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
          WDragCo = WindDrag(WindMag, I)
          IF(LoadDirEffRLen)THEN
             CALL ApplyDirectionalWindReduction(I, WDragCo,
     &         WindMag, DP(I), ETA2(I), WindX, WindY)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
          ENDIF
          IF(LoadCanopyCoef)THEN
            CALL ApplyCanopyCoefficient(I,WindX,WindY)
          ENDIF
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
          WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
          WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
          PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &        -PRBCKGRND_MH2O)
          PR2(I) = PRBCKGRND_MH2O + PRDIFF
          WVNXOUT(I)=RampMete2*WINDX
          WVNYOUT(I)=RampMete2*WINDY
         ENDDO
      ENDIF
#endif

C
C     jgf49.1001 Added NWS29, asymmetric vortex winds (NWS19)
C     embedded in an OWI basin scale gridded met field (NWS12)
C     from NAM data
      ! jgf: Added NWS30 (GAHM+OWI)
      IF ((abs(NWS).EQ.29).or.(abs(NWS).eq.30)) THEN
         ! Now blend in the vortex winds into the background OWI wind field.
         ALLOCATE(vortexWVNX2(NP),vortexWVNY2(NP),vortexPRN2(NP))
         if (abs(nws).eq.29) then
            OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'NWS_19_fort.22',
     &           STATUS='OLD')
            CALL NWS19GET(vortexWVNX2,vortexWVNY2,vortexPRN2,TimeLoc)
         endif
         if (abs(nws).eq.30) then
            OPEN(22,FILE=TRIM(GBLINPUTDIR)//'/'//'NWS_20_fort.22',
     &           STATUS='OLD')
            CALL NWS20GET(vortexWVNX2,vortexWVNY2,vortexPRN2,TimeLoc)
         endif
         DO I=1,NP
            WindX = vortexWVNX2(I)
            WindY = vortexWVNY2(I)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
            CALL getBlendFactor(I, bf)
            WSX2(I) = bf*RampMete2*airwaterdensityrat*WDragCo*WindX*WindMag
     &               +(1.0d0-bf)*WSX2(I)
            WSY2(I) = bf*RampMete2*airwaterdensityrat*WDragCo*WindY*WindMag
     &               +(1.0d0-bf)*WSY2(I)
            vortexPRDIFF = RampMete2*(vortexPRN2(I)-PRBCKGRND_MH2O)
            vortexPRN2(I) = PRBCKGRND_MH2O + vortexPRDIFF
            ! blend vortex into background
            PR2(I) = bf*vortexPRN2(I)+(1.0d0-bf)*PR2(I)
            WVNXOUT(I)=bf*RampMete2*WindX+(1.0d0-bf)*WVNXOUT(I)
            WVNYOUT(I)=bf*RampMete2*WindY+(1.0d0-bf)*WVNYOUT(I)
         ENDDO
      ENDIF
      !
#ifdef DATETIME
#if defined(GRIB2API) || defined(ADCNETCDF)
!     WJP Feb27-2019: Added for NWS=14 for grib2
      IF (ABS(NWS).eq.14) THEN
         ! TimeLoc should be rounded down to nearest WTIMINC
         TimeLoc_NM = nint(TimeLoc/60d0) ! first round to nearest minute
         WTIME1 = WTIMINC*floor(TimeLoc_NM*60d0/WTIMINC)
         WTIME2 = WTIME1 + WTIMINC
         WTRATIO = (TimeLoc-WTIME1)/WTIMINC
         IF(NCICE.NE.0)THEN
             CTRATIO = (TimeLoc-CICE_TIME1)/CICE_TIMINC
         ENDIF
         ! Set datetime and make/check the grib header
         CALL NWS14INIT()
         ! Initialising winds (pressure is already initialized)
         ! Get the first wind snap
         if (NCICE.eq.14) then
            CICE_TIME1 = CICE_TIMINC*floor(TimeLoc_NM*60d0/CICE_TIMINC)
            CICE_TIME2 = CICE_TIME1 + CICE_TIMINC
            CALL NWS14GET(WVNX1,WVNY1,PRN1,CICE1)
         else
            CALL NWS14GET(WVNX1,WVNY1,PRN1)
         endif
         ! Get owi winds (should only overwrite where they exist) 
         if (NWS.eq.-14) then
            CALL NWS12INIT(WVNX1,WVNY1,PRN1,NP,PRBCKGRND)
            CALL NWS12GET(WVNX1,WVNY1,PRN1,NP,PRBCKGRND,HasData)
            WTIME1_12 = WTIME1
            WTIME2_12 = WTIME1_12 + WTIMINC_12
            WTRATIO_12 = (TimeLoc-WTIME1_12)/WTIMINC_12
            WTRATIO_14 = WTRATIO
         endif
         ! Get the second wind snap
         if (NCICE.eq.14) then
            CALL NWS14GET(WVNX2,WVNY2,PRN2,CICE2)
         else
            CALL NWS14GET(WVNX2,WVNY2,PRN2)
         endif
         ! Get owi winds (should only overwrite where they exist) 
         if (NWS.eq.-14) then
            CALL NWS12GET(WVNX2,WVNY2,PRN2,NP,PRBCKGRND,HasData)
         endif
         DO I=1,NP
            IF (NWS.eq.-14) then
               WTRATIO = WTRATIO_14
               IF (HasData(I)) WTRATIO = WTRATIO_12
            ENDIF
            WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
            WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
            WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
            !Casey 110518: Enable Mark Powell's sector-based wind drag.
            WDragCo = WindDrag(WindMag, I)
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               !Casey 110518: Enable Mark Powell's sector-based wind drag.
               WDragCo = WindDrag(WindMag,I)
            ENDIF
            IF (LoadCanopyCoef) then
               CALL ApplyCanopyCoefficient(I,WindX,WindY)
            ENDIF
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               if (NCICE.eq.14) then
                  Cicenow = 100d0*(CICE1(I) + 
     &                             CTRATIO*(CICE2(I)-CICE1(I)))
                  WDragCo = WindIceDrag(WDragCo,Cicenow)
               else
                  WDragCo = WindIceDrag(WDragCo,ciceout(i))
               endif
            ENDIF
            ! WJP: added varying sea surface density
            if (abs(IDEN).ge.5) then
               if (SigTS(I) < -1d2) then
                  airwaterdensityrat = rhoAir/RhoWat0
               else
                  airwaterdensityrat = rhoAir/(SigTS(I) + RhoWat0)
               endif
            endif
            WSX2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDX*WINDMAG
            WSY2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDY*WINDMAG
            PRDIFF  = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                          -PRBCKGRND_MH2O)
#ifdef NOIVB
            PR2(I)  = PRBCKGRND_MH2O
#else
            PR2(I)  = PRBCKGRND_MH2O + PRDIFF
#endif
            WVNXOUT(I) = RampMete2*WINDX
            WVNYOUT(I) = RampMete2*WINDY

          END DO
      ENDIF
#endif
#endif
#ifdef ADCNETCDF
Casey 180318: Added for NWS=13
      IF(ABS(NWS).EQ.13)THEN
        WTIME1 = TIMELOC
        WTIME2 = TIMELOC
        CALL NWS13INIT(TIMELOC,W)
        CALL NWS13GET(TIMELOC,WVNX2,WVNY2,PRN2,
     &         WTIME2,EyeLon,EyeLat,FoundEye,W)
!       IF(TRIM(DragLawString).EQ."Powell".OR.
!    &      TRIM(DragLawString).EQ."POWELL".OR.
!    &      TRIM(DragLawString).EQ."powell")THEN
!         CALL findStormCenter(eyeLatR,eyeLonR,foundEye)
!       ENDIF
        WTRATIO=(TimeLoc-WTIME1)/WTIMINC
        DO I=1,NP
          WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
          WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
          WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
          WDragCo = WindDrag(WindMag, I)
          IF(LoadDirEffRLen)THEN
            CALL ApplyDirectionalWindReduction(I, WDragCo,
     &          WindMag, DP(I), ETA2(I), WindX, WindY)
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)
          ENDIF
          IF(LoadCanopyCoef)THEN
            CALL ApplyCanopyCoefficient(I,WindX,WindY)
          ENDIF
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
          WSX2(I) = RampMete2*0.001293d0*WDRAGCO*WINDX*WINDMAG
          WSY2(I) = RampMete2*0.001293d0*WDRAGCO*WINDY*WINDMAG
          PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &        -PRBCKGRND_MH2O)
          PR2(I) = PRBCKGRND_MH2O + PRDIFF
          WVNXOUT(I)=RampMete2*WINDX
          WVNYOUT(I)=RampMete2*WINDY
         ENDDO
      ENDIF
#endif
      !
      ! jgf50.38.05: Added NWS=15 for reading HWind data
       IF(ABS(NWS).EQ.15) THEN
          CALL NWS15INIT(timeloc)
          CALL NWS15GET(WVNX2,WVNY2,PRN2,timeloc)
          DO I=1,NP
             WDragCo = WindDrag(windSpeeds(i), i)
             windx = wvnx2(i)
             windy = wvny2(i)
             IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(i, WDragCo,
     &              windSpeeds(i), DP(I), ETA2(I),
     &              windx, windy)
                WindMag = SQRT(windx**2+windy**2)
                WDragCo = WindDrag(WindMag, i)
             ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(i,windx,windy)
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
             WSX2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDX*WINDMAG
             WSY2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDY*WINDMAG
             PR2(I) = PRBCKGRND_MH2O + RampMete2*(PRN2(I)-PRBCKGRND_MH2O)
             WVNXOUT(I)=RampMete2*WINDX
             WVNYOUT(I)=RampMete2*WINDY
          END DO
       ENDIF
C
C.....tcm v51.06.02 NWS=16,-16 were added to deal with GFDL Met files.
      IF(NWS.EQ.16) THEN
          WTIME1 = STATIM*86400.D0
          WTIME2 = WTIME1 + WTIMINC
          CALL INIT_GFDL(wtime1)  !this will read in the two time series
          DO IT=1,ITHS
             TIMEIT=IT*DTDPHS + STATIM*86400.D0 ! kmd48.33bc - changed the DTDP to DTDPHS
             IF(TIMEIT.GT.WTIME2) THEN
                WTIME1=WTIME2
                WTIME2=WTIME2+WTIMINC
                DO I=1,NP
                   WVNX1(I)=WVNX2(I)
                   WVNY1(I)=WVNY2(I)
                   PRN1(I)=PRN2(I)
                 END DO
                 CALL NWS16GET(wtime2,WVNX2,WVNY2,PRN2)
              ENDIF
          END DO
          WTRATIO=(TimeLoc-WTIME1)/WTIMINC
          DO I=1,NP
             WINDX = WVNX1(I) + WTRATIO*(WVNX2(I)-WVNX1(I))
             WINDY = WVNY1(I) + WTRATIO*(WVNY2(I)-WVNY1(I))
             WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
             WDragCo = WindDrag(WindMag, I)
             IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY)
                WDragCo = WindDrag(WindMag, I)
             ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
             WSX2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDX*WINDMAG
             WSY2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDY*WINDMAG
             PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
             PR2(I) = PRBCKGRND_MH2O + PRDIFF
             WVNXOUT(I)=RampMete2*WINDX
             WVNYOUT(I)=RampMete2*WINDY
          END DO
       ENDIF
C
       IF(NWS.EQ.-16) THEN
          WTIME1 = TimeLoc
          WTIME2 = WTIME1 + WTIMINC
          call INIT_GFDL(wtime1)  !this reads in the first two time series values
          DO I=1,NP
             WINDX = WVNX1(I)
             WINDY = WVNY1(I)
             WINDMAG = SQRT(WINDX*WINDX+WINDY*WINDY)
             WDragCo = WindDrag(WindMag, I)
             IF (LoadDirEffRLen) THEN
                CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
                WindMag = SQRT(WindX*WindX+WindY*WindY)
                WDragCo = WindDrag(WindMag, I)
             ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
            IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
             WSX2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDX*WINDMAG
             WSY2(I) = RampMete2*airwaterdensityrat*WDRAGCO*WINDY*WINDMAG
             PRDIFF = RampMete2*(PRN1(I)-PRBCKGRND_MH2O)
             PR2(I) = PRBCKGRND_MH2O + PRDIFF
             WVNXOUT(I)=RampMete2*WINDX  ! tcm 20100618 v49.16 changed RampMete to RampMete2
             WVNYOUT(I)=RampMete2*WINDY  ! tcm 20100618 v49.16 changed RampMete to RampMete2
          END DO
       ENDIF

!PV BEG:: NWS=17, there is no NWS=-17 option
!... v55.xx PV -- NWS=17 for NUOPC coupling
      IF(ABS(NWS).EQ.17) THEN ! using ABS(NWS) to make sure that this block is always executed for NWS=17/-17
         WTIME1 = TimeLoc
         WTIME2 = WTIME1 + WTIMINC
         DO IT=1,ITHS
            TIMEIT=IT*DTDPHS + STATIM*86400.D0
            IF(TIMEIT.GT.WTIME2) THEN
               WTIME1=WTIME2
               WTIME2=WTIME2+WTIMINC
               DO I=1,NP
                  WVNX1(I)=WVNX2(I)
                  WVNY1(I)=WVNY2(I)
                  PRN1(I)=PRN2(I)
                END DO
             ENDIF
         END DO

         WTRatio=(TimeLoc-WTIME1)/WTIMINC

         if (WTRatio .GT. 1.0) then
             WTRatio = 1.0
             print *, '[Warning:] WTRatio is .GE. 1.0 and Set to 1.0'
         endif

         DO I=1,NP
            WindX = WVNX1(I) + WTRatio*(WVNX2(I)-WVNX1(I))
            WindY = WVNY1(I) + WTRatio*(WVNY2(I)-WVNY1(I))
            WindMag = SQRT(WindX*WindX+WindY*WindY)
            WDragCo = WindDrag(WindMag, I)

C     jgf46.00 Add directional wind reduction.
            IF (LoadDirEffRLen) THEN
               CALL ApplyDirectionalWindReduction(I, WDragCo,
     &              WindMag, DP(I), ETA2(I), WindX, WindY)
               WindMag = SQRT(WindX*WindX+WindY*WindY)
               WDragCo = WindDrag(WindMag, I)
            ENDIF
            IF (LoadCanopyCoef)
     &         CALL ApplyCanopyCoefficient(I,WindX,WindY)
C    TCM V49.64.01 ADDED ICE EFFECTS ON WIND DRAG COEFF
           IF(NCICE.NE.0) THEN
               WDragCo = WindIceDrag(WDragCo,ciceout(i))
            ENDIF
            WSX2(I) = RampMete2*airwaterdensityrat*WDragCo*WindX*WindMag
            WSY2(I) = RampMete2*airwaterdensityrat*WDragCo*WindY*WindMag
C  tcm v49.16 20100617
            PRDIFF = RampMete2*((PRN1(I)+WTRatio*(PRN2(I)-PRN1(I)))
     &                           -PRBCKGRND_MH2O)
            PR2(I) = PRBCKGRND_MH2O + PRDIFF
            WVNXOUT(I)=RampMete2*WindX
            WVNYOUT(I)=RampMete2*WindY
         END DO
      ENDIF
!PV END:: NWS=17

       !
       ! jgf52.14: Implementing Casey's fix for uninitialized
       ! wvnx1 and wvny1 arrays; they are not used/needed/initialized
       ! by the parametric vortex meteorological models but are used in
       ! the padcswan_init() subroutine to initialize the SWAN wind 
       ! velocities. Failure to initialize them causes NaNs in the 
       ! solution. 
       if ( (abs(nws).eq.8).or.(abs(nws).eq.19).or.(abs(nws).eq.20)) then
          WVNX1 = WVNX2
          WVNY1 = WVNY2
       endif

C      DW/WJP: rotate wind stresses if coordinate transform
       IF ( IFSPROTS .eq. 1 ) THEN
          UVECTMP(1:NP) = WSX2(1:NP); 
          VVECTMP(1:NP) = WSY2(1:NP);
          CALL DRVMAP2DSPVEC( WSX2, WSY2,
     &              UVECTMP, VVECTMP, NP, FWD = .TRUE. ) ;
       ENDIF
      
#if defined(WIND_TRACE) || defined(ALL_TRACE)
       call allMessage(DEBUG,"Return.")
#endif
       call unsetMessageSource()
      
C--------------------END MET FORCING-----------------------------------
      end subroutine hotstartMeteorologicalForcing
C--------------------END MET FORCING-----------------------------------



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

      SUBROUTINE NWS3GET()
      USE SIZES
      USE MESH, ONLY : X,Y,SLAM,SFEA,NP,ICS
      IMPLICIT NONE
      INTEGER, SAVE :: FIRSTCALL = 0
      INTEGER I,J
      INTEGER NScreen,ScreenUnit
      REAL*8 XCOOR,YCOOR
      INTEGER  LATIND1,LATIND2,LONIND1,LONIND2
      REAL(8) WSPEED,WDIR
      REAL(8) WLATM,WLONM,XWRATIO,YWRATIO
      REAL(8),ALLOCATABLE,SAVE :: WVXFN(:,:),WVYFN(:,:)
      INTEGER :: IWYR, IWMO, IWDAY, IWHR

      call setMessageSource("nws3get")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif


C
      IF (FIRSTCALL.EQ.0) THEN
         FIRSTCALL = 1
         ALLOCATE ( WVXFN(NWLAT,NWLON),WVYFN(NWLAT,NWLON) )
      ENDIF
C
      READ(22,*) IWTIME
      IWYR = IWTIME/1000000
      IWMO = IWTIME/10000 - IWYR*100
      IWDAY = IWTIME/100 - IWYR*10000 - IWMO*100
      IWHR = IWTIME - IWYR*1000000 - IWMO*10000 - IWDAY*100
      CALL TIMECONV(IWYR,IWMO,IWDAY,IWHR,0,0.0D0,WTIMED)

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
         IF(ICS.NE.1) THEN
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
         WVNX2(I) = WVXFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVXFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVXFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVXFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         WVNY2(I) = WVYFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVYFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVYFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVYFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
      END DO
      PRN2 = PRBCKGRND_MH2O

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
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

      SUBROUTINE NWS4GET()
      USE SIZES
      USE CONSTANTS, ONLY : RHOWAT0, G
      USE MESH, ONLY : NP
      IMPLICIT NONE
      INTEGER   I,NHG
      REAL(8)  RHOWATG
      CHARACTER*80 PBLJAGF

      call setMessageSource("nws4get")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

C
      RHOWATG=RHOWAT0*G
      DO I=1,NP
        WVNX2(I)=0.d0
        WVNY2(I)=0.d0
        PRN2(I)=101300.d0/RHOWATG
      END DO
      DO WHILE(.TRUE.)
         READ(22,'(A80)') PBLJAGF
         IF(PBLJAGF(2:2).NE.'#') EXIT
      ENDDO
      DO WHILE(.TRUE.)
 171     READ(PBLJAGF,'(I8,5E13.5)') NHG,WVNX2(NHG),WVNY2(NHG),
     &                               PRN2(NHG)
C
C        jgf46.02 From now on, wind files must contain data that are
C        appropriate for the time increment listed within the
C        file. Therefore, the following two lines were commented out.
C        WVNX(NHG)=WVNX(NHG)*1.04d0*0.5144d0 !CONVERT 30-MIN WINDS IN
C        WVNY(NHG)=WVNY(NHG)*1.04d0*0.5144d0 !KNOTS TO 10-MIN WIND IN M/S
C        jgf46.02 Added the following two lines.
         WVNX2(NHG)=WVNX2(NHG)*0.5144d0 !CONVERT KNOTS TO  M/S
         WVNY2(NHG)=WVNY2(NHG)*0.5144d0
         PRN2(NHG)=100.d0*PRN2(NHG)/RHOWATG !CONVERT MILLIBARS TO M OF WATER
         IF(PRN2(NHG).EQ.0.) PRN2(NHG)=101300.d0/RHOWATG                       !RAL0315+ OK
         READ(22,'(A80)') PBLJAGF
         IF(PBLJAGF(2:2).EQ.'#') EXIT
      ENDDO
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
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
      SUBROUTINE NWS6GET()
      USE SIZES
      USE CONSTANTS, ONLY : RHOWAT0, G
      USE MESH, ONLY : X,Y,SLAM,SFEA,NP,ICS
      IMPLICIT NONE
      INTEGER, SAVE :: FIRSTCALL = 0
      INTEGER :: I,J
      REAL(8) :: RHOWATG
      INTEGER  LATIND1,LATIND2,LONIND1,LONIND2
      REAL(8) XWRATIO,YWRATIO
      REAL(8) WLATM,WLONM
      REAL(8) :: XCOOR,YCOOR
      REAL(8),SAVE,ALLOCATABLE :: WVXFN(:,:),WVYFN(:,:),PRN(:,:)

      call setMessageSource("nws6get")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

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
         IF(ICS.NE.1) THEN
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
         WVNX2(I) = WVXFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVXFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVXFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVXFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         WVNY2(I) = WVYFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVYFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVYFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVYFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         PRN2(I) = PRN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + PRN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + PRN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + PRN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         PRN2(I) = PRN2(I)/RHOWATG

      END DO
C
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C     ----------------------------------------------------------------
      END SUBROUTINE NWS6GET
C     ----------------------------------------------------------------


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
      SUBROUTINE NWS7GET()
      USE CONSTANTS, ONLY : RHOWAT0, G
      USE MESH, ONLY : X, Y, SLAM, SFEA, NP, ICS
      IMPLICIT NONE

      INTEGER I                     ! node loop counter
      INTEGER K, J                  ! latitude, longitude loop counters
      REAL(8) RHOWATG             ! ref. dens. of water * grav. constant
      INTEGER  LATIND1,LATIND2,LONIND1,LONIND2
      REAL(8) WLATM,WLONM,XWRATIO,YWRATIO
      REAL(8) XCOOR,YCOOR
      REAL(8), SAVE, ALLOCATABLE :: WVXFN(:,:),WVYFN(:,:),PRN(:,:)
      LOGICAL, SAVE :: MemoryAllocated = .False.

      call setMessageSource("nws7get")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

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
         IF(ICS.NE.1) THEN
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
         WVNX2(I) = WVXFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVXFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVXFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVXFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         WVNY2(I) = WVYFN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + WVYFN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + WVYFN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + WVYFN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         PRN2(I) = PRN(LATIND2,LONIND2)*XWRATIO*YWRATIO
     &        + PRN(LATIND2,LONIND1)*(1.d0-XWRATIO)*YWRATIO
     &        + PRN(LATIND1,LONIND2)*XWRATIO*(1.d0-YWRATIO)
     &        + PRN(LATIND1,LONIND1)*(1.d0-XWRATIO)*(1.d0-YWRATIO)
         PRN2(I) = PRN2(I)/RHOWATG
      END DO
C
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
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
      subroutine nws10Init()
      use global, only : ihot, ncice
      use mesh, only : np, slam, sfea
      implicit none
      real(8) :: gdlon
      integer :: j, jj ! loop counters

      call setMessageSource("nws10init")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      !
      ! gaussian met values
      allocate ( ug(latb*lonb),vg(latb*lonb),pg(latb*lonb) )
      if (ncice.eq.10) allocate ( ig(latb*lonb) )
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
      call g2rini()
      !
      ! Set the counter value. According to the ADCIRC documentation,
      ! fort.200 is read upon hotstart, but not at coldstart.
      if (ihot.eq.0) then
         nws10DatasetCounter = 1
      else
         nws10DatasetCounter = 0
      endif

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

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
      subroutine g2rini()
      use mesh, only : np, slam, sfea
      implicit none
      integer :: N,I,LON,LONP1,LAT,LATP1
      real(8) :: DLAT,DLON,FLONWORK,COLAT,DDLAT,XLAT,DFLAT,DFLAT1,
     &     DDLON,XLON,DFLON,DFLON1

      call setMessageSource("g2rini")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

C
C...Compute estimated DLAT, true DLON for Gaussian grid
      DLAT=PI/REAL(LATB-1,KIND=8)
      DLON=TWOPI/REAL(LONB,KIND=8)
      N=0

C...Loop through all the nodes in the grid to be interpolated onto and
C.....compute the interpolating factors.

      DO I=1,NP
C.....   Compute initial guess of which lon value FLON(I) is in the Gaussian file
C......  Check that this value is reasonable.
         flonwork=slam(i)
         if(flonwork.lt.0.) flonwork=flonwork+twopi
         lon=flonwork/dlon + 1
         lonp1=lon+1
         if(lon.eq.lonb) lonp1=1 !circle condition
         if((lon.lt.1).or.(lon.gt.lonb)) then
            call allMessage(ERROR,'Mesh node outside gaussian grid.')
            write(scratchMessage,
     &       '("i=",i0," lon=",f15.8," dlon=",f15.8," flon(i)=",f15.8)')
     &       i, lon, dlon, slam(i)
            call allMessage(ERROR, scratchMessage)
            call windTerminate()
         endif

C.....   Compute initial guess of which lat value FLAT(I) is in the Gaussian file
C....... Check that this value is reasonable.
         colat=hfpi-sfea(i)
         lat=colat/dlat + 1
         if(lat.eq.latb) lat=lat-1
         latp1=lat+1
         if((lat.lt.1).or.(lat.gt.latb)) then
            call allMessage(ERROR,'Mesh node outside gaussian grid.')
            write(scratchMessage,
     &       '("i=",i0," lat=",f15.8," dlat=",f15.8," flat(i)=",f15.8)')
     &       i, lat, dlat, sfea(i)
            call allMessage(ERROR, scratchMessage)
            call windTerminate()
         endif

         DO WHILE(.TRUE.)
            IF((COLAT.GE.GCLAT(LAT)).AND.(COLAT.LE.GCLAT(LATP1))) EXIT
            IF(COLAT.LT.GCLAT(LAT)) THEN
               LATP1=LAT
               LAT=LAT-1
               IF(LAT.LE.0) THEN
                  LAT=1
                  LATP1=2
                  EXIT
               ENDIF
               CYCLE
            ENDIF
            IF(COLAT.GT.GCLAT(LATP1)) THEN
               LAT=LAT+1
               LATP1=LAT+1
               IF(LAT.GE.LATB ) THEN
                  LAT=LATB-1
                  LATP1=LATB
                  EXIT
               ENDIF
               CYCLE
            ENDIF
            EXIT
         ENDDO
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

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
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
C   WJP 6/29/2020 added ability to read ice concentrations             *
C***********************************************************************

      subroutine nws10get(ull,vll,pll,ill)
      use constants, only : mb2pa
      use sizes, only : gblinputdir
      use mesh, only : np, slam, sfea
      implicit none
      real(8), intent(out) :: ull(:),vll(:),pll(:) ! mesh met values
      real(8), intent(out), optional :: ill(:) ! mesh ice values
      !
      real(8) :: rhowatg ! density of water * gravitational acceleration
      real(8) :: p1,p2,p3,p4 ! pressure at 4 corners of gaussian grid cell
      real(8) :: u1,u2,u3,u4 ! east velocity at 4 corners of gaussian grid cell
      real(8) :: v1,v2,v3,v4 ! north velocity at 4 corners of gaussian grid cell
      real(8) :: i1,i2,i3,i4 ! ice concentration at 4 corners of gaussian grid cell
      character(len=8) :: fname1 ! name of GFS file to read
      integer :: iext        ! numerical extension to file name
      integer :: ios = 0     ! i/o status for opening the fort.xxx files
      integer :: i, n        ! loop counters

      call setMessageSource("nws10get")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

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
!     WJP: With ice
      if (present(ill)) then
         do i=1,lonb*latb
            read(iext,*) pg(i),ug(i),vg(i),ig(i)
         enddo
      else
!     WJP: Without ice
         do i=1,lonb*latb
            read(iext,*) pg(i),ug(i),vg(i)
         enddo
      endif
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

!     WJP: Interpolate ice concentration if present 
      if (present(ill)) then
         do n=1,np
            ! grab the numerical values at the 4 corners of the grid 
            ! cell containing this mesh node
            i1=ig(n00(n))
            i2=ig(n10(n))
            i3=ig(n11(n))
            i4=ig(n01(n))
            ! spatially interpolate from the 4 values to the mesh vertex
            ill(n)=i1*d00(n)+i2*d10(n)+i3*d11(n)+i4*d01(n)
         enddo
      endif
C
      nws10DatasetCounter = nws10DatasetCounter + 1

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
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
      REAL(8) COLRAD(*),WGT(*),WGTCS(*),RCS2(*)
      INTEGER LGGHAF,L2,K,K1,ITER
      REAL(8) SI,SCALEVAL,RL2,DRAD,RAD,P1,P2,EPS,PHI,X,W,SN,RC
      real(8) :: dradz ! initial guess that represents a distance between gridpoints (radian)
      integer :: multiple ! yf51.52.25: reduce initial guess for DRAD for T1534 grid

      call setMessageSource("glats")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

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
         DO WHILE(.TRUE.)
            CALL POLY(L2,RAD,P2)
            DO WHILE(.TRUE.)
               P1 =P2
               ITER=ITER+1
               RAD=RAD+DRAD
               CALL POLY(L2,RAD,P2)
               IF(SIGN(SI,P1).EQ.SIGN(SI,P2)) CYCLE
               EXIT
            ENDDO
            IF(DRAD.LT.EPS) EXIT
            RAD=RAD-DRAD
            DRAD = DRAD * 0.25d0
         ENDDO
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
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
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
      REAL(8) RAD,P,X,Y1,Y2,Y3,G
C
      X = COS(RAD)
      Y1 = 1.0d0
      Y2=X
      DO 1 I=2,N
         G=X*Y2
         Y3=G-Y1+G-(G-Y1)/REAL(I,KIND=8)
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

      SUBROUTINE NWS11GET()
      USE SIZES
      USE CONSTANTS, ONLY : RHOWAT0, G
      USE MESH, ONLY : NP, SLAM, SFEA
      IMPLICIT NONE
      INTEGER,SAVE  ::  ICALL = 0
      INTEGER I,IEXT,IDIG1,IDIG2,IDIG3,KERR,N
      INTEGER IYEAR,IMONTH,IDAY,IHOUR
      REAL*8 RHOWATG100,FLONDEG,FLATDEG
      REAL(8) P1,P2,P3,U1,U2,U3,V1,V2,V3,UE29,VE29,CBETAU,SBETAU
C
      INTEGER,SAVE,ALLOCATABLE ::  N1(:),N2(:),N3(:)
      REAL(8),SAVE,ALLOCATABLE :: D1(:),D2(:),D3(:),BETAU(:)
      REAL(8),SAVE,ALLOCATABLE :: UE(:),VE(:),PE(:)
C
      CHARACTER*1 FNAME2(8)
      CHARACTER*8 FNAME1
      EQUIVALENCE (FNAME1,FNAME2)
      integer :: errorio

      call setMessageSource("nws11get")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

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
         call allMessage(INFO,'Computing ETA29 met field interp factors.')
         DO I=1,NP
            flondeg=rad2deg*slam(i)
            flatdeg=rad2deg*sfea(i)
            CALL E29SEARCH(I,FLONDEG,FLATDEG,N1(I),N2(I),N3(I),
     &           D1(I),D2(I),D3(I),betau(i))
         END DO
#if defined(WIND_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
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
      call openFileForRead(IEXT,FNAME1,errorIO)
      if (errorIO.ne.0) then
         call windTerminate()
      endif

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
         WVNX2(N)=UE29*CBETAU - VE29*SBETAU
         WVNY2(N)=UE29*SBETAU + VE29*CBETAU
         PRN2(N)=P1*D1(N)+P2*D2(N)+P3*D3(N)
         PRN2(N)=PRN2(N)/RHOWATG100
      END DO
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE NWS11GET
C-----------------------------------------------------------------------



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
      real(8) dd1,dd2,dd3,betau,ri,x1,x2,x3,x4,y1,y2,y3,y4
      real(8) aemin,areas,a1,a2,a3,aa,ae,lambda
      real(8) lamda0,phi0,rphi0,cphi0,sphi0,tphi0,dlamda,dphi,rdlamda,
     &       rdphi,rflat,tflat,sflat,cflat,a,rlamar,cphiicrlamda,phiarg,
     &       rphii,rlamda,ri1,ri2,rj,dgtora,flon,flat
      real(8) lamda,lamdaa,lamdab,lamdac,lamdad,lamdae,lamdag
      real(8) phi,phia,phib,phic,phid,phie,phig
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
      real(8) lamda,phi,phii,dlon,dlat,dlnt,arg,betau1,betau2,betau3
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
      real(8) dd1,dd2,dd3,betau
      real(8) lamda,lamdap1,lamdam1,phi,phip1,phim1,dlon,dlat,
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
      SUBROUTINE HollandGet(WVNX,WVNY,PRESS,TIMELOC)
      USE SIZES, ONLY : MyProc, LOCALDIR
      USE MESH, ONLY : X, Y, SLAM, SFEA, NP, ICS
      USE CONSTANTS, ONLY : RHOWAT0, G, mb2pa
      IMPLICIT NONE
      REAL(8), intent(in) :: TIMELOC
      REAL(8), intent(out), dimension(NP) :: WVNX,WVNY
      REAL(8), intent(out), dimension(NP) :: PRESS
      INTEGER I, J
      REAL(8),SAVE,ALLOCATABLE :: RAD(:),DX(:),DY(:)
      REAL(8),SAVE,ALLOCATABLE :: XCOOR(:),YCOOR(:) !lat and lon of nodes
      REAL(8),SAVE,ALLOCATABLE :: V_r(:)            !wind sp. as fn. of dist
      REAL(8),SAVE,ALLOCATABLE :: THETA(:)
      REAL(8) :: TVX,TVY,RRP,RMW,A,B,WTRATIO
      REAL(8) :: TransSpdX, TransSpdY
      REAL(8) :: TM, cpress, lon, lat, spd
      REAL(8) :: ts        ! storm translation speed, m/s
      REAL(8) :: WindMultiplier         ! for storm 2 in LPFS ensemble
      REAL(8) :: centralPressureDeficit ! difference btw ambient and cpress
C
      REAL(8) :: omega, coriolis
      REAL(8) :: mperdeg

      LOGICAL, SAVE :: FIRSTCALL = .True.

      call setMessageSource("hollandget")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

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
      CALL GetHollandStormData(lat,lon,cpress,spd,rrp,rmw,
     &                         tvx,tvy,TIMELOC)
      !
      ! @jasonfleming: If this is a "CALM" period, set winds to zero
      ! velocity and pressure equal to the background pressure and return.
      if ( cpress.lt.0.d0 ) then
         press(:) = PRBCKGRND * mb2pa / rhoWat0g ! convert mb to mH2O
         wvnx(:)=0.d0
         wvny(:)=0.d0
         EyeLat(1) = EyeLat(2)
         EyeLon(1) = EyeLon(2)
         EyeLat(2) = EyeLat(3)
         EyeLon(2) = EyeLon(3)
         EyeLat(3) = lat
         EyeLon(3) = lon
         FoundEye = .true.
#if defined(WIND_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         return
      endif
C
C     jgf50.32: If we are using sector-based wind drag, record the location
C     of the center of the storm.
      IF ((lat.ne.EyeLat(3)).or.(lon.ne.EyeLon(3))) THEN
         EyeLat(1) = EyeLat(2)
         EyeLon(1) = EyeLon(2)
         EyeLat(2) = EyeLat(3)
         EyeLon(2) = EyeLon(3)
         EyeLat(3) = lat
         EyeLon(3) = lon
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
C     RJW v48.45
C     Compute the distances based on haversine formula for distance along a sphere
         rad(i) = sphericalDistance(DX(I),DY(I),lat, ycoor(i))
C     Compute coriolis
         coriolis = 2.0d0 * omega * sin(YCOOR(I)*DEG2RAD)
C     Compute pressure field
!P.V 11/04/2022 to account for Southern Hemisphere
!         PRESS(I)=(cpress+(centralPressureDeficit)*
!     &        EXP(-(RMW*1000.d0/RAD(I))**B)) / (RHOWAT0*G)
         PRESS(I)=cpress+(centralPressureDeficit)*
     &        EXP(-(RMW*1000.d0/RAD(I))**B)
C     Compute velocity field with absolute around 
C     CORIOLIS for the Southern Hempisphere
         V_r(I) = sqrt(
     &        (RMW*1000.d0/RAD(I))**B *
     &        EXP(1.d0-(RMW*1000.d0/RAD(I))**B)*spd**2.d0
     &        + (RAD(I)**2.d0)*(CORIOLIS**2.d0)/4.d0
     &        )
     &        - RAD(I)*ABS(CORIOLIS)/2.d0
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
!P.V 11/04/2022 to account for Southern Hemisphere
        !WVNX(I)=-V_r(I)*SIN(THETA(I))
        !WVNY(I)= V_r(I)*COS(THETA(I))
        ! Find the wind velocity components (caution to SH/NH)
        if(lat.lt.0.d0) then ! SH
         WVNX(I)= V_r(I)*SIN(THETA(I))
         WVNY(I)=-V_r(I)*COS(THETA(I))
        else ! NH
         WVNX(I)=-V_r(I)*SIN(THETA(I))
         WVNY(I)= V_r(I)*COS(THETA(I))
        endif
 
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
!P.V 11/04/2022
        PRESS(I) = max(0.85d5,min(1.1e5,PRESS(I))) !Typhoon Tip 870 hPa ... 12-oct-1979
        WVNX(I)  = max(-200.d0,min(200.d0,WVNX(I)))
        WVNY(I)  = max(-200.d0,min(200.d0,WVNY(I)))
        !-------------------------------------------
        ! Convert atmospheric pressure (Pascals) to
        ! atmospheric pressure-induced water surface
        ! elevation (meters).
        !-------------------------------------------
        PRESS(I) = PRESS(I) / (RHOWAT0*G)

      END DO
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
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
     &     RRPOut,RMWOut,TVXOut,TVYOut,TIMELOC)
      USE SIZES, ONLY : MyProc, GBLINPUTDIR
      USE GLOBAL, ONLY : RNDAY
      USE VORTEX, ONLY : uvtrans
      IMPLICIT NONE
      REAL(8),intent(in) :: TIMELOC
      REAL(8),intent(out) :: LatOut, LonOut, CPressOut
      REAL(8),intent(out) :: SpdOut, RRPOut, RMWOut, TVXOut, TVYOut
C
      INTEGER, ALLOCATABLE, SAVE :: iYear(:),iMth(:),iDay(:),iHr(:)
      INTEGER, ALLOCATABLE, SAVE :: iLat(:),iLon(:)
      CHARACTER(1), ALLOCATABLE, SAVE :: ns(:), ew(:)
      INTEGER, ALLOCATABLE, SAVE :: iSpd(:),iCPress(:),iRRP(:),iRMW(:)
      REAL(8), ALLOCATABLE, SAVE :: Lat(:),Lon(:),Spd(:)
      REAL(8), ALLOCATABLE, SAVE :: CPress(:),RRP(:),RMW(:)
      REAL(8), ALLOCATABLE, SAVE :: TVX(:), TVY(:)      ! jgf46.32jgf9
      CHARACTER(len=4), ALLOCATABLE, SAVE :: CastType(:) !hindcast,forecast
      INTEGER, ALLOCATABLE, SAVE :: iFcstInc(:) ! hours
      REAL(8), ALLOCATABLE, SAVE :: FcstInc(:) ! seconds
      REAL(8),SAVE :: WTRATIO          !time ratio used for interpolation
      REAL(8), ALLOCATABLE, SAVE :: CastTime(:) ! seconds since start of year
      INTEGER,SAVE :: iNowcastCPress, iNowcastRRP, iNowcastRMW
      LOGICAL,SAVE :: FIRSTCALL = .True.
      INTEGER,SAVE :: i  ! Current array counter for fort.22 file
      INTEGER,SAVE :: nl ! Number of lines in the fort.22 file
      INTEGER,SAVE :: pl ! populated length of Holland Data array
      INTEGER :: j ! loop counter
      INTEGER :: ios ! return code for an i/o operation
      call setMessageSource("getHollandStormData")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

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
     &              CastTime(i))
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
     &                 CastTime(i))
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
                  call allMessage(ERROR,
     &            'The storm hindcast/forecast input file (unit 22) '//
     &            'contains invalid data for central pressure or Rmax.')
                  call windTerminate()
               endif
C
      ! @jasonfleming: Adding a new type to allow the analyst to add lines
      ! that do nothing but produce zero winds and background barometric
      ! pressure. These lines can have a date/time like a BEST line or
      ! a date/time and forecast period like an OFCL line.
            CASE("CALM")
               call allMessage(ECHO,
     &            'The fort.22 file contains at least one "CALM" line.')
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

               IF ( iFcstInc(i).eq.0 ) THEN
                  CALL TimeConv(iYear(i),iMth(i),iDay(i),iHr(i),0,0.d0,
     &                 CastTime(i))
                  CastTime(i)=CastTime(i)-WindRefTime
               ELSE
                  FcstInc(i) = FcstInc(i) * 3600.d0 ! convert hours to seconds
                  CastTime(i) = CastTime(i-1) +
     &                 ( FcstInc(i) - FcstInc(i-1) )
               ENDIF

            CASE DEFAULT        ! unrecognized
               call allMessage(ERROR,
     &          'Only "BEST", "OFCL", or "CALM" are allowed '//
     &          'in the 5th column of fort.22.')
               call windTerminate()
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
               write(scratchMessage,2355) FcstInc(i-1)
               call logMessage(DEBUG, scratchMessage)
            endif
 1244       FORMAT('DEBUG HOLLAND: CastTime ',
     &           e16.8,' Lat ',f6.2,' Lon ',f6.2,
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
            IF (ceiling(TIMELOC).ge.CastTime(i-1).and.
     &          floor(TIMELOC).lt.CastTime(i)) THEN
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
      !
      ! @jasonfleming: Check to see if there is enough data to cover
      ! the whole run and bomb out immediately if there isn't.
      if (castTime(pl).lt.RNDAY*86400.d0) then
         call allMessage(ERROR,'The fort.22 file ends before RNDAY.')
         call windTerminate()
      endif

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
C        jgf51.14: Check to see that we haven't gone off the end of
C        meteorological data.
          IF (i.gt.pl) THEN
             call allMessage(ERROR,'The simulation time has extended '//
     &          'beyond the end of the meteorological dataset.')
             call windTerminate()
          ENDIF
      ENDIF
C
C     Interpolate w.r.t. time
      ! @jasonfleming: Add CALM handling.
      LonOut = Lon(i-1) + WTRATIO*(Lon(i)-lon(i-1))
      LatOut = Lat(i-1) + WTRATIO*(Lat(i)-lat(i-1))
      if ( (trim(castType(i)).ne.'CALM').and.
     &      trim(castType(i-1)).ne.'CALM' ) then
         WTRATIO=(TIMELOC-CastTime(i-1))/(CastTime(i)-CastTime(i-1))
         CPressOut = CPress(i-1) + WTRATIO*(CPress(i)-CPress(i-1))
         SpdOut = Spd(i-1) + WTRATIO*(Spd(i)-Spd(i-1))
         RRPOut = RRP(i-1) + WTRATIO*(RRP(i)-RRP(i-1))
         RMWOut = RMW(i-1) + WTRATIO*(RMW(i)-RMW(i-1))
         TVXOut = TVX(i-1) + WTRATIO*(TVX(i)-TVX(i-1))
         TVYOut = TVY(i-1) + WTRATIO*(TVY(i)-TVY(i-1))
      else
         CPressOut = -99999.d0
         SpdOut = -99999.d0
         RRPOut = -99999.d0
         RMWOut = -99999.d0
         TVXOut = -99999.d0
         TVYOut = -99999.d0
      endif
C
C     ----------------------------------------------------------
C     END Code executed on every call to this subroutine
C     ----------------------------------------------------------
 10   format(8x,i4,i2,i2,i2,6x,a4,7x,i3,4x,i3,3x,i3,2x,i4,52x,i3)
 12   format(8x,i4,i2,i2,i2,6x,a4,2x,i3,2x,i3,4x,i3,3x,i3,2x,i4,52x,i3)
 14   format(8x,i4,i2,i2,i2,6x,a4,2x,i3,1x,i4,3x,i4,3x,i3,2x,i4,47x,i3,
     &     2x,i3)

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()

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
      REAL(8) RHOWAT0,G,WTIME1,WTIME2,TIMELOC,WTIMINC
      REAL(8) X(*),Y(*),SLAM(*),SFEA(*)
      REAL(8) WVNX(*),WVNY(*),PRESS(*)
! LOCAL VARIABLES
      REAL(8),SAVE,ALLOCATABLE :: RAD(:),DX(:),DY(:),XCOOR (:),YCOOR(:)
      REAL(8),SAVE,ALLOCATABLE :: EARG(:)
      REAL(8),SAVE,ALLOCATABLE :: VG(:),THETA(:)
      REAL(8),SAVE :: XC1,YC1,DP1,TVX1,TVY1,RMW1,A1,B1
      REAL(8),SAVE :: XC2,YC2,DP2,TVX2,TVY2,RMW2,A2,B2
      REAL(8) :: XC,YC,DP,TVX,TVY,RMW,A,B,WTRATIO
      REAL(8) :: Central_Pressure,T1,T2,T3,D,N,TT,TM
      REAL(8),SAVE :: RHOWATG
      IF (FIRSTCALL.EQ.0) THEN
         FIRSTCALL = 1
         RHOWATG=RHOWAT0*G
         ALLOCATE (RAD(NP),DX(NP),DY(NP),XCOOR(NP),YCOOR(NP),EARG (NP))
         ALLOCATE (VG(NP),THETA(NP))
!         IF(ICS.NE.1) THEN
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
      !    2/19/19     Taylor A.                   Changed ilon to 4 columns to facilitate
      !                                            distant storms.
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
      SUBROUTINE NWS19GET(wvnx,wvny,press,timeloc)
      USE CONSTANTS, ONLY: rhowat0, g, mb2pa
      USE GLOBAL,ONLY: rnday
      USE MESH, ONLY : np, ics, slam, sfea
      USE VORTEX, ONLY : newVortex, setRmaxes, fitRmaxes, uvp,
     &                   getLatestRmax, getLatestAngle
      IMPLICIT NONE
      INTEGER  :: ient,itpc
      INTEGER  :: i, nwi
      INTEGER  :: num_entry
      REAL(8) :: TIMELOC, time1, wtratio
      REAL( 8), SAVE, ALLOCATABLE :: xcoor(:),ycoor(:)
      REAL( 8) :: wvnx(np),wvny(np), press(np)
      REAL(8) :: RhoWatG

      REAL,DIMENSION(:),ALLOCATABLE,SAVE::lat,lon
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE::ilat,ilon,atcfRMW
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: ilat1,ilon1,atcfRMW1
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: ispd,icpress
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: ispd1,icpress1
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: iyear,imth,iday,ihr
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE ::iyear1,imth1,iday1,ihr1
      REAL(8), SAVE ::  lat1, lon1,  spd1,  cpress1,  rmw1
      REAL(8), SAVE ::  dt, told, cLat1, cLon1
      CHARACTER(LEN=4),DIMENSION(:),ALLOCATABLE, SAVE:: type! hindcast/nowcast or forecast?
      CHARACTER(LEN=4),DIMENSION(:),ALLOCATABLE, SAVE:: type1 !hindcast/nowcast or forecast?
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE::iFcstInc1, iFcstInc ! hours between forecasts
      INTEGER , SAVE :: firstCall = 0
      LOGICAL , SAVE :: firstTime = .TRUE.
      REAL(8),DIMENSION(:),ALLOCATABLE,SAVE :: HollB,HollB1
      INTEGER ,DIMENSION(:),ALLOCATABLE,SAVE :: advr,advr1
      REAL(8), SAVE :: Vmax
      REAL(8), SAVE :: Vr
      REAL(8), SAVE, DIMENSION(4) :: r
      REAL(8), SAVE :: Pn
      REAL(8), SAVE :: Pc
      REAL(8), SAVE :: cLat
      REAL(8), SAVE :: cLon
      REAL(8), SAVE :: cHollB
      REAL(8) , DIMENSION(4), SAVE :: crmaxw
      REAL(8), SAVE :: timeOffset
!-------------------------------
      INTEGER,SAVE :: icyc
!------------------------------
      INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: cycle_num
      INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: isotachs_per_cycle
      INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: isotachs_per_cycle1
      INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: ipn,ipn1,cycle_num1

      REAL(8) , DIMENSION(:,:), ALLOCATABLE,SAVE :: rmaxw,rmaxw1
      REAL(8) , DIMENSION(:), ALLOCATABLE,SAVE::CycleTime,CycleTime1
      REAL(8) , DIMENSION(:), ALLOCATABLE,SAVE::uTrans1,uTrans
      REAL(8) , DIMENSION(:), ALLOCATABLE,SAVE::vTrans1,vTrans
      REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: h_speed
      INTEGER , DIMENSION(:,:), ALLOCATABLE,SAVE::quadflag,quadflag1
      INTEGER ,DIMENSION(:),ALLOCATABLE,SAVE :: ivr, dir,speed
      REAL(8), SAVE :: dx,dy
      REAL(8) :: VmaxBL
      REAL(8) :: stormMotionU
      REAL(8) :: stormMotionV
      REAL(8) :: stormMotion
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: ivr1, dir1, speed1
      INTEGER , DIMENSION(:,:),ALLOCATABLE, SAVE :: ir,ir1
      CHARACTER(LEN = 10),DIMENSION(:),ALLOCATABLE,SAVE::name,name1
      REAL(8),SAVE :: timeOld,timeNew
      REAL(8) :: uTransNow, vTransNow ! time-interpolated overland speed, kts
      CHARACTER(1) ew,ns
      INTEGER, SAVE :: num_cycles ! num unique date/times in the fort.22
      REAL(8) :: rndaySeconds ! end time of the run in seconds

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
         timeOffset = nint(TIMELOC/60d0)*60d0 ! round to nearest minute
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
         ! @jasonfleming: compute run end time in seconds to make sure
         ! there are enough data in the file to cover the whole run.
         rndaySeconds = RNDAY * 86400.d0
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

26       FORMAT(3x, i3, 2x, i4, 3i2, 6x, a4, 2(2x,i3),a1, 2x, i4,a1, 2x,
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

         !
         ! @jasonfleming: If there isn't enough data in the fort.22 to
         ! cover the whole run, bomb out immediately.
         if (cycleTime(num_cycles).lt.rndaySeconds) then
            call allMessage(ERROR,"There aren't enough data in the "//
     &      'fort.22 file to cover the entire run time.')
            call windTerminate()
         endif

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
      !
      ! @jasonfleming: Add handling of "CALM"
      if ((trim(type(icyc)).eq.'CALM').or.trim(type(icyc-1)).eq.'CALM') then
         wvnx(:)=0.d0
         wvny(:)=0.d0
         press(:)=PRBCKGRND*mb2pa/(rhowat0*g)
         wtratio=(TIMELOC-CycleTime(icyc-1))/
     &                  (CycleTime(icyc)-CycleTime(icyc-1))
         cLat   = lat(icyc-1)  +
     &                         wtratio * (lat(icyc)-lat(icyc-1))
         cLon  =  lon(icyc-1)  +
     &                         wtratio * (lon(icyc)-lon(icyc-1))
         vortexLat = cLat
         vortexLon = cLon
         vortexRMW = 0.0d0
         EyeLat(1) = EyeLat(2)
         EyeLon(1) = EyeLon(2)
         EyeLat(2) = EyeLat(3)
         EyeLon(2) = EyeLon(3)
         EyeLat(3) = cLat
         EyeLon(3) = cLon
         FoundEye = .true.
#if defined(WIND_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         return
      endif
      !
      ! Interpolate NHC forecast interval in
      ! time to obtain values at model time.
      wtratio=(TIMELOC-CycleTime(icyc-1))/
     &                  (CycleTime(icyc)-CycleTime(icyc-1))
      !
      ! Perform time interpolation, transform variables from integers
      ! to real numbers for hurricane vortex calcualtions.
      cLat   =  lat(icyc-1)  +
     &                         wtratio * (lat(icyc)-lat(icyc-1))
      cLon  =  lon(icyc-1)  +
     &                         wtratio * (lon(icyc)-lon(icyc-1))
      Vmax  =  1.d0*(ispd(icyc-1) +
     &                         wtratio * (ispd(icyc)-ispd(icyc-1)))
      Pn    =  1.d0*(ipn(icyc-1)  +
     &                         wtratio * (ipn(icyc)-ipn(icyc-1)))
      Pc    =  1.d0*(icpress(icyc-1) +
     &                         wtratio*(icpress(icyc)-icpress(icyc-1)))
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
      IF ((cLat.ne.EyeLat(3)).or.(cLon.ne.EyeLon(3))) THEN
         EyeLat(1) = EyeLat(2)
         EyeLon(1) = EyeLon(2)
         EyeLat(2) = EyeLat(3)
         EyeLon(2) = EyeLon(3)
         EyeLat(3) = cLat
         EyeLon(3) = cLon
         FoundEye = .true.
      ENDIF
      !
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
      !    2/19/19     Taylor A.                   Changed ilon to 4 columns to facilitate
      !                                            distant storms.
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
      SUBROUTINE NWS20GET(wvnx, wvny, press, timeloc)
      USE VORTEX, ONLY : setVortex, fitRmaxes4,  uvpr, spInterp,
     &                   Rmaxes4, QuadFlag4, QuadIr4, Bs4, VmBL4  ! jgfdebug
      USE CONSTANTS, ONLY: kt2ms, nm2m, deg2rad, omega, Rearth, m2nm, mb2pa
      USE GLOBAL, ONLY: rnday
      USE MESH, ONLY : slam, sfea, np, ics
      IMPLICIT NONE

      INTEGER  :: ient
      INTEGER  :: i, j
      INTEGER  :: num_entry
      REAL(8) :: TIMELOC, wtratio
      REAL( 8), SAVE, ALLOCATABLE :: xcoor(:),ycoor(:)
      REAL( 8) :: wvnx(np),wvny(np), press(np)
      REAL(8) :: RhoWatG
      REAL(8) :: corio                       ! Coriolis force (1/s)

      !REAL(8),DIMENSION(:),ALLOCATABLE,SAVE::lat,lon
      ! jgf53.dev: FIXME this should be REAL(8) as above
      REAL(8),DIMENSION(:),ALLOCATABLE,SAVE::lat,lon
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE::ilat,ilon,atcfRMW
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: ilat1,ilon1,atcfRMW1
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: ispd,icpress
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: ispd1,icpress1
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: iyear,imth,iday,ihr
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE ::iyear1,imth1,iday1,ihr1
      CHARACTER(LEN=4),DIMENSION(:),ALLOCATABLE, SAVE:: type! hindcast/nowcast or forecast?
      CHARACTER(LEN=4),DIMENSION(:),ALLOCATABLE, SAVE:: type1 !hindcast/nowcast or forecast?
      INTEGER,DIMENSION(:),ALLOCATABLE,SAVE::iFcstInc1, iFcstInc ! hours between forecasts
      INTEGER , SAVE :: firstCall = 0
      LOGICAL , SAVE :: firstTime = .TRUE.
      REAL(8),DIMENSION(:),ALLOCATABLE,SAVE :: HollB,HollB1
      REAL(8),DIMENSION(:,:),ALLOCATABLE,SAVE :: VmaxesBL1,HollBs1
      REAL(8),DIMENSION(:,:,:),ALLOCATABLE,SAVE :: HollBs,VmaxesBL
      INTEGER ,DIMENSION(:),ALLOCATABLE,SAVE :: advr,advr1
      REAL(8), SAVE :: Vr
      REAL(8), SAVE, DIMENSION(4) :: r
      REAL(8), SAVE :: Pn
      REAL(8), SAVE :: Pc
      REAL(8), SAVE :: cLat
      REAL(8), SAVE :: cLon
      REAL(8), SAVE :: told
      ! save parameters for each node
      ! change them to allocatable, save Aug. 2013, Jie
      REAL(8), ALLOCATABLE, SAVE :: cHollBs1(:), cHollBs2(:)
      REAL(8), ALLOCATABLE, SAVE :: cPhiFactor1(:), cPhiFactor2(:)
      REAL(8), ALLOCATABLE, SAVE :: cVmwBL1(:), cVmwBL2(:)
      REAL(8), ALLOCATABLE, SAVE :: crmaxw1(:), crmaxw2(:)
      REAL(8), ALLOCATABLE, SAVE :: crmaxwTrue1(:), crmaxwTrue2(:)
      REAL(8), ALLOCATABLE, SAVE :: cHollBs(:)
      REAL(8), ALLOCATABLE, SAVE :: cPhiFactor(:)
      REAL(8), ALLOCATABLE, SAVE :: cVmwBL(:)
      REAL(8), ALLOCATABLE, SAVE :: crmaxw(:)
      REAL(8), ALLOCATABLE, SAVE :: crmaxwTrue(:)
      REAL(8), ALLOCATABLE, SAVE :: aziangle(:), dist(:)
      REAL(8), SAVE :: timeOffset
!-------------------------------
      INTEGER,SAVE :: icyc
      INTEGER,SAVE :: isot
!------------------------------
      INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: cycle_num
      INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: isotachs_per_cycle
      INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: isotachs_per_cycle1
      INTEGER ,DIMENSION(:),ALLOCATABLE, SAVE :: ipn,ipn1,cycle_num1

      REAL(8) , DIMENSION(:,:), ALLOCATABLE,SAVE :: rmaxw1
      REAL(8) , DIMENSION(:,:,:), ALLOCATABLE,SAVE :: rmaxw
      REAL(8) , DIMENSION(:), ALLOCATABLE,SAVE::CycleTime,CycleTime1
      REAL(8) , DIMENSION(:), ALLOCATABLE,SAVE::uTrans1,uTrans
      REAL(8) , DIMENSION(:), ALLOCATABLE,SAVE::vTrans1,vTrans
      REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: h_speed
      INTEGER , DIMENSION(:,:), ALLOCATABLE,SAVE::quadflag1
      INTEGER , DIMENSION(:,:,:), ALLOCATABLE,SAVE::quadflag
      INTEGER ,DIMENSION(:),ALLOCATABLE,SAVE :: ivr, dir,speed
      REAL(8), SAVE :: dx,dy
      REAL(8) :: VmaxBL
      REAL(8) :: stormMotionU
      REAL(8) :: stormMotionV
      REAL(8) :: stormMotion
      INTEGER,DIMENSION(:),ALLOCATABLE, SAVE :: ivr1, dir1, speed1
      INTEGER , DIMENSION(:,:),ALLOCATABLE, SAVE :: ir1
      INTEGER , DIMENSION(:,:,:),ALLOCATABLE, SAVE :: ir
      CHARACTER(LEN = 10),DIMENSION(:),ALLOCATABLE,SAVE::name,name1
      REAL(8),SAVE :: timeNew
      REAL(8) :: uTransNow, vTransNow ! time-interpolated overland speed, kts
      REAL(8) :: dirNow ! Jie
      REAL(8) :: numRmaxes ! number of Rmax values selected by user
      CHARACTER(1) ew,ns
      INTEGER, SAVE :: num_cycles ! num unique date/times in the fort.22

      REAL(8) :: rndaySeconds ! end time of the run in seconds
      real(8) :: newVortexRMW ! average Rmax at icyc   (next fort.22 time value)
      real(8) :: oldVortexRMW ! average Rmax at icyc-1 (previous fort.22 time value)

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
         
         timeOffset = nint(TIMELOC/60d0)*60d0 ! round initial time to nearest minute
         rndaySeconds = RNDAY * 86400.d0

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
         ir = 0
         rmaxw = 0.0d0
         quadflag = 0
         HollBs = 0.0d0
         VmaxesBL = 0.0d0

         ! Read parsed NHC advisory in modified best-track format

27       FORMAT(3x, i3, 2x, i4, 3i2, 6x, a4, 2(2x,i3),a1, 2x, i4,a1, 2x,
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
               if (isot.gt.4) then
                  write(scratchMessage,'(a,a,i0,a)')
     &             'The GAHM fort.22 input file has more than 4 ',
     &             'isotachs in cycle ',cycle_num1(ient),'.'
                  call allMessage(ERROR,scratchMessage)
                  call windTerminate()
               endif
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
         !
         ! @jasonfleming: If there isn't enough data in the fort.22 to
         ! cover the whole run, bomb out immediately.
         write(scratchMessage,*) 'The last cycle in the fort.22 '//
     &      'file is at ', cycleTime(num_cycles), ' [sec]'
         CALL logMessage(INFO,scratchMessage)
         if (cycleTime(num_cycles).lt.rndaySeconds) then
            call allMessage(ERROR,"There aren't enough data in the "//
     &      'fort.22 file to cover the entire run time.')
            call windTerminate()
         endif

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
      !
      ! @jasonfleming: Add handling of "CALM"
      if ((trim(type(icyc)).eq.'CALM').or.trim(type(icyc-1)).eq.'CALM') then
         wvnx(:)=0.d0
         wvny(:)=0.d0
         press(:)=PRBCKGRND*mb2pa/(rhowat0*g)
         wtratio=(TIMELOC-CycleTime(icyc-1))/
     &                  (CycleTime(icyc)-CycleTime(icyc-1))
         cLat   = lat(icyc-1)  +
     &                         wtratio * (lat(icyc)-lat(icyc-1))
         cLon  =  lon(icyc-1)  +
     &                         wtratio * (lon(icyc)-lon(icyc-1))
         vortexLat = cLat
         vortexLon = cLon
         vortexRMW = 0.0d0
         EyeLat(1) = EyeLat(2)
         EyeLon(1) = EyeLon(2)
         EyeLat(2) = EyeLat(3)
         EyeLon(2) = EyeLon(3)
         EyeLat(3) = cLat
         EyeLon(3) = cLon
         FoundEye = .true.
#if defined(WIND_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG,"Return.")
#endif
         call unsetMessageSource()
         return
      endif

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
      !P.V 11/04/2022
      !Using absolute value for coriolis for Southern Hemisphere
      corio = ABS(2.0d0 * omega * SIN(deg2rad*cLat))

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
      IF ((cLat.ne.EyeLat(3)).or.(cLon.ne.EyeLon(3))) THEN
         EyeLat(1) = EyeLat(2)
         EyeLon(1) = EyeLon(2)
         EyeLat(2) = EyeLat(3)
         EyeLon(2) = EyeLon(3)
         EyeLat(3) = cLat
         EyeLon(3) = cLon
         FoundEye = .true.
      ENDIF

      ! jgf: arithmetic mean of the Rmax in the selected storm quadrants from
      ! the selected isotachs
      numRmaxes = 0.d0 ! count how many Rmax values we have to work with
      do i=1,4  ! loop quads
         do j=1,3 ! loop isotachs
            if (quadflag(icyc-1,i,j)==1) then
               oldVortexRMW = oldVortexRMW + rmaxw(icyc-1,i,j)
               numRmaxes = numRmaxes + 1.d0
            endif
         end do
      enddo
      oldVortexRMW = oldVortexRMW / numRmaxes
      numRmaxes = 0.d0 ! count how many Rmax values we have to work with
      do i=1,4  ! loop quads
         do j=1,3 ! loop isotachs
            if (quadflag(icyc,i,j)==1) then
               newVortexRMW = newVortexRMW + rmaxw(icyc,i,j)
               numRmaxes = numRmaxes + 1.d0
            endif
         end do
      enddo
      newVortexRMW = newVortexRMW / numRmaxes
      vortexRMW = oldVortexRMW + wtratio * ( newVortexRMW - oldVortexRMW )

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


         !P.V 11/04/2022
         press(i) = max(0.85d5,min(1.1e5,press(i)))  ! Typhoon Tip 870 hPa ... 12-oct-1979
         wvnx(i)  = max(-200.d0,min(200.d0,wvnx(i)))
         wvny(i)  = max(-200.d0,min(200.d0,wvny(i)))

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
     &    '("There were ",i0," HWind files specified.")') numFiles
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
     &    '("Cycle time ",i0," must be less cycle time ",i0,".")')i,i+1
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
      EyeLat = 0.d0
      EyeLon = 0.d0

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
      use constants, only : mb2pa, m2nm, nm2m
      use mesh, only : slam, sfea, np, dp
      implicit none
      real(8), intent(out) :: wvnx(:)  ! wind u-velocity (m/s)
      real(8), intent(out) :: wvny(:)  ! wind v-velocity (m/s)
      real(8), intent(out) :: press(:) ! barometric pressure (mH2O)
      real(8), intent(in) :: timeloc    ! adcirc time, seconds since coldstart
C
      real(8) :: wtratio  ! time weighting factor for time interpolation
      real(8) :: cLonNow  ! time interpolated longitude at center of storm
      real(8) :: cLatNow  ! time interpolated latitude at center of storm
      real(8) :: rampValNow  ! time interpolated met ramp value
      real(8) :: k_0      ! constant used in mercator projection
      real(8) :: y_mercator_0 ! mercator y coord at center of storm
      real(8) :: rmaxNow  ! (nm) used in estimating central pressure
      real(8) :: pc       ! (mbar) used to find pressure field
      real(8) :: B        ! Holland cyclone shape parameter
      real(8) :: dist     ! (m) distance from storm center to mesh node
      real(8) :: dxNow,dyNow ! (nm) x and y distance from storm center rmaxNow
      real(8) :: vmaxNow  ! (m/s) time interpolated fulldomain vmax
      real(8) :: oldDX, oldDY ! (nm) distances to rmax in last hwind cycle
      real(8) :: newDX, newDY ! (nm) distances to rmax in next hwind cycle
      real(8) :: angle    ! (nm) angle to location of mesh vmax (on subdomain in parallel)
      real(8) :: maxSpeed ! (m/s) max wind speed on mesh (subdomain in parallel)
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
      if ( (cLatNow.ne.EyeLat(3)).or.(cLonNow.ne.EyeLon(3))) then
         EyeLat(1) = EyeLat(2)
         EyeLon(1) = EyeLon(2)
         EyeLat(2) = EyeLat(3)
         EyeLon(2) = EyeLon(3)
         EyeLat(3) = cLatNow
         EyeLon(3) = cLonNow
         FoundEye = .TRUE.
      endif
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
      use global, only : openFileForRead, nabout
      use constants, only: m2nm
      implicit none
      integer, intent(in) :: cn ! cycle number of file to load
C
      real(8), allocatable    :: windMag(:,:) ! hwind speeds, used to find hwind vmax
      real(8)                 :: dx, dy ! distances (m) to location of vmax
      integer                  :: maxIndices(2) !i,j indices of hwind vmax
      real(8)                 :: throwAwayCoordinate ! dummy variable
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
     &    '("Could not read line ",i0," from the file ",A,".")')
     &    currentLine, trim(hWindFiles(cn)%file_name)
         call allMessage(ERROR,scratchMessage)
      else
         write(scratchMessage,
     &     '("Could not read dataset ",i0," from the file ",A,".")')
     &    dataset, trim(hWindFiles(cn)%file_name)
         call allMessage(ERROR,scratchMessage)
      endif
      write(scratchMessage,'(A,i0,A)')
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
      real(8), intent(out) :: wvnx(:) ! wind u-velocity
      real(8), intent(out) :: wvny(:) ! wind v-velocity
      integer, intent(in)   :: cn      ! cycle number, or hwind file number
      real(8), intent(in)  :: rv      ! ramp value from hwind file
      real(8), intent(in)  :: tw      ! time weighting for time interpolation
C
      real(8)              :: dh      ! readability shorthand for hwind dh
      integer               :: il      ! readability shorthand for hwind integerLimit
      real(8), pointer     :: hwu(:,:)! readability shorthand for hwind u
      real(8), pointer     :: hwv(:,:)! readability shorthand for hwind v
      real(8)              :: oneOverW2 ! 1/w**2 ; const in interp weights
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
       use sizes, only : globaldir 
       use global, only : ihot, nws, nabout, openFileForRead,
     &                setMessageSource,unsetMessageSource,allMessage,
     &                logMessage,scratchMessage,DEBUG,ECHO,INFO,
     &                ERROR,WARNING
      use constants, only: PRBCKGRND,RHOWAT0
      use mesh, only : np, x, y, slam0, sfea0
      implicit none

      real(8), intent(in) :: timeloc ! adcirc time in seconds since coldstart
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
       USE MESH, ONLY: ICS, DRVSPCOORSROTS, CYLINDERMAP
       USE CONSTANTS, ONLY: DEG2RAD,PRBCKGRND,G,RHOWAT0
       USE GLOBAL, ONLY: IFSPROTS,
     &        setMessageSource,unsetMessageSource,allMessage,
     &        logMessage,scratchMessage,DEBUG,ECHO,INFO,ERROR
       USE KDTREE2_MODULE     ! USED FOR FAST SEARCHING

       IMPLICIT NONE
       REAL(8), ALLOCATABLE :: XY_GFDL(:,:),XLOC(:)
       REAL(8), ALLOCATABLE :: PRESS_GFDL(:)
       REAL(8), ALLOCATABLE :: U_WIND_GFDL(:),V_WIND_GFDL(:)
       REAL(8), ALLOCATABLE :: WTS(:)
       INTEGER, ALLOCATABLE :: LOCS(:),II2(:)

       REAL(8) :: LON_GFDL,LAT_GFDL,TOL,DIST, LonR, LatR
       REAL(8) :: RPTS_GFDL,UTMP,VTMP,TTMP,MTMP,SPTMP,PTMP
       REAL(8) :: LONGTMP,LATTMP,HTMP,NTMP,XSTA,YSTA,MAXDIST
       REAL(8), INTENT(IN) :: SLAM0,SFEA0
       REAL(8), INTENT(IN) :: MAX_INTERP_DIST
       REAL(8), INTENT(OUT) :: ADC_U_WIND(NADCPTS),ADC_V_WIND(NADCPTS)
       REAL(8), INTENT(OUT) :: ADC_PRESS(NADCPTS)
       REAL(8), INTENT(IN)  :: ADC_X(NADCPTS),ADC_Y(NADCPTS)
!       REAL(8), PARAMETER ::  PI=3.141592653589793D0
!       REAL(8), PARAMETER ::  DEG2RAD = PI/180.D0  ! DEGREES TO RADIANS

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
        rpts_gfdl = real(npts_gfdl,8)

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
          IF ( IFSPROTS.eq. 1) THEN
             CALL DRVSPCOORSROTS( lonr, latr, LON_GFDL, LAT_GFDL )
          ELSE
             latr = LAT_GFDL ;  lonr = LON_GFDL; 
          END IF
          CALL CYLINDERMAP( XY_GFDL(1,I),XY_GFDL(2,I),
     &           lonr, latr, SLAM0, SFEA0, ICS) ; 
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

       write(scratchMessage,'("Maximum GFDL Interpolation Distance = ",e24.16)') maxdist
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
     &        WARNING
      use constants, only: PRBCKGRND,RHOWAT0,G

      implicit none
      integer :: i,numFiles,currentCycle
      logical :: cycleTimeFound
      real(8), intent(in) :: timeloc   ! adcirc time, seconds since coldstart
      real(8), intent(out) :: wvnx(:)  ! wind u-velocity (m/s)
      real(8), intent(out) :: wvny(:)  ! wind v-velocity (m/s)
      real(8), intent(out) :: press(:)  !barometric pressure (mH20)

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


C----------------------------------------------------------------------
C                S U B R O U T I N E   R S G E T
C***********************************************************************
C   Read onto the ADCIRC grid radiation stress fields in the PBL-JAG   *
C   (hurricane) model format.                                          *
C                           R.L.05/12/99                               *
C     jgf52.19: Added error checking and reporting as well as
C     execution tracing.
C***********************************************************************
      subroutine rsget(rsnx,rsny)
      use global, only : nabout
      use mesh, only : np
      implicit none
      real(8), intent(out), dimension(np) :: rsnx
      real(8), intent(out), dimension(np) :: rsny
      integer :: nhg
      integer, save :: lineNum = 0
      real(8) :: rsnxRaw, rsnyRaw
      character*80 pbljagf
      integer :: errorio
C
      call setMessageSource("rsget")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      rsnx(1:np) = 0.d0
      rsny(1:np) = 0.d0
C
      lineNum = lineNum + 1
      do while(.true.)
         read(23,'(a80)',end=500,err=600,iostat=errorio) pbljagf
         if (pbljagf(1:1).eq.'#') then
            call allMessage(ERROR,'Found "#" in column 1 of fort.23.')
            call windTerminate()
         endif
         if(pbljagf(2:2).ne.'#') exit
      enddo
      lineNum = lineNum + 1
      do while(.true.)
         read(pbljagf,'(i8,5e13.5)',end=500,err=600,iostat=errorio)
     &         nhg, rsnxRaw, rsnyRaw
         rsnx(nhg) = rsnxRaw
         rsny(nhg) = rsnyRaw
         lineNum = lineNum + 1
         read(23,'(a80)') pbljagf
         if (pbljagf(1:1).eq.'#') then
            call allMessage(ERROR,'Found "#" in column 1 of fort.23; '//
     &         ' this character should only appear in column 2.')
            call windTerminate()
         endif
         if(pbljagf(2:2).eq.'#') exit
      enddo

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      !cleans up the setMessageSource set at the beginning of the subroutine
      call unsetMessageSource()
      return
C
C     We jump to this section if there was an error reading a file.
500   call allMessage(ERROR,'Unexpectedly reached end-of-file.') ! END jumps here
600   call allMessage(ERROR,'I/O error during file access.') !  ERR jumps here
      write(scratchMessage,'(a,i0,a)')'Failed to read line ',lineNum,'.'
      call allMessage(ERROR,scratchMessage)
      call allMessage(ERROR,
     &                  'Check the fort.16 file for more information.')
      if (nabout.gt.0) then
         call allMessage(INFO,'Reducing the value of NABOUT to 0'
     &  //' will maximize the information written to the fort.16 file,'
     &  //' which may aid in troubleshooting this issue.')
      endif
      write(scratchMessage,'(a,i0,a)')
     &       'The value of the i/o error flag was ',errorio,'.'
      call allMessage(ERROR,scratchMessage)
      CALL windTerminate()
C----------------------------------------------------------------------
      end subroutine rsget
C----------------------------------------------------------------------


C----------------------------------------------------------------------
C          S U B R O U T I N E
C  I N I T   D Y N A M I C  W A T E R L E V E L  C O R R E C T I O N
C----------------------------------------------------------------------
C      jgf: Initialize i/o for water level correction file.
C
C      file format:
C      ! 1st line is a comment line
C      14400.0       ! 2nd line: time increment between datasets (s)
C      0.0           ! 3rd line: default water level correction value (m)
C      1     0.3048  ! subsequent lines consist of node number and correction value
C      3     0.3048  ! only for non-default values
C      10    0.30
C      115   0.3048
C      307   0.305
C      ##            ! "#" as 2nd character on the line:start of next dataset
C      1     0.2048  ! subsequent lines consist of node number and correction value
C      3     0.2048  ! in second data set ...
C      10    0.20
C      115   0.4048
C      307   0.405
C      ##
C      ...
C
C      If the file consists of only one data set, the data are treated
C      as a constant throughout the run. Time varying ramping of the
C      data are controlled by the dynamicWaterLevelCorrectionRampStart and
C      dynamicWaterLevelCorrectionRampEnd parameters in the fort.15
C      dynamicWaterLevelCorrectionControl namelist.
C
C----------------------------------------------------------------------
      subroutine initDynamicWaterLevelCorrection(timesec)
      use sizes, only : globaldir
      use global, only : nabout, np_g
      use mesh, only : np
      implicit none
      real(8), intent(in) :: timesec ! adcirc model time (s)
      integer :: errorio ! i/o error condition code

      call setMessageSource("initDynamicWaterLevelCorrection")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      !
      ! Check to make sure the correction file actually exists and if so, open it
      errorio = 0
      ! jgf: First try looking in the full domain input file directory
      call openFileForRead(dynamicWaterLevelCorrectionUn,
     &   trim(globaldir)//"/"//trim(dynamicWaterLevelCorrectionFileName), errorio)
      ! If the file was not found, in the directory with the other
      ! full domain input files, perhaps it was provided with a full path?
      if (errorio.gt.0) then
         call openFileForRead(dynamicWaterLevelCorrectionUn,
     &          trim(dynamicWaterLevelCorrectionFileName), errorio)
      endif
      ! give up if the file was still not found
      if (errorio.gt.0) then
         call windTerminate()
      endif
      ! file exists, allocate space to hold the data
#ifdef CMPI
      allocate(dynamicWaterLevelCorrectionFD(np_g))
      dynamicWaterLevelCorrectionFD(:) = 0.d0
#endif
      allocate(dynamicWaterLevelCorrectionsOld_m(np))
      allocate(dynamicWaterLevelCorrectionsNew_m(np))
      !
      ! read header
      dynamicWaterLevelCorrectionCurrentLine = 1
      ! first line is a comment line
      read(dynamicWaterLevelCorrectionUn,'(A)',end=256,err=258,iostat=errorio) dynamicWaterLevelCorrectionComment
      call logMessage(ECHO,"dynamicWaterLevelCorrectionComment: "//trim(dynamicWaterLevelCorrectionComment))
      dynamicWaterLevelCorrectionCurrentLine = dynamicWaterLevelCorrectionCurrentLine + 1
      ! 2nd line is time increment between datasets
      read(dynamicWaterLevelCorrectionUn,fmt=*,end=256,err=258,iostat=errorio) dynamicWaterLevelCorrectionTimeIncrement
      write(scratchMessage,'(a,f15.7)') "dynamicWaterLevelCorrectionTimeIncrement: ",
     &    dynamicWaterLevelCorrectionTimeIncrement
      call logMessage(ECHO,scratchMessage)
      dynamicWaterLevelCorrectionCurrentLine = dynamicWaterLevelCorrectionCurrentLine + 1
      ! tga20190628: this "if" is needed because it could lead to NaNs
      ! if dynamicWaterLevelCorrectionTimeIncrement is zero, plus no reason for negative.
      if (dynamicWaterLevelCorrectionTimeIncrement.le.0.d0) then
         call allMessage(ERROR,'Water level correction time increment ' //
     &       'dynamicWaterLevelCorrectionTimeIncrement must be greater than zero')
         call windTerminate()
      endif
      ! 3nd line is default correction value
      read(dynamicWaterLevelCorrectionUn,fmt=*,end=256,err=258,iostat=errorio) dynamicWaterLevelCorrectionDefaultValue
      write(scratchMessage,'(a,f15.7)') "dynamicWaterLevelCorrectionDefaultValue: ",
     &    dynamicWaterLevelCorrectionDefaultValue
      call logMessage(ECHO,scratchMessage)
      dynamicWaterLevelCorrectionCurrentLine = dynamicWaterLevelCorrectionCurrentLine + 1
      ! initialize array values to zero -- can be used as blank snaps
      dynamicWaterLevelCorrectionsOld_m(:) = 0.d0
      dynamicWaterLevelCorrectionsNew_m(:) = 0.d0
      !
      ! the common case is that dynamicWaterLevelCorrectionSkipSnaps will be 0 and we need
      ! to read the first two datasets
      !
      ! may also need to skip datasets (insert blank datasets if
      ! dynamicWaterLevelCorrectionSkipSnaps is negative)
      !
      ! if there is only one dataset, the correctionIsConstant flag is raised
      dynamicWaterLevelCorrectionDataSetCounter = 0
      dynamicWaterLevelCorrectionBlankSnaps = 2
      do while (dynamicWaterLevelCorrectionDataSetCounter.lt.(dynamicWaterLevelCorrectionSkipSnaps+2))
         if (dynamicWaterLevelCorrectionFileEnded.eqv..true.) then
            exit
         endif
         ! read the next snap from file
         call readDynamicWaterLevelCorrectionDataSet()
         dynamicWaterLevelCorrectionBlankSnaps = dynamicWaterLevelCorrectionBlankSnaps - 1
      end do
      !
      ! set the time corresponding to the two datasets we are
      ! interpolating between
      dynamicWaterLevelCorrectionTimeOld = timesec
      dynamicWaterLevelCorrectionTimeNew = timesec + dynamicWaterLevelCorrectionTimeIncrement
      !
      ! set ramp times relative to specified reference
      select case(trim(dynamicWaterLevelCorrectionRampReferenceTime))
      case('hotstart','Hotstart','HotStart','HOTSTART')
         dynamicWaterLevelCorrectionRampStart = timesec + dynamicWaterLevelCorrectionRampStart
         dynamicWaterLevelCorrectionRampEnd = timesec + dynamicWaterLevelCorrectionRampEnd
      case('coldstart','Coldstart','ColdStart','COLDSTART')
         if (dynamicWaterLevelCorrectionRampStart.lt.0.d0.or.dynamicWaterLevelCorrectionRampEnd.lt.0.d0) then
            call allMessage(ERROR,'Water level correction ramp reference '//
     &      'time was set to coldstart but negative ramp time(s) '//
     &      'were found.')
            call windTerminate()
         endif
      case default
         call allMessage(ERROR,'Water level correction ramp '//
     &   'reference time was set to "'//trim(dynamicWaterLevelCorrectionRampReferenceTime)//
     &   '" which is not a valid value. The dynamicWaterLevelCorrectionRampReferenceTime '//
     &   'parameter must be set to "coldstart" or "hotstart".')
         call windTerminate()
      end select
      !
      ! make sure start of ramp occurs before end of ramp
      if (dynamicWaterLevelCorrectionRampStart.gt.dynamicWaterLevelCorrectionRampEnd) then
         call allMessage(ERROR,
     &   'dynamicWaterLevelCorrectionRampStart is after dynamicWaterLevelCorrectionRampEnd.')
         call windTerminate()
      endif

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return

256   call allMessage(ERROR,'Unexpectedly reached end-of-file.') ! END jumps here
258   call allMessage(ERROR,'I/O error during file access.') !  ERR jumps here
      call allMessage(INFO,
     &   'Check the fort.16 file for more information.')
      if (nabout.gt.0) then
         call allMessage(INFO,
     &  'Reducing the value of NABOUT to 0'
     &  //' will maximize the information written to the fort.16 file,'
     &  //' which may aid in troubleshooting this issue.')
      endif
      write(scratchMessage,'("Could not read line ",i0," of ",a,".")')
     &  dynamicWaterLevelCorrectionCurrentLine, trim(dynamicWaterLevelCorrectionFileName)
      call allMessage(ERROR,scratchMessage)
      write(scratchMessage,'(a,i0,a)')
     &       'The value of the i/o error flag was ',errorio,'.'
      call allMessage(ERROR,scratchMessage)
      CALL windTerminate()


C----------------------------------------------------------------------
      end subroutine initdynamicWaterLevelCorrection
C----------------------------------------------------------------------

C----------------------------------------------------------------------
C      S U B R O U T I N E
C   G E T   D Y N A M I C  W A T E R L E V E L  C O R R E C T I O N S
C----------------------------------------------------------------------
C     jgf: Uses water level correction data that has already been
C     read from file, along with control parameters including
C     dynamicWaterLevelCorrectionMultiplier,
C     dynamicWaterLevelCorrectionRampStart, and
C     dynamicWaterLevelCorrectionRampEnd,
C     as well as the current adcirc time to set the correction values
C     for the current time step.
C
C     If all the data have been read from the file, the correction
C     value remains constant according to the last dataset read
C     from the file.
C----------------------------------------------------------------------
      subroutine getDynamicWaterLevelCorrections(dynamicWaterLevelCorrectionsArray, timesec)
      use mesh, only : np
      implicit none
      real(8), intent(in) :: timesec ! current adcirc time
      real(8), intent(out) :: dynamicWaterLevelCorrectionsArray(np)
      real(8) :: dynamicWaterLevelCorrectionRamp ! linear time varying ramp value
      integer :: i ! jgfdebug
      call setMessageSource("getDynamicWaterLevelCorrections")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      !
      ! hold data constant after last snap
      if ( dynamicWaterLevelCorrectionIsConstant.eqv..false.) then
         if ( (timesec.gt.dynamicWaterLevelCorrectionTimeNew)
     &        .and.(dynamicWaterLevelCorrectionFileEnded.eqv..true.) ) then
            dynamicWaterLevelCorrectionIsConstant = .true.
         endif
      endif
      !
      if ( dynamicWaterLevelCorrectionIsConstant.eqv..true. ) then
         dynamicWaterLevelCorrectionsArray(:) = dynamicWaterLevelCorrectionsNew_m(:)
      else
         ! read new data if it is time
         if ( timesec.gt.dynamicWaterLevelCorrectionTimeNew ) then
            do while (timesec.gt.dynamicWaterLevelCorrectionTimeNew)
               ! insert blanks if needed
               if (dynamicWaterLevelCorrectionSkipSnaps.lt.(-1*dynamicWaterLevelCorrectionBlankSnaps)) then
                  ! both datasets still have all zero values
                  dynamicWaterLevelCorrectionBlankSnaps = dynamicWaterLevelCorrectionBlankSnaps + 1
               else
                  call readDynamicWaterLevelCorrectionDataSet()
               endif
               dynamicWaterLevelCorrectionTimeOld = dynamicWaterLevelCorrectionTimeNew
               dynamicWaterLevelCorrectionTimeNew = dynamicWaterLevelCorrectionTimeNew +
     &                  dynamicWaterLevelCorrectionTimeIncrement
               if ( dynamicWaterLevelCorrectionFileEnded.eqv..true.) then
                  exit
               endif
            end do
         endif
         ! perform time interpolation
         dynamicWaterLevelCorrectionsArray(:) =
     &       (timesec - dynamicWaterLevelCorrectionTimeOld) /
     &              (dynamicWaterLevelCorrectionTimeNew - dynamicWaterLevelCorrectionTimeOld )
     &       * ( dynamicWaterLevelCorrectionsNew_m(:) - dynamicWaterLevelCorrectionsOld_m(:) )
     &          + dynamicWaterLevelCorrectionsOld_m(:)
      endif
      !
      ! set to zero everywhere prior to the start of the ramping period
      if (timesec.lt.dynamicWaterLevelCorrectionRampStart) then
         dynamicWaterLevelCorrectionsArray(:) = 0.d0
      endif
      ! compute and apply linear ramp during ramping period
      if ((timesec.ge.dynamicWaterLevelCorrectionRampStart).and.
     &    (timesec.le.dynamicWaterLevelCorrectionRampEnd)) then
         ! compute linear ramp value
         dynamicWaterLevelCorrectionRamp = ( timesec - dynamicWaterLevelCorrectionRampStart)
     &               / ( dynamicWaterLevelCorrectionRampEnd - dynamicWaterLevelCorrectionRampStart )
         dynamicWaterLevelCorrectionsArray(:) = dynamicWaterLevelCorrectionsArray(:) *
     &      dynamicWaterLevelCorrectionRamp
      endif
      ! just leave the raw values alone after the end of ramp period
      !
      ! apply correction multiplier if any
      dynamicWaterLevelCorrectionsArray(:) = dynamicWaterLevelCorrectionMultiplier *
     &  dynamicWaterLevelCorrectionsArray(:)

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return

C----------------------------------------------------------------------
      end subroutine getdynamicWaterLevelCorrections
C----------------------------------------------------------------------


C----------------------------------------------------------------------
C      S U B R O U T I N E
C R E A D   D Y N A M I C  W A T E R L E V E L  C O R R E C T I O N
C                                                       D A T A S E T
C----------------------------------------------------------------------
C     jgf: Reads a water level correction dataset from file.
C----------------------------------------------------------------------
      subroutine readDynamicWaterLevelCorrectionDataSet()
      use global, only : nabout, np_g, nodes_lg
      use mesh, only : np
      implicit none
      integer :: numCorrectionValues
      integer :: nhg ! fulldomain node number
      real(8) :: correction_val ! water level correction value at that node
      integer :: fd_node_number ! full domain node number
      integer :: sd_node_number ! subdomain node number
      logical :: nextCorrectionDataSetFound ! true if a # is found in 2nd col
      integer :: errorio ! error message code

      call setMessageSource("readDynamicWaterLevelCorrectionDataSet")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      !
      numCorrectionValues = np
#ifdef CMPI
      numCorrectionValues = np_g
#endif
      dynamicWaterLevelCorrectionsOld_m(:) = dynamicWaterLevelCorrectionsNew_m(:)
      dynamicWaterLevelCorrectionsNew_m(:) = dynamicWaterLevelCorrectionDefaultValue
      nextCorrectionDataSetFound = .false.
      do ! loop until we run out of data or hit a "#" in column 2
         read(dynamicWaterLevelCorrectionUn,fmt='(a)',end=133,err=258,iostat=errorio)
     &      scratchMessage
         dynamicWaterLevelCorrectionCurrentLine = dynamicWaterLevelCorrectionCurrentLine + 1
         ! check to see if the next dataset is about to start
         if (scratchMessage(2:2).eq."#") then
            nextCorrectionDataSetFound = .true.
            exit ! break out of do loop
         endif
         ! read the node number and value from the line
         read(scratchMessage,fmt=*,end=256,err=258,iostat=errorio)
     &        nhg, correction_val
         ! check to make sure this is a valid node number
         if ( (nhg.lt.1).or.(nhg.gt.numCorrectionValues) ) then
            write(scratchMessage,'(a,i0,a)')
     &        'The node number on line ',dynamicWaterLevelCorrectionCurrentLine,
     &        ' of the water '//
     &        ' correction file is out of range.'
            call allMessage(ERROR,scratchMessage)
            call windTerminate()
         endif
         ! set the correction value as specified
#ifdef CMPI
         dynamicWaterLevelCorrectionFD(nhg) = correction_val
#else
         dynamicWaterLevelCorrectionsNew_m(nhg) = correction_val
#endif
      end do
 133  continue ! jump to here if we've run out of lines in the correction file
      dynamicWaterLevelCorrectionDataSetCounter = dynamicWaterLevelCorrectionDataSetCounter + 1
      ! write log messages if that was the last dataset
      if (nextCorrectionDataSetFound.eqv..false.) then
         dynamicWaterLevelCorrectionFileEnded = .true.
         call logMessage(INFO,
     &    'Finished reading the final water level correction dataset from "'
     &      //trim(dynamicWaterLevelCorrectionFileName)//'".')
         if ( dynamicWaterLevelCorrectionDataSetCounter.eq.1 ) then
            dynamicWaterLevelCorrectionIsConstant = .true.
            call logMessage(INFO,
     &      'Only one water level correction dataset was found. '//
     &      'The water level correction will therefore be constant in time.')
         endif
      endif

#ifdef CMPI
      ! map the full domain node number to the corresponding
      ! subdomain node number
      DO sd_node_number=1,np
         ! get the corresponding fulldomain node number
         fd_node_number = ABS(NODES_LG(sd_node_number))
         ! fill in the subdomain arrays with the corresponding
         ! fulldomain values
         dynamicWaterLevelCorrectionsNew_m(sd_node_number) = dynamicWaterLevelCorrectionFD(fd_node_number)
      END DO
#endif


#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      return

C
C     We jump to this section if there was an error reading a file.
256   call allMessage(ERROR,'Unexpectedly reached end-of-file.') ! END jumps here
258   call allMessage(ERROR,'I/O error during file access.') !  ERR jumps here
      call allMessage(INFO,
     &   'Check the fort.16 file for more information.')
      if (nabout.gt.0) then
         call allMessage(INFO,
     &  'Reducing the value of NABOUT to 0'
     &  //' will maximize the information written to the fort.16 file,'
     &  //' which may aid in troubleshooting this issue.')
      endif
      write(scratchMessage,'("Could not read line ",i0," of ",a,".")')
     &  dynamicWaterLevelCorrectionCurrentLine, trim(dynamicWaterLevelCorrectionFileName)
      call allMessage(ERROR,scratchMessage)
      write(scratchMessage,'(a,i0,a)')
     &       'The value of the i/o error flag was ',errorio,'.'
      call allMessage(ERROR,scratchMessage)
      CALL windTerminate()

C----------------------------------------------------------------------
      end subroutine readDynamicWaterLevelCorrectionDataSet
C----------------------------------------------------------------------

C----------------------------------------------------------------------
C...  jgf50.38.05: Subroutine to terminate the run cleanly.
C...
C----------------------------------------------------------------------
      SUBROUTINE windTerminate()
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      IMPLICIT NONE
      INTEGER,ALLOCATABLE :: dmy(:)
C
      call setMessageSource("windTerminate")
#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
C
      dmy(10) = 1
      call allMessage(ERROR,"ADCIRC terminating.")
      ! intentionally create seg fault
#ifdef CMPI
      call msg_fini()
#endif
      CALL EXIT(1)
C

#if defined(WIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE windTerminate
C----------------------------------------------------------------------
#ifdef DATETIME
C----------------------------------------------------------------------
      SUBROUTINE NWS14INIT()
C----------------------------------------------------------------------
C.... WJP Adding in subroutines for GRB2 reading below (NWS = 14)
C.... or netcdf (depending on what is present)
C....
C.... CPB 10/2023: reorganized this subroutine to make it more readable
C.... as well as to improve input for TACC
C----------------------------------------------------------------------
      implicit none
      
      integer :: iret, PvarInd
      logical :: Fexists
      
      ! Initialise the datetime
      CurDT = basedatetime + timedelta(minutes=floor(WTIME1/60d0))
      NT = 0; NTC = 0;

      ! check if we are using grib2 or netcdf input files:
      CALL NWS14CHECK_FILETYPE()

      ! set filenames appropriately
      CALL NWS14SET_FILENAMES()

      ! make grib2 .inv files (returns if using netcdf)
      CALL NWS14GRB2_MAKE_INV()

      ! read fort.22 to get netCDF variable info (returns if using
      ! grib2)
      CALL NWS14NC_READ_F22()

      ! Calculate the interpolant weights
      CALL NWS14_CALC_INTERP_WTS()

      call allMessage(ECHO,'Finished init of NWS14')
C----------------------------------------------------------------------
      END SUBROUTINE NWS14INIT
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE NWS14CHECK_FILETYPE()
C----------------------------------------------------------------------
C     CPB 10/2023: This subroutine simply checks whether we are using
C     grib2 or netCDF winds in NWS14 and sets the grb2flag logical
C     appropriately. It does this by:
C        1. Checking if the pressure grib2 file (fort.221.grb2) exists
C        2. If it does not exist or if we are compiled without grib2
C           then it checks if the pressure netCDF file (fort.221.nc)
C           exists.
C        3. If neither exist we abort with an error message.
C----------------------------------------------------------------------
      implicit none
      logical :: Fexists

      Pfile = 'fort.221.grb2'
      INQUIRE(file=Pfile,exist=Fexists)
#ifdef GRIB2API
      IF (Fexists) THEN

         grb2flag = .TRUE.
         RETURN
      ENDIF
#endif
      Pfile = 'fort.221.nc'
      INQUIRE(file=Pfile,exist=Fexists)
      IF (Fexists) THEN
         grb2flag = .FALSE.
         RETURN
      ELSE
         call allMessage(ERROR,
     &       'Neither .grb2 nor .nc wind files exist. Or, if .grb2 '//
     &       'files exist, did you compile with GRIB2 compiler flags?')
         call windTerminate()
      ENDIF

C----------------------------------------------------------------------
      END SUBROUTINE NWS14CHECK_FILETYPE
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE NWS14SET_FILENAMES()
C----------------------------------------------------------------------
C     CPB 10/2023: This subroutine sets the met forcing filenames based
C     on grb2flag as well as the existence or lack thereof of various
C     files. It also opens the netcdf files in read only mode to
C     eliminate opening/closing the files repeatedly
C----------------------------------------------------------------------
#ifdef ADCNETCDF
      use netcdf
#endif
#ifdef CMPI
      USE MESSENGER, ONLY: MSG_BARRIER
#endif
      USE SIZES, ONLY : MYPROC
      USE GLOBAL, ONLY : NCICE
      implicit none
      logical :: Fexists
      ! grib2
      IF (grb2flag) THEN
         Pfile = 'fort.221.grb2'
         Pinv  = 'fort.221.inv'
         Wfile = 'fort.222.grb2'
         Wfile1 = 'fort.222.grb2'
         Winv  = 'fort.222.inv'
         Cfile = 'fort.225.grb2'
         Cinv  = 'fort.225.inv'

         Pvar = 'S'  ! searching for PRMSL or MSLET or PRES:SURFACE
                     ! prev: :PRMSL:mean sea level:'
         Uvar = ':UGRD:10 m above ground:'
         Vvar = ':VGRD:10 m above ground:'
         Cvar = ':ICEC:surface:'
         ! Assign dummy values to ncid variables
         PfileNCID = -9999
         WfileNCID = -9999
         Wfile1NCID = -9999
         IF (NCICE.EQ.14) THEN
            CfileNCID = -9999
         ENDIF
         RETURN
      ELSEIF (.NOT.grb2flag) THEN
         Pfile = 'fort.221.nc'
         Wfile = 'fort.222.nc'
         Wfile1 = 'fort.222.nc'
         Cfile = 'fort.225.nc'
         ! check if we are using two different wind files
         INQUIRE(FILE=Wfile,exist=Fexists)
         IF (.NOT.Fexists) THEN
            ! using 2 wind files
            Wfile = 'fort.222u.nc'
            Wfile1 = 'fort.222v.nc'
         ENDIF
         ! create dummy names for the inventory files
         Pinv = 'dmy'
         Winv = 'dmy'
         Cinv = 'dmy'
         ! create dummy names for the grib2 variables
         Pvar = 'dmy'
         Uvar = 'dmy'
         Vvar = 'dmy'
         Cvar = 'dmy'
#ifdef ADCNETCDF
#ifdef CMPI
         IF (MYPROC.EQ.0) THEN
#endif
            call Check_err(NF90_OPEN(Pfile,nf90_nowrite,PfileNCID))
            call Check_err(NF90_OPEN(Wfile,nf90_nowrite,WfileNCID))
            call Check_err(NF90_OPEN(Wfile1,nf90_nowrite,Wfile1NCID))
            IF (NCICE.EQ.14) THEN
               call Check_err(NF90_OPEN(Cfile,nf90_nowrite,CfileNCID))
            ENDIF
#ifdef CMPI
         ENDIF
         CALL MSG_BARRIER()
#endif
#endif
         RETURN
      ELSE
         ! should be unreachable
         CALL windTerminate()
      ENDIF
C----------------------------------------------------------------------
      END SUBROUTINE NWS14SET_FILENAMES
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE NWS14GRB2_MAKE_INV()
C----------------------------------------------------------------------
C     CPB 10/2023: This subroutine makes the .inv files for grib2
C     meteorological forcing files. If we are not using grib2 then it
C     just returns
C----------------------------------------------------------------------
      USE SIZES, ONLY : MYPROC
#ifdef GRIB2API
      USE wgrib2api
#endif
      USE GLOBAL, ONLY : NCICE
#ifdef CMPI
      USE MESSENGER, ONLY: MSG_BARRIER
#endif
      implicit none
      logical :: Fexists
      integer :: iret
      ! return if we are using netcdf
      IF (.NOT.grb2flag) THEN
         RETURN
      ENDIF
#ifdef GRIB2API
#ifdef CMPI
      ! only do this on proc 0
      IF (MYPROC.EQ.0) THEN
#endif
         ! pressure file (fort.221.grb2)
         inquire(file=Pinv,exist=Fexists)
         if (.not.Fexists) then
            iret = grb2_mk_inv(Pfile, Pinv)
            if (iret.ne.0) then
               call allMessage(ERROR,
     &                        'Fatal error in reading fort.221.grb2.')
               call windTerminate()
            endif
            call allMessage(ECHO,'successfully read fort.221.grb2 '
     &                         //'and wrote out fort.221.inv file')
         endif
         ! wind file (fort.222.grb2)
         inquire(file=Winv,exist=Fexists)
         if (.not.Fexists) then
            iret = grb2_mk_inv(Wfile, Winv)
            if (iret.ne.0) then
               call allMessage(ERROR,
     &                         'Fatal error in reading fort.222.grb2.')
               call windTerminate()
            endif
            call allMessage(ECHO,'successfully read fort.222.grb2 '
     &                         //'and wrote out fort.222.inv file')
         endif
         ! ice file (if we are using it) (fort.225.grb2)
         if (NCICE.eq.14) then
            inquire(file=Cinv,exist=Fexists)
            if (.not.Fexists) then
               iret = grb2_mk_inv(Cfile, Cinv)
               if (iret.ne.0) then
                  call allMessage(ERROR,
     &                         'Fatal error in reading fort.225.grb2.')
                  call windTerminate()
               endif
               call allMessage(ECHO,'successfully read fort.225.grb2 '
     &                         //'and wrote out fort.225.inv file')
            endif
         endif
#ifdef CMPI
      ENDIF
      ! wait until .inv files are all written before proceeding
      CALL MSG_BARRIER()
#endif
#endif
C----------------------------------------------------------------------
      END SUBROUTINE NWS14GRB2_MAKE_INV
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE NWS14NC_READ_F22()
C----------------------------------------------------------------------
C     CPB 10/2023: This subroutine reads the fort.22 file for NWS14 to
C     get netcdf variable names as well as starting indices for reading
C     in the met data. If we are using grib2 it just returns
C----------------------------------------------------------------------
      USE SIZES, ONLY : MYPROC
      USE GLOBAL, ONLY : NCICE
      implicit none
      logical :: Fexists
      integer :: ierr, PvarInd
      ! return if we are using grib2
      IF (grb2flag) THEN
         RETURN
      ENDIF
      ! open the fort.22
      OPEN(22,file='fort.22',ACTION='READ',IOSTAT=ierr)
      ! if we had trouble opening then abort
      IF (ierr.NE.0) THEN
         call allMessage(ERROR,'Unable to open fort.22.')
         call windTerminate()
      ENDIF
      read(22,*) Tdim    ! Time dimension name
      read(22,*) Tvar    ! Time variable name
      read(22,*) Tformat ! Format of time datestr
      read(22,*) Londim  ! Lon dimension name
      read(22,*) Lonvar  ! Lon var name
      read(22,*) Latdim  ! Lat dimension name
      read(22,*) Latvar  ! Lat var name
      read(22,*) Pvar    ! Pressure var name
      read(22,*) Uvar    ! U10 var name
      read(22,*) Vvar    ! V10 Var name
      if (NCICE.eq.14) then
         read(22,*) Cvar ! Ice var name
      endif
      close(22)
      PvarInd = index(Pvar,'HPa')
      if (PvarInd > 0) then
         rhoWat0g = rhoWat0g/100d0
         Pvar = Pvar(1:PvarInd-1)
         write(16,*) 'Pressure is in Hpa, New Pvar = ',trim(Pvar)
      endif
      ! get starting time index for netcdf file
      CALL READNWS14_NC_StaTime(Pfile, PfileNCID)
      RETURN
C----------------------------------------------------------------------
      END SUBROUTINE NWS14NC_READ_F22
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE NWS14_CALC_INTERP_WTS()
C----------------------------------------------------------------------
C     CPB 10/2023: This subroutine calculates the interpolation indices
C     (indp) and weights (weightsp) to interpolate from the grib2 or
C     netcdf input meterological forcing to the ADCIRC grid. NWS = 14
C     was originally set up to use NOAA CFSv2/GFS met forcing products.
C     These meteorological products provide two different resolutions
C     for wind forcing but not for pressure (or ice). To limit
C     pre-processing necessary to use these products, this subroutine
C     checks to see if the wind and pressure have the same dimensions.
C     If so, it only allocates one set of interpolant indices (indp) and
C     weights (weightsp). If they are different then it stores two
C     different sets. If necessary, it also checks the dimensions of the
C     ice forcing and sets the interpolants for ice to be equal to
C     either the wind or pressure forcing depending on the grid.
C----------------------------------------------------------------------
      use mesh, only: NP, SLAM, SFEA, bl_interp, bl_interp2
      use kdtree2_module, only: kdtree2, kdtree2_create
      USE GLOBAL, ONLY : NCICE
      implicit none
      integer :: i, ii, nxp, nyp, indt(4), NTkeep, xi, yi, cc
      real,allocatable,dimension(:,:) :: Pmsl, U10, V10, Cice, lat, lon
      real(8), allocatable ::  lonv(:), latv(:), XY(:,:)
      real(8) :: xx, yy, wt(4) 
      logical :: regular_grid
      type(kdtree2), pointer :: tree
      ! determine if wind/pressure/ice grids match or not.
      ! assume grids are the same to start
      ub = 1
      ubc = 1 ! index where we store the ice data (will match either
              ! wind or pressure). Assume match pressure first.
      ! get pressure dimensions
      CALL READNWS14LatLon(Pfile,lat,lon,Pinv,Pvar,PfileNCID)
      nxp = ubound(lon,1); nyp = ubound(lat,2)
      ! get wind dimensions and see if they match pressure
      CALL READNWS14LatLon(Wfile,lat,lon,Winv,Uvar,WfileNCID)
      if (nxp.NE.ubound(lon,1).OR.nyp.NE.ubound(lat,2)) then
         ! if we don't match both lat and lon then we need to store two
         ! interpolation indices/weights
         ub = 2
         nxp = ubound(lon,1)
         nyp = ubound(lat,2)
      endif
      ! If winds and pressure do not match, then check ice file to see
      ! whether we should use wind or pressure interpolants. If wind and
      ! pressure DO match then we assume ice also matches.
      IF (NCICE.EQ.14.AND.ub.EQ.2) THEN
         CALL READNWS14LatLon(Cfile,lat,lon,Cinv,Cvar,CfileNCID)
         IF (nxp.EQ.ubound(lon,1).AND.nyp.EQ.ubound(lat,2)) THEN
            ! match winds
            ubc = 2
         ENDIF
      ENDIF
      ! Allocate the indices and weights
      allocate(indp(ub,4,np),weightsp(ub,4,np))
      ! Make the interpolant weights. If w 
      do ii = 1,ub
         if (ii.eq.1) then
            CALL READNWS14LatLon(Pfile,lat,lon,Pinv,Pvar,PfileNCID)
         else
            CALL READNWS14LatLon(Wfile,lat,lon,Winv,Uvar,WfileNCID)
         endif
         nxp = ubound(lon,1); nyp = ubound(lat,2)
         ! Check if structured regular..
         if (abs(lon(1,1) - lon(1,nyp)) < 1d-6) then
            regular_grid = .true. 
            allocate(lonv(nxp),latv(nyp))
            lonv = lon(:,1)
            latv = lat(1,:)
         else       
            regular_grid = .false. 
            ! Convert point matrices to 2 X (NXP*NYP) point vector
            allocate(XY(2,nxp*nyp))
            cc = 0
            do xi = 1,nxp
               do yi = 1,nyp
                  cc = cc + 1
                  XY(1,cc) = lon(xi,yi)
                  XY(2,cc) = lat(xi,yi)
               enddo
            enddo
            ! Create the search tree
            tree => kdtree2_create(XY) 
         endif
         ! Loop over each node and interpolate
         do i = 1,np
            xx = rad2deg*slam(i)
            ! Convert our numbers if grids are 0 to 360
            if (maxval(lon).gt.180d0.and.xx < 0d0) xx = xx + 360d0
            yy = rad2deg*sfea(i)
            if (regular_grid) then
               call bl_interp(nxp,lonv,nyp,latv,xx,yy,indt,wt)
            else
               call bl_interp2(nxp,nyp,lon,lat,xx,yy,indt,wt,tree)
            endif
            indp(ii,:,i) = indt
            weightsp(ii,:,i) = wt
         enddo
         if (regular_grid) then
            deallocate(lon,lat,latv,lonv)
         else
            deallocate(lon,lat,XY)
         endif
      enddo
C----------------------------------------------------------------------
      END SUBROUTINE NWS14_CALC_INTERP_WTS
C------------------------------------------------------1----------------
C----------------------------------------------------------------------
      SUBROUTINE NWS14GET(WVNX2,WVNY2,PRN2,CICE2)
      use mesh, only: NP, SLAM, SFEA, bl_interp, bl_interp2
      implicit none

      real(8),intent(out) :: PRN2(NP), WVNX2(NP), WVNY2(NP)
      real(8),intent(out),optional :: CICE2(NP)
      integer :: i, ii, nxp, nyp, indt(4), NTkeep, xi, yi, cc
      real,allocatable,dimension(:,:) :: Pmsl, U10, V10, Cice, lat, lon
      real(8), allocatable ::  lonv(:), latv(:), XY(:,:)
      real(8) :: xx, yy, wt(4) 
      logical :: regular_grid

      ! Show the date in str_date format
      call logMessage(ECHO,'CurDT strng: '//
     &                CurDT%strftime("%Y-%m-%d %H:%M"))
      
      ! Get the Pressure
      !':PRMSL:mean sea level:'
      call READNWS14(Pfile,Pinv,Pvar,Pmsl)

      ! Get the U-winds 
      call READNWS14(Wfile,Winv,Uvar,U10)
      
      ! Get the V-winds 
      call READNWS14(Wfile1,Winv,Vvar,V10)

      ! Get the ice concentration if required 
      if (present(CICE2)) then
         NTkeep = NT; NT = NTC
         call READNWS14(Cfile,Cinv,Cvar,Cice)
         NT = NTkeep
      endif

      ! Doing the interpolation from their grid to ours
      do i = 1,np
         if (indp(1,1,i) < 0) then
            PRN2(i) = PRBCKGRND_MH2O
         else           
            PRN2(i) = ( Pmsl(indp(1,1,i),indp(1,2,i))*weightsp(1,1,i) + 
     &             Pmsl(indp(1,3,i),indp(1,2,i))*weightsp(1,2,i)      +
     &             Pmsl(indp(1,1,i),indp(1,4,i))*weightsp(1,3,i)      +
     &             Pmsl(indp(1,3,i),indp(1,4,i))*weightsp(1,4,i) )  
     &              / rhoWat0g
         endif
         if (indp(ub,1,i) < 0) then
            WVNX2(i) = 0d0; WVNY2(i) = 0d0
         else
            WVNX2(i) = U10(indp(ub,1,i),indp(ub,2,i))*weightsp(ub,1,i) +
     &              U10(indp(ub,3,i),indp(ub,2,i))*weightsp(ub,2,i)    +
     &              U10(indp(ub,1,i),indp(ub,4,i))*weightsp(ub,3,i)    +
     &              U10(indp(ub,3,i),indp(ub,4,i))*weightsp(ub,4,i)
            WVNY2(i) = V10(indp(ub,1,i),indp(ub,2,i))*weightsp(ub,1,i) +
     &              V10(indp(ub,3,i),indp(ub,2,i))*weightsp(ub,2,i)    +
     &              V10(indp(ub,1,i),indp(ub,4,i))*weightsp(ub,3,i)    +
     &              V10(indp(ub,3,i),indp(ub,4,i))*weightsp(ub,4,i)
         endif
#ifdef DEBUG_NWS14
         if (i.gt.1) cycle 
         write(16,*) 'NWS=14 Values check for i = 1'
         write(16,*) 'Pressure: ',PRN2(i)
         write(16,*) 'Winds, U10: ',WVNX2(i),'V10: ',WVNY2(i)
#endif
      enddo 
      deallocate(Pmsl,U10,V10)
      
      ! Add WTIMINC on CurDT for next WTIME  
      CurDT = CurDT + timedelta(minutes=nint(WTIMINC/60d0))
      ! Next time index for the netcdf reading
      NT = NT + 1
      
      ! If ice concentration is present
      if (present(CICE2)) then
         do i = 1,np
            if (indp(ubc,1,i) < 0) then
               CICE2(i) = 0d0;
            else
               CICE2(i) = 
     &             Cice(indp(ubc,1,i),indp(ubc,2,i))*weightsp(ubc,1,i) +
     &             Cice(indp(ubc,3,i),indp(ubc,2,i))*weightsp(ubc,2,i) +
     &             Cice(indp(ubc,1,i),indp(ubc,4,i))*weightsp(ubc,3,i) +
     &             Cice(indp(ubc,3,i),indp(ubc,4,i))*weightsp(ubc,4,i)
            endif
         enddo
         deallocate(Cice)
         ! Next ice time index for the netcdf reading
         NTC = NTC + 1
      endif
      !
C----------------------------------------------------------------------
      END SUBROUTINE NWS14GET
C----------------------------------------------------------------------
      SUBROUTINE READNWS14(fileN,invN,var,data2,data22)
#ifdef GRIB2API
      use wgrib2api
#endif
#ifdef ADCNETCDF
      use netcdf
#endif
      USE GL2LOC_MAPPING ,ONLY : BcastToLocal_2DRealArray, 
     &                           BcastToLocal_Int
      USE SIZES, ONLY : MYPROC
      implicit none
      integer :: iret, iter
      !integer,intent(in) :: startI(4), endI(4)
      character(len=200),intent(in) :: var, fileN, invN
      real,allocatable,intent(out) :: data2(:,:) 
      real,allocatable,intent(out),optional :: data22(:,:) 
      character(len=19),allocatable :: TimeStr(:) 
      character(len=200) :: str_date
      integer :: i, NC_ID, Temp_ID, Lat_ID, Lon_ID, NX, NY, TL, ndims
      ! CPB 4/20/2023
      REAL(8),ALLOCATABLE :: NCTIMES(:)
      ! Getting the date in str_date format
      if (grb2flag) then
         str_date = ':start_FT='//CurDT%strftime("%Y%m%d%H")//'0000:'
      endif      

      if (grb2flag) then
         ! GRIB2
         iret = -1; iter = 0
         do while (iret.le.0.and.iter < 10) 
            if (iter > 0) then
               call logMessage(WARNING,'Trying to read again '//var)
               call sleep(5)
            endif
#ifdef GRIB2API
            iret = grb2_inq(fileN,invN,var,str_date,data2=data2)
            if (iret.gt.1.and.iter.eq.0) then
               ! May have two entries bcause of forecast/DA overlap
               call logMessage(DEBUG,'> 1 msg, trying to read '//var)
               iret = grb2_inq(fileN,invN,var,':anl:',
     &                         str_date,data2=data2)
            endif
#endif
            iter = iter + 1
         enddo
         if (iret.gt.0) then
            call logMessage(DEBUG,'Successfully read '//var)
         else
            call ArnoldSchwarzenegger(iret,var,fileN)
         endif
C..... 
      else
         ! NETCDF
#ifdef ADCNETCDF
         IF (MYPROC.EQ.0) THEN
            call Check_err(NF90_OPEN(fileN,nf90_nowrite,NC_ID))
            call Check_err(NF90_INQ_DIMID(NC_ID,Latdim,Temp_ID))
            call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,len=NY))
            call Check_err(NF90_INQ_DIMID(NC_ID,Londim,Temp_ID))
            call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,len=NX))
         ENDIF
         CALL BcastToLocal_Int(NX)
         CALL BcastToLocal_Int(NY)
         allocate(data2(NX,NY))
         IF (MYPROC.EQ.0) THEN
            call Check_err(NF90_INQ_VARID(NC_ID,var,Temp_ID))
            call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,data2,
     &                     start=[1, 1, NT],count=[NX, NY, 1]))
            call Check_err(NF90_CLOSE(NC_ID))
         ENDIF
         CALL BCastToLocal_2DRealArray(Data2,NX,NY)
#endif
      endif
C----------------------------------------------------------------------
      END SUBROUTINE READNWS14
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE READNWS14LatLon(fileN,lat,lon,invN,var,ncid)
C----------------------------------------------------------------------
C     CPB 10/2023: Added to make NWS = 14 more more readable. Simply 
C     calls the appropriate subroutine based on the filetype we are 
C     using.
C----------------------------------------------------------------------
      implicit none
      CHARACTER(len=200),INTENT(IN) :: fileN
      REAL,ALLOCATABLE,INTENT(OUT) :: lat(:,:),lon(:,:)
      CHARACTER(LEN=200),INTENT(IN),OPTIONAL :: var, invN
      INTEGER,INTENT(IN),OPTIONAL :: ncid
      IF (grb2flag) THEN
         CALL READNWS14LatLon_grib2(fileN,invN,var,lat,lon)
      ELSEIF (.NOT.grb2flag) THEN
         CALL READNWS14LatLon_netCDF(fileN,ncid,lat,lon)
      ELSE
         ! should be unreachable
         CALL windterminate()
      ENDIF
C----------------------------------------------------------------------
      END SUBROUTINE READNWS14LatLon
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE READNWS14LatLon_grib2(fileN,invN,var,lat,lon)
C----------------------------------------------------------------------
C     CPB 10/2023: Reads in the lat and lon from a grib2 format
C     meteorological file. NOTE: reads in on Proc 0 and broadcasts
C----------------------------------------------------------------------
#ifdef GRIB2API
      use wgrib2api
#endif
      USE GL2LOC_MAPPING ,ONLY : BcastToLocal_2DRealArray, 
     &                           BcastToLocal_Int
      USE SIZES, ONLY : MYPROC
      implicit none
      integer :: iret, iter
      character(len=200),intent(in) :: var, fileN, invN
      real,allocatable,intent(out) :: lat(:,:), lon(:,: )
      integer :: NX, NY
      character(len=200) :: str_date
      CALL setMessageSource("READNWS14LatLon_grib2")
      ! Getting the date in str_date format
      str_date = ':start_FT='//CurDT%strftime("%Y%m%d%H")//'0000:'
      ! read in on core 0
      IF (MYPROC.EQ.0) THEN
         iret = -1; iter = 0
         do while (iret.le.0.and.iter < 10) 
            if (iter > 0) then
               call logMessage(WARNING,'Trying to read again '//var)
               call sleep(5)
            endif
#ifdef GRIB2API
            iret = grb2_inq(fileN,invN,var,str_date,
     &                      lat=lat,lon=lon)
#endif
            iter = iter + 1
         enddo
      ENDIF
      ! broadcast status of read
      CALL BcastToLocal_Int(iret)
      if (iret.gt.0) then
         call logMessage(ECHO,'Successfully read LatLon '//var)
      else
         call ArnoldSchwarzenegger(iret,var,fileN)
      endif
      ! broadcast data
      IF (MYPROC.EQ.0) THEN
         NX = UBOUND(lon,1)
         NY = UBOUND(lat,2)
      ENDIF
      CALL BcastToLocaL_Int(NX)
      Call BcastToLocal_Int(NY)
      IF (MYPROC.NE.0) THEN
         ALLOCATE( lat(NX,NY), lon(NX,NY) )
      ENDIF
      ! latitude
      CALL BcastToLocal_2DRealArray(lat,NX,NY)
      ! longitude
      CALL BcastToLocal_2DRealArray(lon,NX,NY)
      CALL unsetMessageSource()
C----------------------------------------------------------------------
      END SUBROUTINE READNWS14LatLon_grib2
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE READNWS14LatLon_netCDF(fileN,ncid,lat,lon)
C----------------------------------------------------------------------
C     CPB 10/2023: Reads in the lat and lon from a netCDF format
C     meteorological forcing file. NOTE: reads in on Proc 0 and
C     broadcasts.
C----------------------------------------------------------------------
#ifdef ADCNETCDF
      use netcdf
#endif
      USE GL2LOC_MAPPING ,ONLY : BcastToLocal_2DRealArray, 
     &                           BcastToLocal_Int
      USE SIZES, ONLY : MYPROC
      implicit none
      INTEGER,INTENT(IN) :: ncid ! already opened ncid tag
      character(len=200),intent(in) :: fileN
      real,allocatable,intent(out) :: lat(:,:), lon(:,:)
      integer :: i, Temp_ID, Lat_ID, Lon_ID, NX, NY, ndims
      CALL setMessageSource("READNWS14LatLon_netCDF")
#ifdef ADCNETCDF
      IF (MYPROC.EQ.0) THEN
         call Check_err(NF90_INQ_DIMID(ncid,Latdim,Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(ncid,Temp_ID,
     &                                         len=NY))
         call Check_err(NF90_INQ_DIMID(ncid,Londim,Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(ncid,Temp_ID,
     &                                         len=NX))
         allocate(lon(NX,NY),lat(NX,NY))
         call Check_err(NF90_INQ_VARID(ncid,Latvar,Lat_ID))
         call Check_err(nf90_inquire_variable(ncid,Lat_ID,
     &               ndims=ndims))
         call Check_err(NF90_INQ_VARID(ncid,Lonvar,Lon_ID))
         if (ndims.eq.3) then
            call Check_err(NF90_GET_VAR(ncid,Lat_ID,lat,
     &                     start=[1, 1, NT],count=[NX, NY, 1]))
            call Check_err(NF90_GET_VAR(ncid,Lon_ID,lon,
     &                     start=[1, 1, NT],count=[NX, NY, 1]))
         elseif (ndims.eq.2) then
            call Check_err(NF90_GET_VAR(ncid,Lat_ID,lat))
            call Check_err(NF90_GET_VAR(ncid,Lon_ID,lon))
         else
            call Check_err(NF90_GET_VAR(ncid,Lat_ID,lat(1,:)))
            call Check_err(NF90_GET_VAR(ncid,Lon_ID,lon(:,1)))
            do i = 2,NX
               lat(i,:) = lat(1,:)
            enddo
            do i = 2,NY
               lon(:,i) = lon(:,1)
            enddo
         endif
      ENDIF
      CALL BcastToLocal_Int(NX)
      CALL BcastToLocal_Int(NY)
      IF (MYPROC.NE.0) THEN
         ALLOCATE( lon(NX,NY), lat(NX,NY) )
      ENDIF
      CALL BcastToLocal_2DRealArray(lon,NX,NY)
      CALL BcastToLocal_2DRealArray(lat,NX,NY)
#endif
      CALL unsetMessageSource()
C----------------------------------------------------------------------
      END SUBROUTINE READNWS14LatLon_netCDF
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE READNWS14_NC_StaTime(FileN, NC_ID)
C----------------------------------------------------------------------
C       CPB 10/2023: Sets the time index from which we start reading in
C       netCDF meteorological forcing. NOTE: reads in on proc 0 and
C       broadcasts.
C----------------------------------------------------------------------
#ifdef ADCNETCDF
      use netcdf
#endif
      USE GL2LOC_MAPPING ,ONLY : BcastToLocal_Int
      USE SIZES, ONLY : MYPROC
      implicit none
      integer :: iret, iter
      INTEGER,INTENT(IN) :: NC_ID
      character(len=200),intent(in) :: fileN
      character(len=19),allocatable :: TimeStr(:) 
      character(len=200) :: str_date
      integer :: i, Temp_ID, TL
      ! CPB 4/20/2023
      REAL(8),ALLOCATABLE :: NCTIMES(:)
#ifdef ADCNETCDF
#ifdef CMPI
      IF (MYPROC.EQ.0) THEN
#endif
         call Check_err(NF90_INQ_DIMID(NC_ID,Tdim,Temp_ID))
         call Check_err(NF90_INQUIRE_DIMENSION(NC_ID,Temp_ID,
     &                                         len=TL))
         IF (Tformat(1:1).eq.'%') THEN
            ! in this case the netCDF files have a datetime variable
            ! so we can determine where to start based off of that
            CALL logMessage(ECHO,"fort.22 indicates that a "
     &          //"datetime variable is provided. Starting index "
     &          //"will be calculated.")
            ALLOCATE( TIMESTR(TL) )
            call Check_err(NF90_INQ_VARID(NC_ID,Tvar,Temp_ID))
            call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,TimeStr))
            ! Determine the first timestep index
            str_date = CurDT%strftime(trim(Tformat)) 
            do NT = 1,TL
               if (trim(str_date).eq.TimeStr(NT)) exit
            enddo
            call logMessage(ECHO,'Starting from '//TimeStr(NT))
            DEALLOCATE( TIMESTR )
         ELSEIF (Tvar.eq.'refdate') THEN
            ! in this case the reference date for the netCDF file is
            ! provided in the fort.22 and the units of the time
            ! variable are in the form of "seconds since refdate".
            ! This allows us to find the start index
            CALL logMessage(ECHO,"fort.22 provides a reference "
     &          //"date. Starting time index will be calculated.")
            ALLOCATE( NCTIMES(TL) )
            refdate = strptime(trim(Tformat),
     &                         "%Y-%m-%dT%H:%M:%S")
            CALL CHECK_ERR(NF90_INQ_VARID(NC_ID,TDIM,TEMP_ID))
            call Check_err(NF90_GET_VAR(NC_ID,Temp_ID,NCTIMES))
            DO NT = 1,TL
               stepdate = refdate +
     &                    timedelta(seconds=INT(NCTIMES(NT)))
               IF (stepdate.EQ.curDT) THEN
                  EXIT
               ENDIF
            ENDDO
            DEALLOCATE( NCTIMES )
         else
            ! If neither of those two things are true we just assume
            ! we start at the beginning of the file
            call logMessage(ECHO,
     &                      'Neither a datetime variable nor a '
     &          //'reference date were provided. Assume met forcing '
     &          //'starts at the beginning of the netCDF file.')
            NT = 1 
         endif
         NTC = NT
#ifdef CMPI
      endif
      CALL BcastToLocal_Int(NT)
      CALL BcastToLocal_Int(NTC)
#endif
#endif
      return
C----------------------------------------------------------------------
      END SUBROUTINE READNWS14_NC_StaTime
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE ArnoldSchwarzenegger(iret,var,fileN)
C----------------------------------------------------------------------
      IMPLICIT NONE
      integer,intent(in) :: iret
      character(len=200),intent(in) :: var, fileN
C
      if (iret.eq.0) then
         call allMessage(ERROR,'Cound not find message when reading '
     &                         //trim(var)//' in '//trim(fileN)//'.')
      elseif (iret.lt.0) then 
         call allMessage(ERROR,'Fatal error when reading '
     &                         //trim(var)//' in '//trim(fileN)//'.')
      elseif (iret.gt.1) then
         call allMessage(ERROR,'Found multiple messages when reading '
     &                         //trim(var)//' in '//trim(fileN)//'.')
      endif
      write(16,*) 'iret =',iret
      call windTerminate()
C----------------------------------------------------------------------
      END SUBROUTINE ArnoldSchwarzenegger
C-----------------------------------------------------------------------
#ifdef ADCNETCDF
C-----------------------------------------------------------------------
C     S U B R O U T I N E   C H E C K  _  E R R
C-----------------------------------------------------------------------
      subroutine check_err(iret)
      USE SIZES, ONLY : myproc
      USE GLOBAL, ONLY : screenUnit
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      USE NETCDF
      IMPLICIT NONE

      INTEGER, intent(in) :: iret
      integer,allocatable :: dmy

      call setMessageSource("check_err")
#ifdef ALL_TRACE
      call allMessage(DEBUG,"Enter.")
#endif
      if (iret .ne. NF90_NOERR) then
         call allMessage(ERROR,nf90_strerror(iret))
         call windTerminate()
      endif
#ifdef ALL_TRACE
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      end subroutine check_err
C-----------------------------------------------------------------------


#endif
C.....END WJP (NWS = 14) ..............................................
#endif

C-----------------------------------------------------------------------
C     S U B R O U T I N E   W I N D L I M I T E R 
C-----------------------------------------------------------------------

C.... DMW 202207 tail off wind forcing in very shallow water using a tanh,
C.... as a function of water column height.

      subroutine windLimiter(H,fwind)
      use global, only : h0
      
      IMPLICIT NONE
      
      REAL(8), intent(in) :: H
      REAL(8), intent(out) :: fwind
      
      call setMessageSource("windLimiter")
#ifdef ALL_TRACE
      call allMessage(DEBUG,"Enter.")
#endif
      fwind = 1.d0
      IF ((H .LE. 4.d0*H0) .AND. (H .GT. H0)) THEN
         fwind = .5d0*tanh(8.d0*(H-(2.5d0*H0))/(4.d0*H0))+.5d0
      ELSEIF (H .LE. H0) THEN
         fwind = 0.d0
      ENDIF

#ifdef ALL_TRACE
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
C-----------------------------------------------------------------------
      end subroutine WindLimiter
C-----------------------------------------------------------------------
C----------------------------------------------------------------------
      END MODULE
C----------------------------------------------------------------------
C----------------------------------------------------------------------
