!-----------------------------------------------------------------------
!  MODULE OWIWIND_NETCDF
!-----------------------------------------------------------------------
!> @author JC Detrich, NC State University
!> @author Alexander Crosby, Oceanweather Inc., alexc@oceanweather.com
!>
!> @copyright Dr. R.A. Luettich and Dr. J.J. Westerink
!>
!> @brief This module handles the I/O and interpolation to the mesh for Oceanweather
!> NetCDF format wind and pressure fields.
!>
!> The module is initialized by calling nws13init and subsequently calling
!> nws13get. The code will place data into the WVNX, WVNY, and PRN arrays, as
!> well as EyeLonR and EyeLatR for Powell drag if FoundEye is True.
!>
!> The fort.22.nc file is a NetCDF file containing groups, each with u/v wind
!> and pressure inputs. Can be named differently if specified in the namelist
!> in the fort.15.
!>
!> An example of the multi-group NetCDF stucture:
!>
!> @verbatim
!>   netcdf fort.22.nc {
!>
!>   // global attributes:
!>   		:group_order = "Main Landfall MovingStorm" ;
!>   		:institution = "Oceanweather Inc. (OWI)" ;
!>   		:conventions = "CF-1.6 OWI-NWS13" ;
!>
!>   group: Main {
!>     dimensions:
!>     	yi = 211 ;
!>     	xi = 221 ;
!>     	time = 1585 ;
!>     variables:
!>     	float lon(yi, xi) ;
!>     		lon:_FillValue = NaNf ;
!>     		lon:units = "degrees_east" ;
!>     		lon:standard_name = "longitude" ;
!>     		lon:axis = "X" ;
!>     		lon:coordinates = "time lat lon" ;
!>     	float lat(yi, xi) ;
!>     		lat:_FillValue = NaNf ;
!>     		lat:units = "degrees_north" ;
!>     		lat:standard_name = "latitude" ;
!>     		lat:axis = "Y" ;
!>     		lat:coordinates = "time lat lon" ;
!>     	int64 time(time) ;
!>     		time:units = "minutes since 1990-01-01T01:00:00" ;
!>     	float U10(time, yi, xi) ;
!>     		U10:_FillValue = NaNf ;
!>     		U10:units = "m s-1" ;
!>     		U10:coordinates = "time lat lon" ;
!>     	float V10(time, yi, xi) ;
!>     		V10:_FillValue = NaNf ;
!>     		V10:units = "m s-1" ;
!>     		V10:coordinates = "time lat lon" ;
!>     	float PSFC(time, yi, xi) ;
!>     		PSFC:_FillValue = NaNf ;
!>     		PSFC:units = "mb" ;
!>     		PSFC:coordinates = "time lat lon" ;
!>
!>     // group attributes:
!>     		:rank = 1 ;
!>     } // group Main
!>
!>   group: Landfall {
!>     dimensions:
!>     	yi = 151 ;
!>     	xi = 151 ;
!>     	time = 1585 ;
!>     variables:
!>     	float lon(yi, xi) ;
!>     		lon:_FillValue = NaNf ;
!>     		lon:units = "degrees_east" ;
!>     		lon:standard_name = "longitude" ;
!>     		lon:axis = "X" ;
!>     		lon:coordinates = "time lat lon" ;
!>     	float lat(yi, xi) ;
!>     		lat:_FillValue = NaNf ;
!>     		lat:units = "degrees_north" ;
!>     		lat:standard_name = "latitude" ;
!>     		lat:axis = "Y" ;
!>     		lat:coordinates = "time lat lon" ;
!>     	int64 time(time) ;
!>     		time:units = "minutes since 1990-01-01T01:00:00" ;
!>     	float U10(time, yi, xi) ;
!>     		U10:_FillValue = NaNf ;
!>     		U10:units = "m s-1" ;
!>     		U10:coordinates = "time lat lon" ;
!>     	float V10(time, yi, xi) ;
!>     		V10:_FillValue = NaNf ;
!>     		V10:units = "m s-1" ;
!>     		V10:coordinates = "time lat lon" ;
!>     	float PSFC(time, yi, xi) ;
!>     		PSFC:_FillValue = NaNf ;
!>     		PSFC:units = "mb" ;
!>     		PSFC:coordinates = "time lat lon" ;
!>
!>     // group attributes:
!>     		:rank = 2 ;
!>     } // group Landfall
!>
!>   group: MovingStorm {
!>     dimensions:
!>     	time = 1585 ;
!>     	yi = 501 ;
!>     	xi = 501 ;
!>     variables:
!>     	int64 time(time) ;
!>     		time:units = "minutes since 1990-01-01T01:00:00" ;
!>     	float lat(time, yi, xi) ;
!>     		lat:_FillValue = NaNf ;
!>     		lat:units = "degrees_north" ;
!>     		lat:standard_name = "latitude" ;
!>     		lat:axis = "Y" ;
!>     		lat:coordinates = "time lat lon" ;
!>     	float lon(time, yi, xi) ;
!>     		lon:_FillValue = NaNf ;
!>     		lon:units = "degrees_east" ;
!>     		lon:standard_name = "longitude" ;
!>     		lon:axis = "X" ;
!>     		lon:coordinates = "time lat lon" ;
!>     	float U10(time, yi, xi) ;
!>     		U10:_FillValue = NaNf ;
!>     		U10:units = "m s-1" ;
!>     		U10:coordinates = "time lat lon" ;
!>     	float V10(time, yi, xi) ;
!>     		V10:_FillValue = NaNf ;
!>     		V10:units = "m s-1" ;
!>     		V10:coordinates = "time lat lon" ;
!>     	float PSFC(time, yi, xi) ;
!>     		PSFC:_FillValue = NaNf ;
!>     		PSFC:units = "mb" ;
!>     		PSFC:coordinates = "time lat lon" ;
!>
!>     // group attributes:
!>     		:rank = 3 ;
!>     } // group MovingStorm
!>   }
!> @endverbatim
!>
!> Wind fields should be provided in groups ranked from coarse (lowest) to fine
!> (higher). The priority with which the wind fields are interpolated to the
!> mesh is determined rank.
!>
!> Wind fields are interpolated onto each node using a bilinear interpolation
!> scheme. If the mesh lies outside of all wind fields, it is set to
!> PRDEFLT (usually 1013mb) and 0.0 m/s velocity. WTIMEINC is used to specify
!> the timestep to perform input-grid(s) to mesh interpolation, for every other
!> timestep where getMeteorologicalForcing (wind.F) is called, temporal
!> interpolation is performed between on-mesh representations of the inputs.
!>
#ifdef ADCNETCDF
!-----------------------------------------------------------------------
      MODULE OWIWIND_NETCDF
!-----------------------------------------------------------------------
      USE SIZES, ONLY: SZ
      USE GLOBAL,ONLY : DEBUG, allMessage, setMessageSource, unsetMessageSource

      CHARACTER(LEN=15)  :: NWS13ColdStartString = "99999999.999999"
      CHARACTER(LEN=1000)  :: NWS13File = "fort.22.nc"
      CHARACTER(LEN=3)   :: NWS13GroupForPowell = "0"
      CHARACTER(LEN=10)  :: NWS13WindMultiplier = "1.0"
      CHARACTER(LEN=100) :: GroupNames(30)

      INTEGER :: ColdStartInMinutesSinceWindRefDateTime
      INTEGER :: NumGroup
      INTEGER :: NWS13GroupForPowellInt

      INTEGER :: NC_DIM
      INTEGER :: NC_ERR
      INTEGER :: NC_ID
      INTEGER :: NC_IDG
      INTEGER :: NC_VAR

      REAL(SZ) :: NWS13WindMultiplierReal

      TYPE NWS13Type
        INTEGER :: InclSnap
        INTEGER :: NextSnap
        INTEGER :: NumSnap
        INTEGER :: PrevSnap
        REAL(SZ) :: NextTime
        REAL(SZ) :: PrevTime
      ENDTYPE
      TYPE(NWS13Type),ALLOCATABLE :: NWS13(:)


      CONTAINS

!-----------------------------------------------------------------------
!> Initializes reading data from the Oceanweather (OWI) NetCDF wind/pre fields
!> @param[in]    timeloc model time
!> @param[inout] w       static weights for stationary "Main" grid
!-----------------------------------------------------------------------
      SUBROUTINE NWS13INIT(TimeLoc,W)
!-----------------------------------------------------------------------
      USE netcdf
      USE SIZES, ONLY: MYPROC
      USE MESH, ONLY: NP

      REAL(SZ),INTENT(IN) :: TimeLoc
      REAL(SZ),INTENT(INOUT),ALLOCATABLE :: W(:,:)

      CHARACTER(LEN=3100) :: GroupOrder
      CHARACTER(LEN=100) :: TimeUnits

      REAL(SZ),ALLOCATABLE :: Lat(:,:)
      REAL(SZ),ALLOCATABLE :: Lon(:,:)

      REAL(SZ) :: ColdStartDateTimeInJulianDays
      REAL(SZ) :: WindRefDateTimeInJulianDays

      INTEGER :: ILat(NP)
      INTEGER :: ILon(NP)
      INTEGER :: TempI(1)

      INTEGER :: GroupNameStart
      INTEGER :: GroupNameEnd
      INTEGER :: NextSpace
      INTEGER :: NumLat
      INTEGER :: NumLon

      INTEGER :: Day
      INTEGER :: Hour
      INTEGER :: Minute
      INTEGER :: Month
      INTEGER :: Second
      INTEGER :: Year

      call setMessageSource("nws13init")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      READ(NWS13WindMultiplier,*) NWS13WindMultiplierReal
      READ(NWS13GroupForPowell,*) NWS13GroupForPowellInt

      ! Open the NetCDF file.
      NC_ERR = NF90_OPEN(TRIM(ADJUSTL(NWS13File)),NF90_NOWRITE,NC_ID)

      ! Develop the number of groups, and the group names.
      NC_ERR = NF90_GET_ATT(NC_ID,NF90_GLOBAL,"group_order",GroupOrder)
      NumGroup = 0
      DO
        DO
          IF(INDEX(TRIM(GroupOrder)," ").EQ.1)THEN
            WRITE(UNIT=GroupOrder,FMT='(A)')
     &          GroupOrder(2:LEN_TRIM(GroupOrder))
          ELSE
            EXIT
          ENDIF
        ENDDO
        NextSpace = INDEX(TRIM(GroupOrder)," ")
        NumGroup = NumGroup + 1
        IF(NextSpace.GT.0)THEN
          GroupNameStart = 1
          GroupNameEnd = NextSpace - 1
        ElSE
          GroupNameStart = 1
          GroupNameEnd = LEN_TRIM(GroupOrder)
        ENDIF
        WRITE(UNIT=GroupNames(NumGroup),FMT='(A)')
     &      GroupOrder(GroupNameStart:GroupNameEnd)
        IF(NextSpace.LE.0)THEN
          EXIT
        ELSE
          WRITE(UNIT=GroupOrder,FMT='(A)')
     &      GroupOrder(NextSpace+1:LEN_TRIM(GroupOrder))
        ENDIF
      ENDDO

      IF(ALLOCATED(NWS13)) DEALLOCATE(NWS13)
      ALLOCATE(NWS13(1:NumGroup))

#ifdef HAVE_NETCDF4
      ! Connect to the first group.
      NC_ERR = NF90_INQ_NCID(NC_ID,TRIM(ADJUSTL(GroupNames(1))),NC_IDG)
#else
      NC_IDG = 1
#endif

      ! Find the starting date/time in YYYY-MM-DDTHH:MM:SS format,
      ! and convert it into Julian format.
      NC_ERR = NF90_INQ_VARID(NC_IDG,"time",NC_VAR)
      NC_ERR = NF90_GET_ATT(NC_IDG,NC_VAR,"units",TimeUnits)
      READ(UNIT=TimeUnits,FMT='(14X,I4,1X,I2,1X,I2,1X,I2,1X,I2,1X,I2)')
     &    Year,Month,Day,Hour,Minute,Second
      WindRefDateTimeInJulianDays = DBLE(JULIAN(Year,Month,Day))
     &    + FRACTIONTIME(Hour,Minute,Second)
      ! Convert the user-specified cold-start time into Julian format.
      READ(UNIT=NWS13ColdStartString,FMT='(I4,I2,I2,1X,I2,I2,I2)')
     &    Year,Month,Day,Hour,Minute,Second
      ColdStartDateTimeInJulianDays = DBLE(JULIAN(Year,Month,Day))
     &    + FRACTIONTIME(Hour,Minute,Second)
      ! Convert the cold-start time into minutes since the starting
      ! date/time specified in the NetCDF file.
      ColdStartInMinutesSinceWindRefDateTime
     &    = NINT(ColdStartDateTimeInJulianDays * 1440
     &    - WindRefDateTimeInJulianDays * 1440)

      ! Find the times associated with the first two time snaps in each
      ! group.  We will use these to build the first time snap.
      DO IG=1,NumGroup
        ! Connect to the group and its time variable.
#ifdef HAVE_NETCDF4
        NC_ERR = NF90_INQ_NCID(NC_ID,TRIM(ADJUSTL(GroupNames(IG))),
     &      NC_IDG)
#else
        NC_IDG = IG
#endif
        NC_ERR = NF90_INQ_DIMID(NC_IDG,"time",NC_DIM)
        NC_ERR = NF90_INQUIRE_DIMENSION(NC_IDG,NC_DIM,
     &      LEN=NWS13(IG)%NumSnap)
        NC_ERR = NF90_INQ_VARID(NC_IDG,"time",NC_VAR)
        ! Pull the time for the first snap ...
        NC_ERR = NF90_GET_VAR(NC_IDG,NC_VAR,TempI(:),
     &      START=(/1/),COUNT=(/1/))
        ! ... and load it into our array.
        NWS13(IG)%PrevTime = DBLE( TempI(1)
     &      - ColdStartInMinutesSinceWindRefDateTime ) * 60.D0
        NWS13(IG)%PrevSnap = 1
        ! Pull the time for the second snap ...
        NC_ERR = NF90_GET_VAR(NC_IDG,NC_VAR,TempI(:),
     &      START=(/2/),COUNT=(/1/))
        ! ... and also load it into our array.
        NWS13(IG)%NextTime = DBLE( TempI(1)
     &      - ColdStartInMinutesSinceWindRefDateTime ) * 60.D0
        NWS13(IG)%NextSnap = 2

        if (IG.EQ.1) then
          ! Read number of cells in longitude, latitude.
          NC_ERR = NF90_INQ_DIMID(NC_IDG,"xi",NC_DIM)
          NC_ERR = NF90_INQUIRE_DIMENSION(NC_IDG,NC_DIM,LEN=NumLon)
          NC_ERR = NF90_INQ_DIMID(NC_IDG,"yi",NC_DIM)
          NC_ERR = NF90_INQUIRE_DIMENSION(NC_IDG,NC_DIM,LEN=NumLat)

          ! Initialize arrays.
          IF(ALLOCATED(Lon)) DEALLOCATE(Lon)
          ALLOCATE(Lon(1:NumLon,1:NumLat))
          IF(ALLOCATED(Lat)) DEALLOCATE(Lat)
          ALLOCATE(Lat(1:NumLon,1:NumLat))
          IF(ALLOCATED(W)) DEALLOCATE(W)
          ALLOCATE(W(1:NP,1:6))

          ! Read mesh.
          NC_ERR = NF90_INQ_VARID(NC_IDG,"lon",NC_VAR)
          NC_ERR = NF90_GET_VAR(NC_IDG,NC_VAR,Lon(:,:),
     &        START=(/1,1,1/),COUNT=(/NumLon,NumLat,1/))
          NC_ERR = NF90_INQ_VARID(NC_IDG,"lat",NC_VAR)
          NC_ERR = NF90_GET_VAR(NC_IDG,NC_VAR,Lat(:,:),
     &        START=(/1,1,1/),COUNT=(/NumLon,NumLat,1/))
          ! get reusable weights for main group
          CALL NWS13INTERP(IG,NumLon,NumLat,Lon,Lat,ILon,ILat,W)
          W(:,5) = ILon
          W(:,6) = ILat
        endif

      ENDDO



      ! Close the NetCDF file.
      NC_ERR = NF90_CLOSE(NC_ID)

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.")
#endif
        call unsetMessageSource()

!-----------------------------------------------------------------------
      ENDSUBROUTINE NWS13INIT
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Reads multi-grid wind and pressure fields from the NetCDF file and
!>        interpolates them to the adcirc mesh
!> @param[inout] timeloc  model time
!> @param[inout] wvnx2    wind speed, x-direction
!> @param[inout] wvny2    wind speed, y-direction
!> @param[inout] prn2     atmospheric pressure
!> @param[inout] wtime2   time of next new field
!> @param[inout] eyelonr  storm center longitude
!> @param[inout] eyelatr  storm center latitude
!> @param[out]   foundeye boolean identifying if file provides storm center
!> @param[in]    w        static weights for stationary "Main" grid
!-----------------------------------------------------------------------
      SUBROUTINE NWS13GET(TimeLoc,WVNX2,WVNY2,PRN2,
     &    WTIME2,EyeLonR,EyeLatR,FoundEye,W)
!-----------------------------------------------------------------------
      USE GLOBAL, ONLY: G,
     &                  RHOWAT0,
     &                  RNDAY
      USE MESH, ONLY: NP
      USE netcdf
      USE SIZES, ONLY: MYPROC
      USE,INTRINSIC :: IEEE_ARITHMETIC

      IMPLICIT NONE

      REAL(SZ),INTENT(INOUT) :: EyeLatR(3)
      REAL(SZ),INTENT(INOUT) :: EyeLonR(3)
      LOGICAL ,INTENT(OUT)   :: FoundEye
      REAL(SZ),INTENT(OUT)   :: PRN2(NP)
      REAL(SZ),INTENT(IN)    :: TimeLoc
      REAL(SZ),INTENT(INOUT)   :: WVNX2(NP)
      REAL(SZ),INTENT(INOUT)   :: WVNY2(NP)
      REAL(SZ),INTENT(INOUT) :: WTIME2
      REAL(SZ),INTENT(IN) :: W(NP,6)

      CHARACTER(LEN=200) :: JunkC
      CHARACTER(LEN=200) :: Line

      INTEGER,SAVE :: PowellGroup = 0
      INTEGER,SAVE :: PowellGroupTemp = 0

      INTEGER :: CurrSnap
      INTEGER :: IG
      INTEGER :: ILat(NP)
      INTEGER :: ILon(NP)
      INTEGER :: IS
      INTEGER :: IV
      INTEGER :: ILT
      INTEGER :: ILN
      INTEGER :: NumLat
      INTEGER :: NumLon
      INTEGER :: TempI(1)

      REAL(SZ),ALLOCATABLE :: Lat(:,:)
      REAL(SZ),ALLOCATABLE :: Lon(:,:)
      REAL(SZ),ALLOCATABLE :: P(:,:)
      REAL(SZ),ALLOCATABLE :: U(:,:)
      REAL(SZ),ALLOCATABLE :: V(:,:)
      REAL(SZ),ALLOCATABLE :: Ws(:,:)

      REAL(SZ),ALLOCATABLE :: NextLat(:,:)
      REAL(SZ),ALLOCATABLE :: NextLon(:,:)
      REAL(SZ),ALLOCATABLE :: NextP(:,:)
      REAL(SZ),ALLOCATABLE :: NextU(:,:)
      REAL(SZ),ALLOCATABLE :: NextV(:,:)
      REAL(SZ),ALLOCATABLE :: NextWs(:,:)

      REAL(SZ),ALLOCATABLE :: PrevLat(:,:)
      REAL(SZ),ALLOCATABLE :: PrevLon(:,:)
      REAL(SZ),ALLOCATABLE :: PrevP(:,:)
      REAL(SZ),ALLOCATABLE :: PrevU(:,:)
      REAL(SZ),ALLOCATABLE :: PrevV(:,:)
      REAL(SZ),ALLOCATABLE :: PrevWs(:,:)

      REAL(SZ) :: CLat(1)
      REAL(SZ) :: CLon(1)
      REAL(SZ) :: NextCLat(1)
      REAL(SZ) :: NextCLon(1)
      REAL(SZ) :: PP(NP)
      REAL(SZ) :: PrevCLat(1)
      REAL(SZ) :: PrevCLon(1)
      REAL(SZ) :: TimeInterpFactor
      REAL(SZ) :: TimeTemp
      REAL(SZ) :: UU(NP)
      REAL(SZ) :: VV(NP)
      REAL(SZ) :: Wind(NP)
      REAL(SZ) :: sWdir
      REAL(SZ) :: W2(NP,4)
      REAL(SZ) :: var_scale_attr
      REAL(SZ) :: var_offset_attr
      REAL(SZ) :: D_QNAN, ZERO

      call setMessageSource("nws13get")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      ZERO = 0.D0
      D_QNAN = ZERO / ZERO
!      WRITE(6, *) D_QNAN

      ! First, fill the ADCIRC wind fields with background values.
      WVNX2 = 0.D0
      WVNY2 = 0.D0
      PRN2  = 101300.D0/RHOWAT0/G

      ! Open a connection to the wind file.
      NC_ERR = NF90_OPEN(TRIM(ADJUSTL(NWS13File)),NF90_NOWRITE,NC_ID)

      ! We need to determine how many groups are active at the current
      ! time in the simulation.  Loop over the groups and find which
      ! groups to include.
      NWS13(:)%InclSnap = 0
      !WTIME1 = WTIME2
      WTIME2 = -1.D0
      DO IG=1,NumGroup
        IF(TimeLoc.LT.NWS13(IG)%PrevTime)THEN
          ! Then we still haven't reached the start of this group,
          ! so we at least need to consider this as the end of the current
          ! interval.  If we don't find any other groups with a closer
          ! time, then this current interval will be zero winds.
          IF(WTIME2.LT.0.D0)THEN
            WTIME2 = NWS13(IG)%PrevTime
          ELSEIF(NWS13(IG)%PrevTime.LT.WTIME2)THEN
            WTIME2 = NWS13(IG)%PrevTime
          ENDIF
          CYCLE
        ENDIF
#ifdef HAVE_NETCDF4
        ! Find the current group and time variable in the file.
        NC_ERR = NF90_INQ_NCID(NC_ID,
     &      TRIM(ADJUSTL(GroupNames(IG))),NC_IDG)
#else
        NC_IDG = IG
#endif
        NC_ERR = NF90_INQ_VARID(NC_IDG,"time",NC_VAR)
        ! Check whether we have progressed past the end of the current
        ! wind snap, and if so, then shift to the next wind snap.
        DO WHILE(TimeLoc.GE.NWS13(IG)%NextTime)
          IF(NWS13(IG)%NextSnap.EQ.NWS13(IG)%NumSnap)THEN
            ! Then we have reached the end of the wind snaps for this
            ! group, and it will be excluded below.
            EXIT
          ENDIF
          NWS13(IG)%PrevSnap = NWS13(IG)%NextSnap
          NWS13(IG)%PrevTime = NWS13(IG)%NextTime
          NWS13(IG)%NextSnap = NWS13(IG)%NextSnap + 1
          ! Get the next time from the file.
          NC_ERR = NF90_GET_VAR(NC_IDG,NC_VAR,TempI(:),
     &        START=(/NWS13(IG)%NextSnap/),COUNT=(/1/))
          NWS13(IG)%NextTime = DBLE( TempI(1)
     &        - ColdStartInMinutesSinceWindRefDateTime ) * 60.D0
        ENDDO
        ! Then check whether we are within the current wind snap.
        IF( NWS13(IG)%PrevTime.LE.TimeLoc .AND.
     &      TimeLoc.LE.NWS13(IG)%NextTime )THEN
          NWS13(IG)%InclSnap = 1
          ! Update our end time if this is the first snap we found ...
          IF(WTIME2.LT.0.D0)THEN
            WTIME2 = NWS13(IG)%NextTime
          ! ... or if it is earlier than what we have found already.
          ELSEIF(NWS13(IG)%NextTime.LT.WTIME2)THEN
            WTIME2 = NWS13(IG)%NextTime
          ENDIF
        ENDIF
      ENDDO
      ! The only way we would not find a next snap is for the scenario
      ! when the simulation has progressed past the end of the wind
      ! file.  In that case, we can use the simulation length as the end
      ! of this next snap.
      IF(WTIME2.LT.0.D0)THEN
        WTIME2 = RNDAY*86400.D0
      ENDIF
      !WTIMINCX = WTIME2 - WTIME1

      ! If we don't need to include winds from any group (because the
      ! simulation has started before or extended after the information
      ! in the wind file), then we can return.
      IF(SUM(NWS13(:)%InclSnap).GT.0)THEN

        IF(MYPROC.EQ.0)THEN
          do ig=1,NumGroup
            IF(NWS13(IG)%InclSnap.EQ.0)THEN
              CYCLE
            ENDIF
            WRITE(Line,'(A,I1,A,A)')
     &          "Using group",ig," to ",
     &          "interpolate pressure and wind fields from time"
            WRITE(JunkC,'(F20.2)') nws13(ig)%PrevTime
            WRITE(Line,'(A,A,A,A)') TRIM(Line)," ",TRIM(ADJUSTL(JunkC)),
     &          " to time"
            WRITE(JunkC,'(F20.2)') nws13(ig)%NextTime
            WRITE(Line,'(A,A,A,A)') TRIM(Line)," ",TRIM(ADJUSTL(JunkC)),"."
            WRITE(*,'(A)') TRIM(Line)
          enddo
        ENDIF

        ! Increment the storm centers.
        EyeLonR(1) = EyeLonR(2)
        EyeLatR(1) = EyeLatR(2)
        EyeLonR(2) = EyeLonR(3)
        EyeLatR(2) = EyeLatR(3)
        EyeLonR(3) = 0.D0
        EyeLatR(3) = 0.D0
        FoundEye   = .FALSE.

        ! Loop over the groups.
        DO IG=1,NumGroup
          IF(NWS13(IG)%InclSnap.EQ.0)THEN
            CYCLE
          ENDIF
#ifdef   HAVE_NETCDF4
          ! Connect to this group.
          NC_ERR = NF90_INQ_NCID(NC_ID,TRIM(ADJUSTL(GroupNames(IG))),
     &        NC_IDG)
#else
          NC_IDG = IG
#endif  
          ! Read number of cells in longitude, latitude.
          NC_ERR = NF90_INQ_DIMID(NC_IDG,"xi",NC_DIM)
          NC_ERR = NF90_INQUIRE_DIMENSION(NC_IDG,NC_DIM,LEN=NumLon)
          NC_ERR = NF90_INQ_DIMID(NC_IDG,"yi",NC_DIM)
          NC_ERR = NF90_INQUIRE_DIMENSION(NC_IDG,NC_DIM,LEN=NumLat)

          ! Initialize arrays.

          ! These are temporary.
          IF(ALLOCATED(Lon)) DEALLOCATE(Lon)
          ALLOCATE(Lon(1:NumLon,1:NumLat))
          IF(ALLOCATED(Lat)) DEALLOCATE(Lat)
          ALLOCATE(Lat(1:NumLon,1:NumLat))
          IF(ALLOCATED(U)) DEALLOCATE(U)
          ALLOCATE(U(1:NumLon,1:NumLat))
          IF(ALLOCATED(V)) DEALLOCATE(V)
          ALLOCATE(V(1:NumLon,1:NumLat))
          IF(ALLOCATED(P)) DEALLOCATE(P)
          ALLOCATE(P(1:NumLon,1:NumLat))

          ! These have information from the previous snap.
          IF(ALLOCATED(PrevLon)) DEALLOCATE(PrevLon)
          ALLOCATE(PrevLon(1:NumLon,1:NumLat))
          IF(ALLOCATED(PrevLat)) DEALLOCATE(PrevLat)
          ALLOCATE(PrevLat(1:NumLon,1:NumLat))
          IF(ALLOCATED(PrevU)) DEALLOCATE(PrevU)
          ALLOCATE(PrevU(1:NumLon,1:NumLat))
          IF(ALLOCATED(PrevV)) DEALLOCATE(PrevV)
          ALLOCATE(PrevV(1:NumLon,1:NumLat))
          IF(ALLOCATED(PrevP)) DEALLOCATE(PrevP)
          ALLOCATE(PrevP(1:NumLon,1:NumLat))

          ! These have information from the next snap.
          IF(ALLOCATED(NextLon)) DEALLOCATE(NextLon)
          ALLOCATE(NextLon(1:NumLon,1:NumLat))
          IF(ALLOCATED(NextLat)) DEALLOCATE(NextLat)
          ALLOCATE(NextLat(1:NumLon,1:NumLat))
          IF(ALLOCATED(NextU)) DEALLOCATE(NextU)
          ALLOCATE(NextU(1:NumLon,1:NumLat))
          IF(ALLOCATED(NextV)) DEALLOCATE(NextV)
          ALLOCATE(NextV(1:NumLon,1:NumLat))
          IF(ALLOCATED(NextP)) DEALLOCATE(NextP)
          ALLOCATE(NextP(1:NumLon,1:NumLat))

          ! Loop over the snaps.
          DO IS=1,2
            IF(IS.EQ.1)THEN
              CurrSnap = NWS13(IG)%PrevSnap
            ELSEIF(IS.EQ.2)THEN
              CurrSnap = NWS13(IG)%NextSnap
            ENDIF

            ! 1. Read this snap from the wind file.

            ! Read mesh.
            NC_ERR = NF90_INQ_VARID(NC_IDG,"lon",NC_VAR)
            NC_ERR = NF90_GET_VAR(NC_IDG,NC_VAR,Lon(:,:),
     &          START=(/1,1,CurrSnap/),COUNT=(/NumLon,NumLat,1/))
            NC_ERR = NF90_INQ_VARID(NC_IDG,"lat",NC_VAR)
            NC_ERR = NF90_GET_VAR(NC_IDG,NC_VAR,Lat(:,:),
     &          START=(/1,1,CurrSnap/),COUNT=(/NumLon,NumLat,1/))

            ! Read winds and pressures.
            ! U-wind component
            NC_ERR = NF90_INQ_VARID(NC_IDG,"U10",NC_VAR)
            NC_ERR = NF90_GET_VAR(NC_IDG,NC_VAR,U(:,:),
     &          START=(/1,1,CurrSnap/),COUNT=(/NumLon,NumLat,1/))
            ! Apply scale/offset if attributes are present
            NC_ERR = NF90_GET_ATT(NC_IDG,NC_VAR, "_FillValue", var_scale_attr)
            IF (NC_ERR .EQ. 0) THEN
              DO ILN=1,NumLon
                DO ILT=1,NumLat
                  IF (U(ILN,ILT).EQ.var_scale_attr) THEN
                    U(ILN,ILT) = D_QNAN
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
            NC_ERR = NF90_GET_ATT(NC_IDG,NC_VAR, "scale_factor", var_scale_attr)
            IF (NC_ERR .EQ. 0) U = U * var_scale_attr
            NC_ERR = NF90_GET_ATT(NC_IDG,NC_VAR, "add_offset", var_offset_attr)
            IF (NC_ERR .EQ. 0) U = U + var_offset_attr
            ! V-wind component
            NC_ERR = NF90_INQ_VARID(NC_IDG,"V10",NC_VAR)
            NC_ERR = NF90_GET_VAR(NC_IDG,NC_VAR,V(:,:),
     &          START=(/1,1,CurrSnap/),COUNT=(/NumLon,NumLat,1/))
            ! Apply scale/offset if attributes are present
            NC_ERR = NF90_GET_ATT(NC_IDG,NC_VAR, "_FillValue", var_scale_attr)
            IF (NC_ERR .EQ. 0) THEN
              DO ILN=1,NumLon
                DO ILT=1,NumLat
                  IF (V(ILN,ILT).EQ.var_scale_attr) THEN
                    V(ILN,ILT) = D_QNAN
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
            NC_ERR = NF90_GET_ATT(NC_IDG,NC_VAR, "scale_factor", var_scale_attr)
            IF (NC_ERR .EQ. 0) V = V * var_scale_attr
            NC_ERR = NF90_GET_ATT(NC_IDG,NC_VAR, "add_offset", var_offset_attr)
            IF (NC_ERR .EQ. 0) V = V + var_offset_attr
            ! Pressure
            NC_ERR = NF90_INQ_VARID(NC_IDG,"PSFC",NC_VAR)
            NC_ERR = NF90_GET_VAR(NC_IDG,NC_VAR,P(:,:),
     &          START=(/1,1,CurrSnap/),COUNT=(/NumLon,NumLat,1/))
            ! Apply scale/offset if attributes are present
            NC_ERR = NF90_GET_ATT(NC_IDG,NC_VAR, "_FillValue", var_scale_attr)
            IF (NC_ERR .EQ. 0) THEN
              DO ILN=1,NumLon
                DO ILT=1,NumLat
                  IF (P(ILN,ILT).EQ.var_scale_attr) THEN
                    P(ILN,ILT) = D_QNAN
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
            NC_ERR = NF90_GET_ATT(NC_IDG,NC_VAR, "scale_factor", var_scale_attr)
            IF (NC_ERR .EQ. 0) P = P * var_scale_attr
            NC_ERR = NF90_GET_ATT(NC_IDG,NC_VAR, "add_offset", var_offset_attr)
            IF (NC_ERR .EQ. 0) P = P + var_offset_attr

            ! Adjust the wind velocity.
            U = U * NWS13WindMultiplierReal
            V = V * NWS13WindMultiplierReal

            ! Convert the pressures into meters of water.
            P = P*100.D0/RHOWAT0/G

            ! Try to read the storm center.
            NC_ERR = NF90_INQ_VARID(NC_IDG,"clon",NC_VAR)
            IF(NC_ERR.EQ.NF90_NOERR)THEN
              NC_ERR = NF90_GET_VAR(NC_IDG,NC_VAR,CLon(:),
     &            START=(/CurrSnap/),COUNT=(/1/))
              FoundEye = .TRUE.
            ENDIF
            NC_ERR = NF90_INQ_VARID(NC_IDG,"clat",NC_VAR)
            IF(NC_ERR.EQ.NF90_NOERR)THEN
              NC_ERR = NF90_GET_VAR(NC_IDG,NC_VAR,CLat(:),
     &            START=(/CurrSnap/),COUNT=(/1/))
              FoundEye = .TRUE.
            ENDIF

            ! Save the snaps.
            IF(IS.EQ.1)THEN
              PrevLon = Lon
              PrevLat = Lat
              PrevU = U
              PrevV = V
              PrevP = P
              IF(FoundEye)THEN
                PrevCLon = CLon
                PrevCLat = CLat
              ENDIF
            ELSEIF(IS.EQ.2)THEN
              NextLon = Lon
              NextLat = Lat
              NextU = U
              NextV = V
              NextP = P
              IF(FoundEye)THEN
                NextCLon = CLon
                NextCLat = CLat
              ENDIF
            ENDIF
          ENDDO

          ! 2. Interpolate the wind snaps to start/end of current interval.

          TimeInterpFactor = ( TimeLoc - NWS13(IG)%PrevTime )
     &             / ( NWS13(IG)%NextTime - NWS13(IG)%PrevTime )
!          IF(MYPROC.EQ.0)THEN
!            write(6,*) TimeInterpFactor,TimeLoc,NWS13(IG)%PrevTime,
!       &               NWS13(IG)%NextTime
!          ENDIF
          Lon  = PrevLon  + ( NextLon  - PrevLon  ) * TimeInterpFactor
          Lat  = PrevLat  + ( NextLat  - PrevLat  ) * TimeInterpFactor
          U    = PrevU    + ( NextU    - PrevU    ) * TimeInterpFactor
          V    = PrevV    + ( NextV    - PrevV    ) * TimeInterpFactor
          P    = PrevP    + ( NextP    - PrevP    ) * TimeInterpFactor
          ! Time Interpolate scalar wind
          PrevWs   = SQRT(PrevU**2+PrevV**2)
          NextWs   = SQRT(NextU**2+NextV**2)
          Ws    = PrevWs    + ( NextWs    - PrevWs    ) * TimeInterpFactor
          IF(FoundEye)THEN
              CLon = PrevCLon + ( NextCLon - PrevCLon ) * TimeInterpFactor
              CLat = PrevCLat + ( NextCLat - PrevCLat ) * TimeInterpFactor
          ENDIF

          ! 3. Interpolate the wind snaps onto the ADCIRC mesh.
          if (IG.GT.1) then
            CALL NWS13INTERP(IG,NumLon,NumLat,Lon,Lat,ILon,ILat,W2)
          else
            W2 = W(:, 1:4)
            ILon = W(:, 5)
            ILat = W(:, 6)
          endif
          DO IV=1,NP
            IF(ILon(IV).GT.0)THEN
              ! do not apply if we pick up a nan
              IF((Ws(ILon(IV)  ,ILat(IV)  )
     &          +  Ws(ILon(IV)+1,ILat(IV)  )
     &          +  Ws(ILon(IV)+1,ILat(IV)+1)
     &          +  Ws(ILon(IV)  ,ILat(IV)+1)).EQ.
     &            (Ws(ILon(IV)  ,ILat(IV)  )
     &          +  Ws(ILon(IV)+1,ILat(IV)  )
     &          +  Ws(ILon(IV)+1,ILat(IV)+1)
     &          +  Ws(ILon(IV)  ,ILat(IV)+1)))THEN
                 UU(IV) = W2(IV,1)*U(ILon(IV)  ,ILat(IV)  )
     &                  + W2(IV,2)*U(ILon(IV)+1,ILat(IV)  )
     &                  + W2(IV,3)*U(ILon(IV)+1,ILat(IV)+1)
     &                  + W2(IV,4)*U(ILon(IV)  ,ILat(IV)+1)
                 VV(IV) = W2(IV,1)*V(ILon(IV)  ,ILat(IV)  )
     &                  + W2(IV,2)*V(ILon(IV)+1,ILat(IV)  )
     &                  + W2(IV,3)*V(ILon(IV)+1,ILat(IV)+1)
     &                  + W2(IV,4)*V(ILon(IV)  ,ILat(IV)+1)
                 Wind(IV) = W2(IV,1)*Ws(ILon(IV)  ,ILat(IV)  )
     &                    + W2(IV,2)*Ws(ILon(IV)+1,ILat(IV)  )
     &                    + W2(IV,3)*Ws(ILon(IV)+1,ILat(IV)+1)
     &                    + W2(IV,4)*Ws(ILon(IV)  ,ILat(IV)+1)
              ENDIF
              ! do not apply if we pick up a nan
              IF((P(ILon(IV)  ,ILat(IV)  )
     &          + P(ILon(IV)+1,ILat(IV)  )
     &          + P(ILon(IV)+1,ILat(IV)+1)
     &          + P(ILon(IV)  ,ILat(IV)+1)).EQ.
     &           (P(ILon(IV)  ,ILat(IV)  )
     &          + P(ILon(IV)+1,ILat(IV)  )
     &          + P(ILon(IV)+1,ILat(IV)+1)
     &          + P(ILon(IV)  ,ILat(IV)+1)))THEN
                 PP(IV) = W2(IV,1)*P(ILon(IV)  ,ILat(IV)  )
     &                  + W2(IV,2)*P(ILon(IV)+1,ILat(IV)  )
     &                  + W2(IV,3)*P(ILon(IV)+1,ILat(IV)+1)
     &                  + W2(IV,4)*P(ILon(IV)  ,ILat(IV)+1)
              ENDIF
	      ! Adjust U/V Based on scalar wind magnitude
              IF ((UU(iV) .eq. 0.) .and.
     &            (VV(iV) .eq. 0.)) then
                sWdir = 0.
              ELSE
                sWdir = ATAN2(UU(iV),VV(iV))
                UU(IV) = -Wind(IV)*SIN(sWDir)
                VV(IV) = -Wind(IV)*COS(sWDir)
              ENDIF
              ! Casey 200528
              PRN2(IV) = PP(IV)
              WVNX2(IV) = UU(IV)
              WVNY2(IV) = VV(IV)
            ENDIF
          ENDDO

          ! Casey 200528
          !       PRN2 = PP
          !       WVNX2 = UU
          !       WVNY2 = VV

          IF(FoundEye)THEN
            IF(      ( NWS13GroupForPowellInt.EQ. 0 )
     &          .OR. ( NWS13GroupForPowellInt.EQ.IG ) )THEN
              ! Assign the storm center to the next spot in the array.
              EyeLonR(3) = CLon(1)
              EyeLatR(3) = CLat(1)
              PowellGroupTemp = IG
            ENDIF
          ENDIF

        ENDDO ! IG=1,NumGroup

      ENDIF

      ! Close the file.
      NC_ERR = NF90_CLOSE(NC_ID)

      IF(PowellGroupTemp.NE.PowellGroup)THEN
        ! If this is a different group then we have been using
        ! for the Powell wind drag scheme, then check whether
        ! the storm center in this new group is "close enough"
        ! to what we have been using.  If not, then reset the
        ! previous storm centers to force the scheme to use
        ! Garratt for the next few snaps.  It will eventually
        ! switch back to Powell for this new group / storm.
        IF( SQRT(   ( EyeLonR(3) - EyeLonR(2) )**2.D0
     &            + ( EyeLatR(3) - EyeLatR(2) )**2.D0 )
     &      .GT. 2.D0 )THEN
          EyeLonR(1) = 0.D0
          EyeLatR(1) = 0.D0
          EyeLonR(2) = 0.D0
          EyeLatR(2) = 0.D0
        ENDIF
        PowellGroup = PowellGroupTemp
      ENDIF

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.")
#endif
        call unsetMessageSource()

!-----------------------------------------------------------------------
      ENDSUBROUTINE NWS13GET
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!> @brief Raw grid to mesh node weight generation with kdtree
!> @param[in]    ig     netcdf group index
!> @param[in]    numlon x-dimension size
!> @param[in]    numlat y-dimension size
!> @param[in]    lon    grid longitudes (2D)
!> @param[in]    lat    grid latitudes (2D)
!> @param[out]   ilon   node positions in grid x-dimension
!> @param[in]    ilat   node positions in grid y-dimension
!> @param[inout] w      node interpolation weights
!-----------------------------------------------------------------------
      SUBROUTINE NWS13INTERP(IG,NumLon,NumLat,Lon,Lat,ILon,ILat,W)
!-----------------------------------------------------------------------
      USE GLOBAL, ONLY: RAD2DEG
      USE kdtree2_module
      USE MESH, ONLY: NP, SLAM, SFEA

      IMPLICIT NONE

      INTEGER,INTENT(IN)  :: IG
      INTEGER,INTENT(OUT) :: ILat(NP)
      INTEGER,INTENT(OUT) :: ILon(NP)
      INTEGER,INTENT(IN)  :: NumLat
      INTEGER,INTENT(IN)  :: NumLon

      REAL(SZ),INTENT(IN)  :: Lon(NumLon,NumLat)
      REAL(SZ),INTENT(IN)  :: Lat(NumLon,NumLat)
      REAL(SZ),INTENT(INOUT) :: W(NP,4)

      INTEGER :: Counter
      INTEGER :: ElemNumber
      INTEGER :: IA
      INTEGER :: IE
      INTEGER :: IO
      INTEGER :: IV
      INTEGER :: NumResult

      LOGICAL :: ElemFound

      REAL(SZ),ALLOCATABLE :: ElemCenter(:,:)
      REAL(SZ) :: AdcLat
      REAL(SZ) :: AdcLon
      REAL(SZ) :: MaxLat
      REAL(SZ) :: MaxLon
      REAL(SZ) :: MinLat
      REAL(SZ) :: MinLon
      REAL(SZ) :: W0
      REAL(SZ) :: W1
      REAL(SZ) :: W2
      REAL(SZ) :: W3
      REAL(SZ) :: W4

      TYPE(KDTREE2),POINTER :: kdTree
      TYPE(KDTREE2_RESULT),ALLOCATABLE :: kdResult(:)

      call setMessageSource("nws13interp")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      NumResult = (NumLon-1)*(NumLat-1)
      NumResult = MIN(NumResult,1000)

      IF(ALLOCATED(ElemCenter)) DEALLOCATE(ElemCenter)
      ALLOCATE(ElemCenter(1:2,1:(NumLon-1)*(NumLat-1)))
      Counter = 0
      DO IO=1,NumLon-1
        DO IA=1,NumLat-1
          Counter = Counter + 1
          ElemCenter(1,Counter) = 0.25D0*(
     &        Lon(IO  ,IA  )+Lon(IO+1,IA  )+
     &        Lon(IO+1,IA+1)+Lon(IO  ,IA+1))
          ElemCenter(2,Counter) = 0.25D0*(
     &        Lat(IO  ,IA  )+Lat(IO+1,IA  )+
     &        Lat(IO+1,IA+1)+Lat(IO  ,IA+1))
        ENDDO
      ENDDO
      kdTree => kdtree2_create(ElemCenter,rearrange=.TRUE.,
     &    sort=.TRUE.)
      IF(ALLOCATED(kdResult)) DEALLOCATE(kdResult)
      ALLOCATE(kdResult(1:NumResult))

      ILon = 0
      ILat = 0
      W = 0.D0

      MaxLat = MAXVAL(Lat)
      MaxLon = MAXVAL(Lon)
      MinLat = MINVAL(Lat)
      MinLon = MINVAL(Lon)
      DO IV=1,NP
        AdcLat = RAD2DEG*SFEA(IV)
        AdcLon = RAD2DEG*SLAM(IV)
        IF(AdcLon.LT.MinLon) CYCLE
        IF(AdcLon.GT.MaxLon) CYCLE
        IF(AdcLat.LT.MinLat) CYCLE
        IF(AdcLat.GT.MaxLat) CYCLE
        ElemFound = .FALSE.
        CALL kdtree2_n_nearest(tp=kdTree,
     &      qv=(/AdcLon,AdcLat/),
     &      nn=NumResult,results=kdResult)
        IE = 1
        DO WHILE((.NOT.ElemFound).AND.(IE.LE.NumResult))
          ElemNumber = kdResult(IE)%idx
          IF(MOD(ElemNumber,NumLat-1).EQ.0)THEN
            IO = ElemNumber/(NumLat-1)
            IA = NumLat-1
          ELSE
            IO = FLOOR(REAL(ElemNumber)/REAL(NumLat-1))+1
            IA = MOD(ElemNumber,NumLat-1)
          ENDIF
          W0 = (   Lon(IO  ,IA  ) * Lat(IO+1,IA  )
     &           - Lon(IO+1,IA  ) * Lat(IO  ,IA  ) )
     &       + (   Lon(IO+1,IA  ) * Lat(IO+1,IA+1)
     &           - Lon(IO+1,IA+1) * Lat(IO+1,IA  ) )
     &       + (   Lon(IO+1,IA+1) * Lat(IO  ,IA+1)
     &           - Lon(IO  ,IA+1) * Lat(IO+1,IA+1) )
     &       + (   Lon(IO  ,IA+1) * Lat(IO  ,IA  )
     &           - Lon(IO  ,IA  ) * Lat(IO  ,IA+1) )
          W0 = ABS(W0)
          W1 = (           AdcLon * Lat(IO+1,IA  )
     &           - Lon(IO+1,IA  ) * AdcLat         )
     &       + (   Lon(IO+1,IA  ) * Lat(IO+1,IA+1)
     &           - Lon(IO+1,IA+1) * Lat(IO+1,IA  ) )
     &       + (   Lon(IO+1,IA+1) * Lat(IO  ,IA+1)
     &           - Lon(IO  ,IA+1) * Lat(IO+1,IA+1) )
     &       + (   Lon(IO  ,IA+1) * AdcLat
     &           -         AdcLon * Lat(IO  ,IA+1) )
          W1 = ABS(W1)
          W2 = (   Lon(IO  ,IA  ) * AdcLat
     &           -         AdcLon * Lat(IO  ,IA  ) )
     &       + (           AdcLon * Lat(IO+1,IA+1)
     &           - Lon(IO+1,IA+1) * AdcLat         )
     &       + (   Lon(IO+1,IA+1) * Lat(IO  ,IA+1)
     &           - Lon(IO  ,IA+1) * Lat(IO+1,IA+1) )
     &       + (   Lon(IO  ,IA+1) * Lat(IO  ,IA  )
     &           - Lon(IO  ,IA  ) * Lat(IO  ,IA+1) )
          W2 = ABS(W2)
          W3 = (   Lon(IO  ,IA  ) * Lat(IO+1,IA  )
     &           - Lon(IO+1,IA  ) * Lat(IO  ,IA  ) )
     &       + (   Lon(IO+1,IA  ) * AdcLat
     &           -         AdcLon * Lat(IO+1,IA  ) )
     &       + (           AdcLon * Lat(IO  ,IA+1)
     &           - Lon(IO  ,IA+1) * AdcLat         )
     &       + (   Lon(IO  ,IA+1) * Lat(IO  ,IA  )
     &           - Lon(IO  ,IA  ) * Lat(IO  ,IA+1) )
          W3 = ABS(W3)
          W4 = (   Lon(IO  ,IA  ) * Lat(IO+1,IA  )
     &           - Lon(IO+1,IA  ) * Lat(IO  ,IA  ) )
     &       + (   Lon(IO+1,IA  ) * Lat(IO+1,IA+1)
     &           - Lon(IO+1,IA+1) * Lat(IO+1,IA  ) )
     &       + (   Lon(IO+1,IA+1) * AdcLat
     &           -         AdcLon * Lat(IO+1,IA+1) )
     &       + (           AdcLon * Lat(IO  ,IA  )
     &           - Lon(IO  ,IA  ) * AdcLat         )
          W4 = ABS(W4)
          IF((W1+W2+W3+W4).LE.(W0*2.001D0))THEN
            ElemFound = .TRUE.
            EXIT
          ENDIF
          IE = IE + 1
        ENDDO
        IF(ElemFound)THEN
          ILon(IV) = IO
          ILat(IV) = IA
          W(IV,1) = W1/(W0*2.D0)
          W(IV,2) = W2/(W0*2.D0)
          W(IV,3) = W3/(W0*2.D0)
          W(IV,4) = W4/(W0*2.D0)
        ENDIF

      ENDDO

      CALL kdtree2_destroy(tp=kdTree)

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
        call allMessage(DEBUG,"Return.")
#endif
        call unsetMessageSource()

!-----------------------------------------------------------------------
      ENDSUBROUTINE NWS13INTERP
!-----------------------------------------------------------------------


      INTEGER FUNCTION JULIAN(YEAR,MONTH,DAY)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: YEAR,MONTH,DAY
      call setMessageSource("julian")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      JULIAN = DAY-32075+1461*(YEAR+4800+(MONTH-14)/12)/4+367
     &    *(MONTH-2-(MONTH-14)/12*12)/12
     &    -3*((YEAR+4900+(MONTH-14)/12)/100)/4
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      END FUNCTION



      REAL(SZ) FUNCTION FRACTIONTIME(HOUR,MINUTE,SECOND)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: HOUR,MINUTE,SECOND
      call setMessageSource("fractiontime")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      FRACTIONTIME = (HOUR + (MINUTE + SECOND/60.D0)/60.D0)/24.D0
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
#endif
      call unsetMessageSource()
      END FUNCTION FRACTIONTIME

      ENDMODULE OWIWIND_NETCDF
#endif
