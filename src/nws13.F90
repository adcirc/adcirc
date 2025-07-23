!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2023 R.A. Luettich, Jr., J.J. Westerink
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------!
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
!>                   :group_order = "Main Landfall MovingStorm" ;
!>                   :institution = "Oceanweather Inc. (OWI)" ;
!>                   :conventions = "CF-1.6 OWI-NWS13" ;
!>
!>   group: Main {
!>     dimensions:
!>             yi = 211 ;
!>             xi = 221 ;
!>             time = 1585 ;
!>     variables:
!>             float lon(yi, xi) ;
!>                     lon:_FillValue = NaNf ;
!>                     lon:units = "degrees_east" ;
!>                     lon:standard_name = "longitude" ;
!>                     lon:axis = "X" ;
!>                     lon:coordinates = "time lat lon" ;
!>             float lat(yi, xi) ;
!>                     lat:_FillValue = NaNf ;
!>                     lat:units = "degrees_north" ;
!>                     lat:standard_name = "latitude" ;
!>                     lat:axis = "Y" ;
!>                     lat:coordinates = "time lat lon" ;
!>             int64 time(time) ;
!>                     time:units = "minutes since 1990-01-01T01:00:00" ;
!>             float U10(time, yi, xi) ;
!>                     U10:_FillValue = NaNf ;
!>                     U10:units = "m s-1" ;
!>                     U10:coordinates = "time lat lon" ;
!>             float V10(time, yi, xi) ;
!>                     V10:_FillValue = NaNf ;
!>                     V10:units = "m s-1" ;
!>                     V10:coordinates = "time lat lon" ;
!>             float PSFC(time, yi, xi) ;
!>                     PSFC:_FillValue = NaNf ;
!>                     PSFC:units = "mb" ;
!>                     PSFC:coordinates = "time lat lon" ;
!>
!>     // group attributes:
!>                     :rank = 1 ;
!>     } // group Main
!>
!>   group: Landfall {
!>     dimensions:
!>             yi = 151 ;
!>             xi = 151 ;
!>             time = 1585 ;
!>     variables:
!>             float lon(yi, xi) ;
!>                     lon:_FillValue = NaNf ;
!>                     lon:units = "degrees_east" ;
!>                     lon:standard_name = "longitude" ;
!>                     lon:axis = "X" ;
!>                     lon:coordinates = "time lat lon" ;
!>             float lat(yi, xi) ;
!>                     lat:_FillValue = NaNf ;
!>                     lat:units = "degrees_north" ;
!>                     lat:standard_name = "latitude" ;
!>                     lat:axis = "Y" ;
!>                     lat:coordinates = "time lat lon" ;
!>             int64 time(time) ;
!>                     time:units = "minutes since 1990-01-01T01:00:00" ;
!>             float U10(time, yi, xi) ;
!>                     U10:_FillValue = NaNf ;
!>                     U10:units = "m s-1" ;
!>                     U10:coordinates = "time lat lon" ;
!>             float V10(time, yi, xi) ;
!>                     V10:_FillValue = NaNf ;
!>                     V10:units = "m s-1" ;
!>                     V10:coordinates = "time lat lon" ;
!>             float PSFC(time, yi, xi) ;
!>                     PSFC:_FillValue = NaNf ;
!>                     PSFC:units = "mb" ;
!>                     PSFC:coordinates = "time lat lon" ;
!>
!>     // group attributes:
!>                     :rank = 2 ;
!>     } // group Landfall
!>
!>   group: MovingStorm {
!>     dimensions:
!>             time = 1585 ;
!>             yi = 501 ;
!>             xi = 501 ;
!>     variables:
!>             int64 time(time) ;
!>                     time:units = "minutes since 1990-01-01T01:00:00" ;
!>             float lat(time, yi, xi) ;
!>                     lat:_FillValue = NaNf ;
!>                     lat:units = "degrees_north" ;
!>                     lat:standard_name = "latitude" ;
!>                     lat:axis = "Y" ;
!>                     lat:coordinates = "time lat lon" ;
!>             float lon(time, yi, xi) ;
!>                     lon:_FillValue = NaNf ;
!>                     lon:units = "degrees_east" ;
!>                     lon:standard_name = "longitude" ;
!>                     lon:axis = "X" ;
!>                     lon:coordinates = "time lat lon" ;
!>             float U10(time, yi, xi) ;
!>                     U10:_FillValue = NaNf ;
!>                     U10:units = "m s-1" ;
!>                     U10:coordinates = "time lat lon" ;
!>             float V10(time, yi, xi) ;
!>                     V10:_FillValue = NaNf ;
!>                     V10:units = "m s-1" ;
!>                     V10:coordinates = "time lat lon" ;
!>             float PSFC(time, yi, xi) ;
!>                     PSFC:_FillValue = NaNf ;
!>                     PSFC:units = "mb" ;
!>                     PSFC:coordinates = "time lat lon" ;
!>
!>     // group attributes:
!>                     :rank = 3 ;
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
module mod_nws13
!-----------------------------------------------------------------------
   use GLOBAL, only: DEBUG, allMessage, setMessageSource, unsetMessageSource

   character(LEN=15)  :: NWS13ColdStartString = "99999999.999999"
   character(LEN=1000)  :: NWS13File = "fort.22.nc"
   character(LEN=3)   :: NWS13GroupForPowell = "0"
   character(LEN=10)  :: NWS13WindMultiplier = "1.0"
   character(LEN=100) :: GroupNames(30)

   integer :: ColdStartInMinutesSinceWindRefDateTime
   integer :: NumGroup
   integer :: NWS13GroupForPowellInt

   integer :: NC_DIM
   integer :: NC_ERR
   integer :: NC_ID
   integer :: NC_IDG
   integer :: NC_VAR

   integer :: PowellGroup = 0
   integer :: PowellGroupTemp = 0

   real(8) :: NWS13WindMultiplierReal

   type NWS13Type
      integer :: InclSnap
      integer :: NextSnap
      integer :: NumSnap
      integer :: PrevSnap
      real(8) :: NextTime
      real(8) :: PrevTime
   end type NWS13Type
   type(NWS13Type), allocatable :: NWS13(:)

   private

   public :: NWS13INIT, NWS13GET, NWS13ColdStartString, NWS13File, &
             NWS13GroupForPowell, NWS13WindMultiplier

contains

!-----------------------------------------------------------------------
!> Initializes reading data from the Oceanweather (OWI) NetCDF wind/pre fields
!> @param[in]    timeloc model time
!> @param[inout] w       static weights for stationary "Main" grid
!-----------------------------------------------------------------------
   subroutine NWS13INIT(W)
!-----------------------------------------------------------------------
      use netcdf, only: NF90_GET_ATT, NF90_OPEN, NF90_NOWRITE, NF90_INQ_NCID, &
                        NF90_INQUIRE_DIMENSION, NF90_INQ_VARID, NF90_GET_VAR, NF90_CLOSE, NF90_INQ_DIMID, &
                        NF90_GLOBAL
      use netcdf_error, only: check_err
      use MESH, only: NP

      implicit none

      real(8), intent(INOUT), allocatable :: W(:, :)

      character(LEN=3100) :: GroupOrder
      character(LEN=100) :: TimeUnits

      real(8), allocatable :: Lat(:, :)
      real(8), allocatable :: Lon(:, :)

      real(8) :: ColdStartDateTimeInJulianDays
      real(8) :: WindRefDateTimeInJulianDays

      integer :: ILat(NP)
      integer :: ILon(NP)
      integer :: TempI(1)

      integer :: GroupNameStart
      integer :: GroupNameEnd
      integer :: NextSpace
      integer :: NumLat
      integer :: NumLon

      integer :: Day
      integer :: Hour
      integer :: Minute
      integer :: Month
      integer :: Second
      integer :: Year
      integer :: IG

      call setMessageSource("nws13init")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      read (NWS13WindMultiplier, *) NWS13WindMultiplierReal
      read (NWS13GroupForPowell, *) NWS13GroupForPowellInt

! Open the NetCDF file.
      call check_err(NF90_OPEN(trim(adjustl(NWS13File)), NF90_NOWRITE, NC_ID))

! Develop the number of groups, and the group names.
      call check_err(NF90_GET_ATT(NC_ID, NF90_GLOBAL, "group_order", GroupOrder))
      NumGroup = 0
      do
         do
            if (index(trim(GroupOrder), " ") == 1) then
               GroupOrder = GroupOrder(2:len_trim(GroupOrder))
            else
               exit
            end if
         end do
         NextSpace = index(trim(GroupOrder), " ")
         NumGroup = NumGroup + 1
         if (NextSpace > 0) then
            GroupNameStart = 1
            GroupNameEnd = NextSpace - 1
         else
            GroupNameStart = 1
            GroupNameEnd = len_trim(GroupOrder)
         end if
         GroupNames(NumGroup) = GroupOrder(GroupNameStart:GroupNameEnd)
         if (NextSpace <= 0) then
            exit
         else
            GroupOrder = GroupOrder(NextSpace + 1:len_trim(GroupOrder))
         end if
      end do

      if (allocated(NWS13)) deallocate (NWS13)
      allocate (NWS13(1:NumGroup))

#ifdef HAVE_NETCDF4
! Connect to the first group.
      call check_err(NF90_INQ_NCID(NC_ID, trim(adjustl(GroupNames(1))), NC_IDG))
#else
      NC_IDG = 1
#endif

! Find the starting date/time in YYYY-MM-DDTHH:MM:SS format,
! and convert it into Julian format.
      call check_err(NF90_INQ_VARID(NC_IDG, "time", NC_VAR))
      call check_err(NF90_GET_ATT(NC_IDG, NC_VAR, "units", TimeUnits))
      read (UNIT=TimeUnits, FMT='(14X,I4,1X,I2,1X,I2,1X,I2,1X,I2,1X,I2)') &
         Year, Month, Day, Hour, Minute, Second
      WindRefDateTimeInJulianDays = dble(JULIAN(Year, Month, Day)) &
                                    + FRACTIONTIME(Hour, Minute, Second)
! Convert the user-specified cold-start time into Julian format.
      read (UNIT=NWS13ColdStartString, FMT='(I4,I2,I2,1X,I2,I2,I2)') &
         Year, Month, Day, Hour, Minute, Second
      ColdStartDateTimeInJulianDays = dble(JULIAN(Year, Month, Day)) &
                                      + FRACTIONTIME(Hour, Minute, Second)
! Convert the cold-start time into minutes since the starting
! date/time specified in the NetCDF file.
      ColdStartInMinutesSinceWindRefDateTime &
         = nint(ColdStartDateTimeInJulianDays*1440 &
                - WindRefDateTimeInJulianDays*1440)

! Find the times associated with the first two time snaps in each
! group.  We will use these to build the first time snap.
      do IG = 1, NumGroup
! Connect to the group and its time variable.
#ifdef HAVE_NETCDF4
         call check_err(NF90_INQ_NCID(NC_ID, trim(adjustl(GroupNames(IG))), &
                                      NC_IDG))
#else
         NC_IDG = IG
#endif
         call check_err(NF90_INQ_DIMID(NC_IDG, "time", NC_DIM))
         call check_err(NF90_INQUIRE_DIMENSION(NC_IDG, NC_DIM, &
                                               LEN=NWS13(IG)%NumSnap))
         call check_err(NF90_INQ_VARID(NC_IDG, "time", NC_VAR))
! Pull the time for the first snap ...
         call check_err(NF90_GET_VAR(NC_IDG, NC_VAR, TempI(:), &
                                     START=[1], COUNT=[1]))
! ... and load it into our array.
         NWS13(IG)%PrevTime = dble(TempI(1) &
                                   - ColdStartInMinutesSinceWindRefDateTime)*60.d0
         NWS13(IG)%PrevSnap = 1
! Pull the time for the second snap ...
         call check_err(NF90_GET_VAR(NC_IDG, NC_VAR, TempI(:), &
                                     START=[2], COUNT=[1]))
! ... and also load it into our array.
         NWS13(IG)%NextTime = dble(TempI(1) &
                                   - ColdStartInMinutesSinceWindRefDateTime)*60.d0
         NWS13(IG)%NextSnap = 2

         if (IG == 1) then
! Read number of cells in longitude, latitude.
            call check_err(NF90_INQ_DIMID(NC_IDG, "xi", NC_DIM))
            call check_err(NF90_INQUIRE_DIMENSION(NC_IDG, NC_DIM, LEN=NumLon))
            call check_err(NF90_INQ_DIMID(NC_IDG, "yi", NC_DIM))
            call check_err(NF90_INQUIRE_DIMENSION(NC_IDG, NC_DIM, LEN=NumLat))

! Initialize arrays.
            if (allocated(Lon)) deallocate (Lon)
            allocate (Lon(1:NumLon, 1:NumLat))
            if (allocated(Lat)) deallocate (Lat)
            allocate (Lat(1:NumLon, 1:NumLat))
            if (allocated(W)) deallocate (W)
            allocate (W(1:NP, 1:6))

! Read mesh.
            call check_err(NF90_INQ_VARID(NC_IDG, "lon", NC_VAR))
            call check_err(NF90_GET_VAR(NC_IDG, NC_VAR, Lon(:, :), &
                                        START=[1, 1, 1], COUNT=[NumLon, NumLat, 1]))
            call check_err(NF90_INQ_VARID(NC_IDG, "lat", NC_VAR))
            call check_err(NF90_GET_VAR(NC_IDG, NC_VAR, Lat(:, :), &
                                        START=[1, 1, 1], COUNT=[NumLon, NumLat, 1]))
! get reusable weights for main group
            call NWS13INTERP(NumLon, NumLat, Lon, Lat, ILon, ILat, W)
            W(:, 5) = dble(ILon)
            W(:, 6) = dble(ILat)
         end if

      end do

! Close the NetCDF file.
      call check_err(NF90_CLOSE(NC_ID))

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()

!-----------------------------------------------------------------------
   end subroutine NWS13INIT
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
   subroutine NWS13GET(TimeLoc, WVNX2, WVNY2, PRN2, &
                       WTIME2, EyeLonR, EyeLatR, FoundEye, W)
!-----------------------------------------------------------------------
      use ADC_CONSTANTS, only: G, RHOWAT0
      use GLOBAL, only: RNDAY
      use MESH, only: NP
      use netcdf, only: NF90_GET_VAR, NF90_INQ_VARID, NF90_INQ_NCID, &
                        NF90_INQUIRE_DIMENSION, NF90_GET_ATT, NF90_OPEN, NF90_NOWRITE, &
                        NF90_CLOSE, NF90_INQ_DIMID, NF90_NOERR
      use SIZES, only: MYPROC
      use netcdf_error, only: check_err

      implicit none

      real(8), parameter :: EPS = epsilon(1d0)

      real(8), intent(INOUT) :: EyeLatR(3)
      real(8), intent(INOUT) :: EyeLonR(3)
      logical, intent(OUT)   :: FoundEye
      real(8), intent(OUT)   :: PRN2(NP)
      real(8), intent(IN)    :: TimeLoc
      real(8), intent(INOUT) :: WVNX2(NP)
      real(8), intent(INOUT) :: WVNY2(NP)
      real(8), intent(OUT)   :: WTIME2
      real(8), intent(IN)    :: W(NP, 6)

      character(LEN=200) :: JunkC
      character(LEN=200) :: Line

      integer :: CurrSnap
      integer :: IG
      integer :: ILat(NP)
      integer :: ILon(NP)
      integer :: IS
      integer :: IV
      integer :: ILT
      integer :: ILN
      integer :: NumLat
      integer :: NumLon
      integer :: TempI(1)

      real(8), allocatable :: Lat(:, :)
      real(8), allocatable :: Lon(:, :)
      real(8), allocatable :: P(:, :)
      real(8), allocatable :: U(:, :)
      real(8), allocatable :: V(:, :)
      real(8), allocatable :: Ws(:, :)

      real(8), allocatable :: NextLat(:, :)
      real(8), allocatable :: NextLon(:, :)
      real(8), allocatable :: NextP(:, :)
      real(8), allocatable :: NextU(:, :)
      real(8), allocatable :: NextV(:, :)
      real(8), allocatable :: NextWs(:, :)

      real(8), allocatable :: PrevLat(:, :)
      real(8), allocatable :: PrevLon(:, :)
      real(8), allocatable :: PrevP(:, :)
      real(8), allocatable :: PrevU(:, :)
      real(8), allocatable :: PrevV(:, :)
      real(8), allocatable :: PrevWs(:, :)

      real(8) :: CLat(1)
      real(8) :: CLon(1)
      real(8) :: NextCLat(1)
      real(8) :: NextCLon(1)
      real(8) :: PP(NP)
      real(8) :: PrevCLat(1)
      real(8) :: PrevCLon(1)
      real(8) :: TimeInterpFactor
      real(8) :: UU(NP)
      real(8) :: VV(NP)
      real(8) :: Wind(NP)
      real(8) :: sWdir
      real(8) :: W2(NP, 4)
      real(8) :: var_scale_attr
      real(8) :: var_offset_attr
      real(8) :: D_QNAN, ZERO

      call setMessageSource("nws13get")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      ZERO = 0.d0
      D_QNAN = ZERO/ZERO
      PrevCLon(1) = 0d0
      PrevCLat(1) = 0d0
      NextCLon(1) = 0d0
      NextCLat(1) = 0d0

! First, fill the ADCIRC wind fields with background values.
      WVNX2 = 0.d0
      WVNY2 = 0.d0
      PRN2 = 101300.d0/RHOWAT0/G

! Open a connection to the wind file.
      call check_err(NF90_OPEN(trim(adjustl(NWS13File)), NF90_NOWRITE, NC_ID))

! We need to determine how many groups are active at the current
! time in the simulation.  Loop over the groups and find which
! groups to include.
      NWS13(:)%InclSnap = 0
!     !WTIME1 = WTIME2
      WTIME2 = -1.d0
      do IG = 1, NumGroup
         if (TimeLoc < NWS13(IG)%PrevTime) then
! Then we still haven't reached the start of this group,
! so we at least need to consider this as the end of the current
! interval.  If we don't find any other groups with a closer
! time, then this current interval will be zero winds.
            if (WTIME2 < 0.d0) then
               WTIME2 = NWS13(IG)%PrevTime
            elseif (NWS13(IG)%PrevTime < WTIME2) then
               WTIME2 = NWS13(IG)%PrevTime
            end if
            cycle
         end if
#ifdef HAVE_NETCDF4
! Find the current group and time variable in the file.
         call check_err(NF90_INQ_NCID(NC_ID, &
                                      trim(adjustl(GroupNames(IG))), NC_IDG))
#else
         NC_IDG = IG
#endif
         call check_err(NF90_INQ_VARID(NC_IDG, "time", NC_VAR))
! Check whether we have progressed past the end of the current
! wind snap, and if so, then shift to the next wind snap.
         do while (TimeLoc >= NWS13(IG)%NextTime)
            if (NWS13(IG)%NextSnap == NWS13(IG)%NumSnap) then
! Then we have reached the end of the wind snaps for this
! group, and it will be excluded below.
               exit
            end if
            NWS13(IG)%PrevSnap = NWS13(IG)%NextSnap
            NWS13(IG)%PrevTime = NWS13(IG)%NextTime
            NWS13(IG)%NextSnap = NWS13(IG)%NextSnap + 1
! Get the next time from the file.
            call check_err(NF90_GET_VAR(NC_IDG, NC_VAR, TempI(:), &
                                        START=[NWS13(IG)%NextSnap], COUNT=[1]))
            NWS13(IG)%NextTime = dble(TempI(1) &
                                      - ColdStartInMinutesSinceWindRefDateTime)*60.d0
         end do
! Then check whether we are within the current wind snap.
         if (NWS13(IG)%PrevTime <= TimeLoc .and. &
             TimeLoc <= NWS13(IG)%NextTime) then
            NWS13(IG)%InclSnap = 1
! Update our end time if this is the first snap we found ...
            if (WTIME2 < 0.d0) then
               WTIME2 = NWS13(IG)%NextTime
! ... or if it is earlier than what we have found already.
            elseif (NWS13(IG)%NextTime < WTIME2) then
               WTIME2 = NWS13(IG)%NextTime
            end if
         end if
      end do
! The only way we would not find a next snap is for the scenario
! when the simulation has progressed past the end of the wind
! file.  In that case, we can use the simulation length as the end
! of this next snap.
      if (WTIME2 < 0.d0) then
         WTIME2 = RNDAY*86400.d0
      end if
!     !WTIMINCX = WTIME2 - WTIME1

! If we don't need to include winds from any group (because the
! simulation has started before or extended after the information
! in the wind file), then we can return.
      if (sum(NWS13(:)%InclSnap) > 0) then

         if (MYPROC == 0) then
            do ig = 1, NumGroup
               if (NWS13(IG)%InclSnap == 0) then
                  cycle
               end if
               write (Line, '(A,I1,A,A)') &
                  "Using group", ig, " to ", &
                  "interpolate pressure and wind fields from time"
               write (JunkC, '(F20.2)') nws13(ig)%PrevTime
               write (Line, '(A,A,A,A)') trim(Line), " ", trim(adjustl(JunkC)), &
                  " to time"
               write (JunkC, '(F20.2)') nws13(ig)%NextTime
               write (Line, '(A,A,A,A)') trim(Line), " ", trim(adjustl(JunkC)), "."
               write (*, '(A)') trim(Line)
            end do
         end if

! Increment the storm centers.
         EyeLonR(1) = EyeLonR(2)
         EyeLatR(1) = EyeLatR(2)
         EyeLonR(2) = EyeLonR(3)
         EyeLatR(2) = EyeLatR(3)
         EyeLonR(3) = 0.d0
         EyeLatR(3) = 0.d0
         FoundEye = .false.

! Loop over the groups.
         do IG = 1, NumGroup
            if (NWS13(IG)%InclSnap == 0) then
               cycle
            end if
#ifdef   HAVE_NETCDF4
! Connect to this group.
            call check_err(NF90_INQ_NCID(NC_ID, trim(adjustl(GroupNames(IG))), &
                                         NC_IDG))
#else
            NC_IDG = IG
#endif
! Read number of cells in longitude, latitude.
            call check_err(NF90_INQ_DIMID(NC_IDG, "xi", NC_DIM))
            call check_err(NF90_INQUIRE_DIMENSION(NC_IDG, NC_DIM, LEN=NumLon))
            call check_err(NF90_INQ_DIMID(NC_IDG, "yi", NC_DIM))
            call check_err(NF90_INQUIRE_DIMENSION(NC_IDG, NC_DIM, LEN=NumLat))

! Initialize arrays.

! These are temporary.
            if (allocated(Lon)) deallocate (Lon)
            allocate (Lon(1:NumLon, 1:NumLat))
            if (allocated(Lat)) deallocate (Lat)
            allocate (Lat(1:NumLon, 1:NumLat))
            if (allocated(U)) deallocate (U)
            allocate (U(1:NumLon, 1:NumLat))
            if (allocated(V)) deallocate (V)
            allocate (V(1:NumLon, 1:NumLat))
            if (allocated(P)) deallocate (P)
            allocate (P(1:NumLon, 1:NumLat))

! These have information from the previous snap.
            if (allocated(PrevLon)) deallocate (PrevLon)
            allocate (PrevLon(1:NumLon, 1:NumLat))
            if (allocated(PrevLat)) deallocate (PrevLat)
            allocate (PrevLat(1:NumLon, 1:NumLat))
            if (allocated(PrevU)) deallocate (PrevU)
            allocate (PrevU(1:NumLon, 1:NumLat))
            if (allocated(PrevV)) deallocate (PrevV)
            allocate (PrevV(1:NumLon, 1:NumLat))
            if (allocated(PrevP)) deallocate (PrevP)
            allocate (PrevP(1:NumLon, 1:NumLat))

! These have information from the next snap.
            if (allocated(NextLon)) deallocate (NextLon)
            allocate (NextLon(1:NumLon, 1:NumLat))
            if (allocated(NextLat)) deallocate (NextLat)
            allocate (NextLat(1:NumLon, 1:NumLat))
            if (allocated(NextU)) deallocate (NextU)
            allocate (NextU(1:NumLon, 1:NumLat))
            if (allocated(NextV)) deallocate (NextV)
            allocate (NextV(1:NumLon, 1:NumLat))
            if (allocated(NextP)) deallocate (NextP)
            allocate (NextP(1:NumLon, 1:NumLat))

! Loop over the snaps.
            do IS = 1, 2
               if (IS == 1) then
                  CurrSnap = NWS13(IG)%PrevSnap
               elseif (IS == 2) then
                  CurrSnap = NWS13(IG)%NextSnap
               end if

! 1. Read this snap from the wind file.

! Read mesh.
               call check_err(NF90_INQ_VARID(NC_IDG, "lon", NC_VAR))
               call check_err(NF90_GET_VAR(NC_IDG, NC_VAR, Lon(:, :), &
                                           START=[1, 1, CurrSnap], COUNT=[NumLon, NumLat, 1]))
               call check_err(NF90_INQ_VARID(NC_IDG, "lat", NC_VAR))
               call check_err(NF90_GET_VAR(NC_IDG, NC_VAR, Lat(:, :), &
                                           START=[1, 1, CurrSnap], COUNT=[NumLon, NumLat, 1]))

! Read winds and pressures.
! U-wind component
               call check_err(NF90_INQ_VARID(NC_IDG, "U10", NC_VAR))
               call check_err(NF90_GET_VAR(NC_IDG, NC_VAR, U(:, :), &
                                           START=[1, 1, CurrSnap], COUNT=[NumLon, NumLat, 1]))
! Apply scale/offset if attributes are present
               NC_ERR = NF90_GET_ATT(NC_IDG, NC_VAR, "_FillValue", var_scale_attr)
               if (NC_ERR == 0) then
                  do ILN = 1, NumLon
                     do ILT = 1, NumLat
                        if (abs(U(ILN, ILT) - var_scale_attr) <= EPS) then
                           U(ILN, ILT) = D_QNAN
                        end if
                     end do
                  end do
               end if
               NC_ERR = NF90_GET_ATT(NC_IDG, NC_VAR, "scale_factor", var_scale_attr)
               if (NC_ERR == 0) U = U*var_scale_attr
               NC_ERR = NF90_GET_ATT(NC_IDG, NC_VAR, "add_offset", var_offset_attr)
               if (NC_ERR == 0) U = U + var_offset_attr
! V-wind component
               call check_err(NF90_INQ_VARID(NC_IDG, "V10", NC_VAR))
               call check_err(NF90_GET_VAR(NC_IDG, NC_VAR, V(:, :), &
                                           START=[1, 1, CurrSnap], COUNT=[NumLon, NumLat, 1]))
! Apply scale/offset if attributes are present
               NC_ERR = NF90_GET_ATT(NC_IDG, NC_VAR, "_FillValue", var_scale_attr)
               if (NC_ERR == 0) then
                  do ILN = 1, NumLon
                     do ILT = 1, NumLat
                        if (abs(V(ILN, ILT) - var_scale_attr) <= eps) then
                           V(ILN, ILT) = D_QNAN
                        end if
                     end do
                  end do
               end if
               NC_ERR = NF90_GET_ATT(NC_IDG, NC_VAR, "scale_factor", var_scale_attr)
               if (NC_ERR == 0) V = V*var_scale_attr
               NC_ERR = NF90_GET_ATT(NC_IDG, NC_VAR, "add_offset", var_offset_attr)
               if (NC_ERR == 0) V = V + var_offset_attr
! Pressure
               call check_err(NF90_INQ_VARID(NC_IDG, "PSFC", NC_VAR))
               call check_err(NF90_GET_VAR(NC_IDG, NC_VAR, P(:, :), &
                                           START=[1, 1, CurrSnap], COUNT=[NumLon, NumLat, 1]))
! Apply scale/offset if attributes are present
               NC_ERR = NF90_GET_ATT(NC_IDG, NC_VAR, "_FillValue", var_scale_attr)
               if (NC_ERR == 0) then
                  do ILN = 1, NumLon
                     do ILT = 1, NumLat
                        if (abs(P(ILN, ILT) - var_scale_attr) < eps) then
                           P(ILN, ILT) = D_QNAN
                        end if
                     end do
                  end do
               end if
               NC_ERR = NF90_GET_ATT(NC_IDG, NC_VAR, "scale_factor", var_scale_attr)
               if (NC_ERR == 0) P = P*var_scale_attr
               NC_ERR = NF90_GET_ATT(NC_IDG, NC_VAR, "add_offset", var_offset_attr)
               if (NC_ERR == 0) P = P + var_offset_attr

! Adjust the wind velocity.
               U = U*NWS13WindMultiplierReal
               V = V*NWS13WindMultiplierReal

! Convert the pressures into meters of water.
               P = P*100.d0/RHOWAT0/G

! Try to read the storm center.
               NC_ERR = NF90_INQ_VARID(NC_IDG, "clon", NC_VAR)
               if (NC_ERR == NF90_NOERR) then
                  NC_ERR = NF90_GET_VAR(NC_IDG, NC_VAR, CLon(:), &
                                        START=[CurrSnap], COUNT=[1])
                  FoundEye = .true.
               end if
               NC_ERR = NF90_INQ_VARID(NC_IDG, "clat", NC_VAR)
               if (NC_ERR == NF90_NOERR) then
                  NC_ERR = NF90_GET_VAR(NC_IDG, NC_VAR, CLat(:), &
                                        START=[CurrSnap], COUNT=[1])
                  FoundEye = .true.
               end if

! Save the snaps.
               if (IS == 1) then
                  PrevLon = Lon
                  PrevLat = Lat
                  PrevU = U
                  PrevV = V
                  PrevP = P
                  if (FoundEye) then
                     PrevCLon = CLon
                     PrevCLat = CLat
                  end if
               elseif (IS == 2) then
                  NextLon = Lon
                  NextLat = Lat
                  NextU = U
                  NextV = V
                  NextP = P
                  if (FoundEye) then
                     NextCLon = CLon
                     NextCLat = CLat
                  end if
               end if
            end do

! 2. Interpolate the wind snaps to start/end of current interval.

            TimeInterpFactor = (TimeLoc - NWS13(IG)%PrevTime) &
                               /(NWS13(IG)%NextTime - NWS13(IG)%PrevTime)
!          IF(MYPROC.EQ.0)THEN
!            write(6,*) TimeInterpFactor,TimeLoc,NWS13(IG)%PrevTime,
!       &               NWS13(IG)%NextTime
!          ENDIF
            Lon = PrevLon + (NextLon - PrevLon)*TimeInterpFactor
            Lat = PrevLat + (NextLat - PrevLat)*TimeInterpFactor
            U = PrevU + (NextU - PrevU)*TimeInterpFactor
            V = PrevV + (NextV - PrevV)*TimeInterpFactor
            P = PrevP + (NextP - PrevP)*TimeInterpFactor
! Time Interpolate scalar wind
            PrevWs = sqrt(PrevU**2 + PrevV**2)
            NextWs = sqrt(NextU**2 + NextV**2)
            Ws = PrevWs + (NextWs - PrevWs)*TimeInterpFactor
            if (FoundEye) then
               CLon = PrevCLon + (NextCLon - PrevCLon)*TimeInterpFactor
               CLat = PrevCLat + (NextCLat - PrevCLat)*TimeInterpFactor
            end if

! 3. Interpolate the wind snaps onto the ADCIRC mesh.
            if (IG > 1) then
               call NWS13INTERP(NumLon, NumLat, Lon, Lat, ILon, ILat, W2)
            else
               W2 = W(:, 1:4)
               ILon = int(W(:, 5))
               ILat = int(W(:, 6))
            end if
            do IV = 1, NP
               if (ILon(IV) > 0) then
! do not apply if we pick up a nan
                  if (abs((Ws(ILon(IV), ILat(IV)) &
                           + Ws(ILon(IV) + 1, ILat(IV)) &
                           + Ws(ILon(IV) + 1, ILat(IV) + 1) &
                           + Ws(ILon(IV), ILat(IV) + 1)) - &
                          (Ws(ILon(IV), ILat(IV)) &
                           + Ws(ILon(IV) + 1, ILat(IV)) &
                           + Ws(ILon(IV) + 1, ILat(IV) + 1) &
                           + Ws(ILon(IV), ILat(IV) + 1))) < eps) then
                     UU(IV) = W2(IV, 1)*U(ILon(IV), ILat(IV)) &
                              + W2(IV, 2)*U(ILon(IV) + 1, ILat(IV)) &
                              + W2(IV, 3)*U(ILon(IV) + 1, ILat(IV) + 1) &
                              + W2(IV, 4)*U(ILon(IV), ILat(IV) + 1)
                     VV(IV) = W2(IV, 1)*V(ILon(IV), ILat(IV)) &
                              + W2(IV, 2)*V(ILon(IV) + 1, ILat(IV)) &
                              + W2(IV, 3)*V(ILon(IV) + 1, ILat(IV) + 1) &
                              + W2(IV, 4)*V(ILon(IV), ILat(IV) + 1)
                     Wind(IV) = W2(IV, 1)*Ws(ILon(IV), ILat(IV)) &
                                + W2(IV, 2)*Ws(ILon(IV) + 1, ILat(IV)) &
                                + W2(IV, 3)*Ws(ILon(IV) + 1, ILat(IV) + 1) &
                                + W2(IV, 4)*Ws(ILon(IV), ILat(IV) + 1)
                  end if
! do not apply if we pick up a nan
                  if (abs((P(ILon(IV), ILat(IV)) &
                           + P(ILon(IV) + 1, ILat(IV)) &
                           + P(ILon(IV) + 1, ILat(IV) + 1) &
                           + P(ILon(IV), ILat(IV) + 1)) - &
                          (P(ILon(IV), ILat(IV)) &
                           + P(ILon(IV) + 1, ILat(IV)) &
                           + P(ILon(IV) + 1, ILat(IV) + 1) &
                           + P(ILon(IV), ILat(IV) + 1))) <= eps) then
                     PP(IV) = W2(IV, 1)*P(ILon(IV), ILat(IV)) &
                              + W2(IV, 2)*P(ILon(IV) + 1, ILat(IV)) &
                              + W2(IV, 3)*P(ILon(IV) + 1, ILat(IV) + 1) &
                              + W2(IV, 4)*P(ILon(IV), ILat(IV) + 1)
                  end if
! Adjust U/V Based on scalar wind magnitude
                  if (abs(UU(iV)) <= eps .and. abs(VV(iV)) <= eps) then
                     sWdir = 0d0
                  else
                     sWdir = atan2(UU(iV), VV(iV))
                     UU(IV) = Wind(IV)*sin(sWDir)
                     VV(IV) = Wind(IV)*cos(sWDir)
                  end if
! Casey 200528
                  PRN2(IV) = PP(IV)
                  WVNX2(IV) = UU(IV)
                  WVNY2(IV) = VV(IV)
               end if
            end do

! Casey 200528
!       PRN2 = PP
!       WVNX2 = UU
!       WVNY2 = VV

            if (FoundEye) then
               if ((NWS13GroupForPowellInt == 0) &
                   .or. (NWS13GroupForPowellInt == IG)) then
! Assign the storm center to the next spot in the array.
                  EyeLonR(3) = CLon(1)
                  EyeLatR(3) = CLat(1)
                  PowellGroupTemp = IG
               end if
            end if

         end do ! IG=1,NumGroup

      end if

! Close the file.
      call check_err(NF90_CLOSE(NC_ID))

      if (PowellGroupTemp /= PowellGroup) then
! If this is a different group then we have been using
! for the Powell wind drag scheme, then check whether
! the storm center in this new group is "close enough"
! to what we have been using.  If not, then reset the
! previous storm centers to force the scheme to use
! Garratt for the next few snaps.  It will eventually
! switch back to Powell for this new group / storm.
         if (sqrt((EyeLonR(3) - EyeLonR(2))**2.d0 &
                  + (EyeLatR(3) - EyeLatR(2))**2.d0) &
             > 2.d0) then
            EyeLonR(1) = 0.d0
            EyeLatR(1) = 0.d0
            EyeLonR(2) = 0.d0
            EyeLatR(2) = 0.d0
         end if
         PowellGroup = PowellGroupTemp
      end if

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()

!-----------------------------------------------------------------------
   end subroutine NWS13GET
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
   subroutine NWS13INTERP(NumLon, NumLat, Lon, Lat, ILon, ILat, W)
!-----------------------------------------------------------------------
      use ADC_CONSTANTS, only: RAD2DEG
      use kdtree2_module, only: KDTREE2, KDTREE2_RESULT, kdtree2_create, &
                                kdtree2_n_nearest, kdtree2_destroy
      use MESH, only: NP, SLAM, SFEA

      implicit none

      integer, intent(OUT) :: ILat(NP)
      integer, intent(OUT) :: ILon(NP)
      integer, intent(IN)  :: NumLat
      integer, intent(IN)  :: NumLon

      real(8), intent(IN)  :: Lon(NumLon, NumLat)
      real(8), intent(IN)  :: Lat(NumLon, NumLat)
      real(8), intent(INOUT) :: W(NP, 4)

      integer :: Counter
      integer :: ElemNumber
      integer :: IA
      integer :: IE
      integer :: IO
      integer :: IV
      integer :: NumResult

      logical :: ElemFound

      real(8), allocatable :: ElemCenter(:, :)
      real(8) :: AdcLat
      real(8) :: AdcLon
      real(8) :: MaxLat
      real(8) :: MaxLon
      real(8) :: MinLat
      real(8) :: MinLon
      real(8) :: W0
      real(8) :: W1
      real(8) :: W2
      real(8) :: W3
      real(8) :: W4

      type(KDTREE2), pointer :: kdTree
      type(KDTREE2_RESULT), allocatable :: kdResult(:)

      call setMessageSource("nws13interp")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      NumResult = (NumLon - 1)*(NumLat - 1)
      NumResult = min(NumResult, 1000)

      if (allocated(ElemCenter)) deallocate (ElemCenter)
      allocate (ElemCenter(1:2, 1:(NumLon - 1)*(NumLat - 1)))
      Counter = 0
      do IO = 1, NumLon - 1
         do IA = 1, NumLat - 1
            Counter = Counter + 1
            ElemCenter(1, Counter) = 0.25d0*( &
                                     Lon(IO, IA) + Lon(IO + 1, IA) + &
                                     Lon(IO + 1, IA + 1) + Lon(IO, IA + 1))
            ElemCenter(2, Counter) = 0.25d0*( &
                                     Lat(IO, IA) + Lat(IO + 1, IA) + &
                                     Lat(IO + 1, IA + 1) + Lat(IO, IA + 1))
         end do
      end do
      kdTree => kdtree2_create(ElemCenter, rearrange=.true., &
                               sort=.true.)
      if (allocated(kdResult)) deallocate (kdResult)
      allocate (kdResult(1:NumResult))

      ILon = 0
      ILat = 0
      W = 0.d0

      MaxLat = maxval(Lat)
      MaxLon = maxval(Lon)
      MinLat = minval(Lat)
      MinLon = minval(Lon)
      do IV = 1, NP
         AdcLat = RAD2DEG*SFEA(IV)
         AdcLon = RAD2DEG*SLAM(IV)
         if (AdcLon < MinLon) cycle
         if (AdcLon > MaxLon) cycle
         if (AdcLat < MinLat) cycle
         if (AdcLat > MaxLat) cycle
         ElemFound = .false.
         call kdtree2_n_nearest(tp=kdTree, &
                                qv=[AdcLon, AdcLat], &
                                nn=NumResult, results=kdResult)
         IE = 1
         do while ((.not. ElemFound) .and. (IE <= NumResult))
            ElemNumber = kdResult(IE)%idx
            if (mod(ElemNumber, NumLat - 1) == 0) then
               IO = ElemNumber/(NumLat - 1)
               IA = NumLat - 1
            else
               IO = floor(real(ElemNumber)/real(NumLat - 1)) + 1
               IA = mod(ElemNumber, NumLat - 1)
            end if
            W0 = (Lon(IO, IA)*Lat(IO + 1, IA) &
                  - Lon(IO + 1, IA)*Lat(IO, IA)) &
                 + (Lon(IO + 1, IA)*Lat(IO + 1, IA + 1) &
                    - Lon(IO + 1, IA + 1)*Lat(IO + 1, IA)) &
                 + (Lon(IO + 1, IA + 1)*Lat(IO, IA + 1) &
                    - Lon(IO, IA + 1)*Lat(IO + 1, IA + 1)) &
                 + (Lon(IO, IA + 1)*Lat(IO, IA) &
                    - Lon(IO, IA)*Lat(IO, IA + 1))
            W0 = abs(W0)
            W1 = (AdcLon*Lat(IO + 1, IA) &
                  - Lon(IO + 1, IA)*AdcLat) &
                 + (Lon(IO + 1, IA)*Lat(IO + 1, IA + 1) &
                    - Lon(IO + 1, IA + 1)*Lat(IO + 1, IA)) &
                 + (Lon(IO + 1, IA + 1)*Lat(IO, IA + 1) &
                    - Lon(IO, IA + 1)*Lat(IO + 1, IA + 1)) &
                 + (Lon(IO, IA + 1)*AdcLat &
                    - AdcLon*Lat(IO, IA + 1))
            W1 = abs(W1)
            W2 = (Lon(IO, IA)*AdcLat &
                  - AdcLon*Lat(IO, IA)) &
                 + (AdcLon*Lat(IO + 1, IA + 1) &
                    - Lon(IO + 1, IA + 1)*AdcLat) &
                 + (Lon(IO + 1, IA + 1)*Lat(IO, IA + 1) &
                    - Lon(IO, IA + 1)*Lat(IO + 1, IA + 1)) &
                 + (Lon(IO, IA + 1)*Lat(IO, IA) &
                    - Lon(IO, IA)*Lat(IO, IA + 1))
            W2 = abs(W2)
            W3 = (Lon(IO, IA)*Lat(IO + 1, IA) &
                  - Lon(IO + 1, IA)*Lat(IO, IA)) &
                 + (Lon(IO + 1, IA)*AdcLat &
                    - AdcLon*Lat(IO + 1, IA)) &
                 + (AdcLon*Lat(IO, IA + 1) &
                    - Lon(IO, IA + 1)*AdcLat) &
                 + (Lon(IO, IA + 1)*Lat(IO, IA) &
                    - Lon(IO, IA)*Lat(IO, IA + 1))
            W3 = abs(W3)
            W4 = (Lon(IO, IA)*Lat(IO + 1, IA) &
                  - Lon(IO + 1, IA)*Lat(IO, IA)) &
                 + (Lon(IO + 1, IA)*Lat(IO + 1, IA + 1) &
                    - Lon(IO + 1, IA + 1)*Lat(IO + 1, IA)) &
                 + (Lon(IO + 1, IA + 1)*AdcLat &
                    - AdcLon*Lat(IO + 1, IA + 1)) &
                 + (AdcLon*Lat(IO, IA) &
                    - Lon(IO, IA)*AdcLat)
            W4 = abs(W4)
            if ((W1 + W2 + W3 + W4) <= (W0*2.001d0)) then
               ElemFound = .true.
               exit
            end if
            IE = IE + 1
         end do
         if (ElemFound) then
            ILon(IV) = IO
            ILat(IV) = IA
            W(IV, 1) = W1/(W0*2.d0)
            W(IV, 2) = W2/(W0*2.d0)
            W(IV, 3) = W3/(W0*2.d0)
            W(IV, 4) = W4/(W0*2.d0)
         end if

      end do

      call kdtree2_destroy(tp=kdTree)

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()

!-----------------------------------------------------------------------
   end subroutine NWS13INTERP
!-----------------------------------------------------------------------

   integer function JULIAN(YEAR, MONTH, DAY)
      implicit none
      integer, intent(IN) :: YEAR, MONTH, DAY
      call setMessageSource("julian")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif
      JULIAN = DAY - 32075 + 1461*(YEAR + 4800 + (MONTH - 14)/12)/4 + 367 &
               *(MONTH - 2 - (MONTH - 14)/12*12)/12 &
               - 3*((YEAR + 4900 + (MONTH - 14)/12)/100)/4
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
   end function JULIAN

   real(8) function FRACTIONTIME(HOUR, MINUTE, SECOND)
      implicit none
      integer, intent(IN) :: HOUR, MINUTE, SECOND
      call setMessageSource("fractiontime")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif
      FRACTIONTIME = (dble(HOUR) + (dble(MINUTE) + dble(SECOND)/60.d0)/60.d0)/24.d0
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
   end function FRACTIONTIME

end module mod_nws13
#endif
