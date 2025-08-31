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
!  MODULE mod_nws13
!-----------------------------------------------------------------------
!> @author JC Detrich, NC State University
!> @author Alexander Crosby, Oceanweather Inc., alexc@oceanweather.com
!> @author Zach Cobell, The Water Institute, zcobell@thewaterinstitute.org
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
!> Wind fields should be provided in groups ranked from lowest priority (1) to
!> highest priority (n). The priority with which the wind fields are interpolated to the
!> mesh is determined by rank.
!>
!> Wind fields are interpolated onto each node using a bilinear interpolation
!> scheme. If the mesh lies outside of all wind fields, it is set to
!> PRDEFLT (usually 1013mb) and 0.0 m/s velocity.

!-----------------------------------------------------------------------
module mod_nws13
!-----------------------------------------------------------------------
   use mod_datetime, only: t_datetime, null_datetime
   use mod_nws13_data, only: t_nws13_group

   implicit none

   private

   !> @brief Path to the NetCDF wind input file (default: fort.22.nc)
   character(LEN=1000) :: NWS13File = "fort.22.nc"

   !> @brief Group number to use for Powell wind drag calculation (0 = auto-select)
   integer :: NWS13GroupForPowell = 0

   !> @brief Currently active group for Powell wind drag scheme
   integer :: PowellGroup = 0

   !> @brief Multiplier applied to wind velocities for sensitivity analysis
   real(8) :: NWS13WindMultiplier = 1.0d0

   real(8), allocatable :: WX1(:), WY1(:), PR1(:)
   real(8), allocatable :: WX2(:), WY2(:), PR2(:)

   !> @brief Reference time for cold start initialization from namelist
   type(t_datetime) :: NWS13ColdStart

   type(t_datetime) :: PrevInterpTime !< Date/time of the prev time we interpolated
   type(t_datetime) :: NextInterpTime !< Date/time of the next time we will interpolate

   !> @brief Array of NWS13 wind groups from NetCDF file
   type(t_nws13_group), allocatable :: wind_groups(:)

   !> @brief Group indices sorted by priority rank
   integer, allocatable :: sorted_group_indices(:)

   public :: NWS13INIT, NWS13GET, nws13_set_namelist_parameters

contains

!-----------------------------------------------------------------------
!> @brief Sets NWS13 module parameters from namelist input values
!> @param[in] NWS13Filename_in        Input wind file path
!> @param[in] NWS13ColdStart_in       Cold start reference time string
!> @param[in] NWS13WindMultiplier_in  Wind velocity multiplier factor
!> @param[in] NWS13GroupForPowell_in  Group ID for Powell wind drag scheme
!-----------------------------------------------------------------------
   subroutine nws13_set_namelist_parameters(NWS13Filename_in, &
                                            NWS13ColdStart_in, &
                                            NWS13WindMultiplier_in, &
                                            NWS13GroupForPowell_in)
      use global, only: logMessage, INFO, scratchMessage
      implicit none

      character(len=*), intent(in) :: NWS13Filename_in
      character(len=*), intent(in) :: NWS13ColdStart_in
      real(8), intent(in) :: NWS13WindMultiplier_in
      integer, intent(in) :: NWS13GroupForPowell_in

      ! Set the parameters from the namelist.
      NWS13File = trim(adjustl(NWS13Filename_in))
      NWS13ColdStart = t_datetime(trim(adjustl(NWS13ColdStart_in)), &
                                  ["%Y%m%d.%H%M%S    ", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S"])
      NWS13WindMultiplier = NWS13WindMultiplier_in
      NWS13GroupForPowell = NWS13GroupForPowell_in

      call logMessage(INFO, "Using NWS13 file: "//trim(adjustl(NWS13File)))
      call logMessage(INFO, "NWS13 Cold Start Date: "//trim(NWS13ColdStart%to_iso_string()))
      write (scratchMessage, '(A,F0.4)') "NWS13 Wind Multiplier: ", NWS13WindMultiplier
      call logMessage(INFO, scratchMessage)

   end subroutine nws13_set_namelist_parameters

!-----------------------------------------------------------------------
!> Initializes reading data from the Oceanweather (OWI) NetCDF wind/pre fields
!-----------------------------------------------------------------------
   subroutine NWS13INIT()
!-----------------------------------------------------------------------
      use netcdf, only: NF90_GET_ATT, NF90_OPEN, NF90_NOWRITE, NF90_INQ_NCID, &
                        NF90_INQUIRE_DIMENSION, NF90_INQ_VARID, NF90_GET_VAR, NF90_CLOSE, NF90_INQ_DIMID, &
                        NF90_INQUIRE_VARIABLE, NF90_INQ_GRPS, NF90_INQ_GRPNAME
      use netcdf_error, only: check_err
      use adc_constants, only: PRBCKGRND, RhoWat0, G
      use mesh, only: np
      use global, only: allMessage, screenMessage, setMessageSource, unsetMessageSource, &
                        logMessage, INFO, ERROR, scratchMessage
#if defined(OWIWIND_TRACE)    || defined(ALL_TRACE)
      use global, only: DEBUG
#endif

#ifdef CMPI
      use messenger, only: msg_fini
#endif

      implicit none

      integer, parameter :: max_groups = 1000
      integer :: NumGroup
      integer :: nc_id
      integer :: IG, IGS
      integer, allocatable :: group_ids(:)

      call setMessageSource("nws13init")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      ! Open the NetCDF file.
      call check_err(NF90_OPEN(trim(adjustl(NWS13File)), NF90_NOWRITE, NC_ID))

      ! Discover groups in the NetCDF file
      allocate (group_ids(max_groups))
      call check_err(NF90_INQ_GRPS(NC_ID, NumGroup, group_ids))
      if (NumGroup == 0) then
         call allMessage(ERROR, "NWS13: No groups found in NetCDF file.")
#ifdef CMPI
         call msg_fini()
#endif
         call exit(1)
      end if

      PrevInterpTime = null_datetime()
      NextInterpTime = null_datetime()

      ! Allocate group_ids array and get the IDs
      allocate (wind_groups(1:NumGroup))
      do IG = 1, NumGroup
         wind_groups(IG) = t_nws13_group(group_ids(IG), IG)
      end do
      deallocate (group_ids)

      allocate (sorted_group_indices(NumGroup))
      sorted_group_indices = sort_groups_by_rank(wind_groups)

      call logMessage(INFO, "Groups will be processed in rank order (highest priority first):")
      do IGS = 1, NumGroup
         IG = sorted_group_indices(IGS)
         write (scratchMessage, '(A,I0,A,A,A,I0,A)') "  Priority ", IGS, ": '", &
            trim(wind_groups(IG)%group_name), "' (rank ", wind_groups(IG)%rank, ")"
         call logMessage(INFO, scratchMessage)
      end do

      call check_err(NF90_CLOSE(NC_ID))

      allocate (wx1(np), wy1(np), pr1(np))
      allocate (wx2(np), wy2(np), pr2(np))

      wx1 = 0d0
      wy1 = 0d0
      pr1 = PRBCKGRND*100d0/RHOWAT0/G
      wx2 = 0d0
      wy2 = 0d0
      pr2 = PRBCKGRND*100d0/RHOWAT0/G

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()

!-----------------------------------------------------------------------
   end subroutine NWS13INIT
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Updates group snapshots and marks active groups
!> @param[in] current_time Current simulation time
!> @param[inout] groups Array of NWS13 groups to update
!> @return Next wind field time as t_datetime
!-----------------------------------------------------------------------
   subroutine update_active_groups(current_time)
      use mod_datetime, only: t_timedelta, operator(+), operator(-), operator(<), operator(==)
      implicit none
      type(t_datetime), intent(in) :: current_time
      integer :: ig

      do ig = 1, size(wind_groups)
         call wind_groups(ig)%advance(current_time)
         if(PrevInterpTime < wind_groups(ig)%PrevFileTime .or. PrevInterpTime == null_datetime())then
             PrevInterpTime = wind_groups(ig)%PrevFileTime
             NextInterpTime = wind_groups(ig)%PrevFileTime + t_timedelta(seconds=600)
         end if
      end do
   end subroutine update_active_groups

   logical pure function nws13_data_started(current_time) result(is_started)
      use mod_datetime, only: operator(>)
      implicit none
      type(t_datetime), intent(in) :: current_time
      integer :: i

      is_started = .false.

      do i = 1, size(wind_groups)
         if (current_time > wind_groups(i)%PrevFileTime) is_started = .true.
      end do
   end function nws13_data_started

   subroutine time_interpolate_wind(CurrentTime, NextInterpTime, WXOUT, WYOUT, PROUT)
      use mod_datetime, only: t_timedelta, operator(-)
      use mesh, only: np
      implicit none
      type(t_datetime), intent(in) :: CurrentTime, NextInterpTime
      real(8), intent(out) :: WXOUT(NP), WYOUT(NP), PROUT(NP)
      real(8) :: wtratio
      type(t_timedelta) :: dt_now, dt_interp

      dt_now = NextInterpTime - CurrentTime
      dt_interp = NextInterpTime - PrevInterpTime
      wtratio = dble(dt_now%total_milliseconds())/dble(dt_interp%total_milliseconds())

      WXOUT = (1.0d0 - wtratio)*WX1 + wtratio*WX2
      WYOUT = (1.0d0 - wtratio)*WY1 + wtratio*WY2
      PROUT = (1.0d0 - wtratio)*PR1 + wtratio*PR2

   end subroutine time_interpolate_wind

   subroutine update_wind_data(CurrentTime, EyeLonR, EyeLatR, FoundEye, PowellGroupTemp)
      use ADC_CONSTANTS, only: G, RHOWAT0, PRBCKGRND
      use netcdf, only: NF90_GET_VAR, NF90_INQ_VARID, NF90_INQ_NCID, &
                        NF90_INQUIRE_DIMENSION, NF90_GET_ATT, NF90_OPEN, NF90_NOWRITE, &
                        NF90_CLOSE, NF90_INQ_DIMID
      use netcdf_error, only: check_err
      use mod_nws13_data, only: null_flag_value
      use mod_datetime, only: t_timedelta, operator(+)
      implicit none

      real(8), parameter :: EPS = epsilon(1d0)

      type(t_datetime), intent(in) :: CurrentTime
      real(8), intent(inout) :: EyeLonR(3), EyeLatR(3)
      logical, intent(inout) :: FoundEye
      integer, intent(inout) :: PowellGroupTemp
      integer :: ig, ig_sorted
      integer :: nc_id

      ! First, initialize the ADCIRC wind fields with flag values.
      ! These will be replaced with actual values or defaults after processing all groups.
      WX2 = null_flag_value
      WY2 = null_flag_value
      PR2 = null_flag_value
      PrevInterpTime = NextInterpTime

      ! Open a connection to the wind file.
      call check_err(NF90_OPEN(trim(adjustl(NWS13File)), NF90_NOWRITE, NC_ID))

      ! Update active groups and get next time from active groups
      call update_active_groups(CurrentTime)

      NextInterpTime = PrevInterpTime + t_timedelta(seconds=600)
      PowellGroupTemp = PowellGroup

      ! Increment the storm centers.
      EyeLonR(1) = EyeLonR(2)
      EyeLatR(1) = EyeLatR(2)
      EyeLonR(2) = EyeLonR(3)
      EyeLatR(2) = EyeLatR(3)
      EyeLonR(3) = 0.d0
      EyeLatR(3) = 0.d0
      FoundEye = .false.

      ! Process each wind group in rank order
      do IG_sorted = 1, size(wind_groups)
         ig = sorted_group_indices(ig_sorted)
         call wind_groups(ig)%process(CurrentTime, WX2, WY2, PR2, &
                                NWS13GroupForPowell, PowellGroupTemp, EyeLonR, &
                                EyeLatR, FoundEye, NWS13WindMultiplier)
      end do

      ! Replace any remaining flag values with default background values
      where (abs(PR2 - null_flag_value) <= eps)
         PR2 = PRBCKGRND*100d0/RHOWAT0/G
         WX2 = 0.d0
         WY2 = 0.d0
      end where

      ! Close the file.
      call check_err(NF90_CLOSE(NC_ID))

   end subroutine update_wind_data

!-----------------------------------------------------------------------
!> @brief Reads multi-grid wind and pressure fields from the NetCDF file and
!>        interpolates them to the adcirc mesh
!> @param[in]    timeloc  model time
!> @param[inout] wx       wind speed, x-direction
!> @param[inout] wy       wind speed, y-direction
!> @param[inout] p        atmospheric pressure
!> @param[inout] wtime2   time of next new field
!> @param[inout] eyelonr  storm center longitude
!> @param[inout] eyelatr  storm center latitude
!> @param[out]   foundeye boolean identifying if file provides storm center
!-----------------------------------------------------------------------
   subroutine NWS13GET(TimeLoc, WXOUT, WYOUT, PROUT, EyeLonR, EyeLatR, FoundEye)
!-----------------------------------------------------------------------
      use MESH, only: NP
      use mod_datetime, only: t_timedelta, null_datetime, operator(+), operator(-), operator(<), operator(==), &
                              operator(>=), operator(<=), operator(/=)
      use global, only: screenMessage, setMessageSource, unsetMessageSource
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      use global, only: allMessage, DEBUG
#endif
      implicit none

      real(8), intent(in)    :: TimeLoc
      real(8), intent(inout) :: EyeLatR(3)
      real(8), intent(inout) :: EyeLonR(3)
      logical, intent(inout) :: FoundEye
      real(8), intent(out) :: PROUT(NP), WXOUT(NP), WYOUT(NP)

      integer :: PowellGroupTemp

      type(t_datetime)  :: CurrentTime

      CurrentTime = NWS13Coldstart + t_timedelta(milliseconds=int(TimeLoc*1000d0))
      if (.not. nws13_data_started(CurrentTime)) then
         return
      elseif (CurrentTime < NextInterpTime) then
         call time_interpolate_wind(CurrentTime, NextInterpTime, WXOUT, WYOUT, PROUT)
      else
         call update_wind_data(CurrentTime, EyeLonR, EyeLatR, FoundEye, PowellGroupTemp)
         call time_interpolate_wind(CurrentTime, NextInterpTime, WXOUT, WYOUT, PROUT)

         if (PowellGroupTemp /= PowellGroup) then
            ! If this is a different group then we have been using
            ! for the Powell wind drag scheme, then check whether
            ! the storm center in this new group is "close enough"
            ! to what we have been using.  If not, then reset the
            ! previous storm centers to force the scheme to use
            ! Garratt for the next few snaps.  It will eventually
            ! switch back to Powell for this new group / storm.
            if (sqrt((EyeLonR(3) - EyeLonR(2))**2.d0 + (EyeLatR(3) - EyeLatR(2))**2.d0) > 2.d0) then
               EyeLonR(1) = 0.d0
               EyeLatR(1) = 0.d0
               EyeLonR(2) = 0.d0
               EyeLatR(2) = 0.d0
            end if
            PowellGroup = PowellGroupTemp
         end if

      end if

!-----------------------------------------------------------------------
   end subroutine NWS13GET
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Sorts group indices by their rank values (descending order)
!> @param[in] groups Array of NWS13Type structures containing rank information
!> @return indices Array of group indices sorted by rank
!-----------------------------------------------------------------------
   pure function sort_groups_by_rank(groups) result(indices)
      type(t_nws13_group), intent(in) :: groups(:)
      integer :: indices(size(groups))
      integer :: n, i, j, temp_idx, temp_rank
      integer :: sorted_ranks(size(groups))

      n = size(groups)

      ! Initialize indices and extract ranks
      do i = 1, n
         indices(i) = i
         sorted_ranks(i) = groups(i)%rank
      end do

      do i = 1, n - 1
         do j = 1, n - i
            if (sorted_ranks(j) < sorted_ranks(j + 1)) then
               ! Swap ranks
               temp_rank = sorted_ranks(j)
               sorted_ranks(j) = sorted_ranks(j + 1)
               sorted_ranks(j + 1) = temp_rank

               ! Swap indices
               temp_idx = indices(j)
               indices(j) = indices(j + 1)
               indices(j + 1) = temp_idx
            end if
         end do
      end do
   end function sort_groups_by_rank

end module mod_nws13
