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
   use mod_nws13_data, only: t_nws13

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

   !> @brief Reference time for cold start initialization from namelist
   type(t_datetime) :: NWS13ColdStart

   !> @brief Array of NWS13 wind groups from NetCDF file
   type(t_nws13), allocatable :: NWS13(:)

   !> @brief Group indices sorted by priority rank
   integer, allocatable :: sorted_group_indices(:)

   public :: NWS13INIT, NWS13GET, nws13_set_namelist_parameters

contains

!-----------------------------------------------------------------------
!> @brief Compares two arrays of t_datetime objects for inequality
!> @param[in] arr1 First datetime array
!> @param[in] arr2 Second datetime array
!> @return differ True if arrays differ in size or any element differs
!-----------------------------------------------------------------------
   pure function datetime_arrays_differ(arr1, arr2) result(differ)
      use mod_datetime, only: operator(/=)
      implicit none
      type(t_datetime), intent(in) :: arr1(:), arr2(:)
      logical :: differ
      integer :: i

      ! Check if sizes differ
      if (size(arr1) /= size(arr2)) then
         differ = .true.
         return
      end if

      ! Check if any elements differ
      differ = .false.
      do i = 1, size(arr1)
         if (arr1(i) /= arr2(i)) then
            differ = .true.
            return
         end if
      end do

   end function datetime_arrays_differ

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
      use global, only: allMessage, screenMessage, setMessageSource, unsetMessageSource, &
                        logMessage, INFO, ERROR, scratchMessage
      use mod_datetime, only: t_datetime, t_timedelta, operator(+), operator(-)
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

      ! Allocate group_ids array and get the IDs
      allocate (NWS13(1:NumGroup))
      do IG = 1, NumGroup
         NWS13(IG) = t_nws13(group_ids(IG), IG)
         if (datetime_arrays_differ(NWS13(IG)%SnapTimes, NWS13(1)%SnapTimes)) then
            call allMessage(ERROR, "Inconsistent timing of wind fields detected.")
#ifdef CMPI
            call msg_fini()
#endif
            call exit(1)
         end if
      end do
      deallocate (group_ids)

      allocate (sorted_group_indices(NumGroup))
      sorted_group_indices = sort_groups_by_rank(NWS13)

      call logMessage(INFO, "Groups will be processed in rank order (highest priority first):")
      do IGS = 1, NumGroup
         IG = sorted_group_indices(IGS)
         write (scratchMessage, '(A,I0,A,A,A,I0,A)') "  Priority ", IGS, ": '", &
            trim(NWS13(IG)%group_name), "' (rank ", NWS13(IG)%rank, ")"
         call logMessage(INFO, scratchMessage)
      end do

      call check_err(NF90_CLOSE(NC_ID))

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()

!-----------------------------------------------------------------------
   end subroutine NWS13INIT
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Determines the next wind field time for groups not yet started
!> @param[in] current_time Current simulation time
!> @param[in] groups Array of NWS13 groups
!> @return Next wind field time as t_datetime
!-----------------------------------------------------------------------
   type(t_datetime) function find_next_wind_time_for_inactive_groups(current_time, groups) result(next_time)
      use mod_datetime, only: t_datetime, null_datetime, operator(<), operator(==)
      implicit none
      type(t_datetime), intent(in) :: current_time
      type(t_nws13), intent(in) :: groups(:)
      integer :: ig

      next_time = null_datetime()

      do ig = 1, size(groups)
         if (current_time < groups(ig)%PrevTime) then
            if (next_time == null_datetime() .or. groups(ig)%PrevTime < next_time) then
               next_time = groups(ig)%PrevTime
            end if
         end if
      end do
   end function find_next_wind_time_for_inactive_groups

!-----------------------------------------------------------------------
!> @brief Updates group snapshots and marks active groups
!> @param[in] current_time Current simulation time
!> @param[inout] groups Array of NWS13 groups to update
!> @return Next wind field time as t_datetime
!-----------------------------------------------------------------------
   type(t_datetime) function update_active_groups(current_time, groups) result(next_time)
      use mod_datetime, only: t_datetime, null_datetime, operator(<=), operator(<), operator(==)
      implicit none
      type(t_datetime), intent(in) :: current_time
      type(t_nws13), intent(inout) :: groups(:)
      integer :: ig

      next_time = null_datetime()
      do ig = 1, size(groups)
         call groups(ig)%advance(current_time, next_time)
      end do
   end function update_active_groups

!-----------------------------------------------------------------------
!> @brief Reads multi-grid wind and pressure fields from the NetCDF file and
!>        interpolates them to the adcirc mesh
!> @param[in]    timeloc  model time
!> @param[inout] wvnx2    wind speed, x-direction
!> @param[inout] wvny2    wind speed, y-direction
!> @param[inout] prn2     atmospheric pressure
!> @param[inout] wtime2   time of next new field
!> @param[inout] eyelonr  storm center longitude
!> @param[inout] eyelatr  storm center latitude
!> @param[out]   foundeye boolean identifying if file provides storm center
!-----------------------------------------------------------------------
   subroutine NWS13GET(TimeLoc, WVNX2, WVNY2, PRN2, &
                       WTIME2, EyeLonR, EyeLatR, FoundEye)
!-----------------------------------------------------------------------
      use ADC_CONSTANTS, only: G, RHOWAT0, PRBCKGRND
      use MESH, only: NP
      use netcdf, only: NF90_GET_VAR, NF90_INQ_VARID, NF90_INQ_NCID, &
                        NF90_INQUIRE_DIMENSION, NF90_GET_ATT, NF90_OPEN, NF90_NOWRITE, &
                        NF90_CLOSE, NF90_INQ_DIMID
      use netcdf_error, only: check_err
      use mod_datetime, only: t_timedelta, null_datetime, operator(+), operator(-), operator(<), operator(==), &
                              operator(>=), operator(<=), operator(/=)
      use mod_nws13_data, only: null_flag_value
      use global, only: screenMessage, setMessageSource, unsetMessageSource, &
                        RNDAY
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      use global, only: allMessage, DEBUG
#endif
      implicit none

      real(8), parameter :: EPS = epsilon(1d0)

      real(8), intent(inout) :: EyeLatR(3)
      real(8), intent(inout) :: EyeLonR(3)
      real(8), intent(out)   :: PRN2(NP)
      real(8), intent(in)    :: TimeLoc
      real(8), intent(inout) :: WVNX2(NP)
      real(8), intent(inout) :: WVNY2(NP)
      real(8), intent(out)   :: WTIME2

      logical, intent(out)   :: FoundEye

      integer :: IG, IG_sorted
      integer :: NC_ID
      integer :: PowellGroupTemp

      type(t_datetime)  :: CurrentTime
      type(t_datetime)  :: wtime2_dt, wtime2_dt_active
      type(t_timedelta) :: dt

      call setMessageSource("nws13get")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      CurrentTime = NWS13Coldstart + t_timedelta(milliseconds=int(TimeLoc*1000d0))
      if (CurrentTime < NWS13(1)%PrevTime) then
         WVNX2 = 0d0
         WVNY2 = 0d0
         PRN2 = PRBCKGRND*100d0/RHOWAT0/G
         call unsetMessageSource()
         return
      end if

      ! First, initialize the ADCIRC wind fields with flag values.
      ! These will be replaced with actual values or defaults after processing all groups.
      WVNX2 = null_flag_value
      WVNY2 = null_flag_value
      PRN2 = null_flag_value
      PowellGroupTemp = PowellGroup

      ! Open a connection to the wind file.
      call check_err(NF90_OPEN(trim(adjustl(NWS13File)), NF90_NOWRITE, NC_ID))

      ! Determine active groups and next wind field time
      CurrentTime = NWS13Coldstart + t_timedelta(milliseconds=int(TimeLoc*1000d0))

      ! First check for groups not yet started
      WTIME2_dt = find_next_wind_time_for_inactive_groups(CurrentTime, NWS13)

      ! Update active groups and get next time from active groups
      wtime2_dt_active = update_active_groups(CurrentTime, NWS13)

      ! Use the earlier of the two times
      if (wtime2_dt_active /= null_datetime()) then
         if (WTIME2_dt == null_datetime() .or. wtime2_dt_active < WTIME2_dt) then
            WTIME2_dt = wtime2_dt_active
         end if
      end if

      ! Convert to seconds or use simulation end time
      if (WTIME2_dt == null_datetime()) then
         WTIME2 = RNDAY*86400.d0
      else
         dt = WTIME2_dt - NWS13ColdStart
         WTIME2 = dble(dt%total_milliseconds())/1000.0d0
      end if

      ! If we don't need to include winds from any group (because the
      ! simulation has started before or extended after the information
      ! in the wind file), then we can return.
      if (sum(NWS13(:)%InclSnap) > 0) then

         ! Increment the storm centers.
         EyeLonR(1) = EyeLonR(2)
         EyeLatR(1) = EyeLatR(2)
         EyeLonR(2) = EyeLonR(3)
         EyeLatR(2) = EyeLatR(3)
         EyeLonR(3) = 0.d0
         EyeLatR(3) = 0.d0
         FoundEye = .false.

         ! Process each wind group in rank order
         do IG_sorted = 1, size(NWS13)
            ig = sorted_group_indices(ig_sorted)
            call NWS13(ig)%process(CurrentTime, WVNX2, WVNY2, PRN2, &
                                   NWS13GroupForPowell, PowellGroupTemp, EyeLonR, &
                                   EyeLatR, FoundEye, NWS13WindMultiplier)
         end do

      end if

      ! Replace any remaining flag values with default background values
      where (abs(PRN2 - null_flag_value) <= eps)
         PRN2 = PRBCKGRND*100d0/RHOWAT0/G
         WVNX2 = 0.d0
         WVNY2 = 0.d0
      end where

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
         if (sqrt((EyeLonR(3) - EyeLonR(2))**2.d0 + (EyeLatR(3) - EyeLatR(2))**2.d0) > 2.d0) then
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
!> @brief Sorts group indices by their rank values (descending order)
!> @param[in] groups Array of NWS13Type structures containing rank information
!> @return indices Array of group indices sorted by rank
!-----------------------------------------------------------------------
   pure function sort_groups_by_rank(groups) result(indices)
      type(t_nws13), intent(in) :: groups(:)
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
