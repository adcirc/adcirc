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
!> @details The module is initialized by creating a t_nws13 object using the constructor,
!> then calling the get() method at each time step. The code interpolates wind and
!> pressure data to mesh nodes, and provides storm center positions for Powell drag
!> calculations when available.
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
!> PRBCKGRND (usually 1013mb) and 0.0 m/s velocity.

!-----------------------------------------------------------------------
module mod_nws13
!-----------------------------------------------------------------------
   use mod_datetime, only: t_datetime, null_datetime
   use mod_nws13_data, only: t_nws13_group

   implicit none

   private

   type t_nws13

      private

      character(1000) :: filename !> Path to the NetCDF wind input file (default: fort.22.nc)
      integer :: current_powell_group = 0 !> Currently active group for Powell wind drag scheme
      integer :: requested_powell_group = 0 !> Group number to use for Powell wind drag calculation (0 = auto-select)
      integer :: n_groups = 0
      integer, allocatable :: sorted_group_indices(:) !> Group indices sorted by priority rank
      real(8) :: wind_multiplier = 1.0d0 !> Multiplier applied to wind velocities for sensitivity analysis
      real(8), allocatable :: wx1(:), wx2(:) !> Wind x-component at previous/next time (m/s)
      real(8), allocatable :: wy1(:), wy2(:) !> Wind y-component at previous/next time (m/s)
      real(8), allocatable :: pr1(:), pr2(:) !> Pressure at previous/next time (m of water)
      type(t_datetime) :: coldstart !> Reference time for cold start initialization from namelist
      type(t_datetime) :: previous_interp_time !> Previous interpolation time
      type(t_datetime) :: next_interp_time !> Next interpolation time
      type(t_nws13_group), allocatable :: wind_groups(:) !> Array of wind groups from NetCDF file
      logical :: is_initialized = .false.

   contains
      procedure, pass(self), public  :: get => t_nws13_get_wind_forcing
      procedure, pass(self), private :: update_active_groups => t_nws13_update_active_groups
      procedure, pass(self), private :: check_data_started => t_nws13_check_data_started
      procedure, pass(self), private :: time_interpolate_wind => t_nws13_time_interpolate_wind
      procedure, pass(self), private :: update_wind_data => t_nws13_update_wind_data
   end type t_nws13

   !> @brief Container for namelist parameters from fort.15
   !> @details Stores NWS13 configuration before object initialization
   type t_nws13_namelist_info
      character(1000) :: filename !> Path to NetCDF wind file
      integer :: powell_group !> Group index for Powell drag (0=auto)
      real(8) :: wind_multiplier !> Wind velocity scaling factor
      type(t_datetime) :: coldstart !> Cold start reference time
      character(100), allocatable :: enabled_groups(:) !> Array of group names to enable (empty=all enabled)
      logical :: initialized = .false. !> Flag indicating parameters are set
   end type t_nws13_namelist_info

   interface t_nws13
      module procedure t_nws13_constructor
   end interface t_nws13

   !> @brief Single source for all user input parameters
   !> @details Stores namelist parameters before NWS13 object creation
   type(t_nws13_namelist_info) :: user_input

   public :: nws13_set_namelist_parameters, t_nws13

contains

!-----------------------------------------------------------------------
!> @brief Sets NWS13 module parameters from namelist input values
!>
!> @details Stores user-specified parameters from the fort.15 namelist for later
!> use during NWS13 initialization. Validates and parses the cold start time string
!> into a datetime object. Logs configuration parameters for debugging.
!>
!> @param[in] NWS13Filename_in        Path to NetCDF wind file (default: fort.22.nc)
!> @param[in] NWS13ColdStart_in       Cold start time string (YYYYMMDD.HHMMSS or ISO format)
!> @param[in] NWS13WindMultiplier_in  Wind velocity scaling factor (typically 1.0)
!> @param[in] NWS13GroupForPowell_in  Group index for Powell drag (0=auto-select)
!> @param[in] NWS13EnabledGroups_in   Array of group names to enable (use default value to enable all)
!-----------------------------------------------------------------------
   subroutine nws13_set_namelist_parameters(NWS13Filename_in, &
                                            NWS13ColdStart_in, &
                                            NWS13WindMultiplier_in, &
                                            NWS13GroupForPowell_in, &
                                            NWS13EnabledGroups_in)
      use global, only: logMessage, INFO, scratchMessage
      implicit none

      character(len=*), intent(in) :: NWS13Filename_in
      character(len=*), intent(in) :: NWS13ColdStart_in
      real(8), intent(in) :: NWS13WindMultiplier_in
      integer, intent(in) :: NWS13GroupForPowell_in
      character(len=*), intent(in) :: NWS13EnabledGroups_in(:)

      integer :: i

      ! Set the parameters from the namelist.
      user_input%filename = trim(adjustl(NWS13Filename_in))
      user_input%coldstart = t_datetime(trim(adjustl(NWS13ColdStart_in)), &
                                        ["%Y%m%d.%H%M%S    ", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S"])
      user_input%wind_multiplier = NWS13WindMultiplier_in
      user_input%powell_group = NWS13GroupForPowell_in

      ! Handle enabled groups selection - check if it has meaningful content
      if (size(NWS13EnabledGroups_in) > 0 .and. trim(NWS13EnabledGroups_in(1)) /= '') then
         allocate (user_input%enabled_groups(size(NWS13EnabledGroups_in)))
         do i = 1, size(NWS13EnabledGroups_in)
            user_input%enabled_groups(i) = trim(adjustl(NWS13EnabledGroups_in(i)))
         end do
      end if

      user_input%initialized = .true.

      call logMessage(INFO, "Using NWS13 file: "//trim(adjustl(user_input%filename)))
      call logMessage(INFO, "NWS13 Cold Start Date: "//trim(user_input%coldstart%to_iso_string()))
      write (scratchMessage, '(A,F0.4)') "NWS13 Wind Multiplier: ", user_input%wind_multiplier
      call logMessage(INFO, scratchMessage)

      if (allocated(user_input%enabled_groups)) then
         write (scratchMessage, '(A,I0,A)') "NWS13 Enabled Groups (", size(user_input%enabled_groups), "):"
         call logMessage(INFO, scratchMessage)
         do i = 1, size(user_input%enabled_groups)
            call logMessage(INFO, "  - "//trim(user_input%enabled_groups(i)))
         end do
      else
         call logMessage(INFO, "NWS13: All groups enabled (default behavior)")
      end if

   end subroutine nws13_set_namelist_parameters

!-----------------------------------------------------------------------
!> @brief Initializes reading data from the Oceanweather (OWI) NetCDF wind/pre fields
!>
!> @details Creates and initializes a new t_nws13 object by reading NetCDF
!> wind group metadata from the specified file. Discovers all groups within
!> the file, sorts them by priority rank, and allocates arrays for wind
!> and pressure data interpolation. Sets up initial values using background
!> pressure and zero wind speeds.
!>
!> @return nws Fully initialized t_nws13 object ready for wind forcing calculations
!-----------------------------------------------------------------------
   type(t_nws13) function t_nws13_constructor() result(nws)
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
      if (.not. user_input%initialized) then
         call allMessage(ERROR, "NWS13 Namelist information not initialized")
#ifdef CMPI
         call msg_fini()
#endif
         call exit(1)
      end if

      nws%filename = user_input%filename
      nws%wind_multiplier = user_input%wind_multiplier
      nws%requested_powell_group = user_input%powell_group
      nws%coldstart = user_input%coldstart

      ! Open the NetCDF file.
      call check_err(NF90_OPEN(trim(adjustl(nws%filename)), NF90_NOWRITE, NC_ID))

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

      nws%n_groups = NumGroup
      nws%previous_interp_time = null_datetime()
      nws%next_interp_time = null_datetime()

      ! Allocate group_ids array and get the IDs
      allocate (nws%wind_groups(1:NumGroup))
      do IG = 1, NumGroup
         nws%wind_groups(IG) = t_nws13_group(group_ids(IG), IG)
         ! Apply group selection logic
         call apply_group_selection(nws%wind_groups(IG))
      end do
      deallocate (group_ids)

      allocate (nws%sorted_group_indices(NumGroup))
      nws%sorted_group_indices = sort_groups_by_rank(nws%wind_groups)

      ! Validate group selection and warn about issues
      call validate_group_selection(nws%wind_groups)

      call logMessage(INFO, "Groups will be processed in rank order (highest priority first):")
      do IGS = 1, NumGroup
         IG = nws%sorted_group_indices(IGS)
         if (nws%wind_groups(IG)%enabled) then
            write (scratchMessage, '(A,I0,A,A,A,I0,A)') "  Priority ", IGS, ": '", &
               trim(nws%wind_groups(IG)%group_name), "' (rank ", nws%wind_groups(IG)%rank, ") [ENABLED]"
         else
            write (scratchMessage, '(A,I0,A,A,A,I0,A)') "  Priority ", IGS, ": '", &
               trim(nws%wind_groups(IG)%group_name), "' (rank ", nws%wind_groups(IG)%rank, ") [DISABLED]"
         end if
         call logMessage(INFO, scratchMessage)
      end do

      call check_err(NF90_CLOSE(NC_ID))

      allocate (nws%wx1(np), nws%wy1(np), nws%pr1(np))
      allocate (nws%wx2(np), nws%wy2(np), nws%pr2(np))

      nws%wx1 = 0d0
      nws%wy1 = 0d0
      nws%pr1 = PRBCKGRND*100d0/RHOWAT0/G
      nws%wx2 = 0d0
      nws%wy2 = 0d0
      nws%pr2 = PRBCKGRND*100d0/RHOWAT0/G

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()

!-----------------------------------------------------------------------
   end function t_nws13_constructor
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Updates group snapshots and marks active groups
!>
!> @details Advances all wind groups to the current simulation time by
!> updating their time snapshots. Sets the previous and next interpolation
!> times based on the earliest group times. This ensures all groups are
!> synchronized for the current time step.
!>
!> @param[inout] self    NWS13 object containing wind groups
!> @param[in] current_time Current simulation time
!-----------------------------------------------------------------------
   subroutine t_nws13_update_active_groups(self, current_time)
      use mod_datetime, only: t_timedelta, operator(+), operator(-), operator(<), operator(==)
      implicit none
      class(t_nws13), intent(inout) :: self
      type(t_datetime), intent(in) :: current_time
      integer :: ig

      do ig = 1, size(self%wind_groups)
         call self%wind_groups(ig)%advance(current_time)
         if (self%previous_interp_time < self%wind_groups(ig)%PrevFileTime .or. &
             self%previous_interp_time == null_datetime()) then
            self%previous_interp_time = self%wind_groups(ig)%PrevFileTime
            self%next_interp_time = self%wind_groups(ig)%PrevFileTime + t_timedelta(seconds=600)
         end if
      end do
   end subroutine t_nws13_update_active_groups

!-----------------------------------------------------------------------
!> @brief Checks if wind data has started for current simulation time
!>
!> @details Determines whether any wind group has data available at the
!> current simulation time. Used to avoid processing wind fields before
!> the first available data snapshot.
!>
!> @param[in] self         NWS13 object containing wind groups
!> @param[in] current_time Current simulation time to check
!> @return    is_started   True if wind data is available, false otherwise
!-----------------------------------------------------------------------
   logical pure function t_nws13_check_data_started(self, current_time) result(is_started)
      use mod_datetime, only: operator(>)
      implicit none
      class(t_nws13), intent(in) :: self
      type(t_datetime), intent(in) :: current_time
      integer :: i

      is_started = .false.
      do i = 1, self%n_groups
         if (current_time > self%wind_groups(i)%PrevFileTime) is_started = .true.
      end do
   end function t_nws13_check_data_started

!-----------------------------------------------------------------------
!> @brief Performs temporal interpolation of wind and pressure fields
!>
!> @details Interpolates wind velocity components and pressure between
!> two time snapshots based on the current simulation time. Uses linear
!> interpolation weighted by the time ratio between snapshots.
!>
!> @param[in]  self        NWS13 object containing wind field data
!> @param[in]  CurrentTime Current simulation time
!> @param[out] WXOUT       Interpolated wind x-component at mesh nodes (m/s)
!> @param[out] WYOUT       Interpolated wind y-component at mesh nodes (m/s)
!> @param[out] PROUT       Interpolated pressure at mesh nodes (m of water)
!-----------------------------------------------------------------------
   subroutine t_nws13_time_interpolate_wind(self, CurrentTime, WXOUT, WYOUT, PROUT)
      use mod_datetime, only: t_timedelta, operator(-)
      use mesh, only: np
      implicit none
      class(t_nws13), intent(in) :: self
      type(t_datetime), intent(in) :: CurrentTime
      real(8), intent(out) :: WXOUT(NP), WYOUT(NP), PROUT(NP)
      real(8) :: wtratio
      type(t_timedelta) :: dt_now, dt_interp

      dt_now = self%next_interp_time - CurrentTime
      dt_interp = self%next_interp_time - self%previous_interp_time
      wtratio = dble(dt_now%total_milliseconds())/dble(dt_interp%total_milliseconds())

      WXOUT = (1.0d0 - wtratio)*self%WX1 + wtratio*self%WX2
      WYOUT = (1.0d0 - wtratio)*self%WY1 + wtratio*self%WY2
      PROUT = (1.0d0 - wtratio)*self%PR1 + wtratio*self%PR2

   end subroutine t_nws13_time_interpolate_wind

!-----------------------------------------------------------------------
!> @brief Updates wind field data for the current time step
!>
!> @details Reads new wind and pressure snapshots from all active groups,
!> processes them in priority order, and updates storm center positions
!> for Powell wind drag calculations. Replaces flag values with background
!> conditions for nodes outside all wind grids.
!>
!> @param[inout] self              NWS13 object to update
!> @param[in]    CurrentTime       Current simulation time
!> @param[inout] EyeLonR           Storm center longitude array (degrees)
!> @param[inout] EyeLatR           Storm center latitude array (degrees)
!> @param[inout] FoundEye          Flag indicating if storm eye was located
!> @param[inout] PowellGroupTemp   Group index used for Powell drag scheme
!-----------------------------------------------------------------------
   subroutine t_nws13_update_wind_data(self, CurrentTime, EyeLonR, EyeLatR, FoundEye, PowellGroupTemp)
      use ADC_CONSTANTS, only: G, RHOWAT0, PRBCKGRND
      use netcdf, only: NF90_GET_VAR, NF90_INQ_VARID, NF90_INQ_NCID, &
                        NF90_INQUIRE_DIMENSION, NF90_GET_ATT, NF90_OPEN, NF90_NOWRITE, &
                        NF90_CLOSE, NF90_INQ_DIMID
      use netcdf_error, only: check_err
      use mod_nws13_data, only: null_flag_value
      use mod_datetime, only: t_timedelta, operator(+)
      implicit none

      real(8), parameter :: EPS = epsilon(1d0)

      class(t_nws13), intent(inout) :: self
      type(t_datetime), intent(in) :: CurrentTime
      real(8), intent(inout) :: EyeLonR(3), EyeLatR(3)
      logical, intent(inout) :: FoundEye
      integer, intent(inout) :: PowellGroupTemp
      integer :: ig, ig_sorted
      integer :: nc_id

      ! First, initialize the ADCIRC wind fields with flag values.
      ! These will be replaced with actual values or defaults after processing all groups.
      self%WX2 = null_flag_value
      self%WY2 = null_flag_value
      self%PR2 = null_flag_value
      self%previous_interp_time = self%next_interp_time

      ! Open a connection to the wind file.
      call check_err(NF90_OPEN(trim(adjustl(self%filename)), NF90_NOWRITE, NC_ID))

      ! Update active groups and get next time from active groups
      call self%update_active_groups(CurrentTime)

      self%next_interp_time = self%previous_interp_time + t_timedelta(seconds=600)
      PowellGroupTemp = self%current_powell_group

      ! Increment the storm centers.
      EyeLonR(1) = EyeLonR(2)
      EyeLatR(1) = EyeLatR(2)
      EyeLonR(2) = EyeLonR(3)
      EyeLatR(2) = EyeLatR(3)
      EyeLonR(3) = 0.d0
      EyeLatR(3) = 0.d0
      FoundEye = .false.

      ! Process each wind group in rank order (only enabled groups)
      do IG_sorted = 1, self%n_groups
         ig = self%sorted_group_indices(ig_sorted)
         if (.not. self%wind_groups(ig)%enabled) cycle
         call self%wind_groups(ig)%process(CurrentTime, self%WX2, self%WY2, self%PR2, &
                                           self%requested_powell_group, PowellGroupTemp, EyeLonR, &
                                           EyeLatR, FoundEye, self%wind_multiplier)
      end do

      ! Replace any remaining flag values with default background values
      where (abs(self%PR2 - null_flag_value) <= eps)
         self%PR2 = PRBCKGRND*100d0/RHOWAT0/G
         self%WX2 = 0.d0
         self%WY2 = 0.d0
      end where

      ! Close the file.
      call check_err(NF90_CLOSE(NC_ID))

   end subroutine t_nws13_update_wind_data

!-----------------------------------------------------------------------
!> @brief Reads multi-grid wind and pressure fields from the NetCDF file and
!>        interpolates them to the ADCIRC mesh
!>
!> @details Main interface for obtaining wind forcing at each time step. Manages
!> temporal interpolation between snapshots and updates wind data when needed.
!> Handles storm center tracking for Powell wind drag scheme, including logic
!> to reset tracking when switching between groups where storm centers may
!> jump large distances or if there is more than one vortex in the basin.
!>
!> @param[inout] self     NWS13 object containing wind field data
!> @param[in]    TimeLoc  Model time in seconds since cold start
!> @param[out]   WXOUT    Wind velocity x-component at mesh nodes (m/s)
!> @param[out]   WYOUT    Wind velocity y-component at mesh nodes (m/s)
!> @param[out]   PROUT    Atmospheric pressure at mesh nodes (m of water)
!> @param[inout] EyeLonR  Storm center longitude history array (degrees)
!> @param[inout] EyeLatR  Storm center latitude history array (degrees)
!> @param[inout] FoundEye Boolean flag indicating if storm center was found
!-----------------------------------------------------------------------
   subroutine t_nws13_get_wind_forcing(self, TimeLoc, WXOUT, WYOUT, PROUT, EyeLonR, EyeLatR, FoundEye)
!-----------------------------------------------------------------------
      use MESH, only: NP
      use mod_datetime, only: t_timedelta, null_datetime, operator(+), operator(-), operator(<), operator(==), &
                              operator(>=), operator(<=), operator(/=)
      use global, only: screenMessage, setMessageSource, unsetMessageSource
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      use global, only: allMessage, DEBUG
#endif
      implicit none
      class(t_nws13), intent(inout) :: self
      real(8), intent(in)    :: TimeLoc
      real(8), intent(inout) :: EyeLatR(3)
      real(8), intent(inout) :: EyeLonR(3)
      logical, intent(inout) :: FoundEye
      real(8), intent(out) :: PROUT(NP), WXOUT(NP), WYOUT(NP)

      integer :: PowellGroupTemp

      type(t_datetime)  :: CurrentTime

      CurrentTime = self%coldstart + t_timedelta(milliseconds=int(TimeLoc*1000d0))
      if (.not. self%check_data_started(CurrentTime)) then
         return
      elseif (CurrentTime < self%next_interp_time) then
         call self%time_interpolate_wind(CurrentTime, WXOUT, WYOUT, PROUT)
      else
         call self%update_wind_data(CurrentTime, EyeLonR, EyeLatR, FoundEye, PowellGroupTemp)
         call self%time_interpolate_wind(CurrentTime, WXOUT, WYOUT, PROUT)

         if (PowellGroupTemp /= self%current_powell_group) then
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
            self%current_powell_group = PowellGroupTemp
         end if

      end if

!-----------------------------------------------------------------------
   end subroutine t_nws13_get_wind_forcing
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Sorts group indices by their rank values (descending order)
!>
!> @details Implements bubble sort to order wind groups by priority rank.
!> Higher rank values indicate higher priority for overlapping grids.
!> The sorted indices are used to process groups in the correct order.
!>
!> @param[in] groups Array of t_nws13_group structures containing rank information
!> @return    indices Array of group indices sorted by rank (highest to lowest)
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

!-----------------------------------------------------------------------
!> @brief Applies group selection logic based on namelist settings
!>
!> @details Sets the enabled flag for a group based on user-specified
!> enabled groups list. If no list is provided, all groups are enabled.
!> If a list is provided, only groups in the list are enabled.
!>
!> @param[inout] group NWS13 group to apply selection logic to
!-----------------------------------------------------------------------
   subroutine apply_group_selection(group)
      implicit none
      class(t_nws13_group), intent(inout) :: group
      integer :: i

      ! Default: all groups enabled if no selection list provided
      if (.not. allocated(user_input%enabled_groups)) then
         group%enabled = .true.
         return
      end if

      ! Check if this group is in the enabled list
      group%enabled = .false.
      do i = 1, size(user_input%enabled_groups)
         if (trim(group%group_name) == trim(user_input%enabled_groups(i))) then
            group%enabled = .true.
            exit
         end if
      end do
   end subroutine apply_group_selection

!-----------------------------------------------------------------------
!> @brief Validates group selection and issues warnings for potential issues
!>
!> @details Checks for common configuration problems including unrecognized
!> group names, all groups disabled, and case sensitivity issues. Issues
!> appropriate warnings to help users debug their configuration.
!>
!> @param[in] groups Array of NWS13 groups to validate
!-----------------------------------------------------------------------
   subroutine validate_group_selection(groups)
      use global, only: logMessage, WARNING, INFO, scratchMessage
      implicit none
      class(t_nws13_group), intent(in) :: groups(:)
      integer :: i, j, enabled_count
      logical :: found_match

      ! If no group selection specified, all groups are enabled - no validation needed
      if (.not. allocated(user_input%enabled_groups)) return

      enabled_count = count(groups(:)%enabled)

      ! Check for no enabled groups
      if (enabled_count == 0) then
         call logMessage(WARNING, "NWS13: No wind groups are enabled! Check group names in namelist.")
         call logMessage(INFO, "Available groups in file:")
         do i = 1, size(groups)
            call logMessage(INFO, "  - '"//trim(groups(i)%group_name)//"'")
         end do
      end if

      ! Check for unrecognized group names in the enabled list
      do i = 1, size(user_input%enabled_groups)
         found_match = .false.
         do j = 1, size(groups)
            if (trim(user_input%enabled_groups(i)) == trim(groups(j)%group_name)) then
               found_match = .true.
               exit
            end if
         end do

         if (.not. found_match) then
            write (scratchMessage, '(A)') "NWS13 WARNING: Group '"// &
               trim(user_input%enabled_groups(i))//"' not found in wind file"
            call logMessage(WARNING, scratchMessage)
         end if
      end do

      ! Summary message
      write (scratchMessage, '(A,I0,A,I0,A)') "NWS13: ", enabled_count, " of ", &
         size(groups), " groups enabled"
      call logMessage(INFO, scratchMessage)
   end subroutine validate_group_selection

end module mod_nws13
