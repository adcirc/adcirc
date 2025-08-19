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
!> @author Zach Cobell, The Water Institute, zcobell@gmail.com
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

   implicit none

   !> @brief Flag value indicating missing or invalid data
   real(8), parameter  :: null_flag_value = -99999999.0d0

   !> @brief Container for meteorological data from NetCDF wind grid
   type t_windData
      real(8), allocatable :: lon(:, :) !< Grid longitude coordinates (degrees)
      real(8), allocatable :: lat(:, :) !< Grid latitude coordinates (degrees)
      real(8), allocatable :: u(:, :) !< Wind velocity, x-direction (m/s)
      real(8), allocatable :: v(:, :) !< Wind velocity, y-direction (m/s)
      real(8), allocatable :: p(:, :) !< Atmospheric pressure (converted to m water)
      real(8), allocatable :: ws(:, :) !< Wind speed magnitude (m/s)
   contains
      final :: wind_data_destructor
   end type t_windData

   !> @brief Stores bilinear interpolation data for grid-to-mesh mapping
   type t_interpolationData
      integer, allocatable :: ilon(:) !< Grid x-indices for each mesh node
      integer, allocatable :: ilat(:) !< Grid y-indices for each mesh node
      real(8), allocatable :: weights(:, :) !< Bilinear interpolation weights (NP, 4)
      logical :: initialized = .false. !< Flag indicating if interpolation is computed
      logical :: is_moving = .false. !< Flag indicating if grid moves with time
      real(8) :: corner_lon(4) = 0.0d0 !< Grid corner longitudes [SW, SE, NE, NW]
      real(8) :: corner_lat(4) = 0.0d0 !< Grid corner latitudes [SW, SE, NE, NW]
   end type t_interpolationData

   !> @brief Stores storm center data for a group to avoid repeated NetCDF reads
   type t_stormCenterData
      logical :: available = .false. !< Flag indicating if storm center data exists in file
      logical :: loaded = .false. !< Flag indicating if data has been loaded into memory
      real(8), allocatable :: clon(:) !< Storm center longitude for each time snapshot
      real(8), allocatable :: clat(:) !< Storm center latitude for each time snapshot
   end type t_stormCenterData

   !> @brief NetCDF wind group metadata and temporal information
   type t_nws13
      integer :: InclSnap !< Flag indicating if group is active for current time step
      integer :: NextSnap !< Index of next time snapshot in group
      integer :: NumSnap !< Total number of time snapshots in group
      integer :: PrevSnap !< Index of previous time snapshot in group
      integer :: rank !< Priority rank for overlapping grids
      integer :: group_id !< NetCDF group identifier
      integer :: NumLon !< Grid x-dimension size
      integer :: NumLat !< Grid y-dimension size
      character(len=100) :: group_name !< NetCDF group name from file
      type(t_datetime) :: NextTime !< Date/time of next snapshot
      type(t_datetime) :: PrevTime !< Date/time of previous snapshot
      type(t_datetime), allocatable :: SnapTimes(:) !< Array of all snapshot times
      type(t_interpolationData) :: interp !< Grid-to-mesh interpolation data
      type(t_stormCenterData) :: storm_center !< Cached storm center data
   end type t_nws13

   interface t_nws13
      module procedure :: nws13_constructor
   end interface

   interface t_stormCenterData
      module procedure :: storm_center_constructor
   end interface

   interface t_windData
      module procedure :: wind_data_constructor
   end interface

   !> @brief Path to the NetCDF wind input file (default: fort.22.nc)
   character(LEN=1000) :: NWS13File = "fort.22.nc"

   !> @brief Group number to use for Powell wind drag calculation (0 = auto-select)
   integer :: NWS13GroupForPowell = 0

   !> @brief Currently active group for Powell wind drag scheme
   integer :: PowellGroup = 0

   !> @brief Temporary group for Powell wind drag scheme validation
   integer :: PowellGroupTemp = 0

   !> @brief Multiplier applied to wind velocities for sensitivity analysis
   real(8) :: NWS13WindMultiplier = 1.0d0

   !> @brief Reference time for cold start initialization from namelist
   type(t_datetime) :: NWS13ColdStart

   !> @brief Array of NWS13 wind groups from NetCDF file
   type(t_nws13), allocatable :: NWS13(:)

   !> @brief Group indices sorted by priority rank
   integer, allocatable :: sorted_group_indices(:)

   private

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
!> @brief Constructor for t_windData structure with optional lon/lat allocation
!> @param[in] ni           Grid x-dimension size
!> @param[in] nj           Grid y-dimension size
!> @param[in] use_lat_lon  Optional flag to allocate lon/lat arrays
!> @return    wind_data    Initialized wind data structure
!-----------------------------------------------------------------------
   pure function wind_data_constructor(ni, nj, use_lat_lon) result(wind_data)
      implicit none
      integer, intent(in) :: ni, nj
      logical, intent(in), optional :: use_lat_lon
      type(t_windData) :: wind_data

      if (present(use_lat_lon)) then
         if (use_lat_lon) then
            allocate (wind_data%lon(ni, nj), wind_data%lat(ni, nj))
         end if
      end if

      allocate (wind_data%u(ni, nj), wind_data%v(ni, nj), wind_data%p(ni, nj), wind_data%ws(ni, nj))

   end function wind_data_constructor

!-----------------------------------------------------------------------
!> @brief Destructor for t_windData structure, deallocates all arrays
!> @param[inout] this  Wind data structure to deallocate
!-----------------------------------------------------------------------
   pure subroutine wind_data_destructor(this)
      type(t_windData), intent(inout) :: this
      if (allocated(this%lon)) deallocate (this%lon, this%lat)
      if (allocated(this%u)) deallocate (this%u, this%v, this%p, this%ws)
   end subroutine wind_data_destructor

!-----------------------------------------------------------------------
!> @brief Constructor for t_stormCenterData structure that reads storm center data if available
!> @param[in] nc_idg NetCDF group ID to read from
!> @param[in] num_snaps Number of time snapshots
!> @return storm_center Initialized storm center data structure with cached data
!-----------------------------------------------------------------------
   function storm_center_constructor(nc_idg, num_snaps) result(storm_center)
      use netcdf, only: NF90_INQ_VARID, NF90_GET_VAR, NF90_NOERR
      use netcdf_error, only: check_err
      implicit none
      integer, intent(in) :: nc_idg, num_snaps
      type(t_stormCenterData) :: storm_center
      integer :: nc_var, nc_err

      ! Allocate arrays
      allocate (storm_center%clon(num_snaps))
      allocate (storm_center%clat(num_snaps))

      ! Check if storm center longitude data exists
      nc_err = NF90_INQ_VARID(nc_idg, "clon", nc_var)
      if (nc_err == NF90_NOERR) then
         storm_center%available = .true.
         call check_err(NF90_GET_VAR(nc_idg, nc_var, storm_center%clon(:), &
                                     START=[1], COUNT=[num_snaps]))

         ! Check if storm center latitude data exists
         nc_err = NF90_INQ_VARID(nc_idg, "clat", nc_var)
         if (nc_err == NF90_NOERR) then
            call check_err(NF90_GET_VAR(nc_idg, nc_var, storm_center%clat(:), &
                                        START=[1], COUNT=[num_snaps]))
            storm_center%loaded = .true.
         else
            storm_center%available = .false.
         end if
      else
         storm_center%available = .false.
      end if

   end function storm_center_constructor

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
         NWS13(IG) = t_nws13(group_ids(IG))
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
!> @brief Constructs a new t_nws13 object from a NetCDF group
!> @param[in]    nc_idg  NetCDF group ID to initialize from
!> @return       grp     Initialized t_nws13 structure
!-----------------------------------------------------------------------
   type(t_nws13) function nws13_constructor(NC_IDG) result(grp)

      use netcdf, only: nf90_inq_dimid, nf90_inq_varid, nf90_inquire_dimension, &
                        nf90_get_var, nf90_inquire_variable, NF90_GLOBAL, nf90_get_att, &
                        nf90_noerr, nf90_inq_grpname
      use netcdf_error, only: check_err
      use global, only: ERROR, allMessage
      use mod_datetime, only: t_timedelta, operator(+)
#ifdef CMPI
      use messenger, only: msg_fini
#endif

      implicit none

      integer, intent(in) :: NC_IDG
      integer :: NC_DIM
      integer :: NC_VAR
      integer :: NC_ERR
      integer :: ndims
      integer :: IGS
      integer, allocatable :: dimids(:)
      integer, allocatable :: TempI(:)
      type(t_datetime) :: WindRefDatetime

      grp%group_id = NC_IDG
      call check_err(NF90_INQ_GRPNAME(NC_IDG, grp%group_name))
      call check_err(NF90_INQ_DIMID(NC_IDG, "time", NC_DIM))
      call check_err(NF90_INQUIRE_DIMENSION(NC_IDG, NC_DIM, LEN=grp%NumSnap))
      call check_err(NF90_INQ_VARID(NC_IDG, "time", NC_VAR))

      ! Read and store grid dimensions
      call check_err(NF90_INQ_DIMID(NC_IDG, "xi", NC_DIM))
      call check_err(NF90_INQUIRE_DIMENSION(NC_IDG, NC_DIM, LEN=grp%NumLon))
      call check_err(NF90_INQ_DIMID(NC_IDG, "yi", NC_DIM))
      call check_err(NF90_INQUIRE_DIMENSION(NC_IDG, NC_DIM, LEN=grp%NumLat))

      WindRefDatetime = read_dataset_reference_time(NC_IDG)

      allocate (grp%SnapTimes(1:grp%NumSnap))
      allocate (TempI(1:grp%NumSnap))
      ! Pull the time for the first snap ...
      call check_err(NF90_GET_VAR(NC_IDG, NC_VAR, TempI(:), START=[1], COUNT=[grp%NumSnap]))
      ! ... and load it into our array.
      do IGS = 1, grp%NumSnap
         grp%SnapTimes(IGS) = WindRefDatetime + t_timedelta(minutes=TempI(IGS))
      end do
      deallocate (TempI)

      grp%PrevTime = grp%SnapTimes(1)
      grp%PrevSnap = 1
      ! Pull the time for the second snap ...
      ! ... and also load it into our array.
      grp%NextTime = grp%SnapTimes(2)
      grp%NextSnap = 2

      ! Read the rank attribute for this group (priority for overlapping grids)
      NC_ERR = NF90_GET_ATT(NC_IDG, NF90_GLOBAL, "rank", grp%rank)
      if (NC_ERR /= NF90_NOERR) then
         call allMessage(ERROR, "Rank attribute not found in nws13 group")
#ifdef CMPI
         call msg_fini()
#endif
         call exit(1)
      end if

      ! Check if this grid moves by examining the dimensionality of lon/lat
      ! Moving grids have 3D arrays (time, yi, xi), stationary have 2D (yi, xi)
      call check_err(NF90_INQ_VARID(NC_IDG, "lon", NC_VAR))
      allocate (dimids(10))
      NC_ERR = NF90_INQUIRE_VARIABLE(NC_IDG, NC_VAR, ndims=ndims, dimids=dimids)
      deallocate (dimids)

      if (ndims == 3) then
         grp%interp%is_moving = .true.
         grp%interp%initialized = .false.
      elseif (ndims == 2) then
         grp%interp = read_2d_grid(NC_IDG)
      else
         call allMessage(ERROR, "Invalid grid dimension for NWS13")
#ifdef CMPI
         call msg_fini()
#endif
         call exit(1)
      end if

      ! Initialize storm center data cache - read all storm center data during startup
      grp%storm_center = t_stormCenterData(NC_IDG, grp%NumSnap)

   end function nws13_constructor

!-----------------------------------------------------------------------
!> @brief Reads the reference date/time from a NetCDF dataset time variable
!> @param[in]    nc_idg  NetCDF group ID
!> @return       reftime Reference date/time as t_datetime object
!-----------------------------------------------------------------------
   type(t_datetime) function read_dataset_reference_time(NC_IDG) result(reftime)
      use netcdf, only: nf90_inq_varid, nf90_get_att
      use netcdf_error, only: check_err
      use global, only: ERROR, allMessage

#ifdef CMPI
      use messenger, only: msg_fini
#endif

      implicit none

      integer, intent(in) :: NC_IDG
      integer :: nc_var
      character(LEN=100) :: TimeUnits
      character(LEN=100) :: TimeString

      ! Find the starting date/time in YYYY-MM-DDTHH:MM:SS format
      call check_err(NF90_INQ_VARID(NC_IDG, "time", NC_VAR))
      call check_err(NF90_GET_ATT(NC_IDG, NC_VAR, "units", TimeUnits))
      TimeString = trim(adjustl(TimeUnits))
      TimeString = TimeUnits(15:len_trim(TimeUnits))

      TimeUnits = TimeUnits(1:14)
      if (TimeUnits /= "minutes since ") then
         call allMessage(ERROR, "NWS13: Invalid time units in NWS13INIT.")
#ifdef CMPI
         call msg_fini()
#endif
         call exit(1)
      end if

      reftime = t_datetime(TimeString, ["%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S"])
      if (.not. reftime%valid()) then
         call allMessage(ERROR, "NWS13: Invalid reference time in NWS13INIT: '"//trim(TimeString)//"'")
#ifdef CMPI
         call msg_fini()
#endif
         call exit(1)
      end if

   end function read_dataset_reference_time

!-----------------------------------------------------------------------
!> @brief Reads static 2D grid coordinates and computes interpolation weights
!> @param[in]    nc_idg      NetCDF group ID
!> @return       interp_data Structure containing interpolation indices and weights
!-----------------------------------------------------------------------
   type(t_interpolationData) function read_2d_grid(NC_IDG) result(interp_data)
      use netcdf, only: NF90_INQ_DIMID, NF90_INQUIRE_DIMENSION, NF90_INQ_VARID, &
                        NF90_GET_VAR
      use netcdf_error, only: check_err

      implicit none

      integer, intent(in) :: NC_IDG
      integer :: NC_DIM
      integer :: NC_VAR
      integer :: NumLon
      integer :: Numlat
      real(8), allocatable :: Lon(:, :)
      real(8), allocatable :: Lat(:, :)

      call check_err(NF90_INQ_DIMID(NC_IDG, "xi", NC_DIM))
      call check_err(NF90_INQUIRE_DIMENSION(NC_IDG, NC_DIM, LEN=NumLon))
      call check_err(NF90_INQ_DIMID(NC_IDG, "yi", NC_DIM))
      call check_err(NF90_INQUIRE_DIMENSION(NC_IDG, NC_DIM, LEN=NumLat))

      allocate (Lon(1:NumLon, 1:NumLat))
      allocate (Lat(1:NumLon, 1:NumLat))

      call check_err(NF90_INQ_VARID(NC_IDG, "lon", NC_VAR))
      call check_err(NF90_GET_VAR(NC_IDG, NC_VAR, Lon(:, :), &
                                  START=[1, 1], COUNT=[NumLon, NumLat]))
      call check_err(NF90_INQ_VARID(NC_IDG, "lat", NC_VAR))
      call check_err(NF90_GET_VAR(NC_IDG, NC_VAR, Lat(:, :), &
                                  START=[1, 1], COUNT=[NumLon, NumLat]))

      interp_data = NWS13INTERP(NumLon, NumLat, Lon, Lat)
      interp_data%is_moving = .false.

      deallocate (Lon, Lat)

   end function read_2d_grid

!-----------------------------------------------------------------------
!> @brief Updates interpolation data if grid has moved. The grid is assumed
!>        not to have moved if the corners are static
!> @param[inout] interp_data  Interpolation structure to update
!> @param[in]    numlon       X-dimension size of grid
!> @param[in]    numlat       Y-dimension size of grid
!> @param[in]    lon          Grid longitudes (2D)
!> @param[in]    lat          Grid latitudes (2D)
!-----------------------------------------------------------------------
   subroutine update_moving_grid_interpolation(interp_data, numlon, numlat, lon, lat)
      implicit none
      type(t_interpolationData), intent(inout) :: interp_data
      integer, intent(in) :: numlon, numlat
      real(8), intent(in) :: lon(numlon, numlat), lat(numlon, numlat)

      real(8) :: new_corners_lon(4), new_corners_lat(4)
      real(8), parameter :: eps = epsilon(1d0)
      logical :: grid_moved

      ! Extract corner coordinates: [SW, SE, NE, NW]
      new_corners_lon(1) = lon(1, 1) ! SW
      new_corners_lon(2) = lon(numlon, 1) ! SE
      new_corners_lon(3) = lon(numlon, numlat) ! NE
      new_corners_lon(4) = lon(1, numlat) ! NW

      new_corners_lat(1) = lat(1, 1) ! SW
      new_corners_lat(2) = lat(numlon, 1) ! SE
      new_corners_lat(3) = lat(numlon, numlat) ! NE
      new_corners_lat(4) = lat(1, numlat) ! NW

      ! Check if any corner has moved significantly
      grid_moved = .false.
      if (interp_data%initialized) then
         grid_moved = any(abs(new_corners_lon - interp_data%corner_lon) > eps) .or. &
                      any(abs(new_corners_lat - interp_data%corner_lat) > eps)
      end if

      ! Update interpolation if grid moved or not initialized
      if (.not. interp_data%initialized .or. grid_moved) then
         interp_data = NWS13INTERP(numlon, numlat, lon, lat)
         interp_data%is_moving = .true.
         interp_data%corner_lon = new_corners_lon
         interp_data%corner_lat = new_corners_lat
      end if

   end subroutine update_moving_grid_interpolation

!-----------------------------------------------------------------------
!> @brief Reads lon/lat coordinates for moving grids and updates interpolation
!> @param[in]    nc_idg    NetCDF group ID
!> @param[in]    currsnap  Current time snapshot index
!> @param[inout] grp       NWS13 group structure to update
!> @param[inout] data      Wind data structure to populate with grid coordinates
!-----------------------------------------------------------------------
   subroutine read_moving_grid(nc_idg, currsnap, grp, data)
      use netcdf, only: NF90_INQ_VARID, NF90_GET_VAR
      use netcdf_error, only: check_err

      implicit none

      integer, intent(in) :: nc_idg, currsnap
      type(t_nws13), intent(inout) :: grp
      type(t_windData), intent(inout) :: data

      integer :: nc_var

      ! Read lon/lat coordinates for this time snapshot
      call check_err(NF90_INQ_VARID(nc_idg, "lon", nc_var))
      call check_err(NF90_GET_VAR(nc_idg, nc_var, data%lon(:, :), &
                                  START=[1, 1, currsnap], &
                                  COUNT=[grp%NumLon, grp%NumLat, 1]))
      call check_err(NF90_INQ_VARID(nc_idg, "lat", nc_var))
      call check_err(NF90_GET_VAR(nc_idg, nc_var, data%lat(:, :), &
                                  START=[1, 1, currsnap], &
                                  COUNT=[grp%NumLon, grp%NumLat, 1]))
      call update_moving_grid_interpolation(grp%interp, grp%NumLon, &
                                            grp%NumLat, data%lon, data%lat)
   end subroutine read_moving_grid

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
                              operator(>=), operator(<=)
      use global, only: screenMessage, setMessageSource, unsetMessageSource, &
                        INFO, RNDAY
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      use global, only: allMessage, DEBUG
#endif
      implicit none

      real(8), parameter :: EPS = epsilon(1d0)

      real(8), intent(inout) :: EyeLatR(3)
      real(8), intent(inout) :: EyeLonR(3)
      logical, intent(out)   :: FoundEye
      real(8), intent(out)   :: PRN2(NP)
      real(8), intent(in)    :: TimeLoc
      real(8), intent(inout) :: WVNX2(NP)
      real(8), intent(inout) :: WVNY2(NP)
      real(8), intent(out)   :: WTIME2

      character(LEN=200) :: Line

      integer :: CurrSnap
      integer :: IG, IG_sorted
      integer :: IS
      integer :: IV
      integer :: NC_ID
      integer :: NC_IDG

      type(t_windData) :: CurrentData, PrevData, NextData

      real(8) :: CLat
      real(8) :: CLon
      real(8) :: NextCLat
      real(8) :: NextCLon
      real(8) :: PP(NP)
      real(8) :: PrevCLat
      real(8) :: PrevCLon
      real(8) :: TimeInterpFactor
      real(8) :: UU(NP)
      real(8) :: VV(NP)
      real(8) :: Wind(NP)
      real(8) :: sWdir

      real(8) :: ws_this(4)
      real(8) :: p_this(4)

      type(t_interpolationData) :: current_interp

      type(t_datetime)  :: CurrentTime
      type(t_datetime)  :: wtime2_dt
      type(t_timedelta) :: dt
      type(t_timedelta) :: TimeInterpNumerator, TimeInterpDenominator

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

      PrevCLon = 0d0
      PrevCLat = 0d0
      NextCLon = 0d0
      NextCLat = 0d0

      ! First, initialize the ADCIRC wind fields with flag values.
      ! These will be replaced with actual values or defaults after processing all groups.
      WVNX2 = null_flag_value
      WVNY2 = null_flag_value
      PRN2 = null_flag_value

      ! Open a connection to the wind file.
      call check_err(NF90_OPEN(trim(adjustl(NWS13File)), NF90_NOWRITE, NC_ID))

      ! We need to determine how many groups are active at the current
      ! time in the simulation.  Loop over the groups and find which
      ! groups to include.
      NWS13(:)%InclSnap = 0
      WTIME2 = -1.d0
      WTIME2_dt = null_datetime()
      CurrentTime = NWS13Coldstart + t_timedelta(milliseconds=int(TimeLoc*1000d0))

      do IG = 1, size(NWS13)
         if (CurrentTime < NWS13(IG)%PrevTime) then
            ! Then we still haven't reached the start of this group,
            ! so we at least need to consider this as the end of the current
            ! interval.  If we don't find any other groups with a closer
            ! time, then this current interval will be zero winds.
            if (WTIME2_dt == null_datetime() .or. NWS13(IG)%PrevTime < WTIME2_dt) then
               WTIME2_dt = NWS13(IG)%PrevTime
               dt = NWS13(IG)%PrevTime - NWS13ColdStart
               WTIME2 = dble(dt%total_milliseconds())/1000.0d0
            end if
            cycle
         end if

         ! Check whether we have progressed past the end of the current
         ! wind snap, and if so, then shift to the next wind snap.
         do while (CurrentTime >= NWS13(IG)%NextTime)
            if (NWS13(IG)%NextSnap == NWS13(IG)%NumSnap) then
               ! Then we have reached the end of the wind snaps for this
               ! group, and it will be excluded below.
               exit
            end if
            NWS13(IG)%PrevSnap = NWS13(IG)%NextSnap
            NWS13(IG)%PrevTime = NWS13(IG)%NextTime
            NWS13(IG)%NextSnap = NWS13(IG)%NextSnap + 1
            NWS13(IG)%NextTime = NWS13(IG)%SnapTimes(NWS13(IG)%NextSnap)
         end do
         ! Then check whether we are within the current wind snap.
         if (NWS13(IG)%PrevTime <= CurrentTime .and. CurrentTime <= NWS13(IG)%NextTime) then
            NWS13(IG)%InclSnap = 1
            ! Update our end time if this is the first snap we found ...
            ! ... or if it is earlier than what we have found already.
            if (WTIME2_dt == null_datetime() .or. (NWS13(IG)%NextTime < WTIME2_dt)) then
               dt = NWS13(IG)%NextTime - NWS13ColdStart
               WTIME2_dt = NWS13(IG)%NextTime
               WTIME2 = dble(dt%total_milliseconds())/1000.0d0
            end if
         end if
      end do
      ! The only way we would not find a next snap is for the scenario
      ! when the simulation has progressed past the end of the wind
      ! file.  In that case, we can use the simulation length as the end
      ! of this next snap.
      if (WTIME2_dt == null_datetime()) then
         WTIME2 = RNDAY*86400.d0
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

         ! Loop over the groups in pre-computed rank order
         do IG_sorted = 1, size(NWS13)
            IG = sorted_group_indices(IG_sorted)
            if (NWS13(IG)%InclSnap == 0) cycle

            write (Line, '(A,A,A,I0,2A)') "Processing group '", trim(NWS13(IG)%group_name), &
               "' with rank ", NWS13(IG)%rank, " for wind/pressure interpolation with time ", &
               trim(NWS13(IG)%NextTime%to_iso_string())
            call screenMessage(INFO, Line)

            ! Connect to this group.
            NC_IDG = NWS13(IG)%group_id

            ! Create wind data objects for this group's dimensions
            CurrentData = t_windData(NWS13(IG)%NumLon, NWS13(IG)%NumLat, NWS13(IG)%interp%is_moving)
            PrevData = t_windData(NWS13(IG)%NumLon, NWS13(IG)%NumLat, NWS13(IG)%interp%is_moving)
            NextData = t_windData(NWS13(IG)%NumLon, NWS13(IG)%NumLat, NWS13(IG)%interp%is_moving)

            ! Loop over the snaps.
            do IS = 1, 2
               if (IS == 1) then
                  CurrSnap = NWS13(IG)%PrevSnap
               elseif (IS == 2) then
                  CurrSnap = NWS13(IG)%NextSnap
               end if

               ! 1. Read this snap from the wind file if it is a moving grid
               if (NWS13(IG)%interp%is_moving) then
                  call read_moving_grid(NC_IDG, CurrSnap, NWS13(IG), CurrentData)
               end if

               ! Read winds and pressures.
               call read_netcdf_var_with_attrs(NC_IDG, "U10", CurrSnap, NWS13(IG)%NumLon, NWS13(IG)%NumLat, CurrentData%U)
               call read_netcdf_var_with_attrs(NC_IDG, "V10", CurrSnap, NWS13(IG)%NumLon, NWS13(IG)%NumLat, CurrentData%V)
               call read_netcdf_var_with_attrs(NC_IDG, "PSFC", CurrSnap, NWS13(IG)%NumLon, NWS13(IG)%NumLat, CurrentData%P)

               ! Adjust the wind velocity.
               CurrentData%U = CurrentData%U*NWS13WindMultiplier
               CurrentData%V = CurrentData%V*NWS13WindMultiplier

               ! Convert the pressures into meters of water.
               CurrentData%P = CurrentData%P*100.d0/RHOWAT0/G
               FoundEye = NWS13(IG)%storm_center%available .and. NWS13(IG)%storm_center%loaded

               ! Save the snaps.
               if (IS == 1) then
                  if (NWS13(IG)%interp%is_moving) then
                     PrevData%Lon = CurrentData%Lon
                     PrevData%Lat = CurrentData%Lat
                  end if
                  PrevData%U = CurrentData%U
                  PrevData%V = CurrentData%V
                  PrevData%P = CurrentData%P
                  if (FoundEye) then
                     PrevCLon = NWS13(IG)%storm_center%clon(CurrSnap)
                     PrevCLat = NWS13(IG)%storm_center%clat(CurrSnap)
                  end if
               elseif (IS == 2) then
                  if (NWS13(IG)%interp%is_moving) then
                     NextData%Lon = CurrentData%Lon
                     NextData%Lat = CurrentData%Lat
                  end if
                  NextData%U = CurrentData%U
                  NextData%V = CurrentData%V
                  NextData%P = CurrentData%P
                  if (FoundEye) then
                     NextCLon = NWS13(IG)%storm_center%clon(CurrSnap)
                     NextCLat = NWS13(IG)%storm_center%clat(CurrSnap)
                  end if
               end if
            end do

            ! 2. Interpolate the wind snaps to start/end of current interval.
            TimeInterpNumerator = CurrentTime - NWS13(IG)%PrevTime
            TimeInterpDenominator = NWS13(IG)%NextTime - NWS13(IG)%PrevTime
            TimeInterpFactor = dble(TimeInterpNumerator%total_milliseconds())/dble(TimeInterpDenominator%total_milliseconds())

            if (NWS13(IG)%interp%is_moving) then
               CurrentData%Lon = PrevData%Lon + (NextData%Lon - PrevData%Lon)*TimeInterpFactor
               CurrentData%Lat = PrevData%Lat + (NextData%Lat - PrevData%Lat)*TimeInterpFactor
            end if
            CurrentData%U = PrevData%U + (NextData%U - PrevData%U)*TimeInterpFactor
            CurrentData%V = PrevData%V + (NextData%V - PrevData%V)*TimeInterpFactor
            CurrentData%P = PrevData%P + (NextData%P - PrevData%P)*TimeInterpFactor

            ! Time Interpolate scalar wind
            PrevData%Ws = sqrt(PrevData%U**2d0 + PrevData%V**2d0)
            NextData%Ws = sqrt(NextData%U**2d0 + NextData%V**2d0)
            CurrentData%Ws = PrevData%Ws + (NextData%Ws - PrevData%Ws)*TimeInterpFactor
            if (FoundEye) then
               CLon = PrevCLon + (NextCLon - PrevCLon)*TimeInterpFactor
               CLat = PrevCLat + (NextCLat - PrevCLat)*TimeInterpFactor
            end if

            ! 3. Use pre-computed interpolation data (updated when grid is read for moving grids)
            current_interp = NWS13(IG)%interp

            do IV = 1, NP
               ! Skip if this node was already assigned by a higher priority grid
               if (abs(PRN2(IV) - null_flag_value) > eps) cycle

               if (current_interp%ilon(IV) > 0) then

                  ws_this(1) = CurrentData%Ws(current_interp%ilon(IV), current_interp%ilat(IV))
                  ws_this(2) = CurrentData%Ws(current_interp%ilon(IV) + 1, current_interp%ilat(IV))
                  ws_this(3) = CurrentData%Ws(current_interp%ilon(IV) + 1, current_interp%ilat(IV) + 1)
                  ws_this(4) = CurrentData%Ws(current_interp%ilon(IV), current_interp%ilat(IV) + 1)
                  p_this(1) = CurrentData%P(current_interp%ilon(IV), current_interp%ilat(IV))
                  p_this(2) = CurrentData%P(current_interp%ilon(IV) + 1, current_interp%ilat(IV))
                  p_this(3) = CurrentData%P(current_interp%ilon(IV) + 1, current_interp%ilat(IV) + 1)
                  p_this(4) = CurrentData%P(current_interp%ilon(IV), current_interp%ilat(IV) + 1)

                  ! do not apply if we pick up a a flag value
                  if (any(abs(ws_this - null_flag_value) > eps) .or. &
                      any(abs(p_this - null_flag_value) > eps)) then
                     UU(IV) = interpolate_bilinear(CurrentData%U, current_interp%weights(IV, :), &
                                                   current_interp%ilon(IV), current_interp%ilat(IV))
                     VV(IV) = interpolate_bilinear(CurrentData%V, current_interp%weights(IV, :), &
                                                   current_interp%ilon(IV), current_interp%ilat(IV))
                     Wind(IV) = interpolate_bilinear(CurrentData%Ws, current_interp%weights(IV, :), &
                                                     current_interp%ilon(IV), current_interp%ilat(IV))
                     PP(IV) = interpolate_bilinear(CurrentData%P, current_interp%weights(IV, :), &
                                                   current_interp%ilon(IV), current_interp%ilat(IV))
                  end if

                  ! Adjust U/V Based on scalar wind magnitude
                  if (abs(UU(IV)) <= eps .and. abs(VV(IV)) <= eps) then
                     sWdir = 0d0
                  else
                     sWdir = atan2(UU(IV), VV(IV))
                     UU(IV) = Wind(IV)*sin(sWDir)
                     VV(IV) = Wind(IV)*cos(sWDir)
                  end if

                  PRN2(IV) = PP(IV)
                  WVNX2(IV) = UU(IV)
                  WVNY2(IV) = VV(IV)
               end if
            end do

            if (FoundEye) then
               if ((NWS13GroupForPowell == 0) .or. (NWS13GroupForPowell == IG)) then
                  ! Assign the storm center to the next spot in the array.
                  EyeLonR(3) = CLon
                  EyeLatR(3) = CLat
                  PowellGroupTemp = IG
               end if
            end if

         end do

      end if

      ! Replace any remaining flag values with default background values
      where (abs(PRN2 - null_flag_value) < eps)
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
!> @brief Computes interpolation weights for mapping grid data to mesh nodes
!> @param[in]    numlon x-dimension size of grid
!> @param[in]    numlat y-dimension size of grid
!> @param[in]    lon    grid longitudes (2D)
!> @param[in]    lat    grid latitudes (2D)
!> @return       interpData Structure containing interpolation indices and weights
!-----------------------------------------------------------------------
   type(t_interpolationData) function NWS13INTERP(NumLon, NumLat, Lon, Lat) result(interpData)
!-----------------------------------------------------------------------
      use ADC_CONSTANTS, only: RAD2DEG
      use kdtree2_module, only: KDTREE2, KDTREE2_RESULT, kdtree2_create, &
                                kdtree2_n_nearest, kdtree2_destroy
      use MESH, only: NP, SLAM, SFEA
      use GLOBAL, only: setMessageSource, unsetMessageSource
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: allMessage, DEBUG
#endif

      implicit none

      integer, intent(in)  :: NumLat
      integer, intent(in)  :: NumLon
      real(8), intent(in)  :: Lon(NumLon, NumLat)
      real(8), intent(in)  :: Lat(NumLon, NumLat)

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
      NumResult = min(NumResult, 20)

      allocate (ElemCenter(1:2, 1:(NumLon - 1)*(NumLat - 1)))
      Counter = 0
      do IO = 1, NumLon - 1
         do IA = 1, NumLat - 1
            Counter = Counter + 1
            ElemCenter(:, Counter) = quad_center_x_y( &
                                     Lon(IO, IA), Lat(IO, IA), &
                                     Lon(IO + 1, IA), Lat(IO + 1, IA), &
                                     Lon(IO + 1, IA + 1), Lat(IO + 1, IA + 1), &
                                     Lon(IO, IA + 1), Lat(IO, IA + 1))
         end do
      end do

      kdTree => kdtree2_create(ElemCenter, rearrange=.true., sort=.true.)
      allocate (kdResult(1:NumResult))

      ! Allocate output arrays
      allocate (interpData%ilon(1:NP))
      allocate (interpData%ilat(1:NP))
      allocate (interpData%weights(1:NP, 1:4))

      ! Initialize
      interpData%ilon = 0
      interpData%ilat = 0
      interpData%weights = 0.d0
      interpData%initialized = .true.

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

            W0 = compute_weight(Lon(IO, IA), Lat(IO, IA), Lon(IO + 1, IA), Lat(IO + 1, IA), &
                                Lon(IO + 1, IA + 1), Lat(IO + 1, IA + 1), Lon(IO, IA + 1), Lat(IO, IA + 1))
            W1 = compute_weight(AdcLon, AdcLat, Lon(IO + 1, IA), Lat(IO + 1, IA), &
                                Lon(IO + 1, IA + 1), Lat(IO + 1, IA + 1), Lon(IO, IA + 1), Lat(IO, IA + 1))
            W2 = compute_weight(Lon(IO, IA), Lat(IO, IA), AdcLon, AdcLat, &
                                Lon(IO + 1, IA + 1), Lat(IO + 1, IA + 1), Lon(IO, IA + 1), Lat(IO, IA + 1))
            W3 = compute_weight(Lon(IO, IA), Lat(IO, IA), Lon(IO + 1, IA), Lat(IO + 1, IA), &
                                AdcLon, AdcLat, Lon(IO, IA + 1), Lat(IO, IA + 1))
            W4 = compute_weight(Lon(IO, IA), Lat(IO, IA), Lon(IO + 1, IA), Lat(IO + 1, IA), &
                                Lon(IO + 1, IA + 1), Lat(IO + 1, IA + 1), AdcLon, AdcLat)

            if ((W1 + W2 + W3 + W4) <= (W0*2.001d0)) then
               ElemFound = .true.
               exit
            end if
            IE = IE + 1
         end do
         if (ElemFound) then
            interpData%ilon(IV) = IO
            interpData%ilat(IV) = IA
            interpData%weights(IV, 1) = W1/(W0*2.d0)
            interpData%weights(IV, 2) = W2/(W0*2.d0)
            interpData%weights(IV, 3) = W3/(W0*2.d0)
            interpData%weights(IV, 4) = W4/(W0*2.d0)
         end if

      end do

      call kdtree2_destroy(tp=kdTree)
      deallocate (kdResult)
      deallocate (ElemCenter)

#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()

   end function NWS13INTERP
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

!-----------------------------------------------------------------------
!> @brief Computes the weight for a point in a quadrilateral
!> @param[in] x1, y1, x2, y2, x3, y3, x4, y4 coordinates of the quadrilateral
!> @return w the weight of the point in the quadrilateral
!-----------------------------------------------------------------------
   real(8) pure function compute_weight(x1, y1, x2, y2, x3, y3, x4, y4) result(w)
      real(8), intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
      w = abs((x1*y2 - x2*y1) + (x2*y3 - x3*y2) + (x3*y4 - x4*y3) + (x4*y1 - x1*y4))
   end function compute_weight

!-----------------------------------------------------------------------
!> @brief Computes the center of a quadrilateral given its corner coordinates
!> @param[in] x1, y1, x2, y2, x3, y3, x4, y4 coordinates of the quadrilateral
!> @return center the coordinates of the center of the quadrilateral
!-----------------------------------------------------------------------
   pure function quad_center_x_y(x1, y1, x2, y2, x3, y3, x4, y4) result(center)
      real(8), intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
      real(8) :: center(2)
      center(1) = 0.25d0*(x1 + x2 + x3 + x4)
      center(2) = 0.25d0*(y1 + y2 + y3 + y4)
   end function quad_center_x_y

!-----------------------------------------------------------------------
!> @brief Bilinear interpolation of a 2D array using weights
!> @param[in] arr 2D array to interpolate
!> @param[in] w   weights for the bilinear interpolation
!> @param[in] i, j indices of the top-left corner of the interpolation
!> @return val the interpolated value
!-----------------------------------------------------------------------
   real(8) pure function interpolate_bilinear(arr, w, i, j) result(val)
      real(8), intent(in) :: arr(:, :)
      real(8), intent(in) :: w(4)
      integer, intent(in) :: i, j
      val = w(1)*arr(i, j) + w(2)*arr(i + 1, j) + w(3)*arr(i + 1, j + 1) + w(4)*arr(i, j + 1)
   end function interpolate_bilinear

!-----------------------------------------------------------------------
!> @brief Reads a NetCDF variable and applies fill value, scale, and offset
!> @param[in]  ncid      NetCDF group ID
!> @param[in]  varname   Name of the variable
!> @param[in]  currsnap  Time index to read
!> @param[in]  numlon    X-dimension size
!> @param[in]  numlat    Y-dimension size
!> @param[out] arr       Output array
!-----------------------------------------------------------------------
   subroutine read_netcdf_var_with_attrs(ncid, varname, currsnap, numlon, numlat, arr)
      use netcdf, only: NF90_INQ_VARID, NF90_GET_VAR, NF90_GET_ATT, NF90_NOERR
      use netcdf_error, only: check_err
      implicit none
      real(8), parameter :: eps = epsilon(1d0)
      integer, intent(in) :: ncid, currsnap, numlon, numlat
      character(len=*), intent(in) :: varname
      real(8), intent(out) :: arr(numlon, numlat)

      integer :: nc_var, nc_err
      real(8) :: var_scale_attr, var_offset_attr, var_fillval

      call check_err(NF90_INQ_VARID(ncid, varname, nc_var))
      call check_err(NF90_GET_VAR(ncid, nc_var, arr(:, :), START=[1, 1, currsnap], COUNT=[numlon, numlat, 1]))

      nc_err = NF90_GET_ATT(ncid, nc_var, "_FillValue", var_fillval)
      if (nc_err == NF90_NOERR) where (abs(arr - var_fillval) <= eps) arr = null_flag_value

      nc_err = NF90_GET_ATT(ncid, nc_var, "scale_factor", var_scale_attr)
      if (nc_err == NF90_NOERR) arr = arr*var_scale_attr
      nc_err = NF90_GET_ATT(ncid, nc_var, "add_offset", var_offset_attr)
      if (nc_err == NF90_NOERR) arr = arr + var_offset_attr
   end subroutine read_netcdf_var_with_attrs

end module mod_nws13
