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

   implicit none

   private

   !> @brief Flag value indicating missing or invalid data
   real(8), parameter  :: null_flag_value = -99999999.0d0

   !> @brief Variable names used in the NetCDF file. Maybe at some point we make these dynamic
   character(2), parameter :: dim_name_lon = "xi" !< Dimension name for longitude
   character(2), parameter :: dim_name_lat = "yi" !< Dimension name for latitude
   character(3), parameter :: var_name_lon = "lon" !< Variable name for longitude
   character(3), parameter :: var_name_lat = "lat" !< Variable name for latitude
   character(3), parameter :: var_name_u_wind = "U10" !< Variable name for U wind component
   character(3), parameter :: var_name_v_wind = "V10" !< Variable name for V wind component
   character(4), parameter :: var_name_pressure = "PSFC" !< Variable name for pressure
   character(4), parameter :: var_name_time = "time" !< Variable name for time
   character(4), parameter :: var_name_storm_center_lon = "clon" !< Variable name for storm track longitude
   character(4), parameter :: var_name_storm_center_lat = "clat" !< Variable name for storm track latitude

   !> @brief Container for meteorological data from NetCDF wind grid
   type t_windData
      real(8), allocatable :: lon(:, :) !< Grid longitude coordinates (degrees)
      real(8), allocatable :: lat(:, :) !< Grid latitude coordinates (degrees)
      real(8), allocatable :: u(:, :) !< Wind velocity, x-direction (m/s)
      real(8), allocatable :: v(:, :) !< Wind velocity, y-direction (m/s)
      real(8), allocatable :: p(:, :) !< Atmospheric pressure (converted to m water)
      real(8), allocatable :: ws(:, :) !< Wind speed magnitude (m/s)
      real(8) :: storm_pos(2) = [0.0d0, 0.0d0] !< Storm center position [lon, lat]
      integer :: file_snap = -1 !< File snap that this data corresponds to (-1 = uninitialized)
   contains
      final :: wind_data_destructor
   end type t_windData

   !> @brief Spatial hash acceleration structure for grid cell location
   type t_spatial_hash
      integer :: nx_hash !< Number of hash cells in x-direction
      integer :: ny_hash !< Number of hash cells in y-direction
      real(8) :: min_lon !< Minimum longitude of grid
      real(8) :: max_lon !< Maximum longitude of grid
      real(8) :: min_lat !< Minimum latitude of grid
      real(8) :: max_lat !< Maximum latitude of grid
      real(8) :: dx_hash !< Hash cell width
      real(8) :: dy_hash !< Hash cell height
      type(t_cell_list), allocatable :: hash_cells(:, :) !< Hash cells containing grid cell lists
   end type t_spatial_hash

   !> @brief List of grid cells overlapping a hash cell
   type t_cell_list
      integer :: n_cells = 0 !< Number of grid cells in list
      integer, allocatable :: i_indices(:) !< Grid i-indices
      integer, allocatable :: j_indices(:) !< Grid j-indices
   end type t_cell_list

   !> @brief Stores bilinear interpolation data for grid-to-mesh mapping
   type t_interpolationData
      integer, allocatable :: ilon(:) !< Grid x-indices for each mesh node
      integer, allocatable :: ilat(:) !< Grid y-indices for each mesh node
      real(8), allocatable :: weights(:, :) !< Bilinear interpolation weights (NP, 4)
      logical :: initialized = .false. !< Flag indicating if interpolation is computed
      logical :: is_moving = .false. !< Flag indicating if grid moves with time
   end type t_interpolationData

   !> @brief Stores storm center data for a group to avoid repeated NetCDF reads
   type t_stormCenterData
      logical :: available = .false. !< Flag indicating if storm center data exists in file
      logical :: loaded = .false. !< Flag indicating if data has been loaded into memory
      real(8), allocatable :: position(:, :) !< Storm center longitude for each time snapshot
   end type t_stormCenterData

   !> @brief Container for interpolated node values
   type t_nodeInterpolationResult
      real(8) :: u = 0.0d0 !< Wind velocity U component
      real(8) :: v = 0.0d0 !< Wind velocity V component
      real(8) :: p = 0.0d0 !< Atmospheric pressure
      logical :: valid = .false. !< Flag indicating if interpolation was successful
   end type t_nodeInterpolationResult

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
      type(t_windData) :: PrevData !< Data corresponding to PrevSnap
      type(t_windData) :: NextData !< Data corresponding to NextSnap
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
   type(t_stormCenterData) function storm_center_constructor(nc_idg, num_snaps) result(storm_center)
      use netcdf, only: NF90_INQ_VARID, NF90_GET_VAR, NF90_NOERR
      use netcdf_error, only: check_err
      implicit none
      integer, intent(in) :: nc_idg, num_snaps
      integer :: nc_var, nc_err

      ! Allocate arrays
      allocate (storm_center%position(num_snaps, 2))

      ! Check if storm center longitude data exists
      nc_err = NF90_INQ_VARID(nc_idg, var_name_storm_center_lon, nc_var)
      if (nc_err == NF90_NOERR) then
         storm_center%available = .true.
         call check_err(NF90_GET_VAR(nc_idg, nc_var, storm_center%position(:, 1), &
                                     START=[1], COUNT=[num_snaps]))

         ! Check if storm center latitude data exists
         nc_err = NF90_INQ_VARID(nc_idg, var_name_storm_center_lat, nc_var)
         if (nc_err == NF90_NOERR) then
            call check_err(NF90_GET_VAR(nc_idg, nc_var, storm_center%position(:, 2), &
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
                        nf90_noerr, nf90_inq_grpname, NF90_MAX_DIMS
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
      call check_err(NF90_INQ_DIMID(NC_IDG, var_name_time, NC_DIM))
      call check_err(NF90_INQUIRE_DIMENSION(NC_IDG, NC_DIM, LEN=grp%NumSnap))
      call check_err(NF90_INQ_VARID(NC_IDG, var_name_time, NC_VAR))

      ! Read and store grid dimensions
      call check_err(NF90_INQ_DIMID(NC_IDG, dim_name_lon, NC_DIM))
      call check_err(NF90_INQUIRE_DIMENSION(NC_IDG, NC_DIM, LEN=grp%NumLon))
      call check_err(NF90_INQ_DIMID(NC_IDG, dim_name_lat, NC_DIM))
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
      call check_err(NF90_INQ_VARID(NC_IDG, var_name_lon, NC_VAR))
      allocate (dimids(NF90_MAX_DIMS))
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
      call check_err(NF90_INQ_VARID(NC_IDG, var_name_time, NC_VAR))
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

      call check_err(NF90_INQ_DIMID(NC_IDG, dim_name_lon, NC_DIM))
      call check_err(NF90_INQUIRE_DIMENSION(NC_IDG, NC_DIM, LEN=NumLon))
      call check_err(NF90_INQ_DIMID(NC_IDG, dim_name_lat, NC_DIM))
      call check_err(NF90_INQUIRE_DIMENSION(NC_IDG, NC_DIM, LEN=NumLat))

      allocate (Lon(1:NumLon, 1:NumLat))
      allocate (Lat(1:NumLon, 1:NumLat))

      call check_err(NF90_INQ_VARID(NC_IDG, var_name_lon, NC_VAR))
      call check_err(NF90_GET_VAR(NC_IDG, NC_VAR, Lon(:, :), &
                                  START=[1, 1], COUNT=[NumLon, NumLat]))
      call check_err(NF90_INQ_VARID(NC_IDG, var_name_lat, NC_VAR))
      call check_err(NF90_GET_VAR(NC_IDG, NC_VAR, Lat(:, :), &
                                  START=[1, 1], COUNT=[NumLon, NumLat]))

      interp_data = NWS13INTERP(NumLon, NumLat, Lon, Lat)
      interp_data%is_moving = .false.

      deallocate (Lon, Lat)

   end function read_2d_grid

!-----------------------------------------------------------------------
!> @brief Reads lon/lat coordinates for moving grids and updates interpolation
!> @param[in]    nc_idg    NetCDF group ID
!> @param[in]    currsnap  Current time snapshot index
!> @param[inout] grp       NWS13 group structure to update
!> @param[inout] data      Wind data structure to populate with grid coordinates
!-----------------------------------------------------------------------
   subroutine read_moving_grid(currsnap, grp, data)
      use netcdf, only: NF90_INQ_VARID, NF90_GET_VAR
      use netcdf_error, only: check_err

      implicit none

      integer, intent(in) :: currsnap
      type(t_nws13), intent(in) :: grp
      type(t_windData), intent(out) :: data

      integer :: nc_var

      ! Read lon/lat coordinates for this time snapshot
      call check_err(NF90_INQ_VARID(grp%group_id, var_name_lon, nc_var))
      call check_err(NF90_GET_VAR(grp%group_id, nc_var, data%lon(:, :), &
                                  START=[1, 1, currsnap], &
                                  COUNT=[grp%NumLon, grp%NumLat, 1]))
      call check_err(NF90_INQ_VARID(grp%group_id, var_name_lat, nc_var))
      call check_err(NF90_GET_VAR(grp%group_id, nc_var, data%lat(:, :), &
                                  START=[1, 1, currsnap], &
                                  COUNT=[grp%NumLon, grp%NumLat, 1]))
   end subroutine read_moving_grid

!-----------------------------------------------------------------------
!> @brief Reads a single wind snapshot from NetCDF file
!> @param[in]    snap_idx  Snapshot index to read
!> @param[in]    grid_info NWS13 grid information
!> @return       t_windData object with read snapshot data
!-----------------------------------------------------------------------
   type(t_windData) function read_wind_snapshot(snap_idx, grid_info) result(wind_data)
      use, intrinsic :: iso_c_binding, only: c_int
      use global, only: screenMessage, INFO
      use ADC_CONSTANTS, only: G, RHOWAT0
      implicit none
      integer, intent(in) :: snap_idx
      type(t_nws13), intent(in) :: grid_info
      character(256) :: message

      ! Log processing message
      write (message, '(A,A,A,I0,2A)') "Processing group '", trim(grid_info%group_name), &
         "' with rank ", grid_info%rank, " for wind/pressure with time ", &
         trim(grid_info%SnapTimes(snap_idx)%to_iso_string())
      call screenMessage(INFO, message)

      ! Initialize wind data object
      wind_data = t_windData(grid_info%NumLon, grid_info%NumLat, grid_info%interp%is_moving)

      ! Read moving grid coordinates if applicable
      if (grid_info%interp%is_moving) then
         call read_moving_grid(snap_idx, grid_info, wind_data)
      end if

      ! Read wind and pressure fields
      call read_netcdf_var_with_attrs(grid_info%group_id, &
                                      var_name_u_wind, snap_idx, grid_info%NumLon, grid_info%NumLat, wind_data%U)
      call read_netcdf_var_with_attrs(grid_info%group_id, &
                                      var_name_v_wind, snap_idx, grid_info%NumLon, grid_info%NumLat, wind_data%V)
      call read_netcdf_var_with_attrs(grid_info%group_id, &
                                      var_name_pressure, snap_idx, grid_info%NumLon, grid_info%NumLat, wind_data%P)

      ! Apply wind multiplier
      wind_data%U = wind_data%U*NWS13WindMultiplier
      wind_data%V = wind_data%V*NWS13WindMultiplier

      ! Convert pressure to meters of water
      wind_data%P = wind_data%P*100.d0/RHOWAT0/G

      ! Get storm center position if available
      if (grid_info%storm_center%available .and. grid_info%storm_center%loaded) then
         wind_data%storm_pos = grid_info%storm_center%position(snap_idx, 1:2)
      else
         wind_data%storm_pos = [0.0d0, 0d0]
      end if

      ! Store the snapshot index this data corresponds to
      wind_data%file_snap = snap_idx

   end function read_wind_snapshot

!-----------------------------------------------------------------------
!> @brief Advances group to the snapshot containing current_time
!> @param[inout] group NWS13 group to advance
!> @param[in] current_time Target time
!-----------------------------------------------------------------------
   subroutine advance_to_current_snapshot(group, current_time)
      use mod_datetime, only: t_datetime, operator(>)
      implicit none
      type(t_nws13), intent(inout) :: group
      type(t_datetime), intent(in) :: current_time

      do while (current_time > group%NextTime .and. group%NextSnap < group%NumSnap)
         group%PrevSnap = group%NextSnap
         group%PrevTime = group%NextTime
         group%NextSnap = group%NextSnap + 1
         group%NextTime = group%SnapTimes(group%NextSnap)
      end do
   end subroutine advance_to_current_snapshot

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
      groups(:)%InclSnap = 0

      do ig = 1, size(groups)
         ! Skip groups that haven't started
         if (current_time < groups(ig)%PrevTime) cycle

         ! Advance to correct snapshot
         call advance_to_current_snapshot(groups(ig), current_time)

         ! Mark active groups and update next time
         if (groups(ig)%PrevTime <= current_time .and. &
             current_time <= groups(ig)%NextTime) then
            groups(ig)%InclSnap = 1
            if (next_time == null_datetime() .or. groups(ig)%NextTime < next_time) then
               next_time = groups(ig)%NextTime
            end if
         end if
      end do
   end function update_active_groups

!-----------------------------------------------------------------------
!> @brief Interpolates wind data between two time snapshots
!> @param[in] grid_info Grid information containing snapshots and timing
!> @param[in] current_time Current simulation time
!> @return Interpolated wind data at current time
!-----------------------------------------------------------------------
   type(t_windData) pure function interpolate_wind_in_time(grid_info, current_time) result(current_data)
      use mod_datetime, only: t_datetime, t_timedelta, operator(-)
      implicit none
      type(t_nws13), intent(in) :: grid_info
      type(t_datetime), intent(in) :: current_time
      type(t_timedelta) :: time_interp_numerator, time_interp_denominator
      real(8) :: time_interp_factor

      ! Initialize current data structure
      current_data = t_windData(grid_info%NumLon, grid_info%NumLat, grid_info%interp%is_moving)

      ! Calculate time interpolation factor
      time_interp_numerator = current_time - grid_info%PrevTime
      time_interp_denominator = grid_info%NextTime - grid_info%PrevTime
      time_interp_factor = dble(time_interp_numerator%total_milliseconds())/ &
                           dble(time_interp_denominator%total_milliseconds())

      ! Interpolate grid coordinates if moving
      if (grid_info%interp%is_moving) then
         current_data%Lon = grid_info%PrevData%Lon + (grid_info%NextData%Lon - grid_info%PrevData%Lon)*time_interp_factor
         current_data%Lat = grid_info%PrevData%Lat + (grid_info%NextData%Lat - grid_info%PrevData%Lat)*time_interp_factor
      end if

      ! Interpolate wind and pressure fields
      current_data%U = grid_info%PrevData%U + (grid_info%NextData%U - grid_info%PrevData%U)*time_interp_factor
      current_data%V = grid_info%PrevData%V + (grid_info%NextData%V - grid_info%PrevData%V)*time_interp_factor
      current_data%P = grid_info%PrevData%P + (grid_info%NextData%P - grid_info%PrevData%P)*time_interp_factor

      ! Calculate and interpolate scalar wind speed
      current_data%Ws = sqrt(grid_info%PrevData%U**2d0 + grid_info%PrevData%V**2d0) + &
                        (sqrt(grid_info%NextData%U**2d0 + grid_info%NextData%V**2d0) - &
                         sqrt(grid_info%PrevData%U**2d0 + grid_info%PrevData%V**2d0))*time_interp_factor

      ! Interpolate storm center position
      current_data%storm_pos = grid_info%PrevData%storm_pos + &
                               (grid_info%NextData%storm_pos - grid_info%PrevData%storm_pos)*time_interp_factor

   end function interpolate_wind_in_time

!-----------------------------------------------------------------------
!> @brief Interpolates wind data to a single mesh node
!> @param[in] node_idx Node index
!> @param[in] wind_data Wind data on grid
!> @param[in] interp Interpolation weights and indices
!> @return t_nodeInterpolationResult containing interpolated values and validity flag
!-----------------------------------------------------------------------
   type(t_nodeInterpolationResult) pure function interpolate_to_node(node_idx, wind_data, interp) result(node_result)
      implicit none
      integer, intent(in) :: node_idx
      type(t_windData), intent(in) :: wind_data
      type(t_interpolationData), intent(in) :: interp
      real(8) :: ws_corners(4), p_corners(4), wind_speed, wind_dir
      real(8), parameter :: eps = epsilon(1d0)
      integer :: ilon, ilat

      ! Initialize result as invalid
      node_result%valid = .false.
      node_result%u = 0.0d0
      node_result%v = 0.0d0
      node_result%p = 0.0d0

      ! Check if this node has valid interpolation data
      if (interp%ilon(node_idx) <= 0) return

      ilon = interp%ilon(node_idx)
      ilat = interp%ilat(node_idx)

      ! Get corner values for validation
      ws_corners(1) = wind_data%Ws(ilon, ilat)
      ws_corners(2) = wind_data%Ws(ilon + 1, ilat)
      ws_corners(3) = wind_data%Ws(ilon + 1, ilat + 1)
      ws_corners(4) = wind_data%Ws(ilon, ilat + 1)

      p_corners(1) = wind_data%P(ilon, ilat)
      p_corners(2) = wind_data%P(ilon + 1, ilat)
      p_corners(3) = wind_data%P(ilon + 1, ilat + 1)
      p_corners(4) = wind_data%P(ilon, ilat + 1)

      ! Skip if any corner has flag values
      if (.not. (all(abs(ws_corners - null_flag_value) > eps) .or. &
                 all(abs(p_corners - null_flag_value) > eps))) return

      ! Perform bilinear interpolation
      node_result%u = interpolate_bilinear(wind_data%U, interp%weights(node_idx, :), ilon, ilat)
      node_result%v = interpolate_bilinear(wind_data%V, interp%weights(node_idx, :), ilon, ilat)
      wind_speed = interpolate_bilinear(wind_data%Ws, interp%weights(node_idx, :), ilon, ilat)
      node_result%p = interpolate_bilinear(wind_data%P, interp%weights(node_idx, :), ilon, ilat)

      ! Adjust U/V based on scalar wind magnitude
      if (abs(node_result%u) <= eps .and. abs(node_result%v) <= eps) then
         wind_dir = 0d0
      else
         wind_dir = atan2(node_result%u, node_result%v)
         node_result%u = wind_speed*sin(wind_dir)
         node_result%v = wind_speed*cos(wind_dir)
      end if

      node_result%valid = .true.

   end function interpolate_to_node

!-----------------------------------------------------------------------
!> @brief Processes a single wind group and applies it to the mesh
!> @param[in] group_idx Group index in sorted order
!> @param[in] current_time Current simulation time
!> @param[inout] wvnx2 Wind X-component array
!> @param[inout] wvny2 Wind Y-component array
!> @param[inout] prn2 Pressure array
!> @param[inout] eye_lon Storm center longitude array
!> @param[inout] eye_lat Storm center latitude array
!> @param[out] found_eye Whether storm eye was found
!-----------------------------------------------------------------------
   subroutine process_wind_group(group_idx, current_time, wvnx2, wvny2, prn2, &
                                 PowellGroupTemp, eye_lon, eye_lat, found_eye)
      use MESH, only: NP
      use mod_datetime, only: t_datetime, t_timedelta, operator(-), operator(==)
      implicit none
      integer, intent(in) :: group_idx
      type(t_datetime), intent(in) :: current_time
      real(8), intent(inout) :: wvnx2(NP), wvny2(NP), prn2(NP)
      integer, intent(inout) :: PowellGroupTemp
      real(8), intent(inout) :: eye_lon(3), eye_lat(3)
      logical, intent(out) :: found_eye

      integer :: ig, iv
      real(8), parameter :: eps = epsilon(1d0)
      type(t_windData) :: prev_data, next_data, current_data
      type(t_interpolationData) :: current_interp
      type(t_nodeInterpolationResult) :: node_result

      ! Get actual group index
      ig = sorted_group_indices(group_idx)

      ! Skip inactive groups
      if (NWS13(ig)%InclSnap == 0) return

      ! Read snapshots (storm position is now included in wind data)
      if (NWS13(ig)%NextData%file_snap /= -1 .and. NWS13(ig)%PrevSnap == NWS13(ig)%NextData%file_snap) then
         NWS13(ig)%PrevData = NWS13(ig)%NextData
      else if (NWS13(ig)%PrevData%file_snap == -1 .or. NWS13(ig)%PrevSnap /= NWS13(ig)%PrevData%file_snap) then
         NWS13(ig)%PrevData = read_wind_snapshot(NWS13(ig)%PrevSnap, NWS13(ig))
      end if

      if (NWS13(ig)%NextData%file_snap == -1 .or. NWS13(ig)%NextData%file_snap /= NWS13(ig)%NextSnap) then
         NWS13(ig)%NextData = read_wind_snapshot(NWS13(ig)%NextSnap, NWS13(ig))
      end if

      ! Check for storm eye data
      found_eye = NWS13(ig)%storm_center%available .and. NWS13(ig)%storm_center%loaded

      ! Time interpolation (includes storm center interpolation)
      current_data = interpolate_wind_in_time(NWS13(ig), current_time)

      ! Get interpolation data
      if (NWS13(ig)%interp%is_moving) then
         current_interp = NWS13INTERP(NWS13(ig)%NumLon, NWS13(ig)%Numlat, &
                                      current_data%lon, current_data%lat)
      else
         current_interp = NWS13(ig)%interp
      end if

      ! Spatial interpolation to mesh nodes
      do iv = 1, NP
         ! Only process nodes not yet assigned
         if (abs(prn2(iv) - null_flag_value) > eps) cycle

         node_result = interpolate_to_node(iv, current_data, current_interp)
         if (node_result%valid) then
            wvnx2(iv) = node_result%u
            wvny2(iv) = node_result%v
            prn2(iv) = node_result%p
         end if
      end do

      ! Update storm eye position for Powell drag
      if (found_eye) then
         if ((NWS13GroupForPowell == 0) .or. (NWS13GroupForPowell == ig)) then
            eye_lon(3) = current_data%storm_pos(1)
            eye_lat(3) = current_data%storm_pos(2)
            PowellGroupTemp = ig
         end if
      end if

   end subroutine process_wind_group

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

      integer :: IG_sorted
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
            call process_wind_group(IG_sorted, CurrentTime, WVNX2, WVNY2, PRN2, &
                                    PowellGroupTemp, EyeLonR, EyeLatR, FoundEye)
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

      integer :: IA
      integer :: IO
      integer :: IV
      integer, dimension(2) :: cell_indices

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

      type(t_spatial_hash) :: hash_accel

      call setMessageSource("nws13interp")
#if defined(OWIWIND_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      ! Build spatial hash acceleration structure
      hash_accel = build_spatial_hash(NumLon, NumLat, Lon, Lat)

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

         ! Find containing cell using spatial hash + walking
         cell_indices = find_containing_cell(hash_accel, AdcLon, AdcLat, &
                                             NumLon, NumLat, Lon, Lat)

         if (cell_indices(1) > 0) then
            IO = cell_indices(1)
            IA = cell_indices(2)

            ! Compute bilinear interpolation weights
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

            interpData%ilon(IV) = IO
            interpData%ilat(IV) = IA
            interpData%weights(IV, 1) = W1/(W0*2.d0)
            interpData%weights(IV, 2) = W2/(W0*2.d0)
            interpData%weights(IV, 3) = W3/(W0*2.d0)
            interpData%weights(IV, 4) = W4/(W0*2.d0)
         end if

      end do

      ! Clean up spatial hash
      call destroy_spatial_hash(hash_accel)

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

!-----------------------------------------------------------------------
!> @brief Builds spatial hash acceleration structure for grid cell location
!> @param[in] nx Grid x-dimension
!> @param[in] ny Grid y-dimension
!> @param[in] lon Grid longitudes
!> @param[in] lat Grid latitudes
!> @return hash_struct Initialized spatial hash structure
!-----------------------------------------------------------------------
   type(t_spatial_hash) function build_spatial_hash(nx, ny, lon, lat) result(hash_struct)
      implicit none
      integer, intent(in) :: nx, ny
      real(8), intent(in) :: lon(nx, ny), lat(nx, ny)
      integer :: i, j, ih, jh, ih_min, ih_max, jh_min, jh_max
      real(8) :: avg_dx, avg_dy, cell_min_lon, cell_max_lon, cell_min_lat, cell_max_lat

      ! Compute grid bounds
      hash_struct%min_lon = minval(lon)
      hash_struct%max_lon = maxval(lon)
      hash_struct%min_lat = minval(lat)
      hash_struct%max_lat = maxval(lat)

      ! Estimate average grid spacing
      avg_dx = (hash_struct%max_lon - hash_struct%min_lon)/real(nx - 1, 8)
      avg_dy = (hash_struct%max_lat - hash_struct%min_lat)/real(ny - 1, 8)

      ! Set hash cell size to ~2.5x average grid spacing for overlap
      hash_struct%dx_hash = avg_dx*2.5d0
      hash_struct%dy_hash = avg_dy*2.5d0

      ! Compute hash grid dimensions
      hash_struct%nx_hash = ceiling((hash_struct%max_lon - hash_struct%min_lon)/hash_struct%dx_hash) + 1
      hash_struct%ny_hash = ceiling((hash_struct%max_lat - hash_struct%min_lat)/hash_struct%dy_hash) + 1

      ! Allocate hash grid
      allocate (hash_struct%hash_cells(hash_struct%nx_hash, hash_struct%ny_hash))

      ! Initialize cell lists with estimated capacity
      do ih = 1, hash_struct%nx_hash
         do jh = 1, hash_struct%ny_hash
            hash_struct%hash_cells(ih, jh)%n_cells = 0
            allocate (hash_struct%hash_cells(ih, jh)%i_indices(16))
            allocate (hash_struct%hash_cells(ih, jh)%j_indices(16))
         end do
      end do

      ! Map grid cells to hash cells
      do i = 1, nx - 1
         do j = 1, ny - 1
            ! Get bounding box of this grid cell
            cell_min_lon = min(lon(i, j), lon(i + 1, j), lon(i + 1, j + 1), lon(i, j + 1))
            cell_max_lon = max(lon(i, j), lon(i + 1, j), lon(i + 1, j + 1), lon(i, j + 1))
            cell_min_lat = min(lat(i, j), lat(i + 1, j), lat(i + 1, j + 1), lat(i, j + 1))
            cell_max_lat = max(lat(i, j), lat(i + 1, j), lat(i + 1, j + 1), lat(i, j + 1))

            ! Find overlapping hash cells
            ih_min = max(1, floor((cell_min_lon - hash_struct%min_lon)/hash_struct%dx_hash) + 1)
            ih_max = min(hash_struct%nx_hash, ceiling((cell_max_lon - hash_struct%min_lon)/hash_struct%dx_hash) + 1)
            jh_min = max(1, floor((cell_min_lat - hash_struct%min_lat)/hash_struct%dy_hash) + 1)
            jh_max = min(hash_struct%ny_hash, ceiling((cell_max_lat - hash_struct%min_lat)/hash_struct%dy_hash) + 1)

            ! Add to all overlapping hash cells
            do ih = ih_min, ih_max
               do jh = jh_min, jh_max
                  call add_to_cell_list(hash_struct%hash_cells(ih, jh), i, j)
               end do
            end do
         end do
      end do

   end function build_spatial_hash

!-----------------------------------------------------------------------
!> @brief Adds a grid cell to a hash cell's list
!> @param[inout] cell_list Hash cell to add to
!> @param[in] i Grid i-index
!> @param[in] j Grid j-index
!-----------------------------------------------------------------------
   subroutine add_to_cell_list(cell_list, i, j)
      implicit none
      type(t_cell_list), intent(inout) :: cell_list
      integer, intent(in) :: i, j
      integer, allocatable :: temp_i(:), temp_j(:)
      integer :: new_size

      cell_list%n_cells = cell_list%n_cells + 1

      ! Resize arrays if needed
      if (cell_list%n_cells > size(cell_list%i_indices)) then
         new_size = size(cell_list%i_indices)*2
         allocate (temp_i(new_size), temp_j(new_size))
         temp_i(1:cell_list%n_cells - 1) = cell_list%i_indices(1:cell_list%n_cells - 1)
         temp_j(1:cell_list%n_cells - 1) = cell_list%j_indices(1:cell_list%n_cells - 1)
         deallocate (cell_list%i_indices, cell_list%j_indices)
         cell_list%i_indices = temp_i
         cell_list%j_indices = temp_j
      end if

      cell_list%i_indices(cell_list%n_cells) = i
      cell_list%j_indices(cell_list%n_cells) = j

   end subroutine add_to_cell_list

!-----------------------------------------------------------------------
!> @brief Destroys spatial hash structure and frees memory
!> @param[inout] hash_struct Spatial hash to destroy
!-----------------------------------------------------------------------
   subroutine destroy_spatial_hash(hash_struct)
      implicit none
      type(t_spatial_hash), intent(inout) :: hash_struct
      integer :: i, j

      if (allocated(hash_struct%hash_cells)) then
         do i = 1, hash_struct%nx_hash
            do j = 1, hash_struct%ny_hash
               if (allocated(hash_struct%hash_cells(i, j)%i_indices)) then
                  deallocate (hash_struct%hash_cells(i, j)%i_indices)
                  deallocate (hash_struct%hash_cells(i, j)%j_indices)
               end if
            end do
         end do
         deallocate (hash_struct%hash_cells)
      end if

   end subroutine destroy_spatial_hash

!-----------------------------------------------------------------------
!> @brief Finds grid cell containing query point using spatial hash + walking
!> @param[in] hash_struct Spatial hash structure
!> @param[in] query_lon Query point longitude
!> @param[in] query_lat Query point latitude
!> @param[in] nx Grid x-dimension
!> @param[in] ny Grid y-dimension
!> @param[in] lon Grid longitudes
!> @param[in] lat Grid latitudes
!> @return indices [i,j] indices of containing cell, or [-1,-1] if not found
!-----------------------------------------------------------------------
   function find_containing_cell(hash_struct, query_lon, query_lat, nx, ny, lon, lat) result(indices)
      implicit none
      type(t_spatial_hash), intent(in) :: hash_struct
      real(8), intent(in) :: query_lon, query_lat
      integer, intent(in) :: nx, ny
      real(8), intent(in) :: lon(nx, ny), lat(nx, ny)
      integer :: indices(2)
      integer :: ih, jh, ic, i, j, best_i, best_j
      real(8) :: dist, min_dist, w0, w1, w2, w3, w4
      logical :: found

      indices = [-1, -1]

      ! Get hash cell for query point
      ih = min(max(1, floor((query_lon - hash_struct%min_lon)/hash_struct%dx_hash) + 1), hash_struct%nx_hash)
      jh = min(max(1, floor((query_lat - hash_struct%min_lat)/hash_struct%dy_hash) + 1), hash_struct%ny_hash)

      ! Check candidate cells from hash
      found = .false.
      min_dist = huge(1.0d0)
      best_i = -1
      best_j = -1

      do ic = 1, hash_struct%hash_cells(ih, jh)%n_cells
         i = hash_struct%hash_cells(ih, jh)%i_indices(ic)
         j = hash_struct%hash_cells(ih, jh)%j_indices(ic)

         ! Quick point-in-quad test using weights
         w0 = compute_weight(lon(i, j), lat(i, j), lon(i + 1, j), lat(i + 1, j), &
                             lon(i + 1, j + 1), lat(i + 1, j + 1), lon(i, j + 1), lat(i, j + 1))
         w1 = compute_weight(query_lon, query_lat, lon(i + 1, j), lat(i + 1, j), &
                             lon(i + 1, j + 1), lat(i + 1, j + 1), lon(i, j + 1), lat(i, j + 1))
         w2 = compute_weight(lon(i, j), lat(i, j), query_lon, query_lat, &
                             lon(i + 1, j + 1), lat(i + 1, j + 1), lon(i, j + 1), lat(i, j + 1))
         w3 = compute_weight(lon(i, j), lat(i, j), lon(i + 1, j), lat(i + 1, j), &
                             query_lon, query_lat, lon(i, j + 1), lat(i, j + 1))
         w4 = compute_weight(lon(i, j), lat(i, j), lon(i + 1, j), lat(i + 1, j), &
                             lon(i + 1, j + 1), lat(i + 1, j + 1), query_lon, query_lat)

         if ((w1 + w2 + w3 + w4) <= (w0*2.001d0)) then
            indices = [i, j]
            found = .true.
            exit
         end if

         ! Track closest cell for walking fallback
         dist = (query_lon - 0.25d0*(lon(i, j) + lon(i + 1, j) + lon(i + 1, j + 1) + lon(i, j + 1)))**2 + &
                (query_lat - 0.25d0*(lat(i, j) + lat(i + 1, j) + lat(i + 1, j + 1) + lat(i, j + 1)))**2
         if (dist < min_dist) then
            min_dist = dist
            best_i = i
            best_j = j
         end if
      end do

      ! If not found in hash, try walking from closest candidate
      if (.not. found .and. best_i > 0) then
         indices = walk_to_containing_cell(best_i, best_j, query_lon, query_lat, nx, ny, lon, lat)
      end if

   end function find_containing_cell

!-----------------------------------------------------------------------
!> @brief Walks grid to find cell containing query point
!> @param[in] start_i Starting i-index
!> @param[in] start_j Starting j-index
!> @param[in] query_lon Query longitude
!> @param[in] query_lat Query latitude
!> @param[in] nx Grid x-dimension
!> @param[in] ny Grid y-dimension
!> @param[in] lon Grid longitudes
!> @param[in] lat Grid latitudes
!> @return indices [i,j] of containing cell or [-1,-1]
!-----------------------------------------------------------------------
   function walk_to_containing_cell(start_i, start_j, query_lon, query_lat, nx, ny, lon, lat) result(indices)
      implicit none
      integer, intent(in) :: start_i, start_j, nx, ny
      real(8), intent(in) :: query_lon, query_lat
      real(8), intent(in) :: lon(nx, ny), lat(nx, ny)
      integer :: indices(2)
      integer :: i, j, di, dj, iter
      real(8) :: w0, w1, w2, w3, w4
      real(8) :: center_lon, center_lat, dx, dy
      logical :: found
      integer, parameter :: max_iter = 20

      i = start_i
      j = start_j
      indices = [-1, -1]
      found = .false.

      do iter = 1, max_iter
         ! Check bounds
         if (i < 1 .or. i >= nx .or. j < 1 .or. j >= ny) exit

         ! Test current cell
         w0 = compute_weight(lon(i, j), lat(i, j), lon(i + 1, j), lat(i + 1, j), &
                             lon(i + 1, j + 1), lat(i + 1, j + 1), lon(i, j + 1), lat(i, j + 1))
         w1 = compute_weight(query_lon, query_lat, lon(i + 1, j), lat(i + 1, j), &
                             lon(i + 1, j + 1), lat(i + 1, j + 1), lon(i, j + 1), lat(i, j + 1))
         w2 = compute_weight(lon(i, j), lat(i, j), query_lon, query_lat, &
                             lon(i + 1, j + 1), lat(i + 1, j + 1), lon(i, j + 1), lat(i, j + 1))
         w3 = compute_weight(lon(i, j), lat(i, j), lon(i + 1, j), lat(i + 1, j), &
                             query_lon, query_lat, lon(i, j + 1), lat(i, j + 1))
         w4 = compute_weight(lon(i, j), lat(i, j), lon(i + 1, j), lat(i + 1, j), &
                             lon(i + 1, j + 1), lat(i + 1, j + 1), query_lon, query_lat)

         if ((w1 + w2 + w3 + w4) <= (w0*2.001d0)) then
            indices = [i, j]
            found = .true.
            exit
         end if

         ! Determine walk direction based on cell center
         center_lon = 0.25d0*(lon(i, j) + lon(i + 1, j) + lon(i + 1, j + 1) + lon(i, j + 1))
         center_lat = 0.25d0*(lat(i, j) + lat(i + 1, j) + lat(i + 1, j + 1) + lat(i, j + 1))
         dx = query_lon - center_lon
         dy = query_lat - center_lat

         ! Walk in the direction of the query point
         di = 0
         dj = 0
         if (abs(dx) > abs(dy)) then
            if (dx > 0 .and. i < nx - 1) di = 1
            if (dx < 0 .and. i > 1) di = -1
         else
            if (dy > 0 .and. j < ny - 1) dj = 1
            if (dy < 0 .and. j > 1) dj = -1
         end if

         ! If no progress possible, try other direction
         if (di == 0 .and. dj == 0) then
            if (abs(dy) > 1d-10) then
               if (dy > 0 .and. j < ny - 1) dj = 1
               if (dy < 0 .and. j > 1) dj = -1
            elseif (abs(dx) > 1d-10) then
               if (dx > 0 .and. i < nx - 1) di = 1
               if (dx < 0 .and. i > 1) di = -1
            else
               exit
            end if
         end if

         i = i + di
         j = j + dj

         ! No movement possible
         if (di == 0 .and. dj == 0) exit
      end do

   end function walk_to_containing_cell

end module mod_nws13
