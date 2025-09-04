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
module mod_nws13_data
   use mod_datetime, only: t_datetime
   implicit none

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

   contains
      procedure, pass(self) :: find => find_containing_cell
      final :: destroy_spatial_hash

   end type t_spatial_hash

   !> @brief List of grid cells overlapping a hash cell
   type t_cell_list
      integer :: n_cells = 0 !< Number of grid cells in list
      integer, allocatable :: i_indices(:) !< Grid i-indices
      integer, allocatable :: j_indices(:) !< Grid j-indices

   contains
      procedure, pass(self) :: add => add_to_cell_list
   end type t_cell_list

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
   end type t_windData

   !> @brief Stores bilinear interpolation data for grid-to-mesh mapping
   type t_interpolationData
      integer, allocatable :: ilon(:) !< Grid x-indices for each mesh node
      integer, allocatable :: ilat(:) !< Grid y-indices for each mesh node
      real(8), allocatable :: weights(:, :) !< Bilinear interpolation weights (NP, 4)
      logical :: initialized = .false. !< Flag indicating if interpolation is computed
      logical :: is_moving = .false. !< Flag indicating if grid moves with time

   contains
      procedure, pass(self) :: interpolate => interpolate_to_node
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
   type t_nws13_group
      integer :: InclSnap !< Flag indicating if group is active for current time step
      integer :: NextSnap !< Index of next time snapshot in group
      integer :: NumSnap !< Total number of time snapshots in group
      integer :: PrevSnap !< Index of previous time snapshot in group
      integer :: rank !< Priority rank for overlapping grids
      integer :: file_index !< Group index in the file
      integer :: group_id !< NetCDF group identifier
      integer :: NumLon !< Grid x-dimension size
      integer :: NumLat !< Grid y-dimension size
      character(len=100) :: group_name !< NetCDF group name from file
      type(t_datetime) :: NextFileTime !< Date/time of next snapshot
      type(t_datetime) :: PrevFileTime !< Date/time of previous snapshot
      type(t_datetime), allocatable :: SnapTimes(:) !< Array of all snapshot times
      type(t_interpolationData) :: interp !< Grid-to-mesh interpolation data
      type(t_stormCenterData) :: storm_center !< Cached storm center data
      type(t_windData) :: PrevData !< Data corresponding to PrevSnap
      type(t_windData) :: NextData !< Data corresponding to NextSnap

   contains
      procedure, pass(self) :: process => t_nws13_process_group
      procedure, pass(self) :: advance => t_nws13_advance_to_time
      procedure, private, pass(self) :: read => t_nws13_read
      procedure, private, pass(self) :: read_snapshot => t_nws13_read_wind_snapshot
      procedure, private, pass(self) :: read_moving_grid => t_nws13_read_moving_grid
      procedure, private, pass(self) :: read_variable => t_nws13_read_netcdf_var_with_attrs
      procedure, private, pass(self) :: time_interpolate => t_nws13_interpolate_wind_in_time
   end type t_nws13_group

   interface t_cell_list
      module procedure t_cell_list_constructor
   end interface t_cell_list

   interface t_spatial_hash
      module procedure t_spatial_hash_constructor
   end interface t_spatial_hash

   interface t_nws13_group
      module procedure :: nws13_constructor
   end interface t_nws13_group

   interface t_stormCenterData
      module procedure :: storm_center_constructor
   end interface t_stormCenterData

   interface t_windData
      module procedure :: wind_data_constructor
   end interface t_windData

   interface t_interoplationData
      module procedure :: t_interpolationData_constructor
   end interface t_interoplationData

   private

   public :: t_nws13_group, null_flag_value

contains

   type(t_cell_list) function t_cell_list_constructor() result(cell_list)
      implicit none
      cell_list%n_cells = 0
      allocate (cell_list%i_indices(16))
      allocate (cell_list%j_indices(16))
   end function t_cell_list_constructor

   !-----------------------------------------------------------------------
   !> @brief Builds spatial hash acceleration structure for grid cell location
   !> @param[in] nx Grid x-dimension
   !> @param[in] ny Grid y-dimension
   !> @param[in] lon Grid longitudes
   !> @param[in] lat Grid latitudes
   !> @return hash_struct Initialized spatial hash structure
   !-----------------------------------------------------------------------
   type(t_spatial_hash) function t_spatial_hash_constructor(nx, ny, lon, lat) result(hash_struct)
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
            hash_struct%hash_cells(ih, jh) = t_cell_list()
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
                  call hash_struct%hash_cells(ih, jh)%add(i, j)
               end do
            end do
         end do
      end do

   end function t_spatial_hash_constructor

   !-----------------------------------------------------------------------
   !> @brief Adds a grid cell to a hash cell's list
   !> @param[inout] cell_list Hash cell to add to
   !> @param[in] i Grid i-index
   !> @param[in] j Grid j-index
   !-----------------------------------------------------------------------
   subroutine add_to_cell_list(self, i, j)
      implicit none
      class(t_cell_list), intent(inout) :: self
      integer, intent(in) :: i, j
      integer, allocatable :: temp_i(:), temp_j(:)
      integer :: new_size

      self%n_cells = self%n_cells + 1

      ! Resize arrays if needed
      if (self%n_cells > size(self%i_indices)) then
         new_size = size(self%i_indices)*2
         allocate (temp_i(new_size), temp_j(new_size))
         temp_i(1:self%n_cells - 1) = self%i_indices(1:self%n_cells - 1)
         temp_j(1:self%n_cells - 1) = self%j_indices(1:self%n_cells - 1)
         deallocate (self%i_indices, self%j_indices)
         self%i_indices = temp_i
         self%j_indices = temp_j
      end if

      self%i_indices(self%n_cells) = i
      self%j_indices(self%n_cells) = j

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
   function find_containing_cell(self, query_lon, query_lat, nx, ny, lon, lat) result(indices)
      implicit none
      class(t_spatial_hash), intent(in) :: self
      real(8), intent(in) :: query_lon, query_lat
      integer, intent(in) :: nx, ny
      real(8), intent(in) :: lon(nx, ny), lat(nx, ny)
      integer :: indices(2)
      integer :: ih, jh, ic, i, j, best_i, best_j
      real(8) :: dist, min_dist, w0, w1, w2, w3, w4
      logical :: found

      indices = [-1, -1]

      ! Get hash cell for query point
      ih = min(max(1, floor((query_lon - self%min_lon)/self%dx_hash) + 1), self%nx_hash)
      jh = min(max(1, floor((query_lat - self%min_lat)/self%dy_hash) + 1), self%ny_hash)

      ! Check candidate cells from hash
      found = .false.
      min_dist = huge(1.0d0)
      best_i = -1
      best_j = -1

      do ic = 1, self%hash_cells(ih, jh)%n_cells
         i = self%hash_cells(ih, jh)%i_indices(ic)
         j = self%hash_cells(ih, jh)%j_indices(ic)

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
   !> @brief Constructs a new t_nws13_group object from a NetCDF group
   !> @param[in]    nc_idg  NetCDF group ID to initialize from
   !> @return       grp     Initialized t_nws13_group structure
   !-----------------------------------------------------------------------
   type(t_nws13_group) function nws13_constructor(NC_IDG, file_index) result(grp)

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
      integer, intent(in) :: file_index
      integer :: NC_DIM
      integer :: NC_VAR
      integer :: NC_ERR
      integer :: ndims
      integer :: IGS
      integer, allocatable :: dimids(:)
      integer, allocatable :: TempI(:)
      type(t_datetime) :: WindRefDatetime

      grp%group_id = NC_IDG
      grp%file_index = file_index
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
      call check_err(NF90_GET_VAR(NC_IDG, NC_VAR, TempI(:), START=[1], COUNT=[grp%NumSnap]))
      do IGS = 1, grp%NumSnap
         grp%SnapTimes(IGS) = WindRefDatetime + t_timedelta(minutes=TempI(IGS))
      end do
      deallocate (TempI)

      grp%PrevFileTime = grp%SnapTimes(1)
      grp%PrevSnap = 1
      grp%NextFileTime = grp%SnapTimes(2)
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

      interp_data = t_interoplationData(NumLon, NumLat, Lon, Lat)
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
   type(t_windData) function t_nws13_read_moving_grid(self, current_snap) result(grid_data)
      use netcdf, only: NF90_INQ_VARID, NF90_GET_VAR
      use netcdf_error, only: check_err

      implicit none

      class(t_nws13_group), intent(in) :: self
      integer, intent(in) :: current_snap

      integer :: nc_var

      ! Allocate the grid data object
      allocate (grid_data%lon(self%NumLon, self%NumLat))
      allocate (grid_data%lat(self%NumLon, self%NumLat))

      ! Read lon/lat coordinates for this time snapshot
      call check_err(NF90_INQ_VARID(self%group_id, var_name_lon, nc_var))
      call check_err(NF90_GET_VAR(self%group_id, nc_var, grid_data%lon(:, :), &
                                  START=[1, 1, current_snap], &
                                  COUNT=[self%NumLon, self%NumLat, 1]))
      call check_err(NF90_INQ_VARID(self%group_id, var_name_lat, nc_var))
      call check_err(NF90_GET_VAR(self%group_id, nc_var, grid_data%lat(:, :), &
                                  START=[1, 1, current_snap], &
                                  COUNT=[self%NumLon, self%NumLat, 1]))
   end function t_nws13_read_moving_grid

   !-----------------------------------------------------------------------
   !> @brief Reads a single wind snapshot from NetCDF file
   !> @param[in]    snap_idx  Snapshot index to read
   !> @param[in]    grid_info NWS13 grid information
   !> @return       t_windData object with read snapshot data
   !-----------------------------------------------------------------------
   type(t_windData) function t_nws13_read_wind_snapshot(self, snap_idx, wind_multiplier) result(wind_data)
      use global, only: screenMessage, INFO
      use ADC_CONSTANTS, only: G, RHOWAT0
      implicit none

      class(t_nws13_group), intent(in) :: self
      integer, intent(in) :: snap_idx
      real(8), intent(in) :: wind_multiplier

      character(256) :: message

      write (message, '(A,A,A,I0,2A)') "Processing group '", trim(self%group_name), &
         "' with rank ", self%rank, " for wind/pressure with time ", &
         trim(self%SnapTimes(snap_idx)%to_iso_string())
      call screenMessage(INFO, message)

      ! Read moving grid coordinates if applicable
      if (self%interp%is_moving) then
         wind_data = self%read_moving_grid(snap_idx)
      else
         wind_data = t_windData(self%NumLon, self%NumLat, .false.)
      end if

      ! Read wind and pressure fields
      wind_data%U = self%read_variable(snap_idx, var_name_u_wind)*wind_multiplier
      wind_data%V = self%read_variable(snap_idx, var_name_v_wind)*wind_multiplier
      wind_data%P = self%read_variable(snap_idx, var_name_pressure)*100d0/RHOWAT0/G

      ! Get storm center position if available
      if (self%storm_center%available .and. self%storm_center%loaded) then
         wind_data%storm_pos = self%storm_center%position(snap_idx, 1:2)
      else
         wind_data%storm_pos = [0.0d0, 0d0]
      end if

      ! Store the snapshot index this data corresponds to
      wind_data%file_snap = snap_idx

   end function t_nws13_read_wind_snapshot

   !-----------------------------------------------------------------------
   !> @brief Advances group to the snapshot containing current_time
   !> @param[inout] group NWS13 group to advance
   !> @param[in] current_time Target time
   !-----------------------------------------------------------------------
   subroutine t_nws13_advance_to_time(self, current_time)
      use mod_datetime, only: t_datetime, operator(>), operator(<=), operator(<), operator(==), null_datetime
      implicit none
      class(t_nws13_group), intent(inout) :: self
      type(t_datetime), intent(in) :: current_time

      ! Skip if this group has not started
      if (current_time < self%PrevFileTime) return

      do while (current_time > self%NextFileTime .and. self%NextSnap < self%NumSnap)
         self%PrevSnap = self%NextSnap
         self%PrevFileTime = self%NextFileTime
         self%NextSnap = self%NextSnap + 1
         self%NextFileTime = self%SnapTimes(self%NextSnap)
      end do

      self%InclSnap = 0

      ! Mark active and update next time
      if (self%PrevFileTime <= current_time .and. current_time <= self%NextFileTime) then
         self%InclSnap = 1
      end if
   end subroutine t_nws13_advance_to_time

   subroutine t_nws13_read(self, wind_multiplier)
      class(t_nws13_group), intent(inout) :: self
      real(8), intent(in) :: wind_multiplier

      if (self%InclSnap /= 0) then
         if (self%NextData%file_snap /= -1 .and. self%PrevSnap == self%NextData%file_snap) then
            self%PrevData = self%NextData
         else if (self%PrevData%file_snap == -1 .or. self%PrevSnap /= self%PrevData%file_snap) then
            self%PrevData = self%read_snapshot(self%PrevSnap, wind_multiplier)
         end if

         if (self%NextData%file_snap == -1 .or. self%NextData%file_snap /= self%NextSnap) then
            self%NextData = self%read_snapshot(self%NextSnap, wind_multiplier)
         end if
      end if

   end subroutine t_nws13_read

   !-----------------------------------------------------------------------
   !> @brief Reads a NetCDF variable and applies fill value, scale, and offset
   !> @param[in]  ncid      NetCDF group ID
   !> @param[in]  varname   Name of the variable
   !> @param[in]  currsnap  Time index to read
   !> @param[in]  numlon    X-dimension size
   !> @param[in]  numlat    Y-dimension size
   !> @param[out] arr       Output array
   !-----------------------------------------------------------------------
   function t_nws13_read_netcdf_var_with_attrs(self, snap_index, varname) result(arr)
      use netcdf, only: NF90_INQ_VARID, NF90_GET_VAR, NF90_GET_ATT, NF90_NOERR
      use netcdf_error, only: check_err
      implicit none
      real(8), parameter :: eps = epsilon(1d0)

      class(t_nws13_group), intent(in) :: self
      integer, intent(in) :: snap_index
      character(len=*), intent(in) :: varname
      real(8) :: arr(self%numlon, self%numlat)

      integer :: nc_var, nc_err
      real(8) :: var_scale_attr, var_offset_attr, var_fillval

      call check_err(NF90_INQ_VARID(self%group_id, varname, nc_var))
      call check_err(NF90_GET_VAR(self%group_id, nc_var, arr(:, :), START=[1, 1, snap_index], COUNT=[self%numlon, self%numlat, 1]))

      nc_err = NF90_GET_ATT(self%group_id, nc_var, "_FillValue", var_fillval)
      if (nc_err == NF90_NOERR) where (abs(arr - var_fillval) <= eps) arr = null_flag_value

      nc_err = NF90_GET_ATT(self%group_id, nc_var, "scale_factor", var_scale_attr)
      if (nc_err == NF90_NOERR) arr = arr*var_scale_attr
      nc_err = NF90_GET_ATT(self%group_id, nc_var, "add_offset", var_offset_attr)
      if (nc_err == NF90_NOERR) arr = arr + var_offset_attr
   end function t_nws13_read_netcdf_var_with_attrs

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
   subroutine t_nws13_process_group(self, current_time, wvnx2, wvny2, prn2, &
                                    PowellGroupIndex, PowellGroupTempIndex, eye_lon, eye_lat, found_eye, wind_multiplier)
      use MESH, only: NP
      use mod_datetime, only: t_datetime, t_timedelta, operator(-), operator(==)
      implicit none
      class(t_nws13_group), intent(inout) :: self
      type(t_datetime), intent(in) :: current_time
      real(8), intent(in) :: wind_multiplier
      real(8), intent(inout) :: wvnx2(NP), wvny2(NP), prn2(NP)
      integer, intent(in) :: PowellGroupIndex
      integer, intent(inout) :: PowellGroupTempIndex
      real(8), intent(inout) :: eye_lon(3), eye_lat(3)
      logical, intent(out) :: found_eye

      integer :: iv
      real(8), parameter :: eps = epsilon(1d0)
      type(t_windData) :: current_data
      type(t_interpolationData) :: current_interp
      type(t_nodeInterpolationResult) :: node_result

      call self%read(wind_multiplier)

      ! Check for storm eye data
      found_eye = self%storm_center%available .and. self%storm_center%loaded

      ! Time interpolation (includes storm center interpolation)
      current_data = self%time_interpolate(current_time)

      ! Get interpolation data
      if (self%interp%is_moving) then
         current_interp = t_interoplationData(self%NumLon, self%Numlat, current_data%lon, current_data%lat)
      else
         current_interp = self%interp
      end if

      ! Spatial interpolation to mesh nodes
      do iv = 1, NP
         ! Only process nodes not yet assigned
         if (abs(prn2(iv) - null_flag_value) > eps) cycle

         node_result = current_interp%interpolate(iv, current_data)
         if (node_result%valid) then
            wvnx2(iv) = node_result%u
            wvny2(iv) = node_result%v
            prn2(iv) = node_result%p
         end if
      end do

      ! Update storm eye position for Powell drag
      if (found_eye) then
         if ((PowellGroupIndex == 0) .or. (PowellGroupIndex == self%file_index)) then
            eye_lon(3) = current_data%storm_pos(1)
            eye_lat(3) = current_data%storm_pos(2)
            PowellGroupTempIndex = self%file_index
         end if
      end if

   end subroutine t_nws13_process_group

   !-----------------------------------------------------------------------
   !> @brief Interpolates wind data between two time snapshots
   !> @param[in] grid_info Grid information containing snapshots and timing
   !> @param[in] current_time Current simulation time
   !> @return Interpolated wind data at current time
   !-----------------------------------------------------------------------
   type(t_windData) pure function t_nws13_interpolate_wind_in_time(self, current_time) result(current_data)
      use mod_datetime, only: t_datetime, t_timedelta, operator(-)
      implicit none
      class(t_nws13_group), intent(in) :: self
      type(t_datetime), intent(in) :: current_time
      type(t_timedelta) :: time_interp_numerator, time_interp_denominator
      real(8) :: time_interp_factor
      real(8), parameter :: eps = epsilon(1d0)

      ! Initialize current data structure
      current_data = t_windData(self%NumLon, self%NumLat, self%interp%is_moving)

      ! Calculate time interpolation factor
      time_interp_numerator = current_time - self%PrevFileTime
      time_interp_denominator = self%NextFileTime - self%PrevFileTime
      time_interp_factor = dble(time_interp_numerator%total_milliseconds())/ &
                           dble(time_interp_denominator%total_milliseconds())

      ! Interpolate grid coordinates if moving
      if (self%interp%is_moving) then
         current_data%Lon = self%PrevData%Lon + (self%NextData%Lon - self%PrevData%Lon)*time_interp_factor
         current_data%Lat = self%PrevData%Lat + (self%NextData%Lat - self%PrevData%Lat)*time_interp_factor
      end if

      ! Initialize fields with null values
      current_data%U = null_flag_value
      current_data%V = null_flag_value
      current_data%P = null_flag_value
      current_data%Ws = null_flag_value

      ! Interpolate all fields where both prev and next values are valid
      where (abs(self%PrevData%U - null_flag_value) > eps .and. &
             abs(self%NextData%U - null_flag_value) > eps .and. &
             abs(self%PrevData%V - null_flag_value) > eps .and. &
             abs(self%NextData%V - null_flag_value) > eps .and. &
             abs(self%PrevData%P - null_flag_value) > eps .and. &
             abs(self%NextData%P - null_flag_value) > eps)

         current_data%U = self%PrevData%U + (self%NextData%U - self%PrevData%U)*time_interp_factor
         current_data%V = self%PrevData%V + (self%NextData%V - self%PrevData%V)*time_interp_factor
         current_data%P = self%PrevData%P + (self%NextData%P - self%PrevData%P)*time_interp_factor
         current_data%Ws = sqrt(self%PrevData%U**2d0 + self%PrevData%V**2d0) + &
                           (sqrt(self%NextData%U**2d0 + self%NextData%V**2d0) - &
                            sqrt(self%PrevData%U**2d0 + self%PrevData%V**2d0))*time_interp_factor
      end where

      ! Interpolate storm center position
      current_data%storm_pos = self%PrevData%storm_pos + &
                               (self%NextData%storm_pos - self%PrevData%storm_pos)*time_interp_factor

   end function t_nws13_interpolate_wind_in_time

   !-----------------------------------------------------------------------
   !> @brief Computes interpolation weights for mapping grid data to mesh nodes
   !> @param[in]    numlon x-dimension size of grid
   !> @param[in]    numlat y-dimension size of grid
   !> @param[in]    lon    grid longitudes (2D)
   !> @param[in]    lat    grid latitudes (2D)
   !> @return       interpData Structure containing interpolation indices and weights
   !-----------------------------------------------------------------------
   type(t_interpolationData) function t_interpolationData_constructor(NumLon, NumLat, Lon, Lat) result(interpData)
      !-----------------------------------------------------------------------
      use ADC_CONSTANTS, only: RAD2DEG
      use MESH, only: NP, SLAM, SFEA

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

      ! Build spatial hash acceleration structure
      hash_accel = t_spatial_hash(NumLon, NumLat, Lon, Lat)

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
         cell_indices = hash_accel%find(AdcLon, AdcLat, NumLon, NumLat, Lon, Lat)

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

   end function t_interpolationData_constructor

   !-----------------------------------------------------------------------
   !> @brief Interpolates wind data to a single mesh node
   !> @param[in] node_idx Node index
   !> @param[in] wind_data Wind data on grid
   !> @param[in] interp Interpolation weights and indices
   !> @return t_nodeInterpolationResult containing interpolated values and validity flag
   !-----------------------------------------------------------------------
   type(t_nodeInterpolationResult) pure function interpolate_to_node(self, node_idx, wind_data) result(node_result)
      implicit none
      class(t_interpolationData), intent(in) :: self
      integer, intent(in) :: node_idx
      type(t_windData), intent(in) :: wind_data
      real(8) :: ws_corners(4), p_corners(4), this_weight(4)
      real(8) :: wind_speed, wind_dir
      real(8), parameter :: eps = epsilon(1d0)
      integer :: ilon, ilat

      ! Initialize result as invalid
      node_result%valid = .false.
      node_result%u = 0.0d0
      node_result%v = 0.0d0
      node_result%p = 0.0d0

      ! Check if this node has valid interpolation data
      if (self%ilon(node_idx) <= 0) return

      ilon = self%ilon(node_idx)
      ilat = self%ilat(node_idx)

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
      this_weight = self%weights(node_idx,:)
      node_result%u = interpolate_bilinear(wind_data%U, this_weight, ilon, ilat)
      node_result%v = interpolate_bilinear(wind_data%V, this_weight, ilon, ilat)
      wind_speed = interpolate_bilinear(wind_data%Ws, this_weight, ilon, ilat)
      node_result%p = interpolate_bilinear(wind_data%P, this_weight, ilon, ilat)

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

end module mod_nws13_data
