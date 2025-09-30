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
module mod_grid_search

   implicit none

   private

   !> @brief Spatial hash acceleration structure for grid cell location
   type t_spatial_hash
      private
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
      procedure, public, pass(self) :: find => find_containing_cell
      final :: destroy_spatial_hash

   end type t_spatial_hash

   !> @brief List of grid cells overlapping a hash cell
   type t_cell_list
      private
      integer :: n_cells = 0 !< Number of grid cells in list
      integer, allocatable :: i_indices(:) !< Grid i-indices
      integer, allocatable :: j_indices(:) !< Grid j-indices

   contains
      procedure, public, pass(self) :: add => add_to_cell_list
   end type t_cell_list

   interface t_cell_list
      module procedure t_cell_list_constructor
   end interface t_cell_list

   interface t_spatial_hash
      module procedure t_spatial_hash_constructor
   end interface t_spatial_hash

   public :: t_spatial_hash, compute_quadrilateral_area, compute_bilinear_weights

contains

   !-----------------------------------------------------------------------
   !> @brief Constructs a new t_cell_list object for spatial hash cells
   !>
   !> @details Creates an empty cell list with initial capacity for storing
   !> grid cell indices. Used internally by the spatial hash acceleration
   !> structure to track which grid cells overlap each hash cell.
   !>
   !> @return cell_list Initialized cell list with pre-allocated arrays
   !-----------------------------------------------------------------------
   pure type(t_cell_list) function t_cell_list_constructor() result(cell_list)
      implicit none
      cell_list%n_cells = 0
      allocate (cell_list%i_indices(16))
      allocate (cell_list%j_indices(16))
   end function t_cell_list_constructor

   !-----------------------------------------------------------------------
   !> @brief Builds spatial hash acceleration structure for grid cell location
   !>
   !> @details Creates a spatial hash structure to accelerate point-in-grid
   !> searches. Divides the spatial domain into hash cells and maps each
   !> grid cell to overlapping hash cells. Hash cell size is set to ~2.5x
   !> the average grid spacing to ensure adequate overlap.
   !>
   !> @param[in] nx         Grid x-dimension
   !> @param[in] ny         Grid y-dimension
   !> @param[in] lon        Grid longitudes (degrees)
   !> @param[in] lat        Grid latitudes (degrees)
   !> @return    hash_struct Initialized spatial hash structure
   !-----------------------------------------------------------------------
   pure type(t_spatial_hash) function t_spatial_hash_constructor(nx, ny, lon, lat) result(hash_struct)
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
   !>
   !> @details Appends a grid cell (i,j) index pair to the cell list.
   !> Automatically resizes the internal arrays when capacity is exceeded
   !> by doubling the array size for efficient growth.
   !>
   !> @param[inout] self Hash cell list to add to
   !> @param[in]    i    Grid i-index
   !> @param[in]    j    Grid j-index
   !-----------------------------------------------------------------------
   pure subroutine add_to_cell_list(self, i, j)
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
   !>
   !> @details Deallocates all arrays associated with the spatial hash
   !> structure including the hash grid and all cell lists. Called
   !> automatically when the hash structure goes out of scope.
   !>
   !> @param[inout] hash_struct Spatial hash to destroy
   !-----------------------------------------------------------------------
   pure subroutine destroy_spatial_hash(hash_struct)
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
   !>
   !> @details Uses a two-stage search algorithm: first checks candidate cells
   !> from the spatial hash, then performs grid walking from the closest
   !> candidate if the initial search fails. Uses area-based containment
   !> test for quadrilateral cells.
   !>
   !> @param[in] self      Spatial hash structure
   !> @param[in] query_lon Query point longitude (degrees)
   !> @param[in] query_lat Query point latitude (degrees)
   !> @param[in] nx        Grid x-dimension
   !> @param[in] ny        Grid y-dimension
   !> @param[in] lon       Grid longitudes (degrees)
   !> @param[in] lat       Grid latitudes (degrees)
   !> @return    indices   [i,j] indices of containing cell, or [-1,-1] if not found
   !-----------------------------------------------------------------------
   pure function find_containing_cell(self, query_lon, query_lat, nx, ny, lon, lat) result(indices)
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

         ! Quick point-in-quad test using areas
         w0 = compute_quadrilateral_area(lon(i, j), lat(i, j), lon(i + 1, j), lat(i + 1, j), &
                                         lon(i + 1, j + 1), lat(i + 1, j + 1), lon(i, j + 1), lat(i, j + 1))
         w1 = compute_quadrilateral_area(query_lon, query_lat, lon(i + 1, j), lat(i + 1, j), &
                                         lon(i + 1, j + 1), lat(i + 1, j + 1), lon(i, j + 1), lat(i, j + 1))
         w2 = compute_quadrilateral_area(lon(i, j), lat(i, j), query_lon, query_lat, &
                                         lon(i + 1, j + 1), lat(i + 1, j + 1), lon(i, j + 1), lat(i, j + 1))
         w3 = compute_quadrilateral_area(lon(i, j), lat(i, j), lon(i + 1, j), lat(i + 1, j), &
                                         query_lon, query_lat, lon(i, j + 1), lat(i, j + 1))
         w4 = compute_quadrilateral_area(lon(i, j), lat(i, j), lon(i + 1, j), lat(i + 1, j), &
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
   !>
   !> @details Performs a directed walk through the grid starting from an
   !> initial cell. At each step, moves toward the query point based on
   !> the relative position to the current cell center. Limited to a maximum
   !> number of iterations to prevent infinite loops.
   !>
   !> @param[in] start_i   Starting i-index
   !> @param[in] start_j   Starting j-index
   !> @param[in] query_lon Query longitude (degrees)
   !> @param[in] query_lat Query latitude (degrees)
   !> @param[in] nx        Grid x-dimension
   !> @param[in] ny        Grid y-dimension
   !> @param[in] lon       Grid longitudes (degrees)
   !> @param[in] lat       Grid latitudes (degrees)
   !> @return    indices   [i,j] of containing cell or [-1,-1]
   !-----------------------------------------------------------------------
   pure function walk_to_containing_cell(start_i, start_j, query_lon, query_lat, nx, ny, lon, lat) result(indices)
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
         w0 = compute_quadrilateral_area(lon(i, j), lat(i, j), lon(i + 1, j), lat(i + 1, j), &
                                         lon(i + 1, j + 1), lat(i + 1, j + 1), lon(i, j + 1), lat(i, j + 1))
         w1 = compute_quadrilateral_area(query_lon, query_lat, lon(i + 1, j), lat(i + 1, j), &
                                         lon(i + 1, j + 1), lat(i + 1, j + 1), lon(i, j + 1), lat(i, j + 1))
         w2 = compute_quadrilateral_area(lon(i, j), lat(i, j), query_lon, query_lat, &
                                         lon(i + 1, j + 1), lat(i + 1, j + 1), lon(i, j + 1), lat(i, j + 1))
         w3 = compute_quadrilateral_area(lon(i, j), lat(i, j), lon(i + 1, j), lat(i + 1, j), &
                                         query_lon, query_lat, lon(i, j + 1), lat(i, j + 1))
         w4 = compute_quadrilateral_area(lon(i, j), lat(i, j), lon(i + 1, j), lat(i + 1, j), &
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
   !> @brief Computes the area of a quadrilateral using the shoelace formula
   !>
   !> @details Calculates the area of a quadrilateral using the shoelace
   !> (surveyor's) formula. Used for point-in-polygon tests and bilinear
   !> interpolation weight calculations. Returns the absolute value of the area.
   !>
   !> @param[in] x1,y1 First vertex coordinates
   !> @param[in] x2,y2 Second vertex coordinates
   !> @param[in] x3,y3 Third vertex coordinates
   !> @param[in] x4,y4 Fourth vertex coordinates
   !> @return    area  Absolute area of the quadrilateral
   !-----------------------------------------------------------------------
   real(8) pure function compute_quadrilateral_area(x1, y1, x2, y2, x3, y3, x4, y4) result(area)
      real(8), intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
      area = abs((x1*y2 - x2*y1) + (x2*y3 - x3*y2) + (x3*y4 - x4*y3) + (x4*y1 - x1*y4))
   end function compute_quadrilateral_area

   !-----------------------------------------------------------------------
   !> @brief Computes bilinear interpolation weights for a quadrilateral cell
   !>
   !> @param[in] query_lon Query point longitude
   !> @param[in] query_lat Query point latitude
   !> @param[in] i         Grid cell i-index
   !> @param[in] j         Grid cell j-index
   !> @param[in] nx        Grid x-dimension
   !> @param[in] ny        Grid y-dimension
   !> @param[in] lon       Grid longitudes (degrees)
   !> @param[in] lat       Grid latitudes (degrees)
   !> @return    weights   [w1,w2,w3,w4] for corners (i,j), (i+1,j), (i+1,j+1), (i,j+1)
   !-----------------------------------------------------------------------
   pure function compute_bilinear_weights(query_lon, query_lat, i, j, nx, ny, lon, lat) result(weights)
      implicit none
      real(8), parameter :: eps = epsilon(1d0)
      real(8), intent(in) :: query_lon, query_lat
      integer, intent(in) :: i, j, nx, ny
      real(8), intent(in) :: lon(nx, ny), lat(nx, ny)
      real(8) :: weights(4)
      real(8) :: total_area

      total_area = compute_quadrilateral_area(lon(i, j), lat(i, j), &
                                              lon(i + 1, j), lat(i + 1, j), &
                                              lon(i + 1, j + 1), lat(i + 1, j + 1), &
                                              lon(i, j + 1), lat(i, j + 1))

      if (total_area <= eps) then
         weights = 0.25d0
      else

         weights(1) = compute_quadrilateral_area(query_lon, query_lat, &
                                                 lon(i + 1, j), lat(i + 1, j), &
                                                 lon(i + 1, j + 1), lat(i + 1, j + 1), &
                                                 lon(i, j + 1), lat(i, j + 1))/(2d0*total_area)

         weights(2) = compute_quadrilateral_area(lon(i, j), lat(i, j), &
                                                 query_lon, query_lat, &
                                                 lon(i + 1, j + 1), lat(i + 1, j + 1), &
                                                 lon(i, j + 1), lat(i, j + 1))/(2d0*total_area)

         weights(3) = compute_quadrilateral_area(lon(i, j), lat(i, j), &
                                                 lon(i + 1, j), lat(i + 1, j), &
                                                 query_lon, query_lat, &
                                                 lon(i, j + 1), lat(i, j + 1))/(2d0*total_area)

         weights(4) = compute_quadrilateral_area(lon(i, j), lat(i, j), &
                                                 lon(i + 1, j), lat(i + 1, j), &
                                                 lon(i + 1, j + 1), lat(i + 1, j + 1), &
                                                 query_lon, query_lat)/(2d0*total_area)

      end if

   end function compute_bilinear_weights

end module mod_grid_search
