!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2025 R.A. Luettich, Jr., J.J. Westerink
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
!===============================================================================
!> \file subgrid.F90
!> \brief Subgrid corrections: lookup table and per-vertex subgrid state.
!>
!> The lookup table (written by the ADCIRC subgrid preprocessor) gives, for each
!> subgrid vertex and each wet-area-fraction level phi_j, the elevation zeta_j at
!> which the level is reached and the averaged quantities at that elevation:
!> grid-averaged depth H_G, wet-averaged depth H_W, friction coefficients (cf at
!> level 0, cmf at level 1) and the advection coefficient C_adv (level 1).
!>
!> The module holds the configuration (namelist &subgridControl), the table and
!> the per-vertex values at the two time levels. The GWCE, momentum, friction
!> and wetting/drying kernels read the per-vertex arrays, which are public and
!> protected; only these procedures write them:
!>   - subgrid_read_namelist, subgrid_validate, subgrid_initialize
!>   - subgrid_update (each time step), subgrid_cold_start, subgrid_hot_start
!>   - subgrid_set_min_depth (drying), subgrid_refresh_vertex (weir
!>     averaging), subgrid_exchange_records (weirs under MPI),
!>     subgrid_finalize
!===============================================================================
#include "logging_macros.h"
module mod_subgrid

   implicit none

   private

   !> Lookup table file name (namelist &subgridControl; adcprep copies the
   !> namelist values into the subdomain fort.15 files)
   character(len=200), public, protected :: subgridFilename = 'null'
   !> Enable the subgrid corrections to wetting/drying, mass and the gravity term
   !> (namelist key level0)
   logical, public, protected :: subgrid_level0 = .false.
   !> Enable the subgrid corrections to bottom friction and advection (namelist
   !> key level1). After subgrid_validate, subgrid_level1 implies subgrid_level0.
   logical, public, protected :: subgrid_level1 = .false.
   !> Include the dphi/dt * dzeta/dt term in the GWCE
   logical, public, protected :: phi_time_derivative = .false.
   !> 1 if phi_time_derivative, 0 otherwise; multiplies the dphi/dt term.
   !> Set by subgrid_read_namelist.
   integer, public, protected :: dphidt = 0
   !> Number of wet-area-fraction levels in the table
   integer :: numPhi = 0

   ! Lookup table (private). The per-vertex arrays are indexed (local
   ! vertex, level).

   !> Wet area fraction of each level (numPhi)
   real(8), allocatable :: setPhi(:)
   !> Friction coefficient (np, numPhi): cmfVertex (with the level 1
   !> correction) if level1, cfVertex otherwise
   real(8), allocatable :: fricVertTab(:, :)
   !> Advection coefficient C_adv (np, numPhi)
   real(8), allocatable :: cadvVertTab(:, :)
   !> Grid-averaged total depth H_G, m (np, numPhi)
   real(8), allocatable :: gridDepthVertTab(:, :)
   !> Wet-averaged total depth H_W, m (np, numPhi)
   real(8), allocatable :: wetDepthVertTab(:, :)
   !> Elevation at which the level is reached, m (np, numPhi)
   real(8), allocatable :: wetFracVertTab(:, :)
   !> Grid-averaged Manning's n, used above the top of the table (np);
   !> allocated only if the table provides it
   real(8), allocatable :: manningAvg(:)

   !> .true. at vertices covered by the table (np)
   logical, allocatable :: isSubgridVert(:)

   ! Per-vertex state at zeta^k (ETA1) and zeta^(k+1) (ETA2). Each record is
   ! the table lookup at the vertex elevation, except after the drying clamp.
   ! subgrid_set_min_depth is depth-only by design: it sets eta and
   ! H_G = H_W = hmin and leaves phi, the friction coefficient and C_adv at
   ! their values for the pre-clamp elevation, so the stored record at a
   ! clamped subgrid vertex is not the table lookup at its eta.
   ! subgrid_refresh_vertex replaces it with the lookup, at the next
   ! subgrid_update or after weir averaging.

   !> Wet area fraction at zeta^k (np)
   real(8), allocatable, public, protected :: wetFracVertETA1(:)
   !> Wet area fraction at zeta^(k+1) (np)
   real(8), allocatable, public, protected :: wetFracVertETA2(:)
   !> Grid-averaged total depth H_G at zeta^k, m (np)
   real(8), allocatable, public, protected :: gridDepthVertETA1(:)
   !> Grid-averaged total depth H_G at zeta^(k+1), m (np)
   real(8), allocatable, public, protected :: gridDepthVertETA2(:)
   !> Wet-averaged total depth H_W at zeta^(k+1), m (np)
   real(8), allocatable, public, protected :: wetDepthVertETA2(:)
   !> Friction coefficient at zeta^(k+1) (np): cmf if level1, cf otherwise.
   !> Read through subgrid_friction.
   real(8), allocatable :: fricVertETA2(:)
   !> Advection coefficient at zeta^k (np), allocated under level1 only
   real(8), allocatable, public, protected :: cadvVertETA1(:)
   !> Advection coefficient at zeta^(k+1) (np), allocated under level1 only
   real(8), allocatable, public, protected :: cadvVertETA2(:)

   !> Table values at one vertex and elevation (result of subgrid_lookup_vertex)
   type :: subgrid_vertex_t
      real(8) :: phi !< wet area fraction
      real(8) :: hg !< grid-averaged total depth H_G, m
      real(8) :: hw !< wet-averaged total depth H_W, m
      real(8) :: fric !< friction coefficient (cmf if level1, cf otherwise)
      real(8) :: cadv = 1.d0 !< advection coefficient C_adv (level1, 1 otherwise)
   end type subgrid_vertex_t

   public :: subgrid_read_namelist
   public :: subgrid_validate
   public :: subgrid_initialize
   public :: subgrid_update
   public :: subgrid_refresh_vertex
   public :: subgrid_exchange_records
   public :: subgrid_cold_start
   public :: subgrid_hot_start
   public :: subgrid_set_min_depth
   public :: subgrid_finalize
   public :: subgrid_has_friction
   public :: subgrid_friction
   public :: subgrid_element_phi

contains

   !=====================================================================
   !> \brief Read the &subgridControl namelist from an open unit.
   !>
   !> A missing namelist (end of file) leaves the defaults in place. Any
   !> other read error, such as a misspelled variable name, stops the run.
   !> Echoes the values, and writes a WARNING if subgridFilename is set
   !> while level0 is off.
   !> \param[in]    unit    unit to read from (positioned by the caller)
   !> \param[out]   ios     iostat of the read
   !> \param[inout] iomsg   iomsg of the read (unchanged if the read
   !>                       succeeds)
   !=====================================================================
   subroutine subgrid_read_namelist(unit, ios, iomsg)
      use global, only: logNamelistReadStatus
      use mod_logging, only: logMessage, allMessage, ECHO, WARNING, t_log_scope, init_log_scope
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none
      integer, intent(in) :: unit
      integer, intent(out) :: ios
      ! allow(C072): the read fills the caller's iomsg buffer, whatever its length
      character(len=*), intent(inout) :: iomsg
      character(len=1024) :: scratchMessage
      logical :: level0, level1
      namelist /subgridControl/ subgridFilename, level0, level1, phi_time_derivative

      LOG_SCOPE_TRACED("subgrid_read_namelist", SUBGRID_TRACING)

      level0 = subgrid_level0
      level1 = subgrid_level1
      read (unit=unit, nml=subgridControl, iostat=ios, iomsg=iomsg)
      call logNamelistReadStatus('subgridControl', ios, iomsg)
      if (ios > 0) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message="Could not read the subgridControl namelist in fort.15: "//trim(iomsg))
      end if
      subgrid_level0 = level0
      subgrid_level1 = level1
      dphidt = merge(1, 0, phi_time_derivative)
      write (scratchMessage, *) "subgridFilename=", trim(adjustl(subgridFilename))
      call logMessage(ECHO, trim(scratchMessage))
      write (scratchMessage, *) "level0=", level0
      call logMessage(ECHO, trim(scratchMessage))
      write (scratchMessage, *) "level1=", level1
      call logMessage(ECHO, trim(scratchMessage))
      write (scratchMessage, *) "phi_time_derivative=", phi_time_derivative
      call logMessage(ECHO, trim(scratchMessage))
      if (subgridFilename /= 'null' .and. .not. subgrid_level0) then
         call allMessage(WARNING, "subgridControl: subgridFilename='"//trim(subgridFilename)// &
                         "' is set but level0=.false.; subgrid corrections are off. "// &
                         "Set level0=.true. to use the table.")
      end if
   end subroutine subgrid_read_namelist

   !=====================================================================
   !> \brief Stop at startup if the subgrid configuration is not supported.
   !>
   !> READ_INPUT calls it once at its end, before the cold start or the
   !> hot-start read.
   !> \param[in] nolifa              finite amplitude and wetting/drying option
   !> \param[in] predictor_corrector .true. if fort.15 gives DTDP < 0
   !> \param[in] load_mannings_n     .true. if mannings_n_at_sea_floor is loaded
   !> \param[in] nddt                time-varying bathymetry option (must be 0)
   !=====================================================================
   subroutine subgrid_validate(nolifa, predictor_corrector, load_mannings_n, nddt)
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      use mod_logging, only: t_log_scope, init_log_scope
      implicit none
      integer, intent(in) :: nolifa
      logical, intent(in) :: predictor_corrector
      logical, intent(in) :: load_mannings_n
      integer, intent(in) :: nddt

      LOG_SCOPE_TRACED("subgrid_validate", SUBGRID_TRACING)

      if (subgrid_level1 .and. .not. subgrid_level0) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message="subgridControl: level1 requires level0; set level0=.true.")
      end if
      if (.not. subgrid_level0) return
#ifndef ADCNETCDF
      call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                     message='Subgrid calculation requires building with netcdf enabled')
#endif
      if (nolifa /= 2) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message="subgrid requires NOLIFA=2 in fort.15.")
      end if
      if (predictor_corrector) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message="subgrid does not support the "// &
                        "predictor-corrector algorithm "// &
                        "(DTDP < 0); set DTDP > 0.")
      end if
      if (.not. load_mannings_n) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message='subgrid requires the mannings_n_at_sea_floor nodal attribute; add it to fort.13 and fort.15.')
      end if
      if (nddt /= 0) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message="subgrid does not support time-varying bathymetry "// &
                        "(NDDT /= 0); set NDDT = 0.")
      end if
   end subroutine subgrid_validate

   !=====================================================================
   !> \brief Read the lookup table on the first call and allocate the
   !>        per-vertex arrays.
   !=====================================================================
   subroutine subgrid_initialize()
      use mesh, only: np
      use mod_logging, only: t_log_scope, init_log_scope
      implicit none

      LOG_SCOPE_TRACED("subgrid_initialize", SUBGRID_TRACING)

      if (.not. allocated(setPhi)) call read_table()
      call allocate_state(np)
   end subroutine subgrid_initialize

   !=====================================================================
   !> \brief Allocate the per-vertex arrays and set their initial values
   !>        (0, and 1 for the advection coefficient).
   !> \param[in] n number of local vertices
   !=====================================================================
   subroutine allocate_state(n)
      use mod_logging, only: t_log_scope, init_log_scope
      implicit none
      integer, intent(in) :: n

      LOG_SCOPE_TRACED("subgrid_allocate_state", SUBGRID_TRACING)

      if (allocated(wetFracVertETA1)) deallocate (wetFracVertETA1)
      allocate (wetFracVertETA1(n))
      if (allocated(gridDepthVertETA1)) deallocate (gridDepthVertETA1)
      allocate (gridDepthVertETA1(n))
      if (allocated(wetFracVertETA2)) deallocate (wetFracVertETA2)
      allocate (wetFracVertETA2(n))
      if (allocated(gridDepthVertETA2)) deallocate (gridDepthVertETA2)
      allocate (gridDepthVertETA2(n))
      if (allocated(wetDepthVertETA2)) deallocate (wetDepthVertETA2)
      allocate (wetDepthVertETA2(n))
      wetFracVertETA1(:) = 0.d0
      wetFracVertETA2(:) = 0.d0
      gridDepthVertETA1(:) = 0.d0
      gridDepthVertETA2(:) = 0.d0
      wetDepthVertETA2(:) = 0.d0
      if (allocated(fricVertETA2)) deallocate (fricVertETA2)
      allocate (fricVertETA2(n))
      fricVertETA2(:) = 0.d0
      if (subgrid_level1) then
         if (allocated(cadvVertETA2)) deallocate (cadvVertETA2)
         allocate (cadvVertETA2(n))
         if (allocated(cadvVertETA1)) deallocate (cadvVertETA1)
         allocate (cadvVertETA1(n))
         ! cadv is 1 (no correction) at vertices outside the subgrid area;
         ! subgrid_refresh_vertex keeps it at 1 there
         cadvVertETA2(:) = 1.d0
         cadvVertETA1(:) = 1.d0
      end if
   end subroutine allocate_state

   !=====================================================================
   !> \brief Read the NetCDF lookup table file subgridFilename.
   !>
   !> Sets numPhi, setPhi, the per-vertex table arrays (cf or cmf and cadv
   !> depending on level1, H_G, H_W, the level elevations), isSubgridVert
   !> and, if present, manningAvg. In parallel, nodes_lg maps the global
   !> vertex numbers to local ones. Stops with the file and variable names if
   !> the file is missing or a read fails.
   !=====================================================================
   subroutine read_table()
#ifdef ADCNETCDF
      use netcdf, only: nf90_open, nf90_close, nf90_inq_dimid, nf90_inquire_dimension, &
                        nf90_inq_varid, nf90_get_var, nf90_strerror, NF90_NOWRITE, NF90_NOERR
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      use mod_logging, only: allMessage, logMessage, WARNING, INFO, t_log_scope, init_log_scope
      use mesh, only: np
      use global, only: np_g
      implicit none
      integer :: numVerts
      integer :: nc_id, nc_var, nc_err
      logical :: file_exists
      real(8), allocatable :: globalVertList(:)
      integer, allocatable :: globalSubVertList(:)
      integer, allocatable :: lg(:)
      character(len=1024) :: scratchMessage

      LOG_SCOPE_TRACED("subgrid_read_table", SUBGRID_TRACING)

      inquire (file=trim(subgridFilename), exist=file_exists)
      if (.not. file_exists) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message="The subgrid table '"//trim(subgridFilename)//"' does not exist; check subgridFilename.")
      end if

      nc_err = nf90_open(trim(subgridFilename), NF90_NOWRITE, nc_id)
      if (nc_err /= NF90_NOERR) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message="Could not open the subgrid table '"//trim(subgridFilename)// &
                        "': "//trim(nf90_strerror(nc_err)))
      end if
      call check_nc(nf90_inq_dimid(nc_id, "numNode", nc_var), "numNode")
      call check_nc(nf90_inquire_dimension(nc_id, nc_var, len=numVerts), "numNode")
      if (numVerts /= np_g) then
         write (scratchMessage, '(a,a,a,i0,a,i0,a)') "The subgrid table '", trim(subgridFilename), &
            "' has ", numVerts, " vertices but the mesh has ", np_g, &
            "; regenerate the table for this fort.14."
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(scratchMessage))
      end if
      call check_nc(nf90_inq_dimid(nc_id, "numPhi", nc_var), "numPhi")
      call check_nc(nf90_inquire_dimension(nc_id, nc_var, len=numPhi), "numPhi")
      if (numPhi < 2) then
         write (scratchMessage, '(a,a,a,i0,a)') "The subgrid table '", trim(subgridFilename), &
            "' has ", numPhi, " phi levels; regenerate the table with at least 2."
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(scratchMessage))
      end if

      ! global vertex number of each local vertex
      allocate (lg(np))
      lg(:) = global_node(np)

      if (allocated(setPhi)) deallocate (setPhi)
      allocate (setPhi(numPhi))
      call check_nc(nf90_inq_varid(nc_id, "phiSet", nc_var), "phiSet")
      call check_nc(nf90_get_var(nc_id, nc_var, setPhi), "phiSet")

      if (subgrid_level1) then
         call read_vertex_table(nc_id, "cmfVertex", numVerts, lg, fricVertTab)
         call read_vertex_table(nc_id, "cadvVertex", numVerts, lg, cadvVertTab)
      else
         call read_vertex_table(nc_id, "cfVertex", numVerts, lg, fricVertTab)
      end if
      call read_vertex_table(nc_id, "gridTotWatDepthVertex", numVerts, lg, gridDepthVertTab)
      call read_vertex_table(nc_id, "wetTotWatDepthVertex", numVerts, lg, wetDepthVertTab)
      call read_vertex_table(nc_id, "wetFractionVertex", numVerts, lg, wetFracVertTab)

      allocate (globalSubVertList(numVerts))
      call check_nc(nf90_inq_varid(nc_id, "binaryVertexList", nc_var), "binaryVertexList")
      call check_nc(nf90_get_var(nc_id, nc_var, globalSubVertList), "binaryVertexList")
      if (allocated(isSubgridVert)) deallocate (isSubgridVert)
      allocate (isSubgridVert(np))
      isSubgridVert(:) = globalSubVertList(lg) == 1
      deallocate (globalSubVertList)

      ! grid-averaged Manning's n is optional
      if (allocated(manningAvg)) deallocate (manningAvg)
      if (nf90_inq_varid(nc_id, "manningAvg", nc_var) == NF90_NOERR) then
         allocate (globalVertList(numVerts))
         call check_nc(nf90_get_var(nc_id, nc_var, globalVertList), "manningAvg")
         allocate (manningAvg(np))
         manningAvg(:) = globalVertList(lg)
         deallocate (globalVertList)
         call logMessage(INFO, "The subgrid table provides manningAvg.")
      else
         call allMessage(WARNING, "subgrid table has no manningAvg; "// &
                         "above the top level the friction stays at the top-level value")
      end if
      call validate_table(lg)
      if (subgrid_level1) call warn_cadv_below_one(lg)

      nc_err = nf90_close(nc_id)
      write (scratchMessage, '(a,a,a,i0,a)') "Read the subgrid table '", trim(subgridFilename), &
         "' with ", numPhi, " phi levels."
      call logMessage(INFO, trim(scratchMessage))
#ifdef CMPI
      write (scratchMessage, '(i0,a)') count(isSubgridVert), &
         " subgrid vertices in this subdomain (local count)."
#else
      write (scratchMessage, '(i0,a)') count(isSubgridVert), " subgrid vertices."
#endif
      call logMessage(INFO, trim(scratchMessage))
#else
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none
      call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                     message='Subgrid calculation requires building with netcdf enabled')
#endif
   end subroutine read_table

#ifdef ADCNETCDF
   !=====================================================================
   !> \brief Stop if the table read by read_table cannot be used.
   !>
   !> The lookups need finite, strictly increasing phi levels and, at every subgrid
   !> vertex, finite and non-decreasing elevations (wetFractionVertex) and
   !> grid-averaged depths (gridTotWatDepthVertex). Equal neighbouring
   !> values are allowed; the interval that the forward and inverse lookups
   !> find still has a nonzero width. The other columns read at subgrid
   !> vertices (wet-averaged depth, friction and, under level1, cadv) must
   !> be finite, in any order. A magnitude of 1e30 or more (the fill value
   !> of the table generator) fails the check in all of these columns.
   !> manningAvg, if present, must be finite and positive at subgrid
   !> vertices. The finiteness test comes before each comparison, so a NaN
   !> cannot raise a floating-point trap before the message is written.
   !> \param[in] lg global vertex number of each local vertex (np)
   !=====================================================================
   subroutine validate_table(lg)
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      use mod_logging, only: t_log_scope, init_log_scope
      implicit none
      integer, intent(in) :: lg(:)
      integer :: i
      logical :: valid
      character(len=9) :: fricName
      character(len=1024) :: scratchMessage

      LOG_SCOPE_TRACED("subgrid_validate_table", SUBGRID_TRACING)

      valid = all(ieee_is_finite(setPhi))
      if (valid) valid = all(setPhi(2:numPhi) > setPhi(1:numPhi - 1))
      if (.not. valid) then
         write (scratchMessage, '(a,a,a)') "The phiSet values in the subgrid table '", &
            trim(subgridFilename), "' are not finite and strictly increasing; regenerate the table."
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(scratchMessage))
      end if
      if (subgrid_level1) then
         fricName = "cmfVertex"
      else
         fricName = "cfVertex"
      end if
      do i = 1, size(isSubgridVert)
         if (.not. isSubgridVert(i)) cycle
         if (.not. valid_row(wetFracVertTab(i, :))) call bad_row(lg(i), "wetFractionVertex")
         if (.not. valid_row(gridDepthVertTab(i, :))) call bad_row(lg(i), "gridTotWatDepthVertex")
         if (.not. finite_row(wetDepthVertTab(i, :))) call bad_value(lg(i), "wetTotWatDepthVertex")
         if (.not. finite_row(fricVertTab(i, :))) call bad_value(lg(i), trim(fricName))
         if (subgrid_level1) then
            if (.not. finite_row(cadvVertTab(i, :))) call bad_value(lg(i), "cadvVertex")
         end if
         if (allocated(manningAvg)) then
            valid = ieee_is_finite(manningAvg(i))
            if (valid) valid = manningAvg(i) > 0.d0
            if (.not. valid) then
               write (scratchMessage, '(a,a,a,i0,a)') "The subgrid table '", trim(subgridFilename), &
                  "' has a non-finite or non-positive manningAvg at node ", lg(i), "; regenerate the table."
               call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(scratchMessage))
            end if
         end if
      end do

   contains

      !> .true. if every value of the row is finite and of magnitude below 1e30
      pure logical function finite_row(row) result(valid)
         implicit none
         real(8), intent(in) :: row(:)

         valid = all(ieee_is_finite(row))
         if (valid) valid = all(abs(row) < 1.d30)
      end function finite_row

      !> .true. if finite_row(row) and the values do not decrease along the row
      pure logical function valid_row(row) result(valid)
         implicit none
         real(8), intent(in) :: row(:)
         integer :: n

         n = size(row)
         valid = finite_row(row)
         if (valid) valid = all(row(2:n) >= row(1:n - 1))
      end function valid_row

      subroutine bad_row(node, name)
         implicit none
         integer, intent(in) :: node
         character(len=*), intent(in) :: name

         write (scratchMessage, '(a,a,a,a,a,i0,a)') "The subgrid table '", trim(subgridFilename), &
            "' has a non-finite, fill-valued (|x| >= 1e30) or decreasing ", name, &
            " row at node ", node, "; regenerate the table."
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(scratchMessage))
      end subroutine bad_row

      subroutine bad_value(node, name)
         implicit none
         integer, intent(in) :: node
         character(len=*), intent(in) :: name

         write (scratchMessage, '(a,a,a,a,a,i0,a)') "The subgrid table '", trim(subgridFilename), &
            "' has a non-finite or fill-valued (|x| >= 1e30) ", name, &
            " row at node ", node, "; regenerate the table."
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, message=trim(scratchMessage))
      end subroutine bad_value
   end subroutine validate_table

   !=====================================================================
   !> \brief Warn if the advection coefficient C_adv is below 1 anywhere in
   !>        the table rows of the subgrid vertices (level1).
   !>
   !> A consistent table has C_adv >= 1, so values below 1 point to a
   !> defect in the table generator. The run continues with the table as
   !> read. Call after validate_table, which ensures the rows are finite.
   !> \param[in] lg global vertex number of each local vertex (np)
   !=====================================================================
   subroutine warn_cadv_below_one(lg)
      use mod_logging, only: allMessage, WARNING, t_log_scope, init_log_scope
      implicit none
      integer, intent(in) :: lg(:)
      integer :: i
      integer :: numBelow
      integer :: worstNode
      real(8) :: rowMin
      real(8) :: worst
      character(len=1024) :: scratchMessage

      LOG_SCOPE_TRACED("subgrid_warn_cadv_below_one", SUBGRID_TRACING)

      numBelow = 0
      worst = 1.d0
      worstNode = 0
      do i = 1, size(isSubgridVert)
         if (.not. isSubgridVert(i)) cycle
         rowMin = minval(cadvVertTab(i, :))
         if (rowMin < 1.d0) then
            numBelow = numBelow + 1
            if (rowMin < worst) then
               worst = rowMin
               worstNode = lg(i)
            end if
         end if
      end do
      if (numBelow == 0) return
#ifdef CMPI
      write (scratchMessage, '(a,a,a,i0,a,g0.8,a,i0,a)') "The subgrid table '", trim(subgridFilename), &
         "' has cadvVertex < 1 at ", numBelow, " subgrid vertices in this subdomain (minimum ", &
         worst, " at node ", worstNode, "); a consistent table has C_adv >= 1."
#else
      write (scratchMessage, '(a,a,a,i0,a,g0.8,a,i0,a)') "The subgrid table '", trim(subgridFilename), &
         "' has cadvVertex < 1 at ", numBelow, " subgrid vertices (minimum ", &
         worst, " at node ", worstNode, "); a consistent table has C_adv >= 1."
#endif
      call allMessage(WARNING, trim(scratchMessage))
   end subroutine warn_cadv_below_one

   !=====================================================================
   !> \brief Read one (numPhi, numNode) table variable into tab(n, numPhi).
   !> \param[in]    nc_id    NetCDF id of the open table file
   !> \param[in]    name     name of the variable
   !> \param[in]    numVerts number of vertices in the file (global)
   !> \param[in]    lg       global vertex number of each local vertex (n)
   !> \param[inout] tab      table, (re)allocated to (n, numPhi)
   !=====================================================================
   subroutine read_vertex_table(nc_id, name, numVerts, lg, tab)
      use netcdf, only: nf90_inq_varid, nf90_get_var
      use mod_logging, only: t_log_scope, init_log_scope
      implicit none
      integer, intent(in) :: nc_id
      character(len=*), intent(in) :: name
      integer, intent(in) :: numVerts
      integer, intent(in) :: lg(:)
      real(8), allocatable, intent(inout) :: tab(:, :)
      real(8), allocatable :: globalVertLookup(:, :)
      integer :: nc_var
      integer :: i

      LOG_SCOPE_TRACED("subgrid_read_vertex_table", SUBGRID_TRACING)

      allocate (globalVertLookup(numPhi, numVerts))
      call check_nc(nf90_inq_varid(nc_id, name, nc_var), name)
      call check_nc(nf90_get_var(nc_id, nc_var, globalVertLookup), name)
      if (allocated(tab)) deallocate (tab)
      allocate (tab(size(lg), numPhi))
      do i = 1, size(lg)
         tab(i, :) = globalVertLookup(:, lg(i))
      end do
      deallocate (globalVertLookup)
   end subroutine read_vertex_table

   !=====================================================================
   !> \brief Stop with the table file name, the item being read and the
   !>        NetCDF message if a NetCDF call failed.
   !> \param[in] status return value of the NetCDF call
   !> \param[in] item   variable or dimension being read
   !=====================================================================
   subroutine check_nc(status, item)
      use netcdf, only: nf90_strerror, NF90_NOERR
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      implicit none
      integer, intent(in) :: status
      character(len=*), intent(in) :: item

      if (status == NF90_NOERR) return
      call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                     message="Error reading '"//item//"' from the subgrid table '"// &
                     trim(subgridFilename)//"': "//trim(nf90_strerror(status)))
   end subroutine check_nc

   !=====================================================================
   !> \brief Global vertex number of each local vertex (identity in serial).
   !>
   !> nodes_lg is negative at ghost vertices; abs gives the global number.
   !> \param[in] n number of local vertices
   !> \return global vertex numbers (n)
   !=====================================================================
   pure function global_node(n) result(lg)
#ifdef CMPI
      use global, only: nodes_lg
#endif
      implicit none
      integer, intent(in) :: n
      integer :: lg(n)
      integer :: i

      do i = 1, n
#ifdef CMPI
         lg(i) = abs(nodes_lg(i))
#else
         lg(i) = i
#endif
      end do
   end function global_node
#endif

   !=====================================================================
   !> \brief Linear interpolation between (x0, y0) and (x1, y1).
   !>
   !> Requires x1 /= x0; validate_table ensures this for the two table
   !> levels that bracket x.
   !> \param[in] x  abscissa to interpolate at
   !> \param[in] x0 lower abscissa
   !> \param[in] x1 upper abscissa
   !> \param[in] y0 value at x0
   !> \param[in] y1 value at x1
   !> \return value at x
   !=====================================================================
   elemental real(8) function subgrid_interpolate(x, x0, x1, y0, y1) result(y)
      implicit none
      real(8), intent(in) :: x
      real(8), intent(in) :: x0
      real(8), intent(in) :: x1
      real(8), intent(in) :: y0
      real(8), intent(in) :: y1

      y = ((x - x0)/(x1 - x0)*(y1 - y0) + y0)
   end function subgrid_interpolate

   !=====================================================================
   !> \brief Table values at subgrid vertex i for elevation eta.
   !>
   !> Below the lowest level the result is the lowest level. Above the top
   !> level it is the top level with H_G and H_W increased by the height of
   !> eta above that level and, if the table has manningAvg, the friction
   !> coefficient computed from Manning's n at that H_G. In between, it
   !> interpolates linearly between the two levels that bracket eta. A
   !> degenerate (equal or decreasing) interval gives the upper level;
   !> validate_table rejects tables that have one.
   !> \param[in] i   local vertex number (isSubgridVert(i))
   !> \param[in] eta elevation, m
   !> \return values at eta: wet area fraction (phi), grid- and wet-averaged
   !>         total depths in m (hg, hw), friction coefficient (fric: cmf if
   !>         level1, cf otherwise) and advection coefficient (cadv: level1,
   !>         1 otherwise)
   !=====================================================================
   pure function subgrid_lookup_vertex(i, eta) result(v)
      use adc_constants, only: g
      implicit none
      integer, intent(in) :: i
      real(8), intent(in) :: eta
      type(subgrid_vertex_t) :: v
      integer :: k, numGreater

      ! number of levels whose elevation is above eta
      numGreater = count(eta < wetFracVertTab(i, :))
      if (numGreater == numPhi) then
         ! dry: below the lowest level
         v = table_row(i, 1)
      elseif (numGreater == 0) then
         ! above the top level: depths grow with eta
         v = table_row(i, numPhi)
         v%hg = v%hg + (eta - wetFracVertTab(i, numPhi))
         v%hw = v%hw + (eta - wetFracVertTab(i, numPhi))
         if (allocated(manningAvg)) then
            v%fric = g*manningAvg(i)**2.d0/((v%hg)**(1.d0/3.d0))
         end if
      else
         k = numPhi - numGreater
         if (wetFracVertTab(i, k + 1) > wetFracVertTab(i, k)) then
            ! eta lies between the elevations of levels k and k+1
            v%phi = subgrid_interpolate(eta, wetFracVertTab(i, k), wetFracVertTab(i, k + 1), &
                                        setPhi(k), setPhi(k + 1))
            v%hg = subgrid_interpolate(eta, wetFracVertTab(i, k), wetFracVertTab(i, k + 1), &
                                       gridDepthVertTab(i, k), gridDepthVertTab(i, k + 1))
            v%hw = subgrid_interpolate(eta, wetFracVertTab(i, k), wetFracVertTab(i, k + 1), &
                                       wetDepthVertTab(i, k), wetDepthVertTab(i, k + 1))
            v%fric = subgrid_interpolate(eta, wetFracVertTab(i, k), wetFracVertTab(i, k + 1), &
                                         fricVertTab(i, k), fricVertTab(i, k + 1))
            if (subgrid_level1) then
               v%cadv = subgrid_interpolate(eta, wetFracVertTab(i, k), wetFracVertTab(i, k + 1), &
                                            cadvVertTab(i, k), cadvVertTab(i, k + 1))
            end if
         else
            ! degenerate interval: use the upper level (validate_table
            ! rejects tables that have one)
            v = table_row(i, k + 1)
         end if
      end if
   end function subgrid_lookup_vertex

   !=====================================================================
   !> \brief Table values of level j at subgrid vertex i.
   !> \param[in] i local vertex number (isSubgridVert(i))
   !> \param[in] j level number (1..numPhi)
   !> \return the table row (cadv = 1 unless level1)
   !=====================================================================
   pure function table_row(i, j) result(v)
      implicit none
      integer, intent(in) :: i
      integer, intent(in) :: j
      type(subgrid_vertex_t) :: v

      v%phi = setPhi(j)
      v%hg = gridDepthVertTab(i, j)
      v%hw = wetDepthVertTab(i, j)
      v%fric = fricVertTab(i, j)
      if (subgrid_level1) v%cadv = cadvVertTab(i, j)
   end function table_row

   !=====================================================================
   !> \brief Shift the zeta^(k+1) values to zeta^k, then look up all vertices
   !>        at the new elevation eta2.
   !> \param[in] eta2 elevation at every local vertex, m (np)
   !=====================================================================
   subroutine subgrid_update(eta2)
      use mesh, only: np
      use mod_logging, only: t_log_scope, init_log_scope
      implicit none
      real(8), intent(in) :: eta2(:)
      integer :: i

      LOG_SCOPE_TRACED("subgrid_update", SUBGRID_TRACING)

      call shift_time_levels()
      do i = 1, np
         call subgrid_refresh_vertex(i, eta2(i))
      end do
   end subroutine subgrid_update

   !=====================================================================
   !> \brief Look up vertex i at elevation eta and store the result as its
   !>        zeta^(k+1) values, without shifting the time levels.
   !>
   !> subgrid_update calls it for every vertex, and wetting/drying calls it
   !> for the two vertices of a weir pair after averaging their elevations.
   !> It replaces the record that subgrid_set_min_depth leaves at a clamped
   !> vertex with the lookup at eta.
   !>
   !> Precondition: level0, and subgrid_initialize has been called.
   !> \param[in] i   local vertex number
   !> \param[in] eta elevation at vertex i, m
   !=====================================================================
   subroutine subgrid_refresh_vertex(i, eta)
      ! no log scope: called per vertex
      use mesh, only: dp
      implicit none
      integer, intent(in) :: i
      real(8), intent(in) :: eta
      type(subgrid_vertex_t) :: v

      if (isSubgridVert(i)) then
         v = subgrid_lookup_vertex(i, eta)
         wetFracVertETA2(i) = v%phi
         gridDepthVertETA2(i) = v%hg
         wetDepthVertETA2(i) = v%hw
         fricVertETA2(i) = v%fric
         if (subgrid_level1) cadvVertETA2(i) = v%cadv
      else
         ! not in the subgrid area: grid total depth
         gridDepthVertETA2(i) = eta + dp(i)
         wetDepthVertETA2(i) = eta + dp(i)
         wetFracVertETA2(i) = 1.d0
         ! no advection correction (the kernels multiply by cadv)
         if (subgrid_level1) cadvVertETA2(i) = 1.d0
      end if
   end subroutine subgrid_refresh_vertex

   !=====================================================================
   !> \brief Give every ghost vertex the zeta^(k+1) values of its owner
   !>        (MPI; no-op in serial).
   !>
   !> Wetting/drying averages the elevations of the two vertices of an
   !> IBTYPE=64 weir pair and refreshes their values on each subdomain
   !> that holds both vertices. When the pair is split between subdomains,
   !> only some subdomains do this. The ETA2 exchange that follows gives
   !> each ghost its owner's elevation, and this exchange gives it the
   !> owner's values. The ghost copies the record instead of repeating the
   !> lookup because the record at a clamped vertex (see
   !> subgrid_set_min_depth) differs from the lookup at its elevation.
   !>
   !> Precondition: level0, and subgrid_initialize has been called.
   !=====================================================================
   subroutine subgrid_exchange_records()
#ifdef CMPI
      use messenger, only: updater
      use global, only: dumy1, dumy2
#endif
      use mod_logging, only: t_log_scope, init_log_scope
      implicit none

      LOG_SCOPE_TRACED("subgrid_exchange_records", SUBGRID_TRACING)

#ifdef CMPI
      call updater(gridDepthVertETA2, wetDepthVertETA2, wetFracVertETA2, 3)
      if (subgrid_level1) then
         call updater(fricVertETA2, cadvVertETA2, dumy1, 2)
      else
         call updater(fricVertETA2, dumy1, dumy2, 1)
      end if
#endif
   end subroutine subgrid_exchange_records

   !=====================================================================
   !> \brief Copy the zeta^(k+1) values that have a zeta^k copy (wet area
   !>        fraction, H_G and, under level1, cadv) to zeta^k.
   !=====================================================================
   subroutine shift_time_levels()
      implicit none

      wetFracVertETA1 = wetFracVertETA2
      gridDepthVertETA1 = gridDepthVertETA2
      if (subgrid_level1) cadvVertETA1 = cadvVertETA2
   end subroutine shift_time_levels

   !=====================================================================
   !> \brief Cold start: raise the dry vertices to the minimum depth and set
   !>        the previous-level values for the first time step.
   !>
   !> Preconditions: subgrid_initialize and subgrid_update(eta2) have been
   !> called with the initial elevation, and nnodecode holds the initial
   !> wet/dry state.
   !> \param[inout] eta2      elevation, m, reset at dry vertices with H_G < hmin (np)
   !> \param[in]    nnodecode wet (1) / dry (0) flag of each vertex (np)
   !> \param[in]    hmin      minimum total depth (0.8*H0), m
   !=====================================================================
   subroutine subgrid_cold_start(eta2, nnodecode, hmin)
      use mesh, only: np
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      use mod_logging, only: logMessage, INFO, t_log_scope, init_log_scope
      implicit none
      real(8), intent(inout) :: eta2(:)
      integer, intent(in) :: nnodecode(:)
      real(8), intent(in) :: hmin
      integer :: i

      LOG_SCOPE_TRACED("subgrid_cold_start", SUBGRID_TRACING)

      if (size(eta2) < np .or. size(nnodecode) < np) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message="subgrid_cold_start: eta2 and nnodecode must have np entries.")
      end if
      do i = 1, np
         if (nnodecode(i) == 0) then
            if (gridDepthVertETA2(i) < hmin) then
               call subgrid_set_min_depth(i, hmin, eta2(i))
            end if
         end if
      end do
      ! the first time step uses the initial values at zeta^k
      call shift_time_levels()
      call logMessage(INFO, "Subgrid cold start: previous-level values copied from the initial ones.")
   end subroutine subgrid_cold_start

   !=====================================================================
   !> \brief Hot start: previous-level values from eta1, current from eta2.
   !>
   !> Precondition: subgrid_initialize has been called. eta1 and eta2 are
   !> the elevations read from the hot-start file.
   !> \param[in] eta1 elevation at the previous time level, m (np)
   !> \param[in] eta2 elevation at the current time level, m (np)
   !=====================================================================
   subroutine subgrid_hot_start(eta1, eta2)
      use mod_logging, only: logMessage, INFO, t_log_scope, init_log_scope
      implicit none
      real(8), intent(in) :: eta1(:)
      real(8), intent(in) :: eta2(:)

      LOG_SCOPE_TRACED("subgrid_hot_start", SUBGRID_TRACING)

      call subgrid_update(eta1)
      ! the second update shifts the eta1 lookups to the previous level
      call subgrid_update(eta2)
      call logMessage(INFO, "Subgrid hot start: lookups set from ETA1 (previous level) and ETA2.")
   end subroutine subgrid_hot_start

   !=====================================================================
   !> \brief Clamp vertex i to the minimum total depth hmin (drying).
   !>
   !> Sets H_G = H_W = hmin and returns in eta the elevation at which the
   !> vertex has H_G = hmin (inverse table lookup; hmin - dp outside the
   !> subgrid area). The clamp is depth-only by design: phi, the friction
   !> coefficient and C_adv keep their values at the pre-clamp elevation,
   !> so the stored record at a clamped subgrid vertex is not the table
   !> lookup at its eta. subgrid_refresh_vertex replaces it with the
   !> lookup, at the next subgrid_update or after weir averaging.
   !>
   !> Precondition: level0, and subgrid_initialize has been called.
   !> \param[in]  i    local vertex number
   !> \param[in]  hmin minimum total depth, m
   !> \param[out] eta  elevation at which vertex i has H_G = hmin, m
   !=====================================================================
   subroutine subgrid_set_min_depth(i, hmin, eta)
      ! no log scope: called per vertex inside the wetting/drying loops
      use mesh, only: dp
      implicit none
      integer, intent(in) :: i
      real(8), intent(in) :: hmin
      real(8), intent(out) :: eta

      gridDepthVertETA2(i) = hmin
      wetDepthVertETA2(i) = hmin
      if (isSubgridVert(i)) then
         eta = subgrid_eta_at_depth(i, hmin)
      else
         eta = hmin - dp(i)
      end if
   end subroutine subgrid_set_min_depth

   !=====================================================================
   !> \brief Elevation at which subgrid vertex i reaches grid-averaged depth h
   !>        (inverse table lookup).
   !>
   !> validate_table guarantees non-decreasing depths, so the bracketing
   !> interval in the interpolation has nonzero width.
   !> \param[in] i local vertex number (isSubgridVert(i))
   !> \param[in] h grid-averaged total depth, m
   !> \return elevation, m
   !=====================================================================
   pure real(8) function subgrid_eta_at_depth(i, h) result(eta)
      implicit none
      integer, intent(in) :: i
      real(8), intent(in) :: h
      integer :: k, numGreater

      numGreater = count(h < gridDepthVertTab(i, :))
      if (numGreater == numPhi) then
         eta = wetFracVertTab(i, 1)
      elseif (numGreater == 0) then
         eta = wetFracVertTab(i, numPhi) + (h - gridDepthVertTab(i, numPhi))
      else
         k = numPhi - numGreater
         eta = subgrid_interpolate(h, gridDepthVertTab(i, k), gridDepthVertTab(i, k + 1), &
                                   wetFracVertTab(i, k), wetFracVertTab(i, k + 1))
      end if
   end function subgrid_eta_at_depth

   !=====================================================================
   !> \brief .true. if the friction at vertex i comes from the subgrid table
   !>        (level0 on and i is a subgrid vertex).
   !> \param[in] i local vertex number
   !> \return .true. if subgrid_friction(i) applies at vertex i
   !=====================================================================
   pure logical function subgrid_has_friction(i) result(has)
      implicit none
      integer, intent(in) :: i

      ! two statements: Fortran does not guarantee that .and. short-circuits,
      ! and isSubgridVert is not allocated when level0 is off
      has = .false.
      if (subgrid_level0) has = isSubgridVert(i)
   end function subgrid_has_friction

   !=====================================================================
   !> \brief Stored zeta^(k+1) friction coefficient at vertex i (cmf if
   !>        level1, cf otherwise). Valid only where subgrid_has_friction(i).
   !> \param[in] i local vertex number
   !> \return friction coefficient
   !=====================================================================
   pure real(8) function subgrid_friction(i) result(f)
      implicit none
      integer, intent(in) :: i

      f = fricVertETA2(i)
   end function subgrid_friction

   !=====================================================================
   !> \brief Element means of the wet area fraction at zeta^k and zeta^(k+1).
   !> \param[in]  n1      first vertex of the element
   !> \param[in]  n2      second vertex of the element
   !> \param[in]  n3      third vertex of the element
   !> \param[out] phiavg1 mean wet area fraction at zeta^k
   !> \param[out] phiavg2 mean wet area fraction at zeta^(k+1)
   !=====================================================================
   pure subroutine subgrid_element_phi(n1, n2, n3, phiavg1, phiavg2)
      implicit none
      integer, intent(in) :: n1
      integer, intent(in) :: n2
      integer, intent(in) :: n3
      real(8), intent(out) :: phiavg1
      real(8), intent(out) :: phiavg2

      phiavg1 = (wetFracVertETA1(n1) + wetFracVertETA1(n2) + wetFracVertETA1(n3))/3.d0
      phiavg2 = (wetFracVertETA2(n1) + wetFracVertETA2(n2) + wetFracVertETA2(n3))/3.d0
   end subroutine subgrid_element_phi

   !=====================================================================
   !> \brief Release the table and the per-vertex arrays, and reset the
   !>        number of table levels.
   !=====================================================================
   subroutine subgrid_finalize()
      use mod_logging, only: t_log_scope, init_log_scope
      implicit none

      LOG_SCOPE_TRACED("subgrid_finalize", SUBGRID_TRACING)
      if (allocated(setPhi)) deallocate (setPhi)
      if (allocated(fricVertTab)) deallocate (fricVertTab)
      if (allocated(cadvVertTab)) deallocate (cadvVertTab)
      if (allocated(gridDepthVertTab)) deallocate (gridDepthVertTab)
      if (allocated(wetDepthVertTab)) deallocate (wetDepthVertTab)
      if (allocated(wetFracVertTab)) deallocate (wetFracVertTab)
      if (allocated(manningAvg)) deallocate (manningAvg)
      if (allocated(isSubgridVert)) deallocate (isSubgridVert)
      if (allocated(wetFracVertETA1)) deallocate (wetFracVertETA1)
      if (allocated(wetFracVertETA2)) deallocate (wetFracVertETA2)
      if (allocated(gridDepthVertETA1)) deallocate (gridDepthVertETA1)
      if (allocated(gridDepthVertETA2)) deallocate (gridDepthVertETA2)
      if (allocated(wetDepthVertETA2)) deallocate (wetDepthVertETA2)
      if (allocated(fricVertETA2)) deallocate (fricVertETA2)
      if (allocated(cadvVertETA1)) deallocate (cadvVertETA1)
      if (allocated(cadvVertETA2)) deallocate (cadvVertETA2)
      numPhi = 0
   end subroutine subgrid_finalize

end module mod_subgrid
