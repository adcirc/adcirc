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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                         M O D U L E  S U B G R I D
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! EXECUTES THE SUBGRID SUBROUTINES THROUGHOUT THE CODE !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FIRST WE CREATE ALL OF THE VARIABLES NEEDED FOR THE SUBGRID ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mod_subgrid

   implicit none

   integer, parameter :: SUBGRID_DISABLED = 100000
   integer, parameter :: SUBGRID_LEVEL_0 = 100001
   integer, parameter :: SUBGRID_LEVEL_1 = 100002

   type t_subgrid
      logical :: subgrid_enabled = .false. !> true if subgrid model is active
      logical :: phi_time_derivative_enabled = .false. !> flag for subgrid phi time derivative
      integer :: subgrid_level = SUBGRID_DISABLED !> subgrid level
      integer :: numPhi = 0 !> number of possible phi values
      real(8) :: dphidt = 0d0 !> used to 0 out the phi time derivative term
      character(len=256) :: filename = "null" !> path to the subgrid lookup file

      integer, allocatable :: subgridVertList(:) !> array of vertex numbers in the subgrid area
      real(8), allocatable :: setPhi(:) !> array of phi values from 0 to 1
      real(8), allocatable :: cfVertTab(:, :) !> averaged bottom friction coefficient without level 1 correction
      real(8), allocatable :: cmfVertTab(:, :) !> averaged bottom friction coefficient with level 1 correction
      real(8), allocatable :: cadvVertTab(:, :) !> averaged advection with level 1 correction
      real(8), allocatable :: gridDepthVertTab(:, :) !> grid averaged total water depths corresponding to each phi value on the vertex
      real(8), allocatable :: wetDepthVertTab(:, :) !> wet averaged total water depths corresponding to each phi value on the vertex
      real(8), allocatable :: wetFracVertTab(:, :) !> depths corresponding to each phi value on the vertex

      real(8), allocatable :: wetFracVertETA1(:) !> vertex averaged wet fractions for zeta at k
      real(8), allocatable :: wetFracVertETA2(:) !> vertex averaged wet fractions for zeta at k+1
      real(8), allocatable :: gridDepthVertETA1(:) !> vertex grid water depths for zeta at k
      real(8), allocatable :: gridDepthVertETA2(:) !> vertex grid water depths for zeta at k+1
      real(8), allocatable :: wetDepthVertETA2(:) !> vertex wet water depths for zeta at k+1
      real(8), allocatable :: cfVertETA2(:) !> vertex bottom friction coefficients for zeta at k+1
      real(8), allocatable :: cmfVertETA2(:) !> vertex corrected bottom friction coefficients for zeta at k+1
      real(8), allocatable :: cadvVertETA2(:) !> vertex corrected advection for zeta at k+1

   contains

      procedure, private, pass(self) :: allocate => allocate_subgrid_variables
      procedure, pass(self) :: active => subgrid_active
      procedure, pass(self) :: level => subgrid_level
      procedure, pass(self) :: compute => compute_subgrid_quantities
      procedure, pass(self) :: levels_greater => subgrid_levels_greater

      ! Interpolation routines for subgrid data
      procedure, private, pass(self) :: subgrid_interpolate_2d
      procedure, private, pass(self) :: subgrid_interpolate_1d
      generic :: interpolate => subgrid_interpolate_2d, subgrid_interpolate_1d

#ifdef ADCNETCDF
      procedure, pass(self) :: read => read_subgrid_lookup_table
#endif

   end type t_subgrid

   interface t_subgrid
      module procedure initialize_subgrid_module
   end interface t_subgrid

   type(t_subgrid) :: subgrid_data !> subgrid data structure

   private

   public :: subgrid_data, t_subgrid, SUBGRID_DISABLED, SUBGRID_LEVEL_0, SUBGRID_LEVEL_1

contains

   function initialize_subgrid_module(np, subgrid_filename, subgrid_level0_active, subgrid_level1_active, &
                                      phi_time_derivative_active) result(instance)
      implicit none
      type(t_subgrid) :: instance
      integer, intent(in) :: np
      logical, intent(in) :: subgrid_level0_active
      logical, intent(in) :: subgrid_level1_active
      logical, intent(in) :: phi_time_derivative_active
      character(len=*), intent(in) :: subgrid_filename

      if (subgrid_level1_active) then
         instance%subgrid_enabled = .true.
         instance%subgrid_level = SUBGRID_LEVEL_1
      else if (subgrid_level0_active) then
         instance%subgrid_enabled = .true.
         instance%subgrid_level = SUBGRID_LEVEL_0
      else
         instance%subgrid_enabled = .false.
         instance%subgrid_level = SUBGRID_DISABLED
      end if

      if (instance%subgrid_enabled) then

         instance%phi_time_derivative_enabled = phi_time_derivative_active
         instance%filename = subgrid_filename

         if (instance%phi_time_derivative_enabled) then
            instance%dphidt = 1.d0
         else
            instance%dphidt = 0.d0
         end if

         call instance%allocate(np)

      end if

   end function initialize_subgrid_module

   pure logical function subgrid_active(self) result(is_active)
      implicit none
      class(t_subgrid), intent(in) :: self
      is_active = self%subgrid_enabled
   end function subgrid_active

   pure integer function subgrid_level(self) result(level)
      implicit none
      class(t_subgrid), intent(in) :: self
      level = self%subgrid_level
   end function subgrid_level

   pure integer function subgrid_levels_greater(self, node_index, water_level) result(numGreater)
      implicit none
      class(t_subgrid), intent(in) :: self
      integer, intent(in) :: node_index !> index of the node to check
      real(8), intent(in) :: water_level !> current water level
      integer :: first_greater_idx(1)

      if (water_level >= self%wetFracVertTab(node_index, self%numPhi)) then
         numGreater = 0
      else if (water_level < self%wetFracVertTab(node_index, 1)) then
         numGreater = self%numPhi
      else
         ! Find the first index where value > water_level
         ! minloc returns the index of the minimum value that satisfies the mask
         first_greater_idx = minloc(self%wetFracVertTab(node_index, 1:self%numPhi), &
                                    mask=self%wetFracVertTab(node_index, 1:self%numPhi) > water_level)
         numGreater = self%numPhi - first_greater_idx(1) + 1
      end if
   end function subgrid_levels_greater

   pure real(8) function subgrid_interpolate_2d(self, node_index, water_level, n_greater, x, y) result(sg)
      implicit none
      class(t_subgrid), intent(in) :: self
      integer, intent(in) :: node_index !> index of the node to interpolate
      real(8), intent(in) :: water_level !> current water level
      integer, intent(in) :: n_greater !> number of phi values greater than the current water level
      real(8), intent(in) :: x(:, :) !> array of quantities to interpolate
      real(8), intent(in) :: y(:, :) !> array of quantities to interpolate
      real(8) :: y0, y1, x0, x1

      if (n_greater == self%numPhi) then
         ! if all phi values are greater than the water level, return the first value
         sg = y(node_index, 1)
      else if (n_greater == 0) then
         ! if no phi values are greater than the water level, return the last value
         sg = y(node_index, self%numPhi)
      else
         ! Otherwise, interpolate
         y0 = y(node_index, self%numPhi - n_greater)
         y1 = y(node_index, self%numPhi - n_greater + 1)
         x0 = x(node_index, self%numPhi - n_greater)
         x1 = x(node_index, self%numPhi - n_greater + 1)
         sg = ((water_level - x0)/(x1 - x0))*(y1 - y0) + y0
      end if

   end function subgrid_interpolate_2d

   pure real(8) function subgrid_interpolate_1d(self, node_index, water_level, n_greater, x, y) result(sg)
      implicit none
      class(t_subgrid), intent(in) :: self
      integer, intent(in) :: node_index !> index of the node to interpolate
      real(8), intent(in) :: water_level !> current water level
      integer, intent(in) :: n_greater !> number of phi values greater than the current water level
      real(8), intent(in) :: x(:, :) !> array of quantities to interpolate
      real(8), intent(in) :: y(:) !> array of quantities to interpolate
      real(8) :: y0, y1, x0, x1

      if (n_greater == self%numPhi) then
         ! if all phi values are greater than the water level, return the first value
         sg = y(1)
      else if (n_greater == 0) then
         ! if no phi values are greater than the water level, return the last value
         sg = y(self%numPhi)
      else
         ! Otherwise, interpolate
         y0 = y(self%numPhi - n_greater)
         y1 = y(self%numPhi - n_greater + 1)
         x0 = x(node_index, self%numPhi - n_greater)
         x1 = x(node_index, self%numPhi - n_greater + 1)
         sg = ((water_level - x0)/(x1 - x0))*(y1 - y0) + y0
      end if

   end function subgrid_interpolate_1d

   !----------------------------------------------------------------------
   !----------------------------------------------------------------------
   !
   !                  INTIALIZE SUBGRID VARIABLES
   !
   !----------------------------------------------------------------------
   !----------------------------------------------------------------------
   ! THIS SECTION OF CODE INITIALIZES THE ARRAYS NEEDED FOR THE LOOKUP
   ! TABLES

   subroutine allocate_subgrid_variables(self, np)
      implicit none
      class(t_subgrid), intent(inout) :: self
      integer, intent(in) :: np

      allocate (self%wetFracVertETA1(NP), source=0d0)
      allocate (self%gridDepthVertETA1(NP), source=0d0)
      allocate (self%wetFracVertETA2(NP), source=0d0)
      allocate (self%gridDepthVertETA2(NP), source=0d0)
      allocate (self%wetDepthVertETA2(NP), source=0d0)

      if (self%subgrid_level == SUBGRID_LEVEL_1) then
         allocate (self%cmfVertETA2(NP), source=0d0)
         allocate (self%cadvVertETA2(NP), source=1d0)
      else
         allocate (self%cfVertETA2(NP), source=0d0)
      end if

   end subroutine allocate_subgrid_variables
   !----------------------------------------------------------------------

   !----------------------------------------------------------------------
   !
   !                  SUBROUTINE READ IN NETCDF FILE
   !
   !----------------------------------------------------------------------
   ! READS IN THE NETCDF LOOKUP TABLE AND POPULATES THE SUBGRID ARRAYS BOTH
   ! FOR SERIAL AND PARALLEL ADCIRC. THIS SUBROUTINE USES A FUNCTION
   ! "CHECK" TO READ IN THE NETCDF FILE. THE "CHECK" SUBROUTINE IS LOCATED
   ! AT THE BOTTOM OF THIS FILE.
#ifdef ADCNETCDF
   subroutine read_subgrid_lookup_table(self)
      use netcdf, only: NF90_OPEN, NF90_NOWRITE, NF90_INQ_DIMID, &
                        NF90_INQUIRE_DIMENSION, NF90_INQ_VARID, &
                        NF90_GET_VAR, NF90_CLOSE
      use global, only: allMessage, ERROR
      use netcdf_error, only: CHECK_ERR
      use mesh, only: np
#ifdef CMPI
      use global, only: nodes_lg
#endif

      implicit none
      class(t_subgrid), intent(inout) :: self

      ! Total number of elements and vertices in mesh (can most
      ! definitely be replaced by variables already in ADCIRC)
      integer :: numVerts
      ! Used for transfering from global to local
      ! integer :: currEle, currNode
      integer :: currNode
      ! Arrays for holding global variables
      real(8), dimension(:, :), allocatable :: globalVertLookup
      integer, dimension(:), allocatable :: globalSubVertList
      ! variables needed to read in main NETCDF lookup table
      integer :: NC_ERR, NC_ID, NC_VAR
      integer :: I
      logical :: file_exists

      inquire (FILE=trim(self%filename), EXIST=file_exists)
      if (.not. file_exists) then
         call allMessage(ERROR, "subgrid lookup file does not exist")
         call terminate()
         call exit(1)
      end if

!JLW: adding subgrid lookup table read in
!JLW: first open subgrid lookup table and read in dimensions
      call CHECK_ERR(NF90_OPEN(trim(self%filename), NF90_NOWRITE, NC_ID))

      call CHECK_ERR(NF90_INQ_DIMID(NC_ID, "numNode", NC_VAR))
      call CHECK_ERR(NF90_INQUIRE_DIMENSION(NC_ID, NC_VAR, len=numVerts))

      call CHECK_ERR(NF90_INQ_DIMID(NC_ID, "numPhi", NC_VAR))
      call CHECK_ERR(NF90_INQUIRE_DIMENSION(NC_ID, NC_VAR, len=self%numPhi))

      allocate (self%setPhi(self%numPhi))

      call check_err(NF90_INQ_VARID(NC_ID, "phiSet", NC_VAR))
      call check_err(NF90_GET_VAR(NC_ID, NC_VAR, self%setPhi))

!JLW: go through and read in each of the vertex averaged variables
!JLW: allocate arrays
      if (self%level() == SUBGRID_LEVEL_1) then
         allocate (self%cmfVertTab(NP, self%numPhi))
         allocate (self%cadvVertTab(NP, self%numPhi))
      else
         allocate (self%cfVertTab(NP, self%numPhi))
      end if
      allocate (self%gridDepthVertTab(NP, self%numPhi))
      allocate (self%wetDepthVertTab(NP, self%numPhi))
      allocate (self%wetFracVertTab(NP, self%numPhi))
      allocate (self%subgridVertList(NP))

!!!!!!!!!!!!!!!!!! bottom friction coefficient !!!!!!!!!!!!!!!!!!!!!!
      allocate (globalVertLookup(self%numPhi, numVerts))

      if (self%level() == SUBGRID_LEVEL_1) then
         call check_err(NF90_INQ_VARID(NC_ID, "cmfVertex", NC_VAR))
         call check_err(NF90_GET_VAR(NC_ID, NC_VAR, globalVertLookup))
         do I = 1, NP
#ifdef CMPI
            currNode = abs(nodes_lg(I))
#else
            currNode = I
#endif
            self%cmfVertTab(I, :) = globalVertLookup(:, currNode)
         end do
      else
         call check_err(NF90_INQ_VARID(NC_ID, "cfVertex", NC_VAR))
         call check_err(NF90_GET_VAR(NC_ID, NC_VAR, globalVertLookup))
         do I = 1, NP
#ifdef CMPI
            currNode = abs(nodes_lg(I))
#else
            currNode = I
#endif
            self%cfVertTab(I, :) = globalVertLookup(:, currNode)
         end do
      end if
!!!!!!!!!!!!!!!!!! ADVECTION CORRECTION VETEX !!!!!!!!!!!!!!!!!!
      if (self%level() == SUBGRID_LEVEL_1) then
         call check_err(NF90_INQ_VARID(NC_ID, "cadvVertex", NC_VAR))
         call check_err(NF90_GET_VAR(NC_ID, NC_VAR, globalVertLookup))
         do I = 1, NP
#ifdef CMPI
            currNode = abs(nodes_lg(I))
#else
            currNode = I
#endif
            self%cadvVertTab(I, :) = globalVertLookup(:, currNode)
         end do
      end if
!!!!!!!!!!!!!!!!! grid depth vertex !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call check_err(NF90_INQ_VARID(NC_ID, "gridTotWatDepthVertex", NC_VAR))
      call check_err(NF90_GET_VAR(NC_ID, NC_VAR, globalVertLookup))
      do I = 1, NP
#ifdef CMPI
         currNode = abs(nodes_lg(I))
#else
         currNode = I
#endif
         self%gridDepthVertTab(I, :) = globalVertLookup(:, currNode)
      end do
!!!!!!!!!!!!!!!!! wet depth vertex !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call check_err(NF90_INQ_VARID(NC_ID, "wetTotWatDepthVertex", NC_VAR))
      call check_err(NF90_GET_VAR(NC_ID, NC_VAR, globalVertLookup))
      do I = 1, NP
#ifdef CMPI
         currNode = abs(nodes_lg(I))
#else
         currNode = I
#endif
         self%wetDepthVertTab(I, :) = globalVertLookup(:, currNode)
      end do
!!!!!!!!!!!!!!!!!!!! wet fraction depths vertex !!!!!!!!!!!!!!!!!!!!!!
      call check_err(NF90_INQ_VARID(NC_ID, "wetFractionVertex", NC_VAR))
      call check_err(NF90_GET_VAR(NC_ID, NC_VAR, globalVertLookup))
      do I = 1, NP
#ifdef CMPI
         currNode = abs(nodes_lg(I))
#else
         currNode = I
#endif
         self%wetFracVertTab(I, :) = globalVertLookup(:, currNode)
      end do
!!!!!!!!!!!!!!!!! subgrid flag for vertex !!!!!!!!!!!!!!!!!!!!!!!!
      allocate (globalSubVertList(numVerts))

      call check_err(NF90_INQ_VARID(NC_ID, "binaryVertexList", NC_VAR))
      call check_err(NF90_GET_VAR(NC_ID, NC_VAR, globalSubVertList))

      do I = 1, NP
#ifdef CMPI
         currNode = abs(nodes_lg(I))
#else
         currNode = I
#endif
         self%subgridVertList(I) = globalSubVertList(currNode)
      end do

!JLW: close netcdf lookup table file
      NC_ERR = NF90_CLOSE(NC_ID)

   end subroutine read_subgrid_lookup_table
#endif
   !----------------------------------------------------------------------

   !---------------------------------------------------------------------
   !
   !                  GET VERTEX AVERAGED VARIABLES
   !
   !---------------------------------------------------------------------
   !---------------------------------------------------------------------
   ! THIS CODE IS VERY SIMILAR TO THE ELEMENTAL SUBROUTINE BUT FOR SUBRID
   ! VARIABLES RELATED TO THE VERTEX AREAS. FOR EACH NEW WATER SURFACE
   ! ELEVATION SUBGRID VERTEX QUANTITIES ARE LOOKED UP.

   subroutine compute_subgrid_quantities(self, np, h0, dp, eta2)
      implicit none
      class(t_subgrid), intent(inout) :: self
      integer, intent(in) :: np !> number of vertices
      real(8), intent(in) :: h0 !> minimum water depth
      real(8), intent(in) :: dp(np) !> nodal depth
      real(8), intent(in) :: eta2(np) !> nodal water surface elevation

      integer :: i
      integer :: numGreater
      real(8) :: HABSMIN

      !JLW: set k = k + 1 for new timestep
      HABSMIN = 0.8d0*H0
      self%wetFracVertETA1 = self%wetFracVertETA2
      self%gridDepthVertETA1 = self%gridDepthVertETA2
      !JLW: now loop through the vertices and lookup vertex averages
      do I = 1, NP
         if (self%subgridVertList(I) == 1) then
            numGreater = self%levels_greater(I, eta2(I))
            self%wetFracVertETA2(I) = self%interpolate(I, eta2(I), numGreater, self%wetFracVertTab, self%setPhi)
            self%gridDepthVertETA2(I) = self%interpolate(I, eta2(I), numGreater, self%wetFracVertTab, self%gridDepthVertTab)
            self%wetDepthVertETA2(I) = self%interpolate(I, eta2(I), numGreater, self%wetFracVertTab, self%wetDepthVertTab)
            if (self%level() == SUBGRID_LEVEL_1) then
               self%cmfVertETA2(I) = self%interpolate(I, eta2(I), numGreater, self%wetFracVertTab, self%cmfVertTab)
               self%cadvVertETA2(I) = self%interpolate(I, eta2(I), numGreater, self%wetFracVertTab, self%cadvVertTab)
            else
               self%cfVertETA2(I) = self%interpolate(I, eta2(I), numGreater, self%wetFracVertTab, self%cfVertTab)
            end if
         else
            !JLW: if not in subgrid area set grid total depth (this was a bug)
            self%gridDepthVertETA2(i) = eta2(i) + dp(i)
            self%wetDepthVertETA2(i) = eta2(i) + dp(i)
            self%wetFracVertETA2(i) = 1.d0
         end if
      end do

   end subroutine compute_subgrid_quantities
   !----------------------------------------------------------------------
#ifdef ADCNETCDF
   subroutine terminate()
#ifdef CMPI
      use MESSENGER, only: MSG_FINI, subdomainFatalError
#endif
      use GLOBAL, only: INFO, allMessage, &
                        setMessageSource, unsetMessageSource
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      use GLOBAL, only: DEBUG
#endif
      implicit none
      call setMessageSource("terminate")
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      call allMessage(INFO, "ADCIRC Terminating.")

#ifdef CMPI
      subdomainFatalError = .true.
      call MSG_FINI()
#endif
      call exit(1)
      !
#if defined(MESH_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.") ! should be unreachable
#endif
      call unsetMessageSource()
   end subroutine terminate
#endif
end module mod_subgrid
