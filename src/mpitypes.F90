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
!> @desription The mpi types module is used to define the custom mpi types and
!> operations that ADCIRC uses for specific parallel operations.
!>
!> @author Zach Cobell
!> @date 2025-04-28
!>
!-------------------------------------------------------------------------------!
module mod_adcirc_mpi_types

   implicit none

   type s_TimestepStatusData
      real(8) :: max_elevation
      real(8) :: max_velocity
      integer :: max_elevation_node
      integer :: max_elevation_rank
      integer :: max_velocity_node
      integer :: max_velocity_rank
   end type s_TimestepStatusData

   integer :: MPI_TimestepStatusData !> Identifier for the MPI type for timestep status data
   integer :: MPI_Op_TimestepDataMax !> Identifier for the reduction performed on timestep status data

   private

   public :: s_TimestepStatusData, create_mpi_types, MPI_TimestepStatusData, MPI_Op_TimestepDataMax

contains

   !--------------------------------------------------------------------
   !> Creates the custom mpi types and operations used in the code
   !> @param[in] myproc calling processor id
   !--------------------------------------------------------------------
   subroutine create_mpi_types(myproc)
      implicit none
      integer, intent(in) :: myproc
      call create_timestep_status_mpi_type(myproc)
   end subroutine create_mpi_types

   !--------------------------------------------------------------------
   !> Checks the MPI return status and writes the specified error if
   !> the status is not MPI_SUCCESS
   !> @param[in] myproc calling processor id
   !> @param[in] message message that is written to screen when there is an error
   !> @param[in] error code returned from the mpi call being checked
   !--------------------------------------------------------------------
   subroutine mpi_check(myproc, message, ierr)
      use mpi, only: MPI_COMM_WORLD, MPI_SUCCESS, MPI_ABORT
      implicit none
      integer, intent(in) :: myproc
      integer, intent(in) :: ierr
      character(len=*), intent(in) :: message
      integer :: mpi_err
      if (ierr /= MPI_SUCCESS) then
         write (*, '(A,I0)') trim(message), ierr
         call MPI_ABORT(MPI_COMM_WORLD, myproc, mpi_err)
      end if
   end subroutine mpi_check

   !--------------------------------------------------------------------
   !> Creates the timestep status mpi objects. This object packs the data
   !> so that we can perform a single reduction among all processors to find
   !> the global maximum values for velocity and water surface elevation and
   !> which processors contained those values
   !--------------------------------------------------------------------
   subroutine create_timestep_status_mpi_type(myproc)
      use mpi, only: MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_TYPE_CREATE_STRUCT, &
                     MPI_TYPE_GET_EXTENT, MPI_TYPE_CREATE_RESIZED, MPI_TYPE_COMMIT, &
                     MPI_OP_CREATE
      implicit none
      integer, intent(in) :: myproc
      integer :: IERR
      integer, parameter :: TS_COUNT = 6
      integer, parameter :: TS_BLOCKLENS(TS_COUNT) = [1, 1, 1, 1, 1, 1]
      integer(8), parameter :: TS_OFFSETS(TS_COUNT) = [0, 8, 16, 20, 24, 28]
      integer, parameter :: TS_TYPES(TS_COUNT) = &
                            [MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, &
                             MPI_INTEGER, MPI_INTEGER, MPI_INTEGER, MPI_INTEGER]
      integer(8) :: LowerBound, Extent
      integer :: TempType

      ! Create a temporary structure and then allow MPI to define a version
      ! that can be used as an array type. We probably won't be broadcasting
      ! this around as an array, however, it would be a nasty bug if someone
      ! tried one day and they had the version that could only be used as a solo
      ! type.
      call MPI_TYPE_CREATE_STRUCT(TS_COUNT, TS_BLOCKLENS, TS_OFFSETS, TS_TYPES, TempType, IERR)
      call MPI_CHECK(myproc, "ERROR: MPI_TYPE_CREATE_STRUCT", IERR)
      call MPI_TYPE_GET_EXTENT(TempType, LowerBound, Extent, IERR)
      call MPI_CHECK(myproc, "ERROR: MPI_TYPE_GET_EXTENT", IERR)
      call MPI_TYPE_CREATE_RESIZED(TempType, LowerBound, Extent, MPI_TimestepStatusData, IERR)
      call MPI_CHECK(myproc, "ERROR: MPI_TYPE_CREATE_RESIZED", IERR)
      call MPI_TYPE_COMMIT(MPI_TimestepStatusData, IERR)
      call MPI_CHECK(myproc, "ERROR: MPI_TYPE_COMMIT", IERR)
      call MPI_OP_CREATE(MPITimestepDataReduce, .true., MPI_OP_TimestepDataMax, IERR)
      call MPI_CHECK(myproc, "ERROR: MPI_OP_CREATE", IERR)

   end subroutine create_timestep_status_mpi_type

   !--------------------------------------------------------------------
   !> Creates the function used to reduce the s_TimstepStatusData type
   !> This is created as a custom MPI operation
   !>
   !> @param[in] invec Input vector of s_TimestepStatusData objects
   !> @param[inout] inoutvec Output vector of s_TimestepSTatusData objects
   !> @param[in] len Length of the input vectors
   !> @param[in] datatype Type of data used in this function
   !> @param[out] ierr Error return value
   !--------------------------------------------------------------------
   subroutine MPITimestepDataReduce(invec, inoutvec, len, datatype)
      use mpi, only: MPI_SUCCESS
      implicit none
      type(s_TimestepStatusData), intent(in) :: invec(len)
      type(s_TimestepStatusData), intent(inout) :: inoutvec(len)
      integer, intent(IN) :: len
      integer, intent(IN) :: datatype
      integer :: i

      if (datatype /= MPI_TimestepStatusData) then
         return
      end if

      do i = 1, len
         if (invec(i)%max_elevation > inoutvec(i)%max_elevation .and. &
             invec(i)%max_elevation_node > 0) then
            inoutvec(i)%max_elevation = invec(i)%max_elevation
            inoutvec(i)%max_elevation_node = invec(i)%max_elevation_node
            inoutvec(i)%max_elevation_rank = invec(i)%max_elevation_rank
         end if
         if (invec(i)%max_velocity > inoutvec(i)%max_velocity .and. &
             invec(i)%max_velocity_node > 0) then
            inoutvec(i)%max_velocity = invec(i)%max_velocity
            inoutvec(i)%max_velocity_node = invec(i)%max_velocity_node
            inoutvec(i)%max_velocity_rank = invec(i)%max_velocity_rank
         end if
      end do

   end subroutine MPITimestepDataReduce

end module mod_adcirc_mpi_types
