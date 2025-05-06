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
module mod_terminate

   implicit none

   private

   public :: terminate

contains

   !-----------------------------------------------------------------------
   ! ADCIRC terminate routine
   ! Needed so that we can clean up mpi when bombing out.
   !-----------------------------------------------------------------------
   subroutine terminate(myproc, error_code_in)
#ifdef CMPI
      use mpi, only: mpi_comm_world, mpi_abort
#endif
      implicit none
      integer, intent(in) :: myproc
      integer, optional, intent(in) :: error_code_in
      integer :: error_code
#ifdef CMPI
      integer :: ierr
#endif

      if (present(error_code_in)) then
         error_code = error_code_in
      else
         error_code = 1
      end if

      if (myproc == 0) then
         write (*, '(A)') "INFO: ADCIRC Terminating"
      end if
      write (16, '(A)') "INFO: ADCIRC Terminating"
      call flush (16)

#ifdef CMPI
      call mpi_abort(MPI_COMM_WORLD, error_code, ierr)
#endif

      call exit(error_code)

   end subroutine terminate
   !-----------------------------------------------------------------------

#ifdef CMPI
   subroutine MSG_ABORT()
      use mpi, only: mpi_comm_world, mpi_abort

      implicit none

      integer, parameter :: myerrcode = 1
      integer :: myerr

      call mpi_abort(mpi_comm_world, myerrcode, myerr)
   end subroutine MSG_ABORT
#endif

end module mod_terminate
