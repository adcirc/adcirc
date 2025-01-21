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
  use sizes, only: myproc

contains

   !-----------------------------------------------------------------------
   ! ADCIRC terminate routine
   ! Needed so that we can clean up mpi when bombing out.
   !-----------------------------------------------------------------------
   subroutine terminate()
      implicit none
#ifdef CMPI
      call msg_abort()
#endif
      if (myproc == 0) then
         write (*, '(A)') "INFO: ADCIRC Terminating"
      end if
      write (16, '(A)') "INFO: ADCIRC Terminating"
      call flush(16)
      call exit(1)

   end subroutine terminate
   !-----------------------------------------------------------------------

#ifdef CMPI
   subroutine MSG_ABORT()
      use mpi, only: mpi_comm_world, mpi_abort
      implicit none
      integer myerr, myerrcode
      call mpi_abort(mpi_comm_world, myerrcode, myerr)
   end subroutine MSG_ABORT
#endif

end module mod_terminate
