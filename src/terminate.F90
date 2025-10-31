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
!-----------------------------------------------------------------------
!  MODULE TERMINATE
!-----------------------------------------------------------------------
!> @brief Module for graceful termination of ADCIRC with proper cleanup
!>
!> This module provides a standardized way to terminate ADCIRC execution
!> with appropriate error codes, message logging, and MPI finalization. It ensures
!> that all processes shut down cleanly and error messages are properly logged.
!-----------------------------------------------------------------------
module mod_terminate
   implicit none

   integer, parameter :: ADCIRC_EXIT_SUCCESS = 0 !< successful exit code
   integer, parameter :: ADCIRC_EXIT_FAILURE = 1 !< failure exit code

   private

   public :: terminate, ADCIRC_EXIT_SUCCESS, ADCIRC_EXIT_FAILURE

contains

!-----------------------------------------------------------------------
!> @brief Terminate ADCIRC execution with proper cleanup and error handling
!>
!> This subroutine gracefully terminates ADCIRC by logging any error
!> messages, setting error flags for parallel execution, finalizing MPI if requested,
!> and calling exit with the specified exit code. Error messages are logged to both
!> screen and log file, while success messages are logged at INFO level.
!>
!> @param[in] exit_code       (optional) exit code to return, defaults to ADCIRC_EXIT_FAILURE
!> @param[in] do_finalize_mpi (optional) whether to call MPI_Finalize, defaults to .true.
!> @param[in] message         (optional) message to log before terminating
!-----------------------------------------------------------------------
   subroutine terminate(exit_code, do_finalize_mpi, message)
#ifdef CMPI
      use MESSENGER, only: subdomainFatalError, msg_fini
#endif
      use global, only: allMessage, INFO, ERROR
      implicit none

      integer, optional, intent(in) :: exit_code
      logical, optional, intent(in) :: do_finalize_mpi
      character(*), optional, intent(in) :: message
      integer :: exit_code_local
      logical :: finalize_mpi
      character(len=256) :: scratchMessage

      if (present(exit_code)) then
         exit_code_local = exit_code
      else
         exit_code_local = ADCIRC_EXIT_FAILURE
      end if

      if (present(do_finalize_mpi)) then
         finalize_mpi = do_finalize_mpi
      else
         finalize_mpi = .true.
      end if

      if (exit_code_local /= ADCIRC_EXIT_SUCCESS) then
         if (present(message)) then
            call allMessage(ERROR, message)
         end if
         write (scratchMessage, "(A,I0)") "ADCIRC terminating with error code: ", exit_code_local
         call allMessage(ERROR, scratchMessage)
      else
         if (present(message)) then
            call allMessage(INFO, message)
         end if
      end if
#ifdef CMPI
      if (exit_code_local /= ADCIRC_EXIT_SUCCESS) subdomainFatalError = .true.
      call msg_fini(finalize_mpi)
#endif
      call exit(exit_code_local)

   end subroutine terminate

end module mod_terminate
