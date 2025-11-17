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
!  MODULE IO
!-----------------------------------------------------------------------
!>
!> @brief I/O utility module for ADCIRC
!>
!> This module contains I/O utility subroutines for ADCIRC including
!> file opening with error checking, file existence verification, and standardized
!> error handling for file operations.
!-----------------------------------------------------------------------

#include "logging_macros.h"

module mod_io
   implicit none

   private

   public :: openFileForRead

contains

!-----------------------------------------------------------------------
!> @brief Open an existing file for reading with error checking
!>
!> This subroutine opens an existing input file for reading and includes
!> comprehensive error checking. It first verifies the file exists, then attempts to
!> open it. If the file is marked as required (default), the program will terminate
!> on failure. If marked as optional, the subroutine returns with errorIO set to non-zero.
!> All operations are logged appropriately.
!>
!> @param[in]  lun      Fortran logical unit number to use for the file
!> @param[in]  filename full pathname of the file to open
!> @param[out] errorIO  zero if the file opened successfully, non-zero otherwise
!> @param[in]  required (optional) .true. = terminate if file not found (default), .false. = return with error
!-----------------------------------------------------------------------
   subroutine openFileForRead(lun, filename, errorIO, required)
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      use mod_logging, only: allMessage, logMessage, INFO, ERROR, t_log_scope, init_log_scope
      integer, intent(in) :: lun ! fortran logical unit number
      character(*), intent(in) :: filename ! full pathname of file
      integer, intent(out) :: errorIO ! zero if the file opened successfully
      logical, intent(in), optional :: required ! .true. = terminate if file not found (default)
      logical :: fileFound ! .true. if the file is present
      logical :: fileRequired ! local copy of required parameter
      character(len=256) :: scratchMessage ! message buffer

      LOG_SCOPE_TRACED("openFileForRead", ALL_TRACING)

      ! Set default for required parameter
      if (present(required)) then
         fileRequired = required
      else
         fileRequired = .true. ! default: file is required
      end if

      errorIO = 0

      ! Check to see if file exists
      write (scratchMessage, 21) lun
21    format("Searching for file to open on unit ", i0, "...")
      call logMessage(INFO, trim(scratchMessage))
      inquire (file=trim(filename), exist=fileFound)
      if (fileFound .eqv. .false.) then
         write (scratchMessage, 23) trim(filename)
23       format("The file '", A, "' was not found.")
         call allMessage(INFO, scratchMessage)
         errorIO = 1
         if (fileRequired) then
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message='Required file not found: '//trim(filename))
         end if
         return ! file not found
      else
         write (scratchMessage, 24) trim(filename)
24       format("The file '", A, "' was found. The file will be opened.")
         call logMessage(INFO, trim(scratchMessage))
      end if

      ! Open existing file
      open (lun, file=trim(filename), status='OLD', &
            action='READ', iostat=errorIO)
      if (errorIO /= 0) then
         write (scratchMessage, 25) trim(filename)
25       format("Could not open the file '", A, "'.")
         call allMessage(ERROR, trim(scratchMessage))
         if (fileRequired) then
            call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                           message='Could not open required file: '//trim(filename))
         end if
         return ! file found but could not be opened
      else
         write (scratchMessage, 26) trim(filename)
26       format("The file '", A, "' was opened successfully.")
         call logMessage(INFO, trim(scratchMessage))
      end if

!-----------------------------------------------------------------------
   end subroutine openFileForRead
!-----------------------------------------------------------------------

end module mod_io
