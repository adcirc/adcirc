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
#ifdef ADCNETCDF

#include "logging_macros.h"

module netcdf_error

   implicit none

   private

   public :: check_err

contains

!-----------------------------------------------------------------------
!     S U B R O U T I N E   C H E C K  _  E R R
!-----------------------------------------------------------------------
!     jgf49.17.02 Checks the return value from netCDF calls; if there
!     was an error, it writes the error message to the screen and to the
!     fort.16 file.
!-----------------------------------------------------------------------
   subroutine check_err(iret)
      use netcdf, only: NF90_NOERR, nf90_strerror
      use mod_terminate, only: terminate, ADCIRC_EXIT_FAILURE
      use mod_logging, only: allMessage, t_log_scope, init_log_scope
      implicit none
      integer, intent(in) :: iret

      LOG_SCOPE_TRACED("check_err", NETCDF_TRACING)

      if (iret /= NF90_NOERR) then
         call terminate(exit_code=ADCIRC_EXIT_FAILURE, &
                        message=nf90_strerror(iret))
      end if
!-----------------------------------------------------------------------
   end subroutine check_err
!-----------------------------------------------------------------------

end module netcdf_error
#endif
