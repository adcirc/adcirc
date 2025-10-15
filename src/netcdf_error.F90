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
module netcdf_error

   implicit none

   private

   public :: check_err, netcdfTerminate

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
      use global, only: ERROR, allMessage, setMessageSource, unsetMessageSource
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      use global, only: DEBUG
#endif
#ifdef CMPI
      use MESSENGER, only: MSG_FINI
#endif
      implicit none
      integer, intent(in) :: iret

      call setMessageSource("check_err")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif
      if (iret /= NF90_NOERR) then
         call allMessage(ERROR, nf90_strerror(iret))
         call netcdfTerminate()
      end if
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
!-----------------------------------------------------------------------
   end subroutine check_err
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     S U B R O U T I N E   N E T C D F   T E R M I N A T E
!-----------------------------------------------------------------------
   subroutine netcdfTerminate(NO_MPI_FINALIZE)
#ifdef CMPI
      use MESSENGER, only: msg_fini, subdomainFatalError
#endif
      use GLOBAL, only: setMessageSource, unsetMessageSource, &
                        allMessage, INFO, allMessage
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      use global, only: DEBUG
#endif
      implicit none
      logical, intent(in), optional :: NO_MPI_FINALIZE
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      real(8), allocatable :: dummy(:)
#endif

      call setMessageSource("netcdfTerminate")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Enter.")
#endif

      call allMessage(INFO, "ADCIRC Terminating.")

#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      ! intentionally create a segmentation fault so that we can get
      ! a stack trace to determine the line number of the netcdf call
      ! that went bad ... this assumes that the code was compiled with
      ! debugging symbols, bounds checking, and stack trace turned on.
      allocate (dummy(1)) ! Allocating (too small) so that -Wuninitialized doesn't complain
      dummy(2) = 99.9d0
#endif

      if (present(NO_MPI_FINALIZE)) then
#ifdef CMPI
         subdomainFatalError = .true.
         call MSG_FINI(NO_MPI_FINALIZE)
#endif
         call exit(1)
      else
#ifdef CMPI
         subdomainFatalError = .true.
         call MSG_FINI()
#endif
         call exit(1)
      end if

#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.") ! should be unreachable
#endif
      call unsetMessageSource()
!-----------------------------------------------------------------------
   end subroutine netcdfTerminate
!-----------------------------------------------------------------------

end module netcdf_error
#endif
