#ifdef ADCNETCDF
module netcdf_error

   use NETCDF

   implicit none

contains

!-----------------------------------------------------------------------
!     S U B R O U T I N E   C H E C K  _  E R R
!-----------------------------------------------------------------------
!     jgf49.17.02 Checks the return value from netCDF calls; if there
!     was an error, it writes the error message to the screen and to the
!     fort.16 file.
!-----------------------------------------------------------------------
   subroutine check_err(iret)
      use mod_logging, only: ERROR, allMessage, &
                             setMessageSource, unsetMessageSource
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      use mod_logging, only: DEBUG
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
#ifdef CMPI
   subroutine netcdfTerminate(NO_MPI_FINALIZE)
      use MESSENGER
#else
      subroutine netcdfTerminate()
#endif
         use mod_logging, only: INFO, allMessage, setMessageSource, unsetMessageSource
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
         use mod_logging, only: DEBUG, ECHO, INFO, WARNING, ERROR
#endif
         implicit none
#ifdef CMPI
         logical, optional :: NO_MPI_FINALIZE
#endif
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

#ifdef CMPI
         subdomainFatalError = .true.
         if (present(NO_MPI_FINALIZE)) then
            call MSG_FINI(NO_MPI_FINALIZE)
         else
            call MSG_FINI()
         end if
#endif
         call exit(1)

#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
         call allMessage(DEBUG, "Return.") ! should be unreachable
#endif
         call unsetMessageSource()
!-----------------------------------------------------------------------
      end subroutine netcdfTerminate
!-----------------------------------------------------------------------

      end module netcdf_error
#endif
