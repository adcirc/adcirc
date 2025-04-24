#ifdef ADCNETCDF
module netcdf_error

   use NETCDF, only: NF90_NOERR, nf90_strerror

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
      use sizes, only: myproc
      use mod_terminate, only: terminate
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
         call terminate(myproc)
      end if
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG, "Return.")
#endif
      call unsetMessageSource()
!-----------------------------------------------------------------------
   end subroutine check_err
!-----------------------------------------------------------------------

      end module netcdf_error
#endif
