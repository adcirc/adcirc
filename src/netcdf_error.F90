#ifdef ADCNETCDF
module netcdf_error
      
    USE NETCDF

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
      USE SIZES, ONLY : myproc
      USE mod_logging, ONLY: screenUnit, ERROR, allMessage, &
          setMessageSource, unsetMessageSource
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      USE mod_logging, ONLY : DEBUG
#endif
#ifdef CMPI
      USE MESSENGER, ONLY : MSG_FINI
#endif
      IMPLICIT NONE
      INTEGER, intent(in) :: iret

      call setMessageSource("check_err")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif
      if (iret .ne. NF90_NOERR) then
         call allMessage(ERROR,nf90_strerror(iret))
         call netcdfTerminate()
      endif
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.")
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
      USE MESSENGER
#endif
      USE mod_logging, ONLY : setMessageSource, unsetMessageSource,&
         allMessage, DEBUG, ECHO, INFO, WARNING, ERROR, allMessage
      IMPLICIT NONE
      LOGICAL, OPTIONAL :: NO_MPI_FINALIZE
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      REAL(8), ALLOCATABLE :: dummy(:)
#endif

      call setMessageSource("netcdfTerminate")
#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Enter.")
#endif

      call allMessage(INFO,"ADCIRC Terminating.")

#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
         ! intentionally create a segmentation fault so that we can get
         ! a stack trace to determine the line number of the netcdf call
         ! that went bad ... this assumes that the code was compiled with
         ! debugging symbols, bounds checking, and stack trace turned on.
         allocate(dummy(1)) ! Allocating (too small) so that -Wuninitialized doesn't complain
         dummy(2) = 99.9d0
#endif

#ifdef CMPI
      subdomainFatalError = .true.
      IF (PRESENT(NO_MPI_FINALIZE)) THEN
        CALL MSG_FINI(NO_MPI_FINALIZE)
      ELSE
        CALL MSG_FINI()
      ENDIF
#endif
      CALL EXIT(1) 

#if defined(NETCDF_TRACE) || defined(ALL_TRACE)
      call allMessage(DEBUG,"Return.") ! should be unreachable
#endif
      call unsetMessageSource()
!-----------------------------------------------------------------------
      end subroutine netcdfTerminate
!-----------------------------------------------------------------------

end module netcdf_error
#endif
