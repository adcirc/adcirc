!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
!  MODULE LOGGING
!-----------------------------------------------------------------------
!>
!> @brief Logging module for ADCIRC message handling
!>
!> This module provides a hierarchical logging system for ADCIRC
!> with configurable logging levels and multiple output destinations. It supports
!> writing to screen, log files (fort.16), and includes call stack
!> tracking when compiled with FULL_STACK.
!>
!> Logging levels (from most to least verbose):
!>   - DEBUG (-1): Write all messages and echo input
!>   - ECHO (0): Echo input plus write all non-debug messages
!>   - INFO (1): Don't echo input, write all non-debug messages
!>   - WARNING (2): Don't echo input, write only warnings and errors
!>   - ERROR (3): Don't echo input, only fatal messages
!-----------------------------------------------------------------------
module mod_logging
   use iso_fortran_env, only: output_unit, error_unit
   implicit none

   private

   type t_log_level
      private
      integer :: log_level !< logging level
      integer :: io_unit !< associated unit number
      character(len=10) :: log_level_name !< name of logging level

   contains
      procedure, public, pass(self) :: level => t_log_level_level
      procedure, public, pass(self) :: unit => t_log_level_unit
      procedure, public, pass(self) :: name => t_log_level_name
   end type t_log_level

   integer, parameter :: LOG_LEVEL_DEBUG = -1
   integer, parameter :: LOG_LEVEL_ECHO = 0
   integer, parameter :: LOG_LEVEL_INFO = 1
   integer, parameter :: LOG_LEVEL_WARNING = 2
   integer, parameter :: LOG_LEVEL_ERROR = 3

   type(t_log_level), parameter :: DEBUG = t_log_level(LOG_LEVEL_DEBUG, output_unit, "DEBUG")
   type(t_log_level), parameter :: ECHO = t_log_level(LOG_LEVEL_ECHO, output_unit, "ECHO")
   type(t_log_level), parameter :: INFO = t_log_level(LOG_LEVEL_INFO, output_unit, "INFO")
   type(t_log_level), parameter :: WARNING = t_log_level(LOG_LEVEL_WARNING, output_unit, "WARNING")
   type(t_log_level), parameter :: ERROR = t_log_level(LOG_LEVEL_ERROR, error_unit, "ERROR")

   integer, parameter :: log_unit = 16 !< log file unit number

   character(len=50), dimension(50) :: messageSources !< subroutine names
   integer :: sourceNumber = 0 !< index into messageSources for current sub

   type(t_log_level) :: nabout = t_log_level(0, output_unit, "INFO") !< logging level from fort.15
   integer :: nscreen = 1 !< screen output destination from fort.15
   integer :: screenUnit !< unit number for screen output

   interface operator(<)
      module procedure t_log_level_lt
   end interface

   interface operator(>)
      module procedure t_log_level_gt
   end interface

   interface operator(<=)
      module procedure t_log_level_lte
   end interface

   interface operator(>=)
      module procedure t_log_level_gte
   end interface

   interface operator(==)
      module procedure t_log_level_eq
   end interface

   interface operator(/=)
      module procedure t_log_level_neq
   end interface

   interface t_log_level
      module procedure t_log_level_constructor
   end interface t_log_level

   public :: DEBUG, ECHO, INFO, WARNING, ERROR, &
             openLogFile, screenMessage, logMessage, allMessage, &
             setMessageSource, unsetMessageSource, screenUnit, &
             nabout, nscreen, log_unit, operator(<), operator(>), &
             operator(<=), operator(>=), operator(==), operator(/=), &
             t_log_level

contains

!-----------------------------------------------------------------------
!> @brief Construct a log level from an integer code
!>
!> Creates a t_log_level object from an integer level code. Valid codes
!> are -1 (DEBUG), 0 (ECHO), 1 (INFO), 2 (WARNING), 3 (ERROR).
!> Invalid codes default to INFO.
!>
!> @param[in] level_code integer code for the logging level
!> @return              log level object corresponding to the code
!-----------------------------------------------------------------------
   type(t_log_level) function t_log_level_constructor(level_code) result(log_level)
      implicit none
      integer, intent(in), optional :: level_code

      select case (level_code)
      case (LOG_LEVEL_DEBUG)
         log_level = DEBUG
      case (LOG_LEVEL_ECHO)
         log_level = ECHO
      case (LOG_LEVEL_INFO)
         log_level = INFO
      case (LOG_LEVEL_WARNING)
         log_level = WARNING
      case (LOG_LEVEL_ERROR)
         log_level = ERROR
      case default
         call allMessage(WARNING, "Invalid log level code. Defaulting to 'INFO'")
         log_level = INFO
      end select

   end function t_log_level_constructor

!-----------------------------------------------------------------------
!> @brief Format a log message with level and source information
!>
!> Pure function that formats a log message according to the configured
!> format (FULL_STACK or standard). Returns the formatted string without
!> performing any I/O operations.
!>
!> @param[in] log_level   logging level
!> @param[in] message     text message to format
!> @return                formatted message string
!-----------------------------------------------------------------------
   pure function formatLogMessage(log_level, message) result(formatted)
      implicit none
      type(t_log_level), intent(in) :: log_level
      character(*), intent(in) :: message
      character(len=:), allocatable :: formatted
      character(len=1024) :: temp
#ifdef FULL_STACK
      integer :: j

      write (temp, '(A, ": ", A, 50(:, "->", A), ": ", A)') &
         trim(log_level%name()), (trim(messageSources(j)), j=1, sourceNumber), trim(message)
      formatted = trim(temp)
#else
      write (temp, '(A, ": ", A, ": ", A)') &
         trim(log_level%name()), trim(messageSources(sourceNumber)), trim(message)
      formatted = trim(temp)
#endif
   end function formatLogMessage

!-----------------------------------------------------------------------
!> @brief Get the integer level value from a log level object
!>
!> @param[in] self  log level object
!> @return          integer level value
!-----------------------------------------------------------------------
   pure integer function t_log_level_level(self) result(level)
      implicit none
      class(t_log_level), intent(in) :: self

      level = self%log_level
   end function t_log_level_level

!-----------------------------------------------------------------------
!> @brief Get the unit number associated with a log level
!>
!> @param[in] self  log level object
!> @return          unit number for output
!-----------------------------------------------------------------------
   pure integer function t_log_level_unit(self) result(unit)
      implicit none
      class(t_log_level), intent(in) :: self

      unit = self%io_unit
   end function t_log_level_unit

!-----------------------------------------------------------------------
!> @brief Get the name string of a log level
!>
!> @param[in] self  log level object
!> @return          name of the log level (e.g., "DEBUG", "INFO")
!-----------------------------------------------------------------------
   pure function t_log_level_name(self) result(name)
      implicit none
      class(t_log_level), intent(in) :: self
      character(len=10) :: name

      name = self%log_level_name
   end function t_log_level_name

!-----------------------------------------------------------------------
!> @brief Less than comparison operator for log levels
!>
!> @param[in] a first log level
!> @param[in] b second log level
!> @return    .true. if a < b
!-----------------------------------------------------------------------
   pure logical function t_log_level_lt(a, b) result(res)
      type(t_log_level), intent(in) :: a, b

      res = a%level() < b%level()
   end function t_log_level_lt

!-----------------------------------------------------------------------
!> @brief Greater than comparison operator for log levels
!>
!> @param[in] a first log level
!> @param[in] b second log level
!> @return    .true. if a > b
!-----------------------------------------------------------------------
   pure logical function t_log_level_gt(a, b) result(res)
      type(t_log_level), intent(in) :: a, b

      res = a%level() > b%level()
   end function t_log_level_gt

!-----------------------------------------------------------------------
!> @brief Less than or equal comparison operator for log levels
!>
!> @param[in] a first log level
!> @param[in] b second log level
!> @return    .true. if a <= b
!-----------------------------------------------------------------------
   pure logical function t_log_level_lte(a, b) result(res)
      type(t_log_level), intent(in) :: a, b

      res = a%level() <= b%level()
   end function t_log_level_lte

!-----------------------------------------------------------------------
!> @brief Greater than or equal comparison operator for log levels
!>
!> @param[in] a first log level
!> @param[in] b second log level
!> @return    .true. if a >= b
!-----------------------------------------------------------------------
   pure logical function t_log_level_gte(a, b) result(res)
      type(t_log_level), intent(in) :: a, b

      res = a%level() >= b%level()
   end function t_log_level_gte

!-----------------------------------------------------------------------
!> @brief Equality comparison operator for log levels
!>
!> @param[in] a first log level
!> @param[in] b second log level
!> @return    .true. if a == b
!-----------------------------------------------------------------------
   pure logical function t_log_level_eq(a, b) result(res)
      type(t_log_level), intent(in) :: a, b

      res = a%level() == b%level()
   end function t_log_level_eq

!-----------------------------------------------------------------------
!> @brief Not equal comparison operator for log levels
!>
!> @param[in] a first log level
!> @param[in] b second log level
!> @return    .true. if a /= b
!-----------------------------------------------------------------------
   pure logical function t_log_level_neq(a, b) result(res)
      type(t_log_level), intent(in) :: a, b

      res = a%level() /= b%level()
   end function t_log_level_neq

!-----------------------------------------------------------------------
!> @brief Open the log file (fort.16)
!>
!> This must be called after make dirname so that we know where to put the log files
!>
!> @param[in] dir (optional) directory where log files will be created, defaults to current directory
!-----------------------------------------------------------------------
   subroutine openLogFile(dir)
      implicit none
      character(*), intent(in), optional :: dir
      character(len=256) :: localdir

      if (present(dir) .eqv. .false.) then
         localdir = '.'
      else
         localdir = dir
      end if

      open (log_unit, FILE=trim(localdir)//'/'//'fort.16', ACTION='WRITE', STATUS='REPLACE')

   end subroutine openLogFile

!-----------------------------------------------------------------------
!> @brief Write a message to the screen with specified logging level
!>
!> General purpose subroutine to write a message to the screen
!> subject to the user's selection of where to write screen output. The logging
!> level is controlled by NABOUT from the fort.15 file. The actual destination
!> of messages written to the screen is controlled by NSCREEN from the fort.15 file.
!> In parallel, only the processor with rank 0 actually writes the message.
!> This subroutine assumes that setMessageSource has been called at the beginning
!> of the subroutine that calls this one, and unsetMessageSource must be called at the end.
!>
!> @param[in] log_level    logging level (DEBUG, ECHO, INFO, WARNING, ERROR)
!> @param[in] message      text message to write
!-----------------------------------------------------------------------
   subroutine screenMessage(log_level, message)
      use iso_fortran_env, only: error_unit
      use sizes, only: myproc
      implicit none
      type(t_log_level), intent(in) :: log_level
      character(*), intent(in) :: message
      integer :: this_unit

      if (myproc == 0 .and. nscreen /= 0 .and. log_level >= nabout) then

         ! Check for the case we are writing to adcirc.log
         if (screenUnit == output_unit) then
            this_unit = log_level%unit()
         else
            this_unit = screenUnit
         end if

         write (this_unit, '(A)') formatLogMessage(log_level, message)
         flush (this_unit)
      end if
   end subroutine screenMessage

!-----------------------------------------------------------------------
!> @brief Write a message to the fort.16 log file
!>
!> General purpose subroutine to write a message to the fort.16 file.
!> In parallel, processors of all ranks will write the message to their own
!> subdomain fort.16 files. This subroutine assumes that setMessageSource has
!> been called at the beginning of the subroutine that calls this one, and
!> unsetMessageSource must be called at the end.
!>
!> @param[in] log_level   logging level (DEBUG, ECHO, INFO, WARNING, ERROR)
!> @param[in] message     text message to write
!-----------------------------------------------------------------------
   subroutine logMessage(log_level, message)
      implicit none
      type(t_log_level), intent(in) :: log_level
      character(*), intent(in) :: message

      if (log_level >= NABOUT) then
         write (log_unit, '(A)') formatLogMessage(log_level, message)
         flush (log_unit)
      end if
   end subroutine logMessage

!-----------------------------------------------------------------------
!> @brief Write a message to both the screen and the fort.16 log file
!>
!> @param[in] level   logging level (DEBUG, ECHO, INFO, WARNING, ERROR)
!> @param[in] message text message to write
!-----------------------------------------------------------------------
   subroutine allMessage(level, message)
      implicit none
      type(t_log_level), intent(in) :: level
      character(*), intent(in) :: message

      call screenMessage(level, message)
      call logMessage(level, message)
   end subroutine allMessage

!-----------------------------------------------------------------------
!> @brief Set the name of the subroutine that is writing messages
!>
!> This must be used at the start of any subroutine that calls
!> screenMessage, logMessage, or allMessage. Creates a call stack for message tracking.
!>
!> @param[in] source name of the calling subroutine
!-----------------------------------------------------------------------
   subroutine setMessageSource(source)
      implicit none
      character(*), intent(in) :: source

      sourceNumber = sourceNumber + 1
      messageSources(sourceNumber) = source
   end subroutine setMessageSource

!-----------------------------------------------------------------------
!> @brief Remove the name of the subroutine from the message source stack
!>
!> This must be used at the end of any subroutine that calls
!> screenMessage, logMessage, or allMessage. Pops the call stack for message tracking.
!-----------------------------------------------------------------------
   subroutine unsetMessageSource()
      implicit none

      sourceNumber = sourceNumber - 1
   end subroutine unsetMessageSource

end module mod_logging
