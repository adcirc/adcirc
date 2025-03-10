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
module mod_logging

   ! Log levels, in order from largest amount of log messages
   ! written (DEBUG) to fewest log messages written (ERROR). Compared
   ! to the value of NABOUT when determining which messages to write
   ! to the screen or to log files.
   integer, parameter :: DEBUG = -1 ! write all messages and echo input
   integer, parameter :: ECHO = 0 ! echo input, plus write all non-debug
   integer, parameter :: INFO = 1 ! don't echo input; write all non-debug
   integer, parameter :: WARNING = 2 ! don't echo input; write only warn/err
   integer, parameter :: ERROR = 3 ! don't echo input; only fatal msgs
   character(len=10), dimension(5), parameter :: logLevelNames = &
                                                 (/"DEBUG  ", "ECHO   ", "INFO   ", "WARNING", "ERROR  "/)

   character(len=50), dimension(50) :: messageSources ! subroutine names
   character(len=1024) :: scratchMessage ! used for formatted messages
   character(len=1024) :: scratchFormat ! used for Fortran format strings
   integer :: sourceNumber ! index into messageSources for current sub
   integer :: nscreen
   integer :: nabout
   integer :: screenUnit = 6

   public :: initLogging, openLogFile, screenMessage, logMessage, &
             allMessage, setMessageSource, unsetMessageSource, &
             DEBUG, ECHO, INFO, WARNING, ERROR, scratchMessage, &
             scratchFormat, screenUnit

   private

contains

   !--------------------------------------------------------------------
   ! S U B R O U T I N E    I N I T   L O G G I N G
   !--------------------------------------------------------------------
   !> Initialize the names for the logging levels and the counter
   !> for the current subroutine.
   !>
   !> @param in_nscreen The value of NSCREEN from the fort.15 file.
   !> @param in_nabout The value of NABOUT from the fort.15 file.
   !--------------------------------------------------------------------
   subroutine initLogging(in_nscreen, in_nabout)
      implicit none

      integer, intent(in) :: in_nscreen, in_nabout

      nscreen = in_nscreen
      nabout = in_nabout
      sourceNumber = 0
   end subroutine initLogging
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! S U B R O U T I N E    O P E N   L O G   F I L E
   !--------------------------------------------------------------------
   !> Open the log file; this must be called after make dirname
   !> so that we know where to put this log file.
   !--------------------------------------------------------------------
   subroutine openLogFile()

#ifdef CMPI
      use sizes, only: localDir
#endif

      implicit none
      !...
      !...  OPEN STATEMENT FOR UNIT 16 OUTPUT FILE (ADCIRC LOG FILE)
      !...
      !     cms51.06: moved fort.33 file from read_input.F to here
#ifdef CMPI
      open (16, FILE=trim(localdir)//'/'//'fort.16', ACTION='WRITE', &
            STATUS='REPLACE')
      open (33, FILE=trim(localdir)//'/'//'fort.33', ACTION='WRITE', &
            STATUS='REPLACE')
#else
      open (16, FILE='fort.16', ACTION='WRITE', STATUS='REPLACE')
      open (33, FILE='fort.33', ACTION='WRITE', STATUS='REPLACE')
#endif

   end subroutine openLogFile
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! S U B R O U T I N E    S C R E E N
   !--------------------------------------------------------------------
   !> General purpose subroutine to write a message to
   !> the screen with a certain "logging level", and subject to the
   !> user's selection of where to write screen output. The logging
   !> level is controlled by NABOUT from the fort.15 file. The actual
   !> destination of messages written to the screen is controlled by
   !> NSCREEN from the fort.15 file.
   !>
   !> In parallel, only the processor with rank 0 actually writes
   !> the message.
   !> This subroutine assumes that the global variable "caller" has
   !> been set to the name of the subroutine calling it. Therefore,
   !> the setMessageSource subroutine must be called at the beginning
   !> of the subroutine that calls this one, and unsetMessageSource
   !> must be called at the end.
   !>
   !> @param level The logging level of the message.
   !> @param message The message to write to the screen.
   !--------------------------------------------------------------------
   subroutine screenMessage(level, message)
      use SIZES, only: myProc

      implicit none

      integer, intent(in) :: level
      character(*), intent(in) :: message

#ifdef FULL_STACK
      integer :: j ! loop counter for stack
#endif

      if (myProc == 0) then
         if (NSCREEN /= 0) then
            if (level >= NABOUT) then
#ifdef FULL_STACK
               write (screenUnit, 331, advance="no") &
                  trim(logLevelNames(level + 2)), &
                  (trim(messageSources(j)), j=1, sourceNumber)
               write (screenUnit, 332) trim(message)
331            format(A, ": ", A, 50(:, "->", A))
332            format(": ", A)
#else
               write (screenUnit, 333) trim(logLevelNames(level + 2)), &
                  trim(messageSources(sourceNumber)), trim(message)
333            format(A, ": ", A, ": ", A)
#endif
               !       !#ifdef FLUSH_MESSAGES
               flush (screenUnit)
               !               !#endif
            end if
         end if
      end if

   end subroutine screenMessage
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! S U B R O U T I N E    L O G   M E S S A G E
   !--------------------------------------------------------------------
   !> General purpose subroutine to write a message to
   !> the fort.16 file. In parallel, processors of all ranks will
   !> write the message to their own subdomain fort.16 files.
   !>
   !> This subroutine assumes that the global variable "caller" has
   !> been set to the name of the subroutine calling it. Therefore,
   !> the setMessageSource subroutine must be called at the beginning
   !> of the subroutine that calls this one, and unsetMessageSource
   !> must be called at the end.
   !>
   !> @param level The logging level of the message.
   !> @param message The message to write to the log file.
   !--------------------------------------------------------------------
   subroutine logMessage(level, message)
      implicit none

      integer, intent(in) :: level
      character(*), intent(in) :: message

#ifdef FULL_STACK
      integer :: j ! loop counter for stack
#endif

      if (level >= NABOUT) then
#ifdef FULL_STACK
         write (16, 331, advance="no") trim(logLevelNames(level + 2)), &
            (trim(messageSources(j)), j=1, sourceNumber)
         write (16, 332) trim(message)
331      format(A, ": ", A, 50(:, "->", A))
332      format(": ", A)
#else
         write (16, 333) trim(logLevelNames(level + 2)), &
            trim(messageSources(sourceNumber)), trim(message)
333      format(A, ": ", A, ": ", A)
#endif
         flush (16)
      end if

   end subroutine logMessage
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! S U B R O U T I N E   A L L    M E S S A G E
   !--------------------------------------------------------------------
   !> General purpose subroutine to write a message to
   !> both the screen and to the fort.16 log file.
   !>
   !> @param level The logging level of the message.
   !> @param message The message to write to the screen and log file.
   !--------------------------------------------------------------------
   subroutine allMessage(level, message)
      implicit none

      integer, intent(in) :: level
      character(*), intent(in) :: message

      call screenMessage(level, message)
      call logMessage(level, message)

   end subroutine allMessage
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! S U B R O U T I N E   S E T   M E S S A G E   S O U R C E
   !--------------------------------------------------------------------
   !> Sets the name of the subroutine that is writing
   !> log and/or screen messages. Must use at the start of any subroutine
   !> that calls screen, logMessage, or allMessage.
   !>
   !> @param source The name of the subroutine writing log and/or screen messages
   !--------------------------------------------------------------------
   subroutine setMessageSource(source)
      implicit none

      character(*), intent(in) :: source

      sourceNumber = sourceNumber + 1
      messageSources(sourceNumber) = source

   end subroutine setMessageSource
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! S U B R O U T I N E   U N S E T   M E S S A G E   S O U R C E
   !--------------------------------------------------------------------
   !> Removes the name of the subroutine that is no longer
   !> writing log and/or screen messages. Must use at the end of
   !> any subroutine that calls screen, logMessage, or allMessage.
   !--------------------------------------------------------------------
   subroutine unsetMessageSource()
      implicit none

      sourceNumber = sourceNumber - 1

   end subroutine unsetMessageSource
   !--------------------------------------------------------------------

end module mod_logging
