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
!> @description The status module is responsible for the functions required to write
!> the model status to the screen and log file as necessary
!>
!> @author Zach Cobell
!> @date 2025-05-01
!-------------------------------------------------------------------------------!
module mod_adcirc_status
   implicit none

   private

   public :: print_model_status

contains

   !--------------------------------------------------------------------
   ! S U B R O U T I N E   P R I N T   M O D E L   S T A T U S
   !--------------------------------------------------------------------
   !> Prints the model status to the screen and to the log file.
   !> This subroutine is called at the end of each time step.
   !>
   !> @param myproc The rank of the processor.
   !> @param iteration The current model time iteration
   !> @param percent_complete The percentage of the current time step
   !> @param num_solver_it The number of solver iterations
   !> @param reference_time The reference time for the model run
   !> @param current_time The current model time in seconds
   !> @param max_elev The maximum elevation
   !> @param max_elev_pos The node number of the maximum elevation
   !> @param max_vel The maximum velocity
   !> @param max_vel_pos The node number of the maximum velocity
   !> @param warning_elev The elevation warning threshold
   !> @param error_elev The elevation error threshold
   !--------------------------------------------------------------------
   subroutine print_model_status(myproc, nscreen, screenUnit, &
                                 model_start_wallclock, num_solver_it, &
                                 iteration, num_hotstart_it, total_it, &
                                 reference_time, current_time, max_elev, &
                                 max_elev_pos, max_elev_rank, max_vel, max_vel_pos, &
                                 max_vel_rank, warning_elev, error_elev, &
                                 warning_vel, error_vel)
      use datetime_module, only: datetime, timedelta
      use mod_date_util, only: get_dt_components, get_dt_string
      implicit none
      integer, intent(in) :: myproc
      integer, intent(in) :: nscreen
      integer, intent(in) :: screenUnit
      real(8), intent(in) :: model_start_wallclock
      integer, intent(in) :: iteration
      integer, intent(in) :: num_solver_it
      integer, intent(in) :: num_hotstart_it
      integer, intent(in) :: total_it
      real(8), intent(in) :: current_time
      real(8), intent(in) :: max_elev
      integer, intent(in) :: max_elev_pos
      real(8), intent(in) :: max_vel
      integer, intent(in) :: max_vel_pos
      integer, intent(in) :: max_elev_rank
      integer, intent(in) :: max_vel_rank
      real(8), intent(in) :: warning_elev
      real(8), intent(in) :: error_elev
      real(8), intent(in) :: warning_vel
      real(8), intent(in) :: error_vel
      type(datetime), intent(in) :: reference_time
      real(8)             :: percent_complete
      character(1024)     :: append_message
      character(2048)     :: message1, message2, message3, message4, message5
      character(22)       :: date_str
      character(128)      :: dt_str
      character(256)      :: elapsed_dt_str, remaining_dt_str
      type(datetime)      :: current_date
      type(timedelta)     :: model_dt, elapsed_dt, remaining_dt

#ifdef ADCIRC_LEGACY_SCREEN_TIME
      logical, parameter :: use_legacy_time_format = .true.
#else
      logical, parameter :: use_legacy_time_format = .false.
#endif

      if (modulo(iteration, nscreen) /= 0 .and. max_elev < warning_elev .and. max_vel < warning_vel) return

      percent_complete = dble(iteration - num_hotstart_it)/dble(total_it - num_hotstart_it)*100d0

      if (max_elev >= error_elev) then
         append_message = "** ERROR: Elevation.gt.ErrorElev, ADCIRC stopping. **"
      elseif (max_elev >= warning_elev) then
         append_message = "** WARNING: Elevation.gt.WarnElev **"
      else
         append_message = ""
      end if

      if (max_vel >= error_vel) then
         append_message = trim(append_message)//"** ERROR: Velocity.gt.ErrorVel, ADCIRC stopping. **"
      elseif (max_vel >= warning_vel) then
         append_message = trim(append_message)//"** WARNING: Velocity.gt.WarnVel **"
      end if

      if (reference_time%isValid() .or. use_legacy_time_format) then
         current_date = reference_time + timedelta(0, 0, 0, int(floor(current_time)))
         model_dt = get_dt_components(current_time)
         date_str = current_date%strftime("%Y-%m-%dT%H:%M:%S")
         dt_str = get_dt_string(model_dt)
         write (message1, 1993) iteration, percent_complete, num_solver_it, date_str, trim(dt_str)
      else
         write (message1, 1992) iteration, percent_complete, num_solver_it, current_time
      end if

      call get_elapsed_remaining_time(percent_complete, model_start_wallclock, elapsed_dt, remaining_dt)

      write (message2, 1994) max_elev, max_elev_pos
      write (message3, 1995) max_vel, max_vel_pos
      write (message4, 1996) max_elev_rank
      write (message5, 1996) max_vel_rank

#ifdef CMPI
      message2 = trim(message2)//trim(message4)
      message3 = trim(message3)//trim(message5)
#endif

      elapsed_dt_str = "    ELAPSED WALL-CLOCK TIME = "//trim(get_dt_string(elapsed_dt))
      remaining_dt_str = "ESTIMATED REMAINING WALL-CLOCK TIME = "//trim(get_dt_string(remaining_dt))

      if (myproc == 0) then
         write (screenUnit, '(A)') trim(message1)
         write (screenUnit, '(A)') trim(message2)
         write (screenUnit, '(A)') trim(message3)
         write (screenUnit, '(A,",",1X,A)') trim(elapsed_dt_str), trim(remaining_dt_str)
         if (len_trim(append_message) > 0) then
            write (screenUnit, '(A)') trim(append_message)
         end if
      end if

      if (max_elev >= warning_elev) then
         write (16, '(A)') trim(message1)
         write (16, '(A)') trim(message2)
         write (16, '(A)') trim(message3)
         if (len_trim(append_message) > 0) then
            write (16, '(A)') trim(append_message)
         end if
      end if

1992  format(1x, 'TIME STEP =', I8, 1x, F6.2, '% COMPLETE', 2x, 'ITERATIONS =', I5, 2x, 'TIME = ', E15.8)
1993  format(1x, 'TIME STEP =', I8, 1x, F6.2, '% COMPLETE', 2x, 'ITERATIONS =', I5, 2x, 'TIME = ', A20, '(', A, ')')
#ifdef CMPI
1994  format(5x, 'WSELMAX = ', 1pe12.4e3, ' AT GLOBAL NODE ', I8)
1995  format(4x, 'SPEEDMAX = ', 1pe12.4e3, ' AT GLOBAL NODE ', I8)
#else
1994  format(5x, 'WSELMAX = ', 1pe12.4e3, ' AT NODE ', I8)
1995  format(4x, 'SPEEDMAX = ', 1pe12.4e3, ' AT NODE ', I8)
#endif
1996  format(1x, 'IN DOMAIN = ', I5)

   end subroutine print_model_status

   !--------------------------------------------------------------------
   !> Generates two timedelta objects containing the elapsed and remaining
   !> time estimate for the model run.
   !>
   !> @param[in] percent_complete The percent complete of the simulation
   !> @param[in] model_start_wallclock Time when the model started running
   !> @param[out] elapsed_dt The elapsed time since the model started
   !> @param[out] remaining_dt The estimated remaining time for the model
   !--------------------------------------------------------------------
   subroutine get_elapsed_remaining_time(percent_complete, model_start_wallclock, elapsed_dt, remaining_dt)
      use datetime_module, only: timedelta
      use mod_date_util, only: get_dt_components
      implicit none
      real(8), intent(in) :: percent_complete
      real(8), intent(in) :: model_start_wallclock
      type(timedelta), intent(out) :: elapsed_dt, remaining_dt
      integer :: current_wallclock_time_int
      integer :: system_clock_rate, system_clock_max
      real(8) :: elapsed_wallclock_time, remaining_wallclock_time, current_wallclock_time

      call system_clock(current_wallclock_time_int, system_clock_rate, system_clock_max)
      current_wallclock_time = dble(current_wallclock_time_int)/dble(system_clock_rate)

      ! Compute the estimated time remaining using the percent complete
      elapsed_wallclock_time = current_wallclock_time - model_start_wallclock
      remaining_wallclock_time = dble(elapsed_wallclock_time)*(100.0d0/percent_complete) - &
                                 dble(elapsed_wallclock_time)

      elapsed_dt = get_dt_components(elapsed_wallclock_time)
      remaining_dt = get_dt_components(remaining_wallclock_time)
   end subroutine get_elapsed_remaining_time

end module mod_adcirc_status
