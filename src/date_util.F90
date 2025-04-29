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
!> @description The date utilities module contains simple functions that are useful
!> for working with dates.
!>
!> @author Zach Cobell
!> @date 2025-05-01
!-------------------------------------------------------------------------------!
module mod_date_util
   implicit none

   private

   public :: get_dt_string, get_dt_components

contains

   character(128) pure function get_dt_string(dt) result(dt_str)
      use datetime_module, only: timedelta
      implicit none
      type(timedelta), intent(in) :: dt

      if (dt%getDays() > 0) then
         write (dt_str, 2000) dt%getDays(), dt%getHours(), dt%getMinutes(), dt%getSeconds()
      else if (dt%getHours() > 0) then
         write (dt_str, 2001) dt%getHours(), dt%getMinutes(), dt%getSeconds()
      else if (dt%getMinutes() > 0) then
         write (dt_str, 2002) dt%getMinutes(), dt%getSeconds()
      else
         write (dt_str, 2003) dt%getSeconds()
      end if

2000  format(I0, 'd', I2.2, 'h', I2.2, 'm', I2.2, 's')
2001  format(I0, 'h', I2.2, 'm', I2.2, 's')
2002  format(I0, 'm', I2.2, 's')
2003  format(I0, 's')

   end function get_dt_string

   type(timedelta) pure function get_dt_components(current_time_in) result(dt)
      use datetime_module, only: timedelta
      implicit none
      real(8), intent(in) :: current_time_in
      real(8) :: current_time
      integer :: days, hours, minutes, seconds

      current_time = current_time_in

      ! Given the current time as seconds since start, compute the number of days
      ! and the time of day in hours, minutes, and seconds.
      days = int(current_time/86400d0)
      current_time = current_time - dble(days)*86400d0
      hours = int(current_time/3600d0)
      current_time = current_time - dble(hours)*3600d0
      minutes = int(current_time/60d0)
      current_time = current_time - dble(minutes)*60d0
      seconds = int(current_time)
      current_time = current_time - dble(seconds)

      dt = timedelta(days, hours, minutes, seconds)

   end function get_dt_components

end module mod_date_util
