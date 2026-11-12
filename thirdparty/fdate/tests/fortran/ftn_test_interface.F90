
module test_fortran_datetime
   use mod_datetime, only: t_datetime, t_timedelta, operator(+)
   implicit none

   public

contains

   subroutine parse_test(date_string)
      implicit none
      character(len=*), intent(in) :: date_string
      type(t_datetime) :: dt
      integer(kind=8) :: timestamp

      dt = t_datetime(date_string)
      timestamp = dt%timestamp()
      write (*, *) timestamp

   end subroutine parse_test

   subroutine arithmetic_test(date_string, days, hours, minutes, seconds)
      implicit none
      character(len=*), intent(in) :: date_string
      integer, intent(in) :: days, hours, minutes, seconds
      type(t_datetime) :: dt
      type(t_timedelta) :: ts

      dt = t_datetime(date_string)
      ts = t_timedelta(days, hours, minutes, seconds)
      dt = dt + ts
      write (*, *) dt%strftime("%Y-%m-%dT%H:%M:%S")

   end subroutine arithmetic_test

end module test_fortran_datetime

program test_interface

   use test_fortran_datetime, only: parse_test, arithmetic_test

   implicit none

   integer :: n_args
   character(len=100) :: date_string
   character(len=20) :: mode
   character(len=10) :: days_str, hours_str, minutes_str, seconds_str
   integer :: days, hours, minutes, seconds

   ! Get the number of command line arguments
   n_args = command_argument_count()
   if (n_args < 1) then
      write (*, *) "Usage: test_interface <mode> <date_string> OR <mode> <date_string> <days> <hours> <minutes> <seconds>"
      call exit(1)
   end if

   call get_command_argument(1, mode)
   if (mode == "parse") then
      call get_command_argument(2, date_string)
      call parse_test(date_string)
   elseif (mode == "arithmetic") then
      call get_command_argument(2, date_string)
      call get_command_argument(3, days_str)
      call get_command_argument(4, hours_str)
      call get_command_argument(5, minutes_str)
      call get_command_argument(6, seconds_str)
      read (days_str, *) days
      read (hours_str, *) hours
      read (minutes_str, *) minutes
      read (seconds_str, *) seconds
      call arithmetic_test(date_string, days, hours, minutes, seconds)
   else
      call exit(1)
   end if

end program test_interface
