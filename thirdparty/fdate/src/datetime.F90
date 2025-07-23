!>
!> FDate - A Fortran Date and Time Library based on C++
!> Copyright (C) 2025 Zach Cobell
!>
!> This program is free software: you can redistribute it and/or modify
!> it under the terms of the GNU General Public License as published by
!> the Free Software Foundation, either version 3 of the License, or
!> (at your option) any later version.
!>
!> This program is distributed in the hope that it will be useful,
!> but WITHOUT ANY WARRANTY; without even the implied warranty of
!> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!> GNU General Public License for more details.
!>
!> You should have received a copy of the GNU General Public License
!> along with this program.  If not, see <https://www.gnu.org/licenses/>.
!>
!> @file datetime_mod.f90
!> @brief Fortran module for datetime handling
!>
!> This module provides a Fortran interface to the C++ datetime library.
!> It defines two derived types, TimeDelta and DateTime, which are implemented
!> using 64-bit integers representing milliseconds. This approach avoids
!> memory management issues while providing a clean Fortran API.

module mod_datetime
   use, intrinsic :: iso_c_binding, only: c_int, c_int64_t, c_char, c_null_char, c_bool, c_ptr, c_null_ptr, c_loc, c_double
   implicit none

   private

   integer, parameter :: DATETIME_STRING_BUFFER_SIZE = 64
   integer, parameter :: TIMEDELTA_STRING_BUFFER_SIZE = 64
   integer, parameter :: DATETIME_MAX_STRING_LENGTH = 10000
   integer, parameter :: DATETIME_MAX_FORMAT_COUNT = 1000

   !> @brief A span of time with various components
   !>
   !> TimeDelta represents a duration that can be expressed in terms of days,
   !> hours, minutes, seconds, and milliseconds. It is stored internally as
   !> an integer number of milliseconds.
   type :: t_timedelta
      private
      integer(kind=c_int64_t) :: ms_count !< Total milliseconds in the timedelta
   contains
      !> @brief Get the days component
      procedure :: days => timedelta_days
      !> @brief Get the hours component
      procedure :: hours => timedelta_hours
      !> @brief Get the minutes component
      procedure :: minutes => timedelta_minutes
      !> @brief Get the seconds component
      procedure :: seconds => timedelta_seconds
      !> @brief Get the milliseconds component
      procedure :: milliseconds => timedelta_milliseconds
      !> @brief Get the total days
      procedure :: total_days => timedelta_total_days
      !> @brief Get the total hours
      procedure :: total_hours => timedelta_total_hours
      !> @brief Get the total minutes
      procedure :: total_minutes => timedelta_total_minutes
      !> @brief Get the total seconds
      procedure :: total_seconds => timedelta_total_seconds
      !> @brief Get the total milliseconds
      procedure :: total_milliseconds => timedelta_total_milliseconds
      !> @brief Convert to a string representation
      procedure :: to_string => timedelta_to_string
   end type t_timedelta

   !> @brief A point in time
   !>
   !> DateTime represents a specific point in time, with year, month, day,
   !> hour, minute, second, and millisecond components. It is stored internally
   !> as milliseconds since the Unix epoch (January 1, 1970, 00:00:00 UTC).
   type :: t_datetime
      private
      integer(kind=c_int64_t) :: timestamp_ms !< Milliseconds since epoch
   contains
      !> @brief Get the year component
      procedure :: year => datetime_year
      !> @brief Get the month component
      procedure :: month => datetime_month
      !> @brief Get the day component
      procedure :: day => datetime_day
      !> @brief Get the hour component
      procedure :: hour => datetime_hour
      !> @brief Get the minute component
      procedure :: minute => datetime_minute
      !> @brief Get the second component
      procedure :: second => datetime_second
      !> @brief Get the millisecond component
      procedure :: millisecond => datetime_millisecond
      !> @brief Get the Julian Day Number
      procedure :: julian_day_number => datetime_julian_day_number
      !> @brief Get the Julian Day (with fractional part)
      procedure :: julian_day => datetime_julian_day
      !> @brief Get the Julian Century (JC)
      procedure :: julian_century => datetime_julian_century
      !> @brief Get the timestamp (milliseconds since epoch)
      procedure :: timestamp => datetime_timestamp
      !> @brief Format the datetime to a string
      procedure :: strftime => datetime_strftime
      !> @brief Convert to ISO 8601 string format
      procedure :: to_iso_string => datetime_to_iso_string
      !> @brief Check if the datetime object is valid
      procedure :: valid => datetime_is_valid
   end type t_datetime

   ! Create a null datetime object
   type(t_datetime), parameter :: null_datetime = t_datetime(0_c_int64_t)

   ! Interface blocks for constructors
   !> @brief Constructor interface for timedelta
   interface t_timedelta
      module procedure :: timedelta_components
   end interface t_timedelta

   !> @brief Constructor interface for datetime
   interface t_datetime
      module procedure :: datetime_default
      module procedure :: datetime_ymd
      module procedure :: datetime_ymd_hms
      module procedure :: datetime_complete
      module procedure :: datetime_from_timestamp
      module procedure :: datetime_strptime
      module procedure :: datetime_strptime_auto
      module procedure :: datetime_strptime_with_formats
   end interface t_datetime

   ! Operator interfaces
   !> @brief Addition operator
   interface operator(+)
      module procedure :: timedelta_add_timedelta
      module procedure :: datetime_add_timedelta
   end interface operator(+)

   !> @brief Subtraction operator
   interface operator(-)
      module procedure :: timedelta_subtract_timedelta
      module procedure :: datetime_subtract_timedelta
      module procedure :: datetime_difference
   end interface operator(-)

   !> @brief Multiplication operator
   interface operator(*)
      module procedure :: timedelta_multiply
   end interface operator(*)

   !> @brief Division operator
   interface operator(/)
      module procedure :: timedelta_divide
   end interface operator(/)

   !> @brief Equality operator
   interface operator(==)
      module procedure :: timedelta_equals
      module procedure :: datetime_equals
   end interface operator(==)

   !> @brief Inequality operator
   interface operator(/=)
      module procedure :: timedelta_not_equals
      module procedure :: datetime_not_equals
   end interface operator(/=)

   !> @brief Less than operator
   interface operator(<)
      module procedure :: timedelta_less_than
      module procedure :: datetime_less_than
   end interface operator(<)

   !> @brief Greater than operator
   interface operator(>)
      module procedure :: timedelta_greater_than
      module procedure :: datetime_greater_than
   end interface operator(>)

   !> @brief Less than or equal operator
   interface operator(<=)
      module procedure :: timedelta_less_equal
      module procedure :: datetime_less_equal
   end interface operator(<=)

   !> @brief Greater than or equal operator
   interface operator(>=)
      module procedure :: timedelta_greater_equal
      module procedure :: datetime_greater_equal
   end interface operator(>=)

   ! C function interfaces
   interface
      !> @brief Create a TimeDelta with components
      pure function f_timedelta_create(days, hours, minutes, seconds, milliseconds) &
         result(ts_ms) bind(C, name="f_timedelta_create")
         import :: c_int, c_int64_t
         implicit none
         integer(c_int), intent(in), value :: days, hours, minutes, seconds, milliseconds
         integer(c_int64_t) :: ts_ms
      end function f_timedelta_create

      !> @brief Get days component from a TimeDelta
      pure function f_timedelta_get_days(ts_ms) result(days) bind(C, name="f_timedelta_get_days")
         import :: c_int, c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: ts_ms
         integer(c_int64_t) :: days
      end function f_timedelta_get_days

      !> @brief Get hours component from a TimeDelta
      pure function f_timedelta_get_hours(ts_ms) result(hours) bind(C, name="f_timedelta_get_hours")
         import :: c_int, c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: ts_ms
         integer(c_int64_t) :: hours
      end function f_timedelta_get_hours

      !> @brief Get minutes component from a TimeDelta
      pure function f_timedelta_get_minutes(ts_ms) result(minutes) bind(C, name="f_timedelta_get_minutes")
         import :: c_int, c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: ts_ms
         integer(c_int64_t) :: minutes
      end function f_timedelta_get_minutes

      !> @brief Get seconds component from a TimeDelta
      pure function f_timedelta_get_seconds(ts_ms) result(seconds) bind(C, name="f_timedelta_get_seconds")
         import :: c_int, c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: ts_ms
         integer(c_int64_t) :: seconds
      end function f_timedelta_get_seconds

      !> @brief Get milliseconds component from a TimeDelta
      pure function f_timedelta_get_milliseconds(ts_ms) result(ms) bind(C, name="f_timedelta_get_milliseconds")
         import :: c_int, c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: ts_ms
         integer(c_int64_t) :: ms
      end function f_timedelta_get_milliseconds

      !> @brief Get total days from a TimeDelta
      pure function f_timedelta_get_total_days(ts_ms) result(days) bind(C, name="f_timedelta_get_total_days")
         import :: c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: ts_ms
         integer(c_int64_t) :: days
      end function f_timedelta_get_total_days

      !> @brief Get total hours from a TimeDelta
      pure function f_timedelta_get_total_hours(ts_ms) result(hours) bind(C, name="f_timedelta_get_total_hours")
         import :: c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: ts_ms
         integer(c_int64_t) :: hours
      end function f_timedelta_get_total_hours

      !> @brief Get total minutes from a TimeDelta
      pure function f_timedelta_get_total_minutes(ts_ms) result(minutes) bind(C, name="f_timedelta_get_total_minutes")
         import :: c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: ts_ms
         integer(c_int64_t) :: minutes
      end function f_timedelta_get_total_minutes

      !> @brief Get total seconds from a TimeDelta
      pure function f_timedelta_get_total_seconds(ts_ms) result(seconds) bind(C, name="f_timedelta_get_total_seconds")
         import :: c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: ts_ms
         integer(c_int64_t) :: seconds
      end function f_timedelta_get_total_seconds

      !> @brief Add two TimeDeltas
      pure function f_timedelta_add(ts1_ms, ts2_ms) result(sum_ms) bind(C, name="f_timedelta_add")
         import :: c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: ts1_ms, ts2_ms
         integer(c_int64_t) :: sum_ms
      end function f_timedelta_add

      !> @brief Subtract one TimeDelta from another
      pure function f_timedelta_subtract(ts1_ms, ts2_ms) result(diff_ms) bind(C, name="f_timedelta_subtract")
         import :: c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: ts1_ms, ts2_ms
         integer(c_int64_t) :: diff_ms
      end function f_timedelta_subtract

      !> @brief Multiply a TimeDelta by a factor
      pure function f_timedelta_multiply(ts_ms, factor) result(result_ms) bind(C, name="f_timedelta_multiply")
         import :: c_int, c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: ts_ms
         integer(c_int), intent(in), value :: factor
         integer(c_int64_t) :: result_ms
      end function f_timedelta_multiply

      !> @brief Divide a TimeDelta by a divisor
      pure function f_timedelta_divide(ts_ms, divisor) result(result_ms) bind(C, name="f_timedelta_divide")
         import :: c_int, c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: ts_ms
         integer(c_int), intent(in), value :: divisor
         integer(c_int64_t) :: result_ms
      end function f_timedelta_divide

      !> @brief Convert a TimeDelta to a string
      subroutine f_timedelta_to_string(ts_ms, buffer, buffer_size) bind(C, name="f_timedelta_to_string")
         import :: c_int64_t, c_int, c_char
         implicit none
         integer(c_int64_t), intent(in), value :: ts_ms
         integer(c_int), intent(in), value :: buffer_size
         character(kind=c_char), intent(in) :: buffer(buffer_size)
      end subroutine f_timedelta_to_string

      !> @brief Check if two TimeDeltas are equal
      pure function f_timedelta_equals(ts1_ms, ts2_ms) result(result) bind(C, name="f_timedelta_equals")
         import :: c_int, c_int64_t, c_bool
         implicit none
         integer(c_int64_t), intent(in), value :: ts1_ms, ts2_ms
         logical(c_bool) :: result
      end function f_timedelta_equals

      !> @brief Check if one TimeDelta is less than another
      pure function f_timedelta_less_than(ts1_ms, ts2_ms) result(result) bind(C, name="f_timedelta_less_than")
         import :: c_int64_t, c_bool
         implicit none
         integer(c_int64_t), intent(in), value :: ts1_ms, ts2_ms
         logical(c_bool) :: result
      end function f_timedelta_less_than

      !> @brief Check if one TimeDelta is greater than another
      pure function f_timedelta_greater_than(ts1_ms, ts2_ms) result(result) bind(C, name="f_timedelta_greater_than")
         import :: c_int64_t, c_bool
         implicit none
         integer(c_int64_t), intent(in), value :: ts1_ms, ts2_ms
         logical(c_bool) :: result
      end function f_timedelta_greater_than

      !> @brief Check if one TimeDelta is less than or equal to another
      pure function f_timedelta_less_equal(ts1_ms, ts2_ms) result(result) bind(C, name="f_timedelta_less_equal")
         import :: c_int64_t, c_bool
         implicit none
         integer(c_int64_t), intent(in), value :: ts1_ms, ts2_ms
         logical(c_bool) :: result
      end function f_timedelta_less_equal

      !> @brief Check if one TimeDelta is greater than or equal to another
      pure function f_timedelta_greater_equal(ts1_ms, ts2_ms) result(result) bind(C, name="f_timedelta_greater_equal")
         import :: c_int64_t, c_bool
         implicit none
         integer(c_int64_t), intent(in), value :: ts1_ms, ts2_ms
         logical(c_bool) :: result
      end function f_timedelta_greater_equal

      ! DateTime C functions
      !> @brief Create a DateTime from components
      pure function f_datetime_create(year, month, day, hour, minute, second, millisecond) &
         result(dt_ms) bind(C, name="f_datetime_create")
         import :: c_int, c_int64_t
         implicit none
         integer(c_int), intent(in), value :: year, month, day, hour, minute, second, millisecond
         integer(c_int64_t) :: dt_ms
      end function f_datetime_create

      !> @brief Get the current DateTime
      function f_datetime_now() result(dt_ms) bind(C, name="f_datetime_now")
         import :: c_int64_t
         implicit none
         integer(c_int64_t) :: dt_ms
      end function f_datetime_now

      !> @brief Parse a DateTime from a string
      function f_datetime_strptime(str, fmt, str_len, format_len) result(dt_ms) &
         bind(C, name="f_datetime_strptime")
         import :: c_int64_t, c_char, c_int
         implicit none
         integer(c_int), intent(in), value :: str_len, format_len
         character(kind=c_char), intent(in) :: fmt(format_len)
         character(kind=c_char), intent(in) :: str(str_len)
         integer(c_int64_t) :: dt_ms
      end function f_datetime_strptime

      !> @brief Get the year component from a DateTime
      pure function f_datetime_get_year(dt_ms) result(year) bind(C, name="f_datetime_get_year")
         import :: c_int, c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: dt_ms
         integer(c_int) :: year
      end function f_datetime_get_year

      !> @brief Get the month component from a DateTime
      pure function f_datetime_get_month(dt_ms) result(month) bind(C, name="f_datetime_get_month")
         import :: c_int, c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: dt_ms
         integer(c_int) :: month
      end function f_datetime_get_month

      !> @brief Get the day component from a DateTime
      pure function f_datetime_get_day(dt_ms) result(day) bind(C, name="f_datetime_get_day")
         import :: c_int, c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: dt_ms
         integer(c_int) :: day
      end function f_datetime_get_day

      !> @brief Get the hour component from a DateTime
      pure function f_datetime_get_hour(dt_ms) result(hour) bind(C, name="f_datetime_get_hour")
         import :: c_int, c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: dt_ms
         integer(c_int) :: hour
      end function f_datetime_get_hour

      !> @brief Get the minute component from a DateTime
      pure function f_datetime_get_minute(dt_ms) result(minute) bind(C, name="f_datetime_get_minute")
         import :: c_int, c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: dt_ms
         integer(c_int) :: minute
      end function f_datetime_get_minute

      !> @brief Get the second component from a DateTime
      pure function f_datetime_get_second(dt_ms) result(second) bind(C, name="f_datetime_get_second")
         import :: c_int, c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: dt_ms
         integer(c_int) :: second
      end function f_datetime_get_second

      !> @brief Get the millisecond component from a DateTime
      pure function f_datetime_get_millisecond(dt_ms) result(ms) bind(C, name="f_datetime_get_millisecond")
         import :: c_int, c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: dt_ms
         integer(c_int) :: ms
      end function f_datetime_get_millisecond

      !> @brief Get the Julian Day Number from a DateTime
      pure function f_datetime_get_julian_day_number(dt_ms) result(jdn) bind(C, name="f_datetime_get_julian_day_number")
         import :: c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: dt_ms
         integer(c_int64_t) :: jdn
      end function f_datetime_get_julian_day_number

      !> @brief Get the Julian Day (with fractional part) from a DateTime
      pure function f_datetime_get_julian_day(dt_ms) result(jd) bind(C, name="f_datetime_get_julian_day")
         import :: c_int64_t, c_double
         implicit none
         integer(c_int64_t), intent(in), value :: dt_ms
         real(c_double) :: jd
      end function f_datetime_get_julian_day

      !> @brief Get the Julian Century from a DateTime
      pure function f_datetime_get_julian_century(dt_ms) result(jc) bind(C, name="f_datetime_get_julian_century")
         import :: c_int64_t, c_double
         implicit none
         integer(c_int64_t), intent(in), value :: dt_ms
         real(c_double) :: jc
      end function f_datetime_get_julian_century

      !> @brief Add a TimeDelta to a DateTime
      pure function f_datetime_add_timedelta(dt_ms, ts_ms) result(result_ms) &
         bind(C, name="f_datetime_add_timedelta")
         import :: c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: dt_ms, ts_ms
         integer(c_int64_t) :: result_ms
      end function f_datetime_add_timedelta

      !> @brief Subtract a TimeDelta from a DateTime
      pure function f_datetime_subtract_timedelta(dt_ms, ts_ms) result(result_ms) &
         bind(C, name="f_datetime_subtract_timedelta")
         import :: c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: dt_ms, ts_ms
         integer(c_int64_t) :: result_ms
      end function f_datetime_subtract_timedelta

      !> @brief Calculate the difference between two DateTimes
      pure function f_datetime_difference(dt1_ms, dt2_ms) result(ts_ms) &
         bind(C, name="f_datetime_difference")
         import :: c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: dt1_ms, dt2_ms
         integer(c_int64_t) :: ts_ms
      end function f_datetime_difference

      !> @brief Format a DateTime to a string
      subroutine f_datetime_strftime(dt_ms, format_str, buffer, format_len, buffer_size) &
         bind(C, name="f_datetime_strftime")
         import :: c_int64_t, c_char, c_int
         implicit none
         integer(c_int64_t), intent(in), value :: dt_ms
         integer(c_int), intent(in), value :: format_len, buffer_size
         character(kind=c_char), intent(in) :: format_str(format_len)
         character(kind=c_char), intent(inout) :: buffer(buffer_size)
      end subroutine f_datetime_strftime

      !> @brief Get the timestamp (milliseconds since epoch) from a DateTime
      subroutine f_datetime_strftime_ms(dt_ms, format_str, buffer, format_len, buffer_size) &
         bind(C, name="f_datetime_strftime_milliseconds")
         import :: c_int64_t, c_char, c_int
         implicit none
         integer(c_int64_t), intent(in), value :: dt_ms
         integer(c_int), intent(in), value :: format_len, buffer_size
         character(kind=c_char), intent(in) :: format_str(format_len)
         character(kind=c_char), intent(inout) :: buffer(buffer_size)
      end subroutine f_datetime_strftime_ms

      !> @brief Convert a DateTime to an ISO 8601 string
      subroutine f_datetime_to_iso_string(dt_ms, buffer, buffer_size) &
         bind(C, name="f_datetime_to_iso_string")
         import :: c_int64_t, c_char, c_int
         implicit none
         integer(c_int64_t), intent(in), value :: dt_ms
         integer(c_int), intent(in), value :: buffer_size
         character(kind=c_char), intent(inout) :: buffer(buffer_size)
      end subroutine f_datetime_to_iso_string

      !> @brief Check if two DateTimes are equal
      pure function f_datetime_equals(dt1_ms, dt2_ms) result(result) bind(C, name="f_datetime_equals")
         import :: c_int64_t, c_bool
         implicit none
         integer(c_int64_t), intent(in), value :: dt1_ms, dt2_ms
         logical(c_bool) :: result
      end function f_datetime_equals

      !> @brief Check if one DateTime is less than another
      pure function f_datetime_less_than(dt1_ms, dt2_ms) result(result) bind(C, name="f_datetime_less_than")
         import :: c_bool, c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: dt1_ms, dt2_ms
         logical(c_bool) :: result
      end function f_datetime_less_than

      !> @brief Check if one DateTime is greater than another
      pure function f_datetime_greater_than(dt1_ms, dt2_ms) result(result) bind(C, name="f_datetime_greater_than")
         import :: c_bool, c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: dt1_ms, dt2_ms
         logical(c_bool) :: result
      end function f_datetime_greater_than

      !> @brief Check if one DateTime is less than or equal to another
      pure function f_datetime_less_equal(dt1_ms, dt2_ms) result(result) bind(C, name="f_datetime_less_equal")
         import :: c_bool, c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: dt1_ms, dt2_ms
         logical(c_bool) :: result
      end function f_datetime_less_equal

      !> @brief Check if one DateTime is greater than or equal to another
      pure function f_datetime_greater_equal(dt1_ms, dt2_ms) result(result) bind(C, name="f_datetime_greater_equal")
         import :: c_bool, c_int64_t
         implicit none
         integer(c_int64_t), intent(in), value :: dt1_ms, dt2_ms
         logical(c_bool) :: result
      end function f_datetime_greater_equal

      !> @brief Check if a DateTime object is valid
      pure function f_datetime_is_valid(dt_ms) result(valid) bind(C, name="f_datetime_is_valid")
         import :: c_int64_t, c_bool
         implicit none
         integer(c_int64_t), intent(in), value :: dt_ms
         logical(c_bool) :: valid
      end function f_datetime_is_valid

      !> @brief Interface to get the invalid DateTime constant
      pure function f_datetime_invalid_timestamp() result(invalid_ts) bind(C, name="f_datetime_invalid_timestamp")
         import :: c_int64_t
         implicit none
         integer(c_int64_t) :: invalid_ts
      end function f_datetime_invalid_timestamp

      !> @brief Parse a DateTime from a string using an array of format options
      function f_datetime_strptime_with_formats(str, formats, str_len, formats_len, num_formats) &
         result(dt_ms) bind(C, name="f_datetime_strptime_with_formats")
         import :: c_int64_t, c_char, c_int, c_ptr
         implicit none
         integer(c_int), intent(in), value :: str_len, num_formats
         character(kind=c_char), intent(in) :: str(str_len)
         type(c_ptr), intent(in), value :: formats
         integer(c_int), intent(in) :: formats_len(num_formats)
         integer(c_int64_t) :: dt_ms
      end function f_datetime_strptime_with_formats

      !> @brief Parse a DateTime with automatic detection and fallback formats
      function f_datetime_strptime_auto_with_fallback(str, fallback_formats, str_len, formats_len, num_formats) &
         result(dt_ms) bind(C, name="f_datetime_strptime_auto_with_fallback")
         import :: c_int64_t, c_char, c_int, c_ptr
         implicit none
         integer(c_int), intent(in), value :: str_len, num_formats
         character(kind=c_char), intent(in) :: str(str_len)
         type(c_ptr), intent(in), value :: fallback_formats
         integer(c_int), intent(in) :: formats_len(num_formats)
         integer(c_int64_t) :: dt_ms
      end function f_datetime_strptime_auto_with_fallback
   end interface

   public :: t_timedelta, t_datetime, now, null_datetime
   public :: operator(+), operator(-), operator(*), operator(/), operator(==)
   public :: operator(/=), operator(<), operator(>), operator(<=), operator(>=)
   public :: datetime_strptime_auto_with_fallback

contains

   !===========================================================================
   ! TimeDelta implementations
   !===========================================================================

   !> @brief Create a TimeDelta from days, hours, minutes, seconds, and milliseconds
   !> @param days Number of days
   !> @param hours Number of hours
   !> @param minutes Number of minutes
   !> @param seconds Number of seconds
   !> @param milliseconds Number of milliseconds
   !> @return TimeDelta with the specified duration
   pure function timedelta_components(days, hours, minutes, seconds, milliseconds) result(ts)
      implicit none
      integer, intent(in), optional :: days, hours, minutes, seconds, milliseconds
      type(t_timedelta) :: ts
      integer :: days_int, hours_int, minutes_int, seconds_int, milliseconds_int

      ! Default values
      if (.not. present(days)) then
         days_int = 0
      else
         days_int = days
      end if

      if (.not. present(hours)) then
         hours_int = 0
      else
         hours_int = hours
      end if

      if (.not. present(minutes)) then
         minutes_int = 0
      else
         minutes_int = minutes
      end if

      if (.not. present(seconds)) then
         seconds_int = 0
      else
         seconds_int = seconds
      end if

      if (.not. present(milliseconds)) then
         milliseconds_int = 0
      else
         milliseconds_int = milliseconds
      end if

      ts%ms_count = f_timedelta_create(days_int, hours_int, minutes_int, seconds_int, milliseconds_int)
   end function timedelta_components

   !> @brief Get the days component from a TimeDelta
   !> @param this TimeDelta object
   !> @return Days component
   pure function timedelta_days(this) result(days)
      implicit none
      class(t_timedelta), intent(in) :: this
      integer(c_int64_t) :: days64
      integer :: days
      days64 = f_timedelta_get_days(this%ms_count)
      days = int(days64)
   end function timedelta_days

   !> @brief Get the hours component from a TimeDelta
   !> @param this TimeDelta object
   !> @return Hours component
   pure function timedelta_hours(this) result(hours)
      implicit none
      class(t_timedelta), intent(in) :: this
      integer(c_int64_t) :: hours64
      integer :: hours
      hours64 = f_timedelta_get_hours(this%ms_count)
      hours = int(hours64)
   end function timedelta_hours

   !> @brief Get the minutes component from a TimeDelta
   !> @param this TimeDelta object
   !> @return Minutes component
   pure function timedelta_minutes(this) result(minutes)
      implicit none
      class(t_timedelta), intent(in) :: this
      integer(c_int64_t) :: minutes64
      integer :: minutes
      minutes64 = f_timedelta_get_minutes(this%ms_count)
      minutes = int(minutes64)
   end function timedelta_minutes

   !> @brief Get the seconds component from a TimeDelta
   !> @param this TimeDelta object
   !> @return Seconds component
   pure function timedelta_seconds(this) result(seconds)
      implicit none
      class(t_timedelta), intent(in) :: this
      integer(c_int64_t) :: seconds64
      integer :: seconds
      seconds64 = f_timedelta_get_seconds(this%ms_count)
      seconds = int(seconds64)
   end function timedelta_seconds

   !> @brief Get the milliseconds component from a TimeDelta
   !> @param this TimeDelta object
   !> @return Milliseconds component
   pure function timedelta_milliseconds(this) result(ms)
      implicit none
      class(t_timedelta), intent(in) :: this
      integer(c_int64_t) :: ms64
      integer :: ms
      ms64 = f_timedelta_get_milliseconds(this%ms_count)
      ms = int(ms64)
   end function timedelta_milliseconds

   !> @brief Get the total days from a TimeDelta
   !> @param this TimeDelta object
   !> @return Total days
   pure function timedelta_total_days(this) result(days)
      implicit none
      class(t_timedelta), intent(in) :: this
      integer(kind=8) :: days
      days = f_timedelta_get_total_days(this%ms_count)
   end function timedelta_total_days

   !> @brief Get the total hours from a TimeDelta
   !> @param this TimeDelta object
   !> @return Total hours
   pure function timedelta_total_hours(this) result(hours)
      implicit none
      class(t_timedelta), intent(in) :: this
      integer(kind=8) :: hours
      hours = f_timedelta_get_total_hours(this%ms_count)
   end function timedelta_total_hours

   !> @brief Get the total minutes from a TimeDelta
   !> @param this TimeDelta object
   !> @return Total minutes
   pure function timedelta_total_minutes(this) result(minutes)
      implicit none
      class(t_timedelta), intent(in) :: this
      integer(kind=8) :: minutes
      minutes = f_timedelta_get_total_minutes(this%ms_count)
   end function timedelta_total_minutes

   !> @brief Get the total seconds from a TimeDelta
   !> @param this TimeDelta object
   !> @return Total seconds
   pure function timedelta_total_seconds(this) result(seconds)
      implicit none
      class(t_timedelta), intent(in) :: this
      integer(kind=8) :: seconds
      seconds = f_timedelta_get_total_seconds(this%ms_count)
   end function timedelta_total_seconds

   !> @brief Get the total milliseconds from a TimeDelta
   !> @param this TimeDelta object
   !> @return Total milliseconds
   pure function timedelta_total_milliseconds(this) result(ms)
      implicit none
      class(t_timedelta), intent(in) :: this
      integer(kind=8) :: ms
      ms = this%ms_count
   end function timedelta_total_milliseconds

   !> @brief Convert a TimeDelta to a string
   !> @param this TimeDelta object
   !> @return String representation of the TimeDelta
   function timedelta_to_string(this) result(str)
      implicit none
      class(t_timedelta), intent(in) :: this
      character(len=TIMEDELTA_STRING_BUFFER_SIZE) :: str
      character(kind=c_char), dimension(TIMEDELTA_STRING_BUFFER_SIZE + 1) :: c_str

      call f_timedelta_to_string(this%ms_count, c_str, TIMEDELTA_STRING_BUFFER_SIZE + 1)
      call c_f_string(c_str, str)
   end function timedelta_to_string

   !> @brief Add two TimeDeltas
   !> @param ts1 First TimeDelta
   !> @param ts2 Second TimeDelta
   !> @return Sum of the two TimeDeltas
   pure function timedelta_add_timedelta(ts1, ts2) result(sum_ts)
      implicit none
      type(t_timedelta), intent(in) :: ts1, ts2
      type(t_timedelta) :: sum_ts

      sum_ts%ms_count = f_timedelta_add(ts1%ms_count, ts2%ms_count)
   end function timedelta_add_timedelta

   !> @brief Subtract one TimeDelta from another
   !> @param ts1 First TimeDelta
   !> @param ts2 Second TimeDelta
   !> @return Difference of the two TimeDeltas
   pure function timedelta_subtract_timedelta(ts1, ts2) result(diff_ts)
      implicit none
      type(t_timedelta), intent(in) :: ts1, ts2
      type(t_timedelta) :: diff_ts

      diff_ts%ms_count = f_timedelta_subtract(ts1%ms_count, ts2%ms_count)
   end function timedelta_subtract_timedelta

   !> @brief Multiply a TimeDelta by a factor
   !> @param ts TimeDelta to multiply
   !> @param factor Multiplication factor
   !> @return Multiplied TimeDelta
   pure function timedelta_multiply(ts, factor) result(product_ts)
      implicit none
      type(t_timedelta), intent(in) :: ts
      integer, intent(in) :: factor
      type(t_timedelta) :: product_ts

      product_ts%ms_count = f_timedelta_multiply(ts%ms_count, factor)
   end function timedelta_multiply

   !> @brief Divide a TimeDelta by a divisor
   !> @param ts TimeDelta to divide
   !> @param divisor Division factor
   !> @return Divided TimeDelta
   pure function timedelta_divide(ts, divisor) result(quotient_ts)
      implicit none
      type(t_timedelta), intent(in) :: ts
      integer, intent(in) :: divisor
      type(t_timedelta) :: quotient_ts

      quotient_ts%ms_count = f_timedelta_divide(ts%ms_count, divisor)
   end function timedelta_divide

   !> @brief Check if two TimeDeltas are equal
   !> @param ts1 First TimeDelta
   !> @param ts2 Second TimeDelta
   !> @return True if equal, False otherwise
   pure function timedelta_equals(ts1, ts2) result(res)
      implicit none
      type(t_timedelta), intent(in) :: ts1, ts2
      logical(c_bool) :: res_t
      logical :: res

      res_t = f_timedelta_equals(ts1%ms_count, ts2%ms_count)
      res = logical(res_t, 4)
   end function timedelta_equals

   !> @brief Check if two TimeDeltas are not equal
   !> @param ts1 First TimeDelta
   !> @param ts2 Second TimeDelta
   !> @return True if not equal, False otherwise
   pure function timedelta_not_equals(ts1, ts2) result(res)
      implicit none
      type(t_timedelta), intent(in) :: ts1, ts2
      logical(c_bool) :: res_t
      logical :: res

      res_t = f_timedelta_equals(ts1%ms_count, ts2%ms_count)
      res = .not. logical(res_t, 4)
   end function timedelta_not_equals

   !> @brief Check if one TimeDelta is less than another
   !> @param ts1 First TimeDelta
   !> @param ts2 Second TimeDelta
   !> @return True if ts1 < ts2, False otherwise
   pure function timedelta_less_than(ts1, ts2) result(res)
      implicit none
      type(t_timedelta), intent(in) :: ts1, ts2
      logical(c_bool) :: res_t
      logical :: res

      res_t = f_timedelta_less_than(ts1%ms_count, ts2%ms_count)
      res = logical(res_t, 4)
   end function timedelta_less_than

   !> @brief Check if one TimeDelta is greater than another
   !> @param ts1 First TimeDelta
   !> @param ts2 Second TimeDelta
   !> @return True if ts1 > ts2, False otherwise
   pure function timedelta_greater_than(ts1, ts2) result(res)
      implicit none
      type(t_timedelta), intent(in) :: ts1, ts2
      logical(c_bool) :: res_t
      logical :: res

      res_t = f_timedelta_greater_than(ts1%ms_count, ts2%ms_count)
      res = logical(res_t, 4)
   end function timedelta_greater_than

   !> @brief Check if one TimeDelta is less than or equal to another
   !> @param ts1 First TimeDelta
   !> @param ts2 Second TimeDelta
   !> @return True if ts1 <= ts2, False otherwise
   pure function timedelta_less_equal(ts1, ts2) result(res)
      implicit none
      type(t_timedelta), intent(in) :: ts1, ts2
      logical(c_bool) :: res_t
      logical :: res

      res_t = f_timedelta_less_equal(ts1%ms_count, ts2%ms_count)
      res = logical(res_t, 4)
   end function timedelta_less_equal

   !> @brief Check if one TimeDelta is greater than or equal to another
   !> @param ts1 First TimeDelta
   !> @param ts2 Second TimeDelta
   !> @return True if ts1 >= ts2, False otherwise
   pure function timedelta_greater_equal(ts1, ts2) result(res)
      implicit none
      type(t_timedelta), intent(in) :: ts1, ts2
      logical(c_bool) :: res_t
      logical :: res

      res_t = f_timedelta_greater_equal(ts1%ms_count, ts2%ms_count)
      res = logical(res_t, 4)
   end function timedelta_greater_equal

   !===========================================================================
   ! DateTime implementations
   !===========================================================================

   !> @brief Create a default DateTime (epoch)
   !> @return DateTime at the epoch
   pure function datetime_default() result(dt)
      implicit none
      type(t_datetime) :: dt
      dt%timestamp_ms = 0_c_int64_t
   end function datetime_default

   !> @brief Create a DateTime from year, month, and day
   !> @param year Year
   !> @param month Month (1-12)
   !> @param day Day (1-31)
   !> @return DateTime with the specified date
   pure function datetime_ymd(year, month, day) result(dt)
      implicit none
      integer, intent(in) :: year, month, day
      type(t_datetime) :: dt

      dt%timestamp_ms = f_datetime_create(year, month, day, 0, 0, 0, 0)
   end function datetime_ymd

   !> @brief Create a DateTime from year, month, day, hour, minute, and second
   !> @param year Year
   !> @param month Month (1-12)
   !> @param day Day (1-31)
   !> @param hour Hour (0-23)
   !> @param minute Minute (0-59)
   !> @param second Second (0-59)
   !> @return DateTime with the specified date and time
   pure function datetime_ymd_hms(year, month, day, hour, minute, second) result(dt)
      implicit none
      integer, intent(in) :: year, month, day, hour, minute, second
      type(t_datetime) :: dt

      dt%timestamp_ms = f_datetime_create(year, month, day, hour, minute, second, 0)
   end function datetime_ymd_hms

   !> @brief Create a DateTime from all components
   !> @param year Year
   !> @param month Month (1-12)
   !> @param day Day (1-31)
   !> @param hour Hour (0-23)
   !> @param minute Minute (0-59)
   !> @param second Second (0-59)
   !> @param millisecond Millisecond (0-999)
   !> @return DateTime with the specified date and time
   pure function datetime_complete(year, month, day, hour, minute, second, millisecond) result(dt)
      implicit none
      integer, intent(in) :: year, month, day, hour, minute, second, millisecond
      type(t_datetime) :: dt

      dt%timestamp_ms = f_datetime_create(year, month, day, hour, minute, second, millisecond)
   end function datetime_complete

   !> @brief Create a DateTime from a timestamp
   !> @param timestamp Milliseconds since the epoch
   !> @return DateTime corresponding to the timestamp
   pure function datetime_from_timestamp(timestamp) result(dt)
      implicit none
      integer(kind=8), intent(in) :: timestamp
      type(t_datetime) :: dt

      dt%timestamp_ms = timestamp
   end function datetime_from_timestamp

   !> @brief Parse a DateTime from a string
   !> @param str String representation of a DateTime
   !> @param format_str Format specification (similar to strftime)
   !> @return DateTime parsed from the string
   function datetime_strptime(str, format_str) result(dt)
      implicit none
      character(len=*), intent(in) :: str, format_str
      type(t_datetime) :: dt

      ! Convert Fortran strings to C-compatible character arrays
      character(kind=c_char, len=1) :: c_str(len_trim(str) + 1), c_format(len_trim(format_str) + 1)
      integer :: i
      logical :: has_milliseconds

      ! If the 4th character from the end of the format string is a dot, it indicates milliseconds
      has_milliseconds = (len_trim(format_str) >= 4 .and. format_str(len_trim(format_str) - 3:len_trim(format_str) - 3) == '.')

      do i = 1, len_trim(str)
         c_str(i) = str(i:i)
      end do
      c_str(len_trim(str) + 1) = c_null_char

      do i = 1, len_trim(format_str)
         c_format(i) = format_str(i:i)
      end do
      c_format(len_trim(format_str) + 1) = c_null_char

      dt%timestamp_ms = f_datetime_strptime(c_str, c_format, len_trim(str), len_trim(format_str))
   end function datetime_strptime

   !> @brief Parse a DateTime from a string with automatic format detection
   !> @param str String representation of a DateTime
   !> @return DateTime parsed from the string
   function datetime_strptime_auto(str) result(dt)
      implicit none
      character(len=*), intent(in) :: str
      type(t_datetime) :: dt
      dt = datetime_strptime(str, "auto")
   end function datetime_strptime_auto

   !> @brief Get the current DateTime
   !> @return DateTime representing the current time
   function now() result(dt)
      implicit none
      type(t_datetime) :: dt

      dt%timestamp_ms = f_datetime_now()
   end function now

   !> @brief Parse a DateTime from a string using an array of format options
   !> @param str String representation of a DateTime
   !> @param formats Array of format strings to try in order
   !> @return DateTime parsed from the string using the first successful format
   function datetime_strptime_with_formats(str, formats) result(dt)
      implicit none
      character(len=*), intent(in) :: str
      character(len=*), intent(in) :: formats(:)
      type(t_datetime) :: dt

      ! Variables for C interface
      character(kind=c_char, len=1) :: c_str(len_trim(str) + 1)
      type(c_ptr), allocatable, target :: formats_ptr_array(:)
      type(c_ptr) :: formats_array_ptr
      integer(c_int), allocatable :: formats_len(:)
      integer :: i, j, max_len_safe
      character(kind=c_char), allocatable, target :: format_chars(:, :)

      ! Input validation
      if (len_trim(str) <= 0 .or. size(formats) <= 0) then
         dt = t_datetime()
         return
      end if

      ! Prevent excessive format counts and string lengths for security
      if (size(formats) > 1000 .or. len_trim(str) > 10000) then
         dt = t_datetime()
         return
      end if

      ! Calculate maximum format length safely (limit individual formats to 1000 chars)
      max_len_safe = 0
      do i = 1, size(formats)
         max_len_safe = max(max_len_safe, min(len_trim(formats(i)), 1000))
      end do

      ! Convert input string to C format
      do i = 1, len_trim(str)
         c_str(i) = str(i:i)
      end do
      c_str(len_trim(str) + 1) = c_null_char

      ! Allocate arrays with safe maximum length
      allocate (formats_ptr_array(size(formats)))
      allocate (formats_len(size(formats)))
      allocate (format_chars(max_len_safe + 1, size(formats)))

      ! Create C-compatible format strings and their pointers
      do i = 1, size(formats)
         formats_len(i) = min(len_trim(formats(i)), 1000) ! Limit format length
         if (formats_len(i) > 0) then
            do j = 1, formats_len(i)
               format_chars(j, i) = formats(i) (j:j)
            end do
            format_chars(formats_len(i) + 1, i) = c_null_char
            formats_ptr_array(i) = c_loc(format_chars(1, i))
         else
            formats_ptr_array(i) = c_null_ptr
         end if
      end do

      ! Get pointer to the array of pointers
      formats_array_ptr = c_loc(formats_ptr_array(1))

      ! Call the C function
      dt%timestamp_ms = f_datetime_strptime_with_formats(c_str, formats_array_ptr, &
                                                         len_trim(str), formats_len, size(formats))

      ! Clean up allocated memory
      deallocate (formats_ptr_array)
      deallocate (formats_len)
      deallocate (format_chars)
   end function datetime_strptime_with_formats

   !> @brief Parse a DateTime with automatic detection and fallback formats
   !> @param str String representation of a DateTime
   !> @param user_formats Optional array of user format strings to try
   !> @return DateTime parsed from the string using auto-detection or fallback formats
   function datetime_strptime_auto_with_fallback(str, user_formats) result(dt)
      implicit none
      character(len=*), intent(in) :: str
      character(len=*), intent(in), optional :: user_formats(:)
      type(t_datetime) :: dt

      ! Variables for C interface
      character(kind=c_char, len=1) :: c_str(len_trim(str) + 1)
      type(c_ptr), allocatable, target :: formats_ptr_array(:)
      type(c_ptr) :: formats_array_ptr
      integer(c_int), allocatable :: formats_len(:)
      integer :: i, j, num_formats, max_len_safe
      character(kind=c_char), allocatable, target :: format_chars(:, :)

      ! Input validation
      if (len_trim(str) <= 0) then
         dt = t_datetime()
         return
      end if

      ! Prevent excessive string length for security
      if (len_trim(str) > DATETIME_MAX_STRING_LENGTH) then
         dt = t_datetime()
         return
      end if

      ! Convert input string to C format
      do i = 1, len_trim(str)
         c_str(i) = str(i:i)
      end do
      c_str(len_trim(str) + 1) = c_null_char

      if (present(user_formats)) then
         num_formats = size(user_formats)

         ! Prevent excessive format counts for security
         if (num_formats > DATETIME_MAX_FORMAT_COUNT) then
            dt = t_datetime()
            return
         end if

         ! Calculate maximum format length safely (limit individual formats to 1000 chars)
         max_len_safe = 0
         do i = 1, num_formats
            max_len_safe = max(max_len_safe, min(len_trim(user_formats(i)), 1000))
         end do

         allocate (formats_ptr_array(num_formats))
         allocate (formats_len(num_formats))
         allocate (format_chars(max_len_safe + 1, num_formats))

         ! Create C-compatible format strings and their pointers
         do i = 1, num_formats
            formats_len(i) = min(len_trim(user_formats(i)), 1000) ! Limit format length
            if (formats_len(i) > 0) then
               do j = 1, formats_len(i)
                  format_chars(j, i) = user_formats(i) (j:j)
               end do
               format_chars(formats_len(i) + 1, i) = c_null_char
               formats_ptr_array(i) = c_loc(format_chars(1, i))
            else
               formats_ptr_array(i) = c_null_ptr
            end if
         end do

         ! Get pointer to the array of pointers
         formats_array_ptr = c_loc(formats_ptr_array(1))

         ! Call the C function with fallback formats
         dt%timestamp_ms = f_datetime_strptime_auto_with_fallback(c_str, formats_array_ptr, &
                                                                  len_trim(str), formats_len, num_formats)

         ! Clean up allocated memory
         deallocate (formats_ptr_array)
         deallocate (formats_len)
         deallocate (format_chars)
      else
         ! Call with null fallback formats (just auto-detection)
         allocate (formats_len(1))
         dt%timestamp_ms = f_datetime_strptime_auto_with_fallback(c_str, c_null_ptr, &
                                                                  len_trim(str), formats_len, 0)
         deallocate (formats_len)
      end if
   end function datetime_strptime_auto_with_fallback

   !> @brief Get the year component from a DateTime
   !> @param this DateTime object
   !> @return Year component
   pure function datetime_year(this) result(year)
      implicit none
      class(t_datetime), intent(in) :: this
      integer :: year

      year = f_datetime_get_year(this%timestamp_ms)
   end function datetime_year

   !> @brief Get the month component from a DateTime
   !> @param this DateTime object
   !> @return Month component (1-12)
   pure function datetime_month(this) result(month)
      implicit none
      class(t_datetime), intent(in) :: this
      integer :: month

      month = f_datetime_get_month(this%timestamp_ms)
   end function datetime_month

   !> @brief Get the day component from a DateTime
   !> @param this DateTime object
   !> @return Day component (1-31)
   pure function datetime_day(this) result(day)
      implicit none
      class(t_datetime), intent(in) :: this
      integer :: day

      day = f_datetime_get_day(this%timestamp_ms)
   end function datetime_day

   !> @brief Get the hour component from a DateTime
   !> @param this DateTime object
   !> @return Hour component (0-23)
   pure function datetime_hour(this) result(hour)
      implicit none
      class(t_datetime), intent(in) :: this
      integer :: hour

      hour = f_datetime_get_hour(this%timestamp_ms)
   end function datetime_hour

   !> @brief Get the minute component from a DateTime
   !> @param this DateTime object
   !> @return Minute component (0-59)
   pure function datetime_minute(this) result(minute)
      implicit none
      class(t_datetime), intent(in) :: this
      integer :: minute

      minute = f_datetime_get_minute(this%timestamp_ms)
   end function datetime_minute

   !> @brief Get the second component from a DateTime
   !> @param this DateTime object
   !> @return Second component (0-59)
   pure function datetime_second(this) result(second)
      implicit none
      class(t_datetime), intent(in) :: this
      integer :: second

      second = f_datetime_get_second(this%timestamp_ms)
   end function datetime_second

   !> @brief Get the millisecond component from a DateTime
   !> @param this DateTime object
   !> @return Millisecond component (0-999)
   pure function datetime_millisecond(this) result(ms)
      implicit none
      class(t_datetime), intent(in) :: this
      integer :: ms

      ms = f_datetime_get_millisecond(this%timestamp_ms)
   end function datetime_millisecond

   !> @brief Get the Julian Day Number from a DateTime
   !> @param this DateTime object
   !> @return Julian Day Number (integer days since JD epoch)
   pure function datetime_julian_day_number(this) result(jdn)
      implicit none
      class(t_datetime), intent(in) :: this
      integer(kind=8) :: jdn

      jdn = f_datetime_get_julian_day_number(this%timestamp_ms)
   end function datetime_julian_day_number

   !> @brief Get the Julian Day (with fractional part) from a DateTime
   !> @param this DateTime object
   !> @return Julian Day with fractional part (days.fraction since JD epoch)
   pure function datetime_julian_day(this) result(jd)
      implicit none
      class(t_datetime), intent(in) :: this
      real(kind=8) :: jd

      jd = f_datetime_get_julian_day(this%timestamp_ms)
   end function datetime_julian_day

   !> @brief Get the Julian Century from a DateTime
   !> @param this DateTime object
   !> @return Julian Century (time unit used in astronomy)
   pure function datetime_julian_century(this) result(jc)
      implicit none
      class(t_datetime), intent(in) :: this
      real(kind=8) :: jc

      jc = f_datetime_get_julian_century(this%timestamp_ms)
   end function datetime_julian_century

   !> @brief Get the timestamp from a DateTime
   !> @param this DateTime object
   !> @return Timestamp (milliseconds since epoch)
   pure function datetime_timestamp(this) result(timestamp)
      implicit none
      class(t_datetime), intent(in) :: this
      integer(kind=8) :: timestamp

      timestamp = this%timestamp_ms
   end function datetime_timestamp

   !> @brief Format a DateTime to a string
   !> @param this DateTime object
   !> @param format Format specification (similar to strftime)
   !> @return String representation of the DateTime
   function datetime_strftime(this, date_format, show_milliseconds) result(str)
      implicit none
      class(t_datetime), intent(in) :: this
      character(len=*), intent(in) :: date_format
      logical, intent(in), optional :: show_milliseconds
      character(len=DATETIME_STRING_BUFFER_SIZE) :: str

      ! Convert Fortran string to C-compatible character array
      character(kind=c_char, len=1) :: c_format(len_trim(date_format) + 1), c_str(65)
      integer :: i
      logical :: show_milliseconds_l

      do i = 1, len_trim(date_format)
         c_format(i) = date_format(i:i)
      end do
      c_format(len_trim(date_format) + 1) = c_null_char

      if (present(show_milliseconds)) then
         show_milliseconds_l = show_milliseconds
      else
         show_milliseconds_l = .false.
      end if

      if (show_milliseconds_l) then
         call f_datetime_strftime_ms(this%timestamp_ms, c_format, c_str, &
                                     len_trim(date_format), DATETIME_STRING_BUFFER_SIZE + 1)
      else
         call f_datetime_strftime(this%timestamp_ms, c_format, c_str, &
                                  len_trim(date_format), DATETIME_STRING_BUFFER_SIZE + 1)
      end if

      call c_f_string(c_str, str)

   end function datetime_strftime

   !> @brief Convert a DateTime to an ISO 8601 string
   !> @param this DateTime object
   !> @param msec Flag to include milliseconds in the output
   !> @return ISO 8601 string representation of the DateTime
   function datetime_to_iso_string(this, msec) result(str)
      implicit none
      class(t_datetime), intent(in) :: this
      character(len=DATETIME_STRING_BUFFER_SIZE) :: str
      logical, intent(in), optional :: msec

      character(kind=c_char, len=1) :: c_str(DATETIME_STRING_BUFFER_SIZE + 1)
      logical :: msec_l

      if (present(msec)) then
         msec_l = msec
      else
         msec_l = .false.
      end if

      if (msec_l) then
         str = datetime_strftime(this, "%Y-%m-%dT%H:%M:%S", .true.)
      else
         str = datetime_strftime(this, "%Y-%m-%dT%H:%M:%S", .false.)
      end if

   end function datetime_to_iso_string

   !> @brief Add a TimeDelta to a DateTime
   !> @param dt DateTime
   !> @param ts TimeDelta
   !> @return DateTime resulting from the addition
   pure function datetime_add_timedelta(dt, ts) result(result_dt)
      implicit none
      type(t_datetime), intent(in) :: dt
      type(t_timedelta), intent(in) :: ts
      type(t_datetime) :: result_dt

      result_dt%timestamp_ms = f_datetime_add_timedelta(dt%timestamp_ms, ts%ms_count)
   end function datetime_add_timedelta

   !> @brief Subtract a TimeDelta from a DateTime
   !> @param dt DateTime
   !> @param ts TimeDelta
   !> @return DateTime resulting from the subtraction
   pure function datetime_subtract_timedelta(dt, ts) result(result_dt)
      implicit none
      type(t_datetime), intent(in) :: dt
      type(t_timedelta), intent(in) :: ts
      type(t_datetime) :: result_dt

      result_dt%timestamp_ms = f_datetime_subtract_timedelta(dt%timestamp_ms, ts%ms_count)
   end function datetime_subtract_timedelta

   !> @brief Calculate the difference between two DateTimes
   !> @param dt1 First DateTime
   !> @param dt2 Second DateTime
   !> @return TimeDelta representing the difference
   pure function datetime_difference(dt1, dt2) result(ts)
      implicit none
      type(t_datetime), intent(in) :: dt1, dt2
      type(t_timedelta) :: ts

      ts%ms_count = f_datetime_difference(dt1%timestamp_ms, dt2%timestamp_ms)
   end function datetime_difference

   !> @brief Check if two DateTimes are equal
   !> @param dt1 First DateTime
   !> @param dt2 Second DateTime
   !> @return True if equal, False otherwise
   pure function datetime_equals(dt1, dt2) result(res)
      implicit none
      type(t_datetime), intent(in) :: dt1, dt2
      logical(c_bool) :: res_t
      logical :: res

      res_t = f_datetime_equals(dt1%timestamp_ms, dt2%timestamp_ms)
      res = logical(res_t, 4)
   end function datetime_equals

   !> @brief Check if two DateTimes are not equal
   !> @param dt1 First DateTime
   !> @param dt2 Second DateTime
   !> @return True if not equal, False otherwise
   pure function datetime_not_equals(dt1, dt2) result(res)
      implicit none
      type(t_datetime), intent(in) :: dt1, dt2
      logical(c_bool) :: res_t
      logical :: res

      res_t = f_datetime_equals(dt1%timestamp_ms, dt2%timestamp_ms)
      res = .not. logical(res_t, 4)
   end function datetime_not_equals

   !> @brief Check if one DateTime is less than another
   !> @param dt1 First DateTime
   !> @param dt2 Second DateTime
   !> @return True if dt1 < dt2, False otherwise
   pure function datetime_less_than(dt1, dt2) result(res)
      implicit none
      type(t_datetime), intent(in) :: dt1, dt2
      logical(c_bool) :: res_t
      logical :: res

      res_t = f_datetime_less_than(dt1%timestamp_ms, dt2%timestamp_ms)
      res = logical(res_t, 4)
   end function datetime_less_than

   !> @brief Check if one DateTime is greater than another
   !> @param dt1 First DateTime
   !> @param dt2 Second DateTime
   !> @return True if dt1 > dt2, False otherwise
   pure function datetime_greater_than(dt1, dt2) result(res)
      implicit none
      type(t_datetime), intent(in) :: dt1, dt2
      logical(c_bool) :: res_t
      logical :: res

      res_t = f_datetime_greater_than(dt1%timestamp_ms, dt2%timestamp_ms)
      res = logical(res_t, 4)
   end function datetime_greater_than

   !> @brief Check if one DateTime is less than or equal to another
   !> @param dt1 First DateTime
   !> @param dt2 Second DateTime
   !> @return True if dt1 <= dt2, False otherwise
   pure function datetime_less_equal(dt1, dt2) result(res)
      implicit none
      type(t_datetime), intent(in) :: dt1, dt2
      logical(c_bool) :: res_t
      logical :: res

      res_t = f_datetime_less_equal(dt1%timestamp_ms, dt2%timestamp_ms)
      res = logical(res_t, 4)
   end function datetime_less_equal

   !> @brief Check if one DateTime is greater than or equal to another
   !> @param dt1 First DateTime
   !> @param dt2 Second DateTime
   !> @return True if dt1 >= dt2, False otherwise
   pure function datetime_greater_equal(dt1, dt2) result(res)
      implicit none
      type(t_datetime), intent(in) :: dt1, dt2
      logical(c_bool) :: res_t
      logical :: res

      res_t = f_datetime_greater_equal(dt1%timestamp_ms, dt2%timestamp_ms)
      res = logical(res_t, 4)
   end function datetime_greater_equal

   !> @brief Check if a given datetime object is valid
   !> @param dt DateTime object
   !> @return True if valid, False otherwise
   pure function datetime_is_valid(dt) result(is_valid)
      implicit none
      class(t_datetime), intent(in) :: dt
      logical :: is_valid
      logical(c_bool) :: is_valid_c

      is_valid_c = f_datetime_is_valid(dt%timestamp_ms)
      is_valid = logical(is_valid_c, 4)
   end function datetime_is_valid

   !> @brief Helper function to convert C string to Fortran string
   !> @param c_string C string
   !> @param f_string Fortran string
   subroutine c_f_string(c_string, f_string)
      implicit none
      character(kind=c_char), dimension(DATETIME_STRING_BUFFER_SIZE + 1), intent(in) :: c_string
      character(len=DATETIME_STRING_BUFFER_SIZE), intent(out) :: f_string
      integer :: i

      f_string = ""
      do i = 1, len(f_string)
         if (c_string(i) == c_null_char) exit
         f_string(i:i) = c_string(i)
      end do
   end subroutine c_f_string

   ! Create a null datetime object
   pure function null_datetime() result(dt)
      type(t_datetime) :: dt
      dt%timestamp_ms = f_datetime_invalid_timestamp()
   end function null_datetime

end module mod_datetime
