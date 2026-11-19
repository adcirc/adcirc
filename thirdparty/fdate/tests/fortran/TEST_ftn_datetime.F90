!> @file test_datetime_mod.f90
!> @brief Custom test suite for datetime_mod module without external dependencies
!>
!> Comprehensive test suite for both TimeDelta and DateTime functionality

module test_utils
   implicit none

   ! Variables to track test results
   integer, save :: tests_run = 0
   integer, save :: tests_failed = 0

   ! Interface for assertion subroutines
   interface assert_equal
      module procedure assert_equal_int
      module procedure assert_equal_int8
      module procedure assert_equal_string
      module procedure assert_equal_real8
   end interface assert_equal

   public

contains

   ! Assertion utilities
   subroutine assert_equal_int(expected, actual, message)
      integer, intent(in) :: expected, actual
      character(len=*), intent(in) :: message

      tests_run = tests_run + 1
      if (expected /= actual) then
         tests_failed = tests_failed + 1
         write (*, '(3A,I0,A,I0)') "FAILED: ", message, &
            " (Expected: ", expected, ", Got: ", actual, ")"
      end if
   end subroutine assert_equal_int

   subroutine assert_equal_int8(expected, actual, message)
      integer(kind=8), intent(in) :: expected, actual
      character(len=*), intent(in) :: message

      tests_run = tests_run + 1
      if (expected /= actual) then
         tests_failed = tests_failed + 1
         write (*, '(3A,I0,A,I0)') "FAILED: ", message, &
            " (Expected: ", expected, ", Got: ", actual, ")"
      end if
   end subroutine assert_equal_int8

   subroutine assert_equal_string(expected, actual, message)
      character(len=*), intent(in) :: expected, actual, message

      tests_run = tests_run + 1
      if (expected /= actual) then
         tests_failed = tests_failed + 1
         write (*, '(7A)') "FAILED: ", message, &
            " (Expected: '", trim(expected), "', Got: '", trim(actual), "')"
      end if
   end subroutine assert_equal_string

   subroutine assert_equal_real8(expected, actual, message)
      real(kind=8), intent(in) :: expected, actual
      character(len=*), intent(in) :: message
      real(kind=8) :: tolerance

      tests_run = tests_run + 1
      tolerance = 1.0d-10
      if (abs(expected - actual) > tolerance) then
         tests_failed = tests_failed + 1
         write (*, '(3A,F0.12,A,F0.12)') "FAILED: ", message, &
            " (Expected: ", expected, ", Got: ", actual, ")"
      end if
   end subroutine assert_equal_real8

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      tests_run = tests_run + 1
      if (.not. condition) then
         tests_failed = tests_failed + 1
         write (*, '(A,A)') "FAILED: ", message
      end if
   end subroutine assert_true

   subroutine assert_false(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      tests_run = tests_run + 1
      if (condition) then
         tests_failed = tests_failed + 1
         write (*, '(A,A)') "FAILED: ", message
      end if
   end subroutine assert_false

   ! Test case runner
   subroutine run_test(test_procedure, test_name)
      interface
         subroutine test_procedure()
            implicit none
         end subroutine test_procedure
      end interface
      character(len=*), intent(in) :: test_name

      integer :: prev_failed, prev_run

      prev_failed = tests_failed
      prev_run = tests_run

      write (*, '(A)') "Running test: "//trim(test_name)
      call test_procedure()

      if (tests_failed == prev_failed) then
         write (*, '(A,I0,A)') "  PASSED (", tests_run - prev_run, " assertions)"
      else
         write (*, '(A,I0,A,I0,A)') "  FAILED (", tests_failed - prev_failed, &
            " of ", tests_run - prev_run, " assertions failed)"
      end if
      write (*, '(A)') ""
   end subroutine run_test

   subroutine print_test_summary()
      write (*, '(A)') "----------------------------------------"
      write (*, '(A)') "Test Summary"
      write (*, '(A)') "----------------------------------------"
      write (*, '(A,I0)') "Total tests run: ", tests_run
      write (*, '(A,I0)') "Tests passed:    ", tests_run - tests_failed
      write (*, '(A,I0)') "Tests failed:    ", tests_failed

      if (tests_failed == 0) then
         write (*, '(A)') "ALL TESTS PASSED!"
      else
         write (*, '(A)') "SOME TESTS FAILED!"
      end if
   end subroutine print_test_summary

end module test_utils

module timedelta_tests
   use test_utils, only: run_test, print_test_summary, &
                         assert_equal, assert_true, assert_false, &
                         tests_run, tests_failed
   use mod_datetime, only: t_timedelta, operator(+), operator(-), &
                           operator(*), operator(/), operator(==), &
                           operator(/=), operator(<), operator(>), &
                           operator(<=), operator(>=)
   implicit none

   public

contains

   subroutine test_timedelta_default()
      use, intrinsic :: iso_fortran_env, only: int64
      type(t_timedelta) :: ts_zero

      ts_zero = t_timedelta()

      call assert_equal(0_int64, ts_zero%total_milliseconds(), "Default TimeDelta should have 0 ms")
      call assert_equal(0, ts_zero%days(), "Default TimeDelta should have 0 days")
      call assert_equal(0, ts_zero%hours(), "Default TimeDelta should have 0 hours")
      call assert_equal(0, ts_zero%minutes(), "Default TimeDelta should have 0 minutes")
      call assert_equal(0, ts_zero%seconds(), "Default TimeDelta should have 0 seconds")
      call assert_equal(0, ts_zero%milliseconds(), "Default TimeDelta should have 0 milliseconds")
   end subroutine test_timedelta_default

   subroutine test_timedelta_components()
      type(t_timedelta) :: ts_components

      ts_components = t_timedelta(2, 3, 4, 5, 6)

      call assert_equal(2, ts_components%days(), "Component days")
      call assert_equal(3, ts_components%hours(), "Component hours")
      call assert_equal(4, ts_components%minutes(), "Component minutes")
      call assert_equal(5, ts_components%seconds(), "Component seconds")
      call assert_equal(6, ts_components%milliseconds(), "Component milliseconds")
   end subroutine test_timedelta_components

   subroutine test_timedelta_total_accessors()
      use, intrinsic :: iso_fortran_env, only: int64
      type(t_timedelta) :: ts_components
      integer(kind=8) :: expected_ms

      ts_components = t_timedelta(2, 3, 4, 5, 6)

      expected_ms = (2_int64*24*60*60*1000) + & ! 2 days
                    (3_int64*60*60*1000) + & ! 3 hours
                    (4_int64*60*1000) + & ! 4 minutes
                    (5_int64*1000) + & ! 5 seconds
                    6_int64 ! 6 milliseconds

      call assert_equal(expected_ms, ts_components%total_milliseconds(), "Total milliseconds")
      call assert_equal(expected_ms/1000, ts_components%total_seconds(), "Total seconds")
      call assert_equal(expected_ms/(60*1000), ts_components%total_minutes(), "Total minutes")
      call assert_equal(expected_ms/(60*60*1000), ts_components%total_hours(), "Total hours")
      call assert_equal(expected_ms/(24*60*60*1000), ts_components%total_days(), "Total days")
   end subroutine test_timedelta_total_accessors

   subroutine test_timedelta_addition()
      type(t_timedelta) :: ts1, ts2, result

      ts1 = t_timedelta(1, 2, 3, 4, 5)
      ts2 = t_timedelta(2, 3, 4, 5, 6)

      result = ts1 + ts2

      call assert_equal(3, result%days(), "Addition - days")
      call assert_equal(5, result%hours(), "Addition - hours")
      call assert_equal(7, result%minutes(), "Addition - minutes")
      call assert_equal(9, result%seconds(), "Addition - seconds")
      call assert_equal(11, result%milliseconds(), "Addition - milliseconds")
   end subroutine test_timedelta_addition

   subroutine test_timedelta_subtraction()
      type(t_timedelta) :: ts1, ts2, result

      ts1 = t_timedelta(3, 4, 5, 6, 7)
      ts2 = t_timedelta(1, 2, 3, 4, 5)

      result = ts1 - ts2

      call assert_equal(2, result%days(), "Subtraction - days")
      call assert_equal(2, result%hours(), "Subtraction - hours")
      call assert_equal(2, result%minutes(), "Subtraction - minutes")
      call assert_equal(2, result%seconds(), "Subtraction - seconds")
      call assert_equal(2, result%milliseconds(), "Subtraction - milliseconds")
   end subroutine test_timedelta_subtraction

   subroutine test_timedelta_multiplication()
      type(t_timedelta) :: ts, result

      ts = t_timedelta(1, 2, 3, 4, 5)
      result = ts*2

      call assert_equal(2, result%days(), "Multiplication - days")
      call assert_equal(4, result%hours(), "Multiplication - hours")
      call assert_equal(6, result%minutes(), "Multiplication - minutes")
      call assert_equal(8, result%seconds(), "Multiplication - seconds")
      call assert_equal(10, result%milliseconds(), "Multiplication - milliseconds")
   end subroutine test_timedelta_multiplication

   subroutine test_timedelta_division()
      type(t_timedelta) :: ts, result

      ts = t_timedelta(2, 4, 6, 8, 10)
      result = ts/2

      call assert_equal(1, result%days(), "Division - days")
      call assert_equal(2, result%hours(), "Division - hours")
      call assert_equal(3, result%minutes(), "Division - minutes")
      call assert_equal(4, result%seconds(), "Division - seconds")
      call assert_equal(5, result%milliseconds(), "Division - milliseconds")
   end subroutine test_timedelta_division

   subroutine test_timedelta_comparisons()
      type(t_timedelta) :: ts1, ts2, ts3

      ts1 = t_timedelta(1, 0, 0, 0, 0)
      ts2 = t_timedelta(2, 0, 0, 0, 0)
      ts3 = t_timedelta(1, 0, 0, 0, 0)

      ! Equality
      call assert_true(ts1 == ts3, "Equality operator for equal timedeltas")
      call assert_false(ts1 == ts2, "Equality operator for different timedeltas")

      ! Inequality
      call assert_true(ts1 /= ts2, "Inequality operator for different timedeltas")
      call assert_false(ts1 /= ts3, "Inequality operator for equal timedeltas")

      ! Less than
      call assert_true(ts1 < ts2, "Less than operator when less")
      call assert_false(ts2 < ts1, "Less than operator when greater")
      call assert_false(ts1 < ts3, "Less than operator when equal")

      ! Greater than
      call assert_true(ts2 > ts1, "Greater than operator when greater")
      call assert_false(ts1 > ts2, "Greater than operator when less")
      call assert_false(ts1 > ts3, "Greater than operator when equal")

      ! Less than or equal
      call assert_true(ts1 <= ts2, "Less than or equal operator when less")
      call assert_true(ts1 <= ts3, "Less than or equal operator when equal")
      call assert_false(ts2 <= ts1, "Less than or equal operator when greater")

      ! Greater than or equal
      call assert_true(ts2 >= ts1, "Greater than or equal operator when greater")
      call assert_true(ts1 >= ts3, "Greater than or equal operator when equal")
      call assert_false(ts1 >= ts2, "Greater than or equal operator when less")
   end subroutine test_timedelta_comparisons

   subroutine test_timedelta_to_string()
      type(t_timedelta) :: ts1, ts2, ts3, ts4

      ts1 = t_timedelta(1, 2, 3, 4, 5)
      ts2 = t_timedelta(0, 2, 3, 4, 5)
      ts3 = t_timedelta(1, 2, 3, 4, 0)
      ts4 = t_timedelta(0, 2, 3, 4, 0)

      call assert_equal("1d 02:03:04.005", ts1%to_string(), "toString with days and ms")
      call assert_equal("02:03:04.005", ts2%to_string(), "toString without days")
      call assert_equal("1d 02:03:04", ts3%to_string(), "toString with days, no ms")
      call assert_equal("02:03:04", ts4%to_string(), "toString without days or ms")
   end subroutine test_timedelta_to_string

   subroutine test_timedelta_negative_duration()
      use, intrinsic :: iso_fortran_env, only: int64
      type(t_timedelta) :: ts_pos, ts_neg

      ts_pos = t_timedelta(1, 2, 3, 4, 5)
      ts_neg = ts_pos*(-1)

      call assert_equal(-1, ts_neg%days(), "Negative timedelta - days")
      call assert_equal(-2, ts_neg%hours(), "Negative timedelta - hours")
      call assert_equal(-3, ts_neg%minutes(), "Negative timedelta - minutes")
      call assert_equal(-4, ts_neg%seconds(), "Negative timedelta - seconds")
      call assert_equal(-5, ts_neg%milliseconds(), "Negative timedelta - milliseconds")

      call assert_equal(-1_int64*ts_pos%total_days(), ts_neg%total_days(), "Negative timedelta - total days")
      call assert_equal(-1_int64*ts_pos%total_hours(), ts_neg%total_hours(), "Negative timedelta - total hours")
      call assert_equal(-1_int64*ts_pos%total_minutes(), ts_neg%total_minutes(), "Negative timedelta - total minutes")
      call assert_equal(-1_int64*ts_pos%total_seconds(), ts_neg%total_seconds(), "Negative timedelta - total seconds")
      call assert_equal(-1_int64*ts_pos%total_milliseconds(), ts_neg%total_milliseconds(), "Negative timedelta - total ms")
   end subroutine test_timedelta_negative_duration

   subroutine test_timedelta_overflow_handling()
      use, intrinsic :: iso_fortran_env, only: int64
      use mod_datetime, only: t_timedelta
      implicit none
      type(t_timedelta) :: ts

      ! 25 hours should become 1 day 1 hour
      ts = t_timedelta(hours=25)

      call assert_equal(1, ts%days(), "Overflow 25 hours - days")
      call assert_equal(1, ts%hours(), "Overflow 25 hours - hours")
      call assert_equal(25_int64, ts%total_hours(), "Overflow 25 hours - total hours")

      ! 60 minutes should become 1 hour
      ts = t_timedelta(minutes=60)

      call assert_equal(0, ts%days(), "Overflow 60 minutes - days")
      call assert_equal(1, ts%hours(), "Overflow 60 minutes - hours")
      call assert_equal(0, ts%minutes(), "Overflow 60 minutes - minutes")
      call assert_equal(60_int64, ts%total_minutes(), "Overflow 60 minutes - total minutes")

      ! 60 seconds should become 1 minute
      ts = t_timedelta(seconds=60)

      call assert_equal(0, ts%days(), "Overflow 60 seconds - days")
      call assert_equal(0, ts%hours(), "Overflow 60 seconds - hours")
      call assert_equal(1, ts%minutes(), "Overflow 60 seconds - minutes")
      call assert_equal(0, ts%seconds(), "Overflow 60 seconds - seconds")
      call assert_equal(60_int64, ts%total_seconds(), "Overflow 60 seconds - total seconds")

      ! 1000 milliseconds should become 1 second
      ts = t_timedelta(milliseconds=1000)

      call assert_equal(0, ts%days(), "Overflow 1000 ms - days")
      call assert_equal(0, ts%hours(), "Overflow 1000 ms - hours")
      call assert_equal(0, ts%minutes(), "Overflow 1000 ms - minutes")
      call assert_equal(1, ts%seconds(), "Overflow 1000 ms - seconds")
      call assert_equal(0, ts%milliseconds(), "Overflow 1000 ms - milliseconds")
      call assert_equal(1000_int64, ts%total_milliseconds(), "Overflow 1000 ms - total milliseconds")
   end subroutine test_timedelta_overflow_handling

end module timedelta_tests

module datetime_tests
   use test_utils, only: run_test, print_test_summary, &
                         assert_equal, assert_true, assert_false, &
                         tests_run, tests_failed
   use mod_datetime, only: t_datetime, t_timedelta, now
   implicit none

   public

contains

   subroutine test_datetime_default()
      use, intrinsic :: iso_fortran_env, only: int64
      type(t_datetime) :: dt_epoch

      dt_epoch = t_datetime()

      call assert_equal(0_int64, dt_epoch%timestamp(), "Default DateTime timestamp")
      call assert_equal(1970, dt_epoch%year(), "Default DateTime year")
      call assert_equal(1, dt_epoch%month(), "Default DateTime month")
      call assert_equal(1, dt_epoch%day(), "Default DateTime day")
      call assert_equal(0, dt_epoch%hour(), "Default DateTime hour")
      call assert_equal(0, dt_epoch%minute(), "Default DateTime minute")
      call assert_equal(0, dt_epoch%second(), "Default DateTime second")
      call assert_equal(0, dt_epoch%millisecond(), "Default DateTime millisecond")
   end subroutine test_datetime_default

   subroutine test_datetime_ymd()
      type(t_datetime) :: dt_ymd

      dt_ymd = t_datetime(2022, 1, 31)

      call assert_equal(2022, dt_ymd%year(), "YMD DateTime year")
      call assert_equal(1, dt_ymd%month(), "YMD DateTime month")
      call assert_equal(31, dt_ymd%day(), "YMD DateTime day")
      call assert_equal(0, dt_ymd%hour(), "YMD DateTime hour")
      call assert_equal(0, dt_ymd%minute(), "YMD DateTime minute")
      call assert_equal(0, dt_ymd%second(), "YMD DateTime second")
      call assert_equal(0, dt_ymd%millisecond(), "YMD DateTime millisecond")
   end subroutine test_datetime_ymd

   subroutine test_datetime_complete()
      type(t_datetime) :: dt_complete

      dt_complete = t_datetime(2022, 1, 31, 12, 34, 56, 789)

      call assert_equal(2022, dt_complete%year(), "Complete DateTime year")
      call assert_equal(1, dt_complete%month(), "Complete DateTime month")
      call assert_equal(31, dt_complete%day(), "Complete DateTime day")
      call assert_equal(12, dt_complete%hour(), "Complete DateTime hour")
      call assert_equal(34, dt_complete%minute(), "Complete DateTime minute")
      call assert_equal(56, dt_complete%second(), "Complete DateTime second")
      call assert_equal(789, dt_complete%millisecond(), "Complete DateTime millisecond")
   end subroutine test_datetime_complete

   subroutine test_datetime_from_timestamp()
      type(t_datetime) :: dt1, dt2

      dt1 = t_datetime(2022, 3, 15, 14, 30, 45, 500)
      ! Create a new datetime from the timestamp of dt1
      dt2 = t_datetime(dt1%timestamp())

      ! Both should represent the same point in time
      call assert_equal(dt1%year(), dt2%year(), "From timestamp - year")
      call assert_equal(dt1%month(), dt2%month(), "From timestamp - month")
      call assert_equal(dt1%day(), dt2%day(), "From timestamp - day")
      call assert_equal(dt1%hour(), dt2%hour(), "From timestamp - hour")
      call assert_equal(dt1%minute(), dt2%minute(), "From timestamp - minute")
      call assert_equal(dt1%second(), dt2%second(), "From timestamp - second")
      call assert_equal(dt1%millisecond(), dt2%millisecond(), "From timestamp - millisecond")
      call assert_equal(dt1%timestamp(), dt2%timestamp(), "From timestamp - timestamp")
   end subroutine test_datetime_from_timestamp

   subroutine test_datetime_strptime()
      type(t_datetime) :: dt1, dt2, dt3

      ! Default format
      dt1 = t_datetime("2022-01-31 12:34:56", "%Y-%m-%d %H:%M:%S")

      call assert_equal(2022, dt1%year(), "Parse default format - year")
      call assert_equal(1, dt1%month(), "Parse default format - month")
      call assert_equal(31, dt1%day(), "Parse default format - day")
      call assert_equal(12, dt1%hour(), "Parse default format - hour")
      call assert_equal(34, dt1%minute(), "Parse default format - minute")
      call assert_equal(56, dt1%second(), "Parse default format - second")
      call assert_equal(0, dt1%millisecond(), "Parse default format - millisecond")

      ! Custom format
      dt2 = t_datetime("31/01/2022 12:34:56", "%d/%m/%Y %H:%M:%S")

      call assert_equal(2022, dt2%year(), "Parse custom format - year")
      call assert_equal(1, dt2%month(), "Parse custom format - month")
      call assert_equal(31, dt2%day(), "Parse custom format - day")
      call assert_equal(12, dt2%hour(), "Parse custom format - hour")
      call assert_equal(34, dt2%minute(), "Parse custom format - minute")
      call assert_equal(56, dt2%second(), "Parse custom format - second")
      call assert_equal(0, dt2%millisecond(), "Parse custom format - millisecond")

      ! With milliseconds
      dt3 = t_datetime("2022-01-31T12:34:56.789", "%Y-%m-%dT%H:%M:%S")

      call assert_equal(2022, dt3%year(), "Parse with ms - year")
      call assert_equal(1, dt3%month(), "Parse with ms - month")
      call assert_equal(31, dt3%day(), "Parse with ms - day")
      call assert_equal(12, dt3%hour(), "Parse with ms - hour")
      call assert_equal(34, dt3%minute(), "Parse with ms - minute")
      call assert_equal(56, dt3%second(), "Parse with ms - second")
      call assert_equal(789, dt3%millisecond(), "Parse with ms - millisecond")
   end subroutine test_datetime_strptime

   subroutine test_datetime_strftime()
      type(t_datetime) :: dt

      dt = t_datetime(2022, 1, 31, 12, 34, 56, 789)

      ! Default format
      call assert_equal("2022-01-31 12:34:56", dt%strftime("%Y-%m-%d %H:%M:%S"), "Format default")

      ! Custom format
      call assert_equal("31/01/2022 12:34:56", dt%strftime("%d/%m/%Y %H:%M:%S"), "Format custom")

      ! Format with milliseconds
      call assert_equal("2022-01-31 12:34:56.789", dt%strftime("%Y-%m-%d %H:%M:%S", .true.), "Format with ms")

      ! ISO string format
      call assert_equal("2022-01-31T12:34:56", dt%to_iso_string(), "ISO string format")

      ! Format with various specifiers
      call assert_equal("2022", dt%strftime("%Y"), "Format year")
      call assert_equal("01", dt%strftime("%m"), "Format month")
      call assert_equal("31", dt%strftime("%d"), "Format day")
      call assert_equal("12", dt%strftime("%H"), "Format hour")
      call assert_equal("34", dt%strftime("%M"), "Format minute")
      call assert_equal("56", dt%strftime("%S"), "Format second")
      call assert_equal("2022-01", dt%strftime("%Y-%m"), "Format year-month")
   end subroutine test_datetime_strftime

   subroutine test_datetime_now()
      type(t_datetime) :: dt_now

      dt_now = now()

      ! Basic validation that now() returns a sensible date (not 1970)
      call assert_true(dt_now%year() >= 2023, "Now() year should be recent")

      ! We can't make exact assertions about the current time since it changes
      ! But we can check that it's within reasonable bounds
      call assert_true(dt_now%month() >= 1 .and. dt_now%month() <= 12, "Now() month in valid range")
      call assert_true(dt_now%day() >= 1 .and. dt_now%day() <= 31, "Now() day in valid range")
      call assert_true(dt_now%hour() >= 0 .and. dt_now%hour() <= 23, "Now() hour in valid range")
      call assert_true(dt_now%minute() >= 0 .and. dt_now%minute() <= 59, "Now() minute in valid range")
      call assert_true(dt_now%second() >= 0 .and. dt_now%second() <= 59, "Now() second in valid range")
      call assert_true(dt_now%millisecond() >= 0 .and. dt_now%millisecond() <= 999, "Now() millisecond in valid range")
   end subroutine test_datetime_now

   subroutine test_datetime_addition_with_timedelta()
      use mod_datetime, only: t_datetime, t_timedelta, operator(+)
      implicit none
      type(t_datetime) :: dt, result
      type(t_timedelta) :: ts

      dt = t_datetime(2022, 1, 15, 12, 0, 0, 0)

      ! Add 10 days
      ts = t_timedelta(days=10)
      result = dt + ts

      call assert_equal(2022, result%year(), "Add 10 days - year")
      call assert_equal(1, result%month(), "Add 10 days - month")
      call assert_equal(25, result%day(), "Add 10 days - day")
      call assert_equal(12, result%hour(), "Add 10 days - hour")

      ! Add 12 hours
      ts = t_timedelta(hours=12)
      result = dt + ts

      call assert_equal(2022, result%year(), "Add 12 hours - year")
      call assert_equal(1, result%month(), "Add 12 hours - month")
      call assert_equal(16, result%day(), "Add 12 hours - day")
      call assert_equal(0, result%hour(), "Add 12 hours - hour")
   end subroutine test_datetime_addition_with_timedelta

   subroutine test_datetime_subtraction_with_timedelta()
      use mod_datetime, only: t_datetime, t_timedelta, operator(-)
      implicit none
      type(t_datetime) :: dt, result
      type(t_timedelta) :: ts

      dt = t_datetime(2022, 1, 15, 12, 0, 0, 0)

      ! Subtract 10 days
      ts = t_timedelta(days=10)
      result = dt - ts

      call assert_equal(2022, result%year(), "Subtract 10 days - year")
      call assert_equal(1, result%month(), "Subtract 10 days - month")
      call assert_equal(5, result%day(), "Subtract 10 days - day")
      call assert_equal(12, result%hour(), "Subtract 10 days - hour")

      ! Subtract 13 hours
      ts = t_timedelta(hours=13)
      result = dt - ts

      call assert_equal(2022, result%year(), "Subtract 13 hours - year")
      call assert_equal(1, result%month(), "Subtract 13 hours - month")
      call assert_equal(14, result%day(), "Subtract 13 hours - day")
      call assert_equal(23, result%hour(), "Subtract 13 hours - hour")
   end subroutine test_datetime_subtraction_with_timedelta

   subroutine test_datetime_difference()
      use, intrinsic :: iso_fortran_env, only: int64
      use mod_datetime, only: t_datetime, t_timedelta, operator(-)
      type(t_datetime) :: dt1, dt2
      type(t_timedelta) :: diff

      dt1 = t_datetime(2022, 1, 15, 12, 0, 0, 0)
      dt2 = t_datetime(2022, 1, 20, 18, 30, 0, 0)

      diff = dt2 - dt1

      call assert_equal(5_int64, diff%total_days(), "DateTime difference - total days")
      call assert_equal(5_int64*24 + 6, diff%total_hours(), "DateTime difference - total hours")
      call assert_equal((5_int64*24 + 6)*60 + 30, diff%total_minutes(), "DateTime difference - total minutes")
   end subroutine test_datetime_difference

   subroutine test_datetime_comparisons()
      use mod_datetime, only: t_datetime, operator(<), operator(>), &
                              operator(==), operator(/=), operator(<=), operator(>=)
      type(t_datetime) :: dt1, dt2, dt3

      dt1 = t_datetime(2022, 1, 15)
      dt2 = t_datetime(2022, 1, 20)
      dt3 = t_datetime(2022, 1, 15)

      ! Equality
      call assert_true(dt1 == dt3, "DateTime equality when equal")
      call assert_false(dt1 == dt2, "DateTime equality when different")

      ! Inequality
      call assert_true(dt1 /= dt2, "DateTime inequality when different")
      call assert_false(dt1 /= dt3, "DateTime inequality when equal")

      ! Less than
      call assert_true(dt1 < dt2, "DateTime less than when less")
      call assert_false(dt2 < dt1, "DateTime less than when greater")
      call assert_false(dt1 < dt3, "DateTime less than when equal")

      ! Greater than
      call assert_true(dt2 > dt1, "DateTime greater than when greater")
      call assert_false(dt1 > dt2, "DateTime greater than when less")
      call assert_false(dt1 > dt3, "DateTime greater than when equal")

      ! Less than or equal
      call assert_true(dt1 <= dt2, "DateTime less than or equal when less")
      call assert_true(dt1 <= dt3, "DateTime less than or equal when equal")
      call assert_false(dt2 <= dt1, "DateTime less than or equal when greater")

      ! Greater than or equal
      call assert_true(dt2 >= dt1, "DateTime greater than or equal when greater")
      call assert_true(dt1 >= dt3, "DateTime greater than or equal when equal")
      call assert_false(dt1 >= dt2, "DateTime greater than or equal when less")
   end subroutine test_datetime_comparisons

   subroutine test_datetime_date_wrapping()
      use mod_datetime, only: t_datetime, t_timedelta, operator(+)
      type(t_datetime) :: dt, result
      type(t_timedelta) :: ts

      ! Create a date and add enough time to wrap to the next month
      dt = t_datetime(2022, 1, 31)
      ts = t_timedelta(days=1)
      result = dt + ts

      call assert_equal(2022, result%year(), "Date wrapping month - year")
      call assert_equal(2, result%month(), "Date wrapping month - month")
      call assert_equal(1, result%day(), "Date wrapping month - day")

      ! Create a date and add enough time to wrap to the next year
      dt = t_datetime(2022, 12, 31)
      ts = t_timedelta(days=1)
      result = dt + ts

      call assert_equal(2023, result%year(), "Date wrapping year - year")
      call assert_equal(1, result%month(), "Date wrapping year - month")
      call assert_equal(1, result%day(), "Date wrapping year - day")
   end subroutine test_datetime_date_wrapping

   subroutine test_datetime_leap_year_handling()
      use mod_datetime, only: t_datetime, t_timedelta, operator(+)
      type(t_datetime) :: leap_day, next_year
      type(t_timedelta) :: ts

      ! Test February 29 in leap year
      leap_day = t_datetime(2020, 2, 29)

      ! Add one year
      ts = t_timedelta(days=366)
      next_year = leap_day + ts

      ! Should be March 1, 2021 (2020 was a leap year)
      call assert_equal(2021, next_year%year(), "Leap year handling - year")
      call assert_equal(3, next_year%month(), "Leap year handling - month")
      call assert_equal(1, next_year%day(), "Leap year handling - day")
   end subroutine test_datetime_leap_year_handling

   subroutine test_datetime_serialization_roundtrip()
      type(t_datetime) :: original, parsed
      character(len=64) :: iso

      original = t_datetime(2022, 3, 15, 14, 30, 45, 500)

      ! Convert to string
      iso = original%to_iso_string(.true.)

      ! Parse back with appropriate format
      parsed = t_datetime(iso, "%Y-%m-%dT%H:%M:%S")

      ! Should be the same timestamp
      call assert_equal(original%timestamp(), parsed%timestamp(), "Serialization roundtrip timestamp")
   end subroutine test_datetime_serialization_roundtrip

   subroutine test_datetime_extreme_values()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime, operator(<), operator(>)
      implicit none
      type(t_datetime) :: ancient, far_future

      ! Far past date
      ancient = t_datetime(1, 1, 1)

      ! Very far future date
      far_future = t_datetime(9999, 12, 31, 23, 59, 59, 999)

      call assert_true(ancient < far_future, "Extreme values comparison")
      call assert_equal(1, ancient%year(), "Ancient year")
      call assert_equal(9999, far_future%year(), "Far future year")
   end subroutine test_datetime_extreme_values

   subroutine test_datetime_invalid()
      use test_utils, only: assert_equal, assert_true, assert_false
      use mod_datetime, only: t_datetime, null_datetime
      implicit none
      type(t_datetime) :: dt_invalid, dt_valid

      ! Parse a date that should be invalid and check for error
      dt_invalid = t_datetime("applesauce")

      ! Check that the date is not valid
      call assert_false(dt_invalid%valid(), "Invalid date should not be valid")

      ! Create a valid date
      dt_valid = t_datetime(2022, 1, 31)
      call assert_true(dt_valid%valid(), "Valid date should be valid")

      dt_valid = t_datetime("2022-01-31 12:34:56", "%Y-%m-%d %H:%M:%S")
      call assert_true(dt_valid%valid(), "Valid date from string should be valid")

      dt_invalid = t_datetime("74/12/2022 12:34:56", "%y/%m/%d %H:%M:%S")
      call assert_false(dt_invalid%valid(), "Invalid date from string should not be valid")

      dt_valid = t_datetime("2025-01-01")
      call assert_true(dt_valid%valid(), "Valid date from string should be valid")

      dt_invalid = null_datetime()
      call assert_false(dt_invalid%valid(), "Null DateTime should not be valid")

   end subroutine test_datetime_invalid

   ! ====================================================
   ! Array-based parsing tests
   ! ====================================================

   subroutine test_datetime_array_parsing_basic()
      use test_utils, only: assert_equal, assert_true, assert_false
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: parsed_dt
      character(len=30) :: formats(3)

      ! Test basic array parsing with first format succeeding
      formats(1) = '%Y-%m-%d %H:%M:%S          '
      formats(2) = '%Y/%m/%d %H:%M:%S          '
      formats(3) = '%d-%m-%Y %H:%M:%S          '

      parsed_dt = t_datetime('2024-03-15 14:30:25', formats)

      call assert_true(parsed_dt%valid(), "Array parsing should succeed with first format")
      call assert_equal(parsed_dt%year(), 2024, "Year should be 2024")
      call assert_equal(parsed_dt%month(), 3, "Month should be 3")
      call assert_equal(parsed_dt%day(), 15, "Day should be 15")
      call assert_equal(parsed_dt%hour(), 14, "Hour should be 14")
      call assert_equal(parsed_dt%minute(), 30, "Minute should be 30")
      call assert_equal(parsed_dt%second(), 25, "Second should be 25")
   end subroutine test_datetime_array_parsing_basic

   subroutine test_datetime_array_parsing_second_format()
      use test_utils, only: assert_equal, assert_true, assert_false
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: parsed_dt
      character(len=30) :: formats(3)

      ! Test array parsing with second format succeeding
      formats(1) = '%Y-%m-%d %H:%M:%S          '
      formats(2) = '%Y/%m/%d %H:%M:%S          '
      formats(3) = '%d-%m-%Y %H:%M:%S          '

      parsed_dt = t_datetime('2024/03/15 14:30:25', formats)

      call assert_true(parsed_dt%valid(), "Array parsing should succeed with second format")
      call assert_equal(parsed_dt%year(), 2024, "Year should be 2024")
      call assert_equal(parsed_dt%month(), 3, "Month should be 3")
      call assert_equal(parsed_dt%day(), 15, "Day should be 15")
      call assert_equal(parsed_dt%hour(), 14, "Hour should be 14")
      call assert_equal(parsed_dt%minute(), 30, "Minute should be 30")
      call assert_equal(parsed_dt%second(), 25, "Second should be 25")
   end subroutine test_datetime_array_parsing_second_format

   subroutine test_datetime_array_parsing_third_format()
      use test_utils, only: assert_equal, assert_true, assert_false
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: parsed_dt
      character(len=30) :: formats(3)

      ! Test array parsing with third format succeeding
      formats(1) = '%Y-%m-%d %H:%M:%S          '
      formats(2) = '%Y/%m/%d %H:%M:%S          '
      formats(3) = '%d-%m-%Y %H:%M:%S          '

      parsed_dt = t_datetime('15-03-2024 14:30:25', formats)

      call assert_true(parsed_dt%valid(), "Array parsing should succeed with third format")
      call assert_equal(parsed_dt%year(), 2024, "Year should be 2024")
      call assert_equal(parsed_dt%month(), 3, "Month should be 3")
      call assert_equal(parsed_dt%day(), 15, "Day should be 15")
      call assert_equal(parsed_dt%hour(), 14, "Hour should be 14")
      call assert_equal(parsed_dt%minute(), 30, "Minute should be 30")
      call assert_equal(parsed_dt%second(), 25, "Second should be 25")
   end subroutine test_datetime_array_parsing_third_format

   subroutine test_datetime_array_parsing_compact_format()
      use test_utils, only: assert_equal, assert_true, assert_false
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: parsed_dt
      character(len=30) :: formats(4)

      ! Test compact format parsing
      formats(1) = '%Y-%m-%d %H:%M:%S          '
      formats(2) = '%Y/%m/%d %H:%M:%S          '
      formats(3) = '%d-%m-%Y %H:%M:%S          '
      formats(4) = '%Y%m%d%H%M%S               '

      parsed_dt = t_datetime('20240315143025', formats)

      call assert_true(parsed_dt%valid(), "Array parsing should succeed with compact format")
      call assert_equal(parsed_dt%year(), 2024, "Year should be 2024")
      call assert_equal(parsed_dt%month(), 3, "Month should be 3")
      call assert_equal(parsed_dt%day(), 15, "Day should be 15")
      call assert_equal(parsed_dt%hour(), 14, "Hour should be 14")
      call assert_equal(parsed_dt%minute(), 30, "Minute should be 30")
      call assert_equal(parsed_dt%second(), 25, "Second should be 25")
   end subroutine test_datetime_array_parsing_compact_format

   subroutine test_datetime_array_parsing_iso_variations()
      use test_utils, only: assert_equal, assert_true, assert_false
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: parsed_dt
      character(len=30) :: formats(4)

      ! Test ISO format variations
      formats(1) = '%Y-%m-%dT%H:%M:%SZ         '
      formats(2) = '%Y-%m-%dT%H:%M:%S          '
      formats(3) = '%Y-%m-%d %H:%M:%S          '
      formats(4) = '%Y/%m/%d %H:%M:%S          '

      parsed_dt = t_datetime('2024-03-15T14:30:25Z', formats)

      call assert_true(parsed_dt%valid(), "Array parsing should succeed with ISO Z format")
      call assert_equal(parsed_dt%year(), 2024, "Year should be 2024")
      call assert_equal(parsed_dt%month(), 3, "Month should be 3")
      call assert_equal(parsed_dt%day(), 15, "Day should be 15")
      call assert_equal(parsed_dt%hour(), 14, "Hour should be 14")
      call assert_equal(parsed_dt%minute(), 30, "Minute should be 30")
      call assert_equal(parsed_dt%second(), 25, "Second should be 25")
   end subroutine test_datetime_array_parsing_iso_variations

   subroutine test_datetime_array_parsing_dot_separated()
      use test_utils, only: assert_equal, assert_true, assert_false
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: parsed_dt
      character(len=30) :: formats(4)

      ! Test dot-separated format parsing
      formats(1) = '%Y-%m-%d %H:%M:%S          '
      formats(2) = '%Y/%m/%d %H:%M:%S          '
      formats(3) = '%d-%m-%Y %H:%M:%S          '
      formats(4) = '%Y.%m.%d %H:%M:%S          '

      parsed_dt = t_datetime('2024.03.15 14:30:25', formats)

      call assert_true(parsed_dt%valid(), "Array parsing should succeed with dot-separated format")
      call assert_equal(parsed_dt%year(), 2024, "Year should be 2024")
      call assert_equal(parsed_dt%month(), 3, "Month should be 3")
      call assert_equal(parsed_dt%day(), 15, "Day should be 15")
      call assert_equal(parsed_dt%hour(), 14, "Hour should be 14")
      call assert_equal(parsed_dt%minute(), 30, "Minute should be 30")
      call assert_equal(parsed_dt%second(), 25, "Second should be 25")
   end subroutine test_datetime_array_parsing_dot_separated

   subroutine test_datetime_array_parsing_failures()
      use test_utils, only: assert_equal, assert_true, assert_false
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: parsed_dt
      character(len=30) :: formats(3)

      ! Test parsing failure when no formats match
      formats(1) = '%Y-%m-%d %H:%M:%S          '
      formats(2) = '%Y/%m/%d %H:%M:%S          '
      formats(3) = '%d-%m-%Y %H:%M:%S          '

      parsed_dt = t_datetime('invalid-date-string', formats)

      call assert_false(parsed_dt%valid(), "Array parsing should fail with invalid date string")
   end subroutine test_datetime_array_parsing_failures

   subroutine test_datetime_array_parsing_single_format()
      use test_utils, only: assert_equal, assert_true, assert_false
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: parsed_dt
      character(len=30) :: formats(1)

      ! Test array parsing with single format
      formats(1) = '%Y-%m-%d %H:%M:%S          '

      parsed_dt = t_datetime('2024-03-15 14:30:25', formats)

      call assert_true(parsed_dt%valid(), "Array parsing should succeed with single format")
      call assert_equal(parsed_dt%year(), 2024, "Year should be 2024")
      call assert_equal(parsed_dt%month(), 3, "Month should be 3")
      call assert_equal(parsed_dt%day(), 15, "Day should be 15")
      call assert_equal(parsed_dt%hour(), 14, "Hour should be 14")
      call assert_equal(parsed_dt%minute(), 30, "Minute should be 30")
      call assert_equal(parsed_dt%second(), 25, "Second should be 25")
   end subroutine test_datetime_array_parsing_single_format

   subroutine test_datetime_array_parsing_date_only()
      use test_utils, only: assert_equal, assert_true, assert_false
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: parsed_dt
      character(len=30) :: formats(2)

      ! Test date-only format parsing (time should default to 00:00:00)
      formats(1) = '%Y-%m-%d %H:%M:%S          ' ! This won't match (has time)
      formats(2) = '%Y-%m-%d                   ' ! This will match

      parsed_dt = t_datetime('2024-03-15', formats)

      call assert_true(parsed_dt%valid(), "Array parsing should succeed with date-only format")
      call assert_equal(parsed_dt%year(), 2024, "Year should be 2024")
      call assert_equal(parsed_dt%month(), 3, "Month should be 3")
      call assert_equal(parsed_dt%day(), 15, "Day should be 15")
      call assert_equal(parsed_dt%hour(), 0, "Hour should default to 0")
      call assert_equal(parsed_dt%minute(), 0, "Minute should default to 0")
      call assert_equal(parsed_dt%second(), 0, "Second should default to 0")
   end subroutine test_datetime_array_parsing_date_only

   subroutine test_datetime_auto_with_fallback()
      use test_utils, only: assert_equal, assert_true, assert_false
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: parsed_dt
      character(len=30) :: fallback_formats(2)

      ! Test auto-detection with fallback formats
      fallback_formats(1) = '%d.%m.%Y %H:%M:%S         '
      fallback_formats(2) = '%Y.%m.%dT%H:%M:%S         '

      ! This format should not be auto-detected, so it should use fallback
      parsed_dt = t_datetime('2024.03.15T14:30:25', fallback_formats)

      call assert_true(parsed_dt%valid(), "Auto with fallback should succeed")
      call assert_equal(parsed_dt%year(), 2024, "Year should be 2024")
      call assert_equal(parsed_dt%month(), 3, "Month should be 3")
      call assert_equal(parsed_dt%day(), 15, "Day should be 15")
      call assert_equal(parsed_dt%hour(), 14, "Hour should be 14")
      call assert_equal(parsed_dt%minute(), 30, "Minute should be 30")
      call assert_equal(parsed_dt%second(), 25, "Second should be 25")
   end subroutine test_datetime_auto_with_fallback

   ! ====================================================
   ! Comprehensive Automatic Format Detection Tests for Fortran
   ! ====================================================

   subroutine test_datetime_auto_format_coverage_1()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: result

      ! Test Format 1: %Y-%m-%d %H:%M:%S (standard format)
      result = t_datetime('2024-03-15 14:30:25')
      call assert_true(result%valid(), "Auto-detection should succeed for standard format")
      call assert_equal(result%year(), 2024, "Year should be 2024")
      call assert_equal(result%month(), 3, "Month should be 3")
      call assert_equal(result%day(), 15, "Day should be 15")
      call assert_equal(result%hour(), 14, "Hour should be 14")
      call assert_equal(result%minute(), 30, "Minute should be 30")
      call assert_equal(result%second(), 25, "Second should be 25")
   end subroutine test_datetime_auto_format_coverage_1

   subroutine test_datetime_auto_format_coverage_2()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: result

      ! Test Format 2: %Y-%m-%dT%H:%M:%SZ (ISO with timezone)
      result = t_datetime('2024-03-15T14:30:25Z')
      call assert_true(result%valid(), "Auto-detection should succeed for ISO Z format")
      call assert_equal(result%year(), 2024, "Year should be 2024")
      call assert_equal(result%month(), 3, "Month should be 3")
      call assert_equal(result%day(), 15, "Day should be 15")
      call assert_equal(result%hour(), 14, "Hour should be 14")
      call assert_equal(result%minute(), 30, "Minute should be 30")
      call assert_equal(result%second(), 25, "Second should be 25")
   end subroutine test_datetime_auto_format_coverage_2

   subroutine test_datetime_auto_format_coverage_3()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: result

      ! Test Format 3: %Y-%m-%dT%H:%M:%S (ISO format)
      result = t_datetime('2024-03-15T14:30:25')
      call assert_true(result%valid(), "Auto-detection should succeed for ISO format")
      call assert_equal(result%year(), 2024, "Year should be 2024")
      call assert_equal(result%month(), 3, "Month should be 3")
      call assert_equal(result%day(), 15, "Day should be 15")
      call assert_equal(result%hour(), 14, "Hour should be 14")
      call assert_equal(result%minute(), 30, "Minute should be 30")
      call assert_equal(result%second(), 25, "Second should be 25")
   end subroutine test_datetime_auto_format_coverage_3

   subroutine test_datetime_auto_format_coverage_4()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: result

      ! Test Format 4: %Y/%m/%d %H:%M:%S (slash format)
      result = t_datetime('2024/03/15 14:30:25')
      call assert_true(result%valid(), "Auto-detection should succeed for slash format")
      call assert_equal(result%year(), 2024, "Year should be 2024")
      call assert_equal(result%month(), 3, "Month should be 3")
      call assert_equal(result%day(), 15, "Day should be 15")
      call assert_equal(result%hour(), 14, "Hour should be 14")
      call assert_equal(result%minute(), 30, "Minute should be 30")
      call assert_equal(result%second(), 25, "Second should be 25")
   end subroutine test_datetime_auto_format_coverage_4

   subroutine test_datetime_auto_format_coverage_5()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: result

      ! Test Format 5: %Y.%m.%d %H:%M:%S (dot format)
      result = t_datetime('2024.03.15 14:30:25')
      call assert_true(result%valid(), "Auto-detection should succeed for dot format")
      call assert_equal(result%year(), 2024, "Year should be 2024")
      call assert_equal(result%month(), 3, "Month should be 3")
      call assert_equal(result%day(), 15, "Day should be 15")
      call assert_equal(result%hour(), 14, "Hour should be 14")
      call assert_equal(result%minute(), 30, "Minute should be 30")
      call assert_equal(result%second(), 25, "Second should be 25")
   end subroutine test_datetime_auto_format_coverage_5

   subroutine test_datetime_auto_format_coverage_6()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: result

      ! Test Format 6: %Y%m%d%H%M%S (compact format)
      result = t_datetime('20240315143025')
      call assert_true(result%valid(), "Auto-detection should succeed for compact format")
      call assert_equal(result%year(), 2024, "Year should be 2024")
      call assert_equal(result%month(), 3, "Month should be 3")
      call assert_equal(result%day(), 15, "Day should be 15")
      call assert_equal(result%hour(), 14, "Hour should be 14")
      call assert_equal(result%minute(), 30, "Minute should be 30")
      call assert_equal(result%second(), 25, "Second should be 25")
   end subroutine test_datetime_auto_format_coverage_6

   subroutine test_datetime_auto_format_coverage_7()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: result

      ! Test Format 7: %Y/%m/%d %H:%M (no seconds)
      result = t_datetime('2024/03/15 14:30')
      call assert_true(result%valid(), "Auto-detection should succeed for format without seconds")
      call assert_equal(result%year(), 2024, "Year should be 2024")
      call assert_equal(result%month(), 3, "Month should be 3")
      call assert_equal(result%day(), 15, "Day should be 15")
      call assert_equal(result%hour(), 14, "Hour should be 14")
      call assert_equal(result%minute(), 30, "Minute should be 30")
      call assert_equal(result%second(), 0, "Second should default to 0")
   end subroutine test_datetime_auto_format_coverage_7

   subroutine test_datetime_auto_format_coverage_8()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: result

      ! Test Format 8: %Y-%m-%d (date only)
      result = t_datetime('2024-03-15')
      call assert_true(result%valid(), "Auto-detection should succeed for date-only format")
      call assert_equal(result%year(), 2024, "Year should be 2024")
      call assert_equal(result%month(), 3, "Month should be 3")
      call assert_equal(result%day(), 15, "Day should be 15")
      call assert_equal(result%hour(), 0, "Hour should default to 0")
      call assert_equal(result%minute(), 0, "Minute should default to 0")
      call assert_equal(result%second(), 0, "Second should default to 0")
   end subroutine test_datetime_auto_format_coverage_8

   subroutine test_datetime_auto_format_coverage_9()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: result

      ! Test Format 9: %Y/%m/%d (date only with slashes)
      result = t_datetime('2024/03/15')
      call assert_true(result%valid(), "Auto-detection should succeed for slash date-only")
      call assert_equal(result%year(), 2024, "Year should be 2024")
      call assert_equal(result%month(), 3, "Month should be 3")
      call assert_equal(result%day(), 15, "Day should be 15")
      call assert_equal(result%hour(), 0, "Hour should default to 0")
      call assert_equal(result%minute(), 0, "Minute should default to 0")
      call assert_equal(result%second(), 0, "Second should default to 0")
   end subroutine test_datetime_auto_format_coverage_9

   subroutine test_datetime_auto_format_coverage_10()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: result

      ! Test Format 10: %Y.%m.%d (date only with dots)
      result = t_datetime('2024.03.15')
      call assert_true(result%valid(), "Auto-detection should succeed for dot date-only")
      call assert_equal(result%year(), 2024, "Year should be 2024")
      call assert_equal(result%month(), 3, "Month should be 3")
      call assert_equal(result%day(), 15, "Day should be 15")
      call assert_equal(result%hour(), 0, "Hour should default to 0")
      call assert_equal(result%minute(), 0, "Minute should default to 0")
      call assert_equal(result%second(), 0, "Second should default to 0")
   end subroutine test_datetime_auto_format_coverage_10

   subroutine test_datetime_auto_format_coverage_11()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: result

      ! Test Format 11: %Y%m%d (compact date only)
      result = t_datetime('20240315')
      call assert_true(result%valid(), "Auto-detection should succeed for compact date-only")
      call assert_equal(result%year(), 2024, "Year should be 2024")
      call assert_equal(result%month(), 3, "Month should be 3")
      call assert_equal(result%day(), 15, "Day should be 15")
      call assert_equal(result%hour(), 0, "Hour should default to 0")
      call assert_equal(result%minute(), 0, "Minute should default to 0")
      call assert_equal(result%second(), 0, "Second should default to 0")
   end subroutine test_datetime_auto_format_coverage_11

   subroutine test_datetime_auto_with_milliseconds()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: result

      ! Test standard format with milliseconds
      result = t_datetime('2024-03-15 14:30:25.123')
      call assert_true(result%valid(), "Auto-detection should succeed with milliseconds")
      call assert_equal(result%year(), 2024, "Year should be 2024")
      call assert_equal(result%month(), 3, "Month should be 3")
      call assert_equal(result%day(), 15, "Day should be 15")
      call assert_equal(result%hour(), 14, "Hour should be 14")
      call assert_equal(result%minute(), 30, "Minute should be 30")
      call assert_equal(result%second(), 25, "Second should be 25")
      call assert_equal(result%millisecond(), 123, "Millisecond should be 123")

      ! Test ISO format with milliseconds
      result = t_datetime('2024-03-15T14:30:25.456')
      call assert_true(result%valid(), "Auto-detection should succeed with ISO milliseconds")
      call assert_equal(result%year(), 2024, "Year should be 2024")
      call assert_equal(result%month(), 3, "Month should be 3")
      call assert_equal(result%day(), 15, "Day should be 15")
      call assert_equal(result%hour(), 14, "Hour should be 14")
      call assert_equal(result%minute(), 30, "Minute should be 30")
      call assert_equal(result%second(), 25, "Second should be 25")
      call assert_equal(result%millisecond(), 456, "Millisecond should be 456")

      ! Test slash format with milliseconds
      result = t_datetime('2024/03/15 14:30:25.789')
      call assert_true(result%valid(), "Auto-detection should succeed with slash milliseconds")
      call assert_equal(result%year(), 2024, "Year should be 2024")
      call assert_equal(result%month(), 3, "Month should be 3")
      call assert_equal(result%day(), 15, "Day should be 15")
      call assert_equal(result%hour(), 14, "Hour should be 14")
      call assert_equal(result%minute(), 30, "Minute should be 30")
      call assert_equal(result%second(), 25, "Second should be 25")
      call assert_equal(result%millisecond(), 789, "Millisecond should be 789")
   end subroutine test_datetime_auto_with_milliseconds

   subroutine test_datetime_auto_format_precedence()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: result

      ! Test that most specific format wins - full datetime vs date-only
      result = t_datetime('2024-03-15 14:30:25')
      call assert_true(result%valid(), "Full datetime should be parsed correctly")
      call assert_equal(result%hour(), 14, "Hour should be parsed from full format")
      call assert_equal(result%minute(), 30, "Minute should be parsed from full format")
      call assert_equal(result%second(), 25, "Second should be parsed from full format")

      ! Test ISO with Z takes precedence
      result = t_datetime('2024-03-15T14:30:25Z')
      call assert_true(result%valid(), "ISO with Z should be parsed correctly")
      call assert_equal(result%year(), 2024, "Year should be 2024")
      call assert_equal(result%month(), 3, "Month should be 3")
      call assert_equal(result%day(), 15, "Day should be 15")

      ! Test minutes-precision format when seconds missing
      result = t_datetime('2024/03/15 14:30')
      call assert_true(result%valid(), "Minutes precision should be parsed correctly")
      call assert_equal(result%year(), 2024, "Year should be 2024")
      call assert_equal(result%month(), 3, "Month should be 3")
      call assert_equal(result%day(), 15, "Day should be 15")
      call assert_equal(result%hour(), 14, "Hour should be 14")
      call assert_equal(result%minute(), 30, "Minute should be 30")
      call assert_equal(result%second(), 0, "Second should default to 0")
   end subroutine test_datetime_auto_format_precedence

   subroutine test_datetime_auto_edge_cases()
      use test_utils, only: assert_equal, assert_true, assert_false
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: result

      ! Test ambiguous date handling - year boundaries
      result = t_datetime('2024-01-01')
      call assert_true(result%valid(), "Year boundary date should parse correctly")
      call assert_equal(result%year(), 2024, "Year should be 2024")
      call assert_equal(result%month(), 1, "Month should be 1")
      call assert_equal(result%day(), 1, "Day should be 1")

      ! Test leap year dates
      result = t_datetime('2024-02-29')
      call assert_true(result%valid(), "Leap year date should parse correctly")
      call assert_equal(result%year(), 2024, "Year should be 2024")
      call assert_equal(result%month(), 2, "Month should be 2")
      call assert_equal(result%day(), 29, "Day should be 29")

      ! Test end of month dates
      result = t_datetime('2024-01-31')
      call assert_true(result%valid(), "End of month date should parse correctly")
      call assert_equal(result%year(), 2024, "Year should be 2024")
      call assert_equal(result%month(), 1, "Month should be 1")
      call assert_equal(result%day(), 31, "Day should be 31")

      ! Test mixed separators should fail
      result = t_datetime('2024-03/15 14:30:25')
      call assert_false(result%valid(), "Mixed separators should fail")

      ! Test partial format strings should fail
      result = t_datetime('2024-03')
      call assert_false(result%valid(), "Partial format should fail")

   end subroutine test_datetime_auto_edge_cases

   subroutine test_datetime_auto_invalid_inputs()
      use test_utils, only: assert_false
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: result

      ! Test completely invalid strings
      result = t_datetime('not a date')
      call assert_false(result%valid(), "Invalid string should fail")

      result = t_datetime('')
      call assert_false(result%valid(), "Empty string should fail")

      result = t_datetime('abc123def')
      call assert_false(result%valid(), "Random string should fail")

      ! Test invalid date values
      result = t_datetime('2024-13-15')
      call assert_false(result%valid(), "Invalid month should fail")

      result = t_datetime('2024-02-30')
      call assert_false(result%valid(), "Invalid day for February should fail")

      result = t_datetime('2024-04-31')
      call assert_false(result%valid(), "Invalid day for April should fail")

      ! Test invalid time value rejection - parsing library should reject out-of-range values
      result = t_datetime('2024-03-15 25:30:25')
      call assert_false(result%valid(), "25 hours should fail to parse")

      result = t_datetime('2024-03-15 14:60:25')
      call assert_false(result%valid(), "60 minutes should fail to parse")

      result = t_datetime('2024-03-15 14:30:60')
      call assert_false(result%valid(), "60 seconds should fail to parse")

      ! Test malformed but partially parseable - use truly malformed strings
      result = t_datetime('2024/03-15')
      call assert_false(result%valid(), "Mixed separators should fail")

      result = t_datetime('2024-')
      call assert_false(result%valid(), "Incomplete date should fail")

      result = t_datetime('202a-03-15')
      call assert_false(result%valid(), "Non-numeric year should fail")
   end subroutine test_datetime_auto_invalid_inputs

   subroutine test_datetime_julian_day_number()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      use, intrinsic :: iso_fortran_env, only: int64
      implicit none
      type(t_datetime) :: dt

      ! Test known Julian Day Numbers for astronomical reference dates

      ! January 1, 2000 12:00 UTC = JDN 2451545
      dt = t_datetime(2000, 1, 1, 12, 0, 0)
      call assert_equal(2451545_int64, dt%julian_day_number(), "JDN for Y2K reference")

      ! Unix epoch: January 1, 1970 00:00 UTC = JDN 2440588
      dt = t_datetime(1970, 1, 1, 0, 0, 0)
      call assert_equal(2440588_int64, dt%julian_day_number(), "JDN for Unix epoch")

      ! GPS epoch: January 6, 1980 00:00 UTC = JDN 2444245
      dt = t_datetime(1980, 1, 6, 0, 0, 0)
      call assert_equal(2444245_int64, dt%julian_day_number(), "JDN for GPS epoch")

      ! Gregorian calendar adoption: October 15, 1582 = JDN 2299161
      dt = t_datetime(1582, 10, 15, 0, 0, 0)
      call assert_equal(2299161_int64, dt%julian_day_number(), "JDN for Gregorian calendar adoption")

      ! Test that time components don't affect JDN (should be same for any time on same date)
      dt = t_datetime(2000, 1, 1, 0, 0, 0)
      call assert_equal(2451545_int64, dt%julian_day_number(), "JDN same for midnight")

      dt = t_datetime(2000, 1, 1, 23, 59, 59)
      call assert_equal(2451545_int64, dt%julian_day_number(), "JDN same for end of day")
   end subroutine test_datetime_julian_day_number

   subroutine test_datetime_julian_day()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: dt
      real(kind=8) :: jd, expected_jd

      ! Test fractional Julian Day calculation

      ! January 1, 2000 12:00 UTC = JD 2451545.0 (noon is 0.0 fractional day)
      dt = t_datetime(2000, 1, 1, 12, 0, 0)
      jd = dt%julian_day()
      call assert_equal(2451545.0d0, jd, "JD for Y2K at noon")

      ! January 1, 2000 00:00 UTC = JD 2451544.5 (midnight is -0.5 fractional day)
      dt = t_datetime(2000, 1, 1, 0, 0, 0)
      jd = dt%julian_day()
      call assert_equal(2451544.5d0, jd, "JD for Y2K at midnight")

      ! January 1, 2000 18:00 UTC = JD 2451545.25 (6 PM is +0.25 fractional day)
      dt = t_datetime(2000, 1, 1, 18, 0, 0)
      jd = dt%julian_day()
      call assert_equal(2451545.25d0, jd, "JD for Y2K at 6 PM")

      ! January 2, 2000 00:00 UTC = JD 2451545.5 (next day midnight)
      dt = t_datetime(2000, 1, 2, 0, 0, 0)
      jd = dt%julian_day()
      call assert_equal(2451545.5d0, jd, "JD for next day midnight")

      ! Test with minutes and seconds
      ! January 1, 2000 12:30:00 UTC = JD 2451545.0 + 30/1440 = 2451545.020833...
      dt = t_datetime(2000, 1, 1, 12, 30, 0)
      jd = dt%julian_day()
      expected_jd = 2451545.0d0 + 30.0d0/1440.0d0
      call assert_equal(expected_jd, jd, "JD for Y2K at 12:30")

      ! Test with milliseconds
      ! January 1, 2000 12:00:00.500 UTC = JD 2451545.0 + 500/(24*60*60*1000)
      dt = t_datetime(2000, 1, 1, 12, 0, 0, 500)
      jd = dt%julian_day()
      expected_jd = 2451545.0d0 + 500.0d0/86400000.0d0
      call assert_equal(expected_jd, jd, "JD for Y2K with milliseconds")
   end subroutine test_datetime_julian_day

   subroutine test_datetime_julian_day_consistency()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      use, intrinsic :: iso_fortran_env, only: int64
      implicit none
      type(t_datetime) :: dt
      integer(kind=8) :: jdn
      real(kind=8) :: jd, fractional_part

      ! Test that julian_day() = julian_day_number() + fractional_part
      dt = t_datetime(2020, 3, 15, 14, 30, 45, 123)

      jdn = dt%julian_day_number()
      jd = dt%julian_day()
      fractional_part = jd - real(jdn, kind=8)

      ! For 14:30:45.123, fractional part should be (14.5-12)/24 + 45/86400 + 123/86400000
      ! = 2.5/24 + 45/86400 + 123/86400000 = 0.104340...
      call assert_true(fractional_part > 0.104d0 .and. fractional_part < 0.105d0, &
                       "Fractional part should be approximately 0.1043")

      ! Test noon gives fractional part of 0.0
      dt = t_datetime(2020, 3, 15, 12, 0, 0, 0)
      jd = dt%julian_day()
      jdn = dt%julian_day_number()
      fractional_part = jd - real(jdn, kind=8)
      call assert_equal(0.0d0, fractional_part, "Noon should give fractional part of 0.0")

      ! Test midnight gives fractional part of -0.5
      dt = t_datetime(2020, 3, 15, 0, 0, 0, 0)
      jd = dt%julian_day()
      jdn = dt%julian_day_number()
      fractional_part = jd - real(jdn, kind=8)
      call assert_equal(-0.5d0, fractional_part, "Midnight should give fractional part of -0.5")
   end subroutine test_datetime_julian_day_consistency

   subroutine test_datetime_julian_day_edge_cases()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      use, intrinsic :: iso_fortran_env, only: int64
      implicit none
      type(t_datetime) :: dt
      integer(kind=8) :: jdn
      real(kind=8) :: jd

      ! Test early date
      dt = t_datetime(1, 1, 1, 12, 0, 0)
      jdn = dt%julian_day_number()
      jd = dt%julian_day()
      call assert_true(jdn > 0, "JDN should be positive for year 1")
      call assert_true(jd > 0.0d0, "JD should be positive for year 1")

      ! Test future date
      dt = t_datetime(3000, 12, 31, 12, 0, 0)
      jdn = dt%julian_day_number()
      jd = dt%julian_day()
      call assert_true(jdn > 2451545_int64, "JDN should be greater than Y2K for year 3000")
      call assert_true(jd > 2451545.0d0, "JD should be greater than Y2K for year 3000")

      ! Test leap year: February 29, 2000
      dt = t_datetime(2000, 2, 29, 12, 0, 0)
      jdn = dt%julian_day_number()
      call assert_true(jdn > 0, "JDN should be valid for leap day")

      ! Test sequential days increment by 1
      dt = t_datetime(2000, 1, 1, 12, 0, 0)
      jdn = dt%julian_day_number()

      dt = t_datetime(2000, 1, 2, 12, 0, 0)
      call assert_equal(jdn + 1_int64, dt%julian_day_number(), "Sequential days should increment JDN by 1")
   end subroutine test_datetime_julian_day_edge_cases

   subroutine test_datetime_julian_century()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: dt
      real(kind=8) :: jc, expected_jc

      ! Test known Julian Century values for astronomical reference dates

      ! J2000.0 epoch: January 1, 2000 12:00 UTC = Julian Century 0.0
      ! (JD 2451545.0 corresponds to J2000.0 epoch)
      dt = t_datetime(2000, 1, 1, 12, 0, 0)
      jc = dt%julian_century()
      call assert_equal(0.0d0, jc, "Julian Century for J2000.0 epoch")

      ! Test one century before J2000.0: January 1, 1900 12:00 UTC = JC -1.0
      dt = t_datetime(1900, 1, 1, 12, 0, 0)
      jc = dt%julian_century()
      call assert_equal(-1.0d0, jc, "Julian Century for 1900.0 (one century before J2000)")

      ! Test one century after J2000.0: January 1, 2100 12:00 UTC = JC 1.0
      dt = t_datetime(2100, 1, 1, 12, 0, 0)
      jc = dt%julian_century()
      call assert_equal(1.0d0, jc, "Julian Century for 2100.0 (one century after J2000)")

      ! Test half century after J2000.0: January 1, 2050 12:00 UTC = JC 0.5
      dt = t_datetime(2050, 1, 1, 12, 0, 0)
      jc = dt%julian_century()
      call assert_equal(0.5d0, jc, "Julian Century for 2050.0 (half century after J2000)")

      ! Test Unix epoch: January 1, 1970 00:00 UTC
      ! JD = 2440587.5, JC = (2440587.5 - 2451545.0) / 36525 = -0.300019...
      dt = t_datetime(1970, 1, 1, 0, 0, 0)
      jc = dt%julian_century()
      expected_jc = (2440587.5d0 - 2451545.0d0)/36525.0d0
      call assert_equal(expected_jc, jc, "Julian Century for Unix epoch")

      ! Test time components affect Julian Century calculation
      ! January 1, 2000 00:00 UTC = JD 2451544.5, JC = (2451544.5 - 2451545.0) / 36525
      dt = t_datetime(2000, 1, 1, 0, 0, 0)
      jc = dt%julian_century()
      expected_jc = (2451544.5d0 - 2451545.0d0)/36525.0d0
      call assert_equal(expected_jc, jc, "Julian Century for J2000 at midnight")

      ! Test time components with more precision
      ! January 1, 2000 18:00 UTC = JD 2451545.25, JC = (2451545.25 - 2451545.0) / 36525
      dt = t_datetime(2000, 1, 1, 18, 0, 0)
      jc = dt%julian_century()
      expected_jc = (2451545.25d0 - 2451545.0d0)/36525.0d0
      call assert_equal(expected_jc, jc, "Julian Century for J2000 at 6 PM")
   end subroutine test_datetime_julian_century

   subroutine test_datetime_julian_century_precision()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: dt
      real(kind=8) :: jc, expected_jc

      ! Test precision with milliseconds
      ! January 1, 2000 12:00:00.500 UTC
      dt = t_datetime(2000, 1, 1, 12, 0, 0, 500)
      jc = dt%julian_century()
      ! JD = 2451545.0 + 500/(24*60*60*1000) = 2451545.0 + 500/86400000
      expected_jc = (2451545.0d0 + 500.0d0/86400000.0d0 - 2451545.0d0)/36525.0d0
      call assert_equal(expected_jc, jc, "Julian Century with millisecond precision")

      ! Test that small time differences result in small JC differences
      ! January 1, 2000 12:30:45.123 UTC
      dt = t_datetime(2000, 1, 1, 12, 30, 45, 123)
      jc = dt%julian_century()
      ! Calculate expected JC for this precise time
      ! JD = 2451545.0 + (30*60 + 45 + 0.123)/(24*60*60) = 2451545.0 + 1845.123/86400
      expected_jc = (2451545.0d0 + 1845.123d0/86400.0d0 - 2451545.0d0)/36525.0d0
      call assert_equal(expected_jc, jc, "Julian Century with precise time components")

      ! Test that the Julian Century is continuous (no jumps between days)
      dt = t_datetime(1999, 12, 31, 23, 59, 59, 999)
      jc = dt%julian_century()
      call assert_true(jc < 0.0d0, "Julian Century should be negative before J2000.0")

      dt = t_datetime(2000, 1, 1, 0, 0, 0, 1)
      jc = dt%julian_century()
      call assert_true(jc < 0.0d0, "Julian Century should be slightly negative just after midnight on J2000.0")
   end subroutine test_datetime_julian_century_precision

   subroutine test_datetime_julian_century_edge_cases()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: dt
      real(kind=8) :: jc

      ! Test very early date
      dt = t_datetime(1, 1, 1, 12, 0, 0)
      jc = dt%julian_century()
      call assert_true(jc < -19.0d0, "Julian Century should be very negative for year 1")

      ! Test far future date
      dt = t_datetime(3000, 12, 31, 12, 0, 0)
      jc = dt%julian_century()
      call assert_true(jc > 10.0d0, "Julian Century should be large positive for year 3000")

      ! Test leap year date: February 29, 2000
      dt = t_datetime(2000, 2, 29, 12, 0, 0)
      jc = dt%julian_century()
      call assert_true(jc > 0.0d0 .and. jc < 0.2d0, "Julian Century should be small positive for leap day 2000")

      ! Test Gregorian calendar adoption: October 15, 1582
      dt = t_datetime(1582, 10, 15, 12, 0, 0)
      jc = dt%julian_century()
      call assert_true(jc < -4.0d0, "Julian Century should be negative for Gregorian calendar adoption")

      ! Test that sequential years show expected JC progression
      dt = t_datetime(2000, 1, 1, 12, 0, 0)
      jc = dt%julian_century()

      dt = t_datetime(2001, 1, 1, 12, 0, 0)
      call assert_true(dt%julian_century() > jc, "Julian Century should increase for next year")
      call assert_true((dt%julian_century() - jc) < 0.02d0, "One year should be approximately 0.01 Julian Century")
   end subroutine test_datetime_julian_century_edge_cases

   subroutine test_datetime_julian_century_consistency()
      use test_utils, only: assert_equal, assert_true
      use mod_datetime, only: t_datetime
      implicit none
      type(t_datetime) :: dt
      real(kind=8) :: jc, jd, expected_jc

      ! Test that julian_century() = (julian_day() - 2451545.0) / 36525.0
      dt = t_datetime(2020, 3, 15, 14, 30, 45, 123)
      jc = dt%julian_century()
      jd = dt%julian_day()
      expected_jc = (jd - 2451545.0d0)/36525.0d0

      call assert_equal(expected_jc, jc, "Julian Century should match formula using Julian Day")

      ! Test multiple dates for consistency
      dt = t_datetime(1995, 7, 4, 6, 12, 30, 456)
      jc = dt%julian_century()
      jd = dt%julian_day()
      expected_jc = (jd - 2451545.0d0)/36525.0d0
      call assert_equal(expected_jc, jc, "Julian Century consistency test for 1995")

      dt = t_datetime(2025, 11, 22, 18, 45, 12, 789)
      jc = dt%julian_century()
      jd = dt%julian_day()
      expected_jc = (jd - 2451545.0d0)/36525.0d0
      call assert_equal(expected_jc, jc, "Julian Century consistency test for 2025")

      ! Test that the relationship holds for the J2000.0 epoch exactly
      dt = t_datetime(2000, 1, 1, 12, 0, 0, 0)
      jc = dt%julian_century()
      jd = dt%julian_day()
      call assert_equal(2451545.0d0, jd, "Julian Day for J2000.0 should be exactly 2451545.0")
      call assert_equal(0.0d0, jc, "Julian Century for J2000.0 should be exactly 0.0")
   end subroutine test_datetime_julian_century_consistency

end module datetime_tests

program test_datetime
   use test_utils, only: run_test, print_test_summary, &
                         assert_equal, assert_true, assert_false, &
                         tests_run, tests_failed
   use timedelta_tests, only: test_timedelta_default, test_timedelta_components, &
                              test_timedelta_total_accessors, test_timedelta_addition, &
                              test_timedelta_subtraction, test_timedelta_multiplication, &
                              test_timedelta_division, test_timedelta_comparisons, &
                              test_timedelta_to_string, test_timedelta_negative_duration, &
                              test_timedelta_overflow_handling
   use datetime_tests, only: test_datetime_default, test_datetime_ymd, &
                             test_datetime_complete, test_datetime_from_timestamp, &
                             test_datetime_strptime, test_datetime_strftime, test_datetime_now, &
                             test_datetime_addition_with_timedelta, test_datetime_subtraction_with_timedelta, &
                             test_datetime_difference, test_datetime_comparisons, &
                             test_datetime_date_wrapping, test_datetime_leap_year_handling, &
                             test_datetime_serialization_roundtrip, test_datetime_extreme_values, &
                             test_datetime_invalid, test_datetime_array_parsing_basic, &
                             test_datetime_array_parsing_second_format, test_datetime_array_parsing_third_format, &
                             test_datetime_array_parsing_compact_format, test_datetime_array_parsing_iso_variations, &
                             test_datetime_array_parsing_dot_separated, test_datetime_array_parsing_failures, &
                             test_datetime_array_parsing_single_format, test_datetime_array_parsing_date_only, &
                             test_datetime_auto_with_fallback, test_datetime_auto_format_coverage_1, &
                             test_datetime_auto_format_coverage_2, test_datetime_auto_format_coverage_3, &
                             test_datetime_auto_format_coverage_4, test_datetime_auto_format_coverage_5, &
                             test_datetime_auto_format_coverage_6, test_datetime_auto_format_coverage_7, &
                             test_datetime_auto_format_coverage_8, test_datetime_auto_format_coverage_9, &
                             test_datetime_auto_format_coverage_10, test_datetime_auto_format_coverage_11, &
                             test_datetime_auto_with_milliseconds, test_datetime_auto_format_precedence, &
                             test_datetime_auto_edge_cases, test_datetime_auto_invalid_inputs, &
                             test_datetime_julian_day_number, test_datetime_julian_day, &
                             test_datetime_julian_day_consistency, test_datetime_julian_day_edge_cases
   implicit none

   integer :: exit_code

   ! Print test header
   write (*, '(A)') "----------------------------------------"
   write (*, '(A)') "Running DateTime Module Tests"
   write (*, '(A)') "----------------------------------------"
   write (*, '(A)') ""

   ! TimeDelta tests
   write (*, '(A)') "=========== TimeDelta Tests ==========="
   call run_test(test_timedelta_default, "TimeDelta Default Constructor")
   call run_test(test_timedelta_components, "TimeDelta Components Constructor")
   call run_test(test_timedelta_total_accessors, "TimeDelta Total Accessors")
   call run_test(test_timedelta_addition, "TimeDelta Addition")
   call run_test(test_timedelta_subtraction, "TimeDelta Subtraction")
   call run_test(test_timedelta_multiplication, "TimeDelta Multiplication")
   call run_test(test_timedelta_division, "TimeDelta Division")
   call run_test(test_timedelta_comparisons, "TimeDelta Comparisons")
   call run_test(test_timedelta_to_string, "TimeDelta ToString")
   call run_test(test_timedelta_negative_duration, "TimeDelta Negative Duration")
   call run_test(test_timedelta_overflow_handling, "TimeDelta Overflow Handling")

   ! DateTime tests
   write (*, '(A)') "=========== DateTime Tests ==========="
   call run_test(test_datetime_default, "DateTime Default Constructor")
   call run_test(test_datetime_ymd, "DateTime YMD Constructor")
   call run_test(test_datetime_complete, "DateTime Complete Constructor")
   call run_test(test_datetime_from_timestamp, "DateTime From Timestamp")
   call run_test(test_datetime_strptime, "DateTime strptime")
   call run_test(test_datetime_strftime, "DateTime strftime")
   call run_test(test_datetime_now, "DateTime Now")
   call run_test(test_datetime_addition_with_timedelta, "DateTime Addition with TimeDelta")
   call run_test(test_datetime_subtraction_with_timedelta, "DateTime Subtraction with TimeDelta")
   call run_test(test_datetime_difference, "DateTime Difference")
   call run_test(test_datetime_comparisons, "DateTime Comparisons")
   call run_test(test_datetime_date_wrapping, "DateTime Date Wrapping")
   call run_test(test_datetime_leap_year_handling, "DateTime Leap Year Handling")
   call run_test(test_datetime_serialization_roundtrip, "DateTime Serialization Roundtrip")
   call run_test(test_datetime_extreme_values, "DateTime Extreme Values")
   call run_test(test_datetime_invalid, "DateTime Invalid")

   ! Array-based parsing tests
   write (*, '(A)') ""
   write (*, '(A)') "====== DateTime Array Parsing Tests ======"
   call run_test(test_datetime_array_parsing_basic, "DateTime Array Parsing Basic")
   call run_test(test_datetime_array_parsing_second_format, "DateTime Array Parsing Second Format")
   call run_test(test_datetime_array_parsing_third_format, "DateTime Array Parsing Third Format")
   call run_test(test_datetime_array_parsing_compact_format, "DateTime Array Parsing Compact Format")
   call run_test(test_datetime_array_parsing_iso_variations, "DateTime Array Parsing ISO Variations")
   call run_test(test_datetime_array_parsing_dot_separated, "DateTime Array Parsing Dot Separated")
   call run_test(test_datetime_array_parsing_failures, "DateTime Array Parsing Failures")
   call run_test(test_datetime_array_parsing_single_format, "DateTime Array Parsing Single Format")
   call run_test(test_datetime_array_parsing_date_only, "DateTime Array Parsing Date Only")
   call run_test(test_datetime_auto_with_fallback, "DateTime Auto with Fallback")

   ! Comprehensive Automatic Format Detection Tests
   write (*, '(A)') ""
   write (*, '(A)') "====== DateTime Auto-Detection Tests ======"
   call run_test(test_datetime_auto_format_coverage_1, "Auto Format Coverage 1 (Standard)")
   call run_test(test_datetime_auto_format_coverage_2, "Auto Format Coverage 2 (ISO Z)")
   call run_test(test_datetime_auto_format_coverage_3, "Auto Format Coverage 3 (ISO)")
   call run_test(test_datetime_auto_format_coverage_4, "Auto Format Coverage 4 (Slash)")
   call run_test(test_datetime_auto_format_coverage_5, "Auto Format Coverage 5 (Dot)")
   call run_test(test_datetime_auto_format_coverage_6, "Auto Format Coverage 6 (Compact)")
   call run_test(test_datetime_auto_format_coverage_7, "Auto Format Coverage 7 (No Seconds)")
   call run_test(test_datetime_auto_format_coverage_8, "Auto Format Coverage 8 (Date Only)")
   call run_test(test_datetime_auto_format_coverage_9, "Auto Format Coverage 9 (Date Only Slash)")
   call run_test(test_datetime_auto_format_coverage_10, "Auto Format Coverage 10 (Date Only Dot)")
   call run_test(test_datetime_auto_format_coverage_11, "Auto Format Coverage 11 (Compact Date)")
   call run_test(test_datetime_auto_with_milliseconds, "Auto Detection With Milliseconds")
   call run_test(test_datetime_auto_format_precedence, "Auto Detection Format Precedence")
   call run_test(test_datetime_auto_edge_cases, "Auto Detection Edge Cases")
   call run_test(test_datetime_auto_invalid_inputs, "Auto Detection Invalid Inputs")

   ! Julian Day tests
   write (*, '(A)') ""
   write (*, '(A)') "========== DateTime Julian Day Tests =========="
   call run_test(test_datetime_julian_day_number, "DateTime Julian Day Number")
   call run_test(test_datetime_julian_day, "DateTime Julian Day")
   call run_test(test_datetime_julian_day_consistency, "DateTime Julian Day Consistency")
   call run_test(test_datetime_julian_day_edge_cases, "DateTime Julian Day Edge Cases")

   ! Print summary
   call print_test_summary()

   ! Set exit code based on test results
   if (tests_failed > 0) then
      exit_code = 1
   else
      exit_code = 0
   end if

   call exit(exit_code)
end program test_datetime
