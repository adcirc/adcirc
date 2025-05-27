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
      use mod_datetime, only: t_datetime
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

   end subroutine test_datetime_invalid

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
                             test_datetime_invalid
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
