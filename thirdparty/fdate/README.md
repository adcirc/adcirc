# FDate Library
[![CI](https://github.com/zcobell/fdate/actions/workflows/push.yaml/badge.svg)](https://github.com/zcobell/fdate/actions/workflows/push.yaml)
[![License: GPLv3](https://img.shields.io/badge/License-GPLv3-yellow.svg)](https://opensource.org/licenses/GPL-3.0)
[![CMake](https://img.shields.io/badge/CMake-3.14%2B-blue.svg)](https://cmake.org/)
[![Fortran](https://img.shields.io/badge/Fortran-2008%2B-blue.svg)](https://fortran-lang.org/)

A modern datetime library for Fortran that provides comprehensive date and time manipulation capabilities through a clean Fortran interface built on top of Howard Hinnant's date library, which was used as the basis for `std::chrono` date operations in C++20.

## Background

The FDate library bridges the gap between modern C++ datetime functionality and Fortran applications. While Fortran excels in scientific computing, it has traditionally lacked robust datetime manipulation capabilities. This library provides:

- Millisecond precision datetime objects
- TimeDelta objects for duration calculations
- Comprehensive formatting and parsing with customizable format strings
- Full operator overloading for natural datetime arithmetic
- Memory-safe design using value types instead of pointers
- Cross-platform compatibility

The library uses a design where both `t_DateTime` and `t_TimeDelta` objects are represented internally as 64-bit integers (milliseconds since epoch for `t_DateTime`, total milliseconds for `t_TimeDelta`), and avoids using any other in-memory structures in the interface between C++ and Fortran.

## Features

### DateTime Operations
- Create dates from components (year, month, day, hour, minute, second, millisecond)
- Parse dates from strings with custom format specifiers
- Parse dates using arrays of format options (tries each format until one succeeds)
- Parse with automatic format detection plus custom fallback formats
- Format dates to strings (including ISO 8601)
- Get current system time
- Extract individual components (year, month, day, etc.)
- Julian Day Number and Julian Day calculations for astronomical applications
- Full comparison operators (`==`, `<`, `>`, `<=`, `>=`, `/=`)
- Add/subtract TimeDeltas to/from DateTime objects

### TimeDelta Operations
- Create time deltas from various units (days, hours, minutes, seconds, milliseconds)
- Extract components and total values
- Arithmetic operations (addition, subtraction, multiplication, division)
- Full comparison operators (`==`, `<`, `>`, `<=`, `>=`, `/=`)
- String representation

## Installation

### Prerequisites
- CMake 3.14 or higher
- Modern C++ compiler supporting at least C++17
- Fortran 2008 compiler (gfortran, ifort, etc.)

### Building
```bash
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install
make && make install
```
This will build a shared library and place it at the specified installation location. From there, you
can use it in your Fortran projects by linking against the installed library and including the Fortran 
module file in the installed include directory.

### CMake Options
- `FDATE_ENABLE_TESTING`: Enable testing suite (default: OFF)

## Fortran Usage Examples

### Basic DateTime Creation and Manipulation

```fortran
program datetime_examples
   use mod_datetime, only: t_datetime, t_timedelta, now
   implicit none
   
   type(t_datetime) :: dt1, dt2, current_time
   type(t_timedelta) :: ts1, ts2, difference
   character(len=64) :: date_string
   
   ! Create a specific datetime: January 15, 2024, 14:30:45.123
   dt1 = t_datetime(2024, 1, 15, 14, 30, 45, 123)
   
   ! Create date-only (time defaults to 00:00:00.000)
   dt2 = t_datetime(2024, 12, 25)
   
   ! Get current system time
   current_time = now()
   
   ! Extract components
   write(*,*) 'Year:', dt1%year()
   write(*,*) 'Month:', dt1%month()
   write(*,*) 'Day:', dt1%day()
   write(*,*) 'Hour:', dt1%hour()
   write(*,*) 'Minute:', dt1%minute()
   write(*,*) 'Second:', dt1%second()
   write(*,*) 'Millisecond:', dt1%millisecond()
   
end program datetime_examples
```

### TimeDelta Operations

```fortran
program timedelta_examples
   use mod_datetime, only: t_timedelta, t_datetime, operator(+), operator(-), operator(*), operator(/)
   implicit none
   
   type(t_timedelta) :: ts1, ts2, ts_result
   type(t_datetime) :: dt_start, dt_end
   
   ! Create timedeltas from different units
   ts1 = t_timedelta(days=2, hours=3, minutes=45, seconds=30, milliseconds=500)
   ts2 = t_timedelta(hours=5)  ! 5 hours
   
   ! Arithmetic operations
   ts_result = ts1 + ts2      ! Addition
   ts_result = ts1 - ts2      ! Subtraction
   ts_result = ts1 * 2        ! Multiplication
   ts_result = ts1 / 2        ! Division
   
   ! Extract components
   write(*,*) 'Days:', ts1%days()
   write(*,*) 'Hours:', ts1%hours()
   write(*,*) 'Total hours:', ts1%total_hours()
   write(*,*) 'Total milliseconds:', ts1%total_milliseconds()
   
   ! Calculate difference between two dates
   dt_start = t_datetime(2024, 1, 1, 9, 0, 0)
   dt_end = t_datetime(2024, 1, 5, 17, 30, 0)
   ts_result = dt_end - dt_start
   
   write(*,*) 'Time difference:', ts_result%to_string()
   
end program timedelta_examples
```

### Date Arithmetic and Comparisons

```fortran
program date_arithmetic
   use mod_datetime, only: t_datetime, t_timedelta, operator(+), operator(-), operator(<), operator(>), operator(==)
   implicit none
   
   type(t_datetime) :: meeting_time, deadline, reminder
   type(t_timedelta) :: one_week, two_hours
   
   ! Set meeting time
   meeting_time = t_datetime(2024, 6, 15, 14, 0, 0)  ! June 15, 2PM
   
   ! Create timedelta objects
   one_week = t_timedelta(days=7)
   two_hours = t_timedelta(hours=2)
   
   ! Calculate deadline (one week after meeting)
   deadline = meeting_time + one_week
   
   ! Set reminder (2 hours before meeting)
   reminder = meeting_time - two_hours
   
   ! Comparisons
   if (reminder < meeting_time) then
      write(*,*) 'Reminder is set before the meeting'
   end if
   
   if (deadline > meeting_time) then
      write(*,*) 'Deadline is after the meeting'
   end if
   
   ! Check if dates are equal
   if (meeting_time == meeting_time) then
      write(*,*) 'Dates are identical'
   end if
   
end program date_arithmetic
```

### String Formatting and Parsing

```fortran
program string_operations
   use mod_datetime, only: t_datetime
   implicit none
   
   type(t_datetime) :: dt, parsed_dt
   character(len=64) :: formatted_string
   
   ! Create a datetime
   dt = t_datetime(2024, 3, 14, 9, 26, 53, 589)
   
   ! Format to different string representations
   formatted_string = dt%strftime('%Y-%m-%d %H:%M:%S')
   write(*,*) 'Standard format: ', trim(formatted_string)
   
   formatted_string = dt%strftime('%B %d, %Y at %I:%M %p')
   write(*,*) 'Verbose format: ', trim(formatted_string)
   
   ! ISO 8601 format
   formatted_string = dt%to_iso_string()
   write(*,*) 'ISO format: ', trim(formatted_string)
   
   ! Some Custom formats
   formatted_string = dt%strftime('%A, %B %d, %Y') ! Full weekday and month names
   write(*,*) 'Custom format: ', trim(formatted_string)
   formatted_string = dt%strftime('%Y-%m-%dT%H:%M:%S',.true.) ! ISO with milliseconds
   write(*,*) 'ISO with milliseconds: ', trim(formatted_string)
   
   ! Parse a date string
   parsed_dt = t_datetime('2024-03-14 09:26:53', 'auto') ! Note, 'auto' as a format string is provided 
                                                         ! as a default if no format is specified 
  
   ! Check if parsing was successful
   if (parsed_dt%valid()) then
       write(*,*) 'Parsed date: ', parsed_dt%strftime('%Y-%m-%d %H:%M:%S')
   else
       write(*,*) 'Failed to parse date string'
   end if
   
end program string_operations
```

### Array-Based Parsing with Multiple Format Options

The library supports parsing with multiple format options, trying each format in order until one succeeds:

```fortran
program array_parsing_examples
   use mod_datetime, only: t_datetime
   implicit none
   
   type(t_datetime) :: parsed_dt
   character(len=30) :: date_strings(3) = [ &
      '2024/03/14 09:26:53        ', &  ! Slash-separated format
      '14-Mar-2024 09:26:53       ', &  ! Month name format  
      '20240314092653             '  ]  ! Compact format
   character(len=30) :: formats(4) = [ &
      '%Y/%m/%d %H:%M:%S          ', &  ! Slash format
      '%d-%b-%Y %H:%M:%S          ', &  ! Month name format
      '%Y%m%d%H%M%S               ', &  ! Compact format
      '%Y-%m-%d %H:%M:%S          '  ]  ! Standard format
   integer :: i
   
   ! Try parsing with an array of format options
   do i = 1, size(date_strings)
      parsed_dt = t_datetime(trim(date_strings(i)), formats)
      
      if (parsed_dt%valid()) then
         write(*,*) 'Successfully parsed: ', trim(date_strings(i))
         write(*,*) 'Result: ', parsed_dt%strftime('%Y-%m-%d %H:%M:%S')
      else
         write(*,*) 'Failed to parse: ', trim(date_strings(i))
      end if
   end do

   ! Parse with auto-detection and custom fallbacks
   parsed_dt = t_datetime('2024.03.14T09:26:53', ['%Y.%m.%dT%H:%M:%S     ', &
                                                   '%Y.%m.%d %H:%M:%S     '])
   
end program array_parsing_examples
```

### Julian Day Calculations

The library provides Julian Day Number (JDN) and Julian Day (JD) calculations for astronomical applications:

```fortran
program julian_day_examples
   use mod_datetime, only: t_datetime
   implicit none
   
   type(t_datetime) :: dt
   integer(kind=8) :: jdn
   real(kind=8) :: jd
   
   ! Create datetime objects for well-known astronomical dates
   
   ! January 1, 2000 at noon UTC (J2000.0 epoch)
   dt = t_datetime(2000, 1, 1, 12, 0, 0)
   jdn = dt%julian_day_number()    ! Returns 2451545
   jd = dt%julian_day()            ! Returns 2451545.0
   
   write(*,*) 'Y2K noon JDN:', jdn
   write(*,*) 'Y2K noon JD:', jd
   
   ! January 1, 2000 at midnight UTC
   dt = t_datetime(2000, 1, 1, 0, 0, 0)
   jd = dt%julian_day()            ! Returns 2451544.5
   write(*,*) 'Y2K midnight JD:', jd
   
   ! Current time with fractional day
   dt = t_datetime(2025, 7, 24, 15, 30, 45, 500)
   jdn = dt%julian_day_number()    ! Integer days since JD epoch
   jd = dt%julian_day()            ! Fractional days including time
   
   write(*,*) 'Current JDN:', jdn
   write(*,*) 'Current JD with time:', jd
   
   ! Julian Day calculations are useful for:
   ! - Astronomical computations
   ! - Converting between different calendar systems  
   ! - High-precision time interval calculations
   ! - Satellite orbital mechanics
   
end program julian_day_examples
```

**Julian Day Notes:**
- Julian Day Number (JDN) is the integer count of days since January 1, 4713 BCE (proleptic Julian calendar)
- Julian Day (JD) includes fractional time, where JD 0.0 begins at noon UTC
- Commonly used in astronomy for precise time calculations
- The library uses overflow-safe 64-bit integer arithmetic for reliability

## Format Specifiers

The library supports typical format specifiers for date formatting:

| Specifier | Description | Example |
|-----------|-------------|---------|
| `%Y` | 4-digit year | 2024 |
| `%m` | Month (01-12) | 03 |
| `%d` | Day of month (01-31) | 14 |
| `%H` | Hour (00-23) | 09 |
| `%M` | Minute (00-59) | 26 |
| `%S` | Second (00-59) | 53 |
| `%B` | Full month name | March |
| `%A` | Full weekday name | Thursday |

## Error Handling

The library handles errors gracefully:
- Invalid date parameters return special error timestamps
- Parsing failures are indicated by returned error values

## Thread Safety

The library is thread-safe for read operations. For write operations in multi-threaded environments, appropriate synchronization should be implemented by the calling application.

## Contributing

Contributions are welcome! Please ensure that:
- All new features include appropriate test cases
- Code follows the existing style conventions
- Documentation is updated for new functionality

## License

This project is licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.

## Support

For questions, bug reports, or feature requests, please create an issue in the project repository.
