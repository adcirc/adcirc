============
DateTime API
============

Overview
========

The ``t_datetime`` type represents a specific point in time with millisecond precision. It supports date/time creation, parsing, formatting, arithmetic operations, and comparisons.

Constructors
============

.. code-block:: fortran

   type(t_datetime) :: dt
   
   ! Default constructor (epoch time - Jan 1, 1970 00:00:00 UTC)
   dt = t_datetime()
   
   ! From year, month, day
   dt = t_datetime(year, month, day)
   
   ! From year, month, day, hour, minute, second
   dt = t_datetime(year, month, day, hour, minute, second)
   
   ! From all components
   dt = t_datetime(year, month, day, hour, minute, second, millisecond)
   
   ! From timestamp (milliseconds since epoch)
   dt = t_datetime(timestamp)
   
   ! From string parsing
   dt = t_datetime_parse(str, format)

Parameters
----------

+-----------------+------------------+----------+---------------------------------------+
| Parameter       | Type             | Optional | Description                           |
+=================+==================+==========+=======================================+
| ``year``        | integer          | No       | Year (e.g., 2025)                     |
+-----------------+------------------+----------+---------------------------------------+
| ``month``       | integer          | No       | Month (1-12)                          |
+-----------------+------------------+----------+---------------------------------------+
| ``day``         | integer          | No       | Day (1-31)                            |
+-----------------+------------------+----------+---------------------------------------+
| ``hour``        | integer          | Yes      | Hour (0-23)                           |
+-----------------+------------------+----------+---------------------------------------+
| ``minute``      | integer          | Yes      | Minute (0-59)                         |
+-----------------+------------------+----------+---------------------------------------+
| ``second``      | integer          | Yes      | Second (0-59)                         |
+-----------------+------------------+----------+---------------------------------------+
| ``millisecond`` | integer          | Yes      | Millisecond (0-999)                   |
+-----------------+------------------+----------+---------------------------------------+
| ``timestamp``   | integer(kind=8)  | No       | Milliseconds since epoch              |
+-----------------+------------------+----------+---------------------------------------+
| ``str``         | character(len=*) | No       | Date/time string to parse             |
+-----------------+------------------+----------+---------------------------------------+
| ``format``      | character(len=*) | No       | Format string (strftime-compatible)   |
+-----------------+------------------+----------+---------------------------------------+

Current Time
============

To get the current date and time:

.. code-block:: fortran

   type(t_datetime) :: current_time
   
   ! Get current time
   current_time = now()

Component Methods
=================

These methods return the individual components of the DateTime:

.. code-block:: fortran

   type(t_datetime) :: dt
   
   ! Get components
   integer :: y, mo, d, h, mi, s, ms
   
   y = dt%year()         ! Year
   mo = dt%month()       ! Month (1-12)
   d = dt%day()          ! Day (1-31)
   h = dt%hour()         ! Hour (0-23)
   mi = dt%minute()      ! Minute (0-59)
   s = dt%second()       ! Second (0-59)
   ms = dt%millisecond() ! Millisecond (0-999)
   
   ! Get timestamp (milliseconds since epoch)
   integer(kind=8) :: ts
   ts = dt%timestamp()

Formatting Methods
==================

.. code-block:: fortran

   type(t_datetime) :: dt
   character(len=64) :: str
   
   ! Format with custom format string (strftime-compatible)
   str = dt%format(date_format, show_milliseconds)
   
   ! Convert to ISO 8601 string
   str = dt%to_iso_string(msec)

Parameters
----------

+-----------------------+-------------------+----------+-----------------------------------------------+
| Parameter             | Type              | Optional | Description                                   |
+=======================+===================+==========+===============================================+
| ``date_format``       | character(len=*)  | No       | Format string (strftime-compatible)           |
+-----------------------+-------------------+----------+-----------------------------------------------+
| ``show_milliseconds`` | logical           | Yes      | Whether to include milliseconds               |
+-----------------------+-------------------+----------+-----------------------------------------------+
| ``msec``              | logical           | Yes      | Whether to include milliseconds in ISO string |
+-----------------------+-------------------+----------+-----------------------------------------------+

Arithmetic Operations
=====================

DateTime supports arithmetic with TimeSpan objects:

.. code-block:: fortran

   type(t_datetime) :: dt1, dt2, result_date
   type(t_timespan) :: diff_ts
   
   ! Add a TimeSpan to a DateTime
   result_date = dt1 + ts
   
   ! Subtract a TimeSpan from a DateTime
   result_date = dt1 - ts
   
   ! Calculate the difference between two DateTimes
   diff_ts = dt1 - dt2

Comparison Operations
=====================

DateTime supports all standard comparison operators:

.. code-block:: fortran

   type(t_datetime) :: dt1, dt2
   logical :: is_equal, is_not_equal
   logical :: is_earlier, is_later
   logical :: is_earlier_or_equal,  is_later_or_equal
   
   ! Equality
   is_equal = dt1 == dt2
   
   ! Inequality
   is_not_equal = dt1 /= dt2
   
   ! Less than (earlier)
   is_earlier = dt1 < dt2
   
   ! Greater than (later)
   is_later = dt1 > dt2
   
   ! Less than or equal
   is_earlier_or_equal = dt1 <= dt2
   
   ! Greater than or equal
   is_later_or_equal = dt1 >= dt2

Examples
========

Creating and Formatting
-----------------------

.. code-block:: fortran

   program datetime_example
      use mod_datetime
      implicit none
      
      type(t_datetime) :: dt
      character(len=64) :: formatted
      
      ! Create a specific date and time
      dt = t_datetime(2025, 5, 16, 14, 30, 45)
      
      ! Format in different ways
      formatted = dt%format("%Y-%m-%d %H:%M:%S")
      print *, "Standard format: ", trim(formatted)  ! "2025-05-16 14:30:45"
      
      formatted = dt%format("%A, %B %d, %Y at %I:%M %p")
      print *, "Custom format: ", trim(formatted)    ! "Friday, May 16, 2025 at 02:30 PM"
      
      formatted = dt%to_iso_string()
      print *, "ISO format: ", trim(formatted)       ! "2025-05-16T14:30:45"
      
   end program datetime_example

Date Calculations
-----------------

.. code-block:: fortran

   program date_calc
      use mod_datetime
      implicit none
      
      type(t_datetime) :: start_date, end_date
      type(t_timespan) :: duration
      
      ! Calculate project timeline
      start_date = t_datetime(2025, 5, 1)
      
      ! Add 45 days to the start date
      duration = t_timespan(days=45)
      end_date = start_date + duration
      
      print *, "Project start: ", start_date%to_iso_string()
      print *, "Project end: ", end_date%to_iso_string()
      print *, "Project duration (days): ", duration%total_days()
      
   end program date_calc

Parsing Dates
-------------

.. code-block:: fortran

   program parse_date
      use mod_datetime
      implicit none
      
      type(t_datetime) :: dt
      character(len=32) :: date_str
      
      ! Parse a date string
      date_str = "2025-05-16 14:30:45"
      dt = t_datetime(date_str, "%Y-%m-%d %H:%M:%S")
      
      print *, "Parsed date components:"
      print *, "Year: ", dt%year()
      print *, "Month: ", dt%month()
      print *, "Day: ", dt%day()
      print *, "Hour: ", dt%hour()
      print *, "Minute: ", dt%minute()
      print *, "Second: ", dt%second()
      
   end program parse_date
