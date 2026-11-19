===============
Getting Started
===============

This guide provides a quick introduction to using FortranDate in your Fortran applications.

Basic Concepts
==============

FortranDate provides two main derived types:

* ``t_datetime`` - Represents a specific point in time with millisecond precision
* ``t_timespan`` - Represents a duration or time interval

Both types support natural arithmetic operations and comparisons through operator overloading.

Using the Module
================

To use FortranDate in your Fortran code, simply include the module:

.. code-block:: fortran

   use mod_datetime

Creating DateTime Objects
=========================

There are several ways to create a ``t_datetime`` object:

.. code-block:: fortran

   ! Current time
   type(t_datetime) :: now_time
   now_time = now()
   
   ! From year, month, day
   type(t_datetime) :: date_only
   date_only = t_datetime(2025, 5, 16)
   
   ! From year, month, day, hour, minute, second
   type(t_datetime) :: date_time
   date_time = t_datetime(2025, 5, 16, 14, 30, 0)
   
   ! With milliseconds
   type(t_datetime) :: precise_time
   precise_time = t_datetime(2025, 5, 16, 14, 30, 0, 500)
   
   ! From a timestamp (milliseconds since epoch)
   type(t_datetime) :: from_timestamp
   from_timestamp = t_datetime(1718558400000_8) ! May 16, 2025 00:00:00
   
   ! Parsing from a string
   type(t_datetime) :: parsed_date
   parsed_date = t_datetime_parse("2025-05-16 14:30:00", "%Y-%m-%d %H:%M:%S")

Accessing DateTime Components
=============================

You can extract individual components from a datetime:

.. code-block:: fortran

   type(t_datetime) :: date_time
   date_time = t_datetime(2025, 5, 16, 14, 30, 45)
   
   print *, "Year:   ", date_time%year()
   print *, "Month:  ", date_time%month()
   print *, "Day:    ", date_time%day()
   print *, "Hour:   ", date_time%hour()
   print *, "Minute: ", date_time%minute()
   print *, "Second: ", date_time%second()
   print *, "Millisecond: ", date_time%millisecond()

Creating TimeSpan Objects
=========================

TimeSpans can be created in several ways:

.. code-block:: fortran

   ! From individual components
   type(t_timespan) :: span1
   span1 = t_timespan(days=1, hours=6, minutes=30)
   
   ! From just one unit
   type(t_timespan) :: one_day, two_hours, thirty_min
   one_day = t_timespan(days=1)
   two_hours = t_timespan(hours=2)
   thirty_min = t_timespan(minutes=30)

DateTime Arithmetic
===================

You can perform various arithmetic operations:

.. code-block:: fortran

   type(t_datetime) :: start_time, end_time
   type(t_timespan) :: duration, extra_time
   
   start_time = t_datetime(2025, 5, 16, 8, 0, 0)
   
   ! Add a timespan to a datetime
   duration = t_timespan(hours=3, minutes=30)
   end_time = start_time + duration
   
   ! Subtract a timespan from a datetime
   start_time = end_time - duration
   
   ! Get duration between two datetimes
   duration = end_time - start_time
   
   ! Timespan arithmetic
   extra_time = duration * 2
   extra_time = duration / 2

Formatting Dates
================

Format dates to strings for output:

.. code-block:: fortran

   type(t_datetime) :: date_time
   character(len=64) :: formatted_date
   
   date_time = t_datetime(2025, 5, 16, 14, 30, 45)
   
   ! Standard format
   formatted_date = date_time%format("%Y-%m-%d %H:%M:%S")
   print *, "Standard format: ", trim(formatted_date)
   
   ! Custom format
   formatted_date = date_time%format("%A, %B %d, %Y at %I:%M %p")
   print *, "Custom format: ", trim(formatted_date)
   
   ! ISO 8601 format
   formatted_date = date_time%to_iso_string()
   print *, "ISO format: ", trim(formatted_date)

Comparing DateTimes
===================

Compare dates using standard operators:

.. code-block:: fortran

   type(t_datetime) :: date1, date2
   logical :: result
   
   date1 = t_datetime(2025, 5, 16)
   date2 = t_datetime(2025, 5, 17)
   
   result = date1 == date2  ! Equality
   result = date1 /= date2  ! Inequality
   result = date1 < date2   ! Less than
   result = date1 > date2   ! Greater than
   result = date1 <= date2  ! Less than or equal
   result = date1 >= date2  ! Greater than or equal

The same comparison operators work with TimeSpan objects as well.

Next Steps
==========

Now that you understand the basics, explore the :doc:`api_reference` for detailed documentation of all available functionality, or check out the :doc:`examples` section for more complex use cases.
