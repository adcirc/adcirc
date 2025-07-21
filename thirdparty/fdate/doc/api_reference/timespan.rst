============
TimeSpan API
============

Overview
========

The ``t_timespan`` type represents a duration or time interval with millisecond precision. It can be used to express periods such as "3 days, 5 hours, and 30 minutes" and to perform duration-based calculations.

Constructor
===========

.. code-block:: fortran

   type(t_timespan) :: ts
   
   ! Default constructor (zero duration)
   ts = t_timespan()
   
   ! From components (all parameters are optional)
   ts = t_timespan(days, hours, minutes, seconds, milliseconds)

Parameters
----------

+------------------+----------+----------+---------------------------------------+
| Parameter        | Type     | Optional | Description                           |
+==================+==========+==========+=======================================+
| ``days``         | integer  | Yes      | Number of days                        |
+------------------+----------+----------+---------------------------------------+
| ``hours``        | integer  | Yes      | Number of hours                       |
+------------------+----------+----------+---------------------------------------+
| ``minutes``      | integer  | Yes      | Number of minutes                     |
+------------------+----------+----------+---------------------------------------+
| ``seconds``      | integer  | Yes      | Number of seconds                     |
+------------------+----------+----------+---------------------------------------+
| ``milliseconds`` | integer  | Yes      | Number of milliseconds                |
+------------------+----------+----------+---------------------------------------+

Component Methods
=================

These methods return the individual components of the TimeSpan:

.. code-block:: fortran

   integer :: d, h, m, s, ms
   type(t_timespan) :: ts
   
   ! Get components
   d = ts%days()        ! Days component (not total days)
   h = ts%hours()       ! Hours component (0-23)
   m = ts%minutes()     ! Minutes component (0-59)
   s = ts%seconds()     ! Seconds component (0-59)
   ms = ts%milliseconds() ! Milliseconds component (0-999)

Total Methods
=============

These methods return the total duration in various units:

.. code-block:: fortran

   integer(kind=8) :: d, h, m, s, ms
   type(t_timespan) :: ts
   
   ! Get total values
   d = ts%total_days()           ! Total number of days
   h = ts%total_hours()          ! Total number of hours
   m = ts%total_minutes()        ! Total number of minutes
   s = ts%total_seconds()        ! Total number of seconds
   ms = ts%total_milliseconds()  ! Total number of milliseconds

String Conversion
=================

.. code-block:: fortran

   character(len=64) :: str
   type(t_timespan) :: ts
   
   ! Convert to string
   str = ts%to_string()

The string format is ``[Nd ]HH:MM:SS[.mmm]`` where parts in brackets are optional:
 * ``N`` days are shown only if non-zero
 * Milliseconds (``.mmm``) are shown only if non-zero

Arithmetic Operations
=====================

TimeSpan supports basic arithmetic operations:

.. code-block:: fortran

   type(t_timespan) :: ts1, ts2, result
   
   ! Addition
   result = ts1 + ts2
   
   ! Subtraction
   result = ts1 - ts2
   
   ! Multiplication by a scalar
   result = ts1 * 2
   
   ! Division by a scalar
   result = ts1 / 2

Comparison Operations
=====================

TimeSpan supports all standard comparison operators:

.. code-block:: fortran

   type(t_timespan) :: ts1, ts2
   logical :: result
   
   ! Equality
   result = ts1 == ts2
   
   ! Inequality
   result = ts1 /= ts2
   
   ! Less than
   result = ts1 < ts2
   
   ! Greater than
   result = ts1 > ts2
   
   ! Less than or equal
   result = ts1 <= ts2
   
   ! Greater than or equal
   result = ts1 >= ts2

Examples
========

Basic Usage
-----------

.. code-block:: fortran

   type(t_timespan) :: duration
   
   ! Create a timespan of 2 days, 5 hours, 30 minutes
   duration = t_timespan(days=2, hours=5, minutes=30)
   
   ! Print the formatted string
   print *, "Duration: ", duration%to_string()  ! "2d 05:30:00"
   
   ! Get hours component (not total hours)
   print *, "Hours component: ", duration%hours()  ! 5
   
   ! Get total hours
   print *, "Total hours: ", duration%total_hours()  ! 53

Arithmetic
----------

.. code-block:: fortran

   type(t_timespan) :: span1, span2, result
   
   ! Define two timespans
   span1 = t_timespan(hours=5)
   span2 = t_timespan(hours=3)
   
   ! Add them
   result = span1 + span2
   print *, "Sum: ", result%total_hours()  ! 8
   
   ! Multiply by 2
   result = span1 * 2
   print *, "Double: ", result%total_hours()  ! 10
