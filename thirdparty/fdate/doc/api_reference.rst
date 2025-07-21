=============
API Reference
=============

This section provides a comprehensive reference of the FortranDate API.

.. toctree::
   :maxdepth: 2
   
   api_reference/timespan
   api_reference/datetime

Overview
========

FortranDate consists of two main types and their associated procedures:

* :doc:`api_reference/timespan` - Represents a duration or time interval
* :doc:`api_reference/datetime` - Represents a specific point in time

Both types support operator overloading for natural syntax in Fortran code, including arithmetic operations and comparisons.

Module Structure
================

All functionality is exposed through the ``mod_datetime`` module. Import it at the beginning of your code:

.. code-block:: fortran

   use mod_datetime

Public and Private Components
=============================

The module exposes the following public entities:

Types:
   * ``t_timespan``
   * ``t_datetime``

Functions:
   * ``now()`` - Get the current time

Operators:
   * ``+`` - Addition
   * ``-`` - Subtraction
   * ``*`` - Multiplication (for TimeSpan)
   * ``/`` - Division (for TimeSpan)
   * ``==`` - Equality
   * ``/=`` - Inequality
   * ``<`` - Less than
   * ``>`` - Greater than
   * ``<=`` - Less than or equal
   * ``>=`` - Greater than or equal

Date and Time Formats
=====================

Both DateTime and TimeSpan types support string conversion with format specifiers.

DateTime Format Specifiers
--------------------------

Format specifiers follow the convention of ``strftime`` from C. Here are the most common:

+-------------+---------------------------------------+---------------+
| Specifier   | Description                           | Example       |
+=============+=======================================+===============+
| ``%Y``      | Year with century (0000-9999)         | 2025          |
+-------------+---------------------------------------+---------------+
| ``%m``      | Month as decimal number (01-12)       | 05            |
+-------------+---------------------------------------+---------------+
| ``%d``      | Day of month (01-31)                  | 16            |
+-------------+---------------------------------------+---------------+
| ``%H``      | Hour in 24-hour format (00-23)        | 14            |
+-------------+---------------------------------------+---------------+
| ``%M``      | Minute (00-59)                        | 30            |
+-------------+---------------------------------------+---------------+
| ``%S``      | Second (00-59)                        | 45            |
+-------------+---------------------------------------+---------------+
| ``%F``      | Short for %Y-%m-%d (ISO date format)  | 2025-05-16    |
+-------------+---------------------------------------+---------------+
| ``%T``      | Short for %H:%M:%S (ISO time format)  | 14:30:45      |
+-------------+---------------------------------------+---------------+
| ``%A``      | Full weekday name                     | Friday        |
+-------------+---------------------------------------+---------------+
| ``%B``      | Full month name                       | May           |
+-------------+---------------------------------------+---------------+
| ``%a``      | Abbreviated weekday name              | Fri           |
+-------------+---------------------------------------+---------------+
| ``%b``      | Abbreviated month name                | May           |
+-------------+---------------------------------------+---------------+
| ``%I``      | Hour in 12-hour format (01-12)        | 02            |
+-------------+---------------------------------------+---------------+
| ``%p``      | AM or PM                              | PM            |
+-------------+---------------------------------------+---------------+
| ``%Z``      | Time zone abbreviation                | UTC           |
+-------------+---------------------------------------+---------------+

TimeSpan String Format
----------------------

TimeSpans are formatted as:
``[Nd ]HH:MM:SS[.mmm]``

Where:
* ``N`` is the number of days (if non-zero)
* ``HH`` is hours (00-23)
* ``MM`` is minutes (00-59)
* ``SS`` is seconds (00-59)
* ``mmm`` is milliseconds (000-999, if non-zero)

Examples:
* ``2d 05:30:00`` - 2 days, 5 hours, 30 minutes
* ``01:15:30.500`` - 1 hour, 15 minutes, 30.5 seconds
