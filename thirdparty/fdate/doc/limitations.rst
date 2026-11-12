==========================
Limitations & Dependencies
==========================

This page outlines the technical limitations, dependencies, and considerations when using FortranDate.

Technical Limitations
=====================

Precision
---------

* FortranDate provides **millisecond precision**, which is sufficient for most scientific and engineering applications
* Sub-millisecond precision is not supported

Date Range
----------

* The library supports dates from January 1, 1970 (Unix epoch) through the year 2262
* Dates outside this range may lead to integer overflow in the internal timestamp representation
* Historical dates before 1970 are represented with negative timestamps

Timezone Handling
-----------------

* All date and time values are in UTC (Coordinated Universal Time)
* The library does not currently provide timezone conversion capabilities
* When working with local times, you'll need to handle timezone conversions manually

Thread Safety
-------------

* The underlying C++ implementation is thread-safe for separate objects
* However, the Fortran interface does not guarantee thread safety when multiple threads manipulate the same object

Leap Seconds
------------

* Leap seconds are not handled explicitly
* When leap seconds occur, they are effectively treated as part of the preceding second

Year 2038 Problem
-----------------

* Unlike some 32-bit time implementations, FortranDate uses 64-bit integers for timestamps
* This avoids the "Year 2038 problem" where 32-bit time values would overflow
* The library can represent dates well beyond 2038 without issues

Dependencies
============

Howard Hinnant's Date Library
-----------------------------

FortranDate relies on Howard Hinnant's date library, a high-quality, header-only C++ library for calendar and time zone operations. This dependency:

* Provides the underlying implementation of calendar calculations
* Is included as a submodule in the FortranDate repository
* Does not require separate installation

C++ Standard Library
--------------------

The C++ portion of FortranDate depends on the C++ standard library, specifically:

* ``<chrono>`` for time handling
* ``<string>`` for string manipulation
* ``<optional>`` for error handling in parsing operations

Compiler Requirements
=====================

C++ Compiler
------------

* C++17 or C++20 compliant compiler
* Tested with:
  * GCC 9.0 or newer
  * Intel oneAPI 2021 or newer
  * NVHPC 21.5 or newer

Fortran Compiler
----------------

* Fortran 2003 or newer for derived type support
* Fortran 2008 or newer preferred for better interoperability features
* Tested with:
  * GFortran 9.0 or newer
  * Intel ifx 2021 or newer
  * NVHPC nvfortran 21.5 or newer

Build System
------------

* CMake 3.20 or newer

Platform Support
================

FortranDate is designed to be cross-platform and has been tested on:

* Linux (various distributions)
* macOS
* Windows (with MSVC, MinGW, or Intel compilers)

Memory Considerations
=====================

FortranDate is designed to be lightweight:

* DateTime objects are 8 bytes (64 bits)
* TimeSpan objects are 8 bytes (64 bits)
* No dynamic memory allocation is used in the core library

Performance Considerations
==========================

FortranDate is designed for efficiency:

* Most operations are O(1) (constant time)
* String parsing and formatting are more expensive operations
* Use timestamp-based comparisons when possible for maximum performance

Known Issues
============

Date Parsing
------------

* Date parsing is forgiving but may accept invalid dates in some edge cases
* Always validate parsed dates in critical applications

Format String Compatibility
---------------------------

* Format strings follow the C `strftime` convention
* Not all format specifiers may be available on all platforms
* Stick to common format specifiers for maximum compatibility

Working with Legacy Codes
=========================

When integrating with legacy Fortran codes:

* Take care with implicit typing rules (use `implicit none`)
* Date formatting may require adjustment of string lengths
* For interoperability with codes using other date representations, conversion functions may be necessary

Future Enhancements
===================

The development team is considering these enhancements for future releases:

* Timezone support
* Higher precision (microsecond and nanosecond)
* Additional parsing formats
* Calendar-specific calculations (business days, holidays)
* Better interface with existing Fortran datetime libraries
