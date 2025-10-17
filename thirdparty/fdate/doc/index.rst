=====
FDate
=====

**A modern, high-precision datetime library for Fortran applications**

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License: MIT

.. image:: https://img.shields.io/badge/C%2B%2B-17%2F20-blue.svg
   :alt: C++ 17/20

Overview
========

FortranDate is a high-performance date and time manipulation library designed specifically for scientific and engineering applications that require precise time handling. It provides an intuitive Fortran interface while leveraging the efficiency and robustness of C++ under the hood.

Key Features
============

* **Millisecond Precision**: All time calculations maintain millisecond precision
* **Modern Fortran Interface**: Clean, object-oriented design with operator overloading
* **Comprehensive Arithmetic**: Add, subtract, compare and manipulate date and time values
* **Format Flexibility**: Parse and format dates with customizable format strings
* **ISO 8601 Support**: Built-in support for ISO 8601 date format

Built for Numerical Modelers
============================

FortranDate was specifically designed with scientific and engineering applications in mind:

* **High Performance**: Optimized for computationally intensive applications
* **Consistency**: Reliable behavior across different systems
* **Precision**: Millisecond resolution for demanding scientific applications
* **Intuitive API**: Familiar, object-oriented interface for Fortran developers

Quick Example
=============

.. code-block:: fortran

   program date_example
      use mod_datetime
      implicit none
      
      type(t_datetime) :: start_time, end_time
      type(t_timespan) :: elapsed
      
      ! Create dates
      start_time = t_datetime(2025, 5, 14, 8, 30, 0)
      end_time = t_datetime(2025, 5, 14, 9, 45, 30)
      
      ! Calculate difference
      elapsed = end_time - start_time
      
      ! Display results
      print *, "Simulation ran for: ", elapsed%total_minutes(), " minutes"
      
   end program date_example

Contents
========

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   installation
   getting_started
   api_reference
   examples
   limitations
   contributing

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
