Fort.16: General Diagnostic Output
==================================

Output file which echo prints information from:

* The Grid and Boundary Information File
* The Model Parameter and Periodic Boundary Condition File

Additionally, this file:

* Provides some processed information
* Prints out error messages from ADCIRC

Notes for fort.16 file:
-----------------------

The Geophysical Fluid Dynamics Laboratory (GFDL) file format was developed to support the output from atmospheric models. The files have the following characteristics:

* Each ASCII GFDL met file contains one or more nested grid data where the nested grids are allowed to change in time.
* Coarse grid data is not stored where finer nest data is given.

The file format of the ADCIRC fort.22 for GFDL is as follows:

Line 1 — ASCII Text Header
Line 2 — Windmultiplier Value (real number)
Line 3 — Maximum Extrapolation distance (m)
Lines 4-end — cycleTime rampValue filename

The file format for each of the actual GFDL files is as follows:

Line 1: Number of grid cells (f10.4) NCELLS
Lines 2-NCELLS+1: Have ten columns of data formatted as 10f10.4

1. u (m/sec)
2. v (m/sec)
3. Temperature (K)
4. mixing ratio(kg/kg)
5. storm accum precipitation (cm)
6. sea level pressure (hPa)
7. longitude (decimal deg)
8. latitude (decimal deg)
9. hurricane hour
10. nest number (this is not always present)

Example
-------

To illustrate the definitions and descriptions provided, a concrete example of an ADCIRC fort.22 file for GFDL is provided as follows:

.. code-block:: none

   ! 1st line is a comment line, max length 1024 characters
   1.0     ! 2nd line is a velocity magnitude multiplier
   10000.0 ! 3rd line: maximum extrapolation distance (m)
   0.0 0.0 "/home/jason/isaac/gfdl/isaac_gfdl_file1" ! time (hours), ramp mult, filename
   6.0 0.5 "/home/jason/isaac/gfdl/isaac_gfdl_file2"
   12.0 1.0 "/home/jason/isaac/gfdl/isaac_gfdl_file3"

Note: When including the path to the data files, if the full path includes forward slashes (as it would if ADCIRC is executing on Unix or Linux), be sure to surround the full path file name with double quotes as shown in the example above. This is required because Fortran treats a bare forward slash in an input file as an end-of-record character. 