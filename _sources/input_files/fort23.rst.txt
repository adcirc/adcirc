Fort.23: Wave Radiation Stress Forcing File
===========================================

The fort.23 file contains wave radiation stress data that can be used to drive ADCIRC either by itself or in combination with other forcing (including winds). This file is read when ABS(:ref:`NWS <NWS>`) >= 100 in the :ref:`Model Parameter and Periodic Boundary Condition File <fort15>`. The format is similar to the meteorological input file used when :ref:`NWS <NWS>` = -4 (i.e. the PBL hurricane model input format following a hot start).

File Structure
--------------

The basic file structure is shown below. Each line of input data is represented by a line containing the input variable name(s). Loops indicate multiple lines of input.

.. parsed-literal::

   :ref:`JN <JN>`, :ref:`RSX(JN) <RSX>`, :ref:`RSY(JN) <RSY>`

Notes
-----

1. At least two datasets must be present in the file to allow for time interpolation. If only one dataset is present, the run will terminate with an unexpected end-of-file error.

2. Radiation stresses are input directly to a subset of nodes in the ADCIRC grid (as specified by the node number :ref:`JN <JN>`).

3. If ADCIRC is cold started, the first set of radiation stress data corresponds to TIME=:ref:`STATIM <STATIM>`. If ADCIRC is hot started, the first set of radiation stress data corresponds to TIME=HOT START TIME. Additional sets of radiation stress data must be provided every :ref:`RSTIMINC <RSTIMINC>`, where :ref:`RSTIMINC <RSTIMINC>` is the radiation stress time interval and is specified in the :ref:`Model Parameter and Periodic Boundary Condition File <fort15>`. Radiation stresses are interpolated in time to the ADCIRC time step.

4. Each data line must have the format I8, 2E13.5. Data input lines are repeated for as many nodes as desired. A line containing the # symbol in column 2 indicates radiation stress data at the next time increment begins on the following line. At each new time, any node that is not specified in the input file is assumed to have zero wave radiation stress.

5. Wave radiation stress must be input in units of velocity squared (consistent with the units of gravity). Stress in these units is obtained by dividing stress in units of force/area by the reference density of water.

6. Data must be provided for the entire model run, otherwise the run will crash!

Example
-------

The following is a simple example of the beginning of a fort.23 file:

.. code-block:: none

   1    0.12345E+00  0.23456E+00
   2    0.34567E+00  0.45678E+00
   3    0.56789E+00  0.67890E+00
   # 
   1    0.12345E+01  0.23456E+01
   2    0.34567E+01  0.45678E+01
   3    0.56789E+01  0.67890E+01
