Fort.24: Self Attraction/Earth Load Tide Forcing File
=====================================================

Self attraction/earth load tides are used to force ADCIRC when NTIP=2 in the :ref:`Model Parameter and Periodic Boundary Condition File <fort15>`.

The format of this file is identical to the "tea" harmonic format with no velocity information included. Entries are grouped by constituents and must be in the same order as the tidal potential terms listed in the :ref:`Model Parameter and Periodic Boundary Condition File <fort15>`. Phases must be in degrees. Amplitudes must be in units compatible with the units of gravity. These values are modified by the nodal factor and equilibrium argument provided for the tidal potential terms.

File Structure
--------------

The basic file structure is shown below. Each line of input data is represented by a line containing the input variable name(s). Loops indicate multiple lines of input.

.. parsed-literal::

   for k=1, :ref:`NTIF <NTIF>`
      Alpha line
      Constituent frequency
      1
      Constituent name (e.g., M2)
      for j=1, :ref:`NP <NP>`
         :ref:`JN <JN>`, :ref:`SALTAMP(k,JN) <SALTAMP>`, :ref:`SALTPHA(k,JN) <SALTPHA>`
      end j loop
   end k loop

Notes
-----

1. The first four lines (Alpha line, Constituent frequency, 1, Constituent name) for each constituent must be present in the file but they are skipped over during the ADCIRC read.

2. Entries must be grouped by constituents and must be in the same order as the tidal potential terms listed in the :ref:`Model Parameter and Periodic Boundary Condition File <fort15>`.

3. Phases must be in degrees.

4. Amplitudes must be in units compatible with the units of gravity.

5. These values are modified by the nodal factor and equilibrium argument provided for the tidal potential terms.

Example
-------

The following is a simple example of the beginning of a fort.24 file:

.. code-block:: none

   Alpha line
   0.0805114017
   1
   M2
   1    0.12345E+00  0.23456E+00
   2    0.34567E+00  0.45678E+00
   3    0.56789E+00  0.67890E+00
