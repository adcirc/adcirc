.. meta::
   :description: Global Storm Tide - Hurricane Katrina in ADCIRC
   :keywords: adcirc, global storm tide - hurricane katrina

Global Storm Tide - Hurricane Katrina
=====================================

This example tests ADCIRC version 55 (and beyond). It tests the simulation of
the storm tides on the spherical Earth under astronomical and atmospheric
forcing during August 2005 when Hurricane Katrina impacted the Gulf of Mexico.
The results of interest are the global elevations, velocities and meteorology.
The test finishes in about 5 minutes in serial ADCIRC for a full month of
simulation. Find the test at the `GitHub test
suite <https://github.com/adcirc/adcirc-cg-testsuite/tree/v55/adcirc/adcirc_global-tide%2Bsurge-2d>`__.

Mesh
----

The mesh is a coarse representation of the spherical Earth with minimum
resolution of approximately 50 km, comprised of 27,330 vertices and 50,859
triangular elements.

.. _optionsfeatures_tested:

Options/Features Tested
-----------------------

-  :ref:`ICS <ICS>` = -22: Uses the Mercator projection with a coordinate
   rotation to remove the pole singularity (need to provide a
   :ref:`fort.rotm <fortrotm>`).
-  :ref:`IM <IM>` = 513113: Uses the fully implicit scheme for the gravity wave
   term (computational time step is 12 minutes).
-  :ref:`NTIP <NTIP>` = 2: Equilibrium tide + self-attraction and loading tide
   (read from a :ref:`fort.24 file <fort24>` forcing for 10 tidal
   constituents.
-  :ref:`NWS <NWS>` = -14: Reads from GRIB2 files that specify the global
   atmospheric forcing (6-hourly CFS reanalysis data) in addition to OWI ASCII
   files that specify the 3-hourly atmospheric forcing in the Hurricane Katrina
   landfall region.
-  :ref:`WTIMINC <WTIMINC>` = 21600, 10800: First value gives the temporal
   interval of the GRIB2 met data (6 hours), second value gives the temporal
   interval of the OWI met data (3 hours) in seconds.
-  :ref:`A00 <A00>`, :ref:`B00 <B00>`, :ref:`C00 <C00>` = 0.5, 0.5, 0:
   Ensures that the fully implicit scheme is stable with a large time step.
-  :ref:`ESLM <ESLM>` = -0.2: Enables the Smagorinsky turbulence closure with a
   coefficient of 0.2.
-  :ref:`NOUTGE <NOUTGE>` = 5: Outputs the global elevations into a netCDF4
   :ref:`fort.63 file <fort63>`.
-  :ref:`NOUTGV <NOUTGV>` = 5: Outputs the global velocities into a netCDF4
   :ref:`fort.64 file <fort64>`.
-  :ref:`NOUTGW <NOUTGW>` = 5: Outputs the global meteorology into a netCDF4
   :ref:`fort.73 file <fort73>` (pressure) and a netCDF4 :ref:`fort.74
   file <fort74>` (velocity).
-  :ref:`internal_tide_friction <internal_tide_friction>`:
   Spatially varying linear wave drag :ref:`fort.13 file <fort13>` attribute
   accounting for energy conversion due to internal tide generation in the deep
   ocean.
-  :ref:`quadratic_friction_coefficient_at_sea_floor <quadratic_friction_coefficient_at_sea_floor>`:
   Spatially varying quadratic bottom friction :ref:`fort.13 file <fort13>`
   attribute.

