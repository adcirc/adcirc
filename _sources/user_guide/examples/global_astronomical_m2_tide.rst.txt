.. meta::
   :description: Global Astronomical M2 Tide in ADCIRC
   :keywords: adcirc, global astronomical m2 tide

Global Astronomical M2 Tide
===========================

This example tests ADCIRC version 55 (and beyond). It tests the simulation of
the astronomical M2 tidal constituent on the spherical Earth under equilibrium
tidal forcing with the inclusion of the self-attraction and loading tide. The
results of interest are the M2 tidal constituent amplitudes and phases of
elevations and velocities from the least-squares harmonic analysis. The test
finishes in about 5 minutes in serial ADCIRC for a full month of simulation.
Find the test at the `GitHub test
suite <https://github.com/adcirc/adcirc-cg-testsuite/tree/v55/adcirc/adcirc_global-tide-2d>`__.

Mesh
----

The mesh is a coarse representation of the spherical Earth with minimum
resolution of approximately 50 km, comprised of 27,330 vertices and 50,859
triangular elements.

Options/Features Tested
-----------------------

-  :ref:`ICS <ICS>` = -22: Uses the Mercator projection with a coordinate
   rotation to remove the pole singularity (need to provide a
   :ref:`fort.rotm <fortrotm>`).
-  :ref:`IM <IM>` = 513113: Uses the fully implicit scheme for the gravity wave
   term (computational time step is 12 minutes).
-  :ref:`NTIP <NTIP>` = 2: equilibrium tide + self-attraction and loading tide
   forcing (read from a :ref:`fort.24 file <fort24>`).
-  :ref:`A00 <A00>`, :ref:`B00 <B00>`, :ref:`C00 <C00>` = 0.5, 0.5, 0:
   Ensures that the fully implicit scheme is stable with a large time step.
-  :ref:`ESLM <ESLM>` = -0.2: enables the Smagorinsky turbulence closure with a
   coefficient of 0.2.
-  :ref:`NHAGE <NHAGE>` = 5: outputs the harmonic constituent elevations into a
   netCDF4 :ref:`fort.53 file <fort53>`.
-  :ref:`NHAGV <NHAGV>` = 5: outputs the harmonic constituent velocities into a
   netCDF4 :ref:`fort.54 file <fort54>`.
-  :ref:`internal_tide_friction <internal_tide_friction>`:
   spatially varying linear wave drag :ref:`fort.13 file <fort13>` attribute
   accounting for energy conversion due to internal tide generation in the deep
   ocean.
-  :ref:`quadratic_friction_coefficient <quadratic_friction_coefficient_at_sea_floor>`:
   spatially varying quadratic bottom friction :ref:`fort.13 file <fort13>`
   attribute.
