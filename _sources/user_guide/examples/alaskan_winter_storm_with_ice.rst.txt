.. meta::
   :description: Alaskan Winter Storm with Ice in ADCIRC
   :keywords: adcirc, alaskan winter storm with ice

Alaskan Winter Storm with Ice
=============================

This example tests ADCIRC version 55 (and beyond). It tests the simulation of
the storm tides in a regional Alaska domain under astronomical and atmospheric
forcing in November 2011 during a strong winter storm in the presence of sea ice
(affecting the surface wind drag) [1]_. The results of interest are the global
elevations, velocities and meteorology. The test finishes in about 5 minutes in
serial ADCIRC for two weeks of simulation. Find the test at the `GitHub test
suite <https://github.com/adcirc/adcirc-cg-testsuite/tree/v55/adcirc/adcirc_alaska_ice-2d>`__.

Mesh
----

The mesh was generated using the OceanMesh2D Alaska
`Example_8_AK.m <https://github.com/CHLNDDEV/OceanMesh2D/blob/Projection/Examples/Example_8_AK.m>`__.
The domain encompasses the Gulf of Alaska, Bering Sea, and Chukchi Sea with a
minimum resolution of 5 km, comprised of 15,876 vertices and 27,757 triangular
elements.

Options/Features Tested
-----------------------

-  :ref:`ICS <ICS>` = 20: Equal-Area cylindrical projection.
-  :ref:`IM <IM>` = 513111: Uses the implicit scheme for the linear component
   of the gravity wave term (computational time step is 4 minutes).
-  :ref:`NTIP <NTIP>` = 2: Equilibrium tide + self-attraction and loading tide
   (read from a :ref:`fort.24 file <fort24>` forcing for 8 tidal
   constituents.
-  :ref:`NWS <NWS>` = 14014: Reads from GRIB2 files that specify the global
   atmospheric forcing and sea-ice concentration (6-hourly CFSv2 reanalysis
   data). Sea-ice concentration affects the wind drag coefficient [1]_.
-  :ref:`WTIMINC <WTIMINC>` = 21600, 21600: First value gives the temporal
   interval of the GRIB2 met data (6 hours), second value gives the temporal
   interval of the GRIB2 ice data (6 hours) - these should always be the same.
-  :ref:`A00 <A00>`, :ref:`B00 <B00>`, :ref:`C00 <C00>` = 0.4, 0.4, 0.2:
   Ensures that the implicit scheme is stable with a fairly large time step.
-  :ref:`ESLM <ESLM>` = -0.2: Enables the Smagorinsky turbulence closure with a
   coefficient of 0.2.
-  :ref:`NOUTGE <NOUTGE>` = 5: Outputs the global elevations into a netCDF4
   :ref:`fort.63 file <fort63>`.
-  :ref:`NOUTGV <NOUTGV>` = 5: Outputs the global velocities into a netCDF4
   :ref:`fort.64 file <fort64>`.
-  :ref:`NOUTGW <NOUTGW>` = 5: Outputs the global meteorology into a netCDF4
   :ref:`fort.73 file <fort73>` (pressure) and a netCDF4 :ref:`fort.74
   file <fort74>` (velocity).
-  :ref:`internal_tide_friction <internal_tide_friction>`: Spatially varying linear wave drag :ref:`fort.13 file <fort13>` attribute
   accounting for energy conversion due to internal tide generation in the deep
   ocean.
-  :ref:`&WarnElevControl namelist <fort15>`: Set
   "WarnElev", the warning elevation level, to 30-m (elevations reach beyond
   20-m [default] but remain below 30-m).
-  `&metControl namelist <Fort.15_file_format#Namelists>`__: Set "rhoAir", to
   1.29193 (density of air at 0 deg C for 1013 mbar); set "WindDragLimit" equal
   to 0.0025; set "invertedBarometerOnElevationBoundary" to true (in Alaska
   extremely large-scale low pressure systems persist and cross over the open
   boundaries, so it is important to have the inverted barometer condition along
   the elevation specified boundary); set "outputWindDrag" to true.

References
----------

.. raw:: html

   <references />

.. [1]
   Joyce, B.R., Pringle, W.J., Wirasaet, D., Westerink, J.J., Van der
   Westhuysen, A.J., Grumbine, R., Feyen, J., 2019. High resolution modeling of
   western Alaskan tides and storm surge under varying sea ice conditions. Ocean
   Model. 141, 101421. doi:10.1016/j.ocemod.2019.101421
