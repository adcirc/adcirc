.. meta::
   :description: Idealized Channel Problem in ADCIRC
   :keywords: adcirc, idealized channel problem

Idealized Channel Problem
=========================

This example tests ADCIRC version 55 (and beyond). It tests the simulation of a
diurnal tide on a sloping beach with a channel along its centerline (adapted
from [1]_). It tests lateral periodic boundary conditions and the
absorption-generation sponge layer [2]_ [3]_. The test finishes in about 8
minutes in parallel ADCIRC (2 processors) for 6 hours of simulation. Note that
the short 6 hour length of the test is chosen only to limit simulation time for
the `GitHub test
suite <https://github.com/adcirc/adcirc-cg-testsuite/tree/v55/adcirc/adcirc_ideal_channel-2d-parallel>`__
where the test case been found. Users may extend the simulation length to
simulate more of the inundating phase of the incoming wave.

Mesh
----

The mesh is comprised of 64,415 vertices and 127,784 triangular elements, with
resolution in the 10-60 m range. The mesh is symmetrical in the east-west
direction so that the east and west lateral boundary vertices match for the
application of the periodic lateral boundary conditions. An elevation specified
boundary condition and absorption-generation sponge layer is prescribed at the
southern end of the domain.

.. figure:: /_static/images/user_guide/examples/idealized_channel_problem/IdealChannel.png
   :width: 1000px

   Left: Mesh triangulation and resolution. Blue line shows the elevation specified boundary 
   condition location, green and yellow lines on the sides show the periodic lateral boundary 
   condition locations. Center: Mesh topo-bathy. Right: The sponge strength coefficients.

.. list-table::
   :class: wrap-table tighter-table
   :widths: 50 50

   * - .. figure:: /_static/images/user_guide/examples/idealized_channel_problem/Channel_Elev.gif

          Elevation time series for the idealized channel problem

     - .. figure:: /_static/images/user_guide/examples/idealized_channel_problem/Channel_Vel.gif

          North-south velocity time series for the idealized channel problem

Options/Features Tested
-----------------------

-  :ref:`IM <IM>` = 111112: Uses the explicit scheme (computational time step
   is 2 seconds).
-  :ref:`A00 <A00>`, :ref:`B00 <B00>`, :ref:`C00 <C00>` = 0.0, 1.0, 0.0:
   Must be used with explicit scheme.
-  :ref:`NOUTGE <NOUTGE>` = 5: Outputs the global elevations into a netCDF4
   :ref:`fort.63 file <fort63>`.
-  :ref:`NOUTGV <NOUTGV>` = 5: Outputs the global velocities into a netCDF4
   :ref:`fort.64 file <fort64>`.
-  :ref:`NOUTGW <NOUTGW>` = 5: Outputs the global meteorology into a netCDF4
   :ref:`fort.73 file <fort73>` (pressure) and a netCDF4 :ref:`fort.74
   file <fort74>` (velocity).
-  :ref:`sponge_generator_layer <sponge_generator_layer>`:
   Applies a sponge layer to absorb outgoing waves while generating incoming
   waves. In this case incoming diurnal tidal waves are generated using the
   ``fort.53001`` and ``fort.54001`` input files.
   :ref:`OceanMesh2D <oceanmesh2d>` functions can be
   used to automatically generate the sponge_generator_layer attribute
   (`Calc_Sponge <https://github.com/CHLNDDEV/OceanMesh2D/blob/Projection/utilities/Calc_Sponge.m>`__)
   and the input files
   (`Make_f5354 <https://github.com/CHLNDDEV/OceanMesh2D/blob/Projection/utilities/Make_f5354.m>`__).
-  :ref:`IBTYPE=94 <ibtype>`: Node pairs are matched along opposite
   lateral boundaries where a periodic (repeating) boundary condition is
   applied.

References
----------

.. raw:: html

   <references />

.. [1]
   Roberts, K.J., Dietrich, J.C., Wirasaet, D., Pringle, W.J., Westerink, J.J.,
   2021. Dynamic Load Balancing for Predictions of Storm Surge and Coastal
   Flooding. Environmental Modelling and Software, 105045.
   https://doi.org/10.1016/j.envsoft.2021.105045

.. [2]
   Pringle, W.J., Wirasaet, D., Suhardjo, A., Meixner, J., Westerink, J.J.,
   Kennedy, A.B., Nong, S., 2018. Finite-Element Barotropic Model for the Indian
   and Western Pacific Oceans: Tidal Model-Data Comparisons and Sensitivities.
   Ocean Model. 129, 13–38. https://doi.org/10.1016/j.ocemod.2018.07.003

.. [3]
   Pringle, W.J., Gonzalez-lopez, J., Joyce, B., Westerink, J.J., van der
   Westhuysen, A.J., 2019. Baroclinic Coupling Improves Depth-Integrated
   Modeling of Coastal Sea Level Variations around Puerto Rico and the U.S.
   Virgin Islands. J. Geophys. Res. Ocean. 124, 2196–2217.
   https://doi.org/10.1029/2018JC014682
