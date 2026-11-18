.. meta::
   :description: Boundary conditions in ADCIRC
   :keywords: adcirc, boundary conditions

.. _boundary_conditions:

Boundary Conditions
===================

**Lateral boundary conditions** allow one to constrain the physics along
boundaries in the model. They are similar to `initial
conditions <initial_conditions>`__, which specify the model state at start time.
As an example, at the walls of a bathtub, the flow normal to the walls is zero
since they are impermeable. For general information on the topic, see the
`Wikipedia article <https://en.wikipedia.org/wiki/Boundary_value_problem>`__.
Since ADCIRC solves for water elevations and velocities (fluxes), typically one
(occasionally both) of these two quantities are constrained at the lateral
boundaries.

Open Boundaries
---------------

Open boundaries are boundaries where the water surface elevation is specified. Locations of the open boundaries are specified in :ref:`fort.14 <fort14>` file. For periodic conditions, the amplitudes, frequencies, and phases are defined in the :ref:`fort.15 <fort15>` file. For non-periodic conditions, time series of the water surface elevations are defined in the :ref:`fort.19 <fort19>` file.

Flux Boundaries
---------------

Flux boundaries are boundaries where the normal flux is specified. Locations of the flux boundaries are specified in :ref:`fort.14 <fort14>` file. For periodic conditions, the amplitudes, frequencies, and phases are defined in the :ref:`fort.15 <fort15>` file. For non-periodic conditions, time series of the normal fluxes are defined in the :ref:`fort.20 <fort20>` file.

By default ADCIRC weakly satisfies the no-flux boundary condition at mesh
boundaries. See :ref:`flux specified boundaries <flux_specified_boundaries>` for
details on fine-grained specification of the flux boundary conditions if
required.

.. _automatic_specific_of_boundary_conditions:

Automatic Specific of Boundary Conditions
-----------------------------------------

The OceanMesh2D [1]_\ `(GitHub
site) <https://github.com/CHLNDDEV/OceanMesh2D>`__ mesh generation toolbox has
the ability to automatically apply the basic no-flux and open ocean elevation
boundary conditions for an ADCIRC mesh. See the **makens** (make node-string)
function using the 'auto' option. Geospatial shoreline data is required to
automatically detect whether a mesh boundary is located in the ocean (applies
open ocean elevation boundary condition) or is along the shoreline/on-land
(applies natural no-flux boundary condition).

The **makens** function also contains other options for manually specifying
other boundary condition types such as rivers and weirs.

References
----------

.. raw:: html

   <references />

.. [1]
   Roberts, K. J., Pringle, W. J., & Westerink, J. J. (2019). OceanMesh2D 1.0:
   MATLAB-based software for two-dimensional unstructured mesh generation in
   coastal ocean modeling. Geoscientific Model Development, 12, 1847–1868.
   https://doi.org/10.5194/gmd-12-1847-2019

