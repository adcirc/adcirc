.. meta::
   :description: ADCIRC Tidal Parameters
   :keywords: adcirc, tides, tidal parameters, tidal forcing

.. _astronomical_tides:

Astronomical Tides
==================

ADCIRC supports two primary methods for implementing tidal forcings in your model: tidal potential terms and open boundary conditions. These can be used independently or in combination to achieve accurate tidal simulations.

Tidal Potential Terms
---------------------

Tidal potential terms represent the direct gravitational effects of celestial bodies (primarily the moon and sun) on the water surface. In ADCIRC, these are specified in the :ref:`fort.15 <fort15>` file through:

* ``NTIP``: Number of tidal potential constituents to be included
* ``TIPOTAG``: Name/tag for each tidal potential constituent
* ``TPK``, ``AMIGT``, ``ETRF``, ``FFT``, ``FACET``: Parameters defining amplitude, frequency, earth tide reduction factor, and nodal factor for each constituent

Open Boundary Conditions
------------------------

Open boundary conditions allow you to specify water surface elevations at the model domain boundaries. These can be implemented in two ways:

1. **Periodic (Harmonic) Conditions**:
   Specified in the :ref:`fort.15 <fort15>` file using:
   
   * ``NBFR``: Number of forcing frequencies on open boundaries
   * ``BOUNTAG``: Name/tag for each boundary forcing constituent
   * ``AMIG``, ``FF``, ``FACE``: Parameters for amplitude, nodal factor, and equilibrium argument
   * ``EMO``, ``EFA``: Elevation amplitude and phase for each constituent at boundary nodes

2. **Non-periodic Conditions**:
   Time series of water surface elevations are defined in the :ref:`fort.19 <fort19>` file.

The locations of open boundaries are defined in the :ref:`fort.14 <fort14>` file through the open boundary node strings.

Internal Tide Energy Conversion
-------------------------------

For large-scale ocean models, internal tide energy conversion should be considered. This represents the conversion of barotropic to baroclinic energy as surface tides flow over deep ocean topography. It can be implemented through:

* A spatially varying nodal attribute called ``internal_tide_friction`` in the fort.13 file
* Typically applied only in deep ocean regions (>100-500m depth)
* Accounts for approximately 30% of global barotropic tidal dissipation

See :ref:`internal_tide_energy_conversion` for more details.

.. toctree::
   :maxdepth: 1
   :caption: Contents:
   :hidden:
   
   internal_tide_energy_conversion

Best Practices
--------------

1. For coastal models, focus on accurate open boundary conditions
2. For large-scale ocean models, include both tidal potential terms and internal tide energy conversion
3. Ensure boundary conditions are properly ramped up using the ``NRAMP`` parameter in fort.15
4. Validate your tidal implementation against known tidal constituents or observations

For detailed implementation guidance, refer to the :ref:`fort.14 <fort14>` and :ref:`fort.15 <fort15>` documentation.

See also:

.. toctree::
   :maxdepth: 1

   ali_dispersion_conrtol
