.. meta::
   :description: WTIMINC in ADCIRC
   :keywords: adcirc, wtiminc

.. _wtiminc_supplemental:

WTIMINC Supplemental
====================

**WTIMINC** is a parameter in the :ref:`fort.15 file <fort15>` that is
the time increment between meteorological forcing data sets (in seconds) in
input files like the :ref:`fort.22 <fort22>` file. This parameter and the :ref:`line on
which it appears :ref:<fort15>` depend on the value of
:ref:`NWS <NWS>`. For a detailed breakdown of what goes on this line in the
:ref:`fort.15 file <fort15>`, see the :ref:`supplemental meteorological/wave/ice
parameters <supplemental_meteorological_wave_ice_parameters>` page. Broadly,
``WTIMINC`` is required for ``ABS(NWS)>1`` unless you are using one of the
parametric vortex models for meteorology. Relatedly, note that when
``NWS=-14``, a second parameter, WTIMINC_12
defines the time increment for meteorological data for the OWI-formatted files.
