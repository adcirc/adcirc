Ramping Meteorological Forcing at Hotstart
==========================================

When running ADCIRC with meteorological forcing, it's generally best to ramp in forcing terms to avoid shocking the system. While parameters like ``DRAMP`` define the start of ramp time from cold-start, there is a special parameter to facilitate ramping of meteorological data at hotstart time.

Configuration
-------------

To enable meteorological ramping at hotstart, you need to:

1. Set ``NRAMP = 8`` in the fort.15 file
2. Supply ``DRAMPUnMete`` on the ``DRAMP`` line, which specifies the number of days at which to start meteorological ramping

Example
-------

Let's say you run a 15-day tide-only simulation and output a hotstart file. To apply a 0.5-day ramp to the meteorological forcing starting at run day 15.0 (the hotstart simulation's start time), your ``DRAMP`` line in fort.15 would look like:

.. code-block:: none

   5 0 0 0 0 0 0.5 0 15 ! DRAMP, DRAMPExtFlux, FluxSettlingTime, DRAMPIntFlux, DRAMPElev, DRAMPTip, DRAMPMete, DRAMPWRad, DRAMPUnMete

Parameter Description
---------------------

The parameters in order are:

- ``DRAMP``: General ramping parameter
- ``DRAMPExtFlux``: External flux ramping
- ``FluxSettlingTime``: Flux settling time
- ``DRAMPIntFlux``: Internal flux ramping
- ``DRAMPElev``: Elevation ramping
- ``DRAMPTip``: Tidal potential ramping
- ``DRAMPMete``: Meteorological forcing ramping (0.5 days in this example)
- ``DRAMPWRad``: Wave radiation stress ramping
- ``DRAMPUnMete``: Days at which to start meteorological ramping (15 days in this example)

Usage Notes
-----------

1. This feature is particularly useful when initializing meteorological forcing at hotstart time
2. The ramping helps prevent numerical instabilities that could arise from sudden introduction of meteorological forcing
3. The ramp duration (``DRAMPMete``) and start time (``DRAMPUnMete``) should be chosen based on your specific simulation needs
4. Make sure the ``DRAMPUnMete`` value matches your hotstart simulation's start time 