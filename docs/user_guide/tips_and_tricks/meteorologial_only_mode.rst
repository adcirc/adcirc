Meteorological-Only Mode
========================

ADCIRC can be run in "meteorological-only" mode, which is a quick and convenient way to read in and write out meteorological data files without performing hydrodynamic calculations. In this mode, only routines necessary for meteorological forcing are called - no water levels or velocities are computed.

This mode is achieved by turning off all ocean state-related output and disabling wave coupling. Since only meteorological data is being output, you can set ``DT`` equal to the time interval at which you want meteorological data to be output, and then set the meteorological output intervals (``NSPOOLM`` and/or ``NSPOOLGW``) to ``1`` (i.e., output every time step).

Required Settings
-----------------

For typical 2D ADCIRC runs, set the following parameters in the fort.15 file:

.. code-block:: none

   NOUTE = 0    ! No elevation station output
   NOUTV = 0    ! No velocity station output  
   NOUTGE = 0   ! No global elevation output
   NOUTGV = 0   ! No global velocity output
   NHSTAR = 0   ! No hotstart output
   NHASE = 0    ! No harmonic analysis of elevations at stations
   NHASV = 0    ! No harmonic analysis of velocities at stations
   NHAGE = 0    ! No harmonic analysis of global elevations
   NHAGV = 0    ! No harmonic analysis of global velocities

Additional Settings
-------------------

For passive scalar transport (``IM = 10``), also set:

.. code-block:: none

   NOUTC = 0    ! No concentration station output
   NOUTGC = 0   ! No global concentration output

For barotropic or baroclinic 3D ADCIRC, also set:

.. code-block:: none

   I3DSD = 0    ! No 3D density station output
   I3DSV = 0    ! No 3D velocity station output
   I3DST = 0    ! No 3D temperature station output
   I3DGD = 0    ! No 3D global density output
   I3DGV = 0    ! No 3D global velocity output
   I3DGT = 0    ! No 3D global temperature output

Usage Notes
-----------

1. This mode is useful when you only need to process meteorological forcing data without running the full hydrodynamic model.
2. The simulation speed can be increased by setting ``DT`` equal to your desired meteorological output interval.
3. Set ``NSPOOLM`` and/or ``NSPOOLGW`` to ``1`` to output meteorological data at every time step.
4. All hydrodynamic calculations are skipped, significantly reducing computational time. 