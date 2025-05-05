Fort.71: Atmospheric Pressure Time Series at Specified Meteorological Recording Stations
========================================================================================

The fort.71 file contains atmospheric pressure time series data at specified meteorological recording stations as defined in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`. This file is generated when meteorological output is enabled.

File Structure
--------------

The basic file structure is shown below. Each line of output is represented by a line containing the output variable name(s). Loops indicate multiple lines of output.

.. parsed-literal::

   :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`
   :ref:`NTRSPM <NTRSPM>`, :ref:`NSTAM <NSTAM>`, :ref:`DTDP <DTDP>`\*:ref:`NSPOOLM <NSPOOLM>`, :ref:`NSPOOLM <NSPOOLM>`, :ref:`IRTYPE <IRTYPE>`
   :ref:`TIME <TIME>`, :ref:`IT <IT>`
   for k=1, :ref:`NSTAM <NSTAM>`
      k, :ref:`RMP00(k) <RMP00>`
   end k loop

   Repeat the block above for each time step until the end of the simulation.

Notes
-----

* Output may be in ascii or binary format depending on how :ref:`NOUTM <NOUTM>` is set in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`
* If binary output is specified, the station number (k) is not included in the output
* The meteorological recording station locations are specified in the fort.15 file 