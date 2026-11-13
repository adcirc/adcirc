Fort.62: Depth-averaged Velocity Time Series at Specified Velocity Recording Stations
=====================================================================================

The fort.62 file contains depth-averaged velocity time series data at specified velocity recording stations as defined in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`. This file is generated when time series output is enabled for velocity stations.

File Structure
--------------

The basic file structure is shown below. Each line of output is represented by a line containing the output variable name(s). Loops indicate multiple lines of output.

.. parsed-literal::

   :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`
   :ref:`NTRSPV <NTRSPV>`, :ref:`NSTAV <NSTAV>`, :ref:`DTDP <DTDP>`\*:ref:`NSPOOLV <NSPOOLV>`, :ref:`NSPOOLV <NSPOOLV>`, :ref:`IRTYPE <IRTYPE>`
   :ref:`TIME <TIME>`, :ref:`IT <IT>`
   for k=1, :ref:`NSTAV <NSTAV>`
      k, :ref:`UU2(k) <UU>`, :ref:`VV2(k) <VV>`
   end k loop

   Repeat the block above for each time step until the end of the simulation.

Notes
-----

* Output may be in ascii or binary format depending on how :ref:`NOUTV <NOUTV>` is set in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`
* If binary output is specified, the station number (k) is not included in the output
* The velocity recording station locations can be specified either directly in the fort.15 file or through a separate vel_stat.151 file when :ref:`NSTAV <NSTAV>` is set to a negative value 