Fort.61: Elevation Time Series at Specified Elevation Recording Stations
========================================================================

Elevation time series output at the elevation recording stations as specified in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`. This file is generated when time series output is enabled for elevation stations.

File Structure
--------------

The basic file structure is shown below. Each line of output is represented by a line containing the output variable name(s). Loops indicate multiple lines of output.

.. parsed-literal::

   :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`
   :ref:`NTRSPE <NTRSPE>`, :ref:`NSTAE <NSTAE>`, :ref:`DTDP <DTDP>`\*:ref:`NSPOOLE <NSPOOLE>`, :ref:`NSPOOLE <NSPOOLE>`, :ref:`IRTYPE <IRTYPE>`
   :ref:`TIME <TIME>`, :ref:`IT <IT>`
   for k=1, :ref:`NSTAE <NSTAE>`
      k, :ref:`ET00(k) <ET00>`
   end k loop

   Repeat the block above for each time step until the end of the simulation.

Notes
-----

* Output may be in ascii or binary format depending on how :ref:`NOUTE <NOUTE>` is set in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`
* If binary output is specified, the station number (k) is not included in the output 