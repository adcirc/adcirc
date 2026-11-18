.. _fort73:

Fort.73: Atmospheric Pressure Time Series at All Nodes in the Model Grid
========================================================================

The fort.73 file contains atmospheric pressure time series data at all nodes in the model grid as defined in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`. This file is generated when meteorological output is enabled for the entire domain.

File Structure
--------------

The basic file structure is shown below. Each line of output is represented by a line containing the output variable name(s). Loops indicate multiple lines of output.

.. parsed-literal::

   :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`
   :ref:`NDSETSW <NDSETSW>`, :ref:`NP <NP>`, :ref:`DTDP <DTDP>` * :ref:`NSPOOLGW <NSPOOLGW>`, :ref:`NSPOOLGW <NSPOOLGW>`, :ref:`IRTYPE <IRTYPE>`
   :ref:`TIME <TIME>`, :ref:`IT <IT>`
   for k=1, :ref:`NP <NP>`
      k, :ref:`PR2(k) <PR2>`
   end k loop

   Repeat the block above for each time step until the end of the simulation.

Notes
-----

* Output may be in ascii or binary format depending on how :ref:`NOUTGW <NOUTGW>` is set in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`
* If binary output is specified, the node number (k) is not included in the output 