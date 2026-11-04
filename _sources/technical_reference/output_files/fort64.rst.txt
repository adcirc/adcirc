.. _fort64:

Fort.64: Depth-averaged Velocity Time Series at All Nodes in the Model Grid
===========================================================================

Depth-averaged velocity time series output at all nodes in the model grid as specified in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`. This file is generated when time series output is enabled for the entire model domain.

File Structure
--------------

The basic file structure is shown below. Each line of output is represented by a line containing the output variable name(s). Loops indicate multiple lines of output.

.. parsed-literal::

   :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`
   :ref:`NDSETSV <NDSETSV>`, :ref:`NP <NP>`, :ref:`DTDP <DTDP>`\*:ref:`NSPOOLGV <NSPOOLGV>`, :ref:`NSPOOLGV <NSPOOLGV>`, :ref:`IRTYPE <IRTYPE>`
   :ref:`TIME <TIME>`, :ref:`IT <IT>`
   for k=1, :ref:`NP <NP>`
      k, :ref:`UU2(k) <UU>`, :ref:`VV2(k) <VV>`
   end k loop

   Repeat the block above for each time step until the end of the simulation.

Notes
-----

* Output may be in ascii or binary format depending on how :ref:`NOUTGV <NOUTGV>` is set in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`
* If binary output is specified, the node number (k) is not included in the output 