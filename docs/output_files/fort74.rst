Fort.74: Wind Stress or VelocityTime Series at All Nodes in the Model Grid
==========================================================================

Wind velocity or stress time series output at all nodes in the model grid as specified in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`.

File Structure
--------------

The basic file structure is shown below. Each line of output is represented by a line containing the output variable name(s). Loops indicate multiple lines of output.

.. parsed-literal::

   :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`
   :ref:`NDSETSW <NDSETSW>`, :ref:`NP <NP>`, :ref:`DTDP <DTDP>`\*:ref:`NSPOOLGW <NSPOOLGW>`, :ref:`NSPOOLGW <NSPOOLGW>`, :ref:`IRTYPE <IRTYPE>`
   :ref:`TIME <TIME>`, :ref:`IT <IT>`
   for k=1, :ref:`NP <NP>`
      k, :ref:`WVNXOUT(k) <WVNXOUT>`, :ref:`WVNYOUT(k) <WVNYOUT>`
   end k loop

   Repeat the block above for each time step until the end of the simulation.

Note
----

* Output may be in ascii or binary format depending on how :ref:`NOUTGW <NOUTGW>` is set in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`
* If binary output is specified, the node number (k) is not included in the output 