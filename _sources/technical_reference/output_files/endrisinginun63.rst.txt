.. _endrisinginun63:

Endrisinginun.63: Inundation Rising at the End of the Run Flag File
===================================================================

The endrisinginun.63 file flags nodes whose inundation depth is rising at the end of the simulation by comparing the water surface elevation on the final time step with the water surface elevation on the previous time step. Nodes with rising inundation levels are flagged with an integer value of 1 and all others are given an integer value of 0.

The writing of the endrisinginun.63 output file is activated when the :ref:`inundationOutput <inundationOutput>` parameter is set to .true. in the optional inundationOutputContol namelist at the bottom of the :doc:`fort.15 <../input_files/fort15>` file.

File Structure
--------------

The basic file structure is shown below. Each line of output is represented by a line containing the output variable name(s). Loops indicate multiple lines of output.

.. parsed-literal::

   :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`
   1, :ref:`NP <NP>`, :ref:`DTDP <DTDP>`\*:ref:`NSPOOLGE <NSPOOLGE>`, :ref:`NSPOOLGE <NSPOOLGE>`, :ref:`IRTYPE <IRTYPE>`
   :ref:`TIME <TIME>`, :ref:`IT <IT>`
   for k=1, :ref:`NP <NP>`
      k, :ref:`endrisinginun(k) <endrisinginun>`
   end k loop

Notes
-----

* Output may be in ascii or netCDF format depending on how :ref:`NOUTGE <NOUTGE>` is set in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`
* The endrisinginun.63 file is written at the very end of the simulation, after timestepping is complete
* The values only reflect the current run, even if the run was hotstarted
* The data only flag rising water surface elevation in areas that are initially dry, according to the :doc:`initiallydry.63 <initiallydry63>` file 