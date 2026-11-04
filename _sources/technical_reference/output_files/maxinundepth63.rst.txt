Maxinundepth.63: Maximum Inundation Depth File
==============================================

The maximum inundation depth (maxinundepth.63) file records the peak inundation depth (in meters) above ground that occurred during the simulation. The data are only recorded in areas that are initially dry according to the :doc:`initiallydry.63 <initiallydry63>` file. The values are -99999.0 if the area was never wet during the simulation.

The writing of the maxinundepth.63 output file is activated when the :ref:`inundationOutput <inundationOutput>` parameter is set to .true. in the optional inundationOutputContol namelist at the bottom of the :doc:`fort.15 <../input_files/fort15>` file.

The file contains two data sets:
1. Peak inundation value
2. Time of occurrence of the peak inundation value in seconds since cold start

File Structure
--------------

The basic file structure is shown below. Each line of output is represented by a line containing the output variable name(s). Loops indicate multiple lines of output.

.. parsed-literal::

   :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`
   2, :ref:`NP <NP>`, :ref:`DTDP <DTDP>`\*:ref:`NSPOOLGE <NSPOOLGE>`, :ref:`NSPOOLGE <NSPOOLGE>`, :ref:`IRTYPE <IRTYPE>`
   :ref:`TIME <TIME>`, :ref:`IT <IT>`
   for k=1, :ref:`NP <NP>`
      k, :ref:`maxinundepth(k) <maxinundepth>`
   end k loop

   :ref:`TIME <TIME>`, :ref:`IT <IT>`
   for k=1, :ref:`NP <NP>`
      k, :ref:`maxinundepth_time(k) <maxinundepth_time>`
   end k loop

Notes
-----

* Output may be in ascii or netCDF format depending on how :ref:`NOUTGE <NOUTGE>` is set in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`
* The maxinundepth.63 file is written at the very end of the simulation, after timestepping is complete
* The values only reflect the current run, even if the run was hotstarted
* The data only record maximum inundation depth in areas that are initially dry, according to the :doc:`initiallydry.63 <initiallydry63>` file 