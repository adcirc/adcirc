Inundationtime.63: Inundation Time File
=======================================

The inundationtime.63 file records the total time that an initially dry area is inundated beyond a certain inundation threshold depth specified in the :doc:`fort.15 <../input_files/fort15>` file in the inundationOutputControl namelist by the parameter :ref:`inunThresh <inunThresh>`.

The writing of the inundationtime.63 output file is activated when the :ref:`inundationOutput <inundationOutput>` parameter is set to .true. in the optional inundationOutputContol namelist at the bottom of the :doc:`fort.15 <../input_files/fort15>` file.

The file contains two data sets:
1. Total accumulated time in seconds that a node was inundated beyond the threshold (periods of inundation are counted toward the total time, even if they are not contiguous)
2. Time of onset of inundation beyond the threshold in seconds since cold start (useful in the context of real time model guidance for decision making)

File Structure
--------------

The basic file structure is shown below. Each line of output is represented by a line containing the output variable name(s). Loops indicate multiple lines of output.

.. parsed-literal::

   :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`
   2, :ref:`NP <NP>`, :ref:`DTDP <DTDP>`\*:ref:`NSPOOLGE <NSPOOLGE>`, :ref:`NSPOOLGE <NSPOOLGE>`, :ref:`IRTYPE <IRTYPE>`
   :ref:`TIME <TIME>`, :ref:`IT <IT>`
   for k=1, :ref:`NP <NP>`
      k, :ref:`inundationtime(k) <inundationtime>`
   end k loop

   :ref:`TIME <TIME>`, :ref:`IT <IT>`
   for k=1, :ref:`NP <NP>`
      k, :ref:`inundationtime_onset(k) <inundationtime_onset>`
   end k loop

Notes
-----

* Output may be in ascii or netCDF format depending on how :ref:`NOUTGE <NOUTGE>` is set in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`
* The inundationtime.63 file is written at the very end of the simulation, after timestepping is complete
* The values only reflect the current run, even if the run was hotstarted
* The data only record maximum inundation depth in areas that are initially dry, according to the :doc:`initiallydry.63 <initiallydry63>` file 