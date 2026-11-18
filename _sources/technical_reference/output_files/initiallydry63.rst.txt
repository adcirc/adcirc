Initiallydry.63: Inundation Data Output File
============================================

The initiallydry.63 file was created as the foundation of inundation output data by specifying the areas that ADCIRC considers wet and dry upon cold start. Those areas that are initially dry can then be considered as inundated areas when they become wet.

When the :ref:`inundationOutput <inundationOutput>` parameter is set to .true. in the optional inundationOutputContol namelist at the bottom of the :doc:`fort.15 <../input_files/fort15>` file, the initiallydry.63 file will be written at the beginning of the simulation whether the run is a cold start or hot start. The nodal values in the initiallydry.63 file are 1 if a node is dry at cold start, and 0 if the node is wet at coldstart. The data in the initiallydry.63 file represent areas that are dry at cold start, even if the run that produced the initiallydry.63 file was hotstarted.

The wet/dry state in the initiallydry.63 file takes into account the :ref:`initial_river_elevation <initial_river_elevation>` nodal attribute, the :ref:`surface_submergence_state <surface_submergence_state>` nodal attribute, and the :ref:`sea_surface_height_above_geoid <sea_surface_height_above_geoid>` nodal attribute from the nodal attributes (:doc:`fort.13 <../input_files/fort13>`) file, as well as the bathymetric depth from the mesh (:doc:`fort.14 <../input_files/fort14>`) file. It also includes the results of the landlocking algorithm, which dries any wet nodes that are completely surrounded by dry nodes. It does not include the effect of any tidal, meteorological, river, or other forcing.

File Structure
--------------

The basic file structure is shown below. Each line of output is represented by a line containing the output variable name(s). Loops indicate multiple lines of output.

.. parsed-literal::

   :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`
   1, :ref:`NP <NP>`, :ref:`DTDP <DTDP>`\*:ref:`NSPOOLGE <NSPOOLGE>`, :ref:`NSPOOLGE <NSPOOLGE>`, :ref:`IRTYPE <IRTYPE>`
   :ref:`TIME <TIME>`, :ref:`IT <IT>`
   for k=1, :ref:`NP <NP>`
      k, :ref:`initiallydry(k) <initiallydry>`
   end k loop

Notes
-----

* Output may be in ascii or netCDF format depending on how :ref:`NOUTGE <NOUTGE>` is set in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`
* The initiallydry.63 is written at the very beginning of a simulation run, before timestepping begins 