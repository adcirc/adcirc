Everdried.63: Dry Node Flagging File
====================================

The everdried.63 file was created to support harmonic analysis by flagging all nodes that had ever become dry during the course of a simulation. These data are useful to harmonic analysis because a node that goes dry for a single time step has a -99999 recorded for its water surface elevation, which contaminates the harmonic analysis solution.

The writing of the everdried.63 output file is activated when the :ref:`inundationOutput <inundationOutput>` parameter is set to .true. in the optional inundationOutputContol namelist at the bottom of the :doc:`fort.15 <../input_files/fort15>` file.

The file contains two data sets:
1. Wet/dry state information: -99999.0 if a node ever went dry during the simulation, 1.0 if it was wet for the entire simulation
2. Total time in seconds that a node was dry during the simulation (0.0 if it was always wet)

File Structure
--------------

The basic file structure is shown below. Each line of output is represented by a line containing the output variable name(s). Loops indicate multiple lines of output.

.. parsed-literal::

   :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`
   2, :ref:`NP <NP>`, :ref:`DTDP <DTDP>`\*:ref:`NSPOOLGE <NSPOOLGE>`, :ref:`NSPOOLGE <NSPOOLGE>`, :ref:`IRTYPE <IRTYPE>`
   :ref:`TIME <TIME>`, :ref:`IT <IT>`
   for k=1, :ref:`NP <NP>`
      k, :ref:`everdried(k) <everdried>`
   end k loop

   :ref:`TIME <TIME>`, :ref:`IT <IT>`
   for k=1, :ref:`NP <NP>`
      k, :ref:`driedtime(k) <driedtime>`
   end k loop

Notes
-----

* Output may be in ascii or netCDF format depending on how :ref:`NOUTGE <NOUTGE>` is set in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`
* The everdried.63 file is written at the very end of the simulation, after timestepping is complete
* The values only reflect the current run, even if the run was hotstarted 