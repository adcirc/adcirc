Nodecode.63: Wet/Dry State of Nodes
===================================

The nodecode.63 file records the wet/dry state of nodes where 1 indicates a node is categorized as wet on the timestep that the dataset was written while a value of 0 indicates that a node is categorized as dry. These data are generally only valuable to ADCIRC developers who are working on experimental wet/dry algorithms.

The writing of the nodecode.63 output file is activated when the :ref:`outputNodeCode <outputNodeCode>` parameter is set to .true. in the optional wetDryContol namelist at the bottom of the :doc:`fort.15 <../input_files/fort15>` file.

Output is only available in the ascii format. The data are produced on the same schedule as the full domain water surface elevation (:doc:`fort.63 <fort63>`) file, i.e., the values of :ref:`TOUTSGE <TOUTSGE>`, :ref:`TOUTFGE <TOUTFGE>`, and :ref:`NSPOOLGE <NSPOOLGE>` are used.

File Structure
--------------

The basic file structure is shown below. Each line of output is represented by a line containing the output variable name(s). Loops indicate multiple lines of output.

.. parsed-literal::

   :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`
   :ref:`NDSETSE <NDSETSE>`, :ref:`NP <NP>`, :ref:`DTDP <DTDP>`\*:ref:`NSPOOLGE <NSPOOLGE>`, :ref:`NSPOOLGE <NSPOOLGE>`, :ref:`IRTYPE <IRTYPE>`
   :ref:`TIME <TIME>`, :ref:`IT <IT>`
   for k=1, :ref:`NP <NP>`
      k, :ref:`nodecode(k) <nodecode>`
   end k loop

Note
----

* The nodecode.63 file contains integer values 