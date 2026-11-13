Noff.100: Wet/Dry Elemental State File
======================================

This file records the wet/dry state of elements in the model grid, where a value of 1 indicates a wet element and 0 indicates a dry element at the time the dataset was written. The data comes from the elemental wet/dry array :ref:`NOFF <NOFF>` in ADCIRC. These data are primarily useful for ADCIRC developers working on experimental wet/dry algorithms.

The file is generated when the `outputNOFF` parameter is set to `.true.` in the optional `wetDryControl` namelist at the bottom of the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`. Output timing follows the same schedule as the :doc:`Water Surface Elevation <fort63>` file, using the parameters :ref:`TOUTSGE <TOUTSGE>`, :ref:`TOUTFGE <TOUTFGE>`, and :ref:`NSPOOLGE <NSPOOLGE>`.

File Structure
--------------

The file is only available in ASCII format. The basic structure is shown below:

.. parsed-literal::

    :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`

    :ref:`NDSETSE <NDSETSE>`, :ref:`NE <NE>`, :ref:`DTDP <DTDP>` * :ref:`NSPOOLGE <NSPOOLGE>`, :ref:`NSPOOLGE <NSPOOLGE>`, :ref:`IRTYPE <IRTYPE>`

    :ref:`TIME <TIME>`, :ref:`IT <IT>`

    for i=1, :ref:`NE <NE>`
        i, :ref:`NOFF(i) <NOFF>`
    end i loop

Notes
-----

* This is the only ADCIRC output file that produces elemental (rather than nodal) data
* Output is available only in ASCII format
* Values are integers: 1 for wet elements, 0 for dry elements
* Output timing matches fort.63 file using :ref:`NSPOOLGE <NSPOOLGE>` parameter
* File generation is controlled by `outputNOFF` in the `wetDryControl` namelist
* Data is primarily intended for wet/dry algorithm development and testing 