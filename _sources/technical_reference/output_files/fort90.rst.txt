.. _fort90:

Fort.90: Primitive Weighting in Continuity Equation
===================================================

This file contains primitive weighting in continuity equation time series data at all nodes in the model grid as specified in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`. The output format and frequency settings are controlled by the same parameters used for the :doc:`Elevation Time Series <fort63>` file.

File Structure
--------------

The file can be written in either ASCII or NetCDF format, depending on how :ref:`NOUTGE <NOUTGE>` is set in the Model Parameter and Periodic Boundary Condition File. The basic structure is shown below:

.. parsed-literal::

    :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`

    :ref:`NDSETSE <NDSETSE>`, :ref:`NP <NP>`, :ref:`DTDP <DTDP>` * :ref:`NSPOOLGE <NSPOOLGE>`, :ref:`NSPOOLGE <NSPOOLGE>`, :ref:`IRTYPE <IRTYPE>`

    :ref:`TIME <TIME>`, :ref:`IT <IT>`

    for k=1, :ref:`NP <NP>`
        k, :ref:`tau0var(k) <tau0var>`
    end k loop

Notes
-----

* Output format (ASCII/NetCDF) is determined by :ref:`NOUTGE <NOUTGE>` in the fort.15 file
* For binary/NetCDF output, the node number (k) is not included in the output
* Time series data is recorded for every node in the model grid
* The output frequency is controlled by :ref:`NSPOOLGE <NSPOOLGE>` parameter, matching fort.63 settings
* The variable tau0var represents the primitive weighting value in the continuity equation at each node 