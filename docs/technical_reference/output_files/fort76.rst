Fort.76: Bathymetry Time Series at All Nodes
============================================

This file contains bathymetry time series data at all nodes in the model grid when ADCIRC's time varying bathymetry feature is activated. The output is controlled by the full domain time varying water surface elevation output parameters (:ref:`NOUTGE <NOUTGE>`, etc.) specified in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`.

File Structure
--------------

The file can be written in either ASCII or binary format, depending on how :ref:`NOUTGE <NOUTGE>` is set in the Model Parameter and Periodic Boundary Condition File. The basic structure is shown below:

.. parsed-literal::

    :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`
    :ref:`NDSETSE <NDSETSE>`, :ref:`NP <NP>`, :ref:`DTDP <DTDP>` * :ref:`NSPOOLGE <NSPOOLGE>`, :ref:`NSPOOLGE <NSPOOLGE>`, :ref:`IRTYPE <IRTYPE>`
    :ref:`TIME <TIME>`, :ref:`IT <IT>`
    for k=1, :ref:`NP <NP>`
        k, :ref:`dp(k) <dp>`
    end k loop

Notes
-----

* Output format (ASCII/binary) is determined by :ref:`NOUTGE <NOUTGE>` in the fort.15 file
* For binary output, the node number (k) is not included in the output
* Time series data is recorded for every node in the model grid
* The output frequency is controlled by :ref:`NSPOOLGE <NSPOOLGE>` parameter 