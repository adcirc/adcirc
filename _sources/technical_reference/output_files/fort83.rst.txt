Fort.83: Scalar Concentration Time Series at All Nodes
======================================================

This file contains scalar concentration time series data at all nodes in the model grid as specified in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`. 

.. note::
   This feature is currently not supported in ADCIRC.

File Structure
--------------

The file can be written in either ASCII or binary format, depending on how :ref:`NOUTGC <NOUTGC>` is set in the Model Parameter and Periodic Boundary Condition File. The basic structure is shown below:

.. parsed-literal::

    :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`

    :ref:`NDSETSC <NDSETSC>`, :ref:`NP <NP>`, :ref:`DTDP <DTDP>` * :ref:`NSPOOLGC <NSPOOLGC>`, :ref:`NSPOOLGC <NSPOOLGC>`, :ref:`IRTYPE <IRTYPE>`

    :ref:`TIME <TIME>`, :ref:`IT <IT>`

    for k=1, :ref:`NP <NP>`
        k, :ref:`C1(k) <C1>`
    end k loop

Notes
-----

* Output format (ASCII/binary) is determined by :ref:`NOUTGC <NOUTGC>` in the fort.15 file
* For binary output, the node number (k) is not included in the output
* Time series data is recorded for every node in the model grid
* The output frequency is controlled by :ref:`NSPOOLGC <NSPOOLGC>` parameter
* This feature is currently under development and not yet supported in ADCIRC 