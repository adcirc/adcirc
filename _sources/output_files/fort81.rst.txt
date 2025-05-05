Fort.81: Scalar Concentration Time Series at Recording Stations
===============================================================

This file contains scalar concentration time series data at specified recording stations as defined in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`. 

.. note::
   This feature is currently not supported in ADCIRC.

File Structure
--------------

The file can be written in either ASCII or binary format, depending on how :ref:`NOUTC <NOUTC>` is set in the Model Parameter and Periodic Boundary Condition File. The basic structure is shown below:

.. parsed-literal::

    :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`

    :ref:`NTRSPC <NTRSPC>`, :ref:`NSTAC <NSTAC>`, :ref:`DTDP <DTDP>` * :ref:`NSPOOLC <NSPOOLC>`, :ref:`NSPOOLC <NSPOOLC>`, :ref:`IRTYPE <IRTYPE>`

    :ref:`TIME <TIME>`, :ref:`IT <IT>`

    for k=1, :ref:`NSTAC <NSTAC>`
        k, :ref:`CC00(k) <CC00>`
    end k loop

Notes
-----

* Output format (ASCII/binary) is determined by :ref:`NOUTC <NOUTC>` in the fort.15 file
* For binary output, the station number (k) is not included in the output
* Time series data is recorded at stations specified for concentration recording
* The output frequency is controlled by :ref:`NSPOOLC <NSPOOLC>` parameter
* This feature is currently under development and not yet supported in ADCIRC 