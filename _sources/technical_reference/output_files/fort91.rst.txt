Fort.91: Ice Coverage Fields at Recording Stations
==================================================

This file contains ice coverage field data at specified meteorological recording stations when ADCIRC's ice coverage input feature is activated. The output is controlled by the meteorological recording station parameters (:ref:`NOUTM <NOUTM>`, etc.) specified in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`.

File Structure
--------------

The file can be written in either ASCII or binary format, depending on how :ref:`NOUTM <NOUTM>` is set in the Model Parameter and Periodic Boundary Condition File. The basic structure is shown below:

.. parsed-literal::

    :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`
    :ref:`NTRSPM <NTRSPM>`, :ref:`NSTAM <NSTAM>`, :ref:`DTDP <DTDP>` * :ref:`NSPOOLM <NSPOOLM>`, :ref:`NSPOOLM <NSPOOLM>`, :ref:`IRTYPE <IRTYPE>`
    :ref:`TIME <TIME>`, :ref:`IT <IT>`
    for k=1, :ref:`NSTAM <NSTAM>`
        k, RMICE00(k)
    end k loop

Notes
-----

* Output format (ASCII/binary) is determined by :ref:`NOUTM <NOUTM>` in the fort.15 file
* For binary output, the station number (k) is not included in the output
* Time series data is recorded at stations specified for meteorological recording
* The output frequency is controlled by :ref:`NSPOOLM <NSPOOLM>` parameter
* The variable RMICE00 represents the ice coverage field value at each recording station 