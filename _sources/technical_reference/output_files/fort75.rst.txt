Fort.75: Bathymetry Time Series at Recording Stations
=====================================================

This file contains bathymetry time series data at specified recording stations when ADCIRC's time varying bathymetry feature is activated. The output is controlled by the water surface elevation recording station parameters (:ref:`NOUTE <NOUTE>`, etc.) specified in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`.

File Structure
--------------

The file can be written in either ASCII or binary format, depending on how :ref:`NOUTE <NOUTE>` is set in the Model Parameter and Periodic Boundary Condition File. The basic structure is shown below:

.. parsed-literal::

    :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`
    :ref:`NTRSPE <NTRSPE>`, :ref:`NSTAE <NSTAE>`, :ref:`DTDP <DTDP>`\*:ref:`NSPOOLE <NSPOOLE>`, :ref:`NSPOOLE <NSPOOLE>`, :ref:`IRTYPE <IRTYPE>`
    :ref:`TIME <TIME>`, :ref:`IT <IT>`
    for k=1, :ref:`NSTAE <NSTAE>`
        k, :ref:`DP00(k) <DP>`
    end k loop

Notes
-----

* Output format (ASCII/binary) is determined by :ref:`NOUTE <NOUTE>` in the fort.15 file
* For binary output, the station number (k) is not included in the output
* Time series data is recorded at stations specified for water surface elevation recording
* The output frequency is controlled by :ref:`NSPOOLE <NSPOOLE>` parameter 