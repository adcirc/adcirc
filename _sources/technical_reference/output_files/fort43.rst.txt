Fort.43: 3D Turbulence at Specified Recording Stations
======================================================

This file contains turbulence time series output at the recording stations as specified in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`.

File Structure
--------------

The file is written in ASCII format. The basic structure is shown below:

.. parsed-literal::

    :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`

    :ref:`NDSET3DST <NDSET3DST>`, :ref:`NSTA3DT <NSTA3DT>`, :ref:`DTDP <DTDP>` * :ref:`NSPO3DST <NSPO3DST>`, :ref:`NSPO3DST <NSPO3DST>`, :ref:`NFEN <NFEN>`, :ref:`IRTYPE <IRTYPE>`

    for k=1, :ref:`NDSET3DST <NDSET3DST>`
        :ref:`TIME <TIME>`, :ref:`IT <IT>`, (:ref:`SIGMA(N) <SIGMA>`, :ref:`SIGMA(N) <SIGMA>`, :ref:`SIGMA(N) <SIGMA>`, N=1, :ref:`NFEN <NFEN>`-1), :ref:`SIGMA(NFEN) <SIGMA>`, :ref:`SIGMA(NFEN) <SIGMA>`

        for j=1, :ref:`NSTA3DT <NSTA3DT>`
            j, (:ref:`q20STA(M) <q20STA>`, :ref:`ISTA(M) <ISTA>`, :ref:`EVSTA(M) <EVSTA>`, M=1, :ref:`NFEN <NFEN>`)
        end j loop
    end k loop

Notes
-----

* Output is only available in ASCII format
* Time series data is recorded at specified recording stations
* The output includes turbulence parameters at multiple vertical levels defined by SIGMA values 