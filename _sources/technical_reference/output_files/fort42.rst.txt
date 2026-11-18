Fort.42: 3D Velocity at Recording Stations
==========================================

This file contains velocity time series output at recording stations as specified in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`.

File Structure
--------------

The file is written in ASCII format. The basic structure is shown below:

.. parsed-literal::

    :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`

    :ref:`NDSET3DSV <NDSET3DSV>`, :ref:`NSTA3DV <NSTA3DV>`, :ref:`DTDP <DTDP>` * :ref:`NSPO3DSV <NSPO3DSV>`, :ref:`NSPO3DSV <NSPO3DSV>`, :ref:`NFEN <NFEN>`, :ref:`IRTYPE <IRTYPE>`

    for k=1, :ref:`NDSET3DSV <NDSET3DSV>`
        :ref:`TIME <TIME>`, :ref:`IT <IT>`, (:ref:`SIGMA(N) <SIGMA>`, :ref:`SIGMA(N) <SIGMA>`, :ref:`SIGMA(N) <SIGMA>`, N=1, :ref:`NFEN <NFEN>`-1), :ref:`SIGMA(NFEN) <SIGMA>`, :ref:`SIGMA(NFEN) <SIGMA>`

        for j=1, :ref:`NSTA3DV <NSTA3DV>`
            j, (:ref:`REAL(QSTA(M)) <REAL(QSTA(k,j))>`, :ref:`AIMAG(QSTA(M)) <AIMAG(QSTA(k,j))>`, :ref:`WZSTA(k,j) <WZSTA>`, M=1, :ref:`NFEN <NFEN>`)
        end j loop
    end k loop

Notes
-----

* Output is only available in ASCII format
* Time series data is recorded at specified recording stations
* Velocity components are stored as:
    * REAL(QSTA): x-component of velocity
    * AIMAG(QSTA): y-component of velocity
    * WZSTA: vertical velocity component
* Data is recorded at multiple vertical levels defined by SIGMA values 