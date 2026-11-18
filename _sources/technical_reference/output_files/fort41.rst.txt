Fort.41: 3D Density, Temperature and/or Salinity at Recording Stations
======================================================================

This file contains density, temperature and salinity time series output at recording stations as specified in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`. The file structure varies based on the value of :ref:`IDEN <IDEN>`.

File Structure
--------------

The file is written in ASCII format. The basic structure varies depending on the value of IDEN:

For IDEN = ±1 (Density Only)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. parsed-literal::

    :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`

    :ref:`NDSET3DSD <NDSET3DSD>`, :ref:`NSTA3DD <NSTA3DD>`, :ref:`DTDP <DTDP>` * :ref:`NSPO3DSD <NSPO3DSD>`, :ref:`NSPO3DSD <NSPO3DSD>`, :ref:`NFEN <NFEN>`, :ref:`IRTYPE <IRTYPE>`

    for k=1, :ref:`NDSET3DSD <NDSET3DSD>`
        :ref:`TIME <TIME>`, :ref:`IT <IT>`, (:ref:`SIGMA(N) <SIGMA>`, N=1, :ref:`NFEN <NFEN>`-1)

        for j=1, :ref:`NSTA3DD <NSTA3DD>`
            j, (:ref:`SIGTSTA(M) <SIGTSTA>`, M=1, :ref:`NFEN <NFEN>`)
        end j loop
    end k loop

For IDEN = ±2 (Density and Salinity)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. parsed-literal::

    :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`

    :ref:`NDSET3DSD <NDSET3DSD>`, :ref:`NSTA3DD <NSTA3DD>`, :ref:`DTDP <DTDP>` * :ref:`NSPO3DSD <NSPO3DSD>`, :ref:`NSPO3DSD <NSPO3DSD>`, :ref:`NFEN <NFEN>`, :ref:`IRTYPE <IRTYPE>`

    for k=1, :ref:`NDSET3DSD <NDSET3DSD>`
        :ref:`TIME <TIME>`, :ref:`IT <IT>`, (:ref:`SIGMA(N) <SIGMA>`, :ref:`SIGMA(N) <SIGMA>`, N=1, :ref:`NFEN <NFEN>`-1), :ref:`SIGMA(NFEN) <SIGMA>`

        for j=1, :ref:`NSTA3DD <NSTA3DD>`
            j, (:ref:`SIGTSTA(M) <SIGTSTA>`, :ref:`SALSTA(M) <SALSTA>`, M=1, :ref:`NFEN <NFEN>`)
        end j loop
    end k loop

For IDEN = ±3 (Density and Temperature)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. parsed-literal::

    :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`

    :ref:`NDSET3DSD <NDSET3DSD>`, :ref:`NSTA3DD <NSTA3DD>`, :ref:`DTDP <DTDP>` * :ref:`NSPO3DSD <NSPO3DSD>`, :ref:`NSPO3DSD <NSPO3DSD>`, :ref:`NFEN <NFEN>`, :ref:`IRTYPE <IRTYPE>`

    for k=1, :ref:`NDSET3DSD <NDSET3DSD>`
        :ref:`TIME <TIME>`, :ref:`IT <IT>`, (:ref:`SIGMA(N) <SIGMA>`, :ref:`SIGMA(N) <SIGMA>`, N=1, :ref:`NFEN <NFEN>`-1), :ref:`SIGMA(NFEN) <SIGMA>`

        for j=1, :ref:`NSTA3DD <NSTA3DD>`
            j, (:ref:`SIGTSTA(M) <SIGTSTA>`, :ref:`TEMPSTA(M) <TEMPSTA>`, M=1, :ref:`NFEN <NFEN>`)
        end j loop
    end k loop

For IDEN = ±4 (Density, Temperature, and Salinity)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. parsed-literal::

    :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`

    :ref:`NDSET3DSD <NDSET3DSD>`, :ref:`NSTA3DD <NSTA3DD>`, :ref:`DTDP <DTDP>` * :ref:`NSPO3DSD <NSPO3DSD>`, :ref:`NSPO3DSD <NSPO3DSD>`, :ref:`NFEN <NFEN>`, :ref:`IRTYPE <IRTYPE>`

    for k=1, :ref:`NDSET3DSD <NDSET3DSD>`
        :ref:`TIME <TIME>`, :ref:`IT <IT>`, (:ref:`SIGMA(N) <SIGMA>`, :ref:`SIGMA(N) <SIGMA>`, :ref:`SIGMA(N) <SIGMA>`, N=1, :ref:`NFEN <NFEN>`-1), :ref:`SIGMA(NFEN) <SIGMA>`, :ref:`SIGMA(NFEN) <SIGMA>`

        for j=1, :ref:`NSTA3DD <NSTA3DD>`
            j, (:ref:`SIGTSTA(M) <SIGTSTA>`, :ref:`TEMPSTA(M) <TEMPSTA>`, :ref:`SALSTA(M) <SALSTA>`, M=1, :ref:`NFEN <NFEN>`)
        end j loop
    end k loop

Notes
-----

* Output is only available in ASCII format
* Time series data is recorded at specified recording stations
* The file structure varies based on the :ref:`IDEN <IDEN>` parameter:
    * IDEN = ±1: Density only
    * IDEN = ±2: Density and salinity
    * IDEN = ±3: Density and temperature
    * IDEN = ±4: Density, temperature, and salinity
* Data is recorded at multiple vertical levels defined by SIGMA values
* SIGTSTA represents density values
* TEMPSTA represents temperature values
* SALSTA represents salinity values 