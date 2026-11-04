.. _fort46:

Fort.46: 3D Turbulence at All Nodes
===================================

This file contains turbulence time series output at all nodes in the model grid as specified in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`.

File Structure
--------------

The file is written in ASCII format. The basic structure is shown below:

.. parsed-literal::

    :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`

    :ref:`NDSET3DGT <NDSET3DGT>`, :ref:`NP <NP>`, :ref:`DTDP <DTDP>` * :ref:`NSPO3DGT <NSPO3DGT>`, :ref:`NSPO3DGT <NSPO3DGT>`, :ref:`NFEN <NFEN>`, :ref:`IRTYPE <IRTYPE>`

    for k=1, :ref:`NDSET3DGT <NDSET3DGT>`
        :ref:`TIME <TIME>`, :ref:`IT <IT>`, (:ref:`SIGMA(N) <SIGMA>`, :ref:`SIGMA(N) <SIGMA>`, :ref:`SIGMA(N) <SIGMA>`, N=1, :ref:`NFEN <NFEN>`-1), :ref:`SIGMA(NFEN) <SIGMA>`, :ref:`SIGMA(NFEN) <SIGMA>`

        for j=1, :ref:`NP <NP>`
            j, (:ref:`q20(j,M) <q20(j,M)>`, :ref:`l(j,M) <l(j,M)>`, :ref:`EV(j,M) <EV(j,M)>`, M=1, :ref:`NFEN <NFEN>`)
        end j loop
    end k loop

Notes
-----

* Output is only available in ASCII format
* Time series data is recorded at all nodes in the model grid
* Turbulence parameters are stored as:
    * q20: turbulent kinetic energy
    * l: turbulent length scale
    * EV: vertical eddy viscosity
* Data is recorded at multiple vertical levels defined by SIGMA values
* Similar to fort.43 but provides data for the entire model domain instead of specific stations 