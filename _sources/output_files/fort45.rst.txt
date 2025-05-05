Fort.45: 3D Velocity at All Nodes
=================================

This file contains velocity time series output at all nodes in the model grid as specified in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`.

File Structure
--------------

The file is written in ASCII format. The basic structure is shown below:

.. parsed-literal::

    :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`

    :ref:`NDSET3DGV <NDSET3DGV>`, :ref:`NP <NP>`, :ref:`DTDP <DTDP>` * :ref:`NSPO3DGV <NSPO3DGV>`, :ref:`NSPO3DGV <NSPO3DGV>`, :ref:`NFEN <NFEN>`, :ref:`IRTYPE <IRTYPE>`

    for k=1, :ref:`NDSET3DGV <NDSET3DGV>`
        :ref:`TIME <TIME>`, :ref:`IT <IT>`, (:ref:`SIGMA(N) <SIGMA>`, :ref:`SIGMA(N) <SIGMA>`, :ref:`SIGMA(N) <SIGMA>`, N=1, :ref:`NFEN <NFEN>`-1), :ref:`SIGMA(NFEN) <SIGMA>`, :ref:`SIGMA(NFEN) <SIGMA>`

        for j=1, :ref:`NP <NP>`
            j, (:ref:`REAL(Q(j,M)) <REAL(Q(j,M))>`, :ref:`AIMAG(Q(j,M)) <AIMAG(Q(j,M))>`, :ref:`WZ(j,M) <WZ(j,M)>`, M=1, :ref:`NFEN <NFEN>`)
        end j loop
    end k loop

Notes
-----

* Output is only available in ASCII format
* Time series data is recorded at all nodes in the model grid
* Velocity components are stored as:
    * REAL(Q): x-component of velocity
    * AIMAG(Q): y-component of velocity
    * WZ: vertical velocity component
* Data is recorded at multiple vertical levels defined by SIGMA values
* Similar to fort.42 but provides data for the entire model domain instead of specific stations 