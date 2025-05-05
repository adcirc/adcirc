.. _fort44:

Fort.44: 3D Density, Temperature and/or Salinity at All Nodes in the Model Grid
===============================================================================

This file contains density, temperature and salinity time series output at all nodes in the model grid as specified in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`.

The file is generated in ASCII format and its structure depends on the value of :ref:`IDEN <iden>`.

File Structure
--------------

The basic file structure is shown below. Each line of output is represented by a line containing the output variable name(s). Loops indicate multiple lines of output.

If :ref:`IDEN <iden>` = 1 or -1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. parsed-literal::

   :ref:`RUNDES <rundes>`, :ref:`RUNID <runid>`, :ref:`AGRID <agrid>`
   :ref:`NDSET3DGD <ndset3dgd>`, :ref:`NP <np>`, :ref:`DTDP <dtdp>` * :ref:`NSPO3DGD <nspo3dgd>`, :ref:`NSPO3DGD <nspo3dgd>`, :ref:`NFEN <nfen>`, :ref:`IDEN <iden>`
   for k = 1 to :ref:`NDSET3DGD <ndset3dgd>`
      :ref:`TIME <time>`, :ref:`IT <it>`, (:ref:`SIGMA <sigma>` (N), N=1,:ref:`NFEN <nfen>`-1)
      for j = 1 to :ref:`NP <np>`
         J, (:ref:`SIGT <sigt>` (M), M=1, :ref:`NFEN <nfen>`)
      end j loop
   end k loop

If :ref:`IDEN <iden>` = 2 or -2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. parsed-literal::

   :ref:`RUNDES <rundes>`, :ref:`RUNID <runid>`, :ref:`AGRID <agrid>`
   :ref:`NDSET3DGD <ndset3dgd>`, :ref:`NP <np>`, :ref:`DTDP <dtdp>` * :ref:`NSPO3DGD <nspo3dgd>`, :ref:`NSPO3DGD <nspo3dgd>`, :ref:`NFEN <nfen>`, :ref:`IRTYPE <irtype>`
   for k = 1 to :ref:`NDSET3DGD <ndset3dgd>`
      :ref:`TIME <time>`, :ref:`IT <it>`, (:ref:`SIGMA <sigma>` (N), N=1,:ref:`NFEN <nfen>`-1), :ref:`SIGMA <sigma>` (:ref:`NFEN <nfen>`)
      for j = 1 to :ref:`NP <np>`
         J, (:ref:`SIGT <sigt>` (M), :ref:`SAL <sal>` (M), M=1, :ref:`NFEN <nfen>`)
      end j loop
   end k loop

If :ref:`IDEN <iden>` = 3 or -3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. parsed-literal::

   :ref:`RUNDES <rundes>`, :ref:`RUNID <runid>`, :ref:`AGRID <agrid>`
   :ref:`NDSET3DGD <ndset3dgd>`, :ref:`NP <np>`, :ref:`DTDP <dtdp>` * :ref:`NSPO3DGD <nspo3dgd>`, :ref:`NSPO3DGD <nspo3dgd>`, :ref:`NFEN <nfen>`, :ref:`IRTYPE <irtype>`
   for k = 1 to :ref:`NDSET3DGD <ndset3dgd>`
      :ref:`TIME <time>`, :ref:`IT <it>`, (:ref:`SIGMA <sigma>` (N), N=1,:ref:`NFEN <nfen>`-1), :ref:`SIGMA <sigma>` (:ref:`NFEN <nfen>`)
      for j = 1 to :ref:`NP <np>`
         J, (:ref:`SIGT <sigt>` (M), :ref:`TEMP <temp>` (M), M=1, :ref:`NFEN <nfen>`)
      end j loop
   end k loop

If :ref:`IDEN <iden>` = 4 or -4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. parsed-literal::

   :ref:`RUNDES <rundes>`, :ref:`RUNID <runid>`, :ref:`AGRID <agrid>`
   :ref:`NDSET3DGD <ndset3dgd>`, :ref:`NP <np>`, :ref:`DTDP <dtdp>` * :ref:`NSPO3DGD <nspo3dgd>`, :ref:`NSPO3DGD <nspo3dgd>`, :ref:`NFEN <nfen>`, :ref:`IRTYPE <irtype>`
   for k = 1 to :ref:`NDSET3DGD <ndset3dgd>`
      :ref:`TIME <time>`, :ref:`IT <it>`, (:ref:`SIGMA <sigma>` (N), N=1,:ref:`NFEN <nfen>`-1), :ref:`SIGMA <sigma>` (:ref:`NFEN <nfen>`), :ref:`SIGMA <sigma>` (:ref:`NFEN <nfen>`)
      for j = 1 to :ref:`NP <np>`
         J, (:ref:`SIGT <sigt>` (M), :ref:`TEMP <temp>` (M), :ref:`SAL <sal>` (M), M=1, :ref:`NFEN <nfen>`)
      end j loop
   end k loop 