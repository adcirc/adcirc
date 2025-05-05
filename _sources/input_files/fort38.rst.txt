Fort.38: Surface Temperature Boundary Values Input File
=======================================================

The Surface Temperature Boundary Values Input File (``fort.38``) is used when the :ref:`RES_BC_FLAG <RES_BC_FLAG>` parameter is set to -3, 3, -4, or 4 in the ``fort.15`` file (i.e., when a lateral temperature boundary condition is being used). The format of this file depends on the :ref:`BCFLAG_TEMP <BCFLAG_TEMP>` parameter in the ``fort.15``, which controls the surface heat flux parameterization in ADCIRC.

File Structure
--------------

The file format varies depending on the value of :ref:`BCFLAG_TEMP <BCFLAG_TEMP>`:

1. If BCFLAG_TEMP=1:

.. parsed-literal::

    for i=1 to numberOfDataSets
        for k=1 to :ref:`NP <NP>`
            k, q_heat(k)
        end k loop
    end i loop

2. If BCFLAG_TEMP=2:

.. parsed-literal::

    for i=1 to numberOfDataSets
        (K, (TMP(K,J), J=1,6),K=1, :ref:`NP <NP>`)
    end i loop

3. If BCFLAG_TEMP=3:

.. parsed-literal::

    for i=1 to numberOfDataSets
        (K, (TMP(K,J), J=1,4),K=1, :ref:`NP <NP>`)
    end i loop

Notes
-----

- This file is only required when running simulations with :ref:`RES_BC_FLAG <RES_BC_FLAG>` = -3, 3, -4, or 4
- NP is the number of nodes in the horizontal mesh (i.e., the 2D fulldomain number of nodes)
- For BCFLAG_TEMP=2 and 3:
  - TMP(K,J) represents the Jth heat flux component for the Kth horizontal mesh node
  - The data are read using an implicit Fortran i/o loop, hence the parentheses
  - BCFLAG_TEMP=2 requires 6 heat flux components
  - BCFLAG_TEMP=3 requires 4 heat flux components
- For BCFLAG_TEMP=1:
  - q_heat(k) represents the heat flux at node k