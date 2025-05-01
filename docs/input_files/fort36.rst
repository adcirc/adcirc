Fort.36: Salinity Boundary Condition Input File
===============================================

The Salinity Boundary Condition Input File (``fort.36``) is used when the :ref:`RES_BC_FLAG <RES_BC_FLAG>` parameter is set to -2, 2, -4, or 4 in the ``fort.15`` file. This file specifies salinity values for ocean boundary nodes at different vertical levels.

File Structure
--------------

The file format is as follows:

.. parsed-literal::

    for i=1 to numberOfDataSets
        comment line (date)
        for k=1 to number_of_ocean_boundary_nodes
            k, (SALBC(k,m), m=1, :ref:`NFEN <NFEN>`)
        end k loop
    end i loop

Notes
-----

- This file is only required when running simulations with :ref:`RES_BC_FLAG <RES_BC_FLAG>` = -2, 2, -4, or 4
- The salinity values should be specified for all ocean boundary nodes
- Each dataset should include a comment line with the date
- SALBC(k,m) represents the salinity value at node k and vertical level m
- NFEN is the number of vertical levels in the model 