Fort.39: Salinity and Temperature River Boundary Values Input File
==================================================================

The Salinity and Temperature River Boundary Values Input File (``fort.39``) is used when the mesh file (``fort.14``) contains a baroclinic river boundary (IBTYPE=122) and :ref:`IDEN <IDEN>` is positive.

File Structure
--------------

The file format is as follows:

.. parsed-literal::

    RIVBCTIMINC, RIVBCSTATIM

    for i=1 to numberOfDataSets
        for j=1 to number_of_river_boundary_nodes
            (bc(j,k), k=1, :ref:`NFEN <NFEN>`)
        end j loop
    end i loop

where:
- RIVBCTIMINC is the time increment (in seconds) between the boundary condition datasets
- RIVBCSTATIM is the time (in seconds) when the boundary condition data start, relative to the cold start time

Notes
-----

- This file is only required when using baroclinic river boundaries (IBTYPE=122) with positive :ref:`IDEN <IDEN>` values
- The bc(j,k) values depend on the value of :ref:`IDEN <IDEN>`:
  - If IDEN=2: bc(j,k) represents salinity boundary condition values
  - If IDEN=3: bc(j,k) represents temperature boundary condition values
  - If IDEN=4: bc(j,k) should be replaced with salbc(j,k),tempbc(j,k)
- NFEN is the number of vertical levels in the model 