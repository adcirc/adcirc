Fort.35: Level of No Motion Boundary Condition Input File
=========================================================

The Level of No Motion Boundary Condition Input File (``fort.35``) is used in 3D baroclinic simulations when the :ref:`BCFLAG_LNM <BCFLAG_LNM>` parameter is set to 1 in the ``fort.15`` file. This file specifies elevation changes for ocean boundary nodes.

File Structure
--------------

The file format is as follows:

.. parsed-literal::

    for i=1 to numberOfDataSets
        comment line (date)
        for k=1 to number_of_ocean_boundary_nodes
            k, elevation_change
        end k loop
    end i loop

Notes
-----

- The elevation changes should be specified for all ocean boundary nodes
- Each dataset should include a comment line with the date
- The elevation changes are used to adjust the level of no motion boundary condition 