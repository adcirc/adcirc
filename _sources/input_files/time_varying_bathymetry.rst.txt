Fort.141: Time Varying Bathymetry Input File
============================================

The ADCIRC Time Varying Bathymetry Input File (``fort.141``) is used by ADCIRC whenever the :ref:`NDDT <NDDT>` value in the ``fort.14`` is nonzero. There are two variations of this file, depending on the value of :ref:`NDDT <NDDT>`.

Full Domain Bathymetry Change
-----------------------------

If :ref:`NDDT <NDDT>` is +/-1, then each bathymetry dataset in the file covers the full domain, and the file format is as follows:

.. parsed-literal::

    for i=1 to numDataSets
        for j=1 to :ref:`NP <NP>`
            j, depth(j)
        end j loop
    end i loop

where:
- :ref:`NP <NP>` is the number of nodes in the horizontal mesh
- j is the node number
- depth(j) has the same meaning as DP in the mesh file (``fort.14``)

Limited Area Bathymetry Change
------------------------------

If :ref:`NDDT <NDDT>` is +/-2, then each dataset in the ``fort.141`` file covers only part of the domain, and the file format is as follows:

.. parsed-literal::

    for i=1 to numDataSets
        "#"
        for j=1 to areaNodes
            j, depth(j)
        end j loop
    end i loop

Notes
-----

- This file is only required when :ref:`NDDT <NDDT>` is nonzero
- For NDDT=+/-1:
  - Each dataset must cover all nodes in the domain
  - The number of records per dataset is fixed (equal to NP)
- For NDDT=+/-2:
  - Each dataset can cover a different subset of nodes
  - The separation between datasets is achieved by placing a hash mark ("#") in the second column of a line
  - This allows for simulations where the number of nodes that change their bathymetry varies over time
- The depth values have the same meaning as the DP values in the mesh file (``fort.14``) 