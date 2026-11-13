.. _fort11:

Fort.11: Density Initial Condition Input File
=============================================

The ``fort.11`` file is used to specify the initial density field for baroclinic runs. Baroclinic 2DDI runs are not yet supported. Baroclinic 3D runs occur if IDEN is not equal to 0.

File Structure
--------------

The file structure varies depending on the value of IDEN in the fort.15 file.

For Baroclinic 2DDI Runs:

If :ref:`IDEN <IDEN>` = 1 or -1 (Density):

.. parsed-literal::

    Header Line 1
    Header Line 2
    :ref:`NVP <NVP>`
    for k=1 to :ref:`NVP <NVP>`
        :ref:`jki <jki>`, :ref:`DASIGT(jki) <DASIGT>`
    end k loop

If :ref:`IDEN <IDEN>` = 2 or -2 (Salinity):

.. parsed-literal::

    Header Line 1
    Header Line 2
    :ref:`NVP <NVP>`
    for k=1 to :ref:`NVP <NVP>`
        :ref:`jki <jki>`, :ref:`DASALT(jki) <DASALT>`
    end k loop

If :ref:`IDEN <IDEN>` = 3 or -3 (Temperature):

.. parsed-literal::

    Header Line 1
    Header Line 2
    :ref:`NVP <NVP>`
    for k=1 to :ref:`NVP <NVP>`
        :ref:`jki <jki>`, :ref:`DATEMP(jki) <DATEMP>`
    end k loop

If :ref:`IDEN <IDEN>` = 4 or -4 (Temperature and Salinity):

.. parsed-literal::

    Header Line 1
    Header Line 2
    :ref:`NVP <NVP>`
    for k=1 to :ref:`NVP <NVP>`
        :ref:`jki <jki>`, :ref:`DATEMP(jki) <DATEMP>`, :ref:`DASALT(jki) <DASALT>`
    end k loop

For Baroclinic 3D Runs:

If :ref:`IDEN <IDEN>` = 1 or -1 (Density):

.. parsed-literal::

    Header Line 1
    Header Line 2
    NVN, NVP
    for k=1 to NVP
        for j=1 to NVN
            k, j, SIGT(NHNN,NVNN)
        end j loop
    end k loop

If :ref:`IDEN <IDEN>` = 2 or -2 (Salinity):

.. parsed-literal::

    Header Line 1
    Header Line 2
    :ref:`NVN <NVN>`, :ref:`NVP <NVP>`
    for k=1 to :ref:`NVP <NVP>`
        for j=1 to :ref:`NVN <NVN>`
            k, j, SAL(k,j)
        end j loop
    end k loop

If :ref:`IDEN <IDEN>` = 3 or -3 (Temperature):

.. parsed-literal::

    Header Line 1
    Header Line 2
    :ref:`NVN <NVN>`, :ref:`NVP <NVP>`
    for k=1 to :ref:`NVP <NVP>`
        for j=1 to :ref:`NVN <NVN>`
            k, j, :ref:`TEMP(k,j) <TEMP>`
        end j loop
    end k loop

If :ref:`IDEN <IDEN>` = 4 or -4 (Temperature and Salinity):

.. parsed-literal::

    Header Line 1
    Header Line 2
    :ref:`NVN <NVN>`, :ref:`NVP <NVP>`
    for k=1 to :ref:`NVP <NVP>`
        for j=1 to :ref:`NVN <NVN>`
            k, j, TEMP(k,j),SAL(k,j)
        end j loop
    end k loop

Notes
-----

- For 3D runs, j=1 represents the bottom layer and j=NVN represents the surface layer
- The file structure depends on the value of IDEN specified in the fort.15 file
- Baroclinic 2DDI runs are not yet supported
- Baroclinic 3D runs occur when IDEN is not equal to 0 