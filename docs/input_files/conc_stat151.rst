Conc_stat.151: Concentration Station Location Input File
========================================================

The reading of the conc_stat.151 file is triggered when the number of concentration recording stations (:ref:`NSTAC`) in the :doc:`fort.15 <fort15>` file is set to a negative value.

File Structure
--------------

.. parsed-literal::

    :ref:`NSTAC2`
    for k=1, :ref:`NSTAC2`
        :ref:`XEC(k) <XEC>`, :ref:`YEC(k) <YEC>`
    end k loop

Notes
-----

If the value of NSTAC2 differs from the value of NSTAC (as read from the fort.15 file) the value of NSTAC2 will be used. If there are fewer than NSTAC2 stations listed in the conc_stat.151 file, ADCIRC will stop with an error. If there are more than NSTAC2 stations listed in the conc_stat.151 file, only the first NSTAC2 of them will be used.

This file is used for concentration station location input and is conditional on the model configuration. 