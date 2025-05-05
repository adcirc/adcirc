Met_stat.151: Meteorological Recording Station Location Input File
==================================================================

The reading of the met_stat.151 file is triggered when the number of meteorological recording stations (:ref:`NSTAM`) in the :doc:`fort.15 <fort15>` file is set to a negative value.

File Structure
--------------

.. parsed-literal::

    :ref:`NSTAM2`
    for k=1, :ref:`NSTAM2`
        :ref:`XEM(k) <XEM>`, :ref:`YEM(k) <YEM>`
    end k loop

Notes
-----

If the value of NSTAM2 differs from the value of NSTAM (as read from the fort.15 file) the value of NSTAM2 will be used. If there are fewer than NSTAM2 stations listed in the met_stat.151 file, ADCIRC will stop with an error. If there are more than NSTAM2 stations listed in the met_stat.151 file, only the first NSTAM2 of them will be used. 