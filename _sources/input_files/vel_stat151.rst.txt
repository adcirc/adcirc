Vel_stat.151: Velocity Station Location Input File
==================================================

The reading of the vel_stat.151 file is triggered when the number of velocity recording stations (:ref:`NSTAV`) in the :doc:`fort.15 <fort15>` file is set to a negative value.

File Structure
--------------

.. parsed-literal::

    :ref:`NSTAV2`
    for k=1, :ref:`NSTAV2`
        :ref:`XEV(k) <XEV>`, :ref:`YEV(k) <YEV>`
    end k loop

Notes
-----

If the value of NSTAV2 differs from the value of NSTAV (as read from the fort.15 file) the value of NSTAV2 will be used. If there are fewer than NSTAV2 stations listed in the vel_stat.151 file, ADCIRC will stop with an error. If there are more than NSTAV2 stations listed in the vel_stat.151 file, only the first NSTAV2 of them will be used.
