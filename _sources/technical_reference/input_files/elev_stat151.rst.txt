.. _elev_stat151:

Elev_stat.151: Elevation Station Location Input File
====================================================

The reading of the elev_stat.151 file is triggered when the number of elevation recording stations (:ref:`NSTAE`) in the :doc:`fort.15 <fort15>` file is set to a negative value.

File Structure
--------------

.. parsed-literal::

    :ref:`NSTAE2`
    for k=1 to :ref:`NSTAE2`
        :ref:`XEL(k) <XEL>`, :ref:`YEL(k) <YEL>`
    end k loop

Notes
-----

If the value of NSTAE2 differs from the value of NSTAE (as read from the fort.15 file) the value of NSTAE2 will be used. If there are fewer than NSTAE2 stations listed in the elev_stat.151 file, ADCIRC will stop with an error. If there are more than NSTAE2 stations listed in the elev_stat.151 file, only the first NSTAE2 of them will be used.
