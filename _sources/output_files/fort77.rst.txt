Fort.77: Time-varying Weir Output
=================================

This file contains time-varying weir elevation data when time-varying weirs are activated and output is requested in ADCIRC. The file records changes in elevation from the original elevations specified in the :doc:`Grid and Boundary Information File <../input_files/fort14>`. When written in sparse format, the data is recorded only for weir nodes, with node numbers referenced to the global node numbering system.

File Structure
--------------

The basic structure is shown below:

.. parsed-literal::

    :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`

    :ref:`NDSETSE <NDSETSE>`, :ref:`NP <NP>`, :ref:`DTDP <DTDP>` * :ref:`NSPOOL_TVW <NSPOOL_TVW>`, :ref:`NSPOOL_TVW <NSPOOL_TVW>`, :ref:`IRTYPE <IRTYPE>`

    :ref:`TIME <TIME>`, :ref:`IT <IT>`

    for k=1, :ref:`NP <NP>`
        k, :ref:`TVW(k) <TVW>`
    end k loop

Notes
-----

* Only elevation changes relative to the original fort.14 elevations are recorded
* When using sparse format, only weir nodes are included in the output
* Node numbers (k) reference the global node numbering system
* The output frequency is controlled by :ref:`NSPOOL_TVW <NSPOOL_TVW>` parameter
* TVW(k) represents the change in elevation at node k from the original fort.14 elevation 