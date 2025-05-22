Fort.93: Ice Coverage Fields at All Nodes
=========================================

This file contains ice coverage field data at all nodes in the model grid. The file follows the same format as the :doc:`Pressure <fort73>` and :doc:`Elevation <fort63>` files. It is generated only when ice fields are activated and global meteorology output is enabled through the parameters (:ref:`NOUTGW <NOUTGW>`, :ref:`TOUTSGW <TOUTSGW>`, :ref:`TOUTFGW <TOUTFGW>`, :ref:`NSPOOLGW <NSPOOLGW>`) in the :doc:`Model Parameter and Periodic Boundary Condition File <../input_files/fort15>`.

File Structure
--------------

The file can be written in either ASCII or binary format, depending on how :ref:`NOUTGW <NOUTGW>` is set in the Model Parameter and Periodic Boundary Condition File. The basic structure is shown below:

.. parsed-literal::

    :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`
    :ref:`NDSETSW <NDSETSW>`, :ref:`NP <NP>`, :ref:`DTDP <DTDP>` * :ref:`NSPOOLGW <NSPOOLGW>`, :ref:`NSPOOLGW <NSPOOLGW>`, :ref:`IRTYPE <IRTYPE>`
    :ref:`TIME <TIME>`, :ref:`IT <IT>`
    for k=1, :ref:`NP <NP>`
        k, iceCoveragePercent(k)
    end k loop

Notes
-----

* Output format (ASCII/binary) is determined by :ref:`NOUTGW <NOUTGW>` in the fort.15 file
* For binary output, the node number (k) is not included in the output
* Time series data is recorded for every node in the model grid
* The output frequency is controlled by :ref:`NSPOOLGW <NSPOOLGW>` parameter
* Output timing can be further controlled using :ref:`TOUTSGW <TOUTSGW>` (start time) and :ref:`TOUTFGW <TOUTFGW>` (end time)
* The variable iceCoveragePercent represents the percentage of ice coverage at each node
* File structure matches the format used in fort.73 (pressure) and fort.63 (elevation) files 