:orphan:

.. _maxmin63:

Maximum and Minimum Value Files (max*.63, min*.63)
==================================================

The following files record extreme values throughout the model domain:

* **maxele.63**: Maximum water surface elevation at each node
* **maxvel.63**: Maximum water velocity magnitude at each node
* **maxwvel.63**: Maximum wind velocity magnitude at each node
* **maxrs.63**: Maximum wave radiation stress at each node (if wave forcing is enabled)
* **minpr.63**: Minimum barometric pressure at each node

These files were developed to efficiently record extreme values at each node in the domain without requiring high-frequency time series output. ADCIRC collects these extreme values at every time step and writes them after the simulation completes.

Hotstart Considerations
-----------------------

When using hotstart files to restart a simulation, ADCIRC will attempt to load existing extreme value files from the run directory at startup. If found, these values are used as initial conditions for the hotstarted run, ensuring continuity of extreme value tracking across multiple simulation segments.

File Structure
--------------

Since ADCIRC version 51, these files contain two full domain datasets:
1. The extreme value at each node
2. The time (in seconds since cold start) when each extreme value occurred

The basic structure is shown below:

.. parsed-literal::

    :ref:`RUNDES <RUNDES>`, :ref:`RUNID <RUNID>`, :ref:`AGRID <AGRID>`

    2, :ref:`NP <NP>`, :ref:`DTDP <DTDP>` * :ref:`NSPOOLGE <NSPOOLGE>`, :ref:`NSPOOLGE <NSPOOLGE>`, :ref:`IRTYPE <IRTYPE>`

    :ref:`TIME <TIME>`, :ref:`IT <IT>`

    for k=1, :ref:`NP <NP>`
        k, extreme_value(k)
    end k loop

    :ref:`TIME <TIME>`, :ref:`IT <IT>`

    for k=1, :ref:`NP <NP>`
        k, time_of_extreme(k)
    end k loop

Notes
-----

* Output format (ASCII/binary) is determined by :ref:`NOUTGE <NOUTGE>` in the fort.15 file
* For binary output, the node number (k) is not included in the output
* Prior to ADCIRC version 51, only the extreme values were included (without timing information)
* Files are written only at the end of the simulation
* Values are tracked continuously throughout the simulation at every time step
* Extreme values from previous runs are considered when using hotstart
* These files help avoid issues with:
    * Large file sizes from full time series output
    * Limited spatial coverage from station-based output
    * Missing peak values due to output frequency limitations 