.. _fort19:

Fort.19: Non-periodic Elevation Boundary Condition File
=======================================================

The fort.19 file contains non-periodic, time varying elevation boundary conditions for "elevation specified" boundary nodes. This file is only read when an "elevation specified" boundary condition has been specified in the :doc:`Grid and Boundary Information File <fort14>` (NOPE>0) and NBFR=0 in the :doc:`Model Parameter and Periodic Boundary Condition File <fort15>`.

File Structure
--------------

The basic file structure is shown below. Each line of input data is represented by a line containing the input variable name(s). Loops indicate multiple lines of input.

.. parsed-literal::

   :ref:`ETIMINC <ETIMINC>`
   for k=1 to :ref:`NETA <NETA>`
      :ref:`ESBIN(k) <ESBIN>`
   end k loop

   Repeat the block above for each time step until the end of the simulation.

Notes
-----

* The first set of elevation values are provided at TIME=STATIM (as specified in fort.15). Additional sets of elevation values are provided every ETIMINC.
* Enough sets of elevation values must be provided to extend for the entire model run, otherwise the run will crash!
* The node order in this file must match the order specified in the elevation specified boundary condition part of the Grid and Boundary Information File.

Example
-------

The following is a simple example of a fort.19 file with three elevation nodes:

.. code-block:: none

   900.0
   0.5
   0.5
   0.5
   0.6
   0.6
   0.6
   0.7
   0.7
   0.7

This example shows:

* A time increment (ETIMINC) of 900.0 seconds (15 minutes)
* Three elevation specified boundary nodes (NETA=3)
* Three time steps of data, with values starting at 0.5 meters and increasing by 0.1 meters for each time step
