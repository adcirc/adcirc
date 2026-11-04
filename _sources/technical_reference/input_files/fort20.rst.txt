.. _fort20:

Fort.20: Non-periodic Normal Flow Boundary Condition File
=========================================================

The fort.20 file contains non-periodic, normal flow boundary conditions for "specified non-zero normal flow" boundary nodes. This file is only read when a "specified non-zero normal flow" boundary condition has been specified in the :doc:`Grid and Boundary Information File <fort14>` (IBTYPE=2, 12, or 22) and NFFR=0 in the :doc:`Model Parameter and Periodic Boundary Condition File <fort15>`.

File Structure
--------------

The basic file structure is shown below. Each line of input data is represented by a line containing the input variable name(s). Loops indicate multiple lines of input.

.. parsed-literal::

   :ref:`FTIMINC <FTIMINC>`
   for k=1 to :ref:`NFLBN <NFLBN>`
      :ref:`QNIN(k) <QNIN>`
   end k loop

   Repeat the block above for each time step until the end of the simulation.

Notes
-----

* The first set of normal flow values is provided at TIME=STATIM (as specified in fort.15). Additional sets of normal flow values are provided every FTIMINC.
* Enough sets of normal flow values must be provided to extend for the entire model run, otherwise the run will crash!
* The node order in this file must match the order specified in the normal flow boundary condition part of the Grid and Boundary Information File.
* A positive flow/unit width is into the domain and a negative flow/unit width is out of the domain.

Example
-------

The following is a simple example of a fort.20 file with two normal flow boundary nodes:

.. code-block:: none

   3600.0
   10.5
   -5.2
   12.0
   -6.3
   15.0
   -8.1

This example shows:

* A time increment (FTIMINC) of 3600.0 seconds (1 hour)
* Two normal flow boundary nodes (NFLBN=2)
* Three time steps of data, with:
   * Node 1: Flow into the domain (positive values) increasing over time
   * Node 2: Flow out of the domain (negative values) increasing in magnitude over time
