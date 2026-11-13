.. _fort14:

Fort.14: Grid and Boundary Information File
===========================================

The fort.14 file contains the finite element grid, the bathymetric data, and the boundary information used by ADCIRC. This file is required to run the ADCIRC model.

File Structure
--------------

This file contains 4 sections: nodal data, element data, open boundary data, and flow boundary data. The basic file structure is shown below. Each line of input data is represented by a line containing the input variable name(s). Loops indicate multiple lines of input. Conditional input is indicated by a clause following the variable name(s).

.. parsed-literal::

   :ref:`AGRID <AGRID>`
   :ref:`NE <NE>`, :ref:`NP <NP>`
   for k=1 to :ref:`NP <NP>`
      :ref:`JN(k) <JN>`, :ref:`X(k) <X>`, :ref:`Y(k) <Y>`, :ref:`DP(k) <DP>`
   end k loop
   for k=1 to :ref:`NE <NE>`
      :ref:`JE(k) <JE>`, :ref:`NHY(k) <NHY>`, :ref:`NM(k,1) <NM>`, :ref:`NM(k,2) <NM>`, :ref:`NM(k,3) <NM>`
   end k loop
   :ref:`NOPE <NOPE>`
   :ref:`NETA <NETA>`
   for k=1 to :ref:`NOPE <NOPE>`
      :ref:`NVDLL(k) <NVDLL>`, :ref:`IBTYPEE(k) <IBTYPEE>`
      for j=1 to :ref:`NVDLL(k) <NVDLL>`
         :ref:`NBDV(k,j) <NBDV>`
      end j loop
   end k loop
   :ref:`NBOU <NBOU>`
   :ref:`NVEL <NVEL>`
   for k=1 to :ref:`NBOU <NBOU>`
      :ref:`NVELL(k) <NVELL>`, :ref:`IBTYPE(k) <IBTYPE>`
      for j=1 to :ref:`NVELL(k) <NVELL>`
         :ref:`NBVV(k,j) <NBVV>`      if IBTYPE(k) = 0, 1, 2, 10, 11, 12, 20, 21, 22, 30
         :ref:`NBVV(k,j) <NBVV>`, :ref:`BARLANHT(k,j) <BARLANHT>`, :ref:`BARLANCFSP(k,j) <BARLANCFSP>`      if IBTYPE(k) = 3, 13, 23
         :ref:`NBVV(k,j) <NBVV>`, :ref:`IBCONN(k,j) <IBCONN>`, :ref:`BARINHT(k,j) <BARINHT>`, :ref:`BARINCFSB(k,j) <BARINCFSB>`, :ref:`BARINCFSP(k,j) <BARINCFSP>`      if IBTYPE(k) = 4, 24, 64
         :ref:`NBVV(k,j) <NBVV>`, :ref:`IBCONN(k,j) <IBCONN>`, :ref:`BARINHT(k,j) <BARINHT>`, :ref:`BARINCFSB(k,j) <BARINCFSB>`, :ref:`BARINCFSP(k,j) <BARINCFSP>`, :ref:`PIPEHT(k,j) <PIPEHT>`, :ref:`PIPECOEF(k,j) <PIPECOEF>`, :ref:`PIPEDIAM(k,j) <PIPEDIAM>`      if IBTYPE(k) = 5, 25
      end j loop
   end k loop

Open Boundary Data
------------------

See :ref:`open boundaries <open_boundaries>` for details of the open boundary data.


Normal Flux Boundary Types
--------------------------

See :ref:`flux specified boundaries <flux_specified_boundaries>` for details of the normal flux boundary types. See also :ref:`IBTYPE <IBTYPE>` for full details of each boundary type.


Example
-------

The following is a simple example of a fort.14 file for a small domain with 3 elements and 4 nodes with an open boundary and no flow boundaries:

.. code-block:: none

   Simple ADCIRC domain
   3 4    ! NE NP
   1 0.0 0.0 -10.0 ! JN X Y DP
   2 1.0 0.0 -10.0
   3 1.0 1.0 -10.0
   4 0.0 1.0 -10.0
   1 3 1 2 3 ! NE NHY NM(1,1) NM(1,2) NM(1,3)
   2 3 1 3 4
   3 3 2 3 1
   1 ! NOPE
   4 ! NETA
   4 0 ! NVDLL(1) IBTYPEE(1)
   1 1 2 ! NBDV(1,1)
   2 1 3 ! NBDV(1,2)
   3 2 3 ! NBDV(1,3)
   4 3 1 ! NBDV(1,4)
   0     ! NBOU