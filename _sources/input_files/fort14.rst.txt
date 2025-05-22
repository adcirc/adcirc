.. _fort14:

Fort.14: Grid and Boundary Information File
===========================================

The fort.14 file contains the finite element grid, the bathymetric data, and the boundary information used by ADCIRC. This file is required to run the ADCIRC model.

File Structure
--------------

The basic file structure is shown below. Each line of input data is represented by a line containing the input variable name(s). Loops indicate multiple lines of input. Conditional input is indicated by a clause following the variable name(s).

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

Normal Flux Boundary Types
--------------------------

The table below describes the valid IBTYPE values and their corresponding boundary condition types. See also :ref:`IBTYPE <IBTYPE>` for full details of each boundary type.

.. list-table:: Normal Flux Boundary Type Values (IBTYPE)
   :widths: 5 8 10 12 12 53
   :width: 100%
   :header-rows: 1
   :class: wrap-table, tight-table

   * - IBTYPE
     - Location
     - Flux Type
     - Num. Impl.
     - Tang. Slip
     - Typical Use
   * - 0
     - External
     - Zero
     - Essential
     - Free
     - Mainland boundaries
   * - 1
     - Internal
     - Zero
     - Essential
     - Free
     - Island boundaries
   * - 2
     - External
     - Nonzero inflow
     - Essential
     - Free
     - River or ocean inflow boundaries; 
       if flux is periodic it is specified in fort.15; 
       if it is time varying and aperiodic, it must be 
       specified in fort.20 file
   * - 3
     - External
     - Outflow
     - Essential
     - Free
     - Flow over a weir out of the domain; 
       must specify levee height; 
       ADCIRC calculates the fluxes
   * - 4
     - Internal
     - Zero or nonzero
     - Essential
     - Free
     - Interior levees; must specify levee height; 
       ADCIRC calculates the fluxes using weir formula
   * - 5
     - Internal
     - Zero or nonzero
     - Essential
     - Free
     - Interior levees with cross-barrier pipes (like a culvert); 
       must specify levee height and other parameters; 
       ADCIRC calculates the fluxes
   * - 10
     - External
     - Zero
     - Essential
     - No slip
     - As ibtype 0 above but no slip
   * - 11
     - Internal
     - Zero
     - Essential
     - No slip
     - As ibtype 1 above but no slip
   * - 12
     - External
     - Nonzero
     - Essential
     - No slip
     - As ibtype 2 above but no slip
   * - 13
     - External
     - Outflow
     - Essential
     - No slip
     - As ibtype 3 above but no slip
   * - 20
     - External
     - Zero (weak)
     - Natural
     - Free
     - As ibtype 0 but natural boundary; 
       preferred over ibtype 0 or 10
   * - 21
     - Internal
     - Zero (weak)
     - Natural
     - Free
     - As ibtype 1 but natural boundary; 
       preferred over ibtype 1 or 11
   * - 22
     - External
     - Nonzero (weak)
     - Natural
     - Free
     - As ibtype 2 but natural boundary; 
       preferred over ibtype 2 or 12
   * - 23
     - External
     - Outflow (weak)
     - Natural
     - Free
     - As ibtype 3 but natural boundary; 
       preferred over ibtype 3 or 13
   * - 24
     - Internal
     - Zero or nonzero (weak)
     - Natural
     - Free
     - As ibtype 4 but natural boundary; 
       preferred over ibtype 4
   * - 25
     - Internal
     - Zero or nonzero (weak)
     - Natural
     - Free
     - As ibtype 5 but natural boundary; 
       preferred over ibtype 5
   * - 64
     - Internal
     - Zero or nonzero (weak)
     - Natural or condensed
     - Free
     - Defines vertical element walls
   * - 102
     - External
     - Nonzero inflow
     - Essential
     - Free
     - As ibtype 2 but baroclinic instead of barotropic; 
       also requires density-related boundary conditions 
       in fort.39 input file
   * - 112
     - External
     - Nonzero
     - Essential
     - No slip
     - As ibtype 12 but baroclinic instead of barotropic; 
       also requires density-related boundary conditions 
       in fort.39 input file
   * - 122
     - External
     - Nonzero (weak)
     - Natural
     - Free
     - As ibtype 22 but baroclinic instead of barotropic; 
       also requires density-related boundary conditions 
       in fort.39 input file

.. note::
   **Abbreviations in column headers:**
   
   * **Num. Impl.** - Numerical Implementation
   * **Tang. Slip** - Tangential Slip

.. raw:: html

   <style>
   .wrap-table th, .wrap-table td {
     white-space: normal !important;
     word-wrap: break-word !important;
     max-width: 100% !important;
     overflow-wrap: break-word !important;
     hyphens: auto !important;
   }
   </style>

.. note::
   Different boundary types require different parameters to be specified in the fort.14 file:
   
   * Basic boundaries (IBTYPE = 0, 1, 2, 10, 11, 12, 20, 21, 22, 102, 112, 122): Node numbers only
   * Weir boundaries (IBTYPE = 3, 13, 23): Node numbers plus barrier height and flow coefficients
   * Interior barrier boundaries (IBTYPE = 4, 24): Node numbers, connected nodes, barrier height, and flow coefficients
   * Pipe boundaries (IBTYPE = 5, 25): Node numbers, connected nodes, barrier parameters, and pipe specifications (height, coefficient, diameter)

Example
-------

The following is a simple example of a fort.14 file for a small domain with 3 elements and 4 nodes:

.. code-block:: none

   Simple ADCIRC domain
   3 4
   1 0.0 0.0 -10.0
   2 1.0 0.0 -10.0
   3 1.0 1.0 -10.0
   4 0.0 1.0 -10.0
   1 3 1 2 3
   2 3 1 3 4
   3 3 2 3 1
   1 4
   1 0
   1
   2
   3
   4
   0 0 