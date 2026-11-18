.. meta::
   :description: Flux specified boundaries in ADCIRC
   :keywords: adcirc, flux specified boundaries

.. _flux_specified_boundaries:

Flux Specified Boundaries
=========================

The table below describes the valid IBTYPE values and their corresponding boundary condition types. See also :ref:`IBTYPE <IBTYPE>` for full details of each boundary type.

.. note::

   * Boundary types 0-5 and 102 above can introduce instabilities into the ADCIRC solution; types 20-25 and 122 are therefore preferred.
   * Avoid boundary types 10-13 and 112 unless very high mesh resolution is provided to allow resolution of the lateral boundary layer.
   * Boundary types 20-25 and 122 are preferred over the corresponding essential boundary conditions (types 0-5 and 102).

.. list-table:: Normal Flux Boundary Type Values (IBTYPE)
   :widths: 10 10 10 10 10 30
   :width: 100%
   :header-rows: 1
   :class: wrap-table, tighter-table

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
   Abbreviations in column headers:

   * **Num. Impl.** - Numerical Implementation
   * **Tang. Slip** - Tangential Slip

.. note::
   Different boundary types require different parameters to be specified in the fort.14 file:
   
   * Basic boundaries (IBTYPE = 0, 1, 2, 10, 11, 12, 20, 21, 22, 102, 112, 122): Node numbers only
   * Weir boundaries (IBTYPE = 3, 13, 23): Node numbers plus barrier height and flow coefficients
   * Interior barrier boundaries (IBTYPE = 4, 24): Node numbers, connected nodes, barrier height, and flow coefficients
   * Pipe boundaries (IBTYPE = 5, 25): Node numbers, connected nodes, barrier parameters, and pipe specifications (height, coefficient, diameter)

**Essential boundaries with free slip** are applied by specifying the
contribution (zero or non-zero) to the normal boundary flux integral in the
continuity equation and by specifying the (zero or non-zero) normal velocity in
the momentum equations. This boundary condition should satisfy the flux balance
in a global sense and the normal flux at each boundary node.

**Essential boundaries with no slip** are applied by specifying the contribution
(zero or non-zero) to the normal boundary flux integral in the continuity
equation and by setting the (zero or non-zero) normal velocity and zero
tangential velocity rather than solving momentum equations along the boundary.
This boundary condition should correctly satisfy the normal flux balance in a
global sense and zero tangential velocity at each boundary node.

**Natural boundaries** are applied by specifying the zero or non-zero
contribution to the normal boundary flux integral in the continuity equation.
There is no constraint on velocity (normal or tangential) in the momentum
equations. This boundary condition should correctly satisfy the flux balance in
a global sense but will only satisfy the normal flow at each boundary node in
the limit of infinite resolution.

.. _general_notes_for_normal_flow_boundary_conditions:

General Notes for Normal Flow Boundary Conditions
-------------------------------------------------

All external (external no normal flow, external with specified normal flow and
external barrier) boundaries should be listed in consecutive order around the
outside of the entire domain before any internal (island with no normal flow or
internal barrier) boundary segments are listed. Internal barrier boundaries that
intersect an external boundary should be specified separately, even though this
will result in some nodes being specified in both boundaries, (see below).

An external no normal flow or specified normal flow boundary that completely
surrounds the domain (e.g., a lake) should be closed by repeating the first node
as the last node.

All no normal flow internal boundaries (e.g., islands) should be closed by
repeating the first node as the last node.

Unless the boundary segment is closed, always start listing the boundary nodes
where two boundaries connect.

When an external specified normal flow or external barrier boundary connects to
an external no normal flow boundary, the initial leg of the external specified
normal flow boundary or external barrier boundary is used to determine the
normal and tangential direction at the node common to both boundaries.

External boundaries with specified (non-zero) normal flow boundary conditions
and external barrier boundaries can not connect. They must be separated by an
external no normal flow boundary or an elevation specified boundary.

An internal barrier boundary can intersect an external no normal flow boundary.
(For example a levee may project out of an external no normal flow boundary in
which case 2 nodes, the front and back node on the internal barrier boundary,
would be common to the external boundary.) However, the common external nodes
must be treated in the weak sense. ADCIRC will automatically accommodate this as
follows:

-  If the external no flow boundary is specified as essential with slip
   (IBTYPE(k)=0) and the internal barrier boundary is specified as essential
   with slip (IBTYPE(k)=4), the common external boundary nodes are automatically
   changed to natural no flow boundary nodes (IBTYPE(k)=20).
-  If the external no flow boundary is specified as essential with no slip
   (IBTYPE(k)=10) and the internal barrier boundary is specified as essential
   with slip (IBTYPE(k)=4), the common external boundary nodes are automatically
   changed to natural no flow boundary nodes (IBTYPE(k)=20).
-  If the external no flow boundary is specified as natural with slip
   (IBTYPE(k)=20) and the internal barrier boundary is specified as essential
   with slip (IBTYPE(k)=4), no changes are made.
-  If the external no flow boundary is specified as essential with slip
   (IBTYPE(k)=0) and the internal barrier boundary is specified as natural with
   slip (IBTYPE(k)=24), no changes are made.
-  If the external no flow boundary is specified as essential with no slip
   (IBTYPE(k)=10) and the internal barrier boundary is specified as natural with
   slip (IBTYPE(k)=24), the common external boundary nodes are automatically
   changed to essential no flow with slip boundary nodes (IBTYPE(k)=0).
-  If the external no flow boundary is specified as natural with slip
   (IBTYPE(k)=20) and the internal barrier boundary is specified as natural with
   slip (IBTYPE(k)=24), no changes are made.

Internal barrier boundaries can not intersect external specified flow boundary
segments, external barrier boundary segments or internal no normal flow
boundaries.

For all normal flow boundaries (i.e. IBTYPE(k) = 0,1,2,3,
4,10,11,12,13,20,21,22,23,24,30), the boundary flux integral in the continuity
equation is evaluated with the appropriate (zero, specified or computed) flux.
This is a natural boundary condition. For natural normal flow boundaries
(IBTYPE(k) = 20,21,22,23,24), this is the only lateral boundary condition that
is used.

For essential normal flow boundaries with tangential slip (IBTYPE(k) =
0,1,2,3,4,10,11,12,13), the normal direction momentum equation (obtained by
re-orienting the x/y momentum equations into normal/tangential directions) is
eliminated and the normal velocity is set by dividing the normal flux per unit
width (zero, specified, or computed) by the total water column height.

For essential normal flow boundaries with no tangential slip (IBTYPE(k) = 10,11,
12,13), both momentum equations are eliminated. The tangential velocity is set
equal to zero and the normal velocity is set by dividing the normal flux per
unit width (zero, specified, or computed) by the total water column height. Use
of this boundary condition requires considerable care since strictly speaking
this type of no slip boundary condition is only mathematically justifiable if
lateral viscous terms are used in the simulation and only physically justifiable
if the lateral boundary layers are sufficiently resolved.

.. _external_barrier_boundary_note:

External Barrier Boundary Note
------------------------------

**(IBTYPE(k) = 3, 13, 23)**

Outward flow per unit width, QN2(k,j), normal to and over an external barrier
boundary is computed as:

   **Case 1** water level below or equal to the barrier height

      QN2(k,j) = 0

   **Case 2** water level above the barrier height

      QN2(k,j) = -(2/3)*BARLANCFSP(k,j)*RBARWL*((2/3)*RBARWL*G)**0.5

      where, RBARWL = ETA1(NBVV(k,j))-BARLANHT(k,j) = water height above the
      barrier
      ETA1(NBVV(k,j)) = water level computed at the previous time step at node
      NBVV(k,j)
      This formula is given by Leendertse (Aspects of SIMSYS2D ? A System for
      Two-Dimensional Flow Computation, Rand/R-3572-USGS, 1987) and is simply
      the formula for a broad crested weir (e.g., see Henderson, Open Channel
      Flow, section 6.6).
      See also General Notes for Normal Flow Boundary Conditions

.. _internal_barrier_boundary_note:

Internal Barrier Boundary Note
------------------------------

**(IBTYPE(k) = 4, 24)**

An internal barrier boundary consists of a long thin island with parallel front
and back faces. Pairs of nodes are placed on either side of the boundary so as
to provide a one to one correspondence between the nodes on the front face and
back faces. Flow is assumed to go across the boundary from one node to its
paired node on the opposite side. The normal flow is equal in magnitude and
opposite in sign on the two sides of the boundary (e.g., outflow on the front
face = inflow on the back face). Normal flow per unit width, QN2(k,j), at
internal barrier boundary node NBVV(k,j) and its paired node IBCONN(k,j) is
computed as:

   **Case 1** water level below or equal to the barrier height on both sides of
   the barrier

      QN2(k,j) = 0

   **Case 2** water level above the barrier height but equal on both sides of
   the barrier

      QN2(k,j) = 0

   *' Case 3*' water level above the barrier height but greater on the front
   side than on the back with subcritical flow across the barrier. Subcritical
   flow from front to back across the barrier occurs if the water level height
   above the barrier on the back side is greater than 2/3 the water level height
   above the barrier on the front side (i.e., RBARWL2 > 0.667*RBARWL1).

      QN2(k,j) = -RAMP*BARINCFSB(k,j)*RBARWL2*(2*G*(RBARWL1-RBARWL2))**0.5

   **Case 4** water level above the barrier height but greater on the front side
   than on the back with supercritical flow across the barrier. Supercritical
   flow from front to back across the barrier occurs if the water level height
   above the barrier on the back side is less than or equal to 2/3 the water
   level height above the barrier on the front side (i.e., RBARWL2 <
   0.667*RBARWL1).

      QN2(k,j) = -(2/3)*RAMP*BARINCFSP(k,j)*RBARWL1*((2/3)*RBARWL1*G)**0.5

   **Case 5** water level above the barrier height but greater on the back side
   than on the front with subcritical flow across the barrier. Subcritical flow
   from back to front across the barrier occurs if the water level height above
   the barrier on the front side is greater than 2/3 the water level height
   above the barrier on the back side (i.e., RBARWL1 > 0.667*RBARWL2).

      QN2(k,j) = RAMP*BARINCFSB(k,j)*RBARWL1*(2*G*(RBARWL2-RBARWL1))**0.5

   **Case 6** water level above the barrier height but greater on the back side
   than on the front with supercritical flow across the barrier. Supercritical
   flow from back to front across the barrier occurs if the water level height
   above the barrier on the front side is less than or equal to 2/3 the water
   level height above the barrier on the back side (i.e., RBARWL1 <
   0.667*RBARWL2).

      QN2(k,j) = (2/3)*RAMP*BARINCFSP(k,j)*RBARWL2*((2/3)*RBARWL2*G)**0.5

         where

            RBARWL1 = ETA1(NBVV(k,j))-BARINHT(k,j) = water height above the
            barrier on the front side of the barrier
            RBARWL2 = ETA1(IBCONN(k,j))- BARINHT(k,j) = water height above the
            barrier on the back side of the barrier
            ETA1(NBVV(k,j)) = water level computed at the previous time step on
            the front side of the barrier
            ETA1(IBCONN(k,j)) = water level computed at the previous time step
            on the back side of the barrier

   These formulae are given by Leendertse (Aspects of SIMSYS2D ? A System for
   Two-Dimensional Flow Computation, Rand/R-3572-USGS, 1987) and are simply the
   formulae for a broad crested weir (e.g., see Henderson, Open Channel Flow,
   section 6.6).
   See also General Notes for Normal Flow Boundary Conditions

.. _internal_barrier_boundary_with_cross_barrier_pipes_note:

Vertical Element Walls Note
---------------------------

**(IBTYPE(k) = 64)**

For IBTYPE(k) = 64, the barrier is a vertical element wall. See :ref:`vertical element walls <special_features_vertical_element_walls>` for details.


Internal Barrier Boundary with Cross Barrier Pipes Note
-------------------------------------------------------

**(IBTYPE(k) = 5,25)**

This type differs from IBTYPE(k)=5 or 25 This type differs from IBTYPE(k)=4 and
24 in that the barrier contains a cross barrier pipe with a specified height,
pipe coefficient, and pipe diameter. The formulation of this boundary type is
described in detail in the Technical Publication “Leaky Internal-Barrier
Normal-Flow Boundaries in the ADCIRC Coastal Hydrodynamics Code“.

Because this boundary type has so much in common with IBTYPE(k)=4 or 24, only
the differences will be described here. Normal flow per unit width, QN2(k,j), at
internal barrier boundary with cross barrier pipes node NBVV(k,j) and its paired
node IBCONN(k,j) is computed as:

   **Case 1**: Water level on both sides of the internal barrier below the
   height of the crown of the cross barrier pipe

      QN2(k,j) = 0

   **Case 2**: Water level on both sides of the internal barrier equal

      QN2(k,j) = 0

   **Case 3:** Water elevation on the front side of the internal barrier greater
   than water elevation on the back side; water elevation on the front side
   greater than the crown height of the cross-barrier pipe; and water elevation
   on the back side below crown height of the pipe

      QN2(k,j) = -RAMP*(0.25*pi*D^2)*(2*G*RBARWL1/(1+PIPECOEFR(k,j)))^0.5

         where

            RBARWL1 = ETA2(NBVV(k,j))-PIPEHTR(k,j)

   **Case 4:** Water elevation on the front side of the internal barrier greater
   than water elevation on the back side; water elevation on the front side
   greater than the crown height of the cross barrier pipe; and water elevation
   on the back side above the crown height of the pipe.

      QN2(k,j) =
      -RAMP*0.25*pi*PIPEDIAMR^2*(2*G*(RBARWL1-RBARWL2)/PIPECOEF(k,j))^0.5

         where

            RBARWL1 = ETA2(NBVV(k,j))-PIPEHTR(k,j)
            RBARWL2 = ETA2(IBCONN(k,j))-PIPEHTR(k,j)

   **Case 5:** Water elevation on the back side of the internal barrier greater
   than water elevation on the front side; water elevation on the back side
   greater than the crown height of the cross-barrier pipe; and water elevation
   on the front side below the crown height of the pipe.

      QN2(I)= RAMP*0.25*pi*(PIPEDIAMR)^2 \* (2*G*RBARWL2 / (1
      +/PIPECOEF(k,j)))^0.5

         where

            RBARWL2 = ETA2(IBCONN(k,j))-PIPEHTR(k,j)

   **Case 6:** Water elevation on the back side of the internal barrier greater
   than water elevation on the front side; water elevation on the back side
   greater than the crown height of the cross-barrier pipe; and water elevation
   on the front side above the crown height of the pipe.

      QN2(I)=
      RAMP**0.25*pi*(PIPEDIAMR(k,j))^2*(2*G*(RBARWL2-RBARWL1)/PIPECOEFR(k,j))^0.5

         where

            RBARWL1 = ETA2(NBVV(k,j))-PIPEHTR(k,j)
            RBARWL2 = ETA2(IBCONN(k,j))-PIPEHTR(k,j)

   Details of the rationale and implementation of this boundary type can be
   found in the following reference:

      Westerink, J.J., R.A. Luettich and A. Militello, 2001, Leaky
      internal-barrier normal-flow boundaries in the ADCIRC coastal
      hydrodynamics code, Coastal and Hydraulic Engineering Technical Note
      ERDC/CHL CHETN-IV-32, U.S. Army Engineer Research and Development Center,
      Vicksburg, MS., February 2001, 28p.
