.. _fort142:

Fort.142: Time Varying Weirs Input File
=======================================

The Time Varying Weirs Input File contains one line for each weir node (or node pair, for internal barrier boundaries) that will change its elevation during the course of the simulation. It also describes time-varying crest :ref:`vertical element walls <special_features_vertical_element_walls>` (``VaryType=4``), which represent overflow gates (crest gates) that rise from the bed to close and lower to open.

The file is named by :ref:`TVW_file <TVW_file>` in the :ref:`TVWControl <TVWControl>` namelist of the fort.15 file.

File Structure
--------------

The first line of the file contains the number of time varying weir nodes (or node pairs) that will appear in the file. Each subsequent line contains a comma separated list of variables that describes how the elevation of that weir node (or node pair) should change during the simulation. The variables on a particular line can appear in any order.

Required Variables
------------------

The following variables are required for each line:

- X1 – x-coordinate location of the first node as part of this weir
- Y1 – y-coordinate location of the first node as part of this weir

For internal boundaries (type 24), two additional parameters are required (but not for ``VaryType=4``; see :ref:`fort142_varytype4`):

- X2 – x-coordinate location of the second node as part of this weir
- Y2 – y-coordinate location of the second node as part of this weir

ADCIRC searches the mesh to associate the X1, Y1 coordinates (and X2, Y2 if specified) with the locations of weir nodes. If a weir node is not found at those coordinates, the code will terminate.

Elevation Change Types
----------------------

The VARYTYPE parameter specifies the trigger or timing of the elevation change for that node:

1. One time elevation change based upon time (VARYTYPE=1)
2. One time elevation change based upon water surface elevation (WSE) (VARYTYPE=2)
3. Schedule based elevation change (VARYTYPE=3)
4. Crest of a vertical element wall following a time series (VARYTYPE=4); see :ref:`fort142_varytype4`

VARYTYPE=1 Parameters
---------------------

The following parameters are all optional except that some combination of them must define both a start time and an end time:

- Hot – 1 if time trigger is relative to hot start time, 0 if relative to cold start time (default: 0)
- TimeStartDay – number of days before start of weir node elevation change
- TimeStartHour – number of hours before start of weir node elevation change
- TimeStartMin – number of minutes before start of weir node elevation change
- TimeStartSec – number of seconds before start of weir node elevation change
- TimeEndDay – number of days before end of weir node elevation change
- TimeEndHour – number of hours before end of weir node elevation change
- TimeEndMin – number of minutes before end of weir node elevation change
- TimeEndSec – number of seconds before end of weir node elevation change

Additional required parameter:
- ZF – final weir elevation in meters relative to the mesh datum (positive up)

VARYTYPE=2 Parameters
---------------------

Required parameters:
- ZF – final weir elevation in meters relative to the mesh datum
- ETA_MAX – water surface elevation that triggers the simulated weir failure

Failure duration parameters (some combination must be specified):
- FailureDurationDay – Time in days for duration of failure
- FailureDurationHour – Time in hours for duration of failure
- FailureDurationMin – Time in minutes for duration of failure
- FailureDurationSec – Time in seconds for duration of failure

VARYTYPE=3 Parameters
---------------------

Required parameter:
- ScheduleFile – name of the Schedule File that describes when and how much the weir nodes should change their height

Optional parameters:
- Loop – Set to 1 to repeat schedule, 0 for no repeat (default), -1 for infinite repeat
- NLoops – Number of times the schedule should repeat
- TimeStartDay – number of days since hotstart when weir node elevation starts to change
- TimeStartHour – number of hours since hotstart when weir node elevation starts to change
- TimeStartMin – number of minutes since hotstart when weir node elevation starts to change
- TimeStartSec – number of seconds since hotstart when weir node elevation starts to change

.. _fort142_varytype4:

VARYTYPE=4: Time-Varying Crest Vertical Element Walls
------------------------------------------------------

``VaryType=4`` moves the crest of an :ref:`IBTYPE <IBTYPE>` = 64 vertical element wall (VEW) pair along a time series, to represent an overflow gate (crest gate) that rises from the bed to close and lowers to open. Underflow gates (flow under the gate) and gates that open sideways are not represented.

Each line locates one VEW pair. All lines that name the same crest table form one gate and move together.

The crest and the wall top move together. In the mesh, the wall top of a VEW pair is the shallower of its two nodes (the wall-top node; the other is the bed node), and the crest (:ref:`BARINHT <BARINHT>` in fort.14) lies a small distance δ above it, e.g. 0.001 m. ADCIRC keeps δ constant: at every time step it sets the crest of both entries of the pair to the crest elevation z\ :sub:`c`\ (t), and the depth of the wall-top node to δ − z\ :sub:`c`\ (t). See :ref:`time_varying_crest_vew` for how the water level is adjusted as the gate moves.

Parameters
~~~~~~~~~~

Required:

- VaryType – 4
- X1, Y1 – search point for the VEW pair, in the coordinates of the fort.14 file (degrees if :ref:`ICS <ICS>` = 2)
- ScheduleFile – name of the crest table file (see below), in the run (fulldomain) directory

Optional:

- SearchRadius – radius of the search in meters, measured in the model's projected coordinates also when ICS = 2. Defaults to :ref:`TVV_SEARCH_RADIUS <TVV_SEARCH_RADIUS>`, or 1e-6 m if neither is given
- Hot – 1 if the times in the crest table are relative to the hot-start time, 0 (default) if relative to cold start. All lines naming the same crest table must use the same value

X2 and Y2 must not be given: a pair is located by X1, Y1 and SearchRadius alone.

Locating the VEW Pair
~~~~~~~~~~~~~~~~~~~~~

The two nodes of a VEW pair may share the same horizontal location or be some distance apart. ADCIRC collects the IBTYPE = 64 boundary nodes within SearchRadius of (X1, Y1), adds the partner of each, and requires that they form exactly one pair:

- If no pair is within the radius, the run stops.
- If more than one pair is within the radius, the run stops and lists the pairs.

For a pair whose nodes share a location, use that location; the default radius is enough, and a larger one tolerates rounding in the coordinates. For a pair whose nodes are apart, use a point between them with a radius just over half their distance. In all cases the radius must be smaller than the distance to the nearest node of any other pair.

In parallel runs, every subdomain reads the ``VaryType=4`` lines from the fulldomain file, and the subdomains check that they located the same pair. adcprep does not distribute these lines; :ref:`use_TVW <use_TVW>` must be ``.true.`` for adcprep to process the file at all.

Crest Table File
~~~~~~~~~~~~~~~~

The first line gives the number of records N (N ≥ 1). Each of the next N lines holds a time in seconds (relative to cold start, or to the hot-start time if Hot = 1) and the crest elevation z\ :sub:`c` in meters relative to the mesh datum (positive up). Anything after the two values on a line is ignored and can be used as a comment.

- Times must increase strictly.
- The crest elevation is interpolated linearly in time and held constant before the first and after the last record.
- The crest is never lowered below the bed elevation plus δ; lower values are clamped, with a warning at startup.
- At cold start, the crest elevation at the start time must equal the fort.14 crest of every pair of the gate, to within 1e-6 m. Write the values with full precision if they must match the fort.14 values exactly.

Startup Checks
~~~~~~~~~~~~~~

The run stops if:

- the GWCE is not lumped (consistent mass matrix, :ref:`ILump <IM>` = 0), or :ref:`TAU0 <TAU0>` is negative (depth-dependent or time-varying);
- δ is not positive, or differs from :ref:`TVV_DELTA <TVV_DELTA>` when that is given;
- a :ref:`condensed_nodes` group contains a wall-top node of a gate together with other nodes, or nodes with different depths or δ;
- ``activateVEW1DChannelWetPerimeter`` is ``.true.`` (its bank elevations are fixed at startup).

Example
~~~~~~~

A wall whose top is a 10 m wide strip between x = 5000 and x = 5010, with VEW pairs on both faces (collocated nodes at y = 0 and y = 100) and the crest 0.001 m above the wall top. The fort.142 file:

.. code-block:: none

   4
   X1=5000.0, Y1=0.0,   VaryType=4, SearchRadius=1.0, ScheduleFile='gate1_crest.txt'
   X1=5000.0, Y1=100.0, VaryType=4, SearchRadius=1.0, ScheduleFile='gate1_crest.txt'
   X1=5010.0, Y1=0.0,   VaryType=4, SearchRadius=1.0, ScheduleFile='gate1_crest.txt'
   X1=5010.0, Y1=100.0, VaryType=4, SearchRadius=1.0, ScheduleFile='gate1_crest.txt'

The crest table ``gate1_crest.txt``: the gate is closed (crest at 0.501 m) until 12 hours, opens to the bed (−4.999 m, with the bed at −5.0 m) over one hour, and closes again from 60 to 61 hours:

.. code-block:: none

   6
        0.0    0.501   closed (matches the fort.14 crest)
    43200.0    0.501
    46800.0   -4.999   open: bed -5.0 plus delta 0.001
   216000.0   -4.999
   219600.0    0.501   closed again
   345600.0    0.501

and in fort.15:

.. code-block:: none

   &TVWControl USE_TVW=.TRUE., TVW_FILE='fort.142', NOUT_TVW=1, TOUTS_TVW=0.0, TOUTF_TVW=4.0, NSPOOL_TVW=300 /

ADCIRC logs one line per gate at startup (number of pairs, crest and δ ranges, wall-top and bed nodes), and one line per completed gate movement with the volume change caused by the crest and water-level update, which should be zero up to round-off.

Notes
-----

- The weir will not decrease below the topographic elevation specified at the node
- For internal weirs with two nodes, the weir height will not decrease below the topographic elevation of either node
- For VaryType 1–3, a hot start simulation will not have any knowledge of previous weir elevation changes. VaryType=4 gates are placed at their crest for the hot-start time, since the crest depends on time only
- The time varying weir input file should not assume anything happens prior to the first ADCIRC time step
- Weir elevation changes prescribed before the start of the current simulation will not be considered (VaryType 1–3)
