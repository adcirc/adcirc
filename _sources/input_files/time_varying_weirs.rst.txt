Fort.142: Time Varying Weirs Input File
=======================================

The Time Varying Weirs Input File contains one line for each weir node (or node pair, for internal barrier boundaries) that will change its elevation during the course of the simulation.

File Structure
--------------

The first line of the file contains the number of time varying weir nodes (or node pairs) that will appear in the file. Each subsequent line contains a comma separated list of variables that describes how the elevation of that weir node (or node pair) should change during the simulation. The variables on a particular line can appear in any order.

Required Variables
------------------

The following variables are required for each line:

- X1 – x-coordinate location of the first node as part of this weir
- Y1 – y-coordinate location of the first node as part of this weir

For internal boundaries (type 24), two additional parameters are required:

- X2 – x-coordinate location of the second node as part of this weir
- Y2 – y-coordinate location of the second node as part of this weir

ADCIRC searches the mesh to associate the X1, Y1 coordinates (and X2, Y2 if specified) with the locations of weir nodes. If a weir node is not found at those coordinates, the code will terminate.

Elevation Change Types
----------------------

The VARYTYPE parameter specifies the trigger or timing of the elevation change for that node:

1. One time elevation change based upon time (VARYTYPE=1)
2. One time elevation change based upon water surface elevation (WSE) (VARYTYPE=2)
3. Schedule based elevation change (VARYTYPE=3)

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

Notes
-----

- The weir will not decrease below the topographic elevation specified at the node
- For internal weirs with two nodes, the weir height will not decrease below the topographic elevation of either node
- A hot start simulation will not have any knowledge of previous weir elevation changes
- The time varying weir input file should not assume anything happens prior to the first ADCIRC time step
- Weir elevation changes prescribed before the start of the current simulation will not be considered 