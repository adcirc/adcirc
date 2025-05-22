Time Varying Weirs Schedule File
================================

The Schedule File is used in conjunction with the Time Varying Weirs Input File when VARYTYPE=3 is specified. The schedule file defines the timing and magnitude of weir height changes for all weir nodes that reference it.

File Structure
--------------

The first line of the schedule file specifies nSections, the number of periods in the schedule where weir height should change. Each subsequent line specifies one section, which includes the starting and ending time of the period, and ZF, the weir height at the end of the period. All weir heights are interpolated linearly in time during the change in height.

Time Parameters
---------------

Use some combination of the following time parameters on each line to specify the start and end of the weir height change period, relative to the start of the schedule:

Start Time Parameters:
- TimeStartDay – number of days since the beginning of the schedule
- TimeStartHour – number of hours since the beginning of the schedule
- TimeStartMin – number of minutes since the beginning of the schedule
- TimeStartSec – number of seconds since the beginning of the schedule

End Time Parameters:
- TimeEndDay – number of days since the beginning of the schedule
- TimeEndHour – number of hours since the beginning of the schedule
- TimeEndMin – number of minutes since the beginning of the schedule
- TimeEndSec – number of seconds since the beginning of the schedule

Weir Height Parameters
----------------------

- ZF – The weir elevation at the end of this change

Special ZF Values:
- -99990 – Make the final weir height equal to the mesh bathy/topo elevation at that node
- -99991 – Add a specified amount (Delta) to the weir height
- -99992 – Subtract a specified amount (Delta) from the weir height
- -99993 – Make the final weir height equal to the original weir height as specified in the fort.14 file

Additional Parameter:
- Delta – Value to add/subtract from the current elevation when using -99991 or -99992 for ZF

Notes
-----

- The entire schedule file applies to all weir nodes that reference it
- Weir heights are interpolated linearly in time during changes
- Time parameters are additive (e.g., TimeStartDay=1 and TimeStartHour=12 represents 1.5 days)
- The schedule can be repeated based on the Loop and NLoops parameters in the Time Varying Weirs Input File 