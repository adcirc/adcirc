.. meta::
   :description: Tips and Tricks in ADCIRC
   :keywords: adcirc, tips and tricks

ADCIRC+SWAN Tips and Tricks
===========================

.. _writing_alternative_station_swan_output:

Writing Alternative Station SWAN Output
---------------------------------------

| SWAN has expansive outputting capabilities, documented in `its
  manual <https://swanmodel.sourceforge.io/online_doc/swanuse/swanuse.html>`__,
  but figuring out how can be a bit difficult. Here's a simple summary of how to
  output various parameters at a user-specified set of locations. Add two lines
  to your :ref:`fort.26 <fort26>` (SWAN control) file. I put these lines after the
  ``NUM`` line in the fort.26 file. The first line declares the existence and
  name of a set of points at which you want outputs:
| ``POINTS 'aStringIdForPointsFile' FILE 'theFileNameOfYourPointsFile'``
| where ``'aStringIdForPointsFile`` is the name you're giving the list of
  locations, and ``'theFileNameOfYourPointsFile`` is the actual file name. This
  file should be placed in your base run directory. The file format is ASCII x-y
  coordinate pairs, no header line, space-delimited, and one point per line,
  e.g.
| x1 y1
| x2 y2
| x3 y3
| ...

| The second line describes the outputs you want:
| ``TABLE 'aStringIdForPointsFile' HEADER 'outputFile' HSIGN HSWELL DIR PDIR TDIR TM02 STEEPNESS OUT 20010101.000000 15 MIN``
| where ``HEADER`` (or ``NOHEADER``) indicates whether you want a header,
  ``'outputFile'`` is the output file name that SWAN writes to,
  ``HSIGN HSWELL...STEEPNESS`` are a list of the parameters to-be-output,
  ``OUT`` denotes the termination of the list of parameters, ``20010101.000000``
  indicates when outputting should start in YYYYMMDD.hhmmss format, and
  ``15 MIN`` indicates a 15-minute output interval. Some notes:

-  The header can be nice but it can also sometimes cause format problems if the
   numeric values being output are too large to fit in the fixed-width columns.
   If you do ``NOHEADER``, the file won't have this issue, though it can be
   harder to inspect visually...it'll be just as easy to read it with a computer
   program or something like Excel.
-  The list of parameters you can output is long, and the definitions are even
   longer. The list is currently on the SWAN website
   `here <https://swanmodel.sourceforge.io/online_doc/swanuse/node32.html>`__
   under the description of the ``BLOCK`` output command, and has `a
   link <https://swanmodel.sourceforge.io/online_doc/swanuse/node35.html#app:defvar>`__
   to the parameter definitions. If that link doesn't work, maybe go to the SWAN
   user manual `base
   page <https://swanmodel.sourceforge.io/online_doc/swanuse/swanuse.html>`__,
   then click on "Write or plot computed quantities".
-  In a parallel run, for model versions later than roughly ADCIRC version v52,
   the output should be a single file in the main run directory. In earlier
   versions, the output was written individually by each subdomain to its PE\*
   folder.

Note that you can add multiple ``POINTS`` and/or ``TABLE`` lines for different
list of coordinates and/or different outputs.

.. _writing_alternative_global_swan_output:

Writing Alternative Global SWAN Output
--------------------------------------

Have you ever wanted to write different SWAN parameters to output in a coupled
ADCIRC+SWAN simulation? If so, these tips and tricks are for you! The current
implementations of the dynamically coupled version of ADCIRC+SWAN only write a
few SWAN parameters to global time-series output files. Those output parameters
and files include:

-  spectrally significant wave height (filename: swan_HS.63)
-  "smoothed" peak wave period (filename: swan_TPS.63)
-  mean absolute wave period (filename: swan_TMM10.63)
-  mean wave direction (filename: swan_DIR.63)
-  and perhaps others like TM01, TM02, and WIND

The SWAN model has many other potentially useful output parameters and, in some
cases, they may be more appropriate for your needs than these default outputs.
For example, what if you preferred seeing the peak wave direction instead of the
mean wave direction? Or maybe you need the swell wave height instead of, or in
addition to, the significant wave height? The two methods below will allow you
to accomplish just this. Granted, this does require that you modify some files
within the SWAN and ADCIRC source codes, and that you recompile the source code
into your desired executables. We hope the instructions below make the task of
modifying the source code a little less daunting. There is one major caveat to
the two methods outlined below: it will not allow you to 'define' your desired
SWAN output in either the fort.26 file OR in namelist form at the bottom of
fort.15. In other words, you will receive whatever outputs you hardwire into
your source code and you cannot easily turn them on or off without again
modifying the source code.

.. _method_1___very_easy:

Method 1 - Very Easy
~~~~~~~~~~~~~~~~~~~~

Use the approach described in the `last
section <Tips_and_Tricks#Writing_Alternative_Station_SWAN_Output>`__ and use the
coordinates of all nodes as your list of stations. It's possible this may not
exactly match what you'd get with other methods, e.g. if SWAN interpolates the
station value even when a point exactly matches a node location then it would be
subject to floating point errors and to effects of wetting/drying at the other
nodes in the triangle.

.. _method_2___easy_but_limited:

Method 2 - Easy but Limited
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is by far the simplest method for modifying SWAN output in a coupled
ADCIRC+SWAN simulation. By following the instructions below, you will be able to
write any SWAN output you desire… but only to the existing SWAN output files.
Existing SWAN output files include swan_HS.63, swan_TPS.63, swan_TMM10.63,
swan_DIR.63, etc. This method will not produce NEW or ADDITIONAL SWAN output
files. Rather, it will simply write the user-specified output to these existing
filenames.

Go to the swan/ folder in your adcirc source code directory tree and open the
file swanout1.f

You are looking for a variable assignment named IVTYPE in the subroutine named
SWOEXA. The IVTYPE value is linked to SWAN output variables. For example, this
is what the first IVTYPE definition looks like in my version of SWAN…

.. code-block:: none

   !
   !       significant wave height
   !
   IVTYPE = 10

Keep searching and note the IVTYPE values corresponding to your desired SWAN
output variable(s). I have included a few below for reference.

| **IVTYPE value / variable name**
| 10 / significant wave height
| 44 / swell wave height
| 12 / peak wave period
| 14 / peak wave direction
| 13 / mean wave direction
| 16 / directional spread
| 15 / energy transport direction
| 17 / average wave length
| 18 / steepness
| With the IVTYPE value(s) noted, proceed to the src/ directory and open the
  file couple2swan.F

Search for the variable IVTYPE. Reassign the IVTYPE values using the values of
your desired output. For example, if you want to write the peak wave direction
to the “swan_DIR.63” file instead of the mean wave direction, you’d want to
change IVTYPE = 13 to IVTYPE = 14. Note that you’ll need to make this change in
approximately five places within the couple2swan.F file. The embedded images
below show comparisons between the original and modified files for this example
(substituting IVTYPE=14 (PDIR) for IVTYPE=13 (DIR)). The change areas appear in
the following subroutines of couple2swan.F:

-  ComputeWaveFrictionProperties
-  SwanOutput

You would of course want to make these IVTYPE reassignments for any and all
existing SWAN output that you'd like to modify. Keep in mind that you only have
four existing files to work with here, so that is ultimately your limit on
output parameters as well. You may want to consider renaming your output files
so that you do not misinterpret your results. The output filenames are found in
src/write_output.F. Using the example IVTYPE reassignment above, you might
consider renaming swan_DIR.63 to swan_PDIR.63. If you need more than four output
parameters and also want to change the filenames, see Method 2 below.

method1_1.png|Change area 1 of 5 (click to enlarge) method1_2.png|Change area 2
of 5 (click to enlarge) method1_3.png|Change area 3 of 5 (click to enlarge)
method1_4.png|Change area 4 of 5 (click to enlarge) method1_5.png|Change area 5
of 5 (click to enlarge)

.. _method_3___moderately_difficult_but_limitless:

Method 3 - Moderately Difficult but Limitless
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This method is more complex but it will allow you to name and create any number
of original SWAN output files. For example, say I wanted to keep the existing
default SWAN output files (referenced above in Method 1) but also write the peak
wave direction output to a brand new file named swan_PDIR.63. You can do that if
you follow the steps below. Note that the example below will also generate a
corresponding swan_PDIR_max.63 file in a manner similar to swan_HS_max.63,
swan_TPS_max.63, etc.

In this more complicated (really just more work) method you will need to perform
four general tasks:

#. Increment an output counter in multiple places (in my example, from 7 to 8
   because I am only adding the PDIR output)
#. Insert new sections of code pertaining to your desired output
#. Create and insert new variable names in GLOBAL and OUTPUT blocks
#. Add new code to direct your parameters to output files

I honestly did this somewhat blindly the first time and it worked just fine (for
me). You do need to be very careful and take your time going through the
couple2swan.F file and all internal subroutines to make sure you are making all
appropriate additions. Like method 1 above, your only modifications will be in
couple2swan.F, but you will need to refer to the IVTYPE values found in
swan/swanout1.f as described above.

As was done above, I am providing screenshots of a "diff" comparison between the
original and modified couple2swan.F files for this simple example of adding the
PDIR output as a new output file. For consistency purposes I do specify PDIR as
an output variable in my fort.26 file, but as Casey and others have mentioned in
the past I don't think that makes any difference whatsoever. If you follow my
example below you would get the PDIR outputs even if you omitted PDIR from the
fort.26 file. The change areas appear in the following subroutines of
couple2swan.F:

-  ComputeWaveFrictionProperties
-  SwanOutput

method2_1.png|Change area 1 of 7 (click to enlarge) method2_2.png|Change area 2
of 7 (click to enlarge) method2_3.png|Change area 3 of 7 (click to enlarge)
method2_4.png|Change area 4 of 7 (click to enlarge) method2_5.png|Change area 5
of 7 (click to enlarge) method2_6.png|Change area 6 of 7 (click to enlarge)
method2_7.png|Change area 7 of 7 (click to enlarge)

