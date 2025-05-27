.. meta::
   :description: Dynamic water level correction in ADCIRC
   :keywords: adcirc, dynamic water level correction

.. _dynamic_water_level_correction:

Dynamic water level correction
==============================

**Dynamic Water Level Correction** is a process by which modeled water levels
are dynamically adjusted by use of a forcing term. The correction can be applied
as constant or varying in space and/or time. The correction is applied as a
forcing term in the momentum equations whose mathematical form is equivalent to
that of an atmospheric pressure term. This means that, for gradually-varying
corrections, corrected water levels should closely follow the input correction,
though these may deviate if a correction is applied very quickly or to an area
that has a very weak connection to an open boundary through which water can
flow. Further discussion is below in the :ref:`FAQ <dwlc_frequently_asked_questions>`.

Overviews and examples of this capability have been provided in multiple
presentations (Luettich et al. 2017 [1]_, Asher et al. 2018 [2]_) and a journal
article (Asher et al. 2019 [3]_) with details and an application to Hurricane
Matthew. Users looking for ways to generate water level correction surfaces can
look to that same article and this digital publication/data repository [4]_,
which holds the code base used in the aforementioned paper.

Version
-------

Early versions of this feature were implemented in v53, however important bug
fixes and a renaming of related variables [Note1]_ were done later, so users are
strongly encouraged to not use these earlier versions.

.. _controlling_water_level_correction:

Controlling Water Level Correction
----------------------------------

| The water level correction feature is triggered by the presence of the
  :ref:`&dynamicWaterLevelCorrectionControl <fort15>` namelist at the bottom of the :ref:`fort.15
  file <fort15>`. Here is an example of how this line is used:

.. code-block:: text

   &dynamicWaterLevelCorrectionControl dynamicWaterLevelCorrectionFileName="hsofs_offset_test_ft.dat" dynamicWaterLevelCorrectionMultiplier=0.3048, dynamicWaterLevelCorrectionRampStart=0.0, dynamicWaterLevelCorrectionRampEnd=259200.0 dynamicWaterLevelCorrectionRampReferenceTime="hotstart", dynamicWaterLevelCorrectionSkipSnaps=0 /

-  :ref:`dynamicWaterLevelCorrectionMultiplier <fort15>` allows one to easily change the
   entire dataset by a single factor, e.g. if the data were provided in feet;
   ADCIRC requires meters.
-  :ref:`dynamicWaterLevelCorrectionRampStart <fort15>` and
   :ref:`dynamicWaterLevelCorrectionRampEnd <fort15>` are the time in seconds for a linear
   ramp up period.
-  :ref:`dynamicWaterLevelCorrectionRampReferenceTime <fort15>` controls whether the ramp
   times are relative to the hotstart time or the coldstart time. Prior to the
   start of the ramp, the correction values are multiplied by 0.0, eliminating
   the correction until the ramp period starts.
-  :ref:`dynamicWaterLevelCorrectionSkipSnaps <fort15>` is used to skip the specified number
   of initial data sets in the corrections data file. It can also be used to
   insert blank snaps at the beginning of the run if it is negative.

.. _water_level_correction_input_file:

Water Level Correction Input File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ADCIRC will begin reading this input file when the run starts, so time snaps are
needed from the start of the run. ADCIRC looks for the water level correction
input file in the same directory as the other full-domain files, both in serial
and in parallel. In parallel, each subdomain does domain decomposition of the
correction data on the fly, so the only change to adcprep is the read/write of
the new namelist to the subdomain fort.15 files. The water level correction
input file starts with the following three header lines:

.. code-block:: text

   # test file for 1ft correction on hsofs mesh during matthew
   14400.0        # time increment between datasets, must be positive -- read but ignored if constant correction
   0.0            # default correction value

The first line is a comment line. The second line is the time interval in
seconds between successive datasets. The third line is the default correction
value; this allows us to create correction input files in a "sparse" format
where only the non-default values are provided.If there is only one dataset in
the correction file, then the values are treated as temporally constant. If
there is more than one dataset, ADCIRC linearly interpolates between them in
time. If ADCIRC runs out of correction data, it continues with the last dataset
as a set of constant values for the remainder of the run.

Because each correction dataset may have a different number of values, the #
(hash, pound, or number) symbol is used as a separator. Specifically, users
should place a line containing # in the first and second columns, i.e. ``##``,
of the file after the end of the current time snap of data. When ADCIRC sees
this, it understands the next line is the start of the next time snap of data.

Output
~~~~~~

If a water level correction input file is specified, then the time interpolated
water level correction values will be written to a file called
`offset.63 <offset.63>`__ on the same schedule as output to the
`fort.63 <fort.63>`__. This file can then be used as a diagnostic to check that
the correction values have been applied as expected.

Compatibility
-------------

Although the correction is implemented as a pseudo barometric pressure, it is
entirely independent of any/all meteorological forcing and meteorological data.
The values do not show up in the fort.71 (barometric pressure stations) or
fort.73 (fulldomain barometric pressure field) output files. The correction can
be used with any meteorological forcing data, or with no meteorological forcing
data at all.

Example
-------

An ADCIRC model is coldstarted with a tide only run that uses a 10 day ramp, and
runs for 15 days. How should the water level correction be configured?

In this case, I suggest that you apply your bias correction after reaching full
strength tidal forcing (i.e., after the 10 day tidal ramp). So to have your bias
correction ramp completed at the end of the 15 day tidal spinup, start the ramp
period at 10.0 days and end it at 13.0 days. This should ensure you clearly see
the correction occurring when you look at the output data, while avoiding
applying the correction too rapidly. In summary, set your ramp period in the
fort.15 file to
``dynamicWaterLevelCorrectionRampStart =864000.0, dynamicWaterLevelCorrectionRampEnd = 1036800.0``

Since this is the model spin-up period, using a time-constant correction is
likely sufficient. That means a single time snap of correction data can be
supplied, for example:

.. code-block:: text

   # bias correction values test for GFS 20161002 cycle 00Z
   99999.9        # time increment in seconds; not used in this case
   0.0            # default nodal value applied to any node not specified below
   234    0.5
   1223  0.6
   1789  -0.2
   4000  0.1
   ##

| After the spinup, in preparation for a forecast run, when you want to migrate
  from one set of bias correction values to another over a 6 hour time period,
  the fort.15 file's ``dynamicWaterLevelCorrectionControl`` line would be
  something like this:

.. code-block:: text

   &dynamicWaterLevelCorrectionControl dynamicWaterLevelCorrectionFileName= "offset_migration_from_00Z_to_06Z.dat", dynamicWaterLevelCorrectionMultiplier = 1.0, dynamicWaterLevelCorrectionRampStart = 1252800.0, dynamicWaterLevelCorrectionRampEnd = 1296000.0, dynamicWaterLevelCorrectionRampReferenceTime = "coldstart", dynamicWaterLevelCorrectionSkipSnaps = 0 /

and the corresponding file would look like something like the following:

.. code-block:: text

   # bias correction values for GFS 20161002 from cycle 00Z to cycle 06Z
   21600.0        # 6 hour time increment in seconds
   0.0            # default nodal value applied to any node not specified below
   234    0.5
   1223  0.6
   1789  -0.2
   4000  0.1
   ##
   512   0.3
   1001  -0.5
   2346  0.74
   4000  -0.1
   ##

At the end of the 15.0 day tidal spinup in this example the bias correction at
node 234 would be 0.5. Three hours later (after the hotstart) it would be 0.25.
Six hours after the hotstart it would be zero.

At the end of the 15.0 day tidal spinup in this example the bias correction at
node 1001 would be 0.0. Three hours later (after the hotstart) it would be
-0.25. Six hours after the hotstart it would be -0.5.

At the end of the 15.0 day tidal spinup in this example the bias correction at
node 4000 would be 0.1. Three hours later (after the hotstart) it would be 0.0.
Six hours after the hotstart it would be -0.1.

.. _dwlc_frequently_asked_questions:

Frequently Asked Questions
--------------------------

Question:
   When writing the &dynamicWaterLevelCorrectionControl namelist,
   should I write the parameters on a single line, comma separated?

Answer:
   I always put fortran namelist parameters on a single line and
   separate the assignments with commas, but other formatting may be possible
   depending on the compiler you are using. Commas has always worked for me, but
   I've not tried any other way.

Question:
   Do I need to ramp the correction, i.e., is it wise to do so for
   some reason?

Answer:
   It should behave the same as a sudden application of atmospheric
   pressure forcing, which could lead to (presumably spurious) waves. Our shortest
   ramp up period is 12 hours but you might get away with 6 hours or even less. We
   have seen that rapidly ramping in the correction can lead to two outcomes that
   may be undesirable. One is a strong geostrophic current response. The other is
   that the water level response can lag behind the correction. This is most common
   for large water bodies with narrow connections to the open ocean, such as
   Pamlico Sound in North Carolina. Further details are in Asher et al. 2019 [3]_.

Question:
   Does the application of a dynamic water level correction surface
   affect the velocity solution (not just the elevation solution as it is primarily
   intended)?

Answer:
   Yes, hopefully in a physically realistic manner. As an example,
   when applying a water level correction to a back bay area with a nearby tidal
   inlet, the water level change in the back bay must be reflected in the velocity
   in the inlet in order to maintain continuity (conservation of mass). In the case
   of an initial application of a temporally steady water level correction, this
   velocity effect should be a one-time occurrence. On the other hand, for a time
   varying water level correction surface, the continuous adjustments to the water
   level will be continuously reflected in the velocity solution. This should be
   physically realistic when the water *should have* come from offshore, but if the
   elevated water level is due to rainfall, then the model is directing flow in the
   wrong direction. Further details are in Asher et al. 2019 [3]_.

Notes
-----

.. [Note1]
   If using an older version, note that in the fort.15 file, instances of
   "dynamicWaterLevelCorrection" should be changed to "offset" in the namelist
   parameter names, e.g. "offsetControl".

References
----------

.. raw:: html

   <references />

.. [1]
   Luettich, R.L., T.G. Asher, B.O. Blanton, J.G. Fleming. Representing Low
   Frequency, Spatially Varying Water Level Anomalies in Storm Surge
   Computations. 2017 American Meteorological Society Annual Meeting. `Link to
   talk <https://ams.confex.com/ams/97Annual/webprogram/Paper316033.html>`__

.. [2]
   Asher, T.G., R.L. Luettich, J.G. Fleming, B.O.Blanton. Assimilation of
   Observed Water Levels into Storm Surge Model Predictions. 2018 American
   Meteorological Society Annual Meeting. `Link to
   talk <https://ams.confex.com/ams/98Annual/webprogram/Paper334044.html>`__

.. [3]
   Asher, T.G., Luettich Jr., R.A., Fleming, J.G., Blanton, B.O., 2019. Low
   frequency water level correction in storm surge models using data
   assimilation. Ocean Modelling 144, 101483.
   https://doi.org/10.1016/j.ocemod.2019.101483

.. [4]
   Asher, T., 2019. Hurricane Matthew (2016) Storm Surge and Wave Simulations
   with Data Assimilation. https://doi.org/10.17603/2Z8H-7K90

