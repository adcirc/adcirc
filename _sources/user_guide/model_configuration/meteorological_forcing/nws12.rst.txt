.. meta::
   :description: NWS12 in ADCIRC
   :keywords: adcirc, nws12

.. _nws12:

NWS12
=====

Oceanweather Inc. (OWI) "WIN" and "PRE" file formats as supported in ADCIRC

fort.22
-------

NWS12 uses the **fort.22** file as a control file because the OWI format
supports individual files for wind and pressure inputs as well as multiple grid
overlays of varying resolutions.

*fort.22 layout*:

| :ref:`NWSET`
| :ref:`NWBS`
| :ref:`DWM`
| `<...to be expanded in new verions>`

Where:

-  **NWSET** is the number of nested grid overlays
-  **NWBS** is the number of "blank" (constant 0 m/s wind speed and 1013mb
   pressure) timesteps to insert
-  **DWM** is a wind multiplier to modify the wind speed magnitudes. Generally
   this should be set to 1 if using inputs sourced from OWI or if winds can be
   considered "marine exposure", with an averaging period between 10-min and
   3-hrs at a height of 10 meters.

NWS12 WIN/PRE text files containing wind and pressure data should be named
fort.22[1-4], as follows in order of increasing precedence in application to the
ADCIRC mesh:

-  *fort.221*: Largest scale pressure (mb) fields
-  *fort.222*: Largest scale wind vector (m/s) fields
-  *fort.223*: Finer scale pressure (mb) fields
-  *fort.224*: Finer scale wind vector (m/s) fields
-  ...to be expanded in new versions

.. _nwbs_details:

NWBS Details
~~~~~~~~~~~~

ADCIRC also has the ability to insert “blank” meteorological data if it is
started before the beginning of the wind data (for a tidal spinup run, for
example). Blank data is characterized by zero wind speed and 1013 mb of
atmospheric pressure. This capability is activated via the NWBS parameter in the
fort.22 file, as mentioned above, and it interacts with the sign of the NWS
value in order to provide full control over the relationship between the ADCIRC
start time (either hot start or cold start) and the beginning of the
meteorological data. This capability and the interaction with the sign of NWS is
described below.

If a NWBS is set to a positive number, NWBS specifies the number of blank snaps
to be inserted before any information is read from fort.22[1-4]. This is useful
when the duration of fort.22[1-4] files does not cover an entire ADCIRC run. Two
examples are as follows:

#. If NWS=+12, and if a hot start run is attempted from the 86400th time step,
   and if the first set, fort.221 and fort.222 files, starts at the 172800th
   time step, and if WTIMINC is set to 900 in UNIT 15, then you should set NWBS
   to (172800-0)/900 = 192. If the start times of the two sets are different,
   use the earlier start time, i.e., the start time of the first set.
#. If NWS=-12, and if a hot start run is attempted from the 86400th time step,
   and if the first set, fort.221 and fort.222 files, starts at the 172800th
   time step, and if WTIMINC is set to 900 in the fort.15, then you should set
   NWBS to (172800 – 86400)/900 = 96. If the start times of the two sets are
   different, use the earlier start time, i.e., the start time of the first set.

If NWBS is set to a negative number, then NWBS determines how many snaps in the
fort.22[1-4] should be skipped before the values in the files are used. This is
useful when the fort.22[1-4] starts earlier than an ADCIRC run.

#. If NWS=+12, and if a hot start run is attempted from the 86400th time step,
   and if the first set, fort.221 and fort.222 files, starts at the -86400th
   time step, and if WTIMINC is set to 900 in the fort.15, then you should set
   NWBS to (-86400-0)/900 = -96. If the start times of the two sets are
   different, use the earlier start time, i.e., the start time of the first set.
#. If NWS=-12, and if a hot start run is attempted from the 86400th time step,
   and if the first set, fort.221 and fort.222 files, starts at the -86400th
   time step, and if WTIMINC is set to 900 in the fort.15, then you should set
   NWBS to (-86400-86400)/900 = -192. If the start times of the two sets are
   different, use the earlier start time, i.e., the start time of the first set.

Notes
~~~~~

#. When the basin and region scale grids are both used in ADCIRC, data from the
   region scale grid takes precedence over data from the basin scale grid in the
   areas where the two grids overlap.
#. If the fort.22[1-4] files are shorter in time than the duration of an ADCIRC
   run, the fort.22[1-4] will run out of data. However, the file reading error
   should be safely caught in ADCIRC and the computation will continue with
   blank snaps.

.. _wtiminc_in_fort.15:

WTIMINC in fort.15
------------------

Although the WIN/PRE files contain datetime and timestepping information for the
input data, *WTIMINC* must be set to the timestep of the input files (in
seconds). This limitation locks all nested grid inputs to the same input
timesteps, and requires that the input timesteps are constant for the duration
of the input files, although the start and end times of the nested input overlay
sets may differ.

.. _file_formatting:

File Formatting
---------------

Winds and pressure data formats are similar. The header format is the same, but
in the wind file the header is followed by U then V components while in the
pressure file the header is followed by just pressures.

The file begins with a header indicating the starting and ending dates and is
followed by a grid/date header for each time step and the u and v components of
the wind in meters/second or pressures in millibars. Starting/Ending dates are
in YYYYMMDDHH format where:

|    YYYY    4-character Year
|    MM      2-character Month
|    DD      2-character Day
|    HH      2-character Hour

*example WIN*:

| OWI WWS Wind Output Ucomp,Vcomp in m/s           Start:1995060600 End:1995060600
| iLat=  67iLong=  67DX= 1.250DY=  .833SWLat=  22.500SWlon= -82.500Dt=199506060000
| -1.16856  -1.06439   -.84875  -1.03460  -1.50047  -2.09462  -2.80243  -3.55863
| -4.24125  -4.84273  -5.59486  -5.37088  -5.30224  -5.12534  -4.89537  -4.67412
| -4.49203  -4.35772  -4.26612  -4.20260  -4.14746  -4.08396  -4.00686  -3.92213
| -3.83615  -3.74765  -3.65182  -3.54998  -3.45299  -3.37660  -3.32959  -3.28037
| -3.10631  -2.67723  -2.08363  -1.53773  -1.12623   -.83526   -.62870   -.47371
| [rest deleted]

-  iLat is the number of parallels
-  iLong is the number of meridians
-  DX is the grid spacing in degrees of longitude
-  DY is the grid spacing in degrees of latitude
-  SWLat is the latitude of the South West corner
-  SWlon is the longitude of the South West corner
-  Dt is the date/time in YYYYMMDDHHmm (same as master header date format but
   with mm Minutes as well)

The number of grid points is iLat*iLong, the u component of the winds in
meters/second is followed by the v component.

Sample code
~~~~~~~~~~~

Sample fortran to read a win file (first time step only):

.. code-block:: fortran

   c     Read in beginning/ending dates of win file
   10    format (t56,i10,t71,i10)
         read (20,10) date1,date2

   c     Read Grid Specifications/Date
   11    format (t6,i4,t16,i4,t23,f6.0,t32,f6.0,t44,f8.0,t58,f8.0,t69,i10,i2)
         read (20,11) iLat, iLong, dx, dy, swlat, swlong, lCYMDH, iMin

   c     Read U/V Components of the wind
   12    format (8f10.0)
         read (20,12) ((uu(i,j),i=1,ilong),j=1,ilat)
         read (20,12) ((vv(i,j),i=1,ilong),j=1,ilat)

Latitude and Longitude for each point can be calculated as follows:

.. code-block:: fortran

         do 20 icnt = 1,iLat
            slat(icnt) = SWlat + (icnt - 1) * DY
   20    continue

         do 30 jcnt = 1,iLong
            slon(jcnt) = SWlong + (jcnt - 1) * DX
   30    continue
