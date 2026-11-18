.. _fort22:

Fort.22: Meteorological Forcing Data
====================================

The fort.22 file contains meteorological forcing data (wind velocity and atmospheric pressure) for the ADCIRC model. This file is read when meteorological forcing has been indicated by the NWS parameter in the :doc:`Model Parameter and Periodic Boundary Condition File <fort15>`.

General Notes
-------------

* The spatial extents of the meteorological data must be consistent with the ADCIRC model domain. For example, if ADCIRC uses negative longitude values to indicate locations W of the Greenwich meridian, the meteorological file must be similarly organized.
* Any grid that crosses the Greenwich Meridian should be organized so that the seam occurs @ 180 deg longitude.
* Meteorological data in formats other than OWI (NWS=12) and GFDL (NWS=16) must be provided for the entire model run, otherwise the run will crash!
* The specific data formats and timing requirements vary significantly based on the NWS parameter value.
* Wind stress and pressure units must be carefully adhered to, according to the requirements of each format.

NWS Parameter Index
-------------------

.. list-table:: NWS Parameter Index
   :widths: 35 65
   :width: 100%
   :header-rows: 1
   :class: wrap-table, tight-table

   * - NWS Value
     - Description
   * - :ref:`NWS = 1 or 101 <nws-1-101>`
     - Meteorological data input directly to all nodes in the ADCIRC grid
   * - :ref:`NWS = 2, -2, 102, or -102 <nws-2-102>`
     - Similar to NWS=1/101, but with different timing options for data input
   * - :ref:`NWS = 3 or 103 <nws-3-103>`
     - Meteorological data input on a longitude-latitude grid and interpolated onto the ADCIRC grid
   * - :ref:`NWS = 4, -4, 104, or -104 <nws-4-104>`
     - Meteorological data from the PBL Hurricane Model input directly to a subset of nodes
   * - :ref:`NWS = 5, -5, 105, or -105 <nws-5-105>`
     - Wind velocity and atmospheric pressure at all grid nodes
   * - :ref:`NWS = 6 or 106 <nws-6-106>`
     - Wind velocity and atmospheric pressure are read in for a rectangular grid
   * - :ref:`NWS = 7 or -7 <nws-7>`
     - Surface stress and pressure values on a regular grid
   * - :ref:`NWS = 8 or 108 <nws-8-108>`
     - Hurricane parameters for the Dynamic Holland model
   * - :ref:`NWS = 10 <nws-10>`
     - Wind velocity (10 m) and atmospheric pressure from a sequence of National Weather Service (NWS) Aviation (AVN) model output files
   * - :ref:`NWS = 11 <nws-11>`
     - Wind velocity (10 m) and atmospheric pressure from a sequence of stripped down National Weather Service (NWS) ETA 29km model output files
   * - :ref:`NWS = 12 or -12 <nws-12>`
     - OWI (Oceanweather Inc.) format with nested meteorological grids in ASCII format
   * - :ref:`NWS = 13 or -13 <nws-12>`
     - OWI (Oceanweather Inc.) format with nested meteorological grids in netCDF format
   * - :ref:`NWS = 14 or -14 <nws-14>`
     - Gridded GRIB2 or netCDF wind and pressure data
   * - :ref:`NWS = 15 <nws-15>`
     - H\*Wind files produced by the NOAA Hurricane Research Division (HRD)
   * - :ref:`NWS = 16 <nws-16>`
     - GFDL model output files produced by the Geophysical Fluid Dynamics Laboratory at NOAA
   * - :ref:`NWS = 19 <nws-19>`
     - Asymmetric hurricane vortex formulation with customizable radius and Holland B parameters
   * - :ref:`NWS = 20 <nws-20>`
     - Generalized Asymmetric Holland Model (GAHM) with theoretical and practical improvements over previous models
   * - :ref:`NWS = 30 <nws-30>`
     - Blended GAHM & Background Gridded Wind and Pressure

.. raw:: html

   <style>
   .wrap-table th, .wrap-table td {
     white-space: normal !important;
     word-wrap: break-word !important;
     max-width: 100% !important;
     overflow-wrap: break-word !important;
     hyphens: auto !important;
     vertical-align: top !important;
   }
   </style>

.. _nws-1-101:

NWS = 1 or 101: Wind Stress & Pressure at All Nodes and Timesteps
-----------------------------------------------------------------

In this format, meteorological data is input directly to all nodes in the ADCIRC grid.

File Structure
^^^^^^^^^^^^^^

.. parsed-literal::

   for k=1 to :ref:`NP <NP>`
      :ref:`JN <JN>`, :ref:`WSX(k) <WSX>`, :ref:`WSY(k) <WSY>`, :ref:`PRN(k,j) <PRN>`
   end k loop

   Repeat the block above for each time step until the end of the simulation.

Format Specifics
^^^^^^^^^^^^^^^^

* Meteorological data is input directly to all nodes in the ADCIRC grid.
* The units for pressure are meters of H₂O.
* The first set of met. data corresponds to TIME=STATIM+DTDP. Additional sets of met data must be provided at every time step, (WTIMINC = DTDP).
* Wind stress must be input in units of velocity squared (consistent with the units of gravity) and surface atmospheric pressure must be input in units of equivalent height of water (e.g., meters of water, feet of water that are consistent with the units of gravity).
* Stress in these units is obtained by dividing stress in units of force/area by the reference density of water.
* Pressure in these units is obtained by dividing pressure in units of force/area by the gravitational constant and the reference density of water.
* For example, 10^5Pa =10^5 N/m^2 =10^5 kg m/(s m)^2 divided by 9.81 m/s^2 and 10^3 kg/m^3 equals 10.2 meters of water.

.. _nws-2-102:

NWS = 2, -2, 102, or -102: Wind Stress & Pressure at All Nodes and Specified Time Interval
------------------------------------------------------------------------------------------

This format is similar to NWS=1/101, but with different timing options for data input.

File Structure
^^^^^^^^^^^^^^

.. parsed-literal::

   :ref:`WTIMINC <WTIMINC>`
   for k=1 to :ref:`NP <NP>`
      :ref:`JN <JN>`, :ref:`WSX(k) <WSX>`, :ref:`WSY(k) <WSY>`, :ref:`PRN(k,j) <PRN>`
   end k loop

   Repeat the block above for each time step until the end of the simulation.

Format Specifics
^^^^^^^^^^^^^^^^

* Meteorological data is input directly to all nodes in the ADCIRC grid.
* The units for pressure are meters of H₂O.
* If NWS = 2 or 102, the first set of met. data corresponds to TIME=STATIM.
* If NWS = -2 or -102, the first set of met data corresponds to TIME=HOT START TIME.
* Additional sets of met. data must be provided every WTIMINC, where WTIMINC is the met. data time interval and is specified in the Model Parameter and Periodic Boundary Condition File.
* Met data is interpolated in time to the ADCIRC time step.
* Wind stress must be input in units of velocity squared (consistent with the units of gravity).
* Surface atmospheric pressure must be input in units of equivalent height of water (e.g., meters of water, feet of water that are consistent with the units of gravity).
* Stress in these units is obtained by dividing stress in units of force/area by the reference density of water.
* Pressure in these units is obtained by dividing pressure in units of force/area by the gravitational constant and the reference density of water.
* For example, 10^5Pa =10^5 N/m^2 =10^5 kg m/(s m)^2 divided by 9.81 m/s^2 and 10^3 kg/m^3 equals 10.2 meters of water.

.. _nws-3-103:

NWS = 3 or 103: Fleet Numeric Format
------------------------------------

In this format, meteorological data is input on a longitude-latitude grid and interpolated onto the ADCIRC grid.

File Structure
^^^^^^^^^^^^^^

.. parsed-literal::

   :ref:`IWTIME <IWTIME>`
   for k=1 to :ref:`NWLAT <NWLAT>`
      for j=1 to :ref:`NWLON <NWLON>`
         :ref:`WSPEED(k,j) <WSPEED>`
      end j loop
   end k loop
   for k=1 to :ref:`NWLAT <NWLAT>`
      for j=1 to :ref:`NWLON <NWLON>`
         :ref:`WDIR(k,j) <WDIR>`
      end j loop
   end k loop

   Repeat the block above for each time step until the end of the simulation.

Format Specifics
^^^^^^^^^^^^^^^^

* Meteorological data is input to a longitude, latitude grid and interpolated in space onto the ADCIRC grid. The ADCIRC grid must be in lon, lat coordinates.
* The first set of met. data must be at or before the date and time listed in the Model Parameter and Periodic Boundary Condition File as the beginning time of the simulation.
* Additional sets of met. data must be provided every WTIMINC, where WTIMINC is the met. data time interval.
* Values for NWLAT, NWLON, WTIMINC, and several other parameters must be set in the Model Parameter and Periodic Boundary Condition File.
* Met data is interpolated in time to the ADCIRC time step.
* Wind velocity (@ 10 m above the water surface) must be input in units of m/s (regardless of the units of gravity).
* The following relations are used to compute wind stress from the input wind velocity
   * WIND_SPEED = magnitude of WIND_VEL
   * DRAG_COEFF = 0.001*(0.75+0.067*WIND_SPEED)
   * If (DRAG_COEFF.gt.0.003) DRAG_COEFF=0.003
   * WIND_STRESS = DRAG_COEFF*0.001293*WIND_VEL*WIND_SPEED


.. _nws-4-104:

NWS = 4, -4, 104, or -104: PBL Hurricane Model format
-----------------------------------------------------

In this format, meteorological data from the PBL Hurricane Model is input directly to a subset of nodes in the ADCIRC grid.

File Structure
^^^^^^^^^^^^^^

.. parsed-literal::

   for k=1 to specified nodes
      :ref:`JN <JN>`, :ref:`WVNX(k) <WVNX>`, :ref:`WVNY(k) <WVNY>`, :ref:`PRN(k) <PRN>`
   end k loop
   #
   for k=1 to specified nodes
      :ref:`JN <JN>`, :ref:`WVNX(k) <WVNX>`, :ref:`WVNY(k) <WVNY>`, :ref:`PRN(k) <PRN>`
   end k loop
   #

   Repeat the block above for each time step until the end of the simulation.

Format Specifics
^^^^^^^^^^^^^^^^

* Meteorological data is input directly to a subset of nodes in the ADCIRC grid (as specified by the node number JN). The ADCIRC grid must be either in lon, lat coordinates or in meter-based Cartesian coordinates.
* If NWS = 4 or 104, the first set of met. data corresponds to TIME=STATIM.
* If NWS = -4 or -104, the first set of met data corresponds to TIME=HOT START TIME.
* Additional sets of met. data must be provided every WTIMINC, where WTIMINC is the met. data time interval and is specified in the Model Parameter and Periodic Boundary Condition File.
* Met data is interpolated in time to the ADCIRC time step.
* Each data line must have the format I8, 3E13.5. Data input lines are repeated for as many nodes as desired.
* A line containing the # symbol in column 2 indicates met data at the next time increment begins on the following line.
* At each new time, any node that is not specified in the input file is assumed to have zero wind velocity and pressure = 1013.
* Wind velocity (assumed to be 10m 10 minute averaged value) must be input in knots and surface atmospheric pressure must be input in hundredths of a millibar.
* The following relations are used to compute wind stress from wind velocity
   * WIND_VEL{m/s @ 10m} = WIND_VEL{knots @ bl average}*0.5144 (In prior ADCIRC versions, an additional factor of 1.04 was included in the formulation to convert from 30 minute avg winds to 10 minute avg winds. This factor was removed and it is currently assumed that the input wind data uses a 10 minute averaging period. Note, this is unrelated to the value of WTIMINC).
   * WIND_SPEED = magnitude of WIND_VEL
   * DRAG_COEFF = 0.001*(0.75+0.067*WIND_SPEED)
   * if(DRAG_COEFF.gt.0.003) DRAG_COEFF=0.003
   * WIND_STRESS = DRAG_COEFF*0.001293*WIND_VEL*WIND_SPEED
* The following relationship is used in ADCIRC to convert to pressure in meters of water from pressure in hundredths of a millibar
   * PRESSURE{m H2O}=PRESSURE{Pa/100}*100/(GRAVITY*DENSITY H2O).


.. _nws-5-105:

NWS = 5, -5, 105, or -105: Wind Velocity & Pressure at All Nodes and Specified Time Interval
--------------------------------------------------------------------------------------------

In this format, meteorological data from a tropical cyclone model (Holland formulation) is input directly to all nodes in the ADCIRC grid.

File Structure
^^^^^^^^^^^^^^

.. parsed-literal::

   for k=1 to NP
      :ref:`JN <JN>`, :ref:`WVX(k) <WVX>`, :ref:`WVY(k) <WVY>`, :ref:`PRN(k,j) <PRN>`
   end k loop

   Repeat the block above for each time step until the end of the simulation.

Format Specifics
^^^^^^^^^^^^^^^^

* Meteorological data is input directly to all nodes in the ADCIRC grid. The ADCIRC grid must be either in lon, lat coordinates or in meter-based Cartesian coordinates.
* If NWS = 5 or 105, the first set of met. data corresponds to TIME=STATIM.
* If NWS = -5 or -105, the first set of met data corresponds to TIME=HOT START TIME.
* Additional sets of met. data must be provided every WTIMINC, where WTIMINC is the met. data time interval and is specified in the Model Parameter and Periodic Boundary Condition File.
* Met data is interpolated in time to the ADCIRC time step.
* Wind velocity (@ 10 m above the water surface) must be input in m/s and surface atmospheric pressure must be input in meters of water.
* The following relations are used to compute wind stress from wind velocity:
   * WIND_SPEED = magnitude of WIND_VEL
   * DRAG_COEFF = 0.001*(0.75+0.067*WIND_SPEED)
   * if(DRAG_COEFF.gt.0.003) DRAG_COEFF=0.003
   * WIND_STRESS = DRAG_COEFF*0.001293*WIND_VEL*WIND_SPEED


.. _nws-6-106:

NWS = 6 or 106: Wind Velocity and Pressure on Rectangular Grid
--------------------------------------------------------------

In this format, meteorological data using a generalized asymmetric Holland formulation is input on a rectangular grid.

File Structure
^^^^^^^^^^^^^^

.. parsed-literal::

   IMDAY, IMMONTH, IMYEAR, IMHOUR, IMMIN, IMSEC
   for k=1 to NWLAT
      for j=1 to NWLON
         :ref:`WVXFN(k,j) <WVXFN>`, :ref:`WVYFN(k,j) <WVYFN>`, :ref:`PRN(k,j) <PRN>`
      end j loop
   end k loop

   Repeat the block above for each time step until the end of the simulation.

Format Specifics
^^^^^^^^^^^^^^^^

* In versions 43 and earlier the format of input was P, U, V. The input has now changed to U, V, P to be consistent with other NWS formats.
* Meteorological data is input on a rectangular grid (either in Longitude, Latitude or Cartesian coordinates, consistent with the grid coordinates) and interpolated in space onto the ADCIRC grid.
* In setting up the meteorological grid it is assumed that y (e.g., latitude) varies from north (k=1) to south (k=NWLAT) and x (e.g., longitude) varies from west (j=1) to east (j=NWLON).
* The spatial extents of the meteorological grid must be consistent with the ADCIRC model domain. For example, if ADCIRC uses negative longitude values to indicate locations W of the Greenwich meridian, the meteorological file must be similarly organized.
* Any grid that crosses the Greenwich Meridian should be organized so that the seam occurs @ 180 deg longitude. Therefore, the meteorological and ADCIRC grids should use negative longitudes W of the Greenwich Meridian and positive longitudes to the E.
* The meteorological grid MUST cover the entire ADCIRC mesh; that is, the ADCIRC mesh must be ENTIRELY within the meteorological grid or an error will result.
* The first set of met. data corresponds to the beginning time of the current simulation.
   * If the model is cold started this corresponds to TIME=STATIM.
   * If the model is hot started, this corresponds to TIME=HOT START TIME.
* Additional sets of met. data must be provided every WTIMINC, where WTIMINC is the met. data time interval.
* Values for NWLAT, NWLON, WTIMINC, and several other parameters must be set in the Model Parameter and Periodic Boundary Condition File.
* Met data is interpolated in time to the ADCIRC time step.
* Wind velocity (@ 10 m above the water surface) must be input in units of m/s and surface atmospheric pressure must be input in units of Pascals = Newtons/square meter.
* The following relations are used to compute wind stress from the input wind velocity:
   * WIND_SPEED = magnitude of WIND_VEL
   * DRAG_COEFF = 0.001*(0.75+0.067*WIND_SPEED)
   * If (DRAG_COEFF.gt.0.003) DRAG_COEFF=0.003
   * WIND_STRESS = DRAG_COEFF*0.001293*WIND_VEL*WIND_SPEED
* The following relationship is used in ADCIRC to convert to pressure in meters of water from pressure in Pascal:
   * PRESSURE{m H2O}=PRESSURE{Pascal}/(GRAVITY*DENSITY H2O).


.. _nws-7:

NWS = 7 or -7: Wind Stress and Pressure on Rectangular Grid
-----------------------------------------------------------

This format has been available since ADCIRC version 55.

File Structure
^^^^^^^^^^^^^^

.. parsed-literal::

   for k=1 to NWLAT
      for j=1 to NWLON
         :ref:`WVXFN(k,j) <WVXFN>`, :ref:`WVYFN(k,j) <WVYFN>`, :ref:`PRN(k,j) <PRN>`
      end j loop
   end k loop

   Repeat the block above for each time step until the end of the simulation.

.. _nws-8-108:

NWS = 8 or 108: Dynamic Symmetric Holland Vortex Model
------------------------------------------------------

In this format, a single tropical cyclone model using the Generalized asymmetric Holland formulation with a moving reference frame and a background wind field is used.

File Structure
^^^^^^^^^^^^^^

The file uses the `ATCF Best Track/Objective Aid/Wind Radii Format <https://science.nrlmry.navy.mil/atcf/docs/database/new/abrdeck.html>`_. Below is a simplified representation of the format structure:

.. code-block:: none

   AL, ##, YYYY, MM/DD/HH24, TECH, TAU, LatN/S, LongE/W, VMAX, MSLP, TY, ...

Format Specifics
^^^^^^^^^^^^^^^^

* Hurricane parameters are read in from the Single File Meteorological Forcing Input File.
* Wind velocity and atmospheric pressure are calculated at every node on the fly by ADCIRC internally using the Dynamic Holland model.
* The input file is fixed width (not comma separated values or csv) and is assumed to correspond to the ATCF Best Track/Objective Aid/Wind Radii Format.
* Historical tracks, real-time hindcast tracks and real-time forecast tracks may be found in this format.
* Selecting NWS = 8 also requires the specification of the cold start time, storm number, and boundary layer adjustment (see YYYY MM DD HH24 StormNumber BLAdj).
* Garret's formula is used to compute wind stress from the wind velocity.
* The symmetric vortex model (NWS=8) in ADCIRC assumes that the longitudes in the fort.22 are west longitude, so it multiplies the longitude values by -1. It ignores the 'E' or 'W' in the longitude column of the fort.22.
* The symmetric vortex model (NWS=8) does not use any of the isotach wind speeds or wind radii data.
* When reading lines labeled "BEST" from the fort.22, it obtains timing information from the year, month, day, and hour in column 3.
* When reading lines labeled "OFCL" from the fort.22, it uses the forecast increment (a.k.a. TAU) from column 6.
* The use of these two different columns by ADCIRC NWS=8 is to maintain consistency with the official file structure for the ATCF file format.
* For NWS=8, ADCIRC knows the current time because the user provides the year, month, day, and hour of cold start on the WTIMINC line in the fort.15 file.
* It also has the time that has elapsed since cold start, because that is provided in the hotstart file, if any.
* It then compares the current time with the date/times in the fort.22 to automatically find the right place to begin reading data from the fort.22.
* If the whole fort.22 consists of "BEST" lines, the symmetric vortex model (NWS=8) only looks at column 3 for time information.
* It automatically knows where to start reading cyclone data, based on the coldstart date/time the user provides in the fort.15 file.


NWS = 9: Asymmetric Holland Vortex Model
----------------------------------------

**This has been deprecated and is no longer available.**

.. _nws-10:

NWS = 10: National Climatic Data Center GFS
-------------------------------------------

In this format, meteorological data is processed in a specified manner.

File Structure
^^^^^^^^^^^^^^

.. code-block:: none

   for k=1, LONB*LATB
   
   PG(k), UG(k), VG(k)
   
   end j loop

.. _nws-11:

NWS = 11: Stripped National Weather Service (NWS) ETA 29km
----------------------------------------------------------

In this format, meteorological data is processed in a specified manner.

File Structure
^^^^^^^^^^^^^^

.. code-block:: none

   for k=1, LONB*LATB
   
   PG(k), UG(k), VG(k)
   
   end j loop

.. _nws-12:

NWS = 12 or -12: Oceanweather ASCII Format Gridded Wind and Pressure
--------------------------------------------------------------------

Oceanweather Inc (OWI) ASCII "WIN"/"PRE" format, details can be found at :ref:`NWS12 <nws12>`.

.. _nws-13:

NWS = 13 or -13: Oceanweather NetCDF Format Gridded Wind and Pressure
---------------------------------------------------------------------

Oceanweather Inc (OWI) NetCDF format, details can be found at :ref:`NWS13 <nws13>`.

.. _nws-14:

NWS = ±14: Gridded GRIB2 or NetCDF Wind and Pressure
----------------------------------------------------

When using NetCDF files (fort.221.nc, fort.222.nc) as met inputs, the fort.22 is required in order to list the relevant variable names so that the internal NetCDF read routines can find the pertinent variables. In the case of GRIB2 input files (fort.221.grb2, fort.222.grb2) the fort.22 is not required because the variable names are standardized and on start-up the internal wgrb2api library prints out inventory look-up files (\*.inv) that it uses to find the information contained within the \*.grb2 files.

This format is available in ADCIRC version 55 and later.

File Structure
^^^^^^^^^^^^^^

.. parsed-literal::

   Temporal dimension name
   Datetime variable name
   Format of the datetime time variable [special note, if first character is not a '%' then it will ignore this variable and will assume to start from the first time snap]
   Zonal (east-west) dimension name
   Longitude variable name
   Meridional (north-south) dimension name
   Latitude variable name
   Sea-level Pressure variable name [plus optional HPa units]
   Zonal 10-m Wind Velocity name
   Meridional 10-m Wind Velocity name
   Ice area-fraction name [optional, only if fort.225.nc present for ice area fraction]

Example (WRF output)
^^^^^^^^^^^^^^^^^^^^

.. parsed-literal::

   Time
   Times [if no datetime variable then set to dummy such as 'none']
   %Y-%m-%d_%H:%M:%S [if no datetime variable then set to dummy such as 'none']
   west_east
   XLONG
   south_north
   XLAT
   PSFC [optional: append HPa to PSFC (i.e., PSFCHPa) to indicate that units of pressure are in hPa/mbar, otherwise assumed to be in Pa]
   U10
   V10
   aice [optional]

Notes
^^^^^

* If the time variable is not in a datetime format (e.g., is a float in minutes since..) the datetime format can be set to a dummy name, e.g., 'minutes'. The code checks to see if the first character of the datetime format is '%'. If not the code will assume to simply begin reading from the first time snap. If the datetime variable is available the code will work out which time snap to start reading from based on the reference date, NCDATE located near or at the bottom of the fort.15 file. If pressure units are in hectopascals, add 'HPa' to the pressure variable line in the fort.22. For example, if the pressure variable is 'mslp', make it 'mslpHPa'.

.. _nws-15:

NWS = 15: H*Wind Gridded Wind and Inferred Pressure
---------------------------------------------------

In this format, ADCIRC uses HWind files produced by the NOAA Hurricane Research Division (HRD).

File Structure
^^^^^^^^^^^^^^

The goal of the implementation of the HWind capability within ADCIRC was to allow HWind files to be used as-is, without resorting to an intermediate format. As a result, the fort.22 file consists of a header line to provide some configuration parameters, and then a list of the filenames of the HWind files to be used in the ADCIRC run. Specifically, the file format is as follows:

.. parsed-literal::

   :ref:`comment line <comment_line>`
   hWindMultiplier
   :ref:`pressureWindRelationship <pressureWindRelationship>`
   for i=1 to numHWindFiles
      :ref:`hours(i) <hours>` :ref:`centralPressure(i) <centralPressure>` :ref:`rampMult(i) <rampMult>` :ref:`filename(i) <filename>`
   end i loop

Format Specifics
^^^^^^^^^^^^^^^^

* HWind files are data assimilated snapshots of the wind velocity fields of tropical cyclones that are produced by the NOAA Hurricane Research Division (HRD). The files have the following characteristics:
   * the format explicitly indicates the center of the storm
   * the (u,v) data are on a regular grid
   * the regular grid is a mercator projection with origin at storm center
   * the mercator grid spacing is in meters and is uniform in x and y (dx=dy)
   * the dimensions (nx,ny) of the mercator grid are equal (nx=ny)
   * the grid dimensions change from snapshot to snapshot; for example, the first shapshot may be 161×161 while the 2nd snapshot may be 121×121 (etc)
   * sequential hwind snapshots will not be evenly spaced in time for a particular storm
   * hwind data do not contain barometric pressure information
   * For the dvorak, knaffzehr, and specifiedPc options, the barometric pressure field is computed by determining the radius to maximum winds Rmax (i.e., the distance of Vmax from the center of the storm), calculating the Holland B parameter, and then using the Holland formulation to calculate barometric pressure as a function of the distance from the center of the storm.

Example
^^^^^^^

The following is an example of a fort.22 file with NWS=15 format:

.. code-block:: none

   ! first line is a comment line, max length 1024 characters
   1.0 ! 2nd line is a velocity magnitude multiplier
   dvorak ! 3rd line: one word for the pressure-wind relationship
   0.0 -1 0.0 "/home/jason/hwind/al092011_0828_1330" ! time (hours), Pc (mb), ramp mult, filename
   6.0 -1 0.5 "/home/jason/hwind/al092011_0828_1930"
   12.0 -1 1.0 "/home/jason/hwind/al092011_0829_0130"

.. _nws-16:

NWS = 16: ASCII NOAA GFDL Gridded Wind and Pressure
---------------------------------------------------

In this format, ADCIRC uses GFDL model output files produced by the Geophysical Fluid Dynamics Laboratory at NOAA.

File Structure
^^^^^^^^^^^^^^

The GFDL input capability uses GFDL model output files as-is; as a result, the fort.22 file consists of a list of GFDL model output files to be used in ADCIRC. The file format is as follows:

.. parsed-literal::

   comment line
   GFDLWindMultiplier
   MaxExtrapolationDistance
   for i=1 to numGFDLFiles
      :ref:`hours(i) <hours>` :ref:`rampMult(i) <rampMult>` :ref:`filename(i) <filename>`
   end i loop

Format Specifics
^^^^^^^^^^^^^^^^

* GFDL model output files are produced by the Geophysical Fluid Dynamics Laboratory at NOAA.
* Each ASCII GFDL model output file contains one or more nested grid dataset where the nested grids are allowed to change in time.
* Coarse grid data is not stored where finer nest data is given.
* The files are formatted as follows:
   * Line 1: Number of grid cells (f10.4) NCELLS
   * Lines 2 through NCELLS+1: Ten columns of data formatted as 10f10.4 as follows:
      * column 1: u (m/s)
      * column 2: v (m/s)
      * column 3: Temperature (K)
      * column 4: mixing ratio(kg/kg)
      * column 5: storm accum precipitation (cm)
      * column 6: sea level pressure (hPa)
      * column 7: longitude (decimal deg)
      * column 8: latitude (decimal deg)
      * column 9: hurricane hour
      * column 10: nest number (this is not always present)
* If the ADCIRC time falls outside the interval of time covered by the GFDL model output files, ADCIRC will insert "blank snaps", i.e., it will set the wind velocity at all mesh vertices to 0.0 m/s and the barometric pressure to a uniform background pressure of 1013mb.

Example
^^^^^^^

The following is an example of a fort.22 file with NWS=16 format:

.. code-block:: none

   ! first line is a comment line, max length 1024 characters
   1.0 ! 2nd line is a velocity magnitude multiplier
   100.0 ! 3rd line: maximum extrapolation distance (m)
   0.0 -1 0.0 "/home/jason/hwind/al092011_0828_1330" ! time (hours), Pc (mb), ramp mult, filename
   6.0 -1 0.5 "/home/jason/hwind/al092011_0828_1930"


.. _nws-19:

NWS = 19: Dynamic Asymmetric Holland Vortex Model
-------------------------------------------------

** Use of this is discouraged. **

In this format, ADCIRC uses an asymmetric hurricane vortex formulation with customizable radius and Holland B parameters.

File Structure
^^^^^^^^^^^^^^

The file needs to be in best track format, see notes below for more information.

Format Specifics
^^^^^^^^^^^^^^^^

* User has the ability to select which Isotach to use in each of the 4 quadrants.
* User also has ability to modify RMAX and Holland's B parameter using the ASWIP program.
* The auxiliary preprocessing program ASWIP.F (located in the /wind directory and executable is created by typing, make aswip, in the work folder after adcirc executable has been generated), will generate the fort.22 input file for NWS=19 from a NWS=9 formatted input file.
* Hurricane parameters are read in from the Single File Meteorological Forcing Input File.
* It is assumed that the line in the fort.22 file with a zero as the forecast increment (i.e., column 6) corresponds to the start of the current simulation run, whether it is a hotstart or cold start.
* In other words, there is no option to set the NWS value negative to indicate that the file starts at the ADCIRC hotstart time. Rather, the forecast increment in hours (column 6) is used to indicate the relationship between the ADCIRC time and the data in the fort.22 file.
* Wind velocity and atmospheric pressure are calculated at exact finite element mesh node locations and directly coupled to ADCIRC at every time step using the asymmetric hurricane vortex formulation (Mattocks et al, 2006; Mattocks and Forbes, 2008) based on the Holland gradient wind model.
* The input file is assumed to correspond to the ATCF Best Track/Objective Aid/Wind Radii Format.
* Historical tracks, real-time hindcast tracks and real-time forecast tracks may be found in this format.
* This option uses the radii at specific wind speeds (34, 50, 64, 100 knots) reported in the four quadrants (NE, SE, SW, NW) of the storm to calculate the radius of maximum winds as a function of the azimuthal angle.
* Garret's formula is used to compute wind stress from the wind velocity.
* The NWS=19 option allows the user to set a value for Rmax and Holland B Parameter.
* Additionally the user can select the isotachs to be used for each of the 4 quadrants.
* The utility program aswip_1.0.3.F located in the /wind folder will generate the NWS=19 formatted file from a NWS=9 formatted fort.22 input file.
* In order to use the NWS=19 option, the file needs to be in best track format.
* The forecast period (column #6) needs to be edited to reflect the time of the forecast/nowcast for each track location (each line) in hours from the start of the simulation (0, 6, 12, 18, etc).
* The original data in that column depends on what type of best track format data is being used. The original data might have 0 or other numbers in that column.
* It is suggested that users change the "BEST" tech type to "ASYM" in column 5 in the fort.22 file to denote that the file has been modified to accommodate the asymmetric wind formulation (the simulation time in hours in the 6th column has been added, etc.) so it will not get confused in the future with a best track file.
* The NWS=19 option requires the following variables in the fort.22 file in a best track format
   * Column 6: Forecast time in hours (enter the time in hours in each record starting at 0)
   * Column 7: Latitude of the eye
   * Column 8: Longitude of the eye
   * Column 9: Maximum sustained wind speed in knots
   * Column 10: Minimum sea level pressure in MB
   * Column 12: Wind intensity in knots of the radii defined in the record (34, 50, 64 or 100 knots)
   * Columns 14, 15, 16, 17: Radius of specified wind intensity for quadrants 1, 2, 3, 4 in NM; ≥ 0
   * Column 18: Background pressure in MB; a standard value of 1013 can be used
   * Column 20: Rmax as reported in the ATCF BEST TRACK file
   * Column 28: Storm Name in ATCF file format
   * Column 29: Time Record number. There can be multiple lines for a given time record depending on the number of isotachs reported in the ATCF File
   * Column 30: Number of isotachs reported in the ATCF file for the corresponding Time record.
   * Columns 31-34: Selection of radii for that particular isotach. 0 indicates do not use this radius, and 1 indicates use this radius and corresponding wind speed.
   * Columns 35-38: Designated Rmax values computed for each of the quadrants selected for each particular isotach.
   * Column 39: Holland B parameter computed using the formulas outlined in the Holland paper, and implemented using the aswip program.
* The format of the file is fixed and users will want to use the aswip program to be sure that the input fort.22 file is properly formatted.
* The command line for NWS=19 is ``./aswip -n 19 -m 2 -z 1``


.. _nws-20:

NWS = 20: Generalized Asymmetric Holland Vortex Model (GAHM)
------------------------------------------------------------

In this format, ADCIRC uses the Generalized Asymmetric Holland Model (GAHM), which provides theoretical and practical improvements over previous parametric vortex models. See :ref:`Generalized Asymmetric Holland Vortex Model <gahm>` for more information.

File Structure
^^^^^^^^^^^^^^

The file format is similar to the NWS = 19 format with 8 additional columns of data, see notes below for more information.

Format Specifics
^^^^^^^^^^^^^^^^

* The Generalized Asymmetric Holland Model (GAHM) provides a set of theoretical and practical improvements over previous parametric meteorological vortex models in ADCIRC. The theory and implementation of the GAHM was initially described at the 2013 ADCIRC Users Group Meeting.
* The NWS=20 option requires the following variables in the fort.22 file in a best track format, which includes the same variables as NWS=19 plus additional columns:
   * Column 6: Forecast time in hours (enter the time in hours in each record starting at 0)
   * Column 7: Latitude of the eye
   * Column 8: Longitude of the eye
   * Column 9: Maximum sustained wind speed in knots
   * Column 10: Minimum sea level pressure in MB
   * Column 12: Wind intensity in knots of the radii defined in the record (34, 50, 64 or 100 knots)
   * Columns 14, 15, 16, 17: Radius of specified wind intensity for quadrants 1, 2, 3, 4 in NM; ≥ 0
   * Column 18: Background pressure in MB; a standard value of 1013 can be used
   * Column 20: Rmax as reported in the ATCF BEST TRACK file
   * Column 28: Storm Name in ATCF file format
   * Column 29: Time Record number. There can be multiple lines for a given time record depending on the number of isotachs reported in the ATCF File
   * Column 30: Number of isotachs reported in the ATCF file for the corresponding Time record.
   * Columns 31-34: Selection of radii for that particular isotach. 0 indicates do not use this radius, and 1 indicates use this radius and corresponding wind speed.
   * Columns 35-38: Designated Rmax values computed for each of the quadrants selected for each particular isotach.
   * Column 39: Holland B parameter computed using the formulas outlined in the Holland paper, and implemented using the aswip program.
   * Columns 40-43: Quadrant-varying Holland B parameter
   * Columns 44-47: Quadrant-varying Vmax calculated at the top of the planetary boundary (a wind reduction factor is applied to reduce the wind speed at the boundary to the 10-m surface)
* The format of the file is fixed and users will want to use the aswip program to be sure that the input fort.22 file is properly formatted.
* Options for the aswip program using NWS = 20 are the following:
   * ``./aswip -n # -m # -z #``
   * -n = nws option
   * -m = methods of selecting isotachs for use in computation of radius/radii to maximum winds
     * 1: always use the 34kt isotach
     * 2: use the highest available isotach in any quadrant each time
     * 3: use the 50kt isotach if it is available; otherwise use the 34kt isotach
     * 4: use all available isotachs (must choose this for GAHM/NWS=20)
   * -z = approaches solving for Rmax, 1 = only rotate wind vectors afterward, 2 = rotate wind vectors before and afterwards (use this for NWS=20)
* The command line for NWS=20 is: ``./aswip -n 20 -m 4 -z 2``

.. _nws-30:

NWS = 30: Blended GAHM & Background Gridded Wind and Pressure
-------------------------------------------------------------

This format is a combination of GAHM (NWS=20) and a gridded background meteorological field (NWS=12). The core of a tropical cyclone is represented using the GAHM model and read in from a file named NWS_20_fort.22, which should be created following the procedure presented in the section on NWS=20. The user must also supply gridded wind and pressure files in the form of an OWI-style fort.22 file, a fort.221 file, and a fort.222 file, details on these are in the section on NWS=12. ADCIRC internally blends the wind and pressure fields from these two sets of inputs over a user-controlled distance. Specifically, the pureVortex and pureBackground inputs in the fort.15 file (for formatting, see info on the meteorological parameter line) are coefficients used to define where transitions between the vortex and background meteorology occur. At or within the distance pureVortex*vortexRMW of the storm's center, the GAHM meteorology is used. At or beyond the distance pureBackground*vortexRMW of the storm's center, the background meteorology is used. In between these, a linear distance-weighted average of the two sets of meteorology is used to define the forcing.
