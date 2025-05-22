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
     - Hurricane parameters for the Dynamic Holland model.
   * - :ref:`NWS = 10 <nws-10>`
     - Wind velocity (10 m) and atmospheric pressure from a sequence of National Weather Service (NWS) Aviation (AVN) model output files.
   * - :ref:`NWS = 11 <nws-11>`
     - Wind velocity (10 m) and atmospheric pressure from a sequence of stripped down National Weather Service (NWS) ETA 29km model output files.
   * - :ref:`NWS = 12 or -12 <nws-12>`
     - Meteorological data uses the OWI (Oceanweather Inc.) format with nested meteorological grids
   * - :ref:`NWS = 15 <nws-15>`
     - HWind files produced by the NOAA Hurricane Research Division (HRD)
   * - :ref:`NWS = 16 <nws-16>`
     - GFDL model output files produced by the Geophysical Fluid Dynamics Laboratory at NOAA
   * - :ref:`NWS = 19 <nws-19>`
     - Asymmetric hurricane vortex formulation with customizable radius and Holland B parameters
   * - :ref:`NWS = 20 <nws-20>`
     - Generalized Asymmetric Holland Model (GAHM) with theoretical and practical improvements over previous models

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

NWS = 1 or 101
--------------

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

Example
^^^^^^^

The following is a simple example of a fort.22 file with NWS=1 format (3 nodes, 2 time steps):

.. code-block:: none

   1 10.0 5.0 10.2
   2 12.0 6.0 10.1
   3 14.0 7.0 10.0
   1 15.0 7.5 10.1
   2 16.0 8.0 10.0
   3 17.0 8.5 9.9

This example shows:

* Three nodes with wind velocity and pressure data (nodes 1, 2, and 3)
* Two time steps of data
* Wind velocity components (WSX, WSY) in m/s and pressure (PRN) in meters of H₂O
* For NWS=1, these wind velocities would be converted to stress using the quadratic drag law

.. _nws-2-102:

NWS = 2, -2, 102, or -102
-------------------------

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

Example
^^^^^^^

The following is a simple example of a fort.22 file with NWS=2 format (3 nodes, 2 time steps):

.. code-block:: none

   3600.0
   1 10.0 5.0 10.2
   2 12.0 6.0 10.1
   3 14.0 7.0 10.0
   1 15.0 7.5 10.1
   2 16.0 8.0 10.0
   3 17.0 8.5 9.9

This example shows:

* A time increment (WTIMINC) of 3600.0 seconds (1 hour)
* Three nodes with wind velocity and pressure data
* Two time steps of data
* Wind velocity components (WSX, WSY) in m/s and pressure (PRN) in meters of H₂O
* For NWS=2, these wind velocities would be converted to stress
* For NWS=102, these would be interpreted as wind stress values directly

.. _nws-3-103:

NWS = 3 or 103
--------------

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

Example
^^^^^^^

The following is a simplified example of a fort.22 file with NWS=3 format (2x3 grid, 1 time step):

.. code-block:: none

   2015083106
   8.2 9.1 8.7
   7.8 8.5 8.2
   270.0 265.0 275.0
   268.0 272.0 270.0

This example shows:

* Timestamp for the data (IWTIME): August 31, 2015, 06:00 UTC
* Wind speed (m/s) for a 2x3 (NWLAT x NWLON) grid
* Wind direction (degrees) for the same grid
* Values would be interpolated to the ADCIRC mesh points

.. _nws-4-104:

NWS = 4, -4, 104, or -104
-------------------------

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

Example
^^^^^^^

The following is a simplified example of a fort.22 file with NWS=4 format (sparse node selection, 2 time steps):

.. code-block:: none

       123     10.50     15.75    1010.25
       456     12.25     14.50    1009.75
       789      9.75     16.25    1011.50
    #
       123     11.00     16.50    1009.50
       456     13.00     15.25    1008.75
       789     10.25     17.00    1010.75
    #

This example shows:

* Three selected nodes (123, 456, and 789) with wind and pressure data
* Two time steps separated by the # symbol
* Wind velocity components (WVNX, WVNY) in knots
* Pressure (PRNP) in hundredths of a millibar
* The specified format I8, 3E13.5 for each data line

.. _nws-5-105:

NWS = 5, -5, 105, or -105
-------------------------

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

Example
^^^^^^^

The following is a simplified example of a fort.22 file with NWS=5 format (3 nodes, 2 time steps):

.. code-block:: none

   1 10.0 5.0 10.2
   2 12.0 6.0 10.1
   3 14.0 7.0 10.0
   1 15.0 7.5 10.1
   2 16.0 8.0 10.0
   3 17.0 8.5 9.9

This example shows:

* Three nodes with wind velocity and pressure data (nodes 1, 2, and 3)
* Two time steps of data
* Wind velocity components (WVX, WVY) in m/s and pressure (PRN) in meters of H₂O
* For NWS=5/-5, these wind velocities would be converted to stress using the formula provided in Format Specifics

.. _nws-6-106:

NWS = 6 or 106
--------------

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

Example
^^^^^^^

The following is a simplified example of a fort.22 file with NWS=6 format (2x3 grid, 1 time step):

.. code-block:: none

   28 09 2022 12 00 00
   10.2 -5.1 101300
   10.8 -5.4 101250
   11.2 -5.8 101200
   9.8 -4.9 101320
   10.5 -5.2 101270
   11.0 -5.6 101220

This example shows:

* Timestamp for the data: September 28, 2022, 12:00:00
* Wind velocity components (WVXFN, WVYFN) in m/s and pressure (PRN) in Pascals for a 2x3 (NWLAT x NWLON) grid
* Each line contains the U-component, V-component, and pressure for a single grid point
* Values would be interpolated to the ADCIRC mesh points

.. _nws-7:

NWS = 7 or -7
-------------

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

NWS = 8 or 108
--------------

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

Example
^^^^^^^

The following is a simplified example of a fort.22 file with NWS=8 format:

.. code-block:: none

   3600.0 2017 09 05 00 09 1.0
   AL, 09, 2017090500, 03, BEST, 000, 167N, 0369W, 125, 948, XX,  34, NEQ,  150,  150,  090,  150,
   AL, 09, 2017090506, 03, BEST, 000, 171N, 0378W, 130, 943, XX,  34, NEQ,  150,  150,  090,  150,
   AL, 09, 2017090512, 03, BEST, 000, 175N, 0386W, 140, 938, XX,  34, NEQ,  160,  160,  090,  160,
   AL, 09, 2017090518, 03, BEST, 000, 180N, 0395W, 155, 925, XX,  34, NEQ,  170,  170,  100,  170,
   AL, 09, 2017090600, 03, BEST, 000, 186N, 0404W, 155, 920, XX,  34, NEQ,  180,  180,  110,  180,

This example shows:

* First line with WTIMINC (time increment) set to 3600.0 seconds (1 hour)
* First line also includes coldstart date (2017-09-05 00:00), storm number (09), and boundary layer adjustment (1.0)
* ATCF format data for Hurricane Irma (AL09 of 2017) at different time steps
* Each line includes:
   * Basin code (AL for Atlantic)
   * Storm number (09)
   * Date/time (YYYYMMDDHH)
   * Forecast source (BEST for best track)
   * Forecast hour (000 for analysis)
   * Position (latN/S, longW/E)
   * Maximum sustained wind (VMAX in knots)
   * Minimum sea level pressure (MSLP in millibars)
   * Storm type (XX)
   * Wind radii information (not used by ADCIRC NWS=8)

.. _nws-10:

NWS = 10
--------

In this format, meteorological data is processed in a specified manner.

File Structure
^^^^^^^^^^^^^^

.. code-block:: none

   for k=1, LONB*LATB
   
   PG(k), UG(k), VG(k)
   
   end j loop

.. _nws-11:

NWS = 11
--------

In this format, meteorological data is processed in a specified manner.

File Structure
^^^^^^^^^^^^^^

.. code-block:: none

   for k=1, LONB*LATB
   
   PG(k), UG(k), VG(k)
   
   end j loop

.. _nws-12:

NWS = 12 or -12
---------------

In this format, meteorological data uses the OWI (Oceanweather Inc.) format with the ability to use nested meteorological grids.

File Structure
^^^^^^^^^^^^^^

Unlike other NWS values, the OWI format uses multiple files:

.. parsed-literal::

   fort.22 (control file):
   :ref:`NWSET <NWSET>`
   :ref:`NWBS <NWBS>`
   :ref:`DWM <DWM>`
   
   fort.221 (basin scale pressure):
   Header line with timing information
   Grid information line
   Pressure data for entire grid
   
   fort.222 (basin scale wind):
   Header line with timing information
   Grid information line
   U-component data for entire grid
   V-component data for entire grid
   
   fort.223 (optional regional scale pressure):
   Same format as fort.221
   
   fort.224 (optional regional scale wind):
   Same format as fort.222

Format Specifics
^^^^^^^^^^^^^^^^

* Most types of meteorological input data for ADCIRC use the fort.22 exclusively; however, the OWI format (activated by setting NWS=12 or -12 in the fort.15) only uses that file for a few control parameters.
* This is because the OWI format is capable of using nested meteorological grids (a large, coarse, basin-scale grid and a smaller, finer region-scale grid), and further specifies that wind and pressure information be stored in separate files, along with information about the size, shape and location of the structured grid(s).
* The basin scale pressure field must be placed in a file called "fort.221", and the basin scale wind field must be placed in a file called "fort.222".
* If regional scale meteorological fields are also used, they must be placed in files called "fort.223" and "fort.224".
* The wind velocity data must be in m/s (10 minute averaged, at 10m) and pressure data must be in millibars.
* The data files are fixed width, meaning that ADCIRC interprets the data according to the exact text column in which it appears.
* The first line in each of the files (fort.221, fort.222, fort.223, fort.224) must be a header line that indicates the type of meteorological data in the file and the starting and ending dates of the meteorological data.
* Each data set starts with a header line with the following information:
   * iLat: number of parallels
   * iLong: number of meridians
   * DX: grid spacing in degrees longitude
   * DY: grid spacing in degrees latitude
   * SWLat: latitude of the South West corner
   * SWLon: longitude of the South West corner
   * Dt: starting time of the wind data, in YYYYMMDDHH24mm format
* Although the data files contain timing information, ADCIRC does not use it. Instead, the user must set WTIMINC in the fort.15 file to the time step of the two gridded meteorological data fields.
* As a result, ADCIRC requires that this time step be of a constant size and that all the OWI meteorological data fields be synchronized to it.
* When the basin and region scale grids are both used in ADCIRC, data from the region scale grid takes precedence over data from the basin scale grid in the areas where the two grids overlap.
* ADCIRC also has the ability to insert "blank" meteorological data if it is started before the beginning of the wind data (for a tidal spinup run, for example).
* Blank data is characterized by zero wind speed and 1013 mb of atmospheric pressure.
* This capability is activated via the NWBS parameter in the fort.22 file and it interacts with the sign of the NWS value in order to provide full control over the relationship between the ADCIRC start time (either hot start or cold start) and the beginning of the meteorological data.
* If NWBS is set to a positive number, NWBS specifies the number of blank snaps to be inserted before any information is read from fort.22[1-4].
* If NWBS is set to a negative number, then NWBS determines how many snaps in the fort.22[1-4] should be skipped before the values in the files are used.
* If the fort.22[1-4] files are shorter in time than the duration of an ADCIRC run, the fort.22[1-4] will run out of data. However, the file reading error should be safely caught in ADCIRC and the computation will continue with blank snaps.

Example
^^^^^^^

The following is a simplified example of the fort.22 file with NWS=12 format:

.. code-block:: none

   192

Example of fort.221 (basin scale pressure):

.. code-block:: none

   OWI WWS Pressure Output in mb        Start:1995060600 End:1995060621
   iLat= 67iLong= 67DX= 1.250DY= .833SWLat= 22.500SWlon= -82.500Dt=199506060000
   1013.0 1013.0 1013.0 1013.0 1013.0 1013.0 1012.9 1012.9 ...
   ... (pressure data for entire grid) ...

Example of fort.222 (basin scale wind):

.. code-block:: none

   OWI WWS Wind Output Ucomp,Vcomp in m/s Start:1995060600 End:1995060621
   iLat= 67iLong= 67DX= 1.250DY= .833SWLat= 22.500SWlon= -82.500Dt=199506060000
   0.0 0.0 0.0 0.0 0.1 0.1 0.2 0.2 ...
   ... (U-component for entire grid) ...
   0.0 0.0 0.0 0.0 0.1 0.1 0.2 0.2 ...
   ... (V-component for entire grid) ...

This example shows:

* fort.22 with NWBS = 192 (indicating 192 blank snaps before reading meteorological data)
* Header lines in fort.221 and fort.222 indicating data type, start time, and end time
* Grid specification lines indicating a 67×67 grid with specific spacing and location
* Actual data values for pressure (in millibars) and wind components (in m/s)
* Regional scale files (fort.223 and fort.224) would follow the same format if used

.. _nws-15:

NWS = 15
--------

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

This example shows:

* A comment line (first line)
* A velocity magnitude multiplier of 1.0 (second line)
* The "dvorak" pressure-wind relationship (third line)
* Three HWind files at different times:
   * The first at 0.0 hours with central pressure automatically calculated (-1), no ramping (0.0)
   * The second at 6.0 hours with central pressure automatically calculated (-1), half ramping (0.5)
   * The third at 12.0 hours with central pressure automatically calculated (-1), full ramping (1.0)
* Each file represents a snapshot of Hurricane Irene (AL09 of 2011) at different times
* The -1 for centralPressure indicates that ADCIRC should calculate the central pressure based on the specified pressure-wind relationship

.. _nws-16:

NWS = 16
--------

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

This example shows:

* A comment line (first line)
* A velocity magnitude multiplier of 1.0 (second line)
* A maximum extrapolation distance of 100.0 meters (third line)
* Two GFDL model output files at different times:
   * The first at 0.0 hours with no ramping (0.0)
   * The second at 6.0 hours with half ramping (0.5)
* Each file represents a snapshot of GFDL model output at the specified time
* The timestamp value corresponds to hours after the cold start time

.. _nws-19:

NWS = 19
--------

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
* The NWS=19 option requires the following variables in the fort.22 file in a best track format:
   * Forecast time in hours (column 6); enter the time in hours in each record starting at 0
   * Latitude of the eye (column 7)
   * Longitude of the eye (column 8)
   * Maximum sustained wind speed in knots (column 9)
   * Minimum sea level pressure in MB (column 10)
   * Wind intensity in knots of the radii defined in the record (34, 50, 64 or 100 knots) (column 12)
   * Radius of specified wind intensity for quadrants 1, 2, 3, 4 in NM (columns 14, 15, 16, 17); ≥ 0
   * Background pressure in MB (column 18); a standard value of 1013 can be used
   * Rmax as reported in the ATCF BEST TRACK file in column 20
   * Storm Name in Column 28 ATCF file format
   * Time Record number in column 29. There can be multiple lines for a given time record depending on the number of isotachs reported in the ATCF File
   * Number of isotachs reported in the ATCF file for the corresponding Time record.
   * Columns 31-34 indicate the selection of radii for that particular isotach. 0 indicates do not use this radius, and 1 indicates use this radius and corresponding wind speed.
   * Columns 35-38 are the designated Rmax values computed for each of the quadrants selected for each particular isotach.
   * Column 39 is the Holland B parameter computed using the formulas outlined in the Holland paper, and implemented using the aswip program.
* The format of the file is fixed and users will want to use the aswip program to be sure that the input fort.22 file is properly formatted.
* The command line for NWS=19 is ``./aswip -n 19 -m 2 -z 1``

Example
^^^^^^^

The following is a simplified example of a fort.22 file with NWS=19 format:

.. code-block:: none

   AL, 09, 2017090500, 03, ASYM, 000, 167N, 0369W, 125, 948, XX, 34, NEQ, 150, 150, 090, 150, 1013, 027, , , , , , , , , IRMA, 1, 4, 0, 0, 1, 0, 000, 000, 025, 000, 1.25, 1.28, 1.26, 1.24, 1.30, 132, 129, 120, 133
   AL, 09, 2017090500, 03, ASYM, 000, 167N, 0369W, 125, 948, XX, 50, NEQ, 075, 075, 050, 075, 1013, 027, , , , , , , , , IRMA, 1, 4, 0, 1, 0, 0, 000, 027, 000, 000, 1.25, 1.28, 1.26, 1.24, 1.30, 132, 129, 120, 133
   AL, 09, 2017090500, 03, ASYM, 000, 167N, 0369W, 125, 948, XX, 64, NEQ, 035, 035, 025, 035, 1013, 027, , , , , , , , , IRMA, 1, 4, 1, 0, 0, 0, 027, 000, 000, 000, 1.25, 1.28, 1.26, 1.24, 1.30, 132, 129, 120, 133
   AL, 09, 2017090500, 03, ASYM, 000, 167N, 0369W, 125, 948, XX, 100, NEQ, 010, 010, 005, 010, 1013, 027, , , , , , , , , IRMA, 1, 4, 0, 0, 0, 1, 000, 000, 000, 027, 1.25, 1.28, 1.26, 1.24, 1.30, 132, 129, 120, 133
   AL, 09, 2017090506, 03, ASYM, 006, 171N, 0378W, 130, 943, XX, 34, NEQ, 160, 160, 100, 160, 1013, 025, , , , , , , , , IRMA, 2, 4, 0, 0, 1, 0, 000, 000, 025, 000, 1.30, 1.32, 1.29, 1.28, 1.33, 138, 135, 128, 140
   AL, 09, 2017090506, 03, ASYM, 006, 171N, 0378W, 130, 943, XX, 50, NEQ, 080, 080, 060, 080, 1013, 025, , , , , , , , , IRMA, 2, 4, 0, 1, 0, 0, 000, 025, 000, 000, 1.30, 1.32, 1.29, 1.28, 1.33, 138, 135, 128, 140
   AL, 09, 2017090506, 03, ASYM, 006, 171N, 0378W, 130, 943, XX, 64, NEQ, 040, 040, 030, 040, 1013, 025, , , , , , , , , IRMA, 2, 4, 1, 0, 0, 0, 025, 000, 000, 000, 1.30, 1.32, 1.29, 1.28, 1.33, 138, 135, 128, 140
   AL, 09, 2017090506, 03, ASYM, 006, 171N, 0378W, 130, 943, XX, 100, NEQ, 012, 012, 008, 012, 1013, 025, , , , , , , , , IRMA, 2, 4, 0, 0, 0, 1, 000, 000, 000, 025, 1.30, 1.32, 1.29, 1.28, 1.33, 138, 135, 128, 140

This example shows:

* Hurricane Irma (AL09 of 2017) with ASYM tech type (modified best track)
* Two time steps (0 and 6 hours) with 4 isotachs (34, 50, 64, 100 knots) at each time
* Radius values in NM for each quadrant (NE, SE, SW, NW) for each isotach
* For each isotach at each time step, different quadrants are selected:
   * 34kt isotach: using SW quadrant (column 33 = 1)
   * 50kt isotach: using SE quadrant (column 32 = 1)
   * 64kt isotach: using NE quadrant (column 31 = 1)
   * 100kt isotach: using NW quadrant (column 34 = 1)
* Rmax values computed for each quadrant based on selected isotachs
* Holland B parameter of 1.25-1.30 computed using aswip program

.. _nws-20:

NWS = 20
--------

In this format, ADCIRC uses the Generalized Asymmetric Holland Model (GAHM), which provides theoretical and practical improvements over previous parametric vortex models.

File Structure
^^^^^^^^^^^^^^

The file format is similar to the NWS = 19 format with 8 additional columns of data, see notes below for more information.

Format Specifics
^^^^^^^^^^^^^^^^

* The Generalized Asymmetric Holland Model (GAHM) provides a set of theoretical and practical improvements over previous parametric meteorological vortex models in ADCIRC. The theory and implementation of the GAHM was initially described at the 2013 ADCIRC Users Group Meeting.
* The NWS=20 option requires the following variables in the fort.22 file in a best track format, which includes the same variables as NWS=19 plus additional columns:
   * Forecast time in hours (column 6); enter the time in hours in each record starting at 0
   * Latitude of the eye (column 7)
   * Longitude of the eye (column 8)
   * Maximum sustained wind speed in knots (column 9)
   * Minimum sea level pressure in MB (column 10)
   * Wind intensity in knots of the radii defined in the record (34, 50, 64 or 100 knots) (column 12)
   * Radius of specified wind intensity for quadrants 1, 2, 3, 4 in NM (columns 14, 15, 16, 17); ≥ 0
   * Background pressure in MB (column 18); a standard value of 1013 can be used
   * Rmax as reported in the ATCF BEST TRACK file in column 20
   * Storm Name in Column 28 ATCF file format
   * Time Record number in column 29. There can be multiple lines for a given time record depending on the number of isotachs reported in the ATCF File
   * Number of isotachs reported in the ATCF file for the corresponding Time record.
   * Columns 31-34 indicate the selection of radii for that particular isotach. 0 indicates do not use this radius, and 1 indicates use this radius and corresponding wind speed.
   * Columns 35-38 are the designated Rmax values computed for each of the quadrants selected for each particular isotach.
   * Column 39 is the Holland B parameter computed using the formulas outlined in the Holland paper, and implemented using the aswip program.
   * Column 40-43 is the quadrant-varying Holland B parameter
   * Column 44-47 are the quadrant-varying Vmax calculated at the top of the planetary boundary (a wind reduction factor is applied to reduce the wind speed at the boundary to the 10-m surface)
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

Example
^^^^^^^

The following is a simplified example of a fort.22 file with NWS=20 format:

.. code-block:: none

   AL, 09, 2017090500, 03, GAHM, 000, 167N, 0369W, 125, 948, XX, 34, NEQ, 150, 150, 090, 150, 1013, 027, , , , , , , , , IRMA, 1, 4, 0, 0, 1, 0, 000, 000, 025, 000, 1.25, 1.28, 1.26, 1.24, 1.30, 132, 129, 120, 133
   AL, 09, 2017090500, 03, GAHM, 000, 167N, 0369W, 125, 948, XX, 50, NEQ, 075, 075, 050, 075, 1013, 027, , , , , , , , , IRMA, 1, 4, 0, 1, 0, 0, 000, 027, 000, 000, 1.25, 1.28, 1.26, 1.24, 1.30, 132, 129, 120, 133
   AL, 09, 2017090500, 03, GAHM, 000, 167N, 0369W, 125, 948, XX, 64, NEQ, 035, 035, 025, 035, 1013, 027, , , , , , , , , IRMA, 1, 4, 1, 0, 0, 0, 027, 000, 000, 000, 1.25, 1.28, 1.26, 1.24, 1.30, 132, 129, 120, 133
   AL, 09, 2017090500, 03, GAHM, 000, 167N, 0369W, 125, 948, XX, 100, NEQ, 010, 010, 005, 010, 1013, 027, , , , , , , , , IRMA, 1, 4, 0, 0, 0, 1, 000, 000, 000, 027, 1.25, 1.28, 1.26, 1.24, 1.30, 132, 129, 120, 133
   AL, 09, 2017090506, 03, GAHM, 006, 171N, 0378W, 130, 943, XX, 34, NEQ, 160, 160, 100, 160, 1013, 025, , , , , , , , , IRMA, 2, 4, 0, 0, 1, 0, 000, 000, 025, 000, 1.30, 1.32, 1.29, 1.28, 1.33, 138, 135, 128, 140
   AL, 09, 2017090506, 03, GAHM, 006, 171N, 0378W, 130, 943, XX, 50, NEQ, 080, 080, 060, 080, 1013, 025, , , , , , , , , IRMA, 2, 4, 0, 1, 0, 0, 000, 025, 000, 000, 1.30, 1.32, 1.29, 1.28, 1.33, 138, 135, 128, 140
   AL, 09, 2017090506, 03, GAHM, 006, 171N, 0378W, 130, 943, XX, 64, NEQ, 040, 040, 030, 040, 1013, 025, , , , , , , , , IRMA, 2, 4, 1, 0, 0, 0, 025, 000, 000, 000, 1.30, 1.32, 1.29, 1.28, 1.33, 138, 135, 128, 140
   AL, 09, 2017090506, 03, GAHM, 006, 171N, 0378W, 130, 943, XX, 100, NEQ, 012, 012, 008, 012, 1013, 025, , , , , , , , , IRMA, 2, 4, 0, 0, 0, 1, 000, 000, 000, 025, 1.30, 1.32, 1.29, 1.28, 1.33, 138, 135, 128, 140

This example shows:

* Hurricane Irma (AL09 of 2017) with GAHM tech type (modified best track)
* Two time steps (0 and 6 hours) with 4 isotachs (34, 50, 64, 100 knots) at each time
* Radius values in NM for each quadrant (NE, SE, SW, NW) for each isotach
* For each isotach at each time step, different quadrants are selected:
   * 34kt isotach: using SW quadrant (column 33 = 1)
   * 50kt isotach: using SE quadrant (column 32 = 1)
   * 64kt isotach: using NE quadrant (column 31 = 1)
   * 100kt isotach: using NW quadrant (column 34 = 1)
* Average Holland B parameter (column 39)
* Quadrant-varying Holland B parameters (columns 40-43)
* Quadrant-varying Vmax values (columns 44-47) at the top of the planetary boundary

Additional NWS Formats
----------------------

ADCIRC supports many other meteorological input formats, including:

* NWS = 7, 107: Tropical cyclone model - Generalized asymmetric Holland formulation with a moving reference frame
* NWS = 9, 109: Multiple tropical cyclone model - Generalized asymmetric Holland formulation
* NWS = 15: HWind format
* NWS = 16: GFDL format

These additional formats will be documented in future updates.
