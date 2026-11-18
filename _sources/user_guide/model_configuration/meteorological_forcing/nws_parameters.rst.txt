.. _nws_parameter:

NWS Parameter
=============

.. contents::
   :local:
   :depth: 2


Introduction
------------

**NWS** is a parameter in the fort.15 file that selects the meteorological forcing input type. The value on the "NWS line" of the fort.15 file also implicitly includes other parameters affecting wave coupling and ice inputs. The NWS parameter controls not just the file type and handling of meteorological data, but also changes what the meteorological parameter line (informally, the ``WTIMINC`` line) looks like in the fort.15 file.

ADCIRC supports a wide range of meteorological input formats, including moving/fixed gridded data in several file formats, tropical cyclone track and parameter data that can be turned into wind/pressure fields via one of several internal vortex models, and direct specification of wind speeds or stresses on nodes.

.. note::

    For combined parameters, the NWS parameter represents:
        - NCICE value: first two digits (ice coverage data format)
        - NRS value: third digit (wave coupling mode)
        - NWS value: last two digits with sign (meteorological forcing mode)

Parameter Summary
-----------------

The following table is a summary of possible ``NWS`` values, their descriptions, and associated meteorological input files (required and optional).

.. note::

   Use of the Dynamic Asymmetric Holland Model (NWS=19) is discouraged. The Generalized Asymmetric Holland Model (NWS=20) provides improved functionality.


.. list-table:: NWS Parameter Summary
   :widths: 15 30 40 15
   :width: 100%
   :header-rows: 1
   :class: wrap-table, tight-table

   * - NWS Value
     - Short-name
     - Description
     - Required Input Files
   * - 0
     - No meteorological forcing
     - No wind, radiation stress or atmospheric pressure forcings are used.
     - n/a
   * - 1
     - Wind stress, every node, every timestep
     - Wind stress and atmospheric pressure are read in at all grid nodes at every model time step from the fort.22 file
     - fort.22
   * - 2
     - Wind stress, every node, every WTIMINC
     - Wind stress and atmospheric pressure are read in at all grid nodes at a time interval that does not equal the model time step from the fort.22 file. Interpolation in time is used to synchronize the wind and pressure information with the model time step.
     - fort.22
   * - 3
     - US Navy Fleet Numeric
     - Wind velocity is read in from a wind file from the fort.22 file in US Navy Fleet Numeric format. This information is interpolated in space onto the ADCIRC grid and in time to synchronize the wind and pressure information with the model time step.
     - fort.22
   * - 4
     - PBL/JAG
     - Wind velocity and atmospheric pressure are read in (PBL/JAG format) at selected ADCIRC grid nodes from the fort.22 file. Interpolation in time is used to synchronize the wind and pressure information with the model time step.
     - fort.22
   * - 5
     - Wind velocity, every node, every WTIMINC
     - Wind velocity and atmospheric pressure are read in at all grid nodes from the fort.22 File. Interpolation in time is used to synchronize the wind and pressure information with the model time step.
     - fort.22
   * - 6
     - Wind velocity, rectangular grid, every WTIMINC
     - Meteorological data (U,V,P) is input on a rectangular grid and interpolated in space onto the ADCIRC grid. The meteorological grid MUST cover the entire ADCIRC mesh.
     - fort.22
   * - 7
     - Wind stress, regular grid, every WTIMINC
     - Surface stress and pressure values are read in on a regular grid from the fort.22 file. Interpolation in time is used to synchronize the wind and pressure information with the model time step.
     - fort.22
   * - 8
     - Symmetric vortex model
     - Hurricane parameters are read in from the fort.22 file. Wind velocity and atmospheric pressure are calculated at every node on the fly by ADCIRC internally using the Dynamic Holland model.
     - fort.22
   * - 10
     - National Weather Service AVN
     - Wind velocity (10 m) and atmospheric pressure are read in from a sequence of National Weather Service (NWS) Aviation (AVN) model output files.
     - fort.200+
   * - 11
     - National Weather Service ETA 29km
     - Wind velocity (10 m) and atmospheric pressure are read in from a sequence of stripped down National Weather Service (NWS) ETA 29km model output files.
     - fort.200+
   * - 12
     - Oceanweather Inc (OWI)
     - Wind velocity (10 minute averaged winds at 10m) and atmospheric pressure are provided in the OWI format on one or two rectangular (lat/lon) grid(s).
     - fort.22, fort.221-224
   * - 13, -13
     - Ramped meteorological forcing
     - Similar to NWS=5/-5, but applies a ramping function to the meteorological inputs over a user-specified time period
     - fort.22
   * - 14, -14
     - GRIB2/NetCDF
     - Wind velocity and atmospheric pressure are read from standardized meteorological formats (GRIB2 or NetCDF)
     - NetCDF/GRIB2 files
   * - 15
     - HWind
     - HWind files are data assimilated snapshots of the wind velocity fields of tropical cyclones that are produced by the NOAA Hurricane Research Division (HRD).
     - fort.22 + HWind files
   * - 19
     - Dynamic Asymmetric Model (deprecated)
     - Wind velocity and atmospheric pressure are calculated directly coupled to ADCIRC at every time step using the asymmetric hurricane vortex formulation based on the Holland gradient wind model.
     - fort.22
   * - 20
     - Generalized Asymmetric Holland Model (GAHM)
     - The GAHM provides a set of theoretical and practical improvements over previous parametric meteorological vortex models in ADCIRC.
     - fort.22

Extended NWS with Ice + Waves
-----------------------------

The following table presents a summary of the extended ``NWS`` values to include ice-coverage and/or wind wave-coupling

.. list-table:: Extended NWS Values with Ice and Waves
   :widths: 30 10 15 15 15 15
   :header-rows: 1
   :class: wrap-table, tight-table

   * - Meteorological Data Format
     - Met. Only
     - Met. plus Waves from fort.23
     - Met. plus Waves SWAN
     - Met. plus Ice Coverage, Waves off
     - Met. plus Ice Coverage OWI-like format plus Waves from SWAN
   * - none
     - 0
     - n/a
     - n/a
     - n/a
     - n/a
   * - wind stress, every node, every timestep
     - 1
     - 101
     - 301
     - n/a
     - 12301
   * - wind stress, every node, every WTIMINC
     - 2
     - 102
     - 302
     - n/a
     - 12302
   * - US Navy Fleet Numeric
     - 3
     - 103
     - 303
     - n/a
     - 12303
   * - PBL/JAG
     - 4
     - 104
     - 304
     - n/a
     - 12304
   * - wind velocity, every node, every WTIMINC
     - 5
     - 105
     - 305
     - n/a
     - 12305
   * - wind velocity, rectangular grid, every WTIMINC
     - 6
     - 106
     - 306
     - n/a
     - 12306
   * - wind stress, regular grid, every WTIMINC
     - 7
     - 107
     - 307
     - n/a
     - 12307
   * - symmetrc vortex model
     - 8
     - 108
     - 308
     - n/a
     - 12308
   * - asymmetric vortex model (no longer available)
     - n/a
     - n/a
     - n/a
     - n/a
     - n/a
   * - National Weather Service AVN
     - 10
     - 110
     - 310
     - 10010
     - 12310
   * - National Weather Service ETA 29km
     - 11
     - 111
     - 311
     - n/a
     - 12311
   * - Oceanweather Inc (OWI)
     - 12
     - 112
     - 312
     - n/a
     - 12312
   * - GRIB2/NetCDF
     - 14
     - 114
     - 314
     - 14014
     - 14314
   * - H*Wind
     - 15
     - 115
     - 315
     - n/a
     - 12315
   * - Dynamic Asymmetric Holland Model (deprecated)
     - 19
     - 119
     - 319
     - n/a
     - 12319
   * - Generalized Asymmetric Holland Model
     - 20
     - 120
     - 320
     - n/a
     - 12320

Detailed NWS Parameter Descriptions
-----------------------------------

**NWS = 0**
    No wind, radiation stress or atmospheric pressure forcings are used.

**NWS = 1**
    Wind stress and atmospheric pressure are read in at all grid nodes at every model time step from the Single File Meteorological Forcing Input File.

**NWS = 2**
    Wind stress and atmospheric pressure are read in at all grid nodes at a time interval that does not equal the model time step from the Single File Meteorological Forcing Input File. Interpolation in time is used to synchronize the wind and pressure information with the model time step. The wind time interval (WTIMINC) is specified below.

**NWS = -2**
    Wind stress and atmospheric pressure are read in at all grid nodes at a time interval that does not equal the model time step from the Single File Meteorological Forcing Input File. It is assumed that the first entry in the Single File Meteorological Forcing Input File corresponds to the time that the current model run is started. Specifically, if the model is hot started, it is assumed that first entry in the Single File Meteorological Forcing Input File corresponds to the model hot start time. Succeeding entries in the Single File Meteorological Forcing Input File occur at the time interval (WTIMINC) specified below. Interpolation in time is used to synchronize the wind and pressure information with the model time step.

**NWS = 3**
    Wind velocity is read in from a wind file from the Single File Meteorological Forcing Input File in US Navy Fleet Numeric format. This information is interpolated in space onto the ADCIRC grid and in time to synchronize the wind and pressure information with the model time step. Garret's formula is used to compute wind stress from the wind velocity. Several parameters (IREFYR, IREFMO, IREFDAY, IREFHR, IREFMIN, REFSEC, NWLAT, NWLON, WLATMAX, WLONMIN, WLATINC, WLONINC, WTIMINC) describing the Fleet Numeric wind file must be specified below.

**NWS = 4**
    Wind velocity and atmospheric pressure are read in (PBL/JAG format) at selected ADCIRC grid nodes from the Single File Meteorological Forcing Input File. It is assumed that the first entry in the Single File Meteorological Forcing Input File corresponds to the beginning of the model run (e.g., the cold start time). Succeeding entries occur at the time interval (WTIMINC) specified below. Thus, if the model is hot started wind data must exist in the fort.22 file dating back to the beginning of the model run so that the model can find its appropriate place in the file. Interpolation in time is used to synchronize the wind and pressure information with the model time step. Garret's formula is used to compute wind stress from wind velocity.

**NWS = -4**
    Wind velocity and atmospheric pressure are read in (PBL/JAG format) at selected ADCIRC grid nodes from the Single File Meteorological Forcing Input File. It is assumed that the first entry in the Single File Meteorological Forcing Input File corresponds to the time that the current model run is started. Specifically, if the model is hot started, it is assumed that first entry in the Single File Meteorological Forcing Input File corresponds to the model hot start time. Succeeding entries in the Single File Meteorological Forcing Input File occur at the time interval (WTIMINC) specified below. Interpolation in time is used to synchronize the wind and pressure information with the model time step. Garret's formula is used to compute wind stress from wind velocity.

**NWS = 5**
    Wind velocity and atmospheric pressure are read in at all grid nodes from the Single File Meteorological Forcing Input File. It is assumed that the first entry in the Single File Meteorological Forcing Input File corresponds to the beginning of the model run (e.g., the cold start time). Succeeding entries occur at the time interval (WTIMINC) specified below. Thus, if the model is hot started wind data must exist in the Single File Meteorological Forcing Input File dating back to the beginning of the model run so that the model can find its appropriate place in the file. Interpolation in time is used to synchronize the wind and pressure information with the model time step. Garret's formula is used to compute wind stress from wind velocity.

**NWS = -5**
    Wind velocity and atmospheric pressure are read in at all grid nodes from the Single File Meteorological Forcing Input File. It is assumed that the first entry in the Single File Meteorological Forcing Input File corresponds to the time that the current model run is started. Specifically, if the model is hot started, it is assumed that first entry in the fort.22 file corresponds to the model hot start time. Succeeding entries in the Single File Meteorological Forcing Input File occur at the time interval (WTIMINC) specified below. Interpolation in time is used to synchronize the wind and pressure information with the model time step. Garret's formula is used to compute wind stress from wind velocity.

**NWS = 6**
    Wind velocity and atmospheric pressure are read in for a rectangular grid (either in Longitude, Latitude or Cartesian coordinates, consistent with the grid coordinates) from the Single File Meteorological Forcing Input File. This information is interpolated in space onto the ADCIRC grid and in time to synchronize the wind and pressure information with the model time step. Garret's formula is used to compute wind stress from the wind velocity. Several parameters describing the rectangular grid and time increment (NWLAT, NWLON, WLATMAX, WLONMIN, WLATINC, WLONINC, WTIMINC) must be specified below. The meterological grid MUST cover the entire ADCIRC mesh; that is, the ADCIRC mesh must be ENTIRELY within the meteorological grid or an error will result.

**NWS = 7**
    Surface stress and pressure values are read in on a regular grid from the Single File Meteorological Forcing Input File. It is assumed that the first entry in the Single File Meteorological Forcing Input File corresponds to the beginning of the model run (e.g., the cold start time). Succeeding entries occur at the time interval (WTIMINC) specified below. Thus, if the model is hot started wind data must exist in the Single File Meteorological Forcing Input File dating back to the beginning of the model run so that the model can find its appropriate place in the file. Interpolation in time is used to synchronize the wind and pressure information with the model time step.

**NWS = -7**
    Surface stress and pressure values are read in on a regular grid from the Single File Meteorological Forcing Input File. It is assumed that the first entry in the Single File Meteorological Forcing Input File corresponds to the time that the current model run is started. Specifically, if the model is hot started, it is assumed that first entry in the fort.22 file corresponds to the model hot start time. Succeeding entries in the Single File Meteorological Forcing Input File occur at the time interval (WTIMINC) specified below. Interpolation in time is used to synchronize the wind and pressure information with the model time step.

**NWS = 8**
    Hurricane parameters are read in from the Single File Meteorological Forcing Input File. Wind velocity and atmospheric pressure are calculated at every node on the fly by ADCIRC internally using the Dynamic Holland model. The input file is assumed to correspond to the ATCF Best Track/Objective Aid/Wind Radii Format. Historical tracks, real-time hindcast tracks and real-time forecast tracks may be found in this format. Selecting NWS = 8 also requires the specification of the cold start time, storm number, and boundary layer adjustment (see YYYY MM DD HH24 StormNumber BLAdj below). Garret's formula is used to compute wind stress from the wind velocity.

**NWS = 9**
    Asymmetric hurricane model, no longer supported

**NWS = 10**
    Wind velocity (10 m) and atmospheric pressure are read in from a sequence of National Weather Service (NWS) Aviation (AVN) model output files. Each AVN file is assumed to contain data on a Gaussian longitude, latitude grid at a single time. Consecutive files in the sequence are separated by N hours in time (where N=WTIMINC/3600 and WTIMINC is read in below). The files are named using the convention: 
    
    * fort.200 – wind & pressure at the time of a model hot start (this file is not used for a cold start)
    * fort.XX1 (where XX1=200+1*N) – wind & pressure N hours after a cold or hot start
    * fort.XX2 (where XX2=200+2*N) – wind & pressure 2N hours after a cold or hot start
    * fort.XX3 (where XX3=200+3*N) – wind & pressure 3N hours after a cold or hot start and so on for all meteorological files
    
    Prior to ADCIRC version 34.05 these files were in binary and created from a larger Grib form file using the program UNPKGRB1. Starting with ADCIRC version 34.05, the files are in ASCII tabular format. If ADCIRC is hot started, it must be done at an even N hour interval so that the hot start time corresponds to the time of a meteorological file. Enough meteorological files must be present to extend through the ending time of the model run. Garret's formula is used to compute wind stress from the wind velocity.

**NWS = 11**
    Wind velocity (10 m) and atmospheric pressure are read in from a sequence of stripped down National Weather Service (NWS) ETA 29km model output files. Each ETA file is assumed to contain data on an E grid for a single day (8 data sets, one every 3 hours, beginning @ 03:00 and continuing through 24:00 of the given day). The files are named using the convention: 
    
    * fort.200 – wind & pressure the day before a model run is hot started. The final data in this file are used as the initial met condition for the hot start. This file is not used for a cold start.
    * fort.201 – wind & pressure during the first day after a cold or hot start.
    * fort.202 – wind & pressure during the second day after a cold or hot start
    * fort.203 – wind & pressure during the third day after a cold or hot start
    
    This sequence continues for all meteorological files. These files are in binary and have the format described below. The wind data is converted to an east-west, north-south coordinate system inside ADCIRC. If the model is hot started, it must be done at an even day interval so that the hot start time corresponds to the time of a meteorological file. Enough meteorological files must be present to extend through the ending time of the model run. Garret's formula is used to compute wind stress from the wind velocity.

**NWS = 12**
    Wind velocity (10 minute averaged winds at 10m) and atmospheric pressure are provided in the OWI format on one or two rectangular (lat/lon) grid(s). If two grids are used, the first is designated as the large ("basin") scale grid, and the second is designated as the small ("region") scale grid. The Single File Meteorological Forcing Input File (fort.22) is only used to specify a few configuration parameters, while the actual wind fields are recorded in files named:
    
    * fort.221
    * fort.222
    * fort.223 (optional)
    * fort.224 (optional)
    
    The time increment of the meteorological forcing is specified through WTIMINC in the fort.15 file. The wind and pressure fields are interpolated in space onto the ADCIRC grid and in time to synchronize the wind and pressure information with the model time step. Garret's formula is used to compute wind stress from wind velocity.

**NWS = 13, -13**
    Similar to NWS = 5/-5, but applies a ramping function to the meteorological inputs. Wind velocity and atmospheric pressure are read in at all grid nodes from the Single File Meteorological Forcing Input File. The positive/negative convention for time reference is the same as with NWS = 5/-5. The meteorological forcing is ramped up over a user-specified time period using the WRAMP parameter.

**NWS = 14, -14**
    Wind velocity and atmospheric pressure are read in from NetCDF or GRIB2 format files. This option allows ADCIRC to directly read standardized meteorological formats. If NWS = 14, ADCIRC assumes data from the cold start time; if NWS = -14, ADCIRC assumes data from the hot start time. The files must include wind velocity components and sea level pressure on a regular grid with proper metadata.

**NWS = 15, -15**
    HWind files are data assimilated snapshots of the wind velocity fields of tropical cyclones that are produced by the NOAA Hurricane Research Division (HRD). 
    
    * If the NWS value is set to +15, the hours column in the associated meterological forcing input file (fort.22) is relative to the cold start time. 
    * If NWS is set to -15, that hours column is relative to the hot start time. 
    
    Please see the documentation of the Single File Meteorological Forcing Input File for complete details.

**NWS = 19**
    User has the ability to select which Isotach to use in each of the 4 quadrants. User also has ability to modify RMAX and Holland's B parameter using the ASWIP program. The auxiliary preprocessing program ASWIP.F (located in the /wind directory and executable is created by typing, make aswip, in the work folder after adcirc executable has been generated), will generate the fort.22 input file for NWS=19 from a NWS=9 formatted input file.

    Hurricane parameters are read in from the Single File Meteorological Forcing Input File. It is assumed that the line in the Single File Meteorological Forcing Input File with a zero as the forecast increment (i.e., column 6) corresponds to the start of the current simulation run, whether it is a hotstart or cold start. In other words, there is no option to set the NWS value negative to indicate that the file starts at the ADCIRC hotstart time. Rather, the forecast increment in hours (column 6) is used to indicate the relationship between the ADCIRC time and the data in the fort.22 file. 
    
    Wind velocity and atmospheric pressure are calculated at exact finite element mesh node locations and directly coupled to ADCIRC at every time step using the asymmetric hurricane vortex formulation (Mattocks et al, 2006; Mattocks and Forbes, 2008) based on the Holland gradient wind model. The input file is assumed to correspond to the ATCF Best Track/Objective Aid/Wind Radii Format. Historical tracks, real-time hindcast tracks and real-time forecast tracks may be found in this format. This option uses the radii at specific wind speeds (34, 50, 64, 100 knots) reported in the four quadrants (NE, SE, SW, NW) of the storm to calculate the radius of maximum winds as a function of the azimuthal angle. Garret's formula is used to compute wind stress from the wind velocity. 
    
    The NWS=19 option allows the user to set a value for Rmax and Holland B Parameter. Additionally the user can select the isotachs to be used for each of the 4 quadrants. The utility program aswip_1.0.3.F located in the /wind folder will generate the NWS=19 fomatted file from a NWS=9 formatted fort.22 input file.

    In order to use the NWS=19 option, the file needs to be in best track format. The forecast period (column #6) needs to be edited to reflect the time of the forecast/nowcast for each track location (each line) in hours from the start of the simulation (0, 6, 12, 18, etc). There is no -19 option to indicate that the hours in column 6 are relative to the hotstart time. For the dynamic asymmetric model (NWS=19), ADCIRC always assumes that hour 0 corresponds to when the model is started, whether that is a cold start or a hot start. Therefore, ADCIRC analysts should not attempt to set NWS to -19.  The original data in that column depends on what type of best track format data is being used. The original data might have 0 or other numbers in that column. See: http://www.nrlmry.navy.mil/atcf_web/docs/database/new/abrdeck.html
    
    It is suggested that users change the "BEST" tech type to "ASYM" in column 5 in the fort.22 file to denote that the file has been modified to accommodate the asymmetric wind formulation (the simulation time in hours in the 6th column has been added, etc.) so it will not get confused in the future with a best track file.

    The NWS=19 option requires the following variables in the fort.22 file in a best track format:

    1. Forecast time in hours (column 6); enter the time in hours in each record starting at 0
    2. Latitude of the eye (column 7)
    3. Longitude of the eye (column 8)
    4. Maximum sustained wind speed in knots (column 9)
    5. Minimum sea level pressure in MB (column 10)
    6. Wind intensity in knots of the radii defined in the record (34, 50, 64 or 100 knots) (column 12)
    7. Radius of specified wind intensity for quadrants 1, 2, 3, 4 in NM (columns 14, 15, 16, 17); ? 0
    8. Background pressure in MB (column 18); a standard value of 1013 can be used
    9. Rmax as reported in the ATCF BEST TRACK file in column 20
    10. Storm Name in Column 28 ATCF file format
    11. Time Record number in column 29. There can be multiple lines for a given time record depending on the number of isotachs reported in the ATCF File
    12. Number of isotachs reported in the ATCF file for the corresponding Time record
    13. Columns 31-34 indicate the selection of radii for that particular isotach. 0 indicates do not use this radius, and 1 indicates use this radius and corresponding wind speed
    14. Columns 35-38 are the designated Rmax values computed for each of the quadrants selected for each particular isotach
    15. Column 39 is the Holland B parameter computed using the formulas outlines in the Holland paper, and implemented using the aswip program

    The format of the file is fixed and users will want to use the aswip program to be sure that the input fort.22 file is properly formatted.

**NWS = 20**
    Generalized Asymmetric Holland Model (GAHM) Format. The track file format is similar to that of the older Dynamic Asymmetric Model (NWS = 19) but with 8 additional columns of data. See notes in the fort.22 file for more information. The theory and implementation of the GAHM was initially described at the 2013 ADCIRC Users Group Meeting.

**NWS = 100, 101, 102, -102, 103, 104, -104, 105, -105, 106, 110, 111**
    Wave radiation stress is applied in addition to meteorological forcing. The meteorological input is specified by: SIGN(NWS)*(ABS(NWS)-100). For example:
    
    * NWS=100 means include wave radiation stress with no meteorological forcing (NWS=0)
    * NWS=101 means include wave radiation stress plus meteorological forcing corresponding to NWS=1
    * NWS=-104 means include wave radiation stress plus meteorological forcing corresponding to NWS=-4, etc.
    
    Wave radiation stress is read from a Wave Radiation Stress Forcing File. The format of this file is similar to the generic meteorological forcing file when NWS=-4 with the exception that no pressure values are read in. The time increment between consecutive radiation stress fields (RSTIMINC) is specified below.

**NWS = 300, 301, 302, -302, 303, 304, -304, 305, -305, 306, 310, 311, 312, -312**
    NWS values in the 300's indicate a SWAN+ADCIRC run. Note: padcswan or adcswan must be precompiled to use this option.

    The SWAN wave model is dynamically coupled to the ADCIRC model. Radiation stresses and currents from the SWAN model are applied in addition to meteorological forcing. The meteorological input is specified by: SIGN(NWS)*(ABS(NWS)-300). For example:
    
    * NWS=300 means include wave radiation stress with no meteorological forcing (NWS=0)
    * NWS=301 means include wave radiation stress plus meteorological forcing corresponding to NWS=1
    * NWS=-304 means include wave radiation stress plus meteorological forcing corresponding to NWS=-4, etc.
    
    Wave radiation stress are computed by the SWAN model every RSTIMINC seconds and passed into ADCIRC. In addition to assigning RSTIMINC the user must have a SWAN input and control file (fort.26) in the same working directory as the fort.15 ADCIRC control file.


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
