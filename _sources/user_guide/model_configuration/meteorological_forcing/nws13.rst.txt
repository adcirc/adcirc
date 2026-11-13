.. meta::
   :description: NWS13 in ADCIRC
   :keywords: adcirc, nws13

.. _nws13:

NWS13
=====

Oceanweather Inc. (OWI) NetCDF Meteorological Input Files

Motivation
----------

A major source of error in ADCIRC can come from the accuracy and representation
of the wind and pressure fields. In order to allow ADCIRC access to high quality
and efficiently stored inputs the following features were added to the existing
:ref:`NWS12 <nws12>` functionality for the newer NetCDF version
:ref:`NWS13 <nws13>`.

.. figure:: /_static/images/user_guide/model_configuration/meteorological_forcing/nws13/overlay_domains.png
   :width: 400px


-  Moving storm-centered grids
-  Grids that can change size
-  Curvilinear grids
-  Arbitrary # of grid overlays
-  New ADCIRC moving interpolation code
-  Arbitrary & irregular timesteps

Some future-proofing attempts:

-  Multi-resolution representation of tropical cyclone wind fields
-  Be able to introduce additional meteorological parameters using the same
   format/NWS input

   -  Possibility to include wind stresses?
   -  Alternative reference heights?
   -  Ice information?
   -  Rainfall?

.. _netcdf_file_schema:

NetCDF File Schema
------------------

This description assumes some familiarity with NetCDF, HDF or other similar
self-describing binary file formats, and in particular features like "groups"
available in HDF5 and NetCDF4. Understanding of NetCDF conventions like Climate
Forecast (CF) may also be helpful. For more information about NetCDF please
visit the links below.

-  https://www.unidata.ucar.edu/software/netcdf/
-  https://www.hdfgroup.org/
-  https://cfconventions.org/
-  https://github.com/Unidata/netcdf4-python
-  https://github.com/pydata/xarray

.. _owi_nws13_convention:

OWI-NWS13 Convention
~~~~~~~~~~~~~~~~~~~~

**NetCDF4 File** contents:

-  Group(s) – 1 group per grid or overlay
-  Variables – U10/V10, PSFC
-  Dimensions – time, yi, xi
-  Attributes – grid rank/priority

**Global Attributes**:

-  *group_order*: space separated list of group names reflecting their
   appropriate rank/order
-  *conventions*: should include “OWI-NWS13“ (Climate Forecast conventions (CF)
   seemingly don’t support NC groups)

**Group Attributes**:

-  *rank*: integer representing the order of overlay/precedence in application
   to nodes

**Group Dimensions**:

-  *time* – length of time dimension
-  *yi* – number of rows in spatial grid description
-  *xi* – number of columns in spatial grid description

**Group Variables**:

-  *U10* (time, yi, xi), U-component of 10m WS (m/s)
-  *V10* (time, yi, xi), V-component of 10m WS (m/s)
-  *PSFC* (time, yi, xi), Surface Pressure (mb)
-  *lon* (yi, xi) or (time, yi, xi), Longitude in Decimal Degrees
-  *lat* (yi, xi) or (time, yi, xi), Latitude in Decimal Degrees
-  *time* (time), Datetime-number with units of “minutes from YYYY-mm-dd
   HH:MM:SS”
-  *clon* and *clat* (time), optional storm center coordinates for Powell drag

Example
~~~~~~~

| netcdf fort.22 {
|     // global attributes:
|     :group_order = "Main JPM0135" ;
|     :institution = "Oceanweather Inc. (OWI)" ;
|     :conventions = "CF-1.6 OWI-NWS13" ;

|  group: Main {
|      dimensions:
|          yi = 211 ;
|          xi = 221 ;
|          time = 133 ;
|      variables:
|          float lon(yi, xi) ;
|              lon:_FillValue = NaNf ;
|              lon:units = "degrees_east" ;
|              lon:standard_name = "longitude" ;
|              lon:axis = "X" ;
|              lon:coordinates = "time lat lon" ;
|          float lat(yi, xi) ;
|              lat:_FillValue = NaNf ;
|              lat:units = "degrees_north" ;
|              lat:standard_name = "latitude" ;
|              lat:axis = "Y" ;
|              lat:coordinates = "time lat lon" ;
|          int64 time(time) ;
|              time:units = "minutes since 1990-01-01T01:00:00" ;
|              time:calendar = "proleptic_gregorian" ;
|          float U10(time, yi, xi) ;
|              U10:_FillValue = NaNf ;
|              U10:units = "m s-1" ;
|              U10:coordinates = "time lat lon" ;
|          float V10(time, yi, xi) ;
|              V10:_FillValue = NaNf ;
|              V10:units = "m s-1" ;
|              V10:coordinates = "time lat lon" ;
|          float PSFC(time, yi, xi) ;
|              PSFC:_FillValue = NaNf ;
|              PSFC:units = "mb" ;
|              PSFC:coordinates = "time lat lon" ;
|  // group attributes:
|      :rank = 1 ;
|  } // group Main

|  group: JPM0135 {
|      dimensions:
|          time = 133 ;
|          yi = 501 ;
|          xi = 501 ;
|      variables:
|          int64 time(time) ;
|              time:units = "minutes since 1990-01-01T01:00:00" ;
|              time:calendar = "proleptic_gregorian" ;
|          float lat(time, yi, xi) ;
|              lat:_FillValue = NaNf ;
|              lat:units = "degrees_north" ;
|              lat:standard_name = "latitude" ;
|              lat:axis = "Y" ;
|              lat:coordinates = "time lat lon" ;
|          float lon(time, yi, xi) ;
|              lon:_FillValue = NaNf ;
|              lon:units = "degrees_east" ;
|              lon:standard_name = "longitude" ;
|              lon:axis = "X" ;
|              lon:coordinates = "time lat lon" ;
|          float U10(time, yi, xi) ;
|              U10:_FillValue = NaNf ;
|              U10:units = "m s-1" ;
|              U10:coordinates = "time lat lon" ;
|          float V10(time, yi, xi) ;
|              V10:_FillValue = NaNf ;
|              V10:units = "m s-1" ;
|              V10:coordinates = "time lat lon" ;
|          float PSFC(time, yi, xi) ;
|              PSFC:_FillValue = NaNf ;
|              PSFC:units = "mb" ;
|              PSFC:coordinates = "time lat lon" ;
|  // group attributes:
|      :rank = 2 ;
|  } // group JPM0135
| }

Notes
~~~~~

-  Decoupling the yi/xi dimensions from lat/lon allows lat and lon to be 2-d
   arrays by depending on both dimensions
-  Regular grids and Curvilinear grids
-  Non-evenly spaced grids as long as they can be expressed in a 2-d "mesh-grid"
-  Grids that change spatial resolution or position in time (but have consistent
   yi/xi array size)
-  Each group/sub-grid can define the timesteps independently, including start
   and stop times
-  Fill Value (and NetCDF packing/compression), and ieee nan floats supported

.. _fort.15_configuration:

fort.15 Configuration
---------------------

**WTIMINC**: configurable grid-to-mesh interpolation timestep in seconds,
separate from the input data timesteps. If the file contains storm-following
grids, there may be some utility in setting this finer than the input
timestep(s).

*&owiWindNetcdf* fort.15 namelist with all inputs quoted as strings

-  *NWS13ColdStartString*: required cold start time of simulation formatted as
   'YYYYMMDD.HHMMSS'
-  *NWS13WindMultiplier*: optional wind speed multiplier (**DWM** from
   `NWS12 <NWS12>`__)
-  *NWS13File*: optional file name for netCDF file ( *fort.22.nc* default )
-  *NWS13GroupForPowell*: optional group # to use for Powell drag

*example*:

 &owiWindNetcdf NWS13File='fort.22.nc' NWS13ColdStartString='20000706.000000' /

Suggestions
-----------

-  Keep the spatial and temporal resolutions to only what is required
-  Add sub-grids, more complex time dimensions, and complex lat/lon grid
   definitions only as necessary by sub-grid/group
-  In most cases 1-hourly specified--if you have a moving-grid overlay and can
   rely on the built-in moving-aware interpolation
-  Avoid discontinuities between overlays/grids
-  Best case: the lowest rank grid should just be a less resolved version of
   exact same higher resolution fields that will be overlaid
-  Lowest rank grid should probably cover the entire domain/mesh, and the entire
   model run time period

.. _interpolation_details:

Interpolation Details
---------------------

.. figure:: /_static/images/user_guide/model_configuration/meteorological_forcing/nws13/Ch_maxmin_interp_timeseries.png
   :width: 300px


.. figure:: /_static/images/user_guide/model_configuration/meteorological_forcing/nws13/Ch_movegrid_interp_spatial.png
   :width: 300px


The NWS13 met module implements a multi-step interpolation
procedure, primarily to minimize aliasing of the wind and pressure fields in
time and space.

The process occurs as follows:

.. _step_1_timespace_on_input_grids:

Step 1: Time/space on input grid(s)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Input fields are interpolated on the input grid (in time, and space if
necessary) to the intermediate timesteps specified in the *fort.15* parameter
**WTIMINC**.

.. _step_2_apply_intermediate_input_grid_to_mesh:

Step 2: Apply intermediate input grid to mesh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The gridded intermediate timesteps are interpolated in space on to the ADCIRC
mesh, appropriately overlaying multiple input grids by rank as supplied. At this
point the inputs are represented on the mesh at **WTIMINC** timesteps.

**NB:** For the lowest ranked overlay, the grid-to-mesh interpolation weights
are pre-computed, which means that this grid cannot be a storm-following moving
grid. For all additional overlays, the weights are dynamically updated.

.. _step_3_time_interpolation_to_model_timestep:

Step 3: Time interpolation to model timestep
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At each model computational timestep (**DTDP**) the inputs at intermediate
timesteps as represented on the mesh are interpolated in time to the appropriate
model timestep node-by-node.

.. figure:: /_static/images/user_guide/model_configuration/meteorological_forcing/nws13/nws13_interp1.png
   :width: 800px


.. figure:: /_static/images/user_guide/model_configuration/meteorological_forcing/nws13/nws13_interp2.png
   :width: 800px




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
