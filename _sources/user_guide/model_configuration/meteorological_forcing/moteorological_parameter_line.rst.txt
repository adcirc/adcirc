.. meta::
   :description: Supplemental meteorological/wave/ice parameters in ADCIRC
   :keywords: adcirc, supplemental meteorological/wave/ice parameters

.. _supplemental_meteorological_wave_ice_parameters:

Supplemental meteorological/wave/ice parameters
===============================================

This table is helpful for understanding file requirements and how the
meteorological parameter line (informally, the
:ref:`WTIMINC` line) should look in the :ref:`fort.15
file <fort15>`. These are principally determined by the value of the
:ref:`NWS <nws_parameter>` line (also in the fort.15 file), though note that the "NWS"
values below are for the full-length value 
in the fort.15 file. Useful information is also contained in the :ref:`fort.22 file
format <fort22>` and :ref:`wind stress <wind_stress>` pages.

.. list-table::
   :header-rows: 1
   :widths: 15 10 10 18 14 8
   :class: wrap-table

   * - Meteorological Data Format
     - :ref:`Forcing <supp_forcing_abbreviations>`
     - :ref:`NWS <nws_parameter>` Value
     - :ref:`WTIMINC <WTIMINC>` Line
     - :ref:`Requirements <supp_forcing_requirements>`
     - :ref:`Notes <supp_forcing_notes>`
   * - **none**
     - nonexistent
     - 0
     - nonexistent
     - none
     - none
   * - **wind stress, every node, every timestep**
     - m
     - 1
     - nonexistent
     - f22
     - [1]_
   * - 
     - m, rs
     - 101
     - ``RSTIMINC``
     - f22, f23
     - [2]_
   * - 
     - m, sw
     - 301
     - ``RSTIMINC``
     - f22, A+S
     - [3]_
   * - 
     - m, st
     - 401
     - ``RSTIMINC``
     - f22, STW
     - [3]_
   * - **wind stress, every node, every ``WTIMINC``**
     - m
     - 2
     - ``WTIMINC``
     - f22
     - [1]_
   * - 
     - m, rs
     - 102
     - ``WTIMINC`` ``RSTIMINC``
     - f22, f23
     - [2]_
   * - 
     - m, sw
     - 302
     - ``WTIMINC`` ``RSTIMINC``
     - f22, A+S
     - [3]_
   * - 
     - m, st
     - 402
     - ``WTIMINC`` ``RSTIMINC``
     - f22, STW
     - [3]_
   * - **US Navy Fleet Numeric**
     - m
     - 3
     - ``irefyr`` ``irefmo`` ``irefday`` ``irefhr`` ``irefmin`` ``refsec``
     - f22
     - none
   * - 
     - m, rs
     - 103
     - ``irefyr`` ``irefmo`` ``irefday`` ``irefhr`` ``irefmin`` ``refsec`` ``RSTIMINC``
     - f22, f23
     - [2]_
   * - 
     - m, sw
     - 303
     - ``irefyr`` ``irefmo`` ``irefday`` ``irefhr`` ``irefmin`` ``refsec`` ``RSTIMINC``
     - f22, A+S
     - [3]_
   * - 
     - m, st
     - 403
     - ``irefyr`` ``irefmo`` ``irefday`` ``irefhr`` ``irefmin`` ``refsec`` ``RSTIMINC``
     - f22, STW
     - [3]_
   * - 
     - m, ic
     - 12003
     - ``irefyr`` ``irefmo`` ``irefday`` ``irefhr`` ``irefmin`` ``refsec`` ``cice_timinc``
     - f22, f25
     - none
   * - 
     - m, rs, ic
     - 12103
     - ``irefyr`` ``irefmo`` ``irefday`` ``irefhr`` ``irefmin`` ``refsec`` ``RSTIMINC`` ``cice_timinc``
     - f22, f23, f25
     - [2]_
   * - 
     - m, sw, ic
     - 12303
     - ``irefyr`` ``irefmo`` ``irefday`` ``irefhr`` ``irefmin`` ``refsec`` ``RSTIMINC`` ``cice_timinc``
     - f22, A+S, f25
     - [3]_
   * - 
     - m, st, ic
     - 12403
     - ``irefyr`` ``irefmo`` ``irefday`` ``irefhr`` ``irefmin`` ``refsec`` ``RSTIMINC`` ``cice_timinc``
     - f22, STW, f25
     - [3]_
   * - **PBL/JAG**
     - m
     - 4
     - ``WTIMINC``
     - f22
     - none
   * - 
     - m, rs
     - 104
     - ``WTIMINC`` ``RSTIMINC``
     - f22, f23
     - [2]_
   * - 
     - m, sw
     - 304
     - ``WTIMINC`` ``RSTIMINC``
     - f22, A+S
     - [3]_
   * - 
     - m, st
     - 404
     - ``WTIMINC`` ``RSTIMINC``
     - f22, STW
     - [3]_
   * - 
     - m, ic
     - 12004
     - ``WTIMINC`` ``cice_timinc``
     - f22, f25
     - none
   * - 
     - m, rs, ic
     - 12104
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - f22, f23, f25
     - [2]_
   * - 
     - m, sw, ic
     - 12304
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - f22, A+S, f25
     - [3]_
   * - 
     - m, st, ic
     - 12404
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - f22, STW, f25
     - [3]_
   * - **wind velocity, every node, every ``WTIMINC``**
     - m
     - 5
     - ``WTIMINC``
     - f22
     - none
   * - 
     - m, rs
     - 105
     - ``WTIMINC`` ``RSTIMINC``
     - f22, f23
     - [2]_
   * - 
     - m, sw
     - 305
     - ``WTIMINC`` ``RSTIMINC``
     - f22, A+S
     - [3]_
   * - 
     - m, st
     - 405
     - ``WTIMINC`` ``RSTIMINC``
     - f22, STW
     - [3]_
   * - 
     - m, ic
     - 12005
     - ``WTIMINC`` ``cice_timinc``
     - f22, f25
     - none
   * - 
     - m, rs, ic
     - 12105
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - f22, f23, f25
     - [2]_
   * - 
     - m, sw, ic
     - 12305
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - f22, A+S, f25
     - [3]_
   * - 
     - m, st, ic
     - 12405
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - f22, STW, f25
     - [3]_
   * - **wind velocity, rectangular grid, every ``WTIMINC``**
     - m
     - 6
     - ``NWLAT NWLON WLATMAX WLONMIN WLATINC WLONINC WTIMINC``
     - f22
     - none
   * - 
     - m, rs
     - 106
     - ``NWLAT NWLON WLATMAX WLONMIN WLATINC WLONINC WTIMINC RSTIMINC``
     - f22, f23
     - [2]_
   * - 
     - m, sw
     - 306
     - ``NWLAT NWLON WLATMAX WLONMIN WLATINC WLONINC WTIMINC RSTIMINC``
     - f22, A+S
     - [3]_
   * - 
     - m, st
     - 406
     - ``NWLAT NWLON WLATMAX WLONMIN WLATINC WLONINC WTIMINC RSTIMINC``
     - f22, STW
     - [3]_
   * - 
     - m, ic
     - 12006
     - ``NWLAT NWLON WLATMAX WLONMIN WLATINC WLONINC WTIMINC cice_timinc``
     - f22, f25
     - none
   * - 
     - m, rs, ic
     - 12106
     - ``NWLAT NWLON WLATMAX WLONMIN WLATINC WLONINC WTIMINC RSTIMINC cice_timinc``
     - f22, f23, f25
     - [2]_
   * - 
     - m, sw, ic
     - 12306
     - ``NWLAT NWLON WLATMAX WLONMIN WLATINC WLONINC WTIMINC RSTIMINC cice_timinc``
     - f22, A+S, f25
     - [3]_
   * - 
     - m, st ,ic
     - 12406
     - ``NWLAT NWLON WLATMAX WLONMIN WLATINC WLONINC WTIMINC RSTIMINC cice_timinc``
     - f22, STW, f25
     - [3]_
   * - **wind stress, regular grid, every ``WTIMINC``**
     - m
     - 7
     - ``WTIMINC``
     - f22
     - [1]_
   * - 
     - m, rs
     - 107
     - ``WTIMINC`` ``RSTIMINC``
     - f22, f23
     - [2]_
   * - 
     - m, sw
     - 307
     - ``WTIMINC`` ``RSTIMINC``
     - f22, A+S
     - [3]_
   * - 
     - m, st
     - 407
     - ``WTIMINC`` ``RSTIMINC``
     - f22, STW
     - [3]_
   * - **symmetric vortex model**
     - m
     - 8
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj``
     - f22
     - none
   * - 
     - m, rs
     - 108
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``RSTIMINC``
     - f22, f23
     - [2]_
   * - 
     - m, sw
     - 308
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``RSTIMINC``
     - f22, A+S
     - [3]_
   * - 
     - m, st
     - 408
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``RSTIMINC``
     - f22, STW
     - [3]_
   * - 
     - m, ic
     - 12008
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``cice_timinc``
     - f22, f25
     - none
   * - 
     - m, rs, ic
     - 12108
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``RSTIMINC`` ``cice_timinc``
     - f22, f23, f25
     - [2]_
   * - 
     - m, sw, ic
     - 12308
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``RSTIMINC`` ``cice_timinc``
     - f22, A+S, f25
     - [3]_
   * - 
     - m, st, ic
     - 12408
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``RSTIMINC`` ``cice_timinc``
     - f22, STW, f25
     - [3]_
   * - **asymmetric vortex model (no longer available)**
     - n/a
     - 9
     - n/a
     - n/a
     - none
   * - **National Climatic Data Center GFS**
     - m
     - 10
     - ``WTIMINC``
     - f2xx+
     - none
   * - 
     - m, rs
     - 110
     - ``WTIMINC`` ``RSTIMINC``
     - f2xx+, f23
     - [2]_
   * - 
     - m, sw
     - 310
     - ``WTIMINC`` ``RSTIMINC``
     - f2xx+, A+S
     - [3]_
   * - 
     - m, st
     - 410
     - ``WTIMINC`` ``RSTIMINC``
     - f2xx+, STW
     - [3]_
   * - 
     - m, ic
     - 12010
     - ``WTIMINC`` ``cice_timinc``
     - f2xx+, f25
     - none
   * - 
     - m, rs, ic
     - 12110
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - f2xx+, f23, f25
     - [2]_
   * - 
     - m, sw, ic
     - 12310
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - f2xx+, A+S, f25
     - [3]_
   * - 
     - m, st, ic
     - 12410
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - f2xx+, STW, f25
     - [3]_
   * - **National Weather Service ETA 29km**
     - m
     - 11
     - nonexistent
     - f22
     - none
   * - 
     - m, rs
     - 111
     - ``RSTIMINC``
     - f22, f23
     - [2]_
   * - 
     - m, sw
     - 311
     - ``RSTIMINC``
     - f22, A+S
     - [3]_
   * - 
     - m, st
     - 411
     - ``RSTIMINC``
     - f22, STW
     - [3]_
   * - 
     - m,ic
     - 12011
     - ``cice_timinc``
     - f22, f25
     - none
   * - 
     - m, rs, ic
     - 12111
     - ``RSTIMINC`` ``cice_timinc``
     - f22, f23, f25
     - [2]_
   * - 
     - m, sw, ic
     - 12311
     - ``RSTIMINC`` ``cice_timinc``
     - f22, A+S, f25
     - [3]_
   * - 
     - m, st, ic
     - 12411
     - ``RSTIMINC`` ``cice_timinc``
     - f22, STW, f25
     - [3]_
   * - **Oceanweather Inc (OWI)**
     - m
     - 12
     - ``WTIMINC``
     - f22x
     - none
   * - 
     - m, rs
     - 112
     - ``WTIMINC`` ``RSTIMINC``
     - f22x, f23
     - [2]_
   * - 
     - m, sw
     - 312
     - ``WTIMINC`` ``RSTIMINC``
     - f22x, A+S
     - [3]_
   * - 
     - m, st
     - 412
     - ``WTIMINC`` ``RSTIMINC``
     - f22x, STW
     - [3]_
   * - 
     - m, ic
     - 12012
     - ``WTIMINC`` ``cice_timinc``
     - f22x, f25
     - none
   * - 
     - m, rs, ic
     - 12112
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - f22x, f23, f25
     - [2]_
   * - 
     - m, sw, ic
     - 12312
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - f22x, A+S, f25
     - [3]_
   * - 
     - m, st, ic
     - 12412
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - f22x, STW, f25
     - [3]_
   * - **Oceanweather Inc (OWI) NetCDF**
     - m
     - 13
     - ``WTIMINC``
     - f22nc
     - `NWS13 <NWS13>`__
   * - 
     - m, rs
     - 113
     - ``WTIMINC`` ``RSTIMINC``
     - f22nc, f23
     - [2]_
   * - 
     - m, sw
     - 313
     - ``WTIMINC`` ``RSTIMINC``
     - f22nc, A+S
     - [3]_
   * - 
     - m, st
     - 413
     - ``WTIMINC`` ``RSTIMINC``
     - f22nc, STW
     - [3]_
   * - 
     - m, ic
     - 12013
     - ``WTIMINC`` ``cice_timinc``
     - f22nc, f25
     - none
   * - 
     - m, rs, ic
     - 12113
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - f22nc, f23, f25
     - [2]_
   * - 
     - m, sw, ic
     - 12313
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - f22nc, A+S, f25
     - [3]_
   * - 
     - m, st, ic
     - 12413
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - f22nc, STW, f25
     - [3]_
   * - **H*Wind**
     - m
     - 15
     - ``WTIMINC``
     - hwind+,f22
     - none
   * - 
     - m, rs
     - 115
     - ``WTIMINC`` ``RSTIMINC``
     - hwind+,f22, f23
     - [2]_
   * - 
     - m, sw
     - 315
     - ``WTIMINC`` ``RSTIMINC``
     - hwind+,f22, A+S
     - [3]_
   * - 
     - m, st
     - 415
     - ``WTIMINC`` ``RSTIMINC``
     - hwind+,f22, STW
     - [3]_
   * - 
     - m, ic
     - 12015
     - ``WTIMINC`` ``cice_timinc``
     - hwind+,f22, f25
     - none
   * - 
     - m, rs, ic
     - 12115
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - hwind+,f22, f23, f25
     - [2]_
   * - 
     - m, sw, ic
     - 12315
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - hwind+,f22, A+S, f25
     - [3]_
   * - 
     - m, st, ic
     - 12415
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - hwind+,f22, STW, f25
     - [3]_
   * - **Geophysical Fluid Dynamics Laboratory (GFDL) Model**
     - m
     - 16
     - ``WTIMINC``
     - gfdl+,f22
     - none
   * - 
     - m, rs
     - 116
     - ``WTIMINC`` ``RSTIMINC``
     - gfdl+,f22, f23
     - [2]_
   * - 
     - m, sw
     - 316
     - ``WTIMINC`` ``RSTIMINC``
     - gfdl+,f22, A+S
     - [3]_
   * - 
     - m, st
     - 416
     - ``WTIMINC`` ``RSTIMINC``
     - gfdl+,f22, STW
     - [3]_
   * - 
     - m, ic
     - 12016
     - ``WTIMINC`` ``cice_timinc``
     - gfdl+,f22, f25
     - none
   * - 
     - m, rs, ic
     - 12116
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - gfdl+,f22, f23, f25
     - [2]_
   * - 
     - m, sw, ic
     - 12316
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - gfdl+,f22, A+S, f25
     - [3]_
   * - 
     - m, st, ic
     - 12416
     - ``WTIMINC`` ``RSTIMINC`` ``cice_timinc``
     - gfdl+,f22, STW, f25
     - [3]_
   * - **Dynamic Asymmetric Model**
     - m
     - 19
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj``
     - pp,f22
     - none
   * - 
     - m, rs
     - 119
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``RSTIMINC``
     - pp, f22, f23
     - [2]_
   * - 
     - m, sw
     - 319
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``RSTIMINC``
     - pp, f22, A+S
     - [3]_
   * - 
     - m, st
     - 419
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``RSTIMINC``
     - pp, f22, STW
     - [3]_
   * - 
     - m, ic
     - 12019
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``cice_timinc``
     - pp, f22, f25
     - none
   * - 
     - m, rs, ic
     - 12119
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``RSTIMINC`` ``cice_timinc``
     - pp, f22, f23, f25
     - [2]_
   * - 
     - m, sw, ic
     - 12319
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``RSTIMINC`` ``cice_timinc``
     - pp, f22, A+S, f25
     - [3]_
   * - 
     - m, st, ic
     - 12419
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``RSTIMINC`` ``cice_timinc``
     - pp, f22, STW, f25
     - [3]_
   * - **Generalized Asymmetric Holland Model (GAHM)**
     - m
     - 20
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``geofactor``
     - pp, f22
     - none
   * - 
     - m, rs
     - 120
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``geofactor`` ``RSTIMINC``
     - pp, f22, f23
     - [2]_
   * - 
     - m, sw
     - 320
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``geofactor`` ``RSTIMINC``
     - pp, f22, A+S
     - [3]_
   * - 
     - m, st
     - 420
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``geofactor`` ``RSTIMINC``
     - pp, f22, STW
     - [3]_
   * - 
     - m, ic
     - 12020
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``geofactor`` ``cice_timinc``
     - pp, f22, f25
     - none
   * - 
     - m, rs, ic
     - 12120
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``geofactor`` ``RSTIMINC`` ``cice_timinc``
     - pp, f22, f23, f25
     - [2]_
   * - 
     - m, sw, ic
     - 12320
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``geofactor`` ``RSTIMINC`` ``cice_timinc``
     - pp, f22, A+S, f25
     - [3]_
   * - 
     - m, st, ic
     - 12420
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``geofactor`` ``RSTIMINC`` ``cice_timinc``
     - pp, f22, STW, f25
     - [3]_
   * - **Blended GAHM and OWI**
     - m
     - 30
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``geofactor`` ``pureVortex`` ``pureBackground``
     - pp, f22, nwsf22
     - none
   * - 
     - m, rs
     - 130
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``geofactor`` ``RSTIMINC`` ``pureVortex`` ``pureBackground``
     - pp, f22, nwsf22, f23
     - [2]_
   * - 
     - m, sw
     - 330
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``geofactor`` ``RSTIMINC`` ``pureVortex`` ``pureBackground``
     - pp, f22, nwsf22, A+S
     - [3]_
   * - 
     - m, st
     - 430
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``geofactor`` ``RSTIMINC`` ``pureVortex`` ``pureBackground``
     - pp, f22, nwsf22, STW
     - [3]_
   * - 
     - m, ic
     - 12030
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``geofactor`` ``cice_timinc`` ``pureVortex`` ``pureBackground``
     - pp, f22, nwsf22, f25
     - none
   * - 
     - m, rs, ic
     - 12130
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``geofactor`` ``RSTIMINC`` ``cice_timinc`` ``pureVortex`` ``pureBackground``
     - pp, f22, nwsf22, f23, f25
     - [2]_
   * - 
     - m, sw, ic
     - 12330
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``geofactor`` ``RSTIMINC`` ``cice_timinc`` ``pureVortex`` ``pureBackground``
     - pp, f22, nwsf22, A+S, f25
     - [3]_
   * - 
     - m, st, ic
     - 12430
     - ``YYYY`` ``MM`` ``DD`` ``HH24`` ``StormNumber`` ``BLAdj`` ``geofactor`` ``RSTIMINC`` ``cice_timinc`` ``pureVortex`` ``pureBackground``
     - pp, f22, nwsf22, STW, f25
     - [3]_

.. _supp_forcing_abbreviations:

Forcing Abbreviations
---------------------

* **m** - meteorology
* **rs** - Wave radiation stress in an analyst-supplied `fort.23 <fort23>` input file
* **sw** - Wave coupling using a simultaneously executing coupled SWAN model
* **st** - Wave coupling using a simultaneously executing ESMF coupled STWAVE model
* **ic** - Ice coverage using an analyst-supplied ice coverage data files

.. _supp_forcing_requirements:

Requirements
------------

* **f22** - Requires analyst-supplied :ref:`fort.22 <fort22>` meteorological forcing file.
* **f2xx+** - Requires at least two fort.2xx meteorological forcing files from
   the GFS model. See fort.22 documentation for file naming and formatting
   details.
* **f22x** - Requires datasets of pressure and wind data in OWI formatted
   :ref:`fort.221 <fort22>` and :ref:`fort.222 <fort22>` files respectively,
   optionally including nested wind and pressure in :ref:`fort.223 <fort22>` and
   :ref:`fort.224 <fort22>` files. See :ref:`fort.22 <fort22>` documentation for details.
* **hwind+** - Requires two or more meteorological forcing files from the
   H*Wind model. See :ref:`fort.22 <fort22>` documentation for details.
* **gfdl+** - Requires two or more meteorological forcing files from the GFDL
   model. See :ref:`fort.22 <fort22>` documentation for details.
* **nwsf22** - Requires analyst-supplied ``NWS_20_fort.22``
   meteorological forcing file. See :ref:`nws-20` for information on the use of ``aswip`` program to generate the ``NWS_20_fort.22`` file.
* **pp** - Requires preprocessing of the ATCF formatted track file using the
   ASymmetric Wind Input Preprocessor (``aswip``) program (distributed with ADCIRC; see :ref:`nws-20`)
   prior to use as input for this parametric vortex model.
* **f23** - Requires analyst-supplied fort.23 wave radiation stress input file.
* **A+S** - Requires ADCIRC+SWAN executable and input files for the coupled
   version of SWAN.
* **STW** - Requires ESMF-coupled ADCIRC and STWAVE and input files for STWAVE.
* **f25** - Requires analyst-supplied :ref:`fort.25 <fort25>` ice coverage input
   file as well as a :ref:`fort.225 <fort25>` basin scale ice coverage file and
   possibly an optional :ref:`fort.227 <fort25>` region scale ice coverage file.

.. _supp_forcing_notes:

Notes
-----

.. raw:: html

   <references group="note" />

.. [1]
   ``NWS`` formats 1, 2, and 7 do not support the use of ice coverage (fort.25)
   files.

.. [2]
   Radiation stress time increment (``RSTIMINC``) represents the time increment
   between the datasets in the fort.23 file (in seconds).

.. [3]
   ``RSTIMINC`` represents the span of ADCIRC simulation time that passes
   between calls to the coupled wave model.
