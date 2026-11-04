.. meta::
   :description: Fort.22x.grb2 file in ADCIRC
   :keywords: adcirc, fort.22x.grb2 file

.. _fort22x_grb2:

Fort.22x.grb2: Meteorological Inputs in GRIB2 Format
====================================================

Fort.22x.grb2 files can be used with ADCIRC version 55 and later.

The fort.22x.grb2 files are meteorological inputs in
`GRIB2 <https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/>`_ binary.
These GRIB2 files are read by ADCIRC through the
`wgrib2api <https://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/wgrib2api.html>`_
static library. The API automatically generates an inventory of the GRIB2 file
at the beginning of simulation to determine variable and valid datetime
information in the files. The user does not need to do anything other than
provide the correct filenames with valid GRIB2 variable names, in addition to
specifying the :ref:`fort.15 file <fort15>` attributes correctly;
:ref:`NWS <nws_parameter>` = 14, correct :ref:`WTIMINC <wtiminc_supplemental>`, and correct
``NCDATE``.

fort.221.grb2
-------------

Contains the surface pressure/sea level pressure data. Valid GRIB2 variable
names for this file are:

-  ``:PRMSL:mean sea level:`` *or*
-  ``:MSLET:mean sea level:`` *or*
-  ``:PRES:surface:``

The mean sea level pressure, MSLET, is known to be a better choice than its
counterpart, PRMSL, in the
`GFS <https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-forcast-system-gfs>`_
model `1 <https://luckgrib.com/tutorials/2018/08/28/gfs-prmsl-vs-mslet.html>`_.
In comparison, the
`CFSv2 <https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/climate-forecast-system-version2-cfsv2>`_
model does not provide MSLET, thus PRMSL should be used. PRES:surface (surface
pressure) should be equivalent to PRMSL/MSLET over the ocean and in general can
be used. However, problems could arise at the land/ocean interface when there is
a significant discrepancy between the ocean/land mask in ADCIRC and the
meteorological model due to resolution differences.

fort.222.grb2
-------------

Contains the wind velocity at 10-m height data. Valid GRIB2 variable names for
this file are:

-  :UGRD:10 m above ground: *and*
-  :VGRD:10 m above ground:

.. _fort.225.grb2_optional:

fort.225.grb2 [optional]
------------------------

Contains the surface ice concentration as an area fraction. Valid GRIB2 variable
names for this file are:

-  :ICEC:surface:

.. _grib2_utility:

GRIB2 utility
-------------

The `WGRIB2 <https://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/>`__ utility
can be used to manipulate and generate GRIB2 files with ease through scripts
such as bash and perl.
