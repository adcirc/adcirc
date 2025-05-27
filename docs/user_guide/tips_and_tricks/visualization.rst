.. meta::
   :description: Visualization in ADCIRC
   :keywords: adcirc, visualization

Visualization
=============

For more general tools and utilities, see :ref:`Tools <tools>`. In the course of
running ADCIRC you may wish to visualize or plot data pertaining to your model.
Visualizing data can help you fine tune the mesh design, conduct QA/QC on
results, and communicate key facets of the modeling to others. Examples of
datasets you might plot include mesh elevation; :ref:`nodal
attributes <nodal_attributes>` such as :ref:`Manning’s
n <mannings_n_at_sea_floor>`; or modeled results such as water surface
elevation. For visualizing on a 2D flexible mesh, ADCIRC users often choose to
plot them as color contours that follow the mesh triangulation.

.. _plotting_considerations:

Plotting Considerations
-----------------------

Options abound for generating plots, and the best tool for the job depends on
project goals and available resources/software. Some factors to consider in
choosing a visualization tool include:

-  Interactive vs. static image – do you intend to zoom and pan around your
   visualized data and/or change contour intervals on the fly, or would you
   prefer to batch process image files at fixed zoom extents and contour limits?
-  Data transfer – do you intend to run simulations on the same computer you
   visualize results?
-  Desired processing speed – does the speed of plotting matter to you? This can
   be important, e.g. when making videos.
-  NetCDF vs. ASCII data – while ADCIRC mesh and nodal attributes input files
   must take ASCII form, outputs can be ASCII or NetCDF. Consider this since
   some plotting tools only work with one or the other. A few considerations:

   -  Multiple applications can read a NetCDF4 file simultaneously, which could
      save considerable time if you can run your plotting tool concurrently on
      multiple processors within one compute node.
   -  NetCDF ADCIRC output files self-contain mesh node, connectivity, and
      boundary information.
   -  The learning curve to coding with netCDF files can be much steeper than
      ASCII.

-  Projection – ADCIRC is often run in geographic coordinate system, which
   references a spherical object; do you need the plotting tool to accurately
   project geographic coordinates, or is distortion acceptable?
-  Which files do you need to visualize? Do you need to plot nodal attribute
   (fort.13) data, or just ADCIRC output?
-  Effort – do you prefer a DIY, high customization approach or a ready-to-run
   program? Somewhat related question: do you want to buy software or leverage
   open source libraries? Common items:

   -  Overlaying a vector field showing wind magnitude and direction or wave
      direction
   -  Adding background aerial imagery
   -  Auto-calculating contour bounds within the specified zoom extents
   -  Overlaying shapefiles showing county boundaries, inshore waterways, storm
      track, etc.

.. _sms_surface_water_modeling_system:

SMS (Surface-water Modeling System)
-----------------------------------

Aquaveo’s SMS software features a suite of many GUI-based pre- and
post-processing tools that allow for interactive visualization. SMS requires a
license, but out of the box you can easily open and view an ADCIRC mesh and any
associated data, including nodal attributes files. SMS plots mesh data as a
triangulated irregular network (TIN) surface, allowing the user to change
contour intervals as desired. You will need to download files from the HPC you
run ADCIRC onto your desktop to plot them with SMS. A few notes about SMS:

-  ADCIRC output NetCDF support is “coming soon” as of version 12.1. SMS
   supports other types of NetCDF data but cannot open a \*.63.nc file.
-  When ADCIRC output files are opened, SMS automatically converts them from
   ASCII to binary HDF5 format. This format allows data to load quickly –
   anecdotal reports estimate this speed to be faster than the same data
   displayed as a TIN in ArcMap.
-  SMS can...

   -  project data in geographic coordinates properly.
   -  download zoom-specific aerial imagery to display as a background layer.
   -  supports vector overlays for fort.74-style ADCIRC data (u v format only,
      so you need to reformat swan_DIR.63 before displaying as vector field).
   -  display shapefiles, but with limited options for changing their symbology.

FigureGen
---------

`FigureGen <https://ccht.ccee.ncsu.edu/figuregen-v-49/>`__ (by Casey Dietrich of
NC State) batch produces images of model data. FigureGen is a Fortran program
that uses the Generic Mapping Tools (GMT) library. To use FigureGen, you must
first install GMT and any of its auxiliary programs that don’t already exist on
your system. These may already be installed if on a supercomputer. When you
compile FigureGen, you can use flags to enable Google Earth output, NetCDF
support, and parallel processing. You’ll then edit a text file of input
parameters, including zoom extent and required output snaps (for time-varying
ADCIRC output file), and FigureGen will generate the series of requested images.
A few notes about FigureGen:

-  FigureGen is a triple threat – it’s fast, easy to use, and free (but please
   cite as appropriate!)
-  You can typically run FigureGen on the same machine you run ADCIRC, no need
   to download large output files.
-  GMT supports map projections so FigureGen plots geographic coordinates
   correctly.
-  Features options for overlaying a vector field
   (`fort.74 <fort.74_file>`__-style data) including automatically finding the
   maximum vector magnitude present within the zoom extent.

   -  Can automatically find the maximum data value present within the zoom
      extent to set the contour maximum.
   -  Can overlay a GSHHS shoreline. If your mesh includes detailed inshore
      waterways, this shoreline may not include sufficient resolution.
   -  Can add a time bar.
   -  Can plot the difference between two (scalar?) files.
   -  Can plot nodal attribute data (except for *elemental_slope_limiter*, but
      that is easy to add yourself)
   -  Can produce several types of image files, and the user can set resolution
      to control file size.

Kalpana
-------

`Kalpana <https://ccht.ccee.ncsu.edu/kalpana/>`__ (by Rosemary Cyriac) produces
either a shapefile or a Google Earth kml file of ADCIRC netCDF output or mesh
bathymetry. Kalpana is a Python script that generates color contours using the
mesh triangulation via the 2D plotting library Matplotlib, then converts the
contours to either shapefile via the Fiona/Shapely libraries, or to kmz via the
Simplekml library. Both the shapefile option and kml option produce polygons
that enclose areas with data values that fall into one contour bin.

-  How-to instructions can be found
   `here <https://ccht.ccee.ncsu.edu/how-to-run-kalpana/>`__
-  Kalpana is free (but cite as appropriate!)
-  No support for plotting nodal attributes parameters yet.
-  Shapely/Fiona and Simplekml support map projections and thus properly plot
   geographic coordinates.
-  Like many Python scripts, Kalpana should be run in a virtual environment.
   Setting up this virtual environment can be simple or frustrating depending on
   existing system installations on your HPC. Jason Fleming wrote a `great post
   on his
   experience <https://ccht.ccee.ncsu.edu/installing-python-modules-for-kalpana/>`__
   in installing the various Python libraries within a virtual environment.
-  Kalpana runs on an HPC, so you don’t need to download data to your desktop to
   plot.
-  There exists a file size and speed trade-off to the two file format options –
   the kml option is slower to generate, but the file size is relatively small;
   the shapefile option generates very quickly (nearly as fast as FigureGen) but
   the file size is on the order of the size of the file it’s plotting. Also
   note that the speed of the kml generation depends on the number of lat/long
   bins you set up in the command line, and the number of bins depends on mesh
   extent and resolution.
-  The shapefiles and kml files produced allow the user to zoom and pan around,
   but not to dynamically change contour intervals on the fly – the polygons
   represent an instance or a snapshot of the data.
-  Both shapefile and kml are popular file formats with the outside world, so
   they’re good choices for communicating your model results to others.
-  Neither file format supports vector fields. If you want vector fields, you
   should edit the Kalpana script to skip the Shapely/Fiona or Simplekml parts
   and use Matplotlib directly to add the vector field and save the Matplotlib
   plot as a png.

ArcGIS
------

Esri’s classic software contains no canned tools to read or write ADCIRC data.
The user must develop scripts to read ADCIRC files via ArcMap's integrated
Python package (ArcPy). Once you’ve imported ADCIRC mesh geometry, nodal
attributes, and model results data into a geodatabase, you can easily create a
TIN using the mesh triangulation as hard breaklines to view your data in a
similar fashion to SMS. Note that the license for ArcMap costs significantly
more than an SMS license. So why the heck would you ever consider using ArcStuff
for your visualization needs?

-  Esri is everywhere – you might already have access to a network license or to
   colleagues with advanced GIS knowledge.
-  While ArcMap doesn’t feature canned tools to read ADCIRC data, it does
   feature canned geoprocessing tools that can aid greatly in assigning nodal
   attribute data or interpolating mesh node elevations. For example, you can
   map data showing land use/land cover in raster or polygon format onto your
   nodes to assign Manning’s *n*, or you can create custom schemes for
   interpolating DEM data onto your mesh nodes.
-  As of version 10.3, ArcPy includes the netCDF4 library.
-  ArcMap allows for a very high level of customization in your image layout.
   North arrows, scales, legends, infinite symbology and overlays – wow!
-  But, you must download ADCIRC output data from your HPC to wherever your
   ArcGIS license exists.

Python
------

See :doc:`Tools <../../tools/index>`.

Matlab
------

See :doc:`Tools <../../tools/index>`.

Paraview
--------

Super powerful, but can be a pain. Tough to generalize here, email the listserv
if interested.

External Links
--------------

-  `SMS by
   Aquaveo <https://www.aquaveo.com/software/sms-surface-water-modeling-system-introduction>`__
-  `FigureGen v49 <https://ccht.ccee.ncsu.edu/figuregen-v-49/>`__ NC State
   Coastal & Computational Hydraulics Team
-  `Kalpana <https://ccht.ccee.ncsu.edu/kalpana/>`__ NC State Coastal &
   Computational Hydraulics Team
-  `ADCIRC.io <http://www.adcirc.io/>`__
-  `ArcMap <http://desktop.arcgis.com/en/arcmap/>`__
-  `Matlab <https://www.mathworks.com/products/matlab.html>`__
-  `Paraview <https://www.paraview.org/>`__
-  `Generic Mapping Tools (GMT) <https://github.com/GenericMappingTools/gmt>`__
-  `Matplotlib <https://matplotlib.org/>`__


