:orphan:

Kalpana
=======

Kalpana is a Python module developed by the Coastal & Computational Hydraulics Team at North Carolina State University for converting ADCIRC model outputs to geospatial vector formats and downscaling maximum water elevations. It enables the visualization of ADCIRC results in common GIS formats and improves inundation representation through high-resolution downscaling.

Features
--------

* **Vector Conversion**: Convert ADCIRC outputs to geospatial vector formats (shapefile or KMZ)
* **Time-Varying Support**: Process both time-varying outputs (fort.63.nc, swan_HS.63.nc) and time-constant outputs (maxele.63.nc)
* **Downscaling**: Downscale maximum water elevations to higher-resolution raster for improved inundation representation
* **Surface Expansion**: Expand water surfaces to intersect with ground surface beyond ADCIRC-predicted extents
* **Multiple Methods**: Support for static and head-loss methods for storm surge expansion and contraction
* **Python Integration**: Fully implemented in Python with modern package structure
* **Visualization Options**: Generate polylines or polygons in multiple formats for visualization

Links
-----

* `GitHub Repository <https://github.com/ccht-ncsu/Kalpana>`_
* `Research Paper <https://link.springer.com/article/10.1007/s11069-021-04634-8>`_ 