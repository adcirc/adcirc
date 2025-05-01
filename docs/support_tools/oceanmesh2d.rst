:orphan:

OceanMesh2D
===========

OceanMesh2D is a two-dimensional triangular mesh generator with pre- and post-processing utilities written in pure MATLAB. It's specifically designed for building models that solve shallow-water equations or wave equations in coastal environments, including support for ADCIRC, FVCOM, WaveWatch3, SWAN, SCHISM, and Telemac models. The software requires no additional MATLAB toolboxes.

Features
--------

* **Mesh Generation**: Creates high-quality triangular meshes for coastal simulations with multi-scale resolution
* **Boundary Handling**: Advanced tools for coastline representation and boundary condition management
* **Attribute Support**: Tools for assigning and manipulating nodal attributes including Manning's coefficients
* **Bathymetry Processing**: Interpolation capabilities for bathymetry and topography data
* **DEM Integration**: Tools for working with Digital Elevation Models of various formats
* **Tidal Data Support**: Generation of tidal boundary conditions from various datasets
* **Mesh Editing**: Capabilities for remeshing, cleaning, and modifying mesh elements
* **Visualization**: Comprehensive plotting functions with support for various visualization types
* **River Forcing**: Support for creating river flow boundary conditions (fort.20 files)
* **Mesh Quality Tools**: Functions for assessing and improving mesh quality

Links
-----

* `GitHub Repository <https://github.com/CHLNDDEV/OceanMesh2D>`_
* `User Guide <https://github.com/CHLNDDEV/OceanMesh2D/tree/master/UserGuide>`_ 