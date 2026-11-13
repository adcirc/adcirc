:orphan:

SubgridADCIRCUtility
====================

SubgridADCIRCUtility is a Python-based toolkit for creating subgrid input files for subgrid-enabled ADCIRC. It enables the representation of high-resolution terrain and landcover features without requiring extremely fine mesh resolution, improving simulation accuracy while maintaining computational efficiency.

Features
--------

* **Subgrid Lookup Tables**: Generate lookup tables that store subgrid information for ADCIRC
* **DEM Integration**: Process high-resolution Digital Elevation Models for subgrid calculations
* **Landcover Processing**: Incorporate landcover data for Manning's n coefficient assignment
* **Water Level Range**: Compute subgrid corrections over user-defined water surface elevation ranges
* **Table Simplification**: Optimize lookup tables for efficiency in ADCIRC simulations
* **Multiple Data Sources**: Support for multiple DEM and landcover files in a single calculation
* **Manning's Value Customization**: Support for both default and custom Manning's n value tables
* **NetCDF Output**: Generate outputs in NetCDF format for direct use with subgrid ADCIRC

Links
-----

* `GitHub Repository <https://github.com/ccht-ncsu/subgridADCIRCUtility>`_ 