Introduction to ADCIRC
======================

What is ADCIRC?
---------------

ADCIRC (ADvanced CIRCulation) is a highly developed computer program for solving time dependent, free surface circulation and transport problems. The model employs the finite element method in space with unstructured grids and uses the Generalized Wave-Continuity Equation (GWCE) formulation. ADCIRC can be run either as a two-dimensional depth integrated (2DDI) model or as a three-dimensional (3D) model. The system has been optimized for parallel computing architectures, making it highly efficient for large-scale simulations.

Model Features
--------------

ADCIRC can be used to simulate:

* Tidal and Wind-Driven Circulation
* Hurricane Storm Surge and Flooding
* Wave-Current Interaction
* Dredging and Material Disposal Studies
* Transport of Materials
* Baroclinic Circulation
* Coastal Inundation and Protection Studies

ADCIRC Programs
---------------

The ADCIRC system provides the following programs:

* **ADCIRC**: Serial version of the core circulation model
* **PADCIRC**: Parallel version of ADCIRC for high-performance computing environments
* **ADCPREP**: Domain decomposition utility for preparing parallel simulations
* **SWAN**: Serial unstructured wave model
* **ADCSWAN**: Serial coupled ADCIRC+SWAN model
* **PADCSWAN**: Parallel coupled ADCIRC+SWAN model
* **PUNSWAN**: Parallel unstructured SWAN model
* **ASWIP**: Preprocessor for converting ATCF formatted hurricane forcing data
* **LIBADC**: Library of ADCIRC subroutines for incorporation into other codes
* **Utilities**: Various utility codes for pre- and post-processing

System Requirements
-------------------

ADCIRC is highly scalable and can be run on a range of platforms, from a standard desktop/laptop computer to high-performance computing (HPC) clusters. The hardware requirements depend on the size and complexity of the computational grid and the length of the simulation.

ADCIRC can be compiled and run on:

* Linux/Unix systems
* Windows (using compatible compilers)
* MacOS

The model is written in Fortran and requires:

* Fortran compiler (gfortran, Intel Fortran, etc.)
* MPI library for parallel execution
* Optional: NetCDF libraries for certain output formats

Getting Started
---------------

To get started with ADCIRC, you'll need to:

1. Install prerequisites (compiler, MPI, etc.)
2. Build ADCIRC using either traditional make or CMake method
3. Prepare input files
4. Run simulations
5. Analyze output

For detailed installation and running instructions, see :doc:`Getting Started <../getting_started/index>`.

ADCIRC Files
------------

ADCIRC uses a set of input files to define the model domain, boundary conditions, and runtime parameters:

* Fort.14: Grid and boundary information
* Fort.15: Model parameters and periodic boundary conditions
* Fort.13: Nodal attributes
* Additional files for meteorological forcing, wave forcing, etc.

The model produces various output files containing water levels, currents, and other variables at specified locations and time intervals.

For detailed descriptions of all input and output files, and parameters, see the :doc:`Input Files <../input_files/index>`, :doc:`Output Files <../output_files/index>`, and :doc:`Parameters <../parameter_definitions/index>` sections of this documentation. 

Contributors
------------

For information about ADCIRC's authors and development team, please refer to the `README <https://github.com/adcirc/adcirc?tab=readme-ov-file>`_ file in the ADCIRC repository.

