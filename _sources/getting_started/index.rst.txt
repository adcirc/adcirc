Getting Started with ADCIRC
===========================

This guide will help you get ADCIRC up and running on your system.

Prerequisites
-------------

Before building ADCIRC, ensure you have the following installed:

* Fortran compiler (gfortran, Intel Fortran)
* MPI implementation (OpenMPI, MPICH)
* NetCDF libraries (optional but recommended)
* CMake (version 3.12 or higher) if using CMake build method
* Perl (required for SWAN coupling)

Building ADCIRC
---------------

ADCIRC can be built using two different methods:

Traditional Build Method
~~~~~~~~~~~~~~~~~~~~~~~~

1. Download the source code::

    git clone https://github.com/adcirc/adcirc.git
    cd adcirc/work

3. Build the executables specigying a compiler and options::

    make adcirc padcirc adcprep compiler=intel NECDF=enable netcdf4=enable

CMake Build Method
~~~~~~~~~~~~~~~~~~

1. Download the source code and create a build directory::

    git clone https://github.com/adcirc/adcirc.git
    mkdir build
    cd build

2. Configure and build using CMake. There are several ways to configure:

   Basic configuration::

    # Configure with defaults
    cmake ..
    
    # Configure with common options
    cmake .. \
        -DBUILD_ADCIRC=ON \
        -DBUILD_PADCIRC=ON \
        -DBUILD_ADCPREP=ON \
        -DENABLE_OUTPUT_NETCDF=ON \
        -DCMAKE_BUILD_TYPE=Release \

   Interactive configuration::

    # Use cmake and ccmake for interactive option selection
    cmake ..
    ccmake ..

    # If you encounter issues with the netcdf-fortran installation, specify the path to the netcdf-fortran installation explicitly:
    ccmake .. -DNETCDF_F90_ROOT=/path/to/netcdf-fortran-install

   Available build options include:

   Executable Options:
   
   * BUILD_ADCIRC: Build serial ADCIRC executable
   * BUILD_PADCIRC: Build parallel ADCIRC executable (requires MPI)
   * BUILD_ADCPREP: Build parallel preprocessor (requires MPI)
   * BUILD_ADCSWAN: Build serial coupled SWAN+ADCIRC (requires Perl)
   * BUILD_PADCSWAN: Build parallel coupled SWAN+ADCIRC (requires MPI and Perl)
   * BUILD_SWAN: Build serial SWAN executable
   * BUILD_PUNSWAN: Build parallel unstructured SWAN
   * BUILD_UTILITIES: Build ADCIRC utility programs

   Output Format Options:
   
   * ENABLE_OUTPUT_NETCDF: Enable NetCDF output format
   * ENABLE_OUTPUT_XDMF: Enable XDMF output format (requires NetCDF)

   Debug Options:
   
   * DEBUG_FULL_STACK: Enable detailed stack trace
   * DEBUG_LOG_LEVEL: Force debug log level for screen messages
   * Various component-specific debug options (e.g., DEBUG_WIND_TRACE, DEBUG_MESH_TRACE)

   Architecture-specific Options:
   
   * Machine-specific optimizations for different platforms (IBM, SGI, SUN, CRAY)
   * VECTOR_COMPUTER: Enable vector computer optimizations

3. Build the executables::

    # Build using 4 cores
    make -j4
    
    # Build using all available cores
    make -j
    
    # Build using single core
    make

Running ADCIRC
--------------

ADCIRC requires several input files to run a simulation:

1. Prepare the required input files:

   * fort.14 (mesh file)
   * fort.15 (model parameters)
   * fort.13 (optional nodal attributes)

2. Run ADCIRC:

   For serial execution::

    ./adcirc

   For parallel execution::

    adcprep --np <number_of_processors> --partmesh   # partition the mesh
    adcprep --np <number_of_processors> --prepall    # prepare all the files
    mpirun -np <number_of_processors> ./padcirc      # run the parallel code

Example Run
~~~~~~~~~~~

Here's a basic example of running a tidal simulation.

1. Prepare the input files::

    # Clone the ADCIRC Test Suite
    git clone https://github.com/adcirc/adcirc-testsuite.git

    # Go to the directory for a quarter annular 2D test case with netcdf format output
    cd adcirc-testsuite/adcirc/adcirc_quarterannular-2d-netcdf

    # Run ADCIRC in serial
    ./adcirc

    # Run ADCIRC in parallel with 4 processors
    adcprep --np 4 --partmesh   # partition the mesh
    adcprep --np 4 --prepall    # prepare all the files
    mpirun -np 4 ./padcirc      # run the parallel code

The simulation will create several output files in the netcdf format including:

* fort.61.nc - elevation time series at specified stations
* fort.63.nc - elevation time series at all nodes
* fort.64.nc - velocity time series at all nodes

For more detailed information on the input and output files, refer to the :doc:`Input Files <../input_files/index>` and :doc:`Output Files <../output_files/index>` sections. 