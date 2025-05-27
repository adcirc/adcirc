.. meta::
   :description: Compiling in ADCIRC
   :keywords: adcirc, compiling

Compiling ADCIRC
================

Compiling is the act of turning the ADCIRC source code into an executable that
can be run. ADCIRC is written in Fortran and C, so Fortran and C compilers are
needed, as a minimum.

Dependencies
------------

The ADCIRC model is minimally dependent on external libraries, though they may
be included at the user's option and enable many important features of the model
code.

.. list-table:: List of dependencies
   :header-rows: 1
   :widths: 20 20 20 20 60
   :class: wrap-table tighter-table

   * - Component
     - Availability
     - Required
     - Source
     - Description
   * - Fortran
     - All versions
     - Required
     - --
     - ADCIRC is mostly written in Fortran, which excels at large matrix operations and numeric operations. Various Fortran compilers have been tested with ADCIRC including the GNU, Intel, and PGI Fortran compilers. The Fortran compiler should support at least Fortran90 standards, though Fortran2008 is recommended.
   * - C
     - All versions
     - Required
     - --
     - ADCIRC uses C code to access both the Metis library and some low-level system functions. As with Fortran, various compilers have been tested. The C compiler should support at least C98 standards.
   * - GNU Make
     - All versions
     - Required
     - --
     - GNU Make is a program used to build source code. It can be used by both CMake and the legacy Makefile, though CMake will allow alternate build programs to be used, such as Ninja. It is natively available on nearly all Linux based systems.
   * - CMake
     - v53.00
     - Required to use CMake build system
     - `Source <https://cmake.org/download/>`__
     - CMake is an open-source, cross-platform family of tools designed to build, test and package software. CMake is used to control the software compilation process using simple platform and compiler independent configuration files, and generate native makefiles and workspaces that can be used in the compiler environment of your choice. The suite of CMake tools were created by Kitware in response to the need for a powerful, cross-platform build environment for open-source projects such as ITK and VTK.
   * - MPI
     - All versions
     - Required for all parallel codes
     - `MPICH <https://www.mpich.org/>`__, `OpenMPI <https://www.open-mpi.org/>`__, `MVAPICH2 <https://mvapich.cse.ohio-state.edu/>`__
     - Message Passing Interface (MPI) is a standardized and portable message-passing standard designed by a group of researchers from academia and industry to function on a wide variety of parallel computing architectures. The standard defines the syntax and semantics of a core of library routines useful to a wide range of users writing portable message-passing programs in C, C++, and Fortran. Different MPI libraries will provide better performance or usability on different systems. As a general rule, OpenMPI is the easiest implementation for beginners to set up on their own and is provided as a package in many common Linux distributions. MVAPICH2 has many performance enhancements that make it ideal for running at very large scales over InfiniBand.
   * - Perl
     - v49.00
     - Required for building the included SWAN model
     - --
     - Perl is a text manipulation language. The SWAN model uses Perl to perform source code preprocessing.
   * - netCDF
     - v51.00
     - Required to write netCDF format outputs, for example, NOUTGE=3,5
     - `Source <https://www.unidata.ucar.edu/software/netcdf/>`__
     - NetCDF (Network Common Data Form) is a set of software libraries and self-describing, machine-independent data formats that support the creation, access, and sharing of array-oriented scientific data. ADCIRC will require both the C and Fortran APIs are installed on the target system. netCDF4, which uses HDF5 as the underlying library for data access and compression, is highly recommended. This will allow the user to write compressed output formats which will save disk space.
   * - XDMF
     - v52.00
     - Required to write XDMF formatted outputs, for example, NOUTGE=7
     - `Source <http://www.xdmf.org/index.php/Main_Page>`__
     - XDMF (eXtensible Data Model and Format) is a library providing a standard way to access data produced by HPC codes. Data format refers to the raw data to be manipulated, the description of the data is separate from the values themselves. It distinguishes the metadata (Light data) and the values themselves (Heavy data). Light data is stored using XML, Heavy data is typically stored using HDF5, so some information is stored redundantly in both XML and HDF5. A Python interface exists for manipulating both Light and Heavy data. ParaView, VisIt and EnSight visualization programs are able to read XDMF.
   * - Grib2
     - v55.00
     - Required for NWS=14 Grib file processing
     - `Source <https://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/>`__
     - GRIB (GRIdded Binary or General Regularly-distributed Information in Binary form[1]) is a concise data format commonly used in meteorology to store historical and forecast weather data. Beginning with ADCIRC v55.00, this data can be directly read and used as forcing for ADCIRC. Using the CMake build system, described below, will automatically build this dependency if selected.
   * - DateTime
     - v55.00
     - Required for NWS=14 Grib file processing
     - `Source <https://github.com/wavebitscientific/datetime-fortran>`__
     - The DateTime library is a Fortran library that provides date manipulation and arithmetic operations. Using the CMake build system, described below, will automatically build this dependency if selected.

.. _cmake_build_system:

CMake Build System
------------------

The CMake build system provides users with a cross-platform (Linux, MacOSX,
Windows) interface for building the ADCIRC model. CMake can be used as either a
command-line utility, a terminal GUI with the Curses library, or a Qt GUI on
Windows. CMake is available on most Linux distributions as either a package or
in the case of shared clusters, as a module. The source code and compiled
executables are available `here <https://cmake.org/download/>`__.

.. _cmake_environment_variables:

CMake Environment Variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ADCIRC CMake build system will search for specific environment variables to
help it find packages in nonstandard locations, as is often the case on shared
systems.

.. list-table::
    :widths: 15 60
    :class: wrap-table tight-table

    * - Variable
      - Description
    * - NETCDFHOME
      - Sets the home path for netCDF C and Fortran libraries. This assumes that the libraries are installed to the same location. This folder should contain the folders "lib" and "include"
    * - XDMFHOME
      - Sets the home path for the XDMF libraries. This folder should contain the folders "lib" and "include"

.. _general_usage:

General Usage
~~~~~~~~~~~~~

CMake is generally used to create a shadow build, that is the build files are
wholely separated from the source code. The simplest way to do this is to create
a folder within the source code folder called something like "build". For the
instructions below, we'll assume that you are inside this directory during the
command examples.

.. _parallel_builds:

Parallel Builds
~~~~~~~~~~~~~~~

CMake can build the ADCIRC model source code using multiple processors. Add the
flag "-j#" where "#" is the number of processors to use to your make command.
This will allow each component of the code to be built separately and may not
provide benefits above 6-7 processors. This cannot be done using the legacy
makefile as it does not contain the appropriate logic.

.. code-block:: bash

   make -j4

.. _cmake_command_line_usage:

CMake Command-Line Usage
~~~~~~~~~~~~~~~~~~~~~~~~

The command-line interface allows the CMake build system to run with minimal
user interaction and can be scripted. CMake building happens in two steps.
First, there is a configuration step where the environment is sensed, tested,
and a Makefile is generated. Second, a build step is triggered which executes
the generated build files. Variables are specified for CMake using the
"-DVARIABLE=VALUE" syntax.

By default, CMake will always select the system default compilers. This can
easily be overridden by specifying the desired compilers on the command line.
For example, to use the Intel compilers:

.. code-block:: bash

   ccmake .. -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort

Note that when specifying the compilers, the MPI library will need to match the
specified compilers (i.e. the mpif90 found in your PATH variable should
correspond to ifort) or parallel code will not be enabled.

When the above command is run, CMake will perform a series of checks, however,
it will not know which programs in the ADCIRC source code suite you want to
build. The following table shows the build options that can be enabled.

.. list-table::
    :widths: 3 2 4
    :class: wrap-table tight-table

    * - Variable
      - Conditional
      - Description
    * - BUILD_ADCIRC
      - --
      - Builds the serial ADCIRC executable
    * - BUILD_ADCPREP
      - adcprep
      - Builds the ADCIRC parallel preprocessor
    * - BUILD_ADCSWAN
      - SWAN code is enabled
      - Builds the serial coupled ADCIRC+SWAN executable
    * - BUILD_ASWIP
      - --
      - Builds the ADCIRC asymmetric wind input preprocessor
    * - BUILD_LIBADCIRC_SHARED
      - --
      - Builds ADCIRC as a shared object library
    * - BUILD_LIBADCIRC_STATIC
      - --
      - Builds ADCIRC as a static library
    * - BUILD_PADCIRC
      - MPI code is enabled
      - Builds the parallel ADCIRC executable
    * - BUILD_PADCSWAN
      - MPI code is enabled, SWAN code is enabled
      - Builds the parallel coupled ADCIRC+SWAN executable
    * - BUILD_UTILITIES
      - --
      - Builds the utility codes in the ADCIRC source code
    * - BUILD_ADCRUN
      - --
      - Builds the ADCIRC executable
    * - BUILD_ADCRUN
      - --
      - Builds the ADCIRC executable
    * - BUILD_ADCRUN
      - --
      - Builds the ADCIRC executable
    * - BUILD_ADCRUN
      - --
      - Builds the ADCIRC executable
    * - BUILD_ADCRUN
      - --
      - Builds the ADCIRC executable
    * - BUILD_ADCRUN
      - --
      - Builds the ADCIRC executable


An example to build adcprep, adcirc, and padcirc might look like:

.. code-block:: bash

   ccmake .. -DBUILD_ADCIRC=ON -DBUILD_PADCIRC=ON -DBUILD_ADCPREP=ON

When CMake runs successfully, you can then simply type "make" to build the
executables.

.. _other_influential_variables:

Other Influential Variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following shows additional variables that can be specified in the command
line which can enable features or optimizations. They are used with the same
syntax as above.

.. list-table::
    :widths: 25 22 30
    :class: wrap-table tighter-table

    * - Variable
      - Description
      - Sample
    * - ENABLE_OUTPUT_NETCDF
      - Enables netCDF formatted output
      - -DENABLE_OUTPUT_NETCDF=ON
    * - ENABLE_OUTPUT_XDMF
      - Enables XDMF formatted output
      - -DENABLE_OUTPUT_XDMF=ON
    * - CMAK
      - Optimization flags
      - -DCMAKE_Fortran_FLAGS_RELEASE="-O3 -funroll-loops --param max-unroll-times=4 -march=native"

.. _cmake_curses_gui:

CMake Curses GUI
~~~~~~~~~~~~~~~~

The Curses GUI available with cmake allows the user to explore and set variables
in a more user-friendly way. It is still recommended to set the compilers on the
command line during the first call of cmake. The executable to use the GUI is
"ccmake" rather than "cmake".

.. code-block:: bash

   ccmake .. -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort

The user will be presented with an interface that looks like:

.. figure:: /_static/images/user_guide/compiling/Adcirc_cmake.JPG.jpg
   :alt: Adcirc_cmake.JPG


The options on the screen can be toggled through and edited interactively using
the keyboard. When you are happy with the specified options, press "c" to
configure. If there are no errors, you will be given the option to press "g" to
generate makefiles. After generating, ccmake will close and you can type "make".

.. figure:: /_static/images/user_guide/compiling/Adcirc_cmake_build.JPG.jpg
   :alt: Adcirc_cmake_build.JPG
