# ADCIRC Build Guide

#### Authors
  * Jason Fleming, <jason.fleming@seahorsecoastal.com>
  * Zach Cobell, <zachary.cobell@arcadis.com>

## Introduction

This document describes how to build ADCIRC and ADCIRC+SWAN using both available build systems currently contained within the source distribution.

There are two build systems available:
1. CMake
2. Traditional GNU Make

The GNU Make suite of tools has been used with ADCIRC for many years and those who have built ADCIRC from its source code in the past will already be familiar with it.

CMake is a new build option and is geared at making building the ADCIRC source code from scratch mode user friendly across many platforms. While the GNU Make approach might work well for seasoned veterns, those unfamiliar with it might struggle, making CMake an attractive option. CMake can also be used in combination with the MinGW suite of tools to build ADCIRC in a Windows environment.

## Building with CMake

CMake is available to use on Windows, Linux, Macintosh, and Cygwin. It provides availability to configure the model using either a GUI interface or the command line.

### Windows Requirements
If you are building on Windows, you will need to install the MinGW64 suite of tools to provide the compilers and make command, A form of Windows Perl, and the CMake interface itself. The following links can provide the packages you need:


| Package                |   Link                          |
|------------------------|---------------------------------|
| CMake                  | https://cmake.org/download/     |
| Strawberry Perl (Perl) | http://strawberryperl.com/      |
| MinGW64                | https://mingw-w64.org/doku.php  |

Make sure that the following commands are available at your command prompt window:

```
perl
gfortran
gcc
mingw32-make
```

If they are not, consult the documentation for the individual packages to insure you've installed them correctly. 

Note that StrawberryPerl installs a number of compilers with it. You should remove them from your Windows Path variable before running CMake. This can cause errors in the build process by having incorrect compilers selected.

### Building Using the Graphical Interface

#### Unix
For Unix, Mac, and Cygwin, you will need to check that the `ccmake` program is installed on your system. From the main source directory, begin by creating a build directory.

```
mkdir build
```

Then, move into the build directory and run the CMake GUI.

```t
cd build
ccmake ..
```

CMake will start. Run the configure section by pressing the `c` key. CMake will then scan your system and determine the build environment. Once it has finished, you will have the opportunity to configure the system to your specifications. Take this opportunity to make sure the selected compilers are the ones you intend to use. CMake will default to the GNU compilers (gcc, gfortran, etc). If you edit the compilers, the configure script will need to re-run and other options besides the compilers will be reset.

Once you've selected the compilers, you can edit the other options. This will include which executables you want to build. Several options labeled `BUILD_[NAME]` will be shown. These can be toggled with your enter key.

There are also several options you can toggle that will enable different portions of the code, such as Powell Wind Drag, SWAN modified friction, and some debugging options. If you do not see some of the specific options you're looking for, type the `t` key to toggle advanced options. For some options, you will need to press the `c` key so that additional options related to the newly enabled option are shown. For example, if netCDF output is turned on, pressing `c` will now display an option to enter the netCDF library path.

If you plan to install the code on your system, be sure to set the `CMAKE_INSTALL_PREFIX` variable appropriately. If you do not have root access on your system, the default option will likely not work and you will need to specify a location where you have write access.

Once you change any options, you will need to re-run the configure portion of the code by pressing the `c` key. When the code has been fully configured, the option to generate makefiles will be enabled. Pressing `g` will write the makefiles and exit CMake. Once the makefiles have been generated, you can build the code with the selected options by running make.

To enter specific optimizations for your system, edit the `CMAKE_Fortran_FLAGS` variable. The specific optimizations required for your system should be available via your system administrator or system documentation.


```
make
```

To clean the build directory, you have two options. You can delete only the files that have been compiled by running:

```
make clean
```

or you can delete all the compiled files including your configuration settings specified in CMake by removing the `build` directory completely.

```
cd ..
rm -r build
```

#### Windows
Begin by running the CMake-GUI application from the start menu. This will bring up the interface. Click the `Browse Source` button and select root directory for the ADCIRC code package. This will contain the file `CMakeLists.txt`. Next, click the `Browse Build` button and select where you'd like to build the executables. This is not the install location; this will just be the directory where files are built. It is recommended that you do this in a folder called `build` in the root directory. It is strongly recommended that you avoid setting the build path to the source root directory. Also, be sure to set the `CMAKE_INSTALL_PREFIX` variable to an appropriate location where you want the executables to be installed.

### Building Using the Command Line

The command line interface for CMake requires that you define a few option flags. These option flags mirror exactly what is found in the GUI interface except they will all be prefixed with `-D`. For example, a build might look like:

```
mkdir build
cd build
cmake .. -DBUILD_ADCIRC=ON -DBUILD_PADCIRC=ON -DBUILD_ADCSWAN=ON \
         -DBUILD_PADCSWAN=ON -DBUILD_ADCPREP=ON -DBUILD_UTILITIES=ON \
         -DBUILD_ASWIP=ON -DBUILD_SWAN=ON -DBUILD_PUNSWAN=ON \
         -DENABLE_OUTPUT_NETCDF=ON -DENABLE_OUTPUT_XDMF=ON \
         -DNETCDFHOME=/usr -DXDMFHOME=/usr -DBUILD_UTILITIES=ON \
         -DCMAKE_Fortran_FLAGS="-mtune=native"
make
```

This process is valid for both Windows and Linux based systems. CMake is also able to build both the serial SWAN and parallel unstructured SWAN using the same build process. There is no need to use SWAN's build system.

Note that CMake will attempt to locate the environment variables "NETCDFHOME" and "XDMFHOME". If these are set in your environment, the netCDF and XDMF package paths will be automatically set, meaning there is no need to set -DXDMFHOME and -DNETCDFHOME.

## Building with GNU make

This section contains details for building ADCIRC executables and
describes the directories and files related to the build process.
For a more complete description of the files and directories that
are included with the ADCIRC source code, see the official ADCIRC
documentation on adcirc.org:

 http://adcirc.org/home/documentation/generic-adcirc-compile-time-operations/

The ADCIRC source tree contains a directory called 'work' where the
ADCIRC and ADCIRC+SWAN executables are built. This directory contains
the following files:

### Definitions
1. config.guess
  * A shell script that is called automatically by make during the build process to guess the hardware and operating system that are being used in the build process. You can execute this shell script yourself on the command line to see what it guesses and confirm that it is guessing correctly. The resulting platform guess is then used automatically in the `cmplrflags.mk` file.

2. cmplrflags.mk
  * Contains sets of compiler flags that are appropriate for various platforms. The build process uses the output from 'config.guess' to select the appropriate set of compiler flags. In the past, in the days of many varieties of proprietary unix, the hardware type and operating system type were sufficient to determine which compiler flags to select. Now that most installations of ADCIRC are on linux, it is often necessary to provide additional hints about which set of compiler flags to use, and these additional hints are provided on the command line when make is executed, as described below.

3. makefile
  * The makefile contains instructions for building `adcirc`, `padcirc`, `padcswan` (if the SWAN source code is also present), `adcprep`, `adcpost`, `hstime`, and `aswip` on Linux or Unix.

When building ADCIRC for the first time on a new platform, begin by going to the ADCIRC 'work' directory

```
cd work
```

and running `config.guess` manually to see what it guesses about the platform architecture.

```
./config.guess
```

Then, open the `cmplrflags.mk` file, to see if the relevant platform and/or compiler is defined in the file, and if it is, how it is triggered on the make command line with the `compiler=` notation.

It is also important to test the availability of the compilers that will be used by executing the `which` command. For example, if using the Intel compilers, type `which ifort`. If running ADCIRC in parallel, also try `which icc` and `which mpif90` to make sure they're all available. If they're not, you should contact your local system administrator.

Looking in the `cmplrflags.mk`, you might find, for example, a set of Intel-appropriate compiler flags by setting `compiler=intel` on the `make` command line. As a result, my first tentative step in building all the ADCIRC executables would be

```
make adcirc compiler=intel
```

Assuming that this completes successfully, try:

```
make padcirc compiler=intel
```

There are a number of additional executables that can be built including:
|Executable|Description|
|----------|-----------|
|adcprep | The parallel ADCIRC simulation preprocessor. Note that if you are building for simulations that will include SWAN, you should add SWAN=enable to your command line command to compile the code |
| aswip | The Asymmetric Wind Input Preprocessor. This utility program is required for the use of an asymmetric vortex model for tropical cyclones in ADCIRC |
| adccmp | Program to compare ADCIRC outputs to a known solution |
| inflate | Program to inflate sparse ASCII ADCIRC output files (NOUT=4) to full format ASCII ADCIRC output files (NOUT=1) |
| hstime | Program to determine the output time for a binary hot start file |

### Building the Tightly Coupled ADCIRC+SWAN
SWAN has a separate build system, since it is normally a separate code. But the SWAN build system must solve the same problems as the ADCIRC build system, that is, detecting the nature of the underlying platform and suggesting an appropriate set of compiler flags. Note that to build SWAN, you will need to have the `perl` command installed on your system.

```
cd ../swan
```

The SWAN build system has a perl script called `platform.pl` that performs much the same function as `config.guess` does for ADCIRC. The difference is that `plaform.pl` picks the compiler flags and writes them to a file called `macros.inc` for use during the build process.

Several sets of `macros.inc` files are already available in the ADCIRC source code for use with certain platforms; type

```
ls macros*
```

to see the currently available choices. To write a set of compiler flags for your platform, type:

```
perl ./platform.pl
```

You may want to have a look at the resulting `macros.inc` file to be sure that it represents the compiler you're using. If it looks ok, then you are ready to compile the coupled ADCIRC+SWAN code.

```
cd ../work
make padcswan compiler=intel
```

Once `padcswan` has been built successfully, you can build `adcprep`, the program that prepares your input files for a parallel run. If you also plan to run ADCIRC+SWAN in parallel, then add the `SWAN=enable` string to the make command line, as follows:

```
make adcprep compiler=intel SWAN=enable
```

You'll see the invocation of your C compiler during the build process for `adcprep`, since `adcprep` uses the metis library (which is written in C) to perform domain decomposition.

Once you've successfully built `adcprep`, the main set of executables is complete.

However, there are a couple utility programs you may like to build as well. The first of these is a very small program called 'hstime' whose only purpose is to read the time in seconds from an ADCIRC hotstart file and write it out to the console window. Use a command like 'make hstime compiler=intel' to build this utility.

### Compiling with netCDF Support

All the above instructions will build the ADCIRC executables without support for netCDF input or output. This was done intentionally, to simplify the build process and avoid complications resulting from the vagaries of building with netCDF.

With that said, the use of netCDF is recommended because of the many advantages of netCDF files. ADCIRC supports the use of three different netCDF capabilities: netCDF3, netCDF4 (classic model), and netCDF4 with internal compression. Your netCDF installation will need to include the base netCDF-C installation as well as netCDF-Fortran. Both are available [here]( http://www.unidata.ucar.edu/downloads/netcdf/index.jsp).

If you've already built all the executables without netCDF, its necessary to get rid of them so that you can build with netCDF.

```
make clobber
```

To enable compilation with netCDF3, the command line listed at the beginning of these instructions for building ADCIRC would be modified as follows:

```
make adcirc compiler=intel NETCDF=enable NETCDFHOME=[PATH/To/NETCDF]
```

with the `NETCDFHOME` variable indicating the directory where netCDF was installed. The ADCIRC build system expects this directory to contain `include/netcdf.inc`.

If netCDF4 without internal compression is installed on your platform, the above would change to the following:

```
make adcirc compiler=intel NETCDF=enable NETCDFHOME=/usr NETCDF4=enable
```

in which case the directory specified in the `NETCDFHOME` variable would also be expected to contain the `include/netcdf.mod` Fortran module file.

Finally, if your netCDF4 installation supports internal compression (only introduced in netCDF4 version 4.1), and you want ADCIRC to automatically make use of internal compression, change the above command line as follows:

```
make adcirc compiler=intel NETCDF=enable NETCDFHOME=/usr NETCDF4=enable \
            NETCDF4_COMPRESSION=enable
```

These changes in the make command line must be made for the build of each of the ADCIRC executables so that they are all consistent.

### Compiling with XDMF Support

Like netCDF, XDMF is another file format option available for use with ADCIRC. Certain postprocessing packages, particularly Paraview, can read the XDMF format and generate visualizations. XDMF requires that your system also supports netCDF4.

To begin an XDMF enabled build, type:

```
make adcirc NETCDF=enable NETCDF4=enable NETCDF4_COMPRESSION=enable \
            XDMF=enable NETCDFHOME=[/path/to/netcdf] XDMFHOME=[/path/to/xdmf]
```

You can get the XDMF source code [here](https://gitlab.kitware.com/xdmf).
