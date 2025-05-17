# WGRIB 2

Provides functionality for interacting with, reading, writing, and
manipulating grib2 files, with a CMake build.

See wgrib2 documentation
[here](https://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/).

See also the [develop branch doxygen
build](https://noaa-emc.github.io/wgrib2).

This release of the wgrib2 package, with CMake build capability, is
part of the [NCEPLIBS](https://github.com/NOAA-EMC/NCEPLIBS) project.

# Users

Wgrib2 is used by the following projects:
* [NOMADS](nomads.ncep.noaa.gov) uses wgrib2 in grib_filter (a wrapper for
   wgrib2) wgrib2 is used to create the *.idx files for random access.
* [global-workflow](https://github.com/NOAA-EMC/global-workflow/),
  including post-processing and downstream products. Downstream jobs
  use wgrib2 to regrid GRIB2 files and extract variables. The new
  ocean/ice components will use wgrib2 to convert netCDF files to
  GRIB2 files.
* [Ensemble Verification
  System(EVS)](https://github.com/NOAA-EMC/EVS), for processing GRIB2
  files.
* [Rapid Refresh Forecast System
  (RRFS)](https://gsl.noaa.gov/focus-areas/unified_forecast_system/rrfs).
* [The Climate Prediction Center
  (CPC)](https://www.cpc.ncep.noaa.gov/) uses wgrib2 extensively to
  process GRIB2 files in their real-time applications.
* The [MET/METplus](https://dtcenter.org/community-code/metplus)
  software uses wgrib2 to extract information from GRIB2 files.
* The [Global Forecast System
  (GFS)](https://www.ncei.noaa.gov/products/weather-climate-models/global-forecast)
  post-processing of the atmosphere component has been using wgrib2 to
  interplate the UPP output master files in Gaussian grid into
  pgrb2/pgrb2b files in lat-lon grid.
* Wgrib2 is the suggested regridding tool for GRIB2 by
  [UPP](https://github.com/NOAA-EMC/UPP).
* [pywgrib2_s](https://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/pywgrib2_s.html)
  and pywgrib2_xy call wgrib2.
* [RNomads](https://cran.r-project.org/web/packages/rNOMADS/rNOMADS.pdf)
  is an R module that can get grib files from nomads.ncep.noaa.gov and
  read them. RNomads will run in windows and linux, and the code uses
  wgrib2.
* wgrib2 is included in anaconda, redhat and other repositories.
* [GrADS](http://cola.gmu.edu/grads/) uses wgrib2 to create control
  files. Also used to create ctl files for the grads dods server.

### Authors

Wesley Ebisuzaki, Reinoud Bokhorst, John Howard, Jaakko Hyvätti, Dusan
Jovic, Daniel Lee, Kristian Nilssen, Karl Pfeiffer, Pablo Romero,
Manfred Schwarb, Gregor Schee, Arlindo da Silva, Niklas Sondell, Sam
Trahan, George Trojan, Sergey Varlamov

CMake Build and testing: Kyle Gerheiser, Alyson Stahl, Edward
Hartnett, Alex Richert

# Prerequisites

Building wgrib2 requires, CMake (version 3.15+), a C compiler, a Fortran compiler
(optional), and various 3rd party libraries depending on the selected build options.

- A Fortran compiler is only required to build the Fortran API (MAKE_FTN_API=ON) or for grid interpolation (USE_IPOLATES=ON).
- wgrib2 optionally depends on netcdf-c to read and write NetCDF files (USE_NETCDF=ON)
- [NCEPLIBS-ip](https://github.com/NOAA-EMC/NCEPLIBS-ip) library (version 5.1.0+) is required for grid interpolation (USE_IPOLATES=ON).
- [libaec.a](https://gitlab.dkrz.de/k202009/libaec) (version 1.0.6+) is required to support files with AEC compression (USE_AEC=ON). 
- [NCEPLIBS-g2c](https://github.com/NOAA-EMC/NCEPLIBS-g2c) library (version 1.9.0+) is required to use wgrib2 with JPEG (USE_JASPER=ON or USE_OPENJPEG=ON) and PNG (USE_PNG=ON) compression. Make sure that the library is built with the appropriate build options enabled (some are not enabled by default). For example, to use wgrib2 with OpenJPEG, g2c must be built with OpenJPEG enabled and Jasper disabled. Please refer to [g2c's readme file](https://github.com/NOAA-EMC/NCEPLIBS-g2c?tab=readme-ov-file) for more detail.


# Installing

Download the tarball from the release page and unpack it, and cd into the main directory of the library. Then run the following commands, substituting your directory locations for the CMAKE_INSTALL_PREFIX (where wgrib2 will be installed), and the CMAKE_PREFIX_PATH (where the build will look for dependencies):

```
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=install -DCMAKE_PREFIX_PATH="/usr/local/NCEPLIBS-g2c;/usr/local/NCEPLIBS-ip"
make
make install
```

The CMake build provided here does not currently support MySQL.

# Using the wgrib2 library

wgrib2 has traditionally been used as a command-line utility, but can now be used as a library.

To use the wgrib2 library, the CMake build offers a package config.

After running `find_package(wgrib2)` CMake will generate the following
targets for use in your project:

```
wgrib2::wgrib2_lib (c library)
wgrib2::wgrib2_api (fortran target)
```

### References

[Wesley Ebisuzaki 20170214 EMC seminar on
wgrib2api](https://www.youtube.com/watch?v=udbxfC1V2DI).

## Disclaimer

The United States Department of Commerce (DOC) GitHub project code is
provided on an "as is" basis and the user assumes responsibility for
its use. DOC has relinquished control of the information and no longer
has responsibility to protect the integrity, confidentiality, or
availability of the information. Any claims against the Department of
Commerce stemming from the use of its GitHub project will be governed
by all applicable Federal law. Any reference to specific commercial
products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of
Commerce. The Department of Commerce seal and logo, or the seal and
logo of a DOC bureau, shall not be used in any manner to imply
endorsement of any commercial product or activity by DOC or the United
States Government.
