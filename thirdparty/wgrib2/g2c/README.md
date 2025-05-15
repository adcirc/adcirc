![Status](https://github.com/NOAA-EMC/NCEPLIBS-g2c/workflows/developer/badge.svg)

# NCEPLIBS-g2c

This library contains C decoder/encoder routines for GRIB edition 2.

GRIdded Binary or General Regularly-distributed Information in Binary
form (GRIB) is a data format for meteorological and forecast data,
standardized by the World Meteorological Organization (WMO). GRIB
edition 2 (GRIB2) was approved by the WMO is 2003.

This library is part of the
[NCEPLIBS](https://github.com/NOAA-EMC/NCEPLIBS) project.

For complete documentation see the latest [NCEPLIBS-g2c
documentation](https://noaa-emc.github.io/NCEPLIBS-g2c/). For more
about GRIB2 see the [NCEP WMO GRIB2
Documentation](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/). For
the WMO GRIB2 templates and tables see [WMO Information Management
GRIB2 GitHub repository](https://github.com/wmo-im/GRIB2).

The NCEPLIBS-g2c library is used by [wgrib2](https://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/),
[grib2io](https://github.com/NOAA-MDL/grib2io), [GRaDS](http://cola.gmu.edu/grads/), and 
[Model Evaluation Tools (MET)](https://metplus.readthedocs.io/en/latest/) projects, among
others.

To submit bug reports, feature requests, or other code-related issues
including installation and usage questions, please create a [GitHub
issue](https://github.com/NOAA-EMC/NCEPLIBS-g2c/issues). For general
NCEPLIBS inquiries, contact [Ed
Hartnett](mailto:edward.hartnett@noaa.gov) (secondary point of contact
[Alex Richert](mailto:alexander.richert@noaa.gov)).

## Related NCEPLIBS Projects

Repository | Notes
-----------|------
[NCEPLIBS-w3emc](https://github.com/NOAA-EMC/NCEPLIBS-w3emc) | Coders/decoders for GRIB1
[NCEPLIBS-g2](https://github.com/NOAA-EMC/NCEPLIBS-g2) | Fortran implementation of the GRIB 2 functions
[NCEPLIBS-grib_util](https://github.com/NOAA-EMC/NCEPLIBS-grib_util) | A collection of GRIB1 and GRIB2 utilities
[NCEPLIBS-g2tmpl](https://github.com/NOAA-EMC/NCEPLIBS-g2tmpl) | Utilities for GRIB2 templates

## Authors

Wesley Ebisuzaki, Eric Engle, Stephen Gilbert, Harry Glahn, Edward
Hartnett, Dusan Jovic, Boi Vuong, other NOAA scientists and engineers.

Code Manager: [Hang Lei](mailto:hang.lei@noaa.gov), [Ed
Hartnett](mailto:edward.hartnett@noaa.gov)

## Prerequisites

- [libjasper.a](http://www.ece.uvic.ca/~mdadams/jasper/) - This
  library is a C implementation of the JPEG-2000 Part-1 standard
  (i.e., ISO/IEC 15444-1). Tested version: jasper-1.900.1. More
  information about JPEG2000 can be found at
  http://www.jpeg.org/JPEG2000.html.

- [libpng.a](http://www.libpng.org/pub/png/libpng.html) - This library
  is a C implementation of the Portable Network Graphics PNG image
  compression format. Tested version: libpng-1.2.44. More information
  about PNG can be found at http://www.libpng.org/pub/png/.

- [libz.a](http://www.gzip.org/zlib/) - This library contains
  compression/decompression routines used by libpng.a for PNG image
  compression support. Tested version: zlib-1.2.6.

- [openjpeg.a](https://www.openjpeg.org/) - OpenJPEG is an open-source
  JPEG 2000 codec written in C language. OpenJPEG is only used if
  CMake build option USE_OpenJPEG is turned on.

- [libaec.a](https://gitlab.dkrz.de/k202009/libaec) - LibAEC is the 
  Adaptive Entropy Coding library.  This library implements 
  extended Golomb-Rice coding as defined in the CCSDS recommended standard [121.0-B-3](https://public.ccsds.org/Pubs/121x0b3.pdf). 
  The library covers the adaptive entropy coder and the preprocessor discussed in
  sections 1 to 5.2.6 of the standard.

## Building

```console
git clone https://github.com/NOAA-EMC/NCEPLIBS-g2c
cmake -S NCEPLIBS-g2c -B NCEPLIBS-g2c/build # -DCMAKE_PREFIX_PATH=/usr/local/jasper-3.0.5 -DCMAKE_INSTALL_PREFIX=/path/to/install/g2c <add'l CMake options>
cmake --build NCEPLIBS-g2c/build --parallel 4
ctest --test-dir NCEPLIBS-g2c/build --parallel 4 # <add'l CTest options>
# Install to CMAKE_INSTALL_PREFIX (/usr/local by default):
cmake --install NCEPLIBS-g2c/build
```

See the [documentation](https://noaa-emc.github.io/NCEPLIBS-g2c/) for a list
of CMake options.

The NCEPLIBS-g2c library supports the PNG and JPEG2000 methods of image compression
algorithms within the GRIB2 standard.

By default the library uses Jasper for JPEG functionality, use 
`-DUSE_OpenJPEG=ON` to use the OpenJPEG library instead.

NCEPLIBS-g2c is also available through [Spack](https://spack.io) as '[g2c](https://github.com/spack/spack/tree/develop/var/spack/repos/builtin/packages/g2c)'.

## References

Hartnett, E., Ator, J, Lei, H., Richert, A., Woollen, J., King, A.,
Hartnett, A., [NCEPLIBS GRIB and BUFR Libraries: Maintaining and
Modernizing NOAA's Libraries for WMO Data
Formats](https://www.researchgate.net/publication/376390180_NCEPLIBS_GRIB_and_BUFR_Libraries_Maintaining_and_Modernizing_NOAA's_Libraries_for_WMO_Data_Formats),
American Geophysical Union (AGU) 2023. (See also
[poster](https://www.researchgate.net/publication/376582005_Poster_-_IN51B-0416_NCEPLIBS_GRIB_and_BUFR_Libraries_Maintaining_and_Modernizing_NOAA's_Libraries_for_WMO_Data_Formats)).

Hartnett, E., Lei, H., Curtis, B, Gerheiser K., [Presentation -
Improving Documentation, Testing, Process, and Code for Legacy NOAA
GRIB2 C Fortran
Libraries](https://www.researchgate.net/publication/360757566_Presentation_-_Improving_Documentation_Testing_Process_and_Code_for_Legacy_NOAA_GRIB2_C_Fortran_Libraries),
NCAR Improving Scientific Software, April 2022.  .

Kumar, V. Krishna, Gilbert, Stephen A., [GRIB2 conversion and its
usage at NCEP](docs/GRIB2_conversion_and_its_usage_at_NCEP.pdf), 14-18
November 2005, 10th Workshop on Meteorological Operational Systems
ECMWF User Orientation, retrieved on July 27, 2021 from
https://www.yumpu.com/en/document/view/11925806/grib2-conversion-and-its-usage-at-ncep.

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
