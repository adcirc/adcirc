# ##############################################################################
# CMake Build File for ADCIRC(+SWAN)
#
# Written By: Zach Cobell
#
# ##############################################################################
#
# The CMake build system enable ADCIRC (and SWAN) to be deployed and built in a
# cross platform environment.
#
# ##############################################################################

# ...Set the default build type
if(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      ${CMAKE_BUILD_TYPE}
      CACHE STRING "Choose the type of
        build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
        Release RelWithDebInfo MinSizeRel.")
else()
  set(CMAKE_BUILD_TYPE
      Release
      CACHE STRING "Choose the type of build,
        options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
        RelWithDebInfo MinSizeRel.")
endif()

# ##############################################################################
# GENERAL OPTIONS
# ##############################################################################
cmake_minimum_required(VERSION 2.8.12)
project(adcirc)

# ##############################################################################

# ##############################################################################
# COMPILERS/LANGUAGE
# ##############################################################################

# ##############################################################################
# ...Perl Perl is required to use SWAN since it generates the source files. All
# SWAN options will be disabled if perl cannot be found. The user is given a
# chance to specify its location
find_program(PERL perl)
if(${PERL} STREQUAL "PERL-NOTFOUND")
  set(PERL_FOUND FALSE)
else()
  set(PERL_FOUND TRUE)
endif()
# ##############################################################################

# ##############################################################################
# ...Language Specifications
enable_language(Fortran)
enable_language(C)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
# ##############################################################################

# ##############################################################################
# ...Put the static libraries in the CMakeFiles folder so they don't contaminate
# the build directory
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles)
# ##############################################################################

# ...Set the version string (used for libraries)
set(ADCIRC_VERSION_MAJOR 55)
set(ADCIRC_VERSION_MINOR 1)
set(ADCIRC_VERSION_PATCH 0)
set(ADCIRC_VERSION_STRING
    ${ADCIRC_VERSION_MAJOR}.${ADCIRC_VERSION_MINOR}.${ADCIRC_VERSION_PATCH})

# ...Include Macros to build components
include(${CMAKE_SOURCE_DIR}/cmake/macros.cmake)

# ...Determine architecture specific parameters
include(${CMAKE_SOURCE_DIR}/cmake/architecture.cmake)

# ...Generate the selection list of options
include(${CMAKE_SOURCE_DIR}/cmake/options_select.cmake)

# ...Translate the selected options list into compiler flags
include(${CMAKE_SOURCE_DIR}/cmake/options_flags.cmake)

# ...Ensure that the netCDF libraries are working
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(${CMAKE_SOURCE_DIR}/cmake/netcdf_check.cmake)

if(ENABLE_GRIB2)
  set(JAS_ENABLE_SHARED
      FALSE
      CACHE BOOL "Enable shared JAS")
  set(ENABLE_DATETIME
      TRUE
      CACHE BOOL "Enable datetime library")
  add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/wgrib2 EXCLUDE_FROM_ALL)
endif()

if(ENABLE_DATETIME)
  add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/datetime_fortran
                   EXCLUDE_FROM_ALL)
endif()

# ...Ensure that the XDMF libraries are working
include(${CMAKE_SOURCE_DIR}/cmake/xdmf_check.cmake)

# ...Build the version module
include(${CMAKE_SOURCE_DIR}/cmake/version.cmake)

# ...Build the mkdir library
include(${CMAKE_SOURCE_DIR}/cmake/mkdir.cmake)

# ...Build ADCIRC
include(${CMAKE_SOURCE_DIR}/cmake/adcirc.cmake)

# ...Build METIS
include(${CMAKE_SOURCE_DIR}/cmake/metis.cmake)

# ...Build ADCPREP
include(${CMAKE_SOURCE_DIR}/cmake/adcprep.cmake)

# ...Build ADCSWAN (Serial ADCIRC+SWAN)
include(${CMAKE_SOURCE_DIR}/cmake/adcswan.cmake)

# ...Build SWAN (Serial structured SWAN)
include(${CMAKE_SOURCE_DIR}/cmake/swan.cmake)

# ...Build PADCIRC (Parallel ADCIRC)
include(${CMAKE_SOURCE_DIR}/cmake/padcirc.cmake)

# ...Build PADCSWAN (Parallel ADCIRC+SWAN)
include(${CMAKE_SOURCE_DIR}/cmake/padcswan.cmake)

# ...Build PUNSWAN (Parallel Unstructured SWAN)
include(${CMAKE_SOURCE_DIR}/cmake/punswan.cmake)

# ...Build ASWIP
include(${CMAKE_SOURCE_DIR}/cmake/aswip.cmake)

# ...Build ADCIRC Utility Codes
include(${CMAKE_SOURCE_DIR}/cmake/utilities.cmake)

# ...Build LIBADC
include(${CMAKE_SOURCE_DIR}/cmake/libadcirc.cmake)
