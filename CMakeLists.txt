# ######################################################################################################################
#
# ADCIRC - The ADvanced CIRCulation model Copyright (C) 1994-2023 R.A. Luettich, Jr., J.J. Westerink
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General
# Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.
#
# ######################################################################################################################
# CMake Build File for ADCIRC(+SWAN)
#
# Written By: Zach Cobell
#
# ######################################################################################################################
#
# The CMake build system enable ADCIRC (and SWAN) to be deployed and built in a cross platform environment.
#
# ######################################################################################################################

# ######################################################################################################################
# GENERAL OPTIONS
# ######################################################################################################################
cmake_minimum_required(VERSION 3.12)
project(adcirc)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
    set(CMAKE_BUILD_TYPE
            RelWithDebInfo
            CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui, ccmake
    set_property(
            CACHE CMAKE_BUILD_TYPE
            PROPERTY STRINGS
            "Debug"
            "Release"
            "MinSizeRel"
            "RelWithDebInfo")
endif()

# ######################################################################################################################

# ######################################################################################################################
# COMPILERS/LANGUAGE
# ######################################################################################################################
enable_language(Fortran)
enable_language(C)
enable_language(CXX)
cmake_policy(SET CMP0069 NEW)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/StandardProjectSettings.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/PreventInSourceBuilds.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/InterproceduralOptimization.cmake)

# ######################################################################################################################
# ...Perl Perl is required to use SWAN since it generates the source files. All SWAN options will be disabled if perl
# cannot be found. The user is given a chance to specify its location
find_program(PERL perl)
if(${PERL} STREQUAL "PERL-NOTFOUND")
  set(PERL_FOUND FALSE)
else()
  set(PERL_FOUND TRUE)
endif()
# ######################################################################################################################

# ######################################################################################################################
# ...Language Specifications
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
# ######################################################################################################################

# ######################################################################################################################
# ...Put the static libraries in the CMakeFiles folder so they don't contaminate the build directory
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles)
# ######################################################################################################################

# ...Set the version string (used for libraries)
set(ADCIRC_VERSION_MAJOR 56)
set(ADCIRC_VERSION_MINOR 0)
set(ADCIRC_VERSION_PATCH 0)
set(ADCIRC_VERSION_STRING ${ADCIRC_VERSION_MAJOR}.${ADCIRC_VERSION_MINOR}.${ADCIRC_VERSION_PATCH})

# ...Include Macros to build components
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/macros.cmake)

# ...Determine architecture specific parameters
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/architecture.cmake)

# ...Generate the selection list of options
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/options_select.cmake)

# ...Translate the selected options list into compiler flags
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/options_flags.cmake)

# ...Ensure that the netCDF libraries are working
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/netcdf_check.cmake)

if(ENABLE_GRIB2)
  set(JAS_ENABLE_SHARED
      FALSE
      CACHE BOOL "Enable shared JAS")
  set(ENABLE_DATETIME
      TRUE
      CACHE BOOL "Enable datetime library")
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/wgrib2 EXCLUDE_FROM_ALL)
endif()

if(ENABLE_DATETIME)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/datetime_fortran EXCLUDE_FROM_ALL)
endif()

# ...Ensure that the XDMF libraries are working
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/xdmf_check.cmake)

# ...Build the version module
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.cmake)

# ...Build the mkdir library
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/mkdir.cmake)

# ...Build ADCIRC
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/adcirc.cmake)

# ...Build METIS
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/metis.cmake)

# ...Build ADCPREP
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/adcprep.cmake)

# ...Build ADCSWAN (Serial ADCIRC+SWAN)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/adcswan.cmake)

# ...Build SWAN (Serial structured SWAN)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/swan.cmake)

# ...Build PADCIRC (Parallel ADCIRC)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/padcirc.cmake)

# ...Build PADCSWAN (Parallel ADCIRC+SWAN)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/padcswan.cmake)

# ...Build PUNSWAN (Parallel Unstructured SWAN)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/punswan.cmake)

# ...Build ASWIP
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/aswip.cmake)

# ...Build ADCIRC Utility Codes
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/utilities.cmake)

# ...Build LIBADC
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/libadcirc.cmake)
