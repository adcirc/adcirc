# ######################################################################################################################
#
# ADCIRC - The ADvanced CIRCulation model Copyright (C) 1994-2025 R.A. Luettich, Jr., J.J. Westerink
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
# GENERAL OPTIONS
# ######################################################################################################################
cmake_minimum_required(VERSION 3.16)
project(adcirc)

# RPATH CONFIGURATION
option(ADCIRC_USE_RPATH "Use RPATH to find shared libraries at runtime" ON)
mark_as_advanced(ADCIRC_USE_RPATH)
if(ADCIRC_USE_RPATH)
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  list(
    FIND
    CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
    "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
    isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
  endif()
endif()
# ######################################################################################################################

# ######################################################################################################################
# COMPILERS/LANGUAGE
# ######################################################################################################################

# ######################################################################################################################
# ...Perl is required to use SWAN since it generates the source files. All SWAN options will be disabled if perl cannot
# be found.
find_package(Perl)
# ######################################################################################################################

# ######################################################################################################################
# ...Language Specifications
enable_language(Fortran)
enable_language(C)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# ######################################################################################################################

# ######################################################################################################################
# ...Put the static libraries in the CMakeFiles folder so they don't contaminate the build directory
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles)
# ######################################################################################################################

# ...Set the version string (used for libraries)
set(ADCIRC_VERSION_STRING "v56.2.0")

# Parse version string into components
string(
  REGEX
  REPLACE "^v"
          ""
          ADCIRC_VERSION_NUMERIC
          "${ADCIRC_VERSION_STRING}")
string(
  REPLACE "."
          ";"
          VERSION_LIST
          ${ADCIRC_VERSION_NUMERIC})
list(
  GET
  VERSION_LIST
  0
  ADCIRC_VERSION_MAJOR)
list(
  GET
  VERSION_LIST
  1
  ADCIRC_VERSION_MINOR)
list(
  GET
  VERSION_LIST
  2
  ADCIRC_VERSION_PATCH)

# Set module path for cleaner include() statements
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ...Include helper functions to build components
include(helper_functions)

# ...Determine architecture specific parameters
include(architecture)

# ...Generate the selection list of options
include(options_select)

# ...Translate the selected options list into compiler flags
include(options_flags)

# ...Ensure that the netCDF libraries are working
include(netcdf_check)

# GRIB2 requires DATETIME, establish dependency
if(ENABLE_GRIB2 AND NOT ENABLE_DATETIME)
  message(STATUS "GRIB2 support requires datetime library, enabling ENABLE_DATETIME")
  set(ENABLE_DATETIME
      TRUE
      CACHE BOOL "Enable datetime library (required by GRIB2)" FORCE)
endif()

# Build datetime_fortran if enabled
if(ENABLE_DATETIME)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/datetime_fortran EXCLUDE_FROM_ALL)
endif()

# Build wgrib2 if GRIB2 is enabled (depends on datetime)
if(ENABLE_GRIB2)
  set(JAS_ENABLE_SHARED
      FALSE
      CACHE BOOL "Enable shared JAS")
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/wgrib2 EXCLUDE_FROM_ALL)
endif()

# ...Include the macros for developer mode
include(developer)

# ...Ensure that the XDMF libraries are working
include(xdmf_check)

# ...Build the version module
include(version)

# ...Build the mkdir library
include(mkdir)

# ...Create interface libraries for compiler flags and dependencies
adcirc_create_interface_libraries()

# ...Build ADCIRC
include(adcirc)

# ...Build METIS
include(metis)

# ...Build ADCPREP
include(adcprep)

# ...Build ADCSWAN (Serial ADCIRC+SWAN)
include(adcswan)

# ...Build SWAN (Serial structured SWAN)
include(swan)

# ...Build PADCIRC (Parallel ADCIRC)
include(padcirc)

# ...Build PADCSWAN (Parallel ADCIRC+SWAN)
include(padcswan)

# ...Build PUNSWAN (Parallel Unstructured SWAN)
include(punswan)

# ...Build ASWIP
include(aswip)

# ...Build ADCIRC Utility Codes
include(utilities)

# ...Build LIBADC
include(libadcirc)
