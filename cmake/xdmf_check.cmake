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

if(NOT ENABLE_OUTPUT_XDMF)
  return()
endif()

include(CheckFortranSourceCompiles)

# ######################################################################################################################
# Find HDF5
# ######################################################################################################################

find_package(HDF5 COMPONENTS C HL)

if(NOT HDF5_FOUND)
  message(
    SEND_ERROR
      "XDMF requires HDF5 but it was not found. Specify HDF5 location or disable XDMF (-DENABLE_OUTPUT_XDMF=OFF).")
  return()
endif()

message(STATUS "HDF5 found: ${HDF5_LIBRARIES}")

# ######################################################################################################################
# Find XDMF libraries
# ######################################################################################################################

set(_xdmf_hints "")
if(XDMFHOME)
  list(
    APPEND
    _xdmf_hints
    "${XDMFHOME}/lib"
    "${XDMFHOME}/lib/x86_64-linux-gnu"
    "${XDMFHOME}/lib64")
endif()

find_library(
  XDMF_LibXdmfCore
  NAMES XdmfCore
  HINTS ${_xdmf_hints}
  DOC "XDMF Core library")

find_library(
  XDMF_LibXdmfUtils
  NAMES XdmfUtils
  HINTS ${_xdmf_hints}
  DOC "XDMF Utils library")

find_library(
  XDMF_LibXdmf
  NAMES Xdmf
  HINTS ${_xdmf_hints}
  DOC "XDMF library")

set(_xdmf_include_hints "")
if(XDMFHOME)
  list(APPEND _xdmf_include_hints "${XDMFHOME}/include")
endif()

find_path(
  XDMF_INCLUDE_DIR
  NAMES Xdmf.f
  HINTS ${_xdmf_include_hints}
  DOC "XDMF include directory")

set(XDMF_AdditionalLibs
    ""
    CACHE STRING "Additional libraries that may be required for XDMF")

if(NOT XDMF_LibXdmfCore
   OR NOT XDMF_LibXdmfUtils
   OR NOT XDMF_LibXdmf
   OR NOT XDMF_INCLUDE_DIR)
  set(_missing_components "")
  if(NOT XDMF_LibXdmfCore)
    list(APPEND _missing_components "XdmfCore library")
  endif()
  if(NOT XDMF_LibXdmfUtils)
    list(APPEND _missing_components "XdmfUtils library")
  endif()
  if(NOT XDMF_LibXdmf)
    list(APPEND _missing_components "Xdmf library")
  endif()
  if(NOT XDMF_INCLUDE_DIR)
    list(APPEND _missing_components "Xdmf.f include file")
  endif()

  set(_error_msg "XDMF components not found.\nMissing: ${_missing_components}")
  if(XDMFHOME)
    set(_error_msg "${_error_msg}\nSearched in: ${XDMFHOME}")
  else()
    set(_error_msg "${_error_msg}\nSpecify -DXDMFHOME=<path> to set XDMF location")
  endif()
  set(_error_msg "${_error_msg}\nOr disable XDMF (-DENABLE_OUTPUT_XDMF=OFF)")

  message(SEND_ERROR "${_error_msg}")
  return()
endif()

message(STATUS "XDMF include: ${XDMF_INCLUDE_DIR}")
message(STATUS "XDMF Core library: ${XDMF_LibXdmfCore}")
message(STATUS "XDMF Utils library: ${XDMF_LibXdmfUtils}")
message(STATUS "XDMF library: ${XDMF_LibXdmf}")

# ######################################################################################################################
# Test that XDMF works
# ######################################################################################################################

set(_cmake_required_libraries_save ${CMAKE_REQUIRED_LIBRARIES})
set(_cmake_required_includes_save ${CMAKE_REQUIRED_INCLUDES})

set(CMAKE_REQUIRED_INCLUDES "${XDMF_INCLUDE_DIR}")
set(CMAKE_REQUIRED_LIBRARIES
    ${XDMF_LibXdmfCore}
    ${XDMF_LibXdmfUtils}
    ${XDMF_LibXdmf}
    ${XDMF_AdditionalLibs}
    ${HDF5_LIBRARIES})

CHECK_Fortran_SOURCE_COMPILES(
  "
  PROGRAM XDMFCHECK
    IMPLICIT NONE
    INCLUDE 'Xdmf.f'
    INTEGER :: xdmfunit
    CALL xdmfInit(xdmfunit)
  END PROGRAM
  "
  XDMF_FORTRAN_WORKS
  SRC_EXT
  F90)

set(CMAKE_REQUIRED_LIBRARIES ${_cmake_required_libraries_save})
set(CMAKE_REQUIRED_INCLUDES ${_cmake_required_includes_save})

if(XDMF_FORTRAN_WORKS)
  set(XDMF_WORKING TRUE)
  message(STATUS "XDMF is compatible with the current compiler")
else()
  set(XDMF_WORKING FALSE)
  message(
    SEND_ERROR "XDMF library is not compatible with ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}.\n"
               "Ensure the XDMF library was compiled with the same Fortran compiler.\n"
               "Try a different XDMFHOME or disable XDMF output (-DENABLE_OUTPUT_XDMF=OFF).")
  return()
endif()

if(XDMF_WORKING AND NOT TARGET XDMF::XDMF)
  add_library(XDMF::XDMF INTERFACE IMPORTED)

  set(_xdmf_all_libs ${XDMF_LibXdmfCore} ${XDMF_LibXdmfUtils} ${XDMF_LibXdmf})
  if(XDMF_AdditionalLibs)
    list(APPEND _xdmf_all_libs ${XDMF_AdditionalLibs})
  endif()
  list(APPEND _xdmf_all_libs ${HDF5_LIBRARIES})

  set_target_properties(
    XDMF::XDMF PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${XDMF_INCLUDE_DIR};${CMAKE_CURRENT_SOURCE_DIR}/src"
                          INTERFACE_LINK_LIBRARIES "${_xdmf_all_libs}")
endif()
