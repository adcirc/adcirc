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

if(NOT ENABLE_OUTPUT_NETCDF)
  return()
endif()

include(CheckFortranSourceCompiles)

# ######################################################################################################################
# Handle NetCDF paths
# ######################################################################################################################

if(NETCDFHOME)
  set(_netcdfhome_value "${NETCDFHOME}")

  message(
    WARNING "NETCDFHOME is deprecated because netcdf-c and netcdf-fortran are often installed separately.\n"
            "Consider using these CMake variables instead:\n"
            "  When both are in the same location:\n"
            "    -DNETCDF_DIR=<path>              : Root directory for both NetCDF-C and NetCDF-Fortran\n"
            "  When installed separately:\n"
            "    -DNETCDF_DIR=<path>              : Root directory for NetCDF-C\n"
            "    -DNETCDF_F90_ROOT=<path>         : Root directory for NetCDF-Fortran\n"
            "  Or set specific paths:\n"
            "    -DNETCDF_INCLUDE_DIR=<path>      : NetCDF-C include directory\n"
            "    -DNETCDF_LIBRARY=<path>          : NetCDF-C library\n"
            "    -DNETCDF_F90_INCLUDE_DIR=<path>  : NetCDF-Fortran include directory\n"
            "    -DNETCDF_F90_LIBRARY=<path>      : NetCDF-Fortran library\n"
            "\nFor backwards compatibility, using NETCDFHOME=${_netcdfhome_value} for both C and Fortran.")

  if(NOT NETCDF_DIR)
    set(NETCDF_DIR "${_netcdfhome_value}")
  endif()
  if(NOT NETCDF_F90_ROOT)
    set(NETCDF_F90_ROOT "${_netcdfhome_value}")
  endif()
endif()

set(NETCDF_F90 "YES")
find_package(NetCDF)

set(NETCDF_AdditionalLibs
    ""
    CACHE STRING "Additional libraries that may be required for netCDF")

if(NOT NETCDF_FOUND)
  message(
    SEND_ERROR "NetCDF not found. Please specify NetCDF paths using:\n"
               "  When both are in the same location:\n"
               "    -DNETCDF_DIR=<path>              : Root directory for both NetCDF-C and NetCDF-Fortran\n"
               "  When installed separately:\n"
               "    -DNETCDF_DIR=<path>              : Root directory for NetCDF-C\n"
               "    -DNETCDF_F90_ROOT=<path>         : Root directory for NetCDF-Fortran\n"
               "  Or set specific paths:\n"
               "    -DNETCDF_INCLUDE_DIR=<path>      : NetCDF-C include directory\n"
               "    -DNETCDF_LIBRARY=<path>          : NetCDF-C library\n"
               "    -DNETCDF_F90_INCLUDE_DIR=<path>  : NetCDF-Fortran include directory\n"
               "    -DNETCDF_F90_LIBRARY=<path>      : NetCDF-Fortran library")
  return()
endif()

message(STATUS "NetCDF-C include: ${NETCDF_C_INCLUDE_DIRS}")
message(STATUS "NetCDF-C library: ${NETCDF_C_LIBRARIES}")
if(NETCDF_F90_INCLUDE_DIRS)
  message(STATUS "NetCDF-Fortran include: ${NETCDF_F90_INCLUDE_DIRS}")
endif()
if(NETCDF_F90_LIBRARIES)
  message(STATUS "NetCDF-Fortran library: ${NETCDF_F90_LIBRARIES}")
endif()

# ######################################################################################################################
# Test that NetCDF-Fortran works
# ######################################################################################################################

set(_cmake_required_libraries_save ${CMAKE_REQUIRED_LIBRARIES})
set(_cmake_required_includes_save ${CMAKE_REQUIRED_INCLUDES})

set(CMAKE_REQUIRED_INCLUDES ${NETCDF_INCLUDE_DIRS})
set(CMAKE_REQUIRED_LIBRARIES ${NETCDF_LIBRARIES} ${NETCDF_AdditionalLibs})

CHECK_Fortran_SOURCE_COMPILES(
  "
  PROGRAM netCDF3Test
    USE NETCDF
    IMPLICIT NONE
    INTEGER :: IERR, NCID
    IERR = NF90_OPEN('test.nc', NF90_NOWRITE, NCID)
  END PROGRAM
  "
  NETCDF3_FORTRAN_WORKS
  SRC_EXT
  F90)

CHECK_Fortran_SOURCE_COMPILES(
  "
  PROGRAM netCDF4Test
    USE NETCDF
    IMPLICIT NONE
    INTEGER :: IERR, NCID, VARID
    IERR = NF90_DEF_VAR_DEFLATE(NCID, VARID, 1, 1, 2)
  END PROGRAM
  "
  NETCDF4_FORTRAN_WORKS
  SRC_EXT
  F90)

set(CMAKE_REQUIRED_LIBRARIES ${_cmake_required_libraries_save})
set(CMAKE_REQUIRED_INCLUDES ${_cmake_required_includes_save})

if(NETCDF3_FORTRAN_WORKS)
  set(NETCDF_WORKING TRUE)
  message(STATUS "NetCDF-Fortran is compatible with the current compiler")

  if(NETCDF4_FORTRAN_WORKS)
    set(NETCDF4_WORKING TRUE)
    message(STATUS "NetCDF4 with compression support detected")
  else()
    set(NETCDF4_WORKING FALSE)
    message(STATUS "NetCDF4 compression not available (NetCDF3 API only)")
  endif()
else()
  set(NETCDF_WORKING FALSE)
  message(
    SEND_ERROR
      "NetCDF-Fortran library is not compatible with ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}.\n"
      "Ensure the NetCDF library was compiled with the same Fortran compiler.\n"
      "Try a different NetCDF installation (see NETCDF_DIR/NETCDF_F90_ROOT above) or disable NetCDF output (-DENABLE_OUTPUT_NETCDF=OFF)."
  )
  return()
endif()

if(NETCDF_WORKING AND NOT TARGET NetCDF::NetCDF)
  add_library(NetCDF::NetCDF INTERFACE IMPORTED)

  set_target_properties(
    NetCDF::NetCDF PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_INCLUDE_DIRS}"
                              INTERFACE_LINK_LIBRARIES "${NETCDF_LIBRARIES};${NETCDF_AdditionalLibs}")
endif()
