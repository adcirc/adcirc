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

# ######################################################################################################################
# ...Output format options
option(ADCIRC_ENABLE_OUTPUT_NETCDF "Turn on netCDF output options" OFF)
if(ADCIRC_ENABLE_OUTPUT_NETCDF)
  option(ADICRC_ENABLE_OUTPUT_XDMF "Turn on XDMF output options" OFF)
endif(ADCIRC_ENABLE_OUTPUT_NETCDF)
# ######################################################################################################################

# ######################################################################################################################
# ...Executables
option(ADCIRC_BUILD_ADCIRC "Build the serial ADCIRC executable" OFF)

if(PERL_FOUND)
  option(ADCIRC_BUILD_ADCSWAN "Build the serial SWAN+ADCIRC executable" OFF)
  option(ADCIRC_BUILD_SWAN "Build the serial SWAN executable" OFF)
  mark_as_advanced(ADCIRC_BUILD_SWAN)
endif(PERL_FOUND)

if(MPI_FOUND)
  option(ADCIRC_BUILD_ADCPREP "Build the MPI parallel ADCIRC preprocessor" OFF)
  option(ADCIRC_BUILD_PADCIRC "Build the MPI parallel ADCIRC executable" OFF)
  option(ADCIRC_BUILD_LIBADCIRC_STATIC "Build the static library version of ADCIRC" OFF)
  option(ADCIRC_BUILD_LIBADCIRC_SHARED "Build the shared library version of ADCIRC" OFF)
  if(PERL_FOUND)
    option(ADCIRC_BUILD_PADCSWAN "Build the MPI parallel SWAN+ADCIRC executable" OFF)
    option(ADCIRC_BUILD_PUNSWAN "Build the MPI parallel unstructured SWAN executable" OFF)
    mark_as_advanced(ADCIRC_BUILD_PUNSWAN)
  endif(PERL_FOUND)
endif(MPI_FOUND)

option(ADCIRC_BUILD_ASWIP "Build ASWIP (ASymmetric Wind Input Preprocessor)")
option(ADCIRC_BUILD_UTILITIES "Build the ADCIRC utility programs" OFF)
mark_as_advanced(ADCIRC_BUILD_UTILITIES)

option(ADCIRC_STRICT_COMPILE "Enable strict compile options for ADCIRC" OFF)

# Options for GRIB2 and DATETIME. Note that if GRIB2 is selected, DATETIME is automatically selected.
option(ADCIRC_ENABLE_GRIB2 "Use GRIB2API static libraries." OFF)
if(ADCIRC_ENABLE_GRIB2)
  set(ADCIRC_ENABLE_DATETIME
      ON
      CACHE BOOL "Use DATETIME static libraries." FORCE)
endif()
option(ADCIRC_ENABLE_DATETIME "Use DATETIME static libraries." ${ADCIRC_ENABLE_DATETIME})

# ######################################################################################################################

# ######################################################################################################################
# ...No need to show CXX compilers on main screen
mark_as_advanced(
  CLEAR
  CMAKE_CXX_COMPILER
  CMAKE_C_COMPILER
  CMAKE_Fortran_COMPILER)
mark_as_advanced(
  CLEAR
  CMAKE_CXX_FLAGS_RELEASE
  CMAKE_C_FLAGS_RELEASE
  CMAKE_Fortran_FLAGS_RELEASE)
mark_as_advanced(
  CLEAR
  CMAKE_CXX_FLAGS_DEBUG
  CMAKE_C_FLAGS_DEBUG
  CMAKE_Fortran_FLAGS_DEBUG)
mark_as_advanced(
  CLEAR
  CMAKE_CXX_FLAGS_RELWITHDEBINFO
  CMAKE_C_FLAGS_RELWITHDEBINFO
  CMAKE_Fortran_FLAGS_RELWITHDEBINFO)
mark_as_advanced(
  CLEAR
  CMAKE_CXX_FLAGS_MINSIZEREL
  CMAKE_C_FLAGS_MINSIZEREL
  CMAKE_Fortran_FLAGS_MINSIZEREL)
# ######################################################################################################################

# ######################################################################################################################
# ...Library paths
if(ADCIRC_ENABLE_OUTPUT_XDMF)
  if(NOT
     ${XDMFHOME}
     STREQUAL
     "")
    set(XDMFHOME
        ${XDMFHOME}
        CACHE STRING "XDMF home path containing lib and include")
  elseif(
    NOT
    $ENV{XDMFHOME}
    STREQUAL
    "")
    set(XDMFHOME
        $ENV{XDMFHOME}
        CACHE STRING "XDMF home path containing lib and include")
  else(
    NOT
    ${XDMFHOME}
    STREQUAL
    "")
    set(XDMFHOME
        "XDMF-NOTFOUND"
        CACHE STRING "XDMF home path containing lib and include")
  endif(
    NOT
    ${XDMFHOME}
    STREQUAL
    "")
else()
  unset(XDMFHOME CACHE)
endif()
# ######################################################################################################################

# ######################################################################################################################
# ...Additional flags
set(ADCIRC_ADDITIONAL_FLAGS_ADCIRC
    ""
    CACHE STRING "Additional flags to compile ADCIRC with")
if(PERL_FOUND)
  set(ADCIRC_ADDITIONAL_FLAGS_SWAN
      ""
      CACHE STRING "Additional flags to compile SWAN with")
endif(PERL_FOUND)
set(ADCIRC_ADDITIONAL_FLAGS_ADCPREP
    ""
    CACHE STRING "Additional flags to compile ADCPREP with")
set(ADCIRC_ADDITIONAL_FLAGS_ASWIP
    ""
    CACHE STRING "Additional flags to compile ASWIP with")
set(ADCIRC_ADDITIONAL_FLAGS_UTLIITIES
    ""
    CACHE STRING "Additional flags for utility programs")
mark_as_advanced(
  ADCIRC_ADDITIONAL_FLAGS_ADCIRC
  ADCIRC_ADDITIONAL_FLAGS_SWAN
  ADCIRC_ADDITIONAL_FLAGS_ADCPREP
  ADCIRC_ADDITIONAL_FLAGS_ASWIP
  ADCIRC_ADDITIONAL_FLAGS_UTLIITIES)
# ######################################################################################################################

# ######################################################################################################################
# ...Options enabled via compiler flags within the code
option(ADCIRC_ENABLE_WARN_ELEV_DEBUG "Enable writing of the fort.69 debug file" OFF)
if(ADCIRC_ENABLE_WARN_ELEV_DEBUG)
  message(
    WARNING
      "The compile time enabled fort.69 file is deprecated and the user should use the &warnElevControl namelist instead. This option will be removed in a future release."
  )
endif()
mark_as_advanced(ADCIRC_ENABLE_WARN_ELEV_DEBUG)

option(ADCIRC_IBM "Format code for IBM based architectures" OFF)
option(ADCIRC_SGI "Format code for SGI based architectures" OFF)
option(ADCIRC_SUN "Format code for SUN based architectures" OFF)
option(ADCIRC_CRAY "Format code for CRAY based architectures" OFF)
option(ADCIRC_CRAYX1 "Format code for CRAYX1 based architectures" OFF)
mark_as_advanced(
  ADCIRC_IBM
  ADCIRC_SGI
  ADCIRC_SUN
  ADCIRC_CRAY
  ADCIRC_CRAYX1)

option(ADCIRC_DEBUG_FULL_STACK "Write the detailed stack trace during debugging" OFF)
option(ADCIRC_DEBUG_FLUSH_MESSAGES "Do not allow caching of screen printed messages" OFF)
option(ADCIRC_DEBUG_LOG_LEVEL "Force debug log level for screen messages" OFF)
option(ADCIRC_DEBUG_ALL_TRACE "Write all tracing debug information to screen output during debugging" OFF)
option(ADCIRC_DEBUG_GLOBALIO_TRACE "Write tracing debug information from the GLOBALIO module" OFF)
option(ADCIRC_DEBUG_WRITER_TRACE "Write tracing debug information for writer processors" OFF)
option(ADCIRC_DEBUG_WRITE_OUTPUT_TRACE "Write tracing debug information for the write output module" OFF)
option(ADCIRC_DEBUG_WIND_TRACE "Write tracing debug information for the wind module" OFF)
option(ADCIRC_DEBUG_WEIR_TRACE "Write tracing debug information for the weir module" OFF)
option(ADCIRC_DEBUG_TVW_TRACE "Write tracing debug information for the time varying weir module" OFF)
option(ADCIRC_DEBUG_VSMY_TRACE "Write tracing debug information for the 3D momentum equation module" OFF)
option(ADCIRC_DEBUG_TIMESTEP_TRACE "Write the tracing debug information for the timestepping module" OFF)
option(ADCIRC_DEBUG_SUBPREP_TRACE "Write the tracing debug information for the subdomain prep module" OFF)
option(ADCIRC_DEBUG_SUBDOMAIN_TRACE "Write the tracing debug information for the subdomain module" OFF)
option(ADCIRC_DEBUG_READ_INPUT_TRACE "Write the tracing debug information for the read_input module" OFF)
option(ADCIRC_DEBUG_OWIWIND_TRACE "Write the tracing debug information for the OWI wind module" OFF)
option(ADCIRC_DEBUG_NODALATTR_TRACE "Write the tracing debug information for the nodal attributes module" OFF)
option(ADCIRC_DEBUG_NETCDF_TRACE "Write the tracing debug information for the netCDF modle" OFF)
option(ADCIRC_DEBUG_MESSENGER_TRACE "Write the tracing debug information for the message passing (MPI) module" OFF)
option(ADCIRC_DEBUG_MESH_TRACE "Write the tracing debug information for the mesh module" OFF)
option(ADCIRC_DEBUG_HOTSTART_TRACE "Write the tracing debug information for the hotstart module" OFF)
option(ADCIRC_DEBUG_GLOBAL_TRACE "Write the tracing debug information for the global module" OFF)
option(ADCIRC_DEBUG_HARM_TRACE "Write the tracing debug information for the harmonics module" OFF)
option(ADCIRC_DEBUG_COLDSTART_TRACE "Write the tracing debug information for the coldstart module" OFF)
option(ADCIRC_DEBUG_COUPLE2SWAN_TRACE "Write the tracing debug information for the couple2swan module" OFF)
option(ADCIRC_DEBUG_ADCIRC_TRACE "Write the tracing debug information for the main ADCIRC module" OFF)
option(ADCIRC_DEBUG_HOLLAND "Write the debugging information for the symmetric Holland model" OFF)
option(ADCIRC_DEBUG_NWS14 "Write the debugging information for NWS=14 interpolation" OFF)
mark_as_advanced(
  ADCIRC_DEBUG_FLUSH_MESSAGES
  ADCIRC_DEBUG_FULL_STACK
  ADCIRC_DEBUG_LOG_LEVEL
  ADCIRC_DEBUG_ALL_TRACE
  ADCIRC_DEBUG_GLOBALIO_TRACE
  ADCIRC_DEBUG_WRITER_TRACE
  ADCIRC_DEBUG_WRITE_OUTPUT_TRACE
  ADCIRC_DEBUG_WIND_TRACE
  ADCIRC_DEBUG_WEIR_TRACE
  ADCIRC_DEBUG_TVW_TRACE
  ADCIRC_DEBUG_VSMY_TRACE
  ADCIRC_DEBUG_TIMESTEP_TRACE
  ADCIRC_DEBUG_SUBPREP_TRACE
  ADCIRC_DEBUG_SUBDOMAIN_TRACE
  ADCIRC_DEBUG_READ_INPUT_TRACE
  ADCIRC_DEBUG_OWIWIND_TRACE
  ADCIRC_DEBUG_NODALATTR_TRACE
  ADCIRC_DEBUG_NETCDF_TRACE
  ADCIRC_DEBUG_MESSENGER_TRACE
  ADCIRC_DEBUG_MESH_TRACE
  ADCIRC_DEBUG_HOTSTART_TRACE
  ADCIRC_DEBUG_GLOBAL_TRACE
  ADCIRC_DEBUG_HARM_TRACE
  ADCIRC_DEBUG_COLDSTART_TRACE
  ADCIRC_DEBUG_COUPLE2SWAN_TRACE
  ADCIRC_DEBUG_ADCIRC_TRACE
  ADCIRC_DEBUG_HOLLAND
  ADCIRC_DEBUG_NWS14
  ADCIRC_IBM)

option(ADCIRC_ENABLE_POWELL
       "Force Powell wind drag to be enabled. Warning: Overrides any other options specified at run time." OFF)
mark_as_advanced(ADCIRC_ENABLE_POWELL)
if(ADCIRC_ENABLE_POWELL)
  message(
    WARNING
      "The compile time enabled Powell wind drag forcing is deprecated and users should switch to the &metControl namelist. This feature will be removed in a future release"
  )
endif()

option(ADCIRC_VECTOR_COMPUTER "Assume the system is a vector computer" OFF)
mark_as_advanced(ADCIRC_VECTOR_COMPUTER)

option(ADCIRC_NOF2008 "The fortran compiler being used does not have F2008 intrinsics" OFF)
mark_as_advanced(ADCIRC_NOF2008)
