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

# ######################################################################################################################
# Create interface libraries for compiler flags, options, and dependencies This function should be called once during
# CMake configuration
# ######################################################################################################################
function(adcirc_create_interface_libraries)
  # ====================================================================================================================
  # adcirc::compiler_flags - Standard ADCIRC Fortran compiler flags
  # ====================================================================================================================
  add_library(adcirc_compiler_flags INTERFACE)
  add_library(adcirc::compiler_flags ALIAS adcirc_compiler_flags)

  if(Fortran_LINELENGTH_FLAG)
    target_compile_options(adcirc_compiler_flags INTERFACE ${Fortran_LINELENGTH_FLAG})
  endif()
  if(Fortran_COMPILER_SPECIFIC_FLAG)
    target_compile_options(adcirc_compiler_flags INTERFACE ${Fortran_COMPILER_SPECIFIC_FLAG})
  endif()

  # ====================================================================================================================
  # adcirc::option_flags - ADCIRC compile-time feature flags
  # ====================================================================================================================
  add_library(adcirc_option_flags INTERFACE)
  add_library(adcirc::option_flags ALIAS adcirc_option_flags)

  if(ADCIRC_OPTION_FLAGS)
    target_compile_definitions(adcirc_option_flags INTERFACE ${ADCIRC_OPTION_FLAGS})
  endif()
  if(WIN32)
    target_compile_definitions(adcirc_option_flags INTERFACE WINDOWS)
  endif()

  # ====================================================================================================================
  # adcirc::swan_compiler_flags - SWAN-specific compiler flags
  # ====================================================================================================================
  add_library(adcirc_swan_compiler_flags INTERFACE)
  add_library(adcirc::swan_compiler_flags ALIAS adcirc_swan_compiler_flags)

  if(Fortran_COMPILER_SPECIFIC_FLAG)
    target_compile_options(adcirc_swan_compiler_flags INTERFACE ${Fortran_COMPILER_SPECIFIC_FLAG})
  endif()

  # SWAN is third-party code, suppress warnings
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
      target_compile_options(adcirc_swan_compiler_flags INTERFACE "-fallow-argument-mismatch" "-w")
    else()
      target_compile_options(adcirc_swan_compiler_flags INTERFACE "-w")
    endif()
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel" OR CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
    target_compile_options(adcirc_swan_compiler_flags INTERFACE "-diag-disable=6843,8291")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "LLVMFlang")
    target_compile_options(adcirc_swan_compiler_flags INTERFACE "-w")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    target_compile_options(adcirc_swan_compiler_flags INTERFACE "-Minform=severe")
  endif()

  # ====================================================================================================================
  # adcirc::link_libraries - All ADCIRC library dependencies with their definitions
  # ====================================================================================================================
  add_library(adcirc_link_libraries INTERFACE)
  add_library(adcirc::link_libraries ALIAS adcirc_link_libraries)

  target_sources(adcirc_link_libraries INTERFACE $<TARGET_OBJECTS:adcirc::version> $<TARGET_OBJECTS:adcirc::mkdir>)

  target_include_directories(
    adcirc_link_libraries INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mod/adcirc_version
                                    ${CMAKE_CURRENT_SOURCE_DIR}/prep ${CMAKE_CURRENT_SOURCE_DIR}/src)

  if(NETCDF_WORKING)
    target_compile_definitions(adcirc_link_libraries INTERFACE ADCNETCDF)
    if(NETCDF4_WORKING)
      target_compile_definitions(adcirc_link_libraries INTERFACE HAVE_NETCDF4 NETCDF_CAN_DEFLATE)
    endif()
    target_link_libraries(adcirc_link_libraries INTERFACE NetCDF::NetCDF)
  endif()

  if(XDMF_WORKING)
    target_compile_definitions(adcirc_link_libraries INTERFACE ADCXDMF)
    target_link_libraries(adcirc_link_libraries INTERFACE XDMF::XDMF)
  endif()

  if(ENABLE_GRIB2)
    target_compile_definitions(adcirc_link_libraries INTERFACE GRIB2API)
    target_link_libraries(adcirc_link_libraries INTERFACE GRIB2::grib2 GRIB2::g2c GRIB2::geo)
  endif()

  if(ENABLE_DATETIME)
    target_compile_definitions(adcirc_link_libraries INTERFACE DATETIME)
    target_include_directories(adcirc_link_libraries
                               INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mod/datetime_fortran)
    target_link_libraries(adcirc_link_libraries INTERFACE datetime)
  endif()

  if(MPI_FOUND)
    add_library(adcirc_mpi_interface INTERFACE)
    add_library(adcirc::mpi ALIAS adcirc_mpi_interface)
    target_link_libraries(adcirc_mpi_interface INTERFACE MPI::MPI_Fortran)
    target_compile_definitions(adcirc_mpi_interface INTERFACE CMPI)
  endif()
endfunction()
