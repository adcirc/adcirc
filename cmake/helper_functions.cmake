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

function(adcirc_set_module_directory TARGET)
  set_target_properties(${TARGET} PROPERTIES Fortran_MODULE_DIRECTORY
                                             ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mod/${TARGET})
endfunction()

function(adcirc_swan_configure_adcswan)
  if(WIN32)
    add_custom_target(
      generate_swan_serial_sources
      COMMAND ${PERL_EXECUTABLE} switch.pl -adcirc -unix *.ftn *.ftn90
      COMMAND if not exist \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_serial_source\" mkdir
              \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_serial_source\"
      COMMAND move /y *.f \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_serial_source/.\"
      COMMAND move /y *.f90 \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_serial_source/.\"
      BYPRODUCTS ${SWAN1SERIAL_SOURCES} ${SWAN2SERIAL_SOURCES}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Serial unSWAN Sources...")
  else()
    add_custom_target(
      generate_swan_serial_sources
      COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_serial_source
      COMMAND ${PERL_EXECUTABLE} switch.pl -adcirc -unix -outdir
              ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_serial_source *.ftn *.ftn90
      BYPRODUCTS ${SWAN1SERIAL_SOURCES} ${SWAN2SERIAL_SOURCES}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Serial unSWAN Sources...")
  endif()

  set_source_files_properties(${SWAN1SERIAL_SOURCES} ${SWAN2SERIAL_SOURCES} PROPERTIES GENERATED TRUE)
endfunction()

function(adcirc_swan_configure_padcswan)
  if(WIN32)
    add_custom_target(
      generate_swan_parallel_sources
      COMMAND ${PERL_EXECUTABLE} switch.pl -pun -adcirc -unix *.ftn *.ftn90
      COMMAND if not exist \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_parallel_source\" mkdir
              \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_parallel_source\"
      COMMAND move /y *.f \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_parallel_source/.\"
      COMMAND move /y *.f90 \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_parallel_source/.\"
      BYPRODUCTS ${SWAN1PARALLEL_SOURCES} ${SWAN2PARALLEL_SOURCES}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Parallel unSWAN Sources...")
  else()
    add_custom_target(
      generate_swan_parallel_sources
      COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_parallel_source
      COMMAND ${PERL_EXECUTABLE} switch.pl -pun -adcirc -unix -outdir
              ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_parallel_source *.ftn *.ftn90
      BYPRODUCTS ${SWAN1PARALLEL_SOURCES} ${SWAN2PARALLEL_SOURCES}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Parallel unSWAN Sources...")
  endif()

  set_source_files_properties(${SWAN1PARALLEL_SOURCES} ${SWAN2PARALLEL_SOURCES} PROPERTIES GENERATED TRUE)
endfunction()

function(adcirc_swan_configure_serial)
  if(WIN32)
    add_custom_target(
      generate_swanonly_serial_sources
      COMMAND ${PERL_EXECUTABLE} switch.pl -unix *.ftn *.ftn90
      COMMAND if not exist \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source\" mkdir
              \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source\"
      COMMAND move /y *.f \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/.\"
      COMMAND move /y *.f90 \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/.\"
      BYPRODUCTS ${SWANONLY_SERIAL_SOURCES}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Serial SWAN stand alone Sources...")
  else()
    add_custom_target(
      generate_swanonly_serial_sources
      COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source
      COMMAND ${PERL_EXECUTABLE} switch.pl -unix -outdir ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source
              *.ftn *.ftn90
      BYPRODUCTS ${SWANONLY_SERIAL_SOURCES}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Serial SWAN stand alone Sources...")
  endif()

  set_source_files_properties(${SWANONLY_SERIAL_SOURCES} PROPERTIES GENERATED TRUE)
endfunction()

function(adcirc_swan_configure_parallel)
  if(WIN32)
    add_custom_target(
      generate_swanonly_parallel_sources
      COMMAND ${PERL_EXECUTABLE} switch.pl -pun -unix *.ftn *.ftn90
      COMMAND if not exist \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source\" mkdir
              \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source\"
      COMMAND move /y *.f \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/.\"
      COMMAND move /y *.f90 \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/.\"
      BYPRODUCTS ${SWANONLY1_PARALLEL_SOURCES} ${SWANONLY2_PARALLEL_SOURCES}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Parallel SWAN Sources...")
  else()
    add_custom_target(
      generate_swanonly_parallel_sources
      COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source
      COMMAND ${PERL_EXECUTABLE} switch.pl -pun -unix -outdir
              ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source *.ftn *.ftn90
      BYPRODUCTS ${SWANONLY1_PARALLEL_SOURCES} ${SWANONLY2_PARALLEL_SOURCES}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Parallel SWAN Sources...")
  endif()

  set_source_files_properties(${SWANONLY1_PARALLEL_SOURCES} ${SWANONLY2_PARALLEL_SOURCES} PROPERTIES GENERATED TRUE)
endfunction()
