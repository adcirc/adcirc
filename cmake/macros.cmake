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
macro(adcirc_add_compiler_flags TARGET)

  set(LOCAL_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG} ${ARGN}")
  set(LOCAL_COMPILER_DEFINITIONS "${ADCIRC_OPTION_FLAGS} ${PRECISION_FLAG} ${MACHINE_FLAG}")

  string(STRIP ${LOCAL_COMPILER_FLAGS} LOCAL_COMPILER_FLAGS)
  separate_arguments(LOCAL_COMPILER_DEFINITIONS)

  set_target_properties(${TARGET} PROPERTIES Fortran_MODULE_DIRECTORY
                                             ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mod/${TARGET})
  set_target_properties(${TARGET} PROPERTIES COMPILE_FLAGS ${LOCAL_COMPILER_FLAGS})
  target_compile_definitions(${TARGET} PRIVATE ${LOCAL_COMPILER_DEFINITIONS})

  adcirc_add_version_definitions(${TARGET})
  adcirc_add_mkdir_definitions(${TARGET})
  adcirc_add_netcdf_definitions(${TARGET})
  adcirc_add_xdmf_definitions(${TARGET})
  adcirc_add_grib2_definitions(${TARGET})
  adcirc_add_datetime_definitions(${TARGET})

  set(LOCAL_COMPILER_FLAGS "")
  set(LOCAL_COMPILER_DEFINITIONS "")

endmacro()

macro(adcirc_add_libraries TARGET)
  adcirc_add_version_library(${TARGET})
  adcirc_add_mkdir_library(${TARGET})
  adcirc_add_datetime_libraries(${TARGET})
  adcirc_add_netcdf_libraries(${TARGET})
  adcirc_add_xdmf_libraries(${TARGET})
  adcirc_add_grib2_libraries(${TARGET})
endmacro()

macro(adcirc_add_compiler_flags_swan TARGET)

  set(LOCAL_COMPILER_FLAGS "${Fortran_COMPILER_SPECIFIC_FLAG} ${ARGN}")

  # SWAN is not under our control, so to get a warning-free build, we need to suppress warning generated in their code
  if(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER 10 OR ${CMAKE_Fortran_COMPILER_VERSION} VERSION_EQUAL 10)
      set(SWAN_WARNING_FLAGS "-fallow-argument-mismatch -w")
      set(LOCAL_COMPILER_FLAGS "${LOCAL_COMPILER_FLAGS} ${SWAN_WARNING_FLAGS}")
    endif()
  elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel" OR ${CMAKE_Fortran_COMPILER_ID} MATCHES "IntelLLVM")
    set(SWAN_WARNING_FLAGS "-diag-disable 6843 -diag-disable 8291")
    set(LOCAL_COMPILER_FLAGS "${LOCAL_COMPILER_FLAGS} ${SWAN_WARNING_FLAGS}")
  elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "LLVMFlang")
    set(SWAN_WARNING_FLAGS "-w")
    set(LOCAL_COMPILER_FLAGS "${LOCAL_COMPILER_FLAGS} ${SWAN_WARNING_FLAGS}")
  elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")
    set(SWAN_WARNING_FLAGS "-Minform=severe")
    set(LOCAL_COMPILER_FLAGS "${LOCAL_COMPILER_FLAGS} ${SWAN_WARNING_FLAGS}")
  endif()

  set_target_properties(${TARGET} PROPERTIES Fortran_MODULE_DIRECTORY
                                             ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mod/${TARGET})
  if(NOT
     "${LOCAL_COMPILER_FLAGS}"
     STREQUAL
     "")
    set_target_properties(${TARGET} PROPERTIES COMPILE_FLAGS ${LOCAL_COMPILER_FLAGS})
  endif(
    NOT
    "${LOCAL_COMPILER_FLAGS}"
    STREQUAL
    "")

  set(LOCAL_COMPILER_FLAGS "")

endmacro()

macro(adcirc_add_grib2_definitions TARGET)
  if(ENABLE_GRIB2)
    target_compile_definitions(${TARGET} PRIVATE GRIB2API)
    target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mod/grib2)
    add_dependencies(
      ${TARGET}
      grib2
      g2c
      geo)
  endif()
endmacro()

macro(adcirc_add_grib2_libraries TARGET)
  if(ENABLE_GRIB2)
    target_link_libraries(
      ${TARGET}
      grib2
      g2c
      geo)
  endif()
endmacro()

macro(adcirc_add_datetime_definitions TARGET)
  target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/mod)
  add_dependencies(${TARGET} fdate::fdate)
endmacro()

macro(adcirc_add_datetime_libraries TARGET)
  target_link_libraries(${TARGET} fdate::fdate)
endmacro()

macro(adcirc_add_netcdf_definitions TARGET)
  if(NETCDF_WORKING)
    target_compile_definitions(${TARGET} PRIVATE ADCNETCDF)
    if(NETCDF4_WORKING)
      target_compile_definitions(${TARGET} PRIVATE HAVE_NETCDF4)
      target_compile_definitions(${TARGET} PRIVATE NETCDF_CAN_DEFLATE)
    endif()
    target_include_directories(${TARGET} PRIVATE ${NETCDF_INCLUDE_DIRS})
  endif()
endmacro()

macro(adcirc_add_netcdf_libraries TARGET)
  if(NETCDF_WORKING)
    target_link_libraries(${TARGET} ${NETCDF_LIBRARIES} ${NETCDF_AdditionalLibs})
  endif()
endmacro()

macro(adcirc_add_xdmf_definitions TARGET)
  if(XDMF_WORKING)
    target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_include_directories(${TARGET} PRIVATE ${XDMFHOME}/include)
    target_compile_definitions(${TARGET} PRIVATE ADCXDMF)
  endif()
endmacro()

macro(adcirc_add_xdmf_libraries TARGET)
  if(XDMF_WORKING)
    target_link_libraries(
      ${TARGET}
      ${XDMF_LibXdmfCore}
      ${XDMF_LibXdmfUtils}
      ${XDMF_LibXdmf}
      ${XDMF_AdditionalLibs})
  endif()
endmacro()

macro(adcirc_add_version_definitions TARGET)
  target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_mod)
  add_dependencies(${TARGET} version)
endmacro()

macro(adcirc_add_version_library TARGET)
  target_link_libraries(${TARGET} version)
endmacro()

macro(adcirc_add_mkdir_definitions TARGET)
  target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/prep)
  add_dependencies(${TARGET} mkdir)
endmacro()

macro(adcirc_add_mkdir_library TARGET)
  target_link_libraries(${TARGET} mkdir)
endmacro()

macro(adcirc_add_metis_library TARGET)
  target_link_libraries(${TARGET} metis)
  add_dependencies(${TARGET} metis)
endmacro()

macro(adcirc_add_mpi TARGET)
  target_include_directories(${TARGET} PRIVATE ${MPI_Fortran_INCLUDE_PATH})
  target_compile_definitions(${TARGET} PRIVATE CMPI)
  target_compile_definitions(${TARGET} PRIVATE ${MPIMOD_FLAG})
  target_link_libraries(${TARGET} ${MPI_Fortran_LIBRARIES})
endmacro()

macro(adcirc_swan_configure_adcswan)
  if(WIN32)
    add_custom_command(
      OUTPUT ${SWAN1SERIAL_SOURCES} ${SWAN2SERIAL_SOURCES}
      COMMAND ${PERL} switch.pl -adcirc -unix *.ftn *.ftn90
      COMMAND if not exist \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_serial_source\" mkdir
              \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_serial_source\"
      COMMAND move /y *.f \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_serial_source/.\"
      COMMAND move /y *.f90 \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_serial_source/.\"
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Serial unSWAN Sources...")
  else(WIN32)
    add_custom_command(
      OUTPUT ${SWAN1SERIAL_SOURCES} ${SWAN2SERIAL_SOURCES}
      COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_serial_source
      COMMAND ${PERL} switch.pl -adcirc -unix -outdir ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_serial_source *.ftn
              *.ftn90
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Serial unSWAN Sources...")
  endif(WIN32)
endmacro()

macro(adcirc_swan_configure_padcswan)
  if(WIN32)
    add_custom_command(
      OUTPUT ${SWAN1PARALLEL_SOURCES} ${SWAN2PARALLEL_SOURCES}
      COMMAND ${PERL} switch.pl -pun -adcirc -unix *.ftn *.ftn90
      COMMAND if not exist \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_parallel_source\" mkdir
              \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_parallel_source\"
      COMMAND move /y *.f \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_parallel_source/.\"
      COMMAND move /y *.f90 \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_parallel_source/.\"
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Serial unSWAN Sources...")
  else(WIN32)
    add_custom_command(
      OUTPUT ${SWAN1PARALLEL_SOURCES} ${SWAN2PARALLEL_SOURCES}
      COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_parallel_source
      COMMAND ${PERL} switch.pl -pun -adcirc -unix -outdir ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_parallel_source
              *.ftn *.ftn90
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Parallel unSWAN Sources...")
  endif(WIN32)
endmacro()

macro(adcirc_swan_configure_serial)
  if(WIN32)
    add_custom_command(
      OUTPUT ${SWANONLY_SERIAL_SOURCES}
      COMMAND ${PERL} switch.pl -unix *.ftn *.ftn90
      COMMAND if not exist \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source\" mkdir
              \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source\"
      COMMAND move /y *.f \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/.\"
      COMMAND move /y *.f90 \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source/.\"
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Serial SWAN stand alone Sources...")
  else(WIN32)
    add_custom_command(
      OUTPUT ${SWANONLY_SERIAL_SOURCES}
      COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source
      COMMAND ${PERL} switch.pl -unix -outdir ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source *.ftn
              *.ftn90
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Serial SWAN stand alone Sources...")
  endif(WIN32)
endmacro()

macro(adcirc_swan_configure_parallel)
  if(WIN32)
    add_custom_command(
      OUTPUT ${SWANONLY1_PARALLEL_SOURCES} ${SWANONLY2_PARALLEL_SOURCES}
      COMMAND ${PERL} switch.pl -pun -unix *.ftn *.ftn90
      COMMAND if not exist \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source\" mkdir
              \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source\"
      COMMAND move /y *.f \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/.\"
      COMMAND move /y *.f90 \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source/.\"
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Parallel SWAN Sources...")
  else(WIN32)
    add_custom_command(
      OUTPUT ${SWANONLY1_PARALLEL_SOURCES} ${SWANONLY2_PARALLEL_SOURCES}
      COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source
      COMMAND ${PERL} switch.pl -pun -unix -outdir ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source *.ftn
              *.ftn90
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Parallel unSWAN Sources...")
  endif(WIN32)
endmacro()
