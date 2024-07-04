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
macro(addCompilerFlags TARGET)

  set(LOCAL_COMPILER_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG}")
  set(LOCAL_COMPILER_DEFINITIONS "${ADCIRC_OPTION_FLAGS} ${PRECISION_FLAG} ${MACHINE_FLAG}")

  string(STRIP ${LOCAL_COMPILER_FLAGS} LOCAL_COMPILER_FLAGS)
  separate_arguments(LOCAL_COMPILER_DEFINITIONS)

  addlibversion(${TARGET})
  addNetCDF(${TARGET})
  addXDMF(${TARGET})
  addGrib2(${TARGET})
  addDatetime(${TARGET})

  set_target_properties(${TARGET} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mod/${TARGET})
  set_target_properties(${TARGET} PROPERTIES COMPILE_FLAGS ${LOCAL_COMPILER_FLAGS})
  target_compile_definitions(${TARGET} PRIVATE ${LOCAL_COMPILER_DEFINITIONS})

  set(LOCAL_COMPILER_FLAGS "")
  set(LOCAL_COMPILER_DEFINITIONS "")

endmacro(addCompilerFlags)

macro(addCompilerFlagsSwan TARGET)

  set(LOCAL_COMPILER_FLAGS "${Fortran_COMPILER_SPECIFIC_FLAG}")

  get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
  if(${Fortran_COMPILER_NAME} MATCHES "gfortran.*")
    if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER 10 OR ${CMAKE_Fortran_COMPILER_VERSION} VERSION_EQUAL 10)
      set(LOCAL_COMPILER_FLAGS "${LOCAL_COMPILER_FLAGS} -fallow-argument-mismatch")
    endif()
  elseif(${Fortran_COMPILER_NAME} MATCHES "ifx.*" OR ${Fortran_COMPILER_NAME} MATCHES "ifort.*")
    set(LOCAL_COMPILER_FLAGS "${LOCAL_COMPILER_FLAGS} -diag-disable 6843 -diag-disable 8291")
  endif()

  set_target_properties(${TARGET} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mod/${TARGET})
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

endmacro(addCompilerFlagsSwan)

macro(addGrib2 TARGET)
  if(ENABLE_GRIB2)
    target_compile_definitions(${TARGET} PRIVATE GRIB2API)
    target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mod/grib2)
    target_link_libraries(
      ${TARGET}
      grib2
      ipolates
      geo)
    add_dependencies(
      ${TARGET}
      grib2
      ipolates
      geo)
  endif(ENABLE_GRIB2)
endmacro(addGrib2)

macro(addDatetime TARGET)
  if(ENABLE_DATETIME)
    target_compile_definitions(${TARGET} PRIVATE DATETIME)
    target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mod/datetime_fortran)
    target_link_libraries(${TARGET} datetime)
    add_dependencies(${TARGET} datetime)
  endif(ENABLE_DATETIME)
endmacro(addDatetime)

macro(addNetCDF TARGET)
  if(NETCDF_WORKING)
    target_compile_definitions(${TARGET} PRIVATE ADCNETCDF)
    if(NETCDF4_WORKING)
      target_compile_definitions(${TARGET} PRIVATE HAVE_NETCDF4)
      target_compile_definitions(${TARGET} PRIVATE NETCDF_CAN_DEFLATE)
    endif(NETCDF4_WORKING)
    target_include_directories(${TARGET} PRIVATE ${NETCDF_INCLUDE_DIRS})
    target_link_libraries(${TARGET} ${NETCDF_LIBRARIES} ${NETCDF_AdditionalLibs})
  endif(NETCDF_WORKING)
endmacro(addNetCDF)

macro(addXDMF TARGET)
  if(XDMF_WORKING)
    target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_include_directories(${TARGET} PRIVATE ${XDMFHOME}/include)
    target_compile_definitions(${TARGET} PRIVATE ADCXDMF)
    target_link_libraries(
      ${TARGET}
      ${XDMF_LibXdmfCore}
      ${XDMF_LibXdmfUtils}
      ${XDMF_LibXdmf}
      ${XDMF_AdditionalLibs})
  endif(XDMF_WORKING)
endmacro(addXDMF)

macro(addLibVersion TARGET)
  target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_mod)
  target_link_libraries(${TARGET} version)
  add_dependencies(${TARGET} version)
endmacro(addLibVersion)

macro(addLibMkdir TARGET)
  target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/prep)
  target_link_libraries(${TARGET} mkdir)
  add_dependencies(${TARGET} mkdir)
endmacro(addLibMkdir)

macro(addLibMetis TARGET)
  target_link_libraries(${TARGET} metis)
  add_dependencies(${TARGET} metis)
endmacro(addLibMetis)

macro(addMPI TARGET)
  target_include_directories(${TARGET} PRIVATE ${MPI_Fortran_INCLUDE_PATH})
  target_compile_definitions(${TARGET} PRIVATE CMPI)
  target_compile_definitions(${TARGET} PRIVATE ${MPIMOD_FLAG})
  target_link_libraries(${TARGET} ${MPI_Fortran_LIBRARIES})
endmacro(addMPI)

macro(swanConfigureAdcswan)
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
      COMMAND ${PERL} switch.pl -adcirc -unix -outdir ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_serial_source *.ftn *.ftn90
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Serial unSWAN Sources...")
  endif(WIN32)
endmacro(swanConfigureAdcswan)

macro(swanConfigurePadcswan)
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
      COMMAND ${PERL} switch.pl -pun -adcirc -unix -outdir ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_parallel_source *.ftn
              *.ftn90
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Parallel unSWAN Sources...")
  endif(WIN32)
endmacro(swanConfigurePadcswan)

macro(swanConfigureSerial)
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
      COMMAND ${PERL} switch.pl -unix -outdir ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_serial_source *.ftn *.ftn90
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Serial SWAN stand alone Sources...")
  endif(WIN32)
endmacro(swanConfigureSerial)

macro(swanConfigureParallel)
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
      COMMAND ${PERL} switch.pl -pun -unix -outdir ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swanonly_parallel_source *.ftn *.ftn90
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/swan
      COMMENT "Generating Parallel unSWAN Sources...")
  endif(WIN32)
endmacro(swanConfigureParallel)
