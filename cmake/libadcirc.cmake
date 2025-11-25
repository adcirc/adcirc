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
set(LIBADC_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/sizes.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/constants.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/terminate.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/io.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/logging.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/KDTREE2/kdtree2.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/global.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/boundaries.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/global_3dvs.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/messenger.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/mesh.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/vew1d.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/harm.F
    ${CMAKE_CURRENT_SOURCE_DIR}/wind/vortex.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/wind.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nws08.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/hashtable.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/owiwind.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/rs2.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/owi_ice.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/itpackv.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/nodalattr.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/globalio.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/subdomain.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/gwce.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/wetdry.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/momentum.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/control.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/writer.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/write_output.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/couple2swan.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/adcirc.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/weir_boundary.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/read_input.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cstart.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/hstart.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/timestep.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/vsmy.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/transport.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/sponge_layer.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/quadrature.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/couple2baroclinic3D.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/gl2loc_mapping.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/internaltide.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/astronomic.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ephemerides.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/tidalpotential.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/sun.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/moon.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/sun_moon_system.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/normal_flow_boundary.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/subgridLookup.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/couple2baroclinic3D.F)

if(NETCDF_WORKING)
  set(LIBADC_SOURCES
      ${LIBADC_SOURCES}
      ${CMAKE_CURRENT_SOURCE_DIR}/src/owiwind_netcdf.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/netcdfio.F90
      ${CMAKE_CURRENT_SOURCE_DIR}/src/netcdf_error.F90)
endif()

if(XDMF_WORKING)
  set(LIBADC_SOURCES ${LIBADC_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/src/xdmfio.F)
endif()

if(BUILD_LIBADCIRC_STATIC AND MPI_FOUND)

  add_library(libadcirc_static STATIC ${LIBADC_SOURCES})
  install(TARGETS libadcirc_static ARCHIVE DESTINATION lib)
  set_target_properties(libadcirc_static PROPERTIES OUTPUT_NAME "adcirc")

  # Configure compiler flags and link libraries
  adcirc_set_module_directory(libadcirc_static)
  target_link_libraries(
    libadcirc_static
    PRIVATE adcirc::compiler_flags
            adcirc::option_flags
            adcirc::link_libraries
            adcirc_mpi_interface)

  if(ADDITIONAL_FLAGS_ADCIRC)
    target_compile_options(libadcirc_static PRIVATE ${ADDITIONAL_FLAGS_ADCIRC})
  endif()
  set_target_properties(libadcirc_static PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  install(
    TARGETS libadcirc_static
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif()

if(BUILD_LIBADCIRC_SHARED AND MPI_FOUND)
  set(LIBADC_SHARED_SOURCES ${LIBADC_SOURCES})
  add_library(libadcirc_shared SHARED ${LIBADC_SHARED_SOURCES})

  # Configure compiler flags and link libraries
  adcirc_set_module_directory(libadcirc_shared)
  target_link_libraries(
    libadcirc_shared
    PRIVATE adcirc::compiler_flags
            adcirc::option_flags
            adcirc::link_libraries
            adcirc::mpi)

  if(ADDITIONAL_FLAGS_ADCIRC)
    target_compile_options(libadcirc_shared PRIVATE ${ADDITIONAL_FLAGS_ADCIRC})
  endif()

  set_target_properties(libadcirc_shared PROPERTIES OUTPUT_NAME "adcirc")

  set_property(TARGET libadcirc_shared PROPERTY POSITION_INDEPENDENT_CODE ON)

  if(APPLE)
    set_property(TARGET libadcirc_shared PROPERTY MACOSX_RPATH ON)
  endif()

  set_target_properties(libadcirc_shared PROPERTIES VERSION ${ADCIRC_VERSION_NUMERIC} SOVERSION ${ADCIRC_VERSION_MAJOR})
  write_basic_package_version_file(
    libadcircConfigVersion.cmake
    VERSION ${ADCIRC_VERSION_NUMERIC}
    COMPATIBILITY SameMajorVersion)

  install(
    TARGETS libadcirc_shared
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libadcircConfigVersion.cmake DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake)

  if(${CMAKE_INSTALL_PREFIX} STREQUAL "/usr/local")
    install(
      CODE "EXECUTE_PROCESS (COMMAND \"${CMAKE_COMMAND}\" -E copy_directory \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mod/libadcirc_shared\" \"${CMAKE_INSTALL_PREFIX}/include/adcirc\")"
    )
  else()
    install(
      CODE "EXECUTE_PROCESS (COMMAND \"${CMAKE_COMMAND}\" -E copy_directory \"${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/mod/libadcirc_shared\" \"${CMAKE_INSTALL_PREFIX}/include\")"
    )
  endif()

  # Conditionally enable strict compiler flags for developers
  enable_developer_mode(${LIBADC_SOURCES})

endif()
