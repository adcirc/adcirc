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
set(LIBADC_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/sizes.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/constants.F
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/KDTREE2/kdtree2.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/global.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/boundaries.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/global_3dvs.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/messenger.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/mesh.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/vew1d.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/harm.F
    ${CMAKE_CURRENT_SOURCE_DIR}/wind/vortex.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/wind.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/wind_modules/nws08.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/hashtable.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/wind_modules/nws12.F90
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
    ${CMAKE_CURRENT_SOURCE_DIR}/src/weir_boundary.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/read_input.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cstart.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/hstart.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/timestep.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/vsmy.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/transport.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/sponge_layer.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/quadrature.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/couple2baroclinic3D.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/gl2loc_mapping.F
    ${CMAKE_CURRENT_SOURCE_DIR}/src/internaltide.F
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
      ${CMAKE_CURRENT_SOURCE_DIR}/src/wind_modules/nws13.F90
      ${CMAKE_CURRENT_SOURCE_DIR}/src/netcdfio.F
      ${CMAKE_CURRENT_SOURCE_DIR}/src/netcdf_error.F90)
endif()

if(XDMF_WORKING)
  set(LIBADC_SOURCES ${LIBADC_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/src/xdmfio.F)
endif()

if(BUILD_LIBADCIRC_STATIC)

  add_library(libadcirc_static STATIC ${LIBADC_SOURCES})
  install(TARGETS libadcirc_static ARCHIVE DESTINATION lib)
  add_dependencies(libadcirc_static version mkdir)
  target_link_libraries(libadcirc_static version mkdir)
  set_target_properties(libadcirc_static PROPERTIES OUTPUT_NAME "adcirc")
  adcirc_add_compiler_flags(libadcirc_static ${ADDITIONAL_FLAGS_ADCIRC})
  adcirc_add_mpi(libadcirc_static)
  set_target_properties(libadcirc_static PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  install(
    TARGETS libadcirc_static
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_LIBADCIRC_STATIC)

if(BUILD_LIBADCIRC_SHARED)
  set(LIBADC_SHARED_SOURCES ${LIBADC_SOURCES})
  add_library(libadcirc_shared SHARED ${LIBADC_SHARED_SOURCES})

  adcirc_add_compiler_flags(libadcirc_shared ${ADDITIONAL_FLAGS_ADCIRC})
  adcirc_add_libraries(libadcirc_shared)
  adcirc_add_mpi(libadcirc_shared)

  set_target_properties(libadcirc_shared PROPERTIES OUTPUT_NAME "adcirc")

  if(APPLE)
    set_property(TARGET libadcirc_shared PROPERTY MACOSX_RPATH ON)
  endif()

  set_target_properties(libadcirc_shared PROPERTIES VERSION ${ADCIRC_VERSION_STRING} SOVERSION ${ADCIRC_VERSION_MAJOR})
  write_basic_package_version_file(
    libadcircConfigVersion.cmake
    VERSION ${ADCIRC_VERSION_STRING}
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

  # Set the linker language to Fortran for the executable
  set_target_properties(libadcirc_shared PROPERTIES LINKER_LANGUAGE Fortran)

endif(BUILD_LIBADCIRC_SHARED)
