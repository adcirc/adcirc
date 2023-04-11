set(LIBADC_SOURCES
    ${CMAKE_SOURCE_DIR}/src/sizes.F
    ${CMAKE_SOURCE_DIR}/thirdparty/KDTREE2/kdtree2.F
    ${CMAKE_SOURCE_DIR}/src/global.F
    ${CMAKE_SOURCE_DIR}/src/boundaries.F
    ${CMAKE_SOURCE_DIR}/src/global_3dvs.F
    ${CMAKE_SOURCE_DIR}/src/messenger.F
    ${CMAKE_SOURCE_DIR}/src/mesh.F
    ${CMAKE_SOURCE_DIR}/src/vew1d.F
    ${CMAKE_SOURCE_DIR}/src/harm.F
    ${CMAKE_SOURCE_DIR}/wind/vortex.F
    ${CMAKE_SOURCE_DIR}/src/wind.F
    ${CMAKE_SOURCE_DIR}/src/hashtable.F
    ${CMAKE_SOURCE_DIR}/src/owiwind.F
    ${CMAKE_SOURCE_DIR}/src/owiwind_netcdf.F
    ${CMAKE_SOURCE_DIR}/src/rs2.F
    ${CMAKE_SOURCE_DIR}/src/owi_ice.F
    ${CMAKE_SOURCE_DIR}/src/itpackv.F
    ${CMAKE_SOURCE_DIR}/src/nodalattr.F
    ${CMAKE_SOURCE_DIR}/src/globalio.F
    ${CMAKE_SOURCE_DIR}/src/subdomain.F
    ${CMAKE_SOURCE_DIR}/src/gwce.F
    ${CMAKE_SOURCE_DIR}/src/wetdry.F
    ${CMAKE_SOURCE_DIR}/src/momentum.F
    ${CMAKE_SOURCE_DIR}/src/netcdfio.F
    ${CMAKE_SOURCE_DIR}/src/control.F
    ${CMAKE_SOURCE_DIR}/src/xdmfio.F
    ${CMAKE_SOURCE_DIR}/src/writer.F
    ${CMAKE_SOURCE_DIR}/src/write_output.F
    ${CMAKE_SOURCE_DIR}/src/couple2swan.F
    ${CMAKE_SOURCE_DIR}/src/adcirc.F
    ${CMAKE_SOURCE_DIR}/src/weir_boundary.F
    ${CMAKE_SOURCE_DIR}/src/read_input.F
    ${CMAKE_SOURCE_DIR}/src/cstart.F
    ${CMAKE_SOURCE_DIR}/src/hstart.F
    ${CMAKE_SOURCE_DIR}/src/timestep.F
    ${CMAKE_SOURCE_DIR}/src/vsmy.F
    ${CMAKE_SOURCE_DIR}/src/transport.F
    ${CMAKE_SOURCE_DIR}/src/sponge_layer.F
    ${CMAKE_SOURCE_DIR}/src/quadrature.F
    ${CMAKE_SOURCE_DIR}/src/couple2baroclinic3D.F)

if(BUILD_LIBADCIRC_STATIC)

  add_library(libadcirc_static STATIC ${LIBADC_SOURCES})
  install(TARGETS libadcirc_static ARCHIVE DESTINATION lib)
  add_dependencies(libadcirc_static version mkdir)
  target_link_libraries(libadcirc_static version mkdir)
  set_target_properties(libadcirc_static PROPERTIES OUTPUT_NAME "adcirc")
  addcompilerflags(libadcirc_static)
  addmpi(libadcirc_static)
  set_target_properties(libadcirc_static PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

  install(
    TARGETS libadcirc_static
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

endif(BUILD_LIBADCIRC_STATIC)

if(BUILD_LIBADCIRC_SHARED)
  set(LIBMKDIR2_SOURCES prep/mkdir.c)
  set(LIBADC_SHARED_SOURCES ${LIBADC_SOURCES} ${CMAKE_BINARY_DIR}/CMakeFiles/version_cmake.F)

  add_library(mkdir2 STATIC ${LIBMKDIR2_SOURCES})
  add_library(libadcirc_shared SHARED ${LIBADC_SHARED_SOURCES})

  addcompilerflags(libadcirc_shared)
  addmpi(libadcirc_shared)

  add_dependencies(libadcirc_shared version mkdir2)
  target_link_libraries(libadcirc_shared mkdir2)

  set_target_properties(libadcirc_shared PROPERTIES OUTPUT_NAME "adcirc")

  set_property(TARGET libadcirc_shared PROPERTY POSITION_INDEPENDENT_CODE ON)
  set_property(TARGET mkdir2 PROPERTY POSITION_INDEPENDENT_CODE ON)

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
      CODE "EXECUTE_PROCESS (COMMAND \"${CMAKE_COMMAND}\" -E copy_directory \"${CMAKE_BINARY_DIR}/CMakeFiles/mod/libadcirc_shared\" \"${CMAKE_INSTALL_PREFIX}/include/adcirc\")"
    )
  else()
    install(
      CODE "EXECUTE_PROCESS (COMMAND \"${CMAKE_COMMAND}\" -E copy_directory \"${CMAKE_BINARY_DIR}/CMakeFiles/mod/libadcirc_shared\" \"${CMAKE_INSTALL_PREFIX}/include\")"
    )
  endif()

endif(BUILD_LIBADCIRC_SHARED)
