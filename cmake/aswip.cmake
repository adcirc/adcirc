if(BUILD_ASWIP)

  set(ASWIP_SOURCES
      ${CMAKE_SOURCE_DIR}/src/sizes.F
      ${CMAKE_SOURCE_DIR}/src/constants.F
      ${CMAKE_SOURCE_DIR}/src/global.F
      ${CMAKE_SOURCE_DIR}/src/global_3dvs.F
      ${CMAKE_SOURCE_DIR}/src/boundaries.F
      ${CMAKE_SOURCE_DIR}/src/hashtable.F
      ${CMAKE_SOURCE_DIR}/src/gl2loc_mapping.F
      ${CMAKE_SOURCE_DIR}/src/internaltide.F
      ${CMAKE_SOURCE_DIR}/src/nodalattr.F
      ${CMAKE_SOURCE_DIR}/src/mesh.F
      ${CMAKE_SOURCE_DIR}/src/wind.F
      ${CMAKE_SOURCE_DIR}/src/owiwind.F
      ${CMAKE_SOURCE_DIR}/src/owiwind_netcdf.F
      ${CMAKE_SOURCE_DIR}/thirdparty/KDTREE2/kdtree2.F
      ${CMAKE_SOURCE_DIR}/src/owi_ice.F
      ${CMAKE_SOURCE_DIR}/wind/vortex.F
      ${CMAKE_SOURCE_DIR}/wind/aswip.F)

  add_executable(aswip ${ASWIP_SOURCES})

  addcompilerflags(aswip ${ADDITIONAL_FLAGS_ASWIP})
  addlibversion(aswip)

  install(TARGETS aswip RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

endif(BUILD_ASWIP)
