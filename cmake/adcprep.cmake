if(BUILD_ADCPREP)

  set(ADCPREP_SOURCES
      ${CMAKE_SOURCE_DIR}/src/sizes.F
      ${CMAKE_SOURCE_DIR}/thirdparty/KDTREE2/kdtree2.F
      ${CMAKE_SOURCE_DIR}/src/global.F
      ${CMAKE_SOURCE_DIR}/src/boundaries.F
      ${CMAKE_SOURCE_DIR}/src/hashtable.F
      ${CMAKE_SOURCE_DIR}/src/mesh.F
      ${CMAKE_SOURCE_DIR}/src/vew1d.F
      ${CMAKE_SOURCE_DIR}/src/global_3dvs.F
      ${CMAKE_SOURCE_DIR}/wind/vortex.F
      ${CMAKE_SOURCE_DIR}/src/owiwind.F
      ${CMAKE_SOURCE_DIR}/src/owiwind_netcdf.F
      ${CMAKE_SOURCE_DIR}/src/rs2.F
      ${CMAKE_SOURCE_DIR}/src/owi_ice.F
      ${CMAKE_SOURCE_DIR}/src/wind.F
      ${CMAKE_SOURCE_DIR}/prep/presizes.F
      ${CMAKE_SOURCE_DIR}/prep/pre_global.F
      ${CMAKE_SOURCE_DIR}/prep/metis.F
      ${CMAKE_SOURCE_DIR}/prep/subprep.F
      ${CMAKE_SOURCE_DIR}/prep/adcprep.F
      ${CMAKE_SOURCE_DIR}/prep/decomp.F
      ${CMAKE_SOURCE_DIR}/prep/prep_weir.F
      ${CMAKE_SOURCE_DIR}/src/itpackv.F
      ${CMAKE_SOURCE_DIR}/src/nodalattr.F
      ${CMAKE_SOURCE_DIR}/src/harm.F
      ${CMAKE_SOURCE_DIR}/prep/read_global.F
      ${CMAKE_SOURCE_DIR}/src/subdomain.F
      ${CMAKE_SOURCE_DIR}/src/gwce.F
      ${CMAKE_SOURCE_DIR}/src/wetdry.F
      ${CMAKE_SOURCE_DIR}/src/momentum.F
      ${CMAKE_SOURCE_DIR}/src/netcdfio.F
      ${CMAKE_SOURCE_DIR}/prep/prep.F
      ${CMAKE_SOURCE_DIR}/prep/interp.F
      ${CMAKE_SOURCE_DIR}/prep/machdep.F
      ${CMAKE_SOURCE_DIR}/src/sponge_layer.F
      ${CMAKE_SOURCE_DIR}/src/quadrature.F
      ${CMAKE_SOURCE_DIR}/src/gl2loc_mapping.F
      ${CMAKE_SOURCE_DIR}/src/couple2baroclinic3D.F
      ${CMAKE_SOURCE_DIR}/src/internaltide.F)

  add_executable(adcprep ${ADCPREP_SOURCES})

  addcompilerflags(adcprep ${ADDITIONAL_FLAGS_ADCPREP})
  addlibmkdir(adcprep)
  addlibmetis(adcprep)

  if(BUILD_PADCSWAN OR BUILD_PUNSWAN)
    target_compile_definitions(adcprep PRIVATE ${PREP_SWAN_FLAG})
  endif(BUILD_PADCSWAN OR BUILD_PUNSWAN)

  target_include_directories(adcprep PRIVATE prep)

  install(TARGETS adcprep RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

endif(BUILD_ADCPREP)
