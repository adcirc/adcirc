IF(BUILD_ADCPREP)
    
    SET( ADCPREP_SOURCES  ${CMAKE_SOURCE_DIR}/src/sizes.F 
                          ${CMAKE_SOURCE_DIR}/thirdparty/KDTREE2/kdtree2.F 
                          ${CMAKE_SOURCE_DIR}/src/global.F 
                          ${CMAKE_SOURCE_DIR}/src/boundaries.F 
                          ${CMAKE_SOURCE_DIR}/src/hashtable.F
                          ${CMAKE_SOURCE_DIR}/src/mesh.F 
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
                          ${CMAKE_SOURCE_DIR}/src/rainfall.F
                          ${CMAKE_SOURCE_DIR}/src/gwce.F
                          ${CMAKE_SOURCE_DIR}/src/wetdry.F 
                          ${CMAKE_SOURCE_DIR}/src/momentum.F 
                          ${CMAKE_SOURCE_DIR}/src/netcdfio.F 
                          ${CMAKE_SOURCE_DIR}/prep/prep.F 
                          ${CMAKE_SOURCE_DIR}/prep/interp.F 
                          ${CMAKE_SOURCE_DIR}/prep/machdep.F
                          ${CMAKE_SOURCE_DIR}/src/sponge_layer.F 
                          ${CMAKE_SOURCE_DIR}/src/quadtrature.F 
                          ${CMAKE_SOURCE_DIR}/src/couple2baroclinic3D.F )

    ADD_EXECUTABLE(adcprep ${ADCPREP_SOURCES})

    addCompilerFlags(adcprep ${ADDITIONAL_FLAGS_ADCPREP})
    addLibMkdir(adcprep)
    addLibMetis(adcprep)

    IF(BUILD_PADCSWAN OR BUILD_PUNSWAN)
        TARGET_COMPILE_DEFINITIONS(adcprep PRIVATE ${PREP_SWAN_FLAG})
    ENDIF(BUILD_PADCSWAN OR BUILD_PUNSWAN)

    TARGET_INCLUDE_DIRECTORIES(adcprep PRIVATE prep)

    INSTALL(TARGETS adcprep RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

ENDIF(BUILD_ADCPREP)
