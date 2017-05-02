IF(BUILD_ADCPREP)
    
    SET( ADCPREP_SOURCES  src/sizes.F KDTREE2/kdtree2.F src/global.F src/boundaries.F src/hashtable.F
                          src/mesh.F src/global_3dvs.F wind/vortex.F src/owiwind.F src/rs2.F
                          src/owi_ice.F src/wind.F prep/presizes.F prep/pre_global.F prep/metis.F
                          prep/subprep.F prep/adcprep.F prep/decomp.F prep/prep_weir.F src/itpackv.F
                          src/nodalattr.F src/harm.F prep/read_global.F src/subdomain.F src/gwce.F
                          src/wetdry.F src/momentum.F src/netcdfio.F prep/prep.F prep/interp.F prep/machdep.F )

    ADD_EXECUTABLE(adcprep ${ADCPREP_SOURCES})

    addCompilerFlags(adcprep ${ADDITIONAL_FLAGS_ADCPREP})
    addLibMkdir(adcprep)
    addLibMetis(adcprep)

    IF(BUILD_PADCSWAN OR BUILD_PUNSWAN)
        TARGET_COMPILE_DEFINITIONS(adcprep PRIVATE ${PREP_SWAN_FLAG})
    ENDIF(BUILD_PADCSWAN OR BUILD_PUNSWAN)

    TARGET_INCLUDE_DIRECTORIES(adcprep PRIVATE prep)

    INSTALL(TARGETS adcprep RUNTIME DESTINATION bin)

ENDIF(BUILD_ADCPREP)
