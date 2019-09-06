
IF(BUILD_ADCIRC)

    SET( ADCIRC_SOURCES  src/sizes.F KDTREE2/kdtree2.F
                          src/global.F src/boundaries.F src/mesh.F
                          src/global_3dvs.F src/harm.F wind/vortex.F
                          src/wind.F src/owiwind.F src/rs2.F
                          src/owi_ice.F src/itpackv.F src/nodalattr.F 
                          src/globalio.F src/subdomain.F src/gwce.F
                          src/wetdry.F src/momentum.F src/netcdfio.F
                          src/control.F src/xdmfio.F src/hashtable.F
                          src/write_output.F src/couple2swan.F src/adcirc.F
                          src/weir_boundary.F src/read_input.F src/cstart.F
                          src/hstart.F src/timestep.F src/vsmy.F
                          src/transport.F src/driver.F )

    ADD_EXECUTABLE(adcirc ${ADCIRC_SOURCES})
    SET(ADCIRC_COMPILER_FLAGS "${ADDITIONAL_FLAGS_ADCIRC} ${ADCIRC_OPTION_FLAGS}")
    addCompilerFlags(adcirc ${ADDITIONAL_FLAGS_ADCIRC})
    addLibVersion(adcirc)
    INSTALL(TARGETS adcirc RUNTIME DESTINATION bin)

ENDIF(BUILD_ADCIRC)
