
IF(BUILD_PADCIRC)
    
    SET( PADCIRC_SOURCES  src/sizes.F KDTREE2/kdtree2.F 
                          src/global.F src/boundaries.F src/global_3dvs.F
                          src/messenger.F src/mesh.F src/harm.F wind/vortex.F 
                          src/wind.F src/hashtable.F src/owiwind.F src/rs2.F src/owi_ice.F 
                          src/itpackv.F src/nodalattr.F src/globalio.F 
                          src/subdomain.F src/gwce.F src/wetdry.F src/momentum.F
                          src/netcdfio.F src/control.F src/xdmfio.F src/writer.F 
                          src/write_output.F src/couple2swan.F src/adcirc.F
                          src/weir_boundary.F src/read_input.F src/cstart.F src/hstart.F 
                          src/timestep.F src/vsmy.F src/transport.F src/driver.F 
                          src/sponge_layer.F src/quadtrature.F src/couple2baroclinic3D.F)

    ADD_EXECUTABLE(padcirc ${PADCIRC_SOURCES})

    addCompilerFlags(padcirc)
    addMPI(padcirc)
    addLibMkdir(padcirc)
    addLibVersion(padcirc)

    ADD_DEPENDENCIES(padcirc version mkdir)
    
    INSTALL(TARGETS padcirc RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

ENDIF(BUILD_PADCIRC)
