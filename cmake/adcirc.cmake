
IF(BUILD_ADCIRC)

    SET( ADCIRC_SOURCES  ${CMAKE_SOURCE_DIR}/src/sizes.F 
                         ${CMAKE_SOURCE_DIR}/thirdparty/KDTREE2/kdtree2.F
                         ${CMAKE_SOURCE_DIR}/src/global.F 
                         ${CMAKE_SOURCE_DIR}/src/boundaries.F 
                         ${CMAKE_SOURCE_DIR}/src/mesh.F
                         ${CMAKE_SOURCE_DIR}/src/global_3dvs.F 
                         ${CMAKE_SOURCE_DIR}/src/harm.F 
                         ${CMAKE_SOURCE_DIR}/wind/vortex.F
                         ${CMAKE_SOURCE_DIR}/src/wind.F 
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
                         ${CMAKE_SOURCE_DIR}/src/hashtable.F
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
                         ${CMAKE_SOURCE_DIR}/src/driver.F
                         ${CMAKE_SOURCE_DIR}/src/sponge_layer.F 
                         ${CMAKE_SOURCE_DIR}/src/quadtrature.F 
                         ${CMAKE_SOURCE_DIR}/src/couple2baroclinic3D.F)

    ADD_EXECUTABLE(adcirc ${ADCIRC_SOURCES})
    SET(ADCIRC_COMPILER_FLAGS "${ADDITIONAL_FLAGS_ADCIRC} ${ADCIRC_OPTION_FLAGS}")
    addCompilerFlags(adcirc ${ADDITIONAL_FLAGS_ADCIRC})
    addLibVersion(adcirc)
    INSTALL(TARGETS adcirc RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

ENDIF(BUILD_ADCIRC)
