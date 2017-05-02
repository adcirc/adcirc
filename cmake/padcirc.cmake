
IF(BUILD_PADCIRC)
    
    SET( PADCIRC1_SOURCES  src/sizes.F KDTREE2/kdtree2.F 
                           src/global.F src/boundaries.F src/global_3dvs.F
                           src/messenger.F )

    SET( PADCIRC2_SOURCES  src/mesh.F src/harm.F wind/vortex.F src/wind.F src/hashtable.F
                           src/owiwind.F src/rs2.F src/owi_ice.F 
                           src/itpackv.F src/nodalattr.F src/globalio.F 
                           src/subdomain.F src/gwce.F src/wetdry.F src/momentum.F
                           src/netcdfio.F src/control.F src/xdmfio.F )

    SET( PADCIRC3_SOURCES  src/writer.F )

    SET( PADCIRC4_SOURCES  src/write_output.F src/couple2swan.F src/adcirc.F
                           src/weir_boundary.F src/read_input.F src/cstart.F src/hstart.F 
                           src/timestep.F src/vsmy.F src/transport.F src/driver.F )

    ADD_LIBRARY(templib_padcirc1 ${PADCIRC1_SOURCES})
    ADD_LIBRARY(templib_padcirc2 ${PADCIRC2_SOURCES})
    ADD_LIBRARY(templib_padcirc3 ${PADCIRC3_SOURCES})
    
    ADD_EXECUTABLE(padcirc ${PADCIRC4_SOURCES})

    addCompilerFlags(templib_padcirc1)
    addCompilerFlags(templib_padcirc2)
    addCompilerFlags(templib_padcirc3)
    addCompilerFlags(padcirc)

    addMPI(templib_padcirc1)
    addMPI(templib_padcirc2)
    addMPI(templib_padcirc3)
    addMPI(padcirc)
    
    addLibMkdir(templib_padcirc1)
    addLibMkdir(templib_padcirc2)
    addLibMkdir(templib_padcirc3)
    addLibMkdir(padcirc)
    addLibVersion(padcirc)

    TARGET_INCLUDE_DIRECTORIES(templib_padcirc2 PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_padcirc1)
    TARGET_INCLUDE_DIRECTORIES(templib_padcirc3 PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_padcirc1)
    TARGET_INCLUDE_DIRECTORIES(padcirc          PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_padcirc1)
    TARGET_INCLUDE_DIRECTORIES(templib_padcirc3 PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_padcirc2)
    TARGET_INCLUDE_DIRECTORIES(padcirc          PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_padcirc2)
    TARGET_INCLUDE_DIRECTORIES(padcirc          PRIVATE ${CMAKE_BINARY_DIR}/CMakeFiles/mod/templib_padcirc3)

    TARGET_LINK_LIBRARIES(padcirc templib_padcirc3 templib_padcirc2 templib_padcirc1)

    ADD_DEPENDENCIES(padcirc          templib_padcirc3)
    ADD_DEPENDENCIES(templib_padcirc2 templib_padcirc1)
    ADD_DEPENDENCIES(templib_padcirc3 templib_padcirc2)
    ADD_DEPENDENCIES(templib_padcirc1 version)
    ADD_DEPENDENCIES(templib_padcirc1 mkdir)
    
    INSTALL(TARGETS padcirc RUNTIME DESTINATION bin)

ENDIF(BUILD_PADCIRC)
